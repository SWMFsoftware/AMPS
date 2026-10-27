//======================================================================================
// Fiber.cpp -- ucontext-based cooperative fibers (see Fiber.h).
//======================================================================================

#include "Fiber.h"

#include <ucontext.h>
#include <cstdint>
#include <sys/mman.h>
#include <unistd.h>

#include <algorithm>
#include <atomic>
#include <chrono>
#include <memory>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <deque>
#include <stdexcept>
#include <string>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

//--------------------------------------------------------------------------------------
// Context switch.  glibc's getcontext()/swapcontext() save and restore the signal mask
// with a system call on every switch; with one fiber per trajectory that is about three
// system calls per trajectory.  On x86-64 the fibers use the user-space switch below
// instead: it saves the callee-saved registers, MXCSR and the x87 control word on the
// current stack, stores the stack pointer, and resumes the target stack the same way.
// Other architectures (or -DGEOSPACE_FIBER_UCONTEXT) keep ucontext.
//--------------------------------------------------------------------------------------
#if defined(__x86_64__) && !defined(GEOSPACE_FIBER_UCONTEXT)
#define GEOSPACE_FIBER_FAST_SWITCH 1
extern "C" void geospace_fiber_switch(void** saveSp,void* loadSp);
__asm__(
    ".text\n"
    ".globl geospace_fiber_switch\n"
    ".type geospace_fiber_switch,@function\n"
    "geospace_fiber_switch:\n"
    "  pushq %rbp\n"
    "  pushq %rbx\n"
    "  pushq %r12\n"
    "  pushq %r13\n"
    "  pushq %r14\n"
    "  pushq %r15\n"
    "  subq $16,%rsp\n"
    "  stmxcsr 8(%rsp)\n"
    "  fnstcw (%rsp)\n"
    "  movq %rsp,(%rdi)\n"
    "  movq %rsi,%rsp\n"
    "  ldmxcsr 8(%rsp)\n"
    "  fldcw (%rsp)\n"
    "  addq $16,%rsp\n"
    "  popq %r15\n"
    "  popq %r14\n"
    "  popq %r13\n"
    "  popq %r12\n"
    "  popq %rbx\n"
    "  popq %rbp\n"
    "  ret\n"
    ".size geospace_fiber_switch,.-geospace_fiber_switch\n");
#else
#define GEOSPACE_FIBER_FAST_SWITCH 0
#endif

namespace Geospace {
namespace GPU {
namespace Fiber {

namespace {

struct JoinGroup {
  long long remaining;
  FiberRec* waiter;              // parent fiber, or nullptr for the top-level loop
  std::exception_ptr firstError;
  bool stop;
  // The parent of a nested loop also executes iterations and may be blocked on a
  // trajectory of its own; it may only be woken by the last child when it is blocked
  // in the join itself.
  bool waiterJoining=true;
};

} // namespace

struct FiberRec {
#if GEOSPACE_FIBER_FAST_SWITCH
  void* sp;                      // saved stack pointer (geospace_fiber_switch)
#else
  ucontext_t ctx;
#endif
  char* stackBase;               // mmap base including the guard page
  std::size_t mapSize;
  std::function<void()> body;
  enum State { Ready, Running, Blocked, Done } state;
  JoinGroup* group;
  bool slab;
};

namespace {

struct Scheduler {
#if GEOSPACE_FIBER_FAST_SWITCH
  void* mainSp=nullptr;
#else
  ucontext_t mainCtx;
#endif
  std::deque<FiberRec*> ready;
  FiberRec* current=nullptr;
  long long live=0;
  std::vector<FiberRec*> freeList;
  Stats stats{};
};

thread_local Scheduler* tSched=nullptr;
bool (*gIdleHook)()=nullptr;
std::function<void()> gServiceHook;
Stats gLastStats{};
int gMaxFibers=-1;

std::size_t StackBytes() {
  static std::size_t bytes=0;
  if (bytes==0) {
    std::size_t kb=1024;
    if (const char* e=std::getenv("AMPS_EARTH_FIBER_STACK_KB")) {
      const long v=std::atol(e);
      if (v>=64) kb=static_cast<std::size_t>(v);
    }
    bytes=kb*1024;
  }
  return bytes;
}

// Stacks are reserved with MAP_NORESERVE; only touched pages use memory.  Up to
// kGuardedLimit live fibers each stack has its own mapping with a PROT_NONE guard page
// (an overflow faults instead of corrupting a neighbour).  Every guarded stack costs
// two kernel memory mappings and vm.max_map_count is typically 65530, so larger fiber
// counts (GPU pools of 10^4-10^5 slots) carve stacks from slabs of kSlab stacks that
// share one mapping and have no guard pages.  AMPS_EARTH_FIBER_GUARD=1 forces guards.
constexpr int kGuardedLimit=16384;
constexpr int kSlab=256;
struct Slab { char* base; std::size_t bytes; };
std::vector<Slab> gSlabs;
std::vector<char*> gSlabFree;
// Fibers (with their stacks) are kept for the life of the process and reused by later
// loops: a stack page is then faulted in once, not once per batch.
std::vector<FiberRec*> gFiberPool;

bool UseGuardPages() {
  static int forced=-1;
  if (forced<0) {
    const char* e=std::getenv("AMPS_EARTH_FIBER_GUARD");
    forced=(e!=nullptr && std::atoi(e)==1) ? 1 : 0;
  }
  return forced==1 || MaxFibers()<=kGuardedLimit;
}

FiberRec* AllocateFiber(Scheduler& s) {
  if (!s.freeList.empty()) {
    FiberRec* f=s.freeList.back();
    s.freeList.pop_back();
    return f;
  }
  if (!gFiberPool.empty()) {
    FiberRec* f=gFiberPool.back();
    gFiberPool.pop_back();
    return f;
  }
  FiberRec* f=new FiberRec();
  const std::size_t page=static_cast<std::size_t>(sysconf(_SC_PAGESIZE));
  f->mapSize=StackBytes()+page;
  if (!UseGuardPages()) {
    if (gSlabFree.empty()) {
      const std::size_t bytes=f->mapSize*kSlab;
      void* p=mmap(nullptr,bytes,PROT_READ|PROT_WRITE,
                   MAP_PRIVATE|MAP_ANONYMOUS|MAP_NORESERVE,-1,0);
      if (p==MAP_FAILED) {
        delete f;
        throw std::runtime_error("Geospace::GPU::Fiber: cannot allocate fiber stack slab");
      }
      gSlabs.push_back(Slab{static_cast<char*>(p),bytes});
      for (int i=kSlab-1;i>=0;--i) gSlabFree.push_back(static_cast<char*>(p)+f->mapSize*i);
    }
    f->stackBase=gSlabFree.back();
    gSlabFree.pop_back();
    f->slab=true;
    return f;
  }
  void* p=mmap(nullptr,f->mapSize,PROT_READ|PROT_WRITE,
               MAP_PRIVATE|MAP_ANONYMOUS|MAP_NORESERVE,-1,0);
  if (p==MAP_FAILED) {
    delete f;
    throw std::runtime_error("Geospace::GPU::Fiber: cannot allocate fiber stack");
  }
  // Guard page at the low end: a stack overflow faults instead of corrupting memory.
  mprotect(p,page,PROT_NONE);
  f->stackBase=static_cast<char*>(p);
  f->slab=false;
  return f;
}

// Switch from the scheduler to fiber f, and from the running fiber f back.
void SwitchToFiber(Scheduler& s,FiberRec* f) {
#if GEOSPACE_FIBER_FAST_SWITCH
  geospace_fiber_switch(&s.mainSp,f->sp);
#else
  swapcontext(&s.mainCtx,&f->ctx);
#endif
}
void SwitchToScheduler(Scheduler& s,FiberRec* f) {
#if GEOSPACE_FIBER_FAST_SWITCH
  geospace_fiber_switch(&f->sp,s.mainSp);
#else
  swapcontext(&f->ctx,&s.mainCtx);
#endif
}

void Trampoline() {
  Scheduler& s=*tSched;
  FiberRec* f=s.current;
  try {
    f->body();
  }
  catch (...) {
    JoinGroup* g=f->group;
    if (!g->firstError) g->firstError=std::current_exception();
    g->stop=true;
  }
  f->state=FiberRec::Done;
  JoinGroup* g=f->group;
  if (--g->remaining==0 && g->waiter!=nullptr && g->waiterJoining) Wake(g->waiter);
  SwitchToScheduler(s,f);
  std::abort(); // a finished fiber is never resumed
}

// Prepare f so that the first switch to it enters Trampoline() on its own stack.
void InitContext(FiberRec* f,std::size_t page) {
#if GEOSPACE_FIBER_FAST_SWITCH
  // Frame consumed by geospace_fiber_switch: x87 control word, MXCSR, r15..r12, rbx,
  // rbp, return address (Trampoline), plus a dummy return address for Trampoline so
  // that its entry sees the ABI alignment (rsp % 16 == 8).
  std::uintptr_t top=reinterpret_cast<std::uintptr_t>(f->stackBase+f->mapSize);
  top&=~static_cast<std::uintptr_t>(15);
  std::uint64_t* sp=reinterpret_cast<std::uint64_t*>(top-80);
  std::uint16_t fcw;
  std::uint32_t mxcsr;
  __asm__ volatile("fnstcw %0" : "=m"(fcw));
  __asm__ volatile("stmxcsr %0" : "=m"(mxcsr));
  for (int i=0;i<10;++i) sp[i]=0;
  std::memcpy(&sp[0],&fcw,sizeof(fcw));
  std::memcpy(reinterpret_cast<char*>(sp)+8,&mxcsr,sizeof(mxcsr));
  sp[8]=reinterpret_cast<std::uint64_t>(&Trampoline);
  sp[9]=0;
  f->sp=sp;
  (void)page;
#else
  getcontext(&f->ctx);
  f->ctx.uc_stack.ss_sp=f->stackBase+page;
  f->ctx.uc_stack.ss_size=f->mapSize-page;
  f->ctx.uc_link=nullptr;
  makecontext(&f->ctx,Trampoline,0);
#endif
}

void StartFiber(Scheduler& s,JoinGroup* group,std::function<void()> body) {
  FiberRec* f=AllocateFiber(s);
  const std::size_t page=static_cast<std::size_t>(sysconf(_SC_PAGESIZE));
  f->body=std::move(body);
  f->group=group;
  f->state=FiberRec::Ready;
  InitContext(f,page);
  ++s.live;
  ++s.stats.fibersCreated;
  s.ready.push_back(f);
}

void Recycle(Scheduler& s,FiberRec* f) {
  f->body=nullptr;
  --s.live;
  s.freeList.push_back(f);
}

// Spawn worker fibers that share one index counter.
long long SpawnWorkers(Scheduler& s,JoinGroup* group,long long begin,long long end,
                       int nWorkers,
                       const std::function<void(long long,int)>& body,
                       std::shared_ptr<std::atomic<long long>> next) {
  for (int w=0;w<nWorkers;++w) {
    StartFiber(s,group,[group,end,w,&body,next]() {
      for (;;) {
        if (group->stop) break;
        const long long i=next->fetch_add(1);
        if (i>=end) break;
        body(i,w);
      }
    });
  }
  (void)begin;
  return nWorkers;
}

void RunScheduler(Scheduler& s,JoinGroup& group) {
  while (group.remaining>0) {
    if (!s.ready.empty()) {
      FiberRec* f=s.ready.front();
      s.ready.pop_front();
      f->state=FiberRec::Running;
      s.current=f;
      ++s.stats.switches;
      const auto tf=std::chrono::steady_clock::now();
      SwitchToFiber(s,f);
      s.stats.fiberSeconds+=std::chrono::duration<double>(std::chrono::steady_clock::now()-tf).count();
      s.current=nullptr;
      if (f->state==FiberRec::Done) Recycle(s,f);
      continue;
    }
    const auto ts=std::chrono::steady_clock::now();
    if (gServiceHook) gServiceHook();
    const auto ti=std::chrono::steady_clock::now();
    s.stats.serviceSeconds+=std::chrono::duration<double>(ti-ts).count();
    ++s.stats.idleCalls;
    const bool ok=(gIdleHook!=nullptr && gIdleHook());
    s.stats.idleSeconds+=std::chrono::duration<double>(std::chrono::steady_clock::now()-ti).count();
    if (!ok) {
      throw std::runtime_error(
          "Geospace::GPU::Fiber: all fibers are blocked and the trajectory engine has no "
          "pending work (deadlock)");
    }
  }
}

} // namespace

int gDefaultMaxFibers=4096;

int MaxFibers() {
  if (gMaxFibers<0) {
    gMaxFibers=gDefaultMaxFibers;
    if (const char* e=std::getenv("AMPS_EARTH_FIBERS")) {
      const int v=std::atoi(e);
      if (v>0) gMaxFibers=v;
    }
  }
  return gMaxFibers;
}

void SetDefaultMaxFibers(int n) {
  if (n>0) {
    gDefaultMaxFibers=n;
    gMaxFibers=-1;
  }
}

bool InFiber() {
  if (tSched==nullptr || tSched->current==nullptr) return false;
#ifdef _OPENMP
  if (omp_in_parallel()) return false;
#endif
  return true;
}

FiberRec* Current() { return (tSched!=nullptr) ? tSched->current : nullptr; }

void Block() {
  Scheduler& s=*tSched;
  FiberRec* f=s.current;
  f->state=FiberRec::Blocked;
  SwitchToScheduler(s,f);
}

void Wake(FiberRec* f) {
  f->state=FiberRec::Ready;
  tSched->ready.push_back(f);
}

void SetIdleHook(bool (*hook)()) { gIdleHook=hook; }
void SetServiceHook(std::function<void()> hook) { gServiceHook=std::move(hook); }

Stats LastStats() { return gLastStats; }

void ParallelFor(long long begin,long long end,int maxWorkers,
                 const std::function<void(long long,int)>& body) {
  if (end<=begin) return;
  const long long n=end-begin;

  if (InFiber()) {
    // Nested loop inside a running fiber.  The calling fiber is one of the workers: it
    // pulls iterations from the shared counter like its children, and before each
    // iteration it spawns additional child workers whenever fiber budget has become
    // available (other loops finishing), up to maxWorkers in total.  A loop started
    // while the budget is exhausted therefore still fans out later instead of running
    // all iterations one at a time.  The parent then waits for its children.
    Scheduler& s=*tSched;
    const long long maxChildren=std::min<long long>(n,std::max(1,maxWorkers))-1;
    JoinGroup group{0,s.current,nullptr,false};
    group.waiterJoining=false;
    auto next=std::make_shared<std::atomic<long long>>(begin);
    long long spawned=0;
    auto spawnAvailable=[&]() {
      while (spawned<maxChildren && next->load()<end &&
             s.live<static_cast<long long>(MaxFibers())) {
        ++group.remaining;
        const int w=static_cast<int>(++spawned);
        JoinGroup* g=&group;
        StartFiber(s,g,[g,end,w,&body,next]() {
          for (;;) {
            if (g->stop) break;
            const long long i=next->fetch_add(1);
            if (i>=end) break;
            body(i,w);
          }
        });
      }
    };
    std::exception_ptr parentError;
    try {
      for (;;) {
        spawnAvailable();
        if (group.stop) break;
        const long long i=next->fetch_add(1);
        if (i>=end) break;
        body(i,0);
      }
    }
    catch (...) {
      parentError=std::current_exception();
      group.stop=true;
    }
    if (group.remaining>0) {
      group.waiterJoining=true;
      Block();
    }
    if (parentError) std::rethrow_exception(parentError);
    if (group.firstError) std::rethrow_exception(group.firstError);
    return;
  }

  if (tSched!=nullptr) {
    // Called from scheduler context (hook) or from an OpenMP region inside a fiber:
    // fibers cannot be nested here, run sequentially with synchronous tracing.
    for (long long i=begin;i<end;++i) body(i,0);
    return;
  }

  Scheduler s;
  tSched=&s;
  long long nWorkers=std::min<long long>(n,std::max(1,maxWorkers));
  nWorkers=std::min<long long>(nWorkers,MaxFibers());
  JoinGroup group{nWorkers,nullptr,nullptr,false};
  auto next=std::make_shared<std::atomic<long long>>(begin);
  try {
    SpawnWorkers(s,&group,begin,end,static_cast<int>(nWorkers),body,next);
    RunScheduler(s,group);
  }
  catch (...) {
    tSched=nullptr;
    throw;
  }
  gLastStats=s.stats;
  gFiberPool.insert(gFiberPool.end(),s.freeList.begin(),s.freeList.end());
  s.freeList.clear();
  tSched=nullptr;
  if (group.firstError) std::rethrow_exception(group.firstError);
}

} // namespace Fiber
} // namespace GPU
} // namespace Geospace
