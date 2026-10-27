//======================================================================================
// EngineDevice.cpp -- device trajectory pool (see EngineDeviceAPI.h).
//
// One source, two targets:
//   * OpenACC (nvc++ -acc=gpu / -acc=multicore): slot data live in device memory
//     (acc_malloc); the advance kernel is an `acc parallel loop` over active slots.
//   * OpenMP host (g++ -fopenmp): the same loop body runs as `omp parallel for` over
//     host memory; device buffers alias the host mirrors, so no copies are made.
//======================================================================================

#include "EngineDeviceAPI.h"
#include "EnginePool.h"
#if !EG_HAVE_OPENACC
#include "HostThreadPool.h"
#endif

#include <algorithm>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <unordered_map>
#include <utility>
#include <new>
#include <stdexcept>
#include <string>
#include <thread>
#include <type_traits>
#include <vector>

#if EG_HAVE_OPENACC
#include <openacc.h>
#endif

using namespace Geospace::GPU;

namespace {

typedef Geospace::GPU::Real R;
typedef Geospace::GPU::Units U;

template <class T>
struct DevBuf {
  std::vector<T> h;
  T* d=nullptr;
  std::size_t n=0;

  void release() {
#if EG_HAVE_OPENACC
    if (d!=nullptr) acc_free(d);
#endif
    // Host builds: d is h.data() or an aliased external array; neither is freed here.
    d=nullptr;
    h.clear();
    h.shrink_to_fit();
    n=0;
  }
  void alloc(std::size_t count) {
    release();
    n=count;
    h.assign(count,T());
#if EG_HAVE_OPENACC
    if (count>0) {
      d=static_cast<T*>(acc_malloc(count*sizeof(T)));
      if (d==nullptr) throw std::bad_alloc();
    }
#else
    d=h.data();
#endif
  }
  void toDevice(std::size_t first,std::size_t count) {
#if EG_HAVE_OPENACC
    if (count>0) acc_memcpy_to_device(d+first,h.data()+first,count*sizeof(T));
#else
    (void)first; (void)count;
#endif
  }
  void fromDevice(std::size_t first,std::size_t count) {
#if EG_HAVE_OPENACC
    if (count>0) acc_memcpy_from_device(h.data()+first,d+first,count*sizeof(T));
#else
    (void)first; (void)count;
#endif
  }
  void allToDevice() { toDevice(0,n); }
  // Host builds only: use the caller's array in place (read-only kernel input).
  void alias(const T* p,std::size_t count) {
    release();
    n=count;
    d=const_cast<T*>(p);
  }
};

struct MeshScalars {
  int nNodes;
  int root;
  int nCells[3];
  R xGlobalMin[3];
  R xGlobalMax[3];
  R dxMaxRefinement[3];
  long long usedLeafBlocks;
  long long interiorCells;
  int nGhost;          // stencil-cell cache: ghost entries per leaf (0: no cache)
};

} // namespace

struct EgDevicePool {
  int capacity=0;
  int K=1024;
  int turnCapacity=0;
  int threads=1;

  DevBuf<EgSlot<R>> slots;
  DevBuf<EgRequest> reqs;
  DevBuf<double> trapD;
  DevBuf<int> trapI;
  DevBuf<EgTurnState> turns;
  long long strideD=0,strideI=0;
  int strideT=0;

  DevBuf<EgContext> ctx;
  int nCtx=0;

  bool haveMesh=false;
  MeshScalars ms{};
  DevBuf<R> mxmin,mxmax,mB;
  DevBuf<int> mixMin,mSize,mUp,mDown,mLevel,mMinL,mMaxL;
  DevBuf<unsigned char> mUsed,mPresence;
  DevBuf<long long> mTempId;
  // Stencil-cell cache (EngineMesh.h BuildStencilCache), built at every set_mesh.
  DevBuf<int> mCacheSlot,mNeib;
  DevBuf<short> mGhostSlot;
  DevBuf<EgGhostRef> mGhost;

  std::vector<int> freeSlots;
  std::vector<int> activeSlots;
  DevBuf<int> activeList;
  DevBuf<int> doneFlags;
  DevBuf<int> gatherList;
  DevBuf<EgResult> packed;
  DevBuf<EgTurnState> packedTurns;
  std::vector<int> finished;   // slot indices finished, waiting for collection

  // Asynchronous launch state.  `inflight` holds the slots of the running launch;
  // slots submitted meanwhile wait in activeSlots for the next launch.
  std::vector<int> inflight;
  bool launchPending=false;

  // Batched submission (OpenACC): eg_dev_submit() only records the new requests on
  // the host; FlushPending() uploads them with one transfer per launch and a scatter
  // kernel writes them, and a reset slot, into their slots.  (Per-request transfers
  // cost one synchronous host-to-device copy each.)
  DevBuf<int> pendIdx;
  DevBuf<EgRequest> pendReq;
  DevBuf<EgSlot<R>> slotTemplate;   // value-initialized slot, done=-1, lastNode=-1
  int nPend=0;

  // Cancellation: tag -> slot of every request not yet collected, and the slots to stop
  // at the next launch (with the tag they had when the cancellation was requested).
  std::unordered_map<long long,int> tagSlot;
  std::vector<std::pair<int,long long>> cancels;
  DevBuf<int> cancelIdx;
};

namespace {

// Advance one active slot by at most K loop iterations.  This is the complete kernel
// body; the OpenACC, OpenMP and thread-pool launchers below differ only in how the
// iterations a=0..nActive-1 are distributed.
EG_HD void AdvanceOne(int a,const int* active,EgSlot<R>* slots,const EgRequest* reqs,
                      const EgContext* ctx,double* td,int* ti,EgTurnState* tt,int* doneFlags,
                      long long sD,long long sI,int sT,int K,const MeshScalars& ms,
                      const R* mxmin,const R* mxmax,const R* mB,const int* mixMin,
                      const int* mSize,const int* mUp,const int* mDown,const int* mLevel,
                      const int* mMinL,const int* mMaxL,const unsigned char* mUsed,
                      const unsigned char* mPresence,const long long* mTempId,
                      const int* mCacheSlot,const short* mGhostSlot,const EgGhostRef* mGhost,
                      const int* mNeib) {
  const int si=active[a];
  EgSlot<R>& s=slots[si];
  const EgRequest& rq=reqs[si];
  const EgContext& cx=ctx[rq.context];

  EgMeshView<R> mv;
  mv.nNodes=ms.nNodes; mv.root=ms.root;
  for (int d=0;d<3;++d) {
    mv.nCells[d]=ms.nCells[d];
    mv.xGlobalMin[d]=ms.xGlobalMin[d];
    mv.xGlobalMax[d]=ms.xGlobalMax[d];
    mv.dxMaxRefinement[d]=ms.dxMaxRefinement[d];
  }
  mv.xmin=mxmin; mv.xmax=mxmax; mv.ixMin=mixMin; mv.sizeIndex=mSize; mv.up=mUp;
  mv.down=mDown; mv.level=mLevel; mv.minNeibLevel=mMinL; mv.maxNeibLevel=mMaxL;
  mv.used=mUsed; mv.tempId=mTempId; mv.usedLeafBlocks=ms.usedLeafBlocks;
  mv.interiorCells=ms.interiorCells; mv.B=mB; mv.presence=mPresence;
  mv.nGhost=ms.nGhost;
  mv.cacheSlot=ms.nGhost>0 ? mCacheSlot : nullptr;
  mv.ghostSlot=ms.nGhost>0 ? mGhostSlot : nullptr;
  mv.ghost=ms.nGhost>0 ? mGhost : nullptr;
  mv.neib=ms.nGhost>0 ? mNeib : nullptr;

  // Row stencils live in thread-local storage of the kernel (GPU local memory).
  EgRowStencil<R> row,fineRow;
  DeviceField<R,U> field(cx,&mv,&row,&fineRow,s.lastNode);
  RunSlot<R,U>(s,rq,cx,field,
               td!=nullptr ? td+(long long)si*sD : nullptr,
               ti!=nullptr ? ti+(long long)si*sI : nullptr,
               tt!=nullptr ? tt+(long long)si*sT : nullptr,
               sT,K);
  doneFlags[a]=s.done;
}

// Launch one advance over p->activeList[0..nActive).  With async=true the call
// returns immediately (OpenACC queue 1 / host worker threads) and WaitKernel() must
// be called before the results are read.
void LaunchKernel(EgDevicePool* p,int nActive,bool async) {
  const int* active=p->activeList.d;
  EgSlot<R>* slots=p->slots.d;
  const EgRequest* reqs=p->reqs.d;
  const EgContext* ctx=p->ctx.d;
  double* td=p->trapD.d;
  int* ti=p->trapI.d;
  EgTurnState* tt=p->turns.d;
  int* doneFlags=p->doneFlags.d;
  const long long sD=p->strideD,sI=p->strideI;
  const int sT=p->strideT;
  const int K=p->K;
  const MeshScalars ms=p->ms;
  const R* mxmin=p->mxmin.d;
  const R* mxmax=p->mxmax.d;
  const R* mB=p->mB.d;
  const int* mixMin=p->mixMin.d;
  const int* mSize=p->mSize.d;
  const int* mUp=p->mUp.d;
  const int* mDown=p->mDown.d;
  const int* mLevel=p->mLevel.d;
  const int* mMinL=p->mMinL.d;
  const int* mMaxL=p->mMaxL.d;
  const unsigned char* mUsed=p->mUsed.d;
  const unsigned char* mPresence=p->mPresence.d;
  const long long* mTempId=p->mTempId.d;
  const int* mCacheSlot=p->mCacheSlot.d;
  const short* mGhostSlot=p->mGhostSlot.d;
  const EgGhostRef* mGhost=p->mGhost.d;
  const int* mNeib=p->mNeib.d;
  const int nThreads=p->threads>0 ? p->threads : 1;

#if EG_HAVE_OPENACC
  (void)nThreads;
  if (async) {
#pragma acc parallel loop gang vector vector_length(64) async(1) \
    deviceptr(active,slots,reqs,ctx,td,ti,tt,doneFlags,mxmin,mxmax,mB,mixMin,mSize,mUp, \
              mDown,mLevel,mMinL,mMaxL,mUsed,mPresence,mTempId,mCacheSlot,mGhostSlot, \
              mGhost,mNeib) firstprivate(ms)
    for (int a=0;a<nActive;++a)
      AdvanceOne(a,active,slots,reqs,ctx,td,ti,tt,doneFlags,sD,sI,sT,K,ms,mxmin,mxmax,mB,
                 mixMin,mSize,mUp,mDown,mLevel,mMinL,mMaxL,mUsed,mPresence,mTempId,
                 mCacheSlot,mGhostSlot,mGhost,mNeib);
  }
  else {
#pragma acc parallel loop gang vector vector_length(64) \
    deviceptr(active,slots,reqs,ctx,td,ti,tt,doneFlags,mxmin,mxmax,mB,mixMin,mSize,mUp, \
              mDown,mLevel,mMinL,mMaxL,mUsed,mPresence,mTempId,mCacheSlot,mGhostSlot, \
              mGhost,mNeib) firstprivate(ms)
    for (int a=0;a<nActive;++a)
      AdvanceOne(a,active,slots,reqs,ctx,td,ti,tt,doneFlags,sD,sI,sT,K,ms,mxmin,mxmax,mB,
                 mixMin,mSize,mUp,mDown,mLevel,mMinL,mMaxL,mUsed,mPresence,mTempId,
                 mCacheSlot,mGhostSlot,mGhost,mNeib);
  }
#else
  auto body=[=](int a) {
    AdvanceOne(a,active,slots,reqs,ctx,td,ti,tt,doneFlags,sD,sI,sT,K,ms,mxmin,mxmax,mB,
               mixMin,mSize,mUp,mDown,mLevel,mMinL,mMaxL,mUsed,mPresence,mTempId,
                 mCacheSlot,mGhostSlot,mGhost,mNeib);
  };
  if (HostRuntimeIsOpenMP()) {
    // OpenMP launches are synchronous; WaitKernel() is then a no-op.
#pragma omp parallel for schedule(dynamic,1) num_threads(nThreads)
    for (int a=0;a<nActive;++a) body(a);
  }
  else if (async) HostParallelForAsync(nActive,nThreads,body);
  else HostParallelFor(nActive,nThreads,body);
#endif
}

void WaitKernel() {
#if EG_HAVE_OPENACC
  acc_wait(1);
#else
  HostParallelWait();
#endif
}

// Write the requests recorded by eg_dev_submit() since the last launch into their
// slots, with the slot reset to slotTemplate (what the host copy of a slot holds).
void FlushPending(EgDevicePool* p) {
#if EG_HAVE_OPENACC
  const int n=p->nPend;
  if (n==0) return;
  p->pendIdx.toDevice(0,n);
  p->pendReq.toDevice(0,n);
  const int* idx=p->pendIdx.d;
  const EgRequest* src=p->pendReq.d;
  const EgSlot<R>* tmpl=p->slotTemplate.d;
  EgSlot<R>* slots=p->slots.d;
  EgRequest* reqs=p->reqs.d;
#pragma acc parallel loop async(1) deviceptr(idx,src,tmpl,slots,reqs)
  for (int k=0;k<n;++k) {
    slots[idx[k]]=tmpl[0];
    reqs[idx[k]]=src[k];
  }
  p->nPend=0;
#else
  (void)p;
#endif
}

// Stop the cancelled slots before the next advance: done=1, termination kTermCancelled,
// result tag = request tag.  A slot whose request changed since the cancellation was
// requested (collected and reused) is left alone.
void ApplyCancels(EgDevicePool* p) {
  if (p->cancels.empty()) return;
  int n=0;
  for (const auto& c:p->cancels) {
    const int si=c.first;
    if (p->reqs.h[si].tag!=c.second) continue;
    if (p->tagSlot.find(c.second)==p->tagSlot.end()) continue;
    p->cancelIdx.h[n++]=si;
  }
  p->cancels.clear();
  if (n==0) return;
  EgSlot<R>* slots=p->slots.d;
  const EgRequest* reqs=p->reqs.d;
#if EG_HAVE_OPENACC
  p->cancelIdx.toDevice(0,n);
  const int* idx=p->cancelIdx.d;
#pragma acc parallel loop async(1) deviceptr(idx,slots,reqs)
#else
  const int* idx=p->cancelIdx.h.data();
#endif
  for (int k=0;k<n;++k) {
    EgSlot<R>& s=slots[idx[k]];
    s.done=1;
    s.res.tag=reqs[idx[k]].tag;
    s.res.termination=kTermCancelled;
  }
}

void GatherKernel(EgDevicePool* p,int nDone) {
  const int* list=p->gatherList.d;
  const EgSlot<R>* slots=p->slots.d;
  EgResult* out=p->packed.d;
  const EgTurnState* tt=p->turns.d;
  EgTurnState* outT=p->packedTurns.d;
  const int sT=p->strideT;
#if EG_HAVE_OPENACC
#pragma acc parallel loop deviceptr(list,slots,out,tt,outT)
#endif
  for (int k=0;k<nDone;++k) {
    const int si=list[k];
    out[k]=slots[si].res;
    if (sT>0 && tt!=nullptr && outT!=nullptr) {
      const int n=slots[si].res.turnStateCount;
      for (int j=0;j<n && j<sT;++j) outT[(long long)k*sT+j]=tt[(long long)si*sT+j];
    }
  }
}

void AllocateArenas(EgDevicePool* p) {
  TrapArenaStride(p->ctx.h.data(),p->nCtx,p->turnCapacity,p->strideD,p->strideI,p->strideT);
  const std::size_t cap=static_cast<std::size_t>(p->capacity);
  if (p->trapD.n!=cap*static_cast<std::size_t>(p->strideD)) p->trapD.alloc(cap*p->strideD);
  if (p->trapI.n!=cap*static_cast<std::size_t>(p->strideI)) p->trapI.alloc(cap*p->strideI);
  if (p->turns.n!=cap*static_cast<std::size_t>(p->strideT)) p->turns.alloc(cap*p->strideT);
  if (p->strideT>0) p->packedTurns.alloc(cap*static_cast<std::size_t>(p->strideT));
  else p->packedTurns.release();
}

template <class S,class T>
void ConvertCopy(DevBuf<T>& dst,const S* src,std::size_t n,double scale) {
  dst.alloc(n);
  (void)scale;
  for (std::size_t i=0;i<n;++i) {
    if constexpr (U::kSI) dst.h[i]=static_cast<T>(src[i]);
    else dst.h[i]=static_cast<T>(src[i]/scale);
  }
  dst.allToDevice();
}

template <class T>
void PlainCopy(DevBuf<T>& dst,const T* src,std::size_t n) {
  dst.alloc(n);
  if (n>0) std::memcpy(dst.h.data(),src,n*sizeof(T));
  dst.allToDevice();
}

// Host address of a mesh buffer: the host mirror, or the aliased caller array.
template <class T>
const T* HostPtr(const DevBuf<T>& b) { return b.h.empty() ? b.d : b.h.data(); }

// Take over a host vector as the buffer contents (device builds: upload, then drop the
// host copy).
template <class T>
void Adopt(DevBuf<T>& dst,std::vector<T>& src) {
  dst.release();
  dst.n=src.size();
  dst.h.swap(src);
#if EG_HAVE_OPENACC
  if (dst.n>0) {
    dst.d=static_cast<T*>(acc_malloc(dst.n*sizeof(T)));
    if (dst.d==nullptr) throw std::bad_alloc();
    dst.allToDevice();
  }
  dst.h.clear();
  dst.h.shrink_to_fit();
#else
  dst.d=dst.h.data();
#endif
}

void ReleaseMeshCache(EgDevicePool* p) {
  p->mCacheSlot.release(); p->mNeib.release(); p->mGhostSlot.release(); p->mGhost.release();
  p->ms.nGhost=0;
}

// Stencil-cell cache of the current mesh (EngineMesh.h BuildStencilCache), computed on
// the host from the same engine-unit arrays the kernels read.  AMPS_EARTH_MESH_CACHE=0
// disables it (uncached path, identical results).
void BuildMeshCache(EgDevicePool* p) {
  ReleaseMeshCache(p);
  const char* env=std::getenv("AMPS_EARTH_MESH_CACHE");
  if (env!=nullptr && std::string(env)=="0") return;
  const auto t0=std::chrono::steady_clock::now();
  EgMeshView<R> v{};
  v.nNodes=p->ms.nNodes; v.root=p->ms.root;
  for (int d=0;d<3;++d) {
    v.nCells[d]=p->ms.nCells[d];
    v.xGlobalMin[d]=p->ms.xGlobalMin[d];
    v.xGlobalMax[d]=p->ms.xGlobalMax[d];
    v.dxMaxRefinement[d]=p->ms.dxMaxRefinement[d];
  }
  v.xmin=HostPtr(p->mxmin); v.xmax=HostPtr(p->mxmax); v.ixMin=HostPtr(p->mixMin);
  v.sizeIndex=HostPtr(p->mSize); v.up=HostPtr(p->mUp); v.down=HostPtr(p->mDown);
  v.level=HostPtr(p->mLevel); v.minNeibLevel=HostPtr(p->mMinL);
  v.maxNeibLevel=HostPtr(p->mMaxL); v.used=HostPtr(p->mUsed); v.tempId=HostPtr(p->mTempId);
  v.usedLeafBlocks=p->ms.usedLeafBlocks; v.interiorCells=p->ms.interiorCells;
  v.B=nullptr; v.presence=nullptr;   // not used by the cache
  const int hw=static_cast<int>(std::thread::hardware_concurrency());
  EgStencilCacheHost c;
  BuildStencilCache(v,c,std::max(1,std::min(hw,16)));
  if (c.nGhost<=0 || c.nSlots<=0) return;
  const double mib=(c.cacheSlot.size()*sizeof(int)+c.ghostSlot.size()*sizeof(short)+
                    c.ghost.size()*sizeof(EgGhostRef)+c.neib.size()*sizeof(int))/1048576.0;
  const long long nSlots=c.nSlots,nRejected=c.nRejected;
  const int nGhost=c.nGhost;
  Adopt(p->mCacheSlot,c.cacheSlot);
  Adopt(p->mGhostSlot,c.ghostSlot);
  Adopt(p->mGhost,c.ghost);
  Adopt(p->mNeib,c.neib);
  p->ms.nGhost=nGhost;
  const char* verbose=std::getenv("AMPS_EARTH_ENGINE_VERBOSE");
  if (verbose!=nullptr && std::string(verbose)=="1") {
    const double dt=std::chrono::duration<double>(std::chrono::steady_clock::now()-t0).count();
    std::printf("[earth-engine] stencil-cell cache: leaves=%lld uncached=%lld ghostCells=%d "
                "memory=%.1f MiB (%.3fs)\n",nSlots,nRejected,nGhost,mib,dt);
    std::fflush(stdout);
  }
}

} // namespace

extern "C" {

const char* eg_dev_backend_name() {
#if EG_HAVE_OPENACC
  return "OpenACC/" EG_REAL_NAME;
#else
  return HostRuntimeIsOpenMP() ? "CPU-OpenMP/" EG_REAL_NAME : "CPU-threads/" EG_REAL_NAME;
#endif
}

int eg_dev_is_accelerator() {
#if EG_HAVE_OPENACC
  return acc_get_device_type()!=acc_device_host ? 1 : 0;
#else
  return 0;
#endif
}

EgDevicePool* eg_dev_create(int capacity,int stepsPerLaunch,int turnCapacity) {
  EgDevicePool* p=new EgDevicePool();
  p->capacity=capacity;
  p->K=stepsPerLaunch;
  p->turnCapacity=turnCapacity;
  p->slots.alloc(capacity);
  p->reqs.alloc(capacity);
  p->activeList.alloc(capacity);
  p->doneFlags.alloc(capacity);
  p->gatherList.alloc(capacity);
  p->packed.alloc(capacity);
  p->cancelIdx.alloc(capacity);
#if EG_HAVE_OPENACC
  p->pendIdx.alloc(capacity);
  p->pendReq.alloc(capacity);
  p->slotTemplate.alloc(1);
  p->slotTemplate.h[0].done=-1;
  p->slotTemplate.h[0].lastNode=-1;
  p->slotTemplate.allToDevice();
#endif
  p->freeSlots.reserve(capacity);
  for (int i=capacity-1;i>=0;--i) p->freeSlots.push_back(i);
  return p;
}

void eg_dev_destroy(EgDevicePool* p) { delete p; }

int eg_dev_set_contexts(EgDevicePool* p,const EgContext* c,int n) {
  if (!p->activeSlots.empty() || !p->finished.empty() || p->launchPending) return -1;
  p->ctx.alloc(n);
  if (n>0) std::memcpy(p->ctx.h.data(),c,n*sizeof(EgContext));
  p->ctx.allToDevice();
  p->nCtx=n;
  AllocateArenas(p);
  return 0;
}

int eg_dev_set_mesh(EgDevicePool* p,const EgMeshHost* m) {
  if (!p->activeSlots.empty() || !p->finished.empty() || p->launchPending) return -1;
  if (m==nullptr) {
    p->haveMesh=false;
    p->ms=MeshScalars{};
    p->mxmin.release(); p->mxmax.release(); p->mB.release();
    p->mixMin.release(); p->mSize.release(); p->mUp.release(); p->mDown.release();
    p->mLevel.release(); p->mMinL.release(); p->mMaxL.release();
    p->mUsed.release(); p->mPresence.release(); p->mTempId.release();
    ReleaseMeshCache(p);
    return 0;
  }
  const std::size_t nn=static_cast<std::size_t>(m->nNodes);
  p->ms.nNodes=m->nNodes;
  p->ms.root=m->root;
  for (int d=0;d<3;++d) {
    p->ms.nCells[d]=m->nCells[d];
    p->ms.xGlobalMin[d]=static_cast<R>(ToLength<U>(m->xGlobalMin[d]));
    p->ms.xGlobalMax[d]=static_cast<R>(ToLength<U>(m->xGlobalMax[d]));
    p->ms.dxMaxRefinement[d]=static_cast<R>(ToLength<U>(m->dxMaxRefinement[d]));
  }
  p->ms.usedLeafBlocks=m->usedLeafBlocks;
  p->ms.interiorCells=m->interiorCells;
#if !EG_HAVE_OPENACC && !defined(GEOSPACE_PRECISION_FLOAT)
  {
    // CPU double/SI build: the kernels run in host memory and read the caller's arrays
    // directly (BatchTracer::SetMesh contract: valid until the next set_mesh), so the
    // compact B array is not duplicated.
    p->mxmin.alias(m->xmin,3*nn); p->mxmax.alias(m->xmax,3*nn);
    p->mB.alias(m->B,3*static_cast<std::size_t>(m->interiorCells));
    p->mixMin.alias(m->ixMin,3*nn); p->mSize.alias(m->sizeIndex,nn); p->mUp.alias(m->up,nn);
    p->mDown.alias(m->down,8*nn); p->mLevel.alias(m->level,nn);
    p->mMinL.alias(m->minNeibLevel,nn); p->mMaxL.alias(m->maxNeibLevel,nn);
    p->mUsed.alias(m->used,nn);
    p->mPresence.alias(m->presence,static_cast<std::size_t>(m->interiorCells));
    p->mTempId.alias(m->tempId,nn);
    BuildMeshCache(p);
    p->haveMesh=true;
    return 0;
  }
#endif
  ConvertCopy(p->mxmin,m->xmin,3*nn,U::kLength_m);
  ConvertCopy(p->mxmax,m->xmax,3*nn,U::kLength_m);
  ConvertCopy(p->mB,m->B,3*static_cast<std::size_t>(m->interiorCells),U::kField_T);
  PlainCopy(p->mixMin,m->ixMin,3*nn);
  PlainCopy(p->mSize,m->sizeIndex,nn);
  PlainCopy(p->mUp,m->up,nn);
  PlainCopy(p->mDown,m->down,8*nn);
  PlainCopy(p->mLevel,m->level,nn);
  PlainCopy(p->mMinL,m->minNeibLevel,nn);
  PlainCopy(p->mMaxL,m->maxNeibLevel,nn);
  PlainCopy(p->mUsed,m->used,nn);
  PlainCopy(p->mPresence,m->presence,static_cast<std::size_t>(m->interiorCells));
  PlainCopy(p->mTempId,m->tempId,nn);
  BuildMeshCache(p);
  p->haveMesh=true;
  return 0;
}

int eg_dev_capacity(EgDevicePool* p) { return p->capacity; }
int eg_dev_active(EgDevicePool* p) {
  return static_cast<int>(p->activeSlots.size()+p->finished.size()+p->inflight.size());
}
int eg_dev_free(EgDevicePool* p) { return static_cast<int>(p->freeSlots.size()); }
void eg_dev_set_steps_per_launch(EgDevicePool* p,int K) { p->K=K; }

void eg_dev_cancel(EgDevicePool* p,long long tag) {
  auto it=p->tagSlot.find(tag);
  if (it==p->tagSlot.end()) return;
  p->cancels.emplace_back(it->second,tag);
}
void eg_dev_set_threads(EgDevicePool* p,int nThreads) { p->threads=nThreads; }

int eg_dev_submit(EgDevicePool* p,const EgRequest* reqs,int n) {
  int accepted=0;
  while (accepted<n && !p->freeSlots.empty()) {
    const int si=p->freeSlots.back();
    p->freeSlots.pop_back();
    p->reqs.h[si]=reqs[accepted];
    EgSlot<R>& s=p->slots.h[si];
    s.done=-1;
    s.lastNode=-1;
#if EG_HAVE_OPENACC
    p->pendIdx.h[p->nPend]=si;
    p->pendReq.h[p->nPend]=reqs[accepted];
    ++p->nPend;
#else
    p->reqs.toDevice(si,1);
    p->slots.toDevice(si,1);
#endif
    p->tagSlot[reqs[accepted].tag]=si;
    p->activeSlots.push_back(si);
    ++accepted;
  }
  return accepted;
}

int eg_dev_launch(EgDevicePool* p) {
  if (p->launchPending) return -1;
  const int nActive=static_cast<int>(p->activeSlots.size());
  if (nActive==0) return 0;
  p->inflight.swap(p->activeSlots);
  p->activeSlots.clear();
  for (int a=0;a<nActive;++a) p->activeList.h[a]=p->inflight[a];
  p->activeList.toDevice(0,nActive);
  FlushPending(p);
  ApplyCancels(p);
  LaunchKernel(p,nActive,true);
  p->launchPending=true;
  return nActive;
}

int eg_dev_wait(EgDevicePool* p) {
  if (!p->launchPending) return static_cast<int>(p->finished.size());
  WaitKernel();
  p->launchPending=false;
  const int nActive=static_cast<int>(p->inflight.size());
  p->doneFlags.fromDevice(0,nActive);
  for (int a=0;a<nActive;++a) {
    if (p->doneFlags.h[a]==1) p->finished.push_back(p->inflight[a]);
    else p->activeSlots.push_back(p->inflight[a]);
  }
  p->inflight.clear();
  return static_cast<int>(p->finished.size());
}

int eg_dev_advance(EgDevicePool* p) {
  if (p->launchPending) eg_dev_wait(p);
  const int nActive=static_cast<int>(p->activeSlots.size());
  if (nActive>0) {
    for (int a=0;a<nActive;++a) p->activeList.h[a]=p->activeSlots[a];
    p->activeList.toDevice(0,nActive);
    FlushPending(p);
    ApplyCancels(p);
    WaitKernel();   // the scatters run on queue 1, the synchronous kernel on the default queue
    LaunchKernel(p,nActive,false);
    p->doneFlags.fromDevice(0,nActive);
    std::vector<int> still;
    still.reserve(nActive);
    for (int a=0;a<nActive;++a) {
      if (p->doneFlags.h[a]==1) p->finished.push_back(p->activeSlots[a]);
      else still.push_back(p->activeSlots[a]);
    }
    p->activeSlots.swap(still);
  }
  return static_cast<int>(p->finished.size());
}

int eg_dev_pending_launch(EgDevicePool* p) { return p->launchPending ? 1 : 0; }

int eg_dev_collect(EgDevicePool* p,EgResult* res,int maxOut,EgTurnState* turns,
                   int turnStride) {
  const int nDone=std::min<int>(maxOut,static_cast<int>(p->finished.size()));
  if (nDone<=0) return 0;
  for (int k=0;k<nDone;++k) p->gatherList.h[k]=p->finished[k];
  p->gatherList.toDevice(0,nDone);
  GatherKernel(p,nDone);
  p->packed.fromDevice(0,nDone);
  if (p->strideT>0) p->packedTurns.fromDevice(0,static_cast<std::size_t>(nDone)*p->strideT);
  for (int k=0;k<nDone;++k) {
    res[k]=p->packed.h[k];
    if (turns!=nullptr && p->strideT>0 && turnStride>0) {
      const int nt=std::min(res[k].turnStateCount,std::min(turnStride,p->strideT));
      for (int j=0;j<nt;++j)
        turns[(long long)k*turnStride+j]=p->packedTurns.h[(long long)k*p->strideT+j];
    }
    p->tagSlot.erase(res[k].tag);
    p->freeSlots.push_back(p->finished[k]);
  }
  p->finished.erase(p->finished.begin(),p->finished.begin()+nDone);
  return nDone;
}

} // extern "C"
