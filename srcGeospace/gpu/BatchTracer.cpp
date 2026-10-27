//======================================================================================
// BatchTracer.cpp -- see BatchTracer.h.
//======================================================================================

#include "BatchTracer.h"
#include "HostThreadPool.h"

#include <algorithm>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <stdexcept>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace Geospace {
namespace GPU {

namespace {

int EnvInt(const char* name,int fallback) {
  if (const char* e=std::getenv(name)) {
    const int v=std::atoi(e);
    if (v>0) return v;
  }
  return fallback;
}

bool IdleHook() { return BatchTracer::Get().Progress(); }

// Results finished and not yet collected (eg_dev_wait returns this count when no
// launch is pending).
int eg_dev_collectable_(EgDevicePool* p) {
  return eg_dev_pending_launch(p) ? 0 : eg_dev_wait(p);
}

} // namespace

BatchTracer& BatchTracer::Get() {
  static BatchTracer instance;
  return instance;
}

bool BatchTracer::Enabled() {
  static int enabled=-1;
  if (enabled<0) {
    const char* e=std::getenv("AMPS_EARTH_ENGINE");
    enabled=(e!=nullptr && (std::strcmp(e,"LEGACY")==0 || std::strcmp(e,"legacy")==0)) ? 0 : 1;
  }
  return enabled==1;
}

BatchTracer::BatchTracer() {
#ifdef _OPENMP
  hostThreads_=omp_get_max_threads();
#endif
  accelerator_=(eg_dev_is_accelerator()!=0);
  turnCapacity_=EnvInt("AMPS_EARTH_TURN_CAPACITY",512);
  {
    const char* e=std::getenv("AMPS_EARTH_ASYNC");
    asyncDevice_=!(e!=nullptr && std::string(e)=="0");
  }
  Fiber::SetIdleHook(&IdleHook);
  // Every fiber has at most one trajectory in flight, so the fiber budget must reach
  // the pool capacity to fill an accelerator pool (CPU pools are smaller than 4096).
  if (accelerator_) Fiber::SetDefaultMaxFibers(EnvInt("AMPS_EARTH_POOL_CAPACITY",65536));
}

void BatchTracer::SetHostThreads(int n) {
  hostThreads_=std::max(1,n);
  if (dev_!=nullptr) eg_dev_set_threads(dev_,hostThreads_);
  if (host_) host_->SetThreads(hostThreads_);
}

int BatchTracer::PreferredFibers() const {
  const int capacity=accelerator_ ? EnvInt("AMPS_EARTH_POOL_CAPACITY",65536)
                                  : EnvInt("AMPS_EARTH_POOL_CAPACITY",
                                           std::max(256,64*hostThreads_));
  return EnvInt("AMPS_EARTH_FIBERS",capacity);
}

std::string BatchTracer::BackendName() const {
  return std::string("device pool: ")+eg_dev_backend_name()+
         (accelerator_ ? " (accelerator)" : " (host)")+"; host pool: "+
         (HostRuntimeIsOpenMP() ? "CPU-OpenMP/" : "CPU-threads/")+EG_REAL_NAME;
}

void BatchTracer::EnsurePools() {
  if (dev_!=nullptr) return;
  capacity_=accelerator_ ? EnvInt("AMPS_EARTH_POOL_CAPACITY",65536)
                         : EnvInt("AMPS_EARTH_POOL_CAPACITY",std::max(256,64*hostThreads_));
  steps_=accelerator_ ? EnvInt("AMPS_EARTH_STEPS_PER_LAUNCH",2048)
                      : EnvInt("AMPS_EARTH_STEPS_PER_LAUNCH",256);
  dev_=eg_dev_create(capacity_,steps_,turnCapacity_);
  eg_dev_set_threads(dev_,hostThreads_);
  const int hostCapacity=EnvInt("AMPS_EARTH_HOST_POOL_CAPACITY",std::max(64,16*hostThreads_));
  host_.reset(new HostPool(hostCapacity,EnvInt("AMPS_EARTH_HOST_STEPS_PER_LAUNCH",64),
                           turnCapacity_,hostThreads_));
  if (Fiber::MaxFibers()<capacity_) {
    // Fiber budget is configured independently (AMPS_EARTH_FIBERS).
  }
}

int BatchTracer::RegisterContext(const EgContext& c) {
  std::lock_guard<std::mutex> lock(mutex_);
  for (std::size_t i=0;i<contexts_.size();++i)
    if (std::memcmp(&contexts_[i],&c,sizeof(EgContext))==0) return static_cast<int>(i);
  contexts_.push_back(c);
  return static_cast<int>(contexts_.size()-1);
}

void BatchTracer::SetMesh(const EgMeshHost& m,unsigned long long generation) {
  EnsurePools();
  if (!inflight_.empty() || !queue_.empty())
    throw std::logic_error("BatchTracer::SetMesh while trajectories are in flight");
  eg_dev_wait(dev_);
  hostMesh_.Set(m);
  if (eg_dev_set_mesh(dev_,&m)!=0)
    throw std::logic_error("BatchTracer::SetMesh: device pool busy");
  meshGeneration_=generation;
  haveMesh_=true;
}

void BatchTracer::ClearMesh() {
  if (dev_!=nullptr) eg_dev_set_mesh(dev_,nullptr);
  hostMesh_.Clear();
  haveMesh_=false;
}

bool BatchTracer::HaveMesh(unsigned long long generation) const {
  return haveMesh_ && meshGeneration_==generation;
}

void BatchTracer::Trace(const EgRequest& rqIn,const HostFieldRef& host,TraceOutput& out) {
  out.turns.clear();
  out.error.clear();
  if (!Fiber::InFiber()) {
    EgContext cx;
    {
      std::lock_guard<std::mutex> lock(mutex_);
      cx=contexts_.at(static_cast<std::size_t>(rqIn.context));
    }
    TraceSync(rqIn,cx,&hostMesh_,host,turnCapacity_,out.res,&out.turns,&out.error);
    out.res.tag=rqIn.tag;
    __atomic_add_fetch(&statSync_,1,__ATOMIC_RELAXED);
    return;
  }
  EnsurePools();
  if (!cancelledFibers_.empty() && cancelledFibers_.count(Fiber::Current())) throw TraceCancelled();
  Pending* p=new Pending();
  p->rq=rqIn;
  p->rq.tag=nextTag_++;
  p->host=host;
  p->out=&out;
  p->fiber=Fiber::Current();
  {
    std::lock_guard<std::mutex> lock(mutex_);
    p->device=(contexts_.at(static_cast<std::size_t>(rqIn.context)).fieldKind!=kFieldHost);
  }
  queue_.push_back(p);
  inflight_[p->rq.tag]=p;
  byFiber_[p->fiber]=p;
  ++statTraced_;
  Fiber::Block();
  byFiber_.erase(p->fiber);
  if (p->cancelled) {
    // Cancel() resumed this fiber; the pool may still return the result later, in
    // which case Deliver() deletes the record.
    p->fiberGone=true;
    if (p->resultArrived) delete p;
    throw TraceCancelled();
  }
  out.res.tag=rqIn.tag;
  delete p;
}

void BatchTracer::Cancel(Fiber::FiberRec* fiber) {
  if (fiber==nullptr) return;
  cancelledFibers_.insert(fiber);
  auto it=byFiber_.find(fiber);
  if (it==byFiber_.end()) return;           // not waiting: the next Trace() throws
  Pending* p=it->second;
  // Already cancelled, or the result was delivered and the fiber is about to resume
  // with it: nothing to cancel.
  if (p->cancelled || p->resultArrived) return;
  p->cancelled=true;
  ++statCancelled_;
  if (!p->submitted) {
    // Still queued: drop it.  Nothing will be delivered for it.
    for (auto q=queue_.begin();q!=queue_.end();++q) if (*q==p) { queue_.erase(q); break; }
    inflight_.erase(p->rq.tag);
    p->resultArrived=true;
  }
  else if (p->device) {
    eg_dev_cancel(dev_,p->rq.tag);          // stopped at the next launch
  }
  Fiber::Wake(fiber);
}

void BatchTracer::ClearCancel(Fiber::FiberRec* fiber) {
  if (!cancelledFibers_.empty()) cancelledFibers_.erase(fiber);
}

void BatchTracer::Deliver(const EgResult& res,std::vector<EgTurnState>* turns,
                          const std::string* err) {
  auto it=inflight_.find(res.tag);
  if (it==inflight_.end()) throw std::logic_error("BatchTracer: unknown result tag");
  Pending* p=it->second;
  inflight_.erase(it);
  if (p->cancelled) {
    // The fiber was already resumed by Cancel(); discard the result.
    p->resultArrived=true;
    if (p->fiberGone) delete p;
    return;
  }
  p->resultArrived=true;
  p->out->res=res;
  if (turns!=nullptr) p->out->turns.swap(*turns);
  if (err!=nullptr) p->out->error=*err;
  ++statDelivered_;
  Fiber::Wake(p->fiber);
}

bool BatchTracer::Progress() {
  EnsurePools();
  bool progressed=false;

  // 1. Upload new contexts to idle pools.
  int nCtx;
  {
    std::lock_guard<std::mutex> lock(mutex_);
    nCtx=static_cast<int>(contexts_.size());
    if (nCtx>devUploaded_ && eg_dev_active(dev_)==0) {
      if (eg_dev_set_contexts(dev_,contexts_.data(),nCtx)==0) devUploaded_=nCtx;
    }
    if (nCtx>hostUploaded_ && host_->Active()==0) {
      host_->SetContexts(contexts_);
      hostUploaded_=nCtx;
    }
  }

  // 2. Submit queued requests (FIFO) whose context is resident.  A request that needs
  //    a newer context stops further submissions to that pool until it drains.
  const auto tsub=std::chrono::steady_clock::now();
  bool devBlocked=false,hostBlocked=false;
  std::deque<Pending*> remaining;
  while (!queue_.empty()) {
    Pending* p=queue_.front();
    queue_.pop_front();
    bool submitted=false;
    if (p->device) {
      if (!devBlocked && p->rq.context<devUploaded_) {
        if (eg_dev_submit(dev_,&p->rq,1)==1) submitted=true;
      }
      else devBlocked=true;
    }
    else {
      if (!hostBlocked && p->rq.context<hostUploaded_) {
        if (host_->Submit(p->rq,p->host)) submitted=true;
      }
      else hostBlocked=true;
    }
    if (submitted) { p->submitted=true; progressed=true; }
    else remaining.push_back(p);
  }
  queue_.swap(remaining);
  statSubmitSeconds_+=std::chrono::duration<double>(std::chrono::steady_clock::now()-tsub).count();

  // 3. Advance pools and deliver finished trajectories.
  //
  // Device pool, asynchronous mode (default): wait for the running launch, start the
  // next one immediately (still-active slots plus requests submitted meanwhile), and
  // only then deliver the finished results.  The resumed fibers execute their search
  // logic and submit new requests while the next launch integrates on the
  // accelerator / worker threads.  Synchronous mode launches and waits in one call.
  const auto t0=std::chrono::steady_clock::now();
  int nDone=0;
  if (asyncDevice_) {
    if (eg_dev_pending_launch(dev_)) {
      const auto tw=std::chrono::steady_clock::now();
      nDone=eg_dev_wait(dev_);
      const auto tl=std::chrono::steady_clock::now();
      statWaitSeconds_+=std::chrono::duration<double>(tl-tw).count();
      statActiveSum_+=eg_dev_active(dev_);
      if (eg_dev_launch(dev_)>0) ++statLaunches_;
      statLaunchSeconds_+=std::chrono::duration<double>(std::chrono::steady_clock::now()-tl).count();
      progressed=true;
    }
    else if (eg_dev_active(dev_)-eg_dev_collectable_(dev_)>0) {
      if (eg_dev_launch(dev_)>0) ++statLaunches_;
      progressed=true;
    }
    if (nDone==0) nDone=eg_dev_collectable_(dev_);
  }
  else if (eg_dev_active(dev_)>0) {
    statActiveSum_+=eg_dev_active(dev_);
    nDone=eg_dev_advance(dev_);
    ++statLaunches_;
    progressed=true;
  }
  if (nDone>0) {
    std::vector<EgResult> res(static_cast<std::size_t>(nDone));
    const int stride=turnCapacity_;
    std::vector<EgTurnState> turns;
    bool anyTurns=false;
    {
      std::lock_guard<std::mutex> lock(mutex_);
      for (int i=0;i<devUploaded_;++i)
        if (contexts_[i].trap.enabled && contexts_[i].trap.driftEnabled &&
            contexts_[i].trap.recordDriftTurnStates) anyTurns=true;
    }
    if (anyTurns) turns.resize(static_cast<std::size_t>(nDone)*stride);
    const auto tc=std::chrono::steady_clock::now();
    const int n=eg_dev_collect(dev_,res.data(),nDone,anyTurns ? turns.data() : nullptr,
                               anyTurns ? stride : 0);
    const auto td=std::chrono::steady_clock::now();
    statCollectSeconds_+=std::chrono::duration<double>(td-tc).count();
    for (int k=0;k<n;++k) {
      std::vector<EgTurnState> t;
      if (anyTurns) {
        const int nt=std::min(res[k].turnStateCount,stride);
        t.assign(turns.begin()+(long long)k*stride,turns.begin()+(long long)k*stride+nt);
      }
      Deliver(res[k],&t,nullptr);
    }
    statDeliverSeconds_+=std::chrono::duration<double>(std::chrono::steady_clock::now()-td).count();
    if (n>0) progressed=true;
  }
  if (host_->Active()>0) {
    host_->Advance();
    ++statLaunches_;
    progressed=true;
    EgResult r;
    std::vector<EgTurnState> t;
    std::string err;
    while (host_->Collect(r,&t,&err)) Deliver(r,&t,&err);
  }
  statKernelSeconds_+=std::chrono::duration<double>(
      std::chrono::steady_clock::now()-t0).count();

  // 4. Pools drained but requests still waiting for a context upload: retry now.
  if (!progressed && !queue_.empty()) {
    std::lock_guard<std::mutex> lock(mutex_);
    if (eg_dev_active(dev_)==0 && static_cast<int>(contexts_.size())>devUploaded_) {
      if (eg_dev_set_contexts(dev_,contexts_.data(),static_cast<int>(contexts_.size()))==0) {
        devUploaded_=static_cast<int>(contexts_.size());
        progressed=true;
      }
    }
    if (host_->Active()==0 && static_cast<int>(contexts_.size())>hostUploaded_) {
      host_->SetContexts(contexts_);
      hostUploaded_=static_cast<int>(contexts_.size());
      progressed=true;
    }
  }
  return progressed;
}

void BatchTracer::PrintStatistics(const char* label) {
  const double avgActive=statLaunches_>0 ? double(statActiveSum_)/double(statLaunches_) : 0.0;
  std::printf("[earth-engine] %s: batched=%lld sync=%lld launches=%lld avgActive=%.1f "
              "kernel+collect=%.3fs backend=%s\n",
              label!=nullptr ? label : "",statTraced_,statSync_,statLaunches_,avgActive,
              statKernelSeconds_,BackendName().c_str());
  std::printf("[earth-engine] %s: progress split: submit=%.3fs wait=%.3fs launch=%.3fs "
              "collect=%.3fs deliver=%.3fs cancelled=%lld\n",label!=nullptr ? label : "",
              statSubmitSeconds_,statWaitSeconds_,statLaunchSeconds_,statCollectSeconds_,
              statDeliverSeconds_,statCancelled_);
  std::fflush(stdout);
}

} // namespace GPU
} // namespace Geospace
