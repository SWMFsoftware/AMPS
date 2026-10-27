//======================================================================================
// EngineHost.cpp -- host executor and synchronous tracer (see EngineHost.h).
//======================================================================================

#include "EngineHost.h"
#include "HostThreadPool.h"

#include <algorithm>
#include <type_traits>
#include <exception>
#include <stdexcept>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace Geospace {
namespace GPU {

namespace {

typedef Real R;
typedef Units U;

// Field functor for host-only evaluators.  Exceptions propagate to the executor.
class HostField {
 public:
  explicit HostField(const HostFieldRef& ref) : ref_(ref) { ResetDipoleErrorStats(stats_); }
  void ResetAvailability() {}
  bool Available() const { return true; }
  int Error() const { return kErrNone; }
  int LastNode() const { return -1; }
  int UnavailableKind() const { return 0; }
  double UnavailableX(int) const { return 0.0; }
  const EgDipoleErrorStats& Stats() const { return stats_; }
  void GetB(const Vec3<R>& x,Vec3<R>& B) {
    double xs[3],bs[3];
    xs[0]=FromLength<U>(x.x); xs[1]=FromLength<U>(x.y); xs[2]=FromLength<U>(x.z);
    ref_.fn(ref_.obj,xs,bs);
    if constexpr (U::kSI) { B.x=bs[0]; B.y=bs[1]; B.z=bs[2]; }
    else { B.x=R(ToField<U>(bs[0])); B.y=R(ToField<U>(bs[1])); B.z=R(ToField<U>(bs[2])); }
  }
 private:
  HostFieldRef ref_;
  EgDipoleErrorStats stats_;
};

template <class T>
void ConvertVec(std::vector<T>& dst,const double* src,std::size_t n,double scale) {
  dst.resize(n);
  for (std::size_t i=0;i<n;++i) {
    if constexpr (U::kSI) { (void)scale; dst[i]=static_cast<T>(src[i]); }
    else dst[i]=static_cast<T>(src[i]/scale);
  }
}

template <class T>
void CopyVec(std::vector<T>& dst,const T* src,std::size_t n) {
  dst.assign(src,src+n);
}

} // namespace

void HostMesh::Set(const EgMeshHost& m) {
  const std::size_t nn=static_cast<std::size_t>(m.nNodes);
#if !defined(GEOSPACE_PRECISION_FLOAT)
  {
    // Double/SI: the engine representation equals the host description, so the view
    // aliases the caller's arrays (they must stay valid until the next Set/Clear; the
    // Mode3D exporter keeps them resident and GlobalMagneticField keeps B until the
    // generation changes).  Saves one copy of the compact B array per rank.
    Clear();
    view_.nNodes=m.nNodes;
    view_.root=m.root;
    for (int d=0;d<3;++d) {
      view_.nCells[d]=m.nCells[d];
      view_.xGlobalMin[d]=m.xGlobalMin[d];
      view_.xGlobalMax[d]=m.xGlobalMax[d];
      view_.dxMaxRefinement[d]=m.dxMaxRefinement[d];
    }
    view_.xmin=m.xmin; view_.xmax=m.xmax; view_.ixMin=m.ixMin; view_.sizeIndex=m.sizeIndex;
    view_.up=m.up; view_.down=m.down; view_.level=m.level; view_.minNeibLevel=m.minNeibLevel;
    view_.maxNeibLevel=m.maxNeibLevel; view_.used=m.used; view_.tempId=m.tempId;
    view_.usedLeafBlocks=m.usedLeafBlocks;
    view_.interiorCells=m.interiorCells;
    view_.B=m.B;
    view_.presence=m.presence;
    valid_=true;
    return;
  }
#endif
  ConvertVec(xmin_,m.xmin,3*nn,U::kLength_m);
  ConvertVec(xmax_,m.xmax,3*nn,U::kLength_m);
  ConvertVec(B_,m.B,3*static_cast<std::size_t>(m.interiorCells),U::kField_T);
  CopyVec(ixMin_,m.ixMin,3*nn);
  CopyVec(size_,m.sizeIndex,nn);
  CopyVec(up_,m.up,nn);
  CopyVec(down_,m.down,8*nn);
  CopyVec(level_,m.level,nn);
  CopyVec(minL_,m.minNeibLevel,nn);
  CopyVec(maxL_,m.maxNeibLevel,nn);
  CopyVec(used_,m.used,nn);
  CopyVec(presence_,m.presence,static_cast<std::size_t>(m.interiorCells));
  CopyVec(tempId_,m.tempId,nn);
  view_.nNodes=m.nNodes;
  view_.root=m.root;
  for (int d=0;d<3;++d) {
    view_.nCells[d]=m.nCells[d];
    view_.xGlobalMin[d]=static_cast<R>(ToLength<U>(m.xGlobalMin[d]));
    view_.xGlobalMax[d]=static_cast<R>(ToLength<U>(m.xGlobalMax[d]));
    view_.dxMaxRefinement[d]=static_cast<R>(ToLength<U>(m.dxMaxRefinement[d]));
  }
  view_.xmin=xmin_.data(); view_.xmax=xmax_.data(); view_.ixMin=ixMin_.data();
  view_.sizeIndex=size_.data(); view_.up=up_.data(); view_.down=down_.data();
  view_.level=level_.data(); view_.minNeibLevel=minL_.data();
  view_.maxNeibLevel=maxL_.data(); view_.used=used_.data(); view_.tempId=tempId_.data();
  view_.usedLeafBlocks=m.usedLeafBlocks;
  view_.interiorCells=m.interiorCells;
  view_.B=B_.data();
  view_.presence=presence_.data();
  valid_=true;
}

void HostMesh::Clear() {
  valid_=false;
  view_=EgMeshView<R>{};
  xmin_.clear(); xmax_.clear(); B_.clear(); ixMin_.clear(); size_.clear(); up_.clear();
  down_.clear(); level_.clear(); minL_.clear(); maxL_.clear(); used_.clear();
  presence_.clear(); tempId_.clear();
}

void TraceSync(const EgRequest& rq,const EgContext& cx,const HostMesh* mesh,
               const HostFieldRef& hostField,int turnCapacity,EgResult& res,
               std::vector<EgTurnState>* turns,std::string* error) {
  EgSlot<R> s;
  s.done=-1;
  s.lastNode=-1;
  const EgTrapLayout layout=MakeTrapLayout(cx.trap,turnCapacity);
  std::vector<double> td(static_cast<std::size_t>(layout.doublesPerSlot));
  std::vector<int> ti(static_cast<std::size_t>(layout.intsPerSlot));
  std::vector<EgTurnState> tt(static_cast<std::size_t>(layout.turnCapacity));
  const int K=1<<30;
  if (cx.fieldKind==kFieldHost) {
    HostField field(hostField);
    try {
      RunSlot<R,U>(s,rq,cx,field,td.data(),ti.data(),tt.data(),layout.turnCapacity,K);
    }
    catch (const std::exception& e) {
      s.res.error=kErrHostField;
      if (error) *error=e.what();
      s.done=1;
    }
  }
  else {
    EgMeshView<R> view{};
    if (mesh!=nullptr && mesh->Valid()) view=mesh->View();
    EgRowStencil<R> row,fineRow;
    DeviceField<R,U> field(cx,&view,&row,&fineRow,-1);
    RunSlot<R,U>(s,rq,cx,field,td.data(),ti.data(),tt.data(),layout.turnCapacity,K);
  }
  res=s.res;
  if (turns!=nullptr) {
    turns->assign(tt.begin(),tt.begin()+std::min<int>(res.turnStateCount,layout.turnCapacity));
  }
}

HostPool::HostPool(int capacity,int stepsPerLaunch,int turnCapacity,int threads)
    : capacity_(capacity),K_(stepsPerLaunch),turnCapacity_(turnCapacity),threads_(threads) {
  slots_.resize(capacity);
  reqs_.resize(capacity);
  fields_.resize(capacity);
  errors_.resize(capacity);
  for (int i=capacity-1;i>=0;--i) free_.push_back(i);
}

HostPool::~HostPool() {}

void HostPool::SetContexts(const std::vector<EgContext>& ctx) {
  if (!active_.empty() || !finished_.empty())
    throw std::logic_error("HostPool::SetContexts called while trajectories are active");
  ctx_=ctx;
  TrapArenaStride(ctx_.data(),static_cast<int>(ctx_.size()),turnCapacity_,
                  strideD_,strideI_,strideT_);
  trapD_.assign(static_cast<std::size_t>(capacity_)*strideD_,0.0);
  trapI_.assign(static_cast<std::size_t>(capacity_)*strideI_,0);
  turns_.assign(static_cast<std::size_t>(capacity_)*strideT_,EgTurnState());
}

bool HostPool::Submit(const EgRequest& rq,const HostFieldRef& field) {
  if (free_.empty()) return false;
  const int si=free_.back();
  free_.pop_back();
  reqs_[si]=rq;
  fields_[si]=field;
  errors_[si].clear();
  slots_[si].done=-1;
  slots_[si].lastNode=-1;
  active_.push_back(si);
  return true;
}

int HostPool::Advance() {
  const int nActive=static_cast<int>(active_.size());
  if (nActive>0) {
    std::vector<int> doneFlag(nActive,0);
    const int nThreads=std::max(1,threads_);
    auto body=[&](int a) {
      const int si=active_[a];
      EgSlot<R>& s=slots_[si];
      const EgRequest& rq=reqs_[si];
      const EgContext& cx=ctx_[rq.context];
      HostField field(fields_[si]);
      try {
        RunSlot<R,U>(s,rq,cx,field,
                     strideD_>0 ? trapD_.data()+(long long)si*strideD_ : nullptr,
                     strideI_>0 ? trapI_.data()+(long long)si*strideI_ : nullptr,
                     strideT_>0 ? turns_.data()+(long long)si*strideT_ : nullptr,
                     strideT_,K_);
      }
      catch (const std::exception& e) {
        s.res.error=kErrHostField;
        errors_[si]=e.what();
        s.done=1;
      }
      catch (...) {
        s.res.error=kErrHostField;
        errors_[si]="unknown exception in host field evaluator";
        s.done=1;
      }
      doneFlag[a]=s.done;
    };
    if (HostRuntimeIsOpenMP()) {
#pragma omp parallel for schedule(dynamic,1) num_threads(nThreads)
      for (int a=0;a<nActive;++a) body(a);
    }
    else HostParallelFor(nActive,nThreads,body);
    std::vector<int> still;
    still.reserve(nActive);
    for (int a=0;a<nActive;++a) {
      if (doneFlag[a]==1) finished_.push_back(active_[a]);
      else still.push_back(active_[a]);
    }
    active_.swap(still);
  }
  return static_cast<int>(finished_.size());
}

bool HostPool::Collect(EgResult& res,std::vector<EgTurnState>* turns,std::string* error) {
  if (finished_.empty()) return false;
  const int si=finished_.back();
  finished_.pop_back();
  res=slots_[si].res;
  if (turns!=nullptr) {
    const int n=std::min(res.turnStateCount,strideT_);
    turns->assign(turns_.begin()+(long long)si*strideT_,
                  turns_.begin()+(long long)si*strideT_+n);
  }
  if (error!=nullptr) *error=errors_[si];
  free_.push_back(si);
  return true;
}

} // namespace GPU
} // namespace Geospace
