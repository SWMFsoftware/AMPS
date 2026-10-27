#ifndef _SRC_GEOSPACE_GPU_ENGINE_HOST_H_
#define _SRC_GEOSPACE_GPU_ENGINE_HOST_H_

//======================================================================================
// EngineHost.h -- host executor of the trajectory engine.
//
//  * HostPool: trajectory pool for field models that cannot run in a device kernel
//    (Geopack IGRF and the Tsyganenko T96/T01/T05/TA15/TA16 Fortran models, and the
//    Mode3D forced-analytic empirical path).  Slots are advanced by OpenMP threads; the
//    field is evaluated through the legacy evaluator object, so results are identical
//    to the legacy tracer.
//  * TraceSync: one trajectory on the calling thread, for any field kind.  Used when a
//    request is made outside a fiber (e.g. inside an OpenMP region or by legacy code
//    paths that were not restructured).
//======================================================================================

#include "EnginePool.h"

#include <memory>
#include <string>
#include <vector>

namespace Geospace {
namespace GPU {

// Evaluate B [T] at x [m] with a host field object.  May throw.
typedef void (*HostFieldEvalFn)(const void* obj,const double x[3],double B[3]);

struct HostFieldRef {
  HostFieldEvalFn fn;
  const void* obj;
};

// Host copy of the flattened mesh in engine units (for TraceSync).
class HostMesh {
 public:
  void Set(const EgMeshHost& m);
  void Clear();
  bool Valid() const { return valid_; }
  const EgMeshView<Real>& View() const { return view_; }
 private:
  bool valid_=false;
  EgMeshView<Real> view_{};
  std::vector<Real> xmin_,xmax_,B_;
  std::vector<int> ixMin_,size_,up_,down_,level_,minL_,maxL_;
  std::vector<unsigned char> used_,presence_;
  std::vector<long long> tempId_;
};

// Single trajectory on the calling thread.  `error` receives an exception message from
// a host field evaluator (res.error==kErrHostField).
void TraceSync(const EgRequest& rq,const EgContext& cx,const HostMesh* mesh,
               const HostFieldRef& hostField,int turnCapacity,EgResult& res,
               std::vector<EgTurnState>* turns,std::string* error);

class HostPool {
 public:
  HostPool(int capacity,int stepsPerLaunch,int turnCapacity,int threads);
  ~HostPool();

  void SetContexts(const std::vector<EgContext>& ctx);
  int Capacity() const { return capacity_; }
  int Active() const { return static_cast<int>(active_.size()+finished_.size()); }
  int Free() const { return static_cast<int>(free_.size()); }
  void SetThreads(int n) { threads_=n; }
  void SetStepsPerLaunch(int K) { K_=K; }

  // Accepts one request; returns false when no slot is free.
  bool Submit(const EgRequest& rq,const HostFieldRef& field);
  int Advance();
  // Collect one finished trajectory; returns false when none.
  bool Collect(EgResult& res,std::vector<EgTurnState>* turns,std::string* error);

 private:
  int capacity_,K_,turnCapacity_,threads_;
  std::vector<EgSlot<Real>> slots_;
  std::vector<EgRequest> reqs_;
  std::vector<HostFieldRef> fields_;
  std::vector<std::string> errors_;
  std::vector<double> trapD_;
  std::vector<int> trapI_;
  std::vector<EgTurnState> turns_;
  long long strideD_=0,strideI_=0;
  int strideT_=0;
  std::vector<EgContext> ctx_;
  std::vector<int> free_,active_,finished_;
};

} // namespace GPU
} // namespace Geospace

#endif
