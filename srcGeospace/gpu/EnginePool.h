#ifndef _SRC_GEOSPACE_GPU_ENGINE_POOL_H_
#define _SRC_GEOSPACE_GPU_ENGINE_POOL_H_

//======================================================================================
// EnginePool.h -- trajectory pool shared by the device and host executors.
//
// A pool owns `capacity` trajectory slots.  The host assigns requests to free slots,
// launches one kernel that advances every active slot by at most K loop iterations
// (each with its own adaptive time step), collects finished slots and refills them.
// Long and short trajectories therefore share the hardware without a global time step
// and without waiting for the slowest trajectory of a batch.
//======================================================================================

#include "EngineStep.h"

namespace Geospace {
namespace GPU {

// Host-side description of the flattened Mode3D mesh (SI, double).  The pool converts
// it to the engine Real/unit system and uploads it.
struct EgMeshHost {
  int nNodes;
  int root;
  int nCells[3];
  double xGlobalMin[3];
  double xGlobalMax[3];
  double dxMaxRefinement[3];
  const double* xmin;
  const double* xmax;
  const int* ixMin;
  const int* sizeIndex;
  const int* up;
  const int* down;
  const int* level;
  const int* minNeibLevel;
  const int* maxNeibLevel;
  const unsigned char* used;
  const long long* tempId;
  long long usedLeafBlocks;
  long long interiorCells;
  const double* B;
  const unsigned char* presence;
};

// Run one slot: TraceBegin when new, then up to K loop iterations.
template <class R,class U,class F>
EG_HD void RunSlot(EgSlot<R>& s,const EgRequest& rq,const EgContext& cx,F& field,
                   double* td,int* ti,EgTurnState* tt,int turnCap,int K) {
  const EgTrapLayout layout=MakeTrapLayout(cx.trap,turnCap);
  EgTrapDetector trap(cx.trap,cx.boxMin_m,cx.boxMax_m,layout,s.trap,
                      layout.doublesPerSlot>0 ? td : nullptr,
                      layout.intsPerSlot>0 ? ti : nullptr,
                      layout.turnCapacity>0 ? tt : nullptr);
  bool finished=false;
  if (s.done<0) finished=Tracer<R,U>::Begin(s,rq,cx,field,trap);
  if (!finished && s.done==0) finished=Tracer<R,U>::Advance(s,rq,cx,field,trap,K);
  (void)finished;
  s.lastNode=field.LastNode();
  Tracer<R,U>::MergeStats(s,field.Stats());
  if (field.UnavailableKind()!=0 && s.res.unavailableKind==0) {
    s.res.unavailableKind=field.UnavailableKind();
    for (int d=0;d<3;++d) s.res.unavailableX_m[d]=field.UnavailableX(d);
  }
  if (field.Error()!=kErrNone) {
    // The legacy tracer terminates the process for these conditions (exit()); the host aborts
    // the run when it receives the result.
    if (s.res.error==kErrNone) s.res.error=field.Error();
    if (s.done!=1) Tracer<R,U>::Finish(s,trap,kTermNumericalFailure);
  }
}

// Per-run trap arena requirement: maximum over all registered contexts.
inline void TrapArenaStride(const EgContext* ctx,int nCtx,int turnCap,
                            long long& strideD,long long& strideI,int& strideT) {
  strideD=0; strideI=0; strideT=0;
  for (int i=0;i<nCtx;++i) {
    const EgTrapLayout l=MakeTrapLayout(ctx[i].trap,turnCap);
    if (l.doublesPerSlot>strideD) strideD=l.doublesPerSlot;
    if (l.intsPerSlot>strideI) strideI=l.intsPerSlot;
    if (l.turnCapacity>strideT) strideT=l.turnCapacity;
  }
}

} // namespace GPU
} // namespace Geospace

#endif
