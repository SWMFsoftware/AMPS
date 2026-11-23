#ifndef SRCSEP_ABSORBING_BRIDGE_H
#define SRCSEP_ABSORBING_BRIDGE_H
#include "sep_transport_coefficients.h"
#include <functional>
namespace SEP { namespace Transport {
// Per-particle, per-physical-step numerical history only. It owns no physical
// state, RNG, field sample or persistent source identity. Reinitializing it at
// the next complete PIC step also makes checkpoint continuation unambiguous.
struct BoundaryIntervalControl {
  double previousAcceptedS=0;
  unsigned successfulSinceGrowth=0;
};
using BoundaryIntervalCheck=std::function<Status(double,bool*)>;
// Refine BEFORE stochastic draws. The callback must query the actual immutable
// geometry and may accept only its unchanged spatial-error criterion. The
// original 64-attempt/minimum-step/progress guards remain authoritative; neither
// this helper nor its optional history advances the physical clock.
Status RefineBoundaryInterval(double beginEpochS,double proposedIntervalS,
    double minimumIntervalS,const BoundaryIntervalCheck&,double* acceptedIntervalS,
    BoundaryIntervalControl* history=nullptr);
struct AbsorbingBridgeResult {
  Status status;
  bool hit=false;
  double fraction=1; // first-hit time / supplied numerical interval, [0,1]
};
// Exact first passage for one Euler Brownian interval against a linear moving
// absorbing boundary, CONDITIONAL on the already sampled endpoint. Distances
// are signed upstream arc distances [m]; diffusionTime=kappa*dt [m^2]. Drift
// and linear boundary speed cancel from the conditional bridge law. This
// consumes the SAME sampled endpoint, never redraws a rejected particle path.
AbsorbingBridgeResult AbsorbingLinearBridge(double distanceBeginM,
    double distanceEndM,double diffusionTimeM2,RandomStream&);
} }
#endif
