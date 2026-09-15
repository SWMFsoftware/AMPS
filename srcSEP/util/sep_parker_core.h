#ifndef SEP_UTIL_SEP_PARKER_CORE_H
#define SEP_UTIL_SEP_PARKER_CORE_H

#include "sep_transport_coefficients.h"

namespace SEP {
namespace Transport {

struct ParkerState {
  ParkerState() {}
  ParkerState(double arcLength, double momentum)
      : arcLengthM(arcLength), momentumKgMPerS(momentum) {}

  double arcLengthM = 0.0;
  double momentumKgMPerS = 0.0;
};

struct ParkerBackground {
  ParkerBackground() {}
  ParkerBackground(double advection, double divergence)
      : plasmaAdvectionMPerS(advection),
        velocityDivergencePerS(divergence) {}

  // U_parallel is convection of the field-line plasma in the chosen positive-s
  // direction.  divU is evaluated from the same immutable background epoch.
  double plasmaAdvectionMPerS = 0.0;
  double velocityDivergencePerS = 0.0;
};

struct ParkerIncrement {
  Status status;
  ParkerState state;
  double displacementM = 0.0;
  double deterministicDisplacementM = 0.0;
  double stochasticDisplacementM = 0.0;
  std::string coefficientProvenance;
};

// One Ito step of the field-aligned Parker pseudo-particle equation:
//   ds=(U_parallel+d(kappa)/ds)dt+sqrt(2 kappa dt)dW,
//   dp=-(p/3) div(U) dt.
// Physical arc length is used here; the PIC adapter applies the flux-tube metric
// exactly once when translating displacement to its non-metric line coordinate.
ParkerIncrement AdvanceParker(const ParkerState& initial,
                              const ParkerBackground& background,
                              double speedMPerS,
                              double dtS,
                              const SpatialDiffusionProvider& provider,
                              RandomStream& random);

}  // namespace Transport
}  // namespace SEP

#endif  // SEP_UTIL_SEP_PARKER_CORE_H
