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
  ParkerBackground(double advection, double divergence,double logAreaGradient=0)
      : plasmaAdvectionMPerS(advection),
        velocityDivergencePerS(divergence),logAreaGradientPerM(logAreaGradient) {}

  // The adapter supplies convection relative to its declared geometry. divU
  // remains the physical divergence from the same immutable background epoch.
  // Compatibility field name: advection relative to the declared geometry
  // [m/s], supplied by the adapter. Zero on plasma-comoving material lines;
  // U dot tangent on fixed geometry. Physical grad(U) is supplied separately.
  double plasmaAdvectionMPerS = 0.0;
  double velocityDivergencePerS = 0.0;
  // d ln(A)/ds [1/m] of the same tube whose physical number is sampled. For
  // A|B|=Phi this equals -d ln|B|/ds. It is a spatial diffusion metric, not
  // solar-wind convection and not part of the physical kappa derivative.
  double logAreaGradientPerM = 0.0;
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
//   ds=(a_coordinate+d(kappa)/ds+kappa*dln(A)/ds)dt+sqrt(2 kappa dt)dW,
//   dp=-(p/3) div(U) dt.
// These are forward physical-number particles: N=4*pi*p^2*A*f per ds dp.
// The A term follows by rewriting (1/A) d_s(A*kappa*d_s f) as a forward
// conservative number equation. Conversion from arc length to segment S is a
// different coordinate operation and cannot supply this missing drift.
ParkerIncrement AdvanceParker(const ParkerState& initial,
                              const ParkerBackground& background,
                              double speedMPerS,
                              double dtS,
                              const SpatialDiffusionProvider& provider,
                              RandomStream& random);

}  // namespace Transport
}  // namespace SEP

#endif  // SEP_UTIL_SEP_PARKER_CORE_H
