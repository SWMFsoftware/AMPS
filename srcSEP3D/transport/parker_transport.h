// ============================================================================
// AMPS-independent three-dimensional Parker transport core.
//
// For a gyrotropic, pitch-angle-averaged distribution with parallel-only
// diffusion, the spatial tensor is K = kappa_parallel b b.  The equivalent
// Ito pseudo-particle equation is
//
//   dX = [U + div(K)] dt + sqrt(2 kappa_parallel) b dW,
//   dP = -(P/3) div(U) dt,
//
// where
//
//   div(K) = b d(kappa_parallel)/ds
//            + kappa_parallel[(b dot grad)b + b div(b)].
//
// Retaining all three terms is essential: omitting the tensor-geometry terms
// breaks uniform-density equilibrium in a curved/nonuniform field.  This
// release intentionally sets perpendicular diffusion and gradient/curvature
// drifts to exactly zero; nonzero values are rejected rather than ignored.
// ============================================================================

#ifndef SEP3D_TRANSPORT_PARKER_TRANSPORT_H
#define SEP3D_TRANSPORT_PARKER_TRANSPORT_H

#include "../core/sep3d_types.h"
#include "keyed_random.h"

namespace SEP3D {
namespace Transport {

struct ParkerParticleState {
  Core::Vec3 positionM;
  double momentumKgMPerS = 0.0;
};

struct ParkerLocalState {
  Core::Vec3 bulkVelocityMPerS;
  Core::Vec3 bHat;
  Core::Vec3 curvaturePerM;
  double divBhatPerM = 0.0;
  double divUPerS = 0.0;
  double kappaParallelM2PerS = 0.0;
  double dKappaParallelDsMPerS = 0.0;

  // Reserved hooks are explicit so a caller cannot accidentally pass a
  // future coefficient to a parallel-only kernel and have it discarded.
  double kappaPerpendicularM2PerS = 0.0;
  Core::Vec3 driftVelocityMPerS;
};

struct ParkerStepResult {
  Core::Status status;
  ParkerParticleState state;
  Core::Tensor3 diffusionTensorM2PerS;
  Core::Vec3 itoDriftMPerS;
  Core::Vec3 deterministicDisplacementM;
  Core::Vec3 stochasticDisplacementM;
};

Core::Tensor3 AssembleParallelDiffusionTensor(double kappaParallelM2PerS,
                                               const Core::Vec3& bHat);
Core::Vec3 ParallelTensorItoDrift(const ParkerLocalState& local);
ParkerStepResult AdvanceParker(const ParkerParticleState& initial,
                               const ParkerLocalState& local,
                               double dtS,
                               KeyedRandomStream* random);

}  // namespace Transport
}  // namespace SEP3D

#endif  // SEP3D_TRANSPORT_PARKER_TRANSPORT_H
