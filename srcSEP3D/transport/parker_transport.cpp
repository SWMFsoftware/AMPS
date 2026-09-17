#include "parker_transport.h"

#include <cmath>
#include <limits>

namespace SEP3D {
namespace Transport {
namespace {

bool Finite(const Core::Vec3& value) {
  return std::isfinite(value.x) && std::isfinite(value.y) &&
         std::isfinite(value.z);
}

Core::Status Invalid(const char* message) {
  return Core::Status(Core::StatusCode::InvalidInput, message);
}

}  // namespace

Core::Tensor3 AssembleParallelDiffusionTensor(
    double kappaParallelM2PerS, const Core::Vec3& bHat) {
  Core::Tensor3 tensor;
  const double b[3] = {bHat.x, bHat.y, bHat.z};
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j)
      tensor(i, j) = kappaParallelM2PerS * b[i] * b[j];
  return tensor;
}

Core::Vec3 ParallelTensorItoDrift(const ParkerLocalState& local) {
  // Product rule for div(kappa b b).  dKappa/ds is b dot grad(kappa),
  // curvature is (b dot grad)b, and divBhat is div(b).
  return local.bHat * local.dKappaParallelDsMPerS +
      (local.curvaturePerM + local.bHat * local.divBhatPerM) *
          local.kappaParallelM2PerS;
}

ParkerStepResult AdvanceParker(const ParkerParticleState& initial,
                               const ParkerLocalState& local,
                               double dtS,
                               KeyedRandomStream* random) {
  ParkerStepResult result;
  result.state = initial;
  const double bNorm = local.bHat.Norm();
  if (!Finite(initial.positionM) ||
      !std::isfinite(initial.momentumKgMPerS) ||
      initial.momentumKgMPerS < 0.0 || !Finite(local.bulkVelocityMPerS) ||
      !Finite(local.bHat) || !Finite(local.curvaturePerM) ||
      !std::isfinite(local.divBhatPerM) || !std::isfinite(local.divUPerS) ||
      !std::isfinite(local.kappaParallelM2PerS) ||
      local.kappaParallelM2PerS < 0.0 ||
      !std::isfinite(local.dKappaParallelDsMPerS) ||
      !std::isfinite(local.kappaPerpendicularM2PerS) ||
      local.kappaPerpendicularM2PerS != 0.0 ||
      !Finite(local.driftVelocityMPerS) ||
      local.driftVelocityMPerS.NormSq() != 0.0 ||
      !std::isfinite(dtS) || dtS < 0.0 ||
      std::fabs(bNorm - 1.0) > 1.0e-12) {
    result.status = Invalid(
        "invalid Parker state or nonzero reserved perpendicular/drift input");
    return result;
  }
  if (local.kappaParallelM2PerS > 0.0 && random == nullptr) {
    result.status = Invalid("diffusive Parker step requires a keyed random stream");
    return result;
  }

  result.diffusionTensorM2PerS = AssembleParallelDiffusionTensor(
      local.kappaParallelM2PerS, local.bHat);
  result.itoDriftMPerS = ParallelTensorItoDrift(local);
  result.deterministicDisplacementM =
      (local.bulkVelocityMPerS + result.itoDriftMPerS) * dtS;
  result.stochasticDisplacementM = local.kappaParallelM2PerS == 0.0
      ? Core::Vec3{}
      : local.bHat *
          (std::sqrt(2.0 * local.kappaParallelM2PerS * dtS) *
           random->Normal01());
  result.state.positionM += result.deterministicDisplacementM +
                            result.stochasticDisplacementM;

  // The exponential is the exact characteristic for frozen div(U).  It is
  // positive by construction and avoids a negative momentum under a large
  // expansion substep; the named cooling limiter still controls accuracy when
  // div(U) varies along the trajectory.
  const double exponent = -local.divUPerS * dtS / 3.0;
  if (!std::isfinite(exponent) ||
      exponent > std::log(std::numeric_limits<double>::max())) {
    result.status = Invalid("Parker cooling exponent is not representable");
    return result;
  }
  result.state.momentumKgMPerS *= std::exp(exponent);
  if (!Finite(result.state.positionM) ||
      !std::isfinite(result.state.momentumKgMPerS)) {
    result.status = Invalid("Parker step produced a non-finite particle state");
    return result;
  }
  result.status = Core::Status::OK();
  return result;
}

}  // namespace Transport
}  // namespace SEP3D
