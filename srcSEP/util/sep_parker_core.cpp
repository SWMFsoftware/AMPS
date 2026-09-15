#include "sep_parker_core.h"

#include <cmath>

namespace SEP {
namespace Transport {

ParkerIncrement AdvanceParker(const ParkerState& initial,
                              const ParkerBackground& background,
                              double speedMPerS,
                              double dtS,
                              const SpatialDiffusionProvider& provider,
                              RandomStream& random) {
  ParkerIncrement result;
  result.state = initial;
  if (!std::isfinite(initial.arcLengthM) ||
      !std::isfinite(initial.momentumKgMPerS) ||
      initial.momentumKgMPerS < 0.0 || !std::isfinite(speedMPerS) ||
      speedMPerS < 0.0 || !std::isfinite(dtS) || dtS < 0.0 ||
      !std::isfinite(background.plasmaAdvectionMPerS) ||
      !std::isfinite(background.velocityDivergencePerS)) {
    result.status = Status::Error(StatusCode::InvalidArgument,
                                  "invalid Parker state, background, or timestep");
    return result;
  }

  const SpatialDiffusionSample coefficient =
      provider.Evaluate(initial.arcLengthM, speedMPerS);
  if (!coefficient.status.ok() ||
      !std::isfinite(coefficient.kappaParallelM2PerS) ||
      !std::isfinite(coefficient.dKappaParallelDsMPerS) ||
      coefficient.kappaParallelM2PerS < 0.0 ||
      coefficient.provenance.empty()) {
    result.status = Status::Error(StatusCode::InvalidCoefficient,
        coefficient.status.message.empty()
            ? "invalid spatial-diffusion provider result"
            : coefficient.status.message);
    return result;
  }

  result.coefficientProvenance = coefficient.provenance;
  result.deterministicDisplacementM =
      (background.plasmaAdvectionMPerS +
       coefficient.dKappaParallelDsMPerS) * dtS;
  result.stochasticDisplacementM =
      coefficient.kappaParallelM2PerS == 0.0
          ? 0.0
          : std::sqrt(2.0 * coefficient.kappaParallelM2PerS * dtS) *
                random.Normal01();
  result.displacementM = result.deterministicDisplacementM +
                         result.stochasticDisplacementM;
  if (!std::isfinite(result.displacementM)) {
    result.status = Status::Error(StatusCode::InvalidParticleState,
                                  "Parker displacement is not finite");
    return result;
  }

  const ScalarResult momentum = ApplyAdiabaticMomentum(
      initial.momentumKgMPerS, background.velocityDivergencePerS, dtS);
  if (!momentum.status.ok()) {
    result.status = momentum.status;
    return result;
  }
  result.state.arcLengthM += result.displacementM;
  result.state.momentumKgMPerS = momentum.value;
  result.status = Status::Ok();
  return result;
}

}  // namespace Transport
}  // namespace SEP
