#include "time_step.h"

#include <algorithm>
#include <cmath>
#include <limits>

namespace SEP3D {
namespace Transport {
namespace {

Core::Status Invalid(const char* message) {
  return Core::Status(Core::StatusCode::InvalidInput, message);
}

double RateLimit(double fraction, double rate) {
  return rate > 0.0 ? fraction / rate
                    : std::numeric_limits<double>::infinity();
}

}  // namespace

const char* Name(StepLimiter limiter) {
  switch (limiter) {
    case StepLimiter::Requested: return "requested";
    case StepLimiter::CellCrossing: return "cell-crossing";
    case StepLimiter::Diffusion: return "diffusion";
    case StepLimiter::Focusing: return "focusing";
    case StepLimiter::Cooling: return "cooling";
    case StepLimiter::FieldVariation: return "field-variation";
    case StepLimiter::ShockCrossing: return "shock-crossing";
    case StepLimiter::SnapshotBoundary: return "snapshot-boundary";
  }
  return "unknown";
}

TimeStepSelection SelectTimeStep(const TimeStepControls& controls,
                                 const TimeStepPhysics& physics) {
  TimeStepSelection result;
  const double controlValues[] = {
      controls.cellCrossingFraction, controls.diffusionFraction,
      controls.focusingFraction, controls.coolingFraction,
      controls.fieldVariationFraction, controls.shockCrossingFraction,
      controls.minimumSubstepS};
  for (double value : controlValues) {
    if (!std::isfinite(value) || value <= 0.0) {
      result.status = Invalid("time-step controls must be finite and positive");
      return result;
    }
  }
  const double physicalValues[] = {
      physics.requestedS, physics.cellSizeM,
      physics.characteristicSpeedMPerS, physics.kappaParallelM2PerS,
      physics.focusingRatePerS, physics.coolingRatePerS,
      physics.fractionalFieldVariationPerS,
      physics.timeToShockCrossingS, physics.timeToSnapshotBoundaryS};
  for (double value : physicalValues) {
    if (!std::isfinite(value) || value < 0.0) {
      result.status = Invalid("time-step physics input is negative or non-finite");
      return result;
    }
  }
  if (physics.requestedS == 0.0 || physics.cellSizeM == 0.0) {
    result.status = Invalid("requested step and cell size must be positive");
    return result;
  }

  const double infinity = std::numeric_limits<double>::infinity();
  result.candidatesS = {{
      physics.requestedS,
      physics.characteristicSpeedMPerS > 0.0
          ? controls.cellCrossingFraction * physics.cellSizeM /
                physics.characteristicSpeedMPerS : infinity,
      physics.kappaParallelM2PerS > 0.0
          ? controls.diffusionFraction * physics.cellSizeM * physics.cellSizeM /
                (2.0 * physics.kappaParallelM2PerS) : infinity,
      RateLimit(controls.focusingFraction, physics.focusingRatePerS),
      RateLimit(controls.coolingFraction, physics.coolingRatePerS),
      RateLimit(controls.fieldVariationFraction,
                physics.fractionalFieldVariationPerS),
      physics.timeToShockCrossingS > 0.0
          ? controls.shockCrossingFraction * physics.timeToShockCrossingS
          : infinity,
      physics.timeToSnapshotBoundaryS > 0.0
          ? physics.timeToSnapshotBoundaryS : infinity}};

  std::size_t selected = 0;
  for (std::size_t i = 1; i < result.candidatesS.size(); ++i) {
    if (result.candidatesS[i] < result.candidatesS[selected]) selected = i;
  }
  result.valueS = result.candidatesS[selected];
  result.limiter = static_cast<StepLimiter>(selected);
  if (!std::isfinite(result.valueS) || result.valueS <= 0.0) {
    result.status = Invalid("no finite positive transport substep is available");
  } else if (result.valueS < controls.minimumSubstepS) {
    result.status = Core::Status(
        Core::StatusCode::StepUnderflow,
        std::string("transport substep underflow at ") + Name(result.limiter));
  } else {
    result.status = Core::Status::OK();
  }
  return result;
}

}  // namespace Transport
}  // namespace SEP3D
