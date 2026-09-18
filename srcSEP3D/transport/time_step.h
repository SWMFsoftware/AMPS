// ============================================================================
// Named transport-substep limits.
//
// Every stability/accuracy restriction remains visible in diagnostics.  The
// selector never replaces an underflowing step with an undocumented floor:
// it returns StepUnderflow together with the exact limiting value and name.
// ============================================================================

#ifndef SEP3D_TRANSPORT_TIME_STEP_H
#define SEP3D_TRANSPORT_TIME_STEP_H

#include "../core/sep3d_types.h"

#include <array>
#include <string>

namespace SEP3D {
namespace Transport {

enum class StepLimiter {
  Requested,
  CellCrossing,
  Diffusion,
  Focusing,
  Cooling,
  FieldVariation,
  ShockCrossing,
  SnapshotBoundary
};

const char* Name(StepLimiter limiter);

struct TimeStepControls {
  double cellCrossingFraction = 0.4;
  double diffusionFraction = 0.2;
  double focusingFraction = 0.2;
  double coolingFraction = 0.2;
  double fieldVariationFraction = 0.2;
  double shockCrossingFraction = 0.5;
  double minimumSubstepS = 1.0e-12;
};

struct TimeStepPhysics {
  double requestedS = 0.0;
  double cellSizeM = 0.0;
  double characteristicSpeedMPerS = 0.0;
  double kappaParallelM2PerS = 0.0;
  // The tensor stability bound is controlled by its largest eigenvalue. For
  // gyrotropic diffusion this is max(kappa_parallel,kappa_perpendicular).
  double kappaPerpendicularM2PerS = 0.0;
  double focusingRatePerS = 0.0;
  double coolingRatePerS = 0.0;
  double fractionalFieldVariationPerS = 0.0;
  double timeToShockCrossingS = 0.0;
  double timeToSnapshotBoundaryS = 0.0;
};

struct TimeStepSelection {
  Core::Status status;
  double valueS = 0.0;
  StepLimiter limiter = StepLimiter::Requested;
  // Fixed order is part of the diagnostic schema and restart reproducibility.
  std::array<double, 8> candidatesS{};
};

TimeStepSelection SelectTimeStep(const TimeStepControls& controls,
                                 const TimeStepPhysics& physics);

}  // namespace Transport
}  // namespace SEP3D

#endif  // SEP3D_TRANSPORT_TIME_STEP_H
