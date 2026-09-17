#include "sep_transport_common.h"

#include <algorithm>
#include <cmath>
#include <limits>

namespace SEP {
namespace Transport {
namespace {

bool Finite(double value) { return std::isfinite(value); }

std::uint64_t Mix(std::uint64_t value) {
  value += UINT64_C(0x9e3779b97f4a7c15);
  value = (value ^ (value >> 30)) * UINT64_C(0xbf58476d1ce4e5b9);
  value = (value ^ (value >> 27)) * UINT64_C(0x94d049bb133111eb);
  return value ^ (value >> 31);
}

NumericalTolerances gNumericalTolerances;

void CountLimiter(const std::string& name, StepDiagnostics* diagnostics) {
  if (!diagnostics) return;
  for (std::size_t i = 0; i < diagnostics->limiterHistogramNames.size(); ++i) {
    if (diagnostics->limiterHistogramNames[i] == name) {
      ++diagnostics->limiterHistogramCounts[i];
      return;
    }
  }
  diagnostics->limiterHistogramNames.push_back(name);
  diagnostics->limiterHistogramCounts.push_back(1);
}

}  // namespace

Status Status::Ok() { return Status(); }

Status Status::Error(StatusCode code, const std::string& message) {
  Status status;
  status.code = code;
  status.message = message;
  return status;
}

Status ValidateNumericalTolerances(const NumericalTolerances& t) {
  if (!Finite(t.geometryFraction) || t.geometryFraction <= 0.0 ||
      t.geometryFraction > 1.0 ||
      !Finite(t.deterministicRelativeTolerance) ||
      t.deterministicRelativeTolerance <= 0.0 ||
      !Finite(t.stochasticPitchRms) || t.stochasticPitchRms <= 0.0 ||
      t.stochasticPitchRms > 1.0 ||
      !Finite(t.coolingLogChange) || t.coolingLogChange <= 0.0 ||
      !Finite(t.focusingPitchChange) || t.focusingPitchChange <= 0.0 ||
      t.focusingPitchChange > 1.0 ||
      !Finite(t.shockFraction) || t.shockFraction <= 0.0 ||
      t.shockFraction > 1.0 ||
      !Finite(t.minimumStepS) || t.minimumStepS < 0.0) {
    return Status::Error(StatusCode::InvalidArgument,
        "transport tolerances require finite positive fractions/tolerances, "
        "geometry, stochastic, focusing and shock fractions no greater than "
        "one, and a non-negative minimum step");
  }
  return Status::Ok();
}

const NumericalTolerances& ActiveNumericalTolerances() {
  return gNumericalTolerances;
}

Status SetActiveNumericalTolerances(const NumericalTolerances& tolerances) {
  const Status status = ValidateNumericalTolerances(tolerances);
  if (status.ok()) gNumericalTolerances = tolerances;
  return status;
}

StepDoublingEstimate EstimateStepDoublingError(
    double fullStepValue, double twoHalfStepValue,
    double absoluteTolerance, double relativeTolerance) {
  StepDoublingEstimate estimate;
  if (!Finite(fullStepValue) || !Finite(twoHalfStepValue) ||
      !Finite(absoluteTolerance) || absoluteTolerance < 0.0 ||
      !Finite(relativeTolerance) || relativeTolerance <= 0.0) {
    estimate.status = Status::Error(StatusCode::InvalidArgument,
        "step-doubling estimator requires finite states and valid tolerances");
    return estimate;
  }
  estimate.absoluteError = std::fabs(twoHalfStepValue - fullStepValue);
  const double scale = absoluteTolerance + relativeTolerance *
      std::max(std::fabs(fullStepValue), std::fabs(twoHalfStepValue));
  if (!(scale > 0.0) || !Finite(scale)) {
    estimate.status = Status::Error(StatusCode::InvalidArgument,
        "step-doubling error scale is zero or non-finite");
    return estimate;
  }
  estimate.normalizedError = estimate.absoluteError / scale;
  estimate.accepted = estimate.normalizedError <= 1.0;
  estimate.status = Status::Ok();
  return estimate;
}

void RecordAcceptedStep(StepDiagnostics* diagnostics) {
  if (diagnostics) ++diagnostics->acceptedSteps;
}

void RecordRejectedStep(StepDiagnostics* diagnostics, bool errorControlled) {
  if (!diagnostics) return;
  ++diagnostics->rejectedSteps;
  if (errorControlled) ++diagnostics->errorControlRejects;
}

Status ValidateParticleState(const ParticleState& state,
                             double speedOfLightMPerS) {
  if (state.species < 0 || state.fieldLineId < 0 ||
      !Finite(state.coordinate) || !Finite(state.vParallelMPerS) ||
      !Finite(state.vNormalMPerS) || !Finite(state.massKg) ||
      state.massKg <= 0.0 || !Finite(speedOfLightMPerS) ||
      speedOfLightMPerS <= 0.0) {
    return Status::Error(StatusCode::InvalidParticleState,
                         "non-finite or non-physical field-line particle state");
  }

  const double speed = std::hypot(state.vParallelMPerS, state.vNormalMPerS);
  if (!Finite(speed) || speed < 0.0 || speed >= speedOfLightMPerS) {
    return Status::Error(StatusCode::InvalidParticleState,
                         "particle speed must be finite and strictly subluminal");
  }
  return Status::Ok();
}

ScalarResult MomentumFromSpeed(double speedMPerS, double massKg,
                               double speedOfLightMPerS) {
  ScalarResult result;
  if (!Finite(speedMPerS) || !Finite(massKg) ||
      !Finite(speedOfLightMPerS) || speedMPerS < 0.0 || massKg <= 0.0 ||
      speedOfLightMPerS <= 0.0 || speedMPerS >= speedOfLightMPerS) {
    result.status = Status::Error(StatusCode::InvalidArgument,
        "speed-to-momentum conversion requires 0 <= v < c and positive mass");
    return result;
  }

  const double beta = speedMPerS / speedOfLightMPerS;
  const double gamma = 1.0 / std::sqrt(1.0 - beta * beta);
  result.value = gamma * massKg * speedMPerS;
  if (!Finite(result.value)) {
    result.status = Status::Error(StatusCode::InvalidParticleState,
                                  "relativistic momentum overflowed");
    return result;
  }
  result.status = Status::Ok();
  return result;
}

ScalarResult SpeedFromMomentum(double momentumKgMPerS, double massKg,
                               double speedOfLightMPerS) {
  ScalarResult result;
  if (!Finite(momentumKgMPerS) || !Finite(massKg) ||
      !Finite(speedOfLightMPerS) || momentumKgMPerS < 0.0 || massKg <= 0.0 ||
      speedOfLightMPerS <= 0.0) {
    result.status = Status::Error(StatusCode::InvalidArgument,
        "momentum-to-speed conversion requires nonnegative momentum and positive mass");
    return result;
  }

  // v = c p/sqrt((mc)^2+p^2) is stable for both non-relativistic and
  // ultra-relativistic finite momentum and is strictly less than c.
  const double mc = massKg * speedOfLightMPerS;
  result.value = speedOfLightMPerS * momentumKgMPerS /
                 std::hypot(mc, momentumKgMPerS);
  if (!Finite(result.value) || result.value >= speedOfLightMPerS) {
    result.status = Status::Error(StatusCode::InvalidParticleState,
                                  "relativistic speed conversion failed");
    return result;
  }
  result.status = Status::Ok();
  return result;
}

ScalarResult FocusingLengthFromDlnBds(double dLnAbsBdsPerM) {
  ScalarResult result;
  if (!Finite(dLnAbsBdsPerM)) {
    result.status = Status::Error(StatusCode::InvalidArgument,
                                  "d ln|B|/ds must be finite");
    return result;
  }
  result.value = dLnAbsBdsPerM == 0.0
      ? std::numeric_limits<double>::infinity()
      : -1.0 / dLnAbsBdsPerM;
  result.status = Status::Ok();
  return result;
}

ScalarResult ApplyAdiabaticMomentum(double momentumKgMPerS,
                                    double velocityDivergencePerS,
                                    double dtS) {
  ScalarResult result;
  if (!Finite(momentumKgMPerS) || momentumKgMPerS < 0.0 ||
      !Finite(velocityDivergencePerS) || !Finite(dtS) || dtS < 0.0) {
    result.status = Status::Error(StatusCode::InvalidArgument,
        "adiabatic update requires p>=0, finite divergence, and dt>=0");
    return result;
  }
  const double exponent = -velocityDivergencePerS * dtS / 3.0;
  if (!Finite(exponent) || exponent > std::log(std::numeric_limits<double>::max())) {
    result.status = Status::Error(StatusCode::InvalidParticleState,
                                  "adiabatic momentum update overflowed");
    return result;
  }
  result.value = momentumKgMPerS * std::exp(exponent);
  if (!Finite(result.value)) {
    result.status = Status::Error(StatusCode::InvalidParticleState,
                                  "adiabatic momentum result is not finite");
    return result;
  }
  result.status = Status::Ok();
  return result;
}

CoordinateAdvance AdvanceCoordinate(double positionM, double displacementM,
                                    double directedSpeedMPerS,
                                    double minimumM, double maximumM,
                                    BoundaryPolicy policy) {
  CoordinateAdvance result;
  result.positionM = positionM;
  result.directedSpeedMPerS = directedSpeedMPerS;
  if (!Finite(positionM) || !Finite(displacementM) ||
      !Finite(directedSpeedMPerS) || !Finite(minimumM) ||
      !Finite(maximumM) || !(minimumM < maximumM) ||
      positionM < minimumM || positionM > maximumM) {
    result.status = Status::Error(StatusCode::InvalidArgument,
                                  "invalid coordinate advance domain or state");
    return result;
  }

  const double attempted = positionM + displacementM;
  if (attempted >= minimumM && attempted <= maximumM) {
    result.positionM = attempted;
    result.status = Status::Ok();
    return result;
  }

  result.crossedBoundary = true;
  if (policy == BoundaryPolicy::Absorb) {
    result.status = Status::Error(StatusCode::OutOfDomain,
                                  "absorbing field-line boundary crossed");
    return result;
  }

  const double width = maximumM - minimumM;
  const double period = 2.0 * width;
  double folded = std::fmod(attempted - minimumM, period);
  if (folded < 0.0) folded += period;
  result.positionM = folded <= width
      ? minimumM + folded
      : maximumM - (folded - width);

  // Count crossings for diagnostics without an O(number-of-reflections) loop.
  // A reflected velocity reverses only after an odd number of boundary hits.
  const double startCell = std::floor((positionM - minimumM) / width);
  const double endCell = std::floor((attempted - minimumM) / width);
  result.reflections = static_cast<unsigned>(std::fabs(endCell - startCell));
  if (result.reflections == 0) result.reflections = 1;
  if ((result.reflections % 2U) != 0U)
    result.directedSpeedMPerS = -result.directedSpeedMPerS;
  result.status = Status::Ok();
  return result;
}

ScalarResult SelectSubstep(double remainingDtS,
                           const std::vector<StepLimit>& limits,
                           double minimumDtS,
                           StepDiagnostics* diagnostics) {
  ScalarResult result;
  if (!Finite(remainingDtS) || remainingDtS <= 0.0 ||
      !Finite(minimumDtS) || minimumDtS < 0.0) {
    result.status = Status::Error(StatusCode::InvalidArgument,
                                  "substep interval must be finite and positive");
    return result;
  }

  result.value = remainingDtS;
  std::string limitingName = "remaining";
  for (const StepLimit& limit : limits) {
    if (limit.name.empty() || !Finite(limit.dtS) || limit.dtS <= 0.0) {
      result.status = Status::Error(StatusCode::InvalidArgument,
                                    "substep limits require a name and dt>0");
      return result;
    }
    if (limit.dtS < result.value) {
      result.value = limit.dtS;
      limitingName = limit.name;
    }
  }

  if (diagnostics) {
    diagnostics->selections++;
    diagnostics->limitingNames.push_back(limitingName);
    CountLimiter(limitingName, diagnostics);
  }
  if (result.value < minimumDtS) {
    if (diagnostics) diagnostics->underflowRejects++;
    result.status = Status::Error(StatusCode::StepUnderflow,
        "stability limit is below the declared minimum substep");
    return result;
  }
  result.status = Status::Ok();
  return result;
}

ScalarResult ComposeSnapshotValidityLimit(double particleStepEpochS,
                                           double elapsedStepS,
                                           double snapshotValidUntilS,
                                           std::uint64_t expectedGeneration,
                                           std::uint64_t sampledGeneration) {
  ScalarResult result;
  if (!std::isfinite(particleStepEpochS) ||
      !std::isfinite(elapsedStepS) || elapsedStepS < 0.0 ||
      std::isnan(snapshotValidUntilS)) {
    result.status = Status::Error(StatusCode::InvalidArgument,
        "snapshot validity limit contains a non-finite time");
    return result;
  }
  if (expectedGeneration != sampledGeneration) {
    result.status = Status::Error(StatusCode::OutOfDomain,
        "background field-line generation changed during the particle step");
    return result;
  }

  result.value = snapshotValidUntilS -
      (particleStepEpochS + elapsedStepS);
  if (result.value < 0.0) {
    result.status = Status::Error(StatusCode::OutOfDomain,
        "particle step requested data beyond the snapshot validity interval");
    return result;
  }
  result.status = Status::Ok();
  return result;
}

KeyedRandomStream::KeyedRandomStream(std::uint64_t campaignSeed,
                                     std::uint64_t particleId,
                                     std::uint64_t operatorId,
                                     std::uint64_t eventIndex)
    : state_(Mix(UINT64_C(0x534550524e474b59))) {
  // Combine fields sequentially with distinct semantic tags.  The former XOR
  // construction was commutative, so swapping (particleId,operatorId) yielded
  // the same stream.  Ordered combination preserves deterministic replay while
  // making every key field independently reachable and collision resistant.
  state_ = Mix(state_ ^ Mix(UINT64_C(1)) ^ Mix(campaignSeed));
  state_ = Mix(state_ ^ Mix(UINT64_C(2)) ^ Mix(particleId));
  state_ = Mix(state_ ^ Mix(UINT64_C(3)) ^ Mix(operatorId));
  state_ = Mix(state_ ^ Mix(UINT64_C(4)) ^ Mix(eventIndex));
}

std::uint64_t KeyedRandomStream::NextU64() {
  state_ += UINT64_C(0x9e3779b97f4a7c15);
  return Mix(state_);
}

double KeyedRandomStream::UniformOpen01() {
  // Retain 53 random bits and offset by half a bin so neither zero nor one is
  // returned.  Logarithms in Box-Muller and event-time kernels are then finite.
  const std::uint64_t mantissa = NextU64() >> 11;
  return (static_cast<double>(mantissa) + 0.5) /
         9007199254740992.0;
}

double KeyedRandomStream::Normal01() {
  if (hasSpareNormal_) {
    hasSpareNormal_ = false;
    return spareNormal_;
  }
  const double u1 = UniformOpen01();
  const double u2 = UniformOpen01();
  const double radius = std::sqrt(-2.0 * std::log(u1));
  const double angle = 6.283185307179586476925286766559 * u2;
  spareNormal_ = radius * std::sin(angle);
  hasSpareNormal_ = true;
  return radius * std::cos(angle);
}

}  // namespace Transport
}  // namespace SEP
