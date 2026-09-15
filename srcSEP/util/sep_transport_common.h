#ifndef SEP_UTIL_SEP_TRANSPORT_COMMON_H
#define SEP_UTIL_SEP_TRANSPORT_COMMON_H

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace SEP {
namespace Transport {

// Status is returned across every numerical boundary.  Production adapters may
// translate a non-OK status into the host application's fatal-error mechanism,
// while component tests can inspect the exact cause without terminating the
// process.  This prevents invalid coefficients and hidden numerical clamps from
// being mistaken for a successful particle step.
enum class StatusCode {
  Ok,
  InvalidArgument,
  InvalidParticleState,
  InvalidCoefficient,
  OutOfDomain,
  StepUnderflow,
  UnsupportedConfiguration,
  UnresolvedCoefficient,
  NondifferentiableCoefficient
};

struct Status {
  StatusCode code = StatusCode::Ok;
  std::string message;

  bool ok() const { return code == StatusCode::Ok; }
  static Status Ok();
  static Status Error(StatusCode code, const std::string& message);
};

struct ScalarResult {
  Status status;
  double value = 0.0;
};

// Minimal state shared by all field-line movers.  coordinate is the host field
// line coordinate (not assumed to be arc length); the PIC-facing adapter is the
// only layer allowed to translate a physical displacement in metres into this
// coordinate.  Velocity components are measured in the local plasma frame.
struct ParticleState {
  int species = -1;
  int fieldLineId = -1;
  double coordinate = 0.0;
  double vParallelMPerS = 0.0;
  double vNormalMPerS = 0.0;
  double massKg = 0.0;
};

Status ValidateParticleState(const ParticleState& state,
                             double speedOfLightMPerS);

// Relativistic conversions use p=gamma*m*v and never clip a luminal input.
// A caller must subcycle or reject an invalid state rather than silently
// replacing it with 0.99c, which changes the represented particle energy.
ScalarResult MomentumFromSpeed(double speedMPerS, double massKg,
                               double speedOfLightMPerS);
ScalarResult SpeedFromMomentum(double momentumKgMPerS, double massKg,
                               double speedOfLightMPerS);

// The focused transport equation uses L=-1/(d ln|B|/ds).  A uniform magnetic
// field has an infinite focusing length and is a valid, zero-focusing case.
ScalarResult FocusingLengthFromDlnBds(double dLnAbsBdsPerM);

// Exact update for dp/dt=-(p/3) div(U), in the declared local plasma frame.
// The exponential form remains positive and avoids a timestep-dependent sign
// reversal that can occur with a forward-Euler momentum update.
ScalarResult ApplyAdiabaticMomentum(double momentumKgMPerS,
                                    double velocityDivergencePerS,
                                    double dtS);

enum class BoundaryPolicy { Absorb, Reflect };

struct CoordinateAdvance {
  Status status;
  double positionM = 0.0;
  double directedSpeedMPerS = 0.0;
  bool crossedBoundary = false;
  unsigned reflections = 0;
};

// Advances a one-dimensional physical arc-length coordinate.  Reflect handles
// arbitrarily long finite steps by repeated mirror images; Absorb reports the
// attempted boundary crossing without inventing an in-domain coordinate.
CoordinateAdvance AdvanceCoordinate(double positionM, double displacementM,
                                    double directedSpeedMPerS,
                                    double minimumM, double maximumM,
                                    BoundaryPolicy policy);

struct StepLimit {
  StepLimit() {}
  StepLimit(const std::string& limitName, double limitDtS)
      : name(limitName), dtS(limitDtS) {}

  std::string name;
  double dtS = 0.0;
};

struct StepDiagnostics {
  std::size_t selections = 0;
  std::size_t acceptedSteps = 0;
  std::size_t rejectedSteps = 0;
  std::size_t underflowRejects = 0;
  std::size_t errorControlRejects = 0;
  std::vector<std::string> limitingNames;

  // Parallel vectors keep the dependency-light C++11 ABI simple while still
  // publishing a histogram of the physical/numerical mechanism that limited
  // each selected step.  A production diagnostic may serialize these pairs
  // without depending on map iteration order.
  std::vector<std::string> limiterHistogramNames;
  std::vector<std::size_t> limiterHistogramCounts;
};

// WP12 replaces mover-local magic numbers with one validated error-control
// contract.  Fractions and tolerances are dimensionless; minimumStepS is in
// seconds.  Defaults reproduce the former conservative limits, but they are
// now named, configurable, fingerprintable, and visible in startup output.
struct NumericalTolerances {
  double geometryFraction = 0.25;
  double deterministicRelativeTolerance = 1.0e-3;
  double stochasticPitchRms = 0.05;
  double coolingLogChange = 0.05;
  double focusingPitchChange = 0.05;
  double shockFraction = 1.0;
  double minimumStepS = 1.0e-12;
};

Status ValidateNumericalTolerances(const NumericalTolerances& tolerances);
const NumericalTolerances& ActiveNumericalTolerances();
Status SetActiveNumericalTolerances(const NumericalTolerances& tolerances);

// The deterministic local-error estimator compares one full step with two
// half steps.  The normalized error is |y_half-y_full| divided by an absolute
// plus relative scale.  It is intentionally independent of a particular mover
// so Parker drift, focusing, and cooling can share exactly one acceptance rule.
struct StepDoublingEstimate {
  Status status;
  double absoluteError = 0.0;
  double normalizedError = 0.0;
  bool accepted = false;
};

StepDoublingEstimate EstimateStepDoublingError(
    double fullStepValue, double twoHalfStepValue,
    double absoluteTolerance, double relativeTolerance);

void RecordAcceptedStep(StepDiagnostics* diagnostics);
void RecordRejectedStep(StepDiagnostics* diagnostics, bool errorControlled);

// ComposeSnapshotValidityLimit converts immutable snapshot metadata into the
// same StepLimit representation used for segment, shock, focusing, and cooling
// limits.  elapsedStepS is measured from the beginning of the current particle
// step; both snapshot times are physical seconds on the authoritative
// simulation clock.  Returning an explicit error for stale generations avoids
// silently sampling a newly published background halfway through a step.
ScalarResult ComposeSnapshotValidityLimit(double particleStepEpochS,
                                           double elapsedStepS,
                                           double snapshotValidUntilS,
                                           std::uint64_t expectedGeneration,
                                           std::uint64_t sampledGeneration);

// Selects the smallest finite positive composable limit.  minimumDtS is an
// error threshold, not a floor: raising a requested timestep would violate the
// very stability condition the selector is intended to enforce.
ScalarResult SelectSubstep(double remainingDtS,
                           const std::vector<StepLimit>& limits,
                           double minimumDtS,
                           StepDiagnostics* diagnostics);

class RandomStream {
 public:
  virtual ~RandomStream() {}
  virtual double UniformOpen01() = 0;
  virtual double Normal01() = 0;
};

// SplitMix64 gives a small, deterministic keyed stream suitable for numerical
// reproducibility tests.  A production campaign should construct one stream
// from its campaign seed, particle identity, operator ID, and event index; MPI
// rank and thread scheduling must not enter those keys.
class KeyedRandomStream : public RandomStream {
 public:
  KeyedRandomStream(std::uint64_t campaignSeed,
                    std::uint64_t particleId,
                    std::uint64_t operatorId,
                    std::uint64_t eventIndex);

  double UniformOpen01() override;
  double Normal01() override;

 private:
  std::uint64_t state_;
  bool hasSpareNormal_ = false;
  double spareNormal_ = 0.0;

  std::uint64_t NextU64();
};

}  // namespace Transport
}  // namespace SEP

#endif  // SEP_UTIL_SEP_TRANSPORT_COMMON_H
