#ifndef SEP_UTIL_SEP_NUMERICAL_EXTENSIONS_H
#define SEP_UTIL_SEP_NUMERICAL_EXTENSIONS_H

#include "sep_transport_common.h"

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace SEP {
namespace NumericalExtensions {

// WP53: a Brownian-tree address is numerical, not physical.  The interval ID
// remains unchanged when a step is subdivided, while depth/index select a node
// in the deterministic bridge.  A parent increment always equals the sum of
// its two children, so refinement compares approximations on the same path.
struct BrownianAddress {
  std::uint64_t campaign = 0;
  std::uint64_t particle = 0;
  std::uint64_t operatorId = 0;
  std::uint64_t physicalInterval = 0;
  unsigned depth = 0;
  std::uint64_t index = 0;
};

Transport::ScalarResult BrownianIncrement(
    const BrownianAddress& address, double rootIntervalS);

struct AdaptiveSdeConfiguration {
  double absoluteTolerance = 1.0e-8;
  double relativeTolerance = 1.0e-4;
  unsigned maximumDepth = 20;
};

struct AdaptiveSdeResult {
  Transport::Status status;
  double value = 0.0;
  unsigned acceptedLeaves = 0;
  unsigned rejectedTrials = 0;
  unsigned maximumDepthReached = 0;
  double maximumNormalizedError = 0.0;
};

// Integrates dX=a*X*dt+b*X*dW using Euler-Maruyama with a Brownian-bridge
// full-step/two-half-step estimator.  State is returned only after the full
// interval succeeds; failure leaves the caller's input untouched.
AdaptiveSdeResult AdvanceAdaptiveMultiplicativeSde(
    double initialValue, double driftPerS, double diffusionPerSqrtS,
    double dtS, const BrownianAddress& root,
    const AdaptiveSdeConfiguration& configuration);

// WP54: first-order upwind remains an explicit baseline.  MUSCL+SSPRK2 is the
// second-order option; limiter selection is named because it changes both
// accuracy and restart compatibility.
enum class AdvectionOrder { FirstOrderUpwind, SecondOrderMuscl };
enum class TvdLimiter { Minmod, MonotonizedCentral };

struct AdvectionResult {
  Transport::Status status;
  std::vector<double> cellAverage;
  std::uint64_t faceUpdates = 0;
  std::uint64_t limiterActivations = 0;
  double conservationResidual = 0.0;
};

AdvectionResult AdvectPeriodic(
    const std::vector<double>& cellAverage,
    double speedMPerS, double dtS, double cellLengthM,
    AdvectionOrder order, TvdLimiter limiter);

// WP55: remap uses a common monotonically increasing physical-volume
// coordinate.  Piecewise-linear reconstruction is limited and integrated over
// exact overlaps; uncovered old/new measure is reported as boundary loss or
// initialization rather than an unexplained numerical correction.
enum class RemapOrder { PiecewiseConstant, PiecewiseLinear };

struct RemapResult {
  Transport::Status status;
  std::vector<double> cellAverage;
  double oldIntegral = 0.0;
  double newIntegral = 0.0;
  double boundaryExchange = 0.0;
  double initializedIntegral = 0.0;
  double numericalResidual = 0.0;
  std::uint64_t limiterActivations = 0;
};

RemapResult RemapConservative(
    const std::vector<double>& oldEdges,
    const std::vector<double>& oldCellAverage,
    const std::vector<double>& newEdges,
    RemapOrder order, double newDomainInitializationValue = 0.0);

// WP56: production shells map typed numerical status to one explicit fate.
// The context fields are intentionally sufficient to replay a failure without
// relying on a transient particle pointer or worker number.
enum class MoverDisposition {
  Completed,
  BoundaryExit,
  RejectedStep,
  QuarantinedParticle,
  InvalidConfiguration,
  FatalInternalError
};
enum class ParticleFate { Retain, DeleteAtBoundary, Quarantine, Unchanged };
enum class FailurePolicy { FailFast, QuarantineParticleLocal };

struct FailureContext {
  std::uint64_t particleId = 0;
  int species = -1;
  int fieldLine = -1;
  int segment = -1;
  double coordinate = 0.0;
  double epochS = 0.0;
  std::string backgroundIdentity;
  std::string coefficientIdentity;
  std::string operatorName;
  std::uint64_t diagnosticCounter = 0;
};

struct MoverResult {
  Transport::Status status;
  MoverDisposition disposition = MoverDisposition::FatalInternalError;
  ParticleFate fate = ParticleFate::Unchanged;
  double committedIntervalS = 0.0;
  FailureContext context;
};

MoverResult ClassifyMoverFailure(const Transport::Status& status,
                                 const FailureContext& context,
                                 FailurePolicy policy);

// A small dependency-free transaction fixture mirrors the production rule:
// particle state is staged locally and is committed only after all associated
// coupling/source records have validated.
class ParticleTransaction {
 public:
  explicit ParticleTransaction(const Transport::ParticleState& initial)
      : initial_(initial), staged_(initial) {}
  Transport::ParticleState* staged() { return &staged_; }
  const Transport::ParticleState& initial() const { return initial_; }
  Transport::Status Commit(Transport::ParticleState* destination);
  void Rollback();
  bool committed() const { return committed_; }

 private:
  Transport::ParticleState initial_;
  Transport::ParticleState staged_;
  bool committed_ = false;
};

}  // namespace NumericalExtensions
}  // namespace SEP

#endif  // SEP_UTIL_SEP_NUMERICAL_EXTENSIONS_H
