#ifndef SEP_UTIL_SEP_TURBULENCE_CORE_H
#define SEP_UTIL_SEP_TURBULENCE_CORE_H

#include "sep_common_header_path.h"
#include SRCSEP_SEP_COMMON_HEADER(sep_transport_common.h)

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace SEP {
namespace Turbulence {

// Source and representation are independent of the particle mover.  The two
// self-consistent values deliberately encode their authoritative storage so a
// run cannot evolve integrated and spectral arrays as competing owners.
enum class Source {
  Prescribed,
  SelfConsistentIntegrated,
  SelfConsistentSpectral,
  SwmfReadOnly,
  SwmfInitialThenEvolveLocal
};
enum class Representation { Integrated, Spectral };
enum class BoundaryPolicy {
  SpecifiedIncomingEnergy,
  SpecifiedIncomingFlux,
  TransparentOutflow,
  FixedReservoir
};
enum class CouplingPolicy { Disabled, StreamingEnergyExchange };
// Lie splitting is deliberately named instead of being an undocumented call
// order.  The production scheme remains first-order in time, but every active
// operator is advanced with the same substep selected from all local limits.
// A future Strang implementation can therefore be added as a new enum value
// without silently changing restart or convergence semantics.
enum class OperatorSplitting { LieFirstOrder };
enum class OperatorPhase {
  Ready,
  SourcesApplied,
  CouplingApplied,
  AdvectionApplied,
  ReflectionApplied,
  CascadeApplied,
  Synchronized
};

const char* SourceName(Source value);
const char* RepresentationName(Representation value);
const char* BoundaryPolicyName(BoundaryPolicy value);
const char* CouplingPolicyName(CouplingPolicy value);
bool ParseSource(const std::string& text, Source* value);
bool ParseBoundaryPolicy(const std::string& text, BoundaryPolicy* value);
bool ParseCouplingPolicy(const std::string& text, CouplingPolicy* value);

struct BoundaryCondition {
  BoundaryPolicy policy = BoundaryPolicy::FixedReservoir;

  // SpecifiedIncomingEnergy and FixedReservoir use joules.  The flux policy
  // uses watts directed into the domain.  TransparentOutflow ignores value.
  double value = 0.0;
};

struct Configuration {
  Source source = Source::SelfConsistentIntegrated;
  Representation representation = Representation::Integrated;
  CouplingPolicy coupling = CouplingPolicy::StreamingEnergyExchange;
  BoundaryCondition innerBoundary;
  BoundaryCondition outerBoundary;
  bool advectionEnabled = true;
  bool reflectionEnabled = true;
  bool cascadeEnabled = true;
  bool shockInjectionEnabled = true;
  bool periodicBoundaries = false;
  double cflSafety = 0.8;
  OperatorSplitting operatorSplitting = OperatorSplitting::LieFirstOrder;
  // Accuracy controls are dimensionless except minimumSubstepS [s].  They are
  // configuration rather than hidden epsilon values, are fingerprinted by the
  // run configuration, and bound source, reflection, and cascade changes in
  // the same manner as the advection CFL bound.
  double operatorAccuracySafety = 0.25;
  double maximumSourceFraction = 0.25;
  double maximumCascadeFraction = 0.20;
  double minimumSubstepS = 1.0e-12;
  std::uint64_t maximumSubsteps = UINT64_C(1000000);
  double reflectionCoefficient = 0.6;
  double cascadeCoefficient = 0.8;
  double perpendicularCorrelationLengthM = 1.0e7;
  double electronHeatingFraction = 0.3;
  double spectralKMinPerM = 1.0e-10;
  double spectralKMaxPerM = 1.0e-2;
  std::size_t spectralBins = 128;
  std::size_t diagnosticCadence = 100;
  double conservationRelativeTolerance = 1.0e-12;
};

Transport::Status ValidateConfiguration(const Configuration& configuration);

// The process-wide configuration is installed once at the driver boundary.
// State objects copy it at initialization/checkpoint time, so later CLI writes
// cannot mutate an in-flight or restarted turbulence generation.
Transport::Status SetActiveConfiguration(const Configuration& configuration);
const Configuration& ActiveConfiguration();
bool EvolvesLocally(Source source);

struct CellState {
  double lengthM = 0.0;
  double volumeM3 = 0.0;
  double plasmaSpeedMPerS = 0.0;
  double alfvenSpeedMPerS = 0.0;
  double dLnAlfvenSpeeddsPerM = 0.0;
  double magneticFieldT = 0.0;
  double massDensityKgPerM3 = 0.0;
  double ePlusJ = 0.0;
  double eMinusJ = 0.0;

  // Spectral storage is branch-major: [E+(k_0)..E+(k_N-1),
  // E-(k_0)..E-(k_N-1)], each in joules per logarithmic bin.
  std::vector<double> spectralEnergyJ;

  // Sources are accumulated before the driver and consumed exactly once by
  // the next locally owned update.  Signs are changes to wave energy [J].
  double pendingParticlePlusJ = 0.0;
  double pendingParticleMinusJ = 0.0;
  double pendingShockPlusJ = 0.0;
  double pendingShockMinusJ = 0.0;
};

struct EnergyLedger {
  double initialEnergyJ = 0.0;
  double innerBoundaryJ = 0.0;
  double outerBoundaryJ = 0.0;
  double shockSourceJ = 0.0;
  double reflectionTransferJ = 0.0;
  double cascadeTransferJ = 0.0;
  double physicalDissipationJ = 0.0;
  double particleExchangeJ = 0.0;
  double limiterCorrectionJ = 0.0;
  // Rejected energy is a request that could not be applied under the declared
  // positivity policy.  It is retained separately from limiterCorrectionJ so
  // global closure can account for the physical request and the applied state.
  double rejectedSourceJ = 0.0;
  double remapCorrectionJ = 0.0;
  double finalEnergyJ = 0.0;
  double closureResidualJ = 0.0;
};

// WP32 classifies every non-ordinary path.  These dispositions are evidence,
// not replacements for Status: fatal invalid states still return a non-OK
// Status, while supported regularization and conservative correction remain
// visible to diagnostics and checkpoint/restart ledgers.
enum class UpdateDisposition {
  Applied,
  PhysicalZero,
  Regularized,
  Corrected,
  Rejected,
  Fatal
};

struct OperatorEvent {
  std::string operatorName;
  UpdateDisposition disposition = UpdateDisposition::Applied;
  double requestedJ = 0.0;
  double appliedJ = 0.0;
  std::string reason;
};

struct Diagnostics {
  std::uint64_t subcycles = 0;
  std::uint64_t advectionSubcycles = 0;
  std::uint64_t reflectionSubcycles = 0;
  std::uint64_t cascadeSubcycles = 0;
  std::uint64_t sourceSubcycles = 0;
  std::uint64_t limiterActivations = 0;
  std::uint64_t rejectedUpdates = 0;
  std::uint64_t physicalZeroEvents = 0;
  std::uint64_t regularizationEvents = 0;
  // Deterministic work counters are independent of wall-clock noise and are
  // the primary WP40 regression signal.  They count model operations, not loop
  // implementation details, so decomposition changes can be compared exactly.
  std::uint64_t advectionFaceUpdates = 0;
  std::uint64_t reflectionCellUpdates = 0;
  std::uint64_t cascadeCellUpdates = 0;
  std::uint64_t spectralBinUpdates = 0;
  std::vector<OperatorEvent> events;
};

struct AdvancePlan {
  Transport::Status status;
  double requestedStepS = 0.0;
  double selectedSubstepS = 0.0;
  std::uint64_t substeps = 0;
  std::string limitingOperator;
  double advectionLimitS = 0.0;
  double sourceLimitS = 0.0;
  double reflectionLimitS = 0.0;
  double cascadeLimitS = 0.0;
};

struct State {
  Configuration configuration;
  std::vector<CellState> cells;
  double epochS = 0.0;
  std::uint64_t fieldLineGeneration = 0;
  std::uint64_t campaignSeed = 0;
  std::uint64_t completedSteps = 0;
  OperatorPhase phase = OperatorPhase::Ready;
  bool handoffCompleted = false;
  double handoffEpochS = 0.0;
  std::string sourceChecksum;
  std::string provenance;
  EnergyLedger accumulatedLedger;
};

struct DerivedCell {
  double wPlusJPerM3 = 0.0;
  double wMinusJPerM3 = 0.0;
  double deltaB2T2 = 0.0;
  double crossHelicity = 0.0;
};

struct CoefficientView {
  Transport::Status status;
  double resonantWaveNumberPerM = 0.0;
  std::size_t resonantBin = 0;
  int resonantBranch = 0;
  double dMuMuPerS = 0.0;
  double lambdaParallelM = 0.0;
  double kappaParallelM2PerS = 0.0;
};

struct StepResult {
  Transport::Status status;
  EnergyLedger ledger;
  Diagnostics diagnostics;
};

Transport::Status InitializeState(State* state);
Transport::Status HandoffImportedState(State* state, double epochS,
                                       const std::string& checksum);
Transport::Status SynchronizeDerivedRepresentation(State* state);
DerivedCell DeriveCell(const CellState& cell);
Transport::ScalarResult MaximumStableAdvectionStepS(const State& state);
// Builds one auditable limit decision for all active operators.  Infinite
// limits mean that the corresponding operator has no active local timescale;
// a finite non-positive or non-finite model input is returned as an error.
AdvancePlan PlanAdvance(const State& state, double dtS);
CoefficientView DeriveCoefficients(const State& state, std::size_t cell,
                                   double speedMPerS, double mu,
                                   double chargeC, double massKg);
StepResult Advance(State* state, double dtS);

// Conservative overlap remap treats each old cell's integrated energy as
// uniformly distributed in physical arc length.  It preserves branch and
// per-bin totals on a covered domain, independent of cell volume changes.
Transport::Status RemapConservatively(const State& oldState,
                                      const std::vector<CellState>& newGeometry,
                                      State* remapped,
                                      EnergyLedger* ledger);

// Checkpoints are deterministic text records with an explicit schema.  They
// include source/ownership, boundary conditions, spectral grid, operator phase,
// pending coupling, RNG metadata, and accumulated ledger totals.
Transport::Status SerializeCheckpoint(const State& state, std::string* text);
Transport::Status DeserializeCheckpoint(const std::string& text, State* state);
std::uint64_t StateHash(const State& state);

}  // namespace Turbulence
}  // namespace SEP

#endif  // SEP_UTIL_SEP_TURBULENCE_CORE_H
