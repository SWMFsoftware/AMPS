#ifndef SEP_TRANSPORT_COMMON_H
#define SEP_TRANSPORT_COMMON_H

#include "sep.h"
#include "util/sep_birth_window.h"
#include "amps_fieldline_frame.h"

#include <cstdint>
#include <array>
#include <string>

namespace SEP {
namespace Transport {
namespace PICAdapter {

enum class StochasticOperator { Parker = 0, Dmumu = 1, MeanFreePath = 2 };

// Only the standalone component-test dispatcher may select this diagnostic
// purpose, before initialization. It permits explicitly labelled numerical
// tests of retained operators; it never certifies the adopted production law.
enum class RunPurpose { Production, ControlledDiagnostic };
void SetStandaloneRunPurpose(RunPurpose purpose);
RunPurpose CurrentRunPurpose();
Status ProductionClosurePreflight();

// Schema 5 retains explicit birth/window history and adds reduced-source
// cohort energy. Older layouts cannot reconstruct that immutable birth value.
// Earlier/unknown records cannot recover omitted history exactly and are rejected.
const std::uint64_t ParticleTransportSchema = UINT64_C(0x5352435345500005);

struct ParticleContext {
  long int pointer = -1;
  PIC::ParticleBuffer::byte* data = NULL;
  PIC::FieldLine::cFieldLineSegment* segment = NULL;
  ParticleState state;
  std::uint64_t stableParticleId = 0;
  std::uint64_t snapshotGeneration = 0;
  double particleStepEpochS = 0.0;
  Birth::Record birth;
  Birth::FirstMovePlan firstMove;
  bool stepPrepared=false,consumesPendingBirth=false,usesDriverWindow=false;
  double moveDurationS=0.0,snapshotFromS=0.0,snapshotUntilS=0.0;
  double statisticalWeight = 0.0;
  OperatorRandomState randomStates[3];
  double mfpOpticalDepth = 0.0;
  std::uint64_t mfpEventIndex = 0;
};

// All three production movers enter through these helpers.  Validation and
// attachment therefore cannot diverge as the mover implementations evolve.
Status LoadParticle(long int pointer, ParticleContext* context);
// Injection sites call this once with a physical source/sequence key and the
// generated three-component momentum [kg m/s].  The resulting ID is persisted
// in the AMPS particle record and remains independent of later buffer moves.
void InitializeParticleTransportState(long int pointer,
                                      std::uint64_t sourceKey,
                                      std::uint64_t sequenceKey,
                                      const double momentumKgMPerS[3]);
// Physical source schedulers assign the ID from persisted source/line/species/
// event/macroparticle keys, independently of the sampled phase-space state.
void InitializeParticleTransportStateWithIdentity(long int pointer,
                                                  std::uint64_t stableId);
// Scheduled source births must supply a declared source window/origin. The
// two-argument overload above is only explicit fresh initial-condition state;
// it cannot reconstruct missing history in an existing record.
Status InitializeParticleTransportStateWithIdentity(long int pointer,
    std::uint64_t stableId,const Birth::Record& birth);
// Collective host boundary, called before PIC::TimeStep changes its clock and
// after joined particle work respectively. Only equal active global species
// timesteps share this first supported coupled clock; unequal dt is rejected.
Status FreezeNativeDriverWindow();
Status CloseNativeDriverWindow();
Status NativeDriverClockCapability();
bool NativeDriverWindowActive();
Status GetNativeDriverWindow(Birth::DriverWindow* window);
Status ValidateNativeDriverSnapshot(const Background::BackgroundSnapshot& snapshot);
Status MakeStepStartBirthRecord(Birth::Record* birth);
Status MakeUniformBirthRecord(double phaseOpen01,Birth::Record* birth);
// Resolve the requested host duration before coefficients/RNG consumption.
// A pending birth gets precisely its persisted remaining duration. All state
// stays private until CommitAndAttach publishes phase space, RNG and birth
// consumption together. Direct diagnostic fixtures without pending births
// may explicitly use their current PIC epoch and requested interval.
Status PrepareParticleStep(ParticleContext* context,double requestedHostDtS,
                            double* effectiveDtS);
ScalarResult ParticleSnapshotRemaining(const ParticleContext& context,
    double elapsedS,std::uint64_t sampledGeneration);
ScalarResult SelectParticleSubstep(const ParticleContext& context,double elapsedS,
    const std::vector<StepLimit>& limits,double minimumAdaptiveStepS,
    StepDiagnostics* diagnostics);
Status AdvanceAlongFieldLine(ParticleContext* context,
                             double displacementM);
// Provider-specific front geometry remains in the background adapter, not in
// a mover. supplied=false retains the legacy datum-based estimate. A supplied
// result is a deterministic proximity limit, not a stochastic crossing law.
Status QueryShockCrossingTime(const ParticleContext&,double relativeStreamingMPerS,
                             double* timeS,bool* supplied);
// Coordinate remains authoritative; preflight its canonical derived Cartesian
// view, then publish PB::X with phase/RNG/birth/list state. Invalid geometry
// returns a typed error before changing any particle bytes or attachment head.
Status CommitAndAttach(const ParticleContext& context);

// Resume/stage a private stream in the context.  Stage performs no shared PIC
// writes; CommitAndAttach publishes phase-space and stochastic history together.
Status ResumeRandomStream(const ParticleContext& context,
                          StochasticOperator operation,
                          KeyedRandomStream* random,
                          std::uint64_t* acceptedIntervals);
Status StageRandomStream(ParticleContext* context,
                         StochasticOperator operation,
                         const KeyedRandomStream& random,
                         std::uint64_t acceptedIntervals);

// Clip an attempted physical displacement to the first absorbing field-line
// boundary.  The returned coordinate is an exact segment endpoint and the
// signed path contains only in-domain work; callers scale interval time by the
// path ratio before enqueuing a BoundaryExit coupling record.
Status ClipPathToFieldLineBoundary(int fieldLineId, double startCoordinate,
                                   double attemptedDisplacementM,
                                   double* boundaryCoordinate,
                                   double* inDomainDisplacementM);

struct LocalBackgroundView {
  // All vectors and scalars are SI and belong to one immutable snapshot.
  // ``tangent`` and ``bUnit`` are kept separately because field-line geometry
  // and magnetic polarity need not share the same orientation.
  double magneticFieldT[3] = {0.0, 0.0, 0.0};
  double magneticFieldMagnitudeT = 0.0;
  double tangent[3] = {0.0, 0.0, 0.0};
  double bUnit[3] = {0.0, 0.0, 0.0};
  double plasmaVelocityMPerS[3] = {0.0, 0.0, 0.0};
  double numberDensityPerM3 = 0.0;
  double massDensityKgPerM3 = 0.0;
  // Provider-defined radial coordinate [m], not an unconditional |global x|.
  // A declared derivative fixture does not replace this frame or validate U.
  // Unknown/undefined radii remain unavailable; consumers must reject rather
  // than use zero as a physical radius or silently choose another origin.
  double heliocentricRadiusM = 0.0;
  bool heliocentricRadiusKnown = false;
  double radiusOriginM[3] = {0.0,0.0,0.0};
  // The same provider translation, in the unchanged Cartesian vector basis.
  // Radius-dependent libraries consume this vector directly; B/U are not
  // rotated. Both vector and radius are usable only when the known flag holds.
  double heliocentricPositionM[3] = {0.0,0.0,0.0};
  std::string radiusFrameIdentity;
  double dLnAbsBdsPerM = 0.0;
  double plasmaAdvectionMPerS = 0.0;
  // Physical U dot tangent above remains available to physical operators.
  // Only this geometry-relative value enters coordinate displacements and
  // segment-crossing limits. A material grid carries its bulk motion itself.
  AMPS::FieldLine::CoordinateFrame particleCoordinateFrame=AMPS::FieldLine::CoordinateFrame::Unknown;
  double coordinateAdvectionMPerS = 0.0;
  double parallelVelocityGradientPerS = 0.0;
  double velocityDivergencePerS = 0.0;
  double fieldAlignedStrainPerS = 0.0;
  std::array<double,9> velocityGradientPerS{};
  double alfvenSpeedMPerS = 0.0;
  int magneticPolarityAlongLine = 1;
  double snapshotSecondsRemaining = 0.0;
  double sampleCoordinate = 0.0;
  std::uint64_t generation = 0;
  std::string identity;
  std::string provenance;
};

typedef LocalBackgroundView LocalBackground;

// Controlled tests that replace a provider's native vertex velocity must
// explicitly declare its Cartesian derivative. Install and remove this guard
// outside the mover read phase; production never installs one.
class ScopedVelocityGradientFixture {
 public:
  explicit ScopedVelocityGradientFixture(
      const std::array<double,9>& gradientPerS);
  // Manufactured transverse tube flow U=(divU/2)*x_perp. The flow vanishes
  // on the line and G=(divU/2)*(I-tt) there; divU is exact even on a curved
  // centerline. Native fixture velocities must all vanish.
  explicit ScopedVelocityGradientFixture(double transverseDivergencePerS);
  ~ScopedVelocityGradientFixture();
  ScopedVelocityGradientFixture(const ScopedVelocityGradientFixture&)=delete;
  ScopedVelocityGradientFixture& operator=(const ScopedVelocityGradientFixture&)=delete;
 private:
  bool previousActive_;
  std::array<double,9> previousGradient_;
  bool previousTransverse_;
  double previousDivergence_;
};

// Reads geometry, flow, and density from one segment while the immutable
// BackgroundSnapshot read phase is active.  The derivative convention is the
// same for Parker and both focused movers.
Status EvaluateLocalBackgroundAt(const ParticleContext& context,
                                 double relativeArcLengthM,
                                 LocalBackgroundView* background);
Status EvaluateLocalBackground(const ParticleContext& context,
                               LocalBackgroundView* background);

extern std::uint64_t CampaignRandomSeed;

enum class CouplingEventType {
  DeterministicInterval = 0,
  ScatteringEvent = 1,
  BoundaryExit = 2
};

struct CouplingRecord {
  // This payload is self-contained by design.  FlushWaveContributions never
  // dereferences a particle-buffer handle, so an absorbing-boundary exit cannot
  // discard already completed in-domain work or leave a dangling pointer.
  std::string turbulenceStateIdentity;
  int fieldLineId = -1;
  int species = -1;
  double statisticalWeight = 0.0;
  std::uint64_t stableParticleId = 0;
  std::uint64_t snapshotGeneration = 0;
  std::uint64_t eventIndex = 0;
  std::uint64_t intervalIndex = 0;
  double dtS = 0.0;
  double preMomentumKgMPerS = 0.0;
  double postMomentumKgMPerS = 0.0;
  double midpointParallelVelocityMPerS = 0.0;
  double midpointNormalVelocityMPerS = 0.0;
  double startCoordinate = 0.0;
  double finishCoordinate = 0.0;
  double signedPathM = 0.0;
  int resonantBranch = 0;
  CouplingEventType eventType = CouplingEventType::DeterministicInterval;
  bool pitchAngleResolved = false;
  bool crossedBoundary = false;
};

// Movers enqueue coupling records in worker-local storage.  The driver calls
// FlushWaveContributions() after PIC::TimeStep(), when no particle thread is
// writing, to sort and apply them deterministically to the legacy G+/G- arrays.
void QueueWaveContribution(const CouplingRecord& record);
void FlushWaveContributions();
// Quiescent driver/checkpoint boundary only: workers must already be joined.
bool PendingWaveContributionsEmpty();

// Install exact-split lifecycle callbacks before stepping. A migration/restart
// byte clone retains history; only the explicit population-control split hook
// creates a new representative. Begin validates current/prospective identities
// across MPI ranks before any controller mutation; it must be called collectively.
void RegisterPopulationControlHooks();
void BeginPopulationControl();

}  // namespace PICAdapter
}  // namespace Transport
}  // namespace SEP

#endif  // SEP_TRANSPORT_COMMON_H
