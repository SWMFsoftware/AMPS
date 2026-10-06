#ifndef SEP_CORONA_SWCME_SHOCK_FRONT_PROVIDER_H
#define SEP_CORONA_SWCME_SHOCK_FRONT_PROVIDER_H

#include "sep_corona_swcme/ambient_state.h"
#include "sep_coronal_cme/mhd_jump_solver.h"

#include <array>
#include <cstdint>
#include <functional>
#include <memory>
#include <string>
#include <vector>

namespace SEP { namespace CoronaSwcme { namespace ShockFront {

using AssetReader = std::function<Core::Result<std::string>(const std::string&)>;

enum class HistoryModel { ConstantSpeed, QuinticPulseThenConstant };
enum class Phase { CoronalHistory, SwcmeOuter };
enum class TrajectoryValidityPolicy { ReportOnly, RequireDeclaredShockCoverage };

// These statuses distinguish a physical absence of a fast shock from a
// numerical inability to reconstruct a requested jump.  Downstream is valid
// only for SolvedFastShock; no zero-filled primitive can masquerade as vacuum.
enum class FrontStatus {
  OutsideFrontSupport,
  BelowPhysicalInnerBoundary,
  AmbientUnavailable,
  NonForwardInflow,
  SubfastFront,
  SolvedFastShock,
  NumericallyUnresolvedWeakShock,
  WrongBranch,
  InvalidJump
};

// Diagnostic validity is deliberately independent of jump validity.  A
// perpendicular fast shock can have a perfectly regular conservative jump
// while the conventional de Hoffmann--Teller (HT) boost is infinite; a
// magnetic null has a hydrodynamic jump limit but no magnetic direction or
// magnetic-compression ratio.  Collapsing those cases into zero-valued
// diagnostics would give downstream clients a plausible but false vector.
enum class HtDiagnosticStatus {
  Valid,
  UpstreamMagneticNull,
  PerpendicularNoFiniteBoost,
  NumericalFailure
};

struct ShockDiagnostics {
  bool magneticCompressionValid = false;
  double magneticCompression = 0.0;
  HtDiagnosticStatus htStatus = HtDiagnosticStatus::NumericalFailure;
  // All velocity vectors use the HCI basis.  htTangentialBoostMPerS is the
  // tangential Galilean boost *from the shock frame to the HT frame*, not the
  // generally larger incident-flow speed in that frame.
  CoronalCME::Vec3 htTangentialBoostMPerS;
  CoronalCME::Vec3 upstreamVelocityHtMPerS;
  CoronalCME::Vec3 downstreamVelocityHtMPerS;
  double incidentHtSpeedMPerS = 0.0;
  double absoluteNormalFieldCosine = 0.0;
  double upstreamElectricCancellation = 0.0;
  double downstreamElectricCancellation = 0.0;
  std::string reason;
};

struct LocalJumpEvaluation {
  FrontStatus status = FrontStatus::InvalidJump;
  CoronalCME::MhdCharacteristicSpeeds characteristics;
  double inflowMPerS = 0.0;
  double fastSpeedMarginMPerS = 0.0;
  double fastMach = 0.0;
  double signedMagneticNormalCosine = 0.0;
  bool magneticDirectionValid = false;
  bool downstreamValid = false;
  CoronalCME::MhdShockSolution jump;
  ShockDiagnostics diagnostics;
  std::string reason;
};

// Evaluate the complete local shock decision used by Provider::Prepare.
// The outward normal points from the represented interior toward upstream and
// normalSpeedMPerS is positive in that direction.  Physical no-shock outcomes
// are successful evaluations with absent downstream state.  Failure to solve
// a positive super-fast candidate remains a typed numerical/model outcome and
// is never silently relabelled sub-fast.
Core::Result<LocalJumpEvaluation> EvaluateLocalJump(
    const CoronalCME::MhdPrimitiveState& upstream,
    CoronalCME::Vec3 outwardNormal,double normalSpeedMPerS,
    double gammaAdiabatic,double weakMachTolerance,
    double residualTolerance,double magneticDirectionThresholdT);

struct Configuration {
  std::string schema;
  std::string profile;
  std::string referenceEpoch;
  std::string timeScale;
  std::string coordinateFrame;
  std::string particleMode;
  AmbientDefinition ambient;
  HistoryModel historyModel = HistoryModel::QuinticPulseThenConstant;
  CoronalCME::Vec3 direction;
  double halfWidthRad = 0.0;
  double physicalInnerRadiusM = 0.0;
  double initialApexRadiusM = 0.0;
  double initialApexSpeedMPerS = 0.0;
  double finalPulseSpeedMPerS = 0.0;
  double accelerationDurationS = 0.0;
  double handoffApexRadiusM = 0.0;
  double dragGammaPerM = 0.0;
  double effectiveTrajectoryWindMPerS = 0.0;
  TrajectoryValidityPolicy validityPolicy = TrajectoryValidityPolicy::ReportOnly;
  double backgroundDtS = 0.0;
  int polarCells = 0;
  int azimuthCells = 0;
  // The topology spelling is part of the event fingerprint.  A front mesh is
  // not a presentation-only choice once its facets can become independently
  // weighted source patches: changing the tessellation changes quadrature,
  // stable labels and therefore future stochastic identities.
  std::string surfaceTopology;
  double weakMachTolerance = 0.0;
  double rhResidualTolerance = 0.0;
  double endpointRadiusM = 0.0;
  CoronalCME::Vec3 observerPositionM;
  bool requireFastShockAtObserver = false;
  std::string harmonicAssetPath;
  std::string harmonicAssetSha256;
  std::string normalizedManifest;
  std::string physicsFingerprint;
};

Core::Result<std::shared_ptr<const Configuration>> ResolveConfiguration(
    const std::string& inputBytes,const AssetReader& reader);

struct TrajectoryState {
  double timeS = 0.0;
  Phase phase = Phase::CoronalHistory;
  double apexRadiusM = 0.0;
  double apexSpeedMPerS = 0.0;
  double apexAccelerationMPerS2 = 0.0;
};

struct GeometrySample {
  std::uint64_t stableId = 0;
  CoronalCME::Vec3 positionM;
  CoronalCME::Vec3 outwardNormal;
  double normalSpeedMPerS = 0.0;
  double areaM2 = 0.0;
  bool supportEdge = false;
};

// Vertices describe only the piecewise-planar embedding used by surface
// consumers.  Plasma and shock classifications intentionally do not live at
// vertices: those fields can be discontinuous where a super-fast portion of
// the prescribed front meets a sub-fast or non-forward portion.
struct SurfaceVertex {
  std::uint64_t stableId = 0;
  CoronalCME::Vec3 positionM;
  bool supportEdge = false;
  bool apex = false;
};

// One triangle is one physical quadrature/source patch.  vertex contains
// zero-based indices into Epoch::vertices.  curvedAreaM2 is the exact area of
// the associated patch in the generating-sphere (mu,phi) chart and is the only
// area used by physics ledgers.  planarAreaM2 is the chord-triangle area; it is
// retained solely to measure geometric approximation error and must never be
// substituted into a source normalization.
struct SurfaceTriangle {
  std::uint64_t stableId = 0;
  std::array<std::uint32_t,3> vertex{{0,0,0}};
  double curvedAreaM2 = 0.0;
  double planarAreaM2 = 0.0;
};

// Complete finite-SSE instantaneous kinematics.  Direction and its derivative
// are HCI unit-vector data with d.ddot=0; the latter changes the physical front
// orientation and is not a tangential relabelling velocity.  The provider's
// first selected profile supplies zero width/direction rates, while this shared
// function keeps the complete level-set derivative testable and prevents a
// future history adapter from silently dropping either term.
struct SseKinematicState {
  double apexRadiusM = 0.0;
  double apexSpeedMPerS = 0.0;
  CoronalCME::Vec3 direction;
  CoronalCME::Vec3 directionRatePerS;
  double halfWidthRad = 0.0;
  double halfWidthRateRadPerS = 0.0;
};
Core::Result<GeometrySample> EvaluateSseRay(
    const SseKinematicState& state,CoronalCME::Vec3 rayDirection,
    std::uint64_t stableId=0);

// Exact constant-coefficient quadratic-drag state matched at (th,Rh,Vh).
// This is a prescribed front proxy, not a global momentum solution.  Keeping
// the closed form public lets tests compare all sign/zero limits against an
// independent ODE integration without duplicating the provider implementation.
Core::Result<TrajectoryState> EvaluateQuadraticDrag(
    double epochS,double handoffTimeS,double handoffRadiusM,
    double handoffSpeedMPerS,double effectiveWindMPerS,double gammaPerM);

struct ShockRecord {
  GeometrySample geometry;
  FrontStatus status = FrontStatus::InvalidJump;
  AmbientPrimitive upstream;
  CoronalCME::MhdCharacteristicSpeeds characteristics;
  double inflowMPerS = 0.0;
  double fastSpeedMarginMPerS = 0.0;
  double fastMach = 0.0;
  double signedMagneticNormalCosine = 0.0;
  bool magneticDirectionValid = false;
  bool downstreamValid = false;
  CoronalCME::MhdShockSolution jump;
  ShockDiagnostics diagnostics;
  std::string reason;
};

struct AreaLedger {
  double geometricSupportM2 = 0.0;
  double belowInnerBoundaryM2 = 0.0;
  double superfastCandidateM2 = 0.0;
  double acceptedShockM2 = 0.0;
  double numericalFailureM2 = 0.0;
};

// Deterministic surface-accounting reduction used by both production epochs
// and independent disconnected/thin-patch fixtures.  Geometric support is the
// denominator; accepted area is not silently renormalized after numerical
// failures.  Zero-area samples contribute nothing and cannot create NaNs.
AreaLedger SummarizeAreas(const std::vector<ShockRecord>& records) noexcept;

struct HandoffReceipt {
  // Exact root and matched HCI state.  The first profile is C1 in front
  // position/normal speed but intentionally permits an acceleration jump.
  double timeS = 0.0;
  double apexRadiusM = 0.0;
  double apexSpeedMPerS = 0.0;
  double coronalAccelerationMPerS2 = 0.0;
  double outerAccelerationMPerS2 = 0.0;
  double maximumSurfacePositionMismatchM = 0.0;
  double maximumNormalSpeedMismatchMPerS = 0.0;
  bool c1Matched = false;
  bool c2Matched = false;
  bool authorityReset = false;
  bool canonicalWindValid = false;
  double canonicalApexWindMPerS = 0.0;
  double effectiveTrajectoryWindMPerS = 0.0;
  double effectiveMinusCanonicalWindMPerS = 0.0;
};

struct Epoch {
  TrajectoryState trajectory;
  std::uint64_t generation = 0;
  std::uint64_t ambientGeneration = 0;
  std::string eventIdentity;
  // The finite SSE cap is a topological disk.  Rings include the exact
  // finite-support boundary and terminate in one shared apex vertex, so the
  // azimuth seam and pole are closed without inventing a rear cap.  records
  // and triangles are one-to-one and have identical stable IDs.
  std::vector<SurfaceVertex> vertices;
  std::vector<SurfaceTriangle> triangles;
  std::vector<ShockRecord> records;
  AreaLedger area;
  HandoffReceipt handoff;
  bool geometricEndpointReached = false;
  bool apexShockAccepted = false;
};

enum class VolumeCapability { AmbientReference, PhysicalDownstreamVolume };

// One provider owns the front through both phases and the canonical ambient
// independently of the phase.  Prepare constructs a private immutable epoch;
// only a completely valid candidate replaces Current().  Physical sub-fast
// records are valid epoch outcomes, while malformed ambient/RH candidates are
// retained as typed records unless the selected coverage policy requires them.
class Provider final {
 public:
  static Core::Result<std::shared_ptr<Provider>> Create(
      std::shared_ptr<const Configuration> configuration);

  Core::Result<TrajectoryState> Trajectory(double epochS) const;
  Core::Result<double> HandoffTimeS() const;
  Core::Result<HandoffReceipt> Handoff() const;
  Core::Result<double> EndpointTimeS() const;
  Core::Result<GeometrySample> EvaluateRay(
      CoronalCME::Vec3 rayDirection,double epochS,std::uint64_t stableId=0) const;
  // Reconstruct and classify an arbitrary point that is asserted to lie on
  // the finite leading surface.  The analytical surface is re-evaluated and a
  // dimensional mismatch is rejected, preventing connectivity/observer code
  // from accepting the generating sphere's rear root.
  Core::Result<ShockRecord> EvaluateFrontPoint(
      CoronalCME::Vec3 positionM,double epochS,std::uint64_t stableId=0) const;
  Core::Result<AmbientState> QueryAmbient(
      CoronalCME::Vec3 positionM,double epochS,std::uint64_t generation) const;
  Core::Status RequireCapability(VolumeCapability capability) const;
  Core::Result<std::shared_ptr<const Epoch>> Prepare(
      double epochS,std::uint64_t generation);

  std::shared_ptr<const Epoch> Current() const { return current_; }
  const Configuration& Event() const { return *configuration_; }
  const AmbientModel& Ambient() const { return *ambient_; }

 private:
  std::shared_ptr<const Configuration> configuration_;
  std::shared_ptr<const AmbientModel> ambient_;
  std::shared_ptr<const Epoch> current_;
};

const char* Name(HistoryModel) noexcept;
const char* Name(Phase) noexcept;
const char* Name(FrontStatus) noexcept;
const char* Name(HtDiagnosticStatus) noexcept;

} } } // namespace SEP::CoronaSwcme::ShockFront

#endif
