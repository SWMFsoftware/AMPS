#ifndef SEP_CORONA_SWCME_CME_EVENT_H
#define SEP_CORONA_SWCME_CME_EVENT_H

#include "sep_coronal_cme/ellipsoid_geometry.h"
#include "sep_status.h"

#include <array>
#include <cstdint>
#include <functional>
#include <map>
#include <memory>
#include <string>
#include <vector>

namespace SEP { namespace CoronaSwcme {

// BG3D-1 owns background configuration only.  No particle/source type appears
// in this API, so a zero-particle background epoch can be prepared without a
// release plan, random stream, wave feedback, or compiled species buffer.
struct AssetIdentity {
  std::string role;
  std::string path;       // provenance; excluded from the physics fingerprint
  std::string sha256;     // identity of the bytes returned by AssetReader
  std::size_t bytes = 0;
};

struct ComponentHistory {
  std::vector<CoronalCME::HermiteKnot> knots;
  Core::Result<CoronalCME::KinematicValue> At(double timeS) const;
};

struct Composition {
  double gammaAdiabatic = 0.0;
  double electronTemperatureK = 0.0;
  double protonTemperatureK = 0.0;
  double alphaTemperatureK = 0.0;
  double alphaToProtonNumberRatio = 0.0;
  bool includeElectronMass = false;
};

// Typed BG3D-1 inputs for the BG3D-2 ambient implementation. Units are SI.
// Harmonic rows remain canonical text here; BG3D-2 constructs and qualifies
// the maintained PFSS/Parker provider from this already checksummed record.
struct AmbientInput {
  double sourceSurfaceRadiusM = 0.0;
  double referenceRadiusM = 0.0;
  double electronDensityAtReferenceM3 = 0.0;
  double closedBaseElectronDensityM3 = 0.0;
  double rotationRateRadPerS = 0.0;
  double minimumMagneticFieldT = 0.0;
  double traceStepM = 0.0;
  double gradientRelativeStep = 0.0;
  int traceMaximumSteps = 0;
  int windTablePoints = 0;
  // Legacy/full-CME profiles default to field-aligned open-tube mass flux.
  // The reduced model's documented radial Parker reference instead conserves
  // rho U r^2 below the source surface as well.  This explicit switch avoids
  // changing BG3D-2 semantics while making the density closure fingerprinted.
  bool sphericalOpenMassFlux = false;
  std::string harmonics;
};

// Frozen closure inputs used by later material-map stages. A vector-potential
// family, rather than a spheromak requirement, owns the reference fluxes.
struct RegionalInput {
  double sheathAdmissionStartS = 0.0;
  // The geometrical fraction belongs to the future ejecta reference body.
  // BG3D-4 must not reinterpret it as a material sheath contact: a surface
  // prescribed as a fixed fraction of a shock generally has nonzero relative
  // mass flux.
  double contactApexFraction = 0.0;
  double ejectaReferenceDensityKgM3 = 0.0;
  double ejectaReferencePressurePa = 0.0;
  double axialFluxWb = 0.0;
  double poloidalFluxWb = 0.0;
  double minimumJacobian = 0.0;
  double maximumIntegratedForceRatio = 0.0;
  double maximumLocalForceRatioP99 = 0.0;
  double maximumForceWorkRatio = 0.0;
  // Contact leakage is graded with an absolute SI term plus a relative term
  // times this finite, independently frozen case scale.  It is not normalized
  // by a local shock flux, which may vanish on sub-fast patches.
  double contactFluxAbsoluteToleranceKgM2S = 0.0;
  double contactFluxRelativeTolerance = 0.0;
  double contactFluxReferenceKgM2S = 0.0;
  // Zero-inventory/startup budgets require a dimensional kg term; dividing by
  // initial mass would be singular for the selected limiting construction.
  double inventoryMassAbsoluteToleranceKg = 0.0;
  double contactNormalVelocityNumericalTolerance = 0.0;
  // Inputs below configure only the retained Level A+ relaxation diagnostic.
  // They are not parameters of the selected production piston closure.
  // Lagrangian post-shock drift L(age)*D_birth uses
  // L'=kappa+(1-kappa)exp(-age/T).  Kappa>0 prevents a singular old-cohort
  // volume while T controls the prescribed relaxation from the exact RH
  // boundary velocity.  These are physical closure parameters, not solvers.
  double sheathDriftAsymptoteFraction = 0.0;
  double sheathDriftRelaxationTimeS = 0.0;
  // Per-ray piston controls live in their own typed record below.  They must
  // not be mapped onto these legacy diagnostic scalars: that would make two
  // physically different closures share an event identity.
  std::string vectorPotentialModel;
  std::string addedHeating;
  std::string sheathStartupModel;
  std::string sheathContactModel;
  std::string sheathReferenceMapModel;
};

// Immutable input for the BG3D-4 Level-B piston/contact authority.  This is
// intentionally separate from ComponentHistory: the latter is the prescribed
// Level-A shock/front used by the legacy RH regression, whereas this record is
// an independently parameterized ejecta boundary which drives the computed
// Level-B compression.  Keeping two typed records prevents a front fit from
// being silently reinterpreted as a material piston.
//
// All values are SI and angles are in the inertial HCI frame.  The selected
// qualification contact has fixed orientation (a C-infinity special case of
// the required C2 orientation law).  Its four geometry rates are multiplied
// by the same quintic startup function, so shape evolution is globally smooth
// rather than blended independently on each ray.
struct PistonContactInput {
  bool enabled = false;
  std::string profile;
  std::string orientationModel;
  std::string attachmentPolicy;
  std::string velocityReduction;
  double startS = 0.0;
  double startupRampDurationS = 0.0;
  double initialCenterDistanceM = 0.0;
  double initialRadialSemiAxisM = 0.0;
  double initialFirstLateralSemiAxisM = 0.0;
  double initialSecondLateralSemiAxisM = 0.0;
  double centerRateMPerS = 0.0;
  double radialRateMPerS = 0.0;
  double firstLateralRateMPerS = 0.0;
  double secondLateralRateMPerS = 0.0;
  double latitudeRad = 0.0;
  double longitudeRad = 0.0;
  double lateralTiltRad = 0.0;
  double handoffApexRadiusM = 0.0;
  double handoffTransitionDurationS = 0.0;
  double outerAmbientSpeedMPerS = 0.0;
  double outerDragCoefficientPerM = 0.0;
  double minimumContactIncidence = 0.0;
  double launchMarginM = 0.0;
  double startupMachTolerance = 0.0;
  double maximumApexAccelerationMPerS2 = 0.0;
  double minimumFinalApexSpeedMPerS = 0.0;
  double maximumFinalApexSpeedMPerS = 0.0;
};

struct PistonRayInput {
  std::uint64_t id = 0;
  CoronalCME::Vec3 direction;
  double solidAngleSr = 0.0;
  // -1 denotes the physical latitude edge of this cell.  Longitude neighbors
  // are periodic and must always be present.  Explicit topology is part of
  // the event identity because future transverse-gradient diagnostics and 3-D
  // assembly must not infer neighbors from MPI ownership or array order.
  int latitudeMinus = -1;
  int latitudePlus = -1;
  int longitudeMinus = -1;
  int longitudePlus = -1;
};

struct PistonRayQuadrature {
  bool enabled = false;
  std::vector<PistonRayInput> rays;
};

// Immutable numerical and finite-domain choices for the Level-B per-ray
// piston.  These values are parsed from the strict event deck and therefore
// participate in the physics fingerprint.  They are background controls only:
// no particle population, source cadence or random-stream setting may alter
// them.  Length/time/pressure thresholds use SI units except the explicitly
// dimensionless CFL, VNR coefficients and relative disturbance/shock ratios.
struct PistonNumericsInput {
  bool enabled = false;
  int initialCells = 0;
  double initialBufferM = 0.0;
  int sourceTablePoints = 0;
  double trajectoryMaximumStepS = 0.0;
  // Absolute event-time cadence for deterministic buffer maintenance.  It is
  // independent of how often an application requests/publishes epochs.
  double bufferCheckIntervalS = 0.0;
  int minimumBufferCells = 0;
  int appendCells = 0;
  double disturbanceRelativeThreshold = 0.0;
  double cfl = 0.0;
  std::string artificialViscosity;
  double quadraticViscosity = 0.0;
  double linearViscosity = 0.0;
  double shockThreshold = 0.0;
  std::string wellBalancedSources;
};

struct HandoffLaw {
  double transitionBeginS = 0.0;
  double transitionEndS = 0.0;
  double ambientSpeedMPerS = 0.0;
  double dragCoefficientPerM = 0.0;
};

struct EventSupport {
  double startS = 0.0;
  double endS = 0.0;
  double solarRadiusM = 0.0;
  double firstValidPlasmaRadiusM = 0.0;
  double coverageRadiusM = 0.0;
  double rootToleranceS = 0.0;
  double maximumApexAccelerationMPerS2 = 0.0;
};

struct RadialExtent {
  double minimumRadiusM = 0.0;
  double maximumRadiusM = 0.0;
  bool intersectsSolarSurface = false;
};

// Exact spatial extrema for a fixed-orientation ellipsoid whose centre is on
// its radial principal axis. Interior extrema are retained; centre +/- radial
// semiaxis alone is not sufficient for a triaxial surface.
Core::Result<RadialExtent> EvaluateRadialExtent(
    const CoronalCME::EllipsoidKinematics& kinematics,double solarRadiusM);

enum class EvolutionPhase { CoronalHistory, HandoffTransition, SwcmeOuter };

struct EventKinematics {
  double timeS = 0.0;
  double handoffWeight = 0.0;
  EvolutionPhase phase = EvolutionPhase::CoronalHistory;
  CoronalCME::EllipsoidKinematics ellipsoid;
  CoronalCME::KinematicValue apexRadiusM;
};

struct EventConfiguration {
  std::string profile;
  std::string coordinateFrame;
  std::string ambientProfile;
  std::string sheathModel;
  std::string ejectaModel;
  std::string outerEvolution;
  std::string attachmentPolicy;
  CoronalCME::RadialPrincipalBasis basis;
  std::array<ComponentHistory,4> components; // center, radial, lateral 1/2
  Composition composition;
  AmbientInput ambient;
  RegionalInput regional;
  PistonContactInput pistonContact;
  PistonRayQuadrature pistonRays;
  PistonNumericsInput pistonNumerics;
  HandoffLaw handoff;
  EventSupport support;
  std::vector<AssetIdentity> assets;
  std::map<std::string,std::string> normalizedAssignments;
  std::string manifest;
  std::string physicsFingerprint;

  Core::Result<EventKinematics> At(double timeS) const;
};

// Applications acquire bytes; the shared resolver owns syntax, content hashes,
// physical validation, and immutable identity.  A reader must return the exact
// requested bytes or a typed failure.  No MPI or filesystem type crosses here.
using AssetReader = std::function<Core::Result<std::string>(const std::string&)>;

// Parse a strict `key=value` background record. Blank lines and lines whose
// first non-space character is `#` are ignored. Keys are unique and closed;
// unknown, duplicate, missing, nonfinite, or unsupported values fail before an
// event object is returned. Asset paths are provenance while SHA-256 content
// checksums participate in the physics identity.
Core::Result<std::shared_ptr<const EventConfiguration>> ResolveEventConfiguration(
    const std::string& inputBytes, const AssetReader& reader);

// Re-run continuous support checks on a typed record.  This is public so a
// coupled host cannot bypass the same history/coverage/closure gate used by the
// text resolver.
Core::Status ValidateEventConfiguration(const EventConfiguration& configuration);

const char* Name(EvolutionPhase phase) noexcept;

} } // namespace SEP::CoronaSwcme

#endif
