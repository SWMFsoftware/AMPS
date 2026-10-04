#ifndef SEP_CORONA_SWCME_CME_EVENT_H
#define SEP_CORONA_SWCME_CME_EVENT_H

#include "sep_coronal_cme/ellipsoid_geometry.h"
#include "sep_status.h"

#include <array>
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
  // Lagrangian post-shock drift L(age)*D_birth uses
  // L'=kappa+(1-kappa)exp(-age/T).  Kappa>0 prevents a singular old-cohort
  // volume while T controls the prescribed relaxation from the exact RH
  // boundary velocity.  These are physical closure parameters, not solvers.
  double sheathDriftAsymptoteFraction = 0.0;
  double sheathDriftRelaxationTimeS = 0.0;
  std::string vectorPotentialModel;
  std::string addedHeating;
  std::string sheathStartupModel;
  std::string sheathContactModel;
  std::string sheathReferenceMapModel;
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
