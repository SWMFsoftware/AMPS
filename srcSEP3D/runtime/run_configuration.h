// ============================================================================
// srcSEP3D immutable run configuration
//
// This is an AMPS-independent R2 contract.  Process arguments and
// AMPS_PARAM.in are translated by a host before this object is constructed;
// neither Runtime nor the coupled entry points parse configuration text.
// ============================================================================

#ifndef SEP3D_RUNTIME_RUN_CONFIGURATION_H
#define SEP3D_RUNTIME_RUN_CONFIGURATION_H

#include "../core/sep3d_types.h"

#include <cstddef>
#include <cstdint>
#include <memory>
#include <string>
#include <vector>

namespace SEP3D {
namespace RuntimeModel {

enum class BackgroundAuthority { AnalyticParker, Swmf };
enum class TurbulenceAuthority { Prescribed, Swmf };
enum class ShockAuthority { None, Swcme };
enum class TransportModel { Parker3D, Focused3D };
// Domain presets are physical choices, not shorthand for a hidden numeric
// default.  ``Earth`` remains an input spelling retained for compatibility;
// normalization maps it to the one-AU preset before fingerprinting.
enum class DomainPreset { Solar, OneAu, Earth, Mars };
enum class OuterRadiusMode { Preset, Explicit };
enum class InnerBoundaryMode { Absorb };
enum class OuterBoundaryMode { Escape, ImportedCoverage };
enum class RefinementProfile { Linear, PowerLaw, Smoothstep };
enum class TubeRadiusMode { PhysicalConstant, ConstantAngularWidth };
enum class RunIntent { TransportOnly, ShockInjection };
enum class MissingTurbulenceMode { Fail, Ballistic };
enum class ResonanceRangeMode { Reject, PowerLawExtension };
enum class PitchAngleSchemeMode { ReflectingMilstein, ReflectingEulerMaruyama };
enum class PerpendicularDiffusionMode { None, Constant, ConstantRatio };
enum class DriftMode { None, GradientB, Curvature, GradientAndCurvature };
enum class ObserverKind {
  FixedCartesian,
  FixedHeliographic,
  MovingCartesian,
  SphericalShell,
  FieldConnected
};
enum class ObserverNormalization { RepresentedParticles, DifferentialIntensity };

const char* Name(BackgroundAuthority value);
const char* Name(TurbulenceAuthority value);
const char* Name(ShockAuthority value);
const char* Name(TransportModel value);
const char* Name(DomainPreset value);
const char* Name(OuterRadiusMode value);
const char* Name(InnerBoundaryMode value);
const char* Name(OuterBoundaryMode value);
const char* Name(RefinementProfile value);
const char* Name(TubeRadiusMode value);
const char* Name(RunIntent value);
const char* Name(MissingTurbulenceMode value);
const char* Name(ResonanceRangeMode value);
const char* Name(PitchAngleSchemeMode value);
const char* Name(PerpendicularDiffusionMode value);
const char* Name(DriftMode value);
const char* Name(ObserverKind value);
const char* Name(ObserverNormalization value);

// C02 typed physical groups.  These records deliberately contain SI values
// only.  The file parser converts unit-bearing text into these records, while
// an SWMF host may populate the same records directly without linking parser
// code into the coupler.
struct ParkerPhysicsOptions {
  double sourceRadiusM = 20.0 * Core::Const::R_sun;
  double sourceLongitudeRad = 0.0;
  double sourceColatitudeRad = 0.5 * Core::Const::kPi;
  double referenceRadiusM = Core::Const::AU;
  double radialFieldAtReferenceT = 3.0e-9;
  double solarRotationRateRadPerS = Core::Const::Omega_sun;
  double solarWindSpeedMPerS = Core::Const::V_sw_default;
  int magneticPolarity = 1;
  double numberDensityAtReferenceM3 = 5.0e6;
  double temperatureK = 1.0e5;
  double validityCadenceS = 3600.0;
  std::string coordinateFrame = "HCI-like-inertial";
};

struct ShockOptions {
  // Times are relative to the host simulation clock.  A shock authority of
  // None ignores the interval; Swcme requires a finite ordered interval and a
  // maximum front radius within the frozen computational domain.
  double activeFromS = 0.0;
  double activeUntilS = 86400.0;
  double initialRadiusM = 20.0 * Core::Const::R_sun;
  double maximumRadiusM = Core::Const::AU;
  double speedMPerS = 1.0e6;
  double compressionRatio = 4.0;
};

struct SourceOptions {
  bool enabled = false;
  // The values in this record apply independently to every species in AMPS'
  // compiled SpeciesList.  physicalParticleRatePerS is therefore a
  // per-species seed-particle rate before injection efficiency and shock-patch
  // partitioning [s^-1], not a total that runtime silently divides among an
  // unknown composition.  Runtime multiplies it by injectionEfficiency and
  // SWCME relative_patch_weight exactly once for each compiled species.
  double physicalParticleRatePerS = 1.0;
  double injectionEfficiency = 1.0e-4;
  // Bounds are total kinetic energy for each particle [J].  The source adapter
  // converts them to momentum separately with that species' immutable AMPS
  // mass, which is essential for mixed electron/ion tables.
  double minimumEnergyJ = 1.0e4 * Core::Const::e;
  double maximumEnergyJ = 1.0e8 * Core::Const::e;
  double spectralIndex = 5.0;
  // Exact number of computational particles injected per active time step,
  // per compiled species, over the complete active shock surface.
  std::uint64_t samplesPerStep = 1000;
};

// Provider-neutral transport record for the standalone [swcme] section.  Key
// ownership and physical validation remain in src/models/swcme; coupled hosts
// may leave this list empty and install their own already-validated provider.
struct SwcmeAssignment {
  std::string key;
  std::string value;
  std::size_t line = 0;
};

struct SpeciesOptions {
  // Species identity is deliberately absent from the post-compile runtime
  // schema.  AMPS fixes the number, order, chemical symbol, mass, and charge
  // when SpeciesList is processed; allowing a second file to restate any of
  // those values would create two conflicting authorities.  This one value is
  // a run-wide numerical policy and is installed for *every* compiled AMPS
  // species.  Individual source plans may apply a conservative per-particle
  // correction, but all block/global base weights begin from this value.
  double macroparticleWeight = 1.0;
};

// AMPS-facing code copies its generated table into these neutral records.  The
// runtime layer can then validate the full compiled table without including
// pic.h or depending on species macros such as _H_PLUS_SPEC_.  ampsIndex is
// required to be the contiguous generated array index; symbol, mass, and
// charge are read-only facts obtained from the AMPS molecular-data API.
struct CompiledSpeciesRecord {
  int ampsIndex = -1;
  std::string symbol;
  double massKg = 0.0;
  double chargeC = 0.0;
};

struct ObserverOptions {
  std::string id;
  Core::Vec3 positionM;
  bool followsTrajectory = false;
  double cadenceS = 60.0;
  unsigned energyBins = 32;
  unsigned pitchAngleBins = 24;
  std::string products = "flux,spectrum";
  ObserverKind kind = ObserverKind::FixedCartesian;
  ObserverNormalization normalization =
      ObserverNormalization::DifferentialIntensity;
  Core::Vec3 velocityMPerS;
  double collectionRadiusM = 0.01 * Core::Const::AU;
  double shellRadiusM = Core::Const::AU;
  double minimumEnergyJ = 1.0e4 * Core::Const::e;
  double maximumEnergyJ = 1.0e8 * Core::Const::e;
  double minimumMu = -1.0;
  double maximumMu = 1.0;
  // `allCompiledSpecies` is resolved only at the AMPS boundary, where the
  // generated SpeciesList count is available.  In the neutral output layer an
  // empty accepted-species vector deliberately means "accept every species";
  // the Boolean distinguishes that valid wildcard from a missing/empty
  // explicit numeric list during configuration validation.
  bool allCompiledSpecies = false;
  std::vector<int> species = {0};
};

struct MemoryModelOptions {
  // Coefficients are explicit because AMPS object sizes are build dependent.
  // They are conservative defaults for dry-run planning, not universal ABI
  // constants.  Native validation may replace them with calibrated values.
  std::size_t baseCellBytes = 256;
  std::size_t baseNodeBytes = 128;
  std::size_t blockStructureBytes = 4096;
  std::size_t communicationBytesPerBlock = 2048;
  std::size_t particleBytes = 160;
  double particlesPerCell = 2.0;
  double haloFraction = 0.20;
  double safetyMarginFraction = 0.25;
};

// Mutable input record used only while the host resolves configuration.  The
// successful factory copies it into a RunConfiguration3D exposed solely
// through const accessors.  Output formatting fields are deliberately present
// so CFG3D03/LIFE3D03 can prove that they do not alter the physics fingerprint.
struct RunConfiguration3DOptions {
  // Human-authored files set this through run.schema_version.  Version 1 is
  // retained for existing campaigns; version 2 additionally requires an
  // explicit finite Parker centreline definition.  Version 3 is the complete
  // standalone initialization contract: every application value, observer,
  // output path, and canonical SWCME3D parameter must be explicit.  Parser-
  // free programmatic SWMF hosts may keep the default and install equivalent
  // validated providers through the typed interface.
  unsigned inputSchemaVersion = 1;
  BackgroundAuthority background = BackgroundAuthority::AnalyticParker;
  TurbulenceAuthority turbulence = TurbulenceAuthority::Prescribed;
  ShockAuthority shock = ShockAuthority::None;
  TransportModel transport = TransportModel::Parker3D;
  DomainPreset domain = DomainPreset::OneAu;
  OuterRadiusMode outerRadiusMode = OuterRadiusMode::Preset;
  InnerBoundaryMode innerBoundary = InnerBoundaryMode::Absorb;
  OuterBoundaryMode outerBoundary = OuterBoundaryMode::Escape;
  Core::Vec3 coordinateOriginM = {0.0, 0.0, 0.0};
  std::string coordinateFrame = "HCI-like-inertial";

  // Finite Parker-line construction used by initialization, diagnostics, and
  // future mesh export.  Coordinates and length are SI.  A zero length/count
  // or a zero initial vector is a programmatic version-1 sentinel; Create()
  // replaces it transactionally with the resolved domain/source defaults.
  Core::Vec3 parkerSpiralOriginM = {0.0, 0.0, 0.0};
  Core::Vec3 parkerSpiralInitialPointM = {0.0, 0.0, 0.0};
  double parkerSpiralLengthM = 0.0;
  std::uint64_t parkerSpiralPointCount = 0;

  double innerRadiusM = 20.0 * Core::Const::R_sun;
  // ``outerRadiusM`` is normalized to a resolved SI value by Create().  The
  // input value is consulted only when outerRadiusMode is Explicit.
  double outerRadiusM = Core::Const::AU;
  double requestedTimeStepS = 1.0;
  std::uint64_t maximumTimeSteps = 100000001;
  std::uint64_t campaignSeed = 1;
  std::uint64_t backgroundCadenceSteps = 1;
  // All runtime schedules are expressed as integer global ticks.  Zero is
  // reserved for a disabled optional event (checkpoint only); physical
  // background and injection schedules always have a positive cadence.
  std::uint64_t injectionCadenceSteps = 1;

  // Phase-M mesh controls.  The radial law and optional Parker tube are
  // evaluated by one AMPS-independent implementation used by both the
  // production localResolution callback and the standalone octree emulator.
  double minimumCellSizeM = 0.01 * Core::Const::AU;
  double backgroundCellSizeM = 0.25 * Core::Const::AU;
  bool enableRadialRefinement = true;
  double solarSurfaceCellSizeM = 0.01 * Core::Const::AU;
  double solarRefinementOuterRadiusM = 0.25 * Core::Const::AU;
  RefinementProfile solarRefinementProfile = RefinementProfile::Smoothstep;
  double solarRefinementExponent = 1.0;
  bool enableTubeRefinement = false;
  double tubeLongitudeRad = 0.0;
  double tubeColatitudeRad = 0.5 * Core::Const::kPi;
  double tubeReferenceRadiusM = Core::Const::AU;
  double tubeRadiusAtReferenceM = 0.03 * Core::Const::AU;
  TubeRadiusMode tubeRadiusMode = TubeRadiusMode::ConstantAngularWidth;
  double tubeCellSizeM = 0.01 * Core::Const::AU;
  RefinementProfile tubeTransverseProfile = RefinementProfile::Smoothstep;
  double tubeTransverseExponent = 1.0;
  unsigned meshCellsPerBlockEdge = 4;
  unsigned maximumMeshLevel = 7;
  std::size_t meshBlockOverheadBytes = 1024;
  std::size_t meshMemoryBudgetBytes =
      std::size_t{4} * 1024 * 1024 * 1024;
  MemoryModelOptions memoryModel;

  // Complete C02 physical inputs.  ``intent`` prevents a transport-only
  // demonstration from being confused with an injection run that silently
  // disabled its source.
  RunIntent intent = RunIntent::TransportOnly;
  ParkerPhysicsOptions parker;
  ShockOptions shockModel;
  SourceOptions source;
  SpeciesOptions species;
  std::vector<SwcmeAssignment> swcmeAssignments;
  // Filled only by the standalone parser after canonical resolution.  These
  // strings enter the immutable physics identity; textual spellings alone do
  // not masquerade as a validated SWCME configuration.
  std::string swcmeConfigurationFingerprint;
  std::string swcmeResolvedManifest;
  std::vector<ObserverOptions> observers = [] {
    ObserverOptions observer;
    observer.id = "default";
    observer.positionM = {0.25 * Core::Const::AU, 0.0, 0.0};
    observer.followsTrajectory = false;
    observer.cadenceS = 60.0;
    observer.energyBins = 32;
    observer.pitchAngleBins = 24;
    observer.products = "flux,spectrum";
    return std::vector<ObserverOptions>{observer};
  }();

  // Phase-T scattering inputs.  Integrated wave amplitude, spectral shape,
  // finite-band behavior, and missing-data behavior are all fingerprinted.
  // SWMF supplies w+/w- amplitudes, but it still uses these declared spectral
  // bounds unless a future coupled interface publishes a resolved spectrum.
  double prescribedDeltaBOverB = 0.3;
  double turbulenceKMinPerM = 1.0e-10;
  double turbulenceKMaxPerM = 1.0e-7;
  double turbulenceSpectralIndex = 5.0 / 3.0;
  double turbulenceCorrelationLengthM = 0.03 * Core::Const::AU;
  MissingTurbulenceMode missingTurbulence = MissingTurbulenceMode::Fail;
  ResonanceRangeMode resonanceRange = ResonanceRangeMode::Reject;

  // Phase-P numerical controls. Each fraction limits the named local time
  // scale; none is an undocumented global safety factor. A selected step below
  // minimumTransportSubstepS returns StepUnderflow and is never silently
  // clamped. The stochastic scheme is fingerprinted because Milstein and
  // Euler-Maruyama do not define identical finite-step trajectories.
  double cellCrossingFraction = 0.4;
  double diffusionFraction = 0.2;
  double focusingFraction = 0.2;
  double coolingFraction = 0.2;
  double fieldVariationFraction = 0.2;
  double shockCrossingFraction = 0.5;
  double minimumTransportSubstepS = 1.0e-12;
  std::uint64_t maximumTransportSubsteps = 100000;
  PitchAngleSchemeMode pitchAngleScheme =
      PitchAngleSchemeMode::ReflectingMilstein;

  // V01 controlled extensions.  The coefficient is isotropic in the plane
  // perpendicular to B. ConstantRatio evaluates k_perp=ratio*k_parallel in
  // each immutable local background; Constant uses the explicit SI value.
  // Drift choices use signed species charge and the relativistic p*v form.
  // Current-sheet drift is intentionally absent: no sheet-geometry contract
  // exists yet, so accepting it would invent coupling physics.
  PerpendicularDiffusionMode perpendicularDiffusion =
      PerpendicularDiffusionMode::None;
  double constantKappaPerpendicularM2PerS = 0.0;
  double kappaPerpendicularToParallelRatio = 0.0;
  DriftMode drift = DriftMode::None;

  // Frozen pre-mesh storage choices.  Offsets are derived by the factory in a
  // canonical order; adapters may not append fields after Configure().
  bool storeMagneticGradient = false;
  bool storeVelocityGradient = false;
  std::size_t samplingBytesPerCell = 0;

  // Output-only controls: these affect products, not particle trajectories.
  std::uint64_t outputCadenceSteps = 1;
  std::uint64_t checkpointCadenceSteps = 0;
  std::string outputDirectory = "output";
  std::string outputPrefix = "sep3d";
  std::string initializationMeshTecplotFile =
      "sep3d-initialization-mesh.dat";
  std::string initializationParkerLineTecplotFile =
      "sep3d-initialization-parker-line.dat";
  std::string restartInputPath;
  std::string restartOutputPath = "restart/sep3d.chk";

  // Deprecated Boolean spellings remain in the typed surface for one source
  // compatibility interval. A true value is rejected with a migration error;
  // text input uses the enum-valued fields above.
  bool enablePerpendicularDiffusion = false;
  bool enableDrifts = false;
  bool enableExternalScriptBackground = false;
  bool enableSelfConsistent3DTurbulence = false;
};

// Validate the complete generated AMPS species table and every observer's
// numeric selection in one fail-closed operation.  The table must contain the
// same positive number of records reported by AMPS, indices must be exactly
// 0..N-1, and symbols must be non-empty and unique.  srcSEP3D's focused SEP
// transport requires finite positive rest mass and non-zero finite charge; a
// neutral compiled species is rejected explicitly because silently applying a
// charged-particle scattering model to it would be physically incorrect.
Core::Status ValidateCompiledSpeciesBinding(
    const RunConfiguration3DOptions& configured, int ampsSpeciesCount,
    const std::vector<CompiledSpeciesRecord>& compiled);

constexpr std::size_t kNoOffset = static_cast<std::size_t>(-1);

struct StorageLayout {
  std::size_t magneticFieldOffset = kNoOffset;       // 3 doubles
  std::size_t bulkVelocityOffset = kNoOffset;        // 3 doubles
  std::size_t numberDensityOffset = kNoOffset;       // 1 double
  std::size_t velocityDivergenceOffset = kNoOffset;  // 1 double
  std::size_t temperatureOffset = kNoOffset;         // 1 double
  std::size_t pressureOffset = kNoOffset;            // 1 double
  std::size_t alfvenSpeedOffset = kNoOffset;         // 1 double
  std::size_t divBhatOffset = kNoOffset;              // 1 double
  std::size_t focusingLengthOffset = kNoOffset;       // 1 double
  std::size_t curvatureOffset = kNoOffset;            // 3 doubles
  std::size_t fieldAlignedStrainOffset = kNoOffset;  // 1 double
  std::size_t magneticGradientOffset = kNoOffset;    // optional 9 doubles
  std::size_t velocityGradientOffset = kNoOffset;    // optional 9 doubles
  // Optional directional magnetic wave variances [T^2], along/against +B.
  // The historic member name is retained as a storage-ABI label.
  std::size_t waveEnergyOffset = kNoOffset;          // optional 2 doubles
  std::size_t cellAssociatedBytes = 0;
  std::size_t samplingBytesPerCell = 0;
  std::string fingerprint;
};

bool operator==(const StorageLayout& left, const StorageLayout& right);
bool operator!=(const StorageLayout& left, const StorageLayout& right);

class RunConfiguration3D final {
 public:
  static Core::Status Create(
      const RunConfiguration3DOptions& options,
      std::shared_ptr<const RunConfiguration3D>* configuration);

  RunConfiguration3D(const RunConfiguration3D&) = default;
  RunConfiguration3D& operator=(const RunConfiguration3D&) = delete;

  const RunConfiguration3DOptions& options() const { return options_; }
  const StorageLayout& storage_layout() const { return storageLayout_; }
  const std::string& physics_fingerprint() const { return physicsFingerprint_; }
  const std::string& resolved_manifest() const { return resolvedManifest_; }

 private:
  RunConfiguration3D(const RunConfiguration3DOptions& options,
                     const StorageLayout& layout,
                     const std::string& physicsFingerprint,
                     const std::string& resolvedManifest);

  const RunConfiguration3DOptions options_;
  const StorageLayout storageLayout_;
  const std::string physicsFingerprint_;
  const std::string resolvedManifest_;
};

}  // namespace RuntimeModel
}  // namespace SEP3D

#endif  // SEP3D_RUNTIME_RUN_CONFIGURATION_H
