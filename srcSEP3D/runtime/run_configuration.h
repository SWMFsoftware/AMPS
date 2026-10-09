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
#include "../core/domain_geometry.h"

#include <cstddef>
#include <cstdint>
#include <memory>
#include <string>
#include <vector>

namespace SEP3D {
namespace RuntimeModel {

// The Python option is a deliberately reserved standalone source.  Giving it
// a typed value now lets an input deck state its intended authority and fail
// with a precise "not implemented" status; it must never fall through to the
// SWMF branch or silently reuse the analytic Parker provider.
enum class BackgroundAuthority {
  // Append future values: existing persisted tags must retain their meaning.
  AnalyticParker,
  PythonInterpolator,
  Swmf,
  Swcme,
  RuntimeModel // A registered model; AMPS publication is provider-independent.
};
enum class TurbulenceAuthority { Prescribed, Swmf };
// A prescribed authority still needs a spectral closure.  Keeping this choice
// separate from authority distinguishes "who supplies the waves" from "which
// inertial-range slope describes them" and leaves the SWMF/AWSoM import path
// unchanged.
enum class PrescribedTurbulenceModel { PowerLaw, Kolmogorov, Kraichnan };
// The spectral slope and the pre-existing wave-amplitude prescription are
// independent physical choices.  ConstantDeltaBOverB reproduces the original
// local closure deltaB=a|B|.  WaveEnergyPowerLaw instead makes the total
// Alfvén-wave energy density itself authoritative at a declared reference
// radius and applies only the explicitly supplied radial exponent.
enum class PrescribedTurbulenceAmplitudeModel {
  ConstantDeltaBOverB,
  WaveEnergyPowerLaw
};
// A finite Parker centreline can either use independently reviewed Cartesian
// coordinates or be tied explicitly to the canonical SWCME launch-apex point.
// The linked mode is resolved by configuration_io.cpp after SWCME validates
// cme.launch_radius and geometry.cme_direction; no coordinate is inferred from
// a textual assignment before that canonical resolution succeeds.
enum class ParkerSpiralStartMode { Explicit, CmeLaunchPoint };
// Runtime records remain independent of SWCME headers.  configuration_io.cpp
// translates the canonical SWCME enum into this neutral value after SWCME has
// validated the complete input layer.
enum class SolarWindThermodynamicClosure { ProtonOnly, MultiSpecies };
enum class ShockAuthority { None, Swcme };
// The shock provider always supplies geometry and compression, but the seed
// spectrum need not always be the test-particle DSA spectrum implied by that
// compression.  FixedPhaseSpacePowerLaw represents an explicitly prescribed
// isotropic phase-space distribution f(p) proportional to p^(-q).  The
// distinction is essential for event studies such as 2013-04-11, whose
// published low-energy boundary uses q=5 independently of the local shock
// compression.  This selector controls spectral *shape* only; physical source
// rate remains the separately declared number-rate contract below.
enum class SourceSpectrumModel {
  LocalCompressionDsa,
  FixedPhaseSpacePowerLaw
};
// ConstantStatisticalWeight is the first reduced-front production sampler:
// every macro represents the same species-specific physical number W_s and
// the live incident flux changes only the Poisson event rate.  The second
// enumerator reserves the reviewed interface for a log-uniform momentum
// proposal whose individual AMPS weight correction will be implemented and
// qualified separately; selecting it currently produces a typed startup
// failure at the native injection boundary.
enum class SourceWeightingModel {
  ConstantStatisticalWeight,
  LogUniformMomentumImportance
};
// This selector governs only the physical rate used to establish the base
// Monte-Carlo weight. It is distinct from the momentum-spectrum model above.
// AcceptedShockIncidentFlux is the gross upstream population swept through
// accepted fast-shock faces; it contains no hidden acceleration efficiency.
enum class SourceRateNormalizationModel {
  ConfiguredConstant,
  AcceptedShockIncidentFlux
};
// The three production movers correspond to three different stochastic
// equations.  Focused3D is retained as a source-level alias for the released
// continuous-D_mumu mover; text input uses the unambiguous canonical name
// ``focused-diffusion`` in schema 4.
enum class TransportModel {
  Parker3D,
  FocusedDiffusion3D,
  Focused3D = FocusedDiffusion3D,
  FocusedScattering3D
};
// Top-level AMPS mover family selected by the input deck.  ``transport``
// remains the concrete physics-integrator selector because the focused family
// already has continuous-D_mumu and discrete-scattering implementations.  The
// two selectors are validated together; this enum never guesses which focused
// flavor the user intended.
enum class ParticleMoverFamily { Parker, FocusedTransport };
// Domain presets are physical choices, not shorthand for a hidden numeric
// default.  ``Earth`` remains an input spelling retained for compatibility;
// normalization maps it to the one-AU preset before fingerprinting.
enum class DomainPreset { Solar, OneAu, Earth, Mars };
enum class OuterRadiusMode { Preset, Explicit, FieldLineEndpoint };
enum class DomainBoxGeometry {
  SunCenteredCube, FieldLineCornerCube, FieldLineXYCornerCube
};
enum class SolarRefinementAnchor { SourceShell, Photosphere };
enum class InnerBoundaryMode { Absorb };
enum class OuterBoundaryMode { Escape, ImportedCoverage };
enum class RefinementProfile { Linear, PowerLaw, Smoothstep };
enum class TubeRadiusMode { PhysicalConstant, ConstantAngularWidth };
// AMPS deactivates complete leaf blocks rather than individual finite-volume
// cells.  FullDomain preserves the historical enclosing cube; ParkerTube keeps
// only blocks conservatively intersecting a configured transport corridor.
enum class ActiveRegionMode { FullDomain, ParkerTube };

// This identifier is part of the frozen physics fingerprint because mask
// semantics determine which physical cells can contain fields and particles.
// Keep it independent of mesh_model.h: the runtime contract is an L0 input to
// the mesh model and must not acquire a reverse include dependency.
constexpr const char* kActiveRegionAlgorithmName =
    "finite-parker-capsule-topological-solar-boundary-v3";
constexpr const char* kSolarSphereActiveRegionAlgorithmName =
    "finite-parker-capsule-topological-solar-sphere-corner-v4";
enum class PopulationControlMode { Off, SplitMerge };
// Coefficient choices are deliberately independent of mover selection.  A
// compatibility check in RunConfiguration3D::Create rejects circular or unused
// combinations instead of silently selecting a coefficient for the user.
enum class SpatialDiffusionModel {
  CorrelationMeanFreePath,
  // Canonical schema-4 name.  The older enumerator remains an equal-valued
  // source alias so typed hosts compiled against the first schema-4 draft do
  // not change ABI or behavior.
  MeanFreePath = CorrelationMeanFreePath,
  PitchAngleIntegral,
  // Schema 5 selects the shared, model-specific library through this distinct
  // value. It is deliberately not an alias for either legacy path: mapping an
  // old name to a scientifically different closure would let an unchanged
  // input deck acquire different transport physics.
  ParallelDiffusionLibrary
};
enum class PitchAngleDiffusionModel {
  Jokipii1966,
  Florinskiy,
  Constant
};
// RadialRigidityPowerLaw is the species-general form
//   lambda=lambda_ref*(r/r_ref)^a*(R/R_ref)^b,
// where rigidity R=p*c/|q| is expressed in volts.  For protons the published
// M-FLAMPA law (pc/1 GeV)^(1/3) is recovered by R_ref=1 GV and b=1/3.
enum class MeanFreePathModel {
  Correlation,
  Constant,
  RadialRigidityPowerLaw
};
enum class FocusedScatteringFrame {
  PlasmaFrameIsotropic,
  AlfvenWaveFrameIsotropic
};
// Propagation observes the canonical shock in a real AMPS lifecycle without
// creating source patches or particles. It is distinct from particle transport.
enum class RunIntent { TransportOnly, ShockInjection, ShockPropagation };
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
// Energy-channel spacing is part of the observer definition rather than a
// hard-coded sampling implementation detail. Logarithmic channels resolve
// power-law SEP spectra; linear channels support narrow-band diagnostics.
enum class EnergyChannelSpacing { Logarithmic, Linear };

const char* Name(BackgroundAuthority value);
const char* Name(TurbulenceAuthority value);
const char* Name(PrescribedTurbulenceModel value);
const char* Name(PrescribedTurbulenceAmplitudeModel value);
const char* Name(ParkerSpiralStartMode value);
const char* Name(SolarWindThermodynamicClosure value);
const char* Name(ShockAuthority value);
const char* Name(SourceSpectrumModel value);
const char* Name(SourceWeightingModel value);
const char* Name(SourceRateNormalizationModel value);
const char* Name(TransportModel value);
const char* Name(ParticleMoverFamily value);
const char* Name(DomainPreset value);
const char* Name(OuterRadiusMode value);
const char* Name(DomainBoxGeometry value);
const char* Name(SolarRefinementAnchor value);
const char* Name(InnerBoundaryMode value);
const char* Name(OuterBoundaryMode value);
const char* Name(RefinementProfile value);
const char* Name(TubeRadiusMode value);
const char* Name(ActiveRegionMode value);
const char* Name(PopulationControlMode value);
const char* Name(SpatialDiffusionModel value);
const char* Name(PitchAngleDiffusionModel value);
const char* Name(MeanFreePathModel value);
const char* Name(FocusedScatteringFrame value);
const char* Name(RunIntent value);
const char* Name(MissingTurbulenceMode value);
const char* Name(ResonanceRangeMode value);
const char* Name(PitchAngleSchemeMode value);
const char* Name(PerpendicularDiffusionMode value);
const char* Name(DriftMode value);
const char* Name(ObserverKind value);
const char* Name(ObserverNormalization value);
const char* Name(EnergyChannelSpacing value);

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
  double densityReferenceRadiusM = Core::Const::AU;
  double temperatureK = 1.0e5;
  // These values complete the canonical SWCME ambient thermodynamic state.
  // numberDensityAtReferenceM3 is the electron density used to normalize the
  // Leblanc profile at densityReferenceRadiusM. This radius is deliberately
  // separate from the magnetic Br reference: SWCME declares density at one AU
  // even when a user chooses another convenient radius for Br. ProtonOnly
  // reproduces the historical p=n_p k_B T_p closure. MultiSpecies enforces
  // charge neutrality and adds electron/alpha pressure and alpha mass exactly
  // as SWCME does.
  double adiabaticIndex = 5.0 / 3.0;
  SolarWindThermodynamicClosure thermodynamicClosure =
      SolarWindThermodynamicClosure::ProtonOnly;
  double alphaToProtonRatio = 0.0;
  double electronTemperatureK = 1.0e5;
  double alphaTemperatureK = 1.0e5;
  // SWCME defines its public one-AU magnetic magnitude at this reference
  // sin(colatitude).  The local 3-D Parker pitch still uses the actual point.
  double referenceSinColatitude = 1.0;
  Core::Vec3 rotationAxis = {0.0, 0.0, 1.0};
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
  // LocalCompressionDsa obtains q=3r/(r-1) from every canonical SWCME shock
  // patch.  FixedPhaseSpacePowerLaw instead uses fixedPhaseSpacePowerIndex for
  // f(p) proportional to p^(-q).  A zero q is the mandatory inactive sentinel
  // in local-DSA mode, preventing a dormant value from silently becoming
  // physics after a later mode edit.
  SourceSpectrumModel spectrumModel =
      SourceSpectrumModel::LocalCompressionDsa;
  double fixedPhaseSpacePowerIndex = 0.0;
  SourceWeightingModel weightingModel =
      SourceWeightingModel::ConstantStatisticalWeight;
  // Legacy schema-1/2 dN/dp exponent.  Schema 3 and later use the unambiguous
  // spectrumModel/fixedPhaseSpacePowerIndex contract above and normalize this
  // old field to zero before fingerprinting.
  double spectralIndex = 5.0;
  // Exact number of computational particles injected per active time step,
  // per compiled species, over the complete active shock surface.
  std::uint64_t samplesPerStep = 1000;
  // A Poisson process is intentionally unbounded.  This input is therefore
  // a fail-closed memory/runaway guard, never a cap: reaching it aborts the
  // step instead of biasing the represented source or changing W_s.
  std::uint64_t maximumMacroparticlesPerSpeciesPerStep = 1000000;
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

// Derived normalization for one post-compile AMPS species. The input file
// does not restate this identity: index/symbol/mass/charge come from the
// generated SpeciesList, while rate and weight are calculated after the mesh
// and reduced front are available. Retaining the result in the immutable run
// record makes native output/restart provenance agree with the values actually
// installed in AMPS' global and block-local arrays.
struct SpeciesParticleNormalization {
  int ampsIndex = -1;
  std::string symbol;
  double physicalSourceRatePerS = 0.0;
  double macroparticleWeight = 0.0;
};

struct ParticleNumericsOptions {
  bool deriveFromMeshAndShock = false;
  SourceRateNormalizationModel sourceRateModel =
      SourceRateNormalizationModel::ConfiguredConstant;
  double maximumParticleSpeedMPerS = 0.0;
  double timeStepMarginFactor = 0.0;
  // No physical default exists for this radius. It remains the zero inactive
  // sentinel unless the selected derived normalization requires it.
  double sourceNormalizationRadiusM = 0.0;
  // Zero means the AMPS mesh has not yet been allocated. The final immutable
  // configuration stores the MPI-reduced minimum characteristic cell size so
  // the installed dt can be rechecked from its defining equation.
  double resolvedMinimumCellSizeM = 0.0;
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
  EnergyChannelSpacing energyChannelSpacing =
      EnergyChannelSpacing::Logarithmic;
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
  // output path, and canonical SWCME3D parameter must be explicit. Version 4
  // adds the active Parker corridor, AMPS particle population control, and the
  // complete mover/coefficient selection contract. Version 5 adds the strict
  // [parallel_diffusion] binding. Its values remain strings here because the
  // shared library is the sole authority that parses each model's schema;
  // srcSEP3D must not duplicate or weaken those model-specific readers. Parser-
  // free programmatic SWMF hosts may keep the default and install equivalent
  // validated providers through the typed interface.
  unsigned inputSchemaVersion = 1;
  BackgroundAuthority background = BackgroundAuthority::AnalyticParker;
  // Used only by runtime-model. Register its factory before parsing/initializing.
  // Included in physics identity; a runtime source cannot be selected by an
  // unrecorded callback or silently substituted for a built-in authority.
  std::string backgroundModelId;
  // Checksummed reduced-model event resolved relative to the application deck
  // before immutable construction.  It is active only for the registered
  // sep-corona-swcme shock-front runtime model and enters the physics identity.
  std::string backgroundModelAssetPath;
  // Shared-section input may carry the complete reduced-model assignment
  // layer inline. Exactly one of this text or backgroundModelAssetPath is
  // active. Relative PFSS/magnetic assets resolve from the recorded directory.
  std::string backgroundModelInlineConfiguration;
  std::string backgroundModelAssetDirectory;
  TurbulenceAuthority turbulence = TurbulenceAuthority::Prescribed;
  PrescribedTurbulenceModel prescribedTurbulenceModel =
      PrescribedTurbulenceModel::Kolmogorov;
  PrescribedTurbulenceAmplitudeModel prescribedTurbulenceAmplitudeModel =
      PrescribedTurbulenceAmplitudeModel::ConstantDeltaBOverB;
  ShockAuthority shock = ShockAuthority::None;
  TransportModel transport = TransportModel::Parker3D;
  ParticleMoverFamily particleMover = ParticleMoverFamily::Parker;
  // File parsing sets this bit only for an explicit run.particle_mover key.
  // It is intentionally excluded from the physics identity: after Create()
  // normalization the selected family itself is the complete physical fact.
  // Parser-free legacy callers may omit it and have the family derived from
  // their already explicit concrete transport selector.
  bool particleMoverExplicitlySelected = false;
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
  ParkerSpiralStartMode parkerSpiralStartMode =
      ParkerSpiralStartMode::Explicit;
  // Parser-derived canonical SWCME launch-apex location.  It is populated
  // only for CmeLaunchPoint mode and then revalidated by Create().  These
  // fields let the AMPS-independent immutable factory enforce the linkage
  // without depending on SWCME headers or reparsing raw text.
  bool cmeLaunchPointResolved = false;
  Core::Vec3 cmeLaunchPointM = {0.0, 0.0, 0.0};
  double parkerSpiralLengthM = 0.0;
  std::uint64_t parkerSpiralPointCount = 0;

  // Parker/CME source and custom-transport cutoff.  This is not the solid
  // solar radius: production registers a separate AMPS internal sphere at the
  // invariant Core::Const::R_sun and validation forbids this shell below it.
  double innerRadiusM = 20.0 * Core::Const::R_sun;
  // ``outerRadiusM`` is normalized to a resolved SI value by Create().  The
  // input value is consulted only when outerRadiusMode is Explicit.
  double outerRadiusM = Core::Const::AU;
  // Opt-in corner layout: a full solar neighbourhood remains inside the box.
  // Zero direction components select a corner automatically from the line;
  // +/-1 select the direction into the domain on each Cartesian axis.
  // FieldLineXYCornerCube uses x/y selections and centers z on the Sun;
  // its z direction must be zero because it has no selected z corner face.
  DomainBoxGeometry domainBoxGeometry = DomainBoxGeometry::SunCenteredCube;
  Core::Vec3 domainCornerDirection = {0.0, 0.0, 0.0};
  double domainCornerMarginM = 0.0;
  // Endpoint is heliocentric distance, not arc length. With endpoint mode,
  // outerRadiusM and parkerSpiralLengthM are derived from this single value.
  double parkerSpiralEndRadiusM = 0.0;
  double requestedTimeStepS = 1.0;
  std::uint64_t maximumTimeSteps = 100000001;
  // Zero retains the step budget. A positive target stops after the first
  // completed native tick whose front reaches it; the bracketing row is kept.
  double stopShockRadiusM = 0.0;
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
  SolarRefinementAnchor solarRefinementAnchor =
      SolarRefinementAnchor::SourceShell;
  bool enableTubeRefinement = false;
  double tubeLongitudeRad = 0.0;
  double tubeColatitudeRad = 0.5 * Core::Const::kPi;
  double tubeReferenceRadiusM = Core::Const::AU;
  double tubeRadiusAtReferenceM = 0.03 * Core::Const::AU;
  TubeRadiusMode tubeRadiusMode = TubeRadiusMode::ConstantAngularWidth;
  double tubeCellSizeM = 0.01 * Core::Const::AU;
  RefinementProfile tubeTransverseProfile = RefinementProfile::Smoothstep;
  double tubeTransverseExponent = 1.0;
  // Optional computational-domain pruning.  The radius is evaluated with the
  // same two physical laws as mesh refinement but remains an independent,
  // normally wider corridor.  bufferBlocks retains complete neighboring leaf
  // blocks for coefficient stencils and particle crossings; it is a numerical
  // halo count, not an addition to the declared physical radius.
  ActiveRegionMode activeRegion = ActiveRegionMode::FullDomain;
  double activeTubeReferenceRadiusM = Core::Const::AU;
  double activeTubeRadiusAtReferenceM = 0.0;
  TubeRadiusMode activeTubeRadiusMode =
      TubeRadiusMode::ConstantAngularWidth;
  unsigned activeTubeBufferBlocks = 0;
  // Union of the finite corridor and a complete spherical neighbourhood.
  // Zero disables this addition. The solid photosphere is still excluded by
  // AMPS' internal boundary; this is an allocation radius, not another Sun.
  double activeSolarSphereRadiusM = 0.0;
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
  ParticleNumericsOptions particleNumerics;
  std::vector<SpeciesParticleNormalization> speciesParticleNormalizations;
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
  // Used only by WaveEnergyPowerLaw.  The value is the total directional sum
  // w_+ + w_- [J m^-3] at turbulenceReferenceRadiusM; the radial law is
  // w(r)=w_ref*(r_ref/r)^p.  Complete input sets the inactive normalization to
  // zero so no reviewed number is silently ignored by the selected model.
  double turbulenceWaveEnergyAtReferenceJPerM3 = 0.0;
  double turbulenceWaveEnergyRadialExponent = 0.0;
  // sigma_c=(deltaB_+^2-deltaB_-^2)/(deltaB_+^2+deltaB_-^2).  Requiring this
  // physical imbalance in complete input avoids an undocumented 50/50 split.
  double turbulenceNormalizedCrossHelicity = 0.0;
  double turbulenceReferenceRadiusM = Core::Const::AU;
  double turbulenceKMinPerM = 1.0e-10;
  double turbulenceKMaxPerM = 1.0e-7;
  double turbulenceKMinRadialExponent = 2.0;
  double turbulenceKMaxRadialExponent = 2.0;
  double turbulenceSpectralIndex = 5.0 / 3.0;
  double turbulenceCorrelationLengthM = 0.03 * Core::Const::AU;
  double turbulenceCorrelationLengthRadialExponent = 1.0;
  double turbulenceValidityCadenceS = 60.0;
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

  // Parallel transport coefficients.  The legacy schema-1..3 values reproduce
  // the released implementation: Parker kappa is derived from the correlation
  // mean free path and focused diffusion uses Jokipii-1966 D_mumu.  Schema 4
  // requires every selector and inactive numeric value explicitly.
  SpatialDiffusionModel spatialDiffusionModel =
      SpatialDiffusionModel::CorrelationMeanFreePath;

  struct ParallelDiffusionAssignment {
    // Names and numeric text are preserved exactly as written. Library keys
    // are case-sensitive and encode SI units (for example rigidity0_V), so the
    // normal srcSEP3D lower-casing rule must not be applied to this section.
    std::string name;
    std::string value;
    std::size_t line = 0;
  };
  std::string parallelDiffusionModelId;
  std::size_t parallelDiffusionModelLine = 0;
  std::vector<ParallelDiffusionAssignment> parallelDiffusionParameters;
  // Derived only after the shared parser validates the complete candidate.
  // It is included in physics/restart identity and never accepted from input.
  std::string parallelDiffusionConfigurationFingerprint;
  PitchAngleDiffusionModel pitchAngleDiffusionModel =
      PitchAngleDiffusionModel::Jokipii1966;
  MeanFreePathModel meanFreePathModel = MeanFreePathModel::Correlation;
  double constantDmumuPerS = 0.0;
  double constantMeanFreePathM = 0.0;
  // Parameters used only by RadialRigidityPowerLaw.  Schema 4 requires all
  // five values to be written explicitly and requires zero sentinels when a
  // different model is active; this prevents a later selector edit from
  // reviving an invisible default.
  double meanFreePathReferenceM = 0.0;
  double meanFreePathReferenceRadiusM = 0.0;
  double meanFreePathReferenceRigidityV = 0.0;
  double meanFreePathRadialExponent = 0.0;
  double meanFreePathRigidityExponent = 0.0;
  double spatialQuadratureAbsoluteToleranceM2PerS = 0.0;
  double spatialQuadratureRelativeTolerance = 1.0e-6;
  unsigned spatialQuadratureMaximumRecursion = 20;
  FocusedScatteringFrame focusedScatteringFrame =
      FocusedScatteringFrame::PlasmaFrameIsotropic;
  std::uint64_t maximumScatteringEventsPerSubstep = 100000;

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

  // AMPS population control is applied after shock injection and before
  // observers/checkpoints at the same joined timestep boundary. Limits are
  // per active AMR cell and per compiled species, which matches the granularity
  // of AMPS' linked-list split/merge primitives and avoids transporting
  // particles between unrelated spatial cells merely to meet a global count.
  PopulationControlMode populationControl = PopulationControlMode::Off;
  unsigned minimumParticlesPerCellPerSpecies = 0;
  unsigned targetParticlesPerCellPerSpecies = 0;
  unsigned maximumParticlesPerCellPerSpecies = 0;
  std::uint64_t populationControlCadenceSteps = 1;

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
  // Base path for the native AMPS data-bearing Tecplot product. A
  // multi-species executable writes one file per AMPS species by inserting a
  // deterministic `.species-N` suffix before the extension.
  std::string initializationDataTecplotFile =
      "sep3d-initialization-data.dat";
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
  // Mandatory directional magnetic wave variances [T^2], along/against +B.
  // The historic member name is retained as a storage-ABI label.  Both
  // prescribed and imported authorities reserve these two values so the
  // initialization product always contains the wave state used by scattering.
  std::size_t waveEnergyOffset = kNoOffset;          // 2 doubles
  std::size_t cellAssociatedBytes = 0;
  std::size_t samplingBytesPerCell = 0;
  std::string fingerprint;
};

bool operator==(const StorageLayout& left, const StorageLayout& right);
bool operator!=(const StorageLayout& left, const StorageLayout& right);

// Resolve the Cartesian root independently of AMPS and mesh_model.h. Callers
// pass normalized options; the immutable factory uses this same operation for
// maximum-level feasibility, and mesh/CLI/native initialization use its bounds.
inline Core::Status ResolveDomainBoundsM(
    const RunConfiguration3DOptions& options,
    Core::Vec3* minimumM, Core::Vec3* maximumM) {
  if (options.domainBoxGeometry != DomainBoxGeometry::SunCenteredCube &&
      options.domainBoxGeometry != DomainBoxGeometry::FieldLineCornerCube &&
      options.domainBoxGeometry != DomainBoxGeometry::FieldLineXYCornerCube)
    return Core::Status(Core::StatusCode::InvalidInput, "unknown domain box geometry");
  Core::DomainGeometryParameters p;
  p.cornerCube = options.domainBoxGeometry != DomainBoxGeometry::SunCenteredCube;
  p.centerZ = options.domainBoxGeometry == DomainBoxGeometry::FieldLineXYCornerCube;
  p.originM = options.coordinateOriginM;
  p.cornerDirection = options.domainCornerDirection;
  p.endpointRadiusM = options.outerRadiusM;
  p.solarSphereRadiusM = options.activeSolarSphereRadiusM;
  p.cornerMarginM = options.domainCornerMarginM;
  p.parker.sourceRadiusM = options.innerRadiusM;
  p.parker.sourceLongitudeRad = options.tubeLongitudeRad;
  p.parker.sourceColatitudeRad = options.tubeColatitudeRad;
  p.parker.solarWindSpeedMPerS = options.parker.solarWindSpeedMPerS;
  p.parker.solarRotationRateRadPerS = options.parker.solarRotationRateRadPerS;
  p.parker.rotationAxis = options.parker.rotationAxis;
  const double active = options.activeTubeRadiusAtReferenceM *
      (options.activeTubeRadiusMode == TubeRadiusMode::PhysicalConstant
           ? 1.0 : options.outerRadiusM / options.activeTubeReferenceRadiusM);
  const double refined = options.enableTubeRefinement
      ? options.tubeRadiusAtReferenceM *
          (options.tubeRadiusMode == TubeRadiusMode::PhysicalConstant
               ? 1.0 : options.outerRadiusM / options.tubeReferenceRadiusM)
      : 0.0;
  p.corridorPaddingM = active > refined ? active : refined;
  return Core::BuildDomainBoundsM(p, minimumM, maximumM);
}

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
  // Checkpoint compatibility deliberately excludes filesystem destinations
  // and the path used to locate the checkpoint being read.  Those values
  // must change when a resumed segment writes into a fresh evidence
  // directory, but they cannot change the evolved state.  The manifest still
  // contains every resolved physical input plus the output/checkpoint clocks,
  // so this is a narrower identity than the full provenance manifest rather
  // than a replacement for physics or storage-layout validation.
  const std::string& restart_compatibility_manifest() const {
    return restartCompatibilityManifest_;
  }

 private:
  RunConfiguration3D(const RunConfiguration3DOptions& options,
                     const StorageLayout& layout,
                     const std::string& physicsFingerprint,
                     const std::string& resolvedManifest,
                     const std::string& restartCompatibilityManifest);

  const RunConfiguration3DOptions options_;
  const StorageLayout storageLayout_;
  const std::string physicsFingerprint_;
  const std::string resolvedManifest_;
  const std::string restartCompatibilityManifest_;
};

}  // namespace RuntimeModel
}  // namespace SEP3D

#endif  // SEP3D_RUNTIME_RUN_CONFIGURATION_H
