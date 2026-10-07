#include "run_configuration.h"

#include "../core/parker_geometry.h"
#include "sep_background_snapshot.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <iomanip>
#include <limits>
#include <sstream>

namespace SEP3D {
namespace RuntimeModel {

namespace {

Core::Status Invalid(const std::string& message) {
  return Core::Status(Core::StatusCode::InvalidInput, message);
}

void AppendField(std::size_t count, std::size_t* cursor,
                 std::size_t* offset) {
  *offset = *cursor;
  *cursor += count * sizeof(double);
}

StorageLayout BuildLayout(const RunConfiguration3DOptions& options) {
  StorageLayout layout;
  std::size_t cursor = 0;
  AppendField(3, &cursor, &layout.magneticFieldOffset);
  AppendField(3, &cursor, &layout.bulkVelocityOffset);
  AppendField(1, &cursor, &layout.numberDensityOffset);
  AppendField(1, &cursor, &layout.velocityDivergenceOffset);
  // Phase B requires one complete ambient state in every populated cell.
  // Keeping these fields unconditional prevents the storage ABI from changing
  // when a run switches between Parker and SWMF background authorities.
  AppendField(1, &cursor, &layout.temperatureOffset);
  AppendField(1, &cursor, &layout.pressureOffset);
  AppendField(1, &cursor, &layout.alfvenSpeedOffset);
  AppendField(1, &cursor, &layout.divBhatOffset);
  AppendField(1, &cursor, &layout.focusingLengthOffset);
  AppendField(3, &cursor, &layout.curvatureOffset);
  AppendField(1, &cursor, &layout.fieldAlignedStrainOffset);
  if (options.storeMagneticGradient) {
    AppendField(9, &cursor, &layout.magneticGradientOffset);
  }
  if (options.storeVelocityGradient) {
    AppendField(9, &cursor, &layout.velocityGradientOffset);
  }
  // Every turbulence authority now publishes the same two directional
  // magnetic variances.  Prescribed waves are no less physical than imported
  // AWSoM waves, so reserving storage only for SWMF made the standalone
  // initialization product omit its wave state even though the provider had
  // evaluated it.  The stored values remain deltaB_+^2 and deltaB_-^2 [T^2];
  // Tecplot also emits w_+/-=deltaB_+/-^2/mu0 [J m^-3].
  AppendField(2, &cursor, &layout.waveEnergyOffset);
  layout.cellAssociatedBytes = cursor;
  layout.samplingBytesPerCell = options.samplingBytesPerCell;

  std::ostringstream canonical;
  canonical << "sep3d-storage-v2"
            << ";B=" << layout.magneticFieldOffset
            << ";U=" << layout.bulkVelocityOffset
            << ";n=" << layout.numberDensityOffset
            << ";divU=" << layout.velocityDivergenceOffset
            << ";T=" << layout.temperatureOffset
            << ";p=" << layout.pressureOffset
            << ";vA=" << layout.alfvenSpeedOffset
            << ";divb=" << layout.divBhatOffset
            << ";focus=" << layout.focusingLengthOffset
            << ";curvature=" << layout.curvatureOffset
            << ";strain=" << layout.fieldAlignedStrainOffset
            << ";gradB=" << layout.magneticGradientOffset
            << ";gradU=" << layout.velocityGradientOffset
            << ";waves=" << layout.waveEnergyOffset
            << ";cell_bytes=" << layout.cellAssociatedBytes
            << ";sample_bytes=" << layout.samplingBytesPerCell;
  layout.fingerprint = SEP::Background::FingerprintConfiguration(canonical.str());
  return layout;
}

double PresetOuterRadiusM(DomainPreset preset) {
  switch (preset) {
    case DomainPreset::Solar: return 0.30 * Core::Const::AU;
    case DomainPreset::OneAu:
    case DomainPreset::Earth: return Core::Const::AU;
    case DomainPreset::Mars: return 1.666 * Core::Const::AU;
  }
  return std::numeric_limits<double>::quiet_NaN();
}

bool FiniteVector(const Core::Vec3& value) {
  return std::isfinite(value.x) && std::isfinite(value.y) &&
         std::isfinite(value.z);
}

bool ValidProfileExponent(double value) {
  return std::isfinite(value) && value > 0.0;
}

bool NearlyEqual(double left, double right) {
  return std::fabs(left - right) <=
      128.0 * std::numeric_limits<double>::epsilon() *
      std::max(1.0, std::max(std::fabs(left), std::fabs(right)));
}

std::string NormalizedSpeciesSymbol(std::string symbol) {
  // Case and surrounding whitespace are presentation details, but punctuation
  // is chemical identity: H, H_PLUS, and H+ must never collapse to one value.
  // This normalization is used only for duplicate detection and diagnostics;
  // the original AMPS symbol is retained in the binding table.
  symbol.erase(std::remove_if(symbol.begin(), symbol.end(),
                              [](unsigned char value) {
                                return std::isspace(value) != 0;
                              }),
               symbol.end());
  std::transform(symbol.begin(), symbol.end(), symbol.begin(),
                 [](unsigned char value) {
                   return static_cast<char>(std::toupper(value));
                 });
  return symbol;
}

}  // namespace

const char* Name(BackgroundAuthority value) {
  switch (value) {
    case BackgroundAuthority::AnalyticParker: return "analytic-parker";
    case BackgroundAuthority::PythonInterpolator: return "python-interpolator";
    case BackgroundAuthority::Swmf: return "swmf";
    case BackgroundAuthority::Swcme: return "swcme";
    case BackgroundAuthority::RuntimeModel: return "runtime-model";
  }
  return "unknown";
}

const char* Name(PrescribedTurbulenceModel value) {
  switch (value) {
    case PrescribedTurbulenceModel::PowerLaw: return "power-law";
    case PrescribedTurbulenceModel::Kolmogorov: return "kolmogorov";
    case PrescribedTurbulenceModel::Kraichnan: return "kraichnan";
  }
  return "unknown";
}

const char* Name(PrescribedTurbulenceAmplitudeModel value) {
  switch (value) {
    case PrescribedTurbulenceAmplitudeModel::ConstantDeltaBOverB:
      return "constant-delta-b-over-b";
    case PrescribedTurbulenceAmplitudeModel::WaveEnergyPowerLaw:
      return "wave-energy-power-law";
  }
  return "unknown";
}

const char* Name(ParkerSpiralStartMode value) {
  switch (value) {
    case ParkerSpiralStartMode::Explicit: return "explicit";
    case ParkerSpiralStartMode::CmeLaunchPoint:
      return "cme-launch-point";
  }
  return "unknown";
}

const char* Name(SolarWindThermodynamicClosure value) {
  switch (value) {
    case SolarWindThermodynamicClosure::ProtonOnly: return "proton-only";
    case SolarWindThermodynamicClosure::MultiSpecies: return "multi-species";
  }
  return "unknown";
}

const char* Name(TurbulenceAuthority value) {
  switch (value) {
    case TurbulenceAuthority::Prescribed: return "prescribed";
    case TurbulenceAuthority::Swmf: return "swmf";
  }
  return "unknown";
}

const char* Name(ShockAuthority value) {
  switch (value) {
    case ShockAuthority::None: return "none";
    case ShockAuthority::Swcme: return "swcme";
  }
  return "unknown";
}

const char* Name(SourceSpectrumModel value) {
  switch (value) {
    case SourceSpectrumModel::LocalCompressionDsa:
      return "local-compression-dsa";
    case SourceSpectrumModel::FixedPhaseSpacePowerLaw:
      return "fixed-phase-space-power-law";
  }
  return "unknown";
}

const char* Name(SourceRateNormalizationModel value) {
  switch (value) {
    case SourceRateNormalizationModel::ConfiguredConstant:
      return "configured-constant";
    case SourceRateNormalizationModel::AcceptedShockIncidentFlux:
      return "accepted-shock-incident-flux";
  }
  return "unknown";
}

const char* Name(TransportModel value) {
  switch (value) {
    case TransportModel::Parker3D: return "parker";
    case TransportModel::FocusedDiffusion3D: return "focused-diffusion";
    case TransportModel::FocusedScattering3D: return "focused-scattering";
  }
  return "unknown";
}

const char* Name(DomainPreset value) {
  switch (value) {
    case DomainPreset::Solar: return "solar";
    case DomainPreset::OneAu: return "one-au";
    case DomainPreset::Earth: return "earth";
    case DomainPreset::Mars: return "mars";
  }
  return "unknown";
}

const char* Name(OuterRadiusMode value) {
  switch (value) {
    case OuterRadiusMode::Preset: return "preset";
    case OuterRadiusMode::Explicit: return "explicit";
    case OuterRadiusMode::FieldLineEndpoint: return "field-line-endpoint";
  }
  return "unknown";
}

const char* Name(DomainBoxGeometry value) {
  return value == DomainBoxGeometry::SunCenteredCube ? "sun-centered-cube" :
      (value == DomainBoxGeometry::FieldLineCornerCube ? "field-line-corner-cube" :
       (value == DomainBoxGeometry::FieldLineXYCornerCube ? "field-line-xy-corner-cube" : "unknown"));
}

const char* Name(SolarRefinementAnchor value) {
  return value == SolarRefinementAnchor::SourceShell ? "source-shell" :
      (value == SolarRefinementAnchor::Photosphere ? "photosphere" : "unknown");
}

const char* Name(InnerBoundaryMode value) {
  return value == InnerBoundaryMode::Absorb ? "absorb" : "unknown";
}

const char* Name(OuterBoundaryMode value) {
  switch (value) {
    case OuterBoundaryMode::Escape: return "escape";
    case OuterBoundaryMode::ImportedCoverage: return "imported-coverage";
  }
  return "unknown";
}

const char* Name(RefinementProfile value) {
  switch (value) {
    case RefinementProfile::Linear: return "linear";
    case RefinementProfile::PowerLaw: return "power-law";
    case RefinementProfile::Smoothstep: return "smoothstep";
  }
  return "unknown";
}

const char* Name(TubeRadiusMode value) {
  switch (value) {
    case TubeRadiusMode::PhysicalConstant: return "physical-constant";
    case TubeRadiusMode::ConstantAngularWidth: return "constant-angular-width";
  }
  return "unknown";
}

const char* Name(ActiveRegionMode value) {
  switch (value) {
    case ActiveRegionMode::FullDomain: return "full-domain";
    case ActiveRegionMode::ParkerTube: return "parker-tube";
  }
  return "unknown";
}

const char* Name(PopulationControlMode value) {
  switch (value) {
    case PopulationControlMode::Off: return "off";
    case PopulationControlMode::SplitMerge: return "split-merge";
  }
  return "unknown";
}

const char* Name(SpatialDiffusionModel value) {
  switch (value) {
    case SpatialDiffusionModel::CorrelationMeanFreePath:
      return "mean-free-path";
    case SpatialDiffusionModel::PitchAngleIntegral:
      return "pitch-angle-integral";
  }
  return "unknown";
}

const char* Name(PitchAngleDiffusionModel value) {
  switch (value) {
    case PitchAngleDiffusionModel::Jokipii1966: return "jokipii-1966";
    case PitchAngleDiffusionModel::Florinskiy: return "florinskiy";
    case PitchAngleDiffusionModel::Constant: return "constant";
  }
  return "unknown";
}

const char* Name(MeanFreePathModel value) {
  switch (value) {
    case MeanFreePathModel::Correlation: return "correlation";
    case MeanFreePathModel::Constant: return "constant";
    case MeanFreePathModel::RadialRigidityPowerLaw:
      return "radial-rigidity-power-law";
  }
  return "unknown";
}

const char* Name(FocusedScatteringFrame value) {
  switch (value) {
    case FocusedScatteringFrame::PlasmaFrameIsotropic:
      return "plasma-frame-isotropic";
    case FocusedScatteringFrame::AlfvenWaveFrameIsotropic:
      return "alfven-wave-frame-isotropic";
  }
  return "unknown";
}

const char* Name(RunIntent value) {
  switch (value) {
    case RunIntent::TransportOnly: return "transport-only";
    case RunIntent::ShockInjection: return "shock-injection";
    case RunIntent::ShockPropagation: return "shock-propagation";
  }
  return "unknown";
}

const char* Name(MissingTurbulenceMode value) {
  switch (value) {
    case MissingTurbulenceMode::Fail: return "fail";
    case MissingTurbulenceMode::Ballistic: return "ballistic";
  }
  return "unknown";
}

const char* Name(ResonanceRangeMode value) {
  switch (value) {
    case ResonanceRangeMode::Reject: return "reject";
    case ResonanceRangeMode::PowerLawExtension: return "power-law-extension";
  }
  return "unknown";
}

const char* Name(PitchAngleSchemeMode value) {
  switch (value) {
    case PitchAngleSchemeMode::ReflectingMilstein:
      return "reflecting-milstein";
    case PitchAngleSchemeMode::ReflectingEulerMaruyama:
      return "reflecting-euler-maruyama";
  }
  return "unknown";
}

const char* Name(PerpendicularDiffusionMode value) {
  switch (value) {
    case PerpendicularDiffusionMode::None: return "none";
    case PerpendicularDiffusionMode::Constant: return "constant";
    case PerpendicularDiffusionMode::ConstantRatio: return "constant-ratio";
  }
  return "unknown";
}

const char* Name(DriftMode value) {
  switch (value) {
    case DriftMode::None: return "none";
    case DriftMode::GradientB: return "gradient-b";
    case DriftMode::Curvature: return "curvature";
    case DriftMode::GradientAndCurvature: return "gradient-curvature";
  }
  return "unknown";
}

const char* Name(ObserverKind value) {
  switch (value) {
    case ObserverKind::FixedCartesian: return "fixed-cartesian";
    case ObserverKind::FixedHeliographic: return "fixed-heliographic";
    case ObserverKind::MovingCartesian: return "moving-cartesian";
    case ObserverKind::SphericalShell: return "spherical-shell";
    case ObserverKind::FieldConnected: return "field-connected";
  }
  return "unknown";
}

const char* Name(ObserverNormalization value) {
  switch (value) {
    case ObserverNormalization::RepresentedParticles:
      return "represented-particles";
    case ObserverNormalization::DifferentialIntensity:
      return "differential-intensity";
  }
  return "unknown";
}

const char* Name(EnergyChannelSpacing value) {
  switch (value) {
    case EnergyChannelSpacing::Logarithmic: return "logarithmic";
    case EnergyChannelSpacing::Linear: return "linear";
  }
  return "unknown";
}

bool operator==(const StorageLayout& left, const StorageLayout& right) {
  return left.magneticFieldOffset == right.magneticFieldOffset &&
         left.bulkVelocityOffset == right.bulkVelocityOffset &&
         left.numberDensityOffset == right.numberDensityOffset &&
         left.velocityDivergenceOffset == right.velocityDivergenceOffset &&
         left.temperatureOffset == right.temperatureOffset &&
         left.pressureOffset == right.pressureOffset &&
         left.alfvenSpeedOffset == right.alfvenSpeedOffset &&
         left.divBhatOffset == right.divBhatOffset &&
         left.focusingLengthOffset == right.focusingLengthOffset &&
         left.curvatureOffset == right.curvatureOffset &&
         left.fieldAlignedStrainOffset == right.fieldAlignedStrainOffset &&
         left.magneticGradientOffset == right.magneticGradientOffset &&
         left.velocityGradientOffset == right.velocityGradientOffset &&
         left.waveEnergyOffset == right.waveEnergyOffset &&
         left.cellAssociatedBytes == right.cellAssociatedBytes &&
         left.samplingBytesPerCell == right.samplingBytesPerCell &&
         left.fingerprint == right.fingerprint;
}

bool operator!=(const StorageLayout& left, const StorageLayout& right) {
  return !(left == right);
}

Core::Status ValidateCompiledSpeciesBinding(
    const RunConfiguration3DOptions& configured, int ampsSpeciesCount,
    const std::vector<CompiledSpeciesRecord>& compiled) {
  if (ampsSpeciesCount <= 0 ||
      compiled.size() != static_cast<std::size_t>(ampsSpeciesCount)) {
    return Core::Status(
        Core::StatusCode::ConfigurationConflict,
        "compiled AMPS species table size does not match the generated count");
  }

  std::vector<std::string> symbols;
  symbols.reserve(compiled.size());
  for (std::size_t slot = 0; slot < compiled.size(); ++slot) {
    const CompiledSpeciesRecord& species = compiled[slot];
    const std::string symbol = NormalizedSpeciesSymbol(species.symbol);
    if (species.ampsIndex != static_cast<int>(slot)) {
      return Core::Status(
          Core::StatusCode::ConfigurationConflict,
          "compiled AMPS species indices are not contiguous at slot " +
              std::to_string(slot));
    }
    if (symbol.empty())
      return Invalid("compiled AMPS species has an empty chemical symbol");
    if (std::find(symbols.begin(), symbols.end(), symbol) != symbols.end()) {
      return Core::Status(
          Core::StatusCode::ConfigurationConflict,
          "compiled AMPS species table contains duplicate symbol '" +
              species.symbol + "'");
    }
    symbols.push_back(symbol);
    if (!std::isfinite(species.massKg) || species.massKg <= 0.0) {
      return Invalid("compiled AMPS species '" + species.symbol +
                     "' has a non-positive or non-finite mass");
    }
    if (!std::isfinite(species.chargeC) || species.chargeC == 0.0) {
      return Core::Status(
          Core::StatusCode::ConfigurationConflict,
          "compiled AMPS species '" + species.symbol +
              "' is neutral or has a non-finite charge; the selected SEP "
              "scattering model requires every injected species to be charged");
    }
  }

  if (!configured.speciesParticleNormalizations.empty()) {
    if (configured.speciesParticleNormalizations.size() != compiled.size())
      return Core::Status(Core::StatusCode::ConfigurationConflict,
          "derived particle normalization does not cover every compiled species");
    for (std::size_t slot=0;slot<compiled.size();++slot) {
      const SpeciesParticleNormalization& normalization=
          configured.speciesParticleNormalizations[slot];
      if(normalization.ampsIndex!=compiled[slot].ampsIndex||
          NormalizedSpeciesSymbol(normalization.symbol)!=
              NormalizedSpeciesSymbol(compiled[slot].symbol))
        return Core::Status(Core::StatusCode::ConfigurationConflict,
            "derived particle normalization species identity differs from the generated AMPS table at slot "+
            std::to_string(slot));
    }
  }

  for (const ObserverOptions& observer : configured.observers) {
    // A wildcard has already been represented as an empty accepted-species
    // vector for the sampling layer.  It is valid for every positive compiled
    // table size and must not be expanded to a build-specific list in the
    // post-compile configuration object.
    if (observer.allCompiledSpecies) continue;
    for (int species : observer.species) {
      if (species < 0 || species >= ampsSpeciesCount) {
        return Core::Status(
            Core::StatusCode::ConfigurationConflict,
            "observer '" + observer.id + "' selects AMPS species index " +
                std::to_string(species) + " outside the compiled table [0," +
                std::to_string(ampsSpeciesCount - 1) + "]");
      }
    }
  }
  return Core::Status::OK();
}

RunConfiguration3D::RunConfiguration3D(
    const RunConfiguration3DOptions& options, const StorageLayout& layout,
    const std::string& physicsFingerprint,
    const std::string& resolvedManifest,
    const std::string& restartCompatibilityManifest)
    : options_(options),
      storageLayout_(layout),
      physicsFingerprint_(physicsFingerprint),
      resolvedManifest_(resolvedManifest),
      restartCompatibilityManifest_(restartCompatibilityManifest) {}

Core::Status RunConfiguration3D::Create(
    const RunConfiguration3DOptions& options,
    std::shared_ptr<const RunConfiguration3D>* configuration) {
  if (configuration == nullptr) return Invalid("configuration output is null");
  configuration->reset();

  // Normalize into a private copy.  No caller-owned options object is ever
  // modified, and no AMPS global is touched until this complete transaction
  // succeeds.  Preset resolution therefore becomes part of the immutable
  // configuration rather than a late mesh-builder side effect.
  RunConfiguration3DOptions normalized = options;
  if (normalized.outerRadiusMode != OuterRadiusMode::Preset &&
      normalized.outerRadiusMode != OuterRadiusMode::Explicit &&
      normalized.outerRadiusMode != OuterRadiusMode::FieldLineEndpoint)
    return Invalid("unknown outer radius mode");
  // Schema 3 obtains the momentum index from each canonical MHD shock state.
  // Zero is an explicit "provider-owned DSA" sentinel; it prevents the legacy
  // constant SourceOptions default from entering the physics fingerprint as if
  // it controlled injection.
  if (normalized.inputSchemaVersion >= 3)
    normalized.source.spectralIndex = 0.0;
  if (normalized.outerRadiusMode == OuterRadiusMode::Preset) {
    normalized.outerRadiusM = PresetOuterRadiusM(normalized.domain);
  }
  if (normalized.domain == DomainPreset::Earth) {
    normalized.domain = DomainPreset::OneAu;
  }
  if (normalized.outerRadiusMode == OuterRadiusMode::FieldLineEndpoint) {
    if (!std::isfinite(normalized.parkerSpiralEndRadiusM) ||
        normalized.parkerSpiralEndRadiusM <= normalized.innerRadiusM)
      return Invalid("field-line endpoint radius must exceed the source radius");
    Core::ParkerSpiralGeometry geometry;
    geometry.sourceRadiusM = normalized.innerRadiusM;
    geometry.sourceLongitudeRad = normalized.tubeLongitudeRad;
    geometry.sourceColatitudeRad = normalized.tubeColatitudeRad;
    geometry.solarWindSpeedMPerS = normalized.parker.solarWindSpeedMPerS;
    geometry.solarRotationRateRadPerS = normalized.parker.solarRotationRateRadPerS;
    geometry.rotationAxis = normalized.parker.rotationAxis;
    const double length = Core::ParkerCurveArcLengthM(
        normalized.parkerSpiralEndRadiusM, geometry);
    if (!std::isfinite(length) || length <= 0.0)
      return Invalid("field-line endpoint has invalid Parker arc length");
    // Zero selects derivation in a raw deck. A resolved options record may
    // carry the matching derived value through parser/factory round trips.
    if (normalized.parkerSpiralLengthM != 0.0 &&
        (!std::isfinite(normalized.parkerSpiralLengthM) ||
         std::fabs(normalized.parkerSpiralLengthM - length) > 1.0e-10 * length))
      return Invalid("length_m conflicts with the selected field-line endpoint; use zero for derivation");
    normalized.outerRadiusM = normalized.parkerSpiralEndRadiusM;
    normalized.parkerSpiralLengthM = length;
  } else if (normalized.parkerSpiralEndRadiusM != 0.0) {
    return Invalid("end_radius_m requires outer_radius_mode=field-line-endpoint");
  }
  normalized.parker.sourceRadiusM = normalized.innerRadiusM;
  normalized.parker.sourceLongitudeRad = normalized.tubeLongitudeRad;
  normalized.parker.sourceColatitudeRad = normalized.tubeColatitudeRad;

  // Version-1/programmatic construction did not expose a finite line.  Derive
  // an equivalent, deterministic definition before validation so old files
  // and typed construction still normalize to the same fingerprint.
  if (normalized.parkerSpiralPointCount == 0)
    normalized.parkerSpiralPointCount = 4001;
  if (normalized.parkerSpiralLengthM == 0.0)
    normalized.parkerSpiralLengthM = normalized.outerRadiusM - normalized.innerRadiusM;
  if (normalized.parkerSpiralInitialPointM.Norm() == 0.0) {
    const double sine = std::sin(normalized.tubeColatitudeRad);
    normalized.parkerSpiralInitialPointM = normalized.coordinateOriginM +
        normalized.innerRadiusM * Core::Vec3(
            sine * std::cos(normalized.tubeLongitudeRad),
            sine * std::sin(normalized.tubeLongitudeRad),
            std::cos(normalized.tubeColatitudeRad));
  }

  if (!std::isfinite(normalized.innerRadiusM) ||
      normalized.innerRadiusM < Core::Const::R_sun) {
    return Invalid(
        "innerRadiusM must be finite and at or above the physical solar "
        "surface (Core::Const::R_sun)");
  }
  if (!std::isfinite(normalized.outerRadiusM) ||
      normalized.outerRadiusM <= normalized.innerRadiusM) {
    return Invalid("outerRadiusM must be finite and exceed innerRadiusM");
  }
  if (!FiniteVector(normalized.coordinateOriginM) ||
      normalized.coordinateOriginM.Norm() != 0.0 ||
      normalized.coordinateFrame.empty()) {
    return Invalid("the current Parker/SWMF contract requires a finite heliocentric origin and named frame");
  }
  if (normalized.inputSchemaVersion < 1 || normalized.inputSchemaVersion > 4)
    return Invalid("inputSchemaVersion must be 1, 2, 3, or 4");
  if (normalized.background == BackgroundAuthority::PythonInterpolator) {
    // The provider value is intentionally recognized before any AMPS state is
    // touched, but execution remains fail-closed until the Python process,
    // units, batching, error, and provenance protocol has its own validation
    // gate.  It is never treated as Parker or SWMF by an else branch.
    return Core::Status::Reserved(
        "Python heliospheric-model interpolation background");
  }
  if (!FiniteVector(normalized.parkerSpiralOriginM) ||
      !FiniteVector(normalized.parkerSpiralInitialPointM) ||
      (normalized.parkerSpiralOriginM - normalized.coordinateOriginM).Norm() != 0.0 ||
      !std::isfinite(normalized.parkerSpiralLengthM) ||
      normalized.parkerSpiralLengthM <= 0.0 ||
      normalized.parkerSpiralPointCount < 2 ||
      normalized.parkerSpiralPointCount > 10000000ULL) {
    return Invalid("finite Parker spiral origin, length, or point count is invalid");
  }
  if (normalized.parkerSpiralStartMode != ParkerSpiralStartMode::Explicit &&
      normalized.parkerSpiralStartMode !=
          ParkerSpiralStartMode::CmeLaunchPoint) {
    return Invalid("Parker spiral start mode is unknown");
  }
  if (normalized.parkerSpiralStartMode ==
      ParkerSpiralStartMode::CmeLaunchPoint) {
    if (!normalized.cmeLaunchPointResolved ||
        !FiniteVector(normalized.cmeLaunchPointM)) {
      return Invalid("CME-linked Parker start has no canonically resolved "
                     "SWCME launch-apex point");
    }
    const double linkageScale = std::max(
        1.0, std::max(normalized.innerRadiusM,
                      normalized.cmeLaunchPointM.Norm()));
    if ((normalized.parkerSpiralInitialPointM -
         normalized.cmeLaunchPointM).Norm() > 1.0e-12 * linkageScale) {
      return Invalid("CME-linked Parker start differs from the canonical "
                     "SWCME launch-apex point");
    }
  } else if (normalized.cmeLaunchPointResolved) {
    return Invalid("an explicit Parker start must not carry a resolved CME "
                   "launch-apex linkage");
  }
  const Core::Vec3 lineSource =
      normalized.parkerSpiralInitialPointM - normalized.parkerSpiralOriginM;
  const double sourceRadius = lineSource.Norm();
  const double geometryTolerance = 1.0e-10 * normalized.innerRadiusM;
  Core::ParkerSpiralGeometry declaredGeometry;
  declaredGeometry.sourceRadiusM = normalized.innerRadiusM;
  declaredGeometry.sourceLongitudeRad = normalized.tubeLongitudeRad;
  declaredGeometry.sourceColatitudeRad = normalized.tubeColatitudeRad;
  declaredGeometry.solarWindSpeedMPerS = normalized.parker.solarWindSpeedMPerS;
  declaredGeometry.solarRotationRateRadPerS = normalized.parker.solarRotationRateRadPerS;
  declaredGeometry.rotationAxis = normalized.parker.rotationAxis;
  const Core::Vec3 declaredSource =
      Core::ParkerCurvePoint(normalized.innerRadiusM, declaredGeometry);
  if (std::fabs(sourceRadius - normalized.innerRadiusM) > geometryTolerance ||
      (lineSource - declaredSource).Norm() > geometryTolerance) {
    return Invalid("Parker initial point must lie at inner_radius_m and match the mesh-tube source angles");
  }
  if (normalized.outerBoundary == OuterBoundaryMode::ImportedCoverage &&
      normalized.background != BackgroundAuthority::Swmf) {
    return Core::Status(Core::StatusCode::ConfigurationConflict,
                        "imported-coverage outer boundary requires SWMF background authority");
  }
  if (!std::isfinite(normalized.requestedTimeStepS) ||
      normalized.requestedTimeStepS <= 0.0) {
    return Invalid("requestedTimeStepS must be finite and positive");
  }
  if (normalized.maximumTimeSteps == 0)
    return Invalid("maximumTimeSteps must be positive");
  if (normalized.campaignSeed == 0) return Invalid("campaignSeed zero is reserved");
  if (normalized.backgroundCadenceSteps == 0) {
    return Invalid("backgroundCadenceSteps must be positive");
  }
  if (normalized.injectionCadenceSteps == 0) {
    return Invalid("injectionCadenceSteps must be positive");
  }
  if (normalized.outputCadenceSteps == 0) {
    return Invalid("outputCadenceSteps must be positive");
  }
  if (normalized.outputDirectory.empty() || normalized.outputPrefix.empty()) {
    return Invalid("output directory and prefix must not be empty");
  }
  const bool reducedShockFront=normalized.background==BackgroundAuthority::RuntimeModel&&
      normalized.backgroundModelId=="sep-corona-swcme-shock-front-v1";
  if (normalized.inputSchemaVersion >= 3 && !reducedShockFront&&
      (normalized.initializationMeshTecplotFile.empty() ||
       normalized.initializationParkerLineTecplotFile.empty() ||
       normalized.initializationDataTecplotFile.empty() ||
       normalized.swcmeAssignments.empty() ||
       normalized.swcmeConfigurationFingerprint.empty() ||
       normalized.swcmeResolvedManifest.empty())) {
    return Invalid("schema version 3 requires validated SWCME and initialization Tecplot outputs");
  }
  const bool propagation = normalized.intent == RunIntent::ShockPropagation;
  if (normalized.inputSchemaVersion >= 3 && !reducedShockFront&&
      (normalized.shock != ShockAuthority::Swcme ||
       (!propagation && (normalized.intent != RunIntent::ShockInjection ||
                        !normalized.source.enabled)))) {
    return Invalid("canonical SWCME requires shock-injection or shock-propagation intent");
  }
  if (!std::isfinite(normalized.stopShockRadiusM) || normalized.stopShockRadiusM < 0.0 ||
      (normalized.stopShockRadiusM > 0.0 &&
       (!propagation || normalized.stopShockRadiusM <= normalized.shockModel.initialRadiusM ||
        normalized.stopShockRadiusM > normalized.outerRadiusM))) {
    return Invalid("run.stop_shock_radius_m requires propagation and a target above launch within the outer boundary");
  }
  // A propagation run needs a surface authority, but the authority need not
  // be the legacy particle-facing SWCME adapter.  The reduced composite owns
  // its front inside the runtime background provider and deliberately keeps
  // [shock].authority=none so no particle crossing/source path can mistake
  // its immediate RH limits for a downstream volume.  Admit exactly that
  // frozen pairing; all other runtime-model/shock combinations remain
  // rejected.  Restart remains excluded from this driver intent because its
  // history writer requires a contiguous tick-zero file.
  const bool propagationAuthority =
      (!reducedShockFront && normalized.shock == ShockAuthority::Swcme) ||
      (reducedShockFront && normalized.shock == ShockAuthority::None);
  if (propagation && (normalized.inputSchemaVersion < 4 ||
      !propagationAuthority || normalized.source.enabled ||
      normalized.populationControl != PopulationControlMode::Off ||
      !normalized.restartInputPath.empty())) {
    return Invalid("shock-propagation requires schema 4, a supported front authority, source.enabled=false, population control off and a fresh run");
  }
  if (normalized.inputSchemaVersion >= 3 &&
      normalized.injectionCadenceSteps != 1) {
    return Invalid("schema version 3 requires source injection on every time step");
  }
  if (normalized.inputSchemaVersion >= 3 &&
      ((normalized.background != BackgroundAuthority::AnalyticParker &&
        !(normalized.inputSchemaVersion >= 4 &&
          (normalized.background == BackgroundAuthority::Swcme ||
           normalized.background == BackgroundAuthority::RuntimeModel))) ||
       normalized.turbulence != TurbulenceAuthority::Prescribed)) {
    return Invalid("standalone schema version 3 requires analytic Parker "
                   "background and prescribed turbulence");
  }
  // The SWCME background and shock share one canonical event configuration.
  // Reject a missing canonical shock owner before factories or AMPS allocation.
  if (normalized.background == BackgroundAuthority::Swcme &&
      (normalized.inputSchemaVersion < 4 || normalized.shock != ShockAuthority::Swcme))
    return Invalid("SWCME background requires schema 4 and canonical SWCME runtime");
  // Require exactly the authority/ID pairing, including parser-free hosts.
  // An unused model_id would otherwise make input provenance ambiguous.
  if ((normalized.background == BackgroundAuthority::RuntimeModel) !=
      !normalized.backgroundModelId.empty())
    return Invalid("background.model_id is required only for runtime-model authority");
  const bool fileBackedReducedModel=!normalized.backgroundModelAssetPath.empty();
  const bool inlineReducedModel=!normalized.backgroundModelInlineConfiguration.empty();
  if(reducedShockFront&&(fileBackedReducedModel==inlineReducedModel))
    return Invalid("reduced shock-front runtime requires exactly one file-backed or inline model configuration");
  if(inlineReducedModel&&normalized.backgroundModelAssetDirectory.empty())
    return Invalid("inline reduced shock-front configuration requires its asset directory");
  if(!inlineReducedModel&&!normalized.backgroundModelAssetDirectory.empty())
    return Invalid("background model asset directory is active only for inline configuration");
  if(!reducedShockFront&&(fileBackedReducedModel||inlineReducedModel))
    return Invalid("reduced shock-front model input is active only for its runtime authority");
  const double meshValues[] = {
      normalized.minimumCellSizeM, normalized.backgroundCellSizeM,
      normalized.solarSurfaceCellSizeM,
      normalized.solarRefinementOuterRadiusM,
      normalized.tubeLongitudeRad, normalized.tubeColatitudeRad,
      normalized.tubeReferenceRadiusM,
      normalized.tubeRadiusAtReferenceM, normalized.tubeCellSizeM,
      normalized.activeTubeReferenceRadiusM,
      normalized.activeTubeRadiusAtReferenceM};
  for (double value : meshValues) {
    if (!std::isfinite(value)) return Invalid("mesh option is not finite");
  }
  if (normalized.minimumCellSizeM <= 0.0 ||
      normalized.backgroundCellSizeM < normalized.minimumCellSizeM ||
      normalized.solarSurfaceCellSizeM < normalized.minimumCellSizeM ||
      normalized.solarSurfaceCellSizeM > normalized.backgroundCellSizeM ||
      normalized.solarRefinementOuterRadiusM <=
          (normalized.solarRefinementAnchor == SolarRefinementAnchor::Photosphere
              ? Core::Const::R_sun : normalized.innerRadiusM) ||
      normalized.solarRefinementOuterRadiusM > normalized.outerRadiusM ||
      !ValidProfileExponent(normalized.solarRefinementExponent) ||
      !ValidProfileExponent(normalized.tubeTransverseExponent) ||
      normalized.meshCellsPerBlockEdge == 0 ||
      normalized.maximumMeshLevel > 19 ||
      normalized.meshMemoryBudgetBytes == 0) {
    return Invalid("mesh sizes, level, block width, or memory budget are invalid");
  }
  if (normalized.solarRefinementAnchor != SolarRefinementAnchor::SourceShell &&
      normalized.solarRefinementAnchor != SolarRefinementAnchor::Photosphere)
    return Invalid("unknown solar refinement anchor");
  if (!std::isfinite(normalized.activeSolarSphereRadiusM) ||
      normalized.activeSolarSphereRadiusM < 0.0 ||
      (normalized.activeSolarSphereRadiusM > 0.0 &&
       (normalized.activeRegion != ActiveRegionMode::ParkerTube ||
        normalized.activeSolarSphereRadiusM < normalized.innerRadiusM ||
        normalized.activeSolarSphereRadiusM >= normalized.outerRadiusM)))
    return Invalid("active solar sphere must connect to the corridor source and lie inside the endpoint radius");
  if (!std::isfinite(normalized.domainCornerMarginM) || normalized.domainCornerMarginM < 0.0 ||
      !FiniteVector(normalized.domainCornerDirection))
    return Invalid("corner margin/direction must be finite and margin nonnegative");
  if ((normalized.domainBoxGeometry == DomainBoxGeometry::FieldLineCornerCube ||
       normalized.domainBoxGeometry == DomainBoxGeometry::FieldLineXYCornerCube) &&
      (normalized.outerRadiusMode != OuterRadiusMode::FieldLineEndpoint ||
       normalized.activeRegion != ActiveRegionMode::ParkerTube ||
       normalized.activeSolarSphereRadiusM == 0.0))
    return Invalid("corner cube requires an endpoint-derived active corridor plus solar sphere");
  if (normalized.domainBoxGeometry == DomainBoxGeometry::SunCenteredCube &&
      (normalized.domainCornerMarginM != 0.0 || normalized.domainCornerDirection.Norm() != 0.0))
    return Invalid("corner controls require a field-line corner geometry");
  if (normalized.domainBoxGeometry == DomainBoxGeometry::FieldLineXYCornerCube &&
      normalized.domainCornerDirection.z != 0.0)
    return Invalid("x-y corner geometry requires corner_direction_z=0");
  if (normalized.enableTubeRefinement &&
      (normalized.tubeReferenceRadiusM <= normalized.innerRadiusM ||
       normalized.tubeRadiusAtReferenceM <= 0.0 ||
       normalized.tubeCellSizeM < normalized.minimumCellSizeM ||
       normalized.tubeCellSizeM > normalized.backgroundCellSizeM ||
       normalized.tubeColatitudeRad < 0.0 ||
       normalized.tubeColatitudeRad > Core::Const::kPi)) {
    return Invalid("Parker-tube refinement options are invalid");
  }
  if (normalized.activeRegion != ActiveRegionMode::FullDomain &&
      normalized.activeRegion != ActiveRegionMode::ParkerTube)
    return Invalid("active-region mode is unknown");
  if (normalized.activeTubeRadiusMode != TubeRadiusMode::PhysicalConstant &&
      normalized.activeTubeRadiusMode !=
          TubeRadiusMode::ConstantAngularWidth)
    return Invalid("active-region radius mode is unknown");
  if (normalized.activeRegion == ActiveRegionMode::ParkerTube) {
    if (normalized.activeTubeReferenceRadiusM <= normalized.innerRadiusM ||
        normalized.activeTubeRadiusAtReferenceM <= 0.0 ||
        normalized.activeTubeBufferBlocks == 0) {
      return Invalid("Parker active-region radius, reference radius, and "
                     "buffer_blocks must be positive");
    }
    // A corridor narrower than the requested refinement tube would discard
    // blocks that the mesh explicitly refined for transport. Both supported
    // radius laws are affine in heliocentric radius (constant or proportional
    // to r), so checking both finite-line endpoints proves containment over
    // the complete interval. The former single-reference test was insufficient
    // when the refinement and active tubes selected different radius modes.
    if (normalized.enableTubeRefinement) {
      const double outerLineLengthM = Core::ParkerCurveArcLengthM(
          normalized.outerRadiusM, declaredGeometry);
      const double finiteLineLengthM = std::min(
          normalized.parkerSpiralLengthM, outerLineLengthM);
      double terminalRadiusM = 0.0;
      const Core::Status terminalStatus =
          Core::ParkerCurveRadiusAtArcLengthM(
              finiteLineLengthM, declaredGeometry, &terminalRadiusM);
      if (!terminalStatus.ok()) return terminalStatus;
      const double endpointRadii[] = {
          normalized.innerRadiusM, terminalRadiusM};
      for (double radiusM : endpointRadii) {
        const double activeRadiusM =
            normalized.activeTubeRadiusMode ==
                    TubeRadiusMode::PhysicalConstant
                ? normalized.activeTubeRadiusAtReferenceM
                : normalized.activeTubeRadiusAtReferenceM * radiusM /
                      normalized.activeTubeReferenceRadiusM;
        const double refinementRadiusM =
            normalized.tubeRadiusMode == TubeRadiusMode::PhysicalConstant
                ? normalized.tubeRadiusAtReferenceM
                : normalized.tubeRadiusAtReferenceM * radiusM /
                      normalized.tubeReferenceRadiusM;
        const double toleranceM = 64.0 *
            std::numeric_limits<double>::epsilon() *
            std::max({1.0, activeRadiusM, refinementRadiusM});
        if (activeRadiusM + toleranceM < refinementRadiusM) {
          return Invalid("active Parker corridor is narrower than the "
                         "refined Parker tube on the finite line");
        }
      }
    }
  } else if (normalized.inputSchemaVersion >= 4 &&
             (normalized.activeTubeRadiusAtReferenceM != 0.0 ||
              normalized.activeTubeBufferBlocks != 0)) {
    return Invalid("full-domain active region requires zero inactive tube "
                   "radius and buffer_blocks");
  }
  Core::Vec3 domainMinimumM, domainMaximumM;
  const Core::Status boundsStatus = ResolveDomainBoundsM(
      normalized, &domainMinimumM, &domainMaximumM);
  if (!boundsStatus.ok()) return boundsStatus;
  const double rootCellM = (domainMaximumM.x - domainMinimumM.x) /
      static_cast<double>(normalized.meshCellsPerBlockEdge);
  const double finestAvailableM = rootCellM /
      std::pow(2.0, static_cast<double>(normalized.maximumMeshLevel));
  double finestRequestedM = normalized.backgroundCellSizeM;
  if (normalized.enableRadialRefinement)
    finestRequestedM = std::min(finestRequestedM,
                                normalized.solarSurfaceCellSizeM);
  if (normalized.enableTubeRefinement)
    finestRequestedM = std::min(finestRequestedM,
                                normalized.tubeCellSizeM);
  if (finestAvailableM > finestRequestedM * (1.0 + 1.0e-12)) {
    std::ostringstream message;
    message << "maximumMeshLevel cannot realize requested cell size: available="
            << finestAvailableM << " requested=" << finestRequestedM;
    return Core::Status(Core::StatusCode::LayoutMismatch, message.str());
  }
  if (normalized.memoryModel.baseCellBytes == 0 ||
      normalized.memoryModel.baseNodeBytes == 0 ||
      normalized.memoryModel.blockStructureBytes == 0 ||
      normalized.memoryModel.particleBytes == 0 ||
      !std::isfinite(normalized.memoryModel.particlesPerCell) ||
      normalized.memoryModel.particlesPerCell < 0.0 ||
      !std::isfinite(normalized.memoryModel.haloFraction) ||
      normalized.memoryModel.haloFraction < 0.0 ||
      !std::isfinite(normalized.memoryModel.safetyMarginFraction) ||
      normalized.memoryModel.safetyMarginFraction < 0.0) {
    return Invalid("memory-model coefficients or safety factors are invalid");
  }
  if (normalized.prescribedTurbulenceModel !=
          PrescribedTurbulenceModel::PowerLaw &&
      normalized.prescribedTurbulenceModel !=
          PrescribedTurbulenceModel::Kolmogorov &&
      normalized.prescribedTurbulenceModel !=
          PrescribedTurbulenceModel::Kraichnan) {
    return Invalid("prescribed turbulence model is unknown");
  }
  if (normalized.prescribedTurbulenceAmplitudeModel !=
          PrescribedTurbulenceAmplitudeModel::ConstantDeltaBOverB &&
      normalized.prescribedTurbulenceAmplitudeModel !=
          PrescribedTurbulenceAmplitudeModel::WaveEnergyPowerLaw) {
    return Invalid("prescribed turbulence amplitude model is unknown");
  }
  if (!std::isfinite(normalized.prescribedDeltaBOverB) ||
      !std::isfinite(normalized.turbulenceWaveEnergyAtReferenceJPerM3) ||
      !std::isfinite(normalized.turbulenceWaveEnergyRadialExponent) ||
      !std::isfinite(normalized.turbulenceNormalizedCrossHelicity) ||
      normalized.turbulenceNormalizedCrossHelicity < -1.0 ||
      normalized.turbulenceNormalizedCrossHelicity > 1.0 ||
      !std::isfinite(normalized.turbulenceReferenceRadiusM) ||
      normalized.turbulenceReferenceRadiusM <= 0.0 ||
      !std::isfinite(normalized.turbulenceKMinPerM) ||
      normalized.turbulenceKMinPerM <= 0.0 ||
      !std::isfinite(normalized.turbulenceKMaxPerM) ||
      normalized.turbulenceKMaxPerM <= normalized.turbulenceKMinPerM ||
      !std::isfinite(normalized.turbulenceKMinRadialExponent) ||
      !std::isfinite(normalized.turbulenceKMaxRadialExponent) ||
      !std::isfinite(normalized.turbulenceSpectralIndex) ||
      normalized.turbulenceSpectralIndex <= 1.0 ||
      !std::isfinite(normalized.turbulenceCorrelationLengthM) ||
      normalized.turbulenceCorrelationLengthM <= 0.0 ||
      !std::isfinite(normalized.turbulenceCorrelationLengthRadialExponent) ||
      !std::isfinite(normalized.turbulenceValidityCadenceS) ||
      normalized.turbulenceValidityCadenceS <= 0.0) {
    return Invalid("turbulence amplitude, imbalance, radial scaling, band, "
                   "index, correlation length, or cadence is invalid");
  }
  // Each amplitude closure has exactly one active normalization.  Requiring a
  // zero sentinel for the inactive branch prevents a complete input deck from
  // carrying a plausible-looking number that does not affect the run.
  if (normalized.prescribedTurbulenceAmplitudeModel ==
          PrescribedTurbulenceAmplitudeModel::ConstantDeltaBOverB) {
    if (normalized.prescribedDeltaBOverB <= 0.0 ||
        normalized.turbulenceWaveEnergyAtReferenceJPerM3 != 0.0 ||
        normalized.turbulenceWaveEnergyRadialExponent != 0.0) {
      return Invalid("constant-delta-b-over-b turbulence requires positive "
                     "delta_b_over_b and zero inactive wave-energy inputs");
    }
  } else {
    if (normalized.prescribedDeltaBOverB != 0.0 ||
        normalized.turbulenceWaveEnergyAtReferenceJPerM3 <= 0.0) {
      return Invalid("wave-energy-power-law turbulence requires zero inactive "
                     "delta_b_over_b and positive reference wave energy");
    }
  }
  // Kolmogorov and Kraichnan are named physical closures, not aliases that
  // silently overwrite a contradictory number.  The explicit index remains
  // in complete input and must agree with the selected closure.  power-law is
  // the opt-in route for another validated q>1 slope.
  if ((normalized.prescribedTurbulenceModel ==
           PrescribedTurbulenceModel::Kolmogorov &&
       !NearlyEqual(normalized.turbulenceSpectralIndex, 5.0 / 3.0)) ||
      (normalized.prescribedTurbulenceModel ==
           PrescribedTurbulenceModel::Kraichnan &&
       !NearlyEqual(normalized.turbulenceSpectralIndex, 3.0 / 2.0))) {
    return Invalid("turbulence.spectral_index contradicts the selected named "
                   "turbulence model");
  }
  if (normalized.turbulence == TurbulenceAuthority::Swmf &&
      normalized.background != BackgroundAuthority::Swmf) {
    return Core::Status(
        Core::StatusCode::ConfigurationConflict,
        "SWMF turbulence requires the SWMF background authority");
  }
  const double transportControls[] = {
      normalized.cellCrossingFraction, normalized.diffusionFraction,
      normalized.focusingFraction, normalized.coolingFraction,
      normalized.fieldVariationFraction, normalized.shockCrossingFraction,
      normalized.minimumTransportSubstepS};
  for (double value : transportControls) {
    if (!std::isfinite(value) || value <= 0.0) {
      return Invalid("transport time-step controls must be finite and positive");
    }
  }
  if (normalized.maximumTransportSubsteps == 0) {
    return Invalid("maximumTransportSubsteps must be positive");
  }
  if (normalized.transport != TransportModel::Parker3D &&
      normalized.transport != TransportModel::FocusedDiffusion3D &&
      normalized.transport != TransportModel::FocusedScattering3D)
    return Invalid("transport mover is unknown");
  if (normalized.spatialDiffusionModel !=
          SpatialDiffusionModel::CorrelationMeanFreePath &&
      normalized.spatialDiffusionModel !=
          SpatialDiffusionModel::PitchAngleIntegral)
    return Invalid("spatial-diffusion model is unknown");
  if (normalized.pitchAngleDiffusionModel !=
          PitchAngleDiffusionModel::Jokipii1966 &&
      normalized.pitchAngleDiffusionModel !=
          PitchAngleDiffusionModel::Florinskiy &&
      normalized.pitchAngleDiffusionModel !=
          PitchAngleDiffusionModel::Constant)
    return Invalid("pitch-angle-diffusion model is unknown");
  if (normalized.meanFreePathModel != MeanFreePathModel::Correlation &&
      normalized.meanFreePathModel != MeanFreePathModel::Constant &&
      normalized.meanFreePathModel !=
          MeanFreePathModel::RadialRigidityPowerLaw)
    return Invalid("mean-free-path model is unknown");
  if (normalized.focusedScatteringFrame !=
          FocusedScatteringFrame::PlasmaFrameIsotropic &&
      normalized.focusedScatteringFrame !=
          FocusedScatteringFrame::AlfvenWaveFrameIsotropic)
    return Invalid("focused-scattering frame is unknown");
  if (!std::isfinite(normalized.constantDmumuPerS) ||
      normalized.constantDmumuPerS < 0.0 ||
      !std::isfinite(normalized.constantMeanFreePathM) ||
      normalized.constantMeanFreePathM < 0.0 ||
      !std::isfinite(normalized.meanFreePathReferenceM) ||
      normalized.meanFreePathReferenceM < 0.0 ||
      !std::isfinite(normalized.meanFreePathReferenceRadiusM) ||
      normalized.meanFreePathReferenceRadiusM < 0.0 ||
      !std::isfinite(normalized.meanFreePathReferenceRigidityV) ||
      normalized.meanFreePathReferenceRigidityV < 0.0 ||
      !std::isfinite(normalized.meanFreePathRadialExponent) ||
      !std::isfinite(normalized.meanFreePathRigidityExponent) ||
      !std::isfinite(normalized.spatialQuadratureAbsoluteToleranceM2PerS) ||
      normalized.spatialQuadratureAbsoluteToleranceM2PerS < 0.0 ||
      !std::isfinite(normalized.spatialQuadratureRelativeTolerance) ||
      normalized.spatialQuadratureRelativeTolerance <= 0.0 ||
      normalized.spatialQuadratureMaximumRecursion == 0 ||
      normalized.spatialQuadratureMaximumRecursion > 60 ||
      normalized.maximumScatteringEventsPerSubstep == 0) {
    return Invalid("parallel transport coefficient values or quadrature "
                   "controls are invalid");
  }
  if (normalized.pitchAngleDiffusionModel ==
          PitchAngleDiffusionModel::Constant &&
      normalized.constantDmumuPerS <= 0.0)
    return Invalid("constant pitch-angle diffusion requires "
                   "constant_dmumu_per_s > 0");
  if (normalized.pitchAngleDiffusionModel !=
          PitchAngleDiffusionModel::Constant &&
      normalized.inputSchemaVersion >= 4 && normalized.constantDmumuPerS != 0.0)
    return Invalid("a non-constant pitch-angle model requires zero inactive "
                   "constant_dmumu_per_s");
  if (normalized.meanFreePathModel == MeanFreePathModel::Constant &&
      normalized.constantMeanFreePathM <= 0.0)
    return Invalid("constant mean-free-path model requires "
                   "constant_mean_free_path_m > 0");
  if (normalized.meanFreePathModel != MeanFreePathModel::Constant &&
      normalized.inputSchemaVersion >= 4 &&
      normalized.constantMeanFreePathM != 0.0)
    return Invalid("a non-constant mean-free-path model requires zero inactive "
                   "constant_mean_free_path_m");
  const bool radialRigidity = normalized.meanFreePathModel ==
      MeanFreePathModel::RadialRigidityPowerLaw;
  if (radialRigidity &&
      (normalized.meanFreePathReferenceM <= 0.0 ||
       normalized.meanFreePathReferenceRadiusM <= 0.0 ||
       normalized.meanFreePathReferenceRigidityV <= 0.0)) {
    return Invalid("radial-rigidity-power-law mean free path requires positive "
                   "reference length, radius, and rigidity");
  }
  if (!radialRigidity && normalized.inputSchemaVersion >= 4 &&
      (normalized.meanFreePathReferenceM != 0.0 ||
       normalized.meanFreePathReferenceRadiusM != 0.0 ||
       normalized.meanFreePathReferenceRigidityV != 0.0 ||
       normalized.meanFreePathRadialExponent != 0.0 ||
       normalized.meanFreePathRigidityExponent != 0.0)) {
    return Invalid("a non-power-law mean-free-path model requires zero inactive "
                   "radial/rigidity power-law parameters");
  }
  if (normalized.transport == TransportModel::FocusedScattering3D &&
      normalized.perpendicularDiffusion != PerpendicularDiffusionMode::None) {
    return Invalid("focused-scattering currently requires perpendicular_"
                   "diffusion=none; an event-partition-invariant transverse "
                   "operator has not been validated");
  }

  if (normalized.populationControl != PopulationControlMode::Off &&
      normalized.populationControl != PopulationControlMode::SplitMerge)
    return Invalid("particle population-control mode is unknown");
  if (normalized.populationControl == PopulationControlMode::SplitMerge) {
    if (normalized.minimumParticlesPerCellPerSpecies < 2 ||
        normalized.targetParticlesPerCellPerSpecies <
            normalized.minimumParticlesPerCellPerSpecies ||
        normalized.maximumParticlesPerCellPerSpecies <
            normalized.targetParticlesPerCellPerSpecies ||
        normalized.populationControlCadenceSteps == 0) {
      return Invalid("split-merge population control requires 2 <= minimum "
                     "<= target <= maximum and a positive cadence");
    }
  } else if (normalized.inputSchemaVersion >= 4 &&
             (normalized.minimumParticlesPerCellPerSpecies != 0 ||
              normalized.targetParticlesPerCellPerSpecies != 0 ||
              normalized.maximumParticlesPerCellPerSpecies != 0)) {
    return Invalid("disabled population control requires zero inactive limits");
  }

  const ParkerPhysicsOptions& parker = normalized.parker;
  const double parkerValues[] = {
      parker.sourceRadiusM, parker.sourceLongitudeRad,
      parker.sourceColatitudeRad, parker.referenceRadiusM,
      parker.radialFieldAtReferenceT, parker.solarRotationRateRadPerS,
      parker.solarWindSpeedMPerS, parker.numberDensityAtReferenceM3,
      parker.densityReferenceRadiusM,
      parker.temperatureK, parker.adiabaticIndex,
      parker.alphaToProtonRatio, parker.electronTemperatureK,
      parker.alphaTemperatureK, parker.referenceSinColatitude,
      parker.validityCadenceS};
  for (double value : parkerValues)
    if (!std::isfinite(value)) return Invalid("Parker configuration contains a non-finite value");
  if (parker.thermodynamicClosure !=
          SolarWindThermodynamicClosure::ProtonOnly &&
      parker.thermodynamicClosure !=
          SolarWindThermodynamicClosure::MultiSpecies) {
    return Invalid("Parker thermodynamic closure is unknown");
  }
  if (parker.sourceRadiusM != normalized.innerRadiusM ||
      parker.referenceRadiusM <= parker.sourceRadiusM ||
      parker.radialFieldAtReferenceT <= 0.0 ||
      parker.solarRotationRateRadPerS < 0.0 ||
      parker.solarWindSpeedMPerS <= 0.0 ||
      parker.numberDensityAtReferenceM3 <= 0.0 ||
      parker.densityReferenceRadiusM <= parker.sourceRadiusM ||
      parker.temperatureK <= 0.0 ||
      parker.adiabaticIndex <= 1.0 || parker.alphaToProtonRatio < 0.0 ||
      parker.electronTemperatureK <= 0.0 ||
      parker.alphaTemperatureK <= 0.0 ||
      parker.referenceSinColatitude < 0.0 ||
      parker.referenceSinColatitude > 1.0 ||
      !FiniteVector(parker.rotationAxis) || parker.rotationAxis.Norm() <= 0.0 ||
      parker.validityCadenceS <= 0.0 ||
      (parker.magneticPolarity != 1 && parker.magneticPolarity != -1) ||
      parker.sourceColatitudeRad < 0.0 ||
      parker.sourceColatitudeRad > Core::Const::kPi ||
      parker.coordinateFrame != normalized.coordinateFrame) {
    return Invalid("Parker physical configuration is inconsistent or outside its range");
  }

  if (normalized.intent == RunIntent::ShockInjection) {
    if (normalized.shock != ShockAuthority::Swcme || !normalized.source.enabled) {
      return Core::Status(Core::StatusCode::ConfigurationConflict,
                          "shock-injection intent requires SWCME shock and enabled source");
    }
  } else if (normalized.source.enabled) {
    return Core::Status(Core::StatusCode::ConfigurationConflict,
                        "enabled source requires shock-injection run intent");
  }
  const ShockOptions& shock = normalized.shockModel;
  if (normalized.inputSchemaVersion < 3 &&
      normalized.shock == ShockAuthority::Swcme &&
      (!std::isfinite(shock.activeFromS) ||
       !std::isfinite(shock.activeUntilS) ||
       !std::isfinite(shock.initialRadiusM) ||
       !std::isfinite(shock.maximumRadiusM) ||
       !std::isfinite(shock.speedMPerS) ||
       !std::isfinite(shock.compressionRatio) ||
       shock.activeUntilS <= shock.activeFromS ||
       shock.initialRadiusM < normalized.innerRadiusM ||
       shock.maximumRadiusM < shock.initialRadiusM ||
       shock.maximumRadiusM > normalized.outerRadiusM ||
       shock.speedMPerS <= 0.0 || shock.compressionRatio <= 1.0)) {
    return Invalid("SWCME shock interval, extent, speed, or compression is invalid");
  }
  const SourceOptions& source = normalized.source;
  if (source.enabled &&
      (!std::isfinite(source.physicalParticleRatePerS) ||
       !std::isfinite(source.injectionEfficiency) ||
       !std::isfinite(source.minimumEnergyJ) ||
       !std::isfinite(source.maximumEnergyJ) ||
       !std::isfinite(source.spectralIndex) ||
       source.physicalParticleRatePerS <= 0.0 ||
       source.injectionEfficiency <= 0.0 || source.injectionEfficiency > 1.0 ||
       source.minimumEnergyJ <= 0.0 ||
       source.maximumEnergyJ <= source.minimumEnergyJ ||
       (normalized.inputSchemaVersion < 3 && source.spectralIndex <= 0.0) ||
       source.samplesPerStep == 0)) {
    return Invalid("source spectrum, efficiency, or sampling controls are invalid");
  }
  if (source.spectrumModel != SourceSpectrumModel::LocalCompressionDsa &&
      source.spectrumModel !=
          SourceSpectrumModel::FixedPhaseSpacePowerLaw) {
    return Invalid("source spectrum model is unknown");
  }
  if (source.spectrumModel == SourceSpectrumModel::LocalCompressionDsa) {
    if (source.fixedPhaseSpacePowerIndex != 0.0) {
      return Invalid("local-compression-dsa source requires zero inactive "
                     "phase_space_power_index");
    }
  } else if (!std::isfinite(source.fixedPhaseSpacePowerIndex) ||
             source.fixedPhaseSpacePowerIndex <= 2.0) {
    // The implementation samples dN/dp = 4*pi*p^2*f(p) over finite positive
    // momentum bounds.  Requiring q>2 gives the intended decreasing number
    // spectrum and rejects the common error of entering the signed exponent
    // -5 instead of the positive q in f proportional to p^(-q).
    return Invalid("fixed-phase-space-power-law source requires a finite "
                   "phase_space_power_index greater than two");
  }
  // The post-compile file owns only a common numerical weight.  Species
  // count, identity, mass, and charge are unavailable in this AMPS-independent
  // factory and are validated against the generated table at the AMPS
  // boundary before any mesh storage is allocated.
  if (!std::isfinite(normalized.species.macroparticleWeight) ||
      normalized.species.macroparticleWeight <= 0.0)
    return Invalid("species.macroparticle_weight must be finite and positive");
  const ParticleNumericsOptions& particleNumerics=normalized.particleNumerics;
  if(particleNumerics.deriveFromMeshAndShock) {
    if(!reducedShockFront||
        particleNumerics.sourceRateModel!=
            SourceRateNormalizationModel::AcceptedShockIncidentFlux||
        !std::isfinite(particleNumerics.maximumParticleSpeedMPerS)||
        particleNumerics.maximumParticleSpeedMPerS<=0||
        particleNumerics.maximumParticleSpeedMPerS>Core::Const::c||
        !std::isfinite(particleNumerics.timeStepMarginFactor)||
        particleNumerics.timeStepMarginFactor<=0||
        particleNumerics.timeStepMarginFactor>1||
        !std::isfinite(particleNumerics.sourceNormalizationRadiusM)||
        particleNumerics.sourceNormalizationRadiusM<=0||
        !std::isfinite(particleNumerics.resolvedMinimumCellSizeM)||
        particleNumerics.resolvedMinimumCellSizeM<0)
      return Invalid("mesh/shock-derived particle numerics require a reduced front, accepted incident-flux model, positive SI speed/radius, and margin in (0,1]");
    if(particleNumerics.resolvedMinimumCellSizeM>0) {
      const double derived=particleNumerics.timeStepMarginFactor*
          particleNumerics.resolvedMinimumCellSizeM/
          particleNumerics.maximumParticleSpeedMPerS;
      if(!NearlyEqual(derived,normalized.requestedTimeStepS))
        return Invalid("global time step differs from margin*minimum_cell_size/maximum_particle_speed");
    }
  } else if(particleNumerics.sourceRateModel!=
          SourceRateNormalizationModel::ConfiguredConstant||
      particleNumerics.maximumParticleSpeedMPerS!=0||
      particleNumerics.timeStepMarginFactor!=0||
      particleNumerics.sourceNormalizationRadiusM!=0||
      particleNumerics.resolvedMinimumCellSizeM!=0||
      !normalized.speciesParticleNormalizations.empty()) {
    return Invalid("inactive derived particle numerics must retain zero sentinels and no per-species records");
  }
  std::vector<int> normalizedSpeciesIndices;
  for(const SpeciesParticleNormalization& item:
      normalized.speciesParticleNormalizations) {
    if(item.ampsIndex<0||item.symbol.empty()||
        !std::isfinite(item.physicalSourceRatePerS)||
        item.physicalSourceRatePerS<=0||
        !std::isfinite(item.macroparticleWeight)||
        item.macroparticleWeight<=0||
        std::find(normalizedSpeciesIndices.begin(),normalizedSpeciesIndices.end(),
            item.ampsIndex)!=normalizedSpeciesIndices.end())
      return Invalid("per-species particle normalization is incomplete or duplicated");
    const double expected=item.physicalSourceRatePerS*
        normalized.requestedTimeStepS/
        static_cast<double>(normalized.source.samplesPerStep);
    if(!NearlyEqual(expected,item.macroparticleWeight))
      return Invalid("per-species particle weight differs from rate*dt/particles_per_iteration");
    normalizedSpeciesIndices.push_back(item.ampsIndex);
  }
  std::vector<std::string> observerIds;
  for (const ObserverOptions& observer : normalized.observers) {
    const double radius = observer.positionM.Norm();
    if (observer.id.empty() || !FiniteVector(observer.positionM) ||
        !std::isfinite(observer.cadenceS) || observer.cadenceS <= 0.0 ||
        observer.energyBins == 0 || observer.pitchAngleBins == 0 ||
        observer.products.empty() ||
        !FiniteVector(observer.velocityMPerS) ||
        !std::isfinite(observer.collectionRadiusM) ||
        observer.collectionRadiusM <= 0.0 ||
        !std::isfinite(observer.shellRadiusM) || observer.shellRadiusM <= 0.0 ||
        !std::isfinite(observer.minimumEnergyJ) ||
        !std::isfinite(observer.maximumEnergyJ) ||
        observer.minimumEnergyJ <= 0.0 ||
        observer.maximumEnergyJ <= observer.minimumEnergyJ ||
        (observer.energyChannelSpacing != EnergyChannelSpacing::Logarithmic &&
         observer.energyChannelSpacing != EnergyChannelSpacing::Linear) ||
        !std::isfinite(observer.minimumMu) ||
        !std::isfinite(observer.maximumMu) ||
        observer.minimumMu < -1.0 || observer.maximumMu > 1.0 ||
        observer.maximumMu <= observer.minimumMu ||
        (observer.allCompiledSpecies && !observer.species.empty()) ||
        (!observer.allCompiledSpecies && observer.species.empty()) ||
        (!observer.followsTrajectory &&
         (radius < normalized.innerRadiusM || radius > normalized.outerRadiusM))) {
      return Invalid("observer identity, location, cadence, bins, or products are invalid");
    }
    const double observerTicks =
        observer.cadenceS / normalized.requestedTimeStepS;
    const double nearestTicks = std::round(observerTicks);
    if (nearestTicks < 1.0 ||
        std::fabs(observerTicks - nearestTicks) >
            64.0 * std::numeric_limits<double>::epsilon() *
                std::max(1.0, std::fabs(observerTicks))) {
      return Invalid(
          "observer cadence must be an integer multiple of requestedTimeStepS");
    }
    if (std::find(observerIds.begin(), observerIds.end(), observer.id) !=
        observerIds.end()) return Invalid("observer IDs must be unique");
    observerIds.push_back(observer.id);

    if (normalized.activeRegion == ActiveRegionMode::ParkerTube) {
      // The active mask is a stationary tube frozen before AMPS allocates any
      // block.  A trajectory supplied later by a moving or coupled observer
      // could leave that tube without a corresponding mesh reactivation, so
      // such a combination is rejected instead of silently publishing empty
      // samples.  Fixed and fixed-heliographic observers are checked as
      // finite collection spheres; a spherical-shell observer is resolved to
      // the same directional point used by ObserverRuntime.
      if (observer.followsTrajectory ||
          observer.kind == ObserverKind::MovingCartesian ||
          observer.kind == ObserverKind::FieldConnected) {
        return Invalid("a stationary Parker active corridor requires fixed "
                       "observers; moving/field-connected observers need a "
                       "future dynamic mesh-reactivation contract");
      }
      Core::Vec3 observerPoint = observer.positionM;
      if (observer.kind == ObserverKind::SphericalShell) {
        Core::Vec3 direction = observer.positionM.Normalized();
        if (direction.NormSq() == 0.0)
          direction = Core::Vec3(1.0, 0.0, 0.0);
        observerPoint = normalized.coordinateOriginM +
            observer.shellRadiusM * direction;
      }
      const Core::Vec3 relative =
          observerPoint - normalized.coordinateOriginM;
      const double observerRadiusM = relative.Norm();
      const double outerLineLengthM = Core::ParkerCurveArcLengthM(
          normalized.outerRadiusM, declaredGeometry);
      const double finiteLineLengthM = std::min(
          normalized.parkerSpiralLengthM, outerLineLengthM);
      double terminalRadiusM = 0.0;
      const Core::Status terminalStatus =
          Core::ParkerCurveRadiusAtArcLengthM(
              finiteLineLengthM, declaredGeometry, &terminalRadiusM);
      if (!terminalStatus.ok()) return terminalStatus;

      double distanceToFiniteTubeM = 0.0;
      double tubeRadiusEvaluationM = observerRadiusM;
      if (observerRadiusM <= terminalRadiusM) {
        // Within the finite line's radial range, retain the same-radius
        // transverse metric used for input validation. The whole-mesh planner
        // subsequently applies the stronger exact segment/AABB test per leaf.
        const Core::Vec3 centreline =
            Core::ParkerCurvePoint(observerRadiusM, declaredGeometry);
        const Core::Vec3 observerDirection = relative.Normalized();
        const Core::Vec3 centrelineDirection = centreline.Normalized();
        distanceToFiniteTubeM = observerRadiusM * std::atan2(
            observerDirection.Cross(centrelineDirection).Norm(),
            std::max(-1.0, std::min(
                1.0, observerDirection.Dot(centrelineDirection))));
      } else {
        // The active mask has a finite end cap. An observer radially beyond it
        // can be accepted only when its collection sphere overlaps that cap;
        // comparing with an infinite analytic continuation would place the
        // observer in storage that the mask intentionally deactivates.
        const Core::Vec3 terminalPoint = normalized.coordinateOriginM +
            Core::ParkerCurvePoint(terminalRadiusM, declaredGeometry);
        distanceToFiniteTubeM = (observerPoint - terminalPoint).Norm();
        tubeRadiusEvaluationM = terminalRadiusM;
      }
      const double activeRadiusM =
          normalized.activeTubeRadiusMode == TubeRadiusMode::PhysicalConstant
              ? normalized.activeTubeRadiusAtReferenceM
              : normalized.activeTubeRadiusAtReferenceM *
                    tubeRadiusEvaluationM /
                    normalized.activeTubeReferenceRadiusM;
      const bool intersectsSolarSphere = normalized.activeSolarSphereRadiusM > 0.0 &&
          observerRadiusM <= normalized.activeSolarSphereRadiusM + observer.collectionRadiusM;
      if (!intersectsSolarSphere && (!std::isfinite(distanceToFiniteTubeM) ||
          distanceToFiniteTubeM > activeRadiusM + observer.collectionRadiusM)) {
        return Invalid("observer '" + observer.id +
                       "' does not intersect the configured finite Parker "
                       "active corridor");
      }
    }

    // The generated AMPS table is intentionally not imported here.  Reject
    // negative indices now; ValidateCompiledSpeciesBinding performs the upper
    // bound check against the generated AMPS count during initialization.
    if (!observer.allCompiledSpecies) {
      for (int species : observer.species)
        if (species < 0)
          return Invalid("observer species indices must be non-negative");
    }
  }

  if (normalized.enablePerpendicularDiffusion || normalized.enableDrifts)
    return Invalid("legacy Boolean perpendicular/drift switches are retired; select the named transport models");
  if (!std::isfinite(normalized.constantKappaPerpendicularM2PerS) ||
      normalized.constantKappaPerpendicularM2PerS < 0.0 ||
      !std::isfinite(normalized.kappaPerpendicularToParallelRatio) ||
      normalized.kappaPerpendicularToParallelRatio < 0.0)
    return Invalid("perpendicular diffusion coefficients must be finite and nonnegative");
  if (normalized.perpendicularDiffusion == PerpendicularDiffusionMode::Constant &&
      normalized.constantKappaPerpendicularM2PerS <= 0.0)
    return Invalid("constant perpendicular diffusion requires positive kappa_perpendicular");
  if (normalized.perpendicularDiffusion == PerpendicularDiffusionMode::ConstantRatio &&
      normalized.kappaPerpendicularToParallelRatio <= 0.0)
    return Invalid("constant-ratio perpendicular diffusion requires a positive ratio");
  // Guiding-centre drift needs grad|B|. Force the storage choice before the
  // layout is frozen so analytic and imported backgrounds use one ABI.
  if (normalized.drift != DriftMode::None)
    normalized.storeMagneticGradient = true;

  const char* reserved = nullptr;
  if (normalized.enableExternalScriptBackground) reserved = "external-script background";
  else if (normalized.enableSelfConsistent3DTurbulence) reserved = "self-consistent 3-D turbulence";
  if (reserved != nullptr) return Core::Status::Reserved(reserved);

  const StorageLayout layout = BuildLayout(normalized);
  // Default legacy profiles retain their existing physical identity. The
  // opt-in extension hashes all new controls and actual shifted root bounds.
  std::ostringstream geometryIdentity;
  geometryIdentity << std::setprecision(17) << std::scientific;
  if (normalized.domainBoxGeometry != DomainBoxGeometry::SunCenteredCube ||
      normalized.parkerSpiralEndRadiusM != 0.0 ||
      normalized.activeSolarSphereRadiusM != 0.0 ||
      normalized.solarRefinementAnchor != SolarRefinementAnchor::SourceShell) {
    geometryIdentity << ";domain_geometry=" << Name(normalized.domainBoxGeometry)
          << ";domain_corner_direction=" << normalized.domainCornerDirection.x << ','
          << normalized.domainCornerDirection.y << ',' << normalized.domainCornerDirection.z
          << ";domain_corner_margin_m=" << normalized.domainCornerMarginM
          << ";domain_min_m=" << domainMinimumM.x << ',' << domainMinimumM.y << ',' << domainMinimumM.z
          << ";domain_max_m=" << domainMaximumM.x << ',' << domainMaximumM.y << ',' << domainMaximumM.z
          << ";line_end_radius_m=" << normalized.parkerSpiralEndRadiusM
          << ";active_solar_sphere_m=" << normalized.activeSolarSphereRadiusM
          << ";solar_anchor=" << Name(normalized.solarRefinementAnchor);
  }
  std::ostringstream physics;
  physics << std::setprecision(17) << std::scientific
          << "sep3d-physics-v8"
          << ";intent=" << Name(normalized.intent)
          << ";background=" << Name(normalized.background)
          << ";turbulence=" << Name(normalized.turbulence)
          << ";shock=" << Name(normalized.shock)
          << ";transport=" << Name(normalized.transport)
          << ";domain=" << Name(normalized.domain)
          << ";outer_mode=" << Name(normalized.outerRadiusMode)
          << ";inner_boundary=" << Name(normalized.innerBoundary)
          << ";outer_boundary=" << Name(normalized.outerBoundary)
          // The physical photosphere is an unconditional AMPS mesh boundary,
          // not an input-selectable transport option.  Record both its
          // implementation contract and exact SI radius nonetheless: restart
          // evidence created before sphere registration must never compare as
          // physics-identical to a run containing the solid solar body.
          << ";solar_boundary=amps-absorbing-sphere-v1"
          << ";solar_radius_m=" << Core::Const::R_sun
          << ";frame=" << normalized.coordinateFrame
          << ";parker_line_start_mode="
          << Name(normalized.parkerSpiralStartMode)
          << ";parker_line_origin=" << normalized.parkerSpiralOriginM.x << ','
          << normalized.parkerSpiralOriginM.y << ','
          << normalized.parkerSpiralOriginM.z
          << ";parker_line_initial=" << normalized.parkerSpiralInitialPointM.x << ','
          << normalized.parkerSpiralInitialPointM.y << ','
          << normalized.parkerSpiralInitialPointM.z
          << ";cme_launch_point_resolved="
          << normalized.cmeLaunchPointResolved;
  if (normalized.cmeLaunchPointResolved) {
    physics << ";cme_launch_point_m=" << normalized.cmeLaunchPointM.x << ','
            << normalized.cmeLaunchPointM.y << ','
            << normalized.cmeLaunchPointM.z;
  }
  physics
          << ";parker_line_length_m=" << normalized.parkerSpiralLengthM
          << ";parker_line_points=" << normalized.parkerSpiralPointCount
          << ";inner_m=" << normalized.innerRadiusM
          << ";outer_m=" << normalized.outerRadiusM
          << ";dt_s=" << normalized.requestedTimeStepS
          << ";maximum_steps=" << normalized.maximumTimeSteps
          << ";seed=" << normalized.campaignSeed
          << ";background_cadence=" << normalized.backgroundCadenceSteps
          << ";injection_cadence=" << normalized.injectionCadenceSteps
          << ";mesh_min_m=" << normalized.minimumCellSizeM
          << ";mesh_background_m=" << normalized.backgroundCellSizeM
          << ";mesh_radial=" << normalized.enableRadialRefinement
          << geometryIdentity.str()
          << ";solar_cell_m=" << normalized.solarSurfaceCellSizeM
          << ";solar_transition_m=" << normalized.solarRefinementOuterRadiusM
          << ";solar_profile=" << Name(normalized.solarRefinementProfile)
          << ";solar_exponent=" << normalized.solarRefinementExponent
          << ";mesh_tube=" << normalized.enableTubeRefinement
          << ";tube_lon=" << normalized.tubeLongitudeRad
          << ";tube_colat=" << normalized.tubeColatitudeRad
          << ";tube_reference_m=" << normalized.tubeReferenceRadiusM
          << ";tube_radius_reference_m=" << normalized.tubeRadiusAtReferenceM
          << ";tube_radius_mode=" << Name(normalized.tubeRadiusMode)
          << ";tube_cell_m=" << normalized.tubeCellSizeM
          << ";tube_profile=" << Name(normalized.tubeTransverseProfile)
          << ";tube_exponent=" << normalized.tubeTransverseExponent
          << ";active_region=" << Name(normalized.activeRegion)
          << ";active_tube_reference_m="
          << normalized.activeTubeReferenceRadiusM
          << ";active_tube_radius_reference_m="
          << normalized.activeTubeRadiusAtReferenceM
          << ";active_tube_radius_mode="
          << Name(normalized.activeTubeRadiusMode)
          << ";active_tube_buffer_blocks="
          << normalized.activeTubeBufferBlocks
          // Mask semantics affect the allocated physical domain and therefore
          // restart compatibility even when every user-facing scalar is
          // unchanged. Bump this literal whenever intersection, topology, or
          // cavity rules change; old evidence must not masquerade as the new
          // hole-free algorithm.
          << ";active_mask_algorithm="
          << (normalized.activeSolarSphereRadiusM > 0.0
                  ? kSolarSphereActiveRegionAlgorithmName : kActiveRegionAlgorithmName)
          << ";mesh_cells_per_block=" << normalized.meshCellsPerBlockEdge
          << ";mesh_max_level=" << normalized.maximumMeshLevel
          << ";mesh_block_overhead=" << normalized.meshBlockOverheadBytes
          << ";mesh_memory_budget=" << normalized.meshMemoryBudgetBytes
          << ";memory_base_cell=" << normalized.memoryModel.baseCellBytes
          << ";memory_base_node=" << normalized.memoryModel.baseNodeBytes
          << ";memory_block=" << normalized.memoryModel.blockStructureBytes
          << ";memory_communication=" << normalized.memoryModel.communicationBytesPerBlock
          << ";memory_particle_bytes=" << normalized.memoryModel.particleBytes
          << ";memory_particles_per_cell=" << normalized.memoryModel.particlesPerCell
          << ";memory_halo_fraction=" << normalized.memoryModel.haloFraction
          << ";memory_safety_fraction=" << normalized.memoryModel.safetyMarginFraction
          << ";parker_reference_m=" << parker.referenceRadiusM
          << ";parker_Br_T=" << parker.radialFieldAtReferenceT
          << ";parker_omega_rad_s=" << parker.solarRotationRateRadPerS
          << ";parker_wind_m_s=" << parker.solarWindSpeedMPerS
          << ";parker_polarity=" << parker.magneticPolarity
          << ";parker_density_m-3=" << parker.numberDensityAtReferenceM3
          << ";parker_density_reference_m="
          << parker.densityReferenceRadiusM
          << ";parker_temperature_K=" << parker.temperatureK
          << ";parker_gamma=" << parker.adiabaticIndex
          << ";parker_thermodynamic_closure="
          << Name(parker.thermodynamicClosure)
          << ";parker_alpha_to_proton=" << parker.alphaToProtonRatio
          << ";parker_electron_temperature_K="
          << parker.electronTemperatureK
          << ";parker_alpha_temperature_K=" << parker.alphaTemperatureK
          << ";parker_reference_sin_colatitude="
          << parker.referenceSinColatitude
          << ";parker_rotation_axis=" << parker.rotationAxis.x << ','
          << parker.rotationAxis.y << ',' << parker.rotationAxis.z
          << ";parker_cadence_s=" << parker.validityCadenceS;
  if (normalized.inputSchemaVersion < 3) {
    physics << ";shock_from_s=" << shock.activeFromS
            << ";shock_until_s=" << shock.activeUntilS
            << ";shock_initial_m=" << shock.initialRadiusM
            << ";shock_maximum_m=" << shock.maximumRadiusM
            << ";shock_speed_m_s=" << shock.speedMPerS
            << ";shock_compression=" << shock.compressionRatio;
  } else {
    // Schema 3 has no surrogate constant-speed ShockOptions.  Its complete
    // kinematics, geometry, and compression are represented only by the
    // canonical SWCME fingerprint below; serializing dormant C++ defaults
    // would falsely make them look like reviewed physics inputs.
    physics << ";shock_model=canonical-swcme3d";
  }
  // Extension selection affects physics even before its provider manifest is
  // available; the common publication boundary also checks that manifest.
  if (normalized.background == BackgroundAuthority::RuntimeModel)
    physics << ";background_model_id=" << normalized.backgroundModelId
            << ";background_model_asset=" << normalized.backgroundModelAssetPath
            << ";background_model_inline_fingerprint="
            << SEP::Background::FingerprintConfiguration(
                normalized.backgroundModelInlineConfiguration)
            << ";background_model_asset_directory="
            << normalized.backgroundModelAssetDirectory;
  physics << ";swcme_fingerprint="
          << normalized.swcmeConfigurationFingerprint
          << ";source_enabled=" << source.enabled
          << ";source_rate_s-1=" << source.physicalParticleRatePerS
          << ";source_efficiency=" << source.injectionEfficiency
          << ";source_min_J=" << source.minimumEnergyJ
          << ";source_max_J=" << source.maximumEnergyJ
          << ";source_spectrum_model=" << Name(source.spectrumModel)
          << ";source_fixed_phase_space_q="
          << source.fixedPhaseSpacePowerIndex;
  if (normalized.inputSchemaVersion < 3)
    physics << ";source_index=" << source.spectralIndex;
  else
    physics << ";source_index=canonical-local-compression";
  physics << ";source_samples_per_compiled_species=" << source.samplesPerStep
          << ";compiled_species_authority=AMPS-SpeciesList"
          << ";species_weight=" << normalized.species.macroparticleWeight
          << ";particle_source_rate_model="
          << Name(normalized.particleNumerics.sourceRateModel)
          << ";particle_numerics_derived="
          << normalized.particleNumerics.deriveFromMeshAndShock
          << ";maximum_particle_speed_m_s="
          << normalized.particleNumerics.maximumParticleSpeedMPerS
          << ";time_step_margin_factor="
          << normalized.particleNumerics.timeStepMarginFactor
          << ";source_normalization_radius_m="
          << normalized.particleNumerics.sourceNormalizationRadiusM
          << ";resolved_minimum_cell_size_m="
          << normalized.particleNumerics.resolvedMinimumCellSizeM;
  for(const SpeciesParticleNormalization& item:
      normalized.speciesParticleNormalizations)
    physics << ";species_normalization=" << item.ampsIndex << ','
            << item.symbol << ',' << item.physicalSourceRatePerS << ','
            << item.macroparticleWeight;
  physics << ";prescribed_turbulence_model="
          << Name(normalized.prescribedTurbulenceModel)
          << ";prescribed_turbulence_amplitude_model="
          << Name(normalized.prescribedTurbulenceAmplitudeModel)
          << ";deltaB_over_B=" << normalized.prescribedDeltaBOverB
          << ";turbulence_wave_energy_reference_J_m-3="
          << normalized.turbulenceWaveEnergyAtReferenceJPerM3
          << ";turbulence_wave_energy_radial_exponent="
          << normalized.turbulenceWaveEnergyRadialExponent
          << ";turbulence_sigma_c="
          << normalized.turbulenceNormalizedCrossHelicity
          << ";turbulence_reference_m="
          << normalized.turbulenceReferenceRadiusM
          << ";turbulence_kmin_m-1=" << normalized.turbulenceKMinPerM
          << ";turbulence_kmax_m-1=" << normalized.turbulenceKMaxPerM
          << ";turbulence_kmin_radial_exponent="
          << normalized.turbulenceKMinRadialExponent
          << ";turbulence_kmax_radial_exponent="
          << normalized.turbulenceKMaxRadialExponent
          << ";turbulence_index=" << normalized.turbulenceSpectralIndex
          << ";turbulence_correlation_m="
          << normalized.turbulenceCorrelationLengthM
          << ";turbulence_correlation_radial_exponent="
          << normalized.turbulenceCorrelationLengthRadialExponent
          << ";turbulence_validity_cadence_s="
          << normalized.turbulenceValidityCadenceS
          << ";missing_turbulence=" << Name(normalized.missingTurbulence)
          << ";resonance_range=" << Name(normalized.resonanceRange)
          << ";cell_crossing_fraction=" << normalized.cellCrossingFraction
          << ";diffusion_fraction=" << normalized.diffusionFraction
          << ";focusing_fraction=" << normalized.focusingFraction
          << ";cooling_fraction=" << normalized.coolingFraction
          << ";field_variation_fraction=" << normalized.fieldVariationFraction
          << ";shock_crossing_fraction=" << normalized.shockCrossingFraction
          << ";minimum_substep_s=" << normalized.minimumTransportSubstepS
          << ";maximum_substeps=" << normalized.maximumTransportSubsteps
          << ";pitch_scheme=" << Name(normalized.pitchAngleScheme)
          << ";spatial_diffusion_model="
          << Name(normalized.spatialDiffusionModel)
          << ";pitch_angle_diffusion_model="
          << Name(normalized.pitchAngleDiffusionModel)
          << ";mean_free_path_model="
          << Name(normalized.meanFreePathModel)
          << ";constant_dmumu_s-1=" << normalized.constantDmumuPerS
          << ";constant_mean_free_path_m="
          << normalized.constantMeanFreePathM
          << ";mean_free_path_reference_m="
          << normalized.meanFreePathReferenceM
          << ";mean_free_path_reference_radius_m="
          << normalized.meanFreePathReferenceRadiusM
          << ";mean_free_path_reference_rigidity_v="
          << normalized.meanFreePathReferenceRigidityV
          << ";mean_free_path_radial_exponent="
          << normalized.meanFreePathRadialExponent
          << ";mean_free_path_rigidity_exponent="
          << normalized.meanFreePathRigidityExponent
          << ";spatial_quadrature_absolute_m2_s="
          << normalized.spatialQuadratureAbsoluteToleranceM2PerS
          << ";spatial_quadrature_relative="
          << normalized.spatialQuadratureRelativeTolerance
          << ";spatial_quadrature_recursion="
          << normalized.spatialQuadratureMaximumRecursion
          << ";focused_scattering_frame="
          << Name(normalized.focusedScatteringFrame)
          << ";maximum_scattering_events_per_substep="
          << normalized.maximumScatteringEventsPerSubstep
          << ";perpendicular_diffusion=" << Name(normalized.perpendicularDiffusion)
          << ";kappa_perpendicular_m2_s="
          << normalized.constantKappaPerpendicularM2PerS
          << ";kappa_perpendicular_ratio="
          << normalized.kappaPerpendicularToParallelRatio
          << ";drift=" << Name(normalized.drift)
          << ";population_control=" << Name(normalized.populationControl)
          << ";population_min_per_cell_species="
          << normalized.minimumParticlesPerCellPerSpecies
          << ";population_target_per_cell_species="
          << normalized.targetParticlesPerCellPerSpecies
          << ";population_max_per_cell_species="
          << normalized.maximumParticlesPerCellPerSpecies
          << ";population_cadence_steps="
          << normalized.populationControlCadenceSteps
          << ";layout=" << layout.fingerprint;
  // Preserve legacy injection/transport fingerprints when the new option is
  // inactive. A propagation target is part of its actual physics identity.
  if (propagation) physics << ";stop_shock_radius_m=" << normalized.stopShockRadiusM;
  for (const ObserverOptions& observer : normalized.observers) {
    physics << ";observer=" << observer.id << ',' << observer.positionM.x
            << ',' << observer.positionM.y << ',' << observer.positionM.z
            << ',' << observer.followsTrajectory << ',' << observer.cadenceS
            << ',' << observer.energyBins << ',' << observer.pitchAngleBins
            << ',' << observer.products << ',' << Name(observer.kind)
            << ',' << Name(observer.normalization)
            << ',' << observer.velocityMPerS.x << ',' << observer.velocityMPerS.y
            << ',' << observer.velocityMPerS.z
            << ',' << observer.collectionRadiusM << ',' << observer.shellRadiusM
            << ',' << observer.minimumEnergyJ << ',' << observer.maximumEnergyJ
            << ',' << Name(observer.energyChannelSpacing)
            << ',' << observer.minimumMu << ',' << observer.maximumMu;
    if (observer.allCompiledSpecies) {
      // Preserve wildcard intent in restart/physics identity.  It must not be
      // fingerprinted as an accidentally empty observer selection.
      physics << ",all-compiled-species";
    } else {
      for (int species : observer.species) physics << ',' << species;
    }
  }
  const std::string fingerprint =
      SEP::Background::FingerprintConfiguration(physics.str());

  // A restart segment normally has a different --output-dir and necessarily
  // has a non-empty restart input path.  Comparing the complete provenance
  // manifest would therefore reject every real restart before mesh creation.
  // Keep a separate compatibility identity: resolved physical inputs and the
  // two integer event clocks remain frozen, while paths and output naming are
  // allowed to relocate.  Physics and native storage are independently
  // guarded by physicsFingerprint and StorageLayout::fingerprint.
  std::ostringstream restartCompatibility;
  restartCompatibility << physics.str()
           << ";swcme_resolved_manifest="
           << normalized.swcmeResolvedManifest
           << ";output_cadence=" << normalized.outputCadenceSteps
           << ";checkpoint_cadence=" << normalized.checkpointCadenceSteps;

  std::ostringstream manifest;
  manifest << restartCompatibility.str()
           << ";output_directory=" << normalized.outputDirectory
           << ";output_prefix=" << normalized.outputPrefix
           << ";initialization_mesh_tecplot="
           << normalized.initializationMeshTecplotFile
           << ";initialization_parker_line_tecplot="
           << normalized.initializationParkerLineTecplotFile
           << ";initialization_data_tecplot="
           << normalized.initializationDataTecplotFile
           << ";restart_input=" << normalized.restartInputPath
           << ";restart_output=" << normalized.restartOutputPath;

  configuration->reset(new RunConfiguration3D(
      normalized, layout, fingerprint, manifest.str(),
      restartCompatibility.str()));
  return Core::Status::OK();
}

}  // namespace RuntimeModel
}  // namespace SEP3D
