#include "run_configuration.h"

#include "sep_background_snapshot.h"

#include <algorithm>
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
  if (options.turbulence == TurbulenceAuthority::Swmf) {
    AppendField(2, &cursor, &layout.waveEnergyOffset);
  }
  layout.cellAssociatedBytes = cursor;
  layout.samplingBytesPerCell = options.samplingBytesPerCell;

  std::ostringstream canonical;
  canonical << "sep3d-storage-v1"
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

}  // namespace

const char* Name(BackgroundAuthority value) {
  switch (value) {
    case BackgroundAuthority::AnalyticParker: return "analytic-parker";
    case BackgroundAuthority::Swmf: return "swmf";
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

const char* Name(TransportModel value) {
  switch (value) {
    case TransportModel::Parker3D: return "parker3d";
    case TransportModel::Focused3D: return "focused3d";
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
  }
  return "unknown";
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

const char* Name(RunIntent value) {
  switch (value) {
    case RunIntent::TransportOnly: return "transport-only";
    case RunIntent::ShockInjection: return "shock-injection";
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

RunConfiguration3D::RunConfiguration3D(
    const RunConfiguration3DOptions& options, const StorageLayout& layout,
    const std::string& physicsFingerprint,
    const std::string& resolvedManifest)
    : options_(options),
      storageLayout_(layout),
      physicsFingerprint_(physicsFingerprint),
      resolvedManifest_(resolvedManifest) {}

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
  if (normalized.outerRadiusMode == OuterRadiusMode::Preset) {
    normalized.outerRadiusM = PresetOuterRadiusM(normalized.domain);
  }
  if (normalized.domain == DomainPreset::Earth) {
    normalized.domain = DomainPreset::OneAu;
  }
  normalized.parker.sourceRadiusM = normalized.innerRadiusM;
  normalized.parker.sourceLongitudeRad = normalized.tubeLongitudeRad;
  normalized.parker.sourceColatitudeRad = normalized.tubeColatitudeRad;

  if (!std::isfinite(normalized.innerRadiusM) || normalized.innerRadiusM <= 0.0) {
    return Invalid("innerRadiusM must be finite and positive");
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
  const double meshValues[] = {
      normalized.minimumCellSizeM, normalized.backgroundCellSizeM,
      normalized.solarSurfaceCellSizeM,
      normalized.solarRefinementOuterRadiusM,
      normalized.tubeLongitudeRad, normalized.tubeColatitudeRad,
      normalized.tubeReferenceRadiusM,
      normalized.tubeRadiusAtReferenceM, normalized.tubeCellSizeM};
  for (double value : meshValues) {
    if (!std::isfinite(value)) return Invalid("mesh option is not finite");
  }
  if (normalized.minimumCellSizeM <= 0.0 ||
      normalized.backgroundCellSizeM < normalized.minimumCellSizeM ||
      normalized.solarSurfaceCellSizeM < normalized.minimumCellSizeM ||
      normalized.solarSurfaceCellSizeM > normalized.backgroundCellSizeM ||
      normalized.solarRefinementOuterRadiusM <= normalized.innerRadiusM ||
      normalized.solarRefinementOuterRadiusM > normalized.outerRadiusM ||
      !ValidProfileExponent(normalized.solarRefinementExponent) ||
      !ValidProfileExponent(normalized.tubeTransverseExponent) ||
      normalized.meshCellsPerBlockEdge == 0 ||
      normalized.maximumMeshLevel > 19 ||
      normalized.meshMemoryBudgetBytes == 0) {
    return Invalid("mesh sizes, level, block width, or memory budget are invalid");
  }
  if (normalized.enableTubeRefinement &&
      (normalized.tubeReferenceRadiusM <= normalized.innerRadiusM ||
       normalized.tubeRadiusAtReferenceM <= 0.0 ||
       normalized.tubeCellSizeM < normalized.minimumCellSizeM ||
       normalized.tubeCellSizeM > normalized.backgroundCellSizeM ||
       normalized.tubeColatitudeRad < 0.0 ||
       normalized.tubeColatitudeRad > Core::Const::kPi)) {
    return Invalid("Parker-tube refinement options are invalid");
  }
  const double rootCellM = 2.0 * normalized.outerRadiusM /
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
  if (!std::isfinite(normalized.prescribedDeltaBOverB) ||
      normalized.prescribedDeltaBOverB <= 0.0 ||
      !std::isfinite(normalized.turbulenceKMinPerM) ||
      normalized.turbulenceKMinPerM <= 0.0 ||
      !std::isfinite(normalized.turbulenceKMaxPerM) ||
      normalized.turbulenceKMaxPerM <= normalized.turbulenceKMinPerM ||
      !std::isfinite(normalized.turbulenceSpectralIndex) ||
      normalized.turbulenceSpectralIndex <= 1.0 ||
      !std::isfinite(normalized.turbulenceCorrelationLengthM) ||
      normalized.turbulenceCorrelationLengthM <= 0.0) {
    return Invalid("turbulence amplitude, band, index, or correlation length is invalid");
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

  const ParkerPhysicsOptions& parker = normalized.parker;
  const double parkerValues[] = {
      parker.sourceRadiusM, parker.sourceLongitudeRad,
      parker.sourceColatitudeRad, parker.referenceRadiusM,
      parker.radialFieldAtReferenceT, parker.solarRotationRateRadPerS,
      parker.solarWindSpeedMPerS, parker.numberDensityAtReferenceM3,
      parker.temperatureK, parker.validityCadenceS};
  for (double value : parkerValues)
    if (!std::isfinite(value)) return Invalid("Parker configuration contains a non-finite value");
  if (parker.sourceRadiusM != normalized.innerRadiusM ||
      parker.referenceRadiusM <= parker.sourceRadiusM ||
      parker.radialFieldAtReferenceT <= 0.0 ||
      parker.solarWindSpeedMPerS <= 0.0 ||
      parker.numberDensityAtReferenceM3 <= 0.0 || parker.temperatureK <= 0.0 ||
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
  if (normalized.shock == ShockAuthority::Swcme &&
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
       source.spectralIndex <= 0.0 || source.samplesPerStep == 0)) {
    return Invalid("source spectrum, efficiency, or sampling controls are invalid");
  }
  if (normalized.species.name.empty() ||
      !std::isfinite(normalized.species.massKg) ||
      !std::isfinite(normalized.species.chargeC) ||
      !std::isfinite(normalized.species.macroparticleWeight) ||
      normalized.species.massKg <= 0.0 ||
      normalized.species.chargeC == 0.0 ||
      normalized.species.macroparticleWeight <= 0.0) {
    return Invalid("species name, mass, charge, or particle weight is invalid");
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
        !std::isfinite(observer.minimumMu) ||
        !std::isfinite(observer.maximumMu) ||
        observer.minimumMu < -1.0 || observer.maximumMu > 1.0 ||
        observer.maximumMu <= observer.minimumMu || observer.species.empty() ||
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
    for (int species : observer.species)
      if (species < 0) return Invalid("observer species index is negative");
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
  std::ostringstream physics;
  physics << std::setprecision(17) << std::scientific
          << "sep3d-physics-v3"
          << ";intent=" << Name(normalized.intent)
          << ";background=" << Name(normalized.background)
          << ";turbulence=" << Name(normalized.turbulence)
          << ";shock=" << Name(normalized.shock)
          << ";transport=" << Name(normalized.transport)
          << ";domain=" << Name(normalized.domain)
          << ";outer_mode=" << Name(normalized.outerRadiusMode)
          << ";inner_boundary=" << Name(normalized.innerBoundary)
          << ";outer_boundary=" << Name(normalized.outerBoundary)
          << ";frame=" << normalized.coordinateFrame
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
          << ";parker_temperature_K=" << parker.temperatureK
          << ";parker_cadence_s=" << parker.validityCadenceS
          << ";shock_from_s=" << shock.activeFromS
          << ";shock_until_s=" << shock.activeUntilS
          << ";shock_initial_m=" << shock.initialRadiusM
          << ";shock_maximum_m=" << shock.maximumRadiusM
          << ";shock_speed_m_s=" << shock.speedMPerS
          << ";shock_compression=" << shock.compressionRatio
          << ";source_enabled=" << source.enabled
          << ";source_rate_s-1=" << source.physicalParticleRatePerS
          << ";source_efficiency=" << source.injectionEfficiency
          << ";source_min_J=" << source.minimumEnergyJ
          << ";source_max_J=" << source.maximumEnergyJ
          << ";source_index=" << source.spectralIndex
          << ";source_samples=" << source.samplesPerStep
          << ";species_name=" << normalized.species.name
          << ";species_mass_kg=" << normalized.species.massKg
          << ";species_charge_C=" << normalized.species.chargeC
          << ";species_weight=" << normalized.species.macroparticleWeight
          << ";deltaB_over_B=" << normalized.prescribedDeltaBOverB
          << ";turbulence_kmin_m-1=" << normalized.turbulenceKMinPerM
          << ";turbulence_kmax_m-1=" << normalized.turbulenceKMaxPerM
          << ";turbulence_index=" << normalized.turbulenceSpectralIndex
          << ";turbulence_correlation_m="
          << normalized.turbulenceCorrelationLengthM
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
          << ";perpendicular_diffusion=" << Name(normalized.perpendicularDiffusion)
          << ";kappa_perpendicular_m2_s="
          << normalized.constantKappaPerpendicularM2PerS
          << ";kappa_perpendicular_ratio="
          << normalized.kappaPerpendicularToParallelRatio
          << ";drift=" << Name(normalized.drift)
          << ";layout=" << layout.fingerprint;
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
            << ',' << observer.minimumMu << ',' << observer.maximumMu;
    for (int species : observer.species) physics << ',' << species;
  }
  const std::string fingerprint =
      SEP::Background::FingerprintConfiguration(physics.str());

  std::ostringstream manifest;
  manifest << physics.str()
           << ";output_cadence=" << normalized.outputCadenceSteps
           << ";checkpoint_cadence=" << normalized.checkpointCadenceSteps
           << ";output_directory=" << normalized.outputDirectory
           << ";output_prefix=" << normalized.outputPrefix
           << ";restart_input=" << normalized.restartInputPath
           << ";restart_output=" << normalized.restartOutputPath;

  configuration->reset(new RunConfiguration3D(
      normalized, layout, fingerprint, manifest.str()));
  return Core::Status::OK();
}

}  // namespace RuntimeModel
}  // namespace SEP3D
