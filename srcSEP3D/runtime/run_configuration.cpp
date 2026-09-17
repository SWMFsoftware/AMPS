#include "run_configuration.h"

#include "sep_background_snapshot.h"

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
    case DomainPreset::Earth: return "earth";
    case DomainPreset::Mars: return "mars";
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

  if (!std::isfinite(options.innerRadiusM) || options.innerRadiusM <= 0.0) {
    return Invalid("innerRadiusM must be finite and positive");
  }
  if (!std::isfinite(options.outerRadiusM) ||
      options.outerRadiusM <= options.innerRadiusM) {
    return Invalid("outerRadiusM must be finite and exceed innerRadiusM");
  }
  if (!std::isfinite(options.requestedTimeStepS) ||
      options.requestedTimeStepS <= 0.0) {
    return Invalid("requestedTimeStepS must be finite and positive");
  }
  if (options.campaignSeed == 0) return Invalid("campaignSeed zero is reserved");
  if (options.backgroundCadenceSteps == 0) {
    return Invalid("backgroundCadenceSteps must be positive");
  }
  if (options.outputCadenceSteps == 0) {
    return Invalid("outputCadenceSteps must be positive");
  }
  if (options.outputDirectory.empty() || options.outputPrefix.empty()) {
    return Invalid("output directory and prefix must not be empty");
  }
  const double meshValues[] = {
      options.minimumCellSizeM, options.backgroundCellSizeM,
      options.tubeLongitudeRad, options.tubeColatitudeRad,
      options.tubeCoreRadiusM, options.tubeShoulderRadiusM,
      options.tubeCellSizeM};
  for (double value : meshValues) {
    if (!std::isfinite(value)) return Invalid("mesh option is not finite");
  }
  if (options.minimumCellSizeM <= 0.0 ||
      options.backgroundCellSizeM < options.minimumCellSizeM ||
      options.meshCellsPerBlockEdge == 0 || options.maximumMeshLevel > 19 ||
      options.meshMemoryBudgetBytes == 0) {
    return Invalid("mesh sizes, level, block width, or memory budget are invalid");
  }
  if (options.enableTubeRefinement &&
      ((options.tubePolarity != 1 && options.tubePolarity != -1) ||
       options.tubeCoreRadiusM <= 0.0 ||
       options.tubeShoulderRadiusM < options.tubeCoreRadiusM ||
       options.tubeCellSizeM <= 0.0 ||
       options.tubeColatitudeRad < 0.0 ||
       options.tubeColatitudeRad > Core::Const::kPi)) {
    return Invalid("Parker-tube refinement options are invalid");
  }
  if (!std::isfinite(options.prescribedDeltaBOverB) ||
      options.prescribedDeltaBOverB <= 0.0 ||
      !std::isfinite(options.turbulenceKMinPerM) ||
      options.turbulenceKMinPerM <= 0.0 ||
      !std::isfinite(options.turbulenceKMaxPerM) ||
      options.turbulenceKMaxPerM <= options.turbulenceKMinPerM ||
      !std::isfinite(options.turbulenceSpectralIndex) ||
      options.turbulenceSpectralIndex <= 1.0 ||
      !std::isfinite(options.turbulenceCorrelationLengthM) ||
      options.turbulenceCorrelationLengthM <= 0.0) {
    return Invalid("turbulence amplitude, band, index, or correlation length is invalid");
  }
  if (options.turbulence == TurbulenceAuthority::Swmf &&
      options.background != BackgroundAuthority::Swmf) {
    return Core::Status(
        Core::StatusCode::ConfigurationConflict,
        "SWMF turbulence requires the SWMF background authority");
  }

  const char* reserved = nullptr;
  if (options.enablePerpendicularDiffusion) reserved = "perpendicular diffusion";
  else if (options.enableDrifts) reserved = "gradient/curvature drifts";
  else if (options.enableExternalScriptBackground) reserved = "external-script background";
  else if (options.enableSelfConsistent3DTurbulence) reserved = "self-consistent 3-D turbulence";
  if (reserved != nullptr) return Core::Status::Reserved(reserved);

  const StorageLayout layout = BuildLayout(options);
  std::ostringstream physics;
  physics << std::setprecision(17) << std::scientific
          << "sep3d-physics-v1"
          << ";background=" << Name(options.background)
          << ";turbulence=" << Name(options.turbulence)
          << ";shock=" << Name(options.shock)
          << ";transport=" << Name(options.transport)
          << ";domain=" << Name(options.domain)
          << ";inner_m=" << options.innerRadiusM
          << ";outer_m=" << options.outerRadiusM
          << ";dt_s=" << options.requestedTimeStepS
          << ";seed=" << options.campaignSeed
          << ";background_cadence=" << options.backgroundCadenceSteps
          << ";mesh_min_m=" << options.minimumCellSizeM
          << ";mesh_background_m=" << options.backgroundCellSizeM
          << ";mesh_radial=" << options.enableRadialRefinement
          << ";mesh_tube=" << options.enableTubeRefinement
          << ";tube_lon=" << options.tubeLongitudeRad
          << ";tube_colat=" << options.tubeColatitudeRad
          << ";tube_polarity=" << options.tubePolarity
          << ";tube_core_m=" << options.tubeCoreRadiusM
          << ";tube_shoulder_m=" << options.tubeShoulderRadiusM
          << ";tube_cell_m=" << options.tubeCellSizeM
          << ";mesh_cells_per_block=" << options.meshCellsPerBlockEdge
          << ";mesh_max_level=" << options.maximumMeshLevel
          << ";mesh_block_overhead=" << options.meshBlockOverheadBytes
          << ";mesh_memory_budget=" << options.meshMemoryBudgetBytes
          << ";deltaB_over_B=" << options.prescribedDeltaBOverB
          << ";turbulence_kmin_m-1=" << options.turbulenceKMinPerM
          << ";turbulence_kmax_m-1=" << options.turbulenceKMaxPerM
          << ";turbulence_index=" << options.turbulenceSpectralIndex
          << ";turbulence_correlation_m="
          << options.turbulenceCorrelationLengthM
          << ";missing_turbulence=" << Name(options.missingTurbulence)
          << ";resonance_range=" << Name(options.resonanceRange)
          << ";layout=" << layout.fingerprint;
  const std::string fingerprint =
      SEP::Background::FingerprintConfiguration(physics.str());

  std::ostringstream manifest;
  manifest << physics.str()
           << ";output_cadence=" << options.outputCadenceSteps
           << ";output_directory=" << options.outputDirectory
           << ";output_prefix=" << options.outputPrefix;

  configuration->reset(new RunConfiguration3D(
      options, layout, fingerprint, manifest.str()));
  return Core::Status::OK();
}

}  // namespace RuntimeModel
}  // namespace SEP3D
