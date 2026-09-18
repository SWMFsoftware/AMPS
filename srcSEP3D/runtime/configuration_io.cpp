#include "configuration_io.h"

#include "../mesh/mesh_model.h"

#include <algorithm>
#include <cctype>
#include <cerrno>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <utility>

namespace SEP3D {
namespace RuntimeModel {
namespace {

Core::Status Invalid(const std::string& message) {
  return Core::Status(Core::StatusCode::InvalidInput, message);
}

std::string Trim(const std::string& value) {
  const std::size_t first = value.find_first_not_of(" \t\r\n");
  if (first == std::string::npos) return {};
  const std::size_t last = value.find_last_not_of(" \t\r\n");
  return value.substr(first, last - first + 1);
}

std::string Lower(std::string value) {
  std::transform(value.begin(), value.end(), value.begin(),
                 [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
  return value;
}

bool ParseDouble(const std::string& text, double* value) {
  if (value == nullptr || text.empty()) return false;
  errno = 0;
  char* end = nullptr;
  const double parsed = std::strtod(text.c_str(), &end);
  if (errno == ERANGE || end == text.c_str() || *end != '\0' ||
      !std::isfinite(parsed)) return false;
  *value = parsed;
  return true;
}

bool ParseUnsigned64(const std::string& text, std::uint64_t* value) {
  if (value == nullptr || text.empty() || text[0] == '-') return false;
  errno = 0;
  char* end = nullptr;
  const unsigned long long parsed = std::strtoull(text.c_str(), &end, 10);
  if (errno == ERANGE || end == text.c_str() || *end != '\0') return false;
  *value = static_cast<std::uint64_t>(parsed);
  return true;
}

bool ParseUnsigned(const std::string& text, unsigned* value) {
  std::uint64_t parsed = 0;
  if (!ParseUnsigned64(text, &parsed) ||
      parsed > std::numeric_limits<unsigned>::max()) return false;
  *value = static_cast<unsigned>(parsed);
  return true;
}

bool ParseSize(const std::string& text, std::size_t* value) {
  std::uint64_t parsed = 0;
  if (!ParseUnsigned64(text, &parsed) ||
      parsed > std::numeric_limits<std::size_t>::max()) return false;
  *value = static_cast<std::size_t>(parsed);
  return true;
}

bool ParseSpeciesList(const std::string& text, std::vector<int>* values) {
  if (values == nullptr) return false;
  std::vector<int> parsed;
  std::istringstream input(text);
  std::string token;
  while (std::getline(input, token, ',')) {
    token = Trim(token);
    if (token.empty() || token[0] == '-') return false;
    errno = 0;
    char* end = nullptr;
    const long value = std::strtol(token.c_str(), &end, 10);
    if (errno == ERANGE || end == token.c_str() || *end != '\0' ||
        value > std::numeric_limits<int>::max()) return false;
    parsed.push_back(static_cast<int>(value));
  }
  if (parsed.empty()) return false;
  *values = std::move(parsed);
  return true;
}

bool ParseBool(const std::string& text, bool* value) {
  const std::string normalized = Lower(text);
  if (normalized == "true" || normalized == "yes" || normalized == "on") {
    *value = true;
    return true;
  }
  if (normalized == "false" || normalized == "no" || normalized == "off") {
    *value = false;
    return true;
  }
  return false;
}

bool KnownSection(const std::string& section) {
  static const std::set<std::string> fixed = {
      "run", "domain", "mesh", "mesh.solar", "mesh.tube", "memory",
      "background", "background.parker", "turbulence", "transport",
      "shock", "source", "species", "storage", "output", "restart"};
  return fixed.count(section) != 0 ||
         (section.rfind("observer.", 0) == 0 &&
          section.size() > std::string("observer.").size());
}

template <typename Enum>
bool ParseEnum(const std::string& text,
               const std::map<std::string, Enum>& values, Enum* result) {
  const auto found = values.find(Lower(text));
  if (found == values.end()) return false;
  *result = found->second;
  return true;
}

ObserverOptions* ObserverForSection(
    const std::string& section, RunConfiguration3DOptions* options,
    bool* observerDefaultsCleared) {
  const std::string prefix = "observer.";
  if (section.rfind(prefix, 0) != 0 || section.size() == prefix.size())
    return nullptr;
  if (!*observerDefaultsCleared) {
    options->observers.clear();
    *observerDefaultsCleared = true;
  }
  const std::string id = section.substr(prefix.size());
  for (ObserverOptions& observer : options->observers)
    if (observer.id == id) return &observer;
  ObserverOptions observer;
  observer.id = id;
  options->observers.push_back(observer);
  return &options->observers.back();
}

Core::Status ApplyField(const std::string& section, const std::string& key,
                        const std::string& value,
                        RunConfiguration3DOptions* o,
                        bool* observerDefaultsCleared) {
  const std::string field = section.empty() ? key : section + "." + key;
  auto invalidValue = [&]() {
    return Invalid("invalid value for '" + field + "': " + value);
  };

  if (field == "run.schema_version") {
    std::uint64_t version = 0;
    if (!ParseUnsigned64(value, &version) || version != 1) return invalidValue();
    return Core::Status::OK();
  }
  if (field == "run.intent") {
    if (!ParseEnum(value, {{"transport-only", RunIntent::TransportOnly},
                           {"shock-injection", RunIntent::ShockInjection}},
                   &o->intent)) return invalidValue();
  } else if (field == "run.transport") {
    if (!ParseEnum(value, {{"parker3d", TransportModel::Parker3D},
                           {"focused3d", TransportModel::Focused3D}},
                   &o->transport)) return invalidValue();
  } else if (field == "run.time_step_s") {
    if (!ParseDouble(value, &o->requestedTimeStepS)) return invalidValue();
  } else if (field == "run.maximum_time_steps") {
    if (!ParseUnsigned64(value, &o->maximumTimeSteps)) return invalidValue();
  } else if (field == "run.campaign_seed") {
    if (!ParseUnsigned64(value, &o->campaignSeed)) return invalidValue();
  } else if (field == "run.background_cadence_steps") {
    if (!ParseUnsigned64(value, &o->backgroundCadenceSteps)) return invalidValue();
  } else if (field == "run.injection_cadence_steps") {
    if (!ParseUnsigned64(value, &o->injectionCadenceSteps)) return invalidValue();
  } else if (field == "domain.preset") {
    if (!ParseEnum(value, {{"solar", DomainPreset::Solar},
                           {"one-au", DomainPreset::OneAu},
                           {"earth", DomainPreset::OneAu},
                           {"mars", DomainPreset::Mars}}, &o->domain))
      return invalidValue();
  } else if (field == "domain.inner_radius_m") {
    if (!ParseDouble(value, &o->innerRadiusM)) return invalidValue();
  } else if (field == "domain.inner_boundary") {
    if (!ParseEnum(value, {{"absorb", InnerBoundaryMode::Absorb}},
                   &o->innerBoundary)) return invalidValue();
  } else if (field == "domain.outer_radius_mode") {
    if (!ParseEnum(value, {{"preset", OuterRadiusMode::Preset},
                           {"explicit", OuterRadiusMode::Explicit}},
                   &o->outerRadiusMode)) return invalidValue();
  } else if (field == "domain.outer_radius_m") {
    if (!ParseDouble(value, &o->outerRadiusM)) return invalidValue();
  } else if (field == "domain.outer_boundary") {
    if (!ParseEnum(value, {{"escape", OuterBoundaryMode::Escape},
                           {"imported-coverage", OuterBoundaryMode::ImportedCoverage}},
                   &o->outerBoundary)) return invalidValue();
  } else if (field == "domain.coordinate_frame") {
    o->coordinateFrame = value;
    o->parker.coordinateFrame = value;
  } else if (field == "domain.origin_x_m") {
    if (!ParseDouble(value, &o->coordinateOriginM.x)) return invalidValue();
  } else if (field == "domain.origin_y_m") {
    if (!ParseDouble(value, &o->coordinateOriginM.y)) return invalidValue();
  } else if (field == "domain.origin_z_m") {
    if (!ParseDouble(value, &o->coordinateOriginM.z)) return invalidValue();
  } else if (field == "mesh.global_cell_size_m") {
    if (!ParseDouble(value, &o->backgroundCellSizeM)) return invalidValue();
  } else if (field == "mesh.minimum_cell_size_m") {
    if (!ParseDouble(value, &o->minimumCellSizeM)) return invalidValue();
  } else if (field == "mesh.cells_per_block_edge") {
    if (!ParseUnsigned(value, &o->meshCellsPerBlockEdge)) return invalidValue();
  } else if (field == "mesh.maximum_level") {
    if (!ParseUnsigned(value, &o->maximumMeshLevel)) return invalidValue();
  } else if (field == "mesh.memory_budget_bytes") {
    if (!ParseSize(value, &o->meshMemoryBudgetBytes)) return invalidValue();
  } else if (field == "mesh.block_overhead_bytes") {
    if (!ParseSize(value, &o->meshBlockOverheadBytes)) return invalidValue();
  } else if (field == "mesh.solar.enabled") {
    if (!ParseBool(value, &o->enableRadialRefinement)) return invalidValue();
  } else if (field == "mesh.solar.surface_cell_size_m") {
    if (!ParseDouble(value, &o->solarSurfaceCellSizeM)) return invalidValue();
  } else if (field == "mesh.solar.transition_outer_radius_m") {
    if (!ParseDouble(value, &o->solarRefinementOuterRadiusM)) return invalidValue();
  } else if (field == "mesh.solar.profile") {
    if (!ParseEnum(value, {{"linear", RefinementProfile::Linear},
                           {"power-law", RefinementProfile::PowerLaw},
                           {"smoothstep", RefinementProfile::Smoothstep}},
                   &o->solarRefinementProfile)) return invalidValue();
  } else if (field == "mesh.solar.exponent") {
    if (!ParseDouble(value, &o->solarRefinementExponent)) return invalidValue();
  } else if (field == "mesh.tube.enabled") {
    if (!ParseBool(value, &o->enableTubeRefinement)) return invalidValue();
  } else if (field == "mesh.tube.source_longitude_rad") {
    if (!ParseDouble(value, &o->tubeLongitudeRad)) return invalidValue();
  } else if (field == "mesh.tube.source_colatitude_rad") {
    if (!ParseDouble(value, &o->tubeColatitudeRad)) return invalidValue();
  } else if (field == "mesh.tube.reference_radius_m") {
    if (!ParseDouble(value, &o->tubeReferenceRadiusM)) return invalidValue();
  } else if (field == "mesh.tube.radius_at_reference_m") {
    if (!ParseDouble(value, &o->tubeRadiusAtReferenceM)) return invalidValue();
  } else if (field == "mesh.tube.radius_mode") {
    if (!ParseEnum(value,
        {{"physical-constant", TubeRadiusMode::PhysicalConstant},
         {"constant-angular-width", TubeRadiusMode::ConstantAngularWidth}},
        &o->tubeRadiusMode)) return invalidValue();
  } else if (field == "mesh.tube.center_cell_size_m") {
    if (!ParseDouble(value, &o->tubeCellSizeM)) return invalidValue();
  } else if (field == "mesh.tube.transverse_profile") {
    if (!ParseEnum(value, {{"linear", RefinementProfile::Linear},
                           {"power-law", RefinementProfile::PowerLaw},
                           {"smoothstep", RefinementProfile::Smoothstep}},
                   &o->tubeTransverseProfile)) return invalidValue();
  } else if (field == "mesh.tube.transverse_exponent") {
    if (!ParseDouble(value, &o->tubeTransverseExponent)) return invalidValue();
  } else if (field == "memory.base_cell_bytes") {
    if (!ParseSize(value, &o->memoryModel.baseCellBytes)) return invalidValue();
  } else if (field == "memory.base_node_bytes") {
    if (!ParseSize(value, &o->memoryModel.baseNodeBytes)) return invalidValue();
  } else if (field == "memory.block_structure_bytes") {
    if (!ParseSize(value, &o->memoryModel.blockStructureBytes)) return invalidValue();
  } else if (field == "memory.communication_bytes_per_block") {
    if (!ParseSize(value, &o->memoryModel.communicationBytesPerBlock)) return invalidValue();
  } else if (field == "memory.particle_bytes") {
    if (!ParseSize(value, &o->memoryModel.particleBytes)) return invalidValue();
  } else if (field == "memory.particles_per_cell") {
    if (!ParseDouble(value, &o->memoryModel.particlesPerCell)) return invalidValue();
  } else if (field == "memory.halo_fraction") {
    if (!ParseDouble(value, &o->memoryModel.haloFraction)) return invalidValue();
  } else if (field == "memory.safety_margin_fraction") {
    if (!ParseDouble(value, &o->memoryModel.safetyMarginFraction)) return invalidValue();
  } else if (field == "background.provider") {
    if (!ParseEnum(value,
        {{"analytic-parker", BackgroundAuthority::AnalyticParker},
         {"swmf", BackgroundAuthority::Swmf}}, &o->background))
      return invalidValue();
  } else if (field == "background.external_script") {
    if (!ParseBool(value, &o->enableExternalScriptBackground)) return invalidValue();
  } else if (field == "background.parker.reference_radius_m") {
    if (!ParseDouble(value, &o->parker.referenceRadiusM)) return invalidValue();
  } else if (field == "background.parker.radial_field_at_reference_t") {
    if (!ParseDouble(value, &o->parker.radialFieldAtReferenceT)) return invalidValue();
  } else if (field == "background.parker.solar_rotation_rate_rad_per_s") {
    if (!ParseDouble(value, &o->parker.solarRotationRateRadPerS)) return invalidValue();
  } else if (field == "background.parker.solar_wind_speed_m_per_s") {
    if (!ParseDouble(value, &o->parker.solarWindSpeedMPerS)) return invalidValue();
  } else if (field == "background.parker.magnetic_polarity") {
    double parsed = 0.0;
    if (!ParseDouble(value, &parsed) || (parsed != -1.0 && parsed != 1.0))
      return invalidValue();
    o->parker.magneticPolarity = static_cast<int>(parsed);
  } else if (field == "background.parker.number_density_at_reference_m3") {
    if (!ParseDouble(value, &o->parker.numberDensityAtReferenceM3)) return invalidValue();
  } else if (field == "background.parker.temperature_k") {
    if (!ParseDouble(value, &o->parker.temperatureK)) return invalidValue();
  } else if (field == "background.parker.validity_cadence_s") {
    if (!ParseDouble(value, &o->parker.validityCadenceS)) return invalidValue();
  } else if (field == "turbulence.authority") {
    if (!ParseEnum(value, {{"prescribed", TurbulenceAuthority::Prescribed},
                           {"swmf", TurbulenceAuthority::Swmf}},
                   &o->turbulence)) return invalidValue();
  } else if (field == "turbulence.delta_b_over_b") {
    if (!ParseDouble(value, &o->prescribedDeltaBOverB)) return invalidValue();
  } else if (field == "turbulence.k_min_per_m") {
    if (!ParseDouble(value, &o->turbulenceKMinPerM)) return invalidValue();
  } else if (field == "turbulence.k_max_per_m") {
    if (!ParseDouble(value, &o->turbulenceKMaxPerM)) return invalidValue();
  } else if (field == "turbulence.spectral_index") {
    if (!ParseDouble(value, &o->turbulenceSpectralIndex)) return invalidValue();
  } else if (field == "turbulence.correlation_length_m") {
    if (!ParseDouble(value, &o->turbulenceCorrelationLengthM)) return invalidValue();
  } else if (field == "turbulence.missing_data") {
    if (!ParseEnum(value, {{"fail", MissingTurbulenceMode::Fail},
                           {"ballistic", MissingTurbulenceMode::Ballistic}},
                   &o->missingTurbulence)) return invalidValue();
  } else if (field == "turbulence.resonance_range") {
    if (!ParseEnum(value, {{"reject", ResonanceRangeMode::Reject},
                           {"power-law-extension", ResonanceRangeMode::PowerLawExtension}},
                   &o->resonanceRange)) return invalidValue();
  } else if (field == "turbulence.self_consistent_3d") {
    if (!ParseBool(value, &o->enableSelfConsistent3DTurbulence)) return invalidValue();
  } else if (field == "transport.cell_crossing_fraction") {
    if (!ParseDouble(value, &o->cellCrossingFraction)) return invalidValue();
  } else if (field == "transport.diffusion_fraction") {
    if (!ParseDouble(value, &o->diffusionFraction)) return invalidValue();
  } else if (field == "transport.focusing_fraction") {
    if (!ParseDouble(value, &o->focusingFraction)) return invalidValue();
  } else if (field == "transport.cooling_fraction") {
    if (!ParseDouble(value, &o->coolingFraction)) return invalidValue();
  } else if (field == "transport.field_variation_fraction") {
    if (!ParseDouble(value, &o->fieldVariationFraction)) return invalidValue();
  } else if (field == "transport.shock_crossing_fraction") {
    if (!ParseDouble(value, &o->shockCrossingFraction)) return invalidValue();
  } else if (field == "transport.minimum_substep_s") {
    if (!ParseDouble(value, &o->minimumTransportSubstepS)) return invalidValue();
  } else if (field == "transport.maximum_substeps") {
    if (!ParseUnsigned64(value, &o->maximumTransportSubsteps)) return invalidValue();
  } else if (field == "transport.pitch_angle_scheme") {
    if (!ParseEnum(value,
        {{"reflecting-milstein", PitchAngleSchemeMode::ReflectingMilstein},
         {"reflecting-euler-maruyama",
          PitchAngleSchemeMode::ReflectingEulerMaruyama}},
        &o->pitchAngleScheme)) return invalidValue();
  } else if (field == "transport.perpendicular_diffusion") {
    if (!ParseBool(value, &o->enablePerpendicularDiffusion)) return invalidValue();
  } else if (field == "transport.drifts") {
    if (!ParseBool(value, &o->enableDrifts)) return invalidValue();
  } else if (field == "shock.authority") {
    if (!ParseEnum(value, {{"none", ShockAuthority::None},
                           {"swcme", ShockAuthority::Swcme}}, &o->shock))
      return invalidValue();
  } else if (field == "shock.active_from_s") {
    if (!ParseDouble(value, &o->shockModel.activeFromS)) return invalidValue();
  } else if (field == "shock.active_until_s") {
    if (!ParseDouble(value, &o->shockModel.activeUntilS)) return invalidValue();
  } else if (field == "shock.initial_radius_m") {
    if (!ParseDouble(value, &o->shockModel.initialRadiusM)) return invalidValue();
  } else if (field == "shock.maximum_radius_m") {
    if (!ParseDouble(value, &o->shockModel.maximumRadiusM)) return invalidValue();
  } else if (field == "shock.speed_m_per_s") {
    if (!ParseDouble(value, &o->shockModel.speedMPerS)) return invalidValue();
  } else if (field == "shock.compression_ratio") {
    if (!ParseDouble(value, &o->shockModel.compressionRatio)) return invalidValue();
  } else if (field == "source.enabled") {
    if (!ParseBool(value, &o->source.enabled)) return invalidValue();
  } else if (field == "source.physical_particle_rate_per_s") {
    if (!ParseDouble(value, &o->source.physicalParticleRatePerS)) return invalidValue();
  } else if (field == "source.injection_efficiency") {
    if (!ParseDouble(value, &o->source.injectionEfficiency)) return invalidValue();
  } else if (field == "source.minimum_energy_j") {
    if (!ParseDouble(value, &o->source.minimumEnergyJ)) return invalidValue();
  } else if (field == "source.maximum_energy_j") {
    if (!ParseDouble(value, &o->source.maximumEnergyJ)) return invalidValue();
  } else if (field == "source.spectral_index") {
    if (!ParseDouble(value, &o->source.spectralIndex)) return invalidValue();
  } else if (field == "source.samples_per_step") {
    if (!ParseUnsigned64(value, &o->source.samplesPerStep)) return invalidValue();
  } else if (field == "species.name") {
    o->species.name = value;
  } else if (field == "species.mass_kg") {
    if (!ParseDouble(value, &o->species.massKg)) return invalidValue();
  } else if (field == "species.charge_c") {
    if (!ParseDouble(value, &o->species.chargeC)) return invalidValue();
  } else if (field == "species.macroparticle_weight") {
    if (!ParseDouble(value, &o->species.macroparticleWeight)) return invalidValue();
  } else if (field == "output.cadence_steps") {
    if (!ParseUnsigned64(value, &o->outputCadenceSteps)) return invalidValue();
  } else if (field == "output.checkpoint_cadence_steps") {
    if (!ParseUnsigned64(value, &o->checkpointCadenceSteps)) return invalidValue();
  } else if (field == "output.directory") {
    o->outputDirectory = value;
  } else if (field == "output.prefix") {
    o->outputPrefix = value;
  } else if (field == "storage.magnetic_gradient") {
    if (!ParseBool(value, &o->storeMagneticGradient)) return invalidValue();
  } else if (field == "storage.velocity_gradient") {
    if (!ParseBool(value, &o->storeVelocityGradient)) return invalidValue();
  } else if (field == "storage.sampling_bytes_per_cell") {
    if (!ParseSize(value, &o->samplingBytesPerCell)) return invalidValue();
  } else if (field == "restart.input_path") {
    // Text files cannot encode an empty right-hand side because that usually
    // indicates a truncated edit.  The explicit token "none" denotes a fresh
    // run and normalizes to the same empty string used by typed construction.
    o->restartInputPath = Lower(value) == "none" ? std::string() : value;
  } else if (field == "restart.output_path") {
    o->restartOutputPath = value;
  } else {
    ObserverOptions* observer = ObserverForSection(
        section, o, observerDefaultsCleared);
    if (observer == nullptr) return Invalid("unknown configuration key '" + field + "'");
    if (key == "position_x_m") {
      if (!ParseDouble(value, &observer->positionM.x)) return invalidValue();
    } else if (key == "position_y_m") {
      if (!ParseDouble(value, &observer->positionM.y)) return invalidValue();
    } else if (key == "position_z_m") {
      if (!ParseDouble(value, &observer->positionM.z)) return invalidValue();
    } else if (key == "follows_trajectory") {
      if (!ParseBool(value, &observer->followsTrajectory)) return invalidValue();
    } else if (key == "cadence_s") {
      if (!ParseDouble(value, &observer->cadenceS)) return invalidValue();
    } else if (key == "energy_bins") {
      if (!ParseUnsigned(value, &observer->energyBins)) return invalidValue();
    } else if (key == "pitch_angle_bins") {
      if (!ParseUnsigned(value, &observer->pitchAngleBins)) return invalidValue();
    } else if (key == "products") {
      observer->products = value;
    } else if (key == "kind") {
      if (!ParseEnum(value,
          {{"fixed-cartesian", ObserverKind::FixedCartesian},
           {"fixed-heliographic", ObserverKind::FixedHeliographic},
           {"moving-cartesian", ObserverKind::MovingCartesian},
           {"spherical-shell", ObserverKind::SphericalShell},
           {"field-connected", ObserverKind::FieldConnected}},
          &observer->kind)) return invalidValue();
    } else if (key == "normalization") {
      if (!ParseEnum(value,
          {{"represented-particles", ObserverNormalization::RepresentedParticles},
           {"differential-intensity", ObserverNormalization::DifferentialIntensity}},
          &observer->normalization)) return invalidValue();
    } else if (key == "velocity_x_m_per_s") {
      if (!ParseDouble(value, &observer->velocityMPerS.x)) return invalidValue();
    } else if (key == "velocity_y_m_per_s") {
      if (!ParseDouble(value, &observer->velocityMPerS.y)) return invalidValue();
    } else if (key == "velocity_z_m_per_s") {
      if (!ParseDouble(value, &observer->velocityMPerS.z)) return invalidValue();
    } else if (key == "collection_radius_m") {
      if (!ParseDouble(value, &observer->collectionRadiusM)) return invalidValue();
    } else if (key == "shell_radius_m") {
      if (!ParseDouble(value, &observer->shellRadiusM)) return invalidValue();
    } else if (key == "minimum_energy_j") {
      if (!ParseDouble(value, &observer->minimumEnergyJ)) return invalidValue();
    } else if (key == "maximum_energy_j") {
      if (!ParseDouble(value, &observer->maximumEnergyJ)) return invalidValue();
    } else if (key == "minimum_mu") {
      if (!ParseDouble(value, &observer->minimumMu)) return invalidValue();
    } else if (key == "maximum_mu") {
      if (!ParseDouble(value, &observer->maximumMu)) return invalidValue();
    } else if (key == "species") {
      if (!ParseSpeciesList(value, &observer->species)) return invalidValue();
    } else {
      return Invalid("unknown configuration key '" + field + "'");
    }
  }
  return Core::Status::OK();
}

}  // namespace

const char* Name(LogVerbosity value) {
  switch (value) {
    case LogVerbosity::Quiet: return "quiet";
    case LogVerbosity::Normal: return "normal";
    case LogVerbosity::Verbose: return "verbose";
  }
  return "unknown";
}

Core::Status ParseStandaloneCommandLine(
    int argc, char* const argv[], StandaloneCommandLine* result) {
  if (result == nullptr) return Invalid("command-line output is null");
  StandaloneCommandLine candidate;
  for (int i = 1; i < argc; ++i) {
    const std::string argument(argv[i]);
    auto requireValue = [&](const char* option, std::string* value) {
      if (i + 1 >= argc) return false;
      *value = argv[++i];
      return !value->empty() && value->rfind("--", 0) != 0 && *value != option;
    };
    if (argument == "--input") {
      if (!requireValue("--input", &candidate.inputPath))
        return Invalid("--input requires a path");
    } else if (argument == "--output-dir") {
      if (!requireValue("--output-dir", &candidate.outputDirectoryOverride))
        return Invalid("--output-dir requires a path");
    } else if (argument == "--restart") {
      if (!requireValue("--restart", &candidate.restartPath))
        return Invalid("--restart requires a path");
    } else if (argument == "--dry-run") {
      candidate.dryRun = true;
    } else if (argument == "--list-tests") {
      candidate.listTests = true;
    } else if (argument == "--all-tests") {
      candidate.allTests = true;
    } else if (argument == "--test") {
      std::string id;
      if (!requireValue("--test", &id)) return Invalid("--test requires an ID");
      candidate.tests.push_back(id);
    } else if (argument == "--log-level") {
      std::string value;
      if (!requireValue("--log-level", &value) ||
          !ParseEnum(value, {{"quiet", LogVerbosity::Quiet},
                             {"normal", LogVerbosity::Normal},
                             {"verbose", LogVerbosity::Verbose}},
                     &candidate.verbosity)) {
        return Invalid("--log-level accepts quiet, normal, or verbose");
      }
    } else {
      return Invalid("unknown srcSEP3D option '" + argument + "'");
    }
  }
  const int selectionModes = static_cast<int>(candidate.listTests) +
      static_cast<int>(candidate.allTests) +
      static_cast<int>(!candidate.tests.empty());
  if (selectionModes > 1)
    return Invalid("--list-tests, --all-tests, and --test are mutually exclusive");
  if (candidate.inputPath.empty() && selectionModes == 0)
    return Invalid("a standalone production run requires --input PATH");
  *result = candidate;
  return Core::Status::OK();
}

Core::Status ParseConfigurationText(
    const std::string& text, RunConfiguration3DOptions* result) {
  if (result == nullptr) return Invalid("configuration output is null");
  RunConfiguration3DOptions candidate;
  std::istringstream input(text);
  std::string section;
  std::string line;
  std::set<std::string> assigned;
  std::set<std::string> sections;
  bool observerDefaultsCleared = false;
  bool schemaSeen = false;
  std::size_t lineNumber = 0;
  while (std::getline(input, line)) {
    ++lineNumber;
    const std::size_t comment = line.find('#');
    if (comment != std::string::npos) line.erase(comment);
    line = Trim(line);
    if (line.empty()) continue;
    if (line.front() == '[' && line.back() == ']') {
      section = Lower(Trim(line.substr(1, line.size() - 2)));
      if (section.empty()) return Invalid("empty section at line " +
                                          std::to_string(lineNumber));
      if (!KnownSection(section))
        return Invalid("unknown configuration section '[" + section + "]'");
      if (!sections.insert(section).second)
        return Invalid("duplicate configuration section '[" + section + "]'");
      continue;
    }
    const std::size_t separator = line.find('=');
    if (separator == std::string::npos)
      return Invalid("expected key=value at line " + std::to_string(lineNumber));
    const std::string key = Lower(Trim(line.substr(0, separator)));
    const std::string value = Trim(line.substr(separator + 1));
    if (key.empty() || value.empty())
      return Invalid("empty key or value at line " + std::to_string(lineNumber));
    const std::string qualified = section.empty() ? key : section + "." + key;
    if (!assigned.insert(qualified).second)
      return Invalid("duplicate configuration key '" + qualified + "'");
    const Core::Status status = ApplyField(
        section, key, value, &candidate, &observerDefaultsCleared);
    if (!status.ok())
      return Invalid("line " + std::to_string(lineNumber) + ": " + status.message);
    if (qualified == "run.schema_version") schemaSeen = true;
  }
  if (!schemaSeen) return Invalid("missing required run.schema_version");
  // C01 distinguishes file input from typed coupled construction.  A coupled
  // host may intentionally rely on C++ defaults while assigning its supplied
  // fields.  A human-authored production file must instead contain every
  // contract group, even when a group merely disables a feature.  This makes
  // omissions visible during parsing rather than after AMPS initialization.
  const char* requiredSections[] = {
      "run", "domain", "mesh", "mesh.solar", "mesh.tube", "memory",
      "background", "background.parker", "turbulence", "transport",
      "shock", "source", "species", "storage", "output", "restart"};
  for (const char* required : requiredSections) {
    if (sections.count(required) == 0) {
      return Invalid("missing required configuration section '[" +
                     std::string(required) + "]'");
    }
  }
  bool observerSectionSeen = false;
  for (const std::string& present : sections) {
    if (present.rfind("observer.", 0) == 0) observerSectionSeen = true;
  }
  if (!observerSectionSeen)
    return Invalid("at least one [observer.ID] section is required");
  // Parsing alone does not expose a half-valid options record.  Invoke the
  // immutable factory as the schema's single normalization/validation gate,
  // then copy its resolved options back for callers that want programmatic
  // construction before the final shared_ptr is installed.
  std::shared_ptr<const RunConfiguration3D> validated;
  const Core::Status status = RunConfiguration3D::Create(candidate, &validated);
  if (!status.ok()) return status;
  *result = validated->options();
  return Core::Status::OK();
}

Core::Status LoadConfigurationFile(
    const std::string& path, RunConfiguration3DOptions* result) {
  if (path.empty()) return Invalid("configuration path is empty");
  std::ifstream input(path.c_str());
  if (!input) return Invalid("cannot open configuration file '" + path + "'");
  std::ostringstream text;
  text << input.rdbuf();
  if (!input.good() && !input.eof())
    return Invalid("failed while reading configuration file '" + path + "'");
  return ParseConfigurationText(text.str(), result);
}

Core::Status BuildStandaloneRunRequest(
    int argc, char* const argv[], StandaloneRunRequest* result) {
  if (result == nullptr) return Invalid("standalone request output is null");
  StandaloneRunRequest candidate;
  Core::Status status = ParseStandaloneCommandLine(
      argc, argv, &candidate.commandLine);
  if (!status.ok()) return status;
  if (candidate.commandLine.listTests || candidate.commandLine.allTests ||
      !candidate.commandLine.tests.empty()) {
    *result = candidate;
    return Core::Status::OK();
  }
  RunConfiguration3DOptions options;
  status = LoadConfigurationFile(candidate.commandLine.inputPath, &options);
  if (!status.ok()) return status;
  if (!candidate.commandLine.outputDirectoryOverride.empty())
    options.outputDirectory = candidate.commandLine.outputDirectoryOverride;
  if (!candidate.commandLine.restartPath.empty())
    options.restartInputPath = candidate.commandLine.restartPath;
  status = RunConfiguration3D::Create(options, &candidate.configuration);
  if (!status.ok()) return status;
  *result = candidate;
  return Core::Status::OK();
}

Core::Status BuildDryRunSummary(const RunConfiguration3D& configuration,
                                std::string* summary) {
  if (summary == nullptr) return Invalid("dry-run summary output is null");
  const RunConfiguration3DOptions& options = configuration.options();
  Mesh::ResolutionConfiguration resolution;
  resolution.innerRadiusM = options.innerRadiusM;
  resolution.outerRadiusM = options.outerRadiusM;
  resolution.minimumCellSizeM = options.minimumCellSizeM;
  resolution.backgroundCellSizeM = options.backgroundCellSizeM;
  resolution.enableRadialRefinement = options.enableRadialRefinement;
  resolution.solarSurfaceCellSizeM = options.solarSurfaceCellSizeM;
  resolution.solarRefinementOuterRadiusM = options.solarRefinementOuterRadiusM;
  resolution.solarRefinementProfile = options.solarRefinementProfile;
  resolution.solarRefinementExponent = options.solarRefinementExponent;
  resolution.enableTubeRefinement = options.enableTubeRefinement;
  resolution.tubeLongitudeRad = options.tubeLongitudeRad;
  resolution.tubeColatitudeRad = options.tubeColatitudeRad;
  resolution.tubeReferenceRadiusM = options.tubeReferenceRadiusM;
  resolution.tubeRadiusAtReferenceM = options.tubeRadiusAtReferenceM;
  resolution.tubeRadiusMode = options.tubeRadiusMode;
  resolution.tubeCellSizeM = options.tubeCellSizeM;
  resolution.tubeTransverseProfile = options.tubeTransverseProfile;
  resolution.tubeTransverseExponent = options.tubeTransverseExponent;
  resolution.solarWindSpeedMPerS = options.parker.solarWindSpeedMPerS;
  resolution.solarRotationRateRadPerS = options.parker.solarRotationRateRadPerS;
  resolution.cellsPerBlockEdge = options.meshCellsPerBlockEdge;
  resolution.maximumLevel = options.maximumMeshLevel;
  resolution.blockOverheadBytes = options.meshBlockOverheadBytes;
  resolution.memoryBudgetBytes = options.meshMemoryBudgetBytes;
  resolution.memoryModel = options.memoryModel;
  const Mesh::DomainBounds domain = Mesh::MakeDomain(options);
  Mesh::RefinementPreflight preflight;
  const Core::Status status = Mesh::BuildRefinementPreflight(
      domain, resolution, configuration.storage_layout(), &preflight);
  if (!status.ok()) return status;
  std::ostringstream output;
  output << std::setprecision(17) << std::scientific
         << "srcSEP3D dry-run configuration\n"
         << "physics_fingerprint=" << configuration.physics_fingerprint() << '\n'
         << "domain_preset=" << Name(options.domain) << '\n'
         << "inner_radius_m=" << options.innerRadiusM << '\n'
         << "outer_radius_m=" << options.outerRadiusM << '\n'
         << "minimum_requested_cell_m=" << preflight.minimumRequestedCellM << '\n'
         << "maximum_requested_cell_m=" << preflight.maximumRequestedCellM << '\n'
         << "tube_radius_at_reference_m=" << preflight.tubeRadiusAtReferenceM << '\n'
         << "estimated_memory_bytes=" << preflight.memory.totalBytes << '\n'
         << "memory_safety_margin_bytes=" << preflight.memory.safetyMarginBytes << '\n'
         << "estimated_blocks_by_level=";
  for (std::size_t i = 0; i < preflight.estimatedBlocksByLevel.size(); ++i) {
    if (i != 0) output << ',';
    output << i << ':' << preflight.estimatedBlocksByLevel[i];
  }
  output << '\n' << "resolved_manifest=" << configuration.resolved_manifest() << '\n';
  *summary = output.str();
  return Core::Status::OK();
}

}  // namespace RuntimeModel
}  // namespace SEP3D
