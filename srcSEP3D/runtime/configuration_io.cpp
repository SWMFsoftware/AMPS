#include "configuration_io.h"

#include "../mesh/mesh_model.h"
#include "swcme3d_input.hpp"

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

bool NearlyEqual(double left, double right) {
  return std::fabs(left - right) <=
      128.0 * std::numeric_limits<double>::epsilon() *
      std::max(std::numeric_limits<double>::min(),
               std::max(std::fabs(left), std::fabs(right)));
}

bool KnownSection(const std::string& section) {
  static const std::set<std::string> fixed = {
      "run", "domain", "parker_spiral", "mesh", "mesh.solar", "mesh.tube", "memory",
      "background", "background.parker", "turbulence", "transport",
      "shock", "source", "species", "storage", "output", "restart",
      "swcme"};
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
    if (!ParseUnsigned64(value, &version) || version < 1 || version > 3)
      return invalidValue();
    o->inputSchemaVersion = static_cast<unsigned>(version);
    return Core::Status::OK();
  }
  if (section == "swcme") {
    SwcmeAssignment assignment;
    assignment.key = key;
    assignment.value = value;
    o->swcmeAssignments.push_back(assignment);
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
  } else if (field == "parker_spiral.origin_x_m") {
    if (!ParseDouble(value, &o->parkerSpiralOriginM.x)) return invalidValue();
  } else if (field == "parker_spiral.origin_y_m") {
    if (!ParseDouble(value, &o->parkerSpiralOriginM.y)) return invalidValue();
  } else if (field == "parker_spiral.origin_z_m") {
    if (!ParseDouble(value, &o->parkerSpiralOriginM.z)) return invalidValue();
  } else if (field == "parker_spiral.start_mode") {
    if (!ParseEnum(value,
        {{"explicit", ParkerSpiralStartMode::Explicit},
         {"cme-launch-point", ParkerSpiralStartMode::CmeLaunchPoint}},
        &o->parkerSpiralStartMode)) return invalidValue();
  } else if (field == "parker_spiral.initial_x_m") {
    if (!ParseDouble(value, &o->parkerSpiralInitialPointM.x)) return invalidValue();
  } else if (field == "parker_spiral.initial_y_m") {
    if (!ParseDouble(value, &o->parkerSpiralInitialPointM.y)) return invalidValue();
  } else if (field == "parker_spiral.initial_z_m") {
    if (!ParseDouble(value, &o->parkerSpiralInitialPointM.z)) return invalidValue();
  } else if (field == "parker_spiral.length_m") {
    if (!ParseDouble(value, &o->parkerSpiralLengthM)) return invalidValue();
  } else if (field == "parker_spiral.point_count") {
    if (!ParseUnsigned64(value, &o->parkerSpiralPointCount)) return invalidValue();
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
         {"python-interpolator", BackgroundAuthority::PythonInterpolator},
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
  } else if (field == "background.parker.number_density_at_one_au_m3") {
    // Complete schema-3 input uses the unambiguous SWCME spelling.  The older
    // "at_reference" key remains parseable for version-1/2 files, where the
    // magnetic and density references historically shared one label.
    if (!ParseDouble(value, &o->parker.numberDensityAtReferenceM3))
      return invalidValue();
    o->parker.densityReferenceRadiusM = Core::Const::AU;
  } else if (field == "background.parker.temperature_k") {
    if (!ParseDouble(value, &o->parker.temperatureK)) return invalidValue();
  } else if (field == "background.parker.validity_cadence_s") {
    if (!ParseDouble(value, &o->parker.validityCadenceS)) return invalidValue();
  } else if (field == "turbulence.authority") {
    if (!ParseEnum(value, {{"prescribed", TurbulenceAuthority::Prescribed},
                           {"swmf", TurbulenceAuthority::Swmf}},
                   &o->turbulence)) return invalidValue();
  } else if (field == "turbulence.model") {
    if (!ParseEnum(value,
        {{"power-law", PrescribedTurbulenceModel::PowerLaw},
         {"kolmogorov", PrescribedTurbulenceModel::Kolmogorov},
         {"kraichnan", PrescribedTurbulenceModel::Kraichnan}},
        &o->prescribedTurbulenceModel)) return invalidValue();
  } else if (field == "turbulence.amplitude_model") {
    if (!ParseEnum(value,
        {{"constant-delta-b-over-b",
              PrescribedTurbulenceAmplitudeModel::ConstantDeltaBOverB},
         {"wave-energy-power-law",
              PrescribedTurbulenceAmplitudeModel::WaveEnergyPowerLaw}},
        &o->prescribedTurbulenceAmplitudeModel)) return invalidValue();
  } else if (field == "turbulence.delta_b_over_b") {
    if (!ParseDouble(value, &o->prescribedDeltaBOverB)) return invalidValue();
  } else if (field ==
             "turbulence.wave_energy_density_at_reference_j_per_m3") {
    if (!ParseDouble(value, &o->turbulenceWaveEnergyAtReferenceJPerM3))
      return invalidValue();
  } else if (field == "turbulence.wave_energy_density_radial_exponent") {
    if (!ParseDouble(value, &o->turbulenceWaveEnergyRadialExponent))
      return invalidValue();
  } else if (field == "turbulence.normalized_cross_helicity") {
    if (!ParseDouble(value, &o->turbulenceNormalizedCrossHelicity))
      return invalidValue();
  } else if (field == "turbulence.reference_radius_m") {
    if (!ParseDouble(value, &o->turbulenceReferenceRadiusM))
      return invalidValue();
  } else if (field == "turbulence.k_min_per_m") {
    if (!ParseDouble(value, &o->turbulenceKMinPerM)) return invalidValue();
  } else if (field == "turbulence.k_max_per_m") {
    if (!ParseDouble(value, &o->turbulenceKMaxPerM)) return invalidValue();
  } else if (field == "turbulence.k_min_radial_exponent") {
    if (!ParseDouble(value, &o->turbulenceKMinRadialExponent))
      return invalidValue();
  } else if (field == "turbulence.k_max_radial_exponent") {
    if (!ParseDouble(value, &o->turbulenceKMaxRadialExponent))
      return invalidValue();
  } else if (field == "turbulence.spectral_index") {
    if (!ParseDouble(value, &o->turbulenceSpectralIndex)) return invalidValue();
  } else if (field == "turbulence.correlation_length_m") {
    if (!ParseDouble(value, &o->turbulenceCorrelationLengthM)) return invalidValue();
  } else if (field == "turbulence.correlation_length_radial_exponent") {
    if (!ParseDouble(value, &o->turbulenceCorrelationLengthRadialExponent))
      return invalidValue();
  } else if (field == "turbulence.validity_cadence_s") {
    if (!ParseDouble(value, &o->turbulenceValidityCadenceS))
      return invalidValue();
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
    if (!ParseEnum(value,
        {{"none", PerpendicularDiffusionMode::None},
         {"constant", PerpendicularDiffusionMode::Constant},
         {"constant-ratio", PerpendicularDiffusionMode::ConstantRatio}},
        &o->perpendicularDiffusion)) {
      // Input files from the parallel-only release used "false". Preserve
      // that inert spelling, but reject "true" because it never identified a
      // coefficient closure and therefore cannot be migrated unambiguously.
      bool legacy = false;
      if (!ParseBool(value, &legacy) || legacy) return invalidValue();
      o->perpendicularDiffusion = PerpendicularDiffusionMode::None;
    }
  } else if (field == "transport.constant_kappa_perpendicular_m2_per_s") {
    if (!ParseDouble(value, &o->constantKappaPerpendicularM2PerS))
      return invalidValue();
  } else if (field == "transport.kappa_perpendicular_to_parallel_ratio") {
    if (!ParseDouble(value, &o->kappaPerpendicularToParallelRatio))
      return invalidValue();
  } else if (field == "transport.drifts") {
    if (!ParseEnum(value,
        {{"none", DriftMode::None}, {"gradient-b", DriftMode::GradientB},
         {"curvature", DriftMode::Curvature},
         {"gradient-curvature", DriftMode::GradientAndCurvature}},
        &o->drift)) {
      bool legacy = false;
      if (!ParseBool(value, &legacy) || legacy) return invalidValue();
      o->drift = DriftMode::None;
    }
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
  } else if (field == "species.macroparticle_weight") {
    // This is deliberately the only post-compile species field.  Count,
    // symbols, masses, charges, and indices are immutable products of the AMPS
    // SpeciesList deck and are enumerated through the AMPS molecular-data API.
    if (!ParseDouble(value, &o->species.macroparticleWeight)) return invalidValue();
  } else if (field == "output.cadence_steps") {
    if (!ParseUnsigned64(value, &o->outputCadenceSteps)) return invalidValue();
  } else if (field == "output.checkpoint_cadence_steps") {
    if (!ParseUnsigned64(value, &o->checkpointCadenceSteps)) return invalidValue();
  } else if (field == "output.directory") {
    o->outputDirectory = value;
  } else if (field == "output.prefix") {
    o->outputPrefix = value;
  } else if (field == "output.initialization_mesh_tecplot_file") {
    o->initializationMeshTecplotFile = value;
  } else if (field == "output.initialization_parker_line_tecplot_file") {
    o->initializationParkerLineTecplotFile = value;
  } else if (field == "output.initialization_data_tecplot_file") {
    o->initializationDataTecplotFile = value;
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
    } else if (key == "energy_spacing") {
      if (!ParseEnum(value,
          {{"logarithmic", EnergyChannelSpacing::Logarithmic},
           {"linear", EnergyChannelSpacing::Linear}},
          &observer->energyChannelSpacing)) return invalidValue();
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
      // Species identities do not exist in the post-compile input schema.
      // `all` therefore means every entry in the immutable AMPS SpeciesList,
      // whatever count that particular executable was generated with.  An
      // explicit comma-separated list remains available for a deliberate
      // observer subset and is upper-bound checked after AMPS exposes its
      // generated table.
      if (Lower(value) == "all") {
        observer->allCompiledSpecies = true;
        observer->species.clear();
      } else {
        if (!ParseSpeciesList(value, &observer->species)) return invalidValue();
        observer->allCompiledSpecies = false;
      }
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
    } else if (argument == "--initialization-only") {
      candidate.initializationOnly = true;
    } else if (argument == "--initialization-output-dir") {
      if (!requireValue("--initialization-output-dir",
                        &candidate.initializationOutputDirectory)) {
        return Invalid("--initialization-output-dir requires a path");
      }
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
  if (candidate.initializationOnly && candidate.dryRun) {
    return Invalid("--initialization-only and --dry-run are mutually exclusive: "
                   "the former builds the AMPS mesh, while the latter forbids "
                   "AMPS allocation");
  }
  if (!candidate.initializationOutputDirectory.empty() &&
      !candidate.initializationOnly) {
    return Invalid("--initialization-output-dir requires --initialization-only");
  }
  if (candidate.initializationOnly && selectionModes != 0) {
    return Invalid("--initialization-only cannot be combined with test selection");
  }
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
    // A bare assignment used to fall through to ApplyField(), which reported
    // only "unknown configuration key".  That message made a legacy flat deck
    // look like a misspelled schema-3 field.  The production grammar is strict
    // INI, so identify the actual structural error and point to the canonical
    // first section without guessing which physics group owned the old key.
    if (section.empty()) {
      return Invalid("line " + std::to_string(lineNumber) + ": key '" + key +
                     "' appears before any [section]; srcSEP3D input uses the "
                     "strict INI schema and must begin with [run] (see "
                     "srcSEP3D/examples/sep3d_analytic_parker.in)");
    }
    const std::string qualified = section.empty() ? key : section + "." + key;
    if (!assigned.insert(qualified).second)
      return Invalid("duplicate configuration key '" + qualified + "'");
    const Core::Status status = ApplyField(
        section, key, value, &candidate, &observerDefaultsCleared);
    if (!status.ok())
      return Invalid("line " + std::to_string(lineNumber) + ": " + status.message);
    if (qualified == "run.schema_version") schemaSeen = true;
    if (section == "swcme" && !candidate.swcmeAssignments.empty())
      candidate.swcmeAssignments.back().line = lineNumber;
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
  if (candidate.inputSchemaVersion >= 2) {
    if (sections.count("parker_spiral") == 0)
      return Invalid("schema version 2 requires [parker_spiral]");
    const char* requiredLineFields[] = {
        "parker_spiral.origin_x_m", "parker_spiral.origin_y_m",
        "parker_spiral.origin_z_m", "parker_spiral.initial_x_m",
        "parker_spiral.initial_y_m", "parker_spiral.initial_z_m",
        "parker_spiral.length_m", "parker_spiral.point_count"};
    for (const char* required : requiredLineFields) {
      if (assigned.count(required) == 0)
        return Invalid("schema version 2 is missing required key '" +
                       std::string(required) + "'");
    }
  }
  if (candidate.inputSchemaVersion >= 3) {
    if (sections.count("swcme") == 0)
      return Invalid("schema version 3 requires a complete [swcme] section");
    // Schema 3 is the no-assumptions production surface.  Every application-
    // owned value is present even if its selected mode makes it inactive; this
    // prevents a later mode edit from reviving a C++ default that never
    // appeared in the reviewed input deck.
    const char* requiredFields[] = {
        "run.schema_version", "run.intent", "run.transport",
        "run.time_step_s", "run.maximum_time_steps", "run.campaign_seed",
        "run.background_cadence_steps", "run.injection_cadence_steps",
        "domain.preset", "domain.inner_radius_m", "domain.inner_boundary",
        "domain.outer_radius_mode", "domain.outer_radius_m",
        "domain.outer_boundary", "domain.coordinate_frame",
        "domain.origin_x_m", "domain.origin_y_m", "domain.origin_z_m",
        "parker_spiral.origin_x_m", "parker_spiral.origin_y_m",
        "parker_spiral.origin_z_m", "parker_spiral.initial_x_m",
        "parker_spiral.initial_y_m", "parker_spiral.initial_z_m",
        "parker_spiral.start_mode", "parker_spiral.length_m",
        "parker_spiral.point_count",
        "mesh.global_cell_size_m", "mesh.minimum_cell_size_m",
        "mesh.cells_per_block_edge", "mesh.maximum_level",
        "mesh.memory_budget_bytes", "mesh.block_overhead_bytes",
        "mesh.solar.enabled", "mesh.solar.surface_cell_size_m",
        "mesh.solar.transition_outer_radius_m", "mesh.solar.profile",
        "mesh.solar.exponent", "mesh.tube.enabled",
        "mesh.tube.source_longitude_rad",
        "mesh.tube.source_colatitude_rad",
        "mesh.tube.reference_radius_m", "mesh.tube.radius_at_reference_m",
        "mesh.tube.radius_mode", "mesh.tube.center_cell_size_m",
        "mesh.tube.transverse_profile", "mesh.tube.transverse_exponent",
        "memory.base_cell_bytes", "memory.base_node_bytes",
        "memory.block_structure_bytes",
        "memory.communication_bytes_per_block", "memory.particle_bytes",
        "memory.particles_per_cell", "memory.halo_fraction",
        "memory.safety_margin_fraction", "background.provider",
        "background.external_script", "background.parker.reference_radius_m",
        "background.parker.radial_field_at_reference_t",
        "background.parker.solar_rotation_rate_rad_per_s",
        "background.parker.solar_wind_speed_m_per_s",
        "background.parker.magnetic_polarity",
        "background.parker.number_density_at_one_au_m3",
        "background.parker.temperature_k",
        "background.parker.validity_cadence_s", "turbulence.authority",
        "turbulence.model", "turbulence.amplitude_model",
        "turbulence.delta_b_over_b",
        "turbulence.wave_energy_density_at_reference_j_per_m3",
        "turbulence.wave_energy_density_radial_exponent",
        "turbulence.normalized_cross_helicity",
        "turbulence.reference_radius_m", "turbulence.k_min_per_m",
        "turbulence.k_max_per_m", "turbulence.k_min_radial_exponent",
        "turbulence.k_max_radial_exponent", "turbulence.spectral_index",
        "turbulence.correlation_length_m",
        "turbulence.correlation_length_radial_exponent",
        "turbulence.validity_cadence_s", "turbulence.missing_data",
        "turbulence.resonance_range", "turbulence.self_consistent_3d",
        "transport.cell_crossing_fraction", "transport.diffusion_fraction",
        "transport.focusing_fraction", "transport.cooling_fraction",
        "transport.field_variation_fraction",
        "transport.shock_crossing_fraction", "transport.minimum_substep_s",
        "transport.maximum_substeps", "transport.pitch_angle_scheme",
        "transport.perpendicular_diffusion",
        "transport.constant_kappa_perpendicular_m2_per_s",
        "transport.kappa_perpendicular_to_parallel_ratio",
        "transport.drifts", "shock.authority", "source.enabled",
        "source.physical_particle_rate_per_s",
        "source.injection_efficiency", "source.minimum_energy_j",
        "source.maximum_energy_j", "source.samples_per_step",
        "species.macroparticle_weight",
        "storage.magnetic_gradient", "storage.velocity_gradient",
        "storage.sampling_bytes_per_cell", "output.cadence_steps",
        "output.checkpoint_cadence_steps", "output.directory",
        "output.prefix", "output.initialization_mesh_tecplot_file",
        "output.initialization_parker_line_tecplot_file",
        "output.initialization_data_tecplot_file",
        "restart.input_path", "restart.output_path"};
    for (const char* required : requiredFields)
      if (assigned.count(required) == 0)
        return Invalid("schema version 3 is missing required key '" +
                       std::string(required) + "'");
    const char* observerFields[] = {
        "kind", "normalization", "position_x_m", "position_y_m",
        "position_z_m", "follows_trajectory", "velocity_x_m_per_s",
        "velocity_y_m_per_s", "velocity_z_m_per_s",
        "collection_radius_m", "shell_radius_m", "cadence_s",
        "energy_bins", "energy_spacing", "pitch_angle_bins", "minimum_energy_j",
        "maximum_energy_j", "minimum_mu", "maximum_mu", "species",
        "products"};
    for (const std::string& presentSection : sections) {
      if (presentSection.rfind("observer.", 0) != 0) continue;
      for (const char* suffix : observerFields) {
        const std::string required = presentSection + "." + suffix;
        if (assigned.count(required) == 0)
          return Invalid("schema version 3 is missing required key '" +
                         required + "'");
      }
    }
    if (candidate.injectionCadenceSteps != 1)
      return Invalid("schema version 3 defines source.samples_per_step as "
                     "an exact count at every simulation step; therefore "
                     "run.injection_cadence_steps must equal 1");
    if (candidate.background == BackgroundAuthority::PythonInterpolator)
      return Core::Status::Reserved(
          "Python heliospheric-model interpolation background");
    if (candidate.background != BackgroundAuthority::AnalyticParker ||
        candidate.turbulence != TurbulenceAuthority::Prescribed)
      return Invalid("standalone schema version 3 requires analytic-parker "
                     "background and prescribed turbulence; coupled SWMF "
                     "hosts must use the parser-free typed interface");

    std::vector<swcme::input3d::Assignment> assignments;
    assignments.reserve(candidate.swcmeAssignments.size());
    for (const SwcmeAssignment& raw : candidate.swcmeAssignments) {
      swcme::input3d::Assignment assignment;
      assignment.key = raw.key;
      assignment.value = raw.value;
      assignment.origin = "srcSEP3D input";
      assignment.line = raw.line;
      assignments.push_back(assignment);
    }
    const swcme::input3d::ResolveResult resolved =
        swcme::input3d::Resolve(assignments);
    if (!resolved.ok()) {
      std::ostringstream message;
      message << "invalid complete SWCME3D configuration key='"
              << resolved.status.key << "'";
      if (resolved.status.line != 0)
        message << " line=" << resolved.status.line;
      if (!resolved.status.message.empty())
        message << ": " << resolved.status.message;
      return Invalid(message.str());
    }
    const swcme3d::Params& model = resolved.configuration.model;
    const swcme::sep::SpectrumConfig& spectrum =
        resolved.configuration.spectrum;
    if (candidate.intent != RunIntent::ShockInjection ||
        candidate.shock != ShockAuthority::Swcme ||
        !candidate.source.enabled)
      return Invalid("schema version 3 requires shock-injection intent, "
                     "shock.authority=swcme, and source.enabled=true");
    if (model.shape != swcme3d::ShockShape::Sphere ||
        model.region_mode != swcme::regions::Mode::ShockOnly ||
        model.shock_acceleration_mode != swcme::acceleration::Mode::Source)
      return Invalid("srcSEP3D currently supports canonical SWCME only as a "
                     "spherical SHOCK_ONLY/SOURCE provider; other geometries "
                     "would require a non-spherical AMPS crossing operator");
    const char* retiredShockFields[] = {
        "shock.active_from_s", "shock.active_until_s",
        "shock.initial_radius_m", "shock.maximum_radius_m",
        "shock.speed_m_per_s", "shock.compression_ratio"};
    for (const char* field : retiredShockFields)
      if (assigned.count(field) != 0)
        return Invalid(std::string("schema version 3 rejects duplicate legacy '") +
                       field + "'; [swcme] is authoritative");
    if (assigned.count("source.spectral_index") != 0)
      return Invalid("schema version 3 rejects source.spectral_index; "
                     "canonical local SWCME compression owns the DSA index");
    if (assigned.count(
            "background.parker.number_density_at_reference_m3") != 0) {
      return Invalid(
          "schema version 3 rejects the ambiguous legacy key "
          "background.parker.number_density_at_reference_m3; use "
          "number_density_at_one_au_m3 because SWCME density is normalized "
          "at one AU independently of the magnetic reference radius");
    }

    // `cme-launch-point` is an explicit cross-model constraint, not a loose
    // suggestion.  The canonical SWCME resolver owns both the launch radius
    // and the apex direction; normalize the latter exactly as SWCME does and
    // require the reviewed Parker and mesh values to identify that same point.
    // The current finite Parker line begins on the AMPS inner sphere, so a
    // linked launch radius outside that sphere is rejected rather than moving
    // the mesh boundary or shortening the line implicitly.
    if (candidate.parkerSpiralStartMode ==
        ParkerSpiralStartMode::CmeLaunchPoint) {
      const double directionNorm = std::sqrt(
          model.cme_dir[0] * model.cme_dir[0] +
          model.cme_dir[1] * model.cme_dir[1] +
          model.cme_dir[2] * model.cme_dir[2]);
      if (!std::isfinite(directionNorm) || directionNorm <= 0.0)
        return Invalid("canonical [swcme] CME direction cannot define its "
                       "launch-apex point");
      const Core::Vec3 direction(
          model.cme_dir[0] / directionNorm,
          model.cme_dir[1] / directionNorm,
          model.cme_dir[2] / directionNorm);
      const double launchRadiusM = model.r0_Rs *
          swcme::constants::SOLAR_RADIUS_M;
      const Core::Vec3 launchPoint =
          candidate.coordinateOriginM + launchRadiusM * direction;
      const double linkageTolerance = 1.0e-10 *
          std::max(1.0, std::max(candidate.innerRadiusM, launchRadiusM));
      if (std::fabs(launchRadiusM - candidate.innerRadiusM) >
          linkageTolerance) {
        return Invalid("parker_spiral.start_mode=cme-launch-point requires "
                       "SWCME cme.launch_radius to equal domain.inner_radius_m");
      }
      if ((candidate.parkerSpiralInitialPointM - launchPoint).Norm() >
          linkageTolerance) {
        return Invalid("[parker_spiral] initial point differs from the "
                       "canonical SWCME launch-apex point");
      }
      const double sine = std::sin(candidate.tubeColatitudeRad);
      const Core::Vec3 tubeDirection(
          sine * std::cos(candidate.tubeLongitudeRad),
          sine * std::sin(candidate.tubeLongitudeRad),
          std::cos(candidate.tubeColatitudeRad));
      if ((tubeDirection - direction).Norm() > 1.0e-12) {
        return Invalid("[mesh.tube] source direction differs from the "
                       "canonical SWCME CME launch-apex direction");
      }
      candidate.cmeLaunchPointResolved = true;
      candidate.cmeLaunchPointM = launchPoint;
      // Replace agreeing decimal spellings with the canonical derived point
      // so fingerprints and geometry use one exact binary authority.
      candidate.parkerSpiralInitialPointM = launchPoint;
    }

    // The AMR Parker field, SWCME upstream state, and source energy interval
    // must describe the same physical system.  Species identity is purposely
    // not compared here: the generated AMPS table is the sole authority and
    // injection converts these declared kinetic-energy bounds to momentum
    // independently with every compiled species' AMPS mass.
    const double parkerSourceM = model.parker_source_radius_Rs *
        swcme::constants::SOLAR_RADIUS_M;
    // SWCME's B1AU_nT is the *total* Parker magnitude at one AU and at the
    // configured reference colatitude.  The analytic background instead
    // accepts the radial component at an arbitrary declared reference radius.
    // First remove the one-AU azimuthal winding, then apply magnetic-flux
    // conservation Br proportional to r^-2.  Using the configured reference
    // radius inside the magnitude conversion would be physically wrong for
    // every otherwise-valid deck whose reference is not exactly one AU.
    const double oneAuWinding = model.solar_rotation_rate_rad_s *
        (Core::Const::AU - parkerSourceM) /
        (model.V_sw_kms * 1.0e3) * model.sin_theta;
    const double radialFieldAtOneAuT = model.B1AU_nT * 1.0e-9 /
        std::sqrt(1.0 + oneAuWinding * oneAuWinding);
    const double radialScale =
        Core::Const::AU / candidate.parker.referenceRadiusM;
    const double canonicalRadialFieldT =
        radialFieldAtOneAuT * radialScale * radialScale;
    const bool sameBackground =
        NearlyEqual(model.V_sw_kms * 1.0e3,
                    candidate.parker.solarWindSpeedMPerS) &&
        NearlyEqual(model.solar_rotation_rate_rad_s,
                    candidate.parker.solarRotationRateRadPerS) &&
        // RunConfiguration3D::Create derives parker.sourceRadiusM from the
        // explicit domain inner radius.  Compare against that input here as
        // well; comparing with the pre-factory C++ default would silently
        // restrict complete schema-3 files to 20 solar radii.
        NearlyEqual(parkerSourceM, candidate.innerRadiusM) &&
        NearlyEqual(model.n1AU_cm3 * 1.0e6,
                    candidate.parker.numberDensityAtReferenceM3) &&
        NearlyEqual(canonicalRadialFieldT,
                    candidate.parker.radialFieldAtReferenceT) &&
        NearlyEqual(model.T_K, candidate.parker.temperatureK) &&
        model.parker_radial_polarity == candidate.parker.magneticPolarity;
    const bool sameSource =
        NearlyEqual(spectrum.kinetic_energy_min_MeV * 1.0e6 * Core::Const::e,
                    candidate.source.minimumEnergyJ) &&
        NearlyEqual(spectrum.kinetic_energy_max_MeV * 1.0e6 * Core::Const::e,
                    candidate.source.maximumEnergyJ) &&
        NearlyEqual(resolved.configuration.injection_efficiency,
                    candidate.source.injectionEfficiency) &&
        spectrum.normalization == swcme::sep::NormalizationMode::RelativeOnly;
    const bool sameParkerGeometry =
        NearlyEqual(model.sin_theta,
                    std::sin(candidate.tubeColatitudeRad)) &&
        NearlyEqual(model.solar_rotation_axis[0], 0.0) &&
        NearlyEqual(model.solar_rotation_axis[1], 0.0) &&
        NearlyEqual(model.solar_rotation_axis[2], 1.0);
    if (!sameBackground)
      return Invalid("[background.parker] differs from the canonical "
                     "[swcme] wind, Parker field, density, or temperature");
    if (!sameSource)
      return Invalid("[source] energy/efficiency differs from canonical "
                     "[swcme], or SWCME normalization is not relative_only");
    if (!sameParkerGeometry)
      return Invalid("[mesh.tube] colatitude or analytic +Z rotation axis "
                     "differs from canonical [swcme] Parker geometry");

    // From this point onward the typed Parker provider receives values copied
    // from the canonical SWCME resolver, not a second independently parsed
    // solar-wind model.  The [background.parker] values above are retained as
    // fail-closed cross-checks for human review; after agreement they are
    // normalized to the exact SWCME binary values so spelling/rounding cannot
    // create two ambient states with one configuration fingerprint.
    candidate.parker.radialFieldAtReferenceT = canonicalRadialFieldT;
    candidate.parker.solarWindSpeedMPerS = model.V_sw_kms * 1.0e3;
    candidate.parker.solarRotationRateRadPerS =
        model.solar_rotation_rate_rad_s;
    candidate.parker.magneticPolarity = model.parker_radial_polarity;
    candidate.parker.numberDensityAtReferenceM3 = model.n1AU_cm3 * 1.0e6;
    candidate.parker.densityReferenceRadiusM = Core::Const::AU;
    candidate.parker.temperatureK = model.T_K;
    candidate.parker.adiabaticIndex = model.gamma_ad;
    candidate.parker.thermodynamicClosure =
        model.thermodynamic_closure ==
                swcme::solarwind::ThermodynamicClosure::MultiSpecies
            ? SolarWindThermodynamicClosure::MultiSpecies
            : SolarWindThermodynamicClosure::ProtonOnly;
    candidate.parker.alphaToProtonRatio = model.alpha_to_proton_ratio;
    candidate.parker.electronTemperatureK = model.electron_T_K;
    candidate.parker.alphaTemperatureK = model.alpha_T_K;
    candidate.parker.referenceSinColatitude = model.sin_theta;
    candidate.parker.rotationAxis = {
        model.solar_rotation_axis[0], model.solar_rotation_axis[1],
        model.solar_rotation_axis[2]};

    const long double representedPerEvent =
        static_cast<long double>(candidate.source.physicalParticleRatePerS) *
        static_cast<long double>(candidate.source.injectionEfficiency) *
        static_cast<long double>(model.relative_source_weight_per_area) *
        static_cast<long double>(candidate.requestedTimeStepS) *
        static_cast<long double>(candidate.injectionCadenceSteps);
    const double derivedWeight = static_cast<double>(
        representedPerEvent /
        static_cast<long double>(candidate.source.samplesPerStep));
    if (!std::isfinite(derivedWeight) || derivedWeight <= 0.0 ||
        !NearlyEqual(derivedWeight, candidate.species.macroparticleWeight))
      return Invalid("species.macroparticle_weight must equal "
                     "rate*efficiency*relative_source_weight*cadence_dt/"
                     "source.samples_per_step");
    candidate.swcmeConfigurationFingerprint =
        resolved.configuration.fingerprint;
    candidate.swcmeResolvedManifest =
        resolved.configuration.normalized_manifest;
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

Core::Status ApplyInitializationOutputDirectory(
    const std::string& directory, RunConfiguration3DOptions* options) {
  if (options == nullptr) return Invalid("initialization options output is null");
  if (directory.empty())
    return Invalid("initialization output directory is empty");

  // Retain all three reviewed leaf names from [output]. Only their parent directory
  // is a command-line concern.  Treat both slash spellings as separators so a
  // deck copied between systems does not accidentally embed its former parent
  // beneath the requested preview directory.
  const auto leafName = [](const std::string& path) {
    const std::size_t separator = path.find_last_of("/\\");
    return separator == std::string::npos ? path : path.substr(separator + 1);
  };
  const std::string meshLeaf = leafName(options->initializationMeshTecplotFile);
  const std::string lineLeaf =
      leafName(options->initializationParkerLineTecplotFile);
  const std::string dataLeaf =
      leafName(options->initializationDataTecplotFile);
  if (meshLeaf.empty() || lineLeaf.empty() || dataLeaf.empty() ||
      meshLeaf == "." || meshLeaf == ".." || lineLeaf == "." ||
      lineLeaf == ".." || dataLeaf == "." || dataLeaf == "..") {
    return Invalid("initialization Tecplot paths must end in file names before "
                   "--initialization-output-dir can be applied");
  }
  const std::string separator =
      (!directory.empty() && directory.back() == '/') ? "" : "/";
  options->initializationMeshTecplotFile = directory + separator + meshLeaf;
  options->initializationParkerLineTecplotFile = directory + separator + lineLeaf;
  options->initializationDataTecplotFile = directory + separator + dataLeaf;
  return Core::Status::OK();
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
  if (!candidate.commandLine.initializationOutputDirectory.empty()) {
    status = ApplyInitializationOutputDirectory(
        candidate.commandLine.initializationOutputDirectory, &options);
    if (!status.ok()) return status;
  }
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
  resolution.originM = options.coordinateOriginM;
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
  resolution.parkerInitialPointM = options.parkerSpiralInitialPointM;
  resolution.parkerLengthM = options.parkerSpiralLengthM;
  resolution.parkerPointCount = options.parkerSpiralPointCount;
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
  std::vector<Core::Vec3> centreline;
  const Core::Status lineStatus =
      Mesh::BuildParkerCenterline(resolution, &centreline);
  if (!lineStatus.ok()) return lineStatus;
  const Core::Vec3& lineEnd = centreline.back();
  std::ostringstream output;
  output << std::setprecision(17) << std::scientific
         << "srcSEP3D dry-run configuration\n"
         << "physics_fingerprint=" << configuration.physics_fingerprint() << '\n'
         << "domain_preset=" << Name(options.domain) << '\n'
         << "inner_radius_m=" << options.innerRadiusM << '\n'
         << "outer_radius_m=" << options.outerRadiusM << '\n'
         << "background_provider=" << Name(options.background) << '\n'
         << "solar_wind_model=canonical-swcme-parker-leblanc\n"
         << "solar_wind_thermodynamic_closure="
         << Name(options.parker.thermodynamicClosure) << '\n'
         << "prescribed_turbulence_model="
         << Name(options.prescribedTurbulenceModel) << '\n'
         << "prescribed_turbulence_amplitude_model="
         << Name(options.prescribedTurbulenceAmplitudeModel) << '\n'
         << "turbulence_wave_energy_reference_j_per_m3="
         << options.turbulenceWaveEnergyAtReferenceJPerM3 << '\n'
         << "turbulence_normalized_cross_helicity="
         << options.turbulenceNormalizedCrossHelicity << '\n'
         << "parker_spiral_start_mode="
         << Name(options.parkerSpiralStartMode) << '\n'
         << "parker_spiral_point_count=" << centreline.size() << '\n'
         << "parker_spiral_length_m=" << options.parkerSpiralLengthM << '\n'
         << "parker_spiral_end_m=" << lineEnd.x << ',' << lineEnd.y << ','
         << lineEnd.z << '\n'
         << "compiled_species_authority=AMPS-SpeciesList\n"
         << "time_step_s=" << options.requestedTimeStepS << '\n'
         << "base_particle_weight=" << options.species.macroparticleWeight << '\n'
         << "observer_count=" << options.observers.size() << '\n'
         << "initialization_data_tecplot_base="
         << options.initializationDataTecplotFile << '\n'
         << "source_samples_per_compiled_species="
         << options.source.samplesPerStep << '\n'
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
