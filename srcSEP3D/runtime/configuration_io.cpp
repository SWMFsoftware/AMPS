#include "configuration_io.h"

#include "../mesh/mesh_model.h"
#include "swcme3d_input.hpp"
#include "background_factory.h"
#include "parallel_diffusion/parallel_diffusion.h"
#include "amps_input.h"
#include "sep_input_layout.h"

#include <algorithm>
#include <cctype>
#include <cerrno>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <filesystem>
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

// Apply one srcSEP3D key.  'spectrum' collects the source: shock reference-
// spectrum settings (decision D4) that are not srcSEP3D options themselves
// but become the canonical SWCME 'source.*' assignments once the whole
// document is known (SEP::InputLayout::SynthesizeSwcmeSource()).
Core::Status ApplyField(const std::string& section, const std::string& key,
                        const std::string& value,
                        RunConfiguration3DOptions* o,
                        bool* observerDefaultsCleared,
                        SEP::InputLayout::SwcmeSourceSettings* spectrum) {
  const std::string field = section.empty() ? key : section + "." + key;
  auto invalidValue = [&]() {
    return Invalid("invalid value for '" + field + "': " + value);
  };

  if (field == "application.sep.schema_version") {
    // Decision D1 (sep_input_layout.h): srcSEP and srcSEP3D share content
    // schema 6.  The srcSEP3D INI-era schemas 1-5 are retired without a
    // compatibility path; naming one is a typed error that points at the
    // converter rather than a silent reinterpretation.  The typed
    // RunConfiguration3DOptions::inputSchemaVersion keeps its contract levels
    // for parser-free hosts; a file always installs level 6, which satisfies
    // every "level >= N" completeness rule of the typed factory.
    std::uint64_t version = 0;
    if (!ParseUnsigned64(value, &version)) return invalidValue();
    if (version != SEP::InputLayout::kSchemaVersion)
      return Invalid("application.sep.schema_version = " + value +
                     " is retired; srcSEP and srcSEP3D read only the shared "
                     "schema 6 (convert with tools/convert_amps_input.py "
                     "sep6-from-srcsep3d)");
    o->inputSchemaVersion = static_cast<unsigned>(version);
    return Core::Status::OK();
  }
  // The model blocks physics: swcme, physics: parallel_diffusion and
  // physics: shock_front never reach ApplyField(): WalkBlocks() hands them to
  // the shared libraries that own their vocabularies.
  if (field == "application.sep.intent") {
    if (!ParseEnum(value, {{"transport-only", RunIntent::TransportOnly},
                           {"shock-injection", RunIntent::ShockInjection},
                           {"shock-propagation", RunIntent::ShockPropagation}},
                   &o->intent)) return invalidValue();
  } else if (field == "application.sep.transport") {
    // One shared key selects the transport formulation in both applications
    // (parker | focused-diffusion | focused-scattering).  The broad mover
    // family follows from it uniquely -- Parker for parker, the focused-
    // transport family for both focused formulations -- so the former
    // schema-5 key application.sep3d.particle_mover is retired and the
    // family is recorded here as an explicit (not defaulted) selection.
    if (!ParseEnum(value, {{"parker", TransportModel::Parker3D},
                           {"focused-diffusion",
                            TransportModel::FocusedDiffusion3D},
                           {"focused-scattering",
                            TransportModel::FocusedScattering3D}},
                   &o->transport)) return invalidValue();
    o->particleMover = o->transport == TransportModel::Parker3D
        ? ParticleMoverFamily::Parker : ParticleMoverFamily::FocusedTransport;
    o->particleMoverExplicitlySelected = true;
  } else if (field == "run.time_step_s") {
    if (!ParseDouble(value, &o->requestedTimeStepS)) return invalidValue();
  } else if (field == "run.time_step_mode") {
    // fixed: run.time_step_s is the global step.  derived-from-mesh: the step
    // is dt = time_step_margin_factor * h_min / maximum_particle_speed_m_s,
    // computed after AMPS allocates the mesh (h_min = smallest allocated
    // characteristic cell length).  Requires the incident-flux normalization.
    if (Lower(value) == "fixed") o->particleNumerics.deriveFromMeshAndShock = false;
    else if (Lower(value) == "derived-from-mesh")
      o->particleNumerics.deriveFromMeshAndShock = true;
    else return invalidValue();
  } else if (field == "run.maximum_particle_speed_m_s") {
    // v_max of the derived step [m/s], in (0, c].
    if (!ParseDouble(value, &o->particleNumerics.maximumParticleSpeedMPerS))
      return invalidValue();
  } else if (field == "run.time_step_margin_factor") {
    // Dimensionless CFL margin f in (0, 1] of the derived step.
    if (!ParseDouble(value, &o->particleNumerics.timeStepMarginFactor))
      return invalidValue();
  } else if (field == "run.maximum_time_steps") {
    if (!ParseUnsigned64(value, &o->maximumTimeSteps)) return invalidValue();
  } else if (field == "run.stop_shock_radius_m") {
    if (!ParseDouble(value, &o->stopShockRadiusM)) return invalidValue();
  } else if (field == "run.campaign_seed") {
    if (!ParseUnsigned64(value, &o->campaignSeed)) return invalidValue();
  } else if (field == "application.sep.background_cadence_steps") {
    if (!ParseUnsigned64(value, &o->backgroundCadenceSteps)) return invalidValue();
  } else if (field == "application.sep.injection_cadence_steps") {
    if (!ParseUnsigned64(value, &o->injectionCadenceSteps)) return invalidValue();
  } else if (field == "mesh.domain.preset") {
    if (!ParseEnum(value, {{"solar", DomainPreset::Solar},
                           {"one-au", DomainPreset::OneAu},
                           {"earth", DomainPreset::OneAu},
                           {"mars", DomainPreset::Mars}}, &o->domain))
      return invalidValue();
  } else if (field == "mesh.domain.inner_radius_m") {
    if (!ParseDouble(value, &o->innerRadiusM)) return invalidValue();
  } else if (field == "mesh.domain.inner_boundary") {
    if (!ParseEnum(value, {{"absorb", InnerBoundaryMode::Absorb}},
                   &o->innerBoundary)) return invalidValue();
  } else if (field == "mesh.domain.outer_radius_mode") {
    if (!ParseEnum(value, {{"preset", OuterRadiusMode::Preset},
                           {"explicit", OuterRadiusMode::Explicit},
                           {"field-line-endpoint", OuterRadiusMode::FieldLineEndpoint}},
                   &o->outerRadiusMode)) return invalidValue();
  } else if (field == "mesh.domain.outer_radius_m") {
    if (!ParseDouble(value, &o->outerRadiusM)) return invalidValue();
  } else if (field == "mesh.domain.box_geometry") {
    if (!ParseEnum(value, {{"sun-centered-cube", DomainBoxGeometry::SunCenteredCube},
                           {"field-line-corner-cube", DomainBoxGeometry::FieldLineCornerCube},
                           {"field-line-xy-corner-cube", DomainBoxGeometry::FieldLineXYCornerCube}},
                   &o->domainBoxGeometry)) return invalidValue();
  } else if (field == "mesh.domain.corner_direction_x") {
    if (!ParseDouble(value, &o->domainCornerDirection.x)) return invalidValue();
  } else if (field == "mesh.domain.corner_direction_y") {
    if (!ParseDouble(value, &o->domainCornerDirection.y)) return invalidValue();
  } else if (field == "mesh.domain.corner_direction_z") {
    if (!ParseDouble(value, &o->domainCornerDirection.z)) return invalidValue();
  } else if (field == "mesh.domain.corner_margin_m") {
    if (!ParseDouble(value, &o->domainCornerMarginM)) return invalidValue();
  } else if (field == "mesh.field_line.end_radius_m") {
    if (!ParseDouble(value, &o->parkerSpiralEndRadiusM)) return invalidValue();
  } else if (field == "mesh.domain.outer_boundary") {
    if (!ParseEnum(value, {{"escape", OuterBoundaryMode::Escape},
                           {"imported-coverage", OuterBoundaryMode::ImportedCoverage}},
                   &o->outerBoundary)) return invalidValue();
  } else if (field == "mesh.domain.coordinate_frame") {
    o->coordinateFrame = value;
    o->parker.coordinateFrame = value;
  } else if (field == "mesh.domain.origin_x_m") {
    if (!ParseDouble(value, &o->coordinateOriginM.x)) return invalidValue();
  } else if (field == "mesh.domain.origin_y_m") {
    if (!ParseDouble(value, &o->coordinateOriginM.y)) return invalidValue();
  } else if (field == "mesh.domain.origin_z_m") {
    if (!ParseDouble(value, &o->coordinateOriginM.z)) return invalidValue();
  } else if (field == "mesh.field_line.origin_x_m") {
    if (!ParseDouble(value, &o->parkerSpiralOriginM.x)) return invalidValue();
  } else if (field == "mesh.field_line.origin_y_m") {
    if (!ParseDouble(value, &o->parkerSpiralOriginM.y)) return invalidValue();
  } else if (field == "mesh.field_line.origin_z_m") {
    if (!ParseDouble(value, &o->parkerSpiralOriginM.z)) return invalidValue();
  } else if (field == "mesh.field_line.start_mode") {
    if (!ParseEnum(value,
        {{"explicit", ParkerSpiralStartMode::Explicit},
         {"cme-launch-point", ParkerSpiralStartMode::CmeLaunchPoint}},
        &o->parkerSpiralStartMode)) return invalidValue();
  } else if (field == "mesh.field_line.initial_x_m") {
    if (!ParseDouble(value, &o->parkerSpiralInitialPointM.x)) return invalidValue();
  } else if (field == "mesh.field_line.initial_y_m") {
    if (!ParseDouble(value, &o->parkerSpiralInitialPointM.y)) return invalidValue();
  } else if (field == "mesh.field_line.initial_z_m") {
    if (!ParseDouble(value, &o->parkerSpiralInitialPointM.z)) return invalidValue();
  } else if (field == "mesh.field_line.length_m") {
    if (!ParseDouble(value, &o->parkerSpiralLengthM)) return invalidValue();
  } else if (field == "mesh.field_line.point_count") {
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
  } else if (field == "mesh.refinement.solar.enabled") {
    if (!ParseBool(value, &o->enableRadialRefinement)) return invalidValue();
  } else if (field == "mesh.refinement.solar.surface_cell_size_m") {
    if (!ParseDouble(value, &o->solarSurfaceCellSizeM)) return invalidValue();
  } else if (field == "mesh.refinement.solar.transition_outer_radius_m") {
    if (!ParseDouble(value, &o->solarRefinementOuterRadiusM)) return invalidValue();
  } else if (field == "mesh.refinement.solar.profile") {
    if (!ParseEnum(value, {{"linear", RefinementProfile::Linear},
                           {"power-law", RefinementProfile::PowerLaw},
                           {"smoothstep", RefinementProfile::Smoothstep}},
                   &o->solarRefinementProfile)) return invalidValue();
  } else if (field == "mesh.refinement.solar.anchor") {
    if (!ParseEnum(value, {{"source-shell", SolarRefinementAnchor::SourceShell},
                           {"photosphere", SolarRefinementAnchor::Photosphere}},
                   &o->solarRefinementAnchor)) return invalidValue();
  } else if (field == "mesh.refinement.solar.exponent") {
    if (!ParseDouble(value, &o->solarRefinementExponent)) return invalidValue();
  } else if (field == "mesh.refinement.tube.enabled") {
    if (!ParseBool(value, &o->enableTubeRefinement)) return invalidValue();
  } else if (field == "mesh.refinement.tube.source_longitude_rad") {
    if (!ParseDouble(value, &o->tubeLongitudeRad)) return invalidValue();
  } else if (field == "mesh.refinement.tube.source_colatitude_rad") {
    if (!ParseDouble(value, &o->tubeColatitudeRad)) return invalidValue();
  } else if (field == "mesh.refinement.tube.reference_radius_m") {
    if (!ParseDouble(value, &o->tubeReferenceRadiusM)) return invalidValue();
  } else if (field == "mesh.refinement.tube.radius_at_reference_m") {
    if (!ParseDouble(value, &o->tubeRadiusAtReferenceM)) return invalidValue();
  } else if (field == "mesh.refinement.tube.radius_mode") {
    if (!ParseEnum(value,
        {{"physical-constant", TubeRadiusMode::PhysicalConstant},
         {"constant-angular-width", TubeRadiusMode::ConstantAngularWidth}},
        &o->tubeRadiusMode)) return invalidValue();
  } else if (field == "mesh.refinement.tube.center_cell_size_m") {
    if (!ParseDouble(value, &o->tubeCellSizeM)) return invalidValue();
  } else if (field == "mesh.refinement.tube.transverse_profile") {
    if (!ParseEnum(value, {{"linear", RefinementProfile::Linear},
                           {"power-law", RefinementProfile::PowerLaw},
                           {"smoothstep", RefinementProfile::Smoothstep}},
                   &o->tubeTransverseProfile)) return invalidValue();
  } else if (field == "mesh.refinement.tube.transverse_exponent") {
    if (!ParseDouble(value, &o->tubeTransverseExponent)) return invalidValue();
  } else if (field == "mesh.active_region.mode") {
    if (!ParseEnum(value,
        {{"full-domain", ActiveRegionMode::FullDomain},
         {"parker-tube", ActiveRegionMode::ParkerTube}},
        &o->activeRegion)) return invalidValue();
  } else if (field == "mesh.active_region.reference_radius_m") {
    if (!ParseDouble(value, &o->activeTubeReferenceRadiusM))
      return invalidValue();
  } else if (field == "mesh.active_region.radius_at_reference_m") {
    if (!ParseDouble(value, &o->activeTubeRadiusAtReferenceM))
      return invalidValue();
  } else if (field == "mesh.active_region.radius_mode") {
    if (!ParseEnum(value,
        {{"physical-constant", TubeRadiusMode::PhysicalConstant},
         {"constant-angular-width", TubeRadiusMode::ConstantAngularWidth}},
        &o->activeTubeRadiusMode)) return invalidValue();
  } else if (field == "mesh.active_region.buffer_blocks") {
    if (!ParseUnsigned(value, &o->activeTubeBufferBlocks))
      return invalidValue();
  } else if (field == "mesh.active_region.solar_sphere_radius_m") {
    if (!ParseDouble(value, &o->activeSolarSphereRadiusM)) return invalidValue();
  } else if (field == "mesh.memory.base_cell_bytes") {
    if (!ParseSize(value, &o->memoryModel.baseCellBytes)) return invalidValue();
  } else if (field == "mesh.memory.base_node_bytes") {
    if (!ParseSize(value, &o->memoryModel.baseNodeBytes)) return invalidValue();
  } else if (field == "mesh.memory.block_structure_bytes") {
    if (!ParseSize(value, &o->memoryModel.blockStructureBytes)) return invalidValue();
  } else if (field == "mesh.memory.communication_bytes_per_block") {
    if (!ParseSize(value, &o->memoryModel.communicationBytesPerBlock)) return invalidValue();
  } else if (field == "mesh.memory.particle_bytes") {
    if (!ParseSize(value, &o->memoryModel.particleBytes)) return invalidValue();
  } else if (field == "mesh.memory.particles_per_cell") {
    if (!ParseDouble(value, &o->memoryModel.particlesPerCell)) return invalidValue();
  } else if (field == "mesh.memory.halo_fraction") {
    if (!ParseDouble(value, &o->memoryModel.haloFraction)) return invalidValue();
  } else if (field == "mesh.memory.safety_margin_fraction") {
    if (!ParseDouble(value, &o->memoryModel.safetyMarginFraction)) return invalidValue();
  } else if (field == "physics.background.provider") {
    if (!ParseEnum(value,
        {{"analytic-parker", BackgroundAuthority::AnalyticParker},
         {"python-interpolator", BackgroundAuthority::PythonInterpolator},
         {"swmf", BackgroundAuthority::Swmf},
         {"swcme", BackgroundAuthority::Swcme},
         {"runtime-model", BackgroundAuthority::RuntimeModel}}, &o->background))
      return invalidValue();
  } else if (field == "physics.background.model_id") {
    // Preserve the exact factory key; immutable configuration validation later
    // requires it for runtime-model and forbids it for every built-in authority.
    o->backgroundModelId = value;
  } else if (field == "physics.background.parker.reference_radius_m") {
    if (!ParseDouble(value, &o->parker.referenceRadiusM)) return invalidValue();
  } else if (field == "physics.background.parker.radial_field_at_reference_t") {
    if (!ParseDouble(value, &o->parker.radialFieldAtReferenceT)) return invalidValue();
  } else if (field == "physics.background.parker.solar_rotation_rate_rad_per_s") {
    if (!ParseDouble(value, &o->parker.solarRotationRateRadPerS)) return invalidValue();
  } else if (field == "physics.background.parker.solar_wind_speed_m_per_s") {
    if (!ParseDouble(value, &o->parker.solarWindSpeedMPerS)) return invalidValue();
  } else if (field == "physics.background.parker.magnetic_polarity") {
    double parsed = 0.0;
    if (!ParseDouble(value, &parsed) || (parsed != -1.0 && parsed != 1.0))
      return invalidValue();
    o->parker.magneticPolarity = static_cast<int>(parsed);
  } else if (field == "physics.background.parker.number_density_at_one_au_m3") {
    // Electron number density at 1 AU [m^-3]; the density reference is fixed
    // at 1 AU independently of reference_radius_m (which expresses B_r).
    if (!ParseDouble(value, &o->parker.numberDensityAtReferenceM3))
      return invalidValue();
    o->parker.densityReferenceRadiusM = Core::Const::AU;
  } else if (field == "physics.background.parker.temperature_k") {
    if (!ParseDouble(value, &o->parker.temperatureK)) return invalidValue();
  } else if (field == "physics.background.parker.validity_cadence_s") {
    if (!ParseDouble(value, &o->parker.validityCadenceS)) return invalidValue();
  } else if (field == "application.sep.turbulence.authority") {
    if (!ParseEnum(value, {{"prescribed", TurbulenceAuthority::Prescribed},
                           {"swmf", TurbulenceAuthority::Swmf}},
                   &o->turbulence)) return invalidValue();
  } else if (field == "application.sep.turbulence.model") {
    if (!ParseEnum(value,
        {{"power-law", PrescribedTurbulenceModel::PowerLaw},
         {"kolmogorov", PrescribedTurbulenceModel::Kolmogorov},
         {"kraichnan", PrescribedTurbulenceModel::Kraichnan}},
        &o->prescribedTurbulenceModel)) return invalidValue();
  } else if (field == "application.sep.turbulence.amplitude_model") {
    if (!ParseEnum(value,
        {{"constant-delta-b-over-b",
              PrescribedTurbulenceAmplitudeModel::ConstantDeltaBOverB},
         {"wave-energy-power-law",
              PrescribedTurbulenceAmplitudeModel::WaveEnergyPowerLaw}},
        &o->prescribedTurbulenceAmplitudeModel)) return invalidValue();
  } else if (field == "application.sep.turbulence.delta_b_over_b") {
    if (!ParseDouble(value, &o->prescribedDeltaBOverB)) return invalidValue();
  } else if (field ==
             "application.sep.turbulence.wave_energy_density_at_reference_j_per_m3") {
    if (!ParseDouble(value, &o->turbulenceWaveEnergyAtReferenceJPerM3))
      return invalidValue();
  } else if (field == "application.sep.turbulence.wave_energy_density_radial_exponent") {
    if (!ParseDouble(value, &o->turbulenceWaveEnergyRadialExponent))
      return invalidValue();
  } else if (field == "application.sep.turbulence.normalized_cross_helicity") {
    if (!ParseDouble(value, &o->turbulenceNormalizedCrossHelicity))
      return invalidValue();
  } else if (field == "application.sep.turbulence.reference_radius_m") {
    if (!ParseDouble(value, &o->turbulenceReferenceRadiusM))
      return invalidValue();
  } else if (field == "application.sep.turbulence.k_min_per_m") {
    if (!ParseDouble(value, &o->turbulenceKMinPerM)) return invalidValue();
  } else if (field == "application.sep.turbulence.k_max_per_m") {
    if (!ParseDouble(value, &o->turbulenceKMaxPerM)) return invalidValue();
  } else if (field == "application.sep.turbulence.k_min_radial_exponent") {
    if (!ParseDouble(value, &o->turbulenceKMinRadialExponent))
      return invalidValue();
  } else if (field == "application.sep.turbulence.k_max_radial_exponent") {
    if (!ParseDouble(value, &o->turbulenceKMaxRadialExponent))
      return invalidValue();
  } else if (field == "application.sep.turbulence.spectral_index") {
    if (!ParseDouble(value, &o->turbulenceSpectralIndex)) return invalidValue();
  } else if (field == "application.sep.turbulence.correlation_length_m") {
    if (!ParseDouble(value, &o->turbulenceCorrelationLengthM)) return invalidValue();
  } else if (field == "application.sep.turbulence.correlation_length_radial_exponent") {
    if (!ParseDouble(value, &o->turbulenceCorrelationLengthRadialExponent))
      return invalidValue();
  } else if (field == "application.sep.turbulence.validity_cadence_s") {
    if (!ParseDouble(value, &o->turbulenceValidityCadenceS))
      return invalidValue();
  } else if (field == "application.sep.turbulence.missing_data") {
    if (!ParseEnum(value, {{"fail", MissingTurbulenceMode::Fail},
                           {"ballistic", MissingTurbulenceMode::Ballistic}},
                   &o->missingTurbulence)) return invalidValue();
  } else if (field == "application.sep.turbulence.resonance_range") {
    if (!ParseEnum(value, {{"reject", ResonanceRangeMode::Reject},
                           {"power-law-extension", ResonanceRangeMode::PowerLawExtension}},
                   &o->resonanceRange)) return invalidValue();
  } else if (field == "application.sep.turbulence.self_consistent_3d") {
    if (!ParseBool(value, &o->enableSelfConsistent3DTurbulence)) return invalidValue();
  } else if (field == "application.sep.transport.cell_crossing_fraction") {
    if (!ParseDouble(value, &o->cellCrossingFraction)) return invalidValue();
  } else if (field == "application.sep.transport.diffusion_fraction") {
    if (!ParseDouble(value, &o->diffusionFraction)) return invalidValue();
  } else if (field == "application.sep.transport.focusing_fraction") {
    if (!ParseDouble(value, &o->focusingFraction)) return invalidValue();
  } else if (field == "application.sep.transport.cooling_fraction") {
    if (!ParseDouble(value, &o->coolingFraction)) return invalidValue();
  } else if (field == "application.sep.transport.field_variation_fraction") {
    if (!ParseDouble(value, &o->fieldVariationFraction)) return invalidValue();
  } else if (field == "application.sep.transport.shock_crossing_fraction") {
    if (!ParseDouble(value, &o->shockCrossingFraction)) return invalidValue();
  } else if (field == "application.sep.transport.minimum_substep_s") {
    if (!ParseDouble(value, &o->minimumTransportSubstepS)) return invalidValue();
  } else if (field == "application.sep.transport.maximum_substeps") {
    if (!ParseUnsigned64(value, &o->maximumTransportSubsteps)) return invalidValue();
  } else if (field == "application.sep.transport.pitch_angle_scheme") {
    if (!ParseEnum(value,
        {{"reflecting-milstein", PitchAngleSchemeMode::ReflectingMilstein},
         {"reflecting-euler-maruyama",
          PitchAngleSchemeMode::ReflectingEulerMaruyama}},
        &o->pitchAngleScheme)) return invalidValue();
  } else if (field == "application.sep.transport.spatial_diffusion_model") {
    if (!ParseEnum(value,
        {{"mean-free-path",
              SpatialDiffusionModel::CorrelationMeanFreePath},
         {"correlation-mean-free-path",
              SpatialDiffusionModel::CorrelationMeanFreePath},
         {"pitch-angle-integral",
              SpatialDiffusionModel::PitchAngleIntegral},
         {"parallel-diffusion-library",
              SpatialDiffusionModel::ParallelDiffusionLibrary}},
        &o->spatialDiffusionModel)) return invalidValue();
  } else if (field == "application.sep.transport.pitch_angle_diffusion_model") {
    if (!ParseEnum(value,
        {{"jokipii-1966", PitchAngleDiffusionModel::Jokipii1966},
         {"florinskiy", PitchAngleDiffusionModel::Florinskiy},
         {"constant", PitchAngleDiffusionModel::Constant}},
        &o->pitchAngleDiffusionModel)) return invalidValue();
  } else if (field == "application.sep.transport.mean_free_path_model") {
    if (!ParseEnum(value,
        {{"correlation", MeanFreePathModel::Correlation},
         {"constant", MeanFreePathModel::Constant},
         {"radial-rigidity-power-law",
              MeanFreePathModel::RadialRigidityPowerLaw}},
        &o->meanFreePathModel)) return invalidValue();
  } else if (field == "application.sep.transport.constant_dmumu_per_s") {
    if (!ParseDouble(value, &o->constantDmumuPerS)) return invalidValue();
  } else if (field == "application.sep.transport.constant_mean_free_path_m") {
    if (!ParseDouble(value, &o->constantMeanFreePathM)) return invalidValue();
  } else if (field == "application.sep.transport.mean_free_path_reference_m") {
    if (!ParseDouble(value, &o->meanFreePathReferenceM)) return invalidValue();
  } else if (field == "application.sep.transport.mean_free_path_reference_radius_m") {
    if (!ParseDouble(value, &o->meanFreePathReferenceRadiusM))
      return invalidValue();
  } else if (field == "application.sep.transport.mean_free_path_reference_rigidity_v") {
    if (!ParseDouble(value, &o->meanFreePathReferenceRigidityV))
      return invalidValue();
  } else if (field == "application.sep.transport.mean_free_path_radial_exponent") {
    if (!ParseDouble(value, &o->meanFreePathRadialExponent))
      return invalidValue();
  } else if (field == "application.sep.transport.mean_free_path_rigidity_exponent") {
    if (!ParseDouble(value, &o->meanFreePathRigidityExponent))
      return invalidValue();
  } else if (field ==
             "application.sep.transport.spatial_quadrature_absolute_tolerance_m2_per_s") {
    if (!ParseDouble(value, &o->spatialQuadratureAbsoluteToleranceM2PerS))
      return invalidValue();
  } else if (field == "application.sep.transport.spatial_quadrature_relative_tolerance") {
    if (!ParseDouble(value, &o->spatialQuadratureRelativeTolerance))
      return invalidValue();
  } else if (field == "application.sep.transport.spatial_quadrature_maximum_recursion") {
    if (!ParseUnsigned(value, &o->spatialQuadratureMaximumRecursion))
      return invalidValue();
  } else if (field == "application.sep.transport.focused_scattering_frame") {
    if (!ParseEnum(value,
        {{"plasma-frame-isotropic",
              FocusedScatteringFrame::PlasmaFrameIsotropic},
         {"alfven-wave-frame-isotropic",
              FocusedScatteringFrame::AlfvenWaveFrameIsotropic}},
        &o->focusedScatteringFrame)) return invalidValue();
  } else if (field == "application.sep.transport.maximum_scattering_events_per_substep") {
    if (!ParseUnsigned64(value, &o->maximumScatteringEventsPerSubstep))
      return invalidValue();
  } else if (field == "application.sep.transport.perpendicular_diffusion") {
    if (!ParseEnum(value,
        {{"none", PerpendicularDiffusionMode::None},
         {"constant", PerpendicularDiffusionMode::Constant},
         {"constant-ratio", PerpendicularDiffusionMode::ConstantRatio}},
        &o->perpendicularDiffusion)) return invalidValue();
  } else if (field == "application.sep.transport.constant_kappa_perpendicular_m2_per_s") {
    if (!ParseDouble(value, &o->constantKappaPerpendicularM2PerS))
      return invalidValue();
  } else if (field == "application.sep.transport.kappa_perpendicular_to_parallel_ratio") {
    if (!ParseDouble(value, &o->kappaPerpendicularToParallelRatio))
      return invalidValue();
  } else if (field == "application.sep.transport.drifts") {
    if (!ParseEnum(value,
        {{"none", DriftMode::None}, {"gradient-b", DriftMode::GradientB},
         {"curvature", DriftMode::Curvature},
         {"gradient-curvature", DriftMode::GradientAndCurvature}},
        &o->drift)) return invalidValue();
  } else if (field == "application.sep.population_control.mode") {
    if (!ParseEnum(value,
        {{"off", PopulationControlMode::Off},
         {"split-merge", PopulationControlMode::SplitMerge}},
        &o->populationControl)) return invalidValue();
  } else if (field ==
             "application.sep.population_control.minimum_particles_per_cell_per_species") {
    if (!ParseUnsigned(value, &o->minimumParticlesPerCellPerSpecies))
      return invalidValue();
  } else if (field ==
             "application.sep.population_control.target_particles_per_cell_per_species") {
    if (!ParseUnsigned(value, &o->targetParticlesPerCellPerSpecies))
      return invalidValue();
  } else if (field ==
             "application.sep.population_control.maximum_particles_per_cell_per_species") {
    if (!ParseUnsigned(value, &o->maximumParticlesPerCellPerSpecies))
      return invalidValue();
  } else if (field == "application.sep.population_control.cadence_steps") {
    if (!ParseUnsigned64(value, &o->populationControlCadenceSteps))
      return invalidValue();
  } else if (field == "physics.shock.authority") {
    if (!ParseEnum(value, {{"none", ShockAuthority::None},
                           {"swcme", ShockAuthority::Swcme}}, &o->shock))
      return invalidValue();
  } else if (field == "source.shock.enabled") {
    if (!ParseBool(value, &o->source.enabled)) return invalidValue();
  } else if (field == "source.shock.physical_particle_rate_per_s") {
    if (!ParseDouble(value, &o->source.physicalParticleRatePerS)) return invalidValue();
  } else if (field == "source.shock.injection_efficiency") {
    // Dimensionless fraction; with physics: swcme it is also the SWCME
    // source.injection_efficiency (D4).
    if (!ParseDouble(value, &o->source.injectionEfficiency)) return invalidValue();
    spectrum->injectionEfficiency = o->source.injectionEfficiency;
  } else if (field == "source.shock.minimum_energy_j") {
    // Kinetic-energy limits [J] of the injected spectrum.  With physics:
    // swcme they are also the SWCME source.energy_min/max (D4).
    if (!ParseDouble(value, &o->source.minimumEnergyJ)) return invalidValue();
    spectrum->minimumEnergyJ = o->source.minimumEnergyJ;
  } else if (field == "source.shock.maximum_energy_j") {
    if (!ParseDouble(value, &o->source.maximumEnergyJ)) return invalidValue();
    spectrum->maximumEnergyJ = o->source.maximumEnergyJ;
  } else if (field == "source.shock.spectrum_reference_energy_j") {
    // Reference spectrum handed to the canonical SWCME source (D4).  srcSEP3D
    // never takes species identity from it: injection converts the energy
    // limits to momentum with each compiled AMPS species mass.
    if (!ParseDouble(value, &spectrum->referenceEnergyJ)) return invalidValue();
  } else if (field == "source.shock.spectrum_particle_mass_kg") {
    if (!ParseDouble(value, &spectrum->particleMassKg)) return invalidValue();
  } else if (field == "source.shock.spectrum_charge_number") {
    std::uint64_t charge = 0;
    if (!ParseUnsigned64(value, &charge) || charge == 0 ||
        charge > static_cast<std::uint64_t>(std::numeric_limits<long>::max()))
      return invalidValue();
    spectrum->chargeNumber = static_cast<long>(charge);
  } else if (field == "source.shock.spectrum_normalization") {
    const std::string normalized = Lower(value);
    if (normalized != "relative-only" &&
        normalized != "reference-differential-intensity") return invalidValue();
    spectrum->normalization = normalized;
  } else if (field ==
             "source.shock.spectrum_reference_intensity_per_m2_s_sr_j") {
    if (!ParseDouble(value, &spectrum->referenceIntensitySi)) return invalidValue();
  } else if (field == "source.shock.spectrum_model") {
    if (!ParseEnum(value,
        {{"local-compression-dsa",
          SourceSpectrumModel::LocalCompressionDsa},
         {"fixed-phase-space-power-law",
          SourceSpectrumModel::FixedPhaseSpacePowerLaw}},
        &o->source.spectrumModel)) return invalidValue();
  } else if (field == "source.shock.phase_space_power_index") {
    if (!ParseDouble(value, &o->source.fixedPhaseSpacePowerIndex))
      return invalidValue();
  } else if (field == "source.shock.samples_per_step") {
    if (!ParseUnsigned64(value, &o->source.samplesPerStep)) return invalidValue();
  } else if (field == "source.shock.normalization") {
    // configured-rate: physical_particle_rate_per_s * injection_efficiency
    // and species.macroparticle_weight define the source (SWCME source).
    // accepted-shock-incident-flux: Ndot_s = sum over accepted reduced-front
    // faces of n_s [V_n - U_1.n]_+ A, W_s = Ndot_s dt / samples_per_step.
    if (!ParseEnum(value,
        {{"configured-rate", SourceRateNormalizationModel::ConfiguredConstant},
         {"accepted-shock-incident-flux",
          SourceRateNormalizationModel::AcceptedShockIncidentFlux}},
        &o->particleNumerics.sourceRateModel)) return invalidValue();
  } else if (field == "source.shock.normalization_radius_m") {
    // Heliocentric apex radius [m] at which Ndot_s (and so W_s) is evaluated.
    if (!ParseDouble(value, &o->particleNumerics.sourceNormalizationRadiusM))
      return invalidValue();
  } else if (field == "source.shock.weight_model") {
    // constant-statistical-weight is implemented; the importance mode is a
    // typed reserved feature that stops startup until it is qualified.
    if (!ParseEnum(value,
        {{"constant-statistical-weight",
          SourceWeightingModel::ConstantStatisticalWeight},
         {"log-uniform-momentum-importance",
          SourceWeightingModel::LogUniformMomentumImportance}},
        &o->source.weightingModel)) return invalidValue();
  } else if (field == "source.shock.maximum_events_per_species_per_step") {
    // Fatal Poisson runaway guard; never truncates or renormalizes a draw.
    if (!ParseUnsigned64(value, &o->source.maximumMacroparticlesPerSpeciesPerStep))
      return invalidValue();
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
  } else if (field == "output.initialization_field_line_tecplot_file") {
    o->initializationParkerLineTecplotFile = value;
  } else if (field == "output.initialization_data_tecplot_file") {
    o->initializationDataTecplotFile = value;
  } else if (field == "application.sep.storage.magnetic_gradient") {
    if (!ParseBool(value, &o->storeMagneticGradient)) return invalidValue();
  } else if (field == "application.sep.storage.velocity_gradient") {
    if (!ParseBool(value, &o->storeVelocityGradient)) return invalidValue();
  } else if (field == "application.sep.storage.sampling_bytes_per_cell") {
    if (!ParseSize(value, &o->samplingBytesPerCell)) return invalidValue();
  } else if (field == "run.restart.input_path") {
    // Text files cannot encode an empty right-hand side because that usually
    // indicates a truncated edit.  The explicit token "none" denotes a fresh
    // run and normalizes to the same empty string used by typed construction.
    o->restartInputPath = Lower(value) == "none" ? std::string() : value;
  } else if (field == "run.restart.output_path") {
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
      // field-line-radius is the shared-vocabulary kind of srcSEP's one-
      // dimensional observers (decision D5); srcSEP3D samples in space.
      if (Lower(value) == "field-line-radius")
        return Invalid(SEP::InputLayout::NotSupportedMessage(
            "observer kind 'field-line-radius' ('" + field + "')",
            SEP::InputLayout::kSrcSep3d));
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

// ---------------------------------------------------------------------------
// Generic AMPS runtime input -> srcSEP3D options
// ---------------------------------------------------------------------------
//
// The srcSEP3D input is written in the generic AMPS runtime-input format
// (src/INPUT_FILE.md; reader src/general/input/amps_input.h) using the SEP
// content schema 6 that srcSEP3D shares with srcSEP.  The complete block and
// key vocabulary of that schema is the table in
// src/models/sep_common/sep_input_layout.h; this walker consults it for every
// block and key, so the two applications accept exactly the same spellings:
//
//   * a block or key that is not in the table is unknown and fatal at its
//     '#begin' or assignment line (unknown input is never ignored);
//   * a block or key that the table marks srcSEP-only is a typed
//     "not supported by srcSEP3D" error (decision D3), e.g. an observer's
//     heliocentric_radius_m;
//   * everything else reaches ApplyField() under its qualified name.
//
// Qualified names.  A key's qualified name is the dot-joined path of the
// blocks that contain it, where a named block '#begin kind: name' contributes
// 'kind.name': 'mesh.domain.inner_radius_m', 'mesh.refinement.solar.profile',
// 'application.sep.turbulence.k_min_per_m', 'observer.earth.cadence_s'.
// ApplyField() and every diagnostic use these names, so an error message
// names exactly the block path a user has to edit.
//
// Model blocks.  Three physics blocks belong to shared model libraries and
// keep those libraries' own vocabularies:
//   physics: swcme              canonical SWCME assignments (dotted keys with
//                               units in the value, e.g. "400 km/s"),
//                               resolved by swcme::input3d::Resolve; its
//                               'source.*' keys are owned by source: shock
//                               (decision D4) and are synthesized below;
//   physics: parallel_diffusion parsed by SEP::ParallelDiffusion::ParseSection;
//   physics: shock_front        reduced shock-front event, resolved by
//                               SEP::CoronaSwcme::ShockFront (usually written
//                               in its own file and #included).
// All other keys are plain identifiers; a dotted key in an application block
// is rejected because nested blocks express that path.

namespace Layout = SEP::InputLayout;

// Accumulated state of one document walk.
struct WalkState {
  RunConfiguration3DOptions* options = nullptr;
  std::set<std::string> assigned;  // qualified keys
  std::set<std::string> sections;  // qualified block paths
  std::map<std::string, std::size_t> keyLines;  // qualified key -> input line
  bool observerDefaultsCleared = false;
  std::vector<SEP::ParallelDiffusion::SectionLine> parallelDiffusionLines;
  std::string parallelDiffusionFile;
  const AMPS::Input::Block* shockFrontBlock = nullptr;
  // source: shock reference spectrum (D4), turned into SWCME source.* keys
  // after the walk when physics: swcme is present.
  Layout::SwcmeSourceSettings spectrum;
  std::size_t physicsSwcmeLine = 0;  // '#begin physics: swcme' line, 0 = absent
};

Core::Status WalkBlocks(const AMPS::Input::Block& parent,
                        const std::string& parentPath, WalkState* state) {
  for (const AMPS::Input::Block& block : parent.blocks) {
    const Layout::BlockRule* rule =
        Layout::FindBlock(parentPath, block.kind, block.name);
    if (rule == nullptr) {
      const std::string where = parentPath.empty()
          ? std::string("at the top level") : "inside '" + parentPath + "'";
      return Invalid(AMPS::Input::ErrorAt(block.begin, "block '" +
          block.Label() + "' is not part of the SEP input " + where +
          " (allowed: " + Layout::AllowedChildren(parentPath) + ")").message);
    }
    const std::string path =
        Layout::QualifiedPath(parentPath, block.kind, block.name);
    if ((rule->support & Layout::kSrcSep3d) == 0)
      return Invalid(AMPS::Input::ErrorAt(block.begin,
          Layout::NotSupportedMessage("block '" + path + "'",
                                      Layout::kSrcSep3d)).message);
    state->sections.insert(path);

    if (path == "physics.swcme") {
      // SWCME owns its dotted, unit-bearing vocabulary; nested blocks are the
      // same key space (AMPS::Input::Flatten).  Locations are kept for the
      // canonical resolver's key/line diagnostics.  Its source.* keys are
      // written once, in source: shock (D4).
      state->physicsSwcmeLine = block.begin.line;
      for (const AMPS::Input::Assignment& item : AMPS::Input::Flatten(block)) {
        if (const char* owner = Layout::SourceShockKeyForSwcmeKey(item.key))
          return Invalid(AMPS::Input::ErrorAt(item.where, "SWCME key '" +
              item.key + "' is written in source: shock as '" + owner +
              "' (each quantity is set in one place, src/INPUT_FILE.md "
              "section 5)", item.text).message);
        SwcmeAssignment assignment;
        assignment.key = item.key;
        assignment.value = item.value;
        assignment.line = item.where.line;
        state->options->swcmeAssignments.push_back(assignment);
        state->assigned.insert(path + "." + item.key);
      }
      continue;
    }
    if (path == "physics.parallel_diffusion") {
      // The library parser owns this block's grammar and model schema; it
      // receives "key = value" records with their physical line numbers.
      if (!block.blocks.empty())
        return Invalid(AMPS::Input::ErrorAt(block.blocks.front().begin,
            "physics: parallel_diffusion holds only key = value records").message);
      for (const AMPS::Input::Assignment& item : block.assignments)
        state->parallelDiffusionLines.push_back(
            {item.where.line, item.key + " = " + item.value});
      state->parallelDiffusionFile = block.begin.file;
      continue;
    }
    if (path == "physics.shock_front") {
      // Resolved later by the shared reduced shock-front library from the
      // canonical text of this block (see ParseConfigurationDocument).
      state->shockFrontBlock = &block;
      continue;
    }

    for (const AMPS::Input::Assignment& item : block.assignments) {
      if (item.key.find('.') != std::string::npos)
        return Invalid(AMPS::Input::ErrorAt(item.where, "key '" + item.key +
            "' contains '.'; SEP keys are plain names and a nested "
            "path is written as nested blocks", item.text).message);
      const std::string qualified = path + "." + item.key;
      const Layout::KeyRule* key = Layout::FindKey(path, item.key);
      if (key == nullptr)
        return Invalid(AMPS::Input::ErrorAt(item.where,
            "unknown configuration key '" + qualified + "'", item.text).message);
      if ((key->support & Layout::kSrcSep3d) == 0)
        return Invalid(AMPS::Input::ErrorAt(item.where,
            Layout::NotSupportedMessage("'" + qualified + "'",
                                        Layout::kSrcSep3d), item.text).message);
      state->assigned.insert(qualified);
      state->keyLines[qualified] = item.where.line;
      const Core::Status status = ApplyField(path, item.key, item.value,
          state->options, &state->observerDefaultsCleared, &state->spectrum);
      if (!status.ok())
        return Invalid(AMPS::Input::ErrorAt(item.where, status.message,
                                            item.text).message);
    }
    const Core::Status nested = WalkBlocks(block, path, state);
    if (!nested.ok()) return nested;
  }
  return Core::Status::OK();
}

// Canonical generic-input text of the physics: shock_front block, passed to
// the shared resolver.  Its file:line provenance is replaced by key names in
// the resolver's diagnostics; the identity is the resolver's normalized
// manifest, so comments, layout and include structure do not change it.
std::string ShockFrontCanonicalText(const AMPS::Input::Block& block) {
  AMPS::Input::Block root;
  root.blocks.push_back(block);
  return AMPS::Input::ToCanonicalText(root);
}

}  // namespace

static Core::Status ParseConfigurationDocument(
    const AMPS::Input::Document& document, RunConfiguration3DOptions* result) {
  if (result == nullptr) return Invalid("configuration output is null");
  RunConfiguration3DOptions candidate;
  WalkState walk;
  walk.options = &candidate;
  {
    const Core::Status status = WalkBlocks(document.root, "", &walk);
    if (!status.ok()) return status;
  }
  const std::set<std::string>& assigned = walk.assigned;
  const std::set<std::string>& sections = walk.sections;
  if (assigned.count("application.sep.schema_version") == 0)
    return Invalid("missing required key application.sep.schema_version "
                   "(write schema_version = 6 in '#begin application: sep')");

  // ---- Particle numerics derived from the mesh and the reduced shock ----
  //
  // run.time_step_mode = derived-from-mesh selects dt = f h_min / v_max,
  // evaluated after AMPS allocates the mesh, and source.shock.normalization
  // = accepted-shock-incident-flux selects W_s = Ndot_s dt / N per compiled
  // species from the accepted upstream incident flux of the reduced front
  // (src/INPUT_FILE.md, "Derived particle numerics").  The keys that the
  // derivation replaces must then be absent, so that no written value can be
  // mistaken for the one the run uses; conversely the derivation inputs are
  // rejected in the configured modes.
  const bool derivedTimeStep = candidate.particleNumerics.deriveFromMeshAndShock;
  const bool incidentFlux = candidate.particleNumerics.sourceRateModel ==
      SourceRateNormalizationModel::AcceptedShockIncidentFlux;
  const std::set<std::string> derivedOwned = {
      "run.time_step_s", "source.shock.physical_particle_rate_per_s",
      "source.shock.injection_efficiency", "species.macroparticle_weight"};
  if (derivedTimeStep != incidentFlux)
    return Invalid("run.time_step_mode=derived-from-mesh and "
                   "source.shock.normalization=accepted-shock-incident-flux "
                   "must be selected together");
  if (incidentFlux) {
    for (const std::string& field : derivedOwned)
      if (assigned.count(field) != 0)
        return Invalid("'" + field + "' is derived from the mesh and the "
                       "accepted shock flux in this mode; remove it");
    const char* inputs[] = {
        "run.maximum_particle_speed_m_s", "run.time_step_margin_factor",
        "source.shock.normalization_radius_m", "source.shock.weight_model",
        "source.shock.maximum_events_per_species_per_step"};
    for (const char* field : inputs)
      if (assigned.count(field) == 0)
        return Invalid(std::string("derived particle numerics require '") +
                       field + "'");
    // The incident flux is itself the declared physical seed rate; no
    // hidden efficiency multiplies it (see particle_normalization.h).
    candidate.source.injectionEfficiency = 1.0;
  } else {
    const char* derivedOnly[] = {
        "run.maximum_particle_speed_m_s", "run.time_step_margin_factor",
        "source.shock.normalization_radius_m", "source.shock.weight_model",
        "source.shock.maximum_events_per_species_per_step"};
    for (const char* field : derivedOnly)
      if (assigned.count(field) != 0)
        return Invalid(std::string("'") + field + "' is used only with "
                       "run.time_step_mode=derived-from-mesh");
  }

  // ---- Schema-6 completeness ----------------------------------------------
  //
  // C01 distinguishes file input from typed coupled construction.  A coupled
  // host may intentionally rely on C++ defaults while assigning its supplied
  // fields.  A human-authored file must instead contain every contract block
  // and every srcSEP3D-owned value, even when the selected mode makes a value
  // inactive: a later mode edit must never revive a C++ default that did not
  // appear in the reviewed input.  The required keys are derived from the
  // shared vocabulary table (sep_input_layout.h) -- every srcSEP3D key of a
  // required block -- minus the explicitly conditional keys below, each of
  // which is checked by its own rule.
  const bool swcmeBlock = walk.physicsSwcmeLine != 0;
  std::vector<std::string> requiredSections = {
      "run", "run.restart", "application.sep", "application.sep.transport",
      "application.sep.population_control", "application.sep.storage",
      "application.sep.turbulence", "mesh", "mesh.domain", "mesh.field_line",
      "mesh.refinement.solar", "mesh.refinement.tube", "mesh.active_region",
      "mesh.memory", "physics.background", "physics.background.parker",
      "physics.shock", "source.shock", "output"};
  if (!incidentFlux) requiredSections.push_back("species");
  for (const std::string& required : requiredSections)
    if (sections.count(required) == 0)
      return Invalid("missing required block '" + required + "'");
  const std::set<std::string> conditionalKeys = {
      // Derived particle numerics (checked above): the step and the weights
      // are either written or derived, never both.
      "run.time_step_s", "run.time_step_mode", "run.maximum_particle_speed_m_s",
      "run.time_step_margin_factor", "source.shock.physical_particle_rate_per_s",
      "source.shock.injection_efficiency", "source.shock.normalization",
      "source.shock.normalization_radius_m", "source.shock.weight_model",
      "source.shock.maximum_events_per_species_per_step",
      "species.macroparticle_weight",
      // Optional stop condition (absent = step budget only).
      "run.stop_shock_radius_m",
      // Corner-cube geometry: required as a group by the corner rule below;
      // the default box is the Sun-centred cube.
      "mesh.domain.box_geometry", "mesh.domain.corner_direction_x",
      "mesh.domain.corner_direction_y", "mesh.domain.corner_direction_z",
      "mesh.domain.corner_margin_m", "mesh.field_line.end_radius_m",
      "mesh.refinement.solar.anchor", "mesh.active_region.solar_sphere_radius_m",
      // Required exactly for provider = runtime-model (typed factory rule).
      "physics.background.model_id"};
  for (const Layout::KeyRule& rule : Layout::kKeys) {
    if ((rule.support & Layout::kSrcSep3d) == 0) continue;
    const std::string block = rule.block;
    if (block == "observer.*") continue;  // per-observer rule below
    if (std::find(requiredSections.begin(), requiredSections.end(), block) ==
        requiredSections.end()) continue;
    const std::string qualified = block + "." + rule.key;
    if (conditionalKeys.count(qualified) != 0) continue;
    if (block == "physics.background.parker" &&
        Layout::SwcmeOwnedParkerKey(rule.key) != nullptr) continue;  // D4 rule
    if (block == "source.shock" &&
        std::find(Layout::SwcmeSourceKeys().begin(),
                  Layout::SwcmeSourceKeys().end(), std::string(rule.key)) !=
            Layout::SwcmeSourceKeys().end()) continue;  // D4 rule
    if (assigned.count(qualified) == 0)
      return Invalid("missing required key '" + qualified + "'");
  }
  if (!incidentFlux) {
    const char* writtenNumerics[] = {
        "run.time_step_s", "source.shock.physical_particle_rate_per_s",
        "source.shock.injection_efficiency", "species.macroparticle_weight"};
    for (const char* field : writtenNumerics)
      if (assigned.count(field) == 0)
        return Invalid(std::string("missing required key '") + field + "'");
  }

  // ---- Decision D4: one owner per physical quantity -------------------------
  //
  // With physics: swcme the canonical SWCME ambient owns the solar wind seen
  // by the analytic Parker provider, so its seven Parker counterparts must be
  // absent (they are copied from the resolved SWCME model below); without it
  // the Parker record is the only description and must be complete.  The
  // reference spectrum of source: shock feeds the SWCME source and exists only
  // together with physics: swcme.
  for (const Layout::KeyRule& rule : Layout::kKeys) {
    if (std::string(rule.block) != "physics.background.parker") continue;
    const char* owner = Layout::SwcmeOwnedParkerKey(rule.key);
    if (owner == nullptr) continue;
    const std::string qualified =
        std::string("physics.background.parker.") + rule.key;
    if (swcmeBlock && assigned.count(qualified) != 0)
      return Invalid("line " + std::to_string(walk.keyLines[qualified]) +
                     ": '" + qualified + "' is owned by physics: swcme (" +
                     owner + "); remove it (each quantity is set in one "
                     "place, src/INPUT_FILE.md section 5)");
    if (!swcmeBlock && assigned.count(qualified) == 0)
      return Invalid("missing required key '" + qualified + "' (without "
                     "physics: swcme the Parker record describes the ambient "
                     "solar wind)");
  }
  for (const std::string& key : Layout::SwcmeSourceKeys()) {
    const std::string qualified = "source.shock." + key;
    if (swcmeBlock && assigned.count(qualified) == 0)
      return Invalid("missing required key '" + qualified + "' (it defines "
                     "the canonical SWCME source spectrum)");
    if (!swcmeBlock && assigned.count(qualified) != 0)
      return Invalid("line " + std::to_string(walk.keyLines[qualified]) +
                     ": '" + qualified + "' is used only with physics: swcme");
  }
  if (swcmeBlock) {
    if (incidentFlux)
      return Invalid("source.shock.normalization=accepted-shock-incident-flux "
                     "requires the reduced shock front, not physics: swcme");
    // The SWCME source.* assignments are synthesized from source: shock and
    // carry the line of the source: shock key that supplies each of them,
    // so the canonical resolver's diagnostics still point at an editable line.
    struct SourceLine { const char* swcme; const char* sourceKey; };
    const SourceLine sourceLines[] = {
        {"source.particle_mass", "spectrum_particle_mass_kg"},
        {"source.charge_number", "spectrum_charge_number"},
        {"source.energy_min", "minimum_energy_j"},
        {"source.energy_max", "maximum_energy_j"},
        {"source.reference_energy", "spectrum_reference_energy_j"},
        {"source.injection_efficiency", "injection_efficiency"},
        {"source.normalization", "spectrum_normalization"},
        {"source.reference_intensity_si",
         "spectrum_reference_intensity_per_m2_s_sr_j"}};
    for (const Layout::SynthesizedAssignment& synthesized :
         Layout::SynthesizeSwcmeSource(walk.spectrum)) {
      SwcmeAssignment assignment;
      assignment.key = synthesized.key;
      assignment.value = synthesized.value;
      for (const SourceLine& line : sourceLines)
        if (synthesized.key == line.swcme)
          assignment.line =
              walk.keyLines[std::string("source.shock.") + line.sourceKey];
      candidate.swcmeAssignments.push_back(assignment);
    }
    // With SWCME present the Parker record expresses B_r at 1 AU; the value
    // itself is derived from the canonical total |B| below.
    candidate.parker.referenceRadiusM = Core::Const::AU;
  }
  const bool parallelSection = sections.count("physics.parallel_diffusion") != 0;
  const bool parallelSelected = candidate.spatialDiffusionModel ==
      SpatialDiffusionModel::ParallelDiffusionLibrary;
  if (parallelSelected && !parallelSection)
    return Invalid("parallel-diffusion-library requires a physics: parallel_diffusion block");
  if (!parallelSelected && parallelSection)
    return Invalid("block physics: parallel_diffusion is present but the library is not selected");
  if (parallelSection) {
    // The library parses the block body: the 'model' selector,
    // duplicate/unknown/missing keys, and every model-specific value and
    // domain.  Its diagnostics carry "line N:" from the physical line numbers
    // recorded by the walker; the file name is prefixed here.  The resulting text and line provenance are
    // copied into the options so RunConfiguration3D::Create() (which typed,
    // parser-free hosts also call) and Runtime installation reproduce the
    // identical configuration and fingerprint.
    SEP::ParallelDiffusion::ParsedSection parsed;
    const SEP::ParallelDiffusion::Status sectionStatus =
        SEP::ParallelDiffusion::ParseSection(walk.parallelDiffusionLines, &parsed);
    if (!sectionStatus.ok())
      return Invalid(walk.parallelDiffusionFile +
                     ": invalid physics: parallel_diffusion block: " +
                     sectionStatus.detail);
    candidate.parallelDiffusionModelId = parsed.modelId;
    candidate.parallelDiffusionModelLine = parsed.modelLineNumber;
    candidate.parallelDiffusionParameters.clear();
    for (const SEP::ParallelDiffusion::SectionParameter& parameter :
         parsed.parameters) {
      RunConfiguration3DOptions::ParallelDiffusionAssignment assignment;
      assignment.name = parameter.name;
      assignment.value = parameter.value;
      assignment.line = parameter.lineNumber;
      candidate.parallelDiffusionParameters.push_back(std::move(assignment));
    }
  }
  if (candidate.domainBoxGeometry == DomainBoxGeometry::FieldLineCornerCube ||
      candidate.domainBoxGeometry == DomainBoxGeometry::FieldLineXYCornerCube) {
    // Additive opt-in contract: legacy decks retain their original defaults.
    // A corner deck explicitly owns every new geometry/refinement control.
    const char* cornerFields[] = {"mesh.domain.outer_radius_mode", "mesh.domain.box_geometry",
        "mesh.domain.corner_direction_x", "mesh.domain.corner_direction_y",
        "mesh.domain.corner_direction_z", "mesh.domain.corner_margin_m",
        "mesh.field_line.end_radius_m", "mesh.field_line.length_m",
        "mesh.refinement.solar.anchor", "mesh.active_region.solar_sphere_radius_m"};
    for (const char* field : cornerFields)
      if (assigned.count(field) == 0)
        return Invalid("corner geometry is missing explicit key '" + std::string(field) + "'");
  }
  {
    const bool reducedShockFront=candidate.background==BackgroundAuthority::RuntimeModel&&
        candidate.backgroundModelId=="sep-corona-swcme-shock-front-v1";
    if (!reducedShockFront&&!swcmeBlock)
      return Invalid("a complete physics: swcme block is required unless "
                     "physics.background selects the reduced shock front");
    // Every observer block is complete: all 21 srcSEP3D observer keys are
    // present even when the selected kind makes some of them inactive.
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
          return Invalid("missing required key '" + required + "'");
      }
    }
    if (candidate.injectionCadenceSteps != 1)
      return Invalid("source.shock.samples_per_step is an exact count at "
                     "every simulation step; therefore "
                     "application.sep.injection_cadence_steps must equal 1");
    if (candidate.background == BackgroundAuthority::PythonInterpolator)
      return Core::Status::Reserved(
          "Python heliospheric-model interpolation background");
    // A standalone file selects the analytic Parker, SWCME or runtime-model
    // background with prescribed turbulence.  SWMF requires its typed host
    // rather than a standalone input alias.
    if ((candidate.background != BackgroundAuthority::AnalyticParker &&
         candidate.background != BackgroundAuthority::Swcme &&
         candidate.background != BackgroundAuthority::RuntimeModel) ||
        candidate.turbulence != TurbulenceAuthority::Prescribed)
      return Invalid("a standalone srcSEP3D input requires the analytic-parker, "
                     "swcme or runtime-model background and prescribed "
                     "turbulence; coupled SWMF hosts must use the parser-free "
                     "typed interface");

    if(!reducedShockFront) {
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
    const bool propagation = candidate.intent == RunIntent::ShockPropagation;
    if (candidate.shock != ShockAuthority::Swcme ||
        (propagation ? candidate.source.enabled
                     : (candidate.intent != RunIntent::ShockInjection || !candidate.source.enabled)))
      return Invalid("canonical SWCME requires injection with source enabled, or shock-propagation with source disabled");
    if (propagation) {
      // Launch-to-arrival telemetry must contain tick zero and have declared
      // temporal coverage through the step budget, even if a radius stop fires.
      const long double horizon = static_cast<long double>(candidate.requestedTimeStepS) * candidate.maximumTimeSteps;
      if (resolved.configuration.launch_epoch_s != 0.0 ||
          resolved.configuration.valid_from_s != 0.0 ||
          resolved.configuration.valid_until_s < horizon)
        return Invalid("shock-propagation requires event.launch_epoch=0, event.valid_from=0 and validity through the complete run horizon");
    }
    // A resolved compression is carried by the mesh. The existing DSA
    // particle source must remain off to avoid counting the same shock twice.
    // This exception is deliberately confined to the source-free propagation
    // control. Finite SSE crossing uses its own moving generating sphere and
    // cap test; this does not qualify an under-resolved compression layer.
    const bool resolvedMesh = candidate.background == BackgroundAuthority::Swcme &&
        propagation && !candidate.source.enabled &&
        model.region_mode == swcme::regions::Mode::FullICME &&
        model.shock_acceleration_mode == swcme::acceleration::Mode::ResolvedCompression;
    // The canonical library validates the SSE axis and angular width. The
    // application now preserves that geometry throughout the mover boundary.
    // Ellipsoid remains reserved until an exact ellipsoid crossing is supplied.
    if ((model.shape != swcme3d::ShockShape::Sphere &&
         model.shape != swcme3d::ShockShape::SSE) ||
        (!resolvedMesh && (model.region_mode != swcme::regions::Mode::ShockOnly ||
         model.shock_acceleration_mode != swcme::acceleration::Mode::Source)))
      return Invalid("srcSEP3D currently supports canonical SWCME only as a "
                     "Sphere/SSE SHOCK_ONLY/SOURCE, or source-free FULL_ICME/"
                     "RESOLVED_COMPRESSION mesh background; ellipsoid crossing is unsupported");
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
        return Invalid("canonical physics: swcme CME direction cannot define its "
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
        return Invalid("mesh.field_line.start_mode=cme-launch-point requires "
                       "SWCME cme.launch_radius to equal mesh.domain.inner_radius_m");
      }
      if ((candidate.parkerSpiralInitialPointM - launchPoint).Norm() >
          linkageTolerance) {
        return Invalid("mesh.field_line initial point differs from the "
                       "canonical SWCME launch-apex point");
      }
      const double sine = std::sin(candidate.tubeColatitudeRad);
      const Core::Vec3 tubeDirection(
          sine * std::cos(candidate.tubeLongitudeRad),
          sine * std::sin(candidate.tubeLongitudeRad),
          std::cos(candidate.tubeColatitudeRad));
      if ((tubeDirection - direction).Norm() > 1.0e-12) {
        return Invalid("mesh.refinement.tube source direction differs from the "
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
    // Decision D4: the SWCME ambient is the only written description of the
    // solar wind, so the former Parker-versus-SWCME agreement test reduces to
    // the one remaining cross-model constraint: the SWCME Parker source
    // (corotation) radius must be the mesh inner boundary, from which
    // RunConfiguration3D::Create derives parker.sourceRadiusM.
    const bool sameBackground =
        NearlyEqual(parkerSourceM, candidate.innerRadiusM);
    // Internal consistency assertion: the SWCME source.* values were
    // synthesized from source: shock above (J -> MeV inside the resolver), so
    // they agree to round-off; a failure here is a programming error.
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
      return Invalid("physics: swcme parker.source_radius differs from "
                     "mesh.domain.inner_radius_m");
    if (!sameSource)
      return Invalid("source.shock.spectrum_normalization must be relative-only "
                     "for srcSEP3D (the source rate is set by "
                     "physical_particle_rate_per_s and injection_efficiency)");
    if (!sameParkerGeometry)
      return Invalid("mesh.refinement.tube colatitude or analytic +Z rotation axis "
                     "differs from canonical physics: swcme Parker geometry");

    // From this point onward the typed Parker provider receives values copied
    // from the canonical SWCME resolver, not a second independently parsed
    // solar-wind model (decision D4: the Parker block holds no ambient keys
    // when physics: swcme is present).  B_r is expressed at the 1-AU
    // reference installed after the walk.
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

    if (candidate.source.enabled) {
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
                       "source.shock.samples_per_step");
    }
    // Keep the legacy typed diagnostic consistent with the sole canonical
    // launch authority; radius-stop validation must use the actual launch.
    if (propagation)
      candidate.shockModel.initialRadiusM = model.r0_Rs * swcme::constants::SOLAR_RADIUS_M;
    candidate.swcmeConfigurationFingerprint =
        resolved.configuration.fingerprint;
    candidate.swcmeResolvedManifest =
        resolved.configuration.normalized_manifest;
    } else {
      // The reduced provider owns its finite front; the legacy particle-facing
      // shock authority stays `none` in every reduced mode.  Two modes exist:
      //  * zero-particle (transport-only or shock-propagation): the source is
      //    disabled and the front is a background/diagnostic authority;
      //  * reduced-front injection (shock-injection): particles are injected
      //    from accepted triangles of the same front with the derived particle
      //    numerics above (AGENTS.md "Authorized srcSEP3D reduced-front source
      //    extension").  This was formerly reachable only through the retired
      //    one-dash section input; it is now an explicit, validated mode.
      if(walk.shockFrontBlock==nullptr)
        return Invalid("physics.background.model_id=sep-corona-swcme-shock-front-v1 "
                       "requires a physics: shock_front block (usually "
                       "#included from its event file)");
      if(candidate.shock!=ShockAuthority::None)
        return Invalid("reduced shock-front runtime requires "
                       "physics.shock.authority=none");
      if(incidentFlux) {
        if(candidate.intent!=RunIntent::ShockInjection||
            !candidate.source.enabled)
          return Invalid("reduced-front injection requires "
                         "application.sep.intent=shock-injection and "
                         "source.shock.enabled=true");
      } else {
        const bool reducedIntent=candidate.intent==RunIntent::TransportOnly||
            candidate.intent==RunIntent::ShockPropagation;
        if(!reducedIntent||candidate.source.enabled)
          return Invalid("zero-particle reduced shock-front runtime requires "
                         "transport-only or shock-propagation and "
                         "source.shock.enabled=false (select "
                         "source.shock.normalization=accepted-shock-incident-flux "
                         "for reduced-front injection)");
      }
      if(!candidate.swcmeAssignments.empty())
        return Invalid("reduced shock-front runtime rejects inactive physics: swcme assignments");
      if(candidate.coordinateFrame!="HCI")
        return Invalid("reduced shock-front runtime requires mesh.domain.coordinate_frame=HCI");
      // The shared resolver receives the canonical block text; relative assets
      // (assets.harmonics_file) resolve from the file that holds the block.
      candidate.backgroundModelAssetPath.clear();
      candidate.backgroundModelInlineConfiguration=
          ShockFrontCanonicalText(*walk.shockFrontBlock);
      std::filesystem::path assetDirectory=
          std::filesystem::path(walk.shockFrontBlock->begin.file).parent_path();
      if(assetDirectory.empty())assetDirectory=".";
      std::error_code error;
      const std::filesystem::path resolvedDirectory=
          std::filesystem::weakly_canonical(assetDirectory,error);
      candidate.backgroundModelAssetDirectory=
          error?assetDirectory.string():resolvedDirectory.string();
    }
  }
  if(walk.shockFrontBlock!=nullptr&&
      !(candidate.background==BackgroundAuthority::RuntimeModel&&
        candidate.backgroundModelId=="sep-corona-swcme-shock-front-v1"))
    return Invalid(AMPS::Input::ErrorAt(walk.shockFrontBlock->begin,
        "physics: shock_front is present but physics.background does not "
        "select provider=runtime-model with model_id="
        "sep-corona-swcme-shock-front-v1").message);
  if(incidentFlux&&walk.shockFrontBlock==nullptr)
    return Invalid("source.shock.normalization=accepted-shock-incident-flux "
                   "requires the reduced shock front (physics: shock_front)");
  bool observerSectionSeen = false;
  for (const std::string& present : sections) {
    if (present.rfind("observer.", 0) == 0) observerSectionSeen = true;
  }
  if (!observerSectionSeen && candidate.intent == RunIntent::ShockPropagation)
    candidate.observers.clear();  // No energetic-particle products in an empty run.
  else if (!observerSectionSeen)
    return Invalid("at least one '#begin observer: ID' block is required");
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

Core::Status ParseConfigurationText(const std::string& text,
                                   RunConfiguration3DOptions* result,
                                   const std::string& baseDirectory) {
  // In-memory input (tests, generated decks).  The display name places the
  // text in baseDirectory so that relative #include paths and the assets of
  // an included physics: shock_front block resolve from there.
  const std::string displayName =
      (std::filesystem::path(baseDirectory) / "<srcSEP3D input text>").string();
  AMPS::Input::Document document;
  const AMPS::Input::Status parsed =
      AMPS::Input::ParseText(text, displayName, baseDirectory, &document);
  if (!parsed.ok()) return Invalid(parsed.message);
  return ParseConfigurationDocument(document, result);
}

Core::Status LoadConfigurationFile(
    const std::string& path, RunConfiguration3DOptions* result) {
  if (path.empty()) return Invalid("configuration path is empty");
  AMPS::Input::Document document;
  const AMPS::Input::Status parsed = AMPS::Input::ParseFile(path, &document);
  if (!parsed.ok()) return Invalid(parsed.message);
  return ParseConfigurationDocument(document, result);
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
  if (candidate.commandLine.listTests) {
    *result = candidate;
    return Core::Status::OK();
  }
  const bool nativeTestRun = candidate.commandLine.allTests ||
      !candidate.commandLine.testSuite.empty() ||
      !candidate.commandLine.tests.empty();
  RunConfiguration3DOptions options;
  status = LoadConfigurationFile(
      nativeTestRun ? candidate.commandLine.testInputPath
                    : candidate.commandLine.inputPath,
      &options);
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
  if (options.background==BackgroundAuthority::Swcme ||
      options.background==BackgroundAuthority::RuntimeModel) {
    // Resolve/prepare the selected runtime before expensive AMR allocation.
    // This temporary provider proves launch-time metadata/coverage only; the
    // production owner-cell provider is created/installed at acquisition and
    // must pass the collective fill gate with its actual samples again.
    std::shared_ptr<Background::BackgroundProvider> provider;
    status=CreateBackgroundProvider(*candidate.configuration,&provider);
    if (status.ok())status=provider->Prepare(0.0);
    if (status.ok() && (!provider->PreparedMetadata() ||
        !BackgroundAuthorityMatches(options.background,provider->PreparedMetadata()->provider) ||
        provider->PreparedMetadata()->coordinateFrame!=options.coordinateFrame))
      return Invalid("runtime background factory returned incompatible provider metadata");
    if (!status.ok())return status;
  }
  if (nativeTestRun &&
      candidate.commandLine.testSteps > options.maximumTimeSteps)
    return Invalid("--test-steps exceeds run.maximum_time_steps in the "
                   "immutable test input");
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
  resolution.solarRefinementAnchor = options.solarRefinementAnchor;
  resolution.activeSolarSphereRadiusM = options.activeSolarSphereRadiusM;
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
  resolution.activeRegion = options.activeRegion;
  resolution.activeTubeReferenceRadiusM = options.activeTubeReferenceRadiusM;
  resolution.activeTubeRadiusAtReferenceM =
      options.activeTubeRadiusAtReferenceM;
  resolution.activeTubeRadiusMode = options.activeTubeRadiusMode;
  resolution.activeTubeBufferBlocks = options.activeTubeBufferBlocks;
  resolution.solarWindSpeedMPerS = options.parker.solarWindSpeedMPerS;
  resolution.solarRotationRateRadPerS = options.parker.solarRotationRateRadPerS;
  resolution.rotationAxis = options.parker.rotationAxis;
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
         << "run_intent=" << Name(options.intent) << '\n'
         << "source_enabled=" << (options.source.enabled ? "true" : "false") << '\n'
         << "stop_shock_radius_m=" << options.stopShockRadiusM << '\n'
         << "maximum_time_steps=" << options.maximumTimeSteps << '\n'
         << "physics_fingerprint=" << configuration.physics_fingerprint() << '\n'
         << "domain_preset=" << Name(options.domain) << '\n'
         << "domain_box_geometry=" << Name(options.domainBoxGeometry) << '\n'
         << "domain_min_m=" << domain.minimumM.x << ',' << domain.minimumM.y << ',' << domain.minimumM.z << '\n'
         << "domain_max_m=" << domain.maximumM.x << ',' << domain.maximumM.y << ',' << domain.maximumM.z << '\n'
         << "domain_side_m=" << domain.maximumM.x - domain.minimumM.x << '\n'
         << "solar_sphere_radius_m=" << options.activeSolarSphereRadiusM << '\n'
         << "solar_refinement_anchor=" << Name(options.solarRefinementAnchor) << '\n'
         << "inner_radius_m=" << options.innerRadiusM << '\n'
         << "outer_radius_m=" << options.outerRadiusM << '\n'
         << "background_provider=" << Name(options.background) << '\n'
         // Distinguish an evolving regional mesh source from a Parker ambient
         // source even when both use SWCME for shock geometry.
         << "solar_wind_model=" << (options.background==BackgroundAuthority::Swcme
             ? "canonical-swcme3d-regional-background" : "canonical-swcme-parker-leblanc") << '\n'
         << "background_model_id=" << options.backgroundModelId << '\n' 
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
         << "parker_spiral_end_radius_m=" << options.parkerSpiralEndRadiusM << '\n'
         << "parker_spiral_end_m=" << lineEnd.x << ',' << lineEnd.y << ','
         << lineEnd.z << '\n'
         << "compiled_species_authority=AMPS-SpeciesList\n"
         << "time_step_s=" << options.requestedTimeStepS << '\n'
         << "base_particle_weight=" << options.species.macroparticleWeight << '\n'
         << "particle_mover=" << Name(options.particleMover) << '\n'
         << "transport_model=" << Name(options.transport) << '\n'
         << "spatial_diffusion_model="
         << Name(options.spatialDiffusionModel) << '\n'
         << "pitch_angle_diffusion_model="
         << Name(options.pitchAngleDiffusionModel) << '\n'
         << "mean_free_path_model=" << Name(options.meanFreePathModel) << '\n'
         << "mean_free_path_reference_m="
         << options.meanFreePathReferenceM << '\n'
         << "mean_free_path_reference_radius_m="
         << options.meanFreePathReferenceRadiusM << '\n'
         << "mean_free_path_reference_rigidity_v="
         << options.meanFreePathReferenceRigidityV << '\n'
         << "mean_free_path_radial_exponent="
         << options.meanFreePathRadialExponent << '\n'
         << "mean_free_path_rigidity_exponent="
         << options.meanFreePathRigidityExponent << '\n'
         << "focused_scattering_frame="
         << Name(options.focusedScatteringFrame) << '\n'
         << "active_region_mode=" << Name(options.activeRegion) << '\n'
         << "active_tube_reference_radius_m="
         << options.activeTubeReferenceRadiusM << '\n'
         << "active_tube_radius_at_reference_m="
         << options.activeTubeRadiusAtReferenceM << '\n'
         << "active_tube_radius_mode="
         << Name(options.activeTubeRadiusMode) << '\n'
         << "active_tube_buffer_blocks="
         << options.activeTubeBufferBlocks << '\n'
         << "population_control_mode="
         << Name(options.populationControl) << '\n'
         << "population_control_minimum_per_cell_species="
         << options.minimumParticlesPerCellPerSpecies << '\n'
         << "population_control_target_per_cell_species="
         << options.targetParticlesPerCellPerSpecies << '\n'
         << "population_control_maximum_per_cell_species="
         << options.maximumParticlesPerCellPerSpecies << '\n'
         << "population_control_cadence_steps="
         << options.populationControlCadenceSteps << '\n'
         << "observer_count=" << options.observers.size() << '\n'
         << "initialization_data_tecplot_base="
         << options.initializationDataTecplotFile << '\n'
         << "source_samples_per_compiled_species="
         << options.source.samplesPerStep << '\n'
         << "source_spectrum_model="
         << Name(options.source.spectrumModel) << '\n'
         << "source_phase_space_power_index="
         << options.source.fixedPhaseSpacePowerIndex << '\n'
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
