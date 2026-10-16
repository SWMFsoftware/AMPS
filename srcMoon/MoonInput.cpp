#include "MoonInput.h"

// Use the repository-wide runtime-input grammar.  This header owns only the
// syntax tree (#begin/#end/#include, duplicate detection, and source
// locations); srcMoon remains the owner of the application: moon keys and of
// every physical/domain validation performed below.
#include "../src/general/input/amps_input.h"

#include <algorithm>
#include <array>
#include <cerrno>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <iterator>
#include <limits>
#include <map>
#include <sstream>
#include <utility>

// Runtime parser for the application-owned portion of a shared AMPS input
// file.  Parsing happens before PIC/MPI initialization, performs no global
// mutation, and commits the completed Configuration only after every syntax,
// value, and path check succeeds.
namespace Moon {
namespace Runtime {
namespace {

namespace fs = std::filesystem;

// The executable has one immutable application configuration.  These objects
// are written only by InstallConfiguration(), before amps_init(), and are read
// thereafter by the production initialization path.
bool gInstalled = false;
Configuration gConfiguration;

struct Setting {
  // Preserve the unparsed value and its complete source location.  Generic
  // input files may obtain assignments through #include, so retaining only a
  // line number from the root file would point owner-level diagnostics at the
  // wrong file.
  std::string value;
  std::string file;
  std::size_t line = 0;
};

// Input syntax is ASCII-oriented, but the unsigned-char conversion avoids the
// undefined behavior of passing a negative signed char to std::isspace().
std::string Trim(const std::string& text) {
  const auto first = std::find_if_not(text.begin(), text.end(),
      [](unsigned char c) { return std::isspace(c) != 0; });
  if (first == text.end()) return std::string();
  const auto last = std::find_if_not(text.rbegin(), text.rend(),
      [](unsigned char c) { return std::isspace(c) != 0; }).base();
  return std::string(first, last);
}

std::string Lower(std::string text) {
  std::transform(text.begin(), text.end(), text.begin(),
      [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
  return text;
}

// Quotes are optional around scalar/path values.  Deliberately do not perform
// shell expansion, environment substitution, or escape processing: the run
// receipt must identify exactly the path supplied by the user.
std::string Unquote(const std::string& value) {
  if (value.size() >= 2 &&
      ((value.front() == '"' && value.back() == '"') ||
       (value.front() == '\'' && value.back() == '\''))) {
    return value.substr(1, value.size() - 2);
  }
  return value;
}

// Reject partial parses (for example "100 km"), overflow, NaN/infinity, zero,
// and negative lengths.  surface_mesh_resolution_m is an SI length.
bool ParsePositiveFinite(const std::string& text, double* result) {
  const std::string scalar = Unquote(Trim(text));
  errno = 0;
  char* end = nullptr;
  const double value = std::strtod(scalar.c_str(), &end);
  if (errno == ERANGE || end == scalar.c_str() || *end != '\0' ||
      !std::isfinite(value) || value <= 0.0) return false;
  *result = value;
  return true;
}

// Source-rate receipts must never interpret a negative number as a physical
// sink. Keep nonnegative parsing separate from the generic positive parser;
// the model-specific block below applies the stricter positive-rate rule.
bool ParseNonnegativeFinite(const std::string& text, double* result) {
  const std::string scalar = Unquote(Trim(text));
  errno = 0;
  char* end = nullptr;
  const double value = std::strtod(scalar.c_str(), &end);
  if (errno == ERANGE || end == scalar.c_str() || *end != '\0' ||
      !std::isfinite(value) || value < 0.0) return false;
  *result = value;
  return true;
}

// Counts such as logarithmic points per decade must be exact positive base-10
// integers.  Rejecting a floating representation ("8.0") keeps receipts
// canonical and prevents an unnoticed truncation from changing the grid.
bool ParsePositiveInteger(const std::string& text, int* result) {
  const std::string scalar = Unquote(Trim(text));
  errno = 0;
  char* end = nullptr;
  const long value = std::strtol(scalar.c_str(), &end, 10);
  if (errno == ERANGE || end == scalar.c_str() || *end != '\0' ||
      value <= 0 || value > std::numeric_limits<int>::max()) return false;
  *result = static_cast<int>(value);
  return true;
}

bool ParseBoolean(const std::string& text, bool* result) {
  const std::string scalar = Lower(Unquote(Trim(text)));
  if (scalar == "true") {
    *result = true;
    return true;
  }
  if (scalar == "false") {
    *result = false;
    return true;
  }
  return false;
}

// Parse a comma-separated, strictly increasing list of SZA values.  Requiring
// explicit values in [0,180] degrees avoids an implicit angular convention and
// makes duplicated profiles a startup error instead of wasted campaign work.
bool ParseSolarZenithAngles(const std::string& text,
                            std::vector<double>* result) {
  std::stringstream stream(Unquote(Trim(text)));
  std::vector<double> values;
  std::string field;
  while (std::getline(stream, field, ',')) {
    double value = 0.0;
    const std::string scalar = Trim(field);
    if (scalar.empty()) return false;

    errno = 0;
    char* end = nullptr;
    value = std::strtod(scalar.c_str(), &end);
    if (errno == ERANGE || end == scalar.c_str() || *end != '\0' ||
        !std::isfinite(value) || value < 0.0 || value > 180.0 ||
        (!values.empty() && value <= values.back())) return false;
    values.push_back(value);
  }
  if (values.empty()) return false;
  *result = std::move(values);
  return true;
}

// Generic comma-separated finite scalar list used for the five heat-capacity
// polynomial coefficients.  Unlike a spatial grid, coefficients may be zero
// or negative, so this parser deliberately imposes no sign/order policy.
bool ParseFiniteList(const std::string& text, std::vector<double>* result) {
  std::stringstream stream(Unquote(Trim(text)));
  std::vector<double> values;
  std::string field;
  while (std::getline(stream,field,',')) {
    const std::string scalar=Trim(field);
    errno=0;
    char* end=nullptr;
    const double value=std::strtod(scalar.c_str(),&end);
    if (errno==ERANGE || end==scalar.c_str() || *end!='\0' ||
        !std::isfinite(value)) return false;
    values.push_back(value);
  }
  if (values.empty()) return false;
  *result=std::move(values);
  return true;
}

// Parse a temperature:probability table.  Temperatures are strictly
// increasing and positive; probabilities are bounded.  Keeping pairs in one
// field prevents independently edited lists from assigning a probability to
// the wrong laboratory temperature.
bool ParseProbabilityTable(const std::string& text,
                           std::vector<double>* temperatureK,
                           std::vector<double>* probability) {
  std::stringstream stream(Unquote(Trim(text)));
  std::vector<double> temperatures,probabilities;
  std::string field;
  while (std::getline(stream,field,',')) {
    const std::size_t colon=field.find(':');
    if (colon==std::string::npos || field.find(':',colon+1)!=std::string::npos)
      return false;
    double temperature=0.0,value=0.0;
    if (!ParsePositiveFinite(field.substr(0,colon),&temperature)) return false;
    const std::string valueText=Trim(field.substr(colon+1));
    errno=0;
    char* end=nullptr;
    value=std::strtod(valueText.c_str(),&end);
    if (errno==ERANGE || end==valueText.c_str() || *end!='\0' ||
        !std::isfinite(value) || value<0.0 || value>1.0 ||
        (!temperatures.empty() && temperature<=temperatures.back()))
      return false;
    temperatures.push_back(temperature);
    probabilities.push_back(value);
  }
  if (temperatures.size()<2) return false;
  *temperatureK=std::move(temperatures);
  *probability=std::move(probabilities);
  return true;
}

// Parse the M1 angular-refinement receipt as ``SZA_DEG:RESOLUTION_M`` pairs.
// Keeping angles and lengths in one ordered field prevents a truncated or
// reordered pair of parallel lists from silently assigning a resolution to
// the wrong solar-zenith angle.  Endpoints and the terminator are mandatory
// because I12 evaluates those physical locations explicitly.
bool ParseAngularResolutionTable(const std::string& text,
                                 std::vector<double>* anglesDeg,
                                 std::vector<double>* resolutionsM) {
  std::stringstream stream(Unquote(Trim(text)));
  std::vector<double> parsedAngles, parsedResolutions;
  std::string field;
  while (std::getline(stream, field, ',')) {
    const std::string pair = Trim(field);
    const std::size_t colon = pair.find(':');
    if (colon == std::string::npos || pair.find(':', colon + 1) !=
        std::string::npos) return false;

    const std::string angleText = Trim(pair.substr(0, colon));
    const std::string resolutionText = Trim(pair.substr(colon + 1));
    errno = 0;
    char* end = nullptr;
    const double angle = std::strtod(angleText.c_str(), &end);
    if (errno == ERANGE || end == angleText.c_str() || *end != '\0' ||
        !std::isfinite(angle) || angle < 0.0 || angle > 180.0 ||
        (!parsedAngles.empty() && angle <= parsedAngles.back())) return false;
    double resolution = 0.0;
    if (!ParsePositiveFinite(resolutionText, &resolution)) return false;
    parsedAngles.push_back(angle);
    parsedResolutions.push_back(resolution);
  }

  const auto terminator = std::find(parsedAngles.begin(), parsedAngles.end(),
                                    90.0);
  if (parsedAngles.size() < 3 || parsedAngles.front() != 0.0 ||
      parsedAngles.back() != 180.0 || terminator == parsedAngles.end()) {
    return false;
  }
  *anglesDeg = std::move(parsedAngles);
  *resolutionsM = std::move(parsedResolutions);
  return true;
}

bool IsIntegralNumberOfSteps(double span, double step) {
  const double count = span / step;
  return std::abs(count - std::round(count)) <=
      32.0 * std::numeric_limits<double>::epsilon() *
          std::max(1.0, std::abs(count));
}

// Construct the exact shared comparison grid in SI units.  Linear endpoints
// are required to be integral step counts, so they cannot be reached by an
// undocumented short final bin.  The far region advances by 10^(1/N); when an
// adjustable far endpoint is not itself on that logarithmic sequence, it is
// appended explicitly and therefore remains represented without extrapolation.
bool BuildProfileAltitudeGrid(Configuration* configuration,
                              std::string* reason) {
  if (configuration->profileAltitudeNearMaximumM >=
          configuration->profileAltitudeMiddleMaximumM ||
      configuration->profileAltitudeMiddleMaximumM >=
          configuration->profileAltitudeFarMaximumM) {
    *reason = "profile altitude maxima must satisfy near < middle < far";
    return false;
  }
  if (!IsIntegralNumberOfSteps(configuration->profileAltitudeNearMaximumM,
                               configuration->profileAltitudeNearStepM)) {
    *reason = "profile_altitude_near_maximum_m must be an integral multiple "
              "of profile_altitude_near_step_m";
    return false;
  }
  if (!IsIntegralNumberOfSteps(
          configuration->profileAltitudeMiddleMaximumM -
              configuration->profileAltitudeNearMaximumM,
          configuration->profileAltitudeMiddleStepM)) {
    *reason = "middle-minus-near altitude span must be an integral multiple "
              "of profile_altitude_middle_step_m";
    return false;
  }

  const long double nearPointCount =
      configuration->profileAltitudeNearMaximumM /
          configuration->profileAltitudeNearStepM + 1.0;
  const long double middlePointCount =
      (configuration->profileAltitudeMiddleMaximumM -
       configuration->profileAltitudeNearMaximumM) /
          configuration->profileAltitudeMiddleStepM;
  const long double farPointCount = std::ceil(
      std::log10(configuration->profileAltitudeFarMaximumM /
                 configuration->profileAltitudeMiddleMaximumM) *
      configuration->profileAltitudeFarPointsPerDecade);

  // M1 deliberately imposes no policy limit on profile points, wall time, or
  // memory.  This check is only a representation guard: the loop counters and
  // std::vector must be able to name the requested elements on this build.
  // A representable request is attempted in full; if the host cannot provide
  // the memory, main()'s new-handler terminates through AMPS exit() with the
  // required resource diagnosis and guidance.
  const long double totalPointCount = std::ceil(
      nearPointCount + middlePointCount + farPointCount);
  const long double representablePointCount = std::min(
      static_cast<long double>(std::numeric_limits<long>::max()),
      static_cast<long double>(std::vector<double>().max_size()));
  if (!std::isfinite(totalPointCount) ||
      totalPointCount > representablePointCount) {
    *reason = "profile altitude grid cannot be represented by this executable";
    return false;
  }

  std::vector<double> grid;
  // Reserve once so genuine resource exhaustion is diagnosed at the point
  // where the requested input grid is materialized, rather than after a long
  // sequence of reallocations with an ambiguous final failure site.
  grid.reserve(static_cast<std::size_t>(totalPointCount));
  const long nearSteps = std::lround(
      configuration->profileAltitudeNearMaximumM /
      configuration->profileAltitudeNearStepM);
  for (long i = 0; i <= nearSteps; ++i) {
    grid.push_back(i * configuration->profileAltitudeNearStepM);
  }

  const long middleSteps = std::lround(
      (configuration->profileAltitudeMiddleMaximumM -
       configuration->profileAltitudeNearMaximumM) /
      configuration->profileAltitudeMiddleStepM);
  for (long i = 1; i <= middleSteps; ++i) {
    grid.push_back(configuration->profileAltitudeNearMaximumM +
                   i * configuration->profileAltitudeMiddleStepM);
  }

  const double ratio = std::pow(10.0,
      1.0 / configuration->profileAltitudeFarPointsPerDecade);
  double altitude = configuration->profileAltitudeMiddleMaximumM;
  while (altitude * ratio < configuration->profileAltitudeFarMaximumM *
      (1.0 - 32.0 * std::numeric_limits<double>::epsilon())) {
    altitude *= ratio;
    grid.push_back(altitude);
  }
  if (grid.back() < configuration->profileAltitudeFarMaximumM) {
    grid.push_back(configuration->profileAltitudeFarMaximumM);
  }

  configuration->profileAltitudeGridM = std::move(grid);
  reason->clear();
  return true;
}

// Give all parser diagnostics a uniform file:line prefix.  Line zero denotes
// a whole-file error such as a missing required section or key.
std::string At(const fs::path& file, std::size_t line,
               const std::string& message) {
  std::ostringstream out;
  out << file.string();
  if (line != 0) out << ':' << line;
  out << ": srcMoon input error: " << message;
  return out.str();
}

// Owner-level value diagnostics use the assignment's native source file,
// including when the generic parser expanded that assignment from #include.
std::string At(const Setting& setting, const std::string& message) {
  return At(fs::path(setting.file), setting.line, message);
}

// Resolve relative paths against the input file—not the launch directory—so
// rerunning the same receipt from a different output directory is invariant.
bool AbsolutePath(const Setting& setting, const std::string& key,
                  fs::path* output,
                  std::string* error) {
  fs::path value(Unquote(Trim(setting.value)));
  if (value.empty()) {
    *error = At(setting, key + " cannot be empty");
    return false;
  }
  // Match the generic #include contract: a relative value is anchored at the
  // file containing that assignment, not at the top-level receipt or launch
  // directory.
  const fs::path inputDirectory = fs::path(setting.file).parent_path();
  if (value.is_relative()) value = inputDirectory / value;
  std::error_code ec;
  value = fs::absolute(value, ec).lexically_normal();
  if (ec) {
    *error = At(setting,
        "cannot resolve " + key + ": " + ec.message());
    return false;
  }
  *output = value;
  return true;
}

}  // namespace

const char* SurfaceGeometryName(SurfaceGeometry geometry) {
  return geometry == SurfaceGeometry::Sphere ? "sphere" : "lola";
}

bool ParseApplicationInput(const std::string& path, Configuration* result,
                           std::string* error) {
  if (result == nullptr || error == nullptr) return false;
  error->clear();

  std::error_code ec;
  fs::path inputFile = fs::absolute(fs::path(path), ec).lexically_normal();
  if (ec) {
    *error = path + ": srcMoon input error: cannot resolve input path: " +
        ec.message();
    return false;
  }

  // Parse the complete document with the common AMPS grammar.  The generic
  // reader validates structure, include balance, duplicate keys/blocks, and
  // source locations.  The Moon owner then requires one flat
  // ``application: moon`` block; accepting an unowned root block would violate
  // the repository rule that unknown configuration is fatal rather than
  // silently ignored.
  AMPS::Input::Document document;
  AMPS::Input::Status status =
      AMPS::Input::ParseFile(inputFile.string(), &document);
  if (!status.ok()) {
    *error = status.message;
    return false;
  }
  status = AMPS::Input::RequireNoAssignments(document.root);
  if (!status.ok()) {
    *error = status.message;
    return false;
  }
  status = AMPS::Input::RequireOnlyBlocks(
      document.root, {"application"}, "the srcMoon executable");
  if (!status.ok()) {
    *error = status.message;
    return false;
  }

  const AMPS::Input::Block* moon = nullptr;
  for (const AMPS::Input::Block& block : document.root.blocks) {
    if (block.name != "moon") {
      *error = AMPS::Input::ErrorAt(
          block.begin, "application block '" + block.Label() +
          "' is not handled by the srcMoon executable").message;
      return false;
    }
    moon = &block;
  }
  if (moon == nullptr) {
    *error = At(inputFile, 0,
        "required '#begin application: moon' block was not found");
    return false;
  }
  status = AMPS::Input::RequireOnlyBlocks(
      *moon, {}, "the srcMoon application: moon owner");
  if (!status.ok()) {
    *error = status.message;
    return false;
  }

  // AMPS generic keys are case-sensitive.  The Moon vocabulary deliberately
  // uses the exact lower-case spellings documented in src/INPUT_FILE.md, so a
  // differently cased key is diagnosed as unknown rather than normalized.
  std::map<std::string, Setting> settings;
  for (const AMPS::Input::Assignment& assignment : moon->assignments) {
    settings.emplace(assignment.key,
        Setting{assignment.value, assignment.where.file,
                assignment.where.line});
  }

  // Require a complete run receipt in both sphere and LOLA modes.  This makes
  // changing surface_geometry explicit and prevents a later LOLA run from
  // inheriting undeclared paths from compiled defaults or process state.
  const char* required[] = {
      "spice_path", "surface_geometry", "surface_mesh_resolution_m",
      "mesh_profile", "mesh_near_surface_cell_size_m",
      "mesh_near_surface_extent_m", "mesh_far_field_cell_size_m",
      "mesh_far_field_extent_m", "mesh_outer_cell_size_m",
      "mesh_output_tecplot",
      "mesh_angular_resolution_m_by_sza",
      "mesh_terminator_half_width_deg",
      "mesh_polar_cap_start_latitude_deg", "mesh_polar_resolution_m",
      "lola_product_id", "lola_image_file", "lola_label_file",
      "surface_cea_file", "surface_tecplot_file",
      "m2_environment",
      "simulation_start_utc", "simulation_environment_mode",
      "random_seed",
      "profile_altitude_near_maximum_m", "profile_altitude_near_step_m",
      "profile_altitude_middle_maximum_m", "profile_altitude_middle_step_m",
      "profile_altitude_far_maximum_m",
      "profile_altitude_far_points_per_decade",
      "profile_solar_zenith_angles_deg",
      "profile_include_body_fixed_poles", "profile_include_shell_average",
      "profile_shell_angular_weighting",
      "limb_observer", "limb_plane", "limb_sides",
      "limb_branch_reference",
      "limb_tangent_maximum_altitude_m",
      "tail_center_direction", "tail_map_coordinate_system",
      "tail_map_display_orientation", "tail_map_half_width_deg",
      "tail_beam_diameter_deg", "tail_exposure_duration_s",
      "tail_beam_response", "tail_temporal_operator"};
  const char* m2Required[] = {
      "surface_temperature_mode", "diviner_temperature_table_file",
      "diviner_eclipse_history_file",
      "diviner_eclipse_selection_contract_file",
      "illumination_ray_origin_offset_m",
      "surface_collision_maximum_segment_length_m",
      "surface_collision_time_tolerance_s",
      "surface_collision_distance_tolerance_m",
      "surface_collision_maximum_iterations",
      "thermal_solar_constant_w_m2",
      "thermal_infrared_emissivity", "thermal_normal_bond_albedo",
      "thermal_albedo_coefficient_a", "thermal_albedo_coefficient_b",
      "thermal_surface_density_kg_m3", "thermal_deep_density_kg_m3",
      "thermal_density_scale_height_m",
      "thermal_surface_conductivity_w_m_k",
      "thermal_deep_conductivity_w_m_k",
      "thermal_radiative_conductivity_parameter",
      "thermal_heat_capacity_coefficients_j_kg_k",
      "thermal_geothermal_flux_w_m2", "thermal_bottom_depth_m",
      "thermal_layer_count", "thermal_time_step_s",
      "thermal_initial_temperature_k", "thermal_nonlinear_tolerance_k",
      "thermal_nonlinear_maximum_iterations",
      "thermal_eclipse_direct_irradiance_model",
      "thermal_solar_disc_radiance_model",
      "thermal_earth_atmosphere_transmission_model",
      "thermal_solar_disc_radius_m", "thermal_opaque_earth_radius_m",
      "thermal_spinup_cycle_start_utc", "thermal_spinup_cycle_end_utc",
      "thermal_spinup_temperature_tolerance_k",
      "thermal_spinup_energy_tolerance_j_m2",
      "sodium_activation_energy_ev", "sodium_attempt_frequency_hz",
      "sodium_sticking_model",
      "sodium_sticking_probability_by_temperature_k",
      "sodium_accommodation_speed_model",
      "sodium_accommodation_coefficient",
      "sodium_surface_interaction_reference",
      "argon_desorption_energy_j_mol", "argon_sticking_model",
      "argon_sticking_probability_by_temperature_k",
      "argon_accommodation_speed_model",
      "argon_accommodation_coefficient",
      "argon_surface_interaction_reference",
      "argon_source_enabled",
      "argon_source_spatial_model", "argon_source_rate_atoms_s",
      "argon_source_time_model", "argon_source_map_file",
      "argon_source_reference",
      "ion_surface_interaction_mode", "cold_trap_retention_timescale_s",
      "surface_inventory_initial_state",
      "surface_inventory_checkpoint_file",
      "surface_inventory_checkpoint_interval_s"};
  for (const char* key : required) {
    if (settings.find(key) == settings.end()) {
      *error = At(inputFile, 0,
                  std::string("application: moon is missing required ") + key);
      return false;
    }
  }
  for (const auto& item : settings) {
    if (std::find(std::begin(required), std::end(required), item.first) ==
            std::end(required) &&
        std::find(std::begin(m2Required), std::end(m2Required), item.first) ==
            std::end(m2Required)) {
      *error = At(item.second,
                  "unrecognized application: moon setting '" +
                  item.first + "'");
      return false;
    }
  }

  // Build a local candidate and assign *result only at the end.  Callers can
  // safely retain their previous object if any validation below fails.
  Configuration candidate;
  candidate.inputFile = inputFile.string();
  const std::string m2Mode=Lower(Unquote(Trim(
      settings["m2_environment"].value)));
  if (m2Mode!="disabled" && m2Mode!="enabled") {
    *error=At(settings["m2_environment"],
        "m2_environment must be disabled or enabled");
    return false;
  }
  candidate.m2EnvironmentEnabled=m2Mode=="enabled";
  if (candidate.m2EnvironmentEnabled) {
    for (const char* key:m2Required) if (settings.find(key)==settings.end()) {
      *error=At(inputFile,0,
          std::string("enabled M2 environment is missing required ")+key);
      return false;
    }
  }
  else {
    // A disabled M2 receipt must not carry dormant physics values.  Rejecting
    // partial groups prevents users from believing that an ignored parameter
    // was active in an M0/M1 run.
    for (const char* key:m2Required) if (settings.find(key)!=settings.end()) {
      *error=At(settings[key],std::string(key)+
          " is not allowed while m2_environment=disabled");
      return false;
    }
  }
  const std::string geometry = Lower(Unquote(settings["surface_geometry"].value));
  if (geometry == "sphere") candidate.surfaceGeometry = SurfaceGeometry::Sphere;
  else if (geometry == "lola") candidate.surfaceGeometry = SurfaceGeometry::Lola;
  else {
    *error = At(settings["surface_geometry"],
                "surface_geometry must be sphere or lola");
    return false;
  }

  if (!ParsePositiveFinite(settings["surface_mesh_resolution_m"].value,
                           &candidate.surfaceMeshResolutionM)) {
    *error = At(settings["surface_mesh_resolution_m"],
                "surface_mesh_resolution_m must be finite and positive");
    return false;
  }

  candidate.meshProfile = Lower(Unquote(Trim(
      settings["mesh_profile"].value)));
  if (candidate.meshProfile != "quick" &&
      candidate.meshProfile != "full" &&
      candidate.meshProfile != "custom") {
    *error = At(settings["mesh_profile"],
        "mesh_profile must be quick, full, or custom");
    return false;
  }
  const char* meshLengthKeys[] = {
      "mesh_near_surface_cell_size_m", "mesh_near_surface_extent_m",
      "mesh_far_field_cell_size_m", "mesh_far_field_extent_m",
      "mesh_outer_cell_size_m", "mesh_terminator_half_width_deg",
      "mesh_polar_cap_start_latitude_deg", "mesh_polar_resolution_m"};
  double* meshLengthTargets[] = {
      &candidate.meshNearSurfaceCellSizeM,
      &candidate.meshNearSurfaceExtentM,
      &candidate.meshFarFieldCellSizeM,
      &candidate.meshFarFieldExtentM,
      &candidate.meshOuterCellSizeM,
      &candidate.meshTerminatorHalfWidthDeg,
      &candidate.meshPolarCapStartLatitudeDeg,
      &candidate.meshPolarResolutionM};
  for (std::size_t i = 0; i < std::size(meshLengthKeys); ++i) {
    if (!ParsePositiveFinite(settings[meshLengthKeys[i]].value,
                             meshLengthTargets[i])) {
      *error = At(settings[meshLengthKeys[i]],
          std::string(meshLengthKeys[i]) + " must be finite and positive");
      return false;
    }
  }
  if (candidate.meshFarFieldCellSizeM <
      candidate.meshNearSurfaceCellSizeM) {
    *error = At(settings["mesh_far_field_cell_size_m"],
        "mesh_far_field_cell_size_m cannot be smaller than "
        "mesh_near_surface_cell_size_m");
    return false;
  }
  if (candidate.meshFarFieldExtentM <=
      candidate.meshNearSurfaceExtentM) {
    *error = At(settings["mesh_far_field_extent_m"],
        "mesh_far_field_extent_m must be greater than "
        "mesh_near_surface_extent_m");
    return false;
  }
  if (candidate.meshOuterCellSizeM <
      candidate.meshFarFieldCellSizeM) {
    *error = At(settings["mesh_outer_cell_size_m"],
        "mesh_outer_cell_size_m cannot be smaller than "
        "mesh_far_field_cell_size_m");
    return false;
  }
  if (!ParseBoolean(settings["mesh_output_tecplot"].value,
                    &candidate.meshOutputTecplot)) {
    *error = At(settings["mesh_output_tecplot"],
        "mesh_output_tecplot must be true or false");
    return false;
  }
  if (candidate.meshTerminatorHalfWidthDeg >= 90.0) {
    *error = At(settings["mesh_terminator_half_width_deg"],
        "mesh_terminator_half_width_deg must be less than 90 degrees");
    return false;
  }
  if (candidate.meshPolarCapStartLatitudeDeg >= 90.0) {
    *error = At(settings["mesh_polar_cap_start_latitude_deg"],
        "mesh_polar_cap_start_latitude_deg must be less than 90 degrees");
    return false;
  }
  if (!ParseAngularResolutionTable(
          settings["mesh_angular_resolution_m_by_sza"].value,
          &candidate.meshAngularSzaDeg,
          &candidate.meshAngularResolutionM)) {
    *error = At(settings["mesh_angular_resolution_m_by_sza"],
        "mesh_angular_resolution_m_by_sza must be an increasing "
        "SZA_DEG:RESOLUTION_M table beginning at 0, containing 90, and "
        "ending at 180 degrees");
    return false;
  }

  // Parse the complete M1 sampling receipt before paths are touched.  Every
  // dimensional input names its SI unit; there are no implicit kilometres or
  // radii.  The explicit checks below prevent a plausible-looking but
  // discontinuous grid from reaching an expensive linked campaign.
  const char* positiveDoubleKeys[] = {
      "profile_altitude_near_maximum_m", "profile_altitude_near_step_m",
      "profile_altitude_middle_maximum_m", "profile_altitude_middle_step_m",
      "profile_altitude_far_maximum_m",
      "limb_tangent_maximum_altitude_m", "tail_map_half_width_deg",
      "tail_beam_diameter_deg", "tail_exposure_duration_s"};
  double* positiveDoubleTargets[] = {
      &candidate.profileAltitudeNearMaximumM,
      &candidate.profileAltitudeNearStepM,
      &candidate.profileAltitudeMiddleMaximumM,
      &candidate.profileAltitudeMiddleStepM,
      &candidate.profileAltitudeFarMaximumM,
      &candidate.limbTangentMaximumAltitudeM,
      &candidate.tailMapHalfWidthDeg,
      &candidate.tailBeamDiameterDeg,
      &candidate.tailExposureDurationS};
  for (std::size_t i = 0; i < std::size(positiveDoubleKeys); ++i) {
    if (!ParsePositiveFinite(settings[positiveDoubleKeys[i]].value,
                             positiveDoubleTargets[i])) {
      *error = At(settings[positiveDoubleKeys[i]],
          std::string(positiveDoubleKeys[i]) +
          " must be finite and positive");
      return false;
    }
  }
  if (!ParsePositiveInteger(
          settings["profile_altitude_far_points_per_decade"].value,
          &candidate.profileAltitudeFarPointsPerDecade)) {
    *error = At(settings["profile_altitude_far_points_per_decade"],
        "profile_altitude_far_points_per_decade must be a positive integer");
    return false;
  }
  if (!ParseSolarZenithAngles(
          settings["profile_solar_zenith_angles_deg"].value,
          &candidate.profileSolarZenithAnglesDeg)) {
    *error = At(settings["profile_solar_zenith_angles_deg"],
        "profile_solar_zenith_angles_deg must be a nonempty, strictly "
        "increasing comma-separated list in [0,180]");
    return false;
  }
  if (!ParseBoolean(settings["profile_include_body_fixed_poles"].value,
                    &candidate.profileIncludeBodyFixedPoles)) {
    *error = At(settings["profile_include_body_fixed_poles"],
        "profile_include_body_fixed_poles must be true or false");
    return false;
  }
  if (!ParseBoolean(settings["profile_include_shell_average"].value,
                    &candidate.profileIncludeShellAverage)) {
    *error = At(settings["profile_include_shell_average"],
        "profile_include_shell_average must be true or false");
    return false;
  }

  candidate.simulationStartUtc =
      Unquote(Trim(settings["simulation_start_utc"].value));
  if (candidate.simulationStartUtc.empty()) {
    *error = At(settings["simulation_start_utc"],
                "simulation_start_utc cannot be empty");
    return false;
  }
  candidate.simulationEnvironmentMode=Lower(Unquote(Trim(
      settings["simulation_environment_mode"].value)));
  if (candidate.simulationEnvironmentMode!="frozen" &&
      candidate.simulationEnvironmentMode!="evolving") {
    *error=At(settings["simulation_environment_mode"],
        "simulation_environment_mode must be frozen or evolving");
    return false;
  }
  if (!ParsePositiveInteger(settings["random_seed"].value,
                            &candidate.randomSeed)) {
    *error=At(settings["random_seed"],
        "random_seed must be a positive base-10 integer");
    return false;
  }

  candidate.limbObserver =
      Lower(Unquote(Trim(settings["limb_observer"].value)));
  candidate.limbPlane = Lower(Unquote(Trim(settings["limb_plane"].value)));
  candidate.limbSides = Lower(Unquote(Trim(settings["limb_sides"].value)));
  candidate.limbBranchReference =
      Lower(Unquote(Trim(settings["limb_branch_reference"].value)));
  candidate.profileShellAngularWeighting = Lower(Unquote(
      Trim(settings["profile_shell_angular_weighting"].value)));
  candidate.tailCenterDirection =
      Lower(Unquote(Trim(settings["tail_center_direction"].value)));
  candidate.tailMapCoordinateSystem = Lower(Unquote(
      Trim(settings["tail_map_coordinate_system"].value)));
  candidate.tailMapDisplayOrientation = Lower(Unquote(
      Trim(settings["tail_map_display_orientation"].value)));
  candidate.tailBeamResponse =
      Lower(Unquote(Trim(settings["tail_beam_response"].value)));
  candidate.tailTemporalOperator =
      Lower(Unquote(Trim(settings["tail_temporal_operator"].value)));
  if (candidate.profileShellAngularWeighting != "solid-angle") {
    *error = At(settings["profile_shell_angular_weighting"],
        "profile_shell_angular_weighting currently supports only "
        "solid-angle; AMR-cell or unweighted-facet averages are biased");
    return false;
  }
  if (candidate.limbObserver != "earth") {
    *error = At(settings["limb_observer"],
        "limb_observer currently supports only earth");
    return false;
  }
  if (candidate.limbPlane != "sun-moon-earth") {
    *error = At(settings["limb_plane"],
        "limb_plane currently supports only sun-moon-earth");
    return false;
  }
  if (candidate.limbSides != "both" &&
      candidate.limbSides != "toward-projected-sun" &&
      candidate.limbSides != "away-from-projected-sun") {
    *error = At(settings["limb_sides"],
        "limb_sides must be both, toward-projected-sun, or "
        "away-from-projected-sun");
    return false;
  }
  if (candidate.limbBranchReference != "projected-sun") {
    *error = At(settings["limb_branch_reference"],
        "limb_branch_reference currently supports only projected-sun");
    return false;
  }
  if (candidate.tailCenterDirection != "anti-lunar") {
    *error = At(settings["tail_center_direction"],
        "tail_center_direction currently supports only anti-lunar");
    return false;
  }
  if (candidate.tailMapCoordinateSystem != "j2000-ra-dec" &&
      candidate.tailMapCoordinateSystem != "eclipj2000-lon-lat") {
    *error = At(settings["tail_map_coordinate_system"],
        "tail_map_coordinate_system must be j2000-ra-dec or "
        "eclipj2000-lon-lat");
    return false;
  }
  if (candidate.tailMapDisplayOrientation != "west-right" &&
      candidate.tailMapDisplayOrientation != "east-right") {
    *error = At(settings["tail_map_display_orientation"],
        "tail_map_display_orientation must be west-right or east-right");
    return false;
  }
  // D14 identifies a circular field of view but contains no authoritative
  // radial response function.  Fail closed on any assumed response law.  A
  // future tabulated response can be added only with a defined file schema,
  // provenance, normalization, and production convolution implementation.
  if (candidate.tailBeamResponse != "unresolved") {
    *error = At(settings["tail_beam_response"],
        "tail_beam_response must remain unresolved until an authoritative "
        "radial response and a production convolution are available");
    return false;
  }
  if (candidate.tailTemporalOperator != "exposure-average") {
    *error = At(settings["tail_temporal_operator"],
        "tail_temporal_operator currently supports only exposure-average");
    return false;
  }
  if (candidate.tailMapHalfWidthDeg >= 180.0) {
    *error = At(settings["tail_map_half_width_deg"],
        "tail_map_half_width_deg must be less than 180 degrees");
    return false;
  }
  if (candidate.tailBeamDiameterDeg >
      2.0 * candidate.tailMapHalfWidthDeg) {
    *error = At(settings["tail_beam_diameter_deg"],
        "tail_beam_diameter_deg cannot exceed the full map width");
    return false;
  }

  std::string gridReason;
  if (!BuildProfileAltitudeGrid(&candidate, &gridReason)) {
    *error = At(inputFile, 0, gridReason);
    return false;
  }
  if (candidate.limbTangentMaximumAltitudeM >
      candidate.profileAltitudeFarMaximumM) {
    *error = At(settings["limb_tangent_maximum_altitude_m"],
        "limb tangent maximum cannot exceed the profile far maximum");
    return false;
  }
  const auto tangentNode = std::find_if(
      candidate.profileAltitudeGridM.begin(),
      candidate.profileAltitudeGridM.end(),
      [&candidate](double altitude) {
        return std::abs(altitude - candidate.limbTangentMaximumAltitudeM) <=
            32.0 * std::numeric_limits<double>::epsilon() *
            std::max(1.0, candidate.limbTangentMaximumAltitudeM);
      });
  if (tangentNode == candidate.profileAltitudeGridM.end()) {
    *error = At(settings["limb_tangent_maximum_altitude_m"],
        "limb_tangent_maximum_altitude_m must coincide with a profile "
        "altitude grid point");
    return false;
  }

  candidate.lolaProductId = Unquote(settings["lola_product_id"].value);
  if (candidate.lolaProductId.empty()) {
    *error = At(settings["lola_product_id"],
                "lola_product_id cannot be empty");
    return false;
  }

  // Parse the M2 group only after all M1 scalar contracts have succeeded.  The
  // group is nevertheless committed transactionally with the whole receipt.
  fs::path divinerTable,divinerEclipseHistory,divinerEclipseContract,
      inventoryCheckpoint;
  if (candidate.m2EnvironmentEnabled) {
    candidate.surfaceTemperatureMode=Lower(Unquote(Trim(
        settings["surface_temperature_mode"].value)));
    if (candidate.surfaceTemperatureMode!="simple-cosine" &&
        candidate.surfaceTemperatureMode!="diviner-map" &&
        candidate.surfaceTemperatureMode!="thermal-model") {
      *error=At(settings["surface_temperature_mode"],
          "surface_temperature_mode must be simple-cosine, diviner-map, or "
          "thermal-model");
      return false;
    }
    const char* positiveM2Keys[]={
        "illumination_ray_origin_offset_m",
        "surface_collision_maximum_segment_length_m",
        "surface_collision_time_tolerance_s",
        "surface_collision_distance_tolerance_m",
        "thermal_solar_constant_w_m2",
        "thermal_infrared_emissivity", "thermal_surface_density_kg_m3",
        "thermal_deep_density_kg_m3", "thermal_density_scale_height_m",
        "thermal_surface_conductivity_w_m_k",
        "thermal_deep_conductivity_w_m_k",
        "thermal_radiative_conductivity_parameter",
        "thermal_geothermal_flux_w_m2", "thermal_bottom_depth_m",
        "thermal_time_step_s", "thermal_initial_temperature_k",
        "thermal_nonlinear_tolerance_k",
        "thermal_solar_disc_radius_m", "thermal_opaque_earth_radius_m",
        "thermal_spinup_temperature_tolerance_k",
        "thermal_spinup_energy_tolerance_j_m2",
        "sodium_activation_energy_ev",
        "sodium_attempt_frequency_hz", "argon_desorption_energy_j_mol",
        "cold_trap_retention_timescale_s",
        "surface_inventory_checkpoint_interval_s"};
    double* positiveM2Targets[]={
        &candidate.illuminationRayOriginOffsetM,
        &candidate.surfaceCollisionMaximumSegmentLengthM,
        &candidate.surfaceCollisionTimeToleranceS,
        &candidate.surfaceCollisionDistanceToleranceM,
        &candidate.thermal.solarConstantWm2,
        &candidate.thermal.infraredEmissivity,
        &candidate.thermal.surfaceDensityKgM3,
        &candidate.thermal.deepDensityKgM3,
        &candidate.thermal.densityScaleHeightM,
        &candidate.thermal.surfaceConductivityWmK,
        &candidate.thermal.deepConductivityWmK,
        &candidate.thermal.radiativeConductivityParameter,
        &candidate.thermal.geothermalFluxWm2,
        &candidate.thermal.bottomDepthM,&candidate.thermalTimeStepS,
        &candidate.thermalInitialTemperatureK,
        &candidate.thermal.nonlinearToleranceK,
        &candidate.thermalSolarDiscRadiusM,
        &candidate.thermalOpaqueEarthRadiusM,
        &candidate.thermalSpinupTemperatureToleranceK,
        &candidate.thermalSpinupEnergyToleranceJm2,
        &candidate.sodiumActivationEnergyEv,
        &candidate.sodiumAttemptFrequencyHz,
        &candidate.argonDesorptionEnergyJmol,
        &candidate.coldTrapRetentionTimescaleS,
        &candidate.surfaceInventoryCheckpointIntervalS};
    for (std::size_t index=0;index<std::size(positiveM2Keys);++index)
      if (!ParsePositiveFinite(settings[positiveM2Keys[index]].value,
                               positiveM2Targets[index])) {
        *error=At(settings[positiveM2Keys[index]],
            std::string(positiveM2Keys[index])+" must be finite and positive");
        return false;
      }
    if (!ParsePositiveInteger(settings["thermal_layer_count"].value,
                              &candidate.thermal.layerCount) ||
        candidate.thermal.layerCount<2 ||
        !ParsePositiveInteger(
            settings["thermal_nonlinear_maximum_iterations"].value,
            &candidate.thermal.nonlinearMaximumIterations)) {
      *error=At(inputFile,0,
          "thermal_layer_count must be at least 2 and the nonlinear "
          "iteration limit must be a positive integer");
      return false;
    }
    if (!ParsePositiveInteger(
            settings["surface_collision_maximum_iterations"].value,
            &candidate.surfaceCollisionMaximumIterations)) {
      *error=At(settings["surface_collision_maximum_iterations"],
          "surface_collision_maximum_iterations must be a positive integer");
      return false;
    }
    const char* boundedKeys[]={"thermal_normal_bond_albedo",
        "sodium_accommodation_coefficient",
        "argon_accommodation_coefficient"};
    double* boundedTargets[]={&candidate.thermal.normalBondAlbedo,
        &candidate.sodiumAccommodationCoefficient,
        &candidate.argonAccommodationCoefficient};
    for (std::size_t index=0;index<std::size(boundedKeys);++index) {
      const std::string scalar=Unquote(Trim(settings[boundedKeys[index]].value));
      errno=0;
      char* end=nullptr;
      const double value=std::strtod(scalar.c_str(),&end);
      if (errno==ERANGE || end==scalar.c_str() || *end!='\0' ||
          !std::isfinite(value) || value<0.0 || value>1.0) {
        *error=At(settings[boundedKeys[index]],
            std::string(boundedKeys[index])+" must be in [0,1]");
        return false;
      }
      *boundedTargets[index]=value;
    }
    const char* finiteKeys[]={"thermal_albedo_coefficient_a",
        "thermal_albedo_coefficient_b"};
    double* finiteTargets[]={&candidate.thermal.albedoCoefficientA,
        &candidate.thermal.albedoCoefficientB};
    for (std::size_t index=0;index<std::size(finiteKeys);++index) {
      std::vector<double> one;
      if (!ParseFiniteList(settings[finiteKeys[index]].value,&one) ||
          one.size()!=1) {
        *error=At(settings[finiteKeys[index]],
            std::string(finiteKeys[index])+" must be one finite scalar");
        return false;
      }
      *finiteTargets[index]=one.front();
    }
    if (!ParseFiniteList(
            settings["thermal_heat_capacity_coefficients_j_kg_k"].value,
            &candidate.thermal.heatCapacityPolynomial) ||
        candidate.thermal.heatCapacityPolynomial.size()!=5) {
      *error=At(settings["thermal_heat_capacity_coefficients_j_kg_k"],
          "thermal heat capacity requires exactly c0,c1,c2,c3,c4");
      return false;
    }
    if (!ParseProbabilityTable(
            settings["sodium_sticking_probability_by_temperature_k"].value,
            &candidate.sodiumStickingTemperatureK,
            &candidate.sodiumStickingProbability)) {
      *error=At(settings["sodium_sticking_probability_by_temperature_k"],
          "sodium sticking table must be increasing TEMPERATURE_K:PROBABILITY "
          "pairs with probabilities in [0,1]");
      return false;
    }
    if (!ParseProbabilityTable(
            settings["argon_sticking_probability_by_temperature_k"].value,
            &candidate.argonStickingTemperatureK,
            &candidate.argonStickingProbability)) {
      *error=At(settings["argon_sticking_probability_by_temperature_k"],
          "argon sticking table must be increasing TEMPERATURE_K:PROBABILITY "
          "pairs with probabilities in [0,1]");
      return false;
    }

    // Model names are required even though the numerical table/factor is also
    // explicit.  This makes the mathematical interpretation reviewable in a
    // run manifest: equal table values do not silently imply a constant law,
    // and an accommodation factor does not silently choose speed versus
    // kinetic-energy interpolation.
    candidate.sodiumStickingModel=Lower(Unquote(Trim(
        settings["sodium_sticking_model"].value)));
    candidate.argonStickingModel=Lower(Unquote(Trim(
        settings["argon_sticking_model"].value)));
    for (const auto& selection : std::array<std::pair<const char*,
             const std::string*>,2>{{
          {"sodium_sticking_model",&candidate.sodiumStickingModel},
          {"argon_sticking_model",&candidate.argonStickingModel}}}) {
      if (*selection.second!="temperature-table" &&
          *selection.second!="constant") {
        *error=At(settings[selection.first],std::string(selection.first)+
            " must be temperature-table or constant");
        return false;
      }
    }
    const auto constantTable=[](const std::vector<double>& probability) {
      return std::all_of(probability.begin(),probability.end(),
          [&probability](double value) { return value==probability.front(); });
    };
    if (candidate.sodiumStickingModel=="constant" &&
        !constantTable(candidate.sodiumStickingProbability)) {
      *error=At(settings["sodium_sticking_probability_by_temperature_k"],
          "constant sodium sticking model requires identical probabilities "
          "at every table node");
      return false;
    }
    if (candidate.argonStickingModel=="constant" &&
        !constantTable(candidate.argonStickingProbability)) {
      *error=At(settings["argon_sticking_probability_by_temperature_k"],
          "constant argon sticking model requires identical probabilities "
          "at every table node");
      return false;
    }
    candidate.sodiumAccommodationSpeedModel=Lower(Unquote(Trim(
        settings["sodium_accommodation_speed_model"].value)));
    candidate.argonAccommodationSpeedModel=Lower(Unquote(Trim(
        settings["argon_accommodation_speed_model"].value)));
    for (const auto& selection : std::array<std::pair<const char*,
             const std::string*>,2>{{
          {"sodium_accommodation_speed_model",
              &candidate.sodiumAccommodationSpeedModel},
          {"argon_accommodation_speed_model",
              &candidate.argonAccommodationSpeedModel}}}) {
      if (*selection.second!="linear-speed" &&
          *selection.second!="kinetic-energy") {
        *error=At(settings[selection.first],std::string(selection.first)+
            " must be linear-speed or kinetic-energy");
        return false;
      }
    }
    candidate.sodiumSurfaceInteractionReference=Unquote(Trim(
        settings["sodium_surface_interaction_reference"].value));
    candidate.argonSurfaceInteractionReference=Unquote(Trim(
        settings["argon_surface_interaction_reference"].value));
    if (candidate.sodiumSurfaceInteractionReference.empty() ||
        candidate.argonSurfaceInteractionReference.empty()) {
      *error=At(inputFile,0,
          "each configured neutral surface law requires a non-empty "
          "surface_interaction_reference receipt");
      return false;
    }

    // A source receipt is separate from the adsorption/residence receipt.
    // Current literature supports a continuous global baseline rate, while
    // the available D10 product is regional and its QA status is incomplete.
    // Reject potassium-map selections rather than filling unmeasured regions
    // or treating a regional normalized proxy as a global source.
    candidate.argonSourceSpatialModel=Lower(Unquote(Trim(
        settings["argon_source_spatial_model"].value)));
    if (!ParseBoolean(settings["argon_source_enabled"].value,
                      &candidate.argonSourceEnabled)) {
      *error=At(settings["argon_source_enabled"],
          "argon_source_enabled must be true or false");
      return false;
    }
    candidate.argonSourceTimeModel=Lower(Unquote(Trim(
        settings["argon_source_time_model"].value)));
    candidate.argonSourceMapFile=Unquote(Trim(
        settings["argon_source_map_file"].value));
    candidate.argonSourceReference=Unquote(Trim(
        settings["argon_source_reference"].value));
    if (!ParseNonnegativeFinite(settings["argon_source_rate_atoms_s"].value,
                                &candidate.argonSourceRateAtomsS)) {
      *error=At(settings["argon_source_rate_atoms_s"],
          "argon_source_rate_atoms_s must be finite and non-negative");
      return false;
    }
    if (candidate.argonSourceSpatialModel!="uniform-global") {
      *error=At(settings["argon_source_spatial_model"],
          "argon_source_spatial_model currently supports only uniform-global; "
          "potassium-weighted production is unavailable while D10 remains "
          "regional and QA-incomplete");
      return false;
    }
    if (candidate.argonSourceTimeModel!="continuous") {
      *error=At(settings["argon_source_time_model"],
          "argon_source_time_model currently supports only continuous");
      return false;
    }
    if (Lower(candidate.argonSourceMapFile)!="none") {
      *error=At(settings["argon_source_map_file"],
          "argon_source_map_file must be none until a global QA-clean D10 "
          "source map is available");
      return false;
    }
    if (!(candidate.argonSourceRateAtomsS>0.0) ||
        candidate.argonSourceReference.empty() ||
        Lower(candidate.argonSourceReference)=="none") {
      *error=At(inputFile,0,
          "uniform-global argon source parameters require a positive rate "
          "and a non-empty literature reference, even when injection is "
          "disabled; this preserves an auditable candidate receipt");
      return false;
    }
    if (candidate.thermal.infraredEmissivity>1.0 ||
        candidate.thermal.normalBondAlbedo>=1.0 ||
        candidate.thermal.deepDensityKgM3<=
            candidate.thermal.surfaceDensityKgM3 ||
        candidate.thermal.deepConductivityWmK<
            candidate.thermal.surfaceConductivityWmK) {
      *error=At(inputFile,0,
          "M2 thermal emissivity/albedo/density/conductivity ordering is "
          "not physical");
      return false;
    }

    // The first I17 implementation has one deliberately narrow physical
    // contract.  Reject aliases and future-looking names so a receipt can
    // never imply atmospheric transmission or limb darkening that the
    // production forcing does not actually calculate.
    candidate.thermalEclipseDirectIrradianceModel=Lower(Unquote(Trim(
        settings["thermal_eclipse_direct_irradiance_model"].value)));
    candidate.thermalSolarDiscRadianceModel=Lower(Unquote(Trim(
        settings["thermal_solar_disc_radiance_model"].value)));
    candidate.thermalEarthAtmosphereTransmissionModel=Lower(Unquote(Trim(
        settings["thermal_earth_atmosphere_transmission_model"].value)));
    if (candidate.thermalEclipseDirectIrradianceModel!=
            "apparent-disc-overlap" ||
        candidate.thermalSolarDiscRadianceModel!="uniform" ||
        candidate.thermalEarthAtmosphereTransmissionModel!="disabled") {
      *error=At(inputFile,0,
          "the approved experimental I17 forcing requires "
          "thermal_eclipse_direct_irradiance_model=apparent-disc-overlap, "
          "thermal_solar_disc_radiance_model=uniform, and "
          "thermal_earth_atmosphere_transmission_model=disabled; no other "
          "eclipse radiance model is implemented");
      return false;
    }
    candidate.thermalSpinupCycleStartUtc=Unquote(Trim(
        settings["thermal_spinup_cycle_start_utc"].value));
    candidate.thermalSpinupCycleEndUtc=Unquote(Trim(
        settings["thermal_spinup_cycle_end_utc"].value));
    if (candidate.thermalSpinupCycleStartUtc.empty() ||
        candidate.thermalSpinupCycleEndUtc.empty() ||
        Lower(candidate.thermalSpinupCycleStartUtc)=="none" ||
        Lower(candidate.thermalSpinupCycleEndUtc)=="none" ||
        candidate.thermalSpinupCycleStartUtc==
            candidate.thermalSpinupCycleEndUtc) {
      *error=At(inputFile,0,
          "thermal spin-up requires distinct, non-empty start/end UTC "
          "values; CSPICE validates their time syntax and ordering before "
          "the first spin-up cycle");
      return false;
    }
    candidate.ionSurfaceInteractionMode=Lower(Unquote(Trim(
        settings["ion_surface_interaction_mode"].value)));
    if (candidate.ionSurfaceInteractionMode!="unconfigured") {
      *error=At(settings["ion_surface_interaction_mode"],
          "ion_surface_interaction_mode must remain unconfigured until "
          "authoritative return/neutralization branching data are supplied");
      return false;
    }
    candidate.surfaceInventoryInitialState=Lower(Unquote(Trim(
        settings["surface_inventory_initial_state"].value)));
    if (candidate.surfaceInventoryInitialState!="zero" &&
        candidate.surfaceInventoryInitialState!="checkpoint") {
      *error=At(settings["surface_inventory_initial_state"],
          "surface_inventory_initial_state must be zero or checkpoint");
      return false;
    }
    const std::string divinerValue=Lower(Unquote(Trim(
        settings["diviner_temperature_table_file"].value)));
    if (divinerValue!="none" &&
        !AbsolutePath(settings["diviner_temperature_table_file"],
                      "diviner_temperature_table_file",&divinerTable,error))
      return false;
    if (candidate.surfaceTemperatureMode=="diviner-map" &&
        (divinerValue=="none" || !fs::is_regular_file(divinerTable))) {
      *error=At(settings["diviner_temperature_table_file"],
          "diviner-map mode requires an existing processed D03 table");
      return false;
    }
    candidate.divinerTemperatureTableFile=
        divinerValue=="none" ? std::string() : divinerTable.string();

    // Observation files are optional for ordinary M2 transport runs but must
    // be explicit.  I17 derives a run receipt that replaces both ``none``
    // values; requiring an existing pair here prevents half-configured
    // evidence from reaching the linked executable.
    const std::string eclipseHistoryValue=Lower(Unquote(Trim(
        settings["diviner_eclipse_history_file"].value)));
    const std::string eclipseContractValue=Lower(Unquote(Trim(
        settings["diviner_eclipse_selection_contract_file"].value)));
    if ((eclipseHistoryValue=="none")!=(eclipseContractValue=="none")) {
      *error=At(inputFile,0,
          "Diviner eclipse history and selection contract must both be none "
          "or both name existing files");
      return false;
    }
    if (eclipseHistoryValue!="none") {
      if (!AbsolutePath(settings["diviner_eclipse_history_file"],
              "diviner_eclipse_history_file",&divinerEclipseHistory,error) ||
          !AbsolutePath(settings["diviner_eclipse_selection_contract_file"],
              "diviner_eclipse_selection_contract_file",
              &divinerEclipseContract,error)) return false;
      if (!fs::is_regular_file(divinerEclipseHistory) ||
          !fs::is_regular_file(divinerEclipseContract)) {
        *error=At(inputFile,0,
            "configured Diviner eclipse history/contract files must exist");
        return false;
      }
      candidate.divinerEclipseHistoryFile=divinerEclipseHistory.string();
      candidate.divinerEclipseSelectionContractFile=
          divinerEclipseContract.string();
    }
    const std::string checkpointValue=Lower(Unquote(Trim(
        settings["surface_inventory_checkpoint_file"].value)));
    if (checkpointValue!="none" &&
        !AbsolutePath(settings["surface_inventory_checkpoint_file"],
                      "surface_inventory_checkpoint_file",
                      &inventoryCheckpoint,error)) return false;
    if (candidate.surfaceInventoryInitialState=="checkpoint") {
      // The setting names a checkpoint *base*.  M2 writes one independently
      // versioned file per modeled neutral species so a partial restart can
      // never be mistaken for a complete Na+Ar surface state.  Validate the
      // concrete files rather than the base path, which is not itself a file.
      if (checkpointValue=="none" ||
          !fs::is_regular_file(inventoryCheckpoint.string()+".Na") ||
          !fs::is_regular_file(inventoryCheckpoint.string()+".Ar")) {
        *error=At(settings["surface_inventory_checkpoint_file"],
            "checkpoint initial state requires existing <base>.Na and "
            "<base>.Ar files");
        return false;
      }
    }
    candidate.surfaceInventoryCheckpointFile=checkpointValue=="none" ?
        std::string() : inventoryCheckpoint.string();
  }

  fs::path spiceRoot, image, label, cea, tecplot;
  if (!AbsolutePath(settings["spice_path"], "spice_path",
                    &spiceRoot, error) ||
      !AbsolutePath(settings["lola_image_file"],
                    "lola_image_file", &image, error) ||
      !AbsolutePath(settings["lola_label_file"],
                    "lola_label_file", &label, error) ||
      !AbsolutePath(settings["surface_cea_file"],
                    "surface_cea_file", &cea, error) ||
      !AbsolutePath(settings["surface_tecplot_file"],
                    "surface_tecplot_file", &tecplot, error)) {
    return false;
  }

  // The requested root follows /home/vtenishe/SPICE: toolkit files live in
  // cspice/, while mission and generic kernels live in Kernels/.
  const fs::path toolkit = spiceRoot / "cspice";
  const fs::path kernels = spiceRoot / "Kernels";
  if (!fs::is_regular_file(toolkit / "include" / "SpiceUsr.h") ||
      !fs::is_regular_file(toolkit / "lib" / "cspice.a") ||
      !fs::is_directory(kernels)) {
    *error = At(settings["spice_path"],
        "spice_path must contain cspice/include/SpiceUsr.h, "
        "cspice/lib/cspice.a, and Kernels/");
    return false;
  }
  // Validate inputs now; output parents are created later by the atomic mesh
  // writers.  std::ifstream still provides the final readable-file check when
  // the DEM is loaded, including permission and concurrent-removal failures.
  if (!fs::is_regular_file(image) || !fs::is_regular_file(label)) {
    *error = At(inputFile, 0,
        "lola_image_file and lola_label_file must name readable files");
    return false;
  }
  if (cea == tecplot) {
    *error = At(inputFile, 0,
                "surface_cea_file and surface_tecplot_file must differ");
    return false;
  }

  candidate.spiceRoot = spiceRoot.string();
  candidate.spiceToolkitDirectory = toolkit.string();
  candidate.spiceKernelDirectory = kernels.string();
  candidate.lolaImageFile = image.string();
  candidate.lolaLabelFile = label.string();
  candidate.surfaceCeaFile = cea.string();
  candidate.surfaceTecplotFile = tecplot.string();
  *result = std::move(candidate);
  return true;
}

double EvaluateAngularMeshResolutionM(const Configuration& configuration,
                                      double solarZenithAngleDeg,
                                      double absoluteBodyLatitudeDeg) {
  // ParseApplicationInput guarantees this invariant for installed production
  // configurations.  Retain a fatal guard here because the function is also a
  // public production kernel used by thin verification adapters.
  if (configuration.meshAngularSzaDeg.size()<3 ||
      configuration.meshAngularSzaDeg.size()!=
          configuration.meshAngularResolutionM.size() ||
      !std::isfinite(solarZenithAngleDeg) || solarZenithAngleDeg<0.0 ||
      solarZenithAngleDeg>180.0 || !std::isfinite(absoluteBodyLatitudeDeg) ||
      absoluteBodyLatitudeDeg<0.0 || absoluteBodyLatitudeDeg>90.0) {
    return std::numeric_limits<double>::quiet_NaN();
  }

  const std::vector<double>& angles=configuration.meshAngularSzaDeg;
  const std::vector<double>& resolutions=
      configuration.meshAngularResolutionM;
  const auto upper=std::upper_bound(angles.begin(),angles.end(),
                                    solarZenithAngleDeg);
  double requestedResolutionM;
  if (upper==angles.begin()) requestedResolutionM=resolutions.front();
  else if (upper==angles.end()) requestedResolutionM=resolutions.back();
  else {
    const std::size_t right=static_cast<std::size_t>(upper-angles.begin());
    const std::size_t left=right-1;
    const double fraction=(solarZenithAngleDeg-angles[left])/
        (angles[right]-angles[left]);
    requestedResolutionM=resolutions[left]+fraction*
        (resolutions[right]-resolutions[left]);
  }

  if (std::abs(solarZenithAngleDeg-90.0)<=
      configuration.meshTerminatorHalfWidthDeg) {
    const auto terminator=std::lower_bound(angles.begin(),angles.end(),90.0);
    if (terminator==angles.end() || *terminator!=90.0)
      return std::numeric_limits<double>::quiet_NaN();
    requestedResolutionM=std::min(requestedResolutionM,
        resolutions[static_cast<std::size_t>(terminator-angles.begin())]);
  }
  if (absoluteBodyLatitudeDeg>=
      configuration.meshPolarCapStartLatitudeDeg) {
    requestedResolutionM=std::min(requestedResolutionM,
                                   configuration.meshPolarResolutionM);
  }
  return requestedResolutionM;
}

double EvaluateVolumeMeshResolutionM(const Configuration& configuration,
                                     double radiusM,
                                     double referenceRadiusM,
                                     double outerRadiusM) {
  if (!std::isfinite(radiusM) || !std::isfinite(referenceRadiusM) ||
      !std::isfinite(outerRadiusM) || radiusM<0.0 ||
      !(referenceRadiusM>0.0) || !(outerRadiusM>referenceRadiusM) ||
      !(configuration.meshNearSurfaceCellSizeM>0.0) ||
      !(configuration.meshNearSurfaceExtentM>0.0) ||
      !(configuration.meshFarFieldExtentM>
          configuration.meshNearSurfaceExtentM) ||
      configuration.meshFarFieldCellSizeM<
          configuration.meshNearSurfaceCellSizeM ||
      configuration.meshOuterCellSizeM<
          configuration.meshFarFieldCellSizeM) {
    return std::numeric_limits<double>::quiet_NaN();
  }
  const double altitude=std::max(0.0,radiusM-referenceRadiusM);
  if (altitude<=configuration.meshNearSurfaceExtentM)
    return configuration.meshNearSurfaceCellSizeM;

  const double nearTransitionRadius=referenceRadiusM+
      configuration.meshNearSurfaceExtentM;
  const double farTransitionRadius=referenceRadiusM+
      configuration.meshFarFieldExtentM;
  if (!(farTransitionRadius>nearTransitionRadius) ||
      !(outerRadiusM>farTransitionRadius))
    return std::numeric_limits<double>::quiet_NaN();

  // Within the declared M1 science domain, transition from the near-surface
  // target to the 40 km far-field target in dimensionless log(radius).  The
  // endpoints are exact, so changing the outer computational box cannot alter
  // the qualified 0--10,000 km spatial contract.
  if (radiusM<=farTransitionRadius) {
    const double fraction=std::log(radiusM/nearTransitionRadius)/
        std::log(farTransitionRadius/nearTransitionRadius);
    return configuration.meshNearSurfaceCellSizeM+fraction*
        (configuration.meshFarFieldCellSizeM-
         configuration.meshNearSurfaceCellSizeM);
  }

  // The legacy Moon domain extends to 400 lunar radii along the tail.  A
  // second explicit profile coarsens outside the M1 science region, avoiding
  // an implicit hidden constant while retaining monotonic resolution.
  if (radiusM>=outerRadiusM) return configuration.meshOuterCellSizeM;
  const double outerFraction=std::log(radiusM/farTransitionRadius)/
      std::log(outerRadiusM/farTransitionRadius);
  return configuration.meshFarFieldCellSizeM+outerFraction*
      (configuration.meshOuterCellSizeM-
       configuration.meshFarFieldCellSizeM);
}

bool InstallConfiguration(const Configuration& configuration,
                          std::string* error) {
  if (error == nullptr) return false;
  if (gInstalled) {
    *error = "srcMoon runtime configuration was already installed";
    return false;
  }
  gConfiguration = configuration;
  gInstalled = true;
  error->clear();
  return true;
}

bool HasConfiguration() { return gInstalled; }

const Configuration& GetConfiguration() {
  // The no-argument regression predates application input.  Preserve it by
  // returning an immutable empty configuration, while production callers use
  // HasConfiguration() to decide whether any fields are meaningful.
  if (!gInstalled) {
    static const Configuration legacy;
    return legacy;
  }
  return gConfiguration;
}

}  // namespace Runtime
}  // namespace Moon
