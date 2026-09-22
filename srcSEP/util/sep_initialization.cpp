#include "sep_initialization.h"

#include <algorithm>
#include <cctype>
#include <cerrno>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <limits>
#include <map>
#include <memory>
#include <set>
#include <sstream>

namespace SEP {
namespace Initialization {
namespace {

std::unique_ptr<const Configuration> gActive;

Transport::Status Error(const std::string& message) {
  return Transport::Status::Error(Transport::StatusCode::InvalidArgument,
                                  message);
}

std::string Trim(const std::string& value) {
  const std::size_t first = value.find_first_not_of(" \t\r\n");
  if (first == std::string::npos) return std::string();
  const std::size_t last = value.find_last_not_of(" \t\r\n");
  return value.substr(first, last - first + 1);
}

std::string Lower(std::string value) {
  std::transform(value.begin(), value.end(), value.begin(),
                 [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
  return value;
}

bool ParseDouble(const std::string& text, double* value) {
  errno = 0;
  char* end = NULL;
  const double parsed = std::strtod(text.c_str(), &end);
  if (text.empty() || end == text.c_str() || *end != '\0' || errno == ERANGE ||
      !std::isfinite(parsed)) return false;
  *value = parsed;
  return true;
}

bool ParseUnsigned64(const std::string& text, std::uint64_t* value) {
  errno = 0;
  char* end = NULL;
  if (text.empty() || text[0] == '-') return false;
  const unsigned long long parsed = std::strtoull(text.c_str(), &end, 10);
  if (end == text.c_str() || *end != '\0' || errno == ERANGE) return false;
  *value = static_cast<std::uint64_t>(parsed);
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

double Dot(const Vec3& left, const Vec3& right) {
  return left.x * right.x + left.y * right.y + left.z * right.z;
}
Vec3 Add(const Vec3& left, const Vec3& right) {
  return {left.x + right.x, left.y + right.y, left.z + right.z};
}
Vec3 Subtract(const Vec3& left, const Vec3& right) {
  return {left.x - right.x, left.y - right.y, left.z - right.z};
}
Vec3 Scale(double value, const Vec3& vector) {
  return {value * vector.x, value * vector.y, value * vector.z};
}
Vec3 Cross(const Vec3& left, const Vec3& right) {
  return {left.y * right.z - left.z * right.y,
          left.z * right.x - left.x * right.z,
          left.x * right.y - left.y * right.x};
}
double Norm(const Vec3& vector) { return std::sqrt(Dot(vector, vector)); }
Vec3 Unit(const Vec3& vector) {
  const double length = Norm(vector);
  return length > 0.0 ? Scale(1.0 / length, vector) : Vec3();
}
bool Finite(const Vec3& vector) {
  return std::isfinite(vector.x) && std::isfinite(vector.y) &&
         std::isfinite(vector.z);
}

double Clamp(double value, double lower, double upper) {
  return std::max(lower, std::min(value, upper));
}

double Profile(double coordinate, RefinementProfile profile, double exponent) {
  const double x = Clamp(coordinate, 0.0, 1.0);
  if (profile == RefinementProfile::Linear) return x;
  if (profile == RefinementProfile::PowerLaw) return std::pow(x, exponent);
  const double smooth = x * x * (3.0 - 2.0 * x);
  return std::pow(smooth, exponent);
}

// Complete set of physically effective canonical SWCME1D inputs.  Data knots
// are conditional because an empty list is the only physically correct value
// for ballistic/DBM, whereas data-driven kinematics requires both lists.  The
// deprecated sheath.compression_floor is intentionally absent: canonical RH
// physics rejects it because it has no effect.
const std::set<std::string>& RequiredSwcmeKeys() {
  static const std::set<std::string> keys = {
      "preset",
      "ambient.wind_speed", "ambient.density_1au",
      "ambient.magnetic_field_1au", "ambient.proton_temperature",
      "ambient.adiabatic_index", "ambient.alpha_to_proton_ratio",
      "ambient.electron_temperature", "ambient.alpha_temperature",
      "ambient.thermodynamic_closure", "parker.radial_polarity",
      "parker.sin_theta", "parker.source_radius", "cme.kinematics",
      "cme.launch_radius", "cme.launch_speed", "cme.drag_coefficient",
      "cme.extrapolation", "shock.region_mode",
      "shock.acceleration_mode", "shock.relative_source_weight_per_area",
      "geometry.sheath_thickness_1au", "geometry.ejecta_thickness_1au",
      "smoothing.shock_width_1au", "smoothing.leading_edge_width_1au",
      "smoothing.trailing_edge_width_1au", "sheath.ramp_power",
      "sheath.leading_edge_speed_factor", "ejecta.density_factor",
      "ejecta.speed_factor", "event.launch_epoch", "event.valid_from",
      "event.valid_until", "source.particle_mass", "source.charge_number",
      "source.energy_min", "source.energy_max", "source.reference_energy",
      "source.injection_efficiency", "source.normalization",
      "source.reference_intensity_si"};
  return keys;
}

// Positive-polarity outward tangent. Magnetic polarity changes B, not this
// geometric curve, so it is deliberately absent from initialization input.
Vec3 ParkerTangent(const Vec3& relativePosition,
                   const Configuration& configuration) {
  const double radius = Norm(relativePosition);
  if (!(radius > 0.0)) return Vec3();
  const Vec3 radial = Scale(1.0 / radius, relativePosition);
  const Vec3 axis = {0.0, 0.0, 1.0};
  const double sourceRadius =
      Norm(Subtract(configuration.parkerInitialPointM,
                    configuration.parkerOriginM));
  const double winding = configuration.solarRotationRateRadPerS *
      std::max(0.0, radius - sourceRadius) /
      configuration.solarWindSpeedMPerS;
  return Unit(Subtract(radial, Scale(winding, Cross(axis, radial))));
}

Vec3 ParkerDirectionAtRadius(double radius,
                             const Configuration& configuration) {
  const Vec3 source = Unit(Subtract(configuration.parkerInitialPointM,
                                    configuration.parkerOriginM));
  const double sourceRadius = configuration.innerRadiusM;
  const double angle = -configuration.solarRotationRateRadPerS *
      std::max(0.0, radius - sourceRadius) /
      configuration.solarWindSpeedMPerS;
  const double cosine = std::cos(angle);
  const double sine = std::sin(angle);
  // Rodrigues rotation about +Z, written explicitly to keep this small layer
  // independent of a general matrix library.
  return {cosine * source.x - sine * source.y,
          sine * source.x + cosine * source.y, source.z};
}

ObserverConfiguration* ObserverForSection(const std::string& section,
                                          Configuration* configuration) {
  const std::string prefix = "observer.";
  if (section.rfind(prefix, 0) != 0 || section.size() == prefix.size())
    return NULL;
  const std::string id = section.substr(prefix.size());
  for (char value : id)
    if (!std::isalnum(static_cast<unsigned char>(value)) && value != '_' &&
        value != '-') return NULL;
  for (ObserverConfiguration& observer : configuration->observers)
    if (observer.id == id) return &observer;
  ObserverConfiguration observer;
  observer.id = id;
  configuration->observers.push_back(observer);
  return &configuration->observers.back();
}

Transport::Status Apply(const std::string& section, const std::string& key,
                        const std::string& value, Configuration* c) {
  const std::string field = section + "." + key;
  auto invalid = [&]() { return Error("invalid value for '" + field + "': " + value); };
  if (field == "run.schema_version") {
    std::uint64_t parsed = 0;
    if (!ParseUnsigned64(value, &parsed) ||
        (parsed != 1 && parsed != 2 && parsed != 3))
      return invalid();
    c->schemaVersion = static_cast<unsigned>(parsed);
  } else if (field == "run.time_step_s") {
    if (!ParseDouble(value, &c->timeStepS)) return invalid();
  } else if (field == "injection.macroparticles_per_step") {
    if (!ParseUnsigned64(value, &c->macroparticlesPerStep)) return invalid();
  } else if (field == "species.particle_weight") {
    if (!ParseDouble(value, &c->particleWeight)) return invalid();
  } else if (field == "observer.heliocentric_radius_m") {
    if (!ParseDouble(value, &c->observerHeliocentricRadiusM)) return invalid();
  } else if (field == "output.mesh_tecplot_file") {
    c->meshTecplotFile = value;
  } else if (field == "output.field_line_tecplot_file") {
    c->fieldLineTecplotFile = value;
  } else if (field == "output.data_tecplot_file") {
    c->dataTecplotFile = value;
  } else if (field == "parker_spiral.origin_x_m") {
    if (!ParseDouble(value, &c->parkerOriginM.x)) return invalid();
  } else if (field == "parker_spiral.origin_y_m") {
    if (!ParseDouble(value, &c->parkerOriginM.y)) return invalid();
  } else if (field == "parker_spiral.origin_z_m") {
    if (!ParseDouble(value, &c->parkerOriginM.z)) return invalid();
  } else if (field == "parker_spiral.initial_x_m") {
    if (!ParseDouble(value, &c->parkerInitialPointM.x)) return invalid();
  } else if (field == "parker_spiral.initial_y_m") {
    if (!ParseDouble(value, &c->parkerInitialPointM.y)) return invalid();
  } else if (field == "parker_spiral.initial_z_m") {
    if (!ParseDouble(value, &c->parkerInitialPointM.z)) return invalid();
  } else if (field == "parker_spiral.length_m") {
    if (!ParseDouble(value, &c->parkerLengthM)) return invalid();
  } else if (field == "parker_spiral.point_count") {
    if (!ParseUnsigned64(value, &c->parkerPointCount)) return invalid();
  } else if (field == "domain.inner_radius_m") {
    if (!ParseDouble(value, &c->innerRadiusM)) return invalid();
  } else if (field == "domain.outer_radius_m") {
    if (!ParseDouble(value, &c->outerRadiusM)) return invalid();
  } else if (field == "mesh.global_cell_size_m") {
    if (!ParseDouble(value, &c->globalCellSizeM)) return invalid();
  } else if (field == "mesh.minimum_cell_size_m") {
    if (!ParseDouble(value, &c->minimumCellSizeM)) return invalid();
  } else if (field == "mesh.maximum_level") {
    std::uint64_t parsed = 0;
    if (!ParseUnsigned64(value, &parsed) ||
        parsed > std::numeric_limits<unsigned>::max()) return invalid();
    c->maximumMeshLevel = static_cast<unsigned>(parsed);
  } else if (field == "mesh.solar.enabled") {
    if (!ParseBool(value, &c->solarRefinementEnabled)) return invalid();
  } else if (field == "mesh.solar.surface_cell_size_m") {
    if (!ParseDouble(value, &c->solarSurfaceCellSizeM)) return invalid();
  } else if (field == "mesh.solar.transition_outer_radius_m") {
    if (!ParseDouble(value, &c->solarTransitionOuterRadiusM)) return invalid();
  } else if (field == "mesh.solar.profile") {
    const std::string normalized = Lower(value);
    if (normalized == "linear") c->solarProfile = RefinementProfile::Linear;
    else if (normalized == "power-law") c->solarProfile = RefinementProfile::PowerLaw;
    else if (normalized == "smoothstep") c->solarProfile = RefinementProfile::Smoothstep;
    else return invalid();
  } else if (field == "mesh.solar.exponent") {
    if (!ParseDouble(value, &c->solarExponent)) return invalid();
  } else if (field == "mesh.tube.enabled") {
    if (!ParseBool(value, &c->tubeRefinementEnabled)) return invalid();
  } else if (field == "mesh.tube.reference_radius_m") {
    if (!ParseDouble(value, &c->tubeReferenceRadiusM)) return invalid();
  } else if (field == "mesh.tube.radius_at_reference_m") {
    if (!ParseDouble(value, &c->tubeRadiusAtReferenceM)) return invalid();
  } else if (field == "mesh.tube.radius_mode") {
    const std::string normalized = Lower(value);
    if (normalized == "physical-constant")
      c->tubeRadiusMode = TubeRadiusMode::PhysicalConstant;
    else if (normalized == "constant-angular-width")
      c->tubeRadiusMode = TubeRadiusMode::ConstantAngularWidth;
    else return invalid();
  } else if (field == "mesh.tube.center_cell_size_m") {
    if (!ParseDouble(value, &c->tubeCenterCellSizeM)) return invalid();
  } else if (field == "mesh.tube.transverse_profile") {
    const std::string normalized = Lower(value);
    if (normalized == "linear") c->tubeProfile = RefinementProfile::Linear;
    else if (normalized == "power-law") c->tubeProfile = RefinementProfile::PowerLaw;
    else if (normalized == "smoothstep") c->tubeProfile = RefinementProfile::Smoothstep;
    else return invalid();
  } else if (field == "mesh.tube.transverse_exponent") {
    if (!ParseDouble(value, &c->tubeExponent)) return invalid();
  } else if (field == "background.parker.solar_wind_speed_m_per_s") {
    if (!ParseDouble(value, &c->solarWindSpeedMPerS)) return invalid();
  } else if (field == "background.parker.solar_rotation_rate_rad_per_s") {
    if (!ParseDouble(value, &c->solarRotationRateRadPerS)) return invalid();
  } else {
    ObserverConfiguration* observer = ObserverForSection(section, c);
    if (observer == NULL)
      return Error("unknown configuration key '" + field + "'");
    if (key == "heliocentric_radius_m") {
      if (!ParseDouble(value, &observer->heliocentricRadiusM)) return invalid();
    } else if (key == "minimum_energy_j") {
      if (!ParseDouble(value, &observer->minimumEnergyJ)) return invalid();
    } else if (key == "maximum_energy_j") {
      if (!ParseDouble(value, &observer->maximumEnergyJ)) return invalid();
    } else if (key == "energy_channels") {
      std::uint64_t parsed = 0;
      if (!ParseUnsigned64(value, &parsed) ||
          parsed > std::numeric_limits<unsigned>::max()) return invalid();
      observer->energyChannels = static_cast<unsigned>(parsed);
    } else if (key == "energy_spacing") {
      const std::string normalized = Lower(value);
      if (normalized == "logarithmic")
        observer->energySpacing = EnergyChannelSpacing::Logarithmic;
      else if (normalized == "linear")
        observer->energySpacing = EnergyChannelSpacing::Linear;
      else return invalid();
    } else if (key == "pitch_angle_bins") {
      std::uint64_t parsed = 0;
      if (!ParseUnsigned64(value, &parsed) ||
          parsed > std::numeric_limits<unsigned>::max()) return invalid();
      observer->pitchAngleBins = static_cast<unsigned>(parsed);
    } else {
      return Error("unknown configuration key '" + field + "'");
    }
  }
  return Transport::Status::Ok();
}

}  // namespace

const char* Name(RefinementProfile profile) {
  if (profile == RefinementProfile::Linear) return "linear";
  if (profile == RefinementProfile::PowerLaw) return "power-law";
  return "smoothstep";
}
const char* Name(TubeRadiusMode mode) {
  return mode == TubeRadiusMode::PhysicalConstant
      ? "physical-constant" : "constant-angular-width";
}
const char* Name(EnergyChannelSpacing spacing) {
  return spacing == EnergyChannelSpacing::Logarithmic
      ? "logarithmic" : "linear";
}

Transport::Status Validate(const Configuration& c) {
  if (!Finite(c.parkerOriginM) || !Finite(c.parkerInitialPointM))
    return Error("Parker origin or initial point is not finite");
  const double values[] = {
      c.parkerOriginM.x, c.parkerOriginM.y, c.parkerOriginM.z,
      c.parkerInitialPointM.x, c.parkerInitialPointM.y,
      c.parkerInitialPointM.z, c.parkerLengthM, c.solarWindSpeedMPerS,
      c.solarRotationRateRadPerS, c.innerRadiusM, c.outerRadiusM,
      c.globalCellSizeM, c.minimumCellSizeM, c.solarSurfaceCellSizeM,
      c.solarTransitionOuterRadiusM, c.solarExponent,
      c.tubeReferenceRadiusM, c.tubeRadiusAtReferenceM,
      c.tubeCenterCellSizeM, c.tubeExponent};
  for (double value : values)
    if (!std::isfinite(value)) return Error("configuration contains a non-finite value");
  if ((c.schemaVersion != 1 && c.schemaVersion != 2 && c.schemaVersion != 3) ||
      c.parkerPointCount < 2 ||
      c.parkerPointCount > 10000000ULL || c.parkerLengthM <= 0.0 ||
      c.solarWindSpeedMPerS <= 0.0 || c.solarRotationRateRadPerS < 0.0 ||
      c.innerRadiusM <= 0.0 || c.outerRadiusM <= c.innerRadiusM ||
      c.minimumCellSizeM <= 0.0 || c.globalCellSizeM < c.minimumCellSizeM ||
      c.maximumMeshLevel > 19 || c.solarExponent <= 0.0 ||
      c.tubeExponent <= 0.0) {
    return Error("Parker, domain, or global mesh values are outside their supported range");
  }
  if (c.schemaVersion >= 2 &&
      (!std::isfinite(c.timeStepS) || c.timeStepS <= 0.0 ||
       c.macroparticlesPerStep == 0 ||
       c.macroparticlesPerStep >
           static_cast<std::uint64_t>(std::numeric_limits<int>::max()) ||
       !std::isfinite(c.particleWeight) || c.particleWeight <= 0.0 ||
       c.meshTecplotFile.empty() || c.fieldLineTecplotFile.empty())) {
    return Error("time step, particle sampling, observer, or Tecplot output "
                 "configuration is invalid");
  }
  if (c.schemaVersion == 2 &&
      (!std::isfinite(c.observerHeliocentricRadiusM) ||
       c.observerHeliocentricRadiusM < c.innerRadiusM ||
       c.observerHeliocentricRadiusM > c.outerRadiusM))
    return Error("schema version 2 observer radius is outside the domain");
  if (c.schemaVersion >= 3) {
    if (c.observers.empty() || c.dataTecplotFile.empty())
      return Error("schema version 3 requires observers and AMPS data output");
    std::set<std::string> observerIds;
    for (const ObserverConfiguration& observer : c.observers) {
      if (observer.id.empty() || observer.id.size() > 48 ||
          !observerIds.insert(observer.id).second ||
          !std::isfinite(observer.heliocentricRadiusM) ||
          observer.heliocentricRadiusM < c.innerRadiusM ||
          observer.heliocentricRadiusM > c.outerRadiusM ||
          !std::isfinite(observer.minimumEnergyJ) ||
          !std::isfinite(observer.maximumEnergyJ) ||
          observer.minimumEnergyJ <= 0.0 ||
          observer.maximumEnergyJ <= observer.minimumEnergyJ ||
          observer.energyChannels == 0 || observer.energyChannels > 1000000U ||
          observer.pitchAngleBins == 0 || observer.pitchAngleBins > 1000000U)
        return Error("one-dimensional observer identity, location, or spectrum is invalid");
    }
  }
  const double sourceRadius = Norm(Subtract(c.parkerInitialPointM, c.parkerOriginM));
  if (std::fabs(sourceRadius - c.innerRadiusM) > 1.0e-10 * c.innerRadiusM)
    return Error("Parker initial point must lie on domain.inner_radius_m");
  if (c.solarRefinementEnabled &&
      (c.solarSurfaceCellSizeM < c.minimumCellSizeM ||
       c.solarSurfaceCellSizeM > c.globalCellSizeM ||
       c.solarTransitionOuterRadiusM <= c.innerRadiusM ||
       c.solarTransitionOuterRadiusM > c.outerRadiusM)) {
    return Error("near-Sun refinement sizes/radii are inconsistent");
  }
  if (c.tubeRefinementEnabled &&
      (c.tubeReferenceRadiusM <= c.innerRadiusM ||
       c.tubeRadiusAtReferenceM <= 0.0 ||
       c.tubeCenterCellSizeM < c.minimumCellSizeM ||
       c.tubeCenterCellSizeM > c.globalCellSizeM)) {
    return Error("magnetic-tube refinement sizes/radii are inconsistent");
  }
  return Transport::Status::Ok();
}

Transport::Status ParseText(const std::string& text, Configuration* result) {
  if (result == NULL) return Error("configuration output is null");
  Configuration candidate;
  std::istringstream input(text);
  std::set<std::string> sections;
  std::set<std::string> assigned;
  std::string section;
  std::string line;
  std::size_t lineNumber = 0;
  const std::set<std::string> known = {
      "run", "parker_spiral", "domain", "mesh", "mesh.solar",
      "mesh.tube", "background.parker", "injection", "species",
      "observer", "output", "swcme"};
  while (std::getline(input, line)) {
    ++lineNumber;
    const std::size_t comment = line.find('#');
    if (comment != std::string::npos) line.erase(comment);
    line = Trim(line);
    if (line.empty()) continue;
    if (line.front() == '[' && line.back() == ']') {
      section = Lower(Trim(line.substr(1, line.size() - 2)));
      const bool namedObserver = section.rfind("observer.", 0) == 0 &&
          section.size() > std::string("observer.").size();
      if (known.count(section) == 0 && !namedObserver)
        return Error("unknown section at line " + std::to_string(lineNumber));
      if (!sections.insert(section).second)
        return Error("duplicate section '" + section + "'");
      continue;
    }
    const std::size_t separator = line.find('=');
    if (section.empty() || separator == std::string::npos)
      return Error("expected section key=value at line " + std::to_string(lineNumber));
    const std::string key = Lower(Trim(line.substr(0, separator)));
    const std::string value = Trim(line.substr(separator + 1));
    const std::string qualified = section + "." + key;
    if (key.empty() || value.empty() || !assigned.insert(qualified).second)
      return Error("empty or duplicate key '" + qualified + "'");
    // [swcme] is a transport envelope.  Model keys commonly contain dots and
    // are intentionally not duplicated in this mesh parser; main.cpp forwards
    // them to the canonical model-owned resolver, which rejects unknown keys,
    // invalid units, duplicates, and inconsistent combinations.
    if (section == "swcme") {
      SwcmeAssignment assignment;
      assignment.key = key;
      assignment.value = value;
      assignment.line = lineNumber;
      candidate.swcmeAssignments.push_back(assignment);
      continue;
    }
    const Transport::Status status = Apply(section, key, value, &candidate);
    if (!status.ok())
      return Error("line " + std::to_string(lineNumber) + ": " + status.message);
  }
  const std::set<std::string> versionOneSections = {
      "run", "parker_spiral", "domain", "mesh", "mesh.solar",
      "mesh.tube", "background.parker"};
  const std::set<std::string> versionTwoSections = known;
  const std::set<std::string> versionThreeSections = {
      "run", "parker_spiral", "domain", "mesh", "mesh.solar",
      "mesh.tube", "background.parker", "injection", "species",
      "output", "swcme"};
  const std::set<std::string>& requiredSections =
      candidate.schemaVersion >= 3 ? versionThreeSections :
      (candidate.schemaVersion >= 2 ? versionTwoSections : versionOneSections);
  for (const std::string& required : requiredSections)
    if (sections.count(required) == 0)
      return Error("missing required section '[" + required + "]'");
  // Every supported field is mandatory. This fail-closed count is paired with
  // duplicate/unknown rejection above, so a future key cannot silently inherit
  // a zero/default value without updating the schema and its tests.
  if (candidate.schemaVersion == 1 && assigned.size() != 28)
    return Error("configuration must assign all 28 version-1 keys exactly once");
  if (candidate.schemaVersion >= 2) {
    // There are 34 application-owned scalar/string keys in version 2.  SWCME
    // keys are counted separately because their names belong to the provider.
    const std::size_t applicationKeyCount =
        assigned.size() - candidate.swcmeAssignments.size();
    if (candidate.schemaVersion == 2 && applicationKeyCount != 34)
      return Error("configuration must assign all 34 version-2 application "
                   "keys exactly once");
    if (candidate.schemaVersion >= 3) {
      if (sections.count("observer") != 0)
        return Error("schema version 3 uses repeatable [observer.ID] sections");
      if (candidate.observers.empty())
        return Error("schema version 3 requires at least one [observer.ID]");
      const std::size_t expected = 34 + 6 * candidate.observers.size();
      if (applicationKeyCount != expected)
        return Error("schema version 3 requires every observer spectrum field "
                     "and all 34 non-observer application keys");
      const char* observerFields[] = {
          "heliocentric_radius_m", "minimum_energy_j", "maximum_energy_j",
          "energy_channels", "energy_spacing", "pitch_angle_bins"};
      for (const ObserverConfiguration& observer : candidate.observers)
        for (const char* key : observerFields)
          if (assigned.count("observer." + observer.id + "." + key) == 0)
            return Error("missing required observer key 'observer." +
                         observer.id + "." + key + "'");
    }
    if (candidate.swcmeAssignments.empty())
      return Error("schema version 2 requires complete [swcme] assignments");
    std::set<std::string> swcmeKeys;
    std::string kinematics;
    for (const SwcmeAssignment& assignment : candidate.swcmeAssignments) {
      swcmeKeys.insert(assignment.key);
      if (assignment.key == "cme.kinematics")
        kinematics = Lower(Trim(assignment.value));
      if (RequiredSwcmeKeys().count(assignment.key) == 0 &&
          assignment.key != "cme.data_times" &&
          assignment.key != "cme.data_radii")
        return Error("unknown SWCME configuration key '" + assignment.key +
                     "'");
    }
    for (const std::string& required : RequiredSwcmeKeys())
      if (swcmeKeys.count(required) == 0)
        return Error("schema version 2 is missing required SWCME key '" +
                     required + "'");
    const bool hasDataTimes = swcmeKeys.count("cme.data_times") != 0;
    const bool hasDataRadii = swcmeKeys.count("cme.data_radii") != 0;
    if (kinematics == "data_driven") {
      if (!hasDataTimes || !hasDataRadii)
        return Error("data_driven CME kinematics requires cme.data_times and "
                     "cme.data_radii");
    }
    else if (hasDataTimes || hasDataRadii) {
      return Error("CME data knots are legal only for data_driven kinematics");
    }
  }
  const Transport::Status valid = Validate(candidate);
  if (!valid.ok()) return valid;
  *result = candidate;
  return Transport::Status::Ok();
}

Transport::Status LoadFile(const std::string& path, Configuration* result) {
  std::ifstream input(path.c_str());
  if (!input.good()) return Error("cannot open initialization input file '" + path + "'");
  std::ostringstream text;
  text << input.rdbuf();
  if (!input.eof() && input.fail())
    return Error("failed while reading initialization input file '" + path + "'");
  return ParseText(text.str(), result);
}

Transport::Status ApplyOutputDirectoryOverride(
    const std::string& directory, Configuration* configuration) {
  if (configuration == NULL) return Error("initialization configuration is null");
  if (directory.empty()) return Error("initialization output directory is empty");

  // Preserve the reviewed product names and remove only their former parent.
  // Both separator spellings are recognized because input decks are commonly
  // copied between workstation and HPC filesystems before a production run.
  const auto leafName = [](const std::string& path) {
    const std::size_t separator = path.find_last_of("/\\");
    return separator == std::string::npos ? path : path.substr(separator + 1);
  };
  const std::string meshLeaf = leafName(configuration->meshTecplotFile);
  const std::string lineLeaf = leafName(configuration->fieldLineTecplotFile);
  const std::string dataLeaf = leafName(configuration->dataTecplotFile);
  if (meshLeaf.empty() || lineLeaf.empty() ||
      (configuration->schemaVersion >= 3 && dataLeaf.empty()) ||
      meshLeaf == "." || meshLeaf == ".." || lineLeaf == "." ||
      lineLeaf == ".." || dataLeaf == "." || dataLeaf == "..") {
    return Error("initialization Tecplot paths must end in file names before "
                 "the output-directory override can be applied");
  }
  const std::string separator = directory.back() == '/' ? "" : "/";
  configuration->meshTecplotFile = directory + separator + meshLeaf;
  configuration->fieldLineTecplotFile = directory + separator + lineLeaf;
  if (configuration->schemaVersion >= 3)
    configuration->dataTecplotFile = directory + separator + dataLeaf;
  return Validate(*configuration);
}

Transport::Status Install(const Configuration& configuration) {
  if (gActive.get() != NULL) return Error("initialization configuration is already installed");
  const Transport::Status valid = Validate(configuration);
  if (!valid.ok()) return valid;
  gActive.reset(new Configuration(configuration));
  return Transport::Status::Ok();
}

bool HasActive() { return gActive.get() != NULL; }
const Configuration& Active() { return *gActive; }

std::string Fingerprint(const Configuration& c) {
  std::ostringstream canonical;
  canonical << std::setprecision(17) << std::scientific
      << "srcsep-initialization-v" << c.schemaVersion;

  // Schema 1 fingerprints were already published before the complete
  // initialization contract existed.  Append version-2-only state only for a
  // version-2 document so restarting or comparing a legacy campaign retains
  // its exact historical identity.
  if (c.schemaVersion >= 2) {
    canonical << ";dt=" << c.timeStepS
        << ";macro_per_step=" << c.macroparticlesPerStep
        << ";particle_weight=" << c.particleWeight
        << ";mesh_output=" << c.meshTecplotFile
        << ";line_output=" << c.fieldLineTecplotFile;
    if (c.schemaVersion == 2)
      canonical << ";observer_radius=" << c.observerHeliocentricRadiusM;
    if (c.schemaVersion >= 3) {
      canonical << ";data_output=" << c.dataTecplotFile;
      for (const ObserverConfiguration& observer : c.observers)
        canonical << ";observer=" << observer.id << ','
            << observer.heliocentricRadiusM << ',' << observer.minimumEnergyJ
            << ',' << observer.maximumEnergyJ << ',' << observer.energyChannels
            << ',' << Name(observer.energySpacing) << ','
            << observer.pitchAngleBins;
    }
  }

  canonical << ";origin=" << c.parkerOriginM.x << ','
      << c.parkerOriginM.y << ',' << c.parkerOriginM.z << ";initial="
      << c.parkerInitialPointM.x << ',' << c.parkerInitialPointM.y << ','
      << c.parkerInitialPointM.z << ";length=" << c.parkerLengthM
      << ";points=" << c.parkerPointCount << ";wind=" << c.solarWindSpeedMPerS
      << ";omega=" << c.solarRotationRateRadPerS << ";inner=" << c.innerRadiusM
      << ";outer=" << c.outerRadiusM << ";global=" << c.globalCellSizeM
      << ";minimum=" << c.minimumCellSizeM << ";level=" << c.maximumMeshLevel
      << ";solar=" << c.solarRefinementEnabled << ',' << c.solarSurfaceCellSizeM
      << ',' << c.solarTransitionOuterRadiusM << ',' << Name(c.solarProfile)
      << ',' << c.solarExponent << ";tube=" << c.tubeRefinementEnabled << ','
      << c.tubeReferenceRadiusM << ',' << c.tubeRadiusAtReferenceM << ','
      << Name(c.tubeRadiusMode) << ',' << c.tubeCenterCellSizeM << ','
      << Name(c.tubeProfile) << ',' << c.tubeExponent;

  if (c.schemaVersion >= 2) {
    // INI field order has no physical meaning.  Sort a private copy so two
    // complete decks with the same SWCME assignments have the same startup
    // fingerprint.  The canonical SWCME resolver separately normalizes units
    // and produces the model-owned physical fingerprint.
    std::vector<SwcmeAssignment> assignments = c.swcmeAssignments;
    std::sort(assignments.begin(), assignments.end(),
        [](const SwcmeAssignment& left, const SwcmeAssignment& right) {
          if (left.key != right.key) return left.key < right.key;
          return left.value < right.value;
        });
    for (const SwcmeAssignment& assignment : assignments)
      canonical << ";swcme." << assignment.key << '=' << assignment.value;
  }
  // FNV-1a is used as a compact deterministic identity, not as a security
  // primitive. The canonical manifest above retains full scientific meaning.
  std::uint64_t hash = 1469598103934665603ULL;
  const std::string bytes = canonical.str();
  for (unsigned char byte : bytes) {
    hash ^= byte;
    hash *= 1099511628211ULL;
  }
  std::ostringstream result;
  result << std::hex << std::setfill('0') << std::setw(16) << hash;
  return result.str();
}

Transport::Status BuildParkerLine(const Configuration& c,
                                  std::vector<Vec3>* points) {
  if (points == NULL) return Error("Parker line output is null");
  const Transport::Status valid = Validate(c);
  if (!valid.ok()) return valid;
  std::vector<Vec3> candidate;
  candidate.reserve(static_cast<std::size_t>(c.parkerPointCount));
  Vec3 point = c.parkerInitialPointM;
  candidate.push_back(point);
  const double step = c.parkerLengthM /
      static_cast<double>(c.parkerPointCount - 1);
  for (std::uint64_t index = 1; index < c.parkerPointCount; ++index) {
    const Vec3 relative = Subtract(point, c.parkerOriginM);
    const Vec3 first = ParkerTangent(relative, c);
    const Vec3 midpoint = Add(relative, Scale(0.5 * step, first));
    const Vec3 tangent = ParkerTangent(midpoint, c);
    if (!(Norm(first) > 0.0) || !(Norm(tangent) > 0.0))
      return Error("Parker tangent vanished during line generation");
    point = Add(point, Scale(step, tangent));
    candidate.push_back(point);
  }
  points->swap(candidate);
  return Transport::Status::Ok();
}

Transport::Status WriteParkerLineTecplot(const Configuration& c,
                                         const std::string& path) {
  if (path.empty()) return Error("Parker-line Tecplot path is empty");
  std::vector<Vec3> points;
  const Transport::Status built = BuildParkerLine(c, &points);
  if (!built.ok()) return built;

  // Open only after all geometry succeeds.  A stream failure is returned to
  // the caller and is fatal during initialization; visualization is part of
  // the declared startup contract rather than a best-effort diagnostic.
  std::ofstream output(path.c_str(), std::ios::out | std::ios::trunc);
  if (!output.good())
    return Error("cannot open Parker-line Tecplot file '" + path + "'");
  output << "TITLE=\"srcSEP initialized Parker field line\"\n"
         << "VARIABLES=\"arc_length_m\",\"x_m\",\"y_m\",\"z_m\","
            "\"heliocentric_radius_m\",\"requested_cell_size_m\"\n"
         << "ZONE T=\"field-line\", I=" << points.size()
         << ", F=POINT\n" << std::scientific << std::setprecision(17);
  double arcLength = 0.0;
  for (std::size_t i = 0; i < points.size(); ++i) {
    if (i != 0) arcLength += Norm(Subtract(points[i], points[i - 1]));
    output << arcLength << ' ' << points[i].x << ' ' << points[i].y << ' '
           << points[i].z << ' '
           << Norm(Subtract(points[i], c.parkerOriginM)) << ' '
           << RequestedCellSizeM(points[i], c) << '\n';
  }
  output.flush();
  if (!output.good())
    return Error("failed while writing Parker-line Tecplot file '" + path + "'");
  return Transport::Status::Ok();
}

double TubeRadiusM(double radius, const Configuration& c) {
  if (!(radius > 0.0) || !(c.tubeReferenceRadiusM > 0.0))
    return std::numeric_limits<double>::quiet_NaN();
  return c.tubeRadiusMode == TubeRadiusMode::PhysicalConstant
      ? c.tubeRadiusAtReferenceM
      : c.tubeRadiusAtReferenceM * radius / c.tubeReferenceRadiusM;
}

double TubeDistanceM(const Vec3& position, const Configuration& c) {
  const Vec3 relative = Subtract(position, c.parkerOriginM);
  const double radius = Norm(relative);
  if (!(radius > 0.0)) return std::numeric_limits<double>::infinity();
  const Vec3 direction = Unit(relative);
  const Vec3 centre = ParkerDirectionAtRadius(radius, c);
  const double sine = Norm(Cross(direction, centre));
  const double cosine = Clamp(Dot(direction, centre), -1.0, 1.0);
  return radius * std::atan2(sine, cosine);
}

double RequestedCellSizeM(const Vec3& position, const Configuration& c) {
  const double radius = Norm(Subtract(position, c.parkerOriginM));
  double requested = c.globalCellSizeM;
  if (c.solarRefinementEnabled) {
    const double fraction = (radius - c.innerRadiusM) /
        (c.solarTransitionOuterRadiusM - c.innerRadiusM);
    requested = std::min(requested, c.solarSurfaceCellSizeM +
        Profile(fraction, c.solarProfile, c.solarExponent) *
        (c.globalCellSizeM - c.solarSurfaceCellSizeM));
  }
  if (c.tubeRefinementEnabled) {
    const double distance = TubeDistanceM(position, c);
    const double fraction = distance / TubeRadiusM(radius, c);
    requested = std::min(requested, c.tubeCenterCellSizeM +
        Profile(fraction, c.tubeProfile, c.tubeExponent) *
        (c.globalCellSizeM - c.tubeCenterCellSizeM));

    // The production AMR builder samples this point function on a Cartesian
    // lattice. The closest lattice point to a curve crossing a cell is within
    // half a cell diagonal, so max(h_tube,2*d) forces progressive capture of
    // the centreline even when the physical tube is initially sub-cell wide.
    constexpr double kStrictRefinement = 2.0 * (1.0 - 1.0e-12);
    const double capture = std::max(
        c.tubeCenterCellSizeM, kStrictRefinement * distance);
    requested = std::min(requested, capture);
  }
  return Clamp(requested, c.minimumCellSizeM, c.globalCellSizeM);
}

}  // namespace Initialization
}  // namespace SEP
