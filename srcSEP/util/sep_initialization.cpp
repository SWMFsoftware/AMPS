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

Transport::Status Apply(const std::string& section, const std::string& key,
                        const std::string& value, Configuration* c) {
  const std::string field = section + "." + key;
  auto invalid = [&]() { return Error("invalid value for '" + field + "': " + value); };
  if (field == "run.schema_version") {
    std::uint64_t parsed = 0;
    if (!ParseUnsigned64(value, &parsed) || parsed != 1) return invalid();
    c->schemaVersion = static_cast<unsigned>(parsed);
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
    return Error("unknown configuration key '" + field + "'");
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
  if (c.schemaVersion != 1 || c.parkerPointCount < 2 ||
      c.parkerPointCount > 10000000ULL || c.parkerLengthM <= 0.0 ||
      c.solarWindSpeedMPerS <= 0.0 || c.solarRotationRateRadPerS < 0.0 ||
      c.innerRadiusM <= 0.0 || c.outerRadiusM <= c.innerRadiusM ||
      c.minimumCellSizeM <= 0.0 || c.globalCellSizeM < c.minimumCellSizeM ||
      c.maximumMeshLevel > 19 || c.solarExponent <= 0.0 ||
      c.tubeExponent <= 0.0) {
    return Error("Parker, domain, or global mesh values are outside their supported range");
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
      "mesh.tube", "background.parker"};
  while (std::getline(input, line)) {
    ++lineNumber;
    const std::size_t comment = line.find('#');
    if (comment != std::string::npos) line.erase(comment);
    line = Trim(line);
    if (line.empty()) continue;
    if (line.front() == '[' && line.back() == ']') {
      section = Lower(Trim(line.substr(1, line.size() - 2)));
      if (known.count(section) == 0)
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
    const Transport::Status status = Apply(section, key, value, &candidate);
    if (!status.ok())
      return Error("line " + std::to_string(lineNumber) + ": " + status.message);
  }
  for (const std::string& required : known)
    if (sections.count(required) == 0)
      return Error("missing required section '[" + required + "]'");
  // Every supported field is mandatory. This fail-closed count is paired with
  // duplicate/unknown rejection above, so a future key cannot silently inherit
  // a zero/default value without updating the schema and its tests.
  if (assigned.size() != 28)
    return Error("configuration must assign all 28 version-1 keys exactly once");
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
      << "srcsep-initialization-v1;origin=" << c.parkerOriginM.x << ','
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
    const double fraction = TubeDistanceM(position, c) / TubeRadiusM(radius, c);
    requested = std::min(requested, c.tubeCenterCellSizeM +
        Profile(fraction, c.tubeProfile, c.tubeExponent) *
        (c.globalCellSizeM - c.tubeCenterCellSizeM));
  }
  return Clamp(requested, c.minimumCellSizeM, c.globalCellSizeM);
}

}  // namespace Initialization
}  // namespace SEP
