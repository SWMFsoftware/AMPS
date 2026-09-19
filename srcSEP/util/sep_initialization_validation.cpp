#include "sep_initialization_validation.h"

#include "sep_initialization.h"

#include <algorithm>
#include <cmath>
#include <sstream>

namespace {

std::string CompleteInput() {
  return R"SRCSEP(
[run]
schema_version = 1

[parker_spiral]
origin_x_m = 0
origin_y_m = 0
origin_z_m = 0
initial_x_m = 1.3914e10
initial_y_m = 0
initial_z_m = 0
length_m = 2.0e11
point_count = 401

[domain]
inner_radius_m = 1.3914e10
outer_radius_m = 1.495978707e11

[mesh]
global_cell_size_m = 3.7399467675e10
minimum_cell_size_m = 1.495978707e9
maximum_level = 7

[mesh.solar]
enabled = true
surface_cell_size_m = 1.495978707e9
transition_outer_radius_m = 3.7399467675e10
profile = smoothstep
exponent = 1

[mesh.tube]
enabled = true
reference_radius_m = 1.495978707e11
radius_at_reference_m = 4.487936121e9
radius_mode = constant-angular-width
center_cell_size_m = 1.495978707e9
transverse_profile = smoothstep
transverse_exponent = 1

[background.parker]
solar_wind_speed_m_per_s = 400000
solar_rotation_rate_rad_per_s = 2.865e-6
)SRCSEP";
}

SEP::Testing::Result Pass(const std::string& message) {
  SEP::Testing::Result result;
  result.status = SEP::Testing::Status::Pass;
  result.message = message;
  result.metrics.push_back({"assertion_failures", 0.0, 0.0, "<=", "count"});
  return result;
}

SEP::Testing::Result Fail(const std::string& message) {
  SEP::Testing::Result result;
  result.status = SEP::Testing::Status::Fail;
  result.message = message;
  result.metrics.push_back({"assertion_failures", 1.0, 0.0, "<=", "count"});
  return result;
}

SEP::Testing::Result RunINIT01() {
  SEP::Initialization::Configuration configuration;
  if (!SEP::Initialization::ParseText(CompleteInput(), &configuration).ok() ||
      configuration.parkerPointCount != 401 ||
      configuration.maximumMeshLevel != 7 ||
      configuration.parkerInitialPointM.x != configuration.innerRadiusM ||
      SEP::Initialization::Fingerprint(configuration).size() != 16) {
    return Fail("complete SI initialization input did not parse/fingerprint exactly");
  }
  std::string missing = CompleteInput();
  const std::string required = "point_count = 401\n";
  missing.erase(missing.find(required), required.size());
  if (SEP::Initialization::ParseText(missing, &configuration).ok())
    return Fail("parser accepted a missing required field");
  std::string unknown = CompleteInput();
  unknown += "\n[unknown]\nvalue=1\n";
  if (SEP::Initialization::ParseText(unknown, &configuration).ok())
    return Fail("parser accepted an unknown section");
  return Pass("complete SI input parses while missing and unknown fields fail closed");
}

SEP::Testing::Result RunINIT02() {
  SEP::Initialization::Configuration configuration;
  if (!SEP::Initialization::ParseText(CompleteInput(), &configuration).ok())
    return Fail("geometry fixture did not parse");
  std::vector<SEP::Initialization::Vec3> points;
  if (!SEP::Initialization::BuildParkerLine(configuration, &points).ok() ||
      points.size() != configuration.parkerPointCount)
    return Fail("Parker generator did not return the requested point count");
  double length = 0.0;
  for (std::size_t i = 1; i < points.size(); ++i) {
    const double dx = points[i].x - points[i - 1].x;
    const double dy = points[i].y - points[i - 1].y;
    const double dz = points[i].z - points[i - 1].z;
    length += std::sqrt(dx*dx + dy*dy + dz*dz);
  }
  const double surface = SEP::Initialization::RequestedCellSizeM(
      configuration.parkerInitialPointM, configuration);
  const double tubeReference = SEP::Initialization::TubeRadiusM(
      configuration.tubeReferenceRadiusM, configuration);
  if (std::fabs(length - configuration.parkerLengthM) >
          1.0e-12 * configuration.parkerLengthM ||
      surface != configuration.minimumCellSizeM ||
      tubeReference != configuration.tubeRadiusAtReferenceM) {
    return Fail("arc length or closed-form refinement identity changed");
  }
  SEP::Testing::Result result = Pass(
      "finite Parker line preserves count/arc length and mesh profiles satisfy their exact endpoints");
  result.metrics.push_back({"polyline_length_relative_error",
      std::fabs(length - configuration.parkerLengthM) /
          configuration.parkerLengthM, 1.0e-12, "<=", ""});
  return result;
}

}  // namespace

std::vector<SEP::Testing::Descriptor>
SEP::Testing::InitializationDescriptors() {
  auto make = [](const char* id, const char* name, const char* description,
                 SEP::Testing::TestCallback callback) {
    SEP::Testing::Descriptor descriptor;
    descriptor.id = id;
    descriptor.name = name;
    descriptor.group = "initialization";
    descriptor.description = description;
    descriptor.initialization = SEP::Testing::InitializationLevel::None;
    descriptor.supportedBuildModes = "standalone-no-AMPS and linked AMPS";
    descriptor.runtime = SEP::Testing::RuntimeClass::Routine;
    descriptor.seedPolicy = "deterministic; no RNG";
    descriptor.stateIsolation = "test-local configuration; no active install";
    descriptor.callback = callback;
    return descriptor;
  };
  return {
      make("INIT01", "Initialization schema",
           "Strict complete-input parsing and fingerprinting.", RunINIT01),
      make("INIT02", "Parker mesh initialization",
           "Finite line and radial/tube refinement identities.", RunINIT02)};
}
