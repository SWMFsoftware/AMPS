#include "sep_initialization_validation.h"

#include "sep_initialization.h"

#include <algorithm>
#include <cstdio>
#include <cmath>
#include <fstream>
#include <sstream>

namespace {

std::string CompleteInput() {
  return R"SRCSEP(
[run]
schema_version = 2
time_step_s = 2

[injection]
macroparticles_per_step = 1024

[species]
particle_weight = 1e25

[observer]
heliocentric_radius_m = 1.0e11

[output]
mesh_tecplot_file = initialization-mesh.dat
field_line_tecplot_file = initialization-field-line.dat

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
solar_rotation_rate_rad_per_s = 2.86533e-6

[swcme]
preset = fast
ambient.wind_speed = 400 km/s
ambient.density_1au = 6 cm^-3
ambient.magnetic_field_1au = 5 nT
ambient.proton_temperature = 1.2e5 K
ambient.adiabatic_index = 1.6666666666666667
ambient.alpha_to_proton_ratio = 0.04
ambient.electron_temperature = 1.2e5 K
ambient.alpha_temperature = 4.8e5 K
ambient.thermodynamic_closure = proton_only
parker.radial_polarity = 1
parker.sin_theta = 1
parker.source_radius = 20 Rs
cme.kinematics = dbm
cme.launch_radius = 20 Rs
cme.launch_speed = 1900 km/s
cme.drag_coefficient = 8e-8 1/km
cme.extrapolation = outside_time
shock.region_mode = shock_only
shock.acceleration_mode = source
shock.relative_source_weight_per_area = 1
geometry.sheath_thickness_1au = 0.12 AU
geometry.ejecta_thickness_1au = 0.22 AU
smoothing.shock_width_1au = 0.010 AU
smoothing.leading_edge_width_1au = 0.020 AU
smoothing.trailing_edge_width_1au = 0.030 AU
sheath.ramp_power = 2
sheath.leading_edge_speed_factor = 1.12
ejecta.density_factor = 0.50
ejecta.speed_factor = 0.80
event.launch_epoch = 0 s
event.valid_from = 0 s
event.valid_until = 604800 s
source.particle_mass = 1 mp
source.charge_number = 1
source.energy_min = 0.1 MeV
source.energy_max = 500 MeV
source.reference_energy = 10 MeV
source.injection_efficiency = 3.4e-4
source.normalization = relative_only
source.reference_intensity_si = 1 si
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
      configuration.timeStepS != 2.0 ||
      configuration.macroparticlesPerStep != 1024 ||
      configuration.swcmeAssignments.size() != 41 ||
      configuration.parkerInitialPointM.x != configuration.innerRadiusM ||
      SEP::Initialization::Fingerprint(configuration).size() != 16) {
    return Fail("complete SI initialization input did not parse/fingerprint exactly");
  }

  // Assignment order is presentation, not physics.  Guard the startup
  // identity against accidental dependence on the order of [swcme] records.
  SEP::Initialization::Configuration reordered = configuration;
  std::reverse(reordered.swcmeAssignments.begin(),
               reordered.swcmeAssignments.end());
  if (SEP::Initialization::Fingerprint(configuration) !=
      SEP::Initialization::Fingerprint(reordered))
    return Fail("schema-2 fingerprint depends on SWCME assignment order");

  // Version-2 state must never contaminate a legacy identity.  This catches
  // regressions in the published schema-1 fingerprint serialization without
  // coupling the test to any particular version-2 numerical/source choice.
  SEP::Initialization::Configuration legacy = configuration;
  legacy.schemaVersion = 1;
  SEP::Initialization::Configuration legacyWithV2Noise = legacy;
  legacyWithV2Noise.timeStepS = 123.0;
  legacyWithV2Noise.macroparticlesPerStep = 9876;
  legacyWithV2Noise.particleWeight = 4.5e27;
  legacyWithV2Noise.observerHeliocentricRadiusM *= 0.75;
  legacyWithV2Noise.meshTecplotFile = "different-mesh.dat";
  legacyWithV2Noise.fieldLineTecplotFile = "different-line.dat";
  legacyWithV2Noise.swcmeAssignments.clear();
  if (SEP::Initialization::Fingerprint(legacy) !=
      SEP::Initialization::Fingerprint(legacyWithV2Noise))
    return Fail("schema-1 fingerprint includes schema-2-only state");

  std::string missing = CompleteInput();
  const std::string required = "point_count = 401\n";
  missing.erase(missing.find(required), required.size());
  if (SEP::Initialization::ParseText(missing, &configuration).ok())
    return Fail("parser accepted a missing required field");
  std::string unknown = CompleteInput();
  unknown += "\n[unknown]\nvalue=1\n";
  if (SEP::Initialization::ParseText(unknown, &configuration).ok())
    return Fail("parser accepted an unknown section");
  std::string missingSwcme = CompleteInput();
  const std::string swcmeRequired = "cme.launch_speed = 1900 km/s\n";
  missingSwcme.erase(missingSwcme.find(swcmeRequired), swcmeRequired.size());
  if (SEP::Initialization::ParseText(missingSwcme, &configuration).ok())
    return Fail("parser accepted an incomplete SWCME model setting");
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

SEP::Testing::Result RunINIT03() {
  SEP::Initialization::Configuration configuration;
  if (!SEP::Initialization::ParseText(CompleteInput(), &configuration).ok())
    return Fail("Tecplot fixture did not parse");
  const std::string path = "initialization-field-line-test.dat";
  const SEP::Transport::Status written =
      SEP::Initialization::WriteParkerLineTecplot(configuration, path);
  std::ifstream input(path.c_str());
  std::ostringstream payload;
  payload << input.rdbuf();
  const bool readOk = input.is_open() && !input.bad();
  std::remove(path.c_str());
  const std::string text = payload.str();
  if (!written.ok() || !readOk ||
      text.find("VARIABLES=\"arc_length_m\"") == std::string::npos ||
      text.find("I=401, F=POINT") == std::string::npos)
    return Fail("field-line Tecplot writer did not emit the declared ordered zone");
  return Pass("initialized finite field line is emitted as a unit-qualified Tecplot zone");
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
           "Finite line and radial/tube refinement identities.", RunINIT02),
      make("INIT03", "Initialization Tecplot output",
           "Validated field-line geometry is serialized after initialization.",
           RunINIT03)};
}
