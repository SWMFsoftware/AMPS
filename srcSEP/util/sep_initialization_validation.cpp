#include "sep_initialization_validation.h"

#include "sep_initialization.h"

#include <algorithm>
#include <cstdio>
#include <cmath>
#include <fstream>
#include <sstream>

namespace {

// Complete srcSEP input in the shared SEP content schema 6 (src/INPUT_FILE.md
// section 5).  Values are independent of the shipped examples: the test
// asserts them literally below.
std::string CompleteInput() {
  return R"SRCSEP(
#begin run
  time_step_s = 2
#end run

#begin application: sep
  schema_version = 6
  transport = focused-diffusion

  #begin turbulence
    source = self-consistent-integrated
    representation = integrated
    coupling = streaming-energy-exchange
    advection = true
    reflection = true
    cascade = true
    shock_injection = true
    inner_boundary = fixed-reservoir
    inner_boundary_value_si = 0
    outer_boundary = fixed-reservoir
    outer_boundary_value_si = 0
    reflection_coefficient = 0.6
    cascade_coefficient = 0.8
    perpendicular_correlation_length_m = 1e7
    electron_heating_fraction = 0.3
    spectral_k_min_per_m = 1e-10
    spectral_k_max_per_m = 1e-2
    spectral_bins = 128
    cfl_safety = 0.8
    operator_accuracy_safety = 0.25
    maximum_source_fraction = 0.25
    maximum_cascade_fraction = 0.2
    minimum_substep_s = 1e-12
    maximum_substeps = 1000000
    conservation_relative_tolerance = 1e-12
    shock_source_efficiency = 0.02
    shock_source_plus_fraction = 0.5
    diagnostic_cadence_steps = 100
  #end turbulence
#end application

#begin mesh
  global_cell_size_m = 3.7399467675e10
  minimum_cell_size_m = 1.495978707e9
  maximum_level = 7

  #begin field_line
    origin_x_m = 0
    origin_y_m = 0
    origin_z_m = 0
    initial_x_m = 1.3914e10
    initial_y_m = 0
    initial_z_m = 0
    length_m = 2.0e11
    point_count = 401
  #end field_line

  #begin domain
    inner_radius_m = 1.3914e10
    outer_radius_m = 1.495978707e11
  #end domain

  #begin refinement: solar
    enabled = true
    surface_cell_size_m = 1.495978707e9
    transition_outer_radius_m = 3.7399467675e10
    profile = smoothstep
    exponent = 1
  #end refinement

  #begin refinement: tube
    enabled = true
    reference_radius_m = 1.495978707e11
    radius_at_reference_m = 4.487936121e9
    radius_mode = constant-angular-width
    center_cell_size_m = 1.495978707e9
    transverse_profile = smoothstep
    transverse_exponent = 1
  #end refinement
#end mesh

#begin physics: background
  provider = swcme
#end physics

#begin physics: shock
  authority = swcme
#end physics

#begin physics: swcme
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
#end physics

#begin source: shock
  enabled = true
  samples_per_step = 1024

  minimum_energy_j = 1.602176634e-14
  maximum_energy_j = 8.01088317e-11
  injection_efficiency = 3.4e-4
  spectrum_reference_energy_j = 1.602176634e-12
  spectrum_particle_mass_kg = 1.67262192595e-27
  spectrum_charge_number = 1
  spectrum_normalization = relative-only
  spectrum_reference_intensity_per_m2_s_sr_j = 1
#end source

#begin species
  macroparticle_weight = 1e25
#end species

#begin observer: earth
  kind = field-line-radius
  heliocentric_radius_m = 1.0e11
  minimum_energy_j = 1.602176634e-15
  maximum_energy_j = 8.01088317e-11
  energy_bins = 70
  energy_spacing = logarithmic
  pitch_angle_bins = 24
#end observer

#begin output
  initialization_mesh_tecplot_file = initialization-mesh.dat
  initialization_field_line_tecplot_file = initialization-field-line.dat
  initialization_data_tecplot_file = initialization-data.dat
#end output
)SRCSEP";
}

// Replace the first occurrence of 'from' (which must exist) by 'to'.
std::string Edited(std::string text, const std::string& from,
                   const std::string& to) {
  const std::size_t at = text.find(from);
  if (at == std::string::npos) return std::string();
  return text.replace(at, from.size(), to);
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
  // identity against accidental dependence on the order of physics: swcme
  // records.
  SEP::Initialization::Configuration reordered = configuration;
  std::reverse(reordered.swcmeAssignments.begin(),
               reordered.swcmeAssignments.end());
  if (SEP::Initialization::Fingerprint(configuration) !=
      SEP::Initialization::Fingerprint(reordered))
    return Fail("startup fingerprint depends on SWCME assignment order");

  // Decision D4: the Parker winding takes the SWCME ambient wind (400 km/s)
  // and the fixed SWCME1D Carrington rate; the SWCME source.* layer is
  // synthesized from source: shock in exact SI decimal text.
  bool synthesizedEnergy = false;
  for (const SEP::Initialization::SwcmeAssignment& assignment :
       configuration.swcmeAssignments)
    if (assignment.key == "source.energy_min" &&
        assignment.value == "1.602176634e-14 J") synthesizedEnergy = true;
  if (configuration.solarWindSpeedMPerS != 4.0e5 ||
      configuration.solarRotationRateRadPerS != 2.86533e-6 ||
      !synthesizedEnergy || !configuration.particleMoverSpecified ||
      configuration.particleMover !=
          SEP::Mover::ProductionMover::FocusedTransportDiffusion ||
      !configuration.turbulenceSpecified ||
      configuration.turbulence.spectralBins != 128 ||
      configuration.shockSourceEfficiency != 0.02 ||
      configuration.diagnosticCadenceSteps != 100)
    return Fail("schema-6 SWCME ownership, mover or turbulence values were "
                "not installed exactly");

  // Negative controls: each must fail closed with its own diagnosis.
  struct Rejection { const char* label; std::string text; const char* expect; };
  const std::string base = CompleteInput();
  const Rejection rejections[] = {
      {"retired schema", Edited(base, "schema_version = 6", "schema_version = 3"),
       "retired"},
      {"srcSEP3D-only block (D3)", Edited(base, "#end mesh",
           "  #begin memory\n    particle_bytes = 160\n  #end memory\n#end mesh"),
       "not supported by srcSEP"},
      {"srcSEP3D-only key (D3)", Edited(base, "  time_step_s = 2",
           "  time_step_s = 2\n  campaign_seed = 1"),
       "not supported by srcSEP"},
      {"srcSEP3D-only provider value (D3)",
       Edited(base, "provider = swcme", "provider = analytic-parker"),
       "not supported by srcSEP"},
      {"srcSEP3D-only observer kind (D3)",
       Edited(base, "kind = field-line-radius", "kind = fixed-cartesian"),
       "not supported by srcSEP"},
      {"SWCME source key in physics: swcme (D4)",
       Edited(base, "  event.valid_until = 604800 s",
              "  event.valid_until = 604800 s\n  source.energy_min = 0.1 MeV"),
       "written in source: shock"},
      {"Parker record (D4: SWCME owns the ambient)",
       Edited(base, "  provider = swcme",
              "  provider = swcme\n  #begin parker\n    "
              "solar_wind_speed_m_per_s = 400000\n  #end parker"),
       "not supported by srcSEP"},
      {"3-D SWCME geometry key (D3)",
       Edited(base, "  event.valid_until = 604800 s",
              "  event.valid_until = 604800 s\n  geometry.shape = sphere"),
       "not supported by srcSEP"},
      {"3-D SWCME rotation different from SWCME1D",
       Edited(base, "  event.valid_until = 604800 s",
              "  event.valid_until = 604800 s\n"
              "  parker.solar_rotation_rate_rad_per_s = 2.865e-6"),
       "not supported by srcSEP"},
      {"retired particle_mover key",
       Edited(base, "  transport = focused-diffusion",
              "  transport = focused-diffusion\n  particle_mover = parker"),
       "unknown configuration key"},
      {"incomplete turbulence block",
       Edited(base, "    spectral_bins = 128\n", ""), "spectral_bins"}};
  for (const Rejection& rejection : rejections) {
    SEP::Initialization::Configuration rejected;
    const SEP::Transport::Status status =
        SEP::Initialization::ParseText(rejection.text, &rejected);
    if (rejection.text.empty() || status.ok() ||
        status.message.find(rejection.expect) == std::string::npos)
      return Fail(std::string("negative control '") + rejection.label +
                  "' was not rejected with '" + rejection.expect + "': " +
                  status.message);
  }
  // The 3-D SWCME rotation key is honoured when it states the fixed SWCME1D
  // rate (D3: reject only what srcSEP cannot honour).
  SEP::Initialization::Configuration honoured;
  if (!SEP::Initialization::ParseText(
          Edited(base, "  event.valid_until = 604800 s",
                 "  event.valid_until = 604800 s\n"
                 "  parker.solar_rotation_rate_rad_per_s = 2.86533e-6"),
          &honoured).ok() ||
      SEP::Initialization::Fingerprint(honoured) !=
          SEP::Initialization::Fingerprint(configuration))
    return Fail("the SWCME1D rotation rate written as the 3-D key was not "
                "accepted as an identical setting");

  std::string missing = CompleteInput();
  const std::string required = "point_count = 401\n";
  missing.erase(missing.find(required), required.size());
  if (SEP::Initialization::ParseText(missing, &configuration).ok())
    return Fail("parser accepted a missing required field");
  std::string unknown = CompleteInput();
  unknown += "\n#begin unknown\nvalue=1\n#end unknown\n";
  if (SEP::Initialization::ParseText(unknown, &configuration).ok())
    return Fail("parser accepted an unknown block");
  // The retired INI syntax is not part of the generic grammar.
  if (SEP::Initialization::ParseText("[run]\ntime_step_s = 1\n",
                                     &configuration).ok())
    return Fail("parser accepted a retired INI section header");
  std::string missingSwcme = CompleteInput();
  const std::string swcmeRequired = "cme.launch_speed = 1900 km/s\n";
  missingSwcme.erase(missingSwcme.find(swcmeRequired), swcmeRequired.size());
  if (SEP::Initialization::ParseText(missingSwcme, &configuration).ok())
    return Fail("parser accepted an incomplete SWCME model setting");
  return Pass("complete schema-6 input parses; retired, unsupported, "
              "duplicated, missing and unknown settings fail closed");
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
