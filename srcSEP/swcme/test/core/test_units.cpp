#include "test_framework.hpp"

#include <swcme1d.hpp>
#include <swcme3d.hpp>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>

namespace {

// CFG02 references are intentionally independent numerical definitions. The
// AU and nominal solar radius are exact adopted SWCME conversion constants;
// decimal SI prefixes and the hour are exact definitions. Pi is written to
// higher precision than double so degree conversion is checked only at the
// rounding imposed by the target floating-point type.
constexpr double REF_AU_M = 149597870700.0;
constexpr double REF_SOLAR_RADIUS_M = 695700000.0;
constexpr double REF_PI = 3.141592653589793238462643383279502884;
constexpr double REF_PROTON_MASS_KG = 1.67262192595e-27;
constexpr double REF_MU0_N_A2 = 1.25663706127e-6;

constexpr double KM_S_TO_M_S = 1000.0;
constexpr double CM3_TO_M3 = 1000000.0;
constexpr double NT_TO_T = 1.0e-9;
constexpr double HOUR_TO_S = 3600.0;
constexpr double INVERSE_KM_TO_INVERSE_M = 1.0e-3;

// Sixty-four ulps at unit scale is the framework's roundoff allowance for a
// short chain of double operations. Exact integer/power-of-ten paths use zero
// whenever the selected input and result have identical binary values.
constexpr double ROUNDOFF_REL_TOL =
    64.0 * std::numeric_limits<double>::epsilon();

// The 1-D writers serialize with "% .9e". Half a unit in the last printed
// decimal corresponds to 5e-10 relative error; the small epsilon margin only
// covers parsing roundoff and is not a physics tolerance.
constexpr double OUTPUT_TEXT_REL_TOL =
    5.0e-10 + 4.0 * std::numeric_limits<double>::epsilon();

struct Metrics {
  double maximum_relative_conversion_error = 0.0;
  double maximum_roundtrip_error = 0.0;
  double maximum_1d_3d_difference = 0.0;
  double alfven_speed_relative_error = 0.0;
  bool any_roundtrip_executed = false;
};

double relative_error(double actual, double reference) {
  const double absolute = std::abs(actual - reference);
  return reference == 0.0 ? absolute : absolute / std::abs(reference);
}

// Print and record a required forward conversion. The factor column is the
// exact/reference SI factor, while `converted` comes from a production model
// state, field evaluator, or exposed adopted conversion constant.
void check_conversion(swcme_test::Context& context, Metrics& metrics,
                      const std::string& quantity, double external,
                      const char* external_unit, double converted,
                      double reference, double factor, double tolerance,
                      bool exact) {
  const double absolute = std::abs(converted - reference);
  const double relative = relative_error(converted, reference);
  const bool pass = std::isfinite(converted) &&
                    (exact ? absolute == 0.0 : relative <= tolerance);
  metrics.maximum_relative_conversion_error =
      std::max(metrics.maximum_relative_conversion_error, relative);

  std::cout << std::left << std::setw(30) << quantity << std::right
            << " external=" << std::scientific << std::setprecision(12)
            << external << ' ' << external_unit
            << " si=" << converted
            << " ref=" << reference
            << " factor=" << factor
            << " abs_err=" << absolute
            << " rel_err=" << relative
            << " tol=" << (exact ? 0.0 : tolerance)
            << ' ' << (pass ? "PASS" : "FAIL") << '\n';
  context.record_result(pass);
}

// Separate pair checks make duplicated 1-D and 3-D conversion paths visible.
// They use the same roundoff policy as forward calculations and contribute to
// the maximum cross-model difference structured metric.
void check_model_consistency(swcme_test::Context& context, Metrics& metrics,
                             const std::string& quantity, double value_1d,
                             double value_3d) {
  const double difference = relative_error(value_1d, value_3d);
  const bool pass = std::isfinite(value_1d) && std::isfinite(value_3d) &&
                    difference <= ROUNDOFF_REL_TOL;
  metrics.maximum_1d_3d_difference =
      std::max(metrics.maximum_1d_3d_difference, difference);
  std::cout << std::left << std::setw(30) << quantity << std::right
            << " 1D=" << std::scientific << std::setprecision(12) << value_1d
            << " 3D=" << value_3d
            << " rel_diff=" << difference
            << " tol=" << ROUNDOFF_REL_TOL
            << ' ' << (pass ? "PASS" : "FAIL") << '\n';
  context.record_result(pass);
}

void print_skip(swcme_test::Context& context, const std::string& quantity,
                const std::string& reason) {
  std::cout << std::left << std::setw(30) << quantity
            << " SKIP — " << reason << '\n';
  context.record_skip();
}

void check_roundtrip(swcme_test::Context& context, Metrics& metrics,
                     const std::string& quantity, double external,
                     double recovered) {
  const double error = relative_error(recovered, external);
  const bool pass = std::isfinite(recovered) && error <= OUTPUT_TEXT_REL_TOL;
  metrics.any_roundtrip_executed = true;
  metrics.maximum_roundtrip_error =
      std::max(metrics.maximum_roundtrip_error, error);
  std::cout << std::left << std::setw(30) << quantity << std::right
            << " original=" << std::scientific << std::setprecision(12)
            << external << " recovered=" << recovered
            << " rel_err=" << error
            << " tol=" << OUTPUT_TEXT_REL_TOL
            << ' ' << (pass ? "PASS" : "FAIL") << '\n';
  context.record_result(pass);
}

// Exercise the actual inverse distance conversions in the 1-D profile writer.
// The temporary file is outside the repository and is removed immediately so
// validation never leaves generated fixtures or intentional failures behind.
bool read_1d_profile_units(double radius_m, double& radius_au,
                           double& radius_rs) {
  const char* path = "/tmp/test_swcme_cfg02_profile.dat";
  const swcme1d::Model model;
  const swcme1d::StepState state = model.prepare_step(0.0);
  double density = 0.0, velocity = 0.0, br = 0.0, bphi = 0.0;
  double bmag = 0.0, div_v = 0.0;
  model.evaluate_radii_with_B_div(state, &radius_m, &density, &velocity,
                                  &br, &bphi, &bmag, &div_v, 1);
  const bool written = model.write_tecplot_radial_profile(
      state, &radius_m, &density, &velocity, &br, &bphi, &bmag, &div_v,
      1, path);

  std::ifstream input(path);
  std::string header;
  for (int i = 0; i < 3 && std::getline(input, header); ++i) {
    // The first three lines are the documented Tecplot headers.
  }
  double parsed_radius_m = 0.0;
  const bool parsed = written &&
      static_cast<bool>(input >> parsed_radius_m >> radius_au >> radius_rs);
  input.close();
  std::remove(path);
  return parsed;
}

// Exercise the 1-D shock writer's production m/s-to-km/s output path. At t=0
// the chosen launch speed is exactly the external value before serialization.
bool read_1d_shock_velocity(double velocity_km_s,
                            double& recovered_velocity_km_s) {
  const char* path = "/tmp/test_swcme_cfg02_shock.dat";
  swcme1d::Params parameters;
  parameters.V_sw_kms = velocity_km_s;
  parameters.V0_sh_kms = velocity_km_s;
  const swcme1d::Model model(parameters);
  const bool written = model.write_tecplot_shock_vs_time(0.0, 1, path);

  std::ifstream input(path);
  std::string header;
  for (int i = 0; i < 3 && std::getline(input, header); ++i) {
    // The first three lines are the documented Tecplot headers.
  }
  double time_s = 0.0, shock_radius_rs = 0.0, compression = 0.0;
  const bool parsed = written && static_cast<bool>(
      input >> time_s >> shock_radius_rs >> recovered_velocity_km_s >> compression);
  input.close();
  std::remove(path);
  return parsed;
}

double velocity_1d(double velocity_km_s) {
  swcme1d::Params parameters;
  parameters.V_sw_kms = velocity_km_s;
  // A positive launch speed and drag keep unrelated DBM fields representable;
  // CFG02 reads only the converted ambient speed from the production state.
  parameters.V0_sh_kms = std::max(1500.0, velocity_km_s + 1.0);
  return swcme1d::Model(parameters).prepare_step(0.0).V_up_ms;
}

double velocity_3d(double velocity_km_s) {
  swcme3d::Params parameters;
  parameters.V_sw_kms = velocity_km_s;
  parameters.V0_sh_kms = std::max(1500.0, velocity_km_s + 1.0);
  return swcme3d::Model(parameters).prepare_step(0.0).V_sw_ms;
}

struct AmbientFields {
  double density_m3;
  double magnetic_field_T;
};

AmbientFields ambient_1d(double density_cm3, double magnetic_field_nT) {
  swcme1d::Params parameters;
  parameters.n1AU_cm3 = density_cm3;
  parameters.B1AU_nT = magnetic_field_nT;
  const swcme1d::Model model(parameters);
  const swcme1d::StepState state = model.prepare_step(0.0);
  const double radius[] = {swcme1d::AU};
  double density[] = {0.0}, velocity[] = {0.0};
  double br[] = {0.0}, bphi[] = {0.0}, bmag[] = {0.0};
  model.evaluate_radii_with_B_div(state, radius, density, velocity, br, bphi,
                                  bmag, nullptr, 1);
  return {density[0], bmag[0]};
}

AmbientFields ambient_3d(double density_cm3, double magnetic_field_nT) {
  swcme3d::Params parameters;
  parameters.n1AU_cm3 = density_cm3;
  parameters.B1AU_nT = magnetic_field_nT;
  const swcme3d::Model model(parameters);
  const swcme3d::StepState state = model.prepare_step(0.0);
  const double x[] = {swcme3d::AU}, y[] = {0.0}, z[] = {0.0};
  double density[] = {0.0}, vx[] = {0.0}, vy[] = {0.0}, vz[] = {0.0};
  double bx[] = {0.0}, by[] = {0.0}, bz[] = {0.0};
  model.evaluate_cartesian_with_B(state, x, y, z, density, vx, vy, vz,
                                  bx, by, bz, 1);
  return {density[0],
          std::sqrt(bx[0] * bx[0] + by[0] * by[0] + bz[0] * bz[0])};
}

double inferred_gamma_1d(double gamma_km_inverse) {
  swcme1d::Params parameters;
  parameters.V_sw_kms = 321.0;
  parameters.V0_sh_kms = 654.0;
  parameters.Gamma_kmInv = gamma_km_inverse;
  const swcme1d::Model model(parameters);
  const swcme1d::StepState state = model.prepare_step(10.0);
  const double u0 = 333000.0;
  const double u = state.V_sh_ms - 321000.0;
  // This algebraically inverts the independently documented DBM speed law;
  // it does not repeat the production km^-1 conversion expression.
  return (u0 / u - 1.0) / (u0 * 10.0);
}

double inferred_gamma_3d(double gamma_km_inverse) {
  swcme3d::Params parameters;
  parameters.V_sw_kms = 321.0;
  parameters.V0_sh_kms = 654.0;
  parameters.Gamma_kmInv = gamma_km_inverse;
  const swcme3d::Model model(parameters);
  const swcme3d::StepState state = model.prepare_step(10.0);
  const double u0 = 333000.0;
  const double u = state.V_sh_ms - 321000.0;
  return (u0 / u - 1.0) / (u0 * 10.0);
}

void check_forward_conversions(swcme_test::Context& context, Metrics& metrics) {
  std::cout << "CFG02 forward conversions\n";

  const double velocities[] = {400.0, 0.0, 1.0, 321.0, 1500.0, 3000.0};
  for (double external : velocities) {
    const double reference = external * KM_S_TO_M_S;
    const double converted_1d = velocity_1d(external);
    const double converted_3d = velocity_3d(external);
    check_conversion(context, metrics, "1D velocity", external, "km/s",
                     converted_1d, reference, KM_S_TO_M_S, 0.0, true);
    check_conversion(context, metrics, "3D velocity", external, "km/s",
                     converted_3d, reference, KM_S_TO_M_S, 0.0, true);
  }

  const double magnetic_fields[] = {5.0, 1.0, 7.0, 100.0, 0.1, 1000.0};
  for (double external : magnetic_fields) {
    const AmbientFields fields_1d = ambient_1d(5.0, external);
    const AmbientFields fields_3d = ambient_3d(5.0, external);
    const double reference = external * NT_TO_T;
    check_conversion(context, metrics, "1D magnetic field", external, "nT",
                     fields_1d.magnetic_field_T, reference, NT_TO_T,
                     ROUNDOFF_REL_TOL, false);
    check_conversion(context, metrics, "3D magnetic field", external, "nT",
                     fields_3d.magnetic_field_T, reference, NT_TO_T,
                     ROUNDOFF_REL_TOL, false);
  }

  const double densities[] = {5.0, 1.0, 8.0, 100.0, 0.01, 1.0e4};
  for (double external : densities) {
    const AmbientFields fields_1d = ambient_1d(external, 5.0);
    const AmbientFields fields_3d = ambient_3d(external, 5.0);
    const double reference = external * CM3_TO_M3;
    check_conversion(context, metrics, "1D number density", external, "cm^-3",
                     fields_1d.density_m3, reference, CM3_TO_M3,
                     ROUNDOFF_REL_TOL, false);
    check_conversion(context, metrics, "3D number density", external, "cm^-3",
                     fields_3d.density_m3, reference, CM3_TO_M3,
                     ROUNDOFF_REL_TOL, false);
  }

  const double distances_au[] = {1.0, 0.1, 0.5, 2.0, 5.0};
  for (double external : distances_au) {
    // SWCME exposes adopted AU constants rather than a callable conversion
    // function; multiplying by each model's public constant is its real path.
    const double reference = external * REF_AU_M;
    check_conversion(context, metrics, "1D astronomical distance", external,
                     "AU", external * swcme1d::AU, reference, REF_AU_M,
                     ROUNDOFF_REL_TOL, false);
    check_conversion(context, metrics, "3D astronomical distance", external,
                     "AU", external * swcme3d::AU, reference, REF_AU_M,
                     ROUNDOFF_REL_TOL, false);
  }

  const double au_in_solar_radii = REF_AU_M / REF_SOLAR_RADIUS_M;
  const double distances_rs[] = {1.0, 2.0, 20.0, au_in_solar_radii};
  for (double external : distances_rs) {
    const double reference = external * REF_SOLAR_RADIUS_M;
    check_conversion(context, metrics, "1D solar-radius distance", external,
                     "Rs", external * swcme1d::Rs, reference,
                     REF_SOLAR_RADIUS_M, ROUNDOFF_REL_TOL, false);
    check_conversion(context, metrics, "3D solar-radius distance", external,
                     "Rs", external * swcme3d::Rs, reference,
                     REF_SOLAR_RADIUS_M, ROUNDOFF_REL_TOL, false);
  }

  // Production APIs accept seconds directly; they expose no hours-to-seconds
  // configuration path. Required hour cases are visible SKIPs rather than a
  // test-side duplicate presented as production behavior.
  const double hours[] = {1.0, 0.5, 24.0, 72.0};
  for (double external : hours) {
    std::cout << "reference time conversion: " << external << " hr = "
              << external * HOUR_TO_S << " s\n";
    print_skip(context, "1D time conversion",
               "production API accepts seconds; hours are not exposed");
    print_skip(context, "3D time conversion",
               "production API accepts seconds; hours are not exposed");
  }

  // The 3-D default is the only production degree conversion. Its exposed
  // 40-degree result supplies the actual production factor for all requested
  // angle probes; 1-D accepts radians and has no degree-valued input.
  const double production_degree_factor =
      swcme3d::Params{}.half_width_rad / 40.0;
  const double degrees[] = {40.0, 90.0, 180.0};
  for (double external : degrees) {
    const double reference = external * REF_PI / 180.0;
    print_skip(context, "1D angle conversion",
               "production 1-D API accepts radians; degrees are not exposed");
    check_conversion(context, metrics, "3D angle", external, "deg",
                     external * production_degree_factor, reference,
                     REF_PI / 180.0, ROUNDOFF_REL_TOL, false);
  }

  const double gamma_external = 2.5e-5;
  check_conversion(context, metrics, "1D inverse length", gamma_external,
                   "km^-1", inferred_gamma_1d(gamma_external), 2.5e-8,
                   INVERSE_KM_TO_INVERSE_M, ROUNDOFF_REL_TOL, false);
  check_conversion(context, metrics, "3D inverse length", gamma_external,
                   "km^-1", inferred_gamma_3d(gamma_external), 2.5e-8,
                   INVERSE_KM_TO_INVERSE_M, ROUNDOFF_REL_TOL, false);
}

void check_roundtrip_support(swcme_test::Context& context, Metrics& metrics) {
  std::cout << "CFG02 round-trip conversions\n";
  // The 1-D writers expose selected output conversions. Parse their generated
  // columns so the inverse path is production code, not a test-side division.
  double recovered_velocity = 0.0;
  const bool velocity_parsed =
      read_1d_shock_velocity(400.0, recovered_velocity);
  context.record_result(velocity_parsed);
  if (velocity_parsed) {
    check_roundtrip(context, metrics, "1D velocity round trip", 400.0,
                    recovered_velocity);
  } else {
    std::cout << "1D velocity writer parse       FAIL\n";
  }
  print_skip(context, "3D velocity round trip",
             "inverse production conversion not exposed");

  double recovered_au = 0.0, recovered_rs = 0.0;
  const bool distance_parsed = read_1d_profile_units(
      0.5 * swcme1d::AU, recovered_au, recovered_rs);
  context.record_result(distance_parsed);
  if (distance_parsed) {
    check_roundtrip(context, metrics, "1D AU round trip", 0.5, recovered_au);
    check_roundtrip(context, metrics, "1D Rs round trip",
                    0.5 * swcme1d::AU / swcme1d::Rs, recovered_rs);
  } else {
    std::cout << "1D distance writer parse       FAIL\n";
  }
  print_skip(context, "3D AU/Rs round trip",
             "inverse production conversion not exposed");

  // The remaining quantities are emitted in SI or have no production inverse.
  print_skip(context, "magnetic-field round trip",
             "inverse production conversion not exposed");
  print_skip(context, "density round trip",
             "inverse production conversion not exposed");
  print_skip(context, "time round trip",
             "forward/inverse production hour conversion not exposed");
  print_skip(context, "angle round trip",
             "inverse production conversion not exposed");
  print_skip(context, "inverse-length round trip",
             "inverse production conversion not exposed");
}

void check_1d_3d_consistency(swcme_test::Context& context, Metrics& metrics) {
  std::cout << "CFG02 1D/3D consistency\n";
  const double velocities[] = {0.0, 1.0, 321.0, 400.0, 1500.0, 3000.0};
  for (double value : velocities) {
    check_model_consistency(context, metrics, "velocity " + std::to_string(value),
                            velocity_1d(value), velocity_3d(value));
  }

  const double magnetic_fields[] = {0.1, 1.0, 5.0, 7.0, 100.0, 1000.0};
  for (double value : magnetic_fields) {
    check_model_consistency(context, metrics,
                            "magnetic field " + std::to_string(value),
                            ambient_1d(5.0, value).magnetic_field_T,
                            ambient_3d(5.0, value).magnetic_field_T);
  }

  const double densities[] = {0.01, 1.0, 5.0, 8.0, 100.0, 1.0e4};
  for (double value : densities) {
    check_model_consistency(context, metrics,
                            "density " + std::to_string(value),
                            ambient_1d(value, 5.0).density_m3,
                            ambient_3d(value, 5.0).density_m3);
  }

  const double distances[] = {0.1, 0.5, 1.0, 2.0, 5.0};
  for (double value : distances) {
    check_model_consistency(context, metrics,
                            "AU distance " + std::to_string(value),
                            value * swcme1d::AU, value * swcme3d::AU);
  }
  check_model_consistency(context, metrics, "inverse length",
                          inferred_gamma_1d(2.5e-5),
                          inferred_gamma_3d(2.5e-5));

  print_skip(context, "time consistency",
             "neither model exposes an hours conversion path");
  print_skip(context, "angle consistency",
             "only 3-D exposes a degree-derived default");
}

void check_equivalent_configuration(swcme_test::Context& context,
                                    Metrics& metrics) {
  std::cout << "CFG02 equivalent-configuration test\n";

  // This is the requested normal heliophysics configuration. Time and sampled
  // radius enter production in SI because Params has no hour or AU fields for
  // those quantities; all other external values use the documented Params path.
  swcme1d::Params parameters_1d;
  parameters_1d.V_sw_kms = 400.0;
  parameters_1d.B1AU_nT = 5.0;
  parameters_1d.n1AU_cm3 = 5.0;
  const swcme1d::Model model_1d(parameters_1d);
  const swcme1d::StepState state_1d = model_1d.prepare_step(86400.0);
  const AmbientFields fields_1d = ambient_1d(5.0, 5.0);

  check_conversion(context, metrics, "1D config velocity", 400.0, "km/s",
                   state_1d.V_up_ms, 400000.0, KM_S_TO_M_S, 0.0, true);
  check_conversion(context, metrics, "1D config magnetic field", 5.0, "nT",
                   fields_1d.magnetic_field_T, 5.0e-9, NT_TO_T,
                   ROUNDOFF_REL_TOL, false);
  check_conversion(context, metrics, "1D config density", 5.0, "cm^-3",
                   fields_1d.density_m3, 5.0e6, CM3_TO_M3,
                   ROUNDOFF_REL_TOL, false);
  check_conversion(context, metrics, "1D config radius", 1.0, "AU",
                   swcme1d::AU, REF_AU_M, REF_AU_M, 0.0, true);
  check_conversion(context, metrics, "1D config time", 86400.0, "s",
                   state_1d.time_s, 86400.0, 1.0, 0.0, true);
  print_skip(context, "1D direct-SI configuration",
             "no unified direct-SI Params path exists for field and density");

  swcme3d::Params parameters_3d;
  parameters_3d.V_sw_kms = 400.0;
  parameters_3d.B1AU_nT = 5.0;
  parameters_3d.n1AU_cm3 = 5.0;
  parameters_3d.half_width_rad = 40.0 * REF_PI / 180.0;
  const swcme3d::Model model_3d(parameters_3d);
  const swcme3d::StepState state_3d = model_3d.prepare_step(86400.0);
  const AmbientFields fields_3d = ambient_3d(5.0, 5.0);

  check_conversion(context, metrics, "3D config velocity", 400.0, "km/s",
                   state_3d.V_sw_ms, 400000.0, KM_S_TO_M_S, 0.0, true);
  check_conversion(context, metrics, "3D config magnetic field", 5.0, "nT",
                   fields_3d.magnetic_field_T, 5.0e-9, NT_TO_T,
                   ROUNDOFF_REL_TOL, false);
  check_conversion(context, metrics, "3D config density", 5.0, "cm^-3",
                   fields_3d.density_m3, 5.0e6, CM3_TO_M3,
                   ROUNDOFF_REL_TOL, false);
  check_conversion(context, metrics, "3D config radius", 1.0, "AU",
                   swcme3d::AU, REF_AU_M, REF_AU_M, 0.0, true);
  check_conversion(context, metrics, "3D config angle", 40.0, "deg",
                   parameters_3d.half_width_rad,
                   0.69813170079773183077, REF_PI / 180.0,
                   ROUNDOFF_REL_TOL, false);
  print_skip(context, "3D prepared-state time",
             "StepState does not retain the prepare_step time argument");
  print_skip(context, "3D direct-SI configuration",
             "no unified direct-SI Params path exists for field and density");
}

void check_alfven_speed(swcme_test::Context& context, Metrics& metrics) {
  std::cout << "CFG02 Alfven-speed dimensional smoke test\n";
  const AmbientFields fields_1d = ambient_1d(5.0, 5.0);
  const AmbientFields fields_3d = ambient_3d(5.0, 5.0);
  const double reference_density_m3 = 5.0e6;
  const double reference_field_T = 5.0e-9;
  const double reference_rho = reference_density_m3 * REF_PROTON_MASS_KG;
  const double reference_va =
      reference_field_T / std::sqrt(REF_MU0_N_A2 * reference_rho);

  const AmbientFields fields[] = {fields_1d, fields_3d};
  const char* labels[] = {"1D", "3D"};
  for (int i = 0; i < 2; ++i) {
    const double external_rho = fields[i].density_m3 *
                                swcme::constants::PROTON_MASS_KG;
    const double external_va = swcme::physics::alfven_speed_m_s(
        fields[i].magnetic_field_T, external_rho);
    const double direct_si_va = swcme::physics::alfven_speed_m_s(
        reference_field_T,
        reference_density_m3 * swcme::constants::PROTON_MASS_KG);
    const double reference_error = relative_error(external_va, reference_va);
    const double path_error = relative_error(external_va, direct_si_va);
    metrics.alfven_speed_relative_error = std::max(
        metrics.alfven_speed_relative_error,
        std::max(reference_error, path_error));

    std::cout << labels[i]
              << " B[T]=" << std::scientific << std::setprecision(12)
              << fields[i].magnetic_field_T
              << " n[m^-3]=" << fields[i].density_m3
              << " rho[kg/m^3]=" << external_rho
              << " mu0=" << swcme::constants::VACUUM_PERMEABILITY_N_A2
              << " SWCME_VA[m/s]=" << external_va
              << " direct_SI_VA[m/s]=" << direct_si_va
              << " reference_VA[m/s]=" << reference_va
              << " rel_err=" << reference_error
              << " tol=" << ROUNDOFF_REL_TOL
              << ' ' << ((reference_error <= ROUNDOFF_REL_TOL &&
                           path_error <= ROUNDOFF_REL_TOL) ? "PASS" : "FAIL")
              << '\n';
    context.record_result(reference_error <= ROUNDOFF_REL_TOL &&
                          path_error <= ROUNDOFF_REL_TOL);
  }
}

void print_metrics(const swcme_test::Context& context, const Metrics& metrics) {
  std::cout << "CFG02 structured metrics\n"
            << "number_of_checks=" << context.checks() << '\n'
            << "number_passed=" << context.passes() << '\n'
            << "number_failed=" << context.failures() << '\n'
            << "number_skipped=" << context.skips() << '\n'
            << "maximum_relative_conversion_error=" << std::scientific
            << std::setprecision(12)
            << metrics.maximum_relative_conversion_error << '\n';
  if (metrics.any_roundtrip_executed) {
    std::cout << "maximum_roundtrip_error="
              << metrics.maximum_roundtrip_error << '\n';
  } else {
    std::cout << "maximum_roundtrip_error=N/A (no production inverse path)\n";
  }
  std::cout << "maximum_1d_3d_difference="
            << metrics.maximum_1d_3d_difference << '\n'
            << "alfven_speed_relative_error="
            << metrics.alfven_speed_relative_error << '\n';
}

}  // namespace

void test_cfg02(swcme_test::Context& context) {
  Metrics metrics;
  std::cout << "CFG02 diagnostic: independent 1-D and 3-D production modules "
               "still duplicate km/s x1e3, cm^-3 x1e6, nT x1e-9, and "
               "km^-1 /1e3 conversion factors.\n";
  check_forward_conversions(context, metrics);
  check_roundtrip_support(context, metrics);
  check_1d_3d_consistency(context, metrics);
  check_equivalent_configuration(context, metrics);
  check_alfven_speed(context, metrics);
  print_metrics(context, metrics);
}
