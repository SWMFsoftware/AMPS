#include "test_framework.hpp"

#include <swcme1d.hpp>
#include <swcme3d.hpp>
#include <swcme_constants.hpp>
#include <swcme_config.hpp>
#include <swcme_units.hpp>

#include <cmath>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>

namespace {

// Independent references for the immutable constant baseline.  CFG01 is now
// primarily a configuration-rejection test, but retaining this compact check
// guards the shared constants used by both public model interfaces.
constexpr double REF_AU_M = 149597870700.0;
constexpr double REF_SOLAR_RADIUS_M = 695700000.0;
constexpr double REF_PROTON_MASS_KG = 1.67262192595e-27;
constexpr double REF_MU0_N_A2 = 1.25663706127e-6;
constexpr double REF_BOLTZMANN_J_K = 1.380649e-23;

void check_equal(swcme_test::Context& context, const std::string& label,
                 double actual, double expected, double tolerance = 0.0) {
  const double error = std::abs(actual - expected);
  const bool pass = std::isfinite(actual) && error <= tolerance;
  std::cout << std::left << std::setw(40) << label << std::right
            << " actual=" << std::scientific << std::setprecision(12) << actual
            << " expected=" << expected << " abs_err=" << error
            << " tol=" << tolerance << ' ' << (pass ? "PASS" : "FAIL") << '\n';
  context.record_result(pass);
}

bool contains_field(const swcme::config::ValidationResult& result,
                    const std::string& field) {
  for (const auto& issue : result.issues) {
    if (issue.field == field) return true;
  }
  return false;
}

void print_validation(const swcme::config::ValidationResult& result) {
  if (result.ok()) {
    std::cout << "  validation=OK\n";
    return;
  }
  for (const auto& issue : result.issues) {
    std::cout << "  field=" << issue.field
              << " code=" << swcme::config::code_name(issue.code)
              << " requirement=" << issue.requirement << '\n';
  }
}

// Check both the side-effect-free validation API and the prepare_step() guard.
// A model is considered correctly protected only if invalid input is visible
// before execution *and* cannot accidentally enter the physics path.
void expect_invalid_1d(swcme_test::Context& context, const std::string& label,
                       const swcme1d::Params& parameters,
                       const std::string& expected_field) {
  const swcme1d::Model model(parameters);
  const auto result = model.validate();
  const bool validation_pass = !result.ok() && contains_field(result, expected_field);
  std::cout << "CFG01 1D " << label << '\n';
  print_validation(result);
  context.record_result(validation_pass);

  bool threw = false;
  try {
    (void)model.prepare_step(0.0);
  } catch (const std::invalid_argument&) {
    threw = true;
  }
  std::cout << "  prepare_step rejection=" << (threw ? "PASS" : "FAIL") << '\n';
  context.record_result(threw);
}

void expect_invalid_3d(swcme_test::Context& context, const std::string& label,
                       const swcme3d::Params& parameters,
                       const std::string& expected_field) {
  const swcme3d::Model model(parameters);
  const auto result = model.validate();
  const bool validation_pass = !result.ok() && contains_field(result, expected_field);
  std::cout << "CFG01 3D " << label << '\n';
  print_validation(result);
  context.record_result(validation_pass);

  bool threw = false;
  try {
    (void)model.prepare_step(0.0);
  } catch (const std::invalid_argument&) {
    threw = true;
  }
  std::cout << "  prepare_step rejection=" << (threw ? "PASS" : "FAIL") << '\n';
  context.record_result(threw);
}

}  // namespace

void test_cfg01(swcme_test::Context& context) {
  std::cout << "CFG01 shared constants baseline\n";
  check_equal(context, "astronomical unit [m]", swcme::constants::AU_M,
              REF_AU_M);
  check_equal(context, "solar radius [m]", swcme::constants::SOLAR_RADIUS_M,
              REF_SOLAR_RADIUS_M);
  check_equal(context, "proton mass [kg]", swcme::constants::PROTON_MASS_KG,
              REF_PROTON_MASS_KG);
  check_equal(context, "vacuum permeability [N/A^2]",
              swcme::constants::VACUUM_PERMEABILITY_N_A2, REF_MU0_N_A2);
  check_equal(context, "Boltzmann constant [J/K]",
              swcme::constants::BOLTZMANN_J_K, REF_BOLTZMANN_J_K);

  std::cout << "CFG01 valid-default acceptance\n";
  const swcme1d::Model valid_1d;
  const swcme3d::Model valid_3d(swcme3d::Params{});
  context.record_result(valid_1d.validate().ok());
  context.record_result(valid_3d.validate().ok());

  // ---------------------------- Common invalid cases -----------------------
  // Mutate one field at a time from an otherwise valid default configuration.
  // This makes every failure diagnostic attributable to one documented rule.
  {
    swcme1d::Params p; p.V_sw_kms = 0.0;
    expect_invalid_1d(context, "zero solar-wind speed", p, "V_sw_kms");
  }
  {
    swcme3d::Params p; p.V_sw_kms = 0.0;
    expect_invalid_3d(context, "zero solar-wind speed", p, "V_sw_kms");
  }
  {
    swcme1d::Params p; p.n1AU_cm3 = -1.0;
    expect_invalid_1d(context, "negative reference density", p, "n1AU_cm3");
  }
  {
    swcme3d::Params p; p.n1AU_cm3 = -1.0;
    expect_invalid_3d(context, "negative reference density", p, "n1AU_cm3");
  }
  {
    swcme1d::Params p; p.Gamma_kmInv = -1.0e-7;
    expect_invalid_1d(context, "negative DBM Gamma", p, "Gamma_kmInv");
  }
  {
    swcme3d::Params p; p.Gamma_kmInv = -1.0e-7;
    expect_invalid_3d(context, "negative DBM Gamma", p, "Gamma_kmInv");
  }
  {
    swcme1d::Params p; p.sin_theta = 1.1;
    expect_invalid_1d(context, "sin(theta) above one", p, "sin_theta");
  }
  {
    swcme3d::Params p; p.sin_theta = -0.1;
    expect_invalid_3d(context, "negative reference sin(theta)", p, "sin_theta");
  }
  {
    swcme1d::Params p; p.T_K = std::numeric_limits<double>::quiet_NaN();
    expect_invalid_1d(context, "NaN temperature", p, "T_K");
  }
  {
    swcme3d::Params p; p.edge_smooth_le_AU_at1AU = -0.01;
    expect_invalid_3d(context, "negative smoothing width", p,
                      "edge_smooth_le_AU_at1AU");
  }

  // Malformed DATA_DRIVEN tables must be rejected at the public configuration
  // boundary, before the SI kinematics component is called.
  {
    swcme1d::Params p;
    p.kinematics_mode = swcme::kinematics::Mode::DataDriven;
    p.data_time_s = {0.0, 10.0, 10.0};
    p.data_radius_Rs = {20.0, 21.0, 22.0};
    expect_invalid_1d(context, "duplicate data-driven time", p, "data_time_s");
  }
  {
    swcme3d::Params p;
    p.kinematics_mode = swcme::kinematics::Mode::DataDriven;
    p.data_time_s = {0.0, 10.0, 20.0};
    p.data_radius_Rs = {20.0, 22.0, 21.0};
    expect_invalid_3d(context, "decreasing data-driven radius", p,
                      "data_radius_Rs");
  }

  // ----------------------------- 3-D-only cases ----------------------------
  {
    swcme3d::Params p; p.cme_dir[0]=0.0; p.cme_dir[1]=0.0; p.cme_dir[2]=0.0;
    expect_invalid_3d(context, "zero CME direction", p, "cme_dir");
  }
  {
    swcme3d::Params p;
    p.solar_rotation_axis[0]=0.0; p.solar_rotation_axis[1]=0.0;
    p.solar_rotation_axis[2]=0.0;
    expect_invalid_3d(context, "zero solar rotation axis", p,
                      "solar_rotation_axis");
  }
  {
    swcme3d::Params p; p.shape=swcme3d::ShockShape::SSE; p.half_width_rad=0.0;
    expect_invalid_3d(context, "zero SSE half width", p, "half_width_rad");
  }
  {
    swcme3d::Params p; p.shape=swcme3d::ShockShape::SSE;
    p.half_width_rad=0.5*swcme::constants::PI + 1.0e-3;
    expect_invalid_3d(context, "SSE half width above pi/2", p,
                      "half_width_rad");
  }
  {
    swcme3d::Params p; p.shape=swcme3d::ShockShape::Ellipsoid; p.axis_ratio_y=0.0;
    expect_invalid_3d(context, "zero ellipsoid axis ratio", p, "axis_ratio_y");
  }
  {
    swcme3d::Params p; p.solar_rotation_rate_rad_s=-1.0;
    expect_invalid_3d(context, "negative solar rotation rate", p,
                      "solar_rotation_rate_rad_s");
  }

  // The validator reports all independent defects in one pass.  This is useful
  // for event configuration files where correcting one field at a time would
  // otherwise require repeated expensive setup attempts.
  {
    swcme3d::Params p;
    p.V_sw_kms = 0.0;
    p.n1AU_cm3 = -1.0;
    p.cme_dir[0]=p.cme_dir[1]=p.cme_dir[2]=0.0;
    const auto result = swcme3d::Model(p).validate();
    const bool pass = result.issues.size() >= 3 &&
                      contains_field(result,"V_sw_kms") &&
                      contains_field(result,"n1AU_cm3") &&
                      contains_field(result,"cme_dir");
    std::cout << "CFG01 multiple-field diagnostics count=" << result.issues.size()
              << ' ' << (pass ? "PASS" : "FAIL") << '\n';
    context.record_result(pass);
  }
}
