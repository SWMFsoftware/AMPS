#include "test_framework.hpp"

#include <swcme1d.hpp>

// Test-only inclusion exposes the file-local constants actually used by the
// 3-D implementation without changing the production API or source.
#if defined(__GNUC__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wformat-security"
#endif
#include "../../swcme3d.cpp"
#if defined(__GNUC__)
#pragma GCC diagnostic pop
#endif

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>

namespace {

// Independent references:
//   AU: IAU 2012 Resolution B2 (exact).
//   Solar radius: IAU 2015 Resolution B3 nominal conversion (exact).
//   mp, mu0, e, kB: 2022 CODATA recommended values (e and kB are exact).
//   Solar rotation: Carrington sidereal period of 25.38 days.
constexpr double REF_AU_M = 149597870700.0;
constexpr double REF_SOLAR_RADIUS_M = 695700000.0;
constexpr double REF_PROTON_MASS_KG = 1.67262192595e-27;
constexpr double REF_PROTON_MASS_STD_UNCERTAINTY_KG = 0.00000000052e-27;
constexpr double REF_ELEMENTARY_CHARGE_C = 1.602176634e-19;
constexpr double REF_MU0_N_A2 = 1.25663706127e-6;
constexpr double REF_MU0_STD_UNCERTAINTY_N_A2 = 0.00000000020e-6;
constexpr double REF_BOLTZMANN_J_K = 1.380649e-23;
constexpr double REF_PI = 3.141592653589793238462643383279502884;
constexpr double REF_SOLAR_ROTATION_RAD_S =
    REF_PI * 2.0 / (25.38 * 86400.0);

// Physical-constant tolerances account for two CODATA standard uncertainties
// plus half a unit in the final decimal place stored by SWCME. Exact constants
// use zero tolerance. Calculated conversions allow only floating-point roundoff.
constexpr double PROTON_MASS_ABS_TOL =
    2.0 * REF_PROTON_MASS_STD_UNCERTAINTY_KG + 0.5e-38;
constexpr double MU0_ABS_TOL =
    2.0 * REF_MU0_STD_UNCERTAINTY_N_A2 + 0.5e-21;
constexpr double SOLAR_ROTATION_ABS_TOL = 0.5e-11;
constexpr double ROUNDOFF_REL_TOL =
    64.0 * std::numeric_limits<double>::epsilon();

enum class ToleranceKind { Absolute, Relative };

void check_value(swcme_test::Context& context, const std::string& quantity,
                 double actual, double reference, double tolerance,
                 ToleranceKind tolerance_kind) {
  const double absolute_error = std::abs(actual - reference);
  const double relative_error = reference == 0.0
                                    ? absolute_error
                                    : absolute_error / std::abs(reference);
  const bool pass = std::isfinite(actual) &&
                    (tolerance_kind == ToleranceKind::Absolute
                         ? absolute_error <= tolerance
                         : relative_error <= tolerance);

  std::cout << std::left << std::setw(35) << quantity << std::right
            << " swcme=" << std::scientific << std::setprecision(12) << actual
            << " ref=" << reference
            << " rel_err=" << relative_error
            << " abs_err=" << absolute_error
            << " tol(" << (tolerance_kind == ToleranceKind::Absolute ? "abs" : "rel")
            << ")=" << tolerance
            << ' ' << (pass ? "PASS" : "FAIL") << '\n';

  context.record_result(pass);
}

double inferred_gamma(double upstream_speed_ms, double initial_speed_ms,
                      double evaluated_speed_ms, double time_s) {
  const double initial_excess = initial_speed_ms - upstream_speed_ms;
  const double evaluated_excess = evaluated_speed_ms - upstream_speed_ms;
  return (initial_excess / evaluated_excess - 1.0) /
         (initial_excess * time_s);
}

}  // namespace

void test_cfg01(swcme_test::Context& context) {
  std::cout << "CFG01 quantity checks\n";

  check_value(context, "1D astronomical unit [m]", swcme1d::AU, REF_AU_M, 0.0,
              ToleranceKind::Absolute);
  check_value(context, "3D astronomical unit [m]", swcme3d::AU, REF_AU_M, 0.0,
              ToleranceKind::Absolute);
  check_value(context, "1D solar radius [m]", swcme1d::Rs,
              REF_SOLAR_RADIUS_M, 0.0, ToleranceKind::Absolute);
  check_value(context, "3D solar radius [m]", swcme3d::Rs,
              REF_SOLAR_RADIUS_M, 0.0, ToleranceKind::Absolute);
  check_value(context, "1D pi", swcme1d::PI, REF_PI, 0.0,
              ToleranceKind::Absolute);
  check_value(context, "3D pi", swcme3d::PI, REF_PI, 0.0,
              ToleranceKind::Absolute);
  check_value(context, "1D proton mass [kg]", swcme1d::MP,
              REF_PROTON_MASS_KG, PROTON_MASS_ABS_TOL,
              ToleranceKind::Absolute);
  check_value(context, "3D proton mass [kg]", ::MP, REF_PROTON_MASS_KG,
              PROTON_MASS_ABS_TOL, ToleranceKind::Absolute);
  check_value(context, "1D vacuum permeability [N/A^2]", swcme1d::MU0,
              REF_MU0_N_A2, MU0_ABS_TOL, ToleranceKind::Absolute);
  check_value(context, "3D vacuum permeability [N/A^2]", ::MU0,
              REF_MU0_N_A2, MU0_ABS_TOL, ToleranceKind::Absolute);
  check_value(context, "1D Boltzmann constant [J/K]", swcme1d::KB,
              REF_BOLTZMANN_J_K, 0.0, ToleranceKind::Absolute);
  check_value(context, "3D Boltzmann constant [J/K]", ::KB,
              REF_BOLTZMANN_J_K, 0.0, ToleranceKind::Absolute);
  check_value(context, "1D solar rotation [rad/s]", swcme1d::OMEGA_SUN,
              REF_SOLAR_ROTATION_RAD_S, SOLAR_ROTATION_ABS_TOL,
              ToleranceKind::Absolute);
  check_value(context, "3D solar rotation [rad/s]", ::OMEGA_SUN,
              REF_SOLAR_ROTATION_RAD_S, SOLAR_ROTATION_ABS_TOL,
              ToleranceKind::Absolute);

  std::cout << "CFG01 1-D versus 3-D consistency checks\n";
  check_value(context, "1D/3D astronomical unit", swcme1d::AU, swcme3d::AU,
              0.0, ToleranceKind::Absolute);
  check_value(context, "1D/3D solar radius", swcme1d::Rs, swcme3d::Rs,
              0.0, ToleranceKind::Absolute);
  check_value(context, "1D/3D pi", swcme1d::PI, swcme3d::PI,
              0.0, ToleranceKind::Absolute);
  check_value(context, "1D/3D proton mass", swcme1d::MP, ::MP,
              0.0, ToleranceKind::Absolute);
  check_value(context, "1D/3D vacuum permeability", swcme1d::MU0, ::MU0,
              0.0, ToleranceKind::Absolute);
  check_value(context, "1D/3D Boltzmann constant", swcme1d::KB, ::KB,
              0.0, ToleranceKind::Absolute);
  check_value(context, "1D/3D solar rotation", swcme1d::OMEGA_SUN,
              ::OMEGA_SUN, 0.0, ToleranceKind::Absolute);

  std::cout << std::left << std::setw(35) << "elementary charge [C]" << std::right
            << " swcme=N/A ref=" << std::scientific << std::setprecision(12)
            << REF_ELEMENTARY_CHARGE_C
            << " rel_err=N/A abs_err=N/A tol=N/A NOT PRESENT (not used)\n";

  swcme1d::Params parameters_1d;
  parameters_1d.r0_Rs = 2.0;
  parameters_1d.V_sw_kms = 321.0;
  parameters_1d.V0_sh_kms = 654.0;
  parameters_1d.Gamma_kmInv = 2.5e-5;
  parameters_1d.n1AU_cm3 = 8.0;
  parameters_1d.B1AU_nT = 7.0;
  const swcme1d::Model model_1d(parameters_1d);
  const swcme1d::StepState state_1d = model_1d.prepare_step(0.0);
  const swcme1d::StepState evolved_1d = model_1d.prepare_step(10.0);

  const double radius_1d[] = {swcme1d::AU};
  double density_1d[] = {0.0};
  double speed_1d[] = {0.0};
  double br_1d[] = {0.0};
  double bphi_1d[] = {0.0};
  double bmag_1d[] = {0.0};
  model_1d.evaluate_radii_with_B_div(state_1d, radius_1d, density_1d,
                                     speed_1d, br_1d, bphi_1d, bmag_1d,
                                     nullptr, 1);

  check_value(context, "1D solar-radius to metre", state_1d.r0_m,
              1391400000.0, 0.0, ToleranceKind::Absolute);
  check_value(context, "1D km/s to m/s", state_1d.V_up_ms, 321000.0, 0.0,
              ToleranceKind::Absolute);
  check_value(context, "1D cm^-3 to m^-3", density_1d[0], 8000000.0,
              ROUNDOFF_REL_TOL, ToleranceKind::Relative);
  check_value(context, "1D nT to T", bmag_1d[0], 7.0e-9,
              ROUNDOFF_REL_TOL, ToleranceKind::Relative);
  check_value(context, "1D km^-1 to m^-1",
              inferred_gamma(state_1d.V_up_ms, 654000.0,
                             evolved_1d.V_sh_ms, 10.0),
              2.5e-8, ROUNDOFF_REL_TOL, ToleranceKind::Relative);

  swcme3d::Params parameters_3d;
  parameters_3d.r0_Rs = 2.0;
  parameters_3d.V_sw_kms = 321.0;
  parameters_3d.V0_sh_kms = 654.0;
  parameters_3d.Gamma_kmInv = 2.5e-5;
  parameters_3d.n1AU_cm3 = 8.0;
  parameters_3d.B1AU_nT = 7.0;
  const swcme3d::Model model_3d(parameters_3d);
  const swcme3d::StepState state_3d = model_3d.prepare_step(0.0);
  const swcme3d::StepState evolved_3d = model_3d.prepare_step(10.0);

  const double x_3d[] = {swcme3d::AU};
  const double y_3d[] = {0.0};
  const double z_3d[] = {0.0};
  double density_3d[] = {0.0};
  double vx_3d[] = {0.0};
  double vy_3d[] = {0.0};
  double vz_3d[] = {0.0};
  double bx_3d[] = {0.0};
  double by_3d[] = {0.0};
  double bz_3d[] = {0.0};
  model_3d.evaluate_cartesian_with_B(state_3d, x_3d, y_3d, z_3d,
                                     density_3d, vx_3d, vy_3d, vz_3d,
                                     bx_3d, by_3d, bz_3d, 1);
  const double bmag_3d =
      std::sqrt(bx_3d[0] * bx_3d[0] + by_3d[0] * by_3d[0] +
                bz_3d[0] * bz_3d[0]);

  check_value(context, "3D solar-radius to metre", state_3d.r_sh_m,
              1391400000.0, 0.0, ToleranceKind::Absolute);
  check_value(context, "3D km/s to m/s", state_3d.V_sw_ms, 321000.0, 0.0,
              ToleranceKind::Absolute);
  check_value(context, "3D cm^-3 to m^-3", density_3d[0], 8000000.0,
              ROUNDOFF_REL_TOL, ToleranceKind::Relative);
  check_value(context, "3D nT to T", bmag_3d, 7.0e-9,
              ROUNDOFF_REL_TOL, ToleranceKind::Relative);
  check_value(context, "3D km^-1 to m^-1",
              inferred_gamma(state_3d.V_sw_ms, 654000.0,
                             evolved_3d.V_sh_ms, 10.0),
              2.5e-8, ROUNDOFF_REL_TOL, ToleranceKind::Relative);
  check_value(context, "3D degree to radian (40 deg)",
              swcme3d::Params{}.half_width_rad,
              0.69813170079773183077, ROUNDOFF_REL_TOL,
              ToleranceKind::Relative);
}
