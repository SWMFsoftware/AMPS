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

#include <cmath>
#include <iomanip>
#include <iostream>
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
// use zero tolerance. These tolerances remain part of the CFG01 constants
// baseline even though configuration rejection itself is not yet exposed.
constexpr double PROTON_MASS_ABS_TOL =
    2.0 * REF_PROTON_MASS_STD_UNCERTAINTY_KG + 0.5e-38;
constexpr double MU0_ABS_TOL =
    2.0 * REF_MU0_STD_UNCERTAINTY_N_A2 + 0.5e-21;
constexpr double SOLAR_ROTATION_ABS_TOL = 0.5e-11;

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

}  // namespace

void test_cfg01(swcme_test::Context& context) {
  // CFG01 is intended to exercise explicit configuration rejection. Neither
  // model currently exposes such an API, so the pre-existing constants checks
  // are retained here and the missing rejection path is reported as a SKIP.
  // Unit-conversion checks formerly below this section now live in CFG02.
  std::cout << "CFG01 constants baseline (configuration rejection pending)\n";

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
  context.record_skip();

  std::cout << "configuration rejection API         SKIP — production SWCME "
               "does not expose validation/rejection results\n";
  context.record_skip();
}
