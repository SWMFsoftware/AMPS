#include "test_framework.hpp"

#include <swcme1d.hpp>
#include <swcme3d.hpp>
#include <swcme_sep_interface.hpp>

#include <algorithm>
#include <cmath>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>

namespace {

// DEN01 deliberately repeats the published Leblanc analytical coefficients in
// test code. This is an independent physical reference, not a copy used in
// place of production: every value under test is obtained through a production
// Model, StepState, and field evaluator.
constexpr double REF_LEBLANC_A_CM3 = 3.3e5;
constexpr double REF_LEBLANC_B_CM3 = 4.1e6;
constexpr double REF_LEBLANC_C_CM3 = 8.0e7;

// AU and the nominal solar radius are the documented, adopted SWCME constants
// validated by CFG01/CFG02. Their ratio is evaluated here instead of storing a
// rounded AU-in-solar-radii conversion that could silently drift.
constexpr double REF_AU_M = 149597870700.0;
constexpr double REF_SOLAR_RADIUS_M = 695700000.0;
constexpr double CM3_TO_M3 = 1.0e6;

// The normalization criterion is the validation specification's mandatory
// physics tolerance. It is intentionally much stricter than a typical model
// comparison because normalization at the defining radius is algebraic.
constexpr double NORMALIZATION_REL_TOL = 1.0e-12;

// Cross-model, coefficient, and scale-factor comparisons contain only short
// double-precision arithmetic chains. Sixty-four machine epsilons allows for
// operation ordering without masking a physical unit or normalization error.
constexpr double ROUNDOFF_REL_TOL =
    64.0 * std::numeric_limits<double>::epsilon();

struct ReferenceLeblanc {
  double term_a_cm3;
  double term_b_cm3;
  double term_c_cm3;
  double total_cm3;
};

struct Metrics {
  int number_of_fixtures = 0;
  int number_passed = 0;
  int number_failed = 0;
  double max_normalization_relative_error = 0.0;
  double max_1d_3d_relative_difference = 0.0;
  bool cross_instance_contamination = false;
  bool coefficients_mutated = false;
};

// This small scalar implementation is intentionally independent of the
// production density evaluator. The radius is dimensionless in nominal solar
// radii and the three returned terms and total are in cm^-3.
ReferenceLeblanc reference_leblanc_unscaled_cm3(double radius_rs) {
  const double inv_r2 = 1.0 / (radius_rs * radius_rs);
  const double inv_r4 = inv_r2 * inv_r2;
  const double inv_r6 = inv_r4 * inv_r2;
  const double term_a = REF_LEBLANC_A_CM3 * inv_r2;
  const double term_b = REF_LEBLANC_B_CM3 * inv_r4;
  const double term_c = REF_LEBLANC_C_CM3 * inv_r6;
  return {term_a, term_b, term_c, term_a + term_b + term_c};
}

double relative_error(double actual, double expected) {
  const double absolute = std::abs(actual - expected);
  return expected == 0.0 ? absolute : absolute / std::abs(expected);
}

// Record a generic roundoff comparison while retaining a diagnostic for the
// exact source of a mismatch. This is separate from the primary 1e-12 density
// normalization residual used by each physical fixture.
bool check_roundoff(swcme_test::Context& context, const std::string& label,
                    double actual, double expected) {
  const double error = relative_error(actual, expected);
  const bool pass = std::isfinite(actual) && error <= ROUNDOFF_REL_TOL;
  std::cout << "  " << std::left << std::setw(39) << label << std::right
            << " actual=" << std::scientific << std::setprecision(12)
            << actual << " expected=" << expected
            << " rel_err=" << error << " tol=" << ROUNDOFF_REL_TOL
            << ' ' << (pass ? "PASS" : "FAIL") << '\n';
  context.record_result(pass);
  return pass;
}

// Evaluate the public 1-D field path. At t=0 the shock is close to the Sun, so
// 1 AU is in the upstream region and the returned density is the production
// Leblanc density without a CME/sheath contribution.
double evaluate_1d_density_m3(const swcme1d::Model& model,
                              const swcme1d::StepState& state,
                              double radius_m) {
  double density_m3 = 0.0;
  double velocity_m_s = 0.0;
  model.evaluate_radii_fast(state, &radius_m, &density_m3, &velocity_m_s, 1);
  return density_m3;
}

// Evaluate the public 3-D field path at an explicitly supplied Cartesian
// point. The default spherical shock is well inside 1 AU at t=0, isolating the
// radial upstream Leblanc profile from downstream CME structure.
double evaluate_3d_density_m3(const swcme3d::Model& model,
                              const swcme3d::StepState& state,
                              double x_m, double y_m, double z_m) {
  double density_m3 = 0.0;
  double vx_m_s = 0.0;
  double vy_m_s = 0.0;
  double vz_m_s = 0.0;
  model.evaluate_cartesian_fast(state, &x_m, &y_m, &z_m, &density_m3,
                                &vx_m_s, &vy_m_s, &vz_m_s, 1);
  return density_m3;
}

// Print the complete independent reference trace once for each model/n_ref
// fixture. The individual 3-D direction results follow beneath this header.
void print_reference_trace(const char* model_name, double n_ref_cm3,
                           double radius_rs, const ReferenceLeblanc& reference,
                           double expected_scale) {
  std::cout << "\nDEN01 fixture: model=" << model_name
            << " n_ref=" << std::fixed << std::setprecision(6) << n_ref_cm3
            << " cm^-3\n"
            << std::scientific << std::setprecision(12)
            << "  r_ref [m]                       = " << REF_AU_M << '\n'
            << "  r_ref [Rs]                      = " << radius_rs << '\n'
            << "  Leblanc nominal terms:\n"
            << "    A/r^2 [cm^-3]                 = " << reference.term_a_cm3 << '\n'
            << "    B/r^4 [cm^-3]                 = " << reference.term_b_cm3 << '\n'
            << "    C/r^6 [cm^-3]                 = " << reference.term_c_cm3 << '\n'
            << "  unscaled total [cm^-3]          = " << reference.total_cm3 << '\n'
            << "  expected scale factor           = " << expected_scale << '\n';
}

// The defining DEN01 check uses the scale-independent residual
// abs(n_model/n_ref - 1). It is strict '< 1e-12', exactly as specified, and
// therefore remains meaningful over the full two-decade fixture range.
bool check_normalization(swcme_test::Context& context, Metrics& metrics,
                         const std::string& label, double production_density_m3,
                         double n_ref_cm3) {
  const double production_cm3 = production_density_m3 / CM3_TO_M3;
  const double absolute_error = std::abs(production_cm3 - n_ref_cm3);
  const double residual = std::abs(production_cm3 / n_ref_cm3 - 1.0);
  const bool pass = std::isfinite(production_cm3) &&
                    residual < NORMALIZATION_REL_TOL;

  ++metrics.number_of_fixtures;
  if (pass) {
    ++metrics.number_passed;
  } else {
    ++metrics.number_failed;
  }
  metrics.max_normalization_relative_error =
      std::max(metrics.max_normalization_relative_error, residual);

  std::cout << "  " << label << '\n'
            << "    production density [cm^-3]    = " << std::scientific
            << std::setprecision(12) << production_cm3 << '\n'
            << "    expected density [cm^-3]      = " << n_ref_cm3 << '\n'
            << "    absolute error [cm^-3]        = " << absolute_error << '\n'
            << "    relative normalization error  = " << residual << '\n'
            << "    tolerance (strictly less than)= " << NORMALIZATION_REL_TOL
            << '\n'
            << "    " << (pass ? "PASS" : "FAIL") << '\n';
  context.record_result(pass);
  return pass;
}

double inferred_scale_from_c2(double coefficient_c2, double radius_m) {
  return coefficient_c2 /
         (REF_LEBLANC_A_CM3 * CM3_TO_M3 * radius_m * radius_m);
}

// Recover each nominal coefficient from the exposed per-step SI cache using
// the independent expected scale. This detects a changed coefficient, power,
// density conversion, or radius conversion without calling the production
// Leblanc evaluator as the reference.
bool check_cached_coefficients(swcme_test::Context& context,
                               const std::string& label, double c2, double c4,
                               double c6, double expected_scale) {
  const double rs2 = REF_SOLAR_RADIUS_M * REF_SOLAR_RADIUS_M;
  const double rs4 = rs2 * rs2;
  const double rs6 = rs4 * rs2;
  const double recovered_a = c2 / (expected_scale * CM3_TO_M3 * rs2);
  const double recovered_b = c4 / (expected_scale * CM3_TO_M3 * rs4);
  const double recovered_c = c6 / (expected_scale * CM3_TO_M3 * rs6);
  bool pass = true;
  pass &= check_roundoff(context, label + " recovered A", recovered_a,
                         REF_LEBLANC_A_CM3);
  pass &= check_roundoff(context, label + " recovered B", recovered_b,
                         REF_LEBLANC_B_CM3);
  pass &= check_roundoff(context, label + " recovered C", recovered_c,
                         REF_LEBLANC_C_CM3);
  return pass;
}

swcme1d::Params make_1d_params(double n_ref_cm3) {
  swcme1d::Params parameters;
  parameters.n1AU_cm3 = n_ref_cm3;
  return parameters;
}

swcme3d::Params make_3d_params(double n_ref_cm3) {
  swcme3d::Params parameters;
  parameters.n1AU_cm3 = n_ref_cm3;
  return parameters;
}

void run_reference_fixtures(swcme_test::Context& context, Metrics& metrics) {
  const double radius_rs = REF_AU_M / REF_SOLAR_RADIUS_M;
  const ReferenceLeblanc reference =
      reference_leblanc_unscaled_cm3(radius_rs);
  const double densities_cm3[] = {1.0, 5.0, 8.0, 20.0, 100.0};

  // Unit vectors include all Cartesian axes and a non-axis-aligned direction.
  // Multiplying each by the adopted AU checks that the 3-D result depends only
  // on heliocentric radius rather than direction.
  const double inv_sqrt_14 = 1.0 / std::sqrt(14.0);
  const double directions[][3] = {
      {1.0, 0.0, 0.0},
      {0.0, 1.0, 0.0},
      {0.0, 0.0, 1.0},
      {inv_sqrt_14, 2.0 * inv_sqrt_14, 3.0 * inv_sqrt_14},
  };
  const char* direction_names[] = {"+X", "+Y", "+Z", "normalized (1,2,3)"};

  std::cout << "DEN01 normalization fixtures\n";
  for (double n_ref_cm3 : densities_cm3) {
    const double expected_scale = n_ref_cm3 / reference.total_cm3;

    const swcme1d::Model model_1d(make_1d_params(n_ref_cm3));
    const swcme1d::StepState state_1d = model_1d.prepare_step(0.0);
    const double density_1d_m3 =
        evaluate_1d_density_m3(model_1d, state_1d, swcme1d::AU);

    print_reference_trace("1D", n_ref_cm3, radius_rs, reference,
                          expected_scale);
    check_normalization(context, metrics, "r=production AU",
                        density_1d_m3, n_ref_cm3);
    const double scale_1d =
        inferred_scale_from_c2(state_1d.C2, swcme1d::Rs);
    check_roundoff(context, "1D exposed/inferred scale", scale_1d,
                    expected_scale);
    check_cached_coefficients(context, "1D", state_1d.C2, state_1d.C4,
                              state_1d.C6, expected_scale);

    const swcme3d::Model model_3d(make_3d_params(n_ref_cm3));
    const swcme3d::StepState state_3d = model_3d.prepare_step(0.0);
    print_reference_trace("3D", n_ref_cm3, radius_rs, reference,
                          expected_scale);

    double density_3d_x_m3 = 0.0;
    for (std::size_t i = 0; i < sizeof(directions) / sizeof(directions[0]);
         ++i) {
      const double x_m = swcme3d::AU * directions[i][0];
      const double y_m = swcme3d::AU * directions[i][1];
      const double z_m = swcme3d::AU * directions[i][2];
      const double magnitude_m = std::sqrt(x_m * x_m + y_m * y_m + z_m * z_m);
      const double density_m3 = evaluate_3d_density_m3(
          model_3d, state_3d, x_m, y_m, z_m);
      if (i == 0) density_3d_x_m3 = density_m3;

      std::cout << "  direction=" << direction_names[i]
                << " x,y,z[m]=(" << std::scientific << std::setprecision(12)
                << x_m << ',' << y_m << ',' << z_m << ") |x|[m]="
                << magnitude_m << '\n';
      check_normalization(context, metrics,
                          std::string("3D direction ") + direction_names[i],
                          density_m3, n_ref_cm3);
    }

    const double scale_3d =
        inferred_scale_from_c2(state_3d.C2, swcme3d::Rs);
    check_roundoff(context, "3D exposed/inferred scale", scale_3d,
                    expected_scale);
    check_cached_coefficients(context, "3D", state_3d.C2, state_3d.C4,
                              state_3d.C6, expected_scale);

    const double model_difference =
        relative_error(density_1d_m3, density_3d_x_m3);
    metrics.max_1d_3d_relative_difference =
        std::max(metrics.max_1d_3d_relative_difference, model_difference);
    check_roundoff(context, "1D/3D density at 1 AU", density_1d_m3,
                    density_3d_x_m3);
    check_roundoff(context, "1D/3D scale factor", scale_1d, scale_3d);

    // Both APIs accept SI meters, so the second input reconstructs one AU via
    // the documented AU/Rs ratio and solar-radius path. This checks coordinate
    // equivalence without duplicating CFG02's broader conversion suite.
    const double radius_via_rs_m = radius_rs * swcme1d::Rs;
    const double density_1d_rs_m3 =
        evaluate_1d_density_m3(model_1d, state_1d, radius_via_rs_m);
    check_roundoff(context, "1D AU versus AU/Rs coordinate path",
                    density_1d_rs_m3, density_1d_m3);
    const double density_3d_rs_m3 = evaluate_3d_density_m3(
        model_3d, state_3d, radius_via_rs_m, 0.0, 0.0);
    check_roundoff(context, "3D AU versus AU/Rs coordinate path",
                    density_3d_rs_m3, density_3d_x_m3);
  }
}

void run_cross_instance_test_1d(swcme_test::Context& context,
                                Metrics& metrics) {
  std::cout << "\nDEN01 1D cross-instance contamination sequence\n";
  const swcme1d::Model model_a(make_1d_params(5.0));
  const swcme1d::StepState state_a = model_a.prepare_step(0.0);
  const double a_c2_before = state_a.C2;
  const double a_c4_before = state_a.C4;
  const double a_c6_before = state_a.C6;
  const double a_density_before = evaluate_1d_density_m3(
      model_a, state_a, swcme1d::AU);
  check_normalization(
      context, metrics, "Model A before B", evaluate_1d_density_m3(
          model_a, state_a, swcme1d::AU), 5.0);

  const swcme1d::Model model_b(make_1d_params(20.0));
  const swcme1d::StepState state_b = model_b.prepare_step(0.0);
  const double b_density_before = evaluate_1d_density_m3(
      model_b, state_b, swcme1d::AU);
  check_normalization(context, metrics, "Model B before C",
                      b_density_before, 20.0);
  const double a_density_after_b = evaluate_1d_density_m3(
      model_a, state_a, swcme1d::AU);
  check_normalization(context, metrics, "Model A after B",
                      a_density_after_b, 5.0);

  const swcme1d::Model model_c(make_1d_params(1.0));
  const swcme1d::StepState state_c = model_c.prepare_step(0.0);
  check_normalization(context, metrics, "Model C",
                      evaluate_1d_density_m3(model_c, state_c, swcme1d::AU),
                      1.0);
  const double b_density_after_c = evaluate_1d_density_m3(
      model_b, state_b, swcme1d::AU);
  const double a_density_after_c = evaluate_1d_density_m3(
      model_a, state_a, swcme1d::AU);
  check_normalization(context, metrics, "Model B after C",
                      b_density_after_c, 20.0);
  check_normalization(context, metrics, "Model A after C",
                      a_density_after_c, 5.0);

  // Exact cache identity is appropriate here: these are snapshots in a
  // value-type state and no arithmetic should be applied to them by creating
  // another model. A change would demonstrate shared mutable state directly.
  const bool state_unchanged = state_a.C2 == a_c2_before &&
                               state_a.C4 == a_c4_before &&
                               state_a.C6 == a_c6_before;
  context.record_result(state_unchanged);
  std::cout << "  Model A cached coefficients C2/C4/C6 unchanged "
            << (state_unchanged ? "PASS" : "FAIL") << '\n';
  const bool behavior_unchanged =
      relative_error(a_density_after_b, a_density_before) <= ROUNDOFF_REL_TOL &&
      relative_error(a_density_after_c, a_density_before) <= ROUNDOFF_REL_TOL &&
      relative_error(b_density_after_c, b_density_before) <= ROUNDOFF_REL_TOL;
  context.record_result(behavior_unchanged);
  std::cout << "  Model A/B before-versus-after behavior unchanged "
            << (behavior_unchanged ? "PASS" : "FAIL") << '\n';
  metrics.cross_instance_contamination |= !behavior_unchanged;
  metrics.coefficients_mutated |= !state_unchanged;
}

void run_cross_instance_test_3d(swcme_test::Context& context,
                                Metrics& metrics) {
  std::cout << "\nDEN01 3D cross-instance contamination sequence\n";
  const swcme3d::Model model_a(make_3d_params(5.0));
  const swcme3d::StepState state_a = model_a.prepare_step(0.0);
  const double a_c2_before = state_a.C2;
  const double a_c4_before = state_a.C4;
  const double a_c6_before = state_a.C6;
  const double a_density_before = evaluate_3d_density_m3(
      model_a, state_a, swcme3d::AU, 0.0, 0.0);
  check_normalization(context, metrics, "Model A before B",
                      a_density_before, 5.0);

  const swcme3d::Model model_b(make_3d_params(20.0));
  const swcme3d::StepState state_b = model_b.prepare_step(0.0);
  const double b_density_before = evaluate_3d_density_m3(
      model_b, state_b, swcme3d::AU, 0.0, 0.0);
  check_normalization(context, metrics, "Model B before C",
                      b_density_before, 20.0);
  const double a_density_after_b = evaluate_3d_density_m3(
      model_a, state_a, swcme3d::AU, 0.0, 0.0);
  check_normalization(context, metrics, "Model A after B",
                      a_density_after_b, 5.0);

  const swcme3d::Model model_c(make_3d_params(1.0));
  const swcme3d::StepState state_c = model_c.prepare_step(0.0);
  check_normalization(context, metrics, "Model C",
                      evaluate_3d_density_m3(
                          model_c, state_c, swcme3d::AU, 0.0, 0.0), 1.0);
  const double b_density_after_c = evaluate_3d_density_m3(
      model_b, state_b, swcme3d::AU, 0.0, 0.0);
  const double a_density_after_c = evaluate_3d_density_m3(
      model_a, state_a, swcme3d::AU, 0.0, 0.0);
  check_normalization(context, metrics, "Model B after C",
                      b_density_after_c, 20.0);
  check_normalization(context, metrics, "Model A after C",
                      a_density_after_c, 5.0);

  const bool state_unchanged = state_a.C2 == a_c2_before &&
                               state_a.C4 == a_c4_before &&
                               state_a.C6 == a_c6_before;
  context.record_result(state_unchanged);
  std::cout << "  Model A cached coefficients C2/C4/C6 unchanged "
            << (state_unchanged ? "PASS" : "FAIL") << '\n';
  const bool behavior_unchanged =
      relative_error(a_density_after_b, a_density_before) <= ROUNDOFF_REL_TOL &&
      relative_error(a_density_after_c, a_density_before) <= ROUNDOFF_REL_TOL &&
      relative_error(b_density_after_c, b_density_before) <= ROUNDOFF_REL_TOL;
  context.record_result(behavior_unchanged);
  std::cout << "  Model A/B before-versus-after behavior unchanged "
            << (behavior_unchanged ? "PASS" : "FAIL") << '\n';
  metrics.cross_instance_contamination |= !behavior_unchanged;
  metrics.coefficients_mutated |= !state_unchanged;
}

void run_radial_sanity(swcme_test::Context& context) {
  std::cout << "\nDEN01 optional radial sanity (DEN02 owns radial accuracy)\n";
  const swcme1d::Model model_1d(make_1d_params(5.0));
  const swcme1d::StepState state_1d = model_1d.prepare_step(0.0);
  const swcme3d::Model model_3d(make_3d_params(5.0));
  const swcme3d::StepState state_3d = model_3d.prepare_step(0.0);
  const double radii_au[] = {0.5, 2.0};
  for (double radius_au : radii_au) {
    const double radius_m = radius_au * REF_AU_M;
    const double density_1d =
        evaluate_1d_density_m3(model_1d, state_1d, radius_m);
    const double density_3d =
        evaluate_3d_density_m3(model_3d, state_3d, radius_m, 0.0, 0.0);
    const bool pass_1d = std::isfinite(density_1d) && density_1d > 0.0;
    const bool pass_3d = std::isfinite(density_3d) && density_3d > 0.0;
    std::cout << "  r=" << radius_au << " AU 1D=" << std::scientific
              << density_1d << " m^-3 " << (pass_1d ? "PASS" : "FAIL")
              << " 3D=" << density_3d << " m^-3 "
              << (pass_3d ? "PASS" : "FAIL") << '\n';
    context.record_result(pass_1d);
    context.record_result(pass_3d);
  }
}

void print_metrics(const swcme_test::Context& context, const Metrics& metrics) {
  std::cout << "\nDEN01 structured metrics\n"
            << "test_id=DEN01\n"
            << "classification=COMMON\n"
            << "number_of_fixtures=" << metrics.number_of_fixtures << '\n'
            << "number_passed=" << metrics.number_passed << '\n'
            << "number_failed=" << metrics.number_failed << '\n'
            << "number_skipped=" << context.skips() << '\n'
            << "max_normalization_relative_error=" << std::scientific
            << std::setprecision(12)
            << metrics.max_normalization_relative_error << '\n'
            << "max_1d_3d_relative_difference="
            << metrics.max_1d_3d_relative_difference << '\n'
            << "cross_instance_contamination="
            << (metrics.cross_instance_contamination ? "true" : "false") << '\n'
            << "coefficients_mutated="
            << (metrics.coefficients_mutated ? "true" : "false") << '\n'
            << "total_subchecks=" << context.checks() << '\n'
            << "total_subchecks_passed=" << context.passes() << '\n'
            << "total_subchecks_failed=" << context.failures() << '\n';
}

}  // namespace

void test_den01(swcme_test::Context& context) {
  Metrics metrics;
  std::cout << "DEN01 diagnostic: 1-D and 3-D independently implement and "
               "cache the same Leblanc coefficients and 1-AU normalization; "
               "neither uses mutable global normalization state.\n";
  run_reference_fixtures(context, metrics);
  run_cross_instance_test_1d(context, metrics);
  run_cross_instance_test_3d(context, metrics);
  run_radial_sanity(context);

  // Neither Params type exposes a configurable normalization radius. This is
  // capability reporting rather than a failed physics check: both production
  // models intentionally fix the reference distance at one adopted AU.
  std::cout << "\nDEN01 generalized normalization radius: SKIP — production "
               "SWCME currently fixes density normalization at 1 AU\n";
  context.record_skip();
  print_metrics(context, metrics);
}

void test_den04(swcme_test::Context& context) {
  std::cout << "DEN04 independent pressure and sound-speed closure\n";

  // These literals deliberately do not reference the production constants.
  // The test therefore detects accidental changes to a mass or to the closure
  // algebra instead of merely exercising both sides of one shared helper.
  constexpr long double kB = 1.380649e-23L;
  constexpr long double mp = 1.67262192595e-27L;
  constexpr long double ma = 6.6446573450e-27L;
  const double gammas[] = {1.2, 1.4, 1.5, 5.0 / 3.0};
  double max_relative_error = 0.0;
  std::size_t fixture_count = 0;

  for (double gamma : gammas) {
    for (int fixture = 0; fixture < 4; ++fixture) {
      swcme::solarwind::ConfigSI cfg;
      cfg.V_sw_m_s = 4.0e5;
      cfg.n1AU_m3 = 5.0e6;
      cfg.B1AU_T = 5.0e-9;
      cfg.T_K = fixture == 3 ? 1.0e-6 : (fixture == 2 ? 8.0e4 : 1.2e5);
      cfg.gamma_ad = gamma;
      cfg.reference_sin_theta = 1.0;
      cfg.electron_T_K = fixture == 2 ? 2.7e5 : cfg.T_K;
      cfg.alpha_T_K = fixture == 2 ? 6.1e5 : cfg.T_K;
      cfg.alpha_to_proton_ratio = fixture == 1 ? 0.0 : 0.05;
      cfg.thermodynamic_closure = fixture == 0
          ? swcme::solarwind::ThermodynamicClosure::ProtonOnly
          : swcme::solarwind::ThermodynamicClosure::MultiSpecies;
      const auto prepared = swcme::solarwind::prepare(cfg);
      const long double ne = 8.25e6L;
      const auto actual = swcme::solarwind::thermodynamic_state(
          prepared, static_cast<double>(ne));

      long double np = ne;
      long double na = 0.0L;
      long double rho = mp * ne;
      long double pressure = ne * kB * cfg.T_K;
      if (cfg.thermodynamic_closure ==
          swcme::solarwind::ThermodynamicClosure::MultiSpecies) {
        np = ne / (1.0L + 2.0L * cfg.alpha_to_proton_ratio);
        na = cfg.alpha_to_proton_ratio * np;
        rho = mp * np + ma * na;
        pressure = kB * (np * cfg.T_K + ne * cfg.electron_T_K +
                          na * cfg.alpha_T_K);
      }
      const long double sound = std::sqrt(
          static_cast<long double>(gamma) * pressure / rho);
      const double errors[] = {
          relative_error(actual.proton_density_m3, static_cast<double>(np)),
          relative_error(actual.alpha_density_m3, static_cast<double>(na)),
          relative_error(actual.mass_density_kg_m3, static_cast<double>(rho)),
          relative_error(actual.pressure_Pa, static_cast<double>(pressure)),
          relative_error(actual.sound_speed_m_s, static_cast<double>(sound))};
      for (double error : errors) {
        max_relative_error = std::max(max_relative_error, error);
        context.expect_true(error < 1.0e-13,
                            "thermodynamic closure residual exceeds 1e-13");
      }
      ++fixture_count;
    }
  }

  // Verify that both dimensional shock APIs consume mass and pressure from
  // the selected closure, rather than computing a hidden proton-only value.
  swcme1d::Params p1;
  p1.thermodynamic_closure =
      swcme::solarwind::ThermodynamicClosure::MultiSpecies;
  p1.alpha_to_proton_ratio = 0.08;
  p1.electron_T_K = 2.0e5;
  p1.alpha_T_K = 5.0e5;
  const swcme1d::Model m1(p1);
  const auto s1 = m1.prepare_step(0.0);
  const auto t1 = swcme::solarwind::thermodynamic_state(
      s1.common.solar_wind, s1.n_up_shock);
  context.expect_true(relative_error(s1.shock_jump.upstream.rho_kg_m3,
                                     t1.mass_density_kg_m3) < 1.0e-13,
                      "1-D shock did not use selected closure mass density");
  context.expect_true(relative_error(s1.shock_jump.upstream.pressure_Pa,
                                     t1.pressure_Pa) < 1.0e-13,
                      "1-D shock did not use selected closure pressure");

  swcme3d::Params p3;
  p3.thermodynamic_closure = p1.thermodynamic_closure;
  p3.alpha_to_proton_ratio = p1.alpha_to_proton_ratio;
  p3.electron_T_K = p1.electron_T_K;
  p3.alpha_T_K = p1.alpha_T_K;
  const swcme3d::Model m3(p3);
  const auto s3 = m3.prepare_step(0.0);
  const double direction[3] = {1.0, 0.0, 0.0};
  swcme3d::LocalShockState shock3;
  context.expect_true(m3.shock_state_direction(s3, direction, shock3),
                      "3-D apex shock state must exist");
  const double n3 = swcme::solarwind::density_m3(
      s3.common.solar_wind, s3.r_sh_m);
  const auto t3 = swcme::solarwind::thermodynamic_state(
      s3.common.solar_wind, n3);
  context.expect_true(relative_error(shock3.upstream.rho_kg_m3,
                                     t3.mass_density_kg_m3) < 1.0e-13,
                      "3-D shock did not use selected closure mass density");
  context.expect_true(relative_error(shock3.upstream.pressure_Pa,
                                     t3.pressure_Pa) < 1.0e-13,
                      "3-D shock did not use selected closure pressure");

  // Invalid closures and composition parameters must fail at configuration
  // validation, before preparation can create NaN pressure or sound speed.
  swcme1d::Params invalid;
  invalid.thermodynamic_closure =
      static_cast<swcme::solarwind::ThermodynamicClosure>(99);
  context.expect_true(!swcme1d::validate_params(invalid).ok(),
                      "unknown closure enum must be rejected");
  invalid = swcme1d::Params{};
  invalid.alpha_to_proton_ratio = -0.01;
  context.expect_true(!swcme1d::validate_params(invalid).ok(),
                      "negative alpha abundance must be rejected");
  invalid = swcme1d::Params{};
  invalid.electron_T_K = std::numeric_limits<double>::infinity();
  context.expect_true(!swcme1d::validate_params(invalid).ok(),
                      "non-finite electron temperature must be rejected");
  invalid = swcme1d::Params{};
  invalid.alpha_T_K = 0.0;
  context.expect_true(!swcme1d::validate_params(invalid).ok(),
                      "non-positive alpha temperature must be rejected");

  std::cout << "  fixtures=" << fixture_count
            << " max_relative_error=" << std::scientific
            << max_relative_error << "\n";
}

void test_den02(swcme_test::Context& context) {
  std::cout << "DEN02 independent Leblanc radial terms and asymptote\n";
  constexpr long double au = 149597870700.0L;
  constexpr long double rs = 695700000.0L;
  constexpr long double cm3_to_m3 = 1.0e6L;
  constexpr long double A = 3.3e5L;
  constexpr long double B = 4.1e6L;
  constexpr long double C = 8.0e7L;
  constexpr long double configured_n1au = 7.3e6L;

  // Normalize from the independently pinned coefficients at the exact IAU AU;
  // no production coefficient, cache, or density helper enters this oracle.
  const long double q = rs / au;
  const long double q2 = q * q;
  const long double q4 = q2 * q2;
  const long double q6 = q4 * q2;
  const long double scale = configured_n1au /
      (cm3_to_m3 * (A * q2 + B * q4 + C * q6));

  swcme1d::Params p1;
  p1.n1AU_cm3 = static_cast<double>(configured_n1au / cm3_to_m3);
  const swcme1d::Model m1(p1);
  const auto s1 = m1.prepare_step(0.0);
  swcme3d::Params p3;
  p3.n1AU_cm3 = p1.n1AU_cm3;
  const swcme3d::Model m3(p3);
  const auto s3 = m3.prepare_step(0.0);

  double max_total_error = 0.0;
  double max_term_error = 0.0;
  double prior_scaled_excess = std::numeric_limits<double>::infinity();
  std::size_t points = 0;
  for (int i = 0; i <= 40; ++i) {
    // Forty logarithmically spaced points span 0.01--100 AU; the explicitly
    // requested 0.5/1/2/5-AU checkpoints are added below through the same path.
    const long double r_au = std::pow(10.0L, -2.0L + 0.1L * i);
    const long double radius = r_au * au;
    const long double x = rs / radius;
    const long double x2 = x * x;
    const long double ref2 = scale * cm3_to_m3 * A * x2;
    const long double ref4 = scale * cm3_to_m3 * B * x2 * x2;
    const long double ref6 = scale * cm3_to_m3 * C * x2 * x2 * x2;
    const long double reference = ref2 + ref4 + ref6;

    const double radius_d = static_cast<double>(radius);
    const double actual1 = evaluate_1d_density_m3(m1, s1, radius_d);
    const double actual3 = evaluate_3d_density_m3(m3, s3, radius_d, 0.0, 0.0);
    const double cached_terms[] = {
        s1.C2 / (radius_d * radius_d),
        s1.C4 / std::pow(radius_d, 4),
        s1.C6 / std::pow(radius_d, 6)};
    const long double reference_terms[] = {ref2, ref4, ref6};
    for (int term = 0; term < 3; ++term) {
      const double error = relative_error(
          cached_terms[term], static_cast<double>(reference_terms[term]));
      max_term_error = std::max(max_term_error, error);
      context.expect_true(error < 1.0e-12,
                          "Leblanc radial term residual exceeds 1e-12");
    }
    for (double actual : {actual1, actual3}) {
      const double error = relative_error(actual, static_cast<double>(reference));
      max_total_error = std::max(max_total_error, error);
      context.expect_true(error < 1.0e-12,
                          "Leblanc total density residual exceeds 1e-12");
    }

    // r^2*n tends to the r^-2 coefficient.  The independently predicted
    // excess consists exactly of the r^-4 and r^-6 contributions after the
    // multiplication; it must shrink monotonically across this log sweep.
    const long double asymptotic = scale * cm3_to_m3 * A * rs * rs;
    const long double scaled = static_cast<long double>(actual1) * radius * radius;
    const double scaled_excess = static_cast<double>((scaled - asymptotic) /
                                                      asymptotic);
    const double predicted_excess = static_cast<double>((ref4 + ref6) / ref2);
    // At 100 AU the excess is only O(1e-8); forming it subtracts two nearly
    // equal doubles. A 1e-7 relative comparison is still a few femto-units in
    // absolute scaled density and reflects this unavoidable cancellation.
    context.expect_true(relative_error(scaled_excess, predicted_excess) < 1.0e-7,
                        "r^2 n excess is not explained by higher-order terms");
    context.expect_true(scaled_excess >= 0.0 &&
                            scaled_excess <= prior_scaled_excess,
                        "r^2 n must approach its positive asymptote monotonically");
    prior_scaled_excess = scaled_excess;
    ++points;
  }

  const double checkpoints[] = {0.5, 1.0, 2.0, 5.0};
  for (double r_au : checkpoints) {
    const long double radius = r_au * au;
    const long double x = rs / radius;
    const long double x2 = x * x;
    const long double reference = scale * cm3_to_m3 *
        (A * x2 + B * x2 * x2 + C * x2 * x2 * x2);
    const double actual = evaluate_1d_density_m3(
        m1, s1, static_cast<double>(radius));
    context.expect_true(relative_error(actual, static_cast<double>(reference)) <
                            1.0e-12,
                        "explicit radial checkpoint exceeds 1e-12");
  }
  std::cout << "  logarithmic_points=" << points
            << " max_term_error=" << std::scientific << max_term_error
            << " max_total_error=" << max_total_error
            << " far_field_relative_excess=" << prior_scaled_excess << '\n';
}

void test_den05(swcme_test::Context& context) {
  std::cout << "DEN05 exact default-closure compatibility\n";
  const auto same_bits = [](double left, double right) {
    return std::memcmp(&left, &right, sizeof(double)) == 0;
  };

  swcme1d::Params legacy1;
  swcme1d::Params explicit1;
  explicit1.thermodynamic_closure =
      swcme::solarwind::ThermodynamicClosure::ProtonOnly;
  explicit1.alpha_to_proton_ratio = swcme::defaults::ALPHA_TO_PROTON_RATIO;
  explicit1.electron_T_K = swcme::defaults::ELECTRON_T_K;
  explicit1.alpha_T_K = swcme::defaults::ALPHA_T_K;
  const swcme1d::Model legacy_model1(legacy1);
  const swcme1d::Model explicit_model1(explicit1);
  const auto legacy_state1 = legacy_model1.prepare_step(7200.0);
  const auto explicit_state1 = explicit_model1.prepare_step(7200.0);
  const double radii[] = {0.2 * REF_AU_M, REF_AU_M, 3.0 * REF_AU_M};
  double ln[3], lv[3], lbr[3], lbp[3], lbm[3], ldiv[3];
  double en[3], ev[3], ebr[3], ebp[3], ebm[3], ediv[3];
  legacy_model1.evaluate_radii_with_B_div(legacy_state1, radii, ln, lv, lbr,
                                           lbp, lbm, ldiv, 3);
  explicit_model1.evaluate_radii_with_B_div(explicit_state1, radii, en, ev,
                                             ebr, ebp, ebm, ediv, 3);
  for (int i = 0; i < 3; ++i) {
    context.expect_true(same_bits(ln[i], en[i]) && same_bits(lv[i], ev[i]) &&
                            same_bits(lbr[i], ebr[i]) &&
                            same_bits(lbp[i], ebp[i]) &&
                            same_bits(lbm[i], ebm[i]) &&
                            same_bits(ldiv[i], ediv[i]),
                        "default and explicit ProtonOnly 1-D fields differ");
  }
  context.expect_true(same_bits(legacy_state1.shock_jump.upstream.rho_kg_m3,
                                explicit_state1.shock_jump.upstream.rho_kg_m3) &&
                          same_bits(legacy_state1.shock_jump.upstream.pressure_Pa,
                                    explicit_state1.shock_jump.upstream.pressure_Pa) &&
                          same_bits(legacy_state1.rc, explicit_state1.rc),
                      "default and explicit ProtonOnly 1-D shock differ");

  swcme3d::Params legacy3;
  legacy3.shape = swcme3d::ShockShape::Sphere;
  swcme3d::Params explicit3 = legacy3;
  explicit3.thermodynamic_closure =
      swcme::solarwind::ThermodynamicClosure::ProtonOnly;
  explicit3.alpha_to_proton_ratio = swcme::defaults::ALPHA_TO_PROTON_RATIO;
  explicit3.electron_T_K = swcme::defaults::ELECTRON_T_K;
  explicit3.alpha_T_K = swcme::defaults::ALPHA_T_K;
  const swcme3d::Model legacy_model3(legacy3);
  const swcme3d::Model explicit_model3(explicit3);
  const auto legacy_state3 = legacy_model3.prepare_step(7200.0);
  const auto explicit_state3 = explicit_model3.prepare_step(7200.0);
  const double x = 0.8 * REF_AU_M, y = 0.1 * REF_AU_M,
               z = -0.05 * REF_AU_M;
  double lo[8] = {}, eo[8] = {};
  legacy_model3.evaluate_cartesian_with_B_div(
      legacy_state3, &x, &y, &z, &lo[0], &lo[1], &lo[2], &lo[3],
      &lo[4], &lo[5], &lo[6], &lo[7], 1);
  explicit_model3.evaluate_cartesian_with_B_div(
      explicit_state3, &x, &y, &z, &eo[0], &eo[1], &eo[2], &eo[3],
      &eo[4], &eo[5], &eo[6], &eo[7], 1);
  for (int i = 0; i < 8; ++i) {
    context.expect_true(same_bits(lo[i], eo[i]),
                        "default and explicit ProtonOnly 3-D field differ");
  }
  const double u[3] = {1.0, 0.0, 0.0};
  swcme3d::LocalShockState ls, es;
  context.expect_true(legacy_model3.shock_state_direction(legacy_state3, u, ls) &&
                          explicit_model3.shock_state_direction(explicit_state3, u, es),
                      "default-compatibility shock surfaces must exist");
  context.expect_true(same_bits(ls.upstream.rho_kg_m3, es.upstream.rho_kg_m3) &&
                          same_bits(ls.upstream.pressure_Pa,
                                    es.upstream.pressure_Pa) &&
                          same_bits(ls.compression, es.compression),
                      "default and explicit ProtonOnly 3-D shock differ");

  // Adapter background and source serialization cover the AMPS-facing public
  // contract without comparing model-owner tokens, which intentionally differ
  // between independently constructed model instances.
  swcme::sep::Interface1D legacy_adapter(legacy1);
  swcme::sep::Interface1D explicit_adapter(explicit1);
  const auto la_step = legacy_adapter.prepare(7200.0);
  const auto ea_step = explicit_adapter.prepare(7200.0);
  swcme::sep::BackgroundState lb, eb;
  context.expect_true(legacy_adapter.evaluate_background(la_step, REF_AU_M, lb).ok() &&
                          explicit_adapter.evaluate_background(ea_step, REF_AU_M, eb).ok(),
                      "default-compatibility adapter backgrounds succeed");
  context.expect_true(same_bits(lb.density_m3, eb.density_m3) &&
                          same_bits(lb.pressure_Pa, eb.pressure_Pa) &&
                          same_bits(lb.focusing_length_m, eb.focusing_length_m),
                      "default and explicit ProtonOnly adapter background differ");
  swcme::sep::SEPSourceState lsource, esource;
  context.expect_true(legacy_adapter.source_at_shock(la_step, lsource).ok() &&
                          explicit_adapter.source_at_shock(ea_step, esource).ok(),
                      "default-compatibility adapter sources succeed");
  context.expect_true(swcme::sep::serialize_source_csv(lsource) ==
                          swcme::sep::serialize_source_csv(esource),
                      "default and explicit ProtonOnly source records differ");

  const std::string default_manifest =
      swcme1d::resolved_configuration_manifest(legacy1);
  context.expect_true(default_manifest.find("thermodynamic_closure=PROTON_ONLY") !=
                          std::string::npos,
                      "default manifest must record ProtonOnly closure");
  swcme1d::Params multi = legacy1;
  multi.thermodynamic_closure =
      swcme::solarwind::ThermodynamicClosure::MultiSpecies;
  multi.alpha_to_proton_ratio = 0.05;
  const std::string multi_manifest =
      swcme1d::resolved_configuration_manifest(multi);
  context.expect_true(multi_manifest.find("thermodynamic_closure=MULTI_SPECIES") !=
                          std::string::npos &&
                          swcme1d::configuration_digest(multi) !=
                              swcme1d::configuration_digest(legacy1),
                      "MultiSpecies must be explicit and change provenance");
  std::cout << "  1D/3D field, shock, adapter, source, and manifest compatibility PASS\n";
}
