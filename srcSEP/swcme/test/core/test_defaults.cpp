#include "test_framework.hpp"

#include "../../swcme_defaults.hpp"
#include "../../swcme1d.hpp"
#include "../../swcme3d.hpp"

#include <cmath>
#include <iostream>
#include <string>
#include <vector>

namespace {

void expect_same(swcme_test::Context& context, double a, double b,
                 const char* label) {
  const double scale = std::max(1.0, std::max(std::abs(a), std::abs(b)));
  context.expect_near(a, b, 8.0e-15 * scale, label);
}

bool contains_key(const std::string& text, const std::string& key) {
  return text.find(key + "=") != std::string::npos;
}

}  // namespace

void test_def01(swcme_test::Context& context) {
  std::cout << "DEF01 canonical 1-D/3-D default equivalence\n";

  const swcme1d::Params p1;
  const swcme3d::Params p3;

  // Every dimensionality-independent physical default must be identical.  This
  // is intentionally explicit rather than comparing one derived state: it
  // catches future drift even when two different inputs happen to produce a
  // similar value at one verification point.
  expect_same(context, p1.V_sw_kms, p3.V_sw_kms, "V_sw default");
  expect_same(context, p1.n1AU_cm3, p3.n1AU_cm3, "n1AU default");
  expect_same(context, p1.B1AU_nT, p3.B1AU_nT, "B1AU default");
  expect_same(context, p1.T_K, p3.T_K, "T default");
  expect_same(context, p1.gamma_ad, p3.gamma_ad, "gamma default");
  expect_same(context, p1.sin_theta, p3.sin_theta,
              "Parker normalization latitude default");
  expect_same(context, p1.r0_Rs, p3.r0_Rs, "DBM r0 default");
  expect_same(context, p1.V0_sh_kms, p3.V0_sh_kms, "V0 default");
  expect_same(context, p1.Gamma_kmInv, p3.Gamma_kmInv, "DBM Gamma default");
  expect_same(context, p1.relative_source_weight_per_area,
              p3.relative_source_weight_per_area, "source weight default");
  expect_same(context, p1.sheath_thick_AU_at1AU,
              p3.sheath_thick_AU_at1AU, "sheath thickness default");
  expect_same(context, p1.ejecta_thick_AU_at1AU,
              p3.ejecta_thick_AU_at1AU, "ejecta thickness default");
  expect_same(context, p1.edge_smooth_shock_AU_at1AU,
              p3.edge_smooth_shock_AU_at1AU, "shock width default");
  expect_same(context, p1.edge_smooth_le_AU_at1AU,
              p3.edge_smooth_le_AU_at1AU, "LE width default");
  expect_same(context, p1.edge_smooth_te_AU_at1AU,
              p3.edge_smooth_te_AU_at1AU, "TE width default");
  expect_same(context, p1.sheath_comp_floor, p3.sheath_comp_floor,
              "deprecated sheath floor neutral default");
  expect_same(context, p1.sheath_ramp_power, p3.sheath_ramp_power,
              "sheath ramp default");
  expect_same(context, p1.V_sheath_LE_factor, p3.V_sheath_LE_factor,
              "sheath LE speed default");
  expect_same(context, p1.f_ME, p3.f_ME, "ejecta density default");
  expect_same(context, p1.V_ME_factor, p3.V_ME_factor,
              "ejecta speed default");

  context.expect_true(p1.kinematics_mode == p3.kinematics_mode,
                      "1-D/3-D kinematic mode defaults differ");
  context.expect_true(p1.data_extrapolation == p3.data_extrapolation,
                      "1-D/3-D extrapolation defaults differ");
  context.expect_true(p1.region_mode == p3.region_mode,
                      "1-D/3-D region defaults differ");
  context.expect_true(p1.shock_acceleration_mode == p3.shock_acceleration_mode,
                      "1-D/3-D acceleration defaults differ");

  context.expect_true(p1.n1AU_cm3 == swcme::defaults::N1AU_CM3,
                      "public defaults must come from swcme_defaults.hpp");
  context.expect_true(p1.r0_Rs == swcme::defaults::DBM_R0_RS,
                      "DBM start radius is not canonical");

  // Also pass both default bundles through the production preparation path.
  // This catches a wrapper that overrides one canonical default while copying
  // Params into the common SI core.
  const swcme1d::Model m1(p1);
  const swcme3d::Model m3(p3);
  const auto s1 = m1.prepare_step(3600.0);
  const auto s3 = m3.prepare_step(3600.0);
  expect_same(context, s1.common.solar_wind.V_sw_m_s,
              s3.common.solar_wind.V_sw_m_s, "prepared V_sw identity");
  expect_same(context, s1.common.solar_wind.Br1AU_T,
              s3.common.solar_wind.Br1AU_T, "prepared Br1AU identity");
  expect_same(context, s1.common.apex.radius_m, s3.common.apex.radius_m,
              "prepared apex radius identity");
  expect_same(context, s1.common.apex.speed_m_s, s3.common.apex.speed_m_s,
              "prepared apex speed identity");
}

void test_def02(swcme_test::Context& context) {
  std::cout << "DEF02 science-scope and geometry conventions\n";

  const swcme1d::Params p1;
  const swcme3d::Params p3;
  const swcme1d::Model m1(p1);
  const swcme3d::Model m3(p3);

  context.expect_true(p1.region_mode == swcme::regions::Mode::ShockOnly,
                      "default region mode must be SHOCK_ONLY");
  context.expect_true(p1.shock_acceleration_mode ==
                          swcme::acceleration::Mode::Source,
                      "default acceleration mode must be SOURCE");
  context.expect_true(m1.model_scope() ==
                          swcme::defaults::ModelScope::ControlledSEPPreShock,
                      "1-D default scope is not CONTROLLED_SEP_PRE_SHOCK");
  context.expect_true(m3.model_scope() ==
                          swcme::defaults::ModelScope::ControlledSEPPreShock,
                      "3-D default scope is not CONTROLLED_SEP_PRE_SHOCK");

  // The science geometry is finite SSE. Sphere remains a verification geometry
  // and must therefore be selected explicitly by tests or callers that need it.
  context.expect_true(p3.shape == swcme3d::ShockShape::SSE,
                      "default 3-D science geometry must be finite SSE");
  context.expect_near(p3.half_width_rad, swcme::defaults::SSE_HALF_WIDTH_RAD,
                      0.0, "SSE half-width default");

  context.expect_true(swcme::defaults::PARKER_RADIAL_POLARITY == +1,
                      "baseline Parker radial polarity must be outward");
  context.expect_true(std::string(swcme::defaults::FRAME_NAME) ==
                          "HCI_like_inertial",
                      "baseline frame convention changed unexpectedly");
  context.expect_true(std::string(swcme::defaults::PARKER_NORMALIZATION_CONVENTION) ==
                          "TOTAL_B_AT_1AU_REFERENCE_LATITUDE",
                      "Parker normalization convention changed unexpectedly");
  context.expect_true(p3.sin_theta == 1.0,
                      "3-D normalization latitude must remain equatorial");
}

void test_def03(swcme_test::Context& context) {
  std::cout << "DEF03 observer-local model-scope gating\n";

  swcme1d::Params p1;
  p1.kinematics_mode = swcme::kinematics::Mode::Ballistic;
  p1.V0_sh_kms = 1000.0;
  const swcme1d::Model m1(p1);

  const double observer_r = swcme::constants::AU_M;
  const auto s1_pre = m1.prepare_step(0.0);
  const auto g1_pre = m1.observer_scope_status(s1_pre, observer_r);
  context.expect_true(g1_pre.valid_input && g1_pre.within_declared_scope &&
                          !g1_pre.shock_has_reached_observer,
                      "1-D pre-shock observer must be inside controlled scope");

  const auto s1_post = m1.prepare_step(2.0e5);
  const auto g1_post = m1.observer_scope_status(s1_post, observer_r);
  context.expect_true(g1_post.valid_input && !g1_post.within_declared_scope &&
                          g1_post.shock_has_reached_observer,
                      "1-D post-arrival Parker background must be flagged out of scope");

  swcme3d::Params p3;
  p3.kinematics_mode = swcme::kinematics::Mode::Ballistic;
  p3.V0_sh_kms = 1000.0;
  p3.shape = swcme3d::ShockShape::SSE;
  const swcme3d::Model m3(p3);
  const double obs_apex[3] = {observer_r, 0.0, 0.0};
  const double obs_outside_cap[3] = {0.0, observer_r, 0.0};

  const auto s3_pre = m3.prepare_step(0.0);
  const auto g3_pre = m3.observer_scope_status(s3_pre, obs_apex);
  context.expect_true(g3_pre.valid_input && g3_pre.within_declared_scope &&
                          g3_pre.shock_surface_on_observer_ray &&
                          !g3_pre.shock_has_reached_observer,
                      "3-D apex observer pre-shock scope classification failed");

  const auto s3_post = m3.prepare_step(2.0e5);
  const auto g3_post = m3.observer_scope_status(s3_post, obs_apex);
  context.expect_true(g3_post.valid_input && !g3_post.within_declared_scope &&
                          g3_post.shock_has_reached_observer,
                      "3-D post-arrival observer must be flagged out of controlled scope");

  // A finite front outside the observer direction is not allowed to trigger a
  // false "shock arrived" scope failure merely because its apex passed 1 AU.
  const auto g3_outside = m3.observer_scope_status(s3_post, obs_outside_cap);
  context.expect_true(g3_outside.valid_input && g3_outside.within_declared_scope &&
                          !g3_outside.shock_surface_on_observer_ray &&
                          !g3_outside.shock_has_reached_observer,
                      "finite-SSE out-of-cap observer scope classification failed");
}

void test_def04(swcme_test::Context& context) {
  std::cout << "DEF04 resolved-configuration manifest completeness\n";

  swcme1d::Params p1;
  p1.n1AU_cm3 = 7.25;  // explicit event override must be visible in metadata
  p1.data_time_s = {0.0, 10.0};
  p1.data_radius_Rs = {20.0, 21.0};
  const std::string m1 = swcme1d::resolved_configuration_manifest(p1);

  swcme3d::Params p3;
  p3.n1AU_cm3 = 7.25;
  p3.cme_dir[0] = 0.5; p3.cme_dir[1] = std::sqrt(0.75); p3.cme_dir[2] = 0.0;
  p3.data_time_s = {0.0, 10.0};
  p3.data_radius_Rs = {20.0, 21.0};
  const std::string m3 = swcme3d::resolved_configuration_manifest(p3);

  const std::vector<std::string> common_keys = {
      "swcme_config_version", "model_scope", "frame", "parker_normalization",
      "parker_radial_polarity", "solar_rotation_rate_rad_s", "V_sw_kms", "n1AU_cm3", "B1AU_nT", "T_K",
      "gamma_ad", "thermodynamic_closure", "alpha_to_proton_ratio",
      "electron_T_K", "alpha_T_K", "kinematics_mode", "r0_Rs",
      "V0_sh_kms", "Gamma_kmInv", "parker_source_radius_Rs",
      "data_extrapolation", "data_time_s.count", "data_radius_Rs.count",
      "region_mode", "shock_acceleration_mode", "relative_source_weight_per_area",
      "sheath_thick_AU_at1AU", "ejecta_thick_AU_at1AU",
      "edge_smooth_shock_AU_at1AU", "edge_smooth_le_AU_at1AU",
      "edge_smooth_te_AU_at1AU", "sheath_comp_floor", "sheath_ramp_power",
      "V_sheath_LE_factor", "f_ME", "V_ME_factor"};
  for (const std::string& key : common_keys) {
    context.expect_true(contains_key(m1, key), "1-D manifest missing " + key);
    context.expect_true(contains_key(m3, key), "3-D manifest missing " + key);
  }

  const std::vector<std::string> geometry_keys = {
      "shape", "axis_ratio_y", "axis_ratio_z", "half_width_rad",
      "flank_slowdown_m", "cme_dir[0]", "cme_dir[1]", "cme_dir[2]",
      "solar_rotation_axis[0]", "solar_rotation_axis[1]",
      "solar_rotation_axis[2]", "solar_rotation_rate_rad_s"};
  for (const std::string& key : geometry_keys) {
    context.expect_true(contains_key(m3, key), "3-D manifest missing " + key);
  }

  context.expect_true(m1.find("n1AU_cm3=7.25000000000000000e+00") !=
                          std::string::npos,
                      "1-D event override missing from resolved manifest");
  context.expect_true(m3.find("n1AU_cm3=7.25000000000000000e+00") !=
                          std::string::npos,
                      "3-D event override missing from resolved manifest");
  context.expect_true(m3.find("model_scope=CONTROLLED_SEP_PRE_SHOCK") !=
                          std::string::npos,
                      "manifest must record the derived model scope");

  // Serialization must be deterministic so later campaign manifests can be
  // diffed and hashed without order-dependent metadata noise.
  context.expect_true(m1 == swcme1d::resolved_configuration_manifest(p1),
                      "1-D resolved manifest is not deterministic");
  context.expect_true(m3 == swcme3d::resolved_configuration_manifest(p3),
                      "3-D resolved manifest is not deterministic");
}
