#include "test_framework.hpp"

#include <swcme1d.hpp>
#include <swcme3d.hpp>
#include <swcme_config.hpp>
#include <swcme_regions.hpp>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>

namespace {

// Use one geometry in both dimensions so every expected limit is derived from
// the shared policy rather than duplicated as a decimal in the test.  FULL_ICME
// plus RESOLVED_COMPRESSION makes all three widths physically active.
void equivalent_policy_params(swcme1d::Params& one,swcme3d::Params& three) {
  one.region_mode=swcme::regions::Mode::FullICME;
  one.shock_acceleration_mode=
      swcme::acceleration::Mode::ResolvedCompression;
  one.sheath_thick_AU_at1AU=0.20;
  one.ejecta_thick_AU_at1AU=0.30;

  three.region_mode=one.region_mode;
  three.shock_acceleration_mode=one.shock_acceleration_mode;
  three.sheath_thick_AU_at1AU=one.sheath_thick_AU_at1AU;
  three.ejecta_thick_AU_at1AU=one.ejecta_thick_AU_at1AU;
}

bool has_issue(const swcme::config::ValidationResult& result,
               const std::string& field,swcme::config::Code code) {
  for (const swcme::config::Issue& issue : result.issues) {
    if (issue.field==field && issue.code==code) return true;
  }
  return false;
}

// CFG03 requires rejection during setup, not only a diagnostic from validate().
// These helpers prove prepare_step() is routed through the identical central
// validator and never reaches make_boundaries() with an invalid request.
bool prepare_rejects(const swcme1d::Params& params) {
  try {
    const swcme1d::Model model(params);
    (void)model.prepare_step(3600.0);
  } catch (const std::exception&) {
    return true;
  }
  return false;
}

bool prepare_rejects(const swcme3d::Params& params) {
  try {
    const swcme3d::Model model(params);
    (void)model.prepare_step(3600.0);
  } catch (const std::exception&) {
    return true;
  }
  return false;
}

void expect_rejected_pair(swcme_test::Context& context,
                          const swcme1d::Params& one,
                          const swcme3d::Params& three,
                          const char* field,swcme::config::Code code,
                          const std::string& label) {
  const swcme::config::ValidationResult one_result=
      swcme1d::Model(one).validate();
  const swcme::config::ValidationResult three_result=
      swcme3d::Model(three).validate();
  context.expect_true(has_issue(one_result,field,code),
                      "1-D "+label+" identifies "+field);
  context.expect_true(has_issue(three_result,field,code),
                      "3-D "+label+" identifies "+field);
  context.expect_true(prepare_rejects(one),
                      "1-D "+label+" is rejected before preparation");
  context.expect_true(prepare_rejects(three),
                      "3-D "+label+" is rejected before preparation");
}

void expect_exact_effective_widths(swcme_test::Context& context,
                                   const swcme::regions::Config& config,
                                   const swcme::regions::Boundaries& boundaries,
                                   double radius,const std::string& label) {
  // One direct floating-point multiplication can differ by a final rounding
  // bit when a cached radius/fraction is reloaded under another optimization
  // mode.  A small epsilon-scaled bound accepts only that representation noise;
  // it remains many orders tighter than any historical 90-percent clipping.
  const auto tolerance=[](double actual,double expected) {
    return 8.0*std::numeric_limits<double>::epsilon()*
           std::max({std::abs(actual),std::abs(expected),1.0});
  };
  const double expected_shock=config.shock_smooth_fraction*radius;
  const double expected_leading=config.leading_smooth_fraction*radius;
  const double expected_trailing=config.trailing_smooth_fraction*radius;
  context.expect_near(boundaries.smooth_shock_width_m,
                      expected_shock,
                      tolerance(boundaries.smooth_shock_width_m,expected_shock),
                      label+" preserves shock width exactly");
  context.expect_near(boundaries.smooth_le_width_m,
                      expected_leading,
                      tolerance(boundaries.smooth_le_width_m,expected_leading),
                      label+" preserves leading width exactly");
  context.expect_near(boundaries.smooth_te_width_m,
                      expected_trailing,
                      tolerance(boundaries.smooth_te_width_m,expected_trailing),
                      label+" preserves trailing width exactly");
}

}  // namespace

void test_cfg03(swcme_test::Context& context) {
  std::cout << "CFG03 smoothing width rejection and exact preservation\n";

  swcme1d::Params base_one;
  swcme3d::Params base_three;
  equivalent_policy_params(base_one,base_three);
  const swcme::regions::SmoothingFractionLimits limits=
      swcme::regions::smoothing_fraction_limits(
          base_one.sheath_thick_AU_at1AU,
          base_one.ejecta_thick_AU_at1AU);

  // Representative interior values exercise the ordinary production path,
  // separately from the zero and exact-limit boundaries below.
  base_one.edge_smooth_shock_AU_at1AU=0.04;
  base_one.edge_smooth_le_AU_at1AU=0.06;
  base_one.edge_smooth_te_AU_at1AU=0.08;
  base_three.edge_smooth_shock_AU_at1AU=
      base_one.edge_smooth_shock_AU_at1AU;
  base_three.edge_smooth_le_AU_at1AU=base_one.edge_smooth_le_AU_at1AU;
  base_three.edge_smooth_te_AU_at1AU=base_one.edge_smooth_te_AU_at1AU;
  const swcme1d::Model normal_model_one(base_one);
  const swcme3d::Model normal_model_three(base_three);
  context.expect_true(normal_model_one.validate().ok() &&
                          normal_model_three.validate().ok(),
                      "ordinary in-limit smoothing configurations are valid");
  const swcme1d::StepState normal_state_one=
      normal_model_one.prepare_step(3600.0);
  const swcme3d::StepState normal_state_three=
      normal_model_three.prepare_step(3600.0);
  expect_exact_effective_widths(context,normal_state_one.region_config,
                                normal_state_one.region_boundaries,
                                normal_state_one.r_sh_m,
                                "1-D ordinary state");
  expect_exact_effective_widths(context,normal_state_three.region_config,
                                normal_state_three.apex_regions,
                                normal_state_three.r_sh_m,
                                "3-D ordinary state");

  // Zero is a valid smoothing request when no resolved shock is selected.  It
  // must remain exactly zero in both prepared records rather than being
  // replaced by a default or minimum numerical layer.
  swcme1d::Params zero_one;
  swcme3d::Params zero_three;
  zero_one.edge_smooth_shock_AU_at1AU=0.0;
  zero_one.edge_smooth_le_AU_at1AU=0.0;
  zero_one.edge_smooth_te_AU_at1AU=0.0;
  zero_three.edge_smooth_shock_AU_at1AU=0.0;
  zero_three.edge_smooth_le_AU_at1AU=0.0;
  zero_three.edge_smooth_te_AU_at1AU=0.0;
  const swcme1d::Model zero_model_one(zero_one);
  const swcme3d::Model zero_model_three(zero_three);
  context.expect_true(zero_model_one.validate().ok() &&
                          zero_model_three.validate().ok(),
                      "zero-width SOURCE configurations are valid");
  const swcme1d::StepState zero_state_one=
      zero_model_one.prepare_step(3600.0);
  const swcme3d::StepState zero_state_three=
      zero_model_three.prepare_step(3600.0);
  context.expect_true(zero_state_one.w_sh_m==0.0 &&
                          zero_state_one.w_le_m==0.0 &&
                          zero_state_one.w_te_m==0.0 &&
                          zero_state_three.w_shock_m==0.0 &&
                          zero_state_three.w_le_m==0.0 &&
                          zero_state_three.w_te_m==0.0,
                      "zero requests produce exactly zero effective widths");

  // Equality at each declared 90-percent limit is valid.  Exercising all
  // limits together is the most restrictive non-overlapping configuration and
  // proves that the comparisons are inclusive rather than approximate.
  swcme1d::Params edge_one=base_one;
  swcme3d::Params edge_three=base_three;
  edge_one.edge_smooth_shock_AU_at1AU=limits.shock;
  edge_one.edge_smooth_le_AU_at1AU=limits.leading;
  edge_one.edge_smooth_te_AU_at1AU=limits.trailing;
  edge_three.edge_smooth_shock_AU_at1AU=limits.shock;
  edge_three.edge_smooth_le_AU_at1AU=limits.leading;
  edge_three.edge_smooth_te_AU_at1AU=limits.trailing;
  const swcme1d::Model edge_model_one(edge_one);
  const swcme3d::Model edge_model_three(edge_three);
  context.expect_true(edge_model_one.validate().ok() &&
                          edge_model_three.validate().ok(),
                      "exact smoothing limits are accepted in both models");
  const swcme1d::StepState edge_state_one=
      edge_model_one.prepare_step(3600.0);
  const swcme3d::StepState edge_state_three=
      edge_model_three.prepare_step(3600.0);
  expect_exact_effective_widths(context,edge_state_one.region_config,
                                edge_state_one.region_boundaries,
                                edge_state_one.r_sh_m,"1-D exact-limit state");
  expect_exact_effective_widths(context,edge_state_three.region_config,
                                edge_state_three.apex_regions,
                                edge_state_three.r_sh_m,"3-D apex exact-limit state");

  // The same accepted fractions must scale without alteration at an arbitrary
  // local 3-D flank radius, not only at the cached apex diagnostic.
  const double local_radius=0.73*edge_state_three.r_sh_m;
  const swcme::regions::Boundaries local_boundaries=
      swcme::regions::make_boundaries(
          local_radius,edge_state_three.region_config);
  expect_exact_effective_widths(context,edge_state_three.region_config,
                                local_boundaries,local_radius,
                                "3-D local exact-limit state");

  // At the simultaneous limits the transition intervals retain finite gaps;
  // their centers must still produce the canonical smoothstep value 1/2.
  const swcme::regions::Boundaries& b=edge_state_one.region_boundaries;
  context.expect_true(
      b.R_sh_m-0.5*b.smooth_shock_width_m>
          b.R_le_m+0.5*b.smooth_le_width_m &&
      b.R_le_m-0.5*b.smooth_le_width_m>
          b.R_te_m+0.5*b.smooth_te_width_m,
      "exact-limit transition intervals remain separated");
  const swcme::regions::Location shock_center=
      swcme::regions::locate(b.R_sh_m,b);
  const swcme::regions::Location leading_center=
      swcme::regions::locate(b.R_le_m,b);
  const swcme::regions::Location trailing_center=
      swcme::regions::locate(b.R_te_m,b);
  context.expect_true(shock_center.region==
                          swcme::regions::Region::ShockTransition &&
                          std::isfinite(shock_center.blend) &&
                          std::abs(shock_center.blend-0.5)<=
                              8.0*std::numeric_limits<double>::epsilon(),
                      "shock center has finite one-half blend");
  context.expect_true(leading_center.region==
                          swcme::regions::Region::LeadingTransition &&
                          std::isfinite(leading_center.blend) &&
                          std::abs(leading_center.blend-0.5)<=
                              8.0*std::numeric_limits<double>::epsilon(),
                      "leading center has finite one-half blend");
  context.expect_true(trailing_center.region==
                          swcme::regions::Region::TrailingTransition &&
                          std::isfinite(trailing_center.blend) &&
                          std::abs(trailing_center.blend-0.5)<=
                              8.0*std::numeric_limits<double>::epsilon(),
                      "trailing center has finite one-half blend");

  // One ULP above each limit must fail independently.  Using nextafter avoids
  // a tolerance-dependent test and protects the exact public boundary.
  {
    swcme1d::Params one=base_one;
    swcme3d::Params three=base_three;
    one.edge_smooth_shock_AU_at1AU=
        std::nextafter(limits.shock,std::numeric_limits<double>::infinity());
    three.edge_smooth_shock_AU_at1AU=one.edge_smooth_shock_AU_at1AU;
    expect_rejected_pair(context,one,three,
        "edge_smooth_shock_AU_at1AU",swcme::config::Code::OutOfRange,
        "one-ULP-above shock limit");
  }
  {
    swcme1d::Params one=base_one;
    swcme3d::Params three=base_three;
    one.edge_smooth_le_AU_at1AU=
        std::nextafter(limits.leading,std::numeric_limits<double>::infinity());
    three.edge_smooth_le_AU_at1AU=one.edge_smooth_le_AU_at1AU;
    expect_rejected_pair(context,one,three,
        "edge_smooth_le_AU_at1AU",swcme::config::Code::OutOfRange,
        "one-ULP-above leading limit");
  }
  {
    swcme1d::Params one=base_one;
    swcme3d::Params three=base_three;
    one.edge_smooth_te_AU_at1AU=
        std::nextafter(limits.trailing,std::numeric_limits<double>::infinity());
    three.edge_smooth_te_AU_at1AU=one.edge_smooth_te_AU_at1AU;
    expect_rejected_pair(context,one,three,
        "edge_smooth_te_AU_at1AU",swcme::config::Code::OutOfRange,
        "one-ULP-above trailing limit");
  }

  // Non-finite values retain the pre-existing NON_FINITE classification rather
  // than being confused with a finite overlap violation.
  {
    swcme1d::Params one=base_one;
    swcme3d::Params three=base_three;
    one.edge_smooth_shock_AU_at1AU=
        std::numeric_limits<double>::quiet_NaN();
    three.edge_smooth_shock_AU_at1AU=one.edge_smooth_shock_AU_at1AU;
    expect_rejected_pair(context,one,three,
        "edge_smooth_shock_AU_at1AU",swcme::config::Code::NonFinite,
        "non-finite shock width");
  }
  {
    swcme1d::Params one=base_one;
    swcme3d::Params three=base_three;
    one.edge_smooth_le_AU_at1AU=std::numeric_limits<double>::infinity();
    three.edge_smooth_le_AU_at1AU=one.edge_smooth_le_AU_at1AU;
    expect_rejected_pair(context,one,three,
        "edge_smooth_le_AU_at1AU",swcme::config::Code::NonFinite,
        "non-finite leading width");
  }
  {
    swcme1d::Params one=base_one;
    swcme3d::Params three=base_three;
    one.edge_smooth_te_AU_at1AU=-
        std::numeric_limits<double>::infinity();
    three.edge_smooth_te_AU_at1AU=one.edge_smooth_te_AU_at1AU;
    expect_rejected_pair(context,one,three,
        "edge_smooth_te_AU_at1AU",swcme::config::Code::NonFinite,
        "non-finite trailing width");
  }

  // A multi-field conflict must report every width in one validation pass;
  // fixing only the first input and retrying would make campaign setup brittle.
  swcme1d::Params conflict_one=base_one;
  swcme3d::Params conflict_three=base_three;
  conflict_one.edge_smooth_shock_AU_at1AU=
      std::nextafter(limits.shock,std::numeric_limits<double>::infinity());
  conflict_one.edge_smooth_le_AU_at1AU=
      std::nextafter(limits.leading,std::numeric_limits<double>::infinity());
  conflict_one.edge_smooth_te_AU_at1AU=
      std::nextafter(limits.trailing,std::numeric_limits<double>::infinity());
  conflict_three.edge_smooth_shock_AU_at1AU=
      conflict_one.edge_smooth_shock_AU_at1AU;
  conflict_three.edge_smooth_le_AU_at1AU=
      conflict_one.edge_smooth_le_AU_at1AU;
  conflict_three.edge_smooth_te_AU_at1AU=
      conflict_one.edge_smooth_te_AU_at1AU;
  const swcme::config::ValidationResult conflict_result_one=
      swcme1d::Model(conflict_one).validate();
  const swcme::config::ValidationResult conflict_result_three=
      swcme3d::Model(conflict_three).validate();
  context.expect_true(conflict_result_one.issues.size()==3U &&
                          conflict_result_three.issues.size()==3U,
                      "multi-width conflicts report all three fields");
  context.expect_true(
      has_issue(conflict_result_one,"edge_smooth_shock_AU_at1AU",
                swcme::config::Code::OutOfRange) &&
      has_issue(conflict_result_one,"edge_smooth_le_AU_at1AU",
                swcme::config::Code::OutOfRange) &&
      has_issue(conflict_result_one,"edge_smooth_te_AU_at1AU",
                swcme::config::Code::OutOfRange) &&
      has_issue(conflict_result_three,"edge_smooth_shock_AU_at1AU",
                swcme::config::Code::OutOfRange) &&
      has_issue(conflict_result_three,"edge_smooth_le_AU_at1AU",
                swcme::config::Code::OutOfRange) &&
      has_issue(conflict_result_three,"edge_smooth_te_AU_at1AU",
                swcme::config::Code::OutOfRange),
      "multi-width conflict diagnostics remain field-specific");

  // Inactive values are still part of the stored configuration and its digest.
  // The selected project policy therefore rejects an oversized dormant width
  // in SHOCK_ONLY/SOURCE just as it does in FULL_ICME; changing mode during the
  // pre-prepare setup phase cannot reveal a previously accepted bad value.
  swcme1d::Params inactive_one;
  swcme3d::Params inactive_three;
  const swcme::regions::SmoothingFractionLimits inactive_limits=
      swcme::regions::smoothing_fraction_limits(
          inactive_one.sheath_thick_AU_at1AU,
          inactive_one.ejecta_thick_AU_at1AU);
  inactive_one.edge_smooth_le_AU_at1AU=std::nextafter(
      inactive_limits.leading,std::numeric_limits<double>::infinity());
  inactive_three.edge_smooth_le_AU_at1AU=
      inactive_one.edge_smooth_le_AU_at1AU;
  expect_rejected_pair(context,inactive_one,inactive_three,
      "edge_smooth_le_AU_at1AU",swcme::config::Code::OutOfRange,
      "inactive one-ULP-above leading limit");
}
