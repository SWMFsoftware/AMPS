#include "test_framework.hpp"

#include <swcme1d.hpp>
#include <swcme3d.hpp>
#include <swcme_config.hpp>
#include <swcme_solarwind.hpp>

#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

bool has_issue(const swcme::config::ValidationResult& result,
               const std::string& field,swcme::config::Code code) {
  for (const swcme::config::Issue& issue : result.issues) {
    if (issue.field==field && issue.code==code) return true;
  }
  return false;
}

// Configuration validity is a pre-preparation contract.  These overloads
// verify that both public dimensional wrappers route prepare_step() through
// the same CFG04 validator instead of discovering a bad radius in later
// solar-wind, geometry, or interpolation code.
bool prepare_rejects(const swcme1d::Params& params) {
  try {
    const swcme1d::Model model(params);
    (void)model.prepare_step(0.0);
  } catch (const std::invalid_argument&) {
    return true;
  }
  return false;
}

bool prepare_rejects(const swcme3d::Params& params) {
  try {
    const swcme3d::Model model(params);
    (void)model.prepare_step(0.0);
  } catch (const std::invalid_argument&) {
    return true;
  }
  return false;
}

void expect_invalid_pair(swcme_test::Context& context,
                         const swcme1d::Params& one,
                         const swcme3d::Params& three,
                         const std::string& field,
                         swcme::config::Code code,
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

void configure_data_driven(swcme1d::Params& one,swcme3d::Params& three,
                           const std::vector<double>& times,
                           const std::vector<double>& radii) {
  one.kinematics_mode=swcme::kinematics::Mode::DataDriven;
  one.data_time_s=times;
  one.data_radius_Rs=radii;
  three.kinematics_mode=one.kinematics_mode;
  three.data_time_s=times;
  three.data_radius_Rs=radii;
}

}  // namespace

void test_cfg04(swcme_test::Context& context) {
  std::cout << "CFG04 configured radius-domain validation\n";

  const double boundary_rs=swcme::solarwind::MIN_RADIUS_RS;
  const double below_rs=std::nextafter(boundary_rs,0.0);
  const double above_rs=
      std::nextafter(boundary_rs,std::numeric_limits<double>::infinity());

  // The model reference radius uses solar-radius input units.  One ULP below
  // the analytical boundary must fail with the public field name, while exact
  // equality and one ULP above must prepare without radius repair.
  swcme1d::Params r0_one;
  swcme3d::Params r0_three;
  r0_one.r0_Rs=below_rs;
  r0_three.r0_Rs=below_rs;
  expect_invalid_pair(context,r0_one,r0_three,"r0_Rs",
                      swcme::config::Code::OutOfRange,
                      "sub-domain reference radius");

  r0_one.r0_Rs=boundary_rs;
  r0_three.r0_Rs=boundary_rs;
  const swcme1d::StepState exact_r0_one=
      swcme1d::Model(r0_one).prepare_step(0.0);
  const swcme3d::StepState exact_r0_three=
      swcme3d::Model(r0_three).prepare_step(0.0);
  context.expect_true(exact_r0_one.r_sh_m==
                          swcme::solarwind::MIN_RADIUS_M &&
                          exact_r0_three.r_sh_m==
                          swcme::solarwind::MIN_RADIUS_M,
                      "exact r0 boundary prepares without clipping");

  r0_one.r0_Rs=above_rs;
  r0_three.r0_Rs=above_rs;
  context.expect_true(swcme1d::Model(r0_one).validate().ok() &&
                          swcme3d::Model(r0_three).validate().ok(),
                      "one-ULP-above r0 boundary is valid");

  r0_one.r0_Rs=std::numeric_limits<double>::quiet_NaN();
  r0_three.r0_Rs=r0_one.r0_Rs;
  expect_invalid_pair(context,r0_one,r0_three,"r0_Rs",
                      swcme::config::Code::NonFinite,
                      "non-finite reference radius");

  // Every PCHIP radius is independently domain checked.  Put the one-ULP
  // violation at each possible position so validation cannot accidentally
  // inspect only the first or last knot.
  for (std::size_t bad_index=0; bad_index<3; ++bad_index) {
    std::vector<double> radii={boundary_rs,above_rs,2.0};
    radii[bad_index]=below_rs;
    swcme1d::Params one;
    swcme3d::Params three;
    configure_data_driven(one,three,{0.0,1.0,2.0},radii);
    expect_invalid_pair(
        context,one,three,
        "data_radius_Rs["+std::to_string(bad_index)+"]",
        swcme::config::Code::OutOfRange,
        "sub-domain PCHIP knot "+std::to_string(bad_index));
  }

  // Exact-boundary and one-ULP-above knots form a valid nondecreasing table.
  // Preparing at the first knot must retain the exact supplied boundary and
  // therefore demonstrates that validation does not clamp or perturb it.
  swcme1d::Params data_one;
  swcme3d::Params data_three;
  configure_data_driven(data_one,data_three,{0.0,1.0},
                        {boundary_rs,above_rs});
  const swcme1d::StepState data_state_one=
      swcme1d::Model(data_one).prepare_step(0.0);
  const swcme3d::StepState data_state_three=
      swcme3d::Model(data_three).prepare_step(0.0);
  context.expect_true(data_state_one.r_sh_m==
                          swcme::solarwind::MIN_RADIUS_M &&
                          data_state_three.r_sh_m==
                          swcme::solarwind::MIN_RADIUS_M,
                      "exact PCHIP domain boundary prepares unchanged");

  // Mixed invalidity must report each malformed radius by knot index in one
  // pass.  Ordering remains a separate table-level diagnostic because it is a
  // relationship rather than a property of one radius.
  configure_data_driven(
      data_one,data_three,{0.0,1.0,1.0,3.0},
      {below_rs,boundary_rs,
       std::numeric_limits<double>::infinity(),above_rs});
  const swcme::config::ValidationResult mixed_one=
      swcme1d::Model(data_one).validate();
  const swcme::config::ValidationResult mixed_three=
      swcme3d::Model(data_three).validate();
  for (const swcme::config::ValidationResult* result :
       {&mixed_one,&mixed_three}) {
    context.expect_true(
        has_issue(*result,"data_radius_Rs[0]",
                  swcme::config::Code::OutOfRange),
        "mixed table reports its sub-domain radius index");
    context.expect_true(
        has_issue(*result,"data_radius_Rs[2]",
                  swcme::config::Code::NonFinite),
        "mixed table reports its non-finite radius index");
    context.expect_true(
        has_issue(*result,"data_time_s",
                  swcme::config::Code::InvalidKinematicsTable),
        "mixed table reports repeated-time ordering");
  }
  context.expect_true(prepare_rejects(data_one) && prepare_rejects(data_three),
                      "mixed-invalid PCHIP tables cannot prepare");

  // ConnectivityOptions::inner_radius_m is the configured source/search
  // surface.  Its SI boundary is inclusive and it must not exceed the observer.
  // CON10 defines equality as a one-point closed-interval query, leaving only
  // reversed ordering as an invalid connectivity configuration.
  const swcme3d::Model connectivity_model(swcme3d::Params{});
  const swcme3d::StepState connectivity_state=
      connectivity_model.prepare_step(0.0);
  swcme3d::ConnectivityOptions options;
  const double boundary_m=swcme::solarwind::MIN_RADIUS_M;
  const double below_m=std::nextafter(boundary_m,0.0);
  const double above_m=
      std::nextafter(boundary_m,std::numeric_limits<double>::infinity());
  const double observer[3]={2.0*boundary_m,0.0,0.0};

  options.inner_radius_m=below_m;
  context.expect_true(
      connectivity_model.observer_connectivity(
          connectivity_state,observer,options).status==
          swcme3d::ConnectivityStatus::InvalidConfiguration,
      "one-ULP-below source surface is invalid configuration");
  options.inner_radius_m=boundary_m;
  context.expect_true(
      connectivity_model.observer_connectivity(
          connectivity_state,observer,options).status!=
          swcme3d::ConnectivityStatus::InvalidConfiguration,
      "exact source-surface boundary is domain-valid");
  options.inner_radius_m=above_m;
  context.expect_true(
      connectivity_model.observer_connectivity(
          connectivity_state,observer,options).status!=
          swcme3d::ConnectivityStatus::InvalidConfiguration,
      "one-ULP-above source-surface boundary is domain-valid");

  options.inner_radius_m=observer[0];
  context.expect_true(
      connectivity_model.observer_connectivity(
          connectivity_state,observer,options).status!=
          swcme3d::ConnectivityStatus::InvalidConfiguration,
      "source surface equal to observer permits one-point closed interval");
  options.inner_radius_m=std::nextafter(
      observer[0],std::numeric_limits<double>::infinity());
  context.expect_true(
      connectivity_model.observer_connectivity(
          connectivity_state,observer,options).status==
          swcme3d::ConnectivityStatus::InvalidConfiguration,
      "source surface outside observer rejects inverted radial domain");

  // Observer scope and connectivity must agree on the lower model boundary.
  // Exact equality is valid for both scope membership and the connectivity
  // solver's one-point closed interval.
  const swcme1d::Model scope_one;
  const swcme1d::StepState scope_state_one=scope_one.prepare_step(0.0);
  const double observer_below[3]={below_m,0.0,0.0};
  const double observer_exact[3]={boundary_m,0.0,0.0};
  const double observer_above[3]={above_m,0.0,0.0};
  context.expect_true(
      !scope_one.observer_scope_status(scope_state_one,below_m).valid_input &&
      !connectivity_model.observer_scope_status(
          connectivity_state,observer_below).valid_input,
      "one-ULP-below observer is outside both scope APIs");
  context.expect_true(
      scope_one.observer_scope_status(scope_state_one,boundary_m).valid_input &&
      connectivity_model.observer_scope_status(
          connectivity_state,observer_exact).valid_input,
      "exact observer boundary is valid in both scope APIs");
  context.expect_true(
      scope_one.observer_scope_status(scope_state_one,above_m).valid_input &&
      connectivity_model.observer_scope_status(
          connectivity_state,observer_above).valid_input,
      "one-ULP-above observer is valid in both scope APIs");

  options.inner_radius_m=boundary_m;
  context.expect_true(
      connectivity_model.observer_connectivity(
          connectivity_state,observer_below,options).status==
          swcme3d::ConnectivityStatus::InvalidObserver,
      "sub-domain observer is not reported as disconnected");
  context.expect_true(
      connectivity_model.observer_connectivity(
          connectivity_state,observer_exact,options).status!=
          swcme3d::ConnectivityStatus::InvalidConfiguration,
      "boundary observer permits exact one-point connectivity query");
  context.expect_true(
      connectivity_model.observer_connectivity(
          connectivity_state,observer_above,options).status!=
          swcme3d::ConnectivityStatus::InvalidObserver,
      "one-ULP-above observer passes the model-domain check");
}
