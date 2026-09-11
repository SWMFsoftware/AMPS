#include "test_framework.hpp"

#include <swcme1d.hpp>
#include <swcme3d.hpp>

#include <array>
#include <cstring>
#include <iostream>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

namespace {

// Compile-time API inspection documents the two supported compatibility
// surfaces: neither dimension may expose a mutable Params reference.  This is
// essential because a reference acquired before preparation could otherwise
// be retained and used after the model freezes.  GetParams() remains read-only.
template <typename T, typename = void>
struct has_mutable_params : std::false_type {};

template <typename T>
struct has_mutable_params<
    T,std::void_t<decltype(std::declval<T&>().MutableParams())>>
    : std::true_type {};

static_assert(!has_mutable_params<swcme3d::Model>::value,
              "3-D model parameters must remain construction-only");
static_assert(!has_mutable_params<swcme1d::Model>::value,
              "1-D model must not expose a retainable mutable Params alias");
static_assert(std::is_same<
                  decltype(std::declval<const swcme1d::Model&>().GetParams()),
                  const swcme1d::Params&>::value,
              "1-D read accessor must not expose mutable parameters");

using OneSnapshot = std::array<double,18>;
using ThreeSnapshot = std::array<double,24>;

// Serialize the complete 1-D checked output at three fixed radii into one
// contiguous double array.  memcmp on this representation enforces the
// validation plan's preferred bitwise-identity criterion.
OneSnapshot capture_1d(swcme_test::Context& context,
                       const swcme1d::Model& model,
                       const swcme1d::StepState& state,
                       const std::string& label) {
  const double radii[3]={0.5*swcme::constants::AU_M,
                         swcme::constants::AU_M,
                         1.5*swcme::constants::AU_M};
  double density[3]={},velocity[3]={},Br[3]={},Bphi[3]={},Bmag[3]={},divV[3]={};
  const swcme::ModelStatus status=model.evaluate_radii_with_B_div_checked(
      state,radii,density,velocity,Br,Bphi,Bmag,divV,3);
  context.expect_true(status.ok(),label+" evaluates successfully");

  OneSnapshot snapshot{};
  std::size_t offset=0;
  for (const double* values : {density,velocity,Br,Bphi,Bmag,divV})
    for (std::size_t i=0; i<3; ++i) snapshot[offset++]=values[i];
  return snapshot;
}

// Serialize density, velocity, magnetic field, and divergence at three 3-D
// points.  The layout is deliberately fixed and contains no struct padding.
ThreeSnapshot capture_3d(swcme_test::Context& context,
                         const swcme3d::Model& model,
                         const swcme3d::StepState& state,
                         const std::string& label) {
  const double AU=swcme::constants::AU_M;
  const double x[3]={0.5*AU,0.0,1.5*AU};
  const double y[3]={0.0,AU,0.0};
  const double z[3]={0.0,0.0,0.2*AU};
  double density[3]={},vx[3]={},vy[3]={},vz[3]={};
  double bx[3]={},by[3]={},bz[3]={},divV[3]={};
  const swcme::ModelStatus status=model.evaluate_cartesian_with_B_div_checked(
      state,x,y,z,density,vx,vy,vz,bx,by,bz,divV,3);
  context.expect_true(status.ok(),label+" evaluates successfully");

  ThreeSnapshot snapshot{};
  std::size_t offset=0;
  for (const double* values : {density,vx,vy,vz,bx,by,bz,divV})
    for (std::size_t i=0; i<3; ++i) snapshot[offset++]=values[i];
  return snapshot;
}

template <std::size_t N>
bool bitwise_equal(const std::array<double,N>& left,
                   const std::array<double,N>& right) {
  return std::memcmp(left.data(),right.data(),sizeof(double)*N)==0;
}

// Every legacy mutator uses the same explicit logic_error contract.  Checking
// the stable diagnostic prevents a future setter from silently becoming a
// special case that reopens an already prepared model.
template <typename Mutation>
void expect_locked_mutation(swcme_test::Context& context, Mutation&& mutation,
                            const std::string& label) {
  bool caught=false;
  std::string diagnostic;
  try {
    mutation();
  } catch (const std::logic_error& error) {
    caught=true;
    diagnostic=error.what();
  }
  context.expect_true(caught,label+" is rejected after preparation");
  context.expect_true(
      diagnostic.find("immutable after successful prepare_step")!=
          std::string::npos,
      label+" explains the immutable lifecycle");
}

}  // namespace

// PST01: a successfully prepared state cannot be made dependent on later
// model-configuration changes.  All retained 1-D mutators and both dimensions'
// assignment paths are rejected after preparation; replacement configuration
// constructs a new owner.  Original-state results remain bitwise identical.
void test_pst01(swcme_test::Context& context) {
  std::cout << "PST01 prepared-state immutability\n";

  // Configure through the legacy fluent API before preparation to preserve
  // source compatibility for existing callers' setup phase.
  swcme1d::Model one;
  one.SetAmbient(400.0,6.0,5.0,1.2e5,5.0/3.0,1.0)
     .SetCME(40.0,1200.0,0.0)
     .SetKinematicsMode(swcme::kinematics::Mode::Ballistic)
     .SetRegionMode(swcme::regions::Mode::ShockOnly)
     .SetShockAccelerationMode(swcme::acceleration::Mode::Source)
     .SetGeometry(0.10,0.20)
     .SetSmoothing(0.01,0.02,0.03)
     .SetSheathEjecta(1.0,2.0,1.10,0.5,0.8);
  const swcme1d::StepState one_state=one.prepare_step(7200.0);
  context.expect_true(one.configuration_locked(),
                      "successful 1-D preparation freezes configuration");
  const swcme::ModelIdentity one_identity=one.model_identity();
  const swcme::ConfigurationDigest one_digest=
      swcme1d::configuration_digest(one.GetParams());
  const OneSnapshot one_baseline=
      capture_1d(context,one,one_state,"1-D baseline");

  swcme1d::Params replacement_params=one.GetParams();
  replacement_params.sin_theta=0.0;
  swcme1d::Model assignment_source(replacement_params);

  // Exercise every public 1-D mutation entry point that predates PST01.
  // Each attempt is followed by a digest/result check, proving rejection took
  // place before even the first configuration field was modified.
  const auto verify_unchanged=[&](const std::string& label) {
    context.expect_true(one.model_identity()==one_identity,
                        label+" preserves model identity");
    context.expect_true(swcme1d::configuration_digest(one.GetParams())==one_digest,
                        label+" preserves configuration digest");
    context.expect_true(bitwise_equal(
                            one_baseline,capture_1d(context,one,one_state,label)),
                        label+" preserves bitwise state results");
  };

  expect_locked_mutation(context,[&]{ one.SetParams(replacement_params); },
                         "SetParams");
  verify_unchanged("SetParams rejection");
  expect_locked_mutation(context,[&]{ one.SetCME(30.0,900.0,1.0e-8); },
                         "SetCME");
  verify_unchanged("SetCME rejection");
  expect_locked_mutation(context,[&]{
    one.SetKinematicsMode(swcme::kinematics::Mode::DBM);
  },"SetKinematicsMode");
  verify_unchanged("SetKinematicsMode rejection");
  expect_locked_mutation(context,[&]{
    one.SetDataDrivenKinematics({0.0,1.0},{40.0,40.1});
  },"SetDataDrivenKinematics");
  verify_unchanged("SetDataDrivenKinematics rejection");
  expect_locked_mutation(context,[&]{
    one.SetAmbient(450.0,7.0,6.0,1.5e5,1.5,0.0);
  },"SetAmbient/sin_theta");
  verify_unchanged("SetAmbient rejection");
  expect_locked_mutation(context,[&]{
    one.SetRegionMode(swcme::regions::Mode::FullICME);
  },"SetRegionMode");
  verify_unchanged("SetRegionMode rejection");
  expect_locked_mutation(context,[&]{
    one.SetShockAccelerationMode(
        swcme::acceleration::Mode::ResolvedCompression);
  },"SetShockAccelerationMode");
  verify_unchanged("SetShockAccelerationMode rejection");
  expect_locked_mutation(context,[&]{ one.SetGeometry(0.12,0.24); },
                         "SetGeometry");
  verify_unchanged("SetGeometry rejection");
  expect_locked_mutation(context,[&]{ one.SetSmoothing(0.02,0.03,0.04); },
                         "SetSmoothing");
  verify_unchanged("SetSmoothing rejection");
  expect_locked_mutation(context,[&]{
    one.SetSheathEjecta(1.1,2.5,1.15,0.6,0.9);
  },"SetSheathEjecta");
  verify_unchanged("SetSheathEjecta rejection");
  expect_locked_mutation(context,[&]{ one=assignment_source; },
                         "1-D copy assignment");
  verify_unchanged("1-D assignment rejection");

  // The supported replacement path creates a new identity/configuration and
  // can be prepared independently without changing the original snapshot.
  swcme1d::Model one_replacement=one.reconfigured(replacement_params);
  context.expect_true(!one_replacement.configuration_locked(),
                      "1-D replacement begins in configuration phase");
  const swcme1d::StepState replacement_state=
      one_replacement.prepare_step(7200.0);
  const OneSnapshot replacement_output=capture_1d(
      context,one_replacement,replacement_state,"1-D replacement");
  context.expect_true(one_replacement.model_identity()!=one_identity,
                      "1-D replacement receives a new owner identity");
  context.expect_true(replacement_state.configuration_digest!=
                          one_state.configuration_digest,
                      "1-D replacement records its changed configuration");
  context.expect_true(!bitwise_equal(one_baseline,replacement_output),
                      "1-D replacement may produce different physics");
  verify_unchanged("1-D replacement creation");

  // A failed preparation must not freeze the setup phase because it never
  // issued a valid state.  The caller can repair the configuration and retry.
  swcme1d::Params invalid_params=replacement_params;
  invalid_params.V_sw_kms=0.0;
  swcme1d::Model repairable(invalid_params);
  bool invalid_rejected=false;
  try {
    (void)repairable.prepare_step(0.0);
  } catch (const std::invalid_argument&) {
    invalid_rejected=true;
  }
  context.expect_true(invalid_rejected && !repairable.configuration_locked(),
                      "failed preparation leaves 1-D configuration repairable");
  repairable.SetAmbient(400.0,6.0,5.0,1.2e5,5.0/3.0,1.0);
  (void)repairable.prepare_step(0.0);
  context.expect_true(repairable.configuration_locked(),
                      "repaired successful preparation freezes 1-D model");

  // Three-dimensional Params are already private and construction-only.  Test
  // the remaining assignment route plus the same replacement-builder rule.
  swcme3d::Params three_params;
  three_params.shape=swcme3d::ShockShape::Sphere;
  three_params.kinematics_mode=swcme::kinematics::Mode::Ballistic;
  three_params.r0_Rs=40.0;
  three_params.V0_sh_kms=1200.0;
  three_params.V_sw_kms=400.0;
  three_params.region_mode=swcme::regions::Mode::ShockOnly;
  three_params.shock_acceleration_mode=swcme::acceleration::Mode::Source;
  swcme3d::Model three(three_params);
  const swcme3d::StepState three_state=three.prepare_step(7200.0);
  const ThreeSnapshot three_baseline=
      capture_3d(context,three,three_state,"3-D baseline");
  const swcme::ModelIdentity three_identity=three.model_identity();
  context.expect_true(three.configuration_locked(),
                      "successful 3-D preparation freezes configuration");

  swcme3d::Params three_changed=three_params;
  three_changed.B1AU_nT+=2.0;
  swcme3d::Model three_assignment_source(three_changed);
  expect_locked_mutation(context,[&]{ three=three_assignment_source; },
                         "3-D copy assignment");
  context.expect_true(three.model_identity()==three_identity,
                      "rejected 3-D assignment preserves identity");
  context.expect_true(bitwise_equal(
                          three_baseline,
                          capture_3d(context,three,three_state,
                                     "3-D assignment rejection")),
                      "rejected 3-D assignment preserves bitwise results");

  swcme3d::Model three_replacement=three.reconfigured(three_changed);
  const swcme3d::StepState three_replacement_state=
      three_replacement.prepare_step(7200.0);
  const ThreeSnapshot three_replacement_output=capture_3d(
      context,three_replacement,three_replacement_state,"3-D replacement");
  context.expect_true(three_replacement.model_identity()!=three_identity,
                      "3-D replacement receives a new owner identity");
  context.expect_true(three_replacement_state.configuration_digest!=
                          three_state.configuration_digest,
                      "3-D replacement records its changed configuration");
  context.expect_true(!bitwise_equal(three_baseline,three_replacement_output),
                      "3-D replacement may produce different physics");
  context.expect_true(bitwise_equal(
                          three_baseline,
                          capture_3d(context,three,three_state,
                                     "3-D replacement creation")),
                      "3-D replacement leaves original state bitwise identical");
}
