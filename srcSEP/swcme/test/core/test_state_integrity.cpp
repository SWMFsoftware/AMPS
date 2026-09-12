#include "test_framework.hpp"

#include <swcme1d.hpp>
#include <swcme3d.hpp>

#include <functional>
#include <iostream>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

namespace {

// The public integrity accessor returns a value, never a mutable reference.
// These compile-time checks prevent accidental API changes that would let
// compatibility code rewrite a seal after modifying cached fields.
static_assert(!std::is_reference<decltype(
                  std::declval<swcme1d::StepState&>().integrity_digest())>::value,
              "1-D prepared-state seal must not be mutable");
static_assert(!std::is_reference<decltype(
                  std::declval<swcme3d::StepState&>().integrity_digest())>::value,
              "3-D prepared-state seal must not be mutable");
static_assert(std::is_copy_constructible<swcme1d::StepState>::value &&
                  std::is_move_constructible<swcme1d::StepState>::value,
              "1-D prepared state must remain copy/move constructible");
static_assert(std::is_copy_constructible<swcme3d::StepState>::value &&
                  std::is_move_constructible<swcme3d::StepState>::value,
              "3-D prepared state must remain copy/move constructible");

using Mutation1D = std::pair<
    std::string,std::function<void(swcme1d::StepState&)>>;
using Mutation3D = std::pair<
    std::string,std::function<void(swcme3d::StepState&)>>;

// Check the complete PST06 status contract and the required transactional
// behavior of a checked 1-D consumer.
void reject_1d_corruption(swcme_test::Context& context,
                          const swcme1d::Model& model,
                          const swcme1d::StepState& baseline,
                          const Mutation1D& mutation) {
  swcme1d::StepState corrupted=baseline;
  mutation.second(corrupted);
  const double radius=swcme::constants::AU_M;
  double density=801.0,velocity=802.0;
  const swcme::ModelStatus status=model.evaluate_radii_fast_checked(
      corrupted,&radius,&density,&velocity,1);
  context.expect_true(status.code==swcme::StatusCode::StalePreparedState,
                      mutation.first+" reports STALE_PREPARED_STATE");
  context.expect_true(status.has_state_integrity &&
                          status.expected_state_integrity==
                              baseline.integrity_digest() &&
                          status.computed_state_integrity!=
                              baseline.integrity_digest(),
                      mutation.first+" carries seal diagnostics");
  context.expect_true(density==801.0 && velocity==802.0,
                      mutation.first+" preserves 1-D outputs");
}

// Apply the same contract to a 3-D checked field call.  Sentinels establish
// that validation happens before any sample result is allocated or written.
void reject_3d_corruption(swcme_test::Context& context,
                          const swcme3d::Model& model,
                          const swcme3d::StepState& baseline,
                          const Mutation3D& mutation) {
  swcme3d::StepState corrupted=baseline;
  mutation.second(corrupted);
  const double x=swcme::constants::AU_M,y=0.0,z=0.0;
  double density=811.0,vx=812.0,vy=813.0,vz=814.0;
  const swcme::ModelStatus status=model.evaluate_cartesian_fast_checked(
      corrupted,&x,&y,&z,&density,&vx,&vy,&vz,1);
  context.expect_true(status.code==swcme::StatusCode::StalePreparedState,
                      mutation.first+" reports STALE_PREPARED_STATE");
  context.expect_true(status.has_state_integrity &&
                          status.expected_state_integrity==
                              baseline.integrity_digest() &&
                          status.computed_state_integrity!=
                              baseline.integrity_digest(),
                      mutation.first+" carries seal diagnostics");
  context.expect_true(
      density==811.0 && vx==812.0 && vy==813.0 && vz==814.0,
      mutation.first+" preserves 3-D outputs");
}

}  // namespace

// PST06: prepared records retain public compatibility mirrors, but every such
// field is covered by a private seal.  Representative nested canonical fields
// and every top-level legacy mirror are corrupted independently.  Valid copy/
// move operations preserve provenance, seal, and evaluability.
void test_pst06(swcme_test::Context& context) {
  std::cout << "PST06 prepared-state record integrity\n";

  swcme1d::Params p1;
  p1.kinematics_mode=swcme::kinematics::Mode::Ballistic;
  p1.r0_Rs=40.0;
  p1.V0_sh_kms=1200.0;
  p1.V_sw_kms=400.0;
  p1.region_mode=swcme::regions::Mode::ShockOnly;
  p1.shock_acceleration_mode=swcme::acceleration::Mode::Source;
  swcme1d::Model one(p1);
  const swcme1d::StepState one_state=one.prepare_step(3600.0);
  context.expect_true(one_state.integrity_digest()!=0 &&
                          swcme1d::prepared_state_integrity(one_state)==
                              one_state.integrity_digest(),
                      "1-D prepared record has a valid private seal");

  // Copy/move construction and assignment are intentionally supported because
  // AMPS can distribute immutable state values.  All forms must preserve the
  // exact seal and remain accepted by the originating model.
  swcme1d::StepState one_copy(one_state);
  swcme1d::StepState one_moved(std::move(one_copy));
  swcme1d::StepState one_copy_assigned;
  one_copy_assigned=one_state;
  swcme1d::StepState one_move_assigned;
  one_move_assigned=std::move(one_copy_assigned);
  const double radius=swcme::constants::AU_M;
  for (const auto* state : {&one_moved,&one_move_assigned}) {
    double density=0.0,velocity=0.0;
    const swcme::ModelStatus status=one.evaluate_radii_fast_checked(
        *state,&radius,&density,&velocity,1);
    context.expect_true(status.ok() &&
                            state->owner_model_identity==one.model_identity() &&
                            state->integrity_digest()==one_state.integrity_digest(),
                        "valid 1-D copy/move preserves record integrity");
  }

  // Each top-level 1-D compatibility mirror is listed explicitly so a future
  // field addition is visible during test review.  Nested canonical records
  // are sampled separately to prove they are sealed as well.
  const std::vector<Mutation1D> one_mutations={
      {"1-D common solar wind",[](auto& s){ s.common.solar_wind.C2*=1.001; }},
      {"1-D common apex",[](auto& s){ s.common.apex.radius_m+=1.0; }},
      {"1-D common r0",[](auto& s){ s.common.r0_m+=1.0; }},
      {"1-D region config",[](auto& s){ s.region_config.f_ME+=0.1; }},
      {"1-D region boundaries",[](auto& s){ s.region_boundaries.R_le_m+=1.0; }},
      {"1-D acceleration config",[](auto& s){
         s.acceleration_config.relative_source_weight_per_area+=0.1; }},
      {"1-D kinematics mode",[](auto& s){
         s.kinematics_mode=swcme::kinematics::Mode::DBM; }},
      {"1-D time",[](auto& s){ s.time_s+=1.0; }},
      {"1-D r0 mirror",[](auto& s){ s.r0_m+=1.0; }},
      {"1-D shock radius",[](auto& s){ s.r_sh_m+=1.0; }},
      {"1-D shock speed",[](auto& s){ s.V_sh_ms+=1.0; }},
      {"1-D leading radius",[](auto& s){ s.r_le_m+=1.0; }},
      {"1-D trailing radius",[](auto& s){ s.r_te_m+=1.0; }},
      // SOURCE mode stores +0 here; flipping only its sign bit verifies that
      // record integrity is exact even though configuration fingerprints
      // intentionally canonicalize numerically equivalent signed zeros.
      {"1-D shock-width sign bit",[](auto& s){ s.w_sh_m=-0.0; }},
      {"1-D leading width",[](auto& s){ s.w_le_m+=1.0; }},
      {"1-D trailing width",[](auto& s){ s.w_te_m+=1.0; }},
      {"1-D upstream speed",[](auto& s){ s.V_up_ms+=1.0; }},
      {"1-D Br mirror",[](auto& s){ s.Br1AU_T+=1.0e-12; }},
      {"1-D Parker pitch",[](auto& s){ s.k_AU+=0.1; }},
      {"1-D upstream B",[](auto& s){ s.B_up_T+=1.0e-12; }},
      {"1-D C2 mirror",[](auto& s){ s.C2*=1.001; }},
      {"1-D C4 mirror",[](auto& s){ s.C4*=1.001; }},
      {"1-D C6 mirror",[](auto& s){ s.C6*=1.001; }},
      {"1-D shock flag",[](auto& s){ s.has_shock=!s.has_shock; }},
      {"1-D solver flag",[](auto& s){
         s.shock_solver_converged=!s.shock_solver_converged; }},
      {"1-D RH record",[](auto& s){ s.shock_jump.compression+=0.1; }},
      // SHK12 added public convergence diagnostics to JumpResult.  Corrupt one
      // independently so PST06 proves those non-primitive fields are sealed,
      // rather than relying only on its existing compression mutation.
      {"1-D RH bracket diagnostic",[](auto& s){
         s.shock_jump.root_bracket_width+=1.0e-6; }},
      {"1-D compression mirror",[](auto& s){ s.rc+=0.1; }},
      {"1-D shock density",[](auto& s){ s.n_up_shock+=1.0; }},
      {"1-D leading density",[](auto& s){ s.n_up_le+=1.0; }},
      {"1-D downstream speed",[](auto& s){ s.V2_shock_ms+=1.0; }},
      {"1-D leading speed",[](auto& s){ s.V_LE_ms+=1.0; }}
  };
  for (const auto& mutation : one_mutations)
    reject_1d_corruption(context,one,one_state,mutation);

  // One legacy wrapper confirms the exception bridge preserves the explicit
  // stale-state code instead of evaluating a corrupted record.
  swcme1d::StepState corrupt_legacy=one_state;
  corrupt_legacy.C2*=1.001;
  bool legacy_caught=false;
  try {
    double density=0.0,velocity=0.0;
    one.evaluate_radii_fast(
        corrupt_legacy,&radius,&density,&velocity,1);
  } catch (const std::runtime_error& error) {
    legacy_caught=std::string(error.what()).find("STALE_PREPARED_STATE")!=
                  std::string::npos;
  }
  context.expect_true(legacy_caught,
                      "1-D legacy wrapper preserves stale-state diagnostic");

  swcme3d::Params p3;
  p3.shape=swcme3d::ShockShape::Sphere;
  p3.kinematics_mode=swcme::kinematics::Mode::Ballistic;
  p3.r0_Rs=40.0;
  p3.V0_sh_kms=1200.0;
  p3.V_sw_kms=400.0;
  p3.region_mode=swcme::regions::Mode::ShockOnly;
  p3.shock_acceleration_mode=swcme::acceleration::Mode::Source;
  swcme3d::Model three(p3);
  const swcme3d::StepState three_state=three.prepare_step(3600.0);
  context.expect_true(three_state.integrity_digest()!=0 &&
                          swcme3d::prepared_state_integrity(three_state)==
                              three_state.integrity_digest(),
                      "3-D prepared record has a valid private seal");

  swcme3d::StepState three_copy(three_state);
  swcme3d::StepState three_moved(std::move(three_copy));
  swcme3d::StepState three_copy_assigned;
  three_copy_assigned=three_state;
  swcme3d::StepState three_move_assigned;
  three_move_assigned=std::move(three_copy_assigned);
  const double x=radius,y=0.0,z=0.0;
  for (const auto* state : {&three_moved,&three_move_assigned}) {
    double density=0.0,vx=0.0,vy=0.0,vz=0.0;
    const swcme::ModelStatus status=three.evaluate_cartesian_fast_checked(
        *state,&x,&y,&z,&density,&vx,&vy,&vz,1);
    context.expect_true(status.ok() &&
                            state->owner_model_identity==three.model_identity() &&
                            state->integrity_digest()==three_state.integrity_digest(),
                        "valid 3-D copy/move preserves record integrity");
  }

  const std::vector<Mutation3D> three_mutations={
      {"3-D common solar wind",[](auto& s){ s.common.solar_wind.C2*=1.001; }},
      {"3-D common apex",[](auto& s){ s.common.apex.radius_m+=1.0; }},
      {"3-D common r0",[](auto& s){ s.common.r0_m+=1.0; }},
      {"3-D region config",[](auto& s){ s.region_config.f_ME+=0.1; }},
      {"3-D apex boundaries",[](auto& s){ s.apex_regions.R_le_m+=1.0; }},
      {"3-D acceleration config",[](auto& s){
         s.acceleration_config.relative_source_weight_per_area+=0.1; }},
      {"3-D e1 frame",[](auto& s){ s.e1[0]+=0.1; }},
      {"3-D e2 frame",[](auto& s){ s.e2[1]+=0.1; }},
      {"3-D e3 frame",[](auto& s){ s.e3[2]+=0.1; }},
      {"3-D kinematics mode",[](auto& s){
         s.kinematics_mode=swcme::kinematics::Mode::DBM; }},
      {"3-D time",[](auto& s){ s.time_s+=1.0; }},
      {"3-D shock radius",[](auto& s){ s.r_sh_m+=1.0; }},
      {"3-D shock speed",[](auto& s){ s.V_sh_ms+=1.0; }},
      {"3-D reference size",[](auto& s){ s.a_m+=1.0; }},
      {"3-D sheath thickness",[](auto& s){ s.dr_sheath_m+=1.0; }},
      {"3-D ejecta thickness",[](auto& s){ s.dr_me_m+=1.0; }},
      {"3-D shock width",[](auto& s){ s.w_shock_m+=1.0; }},
      {"3-D leading width",[](auto& s){ s.w_le_m+=1.0; }},
      {"3-D trailing width",[](auto& s){ s.w_te_m+=1.0; }},
      {"3-D leading radius",[](auto& s){ s.r_le_m+=1.0; }},
      {"3-D trailing radius",[](auto& s){ s.r_te_m+=1.0; }},
      {"3-D leading speed",[](auto& s){ s.V_sheath_LE_ms+=1.0; }},
      {"3-D ejecta speed",[](auto& s){ s.V_ME_ms+=1.0; }},
      {"3-D downstream speed",[](auto& s){ s.V_dn_ms+=1.0; }},
      {"3-D shock flag",[](auto& s){ s.has_shock=!s.has_shock; }},
      {"3-D compression",[](auto& s){ s.rc+=0.1; }},
      {"3-D inverse sheath",[](auto& s){ s.inv_dr_sheath+=1.0e-12; }},
      {"3-D compression floor",[](auto& s){ s.rc_floor+=0.1; }},
      {"3-D wind speed",[](auto& s){ s.V_sw_ms+=1.0; }},
      {"3-D C2",[](auto& s){ s.C2*=1.001; }},
      {"3-D C4",[](auto& s){ s.C4*=1.001; }},
      {"3-D C6",[](auto& s){ s.C6*=1.001; }},
      {"3-D inverse shock-width sign bit",[](auto& s){ s.inv2w_sh=-0.0; }},
      {"3-D inverse leading width",[](auto& s){ s.inv2w_le+=1.0e-12; }},
      {"3-D inverse trailing width",[](auto& s){ s.inv2w_te+=1.0e-12; }},
      {"3-D solar axis",[](auto& s){ s.solar_axis_hat[2]+=0.1; }},
      {"3-D rotation rate",[](auto& s){ s.solar_rotation_rate_rad_s+=1.0e-8; }},
      {"3-D Parker pitch",[](auto& s){ s.k_AU+=0.1; }},
      {"3-D Br",[](auto& s){ s.Br1AU_T+=1.0e-12; }},
      {"3-D ellipsoid a",[](auto& s){ s.a_e+=1.0; }},
      {"3-D ellipsoid b",[](auto& s){ s.b_e+=1.0; }},
      {"3-D ellipsoid c",[](auto& s){ s.c_e+=1.0; }},
      {"3-D inverse a",[](auto& s){ s.inv_a2+=1.0e-20; }},
      {"3-D inverse b",[](auto& s){ s.inv_b2+=1.0e-20; }},
      {"3-D inverse c",[](auto& s){ s.inv_c2+=1.0e-20; }},
      {"3-D sine half width",[](auto& s){ s.sin_half_width+=0.1; }},
      {"3-D cosine half width",[](auto& s){ s.cos_half_width+=0.1; }},
      {"3-D SSE center",[](auto& s){ s.sse_center_m+=1.0; }},
      {"3-D SSE radius",[](auto& s){ s.sse_radius_m+=1.0; }}
  };
  for (const auto& mutation : three_mutations)
    reject_3d_corruption(context,three,three_state,mutation);

  // A deterministic summary must expose both seal values without revealing a
  // mechanism that can modify the expected private seal.
  swcme3d::StepState summary_corruption=three_state;
  summary_corruption.r_sh_m+=1.0;
  double density=821.0,vx=822.0,vy=823.0,vz=824.0;
  const swcme::ModelStatus summary_status=
      three.evaluate_cartesian_fast_checked(
          summary_corruption,&x,&y,&z,&density,&vx,&vy,&vz,1);
  const std::string summary=summary_status.summary();
  context.expect_true(
      summary.find("expected_state_integrity=0x")!=std::string::npos &&
      summary.find("computed_state_integrity=0x")!=std::string::npos,
      "stale-state summary prints both integrity digests");
}
