#include "test_framework.hpp"

#include "../../swcme1d.hpp"
#include "../../swcme3d.hpp"
#include "../../swcme_shock.hpp"
#include "../../swcme_status.hpp"

#include <array>
#include <cstdio>
#include <limits>
#include <stdexcept>
#include <string>

namespace {

swcme1d::Params params1d() {
  swcme1d::Params p;
  p.r0_Rs=20.0;
  p.V0_sh_kms=1200.0;
  p.V_sw_kms=400.0;
  p.region_mode=swcme::regions::Mode::ShockOnly;
  p.shock_acceleration_mode=swcme::acceleration::Mode::Source;
  return p;
}

swcme3d::Params params3d() {
  swcme3d::Params p;
  p.shape=swcme3d::ShockShape::Sphere;
  p.r0_Rs=20.0;
  p.V0_sh_kms=1200.0;
  p.V_sw_kms=400.0;
  p.region_mode=swcme::regions::Mode::ShockOnly;
  p.shock_acceleration_mode=swcme::acceleration::Mode::Source;
  return p;
}

}  // namespace

// ERR01: a query outside the analytical solar-wind domain must be visible to
// the caller.  The old implementation silently clipped it to 1.05 R_sun and
// returned a plausible density/field for a radius that was never evaluated.
void test_err01(swcme_test::Context& context) {
  swcme1d::Model model(params1d());
  const auto step=model.prepare_step(0.0);
  const double r=0.9*swcme::constants::SOLAR_RADIUS_M;
  double n=123.0,V=456.0;
  const swcme::ModelStatus status=model.evaluate_radii_fast_checked(
      step,&r,&n,&V,1);
  context.expect_true(status.code==swcme::StatusCode::OutsideModelDomain,
                      "1-D below-domain query reports OUTSIDE_MODEL_DOMAIN");
  context.expect_true(status.sample_index==0,
                      "1-D status records failed sample index");
  context.expect_true(n==123.0 && V==456.0,
                      "failed sample is not overwritten by fallback values");

  bool threw=false;
  try { model.evaluate_radii_fast(step,&r,&n,&V,1); }
  catch (const std::runtime_error&) { threw=true; }
  context.expect_true(threw,
                      "legacy 1-D void evaluator throws instead of clipping");
}

// ERR02: Cartesian input failures must not be converted to the former r~0
// pseudo-point/zero-output state.  A checked call returns the exact sample and
// the source-compatible wrapper throws.
void test_err02(swcme_test::Context& context) {
  swcme3d::Model model(params3d());
  const auto step=model.prepare_step(0.0);
  const double x[2]={swcme::constants::AU_M,
                     std::numeric_limits<double>::quiet_NaN()};
  const double y[2]={0.0,0.0},z[2]={0.0,0.0};
  double n[2]={-1.0,-2.0},vx[2]={-1.0,-2.0},vy[2]={-1.0,-2.0},vz[2]={-1.0,-2.0};
  const swcme::ModelStatus status=model.evaluate_cartesian_fast_checked(
      step,x,y,z,n,vx,vy,vz,2);
  context.expect_true(status.code==swcme::StatusCode::NonFiniteInput,
                      "3-D NaN coordinate reports NONFINITE_INPUT");
  context.expect_true(status.sample_index==1,
                      "3-D status records the second sample");
  context.expect_true(n[1]==-2.0 && vx[1]==-2.0,
                      "failed Cartesian sample is not replaced by zeros");

  bool threw=false;
  try { model.evaluate_cartesian_fast(step,x,y,z,n,vx,vy,vz,2); }
  catch (const std::runtime_error&) { threw=true; }
  context.expect_true(threw,
                      "legacy 3-D evaluator throws on checked failure");
}

// ERR03: shock solver outcome is explicit.  INVALID_INPUT, NO_SHOCK,
// NUMERICALLY_UNRESOLVED_WEAK_SHOCK, and SOLVED are distinct states rather
// than all being inferred from compression=1 or one boolean convergence flag.
void test_err03(swcme_test::Context& context) {
  using namespace swcme::shock;
  PrimitiveState up;
  up.rho_kg_m3=5.0e6*swcme::constants::PROTON_MASS_KG;
  up.pressure_Pa=5.0e6*swcme::constants::BOLTZMANN_J_K*1.0e5;
  up.velocity_m_s={{4.0e5,0.0,0.0}};
  up.magnetic_T={{5.0e-9,0.0,0.0}};

  const JumpResult invalid=solve_ideal_mhd_fast_shock(
      up,{{0.0,0.0,0.0}},1.0e6,5.0/3.0);
  context.expect_true(invalid.status==SolveStatus::InvalidInput &&
                      !invalid.solver_converged,
                      "degenerate shock normal reports INVALID_INPUT");

  // Non-finite primitive components used to be especially dangerous because
  // sqrt(max(0,dot(B,B))) could collapse a NaN magnetic norm to zero and let
  // the solver continue as if the field were physically absent.
  PrimitiveState bad_field=up;
  bad_field.magnetic_T[1]=std::numeric_limits<double>::quiet_NaN();
  const JumpResult invalid_field=solve_ideal_mhd_fast_shock(
      bad_field,{{1.0,0.0,0.0}},1.0e6,5.0/3.0);
  context.expect_true(invalid_field.status==SolveStatus::InvalidInput &&
                      !invalid_field.solver_converged,
                      "non-finite upstream magnetic field reports INVALID_INPUT");

  const JumpResult no_shock=solve_ideal_mhd_fast_shock(
      up,{{1.0,0.0,0.0}},4.0e5,5.0/3.0);
  context.expect_true(no_shock.status==SolveStatus::NoShock &&
                      !no_shock.has_shock && no_shock.solver_converged,
                      "sub-fast state reports NO_SHOCK explicitly");

  // A representably super-fast state inside the published binary64 weak-shock
  // resolution is physical (has_shock=true) but has no trustworthy RH jump.
  // It must not be confused with the physical sub-fast case above.
  const JumpResult threshold_probe=solve_ideal_mhd_fast_shock(
      up,{{1.0,0.0,0.0}},8.0e5,5.0/3.0);
  const double unresolved_speed=up.velocity_m_s[0]+threshold_probe.fast_speed_m_s*
      (1.0+0.5*WEAK_SHOCK_MACH_RESOLUTION);
  const JumpResult unresolved=solve_ideal_mhd_fast_shock(
      up,{{1.0,0.0,0.0}},unresolved_speed,5.0/3.0);
  context.expect_true(
      unresolved.status==SolveStatus::NumericallyUnresolvedWeakShock &&
          unresolved.has_shock && !unresolved.solver_converged,
      "sub-resolution super-fast state reports explicit unresolved status");

  const JumpResult solved=solve_ideal_mhd_fast_shock(
      up,{{1.0,0.0,0.0}},1.2e6,5.0/3.0);
  context.expect_true(solved.status==SolveStatus::Solved && solved.has_shock &&
                      solved.solver_converged,
                      "well-conditioned fast shock reports SOLVED");
}

// ERR04: bad direction vectors are not normalized to +X.  The checked shock
// API returns DEGENERATE_VECTOR and the legacy bool API throws instead of
// looking like a legitimate +X surface query.
void test_err04(swcme_test::Context& context) {
  swcme3d::Model model(params3d());
  const auto step=model.prepare_step(0.0);
  const double zero[3]={0.0,0.0,0.0};
  swcme3d::LocalShockState state;
  const swcme::ModelStatus status=model.shock_state_direction_checked(
      step,zero,state);
  context.expect_true(status.code==swcme::StatusCode::DegenerateVector,
                      "zero shock direction reports DEGENERATE_VECTOR");
  context.expect_true(!state.surface_exists,
                      "invalid direction does not manufacture a +X surface");

  bool threw=false;
  try { (void)model.shock_state_direction(step,zero,state); }
  catch (const std::runtime_error&) { threw=true; }
  context.expect_true(threw,
                      "legacy shock-direction query throws on bad vector");
}

// ERR05: output code cannot hide a corrupt physics value by replacing NaN with
// zero.  The checked writer rejects the dataset before file creation and
// returns the originating numerical status.
void test_err05(swcme_test::Context& context) {
  swcme3d::Model model(params3d());
  swcme3d::ShockMesh mesh;
  mesh.x={std::numeric_limits<double>::quiet_NaN(),1.0,0.0};
  mesh.y={0.0,0.0,1.0}; mesh.z={0.0,0.0,0.0};
  mesh.n_hat_x={1.0,1.0,1.0}; mesh.n_hat_y={0.0,0.0,0.0};
  mesh.n_hat_z={0.0,0.0,0.0}; mesh.rc={1.0,1.0,1.0};
  mesh.Vsh_n={0.0,0.0,0.0};
  mesh.tri_i={1}; mesh.tri_j={2}; mesh.tri_k={3};
  swcme3d::TriMetrics metrics;
  const char* path="output/ERR05_should_not_exist.dat";
  std::remove(path);
  const swcme::ModelStatus status=
      model.write_shock_surface_center_metrics_tecplot_checked(mesh,metrics,path);
  context.expect_true(status.code==swcme::StatusCode::NonFiniteResult,
                      "writer reports NONFINITE_RESULT for NaN mesh data");
  std::FILE* f=std::fopen(path,"r");
  context.expect_true(f==nullptr,
                      "writer rejects corrupt dataset before creating output");
  if (f) std::fclose(f);
  std::remove(path);
}
