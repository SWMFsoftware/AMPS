#include "test_framework.hpp"

#include <swcme1d.hpp>
#include <swcme3d.hpp>
#include <swcme_regions.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

namespace {

void expect_rel(swcme_test::Context& context, const std::string& label,
                double actual, double expected, double rel_tol,
                double abs_floor=1.0e-30) {
  const double scale=std::max({std::abs(actual),std::abs(expected),abs_floor});
  context.expect_near(actual,expected,rel_tol*scale,label);
}

// Build a pair of physically identical radial/equatorial models.  A spherical
// 3-D shock on +X is the exact geometric counterpart of the 1-D radial ray.
void equivalent_params(swcme1d::Params& p1, swcme3d::Params& p3,
                       swcme::regions::Mode mode) {
  p1.V_sw_kms=400.0;
  p1.n1AU_cm3=5.0;
  p1.B1AU_nT=5.0;
  p1.T_K=1.2e5;
  p1.sin_theta=1.0;
  p1.r0_Rs=20.0;
  p1.V0_sh_kms=1400.0;
  p1.Gamma_kmInv=5.0e-8;
  p1.region_mode=mode;
  p1.shock_acceleration_mode =
      (mode==swcme::regions::Mode::ShockOnly)
          ? swcme::acceleration::Mode::Source
          : swcme::acceleration::Mode::ResolvedCompression;
  p1.sheath_thick_AU_at1AU=0.10;
  p1.ejecta_thick_AU_at1AU=0.20;
  p1.edge_smooth_le_AU_at1AU=0.02;
  p1.edge_smooth_te_AU_at1AU=0.03;
  p1.sheath_ramp_power=2.0;
  p1.V_sheath_LE_factor=1.05;
  p1.f_ME=0.50;
  p1.V_ME_factor=0.80;

  p3.V_sw_kms=p1.V_sw_kms;
  p3.n1AU_cm3=p1.n1AU_cm3;
  p3.B1AU_nT=p1.B1AU_nT;
  p3.T_K=p1.T_K;
  p3.gamma_ad=p1.gamma_ad;
  p3.sin_theta=p1.sin_theta;
  p3.kinematics_mode=p1.kinematics_mode;
  p3.r0_Rs=p1.r0_Rs;
  p3.V0_sh_kms=p1.V0_sh_kms;
  p3.Gamma_kmInv=p1.Gamma_kmInv;
  p3.region_mode=mode;
  p3.shock_acceleration_mode=p1.shock_acceleration_mode;
  p3.sheath_thick_AU_at1AU=p1.sheath_thick_AU_at1AU;
  p3.ejecta_thick_AU_at1AU=p1.ejecta_thick_AU_at1AU;
  p3.edge_smooth_le_AU_at1AU=p1.edge_smooth_le_AU_at1AU;
  p3.edge_smooth_te_AU_at1AU=p1.edge_smooth_te_AU_at1AU;
  p3.sheath_ramp_power=p1.sheath_ramp_power;
  p3.V_sheath_LE_factor=p1.V_sheath_LE_factor;
  p3.f_ME=p1.f_ME;
  p3.V_ME_factor=p1.V_ME_factor;
  p3.shape=swcme3d::ShockShape::Sphere;
  p3.cme_dir[0]=1.0; p3.cme_dir[1]=0.0; p3.cme_dir[2]=0.0;
  p3.solar_rotation_axis[0]=0.0;
  p3.solar_rotation_axis[1]=0.0;
  p3.solar_rotation_axis[2]=1.0;
}

struct OneDState {
  double n=0.0,V=0.0,Br=0.0,Bphi=0.0,Bmag=0.0,divV=0.0;
};

OneDState eval1d(const swcme1d::Model& m, const swcme1d::StepState& s,
                 double r) {
  OneDState out;
  m.evaluate_radii_with_B_div(s,&r,&out.n,&out.V,&out.Br,&out.Bphi,
                              &out.Bmag,&out.divV,1);
  return out;
}

struct ThreeDState {
  double n=0.0,vx=0.0,vy=0.0,vz=0.0,bx=0.0,by=0.0,bz=0.0;
};

ThreeDState eval3d_x(const swcme3d::Model& m, const swcme3d::StepState& s,
                     double r) {
  ThreeDState out;
  const double x=r,y=0.0,z=0.0;
  m.evaluate_cartesian_with_B(s,&x,&y,&z,&out.n,&out.vx,&out.vy,&out.vz,
                              &out.bx,&out.by,&out.bz,1);
  return out;
}

}  // namespace

void test_reg01(swcme_test::Context& context) {
  std::cout << "REG01 SHOCK_ONLY upstream-field identity\n";

  swcme1d::Params p1;
  swcme3d::Params p3;
  equivalent_params(p1,p3,swcme::regions::Mode::ShockOnly);

  // Deliberately choose extreme FULL_ICME shaping parameters.  SHOCK_ONLY must
  // ignore all of them, proving that phenomenological region configuration
  // cannot leak into the controlled SEP background.
  p1.f_ME=0.15; p1.V_ME_factor=0.25;
  p1.sheath_thick_AU_at1AU=0.25; p1.ejecta_thick_AU_at1AU=0.50;
  p3.f_ME=p1.f_ME; p3.V_ME_factor=p1.V_ME_factor;
  p3.sheath_thick_AU_at1AU=p1.sheath_thick_AU_at1AU;
  p3.ejecta_thick_AU_at1AU=p1.ejecta_thick_AU_at1AU;

  const swcme1d::Model m1(p1);
  const swcme3d::Model m3(p3);
  const auto s1=m1.prepare_step(2.0*3600.0);
  const auto s3=m3.prepare_step(2.0*3600.0);

  const double scales[]={0.45,0.65,0.85,0.99,1.01,1.50};
  for (double q: scales) {
    const double r=q*s1.r_sh_m;
    const OneDState a=eval1d(m1,s1,r);
    const double nref=swcme::solarwind::density_m3(s1.common.solar_wind,r);
    const auto Bref=swcme::solarwind::parker_components(
        s1.common.solar_wind,r,p1.sin_theta);
    expect_rel(context,"1D SHOCK_ONLY density",a.n,nref,2.0e-13,1.0);
    expect_rel(context,"1D SHOCK_ONLY velocity",a.V,s1.V_up_ms,2.0e-13,1.0);
    expect_rel(context,"1D SHOCK_ONLY Br",a.Br,Bref.Br_T,2.0e-13,1.0e-15);
    expect_rel(context,"1D SHOCK_ONLY Bphi",a.Bphi,Bref.Bphi_T,2.0e-13,1.0e-15);

    const ThreeDState b=eval3d_x(m3,s3,r);
    const double u[3]={1.0,0.0,0.0};
    const std::array<double,3> axis={{0.0,0.0,1.0}};
    const std::array<double,3> radial={{1.0,0.0,0.0}};
    const auto B3=swcme::solarwind::parker_field_cartesian(
        s3.common.solar_wind,axis,radial,r);
    expect_rel(context,"3D SHOCK_ONLY density",b.n,nref,2.0e-13,1.0);
    expect_rel(context,"3D SHOCK_ONLY Vx",b.vx,s3.V_sw_ms*u[0],2.0e-13,1.0);
    expect_rel(context,"3D SHOCK_ONLY Vy",b.vy,0.0,2.0e-13,1.0);
    expect_rel(context,"3D SHOCK_ONLY Bx",b.bx,B3[0],2.0e-13,1.0e-15);
    expect_rel(context,"3D SHOCK_ONLY By",b.by,B3[1],2.0e-13,1.0e-15);
  }
}

void test_reg02(swcme_test::Context& context) {
  std::cout << "REG02 FULL_ICME resolved-shock inner RH boundary\n";

  swcme1d::Params p1;
  swcme3d::Params p3;
  equivalent_params(p1,p3,swcme::regions::Mode::FullICME);
  const swcme1d::Model m1(p1);
  const swcme3d::Model m3(p3);
  const auto s1=m1.prepare_step(2.0*3600.0);
  const auto s3=m3.prepare_step(2.0*3600.0);
  context.expect_true(s1.has_shock && s1.shock_solver_converged,
                      "1D fixture has a physical fast shock");

  const double u[3]={1.0,0.0,0.0};
  swcme3d::LocalShockState sh3;
  context.expect_true(m3.shock_state_direction(s3,u,sh3),
                      "3D sphere surface exists");
  context.expect_true(sh3.has_shock && sh3.solver_converged,
                      "3D fixture has a physical fast shock");

  // RESOLVED_COMPRESSION uses a symmetric C1 shock layer.  Its INNER edge is
  // pinned to the exact RH downstream state; the physical discontinuous state
  // remains available independently through shock diagnostics.  This protects
  // the RH boundary while giving the transport solver a finite div(V).
  const double r=s1.r_sh_m-0.5*s1.region_boundaries.smooth_shock_width_m;
  const OneDState a=eval1d(m1,s1,r);
  const double n2=s1.shock_jump.downstream.rho_kg_m3/
                  swcme::constants::PROTON_MASS_KG;
  expect_rel(context,"1D shock-layer inner density equals RH",a.n,n2,2.0e-10,1.0);
  expect_rel(context,"1D shock-layer inner speed equals RH",a.V,
             s1.shock_jump.downstream.velocity_m_s[0],2.0e-10,1.0);
  expect_rel(context,"1D shock-layer inner Br equals RH",a.Br,
             s1.shock_jump.downstream.magnetic_T[0],2.0e-10,1.0e-15);
  expect_rel(context,"1D shock-layer inner Bphi equals RH",a.Bphi,
             s1.shock_jump.downstream.magnetic_T[1],2.0e-10,1.0e-15);

  const ThreeDState b=eval3d_x(m3,s3,r);
  expect_rel(context,"3D shock-layer inner density equals RH",b.n,sh3.downstream_n_m3,
             2.0e-10,1.0);
  const double vv[3]={b.vx,b.vy,b.vz};
  const double bb[3]={b.bx,b.by,b.bz};
  for (int k=0;k<3;++k) {
    expect_rel(context,"3D shock-layer inner velocity component",vv[k],
               sh3.downstream.velocity_m_s[k],2.0e-10,1.0);
    expect_rel(context,"3D shock-layer inner B component",bb[k],
               sh3.downstream.magnetic_T[k],2.0e-10,1.0e-15);
  }
}

void test_reg03(swcme_test::Context& context) {
  std::cout << "REG03 magnetic-ejecta density and velocity factors\n";

  const std::pair<double,double> factors[]={{0.5,0.8},{1.0,1.0},{1.4,1.2}};
  for (const auto& factor: factors) {
    swcme1d::Params p1;
    swcme3d::Params p3;
    equivalent_params(p1,p3,swcme::regions::Mode::FullICME);
    p1.f_ME=factor.first; p1.V_ME_factor=factor.second;
    p3.f_ME=factor.first; p3.V_ME_factor=factor.second;
    const swcme1d::Model m1(p1);
    const swcme3d::Model m3(p3);
    const auto s1=m1.prepare_step(2.0*3600.0);
    const auto s3=m3.prepare_step(2.0*3600.0);

    // Mid-ejecta is well outside the LE/TE transition widths for the fixture,
    // so the configured factors must be honored exactly rather than hidden by
    // an interface blend.
    const double r=0.5*(s1.region_boundaries.R_le_m+s1.region_boundaries.R_te_m);
    const double nup=swcme::solarwind::density_m3(s1.common.solar_wind,r);
    const OneDState a=eval1d(m1,s1,r);
    expect_rel(context,"1D ejecta density factor",a.n,factor.first*nup,2.0e-13,1.0);
    expect_rel(context,"1D ejecta velocity factor",a.V,factor.second*s1.V_up_ms,
               2.0e-13,1.0);

    const ThreeDState b=eval3d_x(m3,s3,r);
    expect_rel(context,"3D ejecta density factor",b.n,factor.first*nup,2.0e-13,1.0);
    expect_rel(context,"3D ejecta velocity factor",b.vx,factor.second*s3.V_sw_ms,
               2.0e-13,1.0);
  }

  // Negative factors are configuration errors, not values to clip inside a hot
  // evaluator.  This permanently guards the repaired no-silent-clipping policy.
  swcme1d::Params invalid;
  invalid.f_ME=-0.1;
  const swcme1d::Model bad(invalid);
  context.expect_true(!bad.validate().ok(),"negative f_ME rejected by validation");
}

void test_reg04(swcme_test::Context& context) {
  std::cout << "REG04 self-similar nested shock/LE/TE surfaces\n";

  swcme3d::Params p;
  p.shape=swcme3d::ShockShape::SSE;
  p.half_width_rad=60.0*swcme::constants::PI/180.0;
  p.cme_dir[0]=1.0; p.cme_dir[1]=0.0; p.cme_dir[2]=0.0;
  p.solar_rotation_axis[0]=0.0; p.solar_rotation_axis[1]=0.0;
  p.solar_rotation_axis[2]=1.0;
  p.region_mode=swcme::regions::Mode::FullICME;
  p.sheath_thick_AU_at1AU=0.12;
  p.ejecta_thick_AU_at1AU=0.28;
  const swcme3d::Model m(p);
  const auto s=m.prepare_step(3600.0);

  const double expected_le=1.0-p.sheath_thick_AU_at1AU;
  const double expected_te=1.0-p.sheath_thick_AU_at1AU-p.ejecta_thick_AU_at1AU;
  const double angles_deg[]={0.0,10.0,25.0,40.0,55.0,59.5};
  for (double deg: angles_deg) {
    const double a=deg*swcme::constants::PI/180.0;
    const double u[3]={std::cos(a),std::sin(a),0.0};
    swcme3d::LocalShockState shock;
    context.expect_true(m.shock_state_direction(s,u,shock),
                        "finite SSE surface exists inside half width");
    const auto b=swcme::regions::make_boundaries(shock.Rdir_m,s.region_config);
    context.expect_true(b.R_sh_m>b.R_le_m && b.R_le_m>b.R_te_m && b.R_te_m>0.0,
                        "local FULL_ICME boundaries remain correctly ordered");
    expect_rel(context,"R_LE/R_sh self-similar",b.R_le_m/b.R_sh_m,
               expected_le,2.0e-14,1.0);
    expect_rel(context,"R_TE/R_sh self-similar",b.R_te_m/b.R_sh_m,
               expected_te,2.0e-14,1.0);
    expect_rel(context,"local sheath fractional thickness",
               b.sheath_thickness_m/b.R_sh_m,p.sheath_thick_AU_at1AU,
               2.0e-14,1.0);
    expect_rel(context,"local ejecta fractional thickness",
               b.ejecta_thickness_m/b.R_sh_m,p.ejecta_thick_AU_at1AU,
               2.0e-14,1.0);
  }

  // An inverted public layer configuration is rejected rather than repaired by
  // sorting/clipping radii at runtime.
  swcme3d::Params invalid=p;
  invalid.sheath_thick_AU_at1AU=0.6;
  invalid.ejecta_thick_AU_at1AU=0.5;
  const swcme3d::Model bad(invalid);
  context.expect_true(!bad.validate().ok(),"inverted layer fractions rejected");
}

void test_reg05(swcme_test::Context& context) {
  std::cout << "REG05 continuity and smoothness at artificial region transitions\n";

  swcme1d::Params p1;
  swcme3d::Params p3;
  equivalent_params(p1,p3,swcme::regions::Mode::FullICME);
  const swcme1d::Model m1(p1);
  const swcme3d::Model m3(p3);
  const auto s1=m1.prepare_step(2.0*3600.0);
  const auto s3=m3.prepare_step(2.0*3600.0);
  const auto& b=s1.region_boundaries;

  // First validate the exact common blend contract itself.  The total widths
  // are centered on LE/TE and smoothstep has zero endpoint derivative, which is
  // the mathematical reason the non-shock artificial transitions are C1.
  for (double center_width_pair : {b.smooth_le_width_m,b.smooth_te_width_m}) {
    context.expect_true(center_width_pair>0.0,"reference transition has nonzero width");
  }
  const double hle=0.5*b.smooth_le_width_m;
  const auto le_outer=swcme::regions::locate(b.R_le_m+hle,b);
  const auto le_inner=swcme::regions::locate(b.R_le_m-hle,b);
  expect_rel(context,"LE outer blend is zero",le_outer.blend,0.0,0.0,1.0);
  expect_rel(context,"LE inner blend is one",le_inner.blend,1.0,0.0,1.0);
  const double hte=0.5*b.smooth_te_width_m;
  const auto te_outer=swcme::regions::locate(b.R_te_m+hte,b);
  const auto te_inner=swcme::regions::locate(b.R_te_m-hte,b);
  expect_rel(context,"TE outer blend is zero",te_outer.blend,0.0,0.0,1.0);
  expect_rel(context,"TE inner blend is one",te_inner.blend,1.0,0.0,1.0);

  // Exercise both public dimensional evaluators at and around every transition
  // endpoint.  Equivalent sphere/+X configurations should match to roundoff,
  // while the small two-sided jumps shrink linearly with h because there is no
  // hidden discontinuity.  The physical shock itself is intentionally excluded.
  const double points[]={b.R_le_m+hle,b.R_le_m-hle,
                         b.R_te_m+hte,b.R_te_m-hte};
  for (double r0: points) {
    const double h=1.0e-6*std::max(b.smooth_le_width_m,b.smooth_te_width_m);
    const OneDState am=eval1d(m1,s1,r0-h);
    const OneDState a0=eval1d(m1,s1,r0);
    const OneDState ap=eval1d(m1,s1,r0+h);
    const ThreeDState c0=eval3d_x(m3,s3,r0);
    expect_rel(context,"REG05 1D/3D density identity",a0.n,c0.n,2.0e-12,1.0);
    expect_rel(context,"REG05 1D/3D radial V identity",a0.V,c0.vx,2.0e-12,1.0);
    expect_rel(context,"REG05 1D/3D Br/Bx identity",a0.Br,c0.bx,2.0e-12,1.0e-15);
    expect_rel(context,"REG05 1D/3D Bphi/By identity",a0.Bphi,c0.by,2.0e-12,1.0e-15);

    // One-sided derivatives should converge toward one another at a C1 edge.
    const double dn_left=(a0.n-am.n)/h;
    const double dn_right=(ap.n-a0.n)/h;
    const double dv_left=(a0.V-am.V)/h;
    const double dv_right=(ap.V-a0.V)/h;
    const double nscale=std::max({std::abs(dn_left),std::abs(dn_right),1.0e-30});
    const double vscale=std::max({std::abs(dv_left),std::abs(dv_right),1.0e-30});
    context.expect_true(std::abs(dn_left-dn_right)/nscale<2.0e-3,
                        "density one-sided derivatives agree at C1 edge");

    // Near a constant-speed segment the true velocity derivative is exactly
    // zero, so a purely relative derivative comparison is ill-conditioned: one
    // side can round to zero while the other contains a ~1e-9 1/s finite-
    // difference remainder.  Combine the convergence-relative criterion with
    // an absolute scale based on V_sw divided by the narrowest transition.
    const double transition_scale=std::max(1.0,
        std::min(b.smooth_le_width_m,b.smooth_te_width_m));
    const double dv_abs_tol=5.0e-6*s1.V_up_ms/transition_scale;
    context.expect_true(std::abs(dv_left-dv_right)<2.0e-3*vscale+dv_abs_tol,
                        "velocity one-sided derivatives agree at C1 edge");
  }
}
