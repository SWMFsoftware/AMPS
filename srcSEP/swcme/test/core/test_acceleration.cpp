#include "test_framework.hpp"

#include <swcme1d.hpp>
#include <swcme3d.hpp>
#include <swcme_acceleration.hpp>
#include <swcme_regions.hpp>
#include <swcme_solarwind.hpp>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <string>

namespace {

void expect_rel(swcme_test::Context& context, const std::string& label,
                double actual, double expected, double rel_tol,
                double floor=1.0e-30) {
  const double scale=std::max({std::abs(actual),std::abs(expected),floor});
  context.expect_near(actual,expected,rel_tol*scale,label);
}

void equivalent_full_icme(swcme1d::Params& p1, swcme3d::Params& p3) {
  p1.V_sw_kms=400.0;
  p1.n1AU_cm3=5.0;
  p1.B1AU_nT=5.0;
  p1.T_K=1.2e5;
  p1.sin_theta=1.0;
  p1.r0_Rs=20.0;
  p1.V0_sh_kms=1400.0;
  p1.Gamma_kmInv=5.0e-8;
  p1.region_mode=swcme::regions::Mode::FullICME;
  p1.shock_acceleration_mode=swcme::acceleration::Mode::ResolvedCompression;
  p1.edge_smooth_shock_AU_at1AU=0.02;
  p1.sheath_thick_AU_at1AU=0.12;
  p1.ejecta_thick_AU_at1AU=0.20;
  p1.edge_smooth_le_AU_at1AU=0.02;
  p1.edge_smooth_te_AU_at1AU=0.03;
  p1.sheath_ramp_power=2.0;
  p1.V_sheath_LE_factor=1.05;
  p1.f_ME=0.5;
  p1.V_ME_factor=0.8;

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
  p3.region_mode=p1.region_mode;
  p3.shock_acceleration_mode=p1.shock_acceleration_mode;
  p3.edge_smooth_shock_AU_at1AU=p1.edge_smooth_shock_AU_at1AU;
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

struct OneD {
  double n=0.0,V=0.0,Br=0.0,Bphi=0.0,Bmag=0.0,divV=0.0;
};

OneD eval1(const swcme1d::Model& m, const swcme1d::StepState& s, double r) {
  OneD a;
  m.evaluate_radii_with_B_div(s,&r,&a.n,&a.V,&a.Br,&a.Bphi,&a.Bmag,&a.divV,1);
  return a;
}

struct ThreeD {
  double n=0.0,vx=0.0,vy=0.0,vz=0.0,bx=0.0,by=0.0,bz=0.0;
};

ThreeD eval3(const swcme3d::Model& m, const swcme3d::StepState& s, double r) {
  ThreeD a;
  const double x=r,y=0.0,z=0.0;
  m.evaluate_cartesian_with_B(s,&x,&y,&z,&a.n,&a.vx,&a.vy,&a.vz,
                              &a.bx,&a.by,&a.bz,1);
  return a;
}

}  // namespace

void test_acc01(swcme_test::Context& context) {
  std::cout << "ACC01 SOURCE uses explicit source and unmodified SHOCK_ONLY flow\n";

  swcme1d::Params p;
  p.region_mode=swcme::regions::Mode::ShockOnly;
  p.shock_acceleration_mode=swcme::acceleration::Mode::Source;
  p.V0_sh_kms=1400.0;
  p.relative_source_weight_per_area=1.0;
  const swcme1d::Model m(p);
  context.expect_true(m.validate().ok(),"SOURCE/SHOCK_ONLY configuration accepted");
  const auto s=m.prepare_step(2.0*3600.0);
  const auto source=m.shock_acceleration_state(s);

  context.expect_true(source.physical_shock,"fixture has physical shock");
  context.expect_true(source.source_enabled,"SOURCE enables prescribed source");
  context.expect_true(!source.resolved_compression_enabled,
                      "SOURCE disables resolved compression");
  const double qref=3.0*source.compression/(source.compression-1.0);
  expect_rel(context,"DSA phase-space slope q",source.dsa_q_phase_space,qref,1.0e-14,1.0);
  expect_rel(context,"constant relative source weight",
             source.relative_source_weight_per_area,1.0,0.0,1.0);

  // Probe on both sides of the mathematical shock. SHOCK_ONLY must remain the
  // exact analytical background: the DSA source is bookkeeping only and no RH
  // velocity jump is available for div(V) to accelerate the same population.
  for (double factor : {0.98,0.999,1.001,1.02}) {
    const double r=factor*s.r_sh_m;
    const OneD a=eval1(m,s,r);
    const double nref=swcme::solarwind::density_m3(s.common.solar_wind,r);
    expect_rel(context,"SOURCE SHOCK_ONLY density",a.n,nref,2.0e-13,1.0);
    expect_rel(context,"SOURCE SHOCK_ONLY velocity",a.V,s.V_up_ms,2.0e-13,1.0);
  }
}

void test_acc02(swcme_test::Context& context) {
  std::cout << "ACC02 incompatible acceleration/region combinations are rejected\n";

  swcme1d::Params source_full;
  source_full.region_mode=swcme::regions::Mode::FullICME;
  source_full.shock_acceleration_mode=swcme::acceleration::Mode::Source;
  context.expect_true(!swcme1d::Model(source_full).validate().ok(),
                      "SOURCE + FULL_ICME rejected");

  swcme3d::Params resolved_shock_only;
  resolved_shock_only.region_mode=swcme::regions::Mode::ShockOnly;
  resolved_shock_only.shock_acceleration_mode=
      swcme::acceleration::Mode::ResolvedCompression;
  context.expect_true(!swcme3d::Model(resolved_shock_only).validate().ok(),
                      "RESOLVED_COMPRESSION + SHOCK_ONLY rejected");

  swcme1d::Params zero_width;
  zero_width.region_mode=swcme::regions::Mode::FullICME;
  zero_width.shock_acceleration_mode=swcme::acceleration::Mode::ResolvedCompression;
  zero_width.edge_smooth_shock_AU_at1AU=0.0;
  context.expect_true(!swcme1d::Model(zero_width).validate().ok(),
                      "resolved compression requires positive shock width");

  swcme1d::Params negative_weight;
  negative_weight.region_mode=swcme::regions::Mode::ShockOnly;
  negative_weight.shock_acceleration_mode=swcme::acceleration::Mode::Source;
  negative_weight.relative_source_weight_per_area=-0.1;
  context.expect_true(!swcme1d::Model(negative_weight).validate().ok(),
                      "negative source weight rejected");
}

void test_acc03(swcme_test::Context& context) {
  std::cout << "ACC03 RESOLVED_COMPRESSION C1 shock smoothing and RH endpoints\n";

  swcme1d::Params p1;
  swcme3d::Params p3;
  equivalent_full_icme(p1,p3);
  const swcme1d::Model m(p1);
  const auto s=m.prepare_step(2.0*3600.0);
  context.expect_true(s.has_shock && s.shock_solver_converged,
                      "resolved-compression fixture has shock");
  const double w=s.region_boundaries.smooth_shock_width_m;
  context.expect_true(w>0.0,"resolved shock width is finite");
  const double h=0.5*w;
  const double ro=s.r_sh_m+h;
  const double ri=s.r_sh_m-h;
  const double rc=s.r_sh_m;

  const OneD outer=eval1(m,s,ro);
  const OneD inner=eval1(m,s,ri);
  const OneD center=eval1(m,s,rc);
  const double nup_outer=swcme::solarwind::density_m3(s.common.solar_wind,ro);
  const double n2=s.shock_jump.downstream.rho_kg_m3/swcme::constants::PROTON_MASS_KG;
  expect_rel(context,"outer endpoint density is upstream",outer.n,nup_outer,2.0e-13,1.0);
  expect_rel(context,"outer endpoint velocity is upstream",outer.V,s.V_up_ms,2.0e-13,1.0);
  expect_rel(context,"inner endpoint density is exact RH",inner.n,n2,2.0e-12,1.0);
  expect_rel(context,"inner endpoint velocity is exact RH",inner.V,
             s.shock_jump.downstream.velocity_m_s[0],2.0e-12,1.0);
  context.expect_true(center.V>std::min(outer.V,inner.V) &&
                      center.V<std::max(outer.V,inner.V),
                      "shock midpoint lies strictly between endpoint velocities");

  // C1 check for the velocity profile.  smoothstep has zero derivative at each
  // endpoint; the inner sheath profile is also constructed with zero slope at
  // its start.  Compare one-sided slopes over a small fraction of the layer.
  const double eps=1.0e-4*w;
  const double Vo_m=eval1(m,s,ro-eps).V;
  const double Vo_p=eval1(m,s,ro+eps).V;
  const double Vi_m=eval1(m,s,ri-eps).V;
  const double Vi_p=eval1(m,s,ri+eps).V;
  const double slope_outer=(Vo_p-Vo_m)/(2.0*eps);
  const double slope_inner=(Vi_p-Vi_m)/(2.0*eps);
  const double characteristic=std::abs(inner.V-outer.V)/w;
  context.expect_true(std::abs(slope_outer)<1.0e-3*characteristic,
                      "outer shock smoothing derivative tends to zero");
  context.expect_true(std::abs(slope_inner)<1.0e-3*characteristic,
                      "inner shock/sheath derivative tends to zero");
}

void test_acc04(swcme_test::Context& context) {
  std::cout << "ACC04 1-D/3-D resolved shock profile identity on spherical +X ray\n";

  swcme1d::Params p1;
  swcme3d::Params p3;
  equivalent_full_icme(p1,p3);
  const swcme1d::Model m1(p1);
  const swcme3d::Model m3(p3);
  const auto s1=m1.prepare_step(2.0*3600.0);
  const auto s3=m3.prepare_step(2.0*3600.0);
  const double w=s1.region_boundaries.smooth_shock_width_m;
  expect_rel(context,"1D/3D shock smoothing width",w,s3.apex_regions.smooth_shock_width_m,
             2.0e-14,1.0);

  for (double q : {0.0,0.25,0.50,0.75,1.0}) {
    // q=0 outer/upstream, q=1 inner/RH.  Both models call the same smoothstep
    // and use the same upstream/downstream states in this exact reduction.
    const double r=s1.r_sh_m+0.5*w-q*w;
    const OneD a=eval1(m1,s1,r);
    const ThreeD b=eval3(m3,s3,r);
    expect_rel(context,"resolved profile density identity",a.n,b.n,3.0e-12,1.0);
    expect_rel(context,"resolved profile radial velocity identity",a.V,b.vx,3.0e-12,1.0);
    expect_rel(context,"resolved profile Br/Bx identity",a.Br,b.bx,3.0e-12,1.0e-15);
    expect_rel(context,"resolved profile Bphi/By identity",a.Bphi,b.by,3.0e-12,1.0e-15);
  }
}

void test_acc05(swcme_test::Context& context) {
  std::cout << "ACC05 resolved mode disables DSA source record\n";

  swcme3d::Params p;
  p.region_mode=swcme::regions::Mode::FullICME;
  p.shock_acceleration_mode=swcme::acceleration::Mode::ResolvedCompression;
  p.V0_sh_kms=1400.0;
  const swcme3d::Model m(p);
  const auto s=m.prepare_step(2.0*3600.0);
  const double u[3]={1.0,0.0,0.0};
  swcme::acceleration::ShockAccelerationState a;
  context.expect_true(m.shock_acceleration_state(s,u,a),"apex source surface exists");
  context.expect_true(a.physical_shock,"fixture has physical shock");
  context.expect_true(!a.source_enabled,"resolved mode disables prescribed source");
  context.expect_true(a.resolved_compression_enabled,
                      "resolved mode enables compression accelerator");
  context.expect_true(!std::isfinite(a.dsa_q_phase_space),
                      "resolved mode does not publish active DSA slope");
  context.expect_true(a.relative_source_weight_per_area==0.0,
                      "resolved mode has zero source weight");
  const std::string text=swcme::acceleration::serialize_csv(a);
  context.expect_true(text.find(",NA,")!=std::string::npos,
                      "serialized resolved record marks DSA slope unavailable");
}
