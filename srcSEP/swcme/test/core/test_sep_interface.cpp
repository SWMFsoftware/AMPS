#include "test_framework.hpp"

#include <swcme_sep_interface.hpp>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <numeric>
#include <string>

namespace {

void expect_rel(swcme_test::Context& context,
                const std::string& label,
                double actual,
                double expected,
                double rel_tol=2.0e-12,
                double floor=1.0e-30) {
  const double scale=std::max({std::abs(actual),std::abs(expected),floor});
  context.expect_near(actual,expected,rel_tol*scale,label);
}

void make_equivalent_source_params(swcme1d::Params& p1, swcme3d::Params& p3) {
  p1.V_sw_kms=410.0;
  p1.n1AU_cm3=5.5;
  p1.B1AU_nT=5.2;
  p1.T_K=1.25e5;
  p1.gamma_ad=5.0/3.0;
  p1.sin_theta=1.0;
  p1.kinematics_mode=swcme::kinematics::Mode::DBM;
  p1.r0_Rs=20.0;
  p1.V0_sh_kms=1450.0;
  p1.Gamma_kmInv=5.0e-8;
  p1.region_mode=swcme::regions::Mode::ShockOnly;
  p1.shock_acceleration_mode=swcme::acceleration::Mode::Source;
  p1.relative_source_weight_per_area=1.0;

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
  p3.relative_source_weight_per_area=p1.relative_source_weight_per_area;
  p3.shape=swcme3d::ShockShape::Sphere;
  p3.cme_dir[0]=1.0; p3.cme_dir[1]=0.0; p3.cme_dir[2]=0.0;
  p3.solar_rotation_axis[0]=0.0;
  p3.solar_rotation_axis[1]=0.0;
  p3.solar_rotation_axis[2]=1.0;
}

}  // namespace

void test_sep01(swcme_test::Context& context) {
  std::cout << "SEP01 AMPS background adapter reproduces direct production queries\n";

  swcme3d::Params p;
  p.shape=swcme3d::ShockShape::Sphere;
  p.region_mode=swcme::regions::Mode::ShockOnly;
  p.shock_acceleration_mode=swcme::acceleration::Mode::Source;
  swcme::sep::Interface3D adapter(p);
  const auto step=adapter.prepare(2.0*3600.0);

  const std::array<double,3> x{{0.8*swcme::constants::AU_M,
                                0.17*swcme::constants::AU_M,
                               -0.08*swcme::constants::AU_M}};
  swcme::sep::BackgroundState wrapped;
  const swcme::ModelStatus status=adapter.evaluate_background(step,x,wrapped);
  context.expect_true(status.ok(),"3-D SEP background query succeeds");

  const double px=x[0],py=x[1],pz=x[2];
  double n=0,vx=0,vy=0,vz=0,bx=0,by=0,bz=0,divv=0;
  const swcme::ModelStatus direct=adapter.model().evaluate_cartesian_with_B_div_checked(
      step,&px,&py,&pz,&n,&vx,&vy,&vz,&bx,&by,&bz,&divv,1);
  context.expect_true(direct.ok(),"direct 3-D query succeeds");
  expect_rel(context,"adapter density",wrapped.density_m3,n);
  expect_rel(context,"adapter Vx",wrapped.velocity_m_s[0],vx);
  expect_rel(context,"adapter Vy",wrapped.velocity_m_s[1],vy);
  expect_rel(context,"adapter Vz",wrapped.velocity_m_s[2],vz);
  expect_rel(context,"adapter Bx",wrapped.magnetic_T[0],bx);
  expect_rel(context,"adapter By",wrapped.magnetic_T[1],by);
  expect_rel(context,"adapter Bz",wrapped.magnetic_T[2],bz);
  expect_rel(context,"adapter divV",wrapped.div_velocity_s_inv,divv);
  expect_rel(context,"adapter |B|",wrapped.magnetic_magnitude_T,
             std::hypot(bx,std::hypot(by,bz)));
}

void test_sep02(swcme_test::Context& context) {
  std::cout << "SEP02 source spectrum units and DSA momentum/intensity conversion\n";

  // Unit conversion is deliberately checked with a value that makes the large
  // SI scale obvious.  Round-trip must be exact to floating-point precision.
  const double j_common=17.25; // (cm^2 s sr MeV)^-1
  const double j_si=swcme::sep::differential_intensity_common_to_SI(j_common);
  expect_rel(context,"intensity unit round trip",
             swcme::sep::differential_intensity_SI_to_common(j_si),j_common,2e-15,1.0);

  // A 1-GeV proton has pc=sqrt(K(K+2mc^2)) ~1.696 GeV, so its rigidity is
  // about 1.696 GV for Z=1.  This is an independent numerical sanity check on
  // the relativistic energy/momentum conversion used by the source adapter.
  const double rigidity=swcme::sep::rigidity_GV_from_kinetic_MeV(
      1000.0,swcme::constants::PROTON_MASS_KG,1);
  context.expect_near(rigidity,1.696,3.0e-3,"1-GeV proton rigidity benchmark");

  swcme1d::Params p;
  p.V0_sh_kms=1450.0;
  p.region_mode=swcme::regions::Mode::ShockOnly;
  p.shock_acceleration_mode=swcme::acceleration::Mode::Source;
  swcme::sep::SpectrumConfig spectrum;
  spectrum.kinetic_energy_min_MeV=1.0;
  spectrum.kinetic_energy_max_MeV=1000.0;
  spectrum.reference_energy_MeV=10.0;
  spectrum.normalization=swcme::sep::NormalizationMode::ReferenceDifferentialIntensity;
  spectrum.reference_differential_intensity_SI=j_si;
  swcme::sep::Interface1D adapter(p,spectrum);
  const auto step=adapter.prepare(2.0*3600.0);
  swcme::sep::SEPSourceState source;
  context.expect_true(adapter.source_at_shock(step,source).ok(),"1-D SEP source builds");
  context.expect_true(source.active,"fixture exposes active SOURCE shock");
  expect_rel(context,"momentum intensity index",source.momentum_intensity_index,
             source.q_phase_space-2.0,2e-15,1.0);
  expect_rel(context,"nonrel energy index",source.nonrel_energy_intensity_index,
             0.5*(source.q_phase_space-2.0),2e-15,1.0);

  double shape_ref=0.0;
  context.expect_true(swcme::sep::relative_intensity_shape(
      source,spectrum.reference_energy_MeV,shape_ref).ok(),"shape at Eref evaluates");
  expect_rel(context,"shape(Eref)=1",shape_ref,1.0,2e-15,1.0);
  double physical_ref=0.0;
  context.expect_true(swcme::sep::differential_intensity_SI(
      source,spectrum.reference_energy_MeV,physical_ref).ok(),"J(Eref) evaluates");
  expect_rel(context,"J(Eref)=configured normalization",physical_ref,j_si,2e-15,1.0);

  double shape_hi=0.0;
  context.expect_true(swcme::sep::relative_intensity_shape(source,100.0,shape_hi).ok(),
                      "high-energy shape evaluates");
  context.expect_true(shape_hi>0.0 && shape_hi<1.0,
                      "compressive DSA spectrum decreases above reference energy");

  // The deterministic handoff record must carry its status and all spectrum
  // normalization inputs; otherwise a CSV can no longer be interpreted
  // independently of the process that wrote it.
  const std::string header=swcme::sep::source_csv_header();
  context.expect_true(header.find("status,")==0,
                      "serialized source header begins with explicit status");
  context.expect_true(header.find("particle_mass_kg")!=std::string::npos,
                      "serialized source header includes particle mass");
  context.expect_true(header.find("charge_number")!=std::string::npos,
                      "serialized source header includes charge number");
  context.expect_true(header.find("relative_source_weight_per_area")!=std::string::npos,
                      "serialized source header includes raw source weighting input");
  const std::string record=swcme::sep::serialize_source_csv(source);
  context.expect_true(record.rfind("OK,",0)==0,
                      "serialized active source begins with explicit OK status");
}

void test_sep03(swcme_test::Context& context) {
  std::cout << "SEP03 1-D/3-D adapters emit identical field-aligned source records\n";

  swcme1d::Params p1;
  swcme3d::Params p3;
  make_equivalent_source_params(p1,p3);
  swcme::sep::SpectrumConfig spectrum;
  spectrum.kinetic_energy_min_MeV=2.0;
  spectrum.kinetic_energy_max_MeV=500.0;
  spectrum.reference_energy_MeV=20.0;

  swcme::sep::Interface1D a1(p1,spectrum);
  swcme::sep::Interface3D a3(p3,spectrum);
  const auto s1=a1.prepare(2.0*3600.0);
  const auto s3=a3.prepare(2.0*3600.0);
  swcme::sep::SEPSourceState q1,q3;
  context.expect_true(a1.source_at_shock(s1,q1).ok(),"1-D source record available");
  context.expect_true(a3.source_at_direction(s3,{{1.0,0.0,0.0}},q3).ok(),
                      "3-D +X source record available");
  context.expect_true(q1.active && q3.active,"both equivalent sources active");
  context.expect_true(swcme::sep::serialize_source_csv(q1)==
                      swcme::sep::serialize_source_csv(q3),
                      "AMPS-facing serialized source records are byte-identical");

  // Background identity through the adapters is a direct integration guard:
  // the AMPS path must not alter units or signs relative to standalone SWCME.
  swcme::sep::BackgroundState b1,b3;
  const double r=1.0*swcme::constants::AU_M;
  context.expect_true(a1.evaluate_background(s1,r,b1).ok(),"1-D background adapter");
  context.expect_true(a3.evaluate_background(s3,{{r,0.0,0.0}},b3).ok(),"3-D background adapter");
  expect_rel(context,"adapter 1D/3D density",b1.density_m3,b3.density_m3);
  expect_rel(context,"adapter 1D/3D Vx",b1.velocity_m_s[0],b3.velocity_m_s[0]);
  expect_rel(context,"adapter 1D/3D Bx",b1.magnetic_T[0],b3.magnetic_T[0]);
  expect_rel(context,"adapter 1D/3D By",b1.magnetic_T[1],b3.magnetic_T[1]);
  expect_rel(context,"adapter 1D/3D divV",b1.div_velocity_s_inv,b3.div_velocity_s_inv);
}

void test_sep04(swcme_test::Context& context) {
  std::cout << "SEP04 shock-surface source patches use physical-area normalization\n";

  swcme3d::Params p;
  p.shape=swcme3d::ShockShape::SSE;
  p.half_width_rad=40.0*swcme::constants::PI/180.0;
  p.V0_sh_kms=1450.0;
  p.region_mode=swcme::regions::Mode::ShockOnly;
  p.shock_acceleration_mode=swcme::acceleration::Mode::Source;
  p.relative_source_weight_per_area=2.5;
  swcme::sep::Interface3D adapter(p);
  const auto step=adapter.prepare(2.0*3600.0);
  swcme::sep::SourceSurface surface;
  const swcme::ModelStatus status=adapter.build_shock_surface_source(step,8,16,surface);
  context.expect_true(status.ok(),"source surface builds");
  context.expect_true(surface.patch_count==surface.patches.size() && surface.patch_count>0,
                      "one source record per triangle");
  context.expect_true(surface.active_patch_count>0,"physical source cells exist");
  context.expect_true(surface.active_surface_area_m2>0.0,"active area positive");
  context.expect_true(surface.total_surface_area_m2>=surface.active_surface_area_m2,
                      "active area cannot exceed mesh area");

  double area_fraction_sum=0.0;
  double relative_weight_sum=0.0;
  std::size_t active_count=0;
  for (std::size_t i=0;i<surface.patches.size();++i) {
    const auto& source=surface.patches[i];
    context.expect_true(source.source_id==i,"deterministic source patch id");
    context.expect_true(source.patch_area_m2>0.0,"source patch physical area positive");
    if (source.active) {
      ++active_count;
      area_fraction_sum+=source.area_fraction;
      relative_weight_sum+=source.relative_patch_weight;
      context.expect_true(source.area_fraction>0.0,"active patch area fraction positive");
    } else {
      context.expect_near(source.area_fraction,0.0,0.0,"inactive patch area fraction zero");
      context.expect_near(source.relative_patch_weight,0.0,0.0,"inactive patch weight zero");
    }
  }
  context.expect_true(active_count==surface.active_patch_count,"active count metadata exact");
  expect_rel(context,"active area fractions sum to one",area_fraction_sum,1.0,2e-13,1.0);
  expect_rel(context,"relative patch weights preserve configured source weight",
             relative_weight_sum,p.relative_source_weight_per_area,2e-13,1.0);
}

void test_sep05(swcme_test::Context& context) {
  std::cout << "SEP05 observer cobpoint source reuses production connectivity/shock state\n";

  swcme3d::Params p;
  p.shape=swcme3d::ShockShape::Sphere;
  p.V0_sh_kms=1450.0;
  p.region_mode=swcme::regions::Mode::ShockOnly;
  p.shock_acceleration_mode=swcme::acceleration::Mode::Source;
  swcme::sep::Interface3D adapter(p);
  const auto step=adapter.prepare(8.0*3600.0);
  const std::array<double,3> observer{{swcme::constants::AU_M,0.0,0.0}};
  swcme::sep::SEPSourceState source;
  swcme3d::ConnectivityState connectivity;
  const swcme::ModelStatus status=adapter.source_at_observer_cobpoint(
      step,observer,source,&connectivity);
  context.expect_true(status.ok(),"observer has Parker/shock connection");
  context.expect_true(connectivity.connected && source.connection_evaluated && source.connected,
                      "connection flags propagate to source record");
  context.expect_true(source.active,"connected physical shock activates source");
  const auto& root=connectivity.roots[connectivity.selected_root];
  // The connectivity root is accepted at the documented sub-1e-8-AU
  // residual, while source_at_direction reprojects that direction onto the
  // exact analytical surface.  Their Cartesian positions therefore need only
  // agree to the root tolerance, not bit-for-bit.
  const double root_position_tolerance=1.0e-8*swcme::constants::AU_M;
  for (int k=0;k<3;++k) {
    context.expect_near(source.position_m[k],root.position_m[k],
                        root_position_tolerance,
                        "cobpoint/source position component");
  }
  expect_rel(context,"cobpoint/source compression",source.compression,
             root.shock.compression,2e-12,1.0);
  expect_rel(context,"cobpoint/source thetaBn",source.theta_Bn_rad,
             root.shock.theta_Bn_rad,2e-12,1.0);
  expect_rel(context,"cobpoint/source fast Mach",source.fast_mach,
             root.shock.fast_mach,2e-12,1.0);
}

void test_sep06(swcme_test::Context& context) {
  std::cout << "SEP06 resolved-compression mode exposes no prescribed SEP source\n";

  swcme3d::Params p;
  p.shape=swcme3d::ShockShape::Sphere;
  p.V0_sh_kms=1450.0;
  p.region_mode=swcme::regions::Mode::FullICME;
  p.shock_acceleration_mode=swcme::acceleration::Mode::ResolvedCompression;
  p.edge_smooth_shock_AU_at1AU=0.01;
  swcme::sep::Interface3D adapter(p);
  const auto step=adapter.prepare(2.0*3600.0);
  swcme::sep::SEPSourceState source;
  context.expect_true(adapter.source_at_direction(step,{{1.0,0.0,0.0}},source).ok(),
                      "resolved-compression diagnostic converts cleanly");
  context.expect_true(!source.active,"resolved-compression source is inactive");
  context.expect_true(!std::isfinite(source.q_phase_space),"inactive source has no DSA q");
  double shape=0.0;
  const swcme::ModelStatus shape_status=
      swcme::sep::relative_intensity_shape(source,10.0,shape);
  context.expect_true(shape_status.code==swcme::StatusCode::SourceInactive,
                      "spectrum request explicitly reports SOURCE_INACTIVE");

  const std::string manifest1=adapter.resolved_manifest();
  const std::string manifest2=adapter.resolved_manifest();
  context.expect_true(manifest1==manifest2,"resolved SEP/model manifest deterministic");
  context.expect_true(manifest1.find("sep_source_contract_version=2")!=std::string::npos,
                      "manifest records SEP source contract version");
  context.expect_true(manifest1.find("shock_acceleration_mode=RESOLVED_COMPRESSION")!=
                      std::string::npos,
                      "manifest records no-source acceleration mode");
}
