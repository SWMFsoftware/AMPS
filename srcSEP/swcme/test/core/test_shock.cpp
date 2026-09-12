#include "test_framework.hpp"
#include "../reference/shk05_oblique_v1.hpp"
#include "../reference/shk12_near_mach_v1.hpp"

#include <swcme1d.hpp>
#include <swcme3d.hpp>
#include <swcme_shock.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

namespace {
using Vec3=swcme::shock::Vec3;
using Primitive=swcme::shock::PrimitiveState;

constexpr double GAMMA=5.0/3.0;
constexpr double MP=swcme::constants::PROTON_MASS_KG;
constexpr double KB=swcme::constants::BOLTZMANN_J_K;
constexpr double MU0=swcme::constants::VACUUM_PERMEABILITY_N_A2;

inline double dot(const Vec3& a,const Vec3& b){return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];}
inline Vec3 sub(const Vec3& a,const Vec3& b){return {{a[0]-b[0],a[1]-b[1],a[2]-b[2]}};}
inline Vec3 scale(const Vec3& a,double s){return {{s*a[0],s*a[1],s*a[2]}};}
inline double norm(const Vec3& a){return std::sqrt(std::max(0.0,dot(a,a)));}

void print_check(const std::string& label,double actual,double expected,double tol,bool pass){
  std::cout<<"  "<<std::left<<std::setw(38)<<label<<std::right
           <<" actual="<<std::scientific<<std::setprecision(12)<<actual
           <<" expected="<<expected<<" tol="<<tol<<' '
           <<(pass?"PASS":"FAIL")<<'\n';
}

void expect_rel(swcme_test::Context& c,const std::string& label,
                double actual,double expected,double rel_tol,double floor=1e-30){
  const double scalev=std::max({std::abs(actual),std::abs(expected),floor});
  const double tol=rel_tol*scalev;
  const bool pass=std::isfinite(actual)&&std::isfinite(expected)&&std::abs(actual-expected)<=tol;
  print_check(label,actual,expected,tol,pass); c.expect_true(pass,label);
}

void expect_abs(swcme_test::Context& c,const std::string& label,
                double actual,double expected,double tol){
  const bool pass=std::isfinite(actual)&&std::isfinite(expected)&&std::abs(actual-expected)<=tol;
  print_check(label,actual,expected,tol,pass); c.expect_true(pass,label);
}

Primitive fixture(double theta_deg=40.0,double n_cm3=5.0,double B_nT=5.0,double T=1.2e5){
  Primitive s;
  const double th=theta_deg*swcme::constants::PI/180.0;
  const double n=n_cm3*1e6;
  s.rho_kg_m3=n*MP;
  s.pressure_Pa=n*KB*T;
  s.velocity_m_s={{4.0e5,0.0,0.0}};
  s.magnetic_T={{B_nT*1e-9*std::cos(th),B_nT*1e-9*std::sin(th),0.0}};
  return s;
}

// Independent fast-mode speed used only to construct threshold fixtures.  It
// intentionally does not call the production helper under test.
double reference_fast_speed(const Primitive& s,const Vec3& n,double gamma=GAMMA){
  const double B2=dot(s.magnetic_T,s.magnetic_T);
  const double vA2=B2/(MU0*s.rho_kg_m3);
  const double cs2=gamma*s.pressure_Pa/s.rho_kg_m3;
  const double cosBn=(B2>0.0)?dot(s.magnetic_T,n)/std::sqrt(B2):0.0;
  const double a=vA2+cs2;
  const double disc=std::max(0.0,a*a-4.0*vA2*cs2*cosBn*cosBn);
  return std::sqrt(0.5*(a+std::sqrt(disc)));
}

// Independent perpendicular-shock reference.  For B_n=0 and zero upstream
// tangential flow, ideal-MHD induction gives B2_t=r B1_t and mass continuity
// gives u2_n=u1_n/r.  Normal momentum then determines p2; the scalar energy
// equation below is solved by bisection.  Keeping this test-specific reduction
// independent of swcme_shock.hpp makes SHK04 a real validation rather than a
// second call to the production solver.
double perpendicular_reference_compression(const Primitive& up,double Vsh){
  const double u=Vsh-up.velocity_m_s[0]; // positive inflow magnitude
  const double B=std::abs(up.magnetic_T[1]);
  const double rmax=(GAMMA+1.0)/(GAMMA-1.0);
  auto q=[&](double r){
    const double p2=up.rho_kg_m3*u*u+up.pressure_Pa+B*B/(2*MU0)
                   -up.rho_kg_m3*u*u/r-r*r*B*B/(2*MU0);
    const double E1=0.5*up.rho_kg_m3*u*u
                   +GAMMA/(GAMMA-1.0)*up.pressure_Pa+B*B/MU0;
    const double E2=(1.0/r)*(0.5*up.rho_kg_m3*u*u/r
                   +GAMMA/(GAMMA-1.0)*p2+r*r*B*B/MU0);
    return (E2-E1)/(r-1.0);
  };
  double lo=1.0+1e-8, flo=q(lo), hi=rmax-1e-8, fhi=q(hi);
  if ((flo<0)==(fhi<0)) return 0.0;
  for(int i=0;i<160;++i){
    const double mid=0.5*(lo+hi), fm=q(mid);
    if ((flo<0)!=(fm<0)){hi=mid;fhi=fm;} else {lo=mid;flo=fm;}
  }
  (void)fhi;
  return 0.5*(lo+hi);
}

swcme::shock::JumpResult oblique_fixture_result(){
  return swcme::shock::solve_ideal_mhd_fast_shock(fixture(40.0),{{1,0,0}},1.2e6,GAMMA);
}

} // namespace

void test_shk01(swcme_test::Context& context){
  std::cout<<"SHK01 fast-shock existence and no-shock threshold\n";
  const Primitive up=fixture(45.0);
  const Vec3 n={{1,0,0}};
  const double cf=reference_fast_speed(up,n);
  const double V1n=up.velocity_m_s[0];

  for(double mf : {0.95,1.0,1.05}){
    const auto r=swcme::shock::solve_ideal_mhd_fast_shock(up,n,V1n+mf*cf,GAMMA);
    const bool expect_shock=mf>1.0;
    context.expect_true(r.has_shock==expect_shock,"fast-shock threshold classification");
    if(!expect_shock){
      expect_abs(context,"no-shock compression",r.compression,1.0,0.0);
      expect_rel(context,"no-shock downstream rho",r.downstream.rho_kg_m3,up.rho_kg_m3,1e-15);
      expect_rel(context,"no-shock downstream pressure",r.downstream.pressure_Pa,up.pressure_Pa,1e-15);
    } else {
      context.expect_true(r.solver_converged,"super-fast fixture must converge");
      context.expect_true(r.compression>1.0,"super-fast fixture must compress");
    }
  }

  // Permanent regression guard for the legacy 1-D bug: an empirical sheath
  // compression floor must not create a shock when V_sh=V_sw.
  swcme1d::Params p1;
  p1.V0_sh_kms=p1.V_sw_kms;
  p1.sheath_comp_floor=1.8;
  const swcme1d::StepState s1=swcme1d::Model(p1).prepare_step(0.0);
  context.expect_true(!s1.has_shock,"1-D floor must not manufacture a shock");
  expect_abs(context,"1-D no-shock rc despite floor",s1.rc,1.0,0.0);

  swcme3d::Params p3;
  p3.shape=swcme3d::ShockShape::Sphere;
  p3.V0_sh_kms=p3.V_sw_kms;
  p3.sheath_comp_floor=1.8;
  const swcme3d::Model m3(p3);
  const auto s3=m3.prepare_step(0.0);
  swcme3d::LocalShockState js;
  const double u[3]={1,0,0};
  context.expect_true(m3.shock_state_direction(s3,u,js),"3-D sphere surface exists");
  context.expect_true(!js.has_shock,"3-D floor must not manufacture a shock");
  expect_abs(context,"3-D no-shock rc despite floor",js.compression,1.0,0.0);
}

void test_shk02(swcme_test::Context& context){
  std::cout<<"SHK02 shock obliquity and polarity invariance\n";
  const Vec3 n={{1,0,0}};
  for(double theta : {0.0,30.0,60.0,90.0}){
    Primitive up=fixture(theta);
    auto a=swcme::shock::solve_ideal_mhd_fast_shock(up,n,1.2e6,GAMMA);
    for(double& b:up.magnetic_T)b=-b;
    auto b=swcme::shock::solve_ideal_mhd_fast_shock(up,n,1.2e6,GAMMA);
    expect_abs(context,"theta_Bn polarity invariant",a.theta_Bn_rad,b.theta_Bn_rad,2e-14);
    expect_abs(context,"theta_Bn known fixture",a.theta_Bn_rad,
               theta*swcme::constants::PI/180.0,2e-14);
  }
}

void test_shk03(swcme_test::Context& context){
  std::cout<<"SHK03 parallel fast-shock analytical limit\n";
  Primitive up=fixture(0.0);
  const Vec3 n={{1,0,0}};
  const double Vsh=9.0e5;
  const double U=Vsh-up.velocity_m_s[0];
  const double cs=std::sqrt(GAMMA*up.pressure_Pa/up.rho_kg_m3);
  const double Ms=U/cs;
  const double rref=((GAMMA+1.0)*Ms*Ms)/((GAMMA-1.0)*Ms*Ms+2.0);
  const double p_ratio=(2.0*GAMMA*Ms*Ms-(GAMMA-1.0))/(GAMMA+1.0);
  const double p2ref=p_ratio*up.pressure_Pa;
  const double V2ref=Vsh-(Vsh-up.velocity_m_s[0])/rref;

  const auto r=swcme::shock::solve_ideal_mhd_fast_shock(up,n,Vsh,GAMMA);
  context.expect_true(r.has_shock&&r.solver_converged,"parallel fast shock must converge");
  expect_rel(context,"parallel compression",r.compression,rref,2e-10);
  expect_rel(context,"parallel downstream pressure",r.downstream.pressure_Pa,p2ref,2e-10);
  expect_rel(context,"parallel downstream normal speed",r.downstream.velocity_m_s[0],V2ref,2e-10);
}

void test_shk04(swcme_test::Context& context){
  std::cout<<"SHK04 perpendicular MHD-shock independent benchmark\n";
  Primitive up=fixture(90.0);
  // Eliminate tiny cos(pi/2) roundoff so this is exactly B_n=0.
  up.magnetic_T={{0.0,5.0e-9,0.0}};
  const double Vsh=8.0e5;
  const double rref=perpendicular_reference_compression(up,Vsh);
  context.expect_true(rref>1.0,"independent perpendicular reference root exists");
  const auto r=swcme::shock::solve_ideal_mhd_fast_shock(up,{{1,0,0}},Vsh,GAMMA);
  context.expect_true(r.has_shock&&r.solver_converged,"perpendicular shock must converge");
  expect_rel(context,"perpendicular compression",r.compression,rref,2e-9);
  expect_rel(context,"perpendicular B_t compression",r.downstream.magnetic_T[1],
             rref*up.magnetic_T[1],2e-9);
  expect_rel(context,"perpendicular density compression",r.downstream.rho_kg_m3,
             rref*up.rho_kg_m3,2e-9);
}

void test_shk05(swcme_test::Context& context){
  std::cout<<"SHK05 independent oblique-MHD benchmark\n";
  namespace reference=swcme_test::shk05_reference_v1;
  constexpr std::size_t case_count=sizeof(reference::CASES)/
                                   sizeof(reference::CASES[0]);

  // These metadata assertions prevent a future edit from quietly replacing
  // the promised high-precision, versioned campaign with a few hand-entered
  // doubles.  The generator solves all eight RH equations with Decimal
  // precision 80 and records how many independent Newton seeds reached the
  // unique evolutionary fast branch.
  context.expect_true(reference::FIXTURE_VERSION==1,
                      "SHK05 fixture schema/version is pinned to v1");
  context.expect_true(reference::DECIMAL_PRECISION_DIGITS>=50,
                      "SHK05 reference solve uses at least 50 decimal digits");
  context.expect_true(case_count>=12,
                      "SHK05 matrix spans at least twelve oblique states");

  for(const reference::Case& expected : reference::CASES){
    Primitive upstream;
    upstream.rho_kg_m3=expected.upstream_rho_kg_m3;
    upstream.pressure_Pa=expected.upstream_pressure_Pa;
    for(int component=0;component<3;++component){
      upstream.velocity_m_s[component]=expected.upstream_velocity_m_s[component];
      upstream.magnetic_T[component]=expected.upstream_magnetic_T[component];
    }
    const Vec3 normal={{1.0,0.0,0.0}};
    const auto actual=swcme::shock::solve_ideal_mhd_fast_shock(
        upstream,normal,expected.shock_speed_m_s,expected.gamma);
    const std::string prefix=std::string("SHK05 ")+expected.id+" ";

    // Fixture-generation diagnostics are release data, not merely comments.
    // Every case must have a unique admissible root reached from multiple
    // initial guesses, an 80-digit equation residual below 1e-50, and a
    // nonsingular final Newton path.  A fixture that fails these checks is too
    // poorly conditioned to distinguish a production defect from reference
    // uncertainty and must be regenerated as a new version.
    context.expect_true(std::string(expected.branch)=="EVOLUTIONARY_FAST",
                        prefix+"reference branch is evolutionary fast");
    context.expect_true(expected.physical_root_count==1 &&
                            expected.physical_seed_count>=2 &&
                            expected.converged_seed_count>=
                                expected.physical_seed_count,
                        prefix+"independent seeds select one physical root");
    context.expect_true(expected.reference_max_residual<1.0e-50 &&
                            expected.minimum_newton_pivot>1.0e-10,
                        prefix+"reference residual and conditioning are qualified");

    context.expect_true(actual.has_shock && actual.solver_converged &&
                            actual.status==swcme::shock::SolveStatus::Solved,
                        prefix+"production solver selects a solved fast shock");
    if(!actual.solver_converged) continue;

    // Compare every downstream primitive, not only compression.  The 1e-7
    // validation-plan requirement is tightened to 2e-9 for these deliberately
    // well-conditioned fixtures.  Component-specific floors keep a true zero
    // from being assigned a meaningless relative error while still detecting
    // spurious cross-plane velocity or magnetic field.
    expect_rel(context,prefix+"upstream fast Mach",actual.fast_mach,
               expected.fast_mach,2.0e-10,1.0);
    expect_rel(context,prefix+"compression",actual.compression,
               expected.compression,2.0e-9,1.0);
    expect_rel(context,prefix+"downstream density",
               actual.downstream.rho_kg_m3,expected.downstream_rho_kg_m3,
               2.0e-9,1.0e-30);
    expect_rel(context,prefix+"downstream pressure",
               actual.downstream.pressure_Pa,expected.downstream_pressure_Pa,
               2.0e-9,1.0e-30);
    for(int component=0;component<3;++component){
      expect_rel(context,prefix+"downstream velocity["+
                     std::to_string(component)+"]",
                 actual.downstream.velocity_m_s[component],
                 expected.downstream_velocity_m_s[component],2.0e-9,1.0);
      expect_rel(context,prefix+"downstream magnetic["+
                     std::to_string(component)+"]",
                 actual.downstream.magnetic_T[component],
                 expected.downstream_magnetic_T[component],2.0e-9,1.0e-20);
    }
    expect_rel(context,prefix+"entropy ratio",actual.entropy_ratio,
               expected.entropy_ratio,2.0e-9,1.0);

    // Independently classify the returned branch from downstream
    // characteristics.  A production result is evolutionary fast only if
    // the upstream normal flow is super-fast, the downstream normal flow is
    // sub-fast but remains super-Alfvenic, and entropy increases.
    const double downstream_fast_speed=reference_fast_speed(
        actual.downstream,normal,expected.gamma);
    const double downstream_normal_flow=std::abs(
        actual.downstream.velocity_m_s[0]-expected.shock_speed_m_s);
    const double downstream_fast_mach=
        downstream_normal_flow/downstream_fast_speed;
    const double downstream_normal_alfven=
        std::abs(actual.downstream.magnetic_T[0])/
        std::sqrt(MU0*actual.downstream.rho_kg_m3);
    expect_rel(context,prefix+"downstream fast Mach",downstream_fast_mach,
               expected.downstream_fast_mach,2.0e-9,1.0);
    context.expect_true(actual.fast_mach>1.0 && downstream_fast_mach<1.0 &&
                            downstream_normal_flow>downstream_normal_alfven &&
                            actual.entropy_ratio>1.0,
                        prefix+"production branch is evolutionary fast");

    // Production conservation diagnostics are recomputed from the accepted
    // primitives by swcme_shock.hpp and are deliberately separate from the
    // frozen reference.  Requiring all five here makes a full-state mismatch
    // retain useful evidence about which invariant failed.
    context.expect_true(actual.mass_residual<=1.0e-9 &&
                            actual.normal_B_residual<=1.0e-10 &&
                            actual.electric_residual<=1.0e-8 &&
                            actual.momentum_residual<=1.0e-8 &&
                            actual.energy_residual<=1.0e-8,
                        prefix+"all production RH residuals pass");
  }
}

void test_shk06(swcme_test::Context& context){
  std::cout<<"SHK06 Rankine-Hugoniot mass-flux conservation\n";
  const auto r=oblique_fixture_result();
  context.expect_true(r.has_shock&&r.solver_converged,"reference oblique shock converged");
  expect_abs(context,"mass residual",r.mass_residual,0.0,1e-10);

  // Integration check: RESOLVED_COMPRESSION reaches the exact RH state at
  // the inner edge of its finite C1 shock layer.  The mathematical shock state
  // itself remains exact through shock_state_direction().
  swcme3d::Params p;
  p.shape=swcme3d::ShockShape::Sphere; p.r0_Rs=20.0; p.V0_sh_kms=1200.0;
  p.V_sw_kms=400.0; p.Gamma_kmInv=1e-8; p.sin_theta=1.0;
  // The canonical default is SOURCE/SHOCK_ONLY after Fix 14. This test is
  // specifically a resolved-compression integration check, so both members of
  // the alternative FULL_ICME/RESOLVED_COMPRESSION pair are explicit.
  p.region_mode=swcme::regions::Mode::FullICME;
  p.shock_acceleration_mode=swcme::acceleration::Mode::ResolvedCompression;
  const swcme3d::Model m(p); const auto s=m.prepare_step(0.0);
  const double u[3]={1,0,0}; swcme3d::LocalShockState js;
  context.expect_true(m.shock_state_direction(s,u,js)&&js.has_shock&&js.solver_converged,
                      "3-D apex shock state exists");
  const auto boundaries=swcme::regions::make_boundaries(js.Rdir_m,s.region_config);
  const double rq=js.Rdir_m-0.5*boundaries.smooth_shock_width_m;
  const double x=rq, y=0.0,z=0.0;
  double nq=0,vx=0,vy=0,vz=0,bx=0,by=0,bz=0;
  m.evaluate_cartesian_with_B(s,&x,&y,&z,&nq,&vx,&vy,&vz,&bx,&by,&bz,1);
  expect_rel(context,"3-D resolved-shock inner density",nq,js.downstream_n_m3,2e-7);
  expect_rel(context,"3-D resolved-shock inner Vx",vx,js.downstream.velocity_m_s[0],2e-7,1.0);
  expect_rel(context,"3-D resolved-shock inner Vy",vy,js.downstream.velocity_m_s[1],2e-7,1.0);
  expect_rel(context,"3-D resolved-shock inner Bx",bx,js.downstream.magnetic_T[0],2e-7,1e-30);
  expect_rel(context,"3-D resolved-shock inner By",by,js.downstream.magnetic_T[1],2e-7,1e-30);
}

void test_shk07(swcme_test::Context& context){
  std::cout<<"SHK07 Rankine-Hugoniot normal magnetic-field continuity\n";
  const auto r=oblique_fixture_result();
  context.expect_true(r.has_shock&&r.solver_converged,"reference oblique shock converged");
  expect_abs(context,"normal-B residual",r.normal_B_residual,0.0,1e-12);
}

void test_shk08(swcme_test::Context& context){
  std::cout<<"SHK08 Rankine-Hugoniot tangential electric-field conservation\n";
  const auto r=oblique_fixture_result();
  context.expect_true(r.has_shock&&r.solver_converged,"reference oblique shock converged");
  expect_abs(context,"tangential-E residual",r.electric_residual,0.0,1e-9);
}

void test_shk09(swcme_test::Context& context){
  std::cout<<"SHK09 Rankine-Hugoniot momentum-flux conservation\n";
  const auto r=oblique_fixture_result();
  context.expect_true(r.has_shock&&r.solver_converged,"reference oblique shock converged");
  expect_abs(context,"momentum residual",r.momentum_residual,0.0,1e-9);
}

void test_shk10(swcme_test::Context& context){
  std::cout<<"SHK10 Rankine-Hugoniot total-energy-flux conservation\n";
  const auto r=oblique_fixture_result();
  context.expect_true(r.has_shock&&r.solver_converged,"reference oblique shock converged");
  expect_abs(context,"energy residual",r.energy_residual,0.0,1e-8);
}

void test_shk11(swcme_test::Context& context){
  std::cout<<"SHK11 physical admissibility and entropy increase\n";
  for(double theta : {0.0,30.0,60.0,90.0}){
    Primitive up=fixture(theta);
    if(theta==90.0) up.magnetic_T={{0.0,5e-9,0.0}};
    const auto r=swcme::shock::solve_ideal_mhd_fast_shock(up,{{1,0,0}},8e5,GAMMA);
    context.expect_true(r.has_shock&&r.solver_converged,"admissible shock converged");
    context.expect_true(r.downstream.rho_kg_m3>up.rho_kg_m3,"rho2>rho1");
    context.expect_true(r.downstream.pressure_Pa>0.0,"p2 positive");
    context.expect_true(r.compression>1.0&&r.compression<=4.0*(1+1e-12),"compression within gamma=5/3 bound");
    context.expect_true(r.entropy_ratio>=1.0-1e-10,"entropy proxy increases");
  }
}

void test_shk12(swcme_test::Context& context){
  std::cout<<"SHK12 two-sided near-Mach-one shock limit\n";
  namespace reference=swcme_test::shk12_reference_v1;
  constexpr std::size_t case_count=sizeof(reference::CASES)/
                                   sizeof(reference::CASES[0]);
  const Vec3 n={{1,0,0}};

  // The frozen values are generated with an 80-digit, eight-variable Newton
  // solve and logarithmic continuation.  Pinning the schema and matrix size
  // prevents a future fixture edit from silently reducing the required beta,
  // obliquity, gamma, or Mach-excess coverage.
  context.expect_true(reference::FIXTURE_VERSION==1,
                      "SHK12 fixture schema/version is pinned to v1");
  context.expect_true(reference::DECIMAL_PRECISION_DIGITS>=50,
                      "SHK12 reference uses at least 50 decimal digits");
  context.expect_true(case_count>=48,
                      "SHK12 reference contains eight six-point families");

  double previous_compression=0.0;
  double previous_excess=0.0;
  for(std::size_t index=0;index<case_count;++index){
    const reference::Case& expected=reference::CASES[index];
    Primitive upstream;
    upstream.rho_kg_m3=expected.upstream_rho_kg_m3;
    upstream.pressure_Pa=expected.upstream_pressure_Pa;
    for(int component=0;component<3;++component){
      upstream.velocity_m_s[component]=expected.upstream_velocity_m_s[component];
      upstream.magnetic_T[component]=expected.upstream_magnetic_T[component];
    }
    const auto actual=swcme::shock::solve_ideal_mhd_fast_shock(
        upstream,n,expected.shock_speed_m_s,expected.gamma);
    const std::string prefix=std::string("SHK12 ")+expected.id+" ";

    // Reference conditioning and residual metadata are executable acceptance
    // criteria.  A changed fixture therefore cannot hide a poorly converged
    // high-precision solve behind rounded binary64 values.
    context.expect_true(expected.reference_max_residual<1.0e-50,
                        prefix+"reference residual is below 1e-50");
    context.expect_true(expected.minimum_newton_pivot>1.0e-12,
                        prefix+"reference Jacobian remains conditioned");
    context.expect_true(actual.status==swcme::shock::SolveStatus::Solved &&
                            actual.has_shock && actual.solver_converged,
                        prefix+"well-conditioned weak branch is solved");

    // The plan requires the full primitive state, not compression alone.  The
    // 2e-7 component tolerance leaves a small binary64 bracket margin while
    // remaining independent of and much tighter than visually smooth trends.
    expect_rel(context,prefix+"Mach excess",actual.fast_mach-1.0,
               expected.mach_excess,2.0e-7,1.0e-15);
    expect_rel(context,prefix+"compression",actual.compression,
               expected.compression,2.0e-7);
    expect_rel(context,prefix+"downstream density",actual.downstream.rho_kg_m3,
               expected.downstream_rho_kg_m3,2.0e-7,1.0e-30);
    expect_rel(context,prefix+"downstream pressure",actual.downstream.pressure_Pa,
               expected.downstream_pressure_Pa,2.0e-7,1.0e-30);
    for(int component=0;component<3;++component){
      expect_rel(context,prefix+"downstream velocity["+
                     std::to_string(component)+"]",
                 actual.downstream.velocity_m_s[component],
                 expected.downstream_velocity_m_s[component],2.0e-7,1.0);
      expect_rel(context,prefix+"downstream magnetic["+
                     std::to_string(component)+"]",
                 actual.downstream.magnetic_T[component],
                 expected.downstream_magnetic_T[component],2.0e-7,1.0e-15);
    }

    // Published iteration/bracket diagnostics prove scalar convergence and
    // make any future weak-limit regression directly reproducible.
    context.expect_true(actual.root_iterations>0 && actual.root_iterations<=120,
                        prefix+"iteration count is bounded");
    context.expect_true(actual.root_bracket_lower_compression<=actual.compression &&
                            actual.compression<=actual.root_bracket_upper_compression &&
                            actual.root_bracket_width>=0.0 &&
                            actual.root_bracket_width<=3.0e-13*
                                std::max(1.0,actual.compression),
                        prefix+"final bracket contains the reported root");

    // Each family is stored from larger to smaller Mach excess.  Compression
    // must decrease toward one without an empirical floor; the final 1e-5
    // reference is close enough to expose even a small imposed compression.
    if(index%6!=0){
      context.expect_true(expected.mach_excess<previous_excess &&
                              actual.compression<previous_compression,
                          prefix+"compression approaches one monotonically");
    }
    if(index%6==5){
      context.expect_true(actual.compression-1.0<2.0e-5,
                          prefix+"smallest resolved shock approaches identity");
    }
    previous_excess=expected.mach_excess;
    previous_compression=actual.compression;
  }

  auto make_sweep_state=[](double beta,double theta_deg){
    Primitive state;
    const double B=5.0e-9;
    const double theta=theta_deg*swcme::constants::PI/180.0;
    const double phi=25.0*swcme::constants::PI/180.0;
    state.rho_kg_m3=5.0e6*MP;
    state.pressure_Pa=beta*B*B/(2.0*MU0);
    // Zero normal flow makes the exactly critical floating-point fixture
    // unambiguous; nonzero tangential flow still exercises the vector jump.
    state.velocity_m_s={{0.0,2.0e4,-1.0e4}};
    state.magnetic_T={{B*std::cos(theta),
                       B*std::sin(theta)*std::cos(phi),
                       B*std::sin(theta)*std::sin(phi)}};
    return state;
  };

  // Sweep both sides of Mfast=1 down to 1e-12 for every beta/angle/gamma
  // family.  Subcritical and exactly critical states are physical NoShock;
  // small positive excesses are physical shocks that are explicitly marked
  // numerically unresolved, never silently folded into NoShock.
  for(double beta : {0.01,0.1,1.0,10.0}){
    for(double theta : {1.0,30.0,60.0,89.0}){
      for(double gamma : {1.4,5.0/3.0}){
        const Primitive upstream=make_sweep_state(beta,theta);
        const double fast=reference_fast_speed(upstream,n,gamma);
        for(double excess : {-0.5,-0.1,-1.0e-2,-1.0e-4,-1.0e-8,-1.0e-12,0.0}){
          const auto result=swcme::shock::solve_ideal_mhd_fast_shock(
              upstream,n,(1.0+excess)*fast,gamma);
          context.expect_true(result.status==swcme::shock::SolveStatus::NoShock &&
                                  !result.has_shock && result.solver_converged &&
                                  result.compression==1.0,
                              "SHK12 subcritical/critical state is NoShock");
        }
        for(double excess : {1.0e-12,1.0e-10,1.0e-8,1.0e-7,1.0e-6}){
          const auto result=swcme::shock::solve_ideal_mhd_fast_shock(
              upstream,n,(1.0+excess)*fast,gamma);
          context.expect_true(
              result.status==swcme::shock::SolveStatus::NumericallyUnresolvedWeakShock &&
                  result.has_shock && !result.solver_converged &&
                  std::isfinite(result.downstream.rho_kg_m3) &&
                  std::isfinite(result.downstream.pressure_Pa),
              "SHK12 sub-resolution supercritical state is explicitly unresolved");
        }
      }
    }
  }

  // These low-beta, nearly parallel points expose a near-singular tangential
  // system.  The old outermost-root rule returned a discontinuous r~4--6
  // state.  They must retain a finite diagnostic candidate but carry the
  // unresolved status so downstream SEP physics cannot consume that branch.
  for(const auto parameters :
      {std::array<double,3>{{0.01,1.4,1.0e-2}},
       std::array<double,3>{{0.01,5.0/3.0,1.0e-1}}}){
    const Primitive upstream=make_sweep_state(parameters[0],1.0);
    const double fast=reference_fast_speed(upstream,n,parameters[1]);
    const auto result=swcme::shock::solve_ideal_mhd_fast_shock(
        upstream,n,(1.0+parameters[2])*fast,parameters[1]);
    context.expect_true(
        result.status==swcme::shock::SolveStatus::NumericallyUnresolvedWeakShock &&
            result.has_shock && !result.solver_converged &&
            std::isfinite(result.compression) &&
            std::isfinite(result.downstream.rho_kg_m3) &&
            std::isfinite(result.downstream.pressure_Pa),
        "SHK12 discontinuous near-singular branch is explicitly unresolved");
  }
}
