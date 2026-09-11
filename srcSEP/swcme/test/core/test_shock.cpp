#include "test_framework.hpp"

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
double reference_fast_speed(const Primitive& s,const Vec3& n){
  const double B2=dot(s.magnetic_T,s.magnetic_T);
  const double vA2=B2/(MU0*s.rho_kg_m3);
  const double cs2=GAMMA*s.pressure_Pa/s.rho_kg_m3;
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
  std::cout<<"SHK05 oblique-MHD benchmark grid and branch continuity\n";
  for(double theta : {15.0,30.0,45.0,60.0,75.0}){
    double previous=1.0;
    for(double Vsh : {5.0e5,6.0e5,8.0e5,1.2e6}){
      const auto r=swcme::shock::solve_ideal_mhd_fast_shock(fixture(theta),{{1,0,0}},Vsh,GAMMA);
      context.expect_true(r.has_shock&&r.solver_converged,"oblique benchmark must converge");
      context.expect_true(r.compression>=previous-1e-10,"compression must vary continuously/monotonically with shock speed");
      context.expect_true(r.energy_residual<1e-8,"oblique benchmark energy residual");
      previous=r.compression;
    }
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
  std::cout<<"SHK12 near-Mach-one weak-shock conditioning\n";
  const Primitive up=fixture(60.0);
  const Vec3 n={{1,0,0}};
  const double cf=reference_fast_speed(up,n), V1=up.velocity_m_s[0];
  double previous=10.0;
  for(double delta : {1e-2,5e-3,1e-3}){
    const auto r=swcme::shock::solve_ideal_mhd_fast_shock(up,n,V1+(1.0+delta)*cf,GAMMA);
    context.expect_true(r.has_shock&&r.solver_converged,"weak super-fast shock converged");
    context.expect_true(r.compression>1.0,"weak shock compression above one");
    context.expect_true(r.compression<previous,"compression approaches unity as M_fast approaches one");
    previous=r.compression;
  }
  context.expect_true(previous<1.02,"near-Mach-one compression remains close to unity without floor");
}
