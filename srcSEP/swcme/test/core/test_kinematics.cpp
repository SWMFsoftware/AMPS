#include "test_framework.hpp"

#include <swcme1d.hpp>
#include <swcme3d.hpp>
#include <swcme_kinematics.hpp>

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

constexpr double RS=swcme::constants::SOLAR_RADIUS_M;

void check_rel(swcme_test::Context& context,const std::string& label,
               double actual,double expected,double rel_tol,double abs_floor=1e-30){
  const double scale=std::max({std::abs(actual),std::abs(expected),abs_floor});
  const double err=std::abs(actual-expected)/scale;
  const bool pass=std::isfinite(actual)&&std::isfinite(expected)&&err<=rel_tol;
  context.expect_true(pass,label);
  std::cout<<"  "<<std::left<<std::setw(44)<<label<<std::right
           <<" actual="<<std::scientific<<std::setprecision(12)<<actual
           <<" ref="<<expected<<" rel_err="<<err<<" tol="<<rel_tol
           <<" "<<(pass?"PASS":"FAIL")<<'\n';
}


// Independent closed-form DBM reference.  This helper intentionally does not
// call swcme::kinematics::dbm_state(); it reproduces the governing analytical
// solution directly so KIN01/KIN02 validate the production implementation.
void reference_dbm(double r0,double V0,double Vsw,double gamma,double t,
                   double& r,double& v){
  const double dv0=V0-Vsw;
  const double a=std::abs(dv0);
  if(gamma==0.0 || a==0.0){r=r0+V0*t;v=V0;return;}
  const double x=gamma*a*t;
  v=Vsw+dv0/(1.0+x);
  const double sign=(dv0>0.0)?1.0:-1.0;
  r=r0+Vsw*t+sign*std::log1p(x)/gamma;
}

swcme::kinematics::Config basic_dbm(double V0_kms,double Vsw_kms,double gamma_km_inv){
  swcme::kinematics::Config c;
  c.mode=swcme::kinematics::Mode::DBM;
  c.r0_m=20.0*RS;
  c.V0_m_s=V0_kms*1e3;
  c.Vsw_m_s=Vsw_kms*1e3;
  c.Gamma_m_inv=gamma_km_inv/1e3;
  return c;
}

void compare_wrappers(swcme_test::Context& context,double t_s,double r_ref,double v_ref,
                      const swcme::kinematics::Config& c){
  // Exercise both public dimensional interfaces, not only the common helper.
  // The wrappers use the same SI core and therefore should agree to roundoff.
  swcme1d::Params p1;
  p1.kinematics_mode=c.mode;
  p1.r0_Rs=c.r0_m/RS;
  p1.V0_sh_kms=c.V0_m_s/1e3;
  p1.V_sw_kms=c.Vsw_m_s/1e3;
  p1.Gamma_kmInv=c.Gamma_m_inv*1e3;
  const auto s1=swcme1d::Model(p1).prepare_step(t_s);

  swcme3d::Params p3;
  p3.kinematics_mode=c.mode;
  p3.r0_Rs=c.r0_m/RS;
  p3.V0_sh_kms=c.V0_m_s/1e3;
  p3.V_sw_kms=c.Vsw_m_s/1e3;
  p3.Gamma_kmInv=c.Gamma_m_inv*1e3;
  const auto s3=swcme3d::Model(p3).prepare_step(t_s);

  check_rel(context,"1-D radius versus reference",s1.r_sh_m,r_ref,1e-12,1.0);
  check_rel(context,"3-D radius versus reference",s3.r_sh_m,r_ref,1e-12,1.0);
  check_rel(context,"1-D speed versus reference",s1.V_sh_ms,v_ref,1e-12,1.0);
  check_rel(context,"3-D speed versus reference",s3.V_sh_ms,v_ref,1e-12,1.0);
  check_rel(context,"1-D versus 3-D radius",s1.r_sh_m,s3.r_sh_m,2e-15,1.0);
  check_rel(context,"1-D versus 3-D speed",s1.V_sh_ms,s3.V_sh_ms,2e-15,1.0);
}

swcme::kinematics::Config data_fixture(){
  swcme::kinematics::Config c;
  c.mode=swcme::kinematics::Mode::DataDriven;
  c.r0_m=20.0*RS;
  c.V0_m_s=1000e3;
  c.Vsw_m_s=400e3;
  c.Gamma_m_inv=1e-10;
  c.data_time_s={0.0,3600.0,9000.0,18000.0,36000.0};
  c.data_radius_m={20.0*RS,25.0*RS,31.0*RS,39.0*RS,50.0*RS};
  return c;
}

} // namespace

void test_kin01(swcme_test::Context& context){
  std::cout<<"KIN01 fast-CME sign-aware DBM closed form\n";
  const auto c=basic_dbm(1500.0,400.0,8e-8);
  for(double t : {0.0,3600.0,6.0*3600.0,24.0*3600.0,72.0*3600.0}){
    double rr=0.0,vr=0.0; reference_dbm(c.r0_m,c.V0_m_s,c.Vsw_m_s,c.Gamma_m_inv,t,rr,vr);
    const auto s=swcme::kinematics::evaluate(c,t);
    context.expect_true(s.status==swcme::kinematics::Status::Ok,"fast DBM status");
    check_rel(context,"common fast DBM radius",s.radius_m,rr,1e-12,1.0);
    check_rel(context,"common fast DBM speed",s.speed_m_s,vr,1e-12,1.0);
    compare_wrappers(context,t,rr,vr,c);
  }
}

void test_kin02(swcme_test::Context& context){
  std::cout<<"KIN02 slow-CME sign-aware DBM branch\n";
  const auto c=basic_dbm(300.0,400.0,8e-8);
  double previous=c.V0_m_s;
  for(double t : {0.0,3600.0,6.0*3600.0,24.0*3600.0,72.0*3600.0}){
    double rr=0.0,vr=0.0; reference_dbm(c.r0_m,c.V0_m_s,c.Vsw_m_s,c.Gamma_m_inv,t,rr,vr);
    const auto s=swcme::kinematics::evaluate(c,t);
    context.expect_true(s.status==swcme::kinematics::Status::Ok,"slow DBM status");
    check_rel(context,"common slow DBM radius",s.radius_m,rr,1e-12,1.0);
    check_rel(context,"common slow DBM speed",s.speed_m_s,vr,1e-12,1.0);
    context.expect_true(s.speed_m_s>=previous-1e-10,"slow CME accelerates monotonically toward Vsw");
    context.expect_true(s.speed_m_s<=c.Vsw_m_s+1e-10,"slow CME does not overshoot Vsw");
    previous=s.speed_m_s;
    compare_wrappers(context,t,rr,vr,c);
  }
}

void test_kin03(swcme_test::Context& context){
  std::cout<<"KIN03 exact zero-drag ballistic limit\n";
  auto c=basic_dbm(1500.0,400.0,0.0);
  for(double t : {0.0,1.0,3600.0,24.0*3600.0}){
    const double rr=c.r0_m+c.V0_m_s*t;
    const auto s=swcme::kinematics::evaluate(c,t);
    context.expect_true(s.status==swcme::kinematics::Status::Ok,"Gamma=0 status");
    check_rel(context,"Gamma=0 radius",s.radius_m,rr,2e-15,1.0);
    check_rel(context,"Gamma=0 speed",s.speed_m_s,c.V0_m_s,0.0,1.0);
    compare_wrappers(context,t,rr,c.V0_m_s,c);
  }
}

void test_kin04(swcme_test::Context& context){
  std::cout<<"KIN04 small-Gamma continuity\n";
  const double t=24.0*3600.0;
  const double V0=1500e3, Vsw=400e3;
  const double a=std::abs(V0-Vsw);

  // The production DBM switches from the direct log1p expression to a Taylor
  // representation at x=Gamma*a*t=1e-8.  KIN04 samples immediately on both
  // sides of that NUMERICAL branch boundary.  It must not compare physically
  // different drag coefficients separated by orders of magnitude; doing so
  // would confuse a real physical Gamma dependence with a numerical jump.
  const double gamma_threshold_m=1.0e-8/(a*t);
  const double gamma_threshold_km=gamma_threshold_m*1.0e3;
  const double eps=1.0e-7;
  const std::vector<double> gamma_km={
      0.0,
      gamma_threshold_km*(1.0-10.0*eps),
      gamma_threshold_km*(1.0-eps),
      gamma_threshold_km*(1.0+eps),
      gamma_threshold_km*(1.0+10.0*eps)};

  std::vector<swcme::kinematics::State> states;
  for(double g:gamma_km){
    const auto c=basic_dbm(1500.0,400.0,g);
    const auto st=swcme::kinematics::evaluate(c,t);
    context.expect_true(st.status==swcme::kinematics::Status::Ok,"small-Gamma status");
    double rr=0.0,vr=0.0;
    reference_dbm(c.r0_m,c.V0_m_s,c.Vsw_m_s,c.Gamma_m_inv,t,rr,vr);
    check_rel(context,"small-Gamma analytical radius",st.radius_m,rr,2e-13,1.0);
    check_rel(context,"small-Gamma analytical speed",st.speed_m_s,vr,2e-13,1.0);
    states.push_back(st);
  }

  // Compare the two points that straddle x=1e-8.  Their Gamma values differ
  // only by 2e-7 relative, so the physically expected change is tiny; any
  // larger discontinuity would identify a mismatch between the series and
  // direct-log branches.
  const auto& below=states[2];
  const auto& above=states[3];
  const double rjump=std::abs(above.radius_m-below.radius_m)/
                     std::max(std::abs(above.radius_m),std::abs(below.radius_m));
  const double vjump=std::abs(above.speed_m_s-below.speed_m_s)/
                     std::max(std::abs(above.speed_m_s),std::abs(below.speed_m_s));
  std::cout<<"  branch-boundary radius relative jump="<<std::scientific<<rjump<<'\n'
           <<"  branch-boundary speed relative jump ="<<vjump<<'\n';
  context.expect_true(rjump<1e-11,"no numerical radius jump across small-Gamma branch");
  context.expect_true(vjump<1e-11,"no numerical speed jump across small-Gamma branch");
}

void test_kin05(swcme_test::Context& context){
  std::cout<<"KIN05 long-time DBM asymptotic stability\n";
  for(double V0 : {300.0,1500.0}){
    const auto c=basic_dbm(V0,400.0,8e-8);
    double prev_r=0.0;
    double prev_abs_dv=std::abs(c.V0_m_s-c.Vsw_m_s)+1.0;
    for(double t : {0.0,1e5,1e6,1e7,1e8,1e9}){
      const auto s=swcme::kinematics::evaluate(c,t);
      context.expect_true(s.status==swcme::kinematics::Status::Ok,"long-time DBM status");
      context.expect_true(std::isfinite(s.radius_m)&&std::isfinite(s.speed_m_s),"long-time state finite");
      context.expect_true(s.radius_m>=prev_r,"radius monotonically increasing");
      const double abs_dv=std::abs(s.speed_m_s-c.Vsw_m_s);
      context.expect_true(abs_dv<=prev_abs_dv+1e-8,"speed approaches Vsw monotonically");
      prev_r=s.radius_m;prev_abs_dv=abs_dv;
    }
    const auto late=swcme::kinematics::evaluate(c,1e12);
    context.expect_true(late.status==swcme::kinematics::Status::Ok,"very-late DBM status");
    context.expect_true(std::abs(late.speed_m_s-c.Vsw_m_s)<1.0,"very-late speed approaches Vsw within 1 m/s");
  }
}

void test_kin06(swcme_test::Context& context){
  std::cout<<"KIN06 data-driven PCHIP knot exactness\n";
  const auto c=data_fixture();
  for(std::size_t i=0;i<c.data_time_s.size();++i){
    const auto s=swcme::kinematics::evaluate(c,c.data_time_s[i]);
    context.expect_true(s.status==swcme::kinematics::Status::Ok,"data knot status");
    check_rel(context,"PCHIP knot radius",s.radius_m,c.data_radius_m[i],1e-13,1.0);
    context.expect_true(s.speed_m_s>=0.0,"PCHIP knot derivative nonnegative");

    // Independent centered/one-sided finite-difference derivative of the
    // production interpolated radius checks that the reported speed is the
    // derivative of the same curve rather than an unrelated secant estimate.
    if(i>0 && i+1<c.data_time_s.size()){
      const double min_gap=std::min(c.data_time_s[i]-c.data_time_s[i-1],
                                    c.data_time_s[i+1]-c.data_time_s[i]);
      const double h=1e-5*min_gap;
      const auto a=swcme::kinematics::evaluate(c,c.data_time_s[i]-h);
      const auto b=swcme::kinematics::evaluate(c,c.data_time_s[i]+h);
      const double fd=(b.radius_m-a.radius_m)/(2*h);
      check_rel(context,"PCHIP knot derivative consistency",s.speed_m_s,fd,2e-5,1e-6);
    }
  }

  // Verify the public wrappers consume the same PCHIP trajectory.
  swcme1d::Params p1; p1.kinematics_mode=swcme::kinematics::Mode::DataDriven;
  p1.data_time_s=c.data_time_s; for(double r:c.data_radius_m)p1.data_radius_Rs.push_back(r/RS);
  swcme3d::Params p3; p3.kinematics_mode=swcme::kinematics::Mode::DataDriven;
  p3.data_time_s=c.data_time_s; for(double r:c.data_radius_m)p3.data_radius_Rs.push_back(r/RS);
  const double tq=9000.0;
  const auto a=swcme1d::Model(p1).prepare_step(tq);
  const auto b=swcme3d::Model(p3).prepare_step(tq);
  check_rel(context,"1-D/3-D data-driven knot radius",a.r_sh_m,b.r_sh_m,2e-15,1.0);
  check_rel(context,"1-D/3-D data-driven knot speed",a.V_sh_ms,b.V_sh_ms,2e-15,1.0);
}

void test_kin07(swcme_test::Context& context){
  std::cout<<"KIN07 data-driven PCHIP monotonicity/no overshoot\n";
  const auto c=data_fixture();
  double prev=-1.0;
  for(int i=0;i<=2000;++i){
    const double f=i/2000.0;
    const double t=c.data_time_s.front()+f*(c.data_time_s.back()-c.data_time_s.front());
    const auto s=swcme::kinematics::evaluate(c,t);
    context.expect_true(s.status==swcme::kinematics::Status::Ok,"PCHIP dense-grid status");
    context.expect_true(s.radius_m>=prev-1e-6,"PCHIP radius monotone");
    context.expect_true(s.speed_m_s>=-1e-10,"PCHIP speed nonnegative");
    // Identify the enclosing knot interval and ensure the interpolation remains
    // inside the endpoint envelope; this directly catches cubic overshoot.
    auto hi=std::lower_bound(c.data_time_s.begin(),c.data_time_s.end(),t);
    std::size_t j=(hi==c.data_time_s.begin())?0:static_cast<std::size_t>(hi-c.data_time_s.begin()-1);
    if(j+1>=c.data_radius_m.size())j=c.data_radius_m.size()-2;
    context.expect_true(s.radius_m>=c.data_radius_m[j]-1e-5 &&
                        s.radius_m<=c.data_radius_m[j+1]+1e-5,
                        "PCHIP stays inside knot envelope");
    prev=s.radius_m;
  }
}

void test_kin08(swcme_test::Context& context){
  std::cout<<"KIN08 explicit data-driven extrapolation policy\n";
  auto c=data_fixture();
  c.extrapolation=swcme::kinematics::ExtrapolationPolicy::OutsideTime;
  const auto before=swcme::kinematics::evaluate(c,-1.0);
  const auto after=swcme::kinematics::evaluate(c,c.data_time_s.back()+1.0);
  context.expect_true(before.status==swcme::kinematics::Status::OutsideTime,
                      "default pre-data query returns OUTSIDE_TIME");
  context.expect_true(after.status==swcme::kinematics::Status::OutsideTime,
                      "default post-data query returns OUTSIDE_TIME");

  c.extrapolation=swcme::kinematics::ExtrapolationPolicy::Ballistic;
  const auto first=swcme::kinematics::evaluate(c,c.data_time_s.front());
  const auto last=swcme::kinematics::evaluate(c,c.data_time_s.back());
  const double dt=600.0;
  const auto b=swcme::kinematics::evaluate(c,c.data_time_s.front()-dt);
  const auto a=swcme::kinematics::evaluate(c,c.data_time_s.back()+dt);
  context.expect_true(b.status==swcme::kinematics::Status::Ok,"ballistic pre-data continuation status");
  context.expect_true(a.status==swcme::kinematics::Status::Ok,"ballistic post-data continuation status");
  check_rel(context,"pre-data ballistic continuation",
            b.radius_m,first.radius_m-first.speed_m_s*dt,2e-15,1.0);
  check_rel(context,"post-data ballistic continuation",
            a.radius_m,last.radius_m+last.speed_m_s*dt,2e-15,1.0);

  // Invalid table inputs must be rejected explicitly; no sorting, duplicate
  // removal or endpoint clamping is performed silently by the production core.
  auto bad=data_fixture(); bad.data_time_s[2]=bad.data_time_s[1];
  context.expect_true(swcme::kinematics::evaluate(bad,3600.0).status==
                      swcme::kinematics::Status::InvalidInput,
                      "duplicate data-driven time rejected");
  bad=data_fixture(); bad.data_radius_m[3]=bad.data_radius_m[2]-RS;
  context.expect_true(swcme::kinematics::evaluate(bad,3600.0).status==
                      swcme::kinematics::Status::InvalidInput,
                      "decreasing data-driven radius rejected");
}
