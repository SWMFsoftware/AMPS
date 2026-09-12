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

void test_kin09(swcme_test::Context& context){
  std::cout<<"KIN09 kinematic extrapolation domain\n";
  using swcme::kinematics::Status;

  // BALLISTIC has no algebraic singularity, so negative time is supported
  // until the trajectory crosses 1.05 R_sun.  Log-spaced offsets on both sides
  // of the reference epoch are compared with an independently evaluated
  // long-double line rather than with another production helper.
  swcme::kinematics::Config ballistic;
  ballistic.mode=swcme::kinematics::Mode::Ballistic;
  ballistic.r0_m=20.0*RS;
  ballistic.V0_m_s=1.0e6;
  ballistic.Vsw_m_s=4.0e5;
  ballistic.Gamma_m_inv=0.0;
  for(double t : {-1.0e4,-1.0e3,-1.0e2,-1.0e1,-1.0,
                   0.0,1.0,1.0e1,1.0e2,1.0e3,1.0e4}) {
    const long double expected=
        static_cast<long double>(ballistic.r0_m)+
        static_cast<long double>(ballistic.V0_m_s)*
        static_cast<long double>(t);
    const swcme::kinematics::State state=
        swcme::kinematics::evaluate(ballistic,t);
    context.expect_true(state.status==Status::Ok,
                        "valid ballistic time has OK status");
    check_rel(context,"ballistic independent radius",state.radius_m,
              static_cast<double>(expected),4.0e-15,1.0);
    check_rel(context,"ballistic independent speed",state.speed_m_s,
              ballistic.V0_m_s,0.0,1.0);
  }
  context.expect_true(
      swcme::kinematics::evaluate(ballistic,-2.0e4).status==
          Status::OutsideDomain,
      "ballistic backward crossing returns OUTSIDE_DOMAIN");
  context.expect_true(
      swcme::kinematics::evaluate(
          ballistic,std::numeric_limits<double>::max()).status==
          Status::OutsideDomain,
      "ballistic overflow returns OUTSIDE_DOMAIN");

  // DBM's Gamma=0 branch and equal-speed branch are algebraically ballistic
  // but remain separate production branches.  Exercise them explicitly so a
  // future fast path cannot bypass the common radial/overflow postcondition.
  swcme::kinematics::Config zero_drag=ballistic;
  zero_drag.mode=swcme::kinematics::Mode::DBM;
  context.expect_true(
      swcme::kinematics::evaluate(zero_drag,-2.0e4).status==
          Status::OutsideDomain,
      "zero-drag DBM backward crossing returns OUTSIDE_DOMAIN");
  swcme::kinematics::Config equal_speed=zero_drag;
  equal_speed.Gamma_m_inv=1.0e-10;
  equal_speed.V0_m_s=equal_speed.Vsw_m_s;
  context.expect_true(
      swcme::kinematics::evaluate(
          equal_speed,std::numeric_limits<double>::max()).status==
          Status::OutsideDomain,
      "equal-speed DBM overflow returns OUTSIDE_DOMAIN");
  const swcme::kinematics::State failed_ballistic=
      swcme::kinematics::evaluate(ballistic,-2.0e4);
  context.expect_true(std::isnan(failed_ballistic.radius_m) &&
                          std::isnan(failed_ballistic.speed_m_s),
                      "domain failure exposes no partial ballistic payload");
  context.expect_true(
      swcme::kinematics::evaluate(
          ballistic,std::numeric_limits<double>::quiet_NaN()).status==
          Status::InvalidInput,
      "non-finite time remains INVALID_INPUT");

  // Exercise both signs of DBM velocity contrast.  The independent expression
  // uses long-double log1p and never calls the production DBM implementation;
  // this detects sign, denominator, and branch-status errors simultaneously.
  const std::vector<swcme::kinematics::Config> dbm_cases={
      basic_dbm(1500.0,400.0,8.0e-8),
      basic_dbm(300.0,400.0,8.0e-8)};
  for(const swcme::kinematics::Config& dbm : dbm_cases) {
    for(double t : {-1.0e3,-1.0e2,-1.0e1,-1.0,
                     0.0,1.0,1.0e1,1.0e2,1.0e3,1.0e4}) {
      const long double dv0=
          static_cast<long double>(dbm.V0_m_s)-dbm.Vsw_m_s;
      const long double a=std::abs(dv0);
      const long double x=static_cast<long double>(dbm.Gamma_m_inv)*a*t;
      const long double sign=dv0>0.0L ? 1.0L : -1.0L;
      const long double expected_speed=
          static_cast<long double>(dbm.Vsw_m_s)+dv0/(1.0L+x);
      const long double expected_radius=
          static_cast<long double>(dbm.r0_m)+
          static_cast<long double>(dbm.Vsw_m_s)*t+
          sign*std::log1p(x)/dbm.Gamma_m_inv;
      const swcme::kinematics::State state=
          swcme::kinematics::evaluate(dbm,t);
      context.expect_true(state.status==Status::Ok,
                          "valid DBM extension has OK status");
      check_rel(context,"DBM independent extrapolated radius",state.radius_m,
                static_cast<double>(expected_radius),2.0e-14,1.0);
      check_rel(context,"DBM independent extrapolated speed",state.speed_m_s,
                static_cast<double>(expected_speed),2.0e-14,1.0);
    }
  }

  // For a slow CME, backward continuation has a finite turning point where
  // outward speed reaches zero before the DBM denominator pole.  Sample on
  // both sides: the outward side is valid, while the sunward branch and the
  // pole are explicit domain failures rather than OK or INVALID_INPUT.
  swcme::kinematics::Config turning=basic_dbm(100.0,400.0,1.0e-7);
  turning.r0_m=40.0*RS;
  const double contrast=std::abs(turning.V0_m_s-turning.Vsw_m_s);
  const double turning_time=
      (contrast/turning.Vsw_m_s-1.0)/
      (turning.Gamma_m_inv*contrast);
  const double pole_time=-1.0/(turning.Gamma_m_inv*contrast);
  context.expect_true(
      swcme::kinematics::evaluate(turning,turning_time+1.0).status==Status::Ok,
      "DBM outward side of turning point remains valid");
  context.expect_true(
      swcme::kinematics::evaluate(turning,turning_time-1.0).status==
          Status::OutsideDomain,
      "DBM sunward side of turning point returns OUTSIDE_DOMAIN");
  context.expect_true(
      swcme::kinematics::evaluate(turning,pole_time).status==
          Status::OutsideDomain,
      "DBM denominator pole returns OUTSIDE_DOMAIN");
  context.expect_true(
      swcme::kinematics::evaluate(
          turning,std::numeric_limits<double>::max()).status==
          Status::OutsideDomain &&
      swcme::kinematics::evaluate(
          turning,-std::numeric_limits<double>::max()).status==
          Status::OutsideDomain,
      "DBM overflow-prone extremes return OUTSIDE_DOMAIN");

  // A two-knot table has an independently known endpoint slope.  Explicit
  // BALLISTIC continuation is accepted at/inside the radial boundary and
  // rejected immediately beyond it; the default OUTSIDE_TIME policy remains
  // unchanged and takes precedence before extrapolation is attempted.
  swcme::kinematics::Config data;
  data.mode=swcme::kinematics::Mode::DataDriven;
  data.r0_m=20.0*RS;
  data.V0_m_s=1.0e6;
  data.Vsw_m_s=4.0e5;
  data.Gamma_m_inv=1.0e-10;
  const double minimum=swcme::solarwind::MIN_RADIUS_M;
  data.data_time_s={0.0,100.0};
  data.data_radius_m={2.0*minimum,3.0*minimum};
  data.extrapolation=swcme::kinematics::ExtrapolationPolicy::Ballistic;
  const double endpoint_speed=minimum/100.0;
  for(double t : {-100.0,-50.0,-10.0,-1.0,101.0,110.0,200.0,1100.0}) {
    const bool before=t<0.0;
    const double anchor_time=before ? 0.0 : 100.0;
    const double anchor_radius=before ? 2.0*minimum : 3.0*minimum;
    const long double expected=
        static_cast<long double>(anchor_radius)+
        static_cast<long double>(endpoint_speed)*(t-anchor_time);
    const swcme::kinematics::State state=
        swcme::kinematics::evaluate(data,t);
    context.expect_true(state.status==Status::Ok,
                        "valid PCHIP ballistic extension has OK status");
    check_rel(context,"PCHIP independent extrapolated radius",state.radius_m,
              static_cast<double>(expected),4.0e-15,1.0);
    check_rel(context,"PCHIP independent endpoint speed",state.speed_m_s,
              endpoint_speed,2.0e-15,1.0);
  }
  context.expect_true(
      swcme::kinematics::evaluate(data,-101.0).status==Status::OutsideDomain,
      "pre-PCHIP sub-domain radius returns OUTSIDE_DOMAIN");
  context.expect_true(
      swcme::kinematics::evaluate(
          data,std::numeric_limits<double>::max()).status==
          Status::OutsideDomain,
      "post-PCHIP overflow returns OUTSIDE_DOMAIN");

  data.extrapolation=swcme::kinematics::ExtrapolationPolicy::OutsideTime;
  context.expect_true(
      swcme::kinematics::evaluate(data,-101.0).status==Status::OutsideTime &&
      swcme::kinematics::evaluate(data,101.0).status==Status::OutsideTime,
      "OUTSIDE_TIME policy still refuses both extensions before domain math");

  // Flat PCHIP endpoints are intentionally supported stationary fronts.  Even
  // a very large finite time cannot overflow a zero-slope continuation, which
  // distinguishes valid zero speed from the unsupported negative-speed branch.
  data.data_radius_m={2.0*minimum,2.0*minimum};
  data.extrapolation=swcme::kinematics::ExtrapolationPolicy::Ballistic;
  const swcme::kinematics::State stationary=
      swcme::kinematics::evaluate(data,std::numeric_limits<double>::max());
  context.expect_true(stationary.status==Status::Ok &&
                          stationary.radius_m==2.0*minimum &&
                          stationary.speed_m_s==0.0,
                      "stationary PCHIP continuation remains in domain");

  // Finally verify status propagation through both public model wrappers.  A
  // valid configured table whose explicit backward extension crosses the
  // radial boundary must make prepare_step() fail with OUTSIDE_DOMAIN rather
  // than exposing an invalid StepState or a generic configuration diagnosis.
  swcme1d::Params one;
  swcme3d::Params three;
  one.kinematics_mode=swcme::kinematics::Mode::DataDriven;
  one.data_time_s={0.0,100.0};
  one.data_radius_Rs={2.0,3.0};
  one.data_extrapolation=swcme::kinematics::ExtrapolationPolicy::Ballistic;
  three.kinematics_mode=one.kinematics_mode;
  three.data_time_s=one.data_time_s;
  three.data_radius_Rs=one.data_radius_Rs;
  three.data_extrapolation=one.data_extrapolation;

  // A modest backward extension remains inside the supported radial domain
  // and must reach both dimensional physics paths, proving the wrappers no
  // longer reject all negative time before consulting common kinematics.
  const swcme1d::StepState valid_backward_one=
      swcme1d::Model(one).prepare_step(-10.0);
  const swcme3d::StepState valid_backward_three=
      swcme3d::Model(three).prepare_step(-10.0);
  context.expect_true(valid_backward_one.r_sh_m>
                          swcme::solarwind::MIN_RADIUS_M &&
                          valid_backward_three.r_sh_m==
                          valid_backward_one.r_sh_m,
                      "1-D/3-D wrappers accept the same valid backward state");
  bool one_domain=false;
  bool three_domain=false;
  try {
    (void)swcme1d::Model(one).prepare_step(-200.0);
  } catch (const std::runtime_error& error) {
    one_domain=std::string(error.what()).find("OUTSIDE_DOMAIN")!=
               std::string::npos;
  }
  try {
    (void)swcme3d::Model(three).prepare_step(-200.0);
  } catch (const std::runtime_error& error) {
    three_domain=std::string(error.what()).find("OUTSIDE_DOMAIN")!=
                 std::string::npos;
  }
  context.expect_true(one_domain && three_domain,
                      "1-D/3-D preparation propagates OUTSIDE_DOMAIN");
}
