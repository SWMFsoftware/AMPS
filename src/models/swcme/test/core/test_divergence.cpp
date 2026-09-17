#include "test_framework.hpp"

#include <swcme1d.hpp>
#include <swcme3d.hpp>
#include <swcme_divergence.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <vector>

namespace {

bool relative_close(double a,double b,double rel) {
  const double scale=std::max({std::abs(a),std::abs(b),1.0e-300});
  return std::isfinite(a) && std::isfinite(b) && std::abs(a-b)<=rel*scale;
}

void source_shock_only(swcme1d::Params& p) {
  p.region_mode=swcme::regions::Mode::ShockOnly;
  p.shock_acceleration_mode=swcme::acceleration::Mode::Source;
}

void source_shock_only(swcme3d::Params& p) {
  p.region_mode=swcme::regions::Mode::ShockOnly;
  p.shock_acceleration_mode=swcme::acceleration::Mode::Source;
}

} // namespace

void test_div01(swcme_test::Context& context) {
  std::cout << "DIV01 analytical divergence of constant radial solar wind\n";

  swcme1d::Params p1;
  swcme3d::Params p3;
  source_shock_only(p1);
  source_shock_only(p3);
  p1.V_sw_kms=437.0;
  p3.V_sw_kms=p1.V_sw_kms;
  const swcme1d::Model m1(p1);
  const swcme3d::Model m3(p3);
  const auto s1=m1.prepare_step(3600.0);
  const auto s3=m3.prepare_step(3600.0);

  const std::array<double,4> radii{{
      2.0*swcme::constants::SOLAR_RADIUS_M,
      0.1*swcme::constants::AU_M,
      0.55*swcme::constants::AU_M,
      swcme::constants::AU_M}};
  const std::array<std::array<double,3>,4> dirs{{
      {{1.0,0.0,0.0}},
      {{0.0,1.0,0.0}},
      {{0.36,-0.48,0.80}},
      {{-0.6,0.0,0.8}}}};

  for (std::size_t ir=0;ir<radii.size();++ir) {
    const double r=radii[ir];
    const double expected=2.0*s1.V_up_ms/r;

    double n=0.0,V=0.0,Br=0.0,Bphi=0.0,Bmag=0.0,div1=0.0;
    const swcme::ModelStatus st1=m1.evaluate_radii_with_B_div_checked(
        s1,&r,&n,&V,&Br,&Bphi,&Bmag,&div1,1);
    context.expect_true(st1.ok(),"1-D analytical divergence query succeeds");
    context.expect_true(relative_close(div1,expected,2.0e-14),
                        "1-D SHOCK_ONLY div(V)=2Vsw/r");

    for (const auto& u : dirs) {
      const double x=r*u[0],y=r*u[1],z=r*u[2];
      double div3=0.0;
      const swcme::ModelStatus st3=m3.compute_divV_checked(
          s3,&x,&y,&z,&div3,1,1.0e-3);
      context.expect_true(st3.ok(),"3-D analytical divergence query succeeds");
      context.expect_true(relative_close(div3,expected,2.0e-14),
                          "3-D SHOCK_ONLY divergence is direction-independent");
    }

    std::cout << "  r/AU=" << std::scientific << r/swcme::constants::AU_M
              << " expected=" << expected << " 1D=" << div1 << '\n';
  }
}

void test_div02(swcme_test::Context& context) {
  std::cout << "DIV02 manufactured radial-flow divergence\n";

  // Manufactured positive radial profile.  x=r/r0 and
  //   Vr = V0 (1 + a x + b x^2)
  // gives an independent symbolic reference
  //   div(V) = V0/r0 [2(1+a x+b x^2)/x + a + 2 b x].
  // The production radial operator receives Vr and the analytical dVr/dr; it
  // must recover the symbolic expression to roundoff over the whole domain.
  const double r0=swcme::constants::AU_M;
  const double V0=410.0e3;
  const double a=0.15;
  const double b=0.03;
  for (double x : {0.03,0.1,0.25,0.5,0.8,1.2,2.0}) {
    const double r=x*r0;
    const double Vr=V0*(1.0+a*x+b*x*x);
    const double dVdr=(V0/r0)*(a+2.0*b*x);
    const auto terms=swcme::divergence::radial_terms(r,Vr,dVdr);
    const double expected=(V0/r0)*(
        2.0*(1.0+a*x+b*x*x)/x + a + 2.0*b*x);
    const double rel=std::abs(terms.total_s_inv-expected)/std::abs(expected);
    std::cout << "  x=" << x << " geom=" << terms.geometric_s_inv
              << " dVdr=" << terms.derivative_s_inv
              << " total=" << terms.total_s_inv
              << " rel_err=" << rel << '\n';
    context.expect_true(rel<1.0e-10,
                        "manufactured radial divergence agrees <1e-10 relative");
  }

  // Also verify that the analytical derivative carried by the common region
  // profile differentiates the actual 1-D transport velocity, rather than a
  // parallel approximation.  Mid-region points avoid mathematical interface
  // endpoints; a symmetric numerical derivative is used only as an independent
  // validation reference, not by production div(V).
  swcme1d::Params p;
  p.region_mode=swcme::regions::Mode::FullICME;
  p.shock_acceleration_mode=swcme::acceleration::Mode::ResolvedCompression;
  p.V0_sh_kms=1400.0;
  const swcme1d::Model m(p);
  const auto s=m.prepare_step(2.0*3600.0);
  const auto& bd=s.region_boundaries;
  const std::array<double,4> samples{{
      bd.R_sh_m,
      0.5*((bd.R_sh_m-0.5*bd.smooth_shock_width_m)+bd.R_le_m),
      bd.R_le_m,
      bd.R_te_m}};

  for (double r : samples) {
    const auto rv=swcme::regions::radial_velocity_state(
        r,bd,s.has_shock,s.V_up_ms,s.V2_shock_ms,s.V_LE_ms,
        s.region_config.V_ME_factor*s.V_up_ms,s.region_config.sheath_ramp_power);
    const double h=1.0e-6*r;
    const double rm=r-h,rp=r+h;
    double nm=0.0,Vm=0.0,np=0.0,Vp=0.0;
    context.expect_true(m.evaluate_radii_fast_checked(s,&rm,&nm,&Vm,1).ok(),
                        "left velocity sample succeeds");
    context.expect_true(m.evaluate_radii_fast_checked(s,&rp,&np,&Vp,1).ok(),
                        "right velocity sample succeeds");
    const double numeric=(Vp-Vm)/(2.0*h);
    const double scale=std::max({1.0e-15,std::abs(numeric),std::abs(rv.d_velocity_dr_s_inv)});
    const double rel=std::abs(numeric-rv.d_velocity_dr_s_inv)/scale;
    std::cout << "  region=" << swcme::regions::region_name(rv.region)
              << " analytic_dVdr=" << rv.d_velocity_dr_s_inv
              << " numeric=" << numeric << " scaled_err=" << rel << '\n';
    context.expect_true(rel<2.0e-4,
                        "common region derivative matches independent centered derivative");
  }
}

void test_div03(swcme_test::Context& context) {
  std::cout << "DIV03 general 3-D Cartesian divergence convergence\n";

  // A cubic manufactured vector field is chosen deliberately: a centered
  // second-order derivative is not exact for cubic terms, so the truncation
  // error must decrease ~h^2 under refinement instead of passing trivially.
  const std::array<double,3> x0{{0.7,-0.4,0.6}};
  auto field=[](const std::array<double,3>& x,std::array<double,3>& v) {
    const double X=x[0],Y=x[1],Z=x[2];
    v[0]=0.7*X*X*X + 0.2*X*Y - 0.1*Z;
    v[1]=-0.4*Y*Y*Y + 0.3*Y*Z + 0.2*X;
    v[2]=0.5*Z*Z*Z - 0.25*Z*X + 0.1*Y;
    return swcme::ModelStatus::success();
  };
  const double expected=2.1*x0[0]*x0[0] - 1.2*x0[1]*x0[1]
      +1.5*x0[2]*x0[2] + 0.2*x0[1] + 0.3*x0[2] - 0.25*x0[0];

  std::vector<double> errors;
  for (double h : {0.20,0.10,0.05,0.025}) {
    double div=0.0;
    const swcme::ModelStatus st=swcme::divergence::cartesian_second_order(
        x0,h,field,div);
    context.expect_true(st.ok(),"manufactured Cartesian divergence succeeds");
    const double error=std::abs(div-expected);
    errors.push_back(error);
    std::cout << "  h=" << h << " div=" << std::setprecision(14) << div
              << " exact=" << expected << " abs_err=" << error << '\n';
  }
  for (std::size_t i=1;i<errors.size();++i) {
    const double order=std::log(errors[i-1]/errors[i])/std::log(2.0);
    std::cout << "    observed order " << i << " = " << order << '\n';
    context.expect_true(order>1.90 && order<2.10,
                        "Cartesian divergence shows second-order convergence");
  }

  // Exercise the same production Cartesian operator through the 3-D model on
  // a constant-magnitude radial wind.  The canonical SHOCK_ONLY path is exact,
  // but this explicit numerical path must converge to that exact value when
  // requested for validation.
  swcme3d::Params p;
  source_shock_only(p);
  p.V_sw_kms=463.0;
  const swcme3d::Model m(p);
  const auto s=m.prepare_step(3600.0);
  const double r=0.73*swcme::constants::AU_M;
  const std::array<double,3> u{{0.36,-0.48,0.80}};
  const double x=r*u[0],y=r*u[1],z=r*u[2];
  const double exact=2.0*s.V_sw_ms/r;
  double previous=std::numeric_limits<double>::infinity();
  for (double frac : {2.0e-2,1.0e-2,5.0e-3,2.5e-3}) {
    double div=0.0;
    const auto st=m.compute_divV_cartesian_checked(s,&x,&y,&z,&div,1,frac);
    context.expect_true(st.ok(),"model Cartesian divergence succeeds");
    const double err=std::abs(div-exact);
    std::cout << "  model h/r=" << frac << " rel_err=" << err/std::abs(exact) << '\n';
    context.expect_true(err<previous,"model Cartesian divergence error decreases");
    previous=err;
  }
  context.expect_true(previous/std::abs(exact)<2.0e-5,
                      "refined model Cartesian divergence approaches analytical radial result");
}
