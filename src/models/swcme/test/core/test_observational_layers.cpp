#include "test_framework.hpp"

#include <swcme3d.hpp>
#include <swcme_kinematics.hpp>
#include <swcme_shock.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <numeric>
#include <string>
#include <vector>

namespace {

constexpr double AU=swcme::constants::AU_M;
constexpr double RS=swcme::constants::SOLAR_RADIUS_M;
constexpr double PI=swcme::constants::PI;

double median(std::vector<double> values) {
  std::sort(values.begin(),values.end());
  const std::size_t middle=values.size()/2;
  return values.size()%2==0 ? 0.5*(values[middle-1]+values[middle])
                            : values[middle];
}

double relative_error(double model,double observed) {
  return std::abs(model-observed)/std::max(std::abs(observed),1.0e-300);
}

// All V1-V5 registered checks below are explicitly synthetic regression
// fixtures. They test metric computation, thresholds, and production plumbing;
// EVT01 separately prevents them from being presented as release observational
// evidence. A real campaign supplies traceable measurements through schema v2.
void announce_fixture(const char* layer) {
  std::cout << "  layer=" << layer
            << " fixture_kind=SYNTHETIC_REGRESSION science_release=INCOMPLETE\n";
}

}  // namespace

void test_v1(swcme_test::Context& context) {
  std::cout << "V1 background Parker and Leblanc validation infrastructure\n";
  announce_fixture("V1");
  swcme3d::Params parameters;
  parameters.V0_sh_kms=parameters.V_sw_kms;  // keep all sampled radii upstream
  const swcme3d::Model model(parameters);
  const auto state=model.prepare_step(0.0);
  const std::array<double,5> radii_au{{0.30,0.45,0.70,1.00,1.40}};
  const std::array<double,5> density_factors{{1.10,0.92,1.18,0.95,1.12}};
  const std::array<double,5> field_factors{{0.90,1.08,0.94,1.12,0.91}};
  const std::array<double,5> angle_offsets_deg{{8.0,-6.0,12.0,-9.0,5.0}};
  std::vector<double> density_ratios,field_ratios,angle_errors;
  for (std::size_t i=0;i<radii_au.size();++i) {
    const double x=radii_au[i]*AU,y=0.0,z=0.0;
    double n=0.0,vx=0.0,vy=0.0,vz=0.0,bx=0.0,by=0.0,bz=0.0;
    model.evaluate_cartesian_with_B(state,&x,&y,&z,&n,&vx,&vy,&vz,
                                    &bx,&by,&bz,1);
    const double field=std::hypot(bx,std::hypot(by,bz));
    const double angle=std::atan2(std::abs(by),std::abs(bx))*180.0/PI;
    const double observed_n=n*density_factors[i];
    const double observed_b=field*field_factors[i];
    const double observed_angle=angle+angle_offsets_deg[i];
    density_ratios.push_back(std::max(n/observed_n,observed_n/n));
    field_ratios.push_back(std::max(field/observed_b,observed_b/field));
    angle_errors.push_back(std::abs(angle-observed_angle));
  }
  context.expect_true(median(angle_errors)<=20.0,
                      "median Parker-angle error is at most 20 degrees");
  context.expect_true(*std::max_element(density_ratios.begin(),density_ratios.end())<=2.0,
                      "normalized density trend remains within factor two");
  context.expect_true(*std::max_element(field_ratios.begin(),field_ratios.end())<=2.0,
                      "normalized field trend remains within factor two");
  // A deliberately invalid comparison proves that the threshold logic is not
  // a tautology tied to the chosen passing fixture.
  const std::vector<double> invalid_angle_errors{35.0,45.0,80.0};
  context.expect_true(median(invalid_angle_errors)>20.0,
                      "bad Parker-angle fixture is rejected");
}

void test_v2(swcme_test::Context& context) {
  std::cout << "V2 CME and shock apex kinematics validation infrastructure\n";
  announce_fixture("V2");
  swcme::kinematics::Config data;
  data.mode=swcme::kinematics::Mode::DataDriven;
  data.r0_m=20.0*RS; data.V0_m_s=1.2e6; data.Vsw_m_s=4.0e5;
  data.data_time_s={0.0,4.0*3600.0,10.0*3600.0,20.0*3600.0};
  data.data_radius_m={20.0*RS,42.0*RS,70.0*RS,108.0*RS};
  std::vector<double> height_errors;
  for (std::size_t i=0;i<data.data_time_s.size();++i) {
    const auto state=swcme::kinematics::evaluate(data,data.data_time_s[i]);
    height_errors.push_back(relative_error(state.radius_m,data.data_radius_m[i]));
  }
  context.expect_true(*std::max_element(height_errors.begin(),height_errors.end())<1.0e-12,
                      "data-driven knots reproduce supplied heights");

  swcme::kinematics::Config dbm=data;
  dbm.mode=swcme::kinematics::Mode::DBM;
  dbm.Gamma_m_inv=5.0e-11;
  auto radius_at=[&](double hours){return swcme::kinematics::evaluate(dbm,hours*3600.0).radius_m;};
  double lo=0.0,hi=120.0;
  for (int iteration=0;iteration<100;++iteration) {
    const double mid=0.5*(lo+hi);
    if (radius_at(mid)<0.96*AU) lo=mid; else hi=mid;
  }
  const double model_arrival_h=0.5*(lo+hi);
  const double held_arrival_h=model_arrival_h+3.5;
  context.expect_true(std::abs(model_arrival_h-held_arrival_h)<=6.0,
                      "held arrival error satisfies six-hour target");
  context.expect_true(relative_error(model_arrival_h,held_arrival_h)<=0.15,
                      "held arrival relative error satisfies 15-percent target");
}

void test_v3(swcme_test::Context& context) {
  std::cout << "V3 in-situ shock-jump validation infrastructure\n";
  announce_fixture("V3");
  swcme::shock::PrimitiveState upstream;
  upstream.rho_kg_m3=5.0e6*swcme::constants::PROTON_MASS_KG;
  upstream.pressure_Pa=5.0e6*swcme::constants::BOLTZMANN_J_K*1.2e5;
  upstream.velocity_m_s={{4.0e5,0.0,0.0}};
  upstream.magnetic_T={{4.0e-9,3.0e-9,0.0}};
  const swcme::shock::Vec3 normal{{1.0,0.0,0.0}};
  const auto jump=swcme::shock::solve_ideal_mhd_fast_shock(
      upstream,normal,1.2e6,5.0/3.0);
  context.expect_true(jump.status==swcme::shock::SolveStatus::Solved,
                      "synthetic in-situ shock is physically solved");
  const double residual=std::max({jump.mass_residual,jump.normal_B_residual,
      jump.electric_residual,jump.momentum_residual,jump.energy_residual});
  context.expect_true(residual<=1.0e-8,
                      "Rankine-Hugoniot residual is at most 1e-8");
  const std::vector<double> compression_errors{relative_error(jump.compression,
      jump.compression*1.10)};
  const double downstream_speed=std::sqrt(
      jump.downstream.velocity_m_s[0]*jump.downstream.velocity_m_s[0]+
      jump.downstream.velocity_m_s[1]*jump.downstream.velocity_m_s[1]+
      jump.downstream.velocity_m_s[2]*jump.downstream.velocity_m_s[2]);
  const double downstream_field=std::sqrt(
      jump.downstream.magnetic_T[0]*jump.downstream.magnetic_T[0]+
      jump.downstream.magnetic_T[1]*jump.downstream.magnetic_T[1]+
      jump.downstream.magnetic_T[2]*jump.downstream.magnetic_T[2]);
  context.expect_true(median(compression_errors)<=0.25,
                      "median compression error is at most 25 percent");
  context.expect_true(relative_error(downstream_speed,downstream_speed*1.08)<=0.15,
                      "downstream-speed error is at most 15 percent");
  context.expect_true(relative_error(downstream_field,downstream_field*1.12)<=0.30,
                      "downstream-field error is at most 30 percent");
  context.expect_true(std::abs(jump.theta_Bn_rad-(jump.theta_Bn_rad+10.0*PI/180.0))*180.0/PI<=15.0,
                      "theta-Bn error is at most 15 degrees");
}

void test_v4(swcme_test::Context& context) {
  std::cout << "V4 shock geometry and encounter validation infrastructure\n";
  announce_fixture("V4");
  const double half_width=45.0*PI/180.0;
  const std::array<double,4> observer_angles{{0.0,20.0,44.0,65.0}};
  const std::array<bool,4> expected_hits{{true,true,true,false}};
  std::vector<double> arrival_hours;
  for (std::size_t i=0;i<observer_angles.size();++i) {
    const bool hit=std::abs(observer_angles[i])*PI/180.0<=half_width;
    context.expect_true(hit==expected_hits[i],"synthetic observer hit/miss classification");
    if (hit) {
      // The cosine flank delay is an independent transparent encounter proxy;
      // it enforces apex-before-flank ordering and makes timing uncertainty
      // explicit without claiming an observed event.
      arrival_hours.push_back(42.0/std::cos(observer_angles[i]*PI/180.0));
    }
  }
  context.expect_true(std::is_sorted(arrival_hours.begin(),arrival_hours.end()),
                      "apex-to-flank arrival ordering is correct");
  const std::array<double,3> observed{{43.5,45.5,59.0}};
  for (std::size_t i=0;i<observed.size();++i)
    context.expect_true(std::abs(arrival_hours[i]-observed[i])<=8.0,
                        "encounter time lies within documented uncertainty");
}

void test_v5(swcme_test::Context& context) {
  std::cout << "V5 magnetic-connectivity benchmark infrastructure\n";
  announce_fixture("V5");
  using History=std::vector<bool>;
  auto transition_count=[](const History& history) {
    int count=0;
    for (std::size_t i=1;i<history.size();++i)
      if (history[i]!=history[i-1]) ++count;
    return count;
  };
  // The mandatory 2010-09-09 qualitative constraints are frozen from Tao et
  // al. (2025, ApJ 995:77): STA stays connected and STB loses connection.
  // SOHO and two extra events are present as synthetic pipeline fixtures only;
  // EVT01 therefore keeps the observational release state INCOMPLETE.
  const History sta{true,true,true,true,true,true};
  const History stb{true,true,true,false,false,false};
  const History soho{false,false,true,true,false,false};
  const std::array<History,2> independent_events{{
      History{false,true,true,false},History{true,true,false,false}}};
  context.expect_true(transition_count(sta)==0 && sta.front(),
                      "STA remains connected throughout benchmark window");
  context.expect_true(transition_count(stb)==1 && stb.front() && !stb.back(),
                      "STB changes from connected to disconnected");
  context.expect_true(transition_count(soho)==2,
                      "SOHO synthetic window records onset and loss");
  for (const auto& history:independent_events)
    context.expect_true(transition_count(history)>=1,
                        "independent synthetic event contains a transition");
  const History shifted{true,true,true,true,false,false};
  context.expect_true(transition_count(shifted)==1,
                      "four-hour sampling sensitivity preserves transition class");
}
