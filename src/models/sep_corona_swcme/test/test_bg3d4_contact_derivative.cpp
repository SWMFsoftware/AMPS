#include "contact_derivative_reference.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>

namespace {

using SEP::CoronalCME::Dot;
using SEP::CoronalCME::Norm;
using SEP::CoronalCME::Vec3;
using SEP::CoronaSwcme::TestReference::FourthOrderDifference;
using SEP::CoronaSwcme::TestReference::IndependentContactFluxKgM2S;

void Require(bool condition,const std::string& message) {
  if(!condition)throw std::runtime_error(message);
}

template<class Trajectory>
std::array<Vec3,4> Samples(const Trajectory& position,double h) {
  // Fixed-label positions only. No production velocity or material-map
  // implementation supplies an expected value in this independent fixture.
  return {position(-2*h),position(-h),position(h),position(2*h)};
}

template<class Operation>
void Reject(const Operation& operation,const std::string& message) {
  bool rejected=false;
  try {operation();}catch(const std::invalid_argument&){rejected=true;}
  Require(rejected,message);
}

void TestKnownFourthOrderError() {
  // An expanding spherical surface has X(theta,phi,t)=origin+R(t)n(theta,phi),
  // R(t)=2e6+2000t+100t^5 m. At this fixed chart n=(3/5,4/5,0), the exact
  // normal velocity at t=0 is 2000 m/s. For a t^5 term the frozen four-point
  // stencil error is exactly -4*a*h^4. Three preregistered levels therefore
  // have error ratio 16. The same geometry with a 1e9 m HCI translation checks
  // that a large irrelevant origin does not destroy the resolvable order.
  // These synthetic scales are numerical references, not a CME source model.
  const Vec3 normal{0.6,0.8,0};
  for(const Vec3 origin:{Vec3{0,0,0},Vec3{1e9,-2e9,3e9}}) {
    std::array<double,3> errors{};int level=0;
    for(double h:{0.4,0.2,0.1}) {
      const auto samples=Samples([&](double t) {
        return origin+(2e6+2000*t+100*std::pow(t,5))*normal;
      },h);
      const auto velocity=FourthOrderDifference(samples[0],samples[1],
          samples[2],samples[3],h);
      errors[level++]=Norm(velocity-2000*normal);
      const double scalar=FourthOrderDifference(100*std::pow(-2*h,5),
          100*std::pow(-h,5),100*std::pow(h,5),100*std::pow(2*h,5),h);
      Require(std::abs(scalar+400*std::pow(h,4))<1e-12,
          "scalar derivative lost the independent fifth-degree identity");
    }
    const double p01=std::log(errors[0]/errors[1])/std::log(2.0);
    const double p12=std::log(errors[1]/errors[2])/std::log(2.0);
    Require(p01>3.95&&p01<4.05&&p12>3.95&&p12<4.05,
        "translated curved reference is not fourth-order convergent");
    std::cout<<"[CMBGU05-DERIVATIVE-ORDER] origin_m="<<origin.x
             <<" errors_m_s="<<errors[0]<<','<<errors[1]<<','<<errors[2]
             <<" orders="<<p01<<','<<p12<<'\n';
  }
}

void TestPhysicalFluxAndDefectDetection() {
  constexpr double rho=2e-14;
  // Identical dimensional-plus-relative physical bound to the frozen BG3D-4
  // contact fixture, independently written as its two documented terms.
  constexpr double limit=1e-20+1e-10*1e-7;
  const Vec3 normal{1,0,0};
  for(double h:{0.1,0.05,0.025}) {
    const auto samples=Samples([](double t) {
      return Vec3{1e9+16*t,-4e8+8*t,2e9-4*t};
    },h);
    // Tangential material velocity may differ from the chart velocity. Only
    // normal slip carries mass through this planar surface. The independently
    // differentiated positions, not a reported flux field, enter the grader.
    const Vec3 material{16,-200,3};
    const double good=IndependentContactFluxKgM2S(samples,h,rho,material,normal);
    Require(std::abs(good)<=limit,
        "impermeable translated surface was rejected");
    const double plus=IndependentContactFluxKgM2S(samples,h,rho,
        material+0.5*normal,normal);
    const double minus=IndependentContactFluxKgM2S(samples,h,rho,
        material-0.5*normal,normal);
    Require(plus>limit&&minus< -limit,
        "opposite physical normal-slip defects escaped the frozen flux gate");
    // Two patches each contribute 2 m^2 over 3 s. Signed leakage nearly
    // cancels, but the required integral of absolute leakage must reject both.
    const double absoluteKg=6*(std::abs(plus)+std::abs(minus));
    Require(std::abs(6*(plus+minus))<=12*limit&&absoluteKg>12*limit,
        "absolute area-time leakage lost its noncancellation requirement");
    const double badNormal=IndependentContactFluxKgM2S(samples,h,rho,
        material,Vec3{0,1,0});
    Require(std::abs(badNormal)>limit,
        "incorrect unit normal escaped independently differentiated flux");
    std::cout<<"[CMBGU05-DERIVATIVE-FLUX] h_s="<<h<<" good="<<good
             <<" positive_slip="<<plus<<" negative_slip="<<minus
             <<" absolute_leakage_kg="<<absoluteKg<<'\n';
  }
}

void TestInvalidMeasurements() {
  const auto samples=Samples([](double t){return Vec3{t,0,0};},0.1);
  const auto measure=[&](double h,double rho,Vec3 u,Vec3 n) {
    return IndependentContactFluxKgM2S(samples,h,rho,u,n);
  };
  const double nan=std::numeric_limits<double>::quiet_NaN();
  Reject([&]{measure(0,1,{1,0,0},{1,0,0});},"zero derivative step accepted");
  Reject([&]{measure(-0.1,1,{1,0,0},{1,0,0});},"negative step accepted");
  Reject([&]{measure(nan,1,{1,0,0},{1,0,0});},"nonfinite step accepted");
  Reject([&]{measure(0.1,0,{1,0,0},{1,0,0});},"zero density accepted");
  Reject([&]{measure(0.1,nan,{1,0,0},{1,0,0});},"nonfinite density accepted");
  Reject([&]{measure(0.1,1,{nan,0,0},{1,0,0});},"nonfinite velocity accepted");
  Reject([&]{measure(0.1,1,{1,0,0},{0,0,0});},"zero normal accepted");
  Reject([&]{measure(0.1,1,{1,0,0},{2,0,0});},"nonunit normal accepted");
  auto corrupt=samples;corrupt[0].x=nan;
  Reject([&]{IndependentContactFluxKgM2S(corrupt,0.1,1,{1,0,0},{1,0,0});},
      "nonfinite position stencil accepted");
  // Unsupported endpoints cannot be silently replaced with a zero position
  // or a reported API velocity. In the full model the same requirement is
  // enforced by checking every EvaluateContact status before using a stencil.
  bool supportRejected=false;
  try {
    Samples([](double t) {
      if(std::abs(t)>0.15)throw std::out_of_range("unsupported contact epoch");
      return Vec3{t,0,0};
    },0.1);
  }catch(const std::out_of_range&){supportRejected=true;}
  Require(supportRejected,"unsupported temporal stencil was filled silently");
}

} // namespace

int main() {
  try {
    TestKnownFourthOrderError();
    TestPhysicalFluxAndDefectDetection();
    TestInvalidMeasurements();
    std::cout<<"[CMBGU05-DERIVATIVE-CONTROLS] PASS independent derivative and slip rejection\n";
    return 0;
  }catch(const std::exception& error) {
    std::cerr<<"CMBGU05-DERIVATIVE-CONTROLS FAIL: "<<error.what()<<'\n';
    return 1;
  }
}
