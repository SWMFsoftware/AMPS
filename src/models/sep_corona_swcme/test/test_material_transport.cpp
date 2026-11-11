#include "sep_corona_swcme/material_transport.h"

#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>

namespace {
namespace M=SEP::CoronaSwcme::Material;
namespace XR=SEP::CoronaSwcme::CrossRay;
using SEP::CoronalCME::Cross;
using SEP::CoronalCME::Dot;
using SEP::CoronalCME::Norm;
using SEP::CoronalCME::Vec3;
const XR::Thermodynamics kGas{5.0/3.0,2};
const XR::Primitive kReference{2,{0,0,0},3,{1,-2,0.5}};

void Require(bool ok,const std::string& message) {
  if(!ok)throw std::runtime_error(message);
}

M::MapJet Affine(Vec3 label,double t) {
  const double a=1+0.1*t,b=1+0.05*t,c=1+0.03*t,s=0.2*t*t,h=-0.1*t;
  M::MapJet out;out.deformation={a,s,0,0,b,h,0,0,c};
  out.deformationRatePerS={0.1,0.4*t,0,0,0.05,-0.1,0,0,0.03};
  out.positionHciM={a*label.x+s*label.y+0.5*t*t,b*label.y+h*label.z-t,c*label.z+0.3*t};
  out.velocityHciMPerS={0.1*label.x+0.4*t*label.y+t,0.05*label.y-0.1*label.z-1,
      0.03*label.z+0.3};
  out.accelerationHciMPerS2={0.4*label.y+1,0,0};return out;
}

XR::Primitive AffineOracle(double t) {
  // Closed-form triangular-map solution, independent of production determinant,
  // inverse, rate or flux code. This is a manufactured physical reference.
  const double a=1+0.1*t,b=1+0.05*t,c=1+0.03*t,j=a*b*c;
  return {2/j,{},3*std::pow(j,-5.0/3.0),
      {(a-0.4*t*t)/j,(-2*b-0.05*t)/j,0.5*c/j}};
}

Vec3 EulerVelocity(Vec3 x,double t) {
  // Independent exact triangular back-substitution, rather than the library's
  // cofactor inverse: a_z, a_y, a_x are solved in that order.
  const double az=(x.z-0.3*t)/(1+0.03*t);
  const double ay=(x.y+t+0.1*t*az)/(1+0.05*t);
  const double ax=(x.x-0.5*t*t-0.2*t*t*ay)/(1+0.1*t);
  return {0.1*ax+0.4*t*ay+t,0.05*ay-0.1*az-1,0.03*az+0.3};
}

void Invariants() {
  for(double t:{0.0,0.2,1.0,3.0}) {
    const auto map=Affine({0.2,-0.3,0.4},t);
    const auto out=M::PushForwardAdiabatic(kReference,map,kGas);Require(out.ok(),out.status.message);
    const auto oracle=AffineOracle(t);const double j=(1+0.1*t)*(1+0.05*t)*(1+0.03*t);
    Require(std::abs(out.value.jacobian-j)<2e-14&&
        std::abs(out.value.primitive.densityKgM3-oracle.densityKgM3)<2e-14&&
        std::abs(out.value.primitive.pressurePa-oracle.pressurePa)<3e-14&&
        Norm(out.value.primitive.magneticFieldT-oracle.magneticFieldT)<3e-14,
        "affine mass, thermodynamics or tilted signed field does not match closed-form reference");
    const Vec3 columns[]={{1+0.1*t,0,0},{0.2*t*t,1+0.05*t,0},{0,-0.1*t,1+0.03*t}};
    const Vec3 area[]={Cross(columns[1],columns[2]),Cross(columns[2],columns[0]),
        Cross(columns[0],columns[1])};
    const double initialFlux[]={1,-2,0.5};
    for(int i=0;i<3;++i)for(double sign:{-1.0,1.0})Require(
        std::abs(Dot(out.value.primitive.magneticFieldT,sign*area[i])-sign*initialFlux[i])<3e-14,
        "actual deformed face loses signed magnetic flux");
    Require(std::abs(out.value.primitive.densityKgM3*j-2)<2e-14&&
        std::abs(out.value.primitive.pressurePa*std::pow(j,5.0/3.0)-3)<3e-14,
        "mass or adiabatic entropy is not conserved");
  }
}

void Rates() {
  for(double t:{0.2,1.0,3.0}) {
    const auto map=Affine({0.2,-0.3,0.4},t);
    const auto out=M::PushForwardAdiabatic(kReference,map,kGas);Require(out.ok(),out.status.message);
    const double h=0.001;
    const auto a=AffineOracle(t-h),b=AffineOracle(t+h),c=AffineOracle(t-2*h),d=AffineOracle(t+2*h);
    const double rhoRate=(8*(b.densityKgM3-a.densityKgM3)-(d.densityKgM3-c.densityKgM3))/(12*h);
    const double pRate=(8*(b.pressurePa-a.pressurePa)-(d.pressurePa-c.pressurePa))/(12*h);
    const Vec3 bRate=(8*(b.magneticFieldT-a.magneticFieldT)-(d.magneticFieldT-c.magneticFieldT))/(12*h);
    Require(std::abs(out.value.densityMaterialRateKgM3S-rhoRate)<2e-11&&
        std::abs(out.value.pressureMaterialRatePaPerS-pRate)<2e-11&&
        Norm(out.value.magneticMaterialRateTPerS-bRate)<2e-11,
        "material continuity/adiabatic/induction rates fail independent time differentiation");
    for(int axis=0;axis<3;++axis) {
      Vec3 p=map.positionHciM,m=p;
      if(axis==0){p.x+=h;m.x-=h;}else if(axis==1){p.y+=h;m.y-=h;}else {p.z+=h;m.z-=h;}
      const auto gradient=(EulerVelocity(p,t)-EulerVelocity(m,t))/(2*h);
      Require(std::abs(out.value.velocityGradientPerS[axis]-gradient.x)<2e-11&&
          std::abs(out.value.velocityGradientPerS[3+axis]-gradient.y)<2e-11&&
          std::abs(out.value.velocityGradientPerS[6+axis]-gradient.z)<2e-11,
          "material Fdot inverse(F) gradient fails independent Eulerian velocity field");
    }
  }
}

M::MapJet NonlinearMap(Vec3 a,double t) {
  // X=(a_x,a_y*exp(k*a_x*t),a_z), k=0.2/(m s). This map has spatially
  // varying volume and real lateral displacement. Its inverse is explicit.
  const double k=0.2,g=std::exp(k*a.x*t);
  M::MapJet out;out.positionHciM={a.x,a.y*g,a.z};
  out.velocityHciMPerS={0,k*a.x*a.y*g,0};
  out.accelerationHciMPerS2={0,k*k*a.x*a.x*a.y*g,0};
  out.deformation={1,0,0,k*t*a.y*g,g,0,0,0,1};
  out.deformationRatePerS={0,0,0,k*a.y*g*(1+k*a.x*t),k*a.x*g,0,0,0,0};return out;
}

Vec3 NonlinearField(Vec3 x,double t) {
  const Vec3 a{x.x,x.y*std::exp(-0.2*x.x*t),x.z};
  const auto out=M::PushForwardAdiabatic(kReference,NonlinearMap(a,t),kGas);
  Require(out.ok(),out.status.message);return out.value.primitive.magneticFieldT;
}

void Nonlinear() {
  constexpr double t=2,h=1e-4;double largest=0;
  long double mass=0;std::size_t cells=0;
  for(double x:{-0.4,-0.1,0.2,0.5})for(double y:{-0.3,0.0,0.3})for(double z:{-0.2,0.2}) {
    const auto map=NonlinearMap({x,y,z},t);
    const auto out=M::PushForwardAdiabatic(kReference,map,kGas);Require(out.ok(),out.status.message);
    const double j=std::exp(0.2*x*t);
    const Vec3 exact{std::exp(-0.2*x*t),0.2*t*y-2,0.5*std::exp(-0.2*x*t)};
    Require(Norm(out.value.primitive.magneticFieldT-exact)<2e-14&&
        std::abs(out.value.jacobian-j)<2e-14,"nonlinear Cauchy field/volume differs from analytic solution");
    mass+=out.value.primitive.densityKgM3*j;++cells;
    const auto position=map.positionHciM;
    const double bx1=NonlinearField(position+Vec3{h,0,0},t).x,
        bx0=NonlinearField(position-Vec3{h,0,0},t).x;
    const double by1=NonlinearField(position+Vec3{0,h,0},t).y,
        by0=NonlinearField(position-Vec3{0,h,0},t).y;
    const double bz1=NonlinearField(position+Vec3{0,0,h},t).z,
        bz0=NonlinearField(position-Vec3{0,0,h},t).z;
    largest=std::max(largest,std::abs((bx1-bx0+by1-by0+bz1-bz0)/(2*h)));
  }
  Require(std::abs(mass-2*cells)<2e-12&&largest<3e-10,
      "global nonlinear material mass or independently differentiated divergence does not close");
  std::cout<<"[CSWC0628-XRAY-MATERIAL-EVIDENCE] nonlinear_samples="<<cells<<
      " maximum_div_b_t_m="<<largest<<'\n';
}

void Rejection() {
  auto map=Affine({0,0,0},1);map.deformation={-1,0,0,0,1,0,0,0,1};
  Require(!M::PushForwardAdiabatic(kReference,map,kGas).ok(),"folded material map was regularized");
  map.deformation={1,1,0,1,1+std::numeric_limits<double>::epsilon(),0,0,0,1};
  Require(!M::PushForwardAdiabatic(kReference,map,kGas).ok(),"unresolved map determinant was accepted");
  map.deformation={1e200,0,0,0,1e200,0,0,0,1e200};
  Require(!M::PushForwardAdiabatic(kReference,map,kGas).ok(),"unrepresentable Jacobian was accepted");
  map=Affine({0,0,0},1);map.deformationRatePerS[0]=std::numeric_limits<double>::quiet_NaN();
  Require(!M::PushForwardAdiabatic(kReference,map,kGas).ok(),"nonfinite material rate was accepted");
  map=Affine({0,0,0},1);map.deformation={1e-12,0,0,0,1e-12,0,0,0,1e-12};
  map.deformationRatePerS={};
  const auto tiny=M::PushForwardAdiabatic(kReference,map,kGas);
  Require(tiny.ok()&&tiny.value.jacobian>0&&tiny.value.jacobian<1e-35,
      "resolved small positive Jacobian was replaced by a numerical floor");
  // This kernel check is not production minimum-J qualification. The frozen
  // physical minimum-J gate must still reject such extreme compression there.
  auto zero=kReference;zero.magneticFieldT={};
  const auto zeroField=M::PushForwardAdiabatic(zero,Affine({0,0,0},1),kGas);
  Require(zeroField.ok()&&Norm(zeroField.value.primitive.magneticFieldT)==0,
      "zero magnetic reference received a field floor");
}
} // namespace

int main() {
  try {
    Invariants();std::cout<<"[CSWC0628-XRAY-MATERIAL-INVARIANTS] PASS independent deformed-face mass/entropy/signed-flux references\n";
    Rates();std::cout<<"[CSWC0628-XRAY-MATERIAL-RATES] PASS independent material and Eulerian derivatives\n";
    Nonlinear();std::cout<<"[CSWC0628-XRAY-MATERIAL-NONLINEAR] PASS nonlinear lateral transport mass and divergence\n";
    Rejection();std::cout<<"[CSWC0628-XRAY-MATERIAL-REJECTION] PASS folded/unresolved maps reject without physical floors\n";
  }catch(const std::exception& e){std::cerr<<"material transport FAIL: "<<e.what()<<'\n';return 1;}
  return 0;
}
