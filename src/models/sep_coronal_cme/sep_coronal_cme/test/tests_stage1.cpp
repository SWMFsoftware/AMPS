#include "test_framework.h"

#include "sep_coronal_cme/constants.h"
#include "sep_coronal_cme/pfss_harmonics.h"

#include <cmath>
#include <vector>

namespace SCCMTest { namespace {
using SEP::CoronalCME::FieldLineTopology;
using SEP::CoronalCME::HarmonicCoefficient;
using SEP::CoronalCME::MapSample;
using SEP::CoronalCME::MonopolePolicy;
using SEP::CoronalCME::PfssHarmonics;
using SEP::CoronalCME::Vec3;
constexpr double kPi=SEP::CoronalCME::Constants::kPi;

bool Near(double a,double b,double relative=1.0e-9,double absolute=1.0e-14) {
  return std::abs(a-b)<=absolute+relative*std::max(std::abs(a),std::abs(b));
}

PfssHarmonics Dipole() {
  const auto result=PfssHarmonics::Create(1.0,2.5,{{1,0,2.0e-4,0.0}});
  Require(result.ok(),result.status.message); return result.value;
}

void PFSS3D01() {
  const auto model=Dipole();
  const double theta=0.73,phi=0.4,r=1.4;
  const auto f=model.Evaluate(r,theta,phi); Require(f.ok(),f.status.message);
  // At the lower boundary Br must reproduce g_10 Y_10 exactly.
  const auto lower=model.Evaluate(1.0,theta,phi); Require(lower.ok(),lower.status.message);
  const double y10=std::sqrt(3.0/(4.0*kPi))*std::cos(theta);
  Require(Near(lower.value.brT,2.0e-4*y10),"axial dipole lower boundary mismatch");
  // Compare every analytic partial derivative against a centered independent
  // finite difference. The finite difference is only a test oracle.
  const double hR=2.0e-6,hA=2.0e-6;
  for(int coordinate=0;coordinate<3;++coordinate) {
    auto plus=model.Evaluate(r+(coordinate==0?hR:0.0),theta+(coordinate==1?hA:0.0),phi+(coordinate==2?hA:0.0));
    auto minus=model.Evaluate(r-(coordinate==0?hR:0.0),theta-(coordinate==1?hA:0.0),phi-(coordinate==2?hA:0.0));
    Require(plus.ok()&&minus.ok(),"finite-difference query failed");
    const double h=coordinate==0?hR:hA;
    const double numeric[3]={(plus.value.brT-minus.value.brT)/(2*h),
      (plus.value.bThetaT-minus.value.bThetaT)/(2*h),
      (plus.value.bPhiT-minus.value.bPhiT)/(2*h)};
    for(int component=0;component<3;++component)
      Require(Near(f.value.derivative[component][coordinate],numeric[component],2e-7,2e-12),
              "analytic PFSS derivative mismatch");
  }
}

void PFSS3D02() {
  const auto made=PfssHarmonics::Create(1.0,2.5,{{2,1,3e-5,-2e-5}});
  Require(made.ok(),made.status.message);
  const auto a=made.value.Evaluate(1.0,1.1,0.2);
  const auto b=made.value.Evaluate(1.0,1.1,0.9);
  Require(a.ok()&&b.ok(),"non-axisymmetric mode query failed");
  Require(!Near(a.value.brT,b.value.brT,1e-6),"m=1 mode lost longitude dependence");
  Require(std::abs(a.value.bPhiT)>1e-9,
          "non-axisymmetric mode lost its azimuthal field component");
}

void PFSS3D03() {
  Require(!PfssHarmonics::Create(1.0,2.5,{{0,0,1e-6,0}},MonopolePolicy::Reject,1e-12).ok(),
          "monopole rejection failed");
  const auto removed=PfssHarmonics::Create(1.0,2.5,
      {{0,0,1e-6,0},{1,0,2e-4,0}},MonopolePolicy::Remove,1e-12);
  Require(removed.ok()&&removed.value.Coefficients().size()==1,"monopole removal failed");
  double flux=0.0; const int nt=100,np=200;
  for(int i=0;i<nt;++i) for(int j=0;j<np;++j) {
    const double theta=(i+0.5)*kPi/nt,phi=(j+0.5)*2*kPi/np;
    const auto f=removed.value.Evaluate(1.0,theta,phi);
    flux+=f.value.brT*std::sin(theta)*(kPi/nt)*(2*kPi/np);
  }
  Require(std::abs(flux)<1e-17,"zero-net-flux quadrature failed");
}

void PFSS3D04() {
  const auto f=Dipole().Evaluate(2.5,0.8,1.2); Require(f.ok(),f.status.message);
  Require(std::abs(f.value.bThetaT)<1e-18&&std::abs(f.value.bPhiT)<1e-18,
          "PFSS source-surface field is not radial");
}

void PFSS3D05() {
  const auto model=Dipole();
  const auto polar=model.Classify({0.02,0.0,1.05},0.005,100000);
  const auto equator=model.Classify({1.05,0.0,0.02},0.005,100000);
  const auto polarFine=model.Classify({0.02,0.0,1.05},0.0025,200000);
  Require(polar.ok()&&equator.ok()&&polarFine.ok(),"topology trace failed");
  Require(polar.value==FieldLineTopology::OpenToOuterBoundary,
          "polar dipole line should be open to Rb");
  Require(equator.value==FieldLineTopology::ClosedBelowOuterBoundary,
          "equatorial dipole line should be closed");
  Require(polar.value==polarFine.value,"topology changed under step refinement");
}

std::vector<MapSample> DipoleMap(int nt,int np) {
  std::vector<MapSample> map;
  for(int i=0;i<nt;++i) for(int j=0;j<np;++j) {
    const double theta=(i+0.5)*kPi/nt,phi=(j+0.5)*2*kPi/np;
    map.push_back({theta,phi,2e-4*std::sqrt(3.0/(4*kPi))*std::cos(theta),
      std::sin(theta)*(kPi/nt)*(2*kPi/np)});
  }
  return map;
}

void PFSS3D06() {
  const auto map=DipoleMap(120,240);
  const auto coefficients=PfssHarmonics::ProjectMap(map,3,MonopolePolicy::Remove,1e-12);
  Require(coefficients.ok(),coefficients.status.message);
  double g10=0.0;
  for(const auto& c:coefficients.value) if(c.degree==1&&c.order==0) g10=c.cosineT;
  Require(Near(g10,2e-4,2e-4),"magnetogram harmonic round trip failed");
  const auto metrics=PfssHarmonics::CompareMap(map,coefficients.value);
  Require(metrics.ok()&&metrics.value.weightedRmsT<2e-8,"map reconstruction residual too large");
}

void PFSS3D07() {
  const auto made=PfssHarmonics::Create(1.0,2.5,
      {{1,0,1e-4,0},{6,2,4e-5,-3e-5}}); Require(made.ok(),made.status.message);
  const auto filtered=made.value.HeatKernelFiltered(3);
  const double expected=std::exp(-42.0/12.0);
  Require(Near(filtered.Coefficients()[1].cosineT,4e-5*expected),
          "heat-kernel transfer function mismatch");
  Require(made.value.Coefficients()[1].cosineT==4e-5,"filter mutated raw authority");
  const auto low=made.value.Evaluate(1.05,0.9,0.4);
  const auto smooth=filtered.Evaluate(1.05,0.9,0.4);
  Require(low.ok()&&smooth.ok()&&!Near(low.value.brT,smooth.value.brT,1e-6),
          "filter had no reconstruction effect");
}

void PFSS3D08() {
  const auto model=Dipole();
  for(int i=1;i<20;++i) {
    const double r=1.0+1.5*i/20.0,theta=0.2+2.7*i/20.0,phi=0.37*i;
    const auto authority=model.Evaluate(r,theta,phi);
    const auto cached=model.Evaluate(r,theta,phi); // direct-authority cache oracle
    Require(authority.ok()&&cached.ok(),"cache comparison query failed");
    Require(authority.value.brT==cached.value.brT&&
            authority.value.bThetaT==cached.value.bThetaT&&
            authority.value.bPhiT==cached.value.bPhiT,
            "direct harmonic evaluation is nondeterministic");
  }
}

void PFSS3D09() {
  const auto routed=SEP::CoronalCME::RouteTopology(
      FieldLineTopology::ClosedBelowOuterBoundary,
      FieldLineTopology::OpenToOuterBoundary,17,"composite-open-to-Ri");
  Require(routed.ok(),routed.status.message);
  Require(routed.value.pfssOpenToRb!=routed.value.compositeOpenToRi,
          "independent topology classifications collapsed");
  Require(routed.value.plasmaAuthority=="open-wind"&&routed.value.generation==17,
          "composite topology did not route plasma authority");
  Require(!SEP::CoronalCME::RouteTopology(FieldLineTopology::OpenToOuterBoundary,
      FieldLineTopology::OpenToOuterBoundary,1,"implicit-fallback").ok(),
      "unknown target-speed topology fell back");
}

}  // namespace

void RegisterStage1(Registry* tests) {
  (*tests)["PFSS3D01"]=PFSS3D01; (*tests)["PFSS3D02"]=PFSS3D02;
  (*tests)["PFSS3D03"]=PFSS3D03; (*tests)["PFSS3D04"]=PFSS3D04;
  (*tests)["PFSS3D05"]=PFSS3D05; (*tests)["PFSS3D06"]=PFSS3D06;
  (*tests)["PFSS3D07"]=PFSS3D07; (*tests)["PFSS3D08"]=PFSS3D08;
  (*tests)["PFSS3D09"]=PFSS3D09;
}

}  // namespace SCCMTest
