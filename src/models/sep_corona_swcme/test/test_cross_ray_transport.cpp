#include "sep_corona_swcme/cross_ray_transport.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {

namespace XR=SEP::CoronaSwcme::CrossRay;
using SEP::CoronalCME::Vec3;
constexpr double kMu0=1.25663706212e-6;
constexpr double kPi=3.141592653589793238462643383279502884;
const XR::Thermodynamics kGas{5.0/3,kMu0};

void Require(bool condition,const std::string& message) {
  if(!condition)throw std::runtime_error(message);
}

bool Close(double a,double b,double relative=2e-13) {
  return std::abs(a-b)<=relative*std::max({std::abs(a),std::abs(b),1.0});
}

std::array<double,8> Values(const XR::ConservedDensity& u) {
  return {u.massDensityKgM3,u.momentumDensityKgM2S.x,u.momentumDensityKgM2S.y,
      u.momentumDensityKgM2S.z,u.totalEnergyDensityJPerM3,
      u.magneticFieldT.x,u.magneticFieldT.y,u.magneticFieldT.z};
}

std::array<double,8> Values(const XR::Flux& f) {
  return {f.massKgM2S,f.momentumPa.x,f.momentumPa.y,f.momentumPa.z,
      f.energyWPerM2,f.inductionTMPerS.x,f.inductionTMPerS.y,f.inductionTMPerS.z};
}

XR::Primitive ReferenceState() {
  // A self-consistent SI reference with magnetic pressure 7 Pa, thermal
  // pressure 5 Pa, kinetic energy 14 J/m^3 and thermal energy 7.5 J/m^3.
  // The resulting hand-calculated fluxes below do not call a second version
  // of the production flux formula to generate their expected values.
  const double root=std::sqrt(kMu0);
  return {2,{2,3,-1},5,{root,2*root,3*root}};
}

void CheckFlux(const XR::Flux& got,const std::array<double,8>& expected) {
  const auto values=Values(got);
  for(std::size_t i=0;i<values.size();++i)
    Require(Close(values[i],expected[i]),"independent SI flux component mismatch");
}

void TestFluxAndWaves() {
  const auto p=ReferenceState();const double root=std::sqrt(kMu0);
  const auto u=XR::ToConserved(p,kGas);
  Require(u.ok()&&Close(u.value.totalEnergyDensityJPerM3,28.5),
      "independent conservative total energy mismatch");
  const auto back=XR::ToPrimitive(u.value,kGas);
  Require(back.ok()&&Close(back.value.pressurePa,5)&&
      Close(back.value.velocityMPerS.y,3),"primitive conversion lost pressure/velocity");
  const auto face=XR::HllFlux(p,p,{{1,0,0},0},kGas);
  Require(face.ok(),face.status.message);
  CheckFlux(face.value.value,{4,19,10,-7,76,0,root,7*root});

  auto gas=p;gas.magneticFieldT={};
  const auto gasFace=XR::HllFlux(gas,gas,{{1,0,0},0},kGas);
  Require(gasFace.ok(),gasFace.status.message);
  CheckFlux(gasFace.value.value,{4,13,12,-4,53,0,0,0});
  // Distinct magnetic states exercise the central HLL branch as well as the
  // identical-state physical-flux limit. Here rho=p=1, v=0, Bn=sqrt(mu0),
  // Bt=+/-sqrt(mu0): E and normal stress are equal, opposite tangential
  // stresses average to zero, and the HLL induction flux is cf*sqrt(mu0).
  // The independently reduced characteristic is cf^2=(11+sqrt(61))/6.
  const XR::Primitive ml{1,{},1,{root,root,0}},mr{1,{},1,{root,-root,0}};
  const auto magneticJump=XR::HllFlux(ml,mr,{{1,0,0},0},kGas);
  Require(magneticJump.ok(),magneticJump.status.message);
  CheckFlux(magneticJump.value.value,{0,1,0,0,0,0,
      root*std::sqrt((11+std::sqrt(61.0))/6),0});
  const auto acoustic=XR::FastSpeedMPerS(gas,{1,0,0},kGas);
  Require(acoustic.ok()&&Close(acoustic.value,std::sqrt(25.0/6)),
      "zero-field fast characteristic did not reduce to sound speed");
  gas.magneticFieldT={0,2*root,0};
  const auto perpendicular=XR::FastSpeedMPerS(gas,{1,0,0},kGas);
  Require(perpendicular.ok()&&Close(perpendicular.value,std::sqrt(37.0/6)),
      "perpendicular canonical fast characteristic mismatch");
  gas.magneticFieldT={4*root,0,0};
  const auto parallel=XR::FastSpeedMPerS(gas,{1,0,0},kGas);
  Require(parallel.ok()&&Close(parallel.value,std::sqrt(8.0)),
      "parallel characteristic lost the normal magnetic component");
  gas.magneticFieldT={std::sqrt(kGas.gammaAdiabatic*gas.pressurePa*kMu0),0,0};
  const auto degenerate=XR::FastSpeedMPerS(gas,{1,0,0},kGas);
  Require(degenerate.ok()&&Close(degenerate.value,std::sqrt(25.0/6)),
      "parallel acoustic/Alfven degeneracy produced an invalid discriminant");

  XR::Primitive left{1,{10,0,0},1,{}},right{2,{11,0,0},2,{}};
  const auto plus=XR::HllFlux(left,right,{{1,0,0},0},kGas);
  Require(plus.ok()&&plus.value.leftWaveMPerS==0,"positive upwind branch not selected");
  CheckFlux(plus.value.value,{10,101,0,0,525,0,0,0});
  left.velocityMPerS.x=-10;right.velocityMPerS.x=-11;
  const auto minus=XR::HllFlux(left,right,{{1,0,0},0},kGas);
  Require(minus.ok()&&minus.value.rightWaveMPerS==0,"negative upwind branch not selected");
  CheckFlux(minus.value.value,{-22,244,0,0,-1386,0,0,0});
}

void TestAleAndOrientation() {
  const auto p=ReferenceState();const double root=std::sqrt(kMu0);
  const auto moving=XR::HllFlux(p,p,{{1,0,0},1.5},kGas);
  Require(moving.ok(),moving.status.message);
  // Full Reynolds F-wn*U values, including magnetic and kinetic energy. The
  // normal induction entry is -wn*Bn, not the stationary-face zero value.
  CheckFlux(moving.value.value,{1,13,1,-4,33.25,-1.5*root,-2*root,2.5*root});
  const auto reversed=XR::HllFlux(p,p,{{-1,0,0},-1.5},kGas);
  Require(reversed.ok(),reversed.status.message);
  const auto f=Values(moving.value.value),g=Values(reversed.value.value);
  for(std::size_t i=0;i<f.size();++i)
    Require(Close(f[i],-g[i]),"shared-face orientation reversal is not antisymmetric");
  auto rotated=p;
  rotated.velocityMPerS={-3,2,-1};
  rotated.magneticFieldT={-2*root,root,3*root};
  const auto rotation=XR::HllFlux(rotated,rotated,{{0,1,0},0},kGas);
  Require(rotation.ok(),rotation.status.message);
  CheckFlux(rotation.value.value,{4,-10,19,-7,76,-root,0,7*root});
}

void TestPairConservation() {
  const auto state=XR::ToConserved(ReferenceState(),kGas);
  Require(state.ok(),state.status.message);
  XR::Cell left{2,state.value},right{3,state.value};
  const auto flux=XR::HllFlux(ReferenceState(),ReferenceState(),{{1,0,0},0},kGas);
  Require(flux.ok(),flux.status.message);
  const auto exchanged=XR::TransferStationaryFace(left,right,flux.value,0.5,0.001,kGas);
  Require(exchanged.ok(),exchanged.status.message);
  const auto l=Values(exchanged.value.left.density),r=Values(exchanged.value.right.density);
  const auto original=Values(state.value);
  for(std::size_t i=0;i<l.size();++i)
    Require(Close(2*l[i]+3*r[i],5*original[i]),
        "shared transfer lost mass/vector momentum/energy/induction inventory");
  Require(Close(exchanged.value.transfer.massKg,0.002)&&
      Close(exchanged.value.transfer.momentumKgMPerS.y,0.005)&&
      Close(exchanged.value.transfer.totalEnergyJ,0.038),
      "shared-face extensive ledger has incorrect area/time/SI normalization");
  Require(Values(left.density)==original&&Values(right.density)==original,
      "private transfer mutated its committed inputs");
}

void TestRejectionAndRollback() {
  auto p=ReferenceState();const auto state=XR::ToConserved(p,kGas);
  Require(state.ok(),state.status.message);
  const XR::Cell left{2,state.value},right{3,state.value};
  const auto original=Values(state.value);
  const auto face=XR::HllFlux(p,p,{{1,0,0},0},kGas);
  Require(face.ok(),face.status.message);
  const auto rejected=XR::TransferStationaryFace(left,right,face.value,1,1000,kGas);
  Require(!rejected.ok()&&rejected.status.code==SEP::Core::StatusCode::NumericalFailure,
      "unqualified timestep did not reject nonpositive candidate density/pressure");
  Require(Values(left.density)==original&&Values(right.density)==original,
      "rejected transfer changed a committed donor/recipient");
  const auto moving=XR::HllFlux(p,p,{{1,0,0},1},kGas);
  Require(moving.ok()&&!XR::TransferStationaryFace(
      left,right,moving.value,1,0.01,kGas).ok(),
      "moving-face flux was accepted without geometric-conservation authority");
  p.pressurePa=-1;
  Require(!XR::ToConserved(p,kGas).ok(),"negative pressure was floored");
  p=ReferenceState();p.densityKgM3=0;
  Require(!XR::ToConserved(p,kGas).ok(),"zero density was floored");
  p=ReferenceState();p.magneticFieldT.x+=1;
  const auto monopole=XR::HllFlux(ReferenceState(),p,{{1,0,0},0},kGas);
  Require(!monopole.ok()&&monopole.status.code==SEP::Core::StatusCode::DataIntegrityFailure,
      "inconsistent normal magnetic field was silently projected/cleaned");
  Require(!XR::HllFlux(p,p,{{0,0,0},0},kGas).ok(),"zero face normal was accepted");
  Require(!XR::HllFlux(p,p,{{2,0,0},0},kGas).ok(),"nonunit face normal was accepted");
  auto badFace=face.value;badFace.value.energyWPerM2=
      std::numeric_limits<double>::quiet_NaN();
  Require(!XR::TransferStationaryFace(left,right,badFace,1,0.01,kGas).ok(),
      "nonfinite shared energy flux was accepted");
  Require(!XR::TransferStationaryFace(left,right,face.value,-1,0.01,kGas).ok(),
      "negative physical face area was accepted");
  Require(!XR::TransferStationaryFace(left,right,face.value,1,-0.01,kGas).ok(),
      "negative transport time step was accepted");
  auto corrupt=state.value;corrupt.totalEnergyDensityJPerM3=0;
  Require(!XR::ToPrimitive(corrupt,kGas).ok(),"unresolved thermal energy was floored");
  auto gas=kGas;gas.gammaAdiabatic=1;
  Require(!XR::ToConserved(ReferenceState(),gas).ok(),"gamma=1 ideal-gas closure accepted");
  p=ReferenceState();p.velocityMPerS.x=std::numeric_limits<double>::quiet_NaN();
  Require(!XR::ToConserved(p,kGas).ok(),"nonfinite primitive state was accepted");
}

double AdvectContact(int count) {
  // Independent periodic density contact: rho=1+0.2 sin(2pi(y-vt)), p=1 Pa,
  // v_y=0.4 m/s, B=0. Exact CELL AVERAGES include sinc(pi*dy). The normal
  // direction is lateral y, not the legacy radial axis; all conserved fields
  // evolve through the production HLL/shared-transfer routines.
  constexpr double speed=0.4,end=0.25,cfl=0.2;
  const double width=1.0/count;
  const double sinc=std::sin(kPi*width)/(kPi*width);
  std::vector<XR::Cell> cells;
  std::array<double,8> initial{};
  for(int i=0;i<count;++i) {
    const double rho=1+0.2*sinc*std::sin(2*kPi*(i+0.5)*width);
    const auto value=XR::ToConserved({rho,{0,speed,0},1,{}},kGas);
    Require(value.ok(),value.status.message);
    cells.push_back({width,value.value});
    const auto data=Values(value.value);
    for(std::size_t j=0;j<data.size();++j)initial[j]+=width*data[j];
  }
  double time=0;int steps=0;
  while(time<end) {
    std::vector<XR::FaceFlux> fluxes;double maximum=0;
    for(int i=0;i<count;++i) {
      const auto l=XR::ToPrimitive(cells[i].density,kGas),
          r=XR::ToPrimitive(cells[(i+1)%count].density,kGas);
      Require(l.ok()&&r.ok(),"periodic contact lost positive primitive state");
      const auto f=XR::HllFlux(l.value,r.value,{{0,1,0},0},kGas);
      Require(f.ok(),f.status.message);
      fluxes.push_back(f.value);maximum=std::max(maximum,f.value.maximumRelativeSignalMPerS);
    }
    const double dt=std::min(end-time,cfl*width/maximum);
    for(int i=0;i<count;++i) {
      const int next=(i+1)%count;
      const auto exchanged=XR::TransferStationaryFace(
          cells[i],cells[next],fluxes[i],1,dt,kGas);
      Require(exchanged.ok(),exchanged.status.message);
      cells[i]=exchanged.value.left;cells[next]=exchanged.value.right;
    }
    time+=dt;++steps;
  }
  double error=0;std::array<double,8> final{};
  for(int i=0;i<count;++i) {
    const double exact=1+0.2*sinc*std::sin(2*kPi*((i+0.5)*width-speed*end));
    error+=width*std::abs(cells[i].density.massDensityKgM3-exact);
    const auto data=Values(cells[i].density);
    for(std::size_t j=0;j<data.size();++j)final[j]+=width*data[j];
  }
  // Bound roundoff by the explicit number of additions/transfers, not a
  // production-generated expected conservation residual or a stochastic retry.
  const double relative=128*std::numeric_limits<double>::epsilon()*count*(steps+1);
  for(std::size_t i=0;i<initial.size();++i)
    Require(Close(initial[i],final[i],relative),
        "periodic cross-ray transport lost a global conservative inventory");
  std::cout<<"[CSWC0628-XRAY-ADVECTION-LEVEL] cells="<<count
           <<" steps="<<steps<<" l1_density_error="<<error<<'\n';
  return error;
}

void TestAdvectionConvergence() {
  const double e0=AdvectContact(32),e1=AdvectContact(64),e2=AdvectContact(128);
  const double p0=std::log(e0/e1)/std::log(2.0),p1=std::log(e1/e2)/std::log(2.0);
  Require(p0>0.8&&p0<1.3&&p1>0.8&&p1<1.3,
      "first-order cross-ray HLL contact transport does not converge");
}

} // namespace

int main() {
  const std::pair<const char*,void(*)()> tests[]={
      {"CSWC0628-XRAY-FLUX",TestFluxAndWaves},
      {"CSWC0628-XRAY-ALE",TestAleAndOrientation},
      {"CSWC0628-XRAY-CONSERVATION",TestPairConservation},
      {"CSWC0628-XRAY-REJECTION",TestRejectionAndRollback},
      {"CSWC0628-XRAY-ADVECTION",TestAdvectionConvergence}};
  int failures=0;
  // Emit every owned identity even if an earlier independent check fails;
  // the global adapter verifies exact membership and exit/status agreement.
  for(const auto& test:tests) {
    try {test.second();std::cout<<'['<<test.first<<"] PASS bounded cross-ray kernel\n";}
    catch(const std::exception& error) {
      ++failures;std::cout<<'['<<test.first<<"] FAIL "<<error.what()<<'\n';
    }
  }
  return failures?1:0;
}
