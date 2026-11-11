#include "sep_corona_swcme/cross_ray_fluid.h"
#include "sep_coronal_cme/configuration_parser.h"

#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>

#include <type_traits>

static_assert(!std::is_default_constructible<SEP::CoronaSwcme::CrossRay::StaticFluidState>::value,
    "invalid empty authority must not bypass its validated factory");

namespace {
namespace XR=SEP::CoronaSwcme::CrossRay;
using SEP::CoronalCME::Cross;
using SEP::CoronalCME::Dot;
using SEP::CoronalCME::Norm;
using SEP::CoronalCME::Vec3;
constexpr double kPi=3.141592653589793238462643383279502884;
const XR::Thermodynamics kGas{5.0/3.0,2};
const Vec3 kB{0.4,-0.3,0.2},kV{0.1,-0.2,0.3};
constexpr double kRho=2,kPressure=3;

void Require(bool ok,const std::string& message) {
  if(!ok)throw std::runtime_error(message);
}

double Energy() {
  return kPressure/(kGas.gammaAdiabatic-1)+0.5*kRho*Dot(kV,kV)+
      Dot(kB,kB)/(2*kGas.magneticPermeabilityNPerA2);
}

std::shared_ptr<const XR::StaticFluidState> State() {
  XR::StaticMeshInput in;in.radiusEdgesM={1,2,4};
  for(int j=0;j<=7;++j)in.muEdges.push_back(-1+2.0*j/7);
  for(int k=0;k<=24;++k)in.phiEdgesRad.push_back(2*kPi*k/24);
  in.periodicPhi=true;in.maximumCells=336;
  const auto mesh=XR::StaticSphericalMesh::Create(in,
      SEP::CoronalCME::ComputeContentChecksum("independent-conservative-fluid-CT-reference-v1"),10,3);
  Require(mesh.ok(),mesh.status.message);
  std::vector<double> values;
  for(const auto& face:mesh.value->Faces())
    values.push_back(Dot(kB,face.measure.outwardAreaVectorM2));
  const auto field=XR::StaticMagneticFlux::FromFaceFlux(mesh.value,values);
  Require(field.ok(),field.status.message);
  const auto state=XR::StaticFluidState::Create(field.value,
      std::vector<XR::FluidDensity>(336,{kRho,kRho*kV,Energy()}),kGas);
  Require(state.ok(),state.status.message);return state.value;
}

XR::FluidFluxInterval Interval(const XR::StaticFluidState& state,double dt) {
  const auto& field=*state.Magnetic();const auto& mesh=*field.Mesh();
  XR::FluidFluxInterval out;
  out.magnetic.meshFingerprint=mesh.Fingerprint();
  out.magnetic.eventFingerprint=mesh.EventFingerprint();
  out.magnetic.fieldGeneration=field.Generation();
  out.magnetic.startS=field.EpochS();out.magnetic.endS=field.EpochS()+dt;
  out.magnetic.valuesVolt.resize(mesh.Edges().size());
  out.faces.resize(mesh.Faces().size());out.sources.resize(mesh.Cells().size());
  return out;
}

void CheckBudget(const XR::FluidStepLedger& ledger) {
  Require(std::abs(ledger.residual.massKg)<2e-12&&
      Norm(ledger.residual.momentumKgMPerS)<2e-12&&
      std::abs(ledger.residual.totalEnergyJ)<4e-12,
      "raw extensive mass/vector-momentum/total-energy budget does not close");
}

void Uniform() {
  auto current=State();const auto& mesh=*current->Magnetic()->Mesh();
  // The exact constant-state flux is linear in the surface normal. Its
  // integral therefore uses the analytic area VECTOR, not scalar curved
  // area multiplied by a representative radial normal. This independent
  // reference exercises pressure, Maxwell stress and Poynting work together.
  for(int step=0;step<12;++step) {
    auto interval=Interval(*current,0.025);
    const double pt=kPressure+Dot(kB,kB)/(2*kGas.magneticPermeabilityNPerA2);
    for(std::size_t f=0;f<mesh.Faces().size();++f) {
      const auto s=mesh.Faces()[f].measure.outwardAreaVectorM2;
      const double vs=Dot(kV,s),bs=Dot(kB,s);
      interval.faces[f]={kRho*vs,kRho*vs*kV+pt*s-
          (bs/kGas.magneticPermeabilityNPerA2)*kB,
          (Energy()+pt)*vs-Dot(kV,kB)*bs/kGas.magneticPermeabilityNPerA2};
    }
    const auto electric=-1.0*Cross(kV,kB);
    for(std::size_t e=0;e<mesh.Edges().size();++e) {
      const auto& edge=mesh.Edges()[e];
      interval.magnetic.valuesVolt[e]=Dot(electric,
          mesh.NodesHciM()[edge.lastNode]-mesh.NodesHciM()[edge.firstNode]);
    }
    const auto next=current->Advance(interval);Require(next.ok(),next.status.message);
    for(const auto& p:next.value->Primitives())Require(
        std::abs(p.densityKgM3-kRho)<2e-13&&Norm(p.velocityMPerS-kV)<2e-13&&
        Norm(p.magneticFieldT-kB)<2e-13&&std::abs(p.pressurePa-kPressure)<2e-13,
        "static curved-mesh ideal-MHD uniform state is not preserved");
    Require(next.value->Magnetic()->Generation()==current->Magnetic()->Generation()+1,
        "fluid and CT generations did not advance together");
    CheckBudget(next.value->StepLedger());current=next.value;
  }
}

void Conservation() {
  const auto state=State();const auto& mesh=*state->Magnetic()->Mesh();
  auto interval=Interval(*state,0.1);std::size_t chosen=XR::NoMeshIndex;
  for(std::size_t f=0;f<mesh.Faces().size();++f) {
    const auto& face=mesh.Faces()[f];
    if(face.axis==XR::MeshAxis::Mu&&face.leftCell!=XR::NoMeshIndex&&
        face.rightCell!=XR::NoMeshIndex){chosen=f;break;}
  }
  Require(chosen!=XR::NoMeshIndex,"no actual cross-ray test interface");
  interval.faces[chosen]={0.01,{0.003,-0.002,0.004},0.02};
  const auto next=state->Advance(interval);Require(next.ok(),next.status.message);
  const auto& face=mesh.Faces()[chosen];
  const double dt=interval.magnetic.endS-interval.magnetic.startS;
  for(std::size_t c=0;c<state->Cells().size();++c) {
    const double signedStep=c==face.leftCell?-dt:(c==face.rightCell?dt:0);
    const double v=mesh.Cells()[c].measure.volumeM3;
    const auto& value=next.value->Cells()[c];
    Require(std::abs(value.massKgM3-(kRho+signedStep*0.01/v))<2e-14&&
        Norm(value.momentumKgM2S-(kRho*kV+(signedStep/v)*Vec3{0.003,-0.002,0.004}))<2e-14&&
        std::abs(value.totalEnergyJPerM3-(Energy()+signedStep*0.02/v))<2e-14,
        "lateral face transfer has wrong sign, volume or Cartesian momentum/energy");
  }
  CheckBudget(next.value->StepLedger());
  Require(next.value->StepLedger().outwardBoundaryTransfer.massKg==0&&
      next.value->StepLedger().volumeAddition.massKg==0,
      "internal lateral mass was mislabeled as an external/source contribution");

  // An inner radial boundary has no left cell. A positive canonical flux
  // enters the domain and must appear as NEGATIVE outward transfer.
  interval=Interval(*state,0.1);
  for(std::size_t f=0;f<mesh.Faces().size();++f)if(mesh.Faces()[f].leftCell==XR::NoMeshIndex) {
    interval.faces[f]={0.01,{0.001,0.002,-0.003},0.02};break;
  }
  const auto supplied=state->Advance(interval);Require(supplied.ok(),supplied.status.message);
  CheckBudget(supplied.value->StepLedger());
  Require(std::abs(supplied.value->StepLedger().outwardBoundaryTransfer.massKg+0.01*dt)<1e-15&&
      std::abs(supplied.value->StepLedger().outwardBoundaryTransfer.totalEnergyJ+0.02*dt)<1e-15,
      "inner-boundary mass/energy ledger has the wrong outward sign");
}

std::vector<double> LinearElectricLines(const XR::StaticSphericalMesh& mesh,Vec3 bdot) {
  // Manufactured EXTERNALLY prescribed E=-bdot cross X/2, curl E=-bdot.
  // It tests the coupling kernel with explicit electromagnetic source work;
  // it is not an ideal-MHD edge solver or evidence for an unforced CME.
  const double x[]={0.18343464249564980494,0.52553240991632898582,
      0.79666647741362673959,0.96028985649753623168};
  const double w[]={0.36268378337836198297,0.31370664587788728734,
      0.22238103445337447054,0.10122853629037625915};
  std::vector<double> result;
  for(std::size_t edge=0;edge<mesh.Edges().size();++edge) {
    long double line=0;
    for(int k=0;k<4;++k)for(double sign:{-1.0,1.0}) {
      const auto sample=mesh.SampleEdge(edge,0.5*(1+sign*x[k]));
      Require(sample.ok(),sample.status.message);
      line-=0.25L*w[k]*Dot(Cross(bdot,sample.value.positionHciM),sample.value.derivativeM);
    }
    result.push_back(static_cast<double>(line));
  }
  return result;
}

void EnergyCoupling() {
  const auto state=State();const auto& mesh=*state->Magnetic()->Mesh();
  auto interval=Interval(*state,0.1);const Vec3 bdot{0.1,0.03,-0.02};
  const double dt=interval.magnetic.endS-interval.magnetic.startS;
  const Vec3 nextB=kB+dt*bdot;
  const double magneticGain=(Dot(nextB,nextB)-Dot(kB,kB))/(2*kGas.magneticPermeabilityNPerA2);
  interval.magnetic.valuesVolt=LinearElectricLines(mesh,bdot);
  const auto fixedE=state->Advance(interval);Require(fixedE.ok(),fixedE.status.message);
  for(std::size_t c=0;c<state->Cells().size();++c)Require(
      fixedE.value->Cells()[c].totalEnergyJPerM3==Energy()&&
      Norm(fixedE.value->Primitives()[c].magneticFieldT-nextB)<3e-13&&
      std::abs(fixedE.value->Primitives()[c].pressurePa-
          (kPressure-(kGas.gammaAdiabatic-1)*magneticGain))<3e-13,
      "CT coupling silently reset total energy or ignored the new magnetic energy");
  CheckBudget(fixedE.value->StepLedger());
  for(std::size_t c=0;c<state->Cells().size();++c)
    interval.sources[c].totalPowerW=magneticGain*mesh.Cells()[c].measure.volumeM3/dt;
  const auto supplied=state->Advance(interval);Require(supplied.ok(),supplied.status.message);
  for(const auto& value:supplied.value->Primitives())Require(
      std::abs(value.pressurePa-kPressure)<3e-13,
      "explicit electromagnetic source power was not applied to conserved total energy");
  Require(supplied.value->StepLedger().volumeAddition.totalEnergyJ>0,
      "electromagnetic source work is missing from the extensive budget");
  CheckBudget(supplied.value->StepLedger());
}

void Rejection() {
  const auto state=State();const auto fingerprint=state->Fingerprint();
  auto bad=Interval(*state,0.1);bad.faces.pop_back();
  Require(!state->Advance(bad).ok(),"missing face flux was interpreted as zero");
  bad=Interval(*state,0.1);bad.sources.pop_back();
  Require(!state->Advance(bad).ok(),"missing volume work was interpreted as zero");
  bad=Interval(*state,0.1);bad.magnetic.fieldGeneration++;
  Require(!state->Advance(bad).ok(),"stale magnetic/fluid authority was accepted");
  bad=Interval(*state,0.1);bad.faces[0].momentumN.y=std::numeric_limits<double>::quiet_NaN();
  Require(!state->Advance(bad).ok(),"nonfinite transverse momentum flux was accepted");
  bad=Interval(*state,0.1);bad.sources[0].totalPowerW=-1e6;
  Require(!state->Advance(bad).ok(),"nonpositive thermal energy was repaired with a floor");
  bad=Interval(*state,0.1);
  bad.magnetic.valuesVolt=LinearElectricLines(*state->Magnetic()->Mesh(),{1000,0,0});
  Require(!state->Advance(bad).ok(),"unbudgeted magnetic energy received an energy reset");
  auto cells=state->Cells();cells[0].massKgM3=0;
  Require(!XR::StaticFluidState::Create(state->Magnetic(),cells,kGas).ok(),
      "zero initial fluid density received a floor");
  Require(!XR::StaticFluidState::Create(state->Magnetic(),state->Cells(),{1,2}).ok(),
      "invalid initial thermodynamics were accepted");
  Require(state->Fingerprint()==fingerprint&&state->Magnetic()->EpochS()==10&&
      state->Magnetic()->Generation()==3,"failed epoch preparation mutated committed state");
}
} // namespace

int main() {
  try {
    Uniform();std::cout<<"[CSWC0628-XRAY-FLUID-UNIFORM] PASS static curved-face MHD/CT conservation kernel\n";
    Conservation();std::cout<<"[CSWC0628-XRAY-FLUID-CONSERVATION] PASS cross-ray and boundary extensive budgets\n";
    EnergyCoupling();std::cout<<"[CSWC0628-XRAY-FLUID-ENERGY] PASS conserved total energy with new CT field and explicit work\n";
    Rejection();std::cout<<"[CSWC0628-XRAY-FLUID-REJECTION] PASS atomic epoch and positive thermal-energy rejection\n";
  }catch(const std::exception& e) {
    std::cerr<<"cross-ray fluid FAIL: "<<e.what()<<'\n';return 1;
  }
  return 0;
}
