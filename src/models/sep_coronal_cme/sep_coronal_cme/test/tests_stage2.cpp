#include "test_framework.h"

#include "sep_coronal_cme/closed_field_plasma.h"
#include "sep_coronal_cme/constants.h"
#include "sep_coronal_cme/empirical_wind.h"
#include "sep_coronal_cme/flux_tube_wind.h"
#include "sep_coronal_cme/interface_balance.h"
#include "sep_coronal_cme/pfss_harmonics.h"
#include "sep_coronal_cme/plasma_eos.h"

#include <algorithm>
#include <cmath>
#include <vector>

namespace SCCMTest { namespace {
using namespace SEP::CoronalCME;

bool Close(double a,double b,double relative=1e-9,double absolute=1e-30) {
  return std::abs(a-b)<=absolute+relative*std::max(std::abs(a),std::abs(b));
}
std::vector<WindPoint> Parker() {
  const auto p=SolveRadialIsothermalParker({0.5,0.8,1.0,1.2,2.0,5.0},2.0,1.0,3.0,4.0);
  Require(p.ok(),p.status.message);return p.value;
}
PolytropicWindSolution Polytrope() {
  // Nondimensional manufactured spherical tube: A=r^2, Phi=-2/r,
  // gamma=1.2, K=1/gamma, rho_c=1. At r_c=1 the two critical conditions
  // give u_c^2=a_c^2=Phi'/(ln A)'=1 exactly.
  std::vector<TubeGeometryPoint> g;
  for(double r:{0.8,0.9,1.0,1.2,1.6})
    g.push_back({r,r,r*r,-2.0/r,2.0/r,2.0/(r*r)});
  const auto result=SolvePolytropicTube(g,1.2,1.0/1.2,1.0,2,1e-12);
  Require(result.ok(),result.status.message);return result.value;
}
std::vector<IonSpecies> Hydrogen() {
  return {{"H+",1.0,1.0,Constants::kProtonMassKg,1.0e6}};
}
InterfaceTolerances Loose() {
  InterfaceTolerances t;t.magneticAbsoluteT=1e-8;t.magneticRelative=1;
  t.speedAbsoluteMPerS=1e3;t.speedMach=1;t.massFluxAbsolute=1;
  t.massFluxRelative=1;t.tractionAbsolutePa=1;t.tractionRelative=1;return t;
}
PrimitiveState State(Vec3 u={0,0,0},Vec3 b={0,1e-5,0}) {
  return {1e-12,1e-4,1e5,u,b};
}

void WND3D01(){auto p=Parker();for(const auto& x:p){const double y=std::pow(x.speedMPerS/2.0,2);const double rhs=4*std::log(x.radiusM)+4/x.radiusM-3;Require(Close(y-std::log(y),rhs,1e-10),"Lambert-W Parker relation mismatch");}}
void WND3D02(){double m=0,b=0,e=0;Require(CheckWindInvariants(Parker(),1.0,1e-10,&m,&b,&e).ok(),"isothermal invariants failed");Require(m<1e-12,"mass flux is not constant");}
void WND3D03(){auto solution=Polytrope();Require(CheckWindInvariants(solution.points,1.2,1e-10).ok(),"polytropic invariants failed");}
void WND3D04(){auto solution=Polytrope();for(std::size_t i=0;i<solution.points.size();++i){const auto& p=solution.points[i];const double mach=p.speedMPerS/std::sqrt(p.soundSpeedSquaredM2S2);if(i<solution.selectedCriticalIndex)Require(mach<1,"inner polytropic branch is not subsonic");else if(i==solution.selectedCriticalIndex)Require(Close(mach,1),"critical point is not sonic");else Require(mach>1,"outer polytropic branch is not supersonic");}}
void WND3D05(){auto c=FindCriticalCandidates({1,2,3,4},{1,-1,2,1},{2,2,6,4});Require(c.ok()&&c.value.size()==3,"multiple nozzle critical candidates lost");Require(c.value.back().globallyAdmissible,"global branch not selected");}
void WND3D06(){auto e=EvaluatePlasmaFromElectronDensity(1e12,1e6,Hydrogen(),false,5.0/3.0,1e-5,0.3);Require(e.ok(),e.status.message);Require(e.value.massDensityKgM3>0&&e.value.pressurePa>0&&e.value.alfvenSpeedMPerS>0&&e.value.fastSpeedMPerS>=e.value.alfvenSpeedMPerS,"composition-consistent characteristic speeds failed");}
void WND3D07(){auto result=InvertTargetSpeed([](double t){return SEP::Core::Result<double>::Success(10.0*std::sqrt(t));},200.0,100.0,1000.0);Require(result.ok()&&Close(result.value,400.0,1e-9),"target-speed inversion did not recover temperature");}
void WND3D08(){auto eta=ResolveMassLoading(2e-15,4e5,2e-5);Require(eta.ok(),eta.status.message);Require(Close(eta.value*2e-5/4e5,2e-15),"mass-flux normalization does not close");}
void WND3D09(){auto a=EvaluatePlasmaFromElectronDensity(1e12,1e6,Hydrogen(),false,1.4,1e-5);auto b=EvaluatePlasmaFromElectronDensity(1e12,1e6,Hydrogen(),false,5.0/3.0,1e-5);Require(a.ok()&&b.ok()&&!Close(a.value.soundSpeedMPerS,b.value.soundSpeedMPerS),"gamma_ad did not change characteristic speed");Require(CheckWindInvariants(Parker(),1.0,1e-10).ok(),"explicit isothermal dispatch failed");}
void WND3D10(){Require(!SolveRadialIsothermalParker({},1,1,1).ok(),"empty wind published");Require(!InvertTargetSpeed([](double t){return SEP::Core::Result<double>::Success(t);},100,1,2).ok(),"unbracketed inversion accepted");Require(!FindCriticalCandidates({1,2},{1,1},{-1,-1}).ok(),"nonexistent critical point accepted");}
void WND3D11(){auto left=ResolveMassLoading(2,3,6),right=ResolveMassLoading(4,1.5,6);Require(left.ok()&&right.ok()&&Close(left.value,right.value),"one-sided kink lost mass loading");const double stress=std::abs(2*9-4*2.25);Require(stress>0,"mandatory kink stress residual hidden");}
void WND3D12(){for(double scale:{1.0,0.5,0.25}){auto eta=ResolveMassLoading(2,3,6);Require(eta.ok()&&Close(eta.value,1.0),"loading changed under tube refinement");Require(Close(FluxTubeArea(scale,6),scale/6),"flux area changed with trace partition");}Require(!ResolveMassLoading(1,1,0).ok(),"zero-field loading accepted");}
void WND3D13(){Require(!CertifiedPositiveProfile::Create({{0,1,0,0}},1).ok(),"incomplete empirical profile accepted");auto configuration=ParseGood(Fixture()).schema5;Require(SEP::CoronalCME::CheckCapabilityAvailability(configuration,1).code==SEP::Core::StatusCode::NotImplemented,"unfinished discriminated branch fell back");}
void WND3D14(){const double omega=2e-6;auto radial=ConvertToCorotatingFieldAlignedSpeed(300,VelocityComponent::Radial,VelocityFrame::Inertial,0.6,0.5,{0,0,omega},{1,0,0},{0.6,0.8,0});Require(radial.ok()&&Close(radial.value,500),"radial-to-field speed conversion failed");Require(!ConvertToCorotatingFieldAlignedSpeed(300,VelocityComponent::Radial,VelocityFrame::Inertial,0.4,0.5,{},{},{1,0,0}).ok(),"projection failure accepted");}
void WND3D15(){std::vector<ConsumerMeasure> c={{"a",2,3,4,5,6,true,""},{"b",1,1,1,1,1,true,""}};auto census=BuildCoverageCensus(c,true);Require(census.ok()&&census.value.coveredFluxFraction==1,"finite complete D7 census failed");}
void WND3D16(){const double k=0.2;auto profile=CertifiedPositiveProfile::Create({{0,1,k,k*k},{1,std::exp(k),k*std::exp(k),k*k*std::exp(k)},{2,std::exp(2*k),k*std::exp(2*k),k*k*std::exp(2*k)}},1);Require(profile.ok(),profile.status.message);for(double x:{0.0,0.3,1.0,1.7,2.0}){auto y=profile.value.Evaluate(x);Require(y.ok()&&Close(y.value.value,std::exp(k*x),1e-11),"exact quintic-Hermite reconstruction failed");}auto blend=BlendTwoZoneWind(2,1,3,4,2,5,0.8,1);Require(blend.ok()&&blend.value.densityKgM3>0&&Close(blend.value.speedMPerS*blend.value.densityKgM3/5,0.8),"C2 two-zone continuity authority failed");Require(!BlendTwoZoneWind(2,3,1,4,2,5,0.8,1).ok(),"reversed joins accepted");}
void WND3D17(){auto direct=ResolveMassLoading(2,3,6);auto radial=ResolveMassLoadingFromRadialFlux(6,6);Require(direct.ok()&&radial.ok()&&Close(direct.value,radial.value),"eta_m authorities disagree");Require(!ResolveMassLoadingFromRadialFlux(6,0).ok(),"mass flux without mapped field accepted");}
void WND3D18(){auto h=EvaluatePlasmaFromElectronDensity(1e12,1e6,Hydrogen(),false,5.0/3.0,1e-5);auto ions=Hydrogen();ions.push_back({"He++",0.05,2,Constants::kAlphaMassKg,1e6});auto he=EvaluatePlasmaFromElectronDensity(1e12,1e6,ions,false,5.0/3.0,1e-5);auto eme=EvaluatePlasmaFromElectronDensity(1e12,1e6,Hydrogen(),true,5.0/3.0,1e-5);Require(h.ok()&&he.ok()&&eme.ok()&&he.value.massDensityKgM3!=h.value.massDensityKgM3&&eme.value.massDensityKgM3>h.value.massDensityKgM3,"composition/electron-mass conversion collapsed");}
void WND3D19(){Vec3 omega{0,0,2},x{1,0,0},t{0,1,0};auto inertial=ConvertToCorotatingFieldAlignedSpeed(5,VelocityComponent::FieldAligned,VelocityFrame::Inertial,0,0,omega,x,t);auto rotating=ConvertToCorotatingFieldAlignedSpeed(3,VelocityComponent::FieldAligned,VelocityFrame::Corotating,0,0,omega,x,t);Require(inertial.ok()&&rotating.ok()&&Close(inertial.value,rotating.value),"frame dispatch applied rotation incorrectly");Require(!ConvertToCorotatingFieldAlignedSpeed(3,VelocityComponent::Radial,VelocityFrame::Corotating,1,0.5,omega,x,t).ok(),"unsupported radial-corotating branch accepted");}
void WND3D20(){std::vector<ConsumerMeasure> c={{"covered",3,3,3,3,3,true,""},{"lost",1,2,0,1,4,false,"projection"}};auto sensitivity=BuildCoverageCensus(c,false);Require(sensitivity.ok()&&Close(sensitivity.value.coveredFluxFraction,0.75)&&sensitivity.value.rejectedIds.size()==1,"physical-measure coverage ledger failed");Require(!BuildCoverageCensus(c,true).ok(),"event-nominal missing support accepted");}

void CLS3D01(){const double rho=2,p=8,delta=3;auto h=IsothermalHydrostatic(rho,p,delta);Require(h.ok()&&Close(h.value.pressurePa/h.value.densityKgM3,p/rho),"isothermal hydrostatic balance failed");}
void CLS3D02(){auto h=PolytropicHydrostatic(2,8,1.5,1);Require(h.ok()&&Close(h.value.pressurePa/8,std::pow(h.value.densityKgM3/2,1.5)),"polytropic hydrostatic relation failed");Require(!PolytropicHydrostatic(2,8,1.5,100).ok(),"nonpositive enthalpy accepted");}
void CLS3D03(){WND3D18();auto a=EvaluatePlasmaFromElectronDensity(1e12,2e6,Hydrogen(),false,1.4,2e-5);Require(a.ok()&&a.value.fastSpeedMPerS>0,"closed composition EOS failed");}
void CLS3D04(){auto route=RouteTopology(FieldLineTopology::OpenToOuterBoundary,FieldLineTopology::ClosedBelowOuterBoundary,2,"composite-open-to-Ri");Require(route.ok()&&route.value.plasmaAuthority=="closed-hydrostatic","topology did not select closed authority");}
void CLS3D05(){Require(!CheckFootpointCompatibility(1,2,0.01).ok(),"incompatible two-footpoint closed state was blended");Require(CheckFootpointCompatibility(1,1.001,0.01).ok(),"compatible footpoints rejected");}
void CLS3D06(){Vec3 x{Constants::kSolarRadiusM,0,0};double noRotation=RotatingEffectivePotential(x,{}),rotating=RotatingEffectivePotential(x,{0,0,2e-6});Require(rotating<noRotation,"centrifugal potential was omitted");const double gravity=Constants::kSolarGravitationalParameterM3PerS2/std::pow(Constants::kSolarRadiusM,2),centrifugal=4e-12*Constants::kSolarRadiusM;Require(centrifugal/gravity<1,"rotation/gravity diagnostic invalid");}
void CLS3D07(){auto a=State({10,0,0},{0,1e-5,0}),b=a;auto balance=EvaluateSharpInterface(InterfaceKind::OpenClosedSeparatrix,BalancePolicy::DiagnosticKinematic,StateOrigin::AnalyticComposite,StateOrigin::AnalyticComposite,a,b,{1,0,0},10,Loose());Require(balance.ok()&&balance.value.kinematicGatePassed,"moving-interface relative-normal gate failed");auto stationary=EvaluateSharpInterface(InterfaceKind::OpenClosedSeparatrix,BalancePolicy::DiagnosticKinematic,StateOrigin::AnalyticComposite,StateOrigin::AnalyticComposite,a,b,{1,0,0},0,InterfaceTolerances{});Require(stationary.ok()&&!stationary.value.kinematicGatePassed,"u dot n was used without interface speed");}
void CLS3D08(){auto a=State(),b=a;Require(!EvaluateSharpInterface(InterfaceKind::OpenClosedSeparatrix,BalancePolicy::StationaryTangentialDiscontinuity,StateOrigin::AnalyticComposite,StateOrigin::AnalyticComposite,a,b,{1,0,0},0,Loose()).ok(),"analytic state claimed stationary TD");Require(!EvaluateSharpInterface(InterfaceKind::OpenClosedSeparatrix,BalancePolicy::BoundedApproximation,StateOrigin::AnalyticComposite,StateOrigin::AnalyticComposite,a,b,{1,0,0},0,Loose(),false).ok(),"bounded policy lacked evidence");Require(EvaluateSharpInterface(InterfaceKind::OpenClosedSeparatrix,BalancePolicy::BoundedApproximation,StateOrigin::AnalyticComposite,StateOrigin::AnalyticComposite,a,b,{1,0,0},0,Loose(),true).ok(),"evidenced bounded policy rejected");}
void CLS3D09(){auto a=State(),b=a;auto balance=EvaluateSharpInterface(InterfaceKind::OpenClosedSeparatrix,BalancePolicy::StationaryTangentialDiscontinuity,StateOrigin::SolverProduced,StateOrigin::Imported,a,b,{1,0,0},0,Loose());Require(balance.ok()&&balance.value.tractionNormPa==0&&balance.value.signedMassFluxJumpKgM2S==0&&balance.value.policyPassed,"manufactured sharp TD did not close");b.velocityMPerS={1,10,0};auto mutated=EvaluateSharpInterface(InterfaceKind::OpenClosedSeparatrix,BalancePolicy::StationaryTangentialDiscontinuity,StateOrigin::SolverProduced,StateOrigin::Imported,a,b,{1,0,0},0,InterfaceTolerances{});Require(mutated.ok()&&!mutated.value.tractionGatePassed&&!mutated.value.massFluxGatePassed,"normal/tangential momentum mutation hidden");}
void CLS3D10(){auto closed=EvaluateVolumeMomentumResidual({1,2,3},{4,5,6},{2,3,4},{3,4,5});Require(closed.ok()&&closed.value.normNPerM3==0,"smooth volume balance did not close");auto missing=EvaluateVolumeMomentumResidual({1,2,3},{4,5,6},{2,3,4},{0,0,0});Require(missing.ok()&&Close(missing.value.normNPerM3,std::sqrt(50.0)),"omitted force not surfaced");}
void STR3D01(){auto a=State({100,20,0},{1e-5,2e-5,0}),b=a;auto sharp=EvaluateSharpInterface(InterfaceKind::GenericOpenOpen,BalancePolicy::DiagnosticKinematic,StateOrigin::AnalyticComposite,StateOrigin::AnalyticComposite,a,b,{1,0,0},0,Loose());Require(sharp.ok()&&sharp.value.policyPassed,"generic open-open sharp interface failed");Require(sharp.value.magneticGatePassed&&sharp.value.kinematicGatePassed,"open-open incorrectly inherited separatrix gates");Require(!EvaluateSharpInterface(InterfaceKind::GenericOpenOpen,BalancePolicy::StationaryTangentialDiscontinuity,StateOrigin::Imported,StateOrigin::Imported,a,b,{1,0,0},0,Loose()).ok(),"open-open claimed stationary TD");auto volume=EvaluateVolumeMomentumResidual({},{1,0,0},{},{1,0,0});Require(volume.ok()&&volume.value.normNPerM3==0,"open-open volume residual failed");}

} // namespace

void RegisterStage2(Registry* t) {
  (*t)["WND3D01"]=WND3D01;(*t)["WND3D02"]=WND3D02;(*t)["WND3D03"]=WND3D03;(*t)["WND3D04"]=WND3D04;(*t)["WND3D05"]=WND3D05;
  (*t)["WND3D06"]=WND3D06;(*t)["WND3D07"]=WND3D07;(*t)["WND3D08"]=WND3D08;(*t)["WND3D09"]=WND3D09;(*t)["WND3D10"]=WND3D10;
  (*t)["WND3D11"]=WND3D11;(*t)["WND3D12"]=WND3D12;(*t)["WND3D13"]=WND3D13;(*t)["WND3D14"]=WND3D14;(*t)["WND3D15"]=WND3D15;
  (*t)["WND3D16"]=WND3D16;(*t)["WND3D17"]=WND3D17;(*t)["WND3D18"]=WND3D18;(*t)["WND3D19"]=WND3D19;(*t)["WND3D20"]=WND3D20;
  (*t)["CLS3D01"]=CLS3D01;(*t)["CLS3D02"]=CLS3D02;(*t)["CLS3D03"]=CLS3D03;(*t)["CLS3D04"]=CLS3D04;(*t)["CLS3D05"]=CLS3D05;
  (*t)["CLS3D06"]=CLS3D06;(*t)["CLS3D07"]=CLS3D07;(*t)["CLS3D08"]=CLS3D08;(*t)["CLS3D09"]=CLS3D09;(*t)["CLS3D10"]=CLS3D10;
  (*t)["STR3D01"]=STR3D01;
}
} // namespace SCCMTest
