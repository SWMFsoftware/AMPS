#include "test_framework.h"
#include "sep_coronal_cme/discontinuity_transport.h"
#include "sep_coronal_cme/constants.h"
#include <algorithm>
#include <cmath>
#include <limits>

namespace SCCMTest { namespace {
using namespace SEP::CoronalCME;
bool Close(double a,double b,double rel=1.0e-8,double abs=1.0e-12) {
  return std::abs(a-b)<=abs+rel*std::max(std::abs(a),std::abs(b));
}
FiniteHcsParameters Hcs() {
  FiniteHcsParameters p;
  p.sheet={ {TransportSurfaceKind::FiniteHcs,"exterior-HCS",3},{},{1,0,0},0,0};
  p.tangent={0,1,0}; p.halfThicknessM=0.1; p.fieldMagnitudeT=1;
  p.massDensityKgM3=1; p.pressurePa=2;
  p.outwardWaveEnergyJPerM3=0.03; p.inwardWaveEnergyJPerM3=0.01;
  return p;
}
SheathParameters Sheath() {
  SheathParameters p;
  p.shock={{TransportSurfaceKind::Shock,"front-A",9},{},{1,0,0},0,3};
  // Normalized manufactured SI state with c_s=1 and weak oblique field.
  p.upstream={1,0.6,{0,0.2,-0.1},
      {0.1*std::sqrt(Constants::kVacuumPermeabilityHPerM),
       0.2*std::sqrt(Constants::kVacuumPermeabilityHPerM),0}};
  p.thicknessAtEpochM=2; p.thicknessRateMPerS=0.1;
  p.beginTimeS=0; p.endTimeS=4; p.auditAreaM2=7;
  p.outwardWaveEnergyJPerM3=1.0e-4; p.inwardWaveEnergyJPerM3=2.0e-4;
  return p;
}
void HCS3D04() {
  const auto sheet=FiniteHcsSheet::Create(Hcs()); Require(sheet.ok(),sheet.status.message);
  double signedFlux=0;
  for (int i=-500;i<=500;++i) {
    const double x=i*0.002;
    const auto state=sheet.value->Evaluate({x,0,0}); Require(state.ok(),"HCS state unavailable");
    const auto& s=state.value;
    Require(Close(Norm(s.plasma.magneticFieldT),1) && s.plasma.magneticFieldT.x==0,
        "finite rotational sheet lost field magnitude or normal flux");
    Require(s.plasma.massDensityKgM3==1 && s.plasma.pressurePa==2 &&
        Close(s.waves.totalJPerM3,0.04),"HCS plasma/wave invariant changed");
    signedFlux+=s.plasma.magneticFieldT.y*0.002;
    // An independent difference estimates J x B; derivatives are not taken
    // across an ideal sign-flip stencil or borrowed PFSS/SCS clearance.
    const double dx=1.0e-6;
    const Vec3 derivative=(sheet.value->Evaluate({x+dx,0,0}).value.plasma.magneticFieldT-
        sheet.value->Evaluate({x-dx,0,0}).value.plasma.magneticFieldT)/(2*dx);
    Require(Norm(Cross(Cross({1,0,0},derivative),s.plasma.magneticFieldT))<2.0e-8,
        "finite sheet is not force free");
  }
  Require(std::abs(signedFlux)<1.0e-12,"signed sheet flux did not cancel");
  Require(!sheet.value->Evaluate({0,0,0}).value.signedWaveLabelsValid,
      "center wave polarity was invented");
  auto p=Hcs(); p.halfThicknessM=0;
  Require(!FiniteHcsSheet::Create(p).ok(),"zero thickness silently enabled finite-HCS transport");
}
void HCS3D05() {
  auto sheet=FiniteHcsSheet::Create(Hcs()); Require(sheet.ok(),sheet.status.message);
  OrbitState initial{{-0.5,0,0},{8,0.7,0.2},0};
  OrbitControls c; c.maximumThicknessFraction=0.001; c.maximumGyroAngleRad=0.001;
  const auto reference=sheet.value->Advance(initial,1,3,0.15,c);
  Require(reference.ok() && reference.value.events.size()==1,"manufactured HCS crossing was not located once");
  double previous=1.0e100;
  for (double fraction:{0.2,0.1,0.05,0.025}) {
    c.maximumThicknessFraction=fraction; c.maximumGyroAngleRad=fraction;
    const auto orbit=sheet.value->Advance(initial,1,3,0.15,c);
    Require(orbit.ok() && orbit.value.events.size()==1,"refined orbit lost HCS crossing");
    Require(Close(Norm(orbit.value.state.momentumKgMPerS),Norm(initial.momentumKgMPerS),1.0e-12),
        "HCS orbit failed scattering-frame energy conservation");
    const double error=Norm(orbit.value.state.positionM-reference.value.state.positionM)+
        std::abs(orbit.value.events[0].timeS-reference.value.events[0].timeS);
    Require(error<previous,"HCS time/drift did not converge with particle step"); previous=error;
  }
  c.maximumThicknessFraction=0.002; c.maximumGyroAngleRad=0.002;
  previous=1.0e100;
  for (double grid:{0.1,0.05,0.025,0.0125}) {
    c.fieldGridSpacingM=grid;
    const auto orbit=sheet.value->Advance(initial,1,3,0.15,c);
    Require(orbit.ok() && orbit.value.events.size()==1,"mesh-refined HCS crossing failed");
    const double error=Norm(orbit.value.state.positionM-reference.value.state.positionM);
    Require(error<previous,"HCS transverse drift failed field-grid convergence"); previous=error;
  }
  c.fieldGridSpacingM=0;
  const auto back=sheet.value->Advance(reference.value.state,1,3,-0.15,c);
  Require(back.ok() && back.value.events.size()==1 &&
      Norm(back.value.state.positionM-initial.positionM)<1.0e-6 &&
      Norm(back.value.state.momentumKgMPerS-initial.momentumKgMPerS)<1.0e-6,
      "reverse HCS crossing failed time-reversal invariant");
  const auto exact=LocatePlaneCrossing(Hcs().sheet,{-1,0,0},{1,0,0},2,4);
  Require(exact.ok() && exact.value.crossed && exact.value.event.timeS==3,
      "signed plane crossing time is wrong");
  Require(!LocatePlaneCrossing(Hcs().sheet,{0,0,0},{1,0,0},3,4).value.crossed,
      "root at segment start counted twice");
  Require(!LocatePlaneCrossing(Hcs().sheet,{1,0,0},{2,0,0},3,4).value.crossed,
      "same-side segment invented a crossing");
}
void HCS3D06() {
  for (double thickness:{0.1,0.01,0.001,1.0e-6}) {
    auto p=Hcs(); p.halfThicknessM=thickness;
    const auto sheet=FiniteHcsSheet::Create(p); Require(sheet.ok(),sheet.status.message);
    const auto minus=sheet.value->Evaluate({-1,0,0}),plus=sheet.value->Evaluate({1,0,0});
    Require(Close(minus.value.plasma.magneticFieldT.y,-1,1.0e-7) &&
        Close(plus.value.plasma.magneticFieldT.y,1,1.0e-7),"thin exterior HCS limit lost one-sided sector field");
  }
  const auto ideal=RestoreSector({0,1,0},MagneticSector::Negative,true,true,InterfaceSide::Minus);
  Require(ideal.ok() && ideal.value.valueT.y==-1,"separate ideal-HCS operator changed");
  auto p=Hcs(); p.sheet.identity.kind=TransportSurfaceKind::PfssScsTransition;
  Require(!FiniteHcsSheet::Create(p).ok(),"exterior HCS silently qualified composite transition");
  p.sheet.identity.kind=TransportSurfaceKind::Separatrix;
  Require(!FiniteHcsSheet::Create(p).ok(),"HCS reused separatrix identity");
}
void SHEATH3D01() {
  const auto sheath=PlanarShockSheath::Create(Sheath()); Require(sheath.ok(),sheath.status.message);
  const auto event=LocatePlaneCrossing(Sheath().shock,{1,0,0},{1,0,0},0,1);
  Require(event.ok() && event.value.crossed && Close(event.value.event.timeS,1.0/3.0),
      "moving front crossing time is incorrect");
  // Refining a straight trajectory's cadence does not alter the physical root.
  for (int steps:{3,7,31,127}) {
    unsigned crossings=0;
    for (int i=0;i<steps;++i) {
      const auto e=LocatePlaneCrossing(Sheath().shock,{1,0,0},{1,0,0},double(i)/steps,double(i+1)/steps);
      Require(e.ok(),"shock event location failed");
      if (e.value.crossed) { ++crossings; Require(Close(e.value.event.timeS,1.0/3.0),"shock event depended on cadence"); }
    }
    Require(crossings==1,"shock generation counted more than once");
  }
  auto parallel=Sheath(); parallel.upstream.velocityMPerS={};
  parallel.upstream.magneticFieldT={0.1*std::sqrt(Constants::kVacuumPermeabilityHPerM),0,0};
  const auto normalSheath=PlanarShockSheath::Create(parallel); Require(normalSheath.ok(),normalSheath.status.message);
  OrbitState initial{{1,0,0},{0.5,0.2,0.1},0};
  CrossingLedger orbitLedger; OrbitControls controls; controls.maximumGyroAngleRad=0.01;
  const auto orbit=normalSheath.value->Advance("orbit",initial,1,1,0.5,&orbitLedger,controls);
  Require(orbit.ok() && orbit.value.events.size()==1 && orbitLedger.Keys().size()==1 &&
      Close(orbit.value.events[0].timeS,0.4,1.0e-8) &&
      Close(Norm(orbit.value.state.momentumKgMPerS),Norm(initial.momentumKgMPerS),1.0e-12),
      "split sheath orbit missed event or changed zero-electric-field energy");
  const auto reverseOrbit=normalSheath.value->Advance("reverse-orbit",orbit.value.state,1,1,-0.5,&orbitLedger,controls);
  Require(reverseOrbit.ok() && reverseOrbit.value.events.size()==1 &&
      Norm(reverseOrbit.value.state.positionM-initial.positionM)<1.0e-8,
      "forward/reverse split sheath orbit did not converge");
  CrossingLedger unpublished;
  const auto failed=normalSheath.value->Advance("rear-exit",initial,1,1,2,&unpublished,controls);
  Require(!failed.ok() && unpublished.Keys().empty(),"failed sheath step published a partial crossing ledger");
  const double mass=Constants::kProtonMassKg,c=Constants::kSpeedOfLightMPerS;
  FourMomentum momentum{std::hypot(mass*c*c,1.0e-19*c),{1.0e-19,0,0}};
  CrossingLedger ledger;
  const auto crossed=sheath.value->Cross("particle-7",event.value.event,momentum,&ledger);
  Require(crossed.ok() && crossed.value.applied,"valid shock event not committed");
  Require(Close(FourMomentumInvariant(crossed.value.localPlasma)/FourMomentumInvariant(momentum),1,1.0e-12),
      "plasma-frame crossing broke relativistic mass shell");
  const auto repeat=sheath.value->Cross("particle-7",event.value.event,momentum,&ledger);
  Require(repeat.ok() && !repeat.value.applied && ledger.Keys().size()==1,"duplicate shock callback reaccelerated particle");
  auto reverse=event.value.event; std::swap(reverse.incoming,reverse.outgoing);
  const auto backward=sheath.value->Cross("reverse-particle",reverse,momentum,&ledger);
  Require(backward.ok() && Close(backward.value.shockFrameEnergyJ,crossed.value.shockFrameEnergyJ,1.0e-12),
      "forward/reverse shock-frame energy changed");
  CrossingLedger restored; Require(restored.Restore(ledger.Keys()).ok(),"crossing ledger restore failed");
  Require(!sheath.value->Cross("particle-7",event.value.event,momentum,&restored).value.applied,
      "restart re-applied committed shock crossing");
  auto magnetic=event.value.event; magnetic.surface.kind=TransportSurfaceKind::FiniteHcs;
  magnetic.surface.stableId="HCS";
  auto excluded=magnetic; excluded.surface.kind=TransportSurfaceKind::PfssScsTransition;
  excluded.surface.stableId="S-tr";
  const auto order=OrderCrossingEvents({event.value.event,magnetic,excluded},1.0e-9);
  Require(order[0].surface.kind==TransportSurfaceKind::PfssScsTransition &&
      order[1].surface.kind==TransportSurfaceKind::FiniteHcs && order[2].surface.kind==TransportSurfaceKind::Shock,
      "coincident exclusion/magnetic/shock order is not explicit");
  momentum.totalEnergyJ=std::numeric_limits<double>::infinity();
  Require(!sheath.value->Cross("bad",event.value.event,momentum,&ledger).ok(),"nonfinite crossing energy was published");
}
void SHEATH3D02() {
  const auto sheath=PlanarShockSheath::Create(Sheath()); Require(sheath.ok(),sheath.status.message);
  for (int cells:{4,16,64,256}) for (double time:{0.0,1.0,2.0,4.0}) {
    const auto audit=sheath.value->Audit(cells,time);
    Require(audit.ok() && audit.value.passed && audit.value.divergenceRelative<1.0e-12 &&
        audit.value.energyRelative<1.0e-8 && audit.value.waveFluxRelative<1.0e-12,
        "finite sheath global divergence/mass/flux/energy audit failed refinement");
    Require(audit.value.massInventoryRateKgPerS>0 && audit.value.energyInventoryRateW>0,
        "time-dependent sheath growth omitted inventory accounting");
  }
  for (auto sector:{MagneticSector::Positive,MagneticSector::Negative}) {
    auto reversed=Sheath(); reversed.magneticSector=sector;
    if (sector==MagneticSector::Negative) reversed.upstream.magneticFieldT=-1.0*reversed.upstream.magneticFieldT;
    const auto qualified=PlanarShockSheath::Create(reversed);
    Require(qualified.ok() && qualified.value->Audit(32,1).value.passed,
        "sheath wave flux inferred magnetic sector from shock normal");
  }
  const auto& p=Sheath();
  Require(!sheath.value->Evaluate({0,0,0},0).ok(),"unsided exact-front state accepted");
  Require(sheath.value->Evaluate({0,0,0},0,true,InterfaceSide::Minus).value.massDensityKgM3>p.upstream.massDensityKgM3,
      "downstream one-sided state was not the solved RH state");
  Require(!sheath.value->Evaluate({-3,0,0},0).ok(),"sheath extrapolated behind its rear boundary");
  Require(!sheath.value->Evaluate({0,0,0},5).ok(),"sheath extrapolated beyond time coverage");
  auto invalid=p; invalid.thicknessRateMPerS=-1;
  Require(!PlanarShockSheath::Create(invalid).ok(),"sheath allowed rear boundary to overtake front");
  invalid=p; invalid.outwardWaveEnergyJPerM3=1.0e6;
  Require(!PlanarShockSheath::Create(invalid).ok(),"large waves silently ignored their MHD backreaction");
  invalid=p; invalid.shock.normalSpeedMPerS=0.1;
  Require(!PlanarShockSheath::Create(invalid).ok(),"subfast state fabricated downstream sheath");
}
void SHEATH3D03() {
  const auto hcs=FiniteHcsSheet::Create(Hcs());
  const auto sheath=PlanarShockSheath::Create(Sheath());
  Require(hcs.ok() && sheath.ok(),"independent provider preparation failed");
  Require(ValidateDiscontinuityCapabilities(true,false,hcs.value,{}).ok(),"HCS incorrectly required sheath");
  Require(ValidateDiscontinuityCapabilities(false,true,{},sheath.value).ok(),"sheath incorrectly required HCS");
  Require(!ValidateDiscontinuityCapabilities(true,true,hcs.value,{}).ok(),"HCS enabled missing sheath");
  Require(!ValidateDiscontinuityCapabilities(true,true,{},sheath.value).ok(),"sheath enabled missing HCS");
  Require(!ValidateDiscontinuityCapabilities(true,false,hcs.value,{},true).ok(),"exterior HCS enabled composite transition");
}
} // namespace
void RegisterStage11(Registry* tests) {
  tests->emplace("HCS3D04",HCS3D04); tests->emplace("HCS3D05",HCS3D05);
  tests->emplace("HCS3D06",HCS3D06); tests->emplace("SHEATH3D01",SHEATH3D01);
  tests->emplace("SHEATH3D02",SHEATH3D02); tests->emplace("SHEATH3D03",SHEATH3D03);
}
} // namespace SCCMTest
