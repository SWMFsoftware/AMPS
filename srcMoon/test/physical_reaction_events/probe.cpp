// Thin adapter for the SAME source-bound production event composer. The
// invented kg/state-energy inputs below are explicitly MATHEMATICAL fixtures,
// not a physical sodium state inventory, observations or an admitted profile.
#include "srcMoon/species/PhysicalReactionEvents.h"
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>

namespace E=Moon::Species::PhysicalReactionEvents;
namespace K=Moon::Species::ReducedKinematics;
namespace R=Moon::Species::ReactionExecution;
namespace P=Moon::Species::Photochemistry;
namespace C=Moon::Species::IonCollector;
using Real=K::Real;
std::ofstream checks,metrics;
unsigned total=0,failed=0;
void Check(const std::string& name,bool passed){
  checks<<name<<','<<(passed?1:0)<<'\n';++total;if(!passed){++failed;std::cerr<<"FAIL "<<name<<'\n';}
}
void Metric(const std::string& group,const std::string& key,Real value){metrics<<group<<','<<key<<','<<value<<'\n';}
void Vector(const std::string& group,const std::string& key,const K::Vector3& value){
  for(unsigned i=0;i<3;++i)Metric(group,key+std::to_string(i),value[i]);
}
K::ResolvedState State(const std::string& id,int charge,Real formation){
  K::StateDefinition d;d.identity=id;d.state="ANALYTICAL-STATE-NOT-A-PHYSICAL-CROSSWALK";
  d.elementaryCharge=charge;if(id!="e")d.nuclei={{"23Na",1}};
  d.formation={formation,"ANALYTICAL-fixture-state-gap-not-NIST"};
  d.excitation={0,"ANALYTICAL-explicit-zero"};
  K::ResolvedState s;std::string error;
  if(!K::ResolveState(d,{{"23Na",{4e-26L,"ANALYTICAL-4e-26kg-not-isotope-data"}}},
      {1e-30L,"ANALYTICAL-1e-30kg-not-electron-data"},&s,&error))throw std::runtime_error(error);
  return s;
}
// I22 atomic pilot (approved 2026-11-26): analytical states that carry one
// declared nucleus. One shared mathematical isotope table keeps the three
// states in a common reference inventory; masses are NOT isotope data.
K::ResolvedState NucleusState(const std::string& id,const std::string& nucleus,int charge,Real formation){
  K::StateDefinition d;d.identity=id;d.state="ANALYTICAL-STATE-NOT-A-PHYSICAL-CROSSWALK";
  d.elementaryCharge=charge;if(!nucleus.empty())d.nuclei={{nucleus,1}};
  d.formation={formation,"ANALYTICAL-fixture-state-gap-not-NIST"};
  d.excitation={0,"ANALYTICAL-explicit-zero"};
  K::ResolvedState s;std::string error;
  const K::MassValue mass{4e-26L,"ANALYTICAL-4e-26kg-not-isotope-data"};
  if(!K::ResolveState(d,{{"4He",mass},{"20Ne",mass},{"40Ar",mass},{"39K",mass}},
      {1e-30L,"ANALYTICAL-1e-30kg-not-electron-data"},&s,&error))throw std::runtime_error(error);
  return s;
}
void Catalog(E::Definition* d){
  d->catalog.resize(3);std::string error;double internal=0;
  auto a=E::ReferenceCatalogEntry(d->parent.state,0,&d->catalog[0],&internal);
  auto b=E::ReferenceCatalogEntry(d->ion.state,1,&d->catalog[1],&internal);
  auto c=E::ReferenceCatalogEntry(d->electron.state,-1,&d->catalog[2],&internal);
  if(!a.ok()||!b.ok()||!c.ok())throw std::runtime_error("fixture state binding failed");
}
E::Definition Definition(const P::Registry& source,P::SolarCase solar,const std::string& fixtureSha){
  std::string error;std::vector<P::Channel> rows;
  if(!source.ChannelsForParent("Na",&rows,&error)||rows.size()!=1)throw std::runtime_error(error);
  E::Definition d;d.identity="ANALYTICAL-composition-controls-only";
  d.sourceProfileIdentity="ANALYTICAL-fixture-not-physical-Na-profile";d.sourceProfileSha256=fixtureSha;
  d.stateCrosswalkIdentity="ANALYTICAL-gap-set-to-D04-threshold-not-physical-evidence";d.stateCrosswalkSha256=fixtureSha;
  d.energyReferenceIdentity="ANALYTICAL-zero-reference";d.sourceStateCrosswalkAvailable=true;
  d.electronVoltJ=1.602176634e-19L;d.electronVoltSourceIdentity="SI-exact-eV-definition";
  d.channelId=rows.front().id;d.solar=solar;
  const Real threshold=static_cast<Real>(rows.front().thresholdEnergyEv)*d.electronVoltJ;
  d.parent.state=State("Na",0,0);d.parent.catalogIndex=0;
  d.ion.state=State("Na+",1,threshold);d.ion.catalogIndex=1;
  d.electron.state=State("e",-1,0);d.electron.catalogIndex=2;
  E::StateBinding* bindings[]={&d.parent,&d.ion,&d.electron};
  for(auto* b:bindings){b->physicalRestMass={b->state.referenceMassKg(),"ANALYTICAL-reference-equals-physical-mass-control"};b->modelErrorLimits={0,0,.001L};}
  Catalog(&d);
  // Predeclared SOFTWARE arithmetic budgets, not physical source uncertainties.
  d.numerical.maximumClassicalSpeedMS=3000000;d.numerical.unitVectorTolerance=1e-15L;
  d.numerical.referenceMassKg={1e-48L,1e-12L};d.numerical.momentumKgMS={1e-38L,1e-12L};
  d.numerical.energyJ={1e-31L,1e-12L};d.numerical.thresholdJ={1e-32L,1e-12L};d.rotation={0,1e-15L};
  d.plannerAuthorization.qualified=true;
  d.plannerAuthorization.identity="ANALYTICAL-planner-storage-control-not-physical-qualification";
  d.plannerAuthorization.evidenceIdentity=fixtureSha;d.plannerAuthorization.eventFrame="ANALYTICAL-event-axes";
  d.plannerAuthorization.maximumClassicalSpeedMS=3000000;
  d.plannerAuthorization.bounds={{1e-46,1e-12},{1e-37,1e-12},{1e-30,1e-12}};
  return d;
}
E::EventInput Input(unsigned mode){
  E::EventInput i;i.parentAtEvent.species=0;i.parentAtEvent.referenceWeight=7;i.parentAtEvent.individualCorrection=.375;
  i.parentAtEvent.positionM={{123,456,789}};i.parentAtEvent.internalEnergyJ=0;
  R::TraceInterval::Create(100,1,.375,&i.interval,nullptr);
  R::OpticalDepthClock::FromOpenUniform(.3,&i.freshIonClock,nullptr);
  i.ionReferenceWeight=.4;i.polarOpenUniform=.5;i.azimuthOpenUniform=.25;
  i.frame.inertialIdentity="ANALYTICAL-inertial-axes";i.frame.eventIdentity="ANALYTICAL-event-axes";
  i.frame.transformSourceIdentity="ANALYTICAL-exact-rotation-and-offset";
  i.frame.epochConventionIdentity="ANALYTICAL-frozen-basis-and-advancing-transport";
  i.frame.physicalEventEpochEtS=100.375;i.frame.geometryBasisEpochEtS=90;
  i.frame.inertialToEventRotation={{1,0,0,0,1,0,0,0,1}};
  i.frame.velocityOffsetSupplied=true;
  if(mode==1)i.parentAtEvent.velocityMS={{10000,-20000,30000}};
  if(mode==2){i.frame.inertialToEventRotation={{0,-1,0,1,0,0,0,0,1}};
    i.frame.velocityOffsetEventMS={{100,-50,25}};i.parentAtEvent.velocityMS={{20100,9950,30025}};}
  if(mode==3){
    // Actual rounded SPICE coefficients and large affine velocity offset from
    // the refused commissioning event. These are frame-arithmetic fixtures,
    // not observations or an invented physical sodium state profile.
    i.frame.inertialToEventRotation={{0.560211314016633821211L,-0.760052659393080420891L,-0.329367938022605444282L,
      0.828347105332689293L,0.515022425867632138718L,0.220438140847663277544L,
      0.00208727926048927385841L,-0.396322918594075035426L,0.918108810252010965058L}};
    i.frame.velocityOffsetEventMS={{-192.106854304402989486L,29386.0300909183624842L,-4.30322444344710675068e-13L}};
    i.parentAtEvent.velocityMS={{0,0,0}};
  }
  i.incomingPhotonDirectionInertial={{0,0,-1}};i.pointSunVisible=1;i.distanceAu=2;
  i.photonGeometryIdentity="ANALYTICAL-propagation-minus-z-not-observation";return i;
}
void Record(const std::string& group,const E::Definition& d,const E::EventInput& input,const E::Composition& value){
  const auto& r=value.receipt;
  Metric(group,"mass-parent",d.parent.state.referenceMassKg());Metric(group,"mass-ion",d.ion.state.referenceMassKg());Metric(group,"mass-electron",d.electron.state.referenceMassKg());
  Metric(group,"threshold-J",r.nativeThresholdJ);Metric(group,"mean-J",r.matchedSourceMeanExcessJ);
  Metric(group,"physical-threshold-J",r.physicalThresholdJ);Metric(group,"physical-minus-source-J",r.physicalMinusSourceThresholdJ);
  Metric(group,"input-physical-threshold-J",d.ion.state.definition().formation.joules);
  Metric(group,"source-Q-J",r.stateQJ);Metric(group,"physical-epoch",input.frame.physicalEventEpochEtS);
  Metric(group,"basis-epoch",input.frame.geometryBasisEpochEtS);Metric(group,"native-epsilon",std::numeric_limits<Real>::epsilon());
  Vector(group,"input-parent-event-v",{{input.parentAtEvent.velocityMS[0],input.parentAtEvent.velocityMS[1],input.parentAtEvent.velocityMS[2]}});
  Vector(group,"input-offset-v",input.frame.velocityOffsetEventMS);Vector(group,"input-photon-n",input.incomingPhotonDirectionInertial);
  for(unsigned k=0;k<9;++k)Metric(group,"input-R"+std::to_string(k),input.frame.inertialToEventRotation[k]);
  Vector(group,"parent-inertial-v",r.parentInertialVelocityMS);
  Metric(group,"photon-rest-J",r.incomingPhoton.restEnergyJ);Metric(group,"photon-inertial-J",r.incomingPhoton.inertialEnergyJ);Metric(group,"doppler",r.incomingPhoton.dopplerFactor);
  Vector(group,"photon-inertial-p",r.incomingPhoton.inertialMomentumKgMS);
  Vector(group,"photon-rest-p",r.incomingPhoton.restMomentumKgMS);
  Vector(group,"ion-inertial-v",r.inertialProducts.products[0].velocityMS);Vector(group,"electron-inertial-v",r.outgoingElectron.inertialVelocityMS);
  const auto& ion=value.plan.products().front();Vector(group,"ion-event-v",{{ion.velocityMS[0],ion.velocityMS[1],ion.velocityMS[2]}});
  Vector(group,"electron-event-v",r.outgoingElectron.eventVelocityMS);
  Vector(group,"bath-inertial-p",r.bathInertialMomentumChangeKgMS);Vector(group,"bath-event-p",r.bathEventMomentumChangeKgMS);
  Metric(group,"bath-inertial-J",r.bathInertialEnergyChangeJ);Metric(group,"bath-event-J",r.bathEventEnergyChangeJ);Metric(group,"bath-mass-kg",r.bathModelMassChangeKg);
  Metric(group,"source-rate-1au",r.sourceRate1AuS);Metric(group,"event-rate",r.eventRateS);
  Metric(group,"ion-reference-weight",ion.referenceWeight);Metric(group,"ion-individual-correction",ion.individualCorrection);Metric(group,"ion-number",ion.representedNumber);
  Metric(group,"ion-internal-J",ion.internalEnergyJ);Metric(group,"remaining-s",ion.interval.remainingS());
}
int main(int argc,char** argv){
  if(argc!=4){std::cerr<<"usage: probe FROZEN_PHIDRATES_DIRECTORY ANALYTICAL_FIXTURE_SHA OUTPUT_PREFIX\n";return 2;}
  try {
    checks.open(std::string(argv[3])+".checks.csv");metrics.open(std::string(argv[3])+".metrics.csv");
    if(!checks||!metrics)throw std::runtime_error("cannot open adapter output files");
    checks<<"name,passed\n";metrics<<"case,key,value\n"<<std::scientific<<std::setprecision(std::numeric_limits<Real>::max_digits10);
    P::Registry source;std::string error;
    if(!source.LoadFrozenD04(P::FrozenV1Sources(argv[1]),&error))throw std::runtime_error("frozen source load: "+error);
    E::SodiumPhotoionization provider;
    for(int convention=0;convention<2;++convention)for(int endpoint=0;endpoint<2;++endpoint){
      auto d=Definition(source,static_cast<P::SolarCase>(endpoint),argv[2]);
      if(convention){
        // Explicit mathematical perturbation, NOT an evaluated sodium energy.
        // The real source photon/rate is unchanged. The independently checked
        // input formation gap now differs; a default strict provider refuses it.
        const Real gap=d.ion.state.definition().formation.joules+1e-6L*d.electronVoltJ;
        d.ion.state=State("Na+",1,gap);Catalog(&d);
        d.photonEnergyPolicy=K::PhotonEnergyPolicy::SourceEffectivePhoton;
        d.physicalThresholdSourceIdentity=d.ion.state.definition().formation.sourceIdentity;
      }
      const auto prepared=E::SodiumPhotoionization::Create(source,d,&provider);
      const std::string prefix=std::string(endpoint==0?"quiet":"active")+(convention?"-crosswalk":"");
      Check(prefix+"-source-admission",prepared.ok());if(!prepared.ok())throw std::runtime_error(prepared.message);
      for(unsigned mode=0;mode<4;++mode){
        const std::string group=prefix+(mode==0?"-rest":mode==1?"-moving":mode==2?"-affine":"-rounded");auto input=Input(mode);E::Composition value;
        const auto status=provider.Compose(input,&value);Check(group+"-compose",status.ok());if(!status.ok())throw std::runtime_error(status.message);
        const auto& plan=value.plan;const auto& p=plan.products().front();
        Check(group+"-conservative-plan",plan.valid()&&plan.ledger().integerConservationVerified&&plan.ledger().kinematicConservationVerified);
        Check(group+"-one-ion-plus-bath-electron",plan.products().size()==1&&p.species==1&&plan.ledger().bathTerms.size()==1&&plan.ledger().bathChargeChangeElementary==-1);
        Check(group+"-number-no-timestep-factor",p.representedNumber==2.625&&p.referenceWeight==.4);
        Check(group+"-original-unused-interval",p.interval.startS()==100&&p.interval.durationS()==1&&p.interval.elapsedS()==.375&&p.interval.remainingS()==.625);
        Check(group+"-fresh-daughter-clock",p.clock.fresh()&&p.clock.valid()&&p.clock.remaining()==input.freshIonClock.remaining());
        Check(group+"-frozen-basis-distinct-physical-epoch",value.receipt.frame.geometryBasisEpochEtS==90&&value.receipt.frame.physicalEventEpochEtS==100.375);
        Check(group+"-source-and-policy-identities",value.receipt.sourceChannel.id==d.channelId&&value.receipt.sourceProfileSha256==argv[2]&&plan.eventFrame()==d.plannerAuthorization.eventFrame);
        Record(group,d,input,value);
      }
    }
    {
      const auto rounded=Input(3);const auto& r=rounded.frame.inertialToEventRotation;
      K::Vector3 inverse{},roundtrip{};
      for(unsigned i=0;i<3;++i)for(unsigned j=0;j<3;++j)inverse[i]-=r[3*j+i]*rounded.frame.velocityOffsetEventMS[j];
      for(unsigned i=0;i<3;++i){for(unsigned j=0;j<3;++j)roundtrip[i]+=r[3*i+j]*inverse[j];roundtrip[i]+=rounded.frame.velocityOffsetEventMS[i];}
      Check("recorded-rounded-old-transpose-violates-existing-momentum-bound",std::fabs(4e-26L*roundtrip[1])>1e-38L);
    }
    auto d=Definition(source,P::SolarCase::Quiet,argv[2]);E::SodiumPhotoionization stable;
    if(!E::SodiumPhotoionization::Create(source,d,&stable).ok())throw std::runtime_error("cannot prepare guard fixture");
    auto refuse=[&](const std::string& name,const E::Definition& bad,E::Code code){
      auto retained=stable;auto status=E::SodiumPhotoionization::Create(source,bad,&retained);
      Check(name,status.code==code&&retained.prepared()&&retained.reaction()->identity==stable.reaction()->identity);
    };
    auto bad=d;bad.sourceStateCrosswalkAvailable=false;refuse("missing-crosswalk-preserves-provider",bad,E::Code::Unavailable);
    bad=d;bad.photonEnergyPolicy=static_cast<K::PhotonEnergyPolicy>(9);refuse("unknown-photon-energy-policy-refused",bad,E::Code::InvalidInput);
    bad=d;bad.photonEnergyPolicy=K::PhotonEnergyPolicy::SourceEffectivePhoton;
    refuse("effective-photon-missing-physical-source-refused",bad,E::Code::Unavailable);
    bad.physicalThresholdSourceIdentity="WRONG-SOURCE";
    refuse("effective-photon-wrong-physical-source-refused",bad,E::Code::SourceMismatch);
    bad=d;bad.photonEnergyPolicy=K::PhotonEnergyPolicy::SourceEffectivePhoton;
    bad.ion.state=State("Na+",1,d.ion.state.definition().formation.joules+10*d.electronVoltJ);
    bad.physicalThresholdSourceIdentity=bad.ion.state.definition().formation.sourceIdentity;Catalog(&bad);
    refuse("effective-photon-negative-translation-refused",bad,E::Code::OutsideDomain);
    bad=d;bad.ion.physicalRestMass={};refuse("missing-physical-ion-mass-refused",bad,E::Code::Unavailable);
    bad=d;bad.electron.modelErrorLimits={};refuse("missing-independent-error-budget-refused",bad,E::Code::Unavailable);
    bad=d;bad.ion.state=State("Na+",1,d.ion.state.definition().formation.joules+1e-22L);Catalog(&bad);refuse("wrong-state-threshold-refused",bad,E::Code::SourceMismatch);
    // A caller may legitimately choose a loose arithmetic diagnostic, but
    // that must NEVER grant permission to change the physical source gap.
    // This small mismatch resembles the unresolved D04/ASD scale and would
    // pass the old user-controlled threshold comparison with these bounds.
    bad=d;bad.ion.state=State("Na+",1,d.ion.state.definition().formation.joules+1e-24L);
    bad.numerical.thresholdJ={1,1};Catalog(&bad);
    refuse("loose-numerical-bound-cannot-admit-physical-gap",bad,E::Code::SourceMismatch);
    bad=d;bad.catalog[1].modelMassKg=std::nextafter(bad.catalog[1].modelMassKg,std::numeric_limits<double>::infinity());refuse("actual-catalog-mass-mismatch-refused",bad,E::Code::SourceMismatch);
    bad=d;bad.catalog[2].compiledSpecies=2;refuse("PIC-electron-substitution-refused",bad,E::Code::SourceMismatch);
    bad=d;bad.solar=static_cast<P::SolarCase>(2);refuse("endpoint-interpolation-refused",bad,E::Code::InvalidInput);
    bad=d;bad.solar=static_cast<P::SolarCase>(-1);refuse("omitted-solar-endpoint-refused",bad,E::Code::InvalidInput);
    bad=d;bad.plannerAuthorization.qualified=false;refuse("missing-caller-admission-refused",bad,E::Code::Unavailable);
    bad=d;bad.electron.modelErrorLimits.maximumRelativeRelativisticKineticError=0;refuse("physical-model-error-budget-not-widened",bad,E::Code::OutsideDomain);
    bad=d;bad.electronVoltJ=1e-19L;refuse("wrong-eV-conversion-refused",bad,E::Code::InvalidInput);
    std::vector<P::Channel> molecular;source.ChannelsForParent("H2O",&molecular,&error);bad=d;bad.channelId=molecular.front().id;refuse("unsupported-molecular-family-not-renormalized",bad,E::Code::SourceMismatch);
    auto input=Input(2);E::Composition retained;if(!stable.Compose(input,&retained).ok())throw std::runtime_error("cannot prepare composition guard fixture");
    const auto identity=retained.plan.reactionIdentity();const auto vx=retained.plan.products()[0].velocityMS[0];
    auto failEvent=[&](const std::string& name,const E::EventInput& event,E::Code code){
      const auto status=stable.Compose(event,&retained);Check(name,status.code==code&&retained.plan.reactionIdentity()==identity&&retained.plan.products()[0].velocityMS[0]==vx);
    };
    auto wrong=input;wrong.pointSunVisible=0;failEvent("event-in-shadow-refused-no-output",wrong,E::Code::SourceMismatch);
    wrong=input;wrong.polarOpenUniform=0;failEvent("closed-polar-draw-refused-no-output",wrong,E::Code::InvalidInput);
    wrong=input;wrong.azimuthOpenUniform=1;failEvent("closed-azimuth-draw-refused-no-output",wrong,E::Code::InvalidInput);
    wrong=input;wrong.incomingPhotonDirectionInertial={{0,0,0}};failEvent("missing-photon-ray-refused-no-output",wrong,E::Code::InvalidInput);
    wrong=input;wrong.frame.inertialToEventRotation={{1,0,0,0,1,0,0,0,-1}};failEvent("improper-frame-refused-no-output",wrong,E::Code::InvalidInput);
    wrong=input;wrong.frame.velocityOffsetSupplied=false;failEvent("omitted-state-transform-derivative-refused",wrong,E::Code::InvalidInput);
    wrong=input;wrong.frame.physicalEventEpochEtS=200;failEvent("inconsistent-physical-event-epoch-refused",wrong,E::Code::SourceMismatch);
    wrong=input;wrong.parentAtEvent.internalEnergyJ=1;failEvent("wrong-parent-state-energy-refused-no-output",wrong,E::Code::SourceMismatch);
    wrong=input;wrong.parentAtEvent.velocityMS={{4000000,0,0}};failEvent("outside-classical-parent-domain-refused",wrong,E::Code::OutsideDomain);
    wrong=input;wrong.distanceAu=0;failEvent("invalid-source-distance-refused",wrong,E::Code::OutsideDomain);
    wrong=input;wrong.freshIonClock={};failEvent("missing-daughter-clock-refused",wrong,E::Code::SourceMismatch);
    double rate=123;const auto partial=stable.RateS(1,2,&rate);Check("finite-disc-or-cosine-not-substituted",partial.code==E::Code::InvalidInput&&rate==123);
    const auto dark=stable.RateS(2,0,&rate);Check("point-shadow-hazard-exact-zero",dark.ok()&&rate==0);
    auto rounded=input;rounded.frame.physicalEventEpochEtS=std::nextafter(std::nextafter(100.375,200.0),200.0);
    E::Composition projected;const auto projection=stable.Compose(rounded,&projected);
    Check("supported-absolute-projection-retains-actual-context",projection.ok()&&projected.receipt.frame.physicalEventEpochEtS==rounded.frame.physicalEventEpochEtS&&projected.receipt.intervalProjectedEventEpochEtS==100.375);
    C::StateDefinition collector;const auto state=stable.IonCollectorState(&collector);
    Check("same-source-ion-collector-state",state.ok()&&collector.material.state==d.catalog[1].state&&collector.material.modelMassKg==d.catalog[1].modelMassKg&&collector.formationEnergyJ==static_cast<double>(d.ion.state.definition().formation.joules)&&collector.energyReferenceIdentity==d.energyReferenceIdentity);
    // I22 atomic pilot target selection. The provider must bind each admitted
    // atom to its own complete D04 family and rate, only under the named
    // source-effective photon, and refuse every other atom or combination.
    // The analytical ion gap is the printed D04 threshold + 1e-4 eV: a
    // mathematical perturbation, not a physical ionization energy.
    auto atomic=[&](const std::string& element,const std::string& parentNucleus,const std::string& ionNucleus){
      auto a=Definition(source,P::SolarCase::Quiet,argv[2]);std::vector<P::Channel> rows;std::string failure;
      if(!source.ChannelsForParent(element,&rows,&failure)||rows.size()!=1)throw std::runtime_error("atomic fixture family: "+failure);
      a.channelId=rows.front().id;
      const Real gap=(static_cast<Real>(rows.front().thresholdEnergyEv)+1e-4L)*a.electronVoltJ;
      a.parent.state=NucleusState(element,parentNucleus,0,0);a.ion.state=NucleusState(element+"+",ionNucleus,1,gap);
      a.electron.state=NucleusState("e","",-1,0);
      for(auto* b:{&a.parent,&a.ion,&a.electron})b->physicalRestMass={b->state.referenceMassKg(),"ANALYTICAL-reference-equals-physical-mass-control"};
      Catalog(&a);a.photonEnergyPolicy=K::PhotonEnergyPolicy::SourceEffectivePhoton;
      a.physicalThresholdSourceIdentity=a.ion.state.definition().formation.sourceIdentity;
      return a;
    };
    for(const auto& target:{std::make_pair(std::string("He"),std::string("4He")),std::make_pair(std::string("Ne"),std::string("20Ne")),
                            std::make_pair(std::string("Ar"),std::string("40Ar"))}){
      const auto& element=target.first;auto a=atomic(element,target.second,target.second);E::AtomicPhotoionization provider;
      const auto created=E::AtomicPhotoionization::Create(source,a,&provider);
      const std::string model="EXPERIMENTAL-"+element+"-D04-source-effective-photon-isotropic-ground-ion-v1";
      Check("atomic-"+element+"-source-effective-admitted",created.ok()&&provider.sourceParent()==element&&
          provider.reaction()->identity.find(model)!=std::string::npos);
      double expected=0,actual=-1;std::string failure;
      const bool table=source.TotalDestructionRateS(element,P::SolarCase::Quiet,1,1,&expected,&failure);
      Check("atomic-"+element+"-rate-is-own-D04-family",created.ok()&&table&&provider.RateS(1,1,&actual).ok()&&actual==expected&&actual>0);
      E::Composition composed;const auto status=created.ok()?provider.Compose(Input(0),&composed):created;
      Check("atomic-"+element+"-ground-ion-receipt",status.ok()&&composed.receipt.modelIdentity==model&&
          composed.receipt.scope.find("moving-parent "+element+" to ground "+element+"+")!=std::string::npos&&
          composed.plan.products().size()==1&&composed.plan.ledger().bathChargeChangeElementary==-1&&
          composed.receipt.photonEnergyPolicy=="source-effective-photon");
      auto strict=a;strict.photonEnergyPolicy=K::PhotonEnergyPolicy::SourceEqualsState;strict.physicalThresholdSourceIdentity.clear();
      refuse("atomic-"+element+"-source-equals-state-refused",strict,E::Code::SourceMismatch);
    }
    refuse("atomic-mixed-isotope-ion-refused",atomic("He","4He","20Ne"),E::Code::SourceMismatch);
    refuse("atomic-unsupported-nucleus-refused",atomic("Ar","39K","39K"),E::Code::SourceMismatch);
    {auto wrongChannel=atomic("He","4He","4He");wrongChannel.channelId=d.channelId;
     refuse("atomic-wrong-element-channel-refused",wrongChannel,E::Code::SourceMismatch);}
    Check("atomic-Na-model-and-scope-unchanged",std::string(stable.sourceParent())=="Na"&&retained.receipt.modelIdentity==
        "EXPERIMENTAL-Na-D04-source-mean-isotropic-moving-parent-v1"&&retained.receipt.scope==
        "EXPERIMENTAL source-mean/isotropic moving-parent Na; not validated kinetics");
    checks.flush();metrics.flush();std::cout<<"SUMMARY checks="<<total<<" failed="<<failed<<'\n';return failed?1:0;
  }catch(const std::exception& e){std::cerr<<"ERROR "<<e.what()<<'\n';return 2;}
}
