#ifndef AMPS_CHEMICAL_DATABASE_RESERVOIR_PRODUCTS_H
#define AMPS_CHEMICAL_DATABASE_RESERVOIR_PRODUCTS_H

#include "ReducedModels.h"
#include <cstdint>
#include <map>
#include <string>

namespace AMPS { namespace ChemicalDatabase { namespace ReservoirProducts {
using Real=long double;
using Vector3=ReducedModels::Vector3;
// This service performs no PIC operation, rate selection, RNG, source download,
// permission approval or model activation. All moments use ONE named inertial
// frame and SI kg,m,s,J,C. State energies share the named reference; kinetic
// masses exclude any binding contribution already represented by state Q.
struct Nucleus {std::string identity;std::int64_t count=0;};
struct State {
  std::string identity,state,sourceIdentity,energyReferenceIdentity;
  std::vector<Nucleus> nuclei;
  std::int64_t chargeElementary=0;
  Real kineticMassKg=-1,physicalRestMassKg=-1,stateEnergyJ=std::numeric_limits<Real>::quiet_NaN();
};
struct NumericalBounds {
  Real massAbsoluteKg=-1,massRelative=-1,timeProjectionAbsoluteS=-1;
  ReducedModels::ConservationBounds conservation{-1,-1,-1,-1,-1};
};
struct Evidence {
  // Caller authenticates these pinned reviewed documents before constructing
  // the request. Nonempty names/hashes alone are NOT a scientific certificate.
  std::string modelIdentity,sourceProfileIdentity,sourceProfileSha256;
  std::string stateCrosswalkIdentity,stateCrosswalkSha256;
  std::string domainEvidenceIdentity,modelErrorEvidenceIdentity,modelErrorEvidenceSha256;
  std::string frameIdentity,timeScaleIdentity,energyReferenceIdentity,observableIdentity;
  Real maximumHeavyVelocityErrorMS=-1,maximumHeavyEnergyErrorJ=-1;
  Real maximumRelativeMassError=-1,maximumClassicalHeavySpeedMS=-1;
};
struct EventContext {
  Vector3 positionM{{0,0,0}},parentVelocityMS{{0,0,0}};
  Real physicalStartS=std::numeric_limits<Real>::quiet_NaN(),elapsedS=-1,remainingS=-1;
  Real physicalEventS=std::numeric_limits<Real>::quiet_NaN(),representedNumber=-1;
  std::string frameIdentity,timeScaleIdentity;
};
struct Moment {
  State sourceState; // exact source mass/state/energy/provenance snapshot
  std::string speciesIdentity,stateIdentity;
  std::int64_t number=0;
  std::vector<Nucleus> nuclei;
  Real massKg=0,chargeC=0,kineticEnergyJ=0,stateEnergyJ=0;
  Vector3 momentumKgMS{{0,0,0}};
};
struct ErrorReport {
  // EI: RMS omitted recoil bound, not an all-events velocity-support bound.
  // CEX: pointwise bound conditional on the EXPLICIT supported COM angle cap.
  std::string interpretation;
  Real heavyVelocityUpperMS=0,heavyKineticMeanAbsoluteUpperJ=0,relativeMassError=0;
};
struct Result {
  std::string scope="REDUCED HEAVY-PRODUCT MOMENTS ONLY; NOT FULL DAUGHTER KINETICS OR PHYSICAL QUALIFICATION";
  Evidence evidence;EventContext event;
  Moment parent,trackedProduct,incomingBath,outgoingBath;
  Vector3 trackedVelocityMS{{0,0,0}},outgoingBathVelocityMS{{0,0,0}};
  Real stateQJ=0,bathNumberChange=0,bathMassChangeKg=0,bathChargeChangeC=0,bathEnergyChangeJ=0;
  Vector3 bathMomentumChangeKgMS{{0,0,0}},momentumResidualKgMS{{0,0,0}};
  Real energyResidualJ=0;
  // Packet moments multiply the SAME reacted N. No daughter timestep factor.
  Vector3 representedBathMomentumChangeKgMS{{0,0,0}};
  Real representedBathEnergyChangeJ=0,representedBathChargeChangeC=0;
  ErrorReport error;
  // Preserve the actual moment prescription as well as its pinned profile.
  std::string incomingMomentSourceIdentity,partialChannelIdentity,angleSupportEvidenceIdentity;
  Real incomingThermalEnergyJ=0,incomingLogarithmicRateDerivative=0,incomingConditionalMeanEnergyJ=0;
  Real incomingRelativeKineticEnergyJ=0,sourceMaximumCOMDeflectionRad=0;
};
struct MaxwellianInput {
  // The identity follows differentiation of the Maxwellian rate integral.
  // A fit-derived derivative yields a MODEL conditional mean, not a measurement.
  std::string derivativeSourceIdentity,temperatureDomainIdentity;
  Real thermalEnergyJ=-1,minimumThermalEnergyJ=-1,maximumThermalEnergyJ=-1;
  Real logarithmicRateDerivative=std::numeric_limits<Real>::quiet_NaN();
  Real thresholdJ=-1;
  Vector3 bulkVelocityParentRestMS{{std::numeric_limits<Real>::quiet_NaN(),std::numeric_limits<Real>::quiet_NaN(),std::numeric_limits<Real>::quiet_NaN()}};
  bool isotropicUnshiftedParentRest=false,crossSectionTemperatureIndependent=false;
  bool collisionWeightRotationallyInvariant=false; // needed for zero first moment, not just energy identity
};
inline Status MaxwellianCollisionMeanEnergy(const MaxwellianInput& in,Real* output) {
  if(!output)return {Code::InvalidInput,"null Maxwellian moment output"};
  if(in.derivativeSourceIdentity.empty()||in.temperatureDomainIdentity.empty()||
     !in.isotropicUnshiftedParentRest||!in.crossSectionTemperatureIndependent||
     std::isnan(in.logarithmicRateDerivative)||!ReducedModels::Finite(in.bulkVelocityParentRestMS))
    return {Code::Unavailable,"source derivative/domain and unshifted isotropic T-independent-sigma assumptions required; drift is unsupported"};
  if(!ReducedModels::Finite(in.logarithmicRateDerivative))return {Code::InvalidInput,"nonfinite supplied rate derivative"};
  if(in.bulkVelocityParentRestMS!=Vector3{{0,0,0}})return {Code::OutsideDomain,"actual nonzero parent-rest Maxwellian drift is unsupported"};
  if(!ReducedModels::Finite(in.thermalEnergyJ)||!ReducedModels::Finite(in.minimumThermalEnergyJ)||
     !ReducedModels::Finite(in.maximumThermalEnergyJ)||!ReducedModels::Finite(in.thresholdJ)||
     in.thermalEnergyJ<=0||in.minimumThermalEnergyJ<=0||in.maximumThermalEnergyJ<in.minimumThermalEnergyJ||in.thresholdJ<0)
    return {Code::InvalidInput,"invalid explicit Maxwellian thermal-energy/threshold domain [J]"};
  if(in.thermalEnergyJ<in.minimumThermalEnergyJ||in.thermalEnergyJ>in.maximumThermalEnergyJ)
    return {Code::OutsideDomain,"temperature outside derivative source domain"};
  const Real mean=in.thermalEnergyJ*(1.5L+in.logarithmicRateDerivative);
  if(!ReducedModels::Finite(mean)||mean<=0)return {Code::NumericalFailure,"collision-conditioned mean energy is not representable/positive"};
  if(mean<in.thresholdJ)return {Code::OutsideDomain,"rate derivative implies mean collision energy below selected threshold"};
  *output=mean;return {};
}
namespace detail {
inline bool Text(const std::string& s){return !s.empty();}
inline bool Hash(const std::string& s){return s.size()==64&&s.find_first_not_of("0123456789abcdef")==std::string::npos;}
inline Real Norm(const Vector3& v){return std::hypot(v[0],v[1],v[2]);}
inline Real Dot(const Vector3& a,const Vector3& b){return ReducedModels::Sum3(a[0]*b[0],a[1]*b[1],a[2]*b[2]);}
inline Status StateValid(const State& s,const Evidence& e) {
  if(!Text(s.identity)||!Text(s.state)||!Text(s.sourceIdentity)||!Text(s.energyReferenceIdentity)||
     !ReducedModels::Finite(s.stateEnergyJ)||!ReducedModels::Finite(s.physicalRestMassKg)||s.physicalRestMassKg<=0)
    return {Code::Unavailable,"explicit source/state/physical mass/common energy reference unavailable"};
  if(s.energyReferenceIdentity!=e.energyReferenceIdentity)return {Code::InvalidInput,"state energy reference contradicts selected profile"};
  if(!ReducedModels::Finite(s.kineticMassKg)||s.kineticMassKg<=0)return {Code::InvalidInput,"positive source kinetic mass required"};
  std::map<std::string,bool> seen;
  for(const auto& n:s.nuclei)if(!Text(n.identity)||n.count<=0||!seen.emplace(n.identity,true).second)
    return {Code::InvalidInput,"nuclei must be distinct explicit positive integer isotope identities"};
  return {};
}
inline Status Common(const Evidence& e,const EventContext& event,const NumericalBounds& b) {
  if(!Text(e.modelIdentity)||!Text(e.sourceProfileIdentity)||e.sourceProfileSha256.empty()||
     !Text(e.stateCrosswalkIdentity)||e.stateCrosswalkSha256.empty()||!Text(e.domainEvidenceIdentity)||
     !Text(e.modelErrorEvidenceIdentity)||e.modelErrorEvidenceSha256.empty()||!Text(e.energyReferenceIdentity)||
     !Text(e.observableIdentity)||!Text(e.frameIdentity)||!Text(e.timeScaleIdentity)||
     e.maximumHeavyVelocityErrorMS<0||e.maximumHeavyEnergyErrorJ<0||e.maximumRelativeMassError<0||e.maximumClassicalHeavySpeedMS<=0)
    return {Code::Unavailable,"reviewed source/state/domain/observable model-error profile and explicit limits required"};
  if(!Hash(e.sourceProfileSha256)||!Hash(e.stateCrosswalkSha256)||!Hash(e.modelErrorEvidenceSha256))
    return {Code::InvalidInput,"source/crosswalk/model-error SHA256 syntax is corrupt"};
  for(Real x:{e.maximumHeavyVelocityErrorMS,e.maximumHeavyEnergyErrorJ,e.maximumRelativeMassError,e.maximumClassicalHeavySpeedMS})
    if(!ReducedModels::Finite(x))return {Code::InvalidInput,"nonfinite model error/domain limit"};
  if(!ReducedModels::Valid(b.conservation)||!ReducedModels::Finite(b.massAbsoluteKg)||!ReducedModels::Finite(b.massRelative)||
     !ReducedModels::Finite(b.timeProjectionAbsoluteS)||b.massAbsoluteKg<0||b.massRelative<0||b.timeProjectionAbsoluteS<0)
    return {Code::InvalidInput,"explicit nonnegative NUMERICAL conservation/time bounds required"};
  if(event.frameIdentity!=e.frameIdentity||event.timeScaleIdentity!=e.timeScaleIdentity||
     !ReducedModels::Finite(event.positionM)||!ReducedModels::Finite(event.parentVelocityMS)||
     !ReducedModels::Finite(event.physicalStartS)||!ReducedModels::Finite(event.elapsedS)||!ReducedModels::Finite(event.remainingS)||
     !ReducedModels::Finite(event.physicalEventS)||!ReducedModels::Finite(event.representedNumber)||
     event.elapsedS<0||event.remainingS<0||event.representedNumber<=0)
    return {Code::InvalidInput,"actual same-frame event state/time/positive represented number required"};
  if(std::fabs(ReducedModels::Sum3(event.physicalStartS,event.elapsedS,-event.physicalEventS))>b.timeProjectionAbsoluteS)
    return {Code::InvalidInput,"physical event epoch contradicts supplied interval projection"};
  if(Norm(event.parentVelocityMS)>e.maximumClassicalHeavySpeedMS)return {Code::OutsideDomain,"actual heavy speed outside declared classical domain"};
  return {};
}
inline Moment Make(const State& s,std::int64_t count,const Vector3& momentum,Real kinetic,Real elementaryChargeC) {
  Moment m;m.sourceState=s;m.speciesIdentity=s.identity;m.stateIdentity=s.state;m.number=count;m.nuclei=s.nuclei;
  for(auto& n:m.nuclei)n.count*=count;
  m.massKg=count*s.kineticMassKg;m.chargeC=count*static_cast<Real>(s.chargeElementary)*elementaryChargeC;
  m.stateEnergyJ=count*s.stateEnergyJ;m.momentumKgMS=momentum;m.kineticEnergyJ=kinetic;return m;
}
inline Moment Moving(const State& s,const Vector3& v,Real elementaryChargeC) {
  Vector3 p{};for(int i=0;i<3;++i)p[i]=s.kineticMassKg*v[i];
  return Make(s,1,p,ReducedModels::Kinetic(p,s.kineticMassKg),elementaryChargeC);
}
inline bool Add(std::int64_t* total,std::int64_t value) {
  if((value>0&&*total>std::numeric_limits<std::int64_t>::max()-value)||
     (value<0&&*total<std::numeric_limits<std::int64_t>::min()-value))return false;
  *total+=value;return true;
}
inline Status Topology(const State& p,const State& in,const State& h,const State& out,std::int64_t outCount) {
  std::map<std::string,std::int64_t> left,right;
  for(const State* s:{&p,&in})for(const auto& n:s->nuclei)if(!Add(&left[n.identity],n.count))return {Code::NumericalFailure,"integer nuclei overflow"};
  for(const auto& n:h.nuclei)if(!Add(&right[n.identity],n.count))return {Code::NumericalFailure,"integer nuclei overflow"};
  for(const auto& n:out.nuclei) {
    if(n.count>std::numeric_limits<std::int64_t>::max()/outCount||!Add(&right[n.identity],n.count*outCount))return {Code::NumericalFailure,"integer product nuclei overflow"};
  }
  if(left!=right)return {Code::InvalidInput,"explicit isotope nuclei do not close"};
  std::int64_t a=p.chargeElementary,c=h.chargeElementary;
  if(!Add(&a,in.chargeElementary))return {Code::NumericalFailure,"integer charge overflow"};
  for(std::int64_t i=0;i<outCount;++i)if(!Add(&c,out.chargeElementary))return {Code::NumericalFailure,"integer charge overflow"};
  if(a!=c)return {Code::InvalidInput,"integer state charges do not close"};
  return {};
}
inline bool SameNuclei(const State& a,const State& b) {
  std::map<std::string,std::int64_t> left,right;
  for(const auto& n:a.nuclei)left.emplace(n.identity,n.count);
  for(const auto& n:b.nuclei)right.emplace(n.identity,n.count);
  return left==right;
}
inline Status MaterialMass(const State& parent,const State& incoming,const State& heavy,const State& outgoing,
    std::int64_t outgoingCount,const NumericalBounds& b) {
  const Real before=parent.kineticMassKg+incoming.kineticMassKg;
  const Real after=heavy.kineticMassKg+outgoingCount*outgoing.kineticMassKg;
  const Real limit=b.massAbsoluteKg+b.massRelative*std::max(before,after);
  if(!ReducedModels::Finite(before)||!ReducedModels::Finite(after)||!ReducedModels::Finite(limit))
    return {Code::NumericalFailure,"reference material mass sum/allowance overflow"};
  if(std::fabs(before-after)>limit)return {Code::InvalidInput,"reference material mass convention does not close"};
  return {};
}
inline Real EnergyResidual(const Moment& parent,const Moment& incoming,const Moment& heavy,const Moment& outgoing) {
  // Sum the EIGHT signed terms before comparing with zero. Subtracting two
  // separately rounded large total energies could hide a small supplied Q.
  Real sum=0,correction=0;
  for(Real value:{heavy.kineticEnergyJ,heavy.stateEnergyJ,outgoing.kineticEnergyJ,outgoing.stateEnergyJ,
                 -parent.kineticEnergyJ,-parent.stateEnergyJ,-incoming.kineticEnergyJ,-incoming.stateEnergyJ}) {
    const Real next=sum+value;
    correction+=(std::fabs(sum)>=std::fabs(value))?(sum-next)+value:(value-next)+sum;
    sum=next;
  }
  return sum+correction;
}
inline Status Finish(Result* r,const Evidence& e,const NumericalBounds& b,Real heavyMass) {
  for(const Moment* m:{&r->parent,&r->trackedProduct,&r->incomingBath,&r->outgoingBath})
    if(!ReducedModels::Finite(m->massKg)||!ReducedModels::Finite(m->chargeC)||!ReducedModels::Finite(m->kineticEnergyJ)||
       !ReducedModels::Finite(m->stateEnergyJ)||!ReducedModels::Finite(m->momentumKgMS)||m->massKg<=0||m->kineticEnergyJ<0)
      return {Code::NumericalFailure,"nonfinite/unrepresentable physical moment"};
  const Real beforeMass=r->parent.massKg+r->incomingBath.massKg,afterMass=r->trackedProduct.massKg+r->outgoingBath.massKg;
  const Real allowedMass=b.massAbsoluteKg+b.massRelative*std::max(beforeMass,afterMass);
  if(!ReducedModels::Finite(allowedMass)||!ReducedModels::Finite(beforeMass)||!ReducedModels::Finite(afterMass)||std::fabs(beforeMass-afterMass)>allowedMass)
    return {Code::InvalidInput,"reference material mass convention does not close within explicit arithmetic bounds"};
  r->bathNumberChange=r->outgoingBath.number-r->incomingBath.number;
  r->bathMassChangeKg=r->outgoingBath.massKg-r->incomingBath.massKg;
  r->bathChargeChangeC=r->outgoingBath.chargeC-r->incomingBath.chargeC;
  r->bathEnergyChangeJ=ReducedModels::Sum4(r->outgoingBath.kineticEnergyJ,r->outgoingBath.stateEnergyJ,-r->incomingBath.kineticEnergyJ,-r->incomingBath.stateEnergyJ);
  Real pScale=0;
  for(int i=0;i<3;++i) {
    r->bathMomentumChangeKgMS[i]=r->outgoingBath.momentumKgMS[i]-r->incomingBath.momentumKgMS[i];
    r->momentumResidualKgMS[i]=ReducedModels::Sum4(r->trackedProduct.momentumKgMS[i],r->outgoingBath.momentumKgMS[i],-r->parent.momentumKgMS[i],-r->incomingBath.momentumKgMS[i]);
    pScale=std::max({pScale,std::fabs(r->parent.momentumKgMS[i]),std::fabs(r->incomingBath.momentumKgMS[i]),std::fabs(r->trackedProduct.momentumKgMS[i]),std::fabs(r->outgoingBath.momentumKgMS[i])});
  }
  const Real beforeE=ReducedModels::Sum4(r->parent.kineticEnergyJ,r->parent.stateEnergyJ,r->incomingBath.kineticEnergyJ,r->incomingBath.stateEnergyJ);
  const Real afterE=ReducedModels::Sum4(r->trackedProduct.kineticEnergyJ,r->trackedProduct.stateEnergyJ,r->outgoingBath.kineticEnergyJ,r->outgoingBath.stateEnergyJ);
  r->energyResidualJ=EnergyResidual(r->parent,r->incomingBath,r->trackedProduct,r->outgoingBath);
  const Real allowedP=b.conservation.momentumAbsoluteKgMS+b.conservation.momentumRelative*pScale;
  const Real allowedE=b.conservation.energyAbsoluteJ+b.conservation.energyRelative*std::max(std::fabs(beforeE),std::fabs(afterE));
  if(!ReducedModels::Finite(allowedP)||!ReducedModels::Finite(allowedE)||!ReducedModels::Finite(beforeE)||!ReducedModels::Finite(afterE)||!ReducedModels::Finite(r->bathEnergyChangeJ)||
     !ReducedModels::Finite(r->bathMomentumChangeKgMS)||!ReducedModels::Finite(r->momentumResidualKgMS)||
     Norm(r->momentumResidualKgMS)>allowedP||
     std::fabs(r->energyResidualJ)>allowedE)
    return {Code::NumericalFailure,"computed reduced event violates numerical moment conservation"};
  if(!ReducedModels::Finite(r->error.heavyVelocityUpperMS)||!ReducedModels::Finite(r->error.heavyKineticMeanAbsoluteUpperJ)||
     !ReducedModels::Finite(r->error.relativeMassError))return {Code::NumericalFailure,"model discrepancy bound overflow"};
  if(r->error.heavyVelocityUpperMS>e.maximumHeavyVelocityErrorMS||r->error.heavyKineticMeanAbsoluteUpperJ>e.maximumHeavyEnergyErrorJ||r->error.relativeMassError>e.maximumRelativeMassError)
    return {Code::OutsideDomain,"derived local heavy-observable discrepancy exceeds declared model budget; no tolerance relaxation"};
  if(Norm(r->trackedVelocityMS)>e.maximumClassicalHeavySpeedMS||Norm(r->outgoingBathVelocityMS)>e.maximumClassicalHeavySpeedMS)
    return {Code::OutsideDomain,"output speed outside declared classical heavy-product domain"};
  for(int i=0;i<3;++i)r->representedBathMomentumChangeKgMS[i]=r->event.representedNumber*r->bathMomentumChangeKgMS[i];
  r->representedBathEnergyChangeJ=r->event.representedNumber*r->bathEnergyChangeJ;
  r->representedBathChargeChangeC=r->event.representedNumber*r->bathChargeChangeC;
  for(int i=0;i<3;++i)if(r->bathMomentumChangeKgMS[i]!=0&&r->representedBathMomentumChangeKgMS[i]==0)
    return {Code::NumericalFailure,"positive represented number erased nonzero momentum; increase arithmetic precision"};
  if((r->bathEnergyChangeJ!=0&&r->representedBathEnergyChangeJ==0)||(r->bathChargeChangeC!=0&&r->representedBathChargeChangeC==0))
    return {Code::NumericalFailure,"positive represented number erased nonzero energy/charge; increase arithmetic precision"};
  if(!ReducedModels::Finite(r->representedBathMomentumChangeKgMS)||!ReducedModels::Finite(r->representedBathEnergyChangeJ)||!ReducedModels::Finite(r->representedBathChargeChangeC)||heavyMass<=0)
    return {Code::NumericalFailure,"represented event moments overflow"};
  return {};
}
} // detail
struct ElectronImpactInput {
  Evidence evidence;NumericalBounds numerical;EventContext event;
  State parent,ion,electron;MaxwellianInput incoming;
  Real elementaryChargeC=-1;
  // A reviewed source-to-state threshold bound, NOT an arithmetic tolerance.
  Real sourceThresholdCrosswalkMaximumDifferenceJ=-1;
};
inline Status ElectronImpactInheritedVelocity(const ElectronImpactInput& in,Result* output) {
  try {
  if(!output)return {Code::InvalidInput,"null electron-impact output"};
  auto status=detail::Common(in.evidence,in.event,in.numerical);if(!status.ok())return status;
  for(const State* state:{&in.parent,&in.ion,&in.electron}){status=detail::StateValid(*state,in.evidence);if(!status.ok())return status;}
  if(!in.electron.nuclei.empty()||in.electron.chargeElementary!=-1||!in.parent.nuclei.size()||
     !ReducedModels::Finite(in.elementaryChargeC)||in.elementaryChargeC<=0)
    return {Code::InvalidInput,"explicit electron state/elementary charge and heavy nucleus identities required"};
  status=detail::Topology(in.parent,in.electron,in.ion,in.electron,2);if(!status.ok())return status;
  status=detail::MaterialMass(in.parent,in.electron,in.ion,in.electron,2,in.numerical);if(!status.ok())return status;
  Real mean=0;status=MaxwellianCollisionMeanEnergy(in.incoming,&mean);if(!status.ok())return status;
  if(!in.incoming.collisionWeightRotationallyInvariant)return {Code::Unavailable,"zero incoming electron momentum requires explicit rotationally invariant collision weighting"};
  const Real q=ReducedModels::Sum4(in.parent.stateEnergyJ,in.electron.stateEnergyJ,-in.ion.stateEnergyJ,-2*in.electron.stateEnergyJ);
  if(!ReducedModels::Finite(in.sourceThresholdCrosswalkMaximumDifferenceJ)||in.sourceThresholdCrosswalkMaximumDifferenceJ<0)
    return {Code::Unavailable,"reviewed source threshold/state-Q discrepancy bound unavailable"};
  if(q>0||std::fabs(-q-in.incoming.thresholdJ)>in.sourceThresholdCrosswalkMaximumDifferenceJ)
    return {Code::OutsideDomain,"selected state-Q and rate threshold violate explicit source crosswalk"};
  if(!ReducedModels::Finite(q))return {Code::NumericalFailure,"state-Q sum overflow"};
  const Real outgoingRestEnergy=mean+q;
  if(!ReducedModels::Finite(outgoingRestEnergy)||outgoingRestEnergy<0)return {Code::OutsideDomain,"negative outgoing aggregate electron kinetic energy"};
  Result r;r.evidence=in.evidence;r.event=in.event;r.stateQJ=q;r.parent=detail::Moving(in.parent,in.event.parentVelocityMS,in.elementaryChargeC);
  r.incomingMomentSourceIdentity=in.incoming.derivativeSourceIdentity;
  r.incomingThermalEnergyJ=in.incoming.thermalEnergyJ;r.incomingLogarithmicRateDerivative=in.incoming.logarithmicRateDerivative;
  r.incomingConditionalMeanEnergyJ=mean;
  r.trackedVelocityMS=in.event.parentVelocityMS;
  r.trackedProduct=detail::Moving(in.ion,r.trackedVelocityMS,in.elementaryChargeC);
  Vector3 incomingP{};for(int i=0;i<3;++i)incomingP[i]=in.electron.kineticMassKg*in.event.parentVelocityMS[i];
  const Real incomingE=mean+ReducedModels::Kinetic(incomingP,in.electron.kineticMassKg);
  r.incomingBath=detail::Make(in.electron,1,incomingP,incomingE,in.elementaryChargeC);
  ReducedModels::InheritedState inherited;
  status=ReducedModels::InheritParentVelocity(in.parent.kineticMassKg,in.ion.kineticMassKg,
      in.parent.chargeElementary*in.elementaryChargeC,in.ion.chargeElementary*in.elementaryChargeC,
      -in.elementaryChargeC,in.event.parentVelocityMS,incomingP,incomingE,q,&inherited);
  if(!status.ok())return status;
  r.outgoingBath=detail::Make(in.electron,2,inherited.signedReservoirMomentumKgMS,inherited.signedReservoirEnergyJ,in.elementaryChargeC);
  // A number2 aggregate must have at least its bulk kinetic energy. We do not
  // select either individual electron direction/energy or infer a joint PDF.
  if(r.outgoingBath.kineticEnergyJ<ReducedModels::Kinetic(r.outgoingBath.momentumKgMS,r.outgoingBath.massKg))
    return {Code::OutsideDomain,"outgoing electron aggregate second moment is not realizable"};
  const Real impulse=std::sqrt(2*in.electron.kineticMassKg*mean)+std::sqrt(4*in.electron.kineticMassKg*outgoingRestEnergy);
  r.error.heavyVelocityUpperMS=impulse/in.ion.kineticMassKg;
  r.error.heavyKineticMeanAbsoluteUpperJ=detail::Norm(in.event.parentVelocityMS)*impulse+impulse*impulse/(2*in.ion.kineticMassKg);
  for(const State* state:{&in.parent,&in.ion,&in.electron})r.error.relativeMassError=std::max(r.error.relativeMassError,std::fabs(state->physicalRestMassKg-state->kineticMassKg)/state->physicalRestMassKg);
  r.error.interpretation="Upper RMS omitted electron-recoil speed and mean absolute heavy kinetic discrepancy; no downstream trajectory or all-event bound";
  status=detail::Finish(&r,in.evidence,in.numerical,in.ion.kineticMassKg);if(!status.ok())return status;
  *output=std::move(r);return {};
  } catch(const std::bad_alloc&) {
    return {Code::ResourceLimit,"reduced event allocation failed; reduce metadata/event batch or provide more memory"};
  }
}
struct ForwardChargeExchangeInput {
  Evidence evidence;NumericalBounds numerical;EventContext event;
  State parent,plasmaIon,heavyIon,fastNeutral;
  Vector3 plasmaVelocityMS{{0,0,0}};
  Real elementaryChargeC=-1;
  Real sourceSupportedMaximumCOMDeflectionRad=-1;
  std::string partialChannelIdentity,angleSupportEvidenceIdentity,sourceEnergyDomainMappingIdentity;
  // Caller converts the exact native domain to incoming COM-relative J using
  // the declared masses/convention; no keV/u or lab/COM assumption is made.
  Real minimumIncomingRelativeKineticEnergyJ=-1,maximumIncomingRelativeKineticEnergyJ=-1;
};
inline Status ForwardPartialChargeExchange(const ForwardChargeExchangeInput& in,Result* output) {
  try {
  if(!output)return {Code::InvalidInput,"null charge-exchange output"};
  auto status=detail::Common(in.evidence,in.event,in.numerical);if(!status.ok())return status;
  for(const State* state:{&in.parent,&in.plasmaIon,&in.heavyIon,&in.fastNeutral}){status=detail::StateValid(*state,in.evidence);if(!status.ok())return status;}
  if(in.partialChannelIdentity.empty()||in.angleSupportEvidenceIdentity.empty()||in.sourceEnergyDomainMappingIdentity.empty()||
     in.sourceSupportedMaximumCOMDeflectionRad<0||in.minimumIncomingRelativeKineticEnergyJ<0||in.maximumIncomingRelativeKineticEnergyJ<0)
    return {Code::Unavailable,"explicit partial-channel and source-supported angular error evidence required"};
  if(!ReducedModels::Finite(in.plasmaVelocityMS)||!ReducedModels::Finite(in.elementaryChargeC)||in.elementaryChargeC<=0||
     !ReducedModels::Finite(in.sourceSupportedMaximumCOMDeflectionRad)||in.sourceSupportedMaximumCOMDeflectionRad>std::acos(-1.0L)||
     !ReducedModels::Finite(in.minimumIncomingRelativeKineticEnergyJ)||!ReducedModels::Finite(in.maximumIncomingRelativeKineticEnergyJ)||
     in.maximumIncomingRelativeKineticEnergyJ<in.minimumIncomingRelativeKineticEnergyJ)
    return {Code::InvalidInput,"invalid source beam/elementary charge/COM angular domain"};
  if(detail::Norm(in.plasmaVelocityMS)>in.evidence.maximumClassicalHeavySpeedMS)return {Code::OutsideDomain,"incoming plasma-ion speed outside declared classical domain"};
  if(in.parent.nuclei.empty()||in.plasmaIon.nuclei.empty()||in.parent.chargeElementary!=0||in.heavyIon.chargeElementary!=1||in.plasmaIon.chargeElementary!=1||in.fastNeutral.chargeElementary!=0||
     !detail::SameNuclei(in.parent,in.heavyIon)||!detail::SameNuclei(in.plasmaIon,in.fastNeutral))
    return {Code::OutsideDomain,"forward neutral-target/singly-charged-projectile CEX requires unchanged isotope nuclei per carrier and one transferred electron"};
  status=detail::Topology(in.parent,in.plasmaIon,in.heavyIon,in.fastNeutral,1);if(!status.ok())return status;
  status=detail::MaterialMass(in.parent,in.plasmaIon,in.heavyIon,in.fastNeutral,1,in.numerical);if(!status.ok())return status;
  Vector3 direction{};for(int i=0;i<3;++i)direction[i]=in.event.parentVelocityMS[i]-in.plasmaVelocityMS[i];
  const Real norm=detail::Norm(direction);
  if(!ReducedModels::Finite(norm)||norm==0)return {Code::OutsideDomain,"zero/undefined incoming relative direction"};
  const Real inputTotalMass=in.parent.kineticMassKg+in.plasmaIon.kineticMassKg;
  const Real inputReduced=(in.parent.kineticMassKg/inputTotalMass)*in.plasmaIon.kineticMassKg;
  const Real incomingRelativeEnergy=(inputReduced*norm/2)*norm;
  if(!ReducedModels::Finite(inputTotalMass)||!ReducedModels::Finite(incomingRelativeEnergy)||incomingRelativeEnergy<=0)
    return {Code::NumericalFailure,"incoming relative energy is not representable"};
  if(incomingRelativeEnergy<in.minimumIncomingRelativeKineticEnergyJ||incomingRelativeEnergy>in.maximumIncomingRelativeKineticEnergyJ)
    return {Code::OutsideDomain,"actual incoming relative energy outside explicit source-mapped partial domain"};
  for(auto& x:direction)x/=norm;
  Result r;r.evidence=in.evidence;r.event=in.event;r.parent=detail::Moving(in.parent,in.event.parentVelocityMS,in.elementaryChargeC);
  r.incomingBath=detail::Moving(in.plasmaIon,in.plasmaVelocityMS,in.elementaryChargeC);
  r.partialChannelIdentity=in.partialChannelIdentity;r.angleSupportEvidenceIdentity=in.angleSupportEvidenceIdentity;
  r.incomingRelativeKineticEnergyJ=incomingRelativeEnergy;r.sourceMaximumCOMDeflectionRad=in.sourceSupportedMaximumCOMDeflectionRad;
  r.stateQJ=ReducedModels::Sum4(in.parent.stateEnergyJ,in.plasmaIon.stateEnergyJ,-in.heavyIon.stateEnergyJ,-in.fastNeutral.stateEnergyJ);
  if(!ReducedModels::Finite(r.stateQJ))return {Code::NumericalFailure,"state-Q sum overflow"};
  Vector3 total{};for(int i=0;i<3;++i)total[i]=r.parent.momentumKgMS[i]+r.incomingBath.momentumKgMS[i];
  ReducedModels::TwoBodyState products;
  status=ReducedModels::TwoBodyProducts(total,r.parent.kineticEnergyJ+r.incomingBath.kineticEnergyJ,r.stateQJ,
      in.heavyIon.kineticMassKg,in.fastNeutral.kineticMassKg,direction,in.numerical.conservation,&products);
  if(!status.ok())return status;
  r.trackedVelocityMS=products.velocityMS[0];r.outgoingBathVelocityMS=products.velocityMS[1];
  r.trackedProduct=detail::Make(in.heavyIon,1,products.momentumKgMS[0],products.kineticEnergyJ[0],in.elementaryChargeC);
  r.outgoingBath=detail::Make(in.fastNeutral,1,products.momentumKgMS[1],products.kineticEnergyJ[1],in.elementaryChargeC);
  const Real totalMass=in.heavyIon.kineticMassKg+in.fastNeutral.kineticMassKg;
  const Real mu=(in.heavyIon.kineticMassKg/totalMass)*in.fastNeutral.kineticMassKg;
  const Real p=std::sqrt(2*mu*products.centreOfMassKineticEnergyJ);
  const Real impulse=2*p*std::sin(in.sourceSupportedMaximumCOMDeflectionRad/2);
  r.error.heavyVelocityUpperMS=impulse/in.heavyIon.kineticMassKg;
  r.error.heavyKineticMeanAbsoluteUpperJ=detail::Norm(r.trackedVelocityMS)*impulse+impulse*impulse/(2*in.heavyIon.kineticMassKg);
  for(const State* state:{&in.parent,&in.plasmaIon,&in.heavyIon,&in.fastNeutral})r.error.relativeMassError=std::max(r.error.relativeMassError,std::fabs(state->physicalRestMassKg-state->kineticMassKg)/state->physicalRestMassKg);
  r.error.interpretation="Pointwise forward-scattering discrepancy conditional on explicit supported maximum COM angle; no full differential law";
  status=detail::Finish(&r,in.evidence,in.numerical,in.heavyIon.kineticMassKg);if(!status.ok())return status;
  *output=std::move(r);return {};
  } catch(const std::bad_alloc&) {
    return {Code::ResourceLimit,"reduced event allocation failed; reduce metadata/event batch or provide more memory"};
  }
}
} } }
#endif
