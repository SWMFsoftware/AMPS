#include "PhysicalReactionEvents.h"
#include "ExactArithmetic.h"
#include "../MoonTransportClock.h" // pure arithmetic epoch-projection envelope

#include <algorithm>
#include <cmath>
#include <exception>
#include <limits>
#include <map>
#include <new>
#include <utility>

namespace Moon { namespace Species { namespace PhysicalReactionEvents {
namespace {
// SI defines one electronvolt by the exact elementary charge. This is a
// dimensional definition, not a fitted atomic threshold or a numerical gate.
constexpr Real ExactElectronVoltJ=1.602176634e-19L;
const char* const Model="EXPERIMENTAL-Na-D04-source-mean-isotropic-moving-parent-v1";
Status Error(Code code,const std::string& message){return {code,message};}
bool Finite(Real x){return std::isfinite(x);}
bool Positive(Real x){return Finite(x)&&x>0;}
bool FiniteVector(const Vector3& v){return Finite(v[0])&&Finite(v[1])&&Finite(v[2]);}
Real Norm(const Vector3& v){return std::hypot(std::hypot(v[0],v[1]),v[2]);}
Real Dot(const Vector3& a,const Vector3& b){return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];}
bool BoundValid(const K::Bound& x){return Finite(x.absolute)&&Finite(x.relative)&&x.absolute>=0&&x.relative>=0;}
bool Close(Real residual,Real scale,const K::Bound& b){
  if(!Finite(residual)||!Finite(scale)||scale<0||!BoundValid(b))return false;
  const Real r=std::fabs(residual);
  return r<=b.absolute||(scale>0&&(r-b.absolute)/scale<=b.relative);
}
bool PlannerBoundValid(const R::ClosureBound& b){return Finite(b.absolute)&&Finite(b.relative)&&b.absolute>=0&&b.relative>=0;}
bool Convert(Real value,double* output){
  if(!Finite(value)||std::fabs(value)>std::numeric_limits<double>::max())return false;
  const double converted=static_cast<double>(value);
  if(!std::isfinite(converted)||(value!=0&&converted==0))return false;
  *output=converted;return true;
}
bool ConvertVector(const Vector3& value,R::Vector3* output){
  R::Vector3 candidate;for(unsigned k=0;k<3;++k)if(!Convert(value[k],&candidate[k]))return false;
  *output=candidate;return true;
}
Vector3 Extended(const R::Vector3& value){return {{value[0],value[1],value[2]}};}
Vector3 Rotate(const std::array<Real,9>& r,const Vector3& v){
  Vector3 result{};for(unsigned i=0;i<3;++i)for(unsigned k=0;k<3;++k)result[i]+=r[3*i+k]*v[k];return result;
}
bool InverseRotate(const std::array<Real,9>& r,const Vector3& v,Vector3* output){
  // SPICE supplies rounded coefficients. Proper rotation within an explicit
  // numerical bound does NOT make transpose the exact inverse of those actual
  // coefficients: at a 30 km/s affine offset its error is observable beside
  // a zero LSO parent momentum. Solve the supplied relation, without changing
  // or orthogonalizing its matrix. Partial pivoting is deterministic and the
  // caller verifies the mass-weighted forward roundtrip under EXISTING bounds.
  std::array<std::array<Real,4>,3> a{};
  for(unsigned i=0;i<3;++i){for(unsigned j=0;j<3;++j)a[i][j]=r[3*i+j];a[i][3]=v[i];}
  for(unsigned k=0;k<3;++k){
    unsigned pivot=k;
    for(unsigned i=k+1;i<3;++i)if(std::fabs(a[i][k])>std::fabs(a[pivot][k]))pivot=i;
    if(!Finite(a[pivot][k])||a[pivot][k]==0)return false;
    if(pivot!=k)std::swap(a[pivot],a[k]);
    for(unsigned i=k+1;i<3;++i){
      const Real factor=a[i][k]/a[k][k];
      if(!Finite(factor))return false;
      a[i][k]=0;
      for(unsigned j=k+1;j<4;++j){a[i][j]-=factor*a[k][j];if(!Finite(a[i][j]))return false;}
    }
  }
  Vector3 result{};
  for(int i=2;i>=0;--i){
    Real rhs=a[static_cast<unsigned>(i)][3];
    for(unsigned j=static_cast<unsigned>(i)+1;j<3;++j)rhs-=a[static_cast<unsigned>(i)][j]*result[j];
    result[static_cast<unsigned>(i)]=rhs/a[static_cast<unsigned>(i)][static_cast<unsigned>(i)];
    if(!Finite(result[static_cast<unsigned>(i)]))return false;
  }
  *output=result;return true;
}
bool Kinetic(const Vector3& velocity,Real mass,Real* output){
  const Real scaled=Norm(velocity)*std::sqrt(mass);*output=(scaled/2)*scaled;
  return Finite(*output)&&(*output>0||Norm(velocity)==0);
}
bool MomentumKinetic(const Vector3& momentum,Real mass,Real* output){
  const Real scaled=Norm(momentum)/std::sqrt(mass);*output=(scaled/2)*scaled;
  return Finite(*output)&&(*output>0||Norm(momentum)==0);
}
// An error-free short expansion retains separately supplied formation and
// excitation terms through cancellation. It is bookkeeping, never an energy
// correction. We also compare the binary64 stored-state Q with this value.
struct Sum {
  std::vector<Real> values;
  bool Add(Real value){
    std::size_t kept=0;
    for(Real other:values){
      if(std::fabs(value)<std::fabs(other))std::swap(value,other);
      const Real high=value+other;if(!Finite(high))return false;
      const Real low=other-(high-value);if(low!=0)values[kept++]=low;value=high;
    }
    values.resize(kept);if(value!=0)values.push_back(value);return true;
  }
  Real Value()const{Real result=0;for(auto it=values.rbegin();it!=values.rend();++it)result+=*it;return result;}
};
bool Internal(const K::ResolvedState& state,double* output){
  if(!state.valid())return false;
  Sum total;
  return total.Add(state.definition().formation.joules)&&total.Add(state.definition().excitation.joules)&&Convert(total.Value(),output);
}
std::map<std::string,std::int64_t> Nuclei(const std::vector<R::NucleusCount>& values){
  std::map<std::string,std::int64_t> result;for(const auto& item:values)result.emplace(item.identity,item.count);return result;
}
bool SameState(const K::ResolvedState& state,const R::SpeciesIdentity& material){
  const auto& s=state.definition();
  return state.valid()&&s.identity==material.identity&&s.state==material.state&&s.elementaryCharge==material.elementaryCharge&&
      Nuclei(s.nuclei)==Nuclei(material.nuclei)&&static_cast<double>(state.referenceMassKg())==material.modelMassKg;
}
bool SodiumNuclei(const K::ResolvedState& state){
  const auto& nuclei=state.definition().nuclei;
  return nuclei.size()==1&&nuclei.front().identity=="23Na"&&nuclei.front().count==1;
}
bool ErrorsDeclared(const K::ErrorLimits& limits){
  const Real x[]={limits.maximumRelativeInertialMassError,limits.maximumRelativeAccelerationError,limits.maximumRelativeRelativisticKineticError};
  for(Real v:x)if(!Finite(v)||v<0||v>=1)return false;
  return true;
}
Status ValidateFrame(const EventFrame& frame,const Definition& d){
  if(frame.inertialIdentity.empty()||frame.eventIdentity!=d.plannerAuthorization.eventFrame||frame.transformSourceIdentity.empty()||
      frame.epochConventionIdentity.empty()||!std::isfinite(frame.physicalEventEpochEtS)||!std::isfinite(frame.geometryBasisEpochEtS)||
      !frame.velocityOffsetSupplied||!FiniteVector(frame.velocityOffsetEventMS))return Error(Code::InvalidInput,"actual event/frame/basis epochs and explicitly supplied affine velocity transform provenance are required");
  const auto& r=frame.inertialToEventRotation;
  for(Real x:r)if(!Finite(x))return Error(Code::InvalidInput,"nonfinite inertial-to-event rotation");
  for(unsigned i=0;i<3;++i)for(unsigned j=0;j<3;++j){
    Real inner=0;for(unsigned k=0;k<3;++k)inner+=r[3*k+i]*r[3*k+j];
    if(!Close(inner-(i==j?1:0),1,d.rotation))return Error(Code::InvalidInput,"event velocity transform is not orthogonal within its declared numerical bound");
  }
  const Real det=r[0]*(r[4]*r[8]-r[5]*r[7])-r[1]*(r[3]*r[8]-r[5]*r[6])+r[2]*(r[3]*r[7]-r[4]*r[6]);
  if(!Positive(det)||!Close(det-1,1,d.rotation))return Error(Code::InvalidInput,"event velocity transform is not a proper rotation");
  return {};
}
}

namespace detail {
struct SodiumPhotoionizationData {
  P::Registry source;
  Definition definition;
  P::Channel channel;
  R::Reaction reaction;
  K::FormationBalance formation;
  std::array<K::ModelErrorReport,3> errors;
  double parentInternalJ=0,ionInternalJ=0,electronInternalJ=0;
  Real thresholdJ=0,meanExcessJ=0;
};
}

const char* CodeName(Code code){
  switch(code){case Code::Composed:return "Composed";case Code::Unavailable:return "Unavailable";
    case Code::InvalidInput:return "InvalidInput";case Code::SourceMismatch:return "SourceMismatch";
    case Code::OutsideDomain:return "OutsideDomain";case Code::NumericalFailure:return "NumericalFailure";
    case Code::ResourceLimit:return "ResourceLimit";}
  return "InvalidCode";
}

Status ReferenceCatalogEntry(const K::ResolvedState& state,int compiled,R::SpeciesIdentity* output,double* internal){
  if(!output||!internal||compiled< -1)return Error(Code::InvalidInput,"catalog/internal output and a valid compiled species index are required");
  if(!state.valid())return Error(Code::Unavailable,"reference state has not been resolved from source data");
  try {
    R::SpeciesIdentity value;const auto& s=state.definition();value.identity=s.identity;value.state=s.state;
    value.compiledSpecies=compiled;value.nuclei=s.nuclei;value.elementaryCharge=s.elementaryCharge;double energy=0;
    if(!Convert(state.referenceMassKg(),&value.modelMassKg)||!Positive(value.modelMassKg)||!Internal(state,&energy))
      return Error(Code::NumericalFailure,"state mass/formation/excitation cannot be represented by the actual binary64 planner storage");
    std::string error;if(!R::ValidateCatalog({value},&error))return Error(Code::InvalidInput,error);
    *output=std::move(value);*internal=energy;return {};
  }catch(const std::bad_alloc&){return Error(Code::ResourceLimit,"catalog binding allocation failed; reduce concurrently retained state profiles or provide more memory");}
  catch(const std::exception& e){return Error(Code::InvalidInput,std::string("catalog binding: ")+e.what());}
}

Status SodiumPhotoionization::Create(const P::Registry& source,const Definition& d,SodiumPhotoionization* output){
  if(!output)return Error(Code::InvalidInput,"provider output storage is required");
  if(!source.loaded())return Error(Code::Unavailable,"complete frozen D04 registry is unavailable");
  if(!d.sourceStateCrosswalkAvailable)return Error(Code::Unavailable,"reviewed Na source/state/threshold crosswalk is unavailable; no transport is admitted");
  if(!d.parent.state.valid()||!d.ion.state.valid()||!d.electron.state.valid())
    return Error(Code::Unavailable,"required Na/Na+/electron physical states are unresolved");
  if(d.identity.empty()||d.sourceProfileIdentity.empty()||d.stateCrosswalkIdentity.empty()||d.energyReferenceIdentity.empty()||
      !ReactionData::DigestSyntax(d.sourceProfileSha256)||!ReactionData::DigestSyntax(d.stateCrosswalkSha256)||
      d.electronVoltJ!=ExactElectronVoltJ||d.electronVoltSourceIdentity.empty())
    return Error(Code::InvalidInput,"explicit pinned profile/crosswalk, common energy reference and exact source-identified SI eV conversion are required");
  const int endpoint=static_cast<int>(d.solar);
  if(endpoint!=0&&endpoint!=1)return Error(Code::InvalidInput,"quiet and active are distinct source endpoints; no interpolation is selected");
  if(d.photonEnergyPolicy!=K::PhotonEnergyPolicy::SourceEqualsState&&
      d.photonEnergyPolicy!=K::PhotonEnergyPolicy::SourceEffectivePhoton)
    return Error(Code::InvalidInput,"unknown photon source/physical-state energy convention");
  if(d.photonEnergyPolicy==K::PhotonEnergyPolicy::SourceEffectivePhoton){
    if(d.physicalThresholdSourceIdentity.empty())
      return Error(Code::Unavailable,"source-effective photon needs the authenticated physical threshold source identity");
    if(d.physicalThresholdSourceIdentity!=d.ion.state.definition().formation.sourceIdentity)
      return Error(Code::SourceMismatch,"source-effective threshold evidence differs from the resolved ion formation-energy source");
  }
  if(!d.plannerAuthorization.qualified||d.plannerAuthorization.identity.empty()||d.plannerAuthorization.evidenceIdentity.empty())
    return Error(Code::Unavailable,"caller transport-policy admission is missing; the composer does not qualify physical kinetics");
  const auto& a=d.plannerAuthorization;
  if(a.eventFrame.empty()||!a.products.empty()||a.bath.supplied||!Positive(a.maximumClassicalSpeedMS)||
      a.maximumClassicalSpeedMS>=K::SpeedOfLightMS||a.maximumClassicalSpeedMS!=static_cast<double>(d.numerical.maximumClassicalSpeedMS)||
      !PlannerBoundValid(a.bounds.modelMassKg)||!PlannerBoundValid(a.bounds.momentumKgMS)||!PlannerBoundValid(a.bounds.energyJ)||
      !BoundValid(d.rotation)||d.rotation.absolute+d.rotation.relative>=1)
    return Error(Code::InvalidInput,"authorization needs empty product/bath slots and matching explicit classical, closure and rotation bounds");
  try {
    std::string error;
    if(!R::ValidateCatalog(d.catalog,&error))return Error(Code::InvalidInput,error);
    const std::size_t ids[]={d.parent.catalogIndex,d.ion.catalogIndex,d.electron.catalogIndex};
    if(ids[0]==ids[1]||ids[0]==ids[2]||ids[1]==ids[2])return Error(Code::InvalidInput,"parent, ion and electron require distinct catalog states");
    const StateBinding* states[]={&d.parent,&d.ion,&d.electron};
    for(unsigned i=0;i<3;++i){
      if(ids[i]>=d.catalog.size()||!SameState(states[i]->state,d.catalog[ids[i]]))
        return Error(Code::SourceMismatch,"catalog identity, state, nuclei, charge or actual binary64 mass differs from its resolved source state");
      if(!Positive(states[i]->physicalRestMass.kg)||states[i]->physicalRestMass.sourceIdentity.empty()||!ErrorsDeclared(states[i]->modelErrorLimits))
        return Error(Code::Unavailable,"independent physical rest masses and reviewed reference-inertia/relativistic error limits are required");
    }
    if(d.catalog[ids[0]].compiledSpecies<0||d.catalog[ids[1]].compiledSpecies<0||d.catalog[ids[2]].compiledSpecies!=-1||
        !SodiumNuclei(d.parent.state)||!SodiumNuclei(d.ion.state)||!d.electron.state.definition().nuclei.empty()||
        d.parent.state.definition().elementaryCharge!=0||d.ion.state.definition().elementaryCharge!=1||d.electron.state.definition().elementaryCharge!=-1)
      return Error(Code::SourceMismatch,"reference scope is tracked23Na -> tracked23Na+ + untracked free electron, with exact integer charge");
    if(d.parent.state.referenceInventoryIdentity()!=d.ion.state.referenceInventoryIdentity()||
        d.parent.state.referenceInventoryIdentity()!=d.electron.state.referenceInventoryIdentity())
      return Error(Code::SourceMismatch,"event states use different reference-mass inventories");
    std::vector<P::Channel> family;
    if(!source.ChannelsForParent("Na",&family,&error))return Error(Code::SourceMismatch,error);
    if(family.size()!=1||!family.front().reactive||family.front().id!=d.channelId)
      return Error(Code::SourceMismatch,"Na provider requires its COMPLETE single reactive D04 family; no branch dropping or renormalization");
    const auto& channel=family.front();
    if(channel.physicalParent!="Na"||!channel.hasThreshold||channel.products.size()!=2)
      return Error(Code::Unavailable,"source Na channel lacks a complete threshold/product contract");
    bool ion=false,electron=false;
    for(const auto& product:channel.products){
      if(product.multiplicity!=1)return Error(Code::SourceMismatch,"unsupported Na source product multiplicity");
      if(product.species=="Na+"&&product.charge==1&&!ion)ion=true;
      else if(product.species=="e"&&product.charge==-1&&!electron)electron=true;
      else return Error(Code::SourceMismatch,"source channel is not complete atomic Na photoionization");
    }
    if(!ion||!electron)return Error(Code::SourceMismatch,"source channel omits ion or reservoir electron");
    auto data=std::make_shared<detail::SodiumPhotoionizationData>();data->source=source;data->definition=d;data->channel=channel;
    data->thresholdJ=static_cast<Real>(channel.thresholdEnergyEv)*d.electronVoltJ;
    data->meanExcessJ=static_cast<Real>(channel.meanExcessEnergyEv[endpoint])*d.electronVoltJ;
    if(!Positive(data->thresholdJ)||!Finite(data->meanExcessJ)||data->meanExcessJ<0)
      return Error(Code::SourceMismatch,"source threshold/selected endpoint mean excess is invalid");
    if(!K::FormationQ({d.parent.state},{d.ion.state,d.electron.state},&data->formation,&error))return Error(Code::SourceMismatch,error);
    // The strict historical policy retains EXACT source/state consistency.
    // The alternative is an explicit source model authenticated by the owner:
    // the archived photon remains I_source+mean, while the FULL physical Q
    // below supplies -I_physical. No caller numerical bound admits a mismatch
    // under the strict policy, and no arbitrary correction is accepted under
    // the source-effective policy. Both energy conventions use the SAME rate.
    if(d.photonEnergyPolicy==K::PhotonEnergyPolicy::SourceEqualsState){
    // Physical source consistency is NOT a user numerical-closure tolerance.
    // Compare every unrounded native formation/excitation component with the
    // exact SI conversion of the stored source eV value. Two adjacent native
    // spacings cover only encoding the SI coefficient and its product into
    // this module's long-double state representation. They are not a source
    // uncertainty. The envelope depends solely on the small source threshold,
    // never on large cancelling energy zeros or caller thresholdJ settings.
    // An unrepresentable crosswalk is refused rather than widened.
    std::vector<Real> sourceGap;
    for(const auto* state:{&d.ion.state,&d.electron.state}) {
      sourceGap.push_back(state->definition().formation.joules);
      sourceGap.push_back(state->definition().excitation.joules);
    }
    sourceGap.push_back(-d.parent.state.definition().formation.joules);
    sourceGap.push_back(-d.parent.state.definition().excitation.joules);
    const Real spacing=std::nextafter(data->thresholdJ,std::numeric_limits<Real>::infinity())-data->thresholdJ;
    bool sourceMatches=false;
    if(!Positive(spacing)||!ExactArithmetic::SignedSumMatchesDecimalScale(sourceGap,
        channel.thresholdEnergyEv,UINT64_C(1602176634),28,2*spacing,&sourceMatches,&error))
      return Error(Code::NumericalFailure,"cannot establish exact source/state representation consistency: "+error);
    if(!sourceMatches)
      return Error(Code::SourceMismatch,"exact source threshold/state gap disagree beyond native representation; caller numerical tolerances cannot admit a different physical threshold");
    if(!Close(data->formation.qJ+data->thresholdJ,std::max(std::fabs(data->formation.qJ),data->thresholdJ),d.numerical.thresholdJ))
      return Error(Code::SourceMismatch,"native D04 threshold and supplied physical state gap disagree; an unresolved physical crosswalk cannot be hidden by a numerical tolerance");
    }else if(!Positive(-data->formation.qJ)){
      return Error(Code::SourceMismatch,"physical Na ionization crosswalk requires a finite positive formation-state gap");
    }
    if(!Internal(d.parent.state,&data->parentInternalJ)||!Internal(d.ion.state,&data->ionInternalJ)||!Internal(d.electron.state,&data->electronInternalJ))
      return Error(Code::NumericalFailure,"source state energies cannot be stored in the planner's binary64 internal-energy field");
    Sum storedQ;storedQ.Add(data->parentInternalJ);storedQ.Add(-static_cast<Real>(data->ionInternalJ));storedQ.Add(-static_cast<Real>(data->electronInternalJ));
    if(!Close(storedQ.Value()-data->formation.qJ,std::max(std::fabs(data->formation.qJ),data->thresholdJ),d.numerical.thresholdJ))
      return Error(Code::NumericalFailure,"binary64 state storage erases a source-resolved reaction Q; choose a supported energy zero/storage rather than weakening the bound");
    // Validate the SAME production two-body controls and additive mass/Q
    // contract at rest before any event/RNG/particle operation is possible.
    K::PhotonTwoBodyInput probe;probe.parent=d.parent.state;probe.products={{d.ion.state,d.electron.state}};
    probe.photonDirection={{1,0,0}};probe.channelThreshold={data->thresholdJ,d.stateCrosswalkIdentity};
    probe.energyPolicy=d.photonEnergyPolicy;
    probe.sourceThresholdCrosswalkIdentity=d.stateCrosswalkIdentity+"@sha256:"+d.stateCrosswalkSha256;
    probe.meanExcess={data->meanExcessJ,channel.rateSourceFile[endpoint]+":"+channel.rateSourceRecord[endpoint]};
    probe.polarOpenUniform=.5L;probe.azimuthOpenUniform=.5L;
    probe.context={"parent-rest-inertial-admission",d.identity};probe.controls=d.numerical;K::FinalState temporary;
    if(!K::PhotonTwoBody(probe,&temporary,&error))return Error(Code::OutsideDomain,error);
    for(unsigned i=0;i<3;++i)if(!K::EvaluateModelErrors(states[i]->state,states[i]->physicalRestMass,d.numerical.maximumClassicalSpeedMS,
        states[i]->modelErrorLimits,&data->errors[i],&error))return Error(Code::OutsideDomain,error);
    data->reaction.identity=channel.id+";"+Model+";"+d.stateCrosswalkIdentity;
    if(d.photonEnergyPolicy==K::PhotonEnergyPolicy::SourceEffectivePhoton)
      data->reaction.identity+=";source-effective-photon";
    data->reaction.reactants={{d.parent.catalogIndex,1,R::Population::Tracked}};
    data->reaction.products={{d.ion.catalogIndex,1,R::Population::Tracked},{d.electron.catalogIndex,1,R::Population::Bath}};
    if(!R::ValidateReaction(d.catalog,data->reaction,&error))return Error(Code::SourceMismatch,error);
    SodiumPhotoionization candidate;candidate.data_=std::move(data);*output=std::move(candidate);return {};
  }catch(const std::bad_alloc&){return Error(Code::ResourceLimit,"source/state binding allocation failed; reduce simultaneously retained profiles or provide more memory");}
  catch(const std::exception& e){return Error(Code::InvalidInput,std::string("Na source binding: ")+e.what());}
}

const std::vector<R::SpeciesIdentity>* SodiumPhotoionization::catalog()const{return data_?&data_->definition.catalog:nullptr;}
const R::Reaction* SodiumPhotoionization::reaction()const{return data_?&data_->reaction:nullptr;}
const P::Channel* SodiumPhotoionization::sourceChannel()const{return data_?&data_->channel:nullptr;}
Status SodiumPhotoionization::RateS(double distance,int visible,double* output)const{
  if(!output)return Error(Code::InvalidInput,"rate output storage is required");
  if(!data_)return Error(Code::Unavailable,"Na source/state provider is not prepared");
  if(visible!=0&&visible!=1)return Error(Code::InvalidInput,"this point-Sun model accepts exactly0/1 gas visibility, not a facet cosine or finite-disc fraction");
  std::string error;double value=0;
  if(!data_->source.TotalDestructionRateS("Na",data_->definition.solar,distance,visible,&value,&error))return Error(Code::OutsideDomain,error);
  *output=value;return {};
}

Status SodiumPhotoionization::Compose(const EventInput& input,Composition* output)const{
  if(!output)return Error(Code::InvalidInput,"event composition output storage is required");
  if(!data_)return Error(Code::Unavailable,"Na source/state provider is not prepared");
  const auto& d=data_->definition;
  if(!input.interval.valid()||input.parentAtEvent.species!=d.parent.catalogIndex||input.parentAtEvent.internalEnergyJ!=data_->parentInternalJ||
      !Positive(input.parentAtEvent.referenceWeight)||!Positive(input.parentAtEvent.individualCorrection)||
      !FiniteVector(Extended(input.parentAtEvent.positionM))||!FiniteVector(Extended(input.parentAtEvent.velocityMS))||
      !Positive(input.ionReferenceWeight)||!input.freshIonClock.valid()||input.freshIonClock.consumed()||!input.freshIonClock.fresh())
    return Error(Code::SourceMismatch,"actual parent state/energy/weight/interval and one fresh ion clock must match the admitted source catalog");
  if(input.photonGeometryIdentity.empty()||input.pointSunVisible!=1)
    return Error(Code::SourceMismatch,"a physical photon event requires actual unobscured point-Sun geometry; no event is created in shadow");
  if(!(input.polarOpenUniform>0&&input.polarOpenUniform<1)||!(input.azimuthOpenUniform>0&&input.azimuthOpenUniform<1))
    return Error(Code::InvalidInput,"isotropic recoil requires two supplied strictly open uniforms; no internal resampling exists");
  const auto frameStatus=ValidateFrame(input.frame,d);if(!frameStatus.ok())return frameStatus;
  // Absolute ET is a projection of the immutable relative history. Reuse the
  // SAME established two-addition binary64 envelope as normal Moon tracing;
  // this neither advances exposure nor introduces a chemical-time epsilon.
  const double projectedEpoch=static_cast<double>(static_cast<Real>(input.interval.startS())+input.interval.elapsedS());
  if(!TransportTime::EpochProjectionAgrees(input.frame.physicalEventEpochEtS,projectedEpoch))
    return Error(Code::SourceMismatch,"actual event epoch disagrees with the immutable interval's supported binary64 projection");
  try {
    Composition value;auto& receipt=value.receipt;
    receipt.modelIdentity=Model;receipt.sourceProfileIdentity=d.sourceProfileIdentity;receipt.sourceProfileSha256=d.sourceProfileSha256;
    receipt.stateCrosswalkIdentity=d.stateCrosswalkIdentity;receipt.stateCrosswalkSha256=d.stateCrosswalkSha256;
    receipt.channelId=data_->channel.id;receipt.energyReferenceIdentity=d.energyReferenceIdentity;
    receipt.photonGeometryIdentity=input.photonGeometryIdentity;receipt.solar=d.solar;
    receipt.sourceFiles=data_->source.sources();receipt.sourceChannel=data_->channel;receipt.frame=input.frame;receipt.interval=input.interval;
    receipt.domainErrorReports=data_->errors;receipt.nativeThresholdJ=data_->thresholdJ;receipt.matchedSourceMeanExcessJ=data_->meanExcessJ;
    receipt.photonEnergyPolicy=d.photonEnergyPolicy==K::PhotonEnergyPolicy::SourceEqualsState?
      "source-equals-state":"source-effective-photon";
    receipt.physicalThresholdJ=-data_->formation.qJ;
    Sum thresholdShift;
    if(!thresholdShift.Add(receipt.physicalThresholdJ)||!thresholdShift.Add(-data_->thresholdJ))
      return Error(Code::NumericalFailure,"physical/source threshold difference is not representable");
    receipt.physicalMinusSourceThresholdJ=thresholdShift.Value();
    receipt.intervalProjectedEventEpochEtS=projectedEpoch;
    receipt.stateQJ=data_->formation.qJ;
    auto status=RateS(input.distanceAu,input.pointSunVisible,&receipt.eventRateS);if(!status.ok())return status;
    std::string error;if(!data_->source.ChannelRate1AuS(d.channelId,d.solar,&receipt.sourceRate1AuS,&error))return Error(Code::SourceMismatch,error);
    Vector3 eventParent=Extended(input.parentAtEvent.velocityMS),shifted;
    for(unsigned k=0;k<3;++k)shifted[k]=eventParent[k]-input.frame.velocityOffsetEventMS[k];
    Vector3 parent;
    if(!InverseRotate(input.frame.inertialToEventRotation,shifted,&parent))
      return Error(Code::NumericalFailure,"actual rounded affine parent transform cannot be inverted");
    Vector3 roundtrip=Rotate(input.frame.inertialToEventRotation,parent);
    for(unsigned k=0;k<3;++k){
      roundtrip[k]+=input.frame.velocityOffsetEventMS[k];
      const Real mass=d.parent.state.referenceMassKg();
      const Real residual=mass*(roundtrip[k]-eventParent[k]);
      const Real scale=mass*std::max(std::fabs(roundtrip[k]),std::fabs(eventParent[k]));
      if(!Finite(roundtrip[k])||!Close(residual,scale,d.numerical.momentumKgMS))
        return Error(Code::NumericalFailure,"actual rounded affine parent roundtrip fails its existing mass-weighted momentum bound");
    }
    receipt.parentInertialVelocityMS=parent;
    if(!FiniteVector(parent)||Norm(parent)>d.numerical.maximumClassicalSpeedMS||Norm(eventParent)>d.numerical.maximumClassicalSpeedMS)
      return Error(Code::OutsideDomain,"parent velocity exceeds the explicitly reviewed classical domain in an actual event/inertial frame");
    auto direction=input.incomingPhotonDirectionInertial;const Real n=Norm(direction);
    if(!FiniteVector(direction)||!Positive(n)||std::fabs(n-1)>d.numerical.unitVectorTolerance)
      return Error(Code::InvalidInput,"incoming photon propagation direction must be an explicit nonzero unit vector");
    for(Real& component:direction)component/=n;
    const Real beta=Norm(parent)/K::SpeedOfLightMS;
    const Real gamma=1/std::sqrt((1-beta)*(1+beta));
    // The direction is photon PROPAGATION (Sun -> particle), opposite a
    // particle -> Sun visibility ray. For this approved representative model
    // E_rest=gamma*E_inertial*(1-v dot n/c). This transforms ONE selected
    // photon, not the solar spectrum or its independently frozen rate law.
    const Real doppler=gamma*(1-Dot(parent,direction)/K::SpeedOfLightMS);
    const Real restPhoton=data_->thresholdJ+data_->meanExcessJ;
    const Real photon=restPhoton/doppler;
    if(!Positive(doppler)||!Positive(photon)||!Close(doppler*photon-restPhoton,restPhoton,d.numerical.energyJ))
      return Error(Code::NumericalFailure,"representative photon Doppler transform is unrepresentable or fails its numerical energy bound");
    receipt.incomingPhoton.restEnergyJ=restPhoton;receipt.incomingPhoton.inertialEnergyJ=photon;
    receipt.incomingPhoton.dopplerFactor=doppler;receipt.incomingPhoton.inertialDirection=direction;
    Vector3 total;
    for(unsigned k=0;k<3;++k){
      receipt.incomingPhoton.inertialMomentumKgMS[k]=(photon/K::SpeedOfLightMS)*direction[k];
      total[k]=d.parent.state.referenceMassKg()*parent[k]+receipt.incomingPhoton.inertialMomentumKgMS[k];
    }
    // Retain the inverse-boosted momentum as well as energy. Substituting
    // gamma^2/(gamma+1)/c^2 for (gamma-1)/v^2 avoids a 0/0 at rest and
    // cancellation in gamma-1 for a slow parent. This is an exact photon
    // Lorentz transform; it supplies no outgoing angular/spectral PDF.
    const Real photonBoost=(gamma*gamma/(gamma+1)*Dot(parent,receipt.incomingPhoton.inertialMomentumKgMS)-gamma*photon)/
        (K::SpeedOfLightMS*K::SpeedOfLightMS);
    for(unsigned k=0;k<3;++k)receipt.incomingPhoton.restMomentumKgMS[k]=receipt.incomingPhoton.inertialMomentumKgMS[k]+photonBoost*parent[k];
    if(!Close(Norm(receipt.incomingPhoton.restMomentumKgMS)*K::SpeedOfLightMS-restPhoton,restPhoton,d.numerical.energyJ))
      return Error(Code::NumericalFailure,"incoming photon four-vector fails the explicit parent-rest massless-energy bound");
    Real beforeK=0,comK=0;const Real productMass=d.ion.state.referenceMassKg()+d.electron.state.referenceMassKg();
    if(!Kinetic(parent,d.parent.state.referenceMassKg(),&beforeK)||!MomentumKinetic(total,productMass,&comK))
      return Error(Code::NumericalFailure,"incoming material/COM kinetic energy is unrepresentable");
    // Keep the source mean separate through cancellation with formation Q.
    // A negative result is forbidden, never clipped or empirically corrected.
    Sum available;
    if(!available.Add(beforeK)||!available.Add(data_->thresholdJ/doppler)||!available.Add(data_->meanExcessJ/doppler))
      return Error(Code::NumericalFailure,"incoming representative energy sum overflowed");
    for(const auto& term:data_->formation.terms)for(Real component:{term.formation.joules,term.excitation.joules})
      if(!available.Add(term.incoming?component:-component))return Error(Code::NumericalFailure,"formation/Q sum overflowed");
    const Real totalAvailableKinetic=available.Value();
    if(!available.Add(-comK)||!Finite(available.Value())||available.Value()<0)
      return Error(Code::OutsideDomain,"source photon and declared states leave negative COM kinetic energy; no clipping is allowed");
    K::TwoBodyInput two;two.products={{d.ion.state,d.electron.state}};two.totalMomentumKgMS=total;
    if(!K::IsotropicDirectionFromOpenUniforms(input.polarOpenUniform,input.azimuthOpenUniform,
        &two.relativeDirection,&error))return Error(Code::InvalidInput,error);
    two.availableCOMKineticJ=available.Value();two.controls=d.numerical;
    two.context={input.frame.inertialIdentity,d.identity+";"+d.stateCrosswalkIdentity};
    if(!K::GeneralTwoBodyFinalState(two,&receipt.inertialProducts,&error))return Error(Code::NumericalFailure,error);
    receipt.inertialProducts.modelIdentity=Model;receipt.inertialProducts.photonIncluded=true;
    receipt.inertialProducts.photonEnergyJ=photon;receipt.inertialProducts.photonMomentumKgMS=receipt.incomingPhoton.inertialMomentumKgMS;
    receipt.inertialProducts.stateQJ=data_->formation.qJ;receipt.inertialProducts.formationBalance=data_->formation;
    receipt.inertialProducts.channelThresholdJ=data_->thresholdJ;receipt.inertialProducts.sourceMeanExcessJ=data_->meanExcessJ;
    receipt.inertialProducts.stateAvailableKineticJ=totalAvailableKinetic;
    const auto& ion=receipt.inertialProducts.products[0];const auto& electron=receipt.inertialProducts.products[1];
    Vector3 ionEvent=Rotate(input.frame.inertialToEventRotation,ion.velocityMS);
    Vector3 electronEvent=Rotate(input.frame.inertialToEventRotation,electron.velocityMS);
    for(unsigned k=0;k<3;++k){ionEvent[k]+=input.frame.velocityOffsetEventMS[k];electronEvent[k]+=input.frame.velocityOffsetEventMS[k];}
    if(!FiniteVector(ionEvent)||!FiniteVector(electronEvent)||Norm(ionEvent)>d.numerical.maximumClassicalSpeedMS||Norm(electronEvent)>d.numerical.maximumClassicalSpeedMS)
      return Error(Code::OutsideDomain,"outgoing event-coordinate speed exceeds the explicitly reviewed classical domain");
    auto& e=receipt.outgoingElectron;e.state=d.electron.state;e.inertialVelocityMS=electron.velocityMS;e.eventVelocityMS=electronEvent;
    e.inertialMomentumKgMS=electron.momentumKgMS;
    Sum electronState;
    if(!electronState.Add(d.electron.state.definition().formation.joules)||
        !electronState.Add(d.electron.state.definition().excitation.joules))
      return Error(Code::NumericalFailure,"source electron state-energy sum is unrepresentable");
    e.stateEnergyJ=electronState.Value();
    if(!Kinetic(e.inertialVelocityMS,d.electron.state.referenceMassKg(),&e.inertialKineticJ)||
        !Kinetic(e.eventVelocityMS,d.electron.state.referenceMassKg(),&e.eventKineticJ))
      return Error(Code::NumericalFailure,"outgoing reservoir electron energy is unrepresentable");
    for(unsigned k=0;k<3;++k)e.eventMomentumKgMS[k]=d.electron.state.referenceMassKg()*electronEvent[k];
    receipt.bathModelMassChangeKg=d.electron.state.referenceMassKg();
    receipt.bathInertialEnergyChangeJ=e.inertialKineticJ+e.stateEnergyJ-photon;
    for(unsigned k=0;k<3;++k)receipt.bathInertialMomentumChangeKgMS[k]=e.inertialMomentumKgMS[k]-receipt.incomingPhoton.inertialMomentumKgMS[k];
    const Vector3 rotatedBath=Rotate(input.frame.inertialToEventRotation,receipt.bathInertialMomentumChangeKgMS);
    Real offsetK=0;
    if(!Kinetic(input.frame.velocityOffsetEventMS,receipt.bathModelMassChangeKg,&offsetK))return Error(Code::NumericalFailure,"affine reservoir kinetic term is unrepresentable");
    // Coordinate transformation of the independently assembled signed bath:
    // dp'=R dp+dm*w; dE'=dE+w dot Rdp+dm*w^2/2. This is the classical
    // EVENT-coordinate ledger, not a rotating-frame photon frequency law.
    for(unsigned k=0;k<3;++k)receipt.bathEventMomentumChangeKgMS[k]=rotatedBath[k]+receipt.bathModelMassChangeKg*input.frame.velocityOffsetEventMS[k];
    receipt.bathEventEnergyChangeJ=receipt.bathInertialEnergyChangeJ+Dot(input.frame.velocityOffsetEventMS,rotatedBath)+offsetK;
    value.policy=d.plannerAuthorization;value.policy.bath.supplied=true;
    if(!Convert(receipt.bathModelMassChangeKg,&value.policy.bath.modelMassChangeKg)||
        !Convert(receipt.bathEventEnergyChangeJ,&value.policy.bath.energyChangeJ)||
        !ConvertVector(receipt.bathEventMomentumChangeKgMS,&value.policy.bath.momentumChangeKgMS))
      return Error(Code::NumericalFailure,"signed physical reservoir moments cannot be represented by the planner");
    R::ProductKinematics product;product.term=0;product.occurrence=0;product.referenceWeight=input.ionReferenceWeight;
    product.internalEnergyJ=data_->ionInternalJ;product.nextClock=input.freshIonClock;
    if(!ConvertVector(ionEvent,&product.velocityMS))return Error(Code::NumericalFailure,"ion velocity cannot be represented by actual particle storage");
    value.policy.products.push_back(product);
    if(!R::MakeEventPlan(d.catalog,data_->reaction,input.parentAtEvent,input.interval,R::PlanScope::Transport,&value.policy,&value.plan,&error))
      return Error(Code::NumericalFailure,"complete source/photon/electron/event-coordinate packet failed planner closure: "+error);
    // A photon is not inserted as an invalid zero-mass material catalog
    // species. Its radiation count is nevertheless explicit: one incident
    // photon per event, and one outgoing reservoir electron per event.
    receipt.representedElectronNumberChange=value.plan.ledger().parentNumber;
    receipt.representedPhotonNumberChange=-value.plan.ledger().parentNumber;
    *output=std::move(value);return {};
  }catch(const std::bad_alloc&){return Error(Code::ResourceLimit,"event composition allocation failed; retain the selected event and reduce simultaneous histories or provide more memory; never discard an electron/ion or redraw products");}
  catch(const std::exception& e){return Error(Code::NumericalFailure,std::string("Na event composition: ")+e.what());}
}

Status SodiumPhotoionization::IonCollectorState(IonCollector::StateDefinition* output)const{
  if(!output)return Error(Code::InvalidInput,"collector-state output storage is required");
  if(!data_)return Error(Code::Unavailable,"Na source/state provider is not prepared");
  try {
    IonCollector::StateDefinition state;const auto& d=data_->definition;
    state.material=d.catalog[d.ion.catalogIndex];
    if(!Convert(d.ion.state.definition().formation.joules,&state.formationEnergyJ)||
        !Convert(d.ion.state.definition().excitation.joules,&state.excitationEnergyJ))
      return Error(Code::NumericalFailure,"ion collector cannot store source formation/excitation components");
    state.energyReferenceIdentity=d.energyReferenceIdentity;
    state.sourceEvidenceIdentity=d.sourceProfileIdentity+";"+d.sourceProfileSha256+";"+d.stateCrosswalkIdentity+";"+d.stateCrosswalkSha256;
    std::string error;if(!IonCollector::ValidateStates({state},&error))return Error(Code::SourceMismatch,error);
    *output=std::move(state);return {};
  }catch(const std::bad_alloc&){return Error(Code::ResourceLimit,"collector-state allocation failed; reduce retained profiles or provide more memory");}
  catch(const std::exception& e){return Error(Code::InvalidInput,std::string("collector binding: ")+e.what());}
}
} } }
