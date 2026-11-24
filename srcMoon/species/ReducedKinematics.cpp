#include "ReducedKinematics.h"
#include <algorithm>
#include <cmath>
#include <exception>
#include <iomanip>
#include <limits>
#include <locale>
#include <map>
#include <set>
#include <sstream>
#include <utility>

namespace Moon { namespace Species { namespace ReducedKinematics {
namespace {
bool Fail(const std::string& text,std::string* error){if(error)*error=text;return false;}
bool Finite(Real x){return std::isfinite(x);}
bool Positive(Real x){return Finite(x)&&x>0;}
bool VectorFinite(const Vector3& x){return Finite(x[0])&&Finite(x[1])&&Finite(x[2]);}
Real Norm(const Vector3& x){return std::hypot(std::hypot(x[0],x[1]),x[2]);}
// A short floating expansion retains separately supplied formation/excitation
// terms through cancellation. Collapsing (1e20+1) before subtracting1e20 would
// erase a real1J Q. Error-free FastTwoSum keeps nonoverlapping residuals;
// there is no physical threshold or ad-hoc energy correction here.
struct EnergySum {
  std::vector<Real> partials;
  bool Add(Real x) {
    std::size_t kept=0;
    for(Real y:partials) {
      if(std::fabs(x)<std::fabs(y))std::swap(x,y);
      const Real high=x+y;
      if(!Finite(high))return false;
      const Real low=y-(high-x);
      if(low!=0)partials[kept++]=low;
      x=high;
    }
    partials.resize(kept);if(x!=0)partials.push_back(x);return true;
  }
  Real Value()const {Real total=0;for(auto it=partials.rbegin();it!=partials.rend();++it)total+=*it;return total;}
};
bool BoundValid(const Bound& b){return Finite(b.absolute)&&Finite(b.relative)&&b.absolute>=0&&b.relative>=0;}
bool Close(Real residual,Real scale,const Bound& b) {
  const Real value=std::fabs(residual);
  if(!Finite(value)||!Finite(scale)||scale<0||!BoundValid(b))return false;
  if(value<=b.absolute)return true;
  return scale>0&&(value-b.absolute)/scale<=b.relative;
}
bool Controls(const Context& context,const NumericControls& c,std::string* error) {
  if(context.frameIdentity.empty()||context.sourceIdentity.empty())return Fail("explicit inertial frame and source/model identity are required",error);
  if(!Positive(c.maximumClassicalSpeedMS)||c.maximumClassicalSpeedMS>=SpeedOfLightMS||
      !Finite(c.unitVectorTolerance)||c.unitVectorTolerance<0||c.unitVectorTolerance>=1||!BoundValid(c.referenceMassKg)||
      !BoundValid(c.momentumKgMS)||!BoundValid(c.energyJ)||!BoundValid(c.thresholdJ))
    return Fail("all numerical bounds and a finite subluminal speed guard must be explicitly supplied",error);
  return true;
}
bool Unit(const Vector3& x,Real tolerance,std::string* error) {
  if(!VectorFinite(x)||!Positive(Norm(x))||std::fabs(Norm(x)-1)>tolerance)return Fail("direction must be an explicit nonzero unit vector within the supplied numerical bound",error);
  return true;
}
bool StateValid(const ResolvedState& s,std::string* error) {
  if(!s.valid()||!Positive(s.referenceMassKg())||s.referenceInventoryIdentity().empty())
    return Fail("product/parent state was not resolved from an explicit source inventory",error);
  return true;
}
Real StateEnergy(const ResolvedState& s){return s.definition().formation.joules+s.definition().excitation.joules;}
void Sources(const ResolvedState& s,std::vector<std::string>* out) {
  out->insert(out->end(),s.sourceIdentities().begin(),s.sourceIdentities().end());
}
bool EnergyValid(const EnergyValue& e,bool nonnegative,std::string* error) {
  if(!Finite(e.joules)||(nonnegative&&e.joules<0)||e.sourceIdentity.empty())
    return Fail("every state/channel energy requires a finite value, explicit source and valid sign",error);
  return true;
}
// Momentum-based kinetic energy avoids p^2 overflow before division by mass.
bool MomentumKinetic(const Vector3& p,Real mass,Real* kinetic) {
  const Real root=Norm(p)/std::sqrt(mass);*kinetic=(root/2)*root;
  return Finite(*kinetic)&&(*kinetic>0||Norm(p)==0);
}
bool VelocityKinetic(const Vector3& v,Real mass,Real* kinetic) {
  const Real root=Norm(v)*std::sqrt(mass);*kinetic=(root/2)*root;
  return Finite(*kinetic)&&(*kinetic>0||Norm(v)==0);
}
bool Finalize(const std::vector<ResolvedState>& states,const std::vector<Vector3>& momenta,
    const Vector3& targetMomentum,Real available,const Context& context,
    const NumericControls& controls,const char* model,FinalState* output,std::string* error) {
  FinalState value;value.modelIdentity=model;value.frameIdentity=context.frameIdentity;
  value.sourceIdentities.push_back(context.sourceIdentity);value.availableCOMKineticJ=available;
  Real totalMass=0,momentumScale=Norm(targetMomentum);Vector3 actual{{0,0,0}};
  for(const auto& s:states)totalMass+=s.referenceMassKg();
  if(!Positive(totalMass))return Fail("combined reference mass is unrepresentable",error);
  Vector3 com;for(int d=0;d<3;++d)com[d]=targetMomentum[d]/totalMass;
  for(std::size_t i=0;i<states.size();++i) {
    Product product;product.state=states[i];product.occurrence=i;product.momentumKgMS=momenta[i];
    Vector3 reconstructed,relative;
    for(int d=0;d<3;++d) {
      product.velocityMS[d]=momenta[i][d]/states[i].referenceMassKg();
      reconstructed[d]=states[i].referenceMassKg()*product.velocityMS[d];
      relative[d]=product.velocityMS[d]-com[d];actual[d]+=reconstructed[d];
    }
    if(!VectorFinite(product.velocityMS)||!VectorFinite(reconstructed))return Fail("product momentum/velocity is unrepresentable",error);
    const Real speed=Norm(product.velocityMS);
    if(!Finite(speed)||speed>controls.maximumClassicalSpeedMS)return Fail("product speed exceeds the explicitly supplied nonrelativistic validity guard [m/s]",error);
    value.maximumProductSpeedMS=std::max(value.maximumProductSpeedMS,speed);
    momentumScale=std::max(momentumScale,Norm(reconstructed));
    Real kinetic=0,relativeK=0;
    if(!VelocityKinetic(product.velocityMS,states[i].referenceMassKg(),&kinetic)||
        !VelocityKinetic(relative,states[i].referenceMassKg(),&relativeK))
      return Fail("product kinetic energy is outside the numerical range",error);
    value.totalKineticJ+=kinetic;value.actualCOMKineticJ+=relativeK;
    value.products.push_back(std::move(product));Sources(states[i],&value.sourceIdentities);
  }
  value.totalMomentumKgMS=actual;
  for(int d=0;d<3;++d)value.momentumResidualKgMS[d]=actual[d]-targetMomentum[d];
  value.energyResidualJ=value.actualCOMKineticJ-available;
  Real comEnergy=0;
  if(!MomentumKinetic(targetMomentum,totalMass,&comEnergy)||!Finite(value.totalKineticJ)||
      !Close(Norm(value.momentumResidualKgMS),momentumScale,controls.momentumKgMS)||
      !Close(value.energyResidualJ,std::max(available,value.actualCOMKineticJ),controls.energyJ)||
      !Close(value.totalKineticJ-(comEnergy+available),std::max(value.totalKineticJ,comEnergy+available),controls.energyJ))
    return Fail("stored product velocity/momentum and COM/total kinetic energy fail the explicit numerical closure bounds",error);
  *output=std::move(value);if(error)error->clear();return true;
}
bool AddInteger(std::int64_t a,std::int64_t b,std::int64_t* result) {
  if((b>0&&a>std::numeric_limits<std::int64_t>::max()-b)||
      (b<0&&a<std::numeric_limits<std::int64_t>::min()-b))return false;
  *result=a+b;return true;
}
template<std::size_t N>
bool Material(const ResolvedState& parent,const std::array<ResolvedState,N>& products,std::string* error) {
  std::map<std::string,std::int64_t> incoming,outgoing;std::int64_t charge=0;
  for(const auto& n:parent.definition().nuclei)incoming[n.identity]=n.count;
  for(const auto& s:products) {
    if(!AddInteger(charge,s.definition().elementaryCharge,&charge))return Fail("integer product charge overflows",error);
    for(const auto& n:s.definition().nuclei)
      if(!AddInteger(outgoing[n.identity],n.count,&outgoing[n.identity]))return Fail("integer product nucleus count overflows",error);
  }
  if(incoming!=outgoing||charge!=parent.definition().elementaryCharge)
    return Fail("complete photon products do not conserve explicit isotope nuclei and integer charge",error);
  return true;
}
}

bool ResolveState(const StateDefinition& definition,const std::vector<IsotopeMass>& isotopes,
    const MassValue& electron,ResolvedState* output,std::string* error) {
  if(!output)return Fail("resolved-state output storage is required",error);
  try {
    if(definition.identity.empty()||definition.state.empty()||!Positive(electron.kg)||electron.sourceIdentity.empty()||
        !EnergyValid(definition.formation,false,error)||!EnergyValid(definition.excitation,true,error))
      return Fail("explicit isotope/state identity, electron reference and formation/excitation sources are required",error);
    std::map<std::string,MassValue> inventory;
    for(const auto& isotope:isotopes)
      if(isotope.identity.empty()||!Positive(isotope.neutralAtom.kg)||isotope.neutralAtom.sourceIdentity.empty()||
          !inventory.emplace(isotope.identity,isotope.neutralAtom).second)
        return Fail("isotope inventory has a missing/duplicate identity, invalid mass or missing source",error);
    if(definition.nuclei.empty()&&definition.elementaryCharge!=-1)
      return Fail("an empty nuclear composition is supported only for an explicitly charged single electron",error);
    ResolvedState value;value.definition_=definition;std::set<std::string> seen;
    for(const auto& nucleus:definition.nuclei) {
      const auto found=inventory.find(nucleus.identity);
      if(nucleus.count<=0||!seen.insert(nucleus.identity).second||found==inventory.end())
        return Fail("state has unknown/duplicate isotope identity or nonpositive integer multiplicity",error);
      value.mass_+=static_cast<Real>(nucleus.count)*found->second.kg;
      value.sources_.push_back(found->second.sourceIdentity);
    }
    value.mass_-=static_cast<Real>(definition.elementaryCharge)*electron.kg;
    if(!Positive(value.mass_)||!Finite(StateEnergy(value)))return Fail("additive reference mass/state energy is unrepresentable or nonpositive",error);
    if(definition.elementaryCharge!=0)value.sources_.push_back(electron.sourceIdentity);
    value.sources_.push_back(definition.formation.sourceIdentity);value.sources_.push_back(definition.excitation.sourceIdentity);
    std::ostringstream identity;identity.imbue(std::locale::classic());identity<<std::hexfloat;
    for(const auto& item:inventory)identity<<item.first.size()<<':'<<item.first<<'='<<item.second.kg<<'/'
        <<item.second.sourceIdentity.size()<<':'<<item.second.sourceIdentity<<';';
    identity<<"electron="<<electron.kg<<'/'<<electron.sourceIdentity.size()<<':'<<electron.sourceIdentity;
    value.inventory_=identity.str();value.valid_=true;*output=std::move(value);if(error)error->clear();return true;
  }
  catch(const std::exception& exception){return Fail(std::string("reference-state staging failed: ")+exception.what(),error);}
}

bool FormationQ(const std::vector<ResolvedState>& incoming,const std::vector<ResolvedState>& outgoing,
    FormationBalance* output,std::string* error) {
  if(!output||incoming.empty()||outgoing.empty())return Fail("formation Q needs explicit incoming/outgoing occurrences and output storage",error);
  try {
    FormationBalance value;const std::string inventory=incoming.front().referenceInventoryIdentity();
    EnergySum inputEnergy,outputEnergy,q;
    for(int side=0;side<2;++side)for(std::size_t index=0;index<(side==0?incoming:outgoing).size();++index) {
      const auto& state=(side==0?incoming:outgoing)[index];
      if(!StateValid(state,error)||state.referenceInventoryIdentity()!=inventory)return Fail("formation Q cannot mix unresolved states or reference inventories",error);
      const auto& definition=state.definition();auto& energy=side==0?inputEnergy:outputEnergy;
      for(Real component:{definition.formation.joules,definition.excitation.joules})
        if(!energy.Add(component)||!q.Add(side==0?component:-component))return Fail("state-energy expansion overflowed before publication",error);
      Sources(state,&value.sourceIdentities);
      value.terms.push_back({side==0,index,definition.identity,definition.state,definition.formation,definition.excitation});
    }
    value.incomingStateEnergyJ=inputEnergy.Value();value.outgoingStateEnergyJ=outputEnergy.Value();value.qJ=q.Value();
    if(!Finite(value.qJ)||!Finite(value.incomingStateEnergyJ)||!Finite(value.outgoingStateEnergyJ))return Fail("formation/excitation sum is unrepresentable",error);
    *output=std::move(value);if(error)error->clear();return true;
  }
  catch(const std::exception& exception){return Fail(std::string("formation bookkeeping staging failed: ")+exception.what(),error);}
}

bool EvaluateModelErrors(const ResolvedState& reference,const MassValue& physical,Real speed,
    const ErrorLimits& limits,ModelErrorReport* output,std::string* error) {
  if(!output||!StateValid(reference,error)||!Positive(physical.kg)||physical.sourceIdentity.empty()||
      !Finite(speed)||speed<0||speed>=SpeedOfLightMS)
    return Fail("model-error evaluation needs a resolved reference, source physical mass and finite subluminal speed",error);
  const Real values[]={limits.maximumRelativeInertialMassError,limits.maximumRelativeAccelerationError,
      limits.maximumRelativeRelativisticKineticError};
  for(Real limit:values)if(!Finite(limit)||limit<0||limit>=1)return Fail("independent model-error limits must be explicit in [0,1)",error);
  try {
    ModelErrorReport value;value.speedMS=speed;const Real difference=std::fabs(reference.referenceMassKg()-physical.kg);
    value.relativeInertialMassError=difference/physical.kg;
    value.relativeAccelerationError=difference/reference.referenceMassKg();
    const Real beta=speed/SpeedOfLightMS,s=std::sqrt((1-beta)*(1+beta));
    value.relativeRelativisticKineticError=beta*beta*(2+s)/(2*(1+s));
    if(speed>0&&value.relativeRelativisticKineticError==0)
      return Fail("positive relativistic model error underflows output precision; zero cannot be certified",error);
    if(!Finite(value.relativeInertialMassError)||!Finite(value.relativeAccelerationError)||
        !Finite(value.relativeRelativisticKineticError)||value.relativeInertialMassError>values[0]||
        value.relativeAccelerationError>values[1]||value.relativeRelativisticKineticError>values[2])
      return Fail("reference-inertia or nonrelativistic kinetic error exceeds an explicitly supplied model-error limit",error);
    Sources(reference,&value.sourceIdentities);value.sourceIdentities.push_back(physical.sourceIdentity);
    *output=std::move(value);if(error)error->clear();return true;
  }
  catch(const std::exception& exception){return Fail(std::string("model-error staging failed: ")+exception.what(),error);}
}

bool GeneralTwoBodyFinalState(const TwoBodyInput& input,FinalState* output,std::string* error) {
  if(!output||!Controls(input.context,input.controls,error)||!VectorFinite(input.totalMomentumKgMS)||
      !Finite(input.availableCOMKineticJ)||input.availableCOMKineticJ<0||
      !Unit(input.relativeDirection,input.controls.unitVectorTolerance,error))
    return Fail("invalid explicit two-body momentum, COM energy, direction or controls",error);
  try {
    for(const auto& s:input.products)if(!StateValid(s,error))return false;
    if(input.products[0].referenceInventoryIdentity()!=input.products[1].referenceInventoryIdentity())
      return Fail("two-body products use incompatible reference inventories",error);
    const Real a=input.products[0].referenceMassKg(),b=input.products[1].referenceMassKg(),mass=a+b;
    if(!Positive(mass))return Fail("two-body total mass is unrepresentable",error);
    const Real reduced=std::min(a,b)*(std::max(a,b)/mass);
    const Real relative=std::sqrt(reduced)*std::sqrt(input.availableCOMKineticJ)*std::sqrt(2.0L);
    if(!Finite(relative)||(relative==0&&input.availableCOMKineticJ>0))return Fail("positive two-body relative momentum is unrepresentable",error);
    std::vector<Vector3> momenta(2);
    const Real directionNorm=Norm(input.relativeDirection);
    for(int d=0;d<3;++d) {
      momenta[0][d]=(a/mass)*input.totalMomentumKgMS[d]+relative*(input.relativeDirection[d]/directionNorm);
      momenta[1][d]=(b/mass)*input.totalMomentumKgMS[d]-relative*(input.relativeDirection[d]/directionNorm);
    }
    return Finalize({input.products[0],input.products[1]},momenta,input.totalMomentumKgMS,
        input.availableCOMKineticJ,input.context,input.controls,"declared-two-body-final-state",output,error);
  }
  catch(const std::exception& exception){return Fail(std::string("two-body staging failed: ")+exception.what(),error);}
}

bool IsotropicDirectionFromOpenUniforms(Real polar,Real azimuthUniform,
    Vector3* output,std::string* error) {
  if(!output||!(polar>0&&polar<1)||!(azimuthUniform>0&&azimuthUniform<1))
    return Fail("isotropic direction requires two supplied strictly open uniforms and output storage",error);
  const Real cosine=2*polar-1,azimuth=2*std::acos(-1.0L)*azimuthUniform;
  // This preserves the original PhotonTwoBody formula/order. Reconstructing
  // sine from a rounded cosine would erase a tiny but representable recoil.
  const Real sine=2*std::sqrt(polar)*std::sqrt(1-polar);
  Vector3 value{{sine*std::cos(azimuth),sine*std::sin(azimuth),cosine}};
  if(!VectorFinite(value))return Fail("isotropic direction is unrepresentable",error);
  *output=value;if(error)error->clear();return true;
}

namespace {
// Both photon kernels compute exactly the same source threshold, state Q
// and incident photon momentum. Only the decomposition of available COM
// energy among the products depends on the selected channel.
template<std::size_t N,typename Solver>
bool ParentRestPhoton(const ParentRestPhotonInput<N>& input,const Solver& solve,
                      FinalState* output,std::string* error) {
  if(!output||!Controls(input.context,input.controls,error)||!StateValid(input.parent,error)||
      !VectorFinite(input.parentMomentumKgMS)||Norm(input.parentMomentumKgMS)!=0||
      !Unit(input.photonDirection,input.controls.unitVectorTolerance,error)||
      !EnergyValid(input.channelThreshold,true,error)||!EnergyValid(input.meanExcess,true,error))
    return Fail("photon model requires parent-rest momentum and explicit source energies and direction",error);
  try {
    for(const auto& s:input.products)if(!StateValid(s,error)||s.referenceInventoryIdentity()!=input.parent.referenceInventoryIdentity())
      return Fail("photon products require the same explicit reference inventory as the parent",error);
    if(!Material(input.parent,input.products,error))return false;
    Real mass=0;
    for(const auto& state:input.products)mass+=state.referenceMassKg();
    if(!Close(mass-input.parent.referenceMassKg(),std::max(mass,input.parent.referenceMassKg()),input.controls.referenceMassKg))
      return Fail("photon reference material mass is inconsistent; no artificial bath mass is supplied",error);
    FormationBalance formation;
    if(!FormationQ({input.parent},std::vector<ResolvedState>(input.products.begin(),input.products.end()),&formation,error))return false;
    if(input.energyPolicy!=PhotonEnergyPolicy::SourceEqualsState&&
       input.energyPolicy!=PhotonEnergyPolicy::SourceEffectivePhoton)
      return Fail("unsupported explicit photon energy policy",error);
    if(input.energyPolicy==PhotonEnergyPolicy::SourceEffectivePhoton&&input.sourceThresholdCrosswalkIdentity.empty())
      return Fail("source-effective photon requires an explicit source/physical-state crosswalk identity",error);
    if(input.energyPolicy==PhotonEnergyPolicy::SourceEqualsState&&
       !Close(-formation.qJ-input.channelThreshold.joules,std::max(std::fabs(formation.qJ),input.channelThreshold.joules),input.controls.thresholdJ))
      return Fail("source threshold and explicit formation/excitation state crosswalk disagree",error);
    const Real photon=input.channelThreshold.joules+input.meanExcess.joules;
    if(!Positive(photon))return Fail("representative incoming photon energy must be finite positive",error);
    const Real photonDirectionNorm=Norm(input.photonDirection);
    Vector3 p;for(int d=0;d<3;++d)p[d]=(photon/SpeedOfLightMS)*(input.photonDirection[d]/photonDirectionNorm);
    Real com=0;if(!MomentumKinetic(p,mass,&com))return Fail("photon COM recoil energy is unrepresentable",error);
    // Evaluate the source-state cancellation BEFORE adding the small mean
    // excess. Subtracting two rounded large photon/state totals can erase a
    // genuine, separately supplied translational energy. The receipt retains
    // those components; a later packet/bath adapter must preserve them too.
    EnergySum sourceAvailable;
    if(!sourceAvailable.Add(input.meanExcess.joules)||!sourceAvailable.Add(input.channelThreshold.joules))
      return Fail("source available-energy expansion overflowed",error);
    for(const auto& term:formation.terms)
      for(Real component:{term.formation.joules,term.excitation.joules})
        if(!sourceAvailable.Add(term.incoming?component:-component))return Fail("source-state available-energy expansion overflowed",error);
    const Real stateAvailable=sourceAvailable.Value();
    const Real available=stateAvailable-com;
    if(!Finite(available)||available<0)return Fail("photon/state energies leave negative or unrepresentable COM translation; no clipping is allowed",error);
    FinalState value;
    if(!solve(p,available,&value))return false;
    value.photonIncluded=true;
    value.photonEnergyJ=photon;value.photonMomentumKgMS=p;value.stateQJ=formation.qJ;
    value.channelThresholdJ=input.channelThreshold.joules;value.sourceMeanExcessJ=input.meanExcess.joules;
    value.stateAvailableKineticJ=stateAvailable;
    value.formationBalance=formation;
    value.sourceIdentities.insert(value.sourceIdentities.end(),formation.sourceIdentities.begin(),formation.sourceIdentities.end());
    value.sourceIdentities.push_back(input.channelThreshold.sourceIdentity);value.sourceIdentities.push_back(input.meanExcess.sourceIdentity);
    if(input.energyPolicy==PhotonEnergyPolicy::SourceEffectivePhoton)
      value.sourceIdentities.push_back(input.sourceThresholdCrosswalkIdentity);
    if(!Close(value.totalKineticJ-stateAvailable,std::max(value.totalKineticJ,std::fabs(stateAvailable)),input.controls.energyJ))
      return Fail("complete photon and explicit state-energy ledger failed the numerical closure bound",error);
    *output=std::move(value);if(error)error->clear();return true;
  }
  catch(const std::exception& exception){return Fail(std::string("photon staging failed: ")+exception.what(),error);}
}
} // namespace

bool PhotonTwoBody(const PhotonTwoBodyInput& input,FinalState* output,std::string* error) {
  auto solve=[&](const Vector3& p,Real available,FinalState* value) {
    TwoBodyInput final;final.products=input.products;final.totalMomentumKgMS=p;
    if(!IsotropicDirectionFromOpenUniforms(input.polarOpenUniform,input.azimuthOpenUniform,
        &final.relativeDirection,error))return false;
    final.availableCOMKineticJ=available;final.context=input.context;final.controls=input.controls;
    if(!GeneralTwoBodyFinalState(final,value,error))return false;
    value->modelIdentity=input.energyPolicy==PhotonEnergyPolicy::SourceEqualsState?
      "source-mean-isotropic-parent-rest-photon-two-body":
      "source-effective-mean-isotropic-parent-rest-photon-two-body";
    return true;
  };
  return ParentRestPhoton<2>(input,solve,output,error);
}

bool PhotonChargedBreakup(const PhotonChargedBreakupInput& input,FinalState* output,std::string* error) {
  auto solve=[&](const Vector3& p,Real available,FinalState* value) {
    ChargedBreakupInput final;
    static_cast<ChargedPartition&>(final)=input.partition;
    final.electron=input.products[0];final.ion=input.products[1];final.neutral=input.products[2];
    final.totalMomentumKgMS=p;final.availableCOMKineticJ=available;
    final.context=input.context;final.controls=input.controls;
    if(!ChargedBreakupFinalState(final,value,error))return false;
    value->modelIdentity=input.energyPolicy==PhotonEnergyPolicy::SourceEqualsState?
      "source-mean-parent-rest-photon-charged-breakup":"source-effective-mean-parent-rest-photon-charged-breakup";
    return true;
  };
  return ParentRestPhoton<3>(input,solve,output,error);
}

bool ChargedBreakupFinalState(const ChargedBreakupInput& input,FinalState* output,std::string* error) {
  if(!output||!Controls(input.context,input.controls,error)||input.partitionSourceIdentity.empty()||
     !VectorFinite(input.totalMomentumKgMS)||!Finite(input.availableCOMKineticJ)||input.availableCOMKineticJ<0||
     !Finite(input.electronRelativeEnergyFraction)||input.electronRelativeEnergyFraction<0||
     input.electronRelativeEnergyFraction>1)
    return Fail("charged breakup requires explicit finite momentum, COM energy, frame and sourced fraction in [0,1]",error);
  try {
    const std::vector<ResolvedState> states{input.electron,input.ion,input.neutral};
    for(const auto& state:states)
      if(!StateValid(state,error)||state.referenceInventoryIdentity()!=input.electron.referenceInventoryIdentity())
        return Fail("charged breakup products must share a resolved reference inventory",error);
    const auto& electron=input.electron.definition();
    const auto& ion=input.ion.definition();
    const auto& neutral=input.neutral.definition();
    if(electron.elementaryCharge!=-1||!electron.nuclei.empty()||ion.elementaryCharge!=1||ion.nuclei.empty()||
       neutral.elementaryCharge!=0||neutral.nuclei.empty())
      return Fail("charged breakup requires ordered electron, singly charged nuclear ion and neutral fragment",error);

    const Real me=input.electron.referenceMassKg(),mi=input.ion.referenceMassKg(),mn=input.neutral.referenceMassKg();
    const Real pair=mi+mn,total=pair+me;
    if(!Positive(total))return Fail("charged breakup total mass is unrepresentable",error);
    const Real fraction=input.electronRelativeEnergyFraction;

    // Jacobi relative momentum has magnitude sqrt(2 mu K f). Keep the square
    // roots separate, so a small fraction multiplied by a small K cannot
    // silently remove a representable relative momentum. At zero mode energy
    // there is no distinguished direction and a finite zero vector is valid.
    auto modeMomentum=[&](Real first,Real second,Real share,const Vector3& direction,Vector3& p) {
      if(!VectorFinite(direction))return Fail("nonfinite charged-breakup mode direction",error);
      if(share==0||input.availableCOMKineticJ==0)return true;
      if(!Unit(direction,input.controls.unitVectorTolerance,error))return false;
      const Real mu=std::min(first,second)*(std::max(first,second)/(first+second));
      const Real magnitude=std::sqrt(mu)*std::sqrt(input.availableCOMKineticJ)*std::sqrt(2*share);
      if(!Positive(magnitude))return Fail("charged-breakup mode momentum is unrepresentable",error);
      const Real directionNorm=Norm(direction);
      for(int d=0;d<3;++d)p[d]=magnitude*(direction[d]/directionNorm);
      return true;
    };
    Vector3 eRelative{},heavyRelative{};
    if(!modeMomentum(me,pair,fraction,input.electronRelativeToPair,eRelative)||
       !modeMomentum(mi,mn,1-fraction,input.ionRelativeToNeutral,heavyRelative))return false;
    std::vector<Vector3> momenta(3);
    for(int d=0;d<3;++d) {
      // The heavy-pair COM carries the electron's opposite momentum in the
      // total COM. Share that COM momentum by mass, then add equal and
      // opposite ion–neutral recoil. All three also receive the bulk boost.
      momenta[0][d]=(me/total)*input.totalMomentumKgMS[d]+eRelative[d];
      const Real pairMomentum=(pair/total)*input.totalMomentumKgMS[d]-eRelative[d];
      momenta[1][d]=(mi/pair)*pairMomentum+heavyRelative[d];
      momenta[2][d]=(mn/pair)*pairMomentum-heavyRelative[d];
    }
    FinalState staged;
    if(!Finalize(states,momenta,input.totalMomentumKgMS,input.availableCOMKineticJ,input.context,
                 input.controls,"charged-breakup-explicit-jacobi", &staged,error))return false;
    // A small individual mode can disappear in a large boost even though
    // the *total* energy and momentum pass. Verify both modes independently
    // from the delivered velocities before publishing any product.
    const auto& ve=staged.products[0].velocityMS;
    const auto& vi=staged.products[1].velocityMS;
    const auto& vn=staged.products[2].velocityMS;
    const Real electronMu=std::min(me,pair)*(std::max(me,pair)/total);
    const Real heavyMu=std::min(mi,mn)*(std::max(mi,mn)/pair);
    Vector3 electronError{},heavyError{};
    for(int d=0;d<3;++d) {
      const Real pairVelocity=(mi/pair)*vi[d]+(mn/pair)*vn[d];
      electronError[d]=electronMu*(ve[d]-pairVelocity)-eRelative[d];
      heavyError[d]=heavyMu*(vi[d]-vn[d])-heavyRelative[d];
    }
    // At exactly zero mode energy there is no angle or fractional scale.
    // Its roundoff scale is the largest delivered material momentum.
    const Real scale=std::max({Norm(momenta[0]),Norm(momenta[1]),Norm(momenta[2])});
    if(!Close(Norm(electronError),Norm(eRelative)>0?Norm(eRelative):scale,input.controls.momentumKgMS)||
       !Close(Norm(heavyError),Norm(heavyRelative)>0?Norm(heavyRelative):scale,input.controls.momentumKgMS))
      return Fail("charged-breakup mode momentum was lost at the requested precision",error);
    staged.sourceIdentities.push_back(input.partitionSourceIdentity);
    *output=std::move(staged);
    if(error)error->clear();
    return true;
  }
  catch(const std::exception& exception){return Fail(std::string("charged breakup staging failed: ")+exception.what(),error);}
}

bool ThreeNeutralMicrocanonical(const ThreeBodyInput& input,FinalState* output,std::string* error) {
  if(!output||!Controls(input.context,input.controls,error)||!VectorFinite(input.totalMomentumKgMS)||
      !Finite(input.availableCOMKineticJ)||input.availableCOMKineticJ<0)
    return Fail("three-body model requires explicit finite total momentum, nonnegative COM energy and controls",error);
  try {
    Real mass=0;std::size_t dependent=0;
    for(std::size_t i=0;i<3;++i) {
      if(!StateValid(input.products[i],error)||input.products[i].definition().elementaryCharge!=0||
          input.products[i].referenceInventoryIdentity()!=input.products[0].referenceInventoryIdentity()||!VectorFinite(input.standardNormals[i]))
        return Fail("microcanonical model requires three neutral same-inventory states and finite supplied normals",error);
      mass+=input.products[i].referenceMassKg();
      if(input.products[i].referenceMassKg()>input.products[dependent].referenceMassKg())dependent=i;
    }
    if(!Positive(mass))return Fail("three-body total reference mass is unrepresentable",error);
    std::vector<Vector3> relative(3);Vector3 drawMomentum{{0,0,0}};
    if(input.availableCOMKineticJ>0) {
      for(std::size_t i=0;i<3;++i)for(int d=0;d<3;++d) {
        relative[i][d]=std::sqrt(input.products[i].referenceMassKg())*input.standardNormals[i][d];
        drawMomentum[d]+=relative[i][d];
      }
      for(std::size_t i=0;i<3;++i)if(i!=dependent)for(int d=0;d<3;++d)
        relative[i][d]-=(input.products[i].referenceMassKg()/mass)*drawMomentum[d];
      relative[dependent]={{0,0,0}};
      for(std::size_t i=0;i<3;++i)if(i!=dependent)for(int d=0;d<3;++d)relative[dependent][d]-=relative[i][d];
      Real sampled=0;
      for(std::size_t i=0;i<3;++i){Real energy=0;if(!MomentumKinetic(relative[i],input.products[i].referenceMassKg(),&energy))return Fail("microcanonical draw energy is unrepresentable",error);sampled+=energy;}
      if(!Positive(sampled))return Fail("positive-energy microcanonical draw is degenerate; caller must not redraw an accepted event",error);
      const Real scale=std::sqrt(input.availableCOMKineticJ)/std::sqrt(sampled);
      for(auto& p:relative)for(auto& component:p)component*=scale;
      relative[dependent]={{0,0,0}};
      for(std::size_t i=0;i<3;++i)if(i!=dependent)for(int d=0;d<3;++d)relative[dependent][d]-=relative[i][d];
    }
    std::vector<Vector3> momenta(3);
    for(std::size_t i=0;i<3;++i)for(int d=0;d<3;++d)
      momenta[i][d]=(input.products[i].referenceMassKg()/mass)*input.totalMomentumKgMS[d]+relative[i][d];
    return Finalize({input.products[0],input.products[1],input.products[2]},momenta,input.totalMomentumKgMS,
        input.availableCOMKineticJ,input.context,input.controls,"mass-weighted-microcanonical-three-neutral",output,error);
  }
  catch(const std::exception& exception){return Fail(std::string("microcanonical staging failed: ")+exception.what(),error);}
}
}}}
