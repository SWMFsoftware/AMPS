#include "SodiumPartialCapture.h"
#include "Photochemistry.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <sstream>

namespace Moon { namespace Species { namespace SodiumPartialCapture {
namespace {
const char* const FrozenNodesSha256="2a887ad43497a7f9e2b3d7b6173ae1b2597cc893d897926762c181946c92db93";
bool Error(std::string* error,const std::string& message){if(error)*error=message;return false;}
bool Positive(Real x){return std::isfinite(x)&&x>0;}
bool ConstantSource(const Constants& c){return Positive(c.atomicMassUnitKg)&&Positive(c.electronVoltJ)&&!c.sourceIdentity.empty();}
bool Near(Real a,Real b,const K::Bound& tolerance){
  return std::isfinite(tolerance.absolute)&&std::isfinite(tolerance.relative)&&
      tolerance.absolute>=0&&tolerance.relative>=0&&
      std::abs(a-b)<=tolerance.absolute+tolerance.relative*std::max(std::abs(a),std::abs(b));
}
bool SameNuclei(const K::ResolvedState& a,const K::ResolvedState& b){
  const auto& x=a.definition().nuclei;const auto& y=b.definition().nuclei;
  if(x.size()!=y.size())return false;
  for(const auto& n:x){auto it=std::find_if(y.begin(),y.end(),[&](const ReactionExecution::NucleusCount& m){return m.identity==n.identity;});if(it==y.end()||it->count!=n.count)return false;}
  return true;
}
}

bool Registry::LoadNativeNodes(const std::string& path,std::string* error){
  // Hash verification precedes parsing; an input-supplied hash cannot silently
  // bless a different channel under this registry's fixed scientific meaning.
  if(!ReactionData::VerifyFileSha256(path,FrozenNodesSha256,error))return false;
  std::ifstream stream(path);std::string line;
  if(!std::getline(stream,line))return Error(error,"missing Na partial capture CSV header");
  // The retained native numeric extraction uses CSV CRLF records. Removing
  // the record terminator during parsing does not modify source bytes/data.
  if(!line.empty()&&line.back()=='\r')line.pop_back();
  if(line!="energy_eV_per_u,cross_section_cm2")
    return Error(error,"Na partial capture native CSV header mismatch");
  std::vector<Node> candidate;
  while(std::getline(stream,line)){
    if(!line.empty()&&line.back()=='\r')line.pop_back();
    if(line.empty())return Error(error,"blank Na partial capture source record");
    const auto comma=line.find(',');
    if(comma==std::string::npos||line.find(',',comma+1)!=std::string::npos)
      return Error(error,"malformed Na partial capture source record");
    std::istringstream a(line.substr(0,comma)),b(line.substr(comma+1));
    Node node;Real sourceCm2=0;
    if(!(a>>node.energyEvPerU)||!(b>>sourceCm2))return Error(error,"nonnumeric Na partial capture source record");
    a>>std::ws;b>>std::ws;
    if(!a.eof()||!b.eof()||!Positive(node.energyEvPerU)||!Positive(sourceCm2))
      return Error(error,"invalid Na partial capture source value");
    // cm² -> m² is an exact dimensional conversion, not a calibration.
    node.crossSectionM2=sourceCm2*1.0e-4L;
    if(!Positive(node.crossSectionM2)||(!candidate.empty()&&node.energyEvPerU<=candidate.back().energyEvPerU))
      return Error(error,"invalid or unordered Na partial capture native node");
    candidate.push_back(node);
  }
  if(!stream.eof()||candidate.size()!=12||candidate.front().energyEvPerU!=4000||candidate.back().energyEvPerU!=25000)
    return Error(error,"incomplete Na ground-product partial capture source");
  nodes_.swap(candidate);source_="IAEA-D119990-native-nodes-sha256:"+std::string(FrozenNodesSha256);
  return true;
}

bool Registry::RateAtNativeNode(Real energy,Real density,const Constants& constants,
    Rate* output,std::string* error)const{
  if(!output||!loaded()||!Positive(energy)||!std::isfinite(density)||density<0||!ConstantSource(constants))
    return Error(error,"partial capture requires loaded source, explicit constants, positive energy and nonnegative proton density");
  const auto it=std::find_if(nodes_.begin(),nodes_.end(),[&](const Node& n){return n.energyEvPerU==energy;});
  if(it==nodes_.end())return Error(error,"partial capture energy is not a native D119990 node; no interpolation or extrapolation is implemented");
  Rate candidate;candidate.constants_=constants;candidate.protonDensityM3_=density;candidate.energyEvPerU_=energy;candidate.crossSectionM2_=it->crossSectionM2;
  // E[eV/u] = (m_projectile*g²/2)/(m_projectile/u)/e,
  // hence g²=2*E[eV/u]*e/u, with e[J/eV] and u[kg].
  // The source uses target-rest projectile energy per u; the projectile mass
  // cancels here but is retained independently in momentum/energy kinematics.
  candidate.relativeSpeedMS_=std::sqrt(2*energy*constants.electronVoltJ/constants.atomicMassUnitKg);
  candidate.rateS_=density*candidate.crossSectionM2_*candidate.relativeSpeedMS_;
  candidate.sourceIdentity_=source_+";constants:"+constants.sourceIdentity;
  if(!Positive(candidate.relativeSpeedMS_)||!std::isfinite(candidate.rateS_))
    return Error(error,"Na partial capture rate overflow or unrepresentable source speed");
  candidate.valid_=true;*output=candidate;return true;
}

bool ForwardProducts(const Registry& registry,const ForwardInput& input,K::FinalState* output,std::string* error){
  const auto& p=input.incomingProton;const auto& a=input.incomingSodium;
  const auto& h=input.outgoingHydrogen;const auto& ion=input.outgoingSodiumIon;
  if(!output||!input.sourceRate.valid()||!p.valid()||!a.valid()||!h.valid()||!ion.valid()||!ConstantSource(input.constants)||
      !Positive(input.sourceRate.relativeSpeedMS())||!Positive(input.sourceRate.energyEvPerU())||
      input.sourceRate.sourceIdentity().empty()||!registry.loaded())return Error(error,"partial capture products require explicit resolved states and native source receipt");
  const auto node=std::find_if(registry.nodes().begin(),registry.nodes().end(),[&](const Node& n){return n.energyEvPerU==input.sourceRate.energyEvPerU();});
  if(input.constants.atomicMassUnitKg!=input.sourceRate.constants().atomicMassUnitKg||
      input.constants.electronVoltJ!=input.sourceRate.constants().electronVoltJ||
      input.constants.sourceIdentity!=input.sourceRate.constants().sourceIdentity||
      input.context.frameIdentity!="source-target-rest-inertial"||node==registry.nodes().end()||node->crossSectionM2!=input.sourceRate.crossSectionM2()||
      input.sourceRate.sourceIdentity()!=registry.sourceIdentity()+";constants:"+input.constants.sourceIdentity||
      !std::isfinite(input.sourceRate.rateS())||input.sourceRate.rateS()<0)
    return Error(error,"partial capture products require an authenticated native-node receipt from the same registry/constants");
  // The physical source states must be supplied, not inferred from a compiled
  // PIC species. Names are a versioned semantic crosswalk for this one channel.
  if(p.definition().identity!="H+"||a.definition().identity!="Na"||h.definition().identity!="H"||ion.definition().identity!="Na+"||
      p.definition().state!="bare-proton"||a.definition().state!="3s"||h.definition().state!="1s"||ion.definition().state!="2p6"||
      p.definition().elementaryCharge!=1||a.definition().elementaryCharge!=0||h.definition().elementaryCharge!=0||ion.definition().elementaryCharge!=1||
      !SameNuclei(p,h)||!SameNuclei(a,ion))return Error(error,"D119990 partial capture requires exactly H+/Na(3s) -> H(1s)/Na+(2p6) with preserved isotope identities");
  if(p.definition().nuclei.size()!=1||a.definition().nuclei.size()!=1||
      p.definition().nuclei.front().identity!="1H"||p.definition().nuclei.front().count!=1||
      a.definition().nuclei.front().identity!="23Na"||a.definition().nuclei.front().count!=1)
    return Error(error,"D119990 reference scope is a proton and 23Na; no isotope substitution is inferred");
  const Real massIn=p.referenceMassKg()+a.referenceMassKg(),massOut=h.referenceMassKg()+ion.referenceMassKg();
  if(!Near(massIn,massOut,input.controls.referenceMassKg))return Error(error,"partial capture reference masses do not close");
  Real norm2=0;for(Real x:input.relativeIncidentDirection){if(!std::isfinite(x))return Error(error,"nonfinite partial capture incident direction");norm2+=x*x;}
  if(!std::isfinite(input.controls.unitVectorTolerance)||input.controls.unitVectorTolerance<0||input.controls.unitVectorTolerance>=1||
      !Positive(norm2)||std::abs(norm2-1)>input.controls.unitVectorTolerance)
    return Error(error,"partial capture incident direction must be an explicit unit vector");
  const Real expectedSpeed=std::sqrt(2*input.sourceRate.energyEvPerU()*input.constants.electronVoltJ/input.constants.atomicMassUnitKg);
  // A finite numerical bound is explicit; a caller may not substitute a
  // convenient arbitrary wind speed for the source's energy-per-u convention.
  const Real speedError=std::abs(expectedSpeed-input.sourceRate.relativeSpeedMS());
  if(!std::isfinite(expectedSpeed)||speedError>input.controls.unitVectorTolerance*expectedSpeed)
    return Error(error,"partial capture rate energy and speed receipts disagree");
  K::FormationBalance balance;
  if(!K::FormationQ({p,a},{h,ion},&balance,error))return false;
  const Real incomingK=.5L*p.referenceMassKg()*expectedSpeed*expectedSpeed;
  // An accepted near-unit input carries numerical orientation error only;
  // normalization must not change the authenticated projectile speed/energy.
  K::Vector3 direction=input.relativeIncidentDirection;
  for(Real& value:direction)value/=std::sqrt(norm2);
  K::TwoBodyInput products;products.products={{h,ion}};
  for(unsigned k=0;k<3;++k)products.totalMomentumKgMS[k]=p.referenceMassKg()*expectedSpeed*direction[k];
  Real totalP2=0;for(Real value:products.totalMomentumKgMS)totalP2+=value*value;
  products.availableCOMKineticJ=incomingK+balance.qJ-totalP2/(2*massOut);
  products.relativeDirection=direction;products.context=input.context;products.controls=input.controls;
  K::FinalState candidate;
  if(!K::GeneralTwoBodyFinalState(products,&candidate,error))return false;
  candidate.modelIdentity="EXPERIMENTAL-D119990-ground-partial-forward-COM";
  candidate.stateQJ=balance.qJ;
  candidate.sourceIdentities.push_back(input.sourceRate.sourceIdentity());
  candidate.sourceIdentities.push_back("Na3s-approved-reference-cohort;native-initial-Na-state-unspecified;physical-state-crosswalk-NOT-QUALIFIED");
  candidate.sourceIdentities.insert(candidate.sourceIdentities.end(),balance.sourceIdentities.begin(),balance.sourceIdentities.end());
  *output=candidate;return true;
}
} } }
