#include "IonCollector.h"
#include <algorithm>
#include <cmath>
#include <exception>
#include <limits>
#include <type_traits>
#include <utility>

namespace Moon { namespace Species { namespace IonCollector {
namespace {
bool Reject(std::string* error,const char* text) { if(error)*error=text;return false; }
bool Finite3(const Vector3& v) { for(double x:v)if(!std::isfinite(x))return false;return true; }
bool FiniteExtended(long double x) { return std::isfinite(x); }
}

bool ValidateStates(const std::vector<StateDefinition>& states,std::string* error) {
  if(!error)return false;
  if(states.empty())return Reject(error,"ideal ion collector requires an explicit nonempty state inventory");
  std::vector<R::SpeciesIdentity> material;
  material.reserve(states.size());
  for(const auto& state:states) {
    if(state.material.compiledSpecies<0 || state.material.elementaryCharge==0 || state.material.nuclei.empty() ||
       !std::isfinite(state.formationEnergyJ) || !std::isfinite(state.excitationEnergyJ) ||
       state.excitationEnergyJ<0 || state.energyReferenceIdentity.empty() || state.sourceEvidenceIdentity.empty())
      return Reject(error,"collector state needs a compiled charged species, explicit finite formation/excitation energies, named energy zero and source evidence");
    if(state.energyReferenceIdentity!=states.front().energyReferenceIdentity)
      return Reject(error,"one collector inventory cannot sum energies with different reference zeros");
    material.push_back(state.material);
  }
  // Reuse the production integer-composition/state validator. A rounded mass
  // is never interpreted as nuclear composition and no default ground state
  // is generated when several catalog states map to one compiled species.
  if(!R::ValidateCatalog(material,error))return false;
  error->clear();return true;
}

bool Prepare(const Incident& incident,Record* output,std::string* error) {
  if(!output || !error)return false;
  if(!ValidateStates({incident.state},error))return false;
  if(!(incident.representedNumber>0) || !std::isfinite(incident.representedNumber) ||
     !(incident.elementaryChargeC>0) || !std::isfinite(incident.elementaryChargeC) ||
     !std::isfinite(incident.impactEpochEtS) || !Finite3(incident.positionLsoM) ||
     !Finite3(incident.velocityLsoMS) || !Finite3(incident.frame.outwardNormalLso) ||
     !Finite3(incident.frame.surfaceVelocityLsoMS) || incident.frame.geometryIdentity.empty() ||
     !std::isfinite(incident.frame.geometryEpochEtS) || !incident.frame.stateTransformsSupplied)
    return Reject(error,"collector impact requires actual finite weight/epoch/position/velocity and explicitly supplied geometry/frame transforms");
  for(double x:incident.frame.lsoToBody36)if(!std::isfinite(x))return Reject(error,"invalid collector impact state transform");
  for(double x:incident.frame.bodyToLso36)if(!std::isfinite(x))return Reject(error,"invalid collector impact inverse state transform");
  // A singular/zero array with a claimed supplied flag is not a transform.
  // Check the supplied inverse pair using a floating-point forward bound,
  // not a tunable physical geometry tolerance. Actual frame provenance and
  // epoch remain caller responsibilities; this is arithmetic validation only.
  const long double roundoff=64*std::numeric_limits<double>::epsilon();
  for(int i=0;i<6;++i)for(int j=0;j<6;++j) {
    long double value=0,magnitude=0;
    for(int k=0;k<6;++k) {
      const long double term=static_cast<long double>(incident.frame.lsoToBody36[6*i+k])*incident.frame.bodyToLso36[6*k+j];
      value+=term;magnitude+=std::fabs(term);
    }
    const long double expected=i==j?1:0;
    if(std::fabs(value-expected)>roundoff*(magnitude+expected))
      return Reject(error,"collector impact state transforms are not an inverse pair");
  }
  long double normal2=0,normalVelocity=0;
  for(int d=0;d<3;++d) {
    const long double n=incident.frame.outwardNormalLso[d];normal2+=n*n;
    normalVelocity+=(static_cast<long double>(incident.velocityLsoMS[d])-incident.frame.surfaceVelocityLsoMS[d])*n;
  }
  if(!(normal2>0) || !FiniteExtended(normal2) || normalVelocity>0)
    return Reject(error,"collector needs a nonzero outward normal and a non-outgoing impact in the actual facet rest frame");
  Record record;record.incident=incident;
  const long double number=incident.representedNumber;
  const long double mass=incident.state.material.modelMassKg;
  record.modelMassKg=number*mass;
  record.representedChargeElementary=number*incident.state.material.elementaryCharge;
  record.chargeC=record.representedChargeElementary*incident.elementaryChargeC;
  for(const auto& n:incident.state.material.nuclei)
    record.nuclei.push_back({n.identity,n.count,number*n.count});
  for(int d=0;d<3;++d) {
    const long double v=incident.velocityLsoMS[d];
    const long double relative=v-incident.frame.surfaceVelocityLsoMS[d];
    record.momentumLsoKgMS[d]=record.modelMassKg*v;
    record.momentumFacetKgMS[d]=record.modelMassKg*relative;
    record.kineticEnergyLsoJ+=record.modelMassKg*v*v/2;
    record.kineticEnergyFacetJ+=record.modelMassKg*relative*relative/2;
  }
  record.formationEnergyJ=number*incident.state.formationEnergyJ;
  record.excitationEnergyJ=number*incident.state.excitationEnergyJ;
  record.totalEnergyFacetJ=record.kineticEnergyFacetJ+record.formationEnergyJ+record.excitationEnergyJ;
  // Long-double accumulation avoids needlessly losing small explicit internal
  // contributions. Refuse overflow rather than silently emitting an infinite
  // transferred material/charge/energy budget. Signed formation energy is valid.
  for(long double x:{record.modelMassKg,record.representedChargeElementary,record.chargeC,
      record.kineticEnergyLsoJ,record.kineticEnergyFacetJ,record.formationEnergyJ,
      record.excitationEnergyJ,record.totalEnergyFacetJ})
    if(!FiniteExtended(x))return Reject(error,"collector transfer overflow; reduce represented packet weight or use a smaller explicit physical input domain");
  for(int d=0;d<3;++d)if(!FiniteExtended(record.momentumLsoKgMS[d]) || !FiniteExtended(record.momentumFacetKgMS[d]))
    return Reject(error,"collector momentum transfer overflow");
  for(const auto& n:record.nuclei)if(!FiniteExtended(n.representedNumber))return Reject(error,"collector represented nuclei overflow");
  *output=std::move(record);error->clear();return true;
}

bool Collect(const Incident& incident,Ledger* ledger,std::string* error) {
  if(!ledger || !error)return false;
  try {
    Record record;
    if(!Prepare(incident,&record,error))return false;
    if(!ledger->records.empty() && ledger->records.front().incident.state.energyReferenceIdentity!=incident.state.energyReferenceIdentity)
      return Reject(error,"collector ledger cannot combine transfers expressed with different formation-energy zeros");
    if(ledger->records.size()==ledger->records.max_size())
      return Reject(error,"ideal collector ledger exhausted; stop the run and reduce particles or export/aggregate bounded records before restarting");
    if(ledger->records.size()==ledger->records.capacity()) {
      // Geometric capacity growth avoids copying the entire impact history
      // for every event (quadratic work). Reservation remains before any
      // physical publication; overflow/allocation failure leaves the incident
      // packet and ledger unchanged. This is storage policy, not a yield or
      // a limit on the physical number of collected particles.
      const auto capacity=ledger->records.capacity();
      const auto maximum=ledger->records.max_size();
      const auto grown=capacity>maximum-capacity?maximum:capacity*2;
      ledger->records.reserve(std::max(ledger->records.size()+1,grown));
    }
    // No model packet mutation occurs here. The PIC adapter returns Delete
    // only after this append succeeded; allocation failure leaves both packet
    // and physical ledger unchanged. Moving an allocated record cannot throw.
    static_assert(std::is_nothrow_move_constructible<Record>::value,
        "collector publication must not allocate after ledger reservation");
    ledger->records.push_back(std::move(record));
  } catch(const std::exception&) {
    return Reject(error,"ideal collector record allocation failed; run terminated to preserve material/charge/energy accounting; reduce packet count or provide more memory");
  }
  error->clear();return true;
}
} } }
