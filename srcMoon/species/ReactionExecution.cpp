#include "ReactionExecution.h"
#include "PiecewiseHazard.h" // immutable exact-location receipt; no solver link dependency

// This upstream header is unit/PIC independent. Only its structural table
// validator is reused: Moon retains integer material terms, explicit baths
// and its own transient-number ledger. In particular we do not misuse the
// upstream diagnostic's documented number-RATE moments as number inventories.
#include "src/models/chemical_solver/ReactionNetwork.h"

#include <algorithm>
#include <cmath>
#include <exception>
#include <limits>
#include <map>
#include <set>
#include <utility>

namespace Moon { namespace Species { namespace ReactionExecution {
namespace {
constexpr double SpeedOfLightMS=299792458.0; // exact SI; not a recoil parameter
bool Fail(const std::string& message,std::string* error) {
  if (error) *error=message;
  return false;
}
bool Finite(double value) { return std::isfinite(value); }
bool Positive(double value) { return Finite(value) && value>0; }
bool VectorFinite(const Vector3& v) {
  return Finite(v[0]) && Finite(v[1]) && Finite(v[2]);
}
bool Add(std::int64_t a,std::int64_t b,std::int64_t* out) {
  const auto hi=std::numeric_limits<std::int64_t>::max();
  const auto lo=std::numeric_limits<std::int64_t>::min();
  if ((b>0 && a>hi-b) || (b<0 && a<lo-b)) return false;
  *out=a+b;return true;
}
bool Multiply(std::int64_t a,std::uint32_t b,std::int64_t* out) {
  if (b==0) {*out=0;return true;}
  const auto hi=std::numeric_limits<std::int64_t>::max();
  const auto lo=std::numeric_limits<std::int64_t>::min();
  if ((a>0 && a>hi/b) || (a<0 && a<lo/b)) return false;
  *out=a*static_cast<std::int64_t>(b);return true;
}
bool Subtract(std::int64_t a,std::int64_t b,std::int64_t* out) {
  const auto hi=std::numeric_limits<std::int64_t>::max();
  const auto lo=std::numeric_limits<std::int64_t>::min();
  if ((b>0 && a<lo+b) || (b<0 && a>hi+b)) return false;
  *out=a-b;return true;
}
using Composition=std::map<std::string,std::int64_t>;
struct Material {
  Composition nuclei;
  std::int64_t charge=0,trackedCharge=0,bathCharge=0;
  long double trackedMass=0,bathMass=0;
};
bool Accumulate(const std::vector<SpeciesIdentity>& catalog,
    const std::vector<Term>& terms,Material* material,std::string* error) {
  Material value;
  for (const Term& term:terms) {
    if (term.species>=catalog.size() || term.multiplicity==0)
      return Fail("reaction term has invalid species or zero integer multiplicity",error);
    if (term.population!=Population::Tracked && term.population!=Population::Bath)
      return Fail("reaction term has an invalid tracked/bath population",error);
    const auto& species=catalog[term.species];
    if (term.population==Population::Tracked && species.compiledSpecies<0)
      return Fail("tracked reaction term has no compiled species: "+species.identity,error);
    for (const auto& nucleus:species.nuclei) {
      std::int64_t product,next;
      if (!Multiply(nucleus.count,term.multiplicity,&product) ||
          !Add(value.nuclei[nucleus.identity],product,&next))
        return Fail("integer nucleus count overflows int64",error);
      value.nuclei[nucleus.identity]=next;
    }
    std::int64_t charge,next;
    if (!Multiply(species.elementaryCharge,term.multiplicity,&charge) ||
        !Add(value.charge,charge,&next))
      return Fail("integer elementary charge overflows int64",error);
    value.charge=next;
    auto& populationCharge=term.population==Population::Tracked ? value.trackedCharge:value.bathCharge;
    if (!Add(populationCharge,charge,&next))
      return Fail("integer population charge overflows int64",error);
    populationCharge=next;
    auto& mass=term.population==Population::Tracked ? value.trackedMass:value.bathMass;
    mass+=static_cast<long double>(term.multiplicity)*species.modelMassKg;
    if (!std::isfinite(mass)) return Fail("reaction model mass is unrepresentable",error);
  }
  *material=value;return true;
}
bool TotalHazard(const std::vector<Hazard>& hazards,long double* result,
                 std::string* error) {
  std::set<std::size_t> seen;long double total=0;
  for (const auto& hazard:hazards) {
    if (!seen.insert(hazard.reaction).second)
      return Fail("duplicate competing-hazard reaction index",error);
    if (!Finite(hazard.rateS) || hazard.rateS<0)
      return Fail("hazard must be finite nonnegative [s^-1]",error);
    total+=hazard.rateS;
  }
  if (!std::isfinite(total) || total>std::numeric_limits<double>::max())
    return Fail("combined hazard is not representable [s^-1]",error);
  *result=total;return true;
}
bool Open(double value) { return Finite(value) && value>0 && value<1; }
bool BoundValid(const ClosureBound& bound) {
  return Finite(bound.absolute) && bound.absolute>=0 &&
      Finite(bound.relative) && bound.relative>=0;
}
bool Close(long double residual,long double scale,const ClosureBound& bound) {
  const long double limit=bound.absolute+static_cast<long double>(bound.relative)*scale;
  return std::isfinite(residual) && std::isfinite(limit) && std::abs(residual)<=limit;
}
bool Convert(long double source,double* output) {
  if (!std::isfinite(source) || std::abs(source)>std::numeric_limits<double>::max()) return false;
  const double value=static_cast<double>(source);
  if (source!=0 && value==0) return false; // Do not silently lose an inventory.
  *output=value;return true;
}
bool Moments(const SpeciesIdentity& species,const Vector3& velocity,double internal,
             double maxSpeed,Vector3* momentum,long double* energy,std::string* error) {
  if (!VectorFinite(velocity) || !Finite(internal))
    return Fail("nonfinite caller-supplied velocity/internal energy",error);
  const double speed=std::hypot(std::hypot(velocity[0],velocity[1]),velocity[2]);
  if (!Finite(speed) || speed>maxSpeed)
    return Fail("kinematics exceeds explicitly qualified classical speed domain [m/s]",error);
  for (int k=0;k<3;++k)
    if (!Convert(static_cast<long double>(species.modelMassKg)*velocity[k],&(*momentum)[k]))
      return Fail("classical momentum is unrepresentable [kg m/s]",error);
  *energy=0.5L*species.modelMassKg*speed*speed+internal;
  if (!std::isfinite(*energy)) return Fail("classical/internal energy is unrepresentable [J]",error);
  return true;
}
}

bool ValidateCatalog(const std::vector<SpeciesIdentity>& catalog,std::string* error) {
  if (catalog.empty() || catalog.size()>static_cast<std::size_t>(std::numeric_limits<int>::max()))
    return Fail("species catalog must have a representable positive size",error);
  std::set<std::pair<std::string,std::string>> states;
  for (const auto& species:catalog) {
    if (species.identity.empty() || species.state.empty() ||
        !states.emplace(species.identity,species.state).second)
      return Fail("species/state identities must be explicit and unique",error);
    if (species.compiledSpecies< -1 || !Positive(species.modelMassKg))
      return Fail("invalid compiled species or nonpositive model mass [kg]",error);
    if (species.nuclei.empty() && species.elementaryCharge==0)
      return Fail("massive neutral species lacks explicit nuclear composition",error);
    std::set<std::string> nuclei;
    for (const auto& nucleus:species.nuclei)
      if (nucleus.identity.empty() || nucleus.count<=0 || !nuclei.insert(nucleus.identity).second)
        return Fail("nuclear identities/counts must be explicit, unique and positive",error);
  }
  if (error) error->clear();
  return true;
}

bool ValidateReaction(const std::vector<SpeciesIdentity>& catalog,
                     const Reaction& reaction,std::string* error) {
  if (!ValidateCatalog(catalog,error)) return false;
  if (reaction.identity.empty() || reaction.reactants.empty())
    return Fail("reaction identity/reactants are missing",error);
  std::size_t trackedParents=0;
  for (const auto& term:reaction.reactants) if (term.population==Population::Tracked) {
    ++trackedParents;
    if (term.multiplicity!=1) return Fail("particle reaction requires one unit tracked parent",error);
  }
  if (trackedParents!=1) return Fail("particle reaction requires exactly one tracked parent term",error);
  Material incoming,outgoing;
  if (!Accumulate(catalog,reaction.reactants,&incoming,error) ||
      !Accumulate(catalog,reaction.products,&outgoing,error)) return false;
  if (incoming.nuclei!=outgoing.nuclei || incoming.charge!=outgoing.charge)
    return Fail("reaction does not conserve exact nuclei/elementary charge",error);
  std::vector<ChemicalSolver::cStoichiometricTerm> terms;
  for (const auto& term:reaction.reactants) terms.push_back({static_cast<int>(term.species),static_cast<double>(term.multiplicity)});
  for (const auto& term:reaction.products) terms.push_back({static_cast<int>(term.species),static_cast<double>(term.multiplicity)});
  // Structural view only. Zero A/excess energy does NOT stand for the actual
  // hazard or physical Q. Signed energy/bath policy is handled separately.
  ChemicalSolver::cReaction entry{reaction.identity.c_str(),0,reaction.reactants.size(),
      reaction.reactants.size(),reaction.products.size(),ChemicalSolver::RateLaw::Constant,0,0,0,1,0};
  ChemicalSolver::cReactionNetwork network(terms.data(),terms.size(),&entry,1,static_cast<int>(catalog.size()));
  if (!network.Validate(error)) return false;
  if (error) error->clear();
  return true;
}

bool TraceInterval::Create(double start,double duration,double elapsed,
                          TraceInterval* output,std::string* error) {
  if (!output || !Finite(start) || !Finite(duration) || !Finite(elapsed) ||
      duration<0 || elapsed<0 || elapsed>duration ||
      !std::isfinite(static_cast<long double>(start)+duration) ||
      std::abs(static_cast<long double>(start)+duration)>std::numeric_limits<double>::max())
    return Fail("invalid physical start/duration/elapsed interval [s]",error);
  TraceInterval value;value.valid_=true;value.start_=start;value.duration_=duration;value.elapsed_=elapsed;
  *output=value;if (error) error->clear();
  return true;
}
bool TraceInterval::valid() const {
  return valid_ && Finite(start_) && Finite(duration_) && Finite(elapsed_) &&
      duration_>=0 && elapsed_>=0 && elapsed_<=duration_;
}
bool TraceInterval::Advance(double seconds,TraceInterval* output,std::string* error) const {
  if (!valid() || !output || !Finite(seconds) || seconds<0 || seconds>remainingS())
    return Fail("local-time advance is outside the retained physical interval [s]",error);
  double next=seconds==remainingS() ? duration_:std::min(duration_,elapsed_+seconds);
  if (seconds>0 && next==elapsed_)
    return Fail("positive local-time advance is not representable [s]",error);
  return Create(start_,duration_,next,output,error);
}

bool OpticalDepthClock::FromOpenUniform(double uniform,OpticalDepthClock* output,
                                      std::string* error) {
  if (!output || !Open(uniform)) return Fail("optical-depth draw requires strictly 0<u<1",error);
  OpticalDepthClock value;value.remaining_=-std::log(uniform);
  if (!Positive(value.remaining_)) return Fail("optical depth is not finite positive",error);
  value.valid_=true;value.fresh_=true;
  *output=value;if (error) error->clear();
  return true;
}
bool OpticalDepthClock::valid() const {
  if (!valid_ || !Finite(remaining_) || !Finite(compensation_)) return false;
  if (consumed_) return remaining_==0 && compensation_==0 && !fresh_;
  const long double pending=static_cast<long double>(remaining_)-compensation_;
  return pending>0 && std::isfinite(pending);
}
void OpticalDepthClock::SubtractExposure(long double exposure) {
  // Kahan-style debit with an explicit rounding residual. Even if adding a
  // very small exposure cannot change remaining_, compensation_ carries it
  // into later segments instead of resetting that probability to zero.
  const long double debit=-exposure-compensation_;
  const double next=static_cast<double>(static_cast<long double>(remaining_)+debit);
  compensation_=static_cast<double>((static_cast<long double>(next)-remaining_)-debit);
  remaining_=next;
}
bool ConsumeExposure(const OpticalDepthClock& clock,double exposure,
                     OpticalDepthClock* output,std::string* error) {
  if (!output || !clock.valid() || clock.consumed_ || !Finite(exposure) || exposure<0)
    return Fail("invalid pending clock/integrated optical exposure",error);
  if (exposure>=static_cast<long double>(clock.remaining_)-clock.compensation_)
    return Fail("integrated exposure reaches an event: locate it on the actual trajectory first",error);
  OpticalDepthClock value=clock;value.SubtractExposure(exposure);value.fresh_=false;
  *output=value;if (error) error->clear();
  return true;
}
bool EventProbability(const std::vector<Hazard>& hazards,double dt,double* output,
                      std::string* error) {
  if (!output || !Finite(dt) || dt<0) return Fail("invalid event-probability duration [s]",error);
  long double total;if (!TotalHazard(hazards,&total,error)) return false;
  const long double exposure=total*dt;
  const long double probability=exposure==0 ? 0:-std::expm1(-exposure);
  if (!std::isfinite(probability) || probability<0 || probability>1)
    return Fail("invalid combined event probability",error);
  *output=static_cast<double>(probability);if (error) error->clear();
  return true;
}
bool ConstantHazardEventDelay(const std::vector<Hazard>& hazards,
    const OpticalDepthClock& clock,long double* delayS,std::string* error) {
  if (!delayS || !clock.valid() || clock.consumed_)
    return Fail("constant-hazard preview requires a valid unconsumed optical clock",error);
  long double total;if (!TotalHazard(hazards,&total,error))return false;
  if (total==0) {
    *delayS=std::numeric_limits<long double>::infinity();
    if(error)error->clear();
    return true;
  }
  const long double pending=static_cast<long double>(clock.remaining_)-clock.compensation_;
  long double exact=pending/total;
  if (!(exact>0) || !std::isfinite(exact))
    return Fail("constant-hazard event delay is not representable as positive finite seconds",error);
  // Directed representational rounding, never an optical-depth tolerance.
  if(total*exact<pending)
    exact=std::nextafter(exact,std::numeric_limits<long double>::infinity());
  *delayS=exact;if(error)error->clear();return true;
}
bool ConstantHazardEventDelay(const std::vector<Hazard>& hazards,
    const OpticalDepthClock& clock,double* delayS,std::string* error) {
  if(!delayS)return Fail("constant-hazard preview requires output storage",error);
  long double exact=0;
  if(!ConstantHazardEventDelay(hazards,clock,&exact,error))return false;
  if(std::isinf(exact)){*delayS=std::numeric_limits<double>::infinity();return true;}
  if(exact>std::numeric_limits<double>::max())
    return Fail("constant-hazard event delay is not representable as positive finite seconds",error);
  double represented=static_cast<double>(exact);
  if (!(represented>0) || !std::isfinite(represented))
    return Fail("positive constant-hazard event delay underflows physical-time storage",error);
  if (static_cast<long double>(represented)<exact)
    represented=std::nextafter(represented,std::numeric_limits<double>::infinity());
  if(!std::isfinite(represented))
    return Fail("upward-rounded constant-hazard event delay overflows seconds",error);
  *delayS=represented;if(error)error->clear();return true;
}
bool SelectHazardChannel(const std::vector<Hazard>& hazards,double uniform,
    std::size_t* reaction,std::string* error) {
  if(!reaction || !Open(uniform))
    return Fail("channel selection requires output storage and strictly 0<u<1",error);
  long double total;if(!TotalHazard(hazards,&total,error))return false;
  if(!(total>0))return Fail("zero instantaneous total rate has no conditional event channel",error);
  const long double target=static_cast<long double>(uniform)*total;
  long double cumulative=0;
  for(const auto& hazard:hazards) {
    cumulative+=hazard.rateS;
    if(hazard.rateS>0 && target<cumulative) {
      *reaction=hazard.reaction;if(error)error->clear();return true;
    }
  }
  return Fail("competing-hazard CDF failed to select a positive channel",error);
}

bool OpticalDepthClock::ExactPendingDepth(long double* output,std::string* error) const {
  if(!output || !valid() || consumed_)
    return Fail("exact exposure preview requires output and a pending optical clock",error);
  // Error-free TwoSum of the two exactly representable binary64 components.
  // A nonzero residual proves that a scalar long-double projection loses bits.
  const long double high=remaining_,low=-static_cast<long double>(compensation_);
  const long double sum=high+low,part=sum-high;
  const long double residual=(high-(sum-part))+(low-part);
  if(residual!=0)
    return Fail("compensated pending optical depth is not exactly representable by this long-double exposure interface; retain origin state and qualify a multi-component depth interface",error);
  if(!(sum>0) || !std::isfinite(sum))return Fail("invalid exact pending optical depth",error);
  *output=sum;if(error)error->clear();return true;
}

bool OpticalDepthClock::PendingDepthComponents(double* high,double* low,std::string* error)const {
  if(!high||!low||!valid()||consumed_)
    return Fail("pending components require output and a valid unconsumed optical clock",error);
  *high=remaining_;*low=-compensation_;if(error)error->clear();return true;
}
bool OpticalDepthClock::SamePendingDepth(const OpticalDepthClock& other)const {
  return valid()&&other.valid()&&consumed_==other.consumed_&&
      remaining_==other.remaining_&&compensation_==other.compensation_;
}

bool ApplyLocatedHazardClock(const OpticalDepthClock& clock,
    const PiecewiseHazard::ClockUpdate& location,const std::vector<Hazard>& rates,
    double uniform,ExposureClockAdvance* output,std::string* error) {
  if(!output||!clock.valid()||clock.consumed_||!location.valid()||
      location.origin().high!=clock.remaining_||location.origin().low!=-clock.compensation_||
      !std::isfinite(location.residualRoundingErrorBound())||location.residualRoundingErrorBound()<0)
    return Fail("located clock publication requires the matching full origin components and an exact locator receipt",error);
  ExposureClockAdvance value;value.clock=clock;value.clock.fresh_=false;
  value.residualDepthRoundingErrorBound=location.residualRoundingErrorBound();
  if(location.event()) {
    if(!SelectHazardChannel(rates,uniform,&value.reaction,error))return false;
    long double total;if(!TotalHazard(rates,&total,error))return false;
    value.totalInstantaneousRateS=static_cast<double>(total);
    if(!std::isfinite(value.totalInstantaneousRateS)||!(value.totalInstantaneousRateS>0))
      return Fail("instantaneous total is outside native event-rate storage [s^-1]",error);
    value.event=true;value.clock.remaining_=0;value.clock.compensation_=0;value.clock.consumed_=true;
  }else {
    // The locator formed the exact rational difference BEFORE rounding its
    // high/low storage components. Never debit a rounded scalar integral or
    // reclassify a near-boundary no-event decision as an event here.
    value.clock.remaining_=location.residual().high;
    value.clock.compensation_=-location.residual().low;
    if(!value.clock.valid())return Fail("located positive residual does not fit the native optical-clock representation",error);
  }
  *output=std::move(value);if(error)error->clear();return true;
}

bool AdvanceIntegratedHazardClock(long double exposure,
    const OpticalDepthClock& clock,const std::vector<Hazard>& rates,
    double uniform,ExposureClockAdvance* output,std::string* error) {
  if(!output || !std::isfinite(exposure) || exposure<0 || !clock.valid() || clock.consumed_)
    return Fail("integrated exposure needs a finite nonnegative depth and pending origin clock",error);
  long double pending;
  if(!clock.ExactPendingDepth(&pending,error))return false;
  ExposureClockAdvance value;value.clock=clock;value.clock.fresh_=false;
  if(exposure<pending) {
    value.clock.SubtractExposure(exposure);
    if(!(value.clock.remainingExtended()>0))
      return Fail("positive residual integrated optical depth is unrepresentable",error);
  }
  else {
    if(!SelectHazardChannel(rates,uniform,&value.reaction,error))return false;
    long double total;if(!TotalHazard(rates,&total,error))return false;
    value.totalInstantaneousRateS=static_cast<double>(total);
    value.event=true;value.clock.remaining_=0;value.clock.compensation_=0;
    value.clock.consumed_=true;
  }
  *output=value;if(error)error->clear();return true;
}

bool AdvanceConstantHazardClock(const std::vector<Hazard>& hazards,
    long double elapsed,const OpticalDepthClock& clock,double channelUniform,
    HazardClockAdvance* output,std::string* error) {
  if(!output || !std::isfinite(elapsed) || elapsed<0 || !clock.valid() || clock.consumed_)
    return Fail("constant-clock reconstruction requires nonnegative finite total seconds and a pending origin clock",error);
  long double total;if(!TotalHazard(hazards,&total,error))return false;
  HazardClockAdvance value;value.clock=clock;value.clock.fresh_=false;
  value.totalRateS=static_cast<double>(total);
  const long double exposure=total*elapsed;
  const long double pending=static_cast<long double>(clock.remaining_)-clock.compensation_;
  if(total==0 || exposure<pending) {
    if(total>0) {
      value.clock.SubtractExposure(exposure);
      if(!Positive(value.clock.remaining()))return Fail("residual optical depth is unrepresentable",error);
    }
  }
  else {
    value.eventOffsetS=pending/total;
    if(!SelectHazardChannel(hazards,channelUniform,&value.reaction,error))return false;
    value.event=true;value.clock.remaining_=0;value.clock.compensation_=0;value.clock.consumed_=true;
  }
  *output=value;if(error)error->clear();return true;
}
bool AdvanceConstantHazards(const std::vector<Hazard>& hazards,
    const TraceInterval& interval,const OpticalDepthClock& clock,
    double channelUniform,HazardAdvance* output,std::string* error) {
  if (!output || !interval.valid())
    return Fail("constant-hazard advance requires a valid interval and pending clock",error);
  HazardClockAdvance decision;
  if(!AdvanceConstantHazardClock(hazards,interval.remainingS(),clock,channelUniform,&decision,error))return false;
  HazardAdvance value;value.interval=interval;value.clock=decision.clock;
  value.totalRateS=decision.totalRateS;value.event=decision.event;value.reaction=decision.reaction;
  if(!decision.event) {
    if(!interval.Advance(interval.remainingS(),&value.interval,error))return false;
  }
  else {
    value.eventOffsetS=std::min(interval.remainingS(),static_cast<double>(decision.eventOffsetS));
    if(!(value.eventOffsetS>0) || !Finite(value.eventOffsetS))
      return Fail("positive event offset is not representable [s]",error);
    if(!interval.Advance(value.eventOffsetS,&value.interval,error))return false;
  }
  *output=value;if(error)error->clear();return true;
}

bool MakeEventPlan(const std::vector<SpeciesIdentity>& catalog,const Reaction& reaction,
    const ParentState& parent,const TraceInterval& interval,PlanScope scope,
    const KinematicPolicy* policy,EventPlan* output,std::string* error) {
  if (!output || !ValidateReaction(catalog,reaction,error)) return false;
  if (!interval.valid() || parent.species>=catalog.size() || !Positive(parent.referenceWeight) ||
      !Positive(parent.individualCorrection) || !VectorFinite(parent.positionM) ||
      !VectorFinite(parent.velocityMS) || !Finite(parent.internalEnergyJ))
    return Fail("invalid event parent/physical interval/weight/phase state",error);
  for (const auto& term:reaction.reactants)
    if (term.population==Population::Tracked && term.species!=parent.species)
      return Fail("event parent species/state does not match reaction",error);
  EventPlan value;value.scope_=scope;value.parent_=parent;value.interval_=interval;
  value.reactionIdentity_=reaction.identity;
  value.parentIdentity_=catalog[parent.species].identity;
  value.parentState_=catalog[parent.species].state;
  value.parentCompiledSpecies_=catalog[parent.species].compiledSpecies;
  Material incoming,outgoing;
  if (!Accumulate(catalog,reaction.reactants,&incoming,error) ||
      !Accumulate(catalog,reaction.products,&outgoing,error)) return false;
  auto& ledger=value.ledger_;ledger.integerConservationVerified=true;
  for (const auto& n:incoming.nuclei) ledger.reactantNuclei.push_back({n.first,n.second});
  for (const auto& n:outgoing.nuclei) ledger.productNuclei.push_back({n.first,n.second});
  ledger.reactantChargeElementary=incoming.charge;ledger.productChargeElementary=outgoing.charge;
  if (!Subtract(outgoing.trackedCharge,incoming.trackedCharge,&ledger.trackedChargeChangeElementary) ||
      !Subtract(outgoing.bathCharge,incoming.bathCharge,&ledger.bathChargeChangeElementary))
    return Fail("integer tracked/bath charge difference overflows",error);
  double represented;
  if (!Convert(static_cast<long double>(parent.referenceWeight)*parent.individualCorrection,&represented) || !Positive(represented))
    return Fail("represented parent number is unrepresentable",error);
  ledger.parentNumber=represented;
  if (!Convert(static_cast<long double>(represented)*catalog[parent.species].modelMassKg,&ledger.parentModelMassKg))
    return Fail("represented parent model mass is unrepresentable",error);
  long double declaredProducts=0;
  for (const auto& term:reaction.products)
    if (term.population==Population::Tracked) declaredProducts+=term.multiplicity;
  if (!Convert(declaredProducts*represented,&ledger.productNumber) ||
      !Convert(outgoing.trackedMass*represented,&ledger.productModelMassKg) ||
      !Convert((outgoing.bathMass-incoming.bathMass)*represented,&ledger.bathModelMassChangeKg))
    return Fail("declared represented product inventory is unrepresentable",error);
  for (std::size_t side=0;side<2;++side) {
    const auto& terms=side==0 ? reaction.reactants:reaction.products;
    for (std::size_t i=0;i<terms.size();++i) if (terms[i].population==Population::Bath) {
      double count;
      if (!Convert(static_cast<long double>(represented)*terms[i].multiplicity,&count))
        return Fail("represented bath material number is unrepresentable",error);
      ledger.bathTerms.push_back({terms[i].species,i,side==0,count});
    }
  }
  if (scope==PlanScope::RateOnly) {
    value.valid_=true;*output=std::move(value);if (error) error->clear();
  return true;
  }
  if (scope!=PlanScope::Transport || !policy || !policy->qualified ||
      policy->identity.empty() || policy->evidenceIdentity.empty() || policy->eventFrame.empty() ||
      !Positive(policy->maximumClassicalSpeedMS) || policy->maximumClassicalSpeedMS>=SpeedOfLightMS ||
      !policy->bath.supplied || !Finite(policy->bath.modelMassChangeKg) ||
      !Finite(policy->bath.energyChangeJ) || !VectorFinite(policy->bath.momentumChangeKgMS) ||
      !BoundValid(policy->bounds.modelMassKg) || !BoundValid(policy->bounds.momentumKgMS) ||
      !BoundValid(policy->bounds.energyJ))
    return Fail("transport requires explicit qualified kinematics/frame/classical domain/bath/closure policy",error);
  value.policyIdentity_=policy->identity;value.frame_=policy->eventFrame;
  std::size_t expected=0;
  for (const auto& term:reaction.products) if (term.population==Population::Tracked) {
    if (term.multiplicity>std::numeric_limits<std::size_t>::max()-expected)
      return Fail("tracked product packet size overflows",error);
    expected+=term.multiplicity;
  }
  if (policy->products.size()!=expected) return Fail("missing/extra tracked product occurrences",error);
  std::set<std::pair<std::size_t,std::uint32_t>> occurrences;
  Vector3 beforeMomentum;long double beforeEnergy;
  if (!Moments(catalog[parent.species],parent.velocityMS,parent.internalEnergyJ,
        policy->maximumClassicalSpeedMS,&beforeMomentum,&beforeEnergy,error)) return false;
  long double afterMass=0,afterEnergy=0,number=0;
  long double weightedMass=0,weightedEnergy=0;
  std::array<long double,3> afterMomentum{{0,0,0}};
  std::array<long double,3> weightedMomentum{{0,0,0}};
  for (const auto& product:policy->products) {
    if (product.term>=reaction.products.size()) return Fail("product term index is out of range",error);
    const auto& term=reaction.products[product.term];
    if (term.population!=Population::Tracked || product.occurrence>=term.multiplicity ||
        !occurrences.emplace(product.term,product.occurrence).second)
      return Fail("invalid/duplicate product occurrence/state",error);
    if (!Positive(product.referenceWeight) || !product.nextClock.valid() ||
        product.nextClock.consumed() || !product.nextClock.fresh())
      return Fail("product needs positive reference weight and fresh, unconsumed Markov clock",error);
    ProductPacket packet;packet.species=term.species;packet.term=product.term;
    packet.compiledSpecies=catalog[term.species].compiledSpecies;
    packet.speciesIdentity=catalog[term.species].identity;
    packet.stateIdentity=catalog[term.species].state;
    packet.occurrence=product.occurrence;packet.targetRepresentedNumber=represented;
    packet.referenceWeight=product.referenceWeight;
    if (!Convert(static_cast<long double>(represented)/product.referenceWeight,&packet.individualCorrection) || !Positive(packet.individualCorrection))
      return Fail("product statistical correction is unrepresentable",error);
    // This is the SAME double multiplication used by the particle adapter.
    // Record the actually represented number, not an ideal quotient whose
    // round-trip could differ. Numerical closure below includes that drift.
    packet.representedNumber=packet.referenceWeight*packet.individualCorrection;
    if (!Positive(packet.representedNumber))
      return Fail("product reference-weight round trip is unrepresentable",error);
    packet.positionM=parent.positionM;packet.velocityMS=product.velocityMS;
    packet.internalEnergyJ=product.internalEnergyJ;packet.interval=interval;packet.clock=product.nextClock;
    Vector3 momentum;long double energy;
    if (!Moments(catalog[term.species],product.velocityMS,product.internalEnergyJ,
        policy->maximumClassicalSpeedMS,&momentum,&energy,error)) return false;
    for (int k=0;k<3;++k) {
      afterMomentum[k]+=momentum[k];
      weightedMomentum[k]+=static_cast<long double>(momentum[k])*packet.representedNumber;
    }
    afterMass+=catalog[term.species].modelMassKg;afterEnergy+=energy;
    weightedMass+=static_cast<long double>(catalog[term.species].modelMassKg)*packet.representedNumber;
    weightedEnergy+=energy*packet.representedNumber;number+=packet.representedNumber;
    value.products_.push_back(packet);
  }
  const long double beforeMass=catalog[parent.species].modelMassKg;
  const long double bathMass=policy->bath.modelMassChangeKg;
  const long double expectedBath=outgoing.bathMass-incoming.bathMass;
  if (!Close(bathMass-expectedBath,std::max(std::abs(expectedBath),std::abs(bathMass)),policy->bounds.modelMassKg) ||
      !Close(afterMass-beforeMass+bathMass,std::max({std::abs(beforeMass),std::abs(afterMass),std::abs(bathMass)}),policy->bounds.modelMassKg) ||
      !Close(weightedMass/represented-beforeMass+bathMass,std::max({std::abs(beforeMass),std::abs(weightedMass/represented),std::abs(bathMass)}),policy->bounds.modelMassKg))
    return Fail("declared model-mass/bath transfer does not close",error);
  for (int k=0;k<3;++k) {
    const long double change=afterMomentum[k]-beforeMomentum[k],bath=policy->bath.momentumChangeKgMS[k];
    const long double weightedChange=weightedMomentum[k]-static_cast<long double>(beforeMomentum[k])*represented;
    if (!Close(change+bath,std::max({std::abs(afterMomentum[k]),std::abs(static_cast<long double>(beforeMomentum[k])),std::abs(bath)}),policy->bounds.momentumKgMS) ||
        !Close(weightedChange/represented+bath,std::max({std::abs(weightedMomentum[k]/represented),std::abs(static_cast<long double>(beforeMomentum[k])),std::abs(bath)}),policy->bounds.momentumKgMS) ||
        !Convert(weightedChange,&ledger.trackedMomentumChangeKgMS[k]) ||
        !Convert(bath*represented,&ledger.bathMomentumChangeKgMS[k]))
      return Fail("momentum/bath transfer does not close or is unrepresentable",error);
  }
  const long double energyChange=afterEnergy-beforeEnergy,bathEnergy=policy->bath.energyChangeJ;
  const long double weightedEnergyChange=weightedEnergy-beforeEnergy*represented;
  if (!Close(energyChange+bathEnergy,std::max({std::abs(beforeEnergy),std::abs(afterEnergy),std::abs(bathEnergy)}),policy->bounds.energyJ) ||
      !Close(weightedEnergyChange/represented+bathEnergy,std::max({std::abs(beforeEnergy),std::abs(weightedEnergy/represented),std::abs(bathEnergy)}),policy->bounds.energyJ) ||
      !Convert(weightedEnergyChange,&ledger.trackedEnergyChangeJ) ||
      !Convert(bathEnergy*represented,&ledger.bathEnergyChangeJ) ||
      !Convert(bathMass*represented,&ledger.bathModelMassChangeKg) ||
      !Convert(weightedMass,&ledger.productModelMassKg) ||
      !Convert(number,&ledger.productNumber))
    return Fail("energy/bath/represented product ledger does not close or is unrepresentable",error);
  ledger.kinematicConservationVerified=true;
  value.valid_=true;*output=std::move(value);if (error) error->clear();
  return true;
}

AttemptStatus PendingEvent::TryCommit(TransactionBackend* backend,std::string* error) {
  if (!backend || committed_ || !plan_.valid() || plan_.scope()!=PlanScope::Transport) {
    Fail("only an uncommitted transport plan can reserve/commit",error);return AttemptStatus::Error;
  }
  try {
    std::string localError;
    auto* diagnostic=error ? error:&localError;
    diagnostic->clear();
    const auto reserved=backend->Reserve(plan_,diagnostic);
    if (reserved!=ReservationStatus::Ready) {
      backend->Rollback();
      if (reserved==ReservationStatus::Deferred) return AttemptStatus::Deferred;
      if (diagnostic->empty()) *diagnostic="reaction backend reservation failed";
      return AttemptStatus::Error;
    }
    if (!backend->Stage(plan_,diagnostic)) {
      backend->Rollback();
      if (diagnostic->empty()) *diagnostic="reaction backend staging failed";
      return AttemptStatus::Error;
    }
    if (!backend->Commit(plan_,diagnostic)) {
      backend->Rollback();
      if (diagnostic->empty()) *diagnostic="reaction backend atomic commit failed";
      return AttemptStatus::Error;
    }
    committed_=true;if (error) error->clear();return AttemptStatus::Committed;
  } catch (const std::exception& exception) {
    backend->Rollback();Fail(std::string("reaction transaction backend exception: ")+exception.what(),error);
  } catch (...) {
    backend->Rollback();Fail("reaction transaction backend nonstandard exception",error);
  }
  return AttemptStatus::Error;
}

} } } // namespace Moon::Species::ReactionExecution
