#include "sep_frozen_collision_transaction_validation.h"
#include "sep_frozen_collision_transaction.h"

#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>
#include <memory>
#include <stdexcept>

namespace SEP { namespace Testing { namespace {
namespace F=Transport::FrozenCollisionTransaction;
namespace RC=Transport::ResonantCollision;
namespace PW=Turbulence::Polarized;
const double C=299792458.0;

Result Finish(bool pass,const char* message) {
  Result result;result.status=pass?Status::Pass:Status::Fail;result.message=message;
  result.metrics.push_back({"assertion_failures",pass?0.0:1.0,0.0,"==","count"});
  result.configuration.push_back("evidence_level=prepare-only-frozen-wave-particle-bookkeeping");
  result.configuration.push_back("coupled_CTMC=false;native_PIC_commit=false;production_qualification=false");
  result.configuration.push_back("RNG=none;accepted_proposals_are_never_redrawn_or_suppressed");
  return result;
}

struct Fixture {
  std::shared_ptr<const PW::Snapshot> before;
  std::vector<F::ParticleRecord> particles;
  std::vector<F::WeightedProposal> events;
  F::PhaseContract phase;
  RC::Generator firstGenerator;
};
Fixture Base() {
  // Explicit synthetic SI values cover bookkeeping, not any adopted WKB or
  // warm-plasma domain. Both directions and both helicities remain independent.
  PW::Initialization input;input.grid.centersPerM={1e-8,1e-7,1e-6};
  input.grid.logWidth=std::log(10.0);
  input.plus.known=input.minus.known=true;input.plus.leftFraction=.2;input.minus.leftFraction=.7;
  input.plus.provenance=input.minus.provenance="explicit frozen transaction mathematical fixture";
  input.vacuumPermeabilityHPerM=1.256637062e-6;input.epochS=0;input.validUntilS=1;
  input.geometryGeneration=3;input.backgroundIdentity="frozen-transaction-analytic:g3";
  input.configurationFingerprint="controlled-transaction-values";
  input.initializerProvenance="explicit four-spectrum neutral transaction fixture";
  PW::InitialCell cell;cell.fieldLine=7;cell.segment=11;cell.volumeM3=1e6;cell.meanMagneticFieldT=5e-9;
  cell.massDensityKgM3=cell.meanMagneticFieldT*cell.meanMagneticFieldT/
      (input.vacuumPermeabilityHPerM*4e5*4e5);
  cell.ePlusJ=4e-7;cell.eMinusJ=2e-7;
  cell.domain.authority=RC::DomainAuthority::ValidationStudy;
  cell.domain.provenance="synthetic frozen-wave transaction; no qualified physical limits";
  input.cells.push_back(cell);
  const auto prepared=PW::PrepareInitialization(input);
  if(!prepared.status.ok())throw std::runtime_error(prepared.status.message);
  Fixture fixture;fixture.before=prepared.snapshot;
  auto& phase=fixture.phase;phase.phase=F::Phase::AcceptedFrozenWaveBatch;
  phase.identity="accepted-frozen-operator-batch:0";phase.epochS=0;
  phase.geometryGeneration=fixture.before->geometryGeneration;phase.snapshotRevision=fixture.before->revision;
  phase.snapshotIdentity=fixture.before->identity;phase.backgroundIdentity=fixture.before->backgroundIdentity;
  phase.relativeTolerance=1e-12;phase.operationProvenance="controlled neutral accepted-proposal transaction";
  return fixture;
}
void Append(Fixture& fixture,std::uint64_t identity,RC::Branch branch,double weight=1e6) {
  const auto view=PW::CollisionView(*fixture.before,0);
  RC::Particle point;point.massKg=1.67262192369e-27;point.chargeC=1.602176634e-19;
  const double speed=2e7,mu=.4;
  const double momentum=point.massKg*speed/std::sqrt((1-speed/C)*(1+speed/C));
  point.parallelMomentumKgMPerS=momentum*mu;
  point.perpendicularMomentumKgMPerS=momentum*std::sqrt(1-mu*mu);
  const auto orbit=RC::PrepareOrbit(view.spectrum,branch,point);
  const auto generator=RC::BuildGenerator(view.spectrum,orbit.orbit,32);
  if(!view.status.ok()||!orbit.status.ok()||!generator.status.ok())
    throw std::runtime_error("cannot construct explicit frozen collision fixture");
  const auto from=RC::CellAt(generator.generator,orbit.orbit.initialWavePitch);
  if(from+1>=generator.generator.centers.size())throw std::runtime_error("fixture has no next pitch cell");
  F::WeightedProposal event;event.key.particleIdentity=identity;event.key.eventIndex=0;
  event.fieldLine=7;event.segment=11;event.statisticalWeight=weight;event.branch=branch;
  event.phaseIdentity=fixture.phase.identity;event.expectedLifecycleIdentity="immutable-owner-token:"+std::to_string(identity);
  event.proposal=RC::ProposeJump(generator.generator,generator.generator.centers[from],
      generator.generator.centers[from+1],weight);
  if(!event.proposal.status.ok()||event.proposal.exchange.empty())
    throw std::runtime_error("fixture did not produce an accepted active neighboring jump");
  F::ParticleRecord particle;particle.identity=identity;particle.fieldLine=7;particle.segment=11;
  particle.statisticalWeight=weight;particle.lifecycleIdentity=event.expectedLifecycleIdentity;
  // The literal fixture's true encoded parent is exactly this prepared phase
  // point. A one-ULP mismatched parent is tested as a rejection below; the
  // transaction must never overwrite it with an approximate inverse transform.
  particle.phasePoint=event.proposal.before;
  fixture.particles.push_back(particle);fixture.events.push_back(event);
  if(fixture.events.size()==1)fixture.firstGenerator=generator.generator;
}
long double Kinetic(const RC::Particle& point) {
  // Independent SI reference from the defining relativistic energy. Long
  // double resolves the nonrelativistic kinetic difference in this fixture;
  // production proposal increments instead use the invariant u*Delta p.
  const long double mc=static_cast<long double>(point.massKg)*C;
  return (std::sqrt(mc*mc+static_cast<long double>(point.parallelMomentumKgMPerS)*point.parallelMomentumKgMPerS+
      static_cast<long double>(point.perpendicularMomentumKgMPerS)*point.perpendicularMomentumKgMPerS)-mc)*C;
}
bool RejectsCleanly(const Fixture& fixture,const std::vector<F::WeightedProposal>& events,
                    const F::PhaseContract& phase) {
  const auto result=F::Prepare(*fixture.before,fixture.particles,events,phase);
  return !result.status.ok()&&!result.candidate&&result.particles.empty()&&result.bins.empty()&&
      result.events.empty()&&
      result.operation.beforeIdentity.empty()&&result.operation.afterIdentity.empty();
}

Result WeightedEnergy() {
  auto fixture=Base();Append(fixture,81,RC::Branch::AlongB);Append(fixture,19,RC::Branch::AgainstB,2e6);
  const auto prepared=F::Prepare(*fixture.before,fixture.particles,fixture.events,fixture.phase);
  if(!prepared.status.ok())return Finish(false,prepared.status.message.c_str());
  long double independent=0,scale=0;
  bool pass=prepared.candidate&&!prepared.qualifiedForProduction&&!prepared.candidate->qualifiedForProduction&&
      prepared.candidate->revision==fixture.before->revision+1&&
      prepared.operation.beforeIdentity==fixture.before->identity&&
      prepared.operation.afterIdentity==prepared.candidate->identity;
  for(const auto& event:fixture.events) {
    const long double delta=event.statisticalWeight*(Kinetic(event.proposal.after)-Kinetic(event.proposal.before));
    independent+=delta;scale+=std::fabs(delta);
    const auto particle=std::find_if(prepared.particles.begin(),prepared.particles.end(),[&](const F::ParticleRecord& record){
      return record.identity==event.key.particleIdentity;
    });
    pass=pass&&particle!=prepared.particles.end()&&particle->nextEventIndex==1&&
        particle->lifecycleIdentity==event.expectedLifecycleIdentity&&
        particle->phasePoint.parallelMomentumKgMPerS==event.proposal.after.parallelMomentumKgMPerS;
  }
  const double error=static_cast<double>(std::fabs(independent-prepared.ledger.particleEnergyChangeJ)/scale);
  pass=pass&&error<1e-10&&prepared.ledger.dimensionlessTransferResidual<=fixture.phase.relativeTolerance&&
      prepared.ledger.dimensionlessActualResidual<=fixture.phase.relativeTolerance&&
      PW::ValidateSnapshot(*prepared.candidate).ok()&&PW::ValidateSnapshot(*fixture.before).ok();
  auto result=Finish(pass,"weighted particle increments agree with independent relativistic energy and twice-magnetic wave transfer on both branches");
  result.metrics.push_back({"independent_weighted_energy_relative_error",error,1e-10,"<","dimensionless"});
  return result;
}
Result RoundoffLedger() {
  auto fixture=Base();Append(fixture,81,RC::Branch::AlongB);
  const auto prepared=F::Prepare(*fixture.before,fixture.particles,fixture.events,fixture.phase);
  bool pass=prepared.status.ok()&&prepared.candidate&&prepared.ledger.requestedMagneticChangeJ!=0&&
      prepared.ledger.actualMagneticChangeJ!=0;
  if(prepared.status.ok())for(const auto& bin:prepared.bins) {
    const long double actual=static_cast<long double>(bin.afterMagneticJ)-bin.beforeMagneticJ;
    pass=pass&&bin.magneticRoundoffJ==actual-bin.requestedMagneticChangeJ&&
        std::fabs(bin.magneticRoundoffJ)<=fixture.phase.relativeTolerance*std::fabs(bin.requestedMagneticChangeJ);
  }
  auto unresolved=Base();Append(unresolved,81,RC::Branch::AlongB,1e-20);
  const auto lost=F::Prepare(*unresolved.before,unresolved.particles,unresolved.events,unresolved.phase);
  pass=pass&&lost.status.code==Transport::StatusCode::UnresolvedCoefficient&&
      RejectsCleanly(unresolved,unresolved.events,unresolved.phase)&&unresolved.particles[0].nextEventIndex==0;
  return Finish(pass,"ordinary bounded transfer rounding is explicit; a sub-ULP debit returns Unresolved with no staged particle or wave candidate");
}
Result InsufficientBudget() {
  auto fixture=Base();Append(fixture,19,RC::Branch::AgainstB);
  Append(fixture,81,RC::Branch::AlongB,1e20);
  const auto inputIdentity=fixture.before->identity;
  const bool pass=fixture.events[1].proposal.particleEnergyChangeJ>0&&
      RejectsCleanly(fixture,fixture.events,fixture.phase)&&fixture.before->identity==inputIdentity&&
      fixture.particles[0].nextEventIndex==0&&PW::ValidateSnapshot(*fixture.before).ok();
  return Finish(pass,"an over-budget accepted event rejects the entire batch with no clipped bin, modified particle, new revision or replacement event");
}
Result ProvenanceNegatives() {
  auto fixture=Base();Append(fixture,81,RC::Branch::AlongB);
  std::vector<std::function<void(F::WeightedProposal&)>> mutations={
    [](F::WeightedProposal& e){++e.proposal.generation;},
    [](F::WeightedProposal& e){e.proposal.epochS+=1;},
    [](F::WeightedProposal& e){e.proposal.spectrumIdentity+="stale";},
    [](F::WeightedProposal& e){++e.fieldLine;},
    [](F::WeightedProposal& e){e.proposal.exchange[0].bin=999;},
    [](F::WeightedProposal& e){e.proposal.exchange[0].helicity=static_cast<RC::Helicity>(7);},
    [](F::WeightedProposal& e){e.proposal.exchange[0].branch=RC::Branch::AgainstB;},
    [](F::WeightedProposal& e){e.statisticalWeight*=2;},
    [](F::WeightedProposal& e){e.phaseIdentity+="wrong";},
    [](F::WeightedProposal& e){e.expectedLifecycleIdentity+="wrong";},
    [](F::WeightedProposal& e){e.proposal.before.parallelMomentumKgMPerS=std::nextafter(e.proposal.before.parallelMomentumKgMPerS,1.0);},
    [](F::WeightedProposal& e){e.proposal.after.perpendicularMomentumKgMPerS*=1.01;},
    [](F::WeightedProposal& e){e.proposal.particleEnergyChangeJ*=2;},
    [](F::WeightedProposal& e){e.proposal.exchange[0].magneticEnergyChangeJ*=2;}
  };
  bool pass=true;unsigned rejected=0;
  for(const auto& mutate:mutations) {
    auto events=fixture.events;mutate(events[0]);const bool clean=RejectsCleanly(fixture,events,fixture.phase);
    if(clean)++rejected;
    pass=clean&&pass;
  }
  auto result=Finish(pass,"stale provenance and tampered phase, weight, exact parent, orbit and branch/helicity/bin transfers fail before any candidate escapes");
  result.metrics.push_back({"rejected_negative_controls",static_cast<double>(rejected),14,"==","count"});
  return result;
}
Result EventLifecycle() {
  auto fixture=Base();Append(fixture,81,RC::Branch::AlongB);
  auto second=fixture.events[0];second.key.eventIndex=1;
  const auto& generator=fixture.firstGenerator;
  const auto firstTarget=RC::CellAt(generator,second.proposal.afterWavePitch);
  if(firstTarget+1>=generator.centers.size())return Finish(false,"fixture lacks a second neighboring pitch cell");
  second.proposal=RC::ProposeJump(generator,second.proposal.afterWavePitch,
      generator.centers[firstTarget+1],second.statisticalWeight);
  auto events=fixture.events;events.push_back(second);
  const auto chained=F::Prepare(*fixture.before,fixture.particles,events,fixture.phase);
  bool pass=second.proposal.status.ok()&&chained.status.ok()&&chained.particles[0].nextEventIndex==2;
  auto duplicate=fixture.events;duplicate.push_back(duplicate[0]);
  pass=pass&&RejectsCleanly(fixture,duplicate,fixture.phase);
  events[1].key.eventIndex=2;pass=pass&&RejectsCleanly(fixture,events,fixture.phase);
  fixture.particles[0].nextEventIndex=std::numeric_limits<std::uint64_t>::max();
  fixture.events[0].key.eventIndex=fixture.particles[0].nextEventIndex;
  pass=pass&&RejectsCleanly(fixture,fixture.events,fixture.phase);
  fixture.particles[0].nextEventIndex=0;fixture.events[0].key.eventIndex=0;
  fixture.particles[0].identity=0;fixture.events[0].key.particleIdentity=0;
  pass=pass&&RejectsCleanly(fixture,fixture.events,fixture.phase);
  return Finish(pass,"same-parent proposals retain exact chain state and consecutive events; duplicate, gap and counter-overflow paths reject atomically");
}
Result CanonicalOrdering() {
  auto fixture=Base();for(unsigned i=0;i<20;++i)Append(fixture,1000+7*i,
      i%2?RC::Branch::AgainstB:RC::Branch::AlongB,1e6+i);
  const auto first=F::Prepare(*fixture.before,fixture.particles,fixture.events,fixture.phase);
  std::reverse(fixture.particles.begin(),fixture.particles.end());std::reverse(fixture.events.begin(),fixture.events.end());
  const auto second=F::Prepare(*fixture.before,fixture.particles,fixture.events,fixture.phase);
  bool pass=first.status.ok()&&second.status.ok()&&first.candidate->identity==second.candidate->identity&&
      first.ledger.particleEnergyChangeJ==second.ledger.particleEnergyChangeJ&&
      first.ledger.actualParticleEnergyChangeJ==second.ledger.actualParticleEnergyChangeJ&&
      first.ledger.particleRoundoffJ==second.ledger.particleRoundoffJ&&
      first.ledger.magneticRoundoffJ==second.ledger.magneticRoundoffJ&&first.bins.size()==second.bins.size();
  if(pass)for(std::size_t index=0;index<first.particles.size();++index)
    pass=pass&&first.particles[index].identity==second.particles[index].identity&&
        first.particles[index].nextEventIndex==second.particles[index].nextEventIndex;
  if(pass)for(std::size_t index=0;index<first.bins.size();++index)
    pass=pass&&first.bins[index].afterMagneticJ==second.bins[index].afterMagneticJ&&
        first.bins[index].requestedMagneticChangeJ==second.bins[index].requestedMagneticChangeJ;
  return Finish(pass,"canonical stable-key and physical-bin reductions give identical staged values and identities under reversed particle/event input order");
}
Result PhaseNegatives() {
  auto fixture=Base();Append(fixture,81,RC::Branch::AlongB);
  const std::vector<std::function<void(F::PhaseContract&)>> mutations={
    [](F::PhaseContract& p){p.phase=F::Phase::Unspecified;},
    [](F::PhaseContract& p){p.identity.clear();},
    [](F::PhaseContract& p){p.epochS=1;},
    [](F::PhaseContract& p){++p.geometryGeneration;},
    [](F::PhaseContract& p){++p.snapshotRevision;},
    [](F::PhaseContract& p){p.snapshotIdentity+="stale";},
    [](F::PhaseContract& p){p.backgroundIdentity+="wrong";},
    [](F::PhaseContract& p){p.relativeTolerance=std::numeric_limits<double>::quiet_NaN();},
    [](F::PhaseContract& p){p.relativeTolerance=-1;},
    [](F::PhaseContract& p){p.operationProvenance.clear();}
  };
  bool pass=true;for(const auto& mutate:mutations) {
    auto phase=fixture.phase;mutate(phase);pass=RejectsCleanly(fixture,fixture.events,phase)&&pass;
  }
  return Finish(pass,"undeclared frozen phase, epoch/generation/revision/background mismatch and invalid numerical tolerance fail without an implicit contract");
}
Result PrivateAtomicBundle() {
  auto fixture=Base();Append(fixture,81,RC::Branch::AlongB);
  struct Bundle {std::shared_ptr<const PW::Snapshot> waves;std::vector<F::ParticleRecord> particles;};
  // This fixture-only pointer is deliberately not a new public state owner.
  // It demonstrates why an outer PIC commit must compare all preconditions
  // before one publication point; actual PIC/RNG/list publication is untested.
  auto initial=std::make_shared<Bundle>();initial->waves=fixture.before;initial->particles=fixture.particles;
  std::shared_ptr<const Bundle> current=initial;
  const auto prepared=F::Prepare(*current->waves,current->particles,fixture.events,fixture.phase);
  if(!prepared.status.ok())return Finish(false,prepared.status.message.c_str());
  auto publish=[&](const std::string& expectedIdentity,std::uint64_t expectedEvent){
    if(current->waves->identity!=expectedIdentity||current->particles[0].nextEventIndex!=expectedEvent)return false;
    auto next=std::make_shared<Bundle>();next->waves=prepared.candidate;next->particles=prepared.particles;
    std::shared_ptr<const Bundle> immutable=next;current.swap(immutable);return true;
  };
  bool pass=!publish("stale",0)&&current.get()==initial.get()&&
      !publish(fixture.before->identity,1)&&current.get()==initial.get();
  pass=pass&&publish(fixture.before->identity,0)&&current->waves.get()==prepared.candidate.get()&&
      current->particles[0].nextEventIndex==1&&initial->particles[0].nextEventIndex==0&&
      !publish(fixture.before->identity,0);
  return Finish(pass,"private fixture rejects stale wave/particle preconditions and swaps one complete immutable bundle; no native commit claim");
}
Result EncodedParticleResolution() {
  auto ordinary=Base();Append(ordinary,81,RC::Branch::AlongB);
  const auto prepared=F::Prepare(*ordinary.before,ordinary.particles,ordinary.events,ordinary.phase);
  bool pass=prepared.status.ok()&&prepared.events.size()==1&&
      prepared.ledger.maximumParticleRelativeRoundoff<=ordinary.phase.relativeTolerance;
  if(prepared.status.ok()) {
    const auto& event=prepared.events[0];
    pass=pass&&event.particleRoundoffJ==event.actualEnergyChangeJ-event.requestedEnergyChangeJ&&
        prepared.ledger.actualParticleEnergyChangeJ==event.actualEnergyChangeJ&&
        prepared.ledger.particleRoundoffJ==event.particleRoundoffJ&&
        prepared.ledger.actualTotalResidualJ==prepared.ledger.actualParticleEnergyChangeJ+
                                            2*prepared.ledger.actualMagneticChangeJ;
  }

  auto unresolved=Base();Append(unresolved,81,RC::Branch::AlongB,1e20);
  auto& wire=unresolved.events[0];const auto& generator=unresolved.firstGenerator;
  // Deliberately create a roundoff-scale wire proposal whose phase is
  // indistinguishable from its parent within the RC orbit validator's error
  // bound. Its finite macro weight still resolves the requested wave debit.
  // This tests a bookkeeping negative control, not a stochastic transition or
  // neighbor-generator acceptance claim. No draw or physical input is guessed.
  double pitch=0;bool active=false;
  for(double center:generator.centers) {
    const auto resonance=RC::EvaluateResonance(generator.spectrum,generator.orbit,center);
    if(resonance.status.ok()&&resonance.active&&(!active||std::fabs(center)<std::fabs(pitch))) {
      pitch=center;active=true;
    }
  }
  if(!active)return Finish(false,"encoded-particle negative fixture lacks an active finite-band phase");
  wire.proposal=RC::ProposeJump(generator,pitch,std::nextafter(pitch,1.0),wire.statisticalWeight);
  wire.proposal.after=wire.proposal.before;
  unresolved.particles[0].phasePoint=wire.proposal.before;
  const auto view=PW::CollisionView(*unresolved.before,0);
  const auto canonical=RC::ValidateProposal(view.spectrum,wire.branch,wire.proposal,wire.statisticalWeight);
  bool waveResolved=wire.proposal.status.ok()&&!wire.proposal.exchange.empty();
  for(const auto& exchange:wire.proposal.exchange) {
    const unsigned spectrum=2*static_cast<unsigned>(exchange.branch)+static_cast<unsigned>(exchange.helicity);
    const double before=unresolved.before->cells[0].magneticEnergyJ[spectrum][exchange.bin];
    const long double requested=exchange.magneticEnergyChangeJ;
    const double after=static_cast<double>(static_cast<long double>(before)+requested);
    const long double actual=static_cast<long double>(after)-before;
    waveResolved=waveResolved&&requested!=0&&after>=0&&actual!=0&&
        std::fabs(actual-requested)<=unresolved.phase.relativeTolerance*std::fabs(requested);
  }
  const auto rejected=F::Prepare(*unresolved.before,unresolved.particles,unresolved.events,unresolved.phase);
  pass=pass&&canonical.ok()&&wire.proposal.particleEnergyChangeJ!=0&&waveResolved&&
      Kinetic(wire.proposal.after)==Kinetic(wire.proposal.before)&&
      rejected.status.code==Transport::StatusCode::UnresolvedCoefficient&&
      rejected.status.message.find("encoded particle transfer")!=std::string::npos&&
      RejectsCleanly(unresolved,unresolved.events,unresolved.phase)&&unresolved.particles[0].nextEventIndex==0;
  auto result=Finish(pass,"actual encoded momentum energy has a signed bounded rounding ledger; a canonical roundoff-scale proposal with zero represented particle transfer rejects unresolved before wave debit");
  result.metrics.push_back({"canonical_negative_wire_accepted",canonical.ok()?1.0:0.0,1.0,"==","count"});
  result.metrics.push_back({"negative_control_wave_debit_resolved",waveResolved?1.0:0.0,1.0,"==","count"});
  return result;
}
} // namespace

std::vector<Descriptor> FrozenCollisionTransactionDescriptors() {
  const char* ids[]={"FCT01","FCT02","FCT03","FCT04","FCT05","FCT06","FCT07","FCT08","FCT09"};
  const char* descriptions[]={"Weighted magnetic/total energy bookkeeping","Explicit magnetic storage roundoff",
      "Whole-batch insufficient budget rejection","Accepted-proposal provenance negative controls",
      "Stable event and exact parent lifecycle","Canonical deterministic batch reduction",
      "Explicit frozen-phase numerical contract","Private atomic bundle publication fixture","Encoded particle energy transfer resolution"};
  const std::function<Result()> callbacks[]={WeightedEnergy,RoundoffLedger,InsufficientBudget,ProvenanceNegatives,
      EventLifecycle,CanonicalOrdering,PhaseNegatives,PrivateAtomicBundle,EncodedParticleResolution};
  std::vector<Descriptor> descriptors;
  for(unsigned index=0;index<9;++index) {
    Descriptor descriptor;descriptor.id=ids[index];descriptor.group="frozen-collision-transaction";
    descriptor.name=descriptions[index];descriptor.description=descriptions[index];descriptor.callback=callbacks[index];
    descriptor.initialization=InitializationLevel::None;descriptor.runtime=RuntimeClass::Routine;
    descriptor.seedCapability=SeedCapability::Deterministic;descriptor.repeatSupported=true;
    descriptor.seedPolicy="deterministic explicit accepted proposals; no RNG access or event redraw";
    descriptor.stateIsolation="local immutable four-state snapshot and private particle records; prepare-only; no PIC/global mutation";
    descriptor.supportedBuildModes="source C++11 and native linked diagnostic; no production qualification";
    descriptors.push_back(descriptor);
  }
  return descriptors;
}
} } // namespace SEP::Testing
