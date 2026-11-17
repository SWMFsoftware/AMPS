#include "sep_frozen_collision_transaction.h"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <limits>
#include <map>
#include <sstream>
#include <tuple>
#include <utility>

namespace SEP { namespace Transport { namespace FrozenCollisionTransaction {
namespace {
namespace PW=Turbulence::Polarized;
namespace RC=ResonantCollision;
using CellKey=std::pair<std::uint64_t,std::uint64_t>;
using BinKey=std::tuple<std::uint64_t,std::uint64_t,unsigned,unsigned,std::size_t>;

Status Error(const std::string& text,StatusCode code) {
  return Status::Error(code,"frozen collision transaction: "+text);
}
Prepared Rejected(const std::string& text,StatusCode code=StatusCode::InvalidArgument) {
  Prepared result;result.status=Error(text,code);return result;
}
bool Same(double a,double b) {
  // Exact encoded state ownership matters here. A loose relative comparison
  // could silently replace the true parent momentum by a rounded inverse
  // Lorentz transform. The producer must bind its proposal to that exact
  // state; an adapter needing a different representation must say so openly.
  return std::memcmp(&a,&b,sizeof(double))==0;
}
bool SameParticle(const RC::Particle& a,const RC::Particle& b) {
  return Same(a.massKg,b.massKg)&&Same(a.chargeC,b.chargeC)&&
      Same(a.parallelMomentumKgMPerS,b.parallelMomentumKgMPerS)&&
      Same(a.perpendicularMomentumKgMPerS,b.perpendicularMomentumKgMPerS);
}
bool ValidParticle(const ParticleRecord& particle) {
  const auto& point=particle.phasePoint;
  return particle.identity!=0&&std::isfinite(particle.statisticalWeight)&&particle.statisticalWeight>0&&
      std::isfinite(point.massKg)&&point.massKg>0&&std::isfinite(point.chargeC)&&point.chargeC!=0&&
      std::isfinite(point.parallelMomentumKgMPerS)&&
      std::isfinite(point.perpendicularMomentumKgMPerS)&&point.perpendicularMomentumKgMPerS>=0&&
      !particle.lifecycleIdentity.empty();
}

// Neumaier reduction in canonical event/bin order limits cancellation loss
// while remaining independent of input traversal, thread order or MPI rank.
// This does not claim exact real arithmetic. Its stored-rounding difference is
// accounted below, rather than hidden behind a positive energy floor.
struct Sum {
  long double main=0.0L,correction=0.0L;
  void Add(long double value) {
    const long double next=main+value;
    correction+=std::fabs(main)>=std::fabs(value)?(main-next)+value:(value-next)+main;
    main=next;
  }
  long double Value()const{return main+correction;}
};
double Dimensionless(long double residual,long double scale) {
  if(!std::isfinite(residual)||!std::isfinite(scale)||scale<0)
    return std::numeric_limits<double>::infinity();
  if(scale==0)return residual==0?0.0:std::numeric_limits<double>::infinity();
  return static_cast<double>(std::fabs(residual)/scale);
}
long double EncodedParticleEnergyChange(const RC::Particle& before,const RC::Particle& after,
                                      double weight) {
  // Evaluate the energy of the actual two encoded phase points, without
  // subtracting rest-mass-dominated total energies or changing either point.
  // E=c*sqrt(m^2*c^2+p_parallel^2+p_perp^2), hence
  // Delta E=c*Delta(p^2)/(sqrt(m^2*c^2+p_after^2)+sqrt(m^2*c^2+p_before^2)).
  // The factored square differences retain tiny represented momentum changes.
  // c is the exact SI speed of light, also used by the canonical RC kernel.
  const long double lightSpeed=299792458.0L,mc=static_cast<long double>(before.massKg)*lightSpeed;
  const long double bp=before.parallelMomentumKgMPerS,bn=before.perpendicularMomentumKgMPerS;
  const long double ap=after.parallelMomentumKgMPerS,an=after.perpendicularMomentumKgMPerS;
  const long double squaredChange=(ap-bp)*(ap+bp)+(an-bn)*(an+bn);
  const long double denominator=std::hypot(mc,std::hypot(ap,an))+std::hypot(mc,std::hypot(bp,bn));
  return static_cast<long double>(weight)*lightSpeed*(squaredChange/denominator);
}
std::string BinContext(const BinKey& key) {
  std::ostringstream out;
  out<<"line="<<std::get<0>(key)<<" segment="<<std::get<1>(key)
      <<" branch="<<std::get<2>(key)<<" helicity="<<std::get<3>(key)
      <<" bin="<<std::get<4>(key);
  return out.str();
}
} // namespace

Prepared Prepare(const PW::Snapshot& before,const std::vector<ParticleRecord>& particles,
                 const std::vector<WeightedProposal>& proposals,const PhaseContract& phase) {
  const auto valid=PW::ValidateSnapshot(before);
  if(!valid.ok())return Rejected("invalid immutable input snapshot: "+valid.message);
  if(phase.phase!=Phase::AcceptedFrozenWaveBatch||phase.identity.empty()||
      phase.operationProvenance.empty()||!std::isfinite(phase.epochS)||
      !std::isfinite(phase.relativeTolerance)||phase.relativeTolerance<0)
    return Rejected("explicit frozen-batch phase, provenance and numerical tolerance are required");
  if(!Same(phase.epochS,before.epochS)||phase.geometryGeneration!=before.geometryGeneration||
      phase.snapshotRevision!=before.revision||phase.snapshotIdentity!=before.identity||
      phase.backgroundIdentity!=before.backgroundIdentity)
    return Rejected("phase epoch/generation/revision/snapshot/background provenance is stale");
  if(proposals.empty())return Rejected("batch contains no accepted proposals");

  // Everything below is private scratch until the final success assignment.
  // Exceptions or failed validators cannot expose a partially advanced parent,
  // a partially depleted bin or a candidate revision to the caller.
  std::vector<ParticleRecord> staged=particles;
  std::sort(staged.begin(),staged.end(),[](const ParticleRecord& a,const ParticleRecord& b){
    return a.identity<b.identity;
  });
  std::map<std::uint64_t,std::size_t> particleIndex;
  for(std::size_t index=0;index<staged.size();++index) {
    if(!ValidParticle(staged[index]))return Rejected("invalid SI particle state, weight or lifecycle identity");
    if(!particleIndex.emplace(staged[index].identity,index).second)
      return Rejected("duplicate particle identity in frozen input");
  }
  std::map<CellKey,std::size_t> cellIndex;
  for(std::size_t index=0;index<before.cells.size();++index)
    cellIndex.emplace(CellKey(before.cells[index].fieldLine,before.cells[index].segment),index);
  for(const auto& particle:staged)
    if(cellIndex.find(CellKey(particle.fieldLine,particle.segment))==cellIndex.end())
      return Rejected("frozen particle record references an absent physical cell");

  std::vector<const WeightedProposal*> events;
  for(const auto& proposal:proposals)events.push_back(&proposal);
  std::sort(events.begin(),events.end(),[](const WeightedProposal* a,const WeightedProposal* b){
    return std::tie(a->key.particleIdentity,a->key.eventIndex)<
        std::tie(b->key.particleIdentity,b->key.eventIndex);
  });
  std::map<BinKey,Sum> requested;
  Sum particleEnergy,actualParticleEnergy,particleRoundoff,absoluteTransferScale;
  std::vector<ParticleLedger> eventLedgers;
  double maximumParticleRelativeRoundoff=0;
  bool havePrevious=false;EventKey previous;
  for(const auto* entry:events) {
    if(havePrevious&&entry->key.particleIdentity==previous.particleIdentity&&entry->key.eventIndex==previous.eventIndex)
      return Rejected("duplicate accepted particle/event key");
    previous=entry->key;havePrevious=true;
    const auto foundParticle=particleIndex.find(entry->key.particleIdentity);
    const auto foundCell=cellIndex.find(CellKey(entry->fieldLine,entry->segment));
    if(foundParticle==particleIndex.end()||foundCell==cellIndex.end())
      return Rejected("accepted proposal references an absent particle or physical cell");
    auto& particle=staged[foundParticle->second];
    if(entry->key.eventIndex!=particle.nextEventIndex||
        particle.nextEventIndex==std::numeric_limits<std::uint64_t>::max())
      return Rejected("accepted event index is stale, nonconsecutive or cannot be advanced without overflow");
    if(entry->fieldLine!=particle.fieldLine||entry->segment!=particle.segment||
        entry->phaseIdentity!=phase.identity||entry->expectedLifecycleIdentity!=particle.lifecycleIdentity||
        !Same(entry->statisticalWeight,particle.statisticalWeight)||
        !SameParticle(entry->proposal.before,particle.phasePoint))
      return Rejected("accepted proposal does not match exact particle/cell/phase/weight/lifecycle state");
    if(entry->branch!=RC::Branch::AlongB&&entry->branch!=RC::Branch::AgainstB)
      return Rejected("accepted proposal has an invalid propagation branch");
    const auto view=PW::CollisionView(before,foundCell->second);
    if(!view.status.ok())return Rejected("cannot obtain canonical frozen collision view: "+view.status.message);
    // The canonical collision kernel owns orbit and helicity/bin partition
    // validation. No second resonance law or approximate redistribution is
    // introduced in the bookkeeping layer, and no random stream is consumed.
    const auto proposalValid=RC::ValidateProposal(view.spectrum,entry->branch,
                                                entry->proposal,entry->statisticalWeight);
    if(!proposalValid.ok())return Rejected("invalid accepted collision proposal: "+proposalValid.message,proposalValid.code);
    ParticleLedger eventLedger;eventLedger.key=entry->key;eventLedger.statisticalWeight=entry->statisticalWeight;
    eventLedger.requestedEnergyChangeJ=entry->proposal.particleEnergyChangeJ;
    eventLedger.actualEnergyChangeJ=EncodedParticleEnergyChange(entry->proposal.before,entry->proposal.after,
                                                             entry->statisticalWeight);
    eventLedger.particleRoundoffJ=eventLedger.actualEnergyChangeJ-eventLedger.requestedEnergyChangeJ;
    eventLedger.dimensionlessRoundoff=Dimensionless(eventLedger.particleRoundoffJ,
                                                   std::fabs(eventLedger.requestedEnergyChangeJ));
    const bool unresolvedSign=eventLedger.requestedEnergyChangeJ!=0&&
        (eventLedger.actualEnergyChangeJ==0||
         (eventLedger.actualEnergyChangeJ>0)!=(eventLedger.requestedEnergyChangeJ>0));
    if(unresolvedSign||!std::isfinite(eventLedger.dimensionlessRoundoff)||
        eventLedger.dimensionlessRoundoff>phase.relativeTolerance) {
      std::ostringstream context;context<<"encoded particle transfer is unresolved at the explicit numerical tolerance: particle="
          <<entry->key.particleIdentity<<" event="<<entry->key.eventIndex;
      return Rejected(context.str(),StatusCode::UnresolvedCoefficient);
    }
    maximumParticleRelativeRoundoff=std::max(maximumParticleRelativeRoundoff,eventLedger.dimensionlessRoundoff);
    actualParticleEnergy.Add(eventLedger.actualEnergyChangeJ);particleRoundoff.Add(eventLedger.particleRoundoffJ);
    eventLedgers.push_back(eventLedger);
    particleEnergy.Add(entry->proposal.particleEnergyChangeJ);
    absoluteTransferScale.Add(std::fabs(static_cast<long double>(entry->proposal.particleEnergyChangeJ)));
    for(const auto& exchange:entry->proposal.exchange) {
      const unsigned branch=static_cast<unsigned>(exchange.branch),helicity=static_cast<unsigned>(exchange.helicity);
      // Bounds are checked again at this ownership boundary, independently of
      // the canonical physics validator, before any array address is formed.
      if(branch>1||helicity>1||exchange.bin>=before.grid.centersPerM.size())
        return Rejected("accepted exchange has an invalid branch/helicity/bin index");
      requested[BinKey(entry->fieldLine,entry->segment,branch,helicity,exchange.bin)].Add(exchange.magneticEnergyChangeJ);
      absoluteTransferScale.Add(2*std::fabs(static_cast<long double>(exchange.magneticEnergyChangeJ)));
    }
    particle.phasePoint=entry->proposal.after;
    ++particle.nextEventIndex;
  }

  auto updatedCells=before.cells;
  std::vector<BinLedger> bins;
  Sum requestedMagnetic,actualMagnetic,roundoff;
  for(const auto& item:requested) {
    const auto& key=item.first;
    const auto foundCell=cellIndex.find(CellKey(std::get<0>(key),std::get<1>(key)));
    auto& cell=updatedCells[foundCell->second];
    auto& stored=cell.magneticEnergyJ[2*std::get<2>(key)+std::get<3>(key)][std::get<4>(key)];
    const long double delta=item.second.Value(),candidate=static_cast<long double>(stored)+delta;
    if(!std::isfinite(delta)||!std::isfinite(candidate))
      return Rejected("nonfinite magnetic wave budget at "+BinContext(key));
    if(candidate<0)
      return Rejected("insufficient magnetic wave budget at "+BinContext(key),StatusCode::UnresolvedCoefficient);
    if(candidate>std::numeric_limits<double>::max())
      return Rejected("magnetic wave budget cannot be represented as double at "+BinContext(key),StatusCode::UnresolvedCoefficient);
    const double represented=static_cast<double>(candidate);
    if(candidate>0&&represented==0)
      return Rejected("positive magnetic wave budget underflows storage at "+BinContext(key),StatusCode::UnresolvedCoefficient);
    BinLedger bin;bin.fieldLine=std::get<0>(key);bin.segment=std::get<1>(key);
    bin.branch=static_cast<RC::Branch>(std::get<2>(key));bin.helicity=static_cast<RC::Helicity>(std::get<3>(key));
    bin.bin=std::get<4>(key);bin.beforeMagneticJ=stored;bin.afterMagneticJ=represented;
    bin.requestedMagneticChangeJ=delta;
    const long double actual=static_cast<long double>(represented)-static_cast<long double>(stored);
    bin.magneticRoundoffJ=actual-delta;
    // A rounding ledger is evidence, not authorization to lose an accepted
    // transfer. Resolve the requested NET debit/credit in its own SI scale,
    // independently of a possibly enormous reservoir. The driver must change
    // representation or batch resolution after UnresolvedCoefficient; this
    // helper cannot erase, suppress, clip or redraw an accepted event.
    const double representationError=Dimensionless(bin.magneticRoundoffJ,std::fabs(delta));
    if((delta!=0&&actual==0)||!std::isfinite(representationError)||representationError>phase.relativeTolerance)
      return Rejected("magnetic transfer is unresolved at the explicit per-bin numerical tolerance: "+
          BinContext(key),StatusCode::UnresolvedCoefficient);
    requestedMagnetic.Add(delta);actualMagnetic.Add(actual);roundoff.Add(bin.magneticRoundoffJ);
    stored=represented;bins.push_back(bin);
  }
  Ledger ledger;
  ledger.particleEnergyChangeJ=particleEnergy.Value();
  ledger.actualParticleEnergyChangeJ=actualParticleEnergy.Value();ledger.particleRoundoffJ=particleRoundoff.Value();
  ledger.maximumParticleRelativeRoundoff=maximumParticleRelativeRoundoff;
  ledger.requestedMagneticChangeJ=requestedMagnetic.Value();
  ledger.actualMagneticChangeJ=actualMagnetic.Value();ledger.magneticRoundoffJ=roundoff.Value();
  ledger.physicalTransferResidualJ=ledger.particleEnergyChangeJ+2*ledger.requestedMagneticChangeJ;
  ledger.actualTotalResidualJ=ledger.actualParticleEnergyChangeJ+2*ledger.actualMagneticChangeJ;
  ledger.dimensionlessTransferResidual=Dimensionless(ledger.physicalTransferResidualJ,absoluteTransferScale.Value());
  ledger.dimensionlessActualResidual=Dimensionless(ledger.actualTotalResidualJ,absoluteTransferScale.Value());
  if(!std::isfinite(ledger.dimensionlessTransferResidual)||
      ledger.dimensionlessTransferResidual>phase.relativeTolerance)
    return Rejected("dimensionless particle/total-wave transfer residual exceeds the explicit numerical tolerance");
  if(!std::isfinite(ledger.dimensionlessActualResidual)||ledger.dimensionlessActualResidual>phase.relativeTolerance)
    return Rejected("actual encoded-particle/stored-wave residual exceeds the explicit numerical tolerance",
                    StatusCode::UnresolvedCoefficient);

  PW::EvolutionMetadata metadata;metadata.epochS=before.epochS;metadata.validUntilS=before.validUntilS;
  metadata.geometryGeneration=before.geometryGeneration;metadata.backgroundIdentity=before.backgroundIdentity;
  metadata.operationProvenance=phase.operationProvenance+"; frozen_batch_phase="+phase.identity+
      "; no-coupled-CTMC-or-PIC-publication-claim";
  const auto candidate=PW::PrepareEvolutionCandidate(before,updatedCells,metadata);
  if(!candidate.status.ok())return Rejected("canonical wave candidate rejected: "+candidate.status.message);
  Prepared result;result.status=Status::Ok();result.candidate=candidate.snapshot;result.operation=candidate.operation;
  result.particles=std::move(staged);result.bins=std::move(bins);result.events=std::move(eventLedgers);result.ledger=ledger;
  return result;
}

} } } // namespace SEP::Transport::FrozenCollisionTransaction
