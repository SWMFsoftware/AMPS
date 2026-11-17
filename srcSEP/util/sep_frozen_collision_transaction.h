#ifndef SRCSEP_UTIL_SEP_FROZEN_COLLISION_TRANSACTION_H
#define SRCSEP_UTIL_SEP_FROZEN_COLLISION_TRANSACTION_H

#include "sep_polarized_wave_state.h"
#include <cstdint>
#include <memory>
#include <limits>
#include <string>
#include <vector>

namespace SEP { namespace Transport { namespace FrozenCollisionTransaction {

// This is a neutral, prepare-only bookkeeping boundary. Proposals have already
// been accepted by an explicitly frozen-wave operator. Their hazards were not
// recomputed against intermediate wave changes, so the batch is NOT an exact
// coupled particle/wave CTMC. A future PIC transaction must also own its actual
// particle lists, RNG/hazard histories, geometry and wave publication together.
enum class Phase { Unspecified, AcceptedFrozenWaveBatch };
struct PhaseContract {
  Phase phase=Phase::Unspecified;
  std::string identity;
  double epochS=0.0;
  std::uint64_t geometryGeneration=0,snapshotRevision=0;
  std::string snapshotIdentity,backgroundIdentity;
  // This is a numerical residual limit, not a physical domain or energy floor.
  // It must be supplied explicitly; neither an absolute SI floor nor a
  // scale=max(1,...) acceptance shortcut is introduced by this helper.
  double relativeTolerance=std::numeric_limits<double>::quiet_NaN();
  std::string operationProvenance;
};

struct EventKey {
  std::uint64_t particleIdentity=0,eventIndex=0;
};
struct ParticleRecord {
  // Identity zero is the native/source lifecycle's missing-identity sentinel.
  // nextEventIndex zero is a valid first accepted event.
  std::uint64_t identity=0,nextEventIndex=0,fieldLine=0,segment=0;
  double statisticalWeight=0.0; // Number of real particles represented.
  ResonantCollision::Particle phasePoint; // Plasma-frame SI momentum.
  // Opaque identity supplied by the enclosing lifecycle owner. This helper
  // never hashes, reseeds, resets, or manufactures an RNG/hazard state. The
  // expected identity in every proposal must match the untouched owner token.
  std::string lifecycleIdentity;
};
struct WeightedProposal {
  EventKey key;
  std::uint64_t fieldLine=0,segment=0;
  double statisticalWeight=0.0;
  // Invalid until explicitly supplied: even a zero-transfer proposal must not
  // acquire an assumed propagation direction from a default constructor.
  ResonantCollision::Branch branch=static_cast<ResonantCollision::Branch>(-1);
  std::string phaseIdentity,expectedLifecycleIdentity;
  // ProposeJump already scales both sides by statisticalWeight. It is checked
  // by the canonical collision validator and must NEVER be weighted again.
  ResonantCollision::Proposal proposal;
};

struct BinLedger {
  std::uint64_t fieldLine=0,segment=0;
  ResonantCollision::Branch branch=ResonantCollision::Branch::AlongB;
  ResonantCollision::Helicity helicity=ResonantCollision::Helicity::Left;
  std::size_t bin=0;
  double beforeMagneticJ=0.0,afterMagneticJ=0.0;
  long double requestedMagneticChangeJ=0.0L;
  // Actual stored change minus requested change. This is representation
  // roundoff, not physical heat, a limiter or a suppressed accepted event.
  // Its magnitude must resolve this bin's requested net transfer at the
  // explicit numerical tolerance; a sub-ULP lost transfer is rejected.
  long double magneticRoundoffJ=0.0L;
};
struct Ledger {
  // Canonical invariant increment requested by the accepted proposals.
  long double particleEnergyChangeJ=0.0L;
  // Independently evaluated change of the actual encoded momentum records.
  // Its difference from the canonical increment is numerical phase rounding;
  // it cannot authorize an unresolved or zero represented particle transfer.
  long double actualParticleEnergyChangeJ=0.0L,particleRoundoffJ=0.0L;
  long double requestedMagneticChangeJ=0.0L;
  long double actualMagneticChangeJ=0.0L;
  long double magneticRoundoffJ=0.0L;
  // Ideal nondispersive Alfvén TOTAL wave energy is twice magnetic energy.
  long double physicalTransferResidualJ=0.0L;
  long double actualTotalResidualJ=0.0L;
  double dimensionlessTransferResidual=0.0;
  double dimensionlessActualResidual=0.0,maximumParticleRelativeRoundoff=0.0;
};
struct ParticleLedger {
  EventKey key;
  double statisticalWeight=0.0;
  long double requestedEnergyChangeJ=0.0L,actualEnergyChangeJ=0.0L;
  long double particleRoundoffJ=0.0L;
  double dimensionlessRoundoff=0.0;
};
struct Prepared {
  Status status=Status::Error(StatusCode::InvalidArgument,"frozen collision batch was not prepared");
  std::shared_ptr<const Turbulence::Polarized::Snapshot> candidate;
  Turbulence::Polarized::EvolutionOperationLedger operation;
  std::vector<ParticleRecord> particles;
  std::vector<BinLedger> bins;
  std::vector<ParticleLedger> events;
  Ledger ledger;
  bool qualifiedForProduction=false;
};

// Inputs are immutable. Events are ordered by (particleIdentity,eventIndex),
// then reduced into the canonical (line,segment,branch,helicity,bin) ledger.
// The whole batch's net wave budgets must remain nonnegative. An insufficient
// budget rejects the complete batch; no clipping, event dropping, redraw, or
// partial candidate escapes on any failure. There is deliberately no Commit
// API or second wave/particle owner here. PrepareEvolutionCandidate alone owns
// the canonical snapshot identity, revision and structural validation.
Prepared Prepare(const Turbulence::Polarized::Snapshot& before,
                 const std::vector<ParticleRecord>& particles,
                 const std::vector<WeightedProposal>& proposals,
                 const PhaseContract& phase);

} } } // namespace SEP::Transport::FrozenCollisionTransaction
#endif
