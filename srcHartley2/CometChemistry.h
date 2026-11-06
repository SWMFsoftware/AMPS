#ifndef _GENERIC_COMET_CHEMISTRY_H_
#define _GENERIC_COMET_CHEMISTRY_H_

#include <cstdint>
#include <string>

#if __has_include("models/chemical_solver/ReactionNetwork.h")
#include "models/chemical_solver/ReactionNetwork.h"
#elif __has_include("../models/chemical_solver/ReactionNetwork.h")
#include "../models/chemical_solver/ReactionNetwork.h"
#else
#error "Cannot locate models/chemical_solver/ReactionNetwork.h"
#endif

template <class T> class cTreeNodeAMR;
namespace PIC { namespace Mesh { class cDataBlockAMR; } }

namespace Comet {
namespace Chemistry {

extern bool Enabled;
extern double HeliocentricDistanceAU;

namespace ResidualSampling {
enum Type {IndependentBernoulli=0,RandomPhaseSystematic=1};
}

// The source default remains independent Bernoulli sampling, so existing
// comet inputs reserve no block storage and retain their exact random-number
// stream. The systematic mode uses a separate counter-based phase stream;
// ResidualPhaseSeed therefore never consumes a physical-event random draw.
extern int ResidualSamplingMode;
extern std::uint64_t ResidualPhaseSeed;

// Stage-6 packet splitting is independently default-off.  When enabled, the
// generated policy supplies one complete-reaction packet count per stable
// reaction ID.  MaximumDaughterPacketCount bounds N; the radial hard product
// limit separately bounds N times the concrete stoichiometric product count.
extern bool DaughterPacketSplittingEnabled;
extern long int MaximumDaughterPacketCount;
long int GetDaughterPacketCount(std::size_t reaction);
int GetDaughterPacketPolicySource(std::size_t reaction);

// Stage-7 channel importance sampling changes only the conditional selection
// of a mutually exclusive channel after the ordinary total-reaction event has
// occurred.  Complete parent families are generated from stable reaction IDs.
// Their requested q* values remain positive and normalized; every selected
// product carrier receives L=q/q*.  Unconfigured families and the global
// default use the exact historical analog selector and L=1.
extern bool ChannelImportanceSamplingEnabled;
bool IsChannelImportanceConfigured(std::size_t reaction);
int GetChannelParentSpecies(std::size_t reaction);
double GetPhysicalChannelProbability(std::size_t reaction);
double GetBiasedChannelProbability(std::size_t reaction);
double GetChannelLikelihoodRatio(std::size_t reaction);

// Stage 8 replaces only the Bernoulli decision for explicitly configured
// parent species. One carrier is split into the exact survivor and reacted
// fractions. The reacted carrier remains a parent until the mover reaches a
// conditionally sampled reaction point, then becomes a product (or a neutral
// sink) and uses the rest of the normalized local-timestep interval. There is
// deliberately no user-supplied biased occurrence probability.
extern bool TotalReactionForcingEnabled;
bool IsTotalReactionForcingConfigured(int species);

// Stage-9 estimator arms use one frozen Stage-6/7/8-enabled executable.  The
// configured mode preserves the input-selected behavior.  Analog keeps the
// identical particle layout and installed mover callbacks but executes N=1,
// physical channel probabilities, unit likelihood, and ordinary stochastic
// reaction occurrence.  Enhanced executes the compiled policies.
namespace EstimatorMode {
enum Type {Configured=0,Analog=1,Enhanced=2};
}
extern int ActiveEstimatorMode;
bool ConfigureEstimatorMode(int mode,std::string *error=NULL);
const char* GetEstimatorModeName();
bool IsDaughterPacketSplittingActive();
bool IsChannelImportanceSamplingActive();
bool IsTotalReactionForcingActive();

// Install the three application callbacks used by the core first-order
// trajectory-event hook. This must run after particle-data offsets have been
// requested but before PIC::Mover::Init validates the all-or-none callback
// contract.
void ConfigureTotalReactionForcingMover();

// The returned object is a non-owning view of static tables generated from
// the selected input file.  Validate must succeed before particle stepping.
const ChemicalSolver::cReactionNetwork& GetNetwork();
bool Validate(std::string *error=NULL);

// Reserve owner-block accumulator storage after core local timestep/weight
// storage has been registered and before any AMR block is allocated.
void RequestBlockData();

// An accepted regular carrier transaction can change chemistry W/dt ratios.
// Close the current bounded-error epoch, then rephase lazily under the new
// carriers. Emergency common-factor scaling deliberately does not call this
// function because it preserves every source/product carrier ratio exactly.
void PrepareForRegularCarrierChange();

// Close only the diagnostic activity interval at a complete sampling break;
// physical phases are retained so systematic accumulation continues across
// output windows.
void CloseResidualDiagnosticWindow();

// Exercise the production residual-resampling and product-velocity helpers for
// every generated neutral channel after AMPS has initialized the configured
// species weights and local timesteps.  The total event count is distributed
// across MPI ranks for each reaction.  The returned pass/fail value is
// collective and identical on every rank.
bool ValidateEnsemble(long int eventsPerReaction);

// Exercise Stage-6 packet kinematics, exact per-product correction splitting,
// and N=1/2/4/8 variance/cost behavior for every generated neutral reaction.
// It also drives one repeated-product reaction through the real executor to
// test parent-slot reuse, linked-list commit, processed flags, represented
// packet moments, and byte-exact particle rollback on capacity deferral. This
// focused mode does not replace a transported comet science pilot.
bool ValidateDaughterPacketSplitting(long int eventsPerReaction);

// Exercise Stage-7 channel selection using deterministic global event IDs.
// The validator compares analog and biased channel ensembles, includes
// product-free sinks, checks likelihood-weighted product moments and family
// ESS, and is decomposition-identical at 1/4/12 ranks.  The separate Stage-6
// validator remains the real K*N particle-transaction check and, when Stage 7
// is enabled, also proves that stored corrections include L before /N.
bool ValidateChannelImportanceSampling(long int eventsPerFamily);

// Exercise the deterministic carrier split, conditional position inversion,
// Stage-7 likelihood/Stage-6 packet composition, rollback contract, and
// mover-event continuation with decomposition-independent fixtures.  The
// optional Stage-9 arm continues from the exact live mover products into the
// production general-NTC kernel; the default preserves the accepted Stage-8
// test and its output schema.
bool ValidateTotalReactionForcing(long int eventsPerFixture,
    bool validateTransportedCollision=false);

// Read-only checkpoint audit used by Stage 9.  The void pointer keeps this
// geometry-neutral interface independent of pic.h; callers pass an AMPS
// particle-data record.  A correction-off/Stage-8-off build has no reserved
// marker and therefore always returns false.
bool HasPendingTotalReactionEvent(const void *particleData);

// Exercise owner-only block state, common packet phases, transaction
// rollback, dynamic MPI block migration, and a real local-weight restart.
// This focused validator replaces the science run; it must not be invoked in
// the middle of a production calculation.
bool ValidateSystematicResidualState();

// Reserve one application byte in every particle before ParticleBuffer::Init.
// Transactional chemistry uses it only while one photochemical stage is in
// progress to distinguish completed parents/new daughters from deferred
// parents that must be retried after a capacity-control checkpoint.
void RequestParticleData();
void BeginBatchedModel();
bool IsBatchedParticleEligible(long int ptr);
long int GetBatchedParticleLimit();
void FinishBatchedPass(long int attemptedParents,long int remainingParents);

// AMPS calls this once for each particle during the photolytic-reaction pass.
// The routine samples one mutually exclusive channel, creates all daughters
// from the same event, and conserves event momentum and excess energy.
bool ProcessParticle(long int ptr,long int& firstParticleCell,
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node);

}
}

#endif
