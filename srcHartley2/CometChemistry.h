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
