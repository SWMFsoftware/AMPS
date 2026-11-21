#ifndef _GENERIC_COMET_CHEMISTRY_H_
#define _GENERIC_COMET_CHEMISTRY_H_

#include <cstdint>
#include <string>
#include <utility>
#include <vector>

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

// A two-layer comet build compiles reaction topology (species, stoichiometry,
// rate-law kind, and excess energy) into the executable but supplies the
// case's rate coefficient at startup.  These setters deliberately permit only
// A to change and only before GetNetwork() constructs its lookup index; a
// runtime file therefore cannot alter compiled storage layout or chemistry
// topology.  Legacy compile-time inputs need not call them.
std::size_t GetCompiledReactionCount();
const char* GetCompiledReactionId(std::size_t reaction);
bool ConfigureReactionRate(const std::string &id,double rateCoefficient,
    std::string *error=NULL);
// Chemistry variants (srcHartley2/CHEMISTRY_VARIANTS.md).  A reaction group
// holds alternative complete definitions of one physical process, for example
// OH_LOSS = STATE_RESOLVED | SA1982_VH; exactly one variant per group is
// active in a run.  GetNetwork() sets the rate coefficient of every reaction
// of an inactive variant to zero when it constructs the network, and the
// per-parent channel cache keeps only positive-rate channels in input order,
// so an inactive variant never takes part in channel selection and the
// active variant's selection thresholds and random-number use are unchanged.
// The compiled tables below are generated from the build input by
// input/CometReactionNetwork.pm.  Ungrouped reactions are always active.
namespace RateRole {
enum Type {None=0,Dissociation=1,Ionization=2};
}
struct cCompiledGroup {
  const char *Id;
  int DefaultVariant;              // index into VariantNames
  int VariantCount;
  const char *const *VariantNames;
};
struct cCompiledRateTable {
  const char *Id;
  const char *Model;               // e.g. "SA1982_OH_LIFETIME"
  double ScaleS;                   // [s] per table unit (model parameter)
  double IonFraction;              // k_ion/k_d (model parameter)
};
struct cCompiledReactionVariant {
  int Group;                       // compiled group index, -1 = ungrouped
  int Variant;                     // index into that group's VariantNames
  int RateTable;                   // compiled table index, -1 = input rate
  int RateRole;                    // RateRole::Type
};
std::size_t GetCompiledGroupCount();
const cCompiledGroup& GetCompiledGroup(std::size_t group);
std::size_t GetCompiledRateTableCount();
const cCompiledRateTable& GetCompiledRateTable(std::size_t table);
const cCompiledReactionVariant& GetCompiledReactionVariant(std::size_t reaction);
// Select the active variant of a group before the network is first accessed
// (runtime input).  source describes where the choice came from, for the
// provenance record and error messages ("command line", "runtime input
// i3.in:123", "compiled default").  Without a call the compiled default is
// active.
bool SetActiveVariant(std::size_t group,int variant,const std::string &source,
    std::string *error=NULL);
int GetActiveVariant(std::size_t group);
const std::string& GetActiveVariantSource(std::size_t group);
// True when the reaction is ungrouped or belongs to its group's active variant.
bool IsReactionVariantActive(std::size_t reaction);
// Provenance of the effective variants, rate tables, and table-derived rates,
// set once by the runtime input and written to the run log and to the
// brightness output headers.  Each entry is (Tecplot aux-data name, value).
void SetVariantProvenance(
    const std::vector<std::pair<std::string,std::string>> &records);
const std::vector<std::pair<std::string,std::string>>& GetVariantProvenance();
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
