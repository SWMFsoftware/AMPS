#ifndef _GENERIC_COMET_CHEMISTRY_H_
#define _GENERIC_COMET_CHEMISTRY_H_

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

// The returned object is a non-owning view of static tables generated from
// the selected input file.  Validate must succeed before particle stepping.
const ChemicalSolver::cReactionNetwork& GetNetwork();
bool Validate(std::string *error=NULL);

// Exercise the production residual-resampling and product-velocity helpers for
// every generated neutral channel after AMPS has initialized the configured
// species weights and local timesteps.  The total event count is distributed
// across MPI ranks for each reaction.  The returned pass/fail value is
// collective and identical on every rank.
bool ValidateEnsemble(long int eventsPerReaction);

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
