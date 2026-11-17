#ifndef _AMPS_PIC_FIELD_LINE_POPULATION_H_
#define _AMPS_PIC_FIELD_LINE_POPULATION_H_

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <set>
#include <vector>

namespace PIC { namespace ParticleSplitting { namespace FledLine {

// These host-only hooks belong to population control, not to CloneParticle.
// A byte-exact clone also represents migration/restart of the *same* physical
// history, whereas an exact split creates a second statistical representative.
// Install hooks before particle stepping and do not change them concurrently
// with a split. Applications without extension semantics may leave hooks NULL.
using fParticleIdentity = std::uint64_t (*)(long int particle);
using fPrepareExactSplit = void (*)(long int parent);
using fPostExactSplit = void (*)(long int parent, long int child);
extern fParticleIdentity ParticleIdentity;
extern fPrepareExactSplit PrepareExactSplit;
extern fPostExactSplit PostExactSplit;

// Optional application-owned conservative merger. A registered controller owns
// every outcome, including explicit no-reduction for sparse/rank-limited bins;
// it must never fall through to the legacy velocity-average merger after such a
// result. Registration is serial before stepping. Species/list/bin arguments
// retain the existing core interface, and CloneParticle remains byte-exact.
using fMergeController = void (*)(int species,long int& head,
    int spatialBins,int parallelVelocityBins,int normalVelocityBins,
    double spatialRange,int minimumParticles,int maximumParticles,int threshold);
extern fMergeController MergeController;

namespace ExactSplit {

struct Candidate {
  long int particle;
  double weight;
  std::uint64_t identity;
};

enum class Status { Success, InvalidWeight, InvalidIdentity, DuplicateIdentity };

// Read/validate every application key before any weight or particle-buffer
// mutation. A supplied key is a persistent physical identity: zero and duplicate
// values cannot define a canonical tie break and are rejected, never replaced by
// pointers, list order, rank, or a newly sampled random value. The key callback
// may also preflight application lifecycle counters before a possible split.
inline Status Prepare(std::vector<Candidate>& candidates,
                      fParticleIdentity identity) {
  std::set<std::uint64_t> identities;
  for (auto& candidate : candidates) {
    if (!std::isnormal(candidate.weight) || candidate.weight <= 0.0)
      return Status::InvalidWeight;
    if (identity) {
      candidate.identity = identity(candidate.particle);
      if (candidate.identity == 0) return Status::InvalidIdentity;
      if (!identities.insert(candidate.identity).second)
        return Status::DuplicateIdentity;
    }
  }
  // A physical key makes equal-weight selection independent of the intrusive
  // list/allocation order. Without a key preserve encounter order explicitly;
  // generic core code cannot invent an application's missing identity contract.
  std::stable_sort(candidates.begin(), candidates.end(),
      [identity](const Candidate& a, const Candidate& b) {
        if (a.weight != b.weight) return a.weight > b.weight;
        return identity && a.identity < b.identity;
      });
  return Status::Success;
}

// This is the exact kernel used at the native PIC split site. Operations supplies
// SetWeight(parent,w), Allocate(), and Clone(child,parent); the native Clone
// preserves destination list linkage while copying every requested byte.
// Halving a positive normal weight is exact in binary floating point unless it
// enters the subnormal range, in which case the existing controller skips it.
// Unchanged phase-space coordinates and child weights summing to the old weight
// preserve every instantaneous weighted moment, including relativistic energy.
// The no-fail post hook receives the true parent before list traversal can lose
// that relationship. It may update application identity/stochastic metadata but
// must leave the two half weights and physical phase-space state unchanged.
template<class Operations>
inline bool Apply(const Candidate& candidate, Operations& operations,
                  fPostExactSplit postSplit,
                  fPrepareExactSplit prepareSplit = nullptr) {
  const double halfWeight = candidate.weight / 2.0;
  if (!std::isnormal(halfWeight) || halfWeight <= 0.0) return false;
  // Application preflight may reject an exhausted lineage counter, colliding
  // child identity or incompatible stochastic state by throwing/failing before
  // either weight is changed or a buffer record is allocated. It must not
  // consume a random variate or mutate the parent; postSplit then cannot fail.
  if (prepareSplit) prepareSplit(candidate.particle);
  operations.SetWeight(candidate.particle, halfWeight);
  const long int child = operations.Allocate();
  operations.Clone(child, candidate.particle);
  if (postSplit) postSplit(candidate.particle, child);
  return true;
}

}  // namespace ExactSplit
}}}  // namespace PIC::ParticleSplitting::FledLine

#endif
