#ifndef SEP_UTIL_SEP_POPULATION_CONTROL_H
#define SEP_UTIL_SEP_POPULATION_CONTROL_H

#include "sep_common_header_path.h"
#include SRCSEP_SEP_COMMON_HEADER(sep_transport_common.h)

#include <cstdint>
#include <string>
#include <vector>

namespace SEP {
namespace PopulationControl {

// This dependency-light record is an observation of a host PIC particle, not
// a second particle implementation.  Extensive quantities are represented by
// macro-particle weight times the per-physical-particle value.  All quantities
// use SI units and species is part of the merge compatibility contract.
struct ParticleRecord {
  std::uint64_t stableId = 0;
  std::uint64_t parentId = 0;
  std::uint64_t lineageGeneration = 0;
  int species = -1;
  double weight = 0.0;
  double chargeC = 0.0;
  double kineticEnergyJ = 0.0;
  double parallelMomentumKgMPerS = 0.0;
  double mu = 0.0;
};

struct Configuration {
  int minimumParticlesPerCell = 600;
  int maximumParticlesPerCell = 1000;
  int spatialBins = 20;
  int momentumBins = 20;
  int pitchBins = 20;
  bool deterministicSelection = true;
  double relativeInvariantTolerance = 1.0e-12;
  std::uint64_t lineageSchema = UINT64_C(1);
};

Transport::Status ValidateConfiguration(const Configuration& configuration);
Transport::Status ValidateParticle(const ParticleRecord& particle);

struct Moments {
  double number = 0.0;
  double chargeC = 0.0;
  double energyJ = 0.0;
  double parallelMomentumKgMPerS = 0.0;
  double muFirst = 0.0;
  double muSecond = 0.0;
};

struct InvariantReport {
  Transport::Status status;
  Moments before;
  Moments after;
  double maximumRelativeError = 0.0;
  std::vector<std::string> failedInvariants;
};

Moments ComputeMoments(const std::vector<ParticleRecord>& particles);
InvariantReport CompareInvariants(const std::vector<ParticleRecord>& before,
                                  const std::vector<ParticleRecord>& after,
                                  double relativeTolerance);

// Child identifiers are a tagged, order-sensitive hash of the complete
// lineage tuple.  Rank and thread are intentionally absent: decomposition is
// computational metadata and must not alter physical ancestry or RNG keys.
std::uint64_t DeriveChildId(std::uint64_t parentId,
                            std::uint64_t lineageGeneration,
                            std::uint64_t childOrdinal,
                            std::uint64_t lineageSchema);

// Splitting clones phase-space state and partitions only statistical weight.
// The final child receives the floating residual, which preserves total weight
// exactly even when the requested child count is not a power of two.
Transport::Status SplitParticle(const ParticleRecord& parent,
                                std::size_t childCount,
                                std::uint64_t lineageSchema,
                                std::vector<ParticleRecord>* children);

// Merging is defined only for a single species and charge state.  Weighted
// first and second pitch moments cannot both be preserved by one representative
// scalar mu in general; this function therefore preserves number, charge,
// energy, parallel momentum, and the first pitch moment exactly and exposes the
// second-moment change through CompareInvariants for policy enforcement.
Transport::Status MergeParticles(const std::vector<ParticleRecord>& parents,
                                  std::uint64_t mergedOrdinal,
                                  std::uint64_t lineageSchema,
                                  ParticleRecord* merged);

}  // namespace PopulationControl
}  // namespace SEP

#endif  // SEP_UTIL_SEP_POPULATION_CONTROL_H
