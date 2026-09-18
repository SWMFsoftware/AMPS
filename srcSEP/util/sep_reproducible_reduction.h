#ifndef SEP_UTIL_SEP_REPRODUCIBLE_REDUCTION_H
#define SEP_UTIL_SEP_REPRODUCIBLE_REDUCTION_H

#include "sep_common_header_path.h"
#include SRCSEP_SEP_COMMON_HEADER(sep_transport_common.h)

#include <atomic>
#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace SEP {
namespace Reproducibility {

// Every additive particle-to-wave contribution carries a physical key.  Rank
// and thread number are deliberately absent: changing decomposition must not
// alter either the ordering or the identity of a contribution.
struct ContributionKey {
  // schema identifies the canonical key layout persisted by a restart.  A
  // restart with a different schema is rejected instead of silently changing
  // reduction order and therefore the rounded wave-energy result.
  std::uint64_t schema = 1;
  std::uint64_t source = 0;
  std::uint64_t fieldLine = 0;
  std::uint64_t segment = 0;
  std::uint64_t branch = 0;
  std::uint64_t spectralBin = 0;
  std::uint64_t species = 0;
  std::uint64_t particle = 0;
  std::uint64_t step = 0;
  std::uint64_t event = 0;
  std::uint64_t interval = 0;
  std::uint64_t purpose = 0;
};

struct Contribution {
  ContributionKey key;
  double waveEnergyJ = 0.0;
  double streaming = 0.0;
  std::uint64_t resonantCount = 0;
};

// A worker owns this append-only buffer during a particle loop.  No shared
// floating-point datum is written until all workers have reached the explicit
// deterministic reduction point.
class ThreadLocalBuffer {
 public:
  void Add(const Contribution& value);
  const std::vector<Contribution>& values() const { return values_; }
  void Clear();

 private:
  std::vector<Contribution> values_;
};

struct SegmentAccumulator {
  std::uint64_t fieldLine = 0;
  std::uint64_t segment = 0;
  std::uint64_t branch = 0;
  std::uint64_t spectralBin = 0;
  double waveEnergyJ = 0.0;
  double streaming = 0.0;
  std::uint64_t resonantCount = 0;
};

static const std::uint64_t ProductionContributionKeySchema = 1;

// CanonicalReduction flattens thread buffers, sorts by the complete physical
// contribution key, and accumulates with long-double intermediates.  The same
// routine is used after MPI gathers; MPI_Allreduce is reserved for scalar
// integer counters and other values whose operation is unambiguously additive.
Transport::Status CanonicalReduction(
    const std::vector<ThreadLocalBuffer>& workers,
    std::vector<SegmentAccumulator>* segments);

// This helper models an MPI gather with an arbitrary number of partitions.  It
// exists both as executable policy documentation and as a decomposition test:
// its output must match a one-partition CanonicalReduction bit for bit.
Transport::Status CanonicalPartitionReduction(
    const std::vector<std::vector<Contribution> >& partitions,
    std::vector<SegmentAccumulator>* segments);

// Purpose-separated random streams make stochastic results independent of
// optional diagnostics and unrelated operators.  The purpose identifier is a
// stable caller-owned enum value, never a thread/rank identifier.
Transport::KeyedRandomStream MakeRandomStream(std::uint64_t campaignSeed,
                                               std::uint64_t particleId,
                                               std::uint64_t step,
                                               std::uint64_t purpose);

// Warning and limiter statistics are integers, so atomic increments have no
// floating-point ordering ambiguity.  Snapshot is intended for diagnostics;
// Reset is called at the beginning of each iteration/checkpoint interval.
class AtomicCounters {
 public:
  void AddWarning(std::uint64_t count = 1);
  void AddLimiterActivation(std::uint64_t count = 1);
  std::uint64_t warnings() const;
  std::uint64_t limiterActivations() const;
  void Reset();

 private:
  std::atomic<std::uint64_t> warnings_{0};
  std::atomic<std::uint64_t> limiterActivations_{0};
};

// FNV-1a over a canonical hexadecimal rendering is a compact, portable piece
// of evidence for PAR05.  It is not a cryptographic integrity mechanism.
std::uint64_t EvidenceHash(const std::vector<SegmentAccumulator>& segments);
std::string EvidenceHashHex(const std::vector<SegmentAccumulator>& segments);

}  // namespace Reproducibility
}  // namespace SEP

#endif  // SEP_UTIL_SEP_REPRODUCIBLE_REDUCTION_H
