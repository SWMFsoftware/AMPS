// ============================================================================
// Counter-based random streams for three-dimensional SEP transport.
//
// A particle history must not depend on its current MPI rank, OpenMP worker,
// or position in a cell-linked list.  Every variate is therefore a pure
// function of a semantic key and a draw counter.  Adding a future stochastic
// process uses a new RandomPurpose and cannot shift the Parker or pitch-angle
// streams.  The counter is serializable, so a restarted particle consumes the
// exact next variate used by an uninterrupted trajectory.
// ============================================================================

#ifndef SEP3D_TRANSPORT_KEYED_RANDOM_H
#define SEP3D_TRANSPORT_KEYED_RANDOM_H

#include <cstdint>

namespace SEP3D {
namespace Transport {

enum class RandomPurpose : std::uint64_t {
  ParkerParallel = 1,
  FocusedPitch = 2,
  SourceCount = 3,
  SourceSpectrum = 4,
  SourcePitch = 5,
  SourceGyrophase = 6,
  SourcePosition = 7,
  // V01 assigns new stable values without renumbering any released stream.
  // Each stochastic operator owns a stream so enabling cross-field diffusion
  // cannot change the Parker-parallel or focused-pitch random history.
  PerpendicularFirst = 8,
  PerpendicularSecond = 9,
  ReservedFuturePhysics = 1024
};

struct RandomKey {
  std::uint64_t campaignSeed = 0;
  std::uint64_t particleId = 0;
  std::uint64_t step = 0;
  std::uint64_t substep = 0;
  RandomPurpose purpose = RandomPurpose::ParkerParallel;
};

struct RandomState {
  RandomKey key;
  std::uint64_t drawCounter = 0;
};

class KeyedRandomStream final {
 public:
  explicit KeyedRandomStream(const RandomKey& key,
                             std::uint64_t drawCounter = 0)
      : state_{key, drawCounter} {}

  // UniformOpen01 excludes both endpoints.  That makes log(u) safe in the
  // Box-Muller transform and prevents a source inverse CDF from selecting a
  // formally open endpoint because an integer converted to exactly zero.
  double UniformOpen01();
  double Normal01();

  const RandomState& state() const { return state_; }
  void Restore(const RandomState& state) { state_ = state; }

  // Exposed for reproducibility tests and artifact fingerprints.  SplitMix64
  // is used as a deterministic bit mixer; this is a simulation generator, not
  // a cryptographic primitive.
  static std::uint64_t Hash(const RandomKey& key, std::uint64_t drawCounter);

 private:
  RandomState state_;
};

}  // namespace Transport
}  // namespace SEP3D

#endif  // SEP3D_TRANSPORT_KEYED_RANDOM_H
