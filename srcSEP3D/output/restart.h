// ============================================================================
// Phase-O complete, canonical restart state.
//
// The keyed stochastic method has no hidden generator object: its future is
// determined by campaignSeed plus each particle's stableId/completedStep/
// substep/purpose tuple. Those fields, all lifecycle/output counters, source
// and snapshot generations, sampling counters, and the particle ledger are
// serialized explicitly in a versioned little-endian format.
// ============================================================================

#ifndef SEP3D_OUTPUT_RESTART_H
#define SEP3D_OUTPUT_RESTART_H

#include "sampling.h"
#include "../runtime/runtime.h"
#include "../adapters/source_runtime.h"

#include <cstdint>
#include <functional>
#include <string>
#include <vector>

namespace SEP3D {
namespace Output {

struct RestartState {
  std::string configurationFingerprint;
  std::string resolvedConfigurationManifest;
  std::string storageLayoutFingerprint;
  std::string codeIdentity;
  std::string snapshotFingerprint;
  RuntimeModel::RuntimeCounters runtimeCounters;
  RuntimeModel::EventSchedule eventSchedule;
  RuntimeModel::SnapshotDescriptor activeSnapshot;
  double baseTimeStepS = 0.0;
  std::uint64_t backgroundGeneration = 0;
  std::uint64_t turbulenceGeneration = 0;
  std::uint64_t sourceGeneration = 0;
  std::uint64_t campaignSeed = 0;
  std::uint64_t nextStableParticleId = 0;
  std::uint64_t savedRankCount = 1;
  Adapters::ShockState shockState;
  SamplingState samplingState;
  std::vector<Adapters::ParticleRecord> particles;
  std::vector<Adapters::LedgerRow> ledgerRows;
  std::vector<Adapters::SourceLedgerRow> sourceLedgerRows;
};

enum class MissingSnapshotPolicy { Reject, Wait };
enum class RepartitionPolicy { RequireSameRankCount, DeterministicByStableId };

struct RestartLoadOptions {
  std::string expectedConfigurationFingerprint;
  std::string expectedCodeIdentity;
  std::string expectedSnapshotFingerprint;
  std::string expectedResolvedConfigurationManifest;
  std::string expectedStorageLayoutFingerprint;
  std::uint64_t availableBackgroundGeneration = 0;
  MissingSnapshotPolicy missingSnapshot = MissingSnapshotPolicy::Reject;
  std::uint64_t waitTimeoutMilliseconds = 0;
  std::uint64_t pollMilliseconds = 1;
  std::uint64_t currentRankCount = 1;
  RepartitionPolicy repartition = RepartitionPolicy::RequireSameRankCount;
  // The coupled host may poll its snapshot store. The callback must be
  // read-only and return true only when the requested immutable generation is
  // fully published. It is used only under MissingSnapshotPolicy::Wait.
  std::function<bool(std::uint64_t)> snapshotAvailable;
};

Core::Status WriteRestart(const std::string& path,
                          const RestartState& state);

// Loading is transactional: output is not modified until magic, schema,
// checksum, fingerprints, ranges, particle identities, ledger closure, and
// snapshot availability all validate.
Core::Status ReadRestart(const std::string& path,
                         const RestartLoadOptions& options,
                         RestartState* output);

}  // namespace Output
}  // namespace SEP3D

#endif  // SEP3D_OUTPUT_RESTART_H
