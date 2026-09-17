// ============================================================================
// srcSEP3D typed Runtime lifecycle
//
// Runtime owns orchestration state only: immutable configuration, the frozen
// pre-mesh layout, active snapshot metadata, and restartable counters.  It does
// not allocate an AMPS mesh, parse files, or evaluate transport physics.  Each
// mutating method validates its complete precondition before changing a field,
// which makes an invalid call transactional and straightforward to test.
// ============================================================================

#ifndef SEP3D_RUNTIME_RUNTIME_H
#define SEP3D_RUNTIME_RUNTIME_H

#include "run_configuration.h"

#include <cstdint>
#include <memory>
#include <string>

namespace SEP3D {
namespace RuntimeModel {

enum class LifecycleState {
  Created,
  Configured,
  MeshReady,
  WaitingForSnapshot,
  SnapshotReady,
  Running,
  Checkpointing,
  Finalized
};

const char* Name(LifecycleState state);

enum class AdapterKind { Standalone, Swmf };

struct MeshBinding {
  // Exact byte contract used by the mesh allocator. Runtime copies it only
  // after equality with RunConfiguration3D::storage_layout() is established.
  StorageLayout layout;
};

struct SnapshotDescriptor {
  BackgroundAuthority authority = BackgroundAuthority::AnalyticParker;
  double epochS = 0.0;       // snapshot epoch [s] in the host simulation clock
  double validFromS = 0.0;   // inclusive validity start [s]
  double validUntilS = 0.0;  // inclusive, finite validity end [s]
  std::uint64_t generation = 0;
  bool complete = false;
  std::string coordinateFrame;
  std::string providerIdentity;
  std::string configurationFingerprint;

  bool Covers(double timeS) const;
};

struct RuntimeCounters {
  // These values are checkpoint state. Output cadence is reconstructed from
  // stepsSinceOutput rather than a file-scope/static modulo counter, so a
  // restarted run produces the same next sequence number as an uninterrupted run.
  std::uint64_t completedSteps = 0;
  std::uint64_t stepsSinceOutput = 0;
  std::uint64_t outputSequence = 0;
  std::uint64_t checkpointSequence = 0;
};

class Runtime final {
 public:
  Runtime() = default;
  Runtime(const Runtime&) = delete;
  Runtime& operator=(const Runtime&) = delete;

  LifecycleState state() const { return state_; }
  const std::shared_ptr<const RunConfiguration3D>& configuration() const {
    return configuration_;
  }
  const RuntimeCounters& counters() const { return counters_; }
  bool has_snapshot() const { return hasSnapshot_; }
  const SnapshotDescriptor* active_snapshot() const {
    return hasSnapshot_ ? &activeSnapshot_ : nullptr;
  }
  std::uint64_t pinned_snapshot_generation() const {
    return pinnedSnapshotGeneration_;
  }

  Core::Status Configure(
      const std::shared_ptr<const RunConfiguration3D>& configuration);
  Core::Status BindMesh(const MeshBinding& binding);
  Core::Status BeginBackgroundAcquisition(AdapterKind adapter);
  Core::Status PublishSnapshot(const SnapshotDescriptor& candidate);
  Core::Status BeginStep(double simulationTimeS);
  Core::Status CompleteStep();
  Core::Status BeginCheckpoint();
  Core::Status CompleteCheckpoint();
  Core::Status RestoreCounters(const RuntimeCounters& counters);
  Core::Status Finalize();

 private:
  Core::Status RequireState(LifecycleState required,
                            const char* operation) const;

  LifecycleState state_ = LifecycleState::Created;
  std::shared_ptr<const RunConfiguration3D> configuration_;
  MeshBinding meshBinding_;
  AdapterKind adapterKind_ = AdapterKind::Standalone;
  SnapshotDescriptor activeSnapshot_;
  bool hasSnapshot_ = false;
  std::uint64_t pinnedSnapshotGeneration_ = 0;
  RuntimeCounters counters_;
};

}  // namespace RuntimeModel
}  // namespace SEP3D

#endif  // SEP3D_RUNTIME_RUNTIME_H
