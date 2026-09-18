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

// Snapshot updates have their own transaction state because transport may
// continue to reference the active generation while a coupled provider fills
// the inactive buffer.  The main lifecycle therefore remains SnapshotReady
// until one collective publication commits the staged descriptor.
enum class SnapshotUpdateState {
  Idle,
  Requested,
  Filling,
  Staged,
  Failed
};

const char* Name(SnapshotUpdateState state);

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
  std::uint64_t currentTick = 0;
};

struct EventSchedule {
  std::uint64_t nextBackgroundTick = 0;
  std::uint64_t nextInjectionTick = 0;
  std::uint64_t nextSamplingTick = 0;
  // UINT64_MAX denotes a disabled periodic checkpoint.
  std::uint64_t nextCheckpointTick = UINT64_MAX;
};

enum class ScheduledEvent { Background, Injection, Sampling, Checkpoint };

struct ClockObservation {
  double picTimeS = 0.0;
  double picTimeStepS = 0.0;
  double snapshotTimeS = 0.0;
  double shockTimeS = 0.0;
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
  SnapshotUpdateState snapshot_update_state() const {
    return snapshotUpdateState_;
  }
  const EventSchedule& event_schedule() const { return eventSchedule_; }
  double CurrentTimeS() const;
  double NextStepEndTimeS() const;
  bool EventDue(ScheduledEvent event) const;

  Core::Status Configure(
      const std::shared_ptr<const RunConfiguration3D>& configuration);
  Core::Status BindMesh(const MeshBinding& binding);
  Core::Status BeginBackgroundAcquisition(AdapterKind adapter);
  Core::Status PublishSnapshot(const SnapshotDescriptor& candidate);
  Core::Status RequestSnapshotUpdate(double epochS,
                                     std::uint64_t generation);
  Core::Status BeginSnapshotFill();
  Core::Status StageSnapshot(const SnapshotDescriptor& candidate);
  Core::Status PublishStagedSnapshot(bool collectiveReady);
  Core::Status FailSnapshotUpdate(const std::string& reason);
  Core::Status AcknowledgeSnapshotFailure();
  Core::Status VerifyClockAgreement(const ClockObservation& observation) const;
  Core::Status BeginStep(double simulationTimeS);
  Core::Status CompleteStep();
  Core::Status BeginCheckpoint();
  Core::Status CompleteCheckpoint();
  // Return from a failed transactional checkpoint without incrementing its
  // sequence. This is the only legal rollback transition in Runtime.
  Core::Status AbortCheckpoint();
  Core::Status RestoreCounters(const RuntimeCounters& counters);
  Core::Status RestoreEventSchedule(const EventSchedule& schedule);
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
  EventSchedule eventSchedule_;
  SnapshotUpdateState snapshotUpdateState_ = SnapshotUpdateState::Idle;
  SnapshotDescriptor stagedSnapshot_;
  bool hasStagedSnapshot_ = false;
  double requestedSnapshotEpochS_ = 0.0;
  std::uint64_t requestedSnapshotGeneration_ = 0;
  std::string snapshotUpdateFailure_;
};

}  // namespace RuntimeModel
}  // namespace SEP3D

#endif  // SEP3D_RUNTIME_RUNTIME_H
