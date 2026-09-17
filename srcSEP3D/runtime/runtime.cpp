#include "runtime.h"

#include <cmath>
#include <sstream>

namespace SEP3D {
namespace RuntimeModel {

namespace {

Core::Status InvalidTransition(LifecycleState actual,
                               const char* operation,
                               LifecycleState required) {
  std::ostringstream message;
  message << operation << " requires " << Name(required)
          << " but Runtime is " << Name(actual);
  return Core::Status(Core::StatusCode::InvalidTransition, message.str());
}

Core::Status InvalidSnapshot(const std::string& message) {
  return Core::Status(Core::StatusCode::SnapshotUnavailable, message);
}

}  // namespace

const char* Name(LifecycleState state) {
  switch (state) {
    case LifecycleState::Created: return "Created";
    case LifecycleState::Configured: return "Configured";
    case LifecycleState::MeshReady: return "MeshReady";
    case LifecycleState::WaitingForSnapshot: return "WaitingForSnapshot";
    case LifecycleState::SnapshotReady: return "SnapshotReady";
    case LifecycleState::Running: return "Running";
    case LifecycleState::Checkpointing: return "Checkpointing";
    case LifecycleState::Finalized: return "Finalized";
  }
  return "Unknown";
}

bool SnapshotDescriptor::Covers(double timeS) const {
  return std::isfinite(timeS) && timeS >= validFromS && timeS <= validUntilS;
}

Core::Status Runtime::RequireState(LifecycleState required,
                                   const char* operation) const {
  if (state_ != required) return InvalidTransition(state_, operation, required);
  return Core::Status::OK();
}

Core::Status Runtime::Configure(
    const std::shared_ptr<const RunConfiguration3D>& configuration) {
  const Core::Status order = RequireState(LifecycleState::Created, "Configure");
  if (!order.ok()) return order;
  if (!configuration) {
    return Core::Status(Core::StatusCode::InvalidInput,
                        "Configure requires an immutable configuration");
  }

  // Commit after all checks.  No layout buffer or provider object is allocated
  // here; the host may inspect the frozen layout before creating an AMPS mesh.
  configuration_ = configuration;
  state_ = LifecycleState::Configured;
  return Core::Status::OK();
}

Core::Status Runtime::BindMesh(const MeshBinding& binding) {
  const Core::Status order = RequireState(LifecycleState::Configured, "BindMesh");
  if (!order.ok()) return order;
  if (binding.layout != configuration_->storage_layout()) {
    return Core::Status(
        Core::StatusCode::LayoutMismatch,
        "mesh storage layout differs from the configuration frozen before allocation");
  }
  meshBinding_ = binding;
  state_ = LifecycleState::MeshReady;
  return Core::Status::OK();
}

Core::Status Runtime::BeginBackgroundAcquisition(AdapterKind adapter) {
  const Core::Status order = RequireState(
      LifecycleState::MeshReady, "BeginBackgroundAcquisition");
  if (!order.ok()) return order;
  const BackgroundAuthority expected = configuration_->options().background;
  const bool matches =
      (expected == BackgroundAuthority::AnalyticParker &&
       adapter == AdapterKind::Standalone) ||
      (expected == BackgroundAuthority::Swmf && adapter == AdapterKind::Swmf);
  if (!matches) {
    return Core::Status(
        Core::StatusCode::ConfigurationConflict,
        "background adapter does not match the configured authority");
  }
  adapterKind_ = adapter;
  state_ = LifecycleState::WaitingForSnapshot;
  return Core::Status::OK();
}

Core::Status Runtime::PublishSnapshot(const SnapshotDescriptor& candidate) {
  if (state_ != LifecycleState::WaitingForSnapshot &&
      state_ != LifecycleState::SnapshotReady) {
    return InvalidTransition(state_, "PublishSnapshot",
                             LifecycleState::WaitingForSnapshot);
  }
  if (!candidate.complete) return InvalidSnapshot("snapshot is incomplete");
  if (!std::isfinite(candidate.epochS) ||
      !std::isfinite(candidate.validFromS) ||
      !std::isfinite(candidate.validUntilS) ||
      candidate.validUntilS < candidate.validFromS ||
      candidate.epochS < candidate.validFromS ||
      candidate.epochS > candidate.validUntilS) {
    return InvalidSnapshot("snapshot validity interval is not finite and ordered");
  }
  if (candidate.generation == 0) {
    return InvalidSnapshot("snapshot generation zero is reserved");
  }
  if (candidate.authority != configuration_->options().background) {
    return InvalidSnapshot("snapshot authority differs from configuration");
  }
  if (candidate.configurationFingerprint !=
      configuration_->physics_fingerprint()) {
    return InvalidSnapshot("snapshot configuration fingerprint mismatch");
  }
  if (candidate.coordinateFrame.empty() || candidate.providerIdentity.empty()) {
    return InvalidSnapshot("snapshot frame and provider identity are required");
  }
  if (hasSnapshot_ && candidate.generation <= activeSnapshot_.generation) {
    return InvalidSnapshot("snapshot generation must increase monotonically");
  }

  activeSnapshot_ = candidate;
  hasSnapshot_ = true;
  state_ = LifecycleState::SnapshotReady;
  return Core::Status::OK();
}

Core::Status Runtime::BeginStep(double simulationTimeS) {
  const Core::Status order = RequireState(LifecycleState::SnapshotReady,
                                          "BeginStep");
  if (!order.ok()) return order;
  if (!hasSnapshot_ || !activeSnapshot_.Covers(simulationTimeS)) {
    return InvalidSnapshot("active snapshot does not cover the requested step time");
  }
  pinnedSnapshotGeneration_ = activeSnapshot_.generation;
  state_ = LifecycleState::Running;
  return Core::Status::OK();
}

Core::Status Runtime::CompleteStep() {
  const Core::Status order = RequireState(LifecycleState::Running,
                                          "CompleteStep");
  if (!order.ok()) return order;

  RuntimeCounters next = counters_;
  ++next.completedSteps;
  ++next.stepsSinceOutput;
  if (next.stepsSinceOutput >= configuration_->options().outputCadenceSteps) {
    ++next.outputSequence;
    next.stepsSinceOutput = 0;
  }
  counters_ = next;
  pinnedSnapshotGeneration_ = 0;
  state_ = LifecycleState::SnapshotReady;
  return Core::Status::OK();
}

Core::Status Runtime::BeginCheckpoint() {
  const Core::Status order = RequireState(LifecycleState::SnapshotReady,
                                          "BeginCheckpoint");
  if (!order.ok()) return order;
  state_ = LifecycleState::Checkpointing;
  return Core::Status::OK();
}

Core::Status Runtime::CompleteCheckpoint() {
  const Core::Status order = RequireState(LifecycleState::Checkpointing,
                                          "CompleteCheckpoint");
  if (!order.ok()) return order;
  ++counters_.checkpointSequence;
  state_ = LifecycleState::SnapshotReady;
  return Core::Status::OK();
}

Core::Status Runtime::RestoreCounters(const RuntimeCounters& counters) {
  const Core::Status order = RequireState(LifecycleState::Configured,
                                          "RestoreCounters");
  if (!order.ok()) return order;
  if (counters.stepsSinceOutput >=
      configuration_->options().outputCadenceSteps) {
    return Core::Status(
        Core::StatusCode::InvalidInput,
        "restored stepsSinceOutput must be below output cadence");
  }
  counters_ = counters;
  return Core::Status::OK();
}

Core::Status Runtime::Finalize() {
  if (state_ == LifecycleState::Created ||
      state_ == LifecycleState::Running ||
      state_ == LifecycleState::Checkpointing ||
      state_ == LifecycleState::Finalized) {
    return Core::Status(
        Core::StatusCode::InvalidTransition,
        std::string("Finalize is not legal from ") + Name(state_));
  }
  state_ = LifecycleState::Finalized;
  pinnedSnapshotGeneration_ = 0;
  return Core::Status::OK();
}

}  // namespace RuntimeModel
}  // namespace SEP3D
