#include "runtime.h"

#include <algorithm>
#include <cmath>
#include <limits>
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

Core::Status ValidateSnapshotCandidate(
    const SnapshotDescriptor& candidate,
    const RunConfiguration3D& configuration,
    const SnapshotDescriptor* active) {
  if (!candidate.complete) return InvalidSnapshot("snapshot is incomplete");
  if (!std::isfinite(candidate.epochS) ||
      !std::isfinite(candidate.validFromS) ||
      !std::isfinite(candidate.validUntilS) ||
      candidate.validUntilS < candidate.validFromS ||
      candidate.epochS < candidate.validFromS ||
      candidate.epochS > candidate.validUntilS)
    return InvalidSnapshot("snapshot validity interval is not finite and ordered");
  if (candidate.generation == 0)
    return InvalidSnapshot("snapshot generation zero is reserved");
  if (candidate.authority != configuration.options().background)
    return InvalidSnapshot("snapshot authority differs from configuration");
  if (candidate.configurationFingerprint !=
      configuration.physics_fingerprint())
    return InvalidSnapshot("snapshot configuration fingerprint mismatch");
  if (candidate.coordinateFrame.empty() || candidate.providerIdentity.empty())
    return InvalidSnapshot("snapshot frame and provider identity are required");
  if (active != nullptr && candidate.generation <= active->generation)
    return InvalidSnapshot("snapshot generation must increase monotonically");
  return Core::Status::OK();
}

bool SameTime(double left, double right) {
  const double scale = std::max(1.0, std::max(std::fabs(left), std::fabs(right)));
  return std::fabs(left - right) <=
      64.0 * std::numeric_limits<double>::epsilon() * scale;
}

std::uint64_t AdvanceEvent(std::uint64_t current, std::uint64_t cadence) {
  if (cadence == 0 || current == UINT64_MAX ||
      current > UINT64_MAX - cadence) return UINT64_MAX;
  return current + cadence;
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

const char* Name(SnapshotUpdateState state) {
  switch (state) {
    case SnapshotUpdateState::Idle: return "Idle";
    case SnapshotUpdateState::Requested: return "Requested";
    case SnapshotUpdateState::Filling: return "Filling";
    case SnapshotUpdateState::Staged: return "Staged";
    case SnapshotUpdateState::Failed: return "Failed";
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
  eventSchedule_.nextBackgroundTick =
      configuration_->options().backgroundCadenceSteps;
  eventSchedule_.nextInjectionTick =
      configuration_->options().injectionCadenceSteps;
  eventSchedule_.nextSamplingTick =
      configuration_->options().outputCadenceSteps;
  eventSchedule_.nextCheckpointTick =
      configuration_->options().checkpointCadenceSteps == 0
          ? UINT64_MAX
          : configuration_->options().checkpointCadenceSteps;
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
  if (state_ != LifecycleState::WaitingForSnapshot) {
    return InvalidTransition(state_, "PublishSnapshot",
                             LifecycleState::WaitingForSnapshot);
  }
  const Core::Status valid = ValidateSnapshotCandidate(
      candidate, *configuration_, hasSnapshot_ ? &activeSnapshot_ : nullptr);
  if (!valid.ok()) return valid;

  activeSnapshot_ = candidate;
  hasSnapshot_ = true;
  state_ = LifecycleState::SnapshotReady;
  return Core::Status::OK();
}

Core::Status Runtime::RequestSnapshotUpdate(
    double epochS, std::uint64_t generation) {
  const Core::Status order = RequireState(LifecycleState::SnapshotReady,
                                          "RequestSnapshotUpdate");
  if (!order.ok()) return order;
  if (snapshotUpdateState_ != SnapshotUpdateState::Idle &&
      snapshotUpdateState_ != SnapshotUpdateState::Failed)
    return Core::Status(Core::StatusCode::InvalidTransition,
                        "a snapshot update transaction is already active");
  if (!std::isfinite(epochS) || generation == 0 ||
      (hasSnapshot_ && generation <= activeSnapshot_.generation))
    return InvalidSnapshot("requested snapshot epoch or generation is invalid");
  requestedSnapshotEpochS_ = epochS;
  requestedSnapshotGeneration_ = generation;
  snapshotUpdateFailure_.clear();
  hasStagedSnapshot_ = false;
  snapshotUpdateState_ = SnapshotUpdateState::Requested;
  return Core::Status::OK();
}

Core::Status Runtime::BeginSnapshotFill() {
  if (state_ != LifecycleState::SnapshotReady ||
      snapshotUpdateState_ != SnapshotUpdateState::Requested)
    return Core::Status(Core::StatusCode::InvalidTransition,
                        "BeginSnapshotFill requires a requested update at a joined boundary");
  snapshotUpdateState_ = SnapshotUpdateState::Filling;
  return Core::Status::OK();
}

Core::Status Runtime::StageSnapshot(const SnapshotDescriptor& candidate) {
  if (state_ != LifecycleState::SnapshotReady ||
      snapshotUpdateState_ != SnapshotUpdateState::Filling)
    return Core::Status(Core::StatusCode::InvalidTransition,
                        "StageSnapshot requires a filling update transaction");
  const Core::Status valid = ValidateSnapshotCandidate(
      candidate, *configuration_, hasSnapshot_ ? &activeSnapshot_ : nullptr);
  if (!valid.ok()) return valid;
  if (candidate.generation != requestedSnapshotGeneration_ ||
      !SameTime(candidate.epochS, requestedSnapshotEpochS_))
    return InvalidSnapshot(
        "staged snapshot does not match the requested epoch and generation");
  stagedSnapshot_ = candidate;
  hasStagedSnapshot_ = true;
  snapshotUpdateState_ = SnapshotUpdateState::Staged;
  return Core::Status::OK();
}

Core::Status Runtime::PublishStagedSnapshot(bool collectiveReady) {
  if (state_ != LifecycleState::SnapshotReady ||
      snapshotUpdateState_ != SnapshotUpdateState::Staged ||
      !hasStagedSnapshot_)
    return Core::Status(Core::StatusCode::InvalidTransition,
                        "PublishStagedSnapshot requires one complete staged generation");
  if (!collectiveReady)
    return InvalidSnapshot(
        "collective snapshot readiness failed; active generation is unchanged");
  activeSnapshot_ = stagedSnapshot_;
  hasSnapshot_ = true;
  hasStagedSnapshot_ = false;
  snapshotUpdateState_ = SnapshotUpdateState::Idle;
  return Core::Status::OK();
}

Core::Status Runtime::FailSnapshotUpdate(const std::string& reason) {
  if (state_ != LifecycleState::SnapshotReady ||
      (snapshotUpdateState_ != SnapshotUpdateState::Requested &&
       snapshotUpdateState_ != SnapshotUpdateState::Filling &&
       snapshotUpdateState_ != SnapshotUpdateState::Staged))
    return Core::Status(Core::StatusCode::InvalidTransition,
                        "FailSnapshotUpdate requires an active update transaction");
  if (reason.empty()) return InvalidSnapshot("snapshot failure reason is empty");
  stagedSnapshot_ = SnapshotDescriptor{};
  hasStagedSnapshot_ = false;
  snapshotUpdateFailure_ = reason;
  snapshotUpdateState_ = SnapshotUpdateState::Failed;
  return Core::Status::OK();
}

Core::Status Runtime::AcknowledgeSnapshotFailure() {
  if (snapshotUpdateState_ != SnapshotUpdateState::Failed)
    return Core::Status(Core::StatusCode::InvalidTransition,
                        "no failed snapshot update is awaiting acknowledgement");
  snapshotUpdateFailure_.clear();
  snapshotUpdateState_ = SnapshotUpdateState::Idle;
  return Core::Status::OK();
}

double Runtime::CurrentTimeS() const {
  return configuration_ ? counters_.currentTick *
      configuration_->options().requestedTimeStepS : 0.0;
}

double Runtime::NextStepEndTimeS() const {
  return configuration_ ? (counters_.currentTick + 1) *
      configuration_->options().requestedTimeStepS : 0.0;
}

bool Runtime::EventDue(ScheduledEvent event) const {
  const std::uint64_t tick = counters_.currentTick;
  if (tick == 0 || !configuration_) return false;
  const auto& options = configuration_->options();
  switch (event) {
    case ScheduledEvent::Background:
      return tick % options.backgroundCadenceSteps == 0;
    case ScheduledEvent::Injection:
      return tick % options.injectionCadenceSteps == 0;
    case ScheduledEvent::Sampling:
      return tick % options.outputCadenceSteps == 0;
    case ScheduledEvent::Checkpoint:
      return options.checkpointCadenceSteps != 0 &&
          tick % options.checkpointCadenceSteps == 0;
  }
  return false;
}

Core::Status Runtime::VerifyClockAgreement(
    const ClockObservation& observation) const {
  if (!configuration_ || !std::isfinite(observation.picTimeS) ||
      !std::isfinite(observation.picTimeStepS) ||
      !std::isfinite(observation.snapshotTimeS) ||
      !std::isfinite(observation.shockTimeS))
    return Core::Status(Core::StatusCode::InvalidInput,
                        "clock observation is incomplete or non-finite");
  const double expected = CurrentTimeS();
  const double dt = configuration_->options().requestedTimeStepS;
  if (!SameTime(observation.picTimeS, expected) ||
      !SameTime(observation.picTimeStepS, dt) ||
      !SameTime(observation.snapshotTimeS, expected) ||
      !SameTime(observation.shockTimeS, expected))
    return Core::Status(Core::StatusCode::ConfigurationConflict,
                        "Runtime, PIC, snapshot, and shock clocks disagree");
  return Core::Status::OK();
}

Core::Status Runtime::BeginStep(double simulationTimeS) {
  const Core::Status order = RequireState(LifecycleState::SnapshotReady,
                                          "BeginStep");
  if (!order.ok()) return order;
  if (!SameTime(simulationTimeS, CurrentTimeS())) {
    return Core::Status(Core::StatusCode::ConfigurationConflict,
                        "host simulation time differs from the authoritative Runtime tick");
  }
  if (!hasSnapshot_ || !activeSnapshot_.Covers(simulationTimeS) ||
      !activeSnapshot_.Covers(NextStepEndTimeS())) {
    return InvalidSnapshot("active snapshot does not cover the complete requested step");
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
  ++next.currentTick;
  ++next.stepsSinceOutput;
  if (next.stepsSinceOutput >= configuration_->options().outputCadenceSteps) {
    ++next.outputSequence;
    next.stepsSinceOutput = 0;
  }
  counters_ = next;
  const auto& options = configuration_->options();
  if (counters_.currentTick == eventSchedule_.nextBackgroundTick)
    eventSchedule_.nextBackgroundTick = AdvanceEvent(
        eventSchedule_.nextBackgroundTick, options.backgroundCadenceSteps);
  if (counters_.currentTick == eventSchedule_.nextInjectionTick)
    eventSchedule_.nextInjectionTick = AdvanceEvent(
        eventSchedule_.nextInjectionTick, options.injectionCadenceSteps);
  if (counters_.currentTick == eventSchedule_.nextSamplingTick)
    eventSchedule_.nextSamplingTick = AdvanceEvent(
        eventSchedule_.nextSamplingTick, options.outputCadenceSteps);
  if (counters_.currentTick == eventSchedule_.nextCheckpointTick)
    eventSchedule_.nextCheckpointTick = AdvanceEvent(
        eventSchedule_.nextCheckpointTick, options.checkpointCadenceSteps);
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

Core::Status Runtime::AbortCheckpoint() {
  const Core::Status order = RequireState(LifecycleState::Checkpointing,
                                          "AbortCheckpoint");
  if (!order.ok()) return order;
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
  if (counters.currentTick != counters.completedSteps) {
    return Core::Status(Core::StatusCode::InvalidInput,
                        "restored currentTick must equal completedSteps");
  }
  counters_ = counters;
  const auto& options = configuration_->options();
  auto nextAfter = [](std::uint64_t tick, std::uint64_t cadence) {
    if (cadence == 0) return UINT64_MAX;
    const std::uint64_t quotient = tick / cadence;
    if (quotient >= UINT64_MAX / cadence) return UINT64_MAX;
    return (quotient + 1) * cadence;
  };
  eventSchedule_.nextBackgroundTick = nextAfter(
      counters.currentTick, options.backgroundCadenceSteps);
  eventSchedule_.nextInjectionTick = nextAfter(
      counters.currentTick, options.injectionCadenceSteps);
  eventSchedule_.nextSamplingTick = nextAfter(
      counters.currentTick, options.outputCadenceSteps);
  eventSchedule_.nextCheckpointTick = nextAfter(
      counters.currentTick, options.checkpointCadenceSteps);
  return Core::Status::OK();
}

Core::Status Runtime::RestoreEventSchedule(const EventSchedule& schedule) {
  const Core::Status order = RequireState(LifecycleState::Configured,
                                          "RestoreEventSchedule");
  if (!order.ok()) return order;
  const std::uint64_t tick = counters_.currentTick;
  const auto& options = configuration_->options();
  auto validNext = [tick](std::uint64_t next, std::uint64_t cadence,
                          bool optional) {
    if (optional && cadence == 0) return next == UINT64_MAX;
    return cadence != 0 && next > tick && next % cadence == 0;
  };
  if (!validNext(schedule.nextBackgroundTick,
                 options.backgroundCadenceSteps, false) ||
      !validNext(schedule.nextInjectionTick,
                 options.injectionCadenceSteps, false) ||
      !validNext(schedule.nextSamplingTick,
                 options.outputCadenceSteps, false) ||
      !validNext(schedule.nextCheckpointTick,
                 options.checkpointCadenceSteps, true))
    return Core::Status(Core::StatusCode::InvalidInput,
                        "restored event schedule is not aligned after the current tick");
  eventSchedule_ = schedule;
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
