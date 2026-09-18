#include "output_coordinator.h"

#include <cmath>

namespace SEP3D {
namespace Output {

DuePublicationResult PublishIfDue(RuntimeModel::Runtime* runtime,
                                  double simulationTimeS,
                                  const SamplingRequest& request,
                                  const std::string& codeIdentity,
                                  const std::string& snapshotFingerprint) {
  DuePublicationResult result;
  if (runtime == nullptr || !runtime->configuration() ||
      runtime->state() != RuntimeModel::LifecycleState::SnapshotReady ||
      !std::isfinite(simulationTimeS)) {
    result.status = Core::Status(Core::StatusCode::InvalidTransition,
                                 "output publication requires SnapshotReady Runtime");
    return result;
  }
  // CompleteStep sets stepsSinceOutput to zero exactly when it increments the
  // output sequence. This check contains no independent cadence counter.
  if (runtime->counters().stepsSinceOutput != 0 ||
      runtime->counters().completedSteps == 0) {
    result.status = Core::Status::OK();
    return result;
  }
  const SamplingSnapshot snapshot = Sample(request);
  if (!snapshot.status.ok()) { result.status = snapshot.status; return result; }
  PublicationMetadata metadata;
  metadata.sequence = runtime->counters().outputSequence;
  metadata.simulationTimeS = simulationTimeS;
  metadata.configurationFingerprint =
      runtime->configuration()->physics_fingerprint();
  metadata.codeIdentity = codeIdentity;
  metadata.snapshotFingerprint = snapshotFingerprint;
  const RuntimeModel::SnapshotDescriptor* active = runtime->active_snapshot();
  if (active == nullptr) {
    result.status = Core::Status(Core::StatusCode::SnapshotUnavailable,
                                 "output requires an active snapshot");
    return result;
  }
  metadata.snapshotGeneration = active->generation;
  const auto& options = runtime->configuration()->options();
  result.publication = Publish(options.outputDirectory, options.outputPrefix,
                               metadata, snapshot);
  result.status = result.publication.status;
  result.published = result.status.ok();
  return result;
}

Core::Status WriteRestartAtBoundary(RuntimeModel::Runtime* runtime,
                                    const std::string& path,
                                    RestartState state) {
  if (runtime == nullptr || !runtime->configuration())
    return Core::Status(Core::StatusCode::InvalidInput,
                        "checkpoint requires a configured Runtime");
  Core::Status status = runtime->BeginCheckpoint();
  if (!status.ok()) return status;
  state.configurationFingerprint =
      runtime->configuration()->physics_fingerprint();
  state.resolvedConfigurationManifest =
      runtime->configuration()->resolved_manifest();
  state.storageLayoutFingerprint =
      runtime->configuration()->storage_layout().fingerprint;
  state.runtimeCounters = runtime->counters();
  state.eventSchedule = runtime->event_schedule();
  state.baseTimeStepS =
      runtime->configuration()->options().requestedTimeStepS;
  if (runtime->active_snapshot() != nullptr)
    state.activeSnapshot = *runtime->active_snapshot();
  ++state.runtimeCounters.checkpointSequence;
  status = WriteRestart(path, state);
  if (!status.ok()) {
    const Core::Status rolledBack = runtime->AbortCheckpoint();
    return rolledBack.ok() ? status : rolledBack;
  }
  status = runtime->CompleteCheckpoint();
  return status;
}

Core::Status RestoreRestartBeforeMesh(RuntimeModel::Runtime* runtime,
                                      const std::string& path,
                                      const RestartLoadOptions& options,
                                      RestartState* output) {
  if (runtime == nullptr || output == nullptr || !runtime->configuration() ||
      runtime->state() != RuntimeModel::LifecycleState::Configured)
    return Core::Status(Core::StatusCode::InvalidTransition,
                        "restart restore requires Configured Runtime before mesh binding");
  RestartState candidate;
  Core::Status status = ReadRestart(path, options, &candidate);
  if (!status.ok()) return status;
  status = runtime->RestoreCounters(candidate.runtimeCounters);
  if (!status.ok()) return status;
  status = runtime->RestoreEventSchedule(candidate.eventSchedule);
  if (!status.ok()) return status;
  *output = candidate;
  return Core::Status::OK();
}

}  // namespace Output
}  // namespace SEP3D
