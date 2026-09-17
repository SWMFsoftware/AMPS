// ============================================================================
// Phase-O lifecycle coordinator.
//
// Pure sampling/publication/restart functions do not mutate Runtime. This
// small coordinator is the explicit orchestration layer that binds them to
// output cadence and checkpoint transitions after the host particle step has
// joined all mover workers.
// ============================================================================

#ifndef SEP3D_OUTPUT_OUTPUT_COORDINATOR_H
#define SEP3D_OUTPUT_OUTPUT_COORDINATOR_H

#include "publication.h"
#include "restart.h"

namespace SEP3D {
namespace Output {

struct DuePublicationResult {
  Core::Status status;
  bool published = false;
  PublicationResult publication;
};

DuePublicationResult PublishIfDue(RuntimeModel::Runtime* runtime,
                                  double simulationTimeS,
                                  const SamplingRequest& request,
                                  const std::string& codeIdentity,
                                  const std::string& snapshotFingerprint);

// WriteRestartAtBoundary records the post-commit checkpoint sequence in the
// file, then commits the same increment to Runtime. A failed write invokes the
// sole rollback transition and leaves the sequence unchanged.
Core::Status WriteRestartAtBoundary(RuntimeModel::Runtime* runtime,
                                    const std::string& path,
                                    RestartState state);

Core::Status RestoreRestartBeforeMesh(RuntimeModel::Runtime* runtime,
                                      const std::string& path,
                                      const RestartLoadOptions& options,
                                      RestartState* output);

}  // namespace Output
}  // namespace SEP3D

#endif  // SEP3D_OUTPUT_OUTPUT_COORDINATOR_H
