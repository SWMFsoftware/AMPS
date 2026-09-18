// ============================================================================
// R06 first-class observer orchestration.
//
// Geometry/configuration is resolved at an integer tick, after particle motion
// and source injection have joined.  PreparePublication is read-only;
// CommitPublication is the only operation that clears a pending window.  A
// failed file publication can therefore call AbortPublication (or simply do
// nothing) and retry the identical deterministic reduction.
// ============================================================================

#ifndef SEP3D_OUTPUT_OBSERVER_RUNTIME_H
#define SEP3D_OUTPUT_OBSERVER_RUNTIME_H

#include "sampling.h"
#include "../runtime/run_configuration.h"

namespace SEP3D {
namespace Output {

Core::Status BuildObserverDefinitions(
    const RuntimeModel::RunConfiguration3D& configuration,
    double simulationTimeS,
    std::vector<VirtualSpacecraftDefinition>* definitions);

class ObserverRuntime final {
 public:
  // Capture one joined-boundary observation window.  A second capture before
  // commit is rejected so the caller cannot silently overwrite unpublished
  // statistics.
  Core::Status Capture(SamplingRequest request);
  SamplingSnapshot PreparePublication() const;
  Core::Status CommitPublication(const SamplingSnapshot& published);
  void AbortPublication() {}

  bool has_pending_window() const { return hasPending_; }
  const SamplingState& state() const { return state_; }
  void RestoreState(const SamplingState& state) { state_ = state; }

 private:
  SamplingState state_;
  SamplingRequest pending_;
  bool hasPending_ = false;
};

}  // namespace Output
}  // namespace SEP3D

#endif  // SEP3D_OUTPUT_OBSERVER_RUNTIME_H
