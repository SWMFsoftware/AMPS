#ifndef SRCSEP_FLUX_TUBE_CHECKPOINT_STATE_H
#define SRCSEP_FLUX_TUBE_CHECKPOINT_STATE_H
#include "util/sep_flux_tube_state.h"
#include <memory>
namespace SEP { namespace FieldLine { namespace FluxTubeGeometry {
class CheckpointRestorePlan {
 public:
  CheckpointRestorePlan();~CheckpointRestorePlan();
  CheckpointRestorePlan(CheckpointRestorePlan&&)noexcept;
  CheckpointRestorePlan& operator=(CheckpointRestorePlan&&)noexcept;
  CheckpointRestorePlan(const CheckpointRestorePlan&)=delete;
  CheckpointRestorePlan& operator=(const CheckpointRestorePlan&)=delete;
 private:
  struct Implementation;std::unique_ptr<Implementation> implementation_;
  friend Transport::Status StageCheckpointRestore(const FluxTubeState::State&,CheckpointRestorePlan*);
  friend Transport::Status ValidateCheckpointReady(const CheckpointRestorePlan&);
  friend bool CommitCheckpointRestore(CheckpointRestorePlan*)noexcept;
};
// Quiescent local owner API. Capture is observational and never initializes a
// lazy flux slot. Opaque callbacks reject rather than serializing code pointers
// or clearing an unapproved normalization. The enclosing checkpoint performs
// cross-rank agreement and keeps native writers quiescent through commit.
bool HasOpaqueAreaProfile();
Transport::Status CaptureCheckpointState(FluxTubeState::State* output);
Transport::Status ValidateCheckpointState(const FluxTubeState::State& state);
Transport::Status StageCheckpointRestore(const FluxTubeState::State& state,CheckpointRestorePlan* plan);
Transport::Status ValidateCheckpointReady(const CheckpointRestorePlan& plan);
[[nodiscard]] bool CommitCheckpointRestore(CheckpointRestorePlan* plan)noexcept;
} } }
#endif
