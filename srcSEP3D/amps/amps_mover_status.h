// ============================================================================
// srcSEP3D/amps/amps_mover_status.h
//
// Single translation boundary between srcSEP3D mover outcomes and the AMPS
// particle-mover ABI.  No numerical kernel may return or copy _PARTICLE_*
// integers directly; doing so was the source of the prototype collision that
// Phase R0 removes.
//
// LAYER: L2 (AMPS adapter).  This header intentionally includes pic.h and must
// never be included by core/ or background/.
// ============================================================================

#ifndef SEP3D_AMPS_MOVER_STATUS_H
#define SEP3D_AMPS_MOVER_STATUS_H

#include "pic.h"

#include "../core/sep3d_types.h"

#include <cstdlib>

namespace SEP3D {
namespace AMPS {

// These assertions make an AMPS ABI change a compilation error at the adapter
// boundary.  The exact values are part of the current AMPS mover contract and
// are also checked from the source tree by BLDL3D03.  The distinctness checks
// document the minimum semantic requirement even if AMPS later renumbers the
// constants and this adapter is deliberately updated at the same time.
static_assert(_PARTICLE_DELETED_ON_THE_FACE_ == 0,
              "AMPS deleted-on-face return code changed; review srcSEP3D adapter");
static_assert(_PARTICLE_LEFT_THE_DOMAIN_ == 2,
              "AMPS left-domain return code changed; review srcSEP3D adapter");
static_assert(_PARTICLE_MOTION_FINISHED_ == 3,
              "AMPS motion-finished return code changed; review srcSEP3D adapter");
static_assert(_PARTICLE_LEFT_THE_DOMAIN_ != _PARTICLE_MOTION_FINISHED_,
              "AMPS mover return codes must be distinct");

// Convert a model-level outcome only at the point where the AMPS mover API
// requires an int.  There is intentionally no default mapping for invalid
// background data, numerical underflow, or internal errors: those are model
// statuses that must be handled explicitly before this function is called.
constexpr int ToAmpsMoverReturnCode(
    Core::ParticleMotionOutcome outcome) noexcept {
  switch (outcome) {
    case Core::ParticleMotionOutcome::Advanced:
      return _PARTICLE_MOTION_FINISHED_;
    case Core::ParticleMotionOutcome::LeftDomain:
      return _PARTICLE_LEFT_THE_DOMAIN_;
  }

  // Reaching this point requires a corrupted enum value.  Returning a valid
  // AMPS code would hide the corruption and could mutate particle lists, so
  // terminate immediately instead of applying a fallback.
  std::abort();
}

} // namespace AMPS
} // namespace SEP3D

#endif // SEP3D_AMPS_MOVER_STATUS_H
