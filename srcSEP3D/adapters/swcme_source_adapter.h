// ============================================================================
// Phase-A SWCME -> shared SEP injection adapter.
//
// SWCME owns shock reconstruction and DSA source parameters. sep_common owns
// the normalized sampling algebra. This adapter is intentionally thin: it
// copies a validated swcme::sep::SEPSourceState into a dimension-independent
// event record, derives the unambiguous dN/dp power law, and samples all 3-D
// launch variables from purpose-separated keyed streams. No shock or spectrum
// formula is duplicated here.
// ============================================================================

#ifndef SEP3D_ADAPTERS_SWCME_SOURCE_ADAPTER_H
#define SEP3D_ADAPTERS_SWCME_SOURCE_ADAPTER_H

#include "transport_adapter.h"

#include "sep_injection_spectrum.h"
#include "swcme_sep_source.hpp"

#include <cstdint>
#include <string>

namespace SEP3D {
namespace Adapters {

// Immutable description of one active shock patch. eventGeneration is the
// coupled snapshot/shock generation, not a loop index; it therefore remains
// stable when a surface is redistributed across MPI ranks.
struct ShockSourceRecord {
  Core::Status status;
  bool active = false;
  std::uint64_t eventGeneration = 0;
  std::uint64_t sourceId = 0;
  Core::Vec3 positionM;
  Core::Vec3 outwardNormal;
  double relativePatchWeight = 0.0;
  double compression = 1.0;
  double shockNormalSpeedMPerS = 0.0;
  SEP::Injection::Configuration injection;
  std::string sourceFingerprint;
};

struct InjectedParticle {
  Core::Status status;
  ParticleRecord particle;
};

// Translate a common SWCME source exactly once. For an isotropic DSA source
// f(p)~p^-q, the probability density per linear momentum is dN/dp~p^(2-q).
// sep_common represents a density as x^(-powerIndex), hence powerIndex=q-2.
ShockSourceRecord MakeShockSourceRecord(
    const swcme::sep::SEPSourceState& source,
    std::uint64_t eventGeneration,
    std::uint64_t campaignSeed,
    std::uint64_t macroparticlesPerEvent,
    double injectionEfficiency);

// Sample one macroparticle with independently keyed momentum, pitch-angle,
// and gyrophase streams. Adding a future position draw cannot perturb the
// existing sequence. statisticalWeight is the patch's relative source weight
// multiplied by efficiency and divided equally over the event population.
InjectedParticle SampleInjectedParticle(const ShockSourceRecord& source,
                                         std::uint64_t macroIndex,
                                         int species);

}  // namespace Adapters
}  // namespace SEP3D

#endif  // SEP3D_ADAPTERS_SWCME_SOURCE_ADAPTER_H
