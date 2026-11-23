#ifndef SEP_ADAPTERS_SWCME1D_MATERIAL_ADAPTER_H
#define SEP_ADAPTERS_SWCME1D_MATERIAL_ADAPTER_H

#include "../util/sep_common_header_path.h"
#include SRCSEP_SEP_COMMON_HEADER(sep_material_background.h)
#include "swcme_material_upstream_fwd.hpp"

namespace SEP { namespace SW1DAdapter {

// Detached Stage 12 export for an explicitly built, immutable canonical event.
// The caller supplies its background/geometry generations. This interval gives
// checked smooth HCI primitives, derivatives, trajectories and RH limits. It
// does not register native geometry or publish a particle BackgroundSnapshot.
Core::Result<std::shared_ptr<const Background::Material::Interval>> CaptureUpstreamMaterial(
    std::shared_ptr<const swcme::material1d::Event> event,
    std::uint64_t background_generation,std::uint64_t geometry_generation);

} } // namespace SEP::SW1DAdapter
#endif
