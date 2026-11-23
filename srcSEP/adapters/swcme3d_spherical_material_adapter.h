#ifndef SRCSEP_SWCME3D_SPHERICAL_MATERIAL_ADAPTER_H
#define SRCSEP_SWCME3D_SPHERICAL_MATERIAL_ADAPTER_H

#include "../util/sep_common_header_path.h"
#include SRCSEP_SEP_COMMON_HEADER(sep_material_background.h)

namespace SEP::CoronaSwcme::SphericalShockOnly {
class Event;
struct Numerics;
}

namespace SEP::SW3DSphericalAdapter {

// Exports one immutable HCI interval. This bridge deliberately has no PIC
// state, publication side effects, particle boundary or source authority.
Core::Result<std::shared_ptr<const Background::Material::Interval>> Capture(
    std::shared_ptr<const CoronaSwcme::SphericalShockOnly::Event> event,
    const CoronaSwcme::SphericalShockOnly::Numerics& numerics,
    double beginS, double endS, std::uint64_t backgroundGeneration,
    std::uint64_t geometryGeneration);

}  // namespace SEP::SW3DSphericalAdapter

#endif
