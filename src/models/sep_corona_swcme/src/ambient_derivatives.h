#ifndef SEP_CORONA_SWCME_AMBIENT_DERIVATIVES_H
#define SEP_CORONA_SWCME_AMBIENT_DERIVATIVES_H

#include <functional>

#include "sep_corona_swcme/ambient_state.h"

namespace SEP::CoronaSwcme::Detail {

// Differentiate an already sampled center and its provider at one epoch.
// The callback owns the physical support checks for each displaced point.
// The stencil also checks magnetic sector and region at every node so that
// both providers use exactly the same one-sided handoff convention.
Core::Result<AmbientState> DifferentiateAmbient(CoronalCME::Vec3 point,
    AmbientState center,double relativeStep,
    const std::function<Core::Result<AmbientPrimitive>(CoronalCME::Vec3)>& evaluate,
    bool requireConvergedNeighbors=false);

}
#endif
