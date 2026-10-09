#ifndef SEP_COMMON_PERPENDICULAR_DIFFUSION_INTERNAL_H
#define SEP_COMMON_PERPENDICULAR_DIFFUSION_INTERNAL_H

#include "perpendicular_diffusion.h"

namespace SEP {
namespace PerpendicularDiffusion {
namespace Internal {

ModelResult EvaluateModel(const ParticleState& particle,
                          const LocalState& local,
                          const ModelConfiguration& configuration);

}  // namespace Internal
}  // namespace PerpendicularDiffusion
}  // namespace SEP

#endif
