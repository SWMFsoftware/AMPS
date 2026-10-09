#ifndef SEP_COMMON_PARALLEL_DIFFUSION_PARALLEL_DIFFUSION_ADVANCED_H
#define SEP_COMMON_PARALLEL_DIFFUSION_PARALLEL_DIFFUSION_ADVANCED_H

// Internal linkage contract between the small explicit-model dispatcher and
// the spectral/nonlinear implementation.  Applications include only
// parallel_diffusion.h; keeping these declarations private prevents numerical
// helpers from becoming an accidental public API.

#include "parallel_diffusion.h"

namespace SEP {
namespace ParallelDiffusion {
namespace Internal {

bool IsAdvancedModel(ModelId model);
Status ValidateAdvancedConfiguration(const ModelConfiguration& configuration);
ParallelResult EvaluateAdvanced(const ParticleState& particle,
                                const LocalState& local,
                                const ModelConfiguration& configuration);
Status BuildAdvancedConfiguration(const std::string& modelId,
                                  const std::vector<InputParameter>& parameters,
                                  ModelConfiguration* configuration);
Status EvaluateAdvancedPitchAngle(double mu,
                                  const ParticleState& particle,
                                  const LocalState& local,
                                  const ModelConfiguration& configuration,
                                  double* dMuMuPerS);

}  // namespace Internal
}  // namespace ParallelDiffusion
}  // namespace SEP

#endif
