#ifndef SRCSEP_ADAPTERS_PARALLEL_DIFFUSION_ADAPTER_H
#define SRCSEP_ADAPTERS_PARALLEL_DIFFUSION_ADAPTER_H

// srcSEP binding of the shared parallel-diffusion coefficient library
// (src/models/sep_common/parallel_diffusion, specification revision 1.4).
//
// Ownership:
//   * the input-file section grammar, model selection, per-model key schema,
//     numerical domains, coefficient formulas, and the mover-facing pointer
//     SEP::ParallelDiffusion::ActiveParallelDiffusion belong to the library;
//   * srcSEP owns only (a) locating the schema-4 [parallel_diffusion] section
//     in its --input file (util/sep_initialization.cpp), (b) the selection
//     cross-check with --spatial-diffusion-provider and the mover, (c) the
//     declaration of which optional library inputs srcSEP can supply, and
//     (d) packing the particle/background state of one field-line sample into
//     the library's SI ParticleState/LocalState.
//
// This header and its implementation contain no PIC, MPI, or AMPS code; the
// PIC-dependent sampling stays in coefficient_providers.cpp.  The header uses
// only C++11 types so C++11 unit-test builds of srcSEP code may include it.

#include "../util/sep_initialization.h"

#include <cstdint>
#include <string>

namespace SEP {
namespace ParallelDiffusionBinding {

// Validate the selection and install the model before AMPS initialization.
//
//   librarySelected true when the resolved --spatial-diffusion-provider is
//                   parallel-diffusion-library (SpatialKind::
//                   ParallelDiffusionLibrary).  The caller (main.cpp) owns
//                   the registry enum; a plain flag keeps this adapter
//                   independent of the coefficient-registry header;
//   moverName       canonical production mover name ("parker", ...);
//   initialization  the parsed --input configuration, or NULL when srcSEP runs
//                   without --input.
//
// The library is selected exactly when librarySelected is true, and then a schema-4 [parallel_diffusion] section is required.  A section
// without that selection, the selection without a section, or the selection
// with a focused-transport mover is a typed error.  When selected, the
// section body is parsed by SEP::ParallelDiffusion::ParseSection, checked
// against the inputs srcSEP can supply (none of the optional ones: no
// nucleon count, no slab/2-D turbulence decomposition, no external
// time/region/radial factors, no effective field), and installed with
// SetActiveConfiguration, which also sets ActiveParallelDiffusion.  The
// resolved model ID and library fingerprint are then written into
// *initialization for the startup fingerprint.  When not selected, nothing is
// installed and the function returns Ok.  Errors leave the library's active
// model unchanged.
Transport::Status Configure(bool librarySelected,
                            const std::string& moverName,
                            Initialization::Configuration* initialization);

// True after Configure() installed a model in this process.
bool IsInstalled();

struct KappaSample {
  Transport::Status status;
  double kappaParallelM2PerS = 0.0;  // [m^2 s^-1]
  double lambdaParallelM = 0.0;      // [m]
  std::string modelId;
  std::string configurationFingerprint;
};

// Evaluate kappa_parallel for one particle at one field-line sample through
// the library's ActiveParallelDiffusion pointer.  All inputs are SI:
//
//   massKg, signedChargeC     species rest mass [kg] and charge [C];
//   momentumKgMPerS           total momentum magnitude [kg m s^-1];
//   positionM[3]              Cartesian position [m] in srcSEP's
//                             heliocentric frame (Sun at (0,0,0); a field
//                             line need not start at the origin), so |x| is
//                             the heliocentric radius;
//   meanFieldT[3]             background magnetic field vector [T];
//   timeS                     simulation time of the background snapshot [s];
//   backgroundRevision        snapshot generation, copied to provenance.
//
// No nucleon count, turbulence decomposition, or external factor is passed;
// Configure() has already rejected every model that would need one.  The
// result fails (never substitutes a value) unless the library returns finite,
// positive kappa and lambda.
KappaSample EvaluateKappa(double massKg, double signedChargeC,
                          double momentumKgMPerS, const double positionM[3],
                          const double meanFieldT[3], double timeS,
                          std::uint64_t backgroundRevision);

}  // namespace ParallelDiffusionBinding
}  // namespace SEP

#endif  // SRCSEP_ADAPTERS_PARALLEL_DIFFUSION_ADAPTER_H
