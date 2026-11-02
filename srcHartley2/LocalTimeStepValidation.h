#ifndef _COMET_LOCAL_TIME_STEP_VALIDATION_H_
#define _COMET_LOCAL_TIME_STEP_VALIDATION_H_

namespace Comet {
namespace LocalTimeStepValidation {

struct cOptions {
  long int ParticlesPerCase=200000;
  bool TestFirstOrder=true;
  bool TestSecondOrder=true;
  bool TestMultiFace=true;
};

// Run controlled transport-only tests on the already initialized DIM1 comet
// mesh.  Single-face cases vary timestep and (when compiled) block-local
// particle-weight ratios independently and audit conservation of W/dt-carried
// number, mass, momentum, and energy flux.  No source, collision, chemistry,
// cooling, or sampling operation is called.  The return value is collective:
// every MPI rank receives the same pass/fail result.
bool Run(const cOptions &options);

}
}

#endif
