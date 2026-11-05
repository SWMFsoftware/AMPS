#ifndef _COMET_LOCAL_TIME_STEP_VALIDATION_H_
#define _COMET_LOCAL_TIME_STEP_VALIDATION_H_

namespace Comet {
namespace LocalTimeStepValidation {

struct cOptions {
  long int ParticlesPerCase=200000;
  bool TestFirstOrder=true;
  bool TestSecondOrder=true;
  bool TestMultiFace=true;
  // Correction-on Stage-4 builds normally exercise the wide 1/16...16
  // pattern.  This option forces every injected particle to c=1 through the
  // same correction-enabled binary so the correction-on identity limit can
  // be compared with the correction-off reference path.
  bool UseUnitCorrections=false;
};

// Run controlled transport-only tests on the already initialized DIM1 comet
// mesh.  Single-face cases vary timestep and (when compiled) block-local
// particle-weight ratios independently and audit conservation of
// W*c/dt-carried number, mass, momentum, kinetic energy, and internal energy.
// Correction-on builds inject a deterministic wide cycle of individual
// corrections; exact-bin histograms and physical-moment checks verify that
// local cloning, MPI packing, and continuation preserve c.  Correction-off
// builds use c=1 and retain the same test path.  No source, collision,
// chemistry, cooling, or sampling operation is called.  The return value is
// collective: every MPI rank receives the same pass/fail result.
bool Run(const cOptions &options);

}
}

#endif
