// ============================================================================
// srcSEP3D/main.cpp
//
// Standard AMPS application driver.  Through Phase R2, configuration is a
// host responsibility and the typed Runtime exists, but mesh/background/mover
// phases are not yet executable.  amps_init_mesh() therefore diagnoses the
// missing host configuration or stops at the Phase-M gate.  The driver never
// interprets argv, preserving a clean coupled-library configuration boundary.
// ============================================================================

#include "pic.h"

#include <cstdio>
#include <cstdlib>
#include <iostream>

void amps_init();
void amps_init_mesh();
int amps_time_step();

int main(int argc, char** argv) {
  (void)argc;
  (void)argv;

  amps_init_mesh();
  amps_init();

  for (long int iteration = 0; iteration < 100000001L; ++iteration) {
    if (amps_time_step() == _PIC_TIMESTEP_RETURN_CODE__END_SIMULATION_) break;
  }

  if (_PIC_NIGHTLY_TEST_MODE_ == _PIC_MODE_ON_) {
    char fileName[400];
    std::snprintf(fileName, sizeof(fileName), "%s/test_SEP3D.dat",
                  PIC::OutputDataFileDirectory);
    PIC::RunTimeSystemState::GetMeanParticleMicroscopicParameters(fileName);
  }

  MPI_Finalize();
  std::cout << "End of the run: " << PIC::nTotalSpecies << '\n';
  return EXIT_SUCCESS;
}
