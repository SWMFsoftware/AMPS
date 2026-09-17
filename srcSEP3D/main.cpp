// ============================================================================
// srcSEP3D/main.cpp
//
// Standard standalone AMPS application driver through Phases M/B/T.
//
// This executable is itself the standalone host, so it constructs a typed
// default configuration directly.  It deliberately does not parse argv or an
// AMPS parameter file.  Coupled SWMF builds call the same library entry points
// with their own immutable configuration and imported provider objects.
// ============================================================================

#include "SEP3D.h"

#include <cstdio>
#include <cstdlib>
#include <iostream>

void amps_init();
void amps_init_mesh();
int amps_time_step();

int main(int argc, char** argv) {
  (void)argc;
  (void)argv;

  SEP3D::RuntimeModel::RunConfiguration3DOptions options;
  options.shock = SEP3D::RuntimeModel::ShockAuthority::None;
  std::shared_ptr<const SEP3D::RuntimeModel::RunConfiguration3D> configuration;
  SEP3D::Core::Status status =
      SEP3D::RuntimeModel::RunConfiguration3D::Create(options, &configuration);
  if (status.ok()) status = SEP3D::ConfigureApplication(configuration);
  if (!status.ok()) {
    std::cerr << "srcSEP3D standalone configuration failed: "
              << status.message << '\n';
    return EXIT_FAILURE;
  }

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
