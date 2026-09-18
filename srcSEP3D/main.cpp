// ============================================================================
// srcSEP3D/main.cpp
//
// Standard standalone AMPS application driver through Phases M/B/T/P/A/O.
//
// This executable is itself the standalone host.  It owns the one permitted
// text boundary: argv selects a versioned input file, configuration_io parses
// and normalizes it, and the AMPS-independent immutable factory validates the
// complete request before any AMPS mesh or MPI lifecycle operation begins.
// Coupled SWMF builds bypass this file parser and construct the same typed
// RunConfiguration3DOptions record directly.
// ============================================================================

#include "SEP3D.h"
#include "output/output_coordinator.h"
#include "runtime/configuration_io.h"

#include <cstdio>
#include <cstdlib>
#include <iostream>

void amps_init();
void amps_init_mesh();
int amps_time_step();

int main(int argc, char** argv) {
  SEP3D::RuntimeModel::StandaloneRunRequest request;
  SEP3D::Core::Status status =
      SEP3D::RuntimeModel::BuildStandaloneRunRequest(argc, argv, &request);
  if (!status.ok()) {
    std::cerr << "srcSEP3D command/configuration error: "
              << status.message << '\n';
    return 2;
  }

  // Test discovery/execution uses test/stage1 for AMPS-independent cases and
  // test/run_tests.py for linked/validation cases.  The parser recognizes the
  // common spelling so accidental use fails before AMPS starts, rather than
  // being mistaken for a production simulation.
  if (request.commandLine.listTests || request.commandLine.allTests ||
      !request.commandLine.tests.empty()) {
    std::cerr << "srcSEP3D test selection is provided by test/stage1 and "
                 "test/run_tests.py; no simulation was started\n";
    return 2;
  }

  if (request.commandLine.dryRun) {
    std::string summary;
    status = SEP3D::RuntimeModel::BuildDryRunSummary(
        *request.configuration, &summary);
    if (!status.ok()) {
      std::cerr << "srcSEP3D dry-run preflight failed: " << status.message
                << '\n';
      return EXIT_FAILURE;
    }
    std::cout << summary;
    return EXIT_SUCCESS;
  }

  status = SEP3D::ConfigureApplication(request.configuration);
  if (!status.ok()) {
    std::cerr << "srcSEP3D standalone configuration failed: "
              << status.message << '\n';
    return EXIT_FAILURE;
  }

  if (!request.configuration->options().restartInputPath.empty()) {
    SEP3D::Output::RestartLoadOptions load;
    load.expectedConfigurationFingerprint =
        request.configuration->physics_fingerprint();
    load.expectedResolvedConfigurationManifest =
        request.configuration->resolved_manifest();
    load.expectedStorageLayoutFingerprint =
        request.configuration->storage_layout().fingerprint;
    load.expectedCodeIdentity = "srcSEP3D-R01-R07";
    // Snapshot identity is read from the checkpoint.  Analytic state is
    // rebuilt at that exact epoch/generation below; coupled state must be
    // supplied by its host before initialization.
    load.missingSnapshot = SEP3D::Output::MissingSnapshotPolicy::Wait;
    load.waitTimeoutMilliseconds = 0;
    load.snapshotAvailable = [](std::uint64_t) { return true; };
    load.repartition =
        SEP3D::Output::RepartitionPolicy::DeterministicByStableId;
    SEP3D::Output::RestartState restored;
    status = SEP3D::Output::RestoreRestartBeforeMesh(
        &SEP3D::ApplicationRuntime(),
        request.configuration->options().restartInputPath, load, &restored);
    if (status.ok()) status = SEP3D::InstallRestartState(restored);
    if (!status.ok()) {
      std::cerr << "srcSEP3D restart validation failed: "
                << status.message << '\n';
      return EXIT_FAILURE;
    }
  }

  amps_init_mesh();
  amps_init();

  const std::uint64_t maximumSteps =
      request.configuration->options().maximumTimeSteps;
  for (std::uint64_t iteration = 0; iteration < maximumSteps; ++iteration) {
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
