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
#include "adapters/source_runtime.h"
#include "output/output_coordinator.h"
#include "runtime/configuration_io.h"
#include "validation/coronal_cme_application_test.h"

#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <memory>
#include <vector>

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

  // Discovery is intentionally allocation-free: CI can query the linked
  // executable without starting MPI or constructing a mesh.  Execution below
  // is different—it follows the exact production initialization path before
  // a read-only native probe observes AMPS state at a collective boundary.
  if (request.commandLine.listTests) {
    for (const auto& test : SEP3D::Validation::CoronalCmeNativeTests()) {
      std::cout << test.id << " | " << test.name << " | "
                << test.description << '\n';
    }
    return EXIT_SUCCESS;
  }

  const bool nativeTestMode = request.commandLine.allTests ||
      !request.commandLine.tests.empty();
  std::vector<SEP3D::Validation::NativeTestDescriptor> nativeTests;
  if (nativeTestMode) {
    status = SEP3D::Validation::SelectCoronalCmeNativeTests(
        request.commandLine.allTests, request.commandLine.tests,
        &nativeTests);
    if (!status.ok()) {
      std::cerr << "srcSEP3D native-test selection failed: "
                << status.message << '\n';
      return 2;
    }
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

  if (request.configuration->options().inputSchemaVersion >= 3 &&
      request.configuration->options().shock ==
          SEP3D::RuntimeModel::ShockAuthority::Swcme) {
    std::shared_ptr<SEP3D::Adapters::ShockProvider> shock;
    status = SEP3D::Adapters::CreateStandaloneSwcmeShockProvider(
        *request.configuration, &shock);
    if (status.ok()) status = SEP3D::InstallShockProvider(shock);
    if (!status.ok()) {
      std::cerr << "srcSEP3D SWCME provider initialization failed: "
                << status.message << '\n';
      return EXIT_FAILURE;
    }
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

  if (nativeTestMode) {
    // --test-steps is a small integration horizon, not a second simulation
    // loop.  It invokes the same amps_time_step() used by a production run and
    // stops immediately if the application's normal termination condition is
    // reached.  Zero is allowed for initialization-only native cases.
    for (std::uint64_t iteration = 0;
         iteration < request.commandLine.testSteps; ++iteration) {
      if (amps_time_step() == _PIC_TIMESTEP_RETURN_CODE__END_SIMULATION_) break;
    }
    MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);

    SEP3D::Validation::NativeApplicationState applicationState;
    status = SEP3D::Validation::CaptureNativeApplicationState(
        request.commandLine.expectedMpiRanks, &applicationState);
    int exitCode = status.ok() ? 0 : 2;
    if (PIC::ThisThread == 0 && status.ok()) {
      namespace fs = std::filesystem;
      std::error_code directoryError;
      fs::create_directories(request.commandLine.testArtifactDirectory,
                             directoryError);
      if (directoryError) {
        std::cerr << "srcSEP3D cannot create native-test artifact directory: "
                  << directoryError.message() << '\n';
        exitCode = 2;
      }

      std::vector<SEP3D::Validation::NativeTestResult> results;
      if (exitCode == 0) {
        results = SEP3D::Validation::EvaluateCoronalCmeNativeTests(
            applicationState, nativeTests);

        // Keep a compact human-readable state record next to the machine JSON.
        // It contains identities and counts only; large Tecplot products stay
        // at their declared paths and are referenced rather than duplicated.
        const fs::path statePath = fs::path(
            request.commandLine.testArtifactDirectory) /
            "coronal-cme-application-state.txt";
        std::ofstream stateFile(statePath);
        if (!stateFile) {
          std::cerr << "srcSEP3D cannot open native-test state artifact\n";
          exitCode = 2;
        } else {
          stateFile << "configuration_fingerprint="
                    << applicationState.configurationFingerprint << '\n'
                    << "mpi_ranks=" << applicationState.mpiRankCount << '\n'
                    << "allocated_blocks="
                    << applicationState.globalAllocatedBlocks << '\n'
                    << "physical_cells="
                    << applicationState.globalPhysicalCells << '\n'
                    << "active_region_mode=" << applicationState.activeRegionMode
                    << '\n' << "active_region_plan_installed="
                    << applicationState.activeMaskInstalled << '\n'
                    << "active_region_pruning_applied="
                    << applicationState.activeRegionPruningApplied << '\n'
                    << "active_region_allocation_verified="
                    << applicationState.activeRegionAllocationVerified << '\n'
                    << "planned_active_leaves="
                    << applicationState.plannedActiveLeaves << '\n'
                    << "planned_inactive_leaves="
                    << applicationState.plannedInactiveLeaves << '\n'
                    << "planned_solar_interior_leaves="
                    << applicationState.plannedSolarInteriorLeaves << '\n'
                    << "background_generation="
                    << applicationState.backgroundGeneration << '\n'
                    << "completed_steps="
                    << applicationState.completedSteps << '\n'
                    << "initialization_mask="
                    << applicationState.initializationMask << '\n';
          stateFile.close();
          if (!stateFile) {
            std::cerr << "srcSEP3D failed to close native-test state artifact\n";
            exitCode = 2;
          } else {
            for (auto& result : results)
              result.artifacts.push_back(statePath.string());
          }
        }

        if (exitCode == 0) {
          status = SEP3D::Validation::WriteNativeTestJson(
              request.commandLine.testJsonPath, applicationState, results);
          if (!status.ok()) {
            std::cerr << "srcSEP3D native-test JSON failed: "
                      << status.message << '\n';
            exitCode = 2;
          } else {
            exitCode = SEP3D::Validation::NativeTestExitCode(results);
            for (const auto& result : results) {
              std::cout << '[' << result.id << "] "
                        << SEP3D::Validation::Name(result.status) << " - "
                        << result.message << '\n';
            }
            std::cout << "native_test_json="
                      << request.commandLine.testJsonPath << '\n';
          }
        }
      }
    } else if (PIC::ThisThread == 0 && !status.ok()) {
      std::cerr << "srcSEP3D native state capture failed: "
                << status.message << '\n';
    }
    MPI_Bcast(&exitCode, 1, MPI_INT, 0, MPI_GLOBAL_COMMUNICATOR);
    MPI_Finalize();
    return exitCode;
  }

  // Initialization-only mode is a completed AMPS initialization, not a dry
  // parser pass: the distributed mesh has been built and decomposed, blocks,
  // particle weights, time steps, providers, and observers have been
  // initialized, and the declared Tecplot products have been closed.  Every
  // rank reaches this collective boundary before MPI is finalized; no rank can
  // enter amps_time_step() while another rank is still writing initialization
  // output.
  if (request.commandLine.initializationOnly) {
    MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);
    if (PIC::ThisThread == 0) {
      const auto& options = request.configuration->options();
      std::cout << "srcSEP3D initialization complete; no time steps executed\n"
                << "initialization_mesh="
                << options.initializationMeshTecplotFile << '\n'
                << "initialization_parker_line="
                << options.initializationParkerLineTecplotFile << '\n'
                << "initialization_data_base="
                << options.initializationDataTecplotFile << '\n';
    }
    MPI_Finalize();
    return EXIT_SUCCESS;
  }

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
