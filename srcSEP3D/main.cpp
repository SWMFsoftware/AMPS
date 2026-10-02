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
#include "output/shock_history.h"
#include "runtime/configuration_io.h"
#include "validation/coronal_cme_application_test.h"

#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <sstream>
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
                << test.description << " | suite=" << test.suite << '\n';
    }
    return EXIT_SUCCESS;
  }

  const bool nativeTestMode = request.commandLine.allTests ||
      !request.commandLine.testSuite.empty() ||
      !request.commandLine.tests.empty();
  std::vector<SEP3D::Validation::NativeTestDescriptor> nativeTests;
  if (nativeTestMode) {
    status = SEP3D::Validation::SelectCoronalCmeNativeTests(
        request.commandLine.allTests, request.commandLine.tests,
        &nativeTests, request.commandLine.testSuite);
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
                    << "input_schema_version=" << applicationState.inputSchemaVersion << '\n'
                    << "background_authority=" << applicationState.backgroundAuthority << '\n'
                    << "shock_authority=" << applicationState.shockAuthority << '\n'
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
            unsigned passed = 0, failed = 0, skipped = 0, errors = 0;
            for (const auto& result : results) {
              using SEP3D::Validation::NativeTestStatus;
              switch (result.status) {
                case NativeTestStatus::Pass: ++passed; break;
                case NativeTestStatus::Fail: ++failed; break;
                case NativeTestStatus::Skip: ++skipped; break;
                case NativeTestStatus::Error: ++errors; break;
              }
            }
            std::cout << "native_test_summary: total=" << results.size()
                      << " pass=" << passed << " fail=" << failed
                      << " skip=" << skipped << " error=" << errors << '\n';
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

  const auto& options = request.configuration->options();
  const bool propagation = options.intent ==
      SEP3D::RuntimeModel::RunIntent::ShockPropagation;
  SEP3D::Output::ShockHistoryWriter history;
  SEP3D::Output::ShockHistorySample sample;
  int propagationExit = 0;
  std::string stopReason = "step-budget";
  const std::filesystem::path nativeDirectory(options.outputDirectory);

  // Capture is collective; only rank zero owns the stream. Broadcast every
  // write result before another step can begin so an I/O failure cannot leave
  // the remaining ranks waiting in a later AMPS collective. Tick zero is an
  // observed provider state, not an independently reconstructed trajectory.
  auto publishHistory = [&]() {
    status = SEP3D::CaptureNativeShockHistorySample(&sample);
    propagationExit = status.ok() ? 0 : 2;
    if (PIC::ThisThread == 0 && status.ok()) {
      status = history.Append(sample);
      propagationExit = status.ok() ? 0 : 2;
    }
    if (PIC::ThisThread == 0 && !status.ok())
      std::cerr << "srcSEP3D propagation capture/write failed: " << status.message << '\n';
    MPI_Bcast(&propagationExit, 1, MPI_INT, 0, MPI_GLOBAL_COMMUNICATOR);
  };
  if (propagation) {
    if (PIC::ThisThread == 0) {
      // Existing telemetry is never silently overwritten, including incomplete
      // telemetry from an interrupted invocation. Use a fresh --output-dir.
      status = history.Open(nativeDirectory / "shock-history.csv", options.requestedTimeStepS);
      propagationExit = status.ok() ? 0 : 2;
      if (!status.ok()) std::cerr << "srcSEP3D propagation output failed: " << status.message << '\n';
    }
    MPI_Bcast(&propagationExit, 1, MPI_INT, 0, MPI_GLOBAL_COMMUNICATOR);
    if (!propagationExit) publishHistory();
  }
  for (std::uint64_t iteration = 0;
       !propagationExit && iteration < options.maximumTimeSteps; ++iteration) {
    const int stepCode = amps_time_step();
    if (propagation) {
      publishHistory();
      if (propagationExit) break;
      if (PIC::ThisThread == 0 && (sample.tick == 1 || sample.tick % 60 == 0))
        std::cout << "[propagation] tick=" << sample.tick << " time_s=" << sample.timeS
                  << " radius_au=" << sample.radiusM / SEP3D::Core::Const::AU
                  << " speed_m_s=" << sample.speedMPerS << " shock_active=" << sample.active
                  << " particles=" << sample.particles << std::endl;
      // Write the first row beyond the target before stopping: consumers can
      // interpolate a crossing inside the last two real native time samples.
      if (options.stopShockRadiusM > 0 && sample.radiusM >= options.stopShockRadiusM) {
        stopReason = "shock-radius";
        break;
      }
    }
    if (stepCode == _PIC_TIMESTEP_RETURN_CODE__END_SIMULATION_) {
      stopReason = "amps-termination";
      break;
    }
  }
  if (propagation && !propagationExit) {
    if (PIC::ThisThread == 0) {
      status = history.Close();
      // This is a completed-runtime record, not an observational evidence
      // manifest. The runner owns executable/input checksums, launcher logs
      // and any independently reviewed absolute launch epoch.
      if (status.ok()) {
        int ranks = 0;
        MPI_Comm_size(MPI_GLOBAL_COMMUNICATOR, &ranks);
        std::ofstream manifest(nativeDirectory / "native-runtime.json");
        manifest << std::setprecision(17)
                 << "{\n  \"schema\": \"srcsep3d-native-shock-runtime-v1\",\n"
                 << "  \"producer\": \"srcSEP3D-native\",\n"
                 << "  \"run_intent\": \"shock-propagation\",\n  \"source_enabled\": false,\n"
                 << "  \"mpi_ranks\": " << ranks << ",\n  \"time_step_s\": " << options.requestedTimeStepS
                 << ",\n  \"maximum_time_steps\": " << options.maximumTimeSteps
                 << ",\n  \"completed_steps\": " << sample.tick << ",\n  \"final_time_s\": " << sample.timeS
                 << ",\n  \"final_radius_m\": " << sample.radiusM
                 << ",\n  \"stop_shock_radius_m\": " << options.stopShockRadiusM
                 << ",\n  \"stop_reason\": \"" << stopReason << "\",\n"
                 << "  \"history_file\": \"shock-history.csv\",\n"
                 << "  \"provider_identity\": " << std::quoted(sample.providerIdentity) << ",\n"
                 << "  \"provider_configuration_fingerprint\": " << std::quoted(sample.configurationFingerprint) << ",\n"
                 // Shock history can advance more often than the mesh cadence.
                 // Record the actual installed descriptor's epoch/generation,
                 // not the final shock time, to preserve that distinction.
                 << "  \"mesh_background_authority\": " << std::quoted(SEP3D::RuntimeModel::Name(options.background)) << ",\n"
                 << "  \"mesh_background_cadence_steps\": " << options.backgroundCadenceSteps << ",\n"
                 << "  \"mesh_background_epoch_s\": " << SEP3D::ApplicationRuntime().active_snapshot()->epochS << ",\n"
                 << "  \"mesh_background_generation\": " << SEP3D::ApplicationRuntime().active_snapshot()->generation << ",\n"
                 << "  \"application_configuration_fingerprint\": " << std::quoted(request.configuration->physics_fingerprint()) << "\n}\n";
        manifest.close();
        if (!manifest) status = SEP3D::Core::Status(SEP3D::Core::StatusCode::Error, "native runtime manifest write/close failed");
      }
      propagationExit = status.ok() ? 0 : 2;
      if (!status.ok()) std::cerr << "srcSEP3D propagation close failed: " << status.message << '\n';
      else std::cout << "native_shock_history=" << (nativeDirectory / "shock-history.csv").string()
                     << "\nnative_runtime_manifest=" << (nativeDirectory / "native-runtime.json").string()
                     << "\npropagation_stop=" << stopReason << " tick=" << sample.tick << std::endl;
    }
    MPI_Bcast(&propagationExit, 1, MPI_INT, 0, MPI_GLOBAL_COMMUNICATOR);
  }
  if (propagationExit) {
    MPI_Finalize();
    return propagationExit;
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
