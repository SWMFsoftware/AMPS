// ============================================================================
// Standalone srcSEP3D input-file and command-line boundary
//
// This is the only layer that interprets text.  It converts a versioned,
// unit-bearing schema into RunConfiguration3DOptions and then invokes the same
// immutable Create() factory used by SWMF hosts.  Physics, mesh, background,
// and transport code therefore never read argv, environment variables, or a
// stream and cannot observe a partially parsed configuration.
// ============================================================================

#ifndef SEP3D_RUNTIME_CONFIGURATION_IO_H
#define SEP3D_RUNTIME_CONFIGURATION_IO_H

#include "run_configuration.h"

#include <memory>
#include <string>
#include <vector>

namespace SEP3D {
namespace RuntimeModel {

enum class LogVerbosity { Quiet, Normal, Verbose };

struct StandaloneCommandLine {
  std::string inputPath;
  std::string outputDirectoryOverride;
  // Initialization-only execution builds the real distributed AMPS mesh,
  // completes model initialization, writes the two declared Tecplot products,
  // and exits collectively before the first particle time step.  This is
  // deliberately different from dryRun, which never allocates an AMPS mesh.
  bool initializationOnly = false;
  // Optional parent directory for the initialization mesh and Parker-line
  // products.  The input-deck filenames are retained, so this switch changes
  // only their location and cannot silently rename either declared artifact.
  std::string initializationOutputDirectory;
  std::string restartPath;
  bool dryRun = false;
  bool listTests = false;
  bool allTests = false;
  std::vector<std::string> tests;
  LogVerbosity verbosity = LogVerbosity::Normal;
};

struct StandaloneRunRequest {
  StandaloneCommandLine commandLine;
  std::shared_ptr<const RunConfiguration3D> configuration;
};

Core::Status ParseStandaloneCommandLine(
    int argc, char* const argv[], StandaloneCommandLine* result);
Core::Status ParseConfigurationText(
    const std::string& text, RunConfiguration3DOptions* result);
Core::Status LoadConfigurationFile(
    const std::string& path, RunConfiguration3DOptions* result);
// Replace only the parent directories of the two initialization Tecplot
// products.  The operation is pure: directory creation is deferred until MPI
// initialization, immediately before the collective mesh writer runs.
Core::Status ApplyInitializationOutputDirectory(
    const std::string& directory, RunConfiguration3DOptions* options);
Core::Status BuildStandaloneRunRequest(
    int argc, char* const argv[], StandaloneRunRequest* result);

// The summary is intentionally complete enough for --dry-run resource review
// but never allocates the AMPS mesh.  It includes normalized preset bounds,
// the physics fingerprint, refinement extrema, estimated level counts, and a
// whole-run memory breakdown.
Core::Status BuildDryRunSummary(const RunConfiguration3D& configuration,
                                std::string* summary);

const char* Name(LogVerbosity value);

}  // namespace RuntimeModel
}  // namespace SEP3D

#endif  // SEP3D_RUNTIME_CONFIGURATION_IO_H
