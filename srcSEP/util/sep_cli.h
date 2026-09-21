#ifndef SEP_UTIL_SEP_CLI_H
#define SEP_UTIL_SEP_CLI_H

#include <iosfwd>
#include <string>
#include <vector>

#include "sep_production_mover.h"
#include "sep_common_header_path.h"
#include SRCSEP_SEP_COMMON_HEADER(sep_coefficient_registry.h)
#include "sep_turbulence_core.h"

namespace SEP {
namespace Util {
namespace CLI {

// Runtime command-line controls for the standalone SEP driver.
//
// These options intentionally mirror the three hard-coded turbulence switches
// that used to be set in srcSEP/main.cpp.  The defaults preserve the previous
// model behavior: particle/turbulence coupling, nonlinear cascade, and wave
// reflection are all enabled unless the user explicitly disables them on the
// command line.
struct Options {
  // Optional runtime-initialization file.  An empty path preserves the legacy
  // hard-coded domain/field-line mesh.  A supplied path is parsed and frozen
  // before AMPS initializes MPI or allocates the mesh.
  std::string inputPath;
  // Build the actual AMPS mesh and complete application initialization, write
  // the declared mesh/field-line Tecplot products, and exit collectively
  // before the first particle step.  This mode requires a versioned --input
  // deck so an apparently successful preview can never omit its products.
  bool initializationOnly = false;
  // Optional parent directory for both initialization products.  The leaf
  // names remain those reviewed in [output]; only their location is changed.
  std::string initializationOutputDirectory;
  bool particleCouplingMode = true;
  bool cascadeActive = true;
  bool reflectionActive = true;

  // Select how Alfvén turbulence energy is represented and evolved.
  // "integrated" is the legacy/default model with only E+ and E- per segment.
  // "wave-number-resolved" stores E+(k_j) and E-(k_j), advects each k-bin,
  // and applies particle growth/damping to the resonant k-bin.
  enum class TurbulenceModel { Integrated, WaveNumberResolved };
  TurbulenceModel turbulenceModel = TurbulenceModel::Integrated;

  // The production registry is deliberately limited to the three supported
  // field-line transport formulations.  The default preserves the former FTE
  // Dmumu implementation while assigning it an explicit canonical name.
  Mover::ProductionMover particleMover =
      Mover::ProductionMover::FocusedTransportDiffusion;
  bool particleMoverProvided = false;
  bool listMovers = false;

  // Step 10 exposes one validated coefficient configuration shared by Parker,
  // diffusive FTE, and event-driven FTE.  The defaults preserve the historical
  // srcSEP choices while making their names and source authority explicit.
  Transport::Coefficient::Configuration coefficients;

  // WP12 exposes the one mover-independent local error budget.  These values
  // are applied only after the full command line and compatibility matrix pass
  // validation, so malformed tolerances cannot partially mutate a run.
  Transport::NumericalTolerances numericalTolerances;

  // The pure constant Dmumu provider uses SI s^-1.  It is stored beside the
  // registry selection so parse-only tests can verify CLI-to-provider intent
  // without linking PIC globals.
  double constantDmumuPerS = 0.0;
  bool constantDmumuProvided = false;

  // Step 11 keeps turbulence-source ownership separate from the exactly three
  // particle movers.  This configuration is the authoritative CLI contract for
  // source, representation, boundary, operator, spectral, cadence, and
  // conservation settings; legacy booleans above remain compatibility aliases.
  Turbulence::Configuration turbulence;

  // Run the standalone SEP TestManager diagnostics.  These diagnostics are
  // useful during development but can be intrusive and expensive in normal
  // production runs, so the command-line default is intentionally OFF.
  bool runTestManager = false;

  // Step 1 adds a test-only run mode to the existing parser rather than
  // introducing a second command-line implementation.  Selectors are retained
  // in parse order here; the registry later resolves them case-insensitively,
  // de-duplicates overlaps, and sorts by stable test ID.
  bool listTests = false;
  bool runAllTests = false;
  std::vector<std::string> testIds;
  std::vector<std::string> testGroups;

  // Structured reports are opt-in and valid only for an executing component
  // test selection.  Empty paths disable the corresponding writer; report
  // creation failure is a test ERROR rather than a warning.
  std::string testJsonPath;
  std::string testJunitPath;

  // End-to-end validation cases can require a reviewed model-input manifest
  // and an isolated artifact directory. These paths are carried by the same
  // production parser as every other test option. They do not affect ordinary
  // component tests and are legal only with exactly one explicit --test ID.
  std::string testInputPath;
  std::string testArtifactDirectory;

  // Frequency, in main-loop iterations, for writing the large Tecplot 2-D
  // wave-number-resolved spectrum diagnostic.  The default of 100 keeps the
  // output volume manageable while still giving useful temporal resolution.
  // A value of 0 disables this diagnostic completely.  The option is ignored
  // unless the wave-number-resolved turbulence model is selected.
  int spectralOutputInterval = 100;

  // Number of SEP macroparticles injected on each active field line injection
  // event.  This maps directly to
  // SEP::FieldLine::InjectionParameters::nParticlesPerIteration, whose legacy
  // hard-coded default is 300 in field_line.cpp.  The option is exposed through
  // the CLI because convergence/noise tests often require changing the number
  // of injected macroparticles without recompiling.  The value must be positive:
  // zero injected particles would make the injection weight correction formulas
  // singular and would silently disable the SEP source.
  int injectionParticlesPerIteration = 300;
  // Preserve command-line precedence over the versioned initialization file.
  // Without this provenance bit, the historical default 300 is
  // indistinguishable from an operator's explicit override.
  bool injectionParticlesProvided = false;

  // WP30 driver controls.  Each "Provided" bit preserves provenance when the
  // value remains a default and lets the final RunConfiguration enforce the
  // documented defaults < input file < command line precedence.
  int totalIterations = 100000001;
  bool totalIterationsProvided = false;
  double fieldLineSeedAreaM2 = 3.14159265358979323846;
  bool fieldLineSeedAreaProvided = false;
  double shockTurbulenceEfficiency = 0.02;
  bool shockTurbulenceEfficiencyProvided = false;
  double shockTurbulencePlusFraction = 0.5;
  bool shockTurbulencePlusFractionProvided = false;
  int mergeMinimum = 600;
  int mergeMaximum = 1000;
  bool mergeMinimumProvided = false;
  bool mergeMaximumProvided = false;
  bool analyticalShock = false;
  bool shockModelProvided = false;
  bool slowCmeScenario = false;
  bool cmeScenarioProvided = false;

  // D01 background-failure controls. Strict is intentionally the default.
  // The fallback values are SI and have no effect unless the user explicitly
  // chooses diagnostic-fallback; all choices are copied into the frozen run
  // configuration and consequently into restart/provenance fingerprints.
  enum class SwcmeFailurePolicy { Strict, ClampRadius, DiagnosticFallback };
  SwcmeFailurePolicy swcmeFailurePolicy = SwcmeFailurePolicy::Strict;
  bool swcmeFailurePolicyProvided = false;
  double swcmeFallbackDensityM3 = 5.0e6;
  double swcmeFallbackSpeedMPerS = 4.0e5;
  double swcmeFallbackDivergencePerS = 0.0;
  bool swcmeFallbackDensityProvided = false;
  bool swcmeFallbackSpeedProvided = false;
  bool swcmeFallbackDivergenceProvided = false;
  // Repeatable canonical D02 overrides. The CLI checks only key=value shape;
  // the canonical SWCME resolver owns key names, units, and physical ranges.
  std::vector<std::string> swcmeOverrides;

  bool printHelp = false;
};

// Print the command-line help for the standalone SEP driver.  This routine does
// not modify the model state; it only describes the available switches.
void PrintHelp(const char* program_name, std::ostream& out);

// Parse argc/argv into Options.  Returns false if an option is malformed or
// unknown.  In that case the caller should terminate before AMPS initialization.
bool ParseCommandLine(int argc, char** argv, Options& options,
                      std::ostream& out, std::ostream& err);

// New registry selectors imply a test-only process: the selected tests run and
// the executable exits before the production timestep loop.  The historical
// --test-manager switch is intentionally excluded because its frozen behavior
// runs diagnostics and then continues into production.
bool IsComponentTestExecutionRequested(const Options& options);

// Apply the parsed options to the turbulence model flags used by the physics
// kernels.  Keeping this in a separate function makes the point where CLI
// options modify the model state explicit in main.cpp.
void ApplyTurbulenceOptions(const Options& options);

// Print a compact summary of the active turbulence switches.  This is useful in
// batch logs because it records the run-time configuration independently of the
// input file and compile-time macros.
void PrintTurbulenceOptions(const Options& options, std::ostream& out);

} // namespace CLI
} // namespace Util
} // namespace SEP

#endif // SEP_UTIL_SEP_CLI_H
