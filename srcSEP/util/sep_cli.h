#ifndef SEP_UTIL_SEP_CLI_H
#define SEP_UTIL_SEP_CLI_H

#include <iosfwd>

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
  bool particleCouplingMode = true;
  bool cascadeActive = true;
  bool reflectionActive = true;

  // Select how Alfvén turbulence energy is represented and evolved.
  // "integrated" is the legacy/default model with only E+ and E- per segment.
  // "wave-number-resolved" stores E+(k_j) and E-(k_j), advects each k-bin,
  // and applies particle growth/damping to the resonant k-bin.
  enum class TurbulenceModel { Integrated, WaveNumberResolved };
  TurbulenceModel turbulenceModel = TurbulenceModel::Integrated;

  // Run the standalone SEP TestManager diagnostics.  These diagnostics are
  // useful during development but can be intrusive and expensive in normal
  // production runs, so the command-line default is intentionally OFF.
  bool runTestManager = false;

  // Frequency, in main-loop iterations, for writing the large Tecplot 2-D
  // wave-number-resolved spectrum diagnostic.  The default of 100 keeps the
  // output volume manageable while still giving useful temporal resolution.
  // A value of 0 disables this diagnostic completely.  The option is ignored
  // unless the wave-number-resolved turbulence model is selected.
  int spectralOutputInterval = 100;

  bool printHelp = false;
};

// Print the command-line help for the standalone SEP driver.  This routine does
// not modify the model state; it only describes the available switches.
void PrintHelp(const char* program_name, std::ostream& out);

// Parse argc/argv into Options.  Returns false if an option is malformed or
// unknown.  In that case the caller should terminate before AMPS initialization.
bool ParseCommandLine(int argc, char** argv, Options& options,
                      std::ostream& out, std::ostream& err);

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
