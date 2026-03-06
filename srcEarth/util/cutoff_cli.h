//======================================================================================
// cutoff_cli.h
//======================================================================================
// PURPOSE
//   Provide a small, robust command line interface (CLI) for running the
//   Earth energetic particle tools in "gridless" mode.
//
// USER REQUIREMENTS
//   - CLI options:
//       -h              : print help and exit
//       -mode 3d|gridless
//       -i <input-file> : AMPS_PARAM style input file
//       -mover <name>   : select particle mover (BORIS | RK2 | RK4 | RK6)
//   - The CLI and parser live in srcEarth/util.
//======================================================================================

#ifndef _SRC_EARTH_UTIL_CUTOFF_CLI_H_
#define _SRC_EARTH_UTIL_CUTOFF_CLI_H_

#include <string>

namespace EarthUtil {

  struct CliOptions {
    bool help{false};
    std::string mode{""};
    std::string inputFile{""};
    // Particle mover selection.
    // NOTE: This is intentionally a *string* here to keep the CLI independent of the
    // gridless integrator implementation. The executable can translate this string into
    // a concrete enum (MoverType) and/or an input-file setting.
    //
    // Supported values (case-insensitive):
    //   BORIS           : classic relativistic Boris pusher (legacy default)
    //   RK2/RK4/RK6   : explicit Runge-Kutta movers of order 2/4/6
    //
    // If empty, the executable should use its default / input-file setting.
    std::string mover{""};
  };

  // Parse argc/argv. Throws std::runtime_error for malformed inputs.
  CliOptions ParseCli(int argc,char** argv);

  // Return formatted help message.
  std::string HelpMessage(const char* progName);

}

#endif
