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
  };

  // Parse argc/argv. Throws std::runtime_error for malformed inputs.
  CliOptions ParseCli(int argc,char** argv);

  // Return formatted help message.
  std::string HelpMessage(const char* progName);

}

#endif
