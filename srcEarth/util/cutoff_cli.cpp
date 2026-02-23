//======================================================================================
// cutoff_cli.cpp
//======================================================================================
// IMPLEMENTATION NOTES
//   - Simple manual parsing keeps behavior consistent with the AMPS codebase.
//======================================================================================

#include "cutoff_cli.h"

#include <sstream>
#include <stdexcept>

namespace EarthUtil {

CliOptions ParseCli(int argc,char** argv) {
  CliOptions opt;
  if (argc<=1) return opt;

  for (int i=1;i<argc;i++) {
    std::string a=argv[i];
    if (a=="-h" || a=="--help") {
      opt.help=true;
    }
    else if (a=="-mode") {
      if (i+1>=argc) throw std::runtime_error("Missing value after -mode");
      opt.mode=argv[++i];
    }
    else if (a=="-i") {
      if (i+1>=argc) throw std::runtime_error("Missing value after -i");
      opt.inputFile=argv[++i];
    }
    else {
      std::ostringstream oss;
      oss << "Unknown CLI token: '" << a << "'. Use -h for help.";
      throw std::runtime_error(oss.str());
    }
  }

  return opt;
}

std::string HelpMessage(const char* progName) {
  std::ostringstream out;
  out << "AMPS Earth: Gridless Cutoff Rigidity (prototype)\n\n";
  out << "Usage:\n";
  out << "  " << progName << " -h\n";
  out << "  " << progName << " -mode gridless -i AMPS_PARAM.in\n";
  out << "  " << progName << " -mode 3d      -i AMPS_PARAM.in   (reserved / legacy path)\n\n";
  out << "Options:\n";
  out << "  -h                 Print this help and exit\n";
  out << "  -mode 3d|gridless   Select solver backend\n";
  out << "  -i <file>           Input file in AMPS_PARAM format\n\n";
  out << "Outputs (gridless):\n";
  out << "  - Summary text to stdout\n";
  out << "  - Tecplot ASCII file(s) in current directory:\n";
  out << "      cutoff_gridless_points.dat\n";
  out << "      cutoff_gridless_shells.dat\n\n";
  out << "Notes:\n";
  out << "  * Prototype implementation; coordinate transforms can be added later.\n";
  return out.str();
}

}
