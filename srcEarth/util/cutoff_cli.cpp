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
    else if (a=="-mover" || a=="--mover") {
      // We store the mover choice as a string. The *executable* should translate this
      // into the internal enum (MoverType) and then apply it to the shared mover
      // module (GridlessParticleMovers). This keeps this CLI file small and avoids
      // pulling gridless numerics headers into util/.
      if (i+1>=argc) throw std::runtime_error("Missing value after -mover");
      opt.mover=argv[++i];
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
  out << "  " << progName << " -mode gridless -i AMPS_PARAM.in -mover RK4\n";
  out << "  " << progName << " -mode 3d      -i AMPS_PARAM.in   (reserved / legacy path)\n\n";
  out << "Options:\n";
  out << "  -h                 Print this help and exit\n";
  out << "  -mode 3d|gridless   Select solver backend\n";
  out << "  -i <file>           Input file in AMPS_PARAM format\n\n";
  out << "  -mover <name>       Select particle mover used for trajectory integration\n";
  out << "                     Supported (case-insensitive): BORIS | RK2 | RK4 | RK6\n";
  out << "                     Default: BORIS\n";
  out << "                     NOTE: When provided, this CLI option should override any\n";
  out << "                     mover selection specified in the input file.\n\n";
  out << "Outputs (gridless):\n";
  out << "  - Summary text to stdout\n";
  out << "  - Tecplot ASCII file(s) in current directory:\n";
  out << "      cutoff_gridless_points.dat\n";
  out << "      cutoff_gridless_shells.dat\n\n";
  out << "Notes:\n";
  out << "  * Coordinate transforms (e.g., SM<->GSM) are handled in the solver code.\n";
  out << "  * This CLI only parses tokens; the executable applies them to the solver.\n";
  return out.str();
}

}
