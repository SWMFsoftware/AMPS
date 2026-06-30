//======================================================================================
// cutoff_cli.cpp
//======================================================================================
//
// See cutoff_cli.h for the full CLI reference (all flags, semantics, usage examples,
// and error-handling contract).
//
//======================================================================================
// IMPLEMENTATION NOTES
//======================================================================================
//
// PARSING STRATEGY
// ----------------
// Manual argv scan rather than getopt/boost::program_options. This keeps the
// CLI independent of third-party libraries (consistent with the rest of the
// gridless utility build path, which links only against the C++ standard library
// and the AMPS framework).
//
// The parser walks argv[1..argc-1] with an index cursor i:
//   for i in 1..argc-1:
//     if argv[i] is a recognised flag expecting a value:
//       if i+1 >= argc: exit(__LINE__,__FILE__,"flag needs argument")
//       store argv[++i] as the value
//     elif argv[i] is a recognised boolean flag (e.g., -h):
//       set the corresponding bool field
//     else:
//       exit(__LINE__,__FILE__,"unrecognised option: <argv[i]>")
//
// Values are stored as strings without case normalisation; the solver
// normalises internally (ToUpper, ParseMoverType, etc.).
//
// HELP TEXT FORMAT
// ----------------
// HelpMessage returns a formatted string that the caller prints to stdout (not stderr).
// The format is:
//
//   Usage: <progName> [options]
//
//   Options:
//     -h | --help              Print this message and exit.
//     -mode <3d|gridless>      Select solver mode.
//     -i <file>                Input parameter file (AMPS_PARAM format).
//     -mover <BORIS|RK2|RK4|RK6|GC2|GC4|GC6|HYBRID>   Particle mover (default: BORIS).
//     -mode3d-output-initialized                         Write amps_3d_initialized.data.dat.
//     -mode3d-field-eval <INTERPOLATION|ANALYTIC>        3D B-field source during tracing.
//     -density-mode <ISOTROPIC|ANISOTROPIC>
//                              Override DS_BOUNDARY_MODE from input file.
//     -density-parallel <OPENMP|THREADS|SERIAL>          Mode3D shared-memory backend.
//     -density-threads <int>                             threads per MPI process.
//     -mode3d-mpi-scheduler <DYNAMIC|BLOCK_CYCLIC|STATIC> MPI-rank scheduler.
//       Gridless aliases: -gridless-mpi-scheduler, -gridless-mpi-dynamic-chunk.
//     -mode3d-mpi-dynamic-chunk <int>                    locations per MPI fetch.
//     -cutoff-search <UPPER_SCAN|BINARY>                 cutoff search in 3d/gridless.
//     -cutoff-upper-scan-n <int>                         log-grid samples for UPPER_SCAN.
//
//   See amps_param_parser.h for a full description of the input file format.
//
// The progName argument is typically argv[0] (the executable path as invoked).
//
//======================================================================================

#include "cutoff_cli.h"
#include "specfunc.h"

#include <sstream>

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
      if (i+1>=argc) exit(__LINE__,__FILE__,"Missing value after -mode");
      opt.mode=argv[++i];
    }
    else if (a=="-i") {
      if (i+1>=argc) exit(__LINE__,__FILE__,"Missing value after -i");
      opt.inputFile=argv[++i];
    }
    else if (a=="-mover" || a=="--mover") {
      // We store the mover choice as a string. The *executable* should translate this
      // into the internal enum (MoverType) and then apply it to the shared mover
      // module (GridlessParticleMovers). This keeps this CLI file small and avoids
      // pulling gridless numerics headers into util/.
      if (i+1>=argc) exit(__LINE__,__FILE__,"Missing value after -mover");
      opt.mover=argv[++i];
    }
    else if (a=="-density-mode" || a=="--density-mode") {
      // Overrides DS_BOUNDARY_MODE from the input file.
      // Supported values (case-insensitive): ISOTROPIC | ANISOTROPIC
      if (i+1>=argc) exit(__LINE__,__FILE__,"Missing value after -density-mode");
      opt.densityMode=argv[++i];
    }
    else if (a=="-density-transmission-mode" || a=="--density-transmission-mode" ||
             a=="-density-transmission" || a=="--density-transmission" ||
             a=="-ds-transmission-mode" || a=="--ds-transmission-mode") {
      // Select how the density/flux solver samples the transmission function T(E,Omega).
      // DIRECT keeps the legacy energy grid.  SCAN and ADAPTIVE use a log-spaced
      // rigidity grid, which is more appropriate near cutoffs and penumbra.
      if (i+1>=argc) exit(__LINE__,__FILE__,"Missing value after -density-transmission-mode");
      opt.densityTransmissionMode=argv[++i];
    }
    else if (a=="-density-transmission-scan-n" || a=="--density-transmission-scan-n" ||
             a=="-density-transmission-n" || a=="--density-transmission-n" ||
             a=="-ds-transmission-scan-n" || a=="--ds-transmission-scan-n") {
      if (i+1>=argc) exit(__LINE__,__FILE__,"Missing value after -density-transmission-scan-n");
      opt.densityTransmissionScanN=std::stoi(argv[++i]);
      if (opt.densityTransmissionScanN < 0)
        exit(__LINE__,__FILE__,"-density-transmission-scan-n must be >= 0");
    }
    else if (a=="-density-transmission-refine-n" || a=="--density-transmission-refine-n" ||
             a=="-ds-transmission-refine-n" || a=="--ds-transmission-refine-n") {
      if (i+1>=argc) exit(__LINE__,__FILE__,"Missing value after -density-transmission-refine-n");
      opt.densityTransmissionRefineN=std::stoi(argv[++i]);
      if (opt.densityTransmissionRefineN < 0)
        exit(__LINE__,__FILE__,"-density-transmission-refine-n must be >= 0");
    }
    else if (a=="-density-transmission-max-n" || a=="--density-transmission-max-n" ||
             a=="-ds-transmission-max-n" || a=="--ds-transmission-max-n") {
      if (i+1>=argc) exit(__LINE__,__FILE__,"Missing value after -density-transmission-max-n");
      opt.densityTransmissionMaxN=std::stoi(argv[++i]);
      if (opt.densityTransmissionMaxN < 0)
        exit(__LINE__,__FILE__,"-density-transmission-max-n must be >= 0");
    }
    else if (a=="-density-parallel" || a=="--density-parallel" ||
             a=="-mode3d-density-parallel" || a=="--mode3d-density-parallel" ||
             a=="-mode3d-density-backend" || a=="--mode3d-density-backend" ||
             a=="-mode3d-parallel" || a=="--mode3d-parallel" ||
             a=="-gridless-parallel" || a=="--gridless-parallel" ||
             a=="-gridless-backend" || a=="--gridless-backend" ||
             a=="-mode3d-backend" || a=="--mode3d-backend" ||
             a=="-backtrack-parallel" || a=="--backtrack-parallel" ||
             a=="-backtrack-backend" || a=="--backtrack-backend") {
      // Select shared-memory backend for Mode3D density backtracking within each
      // MPI process.  Validated in main after parsing.
      if (i+1>=argc) exit(__LINE__,__FILE__,"Missing value after -density-parallel");
      opt.densityParallelBackend=argv[++i];
    }
    else if (a=="-density-threads" || a=="--density-threads" ||
             a=="-mode3d-density-threads" || a=="--mode3d-density-threads" ||
             a=="-mode3d-threads" || a=="--mode3d-threads" ||
             a=="-gridless-threads" || a=="--gridless-threads" ||
             a=="-backtrack-threads" || a=="--backtrack-threads" ||
             a=="-n-density-threads" || a=="--n-density-threads") {
      // Number of shared-memory workers per MPI rank for Mode3D density
      // backtracking.  For OPENMP this maps to omp_set_num_threads(N); for
      // THREADS it is the number of std::thread workers.  0 means automatic.
      if (i+1>=argc) exit(__LINE__,__FILE__,"Missing value after -density-threads");
      opt.densityThreads=std::stoi(argv[++i]);
      if (opt.densityThreads < 0)
        exit(__LINE__,__FILE__,"-density-threads must be >= 0 (0 means automatic)");
    }
    else if (a=="-mode3d-mpi-scheduler" || a=="--mode3d-mpi-scheduler" ||
             a=="-mode3d-mpi-backend" || a=="--mode3d-mpi-backend" ||
             a=="-gridless-mpi-scheduler" || a=="--gridless-mpi-scheduler" ||
             a=="-gridless-mpi-backend" || a=="--gridless-mpi-backend" ||
             a=="-backtrack-mpi-scheduler" || a=="--backtrack-mpi-scheduler" ||
             a=="-backtrack-mpi-backend" || a=="--backtrack-mpi-backend") {
      // Inter-rank scheduler for standalone Mode3D backward products.
      // DYNAMIC activates the MPI one-sided atomic work queue; BLOCK_CYCLIC and
      // STATIC are deterministic fallback/debug schedulers.
      if (i+1>=argc) exit(__LINE__,__FILE__,"Missing value after -mode3d-mpi-scheduler");
      opt.mode3dMpiScheduler=argv[++i];
    }
    else if (a=="-mode3d-mpi-dynamic-chunk" || a=="--mode3d-mpi-dynamic-chunk" ||
             a=="-mode3d-mpi-chunk" || a=="--mode3d-mpi-chunk" ||
             a=="-gridless-mpi-dynamic-chunk" || a=="--gridless-mpi-dynamic-chunk" ||
             a=="-gridless-mpi-chunk" || a=="--gridless-mpi-chunk" ||
             a=="-backtrack-mpi-dynamic-chunk" || a=="--backtrack-mpi-dynamic-chunk" ||
             a=="-backtrack-mpi-chunk" || a=="--backtrack-mpi-chunk") {
      // Number of global observation locations fetched per dynamic MPI request.
      // 0 means automatic; negative values are invalid.
      if (i+1>=argc) exit(__LINE__,__FILE__,"Missing value after -mode3d-mpi-dynamic-chunk");
      opt.mode3dMpiDynamicChunk=std::stoi(argv[++i]);
      if (opt.mode3dMpiDynamicChunk < 0)
        exit(__LINE__,__FILE__,"-mode3d-mpi-dynamic-chunk must be >= 0 (0 means automatic)");
    }
    else if (a=="-cutoff-debug-scan" || a=="--cutoff-debug-scan" ||
             a=="-mode3d-cutoff-debug-scan" || a=="--mode3d-cutoff-debug-scan") {
      // Convenience CLI for the Mode3D cutoff rigidity classification diagnostic.
      // The input-file equivalent is CUTOFF_DEBUG_RIGIDITY_SCAN=T plus the
      // CUTOFF_DEBUG_SCAN_LON/LAT/ALT keywords in #CUTOFF_RIGIDITY.
      if (i+3>=argc)
        exit(__LINE__,__FILE__,"Missing values after -cutoff-debug-scan; expected lon_deg lat_deg alt_km");
      opt.cutoffDebugScan=true;
      opt.cutoffDebugScanLon_deg=std::stod(argv[++i]);
      opt.cutoffDebugScanLat_deg=std::stod(argv[++i]);
      opt.cutoffDebugScanAlt_km=std::stod(argv[++i]);
      if (opt.cutoffDebugScanLat_deg < -90.0 || opt.cutoffDebugScanLat_deg > 90.0)
        exit(__LINE__,__FILE__,"-cutoff-debug-scan latitude must be in [-90,90] degrees");
    }
    else if (a=="-cutoff-debug-scan-n" || a=="--cutoff-debug-scan-n") {
      if (i+1>=argc) exit(__LINE__,__FILE__,"Missing value after -cutoff-debug-scan-n");
      opt.cutoffDebugScanN=std::stoi(argv[++i]);
      if (opt.cutoffDebugScanN < 2)
        exit(__LINE__,__FILE__,"-cutoff-debug-scan-n must be >= 2");
    }
    else if (a=="-cutoff-debug-scan-file" || a=="--cutoff-debug-scan-file") {
      if (i+1>=argc) exit(__LINE__,__FILE__,"Missing value after -cutoff-debug-scan-file");
      opt.cutoffDebugScanFile=argv[++i];
      if (opt.cutoffDebugScanFile.empty())
        exit(__LINE__,__FILE__,"-cutoff-debug-scan-file must not be empty");
    }
    else if (a=="-cutoff-debug-exit" || a=="--cutoff-debug-exit" ||
             a=="-mode3d-cutoff-debug-exit" || a=="--mode3d-cutoff-debug-exit") {
      // Detailed exit-classifier diagnostic for one vertical trajectory point.
      // Use -cutoff-debug-exit-r to trace one specific rigidity, or omit it to
      // write a small rigidity list similar to -cutoff-debug-scan.
      if (i+3>=argc)
        exit(__LINE__,__FILE__,"Missing values after -cutoff-debug-exit; expected lon_deg lat_deg alt_km");
      opt.cutoffDebugExit=true;
      opt.cutoffDebugExitLon_deg=std::stod(argv[++i]);
      opt.cutoffDebugExitLat_deg=std::stod(argv[++i]);
      opt.cutoffDebugExitAlt_km=std::stod(argv[++i]);
      if (opt.cutoffDebugExitLat_deg < -90.0 || opt.cutoffDebugExitLat_deg > 90.0)
        exit(__LINE__,__FILE__,"-cutoff-debug-exit latitude must be in [-90,90] degrees");
    }
    else if (a=="-cutoff-debug-exit-r" || a=="--cutoff-debug-exit-r" ||
             a=="-cutoff-debug-exit-r-gv" || a=="--cutoff-debug-exit-r-gv") {
      if (i+1>=argc) exit(__LINE__,__FILE__,"Missing value after -cutoff-debug-exit-r");
      opt.cutoffDebugExitR_GV=std::stod(argv[++i]);
    }
    else if (a=="-cutoff-debug-exit-n" || a=="--cutoff-debug-exit-n") {
      if (i+1>=argc) exit(__LINE__,__FILE__,"Missing value after -cutoff-debug-exit-n");
      opt.cutoffDebugExitN=std::stoi(argv[++i]);
      if (opt.cutoffDebugExitN < 2)
        exit(__LINE__,__FILE__,"-cutoff-debug-exit-n must be >= 2");
    }
    else if (a=="-cutoff-debug-exit-file" || a=="--cutoff-debug-exit-file") {
      if (i+1>=argc) exit(__LINE__,__FILE__,"Missing value after -cutoff-debug-exit-file");
      opt.cutoffDebugExitFile=argv[++i];
      if (opt.cutoffDebugExitFile.empty())
        exit(__LINE__,__FILE__,"-cutoff-debug-exit-file must not be empty");
    }
    else if (a=="-cutoff-search" || a=="--cutoff-search" ||
             a=="-cutoff-search-algorithm" || a=="--cutoff-search-algorithm" ||
             a=="-mode3d-cutoff-search" || a=="--mode3d-cutoff-search" ||
             a=="-gridless-cutoff-search" || a=="--gridless-cutoff-search" ||
             a=="-backtrack-cutoff-search" || a=="--backtrack-cutoff-search") {
      // Select the scalar cutoff-rigidity search algorithm used by both standalone
      // Mode3D and gridless cutoff.  The option is kept generic because the actual
      // solver mode is chosen later by -mode; here we only store the requested token.
      if (i+1>=argc) exit(__LINE__,__FILE__,"Missing value after -cutoff-search");
      opt.cutoffSearchAlgorithm=argv[++i];
    }
    else if (a=="-cutoff-upper-scan-n" || a=="--cutoff-upper-scan-n" ||
             a=="-cutoff-search-n" || a=="--cutoff-search-n" ||
             a=="-mode3d-cutoff-search-n" || a=="--mode3d-cutoff-search-n" ||
             a=="-gridless-cutoff-search-n" || a=="--gridless-cutoff-search-n" ||
             a=="-gridless-cutoff-upper-scan-n" || a=="--gridless-cutoff-upper-scan-n" ||
             a=="-backtrack-cutoff-search-n" || a=="--backtrack-cutoff-search-n") {
      // Number of log-spaced samples used by the UPPER_SCAN pre-scan before the
      // final local bisection.  Applies to both -mode 3d and -mode gridless.
      if (i+1>=argc) exit(__LINE__,__FILE__,"Missing value after -cutoff-upper-scan-n");
      opt.cutoffUpperScanN=std::stoi(argv[++i]);
      if (opt.cutoffUpperScanN < 2)
        exit(__LINE__,__FILE__,"-cutoff-upper-scan-n must be >= 2");
    }
    else if (a=="-mode3d-output-initialized" || a=="--mode3d-output-initialized" ||
             a=="-output-3d-initialized" || a=="--output-3d-initialized") {
      // Optional Mode3D diagnostic output. Kept as a boolean flag because the
      // default behavior should be no output unless explicitly requested.
      opt.mode3dOutputInitialized=true;
    }
    else if (a=="-mode3d-field-eval" || a=="--mode3d-field-eval" ||
             a=="-mode3d-b-field" || a=="--mode3d-b-field") {
      // Optional Mode3D magnetic-field source during tracing.
      // Validated in main after the input file has been parsed.
      if (i+1>=argc) exit(__LINE__,__FILE__,"Missing value after -mode3d-field-eval");
      opt.mode3dFieldEval=argv[++i];
    }
    else if (a=="-mode3d-analytic-field" || a=="--mode3d-analytic-field") {
      // Convenience alias for the common diagnostic/comparison use case.
      opt.mode3dFieldEval="ANALYTIC";
    }
    else if (a=="-max-trace-distance" || a=="--max-trace-distance") {
      // CLI override for the global cumulative trace-distance cap.
      // Units: Earth radii (Re). The option mirrors the input-file key
      // #NUMERICAL / MAX_TRACE_DISTANCE.
      if (i+1>=argc) exit(__LINE__,__FILE__,"Missing value after -max-trace-distance");
      opt.maxTraceDistance_Re=std::stod(argv[++i]);
      if (opt.maxTraceDistance_Re < 0.0) {
        exit(__LINE__,__FILE__,"-max-trace-distance must be >= 0 (0 means: disable the cap)");
      }
    }
    // ----- 3d_forward specific options -----
    else if (a=="-forward-niter" || a=="--forward-niter") {
      // Override Mode3DForwardOptions::nIterations.
      if (i+1>=argc) exit(__LINE__,__FILE__,"Missing value after -forward-niter");
      opt.forward3dNiter = std::stoi(argv[++i]);
      if (opt.forward3dNiter < 1)
        exit(__LINE__,__FILE__,"-forward-niter must be >= 1");
    }
    else if (a=="-forward-nparticles" || a=="--forward-nparticles") {
      // Override Mode3DForwardOptions::nParticlesPerIter.
      // Sets the number of simulation particles injected at the domain boundary per
      // iteration.  The physical particle weight W is then automatically derived from
      // the boundary spectrum integral and total boundary area:
      //   W = (pi x int_{Emin}^{Emax} J(E) dE x A_boundary x dt) / nParticlesPerIter
      // Increasing this improves statistical sampling of the density distribution at
      // the cost of proportionally higher memory and compute per iteration.
      if (i+1>=argc) exit(__LINE__,__FILE__,"Missing value after -forward-nparticles");
      opt.forward3dNparticles = std::stoi(argv[++i]);
      if (opt.forward3dNparticles < 1)
        exit(__LINE__,__FILE__,"-forward-nparticles must be >= 1");
    }
    else if (a=="-forward-boundary-dist" || a=="--forward-boundary-dist") {
      // Override Mode3DForwardOptions::boundaryDistType.
      // Accepted values (case-insensitive): ISOTROPIC (others reserved for future use).
      if (i+1>=argc) exit(__LINE__,__FILE__,"Missing value after -forward-boundary-dist");
      opt.forward3dBoundaryDist = argv[++i];
    }
    else if (a=="-forward-injection-energy" || a=="--forward-injection-energy" ||
             a=="-forward-energy-sampling" || a=="--forward-energy-sampling") {
      // Select the proposal distribution for injected particle kinetic energy.
      // SPECTRUM preserves the legacy branch: sample directly from J(E)dE.
      // LOG_UNIFORM samples uniformly in log(E) and requires individual particle
      // statistical-weight corrections in Mode3DForward.cpp.
      if (i+1>=argc) exit(__LINE__,__FILE__,"Missing value after -forward-injection-energy");
      opt.forward3dInjectionEnergyDistribution = argv[++i];
    }
    else if (a=="-forward-injection-emin" || a=="--forward-injection-emin" ||
             a=="-forward-emin" || a=="--forward-emin" ||
             a=="-forward-injection-energy-min" || a=="--forward-injection-energy-min") {
      // Lower kinetic-energy limit [MeV/n] for 3d_forward injected particles.
      // This overrides the effective #DENSITY_3D DENS_EMIN value so the source
      // range and density-output energy grid remain synchronized.
      if (i+1>=argc) exit(__LINE__,__FILE__,"Missing value after -forward-injection-emin");
      opt.forward3dInjectionEmin_MeV = std::stod(argv[++i]);
      if (!(opt.forward3dInjectionEmin_MeV > 0.0))
        exit(__LINE__,__FILE__,"-forward-injection-emin must be > 0 (MeV/n)");
    }
    else if (a=="-forward-injection-emax" || a=="--forward-injection-emax" ||
             a=="-forward-emax" || a=="--forward-emax" ||
             a=="-forward-injection-energy-max" || a=="--forward-injection-energy-max") {
      // Upper kinetic-energy limit [MeV/n] for 3d_forward injected particles.
      // This overrides the effective #DENSITY_3D DENS_EMAX value so the source
      // range and density-output energy grid remain synchronized.
      if (i+1>=argc) exit(__LINE__,__FILE__,"Missing value after -forward-injection-emax");
      opt.forward3dInjectionEmax_MeV = std::stod(argv[++i]);
      if (!(opt.forward3dInjectionEmax_MeV > 0.0))
        exit(__LINE__,__FILE__,"-forward-injection-emax must be > 0 (MeV/n)");
    }
    else if (a=="-forward-injection-energy-range" || a=="--forward-injection-energy-range" ||
             a=="-forward-energy-range" || a=="--forward-energy-range") {
      // Convenience form: set both limits in one option.
      // Example: -forward-injection-energy-range 1.0 20000.0
      if (i+2>=argc) exit(__LINE__,__FILE__,"Missing values after -forward-injection-energy-range");
      opt.forward3dInjectionEmin_MeV = std::stod(argv[++i]);
      opt.forward3dInjectionEmax_MeV = std::stod(argv[++i]);
      if (!(opt.forward3dInjectionEmin_MeV > 0.0))
        exit(__LINE__,__FILE__,"-forward-injection-energy-range Emin must be > 0 (MeV/n)");
      if (!(opt.forward3dInjectionEmax_MeV > opt.forward3dInjectionEmin_MeV))
        exit(__LINE__,__FILE__,"-forward-injection-energy-range requires Emax > Emin (MeV/n)");
    }
    else {
      std::ostringstream oss;
      oss << "Unknown CLI token: '" << a << "'. Use -h for help.";
      exit(__LINE__,__FILE__,oss.str().c_str());
    }
  }

  if (opt.forward3dInjectionEmin_MeV > 0.0 &&
      opt.forward3dInjectionEmax_MeV > 0.0 &&
      !(opt.forward3dInjectionEmax_MeV > opt.forward3dInjectionEmin_MeV)) {
    exit(__LINE__, __FILE__,
         "3d_forward injection energy limits require Emax > Emin (MeV/n)");
  }

  return opt;
}

std::string HelpMessage(const char* progName) {
  std::ostringstream out;
  out << "AMPS Earth: Gridless Energetic Particle Solver\n\n";

  // -----------------------------------------------------------------------
  out << "Usage:\n";
  out << "  " << progName << " -h\n";
  out << "  " << progName << " -mode gridless -i AMPS_PARAM.in [options]\n";
  out << "  " << progName << " -mode 3d       -i AMPS_PARAM.in [options]\n";
  out << "  " << progName << " -mode 3d_forward -i AMPS_PARAM.in [options]\n";
  out << "\n";

  // -----------------------------------------------------------------------
  out << "Options:\n\n";

  out << "  -h | --help\n";
  out << "      Print this message and exit. No solver is run.\n\n";

  out << "  -mode <3d|gridless>   (required)\n";
  out << "      Select the solver execution path.\n";
  out << "        3d        Full PIC-backed 3D solver with mesh-interpolated fields.\n";
  out << "        gridless  Direct Tsyganenko/IGRF field evaluation + backtracing.\n";
  out << "                  Supports both cutoff-rigidity and density/spectrum modes.\n";
  out << "        3d_forward  Full PIC-backed 3D forward particle transport.\n";
  out << "                  Injects particles at the outer domain boundary, propagates\n";
  out << "                  them forward in time under the configured B/E fields, and\n";
  out << "                  accumulates 3D volumetric number-density on the AMR mesh.\n";
  out << "                  Requires #DENSITY_3D section in the input file.\n\n";

  out << "  -i <file>   (required)\n";
  out << "      Path to the AMPS_PARAM-format input file. Controls all physics and\n";
  out << "      geometry parameters. See 'Input file sections' below for a summary.\n\n";

  out << "  -mover | --mover <BORIS|RK2|RK4|RK6|GC2|GC4|GC6|HYBRID>   (optional; default: BORIS)\n";
  out << "      Select the particle integration algorithm.\n";
  out << "        BORIS   Relativistic Boris pusher. Volume-preserving (symplectic),\n";
  out << "                time-reversible, one field evaluation per step. Recommended\n";
  out << "                for all production runs. Validated against Stormer analytic\n";
  out << "                cutoffs on dipole fields.\n";
  out << "        RK2     Runge-Kutta 2nd order (Heun). Two field evaluations per step.\n";
  out << "                Useful as a quick sanity check against Boris.\n";
  out << "        RK4     Runge-Kutta 4th order. Four field evaluations per step.\n";
  out << "                Use when orbit reconstruction accuracy matters more than speed.\n";
  out << "        RK6     Runge-Kutta 6th order. Six evaluations per step. Recommended\n";
  out << "                near the inner boundary in strong storm-time T05 fields.\n";
  out << "        GC2/4/6 Guiding-center equations integrated with RK2/RK4/RK6.\n";
  out << "                Useful when fast gyromotion is not the quantity of interest.\n";
  out << "        HYBRID  Per-step switch between RK4 and GC4. The code evaluates a\n";
  out << "                local adiabaticity estimate rho/L_eff: if the gyroradius is\n";
  out << "                small compared with field-gradient/curvature scales it uses\n";
  out << "                GC4, otherwise it falls back to full-orbit RK4. Intended to\n";
  out << "                support a larger but still physically correct dt.\n";
  out << "      In -mode 3d_forward, the AMPS-signature mover dispatcher supports\n";
  out << "      BORIS, RK4, GC/GC4, and HYBRID. RK2/RK6/GC2/GC6 remain available in\n";
  out << "      gridless/backward 3D tracing through the shared gridless mover layer.\n";
  out << "      When provided, overrides any mover setting in the input file.\n\n";

  out << "  -density-mode | --density-mode <ISOTROPIC|ANISOTROPIC>   (optional)\n";
  out << "      Override DS_BOUNDARY_MODE from the input file.\n";
  out << "      Only applies when CALC_TARGET = DENSITY_SPECTRUM.\n\n";
  out << "        ISOTROPIC   (default if flag is absent)\n";
  out << "            Boundary differential intensity J_b(E) is assumed uniform\n";
  out << "            over all arrival directions and all positions on the outer domain\n";
  out << "            boundary. Each allowed trajectory contributes unit weight:\n";
  out << "              T(E;x0) = N_allowed / N_dirs\n";
  out << "              J_loc(E;x0) = T(E;x0) * J_b(E)\n\n";
  out << "        ANISOTROPIC\n";
  out << "            Boundary spectrum is factored as:\n";
  out << "              J_b(E, dir, x) = J_b_iso(E) * f_PAD(cos_alpha) * f_spatial(x)\n";
  out << "            where cos_alpha = v_exit . B_hat(x_exit) is the pitch-angle cosine\n";
  out << "            at the domain boundary crossing for each backtraced trajectory.\n";
  out << "            Transmissivity becomes an anisotropy-weighted sum:\n";
  out << "              T_aniso(E;x0) = (1/N_dirs) * sum_k A_k * f_PAD_k * f_spatial_k\n";
  out << "              J_loc(E;x0) = J_b_iso(E) * T_aniso(E;x0)\n";
  out << "            Requires a #BOUNDARY_ANISOTROPY section in the input file.\n";
  out << "            Physically appropriate for: SEP events (field-aligned PAD),\n";
  out << "            radiation belt pancake distributions, bidirectional streaming,\n";
  out << "            day/night asymmetric CME-driven boundary conditions.\n\n";
  out << "  -density-transmission-mode <DIRECT|SCAN|ADAPTIVE>   (optional)\n";
  out << "      Select how density/flux samples the transmission function T(E).\n";
  out << "      DIRECT keeps the legacy DS energy grid. SCAN and ADAPTIVE evaluate\n";
  out << "      T on a log-spaced rigidity grid converted back to kinetic energy,\n";
  out << "      which better resolves cutoff/penumbra structure before folding with\n";
  out << "      the boundary spectrum. Applies to both gridless and Mode3D density/flux.\n";
  out << "      Optional controls: --density-transmission-scan-n <N>,\n";
  out << "      --density-transmission-refine-n <N>, --density-transmission-max-n <N>.\n";
  out << "      Input-file equivalents are DS_TRANSMISSION_MODE and\n";
  out << "      DS_TRANSMISSION_SCAN_N/REFINE_N/MAX_N.\n\n";

  out << "  -mode3d-output-initialized | --mode3d-output-initialized   (optional; default: off)\n";
  out << "      In -mode 3d, write amps_3d_initialized.data.dat after the AMR mesh\n";
  out << "      cell-centered B/E fields have been initialized. This diagnostic file can\n";
  out << "      be large, so it is skipped unless this flag is present.\n\n";

  out << "  -mode3d-field-eval | --mode3d-field-eval <INTERPOLATION|ANALYTIC>   (optional; default: INTERPOLATION)\n";
  out << "      In -mode 3d, select the magnetic-field source used during backtracing.\n";
  out << "        INTERPOLATION  Use the AMR-aware cell-centered interpolation stencil.\n";
  out << "        ANALYTIC       Call the same background-field evaluator used to prepopulate\n";
  out << "                       the mesh cell centers. Alias: --mode3d-analytic-field.\n\n";

  out << "  -cutoff-debug-scan | --cutoff-debug-scan <lon_deg> <lat_deg> <alt_km>   (optional)\n";
  out << "      In -mode 3d cutoff runs, write a rigidity-classification table for one\n";
  out << "      spherical-shell location before the full cutoff calculation starts. The\n";
  out << "      table lists R_GV, TraceAllowed3D(R), and the analytic Störmer vertical\n";
  out << "      cutoff for DIPOLE fields. It is intended for diagnosing cases where the\n";
  out << "      endpoint test returns the lower rigidity bound on large shell grids.\n";
  out << "      Optional controls: --cutoff-debug-scan-n <N>,\n";
  out << "      --cutoff-debug-scan-file <file>. Input-file equivalents are\n";
  out << "      CUTOFF_DEBUG_RIGIDITY_SCAN, CUTOFF_DEBUG_SCAN_LON/LAT/ALT/N/FILE.\n\n";

  out << "  -cutoff-debug-exit | --cutoff-debug-exit <lon_deg> <lat_deg> <alt_km>   (optional)\n";
  out << "      In -mode 3d cutoff runs, write a trajectory-exit diagnostic for one\n";
  out << "      selected vertical backtrace. The file records terminal reason, raw exit\n";
  out << "      point, reconstructed boundary-crossing point and face, box overshoot,\n";
  out << "      rigidity conservation, and the centered-dipole canonical invariant.\n";
  out << "      Optional controls: --cutoff-debug-exit-r <R_GV>,\n";
  out << "      --cutoff-debug-exit-n <N>, --cutoff-debug-exit-file <file>.\n";
  out << "      Input-file equivalents are CUTOFF_DEBUG_EXIT_TRACE and\n";
  out << "      CUTOFF_DEBUG_EXIT_LON/LAT/ALT/R_GV/N/FILE.\n\n";

  out << "  -cutoff-search | --cutoff-search <UPPER_SCAN|BINARY>   (optional)\n";
  out << "      Select the scalar cutoff-rigidity search algorithm in both standalone\n";
  out << "      -mode 3d and -mode gridless cutoff runs. UPPER_SCAN is the default and\n";
  out << "      is penumbra-safe: it first evaluates a log-spaced rigidity grid from\n";
  out << "      Rmin to Rmax, scans from high to low rigidity for the highest forbidden\n";
  out << "      sample, then bisects only that final forbidden/allowed interval. BINARY\n";
  out << "      restores the legacy endpoint bisection and assumes TraceAllowed(R) is\n";
  out << "      monotonic. Gridless aliases: --gridless-cutoff-search and\n";
  out << "      --gridless-cutoff-search-n. Optional control: --cutoff-upper-scan-n <N>.\n";
  out << "      If omitted, UPPER_SCAN uses CUTOFF_NENERGY samples. Input-file\n";
  out << "      equivalents are CUTOFF_SEARCH_ALGORITHM and CUTOFF_UPPER_SCAN_N.\n\n";

  out << "  -density-parallel | --density-parallel <OPENMP|THREADS|SERIAL>   (optional)\n";
  out << "      Select the Mode3D density-backtracking shared-memory backend inside\n";
  out << "      each MPI process. OPENMP preserves the existing OpenMP loops; THREADS\n";
  out << "      uses direct std::thread workers over observation locations and disables\n";
  out << "      nested OpenMP inside those workers; SERIAL disables intra-rank shared\n";
  out << "      memory parallelism. Aliases: --mode3d-density-parallel,\n";
  out << "      --mode3d-density-backend.\n\n";

  out << "  -density-threads | --density-threads <N>   (optional; 0=automatic)\n";
  out << "      Number of shared-memory workers per MPI process for Mode3D density\n";
  out << "      backtracking. For OPENMP this calls omp_set_num_threads(N). For THREADS\n";
  out << "      this sets the number of std::thread workers. Aliases:\n";
  out << "      --mode3d-density-threads, --n-density-threads.\n\n";
  out << "  -mode3d-mpi-scheduler | --mode3d-mpi-scheduler <DYNAMIC|BLOCK_CYCLIC|STATIC>   (optional)\n";
  out << "      Select the inter-rank scheduler for standalone Mode3D cutoff and\n";
  out << "      density/flux backtracking. DYNAMIC uses an MPI one-sided atomic work\n";
  out << "      queue so ranks fetch chunks of global locations as soon as they become\n";
  out << "      idle. BLOCK_CYCLIC is the previous deterministic rank r, r+Nrank, ...\n";
  out << "      partition. STATIC assigns one contiguous slab per rank. Aliases:\n";
  out << "      --backtrack-mpi-scheduler, --mode3d-mpi-backend.\n\n";

  out << "  -mode3d-mpi-dynamic-chunk | --mode3d-mpi-dynamic-chunk <N>   (optional; 0=automatic)\n";
  out << "      Number of global observation locations fetched per MPI atomic request\n";
  out << "      when the scheduler is DYNAMIC. Smaller values improve load balance;\n";
  out << "      larger values reduce MPI scheduler overhead. Alias: --backtrack-mpi-chunk.\n\n";

  out << "  -max-trace-distance | --max-trace-distance <double>   (optional)\n";
  out << "\n";
  out << "3D Forward mode options (-mode 3d_forward):\n\n";
  out << "  -forward-niter | --forward-niter <int>   (optional)\n";
  out << "      Override FORWARD_N_ITERATIONS from the input file.\n";
  out << "      Total number of forward time-step iterations to run.\n\n";
  out << "  -forward-nparticles | --forward-nparticles <int>   (optional; default: 1000)\n";
  out << "      Number of simulation particles injected at the domain boundary each\n";
  out << "      forward iteration.  Overrides FORWARD_N_PARTICLES from the input file.\n";
  out << "      The physical particle weight is computed automatically from the boundary\n";
  out << "      spectrum, total boundary area, and the time step:\n";
  out << "        W = (pi x integral_{Emin}^{Emax} J(E) dE x A_boundary x dt) / N\n";
  out << "      where N is this count.  Larger N gives finer phase-space sampling at\n";
  out << "      the cost of proportionally more memory and compute per iteration.\n\n";
  out << "  -forward-injection-energy | --forward-injection-energy <SPECTRUM|LOG_UNIFORM>   (optional)\n";
  out << "      Select the injected-particle energy proposal distribution.\n";
  out << "        SPECTRUM    (default): sample E from the physical J(E)dE CDF.\n";
  out << "                     This is the legacy equal-weight branch and is efficient\n";
  out << "                     where the spectrum is largest, usually at low energy.\n";
  out << "        LOG_UNIFORM: sample E uniformly in log(E) and apply an individual\n";
  out << "                     statistical-weight correction q(E)=p_target/p_log.\n";
  out << "                     This improves high-energy statistics for steep SEP spectra\n";
  out << "                     while conserving the total physical source rate per step.\n";
  out << "      Alias: --forward-energy-sampling. A matching input-file key can be added\n";
  out << "      later using Mode3DForwardOptions::injectionEnergyDistribution.\n\n";

  out << "  -forward-injection-emin | --forward-injection-emin <MeV/n>   (optional)\n";
  out << "  -forward-injection-emax | --forward-injection-emax <MeV/n>   (optional)\n";
  out << "      Override the 3d_forward particle-energy limits from the command line.\n";
  out << "      These values update the effective #DENSITY_3D DENS_EMIN/DENS_EMAX range,\n";
  out << "      which controls both the volumetric density energy bins and the forward\n";
  out << "      boundary injection/integration range.  Aliases: --forward-emin,\n";
  out << "      --forward-emax, --forward-injection-energy-min, --forward-injection-energy-max.\n";
  out << "      Convenience form: --forward-injection-energy-range <Emin> <Emax>.\n\n";

  out << "  -forward-boundary-dist | --forward-boundary-dist <ISOTROPIC>   (optional)\n";
  out << "      Override Mode3DForwardOptions::boundaryDistType.\n";
  out << "      ISOTROPIC (default): cos(θ)-weighted hemisphere at each domain-boundary face.\n";
  out << "      Other distributions (field-aligned beam, pancake) are reserved for future use.\n\n";
  out << "  -mode3d-output-initialized\n";
  out << "      Also works in 3d_forward mode: writes amps_3dforward_initialized.data.dat\n";
  out << "      after the AMR mesh B/E fields have been populated.\n\n";
  out << "\n";
  out << "      Override #NUMERICAL / MAX_TRACE_DISTANCE from the input file.\n";
  out << "      Units: Earth radii (Re) of cumulative traced path length.\n";
  out << "        value > 0   enable the cumulative-distance cap\n";
  out << "        value = 0   disable the cap\n";
  out << "      Applies to the shared backtracer, so it affects cutoff-rigidity,\n";
  out << "      density, and spectrum calculations consistently.\n\n";
  out << "      Example: compare both modes on the same input file:\n";
  out << "        " << progName << " -mode gridless -i run.in -density-mode ISOTROPIC\n";
  out << "        " << progName << " -mode gridless -i run.in -density-mode ANISOTROPIC\n\n";

  // -----------------------------------------------------------------------
  out << "Input file sections (AMPS_PARAM format):\n\n";

  out << "  #CALCULATION_MODE\n";
  out << "    CALC_TARGET        CUTOFF_RIGIDITY | DENSITY_SPECTRUM\n";
  out << "    FIELD_EVAL_METHOD  GRIDLESS | GRID_3D\n\n";

  out << "  #DENSITY_SPECTRUM   (required when CALC_TARGET = DENSITY_SPECTRUM)\n";
  out << "    DS_EMIN            <double>   MeV/n; lower energy bound\n";
  out << "    DS_EMAX            <double>   MeV/n; upper energy bound\n";
  out << "    DS_NINTERVALS      <int>      energy intervals (nPoints = nIntervals+1)\n";
  out << "    DS_ENERGY_SPACING  LOG | LINEAR\n";
  out << "    DS_MAX_PARTICLES   <int>      total trajectory cap per obs. point (opt)\n";
  out << "    DS_MAX_TRAJ_TIME   <double>   per-trajectory time cap in seconds (opt)\n";
  out << "    DS_BOUNDARY_MODE   ISOTROPIC | ANISOTROPIC    (overridable by -density-mode)\n\n";

  out << "  #BOUNDARY_ANISOTROPY   (required when DS_BOUNDARY_MODE = ANISOTROPIC)\n";
  out << "    BA_PAD_MODEL       ISOTROPIC | SINALPHA_N | COSALPHA_N | BIDIRECTIONAL\n";
  out << "      ISOTROPIC      f(alpha) = 1                (sanity-check identity)\n";
  out << "      SINALPHA_N     f(alpha) = sin^n(alpha)      pancake; equatorial belt\n";
  out << "      COSALPHA_N     f(alpha) = |cos(alpha)|^n    field-aligned beam; SEPs\n";
  out << "      BIDIRECTIONAL  f(alpha) = |cos(alpha)|^n    symmetric; closed-loop\n";
  out << "    BA_PAD_EXPONENT    <double>   exponent n (default 2.0; 0 = isotropic)\n";
  out << "    BA_SPATIAL_MODEL   UNIFORM | DAYSIDE_NIGHTSIDE\n";
  out << "      UNIFORM          f_spatial = 1 everywhere\n";
  out << "      DAYSIDE_NIGHTSIDE\n";
  out << "                       f_spatial = BA_DAYSIDE_FACTOR   if GSM x_exit > 0\n";
  out << "                                 = BA_NIGHTSIDE_FACTOR  if GSM x_exit <= 0\n";
  out << "    BA_DAYSIDE_FACTOR  <double>   sunward multiplier    (default 1.0)\n";
  out << "    BA_NIGHTSIDE_FACTOR <double>  tailward multiplier   (default 1.0)\n\n";

  out << "  #BACKGROUND_FIELD\n";
  out << "    FIELD_MODEL        T96 | T05 | DIPOLE\n";
  out << "    EPOCH              <ISO datetime>    e.g. 2003-11-20T06:00\n";
  out << "    DST / PDYN / IMF_BY / IMF_BZ         driving parameters [nT / nPa]\n";
  out << "    T05_W1 .. T05_W6   storm-time history integrals (T05 only)\n";
  out << "    DIPOLE_MOMENT      <double>   multiple of Earth dipole M_E (DIPOLE only)\n";
  out << "    DIPOLE_TILT        <double>   tilt from +Z_GSM toward +X_GSM [deg]\n\n";

  out << "  #DOMAIN_BOUNDARY\n";
  out << "    DOMAIN_X/Y/ZMIN/MAX <double>  outer rectangular box [km, GSM]\n";
  out << "    R_INNER             <double>  inner loss sphere radius [km]\n\n";

  out << "  #OUTPUT_DOMAIN\n";
  out << "    OUTPUT_MODE        POINTS | SHELLS\n";
  out << "    SHELL_ALTITUDES    <double ...>   km above Earth's surface\n";
  out << "    SHELL_RES          <double>       lon/lat grid resolution [deg]\n";
  out << "    POINTS_BEGIN       (followed by 'x y z' rows in km GSM)\n";
  out << "    POINTS_END\n\n";

  out << "  #NUMERICAL\n";
  out << "    DT_TRACE           <double>   initial time step [s]\n";
  out << "    MAX_STEPS          <int>      hard step count cap\n";
  out << "    MAX_TRACE_TIME     <double>   hard integration time cap [s]\n";
  out << "    MAX_TRACE_DISTANCE <double>   hard cumulative trace-distance cap [Re]\n\n";

  // -----------------------------------------------------------------------
  out << "Outputs (gridless mode):\n\n";

  out << "  CUTOFF_RIGIDITY:\n";
  out << "    gridless_points_cutoff.dat           (POINTS mode)\n";
  out << "    gridless_shell_<A>km_cutoff.dat      (SHELLS mode, one per altitude A)\n";
  out << "    Variables: X_km Y_km Z_km R_km Lon_deg Lat_deg Rc_GV Emin_MeV\n\n";

  out << "  DENSITY_SPECTRUM:\n";
  out << "    gridless_points_density.dat          total density per observation point\n";
  out << "    gridless_points_spectrum.dat         per-point spectrum (one ZONE each)\n";
  out << "    gridless_shell_<A>km_density_channels.dat  (SHELLS mode)\n";
  out << "    Variables: E_MeV  T  J_boundary_perMeV  J_local_perMeV  N_m3  N_cm3\n\n";

  // -----------------------------------------------------------------------
  out << "Notes:\n";
  out << "  * All input positions assumed in GSM coordinates [km]; no transform applied.\n";
  out << "  * -mover and -density-mode override the corresponding input file settings.\n";
  out << "  * Mode3D writes amps_3d_initialized.data.dat only when\n";
  out << "    -mode3d-output-initialized is supplied.\n";
  out << "  * Mode3D magnetic-field tracing uses interpolation by default; use\n";
  out << "    -mode3d-field-eval ANALYTIC for direct background-field evaluation.\n";
  out << "  * Mode3D density backtracking can use -density-parallel OPENMP or THREADS;\n";
  out << "    use -density-threads N to set workers per MPI process.\n";
  out << "  * MPI: rank 0 is master scheduler; ranks 1..N-1 are workers.\n";
  out << "    Dynamic point-level scheduling; no static decomposition.\n";
  out << "  * See AnisotropicSpectrum.h for the full physics derivation of the\n";
  out << "    ANISOTROPIC mode and guidance on extending the boundary spectrum model.\n";

  return out.str();
}

}
