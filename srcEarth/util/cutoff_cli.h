//======================================================================================
// cutoff_cli.h
//======================================================================================
//
// PURPOSE
// -------
// Command-line interface (CLI) for the Earth energetic particle gridless tools.
// Parses argc/argv into a CliOptions struct that the main executable uses to
// configure the solver run. Intentionally kept independent of the solver and parser
// so it can be unit-tested in isolation.
//
//======================================================================================
// SUPPORTED OPTIONS
//======================================================================================
//
//   -h | --help
//       Print help text and exit. No other processing is done.
//
//   -mode <string>
//       Select the solver execution mode.
//       Recognised values:
//         3d        Run the full PIC-backed 3D solver.
//         gridless  Run the gridless field-evaluation solver (Tsyganenko + Boris).
//       Required. Error if absent and -h is not given.
//
//   -i <path>
//       Path to the AMPS_PARAM-format input file. Parsed by ParseAmpsParamFile.
//       Required unless -h is given.
//
//   -mover <string>
//       Select the particle integration algorithm.
//       Recognised values (case-insensitive):
//         BORIS   Relativistic Boris pusher (default; recommended for all production runs)
//         RK2     Runge-Kutta 2nd order (Heun)
//         RK4     Runge-Kutta 4th order (classical)
//         RK6     Runge-Kutta 6th order
//         GC2     Guiding-center equations integrated with RK2
//         GC4     Guiding-center equations integrated with RK4
//         GC6     Guiding-center equations integrated with RK6
//         HYBRID  Switch per step between RK4 and GC4 when the local motion is
//                 sufficiently adiabatic for guiding-center transport
//       In gridless and backward 3D modes, all listed values are handled by the
//       shared GridlessParticleMovers layer. In 3d_forward mode, the AMPS-signature
//       mover manager currently supports BORIS, RK4, GC/GC4, and HYBRID. HC4 is available in gridless/backward 3-D tracing.
//       If omitted, the default mover (BORIS) is used.
//       The string is stored as-is; translation to the concrete mover is done by the caller.
//
//   -density-mode <string>
//       Override the DS_BOUNDARY_MODE key from the input file.
//       Recognised values (case-insensitive):
//         ISOTROPIC     Use uniform isotropic boundary spectrum (original behavior).
//         ANISOTROPIC   Use pitch-angle- and spatially-weighted boundary spectrum.
//                       Requires a #BOUNDARY_ANISOTROPY section in the input file,
//                       or the solver will throw at startup.
//       If omitted, the input file value (or its default ISOTROPIC) is used.
//       CLI override is useful for:
//         - Comparing isotropic and anisotropic results on the same input file
//           without editing the file: run twice with -density-mode ISOTROPIC and
//           -density-mode ANISOTROPIC.
//         - Automated test scripts that exercise both branches from a single
//           test input file.
//
//   -mode3d-output-initialized
//       In -mode 3d, write amps_3d_initialized.data.dat after mesh field
//       initialization. The default is to skip this potentially large diagnostic file.
//
//   -mode3d-field-eval <INTERPOLATION|ANALYTIC>
//       In -mode 3d, select how the magnetic field is evaluated during tracing.
//       INTERPOLATION (default) uses the AMR cell-centered interpolation stencil.
//       ANALYTIC calls the same background-field function used to initialize the
//       mesh cell centers.
//
//   -density-parallel <OPENMP|THREADS|SERIAL>
//       Select the shared-memory backend used by Mode3D density backtracking inside
//       each MPI process. OPENMP preserves the legacy OpenMP path; THREADS uses a
//       direct std::thread work queue over observation locations; SERIAL disables
//       intra-rank shared-memory parallelism.
//
//   -density-threads <int>
//       Number of shared-memory workers per MPI process for the selected density
//       backend. For OPENMP this calls omp_set_num_threads(N); for THREADS it sets
//       the number of std::thread workers.
//
//   -mode3d-mpi-scheduler <DYNAMIC|BLOCK_CYCLIC|STATIC>
//       Select the inter-rank scheduler for standalone Mode3D and gridless cutoff/density runs.
//
//   -mode3d-mpi-dynamic-chunk <int>
//       Number of Mode3D locations or gridless task chunks per dynamic MPI fetch. 0 means automatic.
//
//   -cutoff-search <UPPER_SCAN|BINARY>
//       Select the scalar cutoff-rigidity search algorithm used by both standalone
//       Mode3D and gridless cutoff.  UPPER_SCAN is penumbra-safe and is the
//       default; BINARY is the legacy endpoint method.  Mode-specific aliases are
//       also accepted, including -mode3d-cutoff-search and -gridless-cutoff-search.
//
//   -cutoff-upper-scan-n <int>
//       Number of log-spaced rigidity samples used by UPPER_SCAN before the final
//       forbidden/allowed bisection.  If omitted, the solver reuses CUTOFF_NENERGY
//       from the input file.  The aliases -mode3d-cutoff-search-n and
//       -gridless-cutoff-search-n are accepted.
//
//   -adaptive-dt <T|F>
//       Override #NUMERICAL ADAPTIVE_DT.  T/default means DT_TRACE is a maximum
//       adaptive-step bound; F means fixed-step tracing with DT_TRACE.
//
//   -max-trace-distance <double>
//       Override #NUMERICAL MAX_TRACE_DISTANCE from the input file.
//       Units: Earth radii (Re) of cumulative traced path length.
//       Semantics:
//         value > 0   enable the hard cumulative-distance cap
//         value = 0   disable the cap
//       This mirrors MAX_TRACE_TIME, but limits total geometric distance traveled
//       by a trajectory rather than elapsed integration time.
//
//   -mode3d-mesh-res-earth-re <double>
//   -mode3d-mesh-res-boundary-re <double>
//   -mode3d-mesh-coarsening <LINEAR|LOG|EXPONENTIAL|POWER|CONSTANT>
//   -mode3d-mesh-exponent <double>
//   -mode3d-mesh-r-boundary-re <double>
//       Optional standalone -mode 3d AMR mesh-resolution profile.  If omitted,
//       the current hard-coded main_lib.cpp localResolution() profile is used.
//
//======================================================================================
// USAGE EXAMPLES
//======================================================================================
//
//   Basic cutoff-rigidity run (isotropic, Boris pusher):
//     ./amps -mode gridless -i run.in
//
//   Density/spectrum, anisotropic boundary, override from command line:
//     ./amps -mode gridless -i run.in -density-mode ANISOTROPIC
//
//   Density/spectrum, RK4 mover for comparison:
//     ./amps -mode gridless -i run.in -mover RK4
//
//   Print help:
//     ./amps -h
//
//======================================================================================
// ERROR HANDLING
//======================================================================================
//
// ParseCli throws std::runtime_error for:
//   - An option flag that requires an argument but none follows (e.g., "-i" at end)
//   - An unrecognised option flag (e.g., "-xyz")
//
// It does NOT throw for:
//   - Missing required options (the main executable checks CliOptions after parsing)
//   - Unrecognised values for -mover or -density-mode (passed through as strings;
//     validated by the solver at startup)
//
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
    // a concrete enum (MoverType) or the 3d_forward manager selection.
    //
    // Supported values (case-insensitive):
    //   BORIS           : classic relativistic Boris pusher (legacy default)
    //   HC4           : fourth-order Higuera-Cary/Boris composition; RK4-like accuracy with Boris-like rigidity conservation
    //   RK2/RK4/RK6   : explicit full-orbit Runge-Kutta movers of order 2/4/6
    //   GC2/GC4/GC6   : guiding-center movers integrated with RK2/RK4/RK6
    //   HYBRID        : per-step switch between RK4 and GC4 using a local
    //                   adiabaticity criterion rho/L_eff
    //
    // If empty, the executable should use its default setting.
    std::string mover{""};

    // -density-mode ISOTROPIC|ANISOTROPIC
    // Overrides DS_BOUNDARY_MODE from the input file.
    // ISOTROPIC   : uniform isotropic boundary (original behavior).
    // ANISOTROPIC : pitch-angle-dependent and spatially non-uniform boundary;
    //               requires a #BOUNDARY_ANISOTROPY section in the input file.
    std::string densityMode{""};

    // -density-transmission-mode DIRECT|SCAN|ADAPTIVE
    // Controls the energy/rigidity nodes on which density/flux transmissivity T is
    // evaluated.  DIRECT keeps the legacy DS energy grid.  SCAN/ADAPTIVE build a
    // log-spaced rigidity grid converted back to kinetic energy for spectrum folding.
    std::string densityTransmissionMode{""};

    // -density-transmission-scan-n <int>
    // Number of rigidity scan points for SCAN/ADAPTIVE mode. 0 means no CLI override.
    int densityTransmissionScanN{0};

    // -density-transmission-refine-n <int> and -density-transmission-max-n <int>
    // Parsed for forward-compatible input decks; current production output uses the
    // fixed scan grid so all locations have the same Tecplot energy axis.
    int densityTransmissionRefineN{0};
    int densityTransmissionMaxN{0};

    // -mode3d-output-initialized
    // Boolean flag. When true, Mode3D::Run writes the initialized AMR mesh fields
    // to amps_3d_initialized.data.dat. The default is false to avoid creating this
    // large diagnostic file unless explicitly requested.
    bool mode3dOutputInitialized{false};

    // -mode3d-field-eval <INTERPOLATION|ANALYTIC>
    // Optional Mode3D magnetic-field evaluation override.
    // Empty or INTERPOLATION uses the AMR interpolation stencil.
    // ANALYTIC calls Earth::Mode3D::EvaluateBackgroundMagneticFieldSI directly.
    std::string mode3dFieldEval{""};

    // -density-parallel <OPENMP|THREADS|SERIAL>
    // Optional Mode3D backward-product shared-memory backend override.
    // Empty string means: use default/environment.
    std::string densityParallelBackend{""};

    // -density-threads <int>
    // Number of shared-memory workers per MPI process. 0 means automatic/default.
    int densityThreads{0};

    // -mode3d-mpi-scheduler <DYNAMIC|BLOCK_CYCLIC|STATIC>
    // Inter-rank scheduler for standalone Mode3D and gridless backtracking products.
    // Empty string means: use input-file/default value.
    std::string mode3dMpiScheduler{""};

    // -mode3d-mpi-dynamic-chunk <int>
    // Number of locations/tasks fetched per MPI atomic request in DYNAMIC mode.
    // 0 means: no CLI override / automatic when used from input defaults.
    int mode3dMpiDynamicChunk{0};

    // -cutoff-debug-scan <lon_deg> <lat_deg> <alt_km>
    // Optional Mode3D cutoff diagnostic.  When enabled, rank 0 writes a rigidity
    // classification table for the selected spherical-shell location before the
    // full cutoff calculation starts.  This is useful for dipole/Störmer tests and
    // for diagnosing endpoint/bracketing problems.
    bool cutoffDebugScan{false};
    double cutoffDebugScanLon_deg{0.0};
    double cutoffDebugScanLat_deg{0.0};
    double cutoffDebugScanAlt_km{-1.0};
    int cutoffDebugScanN{0};
    std::string cutoffDebugScanFile{""};

    // -cutoff-debug-exit <lon_deg> <lat_deg> <alt_km>
    // -cutoff-debug-exit-list <file>
    // Optional Mode3D trajectory-exit diagnostic.  It writes the terminal reason,
    // raw exit point, reconstructed boundary-crossing point, rigidity conservation,
    // and dipole canonical-momentum invariant.  The list form traces many
    // lon/lat/alt/R cases in one AMPS run and still produces one output file.
    bool cutoffDebugExit{false};
    double cutoffDebugExitLon_deg{0.0};
    double cutoffDebugExitLat_deg{0.0};
    double cutoffDebugExitAlt_km{-1.0};
    double cutoffDebugExitR_GV{-1.0};
    int cutoffDebugExitN{0};
    std::string cutoffDebugExitListFile{""};
    std::string cutoffDebugExitFile{""};

    // -cutoff-search <UPPER_SCAN|BINARY>
    // Optional override for #CUTOFF_RIGIDITY / CUTOFF_SEARCH_ALGORITHM.
    // This generic option, together with -mode3d-cutoff-search and
    // -gridless-cutoff-search aliases, is applied in both standalone backward
    // modes by ApplyCommonBackwardCli() in srcEarth/main.cpp.
    // Empty means use the input-file/default value.
    std::string cutoffSearchAlgorithm{""};

    // -cutoff-upper-scan-n <int>
    // Optional override for #CUTOFF_RIGIDITY / CUTOFF_UPPER_SCAN_N.
    // Applies to both Mode3D and gridless UPPER_SCAN searches.
    // 0 means no CLI override.
    int cutoffUpperScanN{0};

    // -adaptive-dt T|F
    // Optional CLI override for #NUMERICAL ADAPTIVE_DT.
    // Sentinel convention:
    //   -1 : no CLI override supplied; use input file/default value
    //    0 : fixed-step tracing; DT_TRACE is used directly
    //    1 : adaptive-step tracing; DT_TRACE is the maximum allowed step
    int adaptiveDt{-1};

    // -max-trace-distance <double>
    // Optional CLI override for #NUMERICAL MAX_TRACE_DISTANCE.
    //
    // Units:
    //   Earth radii (Re) of *cumulative* traced path length.
    //
    // Sentinel convention:
    //   < 0   : no CLI override was supplied; use the input file value
    //   = 0   : explicitly disable the cumulative-distance cap
    //   > 0   : enable/override the cap
    double maxTraceDistance_Re{-1.0};

    // Optional CLI override for the standalone Mode3D AMR mesh-resolution profile.
    // These values are expressed in Earth radii for CLI simplicity.  The parser
    // stores the input-file equivalents internally in km, so main.cpp converts the
    // CLI values before merging them into AmpsParam.
    //
    // Sentinel convention for the numeric fields:
    //   < 0 : no CLI override supplied
    //   > 0 : apply/override that value
    double mode3dMeshResEarth_Re{-1.0};
    double mode3dMeshResBoundary_Re{-1.0};
    double mode3dMeshOuterRadius_Re{-1.0};
    std::string mode3dMeshCoarsening{""};
    double mode3dMeshExponent{-1.0};

    // -----------------------------------------------------------------------
    // -mode 3d_forward specific options
    // -----------------------------------------------------------------------

    // -forward-niter <int>
    //   Override Mode3DForwardOptions::nIterations from the input file.
    //   Sentinel: < 0 means no CLI override (use input file value).
    int forward3dNiter{-1};

    // -forward-nparticles <int>
    //   Override Mode3DForwardOptions::nParticlesPerIter from the input file.
    //   Sets the number of simulation particles injected at the domain boundary each
    //   iteration. The physical particle weight is then automatically derived as:
    //     W = (π × ∫J(E)dE × A_boundary × dt) / nParticlesPerIter
    //   Sentinel: < 0 means no CLI override (use input file value, default 1000).
    int forward3dNparticles{-1};

    // -forward-boundary-dist <ISOTROPIC|...>
    //   Override Mode3DForwardOptions::boundaryDistType.
    //   Empty string means: use the input file default (ISOTROPIC).
    std::string forward3dBoundaryDist{""};

    // -forward-injection-energy <SPECTRUM|LOG_UNIFORM>
    //   Select the energy proposal distribution used by the 3d_forward outer-boundary
    //   particle source.
    //     SPECTRUM    : legacy/default; sample E from the physical J(E)dE CDF.
    //     LOG_UNIFORM : sample E uniformly in log(E) and correct each particle with
    //                   an individual statistical-weight factor q(E).
    //   Empty string means: use the input-file/default Mode3DForwardOptions value.
    std::string forward3dInjectionEnergyDistribution{""};

    // -forward-injection-emin <MeV/n>, -forward-injection-emax <MeV/n>
    //   Optional CLI overrides for the 3d_forward particle-energy limits.  These
    //   limits intentionally update the effective #DENSITY_3D range, because in
    //   forward mode that range now controls both the density-output energy grid and
    //   the boundary injection/integration range.
    //   Sentinel: < 0 means no CLI override.
    double forward3dInjectionEmin_MeV{-1.0};
    double forward3dInjectionEmax_MeV{-1.0};
  };

  // Parse argc/argv. Throws std::runtime_error for malformed inputs.
  CliOptions ParseCli(int argc,char** argv);

  // Return formatted help message.
  std::string HelpMessage(const char* progName);

}

#endif
