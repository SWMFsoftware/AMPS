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
//       See GridlessParticleMovers.h for a full description of each mover.
//       If omitted, the default mover (BORIS) is used.
//       The string is stored as-is; translation to MoverType enum is done by the caller.
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
//   -max-trace-distance <double>
//       Override #NUMERICAL MAX_TRACE_DISTANCE from the input file.
//       Units: Earth radii (Re) of cumulative traced path length.
//       Semantics:
//         value > 0   enable the hard cumulative-distance cap
//         value = 0   disable the cap
//       This mirrors MAX_TRACE_TIME, but limits total geometric distance traveled
//       by a trajectory rather than elapsed integration time.
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
    // a concrete enum (MoverType) and/or an input-file setting.
    //
    // Supported values (case-insensitive):
    //   BORIS           : classic relativistic Boris pusher (legacy default)
    //   RK2/RK4/RK6   : explicit full-orbit Runge-Kutta movers of order 2/4/6
    //   GC2/GC4/GC6   : guiding-center movers integrated with RK2/RK4/RK6
    //   HYBRID        : per-step switch between RK4 and GC4 using a local
    //                   adiabaticity criterion rho/L_eff
    //
    // If empty, the executable should use its default / input-file setting.
    std::string mover{""};

    // -density-mode ISOTROPIC|ANISOTROPIC
    // Overrides DS_BOUNDARY_MODE from the input file.
    // ISOTROPIC   : uniform isotropic boundary (original behavior).
    // ANISOTROPIC : pitch-angle-dependent and spatially non-uniform boundary;
    //               requires a #BOUNDARY_ANISOTROPY section in the input file.
    std::string densityMode{""};

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

    // -----------------------------------------------------------------------
    // #DENSITY_3D overrides (-mode 3d_forward)
    // -----------------------------------------------------------------------
    // These flags override the corresponding keys in the #DENSITY_3D section.
    // Sentinel: double < 0 | int < 0 | string empty => no CLI override.
    //
    //   -density3d-emin    <double>       Override DENS_EMIN  [MeV/n]
    //   -density3d-emax    <double>       Override DENS_EMAX  [MeV/n]
    //   -density3d-nenergy <int>          Override DENS_NENERGY
    //   -density3d-spacing <LOG|LINEAR>   Override DENS_ENERGY_SPACING
    double density3dEmin_MeV{-1.0};   ///< <0 = no override
    double density3dEmax_MeV{-1.0};   ///< <0 = no override
    int    density3dNenergy{-1};       ///< <0 = no override
    std::string density3dSpacing{""};  ///< empty = no override
  };

  // Parse argc/argv. Throws std::runtime_error for malformed inputs.
  CliOptions ParseCli(int argc,char** argv);

  // Return formatted help message.
  std::string HelpMessage(const char* progName);

}

#endif
