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
    // ----- #DENSITY_3D overrides (3d_forward mode) -----
    else if (a=="-density3d-emin" || a=="--density3d-emin") {
      if (i+1>=argc) exit(__LINE__,__FILE__,"Missing value after -density3d-emin");
      opt.density3dEmin_MeV = std::stod(argv[++i]);
      if (opt.density3dEmin_MeV <= 0.0)
        exit(__LINE__,__FILE__,"-density3d-emin must be > 0 MeV/n");
    }
    else if (a=="-density3d-emax" || a=="--density3d-emax") {
      if (i+1>=argc) exit(__LINE__,__FILE__,"Missing value after -density3d-emax");
      opt.density3dEmax_MeV = std::stod(argv[++i]);
      if (opt.density3dEmax_MeV <= 0.0)
        exit(__LINE__,__FILE__,"-density3d-emax must be > 0 MeV/n");
    }
    else if (a=="-density3d-nenergy" || a=="--density3d-nenergy") {
      if (i+1>=argc) exit(__LINE__,__FILE__,"Missing value after -density3d-nenergy");
      opt.density3dNenergy = std::stoi(argv[++i]);
      if (opt.density3dNenergy < 1)
        exit(__LINE__,__FILE__,"-density3d-nenergy must be >= 1");
    }
    else if (a=="-density3d-spacing" || a=="--density3d-spacing") {
      if (i+1>=argc) exit(__LINE__,__FILE__,"Missing value after -density3d-spacing");
      opt.density3dSpacing = argv[++i];
      // Validate at parse time for early error reporting
      std::string up = opt.density3dSpacing;
      for (auto& c : up) c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
      if (up != "LOG" && up != "LINEAR")
        exit(__LINE__,__FILE__,"-density3d-spacing must be LOG or LINEAR");
    }
    else {
      std::ostringstream oss;
      oss << "Unknown CLI token: '" << a << "'. Use -h for help.";
      exit(__LINE__,__FILE__,oss.str().c_str());
    }
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
  out << "  -mode3d-output-initialized | --mode3d-output-initialized   (optional; default: off)\n";
  out << "      In -mode 3d, write amps_3d_initialized.data.dat after the AMR mesh\n";
  out << "      cell-centered B/E fields have been initialized. This diagnostic file can\n";
  out << "      be large, so it is skipped unless this flag is present.\n\n";

  out << "  -mode3d-field-eval | --mode3d-field-eval <INTERPOLATION|ANALYTIC>   (optional; default: INTERPOLATION)\n";
  out << "      In -mode 3d, select the magnetic-field source used during backtracing.\n";
  out << "        INTERPOLATION  Use the AMR-aware cell-centered interpolation stencil.\n";
  out << "        ANALYTIC       Call the same background-field evaluator used to prepopulate\n";
  out << "                       the mesh cell centers. Alias: --mode3d-analytic-field.\n\n";

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
  out << "  -forward-boundary-dist | --forward-boundary-dist <ISOTROPIC>   (optional)\n";
  out << "      Override Mode3DForwardOptions::boundaryDistType.\n";
  out << "      ISOTROPIC (default): cos(θ)-weighted hemisphere at each domain-boundary face.\n";
  out << "      Other distributions (field-aligned beam, pancake) are reserved for future use.\n\n";
  out << "  -mode3d-output-initialized\n";
  out << "      Also works in 3d_forward mode: writes amps_3dforward_initialized.data.dat\n";
  out << "      after the AMR mesh B/E fields have been populated.\n\n";
  out << "\n";

  // -----------------------------------------------------------------------
  out << "Energy-binned density options (-mode 3d_forward, overrides #DENSITY_3D):\n\n";

  out << "  -density3d-emin | --density3d-emin <double>   (optional)\n";
  out << "      Override DENS_EMIN from the #DENSITY_3D input section.\n";
  out << "      Lower bound of the energy grid for density sampling.\n";
  out << "      Units: MeV/n (kinetic energy per nucleon).  Must be > 0.\n";
  out << "      Input file default: 1.0 MeV/n\n\n";

  out << "  -density3d-emax | --density3d-emax <double>   (optional)\n";
  out << "      Override DENS_EMAX from the #DENSITY_3D input section.\n";
  out << "      Upper bound of the energy grid for density sampling.\n";
  out << "      Units: MeV/n (kinetic energy per nucleon).  Must be > 0.\n";
  out << "      Input file default: 20000.0 MeV/n\n\n";

  out << "  -density3d-nenergy | --density3d-nenergy <int>   (optional)\n";
  out << "      Override DENS_NENERGY from the #DENSITY_3D input section.\n";
  out << "      Number of energy bins covering [DENS_EMIN, DENS_EMAX].\n";
  out << "      Must be >= 1.  Larger values give finer energy resolution at the\n";
  out << "      cost of proportionally more sampling-buffer memory per cell.\n";
  out << "      Input file default: 30\n\n";

  out << "  -density3d-spacing | --density3d-spacing <LOG|LINEAR>   (optional)\n";
  out << "      Override DENS_ENERGY_SPACING from the #DENSITY_3D input section.\n";
  out << "        LOG     Logarithmically spaced bins (recommended for broad energy\n";
  out << "                ranges spanning orders of magnitude; GCR/SEP spectra).\n";
  out << "        LINEAR  Uniformly spaced bins (equal width in energy; useful for\n";
  out << "                narrow ranges or when spectral resolution is uniform).\n";
  out << "      Input file default: LOG\n\n";

  out << "  Input file equivalent (#DENSITY_3D section):\n";
  out << "    #DENSITY_3D\n";
  out << "    DENS_EMIN              1.0          ! MeV/n\n";
  out << "    DENS_EMAX              20000.0      ! MeV/n\n";
  out << "    DENS_NENERGY           30           ! energy bins\n";
  out << "    DENS_ENERGY_SPACING    LOG          ! LOG or LINEAR\n\n";

  out << "  Example: override energy grid on the command line:\n";
  out << "    " << progName << " -mode 3d_forward -i run.in \\\n";
  out << "        -density3d-emin 10.0 -density3d-emax 1000.0 \\\n";
  out << "        -density3d-nenergy 20 -density3d-spacing LOG\n\n";
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
  out << "  * MPI: rank 0 is master scheduler; ranks 1..N-1 are workers.\n";
  out << "    Dynamic point-level scheduling; no static decomposition.\n";
  out << "  * See AnisotropicSpectrum.h for the full physics derivation of the\n";
  out << "    ANISOTROPIC mode and guidance on extending the boundary spectrum model.\n";

  return out.str();
}

}
