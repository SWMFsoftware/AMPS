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
  out << "  " << progName << " -mode 3d       -i AMPS_PARAM.in\n";
  out << "\n";

  // -----------------------------------------------------------------------
  out << "Options:\n\n";

  out << "  -h | --help\n";
  out << "      Print this message and exit. No solver is run.\n\n";

  out << "  -mode <3d|gridless>   (required)\n";
  out << "      Select the solver execution path.\n";
  out << "        3d        Full PIC-backed 3D solver with mesh-interpolated fields.\n";
  out << "        gridless  Direct Tsyganenko/IGRF field evaluation + backtracing.\n";
  out << "                  Supports both cutoff-rigidity and density/spectrum modes.\n\n";

  out << "  -i <file>   (required)\n";
  out << "      Path to the AMPS_PARAM-format input file. Controls all physics and\n";
  out << "      geometry parameters. See 'Input file sections' below for a summary.\n\n";

  out << "  -mover <BORIS|RK2|RK4|RK6|GC2|GC4|GC6|HYBRID>   (optional; default: BORIS)\n";
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

  out << "  -density-mode <ISOTROPIC|ANISOTROPIC>   (optional)\n";
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
  out << "    MAX_TRACE_TIME     <double>   hard integration time cap [s]\n\n";

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
  out << "  * MPI: rank 0 is master scheduler; ranks 1..N-1 are workers.\n";
  out << "    Dynamic point-level scheduling; no static decomposition.\n";
  out << "  * See AnisotropicSpectrum.h for the full physics derivation of the\n";
  out << "    ANISOTROPIC mode and guidance on extending the boundary spectrum model.\n";

  return out.str();
}

}
