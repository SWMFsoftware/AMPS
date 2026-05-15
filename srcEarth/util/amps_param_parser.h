//======================================================================================
// amps_param_parser.h
//======================================================================================
//
// PURPOSE
// -------
// Self-contained, dependency-free parser for the AMPS_PARAM.in format used by
// the CCMC Runs-on-Request interface for the Geospace energetic particle tools.
// Populates EarthUtil::AmpsParam from a text file; used by both the gridless
// cutoff-rigidity solver and the gridless density/spectrum solver.
//
// DESIGN PRINCIPLES
// -----------------
//   (1) No PIC framework dependencies. This parser can be built and tested standalone.
//   (2) Permissive: unknown keys are stored in a "raw" map rather than rejected.
//   (3) Forward-compatible: new sections and keys can be added without breaking old runs.
//   (4) Clear unit contract: all geometric lengths are expected in km from the caller;
//       the parser documents this but does not enforce conversions (the solvers do).
//
//======================================================================================
// INPUT FILE FORMAT
//======================================================================================
//
// The file is a sequence of named sections. Each section starts with a keyword line
// beginning with '#'. Lines beginning with '!' outside a section header are comments.
// Blank lines are ignored. Within a section, each non-blank, non-comment line is:
//
//   KEY   VALUE   [! optional comment]
//
// The following sections are recognised (unrecognised sections are skipped silently):
//
//   #RUN_INFO
//     RUN_ID      <string>      ! arbitrary run identifier, stored in AmpsParam.runId
//
//   #CALCULATION_MODE
//     CALC_TARGET             CUTOFF_RIGIDITY | DENSITY_SPECTRUM
//     FIELD_EVAL_METHOD       GRIDLESS | GRID_3D
//
//   #CUTOFF_RIGIDITY
//     CUTOFF_EMIN             <double>   ! MeV; lower rigidity scan bound
//     CUTOFF_EMAX             <double>   ! MeV; upper rigidity scan bound
//     CUTOFF_NENERGY          <int>      ! number of rigidity bisection points
//     CUTOFF_MAX_PARTICLES    <int>      ! per-point trajectory cap (optional)
//     CUTOFF_MAX_TRAJ_TIME    <double>   ! per-trajectory time cap [s] (optional)
//     CUTOFF_SAMPLING         VERTICAL | ISOTROPIC
//     DIRECTIONAL_MAP         T|F        ! enable directional cutoff sky-map output
//     DIRMAP_LON_RES          <double>   ! longitude resolution [deg] for sky-map
//     DIRMAP_LAT_RES          <double>   ! latitude resolution [deg] for sky-map
//
//   #DENSITY_SPECTRUM
//     DS_EMIN                 <double>   ! MeV/n; lower energy bound
//     DS_EMAX                 <double>   ! MeV/n; upper energy bound
//     DS_NINTERVALS           <int>      ! number of energy intervals (nPoints = +1)
//     DS_ENERGY_SPACING       LOG | LINEAR
//     DS_MAX_PARTICLES        <int>      ! total trajectory cap per obs. point (optional)
//     DS_MAX_TRAJ_TIME        <double>   ! per-trajectory time cap [s] (optional)
//     DS_BOUNDARY_MODE        ISOTROPIC | ANISOTROPIC
//       ISOTROPIC (default): T(E;x0) = N_allowed/N_dirs, uniform boundary spectrum
//       ANISOTROPIC:         T_aniso(E;x0) = (1/N_dirs)*sum_k A_k*f_PAD_k*f_spatial_k
//                            requires #BOUNDARY_ANISOTROPY section
//
//   #BOUNDARY_ANISOTROPY     (required when DS_BOUNDARY_MODE = ANISOTROPIC)
//     BA_PAD_MODEL            ISOTROPIC | SINALPHA_N | COSALPHA_N | BIDIRECTIONAL
//     BA_PAD_EXPONENT         <double>   ! n in sin^n or |cos|^n  (default 2.0)
//     BA_SPATIAL_MODEL        UNIFORM | DAYSIDE_NIGHTSIDE
//     BA_DAYSIDE_FACTOR       <double>   ! flux multiplier for GSM x > 0  (default 1.0)
//     BA_NIGHTSIDE_FACTOR     <double>   ! flux multiplier for GSM x <= 0 (default 1.0)
//
//   #PARTICLE_SPECIES
//     SPECIES_NAME            <string>   ! e.g. PROTON, ELECTRON, HE4
//     SPECIES_CHARGE          <int>      ! charge in units of e (sign matters)
//     SPECIES_MASS_AMU        <double>   ! mass in atomic mass units
//
//   #BACKGROUND_FIELD
//     FIELD_MODEL             T96 | T05 | DIPOLE
//     EPOCH                   <ISO datetime>   ! e.g. 2003-11-20T06:00
//     DST                     <double>   ! nT
//     PDYN                    <double>   ! nPa
//     IMF_BY                  <double>   ! nT (GSM Y component of IMF)
//     IMF_BZ                  <double>   ! nT (GSM Z component of IMF)
//     IMF_BX                  <double>   ! nT (reserved; not used by T96/T05)
//     SW_VX                   <double>   ! km/s solar wind x velocity (T05 input)
//     SW_N                    <double>   ! cm^-3 solar wind number density (T05 input)
//     T05_W1 .. T05_W6        <double>   ! T05 storm-time history integrals W1..W6
//     DIPOLE_MOMENT           <double>   ! multiple of Earth dipole moment M_E (DIPOLE model)
//     DIPOLE_TILT             <double>   ! tilt from +Z_GSM toward +X_GSM [deg] (DIPOLE model)
//
//   #DOMAIN_BOUNDARY
//     DOMAIN_XMIN             <double>   ! km (GSM)
//     DOMAIN_XMAX             <double>   ! km
//     DOMAIN_YMIN             <double>   ! km
//     DOMAIN_YMAX             <double>   ! km
//     DOMAIN_ZMIN             <double>   ! km
//     DOMAIN_ZMAX             <double>   ! km
//     R_INNER                 <double>   ! km inner loss sphere radius
//
//   #OUTPUT_DOMAIN
//     OUTPUT_MODE             POINTS | SHELLS
//     COORDS                  GSM | GEO | GSE      ! coordinate label (not transformed)
//     SHELL_ALTITUDES         <double> [<double> ...] ! km above Earth's surface; space-separated
//     SHELL_RES               <double>   ! angular resolution [deg] for lon/lat grid
//     POINTS_BEGIN
//       x1 y1 z1
//       x2 y2 z2
//       ...
//     POINTS_END
//     (coordinates in km, GSM by default)
//
//   #NUMERICAL
//     DT_TRACE                <double>   ! initial time step [s]
//     MAX_STEPS               <int>      ! hard cap on integration steps
//     MAX_TRACE_TIME          <double>   ! hard cap on integration time [s]
//     MAX_TRACE_DISTANCE      <double>   ! hard cap on cumulative trace distance [Re]
//                                        ! 0 or negative => disabled
//
//   #SPECTRUM
//     (stored both as a raw key/value map in AmpsParam.spectrum and as a
//      typed metadata view in AmpsParam.particleSpectrum; the downstream
//      evaluator still uses boundary/spectrum.h InitGlobalSpectrumFromKeyValueMap)
//
//   Any unrecognised section or key is stored in AmpsParam.unknown for diagnostics.
//
//   #ENERGY_CHANNELS    (optional; if present, enables per-channel integral flux output)
//     CH_BEGIN
//       NAME   E1_MeV   E2_MeV   [! optional comment]
//       ...
//     CH_END
//     Example:
//       CH_BEGIN
//         CH1   10.0   100.0    ! 10–100 MeV
//         CH2  100.0  1000.0    ! 100 MeV–1 GeV
//         CH3 1000.0 10000.0    ! 1–10 GeV
//         CH4 10000.0 20000.0   ! 10–20 GeV
//       CH_END
//     Channels may overlap and may extend outside [DS_EMIN, DS_EMAX].
//     If this section is absent, only the total integral flux is computed.
//
//======================================================================================
// TYPE CONVERSION RULES
//======================================================================================
//
// Booleans:     T / TRUE / 1  -> true     F / FALSE / 0 -> false   (case-insensitive)
// Integers:     std::stoi     (raises std::invalid_argument on failure)
// Doubles:      std::stod     (raises std::invalid_argument on failure)
// Strings:      trimmed of leading/trailing whitespace; comment text after '!' removed
// Enumerations: stored as uppercase strings; validated lazily by the consumer
//
//======================================================================================
// UNITS CONTRACT
//======================================================================================
//
// The parser stores distances in whatever unit the file uses (km for Runs-on-Request
// inputs; Re for some legacy files). The dominant convention used in production is km.
// The solvers convert km -> m internally.
//
// IMPORTANT: The parser does NOT transform coordinates between GSM, GEO, and GSE.
// The COORDS keyword is stored as a label only. The gridless solver assumes all
// positions are in GSM and will produce incorrect results if non-GSM coordinates are
// passed without external conversion.
//
//======================================================================================
// ERROR HANDLING
//======================================================================================
//
// ParseAmpsParamFile throws std::runtime_error for:
//   - File not found / cannot open
//   - Malformed numeric values for recognised keys (std::stod/stoi failure)
//   - POINTS_END before POINTS_BEGIN
//
// It does NOT throw for:
//   - Unknown section names (silently skipped)
//   - Unknown key names within a recognised section (stored in raw map)
//   - Missing required sections (defaults in structs are used)
//
// Post-parse validation (e.g., checking that #BOUNDARY_ANISOTROPY is present when
// DS_BOUNDARY_MODE = ANISOTROPIC) is performed in amps_param_parser.cpp after the
// parse loop, and also at solver startup in DensityGridless.cpp.
//
//======================================================================================

#ifndef _SRC_EARTH_UTIL_AMPS_PARAM_PARSER_H_
#define _SRC_EARTH_UTIL_AMPS_PARAM_PARSER_H_

#include <string>
#include <vector>
#include <map>
#include <stdexcept>

namespace EarthUtil {

  // Simple 3-component container used throughout.
  struct Vec3 {
    double x{0.0}, y{0.0}, z{0.0};
  };

  //====================================================================================
  // ParticleSpectrum / SpectrumTablePoint
  //====================================================================================
  // A typed representation of the #SPECTRUM block.  The raw key/value table is kept
  // as well, but this structure makes the most important spectrum choices available
  // during initialization without requiring downstream code to reinterpret strings.
  //
  // Design goal:
  //   - keep the parser independent from the actual injection implementation;
  //   - validate the commonly used POWER_LAW inputs early;
  //   - preserve all original raw keys for forward compatibility.
  struct SpectrumTablePoint {
    double energy_MeV{0.0};
    double intensity{0.0};
  };

  struct ParticleSpectrum {
    enum class Type { UNKNOWN, POWER_LAW, POWER_LAW_CUTOFF, LIS_FORCE_FIELD, BAND, TABLE };
    enum class TableFormat { UNKNOWN, TWO_COLUMN, TIME_DEPENDENT };

    Type type{Type::UNKNOWN};
    std::string typeName{"UNKNOWN"};
    bool parsed{false};
    TableFormat tableFormat{TableFormat::UNKNOWN};

    // Common power-law style parameters.  These are populated when present and are
    // required for POWER_LAW.  For other spectrum types they remain optional metadata
    // unless a downstream module chooses to use them.
    double J0{0.0};
    double gamma{0.0};
    double E0_MeV{0.0};
    double Emin_MeV{0.0};
    double Emax_MeV{0.0};

    // Optional extra parameters for extended spectrum forms.  They are intentionally
    // generic because the exact downstream evaluator lives outside this parser.
    double cutoffEnergy_MeV{0.0};
    double modulationPotential_MV{0.0};
    double alpha{0.0};
    double beta{0.0};
    double breakEnergy_MeV{0.0};
    std::string tableFile;
    std::string resolvedTableFile;
    std::string selectedSnapshotEpochUTC;
    double selectedSnapshotTimeOffset_s{0.0};
    std::vector<SpectrumTablePoint> tablePoints;

    std::map<std::string,std::string> raw;
  };

  struct DomainBox {
    // NOTE ON UNITS
    //   For Runs-on-Request inputs we treat *all geometric distances* as being
    //   provided in **kilometers** (km):
    //     - DOMAIN_* bounds
    //     - R_INNER
    //     - POINT coordinates in POINTS_BEGIN..END
    //   The gridless cutoff solver converts km -> meters -> Re internally.
    double xMin{-60.0}, xMax{15.0};
    double yMin{-25.0}, yMax{25.0};
    double zMin{-20.0}, zMax{20.0};
    double rInner{2.0}; // [km] inner loss sphere radius
  };

  struct CutoffScan {
    double eMin_MeV{1.0};
    double eMax_MeV{1000.0};
    int nEnergy{50};
    int maxParticlesPerPoint{500};

    // Optional cap on how long we integrate a *single* backtraced trajectory
    // during the cutoff scan.
    //
    // Why this exists:
    //   Cutoff rigidity searches can spend a lot of time on quasi-trapped
    //   or long-lived trajectories. For interactive/production runs you may
    //   want a tighter limit than the global #NUMERICAL MAX_TRACE_TIME.
    //
    // Semantics (kept consistent with other sections like #DENSITY_SPECTRUM):
    //   - If maxTrajTime_s > 0: use it as the cutoff-scan trace time cap [s].
    //   - If maxTrajTime_s <= 0 or omitted: fall back to #NUMERICAL
    //     MAX_TRACE_TIME.
    double maxTrajTime_s{0.0}; // CUTOFF_MAX_TRAJ_TIME

    // Cutoff sampling mode.
    //
    // VERTICAL:
    //   Compute cutoff using only the local "vertical" arrival direction.
    //   In our convention, this is the direction pointing *toward Earth*:
    //     d_vertical = - unit(r0)
    //   where r0 is the GSM position vector of the injection point.
    //
    // ISOTROPIC:
    //   Compute an isotropic cutoff by scanning many arrival directions
    //   (a pre-defined direction grid) and taking the minimum Rc.
    //
    // IMPORTANT NOTE ABOUT DEFINITIONS:
    //   "Isotropic" in this implementation means "min over sampled sky",
    //   not "effective cutoff" based on allowed fraction / penumbra.
    std::string sampling{"ISOTROPIC"}; // CUTOFF_SAMPLING

    // Optional: compute a directional cutoff rigidity "sky-map" for each
    // injection point.
    //
    // When enabled, the solver evaluates Rc as a function of arrival
    // direction on a lon/lat grid in a *local, point-centered coordinate
    // system* (see detailed documentation in CutoffRigidityGridless.cpp).
    //
    // The intent is to provide a diagnostic product to visualize the
    // directional dependence (penumbra-like structure) and to support
    // scientific interpretation.
    bool directionalMap{false};          // DIRECTIONAL_MAP
    double dirMapLonRes_deg{10.0};       // DIRMAP_LON_RES
    double dirMapLatRes_deg{10.0};       // DIRMAP_LAT_RES
  };

  //====================================================================================
  // Density + spectrum sampling controls (gridless)
  //====================================================================================
  // The gridless density/spectrum workflow (CALC_TARGET = DENSITY_SPECTRUM) uses a
  // dedicated input section:
  //   #DENSITY_SPECTRUM
  //     DS_EMIN           <min energy>      ! MeV/n
  //     DS_EMAX           <max energy>      ! MeV/n
  //     DS_NINTERVALS     <n>               ! number of energy *intervals*
  //     DS_ENERGY_SPACING LOG|LINEAR
  //
  // Notes:
  // - DS_NINTERVALS is stored as "intervals" (not points) because this is the
  //   most robust way to define a grid: Npoints = Nintervals + 1.
  // - Energies are kinetic energy. For PROTON this is MeV per particle.
  //   For ions, MeV/n is commonly used; we keep the unit label but do not
  //   apply any per-nucleon conversion inside the parser.
  struct DensitySpectrumParam {
    double Emin_MeV{1.0};
    double Emax_MeV{1000.0};
    int nIntervals{50};

    // Optional cap on trajectory work per observation point.
    //
    // DS_MAX_PARTICLES limits the *total* number of backtraced trajectories
    // launched from a single observation point across the entire energy grid.
    // In DensityGridless this is enforced by reducing the number of sampled
    // directions per energy:
    //   nDir(E) = min(nDirDefault, floor(DS_MAX_PARTICLES / NenergyPoints)).
    //
    // If omitted or <= 0, no cap is applied.
    int maxParticlesPerPoint{0}; // DS_MAX_PARTICLES

    // Optional cap on how long we integrate a *single* trajectory.
    //
    // DS_MAX_TRAJ_TIME [s] provides an additional (often tighter) time limit
    // on the backtracing integration used to classify ALLOWED/FORBIDDEN.
    // This is useful for density/spectrum calculations because they can trace
    // many more trajectories than a cutoff scan.
    //
    // Semantics:
    //   - If DS_MAX_TRAJ_TIME > 0: use it as the per-trajectory integration
    //     time limit (in seconds).
    //   - If DS_MAX_TRAJ_TIME <= 0 or omitted: fall back to #NUMERICAL
    //     MAX_TRACE_TIME.
    double maxTrajTime_s{0.0}; // DS_MAX_TRAJ_TIME

    enum class Spacing { LOG, LINEAR };
    Spacing spacing{Spacing::LOG};

    int nPoints() const { return (nIntervals > 0) ? (nIntervals + 1) : 0; }

    // DS_BOUNDARY_MODE selects the density/spectrum solver branch.
    //
    // ISOTROPIC   (default, backward-compatible):
    //   T(E; x0) = N_allowed / N_dirs
    //   J_loc(E; x0) = T(E; x0) * J_b(E)
    //   The boundary spectrum J_b is assumed uniform and isotropic; the exit
    //   position and direction of each allowed trajectory are discarded.
    //
    // ANISOTROPIC:
    //   T_aniso(E; x0) = (1/N_dirs) * sum_k [ A_k * f_PAD(cos_alpha_k) * f_spatial(x_k) ]
    //   J_loc(E; x0) = J_b_iso(E) * T_aniso(E; x0)
    //   where cos_alpha_k = v_exit_k . B_hat(x_exit_k) and x_exit_k is the
    //   GSM position where trajectory k crossed the outer domain boundary.
    //   The PAD and spatial modulation models are controlled by #BOUNDARY_ANISOTROPY.
    //   Requires TraceAllowedSharedEx() to return exit state per trajectory.
    std::string boundaryMode{"ISOTROPIC"}; // DS_BOUNDARY_MODE
  };

  //====================================================================================
  // Anisotropic boundary spectrum parameters (#BOUNDARY_ANISOTROPY section)
  //====================================================================================
  // These parameters control the non-isotropic boundary spectrum used when
  // DS_BOUNDARY_MODE = ANISOTROPIC.
  //
  // The full boundary intensity is factored as:
  //   J_b(E, Omega, x) = J_b_iso(E) * f_PAD(cos_alpha) * f_spatial(x)
  //
  // where:
  //   cos_alpha = v_exit . B_hat(x_exit)   (pitch angle at the domain boundary)
  //   x_exit                               (GSM exit position on the outer boundary)
  //
  // PAD MODELS (BA_PAD_MODEL):
  //   ISOTROPIC      f(alpha) = 1                         (reduces to isotropic)
  //   SINALPHA_N     f(alpha) = sin^n(alpha)              (pancake distribution)
  //   COSALPHA_N     f(alpha) = |cos(alpha)|^n            (field-aligned beam)
  //   BIDIRECTIONAL  f(alpha) = |cos(alpha)|^n            (symmetric about equator;
  //                                                        identical to COSALPHA_N
  //                                                        but documents intent)
  //
  // SPATIAL MODELS (BA_SPATIAL_MODEL):
  //   UNIFORM              f_spatial = 1 everywhere
  //   DAYSIDE_NIGHTSIDE    f_spatial = BA_DAYSIDE_FACTOR  if GSM x > 0
  //                                  = BA_NIGHTSIDE_FACTOR if GSM x <= 0
  //====================================================================================
  struct AnisotropyParam {
    // Pitch angle distribution
    std::string padModel{"ISOTROPIC"};   // BA_PAD_MODEL
    double padExponent{2.0};             // BA_PAD_EXPONENT  (n in sin^n or |cos|^n)

    // Spatial flux modulation
    std::string spatialModel{"UNIFORM"}; // BA_SPATIAL_MODEL
    double daysideFactor{1.0};           // BA_DAYSIDE_FACTOR   (GSM x > 0 multiplier)
    double nightsideFactor{1.0};         // BA_NIGHTSIDE_FACTOR (GSM x <= 0 multiplier)
  };

  struct Species {
    std::string name{"PROTON"};
    int charge_e{1};
    double mass_amu{1.0};
  };

  struct BackgroundField {
    // FIELD_MODEL selects the background magnetic field model evaluated by
    // the gridless tools. Historically only Tsyganenko models were supported
    // (T96/T05). We extend this with an analytic dipole for verification and
    // regression tests.
    //
    // Supported values (aliases are normalised to canonical form by the parser):
    //   Canonical  Accepted aliases        Description
    //   --------   ------------------      ----------------------------------------
    //   "T96"      "TS96", "T96S"          Tsyganenko (1996)
    //   "T05"      "TS05","T05S","T04S",   Tsyganenko-Sitnov (2005)
    //              "TS04"
    //   "T01"      "TS01", "T01S"          Tsyganenko (2001)
    //   "TA15N"                            Tsyganenko-Andreeva (2015) northward
    //   "TA15B"                            Tsyganenko-Andreeva (2015) southward
    //   "TA16"    "TA16RBF"               Tsyganenko-Andreeva (2016)
    //   "DIPOLE"                           Analytic centered dipole (internal only)
    std::string model{"T96"};

    // --- Dipole-only parameters (FIELD_MODEL = DIPOLE) ---
    // DIPOLE_MOMENT: dipole moment magnitude as a multiple of Earth's canonical
    // dipole moment M_E (default 1.0).
    double dipoleMoment_Me{1.0}; // keyword: DIPOLE_MOMENT

    // DIPOLE_TILT [deg]: dipole tilt angle in degrees measured from +Z_GSM
    // toward +X_GSM (rotation about +Y_GSM). Default 0.0 aligns the dipole
    // with the GSM Z-axis.
    double dipoleTilt_deg{0.0};  // keyword: DIPOLE_TILT (alias: DIPOLE_TILT_DEG)

    // Common parameters
    double dst_nT{-50.0};
    double pdyn_nPa{2.0};
    double imfBy_nT{0.0};
    double imfBz_nT{0.0};

    // Present in RoR files (reserved)
    double imfBx_nT{0.0};
    double swVx_kms{-400.0};
    double swN_cm3{5.0};

    // T01 storm-time activity indices
    double g[3]{0,0,0};

    // T05/TA15/TA16 storm-time history integrals
    double w[6]{0,0,0,0,0,0};

    // TA15 running averages of southward IMF Bz magnitude
    double bzAvg[6]{0,0,0,0,0,0};

    // TA15/TA16 optimal solar-wind coupling index.
    // For TA15 this maps to PARMOD(4); for TA16 it maps to PARMOD(3).
    // Defaults to 0.0 when not provided.
    double xind{0.0};

    // TA16 coefficient file path (TA16_RBF.par or equivalent).
    // Empty string means use the Fortran default (TA16_RBF.par in the CWD).
    // Set via keyword TA16_COEFF_FILE in #BACKGROUND_FIELD.
    std::string ta16CoeffFile;

    // Snapshot time
    std::string epoch{"2000-01-01T00:00"};

    // Optional external driver table file specified directly in #BACKGROUND_FIELD
    // via DRIVER_FILE.
    //
    // Why this lives in BackgroundField instead of Temporal:
    //   Historically, AMPS inputs often describe the magnetic-field model and its
    //   driver source together in one section:
    //
    //     #BACKGROUND_FIELD
    //     FIELD_MODEL   TS05
    //     DRIVER_FILE   ts05_driving_....txt
    //
    //   The gridless implementation already has a generic time-dependent driver
    //   table capability, but before this change it was only activated through
    //   #TEMPORAL TS_INPUT_MODE=FILE / TS_INPUT_FILE=...
    //
    //   Storing the filename here lets the parser honor the more natural
    //   BACKGROUND_FIELD-based syntax while still loading the same shared
    //   TsDriverTable used by the trajectory/gridless solvers.
    //
    // Semantics:
    //   - empty string: no driver table requested from #BACKGROUND_FIELD
    //   - non-empty    : parser must load the file after the full input is read
    //                    and attach it to temporal.driverTable
    //
    // Precedence rules are implemented in amps_param_parser.cpp.
    std::string driverFile;

    // Raw key/value store for forward compatibility.
    std::map<std::string,std::string> raw;
  };


  struct ElectricField {
    // Extensible electric-field model selector for the AMPS-backed 3D path.
    //
    // Supported values in this patch:
    //   - NONE
    //   - COROTATION_VOLLAND_STERN
    //
    // The intent is to keep the parser/API stable while making it easy to
    // add other models later (Weimer, user tables, solver-backed fields, etc.).
    std::string model{"NONE"};

    // Multiplier applied to the nominal corotation electric field.
    double corotationScale{1.0};

    // Volland-Stern-type shielded convection potential parameters. The exact
    // volumetric extension used by the 3D initializer is documented in
    // srcEarth/3d/ElectricField.cpp.
    double vsPotential_kV{60.0};
    double vsGamma{2.0};
    double vsReferenceL{10.0};
    double vsScale{1.0};

    // Numerical guards used by the volumetric initializer.
    double rMin_km{20.0};
    double lMin{1.0e-3};

    std::map<std::string,std::string> raw;
  };


  //====================================================================================
  // SpacecraftTrajectoryPoint / SpacecraftTrajectory
  //====================================================================================
  // To unify the treatment of standalone POINTS and sampled spacecraft trajectories,
  // every spatial sample is represented as a trajectory point carrying:
  //   - timeUTC    : ISO-8601 timestamp string associated with the sample
  //   - xGSM_m     : location in GSM Cartesian coordinates [m]
  //
  // UNIT CONTRACT
  //   The position is stored in SI meters so it can be handed directly to the
  //   particle-tracing routines (Boris pusher, field evaluators) without any
  //   additional unit conversion at the call site.
  //
  //   The legacy OutputDomain::points array (used for Tecplot output and the
  //   inherited POINTS code path) remains in km; RebuildFlattenedPointsFromTrajectories
  //   performs the m -> km division when populating it.
  //
  // For legacy OUTPUT_MODE=POINTS runs, we create one synthetic trajectory per point.
  // Each such trajectory contains exactly one sample and therefore reuses the exact
  // same downstream computation kernel as true multi-point trajectories.
  //
  // IMPORTANT:
  //   The solver physics is always evaluated in GSM. If the trajectory file is given
  //   in another frame (for example GEO or SM), the parser converts it to GSM when
  //   loading the file so all downstream code sees a single consistent representation.
  struct SpacecraftTrajectoryPoint {
    std::string timeUTC;
    Vec3 xGSM_m;   // GSM Cartesian position [m]
  };

  class SpacecraftTrajectory {
  public:
    std::string name;
    std::string sourceFrame{"GSM"};
    std::vector<SpacecraftTrajectoryPoint> samples;

    inline bool empty() const { return samples.empty(); }
    inline std::size_t size() const { return samples.size(); }

    inline void clear() {
      name.clear();
      sourceFrame = "GSM";
      samples.clear();
    }

    inline void AddSample(const std::string& timeUTC, const Vec3& xGSM_m) {
      samples.push_back(SpacecraftTrajectoryPoint{timeUTC, xGSM_m});
    }
  };

  struct OutputDomain {
    // OUTPUT_MODE:
    //   POINTS     - explicit list of individual locations
    //   SHELLS     - spherical shell map(s)
    //   TRAJECTORY - time-ordered sequence loaded from TRAJ_FILE
    std::string mode{"POINTS"};

    // Coordinate label for reporting / backward compatibility with older inputs.
    // For trajectories, TRAJ_FRAME is the authoritative input-frame selector.
    std::string coords{"GSM"};

    // Legacy flattened point list used by the existing solver kernels.
    // For OUTPUT_MODE=TRAJECTORY this is populated automatically from trajectories[0].
    std::vector<Vec3> points;

    // Unified representation: each standalone point is packed into a one-sample
    // trajectory so trajectory and point workflows can share the same code path.
    std::vector<SpacecraftTrajectory> trajectories;

    // Trajectory-specific controls.
    std::string trajFrame{"GSM"};      // TRAJ_FRAME
    std::string trajFile;               // TRAJ_FILE
    double fluxDt_min{1.0};             // FLUX_DT  [min]

    // SHELLS: altitudes in km + resolution in degrees.
    std::vector<double> shellAlt_km;
    double shellRes_deg{15.0};

    // Helper utilities used by the parser and the solvers.
    inline bool HasTrajectorySamples() const {
      for (const auto& tr : trajectories) if (!tr.empty()) return true;
      return false;
    }

    inline void ClearFlattenedPoints() { points.clear(); }

    inline void RebuildFlattenedPointsFromTrajectories() {
      points.clear();
      for (const auto& tr : trajectories) {
        for (const auto& s : tr.samples) {
          // xGSM_m is stored in meters; the legacy points array is in km.
          points.push_back({s.xGSM_m.x / 1000.0,
                            s.xGSM_m.y / 1000.0,
                            s.xGSM_m.z / 1000.0});
        }
      }
    }

    std::map<std::string,std::string> raw;
  };

  struct Numerical {
    double dtTrace_s{1.0};
    int maxSteps{300000};
    double maxTraceTime_s{7200.0};

    // MAX_TRACE_DISTANCE [Re]
    // -----------------------
    // Optional hard cap on the *cumulative path length* traveled by a single
    // backtraced particle trajectory, measured in Earth radii.
    //
    // Why cumulative path length instead of distance-from-launch?
    //   The tracing workflows here (cutoff, density, and anisotropic spectrum)
    //   can encounter quasi-trapped or drifting trajectories that remain near the
    //   launch location while still traveling a very long total arc length.
    //   Using cumulative path length makes this control analogous to
    //   MAX_TRACE_TIME: it limits total tracing work irrespective of whether the
    //   orbit loops back near the start point.
    //
    // Semantics:
    //   maxTraceDistance_Re <= 0 : disabled (backward-compatible default)
    //   maxTraceDistance_Re >  0 : stop tracing once the accumulated segment
    //                              lengths exceed this threshold.
    //
    // Units:
    //   Earth radii (Re), not km and not meters. The solver converts Re -> m
    //   internally at runtime using _EARTH__RADIUS_.
    double maxTraceDistance_Re{0.0};
  };

  struct CalcMode {
    std::string target{"CUTOFF_RIGIDITY"};
    std::string fieldEvalMethod{"GRIDLESS"};
  };

  //====================================================================================
  //====================================================================================
  // Density3DParam — 3D volumetric density sampling for Mode3DForward (#DENSITY_3D)
  //====================================================================================
  //
  // Controls the energy grid used to bin sampled particle number density in the
  // volumetric (per-cell) sampling pass of Mode3DForward.
  //
  // Input file section:
  //   #DENSITY_3D
  //     DENS_EMIN              1.0          ! MeV/n  lower energy bound
  //     DENS_EMAX              20000.0      ! MeV/n  upper energy bound
  //     DENS_NENERGY           30           ! number of energy bins
  //     DENS_ENERGY_SPACING    LOG          ! LOG or LINEAR bin spacing
  //
  // Resulting output variable:
  //   n(E_i) [m^-3 J^-1] — differential number density per Joule in each cell
  //   Integrated over ΔE_i to get the bin-total number density [m^-3].
  //
  struct Density3DParam {
    double Emin_MeV{1.0};      // DENS_EMIN [MeV/n]
    double Emax_MeV{20000.0};  // DENS_EMAX [MeV/n]
    int    nEnergyBins{30};    // DENS_NENERGY

    // Energy-bin spacing (shared enum shape with DensitySpectrumParam::Spacing
    // to reuse the ParseEnergySpacingToken helper in the parser).
    enum class Spacing { LOG, LINEAR };
    Spacing spacing{Spacing::LOG};  // DENS_ENERGY_SPACING
  };

  //====================================================================================
  // Mode3DForwardOptions — runtime controls for the 3D forward particle mode
  //====================================================================================
  //
  // These options are populated from the CLI and/or the input file (where applicable).
  //
  // CLI flags:
  //   -mode 3d_forward                              select this mode
  //   -mode3d-output-initialized                    write initialized mesh Tecplot file
  //   -forward-niter <int>                          override iteration count
  //   -forward-nparticles <int>                     override simulation particles per iteration;
  //                                                 physical weight is recomputed automatically:
  //                                                   W = (π×∫J dE×A_boundary×dt)/N
  //   -forward-boundary-dist <ISOTROPIC|...>        override boundary distribution type
  //
  // Input file keywords (parsed from #NUMERICAL for max-step reuse):
  //   FORWARD_N_PARTICLES  <int>   simulation particles injected per iteration
  //
  struct Mode3DForwardOptions {
    // Write amps_3dforward_initialized.data.dat after InitMeshFields().
    bool outputInitializedFile{false};

    // Total number of forward integration iterations to run.
    // Can be overridden by -forward-niter on the CLI.
    int nIterations{10000};

    // Simulation particles injected from the domain boundary each iteration.
    // The physical weight per particle is calculated from the boundary flux,
    // the time step dt, and this count:
    //   W = (π × ∫J dE × A_boundary × dt) / nParticlesPerIter
    int nParticlesPerIter{1000};

    // Angular distribution type at the domain boundary.
    //   ISOTROPIC (default) — uniform-intensity boundary with cos(θ)-weighted
    //                          directions (see BoundaryDistribution.h for physics).
    // Future values: COSALPHA_N, SINALPHA_N, DAYSIDE_NIGHTSIDE
    // Populated from -forward-boundary-dist CLI flag or the input file.
    std::string boundaryDistType{"ISOTROPIC"};

    // Energy proposal distribution used by Mode3DForward boundary injection.
    //   SPECTRUM     — legacy/default behavior: sample E from J(E)dE and use
    //                  nearly equal statistical weights.
    //   LOG_UNIFORM  — sample E uniformly in log(E) and apply an individual
    //                  particle statistical-weight correction q(E).  This improves
    //                  high-energy statistics for steep SEP spectra while keeping
    //                  the total physical source rate identical to the input
    //                  spectrum.
    //
    // The field is currently set by the CLI.  The parser also accepts optional
    // input-file keys FORWARD_INJECTION_ENERGY[_DISTRIBUTION] /
    // FORWARD_ENERGY_SAMPLING for forward compatibility with future input files.
    std::string injectionEnergyDistribution{"SPECTRUM"};
  };

  // Mode3DOptions — command-line controls specific to the PIC-backed 3D workflow
  //====================================================================================
  struct Mode3DOptions {
    // Write amps_3d_initialized.data.dat after InitMeshFields().
    // Default is false because this is a diagnostic file and can be large.
    bool outputInitializedFile{false};

    // Magnetic-field source used by the Mode3D cutoff tracer.
    // false: AMR interpolation from cell-centered data populated by InitMeshFields().
    // true : direct analytic/background evaluator call via EvaluateBackgroundMagneticFieldSI().
    bool forceAnalyticMagneticField{false};
  };

  //====================================================================================
  // EnergyChannel — one user-defined integral-flux channel (#ENERGY_CHANNELS section)
  //====================================================================================
  //
  // PHYSICS
  //   An energy channel [E1, E2] defines the integration range for the omnidirectional
  //   integral particle flux at each observation point x0:
  //
  //     F_ch(x0) = 4π ∫_{E1}^{E2} J_loc(E; x0) dE
  //              = 4π ∫_{E1}^{E2} T(E; x0) · J_b(E) dE          [m^-2 s^-1]
  //
  //   Integral flux differs from number density by the absence of a 1/v(E) weight:
  //
  //     n(x0)    = 4π ∫ J_loc(E; x0) / v(E) dE                   [m^-3]
  //
  //   The 1/v weight makes density sensitive to slow particles. Flux treats all
  //   energies in the channel with equal weight (per unit energy interval).
  //
  // ANALYTIC CLOSED FORM (power-law boundary spectrum, gamma ≠ 1)
  //   For J_b(E) = J0 · (E/E0)^{-gamma}:
  //
  //     F_ch = 4π · T · J0 · E0^gamma / (gamma-1) · ( E1^{1-gamma} − E2^{1-gamma} )
  //
  //   Special case gamma = 2 (the test configuration):
  //     F_ch = 4π · T · J0 · E0² · ( 1/E1 − 1/E2 )
  //
  //   Total flux over the full solver energy range [Emin, Emax] is obtained by
  //   substituting E1=Emin, E2=Emax in the same formula.
  //
  // NUMERICAL IMPLEMENTATION
  //   Channel integrals are computed from the T(E) array already computed by the
  //   density solver. For each channel [E1, E2]:
  //
  //     (1) Find all grid nodes E_i with E1 ≤ E_i ≤ E2. Call this sub-array {E_i}.
  //     (2) Prepend E1 with T linearly interpolated from the two surrounding nodes if
  //         E1 does not coincide with a grid node. Append E2 similarly.
  //     (3) Accumulate F = 4π · Σ_i 0.5·(J_loc_i + J_loc_{i+1})·(E_{i+1} − E_i)
  //         using the trapezoidal rule in MeV (then convert to J internally, or use
  //         a consistent unit system throughout).
  //
  //   Channels may overlap and may extend outside [DS_EMIN, DS_EMAX]; the
  //   integration is silently clipped to the available grid range.
  //   No error is raised for an empty intersection (F_ch = 0 is returned).
  //
  // INPUT FORMAT (#ENERGY_CHANNELS section):
  //   CH_BEGIN
  //     NAME   E1_MeV   E2_MeV   [! optional comment]
  //     ...
  //   CH_END
  //
  //   NAME    — identifier string used in output column headings (F_NAME_m2s1).
  //             Allowed characters: letters, digits, underscores.
  //   E1_MeV  — lower channel boundary [MeV], must be > 0.
  //   E2_MeV  — upper channel boundary [MeV], must be > E1_MeV.
  //====================================================================================
  struct EnergyChannel {
    std::string name;       // identifier  (e.g. "CH1", "P10_100")
    double E1_MeV{0.0};    // lower bound [MeV]
    double E2_MeV{0.0};    // upper bound [MeV]
  };

  //====================================================================================
  // TsDriverRecord — one time-stamped snapshot of Tsyganenko model driving parameters
  //====================================================================================
  // Parsed from the Qin-Denton / ViRBO format (qd_sep2017.txt and equivalents).
  // Column layout (1-based, space-separated):
  //   $1  : ISO-8601 datetime  (e.g. 2017-09-10T00:00:00)
  //   $2–7: year month day hour min sec   (not read; $1 is used for time)
  //   $8  : IMF By  [nT]
  //   $9  : IMF Bz  [nT]
  //   $10 : Vsw     [km/s]
  //   $11 : Den_P   [cm^-3]
  //   $12 : Pdyn    [nPa]
  //   $13–15: G1 G2 G3   (T96 external-field G parameters, not presently used)
  //   $16–23: 8 status flags   (skipped)
  //   $24 : Kp      [dimensionless]
  //   $25 : akp3    [dimensionless]   (not presently used)
  //   $26 : Dst     [nT]
  //   $27–32: Bz1..Bz6  (historic average IMF Bz over 1..6 h; not presently used)
  //   $33–38: W1..W6    (T05 storm-time history integrals)
  //   $39–44: W_status flags  (skipped)
  //
  // Only the physically required fields (By, Bz, Vsw, DenP, Pdyn, Dst, W1..W6)
  // are stored; everything else is parsed and discarded.
  struct TsDriverRecord {
    double et{0.0};          // SPICE ephemeris time of this record [s]
    std::string timeUTC;     // original ISO-8601 string (kept for diagnostics)

    // T96 + T05 shared parameters
    double imfBy_nT{0.0};
    double imfBz_nT{0.0};
    double swVx_kms{-400.0}; // stored as positive magnitude in file; sign applied here
    double swN_cm3{5.0};
    double pdyn_nPa{2.0};
    double dst_nT{0.0};

    // T01 storm-time activity indices G1..G3
    double g[3]{0,0,0};

    // T05/TA15/TA16 storm-time integrals W1..W6
    double w[6]{0,0,0,0,0,0};

    // TA15 running averages BZ1..BZ6
    double bzAvg[6]{0,0,0,0,0,0};

    // TA15 coupling index XIND (optional in generic driver tables).
    double xind{0.0};
  };

  //====================================================================================
  // TsDriverTable — time-ordered table of TsDriverRecords with lookup/interpolation
  //====================================================================================
  //
  // USAGE
  //   After construction the table is populated by ParseAmpsParamFile when
  //   TS_INPUT_MODE = FILE and TS_INPUT_FILE is set in #TEMPORAL.
  //
  //   At runtime (per trajectory point), call:
  //     TsDriverRecord rec = table.Lookup(et);
  //   to get a linearly interpolated snapshot of every driving parameter at the
  //   requested ephemeris time.  The result can be applied directly to
  //   BackgroundField before (re-)initialising Geopack and the Tsyganenko model.
  //
  // INTERPOLATION POLICY
  //   Linear interpolation between the two bracketing records is used for all
  //   continuous parameters (By, Bz, Vsw, DenP, Pdyn, Dst, W1..W6).
  //   If the requested time lies outside the table range, the nearest endpoint
  //   record is returned (clamped) with a one-time stderr warning.
  //
  // THREAD SAFETY
  //   The table is read-only after construction.  Concurrent Lookup() calls from
  //   MPI workers (or OpenMP threads) are safe without any locking.
  class TsDriverTable {
  public:
    TsDriverTable() = default;

    // Returns true if the table contains at least one record.
    inline bool empty() const { return records_.empty(); }
    inline std::size_t size() const { return records_.size(); }

    // Add a record (must be called in chronological order).
    void push_back(const TsDriverRecord& r) { records_.push_back(r); }

    // Linear interpolation at ephemeris time et [s].
    // Clamps to the table endpoints when et is out of range.
    // Exits with a clear error message if the table is empty.
    TsDriverRecord Lookup(double et) const;

    // Apply the parameters from a record into a BackgroundField snapshot.
    // This does NOT update field.epoch (the caller is responsible for setting
    // epoch to the trajectory-point timestamp before calling Geopack::Init).
    static void ApplyToField(const TsDriverRecord& rec, BackgroundField& field);

  private:
    std::vector<TsDriverRecord> records_;
    mutable bool clampWarnedLow_{false};
    mutable bool clampWarnedHigh_{false};
  };



  // Build the 10-element Tsyganenko PARMOD vector plus trailing compatibility slot.
  // Layout depends on the canonical FIELD_MODEL name. Unknown or unsupported
  // models fall back to the shared base set [Pdyn,Dst,By,Bz,0..].
  void BuildTsParmod(const BackgroundField& field,
                     const std::string& model,
                     double parmod[11]);

  //====================================================================================
  // TemporalParam — parameters from the #TEMPORAL section
  //====================================================================================
  // Controls the time-series (trajectory) execution mode: event window, update
  // cadences, and the optional time-varying Tsyganenko driver file.
  struct TemporalParam {
    // TEMPORAL_MODE: TIME_SERIES or SNAPSHOT (default SNAPSHOT = legacy behaviour).
    std::string mode{"SNAPSHOT"};

    // Event window [ISO-8601 UTC strings].
    std::string eventStart;
    std::string eventEnd;

    // FIELD_UPDATE_DT [min]: cadence at which the magnetic field model is updated.
    double fieldUpdateDt_min{5.0};

    // INJECT_DT [min]: cadence at which particles are injected along the trajectory.
    double injectDt_min{30.0};

    // TS_INPUT_MODE: FILE (load from TS_INPUT_FILE) or PARAMS (use #BACKGROUND_FIELD).
    std::string tsInputMode{"PARAMS"};

    // TS_INPUT_FILE: path to a Qin-Denton / ViRBO formatted driver file.
    std::string tsInputFile;

    // Populated by ParseAmpsParamFile when tsInputMode == "FILE".
    TsDriverTable driverTable;
  };

  struct AmpsParam {
    std::string runId{"UNKNOWN"};

    CalcMode calc;
    Mode3DOptions mode3d;
    Mode3DForwardOptions mode3dForward;   // -mode 3d_forward controls
    CutoffScan cutoff;
    DensitySpectrumParam densitySpectrum;
    Density3DParam density3d;             // #DENSITY_3D section
    AnisotropyParam anisotropy;
    Species species;
    BackgroundField field;
    ElectricField efield;
    DomainBox domain;
    OutputDomain output;
    Numerical numerics;
    TemporalParam temporal;   // #TEMPORAL section (time-series mode + driver table)

    std::map<std::string,std::string> spectrum;
    ParticleSpectrum particleSpectrum;  // typed view of #SPECTRUM, built at init time
    std::map<std::string,std::string> outputOptions;

    // User-defined integral-flux channels.  Empty when #ENERGY_CHANNELS is absent;
    // in that case only the total integral flux F_tot is written to the output files.
    std::vector<EnergyChannel> fluxChannels;

    // Holds all unknown keys across sections (for diagnostics / forward compat).
    std::map<std::string,std::string> unknown;
  };

  // Parse and validate the #SPECTRUM section into a typed representation while
  // preserving the original raw key/value table for the downstream spectrum
  // evaluator.
  ParticleSpectrum ParseParticleSpectrum(const std::map<std::string,std::string>& kv);

  // Load a trajectory file and convert every sample into GSM Cartesian coordinates
  // in meters.  This is used by OUTPUT_MODE=TRAJECTORY during initialization.
  SpacecraftTrajectory LoadTrajectoryFileAsGsm(const std::string& fileName,
                                               const std::string& trajFrame);

  // Finalize the point-like output representation after parsing.
  //   POINTS     -> pack standalone points into one-sample synthetic trajectories
  //   TRAJECTORY -> load the trajectory file and flatten it into output.points
  //   SHELLS     -> leave the point/trajectory containers untouched
  void InitializePointLikeOutput(AmpsParam& p);

  // Parse an AMPS_PARAM file. Throws std::runtime_error on hard errors.
  AmpsParam ParseAmpsParamFile(const std::string& fileName);

  // Helper conversions (public because both CLI and solver use them).
  bool ToBool(const std::string& s);
  std::string ToUpper(std::string s);

}

#endif
