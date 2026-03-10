//======================================================================================
// amps_param_parser.h
//======================================================================================
// PURPOSE
//   Lightweight parser for the AMPS_PARAM.in style inputs produced by the CCMC
//   Runs-on-Request interface for the Geospace energetic particle tools.
//
//   This parser is intentionally *self-contained* and does NOT depend on the legacy
//   srcEarth/parser.cpp logic. The goal is to support a standalone "gridless" cutoff
//   rigidity computation path that directly evaluates Tsyganenko field models (T96/T05)
//   without requiring the full PIC runtime configuration.
//
// FILE FORMAT OVERVIEW
//   The input examples provided by the user follow a simple sectioned format:
//     - Sections begin with a token like "#RUN_INFO" on its own line.
//     - Inside a section, each line is either:
//         KEY   VALUE   [! comment]
//       or a block delimiter:
//         POINTS_BEGIN / POINTS_END
//       or blank/comment lines.
//     - Lines beginning with '!' or '#' (outside a section header) are treated
//       as comments.
//
//   The parser is permissive:
//     - Unknown keys are stored into a "raw" map for later extension.
//     - Numeric values are parsed with std::stod/std::stoi.
//     - Boolean flags accept: T/F, TRUE/FALSE, 1/0.
//
// WHAT THIS PARSER POPULATES
//   The gridless cutoff solver currently needs only a subset of fields:
//     - CALCULATION_MODE: FIELD_EVAL_METHOD (GRIDLESS vs GRID_3D)
//     - CUTOFF_RIGIDITY: EMIN/EMAX, NENERGY, MAX_PARTICLES (future)
//     - PARTICLE_SPECIES: charge and mass
//     - BACKGROUND_FIELD: FIELD_MODEL (T96 or T05) + model parameters
//     - DOMAIN_BOUNDARY: bounding box and inner loss sphere
//     - OUTPUT_DOMAIN: POINTS or SHELLS definitions
//     - NUMERICAL: DT_TRACE and (optionally) max steps/stop criteria
//
// EXTENSION NOTES
//   The user intends to later extend this to compute fluxes. For that purpose
//   we keep:
//     - Spectrum section as a raw key/value map.
//     - Output options as a raw key/value map.
//
// THREAD SAFETY
//   This parser has no global state; it is safe to use from any thread.
//======================================================================================

#ifndef _SRC_EARTH_UTIL_AMPS_PARAM_PARSER_H_
#define _SRC_EARTH_UTIL_AMPS_PARAM_PARSER_H_

#include <string>
#include <vector>
#include <map>

namespace EarthUtil {

  // Simple 3-component container used throughout.
  struct Vec3 {
    double x{0.0}, y{0.0}, z{0.0};
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
    // Supported values:
    //   - "T96"    Tsyganenko (1996)
    //   - "T05"    Tsyganenko & Sitnov (2005)
    //   - "DIPOLE" Analytic centered dipole (internal field only)
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

    // T05 storm-time integrals
    double w[6]{0,0,0,0,0,0};

    // Snapshot time
    std::string epoch{"2000-01-01T00:00"};

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

  struct OutputDomain {
    // OUTPUT_MODE: "POINTS" or "SHELLS".
    std::string mode{"POINTS"};

    // Coordinate label for reporting.
    std::string coords{"GSM"};

    // POINTS: user-provided list.
    std::vector<Vec3> points;

    // SHELLS: altitudes in km + resolution in degrees.
    std::vector<double> shellAlt_km;
    double shellRes_deg{15.0};

    std::map<std::string,std::string> raw;
  };

  struct Numerical {
    double dtTrace_s{1.0};
    int maxSteps{300000};
    double maxTraceTime_s{7200.0};
  };

  struct CalcMode {
    std::string target{"CUTOFF_RIGIDITY"};
    std::string fieldEvalMethod{"GRIDLESS"};
  };

  struct AmpsParam {
    std::string runId{"UNKNOWN"};

    CalcMode calc;
    CutoffScan cutoff;
    DensitySpectrumParam densitySpectrum;
    Species species;
    BackgroundField field;
    ElectricField efield;
    DomainBox domain;
    OutputDomain output;
    Numerical numerics;

    // For later (flux) expansion
    std::map<std::string,std::string> spectrum;
    std::map<std::string,std::string> outputOptions;

    // Holds all unknown keys across sections.
    std::map<std::string,std::string> unknown;
  };

  // Parse an AMPS_PARAM file. Throws std::runtime_error on hard errors.
  AmpsParam ParseAmpsParamFile(const std::string& fileName);

  // Helper conversions (public because both CLI and solver use them).
  bool ToBool(const std::string& s);
  std::string ToUpper(std::string s);

}

#endif
