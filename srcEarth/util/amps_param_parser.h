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
  };

  struct Species {
    std::string name{"PROTON"};
    int charge_e{1};
    double mass_amu{1.0};
  };

  struct BackgroundField {
    // FIELD_MODEL: "T96" or "T05".
    std::string model{"T96"};

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
    Species species;
    BackgroundField field;
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
