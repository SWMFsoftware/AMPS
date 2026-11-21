#ifndef _GENERIC_COMET_CHEMISTRY_VARIANTS_H_
#define _GENERIC_COMET_CHEMISTRY_VARIANTS_H_

// Rate tables of alternative reaction definitions ("chemistry variants").
// Design: srcHartley2/CHEMISTRY_VARIANTS.md.
//
// This module is deliberately free of AMPS/PIC/MPI dependencies so that the
// table reader, the interpolation, and the rate formulas can be tested by a
// stand-alone program (srcHartley2/tests/test_chemistry_variants.cpp). The
// runtime-input layer (CometRuntimeInput.cpp) calls it once per run, on every
// MPI rank, before the chemistry network is constructed.
//
// Units are SI unless a name says otherwise: velocities in km/s, the legacy
// OH lifetime table in its own unit (1e5 s before the build-time scale), rates
// in s^-1 at 1 AU. Failures return false and describe the problem, including
// the file name and line number, in *error; nothing is clamped silently except
// the documented velocity clamp, which is reported through cInterpolation.

#include <string>
#include <vector>

namespace Comet {
namespace ChemistryVariants {

// Lower-case hexadecimal SHA-256 (FIPS 180-4) of a byte string. Used to record
// the identity of the table file in the run log and the output headers.
std::string Sha256Hex(const std::string &bytes);

// Velocity-dependent OH lifetime table (model SA1982_OH_LIFETIME). File format,
// see srcHartley2/i3/reference/oh_lifetime_sa1982.csv:
//   '#' comment lines and blank lines are ignored;
//   first record:  tau_fuv_s,<FUV/EUV lifetime at 1 AU in s, > 0>
//   second record: v_h_km_s,ohtaus_1e5_s           (column header)
//   then >= 2 rows: <v_h in km/s>,<OHTAUS > 0>, v_h strictly increasing.
// The grid need not be uniform.
struct cLifetimeTable {
  std::string Path;          // as opened
  std::string Sha256;        // of the raw file bytes
  double TauFuvS=0.0;        // [s] at 1 AU
  std::vector<double> VelocityKmS;
  std::vector<double> Ohtaus;  // legacy table unit (multiply by scale S [s])
};

// Parse table text; name is used in error messages only.
bool ParseLifetimeTable(const std::string &text,const std::string &name,
    cLifetimeTable *table,std::string *error);

// Read and parse a table file and record its SHA-256.
bool ReadLifetimeTable(const std::string &path,cLifetimeTable *table,
    std::string *error);

struct cInterpolation {
  double Value=0.0;
  bool Clamped=false;  // v_h was outside [VelocityKmS.front(), VelocityKmS.back()]
};

// Linear interpolation in v_h; outside the table the end value is returned and
// Clamped is set.
cInterpolation InterpolateOhtaus(const cLifetimeTable &table,double velocityKmS);

// Rates of the SA1982_VH OH-loss variant at 1 AU:
//   k_d   = 1/TauFuvS + 1/(scaleS * OHTAUS(v_h))   [s^-1]
//   k_ion = ionFraction * k_d                       [s^-1]
// The caller applies the (1 AU / r_h)^2 scaling through the heliocentric rate
// law of the compiled reaction.
struct cSA1982Rates {
  double VelocityKmS=0.0;
  double Ohtaus=0.0;
  bool Clamped=false;
  double Dissociation1AU=0.0;
  double Ionization1AU=0.0;
};

bool EvaluateSA1982(const cLifetimeTable &table,double velocityKmS,
    double scaleS,double ionFraction,cSA1982Rates *rates,std::string *error);

}  // namespace ChemistryVariants
}  // namespace Comet

#endif
