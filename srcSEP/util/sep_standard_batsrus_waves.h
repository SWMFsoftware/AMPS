#ifndef SRCSEP_UTIL_SEP_STANDARD_BATSRUS_WAVES_H
#define SRCSEP_UTIL_SEP_STANDARD_BATSRUS_WAVES_H
#include "sep_common_header_path.h"
#include SRCSEP_SEP_COMMON_HEADER(sep_coefficient_physics.h)
#include <array>

namespace SEP { namespace Waves { namespace StandardBatsrus {

// USER-SELECTED standard BATSRUS energy export, not representative amplitudes
// or private major/minor arrays. Each entry is total (magnetic+kinetic) Alfven
// wave energy [J/m^3]. Wave1 propagates along actual B; Wave2 opposite actual B.
// Neither magnetic polarity nor geometric line order swaps these two labels.
struct RawTotalAlfvenEndpoint {
  std::array<double,2> waveJPerM3{{0.0,0.0}};
};
struct MagneticBranchView {
  double plusMagneticJPerM3=0.0,minusMagneticJPerM3=0.0;
  double plusVarianceT2=0.0,minusVarianceT2=0.0,totalVarianceT2=0.0;
  std::uint64_t generation=0,checksum=0;
  bool branchDirectionsKnown=false,circularPolarizationsKnown=false;
};
const char* ContractIdentity();
const char* Representation();

// Validate BOTH endpoints, even at fraction zero/one, then interpolate their
// private magnetic views with W_B=rawTotal/2. The only energy->deltaB^2 formula
// is the canonical magnetic helper. Positive conversion/interpolation underflow
// is explicit; no floors, extrapolation, endpoint masking or L/R invention.
// Inputs are immutable and output replacement is atomic on every failure.
Transport::Status InterpolateMagneticView(const RawTotalAlfvenEndpoint& begin,
    const RawTotalAlfvenEndpoint& end,double fraction,std::uint64_t generation,
    MagneticBranchView* output);

} } }
#endif
