//======================================================================================
// DipoleInterface.cpp
//======================================================================================
// IMPLEMENTATION NOTES
//   The dipole field is analytic and extremely cheap to evaluate. Dipole::gParams is
//   retained for backward-compatible serial callers. Step-3 field snapshots instead
//   own a Params value and use the explicit-parameter overload, so their supposedly
//   immutable state cannot be changed by another worker or provider construction.
//======================================================================================

#include "DipoleInterface.h"

namespace Earth {
namespace GridlessMode {
namespace Dipole {

Params gParams; // default-initialized (Earth-like, zero tilt)

void SetMomentScale(double momentScale_Me) {
  gParams.momentScale_Me = momentScale_Me;
}

void SetTiltDeg(double tilt_deg) {
  // Use the same constructor as immutable snapshots so both paths have exactly the
  // same tilt convention and floating-point operations.
  gParams=MakeParams(gParams.momentScale_Me,tilt_deg);
}

} // namespace Dipole
} // namespace GridlessMode
} // namespace Earth
