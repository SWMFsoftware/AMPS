#ifndef SEP_CORONAL_CME_CONSTANTS_H
#define SEP_CORONAL_CME_CONSTANTS_H

namespace SEP {
namespace CoronalCME {
namespace Constants {

// CODATA/IAU-compatible SI constants used by Stages 0--6.  They are named in
// one header so a test cannot unknowingly compare kernels built with different
// solar radii, gravitational parameters, or electromagnetic conventions.
constexpr double kPi = 3.141592653589793238462643383279502884;
constexpr double kSolarRadiusM = 6.957e8;
constexpr double kSolarGravitationalParameterM3PerS2 = 1.32712440018e20;
constexpr double kBoltzmannJPerK = 1.380649e-23;
constexpr double kVacuumPermeabilityHPerM = 1.25663706212e-6;
constexpr double kProtonMassKg = 1.67262192369e-27;
constexpr double kElectronMassKg = 9.1093837015e-31;
constexpr double kAlphaMassKg = 6.6446573357e-27;
constexpr double kElementaryChargeC = 1.602176634e-19;
constexpr double kSpeedOfLightMPerS = 299792458.0;

}  // namespace Constants
}  // namespace CoronalCME
}  // namespace SEP

#endif  // SEP_CORONAL_CME_CONSTANTS_H
