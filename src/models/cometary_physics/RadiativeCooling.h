// Geometry-independent radiative-cooling laws for cometary neutral gases.
//
// The functions in this file evaluate published rates only.  Selecting the
// species, traversing AMPS particles, and removing internal energy belong to
// the comet application adapter, so the same physics can be used by DIM1,
// DIM2, and DIM3 geometries.

#ifndef _AMPS_COMETARY_PHYSICS_RADIATIVE_COOLING_H_
#define _AMPS_COMETARY_PHYSICS_RADIATIVE_COOLING_H_

#include <cmath>
#include <limits>

namespace CometaryPhysics {
namespace RadiativeCooling {

enum class Model {
  Disabled=0,

  // Analytic approximation to the H2O rotational cooling calculations of
  // Crovisier (1984) and Huebner (1985), as presented by Gombosi, Nagy, and
  // Cravens (1986), Rev. Geophys. 24, 667--700, equation (48b).
  Gombosi1986Water=1
};

// Return the positive energy-loss rate of one H2O molecule in J/s.
//
// Gombosi et al. write the volumetric cooling magnitude in erg cm^-3 s^-1:
//
//   |Q| = A T_rot^p n_H2O exp(-tau_IR),
//   tau_IR = 0.4 n_H2O r sigma_IR,
//   sigma_IR = 4.0e-15 cm^2,
//
// with (A,p)=(4.4e-22,3.35) below 52 K and
// (2.0e-20,2.47) at and above 52 K.  Dividing by n_H2O and converting
// erg to joules produces the per-particle rate used here.  Number density is
// supplied in m^-3 and cometocentric distance in m.
inline double Gombosi1986WaterEnergyLossRate(double rotationalTemperatureK,
    double waterNumberDensityM3,double cometocentricDistanceM) {
  if (!(rotationalTemperatureK>0.0) ||
      !(waterNumberDensityM3>0.0) ||
      !(cometocentricDistanceM>=0.0) ||
      !std::isfinite(rotationalTemperatureK) ||
      !std::isfinite(waterNumberDensityM3) ||
      !std::isfinite(cometocentricDistanceM)) return 0.0;

  const double densityCm3=waterNumberDensityM3*1.0e-6;
  const double radiusCm=cometocentricDistanceM*100.0;
  const double infraredAbsorptionCrossSectionCm2=4.0e-15;
  const double opticalDepth=0.4*densityCm3*radiusCm*
      infraredAbsorptionCrossSectionCm2;

  const double coefficient=(rotationalTemperatureK<52.0) ? 4.4e-22 : 2.0e-20;
  const double exponent=(rotationalTemperatureK<52.0) ? 3.35 : 2.47;
  const double rateErgPerSecond=coefficient*
      std::pow(rotationalTemperatureK,exponent)*std::exp(-opticalDepth);

  return rateErgPerSecond*1.0e-7;
}

inline double EnergyLossRate(Model model,double rotationalTemperatureK,
    double numberDensityM3,double cometocentricDistanceM) {
  switch (model) {
  case Model::Disabled:
    return 0.0;
  case Model::Gombosi1986Water:
    return Gombosi1986WaterEnergyLossRate(rotationalTemperatureK,
        numberDensityM3,cometocentricDistanceM);
  }

  return std::numeric_limits<double>::quiet_NaN();
}

} // namespace RadiativeCooling
} // namespace CometaryPhysics

#endif
