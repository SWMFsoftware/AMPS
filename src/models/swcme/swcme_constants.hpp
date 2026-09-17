#ifndef SWCME_CONSTANTS_HPP
#define SWCME_CONSTANTS_HPP

#include <cmath>

namespace swcme {
namespace constants {

// Astronomical conversion constants.
constexpr double PI = 3.1415926535897932384626433832795;
constexpr double AU_M = 1.495978707e11;       // IAU 2012 Resolution B2 [m]
constexpr double SOLAR_RADIUS_M = 6.957e8;    // IAU 2015 nominal radius [m]

// Fundamental physical constants. CODATA 2022 is the SWCME baseline.
constexpr double PROTON_MASS_KG = 1.67262192595e-27;
// CODATA 2022 alpha-particle mass.  Multi-species thermodynamic closure uses
// this value explicitly rather than the 4*m_p approximation so its mass
// density and sound-speed definitions are unambiguous and reproducible.
constexpr double ALPHA_PARTICLE_MASS_KG = 6.6446573450e-27;
constexpr double VACUUM_PERMEABILITY_N_A2 = 1.25663706127e-6;
constexpr double BOLTZMANN_J_K = 1.380649e-23;
constexpr double SPEED_OF_LIGHT_M_S = 299792458.0;
constexpr double ELEMENTARY_CHARGE_C = 1.602176634e-19;
constexpr double MEV_TO_J = 1.0e6 * ELEMENTARY_CHARGE_C;

// Model convention: Carrington sidereal solar rotation rate [rad/s].
constexpr double SOLAR_ROTATION_RAD_S = 2.86533e-6;

}  // namespace constants

namespace physics {

// Return the Alfvén speed for SI magnetic-field magnitude and mass density.
// Both SWCME models use this shared implementation so CFG02 can verify that
// their external-unit preparation paths reach the same dimensional result.
// The zero result preserves the models' existing handling of a non-positive
// field magnitude or density; the physical formula itself is unchanged.
inline double alfven_speed_m_s(double magnetic_field_T,
                               double mass_density_kg_m3) {
  return (magnetic_field_T > 0.0 && mass_density_kg_m3 > 0.0)
             ? magnetic_field_T /
                   std::sqrt(constants::VACUUM_PERMEABILITY_N_A2 *
                             mass_density_kg_m3)
             : 0.0;
}

}  // namespace physics
}  // namespace swcme

#endif
