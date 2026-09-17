#ifndef SWCME_UNITS_HPP
#define SWCME_UNITS_HPP

// ============================================================================
// swcme_units.hpp
// ----------------------------------------------------------------------------
// Centralized unit conversions for the SWCME public parameter interfaces.
//
// Why this file exists
// --------------------
// SWCME historically repeated conversion factors such as 1e3 (km/s -> m/s),
// 1e6 (cm^-3 -> m^-3), and 1e-9 (nT -> T) independently in the 1-D and 3-D
// implementations.  Even when the factors were numerically identical, the
// duplicated code allowed the two model paths to drift and made it difficult
// for CFG02 to verify a single dimensional contract.  All production code now
// uses the helpers below at the external-to-internal boundary.
//
// Design rule
// -----------
// These functions perform *unit conversion only*.  They deliberately do not
// clamp, validate, or otherwise change the physical value.  For example,
// 0 km/s converts exactly to 0 m/s and -1 km/s converts to -1000 m/s.  Whether
// such a value is physically admissible is the separate responsibility of
// swcme_config.hpp.  Keeping conversion and validation separate prevents the
// legacy failure in which the 1-D path silently converted 0 km/s to 1 m/s.
// ============================================================================

#include "swcme_constants.hpp"

namespace swcme {
namespace units {

constexpr double km_per_s_to_m_per_s(double value) { return value * 1.0e3; }
constexpr double m_per_s_to_km_per_s(double value) { return value * 1.0e-3; }

constexpr double cm3_to_m3(double value) { return value * 1.0e6; }
constexpr double m3_to_cm3(double value) { return value * 1.0e-6; }

constexpr double nT_to_T(double value) { return value * 1.0e-9; }
constexpr double T_to_nT(double value) { return value * 1.0e9; }

constexpr double km_inverse_to_m_inverse(double value) { return value * 1.0e-3; }
constexpr double m_inverse_to_km_inverse(double value) { return value * 1.0e3; }

constexpr double au_to_m(double value) { return value * constants::AU_M; }
constexpr double m_to_au(double value) { return value / constants::AU_M; }

constexpr double solar_radii_to_m(double value) {
  return value * constants::SOLAR_RADIUS_M;
}
constexpr double m_to_solar_radii(double value) {
  return value / constants::SOLAR_RADIUS_M;
}

constexpr double hours_to_seconds(double value) { return value * 3600.0; }
constexpr double seconds_to_hours(double value) { return value / 3600.0; }

constexpr double degrees_to_radians(double value) {
  return value * constants::PI / 180.0;
}
constexpr double radians_to_degrees(double value) {
  return value * 180.0 / constants::PI;
}

}  // namespace units
}  // namespace swcme

#endif
