#pragma once

// ============================================================================
// swcme_solarwind.hpp
// ----------------------------------------------------------------------------
// Common analytical ambient-solar-wind physics shared by the SWCME 1-D and
// 3-D interfaces.
//
// This header intentionally contains the pieces that must not diverge between
// dimensional wrappers:
//   * Leblanc density normalization and evaluation,
//   * Parker radial/azimuthal field components,
//   * Cartesian Parker-vector construction for an arbitrary rotation axis,
//   * proton thermal pressure used by the current MHD shock model.
//
// Keeping these equations in one production component prevents the historical
// failure mode in which a correction was applied in 1-D but not in 3-D (or
// vice versa).  The dimensional wrappers remain responsible only for geometry
// that is genuinely different: a fixed-latitude 1-D ray versus a full 3-D
// position and solar-rotation axis.
// ============================================================================

#include "swcme_constants.hpp"

#include <array>
#include <algorithm>
#include <cmath>
#include <limits>

namespace swcme {
namespace solarwind {

// Published Leblanc, Dulk & Bougeret (1998) coefficients.  Radius is measured
// in nominal solar radii and the resulting electron number density is cm^-3.
// These immutable coefficients are normalized per model instance rather than
// modified globally, which is required for thread-safe simultaneous models.
constexpr double LEBLANC_A_CM3 = 3.3e5;
constexpr double LEBLANC_B_CM3 = 4.1e6;
constexpr double LEBLANC_C_CM3 = 8.0e7;

// The analytical SWCME background is defined only outside 1.05 solar radii.
// This constant is a DOMAIN BOUNDARY, not a clipping value.  Public model
// evaluators now return OUTSIDE_MODEL_DOMAIN for smaller radii.  The low-level
// formulas below assume the caller has already enforced that contract and do
// not silently replace an out-of-domain radius with MIN_RADIUS_M.
constexpr double MIN_RADIUS_M =
    1.05 * swcme::constants::SOLAR_RADIUS_M;

// SI input needed to prepare the ambient analytical state.  Unit conversion
// belongs outside this layer (swcme_units.hpp / swcme_core.hpp); this keeps the
// actual solar-wind physics independent of public input-unit conventions.
struct ConfigSI {
  double V_sw_m_s = 0.0;
  double n1AU_m3 = 0.0;
  double B1AU_T = 0.0;
  double T_K = 0.0;
  double gamma_ad = 5.0 / 3.0;

  // B1AU_T is interpreted as total |B| at this documented reference latitude.
  // The 3-D local Parker winding still uses the local sin(theta); this field
  // affects only the conversion of total B1AU into the common radial Br1AU.
  double reference_sin_theta = 1.0;

  // Explicit solar rotation rate shared by the Parker vector and analytical
  // connectivity mapping.  1-D normally uses the default model convention.
  double solar_rotation_rate_rad_s =
      swcme::constants::SOLAR_ROTATION_RAD_S;
};

// Per-step/per-model immutable ambient cache.  It contains all expensive or
// repeated normalization algebra needed by hot field evaluators.
struct PreparedState {
  double V_sw_m_s = 0.0;
  double T_K = 0.0;
  double gamma_ad = 5.0 / 3.0;
  double solar_rotation_rate_rad_s =
      swcme::constants::SOLAR_ROTATION_RAD_S;
  double reference_sin_theta = 1.0;

  // Parker cache. k_AU_equatorial = Omega*AU/Vsw; local latitude enters only
  // when the field is evaluated. Br1AU_T is the common radial normalization.
  double k_AU_equatorial = 0.0;
  double Br1AU_T = 0.0;

  // Leblanc SI coefficients with powers of solar radius absorbed so that the
  // hot-path evaluation is n = C2/r^2 + C4/r^4 + C6/r^6 [m^-3].
  double C2 = 0.0;
  double C4 = 0.0;
  double C6 = 0.0;
};

struct ParkerComponents {
  double Br_T = 0.0;
  double Bphi_T = 0.0;
  double Bmag_T = 0.0;
};

inline PreparedState prepare(const ConfigSI& cfg) {
  PreparedState state;
  state.V_sw_m_s = cfg.V_sw_m_s;
  state.T_K = cfg.T_K;
  state.gamma_ad = cfg.gamma_ad;
  state.solar_rotation_rate_rad_s = cfg.solar_rotation_rate_rad_s;
  state.reference_sin_theta = cfg.reference_sin_theta;

  // Configuration validation guarantees V_sw>0 before this routine is called.
  // Keeping the formula branch-free here makes any future invalid call fail
  // visibly rather than silently substituting a different physical state.
  state.k_AU_equatorial =
      cfg.solar_rotation_rate_rad_s * swcme::constants::AU_M / cfg.V_sw_m_s;

  const double reference_pitch =
      state.k_AU_equatorial * cfg.reference_sin_theta;
  state.Br1AU_T = cfg.B1AU_T /
      std::sqrt(1.0 + reference_pitch * reference_pitch);

  // Normalize the Leblanc profile exactly at 1 AU.  The calculation is kept in
  // one production location so 1-D and 3-D cannot drift in coefficients,
  // radius powers, or cm^-3 -> m^-3 scaling.
  const double rs_over_au =
      swcme::constants::SOLAR_RADIUS_M / swcme::constants::AU_M;
  const double rs2 = rs_over_au * rs_over_au;
  const double rs4 = rs2 * rs2;
  const double rs6 = rs4 * rs2;
  const double n1au_base_cm3 =
      LEBLANC_A_CM3 * rs2 +
      LEBLANC_B_CM3 * rs4 +
      LEBLANC_C_CM3 * rs6;
  const double n1au_base_m3 = n1au_base_cm3 * 1.0e6;
  const double scale = cfg.n1AU_m3 / n1au_base_m3;

  const double Rs = swcme::constants::SOLAR_RADIUS_M;
  const double Rs2 = Rs * Rs;
  const double Rs4 = Rs2 * Rs2;
  const double Rs6 = Rs4 * Rs2;
  state.C2 = scale * LEBLANC_A_CM3 * 1.0e6 * Rs2;
  state.C4 = scale * LEBLANC_B_CM3 * 1.0e6 * Rs4;
  state.C6 = scale * LEBLANC_C_CM3 * 1.0e6 * Rs6;
  return state;
}

inline double density_m3(const PreparedState& state, double radius_m) {
  const double r = radius_m;
  const double inv2 = 1.0 / (r * r);
  const double inv4 = inv2 * inv2;
  const double inv6 = inv4 * inv2;
  return state.C2 * inv2 + state.C4 * inv4 + state.C6 * inv6;
}

// Return Parker components at a specified local sin(colatitude).  This is the
// exact common scalar physics used by both dimensional models: 1-D supplies
// its fixed ray latitude, while 3-D derives sin(theta) from geometry.
inline ParkerComponents parker_components(const PreparedState& state,
                                           double radius_m,
                                           double sin_theta_local) {
  const double r = radius_m;
  const double r_AU = r / swcme::constants::AU_M;
  const double Br = state.Br1AU_T / (r_AU * r_AU);
  const double Bphi = -Br * state.k_AU_equatorial * r_AU * sin_theta_local;
  return {Br, Bphi, std::sqrt(Br * Br + Bphi * Bphi)};
}

// Construct the Cartesian Parker vector for a normalized radial direction and
// normalized solar-rotation axis.  The azimuthal basis is Omega_hat x e_r.  At
// the poles the cross product vanishes; because B_phi also vanishes there, the
// physical limiting field is purely radial and no arbitrary transverse basis
// is introduced.
inline std::array<double, 3> parker_field_cartesian(
    const PreparedState& state,
    const std::array<double, 3>& solar_axis_hat,
    const std::array<double, 3>& radial_hat,
    double radius_m) {
  const std::array<double, 3> cross = {{
      solar_axis_hat[1] * radial_hat[2] - solar_axis_hat[2] * radial_hat[1],
      solar_axis_hat[2] * radial_hat[0] - solar_axis_hat[0] * radial_hat[2],
      solar_axis_hat[0] * radial_hat[1] - solar_axis_hat[1] * radial_hat[0]}};
  const double sin_theta = std::sqrt(
      cross[0] * cross[0] + cross[1] * cross[1] + cross[2] * cross[2]);

  const ParkerComponents components =
      parker_components(state, radius_m, sin_theta);

  constexpr double AXIS_TOL =
      64.0 * std::numeric_limits<double>::epsilon();
  if (sin_theta <= AXIS_TOL) {
    return {{components.Br_T * radial_hat[0],
             components.Br_T * radial_hat[1],
             components.Br_T * radial_hat[2]}};
  }

  const double inv_sin_theta = 1.0 / sin_theta;
  const std::array<double, 3> ephi = {{
      cross[0] * inv_sin_theta,
      cross[1] * inv_sin_theta,
      cross[2] * inv_sin_theta}};
  return {{components.Br_T * radial_hat[0] + components.Bphi_T * ephi[0],
           components.Br_T * radial_hat[1] + components.Bphi_T * ephi[1],
           components.Br_T * radial_hat[2] + components.Bphi_T * ephi[2]}};
}

// Current SWCME shock physics uses a proton-only thermal closure.  Centralizing
// this small relation guarantees that 1-D and 3-D construct identical upstream
// pressure from the same density and temperature until a richer composition/
// electron-pressure model is intentionally introduced and validated.
inline double proton_pressure_Pa(const PreparedState& state,
                                 double number_density_m3) {
  return number_density_m3 * swcme::constants::BOLTZMANN_J_K * state.T_K;
}

}  // namespace solarwind
}  // namespace swcme
