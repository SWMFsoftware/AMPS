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
//   * Parker field-line distance and magnetic focusing length, and
//   * configurable proton-only or electron/proton/alpha thermodynamics.
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
// Keep the public-input threshold in solar-radius units beside its SI form.
// CFG04 uses MIN_RADIUS_RS while validating r0_Rs and data-driven knots before
// unit conversion; evaluators use MIN_RADIUS_M.  Deriving the latter here
// prevents the two boundaries from drifting apart during future revisions.
constexpr double MIN_RADIUS_RS = 1.05;
constexpr double MIN_RADIUS_M =
    MIN_RADIUS_RS * swcme::constants::SOLAR_RADIUS_M;

// The legacy closure remains the numeric zero/default so aggregate/default
// construction preserves every historical result.  MultiSpecies interprets
// the Leblanc value as electron density and obtains charge-neutral proton and
// alpha populations from f_alpha=n_alpha/n_proton.
enum class ThermodynamicClosure {
  ProtonOnly = 0,
  MultiSpecies = 1
};

inline const char* thermodynamic_closure_name(ThermodynamicClosure closure) {
  switch (closure) {
    case ThermodynamicClosure::ProtonOnly: return "PROTON_ONLY";
    case ThermodynamicClosure::MultiSpecies: return "MULTI_SPECIES";
  }
  return "UNKNOWN";
}

// SI input needed to prepare the ambient analytical state.  Unit conversion
// belongs outside this layer (swcme_units.hpp / swcme_core.hpp); this keeps the
// actual solar-wind physics independent of public input-unit conventions.
struct ConfigSI {
  double V_sw_m_s = 0.0;
  double n1AU_m3 = 0.0;
  double B1AU_T = 0.0;
  double T_K = 0.0;
  double gamma_ad = 5.0 / 3.0;
  ThermodynamicClosure thermodynamic_closure =
      ThermodynamicClosure::ProtonOnly;
  double alpha_to_proton_ratio = 0.0;
  double electron_T_K = 0.0;
  double alpha_T_K = 0.0;
  int parker_radial_polarity = +1;

  // B1AU_T is interpreted as total |B| at this documented reference latitude.
  // The 3-D local Parker winding still uses the local sin(theta); this field
  // affects only the conversion of total B1AU into the common radial Br1AU.
  double reference_sin_theta = 1.0;

  // Explicit solar rotation rate shared by the Parker vector and analytical
  // connectivity mapping.  1-D normally uses the default model convention.
  double solar_rotation_rate_rad_s =
      swcme::constants::SOLAR_ROTATION_RAD_S;
  double parker_source_radius_m = 0.0;
};

// Per-step/per-model immutable ambient cache.  It contains all expensive or
// repeated normalization algebra needed by hot field evaluators.
struct PreparedState {
  double V_sw_m_s = 0.0;
  double T_K = 0.0;
  double gamma_ad = 5.0 / 3.0;
  ThermodynamicClosure thermodynamic_closure =
      ThermodynamicClosure::ProtonOnly;
  double alpha_to_proton_ratio = 0.0;
  double electron_T_K = 0.0;
  double alpha_T_K = 0.0;
  int parker_radial_polarity = +1;
  double solar_rotation_rate_rad_s =
      swcme::constants::SOLAR_ROTATION_RAD_S;
  double reference_sin_theta = 1.0;
  double parker_source_radius_m = 0.0;

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

// Return sin(theta), where theta is the colatitude relative to the solar
// rotation axis.  Both inputs are unit vectors at the public call sites.  The
// clamp protects the exact [0,1] trigonometric range from a few ulps of norm
// drift without replacing a genuinely invalid vector; model validation owns
// rejection of zero/non-finite axes and directions before this hot helper.
inline double parker_sin_colatitude(
    const std::array<double, 3>& solar_axis_hat,
    const std::array<double, 3>& radial_hat) {
  const std::array<double, 3> cross = {{
      solar_axis_hat[1] * radial_hat[2] - solar_axis_hat[2] * radial_hat[1],
      solar_axis_hat[2] * radial_hat[0] - solar_axis_hat[0] * radial_hat[2],
      solar_axis_hat[0] * radial_hat[1] - solar_axis_hat[1] * radial_hat[0]}};
  const double magnitude=std::sqrt(cross[0] * cross[0] +
                                   cross[1] * cross[1] +
                                   cross[2] * cross[2]);
  return std::isfinite(magnitude) ? std::min(1.0,magnitude)
                                  : std::numeric_limits<double>::quiet_NaN();
}

inline PreparedState prepare(const ConfigSI& cfg) {
  PreparedState state;
  state.V_sw_m_s = cfg.V_sw_m_s;
  state.T_K = cfg.T_K;
  state.gamma_ad = cfg.gamma_ad;
  state.thermodynamic_closure = cfg.thermodynamic_closure;
  state.alpha_to_proton_ratio = cfg.alpha_to_proton_ratio;
  state.electron_T_K = cfg.electron_T_K;
  state.alpha_T_K = cfg.alpha_T_K;
  state.parker_radial_polarity = cfg.parker_radial_polarity;
  state.solar_rotation_rate_rad_s = cfg.solar_rotation_rate_rad_s;
  state.reference_sin_theta = cfg.reference_sin_theta;
  state.parker_source_radius_m = cfg.parker_source_radius_m;

  // Configuration validation guarantees V_sw>0 before this routine is called.
  // Keeping the formula branch-free here makes any future invalid call fail
  // visibly rather than silently substituting a different physical state.
  state.k_AU_equatorial =
      cfg.solar_rotation_rate_rad_s * swcme::constants::AU_M / cfg.V_sw_m_s;

  const double reference_pitch = state.k_AU_equatorial *
      (1.0 - cfg.parker_source_radius_m / swcme::constants::AU_M) *
      cfg.reference_sin_theta;
  state.Br1AU_T = static_cast<double>(cfg.parker_radial_polarity) * cfg.B1AU_T /
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
  const double winding_radius_AU =
      (r - state.parker_source_radius_m) / swcme::constants::AU_M;
  const double Bphi =
      -Br * state.k_AU_equatorial * winding_radius_AU * sin_theta_local;
  return {Br, Bphi, std::sqrt(Br * Br + Bphi * Bphi)};
}

// Closed-form distance along the same baseline Parker field used by the
// magnetic evaluator.  With k=Omega*sin(theta)/Vsw, ds/dr=sqrt(1+(kr)^2).
// Keeping its antiderivative here makes the model connectivity and AMPS
// adapter consume one production formula.  The explicit small-k branch avoids
// a removable 0/0 singularity and returns the exact radial-field limit.
inline double parker_path_length_m(const PreparedState& state,
                                   double sin_theta_local,
                                   double radius_a_m,
                                   double radius_b_m) {
  if (!std::isfinite(sin_theta_local) || sin_theta_local<0.0 ||
      sin_theta_local>1.0 || !std::isfinite(radius_a_m) ||
      !std::isfinite(radius_b_m) || !(radius_a_m>0.0) ||
      !(radius_b_m>0.0) || !std::isfinite(state.k_AU_equatorial)) {
    return std::numeric_limits<double>::quiet_NaN();
  }

  const double k=state.k_AU_equatorial*sin_theta_local/
                 swcme::constants::AU_M;
  if (std::abs(k)<=64.0*std::numeric_limits<double>::epsilon()/
                         std::max(radius_a_m,radius_b_m)) {
    return std::abs(radius_b_m-radius_a_m);
  }

  const auto primitive=[k,&state](double radius_m) {
    const double x=radius_m-state.parker_source_radius_m;
    const double kx=k*x;
    return 0.5*(x*std::sqrt(1.0+kx*kx)+std::asinh(kx)/k);
  };
  return std::abs(primitive(radius_b_m)-primitive(radius_a_m));
}

// Magnetic focusing length for an OUTWARD path coordinate s,
//
//   L_B = -1/(d ln|B|/ds)
//       = r [1+(kr)^2]^(3/2) / [2+(kr)^2].
//
// The baseline Parker magnitude decreases outward, hence L_B is positive and
// approaches r/2 in the radial (k -> 0) limit.  This function is intentionally
// colocated with parker_components(): an AMPS adapter must never differentiate
// or reinterpret a second magnetic-field formula to obtain focusing.
inline double parker_focusing_length_m(const PreparedState& state,
                                       double radius_m,
                                       double sin_theta_local) {
  if (!std::isfinite(radius_m) || !(radius_m>0.0) ||
      !std::isfinite(sin_theta_local) || sin_theta_local<0.0 ||
      sin_theta_local>1.0 || !std::isfinite(state.k_AU_equatorial)) {
    return std::numeric_limits<double>::quiet_NaN();
  }
  const double k=state.k_AU_equatorial*sin_theta_local/
                 swcme::constants::AU_M;
  const double x=radius_m-state.parker_source_radius_m;
  const double kx=k*x;
  const double one_plus_kx2=1.0+kx*kx;
  // This follows directly from -1/(d ln|B|/ds), with
  // |B| proportional to r^-2 sqrt(1+k^2(r-rb)^2).
  const double denominator=2.0*one_plus_kx2-k*k*radius_m*x;
  return radius_m*one_plus_kx2*std::sqrt(one_plus_kx2)/denominator;
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
  const double sin_theta=parker_sin_colatitude(solar_axis_hat,radial_hat);

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

struct ThermodynamicState {
  double electron_density_m3 = 0.0;
  double proton_density_m3 = 0.0;
  double alpha_density_m3 = 0.0;
  double mass_density_kg_m3 = 0.0;
  double pressure_Pa = 0.0;
  double sound_speed_m_s = 0.0;
};

// Close the Leblanc electron density into the primitive quantities consumed
// by the MHD shock.  Charge neutrality gives ne=np+2*na and the configured
// abundance gives na=f*np.  The proton-only branch is written separately—not
// as the f=0 multi-species limit—because its historical pressure intentionally
// excludes electron pressure and DEN05 requires exact default compatibility.
inline ThermodynamicState thermodynamic_state(
    const PreparedState& state, double electron_density_m3) {
  ThermodynamicState out;
  out.electron_density_m3 = electron_density_m3;
  if (state.thermodynamic_closure == ThermodynamicClosure::ProtonOnly) {
    out.proton_density_m3 = electron_density_m3;
    out.mass_density_kg_m3 =
        electron_density_m3 * swcme::constants::PROTON_MASS_KG;
    out.pressure_Pa = electron_density_m3 *
        swcme::constants::BOLTZMANN_J_K * state.T_K;
  } else {
    const double denominator = 1.0 + 2.0 * state.alpha_to_proton_ratio;
    out.proton_density_m3 = electron_density_m3 / denominator;
    out.alpha_density_m3 =
        state.alpha_to_proton_ratio * out.proton_density_m3;
    out.mass_density_kg_m3 =
        swcme::constants::PROTON_MASS_KG * out.proton_density_m3 +
        swcme::constants::ALPHA_PARTICLE_MASS_KG * out.alpha_density_m3;
    out.pressure_Pa = swcme::constants::BOLTZMANN_J_K *
        (out.proton_density_m3 * state.T_K +
         out.electron_density_m3 * state.electron_T_K +
         out.alpha_density_m3 * state.alpha_T_K);
  }
  out.sound_speed_m_s =
      std::sqrt(state.gamma_ad * out.pressure_Pa / out.mass_density_kg_m3);
  return out;
}

// Compatibility name retained for existing integrations.  Its result is now
// the pressure selected by the prepared closure, so adapter and shock callers
// cannot accidentally continue using proton-only algebra in MultiSpecies mode.
inline double proton_pressure_Pa(const PreparedState& state,
                                 double electron_density_m3) {
  return thermodynamic_state(state, electron_density_m3).pressure_Pa;
}

}  // namespace solarwind
}  // namespace swcme
