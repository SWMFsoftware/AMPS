#ifndef SWCME_REGIONS_HPP
#define SWCME_REGIONS_HPP

// ============================================================================
// swcme_regions.hpp
// ----------------------------------------------------------------------------
// Shared sheath / magnetic-ejecta region geometry and blending utilities.
//
// Why this file exists
// --------------------
// The original 1-D and 3-D implementations carried independent region logic.
// That duplication produced several scientifically important inconsistencies:
//   * the 1-D ejecta path forced f_ME and V_ME_factor to be >= 1, so density
//     depletion and slower ejecta could not actually be represented;
//   * the 3-D code subtracted apex-scaled *absolute* layer widths from every
//     local shock radius, so narrow flanks were not self-similar and could
//     develop distorted or inverted layers;
//   * leading/trailing-edge smoothing was handled differently in 1-D and 3-D;
//   * there was no explicit SHOCK_ONLY mode for the controlled SEP experiment.
//
// This header provides only dimension-independent region mathematics.  It does
// not know about Cartesian vs radial coordinates and it does not compute MHD
// shock jumps.  The dimensional wrappers provide the local shock state and the
// local Parker/upstream state, then use the common boundaries and blend weights
// below.  Keeping the module header-only preserves the lightweight character of
// SWCME while making the region contract identical in 1-D and 3-D.
// ============================================================================

#include <algorithm>
#include <cmath>

namespace swcme {
namespace regions {

// SHOCK_ONLY is the controlled transport background: the analytical upstream
// Parker/Leblanc solar wind is returned everywhere and the shock exists only as
// explicit geometry/source bookkeeping.  FULL_ICME adds the phenomenological
// sheath and ejecta profiles behind the finite CME/shock surface.
enum class Mode {
  ShockOnly,
  FullICME
};

// Nominal region labels.  Leading/TrailingTransition represent the artificial
// sheath/ejecta and ejecta/post-ICME blends.  ShockTransition is different: it
// exists only for RESOLVED_COMPRESSION and is the numerical representation of
// the physical RH discontinuity.  SOURCE mode sets its width to zero so the
// transport background contains no resolved shock accelerator.
enum class Region {
  Upstream,
  ShockTransition,
  Sheath,
  LeadingTransition,
  Ejecta,
  TrailingTransition,
  PostICME
};

inline const char* mode_name(Mode mode) {
  return mode == Mode::ShockOnly ? "SHOCK_ONLY" : "FULL_ICME";
}

inline const char* region_name(Region region) {
  switch (region) {
    case Region::Upstream: return "UPSTREAM";
    case Region::ShockTransition: return "SHOCK_TRANSITION";
    case Region::Sheath: return "SHEATH";
    case Region::LeadingTransition: return "LEADING_TRANSITION";
    case Region::Ejecta: return "EJECTA";
    case Region::TrailingTransition: return "TRAILING_TRANSITION";
    case Region::PostICME: return "POST_ICME";
  }
  return "UNKNOWN";
}

// Public region parameters converted into dimensionless self-similar fractions.
// Quantities named *_at1AU are numerically both AU and fractions of a 1-AU
// reference radius.  For example, a 0.10-AU sheath at R_sh=1 AU corresponds to
// sheath_fraction=0.10 and therefore to a local sheath thickness 0.10*R_sh at
// EVERY apex/flank direction.  This is the key correction that preserves
// self-similarity on a finite non-spherical shock surface.
struct Config {
  Mode mode = Mode::FullICME;
  double sheath_fraction = 0.10;
  double ejecta_fraction = 0.20;
  // Numerical shock width used ONLY by RESOLVED_COMPRESSION.  SOURCE mode
  // passes zero here, leaving the physical source surface out of the flow.
  double shock_smooth_fraction = 0.01;
  double leading_smooth_fraction = 0.02;
  double trailing_smooth_fraction = 0.03;
  double sheath_ramp_power = 2.0;
  double V_sheath_LE_factor = 1.10;
  double f_ME = 0.50;
  double V_ME_factor = 0.80;
};

// CFG03 makes the historical 90-percent geometry margin an explicit input
// contract instead of a runtime correction.  Keeping the factor and all three
// derived limits in the region module gives configuration validation and
// boundary construction one authoritative definition of the policy.
constexpr double SMOOTHING_LAYER_FRACTION_LIMIT = 0.90;

struct SmoothingFractionLimits {
  double shock = 0.0;
  double leading = 0.0;
  double trailing = 0.0;
};

inline SmoothingFractionLimits smoothing_fraction_limits(
    double sheath_fraction,double ejecta_fraction) noexcept {
  // Invalid/non-finite thicknesses are diagnosed independently by the common
  // configuration validator.  Returning zero here avoids manufacturing a
  // meaningful limit from malformed geometry and, importantly, performs no
  // correction of a user smoothing width.
  const double sheath=(std::isfinite(sheath_fraction) &&
                       sheath_fraction>0.0) ? sheath_fraction : 0.0;
  const double ejecta=(std::isfinite(ejecta_fraction) &&
                       ejecta_fraction>0.0) ? ejecta_fraction : 0.0;
  SmoothingFractionLimits limits;
  limits.shock=SMOOTHING_LAYER_FRACTION_LIMIT*sheath;
  limits.leading=SMOOTHING_LAYER_FRACTION_LIMIT*
                 std::min(sheath,ejecta);
  limits.trailing=SMOOTHING_LAYER_FRACTION_LIMIT*ejecta;
  return limits;
}

// Region boundaries along ONE physical ray.  R_sh is the local surface radius,
// not necessarily the apex radius.  R_le and R_te are scaled from that local
// radius, so their ratios to R_sh are independent of direction.
//
// smooth_le_width / smooth_te_width are TOTAL transition widths centered on the
// nominal boundary.  Central configuration validation guarantees that each is
// no more than 90% of its adjacent finite layer.  Boundary construction then
// preserves the accepted request exactly; it never clips a scientific input.
struct Boundaries {
  double R_sh_m = 0.0;
  double R_le_m = 0.0;
  double R_te_m = 0.0;
  double sheath_thickness_m = 0.0;
  double ejecta_thickness_m = 0.0;
  // Total C1 width centered on the mathematical shock surface.  Validation
  // ensures its inner edge remains inside the sheath and cannot overlap LE.
  double smooth_shock_width_m = 0.0;
  double smooth_le_width_m = 0.0;
  double smooth_te_width_m = 0.0;
};

struct Location {
  Region region = Region::Upstream;
  // Blend weight has a region-dependent interpretation:
  //  ShockTransition:    0 -> upstream, 1 -> exact RH downstream.
  //  LeadingTransition:  0 -> pure sheath, 1 -> pure ejecta.
  //  TrailingTransition: 0 -> pure ejecta, 1 -> pure post-ICME ambient.
  //  Other regions:      0.
  double blend = 0.0;
};

inline double clamp01(double x) {
  return x < 0.0 ? 0.0 : (x > 1.0 ? 1.0 : x);
}

inline double smoothstep01(double x) {
  x = clamp01(x);
  return x*x*(3.0 - 2.0*x);
}

inline double lerp(double a, double b, double w) {
  return a + (b-a)*w;
}

inline double log_lerp_positive(double a, double b, double w) {
  // Density is positive by model construction.  A tiny floor protects the
  // logarithm from underflow without changing any physically valid fixture.
  const double aa = std::max(a, 1.0e-300);
  const double bb = std::max(b, 1.0e-300);
  return std::exp(lerp(std::log(aa), std::log(bb), clamp01(w)));
}

inline Boundaries make_boundaries(double local_shock_radius_m,
                                  const Config& config) {
  Boundaries out;
  out.R_sh_m = local_shock_radius_m;
  if (!(local_shock_radius_m > 0.0) || !std::isfinite(local_shock_radius_m)) {
    return out;
  }

  const double fs = config.sheath_fraction;
  const double fe = config.ejecta_fraction;
  out.R_le_m = local_shock_radius_m * (1.0 - fs);
  out.R_te_m = local_shock_radius_m * (1.0 - fs - fe);
  out.sheath_thickness_m = out.R_sh_m - out.R_le_m;
  out.ejecta_thickness_m = out.R_le_m - out.R_te_m;

  // Configuration validation has already proved each requested fraction lies
  // inside the explicit 90-percent limit.  Applying the values directly is the
  // essential CFG03 behavior: a successful run's effective widths are exactly
  // those recorded in its configuration, at every local 3-D surface radius.
  // SOURCE supplies shock_smooth_fraction=0, so its injection surface still
  // contributes no resolved compression layer.
  out.smooth_shock_width_m=config.shock_smooth_fraction*
                           local_shock_radius_m;
  out.smooth_le_width_m=config.leading_smooth_fraction*
                        local_shock_radius_m;
  out.smooth_te_width_m=config.trailing_smooth_fraction*
                        local_shock_radius_m;
  return out;
}

inline Location locate(double r_m, const Boundaries& b) {
  Location out;
  const double h_sh = 0.5 * b.smooth_shock_width_m;
  const double h_le = 0.5 * b.smooth_le_width_m;
  const double h_te = 0.5 * b.smooth_te_width_m;

  // RESOLVED_COMPRESSION represents the physical discontinuity by one finite
  // C1 layer centered on R_sh.  blend=0 is the outer/upstream endpoint and
  // blend=1 is the inner/exact-RH endpoint.  With zero width (SOURCE or legacy
  // unsmoothed diagnostics) the old mathematical discontinuity is recovered.
  if (b.smooth_shock_width_m > 0.0) {
    if (r_m >= b.R_sh_m + h_sh) {
      out.region = Region::Upstream;
      return out;
    }
    if (r_m >= b.R_sh_m - h_sh) {
      const double q=(b.R_sh_m + h_sh - r_m)/b.smooth_shock_width_m;
      out.region=Region::ShockTransition;
      out.blend=smoothstep01(q);
      return out;
    }
  } else if (r_m >= b.R_sh_m) {
    out.region = Region::Upstream;
    return out;
  }

  if (b.smooth_le_width_m > 0.0 &&
      r_m <= b.R_le_m + h_le && r_m >= b.R_le_m - h_le) {
    // r decreases inward.  At the outer edge the state is pure sheath; at the
    // inner edge it is pure ejecta.  smoothstep gives zero derivative at both
    // endpoints, so the artificial interface is C1 when the base states are
    // finite and smooth.
    const double q = (b.R_le_m + h_le - r_m) / b.smooth_le_width_m;
    out.region = Region::LeadingTransition;
    out.blend = smoothstep01(q);
    return out;
  }

  if (r_m > b.R_le_m + h_le) {
    out.region = Region::Sheath;
    return out;
  }

  if (b.smooth_te_width_m > 0.0 &&
      r_m <= b.R_te_m + h_te && r_m >= b.R_te_m - h_te) {
    // At the outer edge return pure ejecta; at the inner edge return pure
    // post-ICME ambient.  As at the leading edge, the smoothstep is C1.
    const double q = (b.R_te_m + h_te - r_m) / b.smooth_te_width_m;
    out.region = Region::TrailingTransition;
    out.blend = smoothstep01(q);
    return out;
  }

  if (r_m > b.R_te_m + h_te) {
    out.region = Region::Ejecta;
    return out;
  }

  out.region = Region::PostICME;
  return out;
}

// Progress through the sheath profile: 0 immediately downstream of the shock,
// 1 at the nominal leading edge.  The function is deliberately clamped so it
// can be evaluated inside the symmetric leading-edge transition without a
// second extrapolation convention.
inline double sheath_progress(double r_m, const Boundaries& b) {
  // In resolved-compression mode the exact RH state is reached at the INNER
  // edge of the numerical shock layer.  Starting the phenomenological sheath
  // relaxation there gives C0/C1 matching: the shock smoothstep has zero slope
  // at its inner endpoint and the sheath smoothstep has zero slope at s=0.
  const double sheath_start=b.R_sh_m-0.5*b.smooth_shock_width_m;
  const double width = sheath_start - b.R_le_m;
  if (!(width > 0.0)) return 1.0;
  return clamp01((sheath_start - r_m) / width);
}

inline double sheath_profile_weight(double r_m, const Boundaries& b,
                                    double ramp_power) {
  const double s = sheath_progress(r_m, b);
  return smoothstep01(std::pow(s, std::max(1.0, ramp_power)));
}

// Value and radial derivative of one scalar blend.  Derivatives are expressed
// with respect to increasing heliocentric radius r [m], so a transition that
// progresses inward naturally has a negative d(weight)/dr.  These helpers are
// used by the analytical 1-D velocity divergence and intentionally reproduce
// exactly the same smoothstep weights used by locate()/field evaluation.
struct WeightDerivative {
  double value = 0.0;
  double d_dr = 0.0; // derivative per meter
};

inline WeightDerivative inward_transition_weight(double r_m,
                                                  double center_m,
                                                  double total_width_m) {
  WeightDerivative out;
  if (!(total_width_m > 0.0)) {
    out.value = r_m < center_m ? 1.0 : 0.0;
    return out;
  }
  const double q=(center_m+0.5*total_width_m-r_m)/total_width_m;
  if (q<=0.0) return out;
  if (q>=1.0) { out.value=1.0; return out; }
  out.value=smoothstep01(q);
  // d/dq [q^2(3-2q)] = 6 q (1-q), dq/dr = -1/width.
  out.d_dr=-6.0*q*(1.0-q)/total_width_m;
  return out;
}

inline WeightDerivative sheath_profile_weight_with_derivative(
    double r_m, const Boundaries& b, double ramp_power) {
  WeightDerivative out;
  const double sheath_start=b.R_sh_m-0.5*b.smooth_shock_width_m;
  const double width=sheath_start-b.R_le_m;
  if (!(width>0.0)) { out.value=1.0; return out; }

  const double raw_s=(sheath_start-r_m)/width;
  if (raw_s<=0.0) return out;
  if (raw_s>=1.0) { out.value=1.0; return out; }

  const double p=std::max(1.0,ramp_power);
  const double sp=std::pow(raw_s,p);
  out.value=smoothstep01(sp);
  const double dsp_ds = (p==1.0) ? 1.0 : p*std::pow(raw_s,p-1.0);
  const double ds_dr=-1.0/width;
  out.d_dr=6.0*sp*(1.0-sp)*dsp_ds*ds_dr;
  return out;
}

// Exact radial velocity profile and derivative used by the 1-D transport
// background.  This function owns the velocity-side interpretation of the
// common region boundaries, ensuring that V(r) and dV/dr cannot drift apart.
// Density and magnetic-field profiles remain separate because their physical
// closures differ, but the region label is returned for audit/debug output.
struct RadialVelocityState {
  Region region = Region::Upstream;
  double velocity_m_s = 0.0;
  double d_velocity_dr_s_inv = 0.0; // (m/s)/m = 1/s
};

inline RadialVelocityState radial_velocity_state(
    double r_m, const Boundaries& b, bool has_shock,
    double upstream_velocity_m_s, double downstream_velocity_m_s,
    double leading_edge_velocity_m_s, double ejecta_velocity_m_s,
    double sheath_ramp_power) {
  RadialVelocityState out;
  const Location loc=locate(r_m,b);
  out.region=loc.region;
  out.velocity_m_s=upstream_velocity_m_s;

  if (loc.region==Region::Upstream || loc.region==Region::PostICME) return out;

  if (loc.region==Region::ShockTransition) {
    if (!has_shock) return out;
    const WeightDerivative w=inward_transition_weight(
        r_m,b.R_sh_m,b.smooth_shock_width_m);
    out.velocity_m_s=lerp(upstream_velocity_m_s,downstream_velocity_m_s,w.value);
    out.d_velocity_dr_s_inv=(downstream_velocity_m_s-upstream_velocity_m_s)*w.d_dr;
    return out;
  }

  auto sheath=[&](double rr) {
    RadialVelocityState sh;
    sh.region=Region::Sheath;
    if (!has_shock) {
      sh.velocity_m_s=upstream_velocity_m_s;
      return sh;
    }
    const WeightDerivative w=sheath_profile_weight_with_derivative(
        rr,b,sheath_ramp_power);
    sh.velocity_m_s=lerp(downstream_velocity_m_s,leading_edge_velocity_m_s,w.value);
    sh.d_velocity_dr_s_inv=(leading_edge_velocity_m_s-downstream_velocity_m_s)*w.d_dr;
    return sh;
  };

  if (loc.region==Region::Sheath) return sheath(r_m);

  if (loc.region==Region::LeadingTransition) {
    const RadialVelocityState sh=sheath(r_m);
    const WeightDerivative w=inward_transition_weight(
        r_m,b.R_le_m,b.smooth_le_width_m);
    out.velocity_m_s=lerp(sh.velocity_m_s,ejecta_velocity_m_s,w.value);
    out.d_velocity_dr_s_inv=(1.0-w.value)*sh.d_velocity_dr_s_inv
        +(ejecta_velocity_m_s-sh.velocity_m_s)*w.d_dr;
    return out;
  }

  if (loc.region==Region::Ejecta) {
    out.velocity_m_s=ejecta_velocity_m_s;
    return out;
  }

  if (loc.region==Region::TrailingTransition) {
    const WeightDerivative w=inward_transition_weight(
        r_m,b.R_te_m,b.smooth_te_width_m);
    out.velocity_m_s=lerp(ejecta_velocity_m_s,upstream_velocity_m_s,w.value);
    out.d_velocity_dr_s_inv=(upstream_velocity_m_s-ejecta_velocity_m_s)*w.d_dr;
    return out;
  }

  return out;
}

// For a forward shock, the phenomenological sheath target at the leading edge
// should remain between ambient V_sw and the outward radial component of the
// exact RH downstream velocity.  The configurable factor therefore shapes the
// relaxation but is not allowed to violate either physical boundary.
inline double leading_edge_speed(double V_sw_m_s,
                                 double downstream_radial_m_s,
                                 double configured_factor) {
  const double low = V_sw_m_s;
  const double high = std::max(low, downstream_radial_m_s);
  const double requested = configured_factor * V_sw_m_s;
  return std::max(low, std::min(requested, high));
}

}  // namespace regions
}  // namespace swcme

#endif
