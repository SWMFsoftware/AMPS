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

// Nominal region labels.  The two transition labels represent finite C1 blend
// zones around the sheath/ejecta leading edge and ejecta/post-ICME trailing
// edge.  The physical shock itself is deliberately NOT smoothed here; r>=R_sh
// is upstream and r<R_sh approaches the exact RH downstream state.
enum class Region {
  Upstream,
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
  double leading_smooth_fraction = 0.02;
  double trailing_smooth_fraction = 0.03;
  double sheath_ramp_power = 2.0;
  double V_sheath_LE_factor = 1.10;
  double f_ME = 0.50;
  double V_ME_factor = 0.80;
};

// Region boundaries along ONE physical ray.  R_sh is the local surface radius,
// not necessarily the apex radius.  R_le and R_te are scaled from that local
// radius, so their ratios to R_sh are independent of direction.
//
// smooth_le_width / smooth_te_width are TOTAL transition widths centered on the
// nominal boundary.  They are capped at 90% of the adjacent finite layer width
// so the leading- and trailing-edge blend zones cannot overlap or invert the
// nominal region ordering even for aggressive user smoothing parameters.
struct Boundaries {
  double R_sh_m = 0.0;
  double R_le_m = 0.0;
  double R_te_m = 0.0;
  double sheath_thickness_m = 0.0;
  double ejecta_thickness_m = 0.0;
  double smooth_le_width_m = 0.0;
  double smooth_te_width_m = 0.0;
};

struct Location {
  Region region = Region::Upstream;
  // Blend weight has a region-dependent interpretation:
  //  LeadingTransition: 0 -> pure sheath, 1 -> pure ejecta.
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

  // The requested smoothing parameters are specified as AU at a 1-AU shock,
  // so numerically they are the corresponding fractions of local R_sh.  Use a
  // symmetric transition about each nominal interface.  Limiting the TOTAL
  // width to 90% of the adjacent finite layer leaves at least 55% of each layer
  // on either side of its center and prevents leading/trailing blends from
  // crossing each other.
  const double requested_le = std::max(0.0, config.leading_smooth_fraction) *
                              local_shock_radius_m;
  const double le_limit = 0.90 * std::max(0.0,
      std::min(out.sheath_thickness_m, out.ejecta_thickness_m));
  out.smooth_le_width_m = std::min(requested_le, le_limit);

  const double requested_te = std::max(0.0, config.trailing_smooth_fraction) *
                              local_shock_radius_m;
  const double te_limit = 0.90 * std::max(0.0, out.ejecta_thickness_m);
  out.smooth_te_width_m = std::min(requested_te, te_limit);
  return out;
}

inline Location locate(double r_m, const Boundaries& b) {
  Location out;
  if (r_m >= b.R_sh_m) {
    out.region = Region::Upstream;
    return out;
  }

  const double h_le = 0.5 * b.smooth_le_width_m;
  const double h_te = 0.5 * b.smooth_te_width_m;

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
  const double width = b.R_sh_m - b.R_le_m;
  if (!(width > 0.0)) return 1.0;
  return clamp01((b.R_sh_m - r_m) / width);
}

inline double sheath_profile_weight(double r_m, const Boundaries& b,
                                    double ramp_power) {
  const double s = sheath_progress(r_m, b);
  return smoothstep01(std::pow(s, std::max(1.0, ramp_power)));
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
