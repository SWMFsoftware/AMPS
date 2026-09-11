#ifndef SWCME_ACCELERATION_HPP
#define SWCME_ACCELERATION_HPP

// ============================================================================
// swcme_acceleration.hpp
// ----------------------------------------------------------------------------
// Common shock-acceleration representation shared by the 1-D and 3-D SWCME
// interfaces.
//
// Motivation
// ----------
// A focused-transport solver can obtain first-order shock acceleration in two
// physically equivalent ways:
//   (1) inject a spectrum that already represents diffusive shock acceleration
//       (DSA), or
//   (2) resolve the compressive velocity gradient and let the transport
//       equation's -(1/3)(div V) p d f/dp term generate the acceleration.
// Enabling both for the same population double-counts the same shock physics.
// SWCME therefore exposes ONE enum, Mode, as the authoritative switch.  There
// are no independent booleans that can accidentally enable both paths.
//
// SOURCE
//   * the shock is an explicit injection/source surface;
//   * the DSA phase-space momentum slope q=3r/(r-1) is supplied when a physical
//     fast shock exists;
//   * the transport background must use SHOCK_ONLY, so no RH compression layer
//     is present in the velocity field and the source cannot be accelerated a
//     second time by a numerically resolved shock.
//
// RESOLVED_COMPRESSION
//   * no DSA source spectrum is enabled;
//   * the RH upstream/downstream velocity change is represented by one shared
//     finite-width C1 transition in swcme_regions.hpp;
//   * the transport solver may then use div(V) to produce compression
//     acceleration.  The exact discontinuous RH state remains available from
//     the shock diagnostic API and is not replaced by the numerical profile.
//
// The small ShockAccelerationState below is deliberately an acceleration
// contract, not the final AMPS SEP-source API.  A later adapter can add energy
// grids, absolute units/normalization, and serialization metadata without
// changing the mutually exclusive physics decision implemented here.
// ============================================================================

#include <array>
#include <cmath>
#include <iomanip>
#include <limits>
#include <sstream>
#include <string>

namespace swcme {
namespace acceleration {

enum class Mode {
  Source,
  ResolvedCompression
};

inline const char* mode_name(Mode mode) {
  return mode == Mode::Source ? "SOURCE" : "RESOLVED_COMPRESSION";
}

// A normalized source weight is intentionally dimensionless.  The baseline
// connectivity/perpendicular-diffusion experiment uses one constant weight per
// unit sampled shock area; conversion to a physical injection rate is owned by
// the later AMPS-facing source adapter, where the particle/spectral units can
// be defined without ambiguity.
struct Config {
  Mode mode = Mode::ResolvedCompression;
  double relative_source_weight_per_area = 1.0;
};

struct ShockAccelerationState {
  Mode mode = Mode::ResolvedCompression;
  bool surface_exists = false;
  bool physical_shock = false;

  // These two flags are exact complements for a physical shock.  SOURCE never
  // exposes a resolved-compression accelerator, and RESOLVED_COMPRESSION never
  // exposes a prescribed DSA source.
  bool source_enabled = false;
  bool resolved_compression_enabled = false;

  double time_s = 0.0;
  std::array<double,3> position_m{{0.0,0.0,0.0}};
  std::array<double,3> normal{{0.0,0.0,0.0}};
  double radius_m = 0.0;
  double normal_speed_m_s = 0.0;
  double compression = 1.0;
  double theta_Bn_rad = 0.0;
  double fast_mach = 0.0;
  double upstream_density_m3 = 0.0;
  double upstream_B_T = 0.0;

  // Valid only when source_enabled=true.  NaN in RESOLVED_COMPRESSION mode is
  // intentional: downstream code cannot silently consume a DSA slope when the
  // selected representation says acceleration must come from div(V).
  double dsa_q_phase_space = std::numeric_limits<double>::quiet_NaN();

  // Relative, dimensionless weighting used by the controlled baseline.  This
  // is zero whenever source_enabled=false.
  double relative_source_weight_per_area = 0.0;
};

inline double dsa_phase_space_slope(double compression) {
  if (!std::isfinite(compression) || !(compression > 1.0)) {
    return std::numeric_limits<double>::quiet_NaN();
  }
  return 3.0*compression/(compression-1.0);
}

// Construct the common acceleration record from already-validated shock
// diagnostics.  This function contains the ONLY production decision that turns
// a physical shock into either an explicit DSA source or a resolved compression
// accelerator.  The 1-D and 3-D wrappers merely provide geometry/state data.
inline ShockAccelerationState make_state(
    const Config& config,
    bool surface_exists,
    bool physical_shock,
    double time_s,
    const std::array<double,3>& position_m,
    const std::array<double,3>& normal,
    double normal_speed_m_s,
    double compression,
    double theta_Bn_rad,
    double fast_mach,
    double upstream_density_m3,
    double upstream_B_T) {
  ShockAccelerationState out;
  out.mode=config.mode;
  out.surface_exists=surface_exists;
  out.physical_shock=surface_exists && physical_shock;
  out.time_s=time_s;
  out.position_m=position_m;
  out.normal=normal;
  out.radius_m=std::sqrt(position_m[0]*position_m[0] +
                         position_m[1]*position_m[1] +
                         position_m[2]*position_m[2]);
  out.normal_speed_m_s=normal_speed_m_s;
  out.compression=out.physical_shock ? compression : 1.0;
  out.theta_Bn_rad=theta_Bn_rad;
  out.fast_mach=fast_mach;
  out.upstream_density_m3=upstream_density_m3;
  out.upstream_B_T=upstream_B_T;

  if (out.physical_shock && config.mode==Mode::Source) {
    out.source_enabled=true;
    out.resolved_compression_enabled=false;
    out.dsa_q_phase_space=dsa_phase_space_slope(out.compression);
    out.relative_source_weight_per_area=config.relative_source_weight_per_area;
  } else if (out.physical_shock && config.mode==Mode::ResolvedCompression) {
    out.source_enabled=false;
    out.resolved_compression_enabled=true;
    out.dsa_q_phase_space=std::numeric_limits<double>::quiet_NaN();
    out.relative_source_weight_per_area=0.0;
  }
  return out;
}

// Deterministic scientific-notation serialization is used only for regression
// and audit comparisons.  It deliberately contains every acceleration input
// needed to prove that equivalent 1-D/3-D fixtures feed the same transport
// experiment before dimensional transport effects are introduced.
inline std::string serialize_csv(const ShockAccelerationState& s) {
  std::ostringstream out;
  out.setf(std::ios::scientific);
  out << std::setprecision(17)
      << mode_name(s.mode) << ','
      << (s.surface_exists ? 1 : 0) << ','
      << (s.physical_shock ? 1 : 0) << ','
      << (s.source_enabled ? 1 : 0) << ','
      << (s.resolved_compression_enabled ? 1 : 0) << ','
      << s.time_s << ','
      << s.position_m[0] << ',' << s.position_m[1] << ',' << s.position_m[2] << ','
      << s.normal[0] << ',' << s.normal[1] << ',' << s.normal[2] << ','
      << s.radius_m << ',' << s.normal_speed_m_s << ',' << s.compression << ','
      << s.theta_Bn_rad << ',' << s.fast_mach << ','
      << s.upstream_density_m3 << ',' << s.upstream_B_T << ',';
  if (std::isfinite(s.dsa_q_phase_space)) out << s.dsa_q_phase_space;
  else out << "NA";
  out << ',' << s.relative_source_weight_per_area;
  return out.str();
}

}  // namespace acceleration
}  // namespace swcme

#endif
