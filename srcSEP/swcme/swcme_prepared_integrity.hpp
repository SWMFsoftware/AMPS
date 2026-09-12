#pragma once

// ============================================================================
// swcme_prepared_integrity.hpp
// ----------------------------------------------------------------------------
// Explicit serializers for the shared fixed-size records embedded in 1-D and
// 3-D StepState.  PST06 must detect changes to physics fields without hashing
// C++ object memory: raw-memory hashes include padding, depend on ABI layout,
// and can change across a valid member-wise copy.  These routines instead feed
// each declared value to the stable ConfigurationDigestBuilder.
// ============================================================================

#include "swcme_acceleration.hpp"
#include "swcme_core.hpp"
#include "swcme_regions.hpp"
#include "swcme_shock.hpp"
#include "swcme_status.hpp"

namespace swcme {
namespace prepared_integrity {

inline void add_primitive(ConfigurationDigestBuilder& digest,
                          const shock::PrimitiveState& state) noexcept {
  digest.add_double(state.rho_kg_m3);
  digest.add_double(state.pressure_Pa);
  for (double value : state.velocity_m_s) digest.add_double(value);
  for (double value : state.magnetic_T) digest.add_double(value);
}

inline void add_jump(ConfigurationDigestBuilder& digest,
                     const shock::JumpResult& jump) noexcept {
  digest.add_bool(jump.has_shock);
  digest.add_bool(jump.solver_converged);
  digest.add_uint64(static_cast<std::uint64_t>(jump.status));
  digest.add_double(jump.compression);
  digest.add_double(jump.theta_Bn_rad);
  digest.add_double(jump.fast_speed_m_s);
  digest.add_double(jump.fast_mach);
  digest.add_double(jump.shock_normal_speed_m_s);
  digest.add_double(jump.upstream_inflow_normal_m_s);
  digest.add_uint64(static_cast<std::uint64_t>(
      static_cast<std::int64_t>(jump.root_iterations)));
  // Root-bracket diagnostics are part of the public JumpResult record and can
  // affect failure interpretation even though they do not alter the primitive
  // state.  Seal them explicitly so post-preparation tampering cannot forge a
  // reassuring convergence history for an otherwise unchanged shock.
  digest.add_double(jump.root_bracket_lower_compression);
  digest.add_double(jump.root_bracket_upper_compression);
  digest.add_double(jump.root_bracket_width);
  // SHK16 exposes singular-conditioning and evolutionary-branch diagnostics.
  // They are sealed just like the bracket endpoints because changing them
  // after preparation would falsify the documented reason a shock was
  // accepted or rejected even when its primitive payload was untouched.
  digest.add_bool(jump.encountered_tangential_singularity);
  digest.add_double(jump.minimum_tangential_determinant_relative);
  digest.add_double(jump.closest_tangential_determinant_relative);
  digest.add_double(jump.selected_tangential_determinant_relative);
  digest.add_double(jump.downstream_fast_mach);
  digest.add_double(jump.downstream_normal_alfven_mach);
  digest.add_bool(jump.evolutionary_fast_branch);
  add_primitive(digest,jump.upstream);
  add_primitive(digest,jump.downstream);
  digest.add_double(jump.mass_residual);
  digest.add_double(jump.normal_B_residual);
  digest.add_double(jump.electric_residual);
  digest.add_double(jump.momentum_residual);
  digest.add_double(jump.energy_residual);
  digest.add_double(jump.entropy_ratio);
}

inline void add_solar_wind(
    ConfigurationDigestBuilder& digest,
    const solarwind::PreparedState& state) noexcept {
  digest.add_double(state.V_sw_m_s);
  digest.add_double(state.T_K);
  digest.add_double(state.gamma_ad);
  digest.add_uint64(static_cast<std::uint64_t>(state.thermodynamic_closure));
  digest.add_double(state.alpha_to_proton_ratio);
  digest.add_double(state.electron_T_K);
  digest.add_double(state.alpha_T_K);
  digest.add_double(state.solar_rotation_rate_rad_s);
  digest.add_double(state.reference_sin_theta);
  digest.add_double(state.parker_source_radius_m);
  digest.add_double(state.k_AU_equatorial);
  digest.add_double(state.Br1AU_T);
  digest.add_double(state.C2);
  digest.add_double(state.C4);
  digest.add_double(state.C6);
}

inline void add_common(ConfigurationDigestBuilder& digest,
                       const core::PreparedState& state) noexcept {
  add_solar_wind(digest,state.solar_wind);
  digest.add_uint64(static_cast<std::uint64_t>(state.apex.status));
  digest.add_double(state.apex.radius_m);
  digest.add_double(state.apex.speed_m_s);
  digest.add_double(state.r0_m);
}

inline void add_region_config(ConfigurationDigestBuilder& digest,
                              const regions::Config& config) noexcept {
  digest.add_uint64(static_cast<std::uint64_t>(config.mode));
  digest.add_double(config.sheath_fraction);
  digest.add_double(config.ejecta_fraction);
  digest.add_double(config.shock_smooth_fraction);
  digest.add_double(config.leading_smooth_fraction);
  digest.add_double(config.trailing_smooth_fraction);
  digest.add_double(config.sheath_ramp_power);
  digest.add_double(config.V_sheath_LE_factor);
  digest.add_double(config.f_ME);
  digest.add_double(config.V_ME_factor);
}

inline void add_boundaries(ConfigurationDigestBuilder& digest,
                           const regions::Boundaries& boundaries) noexcept {
  digest.add_double(boundaries.R_sh_m);
  digest.add_double(boundaries.R_le_m);
  digest.add_double(boundaries.R_te_m);
  digest.add_double(boundaries.sheath_thickness_m);
  digest.add_double(boundaries.ejecta_thickness_m);
  digest.add_double(boundaries.smooth_shock_width_m);
  digest.add_double(boundaries.smooth_le_width_m);
  digest.add_double(boundaries.smooth_te_width_m);
}

inline void add_acceleration_config(
    ConfigurationDigestBuilder& digest,
    const acceleration::Config& config) noexcept {
  digest.add_uint64(static_cast<std::uint64_t>(config.mode));
  digest.add_double(config.relative_source_weight_per_area);
}

}  // namespace prepared_integrity
}  // namespace swcme
