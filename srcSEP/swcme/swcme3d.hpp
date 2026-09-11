#pragma once

#include <atomic>

#include "swcme_constants.hpp"
#include "swcme_units.hpp"
#include "swcme_defaults.hpp"
#include "swcme_status.hpp"
#include "swcme_divergence.hpp"
#include "swcme_config.hpp"
#include "swcme_regions.hpp"
#include "swcme_acceleration.hpp"
#include "swcme_kinematics.hpp"
#include "swcme_solarwind.hpp"
#include "swcme_core.hpp"
#include "swcme_shock.hpp"
#include "swcme_prepared_integrity.hpp"
#include "swcme_output.hpp"
// ============================================================================
// swcme3d.hpp
// ----------------------------------------------------------------------------
// PUBLIC API for a fast, semi-analytical Solar Wind + CME forward-shock model.
//
// What this header declares
// -------------------------
// • Public constants (AU, Rs, PI) with external linkage.
// • Parameter bundle (Params) with physically-named fields and units.
// • Time-dependent state (StepState) computed from Params at time t.
// • Simple geometric containers: ShockMesh, TriMetrics, BoxSpec.
// • Model class with:
//   - Directional shock geometry: radius and normal.
//   - Oblique ideal-MHD fast-shock classification and Rankine-Hugoniot downstream state.
//   - Field evaluators returning density (n), bulk velocity (V), magnetic
//     field (B) and divergence ∇·V at arbitrary Cartesian points.
//   - Observer-to-shock magnetic connectivity on the analytical Parker line,
//     including all cobpoint roots, selected root, path length, and local
//     production ShockState.
//   - Topologically unique shock-surface meshing + per-triangle metrics and
//     deterministic area-weighted source-cell sampling.
//   - Tecplot writers (surface, volume box, an optional box face).
//
// Units & conventions
// -------------------
// • Positions in SI meters (m). Velocities in m/s. Density in m^-3.
// • B field in Tesla (1 nT = 1e-9 T).
// • Angles in radians. Time in seconds.
// • The Sun is at the origin (0,0,0).
// • The CME apex direction is the unit vector Params::cme_dir (global frame).
//
// Physics model (brief)
// ---------------------
// • Ambient density n_up(r): Leblanc et al. (1998), scaled to n(1 AU).
// • Parker spiral upstream B (Parker 1958), evaluated in a local spherical
//   basis defined by an explicit solar-rotation axis. B1AU_nT specifies the
//   total field magnitude at a documented reference colatitude; the local
//   winding everywhere else is computed from the point's actual latitude.
// • CME apex kinematics via shared BALLISTIC / sign-aware DBM / DATA_DRIVEN modes.
// • Shock shape: sphere, self-similar ellipsoid, or finite true-SSE spherical cap.
// • Local fast-shock state from the ideal-MHD Rankine-Hugoniot conditions,
//   including complete downstream rho, p, V, and B.
// • Sheath / ejecta blends using C^1 smoothsteps and independent edge widths.
// • Tangential B amplified smoothly in the sheath (Bn continuous).
// • ∇·V uses the exact 2 V_sw/r result for SHOCK_ONLY.  FULL_ICME, whose
//   RH/sheath velocity can be non-radial and angle-dependent, uses the full
//   second-order Cartesian Jacobian trace dVx/dx+dVy/dy+dVz/dz.
//
// Output (Tecplot)
// ----------------
// VARIABLES order (used consistently across all zones):
//   1: X [m]   2: Y [m]   3: Z [m]
//   4: n [m^-3]      5..7: Vx,Vy,Vz [m/s]
//   8..10: Bx,By,Bz [T]   11: divVsw [1/s]
//   12: rc [-]        13: Vsh_n [m/s]
//   14..16: nx,ny,nz (geom normal; cell-centered in surface_cells)
//   17: area [m^2]    18: rc_mean [-]  19: Vsh_n_mean [m/s]
//   20..22: tnx,tny,tnz (reserved)     23..25: cx,cy,cz [m] (centroid)
//
//
// References (short list)
// -----------------------
// • Parker, E.N. (1958) ApJ 128:664 — Parker spiral.
// • Leblanc, Y.; Dulk, G.A.; Bougeret, J.-L. (1998) Sol. Phys. 183:165 — n(r).
// • Vršnak, B. et al. (2013) Sol. Phys. 285:295 — DBM CME kinematics.
// • Edmiston, J.P.; Kennel, C.F. (1984) J. Plasma Phys. 32:429 — oblique shocks.
// • Priest, E. (2014) Magnetohydrodynamics of the Sun, CUP.
// • Russell, C.T.; Mulligan, T. (2002) PSS 50:527 — CME sheath/ICME structure.
// • Manchester, W.B. IV et al. (2005) ApJ 622:1225 — CME-driven shocks.
//



//
// PHYSICS CONTENTS (Brief):
//  • Upstream solar wind: Parker-spiral magnetic field (Parker 1958) and
//    Leblanc et al. (1998) empirical density profile, both scaled to 1 AU.
//  • CME shock apex kinematics: Drag-Based Model (DBM; Vršnak et al. 2013).
//  • Shock geometry: Sphere / self-similar Ellipsoid / finite true-SSE spherical cap.
//  • Region structure along field lines (radial sampling from the Sun):
//      upstream → (shock) → sheath → (leading edge) → magnetic ejecta →
//      (trailing edge) → downstream ambient,
//    using C^1 smoothstep blends with independent widths at each interface.
//  • Oblique ideal-MHD fast-shock solver: explicit shock/no-shock state,
//    density compression, and conservative downstream rho,p,V,B.
//  • The immediate post-shock magnetic/velocity state is the RH solution; the
//    interior sheath relaxation remains phenomenological.
//  • Mesh generator for the shock surface + triangle metrics (area, normals,
//    centroids, rc_mean, Vsh_n_mean). All saved to Tecplot.
//  • Structured-box samplers for (n, V, B, ∇·V) with explicit ModelStatus
//    propagation; non-finite physics values are never sanitized into zeros.
//
// EFFICIENCY (Hot path):
//  • StepState caches Parker/Leblanc constants and geometry helpers per time.
//    No std::pow in loops; density via r^-2, r^-4, r^-6 with cached coeffs.
//  • Smoothstep edges use precomputed 1/(2w). Vectorization hint (GCC ivdep).
//  • Evaluate functions are re-entrant and allocation-free.
//
// UNITS (SI unless stated):
//  - Position [m];  AU and Rs provided
//  - Velocity [m s^-1]
//  - Density n [m^-3]
//  - Magnetic field B [Tesla]
//  - Divergence ∇·V [s^-1]
//  - Time [s]
//
// KEY REFERENCES:
//  - Parker, E.N. (1958), ApJ 128, 664  (Parker spiral)
//  - Leblanc, Y., Dulk, G.A., Bougeret, J.-L. (1998), Sol. Phys. 183, 165  (n(r))
//  - Vršnak, B. et al. (2013), Sol. Phys. 285, 295  (DBM apex kinematics)
//  - Edmiston, J.P., Kennel, C.F. (1984), JPP 32, 429  (oblique shocks)
//  - Priest, E. (2014), Magnetohydrodynamics of the Sun, CUP  (MHD shock theory)
//  - Russell, C.T., Mulligan, T. (2002), PSS 50, 527; Manchester, W. et al. (2005),
//    ApJ 622, 1225  (CME sheath & ejecta phenomenology)
//
// BUILD EXAMPLE:
//   g++ -std=c++17 -O3 -march=native demo3d_2.cpp swcme3d.cpp -o demo
// ============================================================================

#include <vector>
#include <cstddef>
#include <string>
#include <sstream>
#include <iomanip>

namespace swcme3d {

// ----------------------------------------------------------------------------
// Physical constants (defined in swcme3d.cpp)
// ----------------------------------------------------------------------------
extern const double AU;  // Astronomical Unit [m]
extern const double Rs;  // Solar radius [m]
extern const double PI;  // π

// ----------------------------------------------------------------------------
// Shock geometry selector
//  Sphere: Sun-centered sphere with radius r_sh(t).  This is retained as an
//          exact verification geometry and as a deliberately simple model.
//  Ellipsoid: Sun-centered, self-similarly expanding ellipsoid with axes
//             (a,b,c) aligned with (e1,e2,e3), b/a=axis_ratio_y and
//             c/a=axis_ratio_z.  All axes scale with the apex distance.
//  SSE: finite self-similar-expansion spherical cap.  The cap is the outward
//       arc of a sphere whose center lies on the CME axis; the observer ray is
//       tangent to the cap at the configured half width.  No shock surface is
//       returned outside that angular width.
//
//  ConeSSE is kept as a source-compatible alias for older callers.  Its
//  semantics are now the physically defined SSE spherical cap, NOT the old
//  R=R_apex*cos(theta)^m cosine-cap approximation.
// ----------------------------------------------------------------------------
enum class ShockShape { Sphere = 0, Ellipsoid = 1, SSE = 2, ConeSSE = SSE };

inline const char* shock_shape_name(ShockShape shape) {
  switch (shape) {
    case ShockShape::Sphere: return "SPHERE";
    case ShockShape::Ellipsoid: return "ELLIPSOID";
    case ShockShape::SSE: return "SSE";
  }
  return "UNKNOWN";
}

// ----------------------------------------------------------------------------
// Parameter pack (set at construction). All are read-only thereafter.
// Values are consumed by prepare_step(t) which builds the StepState cache.
// ----------------------------------------------------------------------------
struct Params {
  // Geometry & orientation
  ShockShape shape = ShockShape::SSE;
  double axis_ratio_y = swcme::defaults::ELLIPSOID_AXIS_RATIO_Y;               // Ellipsoid b/a (e2-axis)
  double axis_ratio_z = swcme::defaults::ELLIPSOID_AXIS_RATIO_Z;               // Ellipsoid c/a (e3-axis)
  // Angular half width of the finite SSE cap [rad].  For the science
  // geometry the supported range is 0 < lambda <= pi/2.  The surface exists
  // only for directions whose angular separation from cme_dir is <= lambda.
  double half_width_rad = swcme::defaults::SSE_HALF_WIDTH_RAD;

  // DEPRECATED compatibility field.  Older ConeSSE code used this exponent in
  // R(theta)=R_apex*cos(theta)^m and also applied it a second time to the flank
  // speed.  The corrected SSE geometry derives both radius and normal speed
  // from self-similar spherical-cap geometry, so this parameter is ignored
  // whenever shape==SSE/ConeSSE.  It remains in Params only to avoid breaking
  // existing input/source code while callers migrate to ShockShape::SSE.
  double flank_slowdown_m = 1.0;

  double cme_dir[3] = {1,0,0};             // Global unit vector for apex direction (e1)

  // Solar-rotation axis used to construct the local Parker spherical basis.
  // The vector is normalized once in prepare_step().  Keeping this axis
  // explicit is essential in 3-D: the Parker azimuthal direction is
  // e_phi ∝ Omega_hat × e_r and the local winding scales with
  // sin(theta)=|Omega_hat × e_r|.  The previous implementation implicitly
  // assumed +Z and used one global latitude factor for every point.
  double solar_rotation_axis[3] = {0,0,1};

  // Solar angular rotation rate used by BOTH the Parker field and analytical
  // connectivity mapping.  The default is the SWCME Carrington/sidereal model
  // convention.  Keeping it explicit allows the Omega->0 radial-field limit to
  // be verified without altering global constants or introducing a separate
  // test-only field-line equation.
  double solar_rotation_rate_rad_s = swcme::defaults::SOLAR_ROTATION_RATE_RAD_S;

  // Legacy/reference normalization latitude for B1AU_nT ONLY.  This value no
  // longer controls the local 3-D Parker pitch.  It specifies the sine of the
  // colatitude at which B1AU_nT is interpreted as the total |B| at 1 AU.
  // Thus existing inputs that used sin_theta=1 for an ecliptic reference keep
  // their 1-AU normalization, while every evaluated point now uses its own
  // geometrically correct local sin(theta).  A future API cleanup may rename
  // this field once backward compatibility is no longer required.
  double sin_theta  = swcme::defaults::PARKER_REFERENCE_SIN_THETA;

  // CME/shock-apex kinematics.  The same shared common implementation is
  // consumed by swcme1d, so identical inputs produce identical apex radius
  // and speed in both dimensional interfaces.  DBM remains the default mode.
  swcme::kinematics::Mode kinematics_mode = swcme::defaults::KINEMATICS_MODE;
  double r0_Rs       = swcme::defaults::DBM_R0_RS;               // DBM/ballistic reference radius [Rs]
  double V0_sh_kms   = swcme::defaults::V0_SH_KMS;               // initial/reference shock speed [km/s]
  double V_sw_kms    = swcme::defaults::V_SW_KMS;                // ambient SW speed [km/s]
  double Gamma_kmInv = swcme::defaults::DBM_GAMMA_KM_INV;               // DBM drag parameter Γ [km^-1], >=0

  // DATA_DRIVEN mode: strictly increasing times [s] and nondecreasing apex
  // radii [Rs].  A monotone PCHIP passes exactly through every knot and its
  // derivative is used as V_apex.  The default policy rejects out-of-range
  // queries; ballistic endpoint continuation must be requested explicitly.
  std::vector<double> data_time_s;
  std::vector<double> data_radius_Rs;
  swcme::kinematics::ExtrapolationPolicy data_extrapolation =
      swcme::defaults::DATA_EXTRAPOLATION;

  // Ambient SW scalings at 1 AU
  double n1AU_cm3 = swcme::defaults::N1AU_CM3;                   // upstream density at 1 AU [cm^-3]
  double B1AU_nT  = swcme::defaults::B1AU_TOTAL_NT;                   // |B|(1 AU) [nT] at reference sin_theta above
  double T_K      = swcme::defaults::T_K;                 // proton temperature [K]
  double gamma_ad = swcme::defaults::GAMMA_AD;               // adiabatic index

  // Region mode.  The canonical science default is SHOCK_ONLY: the
  // Parker/Leblanc background is left unchanged and the shock is used for
  // geometry/connectivity/source bookkeeping. FULL_ICME is an explicit
  // phenomenological diagnostic override.
  swcme::regions::Mode region_mode = swcme::defaults::REGION_MODE;

  // Mutually exclusive shock-acceleration representation. SOURCE is the
  // canonical controlled-SEP default and is validated with SHOCK_ONLY.
  // RESOLVED_COMPRESSION is an explicit FULL_ICME diagnostic override.
  swcme::acceleration::Mode shock_acceleration_mode =
      swcme::defaults::ACCELERATION_MODE;
  double relative_source_weight_per_area = swcme::defaults::RELATIVE_SOURCE_WEIGHT_PER_AREA;

  // Self-similar radial thicknesses. Values are AU at a 1-AU shock and thus
  // dimensionless fractions when scaled to each local shock-surface radius.
  double sheath_thick_AU_at1AU = swcme::defaults::SHEATH_THICK_AU_AT_1AU;      // sheath thickness at 1 AU [AU]
  double ejecta_thick_AU_at1AU = swcme::defaults::EJECTA_THICK_AU_AT_1AU;      // ejecta thickness at 1 AU [AU]

  // Edge smoothing (C^1 smoothstep) widths, scale ∝ local R_sh.  The shock
  // width is used only by RESOLVED_COMPRESSION; SOURCE is paired with
  // SHOCK_ONLY and therefore contributes no resolved compression to V(x,t).
  double edge_smooth_shock_AU_at1AU = swcme::defaults::EDGE_SMOOTH_SHOCK_AU_AT_1AU;
  double edge_smooth_le_AU_at1AU    = swcme::defaults::EDGE_SMOOTH_LE_AU_AT_1AU; // sheath→ejecta leading edge
  double edge_smooth_te_AU_at1AU    = swcme::defaults::EDGE_SMOOTH_TE_AU_AT_1AU; // ejecta→downstream trailing edge

  // Target speeds in regions (relative to upstream V_sw)
  double V_sheath_LE_factor = swcme::defaults::V_SHEATH_LE_FACTOR;        // sheath speed at LE = factor*V_sw
  double V_ME_factor        = swcme::defaults::V_ME_FACTOR;        // ejecta speed (bulk) = factor*V_sw

  // Sheath shaping
  double sheath_ramp_power  = swcme::defaults::SHEATH_RAMP_POWER;        // ≥1, steeper density jump near shock
  // Deprecated compatibility parameter: no longer used by shock or region
  // physics.  The MHD RH solver alone determines physical compression.
  double sheath_comp_floor  = swcme::defaults::SHEATH_COMP_FLOOR_COMPAT;

  // Ejecta density factor (relative to upstream)
  double f_ME = swcme::defaults::F_ME;                      // n_ejecta = f_ME * n_up
};

// Deterministic complete 3-D resolved-configuration record.  It includes
// inactive/deprecated compatibility fields as well as the active science
// options so event-specific overrides remain auditable in validation output.
inline std::string resolved_configuration_manifest(const Params& p) {
  std::ostringstream out;
  out << std::setprecision(17) << std::scientific;
  out << "swcme_config_version=" << swcme::defaults::CONFIG_VERSION << '\n';
  out << "model=3D\n";
  out << "frame=" << swcme::defaults::FRAME_NAME << '\n';
  out << "model_scope=" << swcme::defaults::model_scope_name(
      swcme::defaults::model_scope(p.region_mode, p.shock_acceleration_mode)) << '\n';
  out << "parker_normalization="
      << swcme::defaults::PARKER_NORMALIZATION_CONVENTION << '\n';
  out << "parker_radial_polarity=" << swcme::defaults::PARKER_RADIAL_POLARITY << '\n';
  out << "shape=" << shock_shape_name(p.shape) << '\n';
  out << "axis_ratio_y=" << p.axis_ratio_y << '\n';
  out << "axis_ratio_z=" << p.axis_ratio_z << '\n';
  out << "half_width_rad=" << p.half_width_rad << '\n';
  out << "flank_slowdown_m=" << p.flank_slowdown_m << '\n';
  for (int i=0; i<3; ++i) out << "cme_dir[" << i << "]=" << p.cme_dir[i] << '\n';
  for (int i=0; i<3; ++i) out << "solar_rotation_axis[" << i << "]=" << p.solar_rotation_axis[i] << '\n';
  out << "solar_rotation_rate_rad_s=" << p.solar_rotation_rate_rad_s << '\n';
  out << "sin_theta=" << p.sin_theta << '\n';
  out << "kinematics_mode=" << swcme::defaults::kinematics_mode_name(p.kinematics_mode) << '\n';
  out << "r0_Rs=" << p.r0_Rs << '\n';
  out << "V0_sh_kms=" << p.V0_sh_kms << '\n';
  out << "V_sw_kms=" << p.V_sw_kms << '\n';
  out << "Gamma_kmInv=" << p.Gamma_kmInv << '\n';
  out << "data_extrapolation="
      << swcme::defaults::extrapolation_policy_name(p.data_extrapolation) << '\n';
  out << "data_time_s.count=" << p.data_time_s.size() << '\n';
  for (std::size_t i=0; i<p.data_time_s.size(); ++i)
    out << "data_time_s[" << i << "]=" << p.data_time_s[i] << '\n';
  out << "data_radius_Rs.count=" << p.data_radius_Rs.size() << '\n';
  for (std::size_t i=0; i<p.data_radius_Rs.size(); ++i)
    out << "data_radius_Rs[" << i << "]=" << p.data_radius_Rs[i] << '\n';
  out << "n1AU_cm3=" << p.n1AU_cm3 << '\n';
  out << "B1AU_nT=" << p.B1AU_nT << '\n';
  out << "T_K=" << p.T_K << '\n';
  out << "gamma_ad=" << p.gamma_ad << '\n';
  out << "region_mode=" << swcme::defaults::region_mode_name(p.region_mode) << '\n';
  out << "shock_acceleration_mode="
      << swcme::defaults::acceleration_mode_name(p.shock_acceleration_mode) << '\n';
  out << "relative_source_weight_per_area=" << p.relative_source_weight_per_area << '\n';
  out << "sheath_thick_AU_at1AU=" << p.sheath_thick_AU_at1AU << '\n';
  out << "ejecta_thick_AU_at1AU=" << p.ejecta_thick_AU_at1AU << '\n';
  out << "edge_smooth_shock_AU_at1AU=" << p.edge_smooth_shock_AU_at1AU << '\n';
  out << "edge_smooth_le_AU_at1AU=" << p.edge_smooth_le_AU_at1AU << '\n';
  out << "edge_smooth_te_AU_at1AU=" << p.edge_smooth_te_AU_at1AU << '\n';
  out << "V_sheath_LE_factor=" << p.V_sheath_LE_factor << '\n';
  out << "V_ME_factor=" << p.V_ME_factor << '\n';
  out << "sheath_ramp_power=" << p.sheath_ramp_power << '\n';
  out << "sheath_comp_floor=" << p.sheath_comp_floor << '\n';
  out << "f_ME=" << p.f_ME << '\n';
  return out.str();
}

// Deterministically fingerprint the complete resolved 3-D configuration for
// PST03.  The versioned schema tag and fixed serialization order make values
// reproducible across processes; no addresses, struct padding, or locale-
// dependent text formatting enter the digest.
inline swcme::ConfigurationDigest configuration_digest(
    const Params& p) noexcept {
  swcme::ConfigurationDigestBuilder digest;
  digest.add_string("SWCME_CONFIGURATION_DIGEST_V1");
  digest.add_string("3D");
  digest.add_uint64(static_cast<std::uint64_t>(swcme::defaults::CONFIG_VERSION));
  digest.add_string(swcme::defaults::FRAME_NAME);
  digest.add_string(swcme::defaults::PARKER_NORMALIZATION_CONVENTION);
  digest.add_uint64(static_cast<std::uint64_t>(
      static_cast<std::int64_t>(swcme::defaults::PARKER_RADIAL_POLARITY)));
  digest.add_string("PROTON_ONLY_THERMAL_PRESSURE_CLOSURE");

  // Hash every public Params field in declaration order, including currently
  // inactive/deprecated compatibility values.  This matches the completeness
  // promise of resolved_configuration_manifest() and prevents a mode change
  // from reviving a parameter that was omitted from state ownership checks.
  digest.add_uint64(static_cast<std::uint64_t>(p.shape));
  digest.add_double(p.axis_ratio_y); digest.add_double(p.axis_ratio_z);
  digest.add_double(p.half_width_rad); digest.add_double(p.flank_slowdown_m);
  for (double value : p.cme_dir) digest.add_double(value);
  for (double value : p.solar_rotation_axis) digest.add_double(value);
  digest.add_double(p.solar_rotation_rate_rad_s);
  digest.add_double(p.sin_theta);
  digest.add_uint64(static_cast<std::uint64_t>(p.kinematics_mode));
  digest.add_double(p.r0_Rs); digest.add_double(p.V0_sh_kms);
  digest.add_double(p.V_sw_kms); digest.add_double(p.Gamma_kmInv);
  digest.add_uint64(static_cast<std::uint64_t>(p.data_time_s.size()));
  for (double value : p.data_time_s) digest.add_double(value);
  digest.add_uint64(static_cast<std::uint64_t>(p.data_radius_Rs.size()));
  for (double value : p.data_radius_Rs) digest.add_double(value);
  digest.add_uint64(static_cast<std::uint64_t>(p.data_extrapolation));
  digest.add_double(p.n1AU_cm3); digest.add_double(p.B1AU_nT);
  digest.add_double(p.T_K); digest.add_double(p.gamma_ad);
  digest.add_uint64(static_cast<std::uint64_t>(p.region_mode));
  digest.add_uint64(static_cast<std::uint64_t>(p.shock_acceleration_mode));
  digest.add_double(p.relative_source_weight_per_area);
  digest.add_double(p.sheath_thick_AU_at1AU);
  digest.add_double(p.ejecta_thick_AU_at1AU);
  digest.add_double(p.edge_smooth_shock_AU_at1AU);
  digest.add_double(p.edge_smooth_le_AU_at1AU);
  digest.add_double(p.edge_smooth_te_AU_at1AU);
  digest.add_double(p.V_sheath_LE_factor); digest.add_double(p.V_ME_factor);
  digest.add_double(p.sheath_ramp_power); digest.add_double(p.sheath_comp_floor);
  digest.add_double(p.f_ME);
  return digest.value();
}

// Validate the complete 3-D public parameter pack before any basis
// normalization or field/shock calculation.  Common plasma, kinematic, and
// region rules are delegated to swcme_config.hpp; this wrapper adds only the
// geometry/vector rules that have no 1-D analogue.
inline swcme::config::ValidationResult validate_params(const Params& p) {
  swcme::config::CommonConfigView view;
  view.V_sw_kms=p.V_sw_kms; view.n1AU_cm3=p.n1AU_cm3;
  view.B1AU_nT=p.B1AU_nT; view.T_K=p.T_K; view.gamma_ad=p.gamma_ad;
  view.sin_theta=p.sin_theta; view.kinematics_mode=p.kinematics_mode;
  view.r0_Rs=p.r0_Rs; view.V0_sh_kms=p.V0_sh_kms;
  view.Gamma_kmInv=p.Gamma_kmInv; view.data_time_s=&p.data_time_s;
  view.data_radius_Rs=&p.data_radius_Rs;
  view.region_mode=p.region_mode;
  view.acceleration_mode=p.shock_acceleration_mode;
  view.relative_source_weight_per_area=p.relative_source_weight_per_area;
  view.sheath_thick_AU_at1AU=p.sheath_thick_AU_at1AU;
  view.ejecta_thick_AU_at1AU=p.ejecta_thick_AU_at1AU;
  view.edge_smooth_shock_AU_at1AU=p.edge_smooth_shock_AU_at1AU;
  view.edge_smooth_le_AU_at1AU=p.edge_smooth_le_AU_at1AU;
  view.edge_smooth_te_AU_at1AU=p.edge_smooth_te_AU_at1AU;
  view.sheath_comp_floor=p.sheath_comp_floor;
  view.sheath_ramp_power=p.sheath_ramp_power;
  view.V_sheath_LE_factor=p.V_sheath_LE_factor;
  view.f_ME=p.f_ME; view.V_ME_factor=p.V_ME_factor;

  swcme::config::ValidationResult out=swcme::config::validate_common(view);
  swcme::config::require_nonzero_vector(out,"cme_dir",p.cme_dir);
  swcme::config::require_nonzero_vector(out,"solar_rotation_axis",
                                        p.solar_rotation_axis);
  swcme::config::require_nonnegative(out,"solar_rotation_rate_rad_s",
                                     p.solar_rotation_rate_rad_s);

  if (!swcme::config::finite(p.flank_slowdown_m)) {
    out.add("flank_slowdown_m",swcme::config::Code::NonFinite,
            p.flank_slowdown_m,"deprecated compatibility value must be finite");
  }

  if (p.shape==ShockShape::Ellipsoid) {
    swcme::config::require_positive(out,"axis_ratio_y",p.axis_ratio_y);
    swcme::config::require_positive(out,"axis_ratio_z",p.axis_ratio_z);
  }
  if (p.shape==ShockShape::SSE) {
    swcme::config::require_range(out,"half_width_rad",p.half_width_rad,
                                 0.0,0.5*PI,false,true);
  }
  return out;
}

// ----------------------------------------------------------------------------
// Time-dependent state for a given time t (returned by prepare_step).
// CACHES everything required by hot loops (n,V,B,divV evaluators).
// ----------------------------------------------------------------------------
struct StepState {
  // Process-local identity of the exact Model instance that prepared this
  // cache.  Zero denotes an unprepared/default state.  The token is verified
  // before geometry or physics uses any cached value, preventing a foreign
  // cache from being combined with this model's Params.
  swcme::ModelIdentity owner_model_identity = 0;

  // Immutable snapshot of the resolved Params and global conventions that
  // produced this cache.  The model validates it before using any state data.
  swcme::ConfigurationDigest configuration_digest = 0;

  // Canonical dimensionality-independent prepared state.  3-D keeps several
  // legacy mirror fields below because they are part of the current public
  // StepState API, but their values are copied from this common state instead
  // of being recomputed independently.
  swcme::core::PreparedState common;

  // Common region contract. apex_regions mirrors the nominal boundary radii
  // for diagnostics/output; field evaluators recompute the same self-similar
  // boundaries from each local shock radius so finite-SSE/ellipsoid flanks do
  // not inherit an apex-sized absolute layer thickness.
  swcme::regions::Config region_config;
  swcme::regions::Boundaries apex_regions;
  swcme::acceleration::Config acceleration_config;

  // Apex-aligned orthonormal frame (e1 ≡ apex dir; e2,e3 transverse)
  double e1[3], e2[3], e3[3];

  // Apex kinematics & scale.  kinematics_mode records the common solver that
  // generated the step state and is useful when archiving validation metadata.
  swcme::kinematics::Mode kinematics_mode = swcme::defaults::KINEMATICS_MODE;
  double time_s = 0.0;     // time since launch [s]
  double r_sh_m = 0.0;    // shock apex radius [m]
  double V_sh_ms= 0.0;    // shock apex speed [m/s]
  double a_m    = 0.0;    // reference size (== r_sh_m)

  // Region extents and smoothing widths (in meters)
  double dr_sheath_m = 0.0;   // sheath thickness
  double dr_me_m     = 0.0;   // ejecta thickness
  double w_shock_m   = 0.0;   // smooth width at shock
  double w_le_m      = 0.0;   // smooth width at leading edge
  double w_te_m      = 0.0;   // smooth width at trailing edge
  double r_le_m      = 0.0;   // radius of leading edge (= r_sh - dr_sheath)
  double r_te_m      = 0.0;   // radius of trailing edge (= r_le - dr_me)

  // Target speeds
  double V_sheath_LE_ms = 0.0; // sheath speed at LE
  double V_ME_ms        = 0.0; // ejecta bulk speed
  double V_dn_ms        = 0.0; // downstream ambient (== V_sw)

  // Convenience
  // Apex shock diagnostic.  has_shock is deliberately independent of the
  // existence of the geometric CME front: a surface may exist while its
  // normal motion is sub-fast and therefore does not support a physical fast
  // shock.  rc is exactly 1 whenever has_shock is false.
  bool has_shock = false;
  double rc = 1.0;             // diagnostic: apex density compression
  double inv_dr_sheath = 0.0;  // 1 / dr_sheath
  double rc_floor = 1.0;       // deprecated compatibility mirror; always 1

  // ===== CACHES FOR EFFICIENCY (constant during this step) =====
  // Ambient wind speed
  double V_sw_ms = 0.0;

  // Leblanc density, SI coefficients with r-powers absorbed:
  // n(r) = C2 r^-2 + C4 r^-4 + C6 r^-6 [m^-3]
  double C2 = 0.0, C4 = 0.0, C6 = 0.0;

  // Precompute 1/(2w) to reduce divides in smoothstep arguments
  double inv2w_sh = 0.0, inv2w_le = 0.0, inv2w_te = 0.0;

  // Parker spiral cache for the current step.
  //
  // solar_axis_hat is the normalized solar-rotation axis.  k_AU is now the
  // EQUATORIAL pitch coefficient Ω*AU/V_sw; the local Parker pitch is
  // k_AU * r_AU * sin(theta_local), where sin(theta_local) is computed from
  // solar_axis_hat and the point's radial direction.  This deliberately
  // removes the old incorrect assumption of one fixed latitude throughout
  // the 3-D domain.
  double solar_axis_hat[3] = {0.0,0.0,1.0};
  double solar_rotation_rate_rad_s = swcme::defaults::SOLAR_ROTATION_RATE_RAD_S;
  double k_AU = 0.0;
  double Br1AU_T = 0.0;   // Br at 1 AU [T], normalized using Params::sin_theta as reference latitude

  // Ellipsoid helpers: a,b,c and 1/a^2,1/b^2,1/c^2
  double a_e=0.0, b_e=0.0, c_e=0.0;
  double inv_a2=0.0, inv_b2=0.0, inv_c2=0.0;

  // Finite SSE spherical-cap helpers.  For apex distance R_a and half width
  // lambda, the generating sphere has center distance
  //   c = R_a/(1+sin(lambda))
  // and radius
  //   a = R_a*sin(lambda)/(1+sin(lambda)).
  // These are cached because every directional geometry query reuses them.
  double sin_half_width=0.0;
  double cos_half_width=1.0;
  double sse_center_m=0.0;
  double sse_radius_m=0.0;

  // Read-only exposure supports diagnostics and copy/move validation without
  // allowing compatibility code to bless a modified state record.
  swcme::ConfigurationDigest integrity_digest() const noexcept {
    return integrity_digest_;
  }

 private:
  swcme::ConfigurationDigest integrity_digest_ = 0;
  friend class Model;
};

// Compute the PST06 seal over every 3-D cached value and public compatibility
// mirror.  Arrays are serialized element-by-element and the private seal is
// excluded, avoiding ABI padding and self-referential hashing.
inline swcme::ConfigurationDigest prepared_state_integrity(
    const StepState& state) noexcept {
  swcme::ConfigurationDigestBuilder digest(false);
  digest.add_string("SWCME_3D_PREPARED_STATE_INTEGRITY_V1");
  digest.add_uint64(state.owner_model_identity);
  digest.add_uint64(state.configuration_digest);
  swcme::prepared_integrity::add_common(digest,state.common);
  swcme::prepared_integrity::add_region_config(digest,state.region_config);
  swcme::prepared_integrity::add_boundaries(digest,state.apex_regions);
  swcme::prepared_integrity::add_acceleration_config(
      digest,state.acceleration_config);
  for (double value : state.e1) digest.add_double(value);
  for (double value : state.e2) digest.add_double(value);
  for (double value : state.e3) digest.add_double(value);
  digest.add_uint64(static_cast<std::uint64_t>(state.kinematics_mode));
  digest.add_double(state.time_s);
  digest.add_double(state.r_sh_m);
  digest.add_double(state.V_sh_ms);
  digest.add_double(state.a_m);
  digest.add_double(state.dr_sheath_m);
  digest.add_double(state.dr_me_m);
  digest.add_double(state.w_shock_m);
  digest.add_double(state.w_le_m);
  digest.add_double(state.w_te_m);
  digest.add_double(state.r_le_m);
  digest.add_double(state.r_te_m);
  digest.add_double(state.V_sheath_LE_ms);
  digest.add_double(state.V_ME_ms);
  digest.add_double(state.V_dn_ms);
  digest.add_bool(state.has_shock);
  digest.add_double(state.rc);
  digest.add_double(state.inv_dr_sheath);
  digest.add_double(state.rc_floor);
  digest.add_double(state.V_sw_ms);
  digest.add_double(state.C2);
  digest.add_double(state.C4);
  digest.add_double(state.C6);
  digest.add_double(state.inv2w_sh);
  digest.add_double(state.inv2w_le);
  digest.add_double(state.inv2w_te);
  for (double value : state.solar_axis_hat) digest.add_double(value);
  digest.add_double(state.solar_rotation_rate_rad_s);
  digest.add_double(state.k_AU);
  digest.add_double(state.Br1AU_T);
  digest.add_double(state.a_e);
  digest.add_double(state.b_e);
  digest.add_double(state.c_e);
  digest.add_double(state.inv_a2);
  digest.add_double(state.inv_b2);
  digest.add_double(state.inv_c2);
  digest.add_double(state.sin_half_width);
  digest.add_double(state.cos_half_width);
  digest.add_double(state.sse_center_m);
  digest.add_double(state.sse_radius_m);
  return digest.value();
}


// ----------------------------------------------------------------------------
// Complete local shock state returned by Model::shock_state_direction().
//
// This structure is the public bridge between geometry and shock physics.  It
// exposes both the upstream and downstream primitive MHD states at the actual
// shock surface, rather than forcing callers to reconstruct a jump from a
// compression proxy.  surface_exists and has_shock are intentionally separate:
// a finite CME surface can be present even when the relative normal speed is
// sub-fast and no physical shock/downstream jump exists.
// ----------------------------------------------------------------------------
struct LocalShockState {
  bool surface_exists = false;
  bool has_shock = false;
  bool solver_converged = true;
  swcme::ModelStatus status;

  double Rdir_m = 0.0;
  double normal[3] = {0.0,0.0,0.0};
  double Vsh_n_m_s = 0.0;
  double theta_Bn_rad = 0.0;
  double fast_speed_m_s = 0.0;
  double fast_mach = 0.0;
  double compression = 1.0;

  double upstream_n_m3 = 0.0;
  double downstream_n_m3 = 0.0;
  swcme::shock::PrimitiveState upstream;
  swcme::shock::PrimitiveState downstream;

  double mass_residual = 0.0;
  double normal_B_residual = 0.0;
  double electric_residual = 0.0;
  double momentum_residual = 0.0;
  double energy_residual = 0.0;
  double entropy_ratio = 1.0;
};

// ----------------------------------------------------------------------------
// Observer-to-shock magnetic connectivity state.
//
// The baseline upstream magnetic field is the analytical Parker spiral used by
// the production field evaluator.  A magnetic connection is therefore found
// by tracing the observer's Parker field line inward and intersecting it with
// the same production shock geometry used by shape_radius_normal().  Keeping
// the cobpoint calculation on the production geometry/shock APIs prevents the
// connectivity code from drifting to a second, inconsistent copy of the shock
// equations.
// ----------------------------------------------------------------------------
enum class ConnectivityStatus {
  Connected,
  Disconnected,
  InvalidObserver,
  InvalidConfiguration
};

struct ConnectivityOptions {
  // Inner radius of the field-line search.  The default is the lower radial
  // limit historically used by the analytical solar-wind model.  Science
  // applications may raise this (for example to a data-driven/DBM handoff
  // radius) without changing the connectivity algorithm.
  double inner_radius_m = 1.05 * swcme::constants::SOLAR_RADIUS_M;

  // Minimum number of radial scan intervals.  The implementation can increase
  // this automatically when a tightly wound Parker line requires finer phase
  // sampling.  The scan finds every sign-changing intersection and seeds
  // tangent-root searches; final roots are then refined independently.
  std::size_t scan_intervals = 1024;

  // Radial convergence tolerance for an intersection.  1e-9 AU is well below
  // the 1e-8 AU acceptance target in the validation plan while remaining many
  // orders above floating-point spacing at heliospheric radii.
  double radius_tolerance_m = 1.0e-9 * swcme::constants::AU_M;

  // Maximum |r - R_shock(direction)| accepted as a geometrical root.  A
  // slightly looser value than the bisection tolerance is used so tangent
  // roots (which need not change sign) can be identified robustly.
  double surface_residual_tolerance_m = 5.0e-9 * swcme::constants::AU_M;
};

struct ConnectivityRoot {
  // Radial coordinate and Cartesian cobpoint on the observer Parker line.
  double radius_m = 0.0;
  double position_m[3] = {0.0,0.0,0.0};

  // Arc length from this cobpoint outward to the observer along the analytical
  // Parker field line.  This is the length required by field-aligned SEP
  // transport; it is intentionally not the shorter radial separation.
  double path_length_m = 0.0;

  // Signed geometric residual r-R_shock at the final root.
  double surface_residual_m = 0.0;

  // Complete local shock state evaluated by the production shock solver at the
  // cobpoint direction.  has_shock may be false even though the geometrical CME
  // front is intersected; this preserves the distinction between a front and a
  // physical fast shock introduced by the shock-physics correction.
  LocalShockState shock;
};

struct ConnectivityState {
  ConnectivityStatus status = ConnectivityStatus::Disconnected;
  bool connected = false;

  double observer_position_m[3] = {0.0,0.0,0.0};
  double observer_radius_m = 0.0;

  // All intersections are retained in increasing radial order.  The default
  // physical cobpoint is the outermost root (largest radius), i.e. the first
  // shock surface encountered when tracing inward from the observer.  Retaining
  // all roots makes the selection rule auditable and supports future non-convex
  // geometries without silently discarding intersections.
  std::vector<ConnectivityRoot> roots;
  std::size_t selected_root = 0;
};

struct ConnectivityHistorySample {
  double time_s = 0.0;
  ConnectivityState connectivity;
};

// ----------------------------------------------------------------------------
// Shock surface mesh and per-triangle metrics.
//
// Topology contract (Fix 12)
// --------------------------
// The mesh no longer stores a rectangular theta-phi array with duplicated
// phi=0/2pi seam vertices or a complete ring of coincident apex vertices.
// Instead it is an explicitly triangular manifold:
//   * one unique apex vertex;
//   * nPhi unique vertices on every non-polar ring;
//   * periodic ring connectivity implemented by index wrap, not duplicate nodes;
//   * one unique rear pole for closed Sphere/Ellipsoid surfaces; and
//   * one physical outer boundary ring for a finite SSE cap.
//
// `n_theta_intervals` and `n_phi` describe the requested angular resolution and
// are recorded so consumers can audit the construction without reverse
// engineering node counts.  `closed_surface` distinguishes a closed sphere/
// ellipsoid from an open finite-width SSE cap.  Connectivity remains 1-based
// for Tecplot compatibility.
// ----------------------------------------------------------------------------
struct ShockMesh {
  // Nodal fields (size Nv)
  std::vector<double> x, y, z;                   // vertex positions [m]
  std::vector<double> n_hat_x, n_hat_y, n_hat_z; // outward analytic normals
  std::vector<double> rc;                        // physical density compression
  std::vector<double> Vsh_n;                     // normal shock speed [m/s]

  // Connectivity (1-based triangle indices; all three arrays have size Ne).
  std::vector<int> tri_i, tri_j, tri_k;

  // Construction metadata. These fields do not participate in Tecplot output.
  std::size_t n_theta_intervals = 0;
  std::size_t n_phi = 0;
  bool closed_surface = false;
};

struct TriMetrics {
  // Per-triangle (cell-centered) metrics (size Ne)
  std::vector<double> area;       // [m^2]
  std::vector<double> nx, ny, nz; // outward unit normal from triangle winding
  std::vector<double> cx, cy, cz; // centroid [m]
  std::vector<double> rc_mean;    // mean of nodal rc
  std::vector<double> Vsh_n_mean; // mean of nodal Vsh_n [m/s]
};

// OUT06 output contract: a ShockMesh must contain at least three vertices and
// one triangle; all parallel arrays must have exact Nv/Ne sizes; nodal values
// must be finite with unit normals, rc>=1, and Vsh_n>=0; connectivity must be
// distinct, one-based, and in range; and every triangle must pass the canonical
// degeneracy/orientation checks.  TriMetrics may be entirely empty to request
// canonical computation, or complete and consistent with that computation.

// Deterministic area CDF used by shock-source samplers.  The library does not
// own a random-number generator: callers supply a U[0,1) variate to
// sample_triangle_by_area().  This keeps reproducibility under the caller's
// control while ensuring every source adapter uses the same physical area
// weighting.
struct AreaSamplingTable {
  std::vector<double> cumulative_probability; // strictly increasing; last=1
  double total_area_m2 = 0.0;
};

// ----------------------------------------------------------------------------
// Structured volume box specification for sampling fields in a region.
// The box is centered at (cx,cy,cz), with finite nonnegative half-sizes
// (hx,hy,hz), and samples a regular grid Ni×Nj×Nk with at least two points on
// every axis.  OUT04 additionally requires representable bounds/spans and a
// point-count product that fits size_t.  You may build a validated apex-aligned
// default via Model::default_apex_box(...).
// ----------------------------------------------------------------------------
struct BoxSpec {
  double cx=0, cy=0, cz=0;   // box center [m]
  double hx=0, hy=0, hz=0;   // half-sizes [m]
  int Ni=16, Nj=16, Nk=16;   // resolution (≥2 per axis)
};

// ----------------------------------------------------------------------------
// Model class
// ----------------------------------------------------------------------------
class Model {
public:
  explicit Model(const Params&);

  // Copies are independent model owners even when their Params compare equal.
  // The implementation assigns a fresh identity instead of copying the source
  // token, which makes PST02's instance boundary explicit.
  Model(const Model&);
  Model& operator=(const Model&);

  swcme::ModelIdentity model_identity() const noexcept {
    return model_identity_;
  }

  // Return a separately owned replacement rather than mutating an instance
  // that may already have issued prepared states.  The new model starts in
  // the configuration phase and freezes after its own first successful step.
  Model reconfigured(const Params& p) const { return Model(p); }

  bool configuration_locked() const noexcept {
    return configuration_locked_.load(std::memory_order_acquire);
  }

  // Side-effect-free ownership guard used by direct callers and the SEP
  // adapters.  Call it before modifying outputs so mismatch handling remains
  // transactional and diagnostics retain both the expected and supplied IDs.
  swcme::ModelStatus validate_prepared_state(
      const StepState& S, const char* context) const noexcept {
    const swcme::ConfigurationDigest current=configuration_digest(P_);
    // Instance provenance remains the primary error.  Configuration digests
    // are still attached so a foreign-state diagnostic distinguishes equal
    // models from gamma/geometry/Parker/kinematic/region mismatches.
    if (S.owner_model_identity!=model_identity_)
      return swcme::ModelStatus::state_model_mismatch(
          context,model_identity_,S.owner_model_identity,current,
          S.configuration_digest,true);
    if (S.configuration_digest!=current)
      return swcme::ModelStatus::state_configuration_mismatch(
          context,current,S.configuration_digest);
    const swcme::ConfigurationDigest computed=prepared_state_integrity(S);
    if (computed!=S.integrity_digest())
      return swcme::ModelStatus::stale_prepared_state(
          context,S.integrity_digest(),computed);
    return swcme::ModelStatus::success();
  }

  // Side-effect-free validation entry point used by CFG01 and prepare_step().
  // Callers can inspect all invalid fields before starting a simulation.
  swcme::config::ValidationResult validate() const { return validate_params(P_); }

  // Scope is derived from the region/acceleration pair rather than stored as
  // another independent option.  This prevents metadata from claiming a
  // controlled SOURCE run while the transport field actually contains a
  // resolved FULL_ICME compression.
  swcme::defaults::ModelScope model_scope() const {
    return swcme::defaults::model_scope(P_.region_mode,
                                        P_.shock_acceleration_mode);
  }

  // Observer-local scope check. For the default CONTROLLED_SEP_PRE_SHOCK
  // configuration, the analytical Parker background is in scope only until
  // the finite shock front reaches that observer direction/radius. If a finite
  // SSE cap does not intersect the observer ray, the Parker background remains
  // in scope on that ray.
  swcme::defaults::ObserverScopeStatus observer_scope_status(
      const StepState& S, const double observer_m[3]) const;

  // Build time-dependent state (and caches) for time t_s [s].
  StepState prepare_step(double t_s) const;

  // Shock radius along unit direction u=(ux,uy,uz) and outward normal.
  // Returns true only when the selected geometry intersects the outward ray.
  // Sphere and ellipsoid always return true; the finite SSE cap returns false
  // outside its configured half width.  When false, Rdir_m and n_hat are set
  // to zero and must not be interpreted as a physical shock surface.
  bool shape_radius_normal(const StepState& S,
                           double ux,double uy,double uz,
                           double& Rdir_m,double n_hat[3]) const;

  // Evaluate the complete ideal-MHD shock state along unit direction u.
  // The upstream plasma is ALWAYS sampled at the physical shock position
  // Rdir*u; it is never sampled at an arbitrary query point.  The return value
  // is false only when the selected finite geometry has no surface along u.
  // Checked form distinguishes NO_SURFACE from numerical failure.  New
  // transport/AMPS adapters should prefer this method so failures can be
  // propagated without exceptions.  The source-compatible bool wrapper below
  // throws on numerical failure and returns false only for NO_SURFACE.
  swcme::ModelStatus shock_state_direction_checked(
      const StepState& S, const double u[3], LocalShockState& state) const;
  bool shock_state_direction(const StepState& S, const double u[3],
                             LocalShockState& state) const;

  // Convert one directional physical shock state into the single selected
  // acceleration representation.  The checked form preserves numerical/status
  // failures for the AMPS-facing SEP adapter; the bool wrapper is retained for
  // source compatibility and returns false only for NO_SURFACE.
  swcme::ModelStatus shock_acceleration_state_checked(
      const StepState& S, const double u[3],
      swcme::acceleration::ShockAccelerationState& state) const;
  bool shock_acceleration_state(
      const StepState& S, const double u[3],
      swcme::acceleration::ShockAccelerationState& state) const;

  // Return the point on the observer-anchored analytical Parker field line at
  // a requested heliocentric radius.  The line is generated by rotating the
  // observer radial direction about the configured solar-rotation axis by
  //   Delta phi = -Omega_sun (r-r_obs) / V_sw,
  // exactly matching the Parker field implemented by parker_vec_T_fast().
  // This helper is public so validation and transport adapters can verify the
  // same field-line geometry without maintaining a duplicate tracer.
  bool parker_field_line_point(const StepState& S,
                               const double observer_m[3],
                               double radius_m,
                               double point_m[3]) const;

  // Analytical Parker arc length between two radii on the field line anchored
  // at observer_m.  The colatitude is constant along the Parker line, so the
  // integral is closed-form.  The result is always non-negative.
  double parker_field_line_length(const StepState& S,
                                  const double observer_m[3],
                                  double radius_a_m,
                                  double radius_b_m) const;

  // Intersect an observer's Parker field line with the current shock surface.
  // All geometrical roots between options.inner_radius_m and the observer are
  // returned.  The selected cobpoint is the outermost root, while each root
  // carries the complete production LocalShockState and Parker path length.
  ConnectivityState observer_connectivity(
      const StepState& S, const double observer_m[3],
      const ConnectivityOptions& options = ConnectivityOptions{}) const;

  // Convenience history builder for a stationary observer.  Each requested
  // time is evaluated independently through prepare_step() and the connectivity
  // solver; no hidden temporal hysteresis is introduced.  This makes connection
  // onset/loss and cobpoint migration reproducible and easy to validate.
  std::vector<ConnectivityHistorySample> observer_connectivity_history(
      const std::vector<double>& times_s, const double observer_m[3],
      const ConnectivityOptions& options = ConnectivityOptions{}) const;

  // Backward-compatible scalar diagnostic.  This wrapper is intentionally NOT
  // used by new production internals.  r_eval_m, Rdir_m, and n_hat are retained
  // only so older callers still compile; shock physics is recomputed from the
  // direction through shock_state_direction(), which owns the physical surface
  // location and always samples the upstream plasma there.  Therefore changing
  // an arbitrary field-query radius cannot change rc, Vsh_n, or theta_Bn.
  // New code should use shock_state_direction() directly to obtain the complete
  // upstream/downstream state and explicit surface/shock-existence flags.
  void local_oblique_rc(const StepState& S, const double u[3], const double n_hat[3],
                        double Rdir_m, double r_eval_m,
                        double& rc_out, double& Vsh_n_out, double& thetaBn_out) const;

  // Checked batch evaluators return the first failed sample explicitly and
  // never replace a NaN/Inf/out-of-domain query with zero or ambient flow.
  swcme::ModelStatus evaluate_cartesian_fast_checked(
                               const StepState& S,
                               const double* x_m,const double* y_m,const double* z_m,
                               double* n_m3,double* Vx_ms,double* Vy_ms,double* Vz_ms,
                               std::size_t N) const;

  // Evaluate n and V (radial direction) for arrays of Cartesian points.
  // Legacy source-compatible wrapper: throws std::runtime_error on a failed
  // checked evaluation.
  void evaluate_cartesian_fast(const StepState& S,
                               const double* x_m,const double* y_m,const double* z_m,
                               double* n_m3,double* Vx_ms,double* Vy_ms,double* Vz_ms,
                               std::size_t N) const;

  swcme::ModelStatus evaluate_cartesian_with_B_checked(
                                 const StepState& S,
                                 const double* x_m,const double* y_m,const double* z_m,
                                 double* n_m3,double* Vx_ms,double* Vy_ms,double* Vz_ms,
                                 double* Bx_T,double* By_T,double* Bz_T,
                                 std::size_t N) const;

  // Evaluate n, V, and B (Parker upstream + sheath Bt amplification).
  void evaluate_cartesian_with_B(const StepState& S,
                                 const double* x_m,const double* y_m,const double* z_m,
                                 double* n_m3,double* Vx_ms,double* Vy_ms,double* Vz_ms,
                                 double* Bx_T,double* By_T,double* Bz_T,
                                 std::size_t N) const;

  swcme::ModelStatus evaluate_cartesian_with_B_div_checked(
                                     const StepState& S,
                                     const double* x_m,const double* y_m,const double* z_m,
                                     double* n_m3,double* Vx_ms,double* Vy_ms,double* Vz_ms,
                                     double* Bx_T,double* By_T,double* Bz_T,double* divVsw,
                                     std::size_t N, double dr_frac=1e-3) const;

  // Evaluate n, V, B, and div(V) using the canonical mode-aware divergence:
  // exact 2Vsw/r in SHOCK_ONLY, full Cartesian Jacobian trace in FULL_ICME.
  void evaluate_cartesian_with_B_div(const StepState& S,
                                     const double* x_m,const double* y_m,const double* z_m,
                                     double* n_m3,double* Vx_ms,double* Vy_ms,double* Vz_ms,
                                     double* Bx_T,double* By_T,double* Bz_T,double* divVsw,
                                     std::size_t N, double dr_frac=1e-3) const;

  // Canonical divergence evaluator.  SHOCK_ONLY is analytical (2 V_sw/r).
  // FULL_ICME uses the full Cartesian divergence because its local RH/sheath
  // velocity can contain tangential components and angular gradients.
  swcme::ModelStatus compute_divV_checked(
                           const StepState& S,
                           const double* x_m,const double* y_m,const double* z_m,
                           double* divV,std::size_t N,double dr_frac=1e-3) const;
  void compute_divV(const StepState& S,
                    const double* x_m,const double* y_m,const double* z_m,
                    double* divV,std::size_t N,double dr_frac=1e-3) const;

  // Explicit full-vector second-order Cartesian operator.  This entry point is
  // useful for validation/convergence studies even in SHOCK_ONLY, where the
  // canonical compute_divV() intentionally uses the exact analytical result.
  swcme::ModelStatus compute_divV_cartesian_checked(
                           const StepState& S,
                           const double* x_m,const double* y_m,const double* z_m,
                           double* divV,std::size_t N,double dr_frac=1e-3) const;
  void compute_divV_cartesian(const StepState& S,
                              const double* x_m,const double* y_m,const double* z_m,
                              double* divV,std::size_t N,double dr_frac=1e-3) const;

  // Legacy source-compatible name.  Before Fix 13 this routine always applied
  // the radial formula even to non-radial FULL_ICME flow.  It now delegates to
  // the canonical dimensionality-aware compute_divV_checked() implementation.
  swcme::ModelStatus compute_divV_radial_checked(
                           const StepState& S,
                           const double* x_m,const double* y_m,const double* z_m,
                           double* divV,std::size_t N,double dr_frac=1e-3) const;
  void compute_divV_radial(const StepState& S,
                           const double* x_m,const double* y_m,const double* z_m,
                           double* divV,std::size_t N,double dr_frac=1e-3) const;

  // Quick diagnostic at a direction u (unit). Returns false when a finite
  // geometry (currently SSE) has no shock surface in that direction.  On a
  // false return Rdir=0, n_hat=(0,0,0), rc=1, and Vsh_n=0.
  bool diagnose_direction(const StepState& S, const double u[3],
                          double& Rdir_m,double n_hat[3],
                          double& rc_loc,double& Vsh_n) const;

  // Build a topologically unique triangular shock surface. `nTheta` is the
  // number of polar intervals and `nPhi` the number of unique vertices per
  // non-polar ring.  Periodicity is implemented by wrapped connectivity; no
  // duplicate phi=2pi seam nodes are stored.
  ShockMesh build_shock_mesh(const StepState& S, std::size_t nTheta, std::size_t nPhi) const;

  // Compute per-triangle metrics and reject repeated-index or numerically
  // degenerate cells instead of returning a zero-area/zero-normal fallback.
  void compute_triangle_metrics(const ShockMesh& M, TriMetrics& T) const;

  // Build/select from the canonical area-weighted cell distribution.  This is
  // the production primitive for spatially uniform per-unit-area shock source
  // sampling; choosing triangles uniformly by index is intentionally not
  // provided because it biases a nonuniform angular mesh.
  AreaSamplingTable build_area_sampling_table(const TriMetrics& T) const;
  std::size_t sample_triangle_by_area(const AreaSamplingTable& table,
                                      double unit_uniform) const;

  // Useful default volume box (apex-aligned, shifted outward).  OUT04 rejects
  // non-finite/negative half_AU, N<2, or a generated box whose SI bounds or
  // spans overflow, so the factory never returns a structurally invalid box.
  BoxSpec default_apex_box(const StepState& S,double half_AU,int N) const;

  // --- Tecplot writers (VARIABLES defined exactly below) --------------------
  //
  // VARIABLES (same unit-qualified tokens and order in every 3-D product):
  //   1: X[m],  2: Y[m],  3: Z[m],
  //   4: n[m^-3],  5: Vx[m/s], 6: Vy[m/s], 7: Vz[m/s],
  //   8: Bx[T], 9: By[T], 10: Bz[T], 11: divVsw[s^-1],
  //   12: rc[-], 13: Vsh_n[m/s],
  //   14: nx[-], 15: ny[-], 16: nz[-],              // cell geometric normal
  //   17: area[m^2], 18: rc_mean[-], 19: Vsh_n_mean[m/s],
  //   20: tnx[-], 21: tny[-], 22: tnz[-],           // reserved (e.g., tension dir)
  //   23: cx[m], 24: cy[m], 25: cz[m]               // cell centroid
  //
  // Notes:
  //  • In "surface_cells" (FETRIANGLE, BLOCK), [1–3,12–13] are NODAL, all others
  //    are CELLCENTERED. Nodal arrays are length Nv; cell arrays length Ne.
  //  • In "surface_nodal" (FEPOINT), every line lists all 25 values per node;
  //    rc and Vsh_n are meaningful; nx,ny,nz are filled with nodal normals.
  //  • In "volume_box" (structured POINT), n,V,B,divV are meaningful; surface-
  //    specific quantities (rc, Vsh_n, normals, area, centroids) are zeros.
  //  • Non-finite physics data are rejected. Checked writer APIs return the
  //    failure status; writers never replace a bad value by zero.
  //  • Checked writers distinguish FILE_OPEN_FAILURE from FILE_WRITE_FAILURE
  //    and inspect write, flush, stream-error, and close results.  Their final
  //    optional FileOperations pointer is a deterministic validation seam;
  //    normal callers omit it to select the production stdio backend.
  //  • OUT03 stages each validated product beside its destination and commits
  //    it by atomic rename only after close succeeds. FILE_COMMIT_FAILURE
  //    leaves an existing regular destination unchanged.
  //  • OUT05 scans every surface node and every structured volume/face point
  //    before FileOperations is selected.  A non-finite coordinate or radius
  //    below solarwind::MIN_RADIUS_M therefore returns its model-domain status
  //    and flattened row index without creating a staging file.
  //  • OUT04 first validates the BoxSpec itself: finite centers/extents,
  //    nonnegative half sizes, Ni/Nj/Nk >= 2, representable bounds/spans, and
  //    overflow-safe total cardinality.  Structural rejection precedes OUT05
  //    and has the same no-open/no-staging guarantee.
  //  • OUT06 validates ShockMesh structure, nodal invariants, one-based
  //    connectivity, cell quality, and TriMetrics consistency before OUT05 or
  //    any output access.  A completely empty TriMetrics requests canonical
  //    computation; partial or stale records return INVALID_MESH.
  //  • OUT01 independently parses all three product shapes and freezes titles,
  //    unit-qualified variable tokens/order, zone declarations, row widths,
  //    counts, connectivity, finite data, and exact end-of-file consumption.
  // --------------------------------------------------------------------------
  swcme::ModelStatus write_tecplot_dataset_bundle_checked(
                                    const ShockMesh& M,const TriMetrics& T,
                                    const StepState& S,const BoxSpec& B,
                                    const char* path,
                                    const swcme::output::FileOperations*
                                        file_operations=nullptr) const;
  bool write_tecplot_dataset_bundle(const ShockMesh& M,const TriMetrics& T,
                                    const StepState& S,const BoxSpec& B,
                                    const char* path) const;

  // Standalone 2-D face (min-X plane of a given BoxSpec) in Tecplot.
  swcme::ModelStatus write_box_face_minX_tecplot_structured_checked(
                                              const StepState& S,
                                              const BoxSpec& B,
                                              const char* path,
                                              const swcme::output::FileOperations*
                                                  file_operations=nullptr) const;
  bool write_box_face_minX_tecplot_structured(const StepState& S,
                                              const BoxSpec& B,
                                              const char* path) const;

  // Surface-only writer (cell metrics + nodal rc/Vsh_n) as a single zone.
  swcme::ModelStatus write_shock_surface_center_metrics_tecplot_checked(
                                                  const ShockMesh& M,
                                                  const TriMetrics& T,
                                                  const char* path,
                                                  const swcme::output::FileOperations*
                                                      file_operations=nullptr) const;
  bool write_shock_surface_center_metrics_tecplot(const ShockMesh& M,
                                                  const TriMetrics& T,
                                                  const char* path) const;

private:
  // Copy assignment is the only legacy operation capable of replacing a 3-D
  // model's private Params.  Reject it after preparation with the same PST01
  // lifecycle rule used by the 1-D fluent setters.
  void require_configuration_mutable(const char* operation) const {
    if (configuration_locked_.load(std::memory_order_acquire)) {
      throw std::logic_error(std::string("swcme3d::")+operation+
          ": model configuration is immutable after successful prepare_step(); "
          "construct model.reconfigured(params) instead");
    }
  }

  Params P_;
  // Runtime owner token stamped into each StepState prepared by this object.
  // It is deliberately unrelated to configuration values so independently
  // constructed but numerically identical models remain distinct owners.
  swcme::ModelIdentity model_identity_;
  // The atomic lifecycle flag permits concurrent read-only queries after a
  // successful prepare without exposing any path back to mutable Params.
  mutable std::atomic<bool> configuration_locked_;
};

} // namespace swcme3d
