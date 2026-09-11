// ============================================================================
// swcme3d.cpp
// ----------------------------------------------------------------------------
// SOLAR WIND + CME FORWARD-SHOCK SEMI-ANALYTICAL MODEL
//
// PHYSICS OVERVIEW
// ================
// 1) Upstream density n_up(r): Leblanc–Dulk–Bougeret (1998)
//    n(r) [cm^-3] ≈ 3.3e5 (r/Rs)^(-2) + 4.1e6 (r/Rs)^(-4) + 8.0e7 (r/Rs)^(-6);
//    We rescale these coefficients to match user n(1 AU) and convert to SI.
//    Implementation:
//      n(r) [m^-3] = C2 * r^-2 + C4 * r^-4 + C6 * r^-6,
//    where constants C2,C4,C6 (SI) are cached in StepState.
//
// 2) Upstream Parker spiral B_up(r,θ):
//    Br ∝ r^-2 and Bφ = -Br (Ω r sinθ / Vsw).  In 3-D, θ must be the
//    LOCAL colatitude measured from the solar-rotation axis, not one global
//    parameter.  The local azimuthal direction is
//       e_phi = (Omega_hat × e_r)/|Omega_hat × e_r|,
//    with Bphi -> 0 smoothly on the rotation axis.  StepState therefore caches
//    the normalized solar axis and the equatorial pitch coefficient Ω AU/Vsw.
//    B1AU_nT remains a total-field normalization at the legacy/reference
//    colatitude Params::sin_theta; that parameter no longer controls local
//    winding anywhere in the 3-D domain.
//
// 3) CME apex kinematics: Drag-Based Model (DBM)
//    Let u = V_sh - V_sw be the excess speed. DBM gives
//      ΔV(t)=ΔV0/(1+Γ|ΔV0|t), r(t)=r0+V_sw t+sgn(ΔV0)ln(1+Γ|ΔV0|t)/Γ,
//    where Γ is the drag parameter (converted to SI).
//
// 4) Shock geometry (radius and normal along direction u):
//    • Sphere: Rdir = r_sh, n = u.
//    • Ellipsoid: implicit x^2/a^2 + y^2/b^2 + z^2/c^2 = 1 with a=r_sh;
//      for a ray x=λ u1, y=λ u2, z=λ u3 (in apex frame) => λ = 1 / sqrt(u1^2/a^2+…).
//      The outward normal ∝ (x/a^2, y/b^2, z/c^2); rotate to global.
//    • SSE: finite self-similar-expansion spherical cap.  The cap is generated
//      by a sphere centered on the CME axis and tangent to rays at the configured
//      half width.  Radius and outward normal are analytic, and no surface is
//      returned outside the finite angular extent.
//
// 5) Region structure and smoothing:
//    SHOCK_ONLY leaves the analytical Parker/Leblanc transport background
//    unchanged and is paired with explicit SOURCE acceleration. FULL_ICME is
//    paired with RESOLVED_COMPRESSION: one finite C1 layer maps local upstream
//    plasma to the exact surface-owned RH downstream state, after which the
//    phenomenological sheath, magnetic ejecta, and post-ICME ambient follow.
//    Leading/trailing surfaces are local self-similar fractions of R_sh(direction)
//    and use the same C1 convention, so 1-D/3-D numerical acceleration cannot
//    differ merely because their shock smoothing differs.  The canonical shock
//    diagnostic itself remains the exact mathematical RH discontinuity.
//
// 6) Oblique ideal-MHD shock state (local):
//    Compute θBn, the oblique fast-mode speed, and the normal relative inflow.
//    A geometric front is a physical shock only for M_fast>1.  The shared
//    Rankine-Hugoniot solver then returns compression and the complete
//    conservative downstream rho, p, V, and B state.
//
// 7) Magnetic field jump:
//    The immediate downstream magnetic field is the RH solution: B_n remains
//    continuous while B_t follows the full tangential jump conditions.  The
//    interior sheath then relaxes phenomenologically toward the Parker field.
//
// 8) Divergence of V:
//    SHOCK_ONLY has the exact constant-speed radial result ∇·V=2 V_sw/r.
//    FULL_ICME is generally non-radial after an oblique RH jump and also varies
//    angularly with the finite shock surface, so it uses the full Cartesian
//    Jacobian trace dVx/dx+dVy/dy+dVz/dz with a verified second-order stencil.
//
// ─────────────────────────────────────────────────────────────────────────────
// SOLAR WIND + CME FORWARD-SHOCK SEMI-ANALYTICAL MODEL
// ─────────────────────────────────────────────────────────────────────────────
//
// What this file provides
// -----------------------
// • A fast 3D kinematic/phenomenological model of a CME-driven forward shock
//   propagating into a Parker-spiral solar wind. The model returns:
//     - plasma number density n [m^-3]
//     - bulk velocity V = (Vx,Vy,Vz) [m/s]
//     - magnetic field B = (Bx,By,Bz) [Tesla], upstream is Parker; tangential
//       field immediately behind the shock is the conservative RH downstream
//       state and relaxes phenomenologically through the sheath.
//     - ∇·V [1/s]: analytical in SHOCK_ONLY; full Cartesian second-order
//       divergence in FULL_ICME.
//     - a topologically unique triangular shock-surface mesh with nodal normals,
//       nodal rc, nodal normal shock speed, per-cell metrics, and a canonical
//       physical-area CDF for stochastic source-patch selection.
// • Tecplot dataset writers for the shock surface and for a structured volume
//   box near the apex (plus a 2-D face zone). Non-finite physics values are
//   rejected explicitly; output code does not sanitize them into finite data.
//
// Headline physics approximations
// -------------------------------
// A) Ambient plasma density n_up(r):
//    Uses Leblanc, Dulk & Bougeret (1998) density law
//       n(r) [cm^-3] = A (Rs/r)^2 + B (Rs/r)^4 + C (Rs/r)^6
//    with A=3.3e5, B=4.1e6, C=8.0e7 and Rs the solar radius. We scale this law
//    so that n(1 AU) matches user parameter n1AU_cm3. Output is SI [m^-3].
//
// B) Magnetic field B_up(r,θ):
//    Parker spiral (Parker 1958, ApJ 128:664).
//    In spherical coordinates (r,θ,φ) with θ measured from the configured
//    solar-rotation axis, and assuming azimuthal symmetry of the solar
//    rotation, we use a radial Br and toroidal Bφ component:
//       q(r,θ) = Ω⊙ r sinθ / V_sw,  (dimensionless local pitch)
//       Br(r)   = Br(1 AU) (1 AU / r)^2
//       Bφ(r)   = -Br(r) * q(r,θ)
//    A desired |B|(1 AU) = B1AU_nT is enforced at a documented reference
//    colatitude (Params::sin_theta) by choosing Br(1 AU) accordingly.  At every
//    evaluated point, however, the local sin(theta) and e_phi are computed
//    geometrically from the explicit solar-rotation axis.  This distinction
//    preserves the existing normalization input while making the 3-D vector
//    field physically consistent away from the reference latitude.
//
// C) CME apex kinematics:
//    Shared swcme::kinematics supports BALLISTIC, DBM and DATA_DRIVEN modes.
//    The sign-aware DBM solves
//       d(DeltaV)/dt = -Gamma DeltaV |DeltaV|
//    with a=|DeltaV0|:
//       DeltaV(t)=DeltaV0/(1+Gamma a t)
//       r(t)=r0+Vsw t+sgn(DeltaV0) log(1+Gamma a t)/Gamma.
//    Gamma=0 uses the exact ballistic limit; very small Gamma uses a stable
//    series for the logarithmic distance.  DATA_DRIVEN mode uses monotone
//    PCHIP height-time interpolation and its derivative as the apex speed.

// D) Shock shape & direction-dependent speed:
//    Shape options: Sphere, self-similar Ellipsoid, and a finite true-SSE
//    spherical cap.  Every shape scales linearly with the apex distance, so a
//    surface point on a fixed heliocentric ray has radial speed
//       dR/dt = V_apex * R/R_apex.
//    The physically relevant shock-normal speed is the normal projection
//       V_sh,n = dR/dt * (e_r · n_hat).
//    This replaces the old ad hoc cosine flank-speed factor and also fixes the
//    ellipsoid, whose flanks previously moved at the full apex speed.
//
// E) Local ideal-MHD fast-shock solution:
//    The model computes the upstream normal relative inflow and oblique
//    fast-mode speed.  M_fast<=1 gives an explicit no-shock state with rc=1.
//    For M_fast>1, swcme_shock.hpp solves the ideal-MHD Rankine-Hugoniot
//    conditions and returns the complete downstream primitive state.
//
// F) Sheath / ejecta blending (C^1 smoothing):
//    Three radial transitions along a given direction u:
//      1) Exact upstream/downstream RH discontinuity at the shock.
//      2) Sheath → Ejecta ramp at the leading edge r_le = r_sh - dr_sheath.
//      3) Ejecta → Downstream ambient at the trailing edge r_te = r_le - dr_me.
//    Each edge uses a smoothstep s(x)=x^2(3-2x) with its own width. Inside the
//    sheath, density/velocity/B relax from the exact RH downstream state
//    toward phenomenological leading-edge targets with a user power p.
//    Ejecta density is a fraction f_ME of n_up (simple cavity or enhancement).
//
// G) Divergence of V:
//    In SHOCK_ONLY the production value is exact: ∇·V=2 V_sw/r.  In FULL_ICME
//    the velocity may have tangential RH components and angular gradients, so
//    the production operator evaluates ∂Vx/∂x+∂Vy/∂y+∂Vz/∂z with centered
//    second-order Cartesian differences (and an explicit second-order one-sided
//    stencil only when the inner model boundary prevents centering).
//
// H) Triangulated shock surface & per-cell metrics:
//    The surface is parameterized by (θ,φ) in the apex-aligned frame. We
//    compute nodal positions, nodal unit normals, nodal rc and Vsh_n. Then,
//    per triangle we compute geometric normal, area, centroid, and the mean
//    of (rc, Vsh_n) across the three vertices.

// I) Output (Tecplot)
//    We write a dataset with VARIABLES (common order across all zones):
//      1  X [m], 2 Y [m], 3 Z [m],
//      4  n [m^-3], 5 Vx [m/s], 6 Vy [m/s], 7 Vz [m/s],
//      8  Bx [T], 9 By [T], 10 Bz [T], 11 divVsw [1/s],
//      12 rc [-], 13 Vsh_n [m/s],
//      14 nx [-], 15 ny [-], 16 nz [-],
//      17 area [m^2], 18 rc_mean [-], 19 Vsh_n_mean [m/s],
//      20 tnx [-], 21 tny [-], 22 tnz [-]  (reserved),
//      23 cx [m], 24 cy [m], 25 cz [m]     (triangle centroids)
//
// Tecplot VARIABLES (global order across all zones)
// -------------------------------------------------
//  1  "X"       [m]   : position
//  2  "Y"       [m]
//  3  "Z"       [m]
//  4  "n"     [m^-3]  : density
//  5  "Vx"    [m/s]   : bulk velocity
//  6  "Vy"    [m/s]
//  7  "Vz"    [m/s]
//  8  "Bx"      [T]   : magnetic field (Parker upstream; Bt amplified in sheath)
//  9  "By"      [T]
// 10  "Bz"      [T]
// 11  "divVsw" [1/s]  : divergence of V (analytic or Cartesian)
// 12  "rc"      [-]   : compression ratio (nodal on surface zones)
// 13  "Vsh_n"  [m/s]  : normal shock speed (nodal on surface zones)
// 14  "nx"      [-]   : triangle/geometric normal (cell-centered in surface_cells)
// 15  "ny"      [-]
// 16  "nz"      [-]
// 17  "area"   [m^2]  : triangle area (cell-centered in surface_cells)
// 18  "rc_mean"[-]    : mean nodal rc over each triangle (cell-centered)
// 19  "Vsh_n_mean"[m/s]: mean nodal Vsh_n over each triangle (cell-centered)
// 20  "tnx"     [-]   : reserved (0)
// 21  "tny"     [-]
// 22  "tnz"     [-]
// 23  "cx"      [m]   : triangle centroid X (cell-centered)
// 24  "cy"      [m]
// 25  "cz"      [m]

//    Zones:
//      • surface_cells  : FETRIANGLE/BLOCK, with cell-centered metrics
//                         (area, rc_mean, Vsh_n_mean, nx,ny,nz, cx,cy,cz)
//                         and nodal rc, Vsh_n at vertices. This zone appears
//                         first so Tecplot shows non-zero cell metrics by default.

//      • volume_box     : Structured POINT, 3D grid around shock apex.
//      • box_face_minX  : Structured POINT (2D), the plane x = minX of the box.
//    Physics evaluators propagate ModelStatus on invalid/non-finite states.
//    No numeric output is repaired after the fact.
//
// J) References (short list)
//    • Parker, E.N. (1958), ApJ 128, 664 — Parker spiral magnetic field.
//    • Leblanc, Y.; Dulk, G.A.; Bougeret, J.-L. (1998), Sol. Phys. 183, 165 —
//      empirical coronal/solar-wind density profile.
//    • Vršnak, B. et al. (2013), Sol. Phys. 285, 295 — Drag-based CME propagation.
//    • Edmiston, J.P.; Kennel, C.F. (1984), J. Plasma Phys. 32, 429 — Oblique
//      MHD shock jump conditions (used here as a proxy).
//    • Priest, E. (2014), “Magnetohydrodynamics of the Sun”, CUP — background.
//    • Russell, C.T.; Mulligan, T. (2002), Planet. Space Sci. 50, 527 — CME
//      sheath/ICME structure and in situ signatures.
//    • Manchester, W.B. IV et al. (2005), ApJ 622, 1225 — CME-driven shocks.
//
// Build
// -----
//   g++ -std=c++17 -O3 -march=native demo3d_2.cpp swcme3d.cpp -o demo
//
// Header pairing
// --------------
// This file matches the earlier swcme3d.hpp you’re using. It defines AU, Rs,
// PI in the swcme3d namespace (giving them external linkage so other TUs like
// demo3d_2.cpp can reference them), and implements all methods declared there.
//
// ─────────────────────────────────────────────────────────────────────────────


//
// I/O & VARIABLES (Tecplot):
// --------------------------
// VARIABLES (fixed order across all writers):
//   1:X [m], 2:Y [m], 3:Z [m],
//   4:n [m^-3], 5:Vx [m/s], 6:Vy [m/s], 7:Vz [m/s],
//   8:Bx [T], 9:By [T], 10:Bz [T], 11:divVsw [1/s],
//   12:rc [-], 13:Vsh_n [m/s],
//   14:nx [-], 15:ny [-], 16:nz [-],         // cell geometric normal
//   17:area [m^2], 18:rc_mean [-], 19:Vsh_n_mean [m/s],
//   20:tnx [-], 21:tny [-], 22:tnz [-],      // reserved (unused)
//   23:cx [m], 24:cy [m], 25:cz [m]          // cell centroid
//
// SANITIZATION:
//  • Tecplot writers reject non-finite physics data and checked writer APIs
//    return ModelStatus; they do not substitute zeros/ambient values.
//
// PERFORMANCE NOTES:
//  • StepState caches Leblanc coefficients, Parker constants, and geometry.
//  • No std::pow in hot loops; replaces with multiplies on r^-2, etc.
//  • Smooth edges use cached inv(2w). Vectorization hint via GCC ivdep.
//  • No heap allocations in evaluators; all work is per-point arithmetic.
//
// ----------------------------------------------------------------------------

#include "swcme3d.hpp"

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <vector>
#include <string>
#include <algorithm>
#include <limits>
#include <stdexcept>

// Vectorization hint (safe: no loop-carried deps)
#if defined(__GNUC__)
  #define SWCME_IVDEP _Pragma("GCC ivdep")
#else
  #define SWCME_IVDEP
#endif

namespace swcme3d {
  // Physical constants (exported)
  const double AU = swcme::constants::AU_M;
  const double Rs = swcme::constants::SOLAR_RADIUS_M;
  const double PI = swcme::constants::PI;

  // Small helper (in namespace to avoid header pollution)
  inline double clamp01(double v){ return (v<0.0)?0.0:(v>1.0?1.0:v); }
}

// --- file-local helpers ------------------------------------------------------
// Normalize only when the vector is finite and has a meaningful magnitude.
// The former helper silently replaced a zero/invalid vector by +X, which could
// turn a numerical failure into an apparently valid CME/shock direction.
static inline bool normalize_checked(double v[3]){
  const double m=std::hypot(v[0],std::hypot(v[1],v[2]));
  const double tol=64.0*std::numeric_limits<double>::min();
  if (!std::isfinite(m) || !(m>tol)) return false;
  v[0]/=m; v[1]/=m; v[2]/=m;
  return std::isfinite(v[0]) && std::isfinite(v[1]) && std::isfinite(v[2]);
}

static inline bool finite3(const double v[3]){
  return std::isfinite(v[0]) && std::isfinite(v[1]) && std::isfinite(v[2]);
}
static constexpr double MU0 = swcme::constants::VACUUM_PERMEABILITY_N_A2; // [H/m]
static constexpr double MP  = swcme::constants::PROTON_MASS_KG;           // [kg]
static constexpr double KB  = swcme::constants::BOLTZMANN_J_K;            // [J/K]
static constexpr double OMEGA_SUN = swcme::constants::SOLAR_ROTATION_RAD_S; // legacy/test-visible model default [rad/s]

static inline double smoothstep01(double x){
  if (x<=0) return 0;
  if (x>=1) return 1;
  return x*x*(3-2*x);
}

// Rotate a vector about a UNIT axis using Rodrigues' formula.  Connectivity
// uses this operation to construct the analytical Parker field line in a form
// that works for an arbitrary solar-rotation axis; no global +Z assumption is
// hidden in the field-line tracer.
static inline void rotate_about_unit_axis(const double v[3], const double axis[3],
                                          double angle, double out[3]) {
  const double c=std::cos(angle), s=std::sin(angle);
  const double axv[3]={
      axis[1]*v[2]-axis[2]*v[1],
      axis[2]*v[0]-axis[0]*v[2],
      axis[0]*v[1]-axis[1]*v[0]};
  const double adv=axis[0]*v[0]+axis[1]*v[1]+axis[2]*v[2];
  const double one_minus_c=1.0-c;
  out[0]=c*v[0]+s*axv[0]+one_minus_c*adv*axis[0];
  out[1]=c*v[1]+s*axv[1]+one_minus_c*adv*axis[1];
  out[2]=c*v[2]+s*axv[2]+one_minus_c*adv*axis[2];
}

static inline double norm3(const double v[3]) {
  return std::hypot(v[0],std::hypot(v[1],v[2]));
}

// Evaluate the upstream Parker field using the local 3-D spherical basis.
//
// IMPORTANT PHYSICS NOTE:
// The Parker azimuthal unit vector is parallel to Omega_hat × e_r.  The old
// implementation formed (Omega_hat × e_r) × e_r, which points in the local
// meridional direction and therefore rotated the spiral field into the wrong
// plane.  In addition, the old code multiplied the field everywhere by one
// globally supplied sin(theta).  A true 3-D Parker field instead uses the
// local colatitude at each point:
//
//   sin(theta_local) = |Omega_hat × e_r|,
//   B_r               = B_r(1 AU) / r_AU^2,
//   B_phi             = -B_r * (Omega*AU/V_sw) * r_AU * sin(theta_local).
//
// At the rotation poles sin(theta_local)=0, so B_phi must vanish.  We avoid
// normalizing the undefined azimuthal basis there and return the purely radial
// limit analytically.  This keeps the field finite and direction-independent
// as the pole is approached.
static inline void parker_vec_T_fast(const swcme3d::StepState& S,
                                     const double u[3], double r_m,
                                     double B_out[3]){
  // The complete Parker scalar/vector physics now lives in the common solar-
  // wind core.  This file-local wrapper preserves the existing hot-path call
  // sites while making it impossible for the 3-D formula to drift away from
  // the 1-D Parker components.  Geometry supplies only the normalized solar
  // axis and local radial direction.
  const std::array<double,3> axis={{S.solar_axis_hat[0],
                                    S.solar_axis_hat[1],
                                    S.solar_axis_hat[2]}};
  const std::array<double,3> radial={{u[0],u[1],u[2]}};
  const std::array<double,3> B=swcme::solarwind::parker_field_cartesian(
      S.common.solar_wind,axis,radial,r_m);
  B_out[0]=B[0];
  B_out[1]=B[1];
  B_out[2]=B[2];
}

// ----------------------------------------------------------------------------
// Model implementation
// ----------------------------------------------------------------------------
namespace swcme3d {

Model::Model(const Params& P)
    : P_(P), model_identity_(swcme::next_model_identity()),
      configuration_locked_(false) {}

Model::Model(const Model& other)
    : P_(other.P_), model_identity_(swcme::next_model_identity()),
      configuration_locked_(false) {}

Model& Model::operator=(const Model& other) {
  if (this!=&other) {
    require_configuration_mutable("operator=");
    P_=other.P_;
    // Only an unprepared receiver can reach this branch.  Rotate its identity
    // anyway because assignment replaces the logical model represented by the
    // object and PST02 ownership must not depend on configuration equality.
    model_identity_=swcme::next_model_identity();
    configuration_locked_.store(false,std::memory_order_release);
  }
  return *this;
}

swcme::defaults::ObserverScopeStatus Model::observer_scope_status(
    const StepState& S, const double observer_m[3]) const {
  // This legacy value-returning API has no status channel.  Reject a foreign
  // state by exception before reading geometry or modifying caller-visible
  // state, matching the behavior of other source-compatible wrappers.
  swcme::throw_if_error(validate_prepared_state(
      S,"swcme3d::observer_scope_status"));
  // Keep scope bookkeeping side-effect free and geometry-consistent: the
  // observer direction is tested against the same production surface routine
  // used by connectivity and shock diagnostics.  No angular-width shortcut or
  // apex-radius proxy is used here.
  if (observer_m == nullptr || !finite3(observer_m)) {
    return swcme::defaults::observer_scope_status(
        P_.region_mode, P_.shock_acceleration_mode, false, 0.0,
        std::numeric_limits<double>::quiet_NaN());
  }

  const double r_obs = norm3(observer_m);
  if (!std::isfinite(r_obs) || !(r_obs > 0.0)) {
    return swcme::defaults::observer_scope_status(
        P_.region_mode, P_.shock_acceleration_mode, false, 0.0, r_obs);
  }

  const double inv_r = 1.0 / r_obs;
  const double u[3] = {observer_m[0] * inv_r, observer_m[1] * inv_r,
                       observer_m[2] * inv_r};
  double Rdir = 0.0;
  double n_hat[3] = {0.0, 0.0, 0.0};
  const bool exists = shape_radius_normal(S, u[0], u[1], u[2],
                                          Rdir, n_hat);
  return swcme::defaults::observer_scope_status(
      P_.region_mode, P_.shock_acceleration_mode, exists, Rdir, r_obs);
}

StepState Model::prepare_step(double t_s) const {
  // Reject invalid configuration once, before any vector normalization, unit
  // conversion, or finite-value fallback can transform the caller's input.
  // This replaces the former collection of local guards/clamps with one
  // auditable contract shared with the 1-D model.
  const swcme::config::ValidationResult validation=validate();
  if (!validation.ok()) {
    throw std::invalid_argument(validation.summary("swcme3d"));
  }
  if (!std::isfinite(t_s) || t_s<0.0) {
    throw std::invalid_argument("swcme3d: time must be finite and >= 0");
  }

  StepState S{};
  // Stamp the owner before constructing cached geometry and physics.  A
  // successful return therefore always carries the identity of this exact
  // model; an exception exposes no partially prepared state to the caller.
  S.owner_model_identity=model_identity_;
  // Record the complete, stable configuration snapshot after validation.
  // State-consuming APIs recompute this digest before touching their outputs,
  // which makes cross-configuration rejection transactional.
  S.configuration_digest=configuration_digest(P_);
  S.time_s=t_s;

  // 1) Apex-aligned orthonormal basis (e1 along CME apex direction).  The CME
  // vector is guaranteed non-zero by centralized validation, so normalization
  // cannot silently substitute the legacy +X fallback here.
  double e1[3]={P_.cme_dir[0],P_.cme_dir[1],P_.cme_dir[2]};
  if (!::normalize_checked(e1))
    throw std::runtime_error("swcme3d: validated CME direction became degenerate");
  double tmp[3]={0,0,1}; if (std::fabs(e1[2])>0.9){ tmp[0]=1; tmp[1]=0; tmp[2]=0; }
  double e2[3]={ e1[1]*tmp[2]-e1[2]*tmp[1],
                 e1[2]*tmp[0]-e1[0]*tmp[2],
                 e1[0]*tmp[1]-e1[1]*tmp[0] };
  if (!::normalize_checked(e2))
    throw std::runtime_error("swcme3d: failed to construct CME transverse basis");
  double e3[3]={ e1[1]*e2[2]-e1[2]*e2[1],
                 e1[2]*e2[0]-e1[0]*e2[2],
                 e1[0]*e2[1]-e1[1]*e2[0] };
  if (!::normalize_checked(e3))
    throw std::runtime_error("swcme3d: failed to construct CME orthogonal basis");
  S.e1[0]=e1[0]; S.e1[1]=e1[1]; S.e1[2]=e1[2];
  S.e2[0]=e2[0]; S.e2[1]=e2[1]; S.e2[2]=e2[2];
  S.e3[0]=e3[0]; S.e3[1]=e3[1]; S.e3[2]=e3[2];

  // Normalize and cache the solar-rotation axis independently of the CME
  // propagation frame.  Centralized configuration validation already proved
  // the vector finite/non-zero and the rotation rate finite/non-negative, so
  // this block performs geometry only and contains no hidden repair policy.
  {
    const double ax=P_.solar_rotation_axis[0];
    const double ay=P_.solar_rotation_axis[1];
    const double az=P_.solar_rotation_axis[2];
    const double axis_norm=std::hypot(ax,std::hypot(ay,az));
    S.solar_axis_hat[0]=ax/axis_norm;
    S.solar_axis_hat[1]=ay/axis_norm;
    S.solar_axis_hat[2]=az/axis_norm;
    S.solar_rotation_rate_rad_s=P_.solar_rotation_rate_rad_s;
  }

  // 2) Dimensionality-independent ambient state and apex kinematics.
  //
  // swcme::core::prepare() is now the single production path that converts
  // public heliophysics units, normalizes the Leblanc/Parker background, and
  // evaluates BALLISTIC/DBM/DATA_DRIVEN apex kinematics.  The 3-D wrapper adds
  // only its genuinely geometric information (rotation axis, CME frame, shock
  // surface).  This removes the last independent 3-D copies of the density
  // normalization and kinematic setup.
  swcme::core::CommonConfig common_cfg;
  common_cfg.V_sw_kms=P_.V_sw_kms;
  common_cfg.n1AU_cm3=P_.n1AU_cm3;
  common_cfg.B1AU_nT=P_.B1AU_nT;
  common_cfg.T_K=P_.T_K;
  common_cfg.gamma_ad=P_.gamma_ad;
  common_cfg.parker_reference_sin_theta=P_.sin_theta;
  common_cfg.solar_rotation_rate_rad_s=P_.solar_rotation_rate_rad_s;
  common_cfg.kinematics_mode=P_.kinematics_mode;
  common_cfg.r0_Rs=P_.r0_Rs;
  common_cfg.V0_sh_kms=P_.V0_sh_kms;
  common_cfg.Gamma_kmInv=P_.Gamma_kmInv;
  common_cfg.data_time_s=P_.data_time_s;
  common_cfg.data_radius_Rs=P_.data_radius_Rs;
  common_cfg.data_extrapolation=P_.data_extrapolation;

  S.common=swcme::core::prepare(common_cfg,t_s);
  if (S.common.apex.status!=swcme::kinematics::Status::Ok) {
    throw std::runtime_error(std::string("swcme3d kinematics: ")+
                             swcme::kinematics::status_name(S.common.apex.status));
  }
  if (!std::isfinite(S.common.apex.radius_m) ||
      S.common.apex.radius_m < swcme::solarwind::MIN_RADIUS_M) {
    throw std::runtime_error(
        swcme::ModelStatus::make_value(
            swcme::StatusCode::OutsideModelDomain,
            "swcme3d::prepare_step shock radius",S.common.apex.radius_m)
            .summary());
  }

  // Populate the legacy/public StepState mirrors from the common state.  No
  // equations are repeated here; these assignments exist only to preserve the
  // current public ABI/source expectations of callers and validation tools.
  S.V_sw_ms=S.common.solar_wind.V_sw_m_s;
  S.kinematics_mode=P_.kinematics_mode;
  S.r_sh_m=S.common.apex.radius_m;
  S.V_sh_ms=S.common.apex.speed_m_s;
  S.a_m=S.r_sh_m;
  S.C2=S.common.solar_wind.C2;
  S.C4=S.common.solar_wind.C4;
  S.C6=S.common.solar_wind.C6;
  S.k_AU=S.common.solar_wind.k_AU_equatorial;
  S.Br1AU_T=S.common.solar_wind.Br1AU_T;

  // 3) Common self-similar FULL_ICME region geometry.  Public thickness and
  // smoothing parameters are AU at a 1-AU shock and therefore dimensionless
  // fractions when applied to a local shock radius.  apex_regions is retained
  // for output/backward StepState fields; the field evaluators call the same
  // make_boundaries() routine with each *local* Rdir so finite-SSE/ellipsoid
  // flanks never subtract an apex-sized absolute thickness.
  S.acceleration_config.mode=P_.shock_acceleration_mode;
  S.acceleration_config.relative_source_weight_per_area=
      P_.relative_source_weight_per_area;
  S.region_config.mode=P_.region_mode;
  S.region_config.shock_smooth_fraction=
      (P_.shock_acceleration_mode==swcme::acceleration::Mode::ResolvedCompression)
          ? P_.edge_smooth_shock_AU_at1AU : 0.0;
  S.region_config.sheath_fraction=P_.sheath_thick_AU_at1AU;
  S.region_config.ejecta_fraction=P_.ejecta_thick_AU_at1AU;
  S.region_config.leading_smooth_fraction=P_.edge_smooth_le_AU_at1AU;
  S.region_config.trailing_smooth_fraction=P_.edge_smooth_te_AU_at1AU;
  S.region_config.sheath_ramp_power=P_.sheath_ramp_power;
  S.region_config.V_sheath_LE_factor=P_.V_sheath_LE_factor;
  S.region_config.f_ME=P_.f_ME;
  S.region_config.V_ME_factor=P_.V_ME_factor;
  S.apex_regions=swcme::regions::make_boundaries(S.r_sh_m,S.region_config);

  // Legacy mirrors remain populated for existing output/API users.  The
  // shock width is now meaningful only in RESOLVED_COMPRESSION mode and is
  // generated by the same common boundary helper used by 1-D and every 3-D ray.
  S.dr_sheath_m=S.apex_regions.sheath_thickness_m;
  S.dr_me_m=S.apex_regions.ejecta_thickness_m;
  S.w_shock_m=S.apex_regions.smooth_shock_width_m;
  S.w_le_m=S.apex_regions.smooth_le_width_m;
  S.w_te_m=S.apex_regions.smooth_te_width_m;
  S.r_le_m=S.apex_regions.R_le_m;
  S.r_te_m=S.apex_regions.R_te_m;
  S.inv2w_sh=(S.w_shock_m>0.0)?0.5/S.w_shock_m:0.0;
  S.inv2w_le=(S.w_le_m>0.0)?0.5/S.w_le_m:0.0;
  S.inv2w_te=(S.w_te_m>0.0)?0.5/S.w_te_m:0.0;

  // 4) Region target speeds.  V_sheath_LE_ms is finalized after the apex RH
  // state is known below; V_ME is an exact configured factor and is allowed to
  // be below V_sw.
  S.V_sheath_LE_ms=S.V_sw_ms;
  S.V_ME_ms=P_.V_ME_factor*S.V_sw_ms;
  S.V_dn_ms=S.V_sw_ms;

  // 5) Convenience mirrors.  sheath_comp_floor is no longer part of any
  // physical shock/region calculation; keep rc_floor=1 for source compatibility.
  S.inv_dr_sheath=(S.dr_sheath_m>0.0)?1.0/S.dr_sheath_m:0.0;
  S.rc_floor=1.0;

  // 6-7) Leblanc and Parker caches were prepared by swcme::core above.
  // The numbered placeholder is retained in comments because geometry remains
  // step 8 in older documentation/output traces.

  // 8) Geometry caches
  if (P_.shape==ShockShape::Ellipsoid){
    S.a_e    = S.r_sh_m;
    S.b_e    = S.a_e * P_.axis_ratio_y;
    S.c_e    = S.a_e * P_.axis_ratio_z;
    const double a2=S.a_e*S.a_e, b2=S.b_e*S.b_e, c2=S.c_e*S.c_e;
    S.inv_a2=(a2>0)?1.0/a2:0.0; S.inv_b2=(b2>0)?1.0/b2:0.0; S.inv_c2=(c2>0)?1.0/c2:0.0;
  } else if (P_.shape==ShockShape::SSE){
    // A true SSE front is the outward arc of a sphere that expands
    // self-similarly while preserving its angular half width lambda.  If R_a
    // is the apex distance, the generating sphere center is c along the CME
    // axis and its radius is a:
    //   c = R_a/(1+sin(lambda)),  a = c*sin(lambda).
    // The ray at alpha=lambda is tangent to this sphere; therefore the model
    // has a mathematically finite angular extent without artificial clamping.
    // half_width_rad was already checked against (0,pi/2] by validate().
    S.sin_half_width = std::sin(P_.half_width_rad);
    S.cos_half_width = std::cos(P_.half_width_rad);
    const double denom = 1.0 + S.sin_half_width;
    S.sse_center_m = S.r_sh_m/denom;
    S.sse_radius_m = S.sse_center_m*S.sin_half_width;
  }

  // The apex diagnostic below intentionally reuses guarded public geometry and
  // shock routines.  Install a provisional seal now that all fields those
  // routines read are complete; the few diagnostic mirrors they produce are
  // included when the final seal is written afterward.
  S.integrity_digest_=prepared_state_integrity(S);

  // 9) Apex shock diagnostic.  A geometric CME front and a physical fast
  // shock are not synonymous; cache both the explicit existence flag and the
  // physical compression for quick time-series diagnostics.
  {
    double u_apex[3]={e1[0],e1[1],e1[2]};
    LocalShockState apex;
    const bool surface_exists=shock_state_direction(S,u_apex,apex);
    S.has_shock=surface_exists && apex.has_shock && apex.solver_converged;
    S.rc=S.has_shock ? apex.compression : 1.0;
    if (S.has_shock) {
      const double V2_rad=apex.downstream.velocity_m_s[0]*S.e1[0]
                         +apex.downstream.velocity_m_s[1]*S.e1[1]
                         +apex.downstream.velocity_m_s[2]*S.e1[2];
      S.V_sheath_LE_ms=swcme::regions::leading_edge_speed(
          S.V_sw_ms,V2_rad,P_.V_sheath_LE_factor);
    }
  }
  // Replace the provisional seal with the final record digest after the apex
  // shock mirrors have been populated.  No StepState field changes afterward.
  S.integrity_digest_=prepared_state_integrity(S);
  // Freeze only after every geometry and shock calculation succeeds.  Failed
  // preparation therefore does not strand a model in a locked setup state.
  configuration_locked_.store(true,std::memory_order_release);
  return S;
}

// Radius and normal along direction (ux,uy,uz).
//
// The return value is part of the physical geometry contract.  Infinite
// Sun-centered verification geometries (Sphere/Ellipsoid) intersect every
// outward ray, but the SSE shock is a finite cap: directions beyond the
// configured half width have no shock surface and must not be assigned a
// clamped/fabricated flank radius.
bool Model::shape_radius_normal(const StepState& S,
                                double ux,double uy,double uz,
                                double& Rdir_m,double n_hat[3]) const {
  // Ownership precedes output initialization: a rejected foreign-state call
  // must leave Rdir_m and n_hat exactly as supplied by the caller.
  swcme::throw_if_error(validate_prepared_state(
      S,"swcme3d::shape_radius_normal"));
  double u[3]={ux,uy,uz};
  if (!::normalize_checked(u)) {
    throw std::invalid_argument(
        swcme::ModelStatus::make(swcme::StatusCode::DegenerateVector,
                                 "swcme3d::shape_radius_normal direction").summary());
  }

  // Always initialize the outputs to an explicitly nonphysical state.  This
  // ensures a caller that correctly checks the boolean cannot accidentally
  // reuse a previous finite radius/normal when the finite SSE cap is absent.
  Rdir_m=0.0;
  n_hat[0]=0.0; n_hat[1]=0.0; n_hat[2]=0.0;

  // Decompose the ray in the apex-aligned orthonormal frame.  u1=cos(alpha),
  // where alpha is angular separation from the CME propagation axis.
  const double u1=u[0]*S.e1[0]+u[1]*S.e1[1]+u[2]*S.e1[2];
  const double u2=u[0]*S.e2[0]+u[1]*S.e2[1]+u[2]*S.e2[2];
  const double u3=u[0]*S.e3[0]+u[1]*S.e3[1]+u[2]*S.e3[2];

  switch(P_.shape){
    case ShockShape::Sphere:{
      Rdir_m=S.r_sh_m;
      n_hat[0]=u[0]; n_hat[1]=u[1]; n_hat[2]=u[2];
      return true;
    }

    case ShockShape::Ellipsoid:{
      // All three axes scale with the same apex scale, so this remains a
      // self-similar surface.  The ray/ellipsoid intersection lambda follows
      // directly from the implicit quadratic level set.
      const double denom=(u1*u1)*S.inv_a2+(u2*u2)*S.inv_b2+(u3*u3)*S.inv_c2;
      const double lam=(denom>0)? 1.0/std::sqrt(denom):0.0;
      Rdir_m=lam;
      const double x=lam*u1,y=lam*u2,z=lam*u3;
      double nloc[3]={x*S.inv_a2,y*S.inv_b2,z*S.inv_c2};
      double ng[3]={ nloc[0]*S.e1[0]+nloc[1]*S.e2[0]+nloc[2]*S.e3[0],
                     nloc[0]*S.e1[1]+nloc[1]*S.e2[1]+nloc[2]*S.e3[1],
                     nloc[0]*S.e1[2]+nloc[1]*S.e2[2]+nloc[2]*S.e3[2] };
      if (!::normalize_checked(ng)) {
        throw std::runtime_error(
            swcme::ModelStatus::make(swcme::StatusCode::GeometryFailure,
                                     "swcme3d::shape_radius_normal ellipsoid normal").summary());
      }
      n_hat[0]=ng[0]; n_hat[1]=ng[1]; n_hat[2]=ng[2];
      return true;
    }

    case ShockShape::SSE:{
      // True self-similar-expansion (SSE) spherical-cap geometry.
      //
      // Let c be the center distance of the generating sphere from the Sun,
      // a its radius, and alpha the angle between the query ray and CME axis.
      // Intersecting |R*u-c*e1|^2=a^2 gives the outward root
      //
      //   R = c*cos(alpha) + sqrt(a^2-c^2*sin^2(alpha)).
      //
      // The discriminant is zero at alpha=lambda, exactly the tangent flank.
      // For alpha>lambda there is no physical intersection and we return false
      // instead of extending the shock with the legacy clamped radius.
      const double cos_alpha=std::max(-1.0,std::min(1.0,u1));
      const double alpha=std::acos(cos_alpha);
      const double eps=std::numeric_limits<double>::epsilon();
      const double angle_tol=128.0*eps*std::max(1.0,std::fabs(P_.half_width_rad));
      if (alpha>P_.half_width_rad+angle_tol) return false;

      const double c=S.sse_center_m;
      const double a=S.sse_radius_m;
      const double sin2_alpha=std::max(0.0,1.0-cos_alpha*cos_alpha);
      double discriminant=a*a-c*c*sin2_alpha;

      // At the tangent flank the exact discriminant is zero.  Roundoff in the
      // trigonometric projection can make it slightly negative, so only a
      // scale-aware, machine-level negative residual is clamped to zero.  A
      // materially negative value is treated as no intersection rather than
      // silently manufacturing a point.
      const double disc_scale=std::max({a*a,c*c,1.0});
      const double disc_tol=256.0*eps*disc_scale;
      if (discriminant<0.0) {
        if (discriminant>=-disc_tol) discriminant=0.0;
        else return false;
      }

      Rdir_m=c*cos_alpha+std::sqrt(discriminant);

      // The exact outward normal is the normalized gradient of the spherical
      // level set F(x)=|x-c*e1|^2-a^2.  Since x lies on the sphere, division by
      // a already produces a unit vector up to roundoff; normalize_checked removes
      // the remaining drift and reports a true degeneracy rather than inventing
      // an arbitrary normal.
      double ng[3]={
          Rdir_m*u[0]-c*S.e1[0],
          Rdir_m*u[1]-c*S.e1[1],
          Rdir_m*u[2]-c*S.e1[2]};
      if (a<=0.0) return false;  // guarded in prepare_step(); defensive only
      ng[0]/=a; ng[1]/=a; ng[2]/=a;
      if (!::normalize_checked(ng)) {
        throw std::runtime_error(
            swcme::ModelStatus::make(swcme::StatusCode::GeometryFailure,
                                     "swcme3d::shape_radius_normal SSE normal").summary());
      }
      n_hat[0]=ng[0]; n_hat[1]=ng[1]; n_hat[2]=ng[2];
      return true;
    }
  }

  return false;  // defensive for future enum extensions
}

// Complete local shock-state evaluation.
//
// Geometry and shock physics are deliberately separated.  The caller supplies
// only a direction; this routine locates the physical shock surface, evaluates
// the upstream Parker/Leblanc state AT THAT SURFACE, computes the self-similar
// normal shock speed, and then delegates the jump conditions to the shared
// ideal-MHD Rankine-Hugoniot solver in swcme_shock.hpp.
//
// This corrects two important legacy behaviors:
//   * a finite CME surface no longer implies that a fast shock exists; and
//   * shock strength no longer depends on the radius of an arbitrary field
//     query ahead of the shock.
swcme::ModelStatus Model::shock_state_direction_checked(
    const StepState& S, const double u_in[3], LocalShockState& state) const {
  // Check ownership before resetting `state`; PST02 requires a rejected call
  // to preserve every caller-owned output byte and perform no geometry work.
  const swcme::ModelStatus ownership=validate_prepared_state(
      S,"swcme3d::shock_state_direction");
  if (!ownership.ok()) return ownership;
  state=LocalShockState{};
  if (!u_in || !std::isfinite(u_in[0]) || !std::isfinite(u_in[1]) ||
      !std::isfinite(u_in[2])) {
    state.status=swcme::ModelStatus::make(
        swcme::StatusCode::NonFiniteInput,
        "swcme3d::shock_state_direction direction");
    return state.status;
  }

  double u[3]={u_in[0],u_in[1],u_in[2]};
  if (!::normalize_checked(u)) {
    state.status=swcme::ModelStatus::make(
        swcme::StatusCode::DegenerateVector,
        "swcme3d::shock_state_direction direction");
    return state.status;
  }

  double Rdir=0.0,n_hat[3]={0.0,0.0,0.0};
  try {
    if (!shape_radius_normal(S,u[0],u[1],u[2],Rdir,n_hat)) {
      state.status=swcme::ModelStatus::make(
          swcme::StatusCode::NoSurface,
          "swcme3d::shock_state_direction finite surface");
      return state.status;
    }
  } catch (const std::exception&) {
    state.status=swcme::ModelStatus::make(
        swcme::StatusCode::GeometryFailure,
        "swcme3d::shock_state_direction geometry");
    return state.status;
  }

  if (!std::isfinite(Rdir) || Rdir<swcme::solarwind::MIN_RADIUS_M ||
      !::finite3(n_hat)) {
    state.status=swcme::ModelStatus::make_value(
        Rdir<swcme::solarwind::MIN_RADIUS_M
            ? swcme::StatusCode::OutsideModelDomain
            : swcme::StatusCode::GeometryFailure,
        "swcme3d::shock_state_direction surface radius",Rdir);
    return state.status;
  }

  state.surface_exists=true;
  state.Rdir_m=Rdir;
  state.normal[0]=n_hat[0]; state.normal[1]=n_hat[1]; state.normal[2]=n_hat[2];

  // The upstream state is sampled at the physical shock surface.  No radius
  // clamp is performed here: a front inside the documented solar-wind domain
  // is an explicit OUTSIDE_MODEL_DOMAIN result.
  const double n_up_m3=swcme::solarwind::density_m3(S.common.solar_wind,Rdir);
  if (!std::isfinite(n_up_m3) || !(n_up_m3>0.0)) {
    state.status=swcme::ModelStatus::make(
        swcme::StatusCode::NonFiniteResult,
        "swcme3d::shock_state_direction upstream density");
    return state.status;
  }
  state.upstream_n_m3=n_up_m3;

  double B_up[3]={0.0,0.0,0.0};
  ::parker_vec_T_fast(S,u,Rdir,B_up);
  if (!::finite3(B_up)) {
    state.status=swcme::ModelStatus::make(
        swcme::StatusCode::NonFiniteResult,
        "swcme3d::shock_state_direction upstream magnetic field");
    return state.status;
  }

  const double radial_scale=(S.r_sh_m>0.0)? Rdir/S.r_sh_m : 0.0;
  double normal_projection=n_hat[0]*u[0]+n_hat[1]*u[1]+n_hat[2]*u[2];
  const double projection_tol=128.0*std::numeric_limits<double>::epsilon();
  if (normal_projection<0.0 && normal_projection>=-projection_tol)
    normal_projection=0.0; // roundoff-only clamp of a mathematically nonnegative dot product
  state.Vsh_n_m_s=S.V_sh_ms*radial_scale*normal_projection;
  if (!std::isfinite(state.Vsh_n_m_s)) {
    state.status=swcme::ModelStatus::make(
        swcme::StatusCode::NonFiniteResult,
        "swcme3d::shock_state_direction normal shock speed");
    return state.status;
  }

  swcme::shock::PrimitiveState upstream;
  upstream.rho_kg_m3=n_up_m3*MP;
  upstream.pressure_Pa=swcme::solarwind::proton_pressure_Pa(
      S.common.solar_wind,n_up_m3);
  upstream.velocity_m_s={{S.V_sw_ms*u[0],S.V_sw_ms*u[1],S.V_sw_ms*u[2]}};
  upstream.magnetic_T={{B_up[0],B_up[1],B_up[2]}};
  state.upstream=upstream;

  const swcme::shock::JumpResult jump=swcme::shock::solve_ideal_mhd_fast_shock(
      upstream,{{n_hat[0],n_hat[1],n_hat[2]}},state.Vsh_n_m_s,P_.gamma_ad);

  state.has_shock=jump.has_shock;
  state.solver_converged=jump.solver_converged;
  state.compression=(jump.has_shock && jump.solver_converged)? jump.compression : 1.0;
  state.theta_Bn_rad=jump.theta_Bn_rad;
  state.fast_speed_m_s=jump.fast_speed_m_s;
  state.fast_mach=jump.fast_mach;
  state.downstream=(jump.has_shock && jump.solver_converged)? jump.downstream : upstream;
  state.downstream_n_m3=(jump.has_shock && jump.solver_converged)
                           ? jump.compression*n_up_m3 : n_up_m3;
  state.mass_residual=jump.mass_residual;
  state.normal_B_residual=jump.normal_B_residual;
  state.electric_residual=jump.electric_residual;
  state.momentum_residual=jump.momentum_residual;
  state.energy_residual=jump.energy_residual;
  state.entropy_ratio=jump.entropy_ratio;

  if (jump.has_shock && !jump.solver_converged) {
    state.status=swcme::ModelStatus::make(
        swcme::StatusCode::ShockSolverFailure,
        swcme::shock::solve_status_name(jump.status));
    return state.status;
  }

  state.status=swcme::ModelStatus::success();
  return state.status;
}

bool Model::shock_state_direction(const StepState& S, const double u[3],
                                  LocalShockState& state) const {
  const swcme::ModelStatus status=shock_state_direction_checked(S,u,state);
  if (status.no_surface()) return false;
  swcme::throw_if_error(status);
  return true;
}

swcme::ModelStatus Model::shock_acceleration_state_checked(
    const StepState& S, const double u_in[3],
    swcme::acceleration::ShockAccelerationState& state) const {
  // Guard before clearing the acceleration record.  Direct callers therefore
  // receive the same transactional mismatch behavior as the SEP adapter.
  const swcme::ModelStatus ownership=validate_prepared_state(
      S,"swcme3d::shock_acceleration_state");
  if (!ownership.ok()) return ownership;
  state=swcme::acceleration::ShockAccelerationState{};
  LocalShockState shock;
  const swcme::ModelStatus shock_status=shock_state_direction_checked(S,u_in,shock);
  if (!shock_status.ok()) return shock_status;

  double u[3]={u_in[0],u_in[1],u_in[2]};
  if (!::normalize_checked(u)) {
    return swcme::ModelStatus::make(
        swcme::StatusCode::DegenerateVector,
        "swcme3d::shock_acceleration_state direction");
  }
  const std::array<double,3> position{{
      shock.Rdir_m*u[0],shock.Rdir_m*u[1],shock.Rdir_m*u[2]}};
  const std::array<double,3> normal{{
      shock.normal[0],shock.normal[1],shock.normal[2]}};
  const double Bmag=std::hypot(
      shock.upstream.magnetic_T[0],
      std::hypot(shock.upstream.magnetic_T[1],shock.upstream.magnetic_T[2]));
  if (!std::isfinite(Bmag)) {
    return swcme::ModelStatus::make(
        swcme::StatusCode::NonFiniteResult,
        "swcme3d::shock_acceleration_state upstream |B|");
  }
  state=swcme::acceleration::make_state(
      S.acceleration_config,true,shock.has_shock && shock.solver_converged,
      S.time_s,position,normal,shock.Vsh_n_m_s,shock.compression,
      shock.theta_Bn_rad,shock.fast_mach,shock.upstream_n_m3,Bmag);
  return swcme::ModelStatus::success();
}

bool Model::shock_acceleration_state(
    const StepState& S, const double u_in[3],
    swcme::acceleration::ShockAccelerationState& state) const {
  const swcme::ModelStatus status=
      shock_acceleration_state_checked(S,u_in,state);
  if (status.no_surface()) return false;
  swcme::throw_if_error(status);
  return status.ok();
}


// Return a Cartesian point on the analytical Parker field line anchored at an
// observer.  The production Parker field satisfies
//
//   B_phi/B_r = -Omega r sin(theta)/V_sw,
//
// while a field-line tangent satisfies
//
//   r sin(theta) dphi/dr = B_phi/B_r.
//
// Therefore dphi/dr=-Omega/V_sw and the colatitude is constant.  Integrating
// from the observer radius r_obs to a requested radius r gives
//
//   Delta phi = -Omega (r-r_obs)/V_sw.
//
// Rotating the observer radial unit vector by this angle around the configured
// solar axis gives the exact field line corresponding to parker_vec_T_fast().
// This analytical construction is both faster and more accurate than stepping
// a numerical field-line integrator through a field that already has a closed
// form.
bool Model::parker_field_line_point(const StepState& S,
                                    const double observer_m[3],
                                    double radius_m,
                                    double point_m[3]) const {
  // Validate before writing point_m.  The boolean compatibility API cannot
  // carry a status, so ownership misuse is reported through the same exception
  // policy used by the other legacy wrappers.
  swcme::throw_if_error(validate_prepared_state(
      S,"swcme3d::parker_field_line_point"));
  if (!observer_m || !point_m || !std::isfinite(radius_m) ||
      radius_m<swcme::solarwind::MIN_RADIUS_M ||
      !std::isfinite(S.V_sw_ms) || S.V_sw_ms<=0.0) {
    return false;
  }

  const double r_obs=norm3(observer_m);
  if (!std::isfinite(r_obs) || r_obs<swcme::solarwind::MIN_RADIUS_M) {
    return false;
  }

  const double u_obs[3]={observer_m[0]/r_obs,observer_m[1]/r_obs,
                         observer_m[2]/r_obs};
  const double delta_phi=-S.solar_rotation_rate_rad_s*(radius_m-r_obs)/S.V_sw_ms;
  double u[3]={0.0,0.0,0.0};
  rotate_about_unit_axis(u_obs,S.solar_axis_hat,delta_phi,u);
  if (!::normalize_checked(u)) return false; // orthogonal rotation should preserve norm

  point_m[0]=radius_m*u[0];
  point_m[1]=radius_m*u[1];
  point_m[2]=radius_m*u[2];
  return std::isfinite(point_m[0]) && std::isfinite(point_m[1]) &&
         std::isfinite(point_m[2]);
}

// Closed-form arc length on the same Parker line used above.  With constant
// colatitude theta,
//
//   ds/dr = sqrt(1 + (k r)^2),  k=Omega sin(theta)/V_sw.
//
// An antiderivative is
//
//   F(r)=0.5 [ r sqrt(1+(kr)^2) + asinh(kr)/k ].
//
// The exact k=0 (zero rotation or polar field line) limit is |r_b-r_a|.  This
// quantity is carried with every cobpoint because SEP timing depends on path
// length along the field, not merely on radial separation.
double Model::parker_field_line_length(const StepState& S,
                                       const double observer_m[3],
                                       double radius_a_m,
                                       double radius_b_m) const {
  // A NaN return is reserved for invalid numerical arguments.  A foreign
  // state is a programming/ownership error and is rejected explicitly before
  // any cached Parker coefficient is read.
  swcme::throw_if_error(validate_prepared_state(
      S,"swcme3d::parker_field_line_length"));
  if (!observer_m || !std::isfinite(radius_a_m) || !std::isfinite(radius_b_m) ||
      radius_a_m<=0.0 || radius_b_m<=0.0 || !std::isfinite(S.V_sw_ms) ||
      S.V_sw_ms<=0.0) {
    return std::numeric_limits<double>::quiet_NaN();
  }
  const double r_obs=norm3(observer_m);
  if (!std::isfinite(r_obs) || r_obs<=0.0)
    return std::numeric_limits<double>::quiet_NaN();

  const double u[3]={observer_m[0]/r_obs,observer_m[1]/r_obs,observer_m[2]/r_obs};
  const double cross[3]={
      S.solar_axis_hat[1]*u[2]-S.solar_axis_hat[2]*u[1],
      S.solar_axis_hat[2]*u[0]-S.solar_axis_hat[0]*u[2],
      S.solar_axis_hat[0]*u[1]-S.solar_axis_hat[1]*u[0]};
  const double sin_theta=norm3(cross);
  const double k=S.solar_rotation_rate_rad_s*sin_theta/S.V_sw_ms;

  if (std::abs(k)<=64.0*std::numeric_limits<double>::epsilon()/
                         std::max(radius_a_m,radius_b_m)) {
    return std::abs(radius_b_m-radius_a_m);
  }

  const auto primitive=[k](double r) {
    const double kr=k*r;
    return 0.5*(r*std::sqrt(1.0+kr*kr)+std::asinh(kr)/k);
  };
  return std::abs(primitive(radius_b_m)-primitive(radius_a_m));
}

// Intersect one observer-anchored Parker line with the current production
// shock surface.  The algorithm deliberately separates DISCOVERY from ROOT
// REFINEMENT:
//
//   1. scan the inward field line and evaluate h(r)=r-R_shock[u(r)] wherever
//      the selected finite geometry exists;
//   2. refine every sign-changing bracket by bisection;
//   3. independently minimize |h| around local minima so tangent roots, which
//      need not change sign, are not missed;
//   4. inspect finite-SSE validity transitions because first/last connection
//      can occur exactly at the angular boundary;
//   5. deduplicate roots and retain ALL physical intersections in radial order.
//
// The selected cobpoint is the OUTERMOST root.  It is the first shock surface
// encountered when tracing inward from the observer and is therefore a stable,
// documented default if a future non-convex geometry creates multiple roots.
swcme3d::ConnectivityState Model::observer_connectivity(
    const StepState& S, const double observer_m[3],
    const ConnectivityOptions& options) const {
  // Connectivity currently returns a value object rather than ModelStatus.
  // Reject a foreign state before initializing that object or tracing a single
  // field-line point so mixed-model geometry cannot resemble disconnection.
  swcme::throw_if_error(validate_prepared_state(
      S,"swcme3d::observer_connectivity"));
  ConnectivityState result;
  if (observer_m) {
    result.observer_position_m[0]=observer_m[0];
    result.observer_position_m[1]=observer_m[1];
    result.observer_position_m[2]=observer_m[2];
  }

  if (!observer_m || !std::isfinite(observer_m[0]) ||
      !std::isfinite(observer_m[1]) || !std::isfinite(observer_m[2])) {
    result.status=ConnectivityStatus::InvalidObserver;
    return result;
  }

  const double r_obs=norm3(observer_m);
  result.observer_radius_m=r_obs;
  if (!std::isfinite(r_obs) || r_obs<=0.0) {
    result.status=ConnectivityStatus::InvalidObserver;
    return result;
  }
  if (!std::isfinite(S.V_sw_ms) || S.V_sw_ms<=0.0 ||
      !std::isfinite(options.inner_radius_m) || options.inner_radius_m<=0.0 ||
      !std::isfinite(options.radius_tolerance_m) || options.radius_tolerance_m<=0.0 ||
      !std::isfinite(options.surface_residual_tolerance_m) ||
      options.surface_residual_tolerance_m<=0.0) {
    result.status=ConnectivityStatus::InvalidConfiguration;
    return result;
  }

  const double r_min=options.inner_radius_m;
  if (r_min>=r_obs) {
    result.status=ConnectivityStatus::Disconnected;
    return result;
  }

  struct Sample {
    double r=0.0;
    double h=0.0;
    bool surface=false;
  };

  // Evaluate the intersection residual using ONLY public production geometry.
  // A false surface flag is a geometrical no-surface state (e.g. outside an
  // SSE cap), not a numerical failure and must not be converted to a radius.
  const auto evaluate_residual=[&](double r)->Sample {
    Sample sample;
    sample.r=r;
    double x[3]={0.0,0.0,0.0};
    if (!parker_field_line_point(S,observer_m,r,x)) return sample;
    const double invr=1.0/r;
    const double u[3]={x[0]*invr,x[1]*invr,x[2]*invr};
    double R=0.0,n[3]={0.0,0.0,0.0};
    sample.surface=shape_radius_normal(S,u[0],u[1],u[2],R,n);
    if (sample.surface) sample.h=r-R;
    return sample;
  };

  // A minimum angular sampling of 0.5 degree per Parker rotation step prevents
  // tightly wound field lines from crossing a finite cap between sparse radial
  // samples.  The user value remains a lower bound, and an upper cap prevents
  // accidental pathological allocations for an invalid/extreme wind speed.
  const double total_phase=S.solar_rotation_rate_rad_s*(r_obs-r_min)/S.V_sw_ms;
  const double phase_step=0.5*PI/180.0;
  std::size_t intervals=std::max<std::size_t>(32,options.scan_intervals);
  if (std::isfinite(total_phase) && total_phase>0.0) {
    const std::size_t phase_intervals=static_cast<std::size_t>(
        std::ceil(total_phase/phase_step));
    intervals=std::max(intervals,phase_intervals);
  }
  intervals=std::min<std::size_t>(intervals,200000);

  std::vector<Sample> samples(intervals+1);
  for (std::size_t i=0;i<=intervals;++i) {
    const double f=static_cast<double>(i)/static_cast<double>(intervals);
    samples[i]=evaluate_residual(r_min+(r_obs-r_min)*f);
  }

  std::vector<double> candidate_roots;
  const double rtol=options.radius_tolerance_m;
  const double htol=options.surface_residual_tolerance_m;

  const auto add_candidate=[&](double r) {
    if (!std::isfinite(r) || r<r_min-rtol || r>r_obs+rtol) return;
    for (double old : candidate_roots) {
      if (std::abs(old-r)<=std::max(10.0*rtol,2.0*htol)) return;
    }
    candidate_roots.push_back(r);
  };

  // Refine an ordinary sign-changing root.  Bisection is chosen deliberately:
  // the residual is inexpensive, the interval is already small after scanning,
  // and bisection cannot jump to a different shock branch near a cap boundary.
  const auto refine_bracket=[&](double a,double b,Sample fa,Sample fb) {
    if (!fa.surface || !fb.surface) return;
    if (std::abs(fa.h)<=htol) { add_candidate(a); return; }
    if (std::abs(fb.h)<=htol) { add_candidate(b); return; }
    if (fa.h*fb.h>0.0) return;
    for (int it=0; it<100 && b-a>rtol; ++it) {
      const double m=0.5*(a+b);
      const Sample fm=evaluate_residual(m);
      if (!fm.surface) break;  // should not occur inside a normal valid bracket
      if (std::abs(fm.h)<=htol) { a=b=m; fa=fb=fm; break; }
      if (fa.h*fm.h<=0.0) { b=m; fb=fm; }
      else { a=m; fa=fm; }
    }
    const Sample sa=evaluate_residual(a), sb=evaluate_residual(b);
    if (sa.surface && sb.surface) add_candidate(
        std::abs(sa.h)<=std::abs(sb.h) ? a : b);
  };

  for (std::size_t i=1;i<samples.size();++i) {
    const Sample& a=samples[i-1];
    const Sample& b=samples[i];
    if (a.surface && std::abs(a.h)<=htol) add_candidate(a.r);
    if (b.surface && std::abs(b.h)<=htol) add_candidate(b.r);
    if (a.surface && b.surface && a.h*b.h<0.0)
      refine_bracket(a.r,b.r,a,b);
  }

  // Golden-section minimization of |h| catches a tangent intersection for
  // which h touches zero and returns to the same sign.  The objective assigns
  // an infinite cost where a finite geometry does not exist, so a minimizer
  // cannot fabricate an SSE flank outside its supported angular width.
  const auto refine_tangent=[&](double a,double b) {
    constexpr double GR=0.6180339887498948482;
    double c=b-GR*(b-a), d=a+GR*(b-a);
    auto cost=[&](double r) {
      const Sample q=evaluate_residual(r);
      return q.surface ? std::abs(q.h) : std::numeric_limits<double>::infinity();
    };
    double fc=cost(c), fd=cost(d);
    for (int it=0; it<120 && b-a>rtol; ++it) {
      if (fc<fd) {
        b=d; d=c; fd=fc; c=b-GR*(b-a); fc=cost(c);
      } else {
        a=c; c=d; fc=fd; d=a+GR*(b-a); fd=cost(d);
      }
    }
    const double r=(fc<fd)?c:d;
    const Sample q=evaluate_residual(r);
    if (q.surface && std::abs(q.h)<=htol) add_candidate(r);
  };

  for (std::size_t i=1;i+1<samples.size();++i) {
    if (!samples[i].surface) continue;
    const double ai=std::abs(samples[i].h);
    const double al=samples[i-1].surface ? std::abs(samples[i-1].h)
                                         : std::numeric_limits<double>::infinity();
    const double ar=samples[i+1].surface ? std::abs(samples[i+1].h)
                                         : std::numeric_limits<double>::infinity();
    if (ai<=al && ai<=ar) refine_tangent(samples[i-1].r,samples[i+1].r);
  }

  // A finite SSE connection may be born exactly where the Parker direction
  // crosses the cap boundary.  At such a point the surface-existence flag can
  // switch between adjacent scan samples.  Refine that boolean transition and
  // test the valid-side limiting residual explicitly; this avoids classifying
  // a true first/last connection as an iteration failure.
  for (std::size_t i=1;i<samples.size();++i) {
    if (samples[i-1].surface==samples[i].surface) continue;
    double a=samples[i-1].r, b=samples[i].r;
    bool va=samples[i-1].surface;
    for (int it=0; it<100 && b-a>rtol; ++it) {
      const double m=0.5*(a+b);
      const bool vm=evaluate_residual(m).surface;
      if (vm==va) a=m; else b=m;
    }
    // Evaluate just inside the valid side of the transition.  Using the final
    // bisection endpoint rather than an arbitrary epsilon keeps the tolerance
    // tied to the same radial convergence criterion as ordinary roots.
    const double r_valid=va?a:b;
    const Sample q=evaluate_residual(r_valid);
    if (q.surface && std::abs(q.h)<=htol) add_candidate(r_valid);
  }

  std::sort(candidate_roots.begin(),candidate_roots.end());

  // Convert geometric roots into complete, auditable cobpoint records.  A
  // candidate is retained only if a fresh production geometry/shock query at
  // the refined radius satisfies the surface residual tolerance.
  for (double r : candidate_roots) {
    double x[3]={0.0,0.0,0.0};
    if (!parker_field_line_point(S,observer_m,r,x)) continue;
    const double u[3]={x[0]/r,x[1]/r,x[2]/r};
    LocalShockState shock;
    if (!shock_state_direction(S,u,shock)) continue;
    const double residual=r-shock.Rdir_m;
    if (std::abs(residual)>htol) continue;

    ConnectivityRoot root;
    root.radius_m=r;
    root.position_m[0]=x[0]; root.position_m[1]=x[1]; root.position_m[2]=x[2];
    root.surface_residual_m=residual;
    root.path_length_m=parker_field_line_length(S,observer_m,r,r_obs);
    root.shock=shock;
    result.roots.push_back(root);
  }

  if (result.roots.empty()) {
    result.status=ConnectivityStatus::Disconnected;
    result.connected=false;
    result.selected_root=0;
    return result;
  }

  result.connected=true;
  result.status=ConnectivityStatus::Connected;
  result.selected_root=result.roots.size()-1;  // outermost/observer-nearest root
  return result;
}

// Build a deterministic time history for a stationary observer.  Each sample
// is solved from scratch from the same analytical Parker line and production
// shock geometry; this intentionally avoids hidden state that could create or
// suppress connection hysteresis.  Any temporal continuity seen in the result
// must therefore come from the physical evolution of the shock itself.
std::vector<swcme3d::ConnectivityHistorySample>
Model::observer_connectivity_history(const std::vector<double>& times_s,
                                     const double observer_m[3],
                                     const ConnectivityOptions& options) const {
  std::vector<ConnectivityHistorySample> history;
  history.reserve(times_s.size());
  for (double t : times_s) {
    ConnectivityHistorySample sample;
    sample.time_s=t;
    const StepState S=prepare_step(t);
    sample.connectivity=observer_connectivity(S,observer_m,options);
    history.push_back(std::move(sample));
  }
  return history;
}

// Backward-compatible scalar shock diagnostic.  The legacy r_eval_m argument
// is intentionally ignored: shock physics is a property of the shock surface,
// not of the point where a caller happens to request n/V/B.
void Model::local_oblique_rc(const StepState& S, const double u[3], const double n_hat[3],
                             double Rdir_m, double r_eval_m,
                             double& rc_out, double& Vsh_n_out, double& thetaBn_out) const {
  (void)n_hat;
  (void)Rdir_m;
  (void)r_eval_m;
  LocalShockState state;
  const bool surface_exists=shock_state_direction(S,u,state);
  if (!surface_exists) {
    rc_out=1.0; Vsh_n_out=0.0; thetaBn_out=0.0;
    return;
  }
  rc_out=state.has_shock && state.solver_converged ? state.compression : 1.0;
  Vsh_n_out=state.Vsh_n_m_s;
  thetaBn_out=state.theta_Bn_rad;
}

// Common FULL_ICME field construction for one Cartesian ray sample.
//
// Both public 3-D evaluators call this helper, so n/V and n/V/B can no longer
// drift into different region definitions.  Geometry comes from
// swcme_regions::make_boundaries(local R_sh), which is crucial for finite SSE
// and ellipsoid flanks: the layer thickness is a fraction of the LOCAL front
// radius rather than an apex-sized absolute distance.
struct RegionalSample3D {
  double n_m3 = 0.0;
  double velocity_m_s[3] = {0.0,0.0,0.0};
  double magnetic_T[3] = {0.0,0.0,0.0};
};

static inline void ambient_sample_3d(const swcme3d::StepState& S,
                                     const double u[3], double r_m,
                                     bool need_B, RegionalSample3D& out) {
  out.n_m3=swcme::solarwind::density_m3(S.common.solar_wind,r_m);
  out.velocity_m_s[0]=S.V_sw_ms*u[0];
  out.velocity_m_s[1]=S.V_sw_ms*u[1];
  out.velocity_m_s[2]=S.V_sw_ms*u[2];
  if (need_B) {
    ::parker_vec_T_fast(S,u,r_m,out.magnetic_T);
  }
}

static inline void evaluate_region_sample_3d(
    const swcme3d::Params& P, const swcme3d::StepState& S,
    const swcme3d::LocalShockState& shock, bool surface_exists,
    const double u[3], double r_m, bool need_B, RegionalSample3D& out) {
  ambient_sample_3d(S,u,r_m,need_B,out);

  // SHOCK_ONLY intentionally leaves the transport-facing background untouched
  // everywhere.  Shock geometry/connectivity/source diagnostics are still
  // available independently through shock_state_direction()/connectivity APIs.
  if (S.region_config.mode==swcme::regions::Mode::ShockOnly ||
      !surface_exists) {
    return;
  }

  const swcme::regions::Boundaries b=
      swcme::regions::make_boundaries(shock.Rdir_m,S.region_config);
  const swcme::regions::Location loc=swcme::regions::locate(r_m,b);
  if (loc.region==swcme::regions::Region::Upstream) return;

  if (loc.region==swcme::regions::Region::ShockTransition) {
    // RESOLVED_COMPRESSION uses exactly one common C1 shock layer.  At the
    // outer edge this is the local Parker/Leblanc upstream state; at the inner
    // edge it is the exact surface-owned RH downstream state.  The same blend
    // weight and endpoints are used for n, V and B, and the identical helper is
    // used by 1-D, eliminating dimension-dependent numerical acceleration.
    if (shock.has_shock && shock.solver_converged) {
      out.n_m3=swcme::regions::lerp(out.n_m3,shock.downstream_n_m3,loc.blend);
      for (int k=0;k<3;++k) {
        out.velocity_m_s[k]=swcme::regions::lerp(
            out.velocity_m_s[k],shock.downstream.velocity_m_s[k],loc.blend);
        if (need_B) {
          out.magnetic_T[k]=swcme::regions::lerp(
              out.magnetic_T[k],shock.downstream.magnetic_T[k],loc.blend);
        }
      }
    }
    return;
  }

  auto sheath_state = [&](double rr, RegionalSample3D& state) {
    ambient_sample_3d(S,u,rr,need_B,state);
    // A geometric CME front can persist after the fast shock disappears.  In
    // that case no artificial sheath compression/jump is introduced; the
    // sheath portion remains ambient while the optional ejecta can still exist.
    if (!shock.has_shock || !shock.solver_converged) return;

    const double w=swcme::regions::sheath_profile_weight(
        rr,b,S.region_config.sheath_ramp_power);
    const double n_le=swcme::solarwind::density_m3(
        S.common.solar_wind,b.R_le_m);
    state.n_m3=swcme::regions::log_lerp_positive(
        shock.downstream_n_m3,n_le,w);

    const double V2_rad=shock.downstream.velocity_m_s[0]*u[0]
                       +shock.downstream.velocity_m_s[1]*u[1]
                       +shock.downstream.velocity_m_s[2]*u[2];
    const double Vle_mag=swcme::regions::leading_edge_speed(
        S.V_sw_ms,V2_rad,P.V_sheath_LE_factor);
    const double Vle[3]={Vle_mag*u[0],Vle_mag*u[1],Vle_mag*u[2]};
    for (int k=0;k<3;++k) {
      state.velocity_m_s[k]=swcme::regions::lerp(
          shock.downstream.velocity_m_s[k],Vle[k],w);
    }

    if (need_B) {
      double B_le[3];
      ::parker_vec_T_fast(S,u,b.R_le_m,B_le);
      for (int k=0;k<3;++k) {
        state.magnetic_T[k]=swcme::regions::lerp(
            shock.downstream.magnetic_T[k],B_le[k],w);
      }
    }
  };

  auto ejecta_state = [&](double rr, RegionalSample3D& state) {
    ambient_sample_3d(S,u,rr,need_B,state);
    state.n_m3=S.region_config.f_ME*state.n_m3;
    state.velocity_m_s[0]=S.region_config.V_ME_factor*S.V_sw_ms*u[0];
    state.velocity_m_s[1]=S.region_config.V_ME_factor*S.V_sw_ms*u[1];
    state.velocity_m_s[2]=S.region_config.V_ME_factor*S.V_sw_ms*u[2];
    // The baseline ejecta magnetic field remains Parker in this deliberately
    // simple phenomenology.  A flux-rope/ejecta field is outside the present
    // controlled SEP-study scope and can be added later behind this same API.
  };

  if (loc.region==swcme::regions::Region::Sheath) {
    sheath_state(r_m,out);
    return;
  }
  if (loc.region==swcme::regions::Region::LeadingTransition) {
    RegionalSample3D sheath,ejecta;
    sheath_state(r_m,sheath);
    ejecta_state(r_m,ejecta);
    out.n_m3=swcme::regions::lerp(sheath.n_m3,ejecta.n_m3,loc.blend);
    for (int k=0;k<3;++k) {
      out.velocity_m_s[k]=swcme::regions::lerp(
          sheath.velocity_m_s[k],ejecta.velocity_m_s[k],loc.blend);
      if (need_B) out.magnetic_T[k]=swcme::regions::lerp(
          sheath.magnetic_T[k],ejecta.magnetic_T[k],loc.blend);
    }
    return;
  }
  if (loc.region==swcme::regions::Region::Ejecta) {
    ejecta_state(r_m,out);
    return;
  }
  if (loc.region==swcme::regions::Region::TrailingTransition) {
    RegionalSample3D ejecta,ambient;
    ejecta_state(r_m,ejecta);
    ambient_sample_3d(S,u,r_m,need_B,ambient);
    out.n_m3=swcme::regions::lerp(ejecta.n_m3,ambient.n_m3,loc.blend);
    for (int k=0;k<3;++k) {
      out.velocity_m_s[k]=swcme::regions::lerp(
          ejecta.velocity_m_s[k],ambient.velocity_m_s[k],loc.blend);
      if (need_B) out.magnetic_T[k]=swcme::regions::lerp(
          ejecta.magnetic_T[k],ambient.magnetic_T[k],loc.blend);
    }
  }
  // PostICME and Upstream retain the ambient state initialized above.
}

// n, V evaluator (allocation-free; vectorization-friendly)
//
// shock_state_direction() always describes the exact mathematical RH surface.
// The transport field, however, follows the single selected acceleration mode:
// SOURCE+SHOCK_ONLY leaves V uncompressed, while RESOLVED_COMPRESSION+FULL_ICME
// maps upstream to the exact RH downstream state through the shared finite C1
// layer in swcme_regions.hpp.  The same profile is used by 1-D and 3-D so div(V)
// cannot acquire an artificial dimensionality-dependent shock strength.
swcme::ModelStatus Model::evaluate_cartesian_fast_checked(
    const StepState& S,
    const double* x_m,const double* y_m,const double* z_m,
    double* n_m3,double* Vx_ms,double* Vy_ms,double* Vz_ms,
    std::size_t N) const {
  // Ownership has highest precedence and is checked before N or pointers.
  // Consequently a foreign state cannot produce a vacuous OK result or alter
  // any element of a caller-provided output array.
  const swcme::ModelStatus ownership=validate_prepared_state(
      S,"swcme3d::evaluate_cartesian_fast");
  if (!ownership.ok()) return ownership;
  if (N==0) return swcme::ModelStatus::success();
  if (!x_m || !y_m || !z_m || !n_m3 || !Vx_ms || !Vy_ms || !Vz_ms) {
    return swcme::ModelStatus::make(
        swcme::StatusCode::NullPointer,"swcme3d::evaluate_cartesian_fast");
  }

  for (std::size_t i=0;i<N;++i) {
    const double x=x_m[i],y=y_m[i],z=z_m[i];
    if (!std::isfinite(x) || !std::isfinite(y) || !std::isfinite(z)) {
      return swcme::ModelStatus::make(
          swcme::StatusCode::NonFiniteInput,
          "swcme3d::evaluate_cartesian_fast coordinate",i);
    }
    const double r=std::hypot(x,std::hypot(y,z));
    if (!std::isfinite(r)) {
      return swcme::ModelStatus::make(
          swcme::StatusCode::NonFiniteInput,
          "swcme3d::evaluate_cartesian_fast radius",i);
    }
    if (r<swcme::solarwind::MIN_RADIUS_M) {
      return swcme::ModelStatus::make_value(
          swcme::StatusCode::OutsideModelDomain,
          "swcme3d::evaluate_cartesian_fast radius",r,i);
    }
    const double invr=1.0/r;
    const double u[3]={x*invr,y*invr,z*invr};

    LocalShockState shock;
    const swcme::ModelStatus shock_status=shock_state_direction_checked(S,u,shock);
    const bool surface_exists=shock_status.ok();
    if (shock_status.failure()) {
      swcme::ModelStatus out=shock_status;
      out.sample_index=i;
      return out;
    }

    RegionalSample3D sample;
    evaluate_region_sample_3d(P_,S,shock,surface_exists,u,r,false,sample);
    if (!std::isfinite(sample.n_m3) || sample.n_m3<0.0 ||
        !std::isfinite(sample.velocity_m_s[0]) ||
        !std::isfinite(sample.velocity_m_s[1]) ||
        !std::isfinite(sample.velocity_m_s[2])) {
      return swcme::ModelStatus::make(
          swcme::StatusCode::NonFiniteResult,
          "swcme3d::evaluate_cartesian_fast regional state",i);
    }
    n_m3[i]=sample.n_m3;
    Vx_ms[i]=sample.velocity_m_s[0];
    Vy_ms[i]=sample.velocity_m_s[1];
    Vz_ms[i]=sample.velocity_m_s[2];
  }
  return swcme::ModelStatus::success();
}

void Model::evaluate_cartesian_fast(const StepState& S,
                                    const double* x_m,const double* y_m,const double* z_m,
                                    double* n_m3,double* Vx_ms,double* Vy_ms,double* Vz_ms,
                                    std::size_t N) const {
  swcme::throw_if_error(evaluate_cartesian_fast_checked(
      S,x_m,y_m,z_m,n_m3,Vx_ms,Vy_ms,Vz_ms,N));
}

swcme::ModelStatus Model::evaluate_cartesian_with_B_checked(
    const StepState& S,
    const double* x_m,const double* y_m,const double* z_m,
    double* n_m3,double* Vx_ms,double* Vy_ms,double* Vz_ms,
    double* Bx_T,double* By_T,double* Bz_T,std::size_t N) const {
  // Repeat the outer ownership guard here because callers may invoke this API
  // directly; rejection must occur before any background-field output changes.
  const swcme::ModelStatus ownership=validate_prepared_state(
      S,"swcme3d::evaluate_cartesian_with_B");
  if (!ownership.ok()) return ownership;
  if (N==0) return swcme::ModelStatus::success();
  if (!x_m || !y_m || !z_m || !n_m3 || !Vx_ms || !Vy_ms || !Vz_ms ||
      !Bx_T || !By_T || !Bz_T) {
    return swcme::ModelStatus::make(
        swcme::StatusCode::NullPointer,"swcme3d::evaluate_cartesian_with_B");
  }

  for (std::size_t i=0;i<N;++i) {
    const double x=x_m[i],y=y_m[i],z=z_m[i];
    if (!std::isfinite(x) || !std::isfinite(y) || !std::isfinite(z)) {
      return swcme::ModelStatus::make(
          swcme::StatusCode::NonFiniteInput,
          "swcme3d::evaluate_cartesian_with_B coordinate",i);
    }
    const double r=std::hypot(x,std::hypot(y,z));
    if (!std::isfinite(r)) {
      return swcme::ModelStatus::make(
          swcme::StatusCode::NonFiniteInput,
          "swcme3d::evaluate_cartesian_with_B radius",i);
    }
    if (r<swcme::solarwind::MIN_RADIUS_M) {
      return swcme::ModelStatus::make_value(
          swcme::StatusCode::OutsideModelDomain,
          "swcme3d::evaluate_cartesian_with_B radius",r,i);
    }
    const double invr=1.0/r;
    const double u[3]={x*invr,y*invr,z*invr};

    LocalShockState shock;
    const swcme::ModelStatus shock_status=shock_state_direction_checked(S,u,shock);
    const bool surface_exists=shock_status.ok();
    if (shock_status.failure()) {
      swcme::ModelStatus out=shock_status;
      out.sample_index=i;
      return out;
    }

    RegionalSample3D sample;
    evaluate_region_sample_3d(P_,S,shock,surface_exists,u,r,true,sample);
    if (!std::isfinite(sample.n_m3) || sample.n_m3<0.0 ||
        !std::isfinite(sample.velocity_m_s[0]) ||
        !std::isfinite(sample.velocity_m_s[1]) ||
        !std::isfinite(sample.velocity_m_s[2]) ||
        !std::isfinite(sample.magnetic_T[0]) ||
        !std::isfinite(sample.magnetic_T[1]) ||
        !std::isfinite(sample.magnetic_T[2])) {
      return swcme::ModelStatus::make(
          swcme::StatusCode::NonFiniteResult,
          "swcme3d::evaluate_cartesian_with_B regional state",i);
    }
    n_m3[i]=sample.n_m3;
    Vx_ms[i]=sample.velocity_m_s[0];
    Vy_ms[i]=sample.velocity_m_s[1];
    Vz_ms[i]=sample.velocity_m_s[2];
    Bx_T[i]=sample.magnetic_T[0];
    By_T[i]=sample.magnetic_T[1];
    Bz_T[i]=sample.magnetic_T[2];
  }
  return swcme::ModelStatus::success();
}

void Model::evaluate_cartesian_with_B(const StepState& S,
    const double* x_m,const double* y_m,const double* z_m,
    double* n_m3,double* Vx_ms,double* Vy_ms,double* Vz_ms,
    double* Bx_T,double* By_T,double* Bz_T,std::size_t N) const {
  swcme::throw_if_error(evaluate_cartesian_with_B_checked(
      S,x_m,y_m,z_m,n_m3,Vx_ms,Vy_ms,Vz_ms,Bx_T,By_T,Bz_T,N));
}

swcme::ModelStatus Model::compute_divV_cartesian_checked(
    const StepState& S,
    const double* x_m,const double* y_m,const double* z_m,
    double* divV,std::size_t N,double dr_frac) const {
  // Reject mixed ownership before validating the stencil or touching divV;
  // otherwise an invalid numerical-step status could hide the root API error.
  const swcme::ModelStatus ownership=validate_prepared_state(
      S,"swcme3d::compute_divV_cartesian");
  if (!ownership.ok()) return ownership;
  if (N==0) return swcme::ModelStatus::success();
  if (!x_m || !y_m || !z_m || !divV) {
    return swcme::ModelStatus::make(
        swcme::StatusCode::NullPointer,"swcme3d::compute_divV_cartesian");
  }
  if (!std::isfinite(dr_frac) || dr_frac<=0.0) {
    return swcme::ModelStatus::make_value(
        swcme::StatusCode::InvalidNumericalStep,
        "swcme3d::compute_divV_cartesian dr_frac",dr_frac);
  }

  // The absolute floor prevents a vanishing stencil near the inner boundary;
  // elsewhere the relative step controls spatial resolution.  No radial-flow
  // assumption is made: all three diagonal Jacobian entries are differentiated
  // in Cartesian coordinates from the actual production velocity evaluator.
  constexpr double MIN_CARTESIAN_STEP_M=1.0e3;
  for (std::size_t i=0;i<N;++i) {
    const double x=x_m[i],y=y_m[i],z=z_m[i];
    if (!std::isfinite(x) || !std::isfinite(y) || !std::isfinite(z)) {
      return swcme::ModelStatus::make(
          swcme::StatusCode::NonFiniteInput,
          "swcme3d::compute_divV_cartesian coordinate",i);
    }
    const double r=std::hypot(x,std::hypot(y,z));
    if (!std::isfinite(r) || r<swcme::solarwind::MIN_RADIUS_M) {
      return swcme::ModelStatus::make_value(
          std::isfinite(r) ? swcme::StatusCode::OutsideModelDomain
                           : swcme::StatusCode::NonFiniteInput,
          "swcme3d::compute_divV_cartesian radius",r,i);
    }
    const double h=std::max(MIN_CARTESIAN_STEP_M,dr_frac*r);
    if (!std::isfinite(h) || !(h>0.0)) {
      return swcme::ModelStatus::make(
          swcme::StatusCode::InvalidNumericalStep,
          "swcme3d::compute_divV_cartesian stencil",i);
    }

    const std::array<double,3> point{{x,y,z}};
    auto velocity_evaluator=[&](const std::array<double,3>& q,
                                std::array<double,3>& v) -> swcme::ModelStatus {
      double n=0.0;
      return evaluate_cartesian_fast_checked(
          S,&q[0],&q[1],&q[2],&n,&v[0],&v[1],&v[2],1);
    };

    double value=0.0;
    swcme::ModelStatus status=swcme::divergence::cartesian_second_order(
        point,h,velocity_evaluator,value);
    if (!status.ok()) {
      status.sample_index=i;
      return status;
    }
    divV[i]=value;
  }
  return swcme::ModelStatus::success();
}

void Model::compute_divV_cartesian(const StepState& S,
    const double* x_m,const double* y_m,const double* z_m,
    double* divV,std::size_t N,double dr_frac) const {
  swcme::throw_if_error(compute_divV_cartesian_checked(
      S,x_m,y_m,z_m,divV,N,dr_frac));
}

swcme::ModelStatus Model::compute_divV_checked(
    const StepState& S,
    const double* x_m,const double* y_m,const double* z_m,
    double* divV,std::size_t N,double dr_frac) const {
  // This API has an analytical branch that would otherwise read the foreign
  // region mode and solar-wind speed without entering another checked method.
  const swcme::ModelStatus ownership=validate_prepared_state(
      S,"swcme3d::compute_divV");
  if (!ownership.ok()) return ownership;
  if (N==0) return swcme::ModelStatus::success();
  if (!x_m || !y_m || !z_m || !divV) {
    return swcme::ModelStatus::make(
        swcme::StatusCode::NullPointer,"swcme3d::compute_divV");
  }

  // In the controlled SHOCK_ONLY baseline the velocity is exactly
  // V_sw e_r with constant magnitude, so the spherical divergence 2 V_sw/r is
  // analytical.  Using it directly removes finite-difference noise from the
  // adiabatic-energy-change term and makes the result rotationally invariant
  // to roundoff.  dr_frac is intentionally irrelevant in this branch.
  if (S.region_config.mode==swcme::regions::Mode::ShockOnly) {
    for (std::size_t i=0;i<N;++i) {
      const double x=x_m[i],y=y_m[i],z=z_m[i];
      if (!std::isfinite(x) || !std::isfinite(y) || !std::isfinite(z)) {
        return swcme::ModelStatus::make(
            swcme::StatusCode::NonFiniteInput,
            "swcme3d::compute_divV coordinate",i);
      }
      const double r=std::hypot(x,std::hypot(y,z));
      if (!std::isfinite(r) || r<swcme::solarwind::MIN_RADIUS_M) {
        return swcme::ModelStatus::make_value(
            std::isfinite(r) ? swcme::StatusCode::OutsideModelDomain
                             : swcme::StatusCode::NonFiniteInput,
            "swcme3d::compute_divV radius",r,i);
      }
      const swcme::divergence::RadialTerms terms=
          swcme::divergence::radial_terms(r,S.V_sw_ms,0.0);
      if (!std::isfinite(terms.total_s_inv)) {
        return swcme::ModelStatus::make(
            swcme::StatusCode::NonFiniteResult,
            "swcme3d::compute_divV analytical divergence",i);
      }
      divV[i]=terms.total_s_inv;
    }
    return swcme::ModelStatus::success();
  }

  // FULL_ICME is not generally radial: oblique RH jumps contain tangential
  // velocity and a finite SSE/ellipsoid introduces angular gradients in region
  // location.  The old ray derivative omitted both effects, so the canonical
  // path must use the full Cartesian Jacobian trace.
  return compute_divV_cartesian_checked(S,x_m,y_m,z_m,divV,N,dr_frac);
}

void Model::compute_divV(const StepState& S,
    const double* x_m,const double* y_m,const double* z_m,
    double* divV,std::size_t N,double dr_frac) const {
  swcme::throw_if_error(compute_divV_checked(
      S,x_m,y_m,z_m,divV,N,dr_frac));
}

swcme::ModelStatus Model::compute_divV_radial_checked(
    const StepState& S,
    const double* x_m,const double* y_m,const double* z_m,
    double* divV,std::size_t N,double dr_frac) const {
  // Compatibility alias: preserve the existing symbol while fixing its
  // historical semantics.  Callers that truly need the numerical Cartesian
  // operator for a convergence study should use compute_divV_cartesian_checked.
  return compute_divV_checked(S,x_m,y_m,z_m,divV,N,dr_frac);
}

void Model::compute_divV_radial(const StepState& S,
    const double* x_m,const double* y_m,const double* z_m,
    double* divV,std::size_t N,double dr_frac) const {
  swcme::throw_if_error(compute_divV_radial_checked(
      S,x_m,y_m,z_m,divV,N,dr_frac));
}

swcme::ModelStatus Model::evaluate_cartesian_with_B_div_checked(
    const StepState& S,
    const double* x_m,const double* y_m,const double* z_m,
    double* n_m3,double* Vx_ms,double* Vy_ms,double* Vz_ms,
    double* Bx_T,double* By_T,double* Bz_T,double* divVsw,
    std::size_t N,double dr_frac) const {
  swcme::ModelStatus status=evaluate_cartesian_with_B_checked(
      S,x_m,y_m,z_m,n_m3,Vx_ms,Vy_ms,Vz_ms,Bx_T,By_T,Bz_T,N);
  if (!status.ok()) return status;
  return compute_divV_checked(S,x_m,y_m,z_m,divVsw,N,dr_frac);
}

void Model::evaluate_cartesian_with_B_div(const StepState& S,
    const double* x_m,const double* y_m,const double* z_m,
    double* n_m3,double* Vx_ms,double* Vy_ms,double* Vz_ms,
    double* Bx_T,double* By_T,double* Bz_T,double* divVsw,
    std::size_t N,double dr_frac) const {
  swcme::throw_if_error(evaluate_cartesian_with_B_div_checked(
      S,x_m,y_m,z_m,n_m3,Vx_ms,Vy_ms,Vz_ms,Bx_T,By_T,Bz_T,divVsw,N,dr_frac));
}

bool Model::diagnose_direction(const StepState& S,const double u[3],
  double& Rdir_m,double n_hat[3],double& rc_loc,double& Vsh_n) const {
  // All diagnostics are now obtained from the same canonical surface-owned
  // shock state used by the Cartesian field evaluators and connectivity code.
  // This is intentionally stronger than calling shape_radius_normal() followed
  // by a second, scalar shock calculation: one production query establishes the
  // surface point, samples the upstream plasma at that point, and solves the MHD
  // jump.  Consequently a diagnostic cannot accidentally acquire a different
  // Mach number merely because some caller also supplied a field-sampling radius.
  LocalShockState state;
  if (!shock_state_direction(S,u,state)) {
    Rdir_m=0.0;
    n_hat[0]=n_hat[1]=n_hat[2]=0.0;
    rc_loc=1.0;
    Vsh_n=0.0;
    return false;
  }

  Rdir_m=state.Rdir_m;
  n_hat[0]=state.normal[0];
  n_hat[1]=state.normal[1];
  n_hat[2]=state.normal[2];
  rc_loc=(state.has_shock && state.solver_converged)? state.compression : 1.0;
  Vsh_n=state.Vsh_n_m_s;
  return true;
}

// Build a topologically unique triangular shock surface.
//
// The legacy implementation allocated a rectangular (theta,phi) array with
// both endpoints of the periodic phi interval and with an entire nPhi+1 ring at
// theta=0.  All apex-ring entries were the same physical point and phi=0/2pi
// were duplicate seam points.  Triangulating that array necessarily produced
// zero-area apex cells and a duplicated seam.  Filtering those cells after the
// fact is not acceptable because it leaves ambiguous adjacency and double
// counts nodes in source integration.
//
// The construction below encodes the topology directly:
//   * finite SSE cap: one apex + nTheta periodic rings, the last ring being the
//     physical half-width boundary;
//   * Sphere/Ellipsoid: one north/apex pole + (nTheta-1) periodic rings + one
//     south/rear pole.
// Each periodic ring contains exactly nPhi vertices at phi=2*pi*k/nPhi,
// k=0..nPhi-1.  The seam is closed only through the wrapped index (k+1)%nPhi.
ShockMesh Model::build_shock_mesh(const StepState& S,std::size_t nTheta,std::size_t nPhi) const {
  // No mesh allocation or node mutation occurs until ownership is proven.
  // This prevents a foreign state from creating a geometrically plausible
  // surface using the receiving model's shape parameters.
  swcme::throw_if_error(validate_prepared_state(
      S,"swcme3d::build_shock_mesh"));
  ShockMesh M;
  if (nTheta<3) nTheta=3;
  if (nPhi<3) nPhi=3;
  M.n_theta_intervals=nTheta;
  M.n_phi=nPhi;
  M.closed_surface=(P_.shape!=ShockShape::SSE);

  const double thetaMax=M.closed_surface ? swcme3d::PI : P_.half_width_rad;
  if (!std::isfinite(thetaMax) || !(thetaMax>0.0)) {
    throw std::runtime_error(
        swcme::ModelStatus::make(swcme::StatusCode::GeometryFailure,
                                 "swcme3d::build_shock_mesh angular extent").summary());
  }

  // Convert one local polar direction into a canonical production shock-state
  // node.  Keeping all node fields in this helper guarantees that geometry,
  // normal, compression, and normal speed come from the same surface-owned
  // LocalShockState used by diagnostics/connectivity.
  auto add_node = [&](double theta,double phi)->int {
    const double ct=std::cos(theta), st=std::sin(theta);
    const double cp=std::cos(phi), sp=std::sin(phi);
    const double u_loc[3]={ct,st*cp,st*sp};
    double u[3]={ u_loc[0]*S.e1[0]+u_loc[1]*S.e2[0]+u_loc[2]*S.e3[0],
                  u_loc[0]*S.e1[1]+u_loc[1]*S.e2[1]+u_loc[2]*S.e3[1],
                  u_loc[0]*S.e1[2]+u_loc[1]*S.e2[2]+u_loc[2]*S.e3[2] };
    if (!::normalize_checked(u)) {
      throw std::runtime_error(
          swcme::ModelStatus::make(swcme::StatusCode::GeometryFailure,
                                   "swcme3d::build_shock_mesh direction").summary());
    }

    LocalShockState shock;
    const swcme::ModelStatus status=shock_state_direction_checked(S,u,shock);
    if (!status.ok()) {
      // Every requested mesh node is expected to lie on the configured surface.
      // NO_SURFACE or OUTSIDE_MODEL_DOMAIN is therefore a real construction
      // failure, not a reason to insert a zero-radius placeholder node.
      throw std::runtime_error(status.summary());
    }

    const double Rdir=shock.Rdir_m;
    const double xyz[3]={Rdir*u[0],Rdir*u[1],Rdir*u[2]};
    const double rc=shock.has_shock ? shock.compression : 1.0;
    if (!::finite3(xyz) || !::finite3(shock.normal) || !std::isfinite(rc) ||
        !std::isfinite(shock.Vsh_n_m_s)) {
      throw std::runtime_error(
          swcme::ModelStatus::make(swcme::StatusCode::NonFiniteResult,
                                   "swcme3d::build_shock_mesh node").summary());
    }

    const int index=static_cast<int>(M.x.size());
    M.x.push_back(xyz[0]); M.y.push_back(xyz[1]); M.z.push_back(xyz[2]);
    M.n_hat_x.push_back(shock.normal[0]);
    M.n_hat_y.push_back(shock.normal[1]);
    M.n_hat_z.push_back(shock.normal[2]);
    M.rc.push_back(rc);
    M.Vsh_n.push_back(shock.Vsh_n_m_s);
    return index;
  };

  auto add_triangle = [&](int a,int b,int c) {
    // Connectivity is stored as 1-based Tecplot indices.  Repeated indices are
    // forbidden here and independently rejected by compute_triangle_metrics().
    if (a==b || b==c || c==a) {
      throw std::runtime_error(
          swcme::ModelStatus::make(swcme::StatusCode::InvalidMesh,
                                   "swcme3d::build_shock_mesh repeated index").summary());
    }
    M.tri_i.push_back(a+1);
    M.tri_j.push_back(b+1);
    M.tri_k.push_back(c+1);
  };

  // The apex is a single physical vertex, independent of phi.
  const int apex=add_node(0.0,0.0);
  (void)apex; // documented as index zero, but keep construction explicit.

  // Ring offsets contain the 0-based index of phi=0 on each non-polar ring.
  // No phi=2pi vertex is ever generated.
  std::vector<int> ring_offset;
  const std::size_t last_ring_interval=M.closed_surface ? nTheta-1 : nTheta;
  ring_offset.reserve(last_ring_interval);
  for (std::size_t it=1; it<=last_ring_interval; ++it) {
    const double theta=thetaMax*(static_cast<double>(it)/static_cast<double>(nTheta));
    ring_offset.push_back(static_cast<int>(M.x.size()));
    for (std::size_t ip=0; ip<nPhi; ++ip) {
      const double phi=2.0*swcme3d::PI*(static_cast<double>(ip)/static_cast<double>(nPhi));
      add_node(theta,phi);
    }
  }

  auto ring_index = [&](std::size_t ring,std::size_t ip)->int {
    return ring_offset[ring]+static_cast<int>(ip%nPhi);
  };

  // Apex fan.  The local parameterization has dX/dtheta x dX/dphi pointing
  // outward, so [apex,phi,phi+1] has the desired orientation.
  for (std::size_t ip=0; ip<nPhi; ++ip) {
    const std::size_t jp=(ip+1)%nPhi;
    add_triangle(apex,ring_index(0,ip),ring_index(0,jp));
  }

  // Connect adjacent unique rings.  Two triangles per angular quadrilateral;
  // wrapping jp closes the seam without duplicating either endpoint.
  for (std::size_t ring=0; ring+1<ring_offset.size(); ++ring) {
    for (std::size_t ip=0; ip<nPhi; ++ip) {
      const std::size_t jp=(ip+1)%nPhi;
      const int i00=ring_index(ring,ip);
      const int i01=ring_index(ring,jp);
      const int i10=ring_index(ring+1,ip);
      const int i11=ring_index(ring+1,jp);
      add_triangle(i00,i10,i11);
      add_triangle(i00,i11,i01);
    }
  }

  if (M.closed_surface) {
    // Sphere and origin-centered ellipsoid have one unique rear pole at theta=pi.
    // The last fan uses [ring_i,south,ring_{i+1}], which preserves the same
    // outward winding as the rest of the parameterized surface.
    const int rear=add_node(thetaMax,0.0);
    const std::size_t last=ring_offset.size()-1;
    for (std::size_t ip=0; ip<nPhi; ++ip) {
      const std::size_t jp=(ip+1)%nPhi;
      add_triangle(ring_index(last,ip),rear,ring_index(last,jp));
    }
  }

  return M;
}

void Model::compute_triangle_metrics(const ShockMesh& M, TriMetrics& T) const {
  const std::size_t Nv=M.x.size();
  const std::size_t Ne=M.tri_i.size();
  const bool node_sizes=(M.y.size()==Nv && M.z.size()==Nv &&
                         M.n_hat_x.size()==Nv && M.n_hat_y.size()==Nv &&
                         M.n_hat_z.size()==Nv && M.rc.size()==Nv &&
                         M.Vsh_n.size()==Nv);
  if (Nv<3 || Ne==0 || !node_sizes || M.tri_j.size()!=Ne || M.tri_k.size()!=Ne) {
    throw std::runtime_error(
        swcme::ModelStatus::make(swcme::StatusCode::InvalidMesh,
                                 "swcme3d::compute_triangle_metrics sizes").summary());
  }

  T.area.assign(Ne,0.0); T.nx.assign(Ne,0.0); T.ny.assign(Ne,0.0); T.nz.assign(Ne,0.0);
  T.cx.assign(Ne,0.0); T.cy.assign(Ne,0.0); T.cz.assign(Ne,0.0);
  T.rc_mean.assign(Ne,0.0); T.Vsh_n_mean.assign(Ne,0.0);

  for (std::size_t e=0;e<Ne;++e){
    const int ia=M.tri_i[e]-1, ib=M.tri_j[e]-1, ic=M.tri_k[e]-1;
    if (ia<0 || ib<0 || ic<0 || ia==ib || ib==ic || ic==ia ||
        static_cast<std::size_t>(ia)>=Nv || static_cast<std::size_t>(ib)>=Nv ||
        static_cast<std::size_t>(ic)>=Nv) {
      throw std::runtime_error(
          swcme::ModelStatus::make(swcme::StatusCode::InvalidMesh,
                                   "swcme3d::compute_triangle_metrics connectivity",e).summary());
    }

    const double A[3]={M.x[ia],M.y[ia],M.z[ia]};
    const double B[3]={M.x[ib],M.y[ib],M.z[ib]};
    const double C[3]={M.x[ic],M.y[ic],M.z[ic]};
    if (!::finite3(A) || !::finite3(B) || !::finite3(C)) {
      throw std::runtime_error(
          swcme::ModelStatus::make(swcme::StatusCode::NonFiniteResult,
                                   "swcme3d::compute_triangle_metrics vertex",e).summary());
    }

    const double AB[3]={B[0]-A[0],B[1]-A[1],B[2]-A[2]};
    const double AC[3]={C[0]-A[0],C[1]-A[1],C[2]-A[2]};
    const double BC[3]={C[0]-B[0],C[1]-B[1],C[2]-B[2]};
    double n[3]={AB[1]*AC[2]-AB[2]*AC[1],
                 AB[2]*AC[0]-AB[0]*AC[2],
                 AB[0]*AC[1]-AB[1]*AC[0]};
    const double twiceA=std::hypot(n[0],std::hypot(n[1],n[2]));
    const double lab=std::hypot(AB[0],std::hypot(AB[1],AB[2]));
    const double lac=std::hypot(AC[0],std::hypot(AC[1],AC[2]));
    const double lbc=std::hypot(BC[0],std::hypot(BC[1],BC[2]));
    const double edge_scale=std::max(lab,std::max(lac,lbc));

    // A valid triangle needs area that is resolvable relative to its own edge
    // scale.  This rejects exact duplicate/polar cells and grossly ill-
    // conditioned connectivity without imposing a physics-dependent absolute
    // area cutoff.  The factor is deliberately far above one ulp yet far below
    // the angular aspect ratios used by supported meshes.
    const double area_tol=128.0*std::numeric_limits<double>::epsilon()*
                          edge_scale*edge_scale;
    if (!std::isfinite(twiceA) || !std::isfinite(edge_scale) ||
        !(edge_scale>0.0) || !(0.5*twiceA>area_tol)) {
      throw std::runtime_error(
          swcme::ModelStatus::make_value(swcme::StatusCode::InvalidMesh,
                                         "swcme3d::compute_triangle_metrics degenerate area",
                                         0.5*twiceA,e).summary());
    }

    n[0]/=twiceA; n[1]/=twiceA; n[2]/=twiceA;
    const double area=0.5*twiceA;
    const double cx=(A[0]+B[0]+C[0])/3.0;
    const double cy=(A[1]+B[1]+C[1])/3.0;
    const double cz=(A[2]+B[2]+C[2])/3.0;
    const double rc_mean=(M.rc[ia]+M.rc[ib]+M.rc[ic])/3.0;
    const double v_mean=(M.Vsh_n[ia]+M.Vsh_n[ib]+M.Vsh_n[ic])/3.0;

    // Compare the triangle winding with the average analytical nodal normal.
    // A negative/zero dot is an orientation error; silently flipping it here
    // would hide a broken topology from downstream source integration.
    double navg[3]={M.n_hat_x[ia]+M.n_hat_x[ib]+M.n_hat_x[ic],
                    M.n_hat_y[ia]+M.n_hat_y[ib]+M.n_hat_y[ic],
                    M.n_hat_z[ia]+M.n_hat_z[ib]+M.n_hat_z[ic]};
    if (!::normalize_checked(navg)) {
      throw std::runtime_error(
          swcme::ModelStatus::make(swcme::StatusCode::InvalidMesh,
                                   "swcme3d::compute_triangle_metrics analytic normal",e).summary());
    }
    const double orientation=n[0]*navg[0]+n[1]*navg[1]+n[2]*navg[2];
    if (!std::isfinite(orientation) || !(orientation>0.0)) {
      throw std::runtime_error(
          swcme::ModelStatus::make_value(swcme::StatusCode::InvalidMesh,
                                         "swcme3d::compute_triangle_metrics orientation",
                                         orientation,e).summary());
    }

    if (!std::isfinite(area) || !std::isfinite(cx) || !std::isfinite(cy) ||
        !std::isfinite(cz) || !std::isfinite(rc_mean) || !std::isfinite(v_mean)) {
      throw std::runtime_error(
          swcme::ModelStatus::make(swcme::StatusCode::NonFiniteResult,
                                   "swcme3d::compute_triangle_metrics",e).summary());
    }
    T.area[e]=area;
    T.nx[e]=n[0]; T.ny[e]=n[1]; T.nz[e]=n[2];
    T.cx[e]=cx; T.cy[e]=cy; T.cz[e]=cz;
    T.rc_mean[e]=rc_mean;
    T.Vsh_n_mean[e]=v_mean;
  }
}

AreaSamplingTable Model::build_area_sampling_table(const TriMetrics& T) const {
  AreaSamplingTable table;
  if (T.area.empty()) {
    throw std::runtime_error(
        swcme::ModelStatus::make(swcme::StatusCode::InvalidMesh,
                                 "swcme3d::build_area_sampling_table empty").summary());
  }

  // Long-double accumulation prevents the CDF from losing small positive cells
  // when the mesh has a broad area distribution.  Every cell must have finite,
  // strictly positive physical area; zero-area entries are topology defects and
  // are never silently skipped.
  long double total=0.0L;
  for (std::size_t i=0;i<T.area.size();++i) {
    const double a=T.area[i];
    if (!std::isfinite(a) || !(a>0.0)) {
      throw std::runtime_error(
          swcme::ModelStatus::make_value(swcme::StatusCode::InvalidMesh,
                                         "swcme3d::build_area_sampling_table area",
                                         a,i).summary());
    }
    total+=static_cast<long double>(a);
  }
  if (!(total>0.0L) || !std::isfinite(static_cast<double>(total))) {
    throw std::runtime_error(
        swcme::ModelStatus::make(swcme::StatusCode::NonFiniteResult,
                                 "swcme3d::build_area_sampling_table total area").summary());
  }

  table.total_area_m2=static_cast<double>(total);
  table.cumulative_probability.resize(T.area.size());
  long double cumulative=0.0L;
  for (std::size_t i=0;i<T.area.size();++i) {
    cumulative+=static_cast<long double>(T.area[i]);
    table.cumulative_probability[i]=static_cast<double>(cumulative/total);
    if (i>0 && !(table.cumulative_probability[i]>
                 table.cumulative_probability[i-1])) {
      throw std::runtime_error(
          swcme::ModelStatus::make(swcme::StatusCode::InvalidMesh,
                                   "swcme3d::build_area_sampling_table nonmonotone CDF",i).summary());
    }
  }
  // The exact mathematical endpoint is one.  Force only the final roundoff bit
  // so callers can audit `cdf.back()==1` exactly; no cell probability is changed
  // by more than floating-point normalization error.
  table.cumulative_probability.back()=1.0;
  return table;
}

std::size_t Model::sample_triangle_by_area(const AreaSamplingTable& table,
                                           double unit_uniform) const {
  if (!std::isfinite(unit_uniform) || unit_uniform<0.0 || !(unit_uniform<1.0) ||
      table.cumulative_probability.empty() ||
      table.cumulative_probability.back()!=1.0) {
    throw std::invalid_argument(
        "swcme3d::sample_triangle_by_area requires a valid area CDF and 0<=u<1");
  }
  const auto it=std::upper_bound(table.cumulative_probability.begin(),
                                 table.cumulative_probability.end(),unit_uniform);
  if (it==table.cumulative_probability.end()) {
    // With cdf.back()==1 and u<1 this branch is unreachable unless the table was
    // externally corrupted after construction; surface that corruption rather
    // than biasing the sample toward the last cell.
    throw std::runtime_error(
        swcme::ModelStatus::make(swcme::StatusCode::InvalidMesh,
                                 "swcme3d::sample_triangle_by_area corrupted CDF").summary());
  }
  return static_cast<std::size_t>(it-table.cumulative_probability.begin());
}

// --- Tecplot helpers ---------------------------------------------------------
// OUT06 validates the complete mesh record before any output access.  Size and
// connectivity defects are INVALID_MESH; non-finite model-produced fields are
// NONFINITE_RESULT.  Nodal normals, compression, and normal speed also retain
// their physical representation invariants so a finite but corrupt record is
// not accepted merely because printf could serialize it.
static swcme::ModelStatus validate_output_mesh_record(
    const swcme3d::ShockMesh& M,const char* context) {
  const std::size_t Nv=M.x.size(),Ne=M.tri_i.size();
  if (Nv<3 || Ne==0 || M.y.size()!=Nv || M.z.size()!=Nv ||
      M.n_hat_x.size()!=Nv || M.n_hat_y.size()!=Nv ||
      M.n_hat_z.size()!=Nv || M.rc.size()!=Nv || M.Vsh_n.size()!=Nv ||
      M.tri_j.size()!=Ne || M.tri_k.size()!=Ne) {
    return swcme::ModelStatus::make(swcme::StatusCode::InvalidMesh,context);
  }

  for (std::size_t i=0; i<Nv; ++i) {
    const double values[]={M.x[i],M.y[i],M.z[i],M.n_hat_x[i],M.n_hat_y[i],
                           M.n_hat_z[i],M.rc[i],M.Vsh_n[i]};
    for (double value : values) {
      if (!std::isfinite(value))
        return swcme::ModelStatus::make_value(
            swcme::StatusCode::NonFiniteResult,context,value,i);
    }
    const double normal_magnitude=std::hypot(
        M.n_hat_x[i],std::hypot(M.n_hat_y[i],M.n_hat_z[i]));
    if (std::abs(normal_magnitude-1.0)>1.0e-10)
      return swcme::ModelStatus::make_value(
          swcme::StatusCode::InvalidMesh,context,normal_magnitude,i);
    if (M.rc[i]<1.0)
      return swcme::ModelStatus::make_value(
          swcme::StatusCode::InvalidMesh,context,M.rc[i],i);
    if (M.Vsh_n[i]<0.0)
      return swcme::ModelStatus::make_value(
          swcme::StatusCode::InvalidMesh,context,M.Vsh_n[i],i);
  }

  // Connectivity is stored in Tecplot's one-based convention.  Validate the
  // exact serialized integers, including repeated nodes, instead of relying on
  // a supplied TriMetrics array to imply that topology was once valid.
  for (std::size_t e=0; e<Ne; ++e) {
    const int a=M.tri_i[e],b=M.tri_j[e],c=M.tri_k[e];
    if (a<1 || b<1 || c<1 || a==b || b==c || c==a ||
        static_cast<std::size_t>(a)>Nv ||
        static_cast<std::size_t>(b)>Nv ||
        static_cast<std::size_t>(c)>Nv) {
      return swcme::ModelStatus::make(
          swcme::StatusCode::InvalidMesh,context,e);
    }
  }
  return swcme::ModelStatus::success();
}

static inline bool metrics_all_empty(const swcme3d::TriMetrics& T) {
  return T.area.empty() && T.nx.empty() && T.ny.empty() && T.nz.empty() &&
         T.cx.empty() && T.cy.empty() && T.cz.empty() && T.rc_mean.empty() &&
         T.Vsh_n_mean.empty();
}

static inline bool metrics_complete_size(const swcme3d::TriMetrics& T,
                                         std::size_t Ne) {
  return T.area.size()==Ne && T.nx.size()==Ne && T.ny.size()==Ne &&
         T.nz.size()==Ne && T.cx.size()==Ne && T.cy.size()==Ne &&
         T.cz.size()==Ne && T.rc_mean.size()==Ne && T.Vsh_n_mean.size()==Ne;
}

// Compare caller-supplied metrics with a fresh canonical derivation from the
// exact mesh.  This detects stale metrics after coordinates, connectivity, rc,
// or Vsh_n were changed.  A scaled tolerance permits independently reproduced
// floating-point values while remaining far tighter than Tecplot's nine-digit
// serialization precision.
static inline bool metric_matches(double supplied,double canonical) {
  const double scale=std::max(1.0,std::max(std::abs(supplied),
                                           std::abs(canonical)));
  return std::abs(supplied-canonical)<=
      1024.0*std::numeric_limits<double>::epsilon()*scale;
}

static swcme::ModelStatus validate_output_metrics(
    const swcme3d::TriMetrics& supplied,
    const swcme3d::TriMetrics& canonical,const char* context) {
  const std::size_t Ne=canonical.area.size();
  if (!metrics_complete_size(supplied,Ne))
    return swcme::ModelStatus::make(swcme::StatusCode::InvalidMesh,context);

  for (std::size_t e=0; e<Ne; ++e) {
    const double values[]={supplied.area[e],supplied.nx[e],supplied.ny[e],
        supplied.nz[e],supplied.cx[e],supplied.cy[e],supplied.cz[e],
        supplied.rc_mean[e],supplied.Vsh_n_mean[e]};
    for (double value : values) {
      if (!std::isfinite(value))
        return swcme::ModelStatus::make_value(
            swcme::StatusCode::NonFiniteResult,context,value,e);
    }
    const double reference[]={canonical.area[e],canonical.nx[e],
        canonical.ny[e],canonical.nz[e],canonical.cx[e],canonical.cy[e],
        canonical.cz[e],canonical.rc_mean[e],canonical.Vsh_n_mean[e]};
    for (std::size_t field=0; field<9; ++field) {
      if (!metric_matches(values[field],reference[field]))
        return swcme::ModelStatus::make_value(
            swcme::StatusCode::InvalidMesh,context,values[field],e);
    }
  }
  return swcme::ModelStatus::success();
}

// Produce the metric record that will actually be written.  A completely empty
// TriMetrics requests canonical computation for convenience; a partially
// populated record is ambiguous and rejected.  Even complete supplied metrics
// are checked against a fresh derivation, which also reuses the production
// triangle-quality/orientation checks in compute_triangle_metrics().
static swcme::ModelStatus prepare_output_mesh(
    const swcme3d::Model& model,const swcme3d::ShockMesh& M,
    const swcme3d::TriMetrics& supplied,swcme3d::TriMetrics& output,
    const char* mesh_context,const char* metrics_context) {
  const swcme::ModelStatus mesh_status=
      validate_output_mesh_record(M,mesh_context);
  if (!mesh_status.ok()) return mesh_status;

  swcme3d::TriMetrics canonical;
  try {
    model.compute_triangle_metrics(M,canonical);
  } catch (...) {
    // The public metric builder already checks degeneracy and triangle winding.
    // Convert its exception boundary back into the checked writer's status
    // channel while still guaranteeing that no file has been opened.
    return swcme::ModelStatus::make(
        swcme::StatusCode::InvalidMesh,mesh_context);
  }

  if (metrics_all_empty(supplied)) {
    output=canonical;
    return swcme::ModelStatus::success();
  }
  const swcme::ModelStatus metric_status=
      validate_output_metrics(supplied,canonical,metrics_context);
  if (!metric_status.ok()) return metric_status;
  output=supplied;
  return swcme::ModelStatus::success();
}

// Multiply output-grid dimensions without allowing size_t wraparound.  The
// writers stream rows and therefore do not allocate Ni*Nj*Nk elements, but the
// flattened row counter and any downstream consumer still require that total
// to be representable.  Keeping the check here also makes an enormous invalid
// request fail immediately instead of entering a practically unbounded loop.
static inline bool checked_size_product(std::size_t left,std::size_t right,
                                        std::size_t& product) {
  if (right!=0 && left>std::numeric_limits<std::size_t>::max()/right)
    return false;
  product=left*right;
  return true;
}

// OUT04 defines one structural contract for every use of BoxSpec.  Center and
// half-extent fields must be finite, half extents cannot be negative, every
// documented grid dimension must contain at least two samples, all lower/upper
// bounds and doubled spans must be representable, and the complete point count
// must fit size_t.  This helper intentionally performs no model-domain scan;
// OUT05 owns the later per-point radius check.  Returning a precise status here
// before FileOperations is selected guarantees that malformed boxes are as
// side-effect free as out-of-domain points.
static swcme::ModelStatus validate_box_spec(
    const swcme3d::BoxSpec& B,const char* context) {
  const double centers[]={B.cx,B.cy,B.cz};
  const double extents[]={B.hx,B.hy,B.hz};
  const int dimensions[]={B.Ni,B.Nj,B.Nk};

  for (int axis=0; axis<3; ++axis) {
    if (!std::isfinite(centers[axis]))
      return swcme::ModelStatus::make_value(
          swcme::StatusCode::NonFiniteInput,context,centers[axis]);
    if (!std::isfinite(extents[axis]))
      return swcme::ModelStatus::make_value(
          swcme::StatusCode::NonFiniteInput,context,extents[axis]);
    if (extents[axis]<0.0)
      return swcme::ModelStatus::make_value(
          swcme::StatusCode::InvalidConfiguration,context,extents[axis]);
    if (dimensions[axis]<2)
      return swcme::ModelStatus::make_value(
          swcme::StatusCode::InvalidConfiguration,context,
          static_cast<double>(dimensions[axis]));

    // Checking both bounds is not sufficient: a symmetric interval may have
    // finite endpoints while 2*h overflows in the canonical grid expression.
    // Validate the exact intermediate required by structured_coordinate().
    const double lower=centers[axis]-extents[axis];
    const double upper=centers[axis]+extents[axis];
    const double span=2.0*extents[axis];
    if (!std::isfinite(lower) || !std::isfinite(upper) ||
        !std::isfinite(span)) {
      const double bad=!std::isfinite(lower)
          ? lower : (!std::isfinite(upper) ? upper : span);
      return swcme::ModelStatus::make_value(
          swcme::StatusCode::NonFiniteInput,context,bad);
    }
  }

  std::size_t plane_points=0,total_points=0;
  if (!checked_size_product(static_cast<std::size_t>(B.Ni),
                            static_cast<std::size_t>(B.Nj),plane_points) ||
      !checked_size_product(plane_points,static_cast<std::size_t>(B.Nk),
                            total_points)) {
    return swcme::ModelStatus::make(
        swcme::StatusCode::InvalidConfiguration,context);
  }
  return swcme::ModelStatus::success();
}

// Generate a structured-grid coordinate in one canonical place.  OUT05 uses
// this exact expression during preflight and the writers use it again during
// emission, preventing validation and output from silently sampling different
// points because of duplicated interpolation formulas.
static inline double structured_coordinate(double center,double half_extent,
                                           int index,int count) {
  return center+(-half_extent+(2.0*half_extent)*
      (index/double(std::max(1,count-1))));
}

// Validate one requested Cartesian output location without evaluating or
// modifying model state.  The evaluator contract treats non-finite Cartesian
// coordinates (including overflow in grid generation) as NONFINITE_INPUT and
// a finite radius below the shared solar-wind floor as OUTSIDE_MODEL_DOMAIN.
// Carrying the flattened output row makes the failure directly traceable to a
// Tecplot record before any staging file is created.
static inline swcme::ModelStatus preflight_output_point(
    double x,double y,double z,const char* context,std::size_t row) {
  if (!std::isfinite(x) || !std::isfinite(y) || !std::isfinite(z)) {
    const double bad=!std::isfinite(x) ? x : (!std::isfinite(y) ? y : z);
    return swcme::ModelStatus::make_value(
        swcme::StatusCode::NonFiniteInput,context,bad,row);
  }
  const double radius=std::hypot(x,std::hypot(y,z));
  if (!std::isfinite(radius))
    return swcme::ModelStatus::make_value(
        swcme::StatusCode::NonFiniteInput,context,radius,row);
  if (radius<swcme::solarwind::MIN_RADIUS_M)
    return swcme::ModelStatus::make_value(
        swcme::StatusCode::OutsideModelDomain,context,radius,row);
  return swcme::ModelStatus::success();
}

// Surface coordinates are already-computed model results, but their requested
// output locations must obey the same spatial floor as direct background
// queries.  Validate all nodes, not merely the first, so a malformed late mesh
// node cannot cause a filesystem side effect before rejection.
static swcme::ModelStatus preflight_surface_domain(
    const swcme3d::ShockMesh& M,const char* context) {
  for (std::size_t i=0; i<M.x.size(); ++i) {
    const swcme::ModelStatus status=
        preflight_output_point(M.x[i],M.y[i],M.z[i],context,i);
    if (!status.ok()) return status;
  }
  return swcme::ModelStatus::success();
}

// Scan every structured volume point in the exact K/J/I order emitted by the
// bundle.  A separate helper keeps this potentially large but allocation-free
// pass ahead of FileOperations selection and gives sample_index the same linear
// row meaning it has in the resulting Tecplot zone.
static swcme::ModelStatus preflight_volume_domain(
    const swcme3d::BoxSpec& B,const char* context) {
  std::size_t row=0;
  for (int kk=0; kk<B.Nk; ++kk) {
    const double z=structured_coordinate(B.cz,B.hz,kk,B.Nk);
    for (int jj=0; jj<B.Nj; ++jj) {
      const double y=structured_coordinate(B.cy,B.hy,jj,B.Nj);
      for (int ii=0; ii<B.Ni; ++ii,++row) {
        const double x=structured_coordinate(B.cx,B.hx,ii,B.Ni);
        const swcme::ModelStatus status=
            preflight_output_point(x,y,z,context,row);
        if (!status.ok()) return status;
      }
    }
  }
  return swcme::ModelStatus::success();
}

// OUT04 guarantees at least two BoxSpec points along every axis.  Retaining the
// established max(2,...) expression here keeps the preflight byte-compatible
// with the writer and defensive against any future internal caller that reaches
// this helper without passing through the public structural validator.
static swcme::ModelStatus preflight_min_x_face_domain(
    const swcme3d::BoxSpec& B,const char* context) {
  const int I=std::max(2,B.Nj),J=std::max(2,B.Nk);
  const double x=B.cx-B.hx;
  std::size_t row=0;
  for (int j=0; j<J; ++j) {
    const double z=structured_coordinate(B.cz,B.hz,j,J);
    for (int i=0; i<I; ++i,++row) {
      const double y=structured_coordinate(B.cy,B.hy,i,I);
      const swcme::ModelStatus status=
          preflight_output_point(x,y,z,context,row);
      if (!status.ok()) return status;
    }
  }
  return swcme::ModelStatus::success();
}

// Write one Tecplot BLOCK variable while retaining the element index on every
// formatted write.  Returning immediately on failure prevents later values or
// line breaks from obscuring the exact block element that was truncated.
static inline bool dump_array_block(
    swcme::output::CheckedTextFile& output,
    const std::vector<double>& values,const char* context) {
  int columns=0;
  for (std::size_t i=0; i<values.size(); ++i) {
    if (!output.print(context,i,"%.9e ",values[i])) return false;
    if (++columns==8) {
      if (!output.print(context,i,"\n")) return false;
      columns=0;
    }
  }
  return columns==0 || output.print(
      context,values.empty() ? swcme::ModelStatus::npos : values.size()-1,"\n");
}

// Reserved cell fields are emitted through the same checked path as physical
// values.  Even though their content is constant, each zero is real output and
// must carry an item index if the destination stops accepting bytes.
static inline bool dump_zeros_block(
    swcme::output::CheckedTextFile& output,
    std::size_t count,const char* context) {
  int columns=0;
  for (std::size_t i=0; i<count; ++i) {
    if (!output.print(context,i,"0 ")) return false;
    if (++columns==8) {
      if (!output.print(context,i,"\n")) return false;
      columns=0;
    }
  }
  return columns==0 || output.print(
      context,count==0 ? swcme::ModelStatus::npos : count-1,"\n");
}

// Emit the already validated surface dataset through the shared OUT02 stream
// lifecycle and OUT03 transaction.  Keeping this helper independent of
// validation lets both the legacy bool wrapper and checked status API use
// exactly the same staged bytes and commit behavior.
static swcme::ModelStatus write_surface_output(
    const swcme3d::ShockMesh& M,const swcme3d::TriMetrics& T,
    const char* path,const swcme::output::FileOperations& operations) {
  const std::size_t Nv=M.x.size(),Ne=M.tri_i.size();
  swcme::output::CheckedTextFile output(operations);
  if (!output.open_transactional(path)) return swcme::ModelStatus::make(
      swcme::StatusCode::FileOpenFailure,"3D surface output open");

  output.print("3D surface title",swcme::ModelStatus::npos,
               "TITLE=\"Shock surface (cell metrics + nodal rc)\"\n");
  output.print("3D surface variables",swcme::ModelStatus::npos,
      "VARIABLES=\"X\",\"Y\",\"Z\",\"n\",\"Vx\",\"Vy\",\"Vz\","
      "\"Bx\",\"By\",\"Bz\",\"divVsw\",\"rc\",\"Vsh_n\",\"nx\","
      "\"ny\",\"nz\",\"area\",\"rc_mean\",\"Vsh_n_mean\",\"tnx\","
      "\"tny\",\"tnz\",\"cx\",\"cy\",\"cz\"\n");
  output.print("3D surface zone",swcme::ModelStatus::npos,
      "ZONE T=\"surface_cells\", N=%zu, E=%zu, ZONETYPE=FETRIANGLE, "
      "DATAPACKING=BLOCK,\n",Nv,Ne);
  output.print("3D surface variable locations",swcme::ModelStatus::npos,
      "VARLOCATION=([1-3,12-13]=NODAL, [4-11,14-25]=CELLCENTERED)\n");

  dump_array_block(output,M.x,"3D surface X block");
  dump_array_block(output,M.y,"3D surface Y block");
  dump_array_block(output,M.z,"3D surface Z block");
  for (int field=0; field<8 && output.good(); ++field)
    dump_zeros_block(output,Ne,"3D surface zero field block");
  dump_array_block(output,M.rc,"3D surface compression block");
  dump_array_block(output,M.Vsh_n,"3D surface speed block");
  dump_array_block(output,T.nx,"3D surface normal-X block");
  dump_array_block(output,T.ny,"3D surface normal-Y block");
  dump_array_block(output,T.nz,"3D surface normal-Z block");
  dump_array_block(output,T.area,"3D surface area block");
  dump_array_block(output,T.rc_mean,"3D surface mean-compression block");
  dump_array_block(output,T.Vsh_n_mean,"3D surface mean-speed block");
  for (int field=0; field<3 && output.good(); ++field)
    dump_zeros_block(output,Ne,"3D surface reserved field block");
  dump_array_block(output,T.cx,"3D surface centroid-X block");
  dump_array_block(output,T.cy,"3D surface centroid-Y block");
  dump_array_block(output,T.cz,"3D surface centroid-Z block");

  for (std::size_t e=0; e<Ne && output.good(); ++e)
    output.print("3D surface connectivity row",e,"%d %d %d\n",
                 M.tri_i[e],M.tri_j[e],M.tri_k[e]);

  return output.finish("3D surface output flush",
                       "3D surface output stream error",
                       "3D surface output close",
                       "3D surface output commit");
}

// Surface-only: cell metrics + nodal rc/Vsh_n
bool Model::write_shock_surface_center_metrics_tecplot(
  const ShockMesh& M, const TriMetrics& T_in, const char* path) const {
  // Delegate to the checked API so legacy callers receive the same OUT05
  // no-open rejection and OUT03 transaction as status-aware integrations.
  return write_shock_surface_center_metrics_tecplot_checked(
      M,T_in,path).ok();
}

swcme::ModelStatus Model::write_shock_surface_center_metrics_tecplot_checked(
    const ShockMesh& M,const TriMetrics& T,const char* path,
    const swcme::output::FileOperations* file_operations) const {
  if (!path) return swcme::ModelStatus::make(
      swcme::StatusCode::NullPointer,"write_shock_surface_center_metrics_tecplot path");
  TriMetrics checked;
  const swcme::ModelStatus mesh_status=prepare_output_mesh(
      *this,M,T,checked,"write_shock_surface_center_metrics_tecplot mesh",
      "write_shock_surface_center_metrics_tecplot metrics");
  if (!mesh_status.ok()) return mesh_status;
  const swcme::ModelStatus surface_domain=preflight_surface_domain(
      M,"write_shock_surface_center_metrics_tecplot vertex");
  if (!surface_domain.ok()) return surface_domain;
  const swcme::output::FileOperations& operations=file_operations
      ? *file_operations : swcme::output::stdio_file_operations();
  return write_surface_output(M,checked,path,operations);
}

// Default apex-aligned volume box
BoxSpec Model::default_apex_box(const StepState& S,double half_AU,int N) const {
  // The box center depends on prepared orientation and apex distance; reject
  // foreign caches before combining either value with this model's conventions.
  swcme::throw_if_error(validate_prepared_state(
      S,"swcme3d::default_apex_box"));
  // A factory should never manufacture a BoxSpec that its consumers reject.
  // Validate source units before conversion so negative/NaN requests and an
  // under-resolved grid fail at their origin with no partially valid result.
  if (!std::isfinite(half_AU))
    swcme::throw_if_error(swcme::ModelStatus::make_value(
        swcme::StatusCode::NonFiniteInput,
        "swcme3d::default_apex_box half_AU",half_AU));
  if (half_AU<0.0)
    swcme::throw_if_error(swcme::ModelStatus::make_value(
        swcme::StatusCode::InvalidConfiguration,
        "swcme3d::default_apex_box half_AU",half_AU));
  if (N<2)
    swcme::throw_if_error(swcme::ModelStatus::make_value(
        swcme::StatusCode::InvalidConfiguration,
        "swcme3d::default_apex_box resolution",static_cast<double>(N)));
  BoxSpec B; const double h=half_AU*AU;
  B.hx=h; B.hy=h; B.hz=h;
  const double shift=0.4*h; // move box outward along e1 so shock cuts through
  B.cx=S.a_m*S.e1[0]+shift*S.e1[0];
  B.cy=S.a_m*S.e1[1]+shift*S.e1[1];
  B.cz=S.a_m*S.e1[2]+shift*S.e1[2];
  B.Ni=N; B.Nj=N; B.Nk=N;
  // Conversion and outward shifting can overflow even when half_AU itself is
  // finite.  Reuse the public-writer contract so the factory and consumers
  // cannot drift apart as BoxSpec evolves.
  swcme::throw_if_error(validate_box_spec(
      B,"swcme3d::default_apex_box generated box"));
  return B;
}

// Write a complete four-zone dataset after the caller has validated its model
// state, mesh, metrics, and box.  All output passes through CheckedTextFile so
// a partial raw write and a delayed buffered-stream failure have the same
// explicit FILE_WRITE_FAILURE contract, and neither can replace the previous
// destination before OUT03's final commit.
static swcme::ModelStatus write_bundle_output(
    const swcme3d::Model& model,const swcme3d::ShockMesh& M,
    const swcme3d::TriMetrics& T,const swcme3d::StepState& S,
    const swcme3d::BoxSpec& B,const char* path,
    const swcme::output::FileOperations& operations) {
  const std::size_t Nv=M.x.size(),Ne=M.tri_i.size();
  swcme::output::CheckedTextFile output(operations);
  if (!output.open_transactional(path)) return swcme::ModelStatus::make(
      swcme::StatusCode::FileOpenFailure,"3D dataset bundle open");

  output.print("3D dataset title",swcme::ModelStatus::npos,
               "TITLE = \"SW+CME dataset\"\n");
  output.print("3D dataset variables",swcme::ModelStatus::npos,"VARIABLES = "
    "\"X\",\"Y\",\"Z\","
    "\"n\",\"Vx\",\"Vy\",\"Vz\","
    "\"Bx\",\"By\",\"Bz\",\"divVsw\","
    "\"rc\",\"Vsh_n\","
    "\"nx\",\"ny\",\"nz\",\"area\",\"rc_mean\",\"Vsh_n_mean\","
    "\"tnx\",\"tny\",\"tnz\",\"cx\",\"cy\",\"cz\"\n");

  // Zone 1: surface_cells (FETRIANGLE, BLOCK)
  output.print("3D dataset surface-cell zone",swcme::ModelStatus::npos,
      "ZONE T=\"surface_cells\", N=%zu, E=%zu, ZONETYPE=FETRIANGLE, "
      "DATAPACKING=BLOCK,\n",Nv,Ne);
  output.print("3D dataset surface-cell locations",swcme::ModelStatus::npos,
      "VARLOCATION=([1-3,12-13]=NODAL, [4-11,14-25]=CELLCENTERED)\n");

  dump_array_block(output,M.x,"3D dataset surface X block");
  dump_array_block(output,M.y,"3D dataset surface Y block");
  dump_array_block(output,M.z,"3D dataset surface Z block");
  for (int field=0; field<8 && output.good(); ++field)
    dump_zeros_block(output,Ne,"3D dataset surface zero field block");
  dump_array_block(output,M.rc,"3D dataset surface compression block");
  dump_array_block(output,M.Vsh_n,"3D dataset surface speed block");
  dump_array_block(output,T.nx,"3D dataset surface normal-X block");
  dump_array_block(output,T.ny,"3D dataset surface normal-Y block");
  dump_array_block(output,T.nz,"3D dataset surface normal-Z block");
  dump_array_block(output,T.area,"3D dataset surface area block");
  dump_array_block(output,T.rc_mean,"3D dataset surface mean-compression block");
  dump_array_block(output,T.Vsh_n_mean,"3D dataset surface mean-speed block");
  for (int field=0; field<3 && output.good(); ++field)
    dump_zeros_block(output,Ne,"3D dataset surface reserved field block");
  dump_array_block(output,T.cx,"3D dataset surface centroid-X block");
  dump_array_block(output,T.cy,"3D dataset surface centroid-Y block");
  dump_array_block(output,T.cz,"3D dataset surface centroid-Z block");

  for (std::size_t e=0; e<Ne && output.good(); ++e)
    output.print("3D dataset surface connectivity row",e,"%d %d %d\n",
                 M.tri_i[e],M.tri_j[e],M.tri_k[e]);

  // Zone 2: surface_nodal (FEPOINT)
  output.print("3D dataset surface-nodal zone",swcme::ModelStatus::npos,
      "ZONE T=\"surface_nodal\", N=%zu, E=%zu, F=FEPOINT, ET=TRIANGLE\n",
      Nv,Ne);
  for (std::size_t i=0; i<Nv && output.good(); ++i) {
    output.print("3D dataset surface-nodal row",i,"%.9e %.9e %.9e "
      "%.9e %.9e %.9e %.9e "
      "%.9e %.9e %.9e %.9e "
      "%.9e %.9e "
      "%.9e %.9e %.9e "
      "%.9e %.9e %.9e "
      "%.9e %.9e %.9e "
      "%.9e %.9e %.9e\n",
      M.x[i], M.y[i], M.z[i],
      0.0,0.0,0.0,0.0,
      0.0,0.0,0.0,0.0,
      M.rc[i], M.Vsh_n[i],
      M.n_hat_x[i], M.n_hat_y[i], M.n_hat_z[i],
      0.0,0.0,0.0,
      0.0,0.0,0.0,
      0.0,0.0,0.0
    );
  }
  for (std::size_t e=0; e<Ne && output.good(); ++e)
    output.print("3D dataset nodal connectivity row",e,"%d %d %d\n",
                 M.tri_i[e],M.tri_j[e],M.tri_k[e]);

  // Zone 3: volume_box (structured POINT)
  output.print("3D dataset volume zone",swcme::ModelStatus::npos,
      "ZONE T=\"volume_box\", I=%d, J=%d, K=%d, DATAPACKING=POINT\n",
      B.Ni,B.Nj,B.Nk);
  std::size_t volume_row=0;
  for (int kk=0; kk<B.Nk && output.good(); ++kk){
    const double zk=structured_coordinate(B.cz,B.hz,kk,B.Nk);
    for (int jj=0; jj<B.Nj && output.good(); ++jj){
      const double yj=structured_coordinate(B.cy,B.hy,jj,B.Nj);
      for (int ii=0; ii<B.Ni && output.good(); ++ii,++volume_row){
        const double xi=structured_coordinate(B.cx,B.hx,ii,B.Ni);
        double n,Vx,Vy,Vz,Bx,By,Bz,div;
        model.evaluate_cartesian_with_B(
            S,&xi,&yj,&zk,&n,&Vx,&Vy,&Vz,&Bx,&By,&Bz,1);
        model.compute_divV_radial(S,&xi,&yj,&zk,&div,1,1e-3);
        output.print("3D dataset volume row",volume_row,"%.9e %.9e %.9e "
          "%.9e %.9e %.9e %.9e "
          "%.9e %.9e %.9e %.9e "
          "%.9e %.9e "
          "%.9e %.9e %.9e "
          "%.9e %.9e %.9e "
          "%.9e %.9e %.9e "
          "%.9e %.9e %.9e\n",
          xi, yj, zk,
          n, Vx, Vy, Vz,
          Bx, By, Bz, div,
          0.0,0.0, 0.0,0.0,0.0, 0.0,0.0,0.0, 0.0,0.0,0.0, 0.0,0.0,0.0
        );
      }
    }
  }

  // Zone 4: min-X face (structured 2-D POINT)
  const int I=std::max(2,B.Nj), J=std::max(2,B.Nk);
  const double x0=B.cx-B.hx;
  output.print("3D dataset face zone",swcme::ModelStatus::npos,
      "ZONE T=\"box_face_minX\", I=%d, J=%d, DATAPACKING=POINT\n",I,J);
  std::size_t face_row=0;
  for (int j=0; j<J && output.good(); ++j){
    const double z=structured_coordinate(B.cz,B.hz,j,J);
    for (int i=0; i<I && output.good(); ++i,++face_row){
      const double y=structured_coordinate(B.cy,B.hy,i,I);
      const double x=x0;
      double n,Vx,Vy,Vz,Bx,By,Bz,div;
      model.evaluate_cartesian_with_B(
          S,&x,&y,&z,&n,&Vx,&Vy,&Vz,&Bx,&By,&Bz,1);
      model.compute_divV_radial(S,&x,&y,&z,&div,1,1e-3);
      output.print("3D dataset face row",face_row,"%.9e %.9e %.9e "
        "%.9e %.9e %.9e %.9e "
        "%.9e %.9e %.9e %.9e "
        "%.9e %.9e "
        "%.9e %.9e %.9e "
        "%.9e %.9e %.9e "
        "%.9e %.9e %.9e "
        "%.9e %.9e %.9e\n",
        x, y, z,
        n, Vx, Vy, Vz,
        Bx, By, Bz, div,
        0.0,0.0, 0.0,0.0,0.0, 0.0,0.0,0.0, 0.0,0.0,0.0, 0.0,0.0,0.0
      );
    }
  }

  return output.finish("3D dataset bundle flush",
                       "3D dataset bundle stream error",
                       "3D dataset bundle close",
                       "3D dataset bundle commit");
}

// Bundle writer: surface_cells + surface_nodal + volume_box + minX face
bool Model::write_tecplot_dataset_bundle(const ShockMesh& M,const TriMetrics& T_in,
                                         const StepState& S,const BoxSpec& B,
                                         const char* path) const {
  // The source-compatible boolean writer throws on ownership misuse before it
  // can open/truncate `path`.  I/O failures remain the historical false return.
  swcme::throw_if_error(validate_prepared_state(
      S,"swcme3d::write_tecplot_dataset_bundle"));
  return write_tecplot_dataset_bundle_checked(M,T_in,S,B,path).ok();
}

// Emit the standalone face using the same bytes and field order as zone four
// of the bundle.  The helper is separate so injected write and commit failures
// exercise the compiled 3-D transaction rather than only the generic utility.
static swcme::ModelStatus write_face_output(
    const swcme3d::Model& model,const swcme3d::StepState& S,
    const swcme3d::BoxSpec& B,const char* path,
    const swcme::output::FileOperations& operations) {
  const int I=std::max(2,B.Nj), J=std::max(2,B.Nk);
  const double x0=B.cx-B.hx;

  swcme::output::CheckedTextFile output(operations);
  if (!output.open_transactional(path)) return swcme::ModelStatus::make(
      swcme::StatusCode::FileOpenFailure,"3D box face open");
  output.print("3D box face title",swcme::ModelStatus::npos,
               "TITLE = \"Box face (minX)\"\n");
  output.print("3D box face variables",swcme::ModelStatus::npos,"VARIABLES = "
    "\"X\",\"Y\",\"Z\","
    "\"n\",\"Vx\",\"Vy\",\"Vz\","
    "\"Bx\",\"By\",\"Bz\",\"divVsw\","
    "\"rc\",\"Vsh_n\","
    "\"nx\",\"ny\",\"nz\",\"area\",\"rc_mean\",\"Vsh_n_mean\","
    "\"tnx\",\"tny\",\"tnz\",\"cx\",\"cy\",\"cz\"\n");
  output.print("3D box face zone",swcme::ModelStatus::npos,
      "ZONE T=\"box_face_minX\", I=%d, J=%d, DATAPACKING=POINT\n",I,J);

  std::size_t row=0;
  for (int j=0; j<J && output.good(); ++j){
    const double z=structured_coordinate(B.cz,B.hz,j,J);
    for (int i=0; i<I && output.good(); ++i,++row){
      const double y=structured_coordinate(B.cy,B.hy,i,I);
      const double x=x0;
      double n,Vx,Vy,Vz,Bx,By,Bz,div;
      model.evaluate_cartesian_with_B(
          S,&x,&y,&z,&n,&Vx,&Vy,&Vz,&Bx,&By,&Bz,1);
      model.compute_divV_radial(S,&x,&y,&z,&div,1,1e-3);
      output.print("3D box face row",row,"%.9e %.9e %.9e "
        "%.9e %.9e %.9e %.9e "
        "%.9e %.9e %.9e %.9e "
        "%.9e %.9e "
        "%.9e %.9e %.9e "
        "%.9e %.9e %.9e "
        "%.9e %.9e %.9e "
        "%.9e %.9e %.9e\n",
        x, y, z,
        n, Vx, Vy, Vz,
        Bx, By, Bz, div,
        0.0,0.0, 0.0,0.0,0.0, 0.0,0.0,0.0, 0.0,0.0,0.0, 0.0,0.0,0.0
      );
    }
  }
  return output.finish("3D box face flush","3D box face stream error",
                       "3D box face close","3D box face commit");
}

// Standalone 2-D face writer (min-X plane)
bool Model::write_box_face_minX_tecplot_structured(const StepState& S,
                                                   const BoxSpec& B,
                                                   const char* path) const {
  // Ownership precedes open so rejection preserves an existing output file
  // exactly and cannot leave a partial Tecplot header behind.
  swcme::throw_if_error(validate_prepared_state(
      S,"swcme3d::write_box_face_minX_tecplot_structured"));
  return write_box_face_minX_tecplot_structured_checked(S,B,path).ok();
}


swcme::ModelStatus Model::write_tecplot_dataset_bundle_checked(
    const ShockMesh& M,const TriMetrics& T,const StepState& S,const BoxSpec& B,
    const char* path,
    const swcme::output::FileOperations* file_operations) const {
  // Reject before inspecting unrelated mesh/box arguments and, critically,
  // before any destination file is opened.  This makes ownership the root
  // diagnostic and guarantees byte-preserving failure semantics.
  const swcme::ModelStatus ownership=validate_prepared_state(
      S,"swcme3d::write_tecplot_dataset_bundle");
  if (!ownership.ok()) return ownership;
  if (!path) return swcme::ModelStatus::make(
      swcme::StatusCode::NullPointer,"write_tecplot_dataset_bundle path");
  TriMetrics checked;
  const swcme::ModelStatus mesh_status=prepare_output_mesh(
      *this,M,T,checked,"write_tecplot_dataset_bundle mesh",
      "write_tecplot_dataset_bundle metrics");
  if (!mesh_status.ok()) return mesh_status;
  const swcme::ModelStatus box_status=validate_box_spec(
      B,"write_tecplot_dataset_bundle box");
  if (!box_status.ok()) return box_status;
  const swcme::ModelStatus surface_domain=preflight_surface_domain(
      M,"write_tecplot_dataset_bundle surface vertex");
  if (!surface_domain.ok()) return surface_domain;
  const swcme::ModelStatus volume_domain=preflight_volume_domain(
      B,"write_tecplot_dataset_bundle volume point");
  if (!volume_domain.ok()) return volume_domain;
  const swcme::ModelStatus face_domain=preflight_min_x_face_domain(
      B,"write_tecplot_dataset_bundle face point");
  if (!face_domain.ok()) return face_domain;
  try {
    const swcme::output::FileOperations& operations=file_operations
        ? *file_operations : swcme::output::stdio_file_operations();
    return write_bundle_output(*this,M,checked,S,B,path,operations);
  } catch (...) {
    return swcme::ModelStatus::make(swcme::StatusCode::NonFiniteResult,
                                    "write_tecplot_dataset_bundle evaluation");
  }
}

swcme::ModelStatus Model::write_box_face_minX_tecplot_structured_checked(
    const StepState& S,const BoxSpec& B,const char* path,
    const swcme::output::FileOperations* file_operations) const {
  // Validate the state before box checks or file-system access so a foreign
  // cache cannot be masked by a secondary argument error.
  const swcme::ModelStatus ownership=validate_prepared_state(
      S,"swcme3d::write_box_face_minX_tecplot_structured");
  if (!ownership.ok()) return ownership;
  if (!path) return swcme::ModelStatus::make(
      swcme::StatusCode::NullPointer,"write_box_face_minX_tecplot_structured path");
  const swcme::ModelStatus box_status=validate_box_spec(
      B,"write_box_face_minX_tecplot_structured box");
  if (!box_status.ok()) return box_status;
  const swcme::ModelStatus face_domain=preflight_min_x_face_domain(
      B,"write_box_face_minX_tecplot_structured point");
  if (!face_domain.ok()) return face_domain;
  try {
    const swcme::output::FileOperations& operations=file_operations
        ? *file_operations : swcme::output::stdio_file_operations();
    return write_face_output(*this,S,B,path,operations);
  } catch (...) {
    return swcme::ModelStatus::make(swcme::StatusCode::NonFiniteResult,
                                    "write_box_face_minX_tecplot_structured evaluation");
  }
}

} // namespace swcme3d
