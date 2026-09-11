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
// 5) Region structure and smoothing (radial blends):
//    Regions: upstream → shock → sheath → leading edge → magnetic ejecta → trailing edge.
//    We blend with a C^1 smoothstep s(x)=x^2(3-2x) applied to signed distances,
//    with independent widths (w_shock, w_le, w_te), each self-similar ∝ r_sh.
//    Sheath compression at the shock uses a shape:
//      rc_loc = rc_floor + (rc_oblique - rc_floor) * (1 - ξ)^p,  ξ = (Rdir - r)/dr_sheath,
//    where p = sheath_ramp_power ≥ 1; rc_oblique from an oblique-MHD proxy.
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
//    ∇·V = (1/r^2) d/dr ( r^2 V_r ). We compute it with a robust centered
//    finite difference at r±dr along the same ray, using evaluate_cartesian_fast.
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
//     - ∇·V (divergence of bulk speed) [1/s] via a robust radial finite-diff.
//     - a triangulated shock surface mesh with nodal normals, nodal rc, nodal
//       normal shock speed, and per-cell metrics (area, rc_mean, Vsh_n_mean,
//       centroid, geometric normals).
// • Tecplot dataset writers for the shock surface and for a structured volume
//   box near the apex (plus a 2-D face zone). Files are NaN/Inf-sanitized.
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
//    For each point, we compute ∇·V with a robust radial finite difference
//       ∇·V = (1/r^2) ∂(r^2 V_r)/∂r
//    using ±dr around r with dr = max(dr_min, dr_frac · r). We evaluate V at
//    r±dr along the same unit radial direction.
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
// 11  "divVsw" [1/s]  : divergence of V (radial FD)
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
//    All numeric outputs are passed through finite_or(...) to avoid NaN/Inf.
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
//  • All numeric outputs use finite_or(v, fallback) before printing
//    to avoid NaN/Inf in Tecplot files.
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
static inline double finite_or(double v, double fallback=0.0){
  return std::isfinite(v)? v : fallback;
}
static inline void safe_normalize(double v[3]){
  const double m=std::sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
  if (m>0){ v[0]/=m; v[1]/=m; v[2]/=m; } else { v[0]=1; v[1]=0; v[2]=0; }
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
  return std::sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
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
  using namespace swcme3d;

  const double r_AU = r_m / AU;
  const double Br = S.Br1AU_T / (r_AU*r_AU);

  // Cross product Omega_hat × e_r gives both the azimuthal direction and,
  // through its magnitude, the local sin(colatitude).  S.solar_axis_hat is
  // normalized in prepare_step(), while u is a unit radial direction supplied
  // by the point evaluator.
  const double cross[3] = {
      S.solar_axis_hat[1]*u[2] - S.solar_axis_hat[2]*u[1],
      S.solar_axis_hat[2]*u[0] - S.solar_axis_hat[0]*u[2],
      S.solar_axis_hat[0]*u[1] - S.solar_axis_hat[1]*u[0]};
  const double sin_theta_local =
      std::sqrt(cross[0]*cross[0] + cross[1]*cross[1] + cross[2]*cross[2]);

  // The cross product can vanish exactly at the rotation poles.  In that
  // physical limit the Parker winding is zero, so the transverse field is
  // exactly zero and no arbitrary e_phi direction should be manufactured.
  constexpr double AXIS_TOL = 64.0*std::numeric_limits<double>::epsilon();
  if (sin_theta_local <= AXIS_TOL) {
    B_out[0] = Br*u[0];
    B_out[1] = Br*u[1];
    B_out[2] = Br*u[2];
    return;
  }

  const double inv_sin_theta = 1.0/sin_theta_local;
  const double ephi[3] = {
      cross[0]*inv_sin_theta,
      cross[1]*inv_sin_theta,
      cross[2]*inv_sin_theta};

  // S.k_AU is the equatorial coefficient Omega*AU/V_sw.  Multiplication by
  // the local sin(theta) below supplies the latitude dependence required by
  // the Parker solution instead of reusing one global latitude everywhere.
  const double Bphi =
      -Br * S.k_AU * r_AU * sin_theta_local;

  B_out[0] = Br*u[0] + Bphi*ephi[0];
  B_out[1] = Br*u[1] + Bphi*ephi[1];
  B_out[2] = Br*u[2] + Bphi*ephi[2];
}

// ----------------------------------------------------------------------------
// Model implementation
// ----------------------------------------------------------------------------
namespace swcme3d {

Model::Model(const Params& P): P_(P) {}

StepState Model::prepare_step(double t_s) const {
  StepState S{};

  // 1) Apex-aligned orthonormal basis (e1 along CME apex direction)
  double e1[3]={P_.cme_dir[0],P_.cme_dir[1],P_.cme_dir[2]}; ::safe_normalize(e1);
  double tmp[3]={0,0,1}; if (std::fabs(e1[2])>0.9){ tmp[0]=1; tmp[1]=0; tmp[2]=0; }
  double e2[3]={ e1[1]*tmp[2]-e1[2]*tmp[1],
                 e1[2]*tmp[0]-e1[0]*tmp[2],
                 e1[0]*tmp[1]-e1[1]*tmp[0] }; ::safe_normalize(e2);
  double e3[3]={ e1[1]*e2[2]-e1[2]*e2[1],
                 e1[2]*e2[0]-e1[0]*e2[2],
                 e1[0]*e2[1]-e1[1]*e2[0] }; ::safe_normalize(e3);
  S.e1[0]=e1[0]; S.e1[1]=e1[1]; S.e1[2]=e1[2];
  S.e2[0]=e2[0]; S.e2[1]=e2[1]; S.e2[2]=e2[2];
  S.e3[0]=e3[0]; S.e3[1]=e3[1]; S.e3[2]=e3[2];

  // Normalize and cache the solar-rotation axis independently of the CME
  // propagation frame.  Parker geometry is tied to solar rotation, not to the
  // CME direction, so these two axes must never be conflated.  A zero or
  // non-finite rotation axis makes the 3-D Parker basis undefined; unlike the
  // legacy safe_normalize() helper, we fail explicitly here instead of silently
  // substituting an arbitrary direction.  Full centralized configuration
  // validation is planned separately, but this local guard is required for a
  // physically meaningful Parker field now.
  {
    const double ax=P_.solar_rotation_axis[0];
    const double ay=P_.solar_rotation_axis[1];
    const double az=P_.solar_rotation_axis[2];
    const double axis_norm=std::sqrt(ax*ax+ay*ay+az*az);
    if (!std::isfinite(axis_norm) || axis_norm<=0.0) {
      throw std::invalid_argument("swcme3d: solar_rotation_axis must be finite and non-zero");
    }
    S.solar_axis_hat[0]=ax/axis_norm;
    S.solar_axis_hat[1]=ay/axis_norm;
    S.solar_axis_hat[2]=az/axis_norm;
    if (!std::isfinite(P_.solar_rotation_rate_rad_s) ||
        P_.solar_rotation_rate_rad_s<0.0) {
      throw std::invalid_argument(
          "swcme3d: solar_rotation_rate_rad_s must be finite and non-negative");
    }
    S.solar_rotation_rate_rad_s=P_.solar_rotation_rate_rad_s;
  }

  // 2) Shared CME/shock-apex kinematics.
  //
  // The legacy 3-D DBM used (1+Gamma*u0*t) with signed u0 directly in the
  // denominator/logarithm.  That is valid only for the fast-CME branch and
  // makes a slow CME move away from Vsw.  It also divided by Gamma when
  // Gamma=0 and then used finite_or() to hide the invalid radius.  Both 1-D
  // and 3-D now call the same sign-aware common solver, which implements
  //   DeltaV(t)=DeltaV0/(1+Gamma*|DeltaV0|*t)
  // and an exact Gamma=0 ballistic limit.  DATA_DRIVEN mode is evaluated by
  // the same common component using monotone PCHIP interpolation.
  const double r0_m=P_.r0_Rs*Rs;
  S.V_sw_ms=P_.V_sw_kms*1e3;
  S.kinematics_mode=P_.kinematics_mode;

  swcme::kinematics::Config kin;
  kin.mode=P_.kinematics_mode;
  kin.r0_m=r0_m;
  kin.V0_m_s=P_.V0_sh_kms*1e3;
  kin.Vsw_m_s=S.V_sw_ms;
  kin.Gamma_m_inv=P_.Gamma_kmInv/1e3;
  kin.extrapolation=P_.data_extrapolation;
  kin.data_time_s=P_.data_time_s;
  kin.data_radius_m.reserve(P_.data_radius_Rs.size());
  for (double radius_Rs : P_.data_radius_Rs) {
    kin.data_radius_m.push_back(radius_Rs*Rs);
  }

  const swcme::kinematics::State apex=swcme::kinematics::evaluate(kin,t_s);
  if (apex.status!=swcme::kinematics::Status::Ok) {
    // prepare_step() predates an explicit model-status return channel.  Until
    // that broader API refactor is completed, report an invalid or out-of-time
    // trajectory as an exception rather than converting it to a plausible
    // radius/speed with finite_or().
    throw std::runtime_error(std::string("swcme3d kinematics: ")+
                             swcme::kinematics::status_name(apex.status));
  }
  S.r_sh_m=apex.radius_m;
  S.V_sh_ms=apex.speed_m_s;
  S.a_m=S.r_sh_m;

  // 3) Self-similar region widths & derived radii
  const double scaleR = S.r_sh_m / AU;
  S.dr_sheath_m = finite_or(P_.sheath_thick_AU_at1AU * scaleR * AU);
  S.dr_me_m     = finite_or(P_.ejecta_thick_AU_at1AU * scaleR * AU);
  S.w_shock_m   = finite_or(P_.edge_smooth_shock_AU_at1AU * scaleR * AU);
  S.w_le_m      = finite_or(P_.edge_smooth_le_AU_at1AU    * scaleR * AU);
  S.w_te_m      = finite_or(P_.edge_smooth_te_AU_at1AU    * scaleR * AU);
  S.r_le_m = S.r_sh_m - S.dr_sheath_m;
  S.r_te_m = S.r_le_m - S.dr_me_m;

  // Smoothstep helpers
  S.inv2w_sh = (S.w_shock_m>0.0)? 0.5/S.w_shock_m : 0.0;
  S.inv2w_le = (S.w_le_m   >0.0)? 0.5/S.w_le_m    : 0.0;
  S.inv2w_te = (S.w_te_m   >0.0)? 0.5/S.w_te_m    : 0.0;

  // 4) Region target speeds
  S.V_sheath_LE_ms = finite_or(P_.V_sheath_LE_factor * S.V_sw_ms, S.V_sw_ms);
  S.V_ME_ms        = finite_or(P_.V_ME_factor        * S.V_sw_ms, S.V_sw_ms);
  S.V_dn_ms        = S.V_sw_ms;

  // 5) Convenience
  S.inv_dr_sheath  = (S.dr_sheath_m>0.0)? 1.0/S.dr_sheath_m : 0.0;
  S.rc_floor       = (P_.sheath_comp_floor>1.0)? P_.sheath_comp_floor : 1.0;

  // 6) Leblanc coefficients in SI (cached)
  {
    const double A=3.3e5, B=4.1e6, C=8.0e7; // cm^-3 coefficients
    const double Rs2=Rs*Rs, Rs4=Rs2*Rs2, Rs6=Rs4*Rs2;
    const double sAU=Rs/AU;
    const double n1_base=A*sAU*sAU + B*std::pow(sAU,4) + C*std::pow(sAU,6);
    const double leb_scale=(n1_base>0.0)? (P_.n1AU_cm3/n1_base):1.0;
    const double K = leb_scale*1e6; // cm^-3 → m^-3
    S.C2 = K*A*Rs2;
    S.C4 = K*B*Rs4;
    S.C6 = K*C*Rs6;
  }

  // 7) Parker constants for this step.
  //
  // k_AU is intentionally latitude-independent: it stores Omega*AU/V_sw.
  // The local sin(theta) is evaluated from geometry in parker_vec_T_fast().
  // Params::sin_theta is retained only as the reference colatitude used to
  // interpret the legacy B1AU_nT total-field normalization.  This preserves
  // existing inputs while removing the physically incorrect global-latitude
  // assumption from the 3-D field itself.
  {
    const double B1AU_T = P_.B1AU_nT*1e-9;
    const double reference_sin_theta = P_.sin_theta;
    if (!std::isfinite(reference_sin_theta) ||
        reference_sin_theta<0.0 || reference_sin_theta>1.0) {
      throw std::invalid_argument(
          "swcme3d: sin_theta reference normalization must lie in [0,1]");
    }

    // A zero solar-wind speed is not a valid Parker-spiral configuration.  The
    // broader configuration layer will eventually reject it explicitly.  For
    // now we keep prepare_step() finite so existing unit-conversion tests can
    // still inspect zero-valued inputs without producing Inf/NaN state.
    S.k_AU = (S.V_sw_ms!=0.0) ? (S.solar_rotation_rate_rad_s*AU/S.V_sw_ms) : 0.0;

    const double reference_pitch = S.k_AU*reference_sin_theta;
    S.Br1AU_T = B1AU_T /
        std::sqrt(1.0 + reference_pitch*reference_pitch);
  }

  // 8) Geometry caches
  if (P_.shape==ShockShape::Ellipsoid){
    S.a_e    = S.r_sh_m;
    S.b_e    = S.a_e * std::max(1e-3, P_.axis_ratio_y);
    S.c_e    = S.a_e * std::max(1e-3, P_.axis_ratio_z);
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
    if (!std::isfinite(P_.half_width_rad) || P_.half_width_rad<=0.0 ||
        P_.half_width_rad>0.5*PI) {
      throw std::invalid_argument(
          "swcme3d: SSE half_width_rad must satisfy 0 < half_width_rad <= pi/2");
    }
    S.sin_half_width = std::sin(P_.half_width_rad);
    S.cos_half_width = std::cos(P_.half_width_rad);
    const double denom = 1.0 + S.sin_half_width;
    S.sse_center_m = S.r_sh_m/denom;
    S.sse_radius_m = S.sse_center_m*S.sin_half_width;
  }

  // 9) Apex shock diagnostic.  A geometric CME front and a physical fast
  // shock are not synonymous; cache both the explicit existence flag and the
  // physical compression for quick time-series diagnostics.
  {
    double u_apex[3]={e1[0],e1[1],e1[2]};
    LocalShockState apex;
    const bool surface_exists=shock_state_direction(S,u_apex,apex);
    S.has_shock=surface_exists && apex.has_shock && apex.solver_converged;
    S.rc=S.has_shock ? apex.compression : 1.0;
  }
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
  double u[3]={ux,uy,uz}; ::safe_normalize(u);

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
      ::safe_normalize(ng);
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
      // a already produces a unit vector up to roundoff; safe_normalize removes
      // the remaining floating-point drift without changing its direction.
      double ng[3]={
          Rdir_m*u[0]-c*S.e1[0],
          Rdir_m*u[1]-c*S.e1[1],
          Rdir_m*u[2]-c*S.e1[2]};
      if (a<=0.0) return false;  // guarded in prepare_step(); defensive only
      ng[0]/=a; ng[1]/=a; ng[2]/=a;
      ::safe_normalize(ng);
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
bool Model::shock_state_direction(const StepState& S, const double u_in[3],
                                  LocalShockState& state) const {
  state=LocalShockState{};

  double u[3]={u_in[0],u_in[1],u_in[2]};
  ::safe_normalize(u);
  double Rdir=0.0,n_hat[3]={0.0,0.0,0.0};
  if (!shape_radius_normal(S,u[0],u[1],u[2],Rdir,n_hat)) {
    // No finite shock surface exists in this direction (for example outside
    // the angular support of the SSE cap).  This is not a nonlinear-solver
    // failure, so solver_converged remains true and compression remains one.
    return false;
  }

  state.surface_exists=true;
  state.Rdir_m=Rdir;
  state.normal[0]=n_hat[0]; state.normal[1]=n_hat[1]; state.normal[2]=n_hat[2];

  // Upstream density is evaluated at the shock surface, never at the caller's
  // sample radius.  Using the query point here made the same physical shock
  // acquire different Mach numbers depending on where the model was sampled.
  const double r=std::max(Rdir,1.05*Rs);
  const double r2=r*r, inv2=1.0/r2, inv4=inv2*inv2, inv6=inv4*inv2;
  const double n_up_m3=finite_or(S.C2*inv2 + S.C4*inv4 + S.C6*inv6,0.0);
  state.upstream_n_m3=n_up_m3;

  double B_up[3]={0.0,0.0,0.0};
  ::parker_vec_T_fast(S,u,r,B_up);

  // Every supported geometry is self-similar with the apex radius.  The
  // surface point on a fixed ray therefore moves radially at
  // (Rdir/Rapex)*Vapex; the physical jump condition uses only the projection
  // of that motion onto the local outward surface normal.
  const double radial_scale=(S.r_sh_m>0.0)? Rdir/S.r_sh_m : 0.0;
  double normal_projection=n_hat[0]*u[0]+n_hat[1]*u[1]+n_hat[2]*u[2];
  const double projection_tol=128.0*std::numeric_limits<double>::epsilon();
  if (normal_projection<0.0 && normal_projection>=-projection_tol)
    normal_projection=0.0;
  state.Vsh_n_m_s=finite_or(S.V_sh_ms*radial_scale*normal_projection,0.0);

  // Current SWCME thermodynamics uses a proton-only thermal pressure
  // p=n_p k_B T_p and rho=m_p n_p.  Composition/temperature generalization is
  // a separate work package; using the existing convention here keeps this
  // correction focused on shock existence and MHD conservation.
  swcme::shock::PrimitiveState upstream;
  upstream.rho_kg_m3=std::max(0.0,n_up_m3)*MP;
  upstream.pressure_Pa=std::max(0.0,n_up_m3)*KB*std::max(0.0,P_.T_K);
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
  state.downstream_n_m3=state.has_shock && state.solver_converged
                         ? state.compression*n_up_m3 : n_up_m3;
  state.mass_residual=jump.mass_residual;
  state.normal_B_residual=jump.normal_B_residual;
  state.electric_residual=jump.electric_residual;
  state.momentum_residual=jump.momentum_residual;
  state.energy_residual=jump.energy_residual;
  state.entropy_ratio=jump.entropy_ratio;
  return true;
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
  if (!observer_m || !point_m || !std::isfinite(radius_m) || radius_m<=0.0 ||
      !std::isfinite(S.V_sw_ms) || S.V_sw_ms<=0.0) {
    if (point_m) point_m[0]=point_m[1]=point_m[2]=0.0;
    return false;
  }

  const double r_obs=norm3(observer_m);
  if (!std::isfinite(r_obs) || r_obs<=0.0) {
    point_m[0]=point_m[1]=point_m[2]=0.0;
    return false;
  }

  const double u_obs[3]={observer_m[0]/r_obs,observer_m[1]/r_obs,
                         observer_m[2]/r_obs};
  const double delta_phi=-S.solar_rotation_rate_rad_s*(radius_m-r_obs)/S.V_sw_ms;
  double u[3]={0.0,0.0,0.0};
  rotate_about_unit_axis(u_obs,S.solar_axis_hat,delta_phi,u);
  ::safe_normalize(u);  // removes only roundoff from the orthogonal rotation

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

// n, V evaluator (allocation-free; vectorization-friendly)
//
// The shock itself is treated as a physical discontinuity: r>=R_sh is the
// upstream state and the limit r->R_sh^- is the Rankine-Hugoniot downstream
// state returned by shock_state_direction().  The old implementation blended
// upstream and sheath across the shock and therefore returned V_sw at the
// immediate downstream boundary, violating mass conservation.  Smooth region
// shaping remains inside the sheath/ejecta where it does not alter the jump.
void Model::evaluate_cartesian_fast(const StepState& S,
                                    const double* x_m,const double* y_m,const double* z_m,
                                    double* n_m3,double* Vx_ms,double* Vy_ms,double* Vz_ms,
                                    std::size_t N) const {
  const double Vup=S.V_sw_ms;
  const double dr_sheath=S.dr_sheath_m, dr_me=S.dr_me_m;

  SWCME_IVDEP
  for (std::size_t i=0;i<N;++i){
    const double x=x_m[i], y=y_m[i], z=z_m[i];
    const double r2=std::max(1e-12,x*x+y*y+z*z);
    const double r=std::sqrt(r2), invr=1.0/r;
    const double u[3]={x*invr,y*invr,z*invr};
    const double inv2=1.0/r2, inv4=inv2*inv2, inv6=inv4*inv2;
    const double n_up=finite_or(S.C2*inv2+S.C4*inv4+S.C6*inv6,1e6);

    LocalShockState shock;
    const bool surface_exists=shock_state_direction(S,u,shock);
    if (!surface_exists || !shock.has_shock || !shock.solver_converged || r>=shock.Rdir_m) {
      n_m3[i]=n_up;
      Vx_ms[i]=finite_or(Vup*u[0],0.0);
      Vy_ms[i]=finite_or(Vup*u[1],0.0);
      Vz_ms[i]=finite_or(Vup*u[2],0.0);
      continue;
    }

    const double r_le=shock.Rdir_m-dr_sheath;
    const double r_te=r_le-dr_me;
    if (r>=r_le && dr_sheath>0.0) {
      const double xi=swcme3d::clamp01((shock.Rdir_m-r)/dr_sheath);
      const double power=std::max(1.0,P_.sheath_ramp_power);
      const double blend=smoothstep01(std::pow(xi,power));

      // Use the exact RH downstream state as the sheath's outer boundary.
      // The inner/leading-edge target remains phenomenological and is kept
      // radial; interpolation of the full velocity vector preserves possible
      // tangential velocity generated by an oblique MHD jump near the shock.
      const double r_le_safe=std::max(r_le,1.05*Rs);
      const double inv2le=1.0/(r_le_safe*r_le_safe);
      const double n_up_le=finite_or(S.C2*inv2le+S.C4*inv2le*inv2le+
                                     S.C6*inv2le*inv2le*inv2le,n_up);
      const double n2=std::max(shock.downstream_n_m3,1.0e-300);
      const double nle=std::max(n_up_le,1.0e-300);
      const double n_sheath=std::exp((1.0-blend)*std::log(n2)+blend*std::log(nle));

      const double Vle_mag=std::max(Vup,P_.V_sheath_LE_factor*Vup);
      const double Vle[3]={Vle_mag*u[0],Vle_mag*u[1],Vle_mag*u[2]};
      const double Vx=(1.0-blend)*shock.downstream.velocity_m_s[0]+blend*Vle[0];
      const double Vy=(1.0-blend)*shock.downstream.velocity_m_s[1]+blend*Vle[1];
      const double Vz=(1.0-blend)*shock.downstream.velocity_m_s[2]+blend*Vle[2];

      n_m3[i]=finite_or(n_sheath,n_up);
      Vx_ms[i]=finite_or(Vx,Vup*u[0]);
      Vy_ms[i]=finite_or(Vy,Vup*u[1]);
      Vz_ms[i]=finite_or(Vz,Vup*u[2]);
      continue;
    }

    if (r>=r_te) {
      const double n_ejecta=std::max(0.0,P_.f_ME)*n_up;
      const double V_ejecta=S.V_ME_ms;
      n_m3[i]=finite_or(n_ejecta,n_up);
      Vx_ms[i]=finite_or(V_ejecta*u[0],0.0);
      Vy_ms[i]=finite_or(V_ejecta*u[1],0.0);
      Vz_ms[i]=finite_or(V_ejecta*u[2],0.0);
      continue;
    }

    n_m3[i]=n_up;
    Vx_ms[i]=finite_or(Vup*u[0],0.0);
    Vy_ms[i]=finite_or(Vup*u[1],0.0);
    Vz_ms[i]=finite_or(Vup*u[2],0.0);
  }
}

// n, V, B evaluator.  The plasma jump at the shock uses the complete MHD
// downstream state; the sheath interior then relaxes toward a phenomenological
// leading-edge target.  This guarantees that B_n, tangential electric field,
// mass flux, momentum flux, and total-energy flux are correct at r->R_sh^-.
void Model::evaluate_cartesian_with_B(const StepState& S,
  const double* x_m,const double* y_m,const double* z_m,
  double* n_m3,double* Vx_ms,double* Vy_ms,double* Vz_ms,
  double* Bx_T,double* By_T,double* Bz_T,
  std::size_t N) const {

  const double Vup=S.V_sw_ms;
  const double dr_sheath=S.dr_sheath_m, dr_me=S.dr_me_m;

  SWCME_IVDEP
  for (std::size_t i=0;i<N;++i){
    const double x=x_m[i], y=y_m[i], z=z_m[i];
    const double r2=std::max(1e-12,x*x+y*y+z*z);
    const double r=std::sqrt(r2), invr=1.0/r;
    const double u[3]={x*invr,y*invr,z*invr};
    const double inv2=1.0/r2, inv4=inv2*inv2, inv6=inv4*inv2;
    const double n_up=finite_or(S.C2*inv2+S.C4*inv4+S.C6*inv6,1e6);
    double B_up[3]; ::parker_vec_T_fast(S,u,r,B_up);

    LocalShockState shock;
    const bool surface_exists=shock_state_direction(S,u,shock);
    if (!surface_exists || !shock.has_shock || !shock.solver_converged || r>=shock.Rdir_m) {
      n_m3[i]=n_up;
      Vx_ms[i]=finite_or(Vup*u[0],0.0);
      Vy_ms[i]=finite_or(Vup*u[1],0.0);
      Vz_ms[i]=finite_or(Vup*u[2],0.0);
      Bx_T[i]=finite_or(B_up[0],0.0);
      By_T[i]=finite_or(B_up[1],0.0);
      Bz_T[i]=finite_or(B_up[2],0.0);
      continue;
    }

    const double r_le=shock.Rdir_m-dr_sheath;
    const double r_te=r_le-dr_me;
    if (r>=r_le && dr_sheath>0.0) {
      const double xi=swcme3d::clamp01((shock.Rdir_m-r)/dr_sheath);
      const double power=std::max(1.0,P_.sheath_ramp_power);
      const double blend=smoothstep01(std::pow(xi,power));

      const double r_le_safe=std::max(r_le,1.05*Rs);
      const double inv2le=1.0/(r_le_safe*r_le_safe);
      const double n_up_le=finite_or(S.C2*inv2le+S.C4*inv2le*inv2le+
                                     S.C6*inv2le*inv2le*inv2le,n_up);
      const double n2=std::max(shock.downstream_n_m3,1.0e-300);
      const double nle=std::max(n_up_le,1.0e-300);
      const double n_sheath=std::exp((1.0-blend)*std::log(n2)+blend*std::log(nle));

      const double Vle_mag=std::max(Vup,P_.V_sheath_LE_factor*Vup);
      const double Vle[3]={Vle_mag*u[0],Vle_mag*u[1],Vle_mag*u[2]};
      const double Vx=(1.0-blend)*shock.downstream.velocity_m_s[0]+blend*Vle[0];
      const double Vy=(1.0-blend)*shock.downstream.velocity_m_s[1]+blend*Vle[1];
      const double Vz=(1.0-blend)*shock.downstream.velocity_m_s[2]+blend*Vle[2];

      // At the inner sheath edge return to the local Parker field.  This is a
      // phenomenological relaxation, but the outer boundary is now the exact
      // MHD B2 state rather than the old ad hoc B_t*=r_c amplification.
      double B_le[3]; ::parker_vec_T_fast(S,u,r_le_safe,B_le);
      const double Bx=(1.0-blend)*shock.downstream.magnetic_T[0]+blend*B_le[0];
      const double By=(1.0-blend)*shock.downstream.magnetic_T[1]+blend*B_le[1];
      const double Bz=(1.0-blend)*shock.downstream.magnetic_T[2]+blend*B_le[2];

      n_m3[i]=finite_or(n_sheath,n_up);
      Vx_ms[i]=finite_or(Vx,Vup*u[0]);
      Vy_ms[i]=finite_or(Vy,Vup*u[1]);
      Vz_ms[i]=finite_or(Vz,Vup*u[2]);
      Bx_T[i]=finite_or(Bx,B_up[0]);
      By_T[i]=finite_or(By,B_up[1]);
      Bz_T[i]=finite_or(Bz,B_up[2]);
      continue;
    }

    if (r>=r_te) {
      const double n_ejecta=std::max(0.0,P_.f_ME)*n_up;
      const double V_ejecta=S.V_ME_ms;
      n_m3[i]=finite_or(n_ejecta,n_up);
      Vx_ms[i]=finite_or(V_ejecta*u[0],0.0);
      Vy_ms[i]=finite_or(V_ejecta*u[1],0.0);
      Vz_ms[i]=finite_or(V_ejecta*u[2],0.0);
      Bx_T[i]=finite_or(B_up[0],0.0);
      By_T[i]=finite_or(B_up[1],0.0);
      Bz_T[i]=finite_or(B_up[2],0.0);
      continue;
    }

    n_m3[i]=n_up;
    Vx_ms[i]=finite_or(Vup*u[0],0.0);
    Vy_ms[i]=finite_or(Vup*u[1],0.0);
    Vz_ms[i]=finite_or(Vup*u[2],0.0);
    Bx_T[i]=finite_or(B_up[0],0.0);
    By_T[i]=finite_or(B_up[1],0.0);
    Bz_T[i]=finite_or(B_up[2],0.0);
  }
}

void Model::evaluate_cartesian_with_B_div(const StepState& S,
  const double* x_m,const double* y_m,const double* z_m,
  double* n_m3,double* Vx_ms,double* Vy_ms,double* Vz_ms,
  double* Bx_T,double* By_T,double* Bz_T,double* divVsw,
  std::size_t N,double dr_frac) const {
  evaluate_cartesian_with_B(S,x_m,y_m,z_m,n_m3,Vx_ms,Vy_ms,Vz_ms,Bx_T,By_T,Bz_T,N);
  compute_divV_radial(S,x_m,y_m,z_m,divVsw,N,dr_frac);
}

// Robust radial divergence via (1/r^2) d(r^2 V_r)/dr
void Model::compute_divV_radial(const StepState& S,
  const double* x_m,const double* y_m,const double* z_m,
  double* divV,std::size_t N,double dr_frac) const {
  const double rmin=1.05*Rs, dr_min=1.0e-4*AU;

  SWCME_IVDEP
  for (std::size_t i=0;i<N;++i){
    const double x=x_m[i], y=y_m[i], z=z_m[i];
    const double r2=std::max(1e-12, x*x+y*y+z*z);
    const double r =std::sqrt(r2);
    const double invr=1.0/r;
    const double u[3]={x*invr,y*invr,z*invr};
    const double dr=std::max(dr_min, dr_frac*r);
    const double rp=std::max(rmin,r+dr), rm=std::max(rmin,r-dr);
    const double denom_r=(rp>rm)? (rp-rm): std::max(dr_min,std::abs(dr));

    double xp=rp*u[0], yp=rp*u[1], zp=rp*u[2];
    double xm=rm*u[0], ym=rm*u[1], zm=rm*u[2];
    double n,Vxp,Vyp,Vzp,Vxm,Vym,Vzm;
    evaluate_cartesian_fast(S,&xp,&yp,&zp,&n,&Vxp,&Vyp,&Vzp,1);
    evaluate_cartesian_fast(S,&xm,&ym,&zm,&n,&Vxm,&Vym,&Vzm,1);

    const double Vrp=Vxp*u[0]+Vyp*u[1]+Vzp*u[2];
    const double Vrm=Vxm*u[0]+Vym*u[1]+Vzm*u[2];

    const double num=((rp*rp)*Vrp - (rm*rm)*Vrm)/denom_r;
    const double div_val=num/(r*r);
    divV[i]=finite_or(div_val,0.0);
  }
}

bool Model::diagnose_direction(const StepState& S,const double u[3],
  double& Rdir_m,double n_hat[3],double& rc_loc,double& Vsh_n) const {
  const bool has_surface=shape_radius_normal(S,u[0],u[1],u[2],Rdir_m,n_hat);
  if (!has_surface) {
    // Explicit no-surface diagnostics make finite-width behavior visible to
    // callers instead of reporting a plausible but fabricated flank radius.
    rc_loc=1.0;
    Vsh_n=0.0;
    return false;
  }
  double th=0.0;
  local_oblique_rc(S,u,n_hat,Rdir_m,Rdir_m,rc_loc,Vsh_n,th);
  return true;
}

// Build lat–lon mesh on [0,thetaMax]×[0,2π]
ShockMesh Model::build_shock_mesh(const StepState& S,std::size_t nTheta,std::size_t nPhi) const {
  ShockMesh M; if (nTheta<3) nTheta=3; if (nPhi<3) nPhi=3;
  const double PI=swcme3d::PI;
  const double thetaMax=(P_.shape==ShockShape::SSE)? P_.half_width_rad : PI;

  for (std::size_t it=0; it<=nTheta; ++it){
    const double t=thetaMax*(double(it)/double(nTheta));
    const double ct=std::cos(t), st=std::sin(t);
    for (std::size_t ip=0; ip<=nPhi; ++ip){
      const double p=2.0*PI*(double(ip)/double(nPhi)), cp=std::cos(p), sp=std::sin(p);
      double u_loc[3]={ct, st*cp, st*sp};
      double u[3]={ u_loc[0]*S.e1[0]+u_loc[1]*S.e2[0]+u_loc[2]*S.e3[0],
                    u_loc[0]*S.e1[1]+u_loc[1]*S.e2[1]+u_loc[2]*S.e3[1],
                    u_loc[0]*S.e1[2]+u_loc[1]*S.e2[2]+u_loc[2]*S.e3[2] };
      ::safe_normalize(u);

      double Rdir=0.0,n_hat[3]={0,0,0};
      const bool has_surface=shape_radius_normal(S,u[0],u[1],u[2],Rdir,n_hat);
      // build_shock_mesh samples only the mathematically supported angular
      // interval.  Failure here therefore signals an internal geometry error,
      // not an expected outside-cap query, and should never be silently filled
      // with a zero-radius vertex.
      if (!has_surface) {
        throw std::runtime_error("swcme3d: shock mesh requested a direction outside the supported surface");
      }

      M.x.push_back(finite_or(Rdir*u[0],0.0));
      M.y.push_back(finite_or(Rdir*u[1],0.0));
      M.z.push_back(finite_or(Rdir*u[2],0.0));
      M.n_hat_x.push_back(finite_or(n_hat[0],0.0));
      M.n_hat_y.push_back(finite_or(n_hat[1],0.0));
      M.n_hat_z.push_back(finite_or(n_hat[2],1.0));

      double rc_loc=1.0,Vsh_n=0.0,th=0.0;
      local_oblique_rc(S,u,n_hat,Rdir,Rdir,rc_loc,Vsh_n,th);
      M.rc.push_back(finite_or(rc_loc,1.0));
      M.Vsh_n.push_back(finite_or(Vsh_n,0.0));
    }
  }
  const std::size_t NvPhi=nPhi+1;
  auto idx=[&](std::size_t it,std::size_t ip){ return int(it*NvPhi+ip); };
  for (std::size_t it=0; it<nTheta; ++it){
    for (std::size_t ip=0; ip<nPhi; ++ip){
      int i00=idx(it,ip), i01=idx(it,ip+1), i10=idx(it+1,ip), i11=idx(it+1,ip+1);
      M.tri_i.push_back(i00+1); M.tri_j.push_back(i10+1); M.tri_k.push_back(i11+1);
      M.tri_i.push_back(i00+1); M.tri_j.push_back(i11+1); M.tri_k.push_back(i01+1);
    }
  }
  return M;
}

void Model::compute_triangle_metrics(const ShockMesh& M, TriMetrics& T) const {
  const std::size_t Ne=M.tri_i.size();
  T.area.assign(Ne,0.0); T.nx.assign(Ne,0.0); T.ny.assign(Ne,0.0); T.nz.assign(Ne,1.0);
  T.cx.assign(Ne,0.0); T.cy.assign(Ne,0.0); T.cz.assign(Ne,0.0);
  T.rc_mean.assign(Ne,1.0); T.Vsh_n_mean.assign(Ne,0.0);

  for (std::size_t e=0;e<Ne;++e){
    int ia=M.tri_i[e]-1, ib=M.tri_j[e]-1, ic=M.tri_k[e]-1;
    double Ax=M.x[ia], Ay=M.y[ia], Az=M.z[ia];
    double Bx=M.x[ib], By=M.y[ib], Bz=M.z[ib];
    double Cx=M.x[ic], Cy=M.y[ic], Cz=M.z[ic];

    double ABx=Bx-Ax, ABy=By-Ay, ABz=Bz-Az;
    double ACx=Cx-Ax, ACy=Cy-Ay, ACz=Cz-Az;

    double nx=ABy*ACz - ABz*ACy;
    double ny=ABz*ACx - ABx*ACz;
    double nz=ABx*ACy - ABy*ACx;
    double twiceA=std::sqrt(std::max(0.0, nx*nx+ny*ny+nz*nz));
    double area=0.5*twiceA;
    if (twiceA>0){ nx/=twiceA; ny/=twiceA; nz/=twiceA; }

    T.area[e]=finite_or(area,0.0);
    T.nx[e]=finite_or(nx,0.0); T.ny[e]=finite_or(ny,0.0); T.nz[e]=finite_or(nz,1.0);
    T.cx[e]=finite_or((Ax+Bx+Cx)/3.0,0.0);
    T.cy[e]=finite_or((Ay+By+Cy)/3.0,0.0);
    T.cz[e]=finite_or((Az+Bz+Cz)/3.0,0.0);
    T.rc_mean[e]=finite_or((M.rc[ia]+M.rc[ib]+M.rc[ic])/3.0,1.0);
    T.Vsh_n_mean[e]=finite_or((M.Vsh_n[ia]+M.Vsh_n[ib]+M.Vsh_n[ic])/3.0,0.0);
  }
}

// --- Tecplot helpers ---------------------------------------------------------
static inline void dump_array_block(std::FILE* fp, const std::vector<double>& a){
  int cnt=0; for(double v:a){ std::fprintf(fp,"%.9e ",finite_or(v)); if(++cnt==8){std::fprintf(fp,"\n"); cnt=0;} } if(cnt) std::fprintf(fp,"\n");
}
static inline void dump_zeros_block(std::FILE* fp, std::size_t count){
  int cnt=0; for(std::size_t e=0;e<count;++e){ std::fprintf(fp,"0 "); if(++cnt==8){std::fprintf(fp,"\n"); cnt=0;} } if(cnt) std::fprintf(fp,"\n");
}

// Surface-only: cell metrics + nodal rc/Vsh_n
bool Model::write_shock_surface_center_metrics_tecplot(
  const ShockMesh& M, const TriMetrics& T_in, const char* path) const {

  TriMetrics T=T_in;
  const std::size_t Nv=M.x.size(), Ne=M.tri_i.size();
  if (T.area.size()!=Ne){ compute_triangle_metrics(M,T); }

  std::FILE* fp=std::fopen(path,"w"); if(!fp) return false;

  std::fprintf(fp,"TITLE=\"Shock surface (cell metrics + nodal rc)\"\n");
  std::fprintf(fp,"VARIABLES=\"X\",\"Y\",\"Z\",\"n\",\"Vx\",\"Vy\",\"Vz\",\"Bx\",\"By\",\"Bz\",\"divVsw\","
                  "\"rc\",\"Vsh_n\",\"nx\",\"ny\",\"nz\",\"area\",\"rc_mean\",\"Vsh_n_mean\",\"tnx\",\"tny\",\"tnz\",\"cx\",\"cy\",\"cz\"\n");
  std::fprintf(fp,"ZONE T=\"surface_cells\", N=%zu, E=%zu, ZONETYPE=FETRIANGLE, DATAPACKING=BLOCK,\n", Nv, Ne);
  std::fprintf(fp,"VARLOCATION=([1-3,12-13]=NODAL, [4-11,14-25]=CELLCENTERED)\n");

  dump_array_block(fp,M.x); dump_array_block(fp,M.y); dump_array_block(fp,M.z);
  dump_zeros_block(fp,Ne); dump_zeros_block(fp,Ne); dump_zeros_block(fp,Ne);
  dump_zeros_block(fp,Ne);
  dump_zeros_block(fp,Ne); dump_zeros_block(fp,Ne); dump_zeros_block(fp,Ne);
  dump_zeros_block(fp,Ne);
  dump_array_block(fp,M.rc); dump_array_block(fp,M.Vsh_n);
  dump_array_block(fp,T.nx); dump_array_block(fp,T.ny); dump_array_block(fp,T.nz);
  dump_array_block(fp,T.area); dump_array_block(fp,T.rc_mean); dump_array_block(fp,T.Vsh_n_mean);
  dump_zeros_block(fp,Ne); dump_zeros_block(fp,Ne); dump_zeros_block(fp,Ne);
  dump_array_block(fp,T.cx); dump_array_block(fp,T.cy); dump_array_block(fp,T.cz);

  for (std::size_t e=0;e<Ne;++e)
    std::fprintf(fp,"%d %d %d\n", M.tri_i[e], M.tri_j[e], M.tri_k[e]);

  std::fclose(fp);
  return true;
}

// Default apex-aligned volume box
BoxSpec Model::default_apex_box(const StepState& S,double half_AU,int N) const {
  BoxSpec B; const double h=half_AU*AU;
  B.hx=h; B.hy=h; B.hz=h;
  const double shift=0.4*h; // move box outward along e1 so shock cuts through
  B.cx=S.a_m*S.e1[0]+shift*S.e1[0];
  B.cy=S.a_m*S.e1[1]+shift*S.e1[1];
  B.cz=S.a_m*S.e1[2]+shift*S.e1[2];
  B.Ni=N; B.Nj=N; B.Nk=N; return B;
}

// Bundle writer: surface_cells + surface_nodal + volume_box + minX face
bool Model::write_tecplot_dataset_bundle(const ShockMesh& M,const TriMetrics& T_in,
                                         const StepState& S,const BoxSpec& B,
                                         const char* path) const {
  TriMetrics T=T_in;
  const std::size_t Nv=M.x.size(), Ne=M.tri_i.size();
  if (T.area.size()!=Ne){ compute_triangle_metrics(M,T); }

  std::FILE* fp=std::fopen(path,"w"); if(!fp) return false;
  auto p=[&](const char* fmt, auto... args){ std::fprintf(fp,fmt,args...); };

  p("TITLE = \"SW+CME dataset\"\n");
  p("VARIABLES = "
    "\"X\",\"Y\",\"Z\","
    "\"n\",\"Vx\",\"Vy\",\"Vz\","
    "\"Bx\",\"By\",\"Bz\",\"divVsw\","
    "\"rc\",\"Vsh_n\","
    "\"nx\",\"ny\",\"nz\",\"area\",\"rc_mean\",\"Vsh_n_mean\","
    "\"tnx\",\"tny\",\"tnz\",\"cx\",\"cy\",\"cz\"\n");

  // Zone 1: surface_cells (FETRIANGLE, BLOCK)
  p("ZONE T=\"surface_cells\", N=%zu, E=%zu, ZONETYPE=FETRIANGLE, DATAPACKING=BLOCK,\n", Nv, Ne);
  p("VARLOCATION=([1-3,12-13]=NODAL, [4-11,14-25]=CELLCENTERED)\n");

  dump_array_block(fp,M.x); dump_array_block(fp,M.y); dump_array_block(fp,M.z);
  dump_zeros_block(fp,Ne); dump_zeros_block(fp,Ne); dump_zeros_block(fp,Ne);
  dump_zeros_block(fp,Ne);
  dump_zeros_block(fp,Ne); dump_zeros_block(fp,Ne); dump_zeros_block(fp,Ne);
  dump_zeros_block(fp,Ne);
  dump_array_block(fp,M.rc); dump_array_block(fp,M.Vsh_n);
  dump_array_block(fp,T.nx); dump_array_block(fp,T.ny); dump_array_block(fp,T.nz);
  dump_array_block(fp,T.area); dump_array_block(fp,T.rc_mean); dump_array_block(fp,T.Vsh_n_mean);
  dump_zeros_block(fp,Ne); dump_zeros_block(fp,Ne); dump_zeros_block(fp,Ne);
  dump_array_block(fp,T.cx); dump_array_block(fp,T.cy); dump_array_block(fp,T.cz);

  for (std::size_t e=0;e<Ne;++e)
    p("%d %d %d\n", M.tri_i[e], M.tri_j[e], M.tri_k[e]);

  // Zone 2: surface_nodal (FEPOINT)
  p("ZONE T=\"surface_nodal\", N=%zu, E=%zu, F=FEPOINT, ET=TRIANGLE\n", Nv, Ne);
  for (std::size_t i=0;i<Nv;++i){
    p("%.9e %.9e %.9e "
      "%.9e %.9e %.9e %.9e "
      "%.9e %.9e %.9e %.9e "
      "%.9e %.9e "
      "%.9e %.9e %.9e "
      "%.9e %.9e %.9e "
      "%.9e %.9e %.9e "
      "%.9e %.9e %.9e\n",
      finite_or(M.x[i]), finite_or(M.y[i]), finite_or(M.z[i]),
      0.0,0.0,0.0,0.0,
      0.0,0.0,0.0,0.0,
      finite_or(M.rc[i],1.0), finite_or(M.Vsh_n[i],0.0),
      finite_or(M.n_hat_x[i],0.0), finite_or(M.n_hat_y[i],0.0), finite_or(M.n_hat_z[i],1.0),
      0.0,0.0,0.0,
      0.0,0.0,0.0,
      0.0,0.0,0.0
    );
  }
  for (std::size_t e=0;e<Ne;++e)
    p("%d %d %d\n", M.tri_i[e], M.tri_j[e], M.tri_k[e]);

  // Zone 3: volume_box (structured POINT)
  p("ZONE T=\"volume_box\", I=%d, J=%d, K=%d, DATAPACKING=POINT\n", B.Ni,B.Nj,B.Nk);
  for (int kk=0; kk<B.Nk; ++kk){
    const double zk = B.cz + (-B.hz + (2.0*B.hz) * (kk / double(std::max(1,B.Nk-1))));
    for (int jj=0; jj<B.Nj; ++jj){
      const double yj = B.cy + (-B.hy + (2.0*B.hy) * (jj / double(std::max(1,B.Nj-1))));
      for (int ii=0; ii<B.Ni; ++ii){
        const double xi = B.cx + (-B.hx + (2.0*B.hx) * (ii / double(std::max(1,B.Ni-1))));
        double n,Vx,Vy,Vz,Bx,By,Bz,div;
        evaluate_cartesian_with_B(S,&xi,&yj,&zk,&n,&Vx,&Vy,&Vz,&Bx,&By,&Bz,1);
        compute_divV_radial(S,&xi,&yj,&zk,&div,1,1e-3);
        p("%.9e %.9e %.9e "
          "%.9e %.9e %.9e %.9e "
          "%.9e %.9e %.9e %.9e "
          "%.9e %.9e "
          "%.9e %.9e %.9e "
          "%.9e %.9e %.9e "
          "%.9e %.9e %.9e "
          "%.9e %.9e %.9e\n",
          finite_or(xi), finite_or(yj), finite_or(zk),
          finite_or(n), finite_or(Vx), finite_or(Vy), finite_or(Vz),
          finite_or(Bx), finite_or(By), finite_or(Bz), finite_or(div,0.0),
          0.0,0.0, 0.0,0.0,0.0, 0.0,0.0,0.0, 0.0,0.0,0.0, 0.0,0.0,0.0
        );
      }
    }
  }

  // Zone 4: min-X face (structured 2-D POINT)
  const int I=std::max(2,B.Nj), J=std::max(2,B.Nk);
  const double x0=B.cx-B.hx;
  const double y0=B.cy-B.hy, y1=B.cy+B.hy;
  const double z0=B.cz-B.hz, z1=B.cz+B.hz;
  p("ZONE T=\"box_face_minX\", I=%d, J=%d, DATAPACKING=POINT\n", I, J);
  for (int j=0;j<J;++j){
    const double tz=(J==1)?0.0: double(j)/double(J-1);
    const double z=z0+(z1-z0)*tz;
    for (int i=0;i<I;++i){
      const double ty=(I==1)?0.0: double(i)/double(I-1);
      const double y=y0+(y1-y0)*ty;
      const double x=x0;
      double n,Vx,Vy,Vz,Bx,By,Bz,div;
      evaluate_cartesian_with_B(S,&x,&y,&z,&n,&Vx,&Vy,&Vz,&Bx,&By,&Bz,1);
      compute_divV_radial(S,&x,&y,&z,&div,1,1e-3);
      p("%.9e %.9e %.9e "
        "%.9e %.9e %.9e %.9e "
        "%.9e %.9e %.9e %.9e "
        "%.9e %.9e "
        "%.9e %.9e %.9e "
        "%.9e %.9e %.9e "
        "%.9e %.9e %.9e "
        "%.9e %.9e %.9e\n",
        finite_or(x), finite_or(y), finite_or(z),
        finite_or(n), finite_or(Vx), finite_or(Vy), finite_or(Vz),
        finite_or(Bx), finite_or(By), finite_or(Bz), finite_or(div,0.0),
        0.0,0.0, 0.0,0.0,0.0, 0.0,0.0,0.0, 0.0,0.0,0.0, 0.0,0.0,0.0
      );
    }
  }

  std::fclose(fp);
  return true;
}

// Standalone 2-D face writer (min-X plane)
bool Model::write_box_face_minX_tecplot_structured(const StepState& S,
                                                   const BoxSpec& B,
                                                   const char* path) const {
  const int I=std::max(2,B.Nj), J=std::max(2,B.Nk);
  const double x0=B.cx-B.hx;
  const double y0=B.cy-B.hy, y1=B.cy+B.hy;
  const double z0=B.cz-B.hz, z1=B.cz+B.hz;

  std::FILE* fp=std::fopen(path,"w"); if(!fp) return false;
  auto p=[&](const char* fmt, auto... args){ std::fprintf(fp,fmt,args...); };

  p("TITLE = \"Box face (minX)\"\n");
  p("VARIABLES = "
    "\"X\",\"Y\",\"Z\","
    "\"n\",\"Vx\",\"Vy\",\"Vz\","
    "\"Bx\",\"By\",\"Bz\",\"divVsw\","
    "\"rc\",\"Vsh_n\","
    "\"nx\",\"ny\",\"nz\",\"area\",\"rc_mean\",\"Vsh_n_mean\","
    "\"tnx\",\"tny\",\"tnz\",\"cx\",\"cy\",\"cz\"\n");
  p("ZONE T=\"box_face_minX\", I=%d, J=%d, DATAPACKING=POINT\n", I, J);

  for (int j=0;j<J;++j){
    const double tz=(J==1)?0.0: double(j)/double(J-1);
    const double z=z0+(z1-z0)*tz;
    for (int i=0;i<I;++i){
      const double ty=(I==1)?0.0: double(i)/double(I-1);
      const double y=y0+(y1-y0)*ty;
      const double x=x0;
      double n,Vx,Vy,Vz,Bx,By,Bz,div;
      evaluate_cartesian_with_B(S,&x,&y,&z,&n,&Vx,&Vy,&Vz,&Bx,&By,&Bz,1);
      compute_divV_radial(S,&x,&y,&z,&div,1,1e-3);
      p("%.9e %.9e %.9e "
        "%.9e %.9e %.9e %.9e "
        "%.9e %.9e %.9e %.9e "
        "%.9e %.9e "
        "%.9e %.9e %.9e "
        "%.9e %.9e %.9e "
        "%.9e %.9e %.9e "
        "%.9e %.9e %.9e\n",
        finite_or(x), finite_or(y), finite_or(z),
        finite_or(n), finite_or(Vx), finite_or(Vy), finite_or(Vz),
        finite_or(Bx), finite_or(By), finite_or(Bz), finite_or(div,0.0),
        0.0,0.0, 0.0,0.0,0.0, 0.0,0.0,0.0, 0.0,0.0,0.0, 0.0,0.0,0.0
      );
    }
  }
  std::fclose(fp);
  return true;
}

} // namespace swcme3d
