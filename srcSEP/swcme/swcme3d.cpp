// swcme3d.cpp
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
//       field is amplified smoothly within the sheath according to the local
//       compression proxy.
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
//    In spherical coordinates (r,θ,φ) with θ measured from rotation axis
//    (ẑ), and assuming azimuthal symmetry of the solar rotation, we use a
//    purely radial Br and toroidal Bφ component:
//       k_AU = Ω⊙ r sinθ / V_sw,  (dimensionless pitch parameter at r)
//       Br(r) = Br(1 AU) (1 AU / r)^2
//       Bφ(r) = -Br(r) * k_AU
//    A desired |B|(1 AU) = B1AU_nT is enforced by setting Br(1 AU) so that
//    sqrt(Br^2 + Bφ^2) = B1AU_nT (converted to Tesla). We construct Cartesian
//    vectors using local spherical radial and azimuthal directions and the
//    given CME propagation latitude via the parameter sin_theta.
//
// C) CME apex kinematics:
//    Drag-Based Model (DBM), Vršnak et al. (2013), Sol. Phys. 285:295.
//       d(V - Vsw)/dt = - Γ (V - Vsw) |V - Vsw|
//    For constant Γ, the closed form solution yields the apex distance and
//    speed as functions of time t:
//       u0 = V0 - Vsw
//       u(t) = u0 / (1 + Γ u0 t)
//       r(t) = r0 + Vsw t + (1/Γ) ln(1 + Γ u0 t)
//    Inputs: r0 = r0_Rs·Rs, V0 (km/s), Vsw (km/s), Γ (km^-1).
//
// D) Shock shape & direction-dependent speed:
//    Shape options: Sphere, Ellipsoid (axis ratios), Cone (SSE-like).
//    For cones we allow flank slowdown via exponent on cos(θ) to mimic slower
//    expansion off apex. The apex kinematics set an overall scale; the signed
//    normal component of shock speed Vsh_n is projected along the local
//    outward normal.
//
// E) Local compression proxy rc(u, n̂):
//    We estimate oblique fast-mode Mach number M_f,n from upstream cs and vA
//    (Alfvén speed), using the local Parker B direction (θBn angle) and normal
//    component of shock speed Vsh_n relative to solar wind Vsw. Then a
//    hydrodynamic strong-shock formula (Edmiston & Kennel 1984, JPP 32:429;
//    Priest 2014) maps M_f,n to a proxy compression ratio rc, saturated at 4.
//    This is used to smoothly amplify density across the sheath and to
//    increase tangential B.
//
// F) Sheath / ejecta blending (C^1 smoothing):
//    Three radial transitions along a given direction u:
//      1) Upstream → Sheath ramp at shock (width w_shock_m).
//      2) Sheath → Ejecta ramp at the leading edge r_le = r_sh - dr_sheath.
//      3) Ejecta → Downstream ambient at the trailing edge r_te = r_le - dr_me.
//    Each edge uses a smoothstep s(x)=x^2(3-2x) with its own width. Inside the
//    sheath, density transitions toward rc·n_up with a user power p on the
//    (1-xi) factor to allow sharper/softer fronts. Tangential magnetic field
//    is gradually amplified between 1 and rc across the shock and sheath.
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
//
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

//
//    Zones:
//      • surface_cells  : FETRIANGLE/BLOCK, with cell-centered metrics
//                         (area, rc_mean, Vsh_n_mean, nx,ny,nz, cx,cy,cz)
//                         and nodal rc, Vsh_n at vertices. This zone appears
//                         first so Tecplot shows non-zero cell metrics by default.
//      • surface_nodal  : FEPOINT, carries nodal rc, Vsh_n, nodal normals.
//                         (Cell-centered fields are zero here by design.)
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

#include "swcme3d.hpp"

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <vector>
#include <string>
#include <algorithm>
#include <limits>

// If true, the Tecplot bundle writer appends a 2D box face zone for convenience
static constexpr bool ADD_FACE_ZONE_IN_BUNDLE = true;

namespace swcme3d {
  // Fundamental length scales (SI)
  const double AU = 1.495978707e11; // Astronomical Unit [m]
  const double Rs = 6.957e8;        // Solar radius [m]

  // Mathematical constant (external linkage so header/other TUs can reference)
  const double PI = 3.141592653589793;

  // Small helper clamp in namespace for re-use
  inline double clamp01(double v){ return (v<0.0)?0.0:(v>1.0?1.0:v); }
}

// ─────────────────────────────────────────────────────────────────────────────
// Sanity helpers: NaN/Inf protection and safe normalization
// ─────────────────────────────────────────────────────────────────────────────

static inline double finite_or(double v, double fallback=0.0){
  // Returns v if it is finite; otherwise returns fallback (default: 0).
  // We route *all* API-visible numeric values through this to immunize I/O.
  return std::isfinite(v)? v : fallback;
}

static inline void safe_normalize(double v[3]){
  // Normalizes a 3-vector; if zero-length, set to +x̂ to avoid NaNs.
  const double m = std::sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
  if (m>0){ v[0]/=m; v[1]/=m; v[2]/=m; } else { v[0]=1; v[1]=0; v[2]=0; }
}

// ─────────────────────────────────────────────────────────────────────────────
// Physical constants
// ─────────────────────────────────────────────────────────────────────────────

static const double MU0 = 4.0e-7*swcme3d::PI; // μ0 (= 4π×10^-7 H/m). Not constexpr (uses swcme3d::PI)
static constexpr double MP  = 1.67262192369e-27; // proton mass [kg]
static constexpr double KB  = 1.380649e-23;      // Boltzmann [J/K]
static constexpr double OMEGA_SUN = 2.86533e-6;  // solar rotation [rad/s] (~25.4 d)

// ─────────────────────────────────────────────────────────────────────────────
// Smoothstep utility: s(x) = 3x^2 - 2x^3 with s(0)=0, s(1)=1, s'(0)=s'(1)=0
// ─────────────────────────────────────────────────────────────────────────────
static inline double smoothstep01(double x){
  if (x<=0) return 0; if (x>=1) return 1; return x*x*(3-2*x);
}

// ─────────────────────────────────────────────────────────────────────────────
// Parker spiral (Tesla) at position r along unit radial direction u
// ─────────────────────────────────────────────────────────────────────────────
// Notes:
//  • k_AU here is constructed with sinθ ≈ sin_theta provided via Params.
//  • Br is scaled to satisfy |B|(1 AU) = B1AU_nT.
//  • We build a local azimuthal unit vector φ̂ perpendicular to u and ẑ.
// ─────────────────────────────────────────────────────────────────────────────
static inline void parker_vec_T(const swcme3d::Params& P,
                                const double u[3], double r_m,
                                double Vsw_ms, double B_out[3]){
  using namespace swcme3d;
  const double r_AU = r_m / AU;
  const double k_AU = (OMEGA_SUN * AU * P.sin_theta) / Vsw_ms;

  // Enforce |B|(1 AU) = B1AU_nT by choosing Br(1 AU) accordingly
  const double B1AU_T = P.B1AU_nT*1e-9;
  const double Br1AU  = B1AU_T / std::sqrt(1.0 + k_AU*k_AU);

  // Radial scaling; φ-field increases with r via k_AU r/AU
  const double Br   = Br1AU / (r_AU*r_AU);
  const double Bphi = - Br * k_AU * r_AU;

  // Build φ̂ from (ẑ×u)×u
  const double zhat[3]={0,0,1};
  double zxu[3] = { zhat[1]*u[2]-zhat[2]*u[1],
                    zhat[2]*u[0]-zhat[0]*u[2],
                    zhat[0]*u[1]-zhat[1]*u[0] };
  double ph[3] = { zxu[1]*u[2]-zxu[2]*u[1],
                   zxu[2]*u[0]-zxu[0]*u[2],
                   zxu[0]*u[1]-zxu[1]*u[0] };
  safe_normalize(ph);

  // Compose vector field in Cartesian coords
  B_out[0]=Br*u[0]+Bphi*ph[0];
  B_out[1]=Br*u[1]+Bphi*ph[1];
  B_out[2]=Br*u[2]+Bphi*ph[2];
}

// ─────────────────────────────────────────────────────────────────────────────
// DBM apex solution (closed form) — Vršnak et al. 2013
// ─────────────────────────────────────────────────────────────────────────────
static inline void dbm_position_speed(double r0_m, double V0_ms, double Vsw_ms,
                                      double Gamma_mInv, double t_s,
                                      double& r_m, double& V_ms){
  // Relative speed u(t) = (V-Vsw) decays as 1/(1+Γ u0 t). Position has log term.
  const double u0 = V0_ms - Vsw_ms;
  const double den= 1.0 + Gamma_mInv*u0*t_s;
  const double u  = (den!=0.0)? (u0/den) : 0.0;
  r_m = r0_m + Vsw_ms*t_s + std::log(std::max(den,1e-30))/Gamma_mInv;
  V_ms= Vsw_ms + u;
}

// ============================================================================
//                                  Model
// ============================================================================

namespace swcme3d {

Model::Model(const Params& P): P_(P) {}

// ----------------------------------------------------------------------------
// prepare_step(t)
//   Builds all time-dependent state from parameters at time t [s].
//   Sets shock apex radius & speed, sheath/ejecta widths and smoothing,
//   and computes an apex rc diagnostic for information.
// ----------------------------------------------------------------------------
StepState Model::prepare_step(double t_s) const {
  StepState S{};

  // Build an orthonormal apex-aligned frame {e1,e2,e3} where e1 ≡ CME axis.
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

  // DBM apex kinematics (r_sh, V_sh)
  const double r0_m   = P_.r0_Rs * Rs;
  const double V0_ms  = P_.V0_sh_kms*1e3;
  const double Vsw_ms = P_.V_sw_kms*1e3;
  const double Ginv_m = P_.Gamma_kmInv/1e3; // km^-1 → m^-1
  double r_sh=0.0, V_sh=0.0;
  ::dbm_position_speed(r0_m,V0_ms,Vsw_ms,Ginv_m,t_s,r_sh,V_sh);
  S.r_sh_m = finite_or(r_sh, r0_m);
  S.V_sh_ms= finite_or(V_sh, V0_ms);
  S.a_m    = S.r_sh_m;  // semimajor axis for sphere/ellipsoid

  // Self-similar thicknesses (scale ~ r_sh), and smoothing widths
  const double scaleR = S.r_sh_m / AU;
  S.dr_sheath_m = finite_or(P_.sheath_thick_AU_at1AU * scaleR * AU);
  S.dr_me_m     = finite_or(P_.ejecta_thick_AU_at1AU * scaleR * AU);
  S.w_shock_m   = finite_or(P_.edge_smooth_shock_AU_at1AU * scaleR * AU);
  S.w_le_m      = finite_or(P_.edge_smooth_le_AU_at1AU    * scaleR * AU);
  S.w_te_m      = finite_or(P_.edge_smooth_te_AU_at1AU    * scaleR * AU);
  S.r_le_m = S.r_sh_m - S.dr_sheath_m;
  S.r_te_m = S.r_le_m - S.dr_me_m;

  // Target velocities in regions (absolute magnitudes)
  S.V_sheath_LE_ms = finite_or(P_.V_sheath_LE_factor * Vsw_ms, Vsw_ms);
  S.V_ME_ms        = finite_or(P_.V_ME_factor        * Vsw_ms, Vsw_ms);
  S.V_dn_ms        = Vsw_ms; // ambient downstream target

  // Apex rc diagnostic (not used directly in mixing; for info)
  double n_hat_apex[3]={e1[0],e1[1],e1[2]};
  double u_apex[3]    ={e1[0],e1[1],e1[2]};
  double rc_ap=1.0,Vn=0.0,dum=0.0;
  local_oblique_rc(S,u_apex,n_hat_apex,S.r_sh_m,S.r_sh_m,rc_ap,Vn,dum);
  S.rc = finite_or(rc_ap,1.0);

  // Convenience scalars used in tight loops
  S.inv_dr_sheath = (S.dr_sheath_m>0.0)? 1.0/S.dr_sheath_m : 0.0;
  S.rc_floor      = (P_.sheath_comp_floor>1.0)? P_.sheath_comp_floor : 1.0;
  return S;
}

// ----------------------------------------------------------------------------
// shape_radius_normal(S, u) → Rdir, n̂
//   Eval shock radius and outward unit normal along Cartesian unit direction u.
//   Shape: Sphere (constant), Ellipsoid (quadric solve), or Cone-like (SSE).
// ----------------------------------------------------------------------------
void Model::shape_radius_normal(const StepState& S,
                                double ux,double uy,double uz,
                                double& Rdir_m,double n_hat[3]) const {
  double u[3]={ux,uy,uz}; ::safe_normalize(u);

  // Components in apex frame
  const double e1[3]={S.e1[0],S.e1[1],S.e1[2]};
  const double e2[3]={S.e2[0],S.e2[1],S.e2[2]};
  const double e3[3]={S.e3[0],S.e3[1],S.e3[2]};
  const double u1=u[0]*e1[0]+u[1]*e1[1]+u[2]*e1[2];
  const double u2=u[0]*e2[0]+u[1]*e2[1]+u[2]*e2[2];
  const double u3=u[0]*e3[0]+u[1]*e3[1]+u[2]*e3[2];

  switch(P_.shape){
    case ShockShape::Sphere:{
      // Sphere: constant radius; normal aligns with u
      Rdir_m = S.r_sh_m; n_hat[0]=u[0]; n_hat[1]=u[1]; n_hat[2]=u[2];
    } break;

    case ShockShape::Ellipsoid:{
      // Ellipsoid in apex frame: x^2/a^2 + y^2/b^2 + z^2/c^2 = 1
      const double a=S.r_sh_m;
      const double b=a*std::max(1e-3,P_.axis_ratio_y);
      const double c=a*std::max(1e-3,P_.axis_ratio_z);
      // Find λ such that (λu) lies on surface → λ = 1 / sqrt(sum(u_i^2/a_i^2))
      const double denom=(u1*u1)/(a*a)+(u2*u2)/(b*b)+(u3*u3)/(c*c);
      const double lam=(denom>0)? 1.0/std::sqrt(denom):0.0;
      Rdir_m=lam;

      // Outward normal ∝ (x/a^2, y/b^2, z/c^2) in apex frame, then rotate to global
      const double x=Rdir_m*u1, y=Rdir_m*u2, z=Rdir_m*u3;
      double nloc[3]={x/(a*a), y/(b*b), z/(c*c)};
      double ng[3]={ nloc[0]*e1[0]+nloc[1]*e2[0]+nloc[2]*e3[0],
                     nloc[0]*e1[1]+nloc[1]*e2[1]+nloc[2]*e3[1],
                     nloc[0]*e1[2]+nloc[1]*e2[2]+nloc[2]*e3[2] };
      ::safe_normalize(ng);
      n_hat[0]=ng[0]; n_hat[1]=ng[1]; n_hat[2]=ng[2];
    } break;

    case ShockShape::ConeSSE:{
      // Cone-like (SSE): radius scales with cos(theta)^α inside a half-angle.
      // Slower flanks are mimicked by exponent "flank_slowdown_m".
      const double theta=std::acos(std::max(-1.0,std::min(1.0,u1)));
      double R=0.0;
      if (theta>P_.half_width_rad){
        const double cw=std::cos(P_.half_width_rad);
        R=S.r_sh_m*std::pow(std::max(0.0,cw), std::max(0.0,P_.flank_slowdown_m));
      } else {
        const double c=std::cos(theta);
        R=S.r_sh_m*std::pow(std::max(0.0,c), std::max(0.0,P_.flank_slowdown_m));
      }
      Rdir_m=R; n_hat[0]=u[0]; n_hat[1]=u[1]; n_hat[2]=u[2];
    } break;
  }
}

// ----------------------------------------------------------------------------
// local_oblique_rc: proxy compression & normal shock speed at r_eval_m
//   Inputs:
//     u:     radial unit vector (global, from Sun to point)
//     n_hat: local outward shock normal
//     r_eval_m: radius at which upstream quantities are evaluated
//   Outputs:
//     rc_out    : compression proxy (1..≈4)
//     Vsh_n_out : local normal shock speed [m/s]
//     thetaBn_out: obliquity (optional diag)
// ----------------------------------------------------------------------------
void Model::local_oblique_rc(const StepState& S, const double u[3], const double n_hat[3],
                             double /*Rdir_m*/, double r_eval_m,
                             double& rc_out, double& Vsh_n_out, double& thetaBn_out) const {
  using namespace swcme3d;
  const double Vsw=P_.V_sw_kms*1e3;

  // Upstream density via Leblanc at r_eval_m; scale to match n(1 AU)=n1AU_cm3
  const double A=3.3e5, B=4.1e6, C=8.0e7; // cm^-3 coefficients
  const double r0=std::max(r_eval_m, 1.05*Rs); // avoid inside coronal floor
  const double s=Rs/r0, inv2=s*s, inv4=inv2*inv2, inv6=inv4*inv2;
  const double sAU=Rs/AU;
  const double n1_base=A*sAU*sAU + B*std::pow(sAU,4) + C*std::pow(sAU,6);
  const double leb_scale=(n1_base>0.0)? (P_.n1AU_cm3/n1_base):1.0;
  const double n_up_m3=finite_or((leb_scale*(A*inv2+B*inv4+C*inv6))*1e6,1e6); // to [m^-3]
  const double rho=std::max(1e-12,n_up_m3)*MP;

  // Speeds: sound cs, upstream Parker B → vA, then θBn
  const double cs=std::sqrt(std::max(0.0,P_.gamma_ad)*KB*std::max(0.0,P_.T_K)/MP);
  double B_up[3]; ::parker_vec_T(P_,u,r0,Vsw,B_up);
  const double Bmag=std::sqrt(std::max(0.0,B_up[0]*B_up[0]+B_up[1]*B_up[1]+B_up[2]*B_up[2]));
  const double vA=(rho>0.0)? (Bmag/std::sqrt(MU0*rho)):0.0;

  double b_hat[3]={0,0,0}; if (Bmag>0){ b_hat[0]=B_up[0]/Bmag; b_hat[1]=B_up[1]/Bmag; b_hat[2]=B_up[2]/Bmag; }
  const double cosBn=std::fabs(b_hat[0]*n_hat[0]+b_hat[1]*n_hat[1]+b_hat[2]*n_hat[2]);
  thetaBn_out=finite_or(std::acos(std::max(-1.0,std::min(1.0,cosBn))),0.0);

  // Fast-mode speed c_f (cold limit formula with θBn)
  const double a=vA*vA+cs*cs;
  const double disc=std::max(0.0, a*a - 4.0*cs*cs*vA*vA*cosBn*cosBn);
  const double cf=std::sqrt(0.5*(a+std::sqrt(disc)));

  // Normal component of shock speed (directional)
  double Vsh_dir=S.V_sh_ms;
  if (P_.shape==ShockShape::ConeSSE){
    const double u1=u[0]*S.e1[0]+u[1]*S.e1[1]+u[2]*S.e1[2];
    Vsh_dir=S.V_sh_ms*std::pow(std::max(0.0,u1), std::max(0.0,P_.flank_slowdown_m));
  }
  Vsh_n_out=finite_or(Vsh_dir*(n_hat[0]*u[0]+n_hat[1]*u[1]+n_hat[2]*u[2]),0.0);

  // Upstream normal flow relative to shock
  const double Vsw_n=Vsw*(u[0]*n_hat[0]+u[1]*n_hat[1]+u[2]*n_hat[2]);
  const double U1n = std::max(0.0, Vsh_n_out - Vsw_n);
  const double Mfn = (cf>0.0)? (U1n/cf):0.0;

  // Hydrodynamic oblique proxy for rc (cap at ~4)
  double rc=1.0;
  if (Mfn>1.0){
    const double g=std::max(1.01,P_.gamma_ad), M2=Mfn*Mfn;
    rc=((g+1.0)*M2)/((g-1.0)*M2+2.0);
    if (rc>4.0) rc=4.0;
  }
  rc_out=finite_or(rc,1.0);
}

// ----------------------------------------------------------------------------
// evaluate_cartesian_fast: returns (n, V) for N Cartesian points
//   Uses smooth blends across shock/sheath/ejecta based on the local rc proxy.
//   velocity = V(r) û with piecewise targets (ambient, sheath-LE, ejecta).
// ----------------------------------------------------------------------------
void Model::evaluate_cartesian_fast(const StepState& S,
                                    const double* x_m,const double* y_m,const double* z_m,
                                    double* n_m3,double* Vx_ms,double* Vy_ms,double* Vz_ms,
                                    std::size_t N) const {
  const double Vup=P_.V_sw_kms*1e3;

  // Radial self-similar widths
  const double scaleR=S.r_sh_m/AU;
  const double dr_sheath=P_.sheath_thick_AU_at1AU*scaleR*AU;
  const double dr_me    =P_.ejecta_thick_AU_at1AU*scaleR*AU;

  // Smoothing widths at each interface
  const double w_sh=S.w_shock_m, w_le=S.w_le_m, w_te=S.w_te_m;

  // Target speeds in sheath (at LE) and ejecta
  const double Vshe_LE=S.V_sheath_LE_ms, Vme=S.V_ME_ms;

  // Leblanc coefficients and scaling to 1 AU density
  const double A=3.3e5, B=4.1e6, C=8.0e7, sAU=Rs/AU;
  const double n1_base=A*sAU*sAU+B*std::pow(sAU,4)+C*std::pow(sAU,6);
  const double leb_scale=(n1_base>0.0)? (P_.n1AU_cm3/n1_base):1.0;

  // Sheath compression ramp power and floor
  const double p=(P_.sheath_ramp_power<1.0)? 1.0 : P_.sheath_ramp_power;
  const double rc_floor=(P_.sheath_comp_floor<1.0)? 1.0 : P_.sheath_comp_floor;

  // Helper for edge activations
  auto edge=[&](double rnow,double r0,double w)->double{
    if (w<=0.0) return (rnow<=r0)? 1.0:0.0;
    return ::smoothstep01(0.5+(r0-rnow)/(2.0*w));
  };

  for (std::size_t i=0;i<N;++i){
    // Geometry & direction
    const double x=x_m[i], y=y_m[i], z=z_m[i];
    const double r=std::sqrt(std::max(1e-12,x*x+y*y+z*z));
    const double invr=1.0/r;
    double u[3]={x*invr,y*invr,z*invr}; // unit radial

    // Shock radius and outward normal along u
    double Rdir=0.0,n_hat[3]={0,0,1}; shape_radius_normal(S,u[0],u[1],u[2],Rdir,n_hat);

    // Upstream density along the ray
    const double s=Rs/r, inv2=s*s, inv4=inv2*inv2, inv6=inv4*inv2;
    const double n_up=finite_or((leb_scale*(A*inv2+B*inv4+C*inv6))*1e6,1e6);

    // Local shock metrics (for mixing) evaluated at max(r,Rdir)
    double rc_loc=1.0,Vsh_n=0.0,dum=0.0;
    local_oblique_rc(S,u,n_hat,Rdir,(r>Rdir?r:Rdir),rc_loc,Vsh_n,dum);

    // Position within sheath measured inward from Rdir
    const double xi=(dr_sheath>0.0)? swcme3d::clamp01((Rdir-r)/dr_sheath):0.0;

    // Density target in sheath—smooth ramp from rc_floor to rc_loc
    const double ramp=(p==1.0)? (1.0-xi): std::pow(1.0-xi,p);
    const double Csheath=rc_floor + (rc_loc-rc_floor)*ramp;
    const double n_sheath=Csheath*n_up;

    // Velocity targets in sheath/ejecta regions
    const double V_sheath=Vup + (Vshe_LE - Vup)*xi;
    const double n_ejecta=(P_.f_ME>0.0? P_.f_ME:0.0)*n_up;
    const double V_ejecta=Vme;

    // Blend across interfaces (shock, leading, trailing)
    double a=(rc_loc>1.0)? edge(r,Rdir,w_sh):0.0;
    double n_mix=(1.0-a)*n_up + a*n_sheath;
    double Vmag =(1.0-a)*Vup  + a*V_sheath;

    a=(rc_loc>1.0)? edge(r,Rdir-dr_sheath,w_le):0.0;
    n_mix=(1.0-a)*n_mix + a*n_ejecta;
    Vmag =(1.0-a)*Vmag + a*V_ejecta;

    a=(rc_loc>1.0)? edge(r,Rdir-dr_sheath-dr_me,w_te):0.0;
    n_mix=(1.0-a)*n_mix + a*n_up;
    Vmag =(1.0-a)*Vmag + a*Vup;

    // Store results (finite-guarded)
    n_m3[i]=finite_or(n_mix,n_up);
    Vx_ms[i]=finite_or(Vmag*u[0],0.0);
    Vy_ms[i]=finite_or(Vmag*u[1],0.0);
    Vz_ms[i]=finite_or(Vmag*u[2],0.0);
  }
}

// ----------------------------------------------------------------------------
// evaluate_cartesian_with_B: returns (n, V, B) with sheath Bt amplification
//   B_up is Parker. Within the sheath we amplify Bt from 1 → rc using the same
//   ramp used for density, while keeping Bn continuous (Bn assumed ~ const).
// ----------------------------------------------------------------------------
void Model::evaluate_cartesian_with_B(const StepState& S,
  const double* x_m,const double* y_m,const double* z_m,
  double* n_m3,double* Vx_ms,double* Vy_ms,double* Vz_ms,
  double* Bx_T,double* By_T,double* Bz_T,
  std::size_t N) const {

  const double Vup=P_.V_sw_kms*1e3;
  const double scaleR=S.r_sh_m/AU;
  const double dr_sheath=P_.sheath_thick_AU_at1AU*scaleR*AU;
  const double dr_me    =P_.ejecta_thick_AU_at1AU*scaleR*AU;
  const double w_sh=S.w_shock_m, w_le=S.w_le_m, w_te=S.w_te_m;
  const double Vshe_LE=S.V_sheath_LE_ms, Vme=S.V_ME_ms;

  const double A=3.3e5, B=4.1e6, C=8.0e7, sAU=Rs/AU;
  const double n1_base=A*sAU*sAU+B*std::pow(sAU,4)+C*std::pow(sAU,6);
  const double leb_scale=(n1_base>0.0)? (P_.n1AU_cm3/n1_base):1.0;
  const double p=(P_.sheath_ramp_power<1.0)? 1.0 : P_.sheath_ramp_power;
  const double rc_floor=(P_.sheath_comp_floor<1.0)? 1.0 : P_.sheath_comp_floor;

  auto edge=[&](double rnow,double r0,double w)->double{
    if (w<=0.0) return (rnow<=r0)? 1.0:0.0;
    return ::smoothstep01(0.5+(r0-rnow)/(2.0*w));
  };

  for (std::size_t i=0;i<N;++i){
    const double x=x_m[i], y=y_m[i], z=z_m[i];
    const double r=std::sqrt(std::max(1e-12,x*x+y*y+z*z));
    const double invr=1.0/r;
    double u[3]={x*invr,y*invr,z*invr};

    double Rdir=0.0,n_hat[3]={0,0,1}; shape_radius_normal(S,u[0],u[1],u[2],Rdir,n_hat);

    const double s=Rs/r, inv2=s*s, inv4=inv2*inv2, inv6=inv4*inv2;
    const double n_up=finite_or((leb_scale*(A*inv2+B*inv4+C*inv6))*1e6,1e6);

    // Upstream Parker B at the actual r (not at Rdir)
    double B_up[3]; ::parker_vec_T(P_,u,r,Vup,B_up);

    // Local rc and Vsh_n at max(r, Rdir)
    double rc_loc=1.0,Vsh_n=0.0,dum=0.0;
    local_oblique_rc(S,u,n_hat,Rdir,(r>Rdir?r:Rdir),rc_loc,Vsh_n,dum);

    // Sheath "xi" coordinate and density target in sheath
    const double xi=(dr_sheath>0.0)? swcme3d::clamp01((Rdir-r)/dr_sheath):0.0;
    const double ramp=(p==1.0)? (1.0-xi): std::pow(1.0-xi,p);
    const double Csheath=rc_floor + (rc_loc-rc_floor)*ramp;
    const double n_sheath=Csheath*n_up;

    // Velocity targets
    const double V_sheath=Vup + (Vshe_LE - Vup)*xi;
    const double n_ejecta=(P_.f_ME>0.0? P_.f_ME:0.0)*n_up;
    const double V_ejecta=Vme;

    // Smooth blends across edges (same as n,V only function)
    double a=(rc_loc>1.0)? edge(r,Rdir,w_sh):0.0;
    double n_mix=(1.0-a)*n_up + a*n_sheath;
    double Vmag =(1.0-a)*Vup  + a*V_sheath;

    a=(rc_loc>1.0)? edge(r,Rdir-dr_sheath,w_le):0.0;
    n_mix=(1.0-a)*n_mix + a*n_ejecta;
    Vmag =(1.0-a)*Vmag + a*V_ejecta;

    a=(rc_loc>1.0)? edge(r,Rdir-dr_sheath-dr_me,w_te):0.0;
    n_mix=(1.0-a)*n_mix + a*n_up;
    Vmag =(1.0-a)*Vmag + a*Vup;

    // Magnetic field splitting into normal and tangential to shock
    const double Bn_mag=B_up[0]*n_hat[0]+B_up[1]*n_hat[1]+B_up[2]*n_hat[2];
    double Bn[3]={Bn_mag*n_hat[0],Bn_mag*n_hat[1],Bn_mag*n_hat[2]};
    double Bt[3]={B_up[0]-Bn[0],B_up[1]-Bn[1],B_up[2]-Bn[2]};

    // Amplify tangential B smoothly with the same shock-edge activation and ramp
    const double a_shock=(rc_loc>1.0)? edge(r,Rdir,w_sh):0.0;
    const double amp_t=1.0 + (rc_loc-1.0)*a_shock*ramp;
    double Bv[3]={Bn[0]+amp_t*Bt[0],Bn[1]+amp_t*Bt[1],Bn[2]+amp_t*Bt[2]};

    // Store outputs (finite-protected)
    n_m3[i]=finite_or(n_mix,n_up);
    Vx_ms[i]=finite_or(Vmag*u[0],0.0);
    Vy_ms[i]=finite_or(Vmag*u[1],0.0);
    Vz_ms[i]=finite_or(Vmag*u[2],0.0);
    Bx_T[i] =finite_or(Bv[0],0.0);
    By_T[i] =finite_or(Bv[1],0.0);
    Bz_T[i] =finite_or(Bv[2],0.0);
  }
}

// ----------------------------------------------------------------------------
// evaluate_cartesian_with_B_div: same as above + ∇·V
//   We compute divergence using a radial finite-difference of r^2 V_r.
// ----------------------------------------------------------------------------
void Model::evaluate_cartesian_with_B_div(const StepState& S,
  const double* x_m,const double* y_m,const double* z_m,
  double* n_m3,double* Vx_ms,double* Vy_ms,double* Vz_ms,
  double* Bx_T,double* By_T,double* Bz_T,double* divVsw,
  std::size_t N,double dr_frac) const {
  evaluate_cartesian_with_B(S,x_m,y_m,z_m,n_m3,Vx_ms,Vy_ms,Vz_ms,Bx_T,By_T,Bz_T,N);
  compute_divV_radial(S,x_m,y_m,z_m,divVsw,N,dr_frac);
}

// ----------------------------------------------------------------------------
// compute_divV_radial: robust radial ∇·V = (1/r^2) d(r^2 V_r)/dr
//   We evaluate V at r±dr along the same unit direction.
//   dr = max(dr_min, dr_frac * r) with floor to avoid singular behavior.
// ----------------------------------------------------------------------------
void Model::compute_divV_radial(const StepState& S,
  const double* x_m,const double* y_m,const double* z_m,
  double* divV,std::size_t N,double dr_frac) const {
  const double rmin=1.05*Rs, dr_min=1.0e-4*AU;
  for (std::size_t i=0;i<N;++i){
    const double x=x_m[i], y=y_m[i], z=z_m[i];
    const double r=std::sqrt(std::max(1e-12,x*x+y*y+z*z));
    const double invr=1.0/r;
    const double u[3]={x*invr,y*invr,z*invr};
    const double dr=std::max(dr_min, dr_frac*r);
    const double rp=std::max(rmin,r+dr), rm=std::max(rmin,r-dr);
    const double denom_r=(rp>rm)? (rp-rm): std::max(dr_min,std::abs(dr));

    // Evaluate V at r±dr
    double xp=rp*u[0], yp=rp*u[1], zp=rp*u[2];
    double xm=rm*u[0], ym=rm*u[1], zm=rm*u[2];
    double n,Vxp,Vyp,Vzp,Vxm,Vym,Vzm;
    evaluate_cartesian_fast(S,&xp,&yp,&zp,&n,&Vxp,&Vyp,&Vzp,1);
    evaluate_cartesian_fast(S,&xm,&ym,&zm,&n,&Vxm,&Vym,&Vzm,1);
    const double Vrp=Vxp*u[0]+Vyp*u[1]+Vzp*u[2];
    const double Vrm=Vxm*u[0]+Vym*u[1]+Vzm*u[2];

    // Central difference on r^2 V_r
    const double num=((rp*rp)*Vrp - (rm*rm)*Vrm)/denom_r;
    const double div_val=num/(r*r);
    divV[i]=finite_or(div_val,0.0);
  }
}

// ----------------------------------------------------------------------------
// diagnose_direction: convenience diagnostic (Rdir, n̂, rc, Vsh_n) along u
// ----------------------------------------------------------------------------
void Model::diagnose_direction(const StepState& S,const double u[3],
  double& Rdir_m,double n_hat[3],double& rc_loc,double& Vsh_n) const {
  shape_radius_normal(S,u[0],u[1],u[2],Rdir_m,n_hat);
  double th=0.0; local_oblique_rc(S,u,n_hat,Rdir_m,Rdir_m,rc_loc,Vsh_n,th);
}

// ----------------------------------------------------------------------------
// build_shock_mesh: lat–lon sampling of shock surface + nodal attributes
//   Returns a ShockMesh with vertices (x,y,z), nodal normals, nodal rc, Vsh_n,
//   and triangle connectivity (1-based indices).
// ----------------------------------------------------------------------------
ShockMesh Model::build_shock_mesh(const StepState& S,std::size_t nTheta,std::size_t nPhi) const {
  ShockMesh M; if (nTheta<3) nTheta=3; if (nPhi<3) nPhi=3;
  const double thetaMax=(P_.shape==ShockShape::ConeSSE)? P_.half_width_rad : PI;

  // Sample grid in apex frame then rotate to global
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

      double Rdir=0.0,n_hat[3]={0,0,1}; shape_radius_normal(S,u[0],u[1],u[2],Rdir,n_hat);

      // Nodal positions & normals
      M.x.push_back(finite_or(Rdir*u[0],0.0));
      M.y.push_back(finite_or(Rdir*u[1],0.0));
      M.z.push_back(finite_or(Rdir*u[2],0.0));
      M.n_hat_x.push_back(finite_or(n_hat[0],0.0));
      M.n_hat_y.push_back(finite_or(n_hat[1],0.0));
      M.n_hat_z.push_back(finite_or(n_hat[2],1.0));

      // Nodal rc and Vsh_n (diagnostic values used by writers/metrics)
      double rc_loc=1.0,Vsh_n=0.0,th=0.0;
      local_oblique_rc(S,u,n_hat,Rdir,Rdir,rc_loc,Vsh_n,th);
      M.rc.push_back(finite_or(rc_loc,1.0));
      M.Vsh_n.push_back(finite_or(Vsh_n,0.0));
    }
  }

  // 1-based triangle connectivity with two tris per parametric quad
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

// ----------------------------------------------------------------------------
// compute_triangle_metrics: per-cell area, geometric normal, centroid,
//   and simple per-cell means of nodal rc and Vsh_n.
// ----------------------------------------------------------------------------
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

    // Geometric normal from cross(AB, AC)
    double ABx=Bx-Ax, ABy=By-Ay, ABz=Bz-Az;
    double ACx=Cx-Ax, ACy=Cy-Ay, ACz=Cz-Az;
    double nx=ABy*ACz - ABz*ACy;
    double ny=ABz*ACx - ABx*ACz;
    double nz=ABx*ACy - ABy*ACx;
    double twiceA=std::sqrt(std::max(0.0, nx*nx+ny*ny+nz*nz));
    double area=0.5*twiceA;
    if (twiceA>0){ nx/=twiceA; ny/=twiceA; nz/=twiceA; }

    // Store metrics (finite-guarded)
    T.area[e]=finite_or(area,0.0);
    T.nx[e]=finite_or(nx,0.0); T.ny[e]=finite_or(ny,0.0); T.nz[e]=finite_or(nz,1.0);
    T.cx[e]=finite_or((Ax+Bx+Cx)/3.0,0.0);
    T.cy[e]=finite_or((Ay+By+Cy)/3.0,0.0);
    T.cz[e]=finite_or((Az+Bz+Cz)/3.0,0.0);
    T.rc_mean[e]=finite_or((M.rc[ia]+M.rc[ib]+M.rc[ic])/3.0,1.0);
    T.Vsh_n_mean[e]=finite_or((M.Vsh_n[ia]+M.Vsh_n[ib]+M.Vsh_n[ic])/3.0,0.0);
  }
}

// ─────────────────────────────────────────────────────────────────────────────
// Tecplot writers
// ─────────────────────────────────────────────────────────────────────────────
//
// We output a dataset with a consistent VARIABLES list across all zones:
//
//   "X","Y","Z", "n","Vx","Vy","Vz", "Bx","By","Bz","divVsw",
//   "rc","Vsh_n", "nx","ny","nz","area","rc_mean","Vsh_n_mean",
//   "tnx","tny","tnz", "cx","cy","cz"
//
// Notes:
//  • All values are finite-sanitized via finite_or(...).
//  • "surface_cells" (FETRIANGLE/BLOCK) appears FIRST so that Tecplot shows
//    meaningful cell-centered scalars by default (area, rc_mean, ...).
//  • "surface_nodal" (FEPOINT) is a purely nodal view (cell vars printed as 0).
//  • "volume_box" and optional "box_face_minX" are structured POINT zones.
// ─────────────────────────────────────────────────────────────────────────────

// FE surface writer with cell metrics (BLOCK). Auto-computes T if needed.
bool Model::write_shock_surface_center_metrics_tecplot(
  const ShockMesh& M, const TriMetrics& T_in, const char* path) const {

  TriMetrics T = T_in;
  const std::size_t Nv=M.x.size(), Ne=M.tri_i.size();
  if (T.area.size()!=Ne){ compute_triangle_metrics(M,T); }

  std::FILE* fp=std::fopen(path,"w"); if(!fp) return false;

  std::fprintf(fp,"TITLE=\"Shock surface (cell metrics + nodal rc)\"\n");
  std::fprintf(fp,"VARIABLES=\"X\",\"Y\",\"Z\",\"rc\",\"Vsh_n\","
                  "\"area\",\"rc_mean\",\"Vsh_n_mean\",\"tnx\",\"tny\",\"tnz\",\"cx\",\"cy\",\"cz\"\n");
  std::fprintf(fp,"ZONE T=\"surface_cells\", N=%zu, E=%zu, ZONETYPE=FETRIANGLE, DATAPACKING=BLOCK,\n", Nv, Ne);
  std::fprintf(fp,"VARLOCATION=([1-5]=NODAL, [6-14]=CELLCENTERED)\n");

  auto dumpN=[&](const std::vector<double>& a){ int cnt=0; for(double v:a){ std::fprintf(fp,"%.9e ",finite_or(v)); if(++cnt==8){std::fprintf(fp,"\n"); cnt=0;} } if(cnt) std::fprintf(fp,"\n"); };
  auto dumpE=dumpN;

  dumpN(M.x); dumpN(M.y); dumpN(M.z); dumpN(M.rc); dumpN(M.Vsh_n);
  dumpE(T.area); dumpE(T.rc_mean); dumpE(T.Vsh_n_mean);
  dumpE(T.nx); dumpE(T.ny); dumpE(T.nz);
  dumpE(T.cx); dumpE(T.cy); dumpE(T.cz);

  for (std::size_t e=0;e<Ne;++e)
    std::fprintf(fp,"%d %d %d\n", M.tri_i[e], M.tri_j[e], M.tri_k[e]);

  std::fclose(fp);
  return true;
}

// Default apex-centered box (shifted outward along e1 so it samples shock/sheath)
BoxSpec Model::default_apex_box(const StepState& S,double half_AU,int N) const {
  BoxSpec B; const double h=half_AU*AU;
  B.hx=h; B.hy=h; B.hz=h;

  // Shift the box outward (40% of half-size) so the min-X face can pass the Sun
  const double shift=0.4*h;
  B.cx=S.a_m*S.e1[0]+shift*S.e1[0];
  B.cy=S.a_m*S.e1[1]+shift*S.e1[1];
  B.cz=S.a_m*S.e1[2]+shift*S.e1[2];
  B.Ni=N; B.Nj=N; B.Nk=N; return B;
}

// Full dataset bundle (surface_cells + surface_nodal + volume_box + face)
bool Model::write_tecplot_dataset_bundle(const ShockMesh& M,const TriMetrics& T_in,
                                         const StepState& S,const BoxSpec& B,
                                         const char* path) const {
  TriMetrics T=T_in;
  const std::size_t Nv=M.x.size(), Ne=M.tri_i.size();
  if (T.area.size()!=Ne){ compute_triangle_metrics(M,T); }

  std::FILE* fp=std::fopen(path,"w"); if(!fp) return false;
  auto p=[&](const char* fmt, auto... args){ std::fprintf(fp,fmt,args...); };

  // Shared VARIABLES list for all zones in this dataset
  p("TITLE = \"SW+CME dataset\"\n");
  p("VARIABLES = "
    "\"X\",\"Y\",\"Z\","
    "\"n\",\"Vx\",\"Vy\",\"Vz\","
    "\"Bx\",\"By\",\"Bz\",\"divVsw\","
    "\"rc\",\"Vsh_n\","
    "\"nx\",\"ny\",\"nz\",\"area\",\"rc_mean\",\"Vsh_n_mean\","
    "\"tnx\",\"tny\",\"tnz\",\"cx\",\"cy\",\"cz\"\n");

  // ---- Zone 1: surface_cells (BLOCK with cell-centered fields) ----
  p("ZONE T=\"surface_cells\", N=%zu, E=%zu, ZONETYPE=FETRIANGLE, DATAPACKING=BLOCK,\n", Nv, Ne);
  p("VARLOCATION=([1-3,12-13]=NODAL, [4-11,14-25]=CELLCENTERED)\n");

  auto dumpN=[&](const std::vector<double>& a){ int cnt=0; for(double v:a){ std::fprintf(fp,"%.9e ",finite_or(v)); if(++cnt==8){std::fprintf(fp,"\n"); cnt=0;} } if(cnt) std::fprintf(fp,"\n"); };
  auto dumpE=dumpN;
  auto dumpEzeros=[&](std::size_t count){ int cnt=0; for(std::size_t e=0;e<count;++e){ std::fprintf(fp,"0 "); if(++cnt==8){std::fprintf(fp,"\n"); cnt=0;} } if(cnt) std::fprintf(fp,"\n"); };

  // Block order: 25 arrays
  dumpN(M.x); dumpN(M.y); dumpN(M.z);                    // 1..3 nodal positions
  dumpEzeros(Ne); dumpEzeros(Ne); dumpEzeros(Ne);        // 4..6 n,Vx,Vy (cell-only here → 0)
  dumpEzeros(Ne);                                        // 7 Vz
  dumpEzeros(Ne); dumpEzeros(Ne); dumpEzeros(Ne);        // 8..10 Bx,By,Bz
  dumpEzeros(Ne);                                        // 11 divVsw
  dumpN(M.rc); dumpN(M.Vsh_n);                           // 12..13 nodal shock diag
  dumpE(T.nx); dumpE(T.ny); dumpE(T.nz);                 // 14..16 cell normals
  dumpE(T.area); dumpE(T.rc_mean); dumpE(T.Vsh_n_mean);  // 17..19 cell metrics
  dumpEzeros(Ne); dumpEzeros(Ne); dumpEzeros(Ne);        // 20..22 reserved
  dumpE(T.cx); dumpE(T.cy); dumpE(T.cz);                 // 23..25 centroids

  for (std::size_t e=0;e<Ne;++e)
    p("%d %d %d\n", M.tri_i[e], M.tri_j[e], M.tri_k[e]);

  // ---- Zone 2: surface_nodal (FEPOINT). Cell-centered variables are 0. ----
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
      0.0,0.0,0.0,0.0,          // n,V not defined in this zone
      0.0,0.0,0.0,0.0,          // B,div not defined here
      finite_or(M.rc[i],1.0), finite_or(M.Vsh_n[i],0.0),      // nodal rc, Vsh_n
      finite_or(M.n_hat_x[i],0.0), finite_or(M.n_hat_y[i],0.0), finite_or(M.n_hat_z[i],1.0),
      0.0,0.0,0.0,              // area, rc_mean, Vsh_n_mean (cell vars)
      0.0,0.0,0.0,              // reserved
      0.0,0.0,0.0               // centroid (cell var)
    );
  }
  for (std::size_t e=0;e<Ne;++e)
    p("%d %d %d\n", M.tri_i[e], M.tri_j[e], M.tri_k[e]);

  // ---- Zone 3: volume_box (structured POINT). Writes full fields. ----
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
          0.0,0.0,    // rc, Vsh_n (not defined for volume)
          0.0,0.0,0.0,// nx,ny,nz
          0.0,0.0,0.0,// area, rc_mean, Vsh_n_mean
          0.0,0.0,0.0,// reserved
          0.0,0.0,0.0 // centroid
        );
      }
    }
  }

  // ---- Zone 4: optional min-X face (structured 2-D POINT). ----
  if (ADD_FACE_ZONE_IN_BUNDLE){
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
          0.0,0.0,  // rc, Vsh_n
          0.0,0.0,0.0, // nx,ny,nz
          0.0,0.0,0.0, // area, rc_mean, Vsh_n_mean
          0.0,0.0,0.0, // reserved
          0.0,0.0,0.0  // centroid
        );
      }
    }
  }

  std::fclose(fp);
  return true;
}

// Standalone writer for the min-X face (same variables ordering)
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
        0.0,0.0,  // rc, Vsh_n not defined on this 2D slice
        0.0,0.0,0.0, // nx,ny,nz
        0.0,0.0,0.0, // area, rc_mean, Vsh_n_mean
        0.0,0.0,0.0, // reserved
        0.0,0.0,0.0  // centroid
      );
    }
  }
  std::fclose(fp);
  return true;
}

} // namespace swcme3d

