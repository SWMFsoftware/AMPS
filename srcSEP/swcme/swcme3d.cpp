// swcme3d.cpp
//
// ─────────────────────────────────────────────────────────────────────────────
// Solar-wind + CME forward-shock model with Parker spiral B, DBM kinematics,
// sheath/ejecta phenomenology, oblique MHD shock proxy, ∇·Vsw, and Tecplot I/O.
// Includes NaN/Inf sanitization to prevent Tecplot connectivity/parse errors.
// ─────────────────────────────────────────────────────────────────────────────
//
// TECPlot VARIABLES (consistent order in dataset bundle)
// -----------------------------------------------------
// VARIABLES =
//  1  "X"       [m]   : global Cartesian position
//  2  "Y"       [m]
//  3  "Z"       [m]
//  4  "n"     [m^-3]  : proton number density
//  5  "Vx"    [m/s]   : bulk velocity components of solar-wind/CME model
//  6  "Vy"    [m/s]
//  7  "Vz"    [m/s]
//  8  "Bx"      [T]   : magnetic field components; Parker upstream with tangential
//  9  "By"      [T]     amplification inside sheath toward rc (compressed Bt)
// 10  "Bz"      [T]
// 11  "divVsw" [1/s]  : ∇·Vsw using a radial finite-difference approximation
// 12  "rc"      [-]   : shock compression ratio (NODAL on surface zones; 0 elsewhere)
// 13  "Vsh_n"  [m/s]  : shock propagation speed projected on local normal (surface nodes)
// 14  "nx"      [-]   : nodal outward unit normal (surface_nodal) or triangle normal (surface_cells)
// 15  "ny"      [-]
// 16  "nz"      [-]
// 17  "area"   [m^2]  : per-triangle area (CELLCENTERED; surface_cells only)
// 18  "rc_mean"[-]    : per-triangle mean of rc over its 3 vertices (surface_cells)
// 19  "Vsh_n_mean"[m/s]: per-triangle mean Vsh_n (surface_cells)
// 20  "tnx"     [-]   : reserved triangle normal components (unused → 0 in bundle)
// 21  "tny"     [-]
// 22  "tnz"     [-]
// 23  "cx"      [m]   : triangle centroid X (surface_cells)
// 24  "cy"      [m]   : triangle centroid Y
// 25  "cz"      [m]   : triangle centroid Z
//
// PHYSICS SUMMARY
// ---------------
// • Ambient density n(r): Leblanc et al. (1998) scaled so n(1 AU)=n1AU_cm3.
//   n(r) = A (R_s/r)^2 + B (R_s/r)^4 + C (R_s/r)^6, with A,B,C in cm^-3; converted to m^-3.
// • Parker spiral B(r): Br ∝ r^-2, Bφ = -Br (Ω r sinθ / Vsw), normalized so |B|(1 AU)=B1AU_nT.
// • CME apex kinematics via DBM (Vršnak et al. 2013): dV/dt = -Γ (V - Vsw)|V - Vsw|,
//   closed-form V(t) and r(t).
// • Shock/ejecta structure: piecewise radial blending with C^1 smoothsteps and independent
//   smoothing half-widths at shock, leading edge, and trailing edge; density amplified by a
//   ramped compression ratio rc ≥ floor; Bt (tangential B) amplified in sheath.
// • Oblique shock proxy (Edmiston & Kennel 1984; Priest 2014): upstream Parker B, cs, vA
//   ⇒ fast speed cf ⇒ normal fast Mach Mfn ⇒ compression rc (capped at 4 for γ=5/3).
// • Divergence ∇·Vsw: computed by centered finite difference of r^2 V_r in radius (spherical
//   identity: ∇·(V_r r̂) = (1/r^2) d/dr (r^2 V_r)), with guard rails on Δr.
//
// REFERENCES
// ----------
// Parker (1958) ApJ 128:664; Howard & Tappin (2009) Space Sci. Rev. 147:31;
// Leblanc, Dulk, Bougeret (1998) Sol. Phys. 183:165;
// Vršnak et al. (2013) Sol. Phys. 285:295 (DBM);
// Edmiston & Kennel (1984) J. Plasma Phys. 32:429; Priest (2014) MHD of the Sun;
// Russell & Mulligan (2002) Planet. Space Sci. 50:527; Manchester et al. (2005) ApJ 622:1225.
//
// BUILD
// -----
//   g++ -std=c++17 -O3 -march=native demo3d_2.cpp swcme3d.cpp -o demo
//
// NOTE
// ----
// This .cpp pairs with the header previously provided (swcme3d.hpp) that declares:
//   - extern const double AU, Rs, PI;
//   - Model methods including evaluate_cartesian_with_B_div(...) and
//     write_box_face_minX_tecplot_structured(...).
// No changes to that header are required for this file.
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

// Toggle: also write a structured 2-D min-X face zone in dataset bundle.
// Structured zones do NOT require connectivity; this avoids Tecplot FE errors.
static constexpr bool ADD_FACE_ZONE_IN_BUNDLE = true;

// =========================== Namespace-level constants =========================
namespace swcme3d {
  const double AU = 1.495978707e11;     // Astronomical Unit [m]
  const double Rs = 6.957e8;            // Solar radius [m]
  const constexpr double PI = 3.141592653589793;  // Pi
  // Lightweight clamp (avoid std::clamp in older toolchains)
  inline double clamp01(double v){ return (v < 0.0) ? 0.0 : (v > 1.0 ? 1.0 : v); }
}

// =============================== File-scope helpers ===============================

// Finite-number guards (sanitize NaN/Inf → 0.0 by default)
static inline double finite_or(double v, double fallback = 0.0) {
  return std::isfinite(v) ? v : fallback;
}
static inline double finite_bounded(double v, double fallback = 0.0, double vmax = 1e300) {
  if (!std::isfinite(v)) return fallback;
  if (v >  vmax) return vmax;
  if (v < -vmax) return -vmax;
  return v;
}

static constexpr double MU0        = 4.0e-7 * swcme3d::PI; // vacuum permeability [H/m]
static constexpr double MP         = 1.67262192369e-27;    // proton mass [kg]
static constexpr double KB         = 1.380649e-23;         // Boltzmann [J/K]
static constexpr double OMEGA_SUN  = 2.86533e-6;           // solar rotation [rad/s]

// C^1 smoothstep on [0,1]
static inline double smoothstep01(double x){
  if (x <= 0.0) return 0.0;
  if (x >= 1.0) return 1.0;
  return x*x*(3.0 - 2.0*x);
}

static inline void safe_normalize(double v[3]){
  const double m = std::sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
  if (m > 0.0){ v[0] /= m; v[1] /= m; v[2] /= m; }
  else { v[0] = 1.0; v[1] = 0.0; v[2] = 0.0; }
}

// Parker spiral vector at (r,u), normalized so |B|(1 AU) = P.B1AU_nT.
// u is radial unit vector from Sun to the point.
static inline void parker_vec_T(const swcme3d::Params& P,
                                const double u[3], double r_m,
                                double Vsw_ms, double B_out[3])
{
  using namespace swcme3d;
  const double r_AU = r_m / AU;
  const double k_AU = (OMEGA_SUN * AU * P.sin_theta) / Vsw_ms;

  const double B1AU_T = P.B1AU_nT * 1e-9;
  const double Br1AU  = B1AU_T / std::sqrt(1.0 + k_AU*k_AU);  // so |B|(1 AU) = B1AU_T

  const double Br   = Br1AU / (r_AU*r_AU);
  const double Bphi = - Br * k_AU * r_AU;

  // r-hat = u ; phi-hat ∝ (z × u) × u
  const double zhat[3] = {0,0,1};
  double zxu[3] = { zhat[1]*u[2]-zhat[2]*u[1],
                    zhat[2]*u[0]-zhat[0]*u[2],
                    zhat[0]*u[1]-zhat[1]*u[0] };
  double ph[3]  = { zxu[1]*u[2]-zxu[2]*u[1],
                    zxu[2]*u[0]-zxu[0]*u[2],
                    zxu[0]*u[1]-zxu[1]*u[0] };
  safe_normalize(ph);

  B_out[0] = Br*u[0] + Bphi*ph[0];
  B_out[1] = Br*u[1] + Bphi*ph[1];
  B_out[2] = Br*u[2] + Bphi*ph[2];
}

// DBM apex position/speed (Vršnak+ 2013)
static inline void dbm_position_speed(double r0_m, double V0_ms, double Vsw_ms,
                                      double Gamma_mInv, double t_s,
                                      double& r_m, double& V_ms)
{
  const double u0   = V0_ms - Vsw_ms;
  const double den  = 1.0 + Gamma_mInv * u0 * t_s;
  const double u    = (den != 0.0) ? (u0 / den) : 0.0;
  r_m = r0_m + Vsw_ms*t_s + std::log(std::max(den, 1e-30)) / Gamma_mInv;
  V_ms = Vsw_ms + u;
}

// =============================== Implementation ===============================
namespace swcme3d {

Model::Model(const Params& P) : P_(P) {}

// Build state for time t [s]
StepState Model::prepare_step(double t_s) const
{
  StepState S{};

  // Apex-aligned orthonormal basis (e1 = CME direction)
  double e1[3] = {P_.cme_dir[0], P_.cme_dir[1], P_.cme_dir[2]};
  ::safe_normalize(e1);
  double tmp[3] = {0,0,1};
  if (std::fabs(e1[2]) > 0.9) { tmp[0]=1; tmp[1]=0; tmp[2]=0; }
  double e2[3] = { e1[1]*tmp[2]-e1[2]*tmp[1],
                   e1[2]*tmp[0]-e1[0]*tmp[2],
                   e1[0]*tmp[1]-e1[1]*tmp[0] };
  ::safe_normalize(e2);
  double e3[3] = { e1[1]*e2[2]-e1[2]*e2[1],
                   e1[2]*e2[0]-e1[0]*e2[2],
                   e1[0]*e2[1]-e1[1]*e2[0] };
  ::safe_normalize(e3);
  S.e1[0]=e1[0]; S.e1[1]=e1[1]; S.e1[2]=e1[2];
  S.e2[0]=e2[0]; S.e2[1]=e2[1]; S.e2[2]=e2[2];
  S.e3[0]=e3[0]; S.e3[1]=e3[1]; S.e3[2]=e3[2];

  // DBM apex kinematics
  const double r0_m   = P_.r0_Rs * Rs;
  const double V0_ms  = P_.V0_sh_kms * 1e3;
  const double Vsw_ms = P_.V_sw_kms * 1e3;
  const double Ginv_m = P_.Gamma_kmInv / 1e3;  // 1/km → 1/m
  double r_sh=0.0, V_sh=0.0;
  ::dbm_position_speed(r0_m, V0_ms, Vsw_ms, Ginv_m, t_s, r_sh, V_sh);
  S.r_sh_m  = finite_or(r_sh, r0_m);
  S.V_sh_ms = finite_or(V_sh, V0_ms);
  S.a_m     = S.r_sh_m;

  // Self-similar thickness & smoothing widths
  const double scaleR = S.r_sh_m / AU;
  S.dr_sheath_m = finite_or(P_.sheath_thick_AU_at1AU * scaleR * AU);
  S.dr_me_m     = finite_or(P_.ejecta_thick_AU_at1AU * scaleR * AU);
  S.w_shock_m   = finite_or(P_.edge_smooth_shock_AU_at1AU * scaleR * AU);
  S.w_le_m      = finite_or(P_.edge_smooth_le_AU_at1AU    * scaleR * AU);
  S.w_te_m      = finite_or(P_.edge_smooth_te_AU_at1AU    * scaleR * AU);

  S.r_le_m = S.r_sh_m - S.dr_sheath_m;
  S.r_te_m = S.r_le_m - S.dr_me_m;

  // Target speeds
  S.V_sheath_LE_ms = finite_or(P_.V_sheath_LE_factor * Vsw_ms, Vsw_ms);
  S.V_ME_ms        = finite_or(P_.V_ME_factor        * Vsw_ms, Vsw_ms);
  S.V_dn_ms        = Vsw_ms;

  // Apex rc diagnostic
  double n_hat_apex[3] = {e1[0],e1[1],e1[2]};
  double u_apex[3]     = {e1[0],e1[1],e1[2]};
  double rc_apex=1.0, Vsh_n=0.0, dummy=0.0;
  local_oblique_rc(S, u_apex, n_hat_apex, S.r_sh_m, S.r_sh_m, rc_apex, Vsh_n, dummy);
  S.rc = finite_or(rc_apex, 1.0);

  S.inv_dr_sheath = (S.dr_sheath_m > 0.0) ? 1.0 / S.dr_sheath_m : 0.0;
  S.rc_floor      = (P_.sheath_comp_floor > 1.0) ? P_.sheath_comp_floor : 1.0;

  return S;
}

// Directional shock radius and outward normal
void Model::shape_radius_normal(const StepState& S,
                                double ux, double uy, double uz,
                                double& Rdir_m, double n_hat[3]) const
{
  double u[3] = {ux,uy,uz};
  ::safe_normalize(u);

  const double e1[3] = {S.e1[0],S.e1[1],S.e1[2]};
  const double e2[3] = {S.e2[0],S.e2[1],S.e2[2]};
  const double e3[3] = {S.e3[0],S.e3[1],S.e3[2]};

  const double u1 = u[0]*e1[0] + u[1]*e1[1] + u[2]*e1[2];
  const double u2 = u[0]*e2[0] + u[1]*e2[1] + u[2]*e2[2];
  const double u3 = u[0]*e3[0] + u[1]*e3[1] + u[2]*e3[2];

  switch (P_.shape){
    case ShockShape::Sphere: {
      Rdir_m = S.r_sh_m;
      n_hat[0]=u[0]; n_hat[1]=u[1]; n_hat[2]=u[2];
    } break;

    case ShockShape::Ellipsoid: {
      const double a = S.r_sh_m;
      const double b = a * std::max(1e-3, P_.axis_ratio_y);
      const double c = a * std::max(1e-3, P_.axis_ratio_z);
      const double denom = (u1*u1)/(a*a) + (u2*u2)/(b*b) + (u3*u3)/(c*c);
      const double lam   = (denom > 0.0) ? (1.0 / std::sqrt(denom)) : 0.0;
      Rdir_m = lam;

      const double x = Rdir_m * u1, y = Rdir_m * u2, z = Rdir_m * u3;
      double n_loc[3] = { x/(a*a), y/(b*b), z/(c*c) };
      double ng[3] = {
        n_loc[0]*e1[0] + n_loc[1]*e2[0] + n_loc[2]*e3[0],
        n_loc[0]*e1[1] + n_loc[1]*e2[1] + n_loc[2]*e3[1],
        n_loc[0]*e1[2] + n_loc[1]*e2[2] + n_loc[2]*e3[2]
      };
      ::safe_normalize(ng);
      n_hat[0]=ng[0]; n_hat[1]=ng[1]; n_hat[2]=ng[2];
    } break;

    case ShockShape::ConeSSE: {
      // SSE-like cosine cap with flank slowdown exponent
      const double theta = std::acos(std::max(-1.0, std::min(1.0, u1)));
      double R = 0.0;
      if (theta > P_.half_width_rad){
        const double cw = std::cos(P_.half_width_rad);
        R = S.r_sh_m * std::pow(std::max(0.0, cw), std::max(0.0, P_.flank_slowdown_m));
      }else{
        const double c = std::cos(theta);
        R = S.r_sh_m * std::pow(std::max(0.0, c), std::max(0.0, P_.flank_slowdown_m));
      }
      Rdir_m = R;
      n_hat[0]=u[0]; n_hat[1]=u[1]; n_hat[2]=u[2];
    } break;
  }
}

// Oblique fast-mode Mach proxy → rc and Vsh_n
void Model::local_oblique_rc(const StepState& S, const double u[3], const double n_hat[3],
                             double Rdir_m, double r_eval_m,
                             double& rc_out, double& Vsh_n_out, double& thetaBn_out) const
{
  (void)Rdir_m; // not used in this proxy, reserved for variants
  const double Vsw = P_.V_sw_kms * 1e3;

  // Upstream density at r_eval (Leblanc scaled to n(1 AU))
  const double A = 3.3e5, B = 4.1e6, C = 8.0e7; // cm^-3
  const double s = Rs / std::max(r_eval_m, 1.05*Rs);
  const double inv2 = s*s, inv4 = inv2*inv2, inv6 = inv4*inv2;
  const double sAU = Rs / AU;
  const double n1_base = A*sAU*sAU + B*std::pow(sAU,4) + C*std::pow(sAU,6);
  const double leb_scale = (n1_base > 0.0) ? (P_.n1AU_cm3 / n1_base) : 1.0;
  const double n_up_m3 = finite_or((leb_scale * (A*inv2 + B*inv4 + C*inv6)) * 1e6, 1e6);
  const double rho = std::max(1e-12, n_up_m3) * MP;

  // Sound & Alfven speeds
  const double cs  = std::sqrt(std::max(0.0, P_.gamma_ad) * KB * std::max(0.0, P_.T_K) / MP);
  double B_up[3]; ::parker_vec_T(P_, u, r_eval_m, Vsw, B_up);
  const double Bmag = std::sqrt(std::max(0.0, B_up[0]*B_up[0] + B_up[1]*B_up[1] + B_up[2]*B_up[2]));
  const double vA   = (rho>0.0) ? (Bmag / std::sqrt(MU0 * rho)) : 0.0;

  // Obliquity θBn
  const double b_hat[3] = { (Bmag>0.0)? B_up[0]/Bmag : 0.0,
                            (Bmag>0.0)? B_up[1]/Bmag : 0.0,
                            (Bmag>0.0)? B_up[2]/Bmag : 0.0 };
  const double cosBn = std::fabs(b_hat[0]*n_hat[0] + b_hat[1]*n_hat[1] + b_hat[2]*n_hat[2]);
  const double thBn  = std::acos(std::max(-1.0, std::min(1.0, cosBn)));
  thetaBn_out = finite_or(thBn, 0.0);

  // Oblique fast speed
  const double a = vA*vA + cs*cs;
  const double disc = std::max(0.0, a*a - 4.0*cs*cs*vA*vA*cosBn*cosBn);
  const double cf = std::sqrt(0.5 * (a + std::sqrt(disc)));

  // Directional expansion speed projected on local normal
  double Vsh_dir = S.V_sh_ms;
  if (P_.shape == ShockShape::ConeSSE){
    const double u1 = u[0]*S.e1[0] + u[1]*S.e1[1] + u[2]*S.e1[2];
    Vsh_dir = S.V_sh_ms * std::pow(std::max(0.0, u1), std::max(0.0, P_.flank_slowdown_m));
  }
  Vsh_n_out = finite_or( Vsh_dir * (n_hat[0]*u[0] + n_hat[1]*u[1] + n_hat[2]*u[2]) , 0.0);

  const double Vsw_n = Vsw * (u[0]*n_hat[0] + u[1]*n_hat[1] + u[2]*n_hat[2]);
  const double U1n   = std::max(0.0, Vsh_n_out - Vsw_n);
  const double Mfn   = (cf > 0.0) ? (U1n / cf) : 0.0;

  double rc = 1.0;
  if (Mfn > 1.0){
    const double g  = std::max(1.01, P_.gamma_ad);
    const double M2 = Mfn*Mfn;
    rc = ((g+1.0)*M2) / ((g-1.0)*M2 + 2.0);
    if (rc > 4.0) rc = 4.0; // γ=5/3 limit
  }
  rc_out = finite_or(rc, 1.0);
}

// n and V (no B)
void Model::evaluate_cartesian_fast(const StepState& S,
                                    const double* x_m, const double* y_m, const double* z_m,
                                    double* n_m3, double* Vx_ms, double* Vy_ms, double* Vz_ms,
                                    std::size_t N) const
{
  const double Vup = P_.V_sw_kms * 1e3;
  const double scaleR = S.r_sh_m / AU;
  const double dr_sheath = P_.sheath_thick_AU_at1AU * scaleR * AU;
  const double dr_me     = P_.ejecta_thick_AU_at1AU * scaleR * AU;

  const double w_sh = S.w_shock_m;
  const double w_le = S.w_le_m;
  const double w_te = S.w_te_m;

  const double Vshe_LE = S.V_sheath_LE_ms;
  const double Vme     = S.V_ME_ms;

  const double A = 3.3e5, B = 4.1e6, C = 8.0e7;
  const double sAU = Rs / AU;
  const double n1_base = A*sAU*sAU + B*std::pow(sAU,4) + C*std::pow(sAU,6);
  const double leb_scale = (n1_base > 0.0) ? (P_.n1AU_cm3 / n1_base) : 1.0;

  const double p = (P_.sheath_ramp_power < 1.0) ? 1.0 : P_.sheath_ramp_power;
  const double rc_floor = (P_.sheath_comp_floor < 1.0) ? 1.0 : P_.sheath_comp_floor;

  auto edge = [&](double rnow, double r0, double w)->double {
    if (w <= 0.0) return (rnow <= r0) ? 1.0 : 0.0;
    return ::smoothstep01(0.5 + (r0 - rnow)/(2.0*w));
  };

  for (std::size_t i=0; i<N; ++i){
    const double x=x_m[i], y=y_m[i], z=z_m[i];
    const double r2 = x*x + y*y + z*z;
    const double r  = std::sqrt(std::max(1e-12, r2));
    const double invr = 1.0 / r;
    double u[3] = {x*invr, y*invr, z*invr};

    double Rdir=0.0, n_hat[3]={0,0,1};
    shape_radius_normal(S, u[0], u[1], u[2], Rdir, n_hat);

    // Upstream n at r
    const double s = Rs / r;
    const double inv2 = s*s, inv4 = inv2*inv2, inv6 = inv4*inv2;
    const double n_up = finite_or((leb_scale * (A*inv2 + B*inv4 + C*inv6)) * 1e6, 1e6);

    double rc_loc=1.0, Vsh_n=0.0, dummy=0.0;
    local_oblique_rc(S, u, n_hat, Rdir, (r > Rdir ? r : Rdir), rc_loc, Vsh_n, dummy);

    const double xi   = (dr_sheath > 0.0) ? clamp01((Rdir - r) / dr_sheath) : 0.0;
    const double ramp = (p == 1.0) ? (1.0 - xi) : std::pow(1.0 - xi, p);
    const double Csheath = rc_floor + (rc_loc - rc_floor) * ramp;
    const double n_sheath = Csheath * n_up;
    const double V_sheath = Vup + (Vshe_LE - Vup) * xi;

    const double n_ejecta = (P_.f_ME > 0.0 ? P_.f_ME : 0.0) * n_up;
    const double V_ejecta = Vme;

    double a = (rc_loc > 1.0) ? edge(r, Rdir, w_sh) : 0.0;
    double n_mix = (1.0 - a)*n_up + a*n_sheath;
    double Vmag  = (1.0 - a)*Vup  + a*V_sheath;

    a = (rc_loc > 1.0) ? edge(r, Rdir - dr_sheath, w_le) : 0.0;
    n_mix = (1.0 - a)*n_mix + a*n_ejecta;
    Vmag  = (1.0 - a)*Vmag + a*V_ejecta;

    a = (rc_loc > 1.0) ? edge(r, Rdir - dr_sheath - dr_me, w_te) : 0.0;
    n_mix = (1.0 - a)*n_mix + a*n_up;
    Vmag  = (1.0 - a)*Vmag + a*Vup;

    n_m3[i]  = finite_or(n_mix, n_up);
    Vx_ms[i] = finite_or(Vmag * u[0], 0.0);
    Vy_ms[i] = finite_or(Vmag * u[1], 0.0);
    Vz_ms[i] = finite_or(Vmag * u[2], 0.0);
  }
}

// n, V, B with tangential amplification in sheath
void Model::evaluate_cartesian_with_B(const StepState& S,
                                      const double* x_m, const double* y_m, const double* z_m,
                                      double* n_m3, double* Vx_ms, double* Vy_ms, double* Vz_ms,
                                      double* Bx_T, double* By_T, double* Bz_T,
                                      std::size_t N) const
{
  const double Vup = P_.V_sw_kms * 1e3;

  const double scaleR    = S.r_sh_m / AU;
  const double dr_sheath = P_.sheath_thick_AU_at1AU * scaleR * AU;
  const double dr_me     = P_.ejecta_thick_AU_at1AU * scaleR * AU;
  const double w_sh      = S.w_shock_m;
  const double w_le      = S.w_le_m;
  const double w_te      = S.w_te_m;

  const double Vshe_LE = S.V_sheath_LE_ms;
  const double Vme     = S.V_ME_ms;

  const double A = 3.3e5, B = 4.1e6, C = 8.0e7;
  const double sAU = Rs / AU;
  const double n1_base   = A*sAU*sAU + B*std::pow(sAU,4) + C*std::pow(sAU,6);
  const double leb_scale = (n1_base > 0.0) ? (P_.n1AU_cm3 / n1_base) : 1.0;

  const double p = (P_.sheath_ramp_power < 1.0) ? 1.0 : P_.sheath_ramp_power;
  const double rc_floor = (P_.sheath_comp_floor < 1.0) ? 1.0 : P_.sheath_comp_floor;

  auto edge = [&](double rnow, double r0, double w)->double {
    if (w <= 0.0) return (rnow <= r0) ? 1.0 : 0.0;
    return ::smoothstep01(0.5 + (r0 - rnow)/(2.0*w));
  };

  for (std::size_t i=0; i<N; ++i){
    const double x=x_m[i], y=y_m[i], z=z_m[i];
    const double r2 = x*x + y*y + z*z;
    const double r  = std::sqrt(std::max(1e-12, r2));
    const double invr = 1.0 / r;
    double u[3] = {x*invr, y*invr, z*invr};

    double Rdir=0.0, n_hat[3]={0,0,1};
    shape_radius_normal(S, u[0], u[1], u[2], Rdir, n_hat);

    // Upstream n and B at r
    const double s = Rs / r;
    const double inv2 = s*s, inv4 = inv2*inv2, inv6 = inv4*inv2;
    const double n_up = finite_or((leb_scale * (A*inv2 + B*inv4 + C*inv6)) * 1e6, 1e6);
    double B_up[3]; ::parker_vec_T(P_, u, r, Vup, B_up);

    double rc_loc=1.0, Vsh_n=0.0, dummy=0.0;
    local_oblique_rc(S, u, n_hat, Rdir, (r > Rdir ? r : Rdir), rc_loc, Vsh_n, dummy);

    const double xi   = (dr_sheath > 0.0) ? clamp01((Rdir - r) / dr_sheath) : 0.0;
    const double ramp = (p == 1.0) ? (1.0 - xi) : std::pow(1.0 - xi, p);
    const double Csheath = rc_floor + (rc_loc - rc_floor) * ramp;
    const double n_sheath = Csheath * n_up;
    const double V_sheath = Vup + (Vshe_LE - Vup) * xi;

    const double n_ejecta = (P_.f_ME > 0.0 ? P_.f_ME : 0.0) * n_up;
    const double V_ejecta = Vme;

    double a = (rc_loc > 1.0) ? edge(r, Rdir, w_sh) : 0.0;
    double n_mix = (1.0 - a)*n_up + a*n_sheath;
    double Vmag  = (1.0 - a)*Vup  + a*V_sheath;

    a = (rc_loc > 1.0) ? edge(r, Rdir - dr_sheath, w_le) : 0.0;
    n_mix = (1.0 - a)*n_mix + a*n_ejecta;
    Vmag  = (1.0 - a)*Vmag + a*V_ejecta;

    a = (rc_loc > 1.0) ? edge(r, Rdir - dr_sheath - dr_me, w_te) : 0.0;
    n_mix = (1.0 - a)*n_mix + a*n_up;
    Vmag  = (1.0 - a)*Vmag + a*Vup;

    // B: decompose wrt local normal, amplify Bt in sheath toward rc
    const double Bn_mag = B_up[0]*n_hat[0] + B_up[1]*n_hat[1] + B_up[2]*n_hat[2];
    double Bn[3] = {Bn_mag*n_hat[0], Bn_mag*n_hat[1], Bn_mag*n_hat[2]};
    double Bt[3] = {B_up[0]-Bn[0],   B_up[1]-Bn[1],   B_up[2]-Bn[2]  };
    const double a_shock = (rc_loc > 1.0) ? edge(r, Rdir, w_sh) : 0.0; // 0 upstream → 1 just inside
    const double amp_t   = 1.0 + (rc_loc - 1.0) * a_shock * ramp;
    double Bv[3] = { Bn[0] + amp_t*Bt[0],
                     Bn[1] + amp_t*Bt[1],
                     Bn[2] + amp_t*Bt[2] };

    n_m3[i]  = finite_or(n_mix, n_up);
    Vx_ms[i] = finite_or(Vmag * u[0], 0.0);
    Vy_ms[i] = finite_or(Vmag * u[1], 0.0);
    Vz_ms[i] = finite_or(Vmag * u[2], 0.0);
    Bx_T[i]  = finite_or(Bv[0], 0.0);
    By_T[i]  = finite_or(Bv[1], 0.0);
    Bz_T[i]  = finite_or(Bv[2], 0.0);
  }
}

// Convenience: n, V, B, and div(Vsw)
void Model::evaluate_cartesian_with_B_div(const StepState& S,
                                     const double* x_m, const double* y_m, const double* z_m,
                                     double* n_m3, double* Vx_ms, double* Vy_ms, double* Vz_ms,
                                     double* Bx_T, double* By_T, double* Bz_T,
                                     double* divVsw,
                                     std::size_t N,
                                     double dr_frac) const
{
  evaluate_cartesian_with_B(S, x_m, y_m, z_m,
                            n_m3, Vx_ms, Vy_ms, Vz_ms,
                            Bx_T, By_T, Bz_T, N);
  compute_divV_radial(S, x_m, y_m, z_m, divVsw, N, dr_frac);
}

// Radial (centered) ∇·Vsw estimate with guards
void Model::compute_divV_radial(const StepState& S,
                                const double* x_m, const double* y_m, const double* z_m,
                                double* divV, std::size_t N, double dr_frac) const
{
  const double rmin   = 1.05 * Rs;      // avoid inner boundary
  const double dr_min = 1.0e-4 * AU;    // absolute probe step

  for (std::size_t i=0; i<N; ++i){
    const double x = x_m[i], y = y_m[i], z = z_m[i];
    const double r2 = x*x + y*y + z*z;
    const double r  = std::sqrt(std::max(1e-12, r2));
    const double invr = 1.0 / r;
    const double u[3] = { x*invr, y*invr, z*invr };

    const double dr = std::max(dr_min, dr_frac * r);
    const double rp = std::max(rmin, r + dr);
    const double rm = std::max(rmin, r - dr);
    const double denom_r = (rp > rm) ? (rp - rm) : std::max(dr_min, std::abs(dr));

    double xp = rp*u[0], yp = rp*u[1], zp = rp*u[2];
    double xm = rm*u[0], ym = rm*u[1], zm = rm*u[2];

    double ntmp, Vxp,Vyp,Vzp, Vxm,Vym,Vzm;
    evaluate_cartesian_fast(S, &xp,&yp,&zp, &ntmp, &Vxp,&Vyp,&Vzp, 1);
    evaluate_cartesian_fast(S, &xm,&ym,&zm, &ntmp, &Vxm,&Vym,&Vzm, 1);

    const double Vrp = Vxp*u[0] + Vyp*u[1] + Vzp*u[2];
    const double Vrm = Vxm*u[0] + Vym*u[1] + Vzm*u[2];

    const double num = ( (rp*rp)*Vrp - (rm*rm)*Vrm ) / denom_r;
    const double div_val = num / (r*r);

    divV[i] = finite_or(div_val, 0.0);
  }
}

// Directional diagnostic
void Model::diagnose_direction(const StepState& S, const double u[3],
                               double& Rdir_m, double n_hat[3],
                               double& rc_loc, double& Vsh_n) const
{
  shape_radius_normal(S, u[0], u[1], u[2], Rdir_m, n_hat);
  double dummy=0.0;
  local_oblique_rc(S, u, n_hat, Rdir_m, Rdir_m, rc_loc, Vsh_n, dummy);
}

// Lat–lon triangulation over cap (for ConeSSE: θ ≤ half_width)
ShockMesh Model::build_shock_mesh(const StepState& S,
                                  std::size_t nTheta, std::size_t nPhi) const
{
  ShockMesh M;
  if (nTheta < 3) nTheta = 3;
  if (nPhi   < 3) nPhi   = 3;

  const double thetaMax = (P_.shape==ShockShape::ConeSSE) ? P_.half_width_rad : PI;

  // Vertices
  for (std::size_t it=0; it<=nTheta; ++it){
    const double t  = thetaMax * (double(it) / double(nTheta));
    const double ct = std::cos(t), st = std::sin(t);
    for (std::size_t ip=0; ip<=nPhi; ++ip){
      const double p  = 2.0*PI * (double(ip) / double(nPhi));
      const double cp = std::cos(p), sp = std::sin(p);

      // Direction in apex frame → global
      double u_loc[3] = { ct, st*cp, st*sp };
      double u[3] = {
        u_loc[0]*S.e1[0] + u_loc[1]*S.e2[0] + u_loc[2]*S.e3[0],
        u_loc[0]*S.e1[1] + u_loc[1]*S.e2[1] + u_loc[2]*S.e3[1],
        u_loc[0]*S.e1[2] + u_loc[1]*S.e2[2] + u_loc[2]*S.e3[2]
      };
      ::safe_normalize(u);

      double Rdir=0.0, n_hat[3]={0,0,1};
      shape_radius_normal(S, u[0], u[1], u[2], Rdir, n_hat);

      M.x.push_back(finite_or(Rdir*u[0],0.0));
      M.y.push_back(finite_or(Rdir*u[1],0.0));
      M.z.push_back(finite_or(Rdir*u[2],0.0));

      M.n_hat_x.push_back(finite_or(n_hat[0],0.0));
      M.n_hat_y.push_back(finite_or(n_hat[1],0.0));
      M.n_hat_z.push_back(finite_or(n_hat[2],1.0));

      double rc_loc=1.0, Vsh_n=0.0, dummy=0.0;
      local_oblique_rc(S, u, n_hat, Rdir, Rdir, rc_loc, Vsh_n, dummy);
      M.rc.push_back(finite_or(rc_loc,1.0));
      M.Vsh_n.push_back(finite_or(Vsh_n,0.0));
    }
  }

  // Connectivity (two triangles per cell)
  const std::size_t NvPhi   = nPhi + 1;
  auto idx = [&](std::size_t it, std::size_t ip){ return int(it*NvPhi + ip); };
  for (std::size_t it=0; it<nTheta; ++it){
    for (std::size_t ip=0; ip<nPhi; ++ip){
      int i00 = idx(it,   ip);
      int i01 = idx(it,   ip+1);
      int i10 = idx(it+1, ip);
      int i11 = idx(it+1, ip+1);
      M.tri_i.push_back(i00+1); M.tri_j.push_back(i10+1); M.tri_k.push_back(i11+1);
      M.tri_i.push_back(i00+1); M.tri_j.push_back(i11+1); M.tri_k.push_back(i01+1);
    }
  }
  return M;
}

// Per-triangle metrics (area, centroid, unit normal, mean rc, mean Vsh_n)
void Model::compute_triangle_metrics(const ShockMesh& M, TriMetrics& T) const
{
  const std::size_t Ne = M.tri_i.size();
  T.area.resize(Ne); T.nx.resize(Ne); T.ny.resize(Ne); T.nz.resize(Ne);
  T.cx.resize(Ne);   T.cy.resize(Ne); T.cz.resize(Ne);
  T.rc_mean.resize(Ne); T.Vsh_n_mean.resize(Ne);

  for (std::size_t e=0; e<Ne; ++e){
    int ia = M.tri_i[e]-1, ib = M.tri_j[e]-1, ic = M.tri_k[e]-1;
    double Ax=M.x[ia], Ay=M.y[ia], Az=M.z[ia];
    double Bx=M.x[ib], By=M.y[ib], Bz=M.z[ib];
    double Cx=M.x[ic], Cy=M.y[ic], Cz=M.z[ic];

    double ABx=Bx-Ax, ABy=By-Ay, ABz=Bz-Az;
    double ACx=Cx-Ax, ACy=Cy-Ay, ACz=Cz-Az;

    double nx = ABy*ACz - ABz*ACy;
    double ny = ABz*ACx - ABx*ACz;
    double nz = ABx*ACy - ABy*ACx;
    double twiceA = std::sqrt(std::max(0.0, nx*nx + ny*ny + nz*nz));
    double area = 0.5 * twiceA;
    if (twiceA > 0.0){ nx/=twiceA; ny/=twiceA; nz/=twiceA; }

    T.area[e]=finite_or(area,0.0);
    T.nx[e]=finite_or(nx,0.0); T.ny[e]=finite_or(ny,0.0); T.nz[e]=finite_or(nz,1.0);
    T.cx[e]=finite_or((Ax+Bx+Cx)/3.0,0.0);
    T.cy[e]=finite_or((Ay+By+Cy)/3.0,0.0);
    T.cz[e]=finite_or((Az+Bz+Cz)/3.0,0.0);

    T.rc_mean[e]    = finite_or((M.rc[ia] + M.rc[ib] + M.rc[ic]) / 3.0, 1.0);
    T.Vsh_n_mean[e] = finite_or((M.Vsh_n[ia] + M.Vsh_n[ib] + M.Vsh_n[ic]) / 3.0, 0.0);
  }
}

// FE surface with nodal & cell metrics (BLOCK + connectivity)
bool Model::write_shock_surface_center_metrics_tecplot(
    const ShockMesh& M, const TriMetrics& T, const char* path) const
{
  const std::size_t Nv = M.x.size();
  const std::size_t Ne = M.tri_i.size();

  std::FILE* fp = std::fopen(path,"w");
  if(!fp) return false;

  std::fprintf(fp, "TITLE = \"Shock surface (cell-centered metrics + nodal rc)\"\n");
  std::fprintf(fp, "VARIABLES = "
                   "\"X\",\"Y\",\"Z\",\"rc\",\"Vsh_n\","
                   "\"area\",\"rc_mean\",\"Vsh_n_mean\",\"tnx\",\"tny\",\"tnz\",\"cx\",\"cy\",\"cz\"\n");
  std::fprintf(fp, "ZONE T=\"surface_cells\", N=%zu, E=%zu, ZONETYPE=FETRIANGLE, DATAPACKING=BLOCK,\n", Nv, Ne);
  std::fprintf(fp, "VARLOCATION=([1-5]=NODAL, [6-14]=CELLCENTERED)\n");

  const int PER_LINE=8;
  auto dumpN = [&](const std::vector<double>& a){
    int cnt=0; for(double v: a){ std::fprintf(fp,"%.9e ",finite_or(v)); if(++cnt==PER_LINE){std::fprintf(fp,"\n"); cnt=0;} }
    if (cnt) std::fprintf(fp,"\n");
  };
  auto dumpE = dumpN;

  dumpN(M.x); dumpN(M.y); dumpN(M.z); dumpN(M.rc); dumpN(M.Vsh_n);
  dumpE(T.area); dumpE(T.rc_mean); dumpE(T.Vsh_n_mean);
  dumpE(T.nx); dumpE(T.ny); dumpE(T.nz);
  dumpE(T.cx); dumpE(T.cy); dumpE(T.cz);

  for (std::size_t e=0; e<Ne; ++e)
    std::fprintf(fp, "%d %d %d\n", M.tri_i[e], M.tri_j[e], M.tri_k[e]);

  std::fclose(fp);
  return true;
}

// Outward-shifted apex box (axis-aligned globally)
BoxSpec Model::default_apex_box(const StepState& S, double half_AU, int N) const
{
  BoxSpec B;
  const double h = half_AU * AU;
  B.hx=h; B.hy=h; B.hz=h;
  const double OUTSIDE_BIAS = 0.4;         // place box partly outside shock
  const double shift = OUTSIDE_BIAS * h;
  B.cx = S.a_m*S.e1[0] + shift*S.e1[0];
  B.cy = S.a_m*S.e1[1] + shift*S.e1[1];
  B.cz = S.a_m*S.e1[2] + shift*S.e1[2];
  B.Ni=N; B.Nj=N; B.Nk=N;
  return B;
}

// Combined dataset bundle:
//  1) surface_nodal  (FEPOINT, ET=TRIANGLE)   — with connectivity
//  2) surface_cells  (FETRIANGLE, BLOCK)      — with connectivity
//  3) volume_box     (STRUCTURED POINT)       — no connectivity
//  4) box_face_minX  (STRUCTURED POINT, 2-D)  — optional; no connectivity
bool Model::write_tecplot_dataset_bundle(const ShockMesh& M, const TriMetrics& T,
                                         const StepState& S, const BoxSpec& B,
                                         const char* path) const
{
  const std::size_t Nv = M.x.size();
  const std::size_t Ne = M.tri_i.size();

  std::FILE* fp = std::fopen(path,"w");
  if(!fp) return false;
  auto p=[&](const char* fmt, auto... args){ std::fprintf(fp, fmt, args...); };

  // Keep variables identical across all zones in this dataset
  p("TITLE = \"SW+CME dataset\"\n");
  p("VARIABLES = "
    "\"X\",\"Y\",\"Z\","
    "\"n\",\"Vx\",\"Vy\",\"Vz\","
    "\"Bx\",\"By\",\"Bz\",\"divVsw\","
    "\"rc\",\"Vsh_n\","
    "\"nx\",\"ny\",\"nz\",\"area\",\"rc_mean\",\"Vsh_n_mean\","
    "\"tnx\",\"tny\",\"tnz\",\"cx\",\"cy\",\"cz\"\n");

  // ---- Zone 1: surface_nodal (FEPOINT) ----
  p("ZONE T=\"surface_nodal\", N=%zu, E=%zu, F=FEPOINT, ET=TRIANGLE\n", Nv, Ne);
  for (std::size_t i=0;i<Nv;++i){
    p("%.9e %.9e %.9e "
      "%.9e %.9e %.9e %.9e "
      "%.9e %.9e %.9e %.9e "
      "%.6e %.6e "
      "%.6e %.6e %.6e "
      "%.9e %.9e %.9e "
      "%.6e %.6e %.6e "
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

  // ---- Zone 2: surface_cells (BLOCK) ----
  p("ZONE T=\"surface_cells\", N=%zu, E=%zu, ZONETYPE=FETRIANGLE, DATAPACKING=BLOCK,\n", Nv, Ne);
  p("VARLOCATION=([1-3,12-13]=NODAL, [4-11,14-25]=CELLCENTERED)\n");

  const int PER_LINE=8;
  auto dumpN_wrap = [&](const std::vector<double>& a){
    int cnt=0; for(double v: a){
      std::fprintf(fp,"%.9e ", finite_or(v));
      if(++cnt==PER_LINE){ std::fprintf(fp,"\n"); cnt=0; }
    }
    if (cnt) std::fprintf(fp,"\n");
  };
  auto dumpE_wrap = dumpN_wrap;
  auto dumpEzeros_wrap = [&](std::size_t count){
    int cnt=0; for (std::size_t e=0;e<count;++e){
      std::fprintf(fp,"0 ");
      if(++cnt==PER_LINE){ std::fprintf(fp,"\n"); cnt=0; }
    }
    if (cnt) std::fprintf(fp,"\n");
  };

  // BLOCK order: 25 variables
  dumpN_wrap(M.x); dumpN_wrap(M.y); dumpN_wrap(M.z);              // 1..3 (NODAL)
  dumpEzeros_wrap(Ne); dumpEzeros_wrap(Ne); dumpEzeros_wrap(Ne);  // 4..6 n,Vx,Vy
  dumpEzeros_wrap(Ne);                                            // 7   Vz
  dumpEzeros_wrap(Ne); dumpEzeros_wrap(Ne); dumpEzeros_wrap(Ne);  // 8..10 Bx,By,Bz
  dumpEzeros_wrap(Ne);                                            // 11  divVsw
  dumpN_wrap(M.rc); dumpN_wrap(M.Vsh_n);                          // 12..13 (NODAL)
  dumpE_wrap(T.nx); dumpE_wrap(T.ny); dumpE_wrap(T.nz);           // 14..16
  dumpE_wrap(T.area); dumpE_wrap(T.rc_mean); dumpE_wrap(T.Vsh_n_mean); // 17..19
  dumpEzeros_wrap(Ne); dumpEzeros_wrap(Ne); dumpEzeros_wrap(Ne);  // 20..22 tnx,tny,tnz
  dumpE_wrap(T.cx); dumpE_wrap(T.cy); dumpE_wrap(T.cz);           // 23..25

  for (std::size_t e=0;e<Ne;++e)
    p("%d %d %d\n", M.tri_i[e], M.tri_j[e], M.tri_k[e]);

  // ---- Zone 3: volume_box (STRUCTURED POINT; NO CONNECTIVITY) ----
  p("ZONE T=\"volume_box\", I=%d, J=%d, K=%d, DATAPACKING=POINT\n", B.Ni,B.Nj,B.Nk);
  for (int kk=0; kk<B.Nk; ++kk){
    const double zk = B.cz + (-B.hz + (2.0*B.hz) * (kk / double(std::max(1,B.Nk-1))));
    for (int jj=0; jj<B.Nj; ++jj){
      const double yj = B.cy + (-B.hy + (2.0*B.hy) * (jj / double(std::max(1,B.Nj-1))));
      for (int ii=0; ii<B.Ni; ++ii){
        const double xi = B.cx + (-B.hx + (2.0*B.hx) * (ii / double(std::max(1,B.Ni-1))));

        double n,Vx,Vy,Vz,Bx,By,Bz,div;
        evaluate_cartesian_with_B(S, &xi,&yj,&zk, &n,&Vx,&Vy,&Vz, &Bx,&By,&Bz, 1);
        compute_divV_radial(S, &xi,&yj,&zk, &div, 1, 1e-3);

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
          finite_or(Bx), finite_or(By), finite_or(Bz), finite_or(div, 0.0),
          0.0,0.0,         // rc, Vsh_n
          0.0,0.0,0.0,     // nx, ny, nz
          0.0,0.0,0.0,     // area, rc_mean, Vsh_n_mean
          0.0,0.0,0.0,     // tnx, tny, tnz
          0.0,0.0,0.0      // cx, cy, cz
        );
      }
    }
  }

  // ---- Zone 4: optional structured 2-D min-X face (NO CONNECTIVITY) ----
  if (ADD_FACE_ZONE_IN_BUNDLE){
    const int I = std::max(2, B.Nj);  // sample along Y
    const int J = std::max(2, B.Nk);  // sample along Z
    const double x0 = B.cx - B.hx;
    const double y0 = B.cy - B.hy, y1 = B.cy + B.hy;
    const double z0 = B.cz - B.hz, z1 = B.cz + B.hz;

    p("ZONE T=\"box_face_minX\", I=%d, J=%d, DATAPACKING=POINT\n", I, J);
    for (int j=0; j<J; ++j){
      const double tz = (J==1) ? 0.0 : double(j)/double(J-1);
      const double z  = z0 + (z1 - z0) * tz;
      for (int i=0; i<I; ++i){
        const double ty = (I==1) ? 0.0 : double(i)/double(I-1);
        const double y  = y0 + (y1 - y0) * ty;

        const double x = x0;
        double n,Vx,Vy,Vz,Bx,By,Bz,div;
        evaluate_cartesian_with_B(S, &x,&y,&z, &n,&Vx,&Vy,&Vz, &Bx,&By,&Bz, 1);
        compute_divV_radial(S, &x,&y,&z, &div, 1, 1e-3);

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
          finite_or(Bx), finite_or(By), finite_or(Bz), finite_or(div, 0.0),
          0.0,0.0,
          0.0,0.0,0.0,
          0.0,0.0,0.0,
          0.0,0.0,0.0,
          0.0,0.0,0.0
        );
      }
    }
  }

  std::fclose(fp);
  return true;
}

// Standalone structured 2-D min-X face writer (no connectivity required)
bool Model::write_box_face_minX_tecplot_structured(const StepState& S,
                                                   const BoxSpec& B,
                                                   const char* path) const
{
  const int I = std::max(2, B.Nj);
  const int J = std::max(2, B.Nk);

  const double x0 = B.cx - B.hx;
  const double y0 = B.cy - B.hy, y1 = B.cy + B.hy;
  const double z0 = B.cz - B.hz, z1 = B.cz + B.hz;

  std::FILE* fp = std::fopen(path, "w");
  if (!fp) return false;
  auto p=[&](const char* fmt, auto... args){ std::fprintf(fp, fmt, args...); };

  p("TITLE = \"Box face (minX)\"\n");
  p("VARIABLES = "
    "\"X\",\"Y\",\"Z\","
    "\"n\",\"Vx\",\"Vy\",\"Vz\","
    "\"Bx\",\"By\",\"Bz\",\"divVsw\","
    "\"rc\",\"Vsh_n\","
    "\"nx\",\"ny\",\"nz\",\"area\",\"rc_mean\",\"Vsh_n_mean\","
    "\"tnx\",\"tny\",\"tnz\",\"cx\",\"cy\",\"cz\"\n");
  p("ZONE T=\"box_face_minX\", I=%d, J=%d, DATAPACKING=POINT\n", I, J);

  for (int j=0; j<J; ++j){
    const double tz = (J==1) ? 0.0 : double(j)/double(J-1);
    const double z  = z0 + (z1 - z0) * tz;
    for (int i=0; i<I; ++i){
      const double ty = (I==1) ? 0.0 : double(i)/double(I-1);
      const double y  = y0 + (y1 - y0) * ty;

      const double x = x0;

      double n,Vx,Vy,Vz,Bx,By,Bz,div;
      evaluate_cartesian_with_B(S, &x,&y,&z, &n,&Vx,&Vy,&Vz, &Bx,&By,&Bz, 1);
      compute_divV_radial(S, &x,&y,&z, &div, 1, 1e-3);

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
        finite_or(Bx), finite_or(By), finite_or(Bz), finite_or(div, 0.0),
        0.0,0.0,
        0.0,0.0,0.0,
        0.0,0.0,0.0,
        0.0,0.0,0.0,
        0.0,0.0,0.0
      );
    }
  }

  std::fclose(fp);
  return true;
}

} // namespace swcme3d

