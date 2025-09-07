#ifndef SWCME1D_HPP
#define SWCME1D_HPP
/*
===============================================================================
 swcme1d.hpp — Header-only 1-D Solar Wind + CME (DBM) model
-------------------------------------------------------------------------------
 PURPOSE
   Provide a lightweight, numerically efficient 1-D model of the heliocentric
   solar wind plus a driven CME forward shock and its downstream structure
   (sheath and magnetic ejecta, ME). Designed to:
     • return n(r), V(r) and Parker B(r) with CME-induced modifications,
     • expose ∇·V (for SEP adiabatic cooling/heating),
     • offer narrow, independent smoothing at the three edges
       (shock, sheath→ME leading edge, ME→ambient trailing edge),
     • be fast enough for per-particle queries in SEP solvers.

 PHYSICS SUMMARY (compact)
   • Ambient wind: steady, radial, speed V_sw (km/s), density n(r) follows
     a Leblanc (1998)-shaped r^(-2,-4,-6) profile normalized to n(1 AU).
   • Magnetic field: Parker spiral
        Br(r)   = Br(1 AU) (AU/r)^2
        Bphi(r) = -Br(r) (Ω r sinθ / V_sw)
        |B|(r)  = |Br(r)| sqrt(1 + (k r_AU)^2),  k = Ω AU sinθ / V_sw
     with |B|(1 AU) set to B1AU; we solve Br(1 AU) = B1AU / sqrt(1+k^2).
   • CME apex kinematics (DBM): with drag Γ and u(t)=V_sh−V_sw,
        u(t)   = u0 / (1 + Γ u0 t),                     u0=V0_sh−V_sw
        R_sh(t)= r0 + V_sw t + [ ln(1 + Γ u0 t) ] / Γ
        V_sh(t)= V_sw + u(t)
     (Here Γ is in SI 1/m; user supplies in 1/km and it is converted.)
   • Shock compression ratio rc from a fast-mode Mach proxy:
        c_s = sqrt(γ k_B T / m_p), v_A = B/√(μ0 ρ),
        c_f = sqrt(c_s^2 + v_A^2),  M_f ≈ max( (V_sh−V_sw)/c_f , 1 )
        rc  = ((γ+1) M_f^2) / ( (γ-1) M_f^2 + 2 ), capped ≤ 4 (γ=5/3).
   • Downstream structure:
       – Sheath (R_LE < r < R_sh): n decays from rc·n_up at shock to ~n_up at LE;
         V goes from V_dn (shock) to V_sheath_LE_factor · V_sw at LE;
         B amplification applied primarily to tangential component (Bphi) and
         decays from rc (at shock) to ~1 at LE.
       – Magnetic ejecta (R_TE < r < R_LE): n = f_ME n_up; V = V_ME_factor V_sw.
       – Each of the three edges uses its own C^1 smoothing width (shock/LE/TE),
         scaled ∝ R_sh (self-similar with distance).

// ----------------------------------------------------------------------------
// PHYSICS OVERVIEW (1-D ALONG A HELIOCENTRIC RAY; ORIGIN = SUN CENTER)
// --------------------------------------------------------------------
// Upstream density n(r): Leblanc, Dulk & Bougeret (1998), Solar Phys. 183, 165
//   n[r] ~ A (Rs/r)^2 + B (Rs/r)^4 + C (Rs/r)^6   [cm^-3]  with
//   A=3.3e5, B=4.1e6, C=8.0e7. We scale these to match a user-given n(1 AU)
//   and convert to SI [m^-3].  Implementation detail for speed:
//     n(r) = C2 * r^-2 + C4 * r^-4 + C6 * r^-6   [m^-3],
//   where C2,C4,C6 (SI) are cached in StepState and r^-k use fused multiplies.
//
// Upstream magnetic field B(r): Parker (1958), ApJ 128, 664
//   In equatorial approximation (fixed sinθ), with solar rotation Ω and wind Vsw:
//   Br ∝ r^-2, Bφ = -Br * (Ω r sinθ / Vsw).  We choose Br(1AU) so that
//   |B|(1 AU) equals user-given B1AU.  Implementation caches Br1AU_T and
//   k_AU = Ω AU sinθ / Vsw, so evaluation is a few mults per sample.
//
// CME apex kinematics: Drag-Based Model (DBM): Vršnak & Žic (2007); Vršnak et al. (2013)
//   u(t) = Vsh − Vsw.  With drag Γ,
//     u(t) = u0 / (1 + Γ u0 t),   r(t) = r0 + Vsw t + (ln(1+Γ u0 t))/Γ.
//   We guard logs/denominators and convert Γ from km^-1 to m^-1.
//
// Regions and smoothing (self-similar):
//   upstream → [shock at R_sh] → sheath → [leading edge at R_LE] → magnetic ejecta
//   → [trailing edge at R_TE] → downstream ambient.
//   Each interface uses a C¹ smoothstep s(x)=x^2(3−2x) over width w; widths scale
//   ∝ R_sh (so they grow with distance). Sheath density uses a ramp that is steeper
//   near the shock (power p≥1) and bounded below by a floor (rc_floor≥1).
//
// Oblique-MHD shock proxy (1-D quasi-radial reduction): Edmiston & Kennel (1984);
// Priest (2014, CUP): We approximate the normal Mach number using the fast speed
//   c_f = sqrt(c_s^2 + v_A^2) with upstream sound speed c_s and Alfven speed v_A,
//   then take a hydrodynamic-like compression
//     rc = ((γ+1) M_n^2) / ((γ−1) M_n^2 + 2),   1 ≤ rc ≤ 4.
//   Downstream speed at the shock follows a continuity proxy.
//   This rc is used as: sheath density jump and a proxy for Bt amplification.
//
// Divergence of bulk flow:
//   ∇·V = (1/r^2) d/dr ( r^2 V_r ). We compute it with a centered FD using
//   r±dr along the same ray; dr = max(1e-4 AU, 1e-3 r).  Numerically robust.
//

// PHYSICS OVERVIEW (1-D ALONG A HELIOCENTRIC RAY; ORIGIN = SUN CENTER)
// --------------------------------------------------------------------
// Upstream density n(r): Leblanc, Dulk & Bougeret (1998), Solar Phys. 183, 165
//   n[r] ~ A (Rs/r)^2 + B (Rs/r)^4 + C (Rs/r)^6   [cm^-3]  with
//   A=3.3e5, B=4.1e6, C=8.0e7. We scale these to match a user-given n(1 AU)
//   and convert to SI [m^-3].  Implementation detail for speed:
//     n(r) = C2 * r^-2 + C4 * r^-4 + C6 * r^-6   [m^-3],
//   where C2,C4,C6 (SI) are cached in StepState and r^-k use fused multiplies.
//
// Upstream magnetic field B(r): Parker (1958), ApJ 128, 664
//   In equatorial approximation (fixed sinθ), with solar rotation Ω and wind Vsw:
//   Br ∝ r^-2, Bφ = -Br * (Ω r sinθ / Vsw).  We choose Br(1AU) so that
//   |B|(1 AU) equals user-given B1AU.  Implementation caches Br1AU_T and
//   k_AU = Ω AU sinθ / Vsw, so evaluation is a few mults per sample.
//
// CME apex kinematics: Drag-Based Model (DBM): Vršnak & Žic (2007); Vršnak et al. (2013)
//   u(t) = Vsh − Vsw.  With drag Γ,
//     u(t) = u0 / (1 + Γ u0 t),   r(t) = r0 + Vsw t + (ln(1+Γ u0 t))/Γ.
//   We guard logs/denominators and convert Γ from km^-1 to m^-1.
//
// Regions and smoothing (self-similar):
//   upstream → [shock at R_sh] → sheath → [leading edge at R_LE] → magnetic ejecta
//   → [trailing edge at R_TE] → downstream ambient.
//   Each interface uses a C¹ smoothstep s(x)=x^2(3−2x) over width w; widths scale
//   ∝ R_sh (so they grow with distance). Sheath density uses a ramp that is steeper
//   near the shock (power p≥1) and bounded below by a floor (rc_floor≥1).
//
// Oblique-MHD shock proxy (1-D quasi-radial reduction): Edmiston & Kennel (1984);
// Priest (2014, CUP): We approximate the normal Mach number using the fast speed
//   c_f = sqrt(c_s^2 + v_A^2) with upstream sound speed c_s and Alfven speed v_A,
//   then take a hydrodynamic-like compression
//     rc = ((γ+1) M_n^2) / ((γ−1) M_n^2 + 2),   1 ≤ rc ≤ 4.
//   Downstream speed at the shock follows a continuity proxy.
//   This rc is used as: sheath density jump and a proxy for Bt amplification.
//
// Divergence of bulk flow:
//   ∇·V = (1/r^2) d/dr ( r^2 V_r ). We compute it with a centered FD using
//   r±dr along the same ray; dr = max(1e-4 AU, 1e-3 r).  Numerically robust.
//

// *** POST-SHOCK BEHAVIOR (PHYSICAL SANITY) ***
// --------------------------------------------
// • Right at the shock (immediate downstream), a *compressive fast shock* must have
//     – Density increase: rc = n2/n1 > 1 (≤ 4 for γ=5/3).

//     V2 ≈ V_sh + (V1 − V_sh)/rc,  so V2 lies typically between V_sw and V_sh.
//   You should therefore see a **density peak** and **elevated speed** right behind
//   the shock, followed by a **gradual decrease** through the **sheath**.
//
// • Through the **sheath** (shock → LE): it is physical for both **density** and
//   **speed** to **decline** from their immediate post-shock values as compression
//   relaxes and turbulence/expansion redistribute momentum.
//
// • Inside the **magnetic ejecta (ME)**: density is commonly **below ambient** and
//   speed can be **lower than upstream wind**—a well-known depletion region.
//
// • **Unphysical pattern to avoid**: an *immediate* (next-sample) **drop of density
//   below upstream** right behind the shock. That violates jump conditions for a
//   compressive fast shock. If you see this:
//     – Keep `sheath_comp_floor ≥ 1.0` (our default enforces ≥1).
//     – Use a small shock blend width vs. sheath thickness:
//         `edge_smooth_shock_AU_at1AU  <<  sheath_thick_AU_at1AU`.
//     – Avoid over-lapping large smooth widths at shock/LE/TE.
//     – Choose `V_sheath_LE_factor ≳ 1.05–1.2` so the sheath remains faster than ambient.
//
// • Optional (not enforced in this file): a **monotonicity clamp** within the sheath
//   to guarantee `n ≥ n_up` and `V_sw ≤ V ≤ V_sh`. If desired, we can provide a
//   compile-time or runtime switch; for now we document the physics and parameter
//   guidance above (no behavioral change).


 NUMERICAL NOTES
   • Sanitized outputs: n, V, B, |B|, ∇·V are finite (fallbacks applied).
   • No heap allocs in evaluators; no pow() in hot paths.
   • ∇·V uses 1/r^2 · d/dr(r^2 V) with a small centered difference (dr_frac).
   • Tunable widths per edge; keep shock width ≪ sheath thickness.

 UNITS
   • Inputs: V_sw [km/s], n1AU [cm^-3], B1AU [nT], T [K], r0 [R_sun], V0_sh [km/s], Γ [1/km].
   • Outputs: r [m], n [m^-3], V [m/s], Br/Bphi/|B| [T], divV [1/s].
   • Constants exposed: AU [m], Rs [m], PI, OMEGA_SUN [rad/s].

 EXAMPLE (quick start)
 ------------------------------------------------------------------------------
   #include "swcme1d.hpp"
   using namespace swcme1d;

   Model sw;  // default construct with sensible defaults

   sw.SetAmbient(400, 6, 5, 1.2e5)                  // V_sw[km/s], n1AU[cm^-3], B1AU[nT], T[K]
     .SetCME(1.05, 1800, 8e-8)                      // r0[R_sun], V0_sh[km/s], Γ[1/km]
     .SetGeometry(0.10, 0.20)                       // sheath & ME thickness at 1 AU [AU]
     .SetSmoothing(0.01, 0.02, 0.03)                // edge widths at 1 AU [AU] (shock/LE/TE)
     .SetSheathEjecta(1.2, 2.0, 1.10, 0.5, 0.8);    // rc_floor, ramp_p, Vshe_LE, fME, VME

   StepState S = sw.prepare_step(36.0*3600.0);      // 36 hours after launch

   const int N=3;
   double r[N]   ={0.5*AU, 1.0*AU, 1.5*AU};
   double n[N],V[N],Br[N],Bphi[N],Bmag[N],divV[N];

   sw.evaluate_radii_with_B_div(S, r, n, V, Br, Bphi, Bmag, divV, N);
   sw.write_tecplot_radial_profile(S, r, n, V, Br, Bphi, Bmag, divV, N, "profile_1d.dat");

 PARAMETER GUIDE (cheat sheet)
 ------------------------------------------------------------------------------
   Ambient:
     V_sw_kms   300–800      Upstream wind speed (radial)
     n1AU_cm3   2–10         Density at 1 AU (Leblanc-normalized)
     B1AU_nT    3–8          |B|(1 AU) for Parker normalization
     T_K        8e4–2e5      Proton temperature (for c_s)
     gamma_ad   5/3          Adiabatic index
     sin_theta  0.6–1.0      sin(colat); ≈1 in ecliptic

   CME apex (DBM):
     r0_Rs      1.03–1.07    Launch radius
     V0_sh_kms  800–2500     Initial shock speed
     Gamma_kmInv 1e-8–2e-7   Drag Γ (↑ ⇒ stronger decel)

   Geometry (at 1 AU; scales ∝ R_sh):
     sheath_thick_AU_at1AU   0.05–0.20
     ejecta_thick_AU_at1AU   0.10–0.40

   Edge widths (at 1 AU; scales ∝ R_sh):
     edge_smooth_shock_AU_at1AU  0.005–0.02   (keep smallest)
     edge_smooth_le_AU_at1AU     0.01–0.05
     edge_smooth_te_AU_at1AU     0.02–0.06

   Sheath/ME targets:
     sheath_comp_floor    ≥1.0 (1.1–1.5 typical)
     sheath_ramp_power    1–3  (2 steeper near shock)
     V_sheath_LE_factor   1.05–1.2
     f_ME                 0.3–0.8
     V_ME_factor          0.6–0.95

 REFERENCES
   Parker (1958) ApJ 128, 664 — Solar wind & spiral field.
   Leblanc et al. (1998) Sol. Phys. 183, 165 — n_e(R) ~ R^-2,-4,-6.
   Vršnak & Žic (2007) A&A 472, 937 — Drag-Based CME Model (DBM).
   Priest (2014) CUP — MHD of the Sun (shock relations, wave speeds).
===============================================================================
*/

#include <cmath>
#include <cstdio>
#include <cstddef>
#include <algorithm>

namespace swcme1d {

// -------------------- Physical constants (SI) --------------------
constexpr double PI         = 3.1415926535897932384626433832795;
constexpr double AU         = 1.495978707e11;     // Astronomical Unit [m]
constexpr double Rs         = 6.957e8;            // Solar radius [m]
constexpr double OMEGA_SUN  = 2.86533e-6;         // Solar rotation [rad/s] (~25.4 d sidereal)
constexpr double MU0        = 4.0e-7 * PI;        // Vacuum permeability [N/A^2]
constexpr double MP         = 1.67262192369e-27;  // Proton mass [kg]
constexpr double KB         = 1.380649e-23;       // Boltzmann [J/K]

// -------------------- Helpers (fast, branch-light) --------------------
inline double clamp01(double x){ return (x<0.0?0.0:(x>1.0?1.0:x)); }
inline double clamp(double x, double a, double b){ return (x<a?a:(x>b?b:x)); }

// C^1 cubic smoothstep on [0,1]; returns 0 at 0, 1 at 1, monotone
inline double smoothstep01(double x){ x=clamp01(x); return x*x*(3.0-2.0*x); }

// Asymmetric edge smoother: map signed distance d (m) with width w (m) into [0,1]
// where s=0 upstream side, s=1 downstream side (choose sign accordingly).
inline double edge_blend(double d, double w){
  if (w<=0.0) return (d>0.0)?1.0:0.0;
  return smoothstep01( 0.5 + 0.5 * clamp(d/w, -1.0, 1.0) );
}

// -------------------- User parameters --------------------
struct Params {
  // Ambient & thermodynamics
  double V_sw_kms    = 400.0;  // upstream wind speed [km/s]
  double n1AU_cm3    = 6.0;    // density at 1 AU [cm^-3]
  double B1AU_nT     = 5.0;    // |B|(1 AU) [nT]
  double T_K         = 1.2e5;  // proton temperature [K]
  double gamma_ad    = 5.0/3.0;
  double sin_theta   = 1.0;    // ≈1 ecliptic

  // CME launch & drag (DBM)
  double r0_Rs       = 1.05;   // launch radius [R_sun]
  double V0_sh_kms   = 1800.0; // initial shock speed [km/s]
  double Gamma_kmInv = 8e-8;   // drag parameter Γ [1/km]

  // Region geometry (normalized at 1 AU; scales ∝ R_sh)
  double sheath_thick_AU_at1AU = 0.10;
  double ejecta_thick_AU_at1AU = 0.20;

  // Edge smoothing (widths at 1 AU; scales ∝ R_sh)
  double edge_smooth_shock_AU_at1AU = 0.01;
  double edge_smooth_le_AU_at1AU    = 0.02;
  double edge_smooth_te_AU_at1AU    = 0.03;

  // Sheath/ME shaping
  double sheath_comp_floor   = 1.20; // ≥ 1
  double sheath_ramp_power   = 2.00; // 1..3
  double V_sheath_LE_factor  = 1.10; // × V_sw
  double f_ME                = 0.50; // × n_up
  double V_ME_factor         = 0.80; // × V_sw
};

// -------------------- Time cache returned by prepare_step(t) --------------------
struct StepState {
  // time
  double t_s = 0.0;

  // DBM shock kinematics
  double r_sh_m  = 0.0;   // shock radius [m]
  double V_sh_ms = 0.0;   // shock speed [m/s]

  // Region boundaries (outer→inner): R_sh > R_LE > R_TE
  double r_le_m  = 0.0;   // sheath → ME leading edge [m]
  double r_te_m  = 0.0;   // ME → ambient trailing edge [m]

  // Edge widths (scaled to current R_sh)
  double w_sh_m  = 0.0;   // shock smoothing width [m]
  double w_le_m  = 0.0;   // LE width [m]
  double w_te_m  = 0.0;   // TE width [m]

  // Thicknesses at current distance
  double dr_sheath_m = 0.0;
  double dr_me_m     = 0.0;

  // Upstream ambient at snapshot
  double V_up_ms = 0.0;   // = V_sw
  // Leblanc-like density coefficients in SI: n(r)=C2/r^2 + C4/r^4 + C6/r^6
  double C2=0.0, C4=0.0, C6=0.0;

  // Parker normalization for B
  double Br1AU_T = 0.0;   // Br at 1 AU [T]
  double k_AU    = 0.0;   // k = Ω AU sinθ / V_sw  [dimensionless]

  // Shock compression ratio (cap ≤4)
  double rc      = 1.0;

  // Convenience (ambient B at R_sh)
  double B_up_T  = 0.0;
};

// -------------------- Model class (header-only) --------------------
class Model {
public:
  // Default ctor (sensible defaults above)
  Model() : P{} {}

  // Construct from parameters
  explicit Model(const Params& p) : P(p) {}

  // Replace entire parameter pack
  inline Model& SetParams(const Params& p){ P = p; return *this; }

  // Set only the CME launch & drag (DBM)
  inline Model& SetCME(double r0_Rs, double V0_sh_kms, double Gamma_kmInv){
    P.r0_Rs = r0_Rs; P.V0_sh_kms = V0_sh_kms; P.Gamma_kmInv = Gamma_kmInv; return *this;
  }

  // Ambient & thermodynamics
  inline Model& SetAmbient(double V_sw_kms, double n1AU_cm3,
                           double B1AU_nT, double T_K,
                           double gamma_ad=5.0/3.0, double sin_theta=1.0)
  {
    P.V_sw_kms=V_sw_kms; P.n1AU_cm3=n1AU_cm3; P.B1AU_nT=B1AU_nT; P.T_K=T_K;
    P.gamma_ad=gamma_ad; P.sin_theta=sin_theta; return *this;
  }

  // Geometry
  inline Model& SetGeometry(double sheath_thick_AU_at1AU, double ejecta_thick_AU_at1AU){
    P.sheath_thick_AU_at1AU=sheath_thick_AU_at1AU; P.ejecta_thick_AU_at1AU=ejecta_thick_AU_at1AU; return *this;
  }

  // Edge widths
  inline Model& SetSmoothing(double w_sh_AU_at1AU, double w_le_AU_at1AU, double w_te_AU_at1AU){
    P.edge_smooth_shock_AU_at1AU=w_sh_AU_at1AU; P.edge_smooth_le_AU_at1AU=w_le_AU_at1AU; P.edge_smooth_te_AU_at1AU=w_te_AU_at1AU; return *this;
  }

  // Sheath/ME shaping
  inline Model& SetSheathEjecta(double sheath_comp_floor, double sheath_ramp_power,
                                double V_sheath_LE_factor, double f_ME, double V_ME_factor)
  {
    P.sheath_comp_floor=sheath_comp_floor; P.sheath_ramp_power=sheath_ramp_power;
    P.V_sheath_LE_factor=V_sheath_LE_factor; P.f_ME=f_ME; P.V_ME_factor=V_ME_factor; return *this;
  }

  inline const Params& GetParams() const { return P; }
  inline       Params& MutableParams()   { return P; }

  // Build per-time cache (DBM kinematics, Parker & Leblanc normals, rc, widths)
  StepState prepare_step(double t_s) const {
    StepState S; S.t_s = t_s;

    // ---- Ambient basics ----
    const double Vsw = P.V_sw_kms * 1.0e3; // [m/s]
    S.V_up_ms = Vsw;

    // Parker k, Br(1AU)
    S.k_AU = (Vsw>0.0) ? (OMEGA_SUN * AU * P.sin_theta / Vsw) : 0.0;
    const double B1AU_T  = P.B1AU_nT * 1.0e-9;
    S.Br1AU_T = (B1AU_T>0.0) ? (B1AU_T / std::sqrt(1.0 + S.k_AU*S.k_AU)) : 0.0;

    // Leblanc-like normalization to n(1AU)=n1AU_cm3
    const double n1AU = P.n1AU_cm3 * 1.0e6; // [m^-3]
    const double rAU  = AU;                 // [m]
    // Use Leblanc (1998) relative weights at 1 AU to partition the n1AU
    // (3.3e5 R^-2 + 4.1e6 R^-4 + 8.0e7 R^-6) with R in R_sun.
    const double R_AU_in_Rs = AU / Rs;
    const double w2 = 3.3e5 / (R_AU_in_Rs*R_AU_in_Rs);
    const double w4 = 4.1e6 / (R_AU_in_Rs*R_AU_in_Rs*R_AU_in_Rs*R_AU_in_Rs);
    const double w6 = 8.0e7 / (R_AU_in_Rs*R_AU_in_Rs*R_AU_in_Rs*R_AU_in_Rs*R_AU_in_Rs*R_AU_in_Rs);
    const double wsum = std::max(1e-300, w2+w4+w6);
    const double f2 = w2/wsum, f4 = w4/wsum, f6 = w6/wsum;
    // n(r)=C2/r^2 + C4/r^4 + C6/r^6 → enforce n(AU)=n1AU:
    S.C2 = f2 * n1AU * (rAU*rAU);
    S.C4 = f4 * n1AU * (rAU*rAU*rAU*rAU);
    S.C6 = f6 * n1AU * (rAU*rAU*rAU*rAU*rAU*rAU);

    // ---- DBM shock kinematics ----
    const double r0 = P.r0_Rs * Rs;              // [m]
    const double V0 = P.V0_sh_kms * 1.0e3;       // [m/s]
    const double Gamma = P.Gamma_kmInv * 1.0e-3; // [1/m]
    const double u0 = std::max(0.0, V0 - Vsw);
    const double denom = 1.0 + Gamma * u0 * std::max(0.0, t_s);
    const double u = (denom>0.0) ? (u0 / denom) : 0.0;
    S.V_sh_ms = Vsw + u;
    S.r_sh_m  = r0 + Vsw*t_s + ((Gamma>0.0) ? std::log(std::max(1.0, denom))/Gamma : u0*t_s);

    // ---- Region extents & widths (self-similar with R_sh) ----
    const double scale = S.r_sh_m / AU;  // ~1 at 1 AU
    S.dr_sheath_m = std::max(0.0, P.sheath_thick_AU_at1AU * scale * AU);
    S.dr_me_m     = std::max(0.0, P.ejecta_thick_AU_at1AU * scale * AU);
    S.r_le_m      = S.r_sh_m - S.dr_sheath_m;
    S.r_te_m      = S.r_le_m - S.dr_me_m;

    S.w_sh_m      = std::max(0.0, P.edge_smooth_shock_AU_at1AU * scale * AU);
    S.w_le_m      = std::max(0.0, P.edge_smooth_le_AU_at1AU    * scale * AU);
    S.w_te_m      = std::max(0.0, P.edge_smooth_te_AU_at1AU    * scale * AU);

    // ---- Upstream B (for rc) at shock radius ----
    const double r = std::max(S.r_sh_m, 1.05*Rs);
    const double r_AU = r / AU;
    const double Br   = (S.Br1AU_T) * (AU/r)*(AU/r);
    const double Bphi = -Br * (OMEGA_SUN * r * P.sin_theta / Vsw);
    S.B_up_T = std::sqrt(Br*Br + Bphi*Bphi);

    // ---- Compression ratio rc using fast-mode Mach proxy ----
    const double n_up = density_upstream(S, r);
    const double rho  = std::max(1e-30, MP * n_up);
    const double vA   = (S.B_up_T>0.0) ? (S.B_up_T / std::sqrt(MU0 * rho)) : 0.0;
    const double cs   = std::sqrt(std::max(0.0, P.gamma_ad * KB * P.T_K / MP));
    const double cf   = std::sqrt(cs*cs + vA*vA);
    const double Mf   = std::max(1.0, (S.V_sh_ms - Vsw) / std::max(1.0, cf));
    const double g    = P.gamma_ad;
    double rc = ((g+1.0)*Mf*Mf)/((g-1.0)*Mf*Mf + 2.0);
    rc = clamp(rc, 1.0, 4.0);
    // enforce floor (≥1)
    rc = std::max(rc, std::max(1.0, P.sheath_comp_floor));
    S.rc = rc;

    return S;
  }

  // Fast upstream density evaluator (SI), used internally & by adapters
  static inline double density_upstream(const StepState& S, double r_m){
    const double r = std::max(r_m, 1.05*Rs);
    // n(r) = C2/r^2 + C4/r^4 + C6/r^6
    const double inv2 = 1.0/(r*r);
    const double inv4 = inv2*inv2;
    const double inv6 = inv4*inv2;
    double n = S.C2*inv2 + S.C4*inv4 + S.C6*inv6;
    if (!std::isfinite(n) || n<0.0) n = 0.0;
    return n;
  }

  // -------------------- Evaluators --------------------

  // n(r), V(r) only
  void evaluate_radii_fast(const StepState& S,
                           const double* r_m, double* n_m3, double* V_ms,
                           std::size_t N) const
  {
    const double Vsw = S.V_up_ms;
    const double rc  = S.rc;

    for (std::size_t i=0;i<N;++i){
      const double r = std::max(r_m[i], 1.05*Rs);

      // Upstream ambient
      const double n_up = density_upstream(S, r);
      double n = n_up;
      double V = Vsw;

      // Region logic (outer→inner): shock | sheath | ME | inside TE
      if (r >= S.r_sh_m - S.w_sh_m) {
        // Upstream→just behind shock blend over w_sh
        const double V_dn = S.V_sh_ms + (Vsw - S.V_sh_ms)/rc; // immediate downstream
        const double d = S.r_sh_m - r;   // >0 inside post-shock
        const double s = edge_blend(d, S.w_sh_m);
        n = (1.0 - s)*n_up + s*(rc*n_up);
        V = (1.0 - s)*Vsw  + s*(V_dn);
      }
      else if (r > S.r_le_m + S.w_le_m) {
        // Sheath interior (exclude shock & LE skirts)
        const double xi = clamp( (S.r_sh_m - r) / std::max(1e-12, S.dr_sheath_m), 0.0, 1.0 );
        const double p  = clamp(P.sheath_ramp_power, 1.0, 5.0);
        const double w  = std::pow(1.0 - xi, p); // 1 at shock → 0 at LE
        const double n_target = n_up * (1.0 + (rc - 1.0)*w);
        const double V_le   = P.V_sheath_LE_factor * Vsw;
        const double V_dn   = S.V_sh_ms + (Vsw - S.V_sh_ms)/rc;
        const double V_target = V_le + (V_dn - V_le) * w;
        n = n_target; V = V_target;
      }
      else if (r >= S.r_le_m - S.w_le_m) {
        // Blend Sheath→ME across LE width
        const double s = edge_blend(S.r_le_m - r, S.w_le_m); // 0 at sheath side → 1 in ME
        // sheath-side targets at LE
        const double n_sh = n_up; // decayed to ~ambient at LE
        const double V_sh = P.V_sheath_LE_factor * Vsw;
        const double n_me = P.f_ME * n_up;
        const double V_me = P.V_ME_factor * Vsw;
        n = (1.0 - s)*n_sh + s*n_me;
        V = (1.0 - s)*V_sh + s*V_me;
      }
      else if (r > S.r_te_m + S.w_te_m) {
        // ME interior
        n = P.f_ME * n_up;
        V = P.V_ME_factor * Vsw;
      }
      else if (r >= S.r_te_m - S.w_te_m) {
        // Blend ME→ambient across TE
        const double s = edge_blend(r - S.r_te_m, S.w_te_m); // 0 in ME → 1 ambient
        const double n_me = P.f_ME * n_up;
        const double V_me = P.V_ME_factor * Vsw;
        n = (1.0 - s)*n_me + s*n_up;
        V = (1.0 - s)*V_me + s*Vsw;
      }
      else {
        // Inside trailing side (closer than TE - w_te) → ambient again
        n = n_up; V = Vsw;
      }

      // Sanitization
      if (!std::isfinite(n) || n<0.0) n = 0.0;
      if (!std::isfinite(V))          V = 0.0;

      n_m3[i] = n; V_ms[i] = V;
    }
  }

  // Full evaluator: n, V, Br, Bphi, |B|, divV
  void evaluate_radii_with_B_div(const StepState& S,
                                 const double* r_m,
                                 double* n_m3, double* V_ms,
                                 double* Br_T, double* Bphi_T, double* Bmag_T,
                                 double* divV, std::size_t N,
                                 double dr_frac=1e-3) const
  {
    const double Vsw  = S.V_up_ms;
    const double rc   = S.rc;
    const double kAU  = S.k_AU;
    const double Br1  = S.Br1AU_T;

    for (std::size_t i=0;i<N;++i){
      const double r = std::max(r_m[i], 1.05*Rs);
      const double r_AU = r / AU;

      // Upstream Parker field
      double Br  = (Br1) * (AU/r)*(AU/r);
      double Bph = -Br * (OMEGA_SUN * r * P.sin_theta / Vsw);
      double Bmag = std::sqrt(Br*Br + Bph*Bph);

      // Upstream ambient
      const double n_up = density_upstream(S, r);
      double n = n_up;
      double V = Vsw;

      // Region logic (mirrors fast evaluator), also adjust B in sheath
      if (r >= S.r_sh_m - S.w_sh_m) {
        // Shock skirt: blend upstream with immediate downstream targets
        const double V_dn = S.V_sh_ms + (Vsw - S.V_sh_ms)/rc;
        const double d = S.r_sh_m - r;
        const double s = edge_blend(d, S.w_sh_m);

        // Density, speed
        n = (1.0 - s)*n_up + s*(rc*n_up);
        V = (1.0 - s)*Vsw  + s*(V_dn);

        // B: keep Br continuous; amplify Bphi by factor fB(s) from 1 → rc
        const double fB = 1.0 + (rc - 1.0)*s;
        Bph *= fB;
        Bmag = std::sqrt(Br*Br + Bph*Bph);
      }
      else if (r > S.r_le_m + S.w_le_m) {
        // Sheath interior: decay from rc → 1 with xi
        const double xi = clamp( (S.r_sh_m - r) / std::max(1e-12, S.dr_sheath_m), 0.0, 1.0 );
        const double p  = clamp(P.sheath_ramp_power, 1.0, 5.0);
        const double w  = std::pow(1.0 - xi, p); // 1 at shock → 0 at LE

        const double n_target = n_up * (1.0 + (rc - 1.0)*w);
        const double V_le   = P.V_sheath_LE_factor * Vsw;
        const double V_dn   = S.V_sh_ms + (Vsw - S.V_sh_ms)/rc;
        const double V_target = V_le + (V_dn - V_le) * w;
        n = n_target; V = V_target;

        const double fB = 1.0 + (rc - 1.0)*w; // Bphi decay rc→1
        Bph *= fB; Bmag = std::sqrt(Br*Br + Bph*Bph);
      }
      else if (r >= S.r_le_m - S.w_le_m) {
        // LE blend sheath→ME
        const double s = edge_blend(S.r_le_m - r, S.w_le_m); // 0 sheath → 1 ME
        const double n_sh = n_up;
        const double V_sh = P.V_sheath_LE_factor * Vsw;
        const double n_me = P.f_ME * n_up;
        const double V_me = P.V_ME_factor * Vsw;
        n = (1.0 - s)*n_sh + s*n_me;
        V = (1.0 - s)*V_sh + s*V_me;

        // B: blend Bphi amplification (→1 in ME)
        const double fB_sheath = 1.0; // ~1 near LE on sheath side
        const double fB_me     = 1.0; // no extra amplification in ME
        const double fB        = (1.0 - s)*fB_sheath + s*fB_me;
        Bph *= fB; Bmag = std::sqrt(Br*Br + Bph*Bph);
      }
      else if (r > S.r_te_m + S.w_te_m) {
        // ME interior
        n = P.f_ME * n_up;
        V = P.V_ME_factor * Vsw;
        // B unchanged from Parker in this simple 1-D profile
      }
      else if (r >= S.r_te_m - S.w_te_m) {
        // TE blend ME→ambient
        const double s = edge_blend(r - S.r_te_m, S.w_te_m); // 0 ME → 1 ambient
        const double n_me = P.f_ME * n_up;
        const double V_me = P.V_ME_factor * Vsw;
        n = (1.0 - s)*n_me + s*n_up;
        V = (1.0 - s)*V_me + s*Vsw;
        // B Parker
      }
      else {
        // Inside trailing side → ambient again
        n = n_up; V = Vsw;
      }

      // Sanitization
      if (!std::isfinite(n) || n<0.0) n = 0.0;
      if (!std::isfinite(V))          V = 0.0;
      if (!std::isfinite(Br))         Br = 0.0;
      if (!std::isfinite(Bph))        Bph= 0.0;
      if (!std::isfinite(Bmag))       Bmag=std::sqrt(Br*Br + Bph*Bph);

      n_m3[i]=n; V_ms[i]=V; Br_T[i]=Br; Bphi_T[i]=Bph; Bmag_T[i]=Bmag;

      // Divergence: ∇·V = (1/r^2) d/dr (r^2 V)  (centered FD with small dr)
      if (divV){
        const double h = std::max( std::abs(dr_frac*r), 1.0e3 ); // ≥ 1 km
        const double rm = std::max(1.05*Rs, r - h);
        const double rp = r + h;

        double n_m, V_m, n_p, V_p;
        // reuse fast evaluator for V at rp, rm
        evaluate_radii_fast(S, &rm, &n_m, &V_m, 1);
        evaluate_radii_fast(S, &rp, &n_p, &V_p, 1);

        const double num = ( (rp*rp)*V_p - (rm*rm)*V_m );
        const double den = ( (rp - rm) * r * r );
        double d = (den!=0.0) ? (num / den) : 0.0;
        if (!std::isfinite(d)) d = 0.0;
        divV[i] = d;
      }
    }
  }

  // -------------------- Tecplot writer (POINT) --------------------
  // Writes: r[m], n[m^-3], V[m/s], Br[T], Bphi[T], Bmag[T], divV[s^-1], rc[-], R_sh[m], R_LE[m], R_TE[m]
// -------------------- Tecplot writer (POINT) --------------------
// Writes: x[m], R[AU], rSun[R_sun], n[m^-3], V[m/s], Br[T], Bphi[T], Bmag[T],
//         divV[s^-1], rc[-], R_sh[m], R_LE[m], R_TE[m]
bool write_tecplot_radial_profile(const StepState& S,
                                  const double* r_m,
                                  const double* n_m3, const double* V_ms,
                                  const double* Br_T, const double* Bphi_T, const double* Bmag_T,
                                  const double* divV,
                                  std::size_t N, const char* path,double time_simulation = -1.0) const
{
  if (!path || N==0) return false;
  std::FILE* fp = std::fopen(path, "w");
  if (!fp) return false;

  // TITLE with optional simulation time (Tecplot-safe)
  if (time_simulation >= 0.0) {
    std::fprintf(fp, "TITLE=\"1D SW+CME profile (Simulation time = %.3f s)\"\n", time_simulation);
  } else {
    std::fprintf(fp, "TITLE=\"1D SW+CME profile\"\n");
  }

  std::fprintf(fp,
    "VARIABLES=\"x[m]\",\"R[AU]\",\"rSun[R_sun]\",\"n[m^-3]\",\"V[m/s]\","
    "\"Br[T]\",\"Bphi[T]\",\"Bmag[T]\",\"divV[s^-1]\",\"rc\",\"R_sh[m]\",\"R_LE[m]\",\"R_TE[m]\"\n");
  std::fprintf(fp, "ZONE T=\"snapshot\", I=%zu, F=POINT\n", N);

  for (std::size_t i=0;i<N;++i){
    const double r    = safe_finite(r_m[i]);                     // meters
    const double R_AU = (std::isfinite(r) && AU>0.0) ? r/AU : 0.0;
    const double r_Rs = (std::isfinite(r) && Rs>0.0) ? r/Rs : 0.0;

    const double n    = safe_nonneg(n_m3  ? n_m3[i]   : 0.0);
    const double V    = safe_finite( V_ms ? V_ms[i]   : 0.0);
    const double Br   = safe_finite( Br_T ? Br_T[i]   : 0.0);
    const double Bph  = safe_finite(Bphi_T? Bphi_T[i] : 0.0);
    const double Bm   = safe_finite(Bmag_T? Bmag_T[i] : 0.0);
    const double dv   = safe_finite(divV  ? divV[i]   : 0.0);

    // x[m], R[AU], rSun[R_sun], n, V, Br, Bphi, Bmag, divV, rc, R_sh, R_LE, R_TE
    std::fprintf(fp, "%.9e %.9e %.9e %.9e %.9e %.9e %.9e %.9e %.9e %.6f %.9e %.9e %.9e\n",
                 r, R_AU, r_Rs, n, V, Br, Bph, Bm, dv, S.rc, S.r_sh_m, S.r_le_m, S.r_te_m);
  }

  std::fclose(fp);
  return true;
}

  // Convenience: compute and write from radii only
  bool write_tecplot_radial_profile_from_r(const StepState& S,
                                           const double* r_m, std::size_t N,
                                           const char* path,double time_simulation = -1.0) const
  {
    if (!r_m || N==0) return false;
    // stack buffers for small N; fallback to new[] if large (rare)
    const std::size_t M = N;
    double *n=nullptr,*V=nullptr,*Br=nullptr,*Bphi=nullptr,*Bmag=nullptr,*dv=nullptr;
    n    = new double[M]; V    = new double[M]; Br   = new double[M];
    Bphi = new double[M]; Bmag = new double[M]; dv   = new double[M];

    evaluate_radii_with_B_div(S, r_m, n, V, Br, Bphi, Bmag, dv, M);
    const bool ok = write_tecplot_radial_profile(S, r_m, n, V, Br, Bphi, Bmag, dv, M, path,time_simulation);

    delete [] n; delete [] V; delete [] Br; delete [] Bphi; delete [] Bmag; delete [] dv;
    return ok;
  }

private:
  Params P;

  static inline double safe_finite(double x){ return std::isfinite(x)? x:0.0; }
  static inline double safe_nonneg(double x){ return (std::isfinite(x)&&x>0.0)? x:0.0; }
};

} // namespace swcme1d

#endif // SWCME1D_HPP

