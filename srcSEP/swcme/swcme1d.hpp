#pragma once
// ============================================================================
// swcme1d.hpp
// ----------------------------------------------------------------------------
// FAST 1-D SOLAR WIND + CME FORWARD-SHOCK MODEL  (API-SANITIZED, EFFICIENT)
// ----------------------------------------------------------------------------
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
// NUMERICAL / EFFICIENCY DESIGN
// -----------------------------
// • StepState caches: DBM state (R_sh,V_sh), Leblanc SI coefficients (C2,C4,C6),
//   Parker constants (Br1AU_T, k_AU), thermodynamics (c_s), and region geometry.
// • Hot loops: no std::pow; only multiplies/adds, 1/r^2 then r^-4,r^-6 by fusing.
// • Smoothsteps use pre-computed inv2w = 1/(2w) per interface.  No divisions in loops.
// • API is SANITIZED: every scalar written to user arrays is finite-checked and
//   replaced by a physically sensible fallback (see helpers).
// • Vectorization hint (GCC/Clang): IVDEP pragma; OpenMP toggle optional.
//
// UNITS
// -----
//   r [m], t [s];  n [m^-3];  V [m/s];  Br,Bφ,|B| [Tesla];  divV [s^-1].
//
// BUILD
// -----
//   g++ -std=c++17 -O3 -march=native demo1d.cpp -o demo1d
//   # (Optional OpenMP) add: -fopenmp
// ============================================================================

#include <cstddef>
#include <cmath>
#include <algorithm>
#include <cstdio>
#include <vector>

#ifndef SW1D_RESTRICT
# if defined(__clang__) || defined(__GNUC__)
#   define SW1D_RESTRICT __restrict__
# else
#   define SW1D_RESTRICT
# endif
#endif

#if defined(__GNUC__)
  #define SW1D_IVDEP _Pragma("GCC ivdep")
#else
  #define SW1D_IVDEP
#endif

namespace swcme1d {

// ------------------------------
// Physical constants (SI)
// ------------------------------
constexpr double PI   = 3.14159265358979323846;
constexpr double kB   = 1.380649e-23;
constexpr double mp   = 1.67262192369e-27;
constexpr double mu0  = 4.0e-7 * PI;
constexpr double Rs   = 6.957e8;               // solar radius [m]
constexpr double AU   = 1.495978707e11;        // astronomical unit [m]
constexpr double Omega_sun = 2.86533e-6;       // solar rotation [rad/s]

// ------------------------------
// Optional runtime I/O toggle
// ------------------------------
inline bool& _io_flag(){ static bool f=false; return f; }
inline void EnableIO(){ _io_flag()=true; }
inline void DisableIO(){ _io_flag()=false; }
inline bool IsIOEnabled(){ return _io_flag(); }

// ------------------------------
// Numeric safety helpers (API sanitation)
// ------------------------------
static inline double sw1d_finite_or(double v, double fallback) {
  return std::isfinite(v) ? v : fallback;
}
template <typename T>
static inline T sw1d_clamp_nonneg(T v) { return (v < T(0)) ? T(0) : v; }

// ------------------------------
// User parameters (physics + shaping)
// ------------------------------
struct Params {
  // Ambient @ 1 AU and thermodynamics
  double V_sw_kms    = 400.0;  // upstream solar wind speed [km/s]
  double n1AU_cm3    = 5.0;    // upstream density at 1 AU [cm^-3]
  double B1AU_nT     = 5.0;    // magnetic magnitude at 1 AU [nT]
  double T_K         = 1.2e5;  // proton temperature [K]
  double gamma_ad    = 5.0/3.0;

  // Parker pitch (≈1 in ecliptic); treat sinθ constant along the ray
  double sin_theta   = 1.0;

  // DBM initial conditions (launch at r0)
  double r0_Rs       = 1.05;   // [R_sun] from Sun center (0,0,0)
  double V0_sh_kms   = 1800.0; // initial shock speed [km/s]
  double Gamma_kmInv = 8.0e-8; // drag parameter Γ [1/km]

  // Region geometry (self-similar w.r.t. R_sh)
  double sheath_thick_AU_at1AU = 0.10; // sheath thickness at 1 AU [AU]
  double ejecta_thick_AU_at1AU = 0.20; // ejecta thickness at 1 AU [AU]

  // Independent C¹ smooth widths at (shock, LE, TE), scale ∝ R_sh
  double edge_smooth_shock_AU_at1AU = 0.01;
  double edge_smooth_le_AU_at1AU    = 0.02;
  double edge_smooth_te_AU_at1AU    = 0.03;

  // Sheath/ejecta target states
  double sheath_comp_floor     = 1.2;  // min compression at shock (≥1)
  double sheath_ramp_power     = 2.0;  // steeper near shock if >1
  double V_sheath_LE_factor    = 1.10; // sheath speed at LE = factor * V_sw
  double f_ME                  = 0.5;  // ejecta density = f_ME * n_up
  double V_ME_factor           = 0.8;  // ejecta speed = factor * V_sw
};

// ------------------------------
// Per-time cache (built by prepare_step; reused in hot loops)
// ------------------------------
struct StepState {
  // Time and shock state
  double t_s     = 0.0; // [s]
  double r_sh_m  = 0.0; // shock radius [m]
  double V_sh_ms = 0.0; // shock speed [m/s]
  double rc      = 1.0; // compression at shock (≥1)

  // Upstream / downstream speeds at the shock
  double V_up_ms = 0.0;
  double V_dn_ms = 0.0;

  // Region extents & widths
  double dr_sheath_m = 0.0, dr_me_m = 0.0;
  double r_le_m = 0.0, r_te_m = 0.0;
  double inv_dr_sheath = 0.0;

  // Smoothstep precomputed 1/(2w) to avoid divisions in loops
  double inv2w_sh = 0.0, inv2w_le = 0.0, inv2w_te = 0.0;

  // Targets for inner regions
  double rc_floor = 1.2;
  double V_sheath_LE_ms = 0.0;
  double V_ME_ms = 0.0;

  // Cached Leblanc coefficients in SI with Rs powers absorbed:
  // n(r) = C2 r^-2 + C4 r^-4 + C6 r^-6
  double C2=0.0, C4=0.0, C6=0.0;

  // Thermo/Parker
  double c_s = 0.0;               // sound speed [m/s]
  double Br1AU_T=0.0, k_AU=0.0;   // Parker normalization at 1 AU
};

// ------------------------------
// Model
// ------------------------------
class Model {
public:
  explicit Model(const Params& p): P(p) {}

  // Build caches for time t: DBM R_sh,V_sh; rc via fast Mach; regions & widths;
  // Parker/Leblanc constants.  All expensive algebra done once per tick.
  StepState prepare_step(double t_s) const {
    StepState S{};
    S.t_s=t_s;
    S.V_up_ms = V_sw();

    // --- DBM apex kinematics (analytic; guarded logs) -----------------------
    double r_sh, V_sh; shock_position_speed(t_s, r_sh, V_sh);
    S.r_sh_m  = r_sh; S.V_sh_ms = V_sh;

    // --- Leblanc SI coefficients with Rs powers absorbed (fast in loops) ----
    // Start with baseline cm^-3 profile at 1 AU to compute a scale.
    const double nLeb1AU = leblanc_cm3(AU); // baseline at 1 AU (cm^-3)
    const double scale_cm3 = (nLeb1AU>0.0)? (P.n1AU_cm3 / nLeb1AU) : 1.0;
    // Convert to SI and absorb Rs^k:  n(r) = C2/r^2 + C4/r^4 + C6/r^6
    const double A = 3.3e5*scale_cm3, B = 4.1e6*scale_cm3, C = 8.0e7*scale_cm3; // cm^-3
    const double Rs2=Rs*Rs, Rs4=Rs2*Rs2, Rs6=Rs4*Rs2;
    S.C2 = A*1e6*Rs2; S.C4 = B*1e6*Rs4; S.C6 = C*1e6*Rs6; // SI

    // --- Thermo wave speeds (sound) -----------------------------------------
    S.c_s = std::sqrt(std::max(0.0,P.gamma_ad) * kB * std::max(0.0,P.T_K) / mp);

    // --- Parker normalization at 1 AU ---------------------------------------
    S.k_AU = (Omega_sun * AU * P.sin_theta) / S.V_up_ms;          // dimensionless
    const double B1AU_T = P.B1AU_nT * 1e-9;
    S.Br1AU_T = B1AU_T / std::sqrt(1.0 + S.k_AU*S.k_AU);           // choose Br so |B|(1 AU)=B1AU

    // --- Shock strength (quasi-radial proxy) --------------------------------
    const double cf   = c_fast(S, S.r_sh_m);
    const double Vrel = std::max(0.0, S.V_sh_ms - S.V_up_ms);
    if (cf<=0.0 || Vrel<=0.0) {
      S.rc      = 1.0;
      S.V_dn_ms = S.V_up_ms;
    } else {
      const double Mn = Vrel / cf;
      const double g  = P.gamma_ad;
      double rc = ((g+1.0)*Mn*Mn)/((g-1.0)*Mn*Mn+2.0);
      S.rc = std::min(4.0, std::max(1.0, std::isfinite(rc)?rc:1.0));
      // 1-fluid continuity proxy at the discontinuity:
      S.V_dn_ms = S.V_sh_ms + (S.V_up_ms - S.V_sh_ms)/S.rc;
    }

    // --- Self-similar region geometry ---------------------------------------
    const double scaleR = S.r_sh_m / AU;
    const double dr_sheath = std::max(0.0,P.sheath_thick_AU_at1AU) * scaleR * AU;
    const double dr_me     = std::max(0.0,P.ejecta_thick_AU_at1AU) * scaleR * AU;
    S.dr_sheath_m = dr_sheath;
    S.dr_me_m     = dr_me;
    S.r_le_m      = S.r_sh_m - dr_sheath;
    S.r_te_m      = S.r_le_m - dr_me;
    S.inv_dr_sheath = (dr_sheath>0.0)? 1.0/dr_sheath : 0.0;

    // --- Independent C¹ smooth widths (cache inv(2w)) -----------------------
    const double w_sh = std::max(0.0,P.edge_smooth_shock_AU_at1AU) * scaleR * AU;
    const double w_le = std::max(0.0,P.edge_smooth_le_AU_at1AU)    * scaleR * AU;
    const double w_te = std::max(0.0,P.edge_smooth_te_AU_at1AU)    * scaleR * AU;
    S.inv2w_sh = (w_sh>0.0)? 0.5/w_sh : 0.0;
    S.inv2w_le = (w_le>0.0)? 0.5/w_le : 0.0;
    S.inv2w_te = (w_te>0.0)? 0.5/w_te : 0.0;

    // --- Targets in inner regions -------------------------------------------
    S.rc_floor       = std::max(1.0, std::min(P.sheath_comp_floor, S.rc));
    S.V_sheath_LE_ms = P.V_sheath_LE_factor * S.V_up_ms;
    S.V_ME_ms        = P.V_ME_factor        * S.V_up_ms;

    return S;
  }

  // =======================
  // FAST EVALUATORS (1-D r)
  // =======================

  // n(r), V(r) on an array of radii (API-sanitized outputs).
  void evaluate_radii_fast(const StepState& S,
                           const double* SW1D_RESTRICT r_m,
                           double*       SW1D_RESTRICT n_m3,
                           double*       SW1D_RESTRICT V_ms,
                           std::size_t N) const
  {
    const double Vup=S.V_up_ms, Vdn=S.V_dn_ms;
    const double rsh=S.r_sh_m, rle=S.r_le_m, rte=S.r_te_m;
    const double rc=S.rc, rc_floor=S.rc_floor, inv_drs=S.inv_dr_sheath;
    const double Vshe_LE=S.V_sheath_LE_ms, Vme=S.V_ME_ms;
    const double p=P.sheath_ramp_power;
    const double inv2w_sh=S.inv2w_sh, inv2w_le=S.inv2w_le, inv2w_te=S.inv2w_te;

    auto smooth01 = [](double x){ return (x<=0)?0.0: (x>=1)?1.0 : x*x*(3.0-2.0*x); };
    auto edge = [&](double r,double r0,double inv2w)->double{
      if (inv2w<=0.0) return (r<=r0)?1.0:0.0;
      const double t = 0.5 + (r0 - r)*inv2w;
      return smooth01(t);
    };
    auto ramp_pow = [&](double one_minus_xi)->double{
      // common fast path: p==1 → identity; p==2 → square; else pow
      if      (p==1.0) return one_minus_xi;
      else if (p==2.0) return one_minus_xi*one_minus_xi;
      else             return std::pow(one_minus_xi, p);
    };

    SW1D_IVDEP
    for (std::size_t i=0;i<N;++i){
      const double r  = std::max(r_m[i], 1.05*Rs);
      const double r2 = r*r;
      const double inv2 = 1.0/r2;          // r^-2
      const double inv4 = inv2*inv2;       // r^-4
      const double inv6 = inv4*inv2;       // r^-6

      // Upstream density from cached SI coefficients (no pow)
      const double n_up = S.C2*inv2 + S.C4*inv4 + S.C6*inv6; // [m^-3]

      // Sheath graded compression and speed
      const double xi_raw=(S.dr_sheath_m>0.0)? (rsh - r)*inv_drs : 0.0;
      const double xi = (xi_raw<=0)?0.0:((xi_raw>=1.0)?1.0:xi_raw);
      const double Csheath = rc_floor + (rc - rc_floor)*ramp_pow(1.0 - xi);
      const double n_sheath = Csheath * n_up;
      const double V_sheath = Vdn + (Vshe_LE - Vdn)*xi;

      const double n_ejecta = std::max(0.0, P.f_ME) * n_up;
      const double V_ejecta = Vme;

      double a = (rc>1.0)? edge(r, rsh, inv2w_sh):0.0;
      double n_mix = (1.0-a)*n_up + a*n_sheath;
      double V_mix = (1.0-a)*Vup  + a*V_sheath;

      a = (rc>1.0)? edge(r, rle, inv2w_le):0.0;
      n_mix = (1.0-a)*n_mix + a*n_ejecta;
      V_mix = (1.0-a)*V_mix + a*V_ejecta;

      a = (rc>1.0)? edge(r, rte, inv2w_te):0.0;
      n_mix = (1.0-a)*n_mix + a*n_up;
      V_mix = (1.0-a)*V_mix + a*Vup;

      // ---- API SANITIZATION (finite + nonnegative density) -----------------
      n_m3[i] = sw1d_finite_or(n_mix, n_up);
      V_ms[i] = sw1d_finite_or(V_mix, Vup);
      if (n_m3[i] < 0) n_m3[i] = 0;
    }
  }

  // n, V, Br, Bphi, |B|, divV on radii (single pass; API-sanitized outputs).
  void evaluate_radii_with_B_div(const StepState& S,
                                 const double* SW1D_RESTRICT r_m,
                                 double*       SW1D_RESTRICT n_m3,
                                 double*       SW1D_RESTRICT V_ms,
                                 double*       SW1D_RESTRICT Br_T,
                                 double*       SW1D_RESTRICT Bphi_T,
                                 double*       SW1D_RESTRICT Bmag_T,
                                 double*       SW1D_RESTRICT divV,
                                 std::size_t N,
                                 double dr_frac=1e-3) const
  {
    evaluate_radii_fast(S,r_m,n_m3,V_ms,N);

    // Parker spiral components (few multiplies; sanitize)
    SW1D_IVDEP
    for (std::size_t i=0;i<N;++i){
      const double r    = std::max(r_m[i], 1.05*Rs);
      const double r_AU = r / AU;
      const double Br   = S.Br1AU_T / (r_AU*r_AU);
      const double Bphi = - Br * S.k_AU * r_AU;

      Br_T[i]   = sw1d_finite_or(Br,   0.0);
      Bphi_T[i] = sw1d_finite_or(Bphi, 0.0);

      const double Bm2 = Br_T[i]*Br_T[i] + Bphi_T[i]*Bphi_T[i];
      Bmag_T[i] = sw1d_finite_or(std::sqrt(std::max(0.0, Bm2)), 0.0);
    }

    // Divergence via centered FD of (1/r^2) d(r^2 V)/dr  (guarded radii/steps)
    const double dr_min = 1.0e-4*AU;
    SW1D_IVDEP
    for (std::size_t i=0;i<N;++i){
      const double r  = std::max(r_m[i], 1.05*Rs);
      const double dr = std::max(dr_min, dr_frac*r);
      const double rm = std::max(1.05*Rs, r - dr);
      const double rp = r + dr;

      double Vp, Vm;
      sample_V_only(S, rp, Vp);
      sample_V_only(S, rm, Vm);

      const double num = (rp*rp*Vp - rm*rm*Vm) / std::max(1e-12, rp-rm);
      const double dv  = num / (r*r);

      divV[i] = sw1d_finite_or(dv, 0.0);
    }
  }

  // =======================
  // Writers (NaN-safe)
  // =======================

  // Tecplot 1-D radial profile (POINT). Variables (columns):
  //   R[m], n[m^-3], V[m/s], Br[T], Bphi[T], Bmag[T], divV[1/s],
  //   rc[-], R_sh[m], R_LE[m], R_TE[m]
  bool write_tecplot_radial_profile(const StepState& S,
                                    const double* r_m,
                                    const double* n_m3,
                                    const double* V_ms,
                                    const double* Br_T,
                                    const double* Bphi_T,
                                    const double* Bmag_T,
                                    const double* divV,
                                    std::size_t N,
                                    const char* path) const
  {
    std::FILE* fp=std::fopen(path,"w"); if(!fp) return false;
    auto F=[](double v){ return std::isfinite(v)? v:0.0; };

    std::fprintf(fp,"TITLE=\"SW+CME 1D radial profile\"\n");
    std::fprintf(fp,"VARIABLES=\"R[m]\",\"n[m^-3]\",\"V[m/s]\",\"Br[T]\",\"Bphi[T]\",\"Bmag[T]\",\"divV[1/s]\","
                    "\"rc[-]\",\"R_sh[m]\",\"R_LE[m]\",\"R_TE[m]\"\n");
    std::fprintf(fp,"ZONE T=\"profile\", N=%zu, DATAPACKING=POINT\n", N);
    for (std::size_t i=0;i<N;++i){
      std::fprintf(fp,"%.9e %.9e %.9e %.9e %.9e %.9e %.9e %.6f %.9e %.9e %.9e\n",
        F(r_m[i]), F(n_m3[i]), F(V_ms[i]), F(Br_T[i]), F(Bphi_T[i]), F(Bmag_T[i]), F(divV[i]),
        F(S.rc), F(S.r_sh_m), F(S.r_le_m), F(S.r_te_m));
    }
    std::fclose(fp);
    return true;
  }

  // Convenience: compute arrays internally and write Tecplot line.
  bool write_tecplot_radial_profile_from_r(const StepState& S,
                                           const double* r_m, std::size_t N,
                                           const char* path) const
  {
    std::vector<double> n(N), V(N), Br(N), Bp(N), Bm(N), dv(N);
    evaluate_radii_with_B_div(S, r_m, n.data(), V.data(), Br.data(), Bp.data(), Bm.data(), dv.data(), N);
    return write_tecplot_radial_profile(S, r_m, n.data(), V.data(), Br.data(), Bp.data(), Bm.data(), dv.data(), N, path);
  }

  // =======================
  // Optional CSV logging
  // =======================
  static inline void write_step_csv_header(std::FILE* fp=stdout){
    if(!IsIOEnabled()) return; if(!fp) fp=stdout;
    std::fprintf(fp,"t_s,R_sh_AU,V_sh_km_s,rc,V_up_km_s,V_dn_km_s,R_LE_AU,R_TE_AU\n");
  }
  static inline void write_step_csv_line(const StepState& S, std::FILE* fp=stdout){
    if(!IsIOEnabled()) return; if(!fp) fp=stdout;
    std::fprintf(fp,"%.3f,%.9f,%.6f,%.6f,%.6f,%.6f,%.9f,%.9f\n",
      S.t_s, S.r_sh_m/AU, S.V_sh_ms/1e3, S.rc, S.V_up_ms/1e3, S.V_dn_ms/1e3, S.r_le_m/AU, S.r_te_m/AU);
  }

  // Upstream speed [m/s]
  inline double V_sw() const { return P.V_sw_kms*1e3; }

private:
  Params P;

  // ---- micro-helpers (private) ---------------------------------------------

  // Leblanc baseline (cm^-3) before scaling (used once in prepare_step)
  static inline double leblanc_cm3(double r_m){
    const double s=Rs/r_m; const double s2=s*s, s4=s2*s2, s6=s4*s2;
    return 3.3e5*s2 + 4.1e6*s4 + 8.0e7*s6;
  }

  // Parker magnitude at r from S (uses cached Br1AU_T, k_AU) — not used in loops,
  // but kept for completeness.
  static inline double B_parker_mag(const StepState& S, double r_m){
    const double r_AU = r_m / AU;
    const double Br   = S.Br1AU_T / (r_AU*r_AU);
    const double Bphi = - Br * S.k_AU * r_AU;
    return std::sqrt(Br*Br + Bphi*Bphi);
  }

  // Alfven speed using cached Leblanc SI coefficients (fast)
  static inline double v_A(const StepState& S, double r_m){
    const double r2=r_m*r_m, inv2=1.0/r2, inv4=inv2*inv2, inv6=inv4*inv2;
    const double n_up = S.C2*inv2 + S.C4*inv4 + S.C6*inv6; // m^-3
    const double rho  = mp * std::max(1e-3, n_up);
    const double B    = B_parker_mag(S, r_m);
    return B/std::sqrt(mu0*rho);
  }

  // Fast magnetosonic speed ≈ sqrt(c_s^2 + v_A^2)
  static inline double c_fast(const StepState& S, double r_m){
    const double vA=v_A(S,r_m); return std::sqrt(S.c_s*S.c_s + vA*vA);
  }

  // DBM apex r(t), V(t) (analytic; guards for numerics)
  inline void shock_position_speed(double t, double& r_sh_m, double& V_sh_ms) const {
    const double r0  = P.r0_Rs*Rs;
    const double Vsw = V_sw();
    const double V0  = P.V0_sh_kms*1e3;
    const double G   = P.Gamma_kmInv*1e-3; // [1/m]
    const double u0  = V0 - Vsw;

    if (std::fabs(u0) < 1e-12){ r_sh_m = r0 + Vsw*t; V_sh_ms = Vsw; return; }
    if (u0 > 0.0){
      const double den = 1.0 + G*u0*t;
      const double u   = u0 / den;
      const double I   = (1.0/G) * std::log(std::max(den,1e-30));
      r_sh_m = r0 + Vsw*t + I;
      V_sh_ms= Vsw + u;
    } else {
      double den = 1.0 - G*u0*t; if (den<=0.0) den=1e-30;
      const double u   = u0 / den;
      const double I   = -(1.0/G) * std::log(den);
      r_sh_m = r0 + Vsw*t + I;
      V_sh_ms= Vsw + u;
    }
  }

  // Sample only V(r) using exactly the same blending as evaluate_radii_fast
  inline void sample_V_only(const StepState& S, double r, double& V_out) const {
    r = std::max(r, 1.05*Rs);
    const double Vup=S.V_up_ms, Vdn=S.V_dn_ms;
    const double rsh=S.r_sh_m, rle=S.r_le_m, rte=S.r_te_m;
    const double rc=S.rc, inv_drs=S.inv_dr_sheath;
    const double Vshe_LE=S.V_sheath_LE_ms, Vme=S.V_ME_ms;
    const double p=P.sheath_ramp_power;
    const double inv2w_sh=S.inv2w_sh, inv2w_le=S.inv2w_le, inv2w_te=S.inv2w_te;

    auto smooth01 = [](double x){ return (x<=0)?0.0: (x>=1)?1.0 : x*x*(3.0-2.0*x); };
    auto edge = [&](double r,double r0,double inv2w)->double{
      if (inv2w<=0.0) return (r<=r0)?1.0:0.0;
      const double t = 0.5 + (r0 - r)*inv2w;
      return smooth01(t);
    };
    auto ramp_pow = [&](double one_minus_xi)->double{
      if      (p==1.0) return one_minus_xi;
      else if (p==2.0) return one_minus_xi*one_minus_xi;
      else             return std::pow(one_minus_xi, p);
    };

    const double xi_raw=(S.dr_sheath_m>0.0)? (rsh - r)*inv_drs : 0.0;
    const double xi = (xi_raw<=0)?0.0:((xi_raw>=1.0)?1.0:xi_raw);
    const double Csheath = S.rc_floor + (rc - S.rc_floor)*ramp_pow(1.0 - xi);
    (void)Csheath; // Csheath not needed here; kept to mirror density ramp logic.
    const double V_sheath = Vdn + (Vshe_LE - Vdn)*xi;

    double a = (rc>1.0)? edge(r, rsh, inv2w_sh):0.0;
    double V_mix = (1.0-a)*Vup  + a*V_sheath;

    a = (rc>1.0)? edge(r, rle, inv2w_le):0.0;
    V_mix = (1.0-a)*V_mix + a*Vme;

    a = (rc>1.0)? edge(r, rte, inv2w_te):0.0;
    V_mix = (1.0-a)*V_mix + a*Vup;

    V_out = sw1d_finite_or(V_mix, Vup);
  }
};

} // namespace swcme1d
