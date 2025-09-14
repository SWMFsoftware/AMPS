#ifndef SWCME1D_HPP
#define SWCME1D_HPP
/*
================================================================================
 swcme1d.hpp — Header-only 1-D Solar Wind + CME (DBM) model (clean re-write)
--------------------------------------------------------------------------------
PURPOSE
  Provide a lightweight, numerically efficient 1-D model of heliocentric
  solar wind plus a driven CME forward shock and its downstream structure
  (sheath and magnetic ejecta, ME). Designed to:
    • return n(r), V(r) and Parker B(r) with CME-induced modifications,
    • expose ∇·V = (1/r²) d/dr (r² V) for SEP adiabatic terms,
    • allow independent smoothing at three edges (shock, LE, TE),
    • be fast enough for per-particle queries in SEP solvers.

PHYSICS (compact)
  Upstream density (Leblanc, Dulk & Bougeret 1998, scaled to match n at 1 AU):
    n(r) = A (R☉/r)² + B (R☉/r)⁴ + C (R☉/r)⁶  [cm⁻³]  → convert to [m⁻³].
  Parker field (equatorial; uses sinθ):
    Br(r)  = Br(1 AU) (AU/r)²,
    Bφ(r)  = −Br(r) (Ω r sinθ / Vsw). Given |B|(1 AU), Br(1 AU) = B1AU/sqrt(1+k²),
    where k = Ω AU sinθ / Vsw.
  CME shock (Drag-Based Model, e.g., Vršnak & Žic 2007):
    u ≡ Vsh − Vsw, u(t) = u0/(1+Γ u0 t),  Rsh = r0 + Vsw t + ln(1+Γ u0 t)/Γ.
  Compression ratio proxy from fast-mode Mach number:
    cf = sqrt(cs² + vA²), Mf = max(1, (Vsh−Vsw)/cf),
    rc = ((γ+1)Mf²)/((γ−1)Mf²+2) clamped to [1,4] and to a user floor ≥1.
  Downstream structure:
    Regions along the ray (large r → small r): upstream → [Rsh] → SHEATH → [RLE]
    → ME → [RTE] → ambient. Thicknesses at 1 AU are user inputs and scale ∝ Rsh.
  **Monotone sheath construction (artifact-free):**
    Let s∈[0,1] parametrize distance from shock to leading edge: s=0 at Rsh,
    s=1 at RLE. We construct:
      n(s) = exp((1−s) ln(n2_sh) + s ln(n_up(RLE))), with n2_sh = rc·n_up(Rsh).
      V(s) = C¹ smoothstep from V2_sh to V_LE (V_LE ≥ Vsw).
    We clamp to never go below local upstream: n(r) ≥ n_up(r), V(r) ≥ Vsw.
  Tangential B amplification:
    Apply once, inside the sheath only, tapering from fB≈rc at the shock to 1 at
    the LE. No other multiplications of Bφ.

NUMERICAL
  Radii are clipped to ≥1.05 R☉ for stability. ∇·V via centered finite
  difference with h = max(1e−3 r, 1 km). Per-time quantities are cached.

USAGE (examples at bottom of file)
================================================================================
*/

#include <cmath>
#include <cstddef>
#include <cstdio>
#include <algorithm>

namespace swcme1d {

// --------------------------- Physical constants (SI) ---------------------------
constexpr double PI        = 3.1415926535897932384626433832795;
constexpr double AU        = 1.495978707e11;      // Astronomical Unit [m]
constexpr double Rs        = 6.957e8;             // Solar radius [m]
constexpr double OMEGA_SUN = 2.86533e-6;          // Solar rotation [rad/s]
constexpr double MU0       = 4.0e-7 * PI;         // Vacuum permeability [N/A²]
constexpr double MP        = 1.67262192369e-27;   // Proton mass [kg]
constexpr double KB        = 1.380649e-23;        // Boltzmann [J/K]

// --------------------------------- Helpers -----------------------------------
inline double clamp01(double x){ return x<0.0?0.0:(x>1.0?1.0:x); }
inline double clamp(double x,double a,double b){ return x<a?a:(x>b?b:x); }
inline double smoothstep01(double x){ x=clamp01(x); return x*x*(3.0-2.0*x); } // C¹
inline double lerp(double a,double b,double t){ return a + (b-a)*t; }

// ------------------------------ User parameters ------------------------------
struct Params {
  // Ambient & thermodynamics
  double V_sw_kms    = 400.0;  // upstream wind speed [km/s]
  double n1AU_cm3    = 6.0;    // density at 1 AU [cm⁻³]
  double B1AU_nT     = 5.0;    // |B|(1 AU) [nT]
  double T_K         = 1.2e5;  // proton temperature [K]
  double gamma_ad    = 5.0/3.0;
  double sin_theta   = 1.0;    // sin(colatitude) for Parker Bφ

  // CME launch & drag (DBM)
  double r0_Rs       = 1.05;   // launch radius [R☉]
  double V0_sh_kms   = 1500.0; // initial shock speed [km/s]
  double Gamma_kmInv = 8e-8;   // drag [1/km]

  // Geometry: thicknesses at 1 AU, scale ∝ R_sh
  double sheath_thick_AU_at1AU  = 0.10; // AU at 1 AU
  double ejecta_thick_AU_at1AU  = 0.25; // AU at 1 AU

  // Interface smoothing widths at 1 AU (C¹), scale ∝ R_sh
  double edge_smooth_shock_AU_at1AU = 0.01; // shock skirt (kept small)
  double edge_smooth_le_AU_at1AU    = 0.02; // sheath → ME
  double edge_smooth_te_AU_at1AU    = 0.03; // ME → ambient

  // Sheath / ME shaping
  double sheath_comp_floor   = 1.10; // min rc used to build sheath profile
  double sheath_ramp_power   = 2.0;  // controls steepness near shock (≥1)
  double V_sheath_LE_factor  = 1.10; // V at LE relative to V_sw (≥1)
  double f_ME                = 0.50; // ME density factor vs upstream (<1 typical)
  double V_ME_factor         = 0.80; // ME speed factor vs V_sw (<1 typical)
};

// ------------------------------ Per-time cache -------------------------------
struct StepState {
  // DBM / geometry
  double time_s    = 0.0;     // time since launch [s]
  double r0_m      = 1.05*Rs; // launch radius [m]
  double r_sh_m    = 30*AU;   // shock apex radius [m]
  double V_sh_ms   = 1.0e6;   // shock speed [m/s]
  double r_le_m    = 29*AU;   // leading edge [m]
  double r_te_m    = 28*AU;   // trailing edge [m]
  double w_sh_m    = 0.0;     // smoothing widths [m] (kept modest)
  double w_le_m    = 0.0;
  double w_te_m    = 0.0;

  // Ambient / Parker / Leblanc
  double V_up_ms   = 4.0e5;   // upstream wind speed [m/s]
  double Br1AU_T   = 0.0;     // Br at 1 AU [T]
  double k_AU      = 0.0;     // Ω AU sinθ / V_sw [-]
  double B_up_T    = 0.0;     // |B| at R_sh upstream [T]
  double C2        = 0.0;     // Leblanc scaled SI coefficients: n=C2/r²+C4/r⁴+C6/r⁶
  double C4        = 0.0;
  double C6        = 0.0;

  // Shock compression and cached boundary values
  double rc          = 1.0;     // compression ratio used for sheath profile
  double n_up_shock  = 0.0;     // upstream density at shock radius [m⁻3]
  double n_up_le     = 0.0;     // upstream density at leading edge [m⁻3]
  double V2_shock_ms = 0.0;     // immediate downstream speed at the shock [m/s]
  double V_LE_ms     = 0.0;     // sheath speed at the leading edge [m/s]
};

// --------------------------------- Model -------------------------------------
class Model {
public:
  Model() : P{} {}
  explicit Model(const Params& p) : P(p) {}

  // Parameter setters (fluent)
  Model& SetParams(const Params& p){ P=p; return *this; }
  Model& SetCME(double r0_Rs,double V0_sh_kms,double Gamma_kmInv){
    P.r0_Rs=r0_Rs; P.V0_sh_kms=V0_sh_kms; P.Gamma_kmInv=Gamma_kmInv; return *this; }
  Model& SetAmbient(double V_sw_kms,double n1AU_cm3,double B1AU_nT,double T_K,
                    double gamma_ad=5.0/3.0,double sin_theta=1.0){
    P.V_sw_kms=V_sw_kms; P.n1AU_cm3=n1AU_cm3; P.B1AU_nT=B1AU_nT; P.T_K=T_K;
    P.gamma_ad=gamma_ad; P.sin_theta=sin_theta; return *this; }
  Model& SetGeometry(double sheath_thick_AU_at1AU,double ejecta_thick_AU_at1AU){
    P.sheath_thick_AU_at1AU=sheath_thick_AU_at1AU;
    P.ejecta_thick_AU_at1AU=ejecta_thick_AU_at1AU; return *this; }
  Model& SetSmoothing(double w_sh,double w_le,double w_te){
    P.edge_smooth_shock_AU_at1AU=w_sh;
    P.edge_smooth_le_AU_at1AU   =w_le;
    P.edge_smooth_te_AU_at1AU   =w_te; return *this; }
  Model& SetSheathEjecta(double sheath_comp_floor,double sheath_ramp_power,
                         double V_sheath_LE_factor,double f_ME,double V_ME_factor){
    P.sheath_comp_floor=sheath_comp_floor; P.sheath_ramp_power=sheath_ramp_power;
    P.V_sheath_LE_factor=V_sheath_LE_factor; P.f_ME=f_ME; P.V_ME_factor=V_ME_factor;
    return *this; }

  const Params& GetParams() const { return P; }
        Params& MutableParams()   { return P; }

  // ---------------------------- Build per-time cache -------------------------
  StepState prepare_step(double t_s) const {
    StepState S; S.time_s=t_s; S.r0_m = std::max(1.05*Rs, P.r0_Rs*Rs);

    // Upstream wind
    const double Vsw = std::max(1.0, P.V_sw_kms*1.0e3); // [m/s]
    S.V_up_ms = Vsw;

    // Parker spiral constants from |B|(1 AU)
    const double B1AU_T = std::max(0.0, P.B1AU_nT)*1e-9;
    S.k_AU = (Vsw>0.0) ? (OMEGA_SUN*AU*P.sin_theta / Vsw) : 0.0;
    S.Br1AU_T = (B1AU_T>0.0) ? (B1AU_T / std::sqrt(1.0 + S.k_AU*S.k_AU)) : 0.0;

    // Leblanc coefficients scaled to match n(1 AU)
    const double A = 3.3e5, B = 4.1e6, C = 8.0e7; // [cm⁻³]
    const double n1AU_target = std::max(0.0, P.n1AU_cm3)*1e6; // to m⁻³
    const double rrat = Rs/AU;
    const double n1AU_base_cm3 = A*std::pow(rrat,2) + B*std::pow(rrat,4) + C*std::pow(rrat,6);
    const double n1AU_base = n1AU_base_cm3*1e6; // to m⁻³
    const double scale = (n1AU_base>0.0) ? (n1AU_target / n1AU_base) : 0.0;
    // Write n = C2/r² + C4/r⁴ + C6/r⁶ in SI with r in meters:
    S.C2 = scale * (A*1e6) * (Rs*Rs);
    S.C4 = scale * (B*1e6) * (Rs*Rs*Rs*Rs);
    S.C6 = scale * (C*1e6) * (Rs*Rs*Rs*Rs*Rs*Rs);

    // DBM apex kinematics
    const double r0 = S.r0_m;
    const double u0 = std::max(0.0, P.V0_sh_kms*1e3 - Vsw); // [m/s]
    const double Gamma = std::max(0.0, P.Gamma_kmInv/1e3);  // [m⁻¹]
    const double gtu = Gamma * u0 * t_s;
    const double denom = 1.0 + gtu;
    const double u    = (denom>0.0) ? (u0/denom) : 0.0;
    S.V_sh_ms = Vsw + u;
    S.r_sh_m  = r0 + Vsw*t_s + ((denom>0.0 && Gamma>0.0) ? std::log(denom)/Gamma : u0*t_s);
    S.r_sh_m  = std::max(S.r_sh_m, 1.1*Rs);

    // Geometry (self-similar thickness & blending widths)
    const double scale_R = S.r_sh_m / AU; // dimensionless
    const double d_sheath = std::max(0.0, P.sheath_thick_AU_at1AU) * scale_R * AU;
    const double d_me     = std::max(0.0, P.ejecta_thick_AU_at1AU) * scale_R * AU;
    S.r_le_m = std::max(1.05*Rs, S.r_sh_m - d_sheath);
    S.r_te_m = std::max(1.05*Rs, S.r_le_m - d_me);

    S.w_sh_m = std::max(0.0, P.edge_smooth_shock_AU_at1AU) * scale_R * AU;
    S.w_le_m = std::max(0.0, P.edge_smooth_le_AU_at1AU)    * scale_R * AU;
    S.w_te_m = std::max(0.0, P.edge_smooth_te_AU_at1AU)    * scale_R * AU;

    // Upstream |B| at shock for rc estimate
    const double Br_sh  = (S.Br1AU_T>0.0) ? (S.Br1AU_T * (AU/S.r_sh_m)*(AU/S.r_sh_m)) : 0.0;
    const double Bphi_sh= -Br_sh * (OMEGA_SUN * S.r_sh_m * P.sin_theta / Vsw);
    S.B_up_T = std::sqrt(Br_sh*Br_sh + Bphi_sh*Bphi_sh);

    // Upstream density and fast-mode speed at shock
    const double n_up_sh = density_upstream(S, S.r_sh_m);
    const double rho     = std::max(1e-30, MP * n_up_sh);
    const double vA      = (S.B_up_T>0.0) ? (S.B_up_T/std::sqrt(MU0*rho)) : 0.0;
    const double cs      = std::sqrt(std::max(0.0, P.gamma_ad*KB*P.T_K/MP));
    const double cf      = std::sqrt(cs*cs + vA*vA);

    // Compression ratio (proxy)
    const double Mf = std::max(1.0, (S.V_sh_ms - Vsw) / std::max(1.0, cf));
    double rc = ((P.gamma_ad+1.0)*Mf*Mf)/((P.gamma_ad-1.0)*Mf*Mf + 2.0);
    rc = clamp(rc, 1.0, 4.0);
    rc = std::max(rc, std::max(1.0, P.sheath_comp_floor));
    S.rc = rc;

    // Cache boundary values for a strictly monotone sheath
    S.n_up_shock  = n_up_sh;                        // upstream at shock
    S.n_up_le     = density_upstream(S, S.r_le_m);  // upstream at LE
    S.V2_shock_ms = S.V_sh_ms + (Vsw - S.V_sh_ms)/rc; // RH proxy in Sun frame
    S.V_LE_ms     = std::max(Vsw, P.V_sheath_LE_factor*Vsw);

    return S;
  }

  // Upstream Leblanc density (fast; SI). r is clipped ≥1.05 R☉ for stability
  static inline double density_upstream(const StepState& S, double r_m){
    const double r = std::max(r_m, 1.05*Rs);
    const double inv2 = 1.0/(r*r);
    const double inv4 = inv2*inv2;
    const double inv6 = inv4*inv2;
    double n = S.C2*inv2 + S.C4*inv4 + S.C6*inv6;
    return (std::isfinite(n) && n>0.0) ? n : 0.0;
  }

  // ----------------------------- Fast evaluator ------------------------------
  // Returns n[m⁻³] and V[m/s] at each radius (no B, no ∇·V). O(1) per query.
  void evaluate_radii_fast(const StepState& S,
                           const double* r_m,
                           double* n_m3, double* V_ms,
                           std::size_t N) const {
    if (!r_m || !n_m3 || !V_ms || N==0) return;

    const double Vsw = S.V_up_ms;
    const double rc  = S.rc;

    for (std::size_t i=0;i<N;++i){
      const double r = std::max(r_m[i], 1.05*Rs);
      const double n_up = density_upstream(S, r);
      double n = n_up;
      double V = Vsw;

      if (r > S.r_sh_m){
        // upstream ambient
        n = n_up; V = Vsw;

      } else if (r >= S.r_le_m){
        // ------------------------------ SHEATH ------------------------------
        const double Ls = std::max(1e-6, S.r_sh_m - S.r_le_m);
        const double s  = clamp01( (S.r_sh_m - r)/Ls ); // 0 at shock → 1 at LE

        // Density: tie n2(shock)=rc*n_up(R_sh) to n_up(R_LE) (log-linear)
        const double n2_sh = rc * S.n_up_shock;
        const double ln_n  = (1.0 - s)*std::log(std::max(1e-30, n2_sh))
                           + s        *std::log(std::max(1e-30, S.n_up_le));
        const double n_target = std::exp(ln_n);

        // Speed: C¹ smoothstep from V2(shock) to V_LE (≥ Vsw)
        const double t   = std::pow(s, std::max(1.0, P.sheath_ramp_power));
        const double sC1 = smoothstep01(t);
        const double V_target = lerp(S.V2_shock_ms, S.V_LE_ms, sC1);

        // Final (monotone safety vs local upstream)
        n = std::max(n_target, n_up);
        V = std::max(V_target, Vsw);

      } else if (r >= S.r_te_m){
        // ------------------------------- ME ---------------------------------
        const double n_me = std::max(n_up, P.f_ME * n_up);
        const double V_me = std::max(Vsw, P.V_ME_factor * Vsw);

        // Blend with sheath on the ME side of the LE (C¹ at LE)
        if (S.w_le_m>0.0 && r >= S.r_le_m - S.w_le_m){
          const double y  = clamp01( (S.r_le_m - r)/S.w_le_m ); // 0 at LE → 1 inward
          const double sC1 = smoothstep01(y);
          const double V_sheath_at_LE = S.V_LE_ms; // sheath boundary value
          const double n_sheath_at_LE = S.n_up_le; // sheath boundary value
          n = lerp(n_me, n_sheath_at_LE, sC1);
          V = lerp(V_me, V_sheath_at_LE, sC1);
        } else {
          n = n_me; V = V_me;
        }

      } else {
        // ------------------------------ AMBIENT (after TE) ------------------
        n = n_up; V = Vsw;
        // Blend to ME on the ambient side of the TE to make TE C¹
        if (S.w_te_m>0.0 && r >= S.r_te_m - S.w_te_m){
          const double z  = clamp01( (r - (S.r_te_m - S.w_te_m))/S.w_te_m ); // 0 far → 1 at TE
          const double sC1 = smoothstep01(z);
          const double n_me = std::max(1.0, P.f_ME) * n_up;
          const double V_me = std::max(1.0, P.V_ME_factor) * Vsw;
          n = lerp(n, n_me, sC1);
          V = lerp(V, V_me, sC1);
        }
      }

      if (!std::isfinite(n) || n<0.0) n = 0.0;
      if (!std::isfinite(V))          V = 0.0;
      n_m3[i] = n; V_ms[i] = V;
    }
  }

  // -------------------------- Full evaluator (with B, ∇·V) -------------------
  void evaluate_radii_with_B_div(const StepState& S,
                                 const double* r_m,
                                 double* n_m3, double* V_ms,
                                 double* Br_T, double* Bphi_T, double* Bmag_T,
                                 double* divV, std::size_t N,
                                 double dr_frac=1e-3) const {
    if (!r_m || N==0) return;

    // First get n & V
    evaluate_radii_fast(S, r_m, n_m3, V_ms, N);

    const double Vsw = S.V_up_ms;
    const double Br1 = S.Br1AU_T;

    for (std::size_t i=0;i<N;++i){
      const double r = std::max(r_m[i], 1.05*Rs);
      // Parker field upstream baseline
      double Br  = (Br1) * (AU/r)*(AU/r);
      double Bph = -Br * (OMEGA_SUN * r * P.sin_theta / Vsw);

      // Single amplification site: inside the sheath only (taper rc→1)
      if (r < S.r_sh_m && r >= S.r_le_m){
        const double Ls = std::max(1e-6, S.r_sh_m - S.r_le_m);
        const double s  = clamp01( (S.r_sh_m - r)/Ls ); // 0 at shock → 1 at LE
        const double fB = lerp(std::max(1.0, S.rc), 1.0, smoothstep01(s));
        Bph *= fB; // *** the only Bφ amplification ***
      }

      const double Bmag = std::sqrt(Br*Br + Bph*Bph);
      if (Br_T)   Br_T[i]   = std::isfinite(Br)?Br:0.0;
      if (Bphi_T) Bphi_T[i] = std::isfinite(Bph)?Bph:0.0;
      if (Bmag_T) Bmag_T[i] = std::isfinite(Bmag)?Bmag:0.0;

      // ∇·V via centered finite difference
      if (divV){
        const double h = std::max(1.0e3, std::abs(dr_frac*r)); // ≥1 km
        const double rm = std::max(1.05*Rs, r - h);
        const double rp = r + h;
        double n_m, V_m, n_p, V_p;
        evaluate_radii_fast(S, &rm, &n_m, &V_m, 1);
        evaluate_radii_fast(S, &rp, &n_p, &V_p, 1);
        const double num = (rp*rp*V_p - rm*rm*V_m);
        const double den = (rp - rm) * r * r;
        double d = (den!=0.0) ? (num/den) : 0.0;
        divV[i] = std::isfinite(d) ? d : 0.0;
      }
    }
  }

  // ------------------------------- Tecplot writer ----------------------------
  // Writes: r[m], R[AU], rSun[R_s], n[m^-3], V[m/s], Br[T], Bphi[T], Bmag[T],
  //         divV[s^-1], rc, R_sh[m], R_LE[m], R_TE[m]  (13 columns)
  bool write_tecplot_radial_profile(const StepState& S,
                                    const double* r_m,
                                    const double* n_m3,const double* V_ms,
                                    const double* Br_T,const double* Bphi_T,const double* Bmag_T,
                                    const double* divV,std::size_t N,
                                    const char* path,double time_simulation=-1.0) const {
    if (!r_m || !n_m3 || !V_ms || !Br_T || !Bphi_T || !Bmag_T || !path) return false;
    std::FILE* f = std::fopen(path, "w"); if (!f) return false;

    std::fprintf(f,"TITLE=\"1D SW+CME radial profile\"\n"); 
    std::fprintf(f,"VARIABLES=\"r[m]\",\"R[AU]\",\"rSun[R_s]\",\"n[m^-3]\",\"V[m/s]\",\"Br[T]\",\"Bphi[T]\",\"Bmag[T]\",\"divV[s^-1]\",\"rc\",\"R_sh[m]\",\"R_LE[m]\",\"R_TE[m]\"\n");  
    std::fprintf(f, "ZONE T=\"radial\", I=%zu, F=POINT\n", N);

    for (std::size_t i=0;i<N;++i){
      const double r = r_m[i];
      const double R_AU = r/AU;
      const double Rsun = r/Rs;
      const double dv = divV?divV[i]:0.0;
      std::fprintf(f,
        "% .9e % .9e % .9e % .9e % .9e % .9e % .9e % .9e % .9e % .9e % .9e % .9e % .9e\n",
        r, R_AU, Rsun,
        n_m3[i], V_ms[i], Br_T[i], Bphi_T[i], Bmag_T[i], dv,
        S.rc, S.r_sh_m, S.r_le_m, S.r_te_m);
    }

    if (time_simulation>=0.0){
      std::fprintf(f, "# t = %.3f s\n", time_simulation);
    }

    std::fclose(f);
    return true;
  }

  // Convenience wrapper from radii only
  bool write_tecplot_radial_profile_from_r(const StepState& S,
                                           const double* r_m,std::size_t N,
                                           const char* path,double time_simulation=-1.0) const {
    if (!r_m || N==0) return false;
    double *n=new double[N], *V=new double[N], *Br=new double[N], *Bph=new double[N], *Bm=new double[N], *dv=new double[N];
    evaluate_radii_with_B_div(S, r_m, n, V, Br, Bph, Bm, dv, N);
    const bool ok = write_tecplot_radial_profile(S, r_m, n, V, Br, Bph, Bm, dv, N, path, time_simulation);
    delete [] n; delete [] V; delete [] Br; delete [] Bph; delete [] Bm; delete [] dv; return ok;
  }

private:
  Params P;
};

// -----------------------------------------------------------------------------
// Usage examples (copy-paste)
// -----------------------------------------------------------------------------
/*
Example 1 — basic profile at a single time
-----------------------------------------
  using namespace swcme1d;
  Model sw;
  sw.SetAmbient(450, 5.5, 4.0, 1.5e5)
    .SetCME(1.05, 1800, 7.5e-8)
    .SetGeometry(0.12, 0.25)
    .SetSmoothing(0.008, 0.02, 0.03)
    .SetSheathEjecta(1.2, 2.0, 1.10, 0.5, 0.8);

  auto S = sw.prepare_step(24*3600.0);
  double rA[5] = {0.5*AU, 0.8*AU, 1.0*AU, 1.2*AU, 1.5*AU};
  double n[5], V[5], Br[5], Bp[5], Bm[5], dV[5];
  sw.evaluate_radii_with_B_div(S, rA, n, V, Br, Bp, Bm, dV, 5);

  // Optional: write Tecplot file
  sw.write_tecplot_radial_profile(S, rA, n, V, Br, Bp, Bm, dV, 5, "swcme_profile.dat", 24*3600.0);

Example 2 — per-particle usage (fast evaluator)
-----------------------------------------------
  auto S = sw.prepare_step(t_now);
  for (size_t p=0; p<NPART; ++p){
    double r = particle_r[p];
    double n_loc, V_loc;
    sw.evaluate_radii_fast(S, &r, &n_loc, &V_loc, 1);
    // use n_loc,V_loc for adiabatic losses, scattering rates, etc.
  }
*/

} // namespace swcme1d

#endif // SWCME1D_HPP

