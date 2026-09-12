#ifndef SWCME1D_HPP
#define SWCME1D_HPP
/*
================================================================================
 swcme1d.hpp — Header-only 1-D Solar Wind + CME (DBM) model (artifact-free)
--------------------------------------------------------------------------------
OVERVIEW
  This header implements a fast, numerically robust, *1‑D along a radial ray*
  model of the ambient solar wind plus an outward-propagating CME forward shock
  and its downstream structure (compressed sheath and magnetic ejecta, “ME”).
  It is designed for per-particle queries inside SEP transport codes: given a
  radius r and time t since launch, it returns number density n(r), bulk speed
  V(r), Parker spiral B(r)=(Br,Bφ), and the divergence ∇·V(r). It also exposes
  the instantaneous geometry: shock radius R_sh, sheath→ME leading edge R_LE,
  and ME trailing edge R_TE.

  The implementation is *header-only* and uses no dynamic allocations in the
  hot path. All heavy per-time quantities are cached in a StepState, so the
  evaluators are O(1) in r.

PHYSICAL MODEL (succinct but complete)
  • Ambient density n_up(r): Leblanc, Dulk & Bougeret (1998)
      n_cm³(r) = A (R☉/r)^2 + B (R☉/r)^4 + C (R☉/r)^6,
      with canonical A=3.3×10⁵, B=4.1×10⁶, C=8.0×10⁷ [cm⁻³]. We *scale* these
      to match a user-specified n(1 AU) and convert to SI [m⁻³].

  • Magnetic field: equatorial Parker spiral
      Br(r)  = Br(1 AU) (AU/r)^2,
      Bφ(r)  = −Br(r) (Ω r sinθ / V_sw),
      |B|(1 AU) is provided by the user; we infer Br(1 AU)=B1AU/sqrt(1+k²),
      where k≡Ω AU sinθ / V_sw.

  • CME apex kinematics: shared BALLISTIC / DBM / DATA_DRIVEN engine.
      For DBM, ΔV0=V0−Vsw and a=|ΔV0|:
      ΔV(t)=ΔV0/(1+Γ a t),
      R_sh(t)=r0+Vsw t+sgn(ΔV0) log(1+Γ a t)/Γ.
      This sign-aware form decelerates fast CMEs and accelerates slow CMEs
      toward Vsw.  Γ=0 is handled by the exact ballistic limit.

  • Compression and downstream state come from the shared ideal-MHD fast-shock
      solver. A geometric front is a physical shock only when the normal
      shock-frame inflow is super-fast; otherwise compression is exactly one.

  • Downstream structure is controlled by swcme_regions.hpp. SHOCK_ONLY leaves
      the Parker/Leblanc transport background untouched. FULL_ICME uses
      upstream → [R_sh] → sheath → [R_LE] → magnetic ejecta → [R_TE] → ambient,
      with self-similar local thickness fractions.  In RESOLVED_COMPRESSION mode
      the mathematical RH jump is represented in the transport field by one
      finite C1 shock layer; LE/TE transitions are C1 as well.

  • *Crucial sheath construction (artifact-free):*
      Let s∈[0,1] map R_sh→R_LE with s=0 at the shock. We **pin the boundary
      values** and build:
      – density:    n(s) = exp((1−s) ln n₂ + s ln n_up(R_LE)), with n₂=r_c n_up(R_sh),
      – velocity:   V(s) = smoothstep(s^p; V₂→V_LE),  p≥1,
        where V₂ = V_sh − (V_sh−V_sw)/r_c is the downstream speed from mass-flux
        continuity (RH proxy) and V_LE ≥ V_sw (user-set factor).
      In RESOLVED_COMPRESSION mode the inner edge of the numerical shock layer
      is the exact RH state. The sheath starts there and relaxes toward the
      ambient leading-edge target without an empirical compression floor.

  • Magnetic field in FULL_ICME starts from the exact RH downstream vector and
      relaxes to the Parker field at R_LE. The ejecta field remains Parker in
      this intentionally simple baseline model.

NUMERICAL / IMPLEMENTATION CHOICES
  • r < 1.05 R☉ is outside the analytical model domain and is reported
    explicitly; radii are never silently clipped into the supported domain.
  • Per-time quantities (DBM kinematics, Parker parameters, Leblanc scale,
    r_c, sheath/ME geometric radii & widths) are cached in StepState.
  • Because the 1-D velocity is purely radial, ∇·V is evaluated analytically
    as 2 V_r/r + dV_r/dr.  The derivative comes from the same common C1 region
    profile that returns V(r), so no finite-difference noise enters SEP adiabatic
    energy change.
  • In RESOLVED_COMPRESSION, shock/LE/TE blending windows are C¹ and scale
    ∝ R_sh.  The shock layer blends upstream to the exact RH state, reaching
    that state at its inner edge; the subsequent sheath starts from the same
    endpoint with zero derivative.  SOURCE mode has no shock blend at all.

API & UNITS
  • All inputs/outputs documented per method below; defaults live in Params.
  • Units: radii [m], speeds [m/s], densities [m⁻³], magnetic field [T].
  • Public entry points:
      - Params (user configuration), StepState (per-time cache),
      - Model::prepare_step(t)
      - Model::evaluate_radii_fast(S, r[], n[], V[], N)
      - Model::evaluate_radii_with_B_div(S, r[], n[], V[], Br[], Bφ[], |B|[], divV[], N)
      - Model::write_tecplot_radial_profile(...)

COMMON QUESTIONS
  Q: “Why does the speed drop to ~320 km/s somewhere?”
     A: That’s the **ME** default (V_ME_factor=0.80) applied to V_sw=400 → 320.
        In the **sheath**, V≥V_sw always. Adjust with SetSheathEjecta(..., V_ME_factor).

  Q: “Is it physical for downstream to be below upstream right behind a forward
     shock?”  A: No. In the Sun frame V₂ > V_sw and n₂ > n_up. If you see otherwise,
     it’s a blending/branching bug. This implementation makes that impossible.

VALIDATION IDEAS (quick)
  • Check invariants at the shock: V₂>V_sw, n₂=r_c n_up.
  • Verify sheath monotonicity: n decreases from n₂ to n_up(R_LE); V relaxes
    from V₂ to V_LE≥V_sw.
  • Compare ambient n(r) to Leblanc curve; Parker |B|(r)∝r⁻² far out.

REFERENCES
  – Leblanc, Dulk, Bougeret (1998), Solar Phys., 183, 165–180 — density model.
  – Parker, E. N. (1958), ApJ, 128, 664 — spiral field.
  – Vršnak, B., & Žic, T. (2007), A&A, 472, 937 — drag-based CME model (DBM).
  – Vršnak et al. (2013), Sol. Phys., 285, 295 — DBM extensions & applications.
  – Priest, E. (2014), Magnetohydrodynamics of the Sun — shock & MHD basics.

USAGE SKETCH (more complete examples at bottom)
  using namespace swcme1d;
  Model m;
  m.SetAmbient(400, 6, 5, 1.2e5)            // km/s, cm^-3, nT, K
   .SetCME(20.0, 1500, 8e-8)                // R☉, km/s, 1/km
   .SetGeometry(0.10, 0.25)                 // Δsheath@1AU, ΔME@1AU (AU)
   .SetSmoothing(0.01, 0.02, 0.03)          // shock/LE/TE widths @1AU (AU)
   .SetSheathEjecta(1.15, 2.0, 1.10, 0.5, 1.0);  // keep ME speed ≥ V_sw
  auto S = m.prepare_step(36*3600.0);
  double r[3] = {0.6*AU, 1.0*AU, 1.4*AU};
  double n[3], V[3];
  m.evaluate_radii_fast(S, r, n, V, 3);
================================================================================
*/

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
   • CME apex kinematics use the shared swcme::kinematics engine.  DBM uses
        ΔV(t)=ΔV0/(1+Γ|ΔV0|t)
        R_sh=r0+Vsw t+sgn(ΔV0) log(1+Γ|ΔV0|t)/Γ,
     with an exact Γ=0 ballistic branch.  DATA_DRIVEN mode uses monotone PCHIP
     height-time interpolation and returns its derivative as the apex speed.
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
//     ΔV(t)=ΔV0/(1+Γ|ΔV0|t), r(t)=r0+Vsw t+sgn(ΔV0)ln(1+Γ|ΔV0|t)/Γ.
//   We guard logs/denominators and convert Γ from km^-1 to m^-1.
//
// Regions and smoothing (self-similar):
//   SHOCK_ONLY returns the analytical upstream background everywhere. FULL_ICME
//   uses upstream → [shock] → sheath → [R_LE] → magnetic ejecta → [R_TE] →
//   post-ICME ambient. R_LE/R_sh and R_TE/R_sh are fixed fractions.  With
//   RESOLVED_COMPRESSION the transport-facing RH jump is represented by one C1
//   shock layer; SOURCE mode uses SHOCK_ONLY and leaves the flow uncompressed.
//   The artificial LE/TE interfaces are C1-smoothed in FULL_ICME as well.
//
// Shock physics is provided by the shared ideal-MHD Rankine-Hugoniot solver in
// swcme_shock.hpp; the region model consumes that validated downstream state and
// does not maintain its own compression formula or floor.
//
// Divergence of bulk flow:
//   ∇·V = 2 V_r/r + dV_r/dr.  The common region profile supplies an
//   analytical dV_r/dr, so the 1-D result contains no finite-difference noise.
//   (Legacy text below describing a centered FD has been superseded.)
////   r±dr along the same ray; dr = max(1e-4 AU, 1e-3 r).  Numerically robust.
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
//     ΔV(t)=ΔV0/(1+Γ|ΔV0|t), r(t)=r0+Vsw t+sgn(ΔV0)ln(1+Γ|ΔV0|t)/Γ.
//   We guard logs/denominators and convert Γ from km^-1 to m^-1.
//
// Regions and smoothing (self-similar):
//   SHOCK_ONLY returns the analytical upstream background everywhere. FULL_ICME
//   uses upstream → [shock] → sheath → [R_LE] → magnetic ejecta → [R_TE] →
//   post-ICME ambient. R_LE/R_sh and R_TE/R_sh are fixed fractions.  With
//   RESOLVED_COMPRESSION the transport-facing RH jump is represented by one C1
//   shock layer; SOURCE mode uses SHOCK_ONLY and leaves the flow uncompressed.
//   The artificial LE/TE interfaces are C1-smoothed in FULL_ICME as well.
//
// Shock physics is provided by the shared ideal-MHD Rankine-Hugoniot solver in
// swcme_shock.hpp; the region model consumes that validated downstream state and
// does not maintain its own compression formula or floor.
//
// Divergence of bulk flow:
//   ∇·V = 2 V_r/r + dV_r/dr with an analytical derivative of the exact
//   production V(r) blend; no radial finite-difference stencil is used.
////   r±dr along the same ray; dr = max(1e-4 AU, 1e-3 r).  Numerically robust.
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
// • The canonical shock diagnostic remains the exact mathematical RH jump.
//   Transport fields use exactly one acceleration representation: SOURCE keeps
//   SHOCK_ONLY background fields uncompressed, whereas RESOLVED_COMPRESSION maps
//   upstream to the exact RH downstream state across a finite C1 shock layer.
//
// • Artificial leading/trailing ICME boundaries use symmetric C1 smoothstep
//   transitions.  sheath_comp_floor is retained only for source compatibility;
//   it never alters either the exact RH state or the resolved transport profile.


 NUMERICAL NOTES
   • Explicit numerical status: invalid/out-of-domain/non-finite samples are
    reported and are never replaced by plausible fallback values.
   • No heap allocs in evaluators; no pow() in hot paths.
   • ∇·V uses the exact radial identity 2V/r+dV/dr; dV/dr is analytical.
     `dr_frac` remains in the legacy public signature but is ignored in 1-D.
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
     .SetCME(20.0, 1800, 8e-8)                      // r0[R_sun], V0_sh[km/s], Γ[1/km]
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
     r0_Rs      ~15–20       Recommended DBM start radius in drag-dominated heliosphere
     V0_sh_kms  800–2500     Initial shock speed
     Gamma_kmInv 1e-8–2e-7   Drag Γ (↑ ⇒ stronger decel)

   Geometry (at 1 AU; scales ∝ R_sh):
     sheath_thick_AU_at1AU   0.05–0.20
     ejecta_thick_AU_at1AU   0.10–0.40

   Edge widths (at 1 AU; scales ∝ R_sh):
     edge_smooth_shock_AU_at1AU  0.005–0.02; resolved-compression layer width
     edge_smooth_le_AU_at1AU     0.01–0.05
     edge_smooth_te_AU_at1AU     0.02–0.06

   Sheath/ME targets:
     sheath_comp_floor    legacy compatibility input; ignored by physical shock/region state
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

#include <atomic>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <algorithm>
#include <limits>
#include <stdexcept>
#include <vector>
#include <string>
#include <sstream>
#include <iomanip>
#include <type_traits>
#include <utility>

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

namespace swcme1d {

// --------------------------- Physical constants (SI) ---------------------------
constexpr double PI        = swcme::constants::PI;
constexpr double AU        = swcme::constants::AU_M;                    // [m]
constexpr double Rs        = swcme::constants::SOLAR_RADIUS_M;          // [m]
constexpr double OMEGA_SUN = swcme::defaults::SOLAR_ROTATION_RATE_RAD_S;    // [rad/s]
constexpr double MU0       = swcme::constants::VACUUM_PERMEABILITY_N_A2; // [N/A²]
constexpr double MP        = swcme::constants::PROTON_MASS_KG;          // [kg]
constexpr double KB        = swcme::constants::BOLTZMANN_J_K;           // [J/K]

// --------------------------------- Helpers -----------------------------------
inline double clamp01(double x){ return x<0.0?0.0:(x>1.0?1.0:x); }
inline double clamp(double x,double a,double b){ return x<a?a:(x>b?b:x); }
inline double smoothstep01(double x){ x=clamp01(x); return x*x*(3.0-2.0*x); } // C¹
inline double lerp(double a,double b,double t){ return a + (b-a)*t; }

// ------------------------------ User parameters ------------------------------
/**
 * @brief Tunable physical and geometric parameters of the model.
 *
 * Ambient inputs (V_sw, n1AU, |B|1AU, T) set the Parker spiral and upstream
 * thermodynamics. CME inputs feed the shared ballistic/DBM/data-driven apex kinematics.
 * Geometry and smoothing control sheath/ME sizes and edge widths (all scale
 * self-similarly with R_sh).  The selected acceleration mode determines whether
 * the shock is an explicit SOURCE surface or a resolved C1 compression layer.
 * Sheath/ME shaping controls post-shock relaxation and ejecta density/speed.
 *
 * Units noted per field; high-level units: [km/s], [cm^-3], [nT], [K], [AU].
 */
struct Params {
  // Ambient & thermodynamics
  double V_sw_kms    = swcme::defaults::V_SW_KMS;  // upstream wind speed [km/s]
  double n1AU_cm3    = swcme::defaults::N1AU_CM3;    // density at 1 AU [cm⁻³]
  double B1AU_nT     = swcme::defaults::B1AU_TOTAL_NT;    // |B|(1 AU) [nT]
  double T_K         = swcme::defaults::T_K;  // proton temperature [K]
  double gamma_ad    = swcme::defaults::GAMMA_AD;
  double sin_theta   = swcme::defaults::PARKER_REFERENCE_SIN_THETA; // fixed 1-D ray latitude / Parker normalization

  // CME/shock-apex kinematics.  Both 1-D and 3-D now use the shared
  // swcme::kinematics implementation, so selecting the same mode and inputs
  // produces exactly the same apex radius/speed in both models.
  swcme::kinematics::Mode kinematics_mode = swcme::defaults::KINEMATICS_MODE;
  double r0_Rs       = swcme::defaults::DBM_R0_RS; // DBM/ballistic reference radius [R☉]
  double V0_sh_kms   = swcme::defaults::V0_SH_KMS; // initial/reference apex speed [km/s]
  double Gamma_kmInv = swcme::defaults::DBM_GAMMA_KM_INV; // DBM drag coefficient [1/km], must be >=0

  // DATA_DRIVEN mode: monotonically increasing height-time knots.  Times are
  // seconds in the same launch-relative clock used by prepare_step(); radii
  // are solar radii.  The default OUTSIDE_TIME policy deliberately refuses to
  // extrapolate beyond the observed interval unless the caller explicitly
  // requests ballistic endpoint continuation.
  std::vector<double> data_time_s;
  std::vector<double> data_radius_Rs;
  swcme::kinematics::ExtrapolationPolicy data_extrapolation =
      swcme::defaults::DATA_EXTRAPOLATION;

  // Region mode.  The canonical default is SHOCK_ONLY because the validated
  // science baseline is an upstream Parker/Leblanc SEP experiment. FULL_ICME
  // is an explicit phenomenological diagnostic override.
  swcme::regions::Mode region_mode = swcme::defaults::REGION_MODE;

  // Shock acceleration representation.  Exactly one mechanism is selected.
  // SOURCE is the canonical controlled-SEP default and is paired with
  // SHOCK_ONLY; RESOLVED_COMPRESSION is an explicit FULL_ICME diagnostic mode.
  // The relative source weight is dimensionless
  // until the later AMPS source adapter assigns a physical injection unit.
  swcme::acceleration::Mode shock_acceleration_mode =
      swcme::defaults::ACCELERATION_MODE;
  double relative_source_weight_per_area = swcme::defaults::RELATIVE_SOURCE_WEIGHT_PER_AREA;

  // Geometry: thicknesses at 1 AU. These values are interpreted as
  // self-similar fractions of the local shock radius by swcme_regions.hpp.
  double sheath_thick_AU_at1AU  = swcme::defaults::SHEATH_THICK_AU_AT_1AU; // AU at 1 AU
  double ejecta_thick_AU_at1AU  = swcme::defaults::EJECTA_THICK_AU_AT_1AU; // AU at 1 AU

  // Interface smoothing widths at 1 AU (C¹), scale ∝ R_sh.  The shock width
  // is active only for RESOLVED_COMPRESSION; SOURCE mode validates with
  // SHOCK_ONLY and therefore never places this compression in the flow field.
  double edge_smooth_shock_AU_at1AU = swcme::defaults::EDGE_SMOOTH_SHOCK_AU_AT_1AU;
  double edge_smooth_le_AU_at1AU    = swcme::defaults::EDGE_SMOOTH_LE_AU_AT_1AU; // sheath → ME
  double edge_smooth_te_AU_at1AU    = swcme::defaults::EDGE_SMOOTH_TE_AU_AT_1AU; // ME → ambient

  // Sheath / ME shaping
  // Deprecated compatibility parameter retained for source compatibility.
  // It is ignored by both shock and region physics; the MHD RH solver alone
  // determines physical compression.
  double sheath_comp_floor   = swcme::defaults::SHEATH_COMP_FLOOR_COMPAT;
  double sheath_ramp_power   = swcme::defaults::SHEATH_RAMP_POWER; // controls steepness near shock (≥1)
  double V_sheath_LE_factor  = swcme::defaults::V_SHEATH_LE_FACTOR; // V at LE relative to V_sw (≥1)
  double f_ME                = swcme::defaults::F_ME; // ME density factor vs upstream (<1 typical)
  double V_ME_factor         = swcme::defaults::V_ME_FACTOR; // ME speed factor vs V_sw (<1 typical)
};

// Deterministic, complete resolved-configuration record.  This is intentionally
// plain key=value text so validation/campaign tooling can archive the exact
// model inputs without depending on a JSON/YAML library.  Every public Params
// field is emitted, including currently inactive FULL_ICME parameters and the
// data-driven tables, because an override must never disappear from run
// metadata merely because another mode makes it inactive for this run.
inline std::string resolved_configuration_manifest(const Params& p) {
  std::ostringstream out;
  out << std::setprecision(17) << std::scientific;
  out << "swcme_config_version=" << swcme::defaults::CONFIG_VERSION << '\n';
  out << "model=1D\n";
  out << "frame=" << swcme::defaults::FRAME_NAME << '\n';
  out << "model_scope=" << swcme::defaults::model_scope_name(
      swcme::defaults::model_scope(p.region_mode, p.shock_acceleration_mode)) << '\n';
  out << "parker_normalization="
      << swcme::defaults::PARKER_NORMALIZATION_CONVENTION << '\n';
  out << "parker_radial_polarity=" << swcme::defaults::PARKER_RADIAL_POLARITY << '\n';
  out << "solar_rotation_rate_rad_s="
      << swcme::defaults::SOLAR_ROTATION_RATE_RAD_S << '\n';
  out << "V_sw_kms=" << p.V_sw_kms << '\n';
  out << "n1AU_cm3=" << p.n1AU_cm3 << '\n';
  out << "B1AU_nT=" << p.B1AU_nT << '\n';
  out << "T_K=" << p.T_K << '\n';
  out << "gamma_ad=" << p.gamma_ad << '\n';
  out << "sin_theta=" << p.sin_theta << '\n';
  out << "kinematics_mode=" << swcme::defaults::kinematics_mode_name(p.kinematics_mode) << '\n';
  out << "r0_Rs=" << p.r0_Rs << '\n';
  out << "V0_sh_kms=" << p.V0_sh_kms << '\n';
  out << "Gamma_kmInv=" << p.Gamma_kmInv << '\n';
  out << "data_extrapolation="
      << swcme::defaults::extrapolation_policy_name(p.data_extrapolation) << '\n';
  out << "data_time_s.count=" << p.data_time_s.size() << '\n';
  for (std::size_t i=0; i<p.data_time_s.size(); ++i)
    out << "data_time_s[" << i << "]=" << p.data_time_s[i] << '\n';
  out << "data_radius_Rs.count=" << p.data_radius_Rs.size() << '\n';
  for (std::size_t i=0; i<p.data_radius_Rs.size(); ++i)
    out << "data_radius_Rs[" << i << "]=" << p.data_radius_Rs[i] << '\n';
  out << "region_mode=" << swcme::defaults::region_mode_name(p.region_mode) << '\n';
  out << "shock_acceleration_mode="
      << swcme::defaults::acceleration_mode_name(p.shock_acceleration_mode) << '\n';
  out << "relative_source_weight_per_area=" << p.relative_source_weight_per_area << '\n';
  out << "sheath_thick_AU_at1AU=" << p.sheath_thick_AU_at1AU << '\n';
  out << "ejecta_thick_AU_at1AU=" << p.ejecta_thick_AU_at1AU << '\n';
  out << "edge_smooth_shock_AU_at1AU=" << p.edge_smooth_shock_AU_at1AU << '\n';
  out << "edge_smooth_le_AU_at1AU=" << p.edge_smooth_le_AU_at1AU << '\n';
  out << "edge_smooth_te_AU_at1AU=" << p.edge_smooth_te_AU_at1AU << '\n';
  out << "sheath_comp_floor=" << p.sheath_comp_floor << '\n';
  out << "sheath_ramp_power=" << p.sheath_ramp_power << '\n';
  out << "V_sheath_LE_factor=" << p.V_sheath_LE_factor << '\n';
  out << "f_ME=" << p.f_ME << '\n';
  out << "V_ME_factor=" << p.V_ME_factor << '\n';
  return out.str();
}

// Return the PST03 fingerprint of the complete resolved 1-D configuration.
// The schema tag intentionally makes field order part of a versioned contract:
// adding or reinterpreting a parameter must update the tag, preventing an old
// prepared state from being accepted under new physics.  Compile-time Parker
// polarity/normalization and the proton-only pressure closure are included
// because they influence results even though they are not runtime Params.
inline swcme::ConfigurationDigest configuration_digest(
    const Params& p) noexcept {
  swcme::ConfigurationDigestBuilder digest;
  digest.add_string("SWCME_CONFIGURATION_DIGEST_V1");
  digest.add_string("1D");
  digest.add_uint64(static_cast<std::uint64_t>(swcme::defaults::CONFIG_VERSION));
  digest.add_string(swcme::defaults::FRAME_NAME);
  digest.add_string(swcme::defaults::PARKER_NORMALIZATION_CONVENTION);
  digest.add_uint64(static_cast<std::uint64_t>(
      static_cast<std::int64_t>(swcme::defaults::PARKER_RADIAL_POLARITY)));
  digest.add_string("PROTON_ONLY_THERMAL_PRESSURE_CLOSURE");
  digest.add_double(swcme::defaults::SOLAR_ROTATION_RATE_RAD_S);

  // Hash every public Params field in declaration order.  Vector lengths are
  // included before their values so tables with identical prefixes remain
  // distinguishable and default-empty equals explicitly-empty configuration.
  digest.add_double(p.V_sw_kms); digest.add_double(p.n1AU_cm3);
  digest.add_double(p.B1AU_nT); digest.add_double(p.T_K);
  digest.add_double(p.gamma_ad); digest.add_double(p.sin_theta);
  digest.add_uint64(static_cast<std::uint64_t>(p.kinematics_mode));
  digest.add_double(p.r0_Rs); digest.add_double(p.V0_sh_kms);
  digest.add_double(p.Gamma_kmInv);
  digest.add_uint64(static_cast<std::uint64_t>(p.data_time_s.size()));
  for (double value : p.data_time_s) digest.add_double(value);
  digest.add_uint64(static_cast<std::uint64_t>(p.data_radius_Rs.size()));
  for (double value : p.data_radius_Rs) digest.add_double(value);
  digest.add_uint64(static_cast<std::uint64_t>(p.data_extrapolation));
  digest.add_uint64(static_cast<std::uint64_t>(p.region_mode));
  digest.add_uint64(static_cast<std::uint64_t>(p.shock_acceleration_mode));
  digest.add_double(p.relative_source_weight_per_area);
  digest.add_double(p.sheath_thick_AU_at1AU);
  digest.add_double(p.ejecta_thick_AU_at1AU);
  digest.add_double(p.edge_smooth_shock_AU_at1AU);
  digest.add_double(p.edge_smooth_le_AU_at1AU);
  digest.add_double(p.edge_smooth_te_AU_at1AU);
  digest.add_double(p.sheath_comp_floor);
  digest.add_double(p.sheath_ramp_power);
  digest.add_double(p.V_sheath_LE_factor);
  digest.add_double(p.f_ME); digest.add_double(p.V_ME_factor);
  return digest.value();
}

// Validate the complete 1-D public parameter bundle before any unit conversion
// or physics evaluation.  The common rules live in swcme_config.hpp so the
// equivalent 1-D and 3-D fields are judged by exactly the same contract.
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
  return swcme::config::validate_common(view);
}

// ------------------------------ Per‑time cache -------------------------------
/**
 * @brief Per-time cache. Construct once via prepare_step(t) and reuse for
 *        many radius queries. Keeps all heavy computations out of the hot path.
 *
 * Contains: DBM kinematics & geometry (R_sh, R_LE, R_TE, widths), Parker
 * constants, Leblanc coefficients scaled to match n(1 AU), shock compression
 * ratio, explicit shock-existence state, and *pinned* RH boundary values used to build a strictly monotone
 * sheath (n_up at shock & LE; V2 at shock; V at LE).
 */
struct StepState {
  // PST05 lifetime contract: this record owns every value needed by checked
  // evaluation and never stores a pointer/reference to its preparing Model.
  // It may therefore be copied, moved, inspected, or destroyed after that
  // Model's lifetime ends.  Numerical evaluation still requires a live Model
  // carrying the same logical identity; a newly constructed equal Model has a
  // different identity and deterministically rejects this orphaned snapshot.
  // Identity of the exact Model instance that prepared this cache.  Zero means
  // that the object was default-constructed and was never prepared.  The
  // identity is metadata only; it is checked before physics and is never used
  // to select a numerical branch or alter a physical result.
  swcme::ModelIdentity owner_model_identity = 0;

  // Snapshot of the complete resolved configuration used by prepare_step().
  // It is diagnostic metadata and is checked before any cached physics or
  // caller-owned output is touched.
  swcme::ConfigurationDigest configuration_digest = 0;

  // Canonical dimensionality-independent prepared state.  The fields below
  // mirror selected values for backward/source compatibility with existing
  // callers, but ambient normalization and apex kinematics are computed only
  // once by swcme::core::prepare().
  swcme::core::PreparedState common;

  // Shared phenomenological region configuration and the apex/radial boundary
  // set.  The same swcme::regions contract is used by 3-D, so 1-D and 3-D no
  // longer maintain independent interpretations of sheath/ejecta thicknesses
  // or LE/TE smoothing widths.
  swcme::regions::Config region_config;
  swcme::regions::Boundaries region_boundaries;

  // Shared shock-acceleration choice cached with the step so field/source
  // queries cannot observe a mode different from the one that was validated.
  swcme::acceleration::Config acceleration_config;

  // Apex kinematics / geometry.  kinematics_mode records which shared common
  // solver produced r_sh_m and V_sh_ms for traceable diagnostics.
  swcme::kinematics::Mode kinematics_mode = swcme::defaults::KINEMATICS_MODE;
  double time_s    = 0.0;     // time since launch [s]
  double r0_m      = 20.0*Rs; // default DBM reference radius [m]
  double r_sh_m    = 30*AU;   // shock apex radius [m]
  double V_sh_ms   = 1.0e6;   // shock speed [m/s]
  double r_le_m    = 29*AU;   // leading edge [m]
  double r_te_m    = 28*AU;   // trailing edge [m]
  double w_sh_m    = 0.0;     // smoothing widths [m]
  double w_le_m    = 0.0;
  double w_te_m    = 0.0;

  // Ambient / Parker / Leblanc
  double V_up_ms   = 4.0e5;   // upstream wind speed [m/s]
  double Br1AU_T   = 0.0;     // Br at 1 AU [T]
  double k_AU      = 0.0;     // Ω AU sinθ / V_sw [‑]
  double B_up_T    = 0.0;     // |B| at R_sh upstream [T]
  double C2        = 0.0;     // Leblanc scaled SI coefficients: n=C2/r²+C4/r⁴+C6/r⁶
  double C4        = 0.0;
  double C6        = 0.0;

  // Shock compression and cached boundary values.  has_shock is determined
  // from the shared ideal-MHD fast-shock criterion before any sheath shaping
  // is applied; sheath_comp_floor can no longer manufacture a shock.
  bool has_shock   = false;
  bool shock_solver_converged = true;
  // Preserve the complete shared solver result so cross-dimensional tests and
  // future source coupling can compare the exact same MHD state instead of
  // reconstructing it from a few scalar mirrors.
  swcme::shock::JumpResult shock_jump;
  double rc        = 1.0;     // physical density compression; exactly 1 if no shock
  double n_up_shock = 0.0;    // upstream density at shock radius [m⁻3]
  double n_up_le    = 0.0;    // upstream density at leading edge [m⁻3]
  double V2_shock_ms = 0.0;   // immediate downstream speed at the shock [m/s]
  double V_LE_ms     = 0.0;   // sheath speed at the leading edge [m/s]

  // The seal is intentionally the only private StepState datum.  Public fields
  // above remain transitional compatibility mirrors, but callers cannot
  // recompute/reseal a modified record.  PST06 validation therefore converts
  // any mirror corruption into an explicit failure before physics is used.
  swcme::ConfigurationDigest integrity_digest() const noexcept {
    return integrity_digest_;
  }

 private:
  swcme::ConfigurationDigest integrity_digest_ = 0;
  friend class Model;
};

// Serialize every canonical and compatibility field in a fixed order.  This
// function deliberately excludes the private seal itself; prepare_step() stores
// the returned value and consumers independently recompute it.
inline swcme::ConfigurationDigest prepared_state_integrity(
    const StepState& state) noexcept {
  swcme::ConfigurationDigestBuilder digest(false);
  digest.add_string("SWCME_1D_PREPARED_STATE_INTEGRITY_V1");
  digest.add_uint64(state.owner_model_identity);
  digest.add_uint64(state.configuration_digest);
  swcme::prepared_integrity::add_common(digest,state.common);
  swcme::prepared_integrity::add_region_config(digest,state.region_config);
  swcme::prepared_integrity::add_boundaries(digest,state.region_boundaries);
  swcme::prepared_integrity::add_acceleration_config(
      digest,state.acceleration_config);
  digest.add_uint64(static_cast<std::uint64_t>(state.kinematics_mode));
  digest.add_double(state.time_s);
  digest.add_double(state.r0_m);
  digest.add_double(state.r_sh_m);
  digest.add_double(state.V_sh_ms);
  digest.add_double(state.r_le_m);
  digest.add_double(state.r_te_m);
  digest.add_double(state.w_sh_m);
  digest.add_double(state.w_le_m);
  digest.add_double(state.w_te_m);
  digest.add_double(state.V_up_ms);
  digest.add_double(state.Br1AU_T);
  digest.add_double(state.k_AU);
  digest.add_double(state.B_up_T);
  digest.add_double(state.C2);
  digest.add_double(state.C4);
  digest.add_double(state.C6);
  digest.add_bool(state.has_shock);
  digest.add_bool(state.shock_solver_converged);
  swcme::prepared_integrity::add_jump(digest,state.shock_jump);
  digest.add_double(state.rc);
  digest.add_double(state.n_up_shock);
  digest.add_double(state.n_up_le);
  digest.add_double(state.V2_shock_ms);
  digest.add_double(state.V_LE_ms);
  return digest.value();
}

// --------------------------------- Model -------------------------------------
/**
 * @brief Main model class. Create a Model, set Params (or use defaults), then
 *        call prepare_step(t) to get a StepState. Use evaluators to sample n,V
 *        (and B, ∇·V) at arbitrary radii.
 *
 * Threading: read-only after prepare_step(); safe to call evaluators from many
 * threads with the same StepState.
 */
class Model {
public:
  Model() : P{}, model_identity_(swcme::next_model_identity()),
            configuration_digest_(swcme1d::configuration_digest(P)),
            configuration_locked_(false) {}
  explicit Model(const Params& p)
      : P(p), model_identity_(swcme::next_model_identity()),
        configuration_digest_(swcme1d::configuration_digest(P)),
        configuration_locked_(false) {}

  // A copied Model is a new owner even when its Params are identical.  Giving
  // the copy a fresh identity prevents a StepState produced by the source
  // object from being accepted by the copy merely because the compiler copied
  // a hidden token together with the public configuration.
  Model(const Model& other)
      : P(other.P), model_identity_(swcme::next_model_identity()),
        configuration_digest_(other.configuration_digest_),
        configuration_locked_(false) {}
  Model& operator=(const Model& other) {
    if (this!=&other) {
      require_configuration_mutable("operator=");
      P=other.P;
      // Assignment is allowed only before this object has issued a state, but
      // it still changes the logical model and therefore receives a fresh
      // identity consistent with PST02's instance-provenance contract.
      model_identity_=swcme::next_model_identity();
      configuration_digest_=other.configuration_digest_;
      configuration_locked_.store(false,std::memory_order_release);
    }
    return *this;
  }

  // Moving transfers the logical owner, not merely the numerical Params.
  // Consequently, every StepState prepared before the move remains valid with
  // the destination object.  The source is rotated to a fresh identity and
  // unlocked, so it cannot consume those states and remains safe to destroy or
  // assign.  Params uses standard value members (including default-allocator
  // vectors), making its move construction non-throwing; that guarantee lets
  // std::vector relocate prepared Models without falling back to copy semantics
  // (which intentionally create a distinct PST02 owner).
  Model(Model&& other) noexcept(
      std::is_nothrow_move_constructible<Params>::value)
      : P(std::move(other.P)), model_identity_(other.model_identity_),
        configuration_digest_(other.configuration_digest_),
        configuration_locked_(
            other.configuration_locked_.load(std::memory_order_acquire)) {
    other.model_identity_=swcme::next_model_identity();
    other.configuration_digest_=swcme1d::configuration_digest(other.P);
    other.configuration_locked_.store(false,std::memory_order_release);
  }

  // Move assignment follows PST01 as well as PST05.  A destination that has
  // already issued a state is immutable and is rejected before either object
  // changes.  An unprepared destination assumes the source identity and lock
  // state, while the moved-from source receives a new, unprepared lifetime.
  Model& operator=(Model&& other) {
    if (this!=&other) {
      require_configuration_mutable("move operator=");
      P=std::move(other.P);
      model_identity_=other.model_identity_;
      configuration_digest_=other.configuration_digest_;
      configuration_locked_.store(
          other.configuration_locked_.load(std::memory_order_acquire),
          std::memory_order_release);
      other.model_identity_=swcme::next_model_identity();
      other.configuration_digest_=swcme1d::configuration_digest(other.P);
      other.configuration_locked_.store(false,std::memory_order_release);
    }
    return *this;
  }

  // Parameter setters (fluent)
  Model& SetParams(const Params& p){
    require_configuration_mutable("SetParams"); P=p;
    refresh_configuration_digest(); return *this; }
  Model& SetCME(double r0_Rs,double V0_sh_kms,double Gamma_kmInv){
    require_configuration_mutable("SetCME");
    P.r0_Rs=r0_Rs; P.V0_sh_kms=V0_sh_kms; P.Gamma_kmInv=Gamma_kmInv;
    refresh_configuration_digest(); return *this; }
  Model& SetKinematicsMode(swcme::kinematics::Mode mode){
    require_configuration_mutable("SetKinematicsMode");
    P.kinematics_mode=mode; refresh_configuration_digest(); return *this; }
  Model& SetDataDrivenKinematics(const std::vector<double>& time_s,
                                 const std::vector<double>& radius_Rs,
                                 swcme::kinematics::ExtrapolationPolicy policy=
                                     swcme::defaults::DATA_EXTRAPOLATION){
    require_configuration_mutable("SetDataDrivenKinematics");
    P.kinematics_mode=swcme::kinematics::Mode::DataDriven;
    P.data_time_s=time_s; P.data_radius_Rs=radius_Rs; P.data_extrapolation=policy;
    refresh_configuration_digest();
    return *this; }
  Model& SetAmbient(double V_sw_kms,double n1AU_cm3,double B1AU_nT,double T_K,
                    double gamma_ad=swcme::defaults::GAMMA_AD,
                    double sin_theta=swcme::defaults::PARKER_REFERENCE_SIN_THETA){
    require_configuration_mutable("SetAmbient");
    P.V_sw_kms=V_sw_kms; P.n1AU_cm3=n1AU_cm3; P.B1AU_nT=B1AU_nT; P.T_K=T_K;
    P.gamma_ad=gamma_ad; P.sin_theta=sin_theta;
    refresh_configuration_digest(); return *this; }
  Model& SetRegionMode(swcme::regions::Mode mode){
    require_configuration_mutable("SetRegionMode");
    P.region_mode=mode; refresh_configuration_digest(); return *this; }
  Model& SetShockAccelerationMode(swcme::acceleration::Mode mode){
    require_configuration_mutable("SetShockAccelerationMode");
    P.shock_acceleration_mode=mode; refresh_configuration_digest(); return *this; }
  Model& SetGeometry(double sheath_thick_AU_at1AU,double ejecta_thick_AU_at1AU){
    require_configuration_mutable("SetGeometry");
    P.sheath_thick_AU_at1AU=sheath_thick_AU_at1AU;
    P.ejecta_thick_AU_at1AU=ejecta_thick_AU_at1AU;
    refresh_configuration_digest(); return *this; }
  Model& SetSmoothing(double w_sh,double w_le,double w_te){
    require_configuration_mutable("SetSmoothing");
    P.edge_smooth_shock_AU_at1AU=w_sh;
    P.edge_smooth_le_AU_at1AU   =w_le;
    P.edge_smooth_te_AU_at1AU   =w_te;
    refresh_configuration_digest(); return *this; }
  Model& SetSheathEjecta(double sheath_comp_floor,double sheath_ramp_power,
                         double V_sheath_LE_factor,double f_ME,double V_ME_factor){
    require_configuration_mutable("SetSheathEjecta");
    P.sheath_comp_floor=sheath_comp_floor; P.sheath_ramp_power=sheath_ramp_power;
    P.V_sheath_LE_factor=V_sheath_LE_factor; P.f_ME=f_ME;
    P.V_ME_factor=V_ME_factor; refresh_configuration_digest();
    return *this; }

  const Params& GetParams() const { return P; }
  // A raw Params& could be retained before prepare_step() and used afterward,
  // bypassing every runtime lifecycle check.  PST01 therefore makes the
  // legacy escape hatch explicitly unavailable.  Use the guarded setters in
  // the setup phase, or copy GetParams() and call reconfigured(params).
  Params& MutableParams() = delete;

  // Construct a replacement model instead of reopening a frozen instance.
  // The returned object has a fresh PST02 owner identity and remains mutable
  // until its own first successful prepare_step().  This is the supported
  // builder path for parameter sweeps and event-to-event reconfiguration.
  Model reconfigured(const Params& p) const { return Model(p); }

  // Expose the lifecycle state without exposing the mutable flag itself.
  // Acquire ordering pairs with the successful prepare_step() release store.
  bool configuration_locked() const noexcept {
    return configuration_locked_.load(std::memory_order_acquire);
  }

  swcme::ModelIdentity model_identity() const noexcept {
    return model_identity_;
  }

  // Validate ownership without touching any physics output.  Public adapters
  // call this guard before clearing their destination objects, which preserves
  // caller sentinels and makes a rejected mixed-model call transactional.
  swcme::ModelStatus validate_prepared_state(
      const StepState& S, const char* context) const noexcept {
    // Params is mutable only during setup, where each setter refreshes this
    // cache.  Successful preparation freezes both values, so hot validation
    // reads one scalar instead of serializing the complete configuration.
    const swcme::ConfigurationDigest current=configuration_digest_;
    // Preserve PST02 precedence: a foreign instance is always an ownership
    // error, while its digests explain whether the configurations also differ.
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

  // Public validation entry point used by CFG01 and by prepare_step().  It is
  // intentionally side-effect free so callers can inspect a configuration
  // before launching a time-dependent calculation.
  swcme::config::ValidationResult validate() const { return validate_params(P); }

  // The scope is derived from the validated region/acceleration pair; it is
  // not stored independently and therefore cannot drift out of sync.
  swcme::defaults::ModelScope model_scope() const {
    return swcme::defaults::model_scope(P.region_mode, P.shock_acceleration_mode);
  }

  // Report whether the model is being used within its observer-local declared
  // science scope.  In the default CONTROLLED_SEP_PRE_SHOCK mode the Parker
  // background is valid only while the shock radius remains below the observer
  // radius. FULL_ICME_DIAGNOSTIC remains mathematically evaluable across the
  // modeled regions but is not promoted to a validated global ICME model.
  swcme::defaults::ObserverScopeStatus observer_scope_status(
      const StepState& S, double observer_radius_m) const {
    // This value-returning compatibility API has no ModelStatus channel, so a
    // foreign state is surfaced as an exception before cached radius data are
    // combined with this model's region/acceleration parameters.
    swcme::throw_if_error(validate_prepared_state(
        S,"swcme1d::observer_scope_status"));
    return swcme::defaults::observer_scope_status(
        P.region_mode, P.shock_acceleration_mode, true,
        S.r_sh_m, observer_radius_m);
  }

  // ---------------------------- Build per‑time cache -------------------------
  /**
   * @brief Build a per-time cache (StepState) at time t since CME launch.
   *
   * Steps:
   *  1) Parse ambient inputs; compute Parker constants k and Br(1 AU) from |B|1AU.
   *  2) Scale Leblanc coefficients to match n(1 AU) in SI form n=C2/r^2+C4/r^4+C6/r^6.
   *  3) Evaluate the shared apex kinematics mode (BALLISTIC/DBM/DATA_DRIVEN).
   *  4) Build self-similar geometry and the already validated smoothing widths.
   *  5) Evaluate upstream n and |B| at R_sh; compute c_s, v_A, c_f; estimate r_c.
   *  6) Cache *boundary values* for monotone sheath: n_up(R_sh), n_up(R_LE),
   *     V2 (RH proxy), and V_LE ≥ V_sw.
   */
  StepState prepare_step(double t_s) const {
    // Configuration errors are rejected before any normalization, boundary
    // construction, or unit conversion can hide the supplied value.  CFG01
    // covers scalar admissibility and CFG03 covers smoothing/layer conflicts;
    // invalid input fails once at setup instead of becoming a plausible state.
    const swcme::config::ValidationResult validation=validate();
    if (!validation.ok()) {
      throw std::invalid_argument(validation.summary("swcme1d"));
    }
    if (!std::isfinite(t_s) || t_s<0.0) {
      throw std::invalid_argument("swcme1d: time must be finite and >= 0");
    }

    StepState S;
    // Stamp ownership before filling the expensive cache.  prepare_step()
    // either throws and returns no state, or returns a fully initialized state
    // that can be consumed only by this exact Model instance.
    S.owner_model_identity=model_identity_;
    // Capture the exact configuration only after validation succeeds.  PST01
    // freezes supported mutation when this preparation returns; the digest
    // remains a defense-in-depth check before any consumer modifies outputs.
    S.configuration_digest=configuration_digest_;
    S.time_s=t_s;

    // Build the dimensionality-independent core configuration in public units
    // and prepare it once.  This single call now owns unit conversion, Leblanc
    // normalization, Parker radial-field normalization, and apex kinematics.
    // The 1-D wrapper only supplies its fixed reference latitude and later
    // applies 1-D region geometry.
    swcme::core::CommonConfig common_cfg;
    common_cfg.V_sw_kms=P.V_sw_kms;
    common_cfg.n1AU_cm3=P.n1AU_cm3;
    common_cfg.B1AU_nT=P.B1AU_nT;
    common_cfg.T_K=P.T_K;
    common_cfg.gamma_ad=P.gamma_ad;
    common_cfg.parker_reference_sin_theta=P.sin_theta;
    common_cfg.solar_rotation_rate_rad_s=OMEGA_SUN;
    common_cfg.kinematics_mode=P.kinematics_mode;
    common_cfg.r0_Rs=P.r0_Rs;
    common_cfg.V0_sh_kms=P.V0_sh_kms;
    common_cfg.Gamma_kmInv=P.Gamma_kmInv;
    common_cfg.data_time_s=P.data_time_s;
    common_cfg.data_radius_Rs=P.data_radius_Rs;
    common_cfg.data_extrapolation=P.data_extrapolation;

    S.common=swcme::core::prepare(common_cfg,t_s);
    if (S.common.apex.status!=swcme::kinematics::Status::Ok) {
      throw std::runtime_error(std::string("swcme1d kinematics: ")+
                               swcme::kinematics::status_name(S.common.apex.status));
    }
    if (!std::isfinite(S.common.apex.radius_m) ||
        S.common.apex.radius_m < swcme::solarwind::MIN_RADIUS_M) {
      throw std::runtime_error(
          swcme::ModelStatus::make_value(
              swcme::StatusCode::OutsideModelDomain,
              "swcme1d::prepare_step shock radius", S.common.apex.radius_m)
              .summary());
    }

    // Legacy StepState fields are mirrors only.  Keeping them populated avoids
    // breaking existing callers while guaranteeing that both dimensional
    // models receive identical common-core values.
    S.r0_m=S.common.r0_m;
    S.kinematics_mode=P.kinematics_mode;
    S.V_up_ms=S.common.solar_wind.V_sw_m_s;
    S.Br1AU_T=S.common.solar_wind.Br1AU_T;
    S.k_AU=S.common.solar_wind.k_AU_equatorial*P.sin_theta;
    S.C2=S.common.solar_wind.C2;
    S.C4=S.common.solar_wind.C4;
    S.C6=S.common.solar_wind.C6;
    S.r_sh_m=S.common.apex.radius_m;
    S.V_sh_ms=S.common.apex.speed_m_s;
    const double Vsw=S.V_up_ms;

    // Build the common self-similar sheath/ejecta geometry.  Public thickness
    // and smoothing inputs are AU at a 1-AU shock, i.e. dimensionless fractions
    // of the local shock radius.  The same routine is used by 3-D at every
    // flank direction, eliminating the former apex-width subtraction.
    S.acceleration_config.mode=P.shock_acceleration_mode;
    S.acceleration_config.relative_source_weight_per_area=
        P.relative_source_weight_per_area;
    S.region_config.mode=P.region_mode;
    // The physical shock smoothing width belongs to the selected acceleration
    // representation. SOURCE has no resolved compression layer; the shock is
    // an injection surface only. RESOLVED_COMPRESSION uses the same common C1
    // width in both dimensional models.
    S.region_config.shock_smooth_fraction=
        (P.shock_acceleration_mode==swcme::acceleration::Mode::ResolvedCompression)
            ? P.edge_smooth_shock_AU_at1AU : 0.0;
    S.region_config.sheath_fraction=P.sheath_thick_AU_at1AU;
    S.region_config.ejecta_fraction=P.ejecta_thick_AU_at1AU;
    S.region_config.leading_smooth_fraction=P.edge_smooth_le_AU_at1AU;
    S.region_config.trailing_smooth_fraction=P.edge_smooth_te_AU_at1AU;
    S.region_config.sheath_ramp_power=P.sheath_ramp_power;
    S.region_config.V_sheath_LE_factor=P.V_sheath_LE_factor;
    S.region_config.f_ME=P.f_ME;
    S.region_config.V_ME_factor=P.V_ME_factor;
    S.region_boundaries=swcme::regions::make_boundaries(S.r_sh_m,S.region_config);

    // Legacy StepState mirrors remain populated for source compatibility and
    // Tecplot output. In RESOLVED_COMPRESSION w_sh_m is the common total C1
    // transition width; SOURCE sets it exactly to zero and carries acceleration
    // only through the explicit source record.
    S.r_le_m=S.region_boundaries.R_le_m;
    S.r_te_m=S.region_boundaries.R_te_m;
    S.w_sh_m=S.region_boundaries.smooth_shock_width_m;
    S.w_le_m=S.region_boundaries.smooth_le_width_m;
    S.w_te_m=S.region_boundaries.smooth_te_width_m;

    // Upstream Parker field at the shock.  The 1-D ray is treated as the
    // local shock normal, while B_phi remains a tangential component.  This
    // vector decomposition lets the same ideal-MHD jump solver used by 3-D
    // determine both shock existence and the downstream state.
    const swcme::solarwind::ParkerComponents parker_sh =
        swcme::solarwind::parker_components(
            S.common.solar_wind,S.r_sh_m,P.sin_theta);
    const double Br_sh=parker_sh.Br_T;
    const double Bphi_sh=parker_sh.Bphi_T;
    S.B_up_T=parker_sh.Bmag_T;

    const double n_up_sh = density_upstream(S, S.r_sh_m);
    swcme::shock::PrimitiveState upstream;
    upstream.rho_kg_m3 = std::max(0.0,n_up_sh)*MP;
    upstream.pressure_Pa = swcme::solarwind::proton_pressure_Pa(
        S.common.solar_wind,std::max(0.0,n_up_sh));
    upstream.velocity_m_s = {{Vsw,0.0,0.0}};
    upstream.magnetic_T = {{Br_sh,Bphi_sh,0.0}};

    // A compression floor is intentionally NOT passed to the shock solver.
    // The CME front is a physical fast shock only when its normal relative
    // speed exceeds the upstream fast-mode speed.  This removes the legacy
    // behavior in which sheath_comp_floor could manufacture rc>1 even for
    // V_sh=V_sw or another sub-fast disturbance.
    const swcme::shock::JumpResult jump =
        swcme::shock::solve_ideal_mhd_fast_shock(
            upstream,{{1.0,0.0,0.0}},S.V_sh_ms,P.gamma_ad);
    S.shock_jump=jump;
    S.has_shock = jump.has_shock;
    S.shock_solver_converged = jump.solver_converged;
    // A super-fast state for which the RH solve fails is a numerical failure,
    // not a no-shock state.  Do not fall back to compression=1/ambient flow.
    if (S.has_shock && !S.shock_solver_converged) {
      throw std::runtime_error(
          swcme::ModelStatus::make(
              swcme::StatusCode::ShockSolverFailure,
              "swcme1d::prepare_step Rankine-Hugoniot solve")
              .summary());
    }
    S.rc = S.has_shock ? jump.compression : 1.0;

    // Cache exact boundary values for the sheath.  In the no-shock case the
    // downstream state is intentionally identical to upstream, so no jump is
    // fabricated and the leading-edge target also remains ambient.
    S.n_up_shock = n_up_sh;
    S.n_up_le = density_upstream(S,S.r_le_m);
    S.V2_shock_ms = (S.has_shock && S.shock_solver_converged)
                        ? jump.downstream.velocity_m_s[0] : Vsw;

    if (!S.has_shock || !S.shock_solver_converged) {
      S.V_LE_ms = Vsw;
    } else {
      // The shared region helper constrains the phenomenological LE target to
      // lie between ambient flow and the exact RH downstream radial velocity.
      // Configuration validation already requires the requested factor >=1.
      S.V_LE_ms=swcme::regions::leading_edge_speed(
          Vsw,S.V2_shock_ms,P.V_sheath_LE_factor);
    }

    // Seal the complete record only after every cached value and compatibility
    // mirror has reached its final value.  Copies/moves preserve this private
    // tag automatically; callers cannot update it after modifying a mirror.
    S.integrity_digest_=prepared_state_integrity(S);

    // A successfully returned state freezes this model's configuration.
    // Locking only here leaves a model reusable after validation or numerical
    // preparation throws; no partial/failed state changes its lifecycle.
    configuration_locked_.store(true,std::memory_order_release);
    return S;
  }

  // Return the canonical 1-D shock-acceleration record.  The radial +X basis
  // used by the 1-D jump solver is embedded as a 3-vector so the serialized
  // record can be compared field-by-field with the equivalent 3-D +X case.
  swcme::ModelStatus shock_acceleration_state_checked(
      const StepState& S,
      swcme::acceleration::ShockAccelerationState& state) const {
    const swcme::ModelStatus ownership=validate_prepared_state(
        S,"swcme1d::shock_acceleration_state");
    if (!ownership.ok()) return ownership;

    const bool physical=S.has_shock && S.shock_solver_converged;
    const std::array<double,3> position{{S.r_sh_m,0.0,0.0}};
    const std::array<double,3> normal{{1.0,0.0,0.0}};
    state=swcme::acceleration::make_state(
        S.acceleration_config,true,physical,S.time_s,position,normal,S.V_sh_ms,
        S.shock_jump.compression,S.shock_jump.theta_Bn_rad,
        S.shock_jump.fast_mach,S.n_up_shock,S.B_up_T);
    return swcme::ModelStatus::success();
  }

  // Source-compatible value-returning wrapper.  It preserves the historical
  // signature while ensuring that a foreign state is never converted into a
  // plausible source record; ownership failure is surfaced as an exception.
  swcme::acceleration::ShockAccelerationState shock_acceleration_state(
      const StepState& S) const {
    swcme::acceleration::ShockAccelerationState state;
    swcme::throw_if_error(shock_acceleration_state_checked(S,state));
    return state;
  }

  // Upstream Leblanc density (fast; SI). r is clipped ≥1.05 R☉ for stability
  static inline double density_upstream(const StepState& S, double r_m){
    // The Leblanc profile is evaluated by the common production core.  This
    // wrapper is retained only for source compatibility with existing 1-D
    // code/tests; it no longer carries an independent density equation.
    return swcme::solarwind::density_m3(S.common.solar_wind,r_m);
  }

  // ----------------------------- Fast evaluator ------------------------------
  // Returns n[m⁻³] and V[m/s] at each radius (no B, no ∇·V). Designed for O(1)
  // per query, suitable for per‑particle sampling.
  /**
   * @brief Fast evaluator: n(r), V(r). Suitable for per-particle queries.
   * @param S  StepState built at desired time t.
   * @param r_m  Array of radii [m].
   * @param n_m3 Output array [m^-3].
   * @param V_ms Output array [m/s].
   * @param N    Number of points.
   * @details
   *   Region logic:
   *     (i) r > R_sh  → upstream ambient (n_up,V_sw)
   *     (ii) R_LE ≤ r ≤ R_sh → sheath (strictly monotone between pinned bounds)
   *     (iii) R_TE ≤ r < R_LE → magnetic ejecta (ME)
   *     (iv) r < R_TE → ambient again (with C¹ blend across TE)
   *   Sheath construction:
   *     n(s)=exp((1−s)ln n₂ + s ln n_up(R_LE)), V(s)=smoothstep(s^p; V₂→V_LE),
   *     with safety clamps n≥n_up,V≥V_sw to prevent undershoots.
   */
  // Checked batch evaluator.  No output value is synthesized on failure: the
  // first invalid sample is reported through ModelStatus and the caller can
  // decide whether to abort the AMPS step, log the event, or recover at a
  // higher level.  The legacy void wrapper below converts the same status into
  // an exception, preserving source compatibility without preserving silent
  // fallback semantics.
  swcme::ModelStatus evaluate_radii_fast_checked(
      const StepState& S, const double* r_m, double* n_m3, double* V_ms,
      std::size_t N) const {
    // Ownership is checked before N, pointers, or samples so a foreign state
    // is always diagnosed consistently and no destination element can change.
    const swcme::ModelStatus ownership=validate_prepared_state(
        S,"swcme1d::evaluate_radii_fast");
    if (!ownership.ok()) return ownership;
    return evaluate_radii_fast_after_validation(S,r_m,n_m3,V_ms,N);
  }

private:
  // Configuration digesting is setup work.  Centralizing refreshes here keeps
  // the public setter contract auditable and guarantees that PST03 diagnostics
  // remain exact while PST08 removes repeated parameter serialization from
  // every particle/background query.
  void refresh_configuration_digest() noexcept {
    configuration_digest_=swcme1d::configuration_digest(P);
  }

  // PST08 validation lease for the 1-D batch kernel.  This helper is private,
  // so callers cannot bypass ownership/configuration/integrity checks.  A
  // public composite method validates once and may then reuse this kernel;
  // validation is consequently O(1) per API call instead of being repeated by
  // every nested evaluator.  No state or output pointer is retained here.
  swcme::ModelStatus evaluate_radii_fast_after_validation(
      const StepState& S, const double* r_m, double* n_m3, double* V_ms,
      std::size_t N) const {
    if (N==0) return swcme::ModelStatus::success();
    if (!r_m || !n_m3 || !V_ms) {
      return swcme::ModelStatus::make(
          swcme::StatusCode::NullPointer,"swcme1d::evaluate_radii_fast");
    }
    if (S.has_shock && !S.shock_solver_converged) {
      return swcme::ModelStatus::make(
          swcme::StatusCode::ShockSolverFailure,
          "swcme1d::evaluate_radii_fast shock state");
    }

    const double Vsw=S.V_up_ms;
    for (std::size_t i=0;i<N;++i) {
      const double r=r_m[i];
      if (!std::isfinite(r)) {
        return swcme::ModelStatus::make_value(
            swcme::StatusCode::NonFiniteInput,
            "swcme1d::evaluate_radii_fast radius",r,i);
      }
      if (r<swcme::solarwind::MIN_RADIUS_M) {
        return swcme::ModelStatus::make_value(
            swcme::StatusCode::OutsideModelDomain,
            "swcme1d::evaluate_radii_fast radius",r,i);
      }

      const double n_up=density_upstream(S,r);
      if (!std::isfinite(n_up) || n_up<0.0 || !std::isfinite(Vsw)) {
        return swcme::ModelStatus::make(
            swcme::StatusCode::NonFiniteResult,
            "swcme1d::evaluate_radii_fast ambient state",i);
      }

      if (S.region_config.mode==swcme::regions::Mode::ShockOnly) {
        n_m3[i]=n_up;
        V_ms[i]=Vsw;
        continue;
      }

      const swcme::regions::Boundaries& b=S.region_boundaries;
      const swcme::regions::Location loc=swcme::regions::locate(r,b);
      // Density keeps its own physically appropriate interpolation, while the
      // radial velocity is obtained from one shared value+derivative profile in
      // swcme_regions.hpp.  The same profile is used below for analytical
      // div(V), so the divergence cannot differentiate a different V(r) than
      // the transport evaluator actually returns.
      const auto sheath_density=[&](double rr) {
        const double n_local_up=density_upstream(S,rr);
        if (!S.has_shock) return n_local_up;
        const double w=swcme::regions::sheath_profile_weight(
            rr,b,S.region_config.sheath_ramp_power);
        const double n2=S.shock_jump.downstream.rho_kg_m3/MP;
        const double n_le=density_upstream(S,b.R_le_m);
        return swcme::regions::log_lerp_positive(n2,n_le,w);
      };
      const auto ejecta_density=[&](double rr) {
        return S.region_config.f_ME*density_upstream(S,rr);
      };

      double n=n_up;
      if (loc.region==swcme::regions::Region::ShockTransition) {
        if (S.has_shock) {
          const double n2=S.shock_jump.downstream.rho_kg_m3/MP;
          n=swcme::regions::lerp(n_up,n2,loc.blend);
        }
      } else if (loc.region==swcme::regions::Region::Sheath) {
        n=sheath_density(r);
      } else if (loc.region==swcme::regions::Region::LeadingTransition) {
        n=swcme::regions::lerp(
            sheath_density(r),ejecta_density(r),loc.blend);
      } else if (loc.region==swcme::regions::Region::Ejecta) {
        n=ejecta_density(r);
      } else if (loc.region==swcme::regions::Region::TrailingTransition) {
        n=swcme::regions::lerp(ejecta_density(r),n_up,loc.blend);
      }

      const swcme::regions::RadialVelocityState velocity=
          swcme::regions::radial_velocity_state(
              r,b,S.has_shock,Vsw,S.V2_shock_ms,S.V_LE_ms,
              S.region_config.V_ME_factor*Vsw,
              S.region_config.sheath_ramp_power);
      const double V=velocity.velocity_m_s;

      if (!std::isfinite(n) || n<0.0 || !std::isfinite(V)) {
        return swcme::ModelStatus::make(
            swcme::StatusCode::NonFiniteResult,
            "swcme1d::evaluate_radii_fast regional state",i);
      }
      n_m3[i]=n;
      V_ms[i]=V;
    }
    return swcme::ModelStatus::success();
  }

public:

  void evaluate_radii_fast(const StepState& S, const double* r_m,
                           double* n_m3, double* V_ms,
                           std::size_t N) const {
    swcme::throw_if_error(evaluate_radii_fast_checked(S,r_m,n_m3,V_ms,N));
  }

  // -------------------------- Full evaluator (with B, ∇·V) -------------------
  /**
   * @brief Full evaluator: n,V plus Parker B components and ∇·V.
   * @param dr_frac Retained for source compatibility; 1-D divergence is now
   *        analytical and does not use a finite-difference step.
   * @details The Parker field is computed from Br(1 AU) and k.  For the purely
   * radial 1-D velocity field, div(V)=2Vr/r+dVr/dr is evaluated analytically
   * from the exact same region blend that produces V(r).
   */
  swcme::ModelStatus evaluate_radii_with_B_div_checked(
      const StepState& S, const double* r_m,
      double* n_m3, double* V_ms,
      double* Br_T, double* Bphi_T, double* Bmag_T,
      double* divV, std::size_t N, double dr_frac=1e-3) const {
    // Check ownership at this public boundary before even the N==0 shortcut.
    // This gives direct full-field callers the same deterministic PST02 status
    // and untouched-output guarantee as the fast evaluator.
    const swcme::ModelStatus ownership=validate_prepared_state(
        S,"swcme1d::evaluate_radii_with_B_div");
    if (!ownership.ok()) return ownership;
    if (N==0) return swcme::ModelStatus::success();
    if (!r_m || !n_m3 || !V_ms) {
      return swcme::ModelStatus::make(
          swcme::StatusCode::NullPointer,"swcme1d::evaluate_radii_with_B_div");
    }
    (void)dr_frac; // analytical 1-D divergence; kept in the API for compatibility

    // The outer guard above has already authenticated the complete state.
    // Calling the private kernel avoids a second full digest while preserving
    // the exact argument, domain, numerical-status, and output behavior.
    swcme::ModelStatus status=
        evaluate_radii_fast_after_validation(S,r_m,n_m3,V_ms,N);
    if (!status.ok()) return status;

    for (std::size_t i=0;i<N;++i) {
      const double r=r_m[i];  // domain/finite checks already performed above
      const swcme::solarwind::ParkerComponents parker=
          swcme::solarwind::parker_components(
              S.common.solar_wind,r,P.sin_theta);
      double Br=parker.Br_T;
      double Bph=parker.Bphi_T;

      if (S.region_config.mode==swcme::regions::Mode::FullICME) {
        const swcme::regions::Location loc=
            swcme::regions::locate(r,S.region_boundaries);
        if (loc.region==swcme::regions::Region::ShockTransition && S.has_shock) {
          Br=swcme::regions::lerp(
              Br,S.shock_jump.downstream.magnetic_T[0],loc.blend);
          Bph=swcme::regions::lerp(
              Bph,S.shock_jump.downstream.magnetic_T[1],loc.blend);
        } else if (loc.region==swcme::regions::Region::Sheath ||
                   loc.region==swcme::regions::Region::LeadingTransition) {
          const double w=swcme::regions::sheath_profile_weight(
              r,S.region_boundaries,S.region_config.sheath_ramp_power);
          const swcme::solarwind::ParkerComponents parker_le=
              swcme::solarwind::parker_components(
                  S.common.solar_wind,S.region_boundaries.R_le_m,P.sin_theta);
          double Br_sheath=Br,Bph_sheath=Bph;
          if (S.has_shock) {
            Br_sheath=swcme::regions::lerp(
                S.shock_jump.downstream.magnetic_T[0],parker_le.Br_T,w);
            Bph_sheath=swcme::regions::lerp(
                S.shock_jump.downstream.magnetic_T[1],parker_le.Bphi_T,w);
          }
          if (loc.region==swcme::regions::Region::LeadingTransition) {
            Br=swcme::regions::lerp(Br_sheath,Br,loc.blend);
            Bph=swcme::regions::lerp(Bph_sheath,Bph,loc.blend);
          } else {
            Br=Br_sheath;
            Bph=Bph_sheath;
          }
        }
      }

      const double Bmag=std::hypot(Br,Bph);
      if (!std::isfinite(Br) || !std::isfinite(Bph) || !std::isfinite(Bmag)) {
        return swcme::ModelStatus::make(
            swcme::StatusCode::NonFiniteResult,
            "swcme1d::evaluate_radii_with_B_div magnetic field",i);
      }
      if (Br_T) Br_T[i]=Br;
      if (Bphi_T) Bphi_T[i]=Bph;
      if (Bmag_T) Bmag_T[i]=Bmag;

      if (divV) {
        // The 1-D field is exactly radial, so the spherical identity is exact:
        //     div(V) = 2 Vr/r + dVr/dr.
        // SHOCK_ONLY has dVr/dr=0.  FULL_ICME obtains dVr/dr analytically from
        // the same common smoothstep profile used by evaluate_radii_fast().
        double dV_dr=0.0;
        if (S.region_config.mode==swcme::regions::Mode::FullICME) {
          const swcme::regions::RadialVelocityState velocity=
              swcme::regions::radial_velocity_state(
                  r,S.region_boundaries,S.has_shock,S.V_up_ms,S.V2_shock_ms,
                  S.V_LE_ms,S.region_config.V_ME_factor*S.V_up_ms,
                  S.region_config.sheath_ramp_power);
          dV_dr=velocity.d_velocity_dr_s_inv;
          // The fast evaluator was already called above; this equality check is
          // intentionally numerical rather than an assertion so a future drift
          // is surfaced through the public status contract.
          if (!std::isfinite(velocity.velocity_m_s) ||
              std::abs(velocity.velocity_m_s-V_ms[i]) >
                  64.0*std::numeric_limits<double>::epsilon()*
                  std::max({1.0,std::abs(velocity.velocity_m_s),std::abs(V_ms[i])})) {
            return swcme::ModelStatus::make(
                swcme::StatusCode::NonFiniteResult,
                "swcme1d::evaluate_radii_with_B_div velocity/profile mismatch",i);
          }
        }
        const swcme::divergence::RadialTerms terms=
            swcme::divergence::radial_terms(r,V_ms[i],dV_dr);
        if (!std::isfinite(terms.total_s_inv)) {
          return swcme::ModelStatus::make(
              swcme::StatusCode::NonFiniteResult,
              "swcme1d::evaluate_radii_with_B_div divergence",i);
        }
        divV[i]=terms.total_s_inv;
      }
    }
    return swcme::ModelStatus::success();
  }

  void evaluate_radii_with_B_div(const StepState& S,
                                 const double* r_m,
                                 double* n_m3, double* V_ms,
                                 double* Br_T, double* Bphi_T, double* Bmag_T,
                                 double* divV, std::size_t N,
                                 double dr_frac=1e-3) const {
    swcme::throw_if_error(evaluate_radii_with_B_div_checked(
        S,r_m,n_m3,V_ms,Br_T,Bphi_T,Bmag_T,divV,N,dr_frac));
  }

  // ------------------------------- Tecplot writer ----------------------------
  // Writes: r[m], n[m⁻³], V[m/s], Br[T], Bphi[T], Bmag[T], divV[s⁻¹], rc, R_sh, R_LE, R_TE
  /**
   * @brief Write a Tecplot POINT zone with 13 columns: r[m], R[AU], rSun[R_s],
   *        n[m^-3], V[m/s], Br[T], Bphi[T], Bmag[T], divV[s^-1], rc, R_sh[m],
   *        R_LE[m], R_TE[m].
   * @note The VARIABLES list and row format specifier are kept in sync (13/13).
   */
  bool write_tecplot_radial_profile(const StepState& S,
                                    const double* r_m,
                                    const double* n_m3,const double* V_ms,
                                    const double* Br_T,const double* Bphi_T,const double* Bmag_T,
                                    const double* divV,std::size_t N,
                                    const char* path,double time_simulation=-1.0) const {
    // The legacy boolean writer cannot return ModelStatus.  Throw before
    // opening so a foreign state cannot truncate or partially replace an
    // existing file; checked callers should use the method below instead.
    swcme::throw_if_error(validate_prepared_state(
        S,"swcme1d::write_tecplot_radial_profile"));
    return write_tecplot_radial_profile_checked(
        S,r_m,n_m3,V_ms,Br_T,Bphi_T,Bmag_T,divV,N,path,
        time_simulation).ok();
  }

  // Status-returning output entry point used by AMPS and validation code.  The
  // ownership gate precedes all argument inspection and file-system access,
  // which guarantees that PST02 rejection leaves a pre-existing destination
  // byte-for-byte unchanged.  OUT05 then validates every requested row before
  // OUT03 creates a sibling staging file, and the transaction is committed only
  // after the complete OUT02 lifecycle succeeds.
  swcme::ModelStatus write_tecplot_radial_profile_checked(
      const StepState& S, const double* r_m, const double* n_m3,
      const double* V_ms, const double* Br_T, const double* Bphi_T,
      const double* Bmag_T, const double* divV, std::size_t N,
      const char* path, double time_simulation=-1.0,
      const swcme::output::FileOperations* file_operations=nullptr) const {
    const swcme::ModelStatus ownership=validate_prepared_state(
        S,"swcme1d::write_tecplot_radial_profile");
    if (!ownership.ok()) return ownership;
    if (!r_m || !n_m3 || !V_ms || !Br_T || !Bphi_T || !Bmag_T || !path)
      return swcme::ModelStatus::make(
          swcme::StatusCode::NullPointer,
          "swcme1d::write_tecplot_radial_profile arguments");
    if (N==0)
      return swcme::ModelStatus::make(
          swcme::StatusCode::InvalidConfiguration,
          "swcme1d::write_tecplot_radial_profile sample count");
    if (!std::isfinite(time_simulation))
      return swcme::ModelStatus::make_value(
          swcme::StatusCode::NonFiniteInput,
          "swcme1d::write_tecplot_radial_profile time",time_simulation);

    // OUT05 is deliberately a complete, read-only pass over the caller's
    // dataset.  In particular, this pass occurs before FileOperations is even
    // selected, so a bad late row cannot create a staging file or invoke an
    // injected open callback.  Radius is an input and retains the evaluator's
    // NONFINITE_INPUT / OUTSIDE_MODEL_DOMAIN distinction; precomputed model
    // fields are results and therefore use NONFINITE_RESULT.  The first bad
    // row is retained in sample_index for reproducible campaign diagnostics.
    for (std::size_t i=0; i<N; ++i) {
      if (!std::isfinite(r_m[i]))
        return swcme::ModelStatus::make_value(
            swcme::StatusCode::NonFiniteInput,
            "swcme1d::write_tecplot_radial_profile radius",r_m[i],i);
      if (r_m[i]<swcme::solarwind::MIN_RADIUS_M)
        return swcme::ModelStatus::make_value(
            swcme::StatusCode::OutsideModelDomain,
            "swcme1d::write_tecplot_radial_profile radius",r_m[i],i);
      const double values[]={n_m3[i],V_ms[i],Br_T[i],Bphi_T[i],Bmag_T[i],
                             divV ? divV[i] : 0.0};
      for (double value : values) {
        if (!std::isfinite(value))
          return swcme::ModelStatus::make_value(
              swcme::StatusCode::NonFiniteResult,
              "swcme1d::write_tecplot_radial_profile field",value,i);
      }
    }

    const swcme::output::FileOperations& operations=file_operations
        ? *file_operations : swcme::output::stdio_file_operations();
    swcme::output::CheckedTextFile output(operations);
    if (!output.open_transactional(path))
      return swcme::ModelStatus::make(
          swcme::StatusCode::FileOpenFailure,
          "swcme1d::write_tecplot_radial_profile open");

    // Every call carries a stable phase/row context.  CheckedTextFile stops
    // issuing writes after the first failure; finish() closes and discards the
    // staging file while preserving that original diagnostic.
    output.print("swcme1d radial profile title",swcme::ModelStatus::npos,
                 "TITLE=\"1D SW+CME radial profile\"\n");
    output.print("swcme1d radial profile variables",swcme::ModelStatus::npos,
                 "VARIABLES=\"r[m]\",\"R[AU]\",\"rSun[R_s]\","
                 "\"n[m^-3]\",\"V[m/s]\",\"Br[T]\",\"Bphi[T]\","
                 "\"Bmag[T]\",\"divV[s^-1]\",\"rc\","
                 "\"R_sh[m]\",\"R_LE[m]\",\"R_TE[m]\"\n");
    output.print("swcme1d radial profile zone",swcme::ModelStatus::npos,
                 "ZONE T=\"radial\", I=%zu, F=POINT\n",N);

    for (std::size_t i=0; i<N && output.good(); ++i) {
      const double r=r_m[i];
      const double R_AU=swcme::units::m_to_au(r);
      const double Rsun=swcme::units::m_to_solar_radii(r);
      const double dv=divV ? divV[i] : 0.0;
      output.print(
          "swcme1d radial profile row",i,
          "% .9e % .9e % .9e % .9e % .9e % .9e % .9e % .9e "
          "% .9e % .9e % .9e % .9e % .9e\n",
          r,R_AU,Rsun,n_m3[i],V_ms[i],Br_T[i],Bphi_T[i],Bmag_T[i],dv,
          S.rc,S.r_sh_m,S.r_le_m,S.r_te_m);
    }

    if (time_simulation>=0.0 && output.good())
      output.print("swcme1d radial profile time",swcme::ModelStatus::npos,
                   "# t = %.3f s\n",time_simulation);

    return output.finish(
        "swcme1d radial profile flush",
        "swcme1d radial profile stream error",
        "swcme1d radial profile close",
        "swcme1d radial profile commit");
  }

  // Convenience wrapper from radii only
  /**
   * @brief Convenience wrapper: given radii only, compute fields and write
   *        the Tecplot file.
   */
  bool write_tecplot_radial_profile_from_r(const StepState& S,
                                           const double* r_m,std::size_t N,
                                           const char* path,double time_simulation=-1.0) const {
    // Validate before allocating temporary arrays.  Besides preserving the
    // output file, this prevents a foreign-state exception from bypassing the
    // legacy manual cleanup path below.
    swcme::throw_if_error(validate_prepared_state(
        S,"swcme1d::write_tecplot_radial_profile_from_r"));
    if (!r_m || N==0) return false;
    double *n=new double[N], *V=new double[N], *Br=new double[N], *Bph=new double[N], *Bm=new double[N], *dv=new double[N];
    evaluate_radii_with_B_div(S, r_m, n, V, Br, Bph, Bm, dv, N);
    const bool ok = write_tecplot_radial_profile(S, r_m, n, V, Br, Bph, Bm, dv, N, path, time_simulation);
    delete [] n; delete [] V; delete [] Br; delete [] Bph; delete [] Bm; delete [] dv; return ok;
  }

  /**
 * @brief Write a Tecplot POINT zone with shock kinematics vs time.
 *
 * Columns:
 *   1) t[s]          — simulation time since launch [seconds]
 *   2) R_sh[R_s]     — shock apex heliocentric distance in solar radii
 *   3) V_sh[km/s]    — shock apex speed in km/s
 *   4) rc            — compression ratio proxy used in the model
 *
 * Notes:
 *   • The series starts at t = 0 and ends exactly at t = t_end_s.
 *   • If the model’s gating determines there is no forward shock yet,
 *     prepare_step() will set rc ≈ 1 (and your StepState.has_shock would be false
 *     if you kept that flag). We still report R_sh and V_sh from DBM.
 *
 * @param t_end_s  End time [s] (must be ≥ 0).
 * @param N        Number of samples (rows). If N==1, the single row is at t=t_end_s.
 * @param path     Output path for the Tecplot .dat file.
 * @return true on success; false for invalid arguments, open failure, or any
 *         write/flush/stream/close failure.
 */
bool write_tecplot_shock_vs_time(double t_end_s, std::size_t N,
                                 const char* path) const {
  return write_tecplot_shock_vs_time_checked(t_end_s,N,path).ok();
}

// Status-returning companion to the legacy bool API.  OUT02 requires output
// failures to remain distinguishable from malformed input and from kinematic
// preparation failures; OUT05 prepares and validates the complete time series
// before any file is opened; and OUT03 guarantees that only a complete history
// replaces the destination.  Archival callers should use this form.
swcme::ModelStatus write_tecplot_shock_vs_time_checked(
    double t_end_s, std::size_t N, const char* path,
    const swcme::output::FileOperations* file_operations=nullptr) const {
  if (!path)
    return swcme::ModelStatus::make(
        swcme::StatusCode::NullPointer,"swcme1d shock history path");
  if (N==0 || !std::isfinite(t_end_s) || t_end_s<0.0)
    return swcme::ModelStatus::make_value(
        swcme::StatusCode::InvalidConfiguration,
        "swcme1d shock history range",t_end_s);

  // Prepare every requested time first and retain the immutable records for
  // the subsequent formatting pass.  This avoids both time-of-check/time-of-
  // use drift and the former behavior where an out-of-range late data-driven
  // sample was discovered only after a staging file and partial header existed.
  // Computing t from the normalized sample fraction also makes the final time
  // exactly t_end_s instead of accumulating repeated floating-point additions.
  std::vector<StepState> states;
  states.reserve(N);
  try {
    for (std::size_t i=0; i<N; ++i) {
      const double t=(N>1)
          ? t_end_s*(static_cast<double>(i)/static_cast<double>(N-1))
          : t_end_s;
      const StepState prepared=prepare_step(t);
      const double values[]={prepared.r_sh_m,prepared.V_sh_ms,prepared.rc};
      for (double value : values) {
        if (!std::isfinite(value))
          return swcme::ModelStatus::make_value(
              swcme::StatusCode::NonFiniteResult,
              "swcme1d shock history prepared field",value,i);
      }
      if (prepared.r_sh_m<swcme::solarwind::MIN_RADIUS_M)
        return swcme::ModelStatus::make_value(
            swcme::StatusCode::OutsideModelDomain,
            "swcme1d shock history radius",prepared.r_sh_m,i);
      states.push_back(prepared);
    }
  } catch (...) {
    // No transaction exists at this point.  Preserve the first unprepared
    // sample index so a data-driven campaign can identify the exact requested
    // time that exceeded its validated interpolation/extrapolation domain.
    const std::size_t failed_index=states.size();
    const double failed_time=(N>1)
        ? t_end_s*(static_cast<double>(failed_index)/
                   static_cast<double>(N-1)) : t_end_s;
    return swcme::ModelStatus::make_value(
        swcme::StatusCode::NonFiniteResult,
        "swcme1d shock history preparation",failed_time,failed_index);
  }

  const swcme::output::FileOperations& operations=file_operations
      ? *file_operations : swcme::output::stdio_file_operations();
  swcme::output::CheckedTextFile output(operations);
  if (!output.open_transactional(path))
    return swcme::ModelStatus::make(
        swcme::StatusCode::FileOpenFailure,"swcme1d shock history open");

  output.print("swcme1d shock history header",swcme::ModelStatus::npos,
      "TITLE=\"Shock kinematics vs time\"\n"
      "VARIABLES=\"t[s]\",\"R_sh[R_s]\",\"V_sh[km/s]\",\"rc\"\n");
  output.print("swcme1d shock history zone",swcme::ModelStatus::npos,
               "ZONE T=\"shock_vs_time\", N=%zu, F=POINT\n",N);

  for (std::size_t i=0; i<N && output.good(); ++i) {
    const double t=(N>1)
        ? t_end_s*(static_cast<double>(i)/static_cast<double>(N-1))
        : t_end_s;
    const StepState& S=states[i];
    const double Rsh_Rs=S.r_sh_m/Rs;
    const double Vsh_kms=swcme::units::m_per_s_to_km_per_s(S.V_sh_ms);
    output.print("swcme1d shock history row",i,
                 "% .9e % .9e % .9e % .9e\n",
                 t,Rsh_Rs,Vsh_kms,S.rc);
  }

  return output.finish(
      "swcme1d shock history flush",
      "swcme1d shock history stream error",
      "swcme1d shock history close",
      "swcme1d shock history commit");
}


private:
  // Legacy fluent setters have no status return, so post-prepare mutation is
  // rejected with one consistent exception before any Params field changes.
  // Configuration and preparation are an exclusive setup phase; once frozen,
  // the atomic flag supports concurrent read-only state evaluation safely.
  void require_configuration_mutable(const char* operation) const {
    if (configuration_locked_.load(std::memory_order_acquire)) {
      throw std::logic_error(std::string("swcme1d::")+operation+
          ": model configuration is immutable after successful prepare_step(); "
          "construct model.reconfigured(params) instead");
    }
  }

  Params P;
  // Runtime owner token stamped into every StepState returned by this model.
  // It is intentionally not derived from Params; identical Model instances
  // must remain distinct owners for PST02.
  swcme::ModelIdentity model_identity_;
  // Cached complete Params fingerprint.  It changes only through guarded setup
  // operations and is immutable once a StepState has been issued.
  swcme::ConfigurationDigest configuration_digest_;
  // False during the legacy setup phase and permanently true after the first
  // successful prepare_step().  It is never reset on a live model.
  mutable std::atomic<bool> configuration_locked_;
};

// -----------------------------------------------------------------------------
// Extended usage examples (copy‑paste friendly)
// -----------------------------------------------------------------------------
/*
Example 1 — basic profile at a single time
-----------------------------------------
  using namespace swcme1d;
  Model sw;
  sw.SetAmbient(450, 5.5, 4.0, 1.5e5)
    .SetCME(20.0, 1800, 7.5e-8)
    .SetGeometry(0.12, 0.25)
    .SetSmoothing(0.008, 0.02, 0.03)
    .SetSheathEjecta(1.2, 2.0, 1.10, 0.5, 0.8);

  auto S = sw.prepare_step(24*3600.0);
  double rA[5] = {0.5*AU, 0.8*AU, 1.0*AU, 1.2*AU, 1.5*AU};
  double n[5], V[5], Br[5], Bp[5], Bm[5], dV[5];
  sw.evaluate_radii_with_B_div(S, rA, n, V, Br, Bp, Bm, dV, 5);

  // Optional: write Tecplot file
  sw.write_tecplot_radial_profile(S, rA, n, V, Br, Bp, Bm, dV, 5, "swcme_profile.dat", 24*3600.0);

Example 2 — per‑particle usage (fast evaluator)
-----------------------------------------------
  auto S = sw.prepare_step(t_now);
  for (size_t p=0; p<NPART; ++p){
    double r = particle_r[p];
    double n,V; sw.evaluate_radii_fast(S, &r, &n, &V, 1);
    // use n,V for adiabatic losses, scattering rates, etc.
  }

  Example 3 - output shock properties 

  swcme1d::Model m;
  // (configure m’s Params as desired...)
  const double t_end = 48.0 * 3600.0; // 48 hours
  const std::size_t N = 241;          // e.g., every 12 minutes
  m.write_tecplot_shock_vs_time(t_end, N, "shock_vs_time.dat");

*/

} // namespace swcme1d

#endif // SWCME1D_HPP
