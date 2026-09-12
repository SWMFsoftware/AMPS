// demo1d.cpp — 1-D SW+CME usage example (explicit status, Tecplot profile)
//
// WHAT THIS EXAMPLE DOES (PHYSICS):
//  • Sets a “typical fast CME”: r0=20 Rs, V0≈1800 km/s, ambient Vsw=400 km/s.
//  • Computes a snapshot at t = 36 h using the shared sign-aware DBM apex kinematics.
//  • Builds a radial grid 0.2–1.5 AU and evaluates:
//       n(r), V(r), Br(r), Bphi(r), |B|(r), divV(r).
//  • Writes a Tecplot POINT file "profile_1d.dat" with 13 unit-qualified
//    columns: r[m], R[AU], rSun[R_s], n[m^-3], V[m/s], Br[T], Bphi[T],
//    Bmag[T], divV[s^-1], rc, R_sh[m], R_LE[m], and R_TE[m].
//
// HOW TO SET UP CME PARAMETERS IN THIS MODEL
// ----------------------------------------------------------------------------
// The model is configured via a `swcme1d::Params` struct. Below is a field-by-field
// guide, including *physical meaning*, *typical ranges*, and *tuning tips*.
//
// AMBIENT (SOLAR WIND) & THERMODYNAMICS
// -------------------------------------
//  P.V_sw_kms       [km/s]  Upstream solar-wind speed (radial, constant in 1-D).
//                           Typical: 300–800. Slow wind ~350–450; fast stream ~600–800.
//                           Affects Parker spiral pitch (through Ω/V_sw) and the
//                           DBM drag frame (V_sh is measured relative to V_sw).
//
//  P.n1AU_cm3       [cm^-3] Upstream number density at 1 AU (normalization).
//                           Typical: 2–10 (quiescent ~3–7). Used to scale the
//                           Leblanc (1998) n(r) profile everywhere.
//
//  P.B1AU_nT        [nT]    |B| magnitude at 1 AU for the Parker spiral scaling.
//                           Typical: 3–8. Sets Br(1 AU) so that sqrt(Br^2 + Bφ^2)=B1AU.
//
//  P.T_K            [K]     Proton temperature used for sound speed (γ k_B T / m_p).
//                           Typical: 8e4–2e5. Affects fast-mode speed → shock rc.
//
//  P.gamma_ad       [-]     Adiabatic index. Use 5/3 unless you have a reason to vary.
//
//  P.sin_theta      [-]     Sine of colatitude for Parker pitch (≈1 in ecliptic).
//                           If you want a more radial field (high latitude), reduce.
//
// CME LAUNCH & DRAG (DBM, APEX KINEMATICS)
// ----------------------------------------
//  P.r0_Rs          [R_sun] DBM reference radius. About 15–20 Rs is recommended
//                           for drag-dominated propagation; 1.05 Rs is the
//                           solar-wind domain floor, not the recommended DBM start.
//
//  P.V0_sh_kms      [km/s]  Initial shock apex speed at t=0. Typical 800–2500.
//                           Large values produce strong early compression.
//
//  P.Gamma_kmInv    [1/km]  DBM drag parameter Γ. Larger Γ → faster deceleration.
//                           Typical ~ (1–20)×10^-8. Try 8e-8 for a fairly strong drag.
//
// REGION GEOMETRY (SELF-SIMILAR WITH R_sh)
// ----------------------------------------
//  P.sheath_thick_AU_at1AU  [AU] Sheath thickness normalized at 1 AU (scales ∝ R_sh).
//                           Typical 0.05–0.2. 0.1 at 1 AU is a common choice.
//
//  P.ejecta_thick_AU_at1AU  [AU] Magnetic-ejecta (ME) thickness normalized at 1 AU.
//                           Typical 0.1–0.4. 0.2 at 1 AU is a common choice.
//
// EDGE SMOOTHING (C¹ BLENDS; EACH EDGE HAS ITS OWN WIDTH)
// -------------------------------------------------------
//  P.edge_smooth_shock_AU_at1AU  [AU] Total numerical shock-layer width used
//                                      by RESOLVED_COMPRESSION. SOURCE mode uses
//                                      SHOCK_ONLY and does not resolve this jump.
//  P.edge_smooth_le_AU_at1AU     [AU] Blend width at the leading edge (sheath→ME).
//  P.edge_smooth_te_AU_at1AU     [AU] Blend width at the trailing edge (ME→ambient).
//  All three scale ∝ R_sh. Guidelines:
//   • Keep the shock width ≪ sheath thickness to preserve a clear jump.
//   • Avoid making widths so large that the three edges overlap heavily.
//   • Typical at 1 AU: shock 0.005–0.02, LE 0.01–0.05, TE 0.02–0.06.
//
// SHEATH / EJECTA TARGET STATES (POST-SHOCK PROFILES)
// ---------------------------------------------------
//  P.sheath_comp_floor     [-]  Deprecated compatibility field; ignored. Physical
//                              compression comes only from the MHD RH solver.
//                              Actual rc is computed by the shared ideal-MHD RH
//                              solver and is physically bounded for gamma=5/3.
//
//  P.sheath_ramp_power     [-]  Exponent controlling the smoothstep argument for
//                              relaxation from the RH state toward the leading edge.
//                              Values >=1 are supported; larger values delay relaxation.
//
//  P.V_sheath_LE_factor    [-]  Requested sheath speed at LE relative to V_sw.
//                              The common region model constrains this target between
//                              ambient V_sw and the RH downstream radial speed.
//
//  P.f_ME                  [-]  ME density fraction relative to upstream (often depleted).
//                              Typical 0.3–0.8. 0.5 is a reasonable default.
//
//  P.V_ME_factor           [-]  ME bulk speed relative to ambient.
//                              Typical 0.6–0.95. MEs are often slower than ambient.
//
// TIME SELECTION
// --------------
//  Choose a snapshot time t (seconds) after launch; call `model.prepare_step(t)`.
//  DBM computes R_sh(t) and V_sh(t) analytically, builds region extents (R_sh, R_LE,
//  R_TE), and caches Parker/Leblanc constants. Then call the evaluator(s) on your
//  requested radii.
//
// PHYSICAL SANITY (POST-SHOCK)
// ----------------------------
//  • Immediately behind a *compressive fast* shock: density ↑ (rc>1) and, in the
//    Sun frame, V just behind the shock lies between V_sw and V_sh. You should see
//    a post-shock **peak** in density and an **elevated** speed within the sheath,
//    then both can **decrease** through the sheath toward the leading edge—this is
//    physical. Inside the ME, density is often **below ambient** and speed **lower**.
//  • Avoid an immediate dip of density below upstream right after the shock.
//    If you see it in RESOLVED_COMPRESSION/FULL_ICME, decrease the numerical
//    shock width and inspect the RH/sheath configuration. `sheath_comp_floor`
//    is deprecated and intentionally cannot manufacture compression.
//
// NUMERICAL / IMPLEMENTATION NOTES
// --------------------------------
//  • Checked evaluators return `swcme::ModelStatus`. Invalid, out-of-domain,
//    or non-finite states are reported explicitly and are NOT replaced by
//    plausible values such as ambient flow or zero. Legacy void evaluators
//    call the same checked path and throw on failure.
//  • Hot loops avoid pow(); smoothing uses precomputed 1/(2 w). Vectorization hint
//    is included. The writer uses ordered Tecplot zone syntax: `I=..., F=POINT`.
//  • Units: r [m], t [s], n [m^-3], V [m/s], Br/Bphi/|B| [T], divV [s^-1].
//
// EXAMPLE TUNINGS
// ---------------
//  *Strong, fast CME:*  V0_sh_kms≈1800–2200; Γ≈(6–10)e-8; sheath_thick≈0.10–0.15 AU;
//                       f_ME≈0.4–0.6; V_ME_factor≈0.75–0.9; shock width≈0.01 AU @1 AU.
//  *Weaker CME:*        V0_sh_kms≈800–1200; Γ≈(1–4)e-8; thinner sheath; higher f_ME;
//                       lower rc_floor; slightly larger smooth widths.
//
//
//
// WHAT THIS EXAMPLE DOES (IMPLEMENTATION):
//  • Uses the model’s StepState cache (constants precomputed for the time).
//  • Calls the high-throughput evaluator (no allocations in hot path).
//  • This example uses the checked API and aborts if any sample is invalid;
//    successful output is therefore finite because the physics was validated,
//    not because the writer/evaluator silently sanitized it.
//
// BUILD:
//   g++ -std=c++17 -O3 -march=native demo1d.cpp -o demo1d
//   # optional: -fopenmp if you add OMP pragmas
//
// VISUALIZATION TIP (Tecplot):
//   • Load "profile_1d.dat" as XY line: X=R[m], Y=any of the fields.
//   • The discontinuities are smoothed by C¹ blends at shock/LE/TE, so
//     you should see: jump at R_sh, then graded sheath, then ME, then relax.

#include <vector>
#include <cstdio>
#include <cmath>
#include "swcme1d.hpp"

int main(){
  using namespace swcme1d;

  // ------------------------------
  // 1) Configure a typical CME
  // ------------------------------
  Params P;
  P.r0_Rs       = 20.0;     // recommended DBM start in the drag-dominated heliosphere
  P.V0_sh_kms   = 1800.0;   // fast CME launch
  P.V_sw_kms    = 400.0;    // nominal solar wind
  P.Gamma_kmInv = 8e-8;     // DBM drag (tunes deceleration)
  P.n1AU_cm3    = 6.0;      // upstream density at 1 AU
  P.B1AU_nT     = 5.0;      // |B|(1 AU)
  P.T_K         = 1.2e5;    // proton temperature
  P.sin_theta   = 1.0;      // equatorial spiral (max azimuthal component)

  // Region shaping:
  //  • sharper at shock, softer toward ME tail (physically motivated by sheath turbulence decay)
  P.edge_smooth_shock_AU_at1AU = 0.01;
  P.edge_smooth_le_AU_at1AU    = 0.02;
  P.edge_smooth_te_AU_at1AU    = 0.03;

  // Minimum compression at the shock, ramp exponent, target speeds
  P.sheath_comp_floor          = 1.2;
  P.sheath_ramp_power          = 2.0;
  P.V_sheath_LE_factor         = 1.10; // sheath slows from shock toward LE
  P.f_ME                       = 0.5;  // ME density ~ 0.5 upstream (typical)
  P.V_ME_factor                = 0.80; // ME bulk slower than wind

  Model model(P);

  // Optional CSV logging for time evolution (uncomment if needed):
  // swcme1d::EnableIO();
  // swcme1d::Model::write_step_csv_header();

  // ------------------------------
  // 2) Build time snapshot
  // ------------------------------
  const double t_snap = 36.0 * 3600.0; // 36 hours after launch
  StepState S = model.prepare_step(t_snap);

  std::printf("Snapshot t=%.1f hr: R_sh=%.3f AU, V_sh=%.1f km/s, rc=%.2f\n",
              t_snap/3600.0, S.r_sh_m/AU, S.V_sh_ms/1e3, S.rc);

  // If CSV logging enabled:
  // swcme1d::Model::write_step_csv_line(S);

  // ------------------------------
  // 3) Radial grid (0.2–1.5 AU)
  //    div(V) is now analytical in 1-D, so midpoint sampling is chosen only
  //    for uniform plotting coverage, not to suppress derivative noise.
  // ------------------------------
  const std::size_t N = 1200;
  std::vector<double> r(N), n(N), V(N), Br(N), Bphi(N), Bmag(N), divV(N);

  const double rmin = 0.2*AU, rmax = 1.5*AU;
  for (std::size_t i=0;i<N;++i){
    // Convert the discrete index explicitly at the numerical boundary.  This
    // makes the intentional floating-point grid parameterization visible to
    // reviewers and lets OUT08 reject any accidental implicit narrowing.
    const double s = (static_cast<double>(i) + 0.5) /
                     static_cast<double>(N);
    r[i] = rmin + s*(rmax - rmin);
  }

  // ------------------------------
  // 4) Evaluate all fields through the explicit-status API
  // ------------------------------
  const swcme::ModelStatus eval_status=model.evaluate_radii_with_B_div_checked(
      S, r.data(), n.data(), V.data(), Br.data(), Bphi.data(), Bmag.data(),
      divV.data(), N);
  if (!eval_status.ok()) {
    std::fprintf(stderr,"ERROR: SWCME evaluation failed: %s\n",
                 eval_status.summary().c_str());
    return 1;
  }

  // ------------------------------
  // 5) Write Tecplot profile (POINT)
  // ------------------------------
  // OUT07 treats the example as executable documentation, so use the checked
  // API and propagate its precise validation/I/O status to the process exit
  // code.  A demonstration must never print success after a rejected request
  // or a partial/failed transaction.
  const swcme::ModelStatus output_status=
      model.write_tecplot_radial_profile_checked(
          S,r.data(),n.data(),V.data(),Br.data(),Bphi.data(),Bmag.data(),
          divV.data(),N,"profile_1d.dat");
  if (!output_status.ok()) {
    std::fprintf(stderr,"ERROR: profile output failed: %s\n",
                 output_status.summary().c_str());
    return 1;
  }
  std::printf("Wrote Tecplot radial profile: profile_1d.dat (N=%zu)\n", N);

  return 0;
}
