// demo1d.cpp — 1-D SW+CME usage example (sanitized outputs, Tecplot profile)
//
// WHAT THIS EXAMPLE DOES (PHYSICS):
//  • Sets a “typical fast CME”: r0=1.05 Rs, V0≈1800 km/s, ambient Vsw=400 km/s.
//  • Computes a snapshot at t = 36 h after launch using DBM apex kinematics.
//  • Builds a radial grid 0.2–1.5 AU and evaluates:
//       n(r), V(r), Br(r), Bphi(r), |B|(r), divV(r).
//  • Writes a Tecplot POINT file "profile_1d.dat" with columns:
//       R[m] n V Br Bphi |B| divV rc R_sh R_LE R_TE
//
// WHAT THIS EXAMPLE DOES (IMPLEMENTATION):
//  • Uses the model’s StepState cache (constants precomputed for the time).
//  • Calls the high-throughput evaluator (no allocations in hot path).
//  • All arrays returned are NaN/Inf-free by design (API-level sanitation).
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
  P.r0_Rs       = 1.05;     // start just above the photosphere
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
  //    Midpoint sampling helps suppress FD noise in divV
  // ------------------------------
  const std::size_t N = 1200;
  std::vector<double> r(N), n(N), V(N), Br(N), Bphi(N), Bmag(N), divV(N);

  const double rmin = 0.2*AU, rmax = 1.5*AU;
  for (std::size_t i=0;i<N;++i){
    const double s = (i + 0.5) / double(N);
    r[i] = rmin + s*(rmax - rmin);
  }

  // ------------------------------
  // 4) Evaluate all fields (fast path; API-sanitized)
  // ------------------------------
  model.evaluate_radii_with_B_div(S, r.data(), n.data(), V.data(),
                                  Br.data(), Bphi.data(), Bmag.data(), divV.data(), N);

  // ------------------------------
  // 5) Write Tecplot profile (POINT)
  // ------------------------------
  if (!model.write_tecplot_radial_profile(S, r.data(), n.data(), V.data(),
                                          Br.data(), Bphi.data(), Bmag.data(), divV.data(),
                                          N, "profile_1d.dat"))
  {
    std::fprintf(stderr,"ERROR: failed to write profile_1d.dat\n");
    return 1;
  }
  std::printf("Wrote Tecplot radial profile: profile_1d.dat (N=%zu)\n", N);

  return 0;
}

