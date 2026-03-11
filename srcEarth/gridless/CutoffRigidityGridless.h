//======================================================================================
// CutoffRigidityGridless.h
//======================================================================================
//
// PUBLIC INTERFACE FOR THE GRIDLESS CUTOFF RIGIDITY SOLVER (MPI-CAPABLE)
//
//======================================================================================
// ROLE IN srcEarth
//======================================================================================
//
// This header declares two categories of symbols:
//
//   (1) Earth::GridlessMode::RunCutoffRigidity
//       The top-level solver entry point called by srcEarth/main.cpp when
//       '-mode gridless' is selected. Computes cutoff rigidities at all
//       requested locations (explicit POINTS or SHELLS grid), writes Tecplot
//       ASCII output, and returns.
//
//   (2) Earth::GridlessMode::TraceAllowedShared  /  TraceAllowedSharedEx
//       Low-level single-trajectory classifiers that are also called by
//       DensityGridless and any future gridless module. Exporting them here
//       (rather than keeping them file-local) enforces a single canonical
//       definition of "allowed trajectory" across all physics products.
//
//======================================================================================
// PHYSICS: GEOMAGNETIC CUTOFF RIGIDITY
//======================================================================================
//
// The geomagnetic cutoff rigidity Rc at a location x0 is defined as:
//
//   Rc = the minimum rigidity R such that a particle arriving at x0 from
//        some direction can have originated from outside the magnetosphere.
//
// In the Stormer (1955) approximation for a pure dipole field, Rc has a closed-form
// solution. For realistic magnetospheric field models (IGRF + T96/T05) there is no
// analytic solution; Rc must be determined numerically by backtracing.
//
// The Penumbra
// -----------
// Between the Main Cone (lowest R where all directions are allowed) and the Stormer
// Cone (highest R where all directions are forbidden), there is a "penumbra" of
// alternating allowed and forbidden windows. Our effective cutoff is the minimum R
// above which the specific direction being tested is allowed (i.e., we scan rigidity
// downward and record the first forbidden transition). This matches the definition
// used by Cooke et al. (1991) and implemented in MAGNETOCOSMICS.
//
// Backtracing Principle (Liouville's theorem)
// -------------------------------------------
// By time-reversal symmetry of the Lorentz force law, a particle that can arrive
// at x0 from direction d is equivalent to a test particle launched from x0 in
// direction -d that can escape the magnetosphere. We exploit this:
//
//   A particle of rigidity R arriving at x0 from direction d is ALLOWED if and
//   only if a reversed particle launched from x0 in direction -d with the same
//   rigidity escapes the outer domain box before hitting the inner loss sphere.
//
// This converts the "can it arrive?" question into the computationally tractable
// "can it escape?" question, which requires only a forward integration.
//
//======================================================================================
// SOLVER ALGORITHM (OVERVIEW)
//======================================================================================
//
// For each observation point x0:
//
//   Step 1 -- Direction set.
//     Build a uniform sky grid of N_dirs directions (typically 24 x 48 = 1152).
//     For each direction independently:
//
//   Step 2 -- Rigidity scan.
//     Use a bisection search in [Rmin, Rmax] (from the #CUTOFF_RIGIDITY section)
//     to find the cutoff rigidity for this direction. At each trial rigidity R:
//       (a) Convert R [GV] to SI momentum: p = R*1e9*|q|/c
//       (b) Launch a reversed particle from x0 in direction -d with momentum p.
//       (c) Integrate with Boris pusher + adaptive dt until:
//             i.  Particle escapes the outer domain box -> ALLOWED
//             ii. Particle hits the inner loss sphere   -> FORBIDDEN
//             iii.Integration time/step limit reached   -> FORBIDDEN (conservative)
//       (d) Narrow the bisection interval based on allowed/forbidden outcome.
//
//   Step 3 -- Aggregate.
//     Take Rc = min over all directions of the per-direction cutoff (or the
//     isotropic fraction, depending on CUTOFF_SAMPLING mode).
//     Compute Emin = kinetic energy corresponding to Rc.
//
//   Step 4 -- Output.
//     Rank 0 collects results from all workers and writes Tecplot files.
//
//======================================================================================
// MPI LOAD-BALANCING DESIGN
//======================================================================================
//
// Trajectory integration cost varies by orders of magnitude:
//   - A high-latitude point at low energy: most directions are allowed and escape
//     quickly (few steps). Total cost: ~10 ms.
//   - An equatorial point near the cutoff in disturbed T05 fields: many near-cutoff
//     trajectories perform long quasi-trapped orbits before eventually escaping.
//     Total cost: ~10 s.
//
// A static block decomposition (rank k handles points k, k+N, k+2N, ...) would
// severely load-imbalance the MPI job when the observation grid spans a range of
// latitudes.
//
// We use dynamic master/worker scheduling instead:
//   - Rank 0 is the master: it maintains a work queue of observation-point indices.
//   - Ranks 1..nRanks-1 are workers: each repeatedly requests a task, processes it
//     (all directions and all rigidities for that point), and returns the result.
//   - Tasks are dispatched one point at a time, so slow points do not block fast
//     workers from taking new work.
//   - No assumptions about the cost distribution are made; the scheduler is purely
//     reactive.
//
// Message protocol:
//   TaskMsg:    { int pointIdx }        (master -> worker)
//   ResultMsg:  { int pointIdx, double Rc, double Emin }  (worker -> master)
//   Sentinel:   pointIdx = -1 signals no more work (master -> worker termination).
//
// Serial path (nRanks == 1 or MPI unavailable): rank 0 processes all points
// directly without any message passing.
//
//======================================================================================
// SHARED TRAJECTORY CLASSIFIER: TraceAllowedShared AND TraceAllowedSharedEx
//======================================================================================
//
// WHY THESE FUNCTIONS EXIST AS A PUBLIC API
// ------------------------------------------
// The DensityGridless module needs to classify trajectories using exactly the same
// physics as the cutoff solver: same field model, same mover, same adaptive dt,
// same inner/outer boundary geometry. If DensityGridless had its own private copy
// of the tracing loop, any bug fix or improvement to the cutoff solver would silently
// fail to propagate, causing the two physics products to diverge.
//
// Exporting the classifier here ensures a single source of truth.
//
// TraceAllowedShared
// ------------------
// The basic classifier: returns true (ALLOWED) or false (FORBIDDEN).
// Used by the isotropic density branch, which needs only the outcome.
//
// TraceAllowedSharedEx
// --------------------
// Extended version: same physics, but when the trajectory is ALLOWED also fills
// a TrajectoryExitState struct with:
//   - x_exit_m[3]   : GSM position where the particle crossed the outer domain [m]
//   - v_exit_unit[3]: velocity unit vector at that crossing
//   - cosAlpha      : cos(alpha) = v_exit . B_hat(x_exit), the pitch angle cosine
//
// The extra field evaluation at the exit point (needed for B_hat(x_exit)) costs
// one additional GetB_T call per allowed trajectory.
//
// Used by the ANISOTROPIC density branch to weight each allowed trajectory by
// f_PAD(cosAlpha) * f_spatial(x_exit) instead of the flat weight 1.
//
// If exitState is nullptr, TraceAllowedSharedEx is exactly equivalent to
// TraceAllowedShared (same code path, same cost; the exit state is computed
// but not stored).
//
// TrajectoryExitState fields are only meaningful when the function returns true.
// When it returns false (FORBIDDEN), all fields are zero-initialised and must
// not be read by the caller.
//
//======================================================================================
// INPUT / OUTPUT CONTRACT
//======================================================================================
//
// Input (all through EarthUtil::AmpsParam, parsed from AMPS_PARAM.in):
//   - Domain geometry: outer rectangular box [km] + inner loss sphere radius [km]
//   - Field model: FIELD_MODEL = T96 | T05 | DIPOLE + model parameters
//   - Epoch: ISO datetime string used to initialise Geopack (RECALC)
//   - Species: charge [e], mass [amu]
//   - Numerical: initial dt [s], max steps, max trace time [s]
//   - Output mode: POINTS (explicit list) or SHELLS (lon/lat grid at given altitudes)
//   - Cutoff scan: energy range, sampling strategy (VERTICAL | ISOTROPIC)
//
// Output (Tecplot ASCII, written by rank 0):
//   POINTS mode:
//     gridless_points_cutoff.dat
//       Variables: X_km Y_km Z_km R_GSM_km Lon_deg Lat_deg Rc_GV Emin_MeV
//   SHELLS mode:
//     gridless_shell_Akm_cutoff.dat  (one file per shell altitude A)
//       ZONE per shell; same variables as POINTS.
//
//======================================================================================

#ifndef _SRC_EARTH_GRIDLESSM_CUTOFFRIGIDITYGRIDLESS_H_
#define _SRC_EARTH_GRIDLESSM_CUTOFFRIGIDITYGRIDLESS_H_

#include "util/amps_param_parser.h"
#include "GridlessParticleMovers.h"

namespace Earth {
  namespace GridlessMode {
    // Execute the gridless cutoff-rigidity workflow described above.
    // Returns 0 on success; throws std::runtime_error on invalid input
    // or runtime failures (file I/O, unsupported field model, etc.).
    int RunCutoffRigidity(const EarthUtil::AmpsParam& p);

    // Shared low-level trajectory classifier used by both the cutoff-rigidity and
    // density/spectrum workflows.
    //
    // WHY THIS FUNCTION EXISTS
    // ------------------------
    // Earlier revisions kept the actual particle-tracking loop hidden as a file-
    // local static helper inside CutoffRigidityGridless.cpp and then reimplemented
    // a similar loop inside DensityGridless.cpp.  That duplication turned out to be
    // fragile: once one copy drifted (different dt policy, different mover, missed
    // inner-sphere crossing checks, etc.), the two physics products no longer used
    // the same definition of "allowed trajectory."
    //
    // To prevent that divergence, the tracing logic is now exported through this
    // single helper.  Both workflows call the same routine, so any future fix to
    // the mover / adaptive-step / geometry-classification logic automatically
    // affects both cutoff and density calculations in exactly the same way.
    //
    // ARGUMENTS
    // ---------
    //   p                     : full parsed AMPS/gridless parameter block
    //   x0_m[3]               : starting point in GSM Cartesian coordinates [m]
    //   v0_unit[3]            : starting direction (must be normalized or close)
    //   R_GV                  : particle rigidity [GV]
    //   maxTraceTimeOverride_s: optional per-trajectory time cap.  If >0, it
    //                           overrides the default cutoff-specific cap inside
    //                           the shared tracer.  Density/spectrum uses this to
    //                           apply DS_MAX_TRAJ_TIME while still reusing the
    //                           exact same particle push / escape / loss logic as
    //                           the cutoff solver.
    //
    // RETURN VALUE
    // ------------
    //   true  -> trajectory escapes the outer domain before entering the inner sphere
    //   false -> trajectory hits the inner sphere or fails to escape before limits
    //
    // PERFORMANCE NOTE
    // ----------------
    // The implementation caches the field evaluator per thread/run so repeated
    // calls from the density module do not reconstruct the field context on every
    // trajectory.  That keeps this shared-logic design practical even though the
    // function interface itself stays simple.
    bool TraceAllowedShared(const EarthUtil::AmpsParam& p,
                            const double x0_m[3],
                            const double v0_unit[3],
                            double R_GV,
                            double maxTraceTimeOverride_s=-1.0);

    //-------------------------------------------------------------------------
    // Extended trajectory result returned by TraceAllowedSharedEx.
    //
    // Fields are only meaningful when TraceAllowedSharedEx returns true.
    // When the trajectory is FORBIDDEN (returns false), the fields are
    // left at zero/default and must not be used by the caller.
    //-------------------------------------------------------------------------
    struct TrajectoryExitState {
      double x_exit_m[3];    // GSM position where the trajectory crossed the
                             // outer domain boundary [m].
      double v_exit_unit[3]; // Velocity unit vector at that crossing point.
                             // Points in the direction of particle travel (outward).
      double cosAlpha;       // cos of pitch angle at exit:
                             //   cos(alpha) = v_exit_unit . B_hat(x_exit)
                             // where B_hat is the unit field vector at x_exit.
                             // Range [-1, 1].  Sign convention follows the
                             // standard geophysics convention (cosAlpha > 0 means
                             // particle moves along the field).
    };

    // Identical physics to TraceAllowedShared but also fills *exitState when the
    // trajectory is allowed.  exitState may be nullptr; if so this call is
    // equivalent to TraceAllowedShared (the exit state is computed but not stored,
    // with one extra field evaluation at the exit point).
    //
    // Used by the ANISOTROPIC density branch in DensityGridless to obtain the
    // asymptotic direction and exit location of each allowed trajectory so that
    // J_b(E, Omega_inf, x_inf) can be evaluated per-trajectory rather than
    // using a single isotropic J_b(E).
    bool TraceAllowedSharedEx(const EarthUtil::AmpsParam& p,
                              const double x0_m[3],
                              const double v0_unit[3],
                              double R_GV,
                              TrajectoryExitState* exitState,
                              double maxTraceTimeOverride_s=-1.0);
  }
}

#endif
