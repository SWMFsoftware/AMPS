//======================================================================================
// CutoffRigidityGridless.h
//======================================================================================
// PUBLIC INTERFACE FOR THE GRIDLESS CUTOFF RIGIDITY SOLVER (MPI-CAPABLE)
//
// ROLE IN srcEarth
//   This header declares the entry point used by srcEarth/main.cpp when the
//   command-line option '-mode gridless' is selected. The corresponding .cpp
//   file implements a cutoff-rigidity solver that directly evaluates magnetic
//   field models (IGRF + Tsyganenko T96/T05) during particle tracing rather
//   than interpolating from a prebuilt 3D field grid.
//
// WHAT THE SOLVER COMPUTES
//   For each requested location (explicit points or shell lon/lat grid), the
//   solver estimates the cutoff rigidity Rc by sampling many directions and
//   determining the minimum rigidity that allows backtraced particles to escape
//   the computational domain before reaching an inner loss sphere. It also
//   reports Emin (minimum kinetic energy corresponding to Rc) for convenience.
//
// MPI / LOAD BALANCING DESIGN
//   Orbit tracing cost varies significantly between locations and directions.
//   To handle this efficiently in parallel runs, the implementation uses a
//   dynamic master/worker schedule at the trajectory-task level
//   ((location_id, direction_id), plus optional directional-map cells). This avoids the
//   poor utilization that often occurs with static block decomposition when trajectory
//   runtimes vary strongly.
//
// INPUT/OUTPUT CONTRACT (through EarthUtil::AmpsParam)
//   Input (parser-normalized):
//     - Geometry distances in kilometers (km)
//     - Field model selection and parameters (T96/T05, epoch, PARMOD values)
//     - Species properties (charge, mass)
//     - Numerical controls (dt, max steps, max trace time)
//     - Output mode (POINTS or SHELLS) and requested locations
//   Output:
//     - Tecplot ASCII files written by rank 0
//     - Console summary (rank 0)
//
// MAINTENANCE NOTES
//   The single exported function below intentionally keeps the public API small.
//   Internal helper routines (field evaluation, Boris pusher, MPI worker loops,
//   Tecplot writers) are private to the implementation file so they can evolve
//   without forcing rebuilds of unrelated srcEarth modules.
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
  }
}

#endif
