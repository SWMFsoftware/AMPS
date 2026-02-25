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
//   dynamic master/worker schedule (one location per task). This avoids the
//   poor utilization that often occurs with static block decomposition.
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

namespace Earth {
  namespace GridlessMode {
    // Execute the gridless cutoff-rigidity workflow described above.
    // Returns 0 on success; throws std::runtime_error on invalid input
    // or runtime failures (file I/O, unsupported field model, etc.).
    int RunCutoffRigidity(const EarthUtil::AmpsParam& p);
  }
}

#endif
