//======================================================================================
// CutoffRigidityMode3D.h
//======================================================================================
//
// PUBLIC INTERFACE FOR THE 3D (AMR-MESH) CUTOFF RIGIDITY SOLVER
//
//======================================================================================
// ROLE IN srcEarth
//======================================================================================
//
// This module implements geomagnetic cutoff rigidity computation inside the Mode3D
// (full PIC / AMR-mesh) execution branch.  It is called by Mode3D::Run() after the
// AMPS mesh has been initialised and replicated on every MPI rank.
//
// The gridless solver (srcEarth/gridless/CutoffRigidityGridless.cpp) solves the same
// physics problem but bypasses the AMPS mesh entirely: it evaluates B(x) by calling
// the Tsyganenko Fortran libraries directly at each integration step.  Mode3D differs
// in the following ways:
//
//   1. Domain geometry is taken from the AMPS AMR tree (ParsedDomainMin/Max set by
//      ApplyParsedDomain in Mode3D.cpp), not from the AMPS_PARAM.in km values.
//      This guarantees that the cutoff tracer uses exactly the same bounding box and
//      inner-sphere radius as the rest of the Mode3D framework.
//
//   2. Field evaluation is still performed by direct Tsyganenko / IGRF library calls
//      (exactly as in the gridless solver), NOT by interpolating mesh-stored values.
//      This is intentional: trajectory points may carry individual epoch strings and
//      field-model parameters (TRAJECTORY mode), so pre-storing one set of field
//      values in the mesh cells would be incorrect.  Direct evaluation allows
//      per-point Geopack re-initialisation with zero extra memory.
//
//   3. The AMR mesh initialisation follows a REPLICATED-DOMAIN model:
//      every MPI rank independently owns the COMPLETE mesh (all AMR blocks).
//      This differs from the standard AMPS production mode where blocks are
//      distributed across ranks.  The replicated model is required by the
//      cutoff-rigidity physics because:
//        (a) each trajectory point may carry a different field/epoch configuration,
//            making it impossible to broadcast pre-computed field snapshots, and
//        (b) each backtraced trajectory is an independent calculation that never
//            needs to communicate with another rank.
//
//======================================================================================
// PARALLELISATION STRATEGY
//======================================================================================
//
// Three-level hierarchy (coarsest to finest):
//
//   Level 1 — MPI (inter-process, replicated domain)
//   -------------------------------------------------
//   The full input point list is divided into nRanks contiguous slices.  Rank k
//   processes points [k*stride, (k+1)*stride).  Because every rank has the complete
//   mesh and a private field evaluator, no MPI communication is needed during
//   the computation phase.  Rank 0 gathers scalar results (Rc, Emin) at the end
//   and writes the Tecplot output files.
//
//   Why static slicing rather than master/worker dynamic scheduling?
//     In the gridless solver, individual-trajectory cost dominates the total runtime
//     because the gridless field call (Fortran) is cheap but the step count varies
//     by 4+ orders of magnitude (quick escape vs quasi-trapped).  Dynamic scheduling
//     absorbs the resulting load imbalance.
//
//     In Mode3D the field evaluator is identical (same Fortran calls), so the same
//     load-imbalance argument applies.  However, Mode3D adds OpenMP as a second
//     level of parallelism (Level 2 below).  With OpenMP inside each rank, the
//     per-rank work becomes a bag of many parallel tasks rather than a single long
//     serial task.  OpenMP's own internal scheduling (guided or dynamic) absorbs
//     the intra-rank imbalance, so MPI dynamic scheduling is a lower-priority
//     optimisation; static slicing is correct and sufficient.
//     If load imbalance across ranks is later observed to be a bottleneck, a
//     master/worker protocol identical to the gridless solver can be added here
//     without any change to the physics or OpenMP code.
//
//   Level 2 — OpenMP (intra-process, shared memory)
//   ------------------------------------------------
//   Within each MPI rank, trajectory points are processed in parallel by
//   OpenMP threads.  The parallel loop is over the rank's local point slice.
//
//   Thread safety:
//     - The AMPS mesh tree is READ-ONLY during the cutoff phase (only geometry
//       checks; no field-data writes after InitMeshFields in Mode3D::Run).
//       Concurrent reads from multiple threads are safe.
//     - Each thread owns a PRIVATE cMode3DFieldEval instance so that per-point
//       Geopack re-initialisation (epoch/PARMOD update) never touches shared state.
//     - The Tsyganenko Fortran libraries have global state (common blocks). They
//       are NOT thread-safe by default.  The per-thread cMode3DFieldEval instance
//       therefore acquires a fine-grained mutex only around the Fortran call
//       (see implementation notes in CutoffRigidityMode3D.cpp).
//     - Per-thread result arrays are private; only the final reduction (min Rc over
//       directions) is written to the shared result vector, which is guarded by an
//       omp critical section.
//
//   Level 3 — Bisection search (per direction, per point)
//   ------------------------------------------------------
//   For each (point, direction) pair the cutoff rigidity is found by binary search
//   in [Rmin, Rmax].  Each bisection iteration launches one backtraced trajectory
//   integration.  These are sequential within one OpenMP task.
//
//======================================================================================
// REPLICATED DOMAIN INITIALISATION PROTOCOL
//======================================================================================
//
// Called from Mode3D::Run() before RunCutoffRigidity():
//
//   1. Save  PIC::nTotalThreads  and  PIC::ThisThread.
//   2. Temporarily set both to 1 and 0, respectively.
//   3. Call amps_init_mesh().
//      With nTotalThreads=1 the mesh partitioner assigns ALL AMR blocks to thread 0
//      (= this rank).  Every rank therefore builds the complete tree.
//   4. Restore PIC::nTotalThreads and PIC::ThisThread.
//   5. Call InitMeshFields() (writes B/E into mesh cells for non-cutoff-rigidity
//      Mode3D products; the cutoff rigidity solver ignores stored values).
//   6. Call RunCutoffRigidity().
//
// This protocol requires no changes to the AMPS mesh library itself.
//
//======================================================================================
// PARTICLE TRACING (BACKTRACING)
//======================================================================================
//
// Identical physics to TraceAllowedImpl in CutoffRigidityGridless.cpp:
//
//   - Initial momentum:    p = R_GV * 1e9 * |q| / c  * v_unit
//   - Integration:        relativistic Boris pusher with adaptive dt
//   - dt_gyro = GYRO_ANGLE_LIMIT / omega_c     (limits rotation angle per step)
//   - dt_geo  = TRAVEL_FRACTION * d_nearest / v (limits overshoot near boundaries)
//   - dt      = min(dt_gyro, dt_geo, user cap, time remaining)
//   - Inner boundary: loss sphere of radius _EARTH__RADIUS_ centred at origin
//   - Outer boundary: rectangular box from ParsedDomainMin/Max
//   - ALLOWED if trajectory exits the outer box before hitting the inner sphere
//   - FORBIDDEN otherwise (inner sphere hit, time cap, step cap)
//
//======================================================================================
// INPUT / OUTPUT CONTRACT
//======================================================================================
//
// Input (EarthUtil::AmpsParam, parsed from AMPS_PARAM.in):
//   - #BACKGROUND_FIELD: FIELD_MODEL, EPOCH, field parameters (T96/T05/TA16/DIPOLE)
//   - #PARTICLE_SPECIES: charge, mass
//   - #CUTOFF_RIGIDITY:  energy range, max particles per point, sampling mode, etc.
//   - #OUTPUT_DOMAIN:    OUTPUT_MODE (POINTS | SHELLS | TRAJECTORY), coordinates
//   - #NUMERICAL:        DT_TRACE, MAX_STEPS, MAX_TRACE_TIME
//
//   For TRAJECTORY output mode, each point in prm.output.trajectories[0].samples
//   carries its own timeUTC string; RunCutoffRigidity calls field.ReinitGeopack()
//   with that per-point epoch before tracing.
//
// Output (written by rank 0 only):
//   POINTS / TRAJECTORY mode:
//     cutoff_3d_points.dat
//       TITLE + VARIABLES (id, x_km, y_km, z_km, lon_deg, lat_deg, Rc_GV, Emin_MeV)
//       ZONE per run; one row per observation point.
//   SHELLS mode:
//     cutoff_3d_shells.dat
//       TITLE + VARIABLES (lon_deg, lat_deg, Rc_GV, Emin_MeV)
//       One ZONE per shell altitude.
//
//======================================================================================

#ifndef _SRC_EARTH_3D_CUTOFFRIGIDITYMODE3D_H_
#define _SRC_EARTH_3D_CUTOFFRIGIDITYMODE3D_H_

#include "../util/amps_param_parser.h"

namespace Earth {
namespace Mode3D {

//--------------------------------------------------------------------------------------
// RunCutoffRigidity
//
// Top-level entry point for the Mode3D cutoff rigidity workflow.
// Must be called AFTER amps_init_mesh() has completed in the REPLICATED-DOMAIN
// configuration (nTotalThreads=1 during mesh init).
//
// Returns 0 on success.
// Throws std::runtime_error on invalid input or runtime failures.
//--------------------------------------------------------------------------------------
int RunCutoffRigidity(const EarthUtil::AmpsParam& prm);

} // namespace Mode3D
} // namespace Earth

#endif // _SRC_EARTH_3D_CUTOFFRIGIDITYMODE3D_H_
