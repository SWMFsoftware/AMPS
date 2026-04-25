//======================================================================================
// CutoffRigidityMode3D.cpp
//======================================================================================
//
// See CutoffRigidityMode3D.h for the full design rationale, parallelisation strategy,
// and input/output contract.  This file documents implementation specifics.
//
//======================================================================================
// THREAD SAFETY — TSYGANENKO FORTRAN COMMON BLOCKS
//======================================================================================
//
// T96, T05, T01, TA15, TA16, and Geopack store their state in Fortran common blocks
// (IGRF coefficients, parmod, dipole tilt, etc.).  A common block is a single global
// region of static storage; it is NOT thread-safe.
//
// Field evaluation in Mode3D reads from the AMPS AMR mesh (no Fortran calls,
// field-evaluation call.  The lock is released immediately afterward, so the critical
// section contains only the Fortran call itself (typically a few microseconds).
//
// Why a single mutex rather than per-model mutexes?
//   All models share the same Geopack common block (RECALC state, dipole tilt).
//   Separate per-model mutexes would not protect Geopack from concurrent writes by
//   two threads using different models in a mixed-model run.  A single mutex is the
//   conservative, correct choice.
//
// PERFORMANCE NOTE:
//   In a typical run with 1152 directions per point and 20 bisection steps, one point
//   requires ~23 000 Fortran calls.  At ~5 µs/call, a 16-core node running 16 OpenMP
//   threads in lock-step would spend ~7 ms in the critical section per point, while
//   each point takes ~115 ms total.  Lock contention is therefore ~6% overhead at 16
//   threads.  For heavier models (TA16) the fraction is lower because non-Fortran
//   work (Boris pusher, geometry checks) is also more expensive.
//
//   If Fortran-mutex contention ever becomes the bottleneck, the remedy is to link
//   thread-safe versions of the Fortran libraries (compiled with -frecursive or
//   equivalent) and remove the mutex.  That change is isolated to this file.
//
//======================================================================================
// ADAPTIVE TIME STEP
//======================================================================================
//
// Identical policy to SelectAdaptiveDt in CutoffRigidityGridless.cpp:
//   dt = min(DT_TRACE_max,
//            GYRO_ANGLE_LIMIT / omega_c,
//            TRAVEL_FRACTION * d_nearest / v,
//            time_remaining)
//
// GYRO_ANGLE_LIMIT = 0.15 rad   (limits Boris rotation error per step)
// TRAVEL_FRACTION  = 0.20       (limits overshoot near domain boundaries)
//
//======================================================================================
// DOMAIN GEOMETRY
//======================================================================================
//
// Outer box:    Mode3D::ParsedDomainMin[3] / ParsedDomainMax[3]  (SI meters)
//               Set by ApplyParsedDomain() in Mode3D.cpp from prm.domain (km → m).
//
// Inner sphere: radius = _EARTH__RADIUS_  (Earth radius in SI meters)
//               Centred at origin (GSM).
//
// These match the geometry established by amps_init_mesh() so that the cutoff tracer
// and the AMPS particle mover see the same physical domain.
//
//======================================================================================
// RIGIDITY SEARCH
//======================================================================================
//
// Identical to the gridless solver's CutoffForDirection_GV:
//
//   1. Quick endpoint checks at Rmin and Rmax.
//   2. Binary search while (hi - lo) > max(1e-3 GV, 1e-6 * max(|Rmin|,|Rmax|)).
//   3. Return hi (first allowed rigidity from below).
//
// Per-location upper-bound optimisation:
//   Once a direction with cutoff Rc_found is known for a point, later directions
//   for the SAME point search only up to min(Rmax, Rc_found).  This is safe because
//   the isotropic cutoff is the minimum over directions; any direction with a higher
//   individual cutoff cannot lower the minimum below Rc_found.
//
//======================================================================================
// MPI IMPLEMENTATION NOTES
//======================================================================================
//
// MPI is always compiled in (linked with mpi.h).  We guard MPI_Init with
// MPI_Initialized() so that runs launched under mpirun (which calls MPI_Init
// before main()) are not broken.
//
// Point distribution: static block decomposition.
//   nLocal = nLoc / mpiSize  (floor)
//   rank k owns points [k*nLocal, (k+1)*nLocal)
//   last rank picks up the remainder: [(mpiSize-1)*nLocal, nLoc)
//
// Communication:
//   After the computation phase, MPI_Gatherv collects the per-rank result arrays
//   (Rc, Emin) at rank 0.  No communication during computation.
//
//======================================================================================

#include "CutoffRigidityMode3D.h"
#include "Mode3D.h"
#include "ElectricField.h"

// AMPS framework
#include "pic.h"
#include "Earth.h"
#include "specfunc.h"

// Gridless shared utilities (movers, V3 helpers, IGridlessFieldEvaluator)
#include "../gridless/GridlessParticleMovers.h"

// Shared Störmer formula and DipoleInterface global state (gParams.m_hat).
// CutoffRigidityGridless.h declares StormerVerticalCoeff_GV inline so that
// Mode3D and gridless use a single definition with identical constants.
#include "../gridless/CutoffRigidityGridless.h"
#include "../gridless/DipoleInterface.h"

// Parameter parser
#include "../util/amps_param_parser.h"

// Standard library
#include <cstdio>
#include <cmath>
#include <cstring>
#include <string>
#include <vector>
#include <algorithm>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <iostream>

// MPI (always compiled in; see design notes in header)
#include <mpi.h>

// OpenMP
#ifdef _OPENMP
#include <omp.h>
#endif

// SPICE — only for coordinate frame transforms in LocationToX0m (SHELLS mode).
// NOT used for field evaluation: B is read from the AMPS mesh cells.
#ifndef _NO_SPICE_CALLS_
#include "SpiceUsr.h"
#endif

// Constants
#include "constants.h"
#include "constants.PlanetaryData.h"

// NOTE: Tsyganenko model interfaces (T96Interface.h, T05Interface.h, etc.),
// GeopackInterface.h, and DipoleInterface.h are intentionally NOT included here.
// In Mode3D the magnetic field is read directly from the AMPS AMR mesh cells
// that were populated by InitMeshFields() in Mode3D::Run() before the cutoff
// calculation begins.  There is therefore no need to call the Tsyganenko Fortran
// libraries during particle tracing, and no mutex is required.

namespace {

//======================================================================================
//======================================================================================
// SECTION 1 — THREAD SAFETY
//======================================================================================
//
// The magnetic field is read from the AMPS AMR mesh cells that were populated by
// InitMeshFields() in Mode3D::Run() before the cutoff calculation begins.
//
// findTreeNode() is a pure pointer traversal over the frozen (read-only) AMR tree.
// Concurrent reads from multiple OpenMP threads are safe: no shared state is
// modified and no lock is required.
//
// This eliminates the per-step mutex that the gridless evaluator needs around its
// Tsyganenko Fortran calls, giving substantially better OpenMP scaling.
//======================================================================================

//======================================================================================
// SECTION 2 — SMALL VECTOR HELPERS
//======================================================================================
// V3 and operators are imported from GridlessParticleMovers.h.
// We add only those helpers not already provided there.

static inline double v3norm(const V3& a) {
    return std::sqrt(a.x*a.x + a.y*a.y + a.z*a.z);
}

static inline V3 v3unit(const V3& a) {
    const double n = v3norm(a);
    return (n > 0.0) ? mul(1.0/n, a) : V3{0.0, 0.0, 0.0};
}

//======================================================================================
// SECTION 3 — UNIT CONVERSIONS (identical to gridless solver)
//======================================================================================

static inline double MomentumFromKineticEnergy_MeV(double E_MeV, double m0_kg) {
    const double E_J = E_MeV * 1.0e6 * ElectronCharge;
    return Relativistic::Energy2Momentum(E_J, m0_kg);
}

static inline double KineticEnergyFromMomentum_MeV(double p, double m0_kg) {
    const double E_J = Relativistic::Momentum2Energy(p, m0_kg);
    return E_J / (1.0e6 * ElectronCharge);
}

static inline double MomentumFromRigidity_GV(double R_GV, double q_C_abs) {
    return (R_GV * 1.0e9 * q_C_abs) / SpeedOfLight;
}

static inline double RigidityFromMomentum_GV(double p, double q_C_abs) {
    return (q_C_abs > 0.0) ? (p * SpeedOfLight / q_C_abs / 1.0e9) : 0.0;
}

//======================================================================================
// SECTION 4 — DOMAIN GEOMETRY
//======================================================================================
//
// Outer domain: rectangular box in SI meters, from Mode3D::ParsedDomainMin/Max.
// Inner boundary: loss sphere of radius _EARTH__RADIUS_ centred at origin.
//
// These are set once at the start of RunCutoffRigidity and shared across all threads
// as read-only data (no mutation after initialisation).
//======================================================================================

struct DomainBox3D {
    double xMin, xMax;
    double yMin, yMax;
    double zMin, zMax;
    double rInner; // inner loss sphere radius [m]
};

static inline bool InsideBox3D(const V3& x, const DomainBox3D& b) {
    return (x.x >= b.xMin && x.x <= b.xMax &&
            x.y >= b.yMin && x.y <= b.yMax &&
            x.z >= b.zMin && x.z <= b.zMax);
}

static inline bool LostInnerSphere3D(const V3& x, double rInner) {
    const double r2 = x.x*x.x + x.y*x.y + x.z*x.z;
    return (r2 <= rInner * rInner);
}

//======================================================================================
// SECTION 5 — MESH-BASED FIELD EVALUATOR (per-thread)
//======================================================================================
//
// cMode3DMeshFieldEval reads the magnetic field that InitMeshFields() stored in the
// AMPS AMR cell-centre data buffers.  It is the Mode3D replacement for the gridless
// cFieldEvaluator (which calls Tsyganenko Fortran directly).
//
// Design:
//   - Implements IGridlessFieldEvaluator so all shared mover infrastructure
//     (StepParticleChecked, HybridPrepareStepUseGuidingCenter, etc.) is reused
//     without modification.
//   - Each OpenMP thread holds its own instance; the mutable lastNode_ member
//     caches the most recently found AMR block, accelerating findTreeNode() when
//     consecutive GetB_T calls are within the same block (the common case for a
//     particle advancing in small steps).
//   - No mutex: the AMR tree and cell data are READ-ONLY during particle tracing.
//
// Field data layout (mirrors InitMeshFields write):
//   ptr = center->GetAssociatedDataBufferPointer()
//       + PIC::CPLR::DATAFILE::CenterNodeAssociatedDataOffsetBegin
//       + PIC::CPLR::DATAFILE::MULTIFILE::CurrDataFileOffset
//       + PIC::CPLR::DATAFILE::Offset::MagneticField.RelativeOffset
//       + idim * sizeof(double)      (idim = 0,1,2 for Bx,By,Bz)
//
// Fallback: if the particle is outside all allocated blocks (e.g. just before the
// outer-boundary escape check fires), GetB_T returns B = 0.  The integration loop
// detects the escape on the very next geometry check, so the zero-field step is
// never used to advance the trajectory classification.
//======================================================================================

class cMode3DMeshFieldEval : public IGridlessFieldEvaluator {
public:
    cMode3DMeshFieldEval() : lastNode_(nullptr) {}

    // Evaluate B [Tesla] at position x_m [m] by reading from the AMR mesh cells.
    void GetB_T(const V3& x_m, V3& B_T) const override {
        double xArr[3] = {x_m.x, x_m.y, x_m.z};

        // Locate the leaf block that contains xArr.
        // The lastNode_ hint makes subsequent calls O(1) when the particle stays
        // in the same block between integration steps.
        cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node =
            PIC::Mesh::mesh->findTreeNode(xArr, lastNode_);
        lastNode_ = node;          // update hint for next call

        if (node == nullptr || node->block == nullptr) {
            // Particle is outside the domain or in an unallocated block.
            // Return zero; the caller's boundary check will handle this.
            B_T.x = B_T.y = B_T.z = 0.0;
            return;
        }

        // Determine which cell centre within the block is closest to x.
        // Cell i occupies [xmin + i*dx, xmin + (i+1)*dx] where dx = blockSize/_BLOCK_CELLS_X_.
        const double dx = (node->xmax[0] - node->xmin[0]) / _BLOCK_CELLS_X_;
        const double dy = (node->xmax[1] - node->xmin[1]) / _BLOCK_CELLS_Y_;
        const double dz = (node->xmax[2] - node->xmin[2]) / _BLOCK_CELLS_Z_;

        int i = static_cast<int>((xArr[0] - node->xmin[0]) / dx);
        int j = static_cast<int>((xArr[1] - node->xmin[1]) / dy);
        int k = static_cast<int>((xArr[2] - node->xmin[2]) / dz);

        // Clamp to interior cells [0, _BLOCK_CELLS_? - 1].
        i = std::max(0, std::min(i, _BLOCK_CELLS_X_ - 1));
        j = std::max(0, std::min(j, _BLOCK_CELLS_Y_ - 1));
        k = std::max(0, std::min(k, _BLOCK_CELLS_Z_ - 1));

        const int nd = PIC::Mesh::mesh->getCenterNodeLocalNumber(i, j, k);
        PIC::Mesh::cDataCenterNode* center = node->block->GetCenterNode(nd);

        if (center == nullptr ||
            !PIC::CPLR::DATAFILE::Offset::MagneticField.active) {
            B_T.x = B_T.y = B_T.z = 0.0;
            return;
        }

        // Read Bx, By, Bz from the data buffer using the same offset chain
        // that InitMeshFields() used to write them.
        const char* base = center->GetAssociatedDataBufferPointer()
                         + PIC::CPLR::DATAFILE::CenterNodeAssociatedDataOffsetBegin
                         + PIC::CPLR::DATAFILE::MULTIFILE::CurrDataFileOffset
                         + PIC::CPLR::DATAFILE::Offset::MagneticField.RelativeOffset;

        B_T.x = *reinterpret_cast<const double*>(base + 0*sizeof(double));
        B_T.y = *reinterpret_cast<const double*>(base + 1*sizeof(double));
        B_T.z = *reinterpret_cast<const double*>(base + 2*sizeof(double));
    }

private:
    // Mutable so GetB_T can update the search hint while remaining logically const.
    mutable cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* lastNode_;
};

//======================================================================================
// SECTION 6 — FIBONACCI SPHERE DIRECTION GRID
//======================================================================================
// Identical to BuildUniformSphereDirs in CutoffRigidityGridless.cpp.
// Duplicated here (not shared) to keep this file self-contained and independent of
// the gridless module's anonymous namespace.
//======================================================================================

static std::vector<V3> BuildFibonacciDirs(int nDir) {
    if (nDir <= 0)
        throw std::runtime_error("BuildFibonacciDirs: nDir must be > 0");

    std::vector<V3> dirs;
    dirs.reserve((size_t)nDir);

    const double goldenRatio  = 0.5 * (1.0 + std::sqrt(5.0));
    const double goldenAngle  = 2.0 * M_PI * (1.0 - 1.0 / goldenRatio);

    for (int k = 0; k < nDir; ++k) {
        const double z   = 1.0 - 2.0*(static_cast<double>(k) + 0.5)/static_cast<double>(nDir);
        const double r   = std::sqrt(std::max(0.0, 1.0 - z*z));
        const double phi = goldenAngle * static_cast<double>(k);
        dirs.push_back(v3unit(V3{r*std::cos(phi), r*std::sin(phi), z}));
    }

    return dirs;
}

//======================================================================================
// SECTION 7 — ADAPTIVE TIME STEP SELECTOR
//======================================================================================
//
// Mirrors SelectAdaptiveDt in CutoffRigidityGridless.cpp exactly, differing only in
// the domain-geometry representation: the gridless version works in Earth-radii
// (DomainBoxRe) and converts to metres when computing distances; this version
// receives DomainBox3D whose coordinates are already in SI metres, so no conversion
// factor is needed.
//
// useGuidingCenterForThisStep must be determined by the caller BEFORE calling
// SelectDt3D, for exactly the same reason as in the gridless solver:
//   - GC steps do not advance the fast gyrophase, so the gyro-angle limit is
//     irrelevant and must be skipped to avoid unnecessarily small steps.
//   - The dominant GC transport speed is |v_parallel| (field-aligned motion),
//     not the full particle speed |v|.  Using |v| for the boundary limiter in the
//     GC branch makes steps artificially small and wastes integration budget.
//
// The floor expression is identical to the gridless formula, including the
// 100 km / vForBoundaryLimiter term that guarantees forward progress even when
// the particle is nearly stationary.
//======================================================================================

static inline double SelectDt3D(const EarthUtil::AmpsParam& prm,
                                  const cMode3DMeshFieldEval& field,
                                  const V3& x, const V3& p,
                                  double q_C, double m0_kg,
                                  const DomainBox3D& box,
                                  double timeRemaining_s,
                                  bool useGuidingCenterForThisStep) {
    double dt = prm.numerics.dtTrace_s;
    if (timeRemaining_s < dt) dt = timeRemaining_s;

    // Relativistic kinematics — shared by both branches.
    const double p2    = dot(p, p);
    const double mc    = m0_kg * SpeedOfLight;
    const double gamma = std::sqrt(1.0 + p2 / (mc*mc));
    const double pMag  = std::sqrt(std::max(0.0, p2));
    const double vMag  = (gamma > 0.0 && m0_kg > 0.0) ? pMag / (gamma*m0_kg) : 0.0;

    // Evaluate B once; reused for the gyro limiter (full-orbit) and the
    // |v_parallel| boundary limiter (GC branch).
    V3 B; field.GetB_T(x, B);
    const double Bmag = v3norm(B);

    // (1) Branch-aware gyro-angle constraint.
    //
    // Full-orbit branch (Boris, RK2/4/6): limit the rotation angle per step to
    // GYRO_ANGLE_LIMIT = 0.15 rad.  This is the standard Boris accuracy criterion.
    //
    // GC branch (GC2/4/6, HYBRID in GC mode): skip entirely.  The GC equations
    // do not advance the fast gyrophase explicitly — imposing a gyrofrequency
    // limit would make GC steps unnecessarily small without any accuracy benefit.
    if (!useGuidingCenterForThisStep) {
        const double kGyroAngleLimit = 0.15;
        if (Bmag > 0.0) {
            const double omega = std::fabs(q_C) * Bmag / (gamma * m0_kg);
            if (omega > 0.0) dt = std::min(dt, kGyroAngleLimit / omega);
        }
    }

    // (2) Geometry-aware travel constraint (branch-aware speed selection).
    //
    // Limit the per-step displacement to TRAVEL_FRACTION = 20% of the distance
    // to the nearest stopping surface (outer box face or inner loss sphere).
    // This prevents large overshoot near classification boundaries.
    //
    // The effective "transport speed" differs between branches:
    //   Full-orbit branch: use full particle speed |v|.
    //   GC branch:         use |v_parallel| = |p · B̂| / (γ m₀).
    //     Rationale: in the GC approximation the guiding centre drifts along
    //     the field line at |v_parallel|; the drift velocities (grad-B, curvature)
    //     are generally much smaller and are already captured by the mover.
    //     Using |v| would impose the same small dt as the full-orbit branch and
    //     negate the efficiency gain of the GC representation.
    const double dxp = box.xMax - x.x;  const double dxm = x.x - box.xMin;
    const double dyp = box.yMax - x.y;  const double dym = x.y - box.yMin;
    const double dzp = box.zMax - x.z;  const double dzm = x.z - box.zMin;
    // dBox is the distance to the nearest outer-box face (metres).
    // If the particle is outside the box, all differences are negative; we treat
    // that as "box distance = large" because the escape check fires immediately.
    const double dBox  = std::min({dxp, dxm, dyp, dym, dzp, dzm});
    const double dBoxSafe = (dBox < 0.0) ? 1.0e300 : dBox;

    const double r      = v3norm(x);
    const double dInner = std::max(0.0, r - box.rInner); // dist to inner sphere
    const double dNear  = std::min(dBoxSafe, dInner);

    const double kTravelFraction = 0.20;
    double vForBoundaryLimiter = vMag;
    if (useGuidingCenterForThisStep && Bmag > 0.0) {
        const V3 bHat    = mul(1.0 / Bmag, B);
        const double pPar    = dot(p, bHat);
        const double vParAbs = std::fabs(pPar) / (gamma * m0_kg);
        vForBoundaryLimiter  = std::max(vParAbs, 1.0e-12);
    }
    if (vForBoundaryLimiter > 0.0 && dNear < 1.0e299)
        dt = std::min(dt, kTravelFraction * dNear / vForBoundaryLimiter);

    // Floor — identical to the gridless formula.
    // The 100 km / v term (100.0e3 / vForBoundaryLimiter) guarantees at least one
    // step across a 100 km segment regardless of how small the gyro or travel
    // limits become, preventing stagnation when |B| is very large.
    const double dtFloor = std::max(
        std::max(1.0e-12, 1.0e-9 * std::max(prm.numerics.dtTrace_s, 1.0)),
        100.0e3 / vForBoundaryLimiter);
    dt = std::max(dtFloor, dt);
    if (timeRemaining_s > 0.0) dt = std::min(dt, timeRemaining_s);

    return dt;
}

//======================================================================================
// SECTION 8 — CORE BACKTRACER
//======================================================================================
//
// TraceAllowed3D integrates one reversed particle trajectory and classifies it as
// ALLOWED (escaped outer box) or FORBIDDEN (hit inner sphere or exceeded limits).
//
// This function is a direct port of TraceAllowedImpl from CutoffRigidityGridless.cpp.
// The only structural differences are:
//   - Domain geometry uses DomainBox3D (metres) instead of DomainBoxRe (Earth radii).
//     All distance comparisons are therefore done in metres without any Re↔m
//     conversion factor.
//   - The field evaluator type is cMode3DMeshFieldEval (implements IGridlessFieldEvaluator)
//     rather than the gridless cFieldEvaluator.  Both classes expose the same
//     GetB_T(x, B) virtual method, so all mover calls are identical.
//
// Mover selection is fully shared with the gridless solver:
//   - gDefaultMover     (set by SetDefaultMoverType; honours -mover CLI flag)
//   - GetDefaultMoverType()  (used to decide per-step branch in the integration loop)
//   - HybridPrepareStepUseGuidingCenter()  (HYBRID branch decision; thread-local)
//   - ResetHybridTrajectoryContext()       (resets per-trajectory HYBRID state)
//   - StepParticleChecked()               (unified step + inner-sphere contact test)
//
// Using the same shared mover module guarantees that a run produced with
//   ./amps -mode gridless -mover HYBRID
// and one produced with
//   ./amps -mode 3d      -mover HYBRID
// apply exactly the same integration algorithm to each backtraced trajectory.
//
// Thread safety:
//   - gDefaultMover and GetDefaultMoverType() are read-only after SetDefaultMoverType
//     is called once at startup (before any parallel region).
//   - HybridPrepareStepUseGuidingCenter() reads and writes thread-local storage
//     (the HYBRID trajectory context); each OpenMP thread is safe.
//   - The field evaluator (cMode3DMeshFieldEval) is per-thread; see Section 5.
//   - StepParticleChecked calls field.GetB_T which reads from the frozen AMR mesh (no mutex).
//======================================================================================

static bool TraceAllowed3D(const EarthUtil::AmpsParam& prm,
                            cMode3DMeshFieldEval& field,
                            const V3& x0_m,
                            const V3& v0_unit,
                            double R_GV,
                            double q_C,
                            double m0_kg,
                            const DomainBox3D& box,
                            double maxTime_s = -1.0) {
    const double qabs = std::fabs(q_C);
    const double pMag = MomentumFromRigidity_GV(R_GV, qabs);

    V3 p = mul(pMag, v0_unit);
    V3 x = x0_m;

    // Per-trajectory time limit (mirrors TraceAllowedImpl priority: per-call
    // override > #CUTOFF_RIGIDITY CUTOFF_MAX_TRAJ_TIME > #NUMERICAL MAX_TRACE_TIME).
    const double maxTraceTime_s =
        (maxTime_s > 0.0)
            ? maxTime_s
            : ((prm.cutoff.maxTrajTime_s > 0.0)
               ? prm.cutoff.maxTrajTime_s
               : prm.numerics.maxTraceTime_s);

    // Optional cumulative path-length cap (MAX_TRACE_DISTANCE in Re → metres).
    const double maxDist_m =
        (prm.numerics.maxTraceDistance_Re > 0.0)
            ? prm.numerics.maxTraceDistance_Re * _EARTH__RADIUS_
            : -1.0;

    // Reset thread-local HYBRID context for this trajectory.
    // Every backtraced trajectory must start in the full-orbit (RK) branch;
    // the HYBRID mover only switches to GC after observing consecutive
    // adiabatic steps away from the launch point and the loss sphere.
    // Resetting here prevents state leakage from one bisection step to the next.
    // This call is harmless for all non-HYBRID movers.
    ResetHybridTrajectoryContext(x0_m, box.rInner);

    double tTrace = 0.0;
    int    nSteps = 0;
    double sDist  = 0.0;

    while (nSteps < prm.numerics.maxSteps &&
           tTrace < maxTraceTime_s &&
           (maxDist_m <= 0.0 || sDist < maxDist_m)) {

        // Boundary classification BEFORE the push (consistent with gridless).
        if (LostInnerSphere3D(x, box.rInner)) return false; // FORBIDDEN
        if (!InsideBox3D(x, box))             return true;  // ALLOWED

        const double timeRemaining = maxTraceTime_s - tTrace;

        // ---- Mover-branch decision (mirrors TraceAllowedImpl exactly) ----------
        //
        // The branch flag drives two downstream decisions:
        //   1. SelectDt3D: skip gyro-angle limiter and use |v_parallel| for the
        //      boundary limiter when in GC mode (see Section 7).
        //   2. StepParticleChecked: dispatches to the mover selected by gDefaultMover,
        //      which may be BORIS, RK4, GC4, HYBRID, etc. depending on the -mover flag.
        //
        // Control flow (identical to TraceAllowedImpl):
        //   - GC2/GC4/GC6: always GC → useGuidingCenterForThisStep = true.
        //   - HYBRID:       per-step decision from HybridPrepareStepUseGuidingCenter.
        //                   The function reads local adiabaticity and trajectory history
        //                   from thread-local storage; it returns the branch choice for
        //                   the UPCOMING step and updates the history counters.
        //   - All others (BORIS, RK2, RK4, RK6): full-orbit → flag stays false.
        bool useGuidingCenterForThisStep = false;
        const MoverType mover = GetDefaultMoverType();
        if (mover == MoverType::GC2 ||
            mover == MoverType::GC4 ||
            mover == MoverType::GC6) {
            useGuidingCenterForThisStep = true;
        } else if (mover == MoverType::HYBRID) {
            useGuidingCenterForThisStep =
                HybridPrepareStepUseGuidingCenter(x, p, q_C, field);
        }
        // -------------------------------------------------------------------------

        // Adaptive dt — branch-aware (see Section 7).
        const double dt = SelectDt3D(prm, field, x, p, q_C, m0_kg, box,
                                      timeRemaining, useGuidingCenterForThisStep);

        // Advance one step.  StepParticleChecked dispatches to gDefaultMover and
        // returns false immediately if the inner loss sphere is contacted during any
        // intermediate stage of the integrator (RK sub-steps or Boris half-steps).
        const V3 xPrev = x;
        if (!StepParticleChecked(gDefaultMover, x, p, q_C, m0_kg, dt, field, box.rInner)) {
            return false; // FORBIDDEN (mid-step inner-sphere contact)
        }

        sDist  += v3norm(sub(x, xPrev));
        tTrace += dt;
        ++nSteps;
    }

    // Conservative fallback: time/step/distance cap reached without escape.
    return false;
}

//======================================================================================
// SECTION 9 — CUTOFF SEARCH FOR ONE DIRECTION
//======================================================================================
//
// CutoffForDir_GV finds the cutoff rigidity for a single (point, direction) pair
// using the same binary-search algorithm as CutoffForDirection_GV in the gridless
// solver, including the interval-collapse guard and the per-location upper-bound
// optimisation hook (rHi_GV may be set to the current best Rc for this location).
//======================================================================================

static double CutoffForDir_GV(const EarthUtil::AmpsParam& prm,
                               cMode3DMeshFieldEval& field,
                               const V3& x0_m,
                               const V3& dir_unit,   // ARRIVAL direction (backtraced as -dir)
                               double q_C, double m0_kg,
                               const DomainBox3D& box,
                               double Rmin_GV, double Rmax_GV) {
    // Backtracing convention: launch the reversed particle in direction -dir_unit
    const V3 v0 = mul(-1.0, dir_unit);

    // Degenerate bracket guard
    const double tolAbs = 1.0e-3;
    const double tolRel = 1.0e-6 * std::max(std::fabs(Rmin_GV), std::fabs(Rmax_GV));
    const double tol    = std::max(tolAbs, tolRel);
    if (Rmax_GV < Rmin_GV) return -1.0;

    // Quick endpoint classification
    const bool alo = TraceAllowed3D(prm, field, x0_m, v0, Rmin_GV, q_C, m0_kg, box);
    const bool ahi = TraceAllowed3D(prm, field, x0_m, v0, Rmax_GV, q_C, m0_kg, box);

    if (alo && ahi) return Rmin_GV; // already allowed at minimum → cutoff ≤ Rmin
    if (!ahi)       return -1.0;    // still forbidden at maximum → no cutoff in bracket

    // Binary search while interval is wider than tolerance
    double lo = Rmin_GV, hi = Rmax_GV;
    while ((hi - lo) > tol) {
        const double mid = 0.5*(lo + hi);
        const bool a = TraceAllowed3D(prm, field, x0_m, v0, mid, q_C, m0_kg, box);
        if (a) hi = mid; else lo = mid;
    }
    return hi;
}

//======================================================================================
// SECTION 10 — COMPUTE CUTOFF FOR ONE OBSERVATION POINT
//======================================================================================
//
// ComputeCutoffAtPoint_GV runs the full direction scan (Fibonacci sphere or vertical)
// for one trajectory point and returns the minimum cutoff rigidity over all sampled
// directions.
//
// Parameters:
//   prm        — parsed AMPS parameters
//   field      — per-thread field evaluator (private)
//   x0_m       — observation point in GSM metres
//   dirs       — precomputed Fibonacci direction grid (ignored for VERTICAL sampling)
//   samplingVertical — if true, only the local vertical direction is tested
//   q_C, m0_kg — species charge and mass
//   box        — domain geometry
//   Rmin_GV, Rmax_GV — rigidity scan bracket
//
// Returns the minimum cutoff rigidity in GV, or −1 if no cutoff was found in
// [Rmin_GV, Rmax_GV] for any sampled direction.
//======================================================================================

static double ComputeCutoffAtPoint_GV(const EarthUtil::AmpsParam& prm,
                                       cMode3DMeshFieldEval& field,
                                       const V3& x0_m,
                                       const std::vector<V3>& dirs,
                                       bool samplingVertical,
                                       double q_C, double m0_kg,
                                       const DomainBox3D& box,
                                       double Rmin_GV, double Rmax_GV) {
    double rcMin = -1.0;

    const int nDirs = samplingVertical ? 1 : static_cast<int>(dirs.size());

    for (int d = 0; d < nDirs; ++d) {
        V3 dir;
        if (samplingVertical) {
            // Local vertical: particle arrives from the radially outward direction,
            // so the arrival direction is toward Earth: dir = -unit(x0).
            dir = v3unit(mul(-1.0, x0_m));
        } else {
            dir = dirs[(size_t)d];
        }

        // Per-location upper-bound optimisation: once Rc_found is known for this
        // point, searching higher rigidities cannot lower the isotropic minimum.
        const double rHiThis = (!samplingVertical && rcMin > 0.0)
                               ? std::min(Rmax_GV, rcMin)
                               : Rmax_GV;

        const double rc = CutoffForDir_GV(prm, field, x0_m, dir, q_C, m0_kg,
                                           box, Rmin_GV, rHiThis);
        if (rc > 0.0)
            rcMin = (rcMin < 0.0) ? rc : std::min(rcMin, rc);
    }

    return rcMin;
}

//======================================================================================
// SECTION 11 — EPOCH LOOKUP FOR TRAJECTORY POINTS
//======================================================================================

// Returns the UTC epoch string to use for observation point `loc` in the flattened
// location list.  For TRAJECTORY output mode each sample carries its own timestamp;
// for POINTS/SHELLS the global epoch is reused.
static std::string EpochForLocation(const EarthUtil::AmpsParam& prm, int loc) {
    const std::string mode = EarthUtil::ToUpper(prm.output.mode);
    if ((mode == "TRAJECTORY") && !prm.output.trajectories.empty()) {
        const auto& traj = prm.output.trajectories[0];
        if (loc >= 0 && loc < static_cast<int>(traj.samples.size()))
            return traj.samples[(size_t)loc].timeUTC;
    }
    return prm.field.epoch;
}

//======================================================================================
// SECTION 12 — LOCATION INDEX TO 3D POSITION
//======================================================================================

// LocationToX0m converts a flattened location index to a GSM Cartesian position [m].
// POINTS / TRAJECTORY: prm.output.points[loc] in km → metres.
// SHELLS: lon/lat/alt grid node → GSM metres via SPICE (or identity for DIPOLE).
//
// Logic is identical to the LocationToX0m lambda in CutoffRigidityGridless.cpp.
// Replicated here to avoid a cross-module dependency on a local lambda.
static V3 LocationToX0m(const EarthUtil::AmpsParam& prm, int loc,
                         int nLon, int nLat,
                         double d_deg, int nPtsShell) {
    const std::string mode = EarthUtil::ToUpper(prm.output.mode);
    const bool isPoints = (mode == "POINTS" || mode == "TRAJECTORY");

    if (isPoints) {
        const auto& P = prm.output.points[(size_t)loc];
        return V3{P.x*1000.0, P.y*1000.0, P.z*1000.0};
    }

    // SHELLS
    const int s    = loc / nPtsShell;
    const int k    = loc - s*nPtsShell;
    const int iLon = k % nLon;
    const int jLat = k / nLon;

    double lon = d_deg * iLon;
    double lat = -90.0 + d_deg * jLat;
    if (lat > 90.0) lat = 90.0;

    const double alt_km  = prm.output.shellAlt_km[(size_t)s];
    const double r_m     = (_RADIUS_(_EARTH_) + alt_km*1000.0);
    const double lonRad  = lon * M_PI / 180.0;
    const double latRad  = lat * M_PI / 180.0;
    const double cl      = std::cos(latRad);

    const V3 xCart{r_m*cl*std::cos(lonRad), r_m*cl*std::sin(lonRad), r_m*std::sin(latRad)};

    // Pure DIPOLE: use the lon/lat-derived position directly (no frame rotation).
    if (EarthUtil::ToUpper(prm.field.model) == "DIPOLE") return xCart;

    // External models: rotate from Earth-fixed (ITRF93) to GSM via SPICE.
#ifndef _NO_SPICE_CALLS_
    static std::string cachedEpoch;
    static SpiceDouble rot[3][3];
    if (cachedEpoch != prm.field.epoch) {
        cachedEpoch = prm.field.epoch;
        SpiceDouble et;
        str2et_c(prm.field.epoch.c_str(), &et);
        pxform_c("ITRF93", "GSM", et, rot);
    }
    SpiceDouble xGEO[3] = {xCart.x, xCart.y, xCart.z};
    SpiceDouble xGSM[3];
    mxv_c(rot, xGEO, xGSM);
    return V3{xGSM[0], xGSM[1], xGSM[2]};
#else
    // SPICE unavailable: return un-rotated position (acceptable for DIPOLE; warn
    // for external models, but do not hard-fail to allow non-SPICE builds).
    return xCart;
#endif
}

//======================================================================================
// SECTION 13 — TECPLOT OUTPUT WRITERS
//======================================================================================
//
// The Störmer analytic reference is computed via
// Earth::GridlessMode::StormerVerticalCoeff_GV (declared inline in
// CutoffRigidityGridless.h) and Earth::GridlessMode::Dipole::gParams.m_hat
// (from DipoleInterface.h) — the same functions and state used by the gridless
// solver's WriteTecplotShells_DipoleAnalyticCompare.  This guarantees that the
// analytic reference values in cutoff_3d_shells.dat are numerically identical
// to those that would appear in the gridless dipole-compare output file.

static void WriteTecplot3DPoints(const EarthUtil::AmpsParam& prm,
                                 const std::vector<double>& Rc,
                                 const std::vector<double>& Emin,
                                 int nLoc) {
    FILE* f = std::fopen("cutoff_3d_points.dat", "w");
    if (!f) throw std::runtime_error("Cannot write cutoff_3d_points.dat");

    const bool hasTraj = (EarthUtil::ToUpper(prm.output.mode) == "TRAJECTORY" &&
                          !prm.output.trajectories.empty());

    std::fprintf(f, "TITLE=\"Mode3D Cutoff Rigidity (POINTS/TRAJECTORY)\"\n");
    std::fprintf(f, "VARIABLES=");
    if (hasTraj) std::fprintf(f, "\"TimeUTC\" ");
    std::fprintf(f, "\"id\",\"x_km\",\"y_km\",\"z_km\","
                    "\"lon_deg\",\"lat_deg\",\"Rc_GV\",\"Emin_MeV\"\n");
    std::fprintf(f, "ZONE T=\"%s\" I=%d F=POINT\n",
                 hasTraj ? "trajectory" : "points", nLoc);

    for (int i = 0; i < nLoc; ++i) {
        const auto& P = prm.output.points[(size_t)i];
        const double lon = std::atan2(P.y, P.x) * 180.0 / M_PI;
        const double lat = std::atan2(P.z, std::sqrt(P.x*P.x + P.y*P.y)) * 180.0 / M_PI;

        if (hasTraj) {
            const auto& s = prm.output.trajectories[0].samples[(size_t)i];
            std::fprintf(f, "\"%s\" ", s.timeUTC.c_str());
        }
        std::fprintf(f, "%d %e %e %e %e %e %e %e\n",
                     i, P.x, P.y, P.z, lon, lat, Rc[(size_t)i], Emin[(size_t)i]);
    }
    std::fclose(f);
}

static void WriteTecplot3DShells(const EarthUtil::AmpsParam& prm,
                                 const std::vector<std::vector<double>>& RcShell,
                                 const std::vector<std::vector<double>>& EminShell,
                                 int nLon, int nLat, double d_deg) {
    //--------------------------------------------------------------------------
    // Analytic Störmer reference columns
    //
    // Enabled when FIELD_MODEL = DIPOLE  AND  CUTOFF_SAMPLING = VERTICAL.
    //
    // The vertical cutoff in a centred-dipole field has a closed-form Störmer
    // solution:  Rc(λ,r) = R₀ · cos⁴(λ) / (r/Rₑ)²
    // where λ is the magnetic latitude measured from the dipole axis m̂.
    //
    // When CUTOFF_SAMPLING is VERTICAL the numerical solver traces the same
    // radially inward direction as the Störmer vertical — the comparison is
    // therefore physically meaningful.  For ISOTROPIC sampling the numerical
    // result is the minimum over all directions (generally less than the
    // vertical Störmer value), making the comparison misleading; in that case
    // the extra columns are omitted.
    //
    // Extra VARIABLES written per data row when both conditions are met:
    //   Rc_Stormer_GV  — analytic Störmer vertical cutoff [GV]
    //   rel_err        — (Rc_numeric − Rc_Stormer) / Rc_Stormer
    //--------------------------------------------------------------------------
    const bool isDipole      = (EarthUtil::ToUpper(prm.field.model)    == "DIPOLE");
    const bool isVertical    = (EarthUtil::ToUpper(prm.cutoff.sampling) == "VERTICAL");
    const bool writeAnalytic = isDipole && isVertical;

    // Pre-compute Störmer coefficient R₀ and the dipole axis unit vector m̂.
    //
    // We use Earth::GridlessMode::StormerVerticalCoeff_GV and
    // Earth::GridlessMode::Dipole::gParams.m_hat — exactly the same state and
    // function as WriteTecplotShells_DipoleAnalyticCompare in the gridless solver.
    //
    // SetTiltDeg is called explicitly here because the Mode3D flow does not go
    // through the gridless cFieldEvaluator constructor that normally sets gParams.
    double R0_GV = 0.0;
    double mx = 0.0, my = 0.0, mz = 1.0;
    if (writeAnalytic) {
        Earth::GridlessMode::Dipole::SetTiltDeg(prm.field.dipoleTilt_deg);
        mx = Earth::GridlessMode::Dipole::gParams.m_hat[0];
        my = Earth::GridlessMode::Dipole::gParams.m_hat[1];
        mz = Earth::GridlessMode::Dipole::gParams.m_hat[2];
        R0_GV = Earth::GridlessMode::StormerVerticalCoeff_GV(
                    prm.field.dipoleMoment_Me, _EARTH__RADIUS_);
    }

    FILE* f = std::fopen("cutoff_3d_shells.dat", "w");
    if (!f) throw std::runtime_error("Cannot write cutoff_3d_shells.dat");

    std::fprintf(f, "TITLE=\"Mode3D Cutoff Rigidity (SHELLS)\"\n");
    if (writeAnalytic) {
        std::fprintf(f, "VARIABLES=\"lon_deg\",\"lat_deg\","
                        "\"Rc_GV\",\"Emin_MeV\","
                        "\"Rc_Stormer_GV\",\"rel_err\"\n");
    } else {
        std::fprintf(f, "VARIABLES=\"lon_deg\",\"lat_deg\","
                        "\"Rc_GV\",\"Emin_MeV\"\n");
    }

    for (size_t s = 0; s < prm.output.shellAlt_km.size(); ++s) {
        const double alt_km = prm.output.shellAlt_km[s];
        std::fprintf(f, "ZONE T=\"alt_km=%g\" I=%d J=%d F=POINT\n",
                     alt_km, nLon, nLat);

        // r/Re for this shell — must use the same Earth-radius constant as the
        // AMR domain geometry so the Störmer denominator is consistent.
        const double r_m = _RADIUS_(_EARTH_) + alt_km * 1000.0;
        const double rRe = r_m / _EARTH__RADIUS_;

        const int nPts = nLon * nLat;
        for (int k = 0; k < nPts; ++k) {
            const int    iLon = k % nLon;
            const int    jLat = k / nLon;
            double lat = -90.0 + d_deg * jLat;
            if (lat > 90.0) lat = 90.0;
            const double lon  = d_deg * iLon;
            const double rc   = RcShell[s][(size_t)k];
            const double emin = EminShell[s][(size_t)k];

            if (writeAnalytic) {
                // Unit position vector r̂ from (lon, lat).
                const double lonRad = lon * M_PI / 180.0;
                const double latRad = lat * M_PI / 180.0;
                const double clat   = std::cos(latRad);
                const double rhatx  = clat * std::cos(lonRad);
                const double rhaty  = clat * std::sin(lonRad);
                const double rhatz  = std::sin(latRad);

                // Magnetic latitude λ:  sin(λ) = m̂ · r̂
                const double sinLam = mx*rhatx + my*rhaty + mz*rhatz;
                const double cosLam = std::sqrt(std::max(0.0, 1.0 - sinLam*sinLam));

                // Störmer vertical cutoff — identical formula to
                // WriteTecplotShells_DipoleAnalyticCompare in the gridless solver:
                //   Rc_vert = R0 * dipoleMoment_Me * cos^4(lambda) / (r/Re)^2
                const double Rc_S = R0_GV
                                  * prm.field.dipoleMoment_Me
                                  * std::pow(cosLam, 4) / (rRe * rRe);

                // Relative error: defined only when both values are positive.
                // At the magnetic poles cosLam → 0 and Rc_S → 0; avoid 0/0.
                const double relErr = (Rc_S > 0.0 && rc > 0.0)
                                     ? (rc - Rc_S) / Rc_S
                                     : 0.0;

                std::fprintf(f, "%e %e %e %e %e %e\n",
                             lon, lat, rc, emin, Rc_S, relErr);
            } else {
                std::fprintf(f, "%e %e %e %e\n", lon, lat, rc, emin);
            }
        }
    }
    std::fclose(f);
}

} // end anonymous namespace


//======================================================================================
// SECTION 14 — PUBLIC ENTRY POINT: RunCutoffRigidity
//======================================================================================

namespace Earth {
namespace Mode3D {

int RunCutoffRigidity(const EarthUtil::AmpsParam& prm) {

    //==================================================================================
    // 14.1 — MPI initialisation guard
    //==================================================================================
    // MPI is always compiled in.  Guard MPI_Init so the function is safe to call
    // both under mpirun (MPI already initialised) and in single-process debug runs.
    //==================================================================================

    int mpiInited = 0;
    MPI_Initialized(&mpiInited);
    bool ownedMPI = false;
    if (!mpiInited) {
        int argc_dummy = 0; char** argv_dummy = nullptr;
        MPI_Init(&argc_dummy, &argv_dummy);
        ownedMPI = true;
    }

    int mpiRank = 0, mpiSize = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);

    //==================================================================================
    // 14.2 — Species constants
    //==================================================================================

    const double qabs = std::fabs(prm.species.charge_e * ElectronCharge);
    const double q_C  = prm.species.charge_e * ElectronCharge;
    const double m0   = prm.species.mass_amu * _AMU_;

    if (qabs == 0.0)
        throw std::runtime_error("Mode3D cutoff: particle charge must be non-zero.");

    //==================================================================================
    // 14.3 — Rigidity bracket from energy limits
    //==================================================================================

    const double pMin  = MomentumFromKineticEnergy_MeV(prm.cutoff.eMin_MeV, m0);
    const double pMax  = MomentumFromKineticEnergy_MeV(prm.cutoff.eMax_MeV, m0);
    const double Rmin  = RigidityFromMomentum_GV(pMin, qabs);
    const double Rmax  = RigidityFromMomentum_GV(pMax, qabs);

    if (!(Rmax > Rmin) || !(Rmax > 0.0))
        throw std::runtime_error("Mode3D cutoff: invalid energy bracket "
                                 "(CUTOFF_EMIN/EMAX in #CUTOFF_RIGIDITY).");

    //==================================================================================
    // 14.4 — Domain geometry from the AMPS AMR mesh
    //==================================================================================
    //
    // ParsedDomainMin/Max are set in SI metres by ApplyParsedDomain() in Mode3D.cpp
    // before RunCutoffRigidity is called.  The inner sphere radius equals the Earth
    // radius (loss sphere centred at origin, same as in amps_init_mesh).
    //==================================================================================

    if (!ParsedDomainActive)
        throw std::runtime_error("Mode3D cutoff: domain has not been configured "
                                 "(ParsedDomainActive is false). "
                                 "Ensure Mode3D::Run() calls ApplyParsedDomain before "
                                 "RunCutoffRigidity.");

    DomainBox3D box;
    box.xMin   = ParsedDomainMin[0];  box.xMax   = ParsedDomainMax[0];
    box.yMin   = ParsedDomainMin[1];  box.yMax   = ParsedDomainMax[1];
    box.zMin   = ParsedDomainMin[2];  box.zMax   = ParsedDomainMax[2];
    box.rInner = _EARTH__RADIUS_;     // loss sphere = Earth surface

    //==================================================================================
    // 14.5 — Direction grid (Fibonacci sphere)
    //==================================================================================

    const std::string samplingMode = EarthUtil::ToUpper(prm.cutoff.sampling);
    const bool samplingVertical  = (samplingMode == "VERTICAL");
    const bool samplingIsotropic = (samplingMode == "ISOTROPIC" || samplingMode.empty());

    if (!samplingVertical && !samplingIsotropic)
        throw std::runtime_error("Mode3D cutoff: unsupported CUTOFF_SAMPLING '" +
                                 prm.cutoff.sampling + "'. Use VERTICAL or ISOTROPIC.");

    const int nCutoffDirs = prm.cutoff.maxParticlesPerPoint;
    std::vector<V3> dirs;
    if (!samplingVertical)
        dirs = BuildFibonacciDirs(nCutoffDirs);

    //==================================================================================
    // 14.6 — Location list geometry
    //==================================================================================

    const std::string outputMode = EarthUtil::ToUpper(prm.output.mode);
    const bool isPoints = (outputMode == "POINTS" || outputMode == "TRAJECTORY");
    const bool isShells = (outputMode == "SHELLS");

    if (!isPoints && !isShells)
        throw std::runtime_error("Mode3D cutoff: unsupported OUTPUT_MODE '" +
                                 prm.output.mode + "'. Use POINTS, TRAJECTORY, or SHELLS.");

    const double d_deg    = isShells ? prm.output.shellRes_deg : 0.0;
    const int    nShells  = isShells ? static_cast<int>(prm.output.shellAlt_km.size()) : 0;
    const int    nLon     = isShells ? static_cast<int>(std::floor(360.0/d_deg + 0.5)) : 0;
    const int    nLat     = isShells ? static_cast<int>(std::floor(180.0/d_deg + 0.5)) + 1 : 0;
    const int    nPtsShell = isShells ? nLon * nLat : 0;

    const int nLoc = isPoints
                   ? static_cast<int>(prm.output.points.size())
                   : nShells * nPtsShell;

    if (nLoc == 0)
        throw std::runtime_error("Mode3D cutoff: no observation points/shells defined.");

    //==================================================================================
    // 14.7 — Run summary (rank 0 only)
    //==================================================================================

    if (mpiRank == 0) {
        std::cout
            << "======== Mode3D Cutoff Rigidity ========\n"
            << "Field model    : " << prm.field.model    << "\n"
            << "Epoch          : " << prm.field.epoch    << "\n"
            << "Species        : " << prm.species.name
            << " (q=" << prm.species.charge_e << " e"
            << ", m=" << prm.species.mass_amu << " amu)\n"
            << "Rigidity range : [" << Rmin << ", " << Rmax << "] GV\n"
            << "Sampling       : " << (samplingVertical ? "VERTICAL" : "ISOTROPIC") << "\n";

        if (!samplingVertical)
            std::cout << "N_directions   : " << dirs.size()
                      << " (Fibonacci sphere)\n";

        std::cout
            << "N_locations    : " << nLoc  << "\n"
            << "MPI ranks      : " << mpiSize << " (replicated domain, static partition)\n";

#ifdef _OPENMP
        std::cout << "OpenMP threads : " << omp_get_max_threads()
                  << " per rank\n";
#else
        std::cout << "OpenMP threads : 1 (built without OpenMP)\n";
#endif

        std::cout
            << "Domain [m]     : x[" << box.xMin << "," << box.xMax << "] "
            << "y[" << box.yMin << "," << box.yMax << "] "
            << "z[" << box.zMin << "," << box.zMax << "]\n"
            << "Inner sphere r : " << box.rInner << " m\n"
            << "=========================================\n";
        std::cout.flush();
    }

    //==================================================================================
    // 14.8 — Static MPI point distribution
    //==================================================================================
    //
    // Requirement 1: Each MPI process is INDEPENDENT (no inter-rank communication
    //                during computation).
    // Requirement 3: Parallelization of individual points with OpenMP (within rank).
    //
    // Simple static block decomposition:
    //   rank k owns points [locStart, locEnd).
    //   Last rank picks up any remainder so all points are covered.
    //==================================================================================

    const int nPerRank  = nLoc / mpiSize;
    const int locStart  = mpiRank * nPerRank;
    const int locEnd    = (mpiRank == mpiSize - 1) ? nLoc : locStart + nPerRank;
    const int nLocal    = locEnd - locStart;

    // Local result arrays (indexed 0..nLocal-1)
    std::vector<double> rcLocal(nLocal, -1.0);
    std::vector<double> eminLocal(nLocal, -1.0);

    //==================================================================================
    // 14.9 — OpenMP parallel loop over local points
    //==================================================================================
    //
    // Requirement 3: Parallelization of calculation for individual points with OpenMP.
    //
    // Each OpenMP thread owns a PRIVATE cMode3DMeshFieldEval instance.
    // The evaluator reads from the frozen AMR mesh (no shared mutable state),
    // so no mutex or synchronisation is needed between threads.
    //
    // The mutable lastNode_ cache inside cMode3DMeshFieldEval is per-instance
    // (and therefore per-thread), so findTreeNode() hint updates are naturally
    // thread-safe without any explicit protection.
    //==================================================================================

    #pragma omp parallel
    {
        // Thread-private mesh field evaluator.
        // Constructed once per thread; lastNode_ starts as nullptr and is
        // updated on every GetB_T call as the particle moves through the mesh.
        cMode3DMeshFieldEval threadField;

        #pragma omp for schedule(dynamic, 1)
        for (int localIdx = 0; localIdx < nLocal; ++localIdx) {
            const int globalIdx = locStart + localIdx;

            // Reconstruct the observation position for this point.
            const V3 x0_m = LocationToX0m(prm, globalIdx, nLon, nLat, d_deg, nPtsShell);

            // Compute the minimum cutoff rigidity over all sampled directions.
            // The mesh field is the same for all points (initialized once by
            // InitMeshFields in Mode3D::Run before this parallel section).
            const double rc = ComputeCutoffAtPoint_GV(
                prm, threadField, x0_m, dirs, samplingVertical,
                q_C, m0, box, Rmin, Rmax);

            // Store result — each localIdx is unique to this thread in the
            // dynamic schedule, so the critical section is a no-op in practice
            // but guards against any future collapse/vectorisation reordering.
            #pragma omp critical
            {
                rcLocal[(size_t)localIdx]   = rc;
                if (rc > 0.0) {
                    const double pCut = MomentumFromRigidity_GV(rc, qabs);
                    eminLocal[(size_t)localIdx] = KineticEnergyFromMomentum_MeV(pCut, m0);
                }
            }

        } // omp for
    } // omp parallel

    //==================================================================================
    // 14.10 — MPI gather: collect results from all ranks onto rank 0
    //==================================================================================
    //
    // Each rank contributes nLocal doubles for Rc and Emin.
    // MPI_Gatherv handles the variable block sizes that arise from the remainder
    // clause in the static partition (last rank may have a larger slice).
    //==================================================================================

    // Build displacement and count arrays for MPI_Gatherv on rank 0
    std::vector<int> recvCounts(mpiSize, 0);
    std::vector<int> displs(mpiSize, 0);

    // Each rank announces its nLocal to rank 0
    MPI_Gather(&nLocal, 1, MPI_INT, recvCounts.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);

    std::vector<double> rcAll, eminAll;
    if (mpiRank == 0) {
        int offset = 0;
        for (int r = 0; r < mpiSize; ++r) {
            displs[r] = offset;
            offset   += recvCounts[r];
        }
        rcAll.resize((size_t)nLoc, -1.0);
        eminAll.resize((size_t)nLoc, -1.0);
    }

    MPI_Gatherv(rcLocal.data(),   nLocal, MPI_DOUBLE,
                rcAll.data(),   recvCounts.data(), displs.data(), MPI_DOUBLE,
                0, MPI_COMM_WORLD);

    MPI_Gatherv(eminLocal.data(), nLocal, MPI_DOUBLE,
                eminAll.data(), recvCounts.data(), displs.data(), MPI_DOUBLE,
                0, MPI_COMM_WORLD);

    MPI_Barrier(MPI_COMM_WORLD);

    //==================================================================================
    // 14.11 — Output (rank 0 only)
    //==================================================================================

    if (mpiRank == 0) {
        if (isPoints) {
            // Console summary
            for (int i = 0; i < nLoc; ++i) {
                const auto& P = prm.output.points[(size_t)i];
                std::cout << "Point " << i
                          << " (" << P.x << "," << P.y << "," << P.z << " km)"
                          << " -> Rc=" << rcAll[(size_t)i] << " GV"
                          << "  Emin=" << eminAll[(size_t)i] << " MeV\n";
            }

            WriteTecplot3DPoints(prm, rcAll, eminAll, nLoc);
            std::cout << "Wrote: cutoff_3d_points.dat\n";

        } else {
            // SHELLS: reshape flat rcAll into per-shell arrays
            std::vector<std::vector<double>> RcShell(nShells),
                                              EminShell(nShells);
            for (int s = 0; s < nShells; ++s) {
                RcShell[s].assign((size_t)nPtsShell, -1.0);
                EminShell[s].assign((size_t)nPtsShell, -1.0);
                for (int k = 0; k < nPtsShell; ++k) {
                    const int locId = s*nPtsShell + k;
                    RcShell[s][(size_t)k]   = rcAll[(size_t)locId];
                    EminShell[s][(size_t)k] = eminAll[(size_t)locId];
                }
                std::cout << "Shell alt=" << prm.output.shellAlt_km[(size_t)s]
                          << " km done.\n";
            }

            WriteTecplot3DShells(prm, RcShell, EminShell, nLon, nLat, d_deg);
            std::cout << "Wrote: cutoff_3d_shells.dat\n";
        }
        std::cout.flush();
    }

    //==================================================================================
    // 14.12 — MPI finalise (only if we initialised MPI in this function)
    //==================================================================================

    if (ownedMPI) {
        int fin = 0; MPI_Finalized(&fin);
        if (!fin) MPI_Finalize();
    }

    return 0;
}

} // namespace Mode3D
} // namespace Earth
