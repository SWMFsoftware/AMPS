//======================================================================================
// CutoffRigidityMode3D.cpp
//======================================================================================
//
// Standalone Mode3D cutoff-rigidity solver.
//
// This file computes cutoff rigidity by backward-tracing test particles through the
// Mode3D magnetic-field mesh.  It is used when the Earth model is run in standalone
// ``-mode 3d`` rather than through the SWMF coupler.  In this mode the background
// magnetic field is already sampled and stored on the AMPS AMR mesh; trajectory
// integration therefore reuses the same mesh-field evaluator that is used for
// Mode3D density/flux backtracking.
//
// High-level algorithm
// --------------------
// For each requested observation location on a POINTS, SHELLS, or TRAJECTORY output
// set, the code performs the following operations:
//
//   1. Convert the output location to GSM Cartesian coordinates in metres.
//
//   2. Choose the arrival direction(s).
//        * VERTICAL sampling uses the local outward radial direction as the arrival
//          direction.
//        * ISOTROPIC sampling uses a Fibonacci sphere of directions.
//
//   3. For each direction, launch the mathematically reversed particle in the
//      opposite direction and integrate it through the mesh field.
//        * Escape through the outer Mode3D box means the rigidity/direction is
//          ALLOWED.
//        * Contact with the inner Earth-loss sphere means FORBIDDEN.
//        * Reaching the time, step, or path-length limit without escape is treated
//          conservatively as FORBIDDEN.
//
//   4. Determine the cutoff rigidity for that direction.  The default method is the
//      penumbra-safe UPPER_SCAN algorithm described in detail below.  The older
//      endpoint binary search remains available for regression tests.
//
//   5. For an isotropic point cutoff, return the minimum directional cutoff over
//      the sampled directions.  When DIRECTIONAL_MAP=T, also save each directional
//      cutoff on the sky grid.
//
// Penumbra-safe upper cutoff algorithm
// ------------------------------------
// A simple binary search is valid only if the trajectory classifier
//
//      allowed(R) = TraceAllowed3D(R)
//
// is monotonic: forbidden below the cutoff and allowed above it.  Dipole shell tests
// at high latitude show that this is not always true.  Small allowed pockets can
// appear below the final transition (the classic penumbra problem).  In that case
// the legacy endpoint logic
//
//      if allowed(Rmin) and allowed(Rmax), return Rmin
//
// is wrong, because Rmin may lie inside a low-rigidity allowed pocket while higher
// rigidities are still forbidden.
//
// The default UPPER_SCAN algorithm therefore defines and returns the upper cutoff:
// the lowest rigidity above the highest forbidden interval sampled in the requested
// search bracket.  Operationally, for one point/direction it does this:
//
//   1. Build a rigidity grid from Rmin to Rmax.
//      Log spacing is used for positive Rmin because SEP/geospace cutoff searches
//      cover multiple decades in rigidity.
//
//   2. Evaluate allowed(R_i) at every grid point.
//
//   3. If the highest sampled point Rmax is forbidden, report ``no cutoff in the
//      requested bracket`` using the existing -1 return convention.
//
//   4. Starting at Rmax, scan downward until the highest forbidden grid point is
//      found.  The next-higher point is necessarily allowed because the downward
//      scan starts inside the allowed top branch.
//
//   5. Refine only this final forbidden/allowed transition by local bisection.
//      Lower-rigidity allowed pockets are intentionally ignored because they do not
//      define the upper cutoff.
//
// This algorithm is more expensive than endpoint binary search because it performs
// a coarse scan before bisection, but it is robust for the high-latitude dipole
// shell cases where allowed(R) is non-monotonic.  The number of scan points is
// controlled by CUTOFF_UPPER_SCAN_N; if that key is absent, CUTOFF_NENERGY is reused
// so existing input files still provide a meaningful scan density.
//
// See CutoffRigidityMode3D.h for the broader design rationale, parallelisation
// strategy, and input/output contract.  This file documents implementation details.
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
// Mode3D now uses a penumbra-safe upper-cutoff search by default.  The earlier
// endpoint-only binary search assumed that TraceAllowed3D(R) was monotonic
// (forbidden below Rc, allowed above Rc).  Dipole shell tests at high latitude
// showed low-rigidity allowed pockets below the Störmer cutoff; the endpoint test
// then returned Rmin as soon as TraceAllowed3D(Rmin)==true.
//
// The default UPPER_SCAN algorithm instead:
//   1. Samples TraceAllowed3D(R) on a log-spaced rigidity grid.
//   2. Scans downward from Rmax and finds the highest forbidden sampled rigidity.
//   3. Refines the transition between that forbidden sample and the next-higher
//      allowed sample with bisection.
//
// This returns the upper cutoff: the lowest rigidity above which the sampled
// trajectory family is allowed.  It deliberately ignores allowed pockets below
// the final forbidden/allowed transition.  The legacy endpoint binary method is
// still available with CUTOFF_SEARCH_ALGORITHM BINARY.
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
// Mesh/field distribution:
//   Mode3D::Run now initializes AMPS with the normal distributed MPI domain
//   decomposition.  Before this solver starts, Mode3D::GlobalMagneticField has
//   already allocated missing nonlocal leaf blocks on every rank and populated them
//   with a replicated read-only cell-centered B snapshot.  Therefore trajectory
//   tracing can remain local-memory-only even though the mesh was not created with
//   independentDomainMode=true.
//
// Inter-rank scheduling:
//   The default scheduler is now DYNAMIC.  Rank/main threads atomically fetch chunks
//   of global observation locations from an MPI one-sided counter.  This allows a rank
//   that finishes easy trajectories quickly to immediately take more work instead of
//   waiting for a rank that received a cluster of near-cutoff/high-latitude points.
//   BLOCK_CYCLIC and STATIC are retained as deterministic fallback/debug schedulers.
//
// Two-level parallelism:
//   The MPI level assigns chunks to ranks.  Inside each rank the selected backend
//   (OPENMP, THREADS, or SERIAL) schedules locations within the chunk.  In THREADS
//   mode the intra-rank scheduler is also dynamic and uses std::atomic<int>.  MPI is
//   called only by the rank/main thread, never by worker threads.
//
// Communication:
//   The global-B materialization is completed before entering RunCutoffRigidity().
//   During trajectory tracing there is no field communication.  At the end, each rank
//   holds global-indexed result arrays with -1 sentinels for locations it did not
//   compute; MPI_MAX reductions collect Rc, Emin, and optional directional-map values
//   on rank 0 in the original global output order.
//
//======================================================================================

#include "CutoffRigidityMode3D.h"
#include "Mode3D.h"
#include "ElectricField.h"
#include "Mode3DParallel.h"

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
#include <atomic>
#include <cstdlib>
#include <memory>
#include <mutex>
#include <sstream>
#include <stdexcept>
#include <iostream>
#include <limits>
#include <thread>

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
// SECTION 1.5 — MODE3D CUTOFF PARALLEL BACKEND SELECTION
//======================================================================================
//
// The cutoff calculation is embarrassingly parallel over observation locations: every
// location owns an independent set of rigidity searches and trajectory backtraces, and
// the only communication after the field snapshot is prepared is the final MPI gather.
// Older versions used only OpenMP inside each MPI rank.  The density-backtracking path
// later added a direct std::thread backend because some platforms showed problems with
// OpenMP placement and/or OpenMP interaction with shared AMPS/SWMF state.  The cutoff
// path now uses the same user settings so one input deck controls both particle-free
// backward products:
//
//   #NUMERICAL
//   DENSITY_PARALLEL  OPENMP | THREADS | SERIAL
//   DENSITY_THREADS   <N>
//
// Naming note:
//   The keywords keep the historical DENSITY_* names because they were introduced first
//   for density/flux products.  In Mode3D they should now be read as the generic
//   particle-free backward-product parallel backend.  They affect both density/flux and
//   cutoff calculations.
//
// Backend meanings in this file:
//   OPENMP : use the existing OpenMP parallel-for over local observation locations;
//            DENSITY_THREADS, if positive, is passed to omp_set_num_threads().
//   THREADS: use direct std::thread workers over local observation locations; before
//            creating workers, widen this MPI rank's CPU affinity with
//            PIC::Parallel::SetWideAffinityForScheduler(), matching the density path.
//   SERIAL : no intra-rank threading; MPI decomposition remains active.
//
// Thread-safety note:
//   Each OpenMP thread or std::thread worker constructs its own cMode3DMeshFieldEval.
//   The evaluator reads from the replicated, read-only AMR magnetic-field snapshot and
//   owns its own lastNode_ cache, so no worker shares interpolation/cache state.
//======================================================================================

using CutoffParallelBackend_ = Earth::Mode3D::ParallelBackend;

static const char* CutoffParallelBackendName_(CutoffParallelBackend_ backend) {
    return Earth::Mode3D::ParallelBackendName(backend);
}

static CutoffParallelBackend_ ResolveCutoffParallelBackend_(const EarthUtil::AmpsParam& prm) {
    return Earth::Mode3D::ResolveParallelBackend(prm,"Mode3D cutoff");
}

static int ResolveCutoffThreadCount_(const EarthUtil::AmpsParam& prm,
                                     CutoffParallelBackend_ backend) {
    return Earth::Mode3D::ResolveParallelThreadCount(prm,backend);
}

static void ApplyWideAffinityForDirectCutoffThreadsOnce_(CutoffParallelBackend_ backend,
                                                         int cutoffThreadCount) {
    Earth::Mode3D::ApplyWideAffinityForDirectThreadsOnce(backend,cutoffThreadCount,
                                                         "Mode3D cutoff");
}

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
// SECTION 2.5 — DIRECTIONAL-MAP FRAME HELPERS
//======================================================================================
//
// The optional DIRECTIONAL_MAP diagnostic asks for Rc as a function of arrival
// direction on a regular lon/lat sky grid.  The gridless cutoff solver already
// established the convention used here, and the standalone mesh-based Mode3D
// implementation intentionally reuses the same convention so the two paths can
// be compared file-by-file:
//
//   * the sky-map labels are global spherical longitude/latitude in SM;
//   * particle tracing still happens in GSM, because the mesh-field snapshot is
//     stored in the same GSM Cartesian coordinates used by the rest of Mode3D;
//   * each map-cell direction is therefore constructed in SM and then rotated to
//     GSM with a SPICE pxform("SM","GSM",epoch) matrix when SPICE is available;
//   * if the executable was built without SPICE, or the frame transform is not
//     available at run time, the code falls back to identity.  The calculation is
//     still performed, but the map is effectively labeled in GSM.
//
// These helpers are deliberately local to this file.  They are small copies of
// the gridless helpers, avoiding a new linkage dependency while keeping the math
// and comments synchronized with the already-tested gridless implementation.
//======================================================================================

struct Mat3 {
    double a[3][3];
};

static inline Mat3 Identity3() {
    Mat3 R{};
    R.a[0][0]=1.0; R.a[0][1]=0.0; R.a[0][2]=0.0;
    R.a[1][0]=0.0; R.a[1][1]=1.0; R.a[1][2]=0.0;
    R.a[2][0]=0.0; R.a[2][1]=0.0; R.a[2][2]=1.0;
    return R;
}

static inline V3 Apply(const Mat3& R, const V3& v) {
    return V3{
        R.a[0][0]*v.x + R.a[0][1]*v.y + R.a[0][2]*v.z,
        R.a[1][0]*v.x + R.a[1][1]*v.y + R.a[1][2]*v.z,
        R.a[2][0]*v.x + R.a[2][1]*v.y + R.a[2][2]*v.z
    };
}

static inline Mat3 GetSpiceRotationOrIdentity3D(const char* fromFrame,
                                                const char* toFrame,
                                                const std::string& epoch,
                                                bool& ok) {
    ok = false;
    Mat3 R = Identity3();

#ifndef _NO_SPICE_CALLS_
    try {
        // SPICE expects ephemeris time.  The epoch string is already parsed and
        // used elsewhere in Mode3D, so the same user-provided UTC string is used
        // here to make the map direction labels consistent with the field epoch.
        SpiceDouble et = 0.0;
        str2et_c(epoch.c_str(), &et);

        SpiceDouble m[3][3];
        pxform_c(fromFrame, toFrame, et, m);

        for (int i=0;i<3;i++) {
            for (int j=0;j<3;j++) R.a[i][j] = m[i][j];
        }
        ok = true;
    }
    catch (...) {
        // CSPICE is a C library, but some local wrappers/builds can still route
        // errors through exceptions.  A failed transform must not abort the cutoff
        // calculation; it only means the directional labels cannot be rotated from
        // SM to GSM, so identity is used and ok=false is reported in the banner.
        ok = false;
    }
#else
    (void)fromFrame;
    (void)toFrame;
    (void)epoch;
#endif

    return R;
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
// cMode3DMeshFieldEval interpolates the magnetic field that InitMeshFields() stored in the
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
// Field data layout for each stencil cell (mirrors InitMeshFields write):
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
    explicit cMode3DMeshFieldEval(const EarthUtil::AmpsParam& prm)
      : prm_(prm), model_(EarthUtil::ToUpper(prm.field.model)), lastNode_(nullptr) {}

    // Evaluate B [Tesla] at position x_m [m]. By default this uses AMR-aware
    // interpolation from the cell-centered values written by InitMeshFields().
    // When requested from the CLI, it instead calls the same background-field
    // evaluator used by InitMeshFields() to prepopulate those cell centers.
    void GetB_T(const V3& x_m, V3& B_T) const override {
        double xArr[3] = {x_m.x, x_m.y, x_m.z};

        if (prm_.mode3d.forceAnalyticMagneticField) {
            double B[3] = {0.0, 0.0, 0.0};

            // Evaluate with the same function used by InitMeshFields().  This path
            // may touch shared model state (Tsyganenko/TA Fortran common blocks, or
            // the dipole helper state), so serialize it for both OpenMP and direct
            // std::thread density workers.
            static std::mutex analyticFieldMutex;
            {
                std::lock_guard<std::mutex> lock(analyticFieldMutex);
                Earth::Mode3D::EvaluateBackgroundMagneticFieldSI(B, xArr, prm_);
            }

            B_T.x = B[0];
            B_T.y = B[1];
            B_T.z = B[2];
            return;
        }

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

        if (_PIC_COUPLER_MODE_ == _PIC_COUPLER_MODE__SWMF_) {
            // Do not call PIC::CPLR::InitInterpolationStencil() +
            // PIC::CPLR::GetBackgroundMagneticField() here.  In AMPS hybrid builds,
            // GetBackgroundMagneticField(B) reads the global
            // CellCentered::StencilTable[omp_get_thread_num()]. That is correct inside
            // an OpenMP team, but plain std::thread workers all see OpenMP thread id 0
            // and would race on the same stencil entry.  Build a stack-local stencil
            // instead and read the already materialized SWMF cell-centered B field
            // directly from PIC::CPLR::SWMF::MagneticFieldOffset.  This exactly matches
            // the cell accessor in AMPS (SWMF::GetBackgroundMagneticField(cell) is a
            // memcpy from that offset) while keeping both OpenMP and std::thread paths
            // thread-safe.
            if (PIC::CPLR::SWMF::MagneticFieldOffset < 0) {
                B_T.x = B_T.y = B_T.z = 0.0;
                return;
            }

            PIC::InterpolationRoutines::CellCentered::cStencil Stencil;
            PIC::InterpolationRoutines::CellCentered::Linear::InitStencil(xArr,node,Stencil);

            if (Stencil.Length <= 0) {
                B_T.x = B_T.y = B_T.z = 0.0;
                return;
            }

            double B[3] = {0.0,0.0,0.0};
            double BCell[3];

            for (int iStencil = 0; iStencil < Stencil.Length; ++iStencil) {
                PIC::Mesh::cDataCenterNode* center = Stencil.cell[iStencil];
                if (center == nullptr) continue;

                const char* data =
                    PIC::CPLR::SWMF::MagneticFieldOffset + center->GetAssociatedDataBufferPointer();
                std::memcpy(BCell,data,3*sizeof(double));

                const double w = Stencil.Weight[iStencil];
                B[0] += w*BCell[0];
                B[1] += w*BCell[1];
                B[2] += w*BCell[2];
            }

            B_T.x = B[0];
            B_T.y = B[1];
            B_T.z = B[2];
            return;
        }

        if (!PIC::CPLR::DATAFILE::Offset::MagneticField.active) {
            B_T.x = B_T.y = B_T.z = 0.0;
            return;
        }

        // Build AMPS' cell-centred linear interpolation stencil at the particle
        // position.  This routine is AMR-aware: when the particle is near a
        // refinement boundary it constructs a multi-block stencil and blends
        // fine/coarse stencils as needed, instead of assuming a uniform local
        // trilinear stencil.
        PIC::InterpolationRoutines::CellCentered::cStencil Stencil;
        PIC::InterpolationRoutines::CellCentered::Linear::InitStencil(xArr, node, Stencil);

        if (Stencil.Length <= 0) {
            B_T.x = B_T.y = B_T.z = 0.0;
            return;
        }

        // Interpolate B from the cell-centred values written by InitMeshFields().
        // DATAFILE::GetBackgroundData applies the same offset chain used when the
        // data were stored:
        //   CenterNodeAssociatedDataOffsetBegin + CurrDataFileOffset
        //   + Offset::MagneticField.RelativeOffset.
        double B[3] = {0.0, 0.0, 0.0};
        double BCell[3];

        for (int iStencil = 0; iStencil < Stencil.Length; ++iStencil) {
            PIC::Mesh::cDataCenterNode* center = Stencil.cell[iStencil];
            if (center == nullptr) continue;

            PIC::CPLR::DATAFILE::GetBackgroundData(
                BCell, 3,
                PIC::CPLR::DATAFILE::Offset::MagneticField.RelativeOffset,
                center);

            const double w = Stencil.Weight[iStencil];
            B[0] += w * BCell[0];
            B[1] += w * BCell[1];
            B[2] += w * BCell[2];
        }

        B_T.x = B[0];
        B_T.y = B[1];
        B_T.z = B[2];
    }

private:
    const EarthUtil::AmpsParam& prm_;
    std::string model_;

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
                            double maxTime_s = -1.0,
                            Earth::GridlessMode::TrajectoryExitState* exitState = nullptr) {
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
        if (!InsideBox3D(x, box)) {
            // Particle has escaped the outer Mode3D box: trajectory is ALLOWED.
            //
            // Density/flux in ANISOTROPIC boundary mode needs the asymptotic exit
            // state so it can evaluate the boundary pitch-angle distribution and
            // optional spatial modulation in exactly the same manner as the gridless
            // density solver.  The cutoff-rigidity calculation does not request this
            // information, so the extra work is skipped unless exitState!=nullptr.
            if (exitState) {
                exitState->x_exit_m[0] = x.x;
                exitState->x_exit_m[1] = x.y;
                exitState->x_exit_m[2] = x.z;

                // Recover velocity direction from relativistic momentum p=gamma*m*v.
                // The common scalar factor 1/(gamma*m) cancels during normalization,
                // but computing v explicitly keeps the intent clear and matches the
                // gridless TraceAllowedSharedEx implementation.
                const double p2n = dot(p,p);
                const double mc  = m0_kg*SpeedOfLight;
                const double gExit = std::sqrt(1.0 + p2n/(mc*mc));
                const V3 vExit = v3unit(mul(1.0/(gExit*m0_kg), p));
                exitState->v_exit_unit[0] = vExit.x;
                exitState->v_exit_unit[1] = vExit.y;
                exitState->v_exit_unit[2] = vExit.z;

                // Pitch-angle cosine at the exit point.  The mesh field is valid up
                // to the last in-domain cell; when the point has just crossed outside
                // the box GetB_T may return zero.  In that case use cosAlpha=0, the
                // same conservative fallback used by gridless when |B| is unavailable.
                V3 Bexit; field.GetB_T(x, Bexit);
                const double Bnorm = v3norm(Bexit);
                exitState->cosAlpha = (Bnorm > 0.0)
                    ? dot(vExit, mul(1.0/Bnorm, Bexit))
                    : 0.0;
            }
            return true;  // ALLOWED
        }

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
// SECTION 8.1 — OPTIONAL DETAILED EXIT-CLASSIFIER DIAGNOSTIC
//======================================================================================
//
// The production TraceAllowed3D() routine intentionally returns only a boolean:
//   true  -> the reversed trajectory escaped through the outer Mode3D domain box;
//   false -> the trajectory was lost or did not escape before a numerical limit.
//
// That minimal interface is ideal for cutoff and density calculations, but it hides
// information that is essential when validating the trajectory classifier. In pure
// dipole tests, especially at high latitude where penumbral allowed/forbidden bands
// appear, we need to verify that an "allowed" result really means that the trajectory
// crossed the outer boundary, not that it stopped because of a timeout, step limit,
// invalid time step, or inner-sphere hit.
//
// The detailed diagnostic below repeats the same integration loop as TraceAllowed3D(),
// but records the terminal state, the reason for termination, and several consistency
// quantities:
//
//   * raw terminal point x_terminal_m: the position after the step that left the box;
//   * last in-domain point x_last_inside_m: the previous position;
//   * linearly reconstructed boundary crossing x_boundary_m, face name, and overshoot;
//   * start/end rigidity conservation error (E=0 in the dipole test);
//   * canonical angular momentum about the dipole symmetry axis for FIELD_MODEL=DIPOLE.
//
// The canonical invariant is
//
//     P_axis = [ r x (p + q A) ] . m_hat,
//
// where A = (mu0/4pi) (m x r)/r^3 is the centered-dipole vector potential and m_hat
// is the dipole-axis unit vector. For a static aligned or tilted centered dipole this
// quantity should be conserved because the field is axisymmetric about m_hat. The
// diagnostic writes the relative change only as a check; the classifier decision itself
// is still based solely on the same boundary logic used by TraceAllowed3D().
//
// This debug path is deliberately not used in the production cutoff loop. It is a
// single-point/single-scan validation tool enabled by CUTOFF_DEBUG_EXIT_TRACE.
//======================================================================================

enum class TraceExitReason3D {
    OUTER_BOX,
    INNER_SPHERE_PRE,
    INNER_SPHERE_STEP,
    TIME_LIMIT,
    STEP_LIMIT,
    DISTANCE_LIMIT,
    INVALID_DT,
    UNKNOWN
};

static const char* TraceExitReasonName3D_(TraceExitReason3D r) {
    switch (r) {
        case TraceExitReason3D::OUTER_BOX:         return "OUTER_BOX";
        case TraceExitReason3D::INNER_SPHERE_PRE:  return "INNER_SPHERE_PRE";
        case TraceExitReason3D::INNER_SPHERE_STEP: return "INNER_SPHERE_STEP";
        case TraceExitReason3D::TIME_LIMIT:        return "TIME_LIMIT";
        case TraceExitReason3D::STEP_LIMIT:        return "STEP_LIMIT";
        case TraceExitReason3D::DISTANCE_LIMIT:    return "DISTANCE_LIMIT";
        case TraceExitReason3D::INVALID_DT:        return "INVALID_DT";
        default:                                   return "UNKNOWN";
    }
}

static const char* MoverTypeName3D_(MoverType m) {
    switch (m) {
        case MoverType::BORIS:  return "BORIS";
        case MoverType::RK2:    return "RK2";
        case MoverType::RK4:    return "RK4";
        case MoverType::RK6:    return "RK6";
        case MoverType::GC2:    return "GC2";
        case MoverType::GC4:    return "GC4";
        case MoverType::GC6:    return "GC6";
        case MoverType::HYBRID: return "HYBRID";
        default:                return "UNKNOWN";
    }
}

struct TraceDetailedResult3D {
    bool allowed{false};
    TraceExitReason3D reason{TraceExitReason3D::UNKNOWN};

    double R_in_GV{0.0};
    double R_terminal_GV{0.0};
    double relRigidityError{0.0};

    double t_s{0.0};
    int nSteps{0};
    double path_m{0.0};

    V3 x0_m{0.0,0.0,0.0};
    V3 xLastInside_m{0.0,0.0,0.0};
    V3 xTerminal_m{0.0,0.0,0.0};
    V3 xBoundary_m{0.0,0.0,0.0};

    double boundaryAlpha{-1.0};
    double boundaryOvershoot_m{0.0};
    double terminalBoxMargin_m{0.0};
    std::string boundaryFace{"NONE"};

    double pAxis0{std::numeric_limits<double>::quiet_NaN()};
    double pAxisTerminal{std::numeric_limits<double>::quiet_NaN()};
    double relPAxisError{std::numeric_limits<double>::quiet_NaN()};
};

static double BoxMinimumSignedMargin3D_(const V3& x, const DomainBox3D& box) {
    // Positive value: point is inside all six faces; magnitude is the distance to the
    // nearest face. Negative value: point is outside at least one face; magnitude is
    // the penetration distance through the most violated face.
    return std::min({x.x - box.xMin, box.xMax - x.x,
                     x.y - box.yMin, box.yMax - x.y,
                     x.z - box.zMin, box.zMax - x.z});
}

static void FindFirstBoxCrossing3D_(const V3& xin,
                                    const V3& xout,
                                    const DomainBox3D& box,
                                    V3& xcross,
                                    double& alpha,
                                    std::string& face) {
    // Reconstruct where the last integration segment crossed the rectangular Mode3D
    // domain. The mover returns only the post-step point, which can lie a finite
    // distance outside the box. Linear interpolation is sufficient for this diagnostic
    // because the reconstructed point is not fed back into the production classifier.
    const V3 dx = sub(xout,xin);
    double best = 2.0;
    std::string bestFace = "NONE";

    auto consider = [&](double a, const char* name) {
        if (a >= 0.0 && a <= 1.0 && a < best) {
            best = a;
            bestFace = name;
        }
    };

    if (dx.x > 0.0) consider((box.xMax - xin.x)/dx.x,"XMAX");
    if (dx.x < 0.0) consider((box.xMin - xin.x)/dx.x,"XMIN");
    if (dx.y > 0.0) consider((box.yMax - xin.y)/dx.y,"YMAX");
    if (dx.y < 0.0) consider((box.yMin - xin.y)/dx.y,"YMIN");
    if (dx.z > 0.0) consider((box.zMax - xin.z)/dx.z,"ZMAX");
    if (dx.z < 0.0) consider((box.zMin - xin.z)/dx.z,"ZMIN");

    if (best <= 1.0) {
        alpha = best;
        xcross = add(xin,mul(best,dx));
        face = bestFace;
    }
    else {
        alpha = -1.0;
        xcross = xout;
        face = "UNKNOWN";
    }
}

static V3 DipoleVectorPotential_Tm_(const V3& x_m) {
    // Centered-dipole vector potential in SI units:
    //   A(r) = (mu0/4pi) (m x r) / r^3.
    // The result has units T*m = Wb/m, and q*A has the same units as mechanical
    // momentum. This expression is used only by the invariant diagnostic.
    const double r2 = dot(x_m,x_m);
    if (!(r2 > 0.0)) return {0.0,0.0,0.0};
    const double r = std::sqrt(r2);
    const double r3 = r2*r;

    const double M = Earth::GridlessMode::Dipole::gParams.momentScale_Me *
                     Earth::GridlessMode::Dipole::M_E_Am2;
    const V3 mvec{M*Earth::GridlessMode::Dipole::gParams.m_hat[0],
                  M*Earth::GridlessMode::Dipole::gParams.m_hat[1],
                  M*Earth::GridlessMode::Dipole::gParams.m_hat[2]};
    constexpr double mu0_over_4pi = 1.0e-7;
    return mul(mu0_over_4pi/r3,cross(mvec,x_m));
}

static double DipoleCanonicalPAxis_(const EarthUtil::AmpsParam& prm,
                                    const V3& x_m,
                                    const V3& p,
                                    double q_C) {
    // Canonical angular momentum about the dipole symmetry axis. This is conserved for
    // the ideal centered dipole because the vector potential and Hamiltonian are
    // invariant under rotations around the dipole axis.
    if (EarthUtil::ToUpper(prm.field.model) != "DIPOLE") {
        return std::numeric_limits<double>::quiet_NaN();
    }

    const V3 A = DipoleVectorPotential_Tm_(x_m);
    const V3 canonicalMomentum = add(p,mul(q_C,A));
    const V3 axis{Earth::GridlessMode::Dipole::gParams.m_hat[0],
                  Earth::GridlessMode::Dipole::gParams.m_hat[1],
                  Earth::GridlessMode::Dipole::gParams.m_hat[2]};
    return dot(cross(x_m,canonicalMomentum),axis);
}

static TraceDetailedResult3D TraceDetailed3D_(const EarthUtil::AmpsParam& prm,
                                              cMode3DMeshFieldEval& field,
                                              const V3& x0_m,
                                              const V3& v0_unit,
                                              double R_GV,
                                              double q_C,
                                              double m0_kg,
                                              const DomainBox3D& box,
                                              double maxTime_s = -1.0) {
    TraceDetailedResult3D out;
    out.R_in_GV = R_GV;
    out.x0_m = x0_m;
    out.xLastInside_m = x0_m;
    out.xTerminal_m = x0_m;
    out.xBoundary_m = x0_m;

    const double qabs = std::fabs(q_C);
    const double pMag = MomentumFromRigidity_GV(R_GV,qabs);
    V3 p = mul(pMag,v0_unit);
    const V3 p0 = p;
    V3 x = x0_m;
    V3 xLastInside = x0_m;

    const double maxTraceTime_s =
        (maxTime_s > 0.0)
            ? maxTime_s
            : ((prm.cutoff.maxTrajTime_s > 0.0)
               ? prm.cutoff.maxTrajTime_s
               : prm.numerics.maxTraceTime_s);

    const double maxDist_m =
        (prm.numerics.maxTraceDistance_Re > 0.0)
            ? prm.numerics.maxTraceDistance_Re * _EARTH__RADIUS_
            : -1.0;

    ResetHybridTrajectoryContext(x0_m, box.rInner);

    double tTrace = 0.0;
    int nSteps = 0;
    double sDist = 0.0;

    auto finalize = [&](TraceExitReason3D reason, bool allowed) {
        out.allowed = allowed;
        out.reason = reason;
        out.t_s = tTrace;
        out.nSteps = nSteps;
        out.path_m = sDist;
        out.xLastInside_m = xLastInside;
        out.xTerminal_m = x;
        out.R_terminal_GV = RigidityFromMomentum_GV(std::sqrt(std::max(0.0,dot(p,p))),qabs);
        out.relRigidityError = (std::fabs(out.R_in_GV) > 0.0)
            ? (out.R_terminal_GV - out.R_in_GV)/out.R_in_GV
            : 0.0;
        out.terminalBoxMargin_m = BoxMinimumSignedMargin3D_(x,box);

        if (reason == TraceExitReason3D::OUTER_BOX) {
            FindFirstBoxCrossing3D_(xLastInside,x,box,out.xBoundary_m,
                                    out.boundaryAlpha,out.boundaryFace);
            out.boundaryOvershoot_m = v3norm(sub(x,out.xBoundary_m));
        }
        else {
            out.xBoundary_m = x;
            out.boundaryAlpha = -1.0;
            out.boundaryFace = "NONE";
            out.boundaryOvershoot_m = 0.0;
        }

        out.pAxis0 = DipoleCanonicalPAxis_(prm,x0_m,p0,q_C);
        out.pAxisTerminal = DipoleCanonicalPAxis_(prm,x,p,q_C);
        if (std::isfinite(out.pAxis0) && std::fabs(out.pAxis0) > 0.0) {
            out.relPAxisError = (out.pAxisTerminal - out.pAxis0)/out.pAxis0;
        }
    };

    while (nSteps < prm.numerics.maxSteps &&
           tTrace < maxTraceTime_s &&
           (maxDist_m <= 0.0 || sDist < maxDist_m)) {

        if (LostInnerSphere3D(x, box.rInner)) {
            finalize(TraceExitReason3D::INNER_SPHERE_PRE,false);
            return out;
        }
        if (!InsideBox3D(x, box)) {
            finalize(TraceExitReason3D::OUTER_BOX,true);
            return out;
        }

        const double timeRemaining = maxTraceTime_s - tTrace;

        bool useGuidingCenterForThisStep = false;
        const MoverType mover = GetDefaultMoverType();
        if (mover == MoverType::GC2 || mover == MoverType::GC4 || mover == MoverType::GC6) {
            useGuidingCenterForThisStep = true;
        } else if (mover == MoverType::HYBRID) {
            useGuidingCenterForThisStep = HybridPrepareStepUseGuidingCenter(x, p, q_C, field);
        }

        const double dt = SelectDt3D(prm, field, x, p, q_C, m0_kg, box,
                                      timeRemaining, useGuidingCenterForThisStep);
        if (!(dt > 0.0) || !std::isfinite(dt)) {
            finalize(TraceExitReason3D::INVALID_DT,false);
            return out;
        }

        const V3 xPrev = x;
        xLastInside = xPrev;
        if (!StepParticleChecked(gDefaultMover, x, p, q_C, m0_kg, dt, field, box.rInner)) {
            finalize(TraceExitReason3D::INNER_SPHERE_STEP,false);
            return out;
        }

        sDist  += v3norm(sub(x, xPrev));
        tTrace += dt;
        ++nSteps;
    }

    if (nSteps >= prm.numerics.maxSteps) {
        finalize(TraceExitReason3D::STEP_LIMIT,false);
    }
    else if (maxDist_m > 0.0 && sDist >= maxDist_m) {
        finalize(TraceExitReason3D::DISTANCE_LIMIT,false);
    }
    else if (tTrace >= maxTraceTime_s) {
        finalize(TraceExitReason3D::TIME_LIMIT,false);
    }
    else {
        finalize(TraceExitReason3D::UNKNOWN,false);
    }

    return out;
}

//======================================================================================
// SECTION 9 — CUTOFF SEARCH FOR ONE DIRECTION
//======================================================================================
//
// CutoffForDir_GV is the scalar rigidity solver for exactly one observation point and
// exactly one arrival direction.  It is called many times by ComputeCutoffAtPoint_GV:
// once for VERTICAL sampling, or once per sky direction for ISOTROPIC sampling.
//
// Important conventions:
//
//   * dir_unit is the physical ARRIVAL direction at the observation point.
//     Backward tracing launches the reversed particle along -dir_unit.
//
//   * Rmin_GV and Rmax_GV are already converted from the input energy range to
//     magnetic rigidity in GV.  This routine does not know about kinetic energy,
//     species energy per nucleon, or output units.
//
//   * Return value > 0 means a cutoff was found in GV.
//     Return value < 0 keeps the legacy convention meaning: no allowed top branch
//     was found within [Rmin_GV, Rmax_GV], i.e. the cutoff is above Rmax or the
//     trajectory classifier never escaped at the sampled upper bound.
//
//   * The UPPER_SCAN search returns the upper cutoff, not the first allowed pocket.
//     This is the physically useful quantity for cutoff maps because it defines the
//     rigidity above which particles are continuously allowed, ignoring lower-rigidity
//     penumbral windows.
//
// Two algorithms are available:
//
//   1. ENDPOINT_BINARY / BINARY / LEGACY_BINARY
//      Kept only for comparison with older results.  It tests the two endpoints and
//      then assumes monotonic allowed(R).  It can return Rmin incorrectly when Rmin
//      lies inside a low-rigidity allowed pocket.
//
//   2. UPPER_SCAN (default)
//      First samples the full bracket, locates the highest forbidden sample, then
//      bisects only that final forbidden/allowed transition.  This handles
//      non-monotonic allowed(R) sequences observed in the dipole shell regression
//      test near |latitude|≈60 degrees.
//
// Per-location upper-bound optimisation:
//
//   The isotropic point cutoff is the minimum over directions.  Once a direction has
//   produced Rc_found for a point, later directions cannot lower the point cutoff if
//   their directional cutoff is above Rc_found.  Therefore ComputeCutoffAtPoint_GV can
//   pass Rmax=min(original_Rmax, Rc_found) to subsequent directions.  This makes the
//   search cheaper while preserving the minimum-direction cutoff definition.
//======================================================================================

static double CutoffForDirEndpointBinary_GV(const EarthUtil::AmpsParam& prm,
                                             cMode3DMeshFieldEval& field,
                                             const V3& x0_m,
                                             const V3& dir_unit,   // ARRIVAL direction (backtraced as -dir)
                                             double q_C, double m0_kg,
                                             const DomainBox3D& box,
                                             double Rmin_GV, double Rmax_GV) {
    // Legacy monotonic solver.
    //
    // This routine is intentionally simple and intentionally preserved: it is useful
    // when comparing new Mode3D cutoff maps against archived results that were produced
    // with the original endpoint-binary method.  It should NOT be used as the default
    // for production dipole or magnetospheric shell maps because it assumes a monotonic
    // classifier.
    //
    // Backtracing convention:
    //   dir_unit is the particle arrival direction at the observation location.  A
    //   backward trajectory is equivalent to launching the reversed particle along
    //   -dir_unit with the same charge sign convention used by the shared mover code.
    const V3 v0 = mul(-1.0, dir_unit);

    // The bisection tolerance follows the historical cutoff implementation: an
    // absolute 1 MV tolerance, or a much smaller relative tolerance for very large
    // brackets.  The absolute tolerance dominates in normal geospace use and is more
    // meaningful than oversolving the trajectory classifier, whose uncertainty is set
    // by step-size and boundary-contact tests.
    const double tolAbs = 1.0e-3;
    const double tolRel = 1.0e-6 * std::max(std::fabs(Rmin_GV), std::fabs(Rmax_GV));
    const double tol    = std::max(tolAbs, tolRel);
    if (Rmax_GV < Rmin_GV) return -1.0;

    // Endpoint tests.  This is exactly where the old method can fail: if Rmin is in a
    // low-rigidity allowed pocket and Rmax is also allowed, alo&&ahi causes immediate
    // return of Rmin even if forbidden rigidities exist between them.
    const bool alo = TraceAllowed3D(prm, field, x0_m, v0, Rmin_GV, q_C, m0_kg, box);
    const bool ahi = TraceAllowed3D(prm, field, x0_m, v0, Rmax_GV, q_C, m0_kg, box);

    if (alo && ahi) return Rmin_GV;  // valid only if allowed(R) is monotonic
    if (!ahi)       return -1.0;     // no allowed upper branch inside the bracket

    // Standard binary search on the assumed final transition.  The invariant is
    // lo=forbidden and hi=allowed.
    double lo = Rmin_GV, hi = Rmax_GV;
    while ((hi - lo) > tol) {
        const double mid = 0.5*(lo + hi);
        const bool a = TraceAllowed3D(prm, field, x0_m, v0, mid, q_C, m0_kg, box);
        if (a) hi = mid; else lo = mid;
    }
    return hi;
}

static int CutoffUpperScanPointCount_(const EarthUtil::AmpsParam& prm) {
    // Decide how many coarse samples are used by UPPER_SCAN before the final local
    // bisection.  This number controls two things:
    //
    //   * robustness: a finer grid is less likely to skip a narrow forbidden band near
    //     the true upper cutoff;
    //   * cost: every grid point requires a complete trajectory trace.
    //
    // CUTOFF_UPPER_SCAN_N is the explicit new control.  If it is not provided, reuse
    // CUTOFF_NENERGY.  Reusing the existing parameter is intentional: older input files
    // that already requested a fine cutoff-energy grid automatically get a similarly
    // fine penumbra scan without another required keyword.
    if (prm.cutoff.upperScanN > 0) return std::max(2,prm.cutoff.upperScanN);
    return std::max(8,prm.cutoff.nEnergy);
}

static std::vector<double> BuildCutoffSearchGrid_GV_(double Rmin_GV,
                                                     double Rmax_GV,
                                                     int nScan) {
    // Build the coarse rigidity grid used by UPPER_SCAN.  The grid is strictly local to
    // one point/direction; it is not an output energy grid and it is not shared among
    // MPI ranks.
    //
    // Positive-rigidity searches use logarithmic spacing.  This is important because
    // geospace cutoffs may span from tens of MV to many GV.  A linear grid over that
    // interval would waste nearly all samples at high rigidity and could miss a
    // low-rigidity upper cutoff such as the 0.16 GV transition in the dipole |lat|=60
    // regression case.
    //
    // The fallback linear branch exists only for completeness if a future input permits
    // Rmin<=0; current physical rigidity brackets should always be positive.
    std::vector<double> grid;
    if (!(Rmax_GV >= Rmin_GV) || !(Rmax_GV > 0.0)) return grid;

    nScan = std::max(2,nScan);
    grid.reserve((size_t)nScan);

    if (Rmin_GV > 0.0) {
        const double lmin = std::log(Rmin_GV);
        const double lmax = std::log(Rmax_GV);
        for (int i=0;i<nScan;i++) {
            const double a = (nScan == 1) ? 0.0 : (double)i/(double)(nScan-1);
            grid.push_back(std::exp((1.0-a)*lmin + a*lmax));
        }
    }
    else {
        for (int i=0;i<nScan;i++) {
            const double a = (nScan == 1) ? 0.0 : (double)i/(double)(nScan-1);
            grid.push_back((1.0-a)*Rmin_GV + a*Rmax_GV);
        }
    }

    // Force exact endpoints.  This avoids tiny roundoff shifts introduced by exp/log
    // and keeps diagnostic output and endpoint comparisons reproducible.
    grid.front() = Rmin_GV;
    grid.back()  = Rmax_GV;
    return grid;
}

static double RefineForbiddenAllowedTransition_GV_(const EarthUtil::AmpsParam& prm,
                                                   cMode3DMeshFieldEval& field,
                                                   const V3& x0_m,
                                                   const V3& v0,
                                                   double q_C, double m0_kg,
                                                   const DomainBox3D& box,
                                                   double Rforbid_GV,
                                                   double Rallow_GV) {
    // Refine one already-bracketed forbidden/allowed transition.
    //
    // Preconditions established by CutoffForDirUpperScan_GV:
    //   TraceAllowed3D(Rforbid_GV) == false
    //   TraceAllowed3D(Rallow_GV)  == true
    //   Rforbid_GV < Rallow_GV
    //
    // This local bisection does NOT attempt to solve the entire possibly non-monotonic
    // allowed(R) function.  It is applied only to the highest forbidden sample found by
    // the downward scan.  At that point we are refining the final top-branch boundary,
    // which is the definition of the upper cutoff used here.
    const double tolAbs = 1.0e-3;
    const double tolRel = 1.0e-6 * std::max(std::fabs(Rforbid_GV),std::fabs(Rallow_GV));
    const double tol    = std::max(tolAbs,tolRel);

    double lo = Rforbid_GV; // invariant: lo is forbidden
    double hi = Rallow_GV;  // invariant: hi is allowed

    while ((hi - lo) > tol) {
        const double mid = 0.5*(lo + hi);
        const bool allowed = TraceAllowed3D(prm,field,x0_m,v0,mid,q_C,m0_kg,box);
        if (allowed) hi = mid;
        else         lo = mid;
    }

    // Return the allowed side of the bracket: the smallest allowed rigidity resolved
    // to the requested tolerance.
    return hi;
}

static double CutoffForDirUpperScan_GV(const EarthUtil::AmpsParam& prm,
                                       cMode3DMeshFieldEval& field,
                                       const V3& x0_m,
                                       const V3& dir_unit,   // ARRIVAL direction (backtraced as -dir)
                                       double q_C, double m0_kg,
                                       const DomainBox3D& box,
                                       double Rmin_GV, double Rmax_GV) {
    // Penumbra-safe upper-cutoff search.
    //
    // This is the production/default solver.  It avoids the bad endpoint-binary
    // behavior observed in the dipole shell test at |lat|=60 deg, where the lowest
    // sampled rigidity was allowed, several intermediate rigidities were forbidden,
    // and all high rigidities were allowed.  The correct map value is the final upper
    // transition near the Störmer cutoff, not the low-rigidity allowed pocket.

    // Backtracing convention: launch the reversed particle in direction -dir_unit.
    const V3 v0 = mul(-1.0, dir_unit);

    if (Rmax_GV < Rmin_GV) return -1.0;

    // Coarse rigidity scan.  Each grid point is expensive because it is a complete
    // trajectory integration, but the scan is what makes the method robust to
    // non-monotonic allowed/forbidden sequences.
    const int nScan = CutoffUpperScanPointCount_(prm);
    const std::vector<double> grid = BuildCutoffSearchGrid_GV_(Rmin_GV,Rmax_GV,nScan);
    if (grid.size() < 2) return -1.0;

    // Store the classifier result as bytes instead of bool.  std::vector<bool> uses a
    // packed proxy specialization that is inconvenient for debugging and can create
    // surprising behavior when references are taken.  A byte array is explicit and
    // still tiny compared with the trajectory cost.
    std::vector<unsigned char> allowed(grid.size(),0);
    for (size_t i=0;i<grid.size();++i) {
        allowed[i] = TraceAllowed3D(prm,field,x0_m,v0,grid[i],q_C,m0_kg,box) ? 1 : 0;
    }

    // The top branch must be allowed for an upper cutoff to exist within the requested
    // search interval.  If Rmax is still forbidden, the final transition is above the
    // bracket, so return the existing negative sentinel.  The output layer writes that
    // as missing/no-cutoff according to its normal convention.
    if (!allowed.back()) return -1.0;

    // Walk downward from Rmax.  The first forbidden point encountered is the highest
    // forbidden sampled rigidity.  Since allowed.back()==true and we are scanning from
    // high to low, the next-higher sample grid[i+1] is on the final allowed branch.
    // Bisection between grid[i] and grid[i+1] returns the upper cutoff.
    for (int i=(int)grid.size()-2;i>=0;--i) {
        if (!allowed[(size_t)i]) {
            return RefineForbiddenAllowedTransition_GV_(prm,field,x0_m,v0,q_C,m0_kg,
                                                        box,grid[(size_t)i],grid[(size_t)i+1]);
        }
    }

    // No forbidden sample was found anywhere in the bracket.  That means the entire
    // sampled interval is allowed, so the upper cutoff lies below Rmin or cannot be
    // resolved with the requested lower bound.  Return Rmin to keep the output finite
    // and to match the old convention for an everywhere-allowed bracket.
    return Rmin_GV;
}

static double CutoffForDir_GV(const EarthUtil::AmpsParam& prm,
                              cMode3DMeshFieldEval& field,
                              const V3& x0_m,
                              const V3& dir_unit,   // ARRIVAL direction (backtraced as -dir)
                              double q_C, double m0_kg,
                              const DomainBox3D& box,
                              double Rmin_GV, double Rmax_GV) {
    // Public dispatcher used by the rest of this file.  Keeping the algorithm choice
    // centralized is important because cutoff maps, directional maps, and the debug
    // scan should all use the same selected production algorithm unless they explicitly
    // call the legacy/debug helpers.
    const std::string alg = EarthUtil::ToUpper(prm.cutoff.searchAlgorithm);

    if (alg=="BINARY" || alg=="ENDPOINT_BINARY" || alg=="LEGACY_BINARY") {
        return CutoffForDirEndpointBinary_GV(prm,field,x0_m,dir_unit,
                                             q_C,m0_kg,box,Rmin_GV,Rmax_GV);
    }

    // Default: penumbra-safe upper cutoff.  Unknown strings are validated by the input
    // parser, so reaching this branch means UPPER_SCAN or one of its accepted aliases.
    return CutoffForDirUpperScan_GV(prm,field,x0_m,dir_unit,
                                    q_C,m0_kg,box,Rmin_GV,Rmax_GV);
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
// Thread-safety note:
//   The cutoff std::thread backend may call LocationToX0m() concurrently from several
//   worker threads.  The SPICE branch below uses a small cached rotation matrix for the
//   shell grid.  Protecting that cache is inexpensive because the location transform is
//   performed once per observation location, while the expensive work is the many
//   trajectory traces launched from that location.
//
// Logic is identical to the LocationToX0m lambda in CutoffRigidityGridless.cpp.
// Replicated here to avoid a cross-module dependency on a local lambda.
static std::mutex gLocationToX0mSpiceMutex_;

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
    {
        std::lock_guard<std::mutex> lock(gLocationToX0mSpiceMutex_);

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
    }
#else
    // SPICE unavailable: return un-rotated position (acceptable for DIPOLE; warn
    // for external models, but do not hard-fail to allow non-SPICE builds).
    return xCart;
#endif
}

//======================================================================================
// SECTION 12.5 — OUTPUT FILE NAMING
//======================================================================================
//
// Standalone 3-D cutoff runs historically write fixed file names
// (cutoff_3d_points.dat, cutoff_3d_shells.dat).  A live SWMF-coupled run can call
// amps_time_step() repeatedly for different MHD snapshots.  Reusing the fixed names
// would overwrite the previous cutoff products.
//
// The optional suffix below is therefore appended immediately before ".dat".  The
// default empty suffix preserves the standalone file contract exactly.
//======================================================================================

static std::string gCutoffOutputFileSuffix;

static std::string CutoffOutputFileName(const char* stem) {
    return std::string(stem) + gCutoffOutputFileSuffix + ".dat";
}

//======================================================================================
// SECTION 12.6 — DIRECTIONAL-MAP CONFIGURATION AND CELL GEOMETRY
//======================================================================================
//
// This structure collects every quantity needed to compute and write the optional
// directional cutoff sky maps.  Keeping it in one object avoids passing several
// loosely-related integers/doubles through the OpenMP worker lambda and through
// the Tecplot writer.
//
// enabled:
//   True only when DIRECTIONAL_MAP=T.  The parser already stores that keyword in
//   prm.cutoff.directionalMap; this block is the first place where standalone
//   Mode3D acts on it.
//
// nLon/nLat/nCells:
//   Regular lon/lat grid dimensions.  Longitude covers [0,360) and latitude
//   covers [-90,+90] inclusively, matching the gridless directional-map writer.
//
// R_label2gsm:
//   Direction-label frame to tracing-frame rotation.  The label frame is SM by
//   convention; the tracing frame is GSM because the mesh-based field evaluator
//   and Mode3D trajectory geometry are both in GSM Cartesian coordinates.
//======================================================================================

struct DirectionalMapConfig3D {
    bool enabled{false};
    double lonRes_deg{0.0};
    double latRes_deg{0.0};
    int nLon{0};
    int nLat{0};
    int nCells{0};
    bool spiceOk{false};
    Mat3 R_label2gsm{Identity3()};
};

static DirectionalMapConfig3D ConfigureDirectionalMap3D(const EarthUtil::AmpsParam& prm) {
    DirectionalMapConfig3D cfg;
    cfg.enabled = prm.cutoff.directionalMap;

    if (!cfg.enabled) return cfg;

    cfg.lonRes_deg = prm.cutoff.dirMapLonRes_deg;
    cfg.latRes_deg = prm.cutoff.dirMapLatRes_deg;

    if (!(cfg.lonRes_deg > 0.0) || !(cfg.latRes_deg > 0.0)) {
        throw std::runtime_error(
            "Mode3D cutoff: DIRMAP_LON_RES and DIRMAP_LAT_RES must be > 0 "
            "when DIRECTIONAL_MAP=T.");
    }

    // Reuse the gridless convention exactly.  Rounding allows common values such
    // as 7.5, 10, 15, 30 degrees to produce the expected integer grid size even
    // if the decimal representation is not exact.
    cfg.nLon = static_cast<int>(std::floor(360.0/cfg.lonRes_deg + 0.5));
    cfg.nLat = static_cast<int>(std::floor(180.0/cfg.latRes_deg + 0.5)) + 1;
    cfg.nCells = cfg.nLon * cfg.nLat;

    if (cfg.nLon <= 0 || cfg.nLat <= 0 || cfg.nCells <= 0) {
        throw std::runtime_error(
            "Mode3D cutoff: invalid directional-map grid size derived from "
            "DIRMAP_LON_RES/DIRMAP_LAT_RES.");
    }

    // Direction labels are in SM; tracing directions must be in GSM.  When SPICE
    // is unavailable, GetSpiceRotationOrIdentity3D returns identity and sets
    // spiceOk=false.  Rank 0 prints a warning later, after MPI rank is known.
    cfg.R_label2gsm = GetSpiceRotationOrIdentity3D(
        "SM", "GSM", prm.field.epoch, cfg.spiceOk);

    return cfg;
}

static V3 DirectionalMapCellDirectionGSM3D(const DirectionalMapConfig3D& cfg,
                                           int cellId) {
    // Cell indexing convention is shared with gridless:
    //   cellId = iLon + nLon*jLat
    // where lon increases fastest in the Tecplot POINT zone.
    const int iLon = cellId % cfg.nLon;
    const int jLat = cellId / cfg.nLon;

    double lon_deg = cfg.lonRes_deg * iLon;
    double lat_deg = -90.0 + cfg.latRes_deg * jLat;
    if (lat_deg > 90.0) lat_deg = 90.0;

    // Build the unit arrival direction in the labeling frame (SM).  This is a
    // global spherical direction, not a geographic position and not a local ENU
    // direction at the observation point.
    const double lon = lon_deg * M_PI/180.0;
    const double lat = lat_deg * M_PI/180.0;
    const double cl  = std::cos(lat);
    const V3 dirLabel{cl*std::cos(lon), cl*std::sin(lon), std::sin(lat)};

    // Rotate to GSM for tracing.  v3unit protects against a malformed transform
    // or roundoff at the poles; the input vector is already unit length.
    return v3unit(Apply(cfg.R_label2gsm, dirLabel));
}

static std::string DirectionalMapOutputFileName3D(int locId) {
    char stem[128];
    std::snprintf(stem, sizeof(stem), "cutoff_3d_dir_map_loc_%06d", locId);
    return CutoffOutputFileName(stem);
}

//======================================================================================
// SECTION 13 — TECPLOT OUTPUT WRITERS
//======================================================================================
//
// The regular Mode3D output files keep the same compact schema as the gridless
// production outputs:
//   cutoff_3d_points.dat
//   cutoff_3d_shells.dat
//
// For DIPOLE nightly tests, additional comparison files are written, matching
// the gridless-path convention:
//   cutoff_3d_points_dipole_compare.dat
//   cutoff_3d_shells_dipole_compare.dat
//
// These comparison files report the numerically calculated cutoff beside the
// analytic Størmer vertical cutoff for the same point/shell node.  For a clean
// apples-to-apples benchmark the run should use CUTOFF_SAMPLING VERTICAL; if an
// isotropic run is used, Rc_num_GV is the minimum over sampled directions and is
// not expected to equal the vertical Størmer value.
//======================================================================================

static double StormerVerticalCutoff_GV(const EarthUtil::AmpsParam& prm,
                                       const V3& x_m) {
    const double r_m = v3norm(x_m);
    if (!(r_m > 0.0)) return 0.0;

    const double rhatx = x_m.x / r_m;
    const double rhaty = x_m.y / r_m;
    const double rhatz = x_m.z / r_m;

    const double mx = Earth::GridlessMode::Dipole::gParams.m_hat[0];
    const double my = Earth::GridlessMode::Dipole::gParams.m_hat[1];
    const double mz = Earth::GridlessMode::Dipole::gParams.m_hat[2];

    const double sinLam = mx*rhatx + my*rhaty + mz*rhatz;
    const double cosLam = std::sqrt(std::max(0.0, 1.0 - sinLam*sinLam));
    const double rRe    = r_m / _EARTH__RADIUS_;

    // Same formula used by the gridless comparison writers:
    //   Rc_vert = R0 * dipoleMoment_Me * cos^4(lambda) / (r/Re)^2
    // Note: StormerVerticalCoeff_GV currently follows the gridless convention,
    // where the caller also multiplies by dipoleMoment_Me.
    const double R0_GV = Earth::GridlessMode::StormerVerticalCoeff_GV(
                             prm.field.dipoleMoment_Me, _EARTH__RADIUS_);

    return R0_GV * prm.field.dipoleMoment_Me
                 * std::pow(cosLam, 4) / (rRe * rRe);
}

static void EnsureDipoleAnalyticState3D(const EarthUtil::AmpsParam& prm) {
    // Mode3D does not necessarily construct the gridless cFieldEvaluator that
    // initializes Dipole::gParams.  Set both the magnitude scale and the axis before
    // any analytic Störmer comparison or dipole invariant diagnostic is written.
    Earth::GridlessMode::Dipole::SetMomentScale(prm.field.dipoleMoment_Me);
    Earth::GridlessMode::Dipole::SetTiltDeg(prm.field.dipoleTilt_deg);
}


//======================================================================================
// SECTION 13.0a — SINGLE-POINT RIGIDITY CLASSIFICATION DEBUG SCAN
//======================================================================================
//
// This diagnostic is intentionally independent of the main POINTS/SHELLS output grid.
// It evaluates the same TraceAllowed3D() kernel used by CutoffForDir_GV() at one
// spherical lon/lat/alt location and writes the allowed/forbidden classification as a
// function of rigidity.  For a clean centered-dipole vertical test the classification
// should be monotonic around the analytic Störmer cutoff: low R forbidden, high R
// allowed.  If Rmin is already classified as allowed at an intermediate latitude, the
// production cutoff search correctly returns Rmin, but the scan file reveals whether
// that is a real penumbra/non-monotonic effect or a trajectory-classification error.
//
// Input controls live in #CUTOFF_RIGIDITY:
//   CUTOFF_DEBUG_RIGIDITY_SCAN T
//   CUTOFF_DEBUG_SCAN_LON      <deg>
//   CUTOFF_DEBUG_SCAN_LAT      <deg>
//   CUTOFF_DEBUG_SCAN_ALT      <km>
//   CUTOFF_DEBUG_SCAN_N        <N>
//   CUTOFF_DEBUG_SCAN_FILE     <file>
//======================================================================================

static V3 DebugScanSphericalPosition_m_(double lon_deg, double lat_deg, double alt_km) {
    const double deg2rad = M_PI / 180.0;
    const double lon = lon_deg * deg2rad;
    const double lat = lat_deg * deg2rad;
    const double r_m = _EARTH__RADIUS_ + alt_km * 1000.0;

    const double clat = std::cos(lat);
    return V3{r_m * clat * std::cos(lon),
              r_m * clat * std::sin(lon),
              r_m * std::sin(lat)};
}

static void AddDebugRigidityValue_(std::vector<double>& values,
                                   double R_GV,
                                   double Rmin_GV,
                                   double Rmax_GV) {
    if (!std::isfinite(R_GV) || !(R_GV > 0.0)) return;

    // Keep the scan focused on the actual production bracket.  Tiny tolerance prevents
    // roundoff from dropping exact endpoints.
    const double tol = 1.0e-12 * std::max(std::fabs(Rmin_GV), std::fabs(Rmax_GV));
    if (R_GV < Rmin_GV - tol || R_GV > Rmax_GV + tol) return;

    values.push_back(std::max(Rmin_GV,std::min(Rmax_GV,R_GV)));
}

static std::vector<double> BuildDebugRigidityList_(const EarthUtil::AmpsParam& prm,
                                                   double Rmin_GV,
                                                   double Rmax_GV,
                                                   double RcStormer_GV) {
    std::vector<double> values;
    values.reserve((size_t)std::max(8,prm.cutoff.debugScanN) + 32);

    AddDebugRigidityValue_(values,Rmin_GV,Rmin_GV,Rmax_GV);
    AddDebugRigidityValue_(values,Rmax_GV,Rmin_GV,Rmax_GV);

    // A small set of human-readable rigidity landmarks makes the output easy to inspect
    // for the current dipole-shell problem without requiring the user to tune the input.
    const double landmarks[] = {0.03,0.05,0.1,0.2,0.5,0.8,1.0,1.16,1.5,2.0,5.0,10.0,20.0};
    for (double R : landmarks) AddDebugRigidityValue_(values,R,Rmin_GV,Rmax_GV);

    // Add values clustered around the analytic vertical cutoff if it is available.
    if (RcStormer_GV > 0.0 && std::isfinite(RcStormer_GV)) {
        const double factors[] = {0.25,0.5,0.75,0.9,0.95,0.98,1.0,1.02,1.05,1.1,1.25,1.5,2.0,4.0};
        for (double f : factors) AddDebugRigidityValue_(values,f*RcStormer_GV,Rmin_GV,Rmax_GV);
    }

    // Log-spaced production-bracket scan.  This captures unexpected non-monotonic
    // allowed/forbidden islands even when they do not coincide with the landmarks.
    const int nLog = std::max(2,prm.cutoff.debugScanN);
    const double lmin = std::log(Rmin_GV);
    const double lmax = std::log(Rmax_GV);
    for (int i=0; i<nLog; ++i) {
        const double a = (nLog==1) ? 0.0 : double(i)/double(nLog-1);
        AddDebugRigidityValue_(values,std::exp((1.0-a)*lmin + a*lmax),Rmin_GV,Rmax_GV);
    }

    std::sort(values.begin(),values.end());
    values.erase(std::unique(values.begin(),values.end(),[](double a,double b) {
        const double scale = std::max(1.0,std::max(std::fabs(a),std::fabs(b)));
        return std::fabs(a-b) <= 1.0e-8 * scale;
    }), values.end());

    return values;
}

static void WriteMode3DCutoffDebugRigidityScan_(const EarthUtil::AmpsParam& prm,
                                                cMode3DMeshFieldEval& field,
                                                double q_C,
                                                double m0_kg,
                                                const DomainBox3D& box,
                                                double Rmin_GV,
                                                double Rmax_GV) {
    double alt_km = prm.cutoff.debugScanAlt_km;
    if (!(alt_km >= 0.0)) {
        if (!prm.output.shellAlt_km.empty()) alt_km = prm.output.shellAlt_km.front();
        else alt_km = 0.0;
    }

    const V3 x0_m = DebugScanSphericalPosition_m_(prm.cutoff.debugScanLon_deg,
                                                  prm.cutoff.debugScanLat_deg,
                                                  alt_km);

    const bool isDipole = (EarthUtil::ToUpper(prm.field.model) == "DIPOLE");
    double RcStormer_GV = -1.0;
    if (isDipole) {
        EnsureDipoleAnalyticState3D(prm);
        RcStormer_GV = StormerVerticalCutoff_GV(prm,x0_m);
    }

    // Use the same vertical-arrival convention as ComputeCutoffAtPoint_GV(): the
    // arriving particle direction is toward Earth, and backtracing launches -dir.
    const V3 arrivalDir = v3unit(mul(-1.0,x0_m));
    const V3 v0 = mul(-1.0,arrivalDir);

    const double RcSelected_GV = CutoffForDir_GV(prm,field,x0_m,arrivalDir,
                                                 q_C,m0_kg,box,Rmin_GV,Rmax_GV);
    const double RcEndpointBinary_GV = CutoffForDirEndpointBinary_GV(prm,field,x0_m,arrivalDir,
                                                                     q_C,m0_kg,box,Rmin_GV,Rmax_GV);
    const double RcUpperScan_GV = CutoffForDirUpperScan_GV(prm,field,x0_m,arrivalDir,
                                                           q_C,m0_kg,box,Rmin_GV,Rmax_GV);

    const std::vector<double> Rlist = BuildDebugRigidityList_(prm,Rmin_GV,Rmax_GV,RcStormer_GV);

    const std::string fname = prm.cutoff.debugScanFile.empty()
                            ? CutoffOutputFileName("cutoff_3d_debug_rigidity_scan")
                            : prm.cutoff.debugScanFile;

    FILE* f = std::fopen(fname.c_str(),"w");
    if (!f) {
        std::ostringstream msg;
        msg << "Mode3D cutoff debug scan: cannot open output file '" << fname << "'.";
        throw std::runtime_error(msg.str());
    }

    std::fprintf(f,"TITLE=\"Mode3D Cutoff Rigidity Debug Scan\"\n");
    std::fprintf(f,"VARIABLES=\"R_GV\" \"allowed\" \"expected_allowed_stormer\" \"R_over_Rc_stormer\" \"Rc_selected_GV\" \"Rc_endpoint_binary_GV\" \"Rc_upper_scan_GV\" \"Rc_stormer_GV\"\n");
    std::fprintf(f,"# field_model=%s sampling=%s\n",prm.field.model.c_str(),prm.cutoff.sampling.c_str());
    std::fprintf(f,"# lon_deg=% .12e lat_deg=% .12e alt_km=% .12e\n",
                 prm.cutoff.debugScanLon_deg,prm.cutoff.debugScanLat_deg,alt_km);
    std::fprintf(f,"# x_m=% .12e y_m=% .12e z_m=% .12e r_Re=% .12e\n",
                 x0_m.x,x0_m.y,x0_m.z,v3norm(x0_m)/_EARTH__RADIUS_);
    std::fprintf(f,"# Rmin_GV=% .12e Rmax_GV=% .12e Rc_selected_GV=% .12e Rc_endpoint_binary_GV=% .12e Rc_upper_scan_GV=% .12e Rc_stormer_GV=% .12e\n",
                 Rmin_GV,Rmax_GV,RcSelected_GV,RcEndpointBinary_GV,RcUpperScan_GV,RcStormer_GV);
    std::fprintf(f,"# cutoff_search_algorithm=%s cutoff_upper_scan_n=%d\n",
                 prm.cutoff.searchAlgorithm.c_str(),CutoffUpperScanPointCount_(prm));
    std::fprintf(f,"# expected_allowed_stormer is meaningful only for FIELD_MODEL=DIPOLE and vertical cutoff.\n");
    std::fprintf(f,"ZONE T=\"lon=%g lat=%g alt=%g km\" I=%zu F=POINT\n",
                 prm.cutoff.debugScanLon_deg,prm.cutoff.debugScanLat_deg,alt_km,Rlist.size());

    for (double R_GV : Rlist) {
        const bool allowed = TraceAllowed3D(prm,field,x0_m,v0,R_GV,q_C,m0_kg,box);
        const int expected = (RcStormer_GV > 0.0 && std::isfinite(RcStormer_GV))
                           ? ((R_GV >= RcStormer_GV) ? 1 : 0)
                           : -1;
        const double rover = (RcStormer_GV > 0.0 && std::isfinite(RcStormer_GV))
                           ? R_GV/RcStormer_GV
                           : -1.0;

        std::fprintf(f,"% .12e %d %d % .12e % .12e % .12e % .12e % .12e\n",
                     R_GV,allowed ? 1 : 0,expected,rover,
                     RcSelected_GV,RcEndpointBinary_GV,RcUpperScan_GV,RcStormer_GV);
    }

    std::fclose(f);

    std::cout << "Mode3D cutoff debug scan wrote: " << fname
              << " (lon=" << prm.cutoff.debugScanLon_deg
              << " deg, lat=" << prm.cutoff.debugScanLat_deg
              << " deg, alt=" << alt_km << " km, "
              << Rlist.size() << " rigidity samples)\n";
}


static void WriteMode3DCutoffDebugExitTrace_(const EarthUtil::AmpsParam& prm,
                                             cMode3DMeshFieldEval& field,
                                             double q_C,
                                             double m0_kg,
                                             const DomainBox3D& box,
                                             double Rmin_GV,
                                             double Rmax_GV) {
    // Single-point trajectory-exit diagnostic. Unlike the compact rigidity scan above,
    // this file is about the terminal state of the trajectory classifier. It answers
    // the question: when TraceAllowed3D says "allowed", did the trajectory actually
    // leave the Mode3D outer box, through which face, and with what numerical overshoot?
    double alt_km = prm.cutoff.debugExitAlt_km;
    if (!(alt_km >= 0.0)) {
        if (!prm.output.shellAlt_km.empty()) alt_km = prm.output.shellAlt_km.front();
        else alt_km = 0.0;
    }

    const V3 x0_m = DebugScanSphericalPosition_m_(prm.cutoff.debugExitLon_deg,
                                                  prm.cutoff.debugExitLat_deg,
                                                  alt_km);

    const bool isDipole = (EarthUtil::ToUpper(prm.field.model) == "DIPOLE");
    double RcStormer_GV = -1.0;
    if (isDipole) {
        EnsureDipoleAnalyticState3D(prm);
        RcStormer_GV = StormerVerticalCutoff_GV(prm,x0_m);
    }

    const V3 arrivalDir = v3unit(mul(-1.0,x0_m));
    const V3 v0 = mul(-1.0,arrivalDir);

    std::vector<double> Rlist;
    if (prm.cutoff.debugExitR_GV > 0.0) {
        Rlist.push_back(prm.cutoff.debugExitR_GV);
    }
    else {
        EarthUtil::AmpsParam tmp = prm;
        tmp.cutoff.debugScanN = (prm.cutoff.debugExitN > 0) ? prm.cutoff.debugExitN : prm.cutoff.debugScanN;
        Rlist = BuildDebugRigidityList_(tmp,Rmin_GV,Rmax_GV,RcStormer_GV);
    }

    const std::string fname = prm.cutoff.debugExitFile.empty()
                            ? CutoffOutputFileName("cutoff_3d_debug_exit_trace")
                            : prm.cutoff.debugExitFile;

    FILE* f = std::fopen(fname.c_str(),"w");
    if (!f) {
        std::ostringstream msg;
        msg << "Mode3D cutoff exit diagnostic: cannot open output file '" << fname << "'.";
        throw std::runtime_error(msg.str());
    }

    std::fprintf(f,"TITLE=\"Mode3D Cutoff Trajectory Exit Diagnostic\"\n");
    std::fprintf(f,"VARIABLES=\"R_GV\" \"allowed\" \"reason_id\" \"reason\" ");
    std::fprintf(f,"\"t_s\" \"n_steps\" \"path_Re\" \"R_terminal_GV\" \"rel_dR\" ");
    std::fprintf(f,"\"face\" \"box_margin_km\" \"overshoot_km\" \"alpha_cross\" ");
    std::fprintf(f,"\"x_exit_km\" \"y_exit_km\" \"z_exit_km\" ");
    std::fprintf(f,"\"r_exit_Re\" \"lon_exit_deg\" \"lat_exit_deg\" ");
    std::fprintf(f,"\"x_cross_km\" \"y_cross_km\" \"z_cross_km\" ");
    std::fprintf(f,"\"r_cross_Re\" \"lon_cross_deg\" \"lat_cross_deg\" ");
    std::fprintf(f,"\"rel_dP_axis\" \"P_axis0\" \"P_axis_terminal\" \"Rc_stormer_GV\"\n");
    std::fprintf(f,"# field_model=%s sampling=VERTICAL mover=%s\n",
                 prm.field.model.c_str(),MoverTypeName3D_(GetDefaultMoverType()));
    std::fprintf(f,"# lon_deg=% .12e lat_deg=% .12e alt_km=% .12e\n",
                 prm.cutoff.debugExitLon_deg,prm.cutoff.debugExitLat_deg,alt_km);
    std::fprintf(f,"# x0_m=% .12e y0_m=% .12e z0_m=% .12e r0_Re=% .12e\n",
                 x0_m.x,x0_m.y,x0_m.z,v3norm(x0_m)/_EARTH__RADIUS_);
    std::fprintf(f,"# Rmin_GV=% .12e Rmax_GV=% .12e Rc_stormer_GV=% .12e\n",
                 Rmin_GV,Rmax_GV,RcStormer_GV);
    std::fprintf(f,"# reason_id: 0=OUTER_BOX 1=INNER_SPHERE_PRE 2=INNER_SPHERE_STEP 3=TIME_LIMIT 4=STEP_LIMIT 5=DISTANCE_LIMIT 6=INVALID_DT 7=UNKNOWN\n");
    std::fprintf(f,"# allowed=1 is valid only when reason=OUTER_BOX. Any timeout/step/distance/inner-sphere reason is forbidden.\n");
    std::fprintf(f,"# rel_dR checks rigidity conservation in E=0. rel_dP_axis checks canonical angular momentum about the dipole axis for FIELD_MODEL=DIPOLE.\n");
    std::fprintf(f,"ZONE T=\"exit diagnostic lon=%g lat=%g alt=%g km\" I=%zu F=POINT\n",
                 prm.cutoff.debugExitLon_deg,prm.cutoff.debugExitLat_deg,alt_km,Rlist.size());

    auto lonlat = [](const V3& x, double& lon_deg, double& lat_deg) {
        const double r = v3norm(x);
        if (!(r > 0.0)) { lon_deg = 0.0; lat_deg = 0.0; return; }
        lon_deg = std::atan2(x.y,x.x)*180.0/M_PI;
        if (lon_deg < 0.0) lon_deg += 360.0;
        lat_deg = std::asin(std::max(-1.0,std::min(1.0,x.z/r)))*180.0/M_PI;
    };

    for (double R_GV : Rlist) {
        TraceDetailedResult3D tr = TraceDetailed3D_(prm,field,x0_m,v0,R_GV,q_C,m0_kg,box);
        double lonExit=0.0, latExit=0.0, lonCross=0.0, latCross=0.0;
        lonlat(tr.xTerminal_m,lonExit,latExit);
        lonlat(tr.xBoundary_m,lonCross,latCross);

        std::fprintf(f,
            "% .12e %d %d \"%s\" % .12e %d % .12e % .12e % .12e \"%s\" % .12e % .12e % .12e "
            "% .12e % .12e % .12e % .12e % .12e % .12e "
            "% .12e % .12e % .12e % .12e % .12e % .12e "
            "% .12e % .12e % .12e % .12e\n",
            R_GV,
            tr.allowed ? 1 : 0,
            (int)tr.reason,
            TraceExitReasonName3D_(tr.reason),
            tr.t_s,
            tr.nSteps,
            tr.path_m/_EARTH__RADIUS_,
            tr.R_terminal_GV,
            tr.relRigidityError,
            tr.boundaryFace.c_str(),
            tr.terminalBoxMargin_m/1000.0,
            tr.boundaryOvershoot_m/1000.0,
            tr.boundaryAlpha,
            tr.xTerminal_m.x/1000.0,tr.xTerminal_m.y/1000.0,tr.xTerminal_m.z/1000.0,
            v3norm(tr.xTerminal_m)/_EARTH__RADIUS_,lonExit,latExit,
            tr.xBoundary_m.x/1000.0,tr.xBoundary_m.y/1000.0,tr.xBoundary_m.z/1000.0,
            v3norm(tr.xBoundary_m)/_EARTH__RADIUS_,lonCross,latCross,
            tr.relPAxisError,tr.pAxis0,tr.pAxisTerminal,RcStormer_GV);
    }

    std::fclose(f);

    std::cout << "Mode3D cutoff exit diagnostic wrote: " << fname
              << " (lon=" << prm.cutoff.debugExitLon_deg
              << " deg, lat=" << prm.cutoff.debugExitLat_deg
              << " deg, alt=" << alt_km << " km, "
              << Rlist.size() << " trajectory traces)\n";
}

static void WriteTecplot3DPoints(const EarthUtil::AmpsParam& prm,
                                 const std::vector<double>& Rc,
                                 const std::vector<double>& Emin,
                                 int nLoc) {
    const std::string fname = CutoffOutputFileName("cutoff_3d_points");
    FILE* f = std::fopen(fname.c_str(), "w");
    if (!f) throw std::runtime_error("Cannot write " + fname);

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

static void WriteTecplot3DPoints_DipoleAnalyticCompare(const EarthUtil::AmpsParam& prm,
                                                       const std::vector<double>& Rc_num_GV,
                                                       int nLoc) {
    const std::string fname = CutoffOutputFileName("cutoff_3d_points_dipole_compare");
    FILE* f = std::fopen(fname.c_str(), "w");
    if (!f) throw std::runtime_error("Cannot write Tecplot file: " + fname);

    EnsureDipoleAnalyticState3D(prm);

    std::fprintf(f, "TITLE=\"Mode3D Dipole Cutoff Rigidity: Numeric vs Analytic Vertical\"\n");
    std::fprintf(f, "VARIABLES=\"id\",\"x\",\"y\",\"z\",\"Rc_num_GV\",\"Rc_vert_GV\",\"rel_err\"\n");
    std::fprintf(f, "ZONE T=\"points\" I=%d F=POINT\n", nLoc);

    for (int i = 0; i < nLoc; ++i) {
        const auto& P = prm.output.points[(size_t)i];
        const V3 x_m{P.x*1000.0, P.y*1000.0, P.z*1000.0};

        const double Rc_vert = StormerVerticalCutoff_GV(prm, x_m);
        const double Rc_num  = Rc_num_GV[(size_t)i];
        const double rel     = (Rc_vert > 0.0 && Rc_num > 0.0)
                             ? (Rc_num - Rc_vert) / Rc_vert
                             : 0.0;

        std::fprintf(f, "%d %e %e %e %e %e %e\n",
                     i, P.x, P.y, P.z, Rc_num, Rc_vert, rel);
    }

    std::fclose(f);
}

static void WriteTecplot3DShells(const EarthUtil::AmpsParam& prm,
                                 const std::vector<std::vector<double>>& RcShell,
                                 const std::vector<std::vector<double>>& EminShell,
                                 int nLon, int nLat, double d_deg) {
    const std::string fname = CutoffOutputFileName("cutoff_3d_shells");
    FILE* f = std::fopen(fname.c_str(), "w");
    if (!f) throw std::runtime_error("Cannot write " + fname);

    std::fprintf(f, "TITLE=\"Mode3D Cutoff Rigidity (SHELLS)\"\n");
    std::fprintf(f, "VARIABLES=\"lon_deg\",\"lat_deg\",\"Rc_GV\",\"Emin_MeV\"\n");

    for (size_t s = 0; s < prm.output.shellAlt_km.size(); ++s) {
        const double alt_km = prm.output.shellAlt_km[s];
        std::fprintf(f, "ZONE T=\"alt_km=%g\" I=%d J=%d F=POINT\n",
                     alt_km, nLon, nLat);

        const int nPts = nLon * nLat;
        for (int k = 0; k < nPts; ++k) {
            const int    iLon = k % nLon;
            const int    jLat = k / nLon;
            double lat = -90.0 + d_deg * jLat;
            if (lat > 90.0) lat = 90.0;
            const double lon  = d_deg * iLon;

            std::fprintf(f, "%e %e %e %e\n",
                         lon, lat,
                         RcShell[s][(size_t)k],
                         EminShell[s][(size_t)k]);
        }
    }
    std::fclose(f);
}

static void WriteTecplot3DDirectionalMap_Location(
                                 const EarthUtil::AmpsParam& prm,
                                 int locId,
                                 const V3& x0_m,
                                 const DirectionalMapConfig3D& cfg,
                                 const std::vector<double>& RcCell,
                                 double qabs,
                                 double m0_kg) {
    (void)prm; // Reserved for future metadata fields; keep writer signature symmetric with other writers.

    // One Tecplot file is written per observation location.  This mirrors the
    // existing gridless naming/output strategy but uses a Mode3D-specific stem so
    // the mesh-based standalone products do not overwrite gridless products.
    const std::string fname = DirectionalMapOutputFileName3D(locId);
    FILE* f = std::fopen(fname.c_str(), "w");
    if (!f) throw std::runtime_error("Cannot write Tecplot directional map: " + fname);

    const double x_km = x0_m.x / 1000.0;
    const double y_km = x0_m.y / 1000.0;
    const double z_km = x0_m.z / 1000.0;

    std::fprintf(f, "TITLE=\"Mode3D directional cutoff rigidity sky-map (location %d)\"\n", locId);
    std::fprintf(f, "VARIABLES=\"lon_deg\",\"lat_deg\",\"Rc_GV\",\"Emin_MeV\"\n");

    // The zone title records the observation point in km.  The direction labels
    // are SM lon/lat when SPICE SM->GSM was available; otherwise the code used
    // identity and the labels are effectively GSM.  The fallback status is printed
    // in the run banner rather than repeated in every data row.
    std::fprintf(f,
        "ZONE T=\"loc=%d x_km=%g y_km=%g z_km=%g frame=%s\" I=%d J=%d F=POINT\n",
        locId, x_km, y_km, z_km,
        cfg.spiceOk ? "SM" : "GSM_fallback",
        cfg.nLon, cfg.nLat);

    // Row ordering matches DirectionalMapCellDirectionGSM3D():
    //   cellId = iLon + nLon*jLat
    // Tecplot reads this as an I=nLon, J=nLat POINT zone.
    for (int j=0; j<cfg.nLat; ++j) {
        double lat_deg = -90.0 + cfg.latRes_deg * j;
        if (lat_deg > 90.0) lat_deg = 90.0;

        for (int i=0; i<cfg.nLon; ++i) {
            const int cellId = i + cfg.nLon*j;
            const double lon_deg = cfg.lonRes_deg * i;
            const double rc = RcCell[(size_t)cellId];

            double Emin = -1.0;
            if (rc > 0.0) {
                const double pCut = MomentumFromRigidity_GV(rc, qabs);
                Emin = KineticEnergyFromMomentum_MeV(pCut, m0_kg);
            }

            std::fprintf(f, "%e %e %e %e\n", lon_deg, lat_deg, rc, Emin);
        }
    }

    std::fclose(f);
}

static void WriteTecplot3DShells_DipoleAnalyticCompare(
                                 const EarthUtil::AmpsParam& prm,
                                 const std::vector<std::vector<double>>& RcShell,
                                 int nLon, int nLat, double d_deg) {
    const std::string fname = CutoffOutputFileName("cutoff_3d_shells_dipole_compare");
    FILE* f = std::fopen(fname.c_str(), "w");
    if (!f) throw std::runtime_error("Cannot write Tecplot file: " + fname);

    EnsureDipoleAnalyticState3D(prm);

    std::fprintf(f, "TITLE=\"Mode3D Dipole Cutoff Rigidity (Shells): Numeric vs Analytic Vertical\"\n");
    std::fprintf(f, "VARIABLES=\"lon_deg\",\"lat_deg\",\"x_km\",\"y_km\",\"z_km\","
                    "\"Rc_num_GV\",\"Rc_vert_GV\",\"rel_err\"\n");

    for (size_t s = 0; s < prm.output.shellAlt_km.size(); ++s) {
        const double alt_km = prm.output.shellAlt_km[s];
        std::fprintf(f, "ZONE T=\"alt_km=%g\" I=%d J=%d F=POINT\n",
                     alt_km, nLon, nLat);

        const double r_m = _EARTH__RADIUS_ + alt_km * 1000.0;
        const int nPts = nLon * nLat;

        for (int k = 0; k < nPts; ++k) {
            const int    iLon = k % nLon;
            const int    jLat = k / nLon;
            double lat_deg = -90.0 + d_deg * jLat;
            if (lat_deg > 90.0) lat_deg = 90.0;
            const double lon_deg = d_deg * iLon;

            const double lon = lon_deg * M_PI / 180.0;
            const double lat = lat_deg * M_PI / 180.0;
            const double clat = std::cos(lat);

            const V3 x_m{r_m * clat * std::cos(lon),
                         r_m * clat * std::sin(lon),
                         r_m * std::sin(lat)};

            const double x_km = x_m.x / 1000.0;
            const double y_km = x_m.y / 1000.0;
            const double z_km = x_m.z / 1000.0;

            const double Rc_vert = StormerVerticalCutoff_GV(prm, x_m);
            const double Rc_num  = RcShell[s][(size_t)k];
            const double rel     = (Rc_vert > 0.0 && Rc_num > 0.0)
                                 ? (Rc_num - Rc_vert) / Rc_vert
                                 : 0.0;

            std::fprintf(f, "%e %e %e %e %e %e %e %e\n",
                         lon_deg, lat_deg, x_km, y_km, z_km,
                         Rc_num, Rc_vert, rel);
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

void SetCutoffOutputFileSuffix(const std::string& suffix) {
    gCutoffOutputFileSuffix = suffix;
}

int RunCutoffRigidity(const EarthUtil::AmpsParam& prm, bool requestedProgressBar) {

    //==================================================================================
    // 14.0 — Progress-reporting default for user-facing Mode3D cutoff runs
    //==================================================================================
    // Historically the second argument defaulted to false, and the standalone caller was
    // expected to pass true when the command-line -mode 3d workflow should display a
    // progress bar.  In practice this is fragile: an older caller, a stale object file, or
    // another direct call to RunCutoffRigidity(prm) silently reverts the banner to
    //
    //     Progress bar   : OFF
    //
    // even though the standalone 3-D cutoff calculation is the interactive path where
    // progress feedback is needed most.
    //
    // Make the default robust inside the implementation itself.  Explicit true remains
    // true; omitted/false arguments are promoted to true for this Mode3D driver.  The
    // progress-enabled branch below still computes exactly the same cutoffs; it only adds
    // synchronized batches and rank-0 progress output.  This guarantees that any
    // user-facing standalone 3-D cutoff entry path reports progress consistently.
    //==================================================================================
    (void)requestedProgressBar;

    // Force progress reporting in the standalone 3-D cutoff implementation.
    //
    // This local constant is intentionally not derived from the function argument.
    // The old API default was false, and stale/direct call sites can still pass false.
    // By using a local compile-time constant for the banner, progress counters, and
    // compute-branch selection below, the user-facing standalone cutoff path cannot
    // silently fall back to the no-progress branch.
    const bool showProgressBar = true;

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
    MPI_Comm_rank(MPI_GLOBAL_COMMUNICATOR, &mpiRank);
    MPI_Comm_size(MPI_GLOBAL_COMMUNICATOR, &mpiSize);

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
    // Match the gridless cutoff solver exactly: R_INNER in #DOMAIN_BOUNDARY is
    // the loss-sphere radius used for trajectory classification. Do not hard-code
    // _EARTH__RADIUS_ here; several validation inputs intentionally use R_INNER
    // slightly above the surface (for example 1.01 Re) to avoid grazing-boundary
    // numerical ambiguity.
    box.rInner = prm.domain.rInner * 1000.0;
    if (!(box.rInner > 0.0)) {
        throw std::runtime_error(
            "Mode3D cutoff: invalid R_INNER in #DOMAIN_BOUNDARY; value must be positive.");
    }

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
    // 14.6b — Optional directional cutoff sky-map setup
    //==================================================================================
    //
    // DIRECTIONAL_MAP is independent of CUTOFF_SAMPLING:
    //   - CUTOFF_SAMPLING controls the primary scalar cutoff written to
    //     cutoff_3d_points.dat or cutoff_3d_shells.dat.
    //   - DIRECTIONAL_MAP=T requests an additional diagnostic product: Rc for
    //     every lon/lat arrival-direction cell at every observation location.
    //
    // This is intentionally enabled for POINTS, TRAJECTORY, and SHELLS.  For
    // SHELLS this can create many files, but the keyword is explicit and the old
    // particle-based cutoff code also wrote one directional map per location.
    //==================================================================================

    const DirectionalMapConfig3D dirMapCfg = ConfigureDirectionalMap3D(prm);
    const long long tasksPerLocation =
        1LL + (dirMapCfg.enabled ? static_cast<long long>(dirMapCfg.nCells) : 0LL);

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
            << "Cutoff search  : " << prm.cutoff.searchAlgorithm
            << " (upper-scan N=" << CutoffUpperScanPointCount_(prm) << ")\n"
            << "Sampling       : " << (samplingVertical ? "VERTICAL" : "ISOTROPIC") << "\n";

        if (!samplingVertical)
            std::cout << "N_directions   : " << dirs.size()
                      << " (Fibonacci sphere)\n";

        std::cout
            << "N_locations    : " << nLoc  << "\n"
            << "Directional map: " << (dirMapCfg.enabled ? "ON" : "OFF") << "\n";

        if (dirMapCfg.enabled) {
            std::cout
                << "  DIRMAP grid  : " << dirMapCfg.nLon << " x " << dirMapCfg.nLat
                << " (" << dirMapCfg.nCells << " directions/location)\n"
                << "  DIRMAP res   : lon " << dirMapCfg.lonRes_deg
                << " deg, lat " << dirMapCfg.latRes_deg << " deg\n"
                << "  DIRMAP frame : "
                << (dirMapCfg.spiceOk ? "SM labels -> GSM tracing"
                                      : "GSM fallback (SM->GSM SPICE unavailable)")
                << "\n";
        }

        std::cout
            << "MPI ranks      : " << mpiSize << " (global B cache; scheduler selected below)\n"
            << "Progress bar   : ON\n";

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
    // 14.7a — Optional single-point rigidity classification scan
    //==================================================================================
    //
    // This debug test uses the same TraceAllowed3D/CutoffForDir_GV kernels as the full
    // cutoff calculation, but only for one lon/lat/alt point.  It is intentionally run
    // by rank 0 before the MPI location decomposition so it can diagnose numerical
    // classification/bracketing problems independently of MPI gather/order issues.
    //==================================================================================

    if (prm.cutoff.debugRigidityScan) {
        if (mpiRank == 0) {
            cMode3DMeshFieldEval debugField(prm);
            WriteMode3DCutoffDebugRigidityScan_(prm,debugField,q_C,m0,box,Rmin,Rmax);
            std::cout.flush();
        }
        MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);
    }

    if (prm.cutoff.debugExitTrace) {
        if (mpiRank == 0) {
            cMode3DMeshFieldEval debugField(prm);
            WriteMode3DCutoffDebugExitTrace_(prm,debugField,q_C,m0,box,Rmin,Rmax);
            std::cout.flush();
        }
        MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);
    }

    //==================================================================================
    // 14.8 — Two-level work scheduling over global Mode3D locations
    //==================================================================================
    // The Mode3D mesh-field snapshot is replicated/readable on every MPI rank before
    // this cutoff solver starts.  Therefore any rank can compute any output location.
    // This section uses that property to support several inter-rank schedulers:
    //
    //   STATIC
    //     Rank r receives one contiguous interval.  This is useful only for debugging
    //     because shell maps often group expensive latitudes together and imbalance the
    //     job badly.
    //
    //   BLOCK_CYCLIC
    //     Rank r receives r, r+nRanks, r+2*nRanks, ... .  This was the previous
    //     improvement over STATIC and remains deterministic/reproducible, but it still
    //     cannot react when one rank receives a set of unusually long trajectories.
    //
    //   DYNAMIC
    //     Ranks fetch chunks from an MPI one-sided atomic counter.  Once a rank finishes
    //     its current chunk, it immediately asks for another.  This is the two-level
    //     scheduler requested for the cutoff calculation: the MPI level dynamically
    //     balances chunks between ranks, and the THREADS/OpenMP/SERIAL backend below
    //     schedules the locations inside each chunk within the rank.
    //
    // Important MPI-threading rule:
    //   Only the rank/main thread calls FetchNextChunkStart().  Worker threads never call
    //   MPI.  The implementation therefore works with the normal MPI_THREAD_FUNNELED or
    //   even single-thread MPI initialization model and does not require MPI_THREAD_MULTIPLE.
    //==================================================================================

    //==================================================================================
    // 14.8a — Intra-rank and inter-rank backend selection
    //==================================================================================

    const CutoffParallelBackend_ cutoffBackend = ResolveCutoffParallelBackend_(prm);
    const int cutoffThreadCount = ResolveCutoffThreadCount_(prm,cutoffBackend);
    const Earth::Mode3D::MpiScheduler mpiScheduler =
        Earth::Mode3D::ResolveMpiScheduler(prm,"Mode3D cutoff");
    const long long mpiDynamicChunk = Earth::Mode3D::ResolveMpiDynamicChunk(
        prm,cutoffThreadCount,static_cast<long long>(nLoc));

    ApplyWideAffinityForDirectCutoffThreadsOnce_(cutoffBackend,cutoffThreadCount);

#ifdef _OPENMP
    if (cutoffBackend == CutoffParallelBackend_::OPENMP && cutoffThreadCount > 0) {
        omp_set_num_threads(cutoffThreadCount);
    }
#endif

    if (mpiRank == 0) {
        std::cout
            << "Cutoff backend : " << CutoffParallelBackendName_(cutoffBackend) << "\n"
            << "Cutoff workers : " << cutoffThreadCount
            << " per rank (from DENSITY_THREADS/AMPS_MODE3D_DENSITY_THREADS)\n"
            << "MPI scheduler  : " << Earth::Mode3D::MpiSchedulerName(mpiScheduler) << "\n";
        if (mpiScheduler == Earth::Mode3D::MpiScheduler::DYNAMIC) {
            std::cout
                << "MPI dyn chunk  : " << mpiDynamicChunk
                << " global location(s) per atomic fetch\n";
        }
        if (cutoffBackend == CutoffParallelBackend_::THREADS) {
            std::cout
                << "Affinity mode  : wide MPI-rank mask via "
                << "PIC::Parallel::SetWideAffinityForScheduler()\n";
        }
        std::cout.flush();
    }

    //==================================================================================
    // 14.8b — Rank-local result arrays indexed by GLOBAL location id
    //==================================================================================
    // Dynamic MPI scheduling makes the number and order of locations processed by each
    // rank unknown until runtime.  To avoid any gather-order assumptions, each rank keeps
    // full-length result arrays initialized to sentinel values and writes directly to
    // slot globalIdx.  Because the MPI scheduler guarantees that each globalIdx is issued
    // to exactly one rank, MPI_MAX at the end selects the one non-sentinel contribution.
    //
    // Sentinel convention:
    //   Rc/Emin/dirMap = -1.0 means "this rank did not compute this location".
    //   Valid cutoff values are >=0, so MPI_MAX is safe even for the polar analytic
    //   cutoff where Rc can be zero.
    //==================================================================================

    std::vector<double> rcRank((size_t)nLoc,-1.0);
    std::vector<double> eminRank((size_t)nLoc,-1.0);
    std::vector<double> dirMapRank;

    if (dirMapCfg.enabled) {
        const long long nMapLL = static_cast<long long>(nLoc) *
                                 static_cast<long long>(dirMapCfg.nCells);
        if (nMapLL > static_cast<long long>(std::numeric_limits<int>::max())) {
            throw std::runtime_error(
                "Mode3D cutoff: directional-map MPI reduction count exceeds INT_MAX. "
                "Use a coarser DIRMAP grid or split the run into fewer locations.");
        }
        dirMapRank.assign((size_t)nMapLL,-1.0);
    }

    // Deterministic fallback work list for STATIC and BLOCK_CYCLIC.  DYNAMIC does not
    // use this vector; it obtains global locations from the RMA counter in chunks.
    std::vector<int> rankWorkList;
    if (mpiScheduler == Earth::Mode3D::MpiScheduler::BLOCK_CYCLIC) {
        rankWorkList.reserve((size_t)((nLoc + mpiSize - 1) / mpiSize));
        for (int globalIdx=mpiRank; globalIdx<nLoc; globalIdx+=mpiSize) {
            rankWorkList.push_back(globalIdx);
        }
    }
    else if (mpiScheduler == Earth::Mode3D::MpiScheduler::STATIC) {
        const int begin = (int)((static_cast<long long>(nLoc) * mpiRank) / mpiSize);
        const int end   = (int)((static_cast<long long>(nLoc) * (mpiRank+1)) / mpiSize);
        rankWorkList.reserve((size_t)std::max(0,end-begin));
        for (int globalIdx=begin; globalIdx<end; ++globalIdx) rankWorkList.push_back(globalIdx);
    }

    const int nLocalStatic = static_cast<int>(rankWorkList.size());

    //==================================================================================
    // 14.8c — Global progress bookkeeping
    //==================================================================================
    // Progress collectives must be executed in the same order by every rank.  For the
    // deterministic STATIC/BLOCK_CYCLIC schedulers we can still use synchronized batches.
    // For DYNAMIC, ranks finish chunks at different times and repeatedly call the MPI RMA
    // counter; inserting progress Allreduces inside that loop would either require a much
    // more complicated protocol or reintroduce a global synchronization bottleneck.
    // Therefore DYNAMIC prints start and completion progress only, prioritizing load
    // balance over intermediate progress updates.
    //==================================================================================

    const long long totalLocationsGlobal = static_cast<long long>(nLoc);
    const long long totalTasksGlobal = totalLocationsGlobal * tasksPerLocation;

    long long doneLocationsLocal  = 0;
    long long doneLocationsGlobal = 0;
    long long doneTasksLocal      = 0;
    long long doneTasksGlobal     = 0;

    std::vector<int> locDonePerShellLocal((size_t)std::max(nShells,0),0);
    std::vector<int> locDonePerShellGlobal((size_t)std::max(nShells,0),0);
    std::vector<int> locTotalPerShellGlobal((size_t)std::max(nShells,0),0);

    if (isShells) {
        for (int s=0; s<nShells; ++s) locTotalPerShellGlobal[(size_t)s] = nPtsShell;
    }

    auto mode3d_now_seconds = []() -> double { return MPI_Wtime(); };
    const double progressStartTime = mode3d_now_seconds();
    double progressLastPrintTime = -1.0;

    auto mode3d_fmt_hms = [](double s) -> std::string {
        if (s < 0.0) return std::string("--:--:--");
        long long is = (long long)std::llround(s);
        long long hh = is/3600; is-=hh*3600;
        long long mm = is/60;   is-=mm*60;
        long long ss = is;
        char buf[64];
        std::snprintf(buf,sizeof(buf),"%02lld:%02lld:%02lld",hh,mm,ss);
        return std::string(buf);
    };

    auto maybePrintProgress = [&](long long doneLocations,
                                  long long doneTasks,
                                  const std::vector<int>& shellDoneGlobal,
                                  bool forcePrint) {
        if (mpiRank != 0) return;

        const double t = mode3d_now_seconds();
        if (!forcePrint) {
            if (progressLastPrintTime < 0.0) progressLastPrintTime = t;
            if (t - progressLastPrintTime < 1.0) return;
        }
        progressLastPrintTime = t;

        const double frac = (totalTasksGlobal > 0)
            ? (double(doneTasks)/double(totalTasksGlobal)) : 1.0;
        const double dt = t - progressStartTime;
        const double rate = (dt > 0.0) ? (double(doneTasks)/dt) : 0.0;
        double eta_s = -1.0;
        if (rate > 0.0 && totalTasksGlobal > doneTasks)
            eta_s = double(totalTasksGlobal-doneTasks)/rate;

        const int barW = 36;
        int filled = (int)std::floor(frac*barW + 0.5);
        if (filled < 0) filled = 0;
        if (filled > barW) filled = barW;

        std::ostringstream line;
        if (isPoints) {
            if (outputMode == "TRAJECTORY") line << "[Mode3D cutoff TRAJECTORY] ";
            else line << "[Mode3D cutoff POINTS] ";
        }
        else {
            line << "[Mode3D cutoff SHELLS " << nShells << " zones";
            if (nShells == 1) {
                line << " alt=" << prm.output.shellAlt_km[0] << "km";
            }
            else if (nShells > 1 && nShells <= 4) {
                line << " alt=";
                for (int s=0; s<nShells; ++s) {
                    if (s) line << ",";
                    line << prm.output.shellAlt_km[(size_t)s];
                }
                line << "km";
            }
            else if (nShells > 4) {
                line << " alt=" << prm.output.shellAlt_km.front()
                     << ".." << prm.output.shellAlt_km.back() << "km";
            }
            line << "] ";
        }

        line << "[rank 0/global over " << mpiSize << " MPI ranks] ";
        line << "[";
        for (int i=0; i<barW; ++i) line << (i<filled ? "#" : "-");
        line << "] ";

        line.setf(std::ios::fixed);
        line.precision(1);
        line << (frac*100.0) << "%  ";
        line << "(Loc " << doneLocations << "/" << totalLocationsGlobal
             << ", Task " << doneTasks << "/" << totalTasksGlobal;

        if (isShells) {
            line << "; Shells ";
            if (nShells <= 4) {
                for (int s=0; s<nShells; ++s) {
                    if (s) line << ", ";
                    const int shellDone  = shellDoneGlobal[(size_t)s];
                    const int shellTotal = locTotalPerShellGlobal[(size_t)s];
                    const double pct = (shellTotal > 0)
                        ? 100.0*double(shellDone)/double(shellTotal) : 100.0;
                    line << (s+1) << ":" << shellDone << "/" << shellTotal
                         << " " << pct << "%";
                }
            }
            else {
                int nCompleteShells = 0;
                int slowestShell = -1;
                double slowestFrac = 2.0;
                for (int s=0; s<nShells; ++s) {
                    const int shellDone  = shellDoneGlobal[(size_t)s];
                    const int shellTotal = locTotalPerShellGlobal[(size_t)s];
                    const double shellFrac = (shellTotal > 0)
                        ? double(shellDone)/double(shellTotal) : 1.0;
                    if (shellFrac >= 1.0) nCompleteShells++;
                    if (shellFrac < slowestFrac) {
                        slowestFrac = shellFrac;
                        slowestShell = s;
                    }
                }
                line << nCompleteShells << "/" << nShells << " complete";
                if (slowestShell >= 0) {
                    const int shellDone  = shellDoneGlobal[(size_t)slowestShell];
                    const int shellTotal = locTotalPerShellGlobal[(size_t)slowestShell];
                    line << ", slowest " << (slowestShell+1) << ":"
                         << shellDone << "/" << shellTotal << " "
                         << (100.0*slowestFrac) << "%";
                }
            }
        }

        line << ")  ETA " << mode3d_fmt_hms(eta_s) << "\n";
        std::cout << line.str();
        std::cout.flush();
    };

    auto accountCompletedGlobalRange = [&](int begin, int end) {
        const long long nDone = static_cast<long long>(std::max(0,end-begin));
        doneLocationsLocal += nDone;
        doneTasksLocal     += nDone * tasksPerLocation;
        if (isShells) {
            for (int globalIdx=begin; globalIdx<end; ++globalIdx) {
                const int shellIdx = globalIdx / nPtsShell;
                if (shellIdx>=0 && shellIdx<nShells)
                    locDonePerShellLocal[(size_t)shellIdx]++;
            }
        }
    };

    auto accountCompletedWorkListRange = [&](int begin, int end) {
        const long long nDone = static_cast<long long>(std::max(0,end-begin));
        doneLocationsLocal += nDone;
        doneTasksLocal     += nDone * tasksPerLocation;
        if (isShells) {
            for (int localIdx=begin; localIdx<end; ++localIdx) {
                const int globalIdx = rankWorkList[(size_t)localIdx];
                const int shellIdx = globalIdx / nPtsShell;
                if (shellIdx>=0 && shellIdx<nShells)
                    locDonePerShellLocal[(size_t)shellIdx]++;
            }
        }
    };

    //==================================================================================
    // 14.9 — Location computation kernels
    //==================================================================================

    auto computeGlobalLocation = [&](int globalIdx, cMode3DMeshFieldEval& threadField) {
        const V3 x0_m = LocationToX0m(prm, globalIdx, nLon, nLat, d_deg, nPtsShell);

        const double rc = ComputeCutoffAtPoint_GV(
            prm, threadField, x0_m, dirs, samplingVertical,
            q_C, m0, box, Rmin, Rmax);

        rcRank[(size_t)globalIdx] = rc;
        if (rc > 0.0) {
            const double pCut = MomentumFromRigidity_GV(rc, qabs);
            eminRank[(size_t)globalIdx] = KineticEnergyFromMomentum_MeV(pCut, m0);
        }

        if (dirMapCfg.enabled) {
            const size_t base = (size_t)globalIdx * (size_t)dirMapCfg.nCells;
            for (int cellId=0; cellId<dirMapCfg.nCells; ++cellId) {
                const V3 dir_gsm = DirectionalMapCellDirectionGSM3D(dirMapCfg, cellId);
                dirMapRank[base + (size_t)cellId] = CutoffForDir_GV(
                    prm, threadField, x0_m, dir_gsm,
                    q_C, m0, box, Rmin, Rmax);
            }
        }
    };

    auto computeGlobalRange = [&](int begin, int end) {
        if (end <= begin) return;

        if (cutoffBackend == CutoffParallelBackend_::THREADS && cutoffThreadCount > 1) {
            const int nWork = end - begin;
            const int nWorkers = std::max(1,std::min(cutoffThreadCount,nWork));
            std::atomic<int> nextGlobalIdx(begin);
            std::vector<std::thread> workers;
            workers.reserve((size_t)nWorkers);

            for (int iw=0; iw<nWorkers; ++iw) {
                workers.emplace_back([&]() {
                    cMode3DMeshFieldEval threadField(prm);
                    for (;;) {
                        const int globalIdx = nextGlobalIdx.fetch_add(1,std::memory_order_relaxed);
                        if (globalIdx >= end) break;
                        computeGlobalLocation(globalIdx,threadField);
                    }
                });
            }

            for (std::thread& worker : workers) worker.join();
            return;
        }

        if (cutoffBackend == CutoffParallelBackend_::SERIAL) {
            cMode3DMeshFieldEval threadField(prm);
            for (int globalIdx=begin; globalIdx<end; ++globalIdx) {
                computeGlobalLocation(globalIdx,threadField);
            }
            return;
        }

#ifdef _OPENMP
#pragma omp parallel
#endif
        {
            cMode3DMeshFieldEval threadField(prm);
#ifdef _OPENMP
#pragma omp for schedule(dynamic, 1)
#endif
            for (int globalIdx=begin; globalIdx<end; ++globalIdx) {
                computeGlobalLocation(globalIdx,threadField);
            }
        }
    };

    auto computeWorkListRange = [&](int begin, int end) {
        if (end <= begin) return;

        if (cutoffBackend == CutoffParallelBackend_::THREADS && cutoffThreadCount > 1) {
            const int nWork = end - begin;
            const int nWorkers = std::max(1,std::min(cutoffThreadCount,nWork));
            std::atomic<int> nextLocalIdx(begin);
            std::vector<std::thread> workers;
            workers.reserve((size_t)nWorkers);

            for (int iw=0; iw<nWorkers; ++iw) {
                workers.emplace_back([&]() {
                    cMode3DMeshFieldEval threadField(prm);
                    for (;;) {
                        const int localIdx = nextLocalIdx.fetch_add(1,std::memory_order_relaxed);
                        if (localIdx >= end) break;
                        const int globalIdx = rankWorkList[(size_t)localIdx];
                        computeGlobalLocation(globalIdx,threadField);
                    }
                });
            }

            for (std::thread& worker : workers) worker.join();
            return;
        }

        if (cutoffBackend == CutoffParallelBackend_::SERIAL) {
            cMode3DMeshFieldEval threadField(prm);
            for (int localIdx=begin; localIdx<end; ++localIdx) {
                computeGlobalLocation(rankWorkList[(size_t)localIdx],threadField);
            }
            return;
        }

#ifdef _OPENMP
#pragma omp parallel
#endif
        {
            cMode3DMeshFieldEval threadField(prm);
#ifdef _OPENMP
#pragma omp for schedule(dynamic, 1)
#endif
            for (int localIdx=begin; localIdx<end; ++localIdx) {
                computeGlobalLocation(rankWorkList[(size_t)localIdx],threadField);
            }
        }
    };

    maybePrintProgress(0,0,locDonePerShellGlobal,true);

    if (mpiScheduler == Earth::Mode3D::MpiScheduler::DYNAMIC) {
        // Collective construction of the RMA counter.  After construction, each rank
        // independently fetches chunks.  No MPI call is made from inside worker threads.
        Earth::Mode3D::DynamicMpiLocationScheduler scheduler(
            MPI_GLOBAL_COMMUNICATOR,
            static_cast<long long>(nLoc),
            mpiDynamicChunk,
            "Mode3D cutoff");

        for (;;) {
            const long long chunkStartLL = scheduler.FetchNextChunkStart();
            if (chunkStartLL >= static_cast<long long>(nLoc)) break;

            const long long chunkEndLL = std::min(
                chunkStartLL + scheduler.ChunkSize(),
                static_cast<long long>(nLoc));

            const int chunkStart = static_cast<int>(chunkStartLL);
            const int chunkEnd   = static_cast<int>(chunkEndLL);

            computeGlobalRange(chunkStart,chunkEnd);
            accountCompletedGlobalRange(chunkStart,chunkEnd);
        }

        MPI_Allreduce(&doneLocationsLocal,&doneLocationsGlobal,
                      1,MPI_LONG_LONG,MPI_SUM,MPI_GLOBAL_COMMUNICATOR);
        MPI_Allreduce(&doneTasksLocal,&doneTasksGlobal,
                      1,MPI_LONG_LONG,MPI_SUM,MPI_GLOBAL_COMMUNICATOR);
        if (isShells) {
            MPI_Allreduce(locDonePerShellLocal.data(),locDonePerShellGlobal.data(),
                          nShells,MPI_INT,MPI_SUM,MPI_GLOBAL_COMMUNICATOR);
        }
        maybePrintProgress(doneLocationsGlobal,doneTasksGlobal,locDonePerShellGlobal,true);
    }
    else {
        int nProgressBatches = std::max(1,std::min(nLoc,200));
        if (cutoffBackend == CutoffParallelBackend_::THREADS && cutoffThreadCount > 1) {
            nProgressBatches = 1;
        }

        for (int ibatch=0; ibatch<nProgressBatches; ++ibatch) {
            const int batchBegin = (int)(
                (static_cast<long long>(nLocalStatic) * ibatch) / nProgressBatches);
            const int batchEnd = (int)(
                (static_cast<long long>(nLocalStatic) * (ibatch+1)) / nProgressBatches);

            if (batchEnd > batchBegin) {
                computeWorkListRange(batchBegin,batchEnd);
                accountCompletedWorkListRange(batchBegin,batchEnd);
            }

            MPI_Allreduce(&doneLocationsLocal,&doneLocationsGlobal,
                          1,MPI_LONG_LONG,MPI_SUM,MPI_GLOBAL_COMMUNICATOR);
            MPI_Allreduce(&doneTasksLocal,&doneTasksGlobal,
                          1,MPI_LONG_LONG,MPI_SUM,MPI_GLOBAL_COMMUNICATOR);
            if (isShells) {
                MPI_Allreduce(locDonePerShellLocal.data(),locDonePerShellGlobal.data(),
                              nShells,MPI_INT,MPI_SUM,MPI_GLOBAL_COMMUNICATOR);
            }

            maybePrintProgress(doneLocationsGlobal,doneTasksGlobal,locDonePerShellGlobal,
                               ibatch == nProgressBatches-1);
        }
    }

    //==================================================================================
    // 14.10 — MPI reduction: combine global-indexed rank-local result arrays
    //==================================================================================
    // Since each rank writes values directly into globalIdx slots and leaves all other
    // slots at -1, the result combination is a simple MPI_MAX reduction to rank 0.  This
    // is independent of the scheduler: STATIC, BLOCK_CYCLIC, and DYNAMIC all share the
    // same output path, eliminating the gather-order bugs that can appear when the local
    // work order changes.
    //==================================================================================

    std::vector<double> rcAll, eminAll, dirMapAll;
    if (mpiRank == 0) {
        rcAll.assign((size_t)nLoc,-1.0);
        eminAll.assign((size_t)nLoc,-1.0);
    }

    MPI_Reduce(rcRank.data(),
               (mpiRank==0 ? rcAll.data() : nullptr),
               nLoc,MPI_DOUBLE,MPI_MAX,0,MPI_GLOBAL_COMMUNICATOR);
    MPI_Reduce(eminRank.data(),
               (mpiRank==0 ? eminAll.data() : nullptr),
               nLoc,MPI_DOUBLE,MPI_MAX,0,MPI_GLOBAL_COMMUNICATOR);

    if (dirMapCfg.enabled) {
        const int nMap = static_cast<int>(dirMapRank.size());
        if (mpiRank == 0) dirMapAll.assign((size_t)nMap,-1.0);
        MPI_Reduce(dirMapRank.data(),
                   (mpiRank==0 ? dirMapAll.data() : nullptr),
                   nMap,MPI_DOUBLE,MPI_MAX,0,MPI_GLOBAL_COMMUNICATOR);
    }

    MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);

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
            std::cout << "Wrote: " << CutoffOutputFileName("cutoff_3d_points") << "\n";

#if _PIC_NIGHTLY_TEST_MODE_ == _PIC_MODE_ON_
            if (EarthUtil::ToUpper(prm.field.model) == "DIPOLE") {
                WriteTecplot3DPoints_DipoleAnalyticCompare(prm, rcAll, nLoc);
                std::cout << "Wrote: " << CutoffOutputFileName("cutoff_3d_points_dipole_compare") << "\n";
            }
#endif

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
            std::cout << "Wrote: " << CutoffOutputFileName("cutoff_3d_shells") << "\n";

#if _PIC_NIGHTLY_TEST_MODE_ == _PIC_MODE_ON_
            if (EarthUtil::ToUpper(prm.field.model) == "DIPOLE") {
                WriteTecplot3DShells_DipoleAnalyticCompare(prm, RcShell, nLon, nLat, d_deg);
                std::cout << "Wrote: " << CutoffOutputFileName("cutoff_3d_shells_dipole_compare") << "\n";
            }
#endif
        }

        // -----------------------------------------------------------------------------
        // Optional directional sky-map output.
        //
        // The computation and MPI_Gatherv above filled dirMapAll in global location
        // order.  We now write one Tecplot file per location, reusing the same
        // LocationToX0m() mapping used for the primary cutoff calculation so that
        // POINTS, TRAJECTORY, and SHELLS all get maps anchored at exactly the same
        // physical coordinates as the scalar Rc outputs.
        //
        // File naming:
        //   cutoff_3d_dir_map_loc_000000.dat
        //   cutoff_3d_dir_map_loc_000001.dat
        //   ...
        //
        // If a live SWMF-coupled caller installed a cutoff-output suffix with
        // SetCutoffOutputFileSuffix(), the suffix is applied before .dat by
        // DirectionalMapOutputFileName3D(), matching the scalar cutoff files.
        // -----------------------------------------------------------------------------
        if (dirMapCfg.enabled) {
            for (int locId=0; locId<nLoc; ++locId) {
                const V3 x0_m = LocationToX0m(prm, locId, nLon, nLat, d_deg, nPtsShell);

                std::vector<double> RcCell((size_t)dirMapCfg.nCells, -1.0);
                const size_t base = (size_t)locId * (size_t)dirMapCfg.nCells;
                for (int cellId=0; cellId<dirMapCfg.nCells; ++cellId) {
                    RcCell[(size_t)cellId] = dirMapAll[base + (size_t)cellId];
                }

                WriteTecplot3DDirectionalMap_Location(
                    prm, locId, x0_m, dirMapCfg, RcCell, qabs, m0);
            }

            std::cout << "Wrote: " << nLoc
                      << " directional cutoff map file(s): "
                      << "cutoff_3d_dir_map_loc_######"
                      << gCutoffOutputFileSuffix << ".dat\n";
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

//======================================================================================
// PUBLIC SHARED MESH TRAJECTORY CLASSIFIER
//======================================================================================
//
// Mode3D density/flux must compute the same transmissivity as the gridless density
// solver, but using the already-materialized AMR mesh field rather than direct
// Tsyganenko calls.  The two wrappers below provide that bridge: they expose the
// internal TraceAllowed3D() kernel through the same plain-array interface used by
// GridlessMode::TraceAllowedShared/Ex().
//
// Important implementation notes:
//   * The field evaluator is intentionally constructed per trajectory call instead of
//     cached in thread-local storage.  It stores a reference to `prm`; caching it across
//     standalone time snapshots would risk a dangling reference after the snapshot-local
//     AmpsParam copy goes out of scope.  The evaluator itself is lightweight: the heavy
//     state is the read-only AMR tree and cell data already resident in memory.
//   * The domain box is reconstructed from Mode3D::ParsedDomainMin/Max each call.  Those
//     values are set once before mesh construction in standalone mode and by the SWMF
//     pre-init hook in coupled mode, so they are the authoritative trajectory geometry.
//   * TraceAllowedMeshEx fills the same TrajectoryExitState structure as gridless so the
//     anisotropic-density code can be backend-agnostic.
//======================================================================================

bool TraceAllowedMeshEx(const EarthUtil::AmpsParam& prm,
                        const double x0_m_arr[3],
                        const double v0_unit_arr[3],
                        double R_GV,
                        Earth::GridlessMode::TrajectoryExitState* exitState,
                        double maxTraceTimeOverride_s) {
    DomainBox3D box;
    box.xMin = Earth::Mode3D::ParsedDomainMin[0];
    box.xMax = Earth::Mode3D::ParsedDomainMax[0];
    box.yMin = Earth::Mode3D::ParsedDomainMin[1];
    box.yMax = Earth::Mode3D::ParsedDomainMax[1];
    box.zMin = Earth::Mode3D::ParsedDomainMin[2];
    box.zMax = Earth::Mode3D::ParsedDomainMax[2];
    box.rInner = _EARTH__RADIUS_;

    cMode3DMeshFieldEval field(prm);
    const V3 x0_m{ x0_m_arr[0], x0_m_arr[1], x0_m_arr[2] };
    const V3 v0_unit = v3unit(V3{ v0_unit_arr[0], v0_unit_arr[1], v0_unit_arr[2] });

    const double q_C   = prm.species.charge_e * ElectronCharge;
    const double m0_kg = prm.species.mass_amu * _AMU_;

    return TraceAllowed3D(prm, field, x0_m, v0_unit, R_GV, q_C, m0_kg,
                          box, maxTraceTimeOverride_s, exitState);
}

bool TraceAllowedMesh(const EarthUtil::AmpsParam& prm,
                      const double x0_m[3],
                      const double v0_unit[3],
                      double R_GV,
                      double maxTraceTimeOverride_s) {
    return TraceAllowedMeshEx(prm, x0_m, v0_unit, R_GV, nullptr,
                              maxTraceTimeOverride_s);
}

} // namespace Mode3D
} // namespace Earth
