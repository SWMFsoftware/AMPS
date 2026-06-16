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
// Mesh/field distribution:
//   Mode3D::Run now initializes AMPS with the normal distributed MPI domain
//   decomposition.  Before this solver starts, Mode3D::GlobalMagneticField has
//   already allocated missing nonlocal leaf blocks on every rank and populated them
//   with a replicated read-only cell-centered B snapshot.  Therefore trajectory
//   tracing can remain local-memory-only even though the mesh was not created with
//   independentDomainMode=true.
//
// Point distribution: static block decomposition.
//   nLocal = nLoc / mpiSize  (floor)
//   rank k owns points [k*nLocal, (k+1)*nLocal)
//   last rank picks up the remainder: [(mpiSize-1)*nLocal, nLoc)
//
// Communication:
//   The global-B materialization is completed before entering RunCutoffRigidity().
//   During the computation phase there is no field communication; after the
//   computation phase, MPI_Gatherv collects the per-rank result arrays (Rc, Emin)
//   at rank 0.  If DIRECTIONAL_MAP=T, the same gather pattern is reused for the
//   flattened directional-map array Rc(location, lon/lat sky cell).
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
#include <limits>

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
            // the dipole helper state), so serialize it when OpenMP is enabled.
#ifdef _OPENMP
#pragma omp critical(Mode3DAnalyticMagneticFieldEval)
#endif
            {
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
          double B[3];

          PIC::CPLR::InitInterpolationStencil(xArr,node);
	  PIC::CPLR::GetBackgroundMagneticField(B);

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
    // initializes Dipole::gParams, so set the dipole axis explicitly before any
    // analytic Størmer comparison is written.
    Earth::GridlessMode::Dipole::SetTiltDeg(prm.field.dipoleTilt_deg);
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
            << "MPI ranks      : " << mpiSize << " (global B cache, static work partition)\n"
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

    //==================================================================================
    // 14.8b — Optional global progress reporting
    //==================================================================================
    // The gridless cutoff solver prints one global progress bar from its root
    // scheduler.  For Mode3D the cutoff locations are distributed statically
    // across MPI ranks, so there is no central scheduler.  To keep the same
    // user-facing behaviour without requiring MPI calls from OpenMP worker
    // threads, the progress-enabled path below computes the local work in a
    // sequence of synchronized batches:
    //
    //   1. every MPI rank computes the next batch of its own locations with
    //      OpenMP;
    //   2. all ranks enter the same MPI_Allreduce calls;
    //   3. rank 0 receives the summed number of completed locations/tasks and
    //      prints the single global progress bar.
    //
    // When showProgressBar=false, RunCutoffRigidity keeps the original fast
    // path: one OpenMP loop over this rank's entire local location range and no
    // extra inter-rank synchronization during the compute phase.
    //
    // Progress task accounting:
    //   - the primary scalar cutoff contributes one coarse task per location;
    //   - if DIRECTIONAL_MAP=T, each sky-map cell contributes one additional
    //     trajectory-search task for that same location.
    //
    // The computation is still scheduled by location (for simple MPI/OpenMP
    // ownership), so task completion advances in location-sized chunks.  The
    // Task counter nevertheless reflects the actual added directional-map work,
    // which is important because DIRMAP grids can dominate the runtime.
    //==================================================================================

    const long long totalLocationsGlobal = static_cast<long long>(nLoc);
    const long long totalTasksGlobal = totalLocationsGlobal * tasksPerLocation;

    long long doneLocationsLocal  = 0;
    long long doneLocationsGlobal = 0;
    long long doneTasksLocal      = 0;
    long long doneTasksGlobal     = 0;

    std::vector<int> locDonePerShellLocal((size_t)std::max(nShells,0),0);
    std::vector<int> locDonePerShellGlobal((size_t)std::max(nShells,0),0);
    std::vector<int> locTotalPerShellLocal((size_t)std::max(nShells,0),0);
    std::vector<int> locTotalPerShellGlobal((size_t)std::max(nShells,0),0);

    if (showProgressBar && isShells) {
        for (int globalIdx=locStart; globalIdx<locEnd; ++globalIdx) {
            const int s = globalIdx / nPtsShell;
            if (s>=0 && s<nShells) locTotalPerShellLocal[(size_t)s]++;
        }

        // The shell totals are needed by rank 0 to report per-shell progress.
        // In the MPI-parallel Mode3D algorithm different ranks may work on
        // different shells at the same time.  Therefore there is no single
        // globally meaningful "current shell".  The progress line below reports
        // total job progress and a compact per-shell completion summary instead
        // of inferring a shell label from the global linear task counter.
        // This collective is used only when the progress bar is enabled, so the
        // default no-progress path remains free of extra progress synchronization.
        MPI_Allreduce(locTotalPerShellLocal.data(),
                      locTotalPerShellGlobal.data(),
                      nShells, MPI_INT, MPI_SUM, MPI_GLOBAL_COMMUNICATOR);
    }

    auto mode3d_now_seconds = []() -> double { return MPI_Wtime(); };
    const double progressStartTime = mode3d_now_seconds();
    double progressLastPrintTime = -1.0;

    auto mode3d_fmt_hms = [](double s)->std::string {
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
        if (!showProgressBar) return;
        if (mpiRank != 0) return;

        const double t = mode3d_now_seconds();
        if (!forcePrint) {
            if (progressLastPrintTime < 0.0) progressLastPrintTime = t;
            if (t - progressLastPrintTime < 1.0) return;
        }
        progressLastPrintTime = t;

        const long long totalTasks = totalTasksGlobal;
        const long long totalLocations = totalLocationsGlobal;
        const double frac = (totalTasks>0) ?
            (double(doneTasks)/double(totalTasks)) : 1.0;

        const double dt = t - progressStartTime;
        const double rate = (dt>0.0) ? (double(doneTasks)/dt) : 0.0;
        double eta_s = -1.0;
        if (rate>0.0 && totalTasks>doneTasks)
            eta_s = double(totalTasks-doneTasks)/rate;

        const int barW = 36;
        int filled = (int)std::floor(frac*barW + 0.5);
        if (filled < 0) filled = 0;
        if (filled > barW) filled = barW;

        std::ostringstream line;

        if (isPoints) {
            if (outputMode == "TRAJECTORY") line << "[TRAJECTORY] ";
            else line << "[POINTS] ";
        }
        else {
            // Do not print "SHELL i/N" here.  With MPI domain/work splitting,
            // rank 0, rank 1, ... may simultaneously work on locations that
            // belong to different altitude shells.  A single "current shell"
            // would therefore be misleading.  Label the job as a multi-shell
            // calculation and include the actual per-shell counters below.
            line << "[SHELLS " << nShells << " zones";

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

        // Keep the same visual style as the rank-local bar, but make explicit
        // that this line is the root-rank summary over all MPI processes.
        line << "[rank 0/global over " << mpiSize << " MPI ranks] ";

        line << "[";
        for (int i=0;i<barW;i++) line << (i<filled ? "#" : "-");
        line << "] ";

        line.setf(std::ios::fixed);
        line.precision(1);
        line << (frac*100.0) << "%  ";

        line << "(Loc " << doneLocations << "/" << totalLocations << ", "
             << "Task " << doneTasks << "/" << totalTasks;

        if (isShells && nShells > 0) {
            line << "; Shells ";

            if (nShells <= 4) {
                // For the common small number of altitude shells, show every
                // shell explicitly.  These counters come from MPI_Allreduce and
                // therefore summarize work completed by all ranks, not just
                // rank 0.
                for (int s=0; s<nShells; ++s) {
                    if (s) line << ", ";

                    const int shellDone  = shellDoneGlobal[(size_t)s];
                    const int shellTotal = locTotalPerShellGlobal[(size_t)s];
                    const double shellPct = (shellTotal > 0)
                        ? 100.0*double(shellDone)/double(shellTotal)
                        : 100.0;

                    line << (s+1) << ":"
                         << shellDone << "/" << shellTotal
                         << " " << shellPct << "%";
                }
            }
            else {
                // For many shells, keep the line short: report how many shells
                // are complete and identify the least-complete shell, which is
                // the useful diagnostic when progress is imbalanced.
                int nCompleteShells = 0;
                int slowestShell = -1;
                double slowestFrac = 2.0;

                for (int s=0; s<nShells; ++s) {
                    const int shellDone  = shellDoneGlobal[(size_t)s];
                    const int shellTotal = locTotalPerShellGlobal[(size_t)s];
                    const double shellFrac = (shellTotal > 0)
                        ? double(shellDone)/double(shellTotal)
                        : 1.0;

                    if (shellFrac >= 1.0) nCompleteShells++;
                    if (shellFrac < slowestFrac) {
                        slowestFrac = shellFrac;
                        slowestShell = s;
                    }
                }

                line << nCompleteShells << "/" << nShells << " complete";
                if (slowestShell >= 0) {
                    const int shellDone = shellDoneGlobal[(size_t)slowestShell];
                    const int shellTotal = locTotalPerShellGlobal[(size_t)slowestShell];
                    line << ", slowest " << (slowestShell+1) << ":"
                         << shellDone << "/" << shellTotal
                         << " " << (100.0*slowestFrac) << "%";
                }
            }
        }

        line << ")  ETA " << mode3d_fmt_hms(eta_s) << "\n";

        std::cout << line.str();
        std::cout.flush();
    };

    if (showProgressBar) {
        maybePrintProgress(0,0,locDonePerShellGlobal,true);
    }

    // Local result arrays (indexed 0..nLocal-1)
    std::vector<double> rcLocal(nLocal, -1.0);
    std::vector<double> eminLocal(nLocal, -1.0);

    // Optional directional-map results.  The layout is rank-local but otherwise
    // identical to the final rank-0 layout:
    //
    //   dirMapLocal[ localIdx*nCells + cellId ] = Rc_GV(location, sky cell)
    //
    // localIdx is the same local location index used by rcLocal/eminLocal.
    // cellId follows the Tecplot ordering iLon + nLon*jLat.  This flat layout
    // lets MPI_Gatherv collect the maps with the same displacements as the scalar
    // cutoff arrays, simply multiplied by nCells.
    std::vector<double> dirMapLocal;
    if (dirMapCfg.enabled) {
        dirMapLocal.assign((size_t)nLocal * (size_t)dirMapCfg.nCells, -1.0);
    }

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

    auto computeLocalLocation = [&](int localIdx,
                                    cMode3DMeshFieldEval& threadField) {
        const int globalIdx = locStart + localIdx;

        // Reconstruct the observation position for this point.
        const V3 x0_m = LocationToX0m(prm, globalIdx, nLon, nLat, d_deg, nPtsShell);

        // Compute the minimum cutoff rigidity over all sampled directions.
        // The mesh field is the same for all points (initialized once by
        // InitMeshFields in Mode3D::Run before this parallel section).
        const double rc = ComputeCutoffAtPoint_GV(
            prm, threadField, x0_m, dirs, samplingVertical,
            q_C, m0, box, Rmin, Rmax);

        // localIdx is unique in the OpenMP schedule, so no synchronization is
        // needed for these two stores.
        rcLocal[(size_t)localIdx] = rc;
        if (rc > 0.0) {
            const double pCut = MomentumFromRigidity_GV(rc, qabs);
            eminLocal[(size_t)localIdx] = KineticEnergyFromMomentum_MeV(pCut, m0);
        }

        // Optional directional cutoff map.
        //
        // This loop computes Rc for every requested sky-map direction and stores
        // each directional value, not the minimum over directions.  Therefore the
        // per-location upper-bound optimization used by ComputeCutoffAtPoint_GV
        // must NOT be applied here: a directional map needs the physical cutoff
        // of each cell independently, even if some other direction at the same
        // location has already produced a lower Rc.
        //
        // The same thread-private mesh-field evaluator is reused for all cells of
        // this location.  That preserves the lastNode_ cache across adjacent
        // trajectories and avoids rebuilding any field-management state.
        if (dirMapCfg.enabled) {
            const size_t base = (size_t)localIdx * (size_t)dirMapCfg.nCells;
            for (int cellId=0; cellId<dirMapCfg.nCells; ++cellId) {
                const V3 dir_gsm = DirectionalMapCellDirectionGSM3D(dirMapCfg, cellId);
                dirMapLocal[base + (size_t)cellId] = CutoffForDir_GV(
                    prm, threadField, x0_m, dir_gsm,
                    q_C, m0, box, Rmin, Rmax);
            }
        }
    };

    if (!showProgressBar) {
        // Fast path retained only for source-level symmetry with older versions.
        // It is unreachable in the standalone 3-D cutoff driver because
        // showProgressBar is forced true at the top of this function.
        // Keeping the block avoids a larger refactor while documenting the old
        // no-progress behavior.
        #pragma omp parallel
        {
            // Thread-private mesh field evaluator.
            // Constructed once per thread; lastNode_ starts as nullptr and is
            // updated on every GetB_T call as the particle moves through the mesh.
            cMode3DMeshFieldEval threadField(prm);

            #pragma omp for schedule(dynamic, 1)
            for (int localIdx = 0; localIdx < nLocal; ++localIdx) {
                computeLocalLocation(localIdx,threadField);
            }
        }
    }
    else {
        // Progress path: compute all ranks in the same number of batches.  This
        // guarantees that every rank enters the same MPI_Allreduce sequence,
        // avoiding deadlocks and avoiding MPI calls from OpenMP worker threads.
        const int nProgressBatches = std::max(1,std::min(nLoc,200));

        for (int ibatch=0; ibatch<nProgressBatches; ++ibatch) {
            const int batchBegin = (int)(
                (static_cast<long long>(nLocal) * ibatch) / nProgressBatches);
            const int batchEnd = (int)(
                (static_cast<long long>(nLocal) * (ibatch+1)) / nProgressBatches);

            if (batchEnd > batchBegin) {
                #pragma omp parallel
                {
                    cMode3DMeshFieldEval threadField(prm);

                    #pragma omp for schedule(dynamic, 1)
                    for (int localIdx = batchBegin; localIdx < batchEnd; ++localIdx) {
                        computeLocalLocation(localIdx,threadField);
                    }
                }

                // Update this rank's completed-location counters after the
                // OpenMP workers finish the batch.  These counters are used only
                // for progress reporting; the actual results are already stored
                // in rcLocal/eminLocal.
                const long long nBatchLocations = static_cast<long long>(batchEnd-batchBegin);
                doneLocationsLocal += nBatchLocations;
                doneTasksLocal     += nBatchLocations * tasksPerLocation;

                if (isShells) {
                    for (int localIdx=batchBegin; localIdx<batchEnd; ++localIdx) {
                        const int globalIdx = locStart + localIdx;
                        const int shellIdx = globalIdx / nPtsShell;
                        if (shellIdx>=0 && shellIdx<nShells)
                            locDonePerShellLocal[(size_t)shellIdx]++;
                    }
                }
            }

            // Sum completed work over all ranks.  Only rank 0 renders the line;
            // every rank still participates in the collective so the global
            // number represents the full MPI job, not only the root rank.
            MPI_Allreduce(&doneLocationsLocal,&doneLocationsGlobal,
                          1,MPI_LONG_LONG,MPI_SUM,MPI_GLOBAL_COMMUNICATOR);

            MPI_Allreduce(&doneTasksLocal,&doneTasksGlobal,
                          1,MPI_LONG_LONG,MPI_SUM,MPI_GLOBAL_COMMUNICATOR);

            if (isShells) {
                MPI_Allreduce(locDonePerShellLocal.data(),
                              locDonePerShellGlobal.data(),
                              nShells,MPI_INT,MPI_SUM,MPI_GLOBAL_COMMUNICATOR);
            }

            maybePrintProgress(doneLocationsGlobal,doneTasksGlobal,locDonePerShellGlobal,
                               ibatch == nProgressBatches-1);
        }
    }

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
    MPI_Gather(&nLocal, 1, MPI_INT, recvCounts.data(), 1, MPI_INT, 0, MPI_GLOBAL_COMMUNICATOR);

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
                0, MPI_GLOBAL_COMMUNICATOR);

    MPI_Gatherv(eminLocal.data(), nLocal, MPI_DOUBLE,
                eminAll.data(), recvCounts.data(), displs.data(), MPI_DOUBLE,
                0, MPI_GLOBAL_COMMUNICATOR);

    // Optional directional-map gather.
    //
    // The scalar arrays gather nLocal doubles per rank.  The directional map has
    // nCells doubles per location, so the counts/displacements are the scalar
    // counts/displacements multiplied by nCells.  This reuses the same static
    // location decomposition and keeps rank-0 storage ordered by global location:
    //
    //   dirMapAll[ globalLoc*nCells + cellId ]
    //
    // The gathered array can be large for SHELLS mode, but the allocation happens
    // only when the explicit DIRECTIONAL_MAP=T diagnostic is requested.
    std::vector<double> dirMapAll;
    if (dirMapCfg.enabled) {
        std::vector<int> recvCountsMap(mpiSize, 0);
        std::vector<int> displsMap(mpiSize, 0);

        if (mpiRank == 0) {
            for (int r=0; r<mpiSize; ++r) {
                const long long countMap =
                    static_cast<long long>(recvCounts[(size_t)r]) *
                    static_cast<long long>(dirMapCfg.nCells);
                const long long dispMap =
                    static_cast<long long>(displs[(size_t)r]) *
                    static_cast<long long>(dirMapCfg.nCells);

                if (countMap > static_cast<long long>(std::numeric_limits<int>::max()) ||
                    dispMap   > static_cast<long long>(std::numeric_limits<int>::max())) {
                    throw std::runtime_error(
                        "Mode3D cutoff: directional-map MPI_Gatherv count exceeds INT_MAX. "
                        "Use a coarser DIRMAP grid or fewer locations per run.");
                }

                recvCountsMap[(size_t)r] = static_cast<int>(countMap);
                displsMap[(size_t)r]     = static_cast<int>(dispMap);
            }

            dirMapAll.assign((size_t)nLoc * (size_t)dirMapCfg.nCells, -1.0);
        }

        const long long sendCountMapLL =
            static_cast<long long>(nLocal) * static_cast<long long>(dirMapCfg.nCells);
        if (sendCountMapLL > static_cast<long long>(std::numeric_limits<int>::max())) {
            throw std::runtime_error(
                "Mode3D cutoff: local directional-map MPI_Gatherv count exceeds INT_MAX. "
                "Use a coarser DIRMAP grid or fewer locations per rank.");
        }
        const int sendCountMap = static_cast<int>(sendCountMapLL);

        MPI_Gatherv(dirMapLocal.data(), sendCountMap, MPI_DOUBLE,
                    dirMapAll.data(), recvCountsMap.data(), displsMap.data(), MPI_DOUBLE,
                    0, MPI_GLOBAL_COMMUNICATOR);
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
