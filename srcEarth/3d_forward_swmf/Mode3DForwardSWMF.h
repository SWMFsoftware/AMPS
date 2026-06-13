#ifndef _SRC_EARTH_3D_FORWARD_SWMF_MODE3DFORWARDSWMF_H_
#define _SRC_EARTH_3D_FORWARD_SWMF_MODE3DFORWARDSWMF_H_

#include <string>

//======================================================================================
// Mode3DForwardSWMF.h
//======================================================================================
//
// PURPOSE
// -------
// Thin SWMF-coupled initialization bridge for the 3d_forward energetic-particle
// branch.
//
// LIFECYCLE — TWO-PHASE INITIALIZATION
// -------------------------------------
// The bridge is split into two hooks that must be called in order:
//
//   1. amps_pre_init()  — call BEFORE amps_init_mesh()
//      Sets Earth::ModelMode = BoundaryInjectionMode so that amps_init_mesh()
//      takes the correct MAIN_LIB_GEO path (which calls Earth::Init_AfterParser()),
//      and configures the cDensity3D energy grid before
//      PIC::Mesh::initCellSamplingDataBuffer() sizes the per-cell sampling buffer.
//      Also registers Earth::BoundingBoxInjection::SetPrm so that
//      InitDirectionIMF() (called later by main_lib.cpp::amps_init()) has valid prm.
//
//   2. amps_init()      — call at the END of amps_init() (after amps_init_mesh())
//      Completes mover selection, callback registration, spectrum/weight/time-step
//      setup, and sampler initialization.  Uses the parsed parameters cached by
//      amps_pre_init() — the input file is not re-read.
//
// main_lib.cpp wires both calls:
//   amps_init_mesh():   #if SWMF  →  Mode3DForwardSWMF::amps_pre_init()
//   amps_init() end:   #if SWMF  →  Mode3DForwardSWMF::amps_init()
//
//======================================================================================

namespace Earth {
namespace Mode3DForwardSWMF {

// Phase 1: called BEFORE amps_init_mesh().
// Sets model mode and configures the cDensity3D energy grid.
void amps_pre_init();

// Phase 2: called at the END of amps_init().
// Completes the SWMF-coupled 3d_forward runtime initialization.
void amps_init();

// Return true when AMPS_PARAM.in selects the SWMF-coupled 3-D forward solver.
// This is the historical coupled mode.
bool IsForwardMode();

// Return true when AMPS_PARAM.in selects CALC_TARGET = CUTOFF_RIGIDITY.
// In this mode each amps_time_step() call computes a cutoff-rigidity snapshot
// using the current SWMF-coupled fields and writes timestamped output files.
bool IsCutoffRigidityMode();

// Allocate all AMR leaf blocks on every MPI rank and gather the SWMF-coupled
// cell-centered magnetic field into those blocks.
//
// Why this is needed:
//   In a normal coupled AMPS/SWMF MPI run, only local domain blocks and boundary
//   ghost blocks have node->block allocated on each process.  The SWMF coupler
//   therefore fills B only for that local subset.  Backward cutoff tracing is not
//   local in that sense: one particle trajectory can cross into any AMR block, and
//   the Mode3D mesh-field evaluator expects the requested cell data to exist on
//   the rank performing the tracing.
//
// What this routine does before RunCutoffRigidity():
//   1. gives every used AMR leaf a deterministic dense global ID;
//   2. allocates missing blocks locally on every MPI rank;
//   3. packs owner-rank interior-cell B from PIC::CPLR::SWMF::MagneticFieldOffset;
//   4. uses MPI_Allreduce to replicate the global cell-centered B cache;
//   5. fills all local blocks, including ghost cells, from that replicated cache.
//
// After it returns, the existing Mode3D interpolation code can use the SWMF field
// exactly as if the full mesh had been local from the beginning.  The default
// verbose=true prints one rank-0 diagnostic summary per cutoff snapshot.
void PrepareGlobalSWMFCoupledMagneticFieldForCutoff(bool verbose=true);

// Replace the SWMF-coupled cell-centered magnetic field with the analytic dipole
// configured by DIPOLE_MOMENT and DIPOLE_TILT in AMPS_PARAM.in.  The function
// remains available as a debug override and operates on all currently allocated
// AMR blocks.
void RedefineSWMFCoupledMagneticFieldToAnalyticDipole();

// Per-coupling-call cutoff-rigidity driver used by main_lib.cpp::amps_time_step().
void amps_cutoff_time_step();

// Return a file name that uses the same stamp as the most recent
// SWMF-coupled cutoff-rigidity output.
//
// The stamp is created in amps_cutoff_time_step() from
// PIC::SimulationTime::TimeCounter and passed to
// Earth::Mode3D::SetCutoffOutputFileSuffix() before the cutoff products are
// written.  main_lib.cpp uses this helper immediately afterward to name the
// diagnostic AMPS mesh dump as, for example,
//
//   amps_coupled_data.swmf_n000003_t000600.000s.dat
//
// so it is easy to identify which amps_coupled_data file belongs to which
// cutoff snapshot.
std::string GetLastCutoffOutputFileName(const char* stem,const char* extension);

} // namespace Mode3DForwardSWMF
} // namespace Earth

#endif // _SRC_EARTH_3D_FORWARD_SWMF_MODE3DFORWARDSWMF_H_
