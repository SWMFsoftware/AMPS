#ifndef _SRC_EARTH_3D_FORWARD_SWMF_MODE3DFORWARDSWMF_H_
#define _SRC_EARTH_3D_FORWARD_SWMF_MODE3DFORWARDSWMF_H_

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

// Per-coupling-call cutoff-rigidity driver used by main_lib.cpp::amps_time_step().
void amps_cutoff_time_step();

} // namespace Mode3DForwardSWMF
} // namespace Earth

#endif // _SRC_EARTH_3D_FORWARD_SWMF_MODE3DFORWARDSWMF_H_
