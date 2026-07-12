//======================================================================================
// CutoffRigidityMode3D.h
//======================================================================================
//
// Public interface for mesh-backed SEP/GCR cutoff, transmissivity, density, and flux
// calculations.
//
// Mode3D uses the normal distributed AMPS block decomposition.  Before backtracing,
// owner-rank cell-centered B/E values are assembled into compact global arrays by
// GlobalMagneticField.cpp.  The replicated AMR tree is retained, and each used leaf is
// assigned a deterministic node->Temp_ID.  Field evaluation builds an AMPS
// decomposition-independent cRowStencil and applies its (node,i,j,k,weight) entries to
// the compact arrays.  Remote node->block objects are therefore not required.
//
// MPI ranks dynamically share observation locations.  OpenMP or std::thread workers
// operate within each rank.  During trajectory integration the tree and global fields
// are read-only and no inter-rank field communication is performed.
//======================================================================================

#ifndef _SRC_EARTH_3D_CUTOFFRIGIDITYMODE3D_H_
#define _SRC_EARTH_3D_CUTOFFRIGIDITYMODE3D_H_

#include "../util/amps_param_parser.h"
#include "../gridless/CutoffRigidityGridless.h"
#include <string>

namespace Earth {
namespace Mode3D {

//--------------------------------------------------------------------------------------
// RunCutoffRigidity
//
// Top-level entry point for the Mode3D cutoff rigidity workflow.
// Must be called after mesh initialization and compact global-field assembly.
//
// Returns 0 on success.
// Throws std::runtime_error on invalid input or runtime failures.
//--------------------------------------------------------------------------------------
// If showProgressBar is true, rank 0 prints a time-throttled global progress
// bar while all MPI ranks process synchronized work batches.  The user-facing
// standalone 3-D cutoff path defaults to progress reporting; the argument is retained for API compatibility, but the current standalone
// implementation forces the progress path on internally.
int RunCutoffRigidity(const EarthUtil::AmpsParam& prm, bool requestedProgressBar=true);

// Set an optional suffix appended to all RunCutoffRigidity output files.
//
// Default/standalone behavior is unchanged when the suffix is empty:
//   cutoff_3d_points.dat
//   cutoff_3d_shells.dat
//
// SWMF-coupled runs call amps_time_step() multiple times for successive MHD
// snapshots.  In that case the coupling bridge sets a suffix such as
//   .swmf_n000003_t000600.000s
// before each cutoff calculation so every snapshot writes a distinct file.
void SetCutoffOutputFileSuffix(const std::string& suffix);


// Mesh-backed trajectory classifier shared with Mode3D density/flux.
//
// These functions expose the exact same backward-tracing kernel used internally by
// RunCutoffRigidity(), but keep the public interface close to the gridless shared
// classifier so density/flux code can switch field backends with minimal changes.
//
// Inputs:
//   x0_m      — starting point in GSM [m]
//   v0_unit   — initial BACKTRACED velocity direction, unit vector in GSM
//   R_GV      — rigidity [GV]
//   maxTraceTimeOverride_s > 0 overrides #CUTOFF_RIGIDITY/#NUMERICAL time caps
//
// Return value:
//   true  = trajectory escaped the outer Mode3D domain before hitting the loss sphere
//   false = trajectory hit the inner sphere or exceeded integration limits
//
// TraceAllowedMeshEx additionally fills GridlessMode::TrajectoryExitState for allowed
// trajectories; this is needed by DS_BOUNDARY_MODE=ANISOTROPIC so the boundary PAD and
// spatial weighting can be evaluated at the asymptotic exit location/direction.
bool TraceAllowedMesh(const EarthUtil::AmpsParam& prm,
                      const double x0_m[3],
                      const double v0_unit[3],
                      double R_GV,
                      double maxTraceTimeOverride_s=-1.0);

bool TraceAllowedMeshEx(const EarthUtil::AmpsParam& prm,
                        const double x0_m[3],
                        const double v0_unit[3],
                        double R_GV,
                        Earth::GridlessMode::TrajectoryExitState* exitState,
                        double maxTraceTimeOverride_s=-1.0);


} // namespace Mode3D
} // namespace Earth

#endif // _SRC_EARTH_3D_CUTOFFRIGIDITYMODE3D_H_
