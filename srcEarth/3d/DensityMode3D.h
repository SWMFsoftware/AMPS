//======================================================================================
// DensityMode3D.h
//======================================================================================
//
// Public interface for the Mode3D mesh-field density + integral-flux calculator.
//
// This module intentionally mirrors the physics of gridless/DensityGridless.cpp:
// it computes energetic-particle density and flux from backward trajectory
// transmissivity, not from forward Monte-Carlo residence-time sampling.  The only
// intended difference is the magnetic-field backend: Mode3D interpolates compact
// global cell-centered arrays with a decomposition-independent row stencil, while
// gridless evaluates Tsyganenko/dipole fields directly at every trajectory step.
// Relativistic conversions, energy and angular grids, quadrature, channel clipping,
// and structured unresolved-trajectory bounds are shared in util/FluxNumerics.h.
//
// Output includes nominal, lower, and upper transmission/spectrum/density/flux values.
// Nominal values use only resolved trajectories and are NaN if none resolves; lower and
// upper values retain a conservative physical interval.  Validation decks can make an
// unresolved fraction above DS_UNRESOLVED_TOL fatal with DS_FAIL_ON_UNRESOLVED=T.
//
//======================================================================================

#ifndef _SRC_EARTH_3D_DENSITYMODE3D_H_
#define _SRC_EARTH_3D_DENSITYMODE3D_H_

#include "../util/amps_param_parser.h"
#include <string>

namespace Earth {
namespace Mode3D {

// Run the mesh-backed density + flux calculation for the currently prepared magnetic
// field snapshot.  The caller is responsible for preparing the AMPS mesh field first:
//
//   standalone Mode3D:
//     Mode3DPrepareMagneticFieldSnapshot(...) in Mode3D.cpp
//
//   SWMF-coupled Mode3D:
//     Mode3DForwardSWMF::PrepareGlobalSWMFCoupledMagneticFieldForCutoff(...)
//
// Returns 0 on success.  Fatal input/geometry errors use the project-standard
// exit(line,file,message) helper so diagnostics match the rest of srcEarth.
int RunDensityAndFlux(const EarthUtil::AmpsParam& prm);

// Install an optional suffix appended immediately before the .dat extension of every
// density/flux file written by RunDensityAndFlux().
//
// Empty suffix preserves the standalone single-snapshot names:
//   mode3d_points_density.dat
//   mode3d_points_spectrum.dat
//   mode3d_points_flux.dat
//
// Non-empty suffix is used for time-series and SWMF-coupled snapshots:
//   mode3d_points_density_snapshot_000001_....dat
//   mode3d_points_density.swmf_n000001_t0003600.000s.dat
void SetDensityOutputFileSuffix(const std::string& suffix);

} // namespace Mode3D
} // namespace Earth

#endif // _SRC_EARTH_3D_DENSITYMODE3D_H_
