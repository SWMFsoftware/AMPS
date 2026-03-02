//======================================================================================
// DensityGridless.h
//======================================================================================
// GRIDLESS DENSITY + SPECTRUM SOLVER (POINTS MODE)
//
// This module is intentionally separate from CutoffRigidityGridless.*.
// It reuses the *same physical classification concept* (allowed vs forbidden
// backtracing in an analytic magnetospheric field), but the computational goal
// is different:
//
//   - Cutoff rigidity: find a single cutoff value Rc per location.
//   - Density/spectrum: evaluate transmissivity T(E) over an energy grid,
//     fold it with a boundary spectrum J_b(E), and integrate to obtain a
//     total number density n.
//
// Entry point below is called from srcEarth/main.cpp when:
//   -mode gridless
//   and AMPS_PARAM.in contains:
//     #CALCULATION_MODE
//       CALC_TARGET       DENSITY_SPECTRUM
//     #OUTPUT_DOMAIN
//       OUTPUT_MODE       POINTS
//     #DENSITY_SPECTRUM
//       DS_* parameters
//     #SPECTRUM
//       SPECTRUM_TYPE + parameters
//
// Output:
//   (1) gridless_points_density.dat  : one line per point (density)
//   (2) gridless_points_spectrum.dat : one Tecplot ZONE per point (spectrum)
//======================================================================================

#ifndef _SRC_EARTH_GRIDLESS_DENSITYGRIDLESS_H_
#define _SRC_EARTH_GRIDLESS_DENSITYGRIDLESS_H_

#include "util/amps_param_parser.h"

namespace Earth {
  namespace GridlessMode {
    // Returns 0 on success; throws std::runtime_error on invalid input.
    int RunDensityAndSpectrumPoints(const EarthUtil::AmpsParam& p);
  }
}

#endif
