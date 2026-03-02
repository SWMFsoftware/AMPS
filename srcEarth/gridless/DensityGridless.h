//======================================================================================
// DensityGridless.h
//======================================================================================
/**
 * \file DensityGridless.h
 *
 * Public interface for the GRIDLESS energetic-particle density + spectrum calculator.
 *
 * WHY THIS MODULE EXISTS
 * ----------------------
 * The gridless cutoff-rigidity tool determines, for a given location and direction, whether
 * a particle of a given rigidity can access the location from outside the magnetosphere
 * (ALLOWED) or is blocked/absorbed (FORBIDDEN). That is a *geometric/topological* classifier.
 *
 * For many applications we also need *physics outputs* such as:
 *   - local energetic-particle spectrum at a point inside the magnetosphere, and
 *   - local total number density integrated over an energy band.
 *
 * This module reuses exactly the same backtracing classifier as cutoff rigidity, but instead
 * of searching for a single cutoff rigidity, it evaluates the allowed fraction T(E) across
 * an energy grid and folds it with a boundary spectrum J_b(E).
 *
 * THEORY SUMMARY
 * --------------
 * For each observation point x0 and kinetic energy E:
 *   (1) Backtrace a set of directions and compute transmissivity:
 *         T(E;x0) = N_allowed / N_dirs
 *   (2) Local differential intensity (isotropic assumption):
 *         J_loc(E;x0) = T(E;x0) * J_b(E)
 *   (3) Total number density:
 *         n_tot(x0) = 4*pi * integral_{Emin}^{Emax} [ J_loc(E;x0) / v(E) ] dE
 *
 * where v(E) is the relativistic particle speed at kinetic energy E.
 *
 * INPUT CONTRACT
 * --------------
 * This solver requires:
 *   CALC_TARGET            DENSITY_SPECTRUM
 *   FIELD_EVAL_METHOD      GRIDLESS
 *   OUTPUT_MODE            POINTS
 * plus a #DENSITY_SPECTRUM section defining the energy grid and optional work caps.
 *
 * OUTPUT CONTRACT
 * ---------------
 * Writes two Tecplot files (rank 0):
 *   - gridless_points_density.dat
 *   - gridless_points_spectrum.dat
 *
 * See DensityGridless.cpp for full details and derivations.
 */

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
