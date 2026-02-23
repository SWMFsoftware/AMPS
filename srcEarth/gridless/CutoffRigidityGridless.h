//======================================================================================
// CutoffRigidityGridless.h
//======================================================================================
// PURPOSE
//   Gridless cutoff rigidity solver for srcEarth.
//
//   This module implements a "gridless" pipeline that directly evaluates the
//   Tsyganenko magnetic field models (T96 or T05) and IGRF, and traces test
//   particles to determine cutoff rigidity. It is designed to be called from
//   srcEarth/main.cpp when the CLI option "-mode gridless" is used.
//
// ALGORITHM (PROTOTYPE)
//   - For each injection location (POINTS or SHELLS grid), sample a fixed
//     angular grid of directions.
//   - For each direction, find the minimum rigidity Rc (GV) that allows a
//     particle to escape the domain boundary (rectangular box) without hitting
//     the inner loss sphere (R_INNER).
//   - The location cutoff rigidity is the minimum Rc over directions.
//
// UNITS CONVENTION (IMPORTANT)
//   The RoR input files are interpreted as providing *geometric distances in km*:
//     - DOMAIN_X/Y/Z_MIN/MAX  [km]
//     - R_INNER               [km]
//     - POINT coordinates     [km]
//   Internally we convert km -> meters, and then to Re where required by
//   the Tsyganenko Fortran interfaces.
//
// OUTPUT
//   Tecplot ASCII files:
//     cutoff_gridless_points.dat
//     cutoff_gridless_shells.dat
//
//======================================================================================

#ifndef _SRC_EARTH_GRIDLESSM_CUTOFFRIGIDITYGRIDLESS_H_
#define _SRC_EARTH_GRIDLESSM_CUTOFFRIGIDITYGRIDLESS_H_

#include "util/amps_param_parser.h"

namespace Earth {
  namespace GridlessMode {
    int RunCutoffRigidity(const EarthUtil::AmpsParam& p);
  }
}

#endif
