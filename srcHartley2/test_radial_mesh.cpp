#include "RadialMesh.h"

#include <algorithm>
#include <cmath>
#include <cstdio>

namespace {

bool Near(double actual,double expected,double relativeTolerance=1.0e-12) {
  return std::abs(actual-expected)<=relativeTolerance*
      std::max(std::abs(actual),std::abs(expected));
}

}

int main() {
  Comet::RadialModel::cLogarithmicResolutionProfile profile;
  const char *error=NULL;

  if (!profile.Init(1000.0,11000.0,100.0,50,&error)) {
    std::fprintf(stderr,"unexpected initialization failure: %s\n",error);
    return 1;
  }
  if (!Near(profile.GetResolution(1000.0),100.0)||
      !Near(profile.GetCellBoundary(50),11000.0)||
      !(profile.GetOuterResolution()>profile.GetResolution(1000.0))) return 2;

  for (int cell=0;cell<49;cell++)
    if (!Near(profile.GetCellWidth(cell+1)/profile.GetCellWidth(cell),
        profile.GetExpansionRatio())) return 3;

  if (!profile.Init(1000.0,11000.0,100.0,100,&error)||
      !Near(profile.GetOuterResolution(),100.0)||
      !Near(profile.GetCellBoundary(100),11000.0)) return 4;

  if (profile.Init(1000.0,11000.0,100.0,101,&error)) return 5;
  if (profile.Init(0.0,11000.0,100.0,50,&error)) return 6;
  if (profile.Init(1000.0,1000.0,100.0,50,&error)) return 7;
  if (profile.Init(1000.0,1200.0,100.0,1,&error)) return 8;
  if (!profile.Init(1000.0,1100.0,100.0,1,&error)) return 9;

  return 0;
}
