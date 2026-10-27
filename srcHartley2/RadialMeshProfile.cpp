#include "RadialMesh.h"

#include <algorithm>
#include <cmath>
#include <limits>

namespace {

double GeometricWidthSum(double ratio,int cellNumber) {
  if (std::abs(ratio-1.0)<1.0e-10) return static_cast<double>(cellNumber);

  const double logRatio=std::log(ratio);
  return std::expm1(cellNumber*logRatio)/std::expm1(logRatio);
}

}

bool Comet::RadialModel::cLogarithmicResolutionProfile::Init(
    double nucleusRadius,double outerRadius,double surfaceResolution,
    int requestedCellNumber,const char **errorMessage) {
  auto fail=[&](const char *message) {
    if (errorMessage!=NULL) *errorMessage=message;
    return false;
  };

  if (!std::isfinite(nucleusRadius)||nucleusRadius<=0.0)
    return fail("the nucleus radius must be finite and positive");
  if (!std::isfinite(outerRadius)||outerRadius<=nucleusRadius)
    return fail("the outer domain radius must exceed the nucleus radius");
  if (!std::isfinite(surfaceResolution)||surfaceResolution<=0.0)
    return fail("the nucleus-surface resolution must be finite and positive");
  if (requestedCellNumber<=0)
    return fail("the requested radial cell number must be positive");

  const double extent=outerRadius-nucleusRadius;
  const double normalizedExtent=extent/surfaceResolution;

  // The first cell has the supplied surface width and subsequent cells can
  // only be the same size or larger.
  if (normalizedExtent+64.0*std::numeric_limits<double>::epsilon()<requestedCellNumber)
    return fail("the requested cell number requires decreasing cell size");

  if (requestedCellNumber==1&&
      std::abs(normalizedExtent-1.0)>
      64.0*std::numeric_limits<double>::epsilon()*std::max(normalizedExtent,1.0))
    return fail("one radial cell requires the domain extent to equal the surface resolution");

  InnerRadius=nucleusRadius;
  OuterRadius=outerRadius;
  InnerResolution=surfaceResolution;
  CellNumber=requestedCellNumber;

  if (std::abs(normalizedExtent-requestedCellNumber)<=
      64.0*std::numeric_limits<double>::epsilon()*
      std::max(normalizedExtent,static_cast<double>(requestedCellNumber))) {
    ExpansionRatio=1.0;
    return true;
  }

  double lower=1.0,upper=2.0;
  while (GeometricWidthSum(upper,requestedCellNumber)<normalizedExtent) {
    upper*=2.0;
    if (!std::isfinite(upper)||upper>1.0e100)
      return fail("the requested logarithmic expansion is outside the representable range");
  }

  for (int iteration=0;iteration<100;iteration++) {
    const double middle=0.5*(lower+upper);
    if (GeometricWidthSum(middle,requestedCellNumber)<normalizedExtent)
      lower=middle;
    else
      upper=middle;
  }

  ExpansionRatio=0.5*(lower+upper);
  return true;
}

double Comet::RadialModel::cLogarithmicResolutionProfile::GetResolution(
    double radius) const {
  if (radius<=InnerRadius||ExpansionRatio==1.0) return InnerResolution;
  if (radius>=OuterRadius) return GetOuterResolution();

  // Determine which cell of the ideal geometric mesh contains the requested
  // radius. AMPS then approximates this piecewise-constant target with its
  // factor-of-two block refinement levels.
  const double normalizedRadius=(radius-InnerRadius)/InnerResolution;
  int cell=static_cast<int>(std::floor(
      std::log1p((ExpansionRatio-1.0)*normalizedRadius)/
      std::log(ExpansionRatio)));

  if (cell<0) cell=0;
  if (cell>=CellNumber) cell=CellNumber-1;
  return GetCellWidth(cell);
}

double Comet::RadialModel::cLogarithmicResolutionProfile::GetOuterResolution() const {
  return GetCellWidth(CellNumber-1);
}

double Comet::RadialModel::cLogarithmicResolutionProfile::GetCellWidth(int cell) const {
  if (cell<0||cell>=CellNumber) return -1.0;
  return InnerResolution*std::pow(ExpansionRatio,cell);
}

double Comet::RadialModel::cLogarithmicResolutionProfile::GetCellBoundary(
    int boundary) const {
  if (boundary<0||boundary>CellNumber) return -1.0;
  if (ExpansionRatio==1.0)
    return InnerRadius+boundary*InnerResolution;

  return InnerRadius+InnerResolution*
      std::expm1(boundary*std::log(ExpansionRatio))/(ExpansionRatio-1.0);
}
