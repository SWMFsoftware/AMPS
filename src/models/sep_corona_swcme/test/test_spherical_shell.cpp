#include "sep_corona_swcme/spherical_shell.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <utility>

namespace {

void Require(bool condition,const char* message) {
  if(!condition)throw std::runtime_error(message);
}

void Check(double left,double right) {
  // Independent two-node Gauss--Legendre integration. This rule is exact for
  // degree <=3, so it integrates both r^2 and r^3 moments without using the
  // production factored formula or subtracting nearly equal powers. Extended
  // accumulation isolates the double-precision implementation's rounding.
  // All comparisons use the ACTUAL representable input faces, not an intended
  // width that may have rounded when added to a large heliocentric radius.
  const long double lo=left,hi=right,width=hi-lo,mid=lo+width/2;
  const long double offset=width/(2*std::sqrt(3.0L));
  const long double a=mid-offset,b=mid+offset;
  const long double volume=width*(a*a+b*b)/2;
  const long double firstMoment=width*(a*a*a+b*b*b)/2;
  const long double centroid=firstMoment/volume;
  const double gotVolume=SEP::CoronaSwcme::SphericalShellVolumePerSrM3(left,right);
  const double gotCentroid=SEP::CoronaSwcme::SphericalShellCentroidM(left,right);
  // Eight eps is the preregistered roundoff budget for the small number of
  // positive products/sums/divisions, not a fitted physical tolerance. No
  // discretization error exists in these polynomial reference moments.
  constexpr long double relative=8*std::numeric_limits<double>::epsilon();
  Require(std::isfinite(gotVolume)&&gotVolume>0,
      "positive physical shell lost its finite volume");
  Require(std::abs(gotVolume-volume)<=relative*volume,
      "shell volume disagrees with independent exact Gaussian moment");
  Require(std::isfinite(gotCentroid)&&gotCentroid>=left&&gotCentroid<=right,
      "thin shell centroid left its physical faces");
  Require(std::abs(gotCentroid-centroid)<=relative*centroid,
      "shell centroid disagrees with independent exact Gaussian moment");
}

} // namespace

int main() {
  try {
    int shells=0;
    Check(0,1);++shells;
    for(double left:{1e-6,1.0,695700000.0,149597870700.0}) {
      for(double relativeWidth:{1.0,1e-3,1e-9,1e-13}) {
        Check(left,left*(1+relativeWidth));++shells;
      }
      Check(left,std::nextafter(left,std::numeric_limits<double>::infinity()));
      ++shells;
    }
    const double nan=std::numeric_limits<double>::quiet_NaN();
    for(const auto faces:{std::pair<double,double>{-1,2},{1,1},{2,1},
                         {nan,2},{1,nan},{1e200,2e200}}) {
      Require(std::isnan(SEP::CoronaSwcme::SphericalShellVolumePerSrM3(
          faces.first,faces.second)),"invalid/unrepresentable volume was accepted");
      Require(std::isnan(SEP::CoronaSwcme::SphericalShellCentroidM(
          faces.first,faces.second)),"invalid/unrepresentable centroid was accepted");
    }
    std::cout<<"[CMBGU05-SHELL-GEOMETRY-CONTROLS] PASS "<<shells
             <<" independent broad/thin/one-ulp shell moments and invalid-domain checks\n";
    return 0;
  }catch(const std::exception& error) {
    std::cerr<<"CMBGU05-SHELL-GEOMETRY-CONTROLS FAIL: "<<error.what()<<'\n';
    return 1;
  }
}
