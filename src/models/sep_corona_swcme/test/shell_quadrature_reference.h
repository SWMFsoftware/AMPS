#ifndef SEP_CORONA_SWCME_TEST_SHELL_QUADRATURE_REFERENCE_H
#define SEP_CORONA_SWCME_TEST_SHELL_QUADRATURE_REFERENCE_H

#include <cmath>
#include <stdexcept>

namespace SEP { namespace CoronaSwcme { namespace TestReference {

// Independent exact polynomial quadrature, not the production factored shell
// formula. Two Gauss--Legendre nodes integrate both r^2 and r^3 exactly, so
// <r>=integral(r^3 dr)/integral(r^2 dr) is the physical volume centroid in m.
// Extended accumulation avoids the old oracle's subtraction of nearly equal
// heliocentric fourth/third powers. Those cancelled digits had contaminated
// the initial B_r*r^2 reference when production geometry became more accurate.
// Neither the radial-flux tolerance nor a physical input is changed here.
inline double ShellCentroidByGaussianMoments(double left,double right) {
  if(!(std::isfinite(left)&&std::isfinite(right)&&left>=0&&right>left))
    throw std::invalid_argument("invalid independent shell reference faces");
  const long double lo=left,hi=right,width=hi-lo,mid=lo+width/2;
  const long double offset=width/(2*std::sqrt(3.0L));
  const long double a=mid-offset,b=mid+offset;
  return static_cast<double>((a*a*a+b*b*b)/(a*a+b*b));
}

} } } // namespace SEP::CoronaSwcme::TestReference

#endif
