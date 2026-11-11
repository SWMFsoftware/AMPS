#ifndef SEP_CORONA_SWCME_SPHERICAL_SHELL_H
#define SEP_CORONA_SWCME_SPHERICAL_SHELL_H

#include <cmath>
#include <limits>

namespace SEP { namespace CoronaSwcme {

// Exact spherical-sector moments per unit solid angle. Radii are HCI metres
// and satisfy 0 <= left < right. Factoring r_R^3-r_L^3 removes subtraction of
// nearly equal O(r^3) values in compressed material shells. The result is m^3
// per steradian; callers multiply by the immutable physical ray weight.
// No thin-shell approximation, minimum width or geometric floor is used.
// Invalid/overflowed geometry returns NaN, which the existing solver/provider
// checks reject through their typed numerical/out-of-domain status paths.
inline double SphericalShellVolumePerSrM3(double left,double right) {
  if(!(std::isfinite(left)&&std::isfinite(right)&&left>=0&&right>left))
    return std::numeric_limits<double>::quiet_NaN();
  const double volume=(right-left)*(right*right+right*left+left*left)/3;
  return std::isfinite(volume)&&volume>0?volume:
      std::numeric_limits<double>::quiet_NaN();
}

// Volume-weighted radius <r> = integral r^3 dr / integral r^2 dr,
// exactly (3/4)(r_R^4-r_L^4)/(r_R^3-r_L^3). Cancel the common positive width
// symbolically before floating-point evaluation. Evaluate the centroid as a
// positive offset from the inner face so a shell only a few ulps thick cannot
// acquire a centre outside its faces from multiplying large heliocentric
// radii. For w=r_R-r_L the offset is
// w*(r_L^2/2+2*r_L*w/3+w^2/4)/(r_L^2+r_L*w+w^2/3).
// Besides avoiding four generic
// power calls in a frequently evaluated RK-stage geometry, this retains the
// centre of very thin shells whose expanded power differences lose digits.
// The radial finite-volume equations, cell masses, RK4/CFL controls, source
// samples, viscosity and energy accounting are unchanged. Arithmetic rounding
// can change, so independent moments and all physical regressions must pass.
inline double SphericalShellCentroidM(double left,double right) {
  if(!(std::isfinite(left)&&std::isfinite(right)&&left>=0&&right>left))
    return std::numeric_limits<double>::quiet_NaN();
  const double width=right-left;
  const double radius=left+width*(0.5*left*left+
      (2.0/3)*left*width+0.25*width*width)/
      (left*left+left*width+width*width/3);
  return std::isfinite(radius)?radius:
      std::numeric_limits<double>::quiet_NaN();
}

} } // namespace SEP::CoronaSwcme

#endif
