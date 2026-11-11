#ifndef SEP_CORONA_SWCME_TEST_CONTACT_DERIVATIVE_REFERENCE_H
#define SEP_CORONA_SWCME_TEST_CONTACT_DERIVATIVE_REFERENCE_H

#include "sep_coronal_cme/vector_math.h"

#include <array>
#include <cmath>
#include <stdexcept>

namespace SEP { namespace CoronaSwcme { namespace TestReference {

// Independent fourth-order centred derivative at a fixed material label:
// D_h x = [8(x(t+h)-x(t-h))-(x(t+2h)-x(t-2h))]/(12h).
// Positions are metres in the same HCI Cartesian frame and h is seconds, so
// the result is m/s. Pair antisymmetric samples BEFORE multiplying by eight:
// the original expanded sum formed intermediate O(1e9 m) heliocentric terms
// and needlessly rounded away part of the small displacement. This is the
// identical stencil, order, precision and step size, not a changed accuracy
// criterion. It cannot recover rounding already present in sampled positions.
// In particular, it must never substitute the API's reported velocity for an
// independently differentiated position history.
inline double FourthOrderDifference(double minus2,double minus,
    double plus,double plus2,double step) {
  return (8*(plus-minus)-(plus2-minus2))/(12*step);
}

inline CoronalCME::Vec3 FourthOrderDifference(CoronalCME::Vec3 minus2,
    CoronalCME::Vec3 minus,CoronalCME::Vec3 plus,
    CoronalCME::Vec3 plus2,double step) {
  return (8*(plus-minus)-(plus2-minus2))/(12*step);
}

inline bool Finite(CoronalCME::Vec3 value) {
  return std::isfinite(value.x)&&std::isfinite(value.y)&&std::isfinite(value.z);
}

// Signed physical mass flux rho (U - dx_boundary/dt).n, kg/(m^2 s), with the
// normal directed out of the sheath. Callers take the absolute value BEFORE
// area/time integration: opposite leaks cannot cancel. The four positions
// must belong to one supported material chart at t-2h,t-h,t+h,t+2h; support is
// checked by the caller's actual model queries. Invalid input is a test/setup
// rejection, never a zero flux. The normal must already be a unit geometric
// normal; silently normalising an incorrect normal could hide a geometry bug.
inline double IndependentContactFluxKgM2S(
    const std::array<CoronalCME::Vec3,4>& positions,double step,
    double densityKgM3,CoronalCME::Vec3 materialVelocityMPerS,
    CoronalCME::Vec3 outwardNormal) {
  if(!(std::isfinite(step)&&step>0&&std::isfinite(densityKgM3)&&
      densityKgM3>0&&Finite(materialVelocityMPerS)&&Finite(outwardNormal)&&
      std::abs(CoronalCME::Dot(outwardNormal,outwardNormal)-1)<=1e-12))
    throw std::invalid_argument("invalid independent contact-flux inputs");
  for(const auto position:positions)if(!Finite(position))
    throw std::invalid_argument("nonfinite contact position stencil");
  const auto boundaryVelocity=FourthOrderDifference(positions[0],positions[1],
      positions[2],positions[3],step);
  const double flux=densityKgM3*CoronalCME::Dot(
      materialVelocityMPerS-boundaryVelocity,outwardNormal);
  if(!std::isfinite(flux))
    throw std::invalid_argument("nonfinite independent contact mass flux");
  return flux;
}

} } } // namespace SEP::CoronaSwcme::TestReference

#endif
