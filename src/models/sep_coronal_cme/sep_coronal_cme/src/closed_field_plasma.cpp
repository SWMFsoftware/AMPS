#include "sep_coronal_cme/closed_field_plasma.h"
#include "sep_coronal_cme/constants.h"

#include <algorithm>
#include <cmath>

namespace SEP { namespace CoronalCME {
double RotatingEffectivePotential(Vec3 x,Vec3 omega) {
  const double r=Norm(x);if(!(r>0))return std::numeric_limits<double>::quiet_NaN();
  return -Constants::kSolarGravitationalParameterM3PerS2/r-
      0.5*Dot(Cross(omega,x),Cross(omega,x));
}
Core::Result<HydrostaticState> IsothermalHydrostatic(double rho0,double p0,double delta) {
  if(!(rho0>0&&p0>0&&std::isfinite(delta)))return Core::Result<HydrostaticState>::Failure(
      Core::StatusCode::InvalidConfiguration,"isothermal hydrostatic base state is invalid");
  const double c2=p0/rho0,scale=std::exp(-delta/c2);
  return Core::Result<HydrostaticState>::Success({rho0*scale,p0*scale,1.0});
}
Core::Result<HydrostaticState> PolytropicHydrostatic(double rho0,double p0,double gamma,double delta) {
  if(!(rho0>0&&p0>0&&gamma>1&&std::isfinite(delta)))return Core::Result<HydrostaticState>::Failure(
      Core::StatusCode::InvalidConfiguration,"polytropic hydrostatic base state is invalid");
  const double h0=gamma*p0/((gamma-1)*rho0),h=h0-delta;
  if(!(h>0))return Core::Result<HydrostaticState>::Failure(Core::StatusCode::OutOfDomain,
      "polytropic hydrostatic enthalpy is nonpositive");
  const double q=h/h0;
  return Core::Result<HydrostaticState>::Success({rho0*std::pow(q,1/(gamma-1)),
      p0*std::pow(q,gamma/(gamma-1)),std::pow(q,1.0)});
}
Core::Status CheckFootpointCompatibility(double a,double b,double tolerance) {
  if(!(a>0&&b>0&&tolerance>=0))return Core::Status::Failure(
      Core::StatusCode::InvalidConfiguration,"footpoint comparison inputs are invalid");
  if(std::abs(a-b)>tolerance*std::max(a,b))return Core::Status::Failure(
      Core::StatusCode::DataIntegrityFailure,"closed-loop footpoint pressures are incompatible");
  return Core::Status::Success();
}
} }
