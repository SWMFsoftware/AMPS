#include "RadiativeCooling.h"

#include <cassert>
#include <cmath>

namespace {
bool Near(double a,double b) {
  return std::fabs(a-b)<=1.0e-12*(1.0+std::fabs(a)+std::fabs(b));
}
}

int main() {
  using namespace CometaryPhysics::RadiativeCooling;

  // With zero radius the optical-depth factor is one.  Check both published
  // temperature branches after the erg-to-joule conversion.
  assert(Near(Gombosi1986WaterEnergyLossRate(50.0,1.0e18,0.0),
      4.4e-29*std::pow(50.0,3.35)));
  assert(Near(Gombosi1986WaterEnergyLossRate(52.0,1.0e18,0.0),
      2.0e-27*std::pow(52.0,2.47)));

  // n=1e15 m^-3 and r=1 km give tau=0.16 for sigma_IR=4e-15 cm^2.
  const double opticallyThin=Gombosi1986WaterEnergyLossRate(
      100.0,1.0e15,0.0);
  const double attenuated=Gombosi1986WaterEnergyLossRate(
      100.0,1.0e15,1000.0);
  assert(Near(attenuated/opticallyThin,std::exp(-0.16)));

  assert(Gombosi1986WaterEnergyLossRate(0.0,1.0e18,1000.0)==0.0);
  assert(Gombosi1986WaterEnergyLossRate(100.0,0.0,1000.0)==0.0);
  assert(EnergyLossRate(Model::Disabled,100.0,1.0e18,1000.0)==0.0);
  return 0;
}
