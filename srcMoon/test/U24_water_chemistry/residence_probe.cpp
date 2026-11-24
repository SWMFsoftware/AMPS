// Exercise the same scalar rate and finite-time desorption probability used
// by the lunar surface inventory. Compiled species and sticking tables are
// not consulted by either first-order kernel. This conditional calculation
// starts after adsorption and never labels a sodium or argon particle as OH.
#include "MoonSurfacePhysics.h"

#include <iomanip>
#include <iostream>

int main() {
  double fraction, attemptHz, barrierJ, temperatureK, stepS;
  while (std::cin >> fraction >> attemptHz >> barrierJ >> temperatureK >> stepS) {
    Moon::SurfacePhysics::SpeciesParameters parameters;
    parameters.residenceLaw = Moon::SurfacePhysics::ResidenceLaw::FirstOrderArrhenius;
    parameters.attemptFrequencyHz = fraction * attemptHz;
    parameters.activationEnergyJ = barrierJ;
    std::cout << std::setprecision(17)
              << Moon::SurfacePhysics::ResidenceTimeS(parameters, temperatureK) << ' '
              << Moon::SurfacePhysics::DesorptionProbability(parameters, temperatureK, stepS) << '\n';
  }
  return std::cin.eof() ? 0 : 2;
}
