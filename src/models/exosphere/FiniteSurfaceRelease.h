#ifndef AMPS_FINITE_SURFACE_RELEASE_H
#define AMPS_FINITE_SURFACE_RELEASE_H

#include <cmath>
#include <stdexcept>
#include <vector>

namespace Exosphere { namespace SourceProcesses {
// Exact depletion of a continuous atom inventory by competing FIRST-ORDER
// processes over an interval with constant coefficients [s^-1]. Temperature
// and illumination are held at the production step's supplied values; their
// temporal-discretization error still requires convergence. This kernel is
// inappropriate for population-independent sputtering or replenishment.
struct FiniteRelease {
  double hazardS=0.0, probability=0.0, remainingAtoms=0.0;
  std::vector<double> atoms;
};
inline FiniteRelease IntegrateFiniteRelease(double populationAtoms,
    const std::vector<double>& hazardsPerS,double dtS) {
  if (!std::isfinite(populationAtoms) || populationAtoms<0.0 ||
      !std::isfinite(dtS) || dtS<=0.0)
    throw std::invalid_argument("finite release needs nonnegative atoms and positive finite seconds");
  FiniteRelease result;
  result.atoms.resize(hazardsPerS.size(),0.0);
  for (double hazard:hazardsPerS) {
    if (!std::isfinite(hazard) || hazard<0.0)
      throw std::invalid_argument("finite-release hazard must be finite and nonnegative [s^-1]");
    result.hazardS+=hazard;
  }
  if (!std::isfinite(result.hazardS))
    throw std::invalid_argument("finite-release total hazard overflow [s^-1]");
  if (result.hazardS==0.0) {result.remainingAtoms=populationAtoms;return result;}
  result.probability=-std::expm1(-result.hazardS*dtS);
  const double released=populationAtoms*result.probability;
  result.remainingAtoms=populationAtoms-released;
  for (std::size_t i=0;i<hazardsPerS.size();++i)
    result.atoms[i]=released*(hazardsPerS[i]/result.hazardS);
  return result;
}
} }
#endif
