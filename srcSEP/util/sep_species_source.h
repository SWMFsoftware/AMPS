#ifndef SEP_UTIL_SEP_SPECIES_SOURCE_H
#define SEP_UTIL_SEP_SPECIES_SOURCE_H

#include "sep_coefficient_physics.h"

#include <string>
#include <vector>

namespace SEP {
namespace Transport {
namespace SpeciesSource {

// Source spectra must state whether an energy coordinate means total kinetic
// energy or kinetic energy per nucleon.  Electrons may use only TotalKinetic;
// PerNucleon requires a strictly positive nucleon count in SpeciesProperties.
enum class EnergyConvention { TotalKinetic, PerNucleon };

const char* EnergyConventionName(EnergyConvention convention);

// This record joins the transport species definition to its source model.
// abundanceFraction is normalized across the complete configured source;
// injectionEfficiency and spectralIndex remain species-local physical inputs.
struct Configuration {
  CoefficientPhysics::SpeciesProperties species;
  double abundanceFraction = 1.0;
  double injectionEfficiency = 1.0;
  double spectralIndex = 4.0;
  EnergyConvention energyConvention = EnergyConvention::TotalKinetic;
};

// Validate every species before normalization, reject duplicate PIC species
// identifiers, and return fractions whose sum is unity.  The input is never
// mutated, which makes failed initialization transactional.
Status ValidateAndNormalize(const std::vector<Configuration>& requested,
                            std::vector<Configuration>* normalized);

// Install an immutable run-wide source table during initialization.  An empty
// table intentionally selects the legacy source parameters; once a table is
// installed, every injected PIC species must have an explicit entry.
Status SetActiveConfiguration(
    const std::vector<Configuration>& requested);
bool HasActiveConfiguration();
Status FindActiveConfiguration(int modelSpecies, Configuration* result);

// Convert a configured energy coordinate [J] to total kinetic energy [J].
// This is the only conversion that injection kernels should pass to the
// relativistic momentum routines, preventing MeV/nucleon from being treated as
// total MeV for alpha particles and heavy ions.
ScalarResult TotalKineticEnergyJ(const Configuration& configuration,
                                 double configuredEnergyJ);

}  // namespace SpeciesSource
}  // namespace Transport
}  // namespace SEP

#endif  // SEP_UTIL_SEP_SPECIES_SOURCE_H
