#ifndef SRCMOON_SPECIES_SODIUM_CASE1_H
#define SRCMOON_SPECIES_SODIUM_CASE1_H

#include <string>

class cInternalSphericalData;

namespace Moon { namespace Runtime { struct Configuration; } }

namespace Moon {
namespace Species {
namespace SodiumCase1 {

// Validate and install the immutable Tenishev et al. (2013) Case-1 runtime
// receipt after the logical lunar surface has been allocated.  The routine
// initializes the prescribed reservoir in atoms and fails if the executable's
// generated source switches cannot represent the published model.
bool ConfigureProduction(const Runtime::Configuration& configuration,
                         cInternalSphericalData* logicalSurface,
                         std::string* error);

// Ordinary Exosphere runtime callbacks.  Rates are whole-Moon atoms s^-1,
// temperatures are kelvin, and surface densities are atoms m^-2.
double ImpactVaporizationRate(int species);
double ImpactVaporizationTemperatureK(int species);
bool PrescribedSurfaceDensityM2(int species,int surfaceElement,
                                double* numberDensityM2);

// Illuminated neutral-Na lifetime in seconds.  When the paper mode is not
// enabled this returns the historical Huebner-scaled srcMoon value so U14 and
// the no-input regression retain their established behavior.
double PhotoionizationLifetimeS();

bool ProductionEnabled();

}  // namespace SodiumCase1
}  // namespace Species
}  // namespace Moon

#endif  // SRCMOON_SPECIES_SODIUM_CASE1_H
