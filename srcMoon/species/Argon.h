#ifndef SRCMOON_SPECIES_ARGON_H
#define SRCMOON_SPECIES_ARGON_H

#include <string>
#include <vector>

namespace Moon {
namespace Species {
namespace Argon {

// Species-specific source parameters are kept out of the generic surface and
// thermal services.  Rates are physical 40Ar atoms s^-1 over the whole Moon;
// element areas are m^2.  The rate partition is dimensionless and therefore
// remains valid for either the logical sphere or its LOLA radial adapter.
struct SourceParameters {
  bool enabled = false;
  std::string spatialModel;  // Currently only "uniform-global".
  std::string timeModel;     // Only literature-supported "continuous".
  double totalRateAtomsS = 0.0;
  std::string reference;
};

// Validate the source receipt without consulting AMPS global state.  This is
// used by both the production adapter and U21's independent invariant checks.
bool ValidateSourceParameters(const SourceParameters& parameters,
                              std::string* error);

// Partition a global uniform source in proportion to reference-sphere area.
// The resulting rates sum to totalRateAtomsS to floating-point roundoff.  No
// equal-element assumption is made because the logical latitude/longitude
// elements have latitude-dependent areas.  This represents uniform solid-
// angle geography; it does not turn unresolved terrain slope area into an
// unreviewed source enhancement.
bool BuildUniformElementRates(const std::vector<double>& elementAreasM2,
                              double totalRateAtomsS,
                              std::vector<double>* elementRatesAtomsS,
                              std::string* error);

}  // namespace Argon
}  // namespace Species
}  // namespace Moon

#endif  // SRCMOON_SPECIES_ARGON_H
