#ifndef SEP_COMMON_SEP_INJECTION_SPECTRUM_H
#define SEP_COMMON_SEP_INJECTION_SPECTRUM_H

#include "sep_transport_common.h"

#include <cstdint>
#include <string>

namespace SEP {
namespace Injection {

// The exponent is always defined for the probability density in the named
// coordinate.  This removes the historical ambiguity between phase-space
// density f(p), dN/dp, and dN/dlog(p).
enum class Measure { Momentum, LogMomentum, KineticEnergy, LogKineticEnergy };
enum class AngularDistribution { Isotropic, FieldAlignedPlus, FieldAlignedMinus };

struct Spectrum {
  Measure measure = Measure::Momentum;
  double minimum = 0.0;
  double maximum = 0.0;
  double powerIndex = 0.0;
};

struct Configuration {
  std::uint64_t campaignSeed = 0;
  std::uint64_t macroparticlesPerEvent = 0;
  double injectionEfficiency = 0.0;
  Spectrum spectrum;
  AngularDistribution angular = AngularDistribution::Isotropic;
};

Transport::Status Validate(const Configuration& configuration);
Transport::ScalarResult ProbabilityDensity(const Spectrum& spectrum,
                                           double coordinate);
Transport::ScalarResult InverseCdf(const Spectrum& spectrum, double uOpen01);
Transport::ScalarResult ImportanceWeightForLogUniformProposal(
    const Spectrum& target, double coordinate);

// Purpose-separated source keys are position sensitive: swapping field-line
// and species identifiers must change the stream.  Semantic tags also prevent
// adding a new random draw for direction from perturbing the energy sequence.
enum class RandomPurpose : std::uint64_t {
  EventCount = 1,
  Spectrum = 2,
  PitchAngle = 3,
  Gyrophase = 4,
  Position = 5
};

struct RandomKey {
  std::uint64_t campaign = 0;
  std::uint64_t event = 0;
  std::uint64_t fieldLine = 0;
  std::uint64_t species = 0;
  std::uint64_t macroparticle = 0;
  RandomPurpose purpose = RandomPurpose::Spectrum;
};

std::uint64_t HashRandomKey(const RandomKey& key);
Transport::KeyedRandomStream MakeRandomStream(const RandomKey& key);
std::string Fingerprint(const Configuration& configuration);

}  // namespace Injection
}  // namespace SEP

#endif  // SEP_COMMON_SEP_INJECTION_SPECTRUM_H
