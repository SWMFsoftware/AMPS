#ifndef SEP_COEFFICIENT_PROVIDERS_H
#define SEP_COEFFICIENT_PROVIDERS_H

#include "transport_common.h"

#include "util/sep_coefficient_registry.h"

#include <cstdint>
#include <string>

namespace SEP {
namespace Transport {
namespace PICAdapter {

// Runtime counters make the invalid-coefficient policy observable.  In
// particular, selecting "ballistic" can no longer silently clamp lambda to an
// arbitrary minimum or maximum as the legacy event mover did.
struct CoefficientProviderDiagnostics {
  std::uint64_t invalidSamples = 0;
  std::uint64_t ballisticSubstitutions = 0;
};

CoefficientProviderDiagnostics GetCoefficientProviderDiagnostics();
void ResetCoefficientProviderDiagnostics();

// These adapters are the only production bridge from the pure SI coefficient
// interfaces to legacy PIC field-line storage.  All movers therefore share
// provider selection, validation, provenance, and conversion policy.
class PICSpatialDiffusionProvider : public SpatialDiffusionProvider {
 public:
  PICSpatialDiffusionProvider(const ParticleContext& context,
                              const std::string& turbulenceIdentity);
  SpatialDiffusionSample Evaluate(double sM,
                                  double speedMPerS) const override;

 private:
  const ParticleContext& context_;
  std::string turbulenceIdentity_;
};

class PICPitchAngleDiffusionProvider : public PitchAngleDiffusionProvider {
 public:
  PICPitchAngleDiffusionProvider(const ParticleContext& context,
                                 const std::string& turbulenceIdentity);
  PitchAngleDiffusionSample Evaluate(double sM, double momentumKgMPerS,
                                     double mu) const override;

 private:
  const ParticleContext& context_;
  std::string turbulenceIdentity_;
};

class PICMeanFreePathProvider : public MeanFreePathProvider {
 public:
  PICMeanFreePathProvider(const ParticleContext& context,
                          const std::string& turbulenceIdentity);
  MeanFreePathSample Evaluate(double sM, double momentumKgMPerS,
                              double mu) const override;

 private:
  const ParticleContext& context_;
  std::string turbulenceIdentity_;
};

}  // namespace PICAdapter
}  // namespace Transport
}  // namespace SEP

#endif  // SEP_COEFFICIENT_PROVIDERS_H
