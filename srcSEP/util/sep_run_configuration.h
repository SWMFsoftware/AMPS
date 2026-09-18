#ifndef SEP_UTIL_SEP_RUN_CONFIGURATION_H
#define SEP_UTIL_SEP_RUN_CONFIGURATION_H

#include "sep_common_header_path.h"
#include SRCSEP_SEP_COMMON_HEADER(sep_coefficient_registry.h)
#include SRCSEP_SEP_COMMON_HEADER(sep_injection_spectrum.h)
#include "sep_population_control.h"
#include "sep_production_mover.h"
#include "sep_sampling_products.h"
#include "sep_turbulence_core.h"

#include <cstdint>
#include <string>

namespace SEP {
namespace Run {

enum class ShockModel { Analytical, Swcme1d };
enum class CmeScenario { Fast, Slow };
enum class SwcmeFailurePolicy { Strict, ClampRadius, DiagnosticFallback };
enum class ValueSource { Default, InputFile, CommandLine };

template <class T>
struct LayeredValue {
  T value;
  ValueSource source = ValueSource::Default;
};

struct Configuration {
  LayeredValue<Mover::ProductionMover> mover;
  LayeredValue<ShockModel> shockModel;
  LayeredValue<CmeScenario> scenario;
  // D01 makes background recovery part of the immutable run contract. The
  // three fallback values use SI units and are consumed only when the explicit
  // diagnostic-fallback policy is selected.
  LayeredValue<SwcmeFailurePolicy> swcmeFailurePolicy;
  LayeredValue<double> swcmeFallbackDensityM3;
  LayeredValue<double> swcmeFallbackSpeedMPerS;
  LayeredValue<double> swcmeFallbackDivergencePerS;
  // Fingerprint produced by the canonical src/models/swcme resolver after
  // preset expansion, unit conversion, overrides, and validation.
  LayeredValue<std::string> swcmeConfigurationFingerprint;
  LayeredValue<std::uint64_t> totalIterations;
  LayeredValue<double> fieldLineSeedAreaM2;
  LayeredValue<double> shockTurbulenceEfficiency;
  LayeredValue<double> shockTurbulencePlusFraction;
  LayeredValue<int> mergeMinimum;
  LayeredValue<int> mergeMaximum;
  // The legacy pair above remains the CLI provenance surface.  This complete
  // WP33 record owns the actual population-control metric, lineage, and
  // invariant policy and is synchronized from those aliases before freeze.
  PopulationControl::Configuration populationControl;
  LayeredValue<SamplingCore::InvalidParticlePolicy> invalidSamplingPolicy;
  Transport::Coefficient::Configuration coefficients;
  Transport::NumericalTolerances numericalTolerances;
  Turbulence::Configuration turbulence;
  Injection::Configuration injection;
};

Configuration Defaults();
Transport::Status Validate(const Configuration& configuration);

// Merge applies an overlay only when the overlay field comes from a layer at
// least as authoritative as the base.  This codifies defaults < input < CLI and
// prevents a later driver literal from silently shadowing a user choice.
Configuration Merge(const Configuration& base, const Configuration& overlay);

class FrozenConfiguration {
 public:
  static Transport::Status Create(const Configuration& configuration,
                                  FrozenConfiguration* frozen);
  const Configuration& get() const { return configuration_; }
  const std::string& fingerprint() const { return fingerprint_; }
 private:
  Configuration configuration_;
  std::string fingerprint_;
};

std::string Fingerprint(const Configuration& configuration);
Transport::Status Serialize(const FrozenConfiguration& configuration,
                            std::string* text);
Transport::Status Deserialize(const std::string& text,
                              FrozenConfiguration* configuration);
Transport::Status VerifyRestartCompatibility(
    const FrozenConfiguration& requested,
    const FrozenConfiguration& checkpoint);
Transport::Status InstallActive(const Configuration& configuration);
const FrozenConfiguration& Active();
const char* ValueSourceName(ValueSource source);

}  // namespace Run
}  // namespace SEP

#endif  // SEP_UTIL_SEP_RUN_CONFIGURATION_H
