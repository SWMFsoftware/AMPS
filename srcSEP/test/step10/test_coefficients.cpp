#include "sep_coefficient_registry.h"
#include "sep_cli.h"

#include <cmath>
#include <iostream>
#include <limits>
#include <string>
#include <sstream>
#include <vector>

namespace {

void Require(bool condition, const std::string& message) {
  if (!condition) {
    std::cerr << "FAIL: " << message << '\n';
    std::exit(1);
  }
}

bool Near(double actual, double expected, double relativeTolerance) {
  const double scale = std::max(1.0, std::fabs(expected));
  return std::fabs(actual - expected) <= relativeTolerance * scale;
}

bool Parse(const std::vector<std::string>& arguments,
           SEP::Util::CLI::Options* options, std::string* errors) {
  std::vector<std::string> owned;
  owned.push_back("srcsep-test");
  owned.insert(owned.end(), arguments.begin(), arguments.end());
  std::vector<char*> argv;
  for (std::string& value : owned) argv.push_back(&value[0]);
  std::ostringstream output;
  std::ostringstream error;
  const bool ok = SEP::Util::CLI::ParseCommandLine(
      static_cast<int>(argv.size()), argv.data(), *options, output, error);
  if (errors) *errors = error.str();
  return ok;
}

}  // namespace

int main() {
  namespace B = SEP::Background;
  namespace C = SEP::Transport::Coefficient;
  using SEP::Transport::ScalarResult;

  // COEF01: the registry itself is the unit/schema contract consumed by help,
  // documentation, and providers.  These assertions prevent a renamed model
  // or dimension from silently changing a published configuration.
  Require(C::SpatialRegistry().size() == 2,
          "COEF01 spatial registry size changed");
  Require(C::PitchAngleRegistry().size() == 4,
          "COEF01 pitch-angle registry size changed");
  Require(C::MeanFreePathRegistry().size() == 5,
          "COEF01 mean-free-path registry size changed");
  Require(C::SpatialRegistry()[0].units == "m2/s; m/s",
          "COEF01 spatial SI units missing");
  Require(C::PitchAngleRegistry()[0].units == "1/s; 1/s",
          "COEF01 pitch-angle SI units missing");
  Require(C::MeanFreePathRegistry()[0].units == "m",
          "COEF01 mean-free-path SI units missing");
  std::cout << "PASS COEF01: canonical registries publish SI units and schemas\n";

  // COEF02: reject the only recursive conversion and prevent externally owned
  // turbulence from being relabelled as a direct analytical lambda source.
  C::Configuration configuration;
  Require(C::ValidateConfiguration(configuration).ok(),
          "COEF02 default configuration must remain valid");
  configuration.spatial = C::SpatialKind::FromMeanFreePath;
  configuration.meanFreePath = C::MeanFreePathKind::FromSpatial;
  Require(!C::ValidateConfiguration(configuration).ok(),
          "COEF02 conversion cycle was accepted");
  configuration = C::Configuration();
  configuration.source = C::SourceMode::SelfConsistent;
  configuration.meanFreePath = C::MeanFreePathKind::Qlt;
  configuration.pitchAngle = C::PitchAngleKind::Configured;
  Require(!C::ValidateConfiguration(configuration).ok(),
          "COEF02 coupled source accepted a legacy callback");
  configuration.pitchAngle = C::PitchAngleKind::Jokipii1966;
  Require(C::ValidateConfiguration(configuration).ok(),
          "COEF02 source-bound QLT configuration was rejected");
  SEP::Util::CLI::Options cliOptions;
  std::string cliErrors;
  Require(Parse({"--coefficient-source=swmf",
                 "--spatial-diffusion-provider", "from-mfp",
                 "--pitch-angle-diffusion-provider=configured",
                 "--mean-free-path-provider", "qlt1",
                 "--invalid-coefficient-policy=ballistic"},
                &cliOptions, &cliErrors) == false,
          "COEF02 CLI accepted an incompatible SWMF analytical MFP");
  cliOptions = SEP::Util::CLI::Options();
  // WP35 makes source compatibility a complete cross-product contract.  A
  // valid SWMF coefficient selection therefore names read-only SWMF turbulence
  // and disables particle-wave mutation instead of inheriting the default
  // locally evolved source and coupling policy.
  Require(Parse({"--coefficient-source=swmf",
                 "--turbulence-source=swmf-read-only",
                 "--turbulence-coupling-policy=disabled",
                 "--spatial-diffusion-provider", "from-dmumu",
                 "--pitch-angle-diffusion-provider=jokipii-1966",
                 "--mean-free-path-provider=from-spatial",
                 "--invalid-coefficient-policy", "ballistic"},
                &cliOptions, &cliErrors),
          "COEF02 CLI rejected a valid explicit SWMF configuration: " +
              cliErrors);
  Require(cliOptions.coefficients.source == C::SourceMode::Swmf &&
          cliOptions.coefficients.meanFreePath ==
              C::MeanFreePathKind::FromSpatial &&
          cliOptions.coefficients.invalidPolicy ==
              C::InvalidPolicy::Ballistic,
          "COEF02 CLI did not retain canonical coefficient choices");
  std::cout << "PASS COEF02: incompatible parameter combinations are rejected\n";

  // COEF03: all lambda/kappa and isotropic Dmumu conversions meet at these
  // functions.  Round trips at representative SI values detect factor-of-three
  // and factor-of-two regressions directly.
  const double speedMPerS = 2.0e7;
  const double lambdaM = 7.5e9;
  const double mu = 0.35;
  const ScalarResult kappa = C::KappaFromMeanFreePath(lambdaM, speedMPerS);
  Require(kappa.status.ok(), "COEF03 lambda-to-kappa conversion failed");
  const ScalarResult lambdaRoundTrip =
      C::MeanFreePathFromKappa(kappa.value, speedMPerS);
  Require(lambdaRoundTrip.status.ok() &&
          Near(lambdaRoundTrip.value, lambdaM, 1.0e-14),
          "COEF03 lambda/kappa round trip failed");
  const ScalarResult dmumu =
      C::IsotropicDmumuFromMeanFreePath(lambdaM, speedMPerS, mu);
  const ScalarResult lambdaFromDmumu =
      C::MeanFreePathFromIsotropicDmumu(dmumu.value, speedMPerS, mu);
  Require(dmumu.status.ok() && lambdaFromDmumu.status.ok() &&
          Near(lambdaFromDmumu.value, lambdaM, 1.0e-14),
          "COEF03 isotropic Dmumu round trip failed");
  std::cout << "PASS COEF03: centralized coefficient conversions round-trip\n";

  // COEF04: source authority changes provenance/ownership validation, never
  // the mathematical SI closure.  Prescribed analytic and imported SWMF states
  // therefore produce the same conversion while retaining distinct labels.
  Require(C::ValidateSourceAgainstBackground(
              C::SourceMode::Prescribed, B::Provider::Analytic,
              B::Ownership::ModelOwned).ok(),
          "COEF04 prescribed analytic source rejected");
  Require(C::ValidateSourceAgainstBackground(
              C::SourceMode::Swmf, B::Provider::Swmf,
              B::Ownership::ImportedReadOnly).ok(),
          "COEF04 imported SWMF source rejected");
  Require(std::string(C::SourceName(C::SourceMode::Prescribed)) !=
          C::SourceName(C::SourceMode::Swmf),
          "COEF04 source identities collapsed");
  const ScalarResult importedKappa =
      C::KappaFromMeanFreePath(lambdaM, speedMPerS);
  Require(importedKappa.status.ok() &&
          Near(importedKappa.value, kappa.value, 1.0e-15),
          "COEF04 source-independent SI conversion changed");
  std::cout << "PASS COEF04: analytic/imported sources preserve identity and physics\n";

  // COEF05: invalid domains return structured status instead of NaN, clamping,
  // or process termination.  Production policy is applied only after this
  // explicit error reaches the PIC-facing provider.
  Require(!C::KappaFromMeanFreePath(-1.0, speedMPerS).status.ok(),
          "COEF05 negative lambda was accepted");
  Require(!C::MeanFreePathFromKappa(
              std::numeric_limits<double>::quiet_NaN(), speedMPerS).status.ok(),
          "COEF05 NaN kappa was accepted");
  Require(!C::IsotropicDmumuFromMeanFreePath(
              lambdaM, speedMPerS, 1.01).status.ok(),
          "COEF05 |mu|>1 was accepted");
  Require(!C::ValidateSourceAgainstBackground(
              C::SourceMode::Swmf, B::Provider::Analytic,
              B::Ownership::ModelOwned).ok(),
          "COEF05 mislabeled SWMF state was accepted");
  std::cout << "PASS COEF05: invalid coefficients and source states return errors\n";

  // WP14 CLI reachability: the parser must retain the exact SI value and mark
  // it as an explicit override.  Feeding that retained value into the same
  // pure kernel used by the PIC adapter proves the configuration is not merely
  // printed while a different constant controls scattering.
  cliOptions = SEP::Util::CLI::Options();
  Require(Parse({"--pitch-angle-diffusion-provider=constant",
                 "--constant-dmumu=0.375"},
                &cliOptions, &cliErrors),
          "COEF06 constant Dmumu CLI configuration failed: " + cliErrors);
  const SEP::Transport::CoefficientPhysics::PitchAngleResult constant =
      SEP::Transport::CoefficientPhysics::EvaluateConstantDmumu(
          cliOptions.constantDmumuPerS, 0.2);
  Require(cliOptions.constantDmumuProvided && constant.status.ok() &&
          constant.dMuMuPerS == 0.375 &&
          constant.dDmuMuDmuPerS == 0.0,
          "COEF06 parsed constant did not reach the production kernel value");
  std::cout << "PASS COEF06: CLI constant Dmumu reaches the pure provider value\n";
  return 0;
}
