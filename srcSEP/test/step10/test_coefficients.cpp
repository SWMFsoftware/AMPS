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
  // Three spatial providers: from-dmumu, from-mfp, and the shared
  // parallel-diffusion library (added with the schema-4 [parallel_diffusion]
  // binding; Parker mover only).
  Require(C::SpatialRegistry().size() == 3,
          "COEF01 spatial registry size changed");
  Require(C::SpatialRegistry()[2].canonicalName ==
              "parallel-diffusion-library" &&
          C::SpatialRegistry()[2].units == "m2/s; m/s",
          "COEF01 parallel-diffusion-library spatial provider missing");
  // Five pitch-angle and six mean-free-path providers since the 2013-04-11
  // event test added isotropic-from-mfp and radial-rigidity-power-law (both
  // appended, so the earlier entries and their ordinals are unchanged).
  Require(C::PitchAngleRegistry().size() == 5,
          "COEF01 pitch-angle registry size changed");
  Require(C::PitchAngleRegistry()[4].canonicalName == "isotropic-from-mfp",
          "COEF01 isotropic-from-mfp pitch-angle provider missing");
  Require(C::MeanFreePathRegistry().size() == 6,
          "COEF01 mean-free-path registry size changed");
  Require(C::MeanFreePathRegistry()[5].canonicalName ==
              "radial-rigidity-power-law" &&
          C::MeanFreePathRegistry()[5].units == "m",
          "COEF01 radial-rigidity-power-law provider missing");
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

  // ---------------------------------------------------------------------
  // COEF07-COEF10 cover the providers added for the 2013-04-11 real-event
  // nightly test (srcSEP/test/event_2013_04_11).  Independent references are
  // closed-form evaluations written out here, not calls of the kernel under
  // test, so a factor, unit or exponent regression in the kernel fails.
  namespace CP = SEP::Transport::CoefficientPhysics;
  const double au = 1.495978707e11;          // IAU exact [m]
  const double c = 299792458.0;              // [m/s]
  const double e = 1.602176634e-19;          // [C]
  const double mp = 1.67262192595e-27;       // [kg], CODATA 2022
  CP::SpeciesProperties proton;
  proton.modelSpecies = 0;
  proton.name = "proton";
  proton.signedChargeC = e;
  proton.restMassKg = mp;
  proton.nucleonCount = 1.0;
  // Liu et al. (2025) Eqs. 14-16: lambda = 0.3 AU (r/AU) (pc/1 GeV)^(1/3).
  CP::RadialRigidityPowerLawParameters liu;
  liu.referenceMeanFreePathM = 0.3 * au;
  liu.referenceRadiusM = au;
  liu.referenceRigidityV = 1.0e9;
  liu.radialExponent = 1.0;
  liu.rigidityExponent = 1.0 / 3.0;

  // COEF07: closed-form oracle.  For a 10 MeV proton the exact relativistic
  // pc = sqrt(T^2 + 2 T m c^2); R = pc/e.  At r = 0.5 AU the law gives
  // 0.3 AU * 0.5 * (R / 1 GV)^(1/3).  A 1 GV proton at r_ref returns lambda_ref
  // exactly (to rounding of the log form), and the exponents are checked by
  // doubling r and multiplying R by 8 (factor 2 each).
  {
    const double restJ = mp * c * c;
    const double kineticJ = 10.0e6 * e;
    const double pcJ = std::sqrt(kineticJ * kineticJ + 2.0 * kineticJ * restJ);
    const double p = pcJ / c;
    const double rigidityV = pcJ / e;
    const double expected = 0.3 * au * 0.5 * std::cbrt(rigidityV / 1.0e9);
    const CP::MeanFreePathResult tenMeV =
        CP::EvaluateRadialRigidityPowerLawMeanFreePath(0.5 * au, proton, p, liu, c);
    Require(tenMeV.status.ok() && tenMeV.valueState == CP::ValueState::Finite &&
            Near(tenMeV.lambdaParallelM / expected, 1.0, 1.0e-13),
            "COEF07 10 MeV proton lambda differs from the closed form");
    const double p1GV = 1.0e9 * e / c;
    const CP::MeanFreePathResult reference =
        CP::EvaluateRadialRigidityPowerLawMeanFreePath(au, proton, p1GV, liu, c);
    Require(reference.status.ok() &&
            Near(reference.lambdaParallelM / (0.3 * au), 1.0, 1.0e-13),
            "COEF07 reference point does not return lambda_ref");
    const CP::MeanFreePathResult scaled =
        CP::EvaluateRadialRigidityPowerLawMeanFreePath(2.0 * au, proton,
                                                       8.0 * p1GV, liu, c);
    Require(scaled.status.ok() &&
            Near(scaled.lambdaParallelM / (4.0 * 0.3 * au), 1.0, 1.0e-13),
            "COEF07 radial or rigidity exponent is not applied");
    Require(!CP::EvaluateRadialRigidityPowerLawMeanFreePath(
                au, proton, 0.0, liu, c).status.ok(),
            "COEF07 zero rigidity was accepted");
    Require(!CP::EvaluateRadialRigidityPowerLawMeanFreePath(
                au, proton, p1GV, CP::RadialRigidityPowerLawParameters(), c)
                 .status.ok(),
            "COEF07 default (unstated) power law was accepted");
  }
  std::cout << "PASS COEF07: radial-rigidity power law matches the closed form\n";

  // COEF08: registry admission.  The provider requires its stated law, is a
  // prescribed model, and the isotropic closure cannot pair with from-spatial
  // (which would recurse) or with a coupled source.
  {
    C::Configuration law;
    law.meanFreePath = C::MeanFreePathKind::RadialRigidityPowerLaw;
    Require(!C::ValidateConfiguration(law).ok(),
            "COEF08 power-law provider accepted without parameters");
    law.radialRigidityPowerLaw = liu;
    Require(C::ValidateConfiguration(law).ok(),
            "COEF08 stated power law was rejected");
    law.pitchAngle = C::PitchAngleKind::IsotropicFromMeanFreePath;
    law.spatial = C::SpatialKind::FromMeanFreePath;
    Require(C::ValidateConfiguration(law).ok(),
            "COEF08 shared-lambda configuration was rejected");
    for (const char* mover : {"parker", "fte-dmumu", "fte-mfp"})
      Require(C::ValidateMoverCompatibility(law, mover).ok(),
              std::string("COEF08 shared-lambda configuration rejected for ") +
                  mover);
    C::Configuration cycle = law;
    cycle.meanFreePath = C::MeanFreePathKind::FromSpatial;
    Require(!C::ValidateConfiguration(cycle).ok(),
            "COEF08 isotropic-from-mfp/from-spatial cycle was accepted");
    C::Configuration coupled = law;
    coupled.source = C::SourceMode::SelfConsistent;
    Require(!C::ValidateConfiguration(coupled).ok(),
            "COEF08 prescribed power law accepted a coupled source");
    // Fingerprints: unchanged for configurations that do not select the
    // law (pre-existing published identities), sensitive to it when selected.
    C::Configuration other;
    const std::string before = C::ConfigurationFingerprint(other);
    other.radialRigidityPowerLaw = liu;
    Require(C::ConfigurationFingerprint(other) == before,
            "COEF08 unselected power-law parameters changed a fingerprint");
    C::Configuration changed = law;
    changed.radialRigidityPowerLaw.radialExponent = 2.0;
    Require(C::ConfigurationFingerprint(changed) !=
                C::ConfigurationFingerprint(law),
            "COEF08 power-law parameters are absent from the fingerprint");
  }
  std::cout << "PASS COEF08: power-law and isotropic-closure admission rules\n";

  // COEF09: closure identity.  With D_mumu = v(1-mu^2)/(2 lambda) the
  // production from-dmumu quadrature kappa = (v^2/8) Int (1-mu^2)^2/D dmu
  // must equal v lambda / 3 (the from-mfp closure used by parker), so the
  // three movers see one lambda_parallel.  Exact algebra: Int_{-1}^{1}
  // (1-mu^2) dmu = 4/3.
  {
    const double v = 4.37e7;               // ~10 MeV proton [m/s]
    const double lambda = 0.15 * au;       // [m]
    CP::SpatialQuadratureConfiguration quadrature;
    quadrature.absoluteToleranceM2PerS = 0.0;
    quadrature.relativeTolerance = 1.0e-10;
    quadrature.maximumRecursion = 40;
    const CP::SpatialDiffusionResult kappa = CP::IntegrateSpatialDiffusion(
        v,
        [&](double muValue) {
          CP::PitchAngleResult result;
          const ScalarResult d =
              C::IsotropicDmumuFromMeanFreePath(lambda, v, muValue);
          result.status = d.status;
          result.dMuMuPerS = d.value;
          result.dDmuMuDmuPerS = -v * muValue / lambda;
          result.valueState = CP::ValueState::Finite;
          return result;
        },
        quadrature);
    const ScalarResult closure = C::KappaFromMeanFreePath(lambda, v);
    Require(kappa.status.ok() && closure.status.ok() &&
            Near(kappa.kappaParallelM2PerS / (v * lambda / 3.0), 1.0, 1.0e-8) &&
            Near(closure.value / (v * lambda / 3.0), 1.0, 1.0e-15),
            "COEF09 isotropic Dmumu quadrature differs from v*lambda/3");
  }
  std::cout << "PASS COEF09: isotropic Dmumu reproduces kappa = v lambda / 3\n";

  // COEF10: the CLI reaches the registry and source controls.  The values
  // are the event-test law in SI; missing parameters and malformed source
  // controls must stop the parse before any simulation state is touched.
  {
    SEP::Util::CLI::Options parsed;
    std::string parseErrors;
    Require(Parse({"--mean-free-path-provider=radial-rigidity-power-law",
                   "--pitch-angle-diffusion-provider=isotropic-from-mfp",
                   "--spatial-diffusion-provider=from-mfp",
                   "--mfp-reference-m=4.487936121e10",
                   "--mfp-reference-radius-m=1.495978707e11",
                   "--mfp-reference-rigidity-v=1e9",
                   "--mfp-radial-exponent=1",
                   "--mfp-rigidity-exponent=0.3333333333333333",
                   "--shock-maximum-compression-ratio=none",
                   "--shock-energy-sampling=log-uniform"},
                  &parsed, &parseErrors),
            "COEF10 event-test CLI was rejected: " + parseErrors);
    Require(parsed.coefficients.meanFreePath ==
                C::MeanFreePathKind::RadialRigidityPowerLaw &&
            parsed.coefficients.pitchAngle ==
                C::PitchAngleKind::IsotropicFromMeanFreePath &&
            parsed.coefficients.radialRigidityPowerLaw.referenceMeanFreePathM ==
                4.487936121e10 &&
            parsed.coefficients.radialRigidityPowerLaw.rigidityExponent ==
                0.3333333333333333 &&
            std::isinf(parsed.shockMaximumCompressionRatio) &&
            parsed.shockMaximumCompressionRatioProvided &&
            parsed.shockLogUniformEnergySampling &&
            parsed.shockEnergySamplingProvided,
            "COEF10 parsed values did not reach the configuration");
    SEP::Util::CLI::Options missing;
    Require(!Parse({"--mean-free-path-provider=radial-rigidity-power-law",
                    "--mfp-reference-m=4.487936121e10"},
                   &missing, &parseErrors),
            "COEF10 incomplete power law was accepted");
    SEP::Util::CLI::Options defaults;
    Require(Parse({}, &defaults, &parseErrors) &&
                defaults.shockMaximumCompressionRatio == 3.0 &&
                !defaults.shockMaximumCompressionRatioProvided &&
                !defaults.shockLogUniformEnergySampling &&
                !defaults.shockEnergySamplingProvided,
            "COEF10 default source controls changed");
    for (const char* bad : {"--shock-maximum-compression-ratio=1",
                            "--shock-maximum-compression-ratio=0.5",
                            "--shock-maximum-compression-ratio=abc",
                            "--shock-energy-sampling=uniform"}) {
      SEP::Util::CLI::Options rejected;
      Require(!Parse({bad}, &rejected, &parseErrors),
              std::string("COEF10 malformed source control accepted: ") + bad);
    }
  }
  std::cout << "PASS COEF10: CLI configures the shared lambda and DSA source\n";
  return 0;
}
