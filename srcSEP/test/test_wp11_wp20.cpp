#include "sep_coefficient_physics.h"
#include "sep_coefficient_registry.h"
#include "../util/sep_focused_transport_core.h"
#include "sep_species_source.h"
#include "sep_transport_common.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

namespace {

namespace CP = SEP::Transport::CoefficientPhysics;
namespace C = SEP::Transport::Coefficient;
using SEP::Transport::FocusedTransportBackground;
using SEP::Transport::FocusedTransportIncrement;
using SEP::Transport::FocusedTransportState;
using SEP::Transport::KeyedRandomStream;
using SEP::Transport::PitchAngleDiffusionProvider;
using SEP::Transport::PitchAngleDiffusionSample;
using SEP::Transport::Status;

int failures = 0;

void Require(bool condition, const std::string& id,
             const std::string& message) {
  if (!condition) {
    std::cerr << "FAIL " << id << ": " << message << '\n';
    ++failures;
  }
}

bool RelativeNear(double actual, double expected, double tolerance) {
  return std::fabs(actual - expected) <= tolerance *
      std::max(1.0, std::fabs(expected));
}

CP::SpeciesProperties Proton() {
  CP::SpeciesProperties species;
  species.modelSpecies = 0;
  species.name = "proton";
  species.signedChargeC = 1.602176634e-19;
  species.restMassKg = 1.67262192369e-27;
  species.nucleonCount = 1.0;
  return species;
}

CP::LocalInputView Input(double deltaB2T2 = 2.5e-19) {
  CP::LocalInputView input;
  input.source = "controlled-prescribed";
  input.representation = "analytical-branch-variance";
  input.generation = 17;
  input.checksum = 0x1234;
  input.heliocentricRadiusM = 1.495978707e11;
  input.magneticFieldT = 5.0e-9;
  input.deltaB2T2 = deltaB2T2;
  input.deltaBPlus2T2 = 0.5 * deltaB2T2;
  input.deltaBMinus2T2 = 0.5 * deltaB2T2;
  input.alfvenSpeedMPerS = 5.0e4;
  return input;
}

class IsotropicProvider : public PitchAngleDiffusionProvider {
 public:
  explicit IsotropicProvider(double rate) : rate_(rate) {}

  PitchAngleDiffusionSample Evaluate(double, double, double mu) const override {
    PitchAngleDiffusionSample sample;
    if (mu < -1.0 || mu > 1.0) {
      ++outsideCalls_;
      sample.status = Status::Error(
          SEP::Transport::StatusCode::InvalidArgument,
          "controlled provider observed out-of-domain mu");
      return sample;
    }
    sample.dMuMuPerS = rate_ * (1.0 - mu * mu);
    sample.dDmuMuDmuPerS = -2.0 * rate_ * mu;
    sample.provenance = "controlled-isotropic";
    sample.turbulenceStateIdentity = "controlled-generation-17";
    sample.valueState = CP::ValueState::Finite;
    sample.status = Status::Ok();
    return sample;
  }

  std::size_t outsideCalls() const { return outsideCalls_; }

 private:
  double rate_;
  mutable std::size_t outsideCalls_ = 0;
};

void TestWp11BoundedDiffusion() {
  // Begin with an exactly uniform pitch-angle population, which is the
  // stationary zero-flux solution for D=D0(1-mu^2).  A long ensemble detects
  // endpoint pile-up while exact-boundary particles exercise the degenerate
  // D=0/inward-drift limit.
  IsotropicProvider provider(0.8);
  FocusedTransportBackground background;
  background.pitchAngleScheme =
      SEP::Transport::BoundedPitchAngleScheme::ReflectingMilstein;
  const std::size_t particles = 40000;
  const int steps = 40;
  const double dtS = 0.0025;
  double mean = 0.0;
  double p2 = 0.0;
  std::size_t edgePopulation = 0;
  for (std::size_t i = 0; i < particles; ++i) {
    double mu = -1.0 + 2.0 * (static_cast<double>(i) + 0.5) /
        static_cast<double>(particles);
    if (i == 0) mu = -1.0;
    if (i + 1 == particles) mu = 1.0;
    FocusedTransportState state(0.0, 1.0e-20, mu);
    for (int step = 0; step < steps; ++step) {
      KeyedRandomStream random(1100, i + 1, 11, step);
      const FocusedTransportIncrement increment =
          SEP::Transport::AdvanceFocusedTransportDmumu(
              state, background, 1.67262192369e-27, 2.99792458e8,
              dtS, provider, random, NULL);
      Require(increment.status.ok(), "WP11",
              "bounded Milstein increment returned an error");
      if (!increment.status.ok()) break;
      state = increment.state;
      Require(state.mu >= -1.0 && state.mu <= 1.0, "WP11",
              "pitch angle escaped the physical interval");
    }
    mean += state.mu;
    p2 += 0.5 * (3.0 * state.mu * state.mu - 1.0);
    if (std::fabs(state.mu) > 0.99) ++edgePopulation;
  }
  mean /= particles;
  p2 /= particles;
  const double edgeFraction = static_cast<double>(edgePopulation) / particles;
  Require(std::fabs(mean) < 0.012, "WP11",
          "stationary solution lost mu-reversal symmetry");
  Require(std::fabs(p2) < 0.02, "WP11",
          "stationary isotropy developed a persistent P2 bias");
  Require(edgeFraction < 0.02, "WP11",
          "reflecting boundary produced an endpoint probability layer");
  Require(provider.outsideCalls() == 0, "WP11",
          "coefficient was evaluated outside |mu|<=1");
  std::cout << "PASS WP11: zero-flux bounded Milstein diffusion preserves isotropy and domain\n";
}

void TestWp12ErrorControls() {
  SEP::Transport::NumericalTolerances tolerances;
  Require(SEP::Transport::ValidateNumericalTolerances(tolerances).ok(),
          "WP12", "default tolerances are invalid");
  tolerances.geometryFraction = 1.1;
  Require(!SEP::Transport::ValidateNumericalTolerances(tolerances).ok(),
          "WP12", "out-of-range geometry fraction was accepted");

  const SEP::Transport::StepDoublingEstimate pass =
      SEP::Transport::EstimateStepDoublingError(1.0, 1.0005, 0.0, 1.0e-3);
  const SEP::Transport::StepDoublingEstimate fail =
      SEP::Transport::EstimateStepDoublingError(1.0, 1.01, 0.0, 1.0e-3);
  Require(pass.status.ok() && pass.accepted && fail.status.ok() && !fail.accepted,
          "WP12", "full-step/two-half-step acceptance is incorrect");

  SEP::Transport::StepDiagnostics diagnostics;
  std::vector<SEP::Transport::StepLimit> limits;
  limits.push_back(SEP::Transport::StepLimit("geometry", 0.25));
  limits.push_back(SEP::Transport::StepLimit("cooling", 0.5));
  const SEP::Transport::ScalarResult selected =
      SEP::Transport::SelectSubstep(1.0, limits, 1.0e-12, &diagnostics);
  SEP::Transport::RecordAcceptedStep(&diagnostics);
  Require(selected.status.ok() && selected.value == 0.25 &&
          diagnostics.acceptedSteps == 1 &&
          diagnostics.limiterHistogramNames.size() == 1 &&
          diagnostics.limiterHistogramNames[0] == "geometry" &&
          diagnostics.limiterHistogramCounts[0] == 1,
          "WP12", "limiter histogram or accepted-step counter is wrong");
  std::cout << "PASS WP12: named tolerances, step-doubling, and limiter counters are active\n";
}

void TestWp13SourceBinding() {
  CP::SpectrumParameters spectrum;
  spectrum.kMinAtReferencePerM = 1.0e-12;
  spectrum.kMaxAtReferencePerM = 1.0e-3;
  const CP::SpeciesProperties proton = Proton();
  CP::LocalInputView selected = Input(1.0e-18);
  CP::LocalInputView perturbed = selected;
  perturbed.deltaB2T2 *= 4.0;
  perturbed.deltaBPlus2T2 *= 4.0;
  perturbed.deltaBMinus2T2 *= 4.0;
  perturbed.checksum += 1;
  const double speed = 2.0e7;
  const double mu = 0.4;
  const CP::PitchAngleResult base = CP::EvaluateJokipiiSlab(
      selected, spectrum, proton, speed, mu);
  const CP::PitchAngleResult changed = CP::EvaluateJokipiiSlab(
      perturbed, spectrum, proton, speed, mu);
  CP::LocalInputView irrelevant = selected;
  irrelevant.alfvenSpeedMPerS *= 3.0;
  const CP::PitchAngleResult unchanged = CP::EvaluateJokipiiSlab(
      irrelevant, spectrum, proton, speed, mu);
  Require(base.status.ok() && changed.status.ok() && unchanged.status.ok() &&
          RelativeNear(changed.dMuMuPerS / base.dMuMuPerS, 4.0, 1.0e-12) &&
          RelativeNear(unchanged.dMuMuPerS, base.dMuMuPerS, 1.0e-12),
          "WP13", "selected and unselected source fields are not isolated");
  std::cout << "PASS WP13: pure coefficient input responds only to authoritative source data\n";
}

void TestWp14ConstantProvider() {
  const CP::PitchAngleResult zero = CP::EvaluateConstantDmumu(0.0, -1.0);
  const CP::PitchAngleResult positive = CP::EvaluateConstantDmumu(0.25, 0.3);
  Require(zero.status.ok() && zero.dMuMuPerS == 0.0 &&
          positive.status.ok() && positive.dMuMuPerS == 0.25 &&
          positive.dDmuMuDmuPerS == 0.0,
          "WP14", "constant provider does not return its configured SI value");
  Require(!CP::EvaluateConstantDmumu(-1.0, 0.0).status.ok() &&
          !CP::EvaluateConstantDmumu(
              std::numeric_limits<double>::quiet_NaN(), 0.0).status.ok() &&
          !CP::EvaluateConstantDmumu(
              std::numeric_limits<double>::infinity(), 0.0).status.ok(),
          "WP14", "invalid constant values were accepted");
  std::cout << "PASS WP14: constant Dmumu uses the configured value and rejects invalid input\n";
}

void TestWp15JokipiiDerivative() {
  CP::SpectrumParameters spectrum;
  spectrum.kMinAtReferencePerM = 1.0e-12;
  spectrum.kMaxAtReferencePerM = 1.0e-3;
  const CP::LocalInputView input = Input();
  const CP::SpeciesProperties proton = Proton();
  const double speed = 2.0e7;
  for (double mu : {-0.8, -0.35, 0.2, 0.7}) {
    const double h = 1.0e-6;
    const CP::PitchAngleResult center = CP::EvaluateJokipiiSlab(
        input, spectrum, proton, speed, mu);
    const CP::PitchAngleResult plus = CP::EvaluateJokipiiSlab(
        input, spectrum, proton, speed, mu + h);
    const CP::PitchAngleResult minus = CP::EvaluateJokipiiSlab(
        input, spectrum, proton, speed, mu - h);
    const double oracle = (plus.dMuMuPerS - minus.dMuMuPerS) / (2.0 * h);
    Require(center.status.ok() && plus.status.ok() && minus.status.ok() &&
            RelativeNear(center.dDmuMuDmuPerS, oracle, 2.0e-5),
            "WP15", "analytic Jokipii derivative disagrees with finite difference");
  }
  const CP::PitchAngleResult left = CP::EvaluateJokipiiSlab(
      input, spectrum, proton, speed, -1.0);
  const CP::PitchAngleResult right = CP::EvaluateJokipiiSlab(
      input, spectrum, proton, speed, 1.0);
  Require(left.status.ok() && right.status.ok() &&
          std::isfinite(left.dDmuMuDmuPerS) &&
          std::isfinite(right.dDmuMuDmuPerS),
          "WP15", "Jokipii endpoint derivative is singular");
  std::cout << "PASS WP15: Jokipii D and derivative share one finite piecewise kernel\n";
}

void TestWp16FlorinskiyProvider() {
  CP::LocalInputView plusInput = Input();
  plusInput.deltaBPlus2T2 = plusInput.deltaB2T2;
  plusInput.deltaBMinus2T2 = 0.0;
  CP::LocalInputView minusInput = plusInput;
  minusInput.deltaBPlus2T2 = 0.0;
  minusInput.deltaBMinus2T2 = minusInput.deltaB2T2;
  CP::FlorinskiyParameters parameters;
  const CP::SpeciesProperties proton = Proton();
  const double speed = 2.0e7;
  for (double mu : {-0.7, -0.2, 0.3, 0.8}) {
    const CP::PitchAngleResult plus = CP::EvaluateFlorinskiySlab(
        plusInput, parameters, proton, speed, mu);
    const CP::PitchAngleResult mirrored = CP::EvaluateFlorinskiySlab(
        minusInput, parameters, proton, speed, -mu);
    Require(plus.status.ok() && mirrored.status.ok() &&
            plus.dMuMuPerS >= 0.0 &&
            RelativeNear(plus.dMuMuPerS, mirrored.dMuMuPerS, 2.0e-12),
            "WP16", "plus/minus resonance denominators violate mirror symmetry");
  }
  const CP::PitchAngleResult endpoint = CP::EvaluateFlorinskiySlab(
      plusInput, parameters, proton, speed, 1.0);
  Require(endpoint.status.ok() && endpoint.dMuMuPerS == 0.0 &&
          std::isfinite(endpoint.dDmuMuDmuPerS),
          "WP16", "Florinskiy endpoint output is missing or non-finite");
  std::cout << "PASS WP16: Florinskiy outputs, branch denominators, and bounded derivative are repaired\n";
}

void TestWp17AdaptiveSpatialDiffusion() {
  const double speed = 3.0e6;
  const double d0 = 0.4;
  CP::SpatialQuadratureConfiguration configuration;
  configuration.absoluteToleranceM2PerS = 1.0e-3;
  configuration.relativeTolerance = 1.0e-9;
  const CP::SpatialDiffusionResult finite = CP::IntegrateSpatialDiffusion(
      speed,
      [=](double mu) {
        CP::PitchAngleResult r;
        r.status = Status::Ok();
        r.valueState = CP::ValueState::Finite;
        r.dMuMuPerS = d0 * (1.0 - mu * mu);
        r.dDmuMuDmuPerS = -2.0 * d0 * mu;
        return r;
      }, configuration);
  const double expected = speed * speed / (6.0 * d0);
  Require(finite.status.ok() && finite.valueState == CP::ValueState::Finite &&
          RelativeNear(finite.kappaParallelM2PerS, expected, 2.0e-9) &&
          finite.evaluations > 6,
          "WP17", "adaptive quadrature missed the analytical integral");

  const auto gapped = [](double mu) {
    CP::PitchAngleResult r;
    r.status = Status::Ok();
    r.valueState = CP::ValueState::Finite;
    r.dMuMuPerS = std::fabs(mu) < 0.1 ? 0.0 : 0.5 * (1.0 - mu * mu);
    return r;
  };
  configuration.gapPolicy = CP::ResonanceGapPolicy::Reject;
  const CP::SpatialDiffusionResult rejected = CP::IntegrateSpatialDiffusion(
      speed, gapped, configuration);
  configuration.gapPolicy = CP::ResonanceGapPolicy::Ballistic;
  const CP::SpatialDiffusionResult ballistic = CP::IntegrateSpatialDiffusion(
      speed, gapped, configuration);
  Require(!rejected.status.ok() &&
          rejected.status.code == SEP::Transport::StatusCode::UnresolvedCoefficient &&
          ballistic.status.ok() &&
          ballistic.valueState == CP::ValueState::Ballistic &&
          std::isinf(ballistic.kappaParallelM2PerS),
          "WP17", "resonance gap was silently divided or regularized");
  std::cout << "PASS WP17: adaptive quadrature resolves finite cases and types resonance gaps\n";
}

void TestWp18SpeciesAwareness() {
  CP::SpeciesProperties proton = Proton();
  CP::SpeciesProperties alpha = proton;
  alpha.modelSpecies = 1;
  alpha.name = "alpha";
  alpha.signedChargeC = 2.0 * proton.signedChargeC;
  alpha.restMassKg = 4.0 * proton.restMassKg;
  alpha.nucleonCount = 4.0;
  CP::SpeciesProperties electron = proton;
  electron.modelSpecies = 2;
  electron.name = "electron";
  electron.signedChargeC = -proton.signedChargeC;
  electron.restMassKg = 9.1093837015e-31;
  electron.nucleonCount = 0.0;
  const double b = 5.0e-9;
  const double p = 1.0e-20;
  const SEP::Transport::ScalarResult omegaP = CP::GyrofrequencyRadPerS(b, proton);
  const SEP::Transport::ScalarResult omegaA = CP::GyrofrequencyRadPerS(b, alpha);
  const SEP::Transport::ScalarResult rP = CP::LarmorRadiusM(p, b, proton);
  const SEP::Transport::ScalarResult rA = CP::LarmorRadiusM(p, b, alpha);
  const SEP::Transport::ScalarResult rigidityE =
      CP::RigidityVolt(p, electron, 2.99792458e8);
  Require(omegaP.status.ok() && omegaA.status.ok() &&
          RelativeNear(omegaA.value / omegaP.value, 0.5, 1.0e-14) &&
          RelativeNear(rA.value / rP.value, 0.5, 1.0e-14) &&
          rigidityE.status.ok() && CP::ValidateSpecies(electron).ok(),
          "WP18", "charge/mass species scaling is incorrect");

  // The same species metadata also controls source normalization and the
  // total-energy/per-nucleon boundary.  This catches the particularly subtle
  // error where 10 MeV/nucleon alpha particles are injected as 10 MeV total.
  namespace SS = SEP::Transport::SpeciesSource;
  SS::Configuration protonSource;
  protonSource.species = proton;
  protonSource.abundanceFraction = 9.0;
  protonSource.injectionEfficiency = 3.4e-4;
  SS::Configuration alphaSource;
  alphaSource.species = alpha;
  alphaSource.abundanceFraction = 1.0;
  alphaSource.injectionEfficiency = 1.5e-4;
  alphaSource.energyConvention = SS::EnergyConvention::PerNucleon;
  std::vector<SS::Configuration> normalized;
  Require(SS::ValidateAndNormalize({protonSource, alphaSource},
                                   &normalized).ok() &&
          normalized.size() == 2 &&
          normalized[0].abundanceFraction +
              normalized[1].abundanceFraction == 1.0,
          "WP18", "multi-species source abundance did not close exactly");
  const SEP::Transport::ScalarResult alphaEnergy =
      SS::TotalKineticEnergyJ(alphaSource, 2.0e-12);
  Require(alphaEnergy.status.ok() && alphaEnergy.value == 8.0e-12,
          "WP18", "energy-per-nucleon was not converted to total energy");
  electron.nucleonCount = 0.0;
  SS::Configuration invalidElectronSource;
  invalidElectronSource.species = electron;
  invalidElectronSource.energyConvention = SS::EnergyConvention::PerNucleon;
  Require(!SS::ValidateAndNormalize({invalidElectronSource},
                                    &normalized).ok(),
          "WP18", "electron energy-per-nucleon configuration was accepted");
  Require(SS::SetActiveConfiguration({protonSource, alphaSource}).ok(),
          "WP18", "valid species source table was not installable");
  SS::Configuration missing;
  Require(!SS::FindActiveConfiguration(77, &missing).ok(),
          "WP18", "missing species silently inherited another source entry");
  std::cout << "PASS WP18: species scaling and normalized source conventions are explicit\n";
}

void TestWp19NamedScales() {
  C::Configuration configuration;
  Require(C::ValidateConfiguration(configuration).ok(), "WP19",
          "default named scales are invalid");
  configuration.prescribedDeltaBOverB = -0.1;
  Require(!C::ValidateConfiguration(configuration).ok(), "WP19",
          "negative prescribed turbulence ratio was accepted");
  configuration = C::Configuration();
  configuration.spectrum.kMaxAtReferencePerM =
      configuration.spectrum.kMinAtReferencePerM;
  Require(!C::ValidateConfiguration(configuration).ok(), "WP19",
          "unordered spectrum bounds were accepted");

  // Reachability is stronger than parse coverage: perturb named scales and
  // require both the physical kernel and stable run identity to change.
  C::Configuration base;
  base.pitchAngle = C::PitchAngleKind::Jokipii1966;
  C::Configuration spectrumChanged = base;
  spectrumChanged.spectrum.kMinAtReferencePerM *= 2.0;
  C::Configuration correlationChanged = base;
  correlationChanged.correlationLengthAt1AuM *= 2.0;
  C::Configuration policyChanged = base;
  policyChanged.amplitudePolicy = C::TurbulenceAmplitudePolicy::LimitToMeanField;
  C::Configuration constantChanged = base;
  constantChanged.constantDmumuPerS = 0.25;
  Require(C::ConfigurationFingerprint(base) !=
              C::ConfigurationFingerprint(spectrumChanged) &&
          C::ConfigurationFingerprint(base) !=
              C::ConfigurationFingerprint(correlationChanged) &&
          C::ConfigurationFingerprint(base) !=
              C::ConfigurationFingerprint(policyChanged) &&
          C::ConfigurationFingerprint(base) !=
              C::ConfigurationFingerprint(constantChanged),
          "WP19", "named coefficient scale or policy is absent from fingerprint");

  const CP::PitchAngleResult dBase = CP::EvaluateJokipiiSlab(
      Input(), base.spectrum, Proton(), 2.0e7, 0.4);
  const CP::PitchAngleResult dChanged = CP::EvaluateJokipiiSlab(
      Input(), spectrumChanged.spectrum, Proton(), 2.0e7, 0.4);
  const CP::MeanFreePathResult lambdaBase =
      CP::EvaluateCorrelationMeanFreePath(
          Input(), Proton(), 1.0e-20, base.correlationLengthAt1AuM,
          base.spectrum.referenceRadiusM);
  const CP::MeanFreePathResult lambdaChanged =
      CP::EvaluateCorrelationMeanFreePath(
          Input(), Proton(), 1.0e-20,
          correlationChanged.correlationLengthAt1AuM,
          correlationChanged.spectrum.referenceRadiusM);
  Require(dBase.status.ok() && dChanged.status.ok() &&
          dBase.dMuMuPerS != dChanged.dMuMuPerS &&
          lambdaBase.status.ok() && lambdaChanged.status.ok() &&
          lambdaBase.lambdaParallelM != lambdaChanged.lambdaParallelM,
          "WP19", "named spectrum or correlation scale is not physically reachable");
  std::cout << "PASS WP19: turbulence ratio, correlation length, and spectrum scales are validated\n";
}

void TestWp20BallisticCompatibility() {
  C::Configuration configuration;
  configuration.spatial = C::SpatialKind::FromMeanFreePath;
  configuration.invalidPolicy = C::InvalidPolicy::Ballistic;
  Require(!C::ValidateMoverCompatibility(configuration, "parker").ok(),
          "WP20", "Parker accepted a ballistic spatial-from-MFP combination");
  Require(C::ValidateMoverCompatibility(configuration, "fte-mfp").ok(),
          "WP20", "event-driven MFP rejected its defined ballistic state");

  CP::LocalInputView noTurbulence = Input(0.0);
  const CP::MeanFreePathResult ballistic =
      CP::EvaluateCorrelationMeanFreePath(
          noTurbulence, Proton(), 1.0e-20, 1.0e9, 1.495978707e11);
  Require(ballistic.status.ok() &&
          ballistic.valueState == CP::ValueState::Ballistic &&
          std::isinf(ballistic.lambdaParallelM),
          "WP20", "zero turbulence did not produce typed ballistic lambda");
  Require(!C::KappaFromMeanFreePath(
              std::numeric_limits<double>::infinity(), 1.0e7).status.ok(),
          "WP20", "infinite kappa entered an unsupported diffusion operator");

  // Exhaust the cross-field policy matrix at the registry boundary.  Each
  // configuration that survives basic provider validation must also have a
  // defined result for every canonical mover; assertions below restate the
  // forbidden infinity paths independently of the validator implementation.
  const C::SourceMode sources[] = {
      C::SourceMode::Prescribed, C::SourceMode::SelfConsistent,
      C::SourceMode::Swmf};
  const C::SpatialKind spatialKinds[] = {
      C::SpatialKind::FromPitchAngle, C::SpatialKind::FromMeanFreePath};
  const C::PitchAngleKind pitchKinds[] = {
      C::PitchAngleKind::Configured, C::PitchAngleKind::Constant,
      C::PitchAngleKind::Jokipii1966, C::PitchAngleKind::Florinskiy};
  const C::MeanFreePathKind mfpKinds[] = {
      C::MeanFreePathKind::Qlt, C::MeanFreePathKind::Qlt1,
      C::MeanFreePathKind::Tenishev2005, C::MeanFreePathKind::Chen2024,
      C::MeanFreePathKind::FromSpatial};
  const C::InvalidPolicy invalidPolicies[] = {
      C::InvalidPolicy::Fail, C::InvalidPolicy::Ballistic};
  const CP::ResonanceGapPolicy gapPolicies[] = {
      CP::ResonanceGapPolicy::Reject, CP::ResonanceGapPolicy::Ballistic};
  const char* movers[] = {"parker", "fte-dmumu", "fte-mfp"};
  std::size_t accepted = 0, rejected = 0;
  for (C::SourceMode source : sources)
    for (C::SpatialKind spatial : spatialKinds)
      for (C::PitchAngleKind pitch : pitchKinds)
        for (C::MeanFreePathKind mfp : mfpKinds)
          for (C::InvalidPolicy invalid : invalidPolicies)
            for (CP::ResonanceGapPolicy gap : gapPolicies)
              for (const char* mover : movers) {
                C::Configuration candidate;
                candidate.source = source;
                candidate.spatial = spatial;
                candidate.pitchAngle = pitch;
                candidate.meanFreePath = mfp;
                candidate.invalidPolicy = invalid;
                candidate.resonanceGapPolicy = gap;
                const Status status =
                    C::ValidateMoverCompatibility(candidate, mover);
                if (status.ok()) {
                  ++accepted;
                  Require(!(std::string(mover) == "parker" &&
                            spatial == C::SpatialKind::FromMeanFreePath &&
                            invalid == C::InvalidPolicy::Ballistic),
                          "WP20", "matrix admitted ballistic Parker MFP");
                  Require(!(std::string(mover) == "parker" &&
                            gap == CP::ResonanceGapPolicy::Ballistic),
                          "WP20", "matrix admitted ballistic Parker gap");
                }
                else ++rejected;
              }
  Require(accepted > 0 && rejected > 0, "WP20",
          "compatibility matrix did not exercise both supported and rejected states");

  CP::LocalInputView weakTurbulence = Input(1.0e-40);
  const CP::MeanFreePathResult largeFinite =
      CP::EvaluateCorrelationMeanFreePath(
          weakTurbulence, Proton(), 1.0e-20, 1.0e9, 1.495978707e11);
  Require(largeFinite.status.ok() &&
          largeFinite.valueState == CP::ValueState::Finite &&
          std::isfinite(largeFinite.lambdaParallelM) &&
          largeFinite.lambdaParallelM > 1.0e20,
          "WP20", "large finite MFP was confused with typed ballistic state");
  std::cout << "PASS WP20: ballistic lambda is typed and rejected before unsupported Parker use\n";
}

}  // namespace

int main() {
  TestWp11BoundedDiffusion();
  TestWp12ErrorControls();
  TestWp13SourceBinding();
  TestWp14ConstantProvider();
  TestWp15JokipiiDerivative();
  TestWp16FlorinskiyProvider();
  TestWp17AdaptiveSpatialDiffusion();
  TestWp18SpeciesAwareness();
  TestWp19NamedScales();
  TestWp20BallisticCompatibility();
  if (failures != 0) {
    std::cerr << "WP11-WP20 focused tests: " << failures << " failure(s)\n";
    return EXIT_FAILURE;
  }
  std::cout << "WP11-WP20 focused tests: PASS\n";
  return EXIT_SUCCESS;
}
