#include "coefficient_providers.h"

#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstring>
#include <limits>
#include <sstream>

namespace SEP {
namespace Transport {
namespace PICAdapter {
namespace {

namespace CP = CoefficientPhysics;

std::atomic<std::uint64_t> gInvalidSamples(0);
std::atomic<std::uint64_t> gBallisticSubstitutions(0);
std::atomic<std::uint64_t> gAmplitudeRegularizations(0);
std::atomic<std::uint64_t> gQuadratureEvaluations(0);

Status ValidateConfiguredSource() {
  const Background::BackgroundSnapshot& snapshot =
      Background::SnapshotStore::Instance().AcquireForMover();
  return Coefficient::ValidateSourceAgainstBackground(
      Coefficient::ActiveConfiguration().source, snapshot.provider(),
      snapshot.ownership());
}

ParticleContext AtRelativeArcLength(const ParticleContext& origin,
                                    double relativeArcLengthM,
                                    Status* status) {
  ParticleContext sampled = origin;
  if (!status || !std::isfinite(relativeArcLengthM)) {
    if (status) *status = Status::Error(
        StatusCode::InvalidArgument,
        "coefficient sample location must be finite");
    return sampled;
  }
  sampled.state.coordinate =
      PIC::FieldLine::FieldLinesAll[origin.state.fieldLineId].move(
          origin.state.coordinate, relativeArcLengthM, sampled.segment);
  sampled.segment =
      PIC::FieldLine::FieldLinesAll[origin.state.fieldLineId].GetSegment(
          sampled.state.coordinate);
  *status = sampled.segment
      ? Status::Ok()
      : Status::Error(StatusCode::OutOfDomain,
                      "coefficient sample lies outside the field line");
  return sampled;
}

std::uint64_t HashBytes(std::uint64_t hash, const void* data,
                        std::size_t size) {
  const unsigned char* bytes = static_cast<const unsigned char*>(data);
  for (std::size_t i = 0; i < size; ++i) {
    hash ^= bytes[i];
    hash *= UINT64_C(1099511628211);
  }
  return hash;
}

CP::SpeciesProperties SpeciesFor(const ParticleContext& context) {
  CP::SpeciesProperties species;
  species.modelSpecies = context.state.species;
  std::ostringstream name;
  name << "PIC-species-" << context.state.species;
  species.name = name.str();
  species.signedChargeC =
      PIC::MolecularData::GetElectricCharge(context.state.species);
  species.restMassKg = PIC::MolecularData::GetMass(context.state.species);

  // PIC's public molecular table does not expose an atomic-mass number in all
  // enclosing AMPS builds.  For ions, the nearest proton-mass multiple is the
  // only available deterministic metadata; electrons deliberately carry A=0.
  // A future host that exposes isotope data should replace this adapter value,
  // not reintroduce proton constants into the coefficient kernels.
  const double protonMassKg = 1.67262192369e-27;
  const double ratio = species.restMassKg / protonMassKg;
  species.nucleonCount = ratio >= 0.5 ? std::floor(ratio + 0.5) : 0.0;
  return species;
}

Status BuildSourceBoundInput(const ParticleContext& context,
                             CP::LocalInputView* input) {
  if (!input || !context.segment) {
    return Status::Error(StatusCode::InvalidArgument,
                         "coefficient source adapter requires a segment");
  }
  LocalBackgroundView background;
  Status status = EvaluateLocalBackground(context, &background);
  if (!status.ok()) return status;

  double positionM[3] = {0.0, 0.0, 0.0};
  // A particle coordinate stores segment index plus local fraction, whereas
  // cFieldLineSegment::GetCartesian expects only the local [0,1] fraction.
  // Passing the global coordinate would extrapolate positions on every segment
  // after segment zero and corrupt all radius-scaled coefficient parameters.
  const double localFraction = context.state.coordinate -
      std::floor(context.state.coordinate);
  context.segment->GetCartesian(positionM, localFraction);
  input->heliocentricRadiusM = Vector3D::Length(positionM);
  input->magneticFieldT = background.magneticFieldMagnitudeT;
  input->alfvenSpeedMPerS = background.alfvenSpeedMPerS;
  input->generation = background.generation;
  input->source = Coefficient::SourceName(
      Coefficient::ActiveConfiguration().source);
  input->representation = "integrated-branch-energy";

  const Coefficient::Configuration& configuration =
      Coefficient::ActiveConfiguration();
  if (configuration.source == Coefficient::SourceMode::Prescribed) {
    const double ratio = configuration.prescribedDeltaBOverB;
    input->deltaB2T2 = ratio * ratio * input->magneticFieldT *
        input->magneticFieldT;
    input->deltaBPlus2T2 = 0.5 * input->deltaB2T2;
    input->deltaBMinus2T2 = 0.5 * input->deltaB2T2;
    input->representation = "prescribed-deltaB-over-B";
  }
  else if (configuration.source == Coefficient::SourceMode::SelfConsistent) {
    double* wave = context.segment->GetDatum_ptr(
        SEP::AlfvenTurbulence_Kolmogorov::WaveEnergyDensity);
    if (!wave || !std::isfinite(wave[0]) || !std::isfinite(wave[1]) ||
        wave[0] < 0.0 || wave[1] < 0.0) {
      return Status::Error(StatusCode::InvalidCoefficient,
          "self-consistent coefficient source lacks authoritative E+/E- density");
    }
    // For the model's Alfvén-wave energy convention, magnetic variance is
    // mu0 times energy density.  Both branches are converted before entering
    // the pure SI kernel; the kernel never reads PIC datums or global state.
    input->deltaBPlus2T2 = VacuumPermeability * wave[0];
    input->deltaBMinus2T2 = VacuumPermeability * wave[1];
    input->deltaB2T2 = input->deltaBPlus2T2 + input->deltaBMinus2T2;
  }
  else {
    PIC::FieldLine::cFieldLineVertex* begin = context.segment->GetBegin();
    PIC::FieldLine::cFieldLineVertex* end = context.segment->GetEnd();
    double* wave0 = begin ? begin->GetDatum_ptr(
        PIC::FieldLine::DatumAtVertexPlasmaWaves) : NULL;
    double* wave1 = end ? end->GetDatum_ptr(
        PIC::FieldLine::DatumAtVertexPlasmaWaves) : NULL;
    if (!wave0 || !wave1) {
      return Status::Error(StatusCode::InvalidCoefficient,
          "SWMF coefficient source lacks imported branch wave energy");
    }
    const double fraction = context.state.coordinate -
        std::floor(context.state.coordinate);
    const double plus = (1.0 - fraction) * wave0[0] + fraction * wave1[0];
    const double minus = (1.0 - fraction) * wave0[1] + fraction * wave1[1];
    if (!std::isfinite(plus) || !std::isfinite(minus) ||
        plus < 0.0 || minus < 0.0) {
      return Status::Error(StatusCode::InvalidCoefficient,
                           "SWMF branch wave energy is invalid");
    }
    input->deltaBPlus2T2 = VacuumPermeability * plus;
    input->deltaBMinus2T2 = VacuumPermeability * minus;
    input->deltaB2T2 = input->deltaBPlus2T2 + input->deltaBMinus2T2;
    input->representation = "swmf-read-only-branch-energy";
  }

  const double meanField2 = input->magneticFieldT * input->magneticFieldT;
  if (input->deltaB2T2 > meanField2) {
    if (configuration.amplitudePolicy ==
        Coefficient::TurbulenceAmplitudePolicy::Reject) {
      return Status::Error(StatusCode::InvalidCoefficient,
          "deltaB/B exceeds the provider validity domain; select explicit "
          "limit-to-mean-field regularization to continue");
    }
    // Preserve branch imbalance while limiting only the total variance.  This
    // correction is counted so a run cannot silently hide an out-of-domain
    // turbulence amplitude.
    const double scale = meanField2 / input->deltaB2T2;
    input->deltaBPlus2T2 *= scale;
    input->deltaBMinus2T2 *= scale;
    input->deltaB2T2 = meanField2;
    gAmplitudeRegularizations.fetch_add(1, std::memory_order_relaxed);
  }

  std::uint64_t checksum = UINT64_C(1469598103934665603);
  checksum = HashBytes(checksum, input->source.data(), input->source.size());
  checksum = HashBytes(checksum, input->representation.data(),
                       input->representation.size());
  checksum = HashBytes(checksum, &input->generation, sizeof(input->generation));
  checksum = HashBytes(checksum, &input->heliocentricRadiusM,
                       sizeof(input->heliocentricRadiusM));
  checksum = HashBytes(checksum, &input->magneticFieldT,
                       sizeof(input->magneticFieldT));
  checksum = HashBytes(checksum, &input->deltaBPlus2T2,
                       sizeof(input->deltaBPlus2T2));
  checksum = HashBytes(checksum, &input->deltaBMinus2T2,
                       sizeof(input->deltaBMinus2T2));
  input->checksum = checksum == 0 ? UINT64_C(1) : checksum;
  return CP::ValidateLocalInput(*input);
}

std::string Provenance(const char* quantity, const char* provider,
                       const CP::LocalInputView& input) {
  const Background::BackgroundSnapshot& snapshot =
      Background::SnapshotStore::Instance().AcquireForMover();
  std::ostringstream out;
  out << quantity << '=' << provider
      << ";source=" << input.source
      << ";representation=" << input.representation
      << ";background=" << Background::ProviderName(snapshot.provider())
      << ";fingerprint=" << snapshot.configuration_fingerprint()
      << ";coefficient-fingerprint="
      << Coefficient::ConfigurationFingerprint(
             Coefficient::ActiveConfiguration())
      << ";generation=" << input.generation
      << ";source-checksum=" << input.checksum;
  return out.str();
}

void CopyIdentity(const CP::LocalInputView& input,
                  PitchAngleDiffusionSample* sample) {
  sample->source = input.source;
  sample->representation = input.representation;
  sample->generation = input.generation;
  sample->sourceChecksum = input.checksum;
  const CP::SpectrumParameters& spectrum =
      Coefficient::ActiveConfiguration().spectrum;
  sample->spectrumKMinPerM = spectrum.kMinAtReferencePerM * std::pow(
      spectrum.referenceRadiusM / input.heliocentricRadiusM,
      spectrum.kMinRadialExponent);
  sample->spectrumKMaxPerM = spectrum.kMaxAtReferencePerM * std::pow(
      spectrum.referenceRadiusM / input.heliocentricRadiusM,
      spectrum.kMaxRadialExponent);
}

void CopyIdentity(const CP::LocalInputView& input,
                  MeanFreePathSample* sample) {
  sample->source = input.source;
  sample->representation = input.representation;
  sample->generation = input.generation;
  sample->sourceChecksum = input.checksum;
}

CP::PitchAngleResult EvaluatePitchKernel(
    const ParticleContext& context, const CP::LocalInputView& input,
    const CP::SpeciesProperties& species, double momentumKgMPerS, double mu) {
  const ScalarResult speed = SpeedFromMomentum(
      momentumKgMPerS, species.restMassKg, SpeedOfLight);
  if (!speed.status.ok()) {
    CP::PitchAngleResult result;
    result.status = speed.status;
    return result;
  }
  const Coefficient::Configuration& configuration =
      Coefficient::ActiveConfiguration();
  switch (configuration.pitchAngle) {
    case Coefficient::PitchAngleKind::Constant:
      return CP::EvaluateConstantDmumu(
          configuration.constantDmumuPerS, mu);
    case Coefficient::PitchAngleKind::Jokipii1966:
      return CP::EvaluateJokipiiSlab(
          input, configuration.spectrum, species, speed.value, mu);
    case Coefficient::PitchAngleKind::Florinskiy:
      return CP::EvaluateFlorinskiySlab(
          input, configuration.florinskiy, species, speed.value, mu);
    case Coefficient::PitchAngleKind::Configured:
      break;
  }

  CP::PitchAngleResult result;
  result.valueState = CP::ValueState::Finite;
  if (SEP::Diffusion::GetPitchAngleDiffusionCoefficient == NULL) {
    result.status = Status::Ok();
    return result;
  }
  const double boundedMu = std::max(-1.0, std::min(1.0, mu));
  const double vParallel = speed.value * boundedMu;
  const double vNormal = speed.value *
      std::sqrt(std::max(0.0, 1.0 - boundedMu * boundedMu));
  SEP::Diffusion::GetPitchAngleDiffusionCoefficient(
      result.dMuMuPerS, result.dDmuMuDmuPerS, boundedMu,
      vParallel, vNormal, context.state.species,
      context.state.coordinate, context.segment);

  if (SEP::Diffusion::PitchAngleDifferentialMode ==
      SEP::Diffusion::PitchAngleDifferentialModeNumerical) {
    const double h = SEP::Diffusion::muNumericalDifferentiationStep;
    if (!std::isfinite(h) || h <= 0.0) {
      result.status = Status::Error(StatusCode::InvalidCoefficient,
                                    "Dmumu derivative step must be positive");
      return result;
    }
    const double muMinus = std::max(-1.0, boundedMu - h);
    const double muPlus = std::min(1.0, boundedMu + h);
    double dMinus = 0.0, unusedMinus = 0.0;
    double dPlus = 0.0, unusedPlus = 0.0;
    SEP::Diffusion::GetPitchAngleDiffusionCoefficient(
        dMinus, unusedMinus, muMinus, speed.value * muMinus,
        speed.value * std::sqrt(std::max(0.0, 1.0 - muMinus * muMinus)),
        context.state.species, context.state.coordinate, context.segment);
    SEP::Diffusion::GetPitchAngleDiffusionCoefficient(
        dPlus, unusedPlus, muPlus, speed.value * muPlus,
        speed.value * std::sqrt(std::max(0.0, 1.0 - muPlus * muPlus)),
        context.state.species, context.state.coordinate, context.segment);
    if (!(muPlus > muMinus)) {
      result.status = Status::Error(StatusCode::InvalidCoefficient,
                                    "Dmumu derivative stencil collapsed");
      return result;
    }
    result.dDmuMuDmuPerS = (dPlus - dMinus) / (muPlus - muMinus);
    if (boundedMu == -1.0 || boundedMu == 1.0)
      result.derivativeState = CP::DerivativeState::OneSidedLimit;
  }
  result.status = std::isfinite(result.dMuMuPerS) &&
      std::isfinite(result.dDmuMuDmuPerS) && result.dMuMuPerS >= 0.0
      ? Status::Ok()
      : Status::Error(StatusCode::InvalidCoefficient,
                      "configured pitch-angle callback returned invalid output");
  return result;
}

CP::SpatialDiffusionResult EvaluateAdaptiveKappa(
    const ParticleContext& context, double speedMPerS) {
  CP::SpatialDiffusionResult result;
  CP::LocalInputView input;
  Status status = BuildSourceBoundInput(context, &input);
  if (!status.ok()) { result.status = status; return result; }
  const CP::SpeciesProperties species = SpeciesFor(context);
  status = CP::ValidateSpecies(species);
  if (!status.ok()) { result.status = status; return result; }
  const ScalarResult momentum = MomentumFromSpeed(
      speedMPerS, species.restMassKg, SpeedOfLight);
  if (!momentum.status.ok()) { result.status = momentum.status; return result; }
  CP::SpatialQuadratureConfiguration quadrature =
      Coefficient::ActiveConfiguration().spatialQuadrature;
  quadrature.gapPolicy =
      Coefficient::ActiveConfiguration().resonanceGapPolicy;
  result = CP::IntegrateSpatialDiffusion(
      speedMPerS,
      [&](double mu) {
        return EvaluatePitchKernel(context, input, species, momentum.value, mu);
      }, quadrature);
  gQuadratureEvaluations.fetch_add(result.evaluations,
                                   std::memory_order_relaxed);
  return result;
}

MeanFreePathSample EvaluateMeanFreePath(const ParticleContext& context,
                                        double momentumKgMPerS,
                                        const std::string& identity) {
  MeanFreePathSample sample;
  sample.turbulenceStateIdentity = identity;
  const Status sourceStatus = ValidateConfiguredSource();
  if (!sourceStatus.ok()) { sample.status = sourceStatus; return sample; }

  CP::LocalInputView input;
  Status status = BuildSourceBoundInput(context, &input);
  if (!status.ok()) { sample.status = status; return sample; }
  CopyIdentity(input, &sample);
  const CP::SpeciesProperties species = SpeciesFor(context);
  status = CP::ValidateSpecies(species);
  if (!status.ok()) { sample.status = status; return sample; }
  const ScalarResult speed = SpeedFromMomentum(
      momentumKgMPerS, species.restMassKg, SpeedOfLight);
  if (!speed.status.ok()) { sample.status = speed.status; return sample; }
  if (speed.value == 0.0) {
    sample.valueState = CP::ValueState::Ballistic;
    sample.lambdaParallelM = std::numeric_limits<double>::infinity();
    sample.status = Status::Ok();
    return sample;
  }

  const Coefficient::Configuration& configuration =
      Coefficient::ActiveConfiguration();
  switch (configuration.meanFreePath) {
    case Coefficient::MeanFreePathKind::Qlt: {
      CP::SpatialQuadratureConfiguration q = configuration.spatialQuadrature;
      q.gapPolicy = configuration.resonanceGapPolicy;
      const CP::SpatialDiffusionResult kappa = CP::IntegrateSpatialDiffusion(
          speed.value,
          [&](double mu) {
            return CP::EvaluateJokipiiSlab(
                input, configuration.spectrum, species, speed.value, mu);
          }, q);
      gQuadratureEvaluations.fetch_add(kappa.evaluations,
                                       std::memory_order_relaxed);
      if (!kappa.status.ok()) { sample.status = kappa.status; break; }
      if (kappa.valueState == CP::ValueState::Ballistic) {
        sample.valueState = CP::ValueState::Ballistic;
        sample.lambdaParallelM = std::numeric_limits<double>::infinity();
      }
      else {
        const ScalarResult lambda = Coefficient::MeanFreePathFromKappa(
            kappa.kappaParallelM2PerS, speed.value);
        if (!lambda.status.ok()) { sample.status = lambda.status; break; }
        sample.valueState = CP::ValueState::Finite;
        sample.lambdaParallelM = lambda.value;
      }
      sample.status = Status::Ok();
      break;
    }
    case Coefficient::MeanFreePathKind::Qlt1: {
      const CP::MeanFreePathResult lambda =
          CP::EvaluateCorrelationMeanFreePath(
              input, species, momentumKgMPerS,
              configuration.correlationLengthAt1AuM,
              configuration.spectrum.referenceRadiusM);
      sample.status = lambda.status;
      sample.valueState = lambda.valueState;
      sample.lambdaParallelM = lambda.lambdaParallelM;
      break;
    }
    case Coefficient::MeanFreePathKind::Tenishev2005: {
      const double energyJ = Relativistic::Speed2E(
          speed.value, species.restMassKg);
      sample.lambdaParallelM = SEP::Scattering::Tenishev2005AIAA::lambda0 *
          std::pow(energyJ / GeV2J,
                   SEP::Scattering::Tenishev2005AIAA::alpha) *
          std::pow(input.heliocentricRadiusM / _AU_,
                   SEP::Scattering::Tenishev2005AIAA::beta);
      sample.valueState = CP::ValueState::Finite;
      sample.status = Status::Ok();
      break;
    }
    case Coefficient::MeanFreePathKind::Chen2024: {
      const double energyJ = Relativistic::Speed2E(
          speed.value, species.restMassKg);
      const double kappaM2PerS = SEP::Diffusion::Chen2024AA::GetDxx(
          input.heliocentricRadiusM, energyJ);
      const ScalarResult lambda = Coefficient::MeanFreePathFromKappa(
          kappaM2PerS, speed.value);
      sample.status = lambda.status;
      sample.valueState = lambda.status.ok()
          ? CP::ValueState::Finite : CP::ValueState::Invalid;
      sample.lambdaParallelM = lambda.value;
      break;
    }
    case Coefficient::MeanFreePathKind::FromSpatial: {
      const CP::SpatialDiffusionResult kappa =
          EvaluateAdaptiveKappa(context, speed.value);
      if (!kappa.status.ok()) { sample.status = kappa.status; break; }
      if (kappa.valueState == CP::ValueState::Ballistic) {
        sample.valueState = CP::ValueState::Ballistic;
        sample.lambdaParallelM = std::numeric_limits<double>::infinity();
        sample.status = Status::Ok();
      }
      else {
        const ScalarResult lambda = Coefficient::MeanFreePathFromKappa(
            kappa.kappaParallelM2PerS, speed.value);
        sample.status = lambda.status;
        sample.valueState = lambda.status.ok()
            ? CP::ValueState::Finite : CP::ValueState::Invalid;
        sample.lambdaParallelM = lambda.value;
      }
      break;
    }
  }

  sample.provenance = Provenance(
      "lambda_parallel",
      Coefficient::MeanFreePathName(configuration.meanFreePath), input);
  if (sample.status.ok() &&
      ((sample.valueState == CP::ValueState::Finite &&
        std::isfinite(sample.lambdaParallelM) && sample.lambdaParallelM > 0.0) ||
       (sample.valueState == CP::ValueState::Ballistic &&
        std::isinf(sample.lambdaParallelM) && sample.lambdaParallelM > 0.0))) {
    return sample;
  }

  gInvalidSamples.fetch_add(1, std::memory_order_relaxed);
  // An explicitly rejected resonance gap is not an invalid numerical sample;
  // it is a declared unsupported physics regime and cannot be overridden by
  // the unrelated invalid-value policy.
  if (sample.status.code != StatusCode::UnresolvedCoefficient &&
      configuration.invalidPolicy == Coefficient::InvalidPolicy::Ballistic) {
    sample.valueState = CP::ValueState::Ballistic;
    sample.lambdaParallelM = std::numeric_limits<double>::infinity();
    sample.provenance += ";invalid-policy=ballistic";
    gBallisticSubstitutions.fetch_add(1, std::memory_order_relaxed);
    sample.status = Status::Ok();
  }
  else if (sample.status.ok()) {
    sample.valueState = CP::ValueState::Invalid;
    sample.status = Status::Error(StatusCode::InvalidCoefficient,
        "mean-free-path provider returned an invalid typed value");
  }
  return sample;
}

bool FiniteKappaAt(const ParticleContext& origin, double displacementM,
                   double speedMPerS, double* value,
                   std::size_t* evaluations) {
  Status locationStatus;
  const ParticleContext shifted = AtRelativeArcLength(
      origin, displacementM, &locationStatus);
  if (!locationStatus.ok()) return false;
  CP::SpatialDiffusionResult kappa;
  if (Coefficient::ActiveConfiguration().spatial ==
      Coefficient::SpatialKind::FromPitchAngle) {
    kappa = EvaluateAdaptiveKappa(shifted, speedMPerS);
  }
  else {
    const ScalarResult momentum = MomentumFromSpeed(
        speedMPerS, shifted.state.massKg, SpeedOfLight);
    if (!momentum.status.ok()) return false;
    const MeanFreePathSample lambda = EvaluateMeanFreePath(
        shifted, momentum.value, "spatial-derivative");
    if (!lambda.status.ok() || lambda.valueState != CP::ValueState::Finite)
      return false;
    const ScalarResult converted = Coefficient::KappaFromMeanFreePath(
        lambda.lambdaParallelM, speedMPerS);
    if (!converted.status.ok()) return false;
    kappa.status = Status::Ok();
    kappa.valueState = CP::ValueState::Finite;
    kappa.kappaParallelM2PerS = converted.value;
  }
  if (evaluations) *evaluations += kappa.evaluations;
  if (!kappa.status.ok() || kappa.valueState != CP::ValueState::Finite ||
      !std::isfinite(kappa.kappaParallelM2PerS)) return false;
  *value = kappa.kappaParallelM2PerS;
  return true;
}

}  // namespace

CoefficientProviderDiagnostics GetCoefficientProviderDiagnostics() {
  CoefficientProviderDiagnostics result;
  result.invalidSamples = gInvalidSamples.load(std::memory_order_relaxed);
  result.ballisticSubstitutions =
      gBallisticSubstitutions.load(std::memory_order_relaxed);
  result.amplitudeRegularizations =
      gAmplitudeRegularizations.load(std::memory_order_relaxed);
  result.quadratureEvaluations =
      gQuadratureEvaluations.load(std::memory_order_relaxed);
  return result;
}

void ResetCoefficientProviderDiagnostics() {
  gInvalidSamples.store(0, std::memory_order_relaxed);
  gBallisticSubstitutions.store(0, std::memory_order_relaxed);
  gAmplitudeRegularizations.store(0, std::memory_order_relaxed);
  gQuadratureEvaluations.store(0, std::memory_order_relaxed);
}

PICMeanFreePathProvider::PICMeanFreePathProvider(
    const ParticleContext& context, const std::string& turbulenceIdentity)
    : context_(context), turbulenceIdentity_(turbulenceIdentity) {}

MeanFreePathSample PICMeanFreePathProvider::Evaluate(
    double sM, double momentumKgMPerS, double) const {
  Status locationStatus;
  const ParticleContext sampled = AtRelativeArcLength(
      context_, sM, &locationStatus);
  if (!locationStatus.ok()) {
    MeanFreePathSample result;
    result.status = locationStatus;
    return result;
  }
  MeanFreePathSample result = EvaluateMeanFreePath(
      sampled, momentumKgMPerS, turbulenceIdentity_);
  if (result.status.ok()) {
    const ScalarResult speed = SpeedFromMomentum(
        momentumKgMPerS, sampled.state.massKg, SpeedOfLight);
    if (speed.status.ok()) {
      const double total = result.valueState == CP::ValueState::Ballistic
          ? 0.0 : speed.value / result.lambdaParallelM;
      double plusFraction = 0.5;
      CP::LocalInputView input;
      if (BuildSourceBoundInput(sampled, &input).ok() &&
          input.deltaB2T2 > 0.0) {
        plusFraction = input.deltaBPlus2T2 / input.deltaB2T2;
      }
      result.nuPlusPerS = plusFraction * total;
      result.nuMinusPerS = (1.0 - plusFraction) * total;
      result.hasBranchResolvedRates = true;
    }
  }
  return result;
}

PICSpatialDiffusionProvider::PICSpatialDiffusionProvider(
    const ParticleContext& context, const std::string& turbulenceIdentity)
    : context_(context), turbulenceIdentity_(turbulenceIdentity) {}

SpatialDiffusionSample PICSpatialDiffusionProvider::Evaluate(
    double sM, double speedMPerS) const {
  SpatialDiffusionSample sample;
  Status locationStatus;
  const ParticleContext sampled = AtRelativeArcLength(
      context_, sM, &locationStatus);
  if (!locationStatus.ok()) { sample.status = locationStatus; return sample; }
  const Status sourceStatus = ValidateConfiguredSource();
  if (!sourceStatus.ok()) { sample.status = sourceStatus; return sample; }

  if (Coefficient::ActiveConfiguration().spatial ==
          Coefficient::SpatialKind::FromPitchAngle &&
      Coefficient::ActiveConfiguration().pitchAngle ==
          Coefficient::PitchAngleKind::Configured &&
      SEP::Diffusion::GetPitchAngleDiffusionCoefficient == NULL) {
    // A null legacy callback is the explicit manufactured no-diffusion model,
    // not a resonance gap.  Preserve kappa=0 as a finite coefficient so
    // deterministic Parker tests do not reinterpret it as ballistic transport.
    sample.valueState = CP::ValueState::Finite;
    sample.kappaParallelM2PerS = 0.0;
    sample.dKappaParallelDsMPerS = 0.0;
    CP::LocalInputView input;
    const Status inputStatus = BuildSourceBoundInput(sampled, &input);
    if (!inputStatus.ok()) { sample.status = inputStatus; return sample; }
    sample.provenance = Provenance("kappa_parallel", "from-dmumu:none", input);
    sample.status = Status::Ok();
    return sample;
  }

  std::size_t evaluations = 0;
  if (!FiniteKappaAt(sampled, 0.0, speedMPerS,
                     &sample.kappaParallelM2PerS, &evaluations)) {
    sample.status = Status::Error(StatusCode::UnresolvedCoefficient,
        "spatial diffusion is non-finite or ballistic; the Parker mover has "
        "no infinite-kappa numerical operator");
    return sample;
  }
  sample.valueState = CP::ValueState::Finite;

  // Refine a physical arc-length derivative until two nested stencils agree
  // with the configured deterministic error target.  Centered stencils are
  // preferred; a boundary automatically selects the available one-sided form.
  const NumericalTolerances& tolerances = ActiveNumericalTolerances();
  double hM = tolerances.geometryFraction * sampled.segment->GetLength();
  bool derivativeAccepted = false;
  for (int refinement = 0; refinement < 10; ++refinement) {
    double plus = 0.0, minus = 0.0, plusHalf = 0.0, minusHalf = 0.0;
    const bool havePlus = FiniteKappaAt(
        sampled, hM, speedMPerS, &plus, &evaluations);
    const bool haveMinus = FiniteKappaAt(
        sampled, -hM, speedMPerS, &minus, &evaluations);
    const bool havePlusHalf = FiniteKappaAt(
        sampled, 0.5 * hM, speedMPerS, &plusHalf, &evaluations);
    const bool haveMinusHalf = FiniteKappaAt(
        sampled, -0.5 * hM, speedMPerS, &minusHalf, &evaluations);
    double coarse = 0.0, fine = 0.0;
    bool formed = false;
    if (havePlus && haveMinus && havePlusHalf && haveMinusHalf) {
      coarse = (plus - minus) / (2.0 * hM);
      fine = (plusHalf - minusHalf) / hM;
      formed = true;
    }
    else if (havePlus && havePlusHalf) {
      coarse = (plus - sample.kappaParallelM2PerS) / hM;
      fine = (plusHalf - sample.kappaParallelM2PerS) / (0.5 * hM);
      formed = true;
    }
    else if (haveMinus && haveMinusHalf) {
      coarse = (sample.kappaParallelM2PerS - minus) / hM;
      fine = (sample.kappaParallelM2PerS - minusHalf) / (0.5 * hM);
      formed = true;
    }
    if (formed && std::isfinite(fine)) {
      const StepDoublingEstimate error = EstimateStepDoublingError(
          coarse, fine, 1.0e-12,
          tolerances.deterministicRelativeTolerance);
      if (error.status.ok() && error.accepted) {
        sample.dKappaParallelDsMPerS = fine;
        derivativeAccepted = true;
        break;
      }
    }
    hM *= 0.5;
  }
  if (!derivativeAccepted) {
    sample.status = Status::Error(StatusCode::UnresolvedCoefficient,
        "along-field kappa derivative did not meet the configured refinement tolerance");
    return sample;
  }

  CP::LocalInputView input;
  const Status inputStatus = BuildSourceBoundInput(sampled, &input);
  if (!inputStatus.ok()) { sample.status = inputStatus; return sample; }
  sample.estimatedAbsoluteErrorM2PerS = 0.0;
  sample.provenance = Provenance(
      "kappa_parallel",
      Coefficient::SpatialName(Coefficient::ActiveConfiguration().spatial),
      input);
  sample.status = Status::Ok();
  return sample;
}

PICPitchAngleDiffusionProvider::PICPitchAngleDiffusionProvider(
    const ParticleContext& context, const std::string& turbulenceIdentity)
    : context_(context), turbulenceIdentity_(turbulenceIdentity) {}

PitchAngleDiffusionSample PICPitchAngleDiffusionProvider::Evaluate(
    double sM, double momentumKgMPerS, double mu) const {
  PitchAngleDiffusionSample sample;
  sample.turbulenceStateIdentity = turbulenceIdentity_;
  if (!std::isfinite(mu) || std::fabs(mu) > 1.0) {
    sample.status = Status::Error(StatusCode::InvalidArgument,
        "pitch-angle provider refuses evaluation outside |mu|<=1");
    return sample;
  }
  const Status sourceStatus = ValidateConfiguredSource();
  if (!sourceStatus.ok()) { sample.status = sourceStatus; return sample; }
  Status locationStatus;
  const ParticleContext sampled = AtRelativeArcLength(
      context_, sM, &locationStatus);
  if (!locationStatus.ok()) { sample.status = locationStatus; return sample; }
  CP::LocalInputView input;
  Status status = BuildSourceBoundInput(sampled, &input);
  if (!status.ok()) { sample.status = status; return sample; }
  const CP::SpeciesProperties species = SpeciesFor(sampled);
  status = CP::ValidateSpecies(species);
  if (!status.ok()) { sample.status = status; return sample; }
  const CP::PitchAngleResult result = EvaluatePitchKernel(
      sampled, input, species, momentumKgMPerS, mu);
  sample.status = result.status;
  sample.valueState = result.valueState;
  sample.derivativeState = result.derivativeState;
  sample.dMuMuPerS = result.dMuMuPerS;
  sample.dDmuMuDmuPerS = result.dDmuMuDmuPerS;
  CopyIdentity(input, &sample);
  sample.provenance = Provenance(
      "Dmumu", Coefficient::PitchAngleName(
                    Coefficient::ActiveConfiguration().pitchAngle), input);
  return sample;
}

}  // namespace PICAdapter
}  // namespace Transport
}  // namespace SEP
