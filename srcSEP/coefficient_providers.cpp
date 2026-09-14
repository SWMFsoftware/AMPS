#include "coefficient_providers.h"

#include <algorithm>
#include <atomic>
#include <cmath>
#include <limits>
#include <sstream>

namespace SEP {
namespace Transport {
namespace PICAdapter {
namespace {

std::atomic<std::uint64_t> gInvalidSamples(0);
std::atomic<std::uint64_t> gBallisticSubstitutions(0);

Status ValidateConfiguredSource() {
  const Background::BackgroundSnapshot& snapshot =
      Background::SnapshotStore::Instance().AcquireForMover();
  return Coefficient::ValidateSourceAgainstBackground(
      Coefficient::ActiveConfiguration().source, snapshot.provider(),
      snapshot.ownership());
}

std::string Provenance(const char* quantity, const char* provider) {
  const Background::BackgroundSnapshot& snapshot =
      Background::SnapshotStore::Instance().AcquireForMover();
  std::ostringstream out;
  out << quantity << '=' << provider
      << ";source="
      << Coefficient::SourceName(Coefficient::ActiveConfiguration().source)
      << ";background=" << Background::ProviderName(snapshot.provider())
      << ";fingerprint=" << snapshot.configuration_fingerprint()
      << ";generation=" << snapshot.field_line_generation();
  return out.str();
}

// Evaluate one of the historical analytical lambda models at the particle's
// actual field-line location.  The routine returns the raw model value; policy
// for NaN, zero, or negative output is applied once by the public provider.
double RawMeanFreePath(const ParticleContext& context,
                       double momentumKgMPerS) {
  const ScalarResult speed = SpeedFromMomentum(
      momentumKgMPerS, context.state.massKg, SpeedOfLight);
  if (!speed.status.ok()) return std::numeric_limits<double>::quiet_NaN();
  if (speed.value == 0.0) return std::numeric_limits<double>::infinity();

  double x[3] = {0.0, 0.0, 0.0};
  context.segment->GetCartesian(x, context.state.coordinate);
  const double radiusM = Vector3D::Length(x);
  const double absBT = SEP::FieldLineData::GetAbsB(
      context.state.coordinate, context.segment, context.state.fieldLineId);
  const Coefficient::MeanFreePathKind kind =
      Coefficient::ActiveConfiguration().meanFreePath;

  switch (kind) {
    case Coefficient::MeanFreePathKind::Qlt:
      return QLT::calculateMeanFreePath(radiusM, speed.value);
    case Coefficient::MeanFreePathKind::Qlt1:
      return QLT1::calculateMeanFreePath(radiusM, speed.value, absBT);
    case Coefficient::MeanFreePathKind::Tenishev2005: {
      const double energyJ = Relativistic::Speed2E(
          speed.value, context.state.massKg);
      return SEP::Scattering::Tenishev2005AIAA::lambda0 *
          std::pow(energyJ / GeV2J,
                   SEP::Scattering::Tenishev2005AIAA::alpha) *
          std::pow(radiusM / _AU_,
                   SEP::Scattering::Tenishev2005AIAA::beta);
    }
    case Coefficient::MeanFreePathKind::Chen2024: {
      const double energyJ = Relativistic::Speed2E(
          speed.value, context.state.massKg);
      const double kappaM2PerS =
          SEP::Diffusion::Chen2024AA::GetDxx(radiusM, energyJ);
      const ScalarResult lambda = Coefficient::MeanFreePathFromKappa(
          kappaM2PerS, speed.value);
      return lambda.status.ok()
          ? lambda.value : std::numeric_limits<double>::quiet_NaN();
    }
    case Coefficient::MeanFreePathKind::FromSpatial: {
      // This route intentionally calls the legacy Dmumu integral directly,
      // not PICSpatialDiffusionProvider, so the registry cannot recurse.
      double kappaM2PerS = 0.0;
      double derivativeMPerS = 0.0;
      SEP::Diffusion::GetDxx(kappaM2PerS, derivativeMPerS, speed.value,
                             context.state.species,
                             context.state.coordinate, context.segment,
                             context.state.fieldLineId);
      const ScalarResult lambda = Coefficient::MeanFreePathFromKappa(
          kappaM2PerS, speed.value);
      return lambda.status.ok()
          ? lambda.value : std::numeric_limits<double>::quiet_NaN();
    }
  }
  return std::numeric_limits<double>::quiet_NaN();
}

MeanFreePathSample EvaluateMeanFreePath(const ParticleContext& context,
                                        double momentumKgMPerS,
                                        const std::string& identity) {
  MeanFreePathSample sample;
  sample.turbulenceStateIdentity = identity;
  const Status sourceStatus = ValidateConfiguredSource();
  if (!sourceStatus.ok()) {
    sample.status = sourceStatus;
    return sample;
  }

  sample.lambdaParallelM = RawMeanFreePath(context, momentumKgMPerS);
  sample.provenance = Provenance(
      "lambda_parallel",
      Coefficient::MeanFreePathName(
          Coefficient::ActiveConfiguration().meanFreePath));
  if ((std::isfinite(sample.lambdaParallelM) ||
       std::isinf(sample.lambdaParallelM)) && sample.lambdaParallelM > 0.0) {
    sample.status = Status::Ok();
    return sample;
  }

  gInvalidSamples.fetch_add(1, std::memory_order_relaxed);
  if (Coefficient::ActiveConfiguration().invalidPolicy ==
      Coefficient::InvalidPolicy::Ballistic) {
    // Positive infinity is the exact zero-event-rate state under nu=v/lambda.
    // It is deliberately not a large finite clamp, which would introduce an
    // undocumented nonzero scattering probability.
    sample.lambdaParallelM = std::numeric_limits<double>::infinity();
    sample.provenance += ";invalid-policy=ballistic";
    gBallisticSubstitutions.fetch_add(1, std::memory_order_relaxed);
    sample.status = Status::Ok();
  }
  else {
    sample.status = Status::Error(StatusCode::InvalidCoefficient,
        "mean-free-path provider returned NaN or a non-positive value");
  }
  return sample;
}

}  // namespace

CoefficientProviderDiagnostics GetCoefficientProviderDiagnostics() {
  CoefficientProviderDiagnostics result;
  result.invalidSamples = gInvalidSamples.load(std::memory_order_relaxed);
  result.ballisticSubstitutions =
      gBallisticSubstitutions.load(std::memory_order_relaxed);
  return result;
}

void ResetCoefficientProviderDiagnostics() {
  gInvalidSamples.store(0, std::memory_order_relaxed);
  gBallisticSubstitutions.store(0, std::memory_order_relaxed);
}

PICMeanFreePathProvider::PICMeanFreePathProvider(
    const ParticleContext& context, const std::string& turbulenceIdentity)
    : context_(context), turbulenceIdentity_(turbulenceIdentity) {}

MeanFreePathSample PICMeanFreePathProvider::Evaluate(
    double, double momentumKgMPerS, double) const {
  return EvaluateMeanFreePath(context_, momentumKgMPerS,
                              turbulenceIdentity_);
}

PICSpatialDiffusionProvider::PICSpatialDiffusionProvider(
    const ParticleContext& context, const std::string& turbulenceIdentity)
    : context_(context), turbulenceIdentity_(turbulenceIdentity) {}

SpatialDiffusionSample PICSpatialDiffusionProvider::Evaluate(
    double, double speedMPerS) const {
  SpatialDiffusionSample sample;
  const Status sourceStatus = ValidateConfiguredSource();
  if (!sourceStatus.ok()) {
    sample.status = sourceStatus;
    return sample;
  }

  if (Coefficient::ActiveConfiguration().spatial ==
      Coefficient::SpatialKind::FromPitchAngle) {
    if (SEP::Diffusion::GetPitchAngleDiffusionCoefficient == NULL) {
      // A null callback is an intentional deterministic-transport model used
      // by manufactured tests, not an invalid-value recovery path.
      sample.provenance = Provenance("kappa_parallel", "from-dmumu:none");
      sample.status = Status::Ok();
      return sample;
    }
    SEP::Diffusion::GetDxx(sample.kappaParallelM2PerS,
                           sample.dKappaParallelDsMPerS, speedMPerS,
                           context_.state.species, context_.state.coordinate,
                           context_.segment, context_.state.fieldLineId);
  }
  else {
    const ScalarResult momentum = MomentumFromSpeed(
        speedMPerS, context_.state.massKg, SpeedOfLight);
    if (!momentum.status.ok()) {
      sample.status = momentum.status;
      return sample;
    }
    const MeanFreePathSample lambda = EvaluateMeanFreePath(
        context_, momentum.value, turbulenceIdentity_);
    if (!lambda.status.ok()) {
      sample.status = lambda.status;
      return sample;
    }
    const ScalarResult kappa = Coefficient::KappaFromMeanFreePath(
        lambda.lambdaParallelM, speedMPerS);
    if (!kappa.status.ok()) {
      sample.status = kappa.status;
      return sample;
    }
    sample.kappaParallelM2PerS = kappa.value;

    // The selected analytical lambda models depend on heliocentric position.
    // Use a symmetric physical-distance stencil along the field line; if one
    // side reaches a boundary, the valid one-sided stencil remains explicit.
    const double hM = 0.05 * context_.segment->GetLength();
    bool haveMinus = false, havePlus = false;
    double kappaMinus = 0.0, kappaPlus = 0.0;
    for (int sign = -1; sign <= 1; sign += 2) {
      ParticleContext shifted = context_;
      shifted.state.coordinate =
          PIC::FieldLine::FieldLinesAll[shifted.state.fieldLineId].move(
              shifted.state.coordinate, sign * hM, shifted.segment);
      shifted.segment =
          PIC::FieldLine::FieldLinesAll[shifted.state.fieldLineId].GetSegment(
              shifted.state.coordinate);
      if (!shifted.segment) continue;
      const MeanFreePathSample shiftedLambda = EvaluateMeanFreePath(
          shifted, momentum.value, turbulenceIdentity_);
      const ScalarResult shiftedKappa = Coefficient::KappaFromMeanFreePath(
          shiftedLambda.lambdaParallelM, speedMPerS);
      if (!shiftedLambda.status.ok() || !shiftedKappa.status.ok()) continue;
      if (sign < 0) {
        haveMinus = true;
        kappaMinus = shiftedKappa.value;
      }
      else {
        havePlus = true;
        kappaPlus = shiftedKappa.value;
      }
    }
    if (haveMinus && havePlus)
      sample.dKappaParallelDsMPerS = (kappaPlus - kappaMinus) / (2.0 * hM);
    else if (havePlus)
      sample.dKappaParallelDsMPerS =
          (kappaPlus - sample.kappaParallelM2PerS) / hM;
    else if (haveMinus)
      sample.dKappaParallelDsMPerS =
          (sample.kappaParallelM2PerS - kappaMinus) / hM;
    else {
      sample.status = Status::Error(StatusCode::OutOfDomain,
          "cannot form an along-field kappa derivative at this location");
      return sample;
    }
  }

  if (!std::isfinite(sample.kappaParallelM2PerS) ||
      !std::isfinite(sample.dKappaParallelDsMPerS) ||
      sample.kappaParallelM2PerS < 0.0) {
    sample.status = Status::Error(StatusCode::InvalidCoefficient,
        "spatial provider returned invalid SI coefficients");
    return sample;
  }
  sample.provenance = Provenance(
      "kappa_parallel",
      Coefficient::SpatialName(Coefficient::ActiveConfiguration().spatial));
  sample.status = Status::Ok();
  return sample;
}

PICPitchAngleDiffusionProvider::PICPitchAngleDiffusionProvider(
    const ParticleContext& context, const std::string& turbulenceIdentity)
    : context_(context), turbulenceIdentity_(turbulenceIdentity) {}

PitchAngleDiffusionSample PICPitchAngleDiffusionProvider::Evaluate(
    double, double momentumKgMPerS, double mu) const {
  PitchAngleDiffusionSample sample;
  sample.turbulenceStateIdentity = turbulenceIdentity_;
  const Status sourceStatus = ValidateConfiguredSource();
  if (!sourceStatus.ok()) {
    sample.status = sourceStatus;
    return sample;
  }
  if (SEP::Diffusion::GetPitchAngleDiffusionCoefficient == NULL) {
    sample.provenance = Provenance("Dmumu", "configured:none");
    sample.status = Status::Ok();
    return sample;
  }

  const ScalarResult speed = SpeedFromMomentum(
      momentumKgMPerS, context_.state.massKg, SpeedOfLight);
  if (!speed.status.ok()) {
    sample.status = speed.status;
    return sample;
  }
  const double boundedMu = std::max(-1.0, std::min(1.0, mu));
  const double vParallel = speed.value * boundedMu;
  const double vNormal = speed.value *
      std::sqrt(std::max(0.0, 1.0 - boundedMu * boundedMu));
  SEP::Diffusion::GetPitchAngleDiffusionCoefficient(
      sample.dMuMuPerS, sample.dDmuMuDmuPerS, boundedMu,
      vParallel, vNormal, context_.state.species,
      context_.state.coordinate, context_.segment);

  if (SEP::Diffusion::PitchAngleDifferentialMode ==
      SEP::Diffusion::PitchAngleDifferentialModeNumerical) {
    const double h = SEP::Diffusion::muNumericalDifferentiationStep;
    if (!std::isfinite(h) || h <= 0.0) {
      sample.status = Status::Error(StatusCode::InvalidCoefficient,
                                    "Dmumu derivative step must be positive");
      return sample;
    }
    const double muMinus = std::max(-1.0, boundedMu - h);
    const double muPlus = std::min(1.0, boundedMu + h);
    double dMinus = 0.0, unusedMinus = 0.0;
    double dPlus = 0.0, unusedPlus = 0.0;
    SEP::Diffusion::GetPitchAngleDiffusionCoefficient(
        dMinus, unusedMinus, muMinus, speed.value * muMinus,
        speed.value * std::sqrt(std::max(0.0, 1.0 - muMinus * muMinus)),
        context_.state.species, context_.state.coordinate, context_.segment);
    SEP::Diffusion::GetPitchAngleDiffusionCoefficient(
        dPlus, unusedPlus, muPlus, speed.value * muPlus,
        speed.value * std::sqrt(std::max(0.0, 1.0 - muPlus * muPlus)),
        context_.state.species, context_.state.coordinate, context_.segment);
    if (!(muPlus > muMinus)) {
      sample.status = Status::Error(StatusCode::InvalidCoefficient,
                                    "Dmumu derivative stencil collapsed");
      return sample;
    }
    sample.dDmuMuDmuPerS = (dPlus - dMinus) / (muPlus - muMinus);
  }

  if (!std::isfinite(sample.dMuMuPerS) ||
      !std::isfinite(sample.dDmuMuDmuPerS) || sample.dMuMuPerS < 0.0) {
    sample.status = Status::Error(StatusCode::InvalidCoefficient,
        "pitch-angle provider returned an invalid SI coefficient");
    return sample;
  }
  sample.provenance = Provenance(
      "Dmumu", Coefficient::PitchAngleName(
                    Coefficient::ActiveConfiguration().pitchAngle));
  sample.status = Status::Ok();
  return sample;
}

}  // namespace PICAdapter
}  // namespace Transport
}  // namespace SEP
