#include "turbulence_models.h"

// The coefficient declarations and sep_common headers are implementation-only
// dependencies.  Keeping them out of turbulence_models.h lets AMPS compile the
// provider-facing main_lib.cpp without needing an application-specific -I flag.
#include "coefficient_bridge.h"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <limits>
#include <sstream>

namespace SEP3D {
namespace Turbulence {
namespace {

constexpr double kMu0 = 4.0e-7 * Core::Const::kPi;

Core::Status Invalid(const std::string& message) {
  return Core::Status(Core::StatusCode::InvalidInput, message);
}

Core::Status Unavailable(const std::string& message) {
  return Core::Status(Core::StatusCode::SnapshotUnavailable, message);
}

std::uint64_t Digest64(const std::string& text) {
  // FNV-1a is a deterministic provenance tag, not a security hash.  The
  // manifest string remains the authoritative human-readable configuration.
  std::uint64_t result = UINT64_C(1469598103934665603);
  for (unsigned char value : text) {
    result ^= value;
    result *= UINT64_C(1099511628211);
  }
  return result;
}

std::string HexDigest(const std::string& text) {
  std::ostringstream output;
  output << std::hex << std::setw(16) << std::setfill('0') << Digest64(text);
  return output.str();
}

bool SamePoint(const Core::Vec3& left, const Core::Vec3& right) {
  const double scale = std::max(1.0, std::max(left.Norm(), right.Norm()));
  return (left - right).Norm() <= 1.0e-12 * scale;
}

bool FinitePositive(double value) {
  return std::isfinite(value) && value > 0.0;
}

TurbulenceSample MissingSample(MissingTurbulencePolicy policy,
                               const std::string& message,
                               const TurbulenceMetadata& metadata) {
  TurbulenceSample sample;
  sample.sourceIdentity = metadata.providerIdentity;
  sample.generation = metadata.generation;
  sample.configurationDigest = Digest64(metadata.configurationFingerprint);
  if (policy == MissingTurbulencePolicy::Ballistic) {
    // Ballistic is a valid zero-scattering state only because the user chose
    // it explicitly.  The dedicated flag prevents zero variance from being
    // confused with an ordinary finite turbulence record.
    sample.status = Core::Status::Ballistic(message);
    sample.valid = true;
    sample.ballistic = true;
  } else {
    sample.status = Unavailable(message);
  }
  return sample;
}

SEP::Transport::CoefficientPhysics::PitchAngleResult CoefficientError(
    const std::string& message) {
  SEP::Transport::CoefficientPhysics::PitchAngleResult result;
  result.status = SEP::Transport::Status::Error(
      SEP::Transport::StatusCode::InvalidCoefficient, message);
  result.valueState =
      SEP::Transport::CoefficientPhysics::ValueState::Unavailable;
  return result;
}

}  // namespace

PrescribedKolmogorovProvider::PrescribedKolmogorovProvider(
    const PrescribedKolmogorovConfiguration& configuration)
    : configuration_(configuration) {
  metadata_.source = TurbulenceSource::PrescribedKolmogorov;
  metadata_.ownership = TurbulenceOwnership::ModelOwned;
  metadata_.coordinateFrame = configuration_.coordinateFrame;
  metadata_.providerIdentity = CanonicalName();
  metadata_.configurationFingerprint = HexDigest(ResolvedManifest());
  configurationDigest_ = Digest64(metadata_.configurationFingerprint);
}

Core::Status PrescribedKolmogorovProvider::Validate() const {
  const double values[] = {
      configuration_.deltaBOverB, configuration_.referenceRadiusM,
      configuration_.kMinAtReferencePerM,
      configuration_.kMaxAtReferencePerM,
      configuration_.kMinRadialExponent,
      configuration_.kMaxRadialExponent,
      configuration_.spectralIndex,
      configuration_.parallelCorrelationLengthAtReferenceM,
      configuration_.correlationLengthRadialExponent,
      configuration_.validityCadenceS};
  for (double value : values) {
    if (!std::isfinite(value))
      return Invalid("prescribed turbulence configuration is not finite");
  }
  if (configuration_.deltaBOverB <= 0.0 ||
      configuration_.referenceRadiusM <= 0.0 ||
      configuration_.kMinAtReferencePerM <= 0.0 ||
      configuration_.kMaxAtReferencePerM <=
          configuration_.kMinAtReferencePerM ||
      configuration_.spectralIndex <= 1.0 ||
      configuration_.parallelCorrelationLengthAtReferenceM <= 0.0 ||
      configuration_.validityCadenceS <= 0.0 ||
      configuration_.coordinateFrame.empty()) {
    return Invalid("prescribed turbulence scales or spectral band are invalid");
  }
  return Core::Status::OK();
}

Core::Status PrescribedKolmogorovProvider::Prepare(double epochS) {
  const Core::Status valid = Validate();
  if (!valid.ok()) return valid;
  if (!std::isfinite(epochS)) return Invalid("turbulence epoch is not finite");

  TurbulenceMetadata candidate = metadata_;
  candidate.epochS = epochS;
  candidate.validFromS = epochS;
  candidate.validUntilS = epochS + configuration_.validityCadenceS;
  candidate.generation = metadata_.generation + 1;
  metadata_ = candidate;
  prepared_ = true;
  return Core::Status::OK();
}

const TurbulenceMetadata*
PrescribedKolmogorovProvider::PreparedMetadata() const {
  return prepared_ ? &metadata_ : nullptr;
}

TurbulenceSample PrescribedKolmogorovProvider::Evaluate(
    const Core::Vec3& positionM,
    const Background::BackgroundSample& background) const {
  TurbulenceSample sample;
  if (!prepared_)
    return MissingSample(MissingTurbulencePolicy::Fail,
                         "prescribed turbulence is not prepared", metadata_);
  if (!background.status.ok() || !background.valid ||
      !FinitePositive(background.absB)) {
    return MissingSample(MissingTurbulencePolicy::Fail,
                         "prescribed turbulence requires a valid |B|",
                         metadata_);
  }
  const double radius = positionM.Norm();
  if (!FinitePositive(radius))
    return MissingSample(MissingTurbulencePolicy::Fail,
                         "prescribed turbulence radius is invalid", metadata_);

  const double radiusRatio = configuration_.referenceRadiusM / radius;
  sample.deltaB2T2 = std::pow(configuration_.deltaBOverB * background.absB, 2);
  sample.deltaBPlus2T2 = 0.5 * sample.deltaB2T2;
  sample.deltaBMinus2T2 = 0.5 * sample.deltaB2T2;
  if (background.B.Dot(positionM) >= 0.0) {
    sample.deltaBOutward2T2 = sample.deltaBPlus2T2;
    sample.deltaBInward2T2 = sample.deltaBMinus2T2;
  } else {
    sample.deltaBOutward2T2 = sample.deltaBMinus2T2;
    sample.deltaBInward2T2 = sample.deltaBPlus2T2;
  }
  sample.kMinPerM = configuration_.kMinAtReferencePerM *
      std::pow(radiusRatio, configuration_.kMinRadialExponent);
  sample.kMaxPerM = configuration_.kMaxAtReferencePerM *
      std::pow(radiusRatio, configuration_.kMaxRadialExponent);
  sample.spectralIndex = configuration_.spectralIndex;
  sample.parallelCorrelationLengthM =
      configuration_.parallelCorrelationLengthAtReferenceM *
      std::pow(radius / configuration_.referenceRadiusM,
               configuration_.correlationLengthRadialExponent);
  sample.generation = metadata_.generation;
  sample.configurationDigest = configurationDigest_;
  sample.sourceIdentity = CanonicalName();
  sample.valid = true;
  sample.status = Core::Status::OK();
  return sample;
}

std::string PrescribedKolmogorovProvider::ResolvedManifest() const {
  std::ostringstream output;
  output << std::setprecision(17) << std::scientific
         << "prescribed-kolmogorov-v1"
         << ";deltaB_over_B=" << configuration_.deltaBOverB
         << ";reference_m=" << configuration_.referenceRadiusM
         << ";kmin_ref_m-1=" << configuration_.kMinAtReferencePerM
         << ";kmax_ref_m-1=" << configuration_.kMaxAtReferencePerM
         << ";kmin_exp=" << configuration_.kMinRadialExponent
         << ";kmax_exp=" << configuration_.kMaxRadialExponent
         << ";index=" << configuration_.spectralIndex
         << ";correlation_ref_m="
         << configuration_.parallelCorrelationLengthAtReferenceM
         << ";correlation_exp="
         << configuration_.correlationLengthRadialExponent
         << ";cadence_s=" << configuration_.validityCadenceS
         << ";frame=" << configuration_.coordinateFrame;
  return output.str();
}

AwsomTurbulenceProvider::AwsomTurbulenceProvider(
    MissingTurbulencePolicy policy)
    : missingPolicy_(policy) {}

Core::Status AwsomTurbulenceProvider::Load(
    const AwsomTurbulenceImport& imported) {
  if (!std::isfinite(imported.epochS) ||
      !std::isfinite(imported.validUntilS) ||
      imported.validUntilS < imported.epochS || imported.generation == 0 ||
      imported.coordinateFrame != "HCI-like-inertial" ||
      imported.configurationFingerprint.empty() || imported.records.empty() ||
      !FinitePositive(imported.kMinPerM) ||
      !(imported.kMaxPerM > imported.kMinPerM) ||
      !std::isfinite(imported.spectralIndex) ||
      imported.spectralIndex <= 1.0 ||
      !FinitePositive(imported.parallelCorrelationLengthM)) {
    return Invalid("AWSoM turbulence import metadata or spectrum is invalid");
  }

  for (const AwsomWaveRecord& record : imported.records) {
    if (!std::isfinite(record.positionM.x) ||
        !std::isfinite(record.positionM.y) ||
        !std::isfinite(record.positionM.z) ||
        !std::isfinite(record.epochS) || record.epochS != imported.epochS) {
      return Invalid("AWSoM wave record has an invalid position or epoch");
    }
    // An incomplete record is retained intentionally: Evaluate applies the
    // declared Fail/Ballistic policy at that point.  Complete records, in
    // contrast, must have finite non-negative energies before publication.
    if (record.complete &&
        (!std::isfinite(record.wPlusJPerM3) ||
         !std::isfinite(record.wMinusJPerM3) ||
         record.wPlusJPerM3 < 0.0 || record.wMinusJPerM3 < 0.0)) {
      return Invalid("complete AWSoM wave record has invalid energy density");
    }
  }

  TurbulenceMetadata candidate;
  candidate.source = TurbulenceSource::SwmfAwsom;
  candidate.ownership = TurbulenceOwnership::ImportedReadOnly;
  candidate.epochS = imported.epochS;
  candidate.validFromS = imported.epochS;
  candidate.validUntilS = imported.validUntilS;
  candidate.generation = imported.generation;
  candidate.coordinateFrame = imported.coordinateFrame;
  candidate.providerIdentity = CanonicalName();
  candidate.configurationFingerprint = imported.configurationFingerprint;

  imported_ = imported;
  metadata_ = candidate;
  loaded_ = true;
  prepared_ = false;
  return Core::Status::OK();
}

Core::Status AwsomTurbulenceProvider::Validate() const {
  return loaded_ ? Core::Status::OK()
                 : Unavailable("no AWSoM turbulence generation is loaded");
}

Core::Status AwsomTurbulenceProvider::Prepare(double epochS) {
  const Core::Status valid = Validate();
  if (!valid.ok()) return valid;
  if (!std::isfinite(epochS) || epochS != imported_.epochS)
    return Unavailable("requested epoch does not match the AWSoM wave epoch");
  prepared_ = true;
  return Core::Status::OK();
}

const TurbulenceMetadata* AwsomTurbulenceProvider::PreparedMetadata() const {
  return prepared_ ? &metadata_ : nullptr;
}

TurbulenceSample AwsomTurbulenceProvider::Evaluate(
    const Core::Vec3& positionM,
    const Background::BackgroundSample& background) const {
  if (!prepared_)
    return MissingSample(missingPolicy_, "AWSoM turbulence is not prepared",
                         metadata_);
  if (!background.status.ok() || !background.valid ||
      !FinitePositive(background.absB) || positionM.Norm() == 0.0) {
    return MissingSample(missingPolicy_,
                         "AWSoM direction mapping requires valid B and r",
                         metadata_);
  }

  const AwsomWaveRecord* matched = nullptr;
  for (const AwsomWaveRecord& record : imported_.records) {
    if (SamePoint(record.positionM, positionM)) {
      matched = &record;
      break;
    }
  }
  if (matched == nullptr || !matched->complete) {
    return MissingSample(missingPolicy_,
                         "AWSoM wave record is absent or incomplete",
                         metadata_);
  }

  TurbulenceSample sample;
  // AWSoM w is total Alfvén-wave energy.  Equipartition gives magnetic
  // variance deltaB^2=mu0*w, not 2*mu0*w.  The two propagation directions are
  // converted independently so cross helicity is preserved.
  sample.deltaBPlus2T2 = kMu0 * matched->wPlusJPerM3;
  sample.deltaBMinus2T2 = kMu0 * matched->wMinusJPerM3;
  sample.deltaB2T2 = sample.deltaBPlus2T2 + sample.deltaBMinus2T2;
  if (background.B.Dot(positionM) >= 0.0) {
    sample.deltaBOutward2T2 = sample.deltaBPlus2T2;
    sample.deltaBInward2T2 = sample.deltaBMinus2T2;
  } else {
    sample.deltaBOutward2T2 = sample.deltaBMinus2T2;
    sample.deltaBInward2T2 = sample.deltaBPlus2T2;
  }
  sample.kMinPerM = imported_.kMinPerM;
  sample.kMaxPerM = imported_.kMaxPerM;
  sample.spectralIndex = imported_.spectralIndex;
  sample.parallelCorrelationLengthM = imported_.parallelCorrelationLengthM;
  sample.generation = metadata_.generation;
  sample.configurationDigest = Digest64(metadata_.configurationFingerprint);
  sample.sourceIdentity = CanonicalName();
  sample.valid = true;
  sample.status = Core::Status::OK();
  return sample;
}

std::string AwsomTurbulenceProvider::ResolvedManifest() const {
  std::ostringstream output;
  output << std::setprecision(17) << std::scientific
         << "swmf-awsom-waves-v1"
         << ";policy="
         << (missingPolicy_ == MissingTurbulencePolicy::Fail
                 ? "fail" : "ballistic")
         << ";generation=" << metadata_.generation
         << ";epoch_s=" << metadata_.epochS
         << ";valid_until_s=" << metadata_.validUntilS
         << ";frame=" << metadata_.coordinateFrame
         << ";kmin_m-1=" << imported_.kMinPerM
         << ";kmax_m-1=" << imported_.kMaxPerM
         << ";index=" << imported_.spectralIndex
         << ";correlation_m=" << imported_.parallelCorrelationLengthM;
  return output.str();
}

NormalizedPowerLawSpectrum::NormalizedPowerLawSpectrum(
    const TurbulenceSample& sample)
    : varianceT2_(sample.deltaB2T2),
      kMinPerM_(sample.kMinPerM),
      kMaxPerM_(sample.kMaxPerM),
      index_(sample.spectralIndex) {
  if (Validate().ok()) {
    const double exponent = 1.0 - index_;
    normalization_ = varianceT2_ * exponent /
        (std::pow(kMaxPerM_, exponent) -
         std::pow(kMinPerM_, exponent));
  }
}

Core::Status NormalizedPowerLawSpectrum::Validate() const {
  if (!std::isfinite(varianceT2_) || varianceT2_ < 0.0 ||
      !FinitePositive(kMinPerM_) || !(kMaxPerM_ > kMinPerM_) ||
      !std::isfinite(index_) || index_ <= 1.0) {
    return Invalid("power-law spectrum requires variance>=0, ordered positive "
                   "bounds, and index>1");
  }
  return Core::Status::OK();
}

SpectrumValue NormalizedPowerLawSpectrum::Evaluate(
    double waveNumberPerM, ResonanceRangePolicy policy) const {
  SpectrumValue result;
  const Core::Status valid = Validate();
  if (!valid.ok()) {
    result.status = valid;
    return result;
  }
  if (!FinitePositive(waveNumberPerM)) {
    result.status = Invalid("resonant wave number must be finite and positive");
    return result;
  }
  const bool outside = waveNumberPerM < kMinPerM_ ||
                       waveNumberPerM > kMaxPerM_;
  if (outside && policy == ResonanceRangePolicy::Reject) {
    result.status = Core::Status(
        Core::StatusCode::BackgroundInvalid,
        "resonant wave number lies outside the declared turbulence band");
    return result;
  }
  result.resolvedWaveNumberPerM = waveNumberPerM;
  result.extended = outside;
  result.valueT2M = normalization_ * std::pow(waveNumberPerM, -index_);
  result.status = std::isfinite(result.valueT2M)
      ? Core::Status::OK()
      : Invalid("power-law spectral evaluation overflowed");
  return result;
}

double NormalizedPowerLawSpectrum::AnalyticBandVarianceT2() const {
  if (!Validate().ok()) return std::numeric_limits<double>::quiet_NaN();
  const double exponent = 1.0 - index_;
  return normalization_ / exponent *
      (std::pow(kMaxPerM_, exponent) -
       std::pow(kMinPerM_, exponent));
}

SEP::Transport::CoefficientPhysics::LocalInputView
CoefficientBridge::ToSharedInput(
    const TurbulenceSample& turbulence,
    const Background::BackgroundSample& background,
    double heliocentricRadiusM) {
  SEP::Transport::CoefficientPhysics::LocalInputView input;
  input.source = turbulence.sourceIdentity;
  input.representation = "srcsep3d-phase-t-local-view-v1";
  input.generation = turbulence.generation;
  input.checksum = turbulence.configurationDigest;
  input.heliocentricRadiusM = heliocentricRadiusM;
  input.magneticFieldT = background.absB;
  input.deltaB2T2 = turbulence.deltaB2T2;
  input.deltaBPlus2T2 = turbulence.deltaBPlus2T2;
  input.deltaBMinus2T2 = turbulence.deltaBMinus2T2;
  input.alfvenSpeedMPerS = background.alfvenSpeedMpS;
  return input;
}

SEP::Transport::CoefficientPhysics::PitchAngleResult
CoefficientBridge::JokipiiDmumu(
    const TurbulenceSample& turbulence,
    const Background::BackgroundSample& background,
    double heliocentricRadiusM,
    const SEP::Transport::CoefficientPhysics::SpectrumParameters& spectrum,
    const SEP::Transport::CoefficientPhysics::SpeciesProperties& species,
    double speedMPerS, double mu) {
  namespace CP = SEP::Transport::CoefficientPhysics;
  if (turbulence.ballistic && turbulence.status.ballistic()) {
    CP::PitchAngleResult result;
    result.status = SEP::Transport::Status::Ok();
    result.valueState = CP::ValueState::Ballistic;
    result.dMuMuPerS = 0.0;
    result.dDmuMuDmuPerS = 0.0;
    return result;
  }
  if (!turbulence.status.ok() || !turbulence.valid)
    return CoefficientError("turbulence input is unavailable");
  if (!background.status.ok() || !background.valid)
    return CoefficientError("background input is unavailable");
  return CP::EvaluateJokipiiSlab(
      ToSharedInput(turbulence, background, heliocentricRadiusM),
      spectrum, species, speedMPerS, mu);
}

SEP::Transport::ScalarResult CoefficientBridge::KappaFromMeanFreePath(
    double lambdaM, double speedMPerS) {
  return SEP::Transport::Coefficient::KappaFromMeanFreePath(
      lambdaM, speedMPerS);
}

SEP::Transport::ScalarResult CoefficientBridge::MeanFreePathFromKappa(
    double kappaM2PerS, double speedMPerS) {
  return SEP::Transport::Coefficient::MeanFreePathFromKappa(
      kappaM2PerS, speedMPerS);
}

SEP::Transport::ScalarResult
CoefficientBridge::IsotropicDmumuFromMeanFreePath(
    double lambdaM, double speedMPerS, double mu) {
  return SEP::Transport::Coefficient::IsotropicDmumuFromMeanFreePath(
      lambdaM, speedMPerS, mu);
}

SEP::Transport::ScalarResult
CoefficientBridge::MeanFreePathFromIsotropicDmumu(
    double dmumuPerS, double speedMPerS, double mu) {
  return SEP::Transport::Coefficient::MeanFreePathFromIsotropicDmumu(
      dmumuPerS, speedMPerS, mu);
}

}  // namespace Turbulence
}  // namespace SEP3D
