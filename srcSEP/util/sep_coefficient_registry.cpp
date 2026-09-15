#include "sep_coefficient_registry.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <limits>
#include <sstream>

namespace SEP {
namespace Transport {
namespace Coefficient {
namespace {

std::string Lower(std::string value) {
  std::transform(value.begin(), value.end(), value.begin(),
      [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
  return value;
}

Configuration gActiveConfiguration;

}  // namespace

const std::vector<Descriptor>& SpatialRegistry() {
  static const std::vector<Descriptor> values = {
      {"from-dmumu", "kappa_parallel and derivative", "m2/s; m/s",
       "pitch-angle provider; numerical mu quadrature"},
      {"from-mfp", "kappa_parallel and derivative", "m2/s; m/s",
       "lambda_parallel closure kappa=v*lambda/3"}};
  return values;
}

const std::vector<Descriptor>& PitchAngleRegistry() {
  static const std::vector<Descriptor> values = {
      {"configured", "Dmumu and dDmumu/dmu", "1/s; 1/s",
       "legacy prescribed callback compatibility adapter"},
      {"constant", "Dmumu and dDmumu/dmu", "1/s; 1/s",
       "constant_dmumu_per_s>=0"},
      {"jokipii-1966", "Dmumu and dDmumu/dmu", "1/s; 1/s",
       "source-bound B, deltaB2, species, and finite spectral bounds"},
      {"florinskiy", "Dmumu and dDmumu/dmu", "1/s; 1/s",
       "source-bound branch variances, Alfven speed, and correlation length"}};
  return values;
}

const std::vector<Descriptor>& MeanFreePathRegistry() {
  static const std::vector<Descriptor> values = {
      {"qlt", "lambda_parallel", "m", "Kolmogorov slab QLT parameters"},
      {"qlt1", "lambda_parallel", "m", "deltaB/B and correlation scale"},
      {"tenishev-2005", "lambda_parallel", "m",
       "lambda0_m, radial exponent beta, energy exponent alpha"},
      {"chen-2024", "lambda_parallel", "m",
       "published radial and energy power law"},
      {"from-spatial", "lambda_parallel", "m",
       "isotropic closure lambda=3*kappa/v"}};
  return values;
}

const char* SourceName(SourceMode source) {
  switch (source) {
    case SourceMode::Prescribed: return "prescribed";
    case SourceMode::SelfConsistent: return "self-consistent";
    case SourceMode::Swmf: return "swmf";
  }
  return "unknown";
}

const char* SpatialName(SpatialKind kind) {
  return kind == SpatialKind::FromPitchAngle ? "from-dmumu" : "from-mfp";
}
const char* PitchAngleName(PitchAngleKind kind) {
  switch (kind) {
    case PitchAngleKind::Configured: return "configured";
    case PitchAngleKind::Constant: return "constant";
    case PitchAngleKind::Jokipii1966: return "jokipii-1966";
    case PitchAngleKind::Florinskiy: return "florinskiy";
  }
  return "unknown";
}
const char* MeanFreePathName(MeanFreePathKind kind) {
  switch (kind) {
    case MeanFreePathKind::Qlt: return "qlt";
    case MeanFreePathKind::Qlt1: return "qlt1";
    case MeanFreePathKind::Tenishev2005: return "tenishev-2005";
    case MeanFreePathKind::Chen2024: return "chen-2024";
    case MeanFreePathKind::FromSpatial: return "from-spatial";
  }
  return "unknown";
}
const char* InvalidPolicyName(InvalidPolicy policy) {
  return policy == InvalidPolicy::Fail ? "fail" : "ballistic";
}
const char* ResonanceGapPolicyName(
    CoefficientPhysics::ResonanceGapPolicy policy) {
  return policy == CoefficientPhysics::ResonanceGapPolicy::Reject
      ? "reject" : "ballistic";
}
const char* TurbulenceAmplitudePolicyName(TurbulenceAmplitudePolicy policy) {
  return policy == TurbulenceAmplitudePolicy::Reject
      ? "reject" : "limit-to-mean-field";
}

bool ParseSource(const std::string& text, SourceMode* source) {
  if (!source) return false;
  const std::string value = Lower(text);
  if (value == "prescribed") *source = SourceMode::Prescribed;
  else if (value == "self-consistent") *source = SourceMode::SelfConsistent;
  else if (value == "swmf") *source = SourceMode::Swmf;
  else return false;
  return true;
}
bool ParseSpatial(const std::string& text, SpatialKind* kind) {
  if (!kind) return false;
  const std::string value = Lower(text);
  if (value == "from-dmumu") *kind = SpatialKind::FromPitchAngle;
  else if (value == "from-mfp") *kind = SpatialKind::FromMeanFreePath;
  else return false;
  return true;
}
bool ParsePitchAngle(const std::string& text, PitchAngleKind* kind) {
  if (!kind) return false;
  const std::string value = Lower(text);
  if (value == "configured") *kind = PitchAngleKind::Configured;
  else if (value == "constant") *kind = PitchAngleKind::Constant;
  else if (value == "jokipii-1966" || value == "jokipii")
    *kind = PitchAngleKind::Jokipii1966;
  else if (value == "florinskiy") *kind = PitchAngleKind::Florinskiy;
  else return false;
  return true;
}
bool ParseMeanFreePath(const std::string& text, MeanFreePathKind* kind) {
  if (!kind) return false;
  const std::string value = Lower(text);
  if (value == "qlt") *kind = MeanFreePathKind::Qlt;
  else if (value == "qlt1") *kind = MeanFreePathKind::Qlt1;
  else if (value == "tenishev-2005") *kind = MeanFreePathKind::Tenishev2005;
  else if (value == "chen-2024") *kind = MeanFreePathKind::Chen2024;
  else if (value == "from-spatial") *kind = MeanFreePathKind::FromSpatial;
  else return false;
  return true;
}
bool ParseInvalidPolicy(const std::string& text, InvalidPolicy* policy) {
  if (!policy) return false;
  const std::string value = Lower(text);
  if (value == "fail") *policy = InvalidPolicy::Fail;
  else if (value == "ballistic") *policy = InvalidPolicy::Ballistic;
  else return false;
  return true;
}

bool ParseResonanceGapPolicy(
    const std::string& text, CoefficientPhysics::ResonanceGapPolicy* policy) {
  if (!policy) return false;
  const std::string value = Lower(text);
  if (value == "reject")
    *policy = CoefficientPhysics::ResonanceGapPolicy::Reject;
  else if (value == "ballistic")
    *policy = CoefficientPhysics::ResonanceGapPolicy::Ballistic;
  else return false;
  return true;
}

bool ParseTurbulenceAmplitudePolicy(
    const std::string& text, TurbulenceAmplitudePolicy* policy) {
  if (!policy) return false;
  const std::string value = Lower(text);
  if (value == "reject") *policy = TurbulenceAmplitudePolicy::Reject;
  else if (value == "limit-to-mean-field")
    *policy = TurbulenceAmplitudePolicy::LimitToMeanField;
  else return false;
  return true;
}

Status ValidateConfiguration(const Configuration& configuration) {
  if (configuration.spatial == SpatialKind::FromMeanFreePath &&
      configuration.meanFreePath == MeanFreePathKind::FromSpatial) {
    return Status::Error(StatusCode::InvalidCoefficient,
        "from-mfp spatial and from-spatial MFP providers form a conversion cycle");
  }
  if (configuration.source != SourceMode::Prescribed &&
      (configuration.meanFreePath == MeanFreePathKind::Tenishev2005 ||
       configuration.meanFreePath == MeanFreePathKind::Chen2024)) {
    return Status::Error(StatusCode::InvalidCoefficient,
        "Tenishev-2005 and Chen-2024 are prescribed MFP models and cannot "
        "consume self-consistent or SWMF turbulence");
  }
  if (configuration.source != SourceMode::Prescribed &&
      configuration.pitchAngle == PitchAngleKind::Configured) {
    return Status::Error(StatusCode::InvalidCoefficient,
        "a coupled coefficient source requires an explicit source-bound "
        "pitch-angle provider, not the legacy configured callback");
  }
  if (!std::isfinite(configuration.prescribedDeltaBOverB) ||
      configuration.prescribedDeltaBOverB < 0.0 ||
      !std::isfinite(configuration.constantDmumuPerS) ||
      configuration.constantDmumuPerS < 0.0 ||
      !std::isfinite(configuration.correlationLengthAt1AuM) ||
      configuration.correlationLengthAt1AuM <= 0.0) {
    return Status::Error(StatusCode::InvalidCoefficient,
        "prescribed turbulence ratio and correlation length are invalid");
  }
  const Status spectrum =
      CoefficientPhysics::ValidateSpectrum(configuration.spectrum);
  if (!spectrum.ok()) return spectrum;
  if (!std::isfinite(configuration.florinskiy.spectralIndex) ||
      configuration.florinskiy.spectralIndex <= 1.0 ||
      !std::isfinite(configuration.florinskiy.parallelCorrelationLengthM) ||
      configuration.florinskiy.parallelCorrelationLengthM <= 0.0) {
    return Status::Error(StatusCode::InvalidCoefficient,
        "Florinskiy spectrum parameters are invalid");
  }
  if (!std::isfinite(configuration.spatialQuadrature.absoluteToleranceM2PerS) ||
      configuration.spatialQuadrature.absoluteToleranceM2PerS < 0.0 ||
      !std::isfinite(configuration.spatialQuadrature.relativeTolerance) ||
      configuration.spatialQuadrature.relativeTolerance <= 0.0 ||
      configuration.spatialQuadrature.maximumRecursion < 1 ||
      configuration.spatialQuadrature.maximumRecursion > 60) {
    return Status::Error(StatusCode::InvalidCoefficient,
        "spatial-diffusion quadrature tolerances are invalid");
  }
  return Status::Ok();
}

Status ValidateMoverCompatibility(const Configuration& configuration,
                                  const std::string& moverCanonicalName) {
  const Status basic = ValidateConfiguration(configuration);
  if (!basic.ok()) return basic;
  const std::string mover = Lower(moverCanonicalName);
  if (mover != "parker" && mover != "fte-dmumu" && mover != "fte-mfp") {
    return Status::Error(StatusCode::UnsupportedConfiguration,
                         "unknown mover in coefficient compatibility check");
  }
  if (mover == "parker" &&
      configuration.spatial == SpatialKind::FromMeanFreePath &&
      configuration.invalidPolicy == InvalidPolicy::Ballistic) {
    return Status::Error(StatusCode::UnsupportedConfiguration,
        "parker with spatial-from-mfp does not support ballistic lambda; "
        "select invalid-coefficient-policy=fail or fte-mfp");
  }
  if (mover == "parker" &&
      configuration.resonanceGapPolicy ==
          CoefficientPhysics::ResonanceGapPolicy::Ballistic) {
    return Status::Error(StatusCode::UnsupportedConfiguration,
        "parker has no finite numerical operator for infinite kappa; select "
        "resonance-gap-policy=reject or use fte-mfp");
  }
  return Status::Ok();
}

std::string ConfigurationFingerprint(const Configuration& configuration) {
  // FNV-1a is a compact deterministic identity, not a security primitive.
  // Hash individual fields instead of the struct because padding bytes and
  // enum storage are implementation-defined.
  std::uint64_t hash = UINT64_C(1469598103934665603);
  const auto addBytes = [&hash](const void* bytes, std::size_t count) {
    const unsigned char* value =
        static_cast<const unsigned char*>(bytes);
    for (std::size_t i = 0; i < count; ++i) {
      hash ^= value[i];
      hash *= UINT64_C(1099511628211);
    }
  };
  const std::string version = "srcsep-coefficients-v2";
  addBytes(version.data(), version.size());
  const int enums[] = {
      static_cast<int>(configuration.source),
      static_cast<int>(configuration.spatial),
      static_cast<int>(configuration.pitchAngle),
      static_cast<int>(configuration.meanFreePath),
      static_cast<int>(configuration.invalidPolicy),
      static_cast<int>(configuration.resonanceGapPolicy),
      static_cast<int>(configuration.amplitudePolicy)};
  addBytes(enums, sizeof(enums));
  const double values[] = {
      configuration.prescribedDeltaBOverB,
      configuration.constantDmumuPerS,
      configuration.correlationLengthAt1AuM,
      configuration.spectrum.referenceRadiusM,
      configuration.spectrum.kMinAtReferencePerM,
      configuration.spectrum.kMaxAtReferencePerM,
      configuration.spectrum.kMinRadialExponent,
      configuration.spectrum.kMaxRadialExponent,
      configuration.spectrum.spectralIndex,
      configuration.florinskiy.spectralIndex,
      configuration.florinskiy.parallelCorrelationLengthM,
      configuration.spatialQuadrature.absoluteToleranceM2PerS,
      configuration.spatialQuadrature.relativeTolerance};
  addBytes(values, sizeof(values));
  addBytes(&configuration.spatialQuadrature.maximumRecursion,
           sizeof(configuration.spatialQuadrature.maximumRecursion));
  std::ostringstream output;
  output << std::hex << std::setw(16) << std::setfill('0') << hash;
  return output.str();
}

Status ValidateSourceAgainstBackground(SourceMode source,
                                       Background::Provider provider,
                                       Background::Ownership ownership) {
  if (source == SourceMode::Swmf &&
      (provider != Background::Provider::Swmf ||
       ownership != Background::Ownership::ImportedReadOnly)) {
    return Status::Error(StatusCode::InvalidCoefficient,
        "SWMF coefficient source requires a read-only SWMF background snapshot");
  }
  if (source == SourceMode::SelfConsistent &&
      provider == Background::Provider::Swmf &&
      ownership == Background::Ownership::ImportedReadOnly) {
    return Status::Error(StatusCode::InvalidCoefficient,
        "self-consistent coefficients cannot mutate or relabel read-only SWMF turbulence");
  }
  return Status::Ok();
}

ScalarResult KappaFromMeanFreePath(double lambdaM, double speedMPerS) {
  ScalarResult result;
  if (!std::isfinite(lambdaM) || lambdaM <= 0.0 ||
      !std::isfinite(speedMPerS) || speedMPerS <= 0.0) {
    result.status = Status::Error(StatusCode::InvalidArgument,
                                  "kappa=v*lambda/3 requires finite positive inputs");
    return result;
  }
  result.value = speedMPerS * lambdaM / 3.0;
  result.status = std::isfinite(result.value) ? Status::Ok()
      : Status::Error(StatusCode::InvalidCoefficient,
                      "mean-free-path to kappa conversion overflowed");
  return result;
}

ScalarResult MeanFreePathFromKappa(double kappaM2PerS, double speedMPerS) {
  ScalarResult result;
  if (!std::isfinite(kappaM2PerS) || kappaM2PerS < 0.0 ||
      !std::isfinite(speedMPerS) || speedMPerS <= 0.0) {
    result.status = Status::Error(StatusCode::InvalidArgument,
                                  "lambda=3*kappa/v requires kappa>=0 and v>0");
    return result;
  }
  result.value = 3.0 * kappaM2PerS / speedMPerS;
  result.status = result.value > 0.0 ? Status::Ok()
      : Status::Error(StatusCode::InvalidCoefficient,
                      "zero kappa has no finite event-driven mean free path");
  return result;
}

ScalarResult IsotropicDmumuFromMeanFreePath(double lambdaM,
                                            double speedMPerS, double mu) {
  ScalarResult result;
  if (!std::isfinite(mu) || std::fabs(mu) > 1.0 ||
      !std::isfinite(lambdaM) || lambdaM <= 0.0 ||
      !std::isfinite(speedMPerS) || speedMPerS < 0.0) {
    result.status = Status::Error(StatusCode::InvalidArgument,
                                  "isotropic Dmumu closure received invalid input");
    return result;
  }
  result.value = 0.5 * speedMPerS / lambdaM * (1.0 - mu * mu);
  result.status = Status::Ok();
  return result;
}

ScalarResult MeanFreePathFromIsotropicDmumu(double dMuMuPerS,
                                            double speedMPerS, double mu) {
  ScalarResult result;
  const double shape = 1.0 - mu * mu;
  if (!std::isfinite(dMuMuPerS) || dMuMuPerS <= 0.0 ||
      !std::isfinite(speedMPerS) || speedMPerS <= 0.0 ||
      !std::isfinite(mu) || shape <= 0.0) {
    result.status = Status::Error(StatusCode::InvalidArgument,
        "inverse isotropic Dmumu closure requires D>0, v>0, and |mu|<1");
    return result;
  }
  result.value = 0.5 * speedMPerS * shape / dMuMuPerS;
  result.status = std::isfinite(result.value) ? Status::Ok()
      : Status::Error(StatusCode::InvalidCoefficient,
                      "Dmumu to mean-free-path conversion overflowed");
  return result;
}

Configuration& ActiveConfiguration() { return gActiveConfiguration; }
Status SetActiveConfiguration(const Configuration& configuration) {
  const Status status = ValidateConfiguration(configuration);
  if (status.ok()) gActiveConfiguration = configuration;
  return status;
}

}  // namespace Coefficient
}  // namespace Transport
}  // namespace SEP
