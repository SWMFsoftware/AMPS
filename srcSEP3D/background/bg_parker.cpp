#include "bg_parker.h"

#include "sep_background_snapshot.h"

#include <cmath>
#include <iomanip>
#include <limits>
#include <sstream>

namespace SEP3D {
namespace Background {
namespace {

constexpr double kMu0 = 4.0e-7 * Core::Const::kPi;

Core::Status Invalid(const std::string& message) {
  return Core::Status(Core::StatusCode::InvalidInput, message);
}

std::uint64_t Digest64(const std::string& text) {
  // FNV-1a is used only as a compact in-cell tag; the complete versioned
  // fingerprint string remains in SnapshotMetadata and the run manifest.
  std::uint64_t result = 1469598103934665603ULL;
  for (unsigned char value : text) {
    result ^= value;
    result *= 1099511628211ULL;
  }
  return result;
}

Core::Vec3 AxisCrossBasisColumn(const Core::Vec3& axis, int column) {
  if (column == 0) return {0.0, axis.z, -axis.y};
  if (column == 1) return {-axis.z, 0.0, axis.x};
  return {axis.y, -axis.x, 0.0};
}

double Component(const Core::Vec3& value, int index) {
  return index == 0 ? value.x : (index == 1 ? value.y : value.z);
}

}  // namespace

AnalyticParkerProvider::AnalyticParkerProvider(
    const ParkerConfiguration& configuration)
    : configuration_(configuration) {
  metadata_.provider = ProviderKind::AnalyticParker;
  metadata_.ownership = StorageOwnership::ModelOwned;
  metadata_.coordinateFrame = configuration_.coordinateFrame;
  metadata_.providerIdentity = CanonicalName();
  metadata_.configurationFingerprint =
      SEP::Background::FingerprintConfiguration(ResolvedManifest());
  configurationDigest_ = Digest64(metadata_.configurationFingerprint);
}

Core::Status AnalyticParkerProvider::Validate() const {
  const double values[] = {
      configuration_.sourceRadiusM, configuration_.referenceRadiusM,
      configuration_.radialFieldAtReferenceT,
      configuration_.numberDensityAtReferenceM3, configuration_.temperatureK,
      configuration_.solarWindSpeedMPerS,
      configuration_.solarRotationRateRadPerS,
      configuration_.validityCadenceS};
  for (double value : values) {
    if (!std::isfinite(value)) return Invalid("Parker configuration contains a non-finite value");
  }
  if (configuration_.sourceRadiusM <= 0.0 ||
      configuration_.referenceRadiusM <= configuration_.sourceRadiusM ||
      configuration_.radialFieldAtReferenceT <= 0.0 ||
      configuration_.numberDensityAtReferenceM3 <= 0.0 ||
      configuration_.temperatureK <= 0.0 ||
      configuration_.solarWindSpeedMPerS <= 0.0 ||
      configuration_.validityCadenceS <= 0.0 ||
      (configuration_.magneticPolarity != 1 &&
       configuration_.magneticPolarity != -1) ||
      configuration_.coordinateFrame.empty()) {
    return Invalid("Parker configuration is outside its physical range");
  }
  const double axisNorm = configuration_.rotationAxis.Norm();
  if (!std::isfinite(axisNorm) || axisNorm <= 0.0) {
    return Invalid("Parker rotation axis must be finite and non-zero");
  }
  return Core::Status::OK();
}

Core::Status AnalyticParkerProvider::Prepare(double timeS) {
  const Core::Status valid = Validate();
  if (!valid.ok()) return valid;
  if (!std::isfinite(timeS)) return Invalid("Parker preparation time is not finite");
  SnapshotMetadata candidate = metadata_;
  candidate.epochS = timeS;
  candidate.validFromS = timeS;
  candidate.validUntilS = timeS + configuration_.validityCadenceS;
  candidate.generation = metadata_.generation + 1;
  metadata_ = candidate;
  prepared_ = true;
  return Core::Status::OK();
}

const SnapshotMetadata* AnalyticParkerProvider::PreparedMetadata() const {
  return prepared_ ? &metadata_ : nullptr;
}

BackgroundSample AnalyticParkerProvider::Evaluate(
    const Core::Vec3& positionM) const {
  BackgroundSample sample;
  if (!prepared_) {
    sample.status = Core::Status(Core::StatusCode::SnapshotUnavailable,
                                 "Parker provider is not prepared");
    return sample;
  }
  const double radius = positionM.Norm();
  if (!std::isfinite(radius) || radius < configuration_.sourceRadiusM) {
    sample.status = Core::Status(Core::StatusCode::BackgroundInvalid,
                                 "position is inside the Parker source surface");
    return sample;
  }

  const Core::Vec3 rHat = positionM / radius;
  const Core::Vec3 axis = configuration_.rotationAxis.Normalized();
  const Core::Vec3 axisCrossX = axis.Cross(positionM);
  const double polarity = static_cast<double>(configuration_.magneticPolarity);
  const double coefficient = polarity * configuration_.radialFieldAtReferenceT *
      configuration_.referenceRadiusM * configuration_.referenceRadiusM;
  const double winding = configuration_.solarRotationRateRadPerS /
                         configuration_.solarWindSpeedMPerS;
  const double inverseR3 = 1.0 / (radius * radius * radius);
  const double spiralFactor = winding * (radius - configuration_.sourceRadiusM);
  sample.B = coefficient * inverseR3 *
      (positionM - spiralFactor * axisCrossX);
  sample.absB = sample.B.Norm();
  if (!(sample.absB > 0.0) || !std::isfinite(sample.absB)) {
    sample.status = Core::Status(Core::StatusCode::BackgroundInvalid,
                                 "Parker field magnitude is invalid");
    return sample;
  }
  sample.bHat = sample.B / sample.absB;

  // Closed Cartesian derivative.  Writing the spiral as
  // B=C[x/r^3-k(r-r0)(a x x)/r^3] avoids singular spherical basis vectors at
  // the rotation axis and makes div(B)=0 to round-off.
  const double inverseR5 = inverseR3 / (radius * radius);
  const double f = (radius - configuration_.sourceRadiusM) * inverseR3;
  const double dfScale = -2.0 / std::pow(radius, 4) +
                         3.0 * configuration_.sourceRadiusM * inverseR5;
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      const double xi = Component(positionM, i);
      const double xj = Component(positionM, j);
      const double radialDerivative = (i == j ? inverseR3 : 0.0) -
                                      3.0 * xi * xj * inverseR5;
      const Core::Vec3 crossColumn = AxisCrossBasisColumn(axis, j);
      const double spiralDerivative = dfScale * xj *
          Component(axisCrossX, i) + f * Component(crossColumn, i);
      sample.gradB(i, j) = coefficient *
          (radialDerivative - winding * spiralDerivative);
    }
  }

  Core::Vec3 gradAbsB;
  for (int j = 0; j < 3; ++j) {
    double value = 0.0;
    for (int i = 0; i < 3; ++i)
      value += Component(sample.B, i) * sample.gradB(i, j);
    if (j == 0) gradAbsB.x = value / sample.absB;
    else if (j == 1) gradAbsB.y = value / sample.absB;
    else gradAbsB.z = value / sample.absB;
  }
  sample.divBhat = sample.gradB.Trace() / sample.absB -
      sample.B.Dot(gradAbsB) / (sample.absB * sample.absB);
  const double dLnBds = sample.bHat.Dot(gradAbsB) / sample.absB;
  sample.focusingLenM = std::fabs(dLnBds) > 0.0
      ? -1.0 / dLnBds : std::numeric_limits<double>::infinity();
  Core::Tensor3 gradBhat;
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j)
      gradBhat(i, j) = sample.gradB(i, j) / sample.absB -
          Component(sample.B, i) * Component(gradAbsB, j) /
              (sample.absB * sample.absB);
  sample.curvature = gradBhat.Apply(sample.bHat);

  sample.U = configuration_.solarWindSpeedMPerS * rHat;
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j)
      sample.gradU(i, j) = configuration_.solarWindSpeedMPerS / radius *
          ((i == j ? 1.0 : 0.0) -
           Component(rHat, i) * Component(rHat, j));
  sample.divU = sample.gradU.Trace();
  sample.fieldAlignedStrain =
      sample.gradU.DoubleContract(sample.bHat, sample.bHat);
  const double densityScale = configuration_.referenceRadiusM / radius;
  sample.numberDensityM3 = configuration_.numberDensityAtReferenceM3 *
                           densityScale * densityScale;
  sample.temperatureK = configuration_.temperatureK;
  sample.pressurePa = sample.numberDensityM3 * Core::Const::k_B *
                      sample.temperatureK;
  sample.alfvenSpeedMpS = sample.absB /
      std::sqrt(kMu0 * Core::Const::m_p * sample.numberDensityM3);
  sample.generation = metadata_.generation;
  sample.configurationDigest = configurationDigest_;
  sample.valid = true;
  sample.status = Core::Status::OK();
  return sample;
}

std::string AnalyticParkerProvider::ResolvedManifest() const {
  std::ostringstream out;
  out << std::setprecision(17) << std::scientific
      << "parker-provider-v1"
      << ";source_m=" << configuration_.sourceRadiusM
      << ";reference_m=" << configuration_.referenceRadiusM
      << ";Br_ref_T=" << configuration_.radialFieldAtReferenceT
      << ";n_ref_m-3=" << configuration_.numberDensityAtReferenceM3
      << ";temperature_K=" << configuration_.temperatureK
      << ";wind_m_s=" << configuration_.solarWindSpeedMPerS
      << ";omega_rad_s=" << configuration_.solarRotationRateRadPerS
      << ";axis=" << configuration_.rotationAxis.x << ','
      << configuration_.rotationAxis.y << ',' << configuration_.rotationAxis.z
      << ";polarity=" << configuration_.magneticPolarity
      << ";cadence_s=" << configuration_.validityCadenceS
      << ";frame=" << configuration_.coordinateFrame;
  return out.str();
}

ProviderCapabilities AnalyticParkerProvider::Capabilities() const {
  ProviderCapabilities capabilities;
  capabilities.hasAnalyticGradB = true;
  capabilities.hasAnalyticDivBhat = true;
  capabilities.hasAnalyticCurvature = true;
  capabilities.hasAnalyticDivU = true;
  capabilities.hasFieldAlignedStrain = true;
  capabilities.hasPlasmaState = true;
  capabilities.supportsBatchEval = true;
  return capabilities;
}

}  // namespace Background
}  // namespace SEP3D
