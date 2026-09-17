#include "background_snapshot.h"

#include <cmath>
#include <limits>

namespace SEP3D {
namespace Background {
namespace {

Core::Status Invalid(const std::string& message) {
  return Core::Status(Core::StatusCode::SnapshotUnavailable, message);
}

bool FiniteVec(const Core::Vec3& value) {
  return std::isfinite(value.x) && std::isfinite(value.y) &&
         std::isfinite(value.z);
}

bool FiniteTensor(const Core::Tensor3& value) {
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j)
      if (!std::isfinite(value(i, j))) return false;
  return true;
}

double Interpolate(double left, double right, double fraction) {
  return left + fraction * (right - left);
}

Core::Vec3 Interpolate(const Core::Vec3& left, const Core::Vec3& right,
                       double fraction) {
  return left + fraction * (right - left);
}

Core::Tensor3 Interpolate(const Core::Tensor3& left,
                          const Core::Tensor3& right, double fraction) {
  Core::Tensor3 result;
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j)
      result(i, j) = Interpolate(left(i, j), right(i, j), fraction);
  return result;
}

BackgroundSample InterpolateSample(const BackgroundSample& left,
                                   const BackgroundSample& right,
                                   double fraction,
                                   std::uint64_t generation) {
  BackgroundSample result;
  result.status = Core::Status::OK();
  result.valid = true;
  result.B = Interpolate(left.B, right.B, fraction);
  result.absB = result.B.Norm();
  result.bHat = result.B.Normalized();
  result.gradB = Interpolate(left.gradB, right.gradB, fraction);
  result.divBhat = Interpolate(left.divBhat, right.divBhat, fraction);
  result.focusingLenM = Interpolate(left.focusingLenM,
                                    right.focusingLenM, fraction);
  result.curvature = Interpolate(left.curvature, right.curvature, fraction);
  result.U = Interpolate(left.U, right.U, fraction);
  result.gradU = Interpolate(left.gradU, right.gradU, fraction);
  result.divU = Interpolate(left.divU, right.divU, fraction);
  result.fieldAlignedStrain = Interpolate(
      left.fieldAlignedStrain, right.fieldAlignedStrain, fraction);
  result.numberDensityM3 = Interpolate(left.numberDensityM3,
                                       right.numberDensityM3, fraction);
  result.temperatureK = Interpolate(left.temperatureK,
                                    right.temperatureK, fraction);
  result.pressurePa = Interpolate(left.pressurePa, right.pressurePa, fraction);
  result.alfvenSpeedMpS = Interpolate(left.alfvenSpeedMpS,
                                      right.alfvenSpeedMpS, fraction);
  result.generation = generation;
  result.configurationDigest = left.configurationDigest;
  return result;
}

Core::Status ValidateMetadata(const SnapshotMetadata& metadata) {
  if (!std::isfinite(metadata.epochS) ||
      !std::isfinite(metadata.validFromS) ||
      !std::isfinite(metadata.validUntilS) ||
      metadata.validUntilS < metadata.validFromS ||
      metadata.epochS < metadata.validFromS ||
      metadata.epochS > metadata.validUntilS ||
      metadata.generation == 0 || metadata.coordinateFrame.empty() ||
      metadata.providerIdentity.empty() ||
      metadata.configurationFingerprint.empty()) {
    return Invalid("snapshot metadata is incomplete or inconsistent");
  }
  return Core::Status::OK();
}

}  // namespace

Core::Status ValidateCompleteSample(
    const BackgroundSample& sample,
    const ProviderCapabilities& capabilities) {
  if (!sample.status.ok()) return sample.status;
  if (!sample.valid) return Invalid("background sample validity is false");
  if (!FiniteVec(sample.B) || !std::isfinite(sample.absB) ||
      !(sample.absB > 0.0) || !FiniteVec(sample.bHat) ||
      std::fabs(sample.bHat.Norm() - 1.0) > 1.0e-10 ||
      !FiniteVec(sample.U) || !std::isfinite(sample.numberDensityM3) ||
      !(sample.numberDensityM3 > 0.0) ||
      !std::isfinite(sample.temperatureK) || !(sample.temperatureK > 0.0) ||
      !std::isfinite(sample.pressurePa) || !(sample.pressurePa > 0.0) ||
      !std::isfinite(sample.alfvenSpeedMpS) ||
      !(sample.alfvenSpeedMpS > 0.0) || sample.generation == 0) {
    return Invalid("background sample has a missing/non-finite required primitive");
  }
  if (capabilities.hasAnalyticGradB && !FiniteTensor(sample.gradB))
    return Invalid("background sample magnetic gradient is non-finite");
  if (capabilities.hasAnalyticDivBhat &&
      (!std::isfinite(sample.divBhat) ||
       !std::isfinite(sample.focusingLenM))) {
    return Invalid("background sample focusing quantities are non-finite");
  }
  if (capabilities.hasAnalyticCurvature && !FiniteVec(sample.curvature))
    return Invalid("background sample curvature is non-finite");
  if (capabilities.hasAnalyticDivU &&
      (!FiniteTensor(sample.gradU) || !std::isfinite(sample.divU)))
    return Invalid("background sample velocity derivatives are non-finite");
  if (capabilities.hasFieldAlignedStrain &&
      !std::isfinite(sample.fieldAlignedStrain))
    return Invalid("background sample field-aligned strain is non-finite");
  return Core::Status::OK();
}

BackgroundSnapshot::BackgroundSnapshot(
    const SnapshotMetadata& metadata,
    const ProviderCapabilities& capabilities,
    const std::vector<Core::Vec3>& positions,
    const std::vector<BackgroundSample>& samples)
    : metadata_(metadata), capabilities_(capabilities),
      positions_(positions), samples_(samples) {}

bool BackgroundSnapshot::Covers(double timeS) const {
  return std::isfinite(timeS) && timeS >= metadata_.validFromS &&
         timeS <= metadata_.validUntilS;
}

Core::Status BackgroundSnapshotBuilder::Build(
    const BackgroundProvider& provider,
    const std::vector<Core::Vec3>& positions,
    std::shared_ptr<const BackgroundSnapshot>* output) const {
  if (output == nullptr) return Invalid("snapshot output pointer is null");
  if (positions.empty()) return Invalid("snapshot requires at least one point");
  const SnapshotMetadata* metadata = provider.PreparedMetadata();
  if (metadata == nullptr) return Invalid("background provider is not prepared");
  const Core::Status metadataStatus = ValidateMetadata(*metadata);
  if (!metadataStatus.ok()) return metadataStatus;

  std::vector<double> x(positions.size()), y(positions.size()), z(positions.size());
  for (std::size_t i = 0; i < positions.size(); ++i) {
    if (!FiniteVec(positions[i])) return Invalid("snapshot position is non-finite");
    x[i] = positions[i].x; y[i] = positions[i].y; z[i] = positions[i].z;
  }
  std::vector<BackgroundSample> samples(positions.size());
  std::vector<Core::Status> statuses(positions.size());
  const Core::Status batch = provider.EvaluateBatchDetailed(
      x.data(), y.data(), z.data(), positions.size(), samples.data(),
      statuses.data());
  if (!batch.ok()) return batch;
  const ProviderCapabilities capabilities = provider.Capabilities();
  for (std::size_t i = 0; i < samples.size(); ++i) {
    if (!statuses[i].ok()) return statuses[i];
    const Core::Status complete = ValidateCompleteSample(samples[i], capabilities);
    if (!complete.ok()) return complete;
    if (samples[i].generation != metadata->generation)
      return Invalid("sample generation differs from prepared metadata");
  }

  std::shared_ptr<const BackgroundSnapshot> candidate(
      new BackgroundSnapshot(*metadata, capabilities, positions, samples));
  *output = candidate;
  return Core::Status::OK();
}

Core::Status SnapshotBuffer::PublishCurrent(
    const std::shared_ptr<const BackgroundSnapshot>& snapshot) {
  if (!snapshot) return Invalid("current snapshot is null");
  const Core::Status metadata = ValidateMetadata(snapshot->metadata());
  if (!metadata.ok()) return metadata;
  if (current_ && snapshot->metadata().generation <=
                      current_->metadata().generation) {
    return Invalid("current snapshot generation must increase");
  }
  current_ = snapshot;
  next_.reset();
  return Core::Status::OK();
}

Core::Status SnapshotBuffer::StageNext(
    const std::shared_ptr<const BackgroundSnapshot>& snapshot) {
  if (!current_) return Invalid("stage-next requires a current snapshot");
  if (!snapshot) return Invalid("next snapshot is null");
  const SnapshotMetadata& left = current_->metadata();
  const SnapshotMetadata& right = snapshot->metadata();
  if (right.generation <= left.generation || right.epochS <= left.epochS ||
      right.provider != left.provider || right.ownership != left.ownership ||
      right.coordinateFrame != left.coordinateFrame ||
      right.configurationFingerprint != left.configurationFingerprint ||
      snapshot->positions() != current_->positions() ||
      snapshot->samples().size() != current_->samples().size()) {
    return Invalid("next snapshot identity, epoch, or grid is incompatible");
  }
  next_ = snapshot;
  return Core::Status::OK();
}

Core::Status SnapshotBuffer::SnapshotAt(
    double timeS, std::shared_ptr<const BackgroundSnapshot>* output) const {
  if (output == nullptr) return Invalid("interpolated snapshot output is null");
  if (!current_) return Invalid("no current snapshot is published");
  if (!std::isfinite(timeS)) return Invalid("requested snapshot time is non-finite");
  if (!next_) {
    if (!current_->Covers(timeS)) return Invalid("snapshot extrapolation is forbidden");
    *output = current_;
    return Core::Status::OK();
  }
  const double leftTime = current_->metadata().epochS;
  const double rightTime = next_->metadata().epochS;
  if (timeS < leftTime || timeS > rightTime)
    return Invalid("snapshot interpolation request is outside the epoch bracket");
  if (timeS == leftTime) { *output = current_; return Core::Status::OK(); }
  if (timeS == rightTime) { *output = next_; return Core::Status::OK(); }
  const double fraction = (timeS - leftTime) / (rightTime - leftTime);
  std::vector<BackgroundSample> samples;
  samples.reserve(current_->samples().size());
  for (std::size_t i = 0; i < current_->samples().size(); ++i) {
    if (current_->samples()[i].configurationDigest !=
        next_->samples()[i].configurationDigest) {
      return Invalid("cannot interpolate different provider configurations");
    }
    samples.push_back(InterpolateSample(
        current_->samples()[i], next_->samples()[i], fraction,
        next_->metadata().generation));
  }
  SnapshotMetadata metadata = current_->metadata();
  metadata.epochS = timeS;
  metadata.validFromS = leftTime;
  metadata.validUntilS = rightTime;
  metadata.generation = next_->metadata().generation;
  std::shared_ptr<const BackgroundSnapshot> candidate(new BackgroundSnapshot(
      metadata, current_->capabilities(), current_->positions(), samples));
  *output = candidate;
  return Core::Status::OK();
}

}  // namespace Background
}  // namespace SEP3D
