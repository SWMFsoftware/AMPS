#include "bg_swmf.h"
#include "background_snapshot.h"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <iomanip>
#include <limits>
#include <sstream>

namespace SEP3D {
namespace Background {
namespace {

constexpr double kMu0 = 4.0e-7 * Core::Const::kPi;

Core::Status Invalid(const std::string& message) {
  return Core::Status(Core::StatusCode::SnapshotUnavailable, message);
}

bool SamePoint(const Core::Vec3& left, const Core::Vec3& right) {
  const double scale = std::max(1.0, std::max(left.Norm(), right.Norm()));
  return (left - right).Norm() <= 1.0e-12 * scale;
}

std::uint64_t Digest64(const std::string& text) {
  std::uint64_t result = 1469598103934665603ULL;
  for (unsigned char value : text) {
    result ^= value;
    result *= 1099511628211ULL;
  }
  return result;
}

Core::Status Convert(const SwmfImport& imported, const SwmfRawSample& raw,
                     std::uint64_t digest, SwmfAwsomProvider::Record* record) {
  if (record == nullptr) return Invalid("SWMF conversion output is null");
  if (!raw.complete) return Invalid("SWMF record is missing a required field");
  if (!std::isfinite(raw.epochS) || raw.epochS != imported.epochS)
    return Invalid("SWMF record epoch differs from the import epoch");

  const double length = imported.units == SwmfUnitSystem::SI
                            ? 1.0 : Core::Const::R_sun;
  const double magnetic = imported.units == SwmfUnitSystem::SI ? 1.0 : 1.0e-9;
  const double velocity = imported.units == SwmfUnitSystem::SI ? 1.0 : 1.0e3;
  const double density = imported.units == SwmfUnitSystem::SI ? 1.0 : 1.0e6;
  const double pressure = imported.units == SwmfUnitSystem::SI ? 1.0 : 1.0e-9;

  record->positionM = length * raw.position;
  BackgroundSample sample;
  sample.B = magnetic * raw.magnetic;
  sample.absB = sample.B.Norm();
  sample.bHat = sample.B.Normalized();
  sample.U = velocity * raw.velocity;
  sample.numberDensityM3 = density * raw.numberDensity;
  sample.temperatureK = raw.temperature;
  sample.pressurePa = pressure * raw.pressure;
  sample.divU = raw.velocityDivergence;
  sample.divBhat = raw.divBhat / length;
  sample.focusingLenM = raw.focusingLength * length;
  sample.curvature = raw.curvature / length;
  sample.fieldAlignedStrain = raw.fieldAlignedStrain;
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      sample.gradB(i, j) = raw.magneticGradient(i, j) * magnetic / length;
      sample.gradU(i, j) = raw.velocityGradient(i, j) * velocity / length;
    }
  }
  sample.alfvenSpeedMpS = sample.absB /
      std::sqrt(kMu0 * Core::Const::m_p * sample.numberDensityM3);
  sample.generation = imported.generation;
  sample.configurationDigest = digest;
  sample.status = Core::Status::OK();
  sample.valid = true;
  record->sample = sample;
  return Core::Status::OK();
}

}  // namespace

Core::Status SwmfAwsomProvider::Load(const SwmfImport& candidate) {
  if (candidate.coordinateFrame.empty() ||
      candidate.coordinateFrame != candidate.expectedCoordinateFrame) {
    return Invalid("SWMF coordinate frame does not match the configured frame");
  }
  if (candidate.ownership != StorageOwnership::ImportedReadOnly) {
    return Invalid("SWMF storage must be declared imported read-only");
  }
  if (!std::isfinite(candidate.epochS) ||
      !std::isfinite(candidate.validUntilS) ||
      candidate.validUntilS < candidate.epochS || candidate.generation == 0 ||
      candidate.providerIdentity.empty() ||
      candidate.configurationFingerprint.empty() || candidate.samples.empty()) {
    return Invalid("SWMF import metadata is incomplete");
  }

  const std::uint64_t digest = Digest64(candidate.configurationFingerprint);
  std::vector<Record> records;
  records.reserve(candidate.samples.size());
  for (const SwmfRawSample& raw : candidate.samples) {
    Record record;
    const Core::Status converted = Convert(candidate, raw, digest, &record);
    if (!converted.ok()) return converted;
    // The complete-field validator is intentionally applied before commit.
    // Imported data therefore cannot become active with a NaN that later
    // masquerades as a zero-valued physical field.
    const Core::Status complete = ValidateCompleteSample(
        record.sample, Capabilities());
    if (!complete.ok()) return complete;
    records.push_back(record);
  }

  SnapshotMetadata metadata;
  metadata.provider = ProviderKind::SwmfAwsom;
  metadata.ownership = candidate.ownership;
  metadata.epochS = candidate.epochS;
  metadata.validFromS = candidate.epochS;
  metadata.validUntilS = candidate.validUntilS;
  metadata.generation = candidate.generation;
  metadata.coordinateFrame = candidate.coordinateFrame;
  metadata.providerIdentity = candidate.providerIdentity;
  metadata.configurationFingerprint = candidate.configurationFingerprint;

  // Atomic commit after conversion and complete validation of every record.
  import_ = candidate;
  metadata_ = metadata;
  records_.swap(records);
  loaded_ = true;
  prepared_ = false;
  return Core::Status::OK();
}

Core::Status SwmfAwsomProvider::Validate() const {
  return loaded_ ? Core::Status::OK()
                 : Invalid("no validated SWMF import has been loaded");
}

Core::Status SwmfAwsomProvider::Prepare(double timeS) {
  const Core::Status valid = Validate();
  if (!valid.ok()) return valid;
  if (!std::isfinite(timeS) || timeS != metadata_.epochS)
    return Invalid("SWMF preparation time must equal the imported epoch");
  prepared_ = true;
  return Core::Status::OK();
}

const SnapshotMetadata* SwmfAwsomProvider::PreparedMetadata() const {
  return prepared_ ? &metadata_ : nullptr;
}

BackgroundSample SwmfAwsomProvider::Evaluate(const Core::Vec3& positionM) const {
  BackgroundSample failure;
  if (!prepared_) {
    failure.status = Invalid("SWMF provider is not prepared");
    return failure;
  }
  for (const Record& record : records_) {
    if (SamePoint(positionM, record.positionM)) return record.sample;
  }
  failure.status = Core::Status(Core::StatusCode::NotFound,
                                "SWMF import has no record at this point");
  return failure;
}

Core::Status SwmfAwsomProvider::EvaluateBatchDetailed(
    const double* xM, const double* yM, const double* zM,
    std::size_t count, BackgroundSample* output,
    Core::Status* perSampleStatus) const {
  if ((count != 0) && (xM == nullptr || yM == nullptr || zM == nullptr ||
                       output == nullptr || perSampleStatus == nullptr)) {
    return Core::Status(Core::StatusCode::InvalidInput,
                        "SWMF batch arrays must not be null");
  }
  Core::Status aggregate = Core::Status::OK();
  for (std::size_t i = 0; i < count; ++i) {
    const BackgroundSample candidate = Evaluate({xM[i], yM[i], zM[i]});
    perSampleStatus[i] = candidate.status;
    if (candidate.status.ok() && candidate.valid) output[i] = candidate;
    else if (aggregate.ok()) aggregate = candidate.status;
  }
  return aggregate;
}

std::string SwmfAwsomProvider::ResolvedManifest() const {
  std::ostringstream out;
  out << "swmf-awsom-import-v1;units="
      << (import_.units == SwmfUnitSystem::SI ? "SI" : "AWSOM_COUPLING")
      << ";frame=" << import_.coordinateFrame
      << ";ownership=imported-read-only"
      << ";epoch_s=" << std::setprecision(17) << import_.epochS
      << ";generation=" << import_.generation
      << ";records=" << records_.size();
  return out.str();
}

ProviderCapabilities SwmfAwsomProvider::Capabilities() const {
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
