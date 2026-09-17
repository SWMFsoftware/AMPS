#include "runtime_adapters.h"

#include "../background/background_snapshot.h"

namespace SEP3D {
namespace RuntimeModel {

namespace {

Core::Status MissingRuntime() {
  return Core::Status(Core::StatusCode::InvalidInput,
                      "runtime adapter received a null Runtime");
}

SnapshotDescriptor Candidate(
    const Runtime& runtime, BackgroundAuthority authority, double epochS,
    double validUntilS, std::uint64_t generation, bool complete,
    const char* providerIdentity) {
  SnapshotDescriptor result;
  result.authority = authority;
  result.epochS = epochS;
  result.validFromS = epochS;
  result.validUntilS = validUntilS;
  result.generation = generation;
  result.complete = complete;
  result.coordinateFrame = "HCI-like-inertial";
  result.providerIdentity = providerIdentity;
  if (runtime.configuration()) {
    result.configurationFingerprint =
        runtime.configuration()->physics_fingerprint();
  }
  return result;
}

Core::Status CandidateFromSnapshot(
    const Runtime& runtime,
    const Background::BackgroundSnapshot& snapshot,
    BackgroundAuthority expectedAuthority,
    SnapshotDescriptor* result) {
  if (result == nullptr) {
    return Core::Status(Core::StatusCode::InvalidInput,
                        "snapshot descriptor output is null");
  }
  const Background::SnapshotMetadata& metadata = snapshot.metadata();
  const bool authorityMatches =
      (expectedAuthority == BackgroundAuthority::AnalyticParker &&
       metadata.provider == Background::ProviderKind::AnalyticParker) ||
      (expectedAuthority == BackgroundAuthority::Swmf &&
       metadata.provider == Background::ProviderKind::SwmfAwsom);
  if (!authorityMatches) {
    return Core::Status(Core::StatusCode::ConfigurationConflict,
                        "snapshot provider does not match the Runtime adapter");
  }
  if (snapshot.samples().empty()) {
    return Core::Status(Core::StatusCode::SnapshotUnavailable,
                        "cannot publish an empty physical snapshot");
  }
  *result = Candidate(runtime, expectedAuthority, metadata.epochS,
                      metadata.validUntilS, metadata.generation, true,
                      metadata.providerIdentity.c_str());
  result->validFromS = metadata.validFromS;
  result->coordinateFrame = metadata.coordinateFrame;
  return Core::Status::OK();
}

}  // namespace

Core::Status StandaloneAdapter::Initialize(Runtime* runtime) const {
  if (runtime == nullptr) return MissingRuntime();
  return runtime->BeginBackgroundAcquisition(AdapterKind::Standalone);
}

Core::Status StandaloneAdapter::PublishFrozenParker(
    Runtime* runtime, double epochS, double validUntilS,
    std::uint64_t generation) const {
  if (runtime == nullptr) return MissingRuntime();
  return runtime->PublishSnapshot(Candidate(
      *runtime, BackgroundAuthority::AnalyticParker, epochS, validUntilS,
      generation, true, "analytic-parker-r2"));
}

Core::Status StandaloneAdapter::PublishSnapshot(
    Runtime* runtime,
    const Background::BackgroundSnapshot& snapshot) const {
  if (runtime == nullptr) return MissingRuntime();
  SnapshotDescriptor candidate;
  const Core::Status converted = CandidateFromSnapshot(
      *runtime, snapshot, BackgroundAuthority::AnalyticParker, &candidate);
  return converted.ok() ? runtime->PublishSnapshot(candidate) : converted;
}

Core::Status SwmfAdapter::Initialize(Runtime* runtime) const {
  if (runtime == nullptr) return MissingRuntime();
  return runtime->BeginBackgroundAcquisition(AdapterKind::Swmf);
}

Core::Status SwmfAdapter::PublishImported(
    Runtime* runtime, double epochS, double validUntilS,
    std::uint64_t generation, bool complete) const {
  if (runtime == nullptr) return MissingRuntime();
  return runtime->PublishSnapshot(Candidate(
      *runtime, BackgroundAuthority::Swmf, epochS, validUntilS, generation,
      complete, "swmf-import-r2"));
}

Core::Status SwmfAdapter::PublishSnapshot(
    Runtime* runtime,
    const Background::BackgroundSnapshot& snapshot) const {
  if (runtime == nullptr) return MissingRuntime();
  SnapshotDescriptor candidate;
  const Core::Status converted = CandidateFromSnapshot(
      *runtime, snapshot, BackgroundAuthority::Swmf, &candidate);
  return converted.ok() ? runtime->PublishSnapshot(candidate) : converted;
}

}  // namespace RuntimeModel
}  // namespace SEP3D
