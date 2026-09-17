#include "runtime_adapters.h"

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

}  // namespace RuntimeModel
}  // namespace SEP3D
