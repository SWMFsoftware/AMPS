// ============================================================================
// Phase R2 lifecycle tests (LIFE3D01-LIFE3D04)
//
// These tests link only srcSEP3D's AMPS-independent Runtime and sep_common.
// They are the executable specification of the state machine: every operation
// is tried from every lifecycle state, rejected calls are checked for zero
// mutation, and both analytic and coupled adapters publish through the same
// Runtime boundary.
// ============================================================================

#include "sep3d_test_registry.h"
#include "runtime_adapters.h"

#include <array>
#include <memory>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

namespace {

namespace RM = SEP3D::RuntimeModel;
using Result = SEP3D::Testing::Result;
using TestStatus = SEP3D::Testing::Status;

Result Pass(const std::string& message) {
  Result result;
  result.status = TestStatus::Pass;
  result.message = message;
  return result;
}

Result Fail(const std::string& message) {
  Result result;
  result.status = TestStatus::Fail;
  result.message = message;
  return result;
}

bool MakeConfiguration(
    RM::BackgroundAuthority background,
    std::shared_ptr<const RM::RunConfiguration3D>* configuration,
    std::uint64_t outputCadence = 2) {
  RM::RunConfiguration3DOptions options;
  options.background = background;
  options.turbulence = background == RM::BackgroundAuthority::Swmf
                           ? RM::TurbulenceAuthority::Swmf
                           : RM::TurbulenceAuthority::Prescribed;
  options.outputCadenceSteps = outputCadence;
  options.storeMagneticGradient = true;
  options.samplingBytesPerCell = 24;
  return RM::RunConfiguration3D::Create(options, configuration).ok();
}

RM::SnapshotDescriptor Snapshot(
    const RM::Runtime& runtime, std::uint64_t generation) {
  RM::SnapshotDescriptor candidate;
  candidate.authority = runtime.configuration()->options().background;
  candidate.epochS = 0.0;
  candidate.validFromS = 0.0;
  candidate.validUntilS = 10.0;
  candidate.generation = generation;
  candidate.complete = true;
  candidate.coordinateFrame = "HCI-like-inertial";
  candidate.providerIdentity = "lifecycle-fixture";
  candidate.configurationFingerprint =
      runtime.configuration()->physics_fingerprint();
  return candidate;
}

bool ConfigureAndBind(RM::Runtime* runtime,
                      const std::shared_ptr<const RM::RunConfiguration3D>& cfg) {
  if (!runtime->Configure(cfg).ok()) return false;
  RM::MeshBinding binding;
  binding.layout = cfg->storage_layout();
  return runtime->BindMesh(binding).ok();
}

bool DriveTo(RM::LifecycleState target, RM::Runtime* runtime,
             const std::shared_ptr<const RM::RunConfiguration3D>& cfg) {
  if (target == RM::LifecycleState::Created) return true;
  if (!runtime->Configure(cfg).ok()) return false;
  if (target == RM::LifecycleState::Configured) return true;

  RM::MeshBinding binding;
  binding.layout = cfg->storage_layout();
  if (!runtime->BindMesh(binding).ok()) return false;
  if (target == RM::LifecycleState::MeshReady) return true;

  RM::StandaloneAdapter adapter;
  if (!adapter.Initialize(runtime).ok()) return false;
  if (target == RM::LifecycleState::WaitingForSnapshot) return true;
  if (!adapter.PublishFrozenParker(runtime, 0.0, 10.0, 1).ok()) return false;
  if (target == RM::LifecycleState::SnapshotReady) return true;

  if (target == RM::LifecycleState::Running) {
    return runtime->BeginStep(1.0).ok();
  }
  if (target == RM::LifecycleState::Checkpointing) {
    return runtime->BeginCheckpoint().ok();
  }
  if (target == RM::LifecycleState::Finalized) {
    return runtime->Finalize().ok();
  }
  return false;
}

Result RunLIFE3D01() {
  std::shared_ptr<const RM::RunConfiguration3D> configuration;
  if (!MakeConfiguration(RM::BackgroundAuthority::AnalyticParker,
                         &configuration)) {
    return Fail("could not construct the baseline immutable configuration");
  }

  RM::Runtime runtime;
  if (!ConfigureAndBind(&runtime, configuration)) {
    return Fail("Created -> Configured -> MeshReady failed");
  }
  RM::StandaloneAdapter adapter;
  if (!adapter.Initialize(&runtime).ok() ||
      !adapter.PublishFrozenParker(&runtime, 0.0, 10.0, 1).ok()) {
    return Fail("MeshReady -> WaitingForSnapshot -> SnapshotReady failed");
  }
  if (!runtime.BeginStep(1.0).ok() || !runtime.CompleteStep().ok() ||
      !runtime.BeginStep(2.0).ok() || !runtime.CompleteStep().ok()) {
    return Fail("SnapshotReady -> Running -> SnapshotReady failed");
  }
  if (runtime.counters().completedSteps != 2 ||
      runtime.counters().outputSequence != 1 ||
      runtime.counters().stepsSinceOutput != 0) {
    return Fail("Runtime-owned cadence counters do not match two completed steps");
  }
  if (!runtime.BeginCheckpoint().ok() ||
      !runtime.CompleteCheckpoint().ok() ||
      runtime.counters().checkpointSequence != 1) {
    return Fail("checkpoint lifecycle or restartable sequence counter failed");
  }
  if (!runtime.Finalize().ok() ||
      runtime.state() != RM::LifecycleState::Finalized) {
    return Fail("SnapshotReady -> Finalized failed");
  }
  return Pass("canonical lifecycle reached Finalized with runtime-owned cadence and checkpoint counters");
}

enum class Operation {
  Configure,
  BindMesh,
  BeginAcquisition,
  Publish,
  BeginStep,
  CompleteStep,
  BeginCheckpoint,
  CompleteCheckpoint,
  RestoreCounters,
  Finalize
};

const char* OperationName(Operation operation) {
  switch (operation) {
    case Operation::Configure: return "Configure";
    case Operation::BindMesh: return "BindMesh";
    case Operation::BeginAcquisition: return "BeginAcquisition";
    case Operation::Publish: return "Publish";
    case Operation::BeginStep: return "BeginStep";
    case Operation::CompleteStep: return "CompleteStep";
    case Operation::BeginCheckpoint: return "BeginCheckpoint";
    case Operation::CompleteCheckpoint: return "CompleteCheckpoint";
    case Operation::RestoreCounters: return "RestoreCounters";
    case Operation::Finalize: return "Finalize";
  }
  return "Unknown";
}

bool Allowed(RM::LifecycleState state, Operation operation) {
  switch (state) {
    case RM::LifecycleState::Created:
      return operation == Operation::Configure;
    case RM::LifecycleState::Configured:
      return operation == Operation::BindMesh ||
             operation == Operation::RestoreCounters ||
             operation == Operation::Finalize;
    case RM::LifecycleState::MeshReady:
      return operation == Operation::BeginAcquisition ||
             operation == Operation::Finalize;
    case RM::LifecycleState::WaitingForSnapshot:
      return operation == Operation::Publish ||
             operation == Operation::Finalize;
    case RM::LifecycleState::SnapshotReady:
      return operation == Operation::Publish ||
             operation == Operation::BeginStep ||
             operation == Operation::BeginCheckpoint ||
             operation == Operation::Finalize;
    case RM::LifecycleState::Running:
      return operation == Operation::CompleteStep;
    case RM::LifecycleState::Checkpointing:
      return operation == Operation::CompleteCheckpoint;
    case RM::LifecycleState::Finalized:
      return false;
  }
  return false;
}

RM::LifecycleState SuccessfulDestination(Operation operation) {
  switch (operation) {
    case Operation::Configure: return RM::LifecycleState::Configured;
    case Operation::BindMesh: return RM::LifecycleState::MeshReady;
    case Operation::BeginAcquisition:
      return RM::LifecycleState::WaitingForSnapshot;
    case Operation::Publish: return RM::LifecycleState::SnapshotReady;
    case Operation::BeginStep: return RM::LifecycleState::Running;
    case Operation::CompleteStep: return RM::LifecycleState::SnapshotReady;
    case Operation::BeginCheckpoint:
      return RM::LifecycleState::Checkpointing;
    case Operation::CompleteCheckpoint:
      return RM::LifecycleState::SnapshotReady;
    case Operation::RestoreCounters: return RM::LifecycleState::Configured;
    case Operation::Finalize: return RM::LifecycleState::Finalized;
  }
  return RM::LifecycleState::Finalized;
}

SEP3D::Core::Status Apply(
    Operation operation, RM::Runtime* runtime,
    const std::shared_ptr<const RM::RunConfiguration3D>& configuration) {
  switch (operation) {
    case Operation::Configure:
      return runtime->Configure(configuration);
    case Operation::BindMesh: {
      RM::MeshBinding binding;
      binding.layout = configuration->storage_layout();
      return runtime->BindMesh(binding);
    }
    case Operation::BeginAcquisition: {
      RM::StandaloneAdapter adapter;
      return adapter.Initialize(runtime);
    }
    case Operation::Publish: {
      const std::uint64_t generation = runtime->has_snapshot()
                                           ? runtime->active_snapshot()->generation + 1
                                           : 1;
      // Created has no installed configuration, but PublishSnapshot must still
      // be exercised there to prove the ordering guard fires before candidate
      // validation.  Construct an otherwise valid descriptor from the host's
      // immutable configuration instead of dereferencing Runtime state.
      RM::SnapshotDescriptor candidate;
      if (runtime->configuration()) {
        candidate = Snapshot(*runtime, generation);
      } else {
        candidate.authority = configuration->options().background;
        candidate.epochS = 0.0;
        candidate.validFromS = 0.0;
        candidate.validUntilS = 10.0;
        candidate.generation = generation;
        candidate.complete = true;
        candidate.coordinateFrame = "HCI-like-inertial";
        candidate.providerIdentity = "lifecycle-fixture";
        candidate.configurationFingerprint =
            configuration->physics_fingerprint();
      }
      return runtime->PublishSnapshot(candidate);
    }
    case Operation::BeginStep:
      return runtime->BeginStep(1.0);
    case Operation::CompleteStep:
      return runtime->CompleteStep();
    case Operation::BeginCheckpoint:
      return runtime->BeginCheckpoint();
    case Operation::CompleteCheckpoint:
      return runtime->CompleteCheckpoint();
    case Operation::RestoreCounters: {
      RM::RuntimeCounters counters;
      counters.completedSteps = 4;
      counters.stepsSinceOutput = 1;
      counters.outputSequence = 2;
      counters.checkpointSequence = 3;
      return runtime->RestoreCounters(counters);
    }
    case Operation::Finalize:
      return runtime->Finalize();
  }
  return SEP3D::Core::Status::Error("unknown lifecycle operation");
}

Result RunLIFE3D02() {
  const std::array<RM::LifecycleState, 8> states = {{
      RM::LifecycleState::Created,
      RM::LifecycleState::Configured,
      RM::LifecycleState::MeshReady,
      RM::LifecycleState::WaitingForSnapshot,
      RM::LifecycleState::SnapshotReady,
      RM::LifecycleState::Running,
      RM::LifecycleState::Checkpointing,
      RM::LifecycleState::Finalized}};
  const std::array<Operation, 10> operations = {{
      Operation::Configure, Operation::BindMesh,
      Operation::BeginAcquisition, Operation::Publish,
      Operation::BeginStep, Operation::CompleteStep,
      Operation::BeginCheckpoint, Operation::CompleteCheckpoint,
      Operation::RestoreCounters, Operation::Finalize}};

  std::shared_ptr<const RM::RunConfiguration3D> configuration;
  if (!MakeConfiguration(RM::BackgroundAuthority::AnalyticParker,
                         &configuration)) {
    return Fail("could not construct lifecycle matrix configuration");
  }

  std::size_t checked = 0;
  for (RM::LifecycleState state : states) {
    for (Operation operation : operations) {
      RM::Runtime runtime;
      if (!DriveTo(state, &runtime, configuration)) {
        return Fail(std::string("fixture could not reach ") + RM::Name(state));
      }
      const RM::LifecycleState before = runtime.state();
      const RM::RuntimeCounters countersBefore = runtime.counters();
      const RM::RunConfiguration3D* configurationBefore =
          runtime.configuration().get();
      const bool snapshotBefore = runtime.has_snapshot();
      const std::uint64_t generationBefore = snapshotBefore
                                                  ? runtime.active_snapshot()->generation
                                                  : 0;
      const std::uint64_t pinnedBefore = runtime.pinned_snapshot_generation();
      const SEP3D::Core::Status status =
          Apply(operation, &runtime, configuration);
      const bool expected = Allowed(state, operation);
      if (status.ok() != expected) {
        std::ostringstream error;
        error << RM::Name(state) << " + " << OperationName(operation)
              << " expected " << (expected ? "success" : "rejection")
              << " but returned code " << static_cast<int>(status.code)
              << ": " << status.message;
        return Fail(error.str());
      }
      if (expected && runtime.state() != SuccessfulDestination(operation)) {
        return Fail(std::string("legal ") + OperationName(operation) +
                    " ended in " + RM::Name(runtime.state()) +
                    " instead of " +
                    RM::Name(SuccessfulDestination(operation)));
      }
      if (!expected) {
        const std::uint64_t generationAfter = runtime.has_snapshot()
                                                  ? runtime.active_snapshot()->generation
                                                  : 0;
        if (runtime.state() != before ||
            runtime.counters().completedSteps != countersBefore.completedSteps ||
            runtime.counters().stepsSinceOutput != countersBefore.stepsSinceOutput ||
            runtime.counters().outputSequence != countersBefore.outputSequence ||
            runtime.counters().checkpointSequence != countersBefore.checkpointSequence ||
            runtime.configuration().get() != configurationBefore ||
            runtime.has_snapshot() != snapshotBefore ||
            generationAfter != generationBefore ||
            runtime.pinned_snapshot_generation() != pinnedBefore) {
          return Fail(std::string("illegal ") + OperationName(operation) +
                      " mutated Runtime from " + RM::Name(state));
        }
      }
      ++checked;
    }
  }

  Result result = Pass("all legal and illegal lifecycle operation/state pairs obey the transition table without partial mutation");
  result.metrics.push_back({"transition_pairs", static_cast<double>(checked),
                            80.0, "==", "pairs"});
  return result;
}

Result RunLIFE3D03() {
  // Validate authority conflicts and reserved physics before inspecting a
  // valid layout.  Create() must clear its output first, so a caller cannot
  // accidentally retain and run an older configuration after a failed reload.
  RM::RunConfiguration3DOptions conflict;
  conflict.turbulence = RM::TurbulenceAuthority::Swmf;
  std::shared_ptr<const RM::RunConfiguration3D> rejectedConfiguration;
  SEP3D::Core::Status configurationStatus =
      RM::RunConfiguration3D::Create(conflict, &rejectedConfiguration);
  if (configurationStatus.code !=
          SEP3D::Core::StatusCode::ConfigurationConflict ||
      rejectedConfiguration) {
    return Fail("conflicting background/turbulence authority was not rejected atomically");
  }
  RM::RunConfiguration3DOptions reserved;
  reserved.enableDrifts = true;
  configurationStatus =
      RM::RunConfiguration3D::Create(reserved, &rejectedConfiguration);
  if (configurationStatus.code != SEP3D::Core::StatusCode::ReservedFeature ||
      rejectedConfiguration) {
    return Fail("reserved drift physics did not fail during configuration");
  }

  RM::RunConfiguration3DOptions base;
  base.storeMagneticGradient = true;
  base.storeVelocityGradient = true;
  base.samplingBytesPerCell = 40;
  std::shared_ptr<const RM::RunConfiguration3D> first;
  if (!RM::RunConfiguration3D::Create(base, &first).ok()) {
    return Fail("could not create pre-mesh layout configuration");
  }

  RM::RunConfiguration3DOptions formatting = base;
  formatting.outputDirectory = "another-directory";
  formatting.outputPrefix = "different-prefix";
  formatting.outputCadenceSteps = 17;
  std::shared_ptr<const RM::RunConfiguration3D> second;
  if (!RM::RunConfiguration3D::Create(formatting, &second).ok() ||
      first->physics_fingerprint() != second->physics_fingerprint()) {
    return Fail("output-only changes altered the physics fingerprint");
  }

  RM::RunConfiguration3DOptions physics = base;
  physics.requestedTimeStepS = 2.0;
  std::shared_ptr<const RM::RunConfiguration3D> third;
  if (!RM::RunConfiguration3D::Create(physics, &third).ok() ||
      first->physics_fingerprint() == third->physics_fingerprint()) {
    return Fail("physics-relevant timestep change did not alter the fingerprint");
  }

  const RM::StorageLayout layout = first->storage_layout();
  if (layout.magneticFieldOffset != 0 ||
      layout.bulkVelocityOffset != 3 * sizeof(double) ||
      layout.numberDensityOffset != 6 * sizeof(double) ||
      layout.velocityDivergenceOffset != 7 * sizeof(double) ||
      layout.magneticGradientOffset != 8 * sizeof(double) ||
      layout.velocityGradientOffset != 17 * sizeof(double) ||
      layout.cellAssociatedBytes != 26 * sizeof(double) ||
      layout.samplingBytesPerCell != 40) {
    return Fail("canonical pre-mesh offsets or sizes are incorrect");
  }

  RM::RunConfiguration3DOptions imported;
  imported.background = RM::BackgroundAuthority::Swmf;
  imported.turbulence = RM::TurbulenceAuthority::Swmf;
  std::shared_ptr<const RM::RunConfiguration3D> importedConfiguration;
  if (!RM::RunConfiguration3D::Create(imported, &importedConfiguration).ok() ||
      importedConfiguration->storage_layout().waveEnergyOffset !=
          8 * sizeof(double) ||
      importedConfiguration->storage_layout().cellAssociatedBytes !=
          10 * sizeof(double)) {
    return Fail("SWMF turbulence layout did not reserve two wave-energy doubles");
  }

  RM::Runtime runtime;
  if (!runtime.Configure(first).ok()) return Fail("Runtime configuration failed");
  RM::MeshBinding bad;
  bad.layout = layout;
  ++bad.layout.cellAssociatedBytes;
  const SEP3D::Core::Status rejected = runtime.BindMesh(bad);
  if (rejected.code != SEP3D::Core::StatusCode::LayoutMismatch ||
      runtime.state() != RM::LifecycleState::Configured) {
    return Fail("mismatched allocation was not rejected transactionally");
  }
  RM::MeshBinding good;
  good.layout = layout;
  if (!runtime.BindMesh(good).ok()) return Fail("exact frozen layout was rejected");

  return Pass("pre-mesh offsets/sizes are frozen; output-only options preserve the physics fingerprint and layout mismatches are atomic");
}

Result RunLIFE3D04() {
  std::shared_ptr<const RM::RunConfiguration3D> standaloneConfiguration;
  std::shared_ptr<const RM::RunConfiguration3D> swmfConfiguration;
  if (!MakeConfiguration(RM::BackgroundAuthority::AnalyticParker,
                         &standaloneConfiguration) ||
      !MakeConfiguration(RM::BackgroundAuthority::Swmf,
                         &swmfConfiguration)) {
    return Fail("host-supplied adapter configurations were rejected");
  }

  RM::Runtime standalone;
  RM::Runtime coupled;
  if (!ConfigureAndBind(&standalone, standaloneConfiguration) ||
      !ConfigureAndBind(&coupled, swmfConfiguration)) {
    return Fail("host-supplied configurations did not reach MeshReady");
  }
  RM::StandaloneAdapter standaloneAdapter;
  RM::SwmfAdapter swmfAdapter;
  if (!standaloneAdapter.Initialize(&standalone).ok() ||
      !standaloneAdapter.PublishFrozenParker(&standalone, 0.0, 10.0, 1).ok() ||
      !swmfAdapter.Initialize(&coupled).ok() ||
      !swmfAdapter.PublishImported(&coupled, 0.0, 10.0, 1, true).ok()) {
    return Fail("standalone and SWMF adapters did not use the common lifecycle");
  }
  if (standalone.state() != RM::LifecycleState::SnapshotReady ||
      coupled.state() != RM::LifecycleState::SnapshotReady) {
    return Fail("one adapter did not reach SnapshotReady");
  }

  // No filename, argc/argv, environment variable, or parser object is passed
  // anywhere in this fixture.  Both front ends receive a validated immutable
  // object from the host and make the identical Bind/Initialize/Publish calls.
  return Pass("host-supplied analytic and SWMF configurations reach SnapshotReady through the same Runtime API without internal parsing");
}

}  // namespace

std::vector<SEP3D::Testing::Descriptor> RegisterRuntimeTests() {
  using Descriptor = SEP3D::Testing::Descriptor;
  using Initialization = SEP3D::Testing::InitializationLevel;
  using RuntimeClass = SEP3D::Testing::RuntimeClass;

  auto make = [](const char* id, const char* name, const char* description,
                 SEP3D::Testing::TestCallback callback) {
    Descriptor descriptor;
    descriptor.id = id;
    descriptor.name = name;
    descriptor.group = "LIFE3D";
    descriptor.description = description;
    descriptor.initialization = Initialization::None;
    descriptor.supportedBuildModes = "standalone-no-AMPS";
    descriptor.runtime = RuntimeClass::Routine;
    descriptor.seedPolicy = "deterministic-no-rng";
    descriptor.stateIsolation = "fresh Runtime and immutable configuration per case";
    descriptor.callback = std::move(callback);
    return descriptor;
  };

  return {
      make("LIFE3D01", "Legal lifecycle",
           "Drives the complete canonical state path and verifies Runtime-owned counters.",
           RunLIFE3D01),
      make("LIFE3D02", "Illegal lifecycle matrix",
           "Exercises every operation from every lifecycle state and requires rejected calls to be atomic.",
           RunLIFE3D02),
      make("LIFE3D03", "Pre-mesh layout and fingerprint",
           "Freezes offsets before allocation and separates physics from output-only identity.",
           RunLIFE3D03),
      make("LIFE3D04", "Host-supplied adapter parity",
           "Runs standalone and SWMF front ends through the same Runtime calls without internal parsing.",
           RunLIFE3D04),
  };
}
