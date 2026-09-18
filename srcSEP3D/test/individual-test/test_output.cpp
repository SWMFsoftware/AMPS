#include "sep3d_test_registry.h"

#include "publication.h"
#include "restart.h"
#include "keyed_random.h"

#include <chrono>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <utility>

namespace {

using SEP3D::Testing::Result;
namespace A = SEP3D::Adapters;
namespace C = SEP3D::Core;
namespace O = SEP3D::Output;
namespace T = SEP3D::Transport;
namespace fs = std::filesystem;

Result Pass(const std::string& message) {
  Result result; result.status = SEP3D::Testing::Status::Pass;
  result.message = message; return result;
}
Result Fail(const std::string& message) {
  Result result; result.status = SEP3D::Testing::Status::Fail;
  result.message = message; return result;
}

fs::path UniqueDirectory(const char* tag) {
  const auto tick = std::chrono::high_resolution_clock::now()
      .time_since_epoch().count();
  return fs::temp_directory_path() /
      (std::string("srcsep3d-") + tag + "-" + std::to_string(tick));
}

O::SamplingRequest Request() {
  O::SamplingRequest request;
  request.cells.push_back({1, C::Vec3(10.0, 0.0, 0.0), 2.0});
  O::ParticleObservation a;
  a.stableId = 2; a.cellId = 1; a.species = 0;
  a.positionM = C::Vec3(10.0, 0.0, 0.0);
  a.momentumKgMPerS = 1.0e-19; a.restMassKg = C::Const::m_p;
  a.mu = 0.5; a.statisticalWeight = 3.0;
  O::ParticleObservation b = a;
  b.stableId = 1; b.momentumKgMPerS = 2.0e-19;
  b.mu = -0.25; b.statisticalWeight = 2.0;
  request.particles = {a, b};
  O::VirtualSpacecraftDefinition craft;
  craft.name = "earth"; craft.positionM = a.positionM;
  craft.collectionRadiusM = 1.0;
  craft.kineticEnergyEdgesJ = {0.0, 1.0e-12, 1.0e-10};
  request.spacecraft.push_back(craft);
  O::FieldLineProjectionDefinition line;
  line.name = "parker-0"; line.originM = C::Vec3();
  line.direction = C::Vec3(1.0, 0.0, 0.0);
  line.distanceEdgesM = {0.0, 20.0};
  request.fieldLines.push_back(line);
  A::LedgerRow row; row.key = {3, 0}; row.activeStart = 2;
  row.injected = 1; row.activeEnd = 1; row.escaped = 1;
  row.absorbed = 1; row.advanced = 1; row.shockCrossings = 1; row.closed = true;
  request.ledgerRows.push_back(row);
  request.previousState.completedSamplings = 4;
  request.previousState.observationsProcessed = 9;
  return request;
}

std::string SnapshotText(const O::SamplingSnapshot& s) {
  std::ostringstream out; out.precision(17);
  for (const auto& m : s.cellMoments)
    out << m.cellId << '|' << m.species << '|' << m.representedParticles << '|'
        << m.numberDensityM3 << '|' << m.weightedFluxM2PerS.x << '|'
        << m.kineticEnergyDensityJPerM3 << '|' << m.firstPitchMoment << '\n';
  for (const auto& p : s.spacecraft) {
    out << p.name << '|' << p.species << '|' << p.dipoleAnisotropy;
    for (double value : p.representedParticlesPerJ) out << '|' << value;
    out << '\n';
  }
  for (const auto& p : s.fieldLines) {
    out << p.name << '|' << p.species;
    for (double value : p.representedParticlesPerM) out << '|' << value;
    out << '\n';
  }
  for (const auto& shock : s.shocks)
    out << shock.step << '|' << shock.species << '|' << shock.crossings << '\n';
  out << s.nextState.completedSamplings << '|'
      << s.nextState.observationsProcessed;
  return out.str();
}

O::RestartState State() {
  O::RestartState state;
  state.configurationFingerprint = "cfg-123";
  state.resolvedConfigurationManifest = "resolved-cfg-123";
  state.storageLayoutFingerprint = "layout-123";
  state.codeIdentity = "commit-abc";
  state.snapshotFingerprint = "snapshot-77";
  state.runtimeCounters = {8, 1, 3, 2, 8};
  state.eventSchedule = {9, 9, 10, UINT64_MAX};
  state.activeSnapshot.authority =
      SEP3D::RuntimeModel::BackgroundAuthority::AnalyticParker;
  state.activeSnapshot.epochS = 8.0;
  state.activeSnapshot.validFromS = 0.0;
  state.activeSnapshot.validUntilS = 20.0;
  state.activeSnapshot.generation = 77;
  state.activeSnapshot.complete = true;
  state.activeSnapshot.coordinateFrame = "HCI-like-inertial";
  state.activeSnapshot.providerIdentity = "restart-fixture";
  state.activeSnapshot.configurationFingerprint = "cfg-123";
  state.baseTimeStepS = 1.0;
  state.backgroundGeneration = 77;
  state.turbulenceGeneration = 78;
  state.sourceGeneration = 79;
  state.campaignSeed = 1234;
  state.nextStableParticleId = 100;
  state.savedRankCount = 1;
  state.samplingState = {4, 20, 0, 0, 0.0};
  A::ParticleRecord p;
  p.stableId = 7; p.species = 0; p.positionM = C::Vec3(1.0, 2.0, 3.0);
  p.momentumKgMPerS = 1.0e-19; p.mu = 0.3; p.gyrophaseRad = 0.4;
  p.statisticalWeight = 5.0; p.completedStep = 8; p.substep = 13;
  p.lastShockGeneration = 79;
  A::ParticleRecord q = p; q.stableId = 4; q.mu = -0.2;
  state.particles = {p, q};  // writer must canonicalize this deliberately
  A::LedgerRow row; row.key = {7, 0}; row.activeStart = 2; row.injected = 1;
  row.activeEnd = 1; row.advanced = 1; row.escaped = 1; row.absorbed = 1;
  row.closed = true;
  state.ledgerRows.push_back(row);
  return state;
}

O::RestartLoadOptions LoadOptions() {
  O::RestartLoadOptions options;
  options.expectedConfigurationFingerprint = "cfg-123";
  options.expectedCodeIdentity = "commit-abc";
  options.expectedSnapshotFingerprint = "snapshot-77";
  options.availableBackgroundGeneration = 77;
  return options;
}

Result RunNAT3D06() {
  O::SamplingRequest first = Request();
  const O::SamplingRequest pristine = first;
  const O::SamplingSnapshot ordered = O::Sample(first);
  std::reverse(first.particles.begin(), first.particles.end());
  const O::SamplingSnapshot reversed = O::Sample(first);
  if (!ordered.status.ok() || !reversed.status.ok() ||
      SnapshotText(ordered) != SnapshotText(reversed))
    return Fail("sampling depends on AMPS particle traversal order");
  if (first.particles.size() != pristine.particles.size() ||
      first.particles[0].stableId != pristine.particles[1].stableId ||
      first.particles[1].stableId != pristine.particles[0].stableId)
    return Fail("read-only sampler mutated caller-owned observations");
  const O::SamplingSnapshot repeated = O::Sample(pristine);
  if (SnapshotText(ordered) != SnapshotText(repeated))
    return Fail("repeated sampling changed a deterministic product");
  return Pass("sampling is read-only, stable-ID ordered, and repeatable bit-for-bit");
}

Result RunNAT3D07() {
  const fs::path root = UniqueDirectory("publication");
  const O::SamplingSnapshot snapshot = O::Sample(Request());
  O::PublicationMetadata metadata;
  metadata.sequence = 5; metadata.simulationTimeS = 60.0;
  metadata.snapshotGeneration = 77;
  metadata.configurationFingerprint = "cfg-123";
  metadata.codeIdentity = "commit-abc";
  metadata.snapshotFingerprint = "snapshot-77";
  const O::PublicationResult published =
      O::Publish(root.string(), "sep3d", metadata, snapshot);
  if (!published.status.ok()) { fs::remove_all(root); return Fail(published.status.message); }
  O::ParsedPublication parsed;
  const C::Status verified = O::ParseAndVerifyPublication(
      published.directory, &parsed);
  if (!verified.ok() || parsed.metadata.sequence != 5 ||
      parsed.artifactHashes.size() != 4) {
    fs::remove_all(root); return Fail("independent schema/hash parser rejected output");
  }
  std::ofstream corrupt(fs::path(published.directory) / "cells.csv",
                        std::ios::app | std::ios::binary);
  corrupt << "corrupt\n"; corrupt.close();
  O::ParsedPublication sentinel; sentinel.metadata.sequence = 999;
  const C::Status caught = O::ParseAndVerifyPublication(
      published.directory, &sentinel);
  fs::remove_all(root);
  if (caught.ok() || sentinel.metadata.sequence != 999)
    return Fail("artifact corruption was accepted or parser mutated output");
  return Pass("atomic bundle has unit schemas, identity manifest, and verified artifact hashes");
}

Result RunRST3D01() {
  const fs::path root = UniqueDirectory("restart-roundtrip");
  fs::create_directories(root);
  const fs::path file = root / "checkpoint.bin";
  const O::RestartState state = State();
  const C::Status written = O::WriteRestart(file.string(), state);
  O::RestartState loaded;
  const C::Status read = O::ReadRestart(file.string(), LoadOptions(), &loaded);
  fs::remove_all(root);
  if (!written.ok() || !read.ok() || loaded.particles.size() != 2 ||
      loaded.particles[0].stableId != 4 || loaded.particles[1].stableId != 7 ||
      loaded.runtimeCounters.completedSteps != 8 ||
      loaded.samplingState.observationsProcessed != 20)
    return Fail("complete restart state did not round-trip canonically");
  T::RandomKey beforeKey{state.campaignSeed, 7, 8, 13,
                         T::RandomPurpose::ParkerParallel};
  T::RandomKey afterKey{loaded.campaignSeed, loaded.particles[1].stableId,
                        loaded.particles[1].completedStep,
                        loaded.particles[1].substep,
                        T::RandomPurpose::ParkerParallel};
  T::KeyedRandomStream before(beforeKey), after(afterKey);
  for (int i = 0; i < 8; ++i)
    if (before.Normal01() != after.Normal01())
      return Fail("restored stochastic tuple changed future random draws");

  // Compare one complete next transport step as well as raw variates. This
  // catches an omitted particle field that might not affect the first random
  // number but would still diverge the restarted physical trajectory.
  A::MoverInput uninterrupted;
  uninterrupted.particle = state.particles[0];  // stable ID 7
  uninterrupted.model = SEP3D::RuntimeModel::TransportModel::Parker3D;
  uninterrupted.local.background.status = C::Status::OK();
  uninterrupted.local.background.valid = true;
  uninterrupted.local.background.B = C::Vec3(1.0, 0.0, 0.0);
  uninterrupted.local.background.absB = 1.0;
  uninterrupted.local.background.bHat = C::Vec3(1.0, 0.0, 0.0);
  uninterrupted.local.background.U = C::Vec3(2.0, 0.0, 0.0);
  uninterrupted.local.cellSizeM = 1.0e8;
  uninterrupted.local.kappaParallelM2PerS = 1.0e6;
  uninterrupted.speciesMassKg = C::Const::m_p;
  uninterrupted.requestedDtS = 0.25;
  uninterrupted.innerRadiusM = 0.1;
  uninterrupted.outerRadiusM = 1.0e9;
  uninterrupted.campaignSeed = state.campaignSeed;
  A::MoverInput restarted = uninterrupted;
  restarted.particle = loaded.particles[1];
  restarted.campaignSeed = loaded.campaignSeed;
  const A::MoverResult nextA = A::AdvanceParticle(uninterrupted);
  const A::MoverResult nextB = A::AdvanceParticle(restarted);
  if (!nextA.status.ok() || !nextB.status.ok() ||
      !(nextA.particle.positionM == nextB.particle.positionM) ||
      nextA.particle.momentumKgMPerS != nextB.particle.momentumKgMPerS ||
      nextA.particle.substep != nextB.particle.substep)
    return Fail("uninterrupted and restarted next transport steps diverged");
  return Pass("runtime, sampling, ledger, particle, and next stochastic transport step round-trip exactly");
}

Result RunRST3D02() {
  const fs::path root = UniqueDirectory("restart-reject");
  fs::create_directories(root); const fs::path file = root / "checkpoint.bin";
  if (!O::WriteRestart(file.string(), State()).ok()) {
    fs::remove_all(root); return Fail("could not create restart fixture");
  }
  O::RestartLoadOptions mismatch = LoadOptions();
  mismatch.expectedConfigurationFingerprint = "different";
  O::RestartState sentinel; sentinel.nextStableParticleId = 999;
  const C::Status fingerprint = O::ReadRestart(file.string(), mismatch, &sentinel);
  std::fstream bytes(file, std::ios::in | std::ios::out | std::ios::binary);
  bytes.seekp(20); char value = 0; bytes.read(&value, 1); bytes.seekp(20);
  value ^= 0x1; bytes.write(&value, 1); bytes.close();
  const C::Status checksum = O::ReadRestart(file.string(), LoadOptions(), &sentinel);
  fs::remove_all(root);
  if (fingerprint.code != C::StatusCode::ConfigurationConflict ||
      checksum.ok() || sentinel.nextStableParticleId != 999)
    return Fail("mismatch/corruption did not reject transactionally");
  return Pass("fingerprint and checksum failures leave destination state untouched");
}

Result RunRST3D03() {
  const fs::path root = UniqueDirectory("restart-snapshot");
  fs::create_directories(root); const fs::path file = root / "checkpoint.bin";
  if (!O::WriteRestart(file.string(), State()).ok()) {
    fs::remove_all(root); return Fail("could not create snapshot-policy fixture");
  }
  O::RestartLoadOptions reject = LoadOptions();
  reject.availableBackgroundGeneration = 76;
  O::RestartState output;
  const C::Status rejected = O::ReadRestart(file.string(), reject, &output);
  O::RestartLoadOptions wait = reject;
  wait.missingSnapshot = O::MissingSnapshotPolicy::Wait;
  wait.waitTimeoutMilliseconds = 5;
  wait.snapshotAvailable = [](std::uint64_t generation) {
    return generation == 77;
  };
  const C::Status accepted = O::ReadRestart(file.string(), wait, &output);
  fs::remove_all(root);
  if (rejected.code != C::StatusCode::SnapshotUnavailable ||
      !accepted.ok() || output.backgroundGeneration != 77)
    return Fail("missing-snapshot reject/wait policy was not enforced");
  return Pass("restart rejects a missing generation or waits for explicit immutable publication");
}

}  // namespace

std::vector<SEP3D::Testing::Descriptor> RegisterOutputTests() {
  using D = SEP3D::Testing::Descriptor;
  using I = SEP3D::Testing::InitializationLevel;
  using RC = SEP3D::Testing::RuntimeClass;
  auto make = [](const char* id, const char* group, const char* name,
                 SEP3D::Testing::TestCallback callback) {
    D d; d.id = id; d.name = name; d.group = group;
    d.description = "Phase-O deterministic sampling/output/restart acceptance";
    d.initialization = I::None; d.supportedBuildModes = "standalone-no-AMPS";
    d.runtime = RC::Routine; d.seedPolicy = "restored semantic key tuple";
    d.stateIsolation = "unique temporary publication per callback";
    d.callback = std::move(callback); return d;
  };
  return {
      make("NAT3D06", "NAT3D", "Sampling isolation", RunNAT3D06),
      make("NAT3D07", "NAT3D", "Output schema and publication", RunNAT3D07),
      make("RST3D01", "RST3D", "Complete restart round trip", RunRST3D01),
      make("RST3D02", "RST3D", "Transactional restart rejection", RunRST3D02),
      make("RST3D03", "RST3D", "Snapshot restart policy", RunRST3D03),
  };
}
