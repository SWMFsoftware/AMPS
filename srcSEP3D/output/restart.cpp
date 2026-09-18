#include "restart.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <limits>
#include <thread>

namespace SEP3D {
namespace Output {
namespace {

namespace fs = std::filesystem;
constexpr char kMagic[8] = {'S','E','P','3','D','R','0','2'};
constexpr std::uint32_t kSchema = 2;
constexpr std::uint64_t kMaximumRecords = UINT64_C(1000000000);

Core::Status Error(const std::string& message) {
  return Core::Status(Core::StatusCode::Error, message);
}

std::uint64_t Hash(const unsigned char* bytes, std::size_t size) {
  std::uint64_t value = UINT64_C(14695981039346656037);
  for (std::size_t i = 0; i < size; ++i) {
    value ^= bytes[i]; value *= UINT64_C(1099511628211);
  }
  return value;
}

bool AddWouldOverflow(std::uint64_t left, std::uint64_t right) {
  return right > std::numeric_limits<std::uint64_t>::max() - left;
}

// Writing fields one-by-one, rather than dumping C++ structs, excludes
// padding, enum width, host endian, and compiler ABI from the restart format.
class Writer {
 public:
  void U32(std::uint32_t value) {
    for (unsigned shift = 0; shift < 32; shift += 8)
      bytes.push_back(static_cast<unsigned char>((value >> shift) & 0xffU));
  }
  void U64(std::uint64_t value) {
    for (unsigned shift = 0; shift < 64; shift += 8)
      bytes.push_back(static_cast<unsigned char>((value >> shift) & 0xffU));
  }
  void I32(std::int32_t value) { U32(static_cast<std::uint32_t>(value)); }
  void Double(double value) {
    std::uint64_t bits = 0; std::memcpy(&bits, &value, sizeof(bits)); U64(bits);
  }
  void Bool(bool value) { bytes.push_back(value ? 1U : 0U); }
  void String(const std::string& value) {
    U64(value.size()); bytes.insert(bytes.end(), value.begin(), value.end());
  }
  std::vector<unsigned char> bytes;
};

class Reader {
 public:
  Reader(const unsigned char* data, std::size_t size) : data_(data), size_(size) {}
  bool U32(std::uint32_t* value) {
    if (!Need(4)) return false;
    *value = 0;
    for (unsigned shift = 0; shift < 32; shift += 8)
      *value |= static_cast<std::uint32_t>(data_[offset_++]) << shift;
    return true;
  }
  bool U64(std::uint64_t* value) {
    if (!Need(8)) return false;
    *value = 0;
    for (unsigned shift = 0; shift < 64; shift += 8)
      *value |= static_cast<std::uint64_t>(data_[offset_++]) << shift;
    return true;
  }
  bool I32(std::int32_t* value) {
    std::uint32_t raw = 0; if (!U32(&raw)) return false;
    *value = static_cast<std::int32_t>(raw); return true;
  }
  bool Double(double* value) {
    std::uint64_t bits = 0; if (!U64(&bits)) return false;
    std::memcpy(value, &bits, sizeof(bits)); return true;
  }
  bool Bool(bool* value) {
    if (!Need(1) || data_[offset_] > 1) return false;
    *value = data_[offset_++] != 0; return true;
  }
  bool String(std::string* value) {
    std::uint64_t length = 0;
    if (!U64(&length) || length > size_ || !Need(static_cast<std::size_t>(length)))
      return false;
    value->assign(reinterpret_cast<const char*>(data_ + offset_),
                  static_cast<std::size_t>(length));
    offset_ += static_cast<std::size_t>(length); return true;
  }
  bool Done() const { return offset_ == size_; }
 private:
  bool Need(std::size_t count) const { return count <= size_ - offset_; }
  const unsigned char* data_ = nullptr;
  std::size_t size_ = 0;
  std::size_t offset_ = 0;
};

void WriteParticle(Writer* out, const Adapters::ParticleRecord& p) {
  out->U64(p.stableId); out->I32(p.species);
  out->Double(p.positionM.x); out->Double(p.positionM.y); out->Double(p.positionM.z);
  out->Double(p.momentumKgMPerS); out->Double(p.mu); out->Double(p.gyrophaseRad);
  out->Double(p.statisticalWeight); out->U64(p.completedStep);
  out->U64(p.substep); out->U64(p.lastShockGeneration);
}

bool ReadParticle(Reader* in, Adapters::ParticleRecord* p) {
  std::int32_t species = -1;
  if (!in->U64(&p->stableId) || !in->I32(&species) ||
      !in->Double(&p->positionM.x) || !in->Double(&p->positionM.y) ||
      !in->Double(&p->positionM.z) || !in->Double(&p->momentumKgMPerS) ||
      !in->Double(&p->mu) || !in->Double(&p->gyrophaseRad) ||
      !in->Double(&p->statisticalWeight) || !in->U64(&p->completedStep) ||
      !in->U64(&p->substep) || !in->U64(&p->lastShockGeneration)) return false;
  p->species = species; return true;
}

void WriteLedger(Writer* out, const Adapters::LedgerRow& row) {
  out->U64(row.key.step); out->I32(row.key.species);
  out->U64(row.activeStart); out->U64(row.injected); out->U64(row.advanced);
  out->U64(row.escaped); out->U64(row.absorbed); out->U64(row.failed);
  out->U64(row.shockCrossings); out->U64(row.activeEnd); out->Bool(row.closed);
}

bool ReadLedger(Reader* in, Adapters::LedgerRow* row) {
  std::int32_t species = -1;
  if (!in->U64(&row->key.step) || !in->I32(&species) ||
      !in->U64(&row->activeStart) || !in->U64(&row->injected) ||
      !in->U64(&row->advanced) || !in->U64(&row->escaped) ||
      !in->U64(&row->absorbed) || !in->U64(&row->failed) ||
      !in->U64(&row->shockCrossings) || !in->U64(&row->activeEnd) ||
      !in->Bool(&row->closed)) return false;
  row->key.species = species; return true;
}

void WriteSourceLedger(Writer* out, const Adapters::SourceLedgerRow& row) {
  out->U64(row.step); out->I32(row.species); out->U64(row.shockGeneration);
  out->U64(row.sourceId);
  out->Double(row.representedParticles); out->Double(row.injectedEnergyJ);
  out->Double(row.injectedMomentumKgMPerS.x);
  out->Double(row.injectedMomentumKgMPerS.y);
  out->Double(row.injectedMomentumKgMPerS.z);
  out->U64(row.macroparticles); out->U64(row.rejected); out->U64(row.capped);
  out->U64(row.inactivePatches); out->U64(row.disconnectedPatches);
}

bool ReadSourceLedger(Reader* in, Adapters::SourceLedgerRow* row) {
  std::int32_t species = -1;
  if (!in->U64(&row->step) || !in->I32(&species) ||
      !in->U64(&row->shockGeneration) || !in->U64(&row->sourceId) ||
      !in->Double(&row->representedParticles) ||
      !in->Double(&row->injectedEnergyJ) ||
      !in->Double(&row->injectedMomentumKgMPerS.x) ||
      !in->Double(&row->injectedMomentumKgMPerS.y) ||
      !in->Double(&row->injectedMomentumKgMPerS.z) ||
      !in->U64(&row->macroparticles) || !in->U64(&row->rejected) ||
      !in->U64(&row->capped) || !in->U64(&row->inactivePatches) ||
      !in->U64(&row->disconnectedPatches)) return false;
  row->species = species;
  return true;
}

void WriteSnapshot(Writer* out, const RuntimeModel::SnapshotDescriptor& value) {
  out->U32(static_cast<std::uint32_t>(value.authority));
  out->Double(value.epochS); out->Double(value.validFromS);
  out->Double(value.validUntilS); out->U64(value.generation);
  out->Bool(value.complete); out->String(value.coordinateFrame);
  out->String(value.providerIdentity);
  out->String(value.configurationFingerprint);
}

bool ReadSnapshot(Reader* in, RuntimeModel::SnapshotDescriptor* value) {
  std::uint32_t authority = 0;
  if (!in->U32(&authority) || authority > 1 ||
      !in->Double(&value->epochS) || !in->Double(&value->validFromS) ||
      !in->Double(&value->validUntilS) || !in->U64(&value->generation) ||
      !in->Bool(&value->complete) || !in->String(&value->coordinateFrame) ||
      !in->String(&value->providerIdentity) ||
      !in->String(&value->configurationFingerprint)) return false;
  value->authority = static_cast<RuntimeModel::BackgroundAuthority>(authority);
  return true;
}

void WriteShock(Writer* out, const Adapters::ShockState& value) {
  out->Bool(value.active); out->U64(value.generation);
  out->Double(value.epochS); out->Double(value.validUntilS);
  out->Double(value.centerM.x); out->Double(value.centerM.y);
  out->Double(value.centerM.z); out->Double(value.radiusM);
  out->Double(value.radialSpeedMPerS); out->Double(value.compressionRatio);
  out->String(value.providerIdentity);
  out->String(value.configurationFingerprint);
}

bool ReadShock(Reader* in, Adapters::ShockState* value) {
  if (!in->Bool(&value->active) || !in->U64(&value->generation) ||
      !in->Double(&value->epochS) || !in->Double(&value->validUntilS) ||
      !in->Double(&value->centerM.x) || !in->Double(&value->centerM.y) ||
      !in->Double(&value->centerM.z) || !in->Double(&value->radiusM) ||
      !in->Double(&value->radialSpeedMPerS) ||
      !in->Double(&value->compressionRatio) ||
      !in->String(&value->providerIdentity) ||
      !in->String(&value->configurationFingerprint)) return false;
  value->status = Core::Status::OK();
  return true;
}

Core::Status Validate(const RestartState& state) {
  if (state.configurationFingerprint.empty() ||
      state.resolvedConfigurationManifest.empty() ||
      state.storageLayoutFingerprint.empty() || state.codeIdentity.empty() ||
      state.snapshotFingerprint.empty() || state.backgroundGeneration == 0 ||
      state.campaignSeed == 0 || state.nextStableParticleId == 0 ||
      state.savedRankCount == 0 || !std::isfinite(state.baseTimeStepS) ||
      state.baseTimeStepS <= 0.0 ||
      state.runtimeCounters.currentTick !=
          state.runtimeCounters.completedSteps ||
      !state.activeSnapshot.complete || state.activeSnapshot.generation !=
          state.backgroundGeneration)
    return Error("restart identity or generation is invalid");
  std::uint64_t previousId = 0;
  for (const Adapters::ParticleRecord& p : state.particles) {
    if (p.stableId == 0 || p.stableId <= previousId || p.species < 0 ||
        !std::isfinite(p.positionM.x) || !std::isfinite(p.positionM.y) ||
        !std::isfinite(p.positionM.z) || !std::isfinite(p.momentumKgMPerS) ||
        p.momentumKgMPerS < 0.0 || !std::isfinite(p.mu) || p.mu < -1.0 ||
        p.mu > 1.0 || !std::isfinite(p.statisticalWeight) ||
        p.statisticalWeight <= 0.0)
      return Error("restart particle table is invalid or not canonically sorted");
    previousId = p.stableId;
  }
  if (state.nextStableParticleId <= previousId)
    return Error("restart next stable particle ID collides with active state");
  for (const Adapters::LedgerRow& row : state.ledgerRows) {
    if (AddWouldOverflow(row.activeStart, row.injected) ||
        AddWouldOverflow(row.activeEnd, row.escaped) ||
        AddWouldOverflow(row.activeEnd + row.escaped, row.absorbed) ||
        AddWouldOverflow(row.activeEnd + row.escaped + row.absorbed,
                         row.failed))
      return Error("restart particle ledger counter overflow");
    const std::uint64_t left = row.activeStart + row.injected;
    const std::uint64_t right = row.activeEnd + row.escaped + row.absorbed + row.failed;
    if (!row.closed || row.key.species < 0 || left != right ||
        row.advanced != row.activeEnd)
      return Error("restart contains an open or unbalanced particle ledger row");
  }
  const Adapters::SourceLedgerRow* previousSource = nullptr;
  for (const Adapters::SourceLedgerRow& row : state.sourceLedgerRows) {
    if (row.species < 0 || !std::isfinite(row.representedParticles) ||
        row.representedParticles < 0.0 || !std::isfinite(row.injectedEnergyJ) ||
        row.injectedEnergyJ < 0.0 ||
        !std::isfinite(row.injectedMomentumKgMPerS.x) ||
        !std::isfinite(row.injectedMomentumKgMPerS.y) ||
        !std::isfinite(row.injectedMomentumKgMPerS.z))
      return Error("restart source ledger is invalid");
    if (previousSource != nullptr &&
        previousSource->step == row.step &&
        previousSource->species == row.species &&
        previousSource->shockGeneration == row.shockGeneration &&
        previousSource->sourceId == row.sourceId)
      return Error("restart contains a duplicate physical source ledger row");
    previousSource = &row;
  }
  return Core::Status::OK();
}

}  // namespace

Core::Status WriteRestart(const std::string& path, const RestartState& state) {
  RestartState canonical = state;
  std::sort(canonical.particles.begin(), canonical.particles.end(),
            [](const Adapters::ParticleRecord& left,
               const Adapters::ParticleRecord& right) {
              return left.stableId < right.stableId;
            });
  std::sort(canonical.ledgerRows.begin(), canonical.ledgerRows.end(),
            [](const Adapters::LedgerRow& left,
               const Adapters::LedgerRow& right) {
              return left.key < right.key;
            });
  std::sort(canonical.sourceLedgerRows.begin(), canonical.sourceLedgerRows.end(),
            [](const Adapters::SourceLedgerRow& left,
               const Adapters::SourceLedgerRow& right) {
              if (left.step != right.step) return left.step < right.step;
              if (left.species != right.species)
                return left.species < right.species;
              if (left.shockGeneration != right.shockGeneration)
                return left.shockGeneration < right.shockGeneration;
              return left.sourceId < right.sourceId;
            });
  const Core::Status valid = Validate(canonical);
  if (!valid.ok()) return valid;
  Writer payload;
  payload.U32(kSchema);
  payload.String(canonical.configurationFingerprint);
  payload.String(canonical.resolvedConfigurationManifest);
  payload.String(canonical.storageLayoutFingerprint);
  payload.String(canonical.codeIdentity);
  payload.String(canonical.snapshotFingerprint);
  payload.U64(canonical.runtimeCounters.completedSteps);
  payload.U64(canonical.runtimeCounters.stepsSinceOutput);
  payload.U64(canonical.runtimeCounters.outputSequence);
  payload.U64(canonical.runtimeCounters.checkpointSequence);
  payload.U64(canonical.runtimeCounters.currentTick);
  payload.U64(canonical.eventSchedule.nextBackgroundTick);
  payload.U64(canonical.eventSchedule.nextInjectionTick);
  payload.U64(canonical.eventSchedule.nextSamplingTick);
  payload.U64(canonical.eventSchedule.nextCheckpointTick);
  WriteSnapshot(&payload, canonical.activeSnapshot);
  payload.Double(canonical.baseTimeStepS);
  payload.U64(canonical.backgroundGeneration);
  payload.U64(canonical.turbulenceGeneration);
  payload.U64(canonical.sourceGeneration);
  payload.U64(canonical.campaignSeed);
  payload.U64(canonical.nextStableParticleId);
  payload.U64(canonical.savedRankCount);
  WriteShock(&payload, canonical.shockState);
  payload.U64(canonical.samplingState.completedSamplings);
  payload.U64(canonical.samplingState.observationsProcessed);
  payload.U64(canonical.samplingState.pendingWindows);
  payload.U64(canonical.samplingState.pendingObservations);
  payload.Double(canonical.samplingState.pendingRepresentedParticles);
  payload.U64(canonical.particles.size());
  for (const auto& particle : canonical.particles) WriteParticle(&payload, particle);
  payload.U64(canonical.ledgerRows.size());
  for (const auto& row : canonical.ledgerRows) WriteLedger(&payload, row);
  payload.U64(canonical.sourceLedgerRows.size());
  for (const auto& row : canonical.sourceLedgerRows)
    WriteSourceLedger(&payload, row);

  const fs::path final(path);
  const fs::path staging(path + ".staging");
  std::error_code ec;
  if (!final.parent_path().empty()) fs::create_directories(final.parent_path(), ec);
  if (ec || fs::exists(staging))
    return Error("restart staging path is unavailable");
  std::ofstream out(staging, std::ios::binary | std::ios::trunc);
  out.write(kMagic, sizeof(kMagic));
  Writer header; header.U64(payload.bytes.size());
  out.write(reinterpret_cast<const char*>(header.bytes.data()), header.bytes.size());
  out.write(reinterpret_cast<const char*>(payload.bytes.data()), payload.bytes.size());
  Writer trailer; trailer.U64(Hash(payload.bytes.data(), payload.bytes.size()));
  out.write(reinterpret_cast<const char*>(trailer.bytes.data()), trailer.bytes.size());
  out.close();
  if (!out.good()) { fs::remove(staging, ec); return Error("restart write failed"); }
  fs::rename(staging, final, ec);
  if (ec) { fs::remove(staging, ec); return Error("atomic restart rename failed"); }
  return Core::Status::OK();
}

Core::Status ReadRestart(const std::string& path,
                         const RestartLoadOptions& options,
                         RestartState* output) {
  if (output == nullptr) return Error("restart output is null");
  std::ifstream input(path, std::ios::binary);
  if (!input) return Error("restart file is absent");
  input.seekg(0, std::ios::end);
  const std::streamoff length = input.tellg(); input.seekg(0, std::ios::beg);
  if (length < 24) return Error("restart file is truncated");
  std::vector<unsigned char> bytes(static_cast<std::size_t>(length));
  input.read(reinterpret_cast<char*>(bytes.data()), length);
  if (!input || std::memcmp(bytes.data(), kMagic, sizeof(kMagic)) != 0)
    return Error("restart magic is invalid");
  Reader header(bytes.data() + 8, 8);
  std::uint64_t payloadSize = 0;
  if (!header.U64(&payloadSize) || payloadSize != bytes.size() - 24)
    return Error("restart payload length is invalid");
  Reader trailer(bytes.data() + 16 + payloadSize, 8);
  std::uint64_t expectedHash = 0;
  if (!trailer.U64(&expectedHash) ||
      Hash(bytes.data() + 16, payloadSize) != expectedHash)
    return Error("restart checksum mismatch");

  Reader in(bytes.data() + 16, static_cast<std::size_t>(payloadSize));
  RestartState candidate;
  std::uint32_t schema = 0;
  if (!in.U32(&schema) || schema != kSchema ||
      !in.String(&candidate.configurationFingerprint) ||
      !in.String(&candidate.resolvedConfigurationManifest) ||
      !in.String(&candidate.storageLayoutFingerprint) ||
      !in.String(&candidate.codeIdentity) ||
      !in.String(&candidate.snapshotFingerprint) ||
      !in.U64(&candidate.runtimeCounters.completedSteps) ||
      !in.U64(&candidate.runtimeCounters.stepsSinceOutput) ||
      !in.U64(&candidate.runtimeCounters.outputSequence) ||
      !in.U64(&candidate.runtimeCounters.checkpointSequence) ||
      !in.U64(&candidate.runtimeCounters.currentTick) ||
      !in.U64(&candidate.eventSchedule.nextBackgroundTick) ||
      !in.U64(&candidate.eventSchedule.nextInjectionTick) ||
      !in.U64(&candidate.eventSchedule.nextSamplingTick) ||
      !in.U64(&candidate.eventSchedule.nextCheckpointTick) ||
      !ReadSnapshot(&in, &candidate.activeSnapshot) ||
      !in.Double(&candidate.baseTimeStepS) ||
      !in.U64(&candidate.backgroundGeneration) ||
      !in.U64(&candidate.turbulenceGeneration) ||
      !in.U64(&candidate.sourceGeneration) || !in.U64(&candidate.campaignSeed) ||
      !in.U64(&candidate.nextStableParticleId) ||
      !in.U64(&candidate.savedRankCount) ||
      !ReadShock(&in, &candidate.shockState) ||
      !in.U64(&candidate.samplingState.completedSamplings) ||
      !in.U64(&candidate.samplingState.observationsProcessed) ||
      !in.U64(&candidate.samplingState.pendingWindows) ||
      !in.U64(&candidate.samplingState.pendingObservations) ||
      !in.Double(&candidate.samplingState.pendingRepresentedParticles))
    return Error("restart header schema is invalid or truncated");
  std::uint64_t count = 0;
  if (!in.U64(&count) || count > kMaximumRecords)
    return Error("restart particle count is invalid");
  candidate.particles.resize(static_cast<std::size_t>(count));
  for (auto& particle : candidate.particles)
    if (!ReadParticle(&in, &particle)) return Error("restart particle is truncated");
  if (!in.U64(&count) || count > kMaximumRecords)
    return Error("restart ledger count is invalid");
  candidate.ledgerRows.resize(static_cast<std::size_t>(count));
  for (auto& row : candidate.ledgerRows)
    if (!ReadLedger(&in, &row)) return Error("restart ledger is truncated");
  if (!in.U64(&count) || count > kMaximumRecords)
    return Error("restart source-ledger count is invalid");
  candidate.sourceLedgerRows.resize(static_cast<std::size_t>(count));
  for (auto& row : candidate.sourceLedgerRows)
    if (!ReadSourceLedger(&in, &row))
      return Error("restart source ledger is truncated");
  if (!in.Done()) return Error("restart contains undeclared trailing payload");

  const Core::Status valid = Validate(candidate);
  if (!valid.ok()) return valid;
  if ((!options.expectedConfigurationFingerprint.empty() &&
       candidate.configurationFingerprint !=
           options.expectedConfigurationFingerprint) ||
      (!options.expectedCodeIdentity.empty() &&
       candidate.codeIdentity != options.expectedCodeIdentity) ||
      (!options.expectedSnapshotFingerprint.empty() &&
       candidate.snapshotFingerprint != options.expectedSnapshotFingerprint))
    return Core::Status(Core::StatusCode::ConfigurationConflict,
                        "restart configuration, code, or snapshot fingerprint mismatch");
  if ((!options.expectedResolvedConfigurationManifest.empty() &&
       candidate.resolvedConfigurationManifest !=
           options.expectedResolvedConfigurationManifest) ||
      (!options.expectedStorageLayoutFingerprint.empty() &&
       candidate.storageLayoutFingerprint !=
           options.expectedStorageLayoutFingerprint))
    return Core::Status(Core::StatusCode::ConfigurationConflict,
                        "restart resolved manifest or storage layout mismatch");
  if (options.currentRankCount == 0 ||
      (options.repartition == RepartitionPolicy::RequireSameRankCount &&
       candidate.savedRankCount != options.currentRankCount))
    return Core::Status(Core::StatusCode::ConfigurationConflict,
                        "restart rank count differs and deterministic repartition is disabled");

  if (candidate.backgroundGeneration != options.availableBackgroundGeneration) {
    if (options.missingSnapshot == MissingSnapshotPolicy::Reject)
      return Core::Status(Core::StatusCode::SnapshotUnavailable,
                          "restart background generation is not available");
    if (!options.snapshotAvailable)
      return Core::Status(Core::StatusCode::SnapshotUnavailable,
                          "wait policy requires a snapshot availability callback");
    const auto deadline = std::chrono::steady_clock::now() +
        std::chrono::milliseconds(options.waitTimeoutMilliseconds);
    while (!options.snapshotAvailable(candidate.backgroundGeneration)) {
      if (std::chrono::steady_clock::now() >= deadline)
        return Core::Status(Core::StatusCode::SnapshotUnavailable,
                            "timed out waiting for restart background generation");
      std::this_thread::sleep_for(std::chrono::milliseconds(
          std::max<std::uint64_t>(1, options.pollMilliseconds)));
    }
  }
  *output = candidate;
  return Core::Status::OK();
}

}  // namespace Output
}  // namespace SEP3D
