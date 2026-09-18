#include "sep_runtime_contracts.h"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <limits>
#include <set>
#include <sstream>

namespace SEP {
namespace RuntimeContracts {
namespace {

Transport::Status Error(Transport::StatusCode code,
                        const std::string& message) {
  return Transport::Status::Error(code, message);
}

bool Finite(double value) { return std::isfinite(value); }

bool CellGeometryValid(const Turbulence::CellState& cell) {
  return Finite(cell.lengthM) && cell.lengthM > 0.0 &&
         Finite(cell.volumeM3) && cell.volumeM3 > 0.0 &&
         Finite(cell.plasmaSpeedMPerS) &&
         Finite(cell.alfvenSpeedMPerS) && cell.alfvenSpeedMPerS >= 0.0 &&
         Finite(cell.dLnAlfvenSpeeddsPerM) &&
         Finite(cell.magneticFieldT) && cell.magneticFieldT > 0.0 &&
         Finite(cell.massDensityKgPerM3) && cell.massDensityKgPerM3 > 0.0;
}

bool HasPending(const Turbulence::State& state) {
  for (std::size_t i = 0; i < state.cells.size(); ++i) {
    const Turbulence::CellState& c = state.cells[i];
    if (c.pendingParticlePlusJ != 0.0 || c.pendingParticleMinusJ != 0.0 ||
        c.pendingShockPlusJ != 0.0 || c.pendingShockMinusJ != 0.0)
      return true;
  }
  return false;
}

std::string HexEncode(const std::string& value) {
  static const char digits[] = "0123456789abcdef";
  std::string encoded;
  encoded.reserve(value.size() * 2);
  for (std::size_t i = 0; i < value.size(); ++i) {
    const unsigned char c = static_cast<unsigned char>(value[i]);
    encoded.push_back(digits[c >> 4]);
    encoded.push_back(digits[c & 0x0f]);
  }
  return encoded;
}

int HexDigit(char value) {
  if (value >= '0' && value <= '9') return value - '0';
  if (value >= 'a' && value <= 'f') return value - 'a' + 10;
  if (value >= 'A' && value <= 'F') return value - 'A' + 10;
  return -1;
}

bool HexDecode(const std::string& encoded, std::string* value) {
  if (!value || encoded.size() % 2 != 0) return false;
  value->clear();
  value->reserve(encoded.size() / 2);
  for (std::size_t i = 0; i < encoded.size(); i += 2) {
    const int high = HexDigit(encoded[i]);
    const int low = HexDigit(encoded[i + 1]);
    if (high < 0 || low < 0) return false;
    value->push_back(static_cast<char>((high << 4) | low));
  }
  return true;
}

std::uint64_t HashAppend(std::uint64_t hash, std::uint64_t value) {
  // FNV-1a over an explicit little-endian integer representation avoids
  // dependence on host structure padding and compiler ABI.
  for (unsigned i = 0; i < 8; ++i) {
    hash ^= static_cast<unsigned char>((value >> (8 * i)) & UINT64_C(0xff));
    hash *= UINT64_C(1099511628211);
  }
  return hash;
}

bool TransactionLess(const CouplingTransaction& a,
                     const CouplingTransaction& b) {
  if (a.fieldLineId != b.fieldLineId) return a.fieldLineId < b.fieldLineId;
  if (a.cell != b.cell) return a.cell < b.cell;
  if (a.branch != b.branch) return a.branch < b.branch;
  if (a.spectralBin != b.spectralBin) return a.spectralBin < b.spectralBin;
  if (a.particleId != b.particleId) return a.particleId < b.particleId;
  if (a.event != b.event) return a.event < b.event;
  if (a.interval != b.interval) return a.interval < b.interval;
  return a.transactionId < b.transactionId;
}

Transport::Status ValidateTransaction(const CouplingTransaction& t,
                                      const TurbulenceRuntimeStore& store) {
  if (t.schema != UINT64_C(0x5352435345504301) || t.transactionId == 0 ||
      t.particleId == 0 || t.generation == 0 || t.fieldLineId < 0 ||
      (t.branch != -1 && t.branch != 1) ||
      !Finite(t.particleEnergyChangeJ) || t.turbulenceIdentity.empty())
    return Error(Transport::StatusCode::InvalidArgument,
                 "coupling transaction has invalid typed metadata");
  const TurbulenceLineRecord* line = store.Find(t.fieldLineId);
  if (!line)
    return Error(Transport::StatusCode::OutOfDomain,
                 "coupling transaction names an unknown field line");
  if (line->geometryGeneration != t.generation ||
      line->state.fieldLineGeneration != t.generation)
    return Error(Transport::StatusCode::InvalidParticleState,
                 "coupling transaction belongs to a stale geometry generation");
  std::ostringstream expectedIdentity;
  expectedIdentity << line->configurationFingerprint << ":g"
                   << line->geometryGeneration << ":t"
                   << line->state.epochS;
  if (line->state.provenance.find(t.turbulenceIdentity) == std::string::npos &&
      line->state.sourceChecksum != t.turbulenceIdentity &&
      expectedIdentity.str() != t.turbulenceIdentity)
    return Error(Transport::StatusCode::InvalidParticleState,
                 "coupling transaction turbulence identity does not match owner");
  if (t.cell >= line->state.cells.size())
    return Error(Transport::StatusCode::OutOfDomain,
                 "coupling transaction cell is outside its field line");
  if (line->state.configuration.representation ==
          Turbulence::Representation::Spectral &&
      t.spectralBin >= line->state.configuration.spectralBins)
    return Error(Transport::StatusCode::OutOfDomain,
                 "coupling transaction spectral bin is outside the grid");
  return Transport::Status::Ok();
}

}  // namespace

Transport::Status TurbulenceRuntimeStore::Install(
    const TurbulenceLineRecord& record) {
  if (record.fieldLineId < 0 || record.ownerRank < 0 ||
      record.geometryGeneration == 0 ||
      record.configurationFingerprint.empty() || record.state.cells.empty())
    return Error(Transport::StatusCode::InvalidArgument,
                 "runtime turbulence owner metadata is incomplete");
  if (record.state.fieldLineGeneration != record.geometryGeneration)
    return Error(Transport::StatusCode::InvalidArgument,
                 "runtime store and turbulence generation disagree");
  for (std::size_t i = 0; i < record.state.cells.size(); ++i)
    if (!CellGeometryValid(record.state.cells[i]))
      return Error(Transport::StatusCode::InvalidParticleState,
                   "runtime turbulence state has invalid SI geometry");
  if (lines_.count(record.fieldLineId) != 0)
    return Error(Transport::StatusCode::InvalidArgument,
                 "field line already has an authoritative turbulence owner");
  lines_[record.fieldLineId] = record;
  return Transport::Status::Ok();
}

bool TurbulenceRuntimeStore::Contains(int fieldLineId) const {
  return lines_.count(fieldLineId) != 0;
}

const TurbulenceLineRecord* TurbulenceRuntimeStore::Find(
    int fieldLineId) const {
  std::map<int, TurbulenceLineRecord>::const_iterator found =
      lines_.find(fieldLineId);
  return found == lines_.end() ? NULL : &found->second;
}

TurbulenceLineRecord* TurbulenceRuntimeStore::FindMutable(int fieldLineId) {
  std::map<int, TurbulenceLineRecord>::iterator found = lines_.find(fieldLineId);
  return found == lines_.end() ? NULL : &found->second;
}

Transport::Status TurbulenceRuntimeStore::RefreshBackground(
    int fieldLineId, std::uint64_t geometryGeneration,
    const std::string& configurationFingerprint,
    const std::vector<Turbulence::CellState>& geometryAndBackground,
    PendingRemapPolicy pendingPolicy,
    Turbulence::EnergyLedger* remapLedger) {
  TurbulenceLineRecord* record = FindMutable(fieldLineId);
  if (!record)
    return Error(Transport::StatusCode::OutOfDomain,
                 "cannot refresh an uninitialized turbulence field line");
  if (geometryGeneration == 0 || configurationFingerprint.empty() ||
      geometryAndBackground.empty())
    return Error(Transport::StatusCode::InvalidArgument,
                 "background refresh metadata is incomplete");
  if (record->configurationFingerprint != configurationFingerprint)
    return Error(Transport::StatusCode::UnsupportedConfiguration,
                 "turbulence configuration fingerprint changed in flight");
  for (std::size_t i = 0; i < geometryAndBackground.size(); ++i)
    if (!CellGeometryValid(geometryAndBackground[i]))
      return Error(Transport::StatusCode::InvalidParticleState,
                   "background refresh contains invalid geometry or coefficients");

  if (geometryGeneration == record->geometryGeneration) {
    if (geometryAndBackground.size() != record->state.cells.size())
      return Error(Transport::StatusCode::InvalidParticleState,
                   "same-generation turbulence topology changed size");
    // Only fields derived from the immutable background are refreshed.  Wave
    // energy, spectra, and pending sources remain owned by the runtime record.
    for (std::size_t i = 0; i < geometryAndBackground.size(); ++i) {
      Turbulence::CellState& target = record->state.cells[i];
      const Turbulence::CellState& source = geometryAndBackground[i];
      target.lengthM = source.lengthM;
      target.volumeM3 = source.volumeM3;
      target.plasmaSpeedMPerS = source.plasmaSpeedMPerS;
      target.alfvenSpeedMPerS = source.alfvenSpeedMPerS;
      target.dLnAlfvenSpeeddsPerM = source.dLnAlfvenSpeeddsPerM;
      target.magneticFieldT = source.magneticFieldT;
      target.massDensityKgPerM3 = source.massDensityKgPerM3;
    }
    if (remapLedger) *remapLedger = Turbulence::EnergyLedger();
    return Transport::Status::Ok();
  }
  if (geometryGeneration < record->geometryGeneration)
    return Error(Transport::StatusCode::InvalidParticleState,
                 "turbulence geometry generation moved backwards");
  if (HasPending(record->state) &&
      pendingPolicy == PendingRemapPolicy::RejectPending)
    return Error(Transport::StatusCode::InvalidParticleState,
                 "pending source transaction must be settled before remap");

  Turbulence::State remapped;
  Turbulence::EnergyLedger localLedger;
  Transport::Status status = Turbulence::RemapConservatively(
      record->state, geometryAndBackground, &remapped, &localLedger);
  if (!status.ok()) return status;
  remapped.fieldLineGeneration = geometryGeneration;
  remapped.epochS = record->state.epochS;
  remapped.completedSteps = record->state.completedSteps;
  remapped.handoffCompleted = record->state.handoffCompleted;
  remapped.handoffEpochS = record->state.handoffEpochS;
  remapped.sourceChecksum = record->state.sourceChecksum;
  remapped.provenance = record->state.provenance;
  remapped.accumulatedLedger = record->state.accumulatedLedger;
  record->state = remapped;
  record->geometryGeneration = geometryGeneration;
  if (remapLedger) *remapLedger = localLedger;
  return Transport::Status::Ok();
}

void TurbulenceRuntimeStore::Clear() { lines_.clear(); }

Transport::Status TurbulenceRuntimeStore::Serialize(std::string* text) const {
  if (!text)
    return Error(Transport::StatusCode::InvalidArgument,
                 "runtime checkpoint destination is null");
  std::ostringstream out;
  out << "SEP_TURBULENCE_RUNTIME 1\n" << lines_.size() << '\n';
  for (std::map<int, TurbulenceLineRecord>::const_iterator i = lines_.begin();
       i != lines_.end(); ++i) {
    std::string stateText;
    const Transport::Status status =
        Turbulence::SerializeCheckpoint(i->second.state, &stateText);
    if (!status.ok()) return status;
    out << i->second.fieldLineId << ' ' << i->second.ownerRank << ' '
        << i->second.geometryGeneration << ' '
        << HexEncode(i->second.configurationFingerprint) << ' '
        << HexEncode(stateText) << '\n';
  }
  *text = out.str();
  return Transport::Status::Ok();
}

Transport::Status TurbulenceRuntimeStore::Deserialize(
    const std::string& text) {
  std::istringstream in(text);
  std::string magic;
  int version = 0;
  std::size_t count = 0;
  if (!(in >> magic >> version >> count) ||
      magic != "SEP_TURBULENCE_RUNTIME" || version != 1)
    return Error(Transport::StatusCode::InvalidArgument,
                 "runtime turbulence checkpoint header is invalid");
  TurbulenceRuntimeStore staged;
  for (std::size_t i = 0; i < count; ++i) {
    TurbulenceLineRecord record;
    std::string fingerprintHex;
    std::string stateHex;
    std::string stateText;
    if (!(in >> record.fieldLineId >> record.ownerRank >>
          record.geometryGeneration >> fingerprintHex >> stateHex) ||
        !HexDecode(fingerprintHex, &record.configurationFingerprint) ||
        !HexDecode(stateHex, &stateText))
      return Error(Transport::StatusCode::InvalidArgument,
                   "runtime turbulence checkpoint is truncated");
    Transport::Status status =
        Turbulence::DeserializeCheckpoint(stateText, &record.state);
    if (!status.ok()) return status;
    status = staged.Install(record);
    if (!status.ok()) return status;
  }
  std::string trailing;
  if (in >> trailing)
    return Error(Transport::StatusCode::InvalidArgument,
                 "runtime turbulence checkpoint has trailing data");
  lines_.swap(staged.lines_);
  return Transport::Status::Ok();
}

CouplingBatchResult ApplyCouplingTransactions(
    const std::vector<CouplingTransaction>& transactions,
    InvalidTransactionPolicy policy,
    TurbulenceRuntimeStore* store) {
  CouplingBatchResult result;
  if (!store) {
    result.status = Error(Transport::StatusCode::InvalidArgument,
                          "coupling transaction store is null");
    return result;
  }

  std::vector<CouplingTransaction> accepted;
  accepted.reserve(transactions.size());
  std::set<std::uint64_t> ids;
  for (std::size_t i = 0; i < transactions.size(); ++i) {
    const Transport::Status valid = ValidateTransaction(transactions[i], *store);
    const bool duplicate = !ids.insert(transactions[i].transactionId).second;
    if (!valid.ok() || duplicate) {
      if (policy == InvalidTransactionPolicy::RejectBatch) {
        result.status = duplicate
            ? Error(Transport::StatusCode::InvalidArgument,
                    "coupling batch contains a duplicate transaction identity")
            : valid;
        return result;
      }
      ++result.quarantined;
      result.quarantinedTransactionIds.push_back(
          transactions[i].transactionId);
      continue;
    }
    accepted.push_back(transactions[i]);
  }
  std::sort(accepted.begin(), accepted.end(), TransactionLess);

  // Stage only the affected field-line records.  The production store is not
  // modified until every aggregated delta passes finite-value and positivity
  // checks, providing a true all-or-nothing transaction boundary.
  std::map<int, Turbulence::State> staged;
  for (std::size_t i = 0; i < accepted.size(); ++i) {
    const CouplingTransaction& t = accepted[i];
    if (staged.count(t.fieldLineId) == 0)
      staged[t.fieldLineId] = store->Find(t.fieldLineId)->state;
    Turbulence::CellState& cell = staged[t.fieldLineId].cells[t.cell];
    const double waveDelta = -t.particleEnergyChangeJ;
    double* pending = t.branch > 0
        ? &cell.pendingParticlePlusJ : &cell.pendingParticleMinusJ;
    const double available = t.branch > 0 ? cell.ePlusJ : cell.eMinusJ;
    if (!Finite(waveDelta) || !Finite(*pending + waveDelta) ||
        available + *pending + waveDelta < 0.0) {
      result.status = Error(Transport::StatusCode::InvalidParticleState,
          "coupling batch would make a resonant wave branch negative");
      return result;
    }
    *pending += waveDelta;
    result.particleEnergyChangeJ += t.particleEnergyChangeJ;
    result.waveEnergyChangeJ += waveDelta;
    ++result.accepted;
  }
  for (std::map<int, Turbulence::State>::const_iterator i = staged.begin();
       i != staged.end(); ++i)
    store->FindMutable(i->first)->state = i->second;
  result.status = Transport::Status::Ok();
  return result;
}

VelocityDerivatives ComputeVelocityDerivatives(
    const VelocityGradientInput& input) {
  VelocityDerivatives result;
  double b2 = 0.0;
  for (int i = 0; i < 3; ++i) {
    if (!Finite(input.velocityMPerS[i]) || !Finite(input.bUnit[i]) ||
        !Finite(input.curvaturePerM[i])) {
      result.status = Error(Transport::StatusCode::InvalidArgument,
                            "velocity derivative vector input is non-finite");
      return result;
    }
    b2 += input.bUnit[i] * input.bUnit[i];
    for (int j = 0; j < 3; ++j)
      if (!Finite(input.gradientPerS[i][j])) {
        result.status = Error(Transport::StatusCode::InvalidArgument,
                              "velocity gradient tensor is non-finite");
        return result;
      }
  }
  if (std::fabs(b2 - 1.0) > 1.0e-10) {
    result.status = Error(Transport::StatusCode::InvalidArgument,
                          "magnetic-field direction must be a unit vector");
    return result;
  }
  for (int i = 0; i < 3; ++i) {
    result.divergencePerS += input.gradientPerS[i][i];
    result.curvatureContributionPerS +=
        input.velocityMPerS[i] * input.curvaturePerM[i];
    for (int j = 0; j < 3; ++j)
      result.fieldAlignedStrainPerS +=
          input.bUnit[i] * input.gradientPerS[i][j] * input.bUnit[j];
  }
  result.parallelVelocityGradientPerS =
      result.fieldAlignedStrainPerS + result.curvatureContributionPerS;
  result.method = "provider-gradient-tensor+field-line-curvature";
  result.status = Transport::Status::Ok();
  return result;
}

Transport::ScalarResult ContinuityResidualPerS(
    double materialLogDensityDerivativePerS,
    double divergencePerS) {
  Transport::ScalarResult result;
  if (!Finite(materialLogDensityDerivativePerS) || !Finite(divergencePerS)) {
    result.status = Error(Transport::StatusCode::InvalidArgument,
                          "continuity residual inputs must be finite");
    return result;
  }
  result.value = materialLogDensityDerivativePerS + divergencePerS;
  result.status = Transport::Status::Ok();
  return result;
}

ShockValidation ValidateShockState(const ShockState& state) {
  ShockValidation result;
  double normal2 = 0.0;
  for (int i = 0; i < 3; ++i) {
    if (!Finite(state.normal[i])) {
      result.status = Error(Transport::StatusCode::InvalidArgument,
                            "shock normal is non-finite");
      result.reason = result.status.message;
      return result;
    }
    normal2 += state.normal[i] * state.normal[i];
  }
  if (!Finite(state.epochS) || std::fabs(normal2 - 1.0) > 1.0e-10 ||
      !Finite(state.shockNormalSpeedMPerS) ||
      !Finite(state.upstreamNormalSpeedMPerS) ||
      !Finite(state.downstreamNormalSpeedMPerS) ||
      !Finite(state.upstreamNumberDensityPerM3) ||
      !Finite(state.downstreamNumberDensityPerM3) ||
      state.upstreamNumberDensityPerM3 <= 0.0 ||
      state.downstreamNumberDensityPerM3 <= 0.0 ||
      !Finite(state.compressionRatio) || state.compressionRatio <= 0.0 ||
      !Finite(state.alfvenMach) || !Finite(state.sonicMach) ||
      !Finite(state.obliquityRad) || state.provider.empty() ||
      state.provenance.empty()) {
    result.status = Error(Transport::StatusCode::InvalidArgument,
                          "shock state is incomplete or non-physical");
    result.reason = result.status.message;
    return result;
  }
  if (state.compressionRatio <= 1.0 ||
      state.shockNormalSpeedMPerS <= state.upstreamNormalSpeedMPerS) {
    result.status = Transport::Status::Ok();
    result.quality = ShockQuality::NoCompressiveShock;
    result.reason = "upstream-relative normal flux is not compressive";
    return result;
  }
  const double densityRatio = state.downstreamNumberDensityPerM3 /
                              state.upstreamNumberDensityPerM3;
  if (std::fabs(densityRatio - state.compressionRatio) >
      1.0e-6 * std::max(1.0, state.compressionRatio)) {
    result.status = Error(Transport::StatusCode::InvalidCoefficient,
                          "shock density jump and compression ratio disagree");
    result.reason = result.status.message;
    return result;
  }
  result.dsaPhaseSpaceIndex = 3.0 * state.compressionRatio /
                              (state.compressionRatio - 1.0);
  result.quality = ShockQuality::ValidCompressive;
  result.status = Transport::Status::Ok();
  result.reason = "validated compressive shock";
  return result;
}

Transport::ScalarResult ProcessedUpstreamParticleCount(
    const ShockState& state, double shockAreaM2, double dtS,
    double speciesInjectionEfficiency) {
  Transport::ScalarResult result;
  const ShockValidation valid = ValidateShockState(state);
  if (!valid.status.ok()) {
    result.status = valid.status;
    return result;
  }
  if (!Finite(shockAreaM2) || shockAreaM2 < 0.0 || !Finite(dtS) || dtS < 0.0 ||
      !Finite(speciesInjectionEfficiency) || speciesInjectionEfficiency < 0.0 ||
      speciesInjectionEfficiency > 1.0) {
    result.status = Error(Transport::StatusCode::InvalidArgument,
                          "shock source area, time, or efficiency is invalid");
    return result;
  }
  if (valid.quality == ShockQuality::NoCompressiveShock) {
    result.value = 0.0;
    result.status = Transport::Status::Ok();
    return result;
  }
  const double relativeNormalSpeed = std::max(
      0.0, state.shockNormalSpeedMPerS - state.upstreamNormalSpeedMPerS);
  result.value = state.upstreamNumberDensityPerM3 * relativeNormalSpeed *
                 shockAreaM2 * dtS * speciesInjectionEfficiency;
  result.status = Finite(result.value)
      ? Transport::Status::Ok()
      : Error(Transport::StatusCode::InvalidParticleState,
              "shock processed-particle count overflowed");
  return result;
}

Transport::ScalarResult CompressionFromMach(double mach, double gamma) {
  Transport::ScalarResult result;
  if (!Finite(mach) || !Finite(gamma) || mach <= 0.0 || gamma <= 1.0) {
    result.status = Error(Transport::StatusCode::InvalidArgument,
                          "compression requires positive Mach and gamma>1");
    return result;
  }
  if (mach <= 1.0) {
    result.value = 1.0;
    result.status = Transport::Status::Ok();
    return result;
  }
  result.value = (gamma + 1.0) * mach * mach /
      ((gamma - 1.0) * mach * mach + 2.0);
  result.status = Transport::Status::Ok();
  return result;
}

std::uint64_t StableSourceIdentity(const SourceEventKey& key) {
  std::uint64_t hash = UINT64_C(1469598103934665603);
  hash = HashAppend(hash, key.schema);
  hash = HashAppend(hash, key.campaign);
  hash = HashAppend(hash, key.source);
  hash = HashAppend(hash, key.event);
  hash = HashAppend(hash, key.fieldLine);
  hash = HashAppend(hash, key.species);
  hash = HashAppend(hash, key.ordinal);
  return hash == 0 ? UINT64_C(1) : hash;
}

Transport::KeyedRandomStream SourceRandomStream(
    const SourceEventKey& key, SourceRandomPurpose purpose) {
  return Transport::KeyedRandomStream(
      key.campaign, StableSourceIdentity(key),
      static_cast<std::uint64_t>(purpose), key.event);
}

Transport::Status SourceEventScheduler::Serialize(std::string* text) const {
  if (!text)
    return Error(Transport::StatusCode::InvalidArgument,
                 "source scheduler checkpoint destination is null");
  std::ostringstream out;
  out << "SEP_SOURCE_SCHEDULER 2\n" << campaignSeed_ << ' ' << nextEvent_ << '\n';
  *text = out.str();
  return Transport::Status::Ok();
}

Transport::Status SourceEventScheduler::Deserialize(const std::string& text) {
  std::istringstream in(text);
  std::string magic;
  int version = 0;
  std::uint64_t campaign = 0;
  std::uint64_t next = 0;
  std::string trailing;
  if (!(in >> magic >> version >> campaign >> next) ||
      magic != "SEP_SOURCE_SCHEDULER" || version != 2 || (in >> trailing))
    return Error(Transport::StatusCode::InvalidArgument,
                 "source scheduler checkpoint is invalid");
  campaignSeed_ = campaign;
  nextEvent_ = next;
  return Transport::Status::Ok();
}

}  // namespace RuntimeContracts
}  // namespace SEP
