#include "particle_ledger.h"

#include <limits>

namespace SEP3D {
namespace Adapters {
namespace {

Core::Status Invalid(const char* message) {
  return Core::Status(Core::StatusCode::InvalidInput, message);
}

bool AddWouldOverflow(std::uint64_t left, std::uint64_t right) {
  return right > std::numeric_limits<std::uint64_t>::max() - left;
}

}  // namespace

Core::Status ParticleLedger::Begin(std::uint64_t step, int species,
                                   std::uint64_t activeStart) {
  if (species < 0) return Invalid("ledger species is negative");
  const LedgerKey key{step, species};
  if (rows_.find(key) != rows_.end())
    return Invalid("ledger row already exists");
  LedgerRow row; row.key = key; row.activeStart = activeStart;
  rows_.emplace(key, row);
  return Core::Status::OK();
}

Core::Status ParticleLedger::RecordInjection(std::uint64_t step, int species,
                                             std::uint64_t count) {
  auto found = rows_.find(LedgerKey{step, species});
  if (found == rows_.end() || found->second.closed)
    return Invalid("injection requires an open ledger row");
  if (AddWouldOverflow(found->second.injected, count))
    return Invalid("injected particle counter overflow");
  found->second.injected += count;
  return Core::Status::OK();
}

Core::Status ParticleLedger::RecordMover(std::uint64_t step, int species,
                                         ParticleDisposition disposition,
                                         bool shockCrossed) {
  auto found = rows_.find(LedgerKey{step, species});
  if (found == rows_.end() || found->second.closed)
    return Invalid("mover outcome requires an open ledger row");
  // Update a copy and commit only after every counter (including the optional
  // shock diagnostic) has passed its overflow guard. A failed diagnostic must
  // not leave the primary disposition count partially advanced.
  LedgerRow candidate = found->second;
  switch (disposition) {
    case ParticleDisposition::Active:
      if (candidate.advanced == std::numeric_limits<std::uint64_t>::max())
        return Invalid("advanced particle counter overflow");
      ++candidate.advanced; break;
    case ParticleDisposition::Escaped:
      if (candidate.escaped == std::numeric_limits<std::uint64_t>::max())
        return Invalid("escaped particle counter overflow");
      ++candidate.escaped; break;
    case ParticleDisposition::Absorbed:
      if (candidate.absorbed == std::numeric_limits<std::uint64_t>::max())
        return Invalid("absorbed particle counter overflow");
      ++candidate.absorbed; break;
    case ParticleDisposition::Failed:
      if (candidate.failed == std::numeric_limits<std::uint64_t>::max())
        return Invalid("failed particle counter overflow");
      ++candidate.failed; break;
  }
  if (shockCrossed) {
    if (candidate.shockCrossings == std::numeric_limits<std::uint64_t>::max())
      return Invalid("shock-crossing counter overflow");
    ++candidate.shockCrossings;
  }
  found->second = candidate;
  return Core::Status::OK();
}

Core::Status ParticleLedger::Close(std::uint64_t step, int species,
                                   std::uint64_t activeEnd) {
  auto found = rows_.find(LedgerKey{step, species});
  if (found == rows_.end() || found->second.closed)
    return Invalid("close requires an open ledger row");
  LedgerRow candidate = found->second;
  candidate.activeEnd = activeEnd;
  if (AddWouldOverflow(candidate.activeStart, candidate.injected) ||
      AddWouldOverflow(candidate.activeEnd, candidate.escaped) ||
      AddWouldOverflow(candidate.activeEnd + candidate.escaped,
                       candidate.absorbed) ||
      AddWouldOverflow(candidate.activeEnd + candidate.escaped +
                           candidate.absorbed,
                       candidate.failed))
    return Invalid("particle ledger conservation sum overflow");
  const std::uint64_t left = candidate.activeStart + candidate.injected;
  const std::uint64_t right = candidate.activeEnd + candidate.escaped +
      candidate.absorbed + candidate.failed;
  if (left != right || candidate.advanced != candidate.activeEnd)
    return Core::Status(Core::StatusCode::Error,
                        "particle ledger does not close exactly or active count differs");
  candidate.closed = true;
  found->second = candidate;
  return Core::Status::OK();
}

const LedgerRow* ParticleLedger::Find(std::uint64_t step, int species) const {
  const auto found = rows_.find(LedgerKey{step, species});
  return found == rows_.end() ? nullptr : &found->second;
}

std::vector<LedgerRow> ParticleLedger::Rows() const {
  std::vector<LedgerRow> result;
  result.reserve(rows_.size());
  for (const auto& entry : rows_) result.push_back(entry.second);
  return result;
}

}  // namespace Adapters
}  // namespace SEP3D
