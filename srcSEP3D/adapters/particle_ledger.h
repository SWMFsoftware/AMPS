// ============================================================================
// Exact integer particle accounting for Phase A.
//
// Each species/step row obeys
//   active_start + injected = active_end + escaped + absorbed + failed.
// Advanced particles are counted diagnostically and active_end is supplied by
// the owner after list migration.  Closing a row is transactional: a mismatch
// returns an error and leaves the row open for diagnosis.
// ============================================================================

#ifndef SEP3D_ADAPTERS_PARTICLE_LEDGER_H
#define SEP3D_ADAPTERS_PARTICLE_LEDGER_H

#include "transport_adapter.h"

#include <cstdint>
#include <map>
#include <vector>

namespace SEP3D {
namespace Adapters {

struct LedgerKey {
  std::uint64_t step = 0;
  int species = -1;
  bool operator<(const LedgerKey& other) const {
    return step < other.step || (step == other.step && species < other.species);
  }
};

struct LedgerRow {
  LedgerKey key;
  std::uint64_t activeStart = 0;
  std::uint64_t injected = 0;
  std::uint64_t advanced = 0;
  std::uint64_t escaped = 0;
  std::uint64_t absorbed = 0;
  std::uint64_t failed = 0;
  std::uint64_t shockCrossings = 0;
  std::uint64_t activeEnd = 0;
  bool closed = false;
};

class ParticleLedger final {
 public:
  // Remove the previous rank-local working set before a new AMPS particle
  // phase begins.  Production preserves already reduced, closed rows in its
  // separate history vector; this object is the mover's per-step scratchpad.
  void Clear();
  Core::Status Begin(std::uint64_t step, int species,
                     std::uint64_t activeStart);
  Core::Status RecordInjection(std::uint64_t step, int species,
                               std::uint64_t count = 1);
  Core::Status RecordMover(std::uint64_t step, int species,
                           ParticleDisposition disposition,
                           bool shockCrossed = false);
  Core::Status Close(std::uint64_t step, int species,
                     std::uint64_t activeEnd);
  // Install a row that has already been summed over every MPI rank.  The same
  // exact closure checks used by Close are applied before it becomes visible
  // to sampling or restart.  This avoids pretending rank migration is a local
  // source/sink while retaining a globally exact conservation identity.
  Core::Status ImportClosed(const LedgerRow& row);
  const LedgerRow* Find(std::uint64_t step, int species) const;
  std::vector<LedgerRow> Rows() const;

 private:
  std::map<LedgerKey, LedgerRow> rows_;
};

}  // namespace Adapters
}  // namespace SEP3D

#endif  // SEP3D_ADAPTERS_PARTICLE_LEDGER_H
