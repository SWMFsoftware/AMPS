#ifndef SEP_UTIL_BACKGROUND_SNAPSHOT_H
#define SEP_UTIL_BACKGROUND_SNAPSHOT_H

#include <atomic>
#include <cstdint>
#include <memory>
#include <mutex>
#include <string>

namespace SEP {
namespace Background {

// Identifies the authority that supplied the solar-wind, IMF, shock, and
// field-line state consumed by a particle step.  LocalEvolution is deliberately
// distinct from SWMF: once imported SWMF data are copied for local evolution,
// the copy has a new owner and must never be presented as still SWMF-owned.
enum class Provider {
  Analytic,
  Swcme,
  Swmf,
  LocalEvolution
};

// Describes who is permitted to modify the storage represented by a snapshot.
// ImportedReadOnly is mandatory for SWMF imports.  HandoffCopy records the only
// supported transition from an imported state to state evolved by srcSEP.
enum class Ownership {
  ModelOwned,
  ImportedReadOnly,
  HandoffCopy
};

const char* ProviderName(Provider provider);
const char* OwnershipName(Ownership ownership);

// Immutable metadata binding one particle step to one background realization.
// All times are seconds from the authoritative PIC simulation-time origin.  The
// validity interval is closed because a coupling state may be sampled exactly
// at either reported endpoint.  The object has no assignment operator and is
// only published as shared_ptr<const BackgroundSnapshot>, preventing a mover
// from changing provenance while it is using the associated field-line data.
class BackgroundSnapshot {
 public:
  BackgroundSnapshot(Provider provider,
                     Ownership ownership,
                     double epoch_seconds,
                     double valid_from_seconds,
                     double valid_until_seconds,
                     std::uint64_t field_line_generation,
                     const std::string& configuration_fingerprint,
                     const std::string& provenance,
                     double previous_physical_epoch_seconds = 0.0,
                     double current_physical_epoch_seconds = 0.0);

  BackgroundSnapshot(const BackgroundSnapshot&) = default;
  BackgroundSnapshot& operator=(const BackgroundSnapshot&) = delete;

  Provider provider() const { return provider_; }
  Ownership ownership() const { return ownership_; }
  double epoch_seconds() const { return epoch_seconds_; }
  double valid_from_seconds() const { return valid_from_seconds_; }
  double valid_until_seconds() const { return valid_until_seconds_; }
  std::uint64_t field_line_generation() const {
    return field_line_generation_;
  }
  const std::string& configuration_fingerprint() const {
    return configuration_fingerprint_;
  }
  const std::string& provenance() const { return provenance_; }
  double previous_physical_epoch_seconds() const {
    return previous_physical_epoch_seconds_;
  }
  double current_physical_epoch_seconds() const {
    return current_physical_epoch_seconds_;
  }
  double physical_epoch_interval_seconds() const {
    return current_physical_epoch_seconds_ - previous_physical_epoch_seconds_;
  }

  bool Covers(double simulation_time_seconds) const;

 private:
  const Provider provider_;
  const Ownership ownership_;
  const double epoch_seconds_;
  const double valid_from_seconds_;
  const double valid_until_seconds_;
  const std::uint64_t field_line_generation_;
  const std::string configuration_fingerprint_;
  const std::string provenance_;
  // These epochs identify the two physical background realizations used for
  // time derivatives.  They never inherit a particle substep duration.
  const double previous_physical_epoch_seconds_;
  const double current_physical_epoch_seconds_;
};

class SnapshotStore;

// RAII boundary for the part of a global step in which particle movers may
// read background data.  Publication is rejected while a read phase is active,
// so the AMPS vertex arrays named by the snapshot cannot be replaced halfway
// through a threaded particle step.  The guard is movable to support ordinary
// C++11 return-by-value, but is intentionally not copyable.
class ParticleReadPhase {
 public:
  ParticleReadPhase(ParticleReadPhase&& other) noexcept;
  ParticleReadPhase& operator=(ParticleReadPhase&& other) noexcept;
  ~ParticleReadPhase();

  ParticleReadPhase(const ParticleReadPhase&) = delete;
  ParticleReadPhase& operator=(const ParticleReadPhase&) = delete;

  const BackgroundSnapshot& snapshot() const { return *snapshot_; }

 private:
  friend class SnapshotStore;
  ParticleReadPhase(SnapshotStore* store,
                    const std::shared_ptr<const BackgroundSnapshot>& snapshot);
  void Release();

  SnapshotStore* store_;
  std::shared_ptr<const BackgroundSnapshot> snapshot_;
};

// Process-wide publication point for background metadata.  The store does not
// own or copy the large AMPS field-line arrays; it freezes their identity and
// forbids a writer from publishing a replacement during PIC::TimeStep().  That
// narrow contract avoids a per-particle deep copy while still making provider,
// epoch, generation, validity, and configuration explicit and testable.
class SnapshotStore {
 public:
  static SnapshotStore& Instance();

  void Publish(const BackgroundSnapshot& snapshot);

  // A provider change is never accepted by Publish().  The separate handoff
  // operation makes the transfer auditable and verifies the expected previous
  // provider.  In particular, SWMF -> LocalEvolution requires HandoffCopy.
  void PublishHandoff(const BackgroundSnapshot& snapshot,
                      Provider expected_previous_provider);

  // Background implementations call this before changing their backing arrays
  // or prepared caches.  It closes the gap between immutable metadata and the
  // legacy AMPS/SWCME storage by rejecting writes during particle motion and by
  // preventing one provider from mutating another provider's authoritative
  // realization before publication is attempted.
  void AssertProviderMayWrite(Provider provider) const;

  ParticleReadPhase BeginParticleRead(double simulation_time_seconds);

  // Every configured mover calls this inside an active ParticleReadPhase.  The
  // enclosing guard owns the metadata and publication is forbidden, so a const
  // reference remains valid through all legacy early-return paths.  Avoiding a
  // mutex and shared-pointer reference-count operation per particle keeps this
  // safety boundary out of the mover's inner-loop performance profile.
  const BackgroundSnapshot& AcquireForMover() const;

  std::shared_ptr<const BackgroundSnapshot> Current() const;
  bool HasSnapshot() const;

  // Test-only reset.  Production code must never erase provenance mid-run.
  void ResetForTests();

 private:
  friend class ParticleReadPhase;
  SnapshotStore();

  void PublishLocked(const BackgroundSnapshot& snapshot, bool handoff,
                     Provider expected_previous_provider);
  void EndParticleRead();

  mutable std::mutex mutex_;
  std::shared_ptr<const BackgroundSnapshot> current_;
  std::atomic<unsigned int> active_read_phases_;
};

// Deterministic FNV-1a fingerprint for a canonical configuration string.  This
// is an identity/checking hash, not a cryptographic signature.  Providers must
// build the canonical string without mover selection so changing the transport
// algorithm cannot silently change the physical background.
std::string FingerprintConfiguration(const std::string& canonical_configuration);

}  // namespace Background
}  // namespace SEP

#endif  // SEP_UTIL_BACKGROUND_SNAPSHOT_H
