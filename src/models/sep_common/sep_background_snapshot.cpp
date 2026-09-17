#include "sep_background_snapshot.h"

#include <cmath>
#include <iomanip>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <utility>

namespace SEP {
namespace Background {

namespace {

void ValidateOwnership(Provider provider, Ownership ownership) {
  if (provider == Provider::Swmf && ownership != Ownership::ImportedReadOnly) {
    throw std::invalid_argument(
        "an SWMF background snapshot must be imported read-only");
  }

  if (provider == Provider::LocalEvolution &&
      ownership != Ownership::HandoffCopy) {
    throw std::invalid_argument(
        "a locally evolved imported background must use handoff-copy ownership");
  }

  if ((provider == Provider::Analytic || provider == Provider::Swcme) &&
      ownership != Ownership::ModelOwned) {
    throw std::invalid_argument(
        "analytic and SWCME backgrounds must be model-owned");
  }
}

}  // namespace

const char* ProviderName(Provider provider) {
  switch (provider) {
    case Provider::Analytic: return "analytic";
    case Provider::Swcme: return "swcme";
    case Provider::Swmf: return "swmf";
    case Provider::LocalEvolution: return "local-evolution";
  }

  return "unknown";
}

const char* OwnershipName(Ownership ownership) {
  switch (ownership) {
    case Ownership::ModelOwned: return "model-owned";
    case Ownership::ImportedReadOnly: return "imported-read-only";
    case Ownership::HandoffCopy: return "handoff-copy";
  }

  return "unknown";
}

BackgroundSnapshot::BackgroundSnapshot(
    Provider provider, Ownership ownership, double epoch_seconds,
    double valid_from_seconds, double valid_until_seconds,
    std::uint64_t field_line_generation,
    const std::string& configuration_fingerprint,
    const std::string& provenance,
    double previous_physical_epoch_seconds,
    double current_physical_epoch_seconds)
    : provider_(provider),
      ownership_(ownership),
      epoch_seconds_(epoch_seconds),
      valid_from_seconds_(valid_from_seconds),
      valid_until_seconds_(valid_until_seconds),
      field_line_generation_(field_line_generation),
      configuration_fingerprint_(configuration_fingerprint),
      provenance_(provenance),
      previous_physical_epoch_seconds_(previous_physical_epoch_seconds),
      current_physical_epoch_seconds_(current_physical_epoch_seconds) {
  ValidateOwnership(provider_, ownership_);

  // Epoch and lower validity bounds must be finite.  A positive infinity upper
  // bound is permitted for a static analytic field or the most recent imported
  // coupling state, which remains valid until its provider publishes a newer
  // generation.
  if (!std::isfinite(epoch_seconds_) || !std::isfinite(valid_from_seconds_) ||
      std::isnan(valid_until_seconds_)) {
    throw std::invalid_argument("background snapshot contains a non-finite time");
  }

  if (valid_until_seconds_ < valid_from_seconds_ ||
      epoch_seconds_ < valid_from_seconds_ ||
      epoch_seconds_ > valid_until_seconds_) {
    throw std::invalid_argument(
        "background epoch must lie inside an ordered validity interval");
  }

  // The default (0,0) pair represents a static/frozen realization and keeps
  // the constructor source-compatible with older checkpoints.  A provider
  // publishing time-dependent data must supply finite, ordered physical epochs.
  if (!std::isfinite(previous_physical_epoch_seconds_) ||
      !std::isfinite(current_physical_epoch_seconds_) ||
      current_physical_epoch_seconds_ < previous_physical_epoch_seconds_) {
    throw std::invalid_argument("background physical epochs are invalid");
  }

  if (field_line_generation_ == 0) {
    throw std::invalid_argument("field-line generation zero is reserved");
  }

  if (configuration_fingerprint_.empty()) {
    throw std::invalid_argument("configuration fingerprint must not be empty");
  }

  if (provenance_.empty()) {
    throw std::invalid_argument("background provenance must not be empty");
  }
}

bool BackgroundSnapshot::Covers(double simulation_time_seconds) const {
  return std::isfinite(simulation_time_seconds) &&
         simulation_time_seconds >= valid_from_seconds_ &&
         simulation_time_seconds <= valid_until_seconds_;
}

ParticleReadPhase::ParticleReadPhase(
    SnapshotStore* store,
    const std::shared_ptr<const BackgroundSnapshot>& snapshot)
    : store_(store), snapshot_(snapshot) {}

ParticleReadPhase::ParticleReadPhase(ParticleReadPhase&& other) noexcept
    : store_(other.store_), snapshot_(std::move(other.snapshot_)) {
  other.store_ = NULL;
}

ParticleReadPhase& ParticleReadPhase::operator=(
    ParticleReadPhase&& other) noexcept {
  if (this != &other) {
    Release();
    store_ = other.store_;
    snapshot_ = std::move(other.snapshot_);
    other.store_ = NULL;
  }

  return *this;
}

ParticleReadPhase::~ParticleReadPhase() { Release(); }

void ParticleReadPhase::Release() {
  if (store_ != NULL) {
    store_->EndParticleRead();
    store_ = NULL;
    snapshot_.reset();
  }
}

SnapshotStore::SnapshotStore() : active_read_phases_(0) {}

SnapshotStore& SnapshotStore::Instance() {
  static SnapshotStore store;
  return store;
}

void SnapshotStore::Publish(const BackgroundSnapshot& snapshot) {
  std::lock_guard<std::mutex> lock(mutex_);
  PublishLocked(snapshot, false, snapshot.provider());
}

void SnapshotStore::PublishHandoff(
    const BackgroundSnapshot& snapshot,
    Provider expected_previous_provider) {
  std::lock_guard<std::mutex> lock(mutex_);
  PublishLocked(snapshot, true, expected_previous_provider);
}

void SnapshotStore::AssertProviderMayWrite(Provider provider) const {
  std::lock_guard<std::mutex> lock(mutex_);

  if (active_read_phases_.load(std::memory_order_acquire) != 0) {
    throw std::logic_error(
        "background provider attempted to write during particle motion");
  }

  if (current_ && current_->provider() != provider) {
    throw std::logic_error(
        "background provider attempted to overwrite another provider's state");
  }
}

void SnapshotStore::PublishLocked(
    const BackgroundSnapshot& snapshot, bool handoff,
    Provider expected_previous_provider) {
  if (active_read_phases_.load(std::memory_order_acquire) != 0) {
    throw std::logic_error(
        "cannot publish background state during an active particle read phase");
  }

  if (current_) {
    if (snapshot.field_line_generation() <
        current_->field_line_generation()) {
      throw std::logic_error("field-line generation must not decrease");
    }

    if (snapshot.epoch_seconds() < current_->epoch_seconds()) {
      throw std::logic_error("background epoch must not move backwards");
    }

    const bool provider_changed = snapshot.provider() != current_->provider();
    if (provider_changed && !handoff) {
      throw std::logic_error(
          "background provider changes require an explicit handoff");
    }

    if (!provider_changed && handoff) {
      throw std::logic_error("handoff requested without a provider change");
    }

    if (handoff) {
      if (current_->provider() != expected_previous_provider) {
        throw std::logic_error(
            "background handoff did not match the expected previous provider");
      }

      if (current_->provider() == Provider::Swmf &&
          (snapshot.provider() != Provider::LocalEvolution ||
           snapshot.ownership() != Ownership::HandoffCopy)) {
        throw std::logic_error(
            "SWMF state may be handed off only to a local-evolution copy");
      }

      if (snapshot.field_line_generation() <=
          current_->field_line_generation()) {
        throw std::logic_error(
            "a provider handoff must publish a new field-line generation");
      }
    } else {
      if (snapshot.configuration_fingerprint() !=
          current_->configuration_fingerprint()) {
        throw std::logic_error(
            "configuration changed without an explicit provider handoff");
      }

      // A same-epoch replacement is ambiguous and often means that two owners
      // attempted to publish the same step.  A field-line generation increase
      // is the sole exception because SWMF can replace geometry/state at an
      // unchanged timestamp while retaining an explicit new generation.
      if (snapshot.epoch_seconds() == current_->epoch_seconds() &&
          snapshot.field_line_generation() ==
              current_->field_line_generation()) {
        throw std::logic_error(
            "duplicate background publication at the same epoch and generation");
      }
    }
  } else if (handoff) {
    throw std::logic_error("cannot hand off before an initial publication");
  }

  current_ = std::make_shared<const BackgroundSnapshot>(snapshot);
}

ParticleReadPhase SnapshotStore::BeginParticleRead(
    double simulation_time_seconds) {
  std::lock_guard<std::mutex> lock(mutex_);

  if (!current_) {
    throw std::logic_error(
        "particle step requested before a background snapshot was published");
  }

  if (!current_->Covers(simulation_time_seconds)) {
    std::ostringstream message;
    message << "background snapshot at epoch " << current_->epoch_seconds()
            << " s does not cover particle-step time "
            << simulation_time_seconds << " s";
    throw std::logic_error(message.str());
  }

  // Release publication of the active phase after current_ and its metadata
  // have been read under the mutex.  Worker threads acquire this flag before
  // dereferencing current_, establishing visibility without a per-mover lock.
  active_read_phases_.fetch_add(1, std::memory_order_release);
  return ParticleReadPhase(this, current_);
}

const BackgroundSnapshot& SnapshotStore::AcquireForMover() const {
  if (active_read_phases_.load(std::memory_order_acquire) == 0 || !current_) {
    throw std::logic_error(
        "particle mover entered outside the immutable background read phase");
  }

  // current_ cannot be replaced while active_read_phases_ is nonzero.  Its
  // lifetime is additionally retained by ParticleReadPhase::snapshot_, so this
  // reference is stable until PIC::TimeStep() and all worker movers return.
  return *current_;
}

std::shared_ptr<const BackgroundSnapshot> SnapshotStore::Current() const {
  std::lock_guard<std::mutex> lock(mutex_);
  return current_;
}

bool SnapshotStore::HasSnapshot() const {
  std::lock_guard<std::mutex> lock(mutex_);
  return static_cast<bool>(current_);
}

void SnapshotStore::EndParticleRead() {
  std::lock_guard<std::mutex> lock(mutex_);

  const unsigned int readers =
      active_read_phases_.load(std::memory_order_relaxed);
  if (readers == 0) {
    throw std::logic_error("unbalanced background particle-read phase");
  }

  active_read_phases_.store(readers - 1, std::memory_order_release);
}

void SnapshotStore::ResetForTests() {
  std::lock_guard<std::mutex> lock(mutex_);

  if (active_read_phases_.load(std::memory_order_acquire) != 0) {
    throw std::logic_error(
        "cannot reset background store during an active particle read phase");
  }

  current_.reset();
}

std::string FingerprintConfiguration(
    const std::string& canonical_configuration) {
  // 64-bit FNV-1a is stable across compilers and MPI ranks because it processes
  // the canonical byte string explicitly; std::hash is intentionally avoided
  // because the C++ standard does not promise stable values across executions.
  std::uint64_t value = UINT64_C(14695981039346656037);
  const std::uint64_t prime = UINT64_C(1099511628211);

  for (std::string::const_iterator it = canonical_configuration.begin();
       it != canonical_configuration.end(); ++it) {
    value ^= static_cast<unsigned char>(*it);
    value *= prime;
  }

  std::ostringstream out;
  out << std::hex << std::setfill('0') << std::setw(16) << value;
  return out.str();
}

}  // namespace Background
}  // namespace SEP
