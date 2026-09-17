// ============================================================================
// Immutable Phase-B background snapshots and two-slot publication buffer
//
// A builder evaluates a prepared provider entirely into temporary storage.
// Only after every field and metadata item validates is a shared immutable
// snapshot returned.  SnapshotBuffer applies the same validate-before-commit
// rule, so a failed coupled update cannot partially replace the active state.
// ============================================================================

#ifndef SEP3D_BACKGROUND_SNAPSHOT_H
#define SEP3D_BACKGROUND_SNAPSHOT_H

#include "bg_provider.h"

#include <memory>
#include <vector>

namespace SEP3D {
namespace Background {

Core::Status ValidateCompleteSample(const BackgroundSample& sample,
                                    const ProviderCapabilities& capabilities);

class BackgroundSnapshot final {
 public:
  BackgroundSnapshot(const SnapshotMetadata& metadata,
                     const ProviderCapabilities& capabilities,
                     const std::vector<Core::Vec3>& positions,
                     const std::vector<BackgroundSample>& samples);

  const SnapshotMetadata& metadata() const { return metadata_; }
  const ProviderCapabilities& capabilities() const { return capabilities_; }
  const std::vector<Core::Vec3>& positions() const { return positions_; }
  const std::vector<BackgroundSample>& samples() const { return samples_; }
  bool Covers(double timeS) const;

 private:
  const SnapshotMetadata metadata_;
  const ProviderCapabilities capabilities_;
  const std::vector<Core::Vec3> positions_;
  const std::vector<BackgroundSample> samples_;
};

class BackgroundSnapshotBuilder final {
 public:
  // On failure, *output is intentionally unchanged.  This differs from a
  // conventional factory because callers often pass the currently published
  // snapshot and require atomic candidate replacement.
  Core::Status Build(
      const BackgroundProvider& provider,
      const std::vector<Core::Vec3>& positions,
      std::shared_ptr<const BackgroundSnapshot>* output) const;
};

class SnapshotBuffer final {
 public:
  Core::Status PublishCurrent(
      const std::shared_ptr<const BackgroundSnapshot>& snapshot);
  Core::Status StageNext(
      const std::shared_ptr<const BackgroundSnapshot>& snapshot);

  const std::shared_ptr<const BackgroundSnapshot>& current() const {
    return current_;
  }
  const std::shared_ptr<const BackgroundSnapshot>& next() const {
    return next_;
  }

  // Returns the current snapshot at its epoch, the next snapshot at its
  // epoch, or a newly allocated linear interpolation in between.  No temporal
  // extrapolation is permitted.
  Core::Status SnapshotAt(
      double timeS, std::shared_ptr<const BackgroundSnapshot>* output) const;

 private:
  std::shared_ptr<const BackgroundSnapshot> current_;
  std::shared_ptr<const BackgroundSnapshot> next_;
};

}  // namespace Background
}  // namespace SEP3D

#endif  // SEP3D_BACKGROUND_SNAPSHOT_H
