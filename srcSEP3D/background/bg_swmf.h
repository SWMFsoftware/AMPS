// Phase-B SWMF/AWSoM ambient-field adapter.
//
// This is deliberately an import adapter, not an SWMF reader.  The coupled
// host supplies already received records and declares their units, frame,
// epoch, and ownership.  The adapter validates and converts a complete
// candidate before replacing its previously prepared state.

#ifndef SEP3D_BG_SWMF_H
#define SEP3D_BG_SWMF_H

#include "bg_provider.h"

#include <vector>

namespace SEP3D {
namespace Background {

enum class SwmfUnitSystem {
  SI,
  // Common coupling record: position [R_sun], B [nT], U [km/s], density
  // [cm^-3], pressure [nPa], grad(B) [nT/R_sun], curvature/div(b) [1/R_sun].
  AwsomCoupling
};

struct SwmfRawSample {
  Core::Vec3 position;
  Core::Vec3 magnetic;
  Core::Tensor3 magneticGradient;
  Core::Vec3 velocity;
  Core::Tensor3 velocityGradient;
  double numberDensity = 0.0;
  double temperature = 0.0;
  double pressure = 0.0;
  double velocityDivergence = 0.0;
  double divBhat = 0.0;
  double focusingLength = 0.0;
  Core::Vec3 curvature;
  double fieldAlignedStrain = 0.0;
  double epochS = 0.0;
  bool complete = false;
};

struct SwmfImport {
  SwmfUnitSystem units = SwmfUnitSystem::SI;
  std::string coordinateFrame = "HCI-like-inertial";
  std::string expectedCoordinateFrame = "HCI-like-inertial";
  std::string providerIdentity = "swmf-awsom";
  std::string configurationFingerprint;
  StorageOwnership ownership = StorageOwnership::ImportedReadOnly;
  double epochS = 0.0;
  double validUntilS = 0.0;
  std::uint64_t generation = 0;
  std::vector<SwmfRawSample> samples;
};

class SwmfAwsomProvider final : public BackgroundProvider {
 public:
  // Public only so the translation-unit conversion helper can construct a
  // candidate without mutating provider state; callers receive records solely
  // through BackgroundProvider and cannot access the private vector.
  struct Record {
    Core::Vec3 positionM;
    BackgroundSample sample;
  };

  Core::Status Load(const SwmfImport& candidate);

  const char* CanonicalName() const override { return "swmf-awsom-import-v1"; }
  Core::Status Validate() const override;
  Core::Status Prepare(double timeS) override;
  const SnapshotMetadata* PreparedMetadata() const override;
  BackgroundSample Evaluate(const Core::Vec3& positionM) const override;
  Core::Status EvaluateBatchDetailed(
      const double* xM, const double* yM, const double* zM,
      std::size_t count, BackgroundSample* output,
      Core::Status* perSampleStatus) const override;
  std::string ResolvedManifest() const override;
  ProviderCapabilities Capabilities() const override;

 private:
  SwmfImport import_;
  SnapshotMetadata metadata_;
  std::vector<Record> records_;
  bool loaded_ = false;
  bool prepared_ = false;
};

}  // namespace Background
}  // namespace SEP3D

#endif  // SEP3D_BG_SWMF_H
