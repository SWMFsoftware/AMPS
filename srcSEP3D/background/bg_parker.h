// Analytic Parker-spiral background provider used by standalone Phase-B runs.
// All parameters and outputs are SI.  Prepare() freezes only epoch metadata;
// the analytic state itself is immutable after construction.

#ifndef SEP3D_BG_PARKER_H
#define SEP3D_BG_PARKER_H

#include "../core/parker_geometry.h"
#include "bg_provider.h"

namespace SEP3D {
namespace Background {

struct ParkerConfiguration {
  double sourceRadiusM = 2.5 * Core::Const::R_sun;
  double sourceLongitudeRad = 0.0;
  double sourceColatitudeRad = 0.5 * Core::Const::kPi;
  double referenceRadiusM = Core::Const::AU;
  double radialFieldAtReferenceT = 3.0e-9;
  double numberDensityAtReferenceM3 = 5.0e6;
  double temperatureK = 1.0e5;
  double solarWindSpeedMPerS = Core::Const::V_sw_default;
  double solarRotationRateRadPerS = Core::Const::Omega_sun;
  Core::Vec3 rotationAxis = {0.0, 0.0, 1.0};
  int magneticPolarity = 1;
  double validityCadenceS = 3600.0;
  std::string coordinateFrame = "HCI-like-inertial";
};

class AnalyticParkerProvider final : public BackgroundProvider {
 public:
  explicit AnalyticParkerProvider(const ParkerConfiguration& configuration);

  const char* CanonicalName() const override { return "analytic-parker-v1"; }
  Core::Status Validate() const override;
  Core::Status Prepare(double timeS) override;
  const SnapshotMetadata* PreparedMetadata() const override;
  BackgroundSample Evaluate(const Core::Vec3& positionM) const override;
  std::string ResolvedManifest() const override;
  ProviderCapabilities Capabilities() const override;

  const ParkerConfiguration& configuration() const { return configuration_; }

 private:
  ParkerConfiguration configuration_;
  SnapshotMetadata metadata_;
  bool prepared_ = false;
  std::uint64_t configurationDigest_ = 0;
};

}  // namespace Background
}  // namespace SEP3D

#endif  // SEP3D_BG_PARKER_H
