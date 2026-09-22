// Analytic Parker-spiral background provider used by standalone Phase-B runs.
// All parameters and outputs are SI.  Prepare() freezes only epoch metadata;
// the analytic state itself is immutable after construction.

#ifndef SEP3D_BG_PARKER_H
#define SEP3D_BG_PARKER_H

#include "../core/parker_geometry.h"
#include "bg_provider.h"
#include "swcme_solarwind.hpp"

namespace SEP3D {
namespace Background {

struct ParkerConfiguration {
  double sourceRadiusM = 2.5 * Core::Const::R_sun;
  double sourceLongitudeRad = 0.0;
  double sourceColatitudeRad = 0.5 * Core::Const::kPi;
  double referenceRadiusM = Core::Const::AU;
  double radialFieldAtReferenceT = 3.0e-9;
  double numberDensityAtReferenceM3 = 5.0e6;
  double densityReferenceRadiusM = Core::Const::AU;
  double temperatureK = 1.0e5;
  // The following values are the complete canonical SWCME ambient closure.
  // They are explicit here because pressure, mass density, sound speed, and
  // Alfvén speed are not determined by proton temperature alone when alpha
  // particles and electron pressure are enabled.
  double adiabaticIndex = 5.0 / 3.0;
  swcme::solarwind::ThermodynamicClosure thermodynamicClosure =
      swcme::solarwind::ThermodynamicClosure::ProtonOnly;
  double alphaToProtonRatio = 0.0;
  double electronTemperatureK = 1.0e5;
  double alphaTemperatureK = 1.0e5;
  double referenceSinColatitude = 1.0;
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

  const char* CanonicalName() const override {
    return "analytic-parker-swcme-v2";
  }
  Core::Status Validate() const override;
  Core::Status Prepare(double timeS) override;
  const SnapshotMetadata* PreparedMetadata() const override;
  BackgroundSample Evaluate(const Core::Vec3& positionM) const override;
  std::string ResolvedManifest() const override;
  ProviderCapabilities Capabilities() const override;

  const ParkerConfiguration& configuration() const { return configuration_; }

 private:
  ParkerConfiguration configuration_;
  // PreparedState is the canonical, header-only SWCME ambient cache. It owns
  // the Leblanc normalization, one-AU Parker conversion, and thermodynamic
  // closure used by every Evaluate() call. It is replaced only after Validate
  // succeeds, so a failed preparation cannot publish a partial solar wind.
  swcme::solarwind::PreparedState solarWind_;
  SnapshotMetadata metadata_;
  bool prepared_ = false;
  std::uint64_t configurationDigest_ = 0;
};

}  // namespace Background
}  // namespace SEP3D

#endif  // SEP3D_BG_PARKER_H
