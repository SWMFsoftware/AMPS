// ============================================================================
// srcSEP3D immutable run configuration
//
// This is an AMPS-independent R2 contract.  Process arguments and
// AMPS_PARAM.in are translated by a host before this object is constructed;
// neither Runtime nor the coupled entry points parse configuration text.
// ============================================================================

#ifndef SEP3D_RUNTIME_RUN_CONFIGURATION_H
#define SEP3D_RUNTIME_RUN_CONFIGURATION_H

#include "../core/sep3d_types.h"

#include <cstddef>
#include <cstdint>
#include <memory>
#include <string>

namespace SEP3D {
namespace RuntimeModel {

enum class BackgroundAuthority { AnalyticParker, Swmf };
enum class TurbulenceAuthority { Prescribed, Swmf };
enum class ShockAuthority { None, Swcme };
enum class TransportModel { Parker3D, Focused3D };
enum class DomainPreset { Earth, Mars };

const char* Name(BackgroundAuthority value);
const char* Name(TurbulenceAuthority value);
const char* Name(ShockAuthority value);
const char* Name(TransportModel value);
const char* Name(DomainPreset value);

// Mutable input record used only while the host resolves configuration.  The
// successful factory copies it into a RunConfiguration3D exposed solely
// through const accessors.  Output formatting fields are deliberately present
// so CFG3D03/LIFE3D03 can prove that they do not alter the physics fingerprint.
struct RunConfiguration3DOptions {
  BackgroundAuthority background = BackgroundAuthority::AnalyticParker;
  TurbulenceAuthority turbulence = TurbulenceAuthority::Prescribed;
  ShockAuthority shock = ShockAuthority::Swcme;
  TransportModel transport = TransportModel::Parker3D;
  DomainPreset domain = DomainPreset::Earth;

  double innerRadiusM = 20.0 * Core::Const::R_sun;
  double outerRadiusM = Core::Const::AU;
  double requestedTimeStepS = 1.0;
  std::uint64_t campaignSeed = 1;
  std::uint64_t backgroundCadenceSteps = 1;

  // Frozen pre-mesh storage choices.  Offsets are derived by the factory in a
  // canonical order; adapters may not append fields after Configure().
  bool storeMagneticGradient = false;
  bool storeVelocityGradient = false;
  std::size_t samplingBytesPerCell = 0;

  // Output-only controls: these affect products, not particle trajectories.
  std::uint64_t outputCadenceSteps = 1;
  std::string outputDirectory = "output";
  std::string outputPrefix = "sep3d";

  // Reserved physics must fail during configuration rather than quietly
  // becoming a no-op in a mover.
  bool enablePerpendicularDiffusion = false;
  bool enableDrifts = false;
  bool enableExternalScriptBackground = false;
  bool enableSelfConsistent3DTurbulence = false;
};

constexpr std::size_t kNoOffset = static_cast<std::size_t>(-1);

struct StorageLayout {
  std::size_t magneticFieldOffset = kNoOffset;       // 3 doubles
  std::size_t bulkVelocityOffset = kNoOffset;        // 3 doubles
  std::size_t numberDensityOffset = kNoOffset;       // 1 double
  std::size_t velocityDivergenceOffset = kNoOffset;  // 1 double
  std::size_t magneticGradientOffset = kNoOffset;    // optional 9 doubles
  std::size_t velocityGradientOffset = kNoOffset;    // optional 9 doubles
  std::size_t waveEnergyOffset = kNoOffset;          // optional 2 doubles
  std::size_t cellAssociatedBytes = 0;
  std::size_t samplingBytesPerCell = 0;
  std::string fingerprint;
};

bool operator==(const StorageLayout& left, const StorageLayout& right);
bool operator!=(const StorageLayout& left, const StorageLayout& right);

class RunConfiguration3D final {
 public:
  static Core::Status Create(
      const RunConfiguration3DOptions& options,
      std::shared_ptr<const RunConfiguration3D>* configuration);

  RunConfiguration3D(const RunConfiguration3D&) = default;
  RunConfiguration3D& operator=(const RunConfiguration3D&) = delete;

  const RunConfiguration3DOptions& options() const { return options_; }
  const StorageLayout& storage_layout() const { return storageLayout_; }
  const std::string& physics_fingerprint() const { return physicsFingerprint_; }
  const std::string& resolved_manifest() const { return resolvedManifest_; }

 private:
  RunConfiguration3D(const RunConfiguration3DOptions& options,
                     const StorageLayout& layout,
                     const std::string& physicsFingerprint,
                     const std::string& resolvedManifest);

  const RunConfiguration3DOptions options_;
  const StorageLayout storageLayout_;
  const std::string physicsFingerprint_;
  const std::string resolvedManifest_;
};

}  // namespace RuntimeModel
}  // namespace SEP3D

#endif  // SEP3D_RUNTIME_RUN_CONFIGURATION_H
