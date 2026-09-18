// ============================================================================
// R05 shock/source lifecycle and conservation accounting.
//
// Shock reconstruction is provider-owned; particle allocation is host-owned.
// This AMPS-independent layer joins those responsibilities with one immutable
// time-stamped shock state and one deterministic injection plan.  It never
// allocates PIC particles, so an invalid or capped plan can be rejected before
// any linked-list mutation occurs.
// ============================================================================

#ifndef SEP3D_ADAPTERS_SOURCE_RUNTIME_H
#define SEP3D_ADAPTERS_SOURCE_RUNTIME_H

#include "swcme_source_adapter.h"

#include <cstdint>
#include <string>
#include <vector>

namespace SEP3D {
namespace Adapters {

struct ShockState {
  Core::Status status;
  bool active = false;
  std::uint64_t generation = 0;
  double epochS = 0.0;
  double validUntilS = 0.0;
  Core::Vec3 centerM;
  double radiusM = 0.0;
  double radialSpeedMPerS = 0.0;
  double compressionRatio = 1.0;
  std::string providerIdentity;
  std::string configurationFingerprint;
  std::vector<ShockSourceRecord> patches;

  bool Covers(double timeS) const;
};

class ShockProvider {
 public:
  virtual ~ShockProvider() = default;
  virtual const char* CanonicalName() const = 0;
  virtual ShockState Evaluate(double timeS) const = 0;
};

// Deterministic spherical provider used by standalone verification and by
// production analytic-shock campaigns.  Canonical SWCME and host providers
// publish the same ShockState record through PublishedShockProvider.
class AnalyticSphericalShockProvider final : public ShockProvider {
 public:
  AnalyticSphericalShockProvider(
      double activeFromS, double activeUntilS, const Core::Vec3& centerM,
      double initialRadiusM, double speedMPerS, double compressionRatio,
      std::uint64_t generation, const std::string& fingerprint);
  const char* CanonicalName() const override { return "analytic-spherical"; }
  ShockState Evaluate(double timeS) const override;
 private:
  ShockState base_;
  double activeFromS_ = 0.0;
};

class PublishedShockProvider final : public ShockProvider {
 public:
  explicit PublishedShockProvider(std::string identity);
  Core::Status Publish(const ShockState& state);
  const char* CanonicalName() const override { return identity_.c_str(); }
  ShockState Evaluate(double timeS) const override;
 private:
  std::string identity_;
  ShockState published_;
  bool available_ = false;
};

struct SourceLedgerRow {
  std::uint64_t step = 0;
  int species = -1;
  std::uint64_t shockGeneration = 0;
  // Stable shock-patch identity makes canonical restart sorting independent
  // of MPI gather order when several patches inject on the same tick.
  std::uint64_t sourceId = 0;
  double representedParticles = 0.0;
  double injectedEnergyJ = 0.0;
  Core::Vec3 injectedMomentumKgMPerS;
  std::uint64_t macroparticles = 0;
  std::uint64_t rejected = 0;
  std::uint64_t capped = 0;
  std::uint64_t inactivePatches = 0;
  std::uint64_t disconnectedPatches = 0;
};

struct SourceRequest {
  ShockSourceRecord patch;
  std::uint64_t step = 0;
  int species = -1;
  double speciesMassKg = 0.0;
  double intervalS = 0.0;
  double physicalParticleRatePerS = 0.0;
  double macroparticleWeight = 0.0;
  std::uint64_t maximumMacroparticles = 0;
  bool connected = true;
};

struct InjectionPlan {
  Core::Status status;
  std::vector<InjectedParticle> particles;
  SourceLedgerRow ledger;
};

// Apply stochastic rounding with a semantic (campaign,event,patch,species,
// step) key, then sample a complete position/momentum/pitch record for every
// accepted macroparticle.  A cap is visible in the ledger and the weight is
// renormalized so represented physical number remains conservative.
InjectionPlan BuildInjectionPlan(const SourceRequest& request);

}  // namespace Adapters
}  // namespace SEP3D

#endif  // SEP3D_ADAPTERS_SOURCE_RUNTIME_H
