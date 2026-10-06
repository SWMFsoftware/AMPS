#ifndef SEP3D_ADAPTERS_SHOCK_FRONT_BACKGROUND_ADAPTER_H
#define SEP3D_ADAPTERS_SHOCK_FRONT_BACKGROUND_ADAPTER_H

#include "../background/bg_provider.h"
#include "../runtime/run_configuration.h"

#include <memory>
#include <string>

namespace SEP { namespace CoronaSwcme { namespace ShockFront {
class Provider;
struct Epoch;
} } }

namespace SEP3D { namespace Adapters {

// Thin SEP3D boundary for the shared reduced provider.  It converts neutral SI
// vectors/status into srcSEP3D's provider ABI; all front, ambient and RH
// equations remain under src/models.  Native owner/ghost publication therefore
// reuses the same BackgroundProvider/DATAFILE path as every maintained runtime
// model and does not add a PIC-core coupling mode.
class ShockFrontBackgroundAdapter final : public Background::BackgroundProvider {
 public:
  static Core::Status Create(const RuntimeModel::RunConfiguration3D& configuration,
      std::shared_ptr<Background::BackgroundProvider>* output);

  const char* CanonicalName() const override {
    return "sep-corona-swcme-shock-front-ambient-v1";
  }
  Core::Status Validate() const override;
  Core::Status Prepare(double timeS) override;
  const Background::SnapshotMetadata* PreparedMetadata() const override;
  Core::Status RestorePreparedGeneration(std::uint64_t generation) override;
  Background::BackgroundSample Evaluate(const Core::Vec3& positionM) const override;
  std::string ResolvedManifest() const override;
  Background::ProviderCapabilities Capabilities() const override;

  std::shared_ptr<const SEP::CoronaSwcme::ShockFront::Epoch> FrontEpoch() const;
  // Read-only access for native evidence/observer receipts.  Physics remains
  // owned by the shared provider; the application must never call Prepare()
  // through this handle or construct a second event authority.
  std::shared_ptr<const SEP::CoronaSwcme::ShockFront::Provider>
      SharedProvider() const { return provider_; }

 private:
  std::shared_ptr<SEP::CoronaSwcme::ShockFront::Provider> provider_;
  Background::SnapshotMetadata metadata_;
  Core::Vec3 originM_;
  std::uint64_t configurationDigest_ = 0;
  bool prepared_ = false;
};

} } // namespace SEP3D::Adapters

#endif
