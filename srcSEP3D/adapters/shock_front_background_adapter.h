#ifndef SEP3D_ADAPTERS_SHOCK_FRONT_BACKGROUND_ADAPTER_H
#define SEP3D_ADAPTERS_SHOCK_FRONT_BACKGROUND_ADAPTER_H

#include "../background/bg_provider.h"
#include "../runtime/run_configuration.h"
#include "../core/computation_frame.h"

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
  // Computation frame of the run (core/computation_frame.h).  Evaluate()
  // returns frame quantities; the shared provider is always queried in HCI.
  // main_lib.cpp passes this frame to ReducedFrontShockSurface and uses it
  // to rotate the shock axis into the mover's ExpandingShock record.
  const Core::ComputationFrame& Frame() const { return frame_; }
  // HCI ambient primitive at a frame point, converted to frame vectors (B, U)
  // without derivatives.  Usable before Prepare() (tube tracing at mesh
  // construction); sector is the provider's open-field sector (0 = closed).
  //   framePositionM  point in mesh coordinates [m] (the adapter subtracts the
  //                   coordinate origin and rotates into HCI);
  //   epochS          model elapsed time [s] of the query;
  //   magneticFieldT  B [T], frame vector (may be null);
  //   velocityMPerS   U - Omega x x [m/s], frame velocity (may be null);
  //   sector          open-field sector of the ambient (may be null).
  // The value comes straight from the core's ambient model
  // (Provider::Ambient()): no shock state enters, as in Evaluate().
  // Callers: the observer flux-tube tracer (field at the frame's reference
  // epoch) and the tube-line source (B and U' at the release point for the
  // pitch-angle law and the field orientation of the injected particle).
  Core::Status AmbientAt(const Core::Vec3& framePositionM, double epochS,
                         Core::Vec3* magneticFieldT, Core::Vec3* velocityMPerS,
                         int* sector) const;

 private:
  std::shared_ptr<SEP::CoronaSwcme::ShockFront::Provider> provider_;
  Background::SnapshotMetadata metadata_;
  Core::Vec3 originM_;
  Core::ComputationFrame frame_;   // run's computation frame
  // mesh.domain.coordinate_frame of the run (e.g. HCI-corotating); recorded
  // as the coordinate frame of the published snapshot.
  std::string runFrameName_;
  std::uint64_t configurationDigest_ = 0;
  bool prepared_ = false;
};

} } // namespace SEP3D::Adapters

#endif
