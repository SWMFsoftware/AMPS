// SWCME runtime -> provider-neutral immutable mesh snapshots. No AMPS/MPI API.
// This adapter samples the canonical model; main_lib.cpp alone owns mesh
// offsets, cadence scheduling, collective validation and halo publication.
// See background/README.md for the physical units and preparation contract.
#ifndef SEP3D_BG_SWCME_H
#define SEP3D_BG_SWCME_H
#include "bg_provider.h"
#include <memory>
namespace swcme { namespace input3d { struct ResolvedConfiguration; } }
namespace SEP3D { namespace Background {
struct ParkerConfiguration;
// One instance owns a frozen model configuration and one prepared epoch.
// Call Prepare at a joined update boundary, then reuse Evaluate/Batch for that
// epoch. Evaluation is read-only; Prepare must not run concurrently with it.
// The PIMPL keeps canonical SWCME headers out of public AMPS include chains.
class SwcmeBackgroundProvider final : public BackgroundProvider {
 public:
  // Sphere or finite-SSE model is already canonically resolved. innerAmbient is its matching
  // Parker/Leblanc continuation below SWCME's 1.05-Rs domain minimum. dt and
  // cadence are seconds; originM translates global positions to model-local
  // heliocentric coordinates. No cell storage is allocated by this constructor.
  SwcmeBackgroundProvider(const swcme::input3d::ResolvedConfiguration& model,
      const ParkerConfiguration& innerAmbient, double timeStepS,
      double cadenceS, const Core::Vec3& originM);
  ~SwcmeBackgroundProvider() override;
  const char* CanonicalName() const override { return "swcme-runtime-mesh-v1"; }
  Core::Status Validate() const override;
  // Build a complete candidate state at application time timeS, subtracting
  // the event launch epoch for SWCME. Replace prepared state only on success;
  // reject expired coverage, backwards generations and CME/inner-shell overlap.
  Core::Status Prepare(double timeS) override;
  // Non-owning view, null until successful Prepare and invalidated by the next
  // successful Prepare. SnapshotBuilder copies it into immutable storage.
  const SnapshotMetadata* PreparedMetadata() const override;
  // SI primitives plus derivatives from the actual regional vectors. An
  // invalid single query returns a sample whose status carries the failure.
  BackgroundSample Evaluate(const Core::Vec3& positionM) const override;
  // Parallel x/y/z arrays and caller-owned outputs, each of length count.
  // Successful points are written; rejected points retain their old output
  // bytes and receive individual statuses. The aggregate is the first error.
  // A zero-sized batch is legal after Prepare, including on an empty MPI rank.
  Core::Status EvaluateBatchDetailed(const double*,const double*,const double*,
      std::size_t,BackgroundSample*,Core::Status*) const override;
  std::string ResolvedManifest() const override;
  // Availability flags describe numerical, not analytic, vector derivatives.
  ProviderCapabilities Capabilities() const override;
 private:
  struct Implementation;
  std::unique_ptr<Implementation> implementation_;
};
} }
#endif
