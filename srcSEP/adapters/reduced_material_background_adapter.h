#ifndef SRCSEP_REDUCED_MATERIAL_BACKGROUND_ADAPTER_H
#define SRCSEP_REDUCED_MATERIAL_BACKGROUND_ADAPTER_H
#include "../util/sep_reduced_geometry.h"
namespace SEP { namespace ReducedShock { namespace MaterialAdapter {
namespace M=Background::Material;

// This owner is a value capture of the existing reduced interval and explicit
// numerical controls. It has no dependency on mutable PIC clocks, native
// arrays or a future configured event. Querying it cannot prepare/publish a
// core generation. Model equations remain in the canonical provider.
class Owner final : public M::Owner {
 public:
  Owner(IntervalView view,M::Descriptor descriptor,GeometryControls controls)
      :view_(std::move(view)),descriptor_(std::move(descriptor)),controls_(controls){}
  const M::Descriptor& Identity() const override {return descriptor_;}
  Core::Result<M::State> Sample(const M::Query& query) const override {
    return view_.MaterialStateAt(descriptor_,query);
  }
  Core::Result<M::Motion> Move(const Position& position,double begin,double end) const override {
    using R=Core::Result<M::Motion>;
    // The neutral descriptor identifies this EXACT captured core realization.
    // Validate it before evaluating a trajectory, just as MaterialStateAt does
    // for a sampled state. A manually constructed adapter cannot relabel one
    // event/composition/EOS as another merely by presenting valid metadata.
    const auto actual=view_.MaterialDescriptor(descriptor_.geometryGeneration,
        descriptor_.snapshotGeneration,descriptor_.nativeConfigurationIdentity);
    if(!actual.ok())return R::Failure(actual.status.code,actual.status.message);
    if(M::Payload(actual.value)!=M::Payload(descriptor_))
      return R::Failure(Core::StatusCode::DataIntegrityFailure,
          "material trajectory descriptor differs from the captured core owner");

    // MoveVertex is a generic core diagnostic and can trace the PFSS corona.
    // This NEW neutral profile is qualified only in the Parker exterior: a
    // finite core velocity in an open/closed coronal region is not proof of
    // its material/induction closure. Admit the canonical primitive at EVERY
    // RK query, including reverse-time stages; endpoint checks alone must not
    // permit an intermediate excursion through an unqualified region.
    //
    // Keep the existing IntegrateTrajectory controller, exact endpoints,
    // spatial error allocation, work ceiling and supplied numerical controls.
    // The callback returns the same core SI velocity without changing it.
    // No new trajectory equation, coronal fallback or physical AMPS timestep
    // is introduced, and legacy MoveVertex retains its diagnostic contract.
    const auto velocity=[&](const Position& point,double epoch) {
      const auto state=view_.PrimitiveAt(point,epoch);
      if(!state.ok())return Core::Result<Position>::Failure(state.status.code,state.status.message);
      if(!state.value.parkerExterior)
        return Core::Result<Position>::Failure(Core::StatusCode::UnsupportedCapability,
            "PFSS trajectory is not qualified as material Parker exterior evolution");
      return Core::Result<Position>::Success(state.value.velocityMPerS);
    };
    // Zero duration performs no RK evaluations. Validate its actual initial
    // point explicitly so it cannot grant coronal material authority either.
    const auto initial=velocity(position,begin);
    if(!initial.ok())return R::Failure(initial.status.code,initial.status.message);
    const auto moved=IntegrateTrajectory(velocity,position,begin,end,controls_);
    if(!moved.ok())return R::Failure(moved.status.code,moved.status.message);
    // The accepted RK combination can differ from its last stage position.
    // Validate that final position and export its actual HCI velocity [m/s],
    // rather than assuming the final stage's region or substituting a value.
    const auto finalVelocity=velocity(moved.value,end);
    if(!finalVelocity.ok())return R::Failure(finalVelocity.status.code,finalVelocity.status.message);
    M::Motion result;result.positionM=moved.value;result.velocityMPerS=finalVelocity.value;
    return R::Success(std::move(result));
  }
  Core::Result<M::Surface> Shock(const Position& position,double epoch) const override {
    using R=Core::Result<M::Surface>;
    const auto front=view_.FrontAt(position,epoch);
    if(!front.ok())return R::Failure(front.status.code,front.status.message);
    const auto& f=front.value;M::Surface result;
    result.stableIdentity=descriptor_.eventIdentity+":leading-SSE-interface";
    result.eventIdentity=f.eventIdentity;result.classification=f.classification;
    result.epochS=epoch;result.positionM=position;result.outwardNormal=f.outwardNormal;
    result.normalSpeedMPerS=f.normalSpeedMPerS;
    if(f.upstreamValid) {
      M::Query query;query.positionM=position;query.epochS=epoch;query.region=M::Region::ParkerExterior;
      const auto upstream=Sample(query);if(!upstream.ok())return R::Failure(upstream.status.code,upstream.status.message);
      result.upstream=upstream.value;result.upstream->side=M::Side::UpstreamLimit;
      result.upstream->interfaceIdentity=result.stableIdentity;
    }
    if(f.downstreamValid) {
      // Only the actual canonical RH primitive is available here. Species
      // temperatures, gradients and a downstream VOLUME remain absent; they
      // cannot be reconstructed from a compression ratio without a heating law.
      M::State down;down.eventIdentity=down.physicalConfigurationIdentity=descriptor_.eventIdentity;
      down.backgroundGeneration=descriptor_.backgroundGeneration;down.geometryGeneration=descriptor_.geometryGeneration;
      down.epochS=epoch;down.positionM=position;down.region=M::Region::Sheath;
      down.side=M::Side::DownstreamLimit;down.interfaceIdentity=result.stableIdentity;
      down.massDensityKgM3=f.downstreamMassDensityKgM3;down.pressurePa=f.downstreamPressurePa;
      down.velocityMPerS=f.downstreamVelocityMPerS;down.magneticFieldT=f.downstreamMagneticFieldT;
      result.downstream=std::move(down);
    }
    return R::Success(std::move(result));
  }
 private:
  const IntervalView view_;
  const M::Descriptor descriptor_;
  const GeometryControls controls_;
};
inline Core::Result<std::shared_ptr<const M::Interval>> Capture(
    const IntervalView& view,std::uint64_t geometry,std::uint64_t snapshot,
    const std::string& nativeIdentity,const GeometryControls& controls) {
  using R=Core::Result<std::shared_ptr<const M::Interval>>;
  // Accuracy/work controls are supplied by the existing geometry caller.
  // The neutral interface neither guesses tolerances nor changes physical dt.
  if(!std::isfinite(controls.positionToleranceM)||controls.positionToleranceM<=0||
      !std::isfinite(controls.maximumTimeSubstepS)||controls.maximumTimeSubstepS<=0||
      !std::isfinite(controls.maximumArcSubstepM)||controls.maximumArcSubstepM<=0||controls.maximumAttempts==0)
    return R::Failure(Core::StatusCode::InvalidConfiguration,"explicit material numerical controls are required");
  const auto descriptor=view.MaterialDescriptor(geometry,snapshot,nativeIdentity);
  if(!descriptor.ok())return R::Failure(descriptor.status.code,descriptor.status.message);
  return M::Interval::Create(std::make_shared<const Owner>(view,descriptor.value,controls));
}
} } }
#endif
