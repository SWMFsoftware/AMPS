#include "shock_front_background_adapter.h"

#include "../background/bg_provider.h"
#include "provider.h"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <limits>
#include <sstream>

namespace SEP3D { namespace Adapters { namespace {

Core::Status Invalid(const std::string& message) {
  return Core::Status(Core::StatusCode::BackgroundInvalid,message);
}

std::uint64_t Digest(const std::string& text) {
  // DATAFILE samples carry this compact in-memory tag for fast mixed-buffer
  // rejection.  It is derived from, but does not replace, the shared SHA-256
  // event identity recorded in metadata, MPI evidence, and restart receipts.
  std::uint64_t value=1469598103934665603ULL;
  for(unsigned char byte:text){value^=byte;value*=1099511628211ULL;}
  return value;
}

std::string Read(const std::filesystem::path& path) {
  std::ifstream input(path,std::ios::binary);if(!input)return {};
  std::ostringstream bytes;bytes<<input.rdbuf();
  return input.good()||input.eof()?bytes.str():std::string{};
}

Core::Vec3 Convert(SEP::CoronalCME::Vec3 value) {
  return {value.x,value.y,value.z};
}

} // namespace

Core::Status ShockFrontBackgroundAdapter::Create(
    const RuntimeModel::RunConfiguration3D& configuration,
    std::shared_ptr<Background::BackgroundProvider>* output) {
  if(!output)return Invalid("shock-front adapter output pointer is null");
  const auto& options=configuration.options();
  if(options.backgroundModelId!="sep-corona-swcme-shock-front-v1")
    return Invalid("shock-front adapter requires its frozen model id");
  const bool inlineModel=!options.backgroundModelInlineConfiguration.empty();
  const std::filesystem::path eventPath(options.backgroundModelAssetPath);
  const std::filesystem::path assetDirectory=inlineModel?
      std::filesystem::path(options.backgroundModelAssetDirectory):
      eventPath.parent_path();
  const std::string eventBytes=inlineModel?
      options.backgroundModelInlineConfiguration:Read(eventPath);
  if(eventBytes.empty())return Invalid(inlineModel?
      "inline shock-front model configuration is empty":
      "cannot read shock-front event asset '"+eventPath.string()+"'");
  // The application owns filesystem policy only.  Resolve transitive assets
  // beside the event file (or beside the input file that holds the included
  // physics: shock_front block, whose canonical text is the inline
  // configuration), then hand their
  // bytes to the dependency-light shared resolver. Both paths invoke the same
  // strict schema/checksum/physics gate and therefore produce the same event
  // identity for identical physical assignments and magnetic bytes.
  const auto resolved=SEP::CoronaSwcme::ShockFront::ResolveConfiguration(
      eventBytes,[&](const std::string& relative) {
        const auto path=assetDirectory/relative;
        const std::string bytes=Read(path);
        if(bytes.empty())return SEP::Core::Result<std::string>::Failure(
            SEP::Core::StatusCode::OutOfDomain,
            "cannot read reduced asset '"+path.string()+"'");
        return SEP::Core::Result<std::string>::Success(bytes);
      },inlineModel?std::string("srcSEP3D physics: shock_front block"):
          eventPath.string());
  if(!resolved.ok())return Invalid("cannot resolve shock-front event: "+
      resolved.status.message);
  // Both geometries of the shared core are consumed through the
  // provider-neutral ShockSurface (adapters/shock_surface.h).  The mover-side
  // shock record (ExpandingShock) carries Sphere or FiniteSSE accordingly.
  // Any other core geometry is rejected here rather than reaching consumers
  // that know only these two.
  if(resolved.value->geometryModel!=
         SEP::CoronaSwcme::ShockFront::GeometryModel::FiniteSseFixedWidth&&
     resolved.value->geometryModel!=
         SEP::CoronaSwcme::ShockFront::GeometryModel::HeliocentricSphere)
    return Invalid("srcSEP3D reduced shock-front consumers support "
                   "geometry.model=finite-sse-fixed-width or heliocentric-sphere");
  const Core::ComputationFrame frame=options.ComputationFrameValue();
  if(!frame.IsInertial()) {
    // The corotating computation frame is the frame in which the reduced
    // ambient (a rigidly rotating PFSS/Parker pattern) is static.  Any other
    // rate would make the traced tube and the published field time dependent.
    // The comparison is exact (bitwise equality of the two doubles): the
    // deck must state the event's rate literally (2.86533e-6 rad/s in the
    // 2013-04-11 test).
    const double eventRate=resolved.value->ambient.ambient.rotationRateRadPerS;
    if(frame.RotationRateRadPerS()!=eventRate)
      return Invalid("mesh.domain.frame_rotation_rate_rad_per_s must equal the "
                     "event's ambient.rotation_rate_rad_per_s ("+
                     std::to_string(eventRate)+") so the ambient is static in "
                     "the corotating frame");
  }
  const auto shared=SEP::CoronaSwcme::ShockFront::Provider::Create(resolved.value);
  if(!shared.ok())return Invalid("cannot construct shock-front provider: "+
      shared.status.message);
  std::shared_ptr<ShockFrontBackgroundAdapter> candidate(
      new ShockFrontBackgroundAdapter);
  candidate->provider_=shared.value;candidate->originM_=options.coordinateOriginM;
  // Keep the frame for every conversion below and the run's frame label for
  // the published snapshot metadata.
  candidate->frame_=frame;candidate->runFrameName_=options.coordinateFrame;
  candidate->configurationDigest_=Digest(resolved.value->physicsFingerprint);
  const auto valid=candidate->Validate();if(!valid.ok())return valid;
  *output=std::move(candidate);return Core::Status::OK();
}

Core::Status ShockFrontBackgroundAdapter::Validate() const {
  if(!provider_||provider_->Event().particleMode!="disabled"||
      provider_->Event().coordinateFrame!="HCI")return Invalid(
          "shock-front adapter requires a resolved HCI zero-particle event");
  if(!std::isfinite(originM_.x)||!std::isfinite(originM_.y)||
      !std::isfinite(originM_.z))return Invalid("shock-front origin is nonfinite");
  return Core::Status::OK();
}

Core::Status ShockFrontBackgroundAdapter::Prepare(double timeS) {
  const auto valid=Validate();if(!valid.ok())return valid;
  const auto& event=provider_->Event();
  if(!std::isfinite(timeS)||timeS<event.ambient.support.startS||
      timeS>event.ambient.support.endS)return Invalid(
          "shock-front adapter epoch is outside event coverage");
  // The cadence-derived value is a deterministic lower bound, not the whole
  // generation rule.  A particle CFL step can be much shorter than the
  // provider's nominal background/output cadence, yet every distinct front
  // geometry sampled by the live source must still have a distinct immutable
  // generation.  Therefore a forward request advances by at least one from
  // the committed generation.  Sparse 60/600-s background campaigns retain
  // their historical absolute generations, while a 2-s particle campaign no
  // longer attempts to publish many different epochs under generation one.
  // This is an identity correction only: it neither interpolates nor changes
  // the analytical trajectory, ambient state, or shock admission criteria.
  // Long-double arithmetic keeps the absolute lower bound reproducible at
  // heliospheric times and after restart.
  const long double raw=(static_cast<long double>(timeS)-
      event.ambient.support.startS)/event.backgroundDtS;
  if(raw<0||raw>static_cast<long double>(UINT64_MAX-1))return Invalid(
      "shock-front generation would overflow");
  std::uint64_t generation=1+
      static_cast<std::uint64_t>(std::floor(raw+0.5L));
  if(prepared_) {
    if(timeS<metadata_.epochS)
      return Invalid("shock-front epochs must advance monotonically");
    if(timeS>metadata_.epochS) {
      if(metadata_.generation==UINT64_MAX)
        return Invalid("shock-front generation would overflow");
      generation=std::max(generation,metadata_.generation+1);
    } else {
      // Exact-time preparation is idempotent and retains its identity.  This
      // is useful for read-only rebuilding; it must never consume a generation.
      generation=metadata_.generation;
    }
  }
  // Provider::Prepare is transactional.  Publish the application metadata
  // only after the immutable front/ambient candidate commits, so owner arrays,
  // ghost exchange, and the diagnostic surface cannot advertise different
  // clocks or generations after a rejected candidate.
  const auto candidate=provider_->Prepare(timeS,generation);
  if(!candidate.ok())return Invalid("shock-front candidate rejected: "+
      candidate.status.message);
  Background::SnapshotMetadata next;
  next.provider=Background::ProviderKind::RuntimeModel;
  next.ownership=Background::StorageOwnership::ModelOwned;
  next.epochS=next.validFromS=timeS;
  // In the corotating frame of the event's rotation the published ambient
  // is static (a rigidly rotating pattern): it stays valid to the end of the
  // event support, so no republication is needed.  In HCI the pattern turns,
  // so the snapshot keeps the event's background cadence backgroundDtS.
  next.validUntilS=frame_.IsInertial()?
      std::min(event.ambient.support.endS,timeS+event.backgroundDtS):
      event.ambient.support.endS;
  next.generation=generation;
  // The published cells are in the run's computation frame.
  next.coordinateFrame=runFrameName_.empty()?event.coordinateFrame:runFrameName_;
  next.providerIdentity=CanonicalName();
  next.configurationFingerprint=event.physicsFingerprint;
  metadata_=std::move(next);prepared_=true;return Core::Status::OK();
}

const Background::SnapshotMetadata*
ShockFrontBackgroundAdapter::PreparedMetadata() const {
  return prepared_?&metadata_:nullptr;
}

Core::Status ShockFrontBackgroundAdapter::RestorePreparedGeneration(
    std::uint64_t generation) {
  if(!prepared_||generation==0)return Invalid(
      "shock-front restart generation requires a prepared nonzero epoch");
  metadata_.generation=generation;return Core::Status::OK();
}

Background::BackgroundSample ShockFrontBackgroundAdapter::Evaluate(
    const Core::Vec3& global) const {
  Background::BackgroundSample out;
  if(!prepared_) {out.status=Invalid("shock-front adapter is not prepared");return out;}
  // Runtime coordinates may have an application origin offset, while all
  // shared physics is heliocentric HCI.  Evaluate the canonical undisturbed
  // ambient at that HCI position.  Deliberately do not inspect whether the
  // point lies behind the front: the reduced provider has no physical sheath
  // or ejecta volume, and painting a local RH state into cells would fabricate
  // one.  Immediate downstream exists only in FrontEpoch() records.
  // In a corotating run the cell position is a frame point: rotate it into
  // HCI for the query and every returned vector/tensor back into the frame.
  // B rotates as a true vector; U also loses the frame velocity Omega x x at
  // the HCI point (U' = U - Omega x x, field-aligned for the corotating
  // ambient).  Scalars (density, temperature, pressure, Alfven speed) are
  // frame independent.
  const Core::Vec3 local=global-originM_;
  const double epoch=metadata_.epochS;
  const Core::Vec3 hci=frame_.ToInertialPosition(local,epoch);
  const auto state=provider_->QueryAmbient({hci.x,hci.y,hci.z},
      epoch,metadata_.generation);
  if(!state.ok()) {out.status=Invalid(state.status.message);return out;}
  const auto& value=state.value;
  out.B=frame_.FromInertialVector(Convert(value.primitive.magneticFieldT),epoch);
  out.U=frame_.FromInertialVelocity(Convert(value.primitive.velocityMPerS),
                                    hci,epoch);
  out.absB=out.B.Norm();out.bHat=out.absB>0?out.B/out.absB:Core::Vec3{};
  out.numberDensityM3=value.primitive.plasma.electronNumberDensityM3;
  out.temperatureK=value.primitive.protonTemperatureK;
  out.pressurePa=value.primitive.plasma.pressurePa;
  out.alfvenSpeedMpS=value.primitive.plasma.alfvenSpeedMPerS;
  // The core's gradients are row-major d(component)/d(HCI coordinate).  In
  // the frame they become R G R^T with R = Rz(-theta); grad U also loses the
  // antisymmetric grad(Omega x x), which leaves div U and bb:grad U
  // unchanged (core/computation_frame.h).
  Core::ComputationFrame::Tensor gradB{},gradU{};
  for(int i=0;i<9;++i){gradB[i]=value.gradientB[i];gradU[i]=value.gradientU[i];}
  gradB=frame_.FromInertialTensor(gradB,epoch);
  gradU=frame_.FromInertialVelocityGradient(gradU,epoch);
  for(int component=0;component<3;++component)for(int coordinate=0;
      coordinate<3;++coordinate) {
    out.gradB(component,coordinate)=gradB[3*component+coordinate];
    out.gradU(component,coordinate)=gradU[3*component+coordinate];
  }
  Background::CompleteVectorDerivatives(&out);
  out.generation=metadata_.generation;out.configurationDigest=configurationDigest_;
  out.valid=true;out.status=Core::Status::OK();return out;
}

std::string ShockFrontBackgroundAdapter::ResolvedManifest() const {
  return provider_?provider_->Event().normalizedManifest:std::string{};
}

Background::ProviderCapabilities ShockFrontBackgroundAdapter::Capabilities() const {
  Background::ProviderCapabilities out;out.hasGradB=out.hasDivBhat=
      out.hasCurvature=out.hasGradU=out.hasFieldAlignedStrain=
      out.hasPlasmaState=out.supportsBatchEval=true;return out;
}

// Ambient B, U' and sector at a frame point (see the header).  Unlike
// Evaluate() it needs no prepared snapshot and takes its own epoch, so the
// observer-tube tracer can call it while the mesh is being built.  Errors:
// no provider, or the core's ambient evaluation fails (e.g. a point where
// the open-tube wind has no solution).
Core::Status ShockFrontBackgroundAdapter::AmbientAt(
    const Core::Vec3& framePositionM,double epochS,Core::Vec3* magneticFieldT,
    Core::Vec3* velocityMPerS,int* sector) const {
  if(!provider_)return Invalid("shock-front adapter has no provider");
  const Core::Vec3 hci=frame_.ToInertialPosition(framePositionM-originM_,epochS);
  const auto state=provider_->Ambient().Evaluate({hci.x,hci.y,hci.z},epochS);
  if(!state.ok())return Invalid(state.status.message);
  if(magneticFieldT)*magneticFieldT=frame_.FromInertialVector(
      Convert(state.value.magneticFieldT),epochS);
  if(velocityMPerS)*velocityMPerS=frame_.FromInertialVelocity(
      Convert(state.value.velocityMPerS),hci,epochS);
  if(sector)*sector=state.value.magneticSector;
  return Core::Status::OK();
}

std::shared_ptr<const SEP::CoronaSwcme::ShockFront::Epoch>
ShockFrontBackgroundAdapter::FrontEpoch() const {
  return provider_?provider_->Current():nullptr;
}

} } // namespace SEP3D::Adapters
