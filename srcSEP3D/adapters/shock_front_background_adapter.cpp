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
  if(options.backgroundModelId!="sep-corona-swcme-shock-front-v1"||
      options.backgroundModelAssetPath.empty())return Invalid(
          "shock-front adapter requires its frozen model id and event asset");
  const std::filesystem::path eventPath(options.backgroundModelAssetPath);
  const std::string eventBytes=Read(eventPath);
  if(eventBytes.empty())return Invalid("cannot read shock-front event asset '"+
      eventPath.string()+"'");
  // The application owns filesystem policy only.  Resolve transitive assets
  // beside the event deck, then hand their bytes to the dependency-light
  // shared resolver, which verifies checksums/physics and forms one identity.
  const auto resolved=SEP::CoronaSwcme::ShockFront::ResolveConfiguration(
      eventBytes,[&](const std::string& relative) {
        const auto path=eventPath.parent_path()/relative;
        const std::string bytes=Read(path);
        if(bytes.empty())return SEP::Core::Result<std::string>::Failure(
            SEP::Core::StatusCode::OutOfDomain,
            "cannot read reduced asset '"+path.string()+"'");
        return SEP::Core::Result<std::string>::Success(bytes);
      });
  if(!resolved.ok())return Invalid("cannot resolve shock-front event: "+
      resolved.status.message);
  const auto shared=SEP::CoronaSwcme::ShockFront::Provider::Create(resolved.value);
  if(!shared.ok())return Invalid("cannot construct shock-front provider: "+
      shared.status.message);
  std::shared_ptr<ShockFrontBackgroundAdapter> candidate(
      new ShockFrontBackgroundAdapter);
  candidate->provider_=shared.value;candidate->originM_=options.coordinateOriginM;
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
  // Generation is tied to the declared background cadence rather than an
  // application call count.  Long-double arithmetic avoids losing a cadence
  // at heliospheric absolute times; the half-step rounding admits only times
  // that the host parser has already constrained to its cadence.
  const long double raw=(static_cast<long double>(timeS)-
      event.ambient.support.startS)/event.backgroundDtS;
  if(raw<0||raw>static_cast<long double>(UINT64_MAX-1))return Invalid(
      "shock-front generation would overflow");
  const std::uint64_t generation=1+static_cast<std::uint64_t>(std::floor(raw+0.5L));
  if(prepared_&&generation<=metadata_.generation&&timeS!=metadata_.epochS)
    return Invalid("shock-front epochs must advance monotonically");
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
  next.validUntilS=std::min(event.ambient.support.endS,timeS+event.backgroundDtS);
  next.generation=generation;next.coordinateFrame=event.coordinateFrame;
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
  const Core::Vec3 local=global-originM_;
  const auto state=provider_->QueryAmbient({local.x,local.y,local.z},
      metadata_.epochS,metadata_.generation);
  if(!state.ok()) {out.status=Invalid(state.status.message);return out;}
  const auto& value=state.value;
  out.B=Convert(value.primitive.magneticFieldT);
  out.U=Convert(value.primitive.velocityMPerS);
  out.absB=out.B.Norm();out.bHat=out.absB>0?out.B/out.absB:Core::Vec3{};
  out.numberDensityM3=value.primitive.plasma.electronNumberDensityM3;
  out.temperatureK=value.primitive.protonTemperatureK;
  out.pressurePa=value.primitive.plasma.pressurePa;
  out.alfvenSpeedMpS=value.primitive.plasma.alfvenSpeedMPerS;
  for(int component=0;component<3;++component)for(int coordinate=0;
      coordinate<3;++coordinate) {
    out.gradB(component,coordinate)=value.gradientB[3*component+coordinate];
    out.gradU(component,coordinate)=value.gradientU[3*component+coordinate];
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

std::shared_ptr<const SEP::CoronaSwcme::ShockFront::Epoch>
ShockFrontBackgroundAdapter::FrontEpoch() const {
  return provider_?provider_->Current():nullptr;
}

} } // namespace SEP3D::Adapters
