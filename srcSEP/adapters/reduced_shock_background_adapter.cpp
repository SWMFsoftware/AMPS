#include "reduced_shock_background_adapter.h"

#include "provider.h"

#include <cmath>
#include <filesystem>
#include <fstream>
#include <limits>
#include <sstream>
#include <iomanip>

namespace SEP { namespace ReducedShock { namespace {

struct AdapterState {
  std::shared_ptr<CoronaSwcme::ShockFront::Provider> provider;
  EpochMetadata metadata;
  std::string eventPath;
  bool prepared = false;
};

AdapterState gState;

std::string ReadFile(const std::filesystem::path& path) {
  std::ifstream input(path,std::ios::binary);
  if(!input)return {};
  std::ostringstream bytes;
  bytes<<input.rdbuf();
  return input.good()||input.eof()?bytes.str():std::string{};
}

bool Fail(const std::string& message,std::string* error) {
  if(error)*error=message;
  return false;
}

bool FiniteSample(const AmbientSample& sample) {
  for(double value:sample.magneticFieldT)if(!std::isfinite(value))return false;
  for(double value:sample.velocityMPerS)if(!std::isfinite(value))return false;
  return std::isfinite(sample.numberDensityM3)&&sample.numberDensityM3>0.0&&
      std::isfinite(sample.pressurePa)&&sample.pressurePa>0.0&&
      std::isfinite(sample.protonTemperatureK)&&sample.protonTemperatureK>0.0;
}

// Conversion is a value adapter only: EOS, composition, thermal partition and
// vector basis are supplied by the SAME core authority that solves the jump.
// Returned temporaries give failure atomicity to all public query methods.
Core::Result<AmbientSample> ConvertAmbient(
    const CoronaSwcme::AmbientPrimitive& primitive,
    const CoronaSwcme::ShockFront::Configuration& event,
    double epoch,std::uint64_t generation) {
  AmbientSample out;
  out.magneticFieldT={{primitive.magneticFieldT.x,primitive.magneticFieldT.y,primitive.magneticFieldT.z}};
  out.velocityMPerS={{primitive.velocityMPerS.x,primitive.velocityMPerS.y,primitive.velocityMPerS.z}};
  out.numberDensityM3=out.electronNumberDensityM3=primitive.plasma.electronNumberDensityM3;
  out.protonNumberDensityM3=primitive.plasma.protonNumberDensityM3;
  out.alphaNumberDensityM3=event.ambient.composition.alphaToProtonNumberRatio*out.protonNumberDensityM3;
  out.massDensityKgM3=primitive.plasma.massDensityKgM3;
  out.pressurePa=primitive.plasma.pressurePa;out.protonTemperatureK=primitive.protonTemperatureK;
  out.parkerExterior=primitive.region==SEP::CoronaSwcme::AmbientRegion::ParkerExterior;
  out.ambientRegion=static_cast<int>(primitive.region);out.magneticSector=primitive.magneticSector;
  out.epochS=epoch;out.generation=generation;out.eventIdentity=event.physicsFingerprint;
  if(!FiniteSample(out)||!std::isfinite(out.protonNumberDensityM3)||out.protonNumberDensityM3<=0||
      !std::isfinite(out.alphaNumberDensityM3)||out.alphaNumberDensityM3<0||
      !std::isfinite(out.massDensityKgM3)||out.massDensityKgM3<=0)
    return Core::Result<AmbientSample>::Failure(Core::StatusCode::InvalidState,
        "core reduced ambient conversion found invalid species/thermal state");
  return Core::Result<AmbientSample>::Success(std::move(out));
}

} // namespace

Core::Result<IntervalView> CaptureInterval(double departure,double arrival,
    const std::string& expectedIdentity) {
  using Return=Core::Result<IntervalView>;
  if(!gState.prepared||!gState.provider||expectedIdentity.empty()||
      expectedIdentity!=gState.metadata.eventIdentity)
    return Return::Failure(Core::StatusCode::InvalidState,"reduced interval requires the prepared expected event identity");
  if(!std::isfinite(departure)||!std::isfinite(arrival)||arrival<departure||
      departure<gState.metadata.epochS||arrival>gState.metadata.validUntilS)
    return Return::Failure(Core::StatusCode::OutOfDomain,"reduced interval exceeds its committed physical validity; AMPS dt must not be shortened");
  IntervalView view;view.provider_=gState.provider;view.departureS_=departure;view.arrivalS_=arrival;
  view.generation_=gState.metadata.generation;view.identity_=expectedIdentity;
  return Return::Success(std::move(view));
}
Core::Status IntervalView::ValidateEpoch(double epoch) const {
  if(!provider_||generation_==0||identity_!=provider_->Event().physicsFingerprint)
    return Core::Status::Failure(Core::StatusCode::InvalidState,"unbound or foreign reduced interval");
  if(!std::isfinite(epoch)||epoch<departureS_||epoch>arrivalS_) {
    std::ostringstream why;why<<std::setprecision(17)<<"query epoch "<<epoch
        <<" s is outside its immutable reduced interval ["<<departureS_<<", "<<arrivalS_<<"] s";
    return Core::Status::Failure(Core::StatusCode::OutOfDomain,why.str());
  }
  return Core::Status::Success();
}
Core::Result<IntervalView> CapturePhysicalStep(double departure,double arrival,const std::string& identity) {
  using Return=Core::Result<IntervalView>;
  if(!gState.prepared||!gState.provider||identity!=gState.metadata.eventIdentity)
    return Return::Failure(Core::StatusCode::InvalidState,"physical interval requires its prepared core authority");
  const auto& support=gState.provider->Event().ambient.support;
  if(!std::isfinite(departure)||!std::isfinite(arrival)||arrival<departure||
      departure<support.startS||arrival>support.endS)
    return Return::Failure(Core::StatusCode::OutOfDomain,"full physical step exceeds prescribed core event support");
  IntervalView view;view.provider_=gState.provider;view.departureS_=departure;view.arrivalS_=arrival;
  view.generation_=gState.metadata.generation;view.identity_=identity;
  return Return::Success(std::move(view));
}
bool EnsurePreparedForEpoch(double epoch,std::string* error) {
  if(!gState.provider)return Fail("reduced provider is not configured",error);
  const auto& event=gState.provider->Event();
  if(!std::isfinite(epoch)||epoch<event.ambient.support.startS||epoch>event.ambient.support.endS)
    return Fail("application epoch is outside reduced event support",error);
  // This selects the last cache cadence, not a replacement application time.
  // Ambient/front queries still evaluate their actual requested epoch. A full
  // AMPS step can cross this cache cadence because the analytical functions,
  // unlike imported MHD arrays, are defined continuously over event support.
  const double tick=event.ambient.support.startS+event.backgroundDtS*
      std::floor((epoch-event.ambient.support.startS)/event.backgroundDtS);
  return Prepare(tick,error);
}
double EventEndS() {return gState.provider?gState.provider->Event().ambient.support.endS:0;}
Core::Result<AmbientSample> IntervalView::PrimitiveAt(const Position& x,double epoch) const {
  using Return=Core::Result<AmbientSample>;
  const auto valid=ValidateEpoch(epoch);if(!valid.ok())return Return::Failure(valid.code,valid.message);
  const auto state=provider_->Ambient().Evaluate({x[0],x[1],x[2]},epoch);
  if(!state.ok())return Return::Failure(state.status.code,state.status.message);
  return ConvertAmbient(state.value,provider_->Event(),epoch,generation_);
}
Core::Result<AmbientSample> IntervalView::AmbientAt(const Position& x,double epoch) const {
  using Return=Core::Result<AmbientSample>;
  const auto valid=ValidateEpoch(epoch);if(!valid.ok())return Return::Failure(valid.code,valid.message);
  const auto state=provider_->QueryAmbient({x[0],x[1],x[2]},epoch,generation_);
  if(!state.ok())return Return::Failure(state.status.code,state.status.message);
  if(!state.value.derivativesValid||state.value.epochS!=epoch||
      state.value.generation!=generation_||state.value.eventIdentity!=identity_)
    return Return::Failure(Core::StatusCode::DataIntegrityFailure,"core reduced derivative provenance differs from its captured interval");
  auto out=ConvertAmbient(state.value.primitive,provider_->Event(),epoch,generation_);
  if(!out.ok())return out;
  out.value.velocityGradientPerS=state.value.gradientU;
  out.value.magneticGradientTPerM=state.value.gradientB;
  out.value.logMagneticGradientPerM={{state.value.gradientLogBPerM.x,state.value.gradientLogBPerM.y,state.value.gradientLogBPerM.z}};
  for(int i=0;i<9;++i)if(!std::isfinite(out.value.velocityGradientPerS[i])||!std::isfinite(out.value.magneticGradientTPerM[i]))
    return Return::Failure(Core::StatusCode::NumericalFailure,"core reduced derivative tensor is nonfinite");
  for(int i=0;i<3;++i) {
    if(!std::isfinite(out.value.logMagneticGradientPerM[i]))return Return::Failure(Core::StatusCode::NumericalFailure,"core reduced magnetic logarithmic derivative is nonfinite");
    out.value.derivativeStencils[i]=static_cast<int>(state.value.stencils[i]);
  }
  out.value.derivativesValid=true;return out;
}
Core::Status IntervalView::RequirePhysicalDownstream() const {
  const auto valid=ValidateEpoch(departureS_);if(!valid.ok())return valid;
  return provider_->RequireCapability(CoronaSwcme::ShockFront::VolumeCapability::PhysicalDownstreamVolume);
}
Core::Result<SurfaceSphere> IntervalView::SphereAt(double epoch) const {
  using Return=Core::Result<SurfaceSphere>;
  const auto valid=ValidateEpoch(epoch);if(!valid.ok())return Return::Failure(valid.code,valid.message);
  const auto sphere=provider_->Sphere(epoch);
  if(!sphere.ok())return Return::Failure(sphere.status.code,sphere.status.message);
  return Return::Success({{{sphere.value.centerM.x,sphere.value.centerM.y,sphere.value.centerM.z}},
      {{sphere.value.direction.x,sphere.value.direction.y,sphere.value.direction.z}},sphere.value.radiusM,sphere.value.minimumRayCosine});
}
Core::Result<Position> IntervalView::LeadingPointOnRay(const Position& ray,double epoch) const {
  using Return=Core::Result<Position>;
  const auto valid=ValidateEpoch(epoch);if(!valid.ok())return Return::Failure(valid.code,valid.message);
  const auto geometry=provider_->EvaluateRay({ray[0],ray[1],ray[2]},epoch);
  if(!geometry.ok())return Return::Failure(geometry.status.code,geometry.status.message);
  return Return::Success({{geometry.value.positionM.x,geometry.value.positionM.y,geometry.value.positionM.z}});
}
Core::Result<FrontSample> IntervalView::FrontAt(const Position& x,double epoch) const {
  using Return=Core::Result<FrontSample>;
  const auto valid=ValidateEpoch(epoch);if(!valid.ok())return Return::Failure(valid.code,valid.message);
  // Separate unsupported geometry from unavailable physics on a SUPPORTED
  // point. A failed ambient/jump query must never disappear as "no connection".
  const auto point=LeadingPointOnRay(x,epoch);if(!point.ok())return Return::Failure(point.status.code,point.status.message);
  double mismatch=0,radius=0;for(int d=0;d<3;++d){mismatch=std::hypot(mismatch,x[d]-point.value[d]);radius=std::hypot(radius,x[d]);}
  if(mismatch>256*std::numeric_limits<double>::epsilon()*std::max(1.,radius))
    return Return::Failure(Core::StatusCode::OutOfDomain,"point is not on the core leading surface");
  const auto sample=provider_->EvaluateFrontPoint({x[0],x[1],x[2]},epoch);
  if(!sample.ok())return Return::Failure(Core::StatusCode::NumericalFailure,
      "supported core front state unavailable: "+sample.status.message);
  const auto& s=sample.value;FrontSample out;
  if(s.status==CoronaSwcme::ShockFront::FrontStatus::NumericallyUnresolvedWeakShock||
      s.status==CoronaSwcme::ShockFront::FrontStatus::WrongBranch||
      s.status==CoronaSwcme::ShockFront::FrontStatus::InvalidJump)
    return Return::Failure(Core::StatusCode::NumericalFailure,"supported reduced jump is numerically unresolved: "+s.reason);
  out.positionM=x;out.outwardNormal={{s.geometry.outwardNormal.x,s.geometry.outwardNormal.y,s.geometry.outwardNormal.z}};
  out.normalSpeedMPerS=s.geometry.normalSpeedMPerS;out.supportEdge=s.geometry.supportEdge;
  out.inflowMPerS=s.inflowMPerS;out.fastMach=s.fastMach;
  out.classification=CoronaSwcme::ShockFront::Name(s.status);
  out.acceptedShock=s.status==CoronaSwcme::ShockFront::FrontStatus::SolvedFastShock;
  out.downstreamValid=s.downstreamValid;out.epochS=epoch;out.generation=generation_;out.eventIdentity=identity_;
  out.upstreamValid=s.status!=CoronaSwcme::ShockFront::FrontStatus::BelowPhysicalInnerBoundary;
  if(out.upstreamValid) {
    const auto upstream=ConvertAmbient(s.upstream,provider_->Event(),epoch,generation_);
    if(!upstream.ok())return Return::Failure(upstream.status.code,upstream.status.message);
    out.upstream=upstream.value;
  }
  if(out.downstreamValid) {
    out.compressionRatio=s.jump.compressionRatio;
    out.downstreamMassDensityKgM3=s.jump.downstream.massDensityKgM3;
    out.downstreamPressurePa=s.jump.downstream.pressurePa;
    const auto& b=s.jump.downstream.magneticFieldT;out.downstreamMagneticFieldT={{b.x,b.y,b.z}};
    const auto& u=s.jump.downstream.velocityMPerS;out.downstreamVelocityMPerS={{u.x,u.y,u.z}};
  }
  return Return::Success(std::move(out));
}

bool Configure(const std::string& eventPath,std::string* error) {
  if(eventPath.empty())return Fail("reduced shock event path is empty",error);
  const std::filesystem::path path(eventPath);
  const std::string eventBytes=ReadFile(path);
  if(eventBytes.empty())return Fail("cannot read reduced shock event '"+
      path.string()+"'",error);

  // The application resolves files, but the shared resolver owns schema,
  // checksums, units, physics validation, and the normalized event identity.
  // No application default is allowed to fill an omitted physical parameter.
  const auto resolved=CoronaSwcme::ShockFront::ResolveConfiguration(
      eventBytes,[&](const std::string& relative) {
        const std::filesystem::path asset=path.parent_path()/relative;
        const std::string bytes=ReadFile(asset);
        if(bytes.empty())return Core::Result<std::string>::Failure(
            Core::StatusCode::OutOfDomain,
            "cannot read reduced shock asset '"+asset.string()+"'");
        return Core::Result<std::string>::Success(bytes);
      });
  if(!resolved.ok())return Fail("cannot resolve reduced shock event: "+
      resolved.status.message,error);
  if(resolved.value->coordinateFrame!="HCI"||
      resolved.value->particleMode!="disabled")return Fail(
          "srcSEP reduced coupling requires HCI and particle_mode=disabled",error);
  const auto provider=CoronaSwcme::ShockFront::Provider::Create(resolved.value);
  if(!provider.ok())return Fail("cannot create reduced shock provider: "+
      provider.status.message,error);

  AdapterState candidate;
  candidate.provider=provider.value;
  candidate.eventPath=path.string();
  gState=std::move(candidate);
  if(error)error->clear();
  return true;
}

bool Enabled(){return static_cast<bool>(gState.provider);}

Core::Result<CadencePublication> PrepareCadencePublication(double epoch,bool checkpointRestore) {
  using Return=Core::Result<CadencePublication>;
  if(!gState.provider||!gState.prepared)return Return::Failure(Core::StatusCode::InvalidState,"cadence publication needs a configured initial core epoch");
  const auto& event=gState.provider->Event();
  if(!std::isfinite(epoch)||epoch<event.ambient.support.startS||epoch>event.ambient.support.endS)
    return Return::Failure(Core::StatusCode::OutOfDomain,"cadence publication is outside physical event support");
  const long double tick=std::floor((static_cast<long double>(epoch)-event.ambient.support.startS)/event.backgroundDtS);
  if(tick<0||tick>static_cast<long double>(UINT64_MAX-1))return Return::Failure(Core::StatusCode::InvalidState,"cadence generation is not representable");
  const auto generation=1+static_cast<std::uint64_t>(tick);
  const double time=event.ambient.support.startS+static_cast<double>(tick)*event.backgroundDtS;
  if(time<gState.metadata.epochS&&!checkpointRestore)return Return::Failure(Core::StatusCode::InvalidState,"cadence publication cannot reverse the committed cache");
  CadencePublication publication;publication.expected_=gState.provider;
  publication.expectedEpoch_=gState.provider->Current();publication.metadata_=gState.metadata;
  publication.candidate_=gState.provider;
  if(time!=gState.metadata.epochS) {
    // Provider's value copy shares the immutable resolved configuration and
    // ambient authority; only its private Current() pointer is independent.
    // All geometry/RH work still uses the maintained production Prepare path.
    // This is cache staging, never an alternative field or shock equation.
    publication.candidate_=std::make_shared<CoronaSwcme::ShockFront::Provider>(*gState.provider);
    const auto prepared=publication.candidate_->Prepare(time,generation);
    if(!prepared.ok())return Return::Failure(prepared.status.code,prepared.status.message);
    publication.metadata_={time,std::min(event.ambient.support.endS,time+event.backgroundDtS),generation,event.physicsFingerprint};
  }
  return Return::Success(std::move(publication));
}
Core::Status CommitCadencePublication(CadencePublication* publication) {
  if(!publication||!publication->expected_||publication->expected_!=gState.provider||
      publication->expectedEpoch_!=gState.provider->Current())
    return Core::Status::Failure(Core::StatusCode::InvalidState,"stale or consumed reduced cadence publication");
  // Allocation and physics validation preceded readiness. The joined native
  // writer is the only publisher; move/swap cannot allocate or rebind an old
  // read-only IntervalView to this new cache owner.
  gState.provider.swap(publication->candidate_);std::swap(gState.metadata,publication->metadata_);
  publication->expected_.reset();publication->expectedEpoch_.reset();return Core::Status::Success();
}
Core::Result<IntervalView> CadencePublication::PhysicalInterval(double departure,double arrival) const {
  using Return=Core::Result<IntervalView>;
  if(!expected_||!candidate_)return Return::Failure(Core::StatusCode::InvalidState,"missing staged cadence owner");
  const auto& support=candidate_->Event().ambient.support;
  if(!std::isfinite(departure)||!std::isfinite(arrival)||arrival<departure||
      departure<support.startS||arrival>support.endS)
    return Return::Failure(Core::StatusCode::OutOfDomain,"staged material interval exceeds the saved event support");
  IntervalView view;view.provider_=candidate_;view.departureS_=departure;view.arrivalS_=arrival;
  view.generation_=metadata_.generation;view.identity_=metadata_.eventIdentity;
  return Return::Success(std::move(view));
}

bool Prepare(double epochS,std::string* error) {
  if(!Enabled())return Fail("reduced shock provider is not configured",error);
  const auto& event=gState.provider->Event();
  if(!std::isfinite(epochS)||epochS<event.ambient.support.startS||
      epochS>event.ambient.support.endS)return Fail(
          "reduced shock epoch is outside declared event coverage",error);
  if(gState.prepared&&epochS==gState.metadata.epochS)return true;

  // Convert the physical clock to the stable event generation.  Rounding to
  // the nearest cadence tick tolerates only floating-point clock noise; the
  // explicit distance-to-tick check rejects arbitrary sub-cadence updates.
  const long double raw=(static_cast<long double>(epochS)-
      event.ambient.support.startS)/event.backgroundDtS;
  const long double nearest=std::round(raw);
  const double reconstructed=event.ambient.support.startS+
      static_cast<double>(nearest)*event.backgroundDtS;
  const double clockTolerance=64.0*std::numeric_limits<double>::epsilon()*
      std::max(1.0,std::fabs(epochS));
  if(nearest<0.0L||nearest>static_cast<long double>(UINT64_MAX-1)||
      std::fabs(reconstructed-epochS)>clockTolerance)return Fail(
          "srcSEP time step is not aligned with event background_dt_s",error);
  const std::uint64_t generation=1+static_cast<std::uint64_t>(nearest);
  if(gState.prepared&&(epochS<gState.metadata.epochS||
      generation<=gState.metadata.generation))return Fail(
          "reduced shock epochs/generations must advance monotonically",error);

  // Provider::Prepare builds an immutable candidate first.  The adapter does
  // not publish metadata until that candidate has passed shared front,
  // ambient, and jump validation, preserving the old epoch after any failure.
  const auto epoch=gState.provider->Prepare(epochS,generation);
  if(!epoch.ok())return Fail("reduced shock candidate rejected: "+
      epoch.status.message,error);
  EpochMetadata metadata;
  metadata.epochS=epochS;
  metadata.validUntilS=std::min(event.ambient.support.endS,
      epochS+event.backgroundDtS);
  metadata.generation=generation;
  metadata.eventIdentity=event.physicsFingerprint;
  gState.metadata=std::move(metadata);
  gState.prepared=true;
  if(error)error->clear();
  return true;
}

bool EvaluateAmbient(const std::array<double,3>& positionM,
                     AmbientSample* sample,std::string* error) {
  if(!sample)return Fail("ambient sample output pointer is null",error);
  if(!gState.prepared)return Fail("reduced shock epoch is not prepared",error);
  const auto view=CaptureInterval(gState.metadata.epochS,gState.metadata.validUntilS,gState.metadata.eventIdentity);
  if(!view.ok())return Fail(view.status.message,error);
  const auto state=view.value.AmbientAt(positionM,gState.metadata.epochS);
  if(!state.ok())return Fail("reduced ambient query failed: "+state.status.message,error);
  *sample=state.value;
  if(error)error->clear();
  return true;
}

const EpochMetadata& Metadata(){return gState.metadata;}
double BackgroundCadenceS(){return gState.provider?
    gState.provider->Event().backgroundDtS:0.0;}
const std::string& EventPath(){return gState.eventPath;}
const std::string& EventIdentity(){
  static const std::string empty;
  return gState.provider?gState.provider->Event().physicsFingerprint:empty;
}
std::shared_ptr<const CoronaSwcme::ShockFront::Epoch> FrontEpoch(){
  return gState.provider?gState.provider->Current():nullptr;
}
FrontSummary CurrentFrontSummary(){
  FrontSummary summary;
  const auto epoch=FrontEpoch();
  if(epoch) {
    summary.apexRadiusM=epoch->trajectory.apexRadiusM;
    summary.acceptedShockAreaM2=epoch->area.acceptedShockM2;
    summary.apexShockAccepted=epoch->apexShockAccepted;
  }
  return summary;
}

} } // namespace SEP::ReducedShock
