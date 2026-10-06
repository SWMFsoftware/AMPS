#include "reduced_shock_background_adapter.h"

#include "provider.h"

#include <cmath>
#include <filesystem>
#include <fstream>
#include <limits>
#include <sstream>

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

} // namespace

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
  const auto state=gState.provider->QueryAmbient(
      {positionM[0],positionM[1],positionM[2]},gState.metadata.epochS,
      gState.metadata.generation);
  if(!state.ok())return Fail("reduced ambient query failed: "+
      state.status.message,error);
  AmbientSample candidate;
  candidate.magneticFieldT={{state.value.primitive.magneticFieldT.x,
      state.value.primitive.magneticFieldT.y,
      state.value.primitive.magneticFieldT.z}};
  candidate.velocityMPerS={{state.value.primitive.velocityMPerS.x,
      state.value.primitive.velocityMPerS.y,
      state.value.primitive.velocityMPerS.z}};
  candidate.numberDensityM3=
      state.value.primitive.plasma.electronNumberDensityM3;
  candidate.pressurePa=state.value.primitive.plasma.pressurePa;
  candidate.protonTemperatureK=state.value.primitive.protonTemperatureK;
  if(!FiniteSample(candidate))return Fail(
      "reduced ambient query returned nonphysical primitive state",error);
  *sample=candidate;
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
