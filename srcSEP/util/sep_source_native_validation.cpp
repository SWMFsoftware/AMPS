#include "sep_source_native_validation.h"
#include "../source_capabilities.h"
#include "../sep.h"
#include "../transport_common.h"
#include "sep_source_lifecycle.h"

#include <cstring>
#include <limits>
#include <sstream>

namespace SEP { namespace Testing { namespace {
Result RejectionIsObservational() {
  struct ModeRestore {
    bool value;
    int momentum;
    ~ModeRestore(){SEP::BoundingBoxInjection::BoundaryInjectionMode=value;
      SEP::FieldLine::InjectionParameters::InjectionMomentumModel=momentum;}
  } modeRestore{SEP::BoundingBoxInjection::BoundaryInjectionMode,
      SEP::FieldLine::InjectionParameters::InjectionMomentumModel};
  const auto sourceBefore=SourceLifecycle::CaptureNativeState();
  const long int particles=PIC::ParticleBuffer::GetAllPartNum();
  const auto seed=RandomNumberGenerator::rndLastSeed;
  const auto generator=RandomNumberGenerator::gen;
  const auto distribution=RandomNumberGenerator::dist;
  const bool hasNormal=RandomNumberGenerator::NormalSpareAvailable;
  const double normal=RandomNumberGenerator::NormalSpare;
  unsigned failures=0;
  const bool expectedClock=_PIC_COUPLER_MODE_!=_PIC_COUPLER_MODE__SWMF_&&
      (_SIMULATION_TIME_STEP_MODE_==_SINGLE_GLOBAL_TIME_STEP_||
       (_SIMULATION_TIME_STEP_MODE_==_SPECIES_DEPENDENT_GLOBAL_TIME_STEP_&&
        PIC::nTotalSpecies==1&&_PIC_GLOBAL_TIME_COUNTER_MODE_==_PIC_MODE_ON_));
  if(Transport::PICAdapter::NativeDriverClockCapability().ok()!=expectedClock)++failures;
  for(auto source:{Sources::CartesianSource::GalacticCosmicRayBox,
      Sources::CartesianSource::ThermalOuterBox,Sources::CartesianSource::SphericalInnerBoundary,
      Sources::CartesianSource::SphericalShock}) {
    const auto status=Sources::CartesianBirthPreflight(source);
    if(status.code!=Transport::StatusCode::UnsupportedConfiguration||status.message.empty())++failures;
  }
  const auto invalid=Sources::CartesianBirthPreflight(static_cast<Sources::CartesianSource>(999));
  if(invalid.code!=Transport::StatusCode::InvalidArgument)++failures;
  if(Sources::MaxwellianFieldLinePreflight().code!=Transport::StatusCode::UnsupportedConfiguration)++failures;
  if(Sources::ThermalPopulationPreflight().code!=Transport::StatusCode::UnsupportedConfiguration)++failures;
  if(Sources::ThermalBeginningPreflight().code!=Transport::StatusCode::UnsupportedConfiguration)++failures;
  if(!Sources::MovingShockDeparturePreflight(1.0,1.0).ok()||
      Sources::MovingShockDeparturePreflight(1.0,1.01).code!=Transport::StatusCode::UnsupportedConfiguration||
      Sources::MovingShockDeparturePreflight(std::numeric_limits<double>::quiet_NaN(),1.0).code!=Transport::StatusCode::InvalidArgument)
    ++failures;
  SEP::FieldLine::InjectionParameters::InjectionMomentumModel=
      SEP::FieldLine::InjectionParameters::_tenishev2005aiaa;
  SEP::BoundingBoxInjection::BoundaryInjectionMode=false;
  if(!Sources::ValidateConfiguredSources().ok())++failures;
  SEP::BoundingBoxInjection::BoundaryInjectionMode=true;
  const auto enabled=Sources::ValidateConfiguredSources();
#if _DOMAIN_GEOMETRY_ == _DOMAIN_GEOMETRY_BOX_ && _PIC_EARTH_GCR__MODE_ == _PIC_MODE_ON_
  if(enabled.code!=Transport::StatusCode::UnsupportedConfiguration)++failures;
#else
  if(!enabled.ok())++failures;
#endif
  SEP::BoundingBoxInjection::BoundaryInjectionMode=false;
  SEP::FieldLine::InjectionParameters::InjectionMomentumModel=
      SEP::FieldLine::InjectionParameters::_background_sw_temperature;
  if(Sources::ValidateConfiguredSources().code!=Transport::StatusCode::UnsupportedConfiguration)++failures;
  const auto sourceAfter=SourceLifecycle::CaptureNativeState();
  bool sameSource=sourceBefore.schema==sourceAfter.schema&&sourceBefore.campaign==sourceAfter.campaign&&
      sourceBefore.fixedLineLayout==sourceAfter.fixedLineLayout&&sourceBefore.counters.size()==sourceAfter.counters.size();
  for(std::size_t i=0;sameSource&&i<sourceBefore.counters.size();++i) {
    const auto& a=sourceBefore.counters[i];const auto& b=sourceAfter.counters[i];
    sameSource=a.channel.kind==b.channel.kind&&a.channel.line==b.channel.line&&
        a.channel.species==b.channel.species&&a.next==b.next;
  }
  if(!sameSource)++failures;
  const bool untouched=particles==PIC::ParticleBuffer::GetAllPartNum()&&
      seed==RandomNumberGenerator::rndLastSeed&&generator==RandomNumberGenerator::gen&&
      distribution==RandomNumberGenerator::dist&&hasNormal==RandomNumberGenerator::NormalSpareAvailable&&
      std::memcmp(&normal,&RandomNumberGenerator::NormalSpare,sizeof(normal))==0;
  if(!untouched)++failures;
  unsigned global=0;MPI_Allreduce(&failures,&global,1,MPI_UNSIGNED,MPI_SUM,MPI_GLOBAL_COMMUNICATOR);
  Result result;result.status=global?Status::Fail:Status::Pass;
  result.message=global?"source capability preflight changed native state or returned an incorrect status":
      "unsupported legacy and unresolved Maxwellian sources reject before allocation, random draws, or source-counter reservation";
  result.configuration={"evidence_level=native-PIC-and-core-RNG-preflight",
      "source_semantics=no-fabricated-line-attachment;destination=srcSEP3D-without-implied-port",
      "coverage=typed-preflight-nonmutation-and-source-counter;host-entry-order=source-inspected"};
  result.metrics.push_back({"assertion_failures",static_cast<double>(global),0,"<=","count"});
  return result;
}
} // namespace
std::vector<Descriptor> SourceNativeDescriptors() {
  Descriptor descriptor;descriptor.id="SOURCE01";descriptor.name="Unsupported source preflight is observational";
  descriptor.group="source-native";descriptor.description=descriptor.name;
  descriptor.initialization=InitializationLevel::FieldLineModel;descriptor.runtime=RuntimeClass::Routine;
  descriptor.supportedBuildModes="native srcSEP initialized PIC model";
  descriptor.seedPolicy="snapshot complete existing native/core RNG state; no draws";
  descriptor.stateIsolation="source enable flag restored; particle counts and complete core RNG unchanged";
  descriptor.callback=RejectionIsObservational;return {descriptor};
}
} } // namespace SEP::Testing
