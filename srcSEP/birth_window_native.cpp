#include "transport_common.h"
#include "util/sep_background_runtime.h"
#include <climits>
#include <cmath>
#include <cstring>
#include <exception>
#include <set>
#include <utility>
#if _PIC_COUPLER_MODE_ == _PIC_COUPLER_MODE__SWMF_
#include "amps2swmf.h"
#endif

namespace SEP { namespace Transport { namespace PICAdapter { namespace {
// The host freezes this immutable value before worker creation and clears it
// only after every worker and migration has joined. Particle threads only read
// it. No source/operator may replace its epoch from the mutable PIC counter.
Birth::DriverWindow gDriverWindow;
bool gDriverWindowActive=false;
std::shared_ptr<const Background::BackgroundSnapshot> gDriverSnapshot;
Status Invalid(const char* message) {
  return Status::Error(StatusCode::InvalidParticleState,message);
}
Status CollectiveStatus(const Status& local) {
  // Every host rank participates even if one local preflight fails. Broadcast
  // the first actual status so no rank proceeds to RNG/particle mutation while
  // a peer has already rejected the source clock.
  const int candidate=local.ok()?INT_MAX:PIC::ThisThread;int first=INT_MAX;
  MPI_Allreduce(&candidate,&first,1,MPI_INT,MPI_MIN,MPI_GLOBAL_COMMUNICATOR);
  if(first==INT_MAX)return Status::Ok();
  std::string message=PIC::ThisThread==first?local.message:std::string();
  int length=static_cast<int>(message.size()),code=static_cast<int>(local.code);
  MPI_Bcast(&length,1,MPI_INT,first,MPI_GLOBAL_COMMUNICATOR);
  MPI_Bcast(&code,1,MPI_INT,first,MPI_GLOBAL_COMMUNICATOR);
  message.resize(length);MPI_Bcast(length?&message[0]:nullptr,length,MPI_CHAR,first,MPI_GLOBAL_COMMUNICATOR);
  return Status::Error(static_cast<StatusCode>(code),message);
}
Status JoinedBirthsComplete() {
  // Completed windows contain only live, listed field-line records. Inspect
  // both ownership bookkeeping and the persisted schema before accepting a
  // pending flag; uninitialized bytes cannot masquerade as consumed movement.
  namespace FL=PIC::FieldLine;namespace PB=PIC::ParticleBuffer;
  std::set<long int> seen;
  for(int line=0;line<FL::nFieldLine;++line)
    for(auto* segment=FL::FieldLinesAll[line].GetFirstSegment();segment;segment=segment->GetNext()) {
      if(segment->tempFirstParticleIndex.load()!=-1)
        return Invalid("driver birth window cannot close before native moving lists have joined");
      for(auto particle=segment->FirstParticleIndex;particle!=-1;particle=PB::GetNext(particle)) {
        if(particle<0||particle>=PB::GetMaxNPart()||!PB::IsParticleAllocated(particle)||!seen.insert(particle).second)
          return Invalid("driver birth window encountered an invalid or repeated native particle");
        auto* data=PB::GetParticleDataPointer(particle);std::uint64_t schema=0;
        std::memcpy(&schema,data+Offset::TransportSchema,sizeof(schema));
        if(schema!=ParticleTransportSchema)return Invalid("driver birth window encountered incompatible particle history");
        Birth::Record birth;std::memcpy(&birth,data+Offset::BirthRecord,sizeof(birth));
        const auto status=Birth::ValidateRecord(birth);if(!status.ok())return status;
        if(birth.pendingFirstMove)return Invalid("driver window has an unconsumed scheduled birth; remaining movement cannot be discarded");
      }
    }
  if(seen.size()!=static_cast<std::size_t>(PB::GetAllPartNum()))
    return Invalid("driver birth window refuses unlisted or Cartesian particles");
  return Status::Ok();
}
} // namespace

bool NativeDriverWindowActive(){return gDriverWindowActive;}
Status GetNativeDriverWindow(Birth::DriverWindow* window) {
  if(!window)return Status::Error(StatusCode::InvalidArgument,"native driver window output is required");
  if(!gDriverWindowActive)return Status::Error(StatusCode::UnsupportedConfiguration,
      "scheduled particle birth requires an explicit host window frozen before PIC::TimeStep; the current PIC counter cannot infer it");
  *window=gDriverWindow;return Status::Ok();
}
Status MakeStepStartBirthRecord(Birth::Record* birth) {
  Birth::DriverWindow window;const auto status=GetNativeDriverWindow(&window);
  return status.ok()?Birth::StepStartBirth(window,birth):status;
}
Status MakeUniformBirthRecord(double phaseOpen01,Birth::Record* birth) {
  Birth::DriverWindow window;const auto status=GetNativeDriverWindow(&window);
  return status.ok()?Birth::UniformBirth(window,phaseOpen01,birth):status;
}
Status ValidateNativeDriverSnapshot(const Background::BackgroundSnapshot& snapshot) {
  if(gDriverWindowActive && &snapshot!=gDriverSnapshot.get())
    return Invalid("particle background differs from the immutable realization frozen with its source driver window");
  return Status::Ok();
}
Status NativeDriverClockCapability() {
  Birth::HostClockPolicy policy;policy.speciesCount=PIC::nTotalSpecies;
  policy.counterAtEnd=_PIC_GLOBAL_TIME_COUNTER_MODE_==_PIC_MODE_ON_;
  policy.externallyOwned=_PIC_COUPLER_MODE_==_PIC_COUPLER_MODE__SWMF_;
#if _SIMULATION_TIME_STEP_MODE_ == _SINGLE_GLOBAL_TIME_STEP_
  policy.mode=Birth::HostStepMode::SingleGlobal;
#elif _SIMULATION_TIME_STEP_MODE_ == _SPECIES_DEPENDENT_GLOBAL_TIME_STEP_
  policy.mode=Birth::HostStepMode::SpeciesGlobal;
#else
  return Status::Error(StatusCode::UnsupportedConfiguration,
      "native source birth clock currently requires a declared global timestep");
#endif
#if _PIC_COUPLER_MODE_ == _PIC_COUPLER_MODE__SWMF_
  // Merely compiling SWMF is not a clock capability. Accept only this step's
  // still-planned external interval, before its mutable counter changes phase.
  // The value-layer validator checks the audited macro policy, physical
  // departure and full nominal duration; no scheduler endpoint is inferred.
  if(AMPS2SWMF::CoupledStep.Configured()&&AMPS2SWMF::CoupledStep.Planned()) {
    const auto& host=AMPS2SWMF::CoupledStep.Current();
    Birth::DriverWindow bound;bound.schema=1;bound.startEpochS=host.startEpochS;
    bound.durationS=host.durationS;bound.endEpochS=host.endEpochS;
    return Birth::ValidateExternalHostWindow(policy,bound,
        Background::SimulationTimeSeconds(),PIC::ParticleWeightTimeStep::GlobalTimeStep[0]);
  }
#endif
  return Birth::ValidateHostClockPolicy(policy);
}
Status FreezeNativeDriverWindow() {
  Status local=NativeDriverClockCapability();Birth::DriverWindow candidate;
  std::shared_ptr<const Background::BackgroundSnapshot> candidateSnapshot;
  try {
    if(gDriverWindowActive)local=Invalid("cannot replace an active native source driver window");
#if _PIC_COUPLER_MODE_ == _PIC_COUPLER_MODE__SWMF_
    // The generic PT particle restart does not restore this physical ledger
    // and its matching source/background/sampling interval. Do not load it
    // inside PIC after freezing a fresh host window under the scheduler epoch.
    // An SC/IH restart with a newly initialized PT population is unaffected.
    if(local.ok()&&PIC::Restart::LoadRestartSWMF)
      local=Status::Error(StatusCode::UnsupportedConfiguration,
          "coupled PT restart requires restored physical clock, background and sampling state before source-window freeze");
#endif
#if _SIMULATION_TIME_STEP_MODE_ != _SINGLE_GLOBAL_TIME_STEP_ && _SIMULATION_TIME_STEP_MODE_ != _SPECIES_DEPENDENT_GLOBAL_TIME_STEP_
    if(local.ok())local=Status::Error(StatusCode::UnsupportedConfiguration,
        "native source birth clock currently requires a declared global timestep");
#endif
    if(local.ok()&&!PIC::ParticleWeightTimeStep::GlobalTimeStepInitialized)
      local=Invalid("native driver window requires initialized global timestep storage");
    const double dtS=PIC::ParticleWeightTimeStep::GlobalTimeStep[0];
#if _SIMULATION_TIME_STEP_MODE_ == _SPECIES_DEPENDENT_GLOBAL_TIME_STEP_
    if(local.ok())for(int species=0;species<PIC::nTotalSpecies;++species)
      if(PIC::ParticleWeightTimeStep::GlobalTimeStep[species]!=dtS) {
        local=Status::Error(StatusCode::UnsupportedConfiguration,
            "unequal species-global timesteps have no declared common source/background birth clock");break;
      }
#endif
    if(local.ok())local=Birth::FreezeWindow(Background::SimulationTimeSeconds(),dtS,&candidate);
    const auto snapshot=Background::SnapshotStore::Instance().Current();
    if(local.ok()&&(!snapshot||!snapshot->Covers(candidate.startEpochS)||!snapshot->Covers(candidate.endEpochS)))
      local=Invalid("one immutable background snapshot must cover the complete native source driver window");
    if(local.ok())Background::SnapshotStore::Instance().AssertProviderMayWrite(snapshot->provider());
    if(local.ok())candidateSnapshot=snapshot;
  }catch(const std::exception& error){local=Status::Error(StatusCode::InvalidParticleState,error.what());}
  auto status=CollectiveStatus(local);if(!status.ok())return status;
  const double values[3]={candidate.startEpochS,candidate.durationS,candidate.endEpochS};
  double minimum[3],maximum[3];
  MPI_Allreduce(values,minimum,3,MPI_DOUBLE,MPI_MIN,MPI_GLOBAL_COMMUNICATOR);
  MPI_Allreduce(values,maximum,3,MPI_DOUBLE,MPI_MAX,MPI_GLOBAL_COMMUNICATOR);
  for(int i=0;i<3;++i)if(minimum[i]!=maximum[i])
    return Invalid("MPI ranks disagree on the authoritative native source driver window");
  gDriverWindow=candidate;gDriverSnapshot=std::move(candidateSnapshot);
  gDriverWindowActive=true;return Status::Ok();
}
Status CloseNativeDriverWindow() {
  // This is observational until all-rank agreement. Reader/clock/pending-state
  // failures leave the same frozen snapshot/window active for diagnosis/retry.
  Status local=Status::Ok();
  try {
    if(!gDriverWindowActive)local=Invalid("cannot close a native source window that was not frozen");
    if(local.ok()) {
      const auto snapshot=Background::SnapshotStore::Instance().Current();
      if(snapshot!=gDriverSnapshot)local=Invalid("background realization changed before native birth-window closure");
      else Background::SnapshotStore::Instance().AssertProviderMayWrite(snapshot->provider());
    }
    if(local.ok()&&Background::SimulationTimeSeconds()!=gDriverWindow.endEpochS)
      local=Invalid("native driver clock did not reach its frozen end epoch; this is not a completed source timestep");
    if(local.ok())local=JoinedBirthsComplete();
  }catch(const std::exception& error){local=Status::Error(StatusCode::InvalidParticleState,error.what());}
  const auto status=CollectiveStatus(local);if(!status.ok())return status;
  gDriverWindowActive=false;gDriverSnapshot.reset();return Status::Ok();
}
} } }
