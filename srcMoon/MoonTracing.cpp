#include "pic.h"
#include "Moon.h"
#include "MoonInput.h"
#include "MoonReactionController.h"
#include "MoonSurfaceEnvironment.h"
#include "MoonTransportClock.h"

#include <cmath>
#include <algorithm>
#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <iomanip>
#include <limits>
#include <sstream>
#include <type_traits>

namespace Moon {
namespace {
bool TracingSelected() {
  return Runtime::HasConfiguration() &&
      Runtime::GetConfiguration().particleTracingMover!="legacy";
}
constexpr std::uint64_t CommonClockMagic=UINT64_C(0x4d4f4f4e54494d31);
struct CommonParticleTime {
  std::uint64_t magic=0,completedSequence=0,targetSequence=0;
  long double relativeEpochS=0,targetEndS=0;
  bool targeting=false;
};
static_assert(std::is_trivially_copyable<CommonParticleTime>::value,
    "common Moon time metadata must survive actual MPI particle-byte copies");
long int CommonTimeOffset=-1;
TransportTime::Clock CommonClock;
double CommonDurationS=0;
bool ClockError(std::string* error,const std::string& message) {
  if(error)*error=message;
  return false;
}
[[noreturn]] void FatalClock(const std::string& message) {
  const std::string complete="Moon common-time ERROR: "+message+
    ". The physical step is not accepted. Use independent-history for unrelated "
    "trajectories, or configure one identical constant actual physical dt for "
    "all species/blocks and explicitly initialize every source at its correct "
    "birth time. A complete clock/source/particle checkpoint is not supported; "
    "restart this configuration from its declared initial state.";
  exit(__LINE__,__FILE__,complete.c_str());std::abort();
}
CommonParticleTime ReadCommon(long int particle) {
  CommonParticleTime value;
  std::memcpy(&value,PIC::ParticleBuffer::GetParticleDataPointer(particle)+CommonTimeOffset,sizeof(value));
  return value;
}
void WriteCommon(long int particle,const CommonParticleTime& value) {
  std::memcpy(PIC::ParticleBuffer::GetParticleDataPointer(particle)+CommonTimeOffset,&value,sizeof(value));
}
template<class Function>void ForOwnedParticle(Function function) {
  const int cells=_BLOCK_CELLS_X_*_BLOCK_CELLS_Y_*_BLOCK_CELLS_Z_;
  for(unsigned index=0;index<PIC::DomainBlockDecomposition::nLocalBlocks;++index) {
    auto node=PIC::DomainBlockDecomposition::BlockTable[index];
    if(!node || !node->block || node->Thread!=PIC::ThisThread)continue;
    for(int cell=0;cell<cells;++cell)
      for(long int particle=node->block->FirstCellParticleTable[cell];particle!=-1;
          particle=PIC::ParticleBuffer::GetNext(particle))function(particle);
  }
}
bool AllClockRanks(bool local,const std::string& localMessage,std::string* error) {
  int good=local?1:0,all=0;
  MPI_Allreduce(&good,&all,1,MPI_INT,MPI_MIN,MPI_GLOBAL_COMMUNICATOR);
  if(!all)return ClockError(error,localMessage.empty()?
    "another MPI rank rejected the common time invariant; inspect every retained rank diagnostic":localMessage);
  return true;
}
void PrepareCommonMove(long int particle,double remaining) {
  if(!CommonTransportTimeSelected())return;
  const auto* step=CommonClock.step();
  if(CommonTimeOffset<0 || !step)FatalClock("particle motion/source callback occurred outside a declared common Begin/Commit step");
  TransportTime::BirthTime birth;std::string error;
  if(!CommonClock.BirthFromRemaining(remaining,&birth,&error))FatalClock(error);
  const double actual=PIC::Mover::ForceCutSurface::GetParticleEpochS(particle);
  if(!TransportTime::EpochProjectionAgrees(actual,birth.epochEtS))
    FatalClock("particle birth/arrival epoch plus supplied remaining seconds does not target the current common endpoint");
  auto record=ReadCommon(particle);
  if(record.magic!=0 && record.magic!=CommonClockMagic)FatalClock("corrupt common packet time metadata");
  if(record.magic==CommonClockMagic && record.targeting) {
    if(record.targetSequence!=step->sequence || record.targetEndS!=step->relativeEndS)
      FatalClock("continued packet has a stale/different common target; do not reset MPI/reaction histories");
    return;
  }
  if(record.magic==CommonClockMagic && (record.completedSequence!=CommonClock.completedSteps() ||
      record.relativeEpochS!=step->relativeBeginS))
    FatalClock("existing packet was not completed at this common step beginning");
  record.magic=CommonClockMagic;record.targeting=true;
  record.targetSequence=step->sequence;record.targetEndS=step->relativeEndS;
  // A genuinely new source has no prior completed interval. Its actual birth
  // is explicitly constrained above by its remaining physical seconds. All
  // normal MPI/product copies retain this immutable target in their bytes.
  if(record.completedSequence==0)record.relativeEpochS=birth.relativeS;
  WriteCommon(particle,record);
}
}

bool CommonTransportTimeSelected() {
  return TracingSelected() && Runtime::GetConfiguration().tracingClockMode=="common-time";
}
const TransportTime::Clock* GetCommonTransportClock() {
  return CommonTransportTimeSelected() && CommonClock.initialized()?&CommonClock:nullptr;
}
bool InitializeCommonTransportTime(double origin,std::string* error) {
  if(!error)return false;
  if(!CommonTransportTimeSelected()){error->clear();return true;}
#if defined(_NO_SPICE_CALLS_) || _EXOSPHERE__ORBIT_CALCUALTION__MODE_ != _PIC_MODE_ON_
  return ClockError(error,"normal common-time initialization requires the configured SPICE UTC-to-ET origin and orbit-enabled frozen external state; no implicit J2000/time origin is supported");
#endif
  if(CommonTimeOffset<0 || PIC::ParticleBuffer::ParticleDataBuffer!=nullptr)
    return ClockError(error,"initialize explicit common time after tracing storage registration and before particle-buffer allocation");
  if(_PIC_READ_PARTICLE_DATA_RESTART_FILE__MODE_!=_PIC_READ_PARTICLE_DATA_RESTART_FILE__MODE_OFF_ ||
      _PIC_RECOVER_SAMPLING_DATA_RESTART_FILE__MODE_!=_PIC_RECOVER_SAMPLING_DATA_RESTART_FILE__MODE_OFF_ ||
      _PIC_AUTOSAVE_PARTICLE_DATA_RESTART_FILE__MODE_!=_PIC_AUTOSAVE_PARTICLE_DATA_RESTART_FILE__MODE_OFF_ ||
      PIC::Restart::LoadRestartSWMF)
    return ClockError(error,"common-time clock/source metadata has no qualified checkpoint protocol; disable particle/sampling restart read, autosave and SWMF restart");
  if(!CommonClock.Initialize(origin,error))return false;
  // This is the FROZEN coefficient epoch, before initial source-weight/CDF
  // queries. It is not advanced when the transport clock later commits.
  SurfaceEnvironment::BeginProductionParticleStep(origin);
  return true;
}
bool BeginCommonTransportStep(double externalEpoch,std::string* error) {
  if(!error)return false;
  if(!CommonTransportTimeSelected()){error->clear();return true;}
  // Normal PIC refreshes this table at TimeStep entry. Our pre-sampling
  // admission runs just before that entry, so explicitly refresh it after a
  // preceding load balance instead of inspecting stale ownership pointers.
  PIC::DomainBlockDecomposition::UpdateBlockTable();
  bool valid=CommonClock.initialized() && !CommonClock.active() && CommonTimeOffset>=0 &&
    std::isfinite(externalEpoch) && externalEpoch==CommonClock.originEtS();
  std::string problem=valid?"":"common Begin requires initialized idle clock and the unchanged explicit frozen external epoch";
  const long double localClock[2]={static_cast<long double>(CommonClock.originEtS()),CommonClock.elapsedS()};
  long double clockMinimum[2],clockMaximum[2];
  MPI_Allreduce(localClock,clockMinimum,2,MPI_LONG_DOUBLE,MPI_MIN,MPI_GLOBAL_COMMUNICATOR);
  MPI_Allreduce(localClock,clockMaximum,2,MPI_LONG_DOUBLE,MPI_MAX,MPI_GLOBAL_COMMUNICATOR);
  if(clockMinimum[0]!=clockMaximum[0] || clockMinimum[1]!=clockMaximum[1]) {
    valid=false;problem="MPI ranks disagree on the explicit common time origin or accepted relative elapsed time";
  }
  double localMinimum=std::numeric_limits<double>::max(),localMaximum=0;
  long int localBlocks=0;
  for(unsigned index=0;index<PIC::DomainBlockDecomposition::nLocalBlocks;++index) {
    auto node=PIC::DomainBlockDecomposition::BlockTable[index];
    if(!node || !node->block || node->Thread!=PIC::ThisThread)continue;
    ++localBlocks;
    for(int species=0;species<PIC::nTotalSpecies;++species) {
      const double dt=node->block->GetLocalTimeStep(species);
      if(!std::isfinite(dt) || !(dt>0)){valid=false;problem="an actual block/species dt is nonfinite or nonpositive";continue;}
      localMinimum=std::min(localMinimum,dt);localMaximum=std::max(localMaximum,dt);
    }
  }
  double minimum=0,maximum=0;long int blocks=0;
  MPI_Allreduce(&localMinimum,&minimum,1,MPI_DOUBLE,MPI_MIN,MPI_GLOBAL_COMMUNICATOR);
  MPI_Allreduce(&localMaximum,&maximum,1,MPI_DOUBLE,MPI_MAX,MPI_GLOBAL_COMMUNICATOR);
  MPI_Allreduce(&localBlocks,&blocks,1,MPI_LONG,MPI_SUM,MPI_GLOBAL_COMMUNICATOR);
  if(blocks==0 || minimum!=maximum || (CommonDurationS>0 && minimum!=CommonDurationS)) {
    std::ostringstream text;text<<std::setprecision(17)
      <<"common-time requires identical constant ACTUAL dt for every block/species; observed min="<<minimum
      <<" max="<<maximum<<" s, previously accepted dt="<<CommonDurationS<<" s. Mixed local/species clocks and variable-duration iteration averages are not qualified";
    problem=text.str();valid=false;
  }
  if(PIC::Mover::ForceCutSurface::CountPending()!=0){valid=false;problem="Begin found unfinished physical continuations from the preceding step";}
  TransportTime::Clock candidate=CommonClock;TransportTime::Step token;
  if(valid && !candidate.Begin(minimum,&token,&problem))valid=false;
  if(valid)ForOwnedParticle([&](long int particle) {
    const auto record=ReadCommon(particle);const double epoch=PIC::Mover::ForceCutSurface::GetParticleEpochS(particle);
    const bool newInitial=record.magic==0 && CommonClock.completedSteps()==0;
    if(!std::isfinite(epoch) || !TransportTime::EpochProjectionAgrees(epoch,token.beginEtS) ||
        PIC::Mover::ForceCutSurface::GetTimeContext(particle).RemainingTimeS!=0 ||
        (!newInitial && (record.magic!=CommonClockMagic || record.targeting ||
          record.completedSequence!=CommonClock.completedSteps() || record.relativeEpochS!=token.relativeBeginS))) {
      valid=false;problem="existing initial/surviving population is not at the explicit common step beginning; provide a finite qualified t0 and do not mix independent histories";
    }
  });
  if(!AllClockRanks(valid,problem,error))return false;
  CommonClock=candidate;CommonDurationS=minimum;
  // No operation below can allocate or fail after publishing the clock. The
  // initial population is accepted only after proving its EXPLICIT t0 above.
  ForOwnedParticle([&](long int particle) {
    auto record=ReadCommon(particle);
    if(record.magic==0){record.magic=CommonClockMagic;record.relativeEpochS=token.relativeBeginS;
      record.completedSequence=CommonClock.completedSteps();WriteCommon(particle,record);}
  });
  error->clear();return true;
}
bool CommitCommonTransportStep(std::string* error) {
  if(!error)return false;
  if(!CommonTransportTimeSelected()){error->clear();return true;}
  PIC::DomainBlockDecomposition::UpdateBlockTable();
  const auto* step=CommonClock.step();bool valid=step!=nullptr;
  std::string problem=valid?"":"common Commit has no active Begin token";
  if(PIC::Mover::ForceCutSurface::CountPending()!=0){valid=false;problem="Commit found unfinished MPI/AMR/product physical intervals";}
  if(valid)ForOwnedParticle([&](long int particle) {
    const auto record=ReadCommon(particle);const auto time=PIC::Mover::ForceCutSurface::GetTimeContext(particle);
    const double actual=PIC::Mover::ForceCutSurface::GetParticleEpochS(particle);
    if(record.magic!=CommonClockMagic || !record.targeting || record.targetSequence!=step->sequence ||
        record.targetEndS!=step->relativeEndS || time.RemainingTimeS!=0 ||
        !TransportTime::EpochProjectionAgrees(actual,step->endEtS)) {
      valid=false;problem="Commit found an unprepared/stale/incomplete packet or an endpoint outside the derived binary64 projection bound; no provisional state can be sampled as complete";
    }
  });
  TransportTime::Clock candidate=CommonClock;TransportTime::Step token;
  if(valid){token=*step;if(!candidate.Commit(token,&problem))valid=false;}
  if(!AllClockRanks(valid,problem,error))return false;
  // Validation of EVERY rank precedes physical endpoint publication. This is
  // only canonical representation of an already completed interval, never a
  // late flight, reaction/exposure debit, resampling or timestep adjustment.
  ForOwnedParticle([&](long int particle) {
    auto record=ReadCommon(particle);record.completedSequence=token.sequence;
    record.relativeEpochS=token.relativeEndS;record.targeting=false;
    PIC::Mover::ForceCutSurface::SetParticleEpochS(particle,token.endEtS);WriteCommon(particle,record);
  });
  CommonClock=candidate;error->clear();return true;
}

void RegisterParticleTracing() {
  if (!TracingSelected()) return;
  const auto& input=Runtime::GetConfiguration();
  auto& controls=PIC::Mover::ForceCutSurface::Settings;
  controls.RelativeTolerance=input.tracingRelativeTolerance;
  controls.AbsolutePositionToleranceM=input.tracingAbsolutePositionToleranceM;
  controls.AbsoluteVelocityToleranceMPerS=input.tracingAbsoluteVelocityToleranceMS;
  controls.GeometryToleranceM=input.tracingGeometryToleranceM;
  controls.TimeToleranceS=input.tracingTimeToleranceS;
  controls.MinimumSubstepS=input.tracingMinimumSubstepS;
  controls.MaximumSubstepS=input.tracingMaximumSubstepS;
  controls.Adaptive=input.tracingAdaptive;
  controls.MaximumSubsteps=input.tracingMaximumSubsteps;
  controls.MaximumRefinements=input.tracingMaximumRefinements;
  controls.MaximumInteractions=input.tracingMaximumInteractions;
  controls.MaximumSamePointInteractions=input.tracingMaximumSamePointInteractions;
  // Both orders use the configured Moon total-acceleration macro: lunar and
  // differential Earth/Sun gravity, Na radiation pressure/shadows, q/m Lorentz
  // acceleration and the declared rotating-frame terms. No force is copied
  // into the tracing wrapper or tests. Default NULL selects that production
  // callback; external time-dependent adapters belong to a separate contract.
  PIC::Mover::ForceCutSurface::Acceleration=nullptr;
  PIC::Mover::ForceCutSurface::RegisterStorage();
  if(CommonTransportTimeSelected() && CommonTimeOffset<0)
    PIC::ParticleBuffer::RequestDataStorage(CommonTimeOffset,sizeof(CommonParticleTime),alignof(CommonParticleTime));
}

void ConfigureParticleTracingFields() {
  if (!TracingSelected()) return;
  const auto& input=Runtime::GetConfiguration();
  if (_PIC_COUPLER_MODE_!=_PIC_COUPLER_MODE__OFF_)
    exit(__LINE__,__FILE__,"Moon tracing receipt declares prescribed LSO E/B, "
        "but the compiled coupler is ON; disable the coupler or qualify an "
        "explicit time/frame-aware field provider instead of mixing inputs");
  if (input.tracingFieldFrame!="lso" || input.simulationEnvironmentMode!="frozen")
    exit(__LINE__,__FILE__,"Moon's current tracing adapter requires an explicit "
        "stationary LSO field and frozen environment; use a qualified per-particle "
        "ephemeris provider for evolving local-time trajectories");
}

bool GetParticleTracingFields(double* electric,double* magnetic) {
  if (!TracingSelected()) return false;
  const auto& input=Runtime::GetConfiguration();
  // The existing typical B table is intentionally immutable. Select the
  // explicit source through the same production Lorentz branch instead of
  // const-casting a table or silently deriving E from an assumed wind speed.
  for(int k=0;k<3;++k) {
    electric[k]=input.tracingElectricFieldVM[k];
    magnetic[k]=input.tracingMagneticFieldT[k];
  }
  return true;
}

int MoveParticle(long int ptr,double dt,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node) {
  if (!TracingSelected())
    return PIC::Mover::UniformWeight_UniformTimeStep_SecondOrder(ptr,dt,node);
  if (!PIC::Mover::ForceCutSurface::IsStorageRegistered())
    exit(__LINE__,__FILE__,"Moon tracing storage was not registered before particle allocation");
  if (!std::isfinite(PIC::Mover::ForceCutSurface::GetParticleEpochS(ptr))) {
    // A genuinely new source particle has no trajectory epoch. Source
    // injection already publishes its birth offset through this production
    // service; continuing/MPI-cloned particles retain their own stored epoch.
    // No global step or species-ID change resets an existing history.
    PIC::Mover::ForceCutSurface::SetParticleEpochS(ptr,
        SurfaceEnvironment::ProductionInjectionEpochS());
  }
  PrepareCommonMove(ptr,dt);
  // Open exactly the caller's physical interval only after its birth/arrival
  // epoch is known. The opt-in controller retains optical depth through AMR,
  // adaptive retries and MPI ownership changes; those must not draw new event
  // clocks. Installation requires this tracer and a qualified explicit policy.
  if (Chemistry::Controller::Enabled())
    Chemistry::Controller::PrepareParticle(ptr,dt,node);
  if (Runtime::GetConfiguration().particleTracingMover=="force-cut-surface-first-order")
    return PIC::Mover::ForceCutSurfaceFirstOrder(ptr,dt,node);
  return PIC::Mover::ForceCutSurfaceSecondOrder(ptr,dt,node);
}

int MoveParticleBoundaryInjection(long int ptr,double remaining,
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node,bool firstBoundary,
    CutCell::cTriangleFace* injectionFace) {
  if (!TracingSelected())
    return PIC::Mover::UniformWeight_UniformTimeStep_noForce_TraceTrajectory_BoundaryInjection_SecondOrder(
        ptr,remaining,node,firstBoundary,injectionFace);
  if (!PIC::Mover::ForceCutSurface::IsStorageRegistered() || !node || !node->block ||
      !std::isfinite(remaining) || remaining<0)
    exit(__LINE__,__FILE__,"Moon tracing boundary birth requires registered storage, an actual block and finite remaining physical seconds");
  if(CommonTransportTimeSelected()) {
    const auto* clock=GetCommonTransportClock();TransportTime::BirthTime birth;std::string error;
    if(!clock || !clock->active())FatalClock("source boundary callback occurred outside the common physical step");
    const double offset=PIC::Mover::MovingSurfaceInjectionTimeOffsetS;
    const bool success=std::isfinite(offset)?clock->BirthFromOffset(offset,&birth,&error):
      clock->BirthFromRemaining(remaining,&birth,&error);
    if(!success || (std::isfinite(offset) && birth.remainingS!=remaining))
      FatalClock("source departure/remaining interval does not complement the declared common duration: "+error);
    const double existing=PIC::Mover::ForceCutSurface::GetParticleEpochS(ptr);
    if(std::isfinite(existing) && !TransportTime::EpochProjectionAgrees(existing,birth.epochEtS))
      FatalClock("explicit source birth epoch conflicts with its common-step departure/remaining interval");
    if(!std::isfinite(existing))PIC::Mover::ForceCutSurface::SetParticleEpochS(ptr,birth.epochEtS);
  }
  else if (!std::isfinite(PIC::Mover::ForceCutSurface::GetParticleEpochS(ptr))) {
    // Moving-terrain injection supplies its actual Poisson departure offset.
    // The stationary generic source instead passes an explicitly sampled
    // remaining interval, so its birth is step-start + (species dt-remaining).
    // Do not draw a second offset or label that birth as the step beginning.
    double epoch=SurfaceEnvironment::ProductionInjectionEpochS();
    if (!std::isfinite(PIC::Mover::MovingSurfaceInjectionTimeOffsetS)) {
      const int species=PIC::ParticleBuffer::GetI(ptr);
      if (species<0 || species>=PIC::nTotalSpecies)
        exit(__LINE__,__FILE__,"Moon tracing boundary birth has an invalid compiled species");
      const double duration=node->block->GetLocalTimeStep(species);
      if (!std::isfinite(duration) || duration<=0 || remaining>duration)
        exit(__LINE__,__FILE__,"Moon tracing stationary source needs remaining seconds within its actual local species step; pass a qualified explicit birth epoch for another source contract");
      epoch+=duration-remaining;
    }
    PIC::Mover::ForceCutSurface::SetParticleEpochS(ptr,epoch);
  }
  // The new tracer classifies zero-time contacts by the actual incoming
  // velocity and resolves curved re-entry. It must not skip every internal
  // surface over the first segment, which could hide a gravity-driven return.
  // The legacy-only flag/face are forwarded unchanged in the branch above.
  return MoveParticle(ptr,remaining,node);
}
} // namespace Moon
