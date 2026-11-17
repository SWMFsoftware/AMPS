#include "sep_host_step_native_validation.h"
#include "../sep.h"
#include "../transport_common.h"
#include "../physics_checkpoint_state.h"
#include "sep_analytic_radial_background.h"
#include "sep_background_runtime.h"
#include "sep_initialization.h"
#include "sep_run_configuration.h"
#include "sep_shock_source_core.h"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <limits>
#include <mutex>
#include <set>
#include <sstream>
#include <string>
#include <vector>

void amps_time_step();

namespace SEP { namespace Testing { namespace {
namespace A=Transport::PICAdapter;
namespace B=Transport::Birth;
namespace PB=PIC::ParticleBuffer;
namespace FL=PIC::FieldLine;
namespace S=ParticleSource::ShockWave::Tenishev2005;
namespace W=AlfvenTurbulence_Kolmogorov;

struct Input {
  unsigned steps=0,sourceMacros=0;
  double sourceEnergyMeV=0,sourceMu=0,constantDmumuPerS=0;
};
Transport::Status AgreeStatus(const Transport::Status& local) {
  int candidate=local.ok()?std::numeric_limits<int>::max():PIC::ThisThread,first=0;
  MPI_Allreduce(&candidate,&first,1,MPI_INT,MPI_MIN,MPI_GLOBAL_COMMUNICATOR);
  if(first==std::numeric_limits<int>::max())return Transport::Status::Ok();
  int code=static_cast<int>(local.code);char message[2048]={};
  if(PIC::ThisThread==first)std::snprintf(message,sizeof(message),"host fixture rank %d: %s",first,local.message.c_str());
  MPI_Bcast(&code,1,MPI_INT,first,MPI_GLOBAL_COMMUNICATOR);
  MPI_Bcast(message,static_cast<int>(sizeof(message)),MPI_CHAR,first,MPI_GLOBAL_COMMUNICATOR);
  return Transport::Status::Error(static_cast<Transport::StatusCode>(code),message);
}
Transport::Status ReadInput(Input* output) {
  const auto& context=GetExecutionContext();
  std::ifstream input(context.inputPath.c_str());std::string header;
  if(!input||!std::getline(input,header)||header!="SEP_NATIVE_HOST_STEP 1")
    return Transport::Status::Error(Transport::StatusCode::InvalidArgument,"host-step fixture requires its explicit version-one input protocol");
  Input candidate;std::set<std::string> found;std::string key;
  while(input>>key) {
    if(!found.insert(key).second)return Transport::Status::Error(Transport::StatusCode::InvalidArgument,"duplicate host-step fixture input");
    if(key=="steps")input>>candidate.steps;
    else if(key=="source_macros")input>>candidate.sourceMacros;
    else if(key=="source_energy_mev")input>>candidate.sourceEnergyMeV;
    else if(key=="source_mu")input>>candidate.sourceMu;
    else if(key=="constant_dmumu_per_s")input>>candidate.constantDmumuPerS;
    else return Transport::Status::Error(Transport::StatusCode::InvalidArgument,"unknown host-step fixture input");
    if(!input)return Transport::Status::Error(Transport::StatusCode::InvalidArgument,"malformed host-step fixture input value");
  }
  if(found.size()!=5||!input.eof()||candidate.steps<1||candidate.steps>2||
      candidate.sourceMacros<1||candidate.sourceMacros>8||
      !std::isfinite(candidate.sourceEnergyMeV)||candidate.sourceEnergyMeV<=0||
      !std::isfinite(candidate.sourceMu)||candidate.sourceMu<=0||candidate.sourceMu>1||
      !std::isfinite(candidate.constantDmumuPerS)||candidate.constantDmumuPerS<=0)
    return Transport::Status::Error(Transport::StatusCode::InvalidArgument,"host-step fixture needs all explicit bounded source/numerical controls");
  *output=candidate;return Transport::Status::Ok();
}

struct Trace {
  unsigned step=0;std::uint64_t particle=0;
  double clockS=0,startS=0,endS=0,remainingS=0,requestedS=0,radiusM=0;
  bool injection=false,pending=false,cartesianMatches=true,amrBlockPresent=true;
};
struct Probe {
  std::mutex mutex;
  unsigned step=0,failures=0,sourceCalls=0,firstMoves=0,completedSteps=0;
  unsigned requestedSteps=0,traceSequence=0;
  double launchRadiusM=0;
  std::string message;
  std::vector<Trace> trace;
  std::vector<std::string> traceArtifacts;
  long int (*source)()=nullptr;
  PIC::Mover::fSpeciesDependentParticleMover mover=nullptr;
  void Require(bool condition,const std::string& reason) {
    if(condition)return;
    std::lock_guard<std::mutex> lock(mutex);++failures;
    if(message.size()<2400){if(!message.empty())message+="; ";message+=reason;}
  }
  void Append(const Trace& event){std::lock_guard<std::mutex> lock(mutex);trace.push_back(event);}
};
Probe* activeProbe=nullptr;

Trace ParticleTrace(long particle,double requestedS,bool injection) {
  Trace event;event.step=activeProbe->step;event.injection=injection;event.requestedS=requestedS;
  event.clockS=Background::SimulationTimeSeconds();
  auto* data=PB::GetParticleDataPointer(particle);
  B::Record birth;std::memcpy(&birth,data+Offset::BirthRecord,sizeof(birth));
  std::memcpy(&event.particle,data+Offset::StableParticleId,sizeof(event.particle));
  event.pending=birth.pendingFirstMove!=0;event.startS=birth.parent.startEpochS;
  event.endS=birth.parent.endEpochS;event.remainingS=birth.firstMoveDurationS;
  const int line=PB::GetFieldLineId(data);const double coordinate=PB::GetFieldLineCoord(data);
  double canonical[3];FL::FieldLinesAll[line].GetCartesian(canonical,coordinate);
  event.radiusM=std::hypot(std::hypot(canonical[0],canonical[1]),canonical[2]);
  const double* stored=PB::GetX(data);
  for(unsigned dimension=0;dimension<3;++dimension)
    event.cartesianMatches=event.cartesianMatches&&stored[dimension]==canonical[dimension];
  auto* node=PIC::Mesh::mesh->findTreeNode(canonical,nullptr);
  event.amrBlockPresent=node&&node->block;
  return event;
}

void PersistTrace(bool complete=false) {
  // Each native artifact scope gives this rank its own trusted directory.
  // Save after source dispatch, before the core mover can fail, so a genuine
  // PIC fatal exit still leaves clock/source/AMR-ownership evidence. This is
  // diagnostic output only; no physical clock/state is changed to obtain it.
  // Output::WriteTransactional deliberately forbids overwrites. Publish each
  // cumulative snapshot under a new monotone name; retain the early pre-mover
  // snapshot after fatal exits, while only a complete last snapshot can satisfy
  // the successful campaign reader. No mutable latest-pointer is introduced.
  const auto& context=GetExecutionContext();
  std::ostringstream text;text<<std::setprecision(17);
  text<<"{\"schema\":\"srcsep-native-host-step-trace-v2\",\"case_id\":\""<<context.testId
      <<"\",\"mover\":\""<<Mover::Describe(Mover::CurrentProductionMover()).canonicalName
      <<"\",\"rank\":"<<PIC::ThisThread
      <<",\"snapshot_sequence\":"<<activeProbe->traceSequence
      <<",\"complete\":"<<(complete?"true":"false")
      <<",\"requested_host_steps\":"<<activeProbe->requestedSteps
      <<",\"completed_host_steps\":"<<activeProbe->completedSteps
      <<",\"counter_at_end\":"<<(_PIC_GLOBAL_TIME_COUNTER_MODE_==_PIC_MODE_ON_?"true":"false")
      <<",\"events\":[";
  bool comma=false;
  for(const auto& event:activeProbe->trace) {
    if(comma)text<<',';comma=true;
    text<<"{\"step\":"<<event.step<<",\"particle\":"<<event.particle
        <<",\"phase\":\""<<(event.injection?"source":"mover")<<"\",\"clock_s\":"<<event.clockS
        <<",\"parent_start_s\":"<<event.startS<<",\"parent_end_s\":"<<event.endS
        <<",\"first_duration_s\":"<<event.remainingS<<",\"requested_s\":"<<event.requestedS
        <<",\"source_radius_m\":"<<event.radiusM
        <<",\"pending_first\":"<<(event.pending?"true":"false")
        <<",\"cartesian_matches\":"<<(event.cartesianMatches?"true":"false")
        <<",\"canonical_amr_block_local\":"<<(event.amrBlockPresent?"true":"false")<<'}';
  }
  text<<"]}\n";
  Output::ArtifactMetadata metadata;metadata.schemaVersion="srcsep-native-host-step-trace-v2";
  metadata.recordCount=activeProbe->trace.size();metadata.configurationFingerprint=Background::CurrentConfigurationFingerprint();
  std::ostringstream filename;filename<<"host-step-trace-"<<std::setfill('0')<<std::setw(6)<<activeProbe->traceSequence++<<".json";
  const auto written=Output::WriteTransactional(context.artifactDirectory,filename.str(),text.str(),metadata);
  activeProbe->Require(written.status.ok(),"host trace write failed: "+written.status.message);
  if(written.status.ok())activeProbe->traceArtifacts.push_back(context.artifactDirectory+"/"+filename.str());
}

long ObservedInjection() {
  auto& probe=*activeProbe;++probe.sourceCalls;
  B::DriverWindow window;const auto frozen=A::GetNativeDriverWindow(&window);
  probe.Require(frozen.ok(),"real PIC injection has no frozen native driver window");
  const long injected=probe.source();
  for(int line=0;line<FL::nFieldLine;++line)
    for(auto* segment=FL::FieldLinesAll[line].GetFirstSegment();segment;segment=segment->GetNext())
      for(long particle=segment->FirstParticleIndex;particle!=-1;particle=PB::GetNext(particle)) {
        const auto event=ParticleTrace(particle,0,true);
        if(event.pending) {
          probe.Require(event.startS==window.startEpochS&&event.endS==window.endEpochS&&
                        event.remainingS==window.durationS,"actual source birth differs from the host's frozen step-start interval");
          probe.Require(event.cartesianMatches,"actual source stored PB::X differs from its declared field-line coordinate before first movement");
          const double expectedRadius=probe.launchRadiusM+1800000.0*window.startEpochS;
          probe.Require(std::fabs(event.radiusM-expectedRadius)<=2e-12*expectedRadius,
                        "actual analytic shock birth geometry uses a different epoch from its frozen source departure");
          probe.Append(event);
        }
      }
  PersistTrace();return injected;
}
int ObservedMover(long particle,double dt,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node) {
  auto& probe=*activeProbe;
  const auto event=ParticleTrace(particle,dt,false);
  if(event.pending) {
    probe.Require(event.remainingS==dt,"real PIC mover received a different first-birth duration from its pending record");
    {std::lock_guard<std::mutex> lock(probe.mutex);++probe.firstMoves;}
  }
  probe.Append(event);
  const int result=probe.mover(particle,dt,node);
  if(PB::IsParticleAllocated(particle)) {
    B::Record after;std::memcpy(&after,PB::GetParticleDataPointer(particle)+Offset::BirthRecord,sizeof(after));
    probe.Require(after.pendingFirstMove==0,"actual PIC mover left a surviving newborn pending");
  }
  return result;
}

void DeleteOwnedParticles() {
  // Entry requires an empty population, so all subsequent listed particles are
  // this process fixture's. Teardown removes them before the registry restores
  // its source/RNG metadata; PIC private timestep state is isolated by process
  // termination, not falsely described as restored by this cleanup.
  std::set<long> particles;
  for(int line=0;line<FL::nFieldLine;++line)
    for(auto* segment=FL::FieldLinesAll[line].GetFirstSegment();segment;segment=segment->GetNext()) {
      for(const auto head:{segment->FirstParticleIndex,segment->tempFirstParticleIndex.load()})
        for(long particle=head;particle!=-1;) {
          if(!particles.insert(particle).second)break;
          particle=PB::GetNext(particle);
        }
      segment->FirstParticleIndex=-1;segment->tempFirstParticleIndex=-1;
    }
  for(long particle:particles)PB::DeleteParticle(particle);
}

Result RunHost(Mover::ProductionMover selected) {
  Result result;const auto& context=GetExecutionContext();
  if(context.inputPath.empty()) {
    result.status=Status::Skip;result.message="requires the explicit isolated host-step protocol; use test/run_native_host_step_tests.py";return result;
  }
  Input input;const auto parsed=AgreeStatus(ReadInput(&input));
  if(!parsed.ok()){result.status=Status::Error;result.message=parsed.message;return result;}
  const double controls[5]={static_cast<double>(input.steps),static_cast<double>(input.sourceMacros),
      input.sourceEnergyMeV,input.sourceMu,input.constantDmumuPerS};
  double minimum[5],maximum[5];
  MPI_Allreduce(controls,minimum,5,MPI_DOUBLE,MPI_MIN,MPI_GLOBAL_COMMUNICATOR);
  MPI_Allreduce(controls,maximum,5,MPI_DOUBLE,MPI_MAX,MPI_GLOBAL_COMMUNICATOR);
  for(unsigned i=0;i<5;++i)if(minimum[i]!=maximum[i]) {
    result.status=Status::Error;result.message="host fixture input values differ between MPI ranks";return result;
  }
  const auto wind=Background::AnalyticRadial::DeclaredNativeWind();
  const bool local=Background::ConfiguredProvider()==Background::Provider::Analytic&&
      ShockModelType==cShockModelType::Analytic1D&&PB::GetAllPartNum()==0&&
      FL::nFieldLine==1&&PIC::nTotalSpecies==1&&wind.originM==std::array<double,3>{{0,0,0}}&&
      A::CurrentRunPurpose()==A::RunPurpose::ControlledDiagnostic&&!A::NativeDriverWindowActive()&&
      PIC::BC::UserDefinedParticleInjectionFunction&&PIC::Mover::UserDefinedParticleMover;
  const int readyLocal=local?1:0;int ready=0;
  MPI_Allreduce(&readyLocal,&ready,1,MPI_INT,MPI_MIN,MPI_GLOBAL_COMMUNICATOR);
  if(!ready){result.status=Status::Error;result.message="host fixture requires an empty, fixed zero-origin analytic one-line/one-species diagnostic process";return result;}
  Probe probe;probe.trace.reserve(1024);probe.requestedSteps=input.steps;
  const auto originalRun=Run::Active().get();const auto originalCoefficients=Transport::Coefficient::ActiveConfiguration();
  const auto originalTurbulence=Turbulence::ActiveConfiguration();
  const auto originalMover=Mover::CurrentProductionMover();
  const bool originalActive=W::ActiveFlag,originalCoupling=W::ParticleCouplingMode;
  const int originalLocation=FieldLine::InjectionParameters::InjectLocation;
  const int originalMomentum=FieldLine::InjectionParameters::InjectionMomentumModel;
  const int originalAnalytic=FieldLine::InjectionParameters::UseAnalyticShockModel;
  const int originalMacros=FieldLine::InjectionParameters::nParticlesPerIteration;
  const double originalEnergy=FieldLine::InjectionParameters::ConstEnergyInjectionValue;
  const double originalMu=FieldLine::InjectionParameters::ConstMuInjectionValue;
  probe.source=PIC::BC::UserDefinedParticleInjectionFunction;probe.mover=PIC::Mover::UserDefinedParticleMover;
  // Only a dedicated process may enter this function. These public controls
  // are restored; irreversible PIC::TimeStep statistics/first-call flags are
  // contained by the enforced process boundary and remain explicitly so.
  auto controlled=originalRun;controlled.mover.value=selected;
  controlled.coefficients.source=Transport::Coefficient::SourceMode::Prescribed;
  controlled.coefficients.pitchAngle=Transport::Coefficient::PitchAngleKind::Constant;
  controlled.coefficients.constantDmumuPerS=input.constantDmumuPerS;
  controlled.coefficients.spatial=Transport::Coefficient::SpatialKind::FromPitchAngle;
  controlled.coefficients.meanFreePath=Transport::Coefficient::MeanFreePathKind::FromSpatial;
  controlled.turbulence.source=Turbulence::Source::Prescribed;
  controlled.turbulence.coupling=Turbulence::CouplingPolicy::Disabled;
  controlled.turbulence.advectionEnabled=false;controlled.turbulence.reflectionEnabled=false;
  controlled.turbulence.cascadeEnabled=false;controlled.turbulence.shockInjectionEnabled=false;
  controlled.mergeMinimum.value=0;controlled.mergeMaximum.value=100000;
  controlled.populationControl.minimumParticlesPerCell=0;controlled.populationControl.maximumParticlesPerCell=100000;
  const auto installed=AgreeStatus(Run::Validate(controlled));
  if(!installed.ok()){result.status=Status::Error;result.message=installed.message;return result;}
  struct Cleanup {
    Probe& probe;Run::Configuration run;Transport::Coefficient::Configuration coefficients;
    Turbulence::Configuration turbulence;Mover::ProductionMover mover;
    bool active,coupling;int location,momentum,analytic,macros;double energy,mu;bool restored=false;
    void Restore() {
      if(restored)return;restored=true;
      PIC::BC::UserDefinedParticleInjectionFunction=probe.source;PIC::Mover::SetUserDefinedParticleMover(probe.mover);
      DeleteOwnedParticles();activeProbe=nullptr;Mover::SelectProductionMover(mover);
      FieldLine::InjectionParameters::InjectLocation=location;FieldLine::InjectionParameters::InjectionMomentumModel=momentum;
      FieldLine::InjectionParameters::UseAnalyticShockModel=analytic;
      FieldLine::InjectionParameters::nParticlesPerIteration=macros;FieldLine::InjectionParameters::ConstEnergyInjectionValue=energy;
      FieldLine::InjectionParameters::ConstMuInjectionValue=mu;W::ActiveFlag=active;W::ParticleCouplingMode=coupling;
      const auto coefficientStatus=Transport::Coefficient::SetActiveConfiguration(coefficients);
      const auto turbulenceStatus=Turbulence::SetActiveConfiguration(turbulence);
      const auto runStatus=Run::InstallActive(run);
      probe.Require(coefficientStatus.ok()&&turbulenceStatus.ok()&&runStatus.ok(),"public host fixture configuration restoration failed");
    }
    ~Cleanup(){Restore();}
  } cleanup{probe,originalRun,originalCoefficients,originalTurbulence,originalMover,originalActive,originalCoupling,
      originalLocation,originalMomentum,originalAnalytic,originalMacros,originalEnergy,originalMu};
  activeProbe=&probe;
  probe.Require(Run::InstallActive(controlled).ok(),"cannot install controlled host run configuration");
  probe.Require(Transport::Coefficient::SetActiveConfiguration(controlled.coefficients).ok(),"cannot install explicit controlled scattering");
  probe.Require(Turbulence::SetActiveConfiguration(controlled.turbulence).ok(),"cannot freeze prescribed wave controls");
  Mover::SelectProductionMover(selected);W::ActiveFlag=false;W::ParticleCouplingMode=false;
  FieldLine::InjectionParameters::InjectLocation=FieldLine::InjectionParameters::_InjectShockLocations;
  FieldLine::InjectionParameters::InjectionMomentumModel=FieldLine::InjectionParameters::_const_energy;
  FieldLine::InjectionParameters::UseAnalyticShockModel=FieldLine::InjectionParameters::AnalyticShockModel_Tenishev2005;
  FieldLine::InjectionParameters::nParticlesPerIteration=static_cast<int>(input.sourceMacros);
  FieldLine::InjectionParameters::ConstEnergyInjectionValue=input.sourceEnergyMeV;
  FieldLine::InjectionParameters::ConstMuInjectionValue=input.sourceMu;
  PIC::BC::UserDefinedParticleInjectionFunction=ObservedInjection;PIC::Mover::SetUserDefinedParticleMover(ObservedMover);
  // Explicit test-only metadata replacement describes the controlled settings;
  // it is not a provider fallback. The real host creates/freezes every window
  // and advances PIC's real clock. No SetInitialValue or manual clock update is
  // used anywhere in this process fixture.
  Background::SnapshotStore::Instance().ResetForTests();
  const double start=Background::SimulationTimeSeconds();
  const double dt=PIC::ParticleWeightTimeStep::GlobalTimeStep[0];
  double launchRadius=0;
  for(int line=0;line<FL::nFieldLine;++line) {
    const double* point=FL::FieldLinesAll[line].GetFirstVertex()->GetX();
    const double radius=std::hypot(std::hypot(point[0],point[1]),point[2]);
    if(line==0||radius<launchRadius)launchRadius=radius;
  }
  probe.launchRadiusM=launchRadius;
  // The existing first trajectory interval is constant 1800 km/s below0.1 AU.
  // Restrict this bounded reference to that declared interval; do not use the
  // trajectory routine being checked to generate its own endpoint reference.
  probe.Require(start==0&&launchRadius+1800000.0*input.steps*dt<0.1*1.495978707e11,
                "host smoke reference requires the existing pre-first-knot launch interval");
  if(!probe.failures)for(unsigned step=0;step<input.steps;++step) {
    probe.step=step;const double before=Background::SimulationTimeSeconds();
    amps_time_step();
    ++probe.completedSteps;
    const double after=Background::SimulationTimeSeconds();
    probe.Require(after==before+dt,"real PIC host clock did not advance exactly once");
    probe.Require(!A::NativeDriverWindowActive(),"real host left its birth window open");
    PhysicsCheckpoint::AnalyticShockCache cache;
    const auto captured=S::CaptureAnalyticShockCache(&cache);
    probe.Require(captured.ok()&&cache.initialized,"real host lost initialized analytic trajectory history");
    const double expected=launchRadius+1800000.0*after;
    probe.Require(std::fabs(cache.radiusM-expected)<=2e-12*expected&&cache.trajectory.epochS==after,
                  "completed host shock geometry/trajectory still represents departure instead of the actual end epoch");
    const auto history=PhysicsCheckpoint::CaptureHostStepCache();
    probe.Require(std::fabs(history.previousShockRadiusM-expected)<=2e-12*expected,
                  "post-step host shock ledger endpoint differs from actual completed epoch");
    PersistTrace();
  }
  probe.Require(probe.sourceCalls==input.steps,"real PIC did not dispatch the registered source once per host step");
  unsigned globalFirst=0;MPI_Allreduce(&probe.firstMoves,&globalFirst,1,MPI_UNSIGNED,MPI_SUM,MPI_GLOBAL_COMMUNICATOR);
  probe.Require(globalFirst>0,"real PIC performed no pending first-birth mover calls");
  PersistTrace(true);
  cleanup.Restore();
  result.status=probe.failures?Status::Fail:Status::Pass;result.message=probe.failures?probe.message:
      "isolated actual host/PIC steps observe step-start births, one real clock increment and completed analytic shock endpoints";
  result.metrics.push_back({"assertion_failures",static_cast<double>(probe.failures),0,"<=","count"});
  result.metrics.push_back({"global_first_birth_moves",static_cast<double>(globalFirst),0,">","count"});
  result.metrics.push_back({"completed_host_steps",static_cast<double>(probe.completedSteps),static_cast<double>(input.steps),"==","count"});
  result.configuration={"scope=actual-PIC-TimeStep-fixed-analytic-diagnostic;process-isolation-required=true",
      "production_resume_or_source_frame_or_thermal_or_material_or_polarized_closure_qualified=false",
      "waves=frozen-prescribed;source=existing-analytic-shock-constant-energy-controls;clock=actual-host-only"};
  result.artifacts.insert(result.artifacts.end(),probe.traceArtifacts.begin(),probe.traceArtifacts.end());return result;
}
Result Parker(){return RunHost(Mover::ProductionMover::Parker);}
Result Dmumu(){return RunHost(Mover::ProductionMover::FocusedTransportDiffusion);}
Result Mfp(){return RunHost(Mover::ProductionMover::FocusedTransportMeanFreePath);}
} // namespace
std::vector<Descriptor> HostStepNativeDescriptors() {
  const char* names[]={"Actual Parker host birth/shock clock","Actual Dmumu host birth/shock clock","Actual MFP host birth/shock clock"};
  const TestCallback callbacks[]={Parker,Dmumu,Mfp};std::vector<Descriptor> result;
  for(unsigned index=0;index<3;++index) {
    Descriptor d;d.id="HOST0"+std::to_string(index+1);d.name=names[index];d.description=names[index];
    d.group="host-step-native";d.tags={"native-process-exclusive"};
    d.initialization=InitializationLevel::FieldLineModel;d.runtime=RuntimeClass::Extended;
    d.seedCapability=SeedCapability::FixedContract;d.repeatSupported=false;
    d.seedPolicy="fixed native source campaign; explicit fresh process protocol; no seed override/repeat qualification";
    d.stateIsolation="one selected descriptor per fresh process; public controls restored; PIC private first-call/statistics state isolated by process termination";
    d.supportedBuildModes="fixed analytic zero-origin one-line one-species diagnostic; actual host clock policy reported, no SWMF epoch fallback";
    d.callback=callbacks[index];result.push_back(d);
  }
  return result;
}
} }
