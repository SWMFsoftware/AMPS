#include "sep_birth_native_validation.h"
#include "../transport_common.h"
#include "../coefficient_providers.h"
#include "../physics_checkpoint.h"
#include "sep_background_runtime.h"
#include "sep_source_lifecycle.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <iterator>
#include <limits>
#include <set>
#include <string>
#include <unistd.h>

namespace SEP { namespace Testing { namespace {
namespace T=Transport;
namespace A=T::PICAdapter;
namespace B=T::Birth;
namespace PB=PIC::ParticleBuffer;
namespace FL=PIC::FieldLine;
namespace W=AlfvenTurbulence_Kolmogorov;

struct Checks {
  unsigned failures=0;std::string message;
  double maximumCartesianViewResidualM=0.0;
  void Require(bool condition,const std::string& reason) {
    if(condition)return;
    ++failures;if(message.size()<1800){if(!message.empty())message+="; ";message+=reason;}
  }
  bool Status(const T::Status& status){Require(status.ok(),status.message);return status.ok();}
  bool Agree(bool local) {
    int value=local?1:0,all=0;MPI_Allreduce(&value,&all,1,MPI_INT,MPI_MIN,MPI_GLOBAL_COMMUNICATOR);
    if(!all)Require(false,"birth fixture prerequisite failed on an MPI rank");
    return all!=0;
  }
};

// All tests require an initially empty process. The existing corrected native
// codec saves/restores the audited diagnostic fields, including core
// RNG, clock, geometry buffers, provider metadata and private wave reservoirs.
// Only test-owned particles/windows are removed during cleanup. Closing a
// failed fixture after deletion is teardown, never evidence of physical motion.
class Fixture {
 public:
  Checks& check;
  T::Coefficient::Configuration originalCoefficient=T::Coefficient::ActiveConfiguration(),typed=originalCoefficient;
  Mover::ProductionMover originalMover=Mover::CurrentProductionMover();
  bool originalCoupling=W::ParticleCouplingMode;
  bool saved=false;
  std::string path;
  int segmentIndex=3;
  double initialParallel=1.5e5,initialNormal=2.0e5;
  std::uint64_t identity=0;
  SourceLifecycle::Event event;
  explicit Fixture(Checks& checks):check(checks) {
    path=PIC::ThisThread==0?"/tmp/srcsep-native-birth-"+std::to_string(static_cast<long long>(getpid())):std::string();
    int length=static_cast<int>(path.size());MPI_Bcast(&length,1,MPI_INT,0,MPI_GLOBAL_COMMUNICATOR);
    path.resize(length);MPI_Bcast(&path[0],length,MPI_CHAR,0,MPI_GLOBAL_COMMUNICATOR);
    typed.source=T::Coefficient::SourceMode::Prescribed;
    typed.pitchAngle=T::Coefficient::PitchAngleKind::Constant;typed.constantDmumuPerS=0.2;
    typed.spatial=T::Coefficient::SpatialKind::FromPitchAngle;
    typed.meanFreePath=T::Coefficient::MeanFreePathKind::Tenishev2005;
    if(!check.Status(T::Coefficient::SetActiveConfiguration(typed)))return;
    saved=check.Status(PhysicsCheckpoint::SaveCollective(path,0));
    event.campaign=TestSeed(5);event.index=UINT64_C(0x42495254484e);
    event.channel.kind=SourceLifecycle::Kind::FixedBeginning;
    identity=SourceLifecycle::ParticleIdentity(event,17);
    auto random=SourceLifecycle::Random(event,17,3);
    initialParallel+=1.0e5*random.UniformOpen01();
  }
  void DeleteTestParticles() {
    std::set<long int> particles;
    for(int line=0;line<FL::nFieldLine;++line)
      for(auto* segment=FL::FieldLinesAll[line].GetFirstSegment();segment;segment=segment->GetNext()) {
        for(const auto head:{segment->FirstParticleIndex,segment->tempFirstParticleIndex.load()})
          for(auto particle=head;particle!=-1;) {
            if(!particles.insert(particle).second)break;
            particle=PB::GetNext(particle);
          }
        segment->FirstParticleIndex=-1;segment->tempFirstParticleIndex=-1;
      }
    for(auto particle:particles)PB::DeleteParticle(particle);
  }
  ~Fixture() {
    DeleteTestParticles();
    if(A::NativeDriverWindowActive()) {
      B::DriverWindow window;
      if(check.Status(A::GetNativeDriverWindow(&window))) {
        PIC::SimulationTime::SetInitialValue(window.endEpochS);
        check.Status(A::CloseNativeDriverWindow());
      }
    }
    check.Require(A::PendingWaveContributionsEmpty(),"birth fixture left physical wave work for a second owner");
    Mover::SelectProductionMover(originalMover);W::ParticleCouplingMode=originalCoupling;
    check.Status(T::Coefficient::SetActiveConfiguration(typed));
    if(saved) {std::uint64_t completed=0;check.Status(PhysicsCheckpoint::LoadCollective(path,&completed));}
    check.Status(T::Coefficient::SetActiveConfiguration(originalCoefficient));
    if(PIC::ThisThread==0)std::remove(path.c_str());
    MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);
  }
  bool Prepare(Mover::ProductionMover selected,double epochS,double hostDtS) {
    if(!saved)return false;
    Mover::SelectProductionMover(selected);W::ParticleCouplingMode=false;
    PIC::SimulationTime::SetInitialValue(epochS);
    SetDt(hostDtS);
    // Finite validity makes use of the post-increment mutable PIC clock fail:
    // that clock is deliberately at the endpoint during the pending move.
    Background::PublishModelOwnedSnapshot(Background::Provider::Analytic,epochS,epochS+hostDtS,
        "controlled native birth-window timing fixture; unchanged analytic primitives");
    return true;
  }
  void SetDt(double dtS) {
    for(int species=0;species<PIC::nTotalSpecies;++species)
      PIC::ParticleWeightTimeStep::GlobalTimeStep[species]=dtS;
    PIC::ParticleWeightTimeStep::GlobalTimeStepInitialized=true;
  }
  bool Spawn(const B::Record& birth,std::uint64_t macro=0,int allocatingRank=-1,double weight=1.0) {
    auto* segment=FL::FieldLinesAll[0].GetSegment(segmentIndex);bool valid=true;
    const int allocationOwner=allocatingRank<0?segment->Thread:allocatingRank;
    if(allocationOwner==PIC::ThisThread) {
      const auto particle=PB::GetNewParticle(segment->FirstParticleIndex);
      PB::SetI(0,particle);PB::SetFieldLineId(0,particle);
      PB::SetFieldLineCoord(segmentIndex+0.45,particle);
      // Synthetic input is explicitly local-plasma-frame gyrotropic velocity;
      // this test supplies no Maxwellian/lab-frame distribution or conversion.
      PB::SetVParallel(initialParallel,particle);PB::SetVNormal(initialNormal,particle);
      PB::SetIndividualStatWeightCorrection(weight,particle);
      double x[3];FL::FieldLinesAll[0].GetCartesian(x,PB::GetFieldLineCoord(particle));PB::SetX(x,particle);
      valid=check.Status(A::InitializeParticleTransportStateWithIdentity(particle,
          SourceLifecycle::ParticleIdentity(event,17+macro),birth));
    }
    return check.Agree(valid);
  }
  std::string ParticleImage(bool includeBirthHistory=false)const {
    std::string image;
    for(auto* segment=FL::FieldLinesAll[0].GetFirstSegment();segment;segment=segment->GetNext())
      for(auto particle=segment->FirstParticleIndex;particle!=-1;particle=PB::GetNext(particle)) {
        std::string record(reinterpret_cast<const char*>(PB::GetParticleDataPointer(particle)),PB::ParticleDataLength);
        auto* raw=reinterpret_cast<PB::byte*>(&record[0]);PB::SetNext(-1,raw);PB::SetPrev(-1,raw);
        // Compare identical phase-space/RNG/hazard/identity observables while
        // intentionally distinguishing initial-condition vs source history.
        if(!includeBirthHistory)
          std::fill(record.begin()+Offset::BirthRecord,record.begin()+Offset::BirthRecord+sizeof(B::Record),0);
        image+=record;
      }
    return image;
  }
  bool Move(double readEpochS) {
    auto read=Background::SnapshotStore::Instance().BeginParticleRead(readEpochS);
    PIC::Mover::FieldLine::MoveParticles();
    return true;
  }
  void VerifyCartesianViews(bool requireNonzeroMotion) {
    int moved=0;
    for(auto* segment=FL::FieldLinesAll[0].GetFirstSegment();segment;segment=segment->GetNext())
      for(auto particle=segment->FirstParticleIndex;particle!=-1;particle=PB::GetNext(particle)) {
        const double coordinate=PB::GetFieldLineCoord(particle);
        double derived[3],stored[3];FL::FieldLinesAll[0].GetCartesian(derived,coordinate);PB::GetX(stored,particle);
        for(unsigned component=0;component<3;++component) {
          check.maximumCartesianViewResidualM=std::max(check.maximumCartesianViewResidualM,std::abs(derived[component]-stored[component]));
          check.Require(std::isfinite(stored[component])&&stored[component]==derived[component],
              "native mover committed field-line coordinate with stale derived Cartesian position");
        }
        if(coordinate!=segmentIndex+0.45)moved=1;
      }
    int any=0;MPI_Allreduce(&moved,&any,1,MPI_INT,MPI_MAX,MPI_GLOBAL_COMMUNICATOR);
    if(requireNonzeroMotion)check.Require(any!=0,"Cartesian-view regression requires actual nonzero coordinate motion");
  }
  void RejectNonfiniteCartesianCommit(double readEpochS,double hostDurationS) {
    std::vector<A::ParticleContext> contexts;
    {
      auto read=Background::SnapshotStore::Instance().BeginParticleRead(readEpochS);
      for(auto* segment=FL::FieldLinesAll[0].GetFirstSegment();segment;segment=segment->GetNext())
        for(auto particle=segment->FirstParticleIndex;particle!=-1;particle=PB::GetNext(particle)) {
          A::ParticleContext context;double effective=0.0;
          if(!check.Status(A::LoadParticle(particle,&context))||!check.Status(A::PrepareParticleStep(&context,hostDurationS,&effective)))continue;
          // A valid privately staged phase/RNG/birth candidate makes partial
          // publication observable. No stochastic history or pending flag may
          // escape if its coordinate's derived Cartesian view is invalid.
          context.state.coordinate+=0.01;
          T::KeyedRandomStream random(A::CampaignRandomSeed,context.stableParticleId,7,0);
          std::uint64_t intervals=0;
          if(!check.Status(A::ResumeRandomStream(context,A::StochasticOperator::Parker,&random,&intervals)))continue;
          random.UniformOpen01();
          if(!check.Status(A::StageRandomStream(&context,A::StochasticOperator::Parker,random,intervals+1)))continue;
          contexts.push_back(context);
        }
    }
    // Fault injection is outside every particle-read phase. Only this fixture
    // owns its vertices/particles, and all bytes are restored before another
    // read, source event or checkpoint. It is not a geometry update policy.
    Background::SnapshotStore::Instance().AssertProviderMayWrite(Background::Provider::Analytic);
    for(auto& context:contexts) {
      auto* segment=context.segment;const auto oldTemp=segment->tempFirstParticleIndex.load();
      check.Require(oldTemp==-1,"Cartesian rejection fixture requires an empty temporary attachment list");
      if(oldTemp!=-1)continue;
      double savedVertex[3];segment->GetBegin()->GetX(savedVertex);
      double badVertex[3]={std::numeric_limits<double>::quiet_NaN(),savedVertex[1],savedVertex[2]};
      const std::string before(reinterpret_cast<const char*>(context.data),PB::ParticleDataLength);
      segment->GetBegin()->SetX(badVertex);
      const auto rejected=A::CommitAndAttach(context);
      const bool unchanged=std::memcmp(context.data,before.data(),PB::ParticleDataLength)==0&&segment->tempFirstParticleIndex.load()==oldTemp;
      // A fail-before build still needs safe teardown/continued assertions:
      // restore the deliberately corrupted geometry and the old implementation's
      // observable partial publication, including its attachment head.
      segment->GetBegin()->SetX(savedVertex);
      std::memcpy(context.data,before.data(),PB::ParticleDataLength);segment->tempFirstParticleIndex=oldTemp;
      check.Require(!rejected.ok()&&rejected.code==T::StatusCode::InvalidParticleState&&unchanged,
          "nonfinite derived Cartesian position did not reject before phase/RNG/birth/list publication");
    }
  }
};

Result Finish(const Checks& check,const char* success) {
  unsigned failures=0;MPI_Allreduce(&check.failures,&failures,1,MPI_UNSIGNED,MPI_SUM,MPI_GLOBAL_COMMUNICATOR);
  Result result;result.status=failures?Status::Fail:Status::Pass;
  result.message=failures?(check.message.empty()?"native birth fixture failed on another MPI rank":check.message):success;
  result.hasSeed=true;result.seed=TestSeed(5);
  result.configuration={"evidence_level=native-PIC-three-mover-timing-and-persistence",
      "run_purpose=controlled-software-diagnostic;source_physics_qualified=false",
      "input_velocity_frame=explicit-local-plasma;Maxwellian-sampler=not-exercised",
      "mutable_clock=end-epoch-during-move;departure=explicit-frozen-window;finite-snapshot",
      "compiled_timestep_mode="+std::to_string(_SIMULATION_TIME_STEP_MODE_)+
          ";compiled_counter_mode="+std::to_string(_PIC_GLOBAL_TIME_COUNTER_MODE_)+
          ";compiled_coupler_mode="+std::to_string(_PIC_COUPLER_MODE_)};
  result.metrics.push_back({"assertion_failures",static_cast<double>(failures),0,"<=","count"});
  double largest=0.0;MPI_Allreduce(&check.maximumCartesianViewResidualM,&largest,1,MPI_DOUBLE,MPI_MAX,MPI_GLOBAL_COMMUNICATOR);
  result.metrics.push_back({"maximum_cartesian_view_residual",largest,0,"==","m"});
  return result;
}
bool Supported(Result* result) {
  const auto clock=A::NativeDriverClockCapability();
  if(!clock.ok()){result->status=Status::Skip;result->message=clock.message;return false;}
  if(Background::ConfiguredProvider()!=Background::Provider::Analytic||!W::ActiveFlag||
      !FL::FieldLinesAll||FL::nFieldLine<1||FL::FieldLinesAll[0].GetTotalSegmentNumber()<7) {
    result->status=Status::Skip;result->message="native birth fixture requires a fixed analytic line with at least seven segments and initialized wave storage";return false;
  }
  long int particles=PB::GetAllPartNum(),total=0;
  MPI_Allreduce(&particles,&total,1,MPI_LONG,MPI_SUM,MPI_GLOBAL_COMMUNICATOR);
  if(total||A::NativeDriverWindowActive()) {
    result->status=Status::Skip;result->message="run native birth fixture in a fresh empty test process outside a host window";return false;
  }
  return true;
}
Result RemainingDuration(Mover::ProductionMover selected) {
  Result missing;if(!Supported(&missing))return missing;
  Checks check;double remainingS=0,declaredEpochS=0;
  {
    // Startup already owns its published epoch/generation. Place this static
    // analytic diagnostic one second later rather than attempting a different
    // validity interval under that same immutable publication identity.
    Fixture fixture(check);const double start=Background::SimulationTimeSeconds()+1.0,hostDt=0.01;
    declaredEpochS=start;
    if(fixture.Prepare(selected,start,hostDt)&&check.Status(A::FreezeNativeDriverWindow())) {
      B::Record birth;
      if(check.Status(A::MakeUniformBirthRecord(0.75,&birth))&&fixture.Spawn(birth)) {
        remainingS=birth.firstMoveDurationS;
        // This is the actual counter-OFF ordering: the read phase starts at
        // the frozen departure, then PIC advances its mutable clock before
        // movers execute. The source record retains the missing time origin.
        {auto read=Background::SnapshotStore::Instance().BeginParticleRead(start);
         PIC::SimulationTime::SetInitialValue(start+hostDt);
         PIC::Mover::FieldLine::MoveParticles();}
        check.Status(A::CloseNativeDriverWindow());
        fixture.VerifyCartesianViews(true);
        const auto actual=fixture.ParticleImage();
        for(auto* segment=FL::FieldLinesAll[0].GetFirstSegment();segment;segment=segment->GetNext())
          for(auto particle=segment->FirstParticleIndex;particle!=-1;particle=PB::GetNext(particle)) {
            B::Record consumed;std::memcpy(&consumed,PB::GetParticleDataPointer(particle)+Offset::BirthRecord,sizeof(consumed));
            check.Require(consumed.pendingFirstMove==0&&consumed.firstMoveDurationS==remainingS,
                "native mover lost birth provenance or did not atomically consume first movement");
          }
        fixture.DeleteTestParticles();PIC::SimulationTime::SetInitialValue(start);
        B::Record initial;check.Status(B::InitialCondition(start,&initial));
        fixture.SetDt(remainingS);
        if(fixture.Spawn(initial)){fixture.Move(start);fixture.VerifyCartesianViews(true);}
        const auto expected=fixture.ParticleImage();
        check.Require(actual==expected,"pending native birth differs from identical particle advanced for exactly the independently prescribed remaining duration");
        fixture.DeleteTestParticles();fixture.SetDt(hostDt);
        if(fixture.Spawn(initial)){fixture.Move(start);fixture.VerifyCartesianViews(true);}
        const auto full=fixture.ParticleImage();
        int differs=actual!=full?1:0,any=0;
        MPI_Allreduce(&differs,&any,1,MPI_INT,MPI_MAX,MPI_GLOBAL_COMMUNICATOR);
        check.Require(any!=0,"negative control did not distinguish full-host movement from remaining-duration movement");
      }
    }
  }
  auto result=Finish(check,"actual native mover consumes the prescribed remaining birth duration once, preserving exact phase-space and RNG history against direct-duration reference");
  result.configuration.push_back("fixture_start_epoch_s="+std::to_string(declaredEpochS)+";static_start_delay_s=1;startup_snapshot_restored=true");
  result.metrics.push_back({"first_move_duration",remainingS,0.0025,"==","s"});return result;
}
Result Parker(){return RemainingDuration(Mover::ProductionMover::Parker);}
Result Dmumu(){return RemainingDuration(Mover::ProductionMover::FocusedTransportDiffusion);}
Result Mfp(){return RemainingDuration(Mover::ProductionMover::FocusedTransportMeanFreePath);}

Result AtomicTimingRejection() {
  Result missing;if(!Supported(&missing))return missing;
  Checks check;double lateDuration=0;
  {
    Fixture fixture(check);const double epoch=1.0e12,dt=0.01;
    if(fixture.Prepare(Mover::ProductionMover::FocusedTransportMeanFreePath,epoch,dt)&&
        check.Status(A::FreezeNativeDriverWindow())) {
      B::Record birth;check.Status(A::MakeUniformBirthRecord(std::nextafter(1.0,0.0),&birth));
      lateDuration=birth.firstMoveDurationS;
      check.Require(lateDuration>0&&epoch+dt-lateDuration==epoch+dt,
          "late-birth fixture did not expose absolute-epoch rounding loss");
      if(fixture.Spawn(birth)) {
        fixture.RejectNonfiniteCartesianCommit(epoch,dt);
        const auto before=fixture.ParticleImage(true);
        {auto read=Background::SnapshotStore::Instance().BeginParticleRead(epoch);
         PIC::SimulationTime::SetInitialValue(epoch+dt);
         const auto rejected=A::CloseNativeDriverWindow();
         check.Require(!rejected.ok()&&rejected.message.find("during particle motion")!=std::string::npos&&
             A::NativeDriverWindowActive()&&before==fixture.ParticleImage(true),
             "closing an active particle read phase changed or consumed its source window/history");
         for(auto* segment=FL::FieldLinesAll[0].GetFirstSegment();segment;segment=segment->GetNext())
           for(auto particle=segment->FirstParticleIndex;particle!=-1;particle=PB::GetNext(particle)) {
             auto* data=PB::GetParticleDataPointer(particle);B::Record stored;
             std::memcpy(&stored,data+Offset::BirthRecord,sizeof(stored));
             auto unknown=stored;unknown.schema=0;
             std::memcpy(data+Offset::BirthRecord,&unknown,sizeof(unknown));
             A::ParticleContext incompatible;
             const auto rejectedHistory=A::LoadParticle(particle,&incompatible);
             std::memcpy(data+Offset::BirthRecord,&stored,sizeof(stored));
             check.Require(!rejectedHistory.ok()&&before==fixture.ParticleImage(true),
                 "native load guessed missing birth history or mutated the restored record");
             A::ParticleContext context;check.Status(A::LoadParticle(particle,&context));
             double effective=123;
             const auto wrong=A::PrepareParticleStep(&context,0.5*dt,&effective);
             check.Require(!wrong.ok()&&effective==123&&before==fixture.ParticleImage(true),
                 "incompatible first-duration request changed source history or output");
           }
         PIC::Mover::FieldLine::MoveParticles();}
        check.Status(A::CloseNativeDriverWindow());
        for(auto* segment=FL::FieldLinesAll[0].GetFirstSegment();segment;segment=segment->GetNext())
          for(auto particle=segment->FirstParticleIndex;particle!=-1;particle=PB::GetNext(particle)) {
            B::Record consumed;T::OperatorRandomState random[3];
            const auto* data=PB::GetParticleDataPointer(particle);
            std::memcpy(&consumed,data+Offset::BirthRecord,sizeof(consumed));
            std::memcpy(random,data+Offset::OperatorRandomStates,sizeof(random));
            check.Require(!consumed.pendingFirstMove&&random[2].acceptedIntervals>0,
                "strictly positive late birth was discarded below the adaptive floor");
          }
      }
    }
  }
  auto result=Finish(check,"read-phase close and incompatible requests reject without mutation; a representable late birth still executes an actual MFP interval below the adaptive refinement floor");
  result.configuration.push_back("fixture_start_epoch_s=1000000000000;startup_snapshot_restored=true");
  result.metrics.push_back({"preserved_late_first_move_duration",lateDuration,0,">","s"});return result;
}

Result CloneSplitMigration() {
  Result missing;if(!Supported(&missing))return missing;
  Checks check;std::uint64_t splitPopulation=0;double declaredEpochS=0;
  {
    Fixture fixture(check);const double epoch=Background::SimulationTimeSeconds()+1.0,dt=0.01;
    declaredEpochS=epoch;
    if(fixture.Prepare(Mover::ProductionMover::FocusedTransportMeanFreePath,epoch,dt)&&
        check.Status(A::FreezeNativeDriverWindow())) {
      B::Record birth;check.Status(A::MakeUniformBirthRecord(0.75,&birth));
      auto* segment=FL::FieldLinesAll[0].GetSegment(fixture.segmentIndex);
      const int owner=segment->Thread;
      const int originalRank=PIC::nTotalThreads>1?(owner+1)%PIC::nTotalThreads:owner;
      bool created=true;
      for(unsigned macro=0;macro<5;++macro)
        created=fixture.Spawn(birth,macro,originalRank,macro==0?8.0:1.0)&&created;
      if(created) {
        if(PIC::ThisThread==originalRank) {
          // Generic cloning is a transfer/copy primitive. It must preserve
          // physical birth age and pending time; only the actual statistical
          // split callback may assign a child lineage/future random stream.
          const auto parent=segment->FirstParticleIndex,copy=PB::GetNewParticle();
          PB::CloneParticle(copy,parent);
          const auto* a=PB::GetParticleDataPointer(parent);const auto* b=PB::GetParticleDataPointer(copy);
          check.Require(std::memcmp(a+Offset::BirthRecord,b+Offset::BirthRecord,sizeof(B::Record))==0&&
              std::memcmp(a+Offset::OperatorRandomStates,b+Offset::OperatorRandomStates,3*sizeof(T::OperatorRandomState))==0,
              "generic CloneParticle changed pending physical birth or stochastic history");
          PB::DeleteParticle(copy);
        }
        // Place source records on another rank deliberately, then execute the
        // actual native MPI migration before their first movement. The full
        // particle payload, including all64 birth bytes, must reach the owner.
        PIC::ParallelFieldLines::ExchangeFieldLineParticles();
        int localCount=0;bool migrationValid=true;
        for(auto particle=segment->FirstParticleIndex;particle!=-1;particle=PB::GetNext(particle)) {
          ++localCount;B::Record migrated;
          std::memcpy(&migrated,PB::GetParticleDataPointer(particle)+Offset::BirthRecord,sizeof(migrated));
          migrationValid=migrationValid&&PIC::ThisThread==owner&&std::memcmp(&birth,&migrated,sizeof(birth))==0;
        }
        check.Require(localCount==(PIC::ThisThread==owner?5:0),"native exchange changed source population/owner");
        if(check.Agree(migrationValid)) {
          {auto read=Background::SnapshotStore::Instance().BeginParticleRead(epoch);
           PIC::SimulationTime::SetInitialValue(epoch+dt);PIC::Mover::FieldLine::MoveParticles();}
          check.Status(A::CloseNativeDriverWindow());
          A::RegisterPopulationControlHooks();A::BeginPopulationControl();
          if(PIC::ThisThread==owner) {
            double parallelRange=0,normalRange=0;
            // Existing exact half-weight core controller: five representatives
            // become six. This is a statistical copy after the closed physical
            // step, not a new physical source event or a reset birth clock.
            segment->WeightedParticleSplitting(0,segment->FirstParticleIndex,2,2,2,1.0,
                parallelRange,normalRange,5,7,2);
          }
          auto expected=birth;expected.pendingFirstMove=0;
          for(auto particle=segment->FirstParticleIndex;particle!=-1;particle=PB::GetNext(particle)) {
            ++splitPopulation;B::Record inherited;
            std::memcpy(&inherited,PB::GetParticleDataPointer(particle)+Offset::BirthRecord,sizeof(inherited));
            check.Require(std::memcmp(&inherited,&expected,sizeof(inherited))==0,
                "actual exact split invented a new birth age or lost consumed source-window provenance");
          }
        }
      }
    }
  }
  std::uint64_t total=0;MPI_Allreduce(&splitPopulation,&total,1,MPI_UINT64_T,MPI_SUM,MPI_GLOBAL_COMMUNICATOR);
  check.Require(total==6,"actual exact split did not retain six physical-history representatives");
  auto result=Finish(check,"generic clone, actual native MPI exchange and actual exact half-weight split preserve physical birth provenance and pending/consumed history");
  result.configuration.push_back("fixture_start_epoch_s="+std::to_string(declaredEpochS)+";static_start_delay_s=1;startup_snapshot_restored=true");
  result.configuration.push_back(std::string("actual_interrank_transfer=")+(PIC::nTotalThreads>1?"true":"false;requires-MPI2-for-transfer-gate"));
  result.metrics.push_back({"post_split_representatives",static_cast<double>(total),6,"==","count"});return result;
}

std::string ReadRootFile(const std::string& path) {
  if(PIC::ThisThread!=0)return std::string();
  std::ifstream input(path.c_str(),std::ios::binary);
  return std::string((std::istreambuf_iterator<char>(input)),std::istreambuf_iterator<char>());
}
Result CompletedCheckpointRejection() {
  Result missing;if(!Supported(&missing))return missing;
  Checks check;double declaredEpochS=0;
  {
    Fixture fixture(check);const auto originalFile=ReadRootFile(fixture.path);
    const double epoch=Background::SimulationTimeSeconds()+1.0,dt=0.01;
    declaredEpochS=epoch;
    if(fixture.Prepare(Mover::ProductionMover::FocusedTransportMeanFreePath,epoch,dt)&&
        check.Status(A::FreezeNativeDriverWindow())) {
      B::Record pending;check.Status(A::MakeUniformBirthRecord(0.5,&pending));
      if(fixture.Spawn(pending)) {
        const auto before=fixture.ParticleImage(true);std::uint64_t completed=123;
        const auto activeSave=PhysicsCheckpoint::SaveCollective(fixture.path,0);
        const auto activeLoad=PhysicsCheckpoint::LoadCollective(fixture.path,&completed);
        check.Require(!activeSave.ok()&&!activeLoad.ok()&&completed==123&&
            before==fixture.ParticleImage(true)&&A::NativeDriverWindowActive(),
            "active-window checkpoint saved/loaded partial source movement or mutated physics");
        check.Require(ReadRootFile(fixture.path)==originalFile,"rejected active-window save changed the existing checkpoint destination");
        {auto read=Background::SnapshotStore::Instance().BeginParticleRead(epoch);
         PIC::SimulationTime::SetInitialValue(epoch+dt);PIC::Mover::FieldLine::MoveParticles();}
        check.Status(A::CloseNativeDriverWindow());
        // Fuzz only the stored history after a successfully completed step.
        // This models an incompatible incoming/current record, never a valid
        // source API. The completed-step codec must reject pending state even
        // when there is no active host-window flag to catch it earlier.
        for(auto* segment=FL::FieldLinesAll[0].GetFirstSegment();segment;segment=segment->GetNext())
          for(auto particle=segment->FirstParticleIndex;particle!=-1;particle=PB::GetNext(particle))
            std::memcpy(PB::GetParticleDataPointer(particle)+Offset::BirthRecord,&pending,sizeof(pending));
        const auto pendingImage=fixture.ParticleImage(true);
        const auto pendingSave=PhysicsCheckpoint::SaveCollective(fixture.path,1);
        check.Require(!pendingSave.ok()&&pendingSave.message.find("unconsumed scheduled first move")!=std::string::npos&&
            pendingImage==fixture.ParticleImage(true),"completed checkpoint accepted or consumed pending birth history");
        check.Require(ReadRootFile(fixture.path)==originalFile,"rejected pending-record save changed the checkpoint destination");
        auto consumed=pending;consumed.pendingFirstMove=0;
        for(auto* segment=FL::FieldLinesAll[0].GetFirstSegment();segment;segment=segment->GetNext())
          for(auto particle=segment->FirstParticleIndex;particle!=-1;particle=PB::GetNext(particle))
            std::memcpy(PB::GetParticleDataPointer(particle)+Offset::BirthRecord,&consumed,sizeof(consumed));
      }
    }
  }
  auto result=Finish(check,"completed checkpoint API rejects active windows and unconsumed first-move records without changing particles, caller state or an existing destination");
  result.configuration.push_back("fixture_start_epoch_s="+std::to_string(declaredEpochS)+";static_start_delay_s=1;startup_snapshot_restored=true");
  return result;
}
Descriptor Make(const char* id,const char* name,TestCallback callback) {
  Descriptor descriptor;descriptor.id=id;descriptor.name=name;descriptor.description=name;
  descriptor.group="birth-window-native";descriptor.initialization=InitializationLevel::FieldLineModel;
  descriptor.runtime=RuntimeClass::Routine;descriptor.seedCapability=SeedCapability::Context;
  descriptor.repeatSupported=true;descriptor.seedPolicy="requested seed drives actual synthetic source identity and local-plasma velocity; full original physics restored";
  descriptor.supportedBuildModes="fixed analytic native PIC model with equal global dt and corrected checkpoint capability";
  descriptor.stateIsolation="initially empty process; test records/windows removed; original full native physics checkpoint restored";
  descriptor.callback=callback;return descriptor;
}
} // namespace
std::vector<Descriptor> BirthNativeDescriptors() {
  return {Make("BIRTHN01","Native Parker consumes source remaining duration",Parker),
      Make("BIRTHN02","Native Dmumu consumes source remaining duration",Dmumu),
      Make("BIRTHN03","Native MFP consumes source remaining duration",Mfp),
      Make("BIRTHN04","Native birth timing rejection and late interval",AtomicTimingRejection),
      Make("BIRTHN05","Native clone split and migration retain birth history",CloneSplitMigration),
      Make("BIRTHN06","Completed checkpoint refuses pending birth movement",CompletedCheckpointRejection)};
}
} }
