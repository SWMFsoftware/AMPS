#include "sep_checkpoint_native_validation.h"
#include "../physics_checkpoint.h"
#include "../physics_checkpoint_state.h"
#include "../flux_tube_checkpoint_state.h"
#include "../transport_common.h"
#include "../turbulence_production_adapter.h"
#include "../turbulence_wave_number_resolved/turbulence_wave_number_resolved.h"
#include "../turbulence_advection_kolmogorov.h"
#include "sep_background_runtime.h"
#include "sep_initialization.h"
#include "sep_run_configuration.h"
#include "sep_source_lifecycle.h"
#include "sep_turbulence_core.h"
#include "pic_persistent_format.h"

#include <algorithm>
#include <cmath>
#include <climits>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iterator>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unistd.h>
#include <vector>

namespace SEP { namespace Testing { namespace {
namespace T=Transport;
namespace A=T::PICAdapter;
namespace CP=PhysicsCheckpoint;
namespace PB=PIC::ParticleBuffer;
namespace FL=PIC::FieldLine;
namespace W=SEP::AlfvenTurbulence_Kolmogorov;

struct Checks {
  unsigned failures=0;std::string message;
  void Require(bool ok,const std::string& text) {
    if(ok)return;
    ++failures;
    if(message.size()<1500){if(!message.empty())message+="; ";message+=text;}
  }
  bool Status(const T::Status& status) {Require(status.ok(),status.message);return status.ok();}
};
std::string SharedPrefix() {
  std::string path=PIC::ThisThread==0?"/tmp/srcsep-native-checkpoint-"+std::to_string(static_cast<long long>(::getpid())):std::string();
  int size=static_cast<int>(path.size());MPI_Bcast(&size,1,MPI_INT,0,MPI_GLOBAL_COMMUNICATOR);
  path.resize(size);MPI_Bcast(&path[0],size,MPI_CHAR,0,MPI_GLOBAL_COMMUNICATOR);return path;
}
class Fixture {
 public:
  Checks& checks;
  std::string prefix;
  T::Coefficient::Configuration coefficient=T::Coefficient::ActiveConfiguration();
  Turbulence::Configuration turbulence=Turbulence::ActiveConfiguration();
  Mover::ProductionMover mover=Mover::CurrentProductionMover();
  bool coupling=W::ParticleCouplingMode;
  bool saved=false;
  Fixture(Checks& target):checks(target),prefix(SharedPrefix()) {
    auto typed=coefficient;typed.pitchAngle=T::Coefficient::PitchAngleKind::Constant;
    typed.constantDmumuPerS=0.2;typed.spatial=T::Coefficient::SpatialKind::FromPitchAngle;
    typed.meanFreePath=T::Coefficient::MeanFreePathKind::Tenishev2005;
    typed.source=T::Coefficient::SourceMode::Prescribed;
    if(!checks.Status(T::Coefficient::SetActiveConfiguration(typed)))return;
    saved=checks.Status(CP::SaveCollective(prefix+".original",0));
  }
  ~Fixture() {
    // The production turbulence transaction is the sole queue-flush owner.
    // Cleanup must detect an unfinished step rather than applying new physical
    // exchange while restoring the fixture's original checkpoint.
    checks.Require(A::PendingWaveContributionsEmpty(),
        "checkpoint fixture left uncommitted particle-wave contributions");
    if(saved) {
      Mover::SelectProductionMover(mover);W::ParticleCouplingMode=coupling;
      const auto config=Turbulence::SetActiveConfiguration(turbulence);
      checks.Require(config.ok(),"cannot restore fixture turbulence configuration");
      std::uint64_t steps=0;const auto restored=CP::LoadCollective(prefix+".original",&steps);
      checks.Require(restored.ok(),"cannot restore original native physics image: "+restored.message);
    }
    checks.Require(T::Coefficient::SetActiveConfiguration(coefficient).ok(),"cannot restore coefficient configuration");
    if(PIC::ThisThread==0)for(const char* suffix:{".original",".checkpoint",".bad",".sentinel"})
      std::remove((prefix+suffix).c_str());
    MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);
  }
  // Configure only: fresh-process readers must not initialize wave/flux
  // caches or create source events before loading the recorded realization.
  bool Configure(Mover::ProductionMover selected) {
    if(!saved)return false;
    Mover::SelectProductionMover(selected);W::ParticleCouplingMode=false;
    auto waves=turbulence;waves.source=W::WaveNumberResolved::IsActive()?
        Turbulence::Source::SelfConsistentSpectral:Turbulence::Source::SelfConsistentIntegrated;
    waves.representation=W::WaveNumberResolved::IsActive()?Turbulence::Representation::Spectral:Turbulence::Representation::Integrated;
    waves.advectionEnabled=true;waves.reflectionEnabled=false;waves.cascadeEnabled=false;
    waves.shockInjectionEnabled=false;waves.coupling=Turbulence::CouplingPolicy::Disabled;
    if(!checks.Status(Turbulence::SetActiveConfiguration(waves)))return false;
    for(int species=0;species<PIC::nTotalSpecies;++species)
      PIC::ParticleWeightTimeStep::GlobalTimeStep[species]=0.01;
    PIC::ParticleWeightTimeStep::GlobalTimeStepInitialized=true;
    return true;
  }
  bool Prepare(Mover::ProductionMover selected) {
    if(!saved)return false;
    // Check actual startup storage before this fixture writes any native
    // primitive. In particular this catches the old driver-only density
    // initialization, which native component tests bypassed entirely.
    int localPlasmaReady=1;
    double origin[3]={0.0,0.0,0.0};
    if(Initialization::HasActive()) {
      const auto& center=Initialization::Active().parkerOriginM;
      origin[0]=center.x;origin[1]=center.y;origin[2]=center.z;
    }
    for(int l=0;l<FL::nFieldLine;++l)
      for(auto* vertex=FL::FieldLinesAll[l].GetFirstVertex();vertex;vertex=vertex->GetNext()) {
        const double* x=vertex->GetX();
        const double radiusM=std::sqrt((x[0]-origin[0])*(x[0]-origin[0])+
            (x[1]-origin[1])*(x[1]-origin[1])+(x[2]-origin[2])*(x[2]-origin[2]));
        const double expected=5.0e6*std::pow(1.495978707e11/radiusM,2);
        double current=0.0,previous=0.0;
        vertex->GetDatum(FL::DatumAtVertexPlasmaDensity,&current);
        vertex->GetDatum(FL::DatumAtVertexPrevious::DatumAtVertexPlasmaDensity,&previous);
        if(!std::isfinite(expected)||expected<=0.0||!std::isfinite(current)||current<=0.0||
            !std::isfinite(previous)||previous<=0.0||
            std::fabs(current-expected)>16*std::numeric_limits<double>::epsilon()*expected||
            current!=previous)localPlasmaReady=0;
      }
    int plasmaReady=0;
    MPI_Allreduce(&localPlasmaReady,&plasmaReady,1,MPI_INT,MPI_MIN,MPI_GLOBAL_COMMUNICATOR);
    checks.Require(plasmaReady!=0,"shared analytic startup did not initialize positive current/previous density with the declared inverse-square law");
    if(!plasmaReady)return false;
    if(!Configure(selected))return false;
    int localReady=1;
    for(int l=0;l<FL::nFieldLine;++l)for(int s=0;s<FL::FieldLinesAll[l].GetTotalSegmentNumber();++s) {
      auto& line=FL::FieldLinesAll[l];auto* segment=line.GetSegment(s);double magnetic[3];
      line.GetMagneticField(magnetic,s+0.5);
      const double b2=magnetic[0]*magnetic[0]+magnetic[1]*magnetic[1]+magnetic[2]*magnetic[2];
      const double volume=SEP::FieldLine::FluxTubeGeometry::SegmentVolumeM3(segment,l);
      auto* energy=segment->GetDatum_ptr(W::CellIntegratedWaveEnergy);
      auto* derived=segment->GetDatum_ptr(W::WaveEnergyDensity);
      if(!energy||!derived||!std::isfinite(b2)||b2<=0.0||!std::isfinite(volume)||volume<=0.0) {
        checks.Require(false,"native checkpoint wave fixture requires finite positive B and magnetic tube volume");localReady=0;continue;
      }
      // A declared nonzero magnetic initial condition, deltaB/B=0.1 with
      // E+:E-=3:2. Energies are integrated J; no field or plasma replacement.
      const double density=0.01*b2/(2.0*VacuumPermeability);
      energy[0]=0.6*density*volume;energy[1]=0.4*density*volume;
      derived[0]=0.6*density;derived[1]=0.4*density;derived[2]=0.2;
    }
    int ready=0;MPI_Allreduce(&localReady,&ready,1,MPI_INT,MPI_MIN,MPI_GLOBAL_COMMUNICATOR);
    if(!ready){checks.Require(false,"native wave fixture initialization failed on an MPI rank");return false;}
    W::ResetRightBoundaryEminusInitialCondition();W::CaptureRightBoundaryEminusInitialCondition(W::CellIntegratedWaveEnergy);
    if(W::WaveNumberResolved::IsActive()) {
      W::WaveNumberResolved::InitializeSpectrumFromIntegratedEnergy(W::CellIntegratedWaveEnergy);
      W::WaveNumberResolved::ResetRightBoundarySpectrumInitialCondition();
      W::WaveNumberResolved::CaptureRightBoundarySpectrumInitialCondition();
    }
    return Births();
  }
  bool Births() {
    // Every rank advances the same physical channel. Segment index is the
    // source macro label, so changing which rank owns it cannot duplicate ids.
    const auto event=SourceLifecycle::BeginNativeEvent(SourceLifecycle::Kind::FixedBeginning,0,0);
    auto& line=FL::FieldLinesAll[0];const int count=line.GetTotalSegmentNumber();
    int selected=-1;
    for(int segment=2;segment+2<count;++segment)
      if(line.GetSegment(segment)->Thread==PIC::ThisThread){selected=segment;break;}
    if(selected<0)return true;
    auto* segment=line.GetSegment(selected);
    for(unsigned i=0;i<2;++i) {
      const auto macro=(static_cast<std::uint64_t>(selected)<<32)|i;
      auto random=SourceLifecycle::Random(event,macro,2);
      const auto particle=PB::GetNewParticle(segment->FirstParticleIndex);
      PB::SetI(0,particle);PB::SetFieldLineId(0,particle);
      PB::SetFieldLineCoord(selected+0.35+0.1*i,particle);
      PB::SetVParallel(1.0e5*(1.0+random.UniformOpen01()),particle);
      PB::SetVNormal(2.0e5,particle);PB::SetIndividualStatWeightCorrection(1.0,particle);
      double x[3];line.GetCartesian(x,PB::GetFieldLineCoord(particle));PB::SetX(x,particle);
      A::InitializeParticleTransportStateWithIdentity(particle,SourceLifecycle::ParticleIdentity(event,macro));
      if(Mover::CurrentProductionMover()==Mover::ProductionMover::FocusedTransportMeanFreePath) {
        // A positive conditional residual, not a replacement exponential law.
        // This deterministic fixture starts just before its next physical event
        // so native restart exercises an actual scattering history immediately.
        const double residual=1.0e-12;
        std::memcpy(PB::GetParticleDataPointer(particle)+SEP::Offset::MfpOpticalDepth,&residual,sizeof(residual));
      }
    }
    return true;
  }
  bool Step() {
    if(!Births())return false;
    Background::PrepareSnapshotForParticleStep();
    {
      auto read=Background::SnapshotStore::Instance().BeginParticleRead(Background::SimulationTimeSeconds());
      PIC::Mover::FieldLine::MoveParticles();
    }
    PIC::SimulationTime::SetInitialValue(Background::SimulationTimeSeconds()+0.01);
    const auto status=Turbulence::PICAdapter::Advance(0.01,
        std::numeric_limits<double>::quiet_NaN(),std::numeric_limits<double>::quiet_NaN());
    const int local=checks.Status(status)?1:0;int global=0;
    MPI_Allreduce(&local,&global,1,MPI_INT,MPI_MIN,MPI_GLOBAL_COMMUNICATOR);
    checks.Require(global!=0,"native checkpoint step failed on another MPI rank");
    return global!=0;
  }
};

Result Finish(const Checks& checks,const std::string& success) {
  Result result;int failures=checks.failures,total=0;
  MPI_Allreduce(&failures,&total,1,MPI_INT,MPI_SUM,MPI_GLOBAL_COMMUNICATOR);
  result.status=total?Status::Fail:Status::Pass;result.message=total?checks.message:success;
  if(total&&result.message.empty())result.message="checkpoint assertion failed on another MPI rank";
  result.configuration={"evidence_level=native-PIC-operators;fixed-geometry;analytic-provider",
      "run_purpose=controlled-software-diagnostic;production_physics_qualified=false",
      "MPI_policy=same-rank-count-and-segment-ownership;output_continuation=false",
      "checkpoint_scope=same-process-component-restore;fresh-process-resume=not-qualified;production-driver-resume=not-wired",
      "dt_s=0.01;typed_Dmumu_per_s=0.2;wave_operator=production-advection"};
  result.hasSeed=true;result.seed=A::CampaignRandomSeed;
  result.metrics.push_back({"assertion_failures",static_cast<double>(total),0,"<=","count"});
  return result;
}
bool Supported(Result* result) {
  if(Background::ConfiguredProvider()!=Background::Provider::Analytic||!W::ActiveFlag||
      !FL::FieldLinesAll||FL::nFieldLine<1||FL::FieldLinesAll[0].GetTotalSegmentNumber()<6) {
    result->status=Status::Skip;result->message="native checkpoint fixture requires standalone analytic provider, active waves, and a fixed line of at least six segments";
    return false;
  }
  const long int local=PB::GetAllPartNum();long int global=0;
  MPI_Allreduce(&local,&global,1,MPI_LONG,MPI_SUM,MPI_GLOBAL_COMMUNICATOR);
  if(global!=0){result->status=Status::Skip;result->message="run native checkpoint fixture in a fresh test subprocess with an empty particle population";return false;}
  return true;
}

bool WriteInvalidRawFrame(const std::string& valid,const std::string& bad,
    const std::string& record,std::size_t offset,const std::vector<double>& replacement,
    Checks& checks) {
  int good=1;
  if(PIC::ThisThread==0)try {
    if(offset>record.size()||replacement.size()>(record.size()-offset)/sizeof(double))
      throw std::runtime_error("malformed-frame fixture replacement exceeds pinned native record");
    std::ifstream input(valid.c_str(),std::ios::binary);
    const std::vector<unsigned char> file((std::istreambuf_iterator<char>(input)),std::istreambuf_iterator<char>());
    if(file.size()<=PIC::PersistentFormat::HeaderSize)throw std::runtime_error("cannot read native malformed-frame fixture");
    const auto layout=PIC::PersistentFormat::Get64(file.data()+32);
    std::string payload(file.begin()+PIC::PersistentFormat::HeaderSize,file.end());
    // Match the codec's complete length-prefixed raw string, not an arbitrary
    // floating-point byte pattern elsewhere in metadata or another field.
    std::string framed;
    for(unsigned i=0;i<8;++i)framed.push_back(static_cast<char>(static_cast<std::uint64_t>(record.size())>>(8*i)));
    framed+=record;const auto found=payload.find(framed);
    if(found==std::string::npos)throw std::runtime_error("cannot locate complete length-framed native record in checkpoint");
    std::memcpy(&payload[found+8+offset],replacement.data(),replacement.size()*sizeof(double));
    const std::vector<unsigned char> bytes(payload.begin(),payload.end());
    const auto written=PIC::PersistentFormat::WriteAtomic(bad,PIC::PersistentFormat::Kind::ParticleRestart,
        PIC::PersistentFormat::NativeLayout,layout,bytes);
    if(written!=PIC::PersistentFormat::Status::Success)throw std::runtime_error("cannot publish CRC-valid malformed native frame");
  }catch(const std::exception& e){checks.Require(false,e.what());good=0;}
  MPI_Bcast(&good,1,MPI_INT,0,MPI_GLOBAL_COMMUNICATOR);
  return good!=0;
}
std::string RepresentativeParticleRecord(Checks& checks) {
  std::string record;
  for(int line=0;line<FL::nFieldLine&&record.empty();++line)
    for(auto* segment=FL::FieldLinesAll[line].GetFirstSegment();segment&&record.empty();segment=segment->GetNext())
      if(segment->FirstParticleIndex!=-1) {
        record.assign(reinterpret_cast<const char*>(PB::GetParticleDataPointer(segment->FirstParticleIndex)),PB::ParticleDataLength);
        auto* raw=reinterpret_cast<PB::byte*>(&record[0]);PB::SetNext(-1,raw);PB::SetPrev(-1,raw);
      }
  const int candidate=record.empty()?INT_MAX:PIC::ThisThread;int owner=INT_MAX;
  MPI_Allreduce(&candidate,&owner,1,MPI_INT,MPI_MIN,MPI_GLOBAL_COMMUNICATOR);
  if(owner==INT_MAX){checks.Require(false,"malformed-record fixture has no native particle");return std::string();}
  int length=static_cast<int>(record.size());MPI_Bcast(&length,1,MPI_INT,owner,MPI_GLOBAL_COMMUNICATOR);
  record.resize(length);MPI_Bcast(&record[0],length,MPI_CHAR,owner,MPI_GLOBAL_COMMUNICATOR);
  return record;
}
bool WriteInvalidCoordinate(const std::string& valid,const std::string& bad,Checks& checks) {
  const auto record=RepresentativeParticleRecord(checks);if(record.empty())return false;
  return WriteInvalidRawFrame(valid,bad,record,_PIC_PARTICLE_DATA__FIELD_LINE_COORD_OFFSET_,
      {std::numeric_limits<double>::infinity()},checks);
}
bool WritePendingBirth(const std::string& valid,const std::string& bad,Checks& checks) {
  const auto record=RepresentativeParticleRecord(checks);if(record.empty())return false;
  T::Birth::DriverWindow window;T::Birth::Record birth;
  if(!checks.Status(T::Birth::FreezeWindow(Background::SimulationTimeSeconds(),0.01,&window))||
      !checks.Status(T::Birth::UniformBirth(window,0.5,&birth)))return false;
  // Mutate the complete length-framed value record and rebuild the envelope
  // CRC. These double storage slots are copied as bytes, never interpreted as
  // numeric source data; this tests semantic completed-step rejection.
  std::vector<double> replacement(sizeof(birth)/sizeof(double));
  std::memcpy(replacement.data(),&birth,sizeof(birth));
  return WriteInvalidRawFrame(valid,bad,record,Offset::BirthRecord,replacement,checks);
}
bool WriteInvalidParticleValue(const std::string& valid,const std::string& bad,
    std::size_t offset,const std::vector<double>& replacement,Checks& checks) {
  const auto record=RepresentativeParticleRecord(checks);if(record.empty())return false;
  return WriteInvalidRawFrame(valid,bad,record,offset,replacement,checks);
}
double OpaqueFixtureArea(const double*,int){return 1.0;}
Result FluxOwnerContinuation() {
  Result unsupported;if(!Supported(&unsupported))return unsupported;
  Checks checks;
  {
    Fixture fixture(checks);
    if(fixture.saved) {
      namespace F=FieldLine::FluxTubeGeometry;
      FieldLine::FluxTubeState::State original,saved,restored;
      if(checks.Status(F::CaptureCheckpointState(&original))) {
        const double value=F::GetMagneticFluxWb(0)*1.25;
        F::SetMagneticFluxWb(0,value,"native flux persistence diagnostic",9);
        // Initialize J and boundary reservoirs using the declared nondefault
        // volume, so their density normalization is physically consistent.
        if(fixture.Prepare(Mover::ProductionMover::Parker)) {
        checks.Status(F::CaptureCheckpointState(&saved));
        if(checks.Status(CP::SaveCollective(fixture.prefix+".checkpoint",0))) {
          F::SetMagneticFluxWb(0,value*1.5,"changed after checkpoint",10);
          std::uint64_t completed=99;
          if(checks.Status(CP::LoadCollective(fixture.prefix+".checkpoint",&completed))) {
            checks.Status(F::CaptureCheckpointState(&restored));
            std::string expected,actual;
            checks.Status(FieldLine::FluxTubeState::Serialize(saved,&expected));
            checks.Status(FieldLine::FluxTubeState::Serialize(restored,&actual));
            checks.Require(expected==actual,"checkpoint lost actual flux value/generation/provenance or lazy slot state");
          }
        }
        F::CheckpointRestorePlan plan;
        if(checks.Status(F::StageCheckpointRestore(saved,&plan))&&checks.Status(F::ValidateCheckpointReady(plan))) {
          F::SetMagneticFluxWb(0,value*1.75,"stale plan negative control",11);
          checks.Require(!F::CommitCheckpointRestore(&plan),"flux commit accepted a stale owner after readiness");
          checks.Require(!F::ValidateCheckpointReady(plan).ok(),"flux readiness accepted changed owner serial");
        }
        F::SetExplicitAreaProfile(OpaqueFixtureArea);
        const auto opaque=CP::SaveCollective(fixture.prefix+".checkpoint",0);
        checks.Require(opaque.code==T::StatusCode::UnsupportedConfiguration&&F::HasOpaqueAreaProfile(),
            "opaque area callback was serialized or silently cleared");
        F::ClearExplicitAreaProfile();
        if(checks.Status(F::StageCheckpointRestore(original,&plan))&&checks.Status(F::ValidateCheckpointReady(plan)))
          checks.Require(F::CommitCheckpointRestore(&plan),"valid quiescent flux restore failed to commit");
        }
      }
    }
  }
  return Finish(checks,"actual flux registry roundtrips through native checkpoint; stale-owner commit and opaque callback reject without guessed normalization");
}

bool FreshRequest(std::string* role,std::string* prefix,Checks& checks) {
  int valid=1;
  if(PIC::ThisThread==0) {
    std::ifstream input(GetExecutionContext().inputPath.c_str());
    std::string magic,roleKey,prefixKey;unsigned version=0;
    if(!(input>>magic>>version>>roleKey>>std::quoted(*role)>>prefixKey>>std::quoted(*prefix))||
        magic!="SEP_CHECKPOINT_FRESH_PROCESS"||version!=1||roleKey!="role"||prefixKey!="prefix"||
        (*role!="write"&&*role!="resume")||prefix->empty()||(*prefix)[0]!='/')valid=0;
    input>>std::ws;if(!input.eof())valid=0;
  }
  MPI_Bcast(&valid,1,MPI_INT,0,MPI_GLOBAL_COMMUNICATOR);
  checks.Require(valid!=0,"fresh-process fixture input must declare protocol1, write/resume role and absolute prefix");
  if(!valid)return false;
  for(auto* value:{role,prefix}) {
    int size=static_cast<int>(value->size());MPI_Bcast(&size,1,MPI_INT,0,MPI_GLOBAL_COMMUNICATOR);
    value->resize(size);MPI_Bcast(&(*value)[0],size,MPI_CHAR,0,MPI_GLOBAL_COMMUNICATOR);
  }
  return true;
}
std::string ReadBinary(const std::string& path) {
  std::ifstream input(path.c_str(),std::ios::binary);
  return std::string((std::istreambuf_iterator<char>(input)),std::istreambuf_iterator<char>());
}
Result FreshProcess(Mover::ProductionMover selected) {
  Result missing;
  if(GetExecutionContext().inputPath.empty()) {
    missing.status=Status::Skip;missing.message="fresh-process checkpoint case requires paired --test-input roles; run test/run_fresh_process_checkpoint_tests.py";return missing;
  }
  if(!Supported(&missing))return missing;
  Checks checks;std::string role,prefix;
  if(!FreshRequest(&role,&prefix,checks))return Finish(checks,"fresh-process request validated");
  {
    Fixture fixture(checks);
    if(fixture.saved) {
    namespace F=FieldLine::FluxTubeGeometry;
    double magnetic[3];FL::FieldLinesAll[0].GetFirstVertex()->GetMagneticField(magnetic);
    const double expectedFlux=1.25*Run::Active().get().fieldLineSeedAreaM2.value*Vector3D::Length(magnetic);
    if(role=="write") {
      // Deliberately nondefault physical value and provenance: a fresh reader
      // cannot pass by silently rebuilding Phi from its own startup seed area.
      F::SetMagneticFluxWb(0,expectedFlux,"fresh-process diagnostic explicit flux",9);
      if(fixture.Prepare(selected)&&fixture.Step()) {
        // Make the private trajectory cache nontrivial too. This is an explicit
        // post-step query at the authoritative clock, outside any birth/read
        // window; no moving-shock source placement/quadrature is exercised.
        ParticleSource::ShockWave::Tenishev2005::UpdateShockLocation();
        if(checks.Status(CP::SaveCollective(prefix+".checkpoint",1))&&fixture.Step()&&fixture.Step())
          checks.Status(CP::SaveCollective(prefix+".expected",3));
      }
    } else {
      FieldLine::FluxTubeState::State before,after;
      checks.Status(F::CaptureCheckpointState(&before));
      checks.Require(before.records.empty()||before.records[0].generation!=9||
          before.records[0].provenance!="fresh-process diagnostic explicit flux",
          "resume fixture already had writer's private flux state before file load");
      if(fixture.Configure(selected)) {
        std::uint64_t completed=99;
        if(checks.Status(CP::LoadCollective(prefix+".checkpoint",&completed))) {
          checks.Require(completed==1,"fresh process did not restore completed step identity");
          PhysicsCheckpoint::AnalyticShockCache trajectory;
          checks.Status(ParticleSource::ShockWave::Tenishev2005::CaptureAnalyticShockCache(&trajectory));
          checks.Require(trajectory.initialized,"fresh process lost initialized private shock trajectory history");
          checks.Status(F::CaptureCheckpointState(&after));
          checks.Require(!after.records.empty()&&after.records[0].magneticFluxWb==expectedFlux&&
              after.records[0].generation==9&&after.records[0].provenance=="fresh-process diagnostic explicit flux",
              "fresh process did not restore exact nondefault conserved flux/provenance");
          if(fixture.Step()&&fixture.Step()&&checks.Status(CP::SaveCollective(prefix+".resumed",3))&&PIC::ThisThread==0) {
            const auto expected=ReadBinary(prefix+".expected"),actual=ReadBinary(prefix+".resumed");
            checks.Require(!expected.empty()&&expected==actual,
                "fresh-process resumed native mover/source/wave/flux checkpoint differs from uninterrupted final artifact");
          }
        }
      }
    }
    }
  }
  auto result=Finish(checks,"paired fresh-process API checkpoint role completed with exact mover/source/wave/flux state; production resume CLI is separate");
  result.configuration.erase(std::remove_if(result.configuration.begin(),result.configuration.end(),
      [](const std::string& value){return value.find("checkpoint_scope=")==0;}),result.configuration.end());
  result.configuration.push_back("checkpoint_scope=fresh-process-API-role;paired-runner-required-for-continuation-claim");
  result.configuration.push_back("checkpoint_fixture_role="+role+";fresh_process_API_only=true;production_resume_CLI=false");
  result.artifacts={prefix+".checkpoint",prefix+(role=="write"?".expected":".resumed")};return result;
}
Result FreshParker(){return FreshProcess(Mover::ProductionMover::Parker);}
Result FreshDmumu(){return FreshProcess(Mover::ProductionMover::FocusedTransportDiffusion);}
Result FreshMfp(){return FreshProcess(Mover::ProductionMover::FocusedTransportMeanFreePath);}
Result Continuation(Mover::ProductionMover selected,bool requireWaveChange) {
  Result unsupported;if(!Supported(&unsupported))return unsupported;
  Checks checks;std::uint64_t events=0;double changes=0.0;
  {
    Fixture fixture(checks);
    if(fixture.Prepare(selected)) {
      std::vector<double> initialWaves;
      for(int l=0;l<FL::nFieldLine;++l)for(auto* segment=FL::FieldLinesAll[l].GetFirstSegment();segment;segment=segment->GetNext())
        if(segment->Thread==PIC::ThisThread) {
          const auto* energy=segment->GetDatum_ptr(W::CellIntegratedWaveEnergy);
          initialWaves.push_back(energy[0]);initialWaves.push_back(energy[1]);
        }
      if(fixture.Step()) {
      std::size_t position=0;std::uint64_t localChanges=0;
      for(int l=0;l<FL::nFieldLine;++l)for(auto* segment=FL::FieldLinesAll[l].GetFirstSegment();segment;segment=segment->GetNext())
        if(segment->Thread==PIC::ThisThread) {
          const auto* energy=segment->GetDatum_ptr(W::CellIntegratedWaveEnergy);
          if(initialWaves[position]!=energy[0]||initialWaves[position+1]!=energy[1])++localChanges;
          position+=2;
        }
      std::uint64_t globalChanges=0;MPI_Allreduce(&localChanges,&globalChanges,1,MPI_UINT64_T,MPI_SUM,MPI_GLOBAL_COMMUNICATOR);
      std::string first;checks.Status(CP::LocalPhysicsImage(1,&first));
      const auto waveBefore=Turbulence::PICAdapter::LastStepLedger().advectionFaceUpdates;
      if(checks.Status(CP::SaveCollective(fixture.prefix+".checkpoint",1))) {
        bool advanced=fixture.Step()&&fixture.Step();std::string expected;
        if(advanced)advanced=checks.Status(CP::LocalPhysicsImage(3,&expected));
        std::uint64_t completed=99;
        if(advanced&&checks.Status(CP::LoadCollective(fixture.prefix+".checkpoint",&completed))) {
          checks.Require(completed==1,"checkpoint did not restore driver step identity");
          std::string restored;checks.Status(CP::LocalPhysicsImage(1,&restored));
          checks.Require(restored==first,"native physics state differs immediately after exact checkpoint restore");
          if(fixture.Step()&&fixture.Step()) {
            std::string actual;checks.Status(CP::LocalPhysicsImage(3,&actual));
            checks.Require(actual==expected,"resumed native mover/source/wave physics differs from uninterrupted continuation");
          }
        }
      }
      changes=static_cast<double>(globalChanges);
      if(requireWaveChange)checks.Require(waveBefore>0&&globalChanges>0,"checkpoint fixture did not change native waves through production advection");
      if(selected==Mover::ProductionMover::FocusedTransportMeanFreePath) {
        for(int line=0;line<FL::nFieldLine;++line)
          for(auto* segment=FL::FieldLinesAll[line].GetFirstSegment();segment;segment=segment->GetNext())
            for(auto particle=segment->FirstParticleIndex;particle!=-1;particle=PB::GetNext(particle)) {
              std::uint64_t count=0;std::memcpy(&count,PB::GetParticleDataPointer(particle)+SEP::Offset::MfpEventIndex,sizeof(count));events+=count;
            }
        std::uint64_t global=0;MPI_Allreduce(&events,&global,1,MPI_UINT64_T,MPI_SUM,MPI_GLOBAL_COMMUNICATOR);
        checks.Require(global>0,"MFP native checkpoint fixture did not contain a physical scattering event");events=global;
      }
      }
    }
  }
  auto result=Finish(checks,"native PIC mover, synthetic source events, and wave arrays continue bit-for-bit through same-process corrected checkpoint restore; fresh-process persistence remains unqualified");
  result.metrics.push_back({"changed_native_wave_cells",changes,0,">","count"});
  if(selected==Mover::ProductionMover::FocusedTransportMeanFreePath)
    result.metrics.push_back({"native_scattering_events",static_cast<double>(events),0,">","count"});
  return result;
}
Result Parker(){return Continuation(Mover::ProductionMover::Parker,true);}
Result Dmumu(){return Continuation(Mover::ProductionMover::FocusedTransportDiffusion,true);}
Result Mfp(){return Continuation(Mover::ProductionMover::FocusedTransportMeanFreePath,true);}

Result Rejection() {
  Result unsupported;if(!Supported(&unsupported))return unsupported;
  Checks checks;
  {
    Fixture fixture(checks);
    if(fixture.Prepare(Mover::ProductionMover::Parker)&&fixture.Step()) {
      std::string before;checks.Status(CP::LocalPhysicsImage(1,&before));
      if(checks.Status(CP::SaveCollective(fixture.prefix+".checkpoint",1))) {
        if(PIC::ThisThread==0) {
          std::ofstream bad((fixture.prefix+".bad").c_str(),std::ios::binary);bad<<"truncated";
        }
        MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);std::uint64_t completed=99;
        checks.Require(!CP::LoadCollective(fixture.prefix+".bad",&completed).ok()&&completed==99,
            "truncated checkpoint was accepted or changed caller step state");
        std::string after;checks.Status(CP::LocalPhysicsImage(1,&after));
        checks.Require(before==after,"rejected checkpoint partially changed native physics");
        if(WriteInvalidCoordinate(fixture.prefix+".checkpoint",fixture.prefix+".bad",checks)) {
          checks.Require(!CP::LoadCollective(fixture.prefix+".bad",&completed).ok()&&completed==99,
              "CRC-valid nonfinite field-line coordinate was accepted or changed caller state");
          checks.Status(CP::LocalPhysicsImage(1,&after));
          checks.Require(before==after,"malformed native record partially changed physics before rejection");
        }
        if(WritePendingBirth(fixture.prefix+".checkpoint",fixture.prefix+".bad",checks)) {
          const auto pending=CP::LoadCollective(fixture.prefix+".bad",&completed);
          checks.Require(!pending.ok()&&completed==99&&
              pending.message.find("unconsumed scheduled first move")!=std::string::npos,
              "CRC-valid pending birth was accepted into a completed checkpoint or rejected for an unrelated envelope error");
          checks.Status(CP::LocalPhysicsImage(1,&after));
          checks.Require(before==after,"pending-birth checkpoint rejection changed native physics");
        }
        const long int originalParticleCount=PB::GetAllPartNum();
        auto rejectParticle=[&](std::size_t offset,const std::vector<double>& replacement,
            const std::string& reason) {
          if(!WriteInvalidParticleValue(fixture.prefix+".checkpoint",fixture.prefix+".bad",offset,replacement,checks))return;
          const auto rejected=CP::LoadCollective(fixture.prefix+".bad",&completed);
          checks.Require(!rejected.ok()&&completed==99&&PB::GetAllPartNum()==originalParticleCount&&
              rejected.message.find(reason)!=std::string::npos,
              "CRC-valid particle view/hazard did not reach its atomic semantic rejection: "+rejected.message);
          checks.Status(CP::LocalPhysicsImage(1,&after));
          checks.Require(before==after,"rejected particle view/hazard changed phase/source/RNG/background or allocation state");
        };
        // Rebuild a VALID envelope CRC around bad physical bytes. A restored
        // particle may be consumed by generic PIC code before its first mover,
        // so neither NaN X nor a finite one-ULP stale derived view is repairable
        // by waiting for CommitAndAttach. The incoming coordinate is unchanged.
        rejectParticle(_PIC_PARTICLE_DATA__POSITION_OFFSET_,
            {std::numeric_limits<double>::quiet_NaN()},"Cartesian view");
        const auto particleRecord=RepresentativeParticleRecord(checks);
        if(!particleRecord.empty()) {
          double position[3];PB::GetX(position,reinterpret_cast<PB::byte*>(const_cast<char*>(particleRecord.data())));
          rejectParticle(_PIC_PARTICLE_DATA__POSITION_OFFSET_,
              {std::nextafter(position[0],std::numeric_limits<double>::infinity())},"Cartesian view");
        }
        // NaN is explicit unused fresh MFP history, but a completed positive
        // residual process can never persist zero (including negative zero).
        // LoadParticle already enforces this strict contract; the checkpoint
        // validator must reject it before replacing any live records.
        rejectParticle(Offset::MfpOpticalDepth,{0.0},"optical-depth history");
        rejectParticle(Offset::MfpOpticalDepth,{-0.0},"optical-depth history");
        auto* primitiveVertex=FL::FieldLinesAll[0].GetFirstVertex();
        const std::string vertexRecord(primitiveVertex->GetAssociatedDataBufferPointer(),
            primitiveVertex->AssociatedDataLength());
        const long int particleCount=PB::GetAllPartNum();
        auto rejectPrimitive=[&](const PIC::Datum::cDatumStored& datum,
            const std::vector<double>& values,const std::string& expectedReason) {
          if(datum.offset<0){checks.Require(false,"malformed primitive fixture needs an allocated datum");return;}
          if(!WriteInvalidRawFrame(fixture.prefix+".checkpoint",fixture.prefix+".bad",
              vertexRecord,static_cast<std::size_t>(datum.offset),values,checks))return;
          const auto rejected=CP::LoadCollective(fixture.prefix+".bad",&completed);
          checks.Require(!rejected.ok()&&completed==99&&PB::GetAllPartNum()==particleCount,
              "CRC-valid malformed native primitive was accepted or changed caller/allocation state");
          // This proves rejection reached physical preflight, rather than only
          // detecting an intentionally stale CRC or incompatible envelope.
          checks.Require(rejected.message.find(expectedReason)!=std::string::npos,
              "malformed native primitive did not reach its physical validation gate: "+rejected.message);
          checks.Status(CP::LocalPhysicsImage(1,&after));
          checks.Require(before==after,"malformed native primitive changed background/particle/source/RNG state before rejection");
        };
        for(double invalidDensity:{0.0,-1.0,std::numeric_limits<double>::infinity(),
            std::numeric_limits<double>::quiet_NaN(),std::numeric_limits<double>::denorm_min()})
          rejectPrimitive(FL::DatumAtVertexPlasmaDensity,{invalidDensity},"native vertex density");
        rejectPrimitive(FL::DatumAtVertexMagneticField,{0.0,0.0,0.0},"native vertex magnetic field");
        rejectPrimitive(FL::DatumAtVertexPlasmaVelocity,
            {std::numeric_limits<double>::infinity(),0.0,0.0},"native vertex plasma velocity");
        rejectPrimitive(FL::DatumAtVertexPlasmaVelocity,{SpeedOfLight,0.0,0.0},"native vertex plasma velocity");
        auto* vertex=FL::FieldLinesAll[0].GetFirstVertex();double old[3];vertex->GetX(old);
        double changed[3]={std::nextafter(old[0],std::numeric_limits<double>::infinity()),old[1],old[2]};
        vertex->SetX(changed);
        const auto geometry=CP::LoadCollective(fixture.prefix+".checkpoint",&completed);
        vertex->SetX(old);
        checks.Require(!geometry.ok(),"checkpoint accepted a changed physical source line namespace");
        checks.Status(CP::LocalPhysicsImage(1,&after));checks.Require(before==after,
            "geometry rejection changed particle/provider/source history");
        // Saving under an unsupported opaque callback must preserve any existing
        // destination, proving preflight precedes atomic publication.
        if(PIC::ThisThread==0){std::ofstream sentinel((fixture.prefix+".sentinel").c_str());sentinel<<"retain";}
        MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);
        const auto typed=T::Coefficient::ActiveConfiguration();auto opaque=typed;
        opaque.pitchAngle=T::Coefficient::PitchAngleKind::Configured;
        checks.Status(T::Coefficient::SetActiveConfiguration(opaque));
        checks.Require(!CP::SaveCollective(fixture.prefix+".sentinel",1).ok(),"unsupported callback checkpoint was published");
        checks.Status(T::Coefficient::SetActiveConfiguration(typed));
        if(PIC::ThisThread==0){std::ifstream sentinel((fixture.prefix+".sentinel").c_str());std::string text;sentinel>>text;checks.Require(text=="retain","failed checkpoint replaced previous complete destination");}
      }
    }
  }
  return Finish(checks,"truncation, changed geometry, and unsupported callback checkpoints reject before physical-state or destination mutation");
}
Descriptor Make(const char* id,const char* name,TestCallback callback) {
  Descriptor descriptor;descriptor.id=id;descriptor.name=name;descriptor.group="physics-checkpoint-native";
  descriptor.description=name;descriptor.callback=callback;descriptor.initialization=InitializationLevel::FieldLineModel;
  descriptor.runtime=RuntimeClass::Routine;descriptor.supportedBuildModes="native fixed-line analytic MPI-only model with global timestep/weight";
  descriptor.seedPolicy="active campaign; persisted source channels and particle/operator histories";
  descriptor.stateIsolation="original native physics checkpoint restored; output/diagnostic counters excluded";
  return descriptor;
}
} // namespace
std::vector<Descriptor> CheckpointNativeDescriptors() {
  return {Make("CP05C01","Native Parker physics checkpoint continuation",Parker),
          Make("CP05C02","Native Dmumu physics checkpoint continuation",Dmumu),
          Make("CP05C03","Native MFP physics checkpoint continuation",Mfp),
          Make("CP05C04","Native checkpoint transactional rejection",Rejection),
          Make("CP05C05","Native conserved-flux owner persistence",FluxOwnerContinuation)};
}
std::vector<Descriptor> FreshProcessCheckpointDescriptors() {
  std::vector<Descriptor> result;
  const char* names[]={"Fresh-process Parker API continuation","Fresh-process Dmumu API continuation","Fresh-process MFP API continuation"};
  const TestCallback callbacks[]={FreshParker,FreshDmumu,FreshMfp};
  for(unsigned i=0;i<3;++i) {
    auto descriptor=Make(("CP05F0"+std::to_string(i+1)).c_str(),names[i],callbacks[i]);
    descriptor.group="physics-checkpoint-fresh-process";
    // A write role alone does not establish restart equivalence. The paired
    // external launcher supplies explicit role inputs to separate MPI jobs and
    // compares their final artifacts; ordinary routine selection cannot supply
    // that process protocol, so retain it as a separate extended gate.
    descriptor.runtime=RuntimeClass::Extended;
    descriptor.stateIsolation="paired fresh native processes; test input chooses explicit role/prefix; original local state restored after role";
    descriptor.seedPolicy="fixed active model campaign; no converted test-seed or single-process repeat claim";
    result.push_back(descriptor);
  }
  return result;
}
} } // namespace SEP::Testing
