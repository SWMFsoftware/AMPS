// Actual Moon clock/source-boundary/mover hooks on real PIC lists/allocator.
// Numerical chronology controls ONLY: no measured source, reaction, recoil,
// density/brightness qualification or external observation is supplied.
#include "pic.h"
#include "Moon.h"
#include "MoonInput.h"
#include "MoonTransportClock.h"
#include "MoonSurfaceEnvironment.h"
#include "MoonReactionController.h"
#include <array>
#include <cmath>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>
namespace T=PIC::Mover::ForceCutSurface;
namespace C=Moon::Chemistry::Controller;
namespace R=Moon::Species::ReactionExecution;
using Node=cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>;
namespace {
int Order=2,Checks=0,Failures=0;long int TagOffset=-1;
std::ofstream CheckFile,BirthFile;std::vector<Node*> Leaves;
std::ofstream DaughterFile;
Moon::Chemistry::ReactionLedger ReactionLedger;
double CompositionEventEpoch=NAN,CompositionEventRemaining=NAN;
int CompositionProductDraws=0;
std::array<double,3> Start{{0,0,0}};
void Require(bool valid,const std::string& message){if(!valid)throw std::runtime_error(message);}
void Check(bool valid,const char* name) {
  int local=valid?1:0,all=0;MPI_Allreduce(&local,&all,1,MPI_INT,MPI_MIN,MPI_GLOBAL_COMMUNICATOR);
  ++Checks;if(!all)++Failures;
  if(PIC::ThisThread==0){CheckFile<<name<<','<<(all?"PASS":"FAIL")<<'\n';CheckFile.flush();}
}
double Sum(double value){double out=0;MPI_Allreduce(&value,&out,1,MPI_DOUBLE,MPI_SUM,MPI_GLOBAL_COMMUNICATOR);return out;}
int Tag(long int ptr){int tag=-1;std::memcpy(&tag,PIC::ParticleBuffer::GetParticleDataPointer(ptr)+TagOffset,sizeof(tag));return tag;}
bool ResolveControlState(long int particle,Node*,void*,std::size_t* state,double* internal,std::string*) {
  // All existing chronology cohorts are an EXPLICIT inert numerical state.
  // Only the separate tag7 composition control has a constant reaction. Both
  // map to compiled Na; this resolver never guesses a physical excited state.
  *state=Tag(particle)==7?1:0;*internal=0;return true;
}
bool ControlUniform(C::DrawPurpose purpose,long int,std::size_t,void*,double* uniform,std::string*) {
  if(purpose==C::DrawPurpose::SourceClock)*uniform=std::exp(-.4);
  else if(purpose==C::DrawPurpose::Channel)*uniform=.5;
  else {*uniform=.3;++CompositionProductDraws;}
  return true;
}
bool CompositionKinematics(const C::EventContext& event,R::KinematicPolicy* policy,std::string* error) {
  // Nominal-heavy-mass co-moving-bath zero-Q MECHANICAL control only. The
  // incoming bath velocity is explicitly the parent velocity, BEFORE product
  // moments are computed. No residual is fitted, recoil distribution assumed
  // or real CEX channel/measurement installed by this test.
  policy->qualified=true;policy->identity="COMMON-CLOCK-CLONE-ZERO-Q-CONTROL";
  policy->evidenceIdentity="ANALYTICAL-STORAGE-ONLY-NOT-LUNAR-CEX";policy->eventFrame="LSO";
  policy->maximumClassicalSpeedMS=1e5;policy->bounds.modelMassKg={0,1e-13};
  policy->bounds.momentumKgMS={0,1e-13};policy->bounds.energyJ={0,1e-13};policy->bath.supplied=true;
  const double helium=event.catalog->at(4).modelMassKg;policy->bath.modelMassChangeKg=-helium;
  R::ProductKinematics ion,neutral;ion.term=0;neutral.term=1;
  ion.referenceWeight=event.node->block->GetLocalParticleWeight(_NA_PLUS_SPEC_);
  neutral.referenceWeight=event.node->block->GetLocalParticleWeight(_HE_SPEC_);
  for(int k=0;k<3;++k) {
    const double incoming=event.parent.velocityMS[k];ion.velocityMS[k]=incoming;neutral.velocityMS[k]=incoming;
    policy->bath.momentumChangeKgMS[k]=-helium*incoming;
    policy->bath.energyChangeJ-=.5*helium*incoming*incoming;
  }
  double uniform=0;
  if(!event.openUniform(C::DrawPurpose::ProductClock,event.particle,event.channelIndex,event.userData,&uniform,error) ||
      !R::OpticalDepthClock::FromOpenUniform(uniform,&ion.nextClock,error))return false;
  if(!event.openUniform(C::DrawPurpose::ProductClock,event.particle,event.channelIndex,event.userData,&uniform,error) ||
      !R::OpticalDepthClock::FromOpenUniform(uniform,&neutral.nextClock,error))return false;
  policy->products={ion,neutral};CompositionEventEpoch=event.physicalEpochS;
  CompositionEventRemaining=event.interval.remainingS();return true;
}
void InstallCompositionControl() {
  const double na=PIC::MolecularData::GetMass(_NA_SPEC_),nai=PIC::MolecularData::GetMass(_NA_PLUS_SPEC_);
  const double he=PIC::MolecularData::GetMass(_HE_SPEC_),hei=PIC::MolecularData::GetMass(_HE_PLUS_SPEC_);
  Require(na==nai && he==hei,"analytical clone control requires the declared actual nominal heavy mass convention");
  C::Configuration configuration;configuration.identity="COMMON-CLOCK-SINGLE-EVENT-CLONE-CONTROL";configuration.eventFrame="LSO";
  configuration.resolveInitialState=ResolveControlState;configuration.openUniform=ControlUniform;configuration.ledger=&ReactionLedger;
  configuration.catalog={
    {"Na","inert-numerical-control",_NA_SPEC_,{{"Na",1}},0,na},
    {"Na","reactive-numerical-control",_NA_SPEC_,{{"Na",1}},0,na},
    {"Na+","terminal-numerical-control",_NA_PLUS_SPEC_,{{"Na",1}},1,nai},
    {"He","terminal-numerical-control",_HE_SPEC_,{{"He",1}},0,he},
    {"He+","declared-co-moving-bath",_HE_PLUS_SPEC_,{{"He",1}},1,hei}};
  C::Channel channel;channel.parentCatalogSpecies=1;channel.constantRateS=1;
  channel.scope=C::HazardScope::ConstantHistory;channel.hazardEvidenceIdentity="ANALYTICAL-CONSTANT-LAMBDA-ONE-PER-SECOND";
  channel.kinematicProviderIdentity="COMMON-CLOCK-CLONE-ZERO-Q-CONTROL";channel.kinematics=CompositionKinematics;
  channel.reaction={"COMMON-CLOCK-CLONE-CONTROL",{{1,1,R::Population::Tracked},{4,1,R::Population::Bath}},
      {{2,1,R::Population::Tracked},{3,1,R::Population::Tracked}}};configuration.channels={channel};
  std::string error;Require(C::Install(configuration,&error),error);
}
void Collect(Node* node){if(node->lastBranchFlag()==_BOTTOM_BRANCH_TREE_){Leaves.push_back(node);return;}for(int k=0;k<8;++k)if(node->downNode[k])Collect(node->downNode[k]);}
double Resolution(double*){return 1e6;}
double Load(Node*){return 1;}
void SetAllDt(double dt) {
  for(Node* node:Leaves)if(node->block && node->Thread==PIC::ThisThread)
    for(int species=0;species<PIC::nTotalSpecies;++species)node->block->SetLocalTimeStep(dt,species);
}
template<class Function>void ForParticle(Function function) {
  for(unsigned n=0;n<PIC::DomainBlockDecomposition::nLocalBlocks;++n){Node* node=PIC::DomainBlockDecomposition::BlockTable[n];
    if(!node || !node->block || node->Thread!=PIC::ThisThread)continue;
    for(int cell=0;cell<_BLOCK_CELLS_X_*_BLOCK_CELLS_Y_*_BLOCK_CELLS_Z_;++cell)
      for(long int ptr=node->block->FirstCellParticleTable[cell];ptr!=-1;ptr=PIC::ParticleBuffer::GetNext(ptr))function(ptr,node,cell);
  }
}
long int Source(Node* node,const double* x,int tag) {
  const long int ptr=PIC::ParticleBuffer::GetNewParticle();auto* data=PIC::ParticleBuffer::GetParticleDataPointer(ptr);
  double position[3]={x[0],x[1],x[2]},v[3]={730,-410,95};PIC::ParticleBuffer::SetI(_NA_SPEC_,data);
  PIC::ParticleBuffer::SetX(position,data);PIC::ParticleBuffer::SetV(v,data);
  PIC::ParticleBuffer::SetIndividualStatWeightCorrection(1/node->block->GetLocalParticleWeight(_NA_SPEC_),data);
  std::memcpy(data+TagOffset,&tag,sizeof(tag));T::ResetParticleState(data);
  // No finite epoch is assigned here: the ACTUAL source macro must supply
  // the advancing physical birth through the normal selected Moon wrapper.
  Require(!std::isfinite(T::GetParticleEpochS(ptr)),"new/recycled allocator slot did not reset source epoch");return ptr;
}
void RecordBirth(long int ptr,int tag,double expected,double remaining) {
  const auto* clock=Moon::GetCommonTransportClock();const auto* step=clock->step();
  BirthFile<<tag<<','<<step->sequence<<','<<(expected-1000)<<','<<expected<<','
    <<T::GetTimeContext(ptr).InitialEpochS<<','<<remaining<<','<<PIC::ThisThread<<'\n';BirthFile.flush();
}
void AdvanceNormalMover() {double move=0,exchange=0;PIC::TimeStepInternal::ExecutionTrackDefault(move,exchange);}
bool Population(double epoch,int expectedCount) {
  bool valid=true;double count=0;
  ForParticle([&](long int ptr,Node* node,int){++count;
    const double number=node->block->GetLocalParticleWeight(PIC::ParticleBuffer::GetI(ptr))*PIC::ParticleBuffer::GetIndividualStatWeightCorrection(ptr);
    valid=valid && Moon::TransportTime::EpochProjectionAgrees(T::GetParticleEpochS(ptr),epoch) &&
      T::GetTimeContext(ptr).RemainingTimeS==0 && std::abs(number-1)<1e-13;
  });
  const double total=Sum(count);return valid && total==expectedCount;
}
void Setup() {
  PIC::Init_BeforeParser();Moon::Runtime::Configuration c;std::string error;
  c.particleTracingMover=Order==1?"force-cut-surface-first-order":"force-cut-surface-second-order";
  c.tracingClockMode="common-time";c.simulationEnvironmentMode="frozen";c.tracingFieldFrame="lso";
  c.tracingRelativeTolerance=1e-8;c.tracingAbsolutePositionToleranceM=1e-8;c.tracingAbsoluteVelocityToleranceMS=1e-8;
  c.tracingGeometryToleranceM=1e-8;c.tracingTimeToleranceS=1e-12;c.tracingMinimumSubstepS=1e-12;c.tracingMaximumSubstepS=.125;
  c.tracingMaximumSubsteps=10000;c.tracingMaximumRefinements=64;c.tracingMaximumInteractions=1000;c.tracingMaximumSamePointInteractions=32;c.tracingAdaptive=true;
  Require(Moon::Runtime::InstallConfiguration(c,&error),error);Moon::RegisterParticleTracing();Moon::ConfigureParticleTracingFields();
  Check(Moon::CommonTransportTimeSelected() && Moon::Runtime::GetConfiguration().tracingClockMode=="common-time",
      "explicit_common_mode_and_real_tracer_selected");
  PIC::ParticleBuffer::RequestDataStorage(TagOffset,sizeof(int));
  InstallCompositionControl();
  // Explicit analytical ET is the numerical clock oracle. Normal application
  // initialization instead passes its real CSPICE-converted start UTC.
  Require(Moon::InitializeCommonTransportTime(1000,&error),error);
  Moon::OrbitalMotion::et=1000;
  double lower[3]={1.8e6,-5e5,-5e5},upper[3]={2.8e6,5e5,5e5};PIC::Mesh::mesh->AllowBlockAllocation=false;
  PIC::Mesh::mesh->init(lower,upper,Resolution);Require(PIC::Mesh::mesh->splitTreeNode(PIC::Mesh::mesh->rootTree),"root split failed");
  for(int k=0;k<8;++k)Require(PIC::Mesh::mesh->splitTreeNode(PIC::Mesh::mesh->rootTree->downNode[k]),"child split failed");
  PIC::Mesh::mesh->SetNodeNeibResolutionLevelLimit();PIC::Mesh::mesh->CreateBottomBranchNodeList(PIC::Mesh::mesh->rootTree);
  PIC::Mesh::mesh->SetParallelLoadMeasure(Load);PIC::Mesh::mesh->CreateNewParallelDistributionLists();
  PIC::Mesh::initCellSamplingDataBuffer();PIC::Mesh::mesh->AllowBlockAllocation=true;PIC::Mesh::mesh->AllocateTreeBlocks();
  PIC::Mesh::mesh->InitCellMeasure();PIC::Init_AfterParser();PIC::DomainBlockDecomposition::UpdateBlockTable();Collect(PIC::Mesh::mesh->rootTree);
  for(Node* node:Leaves)if(node->block)for(int species=0;species<PIC::nTotalSpecies;++species) {
    node->block->SetLocalTimeStep(1,species);node->block->SetLocalParticleWeight(2+species+.01*(node->xmin[0]-1.8e6)/1e6,species);
  }
  bool found=false;
  for(Node* a:Leaves)for(Node* b:Leaves)if(!found && a->xmax[0]==b->xmin[0] && a->xmin[1]==b->xmin[1] && a->xmin[2]==b->xmin[2] &&
      (PIC::nTotalThreads==1 || a->Thread!=b->Thread)) {
    Start={{a->xmax[0]-730*.2,.5*(a->xmin[1]+a->xmax[1]),.5*(a->xmin[2]+a->xmax[2])}};found=true;
  }
  Require(found,"fixture has no requested actual ownership-crossing face");PIC::ParticleBuffer::Init(64);
  Exosphere::xObjectRadial=_AU_;const double sun[3]={_AU_,0,0},earth[3]={-3.8e8,2e7,5e6};
  for(int k=0;k<3;++k){Exosphere::xSun_SO[k]=sun[k];Exosphere::xEarth_SO[k]=earth[k];Exosphere::RotationVector_SO_FROZEN[k]=0;Exosphere::vObject_SO_FROZEN[k]=0;}
  T::Acceleration=nullptr;PIC::Mover::SetUserDefinedParticleMover(Moon::MoveParticle);
  Check(Moon::GetCommonTransportClock() && Moon::GetCommonTransportClock()->originEtS()==1000 &&
      Moon::SurfaceEnvironment::ProductionSourceEpochS()==1000 && !Moon::GetCommonTransportClock()->active(),
      "common_storage_origin_and_frozen_coefficients_initialized");
}
void Cases() {
  std::string error;bool localBad=false;
  for(Node* node:Leaves)if(!localBad && node->block && node->Thread==PIC::ThisThread) {
    node->block->SetLocalTimeStep(.5,_NA_SPEC_);localBad=true;
  }
  Check(!Moon::BeginCommonTransportStep(1000,&error) && error.find("identical constant ACTUAL dt")!=std::string::npos &&
      !Moon::GetCommonTransportClock()->active() && Moon::GetCommonTransportClock()->elapsedS()==0,
      "unequal_actual_dt_refused_before_physical_step");
  SetAllDt(1);
  Check(!Moon::BeginCommonTransportStep(1001,&error) && !Moon::GetCommonTransportClock()->active(),
      "changed_external_epoch_refused_before_begin");
  Require(Moon::BeginCommonTransportStep(1000,&error),error);
  Check(Moon::GetCommonTransportClock()->step()->relativeBeginS==0 && Moon::GetCommonTransportClock()->step()->relativeEndS==1,
      "first_common_step_begins_with_immutable_relative_endpoints");
  Check(!Moon::BeginCommonTransportStep(1000,&error) && Moon::GetCommonTransportClock()->step()->sequence==1,
      "nested_begin_cannot_replace_active_step");
  Node* initial=PIC::Mesh::mesh->findTreeNode(Start.data());const int owner=initial->Thread;
  bool sourceValid=true;
  if(owner==PIC::ThisThread) {
    for(int tag:{0,9}) {
      const double offset=tag==0?.25:.75,remaining=1-offset;
      const double x[3]={Start[0],Start[1]+10*tag,Start[2]};const long int ptr=Source(initial,x,tag);
      PIC::Mover::MovingSurfaceInjectionTimeOffsetS=offset;
      sourceValid=sourceValid && Moon::SurfaceEnvironment::ProductionInjectionEpochS()==1000+offset;
      _PIC_PARTICLE_MOVER__MOVE_PARTICLE_BOUNDARY_INJECTION_(ptr,remaining,initial,true);
      sourceValid=sourceValid && T::GetTimeContext(ptr).InitialEpochS==1000+offset;
      RecordBirth(ptr,tag,1000+offset,remaining);PIC::Mover::MovingSurfaceInjectionTimeOffsetS=NAN;
    }
  }
  Check(sourceValid,"scoped_offset_source_uses_actual_advancing_transport_birth");
  AdvanceNormalMover();
  Check(Population(1001,2),"first_normal_mpi_continuations_complete");
  bool corrupted=false;
  ForParticle([&](long int ptr,Node*,int){if(Tag(ptr)==0){T::SetParticleEpochS(ptr,1002);corrupted=true;}});
  Check(!Moon::CommitCommonTransportStep(&error) && Moon::GetCommonTransportClock()->active() &&
      Moon::GetCommonTransportClock()->elapsedS()==0,"bad_endpoint_commit_refused_without_clock_publication");
  // Deliberately corrupted diagnostic bytes are repaired ONLY in this control;
  // production does not retry a failed completion after physical publication.
  if(corrupted)ForParticle([&](long int ptr,Node*,int){if(Tag(ptr)==0)T::SetParticleEpochS(ptr,1001);});
  Require(Moon::CommitCommonTransportStep(&error),error);
  const bool firstPopulation=Population(1001,2);
  Check(Moon::GetCommonTransportClock()->completedSteps()==1 && Moon::GetCommonTransportClock()->elapsedS()==1 &&
      firstPopulation,"first_commit_accepts_only_completed_endpoint_population");
  // Remove one explicitly numerical carrier for an allocator lifecycle check,
  // not a physical loss channel. The same slot is reused by a later source.
  Node* recycledNode=nullptr;long int recycled=-1;std::array<double,3> recycledX{{0,0,0}};
  ForParticle([&](long int ptr,Node* node,int){if(Tag(ptr)==9){recycledNode=node;recycled=ptr;PIC::ParticleBuffer::GetX(recycledX.data(),ptr);}});
  if(recycled!=-1) {
    int i,j,k;PIC::Mesh::mesh->FindCellIndex(recycledX.data(),i,j,k,recycledNode,false);
    long int& head=recycledNode->block->FirstCellParticleTable[i+_BLOCK_CELLS_X_*(j+_BLOCK_CELLS_Y_*k)];
    const long int previous=PIC::ParticleBuffer::GetPrev(recycled),next=PIC::ParticleBuffer::GetNext(recycled);
    if(previous==-1)head=next;else PIC::ParticleBuffer::SetNext(next,previous);
    if(next!=-1)PIC::ParticleBuffer::SetPrev(previous,next);PIC::ParticleBuffer::DeleteParticle(recycled);
  }
  SetAllDt(2);
  Check(!Moon::BeginCommonTransportStep(1000,&error) && !Moon::GetCommonTransportClock()->active() && Moon::GetCommonTransportClock()->elapsedS()==1,
      "variable_common_duration_refused_without_advancing_clock");
  SetAllDt(1);Require(Moon::BeginCommonTransportStep(1000,&error),error);
  Check(Moon::GetCommonTransportClock()->step()->beginEtS==1001 && Moon::GetCommonTransportClock()->step()->relativeEndS==2,
      "second_common_step_begins_after_accepted_transport_time");
  bool recycledSource=true;
  if(recycled!=-1) {
    const long int ptr=Source(recycledNode,recycledX.data(),1);recycledSource=ptr==recycled;
    _PIC_PARTICLE_MOVER__MOVE_PARTICLE_BOUNDARY_INJECTION_FACE_(ptr,.5,recycledNode,true,0);
    recycledSource=recycledSource && T::GetTimeContext(ptr).InitialEpochS==1001.5;
    RecordBirth(ptr,1,1001.5,.5);
  }
  Check(recycledSource,"recycled_later_source_uses_fresh_clock_target_and_face_birth");
  AdvanceNormalMover();Require(Moon::CommitCommonTransportStep(&error),error);
  Check(Population(1002,2) && Moon::GetCommonTransportClock()->elapsedS()==2,
      "second_normal_population_completes_common_endpoint");
  Require(Moon::BeginCommonTransportStep(1000,&error),error);AdvanceNormalMover();Require(Moon::CommitCommonTransportStep(&error),error);
  Check(Moon::GetCommonTransportClock()->completedSteps()==3 && Moon::GetCommonTransportClock()->elapsedS()==3,
      "third_step_preserves_surviving_relative_chronology");
  Check(Population(1003,2),"all_packets_finish_common_physical_endpoint");
  double e[3],b[3];const bool fields=Moon::GetParticleTracingFields(e,b);
  Check(Moon::OrbitalMotion::et==1000 && Moon::SurfaceEnvironment::ProductionSourceEpochS()==1000 && fields &&
      e[0]==0 && e[1]==0 && e[2]==0 && b[0]==0 && b[1]==0 && b[2]==0,
      "external_epoch_and_prescribed_fields_remain_frozen");
  // A separate declared numerical cohort verifies COMPOSITION of common time,
  // conservative controller, real parent-byte clone and daughter queue. The
  // preceding chronology cohorts remain inert; no physical source/loss law is
  // inferred from these administrative fixture deletions.
  PIC::ParticleBuffer::DeleteAllParticles();
  Require(Moon::BeginCommonTransportStep(1000,&error),error);
  if(owner==PIC::ThisThread) {
    const long int source=Source(initial,Start.data(),7);
    PIC::Mover::MovingSurfaceInjectionTimeOffsetS=.25;
    _PIC_PARTICLE_MOVER__MOVE_PARTICLE_BOUNDARY_INJECTION_(source,.75,initial,true);
    PIC::Mover::MovingSurfaceInjectionTimeOffsetS=NAN;
  }
  AdvanceNormalMover();Require(Moon::CommitCommonTransportStep(&error),error);
  bool daughters=true;double daughterCount=0,number=0,ions=0,helium=0;
  ForParticle([&](long int particle,Node* node,int) {
    ++daughterCount;const int species=PIC::ParticleBuffer::GetI(particle);
    ions+=species==_NA_PLUS_SPEC_;helium+=species==_HE_SPEC_;
    const double represented=node->block->GetLocalParticleWeight(species)*PIC::ParticleBuffer::GetIndividualStatWeightCorrection(particle);
    number+=represented;C::History history;Require(C::ReadHistory(particle,&history,&error),error);
    const auto time=T::GetTimeContext(particle);const double epoch=T::GetParticleEpochS(particle);
    daughters=daughters && (history.catalogSpecies==2 || history.catalogSpecies==3) && !history.due &&
      history.clock.valid() && !history.clock.consumed() && history.interval.remainingS()==0 &&
      time.RemainingTimeS==0 && epoch==1004 && std::abs(represented-1)<1e-13;
    DaughterFile<<history.catalogSpecies<<','<<species<<','<<epoch<<','<<time.RemainingTimeS<<','<<represented<<','
      <<CompositionEventEpoch<<','<<CompositionEventRemaining<<'\n';
  });
  bool ledger=true;for(const auto& accepted:ReactionLedger.events)
    ledger=ledger && accepted.ledger.integerConservationVerified && accepted.ledger.kinematicConservationVerified &&
      std::abs(accepted.ledger.parentNumber-1)<1e-13 && std::abs(accepted.ledger.productNumber-2)<1e-13 &&
      Moon::TransportTime::EpochProjectionAgrees(CompositionEventEpoch,1003.65) && std::abs(CompositionEventRemaining-.35)<1e-12;
  // Integer topology is EXACT. These predeclared floating N bounds concern
  // only W*(N/W) reconstruction already exercised by the real backend; they
  // do not change a source law, branching probability or physical inventory.
  const double allDaughters=Sum(daughterCount),allNumber=Sum(number),allEvents=Sum(static_cast<double>(ReactionLedger.events.size())),allClocks=Sum(CompositionProductDraws),allIons=Sum(ions),allHelium=Sum(helium);
  Check(daughters && ledger && allDaughters==2 && allIons==1 && allHelium==1 && std::abs(allNumber-2)<1e-12 && allEvents==1 && allClocks==2 &&
      Moon::GetCommonTransportClock()->completedSteps()==4 && Moon::GetCommonTransportClock()->elapsedS()==4,
      "real_conservative_daughters_inherit_common_target_and_finish_unused_interval");
  Check(T::CountPending()==0,"no_unfinished_physical_continuation");
  PIC::ParticleBuffer::DeleteAllParticles();Check(PIC::ParticleBuffer::GetAllPartNum()==0,"no_native_particle_slot_leak");
}
}
int main(int argc,char** argv) {
  PIC::InitMPI();SetExitErrorCode(__LINE__,1);if(argc!=3){MPI_Finalize();return 2;}
  Order=std::stoi(argv[1]);const std::string prefix=std::string(argv[2])+".rank"+std::to_string(PIC::ThisThread);
  if(PIC::ThisThread==0){CheckFile.open(std::string(argv[2])+".checks.csv");CheckFile<<"check,status\n";}
  BirthFile.open(prefix+".births.csv");BirthFile<<std::setprecision(21)<<"source_tag,step,relative_birth_s,expected_birth_et_s,actual_birth_et_s,remaining_s,source_owner\n";
  DaughterFile.open(prefix+".daughters.csv");DaughterFile<<std::setprecision(21)<<"catalog_state,compiled_species,epoch_s,remaining_s,represented_number,event_epoch_s,event_remaining_s\n";
  int failure=0;try{Setup();Cases();}catch(const std::exception& error){failure=1;std::cerr<<"ERROR native common clock: "<<error.what()<<'\n';}
  int bad=0,failed=0;MPI_Allreduce(&failure,&bad,1,MPI_INT,MPI_MAX,MPI_GLOBAL_COMMUNICATOR);MPI_Allreduce(&Failures,&failed,1,MPI_INT,MPI_MAX,MPI_GLOBAL_COMMUNICATOR);
  BirthFile.flush();DaughterFile.flush();if(PIC::ThisThread==0){CheckFile.flush();std::cout<<"M3 common transport clock order="<<Order<<" checks="<<Checks<<" status="<<(bad?"ERROR":failed?"FAIL":"PASS")<<'\n';}
  MPI_Finalize();return bad?2:failed?1:0;
}
