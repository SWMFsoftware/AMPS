// Actual PIC walker/Moon mover/controller controls. The three source states
// and nominal-heavy-mass zero-Q kinematics are ANALYTICAL numerical fixtures,
// not physical Na/CEX rate/energy/angle evidence and never observations.
#include "pic.h"
#include "Moon.h"
#include "MoonInput.h"
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

namespace C=Moon::Chemistry::Controller;
namespace R=Moon::Species::ReactionExecution;
namespace T=PIC::Mover::ForceCutSurface;
using Node=cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>;
namespace {
int Order=2,Failures=0,Checks=0;bool Capacity=false,FutureClock=false;
long int TagOffset=-1;
Moon::Chemistry::ReactionLedger Ledger;
struct Counters {int resolver=0,source=0,channel=0,product=0;};
Counters Counts;
std::ofstream CheckFile,EventFile,ParticleFile;
std::vector<double> CandidateEpochs,CandidateRemaining;
std::vector<int> CandidateTags,CandidateChannels,CandidateOwners;
std::string Phase="main";
std::array<double,3> Start{{0,0,0}};
std::vector<Node*> Leaves;
void Require(bool value,const std::string& error){if(!value)throw std::runtime_error(error);}
void Check(bool condition,const char* name) {
  int local=condition?1:0,all=0;MPI_Allreduce(&local,&all,1,MPI_INT,MPI_MIN,MPI_GLOBAL_COMMUNICATOR);
  ++Checks;if(!all)++Failures;
  if(PIC::ThisThread==0){CheckFile<<name<<','<<(all?"PASS":"FAIL")<<'\n';CheckFile.flush();}
}
double Sum(double value){double all=0;MPI_Allreduce(&value,&all,1,MPI_DOUBLE,MPI_SUM,MPI_GLOBAL_COMMUNICATOR);return all;}
int Tag(long int ptr){int value=-1;std::memcpy(&value,PIC::ParticleBuffer::GetParticleDataPointer(ptr)+TagOffset,sizeof(value));return value;}
void SetTag(long int ptr,int value){std::memcpy(PIC::ParticleBuffer::GetParticleDataPointer(ptr)+TagOffset,&value,sizeof(value));}
bool Resolve(long int ptr,Node*,void*,std::size_t* state,double* internal,std::string* error) {
  ++Counts.resolver;const int tag=Tag(ptr);if(tag<0 || tag>2){*error="analytical source tag is absent";return false;}
  *state=static_cast<std::size_t>(tag);*internal=0;return true;
}
bool Uniform(C::DrawPurpose purpose,long int ptr,std::size_t,void*,double* value,std::string*) {
  const int tag=Tag(ptr);
  if(purpose==C::DrawPurpose::SourceClock){++Counts.source;*value=FutureClock?
      std::nextafter(std::exp(-1.0),0.0):std::exp(-(tag==0?.4:tag==1?.6:1.0));}
  else if(purpose==C::DrawPurpose::Channel){++Counts.channel;*value=tag==1?.8:.2;}
  else {++Counts.product;*value=Counts.product%2?.2:.8;}
  return true;
}
bool Kinematics(const C::EventContext& event,R::KinematicPolicy* policy,std::string* error) {
  const int tag=Tag(event.particle);const bool reverse=event.channelIndex%2;
  policy->qualified=true;policy->identity="ANALYTICAL-policy-"+std::to_string(event.channelIndex);
  policy->evidenceIdentity="STORAGE-CONTROLLER-ONLY-NOT-CEX";policy->eventFrame="LSO";
  policy->maximumClassicalSpeedMS=1e5;policy->bounds.modelMassKg={0,1e-13};
  policy->bounds.momentumKgMS={0,1e-13};policy->bounds.energyJ={0,1e-13};policy->bath.supplied=true;
  const double sodium=event.catalog->at(event.parent.species).modelMassKg;
  const double helium=event.catalog->at(5).modelMassKg;
  policy->bath.modelMassChangeKg=-helium;
  R::ProductKinematics a,b;a.term=0;b.term=1;
  a.referenceWeight=event.node->block->GetLocalParticleWeight(_NA_PLUS_SPEC_);
  b.referenceWeight=event.node->block->GetLocalParticleWeight(_HE_SPEC_);
  for(int k=0;k<3;++k) {
    // Declare incoming bath velocity FIRST. Its signed moments are not fitted
    // from a residual. Reversing relative motion is an explicit elastic
    // numerical control with nominal masses, not an angular/CEX prescription.
    const double parent=event.parent.velocityMS[k],bath=reverse?-parent:parent;
    policy->bath.momentumChangeKgMS[k]=-helium*bath;
    policy->bath.energyChangeJ-=.5*helium*bath*bath;
    if(reverse){a.velocityMS[k]=((sodium-helium)*parent+2*helium*bath)/(sodium+helium);
      b.velocityMS[k]=(2*sodium*parent+(helium-sodium)*bath)/(sodium+helium);}
    else {a.velocityMS[k]=parent;b.velocityMS[k]=bath;}
  }
  double u=0;
  if(!event.openUniform(C::DrawPurpose::ProductClock,event.particle,event.channelIndex,event.userData,&u,error) ||
      !R::OpticalDepthClock::FromOpenUniform(u,&a.nextClock,error))return false;
  if(!event.openUniform(C::DrawPurpose::ProductClock,event.particle,event.channelIndex,event.userData,&u,error) ||
      !R::OpticalDepthClock::FromOpenUniform(u,&b.nextClock,error))return false;
  policy->products={a,b};
  CandidateEpochs.push_back(event.physicalEpochS);CandidateRemaining.push_back(event.interval.remainingS());
  CandidateTags.push_back(tag);CandidateChannels.push_back(static_cast<int>(event.channelIndex));CandidateOwners.push_back(PIC::ThisThread);
  EventFile<<Phase<<','<<tag<<','<<event.channelIndex<<','<<event.physicalEpochS<<','<<event.interval.remainingS()<<','<<PIC::ThisThread<<'\n';EventFile.flush();
  return true;
}
double UnitLoad(Node*){return 1;}
double Resolution(double*){return 1e6;}
void DeclaredZeroAcceleration(double* acceleration,int,long int,const double*,const double*,Node*,const T::cTimeContext&) {
  // Independent straight-line geometry/endpoint control only. Main, boundary
  // and pending-clock cohorts below use the ACTUAL normal Moon force callback.
  // This exact zero vector is not a duplicate lunar force law or observation.
  for(int k=0;k<3;++k)acceleration[k]=0;
}
void Collect(Node* node){if(node->lastBranchFlag()==_BOTTOM_BRANCH_TREE_){Leaves.push_back(node);return;}
  for(int k=0;k<8;++k)if(node->downNode[k])Collect(node->downNode[k]);}
C::Configuration Configuration() {
  C::Configuration c;c.identity="DECLARED-ANALYTICAL-CONSTANT-HISTORY-CONTROL-v1";c.eventFrame="LSO";
  c.resolveInitialState=Resolve;c.openUniform=Uniform;c.ledger=&Ledger;
  const double na=PIC::MolecularData::GetMass(_NA_SPEC_),nai=PIC::MolecularData::GetMass(_NA_PLUS_SPEC_);
  const double he=PIC::MolecularData::GetMass(_HE_SPEC_),hei=PIC::MolecularData::GetMass(_HE_PLUS_SPEC_);
  Require(na==nai && he==hei,"nominal-heavy-mass control requires declared actual neutral/ion mass equality");
  for(int tag=0;tag<3;++tag)c.catalog.push_back({"Na","source-control-"+std::to_string(tag),_NA_SPEC_,{{"Na",1}},0,na});
  c.catalog.push_back({"Na+","terminal-control",_NA_PLUS_SPEC_,{{"Na",1}},1,nai});
  c.catalog.push_back({"He","terminal-control",_HE_SPEC_,{{"He",1}},0,he});
  c.catalog.push_back({"He+","incoming-control-bath",_HE_PLUS_SPEC_,{{"He",1}},1,hei});
  for(int tag=0;tag<3;++tag)for(int branch=0;branch<2;++branch) {
    C::Channel channel;channel.parentCatalogSpecies=tag;channel.constantRateS=branch?.75:.25;
    channel.scope=C::HazardScope::ConstantHistory;channel.hazardEvidenceIdentity="ANALYTICAL-LAMBDA-NOT-CEX";
    channel.kinematicProviderIdentity="ANALYTICAL-policy-"+std::to_string(2*tag+branch);channel.kinematics=Kinematics;
    channel.reaction={"CONTROL-source"+std::to_string(tag)+"-branch"+std::to_string(branch),
      {{static_cast<std::size_t>(tag),1,R::Population::Tracked},{5,1,R::Population::Bath}},
      {{3,1,R::Population::Tracked},{4,1,R::Population::Tracked}}};c.channels.push_back(channel);
  }
#ifdef MOON_CONTROLLER_CONFIGURATION_CALLBACK
  MOON_CONTROLLER_CONFIGURATION_CALLBACK(c);
#endif
  return c;
}
void Setup() {
  PIC::Init_BeforeParser();Moon::Runtime::Configuration c;std::string error;
  c.particleTracingMover=Order==1?"force-cut-surface-first-order":"force-cut-surface-second-order";
  c.simulationEnvironmentMode="frozen";c.tracingFieldFrame="lso";c.tracingRelativeTolerance=1e-8;
  c.tracingAbsolutePositionToleranceM=1e-8;c.tracingAbsoluteVelocityToleranceMS=1e-8;
  c.tracingGeometryToleranceM=1e-8;c.tracingTimeToleranceS=1e-12;c.tracingMinimumSubstepS=1e-12;
  c.tracingMaximumSubstepS=.125;c.tracingMaximumSubsteps=10000;c.tracingMaximumRefinements=64;
  c.tracingMaximumInteractions=1000;c.tracingMaximumSamePointInteractions=32;c.tracingAdaptive=true;
#ifdef MOON_CONTROLLER_RUNTIME_CONFIGURATION_CALLBACK
  MOON_CONTROLLER_RUNTIME_CONFIGURATION_CALLBACK(c);
#endif
  Require(Moon::Runtime::InstallConfiguration(c,&error),error);Moon::RegisterParticleTracing();Moon::ConfigureParticleTracingFields();
  PIC::ParticleBuffer::RequestDataStorage(TagOffset,sizeof(int));
  auto configuration=Configuration(),bad=configuration;bad.channels[0].scope=C::HazardScope::TrajectoryDependent;bad.exposureHistory=nullptr;
  Check(!C::Install(bad,&error) && !C::Enabled() && !PIC::Mover::QueryTrajectoryEvent,"variable_source_hazard_refused_without_activation");
  bad=configuration;bad.resolveInitialState=nullptr;
  Check(!C::Install(bad,&error) && !C::Enabled(),"missing_explicit_source_state_resolver_refused");
  Require(C::Install(configuration,&error),error);
  Check(C::Enabled() && PIC::Mover::QueryTrajectoryEvent==C::Query && PIC::Mover::AdvanceTrajectoryEvent==C::Advance &&
      PIC::Mover::ProcessTrajectoryEvent==C::Process,"complete_actual_core_hook_installation");
  double lower[3]={1.8e6,-5e5,-5e5},upper[3]={2.8e6,5e5,5e5};PIC::Mesh::mesh->AllowBlockAllocation=false;
  PIC::Mesh::mesh->init(lower,upper,Resolution);Require(PIC::Mesh::mesh->splitTreeNode(PIC::Mesh::mesh->rootTree),"root split failed");
  for(int k=0;k<8;++k)Require(PIC::Mesh::mesh->splitTreeNode(PIC::Mesh::mesh->rootTree->downNode[k]),"child split failed");
  PIC::Mesh::mesh->SetNodeNeibResolutionLevelLimit();PIC::Mesh::mesh->CreateBottomBranchNodeList(PIC::Mesh::mesh->rootTree);
  PIC::Mesh::mesh->SetParallelLoadMeasure(UnitLoad);PIC::Mesh::mesh->CreateNewParallelDistributionLists();
  PIC::Mesh::initCellSamplingDataBuffer();PIC::Mesh::mesh->AllowBlockAllocation=true;PIC::Mesh::mesh->AllocateTreeBlocks();
  PIC::Mesh::mesh->InitCellMeasure();PIC::Init_AfterParser();PIC::DomainBlockDecomposition::UpdateBlockTable();Collect(PIC::Mesh::mesh->rootTree);
  for(Node* node:Leaves)if(node->block)for(int s=0;s<PIC::nTotalSpecies;++s){node->block->SetLocalTimeStep(s==_NA_SPEC_?1:.125/(1+s),s);
    node->block->SetLocalParticleWeight(2+s+.01*(node->xmin[0]-1.8e6)/1e6,s);}
  bool found=false;
  for(Node* a:Leaves)for(Node* b:Leaves)if(!found && a->xmax[0]==b->xmin[0] && a->xmin[1]==b->xmin[1] && a->xmin[2]==b->xmin[2] &&
      (PIC::nTotalThreads==1 || a->Thread!=b->Thread)) {
    Start={{a->xmax[0]-730*.2,.5*(a->xmin[1]+a->xmax[1]),.5*(a->xmin[2]+a->xmax[2])}};found=true;
  }
  Require(found,"fixture has no actual requested ownership-crossing face");PIC::ParticleBuffer::Init(64);
  Exosphere::xObjectRadial=_AU_;
  const double sun[3]={_AU_,0,0},earth[3]={-3.8e8,2e7,5e6};
  for(int k=0;k<3;++k){Exosphere::xSun_SO[k]=sun[k];Exosphere::xEarth_SO[k]=earth[k];Exosphere::RotationVector_SO_FROZEN[k]=0;Exosphere::vObject_SO_FROZEN[k]=0;}
  T::Acceleration=nullptr;PIC::Mover::SetUserDefinedParticleMover(Moon::MoveParticle);
}
long int Source(int tag,Node* node) {
  const long int ptr=PIC::ParticleBuffer::GetNewParticle();auto* data=PIC::ParticleBuffer::GetParticleDataPointer(ptr);
  double x[3]={Start[0],Start[1]+10*tag,Start[2]},v[3]={730,-410,95};
  PIC::ParticleBuffer::SetI(_NA_SPEC_,data);PIC::ParticleBuffer::SetX(x,data);PIC::ParticleBuffer::SetV(v,data);
  PIC::ParticleBuffer::SetIndividualStatWeightCorrection(1/node->block->GetLocalParticleWeight(_NA_SPEC_),data);
  SetTag(ptr,tag);T::ResetParticleState(data);T::SetParticleEpochS(ptr,1000);return ptr;
}
void Cases() {
  Node* initial=PIC::Mesh::mesh->findTreeNode(Start.data());const int owner=initial->Thread;std::string error;
  bool recycled=true,pure=true;
  if(owner==PIC::ThisThread) {
    const long int first=Source(0,initial);C::PrepareParticle(first,1,initial);C::History old;Require(C::ReadHistory(first,&old,&error),error);
    PIC::ParticleBuffer::DeleteParticle(first);const long int second=Source(1,initial);C::PrepareParticle(second,1,initial);C::History fresh;
    Require(C::ReadHistory(second,&fresh,&error),error);recycled=first==second && Counts.source==2 && Counts.resolver==2 && fresh.clock.fresh() && fresh.clock.remaining()!=old.clock.remaining();
    PIC::ParticleBuffer::DeleteParticle(second);
  }
  Check(recycled,"recycled_native_slot_starts_exactly_one_fresh_source_clock");Counts={};
  for(bool face : {false,true}) {
    Phase=face?"boundary-face":"boundary";
    if(owner==PIC::ThisThread) {
      const long int ptr=Source(0,initial);
      // Actual configured source macro, with NO direct controller preparation.
      // This independently catches a legacy boundary path bypassing MoonMove.
      // The source has an explicitly supplied finite physical birth epoch;
      // that epoch must survive either macro and its partial .5s interval.
      if(face) {_PIC_PARTICLE_MOVER__MOVE_PARTICLE_BOUNDARY_INJECTION_FACE_(ptr,.5,initial,true,0);}
      else {_PIC_PARTICLE_MOVER__MOVE_PARTICLE_BOUNDARY_INJECTION_(ptr,.5,initial,true);}
    }
    double moving=0,exchange=0;PIC::TimeStepInternal::ExecutionTrackDefault(moving,exchange);
    double children=0;bool correct=true;
    for(unsigned index=0;index<PIC::DomainBlockDecomposition::nLocalBlocks;++index){Node* node=PIC::DomainBlockDecomposition::BlockTable[index];
      for(int cell=0;cell<_BLOCK_CELLS_X_*_BLOCK_CELLS_Y_*_BLOCK_CELLS_Z_;++cell)
        for(long int ptr=node->block->FirstCellParticleTable[cell];ptr!=-1;ptr=PIC::ParticleBuffer::GetNext(ptr)) {
          ++children;C::History h;Require(C::ReadHistory(ptr,&h,&error),error);
          correct=correct && (h.catalogSpecies==3 || h.catalogSpecies==4) && std::abs(T::GetParticleEpochS(ptr)-1000.5)<1e-10 &&
            h.interval.remainingS()<1e-12;
        }}
    const double globalChildren=Sum(children),accepted=Sum(static_cast<double>(Ledger.events.size()));
    const bool draws=Sum(Counts.source)==1 && Sum(Counts.resolver)==1 && Sum(Counts.channel)==1 && Sum(Counts.product)==2;
    Check(correct && draws && globalChildren==2 && accepted==1,
        face?"actual_face_boundary_macro_prepares_source_and_partial_children":"actual_boundary_macro_prepares_source_and_partial_children");
    PIC::ParticleBuffer::DeleteAllParticles();Ledger.events.clear();Counts={};
    CandidateEpochs.clear();CandidateRemaining.clear();CandidateTags.clear();CandidateChannels.clear();CandidateOwners.clear();
  }
  // Exact event/END/internal-AMR-face coincidence. Geometry has deterministic
  // priority; the cached due event must nevertheless publish at this SAME
  // endpoint on a legitimate current owner, with zero-time children. This
  // control does not prescribe which side owns a child on the exact plane.
  // Choosing an exact1m path at1m/s under an explicit zero force avoids tuning
  // an initial position against numerical lunar-force trajectory residuals.
  Phase="endpoint-face";T::Acceleration=DeclaredZeroAcceleration;
  if(owner==PIC::ThisThread) {
    const long int ptr=Source(2,initial);double x[3]={initial->xmax[0]-1,Start[1]+20,Start[2]},v[3]={1,0,0};
    PIC::ParticleBuffer::SetX(x,ptr);PIC::ParticleBuffer::SetV(v,ptr);C::PrepareParticle(ptr,1,initial);
    int i,j,k;PIC::Mesh::mesh->FindCellIndex(x,i,j,k,initial,false);
    long int& head=initial->block->FirstCellParticleTable[i+_BLOCK_CELLS_X_*(j+_BLOCK_CELLS_Y_*k)];
    Moon::Chemistry::ProcessParticle(ptr,head,initial);
  }
  double tieMoving=0,tieExchange=0;PIC::TimeStepInternal::ExecutionTrackDefault(tieMoving,tieExchange);
  double tieChildren=0;bool tieState=true;
  for(unsigned index=0;index<PIC::DomainBlockDecomposition::nLocalBlocks;++index){Node* node=PIC::DomainBlockDecomposition::BlockTable[index];
    for(int cell=0;cell<_BLOCK_CELLS_X_*_BLOCK_CELLS_Y_*_BLOCK_CELLS_Z_;++cell)
      for(long int ptr=node->block->FirstCellParticleTable[cell];ptr!=-1;ptr=PIC::ParticleBuffer::GetNext(ptr)) {
        ++tieChildren;C::History h;Require(C::ReadHistory(ptr,&h,&error),error);
        tieState=tieState && (h.catalogSpecies==3 || h.catalogSpecies==4) && !h.due &&
          h.interval.remainingS()==0 && T::GetParticleEpochS(ptr)==1001;
      }}
  const double tieCount=Sum(tieChildren),tieEvents=Sum(static_cast<double>(Ledger.events.size())),
    tieSources=Sum(Counts.source),tieChannels=Sum(Counts.channel),tieProducts=Sum(Counts.product);
  bool tieCandidates=true;for(std::size_t k=0;k<CandidateEpochs.size();++k)
    tieCandidates=tieCandidates && CandidateEpochs[k]==1001 && CandidateRemaining[k]==0;
  Check(tieState && tieCandidates && tieCount==2 && tieEvents==1 && tieSources==1 && tieChannels==1 && tieProducts==2,
      "exact_internal_face_endpoint_event_commits_zero_time_children");
  PIC::ParticleBuffer::DeleteAllParticles();Ledger.events.clear();Counts={};T::Acceleration=nullptr;
  CandidateEpochs.clear();CandidateRemaining.clear();CandidateTags.clear();CandidateChannels.clear();CandidateOwners.clear();
  // The immediately future binary64 optical deadline MUST NOT be rounded
  // into the endpoint by a clock/acceptance tolerance. Continue that SAME
  // surviving pending clock into the next outer interval; its tiny positive
  // residual also detects cancellation in delay-duration+remaining arithmetic.
  Phase="future";FutureClock=true;
  if(owner==PIC::ThisThread) {
    const long int ptr=Source(2,initial);C::PrepareParticle(ptr,1,initial);
    int i,j,k;double x[3];PIC::ParticleBuffer::GetX(x,ptr);PIC::Mesh::mesh->FindCellIndex(x,i,j,k,initial,false);
    long int& head=initial->block->FirstCellParticleTable[i+_BLOCK_CELLS_X_*(j+_BLOCK_CELLS_Y_*k)];
    Moon::Chemistry::ProcessParticle(ptr,head,initial);
  }
  double futureMoving=0,futureExchange=0;PIC::TimeStepInternal::ExecutionTrackDefault(futureMoving,futureExchange);
  double survivors=0;bool futureState=true;
  for(unsigned index=0;index<PIC::DomainBlockDecomposition::nLocalBlocks;++index){Node* node=PIC::DomainBlockDecomposition::BlockTable[index];
    for(int cell=0;cell<_BLOCK_CELLS_X_*_BLOCK_CELLS_Y_*_BLOCK_CELLS_Z_;++cell)
      for(long int ptr=node->block->FirstCellParticleTable[cell];ptr!=-1;ptr=PIC::ParticleBuffer::GetNext(ptr)) {
        ++survivors;C::History h;Require(C::ReadHistory(ptr,&h,&error),error);
        futureState=futureState && h.catalogSpecies==2 && !h.due && h.clock.remaining()>0 &&
          h.interval.remainingS()==0 && T::GetParticleEpochS(ptr)==1001;
      }}
  const double futureCount=Sum(survivors),futureEvents=Sum(static_cast<double>(Ledger.events.size())),
    futureSources=Sum(Counts.source),futureChannels=Sum(Counts.channel);
  Check(futureState && futureCount==1 && futureEvents==0 && futureSources==1 && futureChannels==0,
      "strict_nextafter_future_deadline_does_not_fire_at_endpoint");
  PIC::TimeStepInternal::ExecutionTrackDefault(futureMoving,futureExchange);
  double continuedChildren=0;bool continuedState=true;
  for(unsigned index=0;index<PIC::DomainBlockDecomposition::nLocalBlocks;++index){Node* node=PIC::DomainBlockDecomposition::BlockTable[index];
    for(int cell=0;cell<_BLOCK_CELLS_X_*_BLOCK_CELLS_Y_*_BLOCK_CELLS_Z_;++cell)
      for(long int ptr=node->block->FirstCellParticleTable[cell];ptr!=-1;ptr=PIC::ParticleBuffer::GetNext(ptr)) {
        ++continuedChildren;C::History h;Require(C::ReadHistory(ptr,&h,&error),error);
        continuedState=continuedState && (h.catalogSpecies==3 || h.catalogSpecies==4) && h.interval.remainingS()==0 &&
          T::GetParticleEpochS(ptr)==1002;
      }}
  const double continuedCount=Sum(continuedChildren),continuedEvents=Sum(static_cast<double>(Ledger.events.size())),
    continuedSources=Sum(Counts.source),continuedChannels=Sum(Counts.channel),continuedProducts=Sum(Counts.product);
  Check(continuedState && continuedCount==2 && continuedEvents==1 && continuedSources==1 &&
      continuedChannels==1 && continuedProducts==2,
      "surviving_tiny_clock_continues_next_outer_interval_without_source_redraw");
  PIC::ParticleBuffer::DeleteAllParticles();Ledger.events.clear();Counts={};FutureClock=false;
  CandidateEpochs.clear();CandidateRemaining.clear();CandidateTags.clear();CandidateChannels.clear();CandidateOwners.clear();
  Phase="main";
  for(int tag=0;tag<3;++tag)if(owner==PIC::ThisThread) {
    const long int ptr=Source(tag,initial);C::PrepareParticle(ptr,1,initial);C::History before,after;Require(C::ReadHistory(ptr,&before,&error),error);
    double event=0;for(int repeat=0;repeat<100;++repeat)C::Query(ptr,_NA_SPEC_,initial,1,.1,event);
    Require(C::ReadHistory(ptr,&after,&error),error);pure=pure && before.clock.remaining()==after.clock.remaining() && before.clock.fresh()==after.clock.fresh() && Counts.source==tag+1 && Counts.channel==0;
    int i,j,k;double x[3];PIC::ParticleBuffer::GetX(x,ptr);PIC::Mesh::mesh->FindCellIndex(x,i,j,k,initial,false);
    long int& head=initial->block->FirstCellParticleTable[i+_BLOCK_CELLS_X_*(j+_BLOCK_CELLS_Y_*k)];
    Moon::Chemistry::ProcessParticle(ptr,head,initial); // actual legacy bypass
  }
  Check(pure,"repeated_refinement_previews_change_no_clock_or_rng");
  Check(Sum(Counts.source)==3 && Sum(Counts.resolver)==3,"explicit_state_and_source_clock_once_per_birth");
  if(Capacity) {
    Require(PIC::nTotalThreads==1,"resource-negative fixture uses one process");
    while(PIC::ParticleBuffer::GetMaxNPart()-PIC::ParticleBuffer::GetAllPartNum()>1)PIC::ParticleBuffer::GetNewParticle();
  }
  double moving=0,exchange=0;PIC::TimeStepInternal::ExecutionTrackDefault(moving,exchange);
  if(Capacity)throw std::runtime_error("resource-negative controller failed to refuse an incomplete packet");
  Check(Sum(Counts.source)==3 && Sum(Counts.resolver)==3 && Sum(Counts.channel)==3 && Sum(Counts.product)==6,
      "combined_hazard_and_product_draw_counts_no_reseed_after_mpi");
  bool events=CandidateEpochs.size()==Ledger.events.size(),budget=true;double parentN=0,productN=0;
  for(std::size_t i=0;i<CandidateEpochs.size();++i) {
    const int tag=CandidateTags[i],branch=tag==1?1:0;const double delay=tag==0?.4:tag==1?.6:1;
    events=events && CandidateChannels[i]==2*tag+branch && std::abs(CandidateEpochs[i]-(1000+delay))<1e-10 &&
      std::abs(CandidateRemaining[i]-(1-delay))<1e-12 && Ledger.events[i].identity=="CONTROL-source"+std::to_string(tag)+"-branch"+std::to_string(branch);
    const auto& entry=Ledger.events[i].ledger;budget=budget && entry.integerConservationVerified && entry.kinematicConservationVerified;
    parentN+=entry.parentNumber;productN+=entry.productNumber;
  }
  const double accepted=Sum(static_cast<double>(Ledger.events.size())),globalParentN=Sum(parentN),globalProductN=Sum(productN);
  Check(events && accepted==3,"selected_source_actual_event_epoch_and_remaining_interval");
  Check(budget && std::abs(globalParentN-3)<1e-12 && std::abs(globalProductN-6)<1e-12,"accepted_exact_material_and_represented_number_budget");
  double count=0,weight=0,ions=0,helium=0;bool state=true;
  for(unsigned index=0;index<PIC::DomainBlockDecomposition::nLocalBlocks;++index){Node* node=PIC::DomainBlockDecomposition::BlockTable[index];
    for(int cell=0;cell<_BLOCK_CELLS_X_*_BLOCK_CELLS_Y_*_BLOCK_CELLS_Z_;++cell)for(long int ptr=node->block->FirstCellParticleTable[cell];ptr!=-1;ptr=PIC::ParticleBuffer::GetNext(ptr)) {
      ++count;const int species=PIC::ParticleBuffer::GetI(ptr);ions+=species==_NA_PLUS_SPEC_;helium+=species==_HE_SPEC_;
      C::History h;Require(C::ReadHistory(ptr,&h,&error),error);const auto time=T::GetTimeContext(ptr);
      const double n=node->block->GetLocalParticleWeight(species)*PIC::ParticleBuffer::GetIndividualStatWeightCorrection(ptr);weight+=n;
      state=state && (h.catalogSpecies==3 || h.catalogSpecies==4) && !h.due && h.clock.valid() && !h.clock.consumed() &&
        std::abs(h.interval.remainingS())<1e-12 && std::abs(T::GetParticleEpochS(ptr)-1001)<1e-10 && time.RemainingTimeS==0 && std::abs(n-1)<1e-13;
      ParticleFile<<Tag(ptr)<<','<<species<<','<<h.catalogSpecies<<','<<n<<','<<T::GetParticleEpochS(ptr)<<','<<h.interval.remainingS()<<','<<h.clock.remaining()<<','<<PIC::ThisThread<<'\n';
    }}
  const double globalCount=Sum(count),globalIons=Sum(ions),globalHelium=Sum(helium),globalWeight=Sum(weight);
  Check(state && globalCount==6 && globalIons==3 && globalHelium==3 && std::abs(globalWeight-6)<1e-12,
      "all_daughters_complete_original_physical_interval_without_birth_dt_factor");
  double crossed=0;for(int eventOwner:CandidateOwners)crossed+=eventOwner!=owner;
  Check(PIC::nTotalThreads==1 || Sum(crossed)>=1,"actual_pending_clock_crosses_mpi_before_event");
  Check(T::CountPending()==0,"no_unfinished_child_or_parent_continuation");
  PIC::ParticleBuffer::DeleteAllParticles();Check(PIC::ParticleBuffer::GetAllPartNum()==0,"no_parent_double_delete_or_native_slot_leak");
}
}
int main(int argc,char** argv) {
  PIC::InitMPI();SetExitErrorCode(__LINE__,1);if(argc!=4){MPI_Finalize();return 2;}
  Order=std::stoi(argv[1]);Capacity=std::string(argv[2])=="capacity";
  const std::string prefix=std::string(argv[3])+".rank"+std::to_string(PIC::ThisThread);
  if(PIC::ThisThread==0){CheckFile.open(std::string(argv[3])+".checks.csv");CheckFile<<"check,status\n";}
  EventFile.open(prefix+".candidates.csv");EventFile<<std::setprecision(17)<<"case,source_tag,channel,event_epoch_s,remaining_s,owner\n";
  ParticleFile.open(prefix+".particles.csv");ParticleFile<<std::setprecision(17)<<"source_tag,compiled_species,catalog_state,represented_number,epoch_s,remaining_s,optical_depth,owner\n";
  int error=0;try{Setup();Cases();}catch(const std::exception& exception){error=1;std::cerr<<"ERROR native controller: "<<exception.what()<<'\n';}
  int bad=0,failed=0;MPI_Allreduce(&error,&bad,1,MPI_INT,MPI_MAX,MPI_GLOBAL_COMMUNICATOR);MPI_Allreduce(&Failures,&failed,1,MPI_INT,MPI_MAX,MPI_GLOBAL_COMMUNICATOR);
  EventFile.flush();ParticleFile.flush();if(PIC::ThisThread==0){CheckFile.flush();std::cout<<"M3 native reaction controller order="<<Order<<" checks="<<Checks<<" status="<<(bad?"ERROR":failed?"FAIL":"PASS")<<'\n';}
  MPI_Finalize();return bad?2:failed?1:0;
}
