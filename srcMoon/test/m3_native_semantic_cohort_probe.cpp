// Analytical normal-startup/PIC composition controls. Nominal heavy masses
// and the declared zero-Q topology are NOT physical lunar ionization/CEX.
#include "pic.h"
#include "Moon.h"
#include "MoonInput.h"
#include "MoonSemanticRandom.h"
#include "MoonSurfaceEnvironment.h"
#include <fstream>
#include <limits>

namespace SR=Moon::Chemistry::SemanticRandom;
namespace PC=Moon::Chemistry::Controller;
namespace PR=Moon::Species::ReactionExecution;
namespace QR=Moon::Species::QualificationRandom;
namespace {
constexpr double ExactSICharge=1.602176634e-19;
double SelectedCharge=ExactSICharge;
bool FailStaging=false;
std::ofstream KeyFile;
SR::Configuration RandomConfiguration(){SR::Configuration c;c.seed=23;
  c.cohort="DECLARED-SEMANTIC-PIC-COMPOSITION-NOT-LUNAR-PHYSICS";
  c.genealogyCapacity=4;return c;}
bool Phase2Kinematics(const PC::EventContext&,PR::KinematicPolicy*,std::string*);
void Phase2Configuration(PC::Configuration&);
void Phase2RuntimeConfiguration(Moon::Runtime::Configuration&);
}

// Reuse the existing *test infrastructure*: actual AMR setup, MPI ownership,
// ordinary Moon forces, source initialization and named collective checks.
// The included main/cases are not executed. No force or reaction production
// service is copied into this fixture.
#define MOON_CONTROLLER_CONFIGURATION_CALLBACK Phase2Configuration
#define MOON_CONTROLLER_RUNTIME_CONFIGURATION_CALLBACK Phase2RuntimeConfiguration
#define main UnusedOriginalControllerMain
#include "m3_native_reaction_controller_probe.cpp"
#undef main
#undef MOON_CONTROLLER_CONFIGURATION_CALLBACK
#undef MOON_CONTROLLER_RUNTIME_CONFIGURATION_CALLBACK

void amps_init();
namespace {
void Phase2RuntimeConfiguration(Moon::Runtime::Configuration& c){
  // Runtime installation is immutable. Declare the optional artificial
  // collector at the ONE initial setup, then install its explicit catalog
  // later before the collector control. No earlier packet hits the surface.
  c.ionSurfaceBoundary="ideal-collector";
}
bool Phase2Initializer(const PR::EventPlan& plan,std::size_t ordinal,long int parent,
    long int child,void*,std::string* error){
  if(!SR::InitializeStagedProduct(plan,ordinal,parent,child,nullptr,error))return false;
  // A fallible staged application initializer is deliberately refused AFTER
  // actual child bytes exist. It must never publish a partial physical packet.
  if(FailStaging && ordinal==1){*error="declared staged-initializer refusal";return false;}
  return true;
}
void Phase2Configuration(PC::Configuration& c){
  std::string error;
  if(!SR::Registered())Require(SR::RegisterStorage(RandomConfiguration(),&error),error);
  c.identity="ANALYTICAL-SEMANTIC-PIC-COMPOSITION-v1";
  c.elementaryChargeC=SelectedCharge;
  c.elementaryChargeSourceIdentity=SelectedCharge==ExactSICharge?
    "SI-2019-exact-elementary-charge-ANALYTICAL-storage-fixture":"AMPS legacy ElectronCharge";
  c.openUniform=SR::DrawForParticle;
  c.initializeProduct=Phase2Initializer;
  c.productInitializerIdentity="full-genealogy-staged-software-control";
  c.channels.resize(1);c.channels[0].constantRateS=1;
  c.channels[0].kinematics=Phase2Kinematics;
  c.channels[0].kinematicProviderIdentity="ANALYTICAL-semantic-zero-Q-policy";
}
bool Phase2Kinematics(const PC::EventContext& event,PR::KinematicPolicy* policy,std::string* error){
  policy->qualified=true;policy->identity="ANALYTICAL-semantic-zero-Q-policy";
  policy->evidenceIdentity="STORAGE-AND-RNG-CONTROL-NOT-PHYSICAL-CHANNEL";
  policy->eventFrame="LSO";policy->maximumClassicalSpeedMS=1e5;
  policy->bounds.modelMassKg={0,1e-13};policy->bounds.momentumKgMS={0,1e-13};
  policy->bounds.energyJ={0,1e-13};policy->bath.supplied=true;
  const double bathMass=event.catalog->at(5).modelMassKg;
  policy->bath.modelMassChangeKg=-bathMass;
  // Explicit analytical incoming bath velocity equals parent velocity. Both
  // outgoing nominal particles have that same velocity: no empirical recoil,
  // zero-recoil physical assumption, or fitted bath residual is introduced.
  for(int k=0;k<3;++k){policy->bath.momentumChangeKgMS[k]=-bathMass*event.parent.velocityMS[k];
    policy->bath.energyChangeJ-=.5*bathMass*event.parent.velocityMS[k]*event.parent.velocityMS[k];}
  // Exercise purpose separation through actual opaque PIC counters. These
  // uniforms do not define an unapproved physical angular distribution.
  double unused=0;
  for(int draw=0;draw<2;++draw)
    if(!SR::DrawForParticle(PC::DrawPurpose::Kinematics,event.particle,event.channelIndex,nullptr,&unused,error))return false;
  std::vector<PR::ProductKinematics> products;
  for(std::size_t ordinal=0;ordinal<2;++ordinal){PR::ProductKinematics p;p.term=ordinal;
    p.referenceWeight=event.node->block->GetLocalParticleWeight(ordinal?_HE_SPEC_:_NA_PLUS_SPEC_);
    for(int k=0;k<3;++k)p.velocityMS[k]=event.parent.velocityMS[k];
    double u=0;if(!SR::ProspectiveProductClockForParticle(event.particle,ordinal,&u,error) ||
        !PR::OpticalDepthClock::FromOpenUniform(u,&p.nextClock,error))return false;
    products.push_back(p);}
  policy->products=std::move(products);return true;
}
std::vector<unsigned char> Bytes(long int ptr){const auto* p=PIC::ParticleBuffer::GetParticleDataPointer(ptr);
  return std::vector<unsigned char>(p,p+PIC::ParticleBuffer::ParticleDataLength);}
bool SameBytes(long int ptr,const std::vector<unsigned char>& bytes){return PIC::ParticleBuffer::IsParticleAllocated(ptr) &&
  std::memcmp(bytes.data(),PIC::ParticleBuffer::GetParticleDataPointer(ptr),bytes.size())==0;}
bool SameKey(const SR::State& a,const SR::State& b){return a.originalParent==b.originalParent &&
  a.genealogy==b.genealogy && a.acceptedDraws==b.acceptedDraws;}
Node* OwnedNode(){for(Node* node:Leaves)if(node->Thread==PIC::ThisThread && node->block)return node;
  throw std::runtime_error("no actual owned PIC leaf");}
long int NewSemanticSource(Node* node,std::uint64_t id){const long int ptr=Source(0,node);
  double x[3];for(int k=0;k<3;++k)x[k]=.5*(node->xmin[k]+node->xmax[k]);
  PIC::ParticleBuffer::SetX(x,ptr);std::string error;Require(SR::InitializeSource(ptr,id,&error),error);return ptr;}
PR::EventPlan DirectPlan(long int ptr,Node* node,const PC::Configuration& c){
  PC::EventContext event;event.particle=ptr;event.node=node;event.catalog=&c.catalog;event.reaction=&c.channels[0].reaction;
  event.parent.species=0;event.parent.referenceWeight=node->block->GetLocalParticleWeight(_NA_SPEC_);
  event.parent.individualCorrection=PIC::ParticleBuffer::GetIndividualStatWeightCorrection(ptr);
  double x[3],v[3];PIC::ParticleBuffer::GetX(x,ptr);PIC::ParticleBuffer::GetV(v,ptr);
  for(int k=0;k<3;++k){event.parent.positionM[k]=x[k];event.parent.velocityMS[k]=v[k];}
  std::string error;Require(PR::TraceInterval::Create(1000,10,3,&event.interval,&error),error);
  PIC::Mover::ForceCutSurface::SetParticleEpochS(ptr,1003);
  PR::KinematicPolicy policy;Require(Phase2Kinematics(event,&policy,&error),error);
  PR::EventPlan plan;Require(PR::MakeEventPlan(c.catalog,c.channels[0].reaction,event.parent,event.interval,
      PR::PlanScope::Transport,&policy,&plan,&error),error);return plan;
}
void DirectCases(){
  Node* node=OwnedNode();auto configuration=Configuration();std::string error;
  Check(PIC::MolecularData::GetElectricCharge(_NA_PLUS_SPEC_)==SelectedCharge &&
      PIC::MolecularData::GetElectricCharge(_HE_PLUS_SPEC_)==SelectedCharge &&
      PIC::MolecularData::GetElectricCharge(_NE_PLUS_SPEC_)==SelectedCharge &&
      PIC::MolecularData::GetElectricCharge(_AR_PLUS_SPEC_)==SelectedCharge &&
      PIC::MolecularData::GetElectricCharge(_NA_SPEC_)==0 && PIC::MolecularData::GetElectricCharge(_HE_SPEC_)==0 &&
      PIC::MolecularData::GetElectricCharge(_NE_SPEC_)==0 && PIC::MolecularData::GetElectricCharge(_AR_SPEC_)==0,
      "actual_compiled_ions_match_explicit_charge_unit");
  const long int parent=NewSemanticSource(node,73);const auto plan=DirectPlan(parent,node,configuration);
  const auto snapshot=Bytes(parent);SR::State original;Require(SR::Read(parent,&original,&error),error);
  long int head=-1;const auto before=PIC::ParticleBuffer::GetAllPartNum();
  {
    Moon::Chemistry::ParticleTransaction wrong(parent,head,node,configuration.catalog,&Ledger,
      Moon::Chemistry::ProductPlacement::DetachedContinuation,true,
      SelectedCharge==ExactSICharge?ElectronCharge:ExactSICharge,Phase2Initializer);
    Check(wrong.Reserve(plan,&error)==PR::ReservationStatus::Error && SameBytes(parent,snapshot) && head==-1 && Ledger.events.empty(),
        "wrong_charge_unit_refused_before_reservation");
  }
  {
    Moon::Chemistry::ParticleTransaction transaction(parent,head,node,configuration.catalog,&Ledger,
      Moon::Chemistry::ProductPlacement::DetachedContinuation,true,SelectedCharge,Phase2Initializer);
    Require(transaction.Reserve(plan,&error)==PR::ReservationStatus::Ready,error);
    FailStaging=true;const bool refused=!transaction.Stage(plan,&error);FailStaging=false;
    transaction.Rollback();SR::State after;Require(SR::Read(parent,&after,&error),error);
    Check(refused && SameKey(original,after) && SameBytes(parent,snapshot) && head==-1 && Ledger.events.empty() &&
        PIC::ParticleBuffer::GetAllPartNum()==before,"staged_partial_child_failure_rolls_back_complete_packet");
  }
  std::vector<long int> children;
  {
    Moon::Chemistry::ParticleTransaction transaction(parent,head,node,configuration.catalog,&Ledger,
      Moon::Chemistry::ProductPlacement::DetachedContinuation,true,SelectedCharge,Phase2Initializer);
    Require(transaction.Reserve(plan,&error)==PR::ReservationStatus::Ready,error);
    Require(transaction.Stage(plan,&error),error);Require(transaction.Commit(plan,&error),error);
    children=transaction.committedProducts();
  }
  SR::State after;Require(SR::Read(parent,&after,&error),error);
  Check(children.size()==2 && Ledger.events.size()==1 && SameKey(original,after) && SameBytes(parent,snapshot),
      "same_plan_retry_commits_without_parent_rng_redraw");
  bool keys=true,clocks=true,intervals=true,clone=true;
  for(std::size_t ordinal=0;ordinal<children.size();++ordinal){SR::State actual,expected;
    Require(SR::Read(children[ordinal],&actual,&error),error);
    Require(SR::ChildState(RandomConfiguration(),original,ordinal,true,&expected,&error),error);
    keys=keys && SameKey(actual,expected) && actual.acceptedDraws[2]==1 && actual.acceptedDraws[0]==0 && actual.acceptedDraws[3]==0;
    PR::TraceInterval interval;PR::OpticalDepthClock clock;std::size_t state;double energy;
    Require(Moon::Chemistry::ReadReactionState(children[ordinal],&interval,&clock,&state,&energy,&error),error);
    clocks=clocks && clock.fresh() && clock.remainingExtended()==plan.products()[ordinal].clock.remainingExtended();
    const int compiled=PIC::ParticleBuffer::GetI(children[ordinal]);
    const double number=node->block->GetLocalParticleWeight(compiled)*PIC::ParticleBuffer::GetIndividualStatWeightCorrection(children[ordinal]);
    const double expectedNumber=plan.parent().referenceWeight*plan.parent().individualCorrection;
    // One division stores N/W_child and one multiplication reconstructs N.
    // Two binary64 roundings bound this diagnostic; it is not a tunable
    // represented-population or physical conservation acceptance tolerance.
    const double reconstructionBound=2*std::numeric_limits<double>::epsilon()*std::abs(expectedNumber);
    intervals=intervals && interval.startS()==1000 && interval.durationS()==10 && interval.elapsedS()==3 && interval.remainingS()==7 &&
      std::abs(number-expectedNumber)<=reconstructionBound;
    const long int copy=PIC::ParticleBuffer::GetNewParticle();
    PIC::ParticleBuffer::CloneParticle(PIC::ParticleBuffer::GetParticleDataPointer(copy),PIC::ParticleBuffer::GetParticleDataPointer(children[ordinal]));
    SR::State copied;Require(SR::Read(copy,&copied,&error),error);clone=clone && SameKey(actual,copied);
    PIC::ParticleBuffer::DeleteParticle_withoutTrajectoryTermination(copy);
    double u=0;Require(SR::ProspectiveProductClock(RandomConfiguration(),original,ordinal,&u,&error),error);
    KeyFile<<"direct,"<<original.originalParent<<','<<ordinal<<','<<u<<','<<clock.remaining()<<','<<actual.acceptedDraws[2]<<','<<PIC::ThisThread<<'\n';
  }
  Check(keys,"staged_daughters_use_complete_child_genealogy_and_clock_counter_one");
  Check(clocks,"staged_clocks_equal_original_prospective_child_draws");
  Check(intervals,"staged_products_keep_represented_number_and_original_local_interval");
  Check(clone,"real_pic_ownership_clone_preserves_semantic_bytes");
  for(long int ptr:children)PIC::ParticleBuffer::DeleteParticle_withoutTrajectoryTermination(ptr);
  PIC::ParticleBuffer::DeleteParticle_withoutTrajectoryTermination(parent);Ledger.events.clear();
  const long int recycled=NewSemanticSource(node,91);SR::State fresh;Require(SR::Read(recycled,&fresh,&error),error);
  Check(recycled==parent && fresh.originalParent==91 && fresh.genealogy.empty() && fresh.acceptedDraws==std::array<std::uint64_t,4>{{0,0,0,0}},
      "recycled_source_slot_explicitly_starts_fresh_semantic_history");
  PIC::ParticleBuffer::DeleteParticle_withoutTrajectoryTermination(recycled);
  Check(PIC::ParticleBuffer::GetAllPartNum()==0 && head==-1,"transaction_controls_leave_no_native_slots_or_published_list");
}
void CapacityCase(){
  // This is an identity-topology storage chain, NOT an asserted physical
  // repeated sodium process. It reaches the EXACT configured genealogy
  // capacity through real staged PIC records, then requires pre-publication
  // refusal rather than truncating, hash-folding or censoring a daughter.
  Node* node=OwnedNode();auto c=Configuration();std::string error;
  long int parent=NewSemanticSource(node,97);PIC::ParticleBuffer::SetI(_NA_PLUS_SPEC_,parent);
  PIC::ParticleBuffer::SetIndividualStatWeightCorrection(1/node->block->GetLocalParticleWeight(_NA_PLUS_SPEC_),parent);
  const PR::Reaction reaction{"ANALYTICAL-IDENTITY-STORAGE-NOT-REACTION",{{3,1,PR::Population::Tracked}},{{3,1,PR::Population::Tracked}}};
  bool valid=true;
  for(std::size_t generation=0;generation<RandomConfiguration().genealogyCapacity;++generation){
    PR::ParentState p;p.species=3;p.referenceWeight=node->block->GetLocalParticleWeight(_NA_PLUS_SPEC_);
    p.individualCorrection=PIC::ParticleBuffer::GetIndividualStatWeightCorrection(parent);
    double x[3],v[3];PIC::ParticleBuffer::GetX(x,parent);PIC::ParticleBuffer::GetV(v,parent);
    for(int k=0;k<3;++k){p.positionM[k]=x[k];p.velocityMS[k]=v[k];}
    PR::TraceInterval interval;Require(PR::TraceInterval::Create(1000,10,3,&interval,&error),error);
    PIC::Mover::ForceCutSurface::SetParticleEpochS(parent,1003);
    PR::KinematicPolicy policy;policy.qualified=true;policy.identity="ANALYTICAL-identity-storage";
    policy.evidenceIdentity="GENEALOGY-CAPACITY-CONTROL-NOT-PHYSICS";policy.eventFrame="LSO";
    policy.maximumClassicalSpeedMS=1e5;policy.bath.supplied=true;
    PR::ProductKinematics product;product.term=0;product.referenceWeight=p.referenceWeight;
    for(int k=0;k<3;++k)product.velocityMS[k]=v[k];
    double u=0;Require(SR::ProspectiveProductClockForParticle(parent,0,&u,&error),error);
    Require(PR::OpticalDepthClock::FromOpenUniform(u,&product.nextClock,&error),error);policy.products={product};
    PR::EventPlan plan;Require(PR::MakeEventPlan(c.catalog,reaction,p,interval,PR::PlanScope::Transport,&policy,&plan,&error),error);
    long int head=-1;std::vector<long int> children;
    {Moon::Chemistry::ParticleTransaction transaction(parent,head,node,c.catalog,&Ledger,
      Moon::Chemistry::ProductPlacement::DetachedContinuation,true,SelectedCharge,Phase2Initializer);
      Require(transaction.Reserve(plan,&error)==PR::ReservationStatus::Ready,error);
      Require(transaction.Stage(plan,&error),error);Require(transaction.Commit(plan,&error),error);children=transaction.committedProducts();}
    Require(children.size()==1,"identity control must publish one complete child");
    PIC::ParticleBuffer::DeleteParticle_withoutTrajectoryTermination(parent);parent=children[0];
    SR::State state;Require(SR::Read(parent,&state,&error),error);
    valid=valid && state.originalParent==97 && state.genealogy.size()==generation+1 && state.acceptedDraws[2]==1;
  }
  Check(valid,"real_staged_identity_chain_retains_full_genealogy_to_capacity");
  const auto before=Bytes(parent);const auto events=Ledger.events.size();const auto count=PIC::ParticleBuffer::GetAllPartNum();
  double unchanged=-9;const bool refused=!SR::ProspectiveProductClockForParticle(parent,0,&unchanged,&error);
  Check(refused && unchanged==-9 && !error.empty() && SameBytes(parent,before) && Ledger.events.size()==events &&
      PIC::ParticleBuffer::GetAllPartNum()==count,"actual_genealogy_capacity_refuses_before_clock_or_particle_publication");
  PIC::ParticleBuffer::DeleteParticle_withoutTrajectoryTermination(parent);Ledger.events.clear();
}
void MoverCase(){
  std::string error;Node* initial=PIC::Mesh::mesh->findTreeNode(Start.data());
  if(initial->Thread==PIC::ThisThread){const long int ptr=Source(0,initial);Require(SR::InitializeSource(ptr,73,&error),error);
    Moon::MoveParticle(ptr,64,initial);}
  double moving=0,exchange=0;PIC::TimeStepInternal::ExecutionTrackDefault(moving,exchange);
  bool valid=true;double count=0,migrated=0;
  for(unsigned i=0;i<PIC::DomainBlockDecomposition::nLocalBlocks;++i){Node* node=PIC::DomainBlockDecomposition::BlockTable[i];
    for(int cell=0;cell<_BLOCK_CELLS_X_*_BLOCK_CELLS_Y_*_BLOCK_CELLS_Z_;++cell)
      for(long int ptr=node->block->FirstCellParticleTable[cell];ptr!=-1;ptr=PIC::ParticleBuffer::GetNext(ptr)){
        ++count;SR::State state;Require(SR::Read(ptr,&state,&error),error);PC::History h;Require(PC::ReadHistory(ptr,&h,&error),error);
        valid=valid && state.originalParent==73 && state.genealogy.size()==1 && state.genealogy[0]<2 &&
          state.acceptedDraws[0]==0 && state.acceptedDraws[1]==0 && state.acceptedDraws[2]==1 && state.acceptedDraws[3]==0 &&
          h.interval.remainingS()==0 && T::GetParticleEpochS(ptr)==1064;
        migrated+=node->Thread!=initial->Thread;
        KeyFile<<"mover,"<<state.originalParent<<','<<state.genealogy[0]<<",0,"<<h.clock.remaining()<<','<<state.acceptedDraws[2]<<','<<PIC::ThisThread<<'\n';
      }}
  const double globalCount=Sum(count),events=Sum(static_cast<double>(Ledger.events.size()));
  Check(valid && globalCount==2 && events==1 && T::CountPending()==0,"actual_movers_continue_complete_semantic_daughter_packet");
  Check(PIC::nTotalThreads==1 || Sum(migrated)>=1,"actual_mpi_migration_preserves_full_daughter_genealogy");
  PIC::ParticleBuffer::DeleteAllParticles();Ledger.events.clear();
  Check(PIC::ParticleBuffer::GetAllPartNum()==0,"semantic_mover_parent_is_deleted_once_without_slot_leak");
}
bool SICollectorResolver(long int ptr,void*,void*,std::size_t* state,std::string* error){
  if(PIC::ParticleBuffer::GetI(ptr)!=_NA_PLUS_SPEC_){*error="analytical SI collector has one explicitly declared nominal ion";return false;}
  *state=0;return true;
}
void SICollectorCase(){
  namespace I=Moon::Species::IonCollector;
  static I::Ledger collector;std::string error;
  Require(Moon::Runtime::GetConfiguration().ionSurfaceBoundary=="ideal-collector",
    "analytical collector must be selected before the one immutable runtime installation");
  I::StateDefinition state;state.material={"Na+","ANALYTICAL-nominal-SI-storage",_NA_PLUS_SPEC_,{{"Na",1}},1,PIC::MolecularData::GetMass(_NA_PLUS_SPEC_)};
  state.formationEnergyJ=0;state.excitationEnergyJ=0;
  state.energyReferenceIdentity="DECLARED-zero-Q-numerical-storage-reference";
  state.sourceEvidenceIdentity="ANALYTICAL-fixture-NOT-lunar-ion-state-data";
  Check(!Moon::SurfaceEnvironment::InstallProductionIonCollector({state},SICollectorResolver,nullptr,&collector,&error,
      SelectedCharge==ExactSICharge?ElectronCharge:ExactSICharge,"deliberately-incompatible-charge-source") &&
      !Moon::SurfaceEnvironment::ProductionIonCollectorInstalled(),"collector_rejects_mismatched_compiled_charge_unit_before_install");
  Check(Moon::SurfaceEnvironment::InstallProductionIonCollector({state},SICollectorResolver,nullptr,&collector,&error,
      SelectedCharge,"SI-2019-exact-elementary-charge-ANALYTICAL-storage-fixture"),
      "collector_accepts_same_explicit_source_charge_as_controller_and_transaction");
  for(int i=0;i<6;++i)for(int j=0;j<6;++j){
    Moon::OrbitalMotion::SO_to_IAU_TransformationMartix[i][j]=i==j?1:0;
    Moon::OrbitalMotion::IAU_to_SO_TransformationMartix[i][j]=i==j?1:0;}
  Moon::OrbitalMotion::et=2000;Node* node=OwnedNode();const long int ptr=NewSemanticSource(node,109);
  PIC::ParticleBuffer::SetI(_NA_PLUS_SPEC_,ptr);T::SetParticleEpochS(ptr,1003);
  double x[3],v[3];PIC::ParticleBuffer::GetX(x,ptr);const double radius=std::hypot(x[0],x[1],x[2]);
  for(int k=0;k<3;++k)v[k]=-2*x[k]/radius;
  PIC::ParticleBuffer::SetV(v,ptr);const auto before=Bytes(ptr);
  cInternalSphericalData sphere;double center[3]={0,0,0};sphere.SetSphereGeometricalParameters(center,radius);
  double remaining=7;const int result=Moon::IonAwareSphereInteraction(_NA_PLUS_SPEC_,ptr,x,v,remaining,node,&sphere);
  bool valid=result==_PARTICLE_DELETED_ON_THE_FACE_ && collector.records.size()==1 && SameBytes(ptr,before) && remaining==7;
  if(valid){const auto& record=collector.records.front();const double number=node->block->GetLocalParticleWeight(_NA_PLUS_SPEC_)*
      PIC::ParticleBuffer::GetIndividualStatWeightCorrection(ptr);
    valid=record.incident.elementaryChargeC==SelectedCharge && record.incident.representedNumber==number &&
      record.chargeC==static_cast<long double>(number)*SelectedCharge && record.incident.impactEpochEtS==1003 &&
      record.incident.frame.geometryEpochEtS==2000;}
  Check(valid,"actual_sphere_collector_preserves_owner_and_uses_source_SI_charge_ledger");
  PIC::ParticleBuffer::DeleteParticle_withoutTrajectoryTermination(ptr);
}
}

int main(int argc,char** argv){
  if(argc!=5 || std::string(argv[1])!="mechanics")return 2;
  const std::string order(argv[2]),charge(argv[3]);
  if((order!="1" && order!="2") || (charge!="si" && charge!="legacy"))return 2;
  Order=std::stoi(order);SelectedCharge=charge=="si"?ExactSICharge:ElectronCharge;
  PIC::InitMPI();SetExitErrorCode(__LINE__,1);int infrastructure=0;
  try{
    if(PIC::ThisThread==0){CheckFile.open(std::string(argv[4])+".checks.csv");CheckFile<<"check,status\n";}
    KeyFile.open(std::string(argv[4])+".rank"+std::to_string(PIC::ThisThread)+".keys.csv");
    KeyFile<<std::setprecision(17)<<"context,original_parent,product_ordinal,prospective_uniform,optical_depth,product_clock_counter,owner\n";
    Setup();DirectCases();CapacityCase();MoverCase();SICollectorCase();
  }catch(const std::exception& exception){infrastructure=1;std::cerr<<"ERROR semantic cohort: "<<exception.what()<<'\n';}
  int bad=0,failed=0;MPI_Allreduce(&infrastructure,&bad,1,MPI_INT,MPI_MAX,MPI_GLOBAL_COMMUNICATOR);
  MPI_Allreduce(&Failures,&failed,1,MPI_INT,MPI_MAX,MPI_GLOBAL_COMMUNICATOR);
  if(PIC::ThisThread==0){CheckFile.flush();std::cout<<"M3 semantic cohort mode=mechanics order="<<Order<<" checks="<<Checks<<" status="<<(bad?"ERROR":failed?"FAIL":"PASS")<<'\n';}
  KeyFile.flush();MPI_Finalize();return bad?2:failed?1:0;
}
