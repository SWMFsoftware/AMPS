// Native adapter controls, NOT a lunar reaction/kinematic model. The declared
// stationary nominal-heavy-mass Na + He+ -> Na+ + He topology exercises the
// exact production plan/PIC publication path. No measured cross section,
// reaction rate, photon recoil or electron rest-mass convention is invented.
#include "pic.h"
#include "Moon.h"
#include "MoonReactionAdapter.h"
#include "MoonInput.h"
#include "MoonSurfaceEnvironment.h"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <new>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

// Allocation refusal is a local infrastructure fault control, never a
// replacement PIC allocator. All ordinary allocations use malloc/new as
// usual. The flag is enabled only at an explicitly staged commit boundary;
// it is consumed by ONE C++ allocation so error reporting remains possible.
namespace { bool RefuseNextAllocation=false; }
void* operator new(std::size_t bytes) {
  if(RefuseNextAllocation){RefuseNextAllocation=false;throw std::bad_alloc();}
  if(void* pointer=std::malloc(bytes?bytes:1))return pointer;
  throw std::bad_alloc();
}
void* operator new[](std::size_t bytes){return ::operator new(bytes);}
void operator delete(void* pointer) noexcept{std::free(pointer);}
void operator delete[](void* pointer) noexcept{std::free(pointer);}
void operator delete(void* pointer,std::size_t) noexcept{std::free(pointer);}
void operator delete[](void* pointer,std::size_t) noexcept{std::free(pointer);}

#if _INDIVIDUAL_PARTICLE_WEIGHT_MODE_ != _INDIVIDUAL_PARTICLE_WEIGHT_ON_
#error Native reaction adapter controls require individual weights ON
#endif
#if _SIMULATION_TIME_STEP_MODE_ != _SPECIES_DEPENDENT_LOCAL_TIME_STEP_
#error Native reaction adapter controls require actual species-local timesteps
#endif
#if _SIMULATION_PARTICLE_WEIGHT_MODE_ != _SPECIES_DEPENDENT_LOCAL_PARTICLE_WEIGHT_
#error Native reaction adapter controls require actual species-local weights
#endif

namespace R=Moon::Species::ReactionExecution;
namespace C=Moon::Chemistry;
namespace T=PIC::Mover::ForceCutSurface;
using Node=cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>;
using Byte=PIC::ParticleBuffer::byte;
namespace {
std::ofstream Checks,Products;
int failures=0,checks=0;
long int OpaqueOffset=-1;
Node* LocalNode=nullptr;
long int* CellHead=nullptr;

void Require(bool condition,const std::string& error) {
  if (!condition) throw std::runtime_error(error);
}
void Check(bool condition,const char* name) {
  ++checks;if (!condition)++failures;
  Checks<<PIC::ThisThread<<','<<name<<','<<(condition?"PASS":"FAIL")<<'\n';
  Checks.flush();
}
bool Near(double value,double expected) {
  // Explicit floating-point reconstruction allowance, not a science gate.
  return std::abs(value-expected)<=1e-14*std::max(std::abs(value),std::abs(expected));
}
std::vector<Byte> Snapshot(long int particle) {
  const auto* begin=PIC::ParticleBuffer::GetParticleDataPointer(particle);
  return std::vector<Byte>(begin,begin+PIC::ParticleBuffer::ParticleDataLength);
}
bool Same(long int particle,const std::vector<Byte>& bytes) {
  return PIC::ParticleBuffer::IsParticleAllocated(particle) &&
      std::memcmp(bytes.data(),PIC::ParticleBuffer::GetParticleDataPointer(particle),bytes.size())==0;
}
R::OpticalDepthClock Fresh(double uniform) {
  R::OpticalDepthClock clock;std::string error;
  Require(R::OpticalDepthClock::FromOpenUniform(uniform,&clock,&error),error);return clock;
}
R::TraceInterval Interval() {
  R::TraceInterval interval;std::string error;
  Require(R::TraceInterval::Create(1000,10,3,&interval,&error),error);return interval;
}
double UnitLoad(Node*) {return 1;}
double Resolution(double*) {return 1;}
void Collect(Node* node,std::vector<Node*>& leaves) {
  if(node->lastBranchFlag()==_BOTTOM_BRANCH_TREE_){leaves.push_back(node);return;}
  for(int child=0;child<8;++child)if(node->downNode[child])Collect(node->downNode[child],leaves);
}

void Setup() {
  PIC::Init_BeforeParser();
  C::RegisterReactionStorage();
  T::RegisterStorage(); // actual physical-epoch extension, before buffer init
  PIC::ParticleBuffer::RequestDataStorage(OpaqueOffset,16);
  double lower[3]={0,0,0},upper[3]={1,1,1};
  PIC::Mesh::mesh->AllowBlockAllocation=false;
  PIC::Mesh::mesh->init(lower,upper,Resolution);
  Require(PIC::Mesh::mesh->splitTreeNode(PIC::Mesh::mesh->rootTree),"root split failed");
  for(int child=0;child<8;++child)
    Require(PIC::Mesh::mesh->splitTreeNode(PIC::Mesh::mesh->rootTree->downNode[child]),"child split failed");
  PIC::Mesh::mesh->SetNodeNeibResolutionLevelLimit();
  PIC::Mesh::mesh->CreateBottomBranchNodeList(PIC::Mesh::mesh->rootTree);
  PIC::Mesh::mesh->SetParallelLoadMeasure(UnitLoad);
  PIC::Mesh::mesh->CreateNewParallelDistributionLists();
  PIC::Mesh::initCellSamplingDataBuffer();
  PIC::Mesh::mesh->AllowBlockAllocation=true;
  PIC::Mesh::mesh->AllocateTreeBlocks();
  PIC::Mesh::mesh->InitCellMeasure();
  PIC::Init_AfterParser();
  PIC::DomainBlockDecomposition::UpdateBlockTable();
  std::vector<Node*> leaves;Collect(PIC::Mesh::mesh->rootTree,leaves);
  for(Node* node:leaves)if(node->Thread==PIC::ThisThread && node->block) {
    for(int species=0;species<PIC::nTotalSpecies;++species) {
      node->block->SetLocalParticleWeight(7+species,species);
      node->block->SetLocalTimeStep(.25+.125*species,species);
    }
    node->block->SetLocalParticleWeight(2,_NA_SPEC_);
    node->block->SetLocalParticleWeight(3,_NA_PLUS_SPEC_);
    node->block->SetLocalParticleWeight(5,_HE_SPEC_);
    node->block->SetLocalTimeStep(.5,_NA_SPEC_);
    node->block->SetLocalTimeStep(.125,_NA_PLUS_SPEC_);
    node->block->SetLocalTimeStep(.0625,_HE_SPEC_);
    if(!LocalNode)LocalNode=node;
  }
  Require(LocalNode!=nullptr,"rank owns no actual allocated block");
  CellHead=&LocalNode->block->FirstCellParticleTable[0];
  PIC::ParticleBuffer::Init(64);
}

long int Particle(int species,bool linked=false) {
  const long int particle=linked?PIC::ParticleBuffer::GetNewParticle(*CellHead):PIC::ParticleBuffer::GetNewParticle();
  auto* data=PIC::ParticleBuffer::GetParticleDataPointer(particle);
  double position[3];
  // A real point in this block's first cell; no geometric location is inferred
  // from rank or ownership, and all products remain at this event position.
  const int cells[3]={_BLOCK_CELLS_X_,_BLOCK_CELLS_Y_,_BLOCK_CELLS_Z_};
  for(int axis=0;axis<3;++axis)
    position[axis]=LocalNode->xmin[axis]+.5*(LocalNode->xmax[axis]-LocalNode->xmin[axis])/cells[axis];
  double velocity[3]={0,0,0};
  PIC::ParticleBuffer::SetI(species,data);
  PIC::ParticleBuffer::SetX(position,data);PIC::ParticleBuffer::SetV(velocity,data);
  PIC::ParticleBuffer::SetIndividualStatWeightCorrection(.75,data);
  if(!linked){PIC::ParticleBuffer::SetNext(-1,data);PIC::ParticleBuffer::SetPrev(-1,data);}
  std::memset(data+OpaqueOffset,77,16);
  T::ResetParticleState(data);T::SetParticleEpochS(particle,1003);
  return particle;
}

struct Fixture {
  // Catalog identity is explicit even for a bath state. No kinetic meaning is
  // attributed to the compiled species' nominal mass equality.
  std::vector<R::SpeciesIdentity> catalog;
  R::Reaction reaction;
  R::ParentState parent;
  R::KinematicPolicy policy;
  long int particle,sentinel;
  Fixture():particle(Particle(_NA_SPEC_)),sentinel(Particle(_NE_SPEC_,true)) {
    const double na=PIC::MolecularData::GetMass(_NA_SPEC_);
    const double nai=PIC::MolecularData::GetMass(_NA_PLUS_SPEC_);
    const double he=PIC::MolecularData::GetMass(_HE_SPEC_);
    const double hei=PIC::MolecularData::GetMass(_HE_PLUS_SPEC_);
    Require(na==nai && he==hei,"analytical nominal-mass fixture requires actual equal neutral/ion heavy masses");
    catalog={{"Na","nominal-ground",_NA_SPEC_,{{"Na",1}},0,na},
      {"Na+","nominal-ground",_NA_PLUS_SPEC_,{{"Na",1}},1,nai},
      {"He","nominal-ground",_HE_SPEC_,{{"He",1}},0,he},
      {"He+","nominal-bath",_HE_PLUS_SPEC_,{{"He",1}},1,hei}};
    reaction={"ANALYTICAL-NOMINAL-MASS-STORAGE-NOT-CEX",{{0,1,R::Population::Tracked},{3,1,R::Population::Bath}},
      {{1,1,R::Population::Tracked},{2,1,R::Population::Tracked}}};
    parent.species=0;parent.referenceWeight=2;parent.individualCorrection=.75;
    double x[3];PIC::ParticleBuffer::GetX(x,particle);
    for(int axis=0;axis<3;++axis)parent.positionM[axis]=x[axis];
    policy.qualified=true;policy.identity="analytical-stationary-nominal-mass-control";
    policy.evidenceIdentity="NUMERICAL-ADAPTER-CONTROL-NOT-LUNAR-CHANNEL";
    policy.eventFrame="analytical-Cartesian-SI";policy.maximumClassicalSpeedMS=1;
    // Exact stationary momenta/energies; a rounding-level mass bound only.
    policy.bounds.modelMassKg={0,1e-14};policy.bounds.momentumKgMS={0,0};policy.bounds.energyJ={0,0};
    policy.bath.supplied=true;policy.bath.modelMassChangeKg=-hei;
    R::ProductKinematics a,b;a.term=0;a.referenceWeight=3;a.nextClock=Fresh(.2);
    b.term=1;b.referenceWeight=5;b.nextClock=Fresh(.8);policy.products={a,b};
  }
  R::EventPlan Plan(R::PlanScope scope=R::PlanScope::Transport) const {
    R::EventPlan plan;std::string error;
    Require(R::MakeEventPlan(catalog,reaction,parent,Interval(),scope,
        scope==R::PlanScope::Transport?&policy:nullptr,&plan,&error),error);return plan;
  }
};

bool ListValid(std::size_t expected) {
  std::set<long int> seen;long int previous=-1;
  for(long int ptr=*CellHead;ptr!=-1;ptr=PIC::ParticleBuffer::GetNext(ptr)) {
    if(ptr<0 || ptr>=PIC::ParticleBuffer::GetMaxNPart() || !PIC::ParticleBuffer::IsParticleAllocated(ptr) ||
        PIC::ParticleBuffer::GetPrev(ptr)!=previous || !seen.insert(ptr).second)return false;
    previous=ptr;
  }
  return seen.size()==expected;
}
void Clear() {
  // Normal list cleanup owns linked records. Explicitly release only detached
  // control records remaining after a negative case; never double-delete.
  PIC::ParticleBuffer::DeleteAllParticles();
  for(long int ptr=0;ptr<PIC::ParticleBuffer::GetMaxNPart();++ptr)
    if(PIC::ParticleBuffer::IsParticleAllocated(ptr))PIC::ParticleBuffer::DeleteParticle_withoutTrajectoryTermination(ptr);
  Require(PIC::ParticleBuffer::GetAllPartNum()==0 && *CellHead==-1,"fixture cleanup leaked slots or list entries");
}
void VerifyProducts(const R::EventPlan& plan,const std::vector<long int>& products,const char* context) {
  Check(products.size()==2,"complete_product_packet");
  for(std::size_t i=0;i<products.size();++i) {
    const auto& expected=plan.products().at(i);const long int ptr=products[i];
    const int compiled=PIC::ParticleBuffer::GetI(ptr);
    R::TraceInterval interval;R::OpticalDepthClock clock;std::size_t catalogIndex=999;double internal=-1;std::string error;
    const bool state=C::ReadReactionState(ptr,&interval,&clock,&catalogIndex,&internal,&error);
    const double weight=LocalNode->block->GetLocalParticleWeight(compiled);
    const double correction=PIC::ParticleBuffer::GetIndividualStatWeightCorrection(ptr);
    Check(compiled==expected.compiledSpecies && catalogIndex==expected.species,"actual_compiled_and_catalog_species");
    Check(state && interval.startS()==1000 && interval.durationS()==10 && interval.elapsedS()==3 && interval.remainingS()==7,
        "original_local_interval_not_daughter_dt");
    Check(state && clock.valid() && clock.fresh() && clock.remaining()==expected.clock.remaining() && internal==0,
        "fresh_product_clock_and_internal_state");
    Check(Near(weight*correction,1.5) && correction==expected.individualCorrection,"represented_number_without_dt_birth_factor");
    Check(T::GetParticleEpochS(ptr)==1003,"mover_epoch_matches_physical_event");
    double actualX[3],actualV[3];PIC::ParticleBuffer::GetX(actualX,ptr);PIC::ParticleBuffer::GetV(actualV,ptr);
    bool kinematics=true;
    for(int axis=0;axis<3;++axis)
      kinematics=kinematics && actualX[axis]==expected.positionM[axis] && actualV[axis]==expected.velocityMS[axis];
    Check(kinematics,"actual_product_position_and_velocity");
    bool opaque=true;auto* data=PIC::ParticleBuffer::GetParticleDataPointer(ptr);
    for(int j=0;j<16;++j)opaque=opaque && data[OpaqueOffset+j]==77;
    Check(opaque,"opaque_lineage_attributes_preserved");
    Products<<PIC::ThisThread<<','<<context<<','<<i<<','<<compiled<<','<<catalogIndex<<','
      <<PIC::MolecularData::GetMass(compiled)<<','<<PIC::MolecularData::GetElectricCharge(compiled)<<','
      <<weight<<','<<correction<<','<<weight*correction<<','<<LocalNode->block->GetLocalTimeStep(compiled)<<','
      <<interval.startS()<<','<<interval.durationS()<<','<<interval.elapsedS()<<','<<interval.remainingS()<<','<<clock.remaining()<<','
      <<actualX[0]<<','<<actualX[1]<<','<<actualX[2]<<','<<actualV[0]<<','<<actualV[1]<<','<<actualV[2]<<'\n';
  }
}

void Cases() {
  Check(PIC::nTotalSpecies==8 && C::ReactionStorageRegistered() && T::IsStorageRegistered(),"actual_eight_species_and_registered_storage");
  {
    const int neutrals[]={_NA_SPEC_,_HE_SPEC_,_NE_SPEC_,_AR_SPEC_};
    const int ions[]={_NA_PLUS_SPEC_,_HE_PLUS_SPEC_,_NE_PLUS_SPEC_,_AR_PLUS_SPEC_};
    std::set<int> actual;bool valid=true;
    for(int i=0;i<4;++i) {
      actual.insert(neutrals[i]);actual.insert(ions[i]);
      valid=valid && neutrals[i]>=0 && ions[i]>=0 &&
        std::isfinite(PIC::MolecularData::GetMass(neutrals[i])) && PIC::MolecularData::GetMass(neutrals[i])>0 &&
        std::isfinite(PIC::MolecularData::GetMass(ions[i])) && PIC::MolecularData::GetMass(ions[i])>0 &&
        PIC::MolecularData::GetElectricCharge(neutrals[i])==0 &&
        PIC::MolecularData::GetElectricCharge(ions[i])==ElectronCharge;
    }
    Check(valid && actual.size()==8,"actual_atomic_mass_charge_and_species_map");
  }
  Check(LocalNode->Thread==PIC::ThisThread && LocalNode->block && *CellHead==-1,"real_owned_block_and_empty_cell_list");
  Check(LocalNode->block->GetLocalTimeStep(_NA_SPEC_)==.5 &&
      LocalNode->block->GetLocalTimeStep(_NA_PLUS_SPEC_)==.125 &&
      LocalNode->block->GetLocalTimeStep(_HE_SPEC_)==.0625,"actual_unequal_species_local_timesteps");
  {
    Fixture f;const auto plan=f.Plan();R::PendingEvent pending(plan);C::ReactionLedger ledger;
    const auto before=Snapshot(f.particle);const long int head=*CellHead;
    std::vector<long int> filler;
    while(PIC::ParticleBuffer::GetMaxNPart()-PIC::ParticleBuffer::GetAllPartNum()>1)
      filler.push_back(Particle(_AR_SPEC_));
    const long int allocated=PIC::ParticleBuffer::GetAllPartNum();std::string error;
    {
      C::ParticleTransaction backend(f.particle,*CellHead,LocalNode,f.catalog,&ledger);
      Check(pending.TryCommit(&backend,&error)==R::AttemptStatus::Deferred,"complete_packet_capacity_deferral");
      Check(Same(f.particle,before) && *CellHead==head && ledger.events.empty() &&
          PIC::ParticleBuffer::GetAllPartNum()==allocated && ListValid(1),"deferral_preserves_parent_list_count_and_ledger");
      for(long int ptr:filler)PIC::ParticleBuffer::DeleteParticle_withoutTrajectoryTermination(ptr);
      Check(pending.TryCommit(&backend,&error)==R::AttemptStatus::Committed,"same_plan_retry_after_real_capacity_release");
      Check(!PIC::ParticleBuffer::IsParticleAllocated(f.particle) && ListValid(3) &&
          PIC::ParticleBuffer::GetAllPartNum()==3 && ledger.events.size()==1,"atomic_parent_delete_product_list_and_ledger");
      VerifyProducts(pending.plan(),backend.committedProducts(),"capacity-retry");
      Check(ledger.events[0].ledger.integerConservationVerified && ledger.events[0].ledger.kinematicConservationVerified &&
          ledger.events[0].ledger.reactantChargeElementary==ledger.events[0].ledger.productChargeElementary,
          "committed_exact_integer_and_kinematic_ledger");
      Check(pending.TryCommit(&backend,&error)==R::AttemptStatus::Error && ledger.events.size()==1 &&
          PIC::ParticleBuffer::GetAllPartNum()==3,"double_commit_refused_without_publication");
    }
    Check(ListValid(3) && PIC::ParticleBuffer::GetAllPartNum()==3,"destructor_retains_committed_products");
  }Clear();
  {
    Fixture f;C::ReactionLedger ledger;std::string error;
    // Explicit analytical elastic mechanics; NOT an angular/CEX sampler.
    // Define the incoming bath velocity first. Its signed momentum/energy
    // removal is independent input, not a residual derived from daughters.
    const double bathVelocity=-1,parentVelocity=1;
    const double sodium=f.catalog[0].modelMassKg,helium=f.catalog[3].modelMassKg;
    f.parent.velocityMS[0]=parentVelocity;
    double initial[3]={parentVelocity,0,0};PIC::ParticleBuffer::SetV(initial,f.particle);
    f.policy.bath.momentumChangeKgMS[0]=-helium*bathVelocity;
    f.policy.bath.energyChangeJ=-.5*helium*bathVelocity*bathVelocity;
    f.policy.maximumClassicalSpeedMS=10;
    f.policy.bounds.momentumKgMS={0,1e-14};f.policy.bounds.energyJ={0,1e-14};
    f.policy.products[0].velocityMS[0]=((sodium-helium)*parentVelocity+2*helium*bathVelocity)/(sodium+helium);
    f.policy.products[1].velocityMS[0]=(2*sodium*parentVelocity+(helium-sodium)*bathVelocity)/(sodium+helium);
    f.policy.identity="analytical-elastic-nominal-mass-control";
    const auto plan=f.Plan();R::PendingEvent pending(plan);
    C::ParticleTransaction backend(f.particle,*CellHead,LocalNode,f.catalog,&ledger);
    Check(pending.TryCommit(&backend,&error)==R::AttemptStatus::Committed && ledger.events.size()==1 && ListValid(3),
        "native_nonzero_analytical_mechanical_commit");
    VerifyProducts(plan,backend.committedProducts(),"elastic-control");
  }Clear();
  {
    // Both transactions share the real accepted-event ledger and have private
    // native slots. Publishing the second consumes the first's old capacity
    // reservation. A forced allocation failure at the first commit MUST occur
    // before any parent/list/product/physical-ledger mutation.
    Fixture a,b;const auto pa=a.Plan(),pb=b.Plan();C::ReactionLedger ledger;
    const auto before=Snapshot(a.particle);std::string error;
    C::ParticleTransaction first(a.particle,*CellHead,LocalNode,a.catalog,&ledger);
    C::ParticleTransaction second(b.particle,*CellHead,LocalNode,b.catalog,&ledger);
    Require(first.Reserve(pa,&error)==R::ReservationStatus::Ready,error);Require(first.Stage(pa,&error),error);
    Require(second.Reserve(pb,&error)==R::ReservationStatus::Ready,error);Require(second.Stage(pb,&error),error);
    Require(second.Commit(pb,&error),error);
    Require(ledger.events.size()==1 && ledger.events.capacity()==1,
        "native allocation-fault control requires the actual single-entry reserved ledger capacity");
    bool refused=false;RefuseNextAllocation=true;
    try { first.Commit(pa,&error); }
    catch(const std::bad_alloc&) { refused=true; }
    RefuseNextAllocation=false;first.Rollback();
    Check(refused && Same(a.particle,before) && ledger.events.size()==1 && ListValid(4) &&
        PIC::ParticleBuffer::GetAllPartNum()==5,"interleaved_ledger_allocation_failure_before_physical_publish");
    R::PendingEvent retry(pa);
    Check(retry.TryCommit(&first,&error)==R::AttemptStatus::Committed && ledger.events.size()==2 &&
        !PIC::ParticleBuffer::IsParticleAllocated(a.particle) && ListValid(6) && PIC::ParticleBuffer::GetAllPartNum()==6,
        "interleaved_shared_ledger_same_plan_retry");
  }Clear();
  {
    Fixture f;const auto plan=f.Plan();C::ReactionLedger ledger;const auto before=Snapshot(f.particle);std::string error;
    {
      C::ParticleTransaction backend(f.particle,*CellHead,LocalNode,f.catalog,&ledger);
      Check(backend.Reserve(plan,&error)==R::ReservationStatus::Ready && PIC::ParticleBuffer::GetAllPartNum()==4,
          "real_private_slots_reserved_before_publish");
      Check(backend.Stage(plan,&error) && ListValid(1) && ledger.events.empty() && Same(f.particle,before),
          "staged_records_unpublished_parent_unchanged");
      backend.Rollback();
      Check(PIC::ParticleBuffer::GetAllPartNum()==2 && ListValid(1) && ledger.events.empty() && Same(f.particle,before),
          "explicit_rollback_releases_all_private_slots");
    }
    Check(PIC::ParticleBuffer::GetAllPartNum()==2 && Same(f.particle,before),"rollback_destructor_idempotent");
  }Clear();
  {
    Fixture f;const auto plan=f.Plan();C::ReactionLedger ledger;const auto before=Snapshot(f.particle);std::string error;
    C::ParticleTransaction backend(f.particle,*CellHead,LocalNode,f.catalog,&ledger);
    Require(backend.Reserve(plan,&error)==R::ReservationStatus::Ready,error);Require(backend.Stage(plan,&error),error);
    auto* data=PIC::ParticleBuffer::GetParticleDataPointer(f.particle);data[OpaqueOffset]=78;
    Check(!backend.Commit(plan,&error) && ledger.events.empty() && ListValid(1),"changed_opaque_parent_bytes_refuse_commit");
    backend.Rollback();data[OpaqueOffset]=77;
    Check(Same(f.particle,before) && PIC::ParticleBuffer::GetAllPartNum()==2,"failed_commit_rollback_preserves_restored_parent");
    R::PendingEvent pending(plan);
    Check(pending.TryCommit(&backend,&error)==R::AttemptStatus::Committed,"same_plan_retry_after_parent_invariant_restored");
    VerifyProducts(plan,backend.committedProducts(),"opaque-retry");
  }Clear();
  {
    Fixture f;const auto plan=f.Plan();C::ReactionLedger ledger;const auto before=Snapshot(f.particle);std::string error;
    auto bad=f.catalog;bad[1].state="changed-after-planning";
    C::ParticleTransaction backend(f.particle,*CellHead,LocalNode,bad,&ledger);
    Check(backend.Reserve(plan,&error)==R::ReservationStatus::Error && Same(f.particle,before) &&
        PIC::ParticleBuffer::GetAllPartNum()==2 && ledger.events.empty(),"stale_catalog_state_identity_refused");
    LocalNode->block->SetLocalParticleWeight(2.5,_NA_SPEC_);
    C::ParticleTransaction weightBackend(f.particle,*CellHead,LocalNode,f.catalog,&ledger);
    Check(weightBackend.Reserve(plan,&error)==R::ReservationStatus::Error && Same(f.particle,before),"stale_reference_weight_refused");
    LocalNode->block->SetLocalParticleWeight(2,_NA_SPEC_);
    double x[3];PIC::ParticleBuffer::GetX(x,f.particle);const double original=x[0];x[0]+=.001;PIC::ParticleBuffer::SetX(x,f.particle);
    Check(weightBackend.Reserve(plan,&error)==R::ReservationStatus::Error && ledger.events.empty(),"stale_parent_position_refused");
    x[0]=original;PIC::ParticleBuffer::SetX(x,f.particle);
  }Clear();
  {
    Fixture f;const auto plan=f.Plan();C::ReactionLedger ledger;const auto before=Snapshot(f.particle);std::string error;
    std::vector<long int> detached;
    {
      C::ParticleTransaction backend(f.particle,*CellHead,LocalNode,f.catalog,&ledger,C::ProductPlacement::DetachedContinuation,true);
      R::PendingEvent pending(plan);
      Check(pending.TryCommit(&backend,&error)==R::AttemptStatus::Committed && Same(f.particle,before) && ListValid(1),
          "detached_retain_parent_preserves_mover_ownership");
      detached=backend.committedProducts();bool unlinked=detached.size()==2;
      for(long int ptr:detached)unlinked=unlinked && PIC::ParticleBuffer::GetNext(ptr)==-1 && PIC::ParticleBuffer::GetPrev(ptr)==-1;
      Check(unlinked && PIC::ParticleBuffer::GetAllPartNum()==4,"detached_products_allocated_but_unlinked");
      VerifyProducts(plan,detached,"detached");
    }
    Check(Same(f.particle,before) && PIC::ParticleBuffer::GetAllPartNum()==4,"detached_destructor_does_not_delete_parent_or_products");
    PIC::ParticleBuffer::DeleteParticle(f.particle);
    Check(PIC::ParticleBuffer::GetAllPartNum()==3 && ledger.events.size()==1,"mover_owner_deletes_retained_parent_exactly_once");
    for(long int ptr:detached)PIC::ParticleBuffer::DeleteParticle_withoutTrajectoryTermination(ptr);
    Check(PIC::ParticleBuffer::GetAllPartNum()==1 && ListValid(1),"detached_continuation_owner_can_release_products");
  }Clear();
  {
    Fixture f;C::ReactionLedger ledger;std::string error;const auto before=Snapshot(f.particle);
    const auto plan=f.Plan();C::ParticleTransaction invalid(f.particle,*CellHead,LocalNode,f.catalog,&ledger,C::ProductPlacement::RebuiltCellList,true);
    Check(invalid.Reserve(plan,&error)==R::ReservationStatus::Error && Same(f.particle,before),"retain_parent_requires_detached_placement");
    R::PendingEvent rateOnly(f.Plan(R::PlanScope::RateOnly));
    C::ParticleTransaction backend(f.particle,*CellHead,LocalNode,f.catalog,&ledger);
    Check(rateOnly.TryCommit(&backend,&error)==R::AttemptStatus::Error && PIC::ParticleBuffer::GetAllPartNum()==2 &&
        ledger.events.empty() && Same(f.particle,before),"rate_only_cannot_publish_pic_particles");
  }Clear();
  Check(PIC::ParticleBuffer::GetAllPartNum()==0 && *CellHead==-1,"all_native_slots_and_cell_lists_clean");
}
}

namespace {
bool CollectorState(long int particle,void*,void*,std::size_t* state,std::string* error) {
  if(PIC::ParticleBuffer::GetI(particle)!=_NA_PLUS_SPEC_){*error="analytical collector fixture has only its declared Na+ state";return false;}
  *state=0;return true;
}
void CollectorCases(const std::string& prefix) {
  // Actual production dispatcher and PIC weight/epoch adapter, with a named
  // nominal analytical state. This fixture supplies no lunar neutralization,
  // formation-energy measurement or physical retention/release model.
  namespace I=Moon::Species::IonCollector;
  Moon::Runtime::Configuration input;
  input.particleTracingMover="force-cut-surface-first-order";
  input.ionSurfaceBoundary="ideal-collector";
  input.simulationEnvironmentMode="frozen";input.tracingFieldFrame="lso";
  input.tracingRelativeTolerance=1e-10;input.tracingAbsolutePositionToleranceM=1e-10;
  input.tracingAbsoluteVelocityToleranceMS=1e-10;input.tracingGeometryToleranceM=1e-10;
  input.tracingTimeToleranceS=1e-10;input.tracingMinimumSubstepS=1e-6;input.tracingMaximumSubstepS=1;
  input.tracingMaximumSubsteps=1000;input.tracingMaximumRefinements=20;
  input.tracingMaximumInteractions=100;input.tracingMaximumSamePointInteractions=10;
  std::string error;Require(Moon::Runtime::InstallConfiguration(input,&error),error);
  I::StateDefinition state;
  state.material.identity="Na+";state.material.state="analytical-nominal-ion-control";
  state.material.compiledSpecies=_NA_PLUS_SPEC_;state.material.nuclei={{"Na",1}};
  state.material.elementaryCharge=1;state.material.modelMassKg=PIC::MolecularData::GetMass(_NA_PLUS_SPEC_);
  state.formationEnergyJ=3e-20;state.excitationEnergyJ=0;
  state.energyReferenceIdentity="analytical-reference-zero-not-lunar-thermochemistry";
  state.sourceEvidenceIdentity="declared-native-control-not-measured-ion-surface-data";
  I::Ledger ledger;
  auto bad=state;bad.material.modelMassKg*=2;
  Check(!Moon::SurfaceEnvironment::InstallProductionIonCollector({bad},CollectorState,nullptr,&ledger,&error)&&
      !Moon::SurfaceEnvironment::ProductionIonCollectorInstalled(),"collector_rejects_actual_pic_mass_mismatch_before_install");
  Check(Moon::SurfaceEnvironment::InstallProductionIonCollector({state},CollectorState,nullptr,&ledger,&error),
      "collector_installs_explicit_nominal_state_resolver");
  Check(!Moon::SurfaceEnvironment::InstallProductionIonCollector({state},CollectorState,nullptr,&ledger,&error),
      "collector_reinstall_refused");
  // Identity synchronized pose is an analytical geometry control. Frozen
  // geometry ET and physical tracer impact ET deliberately differ.
  for(int i=0;i<6;++i)for(int j=0;j<6;++j){
    Moon::OrbitalMotion::SO_to_IAU_TransformationMartix[i][j]=i==j?1:0;
    Moon::OrbitalMotion::IAU_to_SO_TransformationMartix[i][j]=i==j?1:0;
  }
  Moon::OrbitalMotion::et=2000;
  const long int particle=Particle(_NA_PLUS_SPEC_);double x[3],v[3];
  PIC::ParticleBuffer::GetX(x,particle);double norm2=0;
  for(int k=0;k<3;++k)norm2+=x[k]*x[k];const double radius=std::sqrt(norm2);
  for(int k=0;k<3;++k)v[k]=-2*x[k]/radius;
  PIC::ParticleBuffer::SetV(v,particle);const auto before=Snapshot(particle);
  cInternalSphericalData sphere;double center[3]={0,0,0};sphere.SetSphereGeometricalParameters(center,radius);
  double remaining=7;
  const int outcome=Moon::IonAwareSphereInteraction(_NA_PLUS_SPEC_,particle,x,v,remaining,LocalNode,&sphere);
  Check(outcome==_PARTICLE_DELETED_ON_THE_FACE_&&ledger.records.size()==1,"actual_sphere_dispatch_records_collector_and_returns_delete");
  Check(Same(particle,before)&&remaining==7,"collector_callback_preserves_mover_owned_particle_and_remaining_interval");
  Require(ledger.records.size()==1,"missing collector transfer");const auto& record=ledger.records.front();
  const double represented=LocalNode->block->GetLocalParticleWeight(_NA_PLUS_SPEC_)*.75;
  Check(record.incident.representedNumber==represented&&Near(static_cast<double>(record.modelMassKg),represented*state.material.modelMassKg),
      "collector_uses_actual_local_weight_and_individual_correction");
  Check(record.representedChargeElementary==represented&&record.nuclei.size()==1&&record.nuclei.front().representedNumber==represented,
      "collector_actual_charge_and_nuclear_inventory");
  Check(record.incident.impactEpochEtS==1003&&record.incident.frame.geometryEpochEtS==2000,
      "collector_distinguishes_actual_impact_and_frozen_geometry_epochs");
  Check(record.incident.frame.geometryIdentity=="analytic-sphere-ideal-collector"&&record.incident.frame.stateTransformsSupplied,
      "collector_actual_sphere_frame_receipt");
  const double expectedK=represented*.5*state.material.modelMassKg*4;
  Check(Near(static_cast<double>(record.kineticEnergyFacetJ),expectedK)&&
      Near(static_cast<double>(record.totalEnergyFacetJ),expectedK+represented*state.formationEnergyJ),
      "collector_separates_kinetic_and_formation_energy_without_heat_claim");
  const long int allocated=PIC::ParticleBuffer::GetAllPartNum();PIC::ParticleBuffer::DeleteParticle_withoutTrajectoryTermination(particle);
  Check(PIC::ParticleBuffer::GetAllPartNum()==allocated-1,"collector_caller_deletes_mover_owned_incident_exactly_once");
  const long int neutral=Particle(_NA_SPEC_);const auto neutralBefore=Snapshot(neutral);
  I::ImpactFrame unused;
  Check(!Moon::SurfaceEnvironment::HandleProductionIonImpact(neutral,x,v,LocalNode,unused)&&Same(neutral,neutralBefore)&&ledger.records.size()==1,
      "collector_neutral_bypass_preserves_bytes_and_ledger");
  PIC::ParticleBuffer::DeleteParticle_withoutTrajectoryTermination(neutral);
  Check(PIC::ParticleBuffer::GetAllPartNum()==0&&*CellHead==-1,"collector_native_slots_and_lists_clean");
  std::ofstream output(prefix+".collector.csv");output<<std::setprecision(21)
    <<"rank,represented_number,charge_elementary,nuclei_number,model_mass_kg,kinetic_energy_facet_J,formation_energy_J,total_energy_facet_J,impact_epoch_s,geometry_epoch_s\n"
    <<PIC::ThisThread<<','<<represented<<','<<record.representedChargeElementary<<','<<record.nuclei.front().representedNumber<<','
    <<record.modelMassKg<<','<<record.kineticEnergyFacetJ<<','<<record.formationEnergyJ<<','<<record.totalEnergyFacetJ<<','
    <<record.incident.impactEpochEtS<<','<<record.incident.frame.geometryEpochEtS<<'\n';
  output.flush();Require(output.good(),"failed actual collector evidence output");
}
}

int main(int argc,char** argv) {
  PIC::InitMPI();SetExitErrorCode(__LINE__,1);
  if(argc!=2){MPI_Finalize();return 2;}
  const std::string prefix=std::string(argv[1])+".rank"+std::to_string(PIC::ThisThread);
  Checks.open(prefix+".checks.csv");Products.open(prefix+".products.csv");
  Checks<<"rank,check,status\n";Products<<std::setprecision(17)
    <<"rank,context,product,compiled_species,catalog_species,mass_kg,charge_c,reference_weight,individual_correction,represented_number,local_dt_s,start_s,duration_s,elapsed_s,remaining_s,optical_depth,x_m,y_m,z_m,vx_m_s,vy_m_s,vz_m_s\n";
  int error=0;
  try {
    Require(Checks.good() && Products.good(),"cannot open native evidence streams");Setup();
    std::ofstream species(prefix+".species.csv");species<<std::setprecision(17)<<"rank,compiled_species,mass_kg,charge_c\n";
    for(int id=0;id<PIC::nTotalSpecies;++id)
      species<<PIC::ThisThread<<','<<id<<','<<PIC::MolecularData::GetMass(id)<<','<<PIC::MolecularData::GetElectricCharge(id)<<'\n';
    species.flush();Require(species.good(),"failed species inventory output");Cases();CollectorCases(prefix);
    Checks.flush();Products.flush();Require(Checks.good() && Products.good(),"failed native evidence output");
  } catch(const std::exception& exception) {
    error=1;std::cerr<<"ERROR native adapter rank="<<PIC::ThisThread<<": "<<exception.what()<<'\n';
  }
  int globalError=0,globalFailures=0;
  MPI_Allreduce(&error,&globalError,1,MPI_INT,MPI_MAX,MPI_GLOBAL_COMMUNICATOR);
  MPI_Allreduce(&failures,&globalFailures,1,MPI_INT,MPI_SUM,MPI_GLOBAL_COMMUNICATOR);
  std::cout<<"M3 native particle backend rank="<<PIC::ThisThread<<" checks="<<checks<<" status="
    <<(error?"ERROR":failures?"FAIL":"PASS")<<'\n';
  MPI_Finalize();return globalError?2:globalFailures?1:0;
}
