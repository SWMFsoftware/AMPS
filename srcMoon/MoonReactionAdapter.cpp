#include "MoonReactionAdapter.h"

#include <cmath>
#include <cstring>
#include <limits>
#include <type_traits>
#include <utility>

namespace Moon { namespace Chemistry {
namespace {
long int stateOffset=-1;
struct ParticleReactionState {
  bool valid=false;
  std::size_t catalogSpecies=0;
  double internalEnergyJ=0.0; // Joule per physical product under its state convention.
  Execution::TraceInterval interval;
  Execution::OpticalDepthClock clock;
};
static_assert(std::is_trivially_copyable<ParticleReactionState>::value,
    "MPI particle bytes must contain no owning pointers");
static_assert(std::is_nothrow_move_constructible<CommittedReaction>::value,
    "Publishing a pre-reserved ledger entry must not allocate or throw");

bool Reject(std::string* error,const char* message) {
  if (error) *error=message;
  return false;
}
}

void RegisterReactionStorage() {
  if (stateOffset>=0) return;
  if (PIC::ParticleBuffer::ParticleDataBuffer!=nullptr)
    exit(__LINE__,__FILE__,"register Moon chemistry state before initializing the particle buffer");
  PIC::ParticleBuffer::RequestDataStorage(stateOffset,sizeof(ParticleReactionState));
}
bool ReactionStorageRegistered() { return stateOffset>=0; }

bool ReadReactionState(long int particle,Execution::TraceInterval* interval,
    Execution::OpticalDepthClock* clock,std::size_t* species,
    double* internalEnergyJ,std::string* error) {
  if (!error || !interval || !clock || !species || !internalEnergyJ) return false;
  if (stateOffset<0 || particle<0 || particle>=PIC::ParticleBuffer::GetMaxNPart() ||
      !PIC::ParticleBuffer::IsParticleAllocated(particle))
    return Reject(error,"reaction state needs registered storage and an allocated particle");
  ParticleReactionState state;
  std::memcpy(&state,PIC::ParticleBuffer::GetParticleDataPointer(particle)+stateOffset,sizeof(state));
  if (!state.valid || !state.interval.valid() || !state.clock.valid())
    return Reject(error,"particle has no initialized physical reaction interval/clock");
  *interval=state.interval;*clock=state.clock;*species=state.catalogSpecies;
  *internalEnergyJ=state.internalEnergyJ;
  error->clear();return true;
}

ParticleTransaction::ParticleTransaction(long int parent,long int& head,
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node,
    const std::vector<Execution::SpeciesIdentity>& catalog,ReactionLedger* ledger,
    ProductPlacement placement,bool retainParent,double elementaryChargeC,
    StagedProductInitializer initializeProduct,void* initializerUser):parent_(parent),head_(head),
    node_(node),catalog_(catalog),ledger_(ledger),placement_(placement),retainParent_(retainParent),
    elementaryChargeC_(elementaryChargeC),initializeProduct_(initializeProduct),initializerUser_(initializerUser) {}
ParticleTransaction::~ParticleTransaction() { Rollback(); }

bool ParticleTransaction::ValidateParent(const Execution::EventPlan& plan,
    std::string* error) const {
#if _COMPILATION_MODE_ != _COMPILATION_MODE__MPI_
  return Reject(error,"Moon conservative reaction transactions require MPI-only/OMP=1 execution");
#elif _INDIVIDUAL_PARTICLE_WEIGHT_MODE_ != _INDIVIDUAL_PARTICLE_WEIGHT_ON_
  return Reject(error,"Moon conservative reaction products require individual particle weights ON");
#elif _PIC_PARTICLE_TRACKER_MODE_ != _PIC_MODE_OFF_
  return Reject(error,"Moon atomic reaction publication currently requires particle trajectory tracker OFF");
#elif _PIC_LOCAL_TIME_STEP_FLUX_RESAMPLING_MODE_ == _PIC_MODE_ON_
  return Reject(error,"physical reaction intervals cannot use stationary-flux W/dt resampling");
#else
  if (!std::isfinite(elementaryChargeC_) || elementaryChargeC_<=0 ||
      !plan.valid() || plan.scope()!=Execution::PlanScope::Transport ||
      stateOffset<0 || !node_ || !node_->block || !ledger_ ||
      parent_<0 || parent_>=PIC::ParticleBuffer::GetMaxNPart() || head_==parent_ ||
      !PIC::ParticleBuffer::IsParticleAllocated(parent_) ||
      plan.parent().species>=catalog_.size() ||
      (head_!=-1 && (head_<0 || head_>=PIC::ParticleBuffer::GetMaxNPart() ||
                    !PIC::ParticleBuffer::IsParticleAllocated(head_))) ||
      (retainParent_ && placement_!=ProductPlacement::DetachedContinuation) ||
      (placement_!=ProductPlacement::RebuiltCellList && placement_!=ProductPlacement::DetachedContinuation))
    return Reject(error,"invalid transport plan, storage, detached parent, node or reaction ledger");
  const int species=catalog_[plan.parent().species].compiledSpecies;
  const auto& parentIdentity=catalog_[plan.parent().species];
  if (species<0 || species>=PIC::nTotalSpecies || PIC::ParticleBuffer::GetI(parent_)!=species ||
      species!=plan.parentCompiledSpecies() || parentIdentity.identity!=plan.parentSpeciesIdentity() ||
      parentIdentity.state!=plan.parentStateIdentity() ||
      PIC::MolecularData::GetMass(species)!=catalog_[plan.parent().species].modelMassKg ||
      PIC::MolecularData::GetElectricCharge(species)!=
          catalog_[plan.parent().species].elementaryCharge*elementaryChargeC_ ||
      node_->block->GetLocalParticleWeight(species)!=plan.parent().referenceWeight ||
      PIC::ParticleBuffer::GetIndividualStatWeightCorrection(parent_)!=plan.parent().individualCorrection)
    return Reject(error,"reaction parent species/reference/individual weight changed after planning");
  double x[3],v[3];PIC::ParticleBuffer::GetX(x,parent_);PIC::ParticleBuffer::GetV(v,parent_);
  for (int axis=0;axis<3;++axis)
    if (x[axis]!=plan.parent().positionM[axis] || v[axis]!=plan.parent().velocityMS[axis])
      return Reject(error,"reaction parent position/velocity changed after planning");
  for (const auto& product:plan.products()) {
    if (product.species>=catalog_.size())
      return Reject(error,"product state is outside the approved species catalog");
    const int id=catalog_[product.species].compiledSpecies;
    if (id<0 || id>=PIC::nTotalSpecies ||
        id!=product.compiledSpecies ||
        catalog_[product.species].identity!=product.speciesIdentity ||
        catalog_[product.species].state!=product.stateIdentity ||
        node_->block->GetLocalParticleWeight(id)!=product.referenceWeight ||
        PIC::MolecularData::GetMass(id)!=catalog_[product.species].modelMassKg ||
        PIC::MolecularData::GetElectricCharge(id)!=
            catalog_[product.species].elementaryCharge*elementaryChargeC_ ||
        !std::isfinite(product.interval.startS()+product.interval.elapsedS()))
      return Reject(error,"product species/reference weight is not the actual initialized PIC value");
  }
  return true;
#endif
}

Execution::ReservationStatus ParticleTransaction::Reserve(const Execution::EventPlan& plan,
    std::string* error) {
  if (!error) return Execution::ReservationStatus::Error;
  if (committed_ || !reserved_.empty() || staged_) {
    *error="reaction reservation cannot reuse a staged or committed transaction";
    return Execution::ReservationStatus::Error;
  }
  if (!ValidateParent(plan,error))
    return Execution::ReservationStatus::Error;
  const auto count=plan.products().size();
  const long int available=PIC::ParticleBuffer::GetMaxNPart()-PIC::ParticleBuffer::GetAllPartNum();
  if (available<0 || count>static_cast<std::size_t>(available)) {
    *error="insufficient particle slots for the COMPLETE reaction packet; no physical "
        "state was changed. Increase particle-buffer capacity or reduce the initial "
        "model-particle count/weight resolution; never discard products or alter rates";
    return Execution::ReservationStatus::Deferred;
  }
  // Reserve every potentially allocating C++ container before touching the
  // PIC allocator. PendingEvent catches allocation exceptions and rolls back.
  reserved_.reserve(count);committedProducts_.reserve(count);
  parentSnapshot_.resize(PIC::ParticleBuffer::ParticleDataLength);
  std::memcpy(parentSnapshot_.data(),PIC::ParticleBuffer::GetParticleDataPointer(parent_),parentSnapshot_.size());
  ledger_->events.reserve(ledger_->events.size()+1);
  stagedLedger_.identity=plan.reactionIdentity();stagedLedger_.ledger=plan.ledger();
  if(initializeProduct_)stagedLedger_.initialProducts.resize(count);
  for (std::size_t i=0;i<count;++i) reserved_.push_back(PIC::ParticleBuffer::GetNewParticle());
  error->clear();return Execution::ReservationStatus::Ready;
}

bool ParticleTransaction::Stage(const Execution::EventPlan& plan,std::string* error) {
  if (!error || committed_ || staged_ || reserved_.size()!=plan.products().size() ||
      parentSnapshot_.empty() || !ValidateParent(plan,error)) return false;
  for (std::size_t i=0;i<reserved_.size();++i) {
    const auto& product=plan.products()[i];const long int ptr=reserved_[i];
    auto* bytes=PIC::ParticleBuffer::GetParticleDataPointer(ptr);
    // Preserve source lineage/opaque application attributes, but clear list
    // links and the parent's mover continuation. A daughter is a NEW physical
    // trajectory at this event, whereas ordinary MPI clones preserve time.
    PIC::ParticleBuffer::CloneParticle(bytes,parentSnapshot_.data());
    PIC::ParticleBuffer::SetPrev(-1,bytes);PIC::ParticleBuffer::SetNext(-1,bytes);
    PIC::ParticleBuffer::SetI(catalog_[product.species].compiledSpecies,bytes);
    double x[3],v[3];for(int k=0;k<3;++k){x[k]=product.positionM[k];v[k]=product.velocityMS[k];}
    PIC::ParticleBuffer::SetX(x,bytes);PIC::ParticleBuffer::SetV(v,bytes);
    PIC::ParticleBuffer::SetIndividualStatWeightCorrection(product.individualCorrection,bytes);
    PIC::Mover::ForceCutSurface::ResetParticleState(bytes);
    if (PIC::Mover::ForceCutSurface::IsStorageRegistered())
      PIC::Mover::ForceCutSurface::SetParticleEpochS(ptr,product.interval.startS()+product.interval.elapsedS());
    ParticleReactionState state;state.valid=true;state.catalogSpecies=product.species;
    state.internalEnergyJ=product.internalEnergyJ;
    state.interval=product.interval;state.clock=product.clock;
    std::memcpy(bytes+stateOffset,&state,sizeof(state));
    // Full semantic genealogy/counters belong to a NEW physical history, not
    // the copied parent's key. This remains rollback-safe before Commit.
    if (initializeProduct_) {
      if(!initializeProduct_(plan,i,parent_,ptr,initializerUser_,error)||
          !ValidateInitializedProduct(plan,i,&stagedLedger_.initialProducts[i],error)) return false;
    }
  }
  stagedPlan_=&plan;staged_=true;error->clear();return true;
}

bool ParticleTransaction::ValidateInitializedProduct(const Execution::EventPlan& plan,std::size_t ordinal,
    InitialProductSnapshot* out,std::string* error) const {
  if(!out||ordinal>=reserved_.size()||ordinal>=plan.products().size())return Reject(error,"invalid staged product snapshot index");
  const auto ptr=reserved_[ordinal];const auto& product=plan.products()[ordinal];
  if(!PIC::ParticleBuffer::IsParticleAllocated(ptr)||PIC::ParticleBuffer::GetNext(ptr)!=-1||PIC::ParticleBuffer::GetPrev(ptr)!=-1)
    return Reject(error,"opaque initializer changed staged allocation/list ownership before publication");
  InitialProductSnapshot actual;actual.compiledSpecies=PIC::ParticleBuffer::GetI(ptr);
  if(actual.compiledSpecies<0||actual.compiledSpecies>=PIC::nTotalSpecies||actual.compiledSpecies!=product.compiledSpecies)
    return Reject(error,"opaque staged initializer corrupted compiled product species before publication");
  PIC::ParticleBuffer::GetX(actual.positionM.data(),ptr);PIC::ParticleBuffer::GetV(actual.velocityMS.data(),ptr);
  actual.referenceWeight=node_->block->GetLocalParticleWeight(actual.compiledSpecies);
  actual.individualCorrection=PIC::ParticleBuffer::GetIndividualStatWeightCorrection(ptr);
  if(!ReadReactionState(ptr,&actual.interval,&actual.clock,&actual.catalogSpecies,&actual.internalEnergyJ,error))return false;
  if(!PIC::Mover::ForceCutSurface::IsStorageRegistered())return Reject(error,"staged semantic physical products require registered physical epoch storage");
  actual.physicalEpochS=PIC::Mover::ForceCutSurface::GetParticleEpochS(ptr);
  if(actual.compiledSpecies!=product.compiledSpecies||actual.catalogSpecies!=product.species||
      actual.positionM!=product.positionM||actual.velocityMS!=product.velocityMS||
      actual.referenceWeight!=product.referenceWeight||actual.individualCorrection!=product.individualCorrection||
      actual.internalEnergyJ!=product.internalEnergyJ||actual.physicalEpochS!=product.interval.startS()+product.interval.elapsedS()||
      actual.interval.startS()!=product.interval.startS()||actual.interval.durationS()!=product.interval.durationS()||
      actual.interval.elapsedS()!=product.interval.elapsedS()||!actual.clock.SamePendingDepth(product.clock)||
      actual.clock.fresh()!=product.clock.fresh())
    return Reject(error,"opaque staged initializer altered immutable planned species/position/velocity/weight/state/clock/epoch. No product or event ledger was published; correct the opaque offset/initializer before restarting");
  *out=actual;return true;
}
bool ParticleTransaction::Commit(const Execution::EventPlan& plan,std::string* error) {
  if (!error || committed_ || !staged_ || stagedPlan_!=&plan || !ValidateParent(plan,error))
    return false;
  if (std::memcmp(parentSnapshot_.data(),PIC::ParticleBuffer::GetParticleDataPointer(parent_),parentSnapshot_.size())!=0)
    return Reject(error,"reaction parent bytes changed before atomic commit");
  // Re-read all initialized native values immediately before publication;
  // neither an opaque offset bug nor an intervening caller can corrupt the
  // physical packet while leaving a precomputed ledger looking conserved.
  if(initializeProduct_)for(std::size_t i=0;i<reserved_.size();++i)
    if(!ValidateInitializedProduct(plan,i,&stagedLedger_.initialProducts[i],error))return false;
  // Another staged transaction may have committed to this SAME ledger since
  // Reserve. Re-establish capacity while failure is still rollback-safe;
  // MPI-only execution excludes concurrent mutation inside the publication
  // region below, but does not exclude earlier nested/staged transactions.
  ledger_->events.reserve(ledger_->events.size()+1);
  // No allocations, callbacks, random draws or fallible physics after this
  // point. The MPI-only adapter owns list publication and the pre-reserved
  // ledger entry as one transaction. In-flight parent deletion remains the
  // mover's responsibility when retainParent was explicitly requested.
  if (!retainParent_) PIC::ParticleBuffer::DeleteParticle(parent_);
  if (placement_==ProductPlacement::RebuiltCellList) {
    for (long int ptr:reserved_) {
      PIC::ParticleBuffer::SetNext(head_,ptr);PIC::ParticleBuffer::SetPrev(-1,ptr);
      if (head_!=-1) PIC::ParticleBuffer::SetPrev(ptr,head_);
      head_=ptr;
    }
  }
  committedProducts_.swap(reserved_);
  ledger_->events.emplace_back(std::move(stagedLedger_));
  committed_=true;staged_=false;stagedPlan_=nullptr;parentSnapshot_.clear();
  error->clear();return true;
}

void ParticleTransaction::Rollback() noexcept {
  // Slots were never linked/published. Releasing them cannot alter the parent,
  // science ledger or pending event decisions. Allocator attempt counters may
  // advance; they are work diagnostics, not accepted physical reactions.
  for (long int ptr:reserved_)
    if (PIC::ParticleBuffer::IsParticleAllocated(ptr))
      PIC::ParticleBuffer::DeleteParticle_withoutTrajectoryTermination(ptr);
  reserved_.clear();parentSnapshot_.clear();staged_=false;stagedPlan_=nullptr;
}

} } // namespace Moon::Chemistry
