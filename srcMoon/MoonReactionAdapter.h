#ifndef SRCMOON_MOON_REACTION_ADAPTER_H
#define SRCMOON_MOON_REACTION_ADAPTER_H

#include "pic.h"
#include "species/ReactionExecution.h"

#include <string>
#include <vector>

namespace Moon { namespace Chemistry {
namespace Execution=Species::ReactionExecution;

// Extra bytes exist only for explicitly registered chemistry fixtures/modes.
// The two classes contain no pointers and are copied with memcpy, preserving
// the physical interval and optical depth on ordinary MPI migration. Fresh
// daughters receive their newly planned clock instead of copied parent state.
void RegisterReactionStorage();
bool ReactionStorageRegistered();
bool ReadReactionState(long int particle,Execution::TraceInterval* interval,
    Execution::OpticalDepthClock* clock,std::size_t* catalogSpecies,
    double* internalEnergyJ,std::string* error);

// Actual PIC setters are read back while products are still unlinked. This
// fixed numerical snapshot is NOT a duplicate/high-precision hypothetical
// product; independent verifiers reconstruct budgets from these native bytes.
struct InitialProductSnapshot {
  int compiledSpecies=-1;std::size_t catalogSpecies=0;
  Execution::Vector3 positionM{},velocityMS{};
  double referenceWeight=0,individualCorrection=0,internalEnergyJ=0,physicalEpochS=0;
  Execution::TraceInterval interval;Execution::OpticalDepthClock clock;
};
struct CommittedReaction {
  std::string identity;
  Execution::EventLedger ledger;
  std::vector<InitialProductSnapshot> initialProducts;
  double physicalEpochS=0;
  long double exposureApproximationErrorBound=0,clockStorageErrorBound=0;
};
struct ReactionLedger {
  // Append only accepted events. Reservation/deferral attempts are allocator
  // work, never physical reactions. Callers own MPI aggregation and output.
  std::vector<CommittedReaction> events;
};

enum class ProductPlacement { RebuiltCellList, DetachedContinuation };

// This callback runs only on a RESERVED, UNLINKED product, after its physical
// state/clock is assigned and before atomic publication. It may initialize that
// child's opaque application bytes only. It must not draw randomness, mutate
// parent bytes or append external ledgers. A false result rolls back the WHOLE
// packet, so resource/capacity checks cannot discard a physical daughter.
using StagedProductInitializer=bool(*)(const Execution::EventPlan&,
    std::size_t productOrdinal,long int parent,long int child,void*,std::string*);

// Real PIC allocator/list adapter for the common conservative planner.
// This class defines no rate law, recoil model or transport approximation.
// It accepts only an already validated Transport plan from an explicit caller
// policy; rate-only plans cannot reach a commit through PendingEvent.
//
// Preconditions: MPI-only execution; the old-list walker has saved the parent
// successor and parent is not in rebuiltCellHead. For an in-flight callback
// that retains ownership of the current parent, use retainParent=true: the
// parent is left allocated/unchanged for the mover's documented Delete result.
// Detached products stay UNLINKED until the caller queues each continuation.
class ParticleTransaction final : public Execution::TransactionBackend {
 public:
  ParticleTransaction(long int parent,long int& rebuiltCellHead,
      cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node,
      const std::vector<Execution::SpeciesIdentity>& catalog,
      ReactionLedger* ledger,ProductPlacement placement=ProductPlacement::RebuiltCellList,
      bool retainParent=false,double elementaryChargeC=ElectronCharge,
      StagedProductInitializer initializeProduct=nullptr,void* initializerUser=nullptr);
  ~ParticleTransaction() override;
  Execution::ReservationStatus Reserve(const Execution::EventPlan&,std::string* error) override;
  bool Stage(const Execution::EventPlan&,std::string* error) override;
  bool Commit(const Execution::EventPlan&,std::string* error) override;
  void Rollback() noexcept override;
  const std::vector<long int>& committedProducts() const { return committedProducts_; }

 private:
  bool ValidateParent(const Execution::EventPlan&,std::string* error) const;
  bool ValidateInitializedProduct(const Execution::EventPlan&,std::size_t,
      InitialProductSnapshot*,std::string* error) const;
  long int parent_;
  long int& head_;
  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node_;
  const std::vector<Execution::SpeciesIdentity>& catalog_;
  ReactionLedger* ledger_;
  ProductPlacement placement_;
  bool retainParent_,staged_=false,committed_=false;
  double elementaryChargeC_;
  StagedProductInitializer initializeProduct_;
  void* initializerUser_;
  const Execution::EventPlan* stagedPlan_=nullptr;
  std::vector<long int> reserved_,committedProducts_;
  std::vector<PIC::ParticleBuffer::byte> parentSnapshot_;
  CommittedReaction stagedLedger_;
};

} } // namespace Moon::Chemistry
#endif
