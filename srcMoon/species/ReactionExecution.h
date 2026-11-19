#ifndef SRCMOON_SPECIES_REACTION_EXECUTION_H
#define SRCMOON_SPECIES_REACTION_EXECUTION_H

#include <array>
#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace Moon { namespace Species {
namespace PiecewiseHazard {class ClockUpdate;}
namespace ReactionExecution {

using Vector3=std::array<double,3>;

// Material composition is an integer identity, not an atomic mass inferred
// from a rounded AMPS mass table. Names can distinguish isotopes when the
// registry has authoritative isotope data (e.g. "H" versus "H:2"). A state
// remains explicit even when several states share one compiled PIC species.
struct NucleusCount { std::string identity; std::int64_t count=0; };
struct SpeciesIdentity {
  std::string identity,state;
  int compiledSpecies=-1; // -1 is permitted only for untracked bath terms.
  std::vector<NucleusCount> nuclei;
  std::int64_t elementaryCharge=0;
  double modelMassKg=0.0;
};
enum class Population { Tracked, Bath };
struct Term {
  std::size_t species=0;
  std::uint32_t multiplicity=0;
  Population population=Population::Tracked;
};
struct Reaction {
  std::string identity;
  std::vector<Term> reactants,products;
};

// All bool APIs leave output objects unchanged on failure and return a
// diagnostic in error. They neither terminate AMPS nor perform I/O or RNG.
bool ValidateCatalog(const std::vector<SpeciesIdentity>& catalog,std::string* error);
bool ValidateReaction(const std::vector<SpeciesIdentity>& catalog,
                      const Reaction& reaction,std::string* error);

// Immutable local-time interval. startS is a physical epoch on the caller's
// declared scale; duration and elapsed are seconds. Continuation retains the
// original start/duration rather than replacing them with a product timestep.
class TraceInterval {
 public:
  static bool Create(double startS,double durationS,double elapsedS,
                     TraceInterval* output,std::string* error);
  bool Advance(double seconds,TraceInterval* output,std::string* error) const;
  bool valid() const;
  double startS() const { return start_; }
  double durationS() const { return duration_; }
  double elapsedS() const { return elapsed_; }
  double remainingS() const { return duration_-elapsed_; }
 private:
  bool valid_=false;
  double start_=0,duration_=0,elapsed_=0;
};

// An optical depth is dimensionless. A consumed parent clock cannot be
// assigned to a daughter; fresh() distinguishes a new draw from continuation.
class OpticalDepthClock {
 public:
  static bool FromOpenUniform(double uniform,OpticalDepthClock* output,
                              std::string* error);
  bool valid() const;
  bool consumed() const { return consumed_; }
  bool fresh() const { return fresh_; }
  double remaining() const {
    return static_cast<double>(static_cast<long double>(remaining_)-compensation_);
  }
  // Exposure integrators must inspect the compensated pending depth without
  // first rounding it back to binary64. This is a pure dimensionless query;
  // it neither consumes exposure nor changes the clock's freshness.
  long double remainingExtended() const {
    return static_cast<long double>(remaining_)-compensation_;
  }
  // A compensated binary64 pair can contain more bits than one long double.
  // This legacy scalar query still refuses a lossy projection. Actual
  // variable-history location now uses PendingDepthComponents instead.
  bool ExactPendingDepth(long double* output,std::string* error) const;
  bool PendingDepthComponents(double* high,double* low,std::string* error)const;
  // Strong cache identity: BOTH serialized components and consumed state
  // agree. Equal rounded scalar displays are insufficient. No state/RNG is
  // changed, and invalid clocks never compare equal.
  bool SamePendingDepth(const OpticalDepthClock&)const;
 private:
  friend bool ConsumeExposure(const OpticalDepthClock&,double,OpticalDepthClock*,std::string*);
  friend bool ConstantHazardEventDelay(const std::vector<struct Hazard>&,
      const OpticalDepthClock&,double*,std::string*);
  friend bool ConstantHazardEventDelay(const std::vector<struct Hazard>&,
      const OpticalDepthClock&,long double*,std::string*);
  friend bool AdvanceConstantHazardClock(const std::vector<struct Hazard>&,
      long double,const OpticalDepthClock&,double,struct HazardClockAdvance*,std::string*);
  friend bool AdvanceIntegratedHazardClock(long double,const OpticalDepthClock&,
      const std::vector<struct Hazard>&,double,struct ExposureClockAdvance*,std::string*);
  friend bool ApplyLocatedHazardClock(const OpticalDepthClock&,
      const PiecewiseHazard::ClockUpdate&,const std::vector<struct Hazard>&,
      double,struct ExposureClockAdvance*,std::string*);
  friend bool AdvanceConstantHazards(const std::vector<struct Hazard>&,
      const TraceInterval&,const OpticalDepthClock&,double,struct HazardAdvance*,std::string*);
  bool valid_=false,consumed_=false,fresh_=false;
  // Compensated debit retains sub-ULP optical exposures across short tracer
  // segments; repeated tiny segments must not silently erase their hazard.
  double remaining_=0,compensation_=0;
  void SubtractExposure(long double exposure);
};
struct Hazard { std::size_t reaction=0; double rateS=0; };
// Conditional process probabilities at an already located event use the
// INSTANTANEOUS rates [s^-1], never integrated channel exposures. The common
// selector is shared by constant and nonconstant clock decisions. A zero
// total rate has no conditional channel and is refused transactionally.
bool SelectHazardChannel(const std::vector<Hazard>& instantaneousRates,
    double openUniform,std::size_t* reaction,std::string* error);
struct ExposureClockAdvance {
  bool event=false;
  std::size_t reaction=0;
  double totalInstantaneousRateS=0;
  OpticalDepthClock clock;
  long double residualDepthRoundingErrorBound=0;
};
// Accept a caller-integrated dimensionless TOTAL exposure from an immutable
// origin clock. An event must already be located on the actual supported
// trajectory by the caller; this service does not infer a crossing time.
// Only crossing decisions inspect instantaneous rates or the channel draw.
// Pure preview and retry can evaluate this function without mutating origin.
bool AdvanceIntegratedHazardClock(long double integratedExposure,
    const OpticalDepthClock& originClock,
    const std::vector<Hazard>& instantaneousRates,double channelUniform,
    ExposureClockAdvance* output,std::string* error);
// Publish an immutable exact-location receipt from the SAME origin pair.
// The exact locator owns crossing; a rounded displayed integral never
// reclassifies it. Instantaneous CDF semantics are unchanged. There is no
// new clock field/opaque particle ABI, source choice or numerical tolerance.
bool ApplyLocatedHazardClock(const OpticalDepthClock& originClock,
    const PiecewiseHazard::ClockUpdate& location,
    const std::vector<Hazard>& instantaneousRates,double channelUniform,
    ExposureClockAdvance* output,std::string* error);
// Pure preview of a CONSTANT combined hazard. Reads the compensated pending
// optical depth without selecting a channel or advancing state/RNG. The
// finite delay [s] is rounded upward only when needed: it is the first
// representable double not before exact H/sum(lambda). Zero total rate gives
// +infinity. Nonzero rates with unrepresentable positive delay are refused.
bool ConstantHazardEventDelay(const std::vector<Hazard>& hazards,
    const OpticalDepthClock& clock,double* delayS,std::string* error);
// Extended-precision origin-relative delay for a controller which retains an
// immutable interval origin. Its directed rounding uses the SAME exposure
// comparison as AdvanceConstantHazardClock; it performs no RNG or mutation.
bool ConstantHazardEventDelay(const std::vector<Hazard>& hazards,
    const OpticalDepthClock& clock,long double* delayS,std::string* error);
struct HazardClockAdvance {
  bool event=false;
  std::size_t reaction=0;
  double totalRateS=0;
  long double eventOffsetS=0;
  OpticalDepthClock clock;
};
// Reconstruct an immutable constant-rate clock from its interval-start value
// and TOTAL elapsed seconds, rather than summing rounded accepted substeps.
// The caller retains the physical interval separately. All outputs remain
// unchanged on error; a channel uniform is inspected only upon a crossing.
bool AdvanceConstantHazardClock(const std::vector<Hazard>& hazards,
    long double elapsedS,const OpticalDepthClock& originClock,
    double channelUniform,HazardClockAdvance* output,std::string* error);
struct HazardAdvance {
  bool event=false;
  std::size_t reaction=0;
  double totalRateS=0,eventOffsetS=0;
  TraceInterval interval;
  OpticalDepthClock clock;
};
// Constant-rate analytic event time and combined-channel selection. The
// supplied channel draw is inspected only if an event actually occurs.
// At an exact end-of-interval event, products have zero remaining time.
bool AdvanceConstantHazards(const std::vector<Hazard>& hazards,
    const TraceInterval& interval,const OpticalDepthClock& clock,
    double channelUniform,HazardAdvance* output,std::string* error);
bool EventProbability(const std::vector<Hazard>& hazards,double durationS,
                      double* output,std::string* error);
// Caller-integrated variable-rate exposure can be consumed only while it is
// strictly below the pending clock. Crossing requires the caller's actual
// trajectory/root finder; this function never guesses an event location.
bool ConsumeExposure(const OpticalDepthClock& clock,double integratedHazard,
                     OpticalDepthClock* output,std::string* error);

struct ParentState {
  std::size_t species=0;
  double referenceWeight=0,individualCorrection=0;
  Vector3 positionM{{0,0,0}},velocityMS{{0,0,0}};
  double internalEnergyJ=0; // per physical parent, approved state convention
};
struct ProductKinematics {
  std::size_t term=0;
  std::uint32_t occurrence=0; // one carrier for EACH repeated product
  double referenceWeight=0;
  Vector3 velocityMS{{0,0,0}};
  double internalEnergyJ=0;
  OpticalDepthClock nextClock;
};
struct ClosureBound {
  // Absolute SI units and dimensionless relative bound per physical event.
  double absolute=0,relative=0;
};
struct ConservationBounds {
  ClosureBound modelMassKg,momentumKgMS,energyJ;
};
// Signed CHANGE of the untracked environment per ONE physical parent event,
// so tracked delta + bath delta must close. Material bath counts/identities
// come from the reaction terms; these moments must come from the caller's
// qualified photon/electron/projectile/state model, never from residual fit.
struct BathExchange {
  bool supplied=false;
  double modelMassChangeKg=0,energyChangeJ=0;
  Vector3 momentumChangeKgMS{{0,0,0}};
};
struct KinematicPolicy {
  bool qualified=false;
  std::string identity,evidenceIdentity,eventFrame;
  double maximumClassicalSpeedMS=0; // explicit validated nonrelativistic domain
  ConservationBounds bounds;
  BathExchange bath;
  std::vector<ProductKinematics> products; // tracked product occurrences only
};
enum class PlanScope { RateOnly, Transport };
struct ProductPacket {
  // species is the registry/state catalog index, NOT a PIC species index.
  // Resolve catalog[species].compiledSpecies and persist the state identity.
  std::size_t species=0,term=0;
  // Snapshot the resolved identities too: a later accidental catalog reorder
  // must not silently change a pending packet's PIC species/internal state.
  int compiledSpecies=-1;
  std::string speciesIdentity,stateIdentity;
  std::uint32_t occurrence=0;
  double targetRepresentedNumber=0,representedNumber=0;
  double referenceWeight=0,individualCorrection=0;
  Vector3 positionM{{0,0,0}},velocityMS{{0,0,0}};
  double internalEnergyJ=0;
  TraceInterval interval;
  OpticalDepthClock clock;
};
struct BathTermRecord {
  std::size_t species=0,term=0;
  bool reactant=false;
  double representedNumber=0;
};
struct EventLedger {
  // Exact integer composition per ONE physical event, before statistical
  // weights. Noninteger packet weights do not make stoichiometry fractional.
  bool integerConservationVerified=false;
  bool kinematicConservationVerified=false; // false for an explicit rate-only plan
  std::vector<NucleusCount> reactantNuclei,productNuclei;
  std::int64_t reactantChargeElementary=0,productChargeElementary=0;
  std::int64_t trackedChargeChangeElementary=0,bathChargeChangeElementary=0;
  double parentNumber=0,parentModelMassKg=0;
  double productNumber=0,productModelMassKg=0;
  Vector3 trackedMomentumChangeKgMS{{0,0,0}},bathMomentumChangeKgMS{{0,0,0}};
  double trackedEnergyChangeJ=0,bathEnergyChangeJ=0,bathModelMassChangeKg=0;
  std::vector<BathTermRecord> bathTerms;
};

// Published plans are read-only. A plan owns every physical decision needed
// by a capacity retry; no reaction/channel/kinematic draw is repeated there.
class EventPlan {
 public:
  bool valid() const { return valid_; }
  PlanScope scope() const { return scope_; }
  const std::string& reactionIdentity() const { return reactionIdentity_; }
  const std::string& kinematicPolicyIdentity() const { return policyIdentity_; }
  const std::string& eventFrame() const { return frame_; }
  const ParentState& parent() const { return parent_; }
  int parentCompiledSpecies() const { return parentCompiledSpecies_; }
  const std::string& parentSpeciesIdentity() const { return parentIdentity_; }
  const std::string& parentStateIdentity() const { return parentState_; }
  const TraceInterval& interval() const { return interval_; }
  const std::vector<ProductPacket>& products() const { return products_; }
  const EventLedger& ledger() const { return ledger_; }
 private:
  friend bool MakeEventPlan(const std::vector<SpeciesIdentity>&,const Reaction&,
      const ParentState&,const TraceInterval&,PlanScope,const KinematicPolicy*,
      EventPlan*,std::string*);
  bool valid_=false;
  PlanScope scope_=PlanScope::RateOnly;
  std::string reactionIdentity_,policyIdentity_,frame_;
  std::string parentIdentity_,parentState_;
  int parentCompiledSpecies_=-1;
  ParentState parent_;
  TraceInterval interval_;
  std::vector<ProductPacket> products_;
  EventLedger ledger_;
};
// eventInterval is the caller-located physical event interval; positions are
// taken from the parent snapshot AT that event, not extrapolated here.
// RateOnly validates/computes material bookkeeping but cannot commit particles.
bool MakeEventPlan(const std::vector<SpeciesIdentity>& catalog,
    const Reaction& reaction,const ParentState& parent,
    const TraceInterval& eventInterval,PlanScope scope,
    const KinematicPolicy* policy,EventPlan* output,std::string* error);

enum class ReservationStatus { Ready, Deferred, Error };
class TransactionBackend {
 public:
  virtual ~TransactionBackend()=default;
  // Reserve/stage must not publish products or alter parent/physical ledgers.
  // Rollback must release ALL private state, including partial reservations.
  virtual ReservationStatus Reserve(const EventPlan&,std::string* error)=0;
  virtual bool Stage(const EventPlan&,std::string* error)=0;
  // Commit is atomic with respect to parent/products/physical diagnostics.
  // Failure/exception must remain rollback-safe; the PIC adapter owns this
  // guarantee and must verify actual lists/slots before enabling production.
  virtual bool Commit(const EventPlan&,std::string* error)=0;
  virtual void Rollback() noexcept=0;
};
enum class AttemptStatus { Committed, Deferred, Error };
class PendingEvent {
 public:
  explicit PendingEvent(const EventPlan& plan):plan_(plan) {}
  const EventPlan& plan() const { return plan_; }
  bool committed() const { return committed_; }
  AttemptStatus TryCommit(TransactionBackend* backend,std::string* error);
 private:
  const EventPlan plan_;
  bool committed_=false;
};

} } } // namespace Moon::Species::ReactionExecution
#endif
