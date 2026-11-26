#ifndef SRCMOON_SODIUM_COHORT_H
#define SRCMOON_SODIUM_COHORT_H
#include "MoonTrajectoryExposure.h"
#include "MoonSemanticRandom.h"
#include "species/PhysicalReactionEvents.h"
#include "species/PhysicalStateProfile.h"
#include <array>
#include <memory>
namespace Moon { namespace Chemistry { namespace SodiumCohort {
namespace E=Species::PhysicalReactionEvents;
namespace P=Species::PhysicalStateProfile;
// This is an explicit SOURCE INPUT. No field chooses a missing physical
// isotope/state/threshold or widens an error limit. Omitted values are invalid.
struct Request {
  P::Files profile;
  std::string profileMetadataPath,profileMetadataSha256;
  // Deterministic native AMPS projection binds JSON identity and EXACT ordered
  // CSV/source request. It supplies no physical activation/readiness Boolean.
  std::string profileNativeManifestPath,profileNativeManifestSha256;
  std::string parentStateId,ionStateId,electronStateId;
  int parentCompiledSpecies=-1,ionCompiledSpecies=-1;
  // CSV exact header: channel_id,parent_state_id,ion_state_id,electron_state_id,
  // energy_reference_identity. One reviewed complete atomic source row only.
  // Explicit source-effective-photon input uses the documented eight-column
  // extension, retaining the native threshold and physical source identity.
  std::string crosswalkPath,crosswalkSha256;
  E::Definition controls; // all numeric/error/policy/source fields explicit
  SemanticRandom::Configuration random;
  Species::TrajectoryExposure::Controls exposure;
  std::size_t maximumRootRefinements=0;
  bool earthOccultation=false;
  std::string inertialFrame,inertialOrigin; // origin must explicitly be SUN
  double sourceReferenceDistanceM=0;
  std::string sourceReferenceDistanceIdentity;
  // VERIFICATION ONLY (I21 decision 1, 2026-11-26). When true, the frozen
  // photon field uses this static point-Sun position [m, frozen LSO] instead
  // of the SPICE Sun; earthOccultation must then be false. The identity names
  // the analytical fixture. False (default) keeps the normal SPICE geometry.
  bool prescribedSun=false;
  std::array<double,3> prescribedSunPositionM{};
  std::string prescribedSunIdentity;

};
struct ComposedEvent {
  SemanticRandom::State key;E::Receipt receipt;long double relativeEventTimeS=0;
  Species::ReactionExecution::ParentState parentAtEvent;
  Species::ReactionExecution::ProductPacket ionAtEvent;
  long double acceptedInstantaneousRateS=0;
  std::string acceptedExposureIdentity,rateConventionIdentity;
};
class Prepared : public std::enable_shared_from_this<Prepared> {
 public:
  // Data-first inspection does not require invented numerical/model-error
  // choices merely to report a genuinely missing source state. It authenticates
  // exact metadata projection/files and resolves ALL3states+physical masses.
  static E::Status InspectProfile(const Request&,P::Registry* output);
  static E::Status Create(const Species::Photochemistry::Registry&,const Request&,
      std::shared_ptr<Prepared>* output);
  // Register the optional normal amps_init boundary; no particles yet exist.
  E::Status RegisterForNormalStartup();
  const E::Definition& definition()const{return definition_;}
  const P::Registry& profile()const{return profile_;}
  const Species::TrajectoryExposure::Field& geometry()const{return field_;}
  const ReactionLedger& ledger()const{return ledger_;}
  const Species::IonCollector::Ledger& collectorLedger()const{return collectorLedger_;}
  const std::vector<ComposedEvent>& composedReceipts()const{return composed_;}
  // Called only for explicitly constructed source histories at their actual
  // initial epoch; ownership clones never call this reset/initializer.
  bool InitializeSourceKey(long int,std::uint64_t originalParent,std::string*);
 private:
  Request request_;P::Registry profile_;E::Definition definition_;
  E::SodiumPhotoionization provider_;ReactionLedger ledger_;
  std::shared_ptr<const Exposure::FrozenProductionField> geometry_;
  Species::TrajectoryExposure::Field field_;
  std::string exposureIdentity_; // controller exposure-provider identity (SPICE or prescribed Sun)
  // Preview receipts are not science events. These entries are composed only
  // for accepted Process; native output MUST pair them with committed ledger
  // identity/count before publishing them as accepted event evidence.
  std::vector<ComposedEvent> composed_;
  Species::IonCollector::Ledger collectorLedger_;
  double parentInternalEnergyJ_=0;
  double elementaryChargeC_=0;std::string elementaryChargeIdentity_;
  static bool Install(void*,std::string*);
  static bool CollectorState(long int,void*,void*,std::size_t*,std::string*);
  static bool Resolve(long int,Controller::Node*,void*,std::size_t*,double*,std::string*);
  static bool History(const Controller::ExposureRequest&,void*,Controller::ExposureHistory*,std::string*);
  static bool Kinematics(const Controller::EventContext&,Species::ReactionExecution::KinematicPolicy*,std::string*);
};
}}}
#endif
