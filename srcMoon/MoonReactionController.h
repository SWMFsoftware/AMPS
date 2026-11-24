#ifndef SRCMOON_MOON_REACTION_CONTROLLER_H
#define SRCMOON_MOON_REACTION_CONTROLLER_H

#include "MoonReactionAdapter.h"
#include "species/ReactionExposure.h"
#include "species/TrajectoryExposure.h"

#include <cstddef>
#include <cstdint>
#include <limits>
#include <string>
#include <vector>

namespace Moon { namespace Chemistry { namespace Controller {
namespace R=Species::ReactionExecution;
using Node=cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>;

// ConstantHistory declares an immutable hazard [s^-1] for a whole state
// history. TrajectoryDependent requires the explicit accepted-arc provider
// below; it must never borrow or silently freeze constantRateS instead.
enum class HazardScope { ConstantHistory, TrajectoryDependent };
enum class DrawPurpose { SourceClock, Channel, ProductClock, Kinematics };
using OpenUniformProvider=bool(*)(DrawPurpose,long int,std::size_t,
    void*,double*,std::string*);
using InitialStateResolver=bool(*)(long int,Node*,void*,std::size_t*,double*,std::string*);

// Read-only actual arc supplied to a qualified source/geometry provider. The
// acceleration comes through the SAME callback/macro used by the core tracer;
// providers must not substitute a straight chord or another force law. All
// coordinates are LSO SI and times are actual physical seconds/ET.
struct ExposureRequest {
  long int particle=-1;
  Node* node=nullptr;
  std::size_t parentCatalogSpecies=0;
  int compiledSpecies=-1;
  Species::TrajectoryExposure::Arc arc;
  PIC::Mover::ForceCutSurface::cTimeContext time;
  double physicalEpochS=0,maximumS=0;
};
struct ExposureHistory {
  std::vector<Species::PiecewiseHazard::Segment> segments;
  // Conservative absolute dimensionless exposure error over the COMPLETE
  // preview. Debiting this whole bound for an accepted prefix is safe; the
  // controller must not silently discard the provider's numerical budget.
  long double absoluteDepthErrorBound=std::numeric_limits<long double>::quiet_NaN();
};
using ExposureHistoryProvider=bool(*)(const ExposureRequest&,void*,ExposureHistory*,std::string*);

struct EventContext {
  long int particle=-1;
  Node* node=nullptr;
  std::size_t channelIndex=0;
  const std::vector<R::SpeciesIdentity>* catalog=nullptr;
  const R::Reaction* reaction=nullptr;
  R::ParentState parent;
  R::TraceInterval interval;
  double physicalEpochS=0;
  long double selectedInstantaneousRateS=0; // exact selected accepted LEFT-limit [s^-1]
  std::string acceptedExposureIdentity,rateConventionIdentity;
  long double exposureApproximationErrorBound=0,clockStorageErrorBound=0;
  // A qualified kinetic provider can request its own explicit product-clock
  // uniforms. Every returned draw must satisfy the same strict open interval.
  // None is requested by trajectory preview/refinement or allocation retries.
  OpenUniformProvider openUniform=nullptr;
  void* userData=nullptr;
};
using KinematicProvider=bool(*)(const EventContext&,R::KinematicPolicy*,std::string*);
struct Channel {
  std::size_t parentCatalogSpecies=0;
  R::Reaction reaction;
  double constantRateS=0;
  HazardScope scope=HazardScope::TrajectoryDependent; // explicit opt-in required
  std::string hazardEvidenceIdentity,kinematicProviderIdentity;
  KinematicProvider kinematics=nullptr;
};
struct Configuration {
  std::string identity; // caller's immutable source/policy/evidence identity
  std::string eventFrame; // explicit actual Moon event frame, currently "LSO"
  std::vector<R::SpeciesIdentity> catalog;
  std::vector<Channel> channels;
  InitialStateResolver resolveInitialState=nullptr;
  OpenUniformProvider openUniform=nullptr;
  void* userData=nullptr; // callback-only process-local context; never serialized
  ReactionLedger* ledger=nullptr;
  // Exact initialized PIC charge must agree with this explicit source-bound
  // constant; legacy default preserves existing mathematical fixtures.
  double elementaryChargeC=ElectronCharge;
  std::string elementaryChargeSourceIdentity="AMPS legacy ElectronCharge";
  StagedProductInitializer initializeProduct=nullptr;
  std::string productInitializerIdentity;
  void* productInitializerUser=nullptr;
  // Mandatory only for states with TrajectoryDependent channels. It must
  // provide the COMPLETE family over [0,maximumS], including competing
  // constant channels. Provider source/geometry/error controls are bound by
  // exposureProviderIdentity and the immutable run manifest, never inferred.
  ExposureHistoryProvider exposureHistory=nullptr;
  std::string exposureProviderIdentity;
  Species::PiecewiseHazard::NumericControls exposureControls{0};
  // Process configuration only. Keeping this option last preserves source
  // compatibility for callers which use positional aggregate initialization.
  R::OpticalDepthPrecision clockPrecision=R::OpticalDepthPrecision::Binary64Log;
};

// Installation is explicit, once, on every participating rank, before particle
// buffer allocation. It deep-copies immutable catalog/rates/topology and refuses
// legacy movers, generic electron execution and nonempty existing core hooks.
// No .in setting or physical-rate registry installs this controller implicitly.
bool Install(const Configuration& configuration,std::string* error);
bool Enabled();

// Normal Moon::MoveParticle calls this AFTER assigning a finite particle epoch
// and BEFORE selecting the first/second tracer. A truly new source uses the
// explicit state resolver and ONE clock draw. Surviving/MPI-cloned histories
// keep optical depth/channel and only open the caller's next physical interval.
void PrepareParticle(long int particle,double requestedPhysicalIntervalS,Node* node);
// A provider/source deliberately cloning a NEW physical history must reset its
// optional state. Ordinary ownership/migration clones must NEVER call this.
void ResetForNewParticle(PIC::ParticleBuffer::byte* data);

// Thin diagnostic reader, with no mutation or RNG. selectedChannel is meaningful
// only when due==true; the preserved clock is then consumed. This state is copied
// in ordinary opaque PIC particle bytes, never through pointer-based side tables.
struct History {
  R::TraceInterval interval;
  R::OpticalDepthClock clock;
  std::size_t catalogSpecies=0,selectedChannel=0;
  double internalEnergyJ=0;
  bool due=false;
  // Accepted selection, not a new point/ray evaluation. Persisted through
  // exact endpoint MPI handoff; zero until a genuine crossing is accepted.
  long double selectedInstantaneousRateS=0;
  std::uint64_t configurationFingerprint=0;
  // Persisted on MPI continuation and caller-interval reopening. These are
  // numerical diagnostics, not physical source or observational errors.
  long double exposureApproximationErrorBound=0,clockStorageErrorBound=0;
};
bool ReadHistory(long int particle,History* output,std::string* error);

// The actual core trajectory callbacks. Query is side-effect-free; Advance is
// called only for accepted physical flight; Process queues the COMPLETE packet
// and returns Delete while the mover retains sole ownership of parent deletion.
int Query(long int,int,Node*,double,double,double&);
void Advance(long int,int,Node*,double,double);
int Process(long int,double*,double*,int&,Node*,double,double&);

} } }
#endif
