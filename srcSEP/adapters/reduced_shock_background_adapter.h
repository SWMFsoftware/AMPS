#ifndef SRCSEP_ADAPTERS_REDUCED_SHOCK_BACKGROUND_ADAPTER_H
#define SRCSEP_ADAPTERS_REDUCED_SHOCK_BACKGROUND_ADAPTER_H

#include <array>
#include <cstdint>
#include <iosfwd>
#include <memory>
#include <string>
#include "../util/sep_common_header_path.h"
#include SRCSEP_SEP_COMMON_HEADER(sep_status.h)
#include SRCSEP_SEP_COMMON_HEADER(sep_material_background.h)

namespace SEP { namespace CoronaSwcme { namespace ShockFront {
class Provider;
struct Epoch;
} } }

namespace SEP { namespace ReducedShock {

using Position=std::array<double,3>; // HCI Cartesian position [m].

// AMPS-independent SI sample returned by the shared reduced provider.  This
// type is intentionally small: srcSEP owns the conversion into field-line
// vertex storage, while all coronal/Parker plasma and IMF equations remain in
// src/models/sep_corona_swcme.  ``numberDensityM3`` is electron number density
// and ``pressurePa`` is the provider's total thermal pressure.
struct AmbientSample {
  std::array<double,3> magneticFieldT{{0.0,0.0,0.0}};
  std::array<double,3> velocityMPerS{{0.0,0.0,0.0}};
  double numberDensityM3 = 0.0;
  // Preserve the old electron-density field above for source compatibility.
  // Particle sources and native transport use the explicitly named proton
  // density; ne=np is true only for an explicitly pure-H composition.
  double electronNumberDensityM3 = 0.0;
  double protonNumberDensityM3 = 0.0;
  double alphaNumberDensityM3 = 0.0;
  double pressurePa = 0.0;
  double protonTemperatureK = 0.0;
  double massDensityKgM3 = 0.0;
  // Cartesian row-major G_ij=partial U_i/partial x_j [1/s], from the
  // provider's branch-preserving derivative evaluator.
  std::array<double,9> velocityGradientPerS{};
  std::array<double,9> magneticGradientTPerM{};
  std::array<double,3> logMagneticGradientPerM{};
  // Core stencil/region identifiers are diagnostic provenance, not choices
  // for an application-owned derivative or plasma closure.
  std::array<int,3> derivativeStencils{};
  int ambientRegion=0,magneticSector=0;
  double epochS=0;
  std::uint64_t generation=0;
  std::string eventIdentity;
  // Material geometry is qualified for the core Parker exterior. An open PFSS
  // tube is qualified only when the event selects the v1.3 isothermal
  // flux-tube closure (steady corotating field-aligned wind, continuity per
  // unit flux, isothermal Bernoulli); otherwise PFSS branches remain
  // queryable but do not inherit an unproved induction/continuity contract.
  // Closed PFSS tubes are never material states.
  bool parkerExterior=false;
  bool coronalOpen=false;
  bool materialQualified=false;
  bool derivativesValid = false;
};

struct SurfaceSphere {
  Position centerM{},direction{};
  double radiusM=0,minimumRayCosine=0;
};
struct FrontSample {
  Position positionM{},outwardNormal{};
  double normalSpeedMPerS=0,inflowMPerS=0,fastMach=0,compressionRatio=0;
  AmbientSample upstream;
  bool upstreamValid=false,acceptedShock=false,downstreamValid=false,supportEdge=false;
  Position downstreamMagneticFieldT{},downstreamVelocityMPerS{};
  double downstreamMassDensityKgM3=0,downstreamPressurePa=0;
  std::string classification,eventIdentity;
  double epochS=0;
  std::uint64_t generation=0;
};

// A read-only physical realization over a CLOSED interval [departure,arrival].
// It retains the actual provider, so later configuration/publication cannot
// rebind an in-flight geometry or particle query to a different event. The
// provider's prescribed functions are queried at explicit time, without
// mutating Current() or calling Prepare at an artificial cadence tick.
class IntervalView {
 public:
  double DepartureS() const {return departureS_;}
  double ArrivalS() const {return arrivalS_;}
  std::uint64_t Generation() const {return generation_;}
  const std::string& EventIdentity() const {return identity_;}
  Core::Result<AmbientSample> AmbientAt(const Position&,double epochS) const;
  // Trajectory and B tracing need primitives but no derivative stencil.
  // Avoid evaluating unrelated off-line stencil points during these queries.
  Core::Result<AmbientSample> PrimitiveAt(const Position&,double epochS) const;
  Core::Status RequirePhysicalDownstream() const;
  Core::Result<SurfaceSphere> SphereAt(double epochS) const;
  // Conservative bound on the core plasma speed [m/s] over the whole support
  // (AmbientModel::SpeedBoundMPerS); used only to exclude geometry.
  double SpeedBoundMPerS() const;
  Core::Result<Position> LeadingPointOnRay(const Position&,double epochS) const;
  Core::Result<FrontSample> FrontAt(const Position&,double epochS) const;
  // Stage-11 neutral export from this exact captured core owner. These APIs
  // expose existing SI physics without adding fields to AmbientSample or
  // changing that established record's layout/meaning. Geometry and native
  // snapshot generations are separate identities; neither is cache cadence.
  Core::Result<Background::Material::Descriptor> MaterialDescriptor(
      std::uint64_t geometryGeneration,std::uint64_t snapshotGeneration,
      const std::string& nativeConfigurationIdentity) const;
  Core::Result<Background::Material::State> MaterialStateAt(
      const Background::Material::Descriptor&,const Background::Material::Query&) const;
 private:
  friend class CadencePublication;
  friend Core::Result<IntervalView> CaptureInterval(double,double,const std::string&);
  friend Core::Result<IntervalView> CapturePhysicalStep(double,double,const std::string&);
  Core::Status ValidateEpoch(double) const;
  std::shared_ptr<const CoronaSwcme::ShockFront::Provider> provider_;
  double departureS_=0,arrivalS_=0;
  std::uint64_t generation_=0;
  std::string identity_;
};
Core::Result<IntervalView> CaptureInterval(double departureS,double arrivalS,
                                         const std::string& expectedEventIdentity);
// The prescribed core has explicit continuous-time evaluators throughout its
// declared event support. Capture those functions for one full physical step,
// including a cadence crossing, without interpolating prepared arrays or
// changing the committed front. Committed-cadence queries above stay strict.
Core::Result<IntervalView> CapturePhysicalStep(double departureS,double arrivalS,
                                             const std::string& expectedEventIdentity);
bool EnsurePreparedForEpoch(double epochS,std::string* error);
double EventEndS();

// Metadata that binds native field-line storage to one immutable provider
// epoch.  Generation follows the event's declared background cadence rather
// than the number of calls made by the application.
struct EpochMetadata {
  double epochS = 0.0;
  double validUntilS = 0.0;
  std::uint64_t generation = 0;
  std::string eventIdentity;
};

// Prepared cadence-cache ownership token. A native transaction can build the
// complete front candidate without changing the provider seen by existing
// readers. Commit exchanges an already allocated owner, preserving retained
// IntervalViews and leaving every live cache untouched on preparation failure.
class CadencePublication {
 public:
  std::uint64_t Generation() const {return metadata_.generation;}
  // Bind a full physical query interval to the STAGED cache owner. Checkpoint
  // preparation must not obtain a view from the receiver's still-live cache
  // and subsequently relabel it as the saved realization after commit.
  Core::Result<IntervalView> PhysicalInterval(double departureS,double arrivalS) const;
 private:
  friend Core::Result<CadencePublication> PrepareCadencePublication(double,bool);
  friend Core::Status CommitCadencePublication(CadencePublication*);
  std::shared_ptr<CoronaSwcme::ShockFront::Provider> expected_,candidate_;
  std::shared_ptr<const CoronaSwcme::ShockFront::Epoch> expectedEpoch_;
  EpochMetadata metadata_;
};
// The explicit restore flag is reserved for a validated collective checkpoint
// transaction. Normal cadence publication remains monotonic; restoring a saved
// realization may reconstruct an earlier immutable cache from the same event.
Core::Result<CadencePublication> PrepareCadencePublication(double physicalEpochS,bool checkpointRestore=false);
Core::Status CommitCadencePublication(CadencePublication*);

struct FrontSummary {
  double apexRadiusM = 0.0;
  double acceptedShockAreaM2 = 0.0;
  bool apexShockAccepted = false;
};

// Configure the process-wide srcSEP adapter from one event file and its
// repository-relative assets.  Configuration is transactional: a failed
// candidate does not replace an already configured provider.  The selected
// event must be HCI and particle_mode=disabled: the CORE supplies prescribed
// ambient/front functions and does not run particles. A separately fingerprinted
// application-owned particle consumer may use those functions; it neither
// enables particle transport in the core nor supplies a downstream volume.
bool Configure(const std::string& eventPath, std::string* error);
bool Enabled();

// Prepare a complete front/ambient epoch and commit it only after all shared
// geometry, ambient, and RH checks succeed.  A repeated call at the same epoch
// is idempotent; time reversal and two different epochs mapping to one cadence
// generation are rejected because either would make snapshot provenance
// ambiguous.
bool Prepare(double epochS, std::string* error);

// Query the undisturbed ambient authority at an HCI Cartesian position.  The
// reduced model deliberately has no spatial downstream CME volume, so this
// routine never paints a local RH jump behind the front into the field line.
bool EvaluateAmbient(const std::array<double,3>& positionM,
                     AmbientSample* sample, std::string* error);

const EpochMetadata& Metadata();
double BackgroundCadenceS();
const std::string& EventPath();
const std::string& EventIdentity();
std::shared_ptr<const CoronaSwcme::ShockFront::Epoch> FrontEpoch();
FrontSummary CurrentFrontSummary();
// True when the bound event selects the closed heliocentric-sphere geometry
// (schema v1.3) rather than the finite SSE cap. False when no event is bound.
bool HeliocentricSphereFront();
// Write the active fitted apex at a completed physical epoch. Analytical
// configuration and cached-cadence metadata remain immutable.
Core::Status ReportFittedFront(double epochS,std::ostream* writer);

} } // namespace SEP::ReducedShock

#endif
