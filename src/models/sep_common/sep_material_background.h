#ifndef SEP_COMMON_MATERIAL_BACKGROUND_H
#define SEP_COMMON_MATERIAL_BACKGROUND_H

#include "sep_background_snapshot.h"
#include "sep_status.h"
#include <array>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <limits>
#include <locale>
#include <memory>
#include <optional>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

namespace SEP { namespace Background { namespace Material {
using Vector = std::array<double,3>; // Inertial HCI Cartesian components; SI.
using Tensor = std::array<double,9>; // Row-major partial component_i/partial x_j.

// Mandatory numerical records start INVALID until their physical owner fills
// every component. Zero U, zero B or a zero derivative can be an intentional
// supplied value, so zero-initialization must not impersonate missing data.
// These NaNs are construction guards, never physical replacements/absence:
// validation rejects them. Optional physics continues to use std::optional.
constexpr double Uninitialized=std::numeric_limits<double>::quiet_NaN();
constexpr Vector UninitializedVector{{Uninitialized,Uninitialized,Uninitialized}};
constexpr Tensor UninitializedTensor{{Uninitialized,Uninitialized,Uninitialized,
    Uninitialized,Uninitialized,Uninitialized,Uninitialized,Uninitialized,Uninitialized}};

// These IDs identify audited PROFILES, not a user's assertion of capability.
// Values are serialized explicitly below; append future IDs, never reinterpret
// an existing legacy diagnostic profile as a new conservative material law.
enum class Profile : std::uint32_t {
  ReducedParkerV1=1, Swcme1dShockOnlyV1, Swcme1dFullIcmeDiagnosticV1,
  Swcme3dSphereShockOnlyV1, Swcme3dEllipsoidShockOnlyV1, Swcme3dSseShockOnlyV1,
  Swcme3dSphereFullIcmeDiagnosticV1, Swcme3dEllipsoidFullIcmeDiagnosticV1,
  Swcme3dSseFullIcmeDiagnosticV1, CoronalBundleSnapshotV1,
  CoronalEllipsoidLivePlannedV1, FullCoronaCrossRayPlannedV1
};
enum class Use { SnapshotExchange, MaterialTransport, PhysicalDownstream };
enum class Region { Unspecified, ParkerExterior, CoronalOpen, CoronalClosed,
                    Sheath, Ejecta, Transition, PostEvent };
enum class Side { Smooth, UpstreamLimit, DownstreamLimit };
enum class DensityBasis { Unspecified, CoreNonrelativisticFluid,
                          SpeciesProper, DeclaredFrame };
enum class Capability : std::uint64_t {
  Primitives=1, Species=2, Thermodynamics=4, VelocityGradient=8,
  MagneticGradient=16, DenseTime=32, Trajectory=64, InverseMap=128,
  Deformation=256, ShockSurface=512, ContactSurface=1024, DownstreamVolume=2048
};
constexpr std::uint64_t Bit(Capability capability) {
  return static_cast<std::uint64_t>(capability);
}
// Masks are serialized values; unknown bits reject rather than activating
// a future or provider-specific quantity in an older consumer.
constexpr std::uint64_t AllCapabilities=4095;
constexpr std::uint64_t MaterialRequirements=1|2|4|8|16|32|64;
struct ProfileInfo {
  Profile profile;
  const char* key;
  std::uint64_t possible;
  bool materialQualified;
  const char* missing;
};
inline Core::Result<ProfileInfo> Lookup(Profile profile) {
  // "possible" bounds what an actual owner may expose; it is not proof that
  // an optional field exists at a particular query or that physics is closed.
  // Legacy diagnostic evaluators remain callable, but never grant material use.
  const ProfileInfo profiles[]={
    {Profile::ReducedParkerV1,"reduced-sse-parker-exterior-v1",MaterialRequirements|512,true,""},
    {Profile::Swcme1dShockOnlyV1,"swcme1d-shock-only-legacy-v1",1|2|4|8|512,false,"SWCME1D moving SEP transport is unsupported: Leblanc-density/constant-wind material continuity is unclosed"},
    {Profile::Swcme1dFullIcmeDiagnosticV1,"swcme1d-full-icme-diagnostic-v1",1|2|8|512,false,"SWCME1D moving SEP transport is unsupported: material continuity and regional thermodynamic export are unqualified"},
    {Profile::Swcme3dSphereShockOnlyV1,"swcme3d-sphere-shock-only-legacy-v1",1|2|4|8|16|512,false,"SWCME3D sphere material evolution is unqualified"},
    {Profile::Swcme3dEllipsoidShockOnlyV1,"swcme3d-ellipsoid-shock-only-legacy-v1",1|2|4|8|16|512,false,"SWCME3D ellipsoid material evolution is unqualified"},
    {Profile::Swcme3dSseShockOnlyV1,"swcme3d-sse-shock-only-legacy-v1",1|2|4|8|16|512,false,"SWCME3D finite-SSE material evolution is unqualified"},
    {Profile::Swcme3dSphereFullIcmeDiagnosticV1,"swcme3d-sphere-full-icme-diagnostic-v1",1|2|4|8|16|512,false,"SWCME3D regional material/thermal closure is unqualified"},
    {Profile::Swcme3dEllipsoidFullIcmeDiagnosticV1,"swcme3d-ellipsoid-full-icme-diagnostic-v1",1|2|4|8|16|512,false,"SWCME3D regional material/thermal closure is unqualified"},
    {Profile::Swcme3dSseFullIcmeDiagnosticV1,"swcme3d-sse-full-icme-diagnostic-v1",1|2|4|8|16|512,false,"SWCME3D regional material/thermal closure is unqualified"},
    {Profile::CoronalBundleSnapshotV1,"coronal-bundle-snapshot-v1",1|4,false,"Coronal-CME moving SEP transport is unsupported: schema-1 bundle has no material lineage, temporal primitives or full gradient"},
    {Profile::CoronalEllipsoidLivePlannedV1,"coronal-ellipsoid-live-planned-v1",0,false,"Coronal live material consumer awaits qualified G1/G2 background and architecture reconciliation"},
    {Profile::FullCoronaCrossRayPlannedV1,"full-corona-cross-ray-planned-v1",0,false,"Full corona-SWCME moving SEP transport is unsupported: BG3D-4/5/6 closure and native lifecycle remain unqualified"}
  };
  for(const auto& item:profiles)if(item.profile==profile)
    return Core::Result<ProfileInfo>::Success(item);
  return Core::Result<ProfileInfo>::Failure(Core::StatusCode::InvalidConfiguration,
      "unknown material background profile; no family/default fallback");
}
inline Core::Status Require(Profile profile,Use use,Region region=Region::Unspecified,
                            std::uint32_t revision=1) {
  const auto info=Lookup(profile);
  if(!info.ok())return info.status;
  if(revision!=1)return Core::Status::Failure(Core::StatusCode::UnsupportedCapability,
      "unsupported material profile revision");
  if(use==Use::SnapshotExchange) {
    if(!(info.value.possible&Bit(Capability::Primitives)))
      return Core::Status::Failure(Core::StatusCode::UnsupportedCapability,info.value.missing);
    return Core::Status::Success();
  }
  if(use!=Use::MaterialTransport&&use!=Use::PhysicalDownstream)
    return Core::Status::Failure(Core::StatusCode::InvalidConfiguration,"unknown material background operation");
  if(!info.value.materialQualified)
    return Core::Status::Failure(Core::StatusCode::UnsupportedCapability,info.value.missing);
  // The existing qualification is the Parker EXTERIOR ambient and upstream
  // source/absorbing-return model. Its local RH limits never imply a volume.
  if(use==Use::PhysicalDownstream||region!=Region::ParkerExterior)
    return Core::Status::Failure(Core::StatusCode::UnsupportedCapability,
        "profile has no qualified material evolution in the requested region/downstream volume");
  return Core::Status::Success();
}

// Descriptor is serializable identity/coverage only. Deserializing it cannot
// reconstruct a physical model or admit transport: a real immutable Owner and
// operation-specific admission are mandatory. No particle timestep is stored
// as an evolving background cadence, and no owner address enters the codec.
struct Descriptor {
  Profile profile=static_cast<Profile>(0);
  std::uint32_t revision=0;
  std::uint64_t capabilities=0,backgroundGeneration=0,geometryGeneration=0,snapshotGeneration=0;
  // PIC's initial REGISTERED geometry generation is zero. Presence must be
  // explicit; treating zero as missing either rejects a real initial line or
  // encourages renumbering physical/restart/sampling history to satisfy a
  // metadata check. Never infer registration from the numeric counter alone.
  bool geometryGenerationPresent=false;
  std::string eventIdentity,physicalConfigurationIdentity,nativeConfigurationIdentity;
  std::string materialIdentity,compositionIdentity,eosIdentity,frame,units;
  // Seconds from the canonical event origin, including both endpoints.
  // This is physical query coverage, not cache cadence or particle dt.
  double beginS=std::numeric_limits<double>::quiet_NaN();
  double endS=std::numeric_limits<double>::quiet_NaN();
  DensityBasis densityBasis=DensityBasis::Unspecified;
};
inline Core::Status Validate(const Descriptor& value) {
  const auto profile=Lookup(value.profile);if(!profile.ok())return profile.status;
  if(value.revision!=1||value.frame!="HCI"||value.units!="SI"||
      value.eventIdentity.empty()||value.physicalConfigurationIdentity.empty()||
      value.backgroundGeneration==0||!std::isfinite(value.beginS)||
      !std::isfinite(value.endS)||value.endS<value.beginS)
    return Core::Status::Failure(Core::StatusCode::InvalidConfiguration,"incomplete material identity, frame, units, revision or coverage");
  if((value.capabilities&~AllCapabilities)||(value.capabilities&~profile.value.possible))
    return Core::Status::Failure(Core::StatusCode::UnsupportedCapability,"owner claims an unknown/unqualified profile capability");
  if(value.densityBasis!=DensityBasis::CoreNonrelativisticFluid&&
      value.densityBasis!=DensityBasis::SpeciesProper&&value.densityBasis!=DensityBasis::DeclaredFrame)
    return Core::Status::Failure(Core::StatusCode::InvalidConfiguration,"density basis is absent/unknown; no proper/lab conversion is inferred");
  if((value.capabilities&Bit(Capability::Species))&&value.compositionIdentity.empty())
    return Core::Status::Failure(Core::StatusCode::InvalidConfiguration,"species capability has no composition identity");
  if((value.capabilities&Bit(Capability::Thermodynamics))&&value.eosIdentity.empty())
    return Core::Status::Failure(Core::StatusCode::InvalidConfiguration,"thermal capability has no EOS identity");
  if((value.capabilities&Bit(Capability::Trajectory))&&
      (value.materialIdentity.empty()||!value.geometryGenerationPresent))
    return Core::Status::Failure(Core::StatusCode::InvalidConfiguration,"material trajectory has no lineage/geometry generation");
  if(!value.nativeConfigurationIdentity.empty()&&value.snapshotGeneration==0)
    return Core::Status::Failure(Core::StatusCode::InvalidConfiguration,"native attachment has no explicit snapshot generation");
  return Core::Status::Success();
}
inline std::string Payload(const Descriptor& d) {
  std::ostringstream out;out.imbue(std::locale::classic());
  out<<std::setprecision(std::numeric_limits<double>::max_digits10);
  out<<"SEP_MATERIAL_DESCRIPTOR 1\n"<<static_cast<std::uint32_t>(d.profile)<<' '<<d.revision<<' '
     <<d.capabilities<<' '<<d.backgroundGeneration<<' '<<d.geometryGeneration<<' '<<d.geometryGenerationPresent<<' '<<d.snapshotGeneration<<' '
     <<d.beginS<<' '<<d.endS<<' '<<static_cast<unsigned>(d.densityBasis)<<'\n';
  for(const auto* s:{&d.eventIdentity,&d.physicalConfigurationIdentity,&d.nativeConfigurationIdentity,
      &d.materialIdentity,&d.compositionIdentity,&d.eosIdentity,&d.frame,&d.units})out<<std::quoted(*s)<<'\n';
  return out.str();
}
inline Core::Result<std::string> Encode(const Descriptor& d) {
  const auto valid=Validate(d);
  if(!valid.ok())return Core::Result<std::string>::Failure(valid.code,valid.message);
  const auto bytes=Payload(d);
  // The existing deterministic identity hash detects accidental corruption;
  // it is explicitly not an authentication or observational-provenance claim.
  return Core::Result<std::string>::Success(bytes+"digest "+FingerprintConfiguration(bytes)+"\n");
}
inline Core::Result<Descriptor> Decode(const std::string& bytes) {
  using R=Core::Result<Descriptor>;
  std::istringstream in(bytes);in.imbue(std::locale::classic());
  std::string tag;unsigned schema=0,profile=0,basis=0;Descriptor d;
  if(!(in>>tag>>schema)||tag!="SEP_MATERIAL_DESCRIPTOR"||schema!=1)
    return R::Failure(Core::StatusCode::UnsupportedCapability,"unsupported material descriptor schema; legacy bundle/checkpoint bytes are not material state");
  if(!(in>>profile>>d.revision>>d.capabilities>>d.backgroundGeneration>>d.geometryGeneration>>d.geometryGenerationPresent>>d.snapshotGeneration
      >>d.beginS>>d.endS>>basis))
    return R::Failure(Core::StatusCode::DataIntegrityFailure,"truncated material descriptor");
  d.profile=static_cast<Profile>(profile);d.densityBasis=static_cast<DensityBasis>(basis);
  for(auto* s:{&d.eventIdentity,&d.physicalConfigurationIdentity,&d.nativeConfigurationIdentity,
      &d.materialIdentity,&d.compositionIdentity,&d.eosIdentity,&d.frame,&d.units})
    if(!(in>>std::quoted(*s)))return R::Failure(Core::StatusCode::DataIntegrityFailure,"truncated material identity");
  std::string digest;
  if(!(in>>tag>>digest)||tag!="digest"||digest!=FingerprintConfiguration(Payload(d)))
    return R::Failure(Core::StatusCode::DataIntegrityFailure,"material descriptor digest mismatch");
  in>>std::ws;if(!in.eof())return R::Failure(Core::StatusCode::DataIntegrityFailure,"trailing material descriptor data");
  const auto valid=Validate(d);if(!valid.ok())return R::Failure(valid.code,valid.message);
  return R::Success(std::move(d));
}

struct Query {
  Vector positionM=UninitializedVector;
  double epochS=std::numeric_limits<double>::quiet_NaN();
  Region region=Region::Unspecified;
  Side side=Side::Smooth;
  std::string interfaceIdentity;
  std::uint64_t required=Bit(Capability::Primitives);
};
struct Species {
  std::string identity; // Explicit species, never a native AMPS index.
  double numberDensityM3=std::numeric_limits<double>::quiet_NaN();
  std::optional<double> temperatureK; // Absence is not a cold/zero population.
};
struct Gradient {
  Tensor values=UninitializedTensor;
  std::string branchIdentity;
  // A true declaration applies to analytic derivatives as well as actual
  // one-sided stencils. It never licenses smoothing across a discontinuity.
  bool branchPreserved=false;
};
struct State {
  std::string eventIdentity,physicalConfigurationIdentity;
  std::uint64_t backgroundGeneration=0,geometryGeneration=0;
  double epochS=std::numeric_limits<double>::quiet_NaN();
  Vector positionM=UninitializedVector,velocityMPerS=UninitializedVector,magneticFieldT=UninitializedVector;
  Region region=Region::Unspecified;Side side=Side::Smooth;
  std::string interfaceIdentity;
  std::optional<double> massDensityKgM3,pressurePa,gammaAdiabatic;
  std::optional<std::vector<Species>> species;
  std::optional<Gradient> velocityGradientPerS,magneticGradientTPerM;
};
inline bool Finite(const Vector& v) {for(double x:v)if(!std::isfinite(x))return false;return true;}
inline Core::Status Validate(const Descriptor& d,const Query& q,const State& s) {
  if(s.eventIdentity!=d.eventIdentity||s.physicalConfigurationIdentity!=d.physicalConfigurationIdentity||
      s.backgroundGeneration!=d.backgroundGeneration||s.geometryGeneration!=d.geometryGeneration||
      s.epochS!=q.epochS||s.positionM!=q.positionM||s.region!=q.region||s.side!=q.side||
      s.interfaceIdentity!=q.interfaceIdentity)
    return Core::Status::Failure(Core::StatusCode::DataIntegrityFailure,"material sample belongs to a foreign epoch, region, side or owning realization");
  if(!Finite(s.positionM)||!Finite(s.velocityMPerS)||!Finite(s.magneticFieldT)||
      !s.massDensityKgM3||!std::isfinite(*s.massDensityKgM3)||*s.massDensityKgM3<=0)
    return Core::Status::Failure(Core::StatusCode::InvalidState,"material primitive is absent/nonfinite/nonpositive");
  if((q.required&Bit(Capability::Species))&&(!s.species||s.species->empty()))
    return Core::Status::Failure(Core::StatusCode::UnsupportedCapability,"explicit species densities are unavailable");
  // Validate every supplied field, even if this caller did not request it.
  // Missing optional physics is allowed; a present NaN is never absence.
  if(s.species) {
    std::vector<std::string> seen;
    for(const auto& ion:*s.species) {
      if(ion.identity.empty()||!std::isfinite(ion.numberDensityM3)||ion.numberDensityM3<0)
        return Core::Status::Failure(Core::StatusCode::InvalidState,"invalid explicit species density");
      for(const auto& prior:seen)if(prior==ion.identity)
        return Core::Status::Failure(Core::StatusCode::InvalidState,"duplicate species identity");
      seen.push_back(ion.identity);
      if(ion.temperatureK&&(!std::isfinite(*ion.temperatureK)||*ion.temperatureK<=0))
        return Core::Status::Failure(Core::StatusCode::InvalidState,"invalid supplied species temperature");
    }
  }
  if(q.required&Bit(Capability::Thermodynamics)) {
    if(!s.pressurePa||!s.gammaAdiabatic)return Core::Status::Failure(Core::StatusCode::UnsupportedCapability,"thermal state is unavailable");
    if(s.species)for(const auto& species:*s.species)
      if(species.numberDensityM3>0&&(!species.temperatureK||
          !std::isfinite(*species.temperatureK)||*species.temperatureK<=0))
        return Core::Status::Failure(Core::StatusCode::UnsupportedCapability,"populated species has no valid temperature");
  }
  if((s.pressurePa&&(!std::isfinite(*s.pressurePa)||*s.pressurePa<=0))||
      (s.gammaAdiabatic&&(!std::isfinite(*s.gammaAdiabatic)||*s.gammaAdiabatic<=1)))
    return Core::Status::Failure(Core::StatusCode::InvalidState,"invalid supplied pressure/EOS");
  for(const auto& pair:{std::make_pair(Capability::VelocityGradient,&s.velocityGradientPerS),
                         std::make_pair(Capability::MagneticGradient,&s.magneticGradientTPerM)})
    {
      if((q.required&Bit(pair.first))&&!*pair.second)
        return Core::Status::Failure(Core::StatusCode::UnsupportedCapability,"full branch-preserving derivative tensor is unavailable");
      if(*pair.second) {
        if(!(**pair.second).branchPreserved||(**pair.second).branchIdentity.empty())
          return Core::Status::Failure(Core::StatusCode::InvalidState,"derivative does not preserve its physical branch");
        for(double value:(**pair.second).values)if(!std::isfinite(value))
          return Core::Status::Failure(Core::StatusCode::NumericalFailure,"nonfinite derivative tensor");
      }
    }
  return Core::Status::Success();
}

// Optional material-map and surface records carry their physical units and
// absence explicitly. Stage 11 supplies no unchosen lateral flow, inverse,
// thermal law or downstream volume to populate them.
struct Deformation {Tensor gradient=UninitializedTensor,ratePerS=UninitializedTensor;double jacobian=Uninitialized;double referenceEpochS=Uninitialized;};
struct Motion {Vector positionM=UninitializedVector,velocityMPerS=UninitializedVector;std::optional<Deformation> deformation;};
struct Surface {
  std::string stableIdentity,eventIdentity,classification;
  double epochS=Uninitialized,normalSpeedMPerS=Uninitialized;
  Vector positionM=UninitializedVector,outwardNormal=UninitializedVector;
  std::optional<State> upstream,downstream;
};
class Owner {
 public:
  virtual ~Owner()=default;
  virtual const Descriptor& Identity() const=0;
  virtual Core::Result<State> Sample(const Query&) const=0;
  virtual Core::Result<Motion> Move(const Vector&,double,double) const {
    return Core::Result<Motion>::Failure(Core::StatusCode::UnsupportedCapability,"material trajectory is unavailable");
  }
  virtual Core::Result<Vector> Inverse(const Vector&,double) const {
    return Core::Result<Vector>::Failure(Core::StatusCode::UnsupportedCapability,"material inverse is unavailable");
  }
  virtual Core::Result<Surface> Shock(const Vector&,double) const {
    return Core::Result<Surface>::Failure(Core::StatusCode::UnsupportedCapability,"shock surface is unavailable");
  }
  virtual Core::Result<Surface> Contact(const Vector&,double) const {
    return Core::Result<Surface>::Failure(Core::StatusCode::UnsupportedCapability,"contact surface is unavailable");
  }
};
class Interval {
 public:
  static Core::Result<std::shared_ptr<const Interval>> Create(std::shared_ptr<const Owner> owner) {
    using R=Core::Result<std::shared_ptr<const Interval>>;
    if(!owner)return R::Failure(Core::StatusCode::InvalidConfiguration,"material interval has no physical owner");
    const auto valid=Validate(owner->Identity());if(!valid.ok())return R::Failure(valid.code,valid.message);
    return R::Success(std::shared_ptr<const Interval>(new Interval(std::move(owner))));
  }
  const Descriptor& Identity() const {return descriptor_;}
  Core::Status RequireStep(double begin,double end,double fixedDt,Region region) const {
    const auto stable=Stable();if(!stable.ok())return stable;
    const auto admitted=Require(descriptor_.profile,Use::MaterialTransport,region,descriptor_.revision);
    if(!admitted.ok())return admitted;
    if((descriptor_.capabilities&MaterialRequirements)!=MaterialRequirements)
      return Core::Status::Failure(Core::StatusCode::UnsupportedCapability,"material owner lacks required primitive/derivative/trajectory capabilities");
    if(!std::isfinite(fixedDt)||fixedDt<=0||!std::isfinite(begin)||!std::isfinite(end)||
        end!=begin+fixedDt||begin<descriptor_.beginS||end>descriptor_.endS)
      return Core::Status::Failure(Core::StatusCode::OutOfDomain,"complete fixed AMPS step is not covered; shortening dt is forbidden");
    return Core::Status::Success();
  }
  Core::Result<State> Sample(const Query& query) const {
    using R=Core::Result<State>;const auto valid=Check(query.epochS,query.required|Bit(Capability::Primitives));
    if(!valid.ok())return R::Failure(valid.code,valid.message);
    if(!Finite(query.positionM)||query.region<Region::ParkerExterior||query.region>Region::PostEvent||
        (query.side!=Side::Smooth&&query.side!=Side::UpstreamLimit&&query.side!=Side::DownstreamLimit)||
        (query.side!=Side::Smooth&&query.interfaceIdentity.empty()))
      return R::Failure(Core::StatusCode::InvalidConfiguration,"query lacks physical position/region/one-sided interface identity");
    const auto value=owner_->Sample(query);if(!value.ok())return value;
    const auto stable=Stable();if(!stable.ok())return R::Failure(stable.code,stable.message);
    const auto matched=Validate(descriptor_,query,value.value);
    if(!matched.ok())return R::Failure(matched.code,matched.message);
    return value;
  }
  Core::Result<Motion> Move(const Vector& initial,double begin,double end) const {
    using R=Core::Result<Motion>;
    for(double t:{begin,end}) {const auto s=Check(t,Bit(Capability::Trajectory));if(!s.ok())return R::Failure(s.code,s.message);}
    if(!Finite(initial))return R::Failure(Core::StatusCode::InvalidConfiguration,"invalid material trajectory point");
    const auto value=owner_->Move(initial,begin,end);if(!value.ok())return value;
    const auto stable=Stable();if(!stable.ok())return R::Failure(stable.code,stable.message);
    if(!Finite(value.value.positionM)||!Finite(value.value.velocityMPerS))
      return R::Failure(Core::StatusCode::NumericalFailure,"nonfinite material trajectory");
    if(value.value.deformation)
      return R::Failure(Core::StatusCode::UnsupportedCapability,"this profile has no qualified general deformation export");
    return value;
  }
  Core::Result<Vector> Inverse(const Vector& position,double epoch) const {
    using R=Core::Result<Vector>;const auto valid=Check(epoch,Bit(Capability::InverseMap));
    if(!valid.ok())return R::Failure(valid.code,valid.message);
    if(!Finite(position))return R::Failure(Core::StatusCode::InvalidConfiguration,"invalid inverse material position");
    const auto result=owner_->Inverse(position,epoch);if(!result.ok())return result;
    const auto stable=Stable();if(!stable.ok())return R::Failure(stable.code,stable.message);
    if(!Finite(result.value))return R::Failure(Core::StatusCode::NumericalFailure,"nonfinite inverse material position");
    return result;
  }
  Core::Result<Surface> Shock(const Vector& position,double epoch) const {
    return SurfaceAt(position,epoch,false);
  }
  Core::Result<Surface> Contact(const Vector& position,double epoch) const {
    return SurfaceAt(position,epoch,true);
  }
  Core::Status BindSnapshot(const BackgroundSnapshot& snapshot) const {
    const auto stable=Stable();if(!stable.ok())return stable;
    // Existing SnapshotStore remains the sole native publication/read-phase
    // authority. This is preparation-time validation, not another global cache.
    const auto admitted=Require(descriptor_.profile,Use::MaterialTransport,Region::ParkerExterior,descriptor_.revision);
    if(!admitted.ok())return admitted;
    // A profile's catalog eligibility cannot replace actual owner data.
    // Preparation must reject a sparse snapshot-only realization before the
    // sole native store can attach it to a material transport step.
    if((descriptor_.capabilities&MaterialRequirements)!=MaterialRequirements)
      return Core::Status::Failure(Core::StatusCode::UnsupportedCapability,
          "native material binding lacks required owner capabilities");
    if(snapshot.provider()!=Provider::ReducedShock||snapshot.ownership()!=Ownership::ModelOwned||
        descriptor_.nativeConfigurationIdentity.empty()||
        snapshot.configuration_fingerprint()!=descriptor_.nativeConfigurationIdentity||
        snapshot.field_line_generation()!=descriptor_.snapshotGeneration||
        !snapshot.Covers(descriptor_.beginS)||!snapshot.Covers(descriptor_.endS))
      return Core::Status::Failure(Core::StatusCode::DataIntegrityFailure,"native snapshot and material realization differ");
    return Core::Status::Success();
  }
 private:
  explicit Interval(std::shared_ptr<const Owner> owner):descriptor_(owner->Identity()),owner_(std::move(owner)){}
  Core::Result<Surface> SurfaceAt(const Vector& position,double epoch,bool contact) const {
    using R=Core::Result<Surface>;
    const auto valid=Check(epoch,Bit(contact?Capability::ContactSurface:Capability::ShockSurface));
    if(!valid.ok())return R::Failure(valid.code,valid.message);
    if(!Finite(position))return R::Failure(Core::StatusCode::InvalidConfiguration,"invalid physical surface point");
    const auto value=contact?owner_->Contact(position,epoch):owner_->Shock(position,epoch);
    if(!value.ok())return value;
    const auto stable=Stable();if(!stable.ok())return R::Failure(stable.code,stable.message);
    const auto& s=value.value;
    if(s.stableIdentity.empty()||s.classification.empty()||s.eventIdentity!=descriptor_.eventIdentity||
        s.epochS!=epoch||s.positionM!=position||!Finite(s.outwardNormal)||!std::isfinite(s.normalSpeedMPerS))
      return R::Failure(Core::StatusCode::DataIntegrityFailure,"foreign or invalid surface realization");
    const double norm=std::hypot(s.outwardNormal[0],s.outwardNormal[1],s.outwardNormal[2]);
    if(std::abs(norm-1)>512*std::numeric_limits<double>::epsilon())
      return R::Failure(Core::StatusCode::InvalidState,"surface normal must already be unit length; no direction fallback");
    for(const auto& side:{std::make_pair(&s.upstream,Side::UpstreamLimit),std::make_pair(&s.downstream,Side::DownstreamLimit)})
      if(*side.first) {
        const auto& primitive=**side.first;
        if(primitive.region<Region::ParkerExterior||primitive.region>Region::PostEvent)
          return R::Failure(Core::StatusCode::InvalidConfiguration,"one-sided surface primitive has no known physical region");
        Query q;q.positionM=position;q.epochS=epoch;q.region=primitive.region;
        q.side=side.second;q.interfaceIdentity=s.stableIdentity;
        const auto matched=Validate(descriptor_,q,primitive);
        if(!matched.ok())return R::Failure(matched.code,matched.message);
      }
    // Valid local jump limits remain distinct from PhysicalDownstream use.
    return value;
  }
  Core::Status Stable() const {
    if(Payload(descriptor_)!=Payload(owner_->Identity()))
      return Core::Status::Failure(Core::StatusCode::DataIntegrityFailure,"physical owner changed after immutable interval capture");
    return Core::Status::Success();
  }
  Core::Status Check(double epoch,std::uint64_t required) const {
    const auto stable=Stable();if(!stable.ok())return stable;
    if(!std::isfinite(epoch)||epoch<descriptor_.beginS||epoch>descriptor_.endS)
      return Core::Status::Failure(Core::StatusCode::OutOfDomain,"query outside complete material interval");
    if((required&~AllCapabilities)||(required&descriptor_.capabilities)!=required)
      return Core::Status::Failure(Core::StatusCode::UnsupportedCapability,"requested physical quantity is unavailable");
    return Core::Status::Success();
  }
  const Descriptor descriptor_;
  const std::shared_ptr<const Owner> owner_;
};
} } }
#endif
