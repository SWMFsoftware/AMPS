#include "swcme1d_adapter.h"

#include "../util/sep_common_header_path.h"
#include SRCSEP_SEP_COMMON_HEADER(sep_background_snapshot.h)
#include SRCSEP_SEP_COMMON_HEADER(sep_material_background.h)

// The makefile's adapter-only rule supplies AMPS/src/models/swcme as an
// include root. This is the sole srcSEP production translation unit that
// consumes the canonical 1-D SWCME C++ interface.
#include "swcme1d.hpp"
#include "swcme1d_input.hpp"
#include "swcme_material_upstream.hpp"
#include "swcme1d_material_adapter.h"

#include <algorithm>
#include <atomic>
#include <cmath>
#include <iomanip>
#include <limits>
#include <memory>
#include <sstream>
#include <stdexcept>

namespace {

// Reconfiguration replaces the canonical instance rather than mutating a
// prepared Model, whose canonical lifecycle is intentionally immutable.
std::unique_ptr<swcme1d::Model> model(new swcme1d::Model());
swcme1d::StepState state{};
bool state_prepared = false;
bool clamp_sheath = true;
double prepared_epoch_seconds = 0.0;
std::uint64_t prepared_state_id = 0;
swcme::input1d::ResolvedConfiguration resolved_configuration;
bool configuration_resolved = false;
std::vector<SEP::SW1DAdapter::ParameterAssignment> staged_input_assignments;

SEP::SW1DAdapter::FailurePolicy failure_policy =
    SEP::SW1DAdapter::FailurePolicy::Strict;
SEP::SW1DAdapter::BackgroundSample diagnostic_fallback;
bool diagnostic_fallback_configured = false;

// Particle queries may execute concurrently. Counters are observational only
// and therefore use relaxed atomics; no physical state is synchronized through
// them. State preparation remains a lifecycle operation outside particle
// advancement and is serialized by the srcSEP driver.
std::atomic<std::uint64_t> successful_queries{0};
std::atomic<std::uint64_t> failed_queries{0};
std::atomic<std::uint64_t> radius_clamps{0};
std::atomic<std::uint64_t> diagnostic_fallbacks{0};

void AssertSwcmeMayWrite() {
  SEP::Background::SnapshotStore::Instance().AssertProviderMayWrite(
      SEP::Background::Provider::Swcme);
}

void RequirePrepared(const char* operation) {
  if (!state_prepared) {
    throw std::logic_error(std::string(operation) +
                           ": SWCME state has not been prepared");
  }
}

SEP::SW1DAdapter::StatusCode MapCanonicalStatus(
    const swcme::ModelStatus& status) {
  switch (status.code) {
    case swcme::StatusCode::NonFiniteInput:
      return SEP::SW1DAdapter::StatusCode::NonFiniteInput;
    case swcme::StatusCode::OutsideModelDomain:
      return SEP::SW1DAdapter::StatusCode::RadiusOutsideDomain;
    case swcme::StatusCode::NonFiniteResult:
      return SEP::SW1DAdapter::StatusCode::NonFiniteOutput;
    default:
      return SEP::SW1DAdapter::StatusCode::CanonicalFailure;
  }
}

std::string CanonicalDetail(const swcme::ModelStatus& status) {
  std::ostringstream out;
  out << swcme::status_code_name(status.code);
  if (status.context && status.context[0] != '\0') out << " at " << status.context;
  if (status.sample_index != swcme::ModelStatus::npos)
    out << " sample=" << status.sample_index;
  if (status.has_offending_value)
    out << std::setprecision(17) << " value=" << status.offending_value;
  return out.str();
}

SEP::SW1DAdapter::QueryResult FailureResult(
    SEP::SW1DAdapter::StatusCode code, double requested_radius_m,
    double evaluated_radius_m, const std::string& field,
    const std::string& detail) {
  SEP::SW1DAdapter::QueryResult result;
  result.code = code;
  result.original_code = code;
  result.requested_radius_m = requested_radius_m;
  result.evaluated_radius_m = evaluated_radius_m;
  result.epoch_seconds = prepared_epoch_seconds;
  result.source_state_id = prepared_state_id;
  result.failed_field = field;
  result.detail = detail;
  return result;
}

SEP::SW1DAdapter::QueryResult RecoverOrFail(
    SEP::SW1DAdapter::QueryResult result) {
  if (failure_policy == SEP::SW1DAdapter::FailurePolicy::DiagnosticFallback &&
      diagnostic_fallback_configured) {
    const SEP::SW1DAdapter::StatusCode original = result.code;
    result.code = SEP::SW1DAdapter::StatusCode::RecoveredWithDiagnosticFallback;
    result.original_code = original;
    result.sample = diagnostic_fallback;
    result.detail = "diagnostic fallback used after: " + result.detail;
    diagnostic_fallbacks.fetch_add(1, std::memory_order_relaxed);
    successful_queries.fetch_add(1, std::memory_order_relaxed);
    return result;
  }

  failed_queries.fetch_add(1, std::memory_order_relaxed);
  return result;
}

}  // namespace

namespace SEP {
namespace SW1DAdapter {

const char* StatusCodeName(StatusCode code) {
  switch (code) {
    case StatusCode::Success: return "success";
    case StatusCode::RecoveredWithRadiusClamp: return "recovered-radius-clamp";
    case StatusCode::RecoveredWithDiagnosticFallback:
      return "recovered-diagnostic-fallback";
    case StatusCode::ModelNotPrepared: return "model-not-prepared";
    case StatusCode::NonFiniteInput: return "nonfinite-input";
    case StatusCode::RadiusOutsideDomain: return "radius-outside-domain";
    case StatusCode::NonFiniteOutput: return "nonfinite-output";
    case StatusCode::NonPositiveDensity: return "nonpositive-density";
    case StatusCode::InvalidSpeed: return "invalid-speed";
    case StatusCode::InvalidDivergence: return "invalid-divergence";
    case StatusCode::CanonicalFailure: return "canonical-failure";
    case StatusCode::ConfigurationFailure: return "configuration-failure";
    case StatusCode::PreparationFailure: return "preparation-failure";
    case StatusCode::InvalidDiagnosticFallback:
      return "invalid-diagnostic-fallback";
  }
  return "unknown";
}

const char* FailurePolicyName(FailurePolicy policy) {
  switch (policy) {
    case FailurePolicy::Strict: return "strict";
    case FailurePolicy::ClampRadius: return "clamp-radius";
    case FailurePolicy::DiagnosticFallback: return "diagnostic-fallback";
  }
  return "unknown";
}

Status Configure(const ConfigurationRequest& request) {
  AssertSwcmeMayWrite();

  // Convert the provider-neutral transport records into canonical assignments
  // only inside this private translation unit. Staged PARAM and explicit input
  // assignments share one authority layer, so duplicates are rejected rather
  // than becoming order-dependent. CLI/programmatic values form the later,
  // intentionally overriding layer.
  const auto append_assignments=[](
      const std::vector<ParameterAssignment>& source,
      std::vector<swcme::input1d::Assignment>* destination) {
    for (std::size_t i=0;i<source.size();++i) {
      swcme::input1d::Assignment assignment;
      assignment.key=source[i].key;
      assignment.value=source[i].value;
      assignment.origin=source[i].origin;
      assignment.line=source[i].line;
      destination->push_back(assignment);
    }
  };

  swcme::input1d::Layer input_layer;
  input_layer.name="input-file";
  append_assignments(staged_input_assignments,&input_layer.assignments);
  append_assignments(request.input_assignments,&input_layer.assignments);
  swcme::input1d::Layer command_line_layer;
  command_line_layer.name="command-line/programmatic";
  append_assignments(request.command_line_assignments,
                     &command_line_layer.assignments);

  std::vector<swcme::input1d::Layer> layers;
  if (!input_layer.assignments.empty()) layers.push_back(input_layer);
  if (!command_line_layer.assignments.empty()) layers.push_back(command_line_layer);
  const swcme::input1d::ResolveResult resolved=swcme::input1d::Resolve(
      request.preset==Scenario::Fast ? swcme::input1d::Preset::Fast
                                     : swcme::input1d::Preset::Slow,
      layers);
  if (!resolved.ok()) {
    std::ostringstream detail;
    detail << swcme::input1d::code_name(resolved.status.code)
           << " key=" << resolved.status.key
           << " layer=" << resolved.status.layer;
    if (!resolved.status.origin.empty())
      detail << " origin=" << resolved.status.origin;
    if (resolved.status.line!=0) detail << " line=" << resolved.status.line;
    if (!resolved.status.message.empty())
      detail << " detail=" << resolved.status.message;
    return {StatusCode::ConfigurationFailure,detail.str()};
  }
  if (resolved.configuration.source.spectrum.normalization!=
      swcme::sep::NormalizationMode::RelativeOnly) {
    return {StatusCode::ConfigurationFailure,
            "UNSUPPORTED_FIELD key=source.normalization detail=srcSEP's "
            "swept-volume source supports relative_only; calibrated reference "
            "intensity is retained by canonical SWCME but has no 1-D AMPS "
            "source adapter yet"};
  }

  // Model construction and canonical validation complete before either global
  // record is replaced. A malformed override cannot partially reconfigure a
  // previously valid model.
  std::unique_ptr<swcme1d::Model> candidate(
      new swcme1d::Model(resolved.configuration.model));
  model.swap(candidate);
  resolved_configuration=resolved.configuration;
  configuration_resolved=true;
  clamp_sheath=true;
  state_prepared=false;
  prepared_epoch_seconds=0.0;
  prepared_state_id=0;
  failure_policy=FailurePolicy::Strict;
  diagnostic_fallback_configured=false;
  ResetDiagnostics();
  return {};
}

void Configure(Scenario scenario) {
  ConfigurationRequest request;
  request.preset=scenario;
  const Status status=Configure(request);
  if (!status.ok()) throw std::invalid_argument(status.detail);
}

Background::Material::Profile MaterialProfile() {
  // A caller may query admission before preparing an epoch, but must already
  // have resolved the canonical configuration. Do not infer a region mode
  // from a filename, scenario, function pointer or an uninitialized Model.
  using P=Background::Material::Profile;
  if(!configuration_resolved)return static_cast<P>(0);
  switch(resolved_configuration.model.region_mode) {
    case swcme::regions::Mode::ShockOnly:return P::Swcme1dShockOnlyV1;
    case swcme::regions::Mode::FullICME:return P::Swcme1dFullIcmeDiagnosticV1;
  }
  return static_cast<P>(0);
}

namespace {
namespace Material=Background::Material;
namespace Upstream=swcme::material1d;

template<class T> Core::Result<T> MaterialFailure(const swcme::ModelStatus& status) {
  using C=Core::StatusCode;
  C code=C::InvalidState;
  switch(status.code) {
    case swcme::StatusCode::OutsideModelDomain: code=C::OutOfDomain;break;
    case swcme::StatusCode::NonFiniteInput:
    case swcme::StatusCode::DegenerateVector: code=C::InvalidConfiguration;break;
    case swcme::StatusCode::ShockSolverFailure:
    case swcme::StatusCode::NonFiniteResult: code=C::NumericalFailure;break;
    default: break;
  }
  return Core::Result<T>::Failure(code,CanonicalDetail(status));
}

class UpstreamMaterialOwner final : public Material::Owner {
 public:
  UpstreamMaterialOwner(std::shared_ptr<const Upstream::Event> event,
      std::uint64_t generation,std::uint64_t geometry)
      :event_(std::move(event)),descriptor_(Describe(*event_,generation,geometry)) {}

  const Material::Descriptor& Identity() const override {return descriptor_;}

  Core::Result<Material::State> Sample(const Material::Query& q) const override {
    using R=Core::Result<Material::State>;
    // `ParkerExterior` here is the smooth prescribed ambient. The shock is
    // a separately located surface and the core can furnish its RH limits.
    // Claiming a particular downstream region is outside this profile.
    if(q.region!=Material::Region::ParkerExterior||q.side!=Material::Side::Smooth||
       !q.interfaceIdentity.empty())
      return R::Failure(Core::StatusCode::UnsupportedCapability,"SWCME1D material sample requires smooth prescribed ambient");
    Upstream::Sample value;
    const auto status=event_->sample(q.positionM,q.epochS,value);
    if(!status.ok())return MaterialFailure<Material::State>(status);
    return R::Success(Ambient(value,q.positionM,q.epochS));
  }

  Core::Result<Material::Motion> Move(const Material::Vector& p,double begin,double end) const override {
    Upstream::Motion value;
    const auto status=event_->move(p,begin,end,value);
    if(!status.ok())return MaterialFailure<Material::Motion>(status);
    Material::Motion motion;
    motion.positionM=value.position_m;
    motion.velocityMPerS=value.velocity_m_s;
    // Stage 11's shared interface does not yet admit a general deformation
    // result. The exact Jacobian is available from the canonical upstream
    // trajectory API; a neutral caller gets the checked position/velocity.
    return Core::Result<Material::Motion>::Success(std::move(motion));
  }

  Core::Result<Material::Surface> Shock(const Material::Vector& p,double epoch) const override {
    using R=Core::Result<Material::Surface>;
    Upstream::Shock value;
    const auto status=event_->shock_at(p,epoch,value);
    if(!status.ok())return MaterialFailure<Material::Surface>(status);
    const auto& x=value.position_m;
    // Require the actual inertial front. Allow only roundoff from projecting
    // the canonical radius onto the same Cartesian ray again. This bound
    // represents point-coordinate precision, not a physical sheath width.
    if(std::hypot(x[0]-p[0],x[1]-p[1],x[2]-p[2])>
       32.0*std::numeric_limits<double>::epsilon()*std::hypot(x[0],x[1],x[2]))
      return R::Failure(Core::StatusCode::OutOfDomain,"SWCME1D query is not on the prescribed spherical front");
    Material::Surface surface;
    surface.stableIdentity=descriptor_.eventIdentity+":swcme1d-leading-front";
    surface.eventIdentity=descriptor_.eventIdentity;
    surface.classification=swcme::shock::solve_status_name(value.jump.status);
    surface.epochS=epoch;
    surface.positionM=p;
    surface.outwardNormal=value.outward_normal;
    surface.normalSpeedMPerS=value.jump.shock_normal_speed_m_s;
    surface.upstream=Ambient(value.ambient,p,epoch);
    surface.upstream->side=Material::Side::UpstreamLimit;
    surface.upstream->interfaceIdentity=surface.stableIdentity;
    if(value.jump.has_shock) {
      auto down=IdentityAt(p,epoch);
      down.region=Material::Region::Sheath;
      down.side=Material::Side::DownstreamLimit;
      down.interfaceIdentity=surface.stableIdentity;
      down.massDensityKgM3=value.jump.downstream.rho_kg_m3;
      down.pressurePa=value.jump.downstream.pressure_Pa;
      down.gammaAdiabatic=event_->profile().config().gamma;
      down.velocityMPerS=value.jump.downstream.velocity_m_s;
      down.magneticFieldT=value.jump.downstream.magnetic_T;
      // Total downstream pressure does not fix species temperatures or an
      // evolution law. Export only the real one-sided RH bulk state.
      surface.downstream=std::move(down);
    }
    return R::Success(std::move(surface));
  }

 private:
  static Material::Descriptor Describe(const Upstream::Event& event,
      std::uint64_t generation,std::uint64_t geometry) {
    // Freeze identity, composition, frame and validity bounds while the
    // event is captured. The event and geometry generations distinguish
    // physical records with identical coordinates at different host steps.
    // Only smooth upstream fluid and the one-sided shock surface belong
    // to this descriptor, so no volume capability is advertised behind it.
    Material::Descriptor result;
    result.profile=Material::Profile::Swcme1dAdiabaticOuterWindV1;
    result.revision=1;
    result.capabilities=Material::MaterialRequirements|Material::Bit(Material::Capability::ShockSurface);
    result.backgroundGeneration=generation;
    result.geometryGeneration=geometry;
    result.geometryGenerationPresent=true;
    result.eventIdentity=event.config().identity;
    result.physicalConfigurationIdentity=std::to_string(event.fingerprint());
    result.materialIdentity=result.eventIdentity+":radial-reference";
    result.compositionIdentity="charge-neutral-electrons-protons-no-alpha";
    result.eosIdentity="scalar-multispecies-ideal-gas-adiabatic";
    result.frame="HCI";
    result.units="SI";
    result.beginS=event.config().begin_s;
    result.endS=event.config().end_s;
    result.densityBasis=Material::DensityBasis::CoreNonrelativisticFluid;
    return result;
  }

  Material::State IdentityAt(const Material::Vector& position,double epoch) const {
    Material::State state;
    state.eventIdentity=descriptor_.eventIdentity;
    state.physicalConfigurationIdentity=descriptor_.physicalConfigurationIdentity;
    state.backgroundGeneration=descriptor_.backgroundGeneration;
    state.geometryGeneration=descriptor_.geometryGeneration;
    state.positionM=position;
    state.epochS=epoch;
    return state;
  }

  Material::State Ambient(const Upstream::Sample& value,const Material::Vector& position,double epoch) const {
    auto state=IdentityAt(position,epoch);
    state.region=Material::Region::ParkerExterior;
    state.velocityMPerS=value.velocity_m_s;
    state.magneticFieldT=value.magnetic_field_t;
    state.massDensityKgM3=value.mass_density_kg_m3;
    state.pressurePa=value.pressure_pa;
    state.gammaAdiabatic=event_->profile().config().gamma;
    state.species=std::vector<Material::Species>{{"electron",value.electron_density_m3,value.electron_temperature_k},
        {"proton",value.proton_density_m3,value.proton_temperature_k},{"alpha",0.0,std::nullopt}};
    const std::string branch=Upstream::profile_key();
    state.velocityGradientPerS=Material::Gradient{value.velocity_gradient_s_inv,branch,true};
    state.magneticGradientTPerM=Material::Gradient{value.magnetic_gradient_t_m_inv,branch,true};
    return state;
  }

  const std::shared_ptr<const Upstream::Event> event_;
  const Material::Descriptor descriptor_;
};
} // namespace

Core::Result<std::shared_ptr<const Background::Material::Interval>> CaptureUpstreamMaterial(
    std::shared_ptr<const swcme::material1d::Event> event,
    std::uint64_t background_generation,std::uint64_t geometry_generation) {
  using R=Core::Result<std::shared_ptr<const Background::Material::Interval>>;
  if(!event)
    return R::Failure(Core::StatusCode::InvalidConfiguration,"SWCME1D material export requires a canonical event");
  return Background::Material::Interval::Create(
      std::make_shared<const UpstreamMaterialOwner>(std::move(event),background_generation,geometry_generation));
}

ConfigurationSummary GetConfigurationSummary() {
  if (!configuration_resolved)
    throw std::logic_error("SWCME configuration has not been resolved");
  ConfigurationSummary summary;
  summary.preset=swcme::input1d::preset_name(resolved_configuration.preset);
  summary.fingerprint=resolved_configuration.fingerprint;
  summary.normalized_manifest=resolved_configuration.normalized_manifest;
  summary.launch_epoch_s=resolved_configuration.launch_epoch_s;
  summary.valid_from_s=resolved_configuration.valid_from_s;
  summary.valid_until_s=resolved_configuration.valid_until_s;
  summary.source_particle_mass_kg=
      resolved_configuration.source.spectrum.particle_mass_kg;
  summary.source_charge_number=
      resolved_configuration.source.spectrum.charge_number;
  summary.source_energy_min_MeV=
      resolved_configuration.source.spectrum.kinetic_energy_min_MeV;
  summary.source_energy_max_MeV=
      resolved_configuration.source.spectrum.kinetic_energy_max_MeV;
  summary.source_reference_energy_MeV=
      resolved_configuration.source.spectrum.reference_energy_MeV;
  summary.source_injection_efficiency=
      resolved_configuration.source.injection_efficiency;
  summary.relative_source_weight_per_area=
      resolved_configuration.model.relative_source_weight_per_area;
  summary.ambient_wind_speed_m_per_s=
      resolved_configuration.model.V_sw_kms*1000.0;
  summary.solar_rotation_rate_rad_per_s=
      swcme::defaults::SOLAR_ROTATION_RATE_RAD_S;
  summary.parker_source_radius_m=
      resolved_configuration.model.parker_source_radius_Rs*
      swcme::constants::SOLAR_RADIUS_M;
  summary.parker_reference_sin_theta=resolved_configuration.model.sin_theta;
  summary.parker_radial_polarity=
      resolved_configuration.model.parker_radial_polarity;
  // ambient.magnetic_field_1au is canonical total |B| at the configured
  // reference latitude, not Br.  Remove the Parker azimuthal contribution at
  // one AU and retain polarity in the signed radial normalization consumed by
  // srcSEP's field-line/domain field constructor.
  const double reference_winding=
      summary.solar_rotation_rate_rad_per_s*
      (swcme::constants::AU_M-summary.parker_source_radius_m)/
      summary.ambient_wind_speed_m_per_s*
      summary.parker_reference_sin_theta;
  summary.parker_radial_field_at_one_au_t=
      static_cast<double>(summary.parker_radial_polarity)*
      resolved_configuration.model.B1AU_nT*1.0e-9/
      std::sqrt(1.0+reference_winding*reference_winding);
  return summary;
}

void ClearStagedInputAssignments() { staged_input_assignments.clear(); }

Status StageInputAssignment(const ParameterAssignment& assignment) {
  if (assignment.key.empty() || assignment.value.empty())
    return {StatusCode::ConfigurationFailure,
            "SWCME input assignment requires non-empty key and value"};
  staged_input_assignments.push_back(assignment);
  return {};
}

Status PrepareState(double epoch_seconds) {
  AssertSwcmeMayWrite();
  if (!std::isfinite(epoch_seconds)) {
    return {StatusCode::PreparationFailure,
            "SWCME epoch is not finite"};
  }
  if (!configuration_resolved) {
    return {StatusCode::PreparationFailure,
            "SWCME configuration has not been resolved"};
  }
  if (epoch_seconds<resolved_configuration.valid_from_s ||
      epoch_seconds>resolved_configuration.valid_until_s) {
    std::ostringstream detail;
    detail << std::setprecision(17)
           << "simulation epoch " << epoch_seconds
           << " is outside configured validity ["
           << resolved_configuration.valid_from_s << ','
           << resolved_configuration.valid_until_s << ']';
    return {StatusCode::PreparationFailure,detail.str()};
  }

  // Construct the candidate first. Canonical validation may throw; assigning
  // only after success preserves the previous valid cache and snapshot.
  try {
    const double launch_relative_seconds=
        epoch_seconds-resolved_configuration.launch_epoch_s;
    const swcme1d::StepState candidate =
        model->prepare_step(launch_relative_seconds);
    state = candidate;
    state_prepared = true;
    prepared_epoch_seconds = epoch_seconds;
    ++prepared_state_id;
    return {};
  } catch (const std::exception& exception) {
    return {StatusCode::PreparationFailure, exception.what()};
  } catch (...) {
    return {StatusCode::PreparationFailure,
            "canonical SWCME preparation raised a non-standard exception"};
  }
}

Status SetFailurePolicy(FailurePolicy policy) {
  if (policy == FailurePolicy::DiagnosticFallback &&
      !diagnostic_fallback_configured) {
    return {StatusCode::InvalidDiagnosticFallback,
            "diagnostic-fallback policy requires a validated fallback sample"};
  }
  failure_policy = policy;
  return {};
}

Status ValidateSample(const BackgroundSample& sample) {
  if (!std::isfinite(sample.number_density_m3())) {
    return {StatusCode::NonFiniteOutput,
            "number density is not finite [m^-3]"};
  }
  if (!(sample.number_density_m3() > 0.0)) {
    return {StatusCode::NonPositiveDensity,
            "number density must be strictly positive [m^-3]"};
  }
  if (!std::isfinite(sample.speed_m_s())) {
    return {StatusCode::NonFiniteOutput,
            "radial speed is not finite [m/s]"};
  }
  if (!(sample.speed_m_s() > 0.0)) {
    return {StatusCode::InvalidSpeed,
            "radial speed must be strictly positive [m/s]"};
  }
  // Divergence may be positive, zero, or negative; finiteness is its only
  // universal physical invariant at this provider boundary.
  if (!std::isfinite(sample.divergence_s_inv())) {
    return {StatusCode::InvalidDivergence,
            "velocity divergence is not finite [s^-1]"};
  }
  return {};
}

Status SetDiagnosticFallback(const BackgroundSample& sample) {
  const Status status = ValidateSample(sample);
  if (!status.ok()) {
    return {StatusCode::InvalidDiagnosticFallback,
            "invalid diagnostic fallback: " + status.detail};
  }
  diagnostic_fallback = sample;
  diagnostic_fallback_configured = true;
  return {};
}

FailurePolicy GetFailurePolicy() { return failure_policy; }

void EnableSheathClamp(bool enabled) { clamp_sheath = enabled; }

double ShockRadiusM() {
  RequirePrepared("ShockRadiusM");
  return state.r_sh_m;
}

double ShockSpeedMPerS() {
  RequirePrepared("ShockSpeedMPerS");
  return state.V_sh_ms;
}

double CompressionRatio() {
  RequirePrepared("CompressionRatio");
  return state.rc;
}

double RelativeSourceWeightPerArea() {
  if (!configuration_resolved)
    throw std::logic_error(
        "RelativeSourceWeightPerArea: SWCME configuration is unresolved");
  return resolved_configuration.model.relative_source_weight_per_area;
}

double DlnB_Dr_at_r(double radius_m) {
  RequirePrepared("DlnB_Dr_at_r");
  if (!std::isfinite(radius_m) || radius_m < swcme::solarwind::MIN_RADIUS_M)
    throw std::domain_error("DlnB_Dr_at_r: radius is outside the SWCME domain");

  const double radius_au = radius_m / swcme1d::AU;
  const double k = state.k_AU;
  const double k2r2 = (k * radius_au) * (k * radius_au);
  return (-2.0 / radius_au +
          (k * k * radius_au) / (1.0 + k2r2)) /
         swcme1d::AU;
}

QueryResult QueryAtRadius(double radius_m, bool apply_sheath_clamp) {
  if (!state_prepared) {
    return RecoverOrFail(FailureResult(
        StatusCode::ModelNotPrepared, radius_m, radius_m, "state",
        "SWCME query requested before PrepareState"));
  }
  if (!std::isfinite(radius_m)) {
    return RecoverOrFail(FailureResult(
        StatusCode::NonFiniteInput, radius_m, radius_m, "radius_m",
        "query radius is not finite"));
  }

  double evaluated_radius_m = radius_m;
  bool used_radius_clamp = false;
  if (radius_m < swcme::solarwind::MIN_RADIUS_M) {
    if (failure_policy == FailurePolicy::ClampRadius) {
      evaluated_radius_m = swcme::solarwind::MIN_RADIUS_M;
      used_radius_clamp = true;
    } else {
      return RecoverOrFail(FailureResult(
          StatusCode::RadiusOutsideDomain, radius_m, radius_m, "radius_m",
          "radius is below the canonical 1.05-R_sun inner boundary"));
    }
  }

  // All canonical outputs are local temporaries. Even a future evaluator that
  // writes some fields before returning an error cannot leak a partial sample.
  double number_density_m3 = std::numeric_limits<double>::quiet_NaN();
  double speed_m_s = std::numeric_limits<double>::quiet_NaN();
  double radial_field_t = std::numeric_limits<double>::quiet_NaN();
  double azimuthal_field_t = std::numeric_limits<double>::quiet_NaN();
  double field_magnitude_t = std::numeric_limits<double>::quiet_NaN();
  double divergence_s_inv = std::numeric_limits<double>::quiet_NaN();
  const swcme::ModelStatus canonical_status =
      model->evaluate_radii_with_B_div_checked(
          state, &evaluated_radius_m, &number_density_m3, &speed_m_s,
          &radial_field_t, &azimuthal_field_t, &field_magnitude_t,
          &divergence_s_inv, 1);
  if (!canonical_status.ok()) {
    return RecoverOrFail(FailureResult(
        MapCanonicalStatus(canonical_status), radius_m, evaluated_radius_m,
        "canonical-evaluator", CanonicalDetail(canonical_status)));
  }

  if (clamp_sheath && apply_sheath_clamp && state.rc > 1.0 &&
      evaluated_radius_m <= state.r_sh_m &&
      evaluated_radius_m >= state.r_le_m) {
    const double upstream_density =
        swcme1d::Model::density_upstream(state, evaluated_radius_m);
    if (!std::isfinite(upstream_density) || !(upstream_density > 0.0)) {
      return RecoverOrFail(FailureResult(
          StatusCode::NonFiniteOutput, radius_m, evaluated_radius_m,
          "upstream-density", "sheath stabilization received invalid density"));
    }
    number_density_m3 = std::max(number_density_m3, upstream_density);
    const double minimum_speed = state.V_up_ms;
    const double maximum_speed = std::max(state.V_up_ms, state.V_sh_ms);
    speed_m_s = std::min(std::max(speed_m_s, minimum_speed), maximum_speed);
  }

  // Composition is explicitly frozen in the canonical prepared state.  The
  // single-fluid density compression applies to every species. Regional
  // pressure/heating, however, is not exposed by the canonical 1-D API:
  // advertise thermal fields only in its exact ambient SHOCK_ONLY branch.
  const auto thermodynamics=swcme::solarwind::thermodynamic_state(
      state.common.solar_wind,number_density_m3);
  const bool thermal_available=
      state.region_config.mode==swcme::regions::Mode::ShockOnly;
  const BackgroundSample sample(number_density_m3,speed_m_s,divergence_s_inv,
      radial_field_t,azimuthal_field_t,thermodynamics.mass_density_kg_m3,
      thermal_available,thermal_available?thermodynamics.pressure_Pa:0.0,
      thermal_available?state.common.solar_wind.T_K:0.0);
  const Status validation = ValidateSample(sample);
  if (!validation.ok()) {
    std::string field = "background-sample";
    if (validation.code == StatusCode::NonPositiveDensity) field = "density";
    if (validation.code == StatusCode::InvalidSpeed) field = "speed";
    if (validation.code == StatusCode::InvalidDivergence) field = "divergence";
    if (validation.code == StatusCode::NonFiniteOutput) {
      if (validation.detail.find("density") != std::string::npos) field = "density";
      if (validation.detail.find("speed") != std::string::npos) field = "speed";
    }
    return RecoverOrFail(FailureResult(
        validation.code, radius_m, evaluated_radius_m, field,
        validation.detail));
  }

  QueryResult result;
  result.code = used_radius_clamp ? StatusCode::RecoveredWithRadiusClamp
                                  : StatusCode::Success;
  result.original_code = used_radius_clamp ? StatusCode::RadiusOutsideDomain
                                           : StatusCode::Success;
  result.sample = sample;
  result.requested_radius_m = radius_m;
  result.evaluated_radius_m = evaluated_radius_m;
  result.epoch_seconds = prepared_epoch_seconds;
  result.source_state_id = prepared_state_id;
  if (used_radius_clamp) {
    result.failed_field = "radius_m";
    result.detail = "explicit clamp-radius policy evaluated the inner boundary";
    radius_clamps.fetch_add(1, std::memory_order_relaxed);
  }
  successful_queries.fetch_add(1, std::memory_order_relaxed);
  return result;
}

Diagnostics GetDiagnostics() {
  Diagnostics result;
  result.successful_queries = successful_queries.load(std::memory_order_relaxed);
  result.failed_queries = failed_queries.load(std::memory_order_relaxed);
  result.radius_clamps = radius_clamps.load(std::memory_order_relaxed);
  result.diagnostic_fallbacks =
      diagnostic_fallbacks.load(std::memory_order_relaxed);
  result.prepared_state_id = prepared_state_id;
  result.prepared_epoch_seconds = prepared_epoch_seconds;
  return result;
}

void ResetDiagnostics() {
  successful_queries.store(0, std::memory_order_relaxed);
  failed_queries.store(0, std::memory_order_relaxed);
  radius_clamps.store(0, std::memory_order_relaxed);
  diagnostic_fallbacks.store(0, std::memory_order_relaxed);
}

std::string FormatFailure(const QueryResult& result, int mpi_rank) {
  std::ostringstream out;
  out << std::setprecision(17)
      << "SWCME background query " << StatusCodeName(result.code)
      << " rank=" << mpi_rank
      << " epoch_s=" << result.epoch_seconds
      << " requested_radius_m=" << result.requested_radius_m
      << " evaluated_radius_m=" << result.evaluated_radius_m
      << " field=" << (result.failed_field.empty() ? "none" : result.failed_field)
      << " source_state_id=" << result.source_state_id;
  if (!result.detail.empty()) out << " detail=" << result.detail;
  return out.str();
}

void WriteRadialProfileFromR(const double* radii_m, int count,
                             const char* file_name, double epoch_seconds) {
  if (!state_prepared || radii_m == nullptr || count <= 0 ||
      file_name == nullptr) {
    return;
  }
  model->write_tecplot_radial_profile_from_r(
      state, radii_m, count, file_name, epoch_seconds);
}

void WriteShockVsTime(double duration_seconds, int sample_count,
                      const char* file_name) {
  if (sample_count <= 0 || file_name == nullptr) return;
  model->write_tecplot_shock_vs_time(duration_seconds, sample_count, file_name);
}

}  // namespace SW1DAdapter
}  // namespace SEP
