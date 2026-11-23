#include "swcme3d_spherical_material_adapter.h"

#include "spherical_shock_only.h"

namespace SEP::SW3DSphericalAdapter {
namespace {

namespace M = Background::Material;
namespace S = CoronaSwcme::SphericalShockOnly;

CoronalCME::Vec3 CoreVector(const M::Vector& value) {
  return {value[0], value[1], value[2]};
}

M::Vector SepVector(const CoronalCME::Vec3& value) {
  return {value.x, value.y, value.z};
}

class SphericalOwner final : public M::Owner {
 public:
  SphericalOwner(std::shared_ptr<const S::Event> event, const S::Numerics& numerics,
                 M::Descriptor descriptor)
      : event_(std::move(event)), numerics_(numerics), descriptor_(std::move(descriptor)) {}

  const M::Descriptor& Identity() const override { return descriptor_; }

  Core::Result<M::State> Sample(const M::Query& query) const override {
    using R = Core::Result<M::State>;
    if (query.side != M::Side::Smooth || !query.interfaceIdentity.empty())
      return R::Failure(Core::StatusCode::UnsupportedCapability,
                        "use the spherical front query for a one-sided state");
    const auto sample = event_->Evaluate(CoreVector(query.positionM), query.epochS);
    if (!sample.ok()) return R::Failure(sample.status.code, sample.status.message);
    // Evaluation remains available behind the front for core ambient
    // diagnostics. This detached transport interface must never mistake the
    // undisturbed reference there for evolved downstream plasma.
    if (sample.value.side != S::FrontSide::Upstream)
      return R::Failure(Core::StatusCode::UnsupportedCapability,
                        "smooth spherical transport state is upstream of the prescribed front");
    auto state = Ambient(sample.value.ambient, query.positionM, query.epochS);
    if (state.region != query.region)
      return R::Failure(Core::StatusCode::OutOfDomain, "query region differs from the canonical open field branch");
    return R::Success(std::move(state));
  }

  Core::Result<M::Motion> Move(const M::Vector& point, double beginS, double endS) const override {
    using R = Core::Result<M::Motion>;
    const auto moved = event_->Move(CoreVector(point), beginS, endS, numerics_);
    if (!moved.ok()) return R::Failure(moved.status.code, moved.status.message);
    M::Motion result;
    result.positionM = SepVector(moved.value.positionM);
    result.velocityMPerS = SepVector(moved.value.velocityMPerS);
    // The shared consumer currently checks position/velocity only. The core
    // retains F and J for conservation audits. They are never used to replace
    // the consumer's general deformation or downstream-volume contract.
    return R::Success(std::move(result));
  }

  Core::Result<M::Surface> Shock(const M::Vector& point, double epochS) const override {
    using R = Core::Result<M::Surface>;
    const auto shock = event_->Surface(CoreVector(point), epochS);
    if (!shock.ok()) return R::Failure(shock.status.code, shock.status.message);
    const auto& value = shock.value;
    // The core uses ray direction. The neutral contract takes a point already
    // on the sphere. Admit only the floating-point projection error here.
    const auto x = SepVector(value.positionM);
    const double residual = std::hypot(x[0]-point[0], x[1]-point[1], x[2]-point[2]);
    if (residual > 32 * std::numeric_limits<double>::epsilon() * value.radiusM)
      return R::Failure(Core::StatusCode::OutOfDomain, "point is off the spherical front");

    M::Surface result;
    result.eventIdentity = descriptor_.eventIdentity;
    result.stableIdentity = result.eventIdentity + ":leading-sphere";
    result.positionM = point;
    result.epochS = epochS;
    result.outwardNormal = SepVector(value.outwardNormal);
    result.normalSpeedMPerS = value.normalSpeedMPerS;
    result.classification = value.jump ? value.jump->branch : "no-fast-shock";
    result.upstream = Ambient(value.ambient.ambient, point, epochS);
    result.upstream->side = M::Side::UpstreamLimit;
    result.upstream->interfaceIdentity = result.stableIdentity;
    if (value.jump) {
      // The local RH primitive provides total pressure, never a downstream
      // species temperature partition or a usable spatial downstream state.
      const auto& jump = value.jump->downstream;
      auto limit = StateIdentity(point, epochS);
      limit.region = M::Region::Sheath;
      limit.side = M::Side::DownstreamLimit;
      limit.interfaceIdentity = result.stableIdentity;
      limit.massDensityKgM3 = jump.massDensityKgM3;
      limit.pressurePa = jump.pressurePa;
      limit.gammaAdiabatic = event_->Configuration().ambient.composition.gammaAdiabatic;
      limit.magneticFieldT = SepVector(jump.magneticFieldT);
      limit.velocityMPerS = SepVector(jump.velocityMPerS);
      result.downstream = std::move(limit);
    }
    return R::Success(std::move(result));
  }

 private:
  M::State StateIdentity(const M::Vector& point, double epochS) const {
    M::State result;
    result.eventIdentity = descriptor_.eventIdentity;
    result.physicalConfigurationIdentity = descriptor_.physicalConfigurationIdentity;
    result.backgroundGeneration = descriptor_.backgroundGeneration;
    result.geometryGeneration = descriptor_.geometryGeneration;
    result.positionM = point;
    result.epochS = epochS;
    return result;
  }

  M::State Ambient(const CoronaSwcme::AmbientState& sample, const M::Vector& point,
                   double epochS) const {
    const auto& p = sample.primitive;
    const auto& c = event_->Configuration().ambient.composition;
    auto result = StateIdentity(point, epochS);
    result.region = p.region == CoronaSwcme::AmbientRegion::PfssOpen ?
        M::Region::CoronalOpen : M::Region::ParkerExterior;
    result.magneticFieldT = SepVector(p.magneticFieldT);
    result.velocityMPerS = SepVector(p.velocityMPerS);
    result.massDensityKgM3 = p.plasma.massDensityKgM3;
    result.pressurePa = p.plasma.pressurePa;
    result.gammaAdiabatic = c.gammaAdiabatic;
    // Species densities are number per cubic metre. Charge neutrality
    // determines the electron density and the configured alpha abundance
    // is relative to protons. Use the local solved alpha temperature when
    // alpha particles are present, since their pressure law can differ
    // from both the electron and proton laws after the thermal transition.
    result.species = std::vector<M::Species>{
        {"electron", p.plasma.electronNumberDensityM3, p.electronTemperatureK},
        {"proton", p.plasma.protonNumberDensityM3, p.protonTemperatureK}};
    if (c.alphaToProtonNumberRatio > 0)
      result.species->push_back({"alpha", c.alphaToProtonNumberRatio * p.plasma.protonNumberDensityM3,
                                 p.alphaTemperatureK});
    if (sample.derivativesValid) {
      // The core stencils preserve topology and polarity. Include both in the
      // branch identity so opposite sectors cannot be merged by a consumer.
      const std::string branch = std::string(CoronaSwcme::Name(p.region)) + ":" +
                                 std::to_string(p.magneticSector);
      result.velocityGradientPerS = M::Gradient{sample.gradientU, branch, true};
      result.magneticGradientTPerM = M::Gradient{sample.gradientB, branch, true};
    }
    return result;
  }

  const std::shared_ptr<const S::Event> event_;
  const S::Numerics numerics_;
  const M::Descriptor descriptor_;
};

}  // namespace

Core::Result<std::shared_ptr<const M::Interval>> Capture(
    std::shared_ptr<const S::Event> event, const S::Numerics& numerics,
    double beginS, double endS, std::uint64_t backgroundGeneration,
    std::uint64_t geometryGeneration) {
  using R = Core::Result<std::shared_ptr<const M::Interval>>;
  if (!event) return R::Failure(Core::StatusCode::InvalidConfiguration, "a spherical core event is required");
  const auto& support = event->Configuration().ambient.support;
  if (!std::isfinite(beginS) || !std::isfinite(endS) || beginS < support.startS ||
      endS > support.endS || endS < beginS)
    return R::Failure(Core::StatusCode::OutOfDomain, "spherical interval lies outside the declared event epochs");
  if (numerics.maximumAttempts == 0 || !std::isfinite(numerics.positionToleranceM) ||
      numerics.positionToleranceM <= 0 || !std::isfinite(numerics.deformationTolerance) ||
      numerics.deformationTolerance <= 0 || !std::isfinite(numerics.maximumTimeStepS) ||
      numerics.maximumTimeStepS <= 0 || !std::isfinite(numerics.maximumArcStepM) ||
      numerics.maximumArcStepM <= 0)
    return R::Failure(Core::StatusCode::InvalidConfiguration, "explicit spherical trajectory accuracy and work controls are required");

  M::Descriptor d;
  // Bind the thermal profile and its EOS under a distinct descriptor. The
  // temperature anchors alone do not describe the solved species pressure
  // outside the thermal transition. Both profiles still require separate
  // qualification before this interval can bind a native material step.
  const bool resolved = static_cast<bool>(event->ResolvedWind());
  d.profile = resolved ? M::Profile::Swcme3dSpherePfssParkerThermalV1 :
                         M::Profile::Swcme3dSpherePfssParkerKinematicV1;
  d.revision = 1;
  d.capabilities = M::MaterialRequirements | M::Bit(M::Capability::ShockSurface);
  d.backgroundGeneration = backgroundGeneration;
  d.geometryGeneration = geometryGeneration;
  d.geometryGenerationPresent = true;
  d.eventIdentity = event->Identity();
  d.physicalConfigurationIdentity = d.eventIdentity;
  d.materialIdentity = d.eventIdentity + ":open-tubes";
  d.compositionIdentity = d.eventIdentity + ":explicit-e-p-alpha";
  d.eosIdentity = d.eventIdentity + (resolved ? ":isopoly-ambient-adiabatic-jump" :
                                               ":isothermal-ambient-adiabatic-jump");
  d.frame = "HCI";
  d.units = "SI";
  d.beginS = beginS;
  d.endS = endS;
  d.densityBasis = M::DensityBasis::CoreNonrelativisticFluid;
  return M::Interval::Create(std::make_shared<const SphericalOwner>(std::move(event), numerics, std::move(d)));
}

}  // namespace SEP::SW3DSphericalAdapter
