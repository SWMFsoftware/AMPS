#include "sep_corona_swcme/coronal_event_ambient.h"

#include <cmath>

#include "sep_coronal_cme/constants.h"

namespace SEP {
namespace CoronaSwcme {
namespace {

class CanonicalAmbientEpoch final : public CoronalCME::AmbientEpoch {
 public:
  CanonicalAmbientEpoch(std::shared_ptr<const AmbientModel> ambient,
                        CoronalCME::AmbientEpochIdentity identity)
      : ambient_(std::move(ambient)), identity_(std::move(identity)) {}

  const CoronalCME::AmbientEpochIdentity& Identity() const override {
    return identity_;
  }

  Core::Result<CoronalCME::EventAmbientSample> At(
      CoronalCME::Vec3 positionM, CoronalCME::QuerySide requested) const override {
    using Return = Core::Result<CoronalCME::EventAmbientSample>;
    using Side = CoronalCME::QuerySide;
    if (requested != Side::Ordinary && requested != Side::Minus &&
        requested != Side::Plus)
      return Return::Failure(Core::StatusCode::InvalidConfiguration,
                             "unknown canonical ambient radial trace selector");

    // A physical PFSS/Parker radial trace is independent of whether a
    // Cartesian finite-difference stencil happened to be central/one-sided.
    // At exactly Rss the existing canonical evaluator selects Parker. Its
    // opposite trace is not obtained by relabeling that state or shifting the
    // point by an undeclared epsilon: that unsupported request must fail.
    const double radiusM = CoronalCME::Norm(positionM);
    const double sourceRadiusM = ambient_->Definition().ambient.sourceSurfaceRadiusM;
    const Side actual = radiusM < sourceRadiusM ? Side::Minus : Side::Plus;
    if (requested != Side::Ordinary && requested != actual)
      return Return::Failure(
          Core::StatusCode::OutOfDomain,
          "requested PFSS/Parker trace disagrees with the actual radial branch");

    const auto queried = ambient_->EvaluateWithDerivatives(positionM, identity_.timeS,
                                                           identity_.generation);
    if (!queried.ok())
      return Return::Failure(queried.status.code, queried.status.message);
    const auto& state = queried.value;
    if (state.eventIdentity != identity_.physicsIdentity ||
        state.epochS != identity_.timeS || state.generation != identity_.generation ||
        !state.derivativesValid ||
        !ValidAmbientRegionSector(state.primitive.region,
                                  state.primitive.magneticSector))
      return Return::Failure(
          Core::StatusCode::DataIntegrityFailure,
          "canonical ambient returned stale epoch/identity or invalid sector/tensors");

    const auto& primitive = state.primitive;
    const auto& composition = ambient_->Definition().composition;
    CoronalCME::EventAmbientSample result;
    result.primitive = {primitive.plasma.massDensityKgM3, primitive.plasma.pressurePa,
                        primitive.velocityMPerS, primitive.magneticFieldT};
    result.electronDensityM3 = primitive.plasma.electronNumberDensityM3;
    result.electronTemperatureK = primitive.electronTemperatureK;
    result.ions = {{"H+", 1, 1, CoronalCME::Constants::kProtonMassKg,
                    primitive.protonTemperatureK}};
    if (composition.alphaToProtonNumberRatio > 0)
      result.ions.push_back({"He++", composition.alphaToProtonNumberRatio, 2,
                             CoronalCME::Constants::kAlphaMassKg,
                             composition.alphaTemperatureK});
    result.region = primitive.region == AmbientRegion::PfssOpen
                        ? CoronalCME::AmbientRegion::PfssOpen
                    : primitive.region == AmbientRegion::PfssClosed
                        ? CoronalCME::AmbientRegion::PfssClosed
                        : CoronalCME::AmbientRegion::ParkerExterior;
    result.topology = primitive.region == AmbientRegion::PfssClosed
                          ? CoronalCME::FieldLineTopology::ClosedBelowOuterBoundary
                          : CoronalCME::FieldLineTopology::OpenToOuterBoundary;
    result.sector = primitive.magneticSector;
    result.gradientBPerM = state.gradientB;
    result.gradientUPerS = state.gradientU;
    result.derivativesValid = true;
    result.oneSided = requested != Side::Ordinary;
    result.side = result.oneSided ? actual : Side::Ordinary;
    for (std::size_t dimension = 0; dimension < 3; ++dimension) {
      const auto stencil = state.stencils[dimension];
      result.derivativeStencils[dimension] =
          stencil == DifferenceStencil::CentralSecondOrder
              ? CoronalCME::EventDifferenceStencil::CentralSecondOrder
          : stencil == DifferenceStencil::ForwardSecondOrder
              ? CoronalCME::EventDifferenceStencil::ForwardSecondOrder
              : CoronalCME::EventDifferenceStencil::BackwardSecondOrder;
    }
    return Return::Success(std::move(result));
  }

 private:
  std::shared_ptr<const AmbientModel> ambient_;
  CoronalCME::AmbientEpochIdentity identity_;
};
}  // namespace

Core::Result<std::shared_ptr<const CoronalCME::AmbientEpoch>> BindAmbientEpoch(
    std::shared_ptr<const AmbientModel> ambient, double timeS, std::uint64_t generation,
    const std::string& epochUtc, const std::string& expectedAmbientIdentity) {
  using Return = Core::Result<std::shared_ptr<const CoronalCME::AmbientEpoch>>;
  if (!ambient || !std::isfinite(timeS) || generation == 0 || epochUtc.empty() ||
      expectedAmbientIdentity.empty())
    return Return::Failure(Core::StatusCode::InvalidConfiguration,
                           "canonical ambient binding requires an owner, finite epoch, "
                           "generation and identities");
  const auto& definition = ambient->Definition();
  if (definition.physicsFingerprint != expectedAmbientIdentity)
    return Return::Failure(Core::StatusCode::DataIntegrityFailure,
                           "canonical ambient/event physical identity mismatch");
  if (definition.coordinateFrame != "HCI-like-inertial")
    return Return::Failure(
        Core::StatusCode::UnsupportedCapability,
        "canonical ambient binding supports its declared HCI-like inertial frame only");
  if (timeS < definition.support.startS || timeS > definition.support.endS)
    return Return::Failure(
        Core::StatusCode::OutOfDomain,
        "canonical ambient binding time lies outside the complete event history");
  CoronalCME::AmbientEpochIdentity identity = {
      expectedAmbientIdentity,
      definition.coordinateFrame,
      epochUtc,
      timeS,
      timeS,
      timeS,
      definition.composition.gammaAdiabatic,
      generation,
      definition.composition.includeElectronMass};
  std::shared_ptr<const CoronalCME::AmbientEpoch> result =
      std::make_shared<const CanonicalAmbientEpoch>(std::move(ambient),
                                                    std::move(identity));
  return Return::Success(std::move(result));
}

}  // namespace CoronaSwcme
}  // namespace SEP
