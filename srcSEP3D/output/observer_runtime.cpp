#include "observer_runtime.h"

#include <algorithm>
#include <cmath>

namespace SEP3D {
namespace Output {
namespace {

Core::Status Invalid(const std::string& message) {
  return Core::Status(Core::StatusCode::InvalidInput, message);
}

std::vector<double> LogEdges(double minimum, double maximum, unsigned bins) {
  std::vector<double> result;
  result.reserve(static_cast<std::size_t>(bins) + 1);
  const double ratio = std::pow(maximum / minimum, 1.0 / bins);
  for (unsigned index = 0; index <= bins; ++index)
    result.push_back(index == bins ? maximum : minimum * std::pow(ratio, index));
  return result;
}

}  // namespace

Core::Status BuildObserverDefinitions(
    const RuntimeModel::RunConfiguration3D& configuration,
    double simulationTimeS,
    std::vector<VirtualSpacecraftDefinition>* definitions) {
  if (definitions == nullptr || !std::isfinite(simulationTimeS))
    return Invalid("observer-definition output or time is invalid");
  std::vector<VirtualSpacecraftDefinition> candidate;
  for (const RuntimeModel::ObserverOptions& observer :
       configuration.options().observers) {
    VirtualSpacecraftDefinition resolved;
    resolved.name = observer.id;
    resolved.positionM = observer.positionM;
    if (observer.kind == RuntimeModel::ObserverKind::MovingCartesian)
      resolved.positionM += observer.velocityMPerS * simulationTimeS;
    else if (observer.kind == RuntimeModel::ObserverKind::SphericalShell) {
      Core::Vec3 direction = observer.positionM.Normalized();
      if (direction.NormSq() == 0.0) direction = Core::Vec3(1.0, 0.0, 0.0);
      resolved.positionM = direction * observer.shellRadiusM;
    }
    // A coupled field-line service may update observer.positionM before this
    // call.  The resolved sampling contract is still Cartesian and immutable
    // for the complete reduction window.
    resolved.collectionRadiusM = observer.collectionRadiusM;
    resolved.kineticEnergyEdgesJ = LogEdges(
        observer.minimumEnergyJ, observer.maximumEnergyJ,
        observer.energyBins);
    // Sampling defines an empty acceptedSpecies vector as a wildcard.  The
    // configuration factory permits that representation only when the user
    // explicitly selected `species = all`, so no malformed empty list can
    // accidentally broaden an observer here.
    resolved.acceptedSpecies = observer.species;
    resolved.minimumMu = observer.minimumMu;
    resolved.maximumMu = observer.maximumMu;
    resolved.observerKind = RuntimeModel::Name(observer.kind);
    resolved.normalization = RuntimeModel::Name(observer.normalization);
    candidate.push_back(std::move(resolved));
  }
  *definitions = std::move(candidate);
  return Core::Status::OK();
}

Core::Status ObserverRuntime::Capture(SamplingRequest request) {
  if (hasPending_)
    return Core::Status(Core::StatusCode::InvalidTransition,
                        "an observer window is already pending publication");
  request.previousState = state_;
  state_.pendingWindows = 1;
  state_.pendingObservations = request.particles.size();
  state_.pendingRepresentedParticles = 0.0;
  for (const ParticleObservation& particle : request.particles) {
    if (!std::isfinite(particle.statisticalWeight) ||
        particle.statisticalWeight <= 0.0)
      return Invalid("observer capture contains an invalid particle weight");
    state_.pendingRepresentedParticles += particle.statisticalWeight;
  }
  request.previousState = state_;
  pending_ = std::move(request);
  hasPending_ = true;
  return Core::Status::OK();
}

SamplingSnapshot ObserverRuntime::PreparePublication() const {
  if (!hasPending_) {
    SamplingSnapshot missing;
    missing.status = Core::Status(
        Core::StatusCode::InvalidTransition,
        "observer publication requires a captured joined-boundary window");
    return missing;
  }
  return Sample(pending_);
}

Core::Status ObserverRuntime::CommitPublication(
    const SamplingSnapshot& published) {
  if (!hasPending_ || !published.status.ok())
    return Core::Status(Core::StatusCode::InvalidTransition,
                        "only a successful prepared observer snapshot may commit");
  state_ = published.nextState;
  pending_ = SamplingRequest{};
  hasPending_ = false;
  return Core::Status::OK();
}

}  // namespace Output
}  // namespace SEP3D
