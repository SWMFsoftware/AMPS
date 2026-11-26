#include "observer_runtime.h"

#include "observer_trajectory.h"

#include <algorithm>
#include <cmath>
#include <map>

namespace SEP3D {
namespace Output {
namespace {

Core::Status Invalid(const std::string& message) {
  return Core::Status(Core::StatusCode::InvalidInput, message);
}

std::vector<double> EnergyEdges(
    double minimum, double maximum, unsigned bins,
    RuntimeModel::EnergyChannelSpacing spacing) {
  std::vector<double> result;
  result.reserve(static_cast<std::size_t>(bins) + 1);
  if (spacing == RuntimeModel::EnergyChannelSpacing::Linear) {
    const double width = (maximum - minimum) / bins;
    for (unsigned index = 0; index <= bins; ++index)
      result.push_back(index == bins ? maximum : minimum + index * width);
  } else {
    const double ratio = std::pow(maximum / minimum, 1.0 / bins);
    for (unsigned index = 0; index <= bins; ++index)
      result.push_back(index == bins ? maximum
                                    : minimum * std::pow(ratio, index));
  }
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
    if (observer.kind == RuntimeModel::ObserverKind::Body ||
        (observer.kind == RuntimeModel::ObserverKind::MovingCartesian &&
         !configuration.options().ComputationFrameValue().IsInertial())) {
      // Inertial trajectories (an ephemeris body, or a straight HCI line in a
      // corotating run) in the computation frame (observer_trajectory.h).
      // One trajectory per observer and process: a table is read and its
      // SHA-256 checked once.
      //
      // The collection-sphere centre at simulationTimeS is the trajectory's
      // computation-frame position (HCI rotated by -Omega (t - t_ref) in a
      // corotating run) plus the mesh coordinate origin.  It is the same
      // function that traced the observer flux tube (main_lib.cpp), so the
      // sampled sphere follows the path the tube was built around.  A
      // moving-cartesian observer in an HCI run keeps the original branch
      // below (x0 + v t, no frame rotation).
      //
      // The cache key includes the physics fingerprint, so a different
      // configuration in the same process (e.g. the unit-test runner) builds
      // its own trajectory.  The map is not synchronized.  Its callers in
      // main_lib.cpp run on the application's control path, outside
      // threaded regions.
      static std::map<std::string, ObserverTrajectoryPointer> cache;
      ObserverTrajectoryPointer& trajectory =
          cache[configuration.physics_fingerprint() + "|" + observer.id];
      if (!trajectory) {
        const Core::Status created = CreateObserverTrajectory(observer, &trajectory);
        if (!created.ok()) return created;
      }
      const Core::Status placed = ObserverFramePosition(
          *trajectory, configuration.options().ComputationFrameValue(),
          simulationTimeS, &resolved.positionM);
      // A time outside an ephemeris table or the SPICE coverage fails the
      // publication with the observer named; the position is never clamped
      // or extrapolated.
      if (!placed.ok())
        return Core::Status(placed.code, "observer '" + observer.id + "': " +
                                             placed.message);
      // Sun-centred frame position -> mesh coordinates.
      resolved.positionM += configuration.options().coordinateOriginM;
    } else if (observer.kind == RuntimeModel::ObserverKind::MovingCartesian)
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
    resolved.kineticEnergyEdgesJ = EnergyEdges(
        observer.minimumEnergyJ, observer.maximumEnergyJ,
        observer.energyBins, observer.energyChannelSpacing);
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
