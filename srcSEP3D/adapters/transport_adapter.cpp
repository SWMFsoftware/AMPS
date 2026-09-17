#include "transport_adapter.h"

#include <algorithm>
#include <cmath>
#include <limits>

namespace SEP3D {
namespace Adapters {
namespace {

bool Finite(const Core::Vec3& value) {
  return std::isfinite(value.x) && std::isfinite(value.y) &&
         std::isfinite(value.z);
}

Core::Status Invalid(const std::string& message) {
  return Core::Status(Core::StatusCode::InvalidInput, message);
}

MoverResult Failed(const MoverInput& input, const Core::Status& status) {
  MoverResult result;
  result.particle = input.particle;
  result.disposition = ParticleDisposition::Failed;
  result.status = status;
  return result;
}

}  // namespace

const char* Name(ParticleDisposition disposition) {
  switch (disposition) {
    case ParticleDisposition::Active: return "active";
    case ParticleDisposition::Escaped: return "escaped";
    case ParticleDisposition::Absorbed: return "absorbed";
    case ParticleDisposition::Failed: return "failed";
  }
  return "unknown";
}

ShockIntersection FirstShockIntersection(
    const Core::Vec3& initialPositionM,
    const Core::Vec3& finalPositionM,
    double dtS,
    const ExpandingSphericalShock& shock,
    std::uint64_t lastShockGeneration) {
  ShockIntersection result;
  result.generation = shock.generation;
  if (!shock.active || shock.generation == 0 ||
      shock.generation == lastShockGeneration) {
    result.status = Core::Status::OK();
    return result;
  }
  if (!Finite(initialPositionM) || !Finite(finalPositionM) ||
      !Finite(shock.centerM) || !std::isfinite(dtS) || dtS <= 0.0 ||
      !std::isfinite(shock.radiusAtStepStartM) ||
      shock.radiusAtStepStartM <= 0.0 ||
      !std::isfinite(shock.radialSpeedMPerS)) {
    result.status = Invalid("expanding-shock intersection input is invalid");
    return result;
  }

  const Core::Vec3 relative = initialPositionM - shock.centerM;
  const Core::Vec3 segment = finalPositionM - initialPositionM;
  const double radiusChange = shock.radialSpeedMPerS * dtS;
  // |r+u*d|^2=(R+u*dR)^2 produces a quadratic in the substep fraction u.
  const double a = segment.Dot(segment) - radiusChange * radiusChange;
  const double b = 2.0 * (relative.Dot(segment) -
                          shock.radiusAtStepStartM * radiusChange);
  const double c = relative.Dot(relative) -
                   shock.radiusAtStepStartM * shock.radiusAtStepStartM;
  const double scale = std::max(1.0, std::fabs(b) + std::fabs(c));
  double root = std::numeric_limits<double>::infinity();
  if (std::fabs(a) <= 16.0 * std::numeric_limits<double>::epsilon() * scale) {
    if (b != 0.0) root = -c / b;
  } else {
    const double discriminant = b * b - 4.0 * a * c;
    if (discriminant >= 0.0) {
      const double squareRoot = std::sqrt(discriminant);
      const double first = (-b - squareRoot) / (2.0 * a);
      const double second = (-b + squareRoot) / (2.0 * a);
      if (first >= 0.0 && first <= 1.0) root = first;
      if (second >= 0.0 && second <= 1.0) root = std::min(root, second);
    }
  }
  if (std::isfinite(root)) {
    result.crossed = true;
    result.stepFraction = root;
    result.positionM = initialPositionM + segment * root;
  }
  result.status = Core::Status::OK();
  return result;
}

const std::vector<std::string>& ProductionMoverRegistry::CanonicalNames() {
  static const std::vector<std::string> names = {
      "parker3d-tensor", "focused3d-split"};
  return names;
}

Core::Status ProductionMoverRegistry::Validate(
    RuntimeModel::TransportModel selected) {
  switch (selected) {
    case RuntimeModel::TransportModel::Parker3D:
    case RuntimeModel::TransportModel::Focused3D:
      return Core::Status::OK();
  }
  return Invalid("transport selection is not in the production mover registry");
}

MoverResult AdvanceParticle(const MoverInput& input) {
  const Core::Status registered = ProductionMoverRegistry::Validate(input.model);
  if (!registered.ok()) return Failed(input, registered);
  const ParticleRecord& particle = input.particle;
  const Background::BackgroundSample& background = input.local.background;
  if (particle.stableId == 0 || particle.species < 0 ||
      !Finite(particle.positionM) ||
      !std::isfinite(particle.momentumKgMPerS) ||
      particle.momentumKgMPerS < 0.0 || !std::isfinite(particle.mu) ||
      particle.mu < -1.0 || particle.mu > 1.0 ||
      !std::isfinite(particle.statisticalWeight) ||
      particle.statisticalWeight <= 0.0 ||
      !std::isfinite(input.speciesMassKg) || input.speciesMassKg <= 0.0 ||
      !std::isfinite(input.requestedDtS) || input.requestedDtS <= 0.0 ||
      !std::isfinite(input.innerRadiusM) || input.innerRadiusM <= 0.0 ||
      !std::isfinite(input.outerRadiusM) ||
      input.outerRadiusM <= input.innerRadiusM || input.campaignSeed == 0) {
    return Failed(input, Invalid("particle mover input record is invalid"));
  }
  if (!background.status.ok() || !background.valid ||
      !Finite(background.B) || !Finite(background.U) ||
      !std::isfinite(background.absB) || background.absB <= 0.0 ||
      std::fabs(background.bHat.Norm() - 1.0) > 1.0e-12) {
    return Failed(input, Core::Status(
        Core::StatusCode::BackgroundInvalid,
        "particle mover received an invalid background sample"));
  }

  const double speed = Transport::RelativisticSpeed(
      particle.momentumKgMPerS, input.speciesMassKg);
  Transport::TimeStepPhysics physics;
  physics.requestedS = input.requestedDtS;
  physics.cellSizeM = input.local.cellSizeM;
  physics.characteristicSpeedMPerS = background.U.Norm() + speed;
  physics.kappaParallelM2PerS = input.model == RuntimeModel::TransportModel::Parker3D
      ? input.local.kappaParallelM2PerS : 0.0;
  physics.focusingRatePerS = input.model == RuntimeModel::TransportModel::Focused3D
      ? std::fabs(Transport::FocusedPitchDriftPerS(
            particle.mu, speed, [&]() {
              Transport::FocusedLocalState local;
              local.bHat = background.bHat;
              local.bulkVelocityMPerS = background.U;
              local.divBhatPerM = background.divBhat;
              local.divUPerS = background.divU;
              local.fieldAlignedStrainPerS = background.fieldAlignedStrain;
              return local;
            }())) : 0.0;
  physics.coolingRatePerS = std::fabs(background.divU) / 3.0;
  physics.fractionalFieldVariationPerS =
      input.local.fractionalFieldVariationPerS;
  physics.timeToSnapshotBoundaryS = input.local.timeToSnapshotBoundaryS;
  // The exact expanding-sphere intersection is evaluated after a trial move,
  // but the selector still needs a conservative pre-move shock time scale.
  // Divide the shortest radial gap by the largest possible closing speed: the
  // particle speed, plasma advection, and shock expansion.  This bound may
  // subcycle earlier than necessary for a tangential trajectory, but it can
  // never step deliberately over a nearby moving surface.  A zero value means
  // "not applicable" to SelectTimeStep, as for every other named limiter.
  if (input.shock.active && input.shock.radiusAtStepStartM > 0.0) {
    const double gapM = std::fabs(
        (particle.positionM - input.shock.centerM).Norm() -
        input.shock.radiusAtStepStartM);
    const double closingSpeedMPerS = speed + background.U.Norm() +
        std::fabs(input.shock.radialSpeedMPerS);
    if (gapM > 0.0 && closingSpeedMPerS > 0.0)
      physics.timeToShockCrossingS = gapM / closingSpeedMPerS;
  }
  MoverResult result;
  result.particle = particle;
  result.selectedStep = Transport::SelectTimeStep(input.timeStepControls, physics);
  if (!result.selectedStep.status.ok())
    return Failed(input, result.selectedStep.status);

  Transport::RandomKey key;
  key.campaignSeed = input.campaignSeed;
  key.particleId = particle.stableId;
  key.step = particle.completedStep;
  key.substep = particle.substep;
  const double dtS = result.selectedStep.valueS;
  if (input.model == RuntimeModel::TransportModel::Parker3D) {
    key.purpose = Transport::RandomPurpose::ParkerParallel;
    Transport::KeyedRandomStream random(key);
    Transport::ParkerParticleState state;
    state.positionM = particle.positionM;
    state.momentumKgMPerS = particle.momentumKgMPerS;
    Transport::ParkerLocalState local;
    local.bulkVelocityMPerS = background.U;
    local.bHat = background.bHat;
    local.curvaturePerM = background.curvature;
    local.divBhatPerM = background.divBhat;
    local.divUPerS = background.divU;
    local.kappaParallelM2PerS = input.local.kappaParallelM2PerS;
    local.dKappaParallelDsMPerS = input.local.dKappaParallelDsMPerS;
    const auto moved = Transport::AdvanceParker(state, local, dtS, &random);
    if (!moved.status.ok()) return Failed(input, moved.status);
    result.particle.positionM = moved.state.positionM;
    result.particle.momentumKgMPerS = moved.state.momentumKgMPerS;
  } else {
    key.purpose = Transport::RandomPurpose::FocusedPitch;
    Transport::KeyedRandomStream random(key);
    Transport::FocusedParticleState state;
    state.positionM = particle.positionM;
    state.momentumKgMPerS = particle.momentumKgMPerS;
    state.mu = particle.mu;
    Transport::FocusedLocalState local;
    local.bulkVelocityMPerS = background.U;
    local.bHat = background.bHat;
    local.divBhatPerM = background.divBhat;
    local.divUPerS = background.divU;
    local.fieldAlignedStrainPerS = background.fieldAlignedStrain;
    local.dMuMuPerS = input.local.dMuMuPerS;
    local.dDmuMuDmuPerS = input.local.dDmuMuDmuPerS;
    local.scheme = input.pitchScheme ==
            RuntimeModel::PitchAngleSchemeMode::ReflectingMilstein
        ? Transport::PitchAngleScheme::ReflectingMilstein
        : Transport::PitchAngleScheme::ReflectingEulerMaruyama;
    const auto moved = Transport::AdvanceFocused(
        state, local, input.speciesMassKg, dtS, &random);
    if (!moved.status.ok()) return Failed(input, moved.status);
    result.particle.positionM = moved.state.positionM;
    result.particle.momentumKgMPerS = moved.state.momentumKgMPerS;
    result.particle.mu = moved.state.mu;
  }

  ++result.particle.substep;
  const double radius = result.particle.positionM.Norm();
  if (radius < input.innerRadiusM) {
    result.disposition = ParticleDisposition::Absorbed;
    result.status = Core::Status(Core::StatusCode::InnerBoundary,
                                 "particle crossed the inner boundary");
    return result;
  }
  if (radius > input.outerRadiusM) {
    result.disposition = ParticleDisposition::Escaped;
    result.status = Core::Status(Core::StatusCode::DomainExit,
                                 "particle crossed the outer boundary");
    return result;
  }

  result.shockIntersection = FirstShockIntersection(
      particle.positionM, result.particle.positionM, dtS,
      input.shock, particle.lastShockGeneration);
  if (!result.shockIntersection.status.ok())
    return Failed(input, result.shockIntersection.status);
  if (result.shockIntersection.crossed)
    result.particle.lastShockGeneration = input.shock.generation;
  result.disposition = ParticleDisposition::Active;
  result.status = Core::Status::OK();
  return result;
}

}  // namespace Adapters
}  // namespace SEP3D
