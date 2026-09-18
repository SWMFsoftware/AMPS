#include "source_runtime.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <utility>

namespace SEP3D {
namespace Adapters {
namespace {

Core::Status Invalid(const std::string& message) {
  return Core::Status(Core::StatusCode::InvalidInput, message);
}

bool Finite(const Core::Vec3& value) {
  return std::isfinite(value.x) && std::isfinite(value.y) &&
         std::isfinite(value.z);
}

double KineticEnergy(double momentum, double mass) {
  const double mc = mass * Core::Const::c;
  return (std::sqrt(momentum * momentum + mc * mc) - mc) * Core::Const::c;
}

std::uint64_t InjectionSequence(std::uint64_t generation,
                                std::uint64_t step) {
  // SplitMix-style reversible mixing gives every (shock generation, tick)
  // pair its own semantic source event.  It is a key constructor, not a
  // random draw, so changing the number or ownership of shock patches cannot
  // perturb any other particle stream.
  std::uint64_t value = generation ^
      (step + UINT64_C(0x9e3779b97f4a7c15) + (generation << 6) +
       (generation >> 2));
  value ^= value >> 30; value *= UINT64_C(0xbf58476d1ce4e5b9);
  value ^= value >> 27; value *= UINT64_C(0x94d049bb133111eb);
  value ^= value >> 31;
  return value == 0 ? 1 : value;
}

}  // namespace

bool ShockState::Covers(double timeS) const {
  return status.ok() && std::isfinite(timeS) &&
      timeS >= epochS && timeS <= validUntilS;
}

AnalyticSphericalShockProvider::AnalyticSphericalShockProvider(
    double activeFromS, double activeUntilS, const Core::Vec3& centerM,
    double initialRadiusM, double speedMPerS, double compressionRatio,
    std::uint64_t generation, const std::string& fingerprint)
    : activeFromS_(activeFromS) {
  base_.active = true;
  base_.generation = generation;
  base_.epochS = activeFromS;
  base_.validUntilS = activeUntilS;
  base_.centerM = centerM;
  base_.radiusM = initialRadiusM;
  base_.radialSpeedMPerS = speedMPerS;
  base_.compressionRatio = compressionRatio;
  base_.providerIdentity = CanonicalName();
  base_.configurationFingerprint = fingerprint;
  base_.status = Core::Status::OK();
  if (!std::isfinite(activeFromS) || !std::isfinite(activeUntilS) ||
      activeUntilS <= activeFromS || !Finite(centerM) ||
      !std::isfinite(initialRadiusM) || initialRadiusM <= 0.0 ||
      !std::isfinite(speedMPerS) || speedMPerS <= 0.0 ||
      !std::isfinite(compressionRatio) || compressionRatio <= 1.0 ||
      generation == 0 || fingerprint.empty()) {
    base_.status = Invalid("analytic shock configuration is invalid");
    base_.active = false;
  }
}

ShockState AnalyticSphericalShockProvider::Evaluate(double timeS) const {
  ShockState result = base_;
  if (!base_.status.ok()) return result;
  if (!std::isfinite(timeS)) {
    result.status = Invalid("analytic shock evaluation time is invalid");
    result.active = false;
    return result;
  }
  result.active = timeS >= base_.epochS && timeS <= base_.validUntilS;
  if (result.active)
    result.radiusM = base_.radiusM +
        base_.radialSpeedMPerS * (timeS - activeFromS_);
  result.epochS = timeS;
  result.status = Core::Status::OK();
  return result;
}

PublishedShockProvider::PublishedShockProvider(std::string identity)
    : identity_(std::move(identity)) {
  if (identity_.empty()) identity_ = "host-published-shock";
}

Core::Status PublishedShockProvider::Publish(const ShockState& state) {
  if (!state.status.ok() || state.generation == 0 ||
      !std::isfinite(state.epochS) || !std::isfinite(state.validUntilS) ||
      state.validUntilS < state.epochS || !Finite(state.centerM) ||
      !std::isfinite(state.radiusM) || state.radiusM <= 0.0 ||
      !std::isfinite(state.radialSpeedMPerS) ||
      !std::isfinite(state.compressionRatio) ||
      state.compressionRatio <= 1.0 ||
      state.configurationFingerprint.empty())
    return Invalid("published shock state is incomplete or invalid");
  if (available_ && state.generation <= published_.generation)
    return Invalid("published shock generation must increase monotonically");
  ShockState candidate = state;
  candidate.providerIdentity = identity_;
  published_ = std::move(candidate);
  available_ = true;
  return Core::Status::OK();
}

ShockState PublishedShockProvider::Evaluate(double timeS) const {
  if (!available_) {
    ShockState missing;
    missing.status = Core::Status(Core::StatusCode::SnapshotUnavailable,
                                  "no shock generation has been published");
    return missing;
  }
  ShockState result = published_;
  if (!published_.Covers(timeS)) {
    result.status = Core::Status(Core::StatusCode::SnapshotUnavailable,
                                "published shock does not cover requested time");
    result.active = false;
  }
  return result;
}

InjectionPlan BuildInjectionPlan(const SourceRequest& request) {
  InjectionPlan result;
  result.ledger.step = request.step;
  result.ledger.species = request.species;
  result.ledger.shockGeneration = request.patch.eventGeneration;
  result.ledger.sourceId = request.patch.sourceId;
  if (!request.patch.active) {
    ++result.ledger.inactivePatches;
    result.status = Core::Status::OK();
    return result;
  }
  if (!request.connected) {
    ++result.ledger.disconnectedPatches;
    result.status = Core::Status::OK();
    return result;
  }
  if (!request.patch.status.ok() || request.species < 0 ||
      !std::isfinite(request.speciesMassKg) || request.speciesMassKg <= 0.0 ||
      !std::isfinite(request.intervalS) || request.intervalS <= 0.0 ||
      !std::isfinite(request.physicalParticleRatePerS) ||
      request.physicalParticleRatePerS < 0.0 ||
      !std::isfinite(request.macroparticleWeight) ||
      request.macroparticleWeight <= 0.0 ||
      request.maximumMacroparticles == 0) {
    ++result.ledger.rejected;
    result.status = Invalid("source request is incomplete or invalid");
    return result;
  }

  const double represented =
      request.physicalParticleRatePerS * request.intervalS;
  result.ledger.representedParticles = represented;
  if (represented == 0.0) {
    result.status = Core::Status::OK();
    return result;
  }
  const double expectedMacro = represented / request.macroparticleWeight;
  if (!std::isfinite(expectedMacro) ||
      expectedMacro > static_cast<double>(UINT64_MAX)) {
    ++result.ledger.rejected;
    result.status = Invalid("source macroparticle expectation overflows");
    return result;
  }
  std::uint64_t count = static_cast<std::uint64_t>(std::floor(expectedMacro));
  const double fractional = expectedMacro - static_cast<double>(count);
  Transport::RandomKey roundingKey;
  roundingKey.campaignSeed = request.patch.injection.campaignSeed;
  roundingKey.particleId = request.patch.sourceId;
  roundingKey.step = request.step;
  roundingKey.substep = request.patch.eventGeneration;
  roundingKey.purpose = Transport::RandomPurpose::SourceCount;
  Transport::KeyedRandomStream rounding(roundingKey);
  if (fractional > 0.0 && rounding.UniformOpen01() < fractional) ++count;
  if (count == 0) {
    // One weighted macro represents a non-zero source without silently losing
    // the event.  Its statistical weight is the exact represented number.
    count = 1;
  }
  if (count > request.maximumMacroparticles) {
    result.ledger.capped = count - request.maximumMacroparticles;
    count = request.maximumMacroparticles;
  }

  ShockSourceRecord sampledSource = request.patch;
  sampledSource.injectionSequence = InjectionSequence(
      request.patch.eventGeneration, request.step);
  sampledSource.injection.macroparticlesPerEvent = count;
  result.particles.reserve(static_cast<std::size_t>(count));
  const double exactWeight = represented / static_cast<double>(count);
  for (std::uint64_t index = 0; index < count; ++index) {
    InjectedParticle particle = SampleInjectedParticle(
        sampledSource, index, request.species);
    if (!particle.status.ok()) {
      ++result.ledger.rejected;
      result.particles.clear();
      result.status = particle.status;
      return result;
    }
    particle.particle.statisticalWeight = exactWeight;
    // A particle injected after completion of tick request.step must use the
    // following tick for its first transport key.  Restored and uninterrupted
    // runs then consume the identical (particle,step,substep,purpose) stream.
    particle.particle.completedStep = request.step;
    const double energy = KineticEnergy(
        particle.particle.momentumKgMPerS, request.speciesMassKg);
    result.ledger.injectedEnergyJ += exactWeight * energy;
    result.ledger.injectedMomentumKgMPerS +=
        request.patch.outwardNormal *
        (exactWeight * particle.particle.momentumKgMPerS * particle.particle.mu);
    result.particles.push_back(std::move(particle));
  }
  result.ledger.macroparticles = count;
  result.status = Core::Status::OK();
  return result;
}

}  // namespace Adapters
}  // namespace SEP3D
