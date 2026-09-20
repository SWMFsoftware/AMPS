#include "source_runtime.h"

#include "swcme3d_input.hpp"
#include "swcme_sep_interface.hpp"

#include <algorithm>
#include <cmath>
#include <exception>
#include <limits>
#include <memory>
#include <sstream>
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

class StandaloneSwcmeShockProvider final : public ShockProvider {
 public:
  StandaloneSwcmeShockProvider(
      const swcme::input3d::ResolvedConfiguration& model,
      const RuntimeModel::RunConfiguration3DOptions& application)
      : model_(model), application_(application),
        interface_(new swcme::sep::Interface3D(model.model, model.spectrum)) {}

  const char* CanonicalName() const override {
    return "canonical-swcme3d-standalone";
  }

  ShockState Evaluate(double timeS) const override {
    ShockState result;
    result.providerIdentity = CanonicalName();
    result.configurationFingerprint = model_.fingerprint;
    result.centerM = application_.coordinateOriginM;
    if (!std::isfinite(timeS)) {
      result.status = Invalid("SWCME shock evaluation time is not finite");
      return result;
    }
    if (timeS < model_.valid_from_s || timeS > model_.valid_until_s) {
      result.status = Core::Status::OK();
      result.active = false;
      result.epochS = timeS;
      result.validUntilS = model_.valid_until_s;
      return result;
    }

    try {
      const double modelTimeS = timeS - model_.launch_epoch_s;
      const swcme::sep::Interface3D::PreparedStep step =
          interface_->prepare(modelTimeS);
      if (step.common.apex.status != swcme::kinematics::Status::Ok) {
        result.status = Invalid(
            std::string("SWCME kinematics rejected initialization epoch: ") +
            swcme::kinematics::status_name(step.common.apex.status));
        return result;
      }
      result.epochS = timeS;
      result.validUntilS = std::min(
          model_.valid_until_s, timeS + application_.requestedTimeStepS);
      result.radiusM = step.r_sh_m;
      result.radialSpeedMPerS = step.V_sh_ms;
      result.compressionRatio = step.rc;
      result.active = step.has_shock;
      const long double tick =
          (static_cast<long double>(timeS) - model_.valid_from_s) /
          application_.requestedTimeStepS;
      if (!(tick >= 0.0L) || tick > static_cast<long double>(UINT64_MAX - 1)) {
        result.status = Invalid("SWCME shock generation overflows");
        result.active = false;
        return result;
      }
      result.generation = 1 + static_cast<std::uint64_t>(std::llround(tick));

      swcme::sep::SourceSurface surface;
      const swcme::ModelStatus surfaceStatus =
          interface_->build_shock_surface_source(
              step, model_.surface_theta_intervals,
              model_.surface_phi_points, surface);
      if (!surfaceStatus.ok()) {
        result.status = Invalid(
            std::string("canonical SWCME surface source failed: ") +
            surfaceStatus.summary());
        result.active = false;
        return result;
      }
      result.patches.reserve(surface.active_patch_count);
      for (const swcme::sep::SEPSourceState& source : surface.patches) {
        if (!source.active) continue;
        ShockSourceRecord patch = MakeShockSourceRecord(
            source, result.generation, application_.campaignSeed,
            application_.source.samplesPerStep,
            model_.injection_efficiency);
        if (!patch.status.ok()) {
          result.status = patch.status;
          result.active = false;
          result.patches.clear();
          return result;
        }
        // Canonical SWCME coordinates are heliocentric.  Apply the declared
        // application origin exactly once at the provider boundary.
        patch.positionM += application_.coordinateOriginM;
        result.patches.push_back(std::move(patch));
      }
      result.status = Core::Status::OK();
      return result;
    } catch (const std::exception& exception) {
      result.status = Invalid(
          std::string("canonical SWCME provider exception: ") +
          exception.what());
      result.active = false;
      return result;
    }
  }

 private:
  swcme::input3d::ResolvedConfiguration model_;
  RuntimeModel::RunConfiguration3DOptions application_;
  std::unique_ptr<swcme::sep::Interface3D> interface_;
};

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
    if (request.prescribedMacroparticles != 0) {
      ++result.ledger.rejected;
      result.status = Invalid(
          "cannot prescribe macroparticles for a zero physical source");
      return result;
    }
    result.status = Core::Status::OK();
    return result;
  }
  std::uint64_t count = request.prescribedMacroparticles;
  if (count != 0) {
    // An exact schema-v3 allocation must fail rather than cap: capping would
    // contradict the input contract that exactly samples_per_step are born at
    // every active injection boundary.
    if (count > request.maximumMacroparticles) {
      ++result.ledger.rejected;
      result.status = Invalid(
          "prescribed source count exceeds maximumMacroparticles");
      return result;
    }
  } else {
    const double expectedMacro = represented / request.macroparticleWeight;
    if (!std::isfinite(expectedMacro) ||
        expectedMacro > static_cast<double>(UINT64_MAX)) {
      ++result.ledger.rejected;
      result.status = Invalid("source macroparticle expectation overflows");
      return result;
    }
    count = static_cast<std::uint64_t>(std::floor(expectedMacro));
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
      // One weighted macro represents a non-zero source without silently
      // losing the event.  Its statistical weight is the exact represented
      // number, so this legacy fallback remains conservative.
      count = 1;
    }
    if (count > request.maximumMacroparticles) {
      result.ledger.capped = count - request.maximumMacroparticles;
      count = request.maximumMacroparticles;
    }
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

Core::Status AllocateExactPatchMacroparticles(
    const std::vector<ShockSourceRecord>& patches,
    std::uint64_t totalMacroparticles,
    std::vector<std::uint64_t>* perPatchCounts) {
  if (perPatchCounts == nullptr)
    return Invalid("exact patch-allocation output is null");
  perPatchCounts->clear();
  if (patches.empty()) {
    if (totalMacroparticles == 0) return Core::Status::OK();
    return Invalid("cannot allocate source samples without active patches");
  }
  if (totalMacroparticles < patches.size()) {
    return Invalid(
        "samples_per_step is smaller than the active SWCME patch count; "
        "at least one weighted representative per physical patch is required");
  }

  long double weightSum = 0.0L;
  for (const ShockSourceRecord& patch : patches) {
    if (!patch.status.ok() || !patch.active ||
        !std::isfinite(patch.relativePatchWeight) ||
        patch.relativePatchWeight <= 0.0) {
      return Invalid(
          "exact source allocation requires active patches with positive "
          "finite physical weights");
    }
    weightSum += static_cast<long double>(patch.relativePatchWeight);
  }
  if (!(weightSum > 0.0L) || !std::isfinite(weightSum))
    return Invalid("active SWCME patch-weight sum is invalid");

  perPatchCounts->assign(patches.size(), UINT64_C(1));
  const std::uint64_t remaining =
      totalMacroparticles - static_cast<std::uint64_t>(patches.size());
  if (remaining == 0) return Core::Status::OK();

  struct Remainder {
    long double fraction = 0.0L;
    std::uint64_t sourceId = 0;
    std::size_t index = 0;
  };
  std::vector<Remainder> remainders;
  remainders.reserve(patches.size());
  std::uint64_t apportioned = 0;
  for (std::size_t i = 0; i < patches.size(); ++i) {
    const long double quota = static_cast<long double>(remaining) *
        static_cast<long double>(patches[i].relativePatchWeight) / weightSum;
    if (!(quota >= 0.0L) || !std::isfinite(quota) ||
        quota > static_cast<long double>(UINT64_MAX))
      return Invalid("exact source allocation quota is invalid");
    const std::uint64_t base = static_cast<std::uint64_t>(std::floor(quota));
    if (base > remaining - apportioned)
      return Invalid("exact source allocation lost numerical normalization");
    (*perPatchCounts)[i] += base;
    apportioned += base;
    remainders.push_back(
        {quota - static_cast<long double>(base), patches[i].sourceId, i});
  }
  const std::uint64_t residual = remaining - apportioned;
  if (residual > remainders.size())
    return Invalid("largest-remainder source allocation is inconsistent");
  std::sort(remainders.begin(), remainders.end(),
            [](const Remainder& left, const Remainder& right) {
              if (left.fraction != right.fraction)
                return left.fraction > right.fraction;
              if (left.sourceId != right.sourceId)
                return left.sourceId < right.sourceId;
              return left.index < right.index;
            });
  for (std::uint64_t i = 0; i < residual; ++i)
    ++(*perPatchCounts)[remainders[static_cast<std::size_t>(i)].index];
  return Core::Status::OK();
}

Core::Status CreateStandaloneSwcmeShockProvider(
    const RuntimeModel::RunConfiguration3D& configuration,
    std::shared_ptr<ShockProvider>* provider) {
  if (provider == nullptr) return Invalid("SWCME provider output is null");
  const RuntimeModel::RunConfiguration3DOptions& options =
      configuration.options();
  if (options.inputSchemaVersion < 3 ||
      options.shock != RuntimeModel::ShockAuthority::Swcme ||
      options.swcmeAssignments.empty())
    return Invalid("standalone SWCME provider requires schema version 3 and shock authority swcme");

  std::vector<swcme::input3d::Assignment> assignments;
  assignments.reserve(options.swcmeAssignments.size());
  for (const RuntimeModel::SwcmeAssignment& raw : options.swcmeAssignments) {
    swcme::input3d::Assignment assignment;
    assignment.key = raw.key;
    assignment.value = raw.value;
    assignment.origin = "frozen srcSEP3D configuration";
    assignment.line = raw.line;
    assignments.push_back(assignment);
  }
  const swcme::input3d::ResolveResult resolved =
      swcme::input3d::Resolve(assignments);
  if (!resolved.ok()) {
    std::ostringstream message;
    message << "SWCME3D resolution failed key='" << resolved.status.key << "'";
    if (!resolved.status.message.empty())
      message << ": " << resolved.status.message;
    return Invalid(message.str());
  }
  if (resolved.configuration.fingerprint !=
          options.swcmeConfigurationFingerprint ||
      resolved.configuration.normalized_manifest !=
          options.swcmeResolvedManifest)
    return Core::Status(Core::StatusCode::ConfigurationConflict,
                        "re-resolved SWCME configuration differs from frozen identity");
  std::shared_ptr<ShockProvider> candidate;
  try {
    candidate.reset(new StandaloneSwcmeShockProvider(
        resolved.configuration, options));
  } catch (const std::exception& exception) {
    return Invalid(std::string("cannot construct canonical SWCME provider: ") +
                   exception.what());
  }
  // Preflight the first declared active epoch, rather than merely time zero.
  // A delayed event is legitimately inactive at t=0; accepting that trivial
  // state would postpone a malformed MHD jump or surface mesh until after the
  // expensive AMPS mesh had already been allocated.
  const ShockState initial =
      candidate->Evaluate(resolved.configuration.valid_from_s);
  if (!initial.status.ok()) return initial.status;
  if (!initial.active || initial.patches.empty())
    return Invalid("canonical SWCME source has no active patches at event.valid_from");
  std::vector<std::uint64_t> preflightCounts;
  const Core::Status allocation = AllocateExactPatchMacroparticles(
      initial.patches, options.source.samplesPerStep, &preflightCounts);
  if (!allocation.ok()) return allocation;
  *provider = std::move(candidate);
  return Core::Status::OK();
}

}  // namespace Adapters
}  // namespace SEP3D
