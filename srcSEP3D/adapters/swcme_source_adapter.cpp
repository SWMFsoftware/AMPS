#include "swcme_source_adapter.h"

#include <cmath>
#include <iomanip>
#include <limits>
#include <sstream>

namespace SEP3D {
namespace Adapters {
namespace {

Core::Status Invalid(const std::string& message) {
  return Core::Status(Core::StatusCode::InvalidInput, message);
}

// Produce an application-level status without leaking the SWCME status enum
// through the adapter API. The original symbolic code remains in the message
// so a coupled-run diagnostic can be traced to its provider.
Core::Status Translate(const swcme::ModelStatus& status,
                       const char* operation) {
  if (status.ok()) return Core::Status::OK();
  std::ostringstream out;
  out << operation << " failed with SWCME status "
      << swcme::status_code_name(status.code);
  return Invalid(out.str());
}

Core::Status Translate(const SEP::Transport::Status& status,
                       const char* operation) {
  if (status.ok()) return Core::Status::OK();
  return Invalid(std::string(operation) + ": " + status.message);
}

SEP::Injection::RandomKey Key(const ShockSourceRecord& source,
                              std::uint64_t macroIndex,
                              int species,
                              SEP::Injection::RandomPurpose purpose) {
  SEP::Injection::RandomKey key;
  key.campaign = source.injection.campaignSeed;
  key.event = source.eventGeneration;
  key.fieldLine = source.sourceId;  // semantic slot: shock-patch identity
  key.species = static_cast<std::uint64_t>(species);
  key.macroparticle = macroIndex;
  key.purpose = purpose;
  return key;
}

}  // namespace

ShockSourceRecord MakeShockSourceRecord(
    const swcme::sep::SEPSourceState& source,
    std::uint64_t eventGeneration,
    std::uint64_t campaignSeed,
    std::uint64_t macroparticlesPerEvent,
    double injectionEfficiency) {
  ShockSourceRecord result;
  result.status = Translate(source.status, "SWCME source conversion");
  if (!result.status.ok()) return result;
  if (!source.active) {
    result.status = Invalid("SWCME source is not an active physical shock");
    return result;
  }
  if (eventGeneration == 0 || campaignSeed == 0 ||
      macroparticlesPerEvent == 0 ||
      !std::isfinite(injectionEfficiency) || injectionEfficiency <= 0.0 ||
      injectionEfficiency > 1.0 || !std::isfinite(source.q_phase_space) ||
      source.q_phase_space <= 2.0 ||
      !std::isfinite(source.relative_patch_weight) ||
      source.relative_patch_weight <= 0.0) {
    result.status = Invalid("SWCME injection controls or DSA source are invalid");
    return result;
  }

  result.eventGeneration = eventGeneration;
  result.sourceId = static_cast<std::uint64_t>(source.source_id);
  result.positionM = Core::Vec3(source.position_m.data());
  result.outwardNormal = Core::Vec3(source.normal.data()).Normalized();
  result.relativePatchWeight = source.relative_patch_weight;
  result.compression = source.compression;
  result.shockNormalSpeedMPerS = source.normal_speed_m_s;
  if (result.outwardNormal.NormSq() == 0.0) {
    result.status = Invalid("SWCME shock normal is degenerate");
    return result;
  }

  result.injection.campaignSeed = campaignSeed;
  result.injection.macroparticlesPerEvent = macroparticlesPerEvent;
  result.injection.injectionEfficiency = injectionEfficiency;
  result.injection.angular = SEP::Injection::AngularDistribution::Isotropic;
  result.injection.spectrum.measure = SEP::Injection::Measure::Momentum;
  result.injection.spectrum.minimum =
      swcme::sep::momentum_kg_m_s_from_kinetic_MeV(
          source.spectrum.kinetic_energy_min_MeV,
          source.spectrum.particle_mass_kg);
  result.injection.spectrum.maximum =
      swcme::sep::momentum_kg_m_s_from_kinetic_MeV(
          source.spectrum.kinetic_energy_max_MeV,
          source.spectrum.particle_mass_kg);
  result.injection.spectrum.powerIndex = source.q_phase_space - 2.0;
  const SEP::Transport::Status valid = SEP::Injection::Validate(result.injection);
  result.status = Translate(valid, "shared injection-spectrum validation");
  if (!result.status.ok()) return result;

  // The fingerprint combines the exact common spectrum fingerprint with the
  // provider's generation and source identity. This is diagnostic metadata,
  // not a source of random numbers; keys above remain the stochastic authority.
  std::ostringstream identity;
  identity << SEP::Injection::Fingerprint(result.injection) << ':'
           << eventGeneration << ':' << result.sourceId << ':'
           << std::setprecision(17) << result.relativePatchWeight << ':'
           << source.q_phase_space;
  result.sourceFingerprint = identity.str();
  result.active = true;
  result.status = Core::Status::OK();
  return result;
}

InjectedParticle SampleInjectedParticle(const ShockSourceRecord& source,
                                         std::uint64_t macroIndex,
                                         int species) {
  InjectedParticle result;
  if (!source.status.ok() || !source.active || species < 0 ||
      macroIndex >= source.injection.macroparticlesPerEvent) {
    result.status = Invalid("injection request is outside the active event");
    return result;
  }

  SEP::Transport::KeyedRandomStream momentumRandom =
      SEP::Injection::MakeRandomStream(Key(
          source, macroIndex, species, SEP::Injection::RandomPurpose::Spectrum));
  const SEP::Transport::ScalarResult sampled = SEP::Injection::InverseCdf(
      source.injection.spectrum, momentumRandom.UniformOpen01());
  result.status = Translate(sampled.status, "shared inverse-CDF injection");
  if (!result.status.ok()) return result;

  SEP::Transport::KeyedRandomStream pitchRandom =
      SEP::Injection::MakeRandomStream(Key(
          source, macroIndex, species, SEP::Injection::RandomPurpose::PitchAngle));
  SEP::Transport::KeyedRandomStream phaseRandom =
      SEP::Injection::MakeRandomStream(Key(
          source, macroIndex, species, SEP::Injection::RandomPurpose::Gyrophase));

  result.particle.stableId = SEP::Injection::HashRandomKey(Key(
      source, macroIndex, species, SEP::Injection::RandomPurpose::Position));
  // AMPS reserves zero as "uninitialized" in this adapter. HashRandomKey is
  // effectively never zero, but the explicit remap makes the invariant total.
  if (result.particle.stableId == 0) result.particle.stableId = 1;
  result.particle.species = species;
  result.particle.positionM = source.positionM;
  result.particle.momentumKgMPerS = sampled.value;
  result.particle.mu = 2.0 * pitchRandom.UniformOpen01() - 1.0;
  result.particle.gyrophaseRad =
      2.0 * Core::Const::kPi * phaseRandom.UniformOpen01();
  result.particle.statisticalWeight =
      source.relativePatchWeight * source.injection.injectionEfficiency /
      static_cast<double>(source.injection.macroparticlesPerEvent);
  result.particle.lastShockGeneration = source.eventGeneration;
  result.status = Core::Status::OK();
  return result;
}

}  // namespace Adapters
}  // namespace SEP3D
