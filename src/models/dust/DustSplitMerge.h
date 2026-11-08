/*
 * Dust-specific statistical-particle splitting and merging.
 *
 * AMPS dust grains carry continuous physical radius, per-grain mass and
 * per-grain electric charge in addition to the generic PIC phase-space
 * coordinates.  Generic particle consolidation must therefore not combine
 * arbitrary dust particles: changing radius changes charging, drag,
 * radiation, and sampling physics.  This dependency-light production API
 * defines the admissible dust operation before a PIC list adapter allocates or
 * deletes any particle.  It does not alter generic PIC split/merge behavior.
 */

#ifndef AMPS_SRC_MODELS_DUST_DUST_SPLIT_MERGE_H
#define AMPS_SRC_MODELS_DUST_DUST_SPLIT_MERGE_H

#include "DustCore.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <vector>

namespace ElectricallyChargedDust {
namespace SplitMerge {

// Typed outcomes let a caller distinguish invalid data, cross-owner/cross-
// species requests, size-physics rejection, and an operation that cannot
// reduce particle count.  Rejection never mutates the supplied states.
enum class Status {
  kSuccess = 0,
  kInvalidConfiguration,
  kIndividualWeightUnavailable,
  kInvalidParticle,
  kDifferentOwnershipCell,
  kDifferentCompiledSpecies,
  kDifferentSizeGroup,
  kDifferentIntrinsicGrain,
  kInsufficientParticles,
  kSampledSpeedUnrepresentable
};

inline const char* StatusName(Status status) {
  switch (status) {
    case Status::kSuccess: return "success";
    case Status::kInvalidConfiguration: return "invalid-configuration";
    case Status::kIndividualWeightUnavailable:
      return "individual-weight-unavailable";
    case Status::kInvalidParticle: return "invalid-particle";
    case Status::kDifferentOwnershipCell: return "different-ownership-cell";
    case Status::kDifferentCompiledSpecies: return "different-compiled-species";
    case Status::kDifferentSizeGroup: return "different-size-group";
    case Status::kDifferentIntrinsicGrain: return "different-intrinsic-grain";
    case Status::kInsufficientParticles: return "insufficient-particles";
    case Status::kSampledSpeedUnrepresentable:
      return "sampled-speed-unrepresentable";
  }
  return "unknown";
}

// Size groups are equal-width intervals in ln(a), in metres.  Group i is
// [a_min exp(i Delta), a_min exp((i+1) Delta)); the final group owns a_max.
// This is the same bounded mapping used by dust sampling.  Sharing a group is
// necessary but deliberately not sufficient for merging continuous-radius
// grains: exact intrinsic equality is checked below so a merge cannot silently
// replace one charging/drag law with another.
struct SizeGroups {
  double minimum_radius;
  double maximum_radius;
  int count;
  bool individual_weight_available;

  SizeGroups(double minimum, double maximum, int number,
             bool individual_weight_mode = true)
      : minimum_radius(minimum), maximum_radius(maximum), count(number),
        individual_weight_available(individual_weight_mode) {}
};

// A dependency-light view of one dust macroparticle.  The statistical field
// is always the PIC *individual statistical-weight correction*, never the
// species-global base weight.  A PIC adapter may apply this API only to
// particles with the same base species weight and time-step convention.
// ownership_cell is an opaque stable key supplied by the caller; requiring it
// to match prevents a merge from bypassing native AMR ownership/migration.
struct Particle {
  int compiled_species;
  long long ownership_cell;
  double individual_stat_weight_correction;
  double radius;      // m
  double grain_mass;  // kg per represented physical grain
  double grain_charge;  // C per represented physical grain
  double velocity[3];  // m s^-1

  Particle()
      : compiled_species(-1), ownership_cell(-1),
        individual_stat_weight_correction(0.0), radius(0.0),
        grain_mass(0.0), grain_charge(0.0) {
    velocity[0] = velocity[1] = velocity[2] = 0.0;
  }
};

struct Result {
  Status status;
  std::vector<Particle> particles;

  Result() : status(Status::kInvalidParticle) {}
  bool ok() const { return status == Status::kSuccess; }
};

// Conserved represented quantities.  With individual correction w, these are
// N=sum(w), M=sum(w*m), Q=sum(w*q), P=sum(w*m*v), and
// K=sum(0.5*w*m*|v|^2).  The species-global base weight multiplies every term
// by the same constant and therefore does not affect closure ratios.
struct Budget {
  double represented_number;
  double mass;
  double charge;
  double momentum[3];
  double kinetic_energy;
  Core::WeightedMoments sampled_moments;

  Budget()
      : represented_number(0.0), mass(0.0), charge(0.0),
        kinetic_energy(0.0) {
    momentum[0] = momentum[1] = momentum[2] = 0.0;
  }
};

inline Status ValidateConfiguration(const SizeGroups& groups) {
  // Dust splitting never falls back to an implicit unit weight.  A configured
  // PIC adapter must pass its compile-time individual-weight capability here;
  // an OFF build receives a typed unavailable result before particle access.
  if (!groups.individual_weight_available)
    return Status::kIndividualWeightUnavailable;
  if (!std::isfinite(groups.minimum_radius) ||
      !std::isfinite(groups.maximum_radius) ||
      groups.minimum_radius <= 0.0 ||
      groups.maximum_radius <= groups.minimum_radius || groups.count <= 0)
    return Status::kInvalidConfiguration;
  return Status::kSuccess;
}

inline Status ValidateParticle(const Particle& particle,
                               const SizeGroups& groups) {
  const Status configuration = ValidateConfiguration(groups);
  if (configuration != Status::kSuccess) return configuration;
  if (particle.compiled_species < 0 || particle.ownership_cell < 0 ||
      !std::isfinite(particle.individual_stat_weight_correction) ||
      particle.individual_stat_weight_correction <= 0.0 ||
      !std::isfinite(particle.radius) || particle.radius <= 0.0 ||
      !std::isfinite(particle.grain_mass) || particle.grain_mass <= 0.0 ||
      !std::isfinite(particle.grain_charge))
    return Status::kInvalidParticle;
  for (int component = 0; component != 3; ++component)
    if (!std::isfinite(particle.velocity[component]))
      return Status::kInvalidParticle;
  // Radius outside the declared physical interval is configuration/data
  // mismatch, not a clamped merge candidate.
  if (particle.radius < groups.minimum_radius ||
      particle.radius > groups.maximum_radius)
    return Status::kInvalidParticle;
  return Status::kSuccess;
}

inline int Group(const Particle& particle, const SizeGroups& groups) {
  return Core::LogarithmicBin(particle.radius, groups.minimum_radius,
                              groups.maximum_radius, groups.count);
}

inline Budget EvaluateBudget(const std::vector<Particle>& particles) {
  Budget budget;
  for (std::size_t index = 0; index != particles.size(); ++index) {
    const Particle& particle = particles[index];
    const double weight = particle.individual_stat_weight_correction;
    const double represented_mass = weight * particle.grain_mass;
    double speed_squared = 0.0;
    budget.represented_number += weight;
    budget.mass += represented_mass;
    budget.charge += weight * particle.grain_charge;
    for (int component = 0; component != 3; ++component) {
      budget.momentum[component] +=
          represented_mass * particle.velocity[component];
      speed_squared += particle.velocity[component] *
                       particle.velocity[component];
    }
    budget.kinetic_energy += 0.5 * represented_mass * speed_squared;
    Core::AccumulateWeightedMoments(weight, particle.velocity,
                                    particle.grain_charge,
                                    budget.sampled_moments);
  }
  return budget;
}

// Deterministically split one particle into `part_count` identical phase-space
// representatives.  Only the individual correction changes, from w to w/n.
// Consequently number, mass, charge, momentum, kinetic energy, radius/group,
// and every sampled first moment are conserved to roundoff.  Position and PIC
// list links are intentionally outside this algebraic operation and must be
// copied by the native adapter in the same ownership cell.
inline Result Split(const Particle& particle, int part_count,
                    const SizeGroups& groups) {
  Result result;
  result.status = ValidateParticle(particle, groups);
  if (result.status != Status::kSuccess) return result;
  if (part_count < 2) {
    result.status = Status::kInsufficientParticles;
    return result;
  }

  Particle child = particle;
  child.individual_stat_weight_correction =
      particle.individual_stat_weight_correction / part_count;
  if (!(child.individual_stat_weight_correction > 0.0) ||
      !std::isfinite(child.individual_stat_weight_correction)) {
    result.status = Status::kInvalidParticle;
    return result;
  }
  result.particles.assign(static_cast<std::size_t>(part_count), child);
  result.status = Status::kSuccess;
  return result;
}

inline Status MergeAdmissibility(const std::vector<Particle>& particles,
                                 const SizeGroups& groups) {
  const Status configuration = ValidateConfiguration(groups);
  if (configuration != Status::kSuccess) return configuration;
  // Two outputs are required to conserve an arbitrary input's vector momentum
  // and scalar kinetic energy.  At least three inputs are therefore necessary
  // for this operation to reduce particle count.
  if (particles.size() < 3) return Status::kInsufficientParticles;

  for (std::size_t index = 0; index != particles.size(); ++index) {
    const Status valid = ValidateParticle(particles[index], groups);
    if (valid != Status::kSuccess) return valid;
  }

  const Particle& reference = particles[0];
  const int group = Group(reference, groups);
  for (std::size_t index = 1; index != particles.size(); ++index) {
    const Particle& candidate = particles[index];
    if (candidate.ownership_cell != reference.ownership_cell)
      return Status::kDifferentOwnershipCell;
    if (candidate.compiled_species != reference.compiled_species)
      return Status::kDifferentCompiledSpecies;
    if (Group(candidate, groups) != group)
      return Status::kDifferentSizeGroup;
    // Continuous radius controls multiple dust closures.  Exact equality is
    // intentional: accepting merely close values would make mass/charge
    // conservation tolerance-dependent and alter the represented physics.
    if (candidate.radius != reference.radius ||
        candidate.grain_mass != reference.grain_mass ||
        candidate.grain_charge != reference.grain_charge)
      return Status::kDifferentIntrinsicGrain;
  }
  return Status::kSuccess;
}

// Merge N>=3 admissible particles into two representatives.  Both outputs
// carry W/2 and the common intrinsic grain state.  Their velocities are
// u +/- d e_j, where u is the weight-mean velocity, d^2 is the total weighted
// velocity variance, and j is the component with the largest variance
// (x-before-y-before-z tie break).  Equal output weights make momentum exact;
// |d| preserves total scalar kinetic energy.  The pressure tensor and higher
// velocity moments are not generally preserved, an explicit validity limit of
// the two-particle representation.  The sampled raw-speed moment
// sum(w*|v|) is likewise not generally conserved; number, weighted velocity,
// charge and current are exact.  Setting
// require_sampled_speed_conservation=true transactionally returns
// kSampledSpeedUnrepresentable when the candidate cannot also close that raw
// moment, allowing a sampling-sensitive caller to retain the input population.
inline Result Merge(const std::vector<Particle>& particles,
                    const SizeGroups& groups,
                    bool require_sampled_speed_conservation = false) {
  Result result;
  result.status = MergeAdmissibility(particles, groups);
  if (result.status != Status::kSuccess) return result;

  double total_weight = 0.0;
  double mean_velocity[3] = {0.0, 0.0, 0.0};
  for (std::size_t index = 0; index != particles.size(); ++index) {
    const double weight = particles[index].individual_stat_weight_correction;
    total_weight += weight;
    for (int component = 0; component != 3; ++component)
      mean_velocity[component] += weight * particles[index].velocity[component];
  }
  if (!(total_weight > 0.0) || !std::isfinite(total_weight)) {
    result.status = Status::kInvalidParticle;
    return result;
  }
  for (int component = 0; component != 3; ++component)
    mean_velocity[component] /= total_weight;

  double component_variance[3] = {0.0, 0.0, 0.0};
  for (std::size_t index = 0; index != particles.size(); ++index) {
    const double fraction =
        particles[index].individual_stat_weight_correction / total_weight;
    for (int component = 0; component != 3; ++component) {
      const double residual =
          particles[index].velocity[component] - mean_velocity[component];
      component_variance[component] += fraction * residual * residual;
    }
  }
  const double total_variance = component_variance[0] +
      component_variance[1] + component_variance[2];
  if (!std::isfinite(total_variance) || total_variance < 0.0) {
    result.status = Status::kInvalidParticle;
    return result;
  }
  int axis = 0;
  if (component_variance[1] > component_variance[axis]) axis = 1;
  if (component_variance[2] > component_variance[axis]) axis = 2;

  Particle lower = particles[0];
  Particle upper = particles[0];
  lower.individual_stat_weight_correction = 0.5 * total_weight;
  upper.individual_stat_weight_correction = 0.5 * total_weight;
  for (int component = 0; component != 3; ++component) {
    lower.velocity[component] = mean_velocity[component];
    upper.velocity[component] = mean_velocity[component];
  }
  const double offset = std::sqrt(total_variance);
  lower.velocity[axis] -= offset;
  upper.velocity[axis] += offset;

  result.particles.push_back(lower);
  result.particles.push_back(upper);
  if (require_sampled_speed_conservation) {
    const double input_speed = EvaluateBudget(particles).sampled_moments.speed;
    const double output_speed =
        EvaluateBudget(result.particles).sampled_moments.speed;
    const double scale = std::max(1.0, std::fabs(input_speed));
    if (std::fabs(output_speed - input_speed) >
        64.0 * std::numeric_limits<double>::epsilon() * scale) {
      // No list mutation has occurred: discard the candidate representation
      // and let a sampling-sensitive caller retain the original population.
      result.particles.clear();
      result.status = Status::kSampledSpeedUnrepresentable;
      return result;
    }
  }
  result.status = Status::kSuccess;
  return result;
}

}  // namespace SplitMerge
}  // namespace ElectricallyChargedDust

#endif  // AMPS_SRC_MODELS_DUST_DUST_SPLIT_MERGE_H
