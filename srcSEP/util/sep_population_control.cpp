#include "sep_population_control.h"

#include <algorithm>
#include <cmath>
#include <limits>

namespace SEP {
namespace PopulationControl {
namespace {

Transport::Status Error(const std::string& message) {
  return Transport::Status::Error(Transport::StatusCode::InvalidArgument,
                                  message);
}

void HashWord(std::uint64_t tag, std::uint64_t value, std::uint64_t* hash) {
  // FNV-1a is used for stable identity rather than cryptographic security.
  // Hashing a tag before every value prevents commutative field swaps from
  // producing the same lineage key.
  const std::uint64_t words[] = {tag, value};
  for (std::size_t word = 0; word < 2; ++word) {
    for (unsigned byte = 0; byte < 8; ++byte) {
      *hash ^= static_cast<unsigned char>(words[word] >> (8U * byte));
      *hash *= UINT64_C(1099511628211);
    }
  }
}

double RelativeError(double before, double after) {
  const double scale = std::max(1.0, std::max(std::fabs(before),
                                              std::fabs(after)));
  return std::fabs(after - before) / scale;
}

}  // namespace

Transport::Status ValidateConfiguration(const Configuration& c) {
  if (c.minimumParticlesPerCell < 0 ||
      c.maximumParticlesPerCell < c.minimumParticlesPerCell ||
      c.spatialBins <= 0 || c.momentumBins <= 0 || c.pitchBins <= 0 ||
      !(c.relativeInvariantTolerance >= 0.0) ||
      !std::isfinite(c.relativeInvariantTolerance) || c.lineageSchema == 0)
    return Error("population-control configuration is outside its supported domain");
  return Transport::Status::Ok();
}

Transport::Status ValidateParticle(const ParticleRecord& p) {
  if (p.stableId == 0 || p.species < 0 || !(p.weight > 0.0) ||
      !std::isfinite(p.weight) || !std::isfinite(p.chargeC) ||
      !(p.kineticEnergyJ >= 0.0) || !std::isfinite(p.kineticEnergyJ) ||
      !std::isfinite(p.parallelMomentumKgMPerS) || !std::isfinite(p.mu) ||
      std::fabs(p.mu) > 1.0)
    return Error("population-control particle has invalid identity or SI state");
  return Transport::Status::Ok();
}

Moments ComputeMoments(const std::vector<ParticleRecord>& particles) {
  Moments result;
  for (std::size_t i = 0; i < particles.size(); ++i) {
    const ParticleRecord& p = particles[i];
    result.number += p.weight;
    result.chargeC += p.weight * p.chargeC;
    result.energyJ += p.weight * p.kineticEnergyJ;
    result.parallelMomentumKgMPerS +=
        p.weight * p.parallelMomentumKgMPerS;
    result.muFirst += p.weight * p.mu;
    result.muSecond += p.weight * p.mu * p.mu;
  }
  return result;
}

InvariantReport CompareInvariants(const std::vector<ParticleRecord>& before,
                                  const std::vector<ParticleRecord>& after,
                                  double tolerance) {
  InvariantReport report;
  if (!(tolerance >= 0.0) || !std::isfinite(tolerance)) {
    report.status = Error("population invariant tolerance is invalid");
    return report;
  }
  for (std::size_t i = 0; i < before.size(); ++i) {
    report.status = ValidateParticle(before[i]);
    if (!report.status.ok()) return report;
  }
  for (std::size_t i = 0; i < after.size(); ++i) {
    report.status = ValidateParticle(after[i]);
    if (!report.status.ok()) return report;
  }
  report.before = ComputeMoments(before);
  report.after = ComputeMoments(after);
  const double lhs[] = {report.before.number, report.before.chargeC,
      report.before.energyJ, report.before.parallelMomentumKgMPerS,
      report.before.muFirst, report.before.muSecond};
  const double rhs[] = {report.after.number, report.after.chargeC,
      report.after.energyJ, report.after.parallelMomentumKgMPerS,
      report.after.muFirst, report.after.muSecond};
  const char* names[] = {"number", "charge", "energy", "parallel-momentum",
                         "mu-first", "mu-second"};
  for (std::size_t i = 0; i < 6; ++i) {
    const double error = RelativeError(lhs[i], rhs[i]);
    report.maximumRelativeError = std::max(report.maximumRelativeError, error);
    if (error > tolerance) report.failedInvariants.push_back(names[i]);
  }
  report.status = report.failedInvariants.empty()
      ? Transport::Status::Ok()
      : Transport::Status::Error(Transport::StatusCode::OutOfDomain,
          "population control changed one or more declared invariants");
  return report;
}

std::uint64_t DeriveChildId(std::uint64_t parentId,
                            std::uint64_t generation,
                            std::uint64_t ordinal,
                            std::uint64_t schema) {
  std::uint64_t hash = UINT64_C(14695981039346656037);
  HashWord(UINT64_C(0x534348454d41), schema, &hash);
  HashWord(UINT64_C(0x504152454e54), parentId, &hash);
  HashWord(UINT64_C(0x47454e455241), generation, &hash);
  HashWord(UINT64_C(0x4348494c4421), ordinal, &hash);
  // Zero is reserved as an invalid/uninitialized stable identity.
  return hash == 0 ? UINT64_C(1) : hash;
}

Transport::Status SplitParticle(const ParticleRecord& parent,
                                std::size_t childCount,
                                std::uint64_t schema,
                                std::vector<ParticleRecord>* children) {
  if (!children || childCount < 2 || schema == 0)
    return Error("split requires an output, at least two children, and schema");
  Transport::Status status = ValidateParticle(parent);
  if (!status.ok()) return status;
  children->assign(childCount, parent);
  double assigned = 0.0;
  for (std::size_t i = 0; i < childCount; ++i) {
    ParticleRecord& child = (*children)[i];
    child.parentId = parent.stableId;
    child.lineageGeneration = parent.lineageGeneration + 1;
    child.stableId = DeriveChildId(parent.stableId,
        child.lineageGeneration, static_cast<std::uint64_t>(i), schema);
    child.weight = i + 1 == childCount
        ? parent.weight - assigned
        : parent.weight / static_cast<double>(childCount);
    assigned += child.weight;
    status = ValidateParticle(child);
    if (!status.ok()) return status;
  }
  return Transport::Status::Ok();
}

Transport::Status MergeParticles(const std::vector<ParticleRecord>& parents,
                                  std::uint64_t ordinal,
                                  std::uint64_t schema,
                                  ParticleRecord* merged) {
  if (!merged || parents.size() < 2 || schema == 0)
    return Error("merge requires an output, at least two parents, and schema");
  double weight = 0.0;
  std::uint64_t generation = 0;
  for (std::size_t i = 0; i < parents.size(); ++i) {
    Transport::Status status = ValidateParticle(parents[i]);
    if (!status.ok()) return status;
    if (parents[i].species != parents[0].species ||
        parents[i].chargeC != parents[0].chargeC)
      return Error("merge cannot mix species or charge states");
    weight += parents[i].weight;
    generation = std::max(generation, parents[i].lineageGeneration);
  }
  if (!(weight > 0.0) || !std::isfinite(weight))
    return Error("merged statistical weight is invalid");
  const Moments moments = ComputeMoments(parents);
  *merged = parents[0];
  merged->parentId = parents[0].stableId;
  merged->lineageGeneration = generation + 1;
  merged->stableId = DeriveChildId(parents[0].stableId,
      merged->lineageGeneration, ordinal, schema);
  merged->weight = weight;
  merged->kineticEnergyJ = moments.energyJ / weight;
  merged->parallelMomentumKgMPerS = moments.parallelMomentumKgMPerS / weight;
  merged->mu = moments.muFirst / weight;
  return ValidateParticle(*merged);
}

}  // namespace PopulationControl
}  // namespace SEP
