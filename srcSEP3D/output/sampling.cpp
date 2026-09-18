#include "sampling.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <utility>

namespace SEP3D {
namespace Output {
namespace {

Core::Status Invalid(const char* message) {
  return Core::Status(Core::StatusCode::InvalidInput, message);
}

bool Finite(const Core::Vec3& value) {
  return std::isfinite(value.x) && std::isfinite(value.y) &&
         std::isfinite(value.z);
}

// Neumaier summation gives a deterministic, compensated reduction once the
// observations are in canonical stable-ID order. It materially reduces the
// loss caused by broad SEP statistical-weight distributions.
class Sum {
 public:
  void Add(double value) {
    const double next = total_ + value;
    if (std::fabs(total_) >= std::fabs(value))
      correction_ += (total_ - next) + value;
    else
      correction_ += (value - next) + total_;
    total_ = next;
  }
  double Value() const { return total_ + correction_; }
 private:
  double total_ = 0.0;
  double correction_ = 0.0;
};

double Speed(double momentum, double mass) {
  const double mc = mass * Core::Const::c;
  return Core::Const::c * momentum / std::sqrt(momentum * momentum + mc * mc);
}

double KineticEnergy(double momentum, double mass) {
  const double mc = mass * Core::Const::c;
  return (std::sqrt(momentum * momentum + mc * mc) - mc) * Core::Const::c;
}

bool ValidEdges(const std::vector<double>& edges) {
  if (edges.size() < 2) return false;
  for (std::size_t i = 0; i < edges.size(); ++i) {
    if (!std::isfinite(edges[i]) ||
        (i != 0 && !(edges[i] > edges[i - 1]))) return false;
  }
  return true;
}

bool AcceptsSpecies(const VirtualSpacecraftDefinition& craft, int species) {
  return craft.acceptedSpecies.empty() ||
      std::find(craft.acceptedSpecies.begin(), craft.acceptedSpecies.end(),
                species) != craft.acceptedSpecies.end();
}

std::size_t Bin(const std::vector<double>& edges, double value) {
  if (value < edges.front() || value > edges.back()) return edges.size();
  if (value == edges.back()) return edges.size() - 2;
  return static_cast<std::size_t>(
      std::upper_bound(edges.begin(), edges.end(), value) - edges.begin() - 1);
}

}  // namespace

SamplingSnapshot Sample(const SamplingRequest& request) {
  SamplingSnapshot result;
  std::map<std::uint64_t, CellDefinition> cells;
  for (const CellDefinition& cell : request.cells) {
    if (cell.cellId == 0 || !Finite(cell.centerM) ||
        !std::isfinite(cell.volumeM3) || cell.volumeM3 <= 0.0 ||
        !cells.emplace(cell.cellId, cell).second) {
      result.status = Invalid("cell sampling definitions are invalid or duplicate");
      return result;
    }
  }
  for (const VirtualSpacecraftDefinition& craft : request.spacecraft) {
    if (craft.name.empty() || !Finite(craft.positionM) ||
        !std::isfinite(craft.collectionRadiusM) ||
        craft.collectionRadiusM <= 0.0 ||
        !ValidEdges(craft.kineticEnergyEdgesJ) ||
        !std::isfinite(craft.minimumMu) || !std::isfinite(craft.maximumMu) ||
        craft.minimumMu < -1.0 || craft.maximumMu > 1.0 ||
        craft.maximumMu <= craft.minimumMu || craft.observerKind.empty() ||
        craft.normalization.empty()) {
      result.status = Invalid("virtual-spacecraft definition is invalid");
      return result;
    }
  }
  for (const FieldLineProjectionDefinition& line : request.fieldLines) {
    if (line.name.empty() || !Finite(line.originM) || !Finite(line.direction) ||
        std::fabs(line.direction.Norm() - 1.0) > 1.0e-12 ||
        !ValidEdges(line.distanceEdgesM)) {
      result.status = Invalid("field-line projection definition is invalid");
      return result;
    }
  }

  std::vector<ParticleObservation> particles = request.particles;
  std::sort(particles.begin(), particles.end(),
            [](const ParticleObservation& left,
               const ParticleObservation& right) {
              if (left.stableId != right.stableId)
                return left.stableId < right.stableId;
              if (left.species != right.species)
                return left.species < right.species;
              return left.cellId < right.cellId;
            });
  for (std::size_t i = 0; i < particles.size(); ++i) {
    const ParticleObservation& particle = particles[i];
    if (particle.stableId == 0 || particle.species < 0 ||
        cells.find(particle.cellId) == cells.end() ||
        !Finite(particle.positionM) ||
        !std::isfinite(particle.momentumKgMPerS) ||
        particle.momentumKgMPerS < 0.0 ||
        !std::isfinite(particle.restMassKg) || particle.restMassKg <= 0.0 ||
        !std::isfinite(particle.mu) || particle.mu < -1.0 || particle.mu > 1.0 ||
        !std::isfinite(particle.statisticalWeight) ||
        particle.statisticalWeight <= 0.0 ||
        (i != 0 && particle.stableId == particles[i - 1].stableId)) {
      result.status = Invalid("particle observation is invalid or has a duplicate stable ID");
      return result;
    }
  }

  // Accumulators are keyed by explicit physical identity; map ordering makes
  // the published row order independent of AMR block and MPI traversal order.
  struct MomentAcc { Sum w, fx, fy, fz, energy, wmu; };
  std::map<std::pair<std::uint64_t, int>, MomentAcc> moments;
  for (const ParticleObservation& particle : particles) {
    MomentAcc& acc = moments[{particle.cellId, particle.species}];
    const double speed = Speed(particle.momentumKgMPerS, particle.restMassKg);
    const double energy = KineticEnergy(
        particle.momentumKgMPerS, particle.restMassKg);
    acc.w.Add(particle.statisticalWeight);
    // A gyrotropic record carries only the field-aligned first moment. Store
    // that vector along the radial direction as a coordinate-free diagnostic
    // when no local b-hat is part of the sampling input.
    const Core::Vec3 radial = particle.positionM.Normalized();
    const Core::Vec3 flux = radial *
        (particle.statisticalWeight * particle.mu * speed);
    acc.fx.Add(flux.x); acc.fy.Add(flux.y); acc.fz.Add(flux.z);
    acc.energy.Add(particle.statisticalWeight * energy);
    acc.wmu.Add(particle.statisticalWeight * particle.mu);
  }
  for (const auto& item : moments) {
    const CellDefinition& cell = cells[item.first.first];
    const MomentAcc& acc = item.second;
    CellMoment moment;
    moment.cellId = item.first.first;
    moment.species = item.first.second;
    moment.representedParticles = acc.w.Value();
    moment.numberDensityM3 = acc.w.Value() / cell.volumeM3;
    moment.weightedFluxM2PerS = Core::Vec3(
        acc.fx.Value(), acc.fy.Value(), acc.fz.Value()) / cell.volumeM3;
    moment.kineticEnergyDensityJPerM3 = acc.energy.Value() / cell.volumeM3;
    moment.firstPitchMoment = acc.wmu.Value() / acc.w.Value();
    result.cellMoments.push_back(moment);
  }

  for (const VirtualSpacecraftDefinition& craft : request.spacecraft) {
    struct CraftAcc {
      std::vector<Sum> weight;
      std::vector<Sum> weightSquared;
      Sum totalWeight;
      Sum weightedMu;
      std::uint64_t macroCount = 0;
    };
    std::map<int, CraftAcc> bySpecies;
    for (const ParticleObservation& particle : particles) {
      if ((particle.positionM - craft.positionM).Norm() >
          craft.collectionRadiusM || !AcceptsSpecies(craft, particle.species) ||
          particle.mu < craft.minimumMu || particle.mu > craft.maximumMu)
        continue;
      const double energy = KineticEnergy(
          particle.momentumKgMPerS, particle.restMassKg);
      const std::size_t bin = Bin(craft.kineticEnergyEdgesJ, energy);
      if (bin >= craft.kineticEnergyEdgesJ.size() - 1) continue;
      auto& acc = bySpecies[particle.species];
      if (acc.weight.empty()) {
        acc.weight.resize(craft.kineticEnergyEdgesJ.size() - 1);
        acc.weightSquared.resize(craft.kineticEnergyEdgesJ.size() - 1);
      }
      acc.weight[bin].Add(particle.statisticalWeight);
      acc.weightSquared[bin].Add(
          particle.statisticalWeight * particle.statisticalWeight);
      acc.totalWeight.Add(particle.statisticalWeight);
      acc.weightedMu.Add(particle.statisticalWeight * particle.mu);
      ++acc.macroCount;
    }
    for (const auto& item : bySpecies) {
      VirtualSpacecraftProduct product;
      product.name = craft.name; product.species = item.first;
      product.kineticEnergyEdgesJ = craft.kineticEnergyEdgesJ;
      product.observerKind = craft.observerKind;
      product.normalization = craft.normalization;
      product.acceptedMacroparticles = item.second.macroCount;
      for (std::size_t i = 0; i < item.second.weight.size(); ++i) {
        const double width = craft.kineticEnergyEdgesJ[i + 1] -
            craft.kineticEnergyEdgesJ[i];
        product.representedParticlesPerJ.push_back(
            item.second.weight[i].Value() / width);
        product.standardUncertaintyPerJ.push_back(
            std::sqrt(std::max(0.0,
                item.second.weightSquared[i].Value())) / width);
      }
      product.dipoleAnisotropy = item.second.totalWeight.Value() == 0.0 ? 0.0
          : 3.0 * item.second.weightedMu.Value() /
                item.second.totalWeight.Value();
      result.spacecraft.push_back(product);
    }
  }

  for (const FieldLineProjectionDefinition& line : request.fieldLines) {
    std::map<int, std::vector<Sum>> bySpecies;
    for (const ParticleObservation& particle : particles) {
      const double distance = (particle.positionM - line.originM).Dot(line.direction);
      const std::size_t bin = Bin(line.distanceEdgesM, distance);
      if (bin >= line.distanceEdgesM.size() - 1) continue;
      auto& acc = bySpecies[particle.species];
      if (acc.empty()) acc.resize(line.distanceEdgesM.size() - 1);
      acc[bin].Add(particle.statisticalWeight);
    }
    for (const auto& item : bySpecies) {
      FieldLineProjection projection;
      projection.name = line.name; projection.species = item.first;
      projection.distanceEdgesM = line.distanceEdgesM;
      for (std::size_t i = 0; i < item.second.size(); ++i)
        projection.representedParticlesPerM.push_back(
            item.second[i].Value() /
            (line.distanceEdgesM[i + 1] - line.distanceEdgesM[i]));
      result.fieldLines.push_back(projection);
    }
  }

  for (const Adapters::LedgerRow& row : request.ledgerRows) {
    if (!row.closed) {
      result.status = Invalid("shock diagnostics require closed ledger rows");
      return result;
    }
    ShockDiagnostic diagnostic;
    diagnostic.step = row.key.step; diagnostic.species = row.key.species;
    diagnostic.injected = row.injected; diagnostic.escaped = row.escaped;
    diagnostic.absorbed = row.absorbed; diagnostic.failed = row.failed;
    diagnostic.crossings = row.shockCrossings;
    result.shocks.push_back(diagnostic);
  }
  std::sort(result.shocks.begin(), result.shocks.end(),
            [](const ShockDiagnostic& left, const ShockDiagnostic& right) {
              return left.step < right.step ||
                  (left.step == right.step && left.species < right.species);
            });
  result.nextState = request.previousState;
  ++result.nextState.completedSamplings;
  result.nextState.observationsProcessed += particles.size();
  result.nextState.pendingWindows = 0;
  result.nextState.pendingObservations = 0;
  result.nextState.pendingRepresentedParticles = 0.0;
  result.status = Core::Status::OK();
  return result;
}

}  // namespace Output
}  // namespace SEP3D
