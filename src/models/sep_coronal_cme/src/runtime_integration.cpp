#include "sep_coronal_cme/runtime_integration.h"

#include "sep_coronal_cme/constants.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <queue>
#include <set>
#include <sstream>

namespace SEP { namespace CoronalCME { namespace {

bool Finite(double value) { return std::isfinite(value); }

bool ValidBlock(const AxisAlignedBlock& block) {
  return block.stableId != 0 && block.maximumM.x > block.minimumM.x &&
      block.maximumM.y > block.minimumM.y &&
      block.maximumM.z > block.minimumM.z;
}

bool SegmentIntersectsBox(Vec3 start, Vec3 end, Vec3 minimum, Vec3 maximum) {
  double lower = 0.0, upper = 1.0;
  const Vec3 direction = end - start;
  const double origins[3] = {start.x, start.y, start.z};
  const double slopes[3] = {direction.x, direction.y, direction.z};
  const double minima[3] = {minimum.x, minimum.y, minimum.z};
  const double maxima[3] = {maximum.x, maximum.y, maximum.z};
  for (int axis = 0; axis < 3; ++axis) {
    if (std::abs(slopes[axis]) < 1.0e-30) {
      if (origins[axis] < minima[axis] || origins[axis] > maxima[axis])
        return false;
      continue;
    }
    double first = (minima[axis] - origins[axis]) / slopes[axis];
    double second = (maxima[axis] - origins[axis]) / slopes[axis];
    if (first > second) std::swap(first, second);
    lower = std::max(lower, first);
    upper = std::min(upper, second);
    if (lower > upper) return false;
  }
  return true;
}

bool FaceAdjacent(const AxisAlignedBlock& a, const AxisAlignedBlock& b,
                  double tolerance) {
  const auto overlaps = [tolerance](double a0, double a1,
                                    double b0, double b1) {
    return std::min(a1, b1) - std::max(a0, b0) >= -tolerance;
  };
  const bool xFace = std::abs(a.maximumM.x - b.minimumM.x) <= tolerance ||
      std::abs(b.maximumM.x - a.minimumM.x) <= tolerance;
  const bool yFace = std::abs(a.maximumM.y - b.minimumM.y) <= tolerance ||
      std::abs(b.maximumM.y - a.minimumM.y) <= tolerance;
  const bool zFace = std::abs(a.maximumM.z - b.minimumM.z) <= tolerance ||
      std::abs(b.maximumM.z - a.minimumM.z) <= tolerance;
  return (xFace && overlaps(a.minimumM.y, a.maximumM.y,
                            b.minimumM.y, b.maximumM.y) &&
                   overlaps(a.minimumM.z, a.maximumM.z,
                            b.minimumM.z, b.maximumM.z)) ||
      (yFace && overlaps(a.minimumM.x, a.maximumM.x,
                         b.minimumM.x, b.maximumM.x) &&
                overlaps(a.minimumM.z, a.maximumM.z,
                         b.minimumM.z, b.maximumM.z)) ||
      (zFace && overlaps(a.minimumM.x, a.maximumM.x,
                         b.minimumM.x, b.maximumM.x) &&
                overlaps(a.minimumM.y, a.maximumM.y,
                         b.minimumM.y, b.maximumM.y));
}

double RelativisticKineticEnergy(double massKg, Vec3 velocity) {
  const double speed2 = Dot(velocity, velocity);
  const double c2 = Constants::kSpeedOfLightMPerS *
      Constants::kSpeedOfLightMPerS;
  if (!(massKg > 0.0) || speed2 < 0.0 || speed2 >= c2) return -1.0;
  const double gamma = 1.0 / std::sqrt(1.0 - speed2 / c2);
  return (gamma - 1.0) * massKg * c2;
}

}  // namespace

Core::Status ValidateSphericalDomain(const SphericalDomain& domain) {
  if (!(domain.solarRadiusM > 0.0 &&
        domain.outerRadiusM > domain.solarRadiusM &&
        domain.cartesianHalfExtentM.x >= domain.outerRadiusM &&
        domain.cartesianHalfExtentM.y >= domain.outerRadiusM &&
        domain.cartesianHalfExtentM.z >= domain.outerRadiusM)) {
    return Core::Status::Failure(Core::StatusCode::InvalidConfiguration,
        "Cartesian hierarchy must enclose the exact inner/outer spheres");
  }
  return Core::Status::Success();
}

RadialDomainLocation ClassifyRadialDomain(const SphericalDomain& domain,
                                          Vec3 positionM) {
  const double radius = Norm(positionM);
  if (radius < domain.solarRadiusM) return RadialDomainLocation::SolarInterior;
  if (radius > domain.outerRadiusM) return RadialDomainLocation::Escaped;
  return RadialDomainLocation::PhysicalDomain;
}

Core::Result<double> RequestedCellSize(
    const RefinementControls& controls, double solarDistanceM,
    double tubeDistanceM) {
  if (!(controls.globalCellSizeM > 0.0 && controls.solarCellSizeM > 0.0 &&
        controls.tubeCellSizeM > 0.0 && controls.solarDecayLengthM > 0.0 &&
        controls.tubeDecayLengthM > 0.0 && solarDistanceM >= 0.0 &&
        tubeDistanceM >= 0.0 &&
        controls.solarCellSizeM <= controls.globalCellSizeM &&
        controls.tubeCellSizeM <= controls.globalCellSizeM)) {
    return Core::Result<double>::Failure(
        Core::StatusCode::InvalidConfiguration,
        "refinement sizes/decays/distances must be finite, positive, and ordered");
  }
  const auto request = [global = controls.globalCellSizeM](
      double target, double distance, double decay) {
    return target + (global - target) * (1.0 - std::exp(-distance / decay));
  };
  return Core::Result<double>::Success(std::min(
      request(controls.solarCellSizeM, solarDistanceM,
              controls.solarDecayLengthM),
      request(controls.tubeCellSizeM, tubeDistanceM,
              controls.tubeDecayLengthM)));
}

bool BlockIntersectsBufferedPolyline(const AxisAlignedBlock& block,
                                     const std::vector<Vec3>& centerlineM,
                                     double bufferRadiusM) {
  if (!ValidBlock(block) || centerlineM.empty() || bufferRadiusM < 0.0)
    return false;
  const Vec3 minimum{block.minimumM.x - bufferRadiusM,
                     block.minimumM.y - bufferRadiusM,
                     block.minimumM.z - bufferRadiusM};
  const Vec3 maximum{block.maximumM.x + bufferRadiusM,
                     block.maximumM.y + bufferRadiusM,
                     block.maximumM.z + bufferRadiusM};
  if (centerlineM.size() == 1)
    return SegmentIntersectsBox(centerlineM[0], centerlineM[0], minimum, maximum);
  for (std::size_t i = 1; i < centerlineM.size(); ++i)
    if (SegmentIntersectsBox(centerlineM[i - 1], centerlineM[i], minimum, maximum))
      return true;
  return false;
}

Core::Status ValidateFaceConnectedBlocks(
    const std::vector<AxisAlignedBlock>& blocks,
    const std::vector<std::uint64_t>& requiredStableIds, double tolerance) {
  if (blocks.empty() || requiredStableIds.empty() || tolerance < 0.0)
    return Core::Status::Failure(Core::StatusCode::InvalidConfiguration,
                                "corridor connectivity requires blocks/anchors");
  std::vector<bool> reached(blocks.size(), false);
  std::queue<std::size_t> pending;
  std::set<std::uint64_t> required(requiredStableIds.begin(),
                                   requiredStableIds.end());
  for (std::size_t i = 0; i < blocks.size(); ++i) {
    if (!ValidBlock(blocks[i]))
      return Core::Status::Failure(Core::StatusCode::InvalidState,
                                  "corridor contains an invalid block");
    if (blocks[i].active && blocks[i].stableId == requiredStableIds.front()) {
      reached[i] = true;
      pending.push(i);
    }
  }
  if (pending.empty())
    return Core::Status::Failure(Core::StatusCode::InvalidState,
                                "corridor start anchor is inactive or absent");
  while (!pending.empty()) {
    const std::size_t current = pending.front();
    pending.pop();
    required.erase(blocks[current].stableId);
    for (std::size_t next = 0; next < blocks.size(); ++next) {
      if (!reached[next] && blocks[next].active &&
          FaceAdjacent(blocks[current], blocks[next], tolerance)) {
        reached[next] = true;
        pending.push(next);
      }
    }
  }
  if (!required.empty())
    return Core::Status::Failure(Core::StatusCode::InvalidState,
                                "active corridor is not face connected");
  return Core::Status::Success();
}

Core::Status ValidateOneSidedInterpolationStencil(
    const std::vector<CategoricalFieldIdentity>& identities) {
  if (identities.empty())
    return Core::Status::Failure(Core::StatusCode::InvalidConfiguration,
                                "interpolation stencil is empty");
  const auto& first = identities.front();
  if (!first.mappingValid)
    return Core::Status::Failure(Core::StatusCode::InvalidState,
                                "interpolation mapping is invalid");
  for (const auto& identity : identities) {
    if (!identity.mappingValid || identity.region != first.region ||
        identity.sector != first.sector ||
        identity.interfaceIdentity != first.interfaceIdentity)
      return Core::Status::Failure(Core::StatusCode::InvalidState,
          "interpolation stencil crosses a categorical interface");
  }
  return Core::Status::Success();
}

Core::Result<IdealHcsCrossingState> CrossIdealHcs(
    const IdealHcsCrossingState& before, bool driftRequested) {
  if (driftRequested)
    return Core::Result<IdealHcsCrossingState>::Failure(
        Core::StatusCode::NotImplemented,
        "HCS drift requires a qualified finite-thickness sheet");
  if (before.magneticSector != -1 && before.magneticSector != 1)
    return Core::Result<IdealHcsCrossingState>::Failure(
        Core::StatusCode::InvalidState, "magnetic sector must be +/-1");
  IdealHcsCrossingState after = before;
  after.magneticSector *= -1;
  // Labels are tied to physical outward/inward propagation, not the local
  // sign of B, so they and the physical phase-space state remain unchanged.
  return Core::Result<IdealHcsCrossingState>::Success(after);
}

Core::Result<double> ComputeLocalTimeStep(const TimeStepControls& controls,
                                          const LocalStepInputs& inputs) {
  if (!(controls.fixedUpperBoundS > 0.0 && controls.gyroAccuracy > 0.0 &&
        controls.spatialAccuracy > 0.0 &&
        controls.scatteringAccuracy > 0.0 &&
        inputs.gyrofrequencyRadPerS > 0.0 && inputs.speedMPerS > 0.0 &&
        inputs.cellSizeM > 0.0 && inputs.scatteringTimeS > 0.0)) {
    return Core::Result<double>::Failure(Core::StatusCode::InvalidConfiguration,
        "local time step requires positive bounds, accuracies, and scales");
  }
  const double gyro = controls.gyroAccuracy /
      inputs.gyrofrequencyRadPerS;
  const double crossing = controls.spatialAccuracy * inputs.cellSizeM /
      inputs.speedMPerS;
  const double scattering = controls.scatteringAccuracy *
      inputs.scatteringTimeS;
  return Core::Result<double>::Success(
      std::min({controls.fixedUpperBoundS, gyro, crossing, scattering}));
}

Core::Status ValidateAllSpeciesWeights(
    const std::vector<SpeciesWeight>& weights, int count) {
  if (count <= 0 || weights.size() != static_cast<std::size_t>(count))
    return Core::Status::Failure(Core::StatusCode::InvalidConfiguration,
                                "one weight is required per compiled species");
  std::set<int> slots;
  std::set<std::string> ids;
  for (const auto& weight : weights) {
    if (weight.compiledSlot < 0 || weight.compiledSlot >= count ||
        weight.stableSpeciesId.empty() || !(weight.baseWeight > 0.0) ||
        !slots.insert(weight.compiledSlot).second ||
        !ids.insert(weight.stableSpeciesId).second)
      return Core::Status::Failure(Core::StatusCode::InvalidConfiguration,
                                  "species weights are incomplete or ambiguous");
  }
  return Core::Status::Success();
}

PopulationMoments AccumulatePopulationMoments(
    const std::vector<PopulationParticle>& particles) {
  PopulationMoments total;
  for (const auto& particle : particles) {
    total.representedWeight += particle.representedWeight;
    total.representedMomentumKgMPerS = total.representedMomentumKgMPerS +
        particle.representedWeight * particle.momentumKgMPerS;
    total.representedKineticEnergyJ +=
        particle.representedWeight * particle.kineticEnergyJ;
  }
  return total;
}

Core::Result<std::vector<PopulationParticle>> SplitParticleConservatively(
    const PopulationParticle& particle, int copies) {
  if (copies < 2 || !(particle.representedWeight > 0.0))
    return Core::Result<std::vector<PopulationParticle>>::Failure(
        Core::StatusCode::InvalidConfiguration,
        "particle split requires positive weight and at least two copies");
  std::vector<PopulationParticle> result(static_cast<std::size_t>(copies),
                                         particle);
  for (int i = 0; i < copies; ++i) {
    result[static_cast<std::size_t>(i)].id = particle.id * 1024U +
        static_cast<std::uint64_t>(i + 1);
    result[static_cast<std::size_t>(i)].representedWeight =
        particle.representedWeight / copies;
  }
  return Core::Result<std::vector<PopulationParticle>>::Success(result);
}

Core::Result<PopulationParticle> MergeIdenticalStateParticles(
    const std::vector<PopulationParticle>& particles, std::uint64_t mergedId) {
  if (particles.size() < 2 || mergedId == 0)
    return Core::Result<PopulationParticle>::Failure(
        Core::StatusCode::InvalidConfiguration,
        "merge requires at least two particles and a stable output ID");
  PopulationParticle merged = particles.front();
  merged.id = mergedId;
  merged.representedWeight = 0.0;
  for (const auto& particle : particles) {
    if (particle.speciesSlot != merged.speciesSlot ||
        particle.backgroundGeneration != merged.backgroundGeneration ||
        particle.sourceLabel != merged.sourceLabel ||
        Norm(particle.momentumKgMPerS - merged.momentumKgMPerS) > 0.0 ||
        particle.kineticEnergyJ != merged.kineticEnergyJ ||
        !(particle.representedWeight > 0.0))
      return Core::Result<PopulationParticle>::Failure(
          Core::StatusCode::InvalidState,
          "merge group crosses identity or phase-space state");
    merged.representedWeight += particle.representedWeight;
  }
  return Core::Result<PopulationParticle>::Success(merged);
}

void MarkInitialized(InitializationLedger* ledger,
                     InitializationCondition condition) {
  if (ledger) ledger->completedMask |= static_cast<std::uint32_t>(condition);
}

Core::Status ValidateInitializationForOutput(
    const InitializationLedger& ledger) {
  constexpr std::uint32_t all = (1U << 10) - 1U;
  if (ledger.completedMask != all) {
    // Keep the ten-stage gate strict, but report its exact missing stages.
    // An aggregate failure previously hid the distinction between an
    // unprepared provider and a physically inactive, prepared shock.
    struct StageName {
      InitializationCondition condition;
      const char* name;
    };
    constexpr StageName stages[] = {
        {InitializationCondition::Mesh, "Mesh"},
        {InitializationCondition::Boundary, "Boundary"},
        {InitializationCondition::Background, "Background"},
        {InitializationCondition::Turbulence, "Turbulence"},
        {InitializationCondition::Shock, "Shock"},
        {InitializationCondition::Halo, "Halo"},
        {InitializationCondition::Species, "Species"},
        {InitializationCondition::TimeStep, "TimeStep"},
        {InitializationCondition::Observers, "Observers"},
        {InitializationCondition::OutputDictionary, "OutputDictionary"}};
    std::ostringstream message;
    message << "initialization output requested before all ten stages: missing=";
    bool first = true;
    for (const StageName& stage : stages) {
      if ((ledger.completedMask &
           static_cast<std::uint32_t>(stage.condition)) != 0)
        continue;
      if (!first) message << ',';
      message << stage.name;
      first = false;
    }
    if (first) message << "none";
    message << "; completed_mask=0x" << std::hex << ledger.completedMask
            << "; expected_mask=0x" << all;
    const std::uint32_t unexpected = ledger.completedMask & ~all;
    if (unexpected != 0)
      message << "; unexpected_mask_bits=0x" << unexpected;
    return Core::Status::Failure(Core::StatusCode::InvalidState,
                                message.str());
  }
  if (ledger.backgroundGeneration == 0 ||
      ledger.shockBackgroundGeneration != ledger.backgroundGeneration)
    return Core::Status::Failure(Core::StatusCode::InvalidState,
                                "shock/background generations are stale");
  return Core::Status::Success();
}

FiniteOutputValue MakeFiniteOutput(double value, ValueValidity validity,
                                   bool present, double sentinel) {
  if (validity == ValueValidity::Valid && Finite(value))
    return {value, validity, present};
  return {Finite(sentinel) ? sentinel : 0.0, validity, present};
}

Core::Result<std::vector<double>> BuildEnergyEdges(
    EnergyGrid grid, double minimumJ, double maximumJ, int channels,
    const std::vector<double>& explicitEdges) {
  if (grid == EnergyGrid::Explicit) {
    if (explicitEdges.size() < 2)
      return Core::Result<std::vector<double>>::Failure(
          Core::StatusCode::InvalidConfiguration,
          "explicit energy grid requires at least two edges");
    for (std::size_t i = 0; i < explicitEdges.size(); ++i)
      if (!(explicitEdges[i] > 0.0) ||
          (i && explicitEdges[i] <= explicitEdges[i - 1]))
        return Core::Result<std::vector<double>>::Failure(
            Core::StatusCode::InvalidConfiguration,
            "explicit energy edges must be positive and increasing");
    return Core::Result<std::vector<double>>::Success(explicitEdges);
  }
  if (!(minimumJ > 0.0 && maximumJ > minimumJ && channels > 0))
    return Core::Result<std::vector<double>>::Failure(
        Core::StatusCode::InvalidConfiguration,
        "generated energy grid requires 0<minimum<maximum and channels>0");
  std::vector<double> edges(static_cast<std::size_t>(channels + 1));
  for (int i = 0; i <= channels; ++i) {
    const double f = static_cast<double>(i) / channels;
    edges[static_cast<std::size_t>(i)] = grid == EnergyGrid::Linear
        ? minimumJ + f * (maximumJ - minimumJ)
        : minimumJ * std::pow(maximumJ / minimumJ, f);
  }
  return Core::Result<std::vector<double>>::Success(edges);
}

Core::Result<std::vector<ObserverBin>> AccumulateObserver(
    const ObserverDefinition& observer,
    const std::vector<ObserverParticle>& particles) {
  if (observer.stableId.empty() || !(observer.radiusM > 0.0) ||
      observer.energyEdgesJ.size() < 2)
    return Core::Result<std::vector<ObserverBin>>::Failure(
        Core::StatusCode::InvalidConfiguration,
        "observer requires stable ID, geometry, and energy edges");
  if (observer.geometry == ObserverGeometry::SurfaceDisk &&
      Norm(observer.surfaceNormal) == 0.0)
    return Core::Result<std::vector<ObserverBin>>::Failure(
        Core::StatusCode::InvalidConfiguration,
        "disk observer requires a surface normal");
  std::vector<ObserverBin> bins(observer.energyEdgesJ.size() - 1);
  const double volume = 4.0 * Constants::kPi *
      observer.radiusM * observer.radiusM * observer.radiusM / 3.0;
  const double area = Constants::kPi * observer.radiusM * observer.radiusM;
  for (const auto& particle : particles) {
    const Vec3 relativeVelocity = particle.velocityMPerS - observer.velocityMPerS;
    double energy = RelativisticKineticEnergy(particle.massKg,
                                              relativeVelocity);
    if (observer.energyPerNucleon) {
      if (particle.nucleonCount <= 0)
        return Core::Result<std::vector<ObserverBin>>::Failure(
            Core::StatusCode::InvalidState,
            "per-nucleon observer requires an integer nucleon count");
      energy /= particle.nucleonCount;
    }
    if (!(energy >= observer.energyEdgesJ.front() &&
          energy < observer.energyEdgesJ.back()) ||
        !(particle.representedWeight >= 0.0)) continue;
    const auto upper = std::upper_bound(observer.energyEdgesJ.begin(),
                                        observer.energyEdgesJ.end(), energy);
    const std::size_t index = static_cast<std::size_t>(
        std::distance(observer.energyEdgesJ.begin(), upper) - 1);
    bool accepted = false;
    double estimator = 0.0;
    if (observer.geometry == ObserverGeometry::VolumeSphere) {
      accepted = Norm(particle.positionM - observer.centerM) <= observer.radiusM;
      estimator = particle.representedWeight * particle.residenceTimeS / volume;
    } else {
      accepted = particle.crossedSurface && particle.crossingSense > 0;
      const double geometryArea = observer.geometry == ObserverGeometry::SurfaceDisk
          ? area : 4.0 * area;
      estimator = particle.representedWeight / geometryArea;
    }
    if (accepted) {
      bins[index].weightedCount += particle.representedWeight;
      bins[index].intensity += estimator;
      bins[index].particlesPresent = true;
    }
  }
  return Core::Result<std::vector<ObserverBin>>::Success(bins);
}

double DeterministicPhysicalSum(
    const std::vector<std::pair<std::uint64_t, double>>& contributions) {
  auto sorted = contributions;
  std::sort(sorted.begin(), sorted.end(),
            [](const auto& a, const auto& b) { return a.first < b.first; });
  double sum = 0.0, correction = 0.0;
  for (const auto& contribution : sorted) {
    const double adjusted = contribution.second - correction;
    const double next = sum + adjusted;
    correction = (next - sum) - adjusted;
    sum = next;
  }
  return sum;
}

} }  // namespace SEP::CoronalCME
