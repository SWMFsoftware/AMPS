#include "test_framework.h"

#include "sep_coronal_cme/constants.h"
#include "sep_coronal_cme/runtime_integration.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

namespace SCCMTest { namespace {
using namespace SEP::CoronalCME;

bool Close(double a, double b, double tolerance = 1.0e-10) {
  return std::abs(a - b) <= tolerance *
      std::max({1.0, std::abs(a), std::abs(b)});
}

SphericalDomain Domain() { return {1.0, 10.0, {10.0, 10.0, 10.0}}; }

AxisAlignedBlock Block(std::uint64_t id, double x0, double x1,
                       bool active = true) {
  return {id, {x0, -0.5, -0.5}, {x1, 0.5, 0.5}, active};
}

ObserverDefinition VolumeObserver() {
  ObserverDefinition observer;
  observer.stableId = "earth";
  observer.geometry = ObserverGeometry::VolumeSphere;
  observer.radiusM = 2.0;
  observer.energyEdgesJ = {1.0e-23, 1.0e-20, 1.0e-15};
  return observer;
}

ObserverParticle Particle(Vec3 velocity = {1000.0, 0.0, 0.0}) {
  ObserverParticle particle;
  particle.positionM = {0.0, 0.0, 0.0};
  particle.velocityMPerS = velocity;
  particle.massKg = 1.67262192369e-27;
  particle.representedWeight = 4.0;
  particle.residenceTimeS = 2.0;
  particle.crossedSurface = true;
  particle.crossingSense = 1;
  particle.nucleonCount = 1;
  return particle;
}

void BND3D01() {
  const auto domain = Domain();
  Require(ValidateSphericalDomain(domain).ok() &&
          ClassifyRadialDomain(domain, {0.5, 0.0, 0.0}) ==
              RadialDomainLocation::SolarInterior &&
          ClassifyRadialDomain(domain, {1.0, 0.0, 0.0}) ==
              RadialDomainLocation::PhysicalDomain,
          "solar sphere is not the unique inner absorbing boundary");
}

void BND3D02() {
  auto domain = Domain();
  Require(ValidateSphericalDomain(domain).ok() &&
          ClassifyRadialDomain(domain, {10.0, 10.0, 10.0}) ==
              RadialDomainLocation::Escaped,
          "Cartesian corner was treated as part of the spherical domain");
  domain.cartesianHalfExtentM.x = 9.0;
  Require(!ValidateSphericalDomain(domain).ok(),
          "Cartesian hierarchy accepted an unclipped outer sphere");
}

void MESH3D01() {
  RefinementControls controls{10.0, 1.0, 2.0, 2.0, 8.0};
  auto solar = RequestedCellSize(controls, 0.0, 100.0);
  auto tube = RequestedCellSize(controls, 100.0, 0.0);
  auto global = RequestedCellSize(controls, 100.0, 100.0);
  Require(solar.ok() && tube.ok() && global.ok() &&
          Close(solar.value, 1.0) && Close(tube.value, 2.0) &&
          global.value > tube.value,
          "independent solar/tube refinement fields were not combined by min");
}

void COR3D01() {
  AxisAlignedBlock block{1, {0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}, true};
  Require(BlockIntersectsBufferedPolyline(
              block, {{-2.0, 0.9, 0.9}, {2.0, 0.9, 0.9}}, 0.0),
          "segment/block intersection relied on the block center");
}

void COR3D02() {
  std::vector<AxisAlignedBlock> connected{
      Block(1, 0.0, 1.0), Block(2, 1.0, 2.0), Block(3, 2.0, 3.0)};
  Require(ValidateFaceConnectedBlocks(connected, {1, 3}).ok(),
          "face-connected corridor was rejected");
  connected[1].active = false;
  Require(!ValidateFaceConnectedBlocks(connected, {1, 3}).ok(),
          "corridor hole did not fail connectivity closure");
}

void COR3D03() {
  const AxisAlignedBlock block = Block(1, 0.0, 1.0);
  const std::vector<Vec3> line{{-1.0, 1.2, 0.0}, {2.0, 1.2, 0.0}};
  Require(!BlockIntersectsBufferedPolyline(block, line, 0.69) &&
          BlockIntersectsBufferedPolyline(block, line, 0.71),
          "corridor width refinement is not monotone/conservative");
}

void TIM3D01() {
  LocalStepInputs input{2.0, 10.0, 20.0, 5.0};
  auto coarse = ComputeLocalTimeStep({10.0, 0.4, 0.5, 0.5}, input);
  auto fine = ComputeLocalTimeStep({10.0, 0.2, 0.25, 0.25}, input);
  Require(coarse.ok() && fine.ok() && fine.value <= coarse.value &&
          Close(coarse.value, 0.2) && Close(fine.value, 0.1),
          "time-step accuracy refinement is not monotone");
}

void TIM3D02() {
  auto fixed = ComputeLocalTimeStep(
      {0.01, 1.0, 1.0, 1.0}, {1.0, 1.0, 1.0, 1.0});
  Require(fixed.ok() && Close(fixed.value, 0.01),
          "fixed run step cap was not applied exactly");
}

PopulationParticle Population() {
  return {7, 1, 9, "shock-A", 12.0, {2.0, -1.0, 0.5}, 4.0};
}

void POP3D01() {
  auto split = SplitParticleConservatively(Population(), 4);
  Require(split.ok() && split.value.size() == 4 &&
          Close(AccumulatePopulationMoments(split.value).representedWeight, 12.0),
          "conservative split did not reach the requested population");
}

void POP3D02() {
  const auto original = AccumulatePopulationMoments({Population()});
  auto split = SplitParticleConservatively(Population(), 3);
  auto merged = MergeIdenticalStateParticles(split.value, 99);
  const auto final = AccumulatePopulationMoments({merged.value});
  Require(merged.ok() && Close(original.representedWeight, final.representedWeight) &&
          Norm(original.representedMomentumKgMPerS -
               final.representedMomentumKgMPerS) < 1.0e-12 &&
          Close(original.representedKineticEnergyJ,
                final.representedKineticEnergyJ),
          "split/merge changed weight, momentum, or kinetic energy");
}

void POP3D03() {
  auto particle = Population();
  auto incompatible = particle;
  incompatible.sourceLabel = "shock-B";
  Require(!MergeIdenticalStateParticles({particle, incompatible}, 20).ok(),
          "population control crossed a source-label group");
}

void MPI3D01() {
  const std::vector<std::pair<std::uint64_t, double>> first{
      {3, 1.0e16}, {1, 1.0}, {2, -1.0e16}};
  const std::vector<std::pair<std::uint64_t, double>> second{
      {2, -1.0e16}, {3, 1.0e16}, {1, 1.0}};
  Require(DeterministicPhysicalSum(first) == DeterministicPhysicalSum(second),
          "stable-ID reduction depends on rank/request order");
}

InitializationLedger CompleteLedger() {
  InitializationLedger ledger;
  for (std::uint32_t bit = 0; bit < 10; ++bit)
    ledger.completedMask |= (1U << bit);
  ledger.backgroundGeneration = 4;
  ledger.shockBackgroundGeneration = 4;
  return ledger;
}

void INIT3D01() {
  auto ledger = CompleteLedger();
  Require(ValidateInitializationForOutput(ledger).ok(),
          "complete ordered initialization was rejected");

  // Exercise every incomplete stage independently, including the Shock bit
  // missing from the reported pre-activation native run (0x3ef).  Diagnostics
  // must identify the actual absent stage without relaxing acceptance.
  const char* names[] = {"Mesh", "Boundary", "Background", "Turbulence",
                        "Shock", "Halo", "Species", "TimeStep", "Observers",
                        "OutputDictionary"};
  for (std::uint32_t bit = 0; bit < 10; ++bit) {
    auto incomplete = CompleteLedger();
    incomplete.completedMask &= ~(1U << bit);
    const auto result = ValidateInitializationForOutput(incomplete);
    Require(!result.ok() &&
                result.message.find(std::string("missing=") + names[bit] +
                                    ";") != std::string::npos,
            "missing initialization stage was accepted or not identified");
  }
  auto multiple = CompleteLedger();
  multiple.completedMask &= ~static_cast<std::uint32_t>(
      InitializationCondition::Halo);
  multiple.completedMask &= ~static_cast<std::uint32_t>(
      InitializationCondition::Observers);
  const auto multipleResult = ValidateInitializationForOutput(multiple);
  Require(!multipleResult.ok() &&
              multipleResult.message.find("missing=Halo,Observers;") !=
                  std::string::npos,
          "multiple missing stages were not reported in initialization order");
  auto unknown = CompleteLedger();
  unknown.completedMask |= (1U << 10);
  const auto unknownResult = ValidateInitializationForOutput(unknown);
  Require(!unknownResult.ok() &&
              unknownResult.message.find("unexpected_mask_bits=0x400") !=
                  std::string::npos,
          "unknown initialization bits were accepted or not diagnosed");

  ledger.shockBackgroundGeneration = 3;
  Require(!ValidateInitializationForOutput(ledger).ok(),
          "stale shock generation passed initialization");
}

void INIT3D02() {
  CategoricalFieldIdentity a{1, 1, true, 8};
  Require(ValidateOneSidedInterpolationStencil({a, a}).ok(),
          "one-sided categorical halo stencil was rejected");
  auto b = a;
  b.sector = -1;
  Require(!ValidateOneSidedInterpolationStencil({a, b}).ok(),
          "halo/interpolation stencil crossed an HCS");
}

void INIT3D03() {
  const double sentinel = -9.87654321e299;
  auto empty = MakeFiniteOutput(std::numeric_limits<double>::quiet_NaN(),
      ValueValidity::EmptyPopulation, false, sentinel);
  Require(std::isfinite(empty.value) && empty.value == sentinel &&
          !empty.particlesPresent,
          "empty initialization output emitted NaN or lost validity");
}

void INIT3D04() {
  const auto domain = Domain();
  Require(ClassifyRadialDomain(domain, {0.999, 0.0, 0.0}) ==
              RadialDomainLocation::SolarInterior &&
          ClassifyRadialDomain(domain, {1.001, 0.0, 0.0}) ==
              RadialDomainLocation::PhysicalDomain,
          "internal solar sphere interaction was not sharply registered");
}

void NAT3D13() {
  auto invalid = MakeFiniteOutput(std::numeric_limits<double>::infinity(),
      ValueValidity::Inapplicable, false, -1.0e300);
  Require(std::isfinite(invalid.value) &&
          invalid.validity == ValueValidity::Inapplicable,
          "nullable runtime value did not use finite sentinel plus type");
}

void RUN3D02() {
  auto ledger = CompleteLedger();
  ledger.completedMask &= ~static_cast<std::uint32_t>(
      InitializationCondition::Observers);
  Require(!ValidateInitializationForOutput(ledger).ok(),
          "initialization-only stopped before observer initialization");
  MarkInitialized(&ledger, InitializationCondition::Observers);
  Require(ValidateInitializationForOutput(ledger).ok(),
          "initialization-only did not stop after all ten conditions");
}

void OBS3D01() {
  auto first = VolumeObserver();
  auto second = first;
  second.stableId = "solar-orbiter";
  auto a = AccumulateObserver(first, {Particle()});
  auto b = AccumulateObserver(second, {Particle()});
  Require(a.ok() && b.ok() && a.value[0].weightedCount ==
          b.value[0].weightedCount && first.stableId != second.stableId,
          "multiple stable-ID observers were not independent");
}

void OBS3D02() {
  auto linear = BuildEnergyEdges(EnergyGrid::Linear, 1.0, 9.0, 2);
  auto log = BuildEnergyEdges(EnergyGrid::Logarithmic, 1.0, 100.0, 2);
  auto explicitGrid = BuildEnergyEdges(
      EnergyGrid::Explicit, 0.0, 0.0, 0, {1.0, 3.0, 10.0});
  Require(linear.ok() && log.ok() && explicitGrid.ok() &&
          Close(linear.value[1], 5.0) && Close(log.value[1], 10.0) &&
          explicitGrid.value[1] == 3.0,
          "linear/log/explicit energy edges are incorrect");
}

void OBS3D03() {
  auto bins = AccumulateObserver(VolumeObserver(), {});
  Require(bins.ok() && !bins.value[0].particlesPresent &&
          bins.value[0].weightedCount == 0.0 &&
          std::isfinite(bins.value[0].intensity),
          "empty observer bin is not finite and explicitly empty");
}

void OBS3D04() {
  auto observer = VolumeObserver();
  auto bins = AccumulateObserver(observer, {Particle()});
  const double volume = 4.0 * Constants::kPi * 8.0 / 3.0;
  Require(bins.ok() && Close(bins.value[0].intensity, 8.0 / volume),
          "volume residence estimator has wrong measure/time normalization");
}

void OBS3D05() {
  auto observer = VolumeObserver();
  observer.energyPerNucleon = true;
  auto ion = Particle({1.0e6, 0.0, 0.0});
  ion.nucleonCount = 4;
  auto accepted = AccumulateObserver(observer, {ion});
  ion.nucleonCount = 0;
  auto rejected = AccumulateObserver(observer, {ion});
  Require(accepted.ok() && !rejected.ok(),
          "per-nucleon observer did not require a validated mass number");
}

void OBS3D06() {
  auto stationary = VolumeObserver();
  auto moving = stationary;
  moving.velocityMPerS = {900.0, 0.0, 0.0};
  auto particle = Particle({1000.0, 0.0, 0.0});
  auto a = AccumulateObserver(stationary, {particle});
  auto b = AccumulateObserver(moving, {particle});
  Require(a.ok() && b.ok() && a.value[0].particlesPresent &&
          !b.value[0].particlesPresent,
          "energy binning did not use the detector rest frame");
}

void OBS3D07() {
  auto disk = VolumeObserver();
  disk.geometry = ObserverGeometry::SurfaceDisk;
  disk.surfaceNormal = {1.0, 0.0, 0.0};
  auto sphere = disk;
  sphere.geometry = ObserverGeometry::SurfaceSphere;
  auto diskBins = AccumulateObserver(disk, {Particle()});
  auto sphereBins = AccumulateObserver(sphere, {Particle()});
  Require(diskBins.ok() && sphereBins.ok() &&
          Close(diskBins.value[0].intensity,
                4.0 * sphereBins.value[0].intensity),
          "disk/sphere surface measures were not discriminated");
}

}  // namespace

void RegisterStage8(Registry* tests) {
  (*tests)["BND3D01"] = BND3D01; (*tests)["BND3D02"] = BND3D02;
  (*tests)["MESH3D01"] = MESH3D01;
  (*tests)["COR3D01"] = COR3D01; (*tests)["COR3D02"] = COR3D02;
  (*tests)["COR3D03"] = COR3D03;
  (*tests)["TIM3D01"] = TIM3D01; (*tests)["TIM3D02"] = TIM3D02;
  (*tests)["POP3D01"] = POP3D01; (*tests)["POP3D02"] = POP3D02;
  (*tests)["POP3D03"] = POP3D03; (*tests)["MPI3D01"] = MPI3D01;
  (*tests)["INIT3D01"] = INIT3D01; (*tests)["INIT3D02"] = INIT3D02;
  (*tests)["INIT3D03"] = INIT3D03; (*tests)["INIT3D04"] = INIT3D04;
  (*tests)["NAT3D13"] = NAT3D13; (*tests)["RUN3D02"] = RUN3D02;
  (*tests)["OBS3D01"] = OBS3D01; (*tests)["OBS3D02"] = OBS3D02;
  (*tests)["OBS3D03"] = OBS3D03; (*tests)["OBS3D04"] = OBS3D04;
  (*tests)["OBS3D05"] = OBS3D05; (*tests)["OBS3D06"] = OBS3D06;
  (*tests)["OBS3D07"] = OBS3D07;
}

}  // namespace SCCMTest
