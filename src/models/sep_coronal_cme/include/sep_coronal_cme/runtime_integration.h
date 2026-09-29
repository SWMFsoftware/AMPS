#ifndef SEP_CORONAL_CME_RUNTIME_INTEGRATION_H
#define SEP_CORONAL_CME_RUNTIME_INTEGRATION_H

#include "sep_coronal_cme/vector_math.h"
#include "sep_status.h"

#include <cstdint>
#include <string>
#include <vector>

namespace SEP { namespace CoronalCME {

enum class RadialDomainLocation { SolarInterior, PhysicalDomain, Escaped };

struct SphericalDomain {
  double solarRadiusM = 0.0;
  double outerRadiusM = 0.0;
  Vec3 cartesianHalfExtentM;
};
Core::Status ValidateSphericalDomain(const SphericalDomain& domain);
RadialDomainLocation ClassifyRadialDomain(const SphericalDomain& domain,
                                          Vec3 positionM);

struct RefinementControls {
  double globalCellSizeM = 0.0;
  double solarCellSizeM = 0.0;
  double solarDecayLengthM = 0.0;
  double tubeCellSizeM = 0.0;
  double tubeDecayLengthM = 0.0;
};
// The two independent requested sizes asymptote monotonically to the global
// size; the mesh receives their minimum (the finer physical requirement).
Core::Result<double> RequestedCellSize(
    const RefinementControls& controls, double distanceFromSolarSurfaceM,
    double distanceFromAuthoritativeCenterlineM);

struct AxisAlignedBlock {
  std::uint64_t stableId = 0;
  Vec3 minimumM;
  Vec3 maximumM;
  bool active = false;
};
bool BlockIntersectsBufferedPolyline(const AxisAlignedBlock& block,
                                     const std::vector<Vec3>& centerlineM,
                                     double bufferRadiusM);
Core::Status ValidateFaceConnectedBlocks(
    const std::vector<AxisAlignedBlock>& blocks,
    const std::vector<std::uint64_t>& requiredStableIds,
    double coordinateToleranceM = 0.0);

struct CategoricalFieldIdentity {
  int region = 0;
  int sector = 0;
  bool mappingValid = false;
  std::uint64_t interfaceIdentity = 0;
};
Core::Status ValidateOneSidedInterpolationStencil(
    const std::vector<CategoricalFieldIdentity>& identities);

struct IdealHcsCrossingState {
  Vec3 positionM;
  Vec3 physicalVelocityMPerS;
  Vec3 physicalMomentumKgMPerS;
  int magneticSector = 1;
  int outwardWaveLabel = 1;
  int inwardWaveLabel = -1;
};
Core::Result<IdealHcsCrossingState> CrossIdealHcs(
    const IdealHcsCrossingState& before, bool finiteThicknessDriftRequested);

struct TimeStepControls {
  double fixedUpperBoundS = 0.0;
  double gyroAccuracy = 0.0;
  double spatialAccuracy = 0.0;
  double scatteringAccuracy = 0.0;
};
struct LocalStepInputs {
  double gyrofrequencyRadPerS = 0.0;
  double speedMPerS = 0.0;
  double cellSizeM = 0.0;
  double scatteringTimeS = 0.0;
};
Core::Result<double> ComputeLocalTimeStep(const TimeStepControls& controls,
                                          const LocalStepInputs& inputs);

struct SpeciesWeight {
  std::string stableSpeciesId;
  int compiledSlot = -1;
  double baseWeight = 0.0;
};
Core::Status ValidateAllSpeciesWeights(
    const std::vector<SpeciesWeight>& weights, int compiledSpeciesCount);

struct PopulationParticle {
  std::uint64_t id = 0;
  int speciesSlot = -1;
  std::uint64_t backgroundGeneration = 0;
  std::string sourceLabel;
  double representedWeight = 0.0;
  Vec3 momentumKgMPerS;
  double kineticEnergyJ = 0.0;
};
struct PopulationMoments {
  double representedWeight = 0.0;
  Vec3 representedMomentumKgMPerS;
  double representedKineticEnergyJ = 0.0;
};
PopulationMoments AccumulatePopulationMoments(
    const std::vector<PopulationParticle>& particles);
Core::Result<std::vector<PopulationParticle>> SplitParticleConservatively(
    const PopulationParticle& particle, int copies);
Core::Result<PopulationParticle> MergeIdenticalStateParticles(
    const std::vector<PopulationParticle>& particles,
    std::uint64_t mergedId);

enum class InitializationCondition : std::uint32_t {
  Mesh = 1U << 0,
  Boundary = 1U << 1,
  Background = 1U << 2,
  Turbulence = 1U << 3,
  Shock = 1U << 4,
  Halo = 1U << 5,
  Species = 1U << 6,
  TimeStep = 1U << 7,
  Observers = 1U << 8,
  OutputDictionary = 1U << 9
};
struct InitializationLedger {
  std::uint32_t completedMask = 0;
  std::uint64_t backgroundGeneration = 0;
  std::uint64_t shockBackgroundGeneration = 0;
};
void MarkInitialized(InitializationLedger* ledger,
                     InitializationCondition condition);
Core::Status ValidateInitializationForOutput(
    const InitializationLedger& ledger);

enum class ValueValidity { Valid, EmptyPopulation, Inapplicable };
struct FiniteOutputValue {
  double value = 0.0;
  ValueValidity validity = ValueValidity::Valid;
  bool particlesPresent = false;
};
FiniteOutputValue MakeFiniteOutput(double value, ValueValidity validity,
                                   bool particlesPresent,
                                   double finiteSentinel);

enum class EnergyGrid { Linear, Logarithmic, Explicit };
Core::Result<std::vector<double>> BuildEnergyEdges(
    EnergyGrid grid, double minimumJ, double maximumJ, int channels,
    const std::vector<double>& explicitEdgesJ = {});

enum class ObserverGeometry { VolumeSphere, SurfaceDisk, SurfaceSphere };
struct ObserverDefinition {
  std::string stableId;
  ObserverGeometry geometry = ObserverGeometry::VolumeSphere;
  Vec3 centerM;
  Vec3 velocityMPerS;
  Vec3 surfaceNormal;
  double radiusM = 0.0;
  std::vector<double> energyEdgesJ;
  bool energyPerNucleon = false;
};
struct ObserverParticle {
  Vec3 positionM;
  Vec3 velocityMPerS;
  double massKg = 0.0;
  double representedWeight = 0.0;
  double residenceTimeS = 0.0;
  int nucleonCount = 0;
  bool crossedSurface = false;
  int crossingSense = 1;
};
struct ObserverBin {
  double weightedCount = 0.0;
  double intensity = 0.0;
  bool particlesPresent = false;
  bool valueValid = true;
};
Core::Result<std::vector<ObserverBin>> AccumulateObserver(
    const ObserverDefinition& observer,
    const std::vector<ObserverParticle>& particles);

// Sorted compensated reduction is used by verification adapters to prove
// that rank partition and request order cannot change a physical sum.
double DeterministicPhysicalSum(
    const std::vector<std::pair<std::uint64_t, double>>& stableContributions);

} }  // namespace SEP::CoronalCME

#endif  // SEP_CORONAL_CME_RUNTIME_INTEGRATION_H
