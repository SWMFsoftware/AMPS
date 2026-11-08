/*
 * Native one-rank qualification of the legacy spherical exosphere source.
 *
 * This fixture intentionally drives Exosphere::SourceProcesses::
 * InjectionBoundaryModel through a real 3-D AMR mesh, particle buffer,
 * internal sphere, source sampling buffers, and production mover.  The input
 * deck compiles all target-independent spherical source modes currently
 * exposed by exosphere.pl.  SI units are used throughout.
 */
#include "pic.h"
#include "constants.h"
#include "Exosphere.h"

#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <set>

namespace {

const double kSphereRadiusM = 0.40;
const double kTimeStepS = 1.0e-5;
const double kBaseWeight = 1.0e4;
const double kInitialInventoryPerElement = 2.5e11;

double Resolution(double*) { return 0.25; }
double UnitLoad(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>*) { return 1.0; }

void ConfigureLeafBlocks(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node) {
  if (node == NULL) return;
  if (node->lastBranchFlag() == _BOTTOM_BRANCH_TREE_) {
    if (node->block != NULL) {
      for (int species = 0; species < PIC::nTotalSpecies; ++species) {
        node->block->SetLocalTimeStep(kTimeStepS, species);
        node->block->SetLocalParticleWeight(kBaseWeight, species);
      }
    }
    return;
  }
  for (int child = 0; child < (1 << DIM); ++child)
    ConfigureLeafBlocks(node->downNode[child]);
}

void SetIdentityStateTransforms() {
  // With SPICE disabled, explicitly state that the body-fixed and simulation
  // frames coincide and have zero relative velocity.  AMPS stores both maps
  // as row-major 6x6 state transforms [x;v] -> [x';v'].
  for (int i = 0; i < 6; ++i) {
    for (int j = 0; j < 6; ++j) {
      Exosphere::OrbitalMotion::IAU_to_SO_TransformationMartix[i][j] =
          (i == j) ? 1.0 : 0.0;
      Exosphere::OrbitalMotion::SO_to_IAU_TransformationMartix[i][j] =
          (i == j) ? 1.0 : 0.0;
    }
  }
  Exosphere::OrbitalMotion::SunDirection_IAU_OBJECT[0] = 1.0;
  Exosphere::OrbitalMotion::SunDirection_IAU_OBJECT[1] = 0.0;
  Exosphere::OrbitalMotion::SunDirection_IAU_OBJECT[2] = 0.0;
}

double SumPopulation(cInternalSphericalData* sphere) {
  double result = 0.0;
  for (int element = 0;
       element < PIC::BC::InternalBoundary::Sphere::TotalSurfaceElementNumber;
       ++element)
    result += sphere->SurfaceElementPopulation[0][element];
  return result;
}

double SumSurfaceBuffer(double** buffer) {
  double result = 0.0;
  for (int element = 0;
       element < PIC::BC::InternalBoundary::Sphere::TotalSurfaceElementNumber;
       ++element)
    result += buffer[0][element];
  return result;
}

void CollectCommittedParticles(
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node,
    std::set<long int>& particles) {
  if (node == NULL) return;
  if (node->lastBranchFlag() == _BOTTOM_BRANCH_TREE_) {
    if (node->block == NULL) return;
    const int cells = _BLOCK_CELLS_X_ * _BLOCK_CELLS_Y_ * _BLOCK_CELLS_Z_;
    for (int cell = 0; cell < cells; ++cell) {
      long int particle = node->block->FirstCellParticleTable[cell];
      while (particle != -1) {
        particles.insert(particle);
        particle = PIC::ParticleBuffer::GetNext(particle);
      }
    }
    return;
  }
  for (int child = 0; child < (1 << DIM); ++child)
    CollectCommittedParticles(node->downNode[child], particles);
}

}  // namespace

int main(int, char**) {
  PIC::Init_BeforeParser();
  Exosphere::xObjectRadial = _AU_;
  SetIdentityStateTransforms();
  Exosphere::Init_BeforeParser();
  PIC::Mesh::initCellSamplingDataBuffer();

  cInternalSphericalData::SetGeneralSurfaceMeshParameters(4, 8);
  PIC::BC::InternalBoundary::Sphere::Init();
  cInternalBoundaryConditionsDescriptor descriptor =
      PIC::BC::InternalBoundary::Sphere::RegisterInternalSphere();
  cInternalSphericalData* sphere =
      static_cast<cInternalSphericalData*>(descriptor.BoundaryElement);
  double centre[3] = {0.0, 0.0, 0.0};
  sphere->SetSphereGeometricalParameters(centre, kSphereRadiusM);
  sphere->localResolution = Resolution;
  Exosphere::Planet = sphere;
  sphere->Allocate<cInternalSphericalData>(
      PIC::nTotalSpecies,
      PIC::BC::InternalBoundary::Sphere::TotalSurfaceElementNumber,
      _EXOSPHERE__SOURCE_MAX_ID_VALUE_, sphere);

  for (int element = 0;
       element < PIC::BC::InternalBoundary::Sphere::TotalSurfaceElementNumber;
       ++element)
    sphere->SurfaceElementPopulation[0][element] =
        kInitialInventoryPerElement;

  double xmin[3] = {-1.0, -1.0, -1.0};
  double xmax[3] = { 1.0,  1.0,  1.0};
  PIC::Mesh::mesh->AllowBlockAllocation = false;
  PIC::Mesh::mesh->init(xmin, xmax, Resolution);
  PIC::Mesh::mesh->SetNodeNeibResolutionLevelLimit();
  PIC::Mesh::mesh->CreateBottomBranchNodeList(PIC::Mesh::mesh->rootTree);
  PIC::Mesh::mesh->SetParallelLoadMeasure(UnitLoad);
  PIC::Mesh::mesh->CreateNewParallelDistributionLists();
  PIC::Mesh::mesh->AllowBlockAllocation = true;
  PIC::Mesh::mesh->AllocateTreeBlocks();
  PIC::Mesh::mesh->InitCellMeasure();
  ConfigureLeafBlocks(PIC::Mesh::mesh->rootTree);
  PIC::DomainBlockDecomposition::UpdateBlockTable();

  PIC::Init_AfterParser();
  PIC::Mover::Init();
  PIC::ParticleBuffer::Init(8192);
  PIC::ParticleWeightTimeStep::GlobalParticleWeight[0] = kBaseWeight;
  PIC::ParticleWeightTimeStep::GlobalTimeStep[0] = kTimeStepS;

  // The three surface-weighted source generators require their discrete
  // element distributions to be bound to this sphere.  Production target
  // adapters perform the same binding after constructing their surface.
  const int last_element =
      PIC::BC::InternalBoundary::Sphere::TotalSurfaceElementNumber - 1;
  Exosphere::SourceProcesses::PhotonStimulatedDesorption::
      SurfaceInjectionDistribution[0].SetLimits(
          0, last_element,
          Exosphere::SourceProcesses::PhotonStimulatedDesorption::
              GetSurfaceElementProductionRate);
  Exosphere::SourceProcesses::ThermalDesorption::
      SurfaceInjectionDistribution[0].SetLimits(
          0, last_element,
          Exosphere::SourceProcesses::ThermalDesorption::
              GetSurfaceElementProductionRate);
  Exosphere::SourceProcesses::SolarWindSputtering::
      SurfaceInjectionDistribution[0].SetLimits(
          0, last_element,
          Exosphere::SourceProcesses::SolarWindSputtering::
              GetSurfaceElementProductionRate);

  Exosphere::Init_AfterParser();

  // PSD and thermal totals are ordinarily refreshed by the target's surface
  // exchange step.  Bind them to equal physical rates here while retaining
  // their real element-shape functions, so the test observes every compiled
  // dispatcher branch without making its oracle depend on random counts.
  Exosphere::SourceProcesses::PhotonStimulatedDesorption::SourceRate[0] =
      3.2e10;
  Exosphere::SourceProcesses::ThermalDesorption::SourceRate[0] = 3.2e10;

  std::cout << "[CASE] RUN  CORE-EXOSPHERE-SPHERE-INTEGRATION-001 "
               "legacy spherical sources, native AMR ownership and ledgers"
            << std::endl;

  bool pass = Exosphere::nTotalSourceProcesses == 5;
  const double initial_population = SumPopulation(sphere);
  long int injected = 0;
  // Four independent physical intervals make each approximately equal source
  // mode observable in the fixed stream while remaining far below the 8192
  // record capacity.  Exact source closure below uses realized weights.
  for (int interval = 0; interval < 4; ++interval)
    injected += Exosphere::SourceProcesses::InjectionBoundaryModel(
        0, _INTERNAL_BOUNDARY_TYPE_SPHERE_, sphere);
  PIC::Mover::CommitTempParticleMovingListsToFirstCellParticleTable();

  std::set<long int> committed;
  CollectCommittedParticles(PIC::Mesh::mesh->rootTree, committed);
  pass = pass && injected > 0 &&
      static_cast<long int>(committed.size()) ==
          PIC::ParticleBuffer::GetAllPartNum();

  long int source_counts[1 + _EXOSPHERE__SOURCE_MAX_ID_VALUE_] = {};
  for (std::set<long int>::const_iterator it = committed.begin();
       it != committed.end(); ++it) {
    PIC::ParticleBuffer::byte* data =
        PIC::ParticleBuffer::GetParticleDataPointer(*it);
    double position[3];
    PIC::ParticleBuffer::GetX(position, data);
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* owner =
        PIC::Mesh::mesh->findTreeNode(position, NULL);
    const int source = Exosphere::Sampling::GetParticleSourceID(data);
    const double correction =
        PIC::ParticleBuffer::GetIndividualStatWeightCorrection(data);
    if (owner == NULL || owner->block == NULL ||
        owner->Thread != PIC::Mesh::mesh->ThisThread || source < 0 ||
        source > _EXOSPHERE__SOURCE_MAX_ID_VALUE_ ||
        !(correction > 0.0) || !std::isfinite(correction)) {
      pass = false;
      continue;
    }
    ++source_counts[source];
  }

  const int compiled_sources[] = {
      _EXOSPHERE_SOURCE__ID__IMPACT_VAPORIZATION_,
      _EXOSPHERE_SOURCE__ID__PHOTON_STIMULATED_DESPRPTION_,
      _EXOSPHERE_SOURCE__ID__THERMAL_DESORPTION_,
      _EXOSPHERE_SOURCE__ID__SOLAR_WIND_SPUTTERING_,
      _EXOSPHERE_SOURCE__ID__VERTICAL_INJECTION_};
  for (std::size_t i = 0;
       i < sizeof(compiled_sources) / sizeof(compiled_sources[0]); ++i)
    pass = pass && source_counts[compiled_sources[i]] > 0;

  double calculated_rate = 0.0;
  double sampled_rate = 0.0;
  for (int source = 0; source <= _EXOSPHERE__SOURCE_MAX_ID_VALUE_; ++source) {
    calculated_rate += Exosphere::Sampling::CalculatedSourceRate[0][source];
    for (int element = 0;
         element < PIC::BC::InternalBoundary::Sphere::TotalSurfaceElementNumber;
         ++element)
      sampled_rate += sphere->SampleSpeciesSurfaceSourceRate[0][element][source];
  }
  const double rate_scale = std::fmax(1.0, std::fabs(calculated_rate));
  pass = pass &&
      std::fabs(sampled_rate - calculated_rate) <= 2.0e-13 * rate_scale &&
      std::fabs(PIC::BC::ParticleProductionRate[0] - calculated_rate) <=
          2.0e-13 * rate_scale;

  // Exercise the production accommodation/adsorption path with a committed
  // particle.  The application closure returns sticking=1 and retained
  // adsorption fraction 0.75, giving exact represented-weight expectations.
  const long int impact_particle = committed.empty() ? -1 : *committed.begin();
  double impact_weight = 0.0;
  if (impact_particle >= 0) {
    impact_weight = kBaseWeight *
        PIC::ParticleBuffer::GetIndividualStatWeightCorrection(
            PIC::ParticleBuffer::GetParticleDataPointer(impact_particle));
    double x[3] = {kSphereRadiusM, 0.0, 0.0};
    double v[3] = {-100.0, 0.0, 0.0};
    double remaining_time = 0.0;
    const int interaction =
        Exosphere::SurfaceInteraction::
            ParticleSphereInteraction_SurfaceAccomodation(
                0, impact_particle, x, v, remaining_time,
                PIC::Mesh::mesh->rootTree, sphere);
    pass = pass && interaction == _PARTICLE_DELETED_ON_THE_FACE_;
  }

  const double adsorption =
      SumSurfaceBuffer(sphere->SurfaceElementAdsorptionFluxDOWN);
  const double desorption =
      SumSurfaceBuffer(sphere->SurfaceElementDesorptionFluxUP);
  pass = pass &&
      std::fabs(Exosphere::Sampling::TotalPlanetReturnFlux[0] -
                impact_weight) <= 1.0e-13 * std::fmax(1.0, impact_weight) &&
      std::fabs(Exosphere::Sampling::PlanetSurfaceStickingRate[0] -
                impact_weight) <= 1.0e-13 * std::fmax(1.0, impact_weight) &&
      std::fabs(adsorption - 0.75 * impact_weight) <=
          1.0e-13 * std::fmax(1.0, impact_weight);

  Exosphere::ExchangeSurfaceAreaDensity();
  const double final_population = SumPopulation(sphere);
  const double population_scale = std::fmax(1.0, initial_population);
  pass = pass &&
      std::fabs(final_population -
                (initial_population + adsorption - desorption)) <=
          2.0e-13 * population_scale &&
      SumSurfaceBuffer(sphere->SurfaceElementAdsorptionFluxDOWN) == 0.0 &&
      SumSurfaceBuffer(sphere->SurfaceElementDesorptionFluxUP) == 0.0;

  std::cout << std::setprecision(17)
            << "exosphere_native injected=" << injected
            << " committed=" << committed.size()
            << " calculated_rate_s-1=" << calculated_rate
            << " sampled_rate_s-1=" << sampled_rate
            << " adsorption=" << adsorption
            << " desorption=" << desorption
            << " population_delta=" << final_population - initial_population
            << std::endl;
  std::cout << "exosphere_source_counts impact="
            << source_counts[_EXOSPHERE_SOURCE__ID__IMPACT_VAPORIZATION_]
            << " psd="
            << source_counts[
                   _EXOSPHERE_SOURCE__ID__PHOTON_STIMULATED_DESPRPTION_]
            << " thermal="
            << source_counts[_EXOSPHERE_SOURCE__ID__THERMAL_DESORPTION_]
            << " sputtering="
            << source_counts[_EXOSPHERE_SOURCE__ID__SOLAR_WIND_SPUTTERING_]
            << " vertical="
            << source_counts[_EXOSPHERE_SOURCE__ID__VERTICAL_INJECTION_]
            << std::endl;
  std::cout << "[CORE-EXOSPHERE-SPHERE-INTEGRATION-001] "
            << (pass ? "PASS" : "FAIL") << std::endl;

  MPI_Finalize();
  return pass ? EXIT_SUCCESS : EXIT_FAILURE;
}
