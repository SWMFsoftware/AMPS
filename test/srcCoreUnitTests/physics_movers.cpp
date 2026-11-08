#include "pic.h"
#include "PhotolyticReactions.h"
#include "src/test/native/test_harness.h"
#include "core_test_state.h"
#include "core_test_mesh.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <vector>

namespace {

double UniformResolution(double*) {
  // A constant physical resolution produces a small uniform tree.  The test
  // resolves the owning leaf explicitly, so its oracle is independent of the
  // generator's block-refinement threshold while cut-cell and adaptive
  // resolution physics remain absent.
  return 2.0;
}

double UnitLoad(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>*) { return 1.0; }

int surface_interaction_count = 0;
cInternalSphericalData* test_sphere = NULL;

int SpecularSphereInteraction(
    int, long int, double* position, double* velocity, double&,
    void*, void* sphere_data) {
  // The callback receives a point on the surface in the simulation frame.
  // Specular reflection is v' = v - 2(v.n)n for the outward unit normal n;
  // it conserves speed and reverses only the normal component. The callback
  // does not consume the remaining flight time; the production mover resumes
  // it after the rejected-on-face result.
  cInternalSphericalData* sphere =
      static_cast<cInternalSphericalData*>(sphere_data);
  double* origin = NULL;
  double radius = 0.0;
  sphere->GetSphereGeometricalParameters(origin, radius);
  double normal[3], norm = 0.0;
  for (int dimension = 0; dimension < 3; ++dimension) {
    normal[dimension] = position[dimension] - origin[dimension];
    norm += normal[dimension] * normal[dimension];
  }
  norm = std::sqrt(norm);
  AMPS_EXPECT_TRUE(radius > 0.0 && norm > 0.0);
  double normal_velocity = 0.0;
  for (int dimension = 0; dimension < 3; ++dimension) {
    normal[dimension] /= norm;
    normal_velocity += velocity[dimension] * normal[dimension];
  }
  for (int dimension = 0; dimension < 3; ++dimension)
    velocity[dimension] -= 2.0 * normal_velocity * normal[dimension];
  ++surface_interaction_count;
  return _PARTICLE_REJECTED_ON_THE_FACE_;
}

cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* MoverNode(const double* position) {
  static bool mesh_initialized = false;
  if (!mesh_initialized) {
    const double xmin[3] = {-1.0, -1.0, -1.0};
    const double xmax[3] = {1.0, 1.0, 1.0};
    // Register one analytic surface before mesh construction so the AMR tree
    // carries the same internal-boundary descriptors used by applications.
    // The radius is small enough that the existing off-axis ballistic cases
    // do not intersect it; the dedicated surface case crosses it normally.
    SetGeneralSphereSurfaceMeshParameters(8, 16);
    PIC::BC::InternalBoundary::Sphere::Init();
    cInternalBoundaryConditionsDescriptor descriptor =
        PIC::BC::InternalBoundary::Sphere::RegisterInternalSphere();
    test_sphere = static_cast<cInternalSphericalData*>(descriptor.BoundaryElement);
    double origin[3] = {0.0, 0.0, 0.0};
    test_sphere->SetSphereGeometricalParameters(origin, 0.2);
    test_sphere->ParticleSphereInteraction = SpecularSphereInteraction;
    test_sphere->InjectionRate = NULL;
    test_sphere->InjectionBoundaryCondition = NULL;
    PIC::Mesh::mesh->AllowBlockAllocation = false;
    PIC::Mesh::mesh->init(const_cast<double*>(xmin), const_cast<double*>(xmax),
                          UniformResolution);
    PIC::Mesh::mesh->buildMesh();
    PIC::Mesh::mesh->SetParallelLoadMeasure(UnitLoad);
    PIC::Mesh::mesh->CreateNewParallelDistributionLists();
    PIC::Mesh::mesh->AllowBlockAllocation = true;
    PIC::Mesh::mesh->AllocateTreeBlocks();
    PIC::Mesh::mesh->InitCellMeasure();
    // Whole-mesh PIC modules do not traverse the AMR linked list directly;
    // they use DomainBlockDecomposition::BlockTable.  Native applications
    // refresh that table during initialization/load balancing.  The fixture
    // must do the same after block allocation, otherwise chemistry/collision
    // drivers execute a zero-block no-op while cell-local tests still appear
    // healthy.  This table contains the eight one-rank owner leaves and has no
    // particle or physical-state side effects.
    PIC::DomainBlockDecomposition::UpdateBlockTable();
    mesh_initialized = true;
  }

  // AMR internal nodes do not own data blocks.  The resolution request may
  // refine the root, so always resolve the leaf that contains this particle
  // rather than assuming rootTree is an allocated owner.
  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node =
      PIC::Mesh::mesh->findTreeNode(const_cast<double*>(position));
  AMPS_EXPECT_TRUE(node != NULL && node->block != NULL);
  AMPS_EXPECT_EQ(node->Thread, PIC::ThisThread);
  return node;
}

using AMPS::CoreTests::Fixture::EmptyListHeadScope;
using AMPS::CoreTests::Fixture::MeshCellStateScope;
using AMPS::CoreTests::Fixture::ParticleBufferScope;
using AMPS::CoreTests::Fixture::PhotolysisLifetimeScope;

void RunBallisticMover(bool boris) {
  const double initial_position[3] = {-0.35, 0.65, -0.10};
  // Keep the endpoint away from AMR faces.  A point exactly on a face has two
  // geometrically valid adjacent leaves and is a separate ownership-policy
  // test, not part of this ballistic mover oracle.
  const double velocity[3] = {0.40, -0.30, 0.16};
  const double dt = 0.5;
  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node = MoverNode(initial_position);
  double expected[3];
  for (int dimension = 0; dimension < 3; ++dimension)
    expected[dimension] = initial_position[dimension] + velocity[dimension] * dt;
  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* destination =
      PIC::Mesh::mesh->findTreeNode(expected, node);
  AMPS_EXPECT_TRUE(destination != NULL && destination->block != NULL);
  int i, j, k;
  AMPS_EXPECT_TRUE(PIC::Mesh::mesh->FindCellIndex(
                       expected, i, j, k, destination, false) >= 0);
  long int& destination_head = destination->block->tempParticleMovingListTable[
      i + _BLOCK_CELLS_X_ * (j + _BLOCK_CELLS_Y_ * k)];
  ParticleBufferScope particles(4096);
  EmptyListHeadScope destination_guard(destination_head);
  const long int particle = PIC::ParticleBuffer::GetNewParticle(true);
  PIC::ParticleBuffer::SetI(_H2O_SPEC_, particle);  // neutral: q/m=0
  PIC::ParticleBuffer::SetX(const_cast<double*>(initial_position), particle);
  PIC::ParticleBuffer::SetV(const_cast<double*>(velocity), particle);

  const int result = boris ? PIC::Mover::Boris(particle, dt, node)
                           : PIC::Mover::Simple(particle, dt, node);
  AMPS_EXPECT_EQ(result, _PARTICLE_MOTION_FINISHED_);
  const double* final_position = PIC::ParticleBuffer::GetX(particle);
  const double* final_velocity = PIC::ParticleBuffer::GetV(particle);
  for (int dimension = 0; dimension < 3; ++dimension) {
    AMPS_EXPECT_NEAR(final_position[dimension], expected[dimension], 3.0e-15);
    AMPS_EXPECT_NEAR(final_velocity[dimension], velocity[dimension], 3.0e-15);
  }

  AMPS_EXPECT_EQ(destination_head, particle);
  AMPS_EXPECT_EQ(PIC::ParticleBuffer::GetNext(particle), -1);
}

struct PopulationTotals {
  double momentum[3];
  double translational;
  double rotational;
};

PopulationTotals MeasurePopulation(long int first) {
  PopulationTotals totals = {{0.0, 0.0, 0.0}, 0.0, 0.0};
  for (long int particle = first; particle != -1;
       particle = PIC::ParticleBuffer::GetNext(particle)) {
    const int species = PIC::ParticleBuffer::GetI(particle);
    const double mass = PIC::MolecularData::GetMass(species);
    const double* velocity = PIC::ParticleBuffer::GetV(particle);
    for (int dimension = 0; dimension < 3; ++dimension) {
      totals.momentum[dimension] += mass * velocity[dimension];
      totals.translational +=
          0.5 * mass * velocity[dimension] * velocity[dimension];
    }
    totals.rotational += PIC::IDF::LB::GetRotE(
        PIC::ParticleBuffer::GetParticleDataPointer(particle));
  }
  return totals;
}

int PopulationCell(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node,
                   const double* position, int& i, int& j, int& k) {
  const int local = PIC::Mesh::mesh->FindCellIndex(
      const_cast<double*>(position), i, j, k, node, false);
  AMPS_EXPECT_TRUE(local >= 0);
  return i + _BLOCK_CELLS_X_ * (j + _BLOCK_CELLS_Y_ * k);
}

}  // namespace

namespace AMPS {
namespace CoreTests {

void LinkMoverTests() {}

cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* MoverFixtureNode(
    const double* position) {
  return ::MoverNode(position);
}

cInternalSphericalData* MoverFixtureSphere() {
  // Force lazy initialization before publishing the registered descriptor.
  const double outside[3] = {0.6, 0.6, 0.6};
  ::MoverNode(outside);
  return ::test_sphere;
}

}  // namespace CoreTests
}  // namespace AMPS

AMPS_TEST_CASE("PHYS-MOVER-SIMPLE-001",
               "simple mover reproduces neutral ballistic motion and cell insertion") {
  RunBallisticMover(false);
}

AMPS_TEST_CASE("PHYS-MOVER-BORIS-001",
               "Boris mover reduces to neutral ballistic motion and cell insertion") {
  // With q=0, E=B=0, and no target gravity, the exact equations are dx/dt=v
  // and dv/dt=0.  This analytic limit tests the production Boris dispatch,
  // state update, AMR lookup, and destination-list path without conflating
  // field interpolation error.  Uniform-B convergence remains future work.
  RunBallisticMover(true);
}

AMPS_TEST_CASE("PHYS-MOVER-SURFACE-001",
               "trajectory mover resolves and specularly reflects at a sphere") {
  const double initial_position[3] = {-0.60, 0.0, 0.0};  // m
  const double initial_velocity[3] = {1.0, 0.0, 0.0};    // m/s
  const double dt = 1.0;                                  // s
  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node =
      MoverNode(initial_position);
  AMPS_EXPECT_TRUE(test_sphere != NULL);
  double expected_position[3] = {-0.8, 0.0, 0.0};
  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* destination =
      PIC::Mesh::mesh->findTreeNode(expected_position, node);
  AMPS_EXPECT_TRUE(destination != NULL && destination->block != NULL);
  int i, j, k;
  AMPS_EXPECT_TRUE(PIC::Mesh::mesh->FindCellIndex(
                       expected_position, i, j, k, destination, false) >= 0);
  long int& destination_head = destination->block->tempParticleMovingListTable[
      i + _BLOCK_CELLS_X_ * (j + _BLOCK_CELLS_Y_ * k)];
  ParticleBufferScope particles(4096);
  EmptyListHeadScope destination_guard(destination_head);
  const long int particle = PIC::ParticleBuffer::GetNewParticle(true);
  PIC::ParticleBuffer::SetI(_H2O_SPEC_, particle);
  PIC::ParticleBuffer::SetX(const_cast<double*>(initial_position), particle);
  PIC::ParticleBuffer::SetV(const_cast<double*>(initial_velocity), particle);
  AMPS::CoreTests::Fixture::ValueScope<int> interaction_count_scope(
      surface_interaction_count);
  surface_interaction_count = 0;

  const int result =
      PIC::Mover::UniformWeight_UniformTimeStep_noForce_TraceTrajectory(
          particle, dt, node);
  AMPS_EXPECT_EQ(result, _PARTICLE_MOTION_FINISHED_);
  AMPS_EXPECT_EQ(surface_interaction_count, 1);
  const double* position = PIC::ParticleBuffer::GetX(particle);
  const double* velocity = PIC::ParticleBuffer::GetV(particle);
  // The particle reaches x=-0.2 after 0.4 s, reflects, and covers 0.6 m
  // outward during the remaining interval. This analytic event-time oracle
  // checks intersection location, callback dispatch, and remaining-time use.
  AMPS_EXPECT_NEAR(position[0], -0.8, 2.0e-12);
  AMPS_EXPECT_NEAR(position[1], 0.0, 2.0e-12);
  AMPS_EXPECT_NEAR(position[2], 0.0, 2.0e-12);
  AMPS_EXPECT_NEAR(velocity[0], -1.0, 2.0e-12);
  AMPS_EXPECT_NEAR(velocity[1], 0.0, 2.0e-12);
  AMPS_EXPECT_NEAR(velocity[2], 0.0, 2.0e-12);
  AMPS_EXPECT_EQ(destination_head, particle);
}

AMPS_TEST_CASE("PHYS-3D-COLLISION-RELAXATION-001",
               "populated 3-D cell conserves and equilibrates through NTC/LB collisions") {
  const double position[3] = {0.60, 0.60, 0.60};  // m, outside test sphere
  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node = MoverNode(position);
  int i, j, k;
  const int local_cell = PopulationCell(node, position, i, j, k);
  long int& first = node->block->FirstCellParticleTable[local_cell];
  ParticleBufferScope particles(4096);
  EmptyListHeadScope cell_list_guard(first);
  const int population = 160;
  const double beam_speed = 800.0;  // m/s
  for (int index = 0; index < population; ++index) {
    const long int particle = PIC::ParticleBuffer::GetNewParticle(first);
    const double velocity[3] = {
        (index % 2 == 0) ? beam_speed : -beam_speed, 0.0, 0.0};
    PIC::ParticleBuffer::SetI(_H2O_SPEC_, particle);
    PIC::ParticleBuffer::SetX(const_cast<double*>(position), particle);
    PIC::ParticleBuffer::SetV(const_cast<double*>(velocity), particle);
    PIC::IDF::LB::SetRotE(
        0.0, PIC::ParticleBuffer::GetParticleDataPointer(particle));
  }

  PIC::Mesh::cDataCenterNode* cell =
      node->block->GetCenterNode(_getCenterNodeLocalNumber(i, j, k));
  AMPS_EXPECT_TRUE(cell != NULL && cell->Measure > 0.0);
  MeshCellStateScope cell_state(node->block, cell);
  const double cross_section = 3.4e-19;  // m^2, configured H2O-H2O value
  const double dt = 1.0;                 // s
  const double target_candidates = 30.0;
  // NTC's identical-particle candidate expectation is
  // 0.5*N^2*W*sigma*g_max*dt/V. Choose W from that equation so each call is
  // statistically useful without triggering the collision-limiting branch.
  const double weight = target_candidates * cell->Measure /
      (0.5 * population * population * cross_section *
       (2.0 * beam_speed) * dt);
  node->block->SetLocalTimeStep(dt, _H2O_SPEC_);
  node->block->SetLocalParticleWeight(weight, _H2O_SPEC_);

  const PopulationTotals initial = MeasurePopulation(first);
  rnd_seed(5701);
  for (int iteration = 0; iteration < 500; ++iteration)
    PIC::MolecularCollisions::ParticleCollisionModel::ModelCellCollisions_ntc(
        i, j, k, node);
  const PopulationTotals final = MeasurePopulation(first);

  const double initial_total = initial.translational + initial.rotational;
  const double final_total = final.translational + final.rotational;
  AMPS_EXPECT_TRUE(initial_total > 0.0);
  AMPS_EXPECT_TRUE(std::fabs(final_total - initial_total) / initial_total <
                   5.0e-13);
  const double momentum_scale =
      population * PIC::MolecularData::GetMass(_H2O_SPEC_) * beam_speed;
  for (int dimension = 0; dimension < 3; ++dimension)
    AMPS_EXPECT_TRUE(std::fabs(final.momentum[dimension] -
                               initial.momentum[dimension]) /
                         momentum_scale < 5.0e-13);
  AMPS_EXPECT_TRUE(final.rotational > 0.0);
  // H2O has three translational and three rotational quadratic degrees per
  // molecule. At equilibrium each reservoir therefore receives half of the
  // COM-frame energy (the three conserved COM degrees give only an O(1/N)
  // correction, and here total momentum is exactly zero).
  const double rotational_fraction = final.rotational / final_total;
  AMPS_EXPECT_TRUE(std::fabs(rotational_fraction - 0.5) < 0.10);
}

AMPS_TEST_CASE("PHYS-3D-COLLISION-RATE-001",
               "populated 3-D NTC cell reproduces the Maxwellian collision rate") {
  namespace Collision = PIC::MolecularCollisions::ParticleCollisionModel;
  const double position[3] = {0.60, 0.60, 0.60};
  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node = MoverNode(position);
  int i, j, k;
  const int local_cell = PopulationCell(node, position, i, j, k);
  long int& first = node->block->FirstCellParticleTable[local_cell];
  ParticleBufferScope particles(4096);
  EmptyListHeadScope cell_list_guard(first);
  const int per_species = 240;
  const double temperature = 500.0;  // K
  const int species[2] = {_H2O_SPEC_, _O_SPEC_};
  const double zero_drift[3] = {0.0, 0.0, 0.0};
  rnd_seed(5703);
  for (int group = 0; group < 2; ++group) {
    for (int index = 0; index < per_species; ++index) {
      const long int particle = PIC::ParticleBuffer::GetNewParticle(first);
      double velocity[3];
      PIC::Distribution::MaxwellianVelocityDistribution(
          velocity, zero_drift, temperature, species[group]);
      PIC::ParticleBuffer::SetI(species[group], particle);
      PIC::ParticleBuffer::SetX(const_cast<double*>(position), particle);
      PIC::ParticleBuffer::SetV(velocity, particle);
      if (species[group] == _H2O_SPEC_)
        PIC::IDF::LB::InitRotTemp(
            temperature, PIC::ParticleBuffer::GetParticleDataPointer(particle));
    }
  }

  PIC::Mesh::cDataCenterNode* cell =
      node->block->GetCenterNode(_getCenterNodeLocalNumber(i, j, k));
  AMPS_EXPECT_TRUE(cell != NULL && cell->Measure > 0.0);
  MeshCellStateScope cell_state(node->block, cell);
  const double mass0 = PIC::MolecularData::GetMass(species[0]);
  const double mass1 = PIC::MolecularData::GetMass(species[1]);
  const double reduced_mass = mass0 * mass1 / (mass0 + mass1);
  const double mean_relative_speed =
      std::sqrt(8.0 * Kbol * temperature / (Pi * reduced_mass));
  const double cross_section = 3.4e-19;  // m^2, configured cross pair
  const double dt = 1.0;                 // s
  const double target_candidates = 60.0;
  // A Maxwellian maximum in this bounded population is normally below 3.5
  // times its mean relative speed. This estimate selects a useful statistical
  // weight only; the production NTC kernel evaluates its own majorant.
  const double weight = target_candidates * cell->Measure /
      (per_species * per_species * cross_section *
       3.5 * mean_relative_speed * dt);
  for (int group = 0; group < 2; ++group) {
    node->block->SetLocalTimeStep(dt, species[group]);
    node->block->SetLocalParticleWeight(weight, species[group]);
  }

  char* sampling = cell->GetAssociatedDataBufferPointer() +
      PIC::Mesh::collectingCellSampleDataPointerOffset;
  const int directed_offset =
      Collision::CollsionFrequentcySampling::SamplingBufferOffset +
      sizeof(double) * Collision::CollsionFrequentcySampling::Offset(
          species[0], species[1]);
  *reinterpret_cast<double*>(sampling + directed_offset) = 0.0;
  const int iterations = 1000;
  for (int iteration = 0; iteration < iterations; ++iteration)
    Collision::ModelCellCollisions_ntc(i, j, k, node);
  const double observed_rate =
      *reinterpret_cast<double*>(sampling + directed_offset) / iterations;
  const double number_density = per_species * weight / cell->Measure;
  const double expected_rate = number_density * number_density *
                               cross_section * mean_relative_speed;
  // Finite-population temperature fluctuations, stochastic NTC acceptance,
  // and the sampled initial Maxwellians all contribute. Twelve percent is a
  // conservative gate for this fixed 240+240 particle fixture while still
  // rejecting missing 1/V, dt, weight, or relative-speed factors.
  AMPS_EXPECT_TRUE(std::fabs(observed_rate - expected_rate) / expected_rate <
                   0.12);
}

AMPS_TEST_CASE("PHYS-3D-PHOTODECAY-001",
               "allocated 3-D population follows analytic parent survival") {
  using PIC::ChemicalReactions::PhotolyticReactions::ConstantTotalLifeTime;
  using PIC::ChemicalReactions::PhotolyticReactions::PhotolyticReaction;
  const double position[3] = {0.60, 0.60, 0.60};
  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node = MoverNode(position);
  int i, j, k;
  const int local_cell = PopulationCell(node, position, i, j, k);
  long int& first = node->block->FirstCellParticleTable[local_cell];
  ParticleBufferScope particles(4096);
  EmptyListHeadScope cell_list_guard(first);
  PhotolysisLifetimeScope lifetime_scope;
  const int population = 2400;
  for (int index = 0; index < population; ++index) {
    const long int particle = PIC::ParticleBuffer::GetNewParticle(first);
    PIC::ParticleBuffer::SetI(_H2O_SPEC_, particle);
    PIC::ParticleBuffer::SetX(const_cast<double*>(position), particle);
  }
  // Exercise the configured whole-mesh production driver before applying the
  // generic parent-decay scheduler below.  The core-unit-tests configuration
  // deliberately installs PhotolyticReactionProcessor_default, whose contract
  // is traversal/reinsertion without daughter production.  It reverses each
  // intrusive list; list integrity and population must remain unchanged.
  PIC::ChemicalReactions::PhotolyticReactions::ExecutePhotochemicalModel();
  PIC::ParticleBuffer::CheckParticleList();
  int driver_population = 0;
  long int driver_previous = -1;
  for (long int particle = first; particle != -1;
       particle = PIC::ParticleBuffer::GetNext(particle)) {
    ++driver_population;
    AMPS_EXPECT_EQ(PIC::ParticleBuffer::GetPrev(particle), driver_previous);
    driver_previous = particle;
  }
  AMPS_EXPECT_EQ(driver_population, population);
  for (int species = 0; species < PIC::nTotalSpecies; ++species)
    ConstantTotalLifeTime[species] = -1.0;
  const double lifetime = 2.0;  // s
  const double interval_requested = 1.0;  // s
  ConstantTotalLifeTime[_H2O_SPEC_] = lifetime;

  int decays = 0;
  long int survivors = -1;
  rnd_seed(5702);
  for (long int particle = first; particle != -1;) {
    const long int next = PIC::ParticleBuffer::GetNext(particle);
    int species = PIC::ParticleBuffer::GetI(particle);
    double interval = interval_requested;
    const int result = PhotolyticReaction(
        PIC::ParticleBuffer::GetX(particle), particle, species, interval, node);
    if (result == _PHOTOLYTIC_REACTION_OCCURES_) {
      ++decays;
      PIC::ParticleBuffer::DeleteParticle(particle);
    } else {
      AMPS_EXPECT_EQ(result, _PHOTOLYTIC_REACTIONS_NO_TRANSPHORMATION_);
      PIC::ParticleBuffer::SetPrev(-1, particle);
      PIC::ParticleBuffer::SetNext(survivors, particle);
      if (survivors != -1) PIC::ParticleBuffer::SetPrev(particle, survivors);
      survivors = particle;
    }
    particle = next;
  }
  first = survivors;
  const double probability = 1.0 - std::exp(-interval_requested / lifetime);
  const double sigma = std::sqrt(
      population * probability * (1.0 - probability));
  AMPS_EXPECT_TRUE(std::fabs(decays - population * probability) < 6.0 * sigma);
  int survivor_count = 0;
  long int previous = -1;
  for (long int particle = first; particle != -1;
       particle = PIC::ParticleBuffer::GetNext(particle)) {
    ++survivor_count;
    AMPS_EXPECT_EQ(PIC::ParticleBuffer::GetI(particle), _H2O_SPEC_);
    AMPS_EXPECT_EQ(PIC::ParticleBuffer::GetPrev(particle), previous);
    previous = particle;
  }
  AMPS_EXPECT_EQ(survivor_count, population - decays);
  // CheckParticleList independently recounts every AMR cell and rejects an
  // unallocated entry, a cycle/overrun, a nonempty mover list, or disagreement
  // with NAllPart.  The explicit walk above additionally validates Prev links,
  // which the production checker intentionally does not inspect.
  PIC::ParticleBuffer::CheckParticleList();
}
