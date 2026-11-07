#include "pic.h"
#include "src/test/native/test_harness.h"
#include "core_test_mesh.h"
#include "core_test_state.h"

#include <algorithm>
#include <cmath>
#include <vector>

/* Stage-T5 collision and internal-energy fixtures.
 *
 * Every compiled collision selector is driven by the same populated-cell
 * rate oracle n_i*n_j*<sigma*g>, with the documented directed-event
 * convention for equal species. The fixtures separately exercise zero-rate,
 * empty/singleton, collision-limiting, unequal global W/dt, and supported
 * individual-record-weight domains while checking SI momentum [kg m/s] and
 * total energy [J]. VHS/VSS expectations are reconstructed from literal deck
 * parameters and independent quadrature/angular moments. Larsen--Borgnakke
 * cases use conserved microcanonical energy partitions and explicit
 * monatomic/linear/nonlinear mode counts. Unsupported Yinsi individual
 * weights and VT/VV relaxation remain typed capabilities; the tests do not
 * manufacture a physical law merely to obtain coverage.
 */

namespace AMPS {
namespace CoreTests {

// Archive anchor: the native application is linked from a static archive, so
// main.cpp references one symbol from each test translation unit to guarantee
// that static AMPS_TEST_CASE registrations are retained by the linker.
void LinkT5Tests() {}

}  // namespace CoreTests
}  // namespace AMPS

namespace {

using AMPS::CoreTests::Fixture::EmptyListHeadScope;
using AMPS::CoreTests::Fixture::MeshCellStateScope;
using AMPS::CoreTests::Fixture::ParticleBufferScope;

typedef void (*CellCollisionFunction)(
    int, int, int, cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>*);

struct ConservedState {
  double momentum[3];       // kg m s^-1 per equally weighted model population
  double total_energy;      // J: translation + rotation + vibration
};

ConservedState MeasureConservedState(long int first) {
  ConservedState result = {{0.0, 0.0, 0.0}, 0.0};
  for (long int particle = first; particle != -1;
       particle = PIC::ParticleBuffer::GetNext(particle)) {
    const int species = PIC::ParticleBuffer::GetI(particle);
    const double mass = PIC::MolecularData::GetMass(species);
    const double* velocity = PIC::ParticleBuffer::GetV(particle);
    for (int dimension = 0; dimension < 3; ++dimension) {
      result.momentum[dimension] += mass * velocity[dimension];
      result.total_energy +=
          0.5 * mass * velocity[dimension] * velocity[dimension];
    }
    // The T5 deck enables classical LB storage.  Collision algorithms may
    // move energy between translation and rotation, so a translational-only
    // conservation oracle would falsely report a physically valid RT event.
    PIC::ParticleBuffer::byte* record =
        PIC::ParticleBuffer::GetParticleDataPointer(particle);
    result.total_energy += PIC::IDF::LB::GetRotE(record);
    result.total_energy += PIC::IDF::LB::GetVibE(-1, record);
  }
  return result;
}

void RunCollisionRateFixture(CellCollisionFunction collide,
                             unsigned int seed, bool same_species) {
  namespace Collision = PIC::MolecularCollisions::ParticleCollisionModel;
  const double position[3] = {0.60, 0.60, 0.60};  // m; outside T4 sphere
  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node =
      AMPS::CoreTests::MoverFixtureNode(position);
  int i, j, k;
  const int local = PIC::Mesh::mesh->FindCellIndex(
      const_cast<double*>(position), i, j, k, node, false);
  AMPS_EXPECT_TRUE(local >= 0);
  const int cell_index =
      i + _BLOCK_CELLS_X_ * (j + _BLOCK_CELLS_Y_ * k);
  long int& first = node->block->FirstCellParticleTable[cell_index];

  // Equal statistical weights deliberately isolate each algorithm's physical
  // pair-rate convention.  Unequal individual weights are a separate profile
  // because their conservation measure must include the record correction.
  ParticleBufferScope particles(4096);
  EmptyListHeadScope list_scope(first);
  const int count0 = same_species ? 240 : 160;
  const int count1 = same_species ? 0 : 160;
  const int species0 = _H2O_SPEC_;
  const int species1 = same_species ? _H2O_SPEC_ : _O_SPEC_;
  const double temperature = 500.0;  // K
  const double bulk[3] = {0.0, 0.0, 0.0};
  rnd_seed(seed);
  for (int group = 0; group < (same_species ? 1 : 2); ++group) {
    const int species = group == 0 ? species0 : species1;
    const int count = group == 0 ? count0 : count1;
    for (int index = 0; index < count; ++index) {
      const long int particle = PIC::ParticleBuffer::GetNewParticle(first);
      double velocity[3];
      PIC::Distribution::MaxwellianVelocityDistribution(
          velocity, bulk, temperature, species);
      PIC::ParticleBuffer::SetI(species, particle);
      PIC::ParticleBuffer::SetX(const_cast<double*>(position), particle);
      PIC::ParticleBuffer::SetV(velocity, particle);
      // Initialize rotation at the same 500 K as translation.  Otherwise the
      // repeated RT kernel would deliberately cool translation while filling
      // a cold rotational reservoir, and a stationary-Maxwellian rate would
      // no longer be the correct oracle.  Vibration stays cold because VT is
      // explicitly unqualified/off in this profile; all reservoirs remain in
      // the exact total-energy check above.
      PIC::ParticleBuffer::byte* record =
          PIC::ParticleBuffer::GetParticleDataPointer(particle);
      PIC::IDF::LB::InitRotTemp(temperature, record);
      for (int mode = 0;
           mode < PIC::IDF::nTotalVibtationalModes[species]; ++mode)
        PIC::IDF::LB::SetVibE(0.0, mode, record);
    }
  }

  PIC::Mesh::cDataCenterNode* cell =
      node->block->GetCenterNode(_getCenterNodeLocalNumber(i, j, k));
  AMPS_EXPECT_TRUE(cell != NULL && cell->Measure > 0.0);
  MeshCellStateScope cell_state(node->block, cell);
  PIC::Mesh::flushCollectingSamplingBuffer(cell);

  const double mass0 = PIC::MolecularData::GetMass(species0);
  const double mass1 = PIC::MolecularData::GetMass(species1);
  const double reduced_mass = mass0 * mass1 / (mass0 + mass1);
  // Independent Maxwell relative-speed law:
  // <g>=sqrt(8 k T/(pi mu)), with mu=m0*m1/(m0+m1).
  const double mean_relative_speed =
      std::sqrt(8.0 * Kbol * temperature / (Pi * reduced_mass));
  const double sigma = 3.4e-19;  // m^2, generated constant pair table
  const double dt = 0.5;         // s
  const double pair_factor = same_species
      ? 0.5 * count0 * count0
      : static_cast<double>(count0) * count1;
  const double target_candidates = 24.0;
  // Candidate count scales as N_pair W sigma g_max dt/V.  Four times the
  // Maxwell mean is only a weight-selection estimate; each production kernel
  // independently constructs and applies its actual majorant.
  const double weight = target_candidates * cell->Measure /
      (pair_factor * sigma * (4.0 * mean_relative_speed) * dt);
  node->block->SetLocalParticleWeight(weight, species0);
  node->block->SetLocalTimeStep(dt, species0);
  if (!same_species) {
    node->block->SetLocalParticleWeight(weight, species1);
    node->block->SetLocalTimeStep(dt, species1);
  }

  char* sampling = cell->GetAssociatedDataBufferPointer() +
      PIC::Mesh::collectingCellSampleDataPointerOffset;
  const int offset = Collision::CollsionFrequentcySampling::SamplingBufferOffset +
      sizeof(double) * Collision::CollsionFrequentcySampling::Offset(
          species0, species1);
  *reinterpret_cast<double*>(sampling + offset) = 0.0;
  const ConservedState before = MeasureConservedState(first);
  const int iterations = 700;
  for (int iteration = 0; iteration < iterations; ++iteration)
    collide(i, j, k, node);
  const ConservedState after = MeasureConservedState(first);

  // The sampled quantity is a directed physical collision frequency per
  // volume [m^-3 s^-1].  For unlike species it is n0*n1*<sigma*g>.  AMPS'
  // documented same-species sampling convention records both participants,
  // so its directed value is n^2*<sigma*g>, not the unordered-event value
  // 0.5*n^2*<sigma*g>.
  const double density0 = count0 * weight / cell->Measure;
  const double density1 = same_species
      ? density0 : count1 * weight / cell->Measure;
  const double expected_rate =
      density0 * density1 * sigma * mean_relative_speed;
  const double observed_rate =
      *reinterpret_cast<double*>(sampling + offset) / iterations;
  AMPS_EXPECT_TRUE(observed_rate > 0.0);
  // Normalize before reporting so a failed receipt carries a dimensionless
  // convention error (for example exactly 1/2 or 2), not an unreadable
  // 10^N physical rate.  The 15% gate includes finite-population Maxwellian
  // fluctuations and stochastic majorant acceptance but rejects missing
  // density, pair-factor, time-step, or volume factors.
  AMPS_EXPECT_NEAR(observed_rate / expected_rate, 1.0, 0.15);

  const double energy_scale = std::max(before.total_energy, 1.0e-300);
  AMPS_EXPECT_TRUE(std::fabs(after.total_energy - before.total_energy) /
                       energy_scale < 2.0e-11);
  const double momentum_scale = (count0 + count1) *
      std::max(mass0, mass1) * mean_relative_speed;
  for (int dimension = 0; dimension < 3; ++dimension)
    AMPS_EXPECT_TRUE(std::fabs(after.momentum[dimension] -
                               before.momentum[dimension]) /
                         momentum_scale < 2.0e-12);
}

void RunEmptyAndSingletonFixture(CellCollisionFunction collide,
                                 unsigned int seed) {
  const double position[3] = {0.60, 0.60, 0.60};
  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node =
      AMPS::CoreTests::MoverFixtureNode(position);
  int i, j, k;
  AMPS_EXPECT_TRUE(PIC::Mesh::mesh->FindCellIndex(
                       const_cast<double*>(position), i, j, k, node, false) >= 0);
  long int& first = node->block->FirstCellParticleTable[
      i + _BLOCK_CELLS_X_ * (j + _BLOCK_CELLS_Y_ * k)];
  ParticleBufferScope particles(32);
  EmptyListHeadScope list_scope(first);
  PIC::Mesh::cDataCenterNode* cell =
      node->block->GetCenterNode(_getCenterNodeLocalNumber(i, j, k));
  MeshCellStateScope cell_state(node->block, cell);
  node->block->SetLocalParticleWeight(1.0, _O_SPEC_);
  node->block->SetLocalTimeStep(1.0, _O_SPEC_);

  // Empty and singleton cells have no unordered pair.  Every selector must
  // return without allocating/deleting a record or perturbing its velocity.
  collide(i, j, k, node);
  AMPS_EXPECT_EQ(first, -1);
  rnd_seed(seed);
  const long int particle = PIC::ParticleBuffer::GetNewParticle(first);
  PIC::ParticleBuffer::SetI(_O_SPEC_, particle);
  PIC::ParticleBuffer::SetX(const_cast<double*>(position), particle);
  double velocity[3] = {120.0, -30.0, 8.0};
  PIC::ParticleBuffer::SetV(velocity, particle);
  collide(i, j, k, node);
  AMPS_EXPECT_EQ(first, particle);
  AMPS_EXPECT_TRUE(PIC::ParticleBuffer::IsParticleAllocated(particle));
  for (int dimension = 0; dimension < 3; ++dimension)
    AMPS_EXPECT_NEAR(PIC::ParticleBuffer::GetV(particle)[dimension],
                     velocity[dimension], 0.0);
}

struct CollisionLimitingThresholdScope {
  double saved;
  explicit CollisionLimitingThresholdScope(double value)
      : saved(PIC::MolecularCollisions::ParticleCollisionModel::
                  CollisionLimitingThrehold) {
    // This dimensionless value is the maximum expected number of collision
    // trials per model particle before AMPS thins trials and compensates their
    // sampled physical weight. The guard makes selector-order tests isolated.
    PIC::MolecularCollisions::ParticleCollisionModel::
        CollisionLimitingThrehold = value;
  }
  ~CollisionLimitingThresholdScope() {
    PIC::MolecularCollisions::ParticleCollisionModel::
        CollisionLimitingThrehold = saved;
  }
};

void RunZeroRelativeSpeedFixture(CellCollisionFunction collide,
                                 unsigned int seed) {
  namespace Collision = PIC::MolecularCollisions::ParticleCollisionModel;
  const double position[3] = {0.60, 0.60, 0.60};  // m
  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node =
      AMPS::CoreTests::MoverFixtureNode(position);
  int i, j, k;
  AMPS_EXPECT_TRUE(PIC::Mesh::mesh->FindCellIndex(
      const_cast<double*>(position), i, j, k, node, false) >= 0);
  long int& first = node->block->FirstCellParticleTable[
      i + _BLOCK_CELLS_X_ * (j + _BLOCK_CELLS_Y_ * k)];
  ParticleBufferScope particles(16);
  EmptyListHeadScope list_scope(first);
  PIC::Mesh::cDataCenterNode* cell =
      node->block->GetCenterNode(_getCenterNodeLocalNumber(i, j, k));
  MeshCellStateScope cell_state(node->block, cell);
  PIC::Mesh::flushCollectingSamplingBuffer(cell);
  node->block->SetLocalParticleWeight(2.0, _O_SPEC_);
  node->block->SetLocalTimeStep(0.5, _O_SPEC_);

  // Two identical monatomic records have g=0 and therefore sigma*g=0 for the
  // configured constant cross section. A selector must return a finite
  // no-collision result even when its estimated majorant is exactly zero.
  const double common_velocity[3] = {17.0, -9.0, 4.0};  // m/s
  for (int index = 0; index < 2; ++index) {
    const long int particle = PIC::ParticleBuffer::GetNewParticle(first);
    PIC::ParticleBuffer::SetI(_O_SPEC_, particle);
    PIC::ParticleBuffer::SetX(const_cast<double*>(position), particle);
    PIC::ParticleBuffer::SetV(
        const_cast<double*>(common_velocity), particle);
  }
  char* sampling = cell->GetAssociatedDataBufferPointer() +
      PIC::Mesh::collectingCellSampleDataPointerOffset;
  const int offset = Collision::CollsionFrequentcySampling::SamplingBufferOffset +
      sizeof(double) * Collision::CollsionFrequentcySampling::Offset(
          _O_SPEC_, _O_SPEC_);
  *reinterpret_cast<double*>(sampling + offset) = 0.0;
  rnd_seed(seed);
  collide(i, j, k, node);

  int count = 0;
  for (long int particle = first; particle != -1;
       particle = PIC::ParticleBuffer::GetNext(particle)) {
    ++count;
    const double* velocity = PIC::ParticleBuffer::GetV(particle);
    for (int dimension = 0; dimension < 3; ++dimension)
      AMPS_EXPECT_NEAR(velocity[dimension], common_velocity[dimension], 0.0);
  }
  AMPS_EXPECT_EQ(count, 2);
  AMPS_EXPECT_NEAR(*reinterpret_cast<double*>(sampling + offset), 0.0, 0.0);
}

void RunCollisionLimitingFixture(CellCollisionFunction collide,
                                 unsigned int seed) {
  namespace Collision = PIC::MolecularCollisions::ParticleCollisionModel;
  const double position[3] = {0.60, 0.60, 0.60};  // m
  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node =
      AMPS::CoreTests::MoverFixtureNode(position);
  int i, j, k;
  AMPS_EXPECT_TRUE(PIC::Mesh::mesh->FindCellIndex(
      const_cast<double*>(position), i, j, k, node, false) >= 0);
  long int& first = node->block->FirstCellParticleTable[
      i + _BLOCK_CELLS_X_ * (j + _BLOCK_CELLS_Y_ * k)];
  ParticleBufferScope particles(256);
  EmptyListHeadScope list_scope(first);
  PIC::Mesh::cDataCenterNode* cell =
      node->block->GetCenterNode(_getCenterNodeLocalNumber(i, j, k));
  MeshCellStateScope cell_state(node->block, cell);
  PIC::Mesh::flushCollectingSamplingBuffer(cell);

  const int count = 40;
  const double speed = 600.0;       // m/s; opposing beams have g=1200 m/s
  const double sigma = 3.4e-19;     // m^2; independent deck value for O/O
  const double dt = 0.4;            // s
  const double target_trials = 20000.0;
  // For one species, the NTC prospective unordered trial count is
  // 0.5*N^2*W*sigma*g_max*dt/V. Choosing W from that equation forces the
  // candidate count above N^2 and above threshold*N for every selector.
  const double weight = target_trials * 2.0 * cell->Measure /
      (count * count * sigma * (2.0 * speed) * dt);
  node->block->SetLocalParticleWeight(weight, _O_SPEC_);
  node->block->SetLocalTimeStep(dt, _O_SPEC_);
  std::vector<long int> records;
  for (int index = 0; index < count; ++index) {
    const long int particle = PIC::ParticleBuffer::GetNewParticle(first);
    PIC::ParticleBuffer::SetI(_O_SPEC_, particle);
    PIC::ParticleBuffer::SetX(const_cast<double*>(position), particle);
    const double velocity[3] = {
        index < count / 2 ? speed : -speed, 0.0, 0.0};
    PIC::ParticleBuffer::SetV(const_cast<double*>(velocity), particle);
    records.push_back(particle);
  }

  char* sampling = cell->GetAssociatedDataBufferPointer() +
      PIC::Mesh::collectingCellSampleDataPointerOffset;
  const int offset = Collision::CollsionFrequentcySampling::SamplingBufferOffset +
      sizeof(double) * Collision::CollsionFrequentcySampling::Offset(
          _O_SPEC_, _O_SPEC_);
  *reinterpret_cast<double*>(sampling + offset) = 0.0;
  CollisionLimitingThresholdScope threshold(0.25);
  rnd_seed(seed);
  const int iterations = 300;
  for (int iteration = 0; iteration < iterations; ++iteration) {
    const ConservedState before = MeasureConservedState(first);
    collide(i, j, k, node);
    const ConservedState after = MeasureConservedState(first);
    AMPS_EXPECT_TRUE(std::fabs(after.total_energy - before.total_energy) /
                         before.total_energy < 3.0e-14);
    const double momentum_scale = count *
        PIC::MolecularData::GetMass(_O_SPEC_) * speed;
    for (int dimension = 0; dimension < 3; ++dimension)
      AMPS_EXPECT_TRUE(std::fabs(after.momentum[dimension] -
                                 before.momentum[dimension]) /
                           momentum_scale < 2.0e-13);

    // Reset the two monokinetic beams after each independent collision step.
    // Their discrete relative-speed mean is exactly 600 m/s (half the ordered
    // pairs have g=0, half g=1200), so the rate oracle does not assume that a
    // repeatedly scattered finite population remains Maxwellian.
    for (int index = 0; index < count; ++index) {
      const double velocity[3] = {
          index < count / 2 ? speed : -speed, 0.0, 0.0};
      PIC::ParticleBuffer::SetV(
          const_cast<double*>(velocity), records[index]);
    }
  }

  const double density = count * weight / cell->Measure;
  const double expected_rate = density * density * sigma * speed;
  const double observed_rate =
      *reinterpret_cast<double*>(sampling + offset) / iterations;
  AMPS_EXPECT_TRUE(observed_rate > 0.0);
  // The aggressive limiter deliberately reduces some selectors to only a few
  // weighted trials per step. A 22% statistical gate accepts that variance
  // while rejecting an omitted limiting compensation factor (orders of
  // magnitude) or the unordered-event factor of two.
  AMPS_EXPECT_NEAR(observed_rate / expected_rate, 1.0, 0.22);
}

void RunUnequalWeightTimeFixture(CellCollisionFunction collide,
                                 unsigned int seed) {
  namespace Collision = PIC::MolecularCollisions::ParticleCollisionModel;
  const double position[3] = {0.60, 0.60, 0.60};  // m
  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node =
      AMPS::CoreTests::MoverFixtureNode(position);
  int i, j, k;
  AMPS_EXPECT_TRUE(PIC::Mesh::mesh->FindCellIndex(
      const_cast<double*>(position), i, j, k, node, false) >= 0);
  long int& first = node->block->FirstCellParticleTable[
      i + _BLOCK_CELLS_X_ * (j + _BLOCK_CELLS_Y_ * k)];
  ParticleBufferScope particles(1024);
  EmptyListHeadScope list_scope(first);
  PIC::Mesh::cDataCenterNode* cell =
      node->block->GetCenterNode(_getCenterNodeLocalNumber(i, j, k));
  MeshCellStateScope cell_state(node->block, cell);
  PIC::Mesh::flushCollectingSamplingBuffer(cell);

  const int count0 = 160, count1 = 160;
  const int species0 = _H2O_SPEC_, species1 = _O_SPEC_;
  const double temperature = 500.0;  // K
  const double mass0 = PIC::MolecularData::GetMass(species0);
  const double mass1 = PIC::MolecularData::GetMass(species1);
  const double reduced_mass = mass0 * mass1 / (mass0 + mass1);
  const double mean_relative_speed =
      std::sqrt(8.0 * Kbol * temperature / (Pi * reduced_mass));
  const double sigma = 3.4e-19;  // m^2
  const double dt0 = 0.25, dt1 = 0.70;  // s
  // The advanced selectors form candidates over max(dt) and min(weight), then
  // independently accept updates with probabilities proportional to dt/W.
  // W1=3*W0 and unequal dt exercise all four boolean update outcomes without
  // changing the physical directed rate n0*n1*<sigma*g>.
  const double target_candidates = 36.0;
  const double weight0 = target_candidates * cell->Measure /
      (3.0 * count0 * count1 * sigma *
       (4.0 * mean_relative_speed) * dt1);
  const double weight1 = 3.0 * weight0;
  node->block->SetLocalParticleWeight(weight0, species0);
  node->block->SetLocalParticleWeight(weight1, species1);
  node->block->SetLocalTimeStep(dt0, species0);
  node->block->SetLocalTimeStep(dt1, species1);

  const double bulk[3] = {0.0, 0.0, 0.0};
  rnd_seed(seed);
  for (int group = 0; group < 2; ++group) {
    const int species = group == 0 ? species0 : species1;
    const int count_group = group == 0 ? count0 : count1;
    for (int index = 0; index < count_group; ++index) {
      const long int particle = PIC::ParticleBuffer::GetNewParticle(first);
      double velocity[3];
      PIC::Distribution::MaxwellianVelocityDistribution(
          velocity, bulk, temperature, species);
      PIC::ParticleBuffer::SetI(species, particle);
      PIC::ParticleBuffer::SetX(const_cast<double*>(position), particle);
      PIC::ParticleBuffer::SetV(velocity, particle);
      PIC::ParticleBuffer::byte* record =
          PIC::ParticleBuffer::GetParticleDataPointer(particle);
      PIC::IDF::LB::InitRotTemp(temperature, record);
    }
  }

  char* sampling = cell->GetAssociatedDataBufferPointer() +
      PIC::Mesh::collectingCellSampleDataPointerOffset;
  const int offset01 = Collision::CollsionFrequentcySampling::SamplingBufferOffset +
      sizeof(double) * Collision::CollsionFrequentcySampling::Offset(
          species0, species1);
  const int offset10 = Collision::CollsionFrequentcySampling::SamplingBufferOffset +
      sizeof(double) * Collision::CollsionFrequentcySampling::Offset(
          species1, species0);
  *reinterpret_cast<double*>(sampling + offset01) = 0.0;
  *reinterpret_cast<double*>(sampling + offset10) = 0.0;
  const int iterations = 900;
  for (int iteration = 0; iteration < iterations; ++iteration)
    collide(i, j, k, node);

  const double density0 = count0 * weight0 / cell->Measure;
  const double density1 = count1 * weight1 / cell->Measure;
  const double expected_rate =
      density0 * density1 * sigma * mean_relative_speed;
  const double observed01 =
      *reinterpret_cast<double*>(sampling + offset01) / iterations;
  const double observed10 =
      *reinterpret_cast<double*>(sampling + offset10) / iterations;
  AMPS_EXPECT_NEAR(observed01 / expected_rate, 1.0, 0.18);
  AMPS_EXPECT_NEAR(observed10 / expected_rate, 1.0, 0.18);
  AMPS_EXPECT_NEAR(observed01 / observed10, 1.0, 0.16);

  int final_count = 0;
  for (long int particle = first; particle != -1;
       particle = PIC::ParticleBuffer::GetNext(particle)) {
    ++final_count;
    const double* velocity = PIC::ParticleBuffer::GetV(particle);
    AMPS_EXPECT_TRUE(std::isfinite(velocity[0]) &&
                     std::isfinite(velocity[1]) &&
                     std::isfinite(velocity[2]));
  }
  AMPS_EXPECT_EQ(final_count, count0 + count1);
}

void RunIndividualWeightFixture(CellCollisionFunction collide,
                                unsigned int seed) {
#if _INDIVIDUAL_PARTICLE_WEIGHT_MODE_ == _INDIVIDUAL_PARTICLE_WEIGHT_ON_
  namespace Collision = PIC::MolecularCollisions::ParticleCollisionModel;
  const double position[3] = {0.60, 0.60, 0.60};  // m
  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node =
      AMPS::CoreTests::MoverFixtureNode(position);
  int i, j, k;
  AMPS_EXPECT_TRUE(PIC::Mesh::mesh->FindCellIndex(
      const_cast<double*>(position), i, j, k, node, false) >= 0);
  long int& first = node->block->FirstCellParticleTable[
      i + _BLOCK_CELLS_X_ * (j + _BLOCK_CELLS_Y_ * k)];
  ParticleBufferScope particles(1024);
  EmptyListHeadScope list_scope(first);
  PIC::Mesh::cDataCenterNode* cell =
      node->block->GetCenterNode(_getCenterNodeLocalNumber(i, j, k));
  MeshCellStateScope cell_state(node->block, cell);
  PIC::Mesh::flushCollectingSamplingBuffer(cell);

  const int count0 = 160, count1 = 160;
  const int species0 = _H2O_SPEC_, species1 = _O_SPEC_;
  const double temperature = 500.0;  // K
  const double mass0 = PIC::MolecularData::GetMass(species0);
  const double mass1 = PIC::MolecularData::GetMass(species1);
  const double reduced_mass = mass0 * mass1 / (mass0 + mass1);
  const double mean_relative_speed =
      std::sqrt(8.0 * Kbol * temperature / (Pi * reduced_mass));
  const double sigma = 3.4e-19;  // m^2
  const double dt = 0.5;         // s
  const double target_candidates = 36.0;
  // Corrections alternate around unity, so each species has total physical
  // weight N*W while its minimum record is 0.5W (H2O) or 0.25W (O). The
  // advanced selectors build candidates using that minimum and compensate
  // each real-record update by 1/correction. Choosing W from the conservative
  // 0.25W majorant estimate keeps the fixture stochastic but inexpensive.
  const double base_weight = target_candidates * cell->Measure /
      (4.0 * count0 * count1 * sigma *
       (4.0 * mean_relative_speed) * dt);
  node->block->SetLocalParticleWeight(base_weight, species0);
  node->block->SetLocalParticleWeight(base_weight, species1);
  node->block->SetLocalTimeStep(dt, species0);
  node->block->SetLocalTimeStep(dt, species1);

  const double bulk[3] = {0.0, 0.0, 0.0};
  double correction_sum[2] = {0.0, 0.0};
  rnd_seed(seed);
  for (int group = 0; group < 2; ++group) {
    const int species = group == 0 ? species0 : species1;
    const int count_group = group == 0 ? count0 : count1;
    for (int index = 0; index < count_group; ++index) {
      const long int particle = PIC::ParticleBuffer::GetNewParticle(first);
      double velocity[3];
      PIC::Distribution::MaxwellianVelocityDistribution(
          velocity, bulk, temperature, species);
      PIC::ParticleBuffer::SetI(species, particle);
      PIC::ParticleBuffer::SetX(const_cast<double*>(position), particle);
      PIC::ParticleBuffer::SetV(velocity, particle);
      const double correction = group == 0
          ? (index % 2 == 0 ? 0.5 : 1.5)
          : (index % 2 == 0 ? 0.25 : 1.75);
      PIC::ParticleBuffer::SetIndividualStatWeightCorrection(
          correction, particle);
      correction_sum[group] += correction;
      PIC::ParticleBuffer::byte* record =
          PIC::ParticleBuffer::GetParticleDataPointer(particle);
      PIC::IDF::LB::InitRotTemp(temperature, record);
    }
  }
  AMPS_EXPECT_NEAR(correction_sum[0], count0, 0.0);
  AMPS_EXPECT_NEAR(correction_sum[1], count1, 0.0);

  char* sampling = cell->GetAssociatedDataBufferPointer() +
      PIC::Mesh::collectingCellSampleDataPointerOffset;
  const int offset01 = Collision::CollsionFrequentcySampling::SamplingBufferOffset +
      sizeof(double) * Collision::CollsionFrequentcySampling::Offset(
          species0, species1);
  const int offset10 = Collision::CollsionFrequentcySampling::SamplingBufferOffset +
      sizeof(double) * Collision::CollsionFrequentcySampling::Offset(
          species1, species0);
  *reinterpret_cast<double*>(sampling + offset01) = 0.0;
  *reinterpret_cast<double*>(sampling + offset10) = 0.0;
  const int iterations = 900;
  for (int iteration = 0; iteration < iterations; ++iteration)
    collide(i, j, k, node);

  const double density0 = correction_sum[0] * base_weight / cell->Measure;
  const double density1 = correction_sum[1] * base_weight / cell->Measure;
  const double expected_rate =
      density0 * density1 * sigma * mean_relative_speed;
  const double observed01 =
      *reinterpret_cast<double*>(sampling + offset01) / iterations;
  const double observed10 =
      *reinterpret_cast<double*>(sampling + offset10) / iterations;
  AMPS_EXPECT_NEAR(observed01 / expected_rate, 1.0, 0.20);
  AMPS_EXPECT_NEAR(observed10 / expected_rate, 1.0, 0.20);
  AMPS_EXPECT_NEAR(observed01 / observed10, 1.0, 0.18);

  int final_count = 0;
  for (long int particle = first; particle != -1;
       particle = PIC::ParticleBuffer::GetNext(particle)) {
    ++final_count;
    const double* velocity = PIC::ParticleBuffer::GetV(particle);
    AMPS_EXPECT_TRUE(std::isfinite(velocity[0]) &&
                     std::isfinite(velocity[1]) &&
                     std::isfinite(velocity[2]));
  }
  AMPS_EXPECT_EQ(final_count, count0 + count1);
#else
  (void)collide;
  (void)seed;
  AMPS_EXPECT_TRUE(false);
#endif
}

double ExpectedVariableSphereCrossSection(
    double relative_speed, int species0, int species1) {
  // These are independent fixtures matching the explicit SI values in the
  // VHS/VSS input decks.  Do not read the generated production tables here:
  // doing so would allow a parser error to reproduce itself in the oracle.
  const double diameter0 = species0 == _H2O_SPEC_ ? 3.2e-10 : 2.8e-10;
  const double diameter1 = species1 == _H2O_SPEC_ ? 3.2e-10 : 2.8e-10;
  const double reference_speed0 = species0 == _H2O_SPEC_ ? 500.0 : 700.0;
  const double reference_speed1 = species1 == _H2O_SPEC_ ? 500.0 : 700.0;
  const double power0 = species0 == _H2O_SPEC_ ? 0.25 : 0.15;
  const double power1 = species1 == _H2O_SPEC_ ? 0.25 : 0.15;
  const double effective0 =
      diameter0 * std::pow(reference_speed0 / relative_speed, power0);
  const double effective1 =
      diameter1 * std::pow(reference_speed1 / relative_speed, power1);
  const double effective = species0 == species1
      ? effective0 : 0.5 * (effective0 + effective1);
  return Pi * effective * effective;
}

double MaxwellRelativeSpeedDensity(double speed, double reduced_mass,
                                   double temperature) {
  // Probability density [s/m] for the magnitude of the difference of two
  // independent 3-D Maxwell velocities at the same T:
  // 4*pi*(mu/(2*pi*kT))^(3/2)*g^2*exp[-mu*g^2/(2*kT)].
  const double coefficient = 4.0 * Pi * std::pow(
      reduced_mass / (2.0 * Pi * Kbol * temperature), 1.5);
  return coefficient * speed * speed *
      std::exp(-reduced_mass * speed * speed /
               (2.0 * Kbol * temperature));
}

void ValidateVariableSphereCrossSection(bool expect_vss) {
  namespace Models = PIC::MolecularData::MolecularModels;
  const int expected_model = expect_vss ? _VSS_MOL_MODEL_ : _VHS_MOL_MODEL_;
  AMPS_EXPECT_EQ(Models::CompiledCollisionCrossSectionModel, expected_model);

  const double speeds[3] = {125.0, 500.0, 2000.0};  // m/s
  const int pairs[2][2] = {
      {_H2O_SPEC_, _H2O_SPEC_}, {_H2O_SPEC_, _O_SPEC_}};
  for (int pair = 0; pair < 2; ++pair) {
    for (int sample = 0; sample < 3; ++sample) {
      double velocity0[3] = {0.0, 0.0, 0.0};
      double velocity1[3] = {speeds[sample], 0.0, 0.0};
      const double observed = Models::GetTotalCrossSection(
          pairs[pair][0], velocity0, pairs[pair][1], velocity1);
      const double expected = ExpectedVariableSphereCrossSection(
          speeds[sample], pairs[pair][0], pairs[pair][1]);
      AMPS_EXPECT_TRUE(std::fabs(observed - expected) / expected < 2.0e-14);
    }
  }

  // Independently integrate <sigma*g> for the unequal pair.  A midpoint rule
  // avoids evaluating the integrable VHS singularity at g=0; twelve thermal
  // speeds makes the omitted Maxwell tail negligible in double precision.
  const double temperature = 650.0;  // K
  const double mass0 = PIC::MolecularData::GetMass(_H2O_SPEC_);
  const double mass1 = PIC::MolecularData::GetMass(_O_SPEC_);
  const double reduced_mass = mass0 * mass1 / (mass0 + mass1);
  const double thermal_speed =
      std::sqrt(2.0 * Kbol * temperature / reduced_mass);
  const int intervals = 200000;
  const double maximum_speed = 12.0 * thermal_speed;
  const double delta = maximum_speed / intervals;
  double expected_integral = 0.0, observed_integral = 0.0;
  for (int interval = 0; interval < intervals; ++interval) {
    const double speed = (interval + 0.5) * delta;
    const double density = MaxwellRelativeSpeedDensity(
        speed, reduced_mass, temperature);
    const double expected_sigma = ExpectedVariableSphereCrossSection(
        speed, _H2O_SPEC_, _O_SPEC_);
    double velocity0[3] = {0.0, 0.0, 0.0};
    double velocity1[3] = {speed, 0.0, 0.0};
    const double observed_sigma = Models::GetTotalCrossSection(
        _H2O_SPEC_, velocity0, _O_SPEC_, velocity1);
    expected_integral += expected_sigma * speed * density * delta;
    observed_integral += observed_sigma * speed * density * delta;
  }
  AMPS_EXPECT_TRUE(expected_integral > 0.0);
  AMPS_EXPECT_TRUE(std::fabs(observed_integral - expected_integral) /
                       expected_integral < 2.0e-12);
}

void ValidateVssAngularMoments() {
  namespace Scattering = PIC::MolecularCollisions::VelocityScattering;
  // The deck supplies alpha_H2O=1.70 and alpha_O=2.10.  The production
  // binary-pair closure declares their arithmetic mean, alpha_01=1.90.  Keep
  // these literal fixture values independent of the generated table so an
  // omitted/misordered CHI entry cannot reproduce itself in the oracle.
  const double alpha = 0.5 * (1.70 + 2.10);
  const double expected_first = (alpha - 1.0) / (alpha + 1.0);
  const double expected_second =
      4.0 * alpha / (alpha + 2.0) -
      4.0 * alpha / (alpha + 1.0) + 1.0;
  const int samples = 200000;
  double first_sum = 0.0, second_sum = 0.0;
  const double mass0 = PIC::MolecularData::GetMass(_H2O_SPEC_);
  const double mass1 = PIC::MolecularData::GetMass(_O_SPEC_);
  const double incoming_speed = 900.0;  // m/s

  rnd_seed(5681);
  for (int sample = 0; sample < samples; ++sample) {
    // In the center-of-mass frame, the incoming relative velocity points in
    // +x. The cosine of the physical deflection angle is therefore simply
    // vrel'_x/|vrel'|. No production angular sampler is used in the oracle:
    // for U~Uniform(0,1), mu=2 U^(1/alpha)-1 gives
    // E[mu]=(alpha-1)/(alpha+1) and the second moment above.
    double vcm[3] = {0.0, 0.0, 0.0};
    double relative[3] = {incoming_speed, 0.0, 0.0};
    double velocity0[3] = {
        -mass1 / (mass0 + mass1) * incoming_speed, 0.0, 0.0};
    double velocity1[3] = {
        mass0 / (mass0 + mass1) * incoming_speed, 0.0, 0.0};
    AMPS_EXPECT_TRUE(Scattering::Scatter(
        _H2O_SPEC_, velocity0, mass0, _O_SPEC_, velocity1, mass1,
        vcm, relative));
    double outgoing[3];
    for (int dimension = 0; dimension < 3; ++dimension)
      outgoing[dimension] = velocity1[dimension] - velocity0[dimension];
    const double outgoing_speed = std::sqrt(
        outgoing[0] * outgoing[0] + outgoing[1] * outgoing[1] +
        outgoing[2] * outgoing[2]);
    AMPS_EXPECT_TRUE(std::fabs(outgoing_speed - incoming_speed) /
                         incoming_speed < 3.0e-14);
    const double cosine = outgoing[0] / outgoing_speed;
    first_sum += cosine;
    second_sum += cosine * cosine;
  }

  // Sampling errors scale as N^-1/2. These absolute gates are wider than six
  // standard errors for the configured alpha but much smaller than the
  // difference between isotropic scattering (<mu>=0) and VSS alpha=1.9.
  AMPS_EXPECT_NEAR(first_sum / samples, expected_first, 0.006);
  AMPS_EXPECT_NEAR(second_sum / samples, expected_second, 0.006);
}

void ValidateRotationalRelaxationProbability(double expected_collision_number,
                                              unsigned int seed) {
  AMPS_EXPECT_NEAR(PIC::IDF::RotationZnumber[_H2O_SPEC_],
                   expected_collision_number, 0.0);
  const int trials = 120000;
  const double initial_relative_speed = 1400.0;  // m/s
  std::vector<PIC::ParticleBuffer::byte> first(
      PIC::ParticleBuffer::ParticleDataLength, 0);
  std::vector<PIC::ParticleBuffer::byte> second(
      PIC::ParticleBuffer::ParticleDataLength, 0);
  PIC::ParticleBuffer::SetI(_H2O_SPEC_, &first[0]);
  PIC::ParticleBuffer::SetI(_H2O_SPEC_, &second[0]);
  int transitions = 0;
  rnd_seed(seed);
  for (int trial = 0; trial < trials; ++trial) {
    PIC::IDF::LB::SetRotE(0.0, &first[0]);
    PIC::IDF::LB::SetRotE(0.0, &second[0]);
    double relative_speed = initial_relative_speed;
    bool update[2] = {true, true};
    PIC::IDF::LB::RedistributeEnergy(
        &first[0], &second[0], relative_speed, update, NULL);
    // A successful cold RT transition assigns positive rotational energy and
    // therefore lowers the relative translational speed.  No floating-point
    // equality with the sampled energy is needed for event classification.
    if (relative_speed < initial_relative_speed * (1.0 - 1.0e-14))
      ++transitions;
  }
  const double probability = 1.0 / expected_collision_number;
  const double sigma =
      std::sqrt(trials * probability * (1.0 - probability));
  AMPS_EXPECT_TRUE(std::fabs(transitions - trials * probability) < 6.0 * sigma);
}

void ValidateRtHotColdDirections() {
  const int samples = 90000;
  const double total_energy = 12.0 * Kbol * 400.0;  // J per model pair
  double mean_rotational_fraction[2] = {0.0, 0.0};
  std::vector<PIC::ParticleBuffer::byte> first(
      PIC::ParticleBuffer::ParticleDataLength, 0);
  std::vector<PIC::ParticleBuffer::byte> second(
      PIC::ParticleBuffer::ParticleDataLength, 0);
  PIC::ParticleBuffer::SetI(_H2O_SPEC_, &first[0]);
  PIC::ParticleBuffer::SetI(_H2O_SPEC_, &second[0]);

  for (int direction = 0; direction < 2; ++direction) {
    rnd_seed(5720 + direction);
    double rotational_sum = 0.0;
    for (int sample = 0; sample < samples; ++sample) {
      // direction=0 starts with hot translation and cold rotation;
      // direction=1 starts with zero relative translation and all energy in
      // rotation. Larsen--Borgnakke must forget that preparation and sample
      // the same microcanonical equilibrium partition in both directions.
      const double initial_rotation =
          direction == 0 ? 0.0 : 0.5 * total_energy;
      PIC::IDF::LB::SetRotE(initial_rotation, &first[0]);
      PIC::IDF::LB::SetRotE(initial_rotation, &second[0]);
      double collision_energy = direction == 0 ? total_energy : 0.0;
      bool update[2] = {true, true};
      bool redistributed = false;
      PIC::IDF::LB::CalculateRTTransition(
          &first[0], &second[0], collision_energy, update, redistributed, 0.5);
      const double rotation = PIC::IDF::LB::GetRotE(&first[0]) +
          PIC::IDF::LB::GetRotE(&second[0]);
      AMPS_EXPECT_TRUE(redistributed && collision_energy >= 0.0 &&
                       rotation >= 0.0);
      AMPS_EXPECT_NEAR(collision_energy + rotation, total_energy,
                       2.0e-14 * total_energy);
      rotational_sum += rotation;
    }
    mean_rotational_fraction[direction] =
        rotational_sum / (samples * total_energy);
  }

  // For omega=0.5 and two nonlinear H2O rotors, the production beta shapes
  // are translational 2 and rotational 3, hence <E_rot/E_total>=3/5.
  AMPS_EXPECT_NEAR(mean_rotational_fraction[0], 3.0 / 5.0, 0.006);
  AMPS_EXPECT_NEAR(mean_rotational_fraction[1], 3.0 / 5.0, 0.006);
  AMPS_EXPECT_NEAR(mean_rotational_fraction[0],
                   mean_rotational_fraction[1], 0.006);
}

void ValidateRtUpdateFlags() {
  const double energy_unit = Kbol * 300.0;  // J
  std::vector<PIC::ParticleBuffer::byte> first(
      PIC::ParticleBuffer::ParticleDataLength, 0);
  std::vector<PIC::ParticleBuffer::byte> second(
      PIC::ParticleBuffer::ParticleDataLength, 0);
  PIC::ParticleBuffer::SetI(_H2O_SPEC_, &first[0]);
  PIC::ParticleBuffer::SetI(_H2O_SPEC_, &second[0]);

  for (int mask = 0; mask < 4; ++mask) {
    const double initial0 = 2.0 * energy_unit;
    const double initial1 = 4.0 * energy_unit;
    PIC::IDF::LB::SetRotE(initial0, &first[0]);
    PIC::IDF::LB::SetRotE(initial1, &second[0]);
    double collision_energy = 7.0 * energy_unit;
    bool update[2] = {(mask & 1) != 0, (mask & 2) != 0};
    bool redistributed = false;
    rnd_seed(5730 + mask);
    PIC::IDF::LB::CalculateRTTransition(
        &first[0], &second[0], collision_energy, update, redistributed, 0.5);
    AMPS_EXPECT_TRUE(redistributed && collision_energy >= 0.0);

    // Variable-weight/time-step collision selectors draw the two update flags
    // independently. A false flag represents a virtual collision partner and
    // must leave that record's internal state byte-for-byte unchanged; a true
    // flag receives the sampled nonnegative reservoir. Per-event conservation
    // is required only when both real records update. Asymmetric events
    // conserve the weighted ensemble in expectation through selector-specific
    // update probabilities and are tested by the directed-rate fixture.
    const double rotation0 = PIC::IDF::LB::GetRotE(&first[0]);
    const double rotation1 = PIC::IDF::LB::GetRotE(&second[0]);
    if (update[0])
      AMPS_EXPECT_TRUE(rotation0 >= 0.0 && std::isfinite(rotation0));
    else
      AMPS_EXPECT_NEAR(rotation0, initial0, 0.0);
    if (update[1])
      AMPS_EXPECT_TRUE(rotation1 >= 0.0 && std::isfinite(rotation1));
    else
      AMPS_EXPECT_NEAR(rotation1, initial1, 0.0);
    if (update[0] && update[1])
      AMPS_EXPECT_NEAR(collision_energy + rotation0 + rotation1,
                       13.0 * energy_unit, 3.0e-14 * energy_unit);
  }
}

void ValidateModeFamiliesAndTemperatureLimits() {
  // The generated mode table must distinguish monatomic O, linear OH/CO2,
  // and nonlinear H2O. Rotational equipartition then predicts means of
  // 0, kT, and 3kT/2 respectively, independently of the sampler internals.
  AMPS_EXPECT_EQ(PIC::IDF::nTotalRotationalModes[_O_SPEC_], 0);
  AMPS_EXPECT_EQ(PIC::IDF::nTotalRotationalModes[_OH_SPEC_], 2);
  AMPS_EXPECT_EQ(PIC::IDF::nTotalRotationalModes[_CO2_SPEC_], 2);
  AMPS_EXPECT_EQ(PIC::IDF::nTotalRotationalModes[_H2O_SPEC_], 3);
  AMPS_EXPECT_EQ(PIC::IDF::nTotalVibtationalModes[_O_SPEC_], 0);
  AMPS_EXPECT_EQ(PIC::IDF::nTotalVibtationalModes[_OH_SPEC_], 1);
  AMPS_EXPECT_EQ(PIC::IDF::nTotalVibtationalModes[_CO2_SPEC_], 4);
  AMPS_EXPECT_EQ(PIC::IDF::nTotalVibtationalModes[_H2O_SPEC_], 3);

  const double rotational_temperature = 420.0;  // K
  const int rotational_samples = 80000;
  const int rotating_species[2] = {_OH_SPEC_, _H2O_SPEC_};
  const double factors[2] = {1.0, 1.5};
  for (int family = 0; family < 2; ++family) {
    std::vector<PIC::ParticleBuffer::byte> record(
        PIC::ParticleBuffer::ParticleDataLength, 0);
    PIC::ParticleBuffer::SetI(rotating_species[family], &record[0]);
    rnd_seed(5740 + family);
    double sum = 0.0;
    for (int sample = 0; sample < rotational_samples; ++sample) {
      PIC::IDF::LB::InitRotTemp(rotational_temperature, &record[0]);
      sum += PIC::IDF::LB::GetRotE(&record[0]);
    }
    const double expected =
        factors[family] * Kbol * rotational_temperature;
    AMPS_EXPECT_NEAR((sum / rotational_samples) / expected, 1.0, 0.018);
  }
  std::vector<PIC::ParticleBuffer::byte> monatomic(
      PIC::ParticleBuffer::ParticleDataLength, 0);
  PIC::ParticleBuffer::SetI(_O_SPEC_, &monatomic[0]);
  PIC::IDF::LB::InitRotTemp(rotational_temperature, &monatomic[0]);
  AMPS_EXPECT_NEAR(PIC::IDF::LB::GetRotE(&monatomic[0]), 0.0, 0.0);

  std::vector<PIC::ParticleBuffer::byte> carbon_dioxide(
      PIC::ParticleBuffer::ParticleDataLength, 0);
  PIC::ParticleBuffer::SetI(_CO2_SPEC_, &carbon_dioxide[0]);
  PIC::IDF::LB::InitVibTemp(0.01, &carbon_dioxide[0]);  // K
  for (int mode = 0; mode < 4; ++mode)
    AMPS_EXPECT_NEAR(
        PIC::IDF::LB::GetVibE(mode, &carbon_dioxide[0]), 0.0, 0.0);

  // At high but finite temperature, each harmonic mode has mean energy
  // k*Theta/[exp(Theta/T)-1], approaching kT without adding zero-point energy.
  const double high_temperature = 50000.0;  // K
  const double theta[4] = {800.0, 900.0, 1000.0, 1100.0};  // K
  const int vibrational_samples = 80000;
  double sums[4] = {0.0, 0.0, 0.0, 0.0};
  rnd_seed(5745);
  for (int sample = 0; sample < vibrational_samples; ++sample) {
    PIC::IDF::LB::InitVibTemp(high_temperature, &carbon_dioxide[0]);
    for (int mode = 0; mode < 4; ++mode)
      sums[mode] += PIC::IDF::LB::GetVibE(mode, &carbon_dioxide[0]);
  }
  for (int mode = 0; mode < 4; ++mode) {
    const double expected = Kbol * theta[mode] /
        (std::exp(theta[mode] / high_temperature) - 1.0);
    AMPS_EXPECT_NEAR((sums[mode] / vibrational_samples) / expected,
                     1.0, 0.025);
  }
}

}  // namespace

#define AMPS_T5_COLLISION_CASE(ID, DESCRIPTION, FUNCTION, SEED)             \
  AMPS_TEST_CASE(ID, DESCRIPTION) {                                         \
    RunCollisionRateFixture(FUNCTION, SEED, false);                         \
    RunCollisionRateFixture(FUNCTION, (SEED) + 100, true);                  \
  }

AMPS_T5_COLLISION_CASE(
    "T5-COLLISION-NTC-001",
    "NTC reproduces unlike/same directed rates and conserves pair invariants",
    PIC::MolecularCollisions::ParticleCollisionModel::ModelCellCollisions_ntc,
    5501)
AMPS_T5_COLLISION_CASE(
    "T5-COLLISION-NTC-SIMPLIFIED-001",
    "simplified NTC reproduces unlike/same directed rates and invariants",
    PIC::MolecularCollisions::ParticleCollisionModel::
        ModelCellCollisions_ntc_simplified,
    5502)
AMPS_T5_COLLISION_CASE(
    "T5-COLLISION-MF-001",
    "majorant-frequency reproduces unlike/same directed rates and invariants",
    PIC::MolecularCollisions::ParticleCollisionModel::ModelCellCollisions_mf,
    5503)
AMPS_T5_COLLISION_CASE(
    "T5-COLLISION-MF-SIMPLIFIED-001",
    "simplified majorant-frequency reproduces directed rates and invariants",
    PIC::MolecularCollisions::ParticleCollisionModel::
        ModelCellCollisions_mf_simplified,
    5504)
AMPS_T5_COLLISION_CASE(
    "T5-COLLISION-MF-YINSI-001",
    "Yinsi majorant-frequency reproduces directed rates and invariants",
    PIC::MolecularCollisions::ParticleCollisionModel::ModelCellCollisions_mf_Yinsi,
    5505)

#undef AMPS_T5_COLLISION_CASE

AMPS_TEST_CASE("T5-COLLISION-EDGE-CELLS-001",
               "all collision selectors preserve empty and singleton cells") {
  CellCollisionFunction functions[] = {
      PIC::MolecularCollisions::ParticleCollisionModel::ModelCellCollisions_ntc,
      PIC::MolecularCollisions::ParticleCollisionModel::
          ModelCellCollisions_ntc_simplified,
      PIC::MolecularCollisions::ParticleCollisionModel::ModelCellCollisions_mf,
      PIC::MolecularCollisions::ParticleCollisionModel::
          ModelCellCollisions_mf_simplified,
      PIC::MolecularCollisions::ParticleCollisionModel::ModelCellCollisions_mf_Yinsi};
  for (int algorithm = 0; algorithm < 5; ++algorithm)
    RunEmptyAndSingletonFixture(functions[algorithm], 5600 + algorithm);
}

#define AMPS_T5_ZERO_CASE(ID, DESCRIPTION, FUNCTION, SEED)                 \
  AMPS_TEST_CASE(ID, DESCRIPTION) {                                        \
    RunZeroRelativeSpeedFixture(FUNCTION, SEED);                           \
  }

AMPS_T5_ZERO_CASE(
    "T5-COLLISION-ZERO-NTC-001", "NTC preserves a cold zero-rate pair",
    PIC::MolecularCollisions::ParticleCollisionModel::ModelCellCollisions_ntc,
    5610)
AMPS_T5_ZERO_CASE(
    "T5-COLLISION-ZERO-NTC-SIMPLIFIED-001",
    "simplified NTC preserves a cold zero-rate pair",
    PIC::MolecularCollisions::ParticleCollisionModel::
        ModelCellCollisions_ntc_simplified,
    5611)
AMPS_T5_ZERO_CASE(
    "T5-COLLISION-ZERO-MF-001", "majorant-frequency preserves a cold zero-rate pair",
    PIC::MolecularCollisions::ParticleCollisionModel::ModelCellCollisions_mf,
    5612)
AMPS_T5_ZERO_CASE(
    "T5-COLLISION-ZERO-MF-SIMPLIFIED-001",
    "simplified majorant-frequency preserves a cold zero-rate pair",
    PIC::MolecularCollisions::ParticleCollisionModel::
        ModelCellCollisions_mf_simplified,
    5613)
AMPS_T5_ZERO_CASE(
    "T5-COLLISION-ZERO-MF-YINSI-001",
    "Yinsi majorant-frequency preserves a cold zero-rate pair",
    PIC::MolecularCollisions::ParticleCollisionModel::ModelCellCollisions_mf_Yinsi,
    5614)

#undef AMPS_T5_ZERO_CASE

#define AMPS_T5_LIMIT_CASE(ID, DESCRIPTION, FUNCTION, SEED)                \
  AMPS_TEST_CASE(ID, DESCRIPTION) {                                        \
    RunCollisionLimitingFixture(FUNCTION, SEED);                           \
  }

AMPS_T5_LIMIT_CASE(
    "T5-COLLISION-LIMITING-NTC-001", "NTC limiting preserves rate and invariants",
    PIC::MolecularCollisions::ParticleCollisionModel::ModelCellCollisions_ntc,
    5620)
AMPS_T5_LIMIT_CASE(
    "T5-COLLISION-LIMITING-NTC-SIMPLIFIED-001",
    "simplified NTC limiting preserves rate and invariants",
    PIC::MolecularCollisions::ParticleCollisionModel::
        ModelCellCollisions_ntc_simplified,
    5621)
AMPS_T5_LIMIT_CASE(
    "T5-COLLISION-LIMITING-MF-001",
    "majorant-frequency limiting preserves rate and invariants",
    PIC::MolecularCollisions::ParticleCollisionModel::ModelCellCollisions_mf,
    5622)
AMPS_T5_LIMIT_CASE(
    "T5-COLLISION-LIMITING-MF-SIMPLIFIED-001",
    "simplified majorant-frequency limiting preserves rate and invariants",
    PIC::MolecularCollisions::ParticleCollisionModel::
        ModelCellCollisions_mf_simplified,
    5623)
AMPS_T5_LIMIT_CASE(
    "T5-COLLISION-LIMITING-MF-YINSI-001",
    "Yinsi majorant-frequency limiting preserves rate and invariants",
    PIC::MolecularCollisions::ParticleCollisionModel::ModelCellCollisions_mf_Yinsi,
    5624)

#undef AMPS_T5_LIMIT_CASE

// The two simplified selectors explicitly assume one common W and dt and are
// qualified on that declared domain by their algorithm cases. The advanced
// paths own independent update probabilities and receive asymmetric fixtures.
#define AMPS_T5_ASYMMETRIC_CASE(ID, DESCRIPTION, FUNCTION, SEED)           \
  AMPS_TEST_CASE(ID, DESCRIPTION) {                                        \
    RunUnequalWeightTimeFixture(FUNCTION, SEED);                           \
  }

AMPS_T5_ASYMMETRIC_CASE(
    "T5-COLLISION-ASYMMETRIC-NTC-001",
    "NTC preserves directed rates for unequal global W and dt",
    PIC::MolecularCollisions::ParticleCollisionModel::ModelCellCollisions_ntc,
    5630)
AMPS_T5_ASYMMETRIC_CASE(
    "T5-COLLISION-ASYMMETRIC-MF-001",
    "majorant-frequency preserves directed rates for unequal global W and dt",
    PIC::MolecularCollisions::ParticleCollisionModel::ModelCellCollisions_mf,
    5631)
AMPS_T5_ASYMMETRIC_CASE(
    "T5-COLLISION-ASYMMETRIC-MF-YINSI-001",
    "Yinsi majorant-frequency preserves directed rates for unequal global W and dt",
    PIC::MolecularCollisions::ParticleCollisionModel::ModelCellCollisions_mf_Yinsi,
    5632)

#undef AMPS_T5_ASYMMETRIC_CASE

AMPS_TEST_CASE("T5-COLLISION-WEIGHT-CAPABILITY-001",
               "selectors report supported, fixed-only, and absent record-weight laws") {
  namespace Collision = PIC::MolecularCollisions::ParticleCollisionModel;
  AMPS_EXPECT_EQ(Collision::GetIndividualWeightCapability(
                     Collision::COLLISION_SELECTOR_NTC),
                 Collision::INDIVIDUAL_WEIGHT_SUPPORTED);
  AMPS_EXPECT_EQ(Collision::GetIndividualWeightCapability(
                     Collision::COLLISION_SELECTOR_MF),
                 Collision::INDIVIDUAL_WEIGHT_SUPPORTED);
  AMPS_EXPECT_EQ(Collision::GetIndividualWeightCapability(
                     Collision::COLLISION_SELECTOR_NTC_SIMPLIFIED),
                 Collision::INDIVIDUAL_WEIGHT_REQUIRES_FIXED_RECORD_WEIGHT);
  AMPS_EXPECT_EQ(Collision::GetIndividualWeightCapability(
                     Collision::COLLISION_SELECTOR_MF_SIMPLIFIED),
                 Collision::INDIVIDUAL_WEIGHT_REQUIRES_FIXED_RECORD_WEIGHT);
  AMPS_EXPECT_EQ(Collision::GetIndividualWeightCapability(
                     Collision::COLLISION_SELECTOR_MF_YINSI),
                 Collision::INDIVIDUAL_WEIGHT_NOT_IMPLEMENTED);
}

AMPS_TEST_CASE("T5-COLLISION-INDIVIDUAL-WEIGHT-NTC-001",
               "NTC preserves both directed rates for unequal record weights") {
  RunIndividualWeightFixture(
      PIC::MolecularCollisions::ParticleCollisionModel::ModelCellCollisions_ntc,
      5640);
}

AMPS_TEST_CASE("T5-COLLISION-INDIVIDUAL-WEIGHT-MF-001",
               "majorant-frequency preserves directed rates for unequal record weights") {
  RunIndividualWeightFixture(
      PIC::MolecularCollisions::ParticleCollisionModel::ModelCellCollisions_mf,
      5641);
}

AMPS_TEST_CASE("T5-IDF-ROTZ-CONTRACT-001",
               "LB resolves configured rotational collision numbers without mixing") {
  double collision_number = -1.0;
  AMPS_EXPECT_EQ(PIC::IDF::LB::EvaluateRotationalCollisionNumber(
                     _H2O_SPEC_, _H2O_SPEC_, collision_number),
                 PIC::IDF::LB::ROTATIONAL_COLLISION_NUMBER_OK);
  AMPS_EXPECT_NEAR(collision_number,
                   PIC::IDF::RotationZnumber[_H2O_SPEC_], 0.0);
  AMPS_EXPECT_EQ(PIC::IDF::LB::EvaluateRotationalCollisionNumber(
                     _H2O_SPEC_, _O_SPEC_, collision_number),
                 PIC::IDF::LB::ROTATIONAL_COLLISION_NUMBER_OK);
  AMPS_EXPECT_NEAR(collision_number,
                   PIC::IDF::RotationZnumber[_H2O_SPEC_], 0.0);
  AMPS_EXPECT_EQ(PIC::IDF::LB::EvaluateRotationalCollisionNumber(
                     _O_SPEC_, _O_SPEC_, collision_number),
                 PIC::IDF::LB::ROTATIONAL_COLLISION_NUMBER_NO_ROTATIONAL_MODES);
  AMPS_EXPECT_TRUE(std::isinf(collision_number));
}

AMPS_TEST_CASE("T5-IDF-CAPABILITY-001",
               "LB reports RT supported while VT and VV remain explicitly unqualified") {
  AMPS_EXPECT_EQ(PIC::IDF::LB::GetRelaxationCapability(
                     PIC::IDF::LB::RELAXATION_CHANNEL_RT),
                 PIC::IDF::LB::RELAXATION_CAPABILITY_SUPPORTED);
  AMPS_EXPECT_EQ(PIC::IDF::LB::GetRelaxationCapability(
                     PIC::IDF::LB::RELAXATION_CHANNEL_VT),
                 PIC::IDF::LB::RELAXATION_CAPABILITY_UNVERIFIED_TERMINATES);
  AMPS_EXPECT_EQ(PIC::IDF::LB::GetRelaxationCapability(
                     PIC::IDF::LB::RELAXATION_CHANNEL_VV),
                 PIC::IDF::LB::RELAXATION_CAPABILITY_NOT_IMPLEMENTED);
}

AMPS_TEST_CASE("T5-IDF-RT-HOT-COLD-001",
               "LB reaches the same conserved RT equilibrium from either reservoir") {
  ValidateRtHotColdDirections();
}

AMPS_TEST_CASE("T5-IDF-UPDATE-FLAGS-001",
               "LB honors all real/virtual collision-record update permutations") {
  ValidateRtUpdateFlags();
}

AMPS_TEST_CASE("T5-IDF-MODE-FAMILIES-001",
               "monatomic, linear, nonlinear, and vibrational limits match theory") {
  ValidateModeFamiliesAndTemperatureLimits();
}

AMPS_TEST_CASE("T5-COLLISION-VHS-001",
               "generated VHS tables reproduce sigma(g) and Maxwell quadrature") {
  ValidateVariableSphereCrossSection(false);
  AMPS_EXPECT_EQ(PIC::MolecularData::MolecularModels::
                     GetAngularScatteringCapability(),
                 PIC::MolecularData::MolecularModels::
                     ANGULAR_SCATTERING_ISOTROPIC_SUPPORTED);
}

AMPS_TEST_CASE("T5-COLLISION-VSS-001",
               "VSS sigma, pair alpha, angular moments, and invariants are analytical") {
  ValidateVariableSphereCrossSection(true);
  AMPS_EXPECT_NEAR(PIC::MolecularData::MolecularModels::VSS::GetChi(
                       _H2O_SPEC_), 1.70, 0.0);
  AMPS_EXPECT_NEAR(PIC::MolecularData::MolecularModels::VSS::GetChi(
                       _O_SPEC_), 2.10, 0.0);
  AMPS_EXPECT_EQ(PIC::MolecularData::MolecularModels::
                     GetAngularScatteringCapability(),
                 PIC::MolecularData::MolecularModels::
                     ANGULAR_SCATTERING_VSS_SUPPORTED);
  ValidateVssAngularMoments();
}

AMPS_TEST_CASE("T5-IDF-ROTZ-Z2-001",
               "configured Zrot=2 produces an RT probability of one half") {
  ValidateRotationalRelaxationProbability(2.0, 5702);
}

AMPS_TEST_CASE("T5-IDF-ROTZ-Z10-001",
               "configured Zrot=10 produces an RT probability of one tenth") {
  ValidateRotationalRelaxationProbability(10.0, 5710);
}
