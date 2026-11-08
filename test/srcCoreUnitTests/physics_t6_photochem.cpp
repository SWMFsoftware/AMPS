#include "pic.h"
#include "PhotolyticReactions.h"
#include "src/test/native/test_harness.h"
#include "core_test_chemistry.h"
#include "core_test_mesh.h"
#include "core_test_state.h"

#include <cmath>
#include <limits>

namespace {

double Norm3(const double* value) {
  return std::sqrt(value[0] * value[0] + value[1] * value[1] +
                   value[2] * value[2]);
}

void ReinsertUnchanged(
    long int particle, long int& head,
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node) {
  PIC::ChemicalReactions::PhotolyticReactions::
      PhotolyticReactionProcessor_default(particle, head, node);
}

int CellIndex(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node,
              const double* position, int& i, int& j, int& k) {
  AMPS_EXPECT_TRUE(node != NULL && node->block != NULL);
  AMPS_EXPECT_TRUE(PIC::Mesh::mesh->FindCellIndex(
                       const_cast<double*>(position), i, j, k, node, false) >=
                   0);
  return i + _BLOCK_CELLS_X_ * (j + _BLOCK_CELLS_Y_ * k);
}

double PhysicalWeight(
    long int particle, cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node) {
  const int species = PIC::ParticleBuffer::GetI(particle);
  return node->block->GetLocalParticleWeight(species) *
         PIC::ParticleBuffer::GetIndividualStatWeightCorrection(particle);
}

}  // namespace

namespace AMPS {
namespace CoreTests {

void LinkT6PhotochemistryTests() {}

int T6DaughterProcessorStatusCode = T6_DAUGHTER_PROCESSOR_OK;
int T6DaughterProcessorEventCount = 0;

void ResetT6DaughterProcessorState() {
  T6DaughterProcessorStatusCode = T6_DAUGHTER_PROCESSOR_OK;
  T6DaughterProcessorEventCount = 0;
}

void T6PhotochemistryDaughterProcessor(
    long int parent, long int& first_particle_cell,
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node) {
#if _AMPS_CORE_TEST_T6_DAUGHTER_PROCESSOR_MODE_ == \
    _AMPS_CORE_TEST_T6_DAUGHTER_PROCESSOR_MODE__ON_
#if _INDIVIDUAL_PARTICLE_WEIGHT_MODE_ != _INDIVIDUAL_PARTICLE_WEIGHT_ON_
#error The T6 daughter processor profile requires individual particle weights
#endif
  // This callback represents one already-selected H2O photodissociation event
  //   H2O + photon -> OH + H,  E_excess=4.04 eV.
  // Scheduling is independently validated below. The generic production
  // driver owns old-list traversal; this callback owns daughter creation in
  // the new list and deletion of the admitted parent. Other species survive.
  const int parent_species = PIC::ParticleBuffer::GetI(parent);
  if (parent_species != _H2O_SPEC_) {
    ReinsertUnchanged(parent, first_particle_cell, node);
    return;
  }

  if (node == NULL || node->block == NULL) {
    T6DaughterProcessorStatusCode =
        T6_DAUGHTER_PROCESSOR_INVALID_PARENT_WEIGHT;
    ReinsertUnchanged(parent, first_particle_cell, node);
    return;
  }

  const double parent_base_weight =
      node->block->GetLocalParticleWeight(parent_species);
  const double parent_weight = parent_base_weight *
      PIC::ParticleBuffer::GetIndividualStatWeightCorrection(parent);
  if (!std::isfinite(parent_weight) || parent_weight <= 0.0) {
    T6DaughterProcessorStatusCode =
        T6_DAUGHTER_PROCESSOR_INVALID_PARENT_WEIGHT;
    ReinsertUnchanged(parent, first_particle_cell, node);
    return;
  }

  // Isolate the physical first H2O channel while retaining the production
  // Huebner two-body kernel. Masses are kg and excess energy is joules.
  double rate[1] = {1.0};
  int product_table[2] = {_OH_SPEC_, _H_SPEC_};
  int product_count_table[1] = {2};
  double product_mass[2] = {_OH__MASS_, _H__MASS_};
  double excess_energy[1] = {4.04 * eV2J};
  int products[2] = {-1, -1};
  double relative_velocity[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  int channel = -1, product_count = 0;
  PhotolyticReactions::Huebner1992ASS::GenerateReactionProducts(
      channel, product_count, products, relative_velocity, rate, 1,
      product_table, product_count_table, product_mass, 2, 1.0,
      excess_energy);
  if (channel != 0 || product_count != 2 || products[0] != _OH_SPEC_ ||
      products[1] != _H_SPEC_) {
    T6DaughterProcessorStatusCode = T6_DAUGHTER_PROCESSOR_INVALID_CHANNEL;
    ReinsertUnchanged(parent, first_particle_cell, node);
    return;
  }

  double correction[2];
  for (int product = 0; product < product_count; ++product) {
    const double base_weight =
        node->block->GetLocalParticleWeight(products[product]);
    if (!std::isfinite(base_weight) || base_weight <= 0.0) {
      T6DaughterProcessorStatusCode =
          T6_DAUGHTER_PROCESSOR_INVALID_PRODUCT_WEIGHT;
      ReinsertUnchanged(parent, first_particle_cell, node);
      return;
    }
    // One parent molecule produces one molecule/atom of each daughter.
    // c_s=W_parent/W_s makes each daughter carry exactly W_parent physical
    // events even when species-global statistical weights differ.
    correction[product] = parent_weight / base_weight;
  }

  double parent_velocity[3];
  PIC::ParticleBuffer::GetV(parent_velocity, parent);
  for (int product = 0; product < product_count; ++product) {
    const long int daughter =
        PIC::ParticleBuffer::GetNewParticle(first_particle_cell);
    PIC::ParticleBuffer::CloneParticle(daughter, parent);
    PIC::ParticleBuffer::SetI(products[product], daughter);
    PIC::ParticleBuffer::SetIndividualStatWeightCorrection(
        correction[product], daughter);
    double velocity[3];
    for (int dimension = 0; dimension < 3; ++dimension) {
      velocity[dimension] = parent_velocity[dimension] +
          relative_velocity[3 * product + dimension];
    }
    PIC::ParticleBuffer::SetV(velocity, daughter);
  }

  PIC::ParticleBuffer::DeleteParticle(parent);
  ++T6DaughterProcessorEventCount;
#else
  // The symbol remains linkable in ordinary profiles, but the profile does
  // not redirect the production hook here. Explicit calls fail closed.
  ReinsertUnchanged(parent, first_particle_cell, node);
#endif
}

}  // namespace CoreTests
}  // namespace AMPS

AMPS_TEST_CASE("T6-PHOTOCHEM-RATE-DOMAIN-001",
               "all neutral parent rates are idempotent and obey finite inverse-square dilution") {
  // Independent active-Sun channel sums in s^-1 at 1 AU. Initializing twice
  // detects accidental accumulation rather than only checking one sum.
  const int species[6] = {
      _H2O_SPEC_, _O2_SPEC_, _H_SPEC_, _H2_SPEC_, _O_SPEC_, _OH_SPEC_};
  const double expected[6] = {
      1.76e-5 + 1.48e-6 + 1.91e-6 + 8.28e-7 + 1.51e-7 + 2.21e-8 + 4.07e-8,
      2.15e-7 + 6.47e-6 + 9.35e-8 + 1.18e-6 + 3.47e-7,
      1.72e-7,
      1.09e-7 + 8.21e-8 + 1.15e-7 + 2.79e-8,
      6.43e-7,
      1.37e-5 + 1.76e-5 + 2.11e-6 + 6.43e-7};

  PhotolyticReactions::Init();
  for (int parent = 0; parent < 6; ++parent) {
    AMPS_EXPECT_TRUE(PhotolyticReactions::ModelAvailable(species[parent]));
    AMPS_EXPECT_NEAR(PhotolyticReactions::GetTotalReactionRate(
                         species[parent], 1.0 * _AU_),
                     expected[parent], 4.0e-18);
  }
  PhotolyticReactions::Init();
  for (int parent = 0; parent < 6; ++parent) {
    AMPS_EXPECT_NEAR(PhotolyticReactions::GetTotalReactionRate(
                         species[parent], 0.5 * _AU_),
                     4.0 * expected[parent], 1.6e-17);
    AMPS_EXPECT_NEAR(PhotolyticReactions::GetTotalReactionRate(
                         species[parent], 2.0 * _AU_),
                     0.25 * expected[parent], 1.0e-18);
  }
  AMPS_EXPECT_FALSE(PhotolyticReactions::ModelAvailable(_CO2_SPEC_));

  // Invalid distances return NaN, distinguishable from a valid small rate.
  AMPS_EXPECT_TRUE(std::isnan(
      PhotolyticReactions::H2O::GetTotalReactionRate(0.0)));
  AMPS_EXPECT_TRUE(std::isnan(
      PhotolyticReactions::H2O::GetTotalReactionRate(-_AU_)));
  AMPS_EXPECT_TRUE(std::isnan(PhotolyticReactions::H2O::GetTotalReactionRate(
      std::numeric_limits<double>::infinity())));
  AMPS_EXPECT_TRUE(std::isnan(PhotolyticReactions::H2O::GetTotalReactionRate(
      std::numeric_limits<double>::quiet_NaN())));
}

AMPS_TEST_CASE("T6-PHOTOCHEM-SCHEDULER-DOMAIN-001",
               "scheduler rejects disabled, unavailable, and invalid lifetime/time domains") {
  using namespace PIC::ChemicalReactions::PhotolyticReactions;
  AMPS::CoreTests::Fixture::PhotolysisLifetimeScope lifetime_scope;
  double* const original_table = ConstantTotalLifeTime;
  ConstantTotalLifeTime[_H2O_SPEC_] = 9.0;
  Init();
  AMPS_EXPECT_TRUE(ConstantTotalLifeTime == original_table);
  for (int species = 0; species < PIC::nTotalSpecies; ++species)
    AMPS_EXPECT_NEAR(ConstantTotalLifeTime[species], -1.0, 0.0);

  bool allowed = true;
  AMPS_EXPECT_NEAR(TotalLifeTime_default(NULL, -1, -1, allowed, NULL),
                   -1.0, 0.0);
  AMPS_EXPECT_FALSE(allowed);
  AMPS_EXPECT_NEAR(
      TotalLifeTime_default(NULL, PIC::nTotalSpecies, -1, allowed, NULL),
      -1.0, 0.0);
  AMPS_EXPECT_FALSE(allowed);

  double position[3] = {0.0, 0.0, 0.0};
  ConstantTotalLifeTime[_H2O_SPEC_] = 4.0;  // s
  int species = -1;
  double interval = 1.0;
  AMPS_EXPECT_EQ(PhotolyticReaction(position, -1, species, interval, NULL),
                 _PHOTOLYTIC_REACTIONS_NO_TRANSPHORMATION_);
  AMPS_EXPECT_NEAR(interval, 1.0, 0.0);

  const double invalid_intervals[3] = {
      -1.0, std::numeric_limits<double>::infinity(),
      std::numeric_limits<double>::quiet_NaN()};
  for (int domain = 0; domain < 3; ++domain) {
    species = _H2O_SPEC_;
    interval = invalid_intervals[domain];
    AMPS_EXPECT_EQ(PhotolyticReaction(position, -1, species, interval, NULL),
                   _PHOTOLYTIC_REACTIONS_NO_TRANSPHORMATION_);
    if (std::isnan(invalid_intervals[domain]))
      AMPS_EXPECT_TRUE(std::isnan(interval));
    else
      AMPS_EXPECT_TRUE(interval == invalid_intervals[domain]);
  }

  const double invalid_lifetimes[4] = {
      -1.0, 0.0, std::numeric_limits<double>::infinity(),
      std::numeric_limits<double>::quiet_NaN()};
  for (int domain = 0; domain < 4; ++domain) {
    ConstantTotalLifeTime[_H2O_SPEC_] = invalid_lifetimes[domain];
    species = _H2O_SPEC_;
    interval = 1.0;
    AMPS_EXPECT_EQ(PhotolyticReaction(position, -1, species, interval, NULL),
                   _PHOTOLYTIC_REACTIONS_NO_TRANSPHORMATION_);
    AMPS_EXPECT_NEAR(interval, 1.0, 0.0);
  }
}

AMPS_TEST_CASE("T6-PHOTOCHEM-BRANCH-DISTANCE-001",
               "channel branching is distance invariant and invalid tables fail closed") {
  const double reference_rate[7] = {
      1.76e-5, 1.48e-6, 1.91e-6, 8.28e-7, 1.51e-7, 2.21e-8, 4.07e-8};
  double reference_total = 0.0;
  for (int channel = 0; channel < 7; ++channel)
    reference_total += reference_rate[channel];

  // All product IDs are unavailable so this statistical case isolates the
  // production cumulative channel selector from kinematic sampling cost.
  int absent_products[21];
  double unused_mass[21] = {0.0};
  int channel_products[7] = {2, 2, 3, 2, 3, 3, 3};
  for (int slot = 0; slot < 21; ++slot) absent_products[slot] = -1;
  const double distance_au[3] = {0.5, 1.0, 2.0};
  const int samples = 70000;
  for (int distance = 0; distance < 3; ++distance) {
    const double scale = 1.0 / (distance_au[distance] * distance_au[distance]);
    double rate[7], total = 0.0;
    for (int channel = 0; channel < 7; ++channel) {
      rate[channel] = scale * reference_rate[channel];
      total += rate[channel];
    }
    int histogram[7] = {0, 0, 0, 0, 0, 0, 0};
    rnd_seed(6600 + distance);
    for (int sample = 0; sample < samples; ++sample) {
      int channel = -1, count = -1, products[3];
      double velocity[9];
      PhotolyticReactions::Huebner1992ASS::GenerateReactionProducts(
          channel, count, products, velocity, rate, 7, absent_products,
          channel_products, unused_mass, 3, total, NULL);
      AMPS_EXPECT_TRUE(channel >= 0 && channel < 7);
      AMPS_EXPECT_EQ(count, 0);
      ++histogram[channel];
    }
    for (int channel = 0; channel < 7; ++channel) {
      const double probability = reference_rate[channel] / reference_total;
      const double sigma = std::sqrt(
          samples * probability * (1.0 - probability));
      AMPS_EXPECT_TRUE(std::fabs(histogram[channel] - samples * probability) <
                       7.0 * sigma + 2.0);
    }
  }

  // Negative rates are not probabilities. The selector must report the typed
  // sentinel rather than assigning the invalid residual to the final channel.
  int channel = 7, count = 7, products[3];
  double velocity[9];
  double invalid_rate[2] = {1.0, -0.1};
  int product_table[4] = {-1, -1, -1, -1};
  int product_counts[2] = {2, 2};
  double masses[4] = {1.0, 1.0, 1.0, 1.0};
  PhotolyticReactions::Huebner1992ASS::GenerateReactionProducts(
      channel, count, products, velocity, invalid_rate, 2, product_table,
      product_counts, masses, 2, 0.9, NULL);
  AMPS_EXPECT_EQ(channel, -1);
  AMPS_EXPECT_EQ(count, 0);
  double zero_rate[2] = {0.0, 0.0};
  AMPS_EXPECT_NEAR(
      PhotolyticReactions::Huebner1992ASS::GetSpeciesReactionYield(
          _H_SPEC_, zero_rate, 2, product_table, 2),
      0.0, 0.0);
  AMPS_EXPECT_TRUE(std::isnan(
      PhotolyticReactions::Huebner1992ASS::GetSpeciesReactionYield(
          _H_SPEC_, invalid_rate, 2, product_table, 2)));
}

AMPS_TEST_CASE("T6-PHOTOCHEM-GIVEN-PRODUCT-001",
               "conditioned product selection terminates for unavailable daughters and preserves branching") {
  int channel = 9, count = 9;
  int* products = NULL;
  double* velocity = NULL;
  PhotolyticReactions::H2O::GenerateGivenProducts(
      _CO2_SPEC_, channel, count, products, velocity);
  AMPS_EXPECT_EQ(channel, -1);
  AMPS_EXPECT_EQ(count, 0);

  // OH occurs in H2O channels 0 and 6. Conditioning must use
  // k_i/(k_0+k_6), not unconditional fractions. This covers the supported
  // neutral two-body and ionizing three-body product lists.
  const double rate0 = 1.76e-5, rate6 = 4.07e-8;
  const double probability6 = rate6 / (rate0 + rate6);
  const int samples = 50000;
  int channel6 = 0;
  rnd_seed(6610);
  for (int sample = 0; sample < samples; ++sample) {
    PhotolyticReactions::H2O::GenerateGivenProducts(
        _OH_SPEC_, channel, count, products, velocity);
    AMPS_EXPECT_TRUE(channel == 0 || channel == 6);
    AMPS_EXPECT_EQ(products[0], _OH_SPEC_);
    AMPS_EXPECT_TRUE(count == 2 || count == 3);
    for (int product = 0; product < count; ++product)
      for (int dimension = 0; dimension < 3; ++dimension)
        AMPS_EXPECT_TRUE(std::isfinite(velocity[3 * product + dimension]));
    if (channel == 6) ++channel6;
  }
  const double sigma =
      std::sqrt(samples * probability6 * (1.0 - probability6));
  AMPS_EXPECT_TRUE(std::fabs(channel6 - samples * probability6) <
                   7.0 * sigma + 2.0);
}

AMPS_TEST_CASE("T6-PHOTOCHEM-THREE-BODY-ANGULAR-001",
               "three-body events conserve momentum/energy and have Haar-isotropic orientation") {
  const double mass[3] = {_O__MASS_, _H__MASS_, _H__MASS_};
  const double energy = 0.70 * eV2J;
  const int samples = 60000;
  double mean[3] = {0.0, 0.0, 0.0};
  double second[3] = {0.0, 0.0, 0.0};
  double cross[3] = {0.0, 0.0, 0.0};  // xy,yz,zx
  rnd_seed(6620);
  for (int sample = 0; sample < samples; ++sample) {
    double velocity[9];
    PhotolyticReactions::Huebner1992ASS::GetProductVelocity_3(
        energy, const_cast<double*>(mass), velocity);
    double momentum[3] = {0.0, 0.0, 0.0};
    double momentum_scale = 0.0, kinetic = 0.0;
    for (int product = 0; product < 3; ++product) {
      const double speed = Norm3(velocity + 3 * product);
      momentum_scale += mass[product] * speed;
      kinetic += 0.5 * mass[product] * speed * speed;
      for (int dimension = 0; dimension < 3; ++dimension)
        momentum[dimension] +=
            mass[product] * velocity[3 * product + dimension];
    }
    AMPS_EXPECT_TRUE(Norm3(momentum) / momentum_scale < 3.0e-12);
    AMPS_EXPECT_TRUE(std::fabs(kinetic - energy) / energy < 3.0e-12);

    const double speed = Norm3(velocity);
    const double direction[3] = {
        velocity[0] / speed, velocity[1] / speed, velocity[2] / speed};
    for (int dimension = 0; dimension < 3; ++dimension) {
      mean[dimension] += direction[dimension];
      second[dimension] += direction[dimension] * direction[dimension];
    }
    cross[0] += direction[0] * direction[1];
    cross[1] += direction[1] * direction[2];
    cross[2] += direction[2] * direction[0];
  }
  for (int dimension = 0; dimension < 3; ++dimension) {
    mean[dimension] /= samples;
    second[dimension] /= samples;
    cross[dimension] /= samples;
    AMPS_EXPECT_TRUE(std::fabs(mean[dimension]) < 0.012);
    AMPS_EXPECT_TRUE(std::fabs(second[dimension] - 1.0 / 3.0) < 0.012);
    AMPS_EXPECT_TRUE(std::fabs(cross[dimension]) < 0.012);
  }
}

#if _AMPS_CORE_TEST_T6_DAUGHTER_PROCESSOR_MODE_ == \
    _AMPS_CORE_TEST_T6_DAUGHTER_PROCESSOR_MODE__OFF_
AMPS_TEST_CASE("T6-PHOTOCHEM-DRIVER-NO-DAUGHTER-001",
               "whole-mesh default processor preserves mixed parents and both list links") {
  const double position[3] = {0.60, 0.60, 0.60};  // m
  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node =
      AMPS::CoreTests::MoverFixtureNode(position);
  int i, j, k;
  long int& head = node->block->FirstCellParticleTable[
      CellIndex(node, position, i, j, k)];
  AMPS::CoreTests::Fixture::ParticleBufferScope particles(256);
  AMPS::CoreTests::Fixture::EmptyListHeadScope head_scope(head);
  const int parent_species[3] = {_H2O_SPEC_, _O2_SPEC_, _CO2_SPEC_};
  const int each = 24;
  for (int group = 0; group < 3; ++group) {
    for (int particle = 0; particle < each; ++particle) {
      const long int pointer = PIC::ParticleBuffer::GetNewParticle(head);
      PIC::ParticleBuffer::SetI(parent_species[group], pointer);
      PIC::ParticleBuffer::SetX(const_cast<double*>(position), pointer);
    }
  }

  // Each pass traverses native local blocks and reconstructs every cell list.
  // The default callback has no physical daughter capability: all supported
  // and unsupported parent species must survive, with Next and Prev coherent.
  for (int pass = 0; pass < 2; ++pass) {
    PIC::ChemicalReactions::PhotolyticReactions::ExecutePhotochemicalModel();
    int histogram[PIC::nTotalSpecies] = {0};
    int count = 0;
    long int previous = -1;
    for (long int pointer = head; pointer != -1;
         pointer = PIC::ParticleBuffer::GetNext(pointer)) {
      ++count;
      ++histogram[PIC::ParticleBuffer::GetI(pointer)];
      AMPS_EXPECT_EQ(PIC::ParticleBuffer::GetPrev(pointer), previous);
      previous = pointer;
    }
    AMPS_EXPECT_EQ(count, 3 * each);
    for (int group = 0; group < 3; ++group)
      AMPS_EXPECT_EQ(histogram[parent_species[group]], each);
    AMPS_EXPECT_EQ(PIC::ParticleBuffer::GetAllPartNum(), 3 * each);
    PIC::ParticleBuffer::CheckParticleList();
  }
}
#endif

#if _AMPS_CORE_TEST_T6_DAUGHTER_PROCESSOR_MODE_ == \
    _AMPS_CORE_TEST_T6_DAUGHTER_PROCESSOR_MODE__ON_
AMPS_TEST_CASE("T6-PHOTOCHEM-DRIVER-DAUGHTER-001",
               "whole-mesh hook inserts weighted H2O daughters with closed physical budgets") {
  const double position[3] = {0.60, 0.60, 0.60};  // m
  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node =
      AMPS::CoreTests::MoverFixtureNode(position);
  int i, j, k;
  const int cell_index = CellIndex(node, position, i, j, k);
  long int& head = node->block->FirstCellParticleTable[cell_index];
  PIC::Mesh::cDataCenterNode* cell =
      node->block->GetCenterNode(_getCenterNodeLocalNumber(i, j, k));
  AMPS_EXPECT_TRUE(cell != NULL);
  AMPS::CoreTests::Fixture::ParticleBufferScope particles(64);
  AMPS::CoreTests::Fixture::EmptyListHeadScope head_scope(head);
  AMPS::CoreTests::Fixture::MeshCellStateScope cell_scope(node->block, cell);

  // Unequal base weights prove conservation in physical molecules W_s*c_s,
  // not merely in model-particle count.
  PIC::ParticleWeightTimeStep::GlobalParticleWeight[_H2O_SPEC_] = 10.0;
  PIC::ParticleWeightTimeStep::GlobalParticleWeight[_OH_SPEC_] = 2.0;
  PIC::ParticleWeightTimeStep::GlobalParticleWeight[_H_SPEC_] = 8.0;
  const long int parent = PIC::ParticleBuffer::GetNewParticle(head);
  PIC::ParticleBuffer::SetI(_H2O_SPEC_, parent);
  PIC::ParticleBuffer::SetX(const_cast<double*>(position), parent);
  double parent_velocity[3] = {0.0, 0.0, 0.0};
  PIC::ParticleBuffer::SetV(parent_velocity, parent);
  PIC::ParticleBuffer::SetIndividualStatWeightCorrection(0.4, parent);
  const double parent_physical_weight = 4.0;
  const double parent_mass_budget =
      parent_physical_weight * PIC::MolecularData::GetMass(_H2O_SPEC_);
  const double parent_h_budget = 2.0 * parent_physical_weight;
  const double parent_o_budget = parent_physical_weight;

  AMPS::CoreTests::ResetT6DaughterProcessorState();
  rnd_seed(6630);
  PIC::ChemicalReactions::PhotolyticReactions::ExecutePhotochemicalModel();
  AMPS_EXPECT_EQ(AMPS::CoreTests::T6DaughterProcessorStatusCode,
                 AMPS::CoreTests::T6_DAUGHTER_PROCESSOR_OK);
  AMPS_EXPECT_EQ(AMPS::CoreTests::T6DaughterProcessorEventCount, 1);
  AMPS_EXPECT_EQ(PIC::ParticleBuffer::GetAllPartNum(), 2);

  int daughter_count = 0;
  int species_count[PIC::nTotalSpecies] = {0};
  double h_budget = 0.0, o_budget = 0.0, charge_budget = 0.0;
  double mass_budget = 0.0, kinetic_budget = 0.0;
  double momentum[3] = {0.0, 0.0, 0.0};
  double momentum_scale = 0.0;
  long int previous = -1;
  for (long int daughter = head; daughter != -1;
       daughter = PIC::ParticleBuffer::GetNext(daughter)) {
    ++daughter_count;
    AMPS_EXPECT_EQ(PIC::ParticleBuffer::GetPrev(daughter), previous);
    previous = daughter;
    const int species = PIC::ParticleBuffer::GetI(daughter);
    ++species_count[species];
    const double weight = PhysicalWeight(daughter, node);
    AMPS_EXPECT_NEAR(weight, parent_physical_weight, 2.0e-15);
    if (species == _OH_SPEC_) h_budget += weight, o_budget += weight;
    if (species == _H_SPEC_) h_budget += weight;
    const double mass = PIC::MolecularData::GetMass(species);
    const double* velocity = PIC::ParticleBuffer::GetV(daughter);
    const double speed = Norm3(velocity);
    mass_budget += weight * mass;
    charge_budget +=
        weight * PIC::MolecularData::GetElectricCharge(species);
    kinetic_budget += 0.5 * weight * mass * speed * speed;
    momentum_scale += weight * mass * speed;
    for (int dimension = 0; dimension < 3; ++dimension)
      momentum[dimension] += weight * mass * velocity[dimension];
  }
  AMPS_EXPECT_EQ(daughter_count, 2);
  AMPS_EXPECT_EQ(species_count[_OH_SPEC_], 1);
  AMPS_EXPECT_EQ(species_count[_H_SPEC_], 1);
  AMPS_EXPECT_NEAR(h_budget, parent_h_budget, 2.0e-15);
  AMPS_EXPECT_NEAR(o_budget, parent_o_budget, 2.0e-15);
  AMPS_EXPECT_NEAR(charge_budget, 0.0, 0.0);
  // Legacy constants round OH/H to integer atomic masses while H2O uses
  // 18.01528 u. Element counts close exactly; this bounded mass residual
  // records that database approximation rather than hiding it.
  AMPS_EXPECT_TRUE(std::fabs(mass_budget - parent_mass_budget) /
                       parent_mass_budget <
                   2.0e-3);
  AMPS_EXPECT_TRUE(Norm3(momentum) / momentum_scale < 3.0e-12);
  AMPS_EXPECT_TRUE(std::fabs(kinetic_budget -
                             parent_physical_weight * 4.04 * eV2J) /
                       (parent_physical_weight * 4.04 * eV2J) <
                   3.0e-12);
  PIC::ParticleBuffer::CheckParticleList();
}
#endif
