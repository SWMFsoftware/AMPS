#include "pic.h"
#include "PhotolyticReactions.h"
#include "src/test/native/test_harness.h"
#include "core_test_state.h"

#include <algorithm>
#include <cmath>

namespace {

double Norm3(const double* value) {
  return std::sqrt(value[0] * value[0] + value[1] * value[1] +
                   value[2] * value[2]);
}

void CheckProductConservation(int products, const double* masses,
                              const double* velocity, double excess_energy) {
  // Product velocities are in the parent centre-of-mass frame.  Momentum must
  // consequently close to zero and the translational kinetic energy must equal
  // the channel excess energy.  Normalized residuals avoid a dimensional
  // absolute tolerance that would be meaningless across electron/molecule
  // mass ratios and reaction energies.
  double momentum[3] = {0.0, 0.0, 0.0};
  double momentum_scale = 0.0;
  double kinetic_energy = 0.0;
  for (int product = 0; product < products; ++product) {
    const double speed = Norm3(velocity + 3 * product);
    momentum_scale += masses[product] * speed;
    kinetic_energy += 0.5 * masses[product] * speed * speed;
    for (int dimension = 0; dimension < 3; ++dimension) {
      momentum[dimension] += masses[product] * velocity[3 * product + dimension];
    }
  }
  AMPS_EXPECT_TRUE(momentum_scale > 0.0);
  AMPS_EXPECT_TRUE(Norm3(momentum) / momentum_scale < 2.0e-12);
  AMPS_EXPECT_TRUE(std::fabs(kinetic_energy - excess_energy) / excess_energy <
                   2.0e-12);
}

}  // namespace

namespace AMPS {
namespace CoreTests {

void LinkPhotochemistryTests() {}

}  // namespace CoreTests
}  // namespace AMPS

AMPS_TEST_CASE("PHYS-PHOTOCHEM-RATES-001",
               "Huebner H2O rate sum, inverse-square scaling, and OH yield") {
  // The active-Sun channel values are a hand-transcribed reference from the
  // Huebner table represented by this production model.  Summing that fixture
  // independently catches omitted/doubled channels; checking 2 AU catches a
  // distance-unit or exponent regression in addition to a table regression.
  const double expected_rates[7] = {1.76e-5, 1.48e-6, 1.91e-6, 8.28e-7,
                                    1.51e-7, 2.21e-8, 4.07e-8};
  double expected_total = 0.0;
  for (int channel = 0; channel < 7; ++channel) expected_total += expected_rates[channel];

  const double at_one_au =
      PhotolyticReactions::H2O::GetTotalReactionRate(1.0 * _AU_);
  const double at_two_au =
      PhotolyticReactions::H2O::GetTotalReactionRate(2.0 * _AU_);
  AMPS_EXPECT_NEAR(at_one_au, expected_total, 1.0e-18);
  AMPS_EXPECT_NEAR(at_two_au, expected_total / 4.0, 3.0e-19);

  // OH is produced by active-Sun channels 0 and 6, exactly once in each.
  const double expected_oh_yield = (expected_rates[0] + expected_rates[6]) /
                                   expected_total;
  AMPS_EXPECT_NEAR(
      PhotolyticReactions::H2O::GetSpeciesReactionYield(_OH_SPEC_),
      expected_oh_yield, 2.0e-14);
}

AMPS_TEST_CASE("PHYS-PHOTOCHEM-DECAY-001",
               "exponential parent decay probability and conditional event time") {
  AMPS::CoreTests::Fixture::PhotolysisLifetimeScope lifetime_scope;
  using PIC::ChemicalReactions::PhotolyticReactions::ConstantTotalLifeTime;
  using PIC::ChemicalReactions::PhotolyticReactions::PhotolyticReaction;

  AMPS_EXPECT_TRUE(ConstantTotalLifeTime != NULL);
  for (int species = 0; species < PIC::nTotalSpecies; ++species)
    ConstantTotalLifeTime[species] = -1.0;
  ConstantTotalLifeTime[_H2O_SPEC_] = 4.0;  // seconds

  // An unavailable lifetime and a zero interval are exact boundary cases, not
  // statistical assertions.  The scheduler must neither decay nor alter the
  // caller's remaining interval.
  double position[3] = {0.0, 0.0, 0.0};
  int species = _O_SPEC_;
  double interval = 1.3;
  AMPS_EXPECT_EQ(PhotolyticReaction(position, -1, species, interval, NULL),
                 _PHOTOLYTIC_REACTIONS_NO_TRANSPHORMATION_);
  AMPS_EXPECT_NEAR(interval, 1.3, 0.0);
  species = _H2O_SPEC_;
  interval = 0.0;
  AMPS_EXPECT_EQ(PhotolyticReaction(position, -1, species, interval, NULL),
                 _PHOTOLYTIC_REACTIONS_NO_TRANSPHORMATION_);
  AMPS_EXPECT_NEAR(interval, 0.0, 0.0);

  const int trials = 100000;
  const double requested_interval = 1.3;
  const double probability = 1.0 - std::exp(-requested_interval / 4.0);
  const double sigma = std::sqrt(trials * probability * (1.0 - probability));
  int decays = 0;
  rnd_seed(1701);
  for (int trial = 0; trial < trials; ++trial) {
    species = _H2O_SPEC_;
    interval = requested_interval;
    const int result = PhotolyticReaction(position, -1, species, interval, NULL);
    if (result == _PHOTOLYTIC_REACTION_OCCURES_) {
      ++decays;
      AMPS_EXPECT_TRUE(interval >= 0.0 && interval <= requested_interval);
    } else {
      AMPS_EXPECT_EQ(result, _PHOTOLYTIC_REACTIONS_NO_TRANSPHORMATION_);
      AMPS_EXPECT_NEAR(interval, requested_interval, 0.0);
    }
  }
  // A six-sigma binomial gate has negligible false-failure probability while
  // detecting use of dt/tau, the survival probability, or a wrong lifetime.
  AMPS_EXPECT_TRUE(std::fabs(decays - trials * probability) < 6.0 * sigma);
}

AMPS_TEST_CASE("PHYS-PHOTOCHEM-KINEMATICS-001",
               "two/three-body photoproduct conservation and isotropic direction") {
  const double two_masses[2] = {_OH__MASS_, _H__MASS_};
  const double three_masses[3] = {_O__MASS_, _H__MASS_, _H__MASS_};
  const double two_energy = 4.04 * eV2J;
  const double three_energy = 0.70 * eV2J;
  double velocity[9];

  rnd_seed(1702);
  for (int sample = 0; sample < 2000; ++sample) {
    PhotolyticReactions::Huebner1992ASS::GetProductVelocity_2(
        two_energy, const_cast<double*>(two_masses), velocity);
    CheckProductConservation(2, two_masses, velocity, two_energy);
    PhotolyticReactions::Huebner1992ASS::GetProductVelocity_3(
        three_energy, const_cast<double*>(three_masses), velocity);
    CheckProductConservation(3, three_masses, velocity, three_energy);
  }

  // For uniform solid angle each Cartesian direction cosine has zero mean and
  // second moment 1/3.  This independent moment oracle specifically prevents
  // sampling polar angle uniformly, which overpopulates the two poles.
  const int samples = 120000;
  double mean[3] = {0.0, 0.0, 0.0};
  double second[3] = {0.0, 0.0, 0.0};
  rnd_seed(1703);
  for (int sample = 0; sample < samples; ++sample) {
    PhotolyticReactions::Huebner1992ASS::GetProductVelocity_2(
        two_energy, const_cast<double*>(two_masses), velocity);
    const double speed = Norm3(velocity);
    for (int dimension = 0; dimension < 3; ++dimension) {
      const double direction = velocity[dimension] / speed;
      mean[dimension] += direction;
      second[dimension] += direction * direction;
    }
  }
  for (int dimension = 0; dimension < 3; ++dimension) {
    mean[dimension] /= samples;
    second[dimension] /= samples;
    AMPS_EXPECT_TRUE(std::fabs(mean[dimension]) < 0.008);
    AMPS_EXPECT_TRUE(std::fabs(second[dimension] - 1.0 / 3.0) < 0.008);
  }
}

AMPS_TEST_CASE("PHYS-PHOTOCHEM-BRANCHING-001",
               "Huebner H2O channel sampler reproduces tabulated branching ratios") {
  // The seven expected probabilities are independently normalized from the
  // published active-Sun rates used by the fixture. A per-channel six-sigma
  // binomial bound tests the production cumulative selector. Product count
  // and identity are checked against a hand-entered neutral/ion table so a
  // correct histogram cannot conceal wrong channel products.
  const double rates[7] = {1.76e-5, 1.48e-6, 1.91e-6, 8.28e-7,
                           1.51e-7, 2.21e-8, 4.07e-8};
  const int expected_products[7][3] = {
      {_OH_SPEC_, _H_SPEC_, -1},
      {_O_SPEC_, _H2_SPEC_, -1},
      {_O_SPEC_, _H_SPEC_, _H_SPEC_},
      {_H2O_PLUS_SPEC_, _ELECTRON_SPEC_, -1},
      {_OH_PLUS_SPEC_, _H_SPEC_, _ELECTRON_SPEC_},
      {_O_PLUS_SPEC_, _H2_SPEC_, _ELECTRON_SPEC_},
      {_OH_SPEC_, _H_PLUS_SPEC_, _ELECTRON_SPEC_},
  };
  double total_rate = 0.0;
  for (int channel = 0; channel < 7; ++channel) total_rate += rates[channel];
  const int samples = 240000;
  int histogram[7] = {0, 0, 0, 0, 0, 0, 0};
  rnd_seed(1704);
  AMPS_EXPECT_TRUE(
      PhotolyticReactions::H2O::Huebner1992ASS::ReactionRateTable != NULL);
  // The generic production selector returns before product kinematics when a
  // channel has no available compiled daughters. Use an intentionally absent
  // product table for the high-statistics branching histogram so this case
  // measures the selector rather than spending its budget generating 240,000
  // already-covered three-body velocity samples.
  int absent_products[7 * 3];
  double unused_masses[7 * 3] = {0.0};
  int returned_products[3];
  double returned_velocities[9];
  int channel_product_count[7] = {2, 2, 3, 2, 3, 3, 3};
  for (int slot = 0; slot < 7 * 3; ++slot) absent_products[slot] = -1;
  for (int sample = 0; sample < samples; ++sample) {
    int channel = -1, product_count = -1;
    PhotolyticReactions::Huebner1992ASS::GenerateReactionProducts(
        channel, product_count, returned_products, returned_velocities,
        PhotolyticReactions::H2O::Huebner1992ASS::ReactionRateTable, 7,
        absent_products,
        channel_product_count, unused_masses, 3, total_rate, NULL);
    AMPS_EXPECT_TRUE(channel >= 0 && channel < 7);
    AMPS_EXPECT_EQ(product_count, 0);
    ++histogram[channel];
  }
  for (int channel = 0; channel < 7; ++channel) {
    const double probability = rates[channel] / total_rate;
    const double expected = samples * probability;
    const double sigma = std::sqrt(
        samples * probability * (1.0 - probability));
    // The rarest channel still has an expectation above 230 events, so the
    // Gaussian six-sigma gate is conservative and deterministic here.
    AMPS_EXPECT_TRUE(std::fabs(histogram[channel] - expected) < 6.0 * sigma);
  }

  int product_channel_seen[7] = {0, 0, 0, 0, 0, 0, 0};
  for (int sample = 0; sample < 12000; ++sample) {
    int channel = -1, product_count = -1;
    int* products = NULL;
    double* velocities = NULL;
    PhotolyticReactions::H2O::GenerateReactionProducts(
        channel, product_count, products, velocities);
    AMPS_EXPECT_TRUE(channel >= 0 && channel < 7);
    // This fixture deliberately compiles neutrals only. Species macros for
    // unavailable ions/electrons are negative, and the production generator
    // omits them while compacting the available product list. Count that
    // hand-entered channel row under the same public availability contract;
    // neutral channels 0--2 remain complete, while photoionization channels
    // exercise the configured missing-daughter behavior.
    int expected_count = 0;
    for (int slot = 0; slot < 3; ++slot)
      if (expected_products[channel][slot] >= 0) ++expected_count;
    AMPS_EXPECT_EQ(product_count, expected_count);
    for (int slot = 0, product = 0; slot < 3; ++slot) {
      if (expected_products[channel][slot] < 0) continue;
      AMPS_EXPECT_EQ(products[product], expected_products[channel][slot]);
      ++product;
    }
    AMPS_EXPECT_TRUE(velocities != NULL);
    ++product_channel_seen[channel];
  }
  for (int channel = 0; channel < 7; ++channel) {
    // With the fixed seed the smallest expected count is about twelve. This
    // explicit check prevents the product-contract loop from silently missing
    // one of the ionization-channel availability branches.
    AMPS_EXPECT_TRUE(product_channel_seen[channel] > 0);
  }

  // Regression for configured missing daughters: IDs and velocities must be
  // compacted together. Force a three-body channel whose only available
  // daughter occupies original slot 1, initialize every return velocity to
  // NaN, and require returned slot 0 to be overwritten with finite kinematics.
  double one_rate[1] = {1.0};
  int sparse_products[3] = {-1, _H_SPEC_, -1};
  int sparse_counts[1] = {3};
  double sparse_masses[3] = {_OH__MASS_, _H__MASS_, _ELECTRON__MASS_};
  double sparse_energy[1] = {2.0 * eV2J};
  int compact_products[3] = {-1, -1, -1};
  double compact_velocity[9];
  for (int entry = 0; entry < 9; ++entry)
    compact_velocity[entry] = std::nan("");
  int forced_channel = -1, compact_count = -1;
  PhotolyticReactions::Huebner1992ASS::GenerateReactionProducts(
      forced_channel, compact_count, compact_products, compact_velocity,
      one_rate, 1, sparse_products, sparse_counts, sparse_masses, 3, 1.0,
      sparse_energy);
  AMPS_EXPECT_EQ(forced_channel, 0);
  AMPS_EXPECT_EQ(compact_count, 1);
  AMPS_EXPECT_EQ(compact_products[0], _H_SPEC_);
  for (int dimension = 0; dimension < 3; ++dimension)
    AMPS_EXPECT_TRUE(std::isfinite(compact_velocity[dimension]));
}
