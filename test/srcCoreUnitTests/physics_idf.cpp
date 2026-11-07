#include "pic.h"
#include "src/test/native/test_harness.h"

#include <cmath>
#include <vector>

namespace {

std::vector<PIC::ParticleBuffer::byte> ParticleRecord(int species) {
  std::vector<PIC::ParticleBuffer::byte> record(
      PIC::ParticleBuffer::ParticleDataLength, 0);
  PIC::ParticleBuffer::SetI(species, &record[0]);
  return record;
}

}  // namespace

namespace AMPS {
namespace CoreTests {

void LinkInternalDegreesOfFreedomTests() {}

}  // namespace CoreTests
}  // namespace AMPS

AMPS_TEST_CASE("PHYS-IDF-STORAGE-001",
               "LB rotational and per-mode vibrational particle storage") {
  AMPS_EXPECT_TRUE(PIC::IDF::LB::_ROTATIONAL_ENERGY_OFFSET_ >= 0);
  AMPS_EXPECT_TRUE(PIC::IDF::LB::_VIBRATIONAL_ENERGY_OFFSET_ >
                   PIC::IDF::LB::_ROTATIONAL_ENERGY_OFFSET_);
  AMPS_EXPECT_EQ(PIC::IDF::nTotalRotationalModes[_H2O_SPEC_], 3);
  AMPS_EXPECT_EQ(PIC::IDF::nTotalVibtationalModes[_H2O_SPEC_], 3);
  std::vector<PIC::ParticleBuffer::byte> record = ParticleRecord(_H2O_SPEC_);
  PIC::IDF::LB::SetRotE(1.25, &record[0]);
  AMPS_EXPECT_NEAR(PIC::IDF::LB::GetRotE(&record[0]), 1.25, 0.0);
  double expected_sum = 0.0;
  for (int mode = 0; mode < 3; ++mode) {
    const double energy = 2.0 + 0.75 * mode;
    expected_sum += energy;
    PIC::IDF::LB::SetVibE(energy, mode, &record[0]);
    AMPS_EXPECT_NEAR(PIC::IDF::LB::GetVibE(mode, &record[0]), energy, 0.0);
  }
  AMPS_EXPECT_NEAR(PIC::IDF::LB::GetVibE(-1, &record[0]), expected_sum, 0.0);
}

AMPS_TEST_CASE("PHYS-IDF-ROTATION-001",
               "LB rotational equilibrium mean follows equipartition") {
  // For f classical rotational degrees of freedom, E_rot is Gamma(f/2,kT)
  // and therefore <E_rot>=f kT/2.  H2O has f=3 in the configured fixture.
  const double temperature = 350.0;
  const int samples = 100000;
  std::vector<PIC::ParticleBuffer::byte> record = ParticleRecord(_H2O_SPEC_);
  double sum = 0.0;
  rnd_seed(3701);
  for (int sample = 0; sample < samples; ++sample) {
    PIC::IDF::LB::InitRotTemp(temperature, &record[0]);
    sum += PIC::IDF::LB::GetRotE(&record[0]);
  }
  const double expected = 1.5 * Kbol * temperature;
  AMPS_EXPECT_TRUE(std::fabs(sum / samples - expected) / expected < 0.015);

  // A monatomic species has no rotational storage energy regardless of T.
  record = ParticleRecord(_O_SPEC_);
  PIC::IDF::LB::SetRotE(12.0, &record[0]);
  PIC::IDF::LB::InitRotTemp(temperature, &record[0]);
  AMPS_EXPECT_NEAR(PIC::IDF::LB::GetRotE(&record[0]), 0.0, 0.0);
}

AMPS_TEST_CASE("PHYS-IDF-VIBRATION-001",
               "LB vibrational equilibrium mean follows harmonic-oscillator law") {
  // The classical LB sampler uses Bird's effective vibrational degrees of
  // freedom.  Its per-mode mean must reduce to k*Theta/(exp(Theta/T)-1), the
  // thermal energy above the zero-point level for a harmonic oscillator.
  const double temperature = 350.0;
  const int samples = 100000;
  std::vector<PIC::ParticleBuffer::byte> record = ParticleRecord(_H2O_SPEC_);
  double sum[3] = {0.0, 0.0, 0.0};
  rnd_seed(3702);
  for (int sample = 0; sample < samples; ++sample) {
    PIC::IDF::LB::InitVibTemp(temperature, &record[0]);
    for (int mode = 0; mode < 3; ++mode)
      sum[mode] += PIC::IDF::LB::GetVibE(mode, &record[0]);
  }
  const double theta[3] = {200.0, 300.0, 400.0};
  for (int mode = 0; mode < 3; ++mode) {
    const double expected = Kbol * theta[mode] /
                            (std::exp(theta[mode] / temperature) - 1.0);
    AMPS_EXPECT_TRUE(std::fabs(sum[mode] / samples - expected) / expected < 0.025);
  }
}

AMPS_TEST_CASE("PHYS-IDF-PARTITION-001",
               "LB energy partition conserves energy and follows beta mean") {
  // Microcanonical partition of a fixed energy between f_i internal and f_t
  // translational quadratic degrees of freedom has mean fraction
  // f_i/(f_i+f_t).  This tests both that distribution and exact per-draw
  // conservation; it does not use the production sampler to form its oracle.
  const int samples = 100000;
  const double total = 9.25;
  const double internal_df = 5.0;
  const double translational_df = 3.0;
  double internal_sum = 0.0;
  rnd_seed(3703);
  for (int sample = 0; sample < samples; ++sample) {
    double collision = total;
    double internal = 0.0;
    PIC::IDF::LB::DistributeEnergy(internal, collision, internal_df,
                                   translational_df);
    AMPS_EXPECT_TRUE(internal >= 0.0 && collision >= 0.0);
    AMPS_EXPECT_NEAR(internal + collision, total, 2.0e-15);
    internal_sum += internal;
  }
  const double expected_fraction = internal_df /
                                   (internal_df + translational_df);
  AMPS_EXPECT_TRUE(std::fabs(internal_sum / (samples * total) -
                             expected_fraction) < 0.008);
}

AMPS_TEST_CASE("PHYS-IDF-RT-001",
               "LB rotational-translational transition conserves and partitions energy") {
  // CalculateRTTransition performs the collision-conditioned LB partition in
  // two stages. With omega=0.5, its translational exponent is 1 and the H2O
  // rotational exponent uses both molecules: (f_rot,0+f_rot,1)/2-1=2.
  // Thus the mean total rotational fraction is
  // (2+1)/[(1+1)+(2+1)]=3/5. For two H2O molecules the second beta split is
  // symmetric, so each receives half of that mean.
  const int samples = 120000;
  const double total_energy = 11.0 * Kbol * 300.0;
  const double expected_total_rotational_fraction = 3.0 / 5.0;
  std::vector<PIC::ParticleBuffer::byte> first = ParticleRecord(_H2O_SPEC_);
  std::vector<PIC::ParticleBuffer::byte> second = ParticleRecord(_H2O_SPEC_);
  double rotational_sum[2] = {0.0, 0.0};
  rnd_seed(3704);
  for (int sample = 0; sample < samples; ++sample) {
    PIC::IDF::LB::SetRotE(0.0, &first[0]);
    PIC::IDF::LB::SetRotE(0.0, &second[0]);
    double collision_energy = total_energy;
    bool change[2] = {true, true};
    bool redistributed = false;
    PIC::IDF::LB::CalculateRTTransition(
        &first[0], &second[0], collision_energy, change, redistributed, 0.5);
    const double rotational0 = PIC::IDF::LB::GetRotE(&first[0]);
    const double rotational1 = PIC::IDF::LB::GetRotE(&second[0]);
    AMPS_EXPECT_TRUE(redistributed);
    AMPS_EXPECT_TRUE(collision_energy >= 0.0 && rotational0 >= 0.0 &&
                     rotational1 >= 0.0);
    AMPS_EXPECT_NEAR(collision_energy + rotational0 + rotational1,
                     total_energy, 2.0e-14 * total_energy);
    rotational_sum[0] += rotational0;
    rotational_sum[1] += rotational1;
  }
  const double observed_total_fraction =
      (rotational_sum[0] + rotational_sum[1]) / (samples * total_energy);
  AMPS_EXPECT_TRUE(std::fabs(observed_total_fraction -
                             expected_total_rotational_fraction) < 0.006);
  const double first_share = rotational_sum[0] /
                             (rotational_sum[0] + rotational_sum[1]);
  AMPS_EXPECT_TRUE(std::fabs(first_share - 0.5) < 0.006);

  // The monatomic/molecular branch must put no rotational energy in O. This
  // explicitly replaces the legacy assertion that indexed element 2 of a
  // two-element accumulator and therefore invoked undefined behavior.
  first = ParticleRecord(_O_SPEC_);
  second = ParticleRecord(_H2O_SPEC_);
  PIC::IDF::LB::SetRotE(0.0, &first[0]);
  PIC::IDF::LB::SetRotE(0.0, &second[0]);
  double collision_energy = total_energy;
  bool change[2] = {true, true};
  bool redistributed = false;
  PIC::IDF::LB::CalculateRTTransition(
      &first[0], &second[0], collision_energy, change, redistributed, 0.5);
  AMPS_EXPECT_NEAR(PIC::IDF::LB::GetRotE(&first[0]), 0.0, 0.0);
  AMPS_EXPECT_NEAR(collision_energy + PIC::IDF::LB::GetRotE(&second[0]),
                   total_energy, 2.0e-14 * total_energy);
}

AMPS_TEST_CASE("PHYS-IDF-RT-RATE-001",
               "LB RT relaxation occurrence reproduces its analytic collision number") {
  // The baseline deck generates Z_rot(H2O)=5 and RedistributeEnergy consumes
  // that table, so an eligible collision performs RT redistribution with
  // probability 1/5. Resetting two cold H2O records and the same relative
  // speed on every trial makes a reduced speed an unambiguous event marker;
  // T5's separate Z=2 and Z=10 profiles prove this is not a hidden literal.
  const int trials = 120000;
  const double initial_relative_speed = 1400.0;  // m/s
  std::vector<PIC::ParticleBuffer::byte> first = ParticleRecord(_H2O_SPEC_);
  std::vector<PIC::ParticleBuffer::byte> second = ParticleRecord(_H2O_SPEC_);
  int transitions = 0;
  rnd_seed(3705);
  for (int trial = 0; trial < trials; ++trial) {
    PIC::IDF::LB::SetRotE(0.0, &first[0]);
    PIC::IDF::LB::SetRotE(0.0, &second[0]);
    double relative_speed = initial_relative_speed;
    bool change[2] = {true, true};
    PIC::IDF::LB::RedistributeEnergy(
        &first[0], &second[0], relative_speed, change, NULL);
    if (relative_speed < initial_relative_speed * (1.0 - 1.0e-14))
      ++transitions;
  }
  const double probability = 0.2;
  const double sigma =
      std::sqrt(trials * probability * (1.0 - probability));
  AMPS_EXPECT_TRUE(std::fabs(transitions - trials * probability) < 6.0 * sigma);
}
