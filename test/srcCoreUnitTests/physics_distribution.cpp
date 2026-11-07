#include "pic.h"
#include "src/test/native/test_harness.h"

#include <cmath>

namespace AMPS {
namespace CoreTests {

void LinkDistributionTests() {}

}  // namespace CoreTests
}  // namespace AMPS

AMPS_TEST_CASE("PHYS-DISTRIBUTION-MAXWELLIAN-001",
               "Maxwellian sampler reproduces drift, mean speed, and temperature") {
  // For a three-dimensional Maxwell-Boltzmann distribution with component
  // variance kT/m, E[v]=u, E[|v-u|]=sqrt(8kT/(pi*m)), and
  // 2<E_thermal>/(3k)=T.  These formulas are evaluated independently of the
  // production sampler and test the legacy srcUnitTests distribution case
  // without its GoogleTest or accumulator dependencies.
  const int samples = 160000;
  const double temperature = 425.0;  // K
  const double drift[3] = {125.0, -75.0, 40.0};  // m/s
  const double mass = PIC::MolecularData::GetMass(_H2O_SPEC_);
  const double expected_mean_speed =
      std::sqrt(8.0 * Kbol * temperature / (Pi * mass));
  double component_sum[3] = {0.0, 0.0, 0.0};
  double thermal_speed_sum = 0.0;
  double thermal_energy_sum = 0.0;

  rnd_seed(4701);
  for (int sample = 0; sample < samples; ++sample) {
    double velocity[3];
    PIC::Distribution::MaxwellianVelocityDistribution(
        velocity, drift, temperature, _H2O_SPEC_);
    double thermal_speed_squared = 0.0;
    for (int dimension = 0; dimension < 3; ++dimension) {
      component_sum[dimension] += velocity[dimension];
      const double thermal = velocity[dimension] - drift[dimension];
      thermal_speed_squared += thermal * thermal;
    }
    thermal_speed_sum += std::sqrt(thermal_speed_squared);
    thermal_energy_sum += 0.5 * mass * thermal_speed_squared;
  }

  // The standard error of a component mean is sqrt(kT/(mN)); eight standard
  // errors is a stable deterministic gate while still detecting a missing or
  // sign-reversed bulk velocity.
  const double component_error =
      8.0 * std::sqrt(Kbol * temperature / (mass * samples));
  for (int dimension = 0; dimension < 3; ++dimension)
    AMPS_EXPECT_TRUE(std::fabs(component_sum[dimension] / samples -
                               drift[dimension]) < component_error);
  AMPS_EXPECT_TRUE(std::fabs(thermal_speed_sum / samples -
                             expected_mean_speed) /
                       expected_mean_speed < 0.008);
  const double sampled_temperature =
      2.0 * thermal_energy_sum / (3.0 * samples * Kbol);
  AMPS_EXPECT_TRUE(std::fabs(sampled_temperature - temperature) /
                       temperature < 0.008);
}
