#include "pic.h"
#include "src/test/native/test_harness.h"

#include <cmath>

namespace {

double SquaredNorm(const double* value) {
  return value[0] * value[0] + value[1] * value[1] + value[2] * value[2];
}

void ScatterAndCheck(const double* initial0, double mass0,
                     const double* initial1, double mass1) {
  double velocity0[3], velocity1[3], center_of_mass[3], relative[3];
  double initial_momentum[3];
  for (int dimension = 0; dimension < 3; ++dimension) {
    velocity0[dimension] = initial0[dimension];
    velocity1[dimension] = initial1[dimension];
    initial_momentum[dimension] = mass0 * initial0[dimension] +
                                  mass1 * initial1[dimension];
    center_of_mass[dimension] = initial_momentum[dimension] / (mass0 + mass1);
    relative[dimension] = initial1[dimension] - initial0[dimension];
  }
  const double initial_energy =
      0.5 * mass0 * SquaredNorm(initial0) +
      0.5 * mass1 * SquaredNorm(initial1);
  const double relative_speed_squared = SquaredNorm(relative);

  PIC::MolecularCollisions::VelocityScattering::HS::VelocityAfterCollision(
      velocity0, mass0, velocity1, mass1, center_of_mass, relative);

  const double final_energy =
      0.5 * mass0 * SquaredNorm(velocity0) +
      0.5 * mass1 * SquaredNorm(velocity1);
  double final_relative[3];
  for (int dimension = 0; dimension < 3; ++dimension) {
    const double final_momentum = mass0 * velocity0[dimension] +
                                  mass1 * velocity1[dimension];
    const double scale = std::max(1.0, std::fabs(initial_momentum[dimension]));
    AMPS_EXPECT_TRUE(std::fabs(final_momentum - initial_momentum[dimension]) <
                     3.0e-14 * scale);
    final_relative[dimension] = velocity1[dimension] - velocity0[dimension];
  }
  AMPS_EXPECT_TRUE(std::fabs(final_energy - initial_energy) / initial_energy <
                   3.0e-14);
  AMPS_EXPECT_TRUE(std::fabs(SquaredNorm(final_relative) - relative_speed_squared) /
                       relative_speed_squared <
                   3.0e-14);
}

}  // namespace

namespace AMPS {
namespace CoreTests {

void LinkCollisionTests() {}

}  // namespace CoreTests
}  // namespace AMPS

AMPS_TEST_CASE("PHYS-COLLISION-CONSERVATION-001",
               "hard-sphere scattering conserves momentum and kinetic energy") {
  // Unequal masses and asymmetric velocities expose sign, reduced-mass, and
  // centre-of-mass errors.  The first pair exercises the general rotation
  // branch; the second has relative velocity on x and exercises the axial
  // branch selected when sqrt(vrel_y^2+vrel_z^2) is small.
  const double general0[3] = {2.0, -3.5, 7.25};
  const double general1[3] = {-5.0, 11.0, 0.75};
  const double axial0[3] = {-4.0, 2.0, -1.0};
  const double axial1[3] = {9.0, 2.0, -1.0};
  rnd_seed(2701);
  for (int sample = 0; sample < 5000; ++sample) {
    ScatterAndCheck(general0, 2.0, general1, 7.0);
    ScatterAndCheck(axial0, 11.0, axial1, 3.0);
  }
}

AMPS_TEST_CASE("PHYS-COLLISION-ISOTROPY-001",
               "hard-sphere relative velocity is uniform on the sphere") {
  const int samples = 120000;
  double mean[3] = {0.0, 0.0, 0.0};
  double second[3] = {0.0, 0.0, 0.0};
  const double mass0 = 2.0, mass1 = 5.0;
  rnd_seed(2702);
  for (int sample = 0; sample < samples; ++sample) {
    double velocity0[3] = {-1.0, 0.0, 0.0};
    double velocity1[3] = {1.0, 0.0, 0.0};
    double center_of_mass[3] = {
        (mass0 * velocity0[0] + mass1 * velocity1[0]) / (mass0 + mass1),
        0.0, 0.0};
    double relative[3] = {2.0, 0.0, 0.0};
    PIC::MolecularCollisions::VelocityScattering::HS::VelocityAfterCollision(
        velocity0, mass0, velocity1, mass1, center_of_mass, relative);
    double direction[3];
    for (int dimension = 0; dimension < 3; ++dimension)
      direction[dimension] = (velocity1[dimension] - velocity0[dimension]) / 2.0;
    for (int dimension = 0; dimension < 3; ++dimension) {
      mean[dimension] += direction[dimension];
      second[dimension] += direction[dimension] * direction[dimension];
    }
  }
  // Uniform solid angle gives E[n_i]=0 and E[n_i^2]=1/3.  The 0.008
  // tolerance is over five standard errors at this sample size, so the test
  // is stable while decisively rejecting planar or polar-biased scattering.
  for (int dimension = 0; dimension < 3; ++dimension) {
    mean[dimension] /= samples;
    second[dimension] /= samples;
    AMPS_EXPECT_TRUE(std::fabs(mean[dimension]) < 0.008);
    AMPS_EXPECT_TRUE(std::fabs(second[dimension] - 1.0 / 3.0) < 0.008);
  }
}
