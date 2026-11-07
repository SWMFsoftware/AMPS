#include "pic.h"
#include "src/test/native/test_harness.h"
#include "core_test_fields.h"
#include "core_test_mesh.h"
#include "core_test_state.h"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <limits>
#include <vector>

/* Stage-T4 native physics fixtures.
 *
 * These cases call production movers, boundary dispatch, source samplers, and
 * cell-sampling storage through a deliberately small one-rank AMR problem.
 * Analytical expectations are derived in this file rather than copied from
 * production: constant-field Lorentz solutions, exact moving-surface event
 * roots, Maxwell half-space flux, Poisson moments, uniform-volume CDFs, and
 * raw weighted velocity moments. SI units are used throughout. Statistical
 * gates state their sampling allowance locally; deterministic conservation
 * and event geometry use roundoff-scale gates. Each process-global callback
 * or field vector is guarded and restored so case order cannot provide hidden
 * initial conditions. Distributed ownership and unimplemented surface laws
 * are intentionally outside this one-rank fixture and remain explicit gaps in
 * the stage ledger.
 */

namespace AMPS {
namespace CoreTests {

double UniformElectricField[3] = {0.0, 0.0, 0.0};
double UniformMagneticField[3] = {0.0, 0.0, 0.0};

void BorisUniformFieldAcceleration(
    double* acceleration, double* rotation, int species, long int,
    double*, double*, cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>*) {
  const double charge_to_mass =
      PIC::MolecularData::GetElectricCharge(species) /
      PIC::MolecularData::GetMass(species);
  for (int dimension = 0; dimension < 3; ++dimension) {
    acceleration[dimension] =
        charge_to_mass * UniformElectricField[dimension];
    // The production Boris convention is dv/dt=a+Omega x v with
    // Omega=-(q/m)B, which is algebraically q(E+v x B)/m.
    rotation[dimension] =
        -charge_to_mass * UniformMagneticField[dimension];
  }
}

void RelativisticBorisUniformFields(double* electric_field,
                                    double* magnetic_field) {
  std::copy(UniformElectricField, UniformElectricField + 3, electric_field);
  std::copy(UniformMagneticField, UniformMagneticField + 3, magnetic_field);
}

void LinkT4Tests() {}

}  // namespace CoreTests
}  // namespace AMPS

namespace {

using AMPS::CoreTests::Fixture::EmptyListHeadScope;
using AMPS::CoreTests::Fixture::MeshCellStateScope;
using AMPS::CoreTests::Fixture::ParticleBufferScope;

struct UniformFieldScope {
  double electric[3];
  double magnetic[3];
  UniformFieldScope() {
    std::copy(AMPS::CoreTests::UniformElectricField,
              AMPS::CoreTests::UniformElectricField + 3, electric);
    std::copy(AMPS::CoreTests::UniformMagneticField,
              AMPS::CoreTests::UniformMagneticField + 3, magnetic);
  }
  ~UniformFieldScope() {
    std::copy(electric, electric + 3,
              AMPS::CoreTests::UniformElectricField);
    std::copy(magnetic, magnetic + 3,
              AMPS::CoreTests::UniformMagneticField);
  }
};

struct FlightResult {
  double position[3];
  double velocity[3];
};

enum MoverKind { NONRELATIVISTIC_BORIS, RELATIVISTIC_BORIS };

void DetachMovedParticle(long int particle,
                         cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node) {
  const double* position = PIC::ParticleBuffer::GetX(particle);
  int i, j, k;
  AMPS_EXPECT_TRUE(PIC::Mesh::mesh->FindCellIndex(
                       const_cast<double*>(position), i, j, k, node, false) >= 0);
#if _COMPILATION_MODE_ == _COMPILATION_MODE__MPI_
  long int& head = node->block->tempParticleMovingListTable[
      i + _BLOCK_CELLS_X_ * (j + _BLOCK_CELLS_Y_ * k)];
  AMPS_EXPECT_EQ(head, particle);
  head = PIC::ParticleBuffer::GetNext(particle);
  if (head != -1) PIC::ParticleBuffer::SetPrev(-1, head);
  PIC::ParticleBuffer::SetNext(-1, particle);
  PIC::ParticleBuffer::SetPrev(-1, particle);
#else
#error The T4 one-rank fixture requires the configured MPI particle-list mode
#endif
}

FlightResult AdvanceParticle(int species, const double* initial_position,
                             const double* initial_velocity, double total_time,
                             int steps, MoverKind kind) {
  AMPS_EXPECT_TRUE(steps > 0 && total_time > 0.0);
  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node =
      AMPS::CoreTests::MoverFixtureNode(initial_position);
  ParticleBufferScope particles(64);
  const long int particle = PIC::ParticleBuffer::GetNewParticle(true);
  PIC::ParticleBuffer::SetI(species, particle);
  PIC::ParticleBuffer::SetX(const_cast<double*>(initial_position), particle);
  PIC::ParticleBuffer::SetV(const_cast<double*>(initial_velocity), particle);
  const double dt = total_time / steps;
  for (int step = 0; step < steps; ++step) {
    const int status = kind == NONRELATIVISTIC_BORIS
        ? PIC::Mover::Boris(particle, dt, node)
        : PIC::Mover::Relativistic::Boris(particle, dt, node);
    AMPS_EXPECT_EQ(status, _PARTICLE_MOTION_FINISHED_);
    node = PIC::Mesh::mesh->findTreeNode(PIC::ParticleBuffer::GetX(particle),
                                         node);
    AMPS_EXPECT_TRUE(node != NULL && node->block != NULL);
    DetachMovedParticle(particle, node);
  }
  FlightResult result;
  PIC::ParticleBuffer::GetX(result.position, particle);
  PIC::ParticleBuffer::GetV(result.velocity, particle);
  return result;
}

double Norm2(const double* value) {
  return value[0] * value[0] + value[1] * value[1] +
         value[2] * value[2];
}

double VelocityError(const FlightResult& result, const double* expected) {
  double error = 0.0;
  for (int dimension = 0; dimension < 3; ++dimension)
    error += std::pow(result.velocity[dimension] - expected[dimension], 2);
  return std::sqrt(error);
}

int absorbing_count = 0;
int AbsorbSphere(int, long int, double*, double*, double&, void*, void*) {
  ++absorbing_count;
  return _PARTICLE_DELETED_ON_THE_FACE_;
}

int tangent_count = 0;
int CountAndRejectSphere(int, long int, double*, double*, double&, void*, void*) {
  ++tangent_count;
  return _PARTICLE_REJECTED_ON_THE_FACE_;
}

int moving_finder_count = 0;
int moving_interaction_count = 0;
const double moving_plane_origin = 0.0;
const double moving_plane_speed = 0.2;
const double moving_sphere_origin[3] = {0.0, 0.6, 0.6};  // m at t=0
const double moving_sphere_velocity[3] = {0.1, 0.0, 0.0}; // m/s
const double moving_sphere_radius = 0.2;                   // m

bool FindMovingPlane(long int, const double* position, const double* velocity,
                     const double* acceleration, double elapsed,
                     double maximum_time, double& hit_time,
                     double* hit_position, double* outward_normal,
                     double* surface_velocity) {
  ++moving_finder_count;
  // The exact fixture is inertial and force-free. The plane is x=u*t and its
  // outward normal faces the incident negative-x half-space. Solve
  // x0+v*t = u*(elapsed+t); a non-positive root is not a new event.
  AMPS_EXPECT_NEAR(acceleration[0], 0.0, 0.0);
  const double relative_speed = velocity[0] - moving_plane_speed;
  if (relative_speed == 0.0) return false;
  const double root =
      (moving_plane_origin + moving_plane_speed * elapsed - position[0]) /
      relative_speed;
  if (!(root > 0.0) || root > maximum_time) return false;
  hit_time = root;
  for (int dimension = 0; dimension < 3; ++dimension) {
    hit_position[dimension] = position[dimension] + velocity[dimension] * root;
    outward_normal[dimension] = 0.0;
    surface_velocity[dimension] = 0.0;
  }
  outward_normal[0] = -1.0;
  surface_velocity[0] = moving_plane_speed;
  return true;
}

int ReflectMovingPlane(long int, double*, double* velocity,
                       const double* outward_normal,
                       const double* surface_velocity,
                       cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>*) {
  ++moving_interaction_count;
  double normal_speed = 0.0;
  for (int dimension = 0; dimension < 3; ++dimension)
    normal_speed += (velocity[dimension] - surface_velocity[dimension]) *
                    outward_normal[dimension];
  for (int dimension = 0; dimension < 3; ++dimension)
    velocity[dimension] -= 2.0 * normal_speed * outward_normal[dimension];
  return _PARTICLE_REJECTED_ON_THE_FACE_;
}

bool FindMovingSphere(long int, const double* position, const double* velocity,
                      const double* acceleration, double elapsed,
                      double maximum_time, double& hit_time,
                      double* hit_position, double* outward_normal,
                      double* surface_velocity) {
  ++moving_finder_count;
  // Work in the inertial frame translated with the prescribed sphere centre.
  // With zero acceleration, r(t)=r0+w*t where r0=x-c(elapsed) [m] and
  // w=v-u [m/s]. The first positive root of
  // |r0+w*t|^2-R^2=0 is the physical entry event. A negative discriminant or
  // root outside the remaining substep is a typed no-hit result (false).
  double relative_position[3], relative_velocity[3];
  double a = 0.0, b_half = 0.0, c = -moving_sphere_radius * moving_sphere_radius;
  for (int dimension = 0; dimension < 3; ++dimension) {
    AMPS_EXPECT_NEAR(acceleration[dimension], 0.0, 0.0);
    const double centre = moving_sphere_origin[dimension] +
        moving_sphere_velocity[dimension] * elapsed;
    relative_position[dimension] = position[dimension] - centre;
    relative_velocity[dimension] =
        velocity[dimension] - moving_sphere_velocity[dimension];
    a += relative_velocity[dimension] * relative_velocity[dimension];
    b_half += relative_position[dimension] * relative_velocity[dimension];
    c += relative_position[dimension] * relative_position[dimension];
  }
  if (!(a > 0.0)) return false;
  const double discriminant = b_half * b_half - a * c;
  if (discriminant < 0.0) return false;
  const double root0 = (-b_half - std::sqrt(discriminant)) / a;
  const double root1 = (-b_half + std::sqrt(discriminant)) / a;
  const double root = root0 > 0.0 ? root0 : root1;
  if (!(root > 0.0) || root > maximum_time) return false;

  hit_time = root;
  double normal_norm2 = 0.0;
  for (int dimension = 0; dimension < 3; ++dimension) {
    hit_position[dimension] = position[dimension] + velocity[dimension] * root;
    const double centre = moving_sphere_origin[dimension] +
        moving_sphere_velocity[dimension] * (elapsed + root);
    outward_normal[dimension] = hit_position[dimension] - centre;
    normal_norm2 += outward_normal[dimension] * outward_normal[dimension];
    surface_velocity[dimension] = moving_sphere_velocity[dimension];
  }
  AMPS_EXPECT_TRUE(normal_norm2 > 0.0);
  const double inverse_norm = 1.0 / std::sqrt(normal_norm2);
  for (int dimension = 0; dimension < 3; ++dimension)
    outward_normal[dimension] *= inverse_norm;
  return true;
}

double IndependentInwardMaxwellianFlux(double density, double temperature,
                                       double normal_bulk, int species) {
  const double thermal_sigma = std::sqrt(
      Kbol * temperature / PIC::MolecularData::GetMass(species));
  const double z = normal_bulk / (std::sqrt(2.0) * thermal_sigma);
  return density * (thermal_sigma / std::sqrt(2.0 * Pi) *
      std::exp(-0.5 * std::pow(normal_bulk / thermal_sigma, 2)) -
      0.5 * normal_bulk * std::erfc(z));
}

}  // namespace

AMPS_TEST_CASE("T4-MOVER-ELECTRIC-CONVERGENCE-001",
               "charged Boris reproduces uniform-E velocity and converges to position") {
  UniformFieldScope fields;
  AMPS::CoreTests::UniformElectricField[0] = 1.0e-10;  // V/m
  AMPS::CoreTests::UniformElectricField[1] = 0.0;
  AMPS::CoreTests::UniformElectricField[2] = 0.0;
  std::fill(AMPS::CoreTests::UniformMagneticField,
            AMPS::CoreTests::UniformMagneticField + 3, 0.0);
  const double x0[3] = {0.45, 0.65, 0.55};
  const double v0[3] = {0.08, -0.01, 0.015};
  const double time = 2.0;
  const double acceleration =
      PIC::MolecularData::GetElectricCharge(_H_PLUS_SPEC_) /
      PIC::MolecularData::GetMass(_H_PLUS_SPEC_) *
      AMPS::CoreTests::UniformElectricField[0];
  double expected_velocity[3] = {
      v0[0] + acceleration * time, v0[1], v0[2]};
  double previous_error = std::numeric_limits<double>::infinity();
  for (int level = 0; level < 3; ++level) {
    const int steps = 10 << level;
    const FlightResult result = AdvanceParticle(
        _H_PLUS_SPEC_, x0, v0, time, steps, NONRELATIVISTIC_BORIS);
    AMPS_EXPECT_TRUE(VelocityError(result, expected_velocity) < 2.0e-14);
    const double exact_x = x0[0] + v0[0] * time +
                           0.5 * acceleration * time * time;
    const double error = std::fabs(result.position[0] - exact_x);
    if (level != 0) {
      // AMPS stores velocity at the leapfrog half level and advances x with
      // the updated value. Against a same-time analytic position this is
      // first-order; halving dt must halve the bounded staggering error.
      AMPS_EXPECT_TRUE(previous_error / error > 1.85 &&
                       previous_error / error < 2.15);
    }
    previous_error = error;
  }

  // Field/charge reversal and q/m scaling are checked at the exact force hook
  // independently of integration error.
  double ap[3], rp[3], ae[3], re[3];
  AMPS::CoreTests::BorisUniformFieldAcceleration(
      ap, rp, _H_PLUS_SPEC_, -1, NULL, NULL, NULL);
  AMPS::CoreTests::BorisUniformFieldAcceleration(
      ae, re, _ELECTRON_SPEC_, -1, NULL, NULL, NULL);
  AMPS_EXPECT_TRUE(ap[0] > 0.0 && ae[0] < 0.0);
  const double expected_ratio =
      PIC::MolecularData::GetMass(_H_PLUS_SPEC_) /
      PIC::MolecularData::GetMass(_ELECTRON_SPEC_);
  AMPS_EXPECT_TRUE(std::fabs(std::fabs(ae[0] / ap[0]) - expected_ratio) /
                       expected_ratio < 1.0e-12);
}

AMPS_TEST_CASE("T4-MOVER-MAGNETIC-CONVERGENCE-001",
               "charged Boris conserves energy and converges to uniform-B gyromotion") {
  UniformFieldScope fields;
  std::fill(AMPS::CoreTests::UniformElectricField,
            AMPS::CoreTests::UniformElectricField + 3, 0.0);
  AMPS::CoreTests::UniformMagneticField[0] = 0.0;
  AMPS::CoreTests::UniformMagneticField[1] = 0.0;
  AMPS::CoreTests::UniformMagneticField[2] = 1.0e-8;  // T
  const double x0[3] = {0.55, 0.65, 0.55};
  const double v0[3] = {0.08, 0.0, 0.02};
  const double omega =
      PIC::MolecularData::GetElectricCharge(_H_PLUS_SPEC_) *
      AMPS::CoreTests::UniformMagneticField[2] /
      PIC::MolecularData::GetMass(_H_PLUS_SPEC_);
  const double time = 1.3;
  const double expected_velocity[3] = {
      v0[0] * std::cos(omega * time),
      -v0[0] * std::sin(omega * time), v0[2]};
  double previous_error = std::numeric_limits<double>::infinity();
  for (int level = 0; level < 3; ++level) {
    const FlightResult result = AdvanceParticle(
        _H_PLUS_SPEC_, x0, v0, time, 8 << level,
        NONRELATIVISTIC_BORIS);
    AMPS_EXPECT_NEAR(Norm2(result.velocity), Norm2(v0), 3.0e-16);
    AMPS_EXPECT_NEAR(result.velocity[2], v0[2], 2.0e-15);
    const double error = VelocityError(result, expected_velocity);
    if (level != 0)
      AMPS_EXPECT_TRUE(previous_error / error > 3.5);
    previous_error = error;
  }

  // Reversing B reverses gyrophase while preserving speed.
  AMPS::CoreTests::UniformMagneticField[2] *= -1.0;
  const FlightResult reversed = AdvanceParticle(
      _H_PLUS_SPEC_, x0, v0, time, 64, NONRELATIVISTIC_BORIS);
  AMPS_EXPECT_TRUE(reversed.velocity[1] > 0.0);
  AMPS_EXPECT_NEAR(Norm2(reversed.velocity), Norm2(v0), 3.0e-16);
}

AMPS_TEST_CASE("T4-MOVER-EXB-DRIFT-001",
               "charged Boris converges to the analytic crossed-field drift") {
  UniformFieldScope fields;
  AMPS::CoreTests::UniformElectricField[0] = 2.0e-10;
  AMPS::CoreTests::UniformElectricField[1] = 0.0;
  AMPS::CoreTests::UniformElectricField[2] = 0.0;
  AMPS::CoreTests::UniformMagneticField[0] = 0.0;
  AMPS::CoreTests::UniformMagneticField[1] = 0.0;
  AMPS::CoreTests::UniformMagneticField[2] = 1.0e-8;
  const double drift_y = -AMPS::CoreTests::UniformElectricField[0] /
                          AMPS::CoreTests::UniformMagneticField[2];
  const double omega =
      PIC::MolecularData::GetElectricCharge(_H_PLUS_SPEC_) *
      AMPS::CoreTests::UniformMagneticField[2] /
      PIC::MolecularData::GetMass(_H_PLUS_SPEC_);
  const double time = 2.0 * Pi / omega;
  const double x0[3] = {0.55, 0.75, 0.55};
  const double v0[3] = {0.02, drift_y, 0.0};
  double previous_error = std::numeric_limits<double>::infinity();
  for (int level = 0; level < 3; ++level) {
    const FlightResult result = AdvanceParticle(
        _H_PLUS_SPEC_, x0, v0, time, 64 << level,
        NONRELATIVISTIC_BORIS);
    const double observed_drift = (result.position[1] - x0[1]) / time;
    const double error = std::fabs(observed_drift - drift_y);
    if (level != 0) AMPS_EXPECT_TRUE(error < previous_error * 0.7);
    previous_error = error;
    AMPS_EXPECT_TRUE(error < 2.0e-3);
  }
}

AMPS_TEST_CASE("T4-MOVER-RELATIVISTIC-BORIS-001",
               "relativistic Boris preserves gamma and converges to gyrofrequency") {
  UniformFieldScope fields;
  std::fill(AMPS::CoreTests::UniformElectricField,
            AMPS::CoreTests::UniformElectricField + 3, 0.0);
  AMPS::CoreTests::UniformMagneticField[0] = 0.0;
  AMPS::CoreTests::UniformMagneticField[1] = 0.0;
  AMPS::CoreTests::UniformMagneticField[2] = 4.0;  // T
  const double speed = 0.2 * SpeedOfLight;
  const double gamma = 1.0 / std::sqrt(1.0 - 0.04);
  const double omega =
      PIC::MolecularData::GetElectricCharge(_H_PLUS_SPEC_) *
      AMPS::CoreTests::UniformMagneticField[2] /
      (gamma * PIC::MolecularData::GetMass(_H_PLUS_SPEC_));
  const double time = 0.5 * Pi / omega;
  const double x0[3] = {0.55, 0.55, 0.55};
  const double v0[3] = {0.0, speed, 0.0};
  const double expected_velocity[3] = {speed, 0.0, 0.0};
  double previous_error = std::numeric_limits<double>::infinity();
  for (int level = 0; level < 3; ++level) {
    const FlightResult result = AdvanceParticle(
        _H_PLUS_SPEC_, x0, v0, time, 8 << level, RELATIVISTIC_BORIS);
    AMPS_EXPECT_TRUE(std::fabs(std::sqrt(Norm2(result.velocity)) - speed) /
                         speed < 2.0e-14);
    const double observed_gamma = 1.0 / std::sqrt(
        1.0 - Norm2(result.velocity) / (SpeedOfLight * SpeedOfLight));
    AMPS_EXPECT_NEAR(observed_gamma, gamma, 3.0e-14);
    const double error = VelocityError(result, expected_velocity) / speed;
    if (level != 0) AMPS_EXPECT_TRUE(previous_error / error > 3.0);
    previous_error = error;
  }
}

AMPS_TEST_CASE("T4-GYRO-DRIFTS-001",
               "guiding-center RHS reproduces parallel, ExB, grad-B, and curvature terms") {
  const double mass = PIC::MolecularData::GetMass(_H_PLUS_SPEC_);
  const double charge = PIC::MolecularData::GetElectricCharge(_H_PLUS_SPEC_);
  const double magnetic_strength = 2.0;  // T
  const double grad_magnitude = 0.5;     // T/m, dBz/dx
  const double curvature_gradient = 0.25; // T/m, dBx/dz
  const double mirror_gradient = 0.1;    // T/m, dBz/dz
  const double magnetic_moment =
      0.3 * charge * magnetic_strength / grad_magnitude;
  const double parallel_velocity = std::sqrt(
      0.4 * charge * magnetic_strength * magnetic_strength /
      (mass * curvature_gradient));
  const double electric[3] = {4.0, 0.0, 3.0};
  const double magnetic[3] = {0.0, 0.0, magnetic_strength};
  double gradient[9] = {0.0, 0.0, curvature_gradient,
                        0.0, 0.0, 0.0,
                        grad_magnitude, 0.0, mirror_gradient};
  double abs_b, direction[3], drift[3], parallel_acceleration;
  AMPS_EXPECT_EQ(PIC::GYROKINETIC::EvaluateRhsFromFields(
                     electric, magnetic, gradient, parallel_velocity,
                     magnetic_moment, _H_PLUS_SPEC_, abs_b, direction, drift,
                     parallel_acceleration),
                 PIC::GYROKINETIC::GUIDING_CENTER_RHS_OK);
  AMPS_EXPECT_NEAR(abs_b, magnetic_strength, 0.0);
  AMPS_EXPECT_NEAR(direction[0], 0.0, 0.0);
  AMPS_EXPECT_NEAR(direction[1], 0.0, 0.0);
  AMPS_EXPECT_NEAR(direction[2], 1.0, 0.0);
  // E x B contributes -2 y, grad-B +0.3 y, curvature +0.4 y.
  AMPS_EXPECT_NEAR(drift[0], 0.0, 2.0e-14);
  AMPS_EXPECT_NEAR(drift[1], -1.3, 2.0e-13);
  AMPS_EXPECT_NEAR(drift[2], 0.0, 2.0e-14);
  const double expected_parallel = charge / mass * electric[2] -
      magnetic_moment / mass * mirror_gradient;
  AMPS_EXPECT_TRUE(std::fabs(parallel_acceleration - expected_parallel) /
                       std::fabs(expected_parallel) < 2.0e-14);
}

AMPS_TEST_CASE("T4-GYRO-VALIDITY-001",
               "guiding-center RHS returns typed neutral, zero-B, and nonfinite states") {
  const double electric[3] = {0.0, 0.0, 0.0};
  double magnetic[3] = {0.0, 0.0, 1.0};
  double gradient[9] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                        0.0, 0.0, 0.0};
  double abs_b, direction[3], drift[3], acceleration;
  AMPS_EXPECT_EQ(PIC::GYROKINETIC::EvaluateRhsFromFields(
                     electric, magnetic, gradient, 1.0, 0.0, _H2O_SPEC_,
                     abs_b, direction, drift, acceleration),
                 PIC::GYROKINETIC::GUIDING_CENTER_RHS_NEUTRAL_SPECIES);
  magnetic[2] = 0.0;
  AMPS_EXPECT_EQ(PIC::GYROKINETIC::EvaluateRhsFromFields(
                     electric, magnetic, gradient, 1.0, 0.0, _H_PLUS_SPEC_,
                     abs_b, direction, drift, acceleration),
                 PIC::GYROKINETIC::GUIDING_CENTER_RHS_ZERO_MAGNETIC_FIELD);
  magnetic[2] = 1.0;
  gradient[4] = std::numeric_limits<double>::quiet_NaN();
  AMPS_EXPECT_EQ(PIC::GYROKINETIC::EvaluateRhsFromFields(
                     electric, magnetic, gradient, 1.0, 0.0, _H_PLUS_SPEC_,
                     abs_b, direction, drift, acceleration),
                 PIC::GYROKINETIC::GUIDING_CENTER_RHS_NONFINITE_INPUT);
  AMPS_EXPECT_NEAR(abs_b, 0.0, 0.0);
  AMPS_EXPECT_NEAR(Norm2(drift), 0.0, 0.0);
}

AMPS_TEST_CASE("T4-BOUNDARY-ABSORB-OUTER-001",
               "sphere absorption and outer-domain loss delete exactly one particle") {
  cInternalSphericalData* sphere = AMPS::CoreTests::MoverFixtureSphere();
  AMPS_EXPECT_TRUE(sphere != NULL);
  AMPS::CoreTests::Fixture::ValueScope<decltype(sphere->ParticleSphereInteraction)>
      callback_scope(sphere->ParticleSphereInteraction);
  AMPS::CoreTests::Fixture::ValueScope<int> count_scope(absorbing_count);
  sphere->ParticleSphereInteraction = AbsorbSphere;
  absorbing_count = 0;
  {
    const double x[3] = {-0.6, 0.0, 0.0};
    const double v[3] = {1.0, 0.0, 0.0};
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node =
        AMPS::CoreTests::MoverFixtureNode(x);
    ParticleBufferScope particles(16);
    const long int particle = PIC::ParticleBuffer::GetNewParticle(true);
    PIC::ParticleBuffer::SetI(_H2O_SPEC_, particle);
    PIC::ParticleBuffer::SetX(const_cast<double*>(x), particle);
    PIC::ParticleBuffer::SetV(const_cast<double*>(v), particle);
    AMPS_EXPECT_EQ(PIC::Mover::UniformWeight_UniformTimeStep_noForce_TraceTrajectory(
                       particle, 1.0, node),
                   _PARTICLE_LEFT_THE_DOMAIN_);
    AMPS_EXPECT_EQ(absorbing_count, 1);
    AMPS_EXPECT_TRUE(!PIC::ParticleBuffer::IsParticleAllocated(particle));
  }
  {
    const double x[3] = {0.9, 0.6, 0.6};
    const double v[3] = {0.4, 0.0, 0.0};
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node =
        AMPS::CoreTests::MoverFixtureNode(x);
    ParticleBufferScope particles(16);
    const long int particle = PIC::ParticleBuffer::GetNewParticle(true);
    PIC::ParticleBuffer::SetI(_H2O_SPEC_, particle);
    PIC::ParticleBuffer::SetX(const_cast<double*>(x), particle);
    PIC::ParticleBuffer::SetV(const_cast<double*>(v), particle);
    AMPS_EXPECT_EQ(PIC::Mover::Simple(particle, 1.0, node),
                   _PARTICLE_LEFT_THE_DOMAIN_);
    AMPS_EXPECT_TRUE(!PIC::ParticleBuffer::IsParticleAllocated(particle));
    AMPS_EXPECT_EQ(absorbing_count, 1);
  }
}

AMPS_TEST_CASE("T4-BOUNDARY-TANGENT-001",
               "near-tangent sphere trajectory is a no-hit with intact list state") {
  cInternalSphericalData* sphere = AMPS::CoreTests::MoverFixtureSphere();
  AMPS::CoreTests::Fixture::ValueScope<decltype(sphere->ParticleSphereInteraction)>
      callback_scope(sphere->ParticleSphereInteraction);
  AMPS::CoreTests::Fixture::ValueScope<int> count_scope(tangent_count);
  sphere->ParticleSphereInteraction = CountAndRejectSphere;
  tangent_count = 0;
  const double x[3] = {-0.6, 0.2001, 0.0};
  const double v[3] = {1.0, 0.0, 0.0};
  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node =
      AMPS::CoreTests::MoverFixtureNode(x);
  ParticleBufferScope particles(16);
  const long int particle = PIC::ParticleBuffer::GetNewParticle(true);
  PIC::ParticleBuffer::SetI(_H2O_SPEC_, particle);
  PIC::ParticleBuffer::SetX(const_cast<double*>(x), particle);
  PIC::ParticleBuffer::SetV(const_cast<double*>(v), particle);
  AMPS_EXPECT_EQ(PIC::Mover::UniformWeight_UniformTimeStep_noForce_TraceTrajectory(
                     particle, 1.0, node),
                 _PARTICLE_MOTION_FINISHED_);
  AMPS_EXPECT_EQ(tangent_count, 0);
  // The trajectory tracer deliberately nudges block/surface events by the
  // mesh tolerance to prevent a zero-time rediscovery; this miss remains
  // geometric, not an exact arithmetic free-flight path.
  AMPS_EXPECT_NEAR(PIC::ParticleBuffer::GetX(particle)[0], 0.4, 2.0e-5);
  node = PIC::Mesh::mesh->findTreeNode(PIC::ParticleBuffer::GetX(particle));
  DetachMovedParticle(particle, node);
}

AMPS_TEST_CASE("T4-BOUNDARY-START-INSIDE-001",
               "a particle starting inside a sphere is projected and absorbed once") {
  cInternalSphericalData* sphere = AMPS::CoreTests::MoverFixtureSphere();
  AMPS_EXPECT_TRUE(sphere != NULL);
  AMPS::CoreTests::Fixture::ValueScope<decltype(sphere->ParticleSphereInteraction)>
      callback_scope(sphere->ParticleSphereInteraction);
  AMPS::CoreTests::Fixture::ValueScope<int> count_scope(absorbing_count);
  sphere->ParticleSphereInteraction = AbsorbSphere;
  absorbing_count = 0;

  // The point is strictly inside the 0.2 m test sphere but is not its centre,
  // so the production recovery path has a well-defined outward projection.
  // The callback reports absorption; ownership of ParticleBuffer deletion
  // remains with the mover and is checked explicitly below.
  const double x[3] = {0.10, 0.0, 0.0};
  const double v[3] = {0.20, 0.0, 0.0};
  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node =
      AMPS::CoreTests::MoverFixtureNode(x);
  ParticleBufferScope particles(16);
  const long int particle = PIC::ParticleBuffer::GetNewParticle(true);
  PIC::ParticleBuffer::SetI(_H2O_SPEC_, particle);
  PIC::ParticleBuffer::SetX(const_cast<double*>(x), particle);
  PIC::ParticleBuffer::SetV(const_cast<double*>(v), particle);
  AMPS_EXPECT_EQ(PIC::Mover::UniformWeight_UniformTimeStep_noForce_TraceTrajectory(
                     particle, 0.5, node),
                 _PARTICLE_LEFT_THE_DOMAIN_);
  AMPS_EXPECT_EQ(absorbing_count, 1);
  AMPS_EXPECT_TRUE(!PIC::ParticleBuffer::IsParticleAllocated(particle));
}

AMPS_TEST_CASE("T4-BOUNDARY-MOVING-PLANE-001",
               "moving-plane event time, reflection, remaining step, and callback are exact") {
  AMPS::CoreTests::Fixture::ValueScope<PIC::Mover::fFindMovingSurfaceIntersection>
      finder_scope(PIC::Mover::FindMovingSurfaceIntersection);
  AMPS::CoreTests::Fixture::ValueScope<PIC::Mover::fProcessMovingSurfaceIntersection>
      interaction_scope(PIC::Mover::ProcessMovingSurfaceIntersection);
  AMPS::CoreTests::Fixture::ValueScope<int> finder_count_scope(moving_finder_count);
  AMPS::CoreTests::Fixture::ValueScope<int> interaction_count_scope(
      moving_interaction_count);
  PIC::Mover::FindMovingSurfaceIntersection = FindMovingPlane;
  PIC::Mover::ProcessMovingSurfaceIntersection = ReflectMovingPlane;
  moving_finder_count = 0;
  moving_interaction_count = 0;
  const double x[3] = {-0.6, 0.6, 0.6};
  const double v[3] = {1.0, 0.0, 0.0};
  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node =
      AMPS::CoreTests::MoverFixtureNode(x);
  ParticleBufferScope particles(16);
  const long int particle = PIC::ParticleBuffer::GetNewParticle(true);
  PIC::ParticleBuffer::SetI(_H2O_SPEC_, particle);
  PIC::ParticleBuffer::SetX(const_cast<double*>(x), particle);
  PIC::ParticleBuffer::SetV(const_cast<double*>(v), particle);
  AMPS_EXPECT_EQ(PIC::Mover::UniformWeight_UniformTimeStep_noForce(
                     particle, 1.0, node),
                 _PARTICLE_MOTION_FINISHED_);
  AMPS_EXPECT_EQ(moving_interaction_count, 1);
  AMPS_EXPECT_TRUE(moving_finder_count >= 2);
  AMPS_EXPECT_NEAR(PIC::ParticleBuffer::GetV(particle)[0], -0.6, 2.0e-14);
  // Hit is at t=0.75, x=0.15. The reflected particle travels -0.15 m
  // during the remaining 0.25 s, plus the documented outward mesh epsilon.
  AMPS_EXPECT_NEAR(PIC::ParticleBuffer::GetX(particle)[0], -1.0e-6, 2.0e-12);
  node = PIC::Mesh::mesh->findTreeNode(PIC::ParticleBuffer::GetX(particle));
  DetachMovedParticle(particle, node);
}

AMPS_TEST_CASE("T4-BOUNDARY-MOVING-SPHERE-001",
               "moving-sphere quadratic event and surface-frame reflection are exact") {
  AMPS::CoreTests::Fixture::ValueScope<PIC::Mover::fFindMovingSurfaceIntersection>
      finder_scope(PIC::Mover::FindMovingSurfaceIntersection);
  AMPS::CoreTests::Fixture::ValueScope<PIC::Mover::fProcessMovingSurfaceIntersection>
      interaction_scope(PIC::Mover::ProcessMovingSurfaceIntersection);
  AMPS::CoreTests::Fixture::ValueScope<int> finder_count_scope(moving_finder_count);
  AMPS::CoreTests::Fixture::ValueScope<int> interaction_count_scope(
      moving_interaction_count);
  PIC::Mover::FindMovingSurfaceIntersection = FindMovingSphere;
  PIC::Mover::ProcessMovingSurfaceIntersection = ReflectMovingPlane;
  moving_finder_count = 0;
  moving_interaction_count = 0;

  // The particle and sphere move along +x at 1.0 and 0.1 m/s. Their initial
  // leading-surface gap is 0.4 m, so t_hit=0.4/(1.0-0.1)=4/9 s. Specular
  // reflection in the sphere frame gives v'=2u-v=-0.8 m/s; the remaining 5/9
  // s returns the particle to x=-0.6 m, apart from the mover's face epsilon.
  const double x[3] = {-0.6, 0.6, 0.6};
  const double v[3] = {1.0, 0.0, 0.0};
  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node =
      AMPS::CoreTests::MoverFixtureNode(x);
  ParticleBufferScope particles(16);
  const long int particle = PIC::ParticleBuffer::GetNewParticle(true);
  PIC::ParticleBuffer::SetI(_H2O_SPEC_, particle);
  PIC::ParticleBuffer::SetX(const_cast<double*>(x), particle);
  PIC::ParticleBuffer::SetV(const_cast<double*>(v), particle);
  AMPS_EXPECT_EQ(PIC::Mover::UniformWeight_UniformTimeStep_noForce(
                     particle, 1.0, node),
                 _PARTICLE_MOTION_FINISHED_);
  AMPS_EXPECT_EQ(moving_interaction_count, 1);
  AMPS_EXPECT_TRUE(moving_finder_count >= 2);
  AMPS_EXPECT_NEAR(PIC::ParticleBuffer::GetV(particle)[0], -0.8, 2.0e-14);
  AMPS_EXPECT_NEAR(PIC::ParticleBuffer::GetX(particle)[0], -0.6, 2.0e-6);
  AMPS_EXPECT_NEAR(PIC::ParticleBuffer::GetX(particle)[1], 0.6, 2.0e-14);
  AMPS_EXPECT_NEAR(PIC::ParticleBuffer::GetX(particle)[2], 0.6, 2.0e-14);
  node = PIC::Mesh::mesh->findTreeNode(PIC::ParticleBuffer::GetX(particle));
  DetachMovedParticle(particle, node);
}

AMPS_TEST_CASE("T4-INJECTION-MAXWELLIAN-FLUX-001",
               "half-space Maxwellian flux matches independent zero/inward/outward laws") {
  const double density = 3.0e11;       // m^-3
  const double temperature = 420.0;    // K
  const double normal[3] = {1.0, 0.0, 0.0};
  const double thermal_sigma = std::sqrt(
      Kbol * temperature / PIC::MolecularData::GetMass(_H2O_SPEC_));
  const double normal_bulk[3] = {0.0, -1.5 * thermal_sigma,
                                 1.5 * thermal_sigma};
  for (int fixture = 0; fixture < 3; ++fixture) {
    const double bulk[3] = {normal_bulk[fixture], 0.2 * thermal_sigma, 0.0};
    const double observed = PIC::BC::CalculateInjectionRate_MaxwellianDistribution(
        density, temperature, bulk, const_cast<double*>(normal), _H2O_SPEC_);
    const double expected = IndependentInwardMaxwellianFlux(
        density, temperature, normal_bulk[fixture], _H2O_SPEC_);
    AMPS_EXPECT_TRUE(std::fabs(observed - expected) / expected < 3.0e-14);
  }

  rnd_seed(4401);
  const double bulk[3] = {-0.4 * thermal_sigma, 0.1 * thermal_sigma, 0.0};
  double mean_normal = 0.0;
  const int samples = 100000;
  for (int sample = 0; sample < samples; ++sample) {
    double velocity[3];
    PIC::Distribution::InjectMaxwellianDistribution(
        velocity, bulk, temperature, const_cast<double*>(normal), _H2O_SPEC_);
    AMPS_EXPECT_TRUE(velocity[0] <= 0.0);
    mean_normal += velocity[0];
  }
  AMPS_EXPECT_TRUE(mean_normal / samples < 0.0);
}

AMPS_TEST_CASE("T4-INJECTION-POISSON-001",
               "production exponential waits reproduce Poisson mean and variance") {
  AMPS_EXPECT_TRUE(std::isinf(PIC::BC::SamplePoissonWaitingTime(0.0)));
  AMPS_EXPECT_TRUE(std::isinf(PIC::BC::SamplePoissonWaitingTime(-1.0)));
  AMPS_EXPECT_TRUE(std::isinf(PIC::BC::SamplePoissonWaitingTime(
      std::numeric_limits<double>::quiet_NaN())));
  const double rate = 3.0;       // s^-1
  const double interval = 0.4;   // s
  const int windows = 100000;
  double sum = 0.0, sum2 = 0.0;
  rnd_seed(4402);
  for (int window = 0; window < windows; ++window) {
    double elapsed = 0.0;
    int count = 0;
    while ((elapsed += PIC::BC::SamplePoissonWaitingTime(rate)) < interval)
      ++count;
    sum += count;
    sum2 += count * count;
  }
  const double mean = sum / windows;
  const double variance = sum2 / windows - mean * mean;
  const double expected = rate * interval;
  AMPS_EXPECT_TRUE(std::fabs(mean - expected) < 0.015);
  AMPS_EXPECT_TRUE(std::fabs(variance - expected) < 0.025);
}

AMPS_TEST_CASE("T4-INJECTION-VOLUME-POSITION-001",
               "volume source positions are uniform inside the selected AMR cell") {
#if _PIC_VOLUME_PARTICLE_INJECTION_MODE_ == _PIC_VOLUME_PARTICLE_INJECTION_MODE__ON_
  const double probe[3] = {0.6, 0.6, 0.6};
  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node =
      AMPS::CoreTests::MoverFixtureNode(probe);
  const int i = 1, j = 1, k = 1;
  const double cell_min[3] = {
      node->xmin[0] + i * (node->xmax[0] - node->xmin[0]) / _BLOCK_CELLS_X_,
      node->xmin[1] + j * (node->xmax[1] - node->xmin[1]) / _BLOCK_CELLS_Y_,
      node->xmin[2] + k * (node->xmax[2] - node->xmin[2]) / _BLOCK_CELLS_Z_};
  const double cell_width[3] = {
      (node->xmax[0] - node->xmin[0]) / _BLOCK_CELLS_X_,
      (node->xmax[1] - node->xmin[1]) / _BLOCK_CELLS_Y_,
      (node->xmax[2] - node->xmin[2]) / _BLOCK_CELLS_Z_};
  double mean[3] = {0.0, 0.0, 0.0};
  double second[3] = {0.0, 0.0, 0.0};
  const int samples = 120000;
  rnd_seed(4403);
  for (int sample = 0; sample < samples; ++sample) {
    double position[3];
    PIC::VolumeParticleInjection::GetRandomCellPosition(position, i, j, k, node);
    for (int dimension = 0; dimension < 3; ++dimension) {
      const double unit = (position[dimension] - cell_min[dimension]) /
                          cell_width[dimension];
      AMPS_EXPECT_TRUE(unit >= 0.0 && unit < 1.0);
      mean[dimension] += unit;
      second[dimension] += unit * unit;
    }
  }
  for (int dimension = 0; dimension < 3; ++dimension) {
    mean[dimension] /= samples;
    second[dimension] /= samples;
    AMPS_EXPECT_TRUE(std::fabs(mean[dimension] - 0.5) < 0.0035);
    AMPS_EXPECT_TRUE(std::fabs(second[dimension] - 1.0 / 3.0) < 0.0035);
  }
#else
  AMPS_EXPECT_TRUE(false);
#endif
}

AMPS_TEST_CASE("T4-INJECTION-FACE-POSITION-001",
               "triangle-face source positions are uniform in physical area") {
  // A unit right triangle has x>=0, y>=0, x+y<=1.  Uniform area gives
  // E[x]=E[y]=1/3, E[x^2]=E[y^2]=1/6, and E[xy]=1/12.  These independent
  // simplex moments detect edge-uniform and biased-barycentric samplers.
  CutCell::cTriangleFace face;
  double a[3] = {0.0, 0.0, 0.0};
  double b[3] = {1.0, 0.0, 0.0};
  double c[3] = {0.0, 1.0, 0.0};
  face.SetFaceNodes(a, b, c);
  const int samples = 120000;
  double mean_x = 0.0, mean_y = 0.0;
  double second_x = 0.0, second_y = 0.0, cross = 0.0;
  rnd_seed(4404);
  for (int sample = 0; sample < samples; ++sample) {
    double position[3];
    face.GetRandomPosition(position);
    AMPS_EXPECT_TRUE(position[0] >= 0.0 && position[1] >= 0.0 &&
                     position[0] + position[1] <= 1.0 + 2.0e-15);
    AMPS_EXPECT_NEAR(position[2], 0.0, 0.0);
    mean_x += position[0];
    mean_y += position[1];
    second_x += position[0] * position[0];
    second_y += position[1] * position[1];
    cross += position[0] * position[1];
  }
  mean_x /= samples;
  mean_y /= samples;
  second_x /= samples;
  second_y /= samples;
  cross /= samples;
  AMPS_EXPECT_TRUE(std::fabs(mean_x - 1.0 / 3.0) < 0.0035);
  AMPS_EXPECT_TRUE(std::fabs(mean_y - 1.0 / 3.0) < 0.0035);
  AMPS_EXPECT_TRUE(std::fabs(second_x - 1.0 / 6.0) < 0.0035);
  AMPS_EXPECT_TRUE(std::fabs(second_y - 1.0 / 6.0) < 0.0035);
  AMPS_EXPECT_TRUE(std::fabs(cross - 1.0 / 12.0) < 0.0025);
}

AMPS_TEST_CASE("T4-INJECTION-SPHERE-VOLUME-001",
               "spherical source radii follow the uniform-volume r-cubed CDF") {
  AMPS_EXPECT_TRUE(std::isnan(PIC::BC::SampleUniformSphereVolumeRadius(0.0)));
  AMPS_EXPECT_TRUE(std::isnan(PIC::BC::SampleUniformSphereVolumeRadius(-1.0)));
  AMPS_EXPECT_TRUE(std::isnan(PIC::BC::SampleUniformSphereVolumeRadius(
      std::numeric_limits<double>::infinity())));

  // For X=r/R with CDF X^3, E[X^k]=3/(k+3).  Checking both the first and
  // second moments distinguishes the correct volume law from uniform radius
  // and uniform surface sampling without using production output as oracle.
  const double radius = 7.5;  // m; scale must cancel from normalized moments.
  const int samples = 120000;
  double mean = 0.0, second = 0.0;
  rnd_seed(4405);
  for (int sample = 0; sample < samples; ++sample) {
    const double normalized =
        PIC::BC::SampleUniformSphereVolumeRadius(radius) / radius;
    AMPS_EXPECT_TRUE(normalized >= 0.0 && normalized < 1.0);
    mean += normalized;
    second += normalized * normalized;
  }
  mean /= samples;
  second /= samples;
  AMPS_EXPECT_TRUE(std::fabs(mean - 3.0 / 4.0) < 0.0030);
  AMPS_EXPECT_TRUE(std::fabs(second - 3.0 / 5.0) < 0.0030);
}

AMPS_TEST_CASE("T4-SAMPLING-MOMENTS-001",
               "cell sampling recovers density, flux, pressure, energy, and temperature") {
  const double position[3] = {0.6, 0.6, 0.6};
  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node =
      AMPS::CoreTests::MoverFixtureNode(position);
  int i, j, k;
  const int local = PIC::Mesh::mesh->FindCellIndex(
      const_cast<double*>(position), i, j, k, node, false);
  AMPS_EXPECT_TRUE(local >= 0);
  PIC::Mesh::cDataCenterNode* cell = node->block->GetCenterNode(local);
  long int& first = node->block->FirstCellParticleTable[
      i + _BLOCK_CELLS_X_ * (j + _BLOCK_CELLS_Y_ * k)];
  ParticleBufferScope particles(32);
  EmptyListHeadScope list_scope(first);
  MeshCellStateScope state_scope(node->block, cell);
  const double weight = 2.5e7;
  node->block->SetLocalParticleWeight(weight, _H2O_SPEC_);
  node->block->SetLocalTimeStep(1.0, _H2O_SPEC_);
  const double velocity[3][3] = {
      {1.0, 2.0, 3.0}, {3.0, -1.0, 2.0}, {-2.0, 2.0, 1.0}};
  for (int particle_index = 0; particle_index < 3; ++particle_index) {
    const long int particle = PIC::ParticleBuffer::GetNewParticle(first);
    PIC::ParticleBuffer::SetI(_H2O_SPEC_, particle);
    PIC::ParticleBuffer::SetX(const_cast<double*>(position), particle);
    PIC::ParticleBuffer::SetV(const_cast<double*>(velocity[particle_index]),
                              particle);
  }
  int* counts[1];
  counts[0] = new int[PIC::nTotalSpecies];
  std::fill(counts[0], counts[0] + PIC::nTotalSpecies, 0);
  PIC::Mesh::flushCollectingSamplingBuffer(cell);
  PIC::Sampling::ProcessCell(i, j, k, counts, node, 0);
  delete[] counts[0];
  PIC::Mesh::switchSamplingBuffers();
  AMPS::CoreTests::Fixture::ValueScope<long int> sample_length_scope(
      PIC::LastSampleLength);
  PIC::LastSampleLength = 1;
  const double expected_mean[3] = {2.0 / 3.0, 1.0, 2.0};
  double observed_mean[3];
  cell->GetBulkVelocity(observed_mean, _H2O_SPEC_);
  for (int dimension = 0; dimension < 3; ++dimension)
    AMPS_EXPECT_NEAR(observed_mean[dimension], expected_mean[dimension],
                     2.0e-14);
  AMPS_EXPECT_NEAR(cell->GetNumberDensity(_H2O_SPEC_),
                   3.0 * weight / cell->Measure, 1.0e-12 *
                   (3.0 * weight / cell->Measure));
  const double number_density = 3.0 * weight / cell->Measure;
  AMPS_EXPECT_NEAR(cell->GetDatumAverage(
                       PIC::Mesh::DatumParticleNumber, _H2O_SPEC_),
                   3.0, 0.0);
  AMPS_EXPECT_NEAR(cell->GetCompleteSampleCellParticleWeight(_H2O_SPEC_),
                   3.0 * weight, 0.0);

  double observed_velocity2[3], observed_cross[3];
  cell->GetDatumAverage(PIC::Mesh::DatumParticleVelocity2,
                        observed_velocity2, _H2O_SPEC_);
  cell->GetDatumAverage(PIC::Mesh::DatumParticleVelocity2Tensor,
                        observed_cross, _H2O_SPEC_);
  const double expected_velocity2[3] = {14.0 / 3.0, 3.0, 14.0 / 3.0};
  const double expected_cross[3] = {-5.0 / 3.0, 2.0, 7.0 / 3.0};
  for (int dimension = 0; dimension < 3; ++dimension) {
    AMPS_EXPECT_NEAR(observed_velocity2[dimension],
                     expected_velocity2[dimension], 2.0e-14);
    AMPS_EXPECT_NEAR(observed_cross[dimension],
                     expected_cross[dimension], 2.0e-14);
  }
  const double expected_mean_speed =
      (2.0 * std::sqrt(14.0) + 3.0) / 3.0;
  AMPS_EXPECT_NEAR(cell->GetMeanParticleSpeed(_H2O_SPEC_),
                   expected_mean_speed, 2.0e-14);

  // The raw velocity moments define the physical observables without a
  // production-derived oracle: particle flux is n*u, kinetic-energy density
  // is (m*n/2)<v^2>, and pressure is
  // P_ij=m*n*(<v_i v_j>-u_i*u_j), all in SI units. Exact raw moments above
  // therefore validate flux [m^-2 s^-1], energy [J m^-3], and all six
  // independent components of the symmetric pressure tensor [Pa].
  double particle_flux[3], pressure_diagonal[3], pressure_cross[3];
  const double mass = PIC::MolecularData::GetMass(_H2O_SPEC_);
  for (int dimension = 0; dimension < 3; ++dimension) {
    particle_flux[dimension] = number_density * observed_mean[dimension];
    pressure_diagonal[dimension] = mass * number_density *
        (observed_velocity2[dimension] -
         observed_mean[dimension] * observed_mean[dimension]);
  }
  pressure_cross[0] = mass * number_density *
      (observed_cross[0] - observed_mean[0] * observed_mean[1]);
  pressure_cross[1] = mass * number_density *
      (observed_cross[1] - observed_mean[1] * observed_mean[2]);
  pressure_cross[2] = mass * number_density *
      (observed_cross[2] - observed_mean[2] * observed_mean[0]);
  const double kinetic_energy_density = 0.5 * mass * number_density *
      (observed_velocity2[0] + observed_velocity2[1] +
       observed_velocity2[2]);
  AMPS_EXPECT_NEAR(particle_flux[0], number_density * 2.0 / 3.0,
                   2.0e-14 * number_density);
  AMPS_EXPECT_NEAR(particle_flux[1], number_density,
                   2.0e-14 * number_density);
  AMPS_EXPECT_NEAR(particle_flux[2], 2.0 * number_density,
                   4.0e-14 * number_density);
  AMPS_EXPECT_TRUE(kinetic_energy_density > 0.0);
  for (int component = 0; component < 3; ++component) {
    AMPS_EXPECT_TRUE(std::isfinite(pressure_diagonal[component]) &&
                     pressure_diagonal[component] >= 0.0);
    AMPS_EXPECT_TRUE(std::isfinite(pressure_cross[component]));
  }
  double variance_sum = 0.0;
  for (int particle_index = 0; particle_index < 3; ++particle_index)
    for (int dimension = 0; dimension < 3; ++dimension)
      variance_sum += std::pow(velocity[particle_index][dimension] -
                               expected_mean[dimension], 2);
  const double expected_temperature =
      PIC::MolecularData::GetMass(_H2O_SPEC_) * variance_sum /
      (9.0 * Kbol);
  double observed_temperature = 0.0;
  cell->GetTranslationalTemperature(&observed_temperature, _H2O_SPEC_);
  AMPS_EXPECT_TRUE(std::fabs(observed_temperature - expected_temperature) /
                       expected_temperature < 2.0e-14);
  PIC::Mesh::flushCompletedSamplingBuffer(cell);
  PIC::Mesh::switchSamplingBuffers();
  PIC::Mesh::flushCollectingSamplingBuffer(cell);
}

AMPS_TEST_CASE("T4-SAMPLING-WEIGHT-RESCALING-001",
               "legal particle-count/weight rescaling preserves sampled moments") {
  const double position[3] = {0.6, 0.6, 0.6};
  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node =
      AMPS::CoreTests::MoverFixtureNode(position);
  int i, j, k;
  const int local = PIC::Mesh::mesh->FindCellIndex(
      const_cast<double*>(position), i, j, k, node, false);
  PIC::Mesh::cDataCenterNode* cell = node->block->GetCenterNode(local);
  long int& first = node->block->FirstCellParticleTable[
      i + _BLOCK_CELLS_X_ * (j + _BLOCK_CELLS_Y_ * k)];
  ParticleBufferScope particles(32);
  EmptyListHeadScope list_scope(first);
  MeshCellStateScope state_scope(node->block, cell);
  AMPS::CoreTests::Fixture::ValueScope<long int> sample_length_scope(
      PIC::LastSampleLength);
  PIC::LastSampleLength = 1;
  PIC::Mesh::flushCollectingSamplingBuffer(cell);
  PIC::Mesh::flushCompletedSamplingBuffer(cell);

  const double base_weight = 4.0e6;
  const double states[2][3] = {{-2.0, 1.0, 0.5}, {4.0, -3.0, 2.5}};
  double reference_density = 0.0, reference_mean[3] = {0.0, 0.0, 0.0};
  double reference_temperature = 0.0, reference_speed = 0.0;
  for (int representation = 0; representation < 2; ++representation) {
    // Species-global setters are monotone minimum reducers because different
    // blocks may independently report their stability limit. Exercise the
    // equivalent populations in decreasing-W order: two records at 2W, then
    // four records at W. Raising W in-place would violate the public setter
    // contract and test fixture behavior rather than physical rescaling.
    const int copies = representation == 0 ? 1 : 2;
    const double weight = representation == 0 ? 2.0 * base_weight
                                               : base_weight;
    node->block->SetLocalParticleWeight(weight, _H2O_SPEC_);
    node->block->SetLocalTimeStep(1.0, _H2O_SPEC_);
    std::vector<long int> records;
    for (int state = 0; state < 2; ++state) {
      for (int copy = 0; copy < copies; ++copy) {
        const long int particle = PIC::ParticleBuffer::GetNewParticle(first);
        PIC::ParticleBuffer::SetI(_H2O_SPEC_, particle);
        PIC::ParticleBuffer::SetX(const_cast<double*>(position), particle);
        PIC::ParticleBuffer::SetV(
            const_cast<double*>(states[state]), particle);
        records.push_back(particle);
      }
    }
    int* counts[1];
    counts[0] = new int[PIC::nTotalSpecies];
    std::fill(counts[0], counts[0] + PIC::nTotalSpecies, 0);
    PIC::Sampling::ProcessCell(i, j, k, counts, node, 0);
    delete[] counts[0];
    PIC::Mesh::switchSamplingBuffers();

    const double density = cell->GetNumberDensity(_H2O_SPEC_);
    double mean[3], temperature = 0.0;
    cell->GetBulkVelocity(mean, _H2O_SPEC_);
    cell->GetTranslationalTemperature(&temperature, _H2O_SPEC_);
    const double mean_speed = cell->GetMeanParticleSpeed(_H2O_SPEC_);
    if (representation == 0) {
      reference_density = density;
      std::copy(mean, mean + 3, reference_mean);
      reference_temperature = temperature;
      reference_speed = mean_speed;
    } else {
      AMPS_EXPECT_NEAR(density, reference_density,
                       2.0e-14 * reference_density);
      for (int dimension = 0; dimension < 3; ++dimension)
        AMPS_EXPECT_NEAR(mean[dimension], reference_mean[dimension], 2.0e-14);
      AMPS_EXPECT_NEAR(temperature, reference_temperature,
                       2.0e-14 * reference_temperature);
      AMPS_EXPECT_NEAR(mean_speed, reference_speed, 2.0e-14);
    }

    // Remove the current representation before constructing the rescaled
    // one. Completed data are consumed and cleared, then a buffer swap puts a
    // clean collecting region back in service without stale carry-over.
    for (std::size_t record = 0; record < records.size(); ++record)
      PIC::ParticleBuffer::DeleteParticle(records[record], first);
    AMPS_EXPECT_EQ(first, -1);
    PIC::Mesh::flushCompletedSamplingBuffer(cell);
    PIC::Mesh::switchSamplingBuffers();
    PIC::Mesh::flushCollectingSamplingBuffer(cell);
  }
}

AMPS_TEST_CASE("T4-SAMPLING-BUFFER-LIFECYCLE-001",
               "sampling collect/swap/flush exposes no stale carry-over") {
  const double position[3] = {0.6, 0.6, 0.6};
  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node =
      AMPS::CoreTests::MoverFixtureNode(position);
  int i, j, k;
  const int local = PIC::Mesh::mesh->FindCellIndex(
      const_cast<double*>(position), i, j, k, node, false);
  PIC::Mesh::cDataCenterNode* cell = node->block->GetCenterNode(local);
  MeshCellStateScope state_scope(node->block, cell);
  AMPS_EXPECT_TRUE(PIC::Mesh::sampleSetDataLength >=
                   static_cast<int>(sizeof(double)));
  double* collecting = reinterpret_cast<double*>(
      cell->GetAssociatedDataBufferPointer() +
      PIC::Mesh::collectingCellSampleDataPointerOffset);
  double* completed = reinterpret_cast<double*>(
      cell->GetAssociatedDataBufferPointer() +
      PIC::Mesh::completedCellSampleDataPointerOffset);
  collecting[0] = 17.0;
  completed[0] = 3.0;
  PIC::Mesh::switchSamplingBuffers();
  completed = reinterpret_cast<double*>(cell->GetAssociatedDataBufferPointer() +
      PIC::Mesh::completedCellSampleDataPointerOffset);
  collecting = reinterpret_cast<double*>(cell->GetAssociatedDataBufferPointer() +
      PIC::Mesh::collectingCellSampleDataPointerOffset);
  AMPS_EXPECT_NEAR(completed[0], 17.0, 0.0);
  AMPS_EXPECT_NEAR(collecting[0], 3.0, 0.0);
  PIC::Mesh::flushCollectingSamplingBuffer(cell);
  AMPS_EXPECT_NEAR(collecting[0], 0.0, 0.0);
  PIC::Mesh::flushCompletedSamplingBuffer(cell);
  AMPS_EXPECT_NEAR(completed[0], 0.0, 0.0);
  // Restore the global offset orientation for test-order independence.
  PIC::Mesh::switchSamplingBuffers();
}

#if _PIC_BC__PERIODIC_MODE_ == _PIC_BC__PERIODIC_MODE_ON_
namespace {

double PeriodicFixtureResolution(double*) {
  // The AMR resolution callback is a requested physical cell size [m].  The
  // legacy periodic manager requires its added ghost shell and user domain to
  // land on the same dyadic refinement level.  L/16 satisfies that invariant
  // for the configured two-cell blocks and retains a bounded one-rank mesh.
  return 0.125;
}

}  // namespace

AMPS_TEST_CASE("T4-BOUNDARY-PERIODIC-001",
               "periodic AMR pairs wrap particles and synchronize node data") {
  // This case owns the process mesh.  The periodic profile selects it alone,
  // so no earlier lazy mover fixture can have initialized a nonperiodic root.
  AMPS_EXPECT_TRUE(PIC::Mesh::mesh->rootTree == NULL);
  double xmin[3] = {-1.0, -1.0, -1.0};
  double xmax[3] = {1.0, 1.0, 1.0};
  PIC::Mesh::mesh->AllowBlockAllocation = false;
  PIC::BC::ExternalBoundary::Periodic::Init(
      xmin, xmax, PeriodicFixtureResolution);
  PIC::Mesh::mesh->buildMesh();
  PIC::Mesh::mesh->CreateNewParallelDistributionLists();
  PIC::Mesh::mesh->AllowBlockAllocation = true;
  PIC::Mesh::mesh->AllocateTreeBlocks();
  PIC::Mesh::mesh->InitCellMeasure();
  PIC::BC::ExternalBoundary::Periodic::InitBlockPairTable(true);

  for (int dimension = 0; dimension < 3; ++dimension) {
    AMPS_EXPECT_NEAR(PIC::BC::ExternalBoundary::Periodic::L[dimension],
                     2.0, 0.0);
    AMPS_EXPECT_TRUE(PIC::BC::ExternalBoundary::Periodic::BoundaryDx[dimension]
                     > 0.0);
  }
  AMPS_EXPECT_TRUE(
      PIC::BC::ExternalBoundary::Periodic::BlockPairTableLength > 0);

  PIC::BC::ExternalBoundary::Periodic::cBlockPairTable* selected = NULL;
  for (int pair = 0;
       pair < PIC::BC::ExternalBoundary::Periodic::BlockPairTableLength;
       ++pair) {
    PIC::BC::ExternalBoundary::Periodic::cBlockPairTable& candidate =
        PIC::BC::ExternalBoundary::Periodic::BlockPairTable[pair];
    AMPS_EXPECT_TRUE(candidate.GhostBlock != NULL &&
                     candidate.RealBlock != NULL);
    if (candidate.GhostBlock->Thread == PIC::ThisThread &&
        candidate.RealBlock->Thread == PIC::ThisThread &&
        candidate.GhostBlock->block != NULL &&
        candidate.RealBlock->block != NULL) {
      selected = &candidate;
      break;
    }
  }
  AMPS_EXPECT_TRUE(selected != NULL);

  const int local_cell = 0;
  long int& ghost_head =
      selected->GhostBlock->block->FirstCellParticleTable[local_cell];
  long int& real_head =
      selected->RealBlock->block->FirstCellParticleTable[local_cell];
  ParticleBufferScope particles(32);
  EmptyListHeadScope ghost_list_scope(ghost_head);
  EmptyListHeadScope real_list_scope(real_head);
  const long int particle = PIC::ParticleBuffer::GetNewParticle(ghost_head);
  PIC::ParticleBuffer::SetI(_H2O_SPEC_, particle);
  double before[3];
  for (int dimension = 0; dimension < 3; ++dimension)
    before[dimension] = 0.5 * (selected->GhostBlock->xmin[dimension] +
                               selected->GhostBlock->xmax[dimension]);
  PIC::ParticleBuffer::SetX(before, particle);
  double velocity[3] = {1.0, -2.0, 3.0};
  PIC::ParticleBuffer::SetV(velocity, particle);

  PIC::BC::ExternalBoundary::Periodic::ExchangeParticlesLocal(
      selected->RealBlock, selected->GhostBlock);
  AMPS_EXPECT_EQ(ghost_head, -1);
  AMPS_EXPECT_EQ(real_head, particle);
  const double* after = PIC::ParticleBuffer::GetX(particle);
  for (int dimension = 0; dimension < 3; ++dimension) {
    const double expected = before[dimension] +
        selected->RealBlock->xmin[dimension] -
        selected->GhostBlock->xmin[dimension];
    AMPS_EXPECT_NEAR(after[dimension], expected, 2.0e-14);
    AMPS_EXPECT_TRUE(after[dimension] >= selected->RealBlock->xmin[dimension]);
    AMPS_EXPECT_TRUE(after[dimension] < selected->RealBlock->xmax[dimension]);
    AMPS_EXPECT_NEAR(PIC::ParticleBuffer::GetV(particle)[dimension],
                     velocity[dimension], 0.0);
  }

  // Exercise the same local block-data copy selected by periodic halo
  // synchronization.  Sampling storage is opaque here: only byte-for-byte
  // direction and completeness are asserted, independent of PIC field layout.
  PIC::Mesh::cDataCenterNode* real_center =
      selected->RealBlock->block->GetCenterNode(
          _getCenterNodeLocalNumber(0, 0, 0));
  PIC::Mesh::cDataCenterNode* ghost_center =
      selected->GhostBlock->block->GetCenterNode(
          _getCenterNodeLocalNumber(0, 0, 0));
  AMPS_EXPECT_TRUE(real_center != NULL && ghost_center != NULL);
  AMPS_EXPECT_TRUE(PIC::Mesh::cDataCenterNode_static_data::
                       totalAssociatedDataLength >=
                   static_cast<int>(sizeof(double)));
  const double marker = 6.125;
  const double stale = -9.5;
  std::memcpy(real_center->GetAssociatedDataBufferPointer(), &marker,
              sizeof(marker));
  std::memcpy(ghost_center->GetAssociatedDataBufferPointer(), &stale,
              sizeof(stale));
  PIC::BC::ExternalBoundary::Periodic::ExchangeBlockDataLocal(
      *selected, NULL, NULL);
  double received = 0.0;
  std::memcpy(&received, ghost_center->GetAssociatedDataBufferPointer(),
              sizeof(received));
  AMPS_EXPECT_NEAR(received, marker, 0.0);

  // EmptyListHeadScope restores both list heads before ParticleBufferScope
  // releases the record, preventing a dangling particle index in the mesh.
}
#endif

