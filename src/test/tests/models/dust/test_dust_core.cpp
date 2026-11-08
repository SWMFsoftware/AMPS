#include "test_harness.h"

#include "DustCore.h"
#include "DustSplitMerge.h"

#include <cmath>
#include <limits>
#include <stdexcept>

// This table header is also included by Dust.h inside the production Mie
// namespace.  Include it here under the same namespace so the focused test
// calls the exact production lookup without pulling in PIC or MPI.
namespace ElectricallyChargedDust {
namespace ScatteringEfficientcy {
namespace Mie {
#include "Dust_MieScatteringEfficiency_Fink2012Icarus.h"
}  // namespace Mie
}  // namespace ScatteringEfficientcy
}  // namespace ElectricallyChargedDust

namespace {
template <typename Function>
void ExpectInvalid(Function function) {
  bool rejected = false;
  try {
    function();
  } catch (const std::invalid_argument&) {
    rejected = true;
  }
  AMPS_EXPECT_TRUE(rejected);
}

double RadiusForSizeParameterLog10(double logarithm, double wavelength) {
  const double pi = 3.141592653589793238462643383279502884;
  return wavelength * std::pow(10.0, logarithm) / (2.0 * pi);
}

ElectricallyChargedDust::SplitMerge::Particle MakeDustParticle(
    double weight, double radius, double grain_mass, double grain_charge,
    double vx, double vy = 0.0, double vz = 0.0) {
  ElectricallyChargedDust::SplitMerge::Particle particle;
  particle.compiled_species = 7;
  particle.ownership_cell = 42;
  particle.individual_stat_weight_correction = weight;
  particle.radius = radius;
  particle.grain_mass = grain_mass;
  particle.grain_charge = grain_charge;
  particle.velocity[0] = vx;
  particle.velocity[1] = vy;
  particle.velocity[2] = vz;
  return particle;
}

void ExpectBudgetsNear(
    const ElectricallyChargedDust::SplitMerge::Budget& left,
    const ElectricallyChargedDust::SplitMerge::Budget& right,
    double relative_tolerance) {
  const double number_scale = std::fmax(1.0, std::fabs(left.represented_number));
  const double mass_scale = std::fmax(1.0e-300, std::fabs(left.mass));
  const double charge_scale = std::fmax(1.0e-300, std::fabs(left.charge));
  const double energy_scale = std::fmax(1.0e-300, std::fabs(left.kinetic_energy));
  AMPS_EXPECT_NEAR(right.represented_number, left.represented_number,
                   relative_tolerance * number_scale);
  AMPS_EXPECT_NEAR(right.mass, left.mass, relative_tolerance * mass_scale);
  AMPS_EXPECT_NEAR(right.charge, left.charge,
                   relative_tolerance * charge_scale);
  AMPS_EXPECT_NEAR(right.kinetic_energy, left.kinetic_energy,
                   relative_tolerance * energy_scale);
  for (int component = 0; component != 3; ++component) {
    const double scale = std::fmax(1.0e-300,
        std::fabs(left.momentum[component]));
    AMPS_EXPECT_NEAR(right.momentum[component], left.momentum[component],
                     relative_tolerance * scale);
  }
}
}  // namespace

AMPS_TEST_CASE("CORE-DUST-SIZE-BINS-001",
               "dust radii and speeds map to bounded logarithmic bins") {
  using ElectricallyChargedDust::Core::LogarithmicBin;
  // Four bins span four decades: their interiors are represented by 3e-N,
  // avoiding a test oracle that repeats the production logarithm/index code.
  AMPS_EXPECT_EQ(LogarithmicBin(0.0, 1.0e-7, 1.0e-3, 4), 0);
  AMPS_EXPECT_EQ(LogarithmicBin(3.0e-7, 1.0e-7, 1.0e-3, 4), 0);
  AMPS_EXPECT_EQ(LogarithmicBin(3.0e-6, 1.0e-7, 1.0e-3, 4), 1);
  AMPS_EXPECT_EQ(LogarithmicBin(3.0e-5, 1.0e-7, 1.0e-3, 4), 2);
  AMPS_EXPECT_EQ(LogarithmicBin(3.0e-4, 1.0e-7, 1.0e-3, 4), 3);
  AMPS_EXPECT_EQ(LogarithmicBin(1.0e-3, 1.0e-7, 1.0e-3, 4), 3);
  AMPS_EXPECT_EQ(LogarithmicBin(1.0, 1.0e-7, 1.0e-3, 4), 3);

  ExpectInvalid([]() {
    (void)ElectricallyChargedDust::Core::LogarithmicBin(
        -1.0, 1.0e-7, 1.0e-3, 4);
  });
  ExpectInvalid([]() {
    (void)ElectricallyChargedDust::Core::LogarithmicBin(
        1.0e-6, 1.0e-7, 1.0e-3, 0);
  });
}

AMPS_TEST_CASE("CORE-DUST-SIZE-DISTRIBUTION-001",
               "log-sampled power law has analytic normalization and p=1 limit") {
  using namespace ElectricallyChargedDust::Core;
  // For p=2 on [1,100], integral a^-1 d(ln a) = 1-1/100 = 0.99.
  AMPS_EXPECT_NEAR(PowerLawLogNormalization(1.0, 100.0, 2.0),
                   0.99, 2.0e-15);
  AMPS_EXPECT_NEAR(PowerLawLogDensity(10.0, 1.0, 100.0, 2.0),
                   0.1 / 0.99, 2.0e-15);
  // The removable p=1 singularity is exactly integral 1 d(ln a)=ln(100).
  const double logarithmic_width = std::log(100.0);
  AMPS_EXPECT_NEAR(PowerLawLogNormalization(1.0, 100.0, 1.0),
                   logarithmic_width, 2.0e-15);
  AMPS_EXPECT_NEAR(PowerLawLogDensity(10.0, 1.0, 100.0, 1.0),
                   1.0 / logarithmic_width, 2.0e-15);
  AMPS_EXPECT_NEAR(PowerLawLogMaximum(1.0, 100.0, 0.0), 100.0, 0.0);
  AMPS_EXPECT_NEAR(PowerLawLogMaximum(1.0, 100.0, 2.0), 1.0, 0.0);
}

AMPS_TEST_CASE("CORE-DUST-CHARGING-POTENTIAL-001",
               "isolated-grain charge and radius give the SI capacitance potential") {
  using ElectricallyChargedDust::Core::IsolatedSpherePotential;
  // Independent literal oracle: Coulomb's constant 8.98755179226117e9 SI,
  // q=1e-15 C and a=1e-6 m give phi=k*q/a=8.98755179226117 V.
  const double epsilon0 = 8.8541878128e-12;
  const double potential = IsolatedSpherePotential(1.0e-6, 1.0e-15, epsilon0);
  AMPS_EXPECT_NEAR(potential, 8.98755179226117, 2.0e-10);
  AMPS_EXPECT_NEAR(IsolatedSpherePotential(
      2.0e-6, -1.0e-15, epsilon0), -0.5 * potential, 2.0e-14);
  AMPS_EXPECT_NEAR(IsolatedSpherePotential(
      1.0e-6, 0.0, epsilon0), 0.0, 0.0);
  ExpectInvalid([epsilon0]() {
    (void)IsolatedSpherePotential(0.0, 1.0e-15, epsilon0);
  });
}

AMPS_TEST_CASE("CORE-DUST-FORCE-LIMITS-001",
               "Lorentz acceleration reproduces electric, magnetic, and zero limits") {
  using ElectricallyChargedDust::Core::LorentzAcceleration;
  const double electric[3] = {2.0, 4.0, 6.0};
  const double zero[3] = {0.0, 0.0, 0.0};
  const double velocity[3] = {3.0, 4.0, 0.0};
  double acceleration[3];

  LorentzAcceleration(2.0, 4.0, electric, zero, velocity, acceleration);
  AMPS_EXPECT_NEAR(acceleration[0], 1.0, 0.0);
  AMPS_EXPECT_NEAR(acceleration[1], 2.0, 0.0);
  AMPS_EXPECT_NEAR(acceleration[2], 3.0, 0.0);

  // v=(3,4,0) and B=(0,0,2) give v cross B=(8,-6,0); q/m=1/2.
  const double magnetic[3] = {0.0, 0.0, 2.0};
  LorentzAcceleration(2.0, 4.0, zero, magnetic, velocity, acceleration);
  AMPS_EXPECT_NEAR(acceleration[0], 4.0, 0.0);
  AMPS_EXPECT_NEAR(acceleration[1], -3.0, 0.0);
  AMPS_EXPECT_NEAR(acceleration[2], 0.0, 0.0);

  LorentzAcceleration(0.0, 4.0, electric, magnetic, velocity, acceleration);
  AMPS_EXPECT_NEAR(acceleration[0], 0.0, 0.0);
  AMPS_EXPECT_NEAR(acceleration[1], 0.0, 0.0);
  AMPS_EXPECT_NEAR(acceleration[2], 0.0, 0.0);
  ExpectInvalid([&]() {
    LorentzAcceleration(1.0, 0.0, electric, magnetic, velocity, acceleration);
  });
}

AMPS_TEST_CASE("CORE-DUST-INJECTION-LAW-001",
               "spherical-source direction and grain mass obey analytic laws") {
  using namespace ElectricallyChargedDust::Core;
  double direction[3] = {3.0, 4.0, 0.0};
  NormalizeInjectionDirection(direction);
  AMPS_EXPECT_NEAR(direction[0], 0.6, 1.0e-15);
  AMPS_EXPECT_NEAR(direction[1], 0.8, 1.0e-15);
  AMPS_EXPECT_NEAR(direction[2], 0.0, 0.0);
  AMPS_EXPECT_NEAR(direction[0] * direction[0] +
                   direction[1] * direction[1] +
                   direction[2] * direction[2], 1.0, 2.0e-15);

  // A nonzero centre is a retained literal regression for the production
  // spherical injector: x=x0+R*n and v=u*n, not an origin-centred shortcut.
  const double center[3] = {10.0, -2.0, 7.0};
  double position[3];
  double velocity[3];
  SphericalSourceState(center, 5.0, 20.0, direction, position, velocity);
  AMPS_EXPECT_NEAR(position[0], 13.0, 1.0e-15);
  AMPS_EXPECT_NEAR(position[1], 2.0, 1.0e-15);
  AMPS_EXPECT_NEAR(position[2], 7.0, 1.0e-15);
  AMPS_EXPECT_NEAR(velocity[0], 12.0, 1.0e-15);
  AMPS_EXPECT_NEAR(velocity[1], 16.0, 1.0e-15);
  AMPS_EXPECT_NEAR(velocity[2], 0.0, 0.0);

  // Independent literal: (4*pi/3)*1000*(1e-6)^3 kg.
  AMPS_EXPECT_NEAR(SphericalGrainMass(1.0e-6, 1000.0),
                   4.1887902047863905e-15, 2.0e-30);
  ExpectInvalid([]() {
    double zero_direction[3] = {0.0, 0.0, 0.0};
    ElectricallyChargedDust::Core::NormalizeInjectionDirection(zero_direction);
  });
  ExpectInvalid([&]() {
    const double nonunit[3] = {1.0, 1.0, 0.0};
    ElectricallyChargedDust::Core::SphericalSourceState(
        center, 5.0, 20.0, nonunit, position, velocity);
  });
}

AMPS_TEST_CASE("CORE-DUST-MIE-INTERPOLATION-001",
               "Fink-2012 table clamps and interpolates in log size parameter") {
  namespace Fink = ElectricallyChargedDust::ScatteringEfficientcy::Mie::
      Fink2012Icarus;
  const double wavelength = 1.0;
  AMPS_EXPECT_NEAR(Fink::Get(RadiusForSizeParameterLog10(-2.0, wavelength),
                            wavelength, Fink::n20i0001),
                   0.013, 0.0);
  AMPS_EXPECT_NEAR(Fink::Get(RadiusForSizeParameterLog10(4.0, wavelength),
                            wavelength, Fink::n20i04),
                   1.196, 0.0);

  // Table point 12 is aLog=-1+12*0.08=-0.04 and is a retained literal,
  // independent of the production interpolation arithmetic.
  AMPS_EXPECT_NEAR(Fink::Get(RadiusForSizeParameterLog10(-0.04, wavelength),
                            wavelength, Fink::n20i01),
                   0.836, 2.0e-14);
  // Halfway in log space between table points 10 (0.180) and 11 (0.495).
  AMPS_EXPECT_NEAR(Fink::Get(RadiusForSizeParameterLog10(-0.16, wavelength),
                            wavelength, Fink::n20i0001),
                   0.3375, 2.0e-14);
  ExpectInvalid([&]() {
    (void)Fink::Get(1.0, wavelength, Fink::nCases);
  });
  ExpectInvalid([&]() {
    (void)Fink::Get(1.0, 0.0, Fink::n20i0001);
  });
}

AMPS_TEST_CASE("CORE-DUST-MOMENTS-001",
               "weighted dust count velocity charge and current moments close") {
  using namespace ElectricallyChargedDust::Core;
  WeightedMoments moments;
  const double first_velocity[3] = {3.0, 4.0, 0.0};
  const double second_velocity[3] = {-1.0, 0.0, 0.0};
  AccumulateWeightedMoments(2.0, first_velocity, -2.0, moments);
  AccumulateWeightedMoments(3.0, second_velocity, 1.0, moments);

  // Literal sums: N=2+3; sum(w*v)=(3,8,0); sum(w*|v|)=10+3;
  // sum(w*q)=-4+3; current=(-12,-16,0)+(-3,0,0).  The absolute-current
  // budget deliberately exercises both signs: 2*|-2|*5+3*|1|*1=23.
  AMPS_EXPECT_NEAR(moments.number, 5.0, 0.0);
  AMPS_EXPECT_NEAR(moments.velocity[0], 3.0, 0.0);
  AMPS_EXPECT_NEAR(moments.velocity[1], 8.0, 0.0);
  AMPS_EXPECT_NEAR(moments.velocity[2], 0.0, 0.0);
  AMPS_EXPECT_NEAR(moments.speed, 13.0, 0.0);
  AMPS_EXPECT_NEAR(moments.charge, -1.0, 0.0);
  AMPS_EXPECT_NEAR(moments.current[0], -15.0, 0.0);
  AMPS_EXPECT_NEAR(moments.current[1], -16.0, 0.0);
  AMPS_EXPECT_NEAR(moments.current[2], 0.0, 0.0);
  AMPS_EXPECT_NEAR(moments.absolute_current, 23.0, 0.0);

  ExpectInvalid([&]() {
    AccumulateWeightedMoments(-1.0, first_velocity, 0.0, moments);
  });
}

AMPS_TEST_CASE("CORE-DUST-SPLIT-MERGE-001",
               "dust split/merge preserves budgets inside one intrinsic size group") {
  using namespace ElectricallyChargedDust::SplitMerge;
  const SizeGroups groups(1.0e-7, 1.0e-3, 4);
  const double radius = 2.0e-6;
  const double mass = 3.2e-15;
  const double charge = -4.0e-18;

  // Splitting copies the complete grain state and divides only AMPS's
  // individual statistical-weight correction.  This makes every represented
  // budget and sampled moment invariant without depending on a global weight.
  const Particle parent = MakeDustParticle(
      6.0, radius, mass, charge, 3.0, 4.0, 0.0);
  const std::vector<Particle> parent_vector(1, parent);
  const Budget parent_budget = EvaluateBudget(parent_vector);
  const Result split = Split(parent, 3, groups);
  AMPS_EXPECT_TRUE(split.ok());
  AMPS_EXPECT_EQ(split.particles.size(), static_cast<std::size_t>(3));
  for (std::size_t index = 0; index != split.particles.size(); ++index) {
    AMPS_EXPECT_NEAR(split.particles[index].individual_stat_weight_correction,
                     2.0, 0.0);
    AMPS_EXPECT_NEAR(split.particles[index].radius, radius, 0.0);
    AMPS_EXPECT_EQ(Group(split.particles[index], groups), 1);
  }
  const Budget split_budget = EvaluateBudget(split.particles);
  ExpectBudgetsNear(parent_budget, split_budget, 2.0e-15);
  AMPS_EXPECT_NEAR(split_budget.sampled_moments.number,
                   parent_budget.sampled_moments.number, 2.0e-14);
  AMPS_EXPECT_NEAR(split_budget.sampled_moments.speed,
                   parent_budget.sampled_moments.speed, 2.0e-14);
  AMPS_EXPECT_NEAR(split_budget.sampled_moments.absolute_current,
                   parent_budget.sampled_moments.absolute_current, 2.0e-30);

  // Three collinear, positive-x populations are an independent literal
  // fixture.  The two-output merge preserves sum(w), mass, charge, momentum,
  // and kinetic energy.  Because both reconstructed speeds remain positive in
  // this fixture, the sampled mean-speed moment is conserved as well.
  std::vector<Particle> merge_input;
  merge_input.push_back(MakeDustParticle(1.0, radius, mass, charge, 1.0));
  merge_input.push_back(MakeDustParticle(2.0, radius, mass, charge, 2.0));
  merge_input.push_back(MakeDustParticle(3.0, radius, mass, charge, 3.0));
  const Budget merge_input_budget = EvaluateBudget(merge_input);
  AMPS_EXPECT_NEAR(merge_input_budget.represented_number, 6.0, 0.0);
  AMPS_EXPECT_NEAR(merge_input_budget.momentum[0], 14.0 * mass,
                   2.0e-29);
  AMPS_EXPECT_NEAR(merge_input_budget.kinetic_energy, 18.0 * mass,
                   2.0e-29);

  const Result merged = Merge(merge_input, groups, true);
  AMPS_EXPECT_TRUE(merged.ok());
  AMPS_EXPECT_EQ(merged.particles.size(), static_cast<std::size_t>(2));
  AMPS_EXPECT_NEAR(merged.particles[0].individual_stat_weight_correction,
                   3.0, 0.0);
  AMPS_EXPECT_NEAR(merged.particles[1].individual_stat_weight_correction,
                   3.0, 0.0);
  AMPS_EXPECT_EQ(Group(merged.particles[0], groups), 1);
  AMPS_EXPECT_EQ(Group(merged.particles[1], groups), 1);
  const Budget merged_budget = EvaluateBudget(merged.particles);
  ExpectBudgetsNear(merge_input_budget, merged_budget, 4.0e-15);
  AMPS_EXPECT_NEAR(merged_budget.sampled_moments.number,
                   merge_input_budget.sampled_moments.number, 2.0e-14);
  AMPS_EXPECT_NEAR(merged_budget.sampled_moments.velocity[0],
                   merge_input_budget.sampled_moments.velocity[0], 2.0e-14);
  AMPS_EXPECT_NEAR(merged_budget.sampled_moments.speed,
                   merge_input_budget.sampled_moments.speed, 2.0e-14);
  AMPS_EXPECT_NEAR(merged_budget.sampled_moments.charge,
                   merge_input_budget.sampled_moments.charge, 2.0e-30);
  AMPS_EXPECT_NEAR(merged_budget.sampled_moments.current[0],
                   merge_input_budget.sampled_moments.current[0], 2.0e-29);

  // Repeating an operation gives bitwise-identical ordered outputs: no random
  // direction or container iteration order participates in reconstruction.
  const Result repeated = Merge(merge_input, groups, true);
  AMPS_EXPECT_TRUE(repeated.ok());
  for (std::size_t index = 0; index != merged.particles.size(); ++index) {
    AMPS_EXPECT_NEAR(repeated.particles[index].velocity[0],
                     merged.particles[index].velocity[0], 0.0);
    AMPS_EXPECT_NEAR(repeated.particles[index].velocity[1],
                     merged.particles[index].velocity[1], 0.0);
    AMPS_EXPECT_NEAR(repeated.particles[index].velocity[2],
                     merged.particles[index].velocity[2], 0.0);
  }

  // Typed rejections demonstrate that group equality is mandatory and still
  // insufficient when two continuous radii would change grain physics.
  std::vector<Particle> wrong_group = merge_input;
  wrong_group[2].radius = 2.0e-4;
  AMPS_EXPECT_TRUE(Merge(wrong_group, groups).status ==
                   Status::kDifferentSizeGroup);
  std::vector<Particle> changed_radius = merge_input;
  changed_radius[2].radius = 3.0e-6;  // same logarithmic group, different grain
  AMPS_EXPECT_TRUE(Merge(changed_radius, groups).status ==
                   Status::kDifferentIntrinsicGrain);
  std::vector<Particle> wrong_owner = merge_input;
  wrong_owner[1].ownership_cell = 43;
  AMPS_EXPECT_TRUE(Merge(wrong_owner, groups).status ==
                   Status::kDifferentOwnershipCell);
  std::vector<Particle> wrong_species = merge_input;
  wrong_species[1].compiled_species = 8;
  AMPS_EXPECT_TRUE(Merge(wrong_species, groups).status ==
                   Status::kDifferentCompiledSpecies);
  std::vector<Particle> invalid_weight = merge_input;
  invalid_weight[0].individual_stat_weight_correction = 0.0;
  AMPS_EXPECT_TRUE(Merge(invalid_weight, groups).status ==
                   Status::kInvalidParticle);
  AMPS_EXPECT_TRUE(Merge(std::vector<Particle>(2, parent), groups).status ==
                   Status::kInsufficientParticles);
  const SizeGroups unavailable_weights(1.0e-7, 1.0e-3, 4, false);
  AMPS_EXPECT_TRUE(Split(parent, 2, unavailable_weights).status ==
                   Status::kIndividualWeightUnavailable);

  // Scalar energy and the raw sampled speed are independent constraints.  For
  // symmetric velocities {-1,0,+1}, the conservative two-particle energy
  // shell has speeds sqrt(2/3), whose weighted sum differs from the input's
  // literal value 2.  A sampling-sensitive request is rejected transactionally.
  std::vector<Particle> speed_unrepresentable;
  speed_unrepresentable.push_back(
      MakeDustParticle(1.0, radius, mass, charge, -1.0));
  speed_unrepresentable.push_back(
      MakeDustParticle(1.0, radius, mass, charge, 0.0));
  speed_unrepresentable.push_back(
      MakeDustParticle(1.0, radius, mass, charge, 1.0));
  const Result speed_rejected = Merge(speed_unrepresentable, groups, true);
  AMPS_EXPECT_TRUE(speed_rejected.status ==
                   Status::kSampledSpeedUnrepresentable);
  AMPS_EXPECT_TRUE(speed_rejected.particles.empty());
}
