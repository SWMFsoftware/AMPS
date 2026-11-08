#include "pic.h"
#include "src/test/native/test_harness.h"

#ifndef _AMPS_CORE_TEST_T6_REACTIONS_
#define _AMPS_CORE_TEST_T6_REACTIONS_ 0
#endif

#if _AMPS_CORE_TEST_T6_REACTIONS_ == 1
#include "ChargeExchange.h"
#include "ElectronImpact.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>
#endif

/* Stage-T6 charge-exchange and electron-impact validation.
 *
 * This translation unit is profile-gated because the baseline core-test
 * application intentionally does not link either optional reaction archive.
 * The T6 profile must compile H, H2, O, O2, OH, H2O, CO2, electron, H+, H2+,
 * O+, O2+, OH+, and H2O+ and must enable the Burger electron-impact and both
 * charge-exchange source files.  The tests use SI rates [m^3 s^-1], number
 * density [m^-3], speed [m s^-1], momentum [kg m s^-1], and energy [J]; the
 * only electron-temperature argument is the documented eV input.
 *
 * Charge exchange currently exposes an analytic lifetime only.  It has no
 * reaction-product, statistical-weight, lookup-table, interpolation, or
 * Maxwellian-quadrature API, so this suite does not invent those capabilities.
 * Instead it independently reconstructs the two published analytic fits and
 * exercises both branches of the Heerkhuisen effective-relative-speed law.
 * Electron impact does expose tabulated rates and product generation; those
 * contracts receive knot, interpolation, boundary, branching, product,
 * stoichiometric, charge, momentum, and energy checks below.
 */

namespace AMPS {
namespace CoreTests {

// Static-archive anchor.  main.cpp must reference this symbol when the source
// is added to the configured application's archive; registrations exist only
// in the explicitly enabled T6 profile.
void LinkT6ReactionTests() {}

}  // namespace CoreTests
}  // namespace AMPS

#if _AMPS_CORE_TEST_T6_REACTIONS_ == 1
namespace {

double RelativeError(double value, double reference) {
  return std::fabs(value - reference) /
      std::max(std::fabs(reference), 1.0e-300);
}

double IndependentMeanRelativeSpeed(const double* particle_velocity,
                                    const double* plasma_velocity,
                                    double plasma_temperature) {
  // Heerkhuisen et al. (2006), Eq. 8.  For an H Maxwellian with
  // v_th=sqrt(2 k_B T/m_H), omega=|v_n-v_p|/v_th, the speed averaged over the
  // ion distribution is
  //   <g>/v_th = exp(-omega^2)/sqrt(pi)
  //              +(omega+1/(2 omega))*erf(omega).
  // The O(omega^2) series is evaluated below 1e-2, matching the production
  // validity switch while avoiding division by zero at zero drift.
  const double thermal_speed =
      std::sqrt(2.0 * Kbol * plasma_temperature / _H__MASS_);
  double drift2 = 0.0;
  for (int dimension = 0; dimension < 3; ++dimension) {
    const double delta = particle_velocity[dimension] -
        plasma_velocity[dimension];
    drift2 += delta * delta;
  }
  const double omega = std::sqrt(drift2) / thermal_speed;
  if (omega > 1.0e-2)
    return thermal_speed *
        (std::exp(-omega * omega) / std::sqrt(Pi) +
         (omega + 0.5 / omega) * std::erf(omega));
  return thermal_speed * 2.0 * (1.0 + omega * omega / 3.0) /
      std::sqrt(Pi);
}

double IndependentLindsayLifetime(const double* particle_velocity,
                                  const double* plasma_velocity,
                                  double plasma_temperature,
                                  double plasma_number_density) {
  // Lindsay & Stebbings (2005) H-H+ fit. E is the relative collision energy
  // in keV and sigma is returned in m^2:
  // sigma=(4.15-0.531 ln E)^2 C(E) 1e-20,
  // C=(1-exp(-67.3/E))^4.5 for E>4 keV and C=1 below that asymptote.
  // The event lifetime is tau=[n_p <g> sigma]^-1 [s].
  const double relative_speed = IndependentMeanRelativeSpeed(
      particle_velocity, plasma_velocity, plasma_temperature);
  const double energy_kev = 0.5 * _H__MASS_ * relative_speed * relative_speed /
      eV2J / 1000.0;
  const double correction = energy_kev > 4.0
      ? std::pow(1.0 - std::exp(-67.3 / energy_kev), 4.5) : 1.0;
  const double logarithm = energy_kev > 1.0e-10
      ? std::log(energy_kev) : -10.0 * std::log(10.0);
  const double sigma = std::pow(4.15 - 0.531 * logarithm, 2.0) *
      correction * 1.0e-20;
  return 1.0 / (plasma_number_density * relative_speed * sigma);
}

double IndependentMaherLifetime(const double* particle_velocity,
                                const double* plasma_velocity,
                                double plasma_temperature,
                                double plasma_number_density) {
  // Maher & Tinsley (1977), as used by Heerikhuisen et al. (2006):
  // sigma=[1.6-0.0695 ln(g_cm/s)]^2 1e-18 m^2.  The factor 100 converts
  // g from m/s to cm/s, and tau=[n_p <g> sigma]^-1 [s].
  const double relative_speed = IndependentMeanRelativeSpeed(
      particle_velocity, plasma_velocity, plasma_temperature);
  const double sigma =
      std::pow(1.6 - 0.0695 * std::log(relative_speed * 100.0), 2.0) *
      1.0e-18;
  return 1.0 / (plasma_number_density * relative_speed * sigma);
}

typedef double (*ChargeExchangeLifetime)(int, double*, double*, double,
                                         double);

void CheckChargeExchangeModel(ChargeExchangeLifetime model,
                              bool lindsay_stebbings) {
  const double plasma_velocity[3] = {13000.0, -7000.0, 2000.0};
  double zero_drift[3] = {13000.0, -7000.0, 2000.0};
  double high_drift[3] = {163000.0, 18000.0, -9000.0};
  double mutable_plasma_velocity[3] = {
      plasma_velocity[0], plasma_velocity[1], plasma_velocity[2]};
  const double temperature = 8000.0;       // K
  const double number_density = 4.0e6;     // m^-3

  const double expected_zero = lindsay_stebbings
      ? IndependentLindsayLifetime(zero_drift, plasma_velocity, temperature,
                                   number_density)
      : IndependentMaherLifetime(zero_drift, plasma_velocity, temperature,
                                 number_density);
  const double expected_high = lindsay_stebbings
      ? IndependentLindsayLifetime(high_drift, plasma_velocity, temperature,
                                   number_density)
      : IndependentMaherLifetime(high_drift, plasma_velocity, temperature,
                                 number_density);
  const double observed_zero = model(_H_SPEC_, zero_drift,
      mutable_plasma_velocity, temperature, number_density);
  const double observed_high = model(_H_SPEC_, high_drift,
      mutable_plasma_velocity, temperature, number_density);
  AMPS_EXPECT_TRUE(RelativeError(observed_zero, expected_zero) < 3.0e-13);
  AMPS_EXPECT_TRUE(RelativeError(observed_high, expected_high) < 3.0e-13);

  // tau is exactly inverse in the proton number density.  This deterministic
  // ratio rejects a hidden mass-density conversion or omitted density factor.
  const double doubled_density = model(_H_SPEC_, high_drift,
      mutable_plasma_velocity, temperature, 2.0 * number_density);
  AMPS_EXPECT_TRUE(RelativeError(doubled_density, 0.5 * observed_high) <
                   3.0e-15);

  // The current production fit is explicitly H/H-ENA-only.  CO2 is a compiled
  // species in this profile but has no declared cross section; 1e100 s is the
  // legacy typed numeric sentinel for an unsupported/no-reaction species.
  AMPS_EXPECT_EQ(model(_CO2_SPEC_, high_drift, mutable_plasma_velocity,
                       temperature, number_density), 1.0e100);

  // Invalid scalar/vector domains use the same documented no-reaction
  // sentinel.  This boundary contract prevents sqrt/log/division from
  // manufacturing NaN, infinity, or a negative event lifetime.
  double invalid_velocity[3] = {
      std::numeric_limits<double>::quiet_NaN(), 0.0, 0.0};
  AMPS_EXPECT_EQ(model(_H_SPEC_, high_drift, mutable_plasma_velocity,
                       0.0, number_density), 1.0e100);
  AMPS_EXPECT_EQ(model(_H_SPEC_, high_drift, mutable_plasma_velocity,
                       temperature, 0.0), 1.0e100);
  AMPS_EXPECT_EQ(model(_H_SPEC_, invalid_velocity, mutable_plasma_velocity,
                       temperature, number_density), 1.0e100);
  AMPS_EXPECT_EQ(model(_H_SPEC_, NULL, mutable_plasma_velocity,
                       temperature, number_density), 1.0e100);
}

double LiteralH2ORate(int channel, double log_rate) {
  // Burger's digitized tables store log10(k/[cm^3 s^-1]); 1 cm^3=1e-6 m^3.
  (void)channel;
  return 1.0e-6 * std::pow(10.0, log_rate);
}

std::vector<double> H2ORates(double temperature_ev) {
  std::vector<double> rates(
      ElectronImpact::H2O::Burger2010SSR::nReactionChannels, 0.0);
  ElectronImpact::Burger2010SSR::GetTotalRateCoefficient(
      &rates[0], temperature_ev,
      ElectronImpact::H2O::Burger2010SSR::nReactionChannels,
      &ElectronImpact::H2O::Burger2010SSR::log10RateCoefficientTable[0][0],
      ElectronImpact::H2O::Burger2010SSR::minTemp);
  return rates;
}

double Sum(const std::vector<double>& values) {
  double result = 0.0;
  for (std::size_t index = 0; index < values.size(); ++index)
    result += values[index];
  return result;
}

struct ElementCount {
  int hydrogen;
  int oxygen;
};

ElementCount SpeciesAtoms(int species) {
  // Product identity is checked by an independent elemental map rather than
  // by summing floating molecular masses.  AMPS' historical _H2O__MASS_
  // literal and integer-amu fragment constants differ by about 8.5e-4, so a
  // rest-mass equality would test unrelated constant provenance, not reaction
  // stoichiometry. Electrons carry no nuclei.
  if (species == _H_SPEC_ || species == _H_PLUS_SPEC_)
    return ElementCount{1, 0};
  if (species == _H2_SPEC_ || species == _H2_PLUS_SPEC_)
    return ElementCount{2, 0};
  if (species == _O_SPEC_ || species == _O_PLUS_SPEC_)
    return ElementCount{0, 1};
  if (species == _OH_SPEC_ || species == _OH_PLUS_SPEC_)
    return ElementCount{1, 1};
  if (species == _H2O_SPEC_ || species == _H2O_PLUS_SPEC_)
    return ElementCount{2, 1};
  if (species == _ELECTRON_SPEC_) return ElementCount{0, 0};
  throw AMPS::Test::Failure("unknown H2O reaction product in atomic oracle");
}

int SpeciesChargeNumber(int species) {
  if (species == _ELECTRON_SPEC_) return -1;
  if (species == _H_PLUS_SPEC_ || species == _H2_PLUS_SPEC_ ||
      species == _O_PLUS_SPEC_ || species == _O2_PLUS_SPEC_ ||
      species == _OH_PLUS_SPEC_ || species == _H2O_PLUS_SPEC_)
    return 1;
  return 0;
}

double KinematicMass(int species) {
  // These are the masses explicitly used by the Burger kinematic closure.
  // They are stated here by product identity so the conservation sum is not
  // read back from the mutable production mass table under test.  Ion masses
  // use the corresponding neutral/nuclear mass; electron binding-energy mass
  // defects are outside this nonrelativistic empirical model.
  if (species == _ELECTRON_SPEC_) return _ELECTRON__MASS_;
  if (species == _H_SPEC_ || species == _H_PLUS_SPEC_) return _H__MASS_;
  if (species == _H2_SPEC_ || species == _H2_PLUS_SPEC_) return _H2__MASS_;
  if (species == _O_SPEC_ || species == _O_PLUS_SPEC_) return _O__MASS_;
  if (species == _O2_SPEC_ || species == _O2_PLUS_SPEC_) return _O2__MASS_;
  if (species == _OH_SPEC_ || species == _OH_PLUS_SPEC_) return _OH__MASS_;
  if (species == _H2O_SPEC_ || species == _H2O_PLUS_SPEC_)
    return _H2O__MASS_;
  throw AMPS::Test::Failure("unknown H2O reaction product in mass oracle");
}

void CheckH2OProductEvent(double temperature_ev, int reaction_channel,
                          int number_products, const int* products,
                          const double* velocities) {
  static const int expected_count[5] = {4, 4, 4, 3, 3};
  static const int expected_products[5][4] = {
      {_O_PLUS_SPEC_, _H2_SPEC_, _ELECTRON_SPEC_, _ELECTRON_SPEC_},
      {_OH_SPEC_, _H_PLUS_SPEC_, _ELECTRON_SPEC_, _ELECTRON_SPEC_},
      {_OH_PLUS_SPEC_, _H_SPEC_, _ELECTRON_SPEC_, _ELECTRON_SPEC_},
      {_H2O_PLUS_SPEC_, _ELECTRON_SPEC_, _ELECTRON_SPEC_, -1},
      {_OH_SPEC_, _H_SPEC_, _ELECTRON_SPEC_, -1}};
  static const double threshold_ev[5] = {4.0, 2.8, 2.3, 1.5, 1.0};

  AMPS_EXPECT_TRUE(reaction_channel >= 0 && reaction_channel < 5);
  AMPS_EXPECT_EQ(number_products, expected_count[reaction_channel]);
  int hydrogen = 0, oxygen = 0, charge = 0;
  double momentum[3] = {0.0, 0.0, 0.0};
  double kinetic_energy = 0.0;
  for (int product = 0; product < number_products; ++product) {
    AMPS_EXPECT_EQ(products[product],
                   expected_products[reaction_channel][product]);
    const ElementCount atoms = SpeciesAtoms(products[product]);
    hydrogen += atoms.hydrogen;
    oxygen += atoms.oxygen;
    charge += SpeciesChargeNumber(products[product]);
    const double mass = KinematicMass(products[product]);
    AMPS_EXPECT_TRUE(mass > 0.0 && std::isfinite(mass));
    double speed2 = 0.0;
    for (int dimension = 0; dimension < 3; ++dimension) {
      const double velocity = velocities[3 * product + dimension];
      AMPS_EXPECT_TRUE(std::isfinite(velocity));
      momentum[dimension] += mass * velocity;
      speed2 += velocity * velocity;
    }
    kinetic_energy += 0.5 * mass * speed2;
  }
  AMPS_EXPECT_EQ(hydrogen, 2);
  AMPS_EXPECT_EQ(oxygen, 1);
  // A neutral parent plus one incident electron has net charge -e.  Ionizing
  // branches return one +e ion and two electrons; dissociation returns one
  // electron.  Charge numbers avoid trusting deck initialization as an oracle.
  AMPS_EXPECT_EQ(charge, -1);

  // Product kinematics are generated in the stationary-parent frame.  The
  // incident electron direction is isotropic and not returned, so only the
  // exactly invariant momentum magnitude can be checked:
  // |sum m_i v_i|=sqrt(2 m_e T_e e), while sum(1/2 m_i v_i^2)
  // =(T_e-T_threshold)e.  These equations are independent of the production
  // random decomposition and rotation.
  const double momentum_magnitude = std::sqrt(
      momentum[0] * momentum[0] + momentum[1] * momentum[1] +
      momentum[2] * momentum[2]);
  const double expected_momentum =
      std::sqrt(2.0 * _ELECTRON__MASS_ * temperature_ev * eV2J);
  const double expected_energy =
      (temperature_ev - threshold_ev[reaction_channel]) * eV2J;
  AMPS_EXPECT_TRUE(RelativeError(momentum_magnitude, expected_momentum) <
                   2.0e-10);
  AMPS_EXPECT_TRUE(RelativeError(kinetic_energy, expected_energy) < 2.0e-10);
}

}  // namespace

AMPS_TEST_CASE(
    "T6-CHARGE-EXCHANGE-LINDSAY-001",
    "Lindsay-Stebbings analytic H-H+ lifetime, branches and density scaling") {
  CheckChargeExchangeModel(&ChargeExchange::LindsayStebbings::LifeTime, true);
}

AMPS_TEST_CASE(
    "T6-CHARGE-EXCHANGE-MAHER-001",
    "Maher-Tinsley analytic H-H+ lifetime, branches and density scaling") {
  CheckChargeExchangeModel(&ChargeExchange::MaherTinsley::LifeTime, false);
}

AMPS_TEST_CASE(
    "T6-CHARGE-EXCHANGE-SELECTOR-001",
    "configured charge-exchange wrapper dispatches to the selected model") {
  double particle_velocity[3] = {125000.0, -21000.0, 9000.0};
  double plasma_velocity[3] = {5000.0, -1000.0, 2000.0};
  const double selected = ChargeExchange::LifeTime(
      _H_SPEC_, particle_velocity, plasma_velocity, 12000.0, 7.0e6);
#if _CHARGE_EXCHANGE__MODEL_ == _CHARGE_EXCHANGE__LINDSAY_STEBBINGS_
  const double direct = ChargeExchange::LindsayStebbings::LifeTime(
      _H_SPEC_, particle_velocity, plasma_velocity, 12000.0, 7.0e6);
#elif _CHARGE_EXCHANGE__MODEL_ == _CHARGE_EXCHANGE__MAHER_TINSLEY_
  const double direct = ChargeExchange::MaherTinsley::LifeTime(
      _H_SPEC_, particle_velocity, plasma_velocity, 12000.0, 7.0e6);
#else
#error T6 reaction profile must select a supported charge-exchange model
#endif
  AMPS_EXPECT_NEAR(selected, direct, 0.0);
}

AMPS_TEST_CASE(
    "T6-ELECTRON-IMPACT-RATES-001",
    "Burger H2O thresholds, knots, log interpolation, clamp and invalid input") {
  static const double knot_log_rate[5] = {
      -10.10200, -9.22257, -8.83602, -8.10094, -7.85134};
  static const double next_log_rate[5] = {
      -9.96663, -9.10097, -8.74298, -8.02529, -7.78419};
  static const double final_log_rate[5] = {
      -8.44040, -7.72228, -7.65143, -7.15907, -7.01367};

  std::vector<double> rates = H2ORates(0.999);
  AMPS_EXPECT_NEAR(Sum(rates), 0.0, 0.0);
  rates = H2ORates(10.0);  // log10(T_e)=1.00, table index 25
  for (int channel = 0; channel < 5; ++channel)
    AMPS_EXPECT_TRUE(RelativeError(
        rates[channel], LiteralH2ORate(channel, knot_log_rate[channel])) <
        3.0e-14);

  // The table is linear in log10(k) on a uniform log10(T_e) grid.  At
  // log10(T_e)=1.02 the expected rate is the geometric mean of indices 25/26.
  rates = H2ORates(std::pow(10.0, 1.02));
  for (int channel = 0; channel < 5; ++channel) {
    const double expected = LiteralH2ORate(
        channel, 0.5 * (knot_log_rate[channel] + next_log_rate[channel]));
    AMPS_EXPECT_TRUE(RelativeError(rates[channel], expected) < 3.0e-14);
  }

  rates = H2ORates(1000.0);  // documented high-side constant clamp
  for (int channel = 0; channel < 5; ++channel)
    AMPS_EXPECT_TRUE(RelativeError(
        rates[channel], LiteralH2ORate(channel, final_log_rate[channel])) <
        3.0e-14);

  // One independently transcribed 10 eV knot from every other Burger parent
  // table verifies that namespace wrappers own the correct channel count,
  // flattened table base, and threshold vector.  This guards against testing
  // H2O thoroughly while silently dispatching O2/H2/H/O through the wrong
  // storage.
  const double expected_o2 = LiteralH2ORate(0, -8.79029) +
      LiteralH2ORate(1, -8.15209) + LiteralH2ORate(2, -8.15085);
  const double expected_h2 = LiteralH2ORate(0, -9.82660) +
      LiteralH2ORate(1, -8.25713) + LiteralH2ORate(2, -9.48428);
  AMPS_EXPECT_TRUE(RelativeError(
      ElectronImpact::O2::Burger2010SSR::GetTotalRateCoefficient(10.0),
      expected_o2) < 3.0e-14);
  AMPS_EXPECT_TRUE(RelativeError(
      ElectronImpact::H2::Burger2010SSR::GetTotalRateCoefficient(10.0),
      expected_h2) < 3.0e-14);
  AMPS_EXPECT_TRUE(RelativeError(
      ElectronImpact::H::Burger2010SSR::GetTotalRateCoefficient(10.0),
      LiteralH2ORate(0, -8.30222)) < 3.0e-14);
  AMPS_EXPECT_TRUE(RelativeError(
      ElectronImpact::O::Burger2010SSR::GetTotalRateCoefficient(10.0),
      LiteralH2ORate(0, -8.19213)) < 3.0e-14);

  const double invalid_temperatures[3] = {
      0.0, -1.0, std::numeric_limits<double>::quiet_NaN()};
  for (int invalid = 0; invalid < 3; ++invalid) {
    rates = H2ORates(invalid_temperatures[invalid]);
    for (int channel = 0; channel < 5; ++channel)
      AMPS_EXPECT_NEAR(rates[channel], 0.0, 0.0);
  }
}

AMPS_TEST_CASE(
    "T6-ELECTRON-IMPACT-YIELDS-001",
    "Burger H2O multiplicity-weighted yields and exact threshold capability") {
  const std::vector<double> rates = H2ORates(10.0);
  const double total = Sum(rates);
  const double expected_oh = (rates[1] + rates[4]) / total;
  const double expected_electron =
      (2.0 * (rates[0] + rates[1] + rates[2] + rates[3]) + rates[4]) /
      total;
  AMPS_EXPECT_TRUE(RelativeError(
      ElectronImpact::H2O::Burger2010SSR::GetSpeciesReactionYield(
          _OH_SPEC_, 10.0), expected_oh) < 2.0e-14);
  AMPS_EXPECT_TRUE(RelativeError(
      ElectronImpact::H2O::Burger2010SSR::GetSpeciesReactionYield(
          _ELECTRON_SPEC_, 10.0), expected_electron) < 2.0e-14);
  AMPS_EXPECT_NEAR(
      ElectronImpact::H2O::Burger2010SSR::GetSpeciesReactionYield(
          _CO2_SPEC_, 10.0), 0.0, 0.0);

  AMPS_EXPECT_TRUE(ElectronImpact::ModelAvailable(_H2O_SPEC_));
  AMPS_EXPECT_TRUE(ElectronImpact::ModelAvailable(_O2_SPEC_));
  AMPS_EXPECT_TRUE(ElectronImpact::ModelAvailable(_H2_SPEC_));
  AMPS_EXPECT_TRUE(ElectronImpact::ModelAvailable(_H_SPEC_));
  AMPS_EXPECT_TRUE(ElectronImpact::ModelAvailable(_O_SPEC_));
  AMPS_EXPECT_FALSE(ElectronImpact::ModelAvailable(_CO2_SPEC_));
  // Channel 4 opens at exactly 1 eV; capability and rate APIs must use the
  // same inclusive boundary and must both reject a point below it.
  AMPS_EXPECT_TRUE(ElectronImpact::hasDaughterSpec(
      _H2O_SPEC_, _OH_SPEC_, 1.0));
  AMPS_EXPECT_FALSE(ElectronImpact::hasDaughterSpec(
      _H2O_SPEC_, _OH_SPEC_, std::nextafter(1.0, 0.0)));
}

AMPS_TEST_CASE(
    "T6-ELECTRON-IMPACT-PRODUCTS-001",
    "Burger H2O product branching, identity and physical conservation") {
  const double temperature_ev = 10.0;
  const std::vector<double> rates = H2ORates(temperature_ev);
  const double total = Sum(rates);
  const int samples = 6000;
  int counts[5] = {0, 0, 0, 0, 0};
  rnd_seed(2606101);
  for (int sample = 0; sample < samples; ++sample) {
    int reaction_channel = -2;
    int number_products = -1;
    int* products = NULL;
    double* velocities = NULL;
    ElectronImpact::H2O::Burger2010SSR::GenerateReactionProducts(
        temperature_ev, reaction_channel, number_products, products,
        velocities);
    CheckH2OProductEvent(temperature_ev, reaction_channel, number_products,
                         products, velocities);
    ++counts[reaction_channel];
  }
  for (int channel = 0; channel < 5; ++channel) {
    const double expected = rates[channel] / total;
    const double observed = static_cast<double>(counts[channel]) / samples;
    const double six_sigma =
        6.0 * std::sqrt(expected * (1.0 - expected) / samples);
    AMPS_EXPECT_TRUE(std::fabs(observed - expected) < six_sigma + 0.003);
  }

  int reaction_channel = 44;
  int number_products = 44;
  int* products = NULL;
  double* velocities = NULL;
  ElectronImpact::H2O::Burger2010SSR::GenerateReactionProducts(
      0.5, reaction_channel, number_products, products, velocities);
  AMPS_EXPECT_EQ(reaction_channel, -1);
  AMPS_EXPECT_EQ(number_products, 0);
}

AMPS_TEST_CASE(
    "T6-ELECTRON-IMPACT-CONDITIONED-001",
    "conditioned Burger products sample eligible channels and reject gaps") {
  rnd_seed(2606102);
  for (int sample = 0; sample < 300; ++sample) {
    int reaction_channel = -2;
    int number_products = -1;
    int* products = NULL;
    double* velocities = NULL;
    ElectronImpact::H2O::Burger2010SSR::GenerateGivenProducts(
        _OH_SPEC_, 10.0, reaction_channel, number_products, products,
        velocities);
    AMPS_EXPECT_TRUE(reaction_channel == 1 || reaction_channel == 4);
    CheckH2OProductEvent(10.0, reaction_channel, number_products, products,
                         velocities);
    AMPS_EXPECT_TRUE(std::find(products, products + number_products,
                               _OH_SPEC_) != products + number_products);
  }

  int reaction_channel = 17;
  int number_products = 17;
  int* products = NULL;
  double* velocities = NULL;
  ElectronImpact::H2O::Burger2010SSR::GenerateGivenProducts(
      _CO2_SPEC_, 10.0, reaction_channel, number_products, products,
      velocities);
  AMPS_EXPECT_EQ(reaction_channel, -1);
  AMPS_EXPECT_EQ(number_products, 0);
}

#endif  // _AMPS_CORE_TEST_T6_REACTIONS_ == 1
