#include "test_harness.h"

#include "Sputtering.h"

#include <cmath>
#include <cstddef>
#include <functional>
#include <limits>
#include <stdexcept>

double amps_sputtering_test_random = 0.5;

namespace {
// Independent SI literals used to map a tabulated projectile energy to the
// speed accepted by the production API, v=sqrt(2 E e/m).  Keeping these out of
// the production mass helper means an erroneous production mass selection is
// visible to the table tests.
const double kElectronVoltJoule = 1.602176634e-19;
const double kAtomicMassKilogramOracle = 1.66053906660e-27;
const double kMassH = 1.0 * kAtomicMassKilogramOracle;
const double kMassO = 16.0 * kAtomicMassKilogramOracle;
const double kMassO2 = 32.0 * kAtomicMassKilogramOracle;
const double kMassH2O = 18.0 * kAtomicMassKilogramOracle;

double SpeedFromEnergyEv(double energy_ev, double mass_kg) {
  return std::sqrt(2.0 * energy_ev * kElectronVoltJoule / mass_kg);
}

void ExpectRelativeNear(double value, double reference, double relative) {
  AMPS_EXPECT_NEAR(value, reference,
                   relative * std::max(1.0, std::fabs(reference)));
}

void ExpectInvalid(const std::function<void()>& operation) {
  bool rejected = false;
  try {
    operation();
  } catch (const std::invalid_argument&) {
    rejected = true;
  }
  AMPS_EXPECT_TRUE(rejected);
}

void ExpectOverflow(const std::function<void()>& operation) {
  bool rejected = false;
  try {
    operation();
  } catch (const std::overflow_error&) {
    rejected = true;
  }
  AMPS_EXPECT_TRUE(rejected);
}

double Dot(const std::array<double, 3>& left,
           const std::array<double, 3>& right) {
  return left[0] * right[0] + left[1] * right[1] + left[2] * right[2];
}

double Norm(const std::array<double, 3>& value) {
  return std::sqrt(Dot(value, value));
}

double H2OCumulativeEnergyProbability(double energy_joule) {
  // The retained H2O law is p(E)=2 U E/(E+U)^3 [J^-1], U=0.055 eV.
  // Its independently integrated CDF is [E/(E+U)]^2.  This is used only as
  // an oracle for the sampled relative statistical weights.
  const double binding_joule = 0.055 * kElectronVoltJoule;
  const double ratio = energy_joule / (energy_joule + binding_joule);
  return ratio * ratio;
}

double O2CumulativeEnergyProbability(double energy_joule) {
  // For p(E)=U/(E+U)^2 [J^-1], U=0.015 eV, direct integration from zero gives
  // F(E)=E/(E+U).  Both analytic PDFs integrate to unity on [0,infinity).
  const double binding_joule = 0.015 * kElectronVoltJoule;
  return energy_joule / (energy_joule + binding_joule);
}

double IntegrateReturnedLogWeight(int species, double v_min, double v_max,
                                  double log_pdf_max, int intervals) {
  // The sampler is uniform in ln(v).  If production returns
  // w(v)=p_ln(v)/max(p_ln), midpoint quadrature reconstructs the physical
  // probability over its finite speed window as
  //   integral p_ln(v)dln(v) = Delta_ln(v) max(p_ln) mean(w).
  // This checks normalization without using a production-generated expected
  // value.  Speeds are m/s; p_ln and w are dimensionless.
  double sum = 0.0;
  for (int i = 0; i < intervals; ++i) {
    amps_sputtering_test_random = (i + 0.5) / intervals;
    double speed = 0.0;
    double weight = 0.0;
    Sputtering::GetSputteringSpeed(species, speed, weight);
    AMPS_EXPECT_TRUE(speed >= v_min && speed <= v_max);
    AMPS_EXPECT_TRUE(weight >= 0.0 && weight <= 1.0 + 1.0e-13);
    sum += weight;
  }
  return std::log(v_max / v_min) * log_pdf_max * sum / intervals;
}
}  // namespace

AMPS_TEST_CASE("SPUTTERING-TEOLIS-TABLE-001",
               "Teolis table endpoints, thresholds, and log-log interpolation") {
  // These decimal yields are independently pinned from rows 0, 1, and 99 of
  // the retained Benjamin Teolis table.  This fixture validates table access;
  // it does not assert external laboratory provenance beyond that retained
  // model attribution.
  const double energy0 = 10.0;
  const double energy_max = 6.0e6;
  const double ratio = std::pow(energy_max / energy0, 1.0 / 99.0);
  const double energy1 = energy0 * ratio;
  const double geometric_midpoint = std::sqrt(energy0 * energy1);

  ExpectRelativeNear(Sputtering::Ice::Teolis::SputteringYield(
      SpeedFromEnergyEv(energy0, kMassH), _H_PLUS_SPEC_, _O2_SPEC_),
      0.00903, 2.0e-12);
  ExpectRelativeNear(Sputtering::Ice::Teolis::SputteringYield(
      SpeedFromEnergyEv(energy1, kMassH), _H_PLUS_SPEC_, _O2_SPEC_),
      0.01125, 2.0e-12);
  ExpectRelativeNear(Sputtering::Ice::Teolis::SputteringYield(
      SpeedFromEnergyEv(geometric_midpoint, kMassH), _H_PLUS_SPEC_,
      _O2_SPEC_), std::sqrt(0.00903 * 0.01125), 2.0e-12);

  // The upper endpoint is part of the 100-row grid and must return row 99,
  // not read a nonexistent row 100.  Check a second target/projectile pair so
  // target and projectile indexing cannot be accidentally interchanged.
  ExpectRelativeNear(Sputtering::Ice::Teolis::SputteringYield(
      SpeedFromEnergyEv(energy_max, kMassO2), _O2_PLUS_SPEC_, _H2O_SPEC_),
      9811.88600, 2.0e-12);
  ExpectRelativeNear(Sputtering::Ice::Teolis::SputteringYield(
      SpeedFromEnergyEv(energy_max, kMassO), _O_PLUS_SPEC_, _O2_SPEC_),
      35.72758, 2.0e-12);

  AMPS_EXPECT_EQ(Sputtering::Ice::Teolis::SputteringYield(
      SpeedFromEnergyEv(9.0, kMassH), _H_PLUS_SPEC_, _O2_SPEC_), 0.0);
  AMPS_EXPECT_EQ(Sputtering::Ice::Teolis::SputteringYield(
      SpeedFromEnergyEv(6.1e6, kMassH), _H_PLUS_SPEC_, _O2_SPEC_), 0.0);
  AMPS_EXPECT_EQ(Sputtering::Ice::Teolis::SputteringYield(
      std::numeric_limits<double>::quiet_NaN(), _H_PLUS_SPEC_, _O2_SPEC_),
      0.0);
}

AMPS_TEST_CASE("SPUTTERING-RUBIN-TABLE-001",
               "Rubin table thresholds, boundaries, and log-log interpolation") {
  // Rubin's retained grid has log10(E/eV)=1.47+0.1 i.  Row 0 is explicitly
  // below the model's one-molecule/ion threshold; rows 1 and 2 for H+ are
  // 1.2189 and 1.6346.  Half a log-energy cell must therefore return their
  // geometric, not arithmetic, mean.
  const double below_threshold = std::pow(10.0, 1.52);
  const double row1_energy = std::pow(10.0, 1.57);
  const double midpoint_energy = std::pow(10.0, 1.62);
  AMPS_EXPECT_EQ(Sputtering::Ice::Rubin::SputteringYield(
      SpeedFromEnergyEv(below_threshold, kMassH), _H_PLUS_SPEC_,
      _H2O_SPEC_), 0.0);
  ExpectRelativeNear(Sputtering::Ice::Rubin::SputteringYield(
      SpeedFromEnergyEv(row1_energy, kMassH), _H_PLUS_SPEC_, _H2O_SPEC_),
      1.2189, 3.0e-12);
  ExpectRelativeNear(Sputtering::Ice::Rubin::SputteringYield(
      SpeedFromEnergyEv(midpoint_energy, kMassH), _H_PLUS_SPEC_,
      _H2O_SPEC_), std::sqrt(1.2189 * 1.6346), 3.0e-12);

  // Zero speed formerly reached log10(0) and converted -infinity to an array
  // index.  The documented 1 GeV ceiling and nonfinite inputs must also fail
  // closed with zero yield.
  AMPS_EXPECT_EQ(Sputtering::Ice::Rubin::SputteringYield(
      0.0, _H_PLUS_SPEC_, _H2O_SPEC_), 0.0);
  AMPS_EXPECT_EQ(Sputtering::Ice::Rubin::SputteringYield(
      SpeedFromEnergyEv(1.01e9, kMassO), _O_PLUS_SPEC_, _H2O_SPEC_), 0.0);
  AMPS_EXPECT_EQ(Sputtering::Ice::Rubin::SputteringYield(
      std::numeric_limits<double>::infinity(), _O2_PLUS_SPEC_, _H2O_SPEC_),
      0.0);
}

AMPS_TEST_CASE("SPUTTERING-PRODUCT-CATALOGUE-001",
               "master dispatch preserves product identity and yield ordering") {
  // Verify each supported ice model's compile-time product catalogue before
  // checking the configured master dispatch.  Rubin produces H2O only;
  // Teolis produces O2 and H2O in the same order as the yield table's target
  // axis.  Species IDs are identities, not multiplicative source weights.
  AMPS_EXPECT_EQ(Sputtering::Ice::Rubin::nSpecTarget, 1);
  AMPS_EXPECT_EQ(Sputtering::Ice::Rubin::SpecTarget[0], _H2O_SPEC_);
  AMPS_EXPECT_EQ(Sputtering::Ice::Teolis::nSpecTarget, 2);
  AMPS_EXPECT_EQ(Sputtering::Ice::Teolis::SpecTarget[0], _O2_SPEC_);
  AMPS_EXPECT_EQ(Sputtering::Ice::Teolis::SpecTarget[1], _H2O_SPEC_);

  // The configured default is Teolis ice: O2 then H2O.  Calling Init twice
  // exercises catalogue ownership/replacement as well as the public linkage;
  // neither operation may create a duplicate product or stale pointer state.
  Sputtering::Init();
  Sputtering::Init();
  AMPS_EXPECT_EQ(Sputtering::nSpecTarget, 2);
  AMPS_EXPECT_TRUE(Sputtering::SpecTarget != NULL);
  AMPS_EXPECT_TRUE(Sputtering::Yield != NULL);
  AMPS_EXPECT_EQ(Sputtering::SpecTarget[0], _O2_SPEC_);
  AMPS_EXPECT_EQ(Sputtering::SpecTarget[1], _H2O_SPEC_);

  const double speed = SpeedFromEnergyEv(100.0, kMassH);
  const int* products = NULL;
  double* yields = NULL;
  const int count = Sputtering::GetSputteringYield(
      speed, _H_PLUS_SPEC_, products, yields);
  AMPS_EXPECT_EQ(count, 2);
  AMPS_EXPECT_TRUE(products == Sputtering::SpecTarget);
  AMPS_EXPECT_TRUE(yields == Sputtering::Yield);
  for (int i = 0; i < count; ++i)
    ExpectRelativeNear(yields[i],
        Sputtering::Ice::Teolis::SputteringYield(
            speed, _H_PLUS_SPEC_, products[i]), 2.0e-13);

  // Unsupported projectile identities have no tabulated yield.  The legacy
  // public contract is zero rather than an invented mass or charge mapping.
  Sputtering::GetSputteringYield(speed, 987654, products, yields);
  AMPS_EXPECT_EQ(yields[0], 0.0);
  AMPS_EXPECT_EQ(yields[1], 0.0);
}

AMPS_TEST_CASE("SPUTTERING-SPEED-WEIGHT-001",
               "H2O and O2 log-speed samples reproduce analytic relative weights") {
  // Quantile u=1/2 of a uniform-ln(v) proposal is sqrt(v_min v_max).
  amps_sputtering_test_random = 0.5;
  double speed = 0.0;
  double weight = 0.0;
  Sputtering::GetSputteringSpeed(_H2O_SPEC_, speed, weight);
  ExpectRelativeNear(speed, std::sqrt(10.0 * 6000.0), 2.0e-13);
  const double h2o_u = 0.055 * kElectronVoltJoule;
  const double h2o_energy = 0.5 * kMassH2O * speed * speed;
  const double expected_h2o =
      (4.0 * h2o_u * h2o_energy * h2o_energy /
       std::pow(h2o_energy + h2o_u, 3.0)) / (16.0 / 27.0);
  ExpectRelativeNear(weight, expected_h2o, 2.0e-13);

  amps_sputtering_test_random = 0.5;
  Sputtering::GetSputteringSpeed(_O2_SPEC_, speed, weight);
  ExpectRelativeNear(speed, std::sqrt(10.0 * 2000.0), 2.0e-13);
  const double o2_u = 0.015 * kElectronVoltJoule;
  const double o2_energy = 0.5 * kMassO2 * speed * speed;
  const double expected_o2 =
      (2.0 * o2_u * o2_energy / std::pow(o2_energy + o2_u, 2.0)) / 0.5;
  ExpectRelativeNear(weight, expected_o2, 2.0e-13);

  // The H2O log-PDF reaches 16/27 at E=2U; the O2 log-PDF reaches 1/2 at
  // E=U.  Request those speeds through the sampler and verify that the
  // dimensionless correction is normalized to one at each maximum.
  const double h2o_peak_speed = std::sqrt(4.0 * h2o_u / kMassH2O);
  amps_sputtering_test_random =
      std::log(h2o_peak_speed / 10.0) / std::log(6000.0 / 10.0);
  Sputtering::GetSputteringSpeed(_H2O_SPEC_, speed, weight);
  ExpectRelativeNear(speed, h2o_peak_speed, 2.0e-13);
  AMPS_EXPECT_NEAR(weight, 1.0, 2.0e-13);

  const double o2_peak_speed = std::sqrt(2.0 * o2_u / kMassO2);
  amps_sputtering_test_random =
      std::log(o2_peak_speed / 10.0) / std::log(2000.0 / 10.0);
  Sputtering::GetSputteringSpeed(_O2_SPEC_, speed, weight);
  ExpectRelativeNear(speed, o2_peak_speed, 2.0e-13);
  AMPS_EXPECT_NEAR(weight, 1.0, 2.0e-13);
}

AMPS_TEST_CASE("SPUTTERING-WEIGHT-NORMALIZATION-001",
               "relative weights integrate to each truncated analytic energy law") {
  const int intervals = 8192;
  const double h2o_integral = IntegrateReturnedLogWeight(
      _H2O_SPEC_, 10.0, 6000.0, 16.0 / 27.0, intervals);
  const double h2o_e_min = 0.5 * kMassH2O * 10.0 * 10.0;
  const double h2o_e_max = 0.5 * kMassH2O * 6000.0 * 6000.0;
  const double h2o_expected =
      H2OCumulativeEnergyProbability(h2o_e_max) -
      H2OCumulativeEnergyProbability(h2o_e_min);
  AMPS_EXPECT_NEAR(h2o_integral, h2o_expected, 2.0e-8);

  const double o2_integral = IntegrateReturnedLogWeight(
      _O2_SPEC_, 10.0, 2000.0, 0.5, intervals);
  const double o2_e_min = 0.5 * kMassO2 * 10.0 * 10.0;
  const double o2_e_max = 0.5 * kMassO2 * 2000.0 * 2000.0;
  const double o2_expected = O2CumulativeEnergyProbability(o2_e_max) -
                             O2CumulativeEnergyProbability(o2_e_min);
  AMPS_EXPECT_NEAR(o2_integral, o2_expected, 2.0e-8);
}

AMPS_TEST_CASE("SPUTTERING-ANGLE-SELECTOR-001",
               "default and named emission-angle selectors are explicit and typed") {
  using Sputtering::Source::EmissionAngleLaw;
  AMPS_EXPECT_TRUE(Sputtering::Source::kDefaultEmissionAngleLaw ==
                   EmissionAngleLaw::KnudsenCosine);
  AMPS_EXPECT_TRUE(Sputtering::Source::ParseEmissionAngleLaw(
      "cosine-knudsen") == EmissionAngleLaw::KnudsenCosine);
  AMPS_EXPECT_TRUE(Sputtering::Source::ParseEmissionAngleLaw(
      "uniform-hemisphere") == EmissionAngleLaw::UniformHemisphere);
  AMPS_EXPECT_EQ(std::string(Sputtering::Source::EmissionAngleLawName(
      EmissionAngleLaw::KnudsenCosine)), std::string("cosine-knudsen"));
  AMPS_EXPECT_EQ(std::string(Sputtering::Source::EmissionAngleLawName(
      EmissionAngleLaw::UniformHemisphere)),
      std::string("uniform-hemisphere"));

  ExpectInvalid([]() {
    (void)Sputtering::Source::ParseEmissionAngleLaw("unknown");
  });
  ExpectInvalid([]() {
    (void)Sputtering::Source::DirectionFromUnitSamples(
        std::array<double, 3>{{0.0, 0.0, 1.0}}, 0.5, 0.5,
        static_cast<EmissionAngleLaw>(777));
  });
  ExpectInvalid([]() {
    (void)Sputtering::Source::DirectionFromUnitSamples(
        std::array<double, 3>{{0.0, 0.0, 0.0}}, 0.5, 0.5);
  });
}

AMPS_TEST_CASE("SPUTTERING-ANGLE-MOMENTS-001",
               "cosine and uniform hemisphere samples close analytic angular moments") {
  using Sputtering::Source::EmissionAngleLaw;
  const std::array<double, 3> normal = {{2.0, -1.0, 3.0}};
  const double inverse_normal = 1.0 / Norm(normal);
  const std::array<double, 3> unit_normal = {{
      normal[0] * inverse_normal, normal[1] * inverse_normal,
      normal[2] * inverse_normal}};
  const EmissionAngleLaw laws[2] = {
      EmissionAngleLaw::KnudsenCosine,
      EmissionAngleLaw::UniformHemisphere};
  const int n_mu = 400;
  const int n_phi = 64;

  for (int law_index = 0; law_index < 2; ++law_index) {
    const EmissionAngleLaw law = laws[law_index];
    const Sputtering::Source::AngularMoments analytic =
        Sputtering::Source::GetAngularMoments(law);
    double mean_mu = 0.0;
    double mean_mu_squared = 0.0;
    std::array<double, 3> mean_direction = {{0.0, 0.0, 0.0}};
    for (int i = 0; i < n_mu; ++i) {
      for (int j = 0; j < n_phi; ++j) {
        const std::array<double, 3> direction =
            Sputtering::Source::DirectionFromUnitSamples(
                normal, (i + 0.5) / n_mu, (j + 0.5) / n_phi, law);
        const double mu = Dot(direction, unit_normal);
        AMPS_EXPECT_NEAR(Norm(direction), 1.0, 8.0e-15);
        AMPS_EXPECT_TRUE(mu >= -2.0e-15);
        mean_mu += mu;
        mean_mu_squared += mu * mu;
        for (int d = 0; d < 3; ++d) mean_direction[d] += direction[d];
      }
    }
    const double samples = static_cast<double>(n_mu * n_phi);
    mean_mu /= samples;
    mean_mu_squared /= samples;
    for (int d = 0; d < 3; ++d) mean_direction[d] /= samples;

    // p(mu)=2mu gives <mu>=2/3 and <mu^2>=1/2; p(mu)=1 gives
    // <mu>=1/2 and <mu^2>=1/3.  Stratified midpoint samples remove random
    // seed uncertainty while retaining the production inverse-CDF mapping.
    AMPS_EXPECT_NEAR(mean_mu, analytic.mean_mu, 4.0e-5);
    AMPS_EXPECT_NEAR(mean_mu_squared, analytic.mean_mu_squared, 2.0e-6);
    for (int d = 0; d < 3; ++d)
      AMPS_EXPECT_NEAR(mean_direction[d],
                       analytic.mean_mu * unit_normal[d], 4.0e-5);
  }
}

AMPS_TEST_CASE("SPUTTERING-SOURCE-RATE-001",
               "incident flux, area, yield, and statistical weight close event rate") {
  // Independent dimensional oracle:
  // (2.5e12 projectiles m^-2 s^-1)(4 m^2)(3 products/projectile)
  // =3e13 products/s.  With W=2e8 products/event, lambda=1.5e5 events/s,
  // hence 3e5 expected events in two seconds.
  const Sputtering::Source::SourceRateAccounting accounting =
      Sputtering::Source::ComputeSourceRateAccounting(
          2.5e12, 4.0, 3.0, 2.0e8, 2.0);
  AMPS_EXPECT_NEAR(accounting.incident_projectiles_per_s, 1.0e13, 0.0);
  AMPS_EXPECT_NEAR(accounting.physical_products_per_s, 3.0e13, 0.0);
  AMPS_EXPECT_NEAR(accounting.physical_products_per_numerical_event,
                   2.0e8, 0.0);
  AMPS_EXPECT_NEAR(accounting.numerical_events_per_s, 1.5e5, 0.0);
  AMPS_EXPECT_NEAR(accounting.expected_numerical_events, 3.0e5, 0.0);
  AMPS_EXPECT_NEAR(Sputtering::Source::EffectivePhysicalWeight(4.0e8, 0.5),
                   2.0e8, 0.0);

  ExpectInvalid([]() {
    (void)Sputtering::Source::ComputeSourceRateAccounting(
        -1.0, 1.0, 1.0, 1.0, 1.0);
  });
  ExpectInvalid([]() {
    (void)Sputtering::Source::ComputeSourceRateAccounting(
        1.0, 1.0, 1.0, 0.0, 1.0);
  });
  ExpectInvalid([]() {
    (void)Sputtering::Source::EffectivePhysicalWeight(1.0, 0.0);
  });
  ExpectOverflow([]() {
    (void)Sputtering::Source::ComputeSourceRateAccounting(
        std::numeric_limits<double>::max(), 2.0, 1.0, 1.0, 1.0);
  });
}

AMPS_TEST_CASE("SPUTTERING-PRODUCT-LEDGER-001",
               "product identity, mass, momentum, energy, and rate ledgers close") {
  const int species = _O2_SPEC_;
  const double mass = 5.0e-26;
  Sputtering::Source::ProductLedger ledger(species);
  ledger.RecordEvent(species, mass,
                     std::array<double, 3>{{100.0, 0.0, 0.0}}, 100.0, 0.5);
  ledger.RecordEvent(species, mass,
                     std::array<double, 3>{{-20.0, 30.0, 40.0}}, 25.0, 2.0);

  // Each event represents 50 physical molecules.  The independent extensive
  // sums are m*50[(100,0,0)+(-20,30,40)] and
  // (m*50/2)[100^2+(-20)^2+30^2+40^2].
  AMPS_EXPECT_EQ(ledger.numerical_events(), static_cast<std::uint64_t>(2));
  AMPS_EXPECT_NEAR(ledger.physical_products(), 100.0, 0.0);
  AMPS_EXPECT_NEAR(ledger.mass_kg(), 100.0 * mass, 1.0e-38);
  AMPS_EXPECT_NEAR(ledger.momentum_kg_m_per_s()[0], mass * 4000.0, 1.0e-35);
  AMPS_EXPECT_NEAR(ledger.momentum_kg_m_per_s()[1], mass * 1500.0, 1.0e-35);
  AMPS_EXPECT_NEAR(ledger.momentum_kg_m_per_s()[2], mass * 2000.0, 1.0e-35);
  AMPS_EXPECT_NEAR(ledger.kinetic_energy_joule(), mass * 322500.0, 1.0e-33);

  const Sputtering::Source::ProductRates rates = ledger.Rates(2.0);
  AMPS_EXPECT_EQ(rates.product_species, species);
  AMPS_EXPECT_NEAR(rates.physical_products_per_s, 50.0, 0.0);
  AMPS_EXPECT_NEAR(rates.numerical_events_per_s, 1.0, 0.0);
  AMPS_EXPECT_NEAR(rates.mass_kg_per_s, 50.0 * mass, 1.0e-38);
  AMPS_EXPECT_NEAR(rates.momentum_newton[0], mass * 2000.0, 1.0e-35);
  AMPS_EXPECT_NEAR(rates.kinetic_power_watt, mass * 161250.0, 1.0e-33);

  const double products_before = ledger.physical_products();
  ExpectInvalid([&ledger, mass, species]() {
    ledger.RecordEvent(species + 1, mass,
                       std::array<double, 3>{{1.0, 0.0, 0.0}}, 1.0);
  });
  ExpectInvalid([&ledger, mass, species]() {
    ledger.RecordEvent(species, mass,
        std::array<double, 3>{{std::numeric_limits<double>::infinity(),
                              0.0, 0.0}}, 1.0);
  });
  AMPS_EXPECT_NEAR(ledger.physical_products(), products_before, 0.0);
}

AMPS_TEST_CASE("SPUTTERING-RATE-LEDGER-CLOSURE-001",
               "angular samples close planned product, mass, momentum, and energy rates") {
  using Sputtering::Source::EmissionAngleLaw;
  const double flux = 2.0e10;       // projectiles m^-2 s^-1
  const double area = 0.4;          // m^2
  const double yield = 0.8;         // products/projectile
  const double duration = 1.0;      // s
  const double mass = 32.0 * kAtomicMassKilogramOracle;
  const double speed = 500.0;       // m/s
  const int n_mu = 200;
  const int n_phi = 128;
  const int event_count = n_mu * n_phi;
  const double product_rate = flux * area * yield;
  const double event_weight = product_rate * duration / event_count;
  const std::array<double, 3> normal = {{0.0, 0.0, 1.0}};
  const EmissionAngleLaw laws[2] = {
      EmissionAngleLaw::KnudsenCosine,
      EmissionAngleLaw::UniformHemisphere};

  for (int law_index = 0; law_index < 2; ++law_index) {
    const EmissionAngleLaw law = laws[law_index];
    const Sputtering::Source::SourceRateAccounting plan =
        Sputtering::Source::ComputeSourceRateAccounting(
            flux, area, yield, event_weight, duration);
    AMPS_EXPECT_NEAR(plan.expected_numerical_events, event_count, 2.0e-11);

    Sputtering::Source::ProductLedger ledger(_O2_SPEC_);
    for (int i = 0; i < n_mu; ++i) {
      for (int j = 0; j < n_phi; ++j) {
        const std::array<double, 3> direction =
            Sputtering::Source::DirectionFromUnitSamples(
                normal, (i + 0.5) / n_mu, (j + 0.5) / n_phi, law);
        const std::array<double, 3> velocity = {{
            speed * direction[0], speed * direction[1],
            speed * direction[2]}};
        ledger.RecordEvent(_O2_SPEC_, mass, velocity, event_weight);
      }
    }
    const Sputtering::Source::ProductRates rates = ledger.Rates(duration);
    const double mean_mu = Sputtering::Source::GetAngularMoments(law).mean_mu;
    ExpectRelativeNear(rates.physical_products_per_s, product_rate, 2.0e-13);
    ExpectRelativeNear(rates.mass_kg_per_s, product_rate * mass, 2.0e-13);
    AMPS_EXPECT_NEAR(rates.momentum_newton[0], 0.0, 2.0e-26);
    AMPS_EXPECT_NEAR(rates.momentum_newton[1], 0.0, 2.0e-26);
    ExpectRelativeNear(rates.momentum_newton[2],
        product_rate * mass * speed * mean_mu, 6.0e-5);
    ExpectRelativeNear(rates.kinetic_power_watt,
        product_rate * 0.5 * mass * speed * speed, 2.0e-13);
  }
}
