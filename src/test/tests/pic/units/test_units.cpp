#include "test_harness.h"
#include "pic_units_normalization.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>

namespace {

void ExpectRelative(double value, double reference, double tolerance) {
  AMPS_EXPECT_TRUE(std::isfinite(value));
  AMPS_EXPECT_TRUE(std::isfinite(reference));
  const double scale = std::max(std::fabs(reference), 1.0e-300);
  AMPS_EXPECT_TRUE(std::fabs(value - reference) <= tolerance * scale);
}

picunits::NormScalesSI TestScales() {
  // Deliberately non-decimal-related scales make missing powers of 10 in the
  // SI-to-CGS path visible instead of allowing them to cancel accidentally.
  return picunits::NormScalesSI{2.5e6, 7.5e4, 3.25e-21};
}

void ExpectInvalid(const picunits::NormScalesSI& scales) {
  bool rejected = false;
  try {
    (void)picunits::build(scales);
  } catch (const std::invalid_argument& error) {
    rejected = true;
    AMPS_EXPECT_TRUE(std::string(error.what()).find("finite and positive") !=
                     std::string::npos);
  }
  AMPS_EXPECT_TRUE(rejected);
}

void ExpectUnrepresentable(const picunits::NormScalesSI& scales) {
  bool rejected = false;
  try {
    (void)picunits::build(scales);
  } catch (const std::overflow_error& error) {
    rejected = true;
    AMPS_EXPECT_TRUE(std::string(error.what()).find("not representable") !=
                     std::string::npos);
  }
  AMPS_EXPECT_TRUE(rejected);
}

template <class Callable>
void ExpectInvalidScaleParameter(Callable callable) {
  bool rejected = false;
  try {
    callable();
  } catch (const std::invalid_argument& error) {
    rejected = true;
    AMPS_EXPECT_TRUE(std::string(error.what()).find("finite and positive") !=
                     std::string::npos);
  }
  AMPS_EXPECT_TRUE(rejected);
}

template <class Callable>
void ExpectScaleOverflow(Callable callable) {
  bool rejected = false;
  try {
    callable();
  } catch (const std::overflow_error& error) {
    rejected = true;
    AMPS_EXPECT_TRUE(std::string(error.what()).find("not representable") !=
                     std::string::npos);
  }
  AMPS_EXPECT_TRUE(rejected);
}

}  // namespace

AMPS_TEST_CASE("UNITS-FACTOR-DERIVATION-001",
               "SI base scales reproduce independently derived CGS factors") {
  const picunits::NormScalesSI scales = TestScales();
  const picunits::Factors factors = picunits::build(scales);

  // These equations are evaluated here from the dimensional definitions, not
  // copied from a Factors instance. They independently pin every power of the
  // primary length, speed, and mass scales.
  const double length_cgs = 100.0 * scales.lSI_m;
  const double speed_cgs = 100.0 * scales.uSI_mps;
  const double mass_cgs = 1000.0 * scales.mSI_kg;
  const double density_cgs = mass_cgs /
      (length_cgs * length_cgs * length_cgs);
  const double magnetic_cgs = std::sqrt(density_cgs) * speed_cgs;
  const double charge_cgs = std::sqrt(mass_cgs * length_cgs) * speed_cgs;
  const double pressure_cgs = density_cgs * speed_cgs * speed_cgs;
  const double current_cgs = charge_cgs * speed_cgs /
      (length_cgs * length_cgs * length_cgs);

  ExpectRelative(factors.L0_cm, length_cgs, 2.0e-15);
  ExpectRelative(factors.U0_cms, speed_cgs, 2.0e-15);
  ExpectRelative(factors.M0_g, mass_cgs, 2.0e-15);
  ExpectRelative(factors.rho0, density_cgs, 3.0e-15);
  ExpectRelative(factors.B0_G, magnetic_cgs, 3.0e-15);
  ExpectRelative(factors.Q0, charge_cgs, 3.0e-15);
  ExpectRelative(factors.P0, pressure_cgs, 3.0e-15);
  ExpectRelative(factors.J0, current_cgs, 4.0e-15);
  ExpectRelative(factors.T0_s, scales.lSI_m / scales.uSI_mps, 2.0e-15);
  ExpectRelative(factors.Si2NoRho, 1.0e-3 / density_cgs, 3.0e-15);
  ExpectRelative(factors.Si2NoV, 100.0 / speed_cgs, 2.0e-15);
  ExpectRelative(factors.Si2NoB, 1.0e4 / magnetic_cgs, 3.0e-15);
  ExpectRelative(factors.Si2NoP, 10.0 / pressure_cgs, 3.0e-15);
  ExpectRelative(factors.Si2NoJ,
                 (1.0e-5 * speed_cgs) / current_cgs, 4.0e-15);
  ExpectRelative(factors.Si2NoE,
                 (1.0e6 / speed_cgs) / magnetic_cgs, 4.0e-15);
}

AMPS_TEST_CASE("UNITS-SCALAR-ROUNDTRIP-001",
               "scalar and number-density helpers preserve SI values") {
  const picunits::Factors factors = picunits::build(TestScales());
  const double values[10] = {4.25e-18, -8.75e4, 6.5e-9, 2.1e-7, -3.2e-11,
                             9.4e-5, 7.3e8, 12.75, 2.8e-26, -1.602176634e-19};
  const double normalized[10] = {
      picunits::si2no_rho(values[0], factors),
      picunits::si2no_v(values[1], factors),
      picunits::si2no_B(values[2], factors),
      picunits::si2no_P(values[3], factors),
      picunits::si2no_J(values[4], factors),
      picunits::si2no_E(values[5], factors),
      picunits::si2no_L(values[6], factors),
      picunits::si2no_t(values[7], factors),
      picunits::si2no_m(values[8], factors),
      picunits::si2no_q(values[9], factors)};
  const double recovered[10] = {
      picunits::no2si_rho(normalized[0], factors),
      picunits::no2si_v(normalized[1], factors),
      picunits::no2si_B(normalized[2], factors),
      picunits::no2si_P(normalized[3], factors),
      picunits::no2si_J(normalized[4], factors),
      picunits::no2si_E(normalized[5], factors),
      picunits::no2si_L(normalized[6], factors),
      picunits::no2si_t(normalized[7], factors),
      picunits::no2si_m(normalized[8], factors),
      picunits::no2si_q(normalized[9], factors)};
  for (int quantity = 0; quantity < 10; ++quantity)
    ExpectRelative(recovered[quantity], values[quantity], 4.0e-15);

  const double density_si = 4.5e7;  // m^-3
  const double reference_cm3 = 1250.0;
  const double explicit_no = picunits::si2no_n(density_si, reference_cm3);
  ExpectRelative(explicit_no, 0.036, 2.0e-15);
  ExpectRelative(picunits::no2si_n(explicit_no, reference_cm3), density_si,
                 2.0e-15);
  const double species_mass = 2.65686250656e-26;  // kg, neutral O reference
  const double derived_no = picunits::si2no_n(density_si, factors, species_mass);
  ExpectRelative(picunits::no2si_n(derived_no, factors, species_mass), density_si,
                 3.0e-15);
}

AMPS_TEST_CASE("UNITS-VECTOR-ROUNDTRIP-001",
               "array and pointer vector helpers preserve signs and alias safely") {
  const picunits::Factors factors = picunits::build(TestScales());
  const picunits::Vec3 source = {{-3.25e4, 0.0, 8.75e4}};
  const picunits::Vec3 velocity =
      picunits::no2si_v3(picunits::si2no_v3(source, factors), factors);
  const picunits::Vec3 magnetic =
      picunits::no2si_B3(picunits::si2no_B3(source, factors), factors);
  const picunits::Vec3 electric =
      picunits::no2si_E3(picunits::si2no_E3(source, factors), factors);
  const picunits::Vec3 current =
      picunits::no2si_J3(picunits::si2no_J3(source, factors), factors);
  const picunits::Vec3 length =
      picunits::no2si_L3(picunits::si2no_L3(source, factors), factors);
  for (int dimension = 0; dimension < 3; ++dimension) {
    ExpectRelative(velocity[dimension], source[dimension], 3.0e-15);
    ExpectRelative(magnetic[dimension], source[dimension], 3.0e-15);
    ExpectRelative(electric[dimension], source[dimension], 3.0e-15);
    ExpectRelative(current[dimension], source[dimension], 3.0e-15);
    ExpectRelative(length[dimension], source[dimension], 3.0e-15);
  }

  // Raw-pointer helpers explicitly promise in-place operation. A nonuniform
  // vector catches overwrite-before-read and component-index regressions.
  double in_place[3] = {1.25e-7, -9.5e-8, 3.75e-9};
  const double original[3] = {in_place[0], in_place[1], in_place[2]};
  picunits::si2no_B3(in_place, in_place, factors);
  picunits::no2si_B3(in_place, in_place, factors);
  for (int dimension = 0; dimension < 3; ++dimension)
    ExpectRelative(in_place[dimension], original[dimension], 3.0e-15);
}

AMPS_TEST_CASE("UNITS-DIMENSIONAL-IDENTITIES-001",
               "normalizers preserve time, pressure, Lorentz, and density identities") {
  const picunits::NormScalesSI scales = TestScales();
  const picunits::Factors factors = picunits::build(scales);
  const double inverse_products[] = {
      factors.Si2NoRho * factors.No2SiRho,
      factors.Si2NoV * factors.No2SiV,
      factors.Si2NoB * factors.No2SiB,
      factors.Si2NoP * factors.No2SiP,
      factors.Si2NoJ * factors.No2SiJ,
      factors.Si2NoE * factors.No2SiE,
      factors.Si2NoL * factors.No2SiL,
      factors.Si2NoT * factors.No2SiT,
      factors.Si2NoM * factors.No2SiM,
      factors.Si2NoQ * factors.No2SiQ};
  for (std::size_t index = 0;
       index < sizeof(inverse_products) / sizeof(inverse_products[0]); ++index)
    ExpectRelative(inverse_products[index], 1.0, 3.0e-15);

  ExpectRelative(factors.Si2NoT,
                 factors.Si2NoL / factors.Si2NoV, 2.0e-15);
  ExpectRelative(picunits::si2no_L(scales.lSI_m, factors), 1.0, 2.0e-15);
  ExpectRelative(picunits::si2no_v(scales.uSI_mps, factors), 1.0, 2.0e-15);
  ExpectRelative(picunits::si2no_t(scales.lSI_m / scales.uSI_mps, factors),
                 1.0, 2.0e-15);
  ExpectRelative(picunits::si2no_m(scales.mSI_kg, factors), 1.0, 2.0e-15);

  const double magnetic_si = -7.25e-9;
  const double convective_e_si = scales.uSI_mps * magnetic_si;
  ExpectRelative(picunits::si2no_E(convective_e_si, factors),
                 picunits::si2no_v(scales.uSI_mps, factors) *
                     picunits::si2no_B(magnetic_si, factors),
                 4.0e-15);

  const double reference_density_si = factors.rho0 * 1.0e3;  // g/cm3 -> kg/m3
  const double pressure_si = reference_density_si *
                             scales.uSI_mps * scales.uSI_mps;
  ExpectRelative(picunits::si2no_P(pressure_si, factors), 1.0, 4.0e-15);

  const double species_mass = 2.991507e-26;
  const double number_density = 8.2e11;
  ExpectRelative(picunits::si2no_n(number_density, factors, species_mass),
                 picunits::si2no_rho(number_density * species_mass, factors),
                 4.0e-15);

  const double charge = 1.602176634e-19;
  const double particle_mass = 1.67262192369e-27;
  ExpectRelative(picunits::si2no_q(charge, factors) /
                     picunits::si2no_m(particle_mass, factors),
                 (charge / particle_mass) * factors.No2SiB * factors.No2SiT,
                 4.0e-15);

  // Independently reproduce the scale choice used by the production AMU/e
  // initializer. This both checks the dimensional derivation and proves the
  // deliberately microscopic length scale remains representable after the
  // stricter build-domain validation.
  const double atomic_mass_si = 1.66053906660e-27;
  const double elementary_charge_si = 1.602176634e-19;
  const double amu_length = 1.0e-7 * elementary_charge_si *
                            elementary_charge_si / atomic_mass_si;
  const picunits::Factors amu_factors = picunits::build(
      picunits::NormScalesSI{amu_length, 1.0e6, atomic_mass_si});
  ExpectRelative(amu_factors.No2SiM, atomic_mass_si, 3.0e-15);
  ExpectRelative(amu_factors.No2SiQ, elementary_charge_si, 4.0e-15);
}

AMPS_TEST_CASE("UNITS-INVALID-SCALES-001",
               "invalid and unrepresentable primary/reference scales are rejected") {
  const picunits::NormScalesSI valid = TestScales();
  const double nan = std::numeric_limits<double>::quiet_NaN();
  const double infinity = std::numeric_limits<double>::infinity();
  const double maximum = std::numeric_limits<double>::max();
  const double tiny = std::numeric_limits<double>::denorm_min();
  ExpectInvalid(picunits::NormScalesSI{0.0, valid.uSI_mps, valid.mSI_kg});
  ExpectInvalid(picunits::NormScalesSI{valid.lSI_m, -1.0, valid.mSI_kg});
  ExpectInvalid(picunits::NormScalesSI{valid.lSI_m, valid.uSI_mps, 0.0});
  ExpectInvalid(picunits::NormScalesSI{nan, valid.uSI_mps, valid.mSI_kg});
  ExpectInvalid(picunits::NormScalesSI{valid.lSI_m, infinity, valid.mSI_kg});
  ExpectInvalid(picunits::NormScalesSI{valid.lSI_m, valid.uSI_mps, infinity});

  // These inputs satisfy the finite-positive precondition, but their CGS
  // bases or derived normalizers cannot be represented by a nonzero finite
  // double. They must fail rather than returning a partially valid Factors.
  ExpectUnrepresentable(
      picunits::NormScalesSI{maximum, valid.uSI_mps, valid.mSI_kg});
  ExpectUnrepresentable(
      picunits::NormScalesSI{valid.lSI_m, maximum, valid.mSI_kg});
  ExpectUnrepresentable(
      picunits::NormScalesSI{valid.lSI_m, valid.uSI_mps, maximum});
  ExpectUnrepresentable(
      picunits::NormScalesSI{tiny, valid.uSI_mps, valid.mSI_kg});
  ExpectUnrepresentable(
      picunits::NormScalesSI{valid.lSI_m, tiny, valid.mSI_kg});
  ExpectUnrepresentable(
      picunits::NormScalesSI{valid.lSI_m, valid.uSI_mps, tiny});

  const picunits::Factors factors = picunits::build(valid);
  const double invalid_reference[4] = {0.0, -1.0, nan, infinity};
  for (int i = 0; i < 4; ++i) {
    ExpectInvalidScaleParameter([&invalid_reference, i]() {
      (void)picunits::si2no_n(1.0, invalid_reference[i]);
    });
    ExpectInvalidScaleParameter([&invalid_reference, i]() {
      (void)picunits::no2si_n(1.0, invalid_reference[i]);
    });
    ExpectInvalidScaleParameter([&factors, &invalid_reference, i]() {
      (void)picunits::si2no_n(1.0, factors, invalid_reference[i]);
    });
    ExpectInvalidScaleParameter([&factors, &invalid_reference, i]() {
      (void)picunits::no2si_n(1.0, factors, invalid_reference[i]);
    });
  }

  // A finite positive kg value can still overflow the kg-to-g scale used to
  // derive the reference number density; this has a distinct typed failure.
  ExpectScaleOverflow([&factors, maximum]() {
    (void)picunits::si2no_n(1.0, factors, maximum);
  });
  ExpectScaleOverflow([&factors, maximum]() {
    (void)picunits::no2si_n(1.0, factors, maximum);
  });
}
