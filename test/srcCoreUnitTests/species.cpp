#include "pic.h"
#include "src/test/native/test_harness.h"

#include <algorithm>
#include <cmath>
#include <cstring>

namespace {

void ExpectRelative(double value, double reference, double tolerance) {
  AMPS_EXPECT_TRUE(std::isfinite(value));
  const double scale = std::max(std::fabs(reference), 1.0e-300);
  AMPS_EXPECT_TRUE(std::fabs(value - reference) <= tolerance * scale);
}

struct ExpectedSpecies {
  const char* symbol;
  int generated_index;
  double mass_in_amu;
  double charge_in_elementary_units;
};

/* Literal oracle for input/test/core-unit-tests.input. The first seven
 * entries preserve the neutral chemistry/IDF order. H+ and the electron were
 * appended for T4 charged-mover tests; including them here makes a future
 * SpeciesList change fail the generated-table contract instead of silently
 * reducing this test to a neutral prefix. Masses are deck values in atomic
 * mass units, and charges are dimensionless multiples of +e. */
const ExpectedSpecies kExpectedSpecies[] = {
    {"H2O", _H2O_SPEC_, 18.0, 0.0},
    {"OH", _OH_SPEC_, 17.0, 0.0},
    {"H", _H_SPEC_, 1.0, 0.0},
    {"O", _O_SPEC_, 16.0, 0.0},
    {"H2", _H2_SPEC_, 2.0, 0.0},
    {"O2", _O2_SPEC_, 32.0, 0.0},
    {"CO2", _CO2_SPEC_, 44.0, 0.0},
    {"H_PLUS", _H_PLUS_SPEC_, 1.0, 1.0},
    {"ELECTRON", _ELECTRON_SPEC_, 0.000548579909, -1.0},
};

}  // namespace

namespace AMPS {
namespace CoreTests {

void LinkSpeciesTests() {}

}  // namespace CoreTests
}  // namespace AMPS

AMPS_TEST_CASE("PIC-SPECIES-CONSTANTS-001",
               "fundamental constants agree with independent SI/CODATA values") {
  // The atomic-mass reference is the CODATA 2018 value in kg. Elementary
  // charge and Boltzmann references are exact in the post-2019 SI. AMPS keeps
  // older pre-redefinition literals for the latter two, so the tolerances
  // explicitly bound their known historical differences rather than deriving
  // expected values from production macros.
  const double atomic_mass_reference = 1.66053906660e-27;
  const double elementary_charge_si = 1.602176634e-19;
  const double boltzmann_si = 1.380649e-23;
  ExpectRelative(_AMU_, atomic_mass_reference, 5.0e-12);
  ExpectRelative(ElectronCharge, elementary_charge_si, 5.0e-8);
  ExpectRelative(Kbol, boltzmann_si, 1.0e-6);
  AMPS_EXPECT_NEAR(SpeedOfLight, 299792458.0, 0.0);
  ExpectRelative(Pi, 3.141592653589793238462643383279502884, 1.0e-15);
  // This legacy H2O constant is independently tabulated and intentionally not
  // used as the oracle for the generated deck's `18*_AMU_` mass expression.
  ExpectRelative(_H2O__MASS_, 2.99151e-26, 1.0e-15);
}

AMPS_TEST_CASE("PIC-SPECIES-GENERATED-MAPPING-001",
               "deck order, generated macros, and chemical symbols map one-to-one") {
  const int count = static_cast<int>(
      sizeof(kExpectedSpecies) / sizeof(kExpectedSpecies[0]));
  AMPS_EXPECT_EQ(PIC::nTotalSpecies, count);
  for (int index = 0; index < count; ++index) {
    // The array order is the literal SpeciesList order in
    // input/test/core-unit-tests.input, an oracle outside generated pic.h.
    AMPS_EXPECT_EQ(kExpectedSpecies[index].generated_index, index);
    AMPS_EXPECT_TRUE(std::strcmp(
        PIC::MolecularData::GetChemSymbol(index),
        kExpectedSpecies[index].symbol) == 0);
    char copied[_MAX_STRING_LENGTH_PIC_];
    PIC::MolecularData::GetChemSymbol(copied, index);
    AMPS_EXPECT_TRUE(std::strcmp(copied, kExpectedSpecies[index].symbol) == 0);
  }
}

AMPS_TEST_CASE("PIC-SPECIES-GENERATED-TABLES-001",
               "generated masses, charges, and gas types match the deck") {
  // CODATA 2018 atomic mass is repeated here as test data. Expected masses are
  // computed from the input deck's integer molecular weights, not AMPS' mass
  // macros or the generated MolMass table being tested.
  const double atomic_mass_reference = 1.66053906660e-27;
  const double elementary_charge_reference = 1.602176634e-19;
  const int count = static_cast<int>(
      sizeof(kExpectedSpecies) / sizeof(kExpectedSpecies[0]));
  for (int index = 0; index < count; ++index) {
    const double expected_mass =
        kExpectedSpecies[index].mass_in_amu * atomic_mass_reference;
    ExpectRelative(PIC::MolecularData::GetMass(index), expected_mass, 5.0e-12);
    const double expected_charge =
        kExpectedSpecies[index].charge_in_elementary_units *
        elementary_charge_reference;
    if (expected_charge == 0.0)
      AMPS_EXPECT_NEAR(
          PIC::MolecularData::GetElectricCharge(index), 0.0, 0.0);
    else
      // AMPS retains a pre-2019 elementary-charge literal, so use the same
      // independently documented tolerance as PIC-SPECIES-CONSTANTS-001.
      ExpectRelative(PIC::MolecularData::GetElectricCharge(index),
                     expected_charge, 5.0e-8);
    AMPS_EXPECT_EQ(PIC::MolecularData::GetSpecieType(index),
                   _PIC_SPECIE_TYPE__GAS_);
  }
}

AMPS_TEST_CASE("PIC-SPECIES-LOOKUP-DIAGNOSTICS-001",
               "compiled species are unique and unknown names return a typed miss") {
  const int count = static_cast<int>(
      sizeof(kExpectedSpecies) / sizeof(kExpectedSpecies[0]));
  for (int left = 0; left < count; ++left) {
    char known[_MAX_STRING_LENGTH_PIC_];
    std::strcpy(known, kExpectedSpecies[left].symbol);
    AMPS_EXPECT_EQ(PIC::MolecularData::GetSpecieNumber(known), left);
    for (int right = left + 1; right < count; ++right)
      AMPS_EXPECT_TRUE(std::strcmp(PIC::MolecularData::GetChemSymbol(left),
                                   PIC::MolecularData::GetChemSymbol(right)) != 0);
  }

  char unknown[] = "Xe_NOT_COMPILED";
  char wrong_case[] = "h2o";
  AMPS_EXPECT_EQ(PIC::MolecularData::GetSpecieNumber(unknown), -1);
  AMPS_EXPECT_EQ(PIC::MolecularData::GetSpecieNumber(wrong_case), -1);
}
