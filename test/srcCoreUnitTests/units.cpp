#include "pic.h"
#include "src/test/native/test_harness.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>

namespace {

// PIC::Units::Factors is process-global production state.  Native cases can
// run in normal, reverse, or shuffled order, so every case restores the exact
// state it observed and cannot make another case order-dependent.
class UnitFactorsScope {
 public:
  UnitFactorsScope() : saved_(PIC::Units::Factors) {
    PIC::Units::Factors = picunits::Factors{};
  }

  ~UnitFactorsScope() { PIC::Units::Factors = saved_; }

 private:
  picunits::Factors saved_;
};

void ExpectRelative(double value, double reference, double tolerance) {
  AMPS_EXPECT_TRUE(std::isfinite(value));
  AMPS_EXPECT_TRUE(std::isfinite(reference));
  const double scale = std::max(std::fabs(reference), 1.0e-300);
  AMPS_EXPECT_TRUE(std::fabs(value-reference)<=tolerance*scale);
}

}  // namespace

namespace AMPS {
namespace CoreTests {

// The application is archived before final linkage.  This explicit anchor
// ensures this translation unit and its static test registrations are loaded.
void LinkUnitTests() {}

}  // namespace CoreTests
}  // namespace AMPS

AMPS_TEST_CASE("PIC-UNITS-GLOBAL-INIT-001",
               "global AMU/charge initialization is accurate and idempotent") {
  UnitFactorsScope state;
  PIC::Units::InitializeAMUChargeNormalization(false);
  const picunits::Factors first = PIC::Units::Factors;

  // These are independent SI reference values.  AMPS uses a historical AMU
  // and elementary-charge literal, so the tolerances bound those known source
  // differences while still detecting a broken normalization derivation.
  const double atomic_mass_si = 1.66053906660e-27;
  const double elementary_charge_si = 1.602176634e-19;
  ExpectRelative(first.No2SiM, atomic_mass_si, 5.0e-12);
  ExpectRelative(first.No2SiQ, elementary_charge_si, 5.0e-8);
  ExpectRelative(first.Si2NoM*first.No2SiM, 1.0, 3.0e-15);
  ExpectRelative(first.Si2NoQ*first.No2SiQ, 1.0, 3.0e-15);

  // A repeated startup call must not alter the already committed global unit
  // system.  Exact equality is appropriate because the implementation returns
  // before recomputing or assigning any factor.
  PIC::Units::InitializeAMUChargeNormalization(false);
  AMPS_EXPECT_NEAR(PIC::Units::Factors.L0_cm, first.L0_cm, 0.0);
  AMPS_EXPECT_NEAR(PIC::Units::Factors.U0_cms, first.U0_cms, 0.0);
  AMPS_EXPECT_NEAR(PIC::Units::Factors.M0_g, first.M0_g, 0.0);
  AMPS_EXPECT_NEAR(PIC::Units::Factors.Si2NoM, first.Si2NoM, 0.0);
  AMPS_EXPECT_NEAR(PIC::Units::Factors.No2SiQ, first.No2SiQ, 0.0);
}

AMPS_TEST_CASE("PIC-UNITS-CONVERSION-TABLE-001",
               "global conversion-table wrapper prints active production factors") {
  UnitFactorsScope state;
  PIC::Units::InitializeAMUChargeNormalization(false);

  FILE* output = std::tmpfile();
  AMPS_EXPECT_TRUE(output!=NULL);
  PIC::Units::PrintConversionTable(output);
  AMPS_EXPECT_EQ(std::fflush(output), 0);
  AMPS_EXPECT_EQ(std::fseek(output, 0, SEEK_END), 0);
  const long length = std::ftell(output);
  AMPS_EXPECT_TRUE(length>0);
  AMPS_EXPECT_EQ(std::fseek(output, 0, SEEK_SET), 0);
  std::string text(static_cast<std::size_t>(length), '\0');
  AMPS_EXPECT_EQ(std::fread(&text[0], 1, text.size(), output), text.size());
  std::fclose(output);

  AMPS_EXPECT_TRUE(text.find("PIC / AMPS Unit Conversion Coefficients")!=
                   std::string::npos);
  AMPS_EXPECT_TRUE(text.find("Particle charge")!=std::string::npos);
  AMPS_EXPECT_TRUE(text.find("Number density (p+)")!=std::string::npos);

  // The formatted scientific value binds the wrapper output to the active
  // global Factors object, rather than merely proving that a static heading
  // was printed.
  char expected_mass[64];
  std::snprintf(expected_mass, sizeof(expected_mass), "%.10e",
                PIC::Units::Factors.No2SiM);
  AMPS_EXPECT_TRUE(text.find(expected_mass)!=std::string::npos);
}
