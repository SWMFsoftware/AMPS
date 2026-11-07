#include "test_harness.h"

#include "He.h"
#include "Na.h"
#include "constants.h"

#include <cmath>
#include <limits>
#include <stdexcept>

namespace {
void ExpectInvalid(double (*function)(double, double), double velocity,
                   double distance) {
  bool rejected = false;
  try { (void)function(velocity, distance); }
  catch (const std::invalid_argument&) { rejected = true; }
  AMPS_EXPECT_TRUE(rejected);
}


void ExpectOverflow(double (*function)(double, double), double velocity,
                    double distance) {
  bool rejected = false;
  try { (void)function(velocity, distance); }
  catch (const std::overflow_error&) { rejected = true; }
  AMPS_EXPECT_TRUE(rejected);
}
}

AMPS_TEST_CASE("SPECIES-RADIATION-REFERENCE-001",
               "helium and sodium tables reproduce retained digitized regression points") {
  const double combi_distance = 149598000000.0;
  const double killen_distance = 0.352 * _AU_;
  AMPS_EXPECT_NEAR(SodiumGfactor__D1D2__Combi_1997_icarus(
      -50000.0, combi_distance), 14.55069, 1.0e-13);
  // 0.5 km/s is halfway between independently tabulated values 0.80269 and
  // 0.85768 at 0 and 1 km/s, respectively.
  AMPS_EXPECT_NEAR(SodiumGfactor__D1D2__Combi_1997_icarus(
      500.0, combi_distance), 0.830185, 1.0e-13);
  AMPS_EXPECT_NEAR(SodiumGfactor__5891_58A__Killen_2009_AJSS(
      0.0, killen_distance), 4.24193, 1.0e-13);
  AMPS_EXPECT_NEAR(SodiumGfactor__5897_56A__Killen_2009_AJSS(
      0.0, killen_distance), 2.62726, 1.0e-13);
  AMPS_EXPECT_NEAR(HeliumGfactor__584_35A__Killen_2009_AJSS(
      0.0, killen_distance), 4.64032e-5, 1.0e-17);
}

AMPS_TEST_CASE("SPECIES-RADIATION-SCALING-001",
               "all g-factors obey inverse-square distance scaling") {
  const double base = 0.352 * _AU_;
  const double velocity = 3000.0;
  const double he = HeliumGfactor__584_35A__Killen_2009_AJSS(velocity, base);
  const double d1 = SodiumGfactor__5891_58A__Killen_2009_AJSS(velocity, base);
  const double d2 = SodiumGfactor__5897_56A__Killen_2009_AJSS(velocity, base);
  const double combi = SodiumGfactor__D1D2__Combi_1997_icarus(velocity, base);
  AMPS_EXPECT_NEAR(HeliumGfactor__584_35A__Killen_2009_AJSS(
      velocity, 2.0 * base), 0.25 * he, 1.0e-18);
  AMPS_EXPECT_NEAR(SodiumGfactor__5891_58A__Killen_2009_AJSS(
      velocity, 2.0 * base), 0.25 * d1, 1.0e-13);
  AMPS_EXPECT_NEAR(SodiumGfactor__5897_56A__Killen_2009_AJSS(
      velocity, 2.0 * base), 0.25 * d2, 1.0e-13);
  AMPS_EXPECT_NEAR(SodiumGfactor__D1D2__Combi_1997_icarus(
      velocity, 2.0 * base), 0.25 * combi, 1.0e-13);
}

AMPS_TEST_CASE("SPECIES-RADIATION-ACCELERATION-001",
               "radiation accelerations retain published units and photon momentum") {
  const double combi_distance = 149598000000.0;
  // Combi boundary values are cm/s^2; the production result is m/s^2.
  AMPS_EXPECT_NEAR(SodiumRadiationPressureAcceleration__Combi_1997_icarus(
      -60000.0, combi_distance), 0.464855, 1.0e-14);
  AMPS_EXPECT_NEAR(SodiumRadiationPressureAcceleration__Combi_1997_icarus(
      60000.0, combi_distance), 0.462883, 1.0e-14);
  const double distance = 0.352 * _AU_;
  // Use independent SI constants and the retained 0 km/s digitized g-factor,
  // rather than calling the production g-factor or reusing constants.h in the
  // oracle. h and c are exact SI definitions; the atomic-mass value is the
  // CODATA 2018 value used by the current AMPS constants table.
  const double planck_si = 6.62607015e-34;
  const double speed_of_light_si = 299792458.0;
  const double atomic_mass_si = 1.66053906660e-27;
  const double digitized_g = 4.64032e-5;
  const double expected = planck_si * 5.1303578e15 * digitized_g /
                          (1.00794 * atomic_mass_si * speed_of_light_si);
  // AMPS retains an older CODATA atomic-mass constant, so the independent
  // modern-CODATA oracle uses a 2e-7 relative compatibility tolerance.
  AMPS_EXPECT_NEAR(RadiationPressureAcceleration__Chamberlian_TOPA_Eq7_2_9(
      0.0, distance), expected, 2.0e-7 * std::fabs(expected));
}

AMPS_TEST_CASE("SPECIES-RADIATION-INVALID-001",
               "nonfinite velocity and nonpositive/nonfinite distance are rejected") {
  const double nan = std::numeric_limits<double>::quiet_NaN();
  const double infinity = std::numeric_limits<double>::infinity();
  ExpectInvalid(HeliumGfactor__584_35A__Killen_2009_AJSS, nan, _AU_);
  ExpectInvalid(RadiationPressureAcceleration__Chamberlian_TOPA_Eq7_2_9,
                0.0, 0.0);
  ExpectInvalid(SodiumGfactor__D1D2__Combi_1997_icarus, 0.0, -_AU_);
  ExpectInvalid(SodiumRadiationPressureAcceleration__Combi_1997_icarus,
                0.0, infinity);
  ExpectInvalid(SodiumGfactor__5891_58A__Killen_2009_AJSS, infinity, _AU_);
  ExpectInvalid(SodiumGfactor__5897_56A__Killen_2009_AJSS, 0.0, nan);

  // Finite velocities outside every table must clamp before index conversion;
  // DBL_MAX-to-int was undefined in the former Killen lookup implementation.
  const double largest = std::numeric_limits<double>::max();
  const double killen_distance = 0.352 * _AU_;
  AMPS_EXPECT_NEAR(HeliumGfactor__584_35A__Killen_2009_AJSS(
      -largest, killen_distance), 3.68739e-5, 1.0e-17);
  AMPS_EXPECT_NEAR(HeliumGfactor__584_35A__Killen_2009_AJSS(
      largest, killen_distance), 5.03423e-5, 1.0e-17);
  AMPS_EXPECT_NEAR(SodiumGfactor__5891_58A__Killen_2009_AJSS(
      -largest, killen_distance), 42.33717, 1.0e-13);
  AMPS_EXPECT_NEAR(SodiumGfactor__5891_58A__Killen_2009_AJSS(
      largest, killen_distance), 42.06349, 1.0e-13);
  AMPS_EXPECT_NEAR(SodiumGfactor__5897_56A__Killen_2009_AJSS(
      -largest, killen_distance), 26.13574, 1.0e-13);
  AMPS_EXPECT_NEAR(SodiumGfactor__5897_56A__Killen_2009_AJSS(
      largest, killen_distance), 25.12315, 1.0e-13);

  // Tiny positive distance is in-domain geometrically but its inverse-square
  // value is not representable as a double; each public path must fail with a
  // typed range error instead of returning infinity.
  const double tiny = std::numeric_limits<double>::denorm_min();
  ExpectOverflow(HeliumGfactor__584_35A__Killen_2009_AJSS, 0.0, tiny);
  ExpectOverflow(RadiationPressureAcceleration__Chamberlian_TOPA_Eq7_2_9,
                 0.0, tiny);
  ExpectOverflow(SodiumGfactor__D1D2__Combi_1997_icarus, 0.0, tiny);
  ExpectOverflow(SodiumRadiationPressureAcceleration__Combi_1997_icarus,
                 0.0, tiny);
  ExpectOverflow(SodiumGfactor__5891_58A__Killen_2009_AJSS, 0.0, tiny);
  ExpectOverflow(SodiumGfactor__5897_56A__Killen_2009_AJSS, 0.0, tiny);
}
