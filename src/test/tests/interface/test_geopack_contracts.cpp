#include "test_harness.h"
#include "GeopackInterface.h"

#include <cstring>
#include <stdexcept>

namespace {

double intercepted_position[3] = {0.0, 0.0, 0.0};

void ExpectIdentity(const double matrix[3][3]) {
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j)
      AMPS_EXPECT_NEAR(matrix[i][j], i == j ? 1.0 : 0.0, 0.0);
}

}  // namespace

/*
 * Intercept the exact legacy Fortran symbol called by production
 * Geopack::IGRF::GetMagneticField.  The stub is an observation boundary, not
 * an expected-value implementation: it records arguments and returns a fixed
 * non-symmetric GSM vector.  Hand calculations below remain independent.
 */
extern "C" void igrf_gsw_08_(double* x, double* y, double* z,
                              double* bx, double* by, double* bz) {
  intercepted_position[0] = *x;
  intercepted_position[1] = *y;
  intercepted_position[2] = *z;
  *bx = 7.0;
  *by = 11.0;
  *bz = 13.0;
}

/*
 * The production translation unit declares the complete legacy GEOPACK ABI.
 * Provide inert definitions for entry points that these bounded cases do not
 * call, so the test remains portable across linkers without relying on ELF
 * section garbage collection to discard unresolved functions.
 */
extern "C" void recalc_08_(int*, int*, int*, int*, int*,
                            double*, double*, double*) {}
extern "C" void igrf_prepare_context_08_(
    int*, int*, int*, int*, int*, double*, double*, double*,
    double*, double*, double*, double*, double*) {}
extern "C" void igrf_gsw_08_ctx_(
    const double*, const double*, const double*, const double*,
    double*, double*, double*, double*, double*, double*) {}
extern "C" void sphcar_08_(double*, double*, double*, double*, double*,
                            double*, int*) {}
extern "C" void bspcar_08_(double*, double*, double*, double*, double*,
                            double*, double*, double*) {}
extern "C" void igrf_geo_08_(double*, double*, double*, double*, double*,
                              double*) {}
extern "C" void gswgse_08_(double*, double*, double*, double*, double*,
                            double*, int*) {}

[[noreturn]] void exit(int, const char*, const char* message) {
  throw std::runtime_error(message ? message : "AMPS interface fatal diagnostic");
}

AMPS_TEST_CASE("INTERFACE-GEOPACK-MATRIX-001",
               "row-major Geopack matrix/vector argument ordering") {
  const double matrix[3][3] = {
      {2.0, -1.0, 4.0},
      {0.5, 3.0, -2.0},
      {-3.0, 5.0, 1.0}};
  const double input[3] = {7.0, -2.0, 0.5};
  double output[3] = {0.0, 0.0, 0.0};
  Geopack::MatrixVectorMultiply(matrix, input, output);
  AMPS_EXPECT_NEAR(output[0], 18.0, 1.0e-15);
  AMPS_EXPECT_NEAR(output[1], -3.5, 1.0e-15);
  AMPS_EXPECT_NEAR(output[2], -30.5, 1.0e-15);
}

AMPS_TEST_CASE("INTERFACE-GEOPACK-CONVENTION-001",
               "production user-to-GSM position and GSM-to-user field flow") {
  /* +90 degree active rotation about +Z, and its inverse. */
  const double user_to_gsm[3][3] = {
      {0.0, -1.0, 0.0},
      {1.0, 0.0, 0.0},
      {0.0, 0.0, 1.0}};
  const double gsm_to_user[3][3] = {
      {0.0, 1.0, 0.0},
      {-1.0, 0.0, 0.0},
      {0.0, 0.0, 1.0}};
  std::memcpy(Geopack::UserFrame2GSM, user_to_gsm, sizeof(user_to_gsm));
  std::memcpy(Geopack::GSM2UserFrame, gsm_to_user, sizeof(gsm_to_user));
  Geopack::Rotate2GSM = true;

  const double earth_radius_m = 6371000.0;
  double position_m[3] = {
      2.0 * earth_radius_m, 3.0 * earth_radius_m, 5.0 * earth_radius_m};
  double field_t[3] = {0.0, 0.0, 0.0};
  Geopack::IGRF::GetMagneticField(field_t, position_m);

  /* user (2,3,5) -> GSM (-3,2,5), all in Earth radii. */
  AMPS_EXPECT_NEAR(intercepted_position[0], -3.0, 1.0e-14);
  AMPS_EXPECT_NEAR(intercepted_position[1], 2.0, 1.0e-14);
  AMPS_EXPECT_NEAR(intercepted_position[2], 5.0, 1.0e-14);
  /* fixed GSM (7,11,13) nT -> user (11,-7,13) nT -> tesla. */
  AMPS_EXPECT_NEAR(field_t[0], 11.0e-9, 1.0e-22);
  AMPS_EXPECT_NEAR(field_t[1], -7.0e-9, 1.0e-22);
  AMPS_EXPECT_NEAR(field_t[2], 13.0e-9, 1.0e-22);
}

AMPS_TEST_CASE("INTERFACE-GEOPACK-NOSPICE-001",
               "production no-SPICE fallback is identity with explicit false flag") {
  double forward[3][3];
  double reverse[3][3];
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j) {
      forward[i][j] = 100.0 + 10.0 * i + j;
      reverse[i][j] = -100.0 - 10.0 * i - j;
    }
  bool rotate = true;
  Geopack::SetFrameRotation("2012-03-04T05:06:07", "J2000", "GSM",
                            forward, reverse, rotate);
  AMPS_EXPECT_FALSE(rotate);
  ExpectIdentity(forward);
  ExpectIdentity(reverse);

  const double observed =
      Geopack::DriverTimeTagFromYearDoy(2012, 64, 5, 6, 7);
  const double expected =
      ((((2012.0 * 367.0 + 64.0) * 24.0 + 5.0) * 60.0 + 6.0) * 60.0 + 7.0);
  AMPS_EXPECT_NEAR(observed, expected, 0.0);
}
