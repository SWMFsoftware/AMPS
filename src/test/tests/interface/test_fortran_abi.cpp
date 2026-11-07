#include "test_harness.h"

#include <cstdint>

extern "C" {
/* Historical AMPS wrappers declare legacy Fortran procedures this way. */
void amps_fortran_legacy_echo_(const std::int32_t* tag, const double* scalar,
                               const double vector[3], const double matrix[9],
                               std::int32_t* tag_out, double* scalar_out,
                               double vector_out[3], double probes_out[3]);

/* Standards-defined control supplied by ISO_C_BINDING/BIND(C). */
void amps_fortran_bindc_echo(const std::int32_t* tag, const double* scalar,
                             const double vector[3], const double matrix[9],
                             std::int32_t* tag_out, double* scalar_out,
                             double vector_out[3], double probes_out[3]);
}

namespace {

const std::int32_t kTag = 17;
const double kScalar = -2.5;
const double kVector[3] = {1.25, -3.0, 7.5};

/*
 * C++ storage below is deliberately the Fortran column-major sequence for
 * matrix(3,3): columns are (11,21,31), (12,22,32), (13,23,33).  Non-symmetric
 * values make a row/column transposition observable.
 */
const double kFortranMatrix[9] = {
    11.0, 21.0, 31.0,
    12.0, 22.0, 32.0,
    13.0, 23.0, 33.0};

}  // namespace

AMPS_TEST_CASE("INTERFACE-FORTRAN-LEGACY-001",
               "legacy underscore Fortran ABI and column-major argument order") {
  std::int32_t tag_out = 0;
  double scalar_out = 0.0;
  double vector_out[3] = {0.0, 0.0, 0.0};
  double probes_out[3] = {0.0, 0.0, 0.0};
  amps_fortran_legacy_echo_(&kTag, &kScalar, kVector, kFortranMatrix,
                            &tag_out, &scalar_out, vector_out, probes_out);

  AMPS_EXPECT_EQ(tag_out, 118);
  AMPS_EXPECT_NEAR(scalar_out, -27.0, 1.0e-15);
  AMPS_EXPECT_NEAR(vector_out[0], 24.5, 1.0e-15);
  AMPS_EXPECT_NEAR(vector_out[1], 7.25, 1.0e-15);
  AMPS_EXPECT_NEAR(vector_out[2], -5.5, 1.0e-15);
  AMPS_EXPECT_NEAR(probes_out[0], 3211.0, 1.0e-15);
  AMPS_EXPECT_NEAR(probes_out[1], 11.0, 1.0e-15);
  AMPS_EXPECT_NEAR(probes_out[2], 53.0, 1.0e-15);
}

AMPS_TEST_CASE("INTERFACE-FORTRAN-ISO-C-001",
               "ISO_C_BINDING Fortran ABI and by-reference array round trip") {
  std::int32_t tag_out = 0;
  double scalar_out = 0.0;
  double vector_out[3] = {0.0, 0.0, 0.0};
  double probes_out[3] = {0.0, 0.0, 0.0};
  amps_fortran_bindc_echo(&kTag, &kScalar, kVector, kFortranMatrix,
                          &tag_out, &scalar_out, vector_out, probes_out);

  AMPS_EXPECT_EQ(tag_out, 8);
  AMPS_EXPECT_NEAR(scalar_out, 6.0, 1.0e-15);
  AMPS_EXPECT_NEAR(vector_out[0], -0.5, 1.0e-15);
  AMPS_EXPECT_NEAR(vector_out[1], 8.75, 1.0e-15);
  AMPS_EXPECT_NEAR(vector_out[2], 5.75, 1.0e-15);
  AMPS_EXPECT_NEAR(probes_out[0], 22.0, 1.0e-15);
  AMPS_EXPECT_NEAR(probes_out[1], 151.0, 1.0e-15);
  AMPS_EXPECT_NEAR(probes_out[2], 10.0, 1.0e-15);
}
