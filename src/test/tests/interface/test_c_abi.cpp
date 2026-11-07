#include "test_harness.h"

#include <cstdint>

extern "C" {
void amps_c_abi_echo(const std::int32_t* tag, const double* scale,
                     const double vector[3], const double matrix[6],
                     const unsigned char bytes[4], std::int32_t* tag_out,
                     double vector_out[3], double probes_out[3],
                     std::uint64_t* byte_code_out);
}

AMPS_TEST_CASE("INTERFACE-C-ABI-001",
               "C/C++ scalar, array, byte, and argument-order round trip") {
  const std::int32_t tag = 17;
  const double scale = -2.5;
  const double vector[3] = {1.25, -3.0, 7.5};
  const double matrix[6] = {11.0, 12.0, 13.0, 21.0, 22.0, 23.0};
  const unsigned char bytes[4] = {0x12, 0x34, 0x56, 0x78};
  std::int32_t tag_out = 0;
  double vector_out[3] = {0.0, 0.0, 0.0};
  double probes_out[3] = {0.0, 0.0, 0.0};
  std::uint64_t byte_code = 0;

  amps_c_abi_echo(&tag, &scale, vector, matrix, bytes, &tag_out,
                  vector_out, probes_out, &byte_code);

  AMPS_EXPECT_EQ(tag_out, 54);
  AMPS_EXPECT_NEAR(vector_out[0], 5.0, 1.0e-15);
  AMPS_EXPECT_NEAR(vector_out[1], 7.25, 1.0e-15);
  AMPS_EXPECT_NEAR(vector_out[2], 14.0, 1.0e-15);
  AMPS_EXPECT_NEAR(probes_out[0], 241.0, 1.0e-15);
  AMPS_EXPECT_NEAR(probes_out[1], -9.0, 1.0e-15);
  AMPS_EXPECT_NEAR(probes_out[2], 35.0, 1.0e-15);
  AMPS_EXPECT_EQ(byte_code, UINT64_C(0x12345678));
}
