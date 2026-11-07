#include <stdint.h>
#include <stddef.h>

/*
 * This fixture is intentionally written as C and compiled by the configured C
 * compiler.  The C++ caller must therefore cross a real language-linkage
 * boundary; declaring the implementation in the C++ test would only test C++
 * against itself.  Every output combines nonadjacent inputs so swapped scalar,
 * vector, matrix, or byte-buffer arguments cannot accidentally pass.
 */
void amps_c_abi_echo(const int32_t *tag, const double *scale,
                     const double vector[3], const double matrix[6],
                     const unsigned char bytes[4], int32_t *tag_out,
                     double vector_out[3], double probes_out[3],
                     uint64_t *byte_code_out) {
  *tag_out = *tag + 37;
  vector_out[0] = vector[2] + *scale;
  vector_out[1] = vector[0] - 2.0 * vector[1];
  vector_out[2] = vector[1] + (double)(*tag);

  /* matrix is a two-row, three-column row-major C array. */
  probes_out[0] = matrix[0] + 10.0 * matrix[5];
  probes_out[1] = matrix[1] - matrix[3];
  probes_out[2] = matrix[2] + matrix[4];

  *byte_code_out = ((uint64_t)bytes[0] << 24) |
                   ((uint64_t)bytes[1] << 16) |
                   ((uint64_t)bytes[2] << 8) |
                   (uint64_t)bytes[3];
}
