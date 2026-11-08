/*
 * CUDA compile/runtime boundary for one production AMPS constants header.
 *
 * The kernel checks a simple identity involving SpeedOfLight from
 * src/general/constants.h.  This establishes that the maintained constants
 * header is device-compilable.  It does not claim CUDA PIC support: device
 * availability and kernel execution are a separate case, and absence of a
 * driver/device is reported as SKIP by the Python adapter.
 */

#include "constants.h"

#include <cuda_runtime.h>
#include <cstdio>

__global__ void CheckConstant(int* result) {
  const double light_seconds = SpeedOfLight * 2.0 / SpeedOfLight;
  *result = (light_seconds == 2.0) ? 1 : 0;
}

int main() {
  int device_count = 0;
  const cudaError_t count_status = cudaGetDeviceCount(&device_count);
  if (count_status != cudaSuccess || device_count == 0) {
    std::printf("CUDA_NO_DEVICE status=%s count=%d\n",
                cudaGetErrorString(count_status), device_count);
    return 77;
  }
  int* device_result = NULL;
  if (cudaMalloc(&device_result, sizeof(int)) != cudaSuccess) return 2;
  CheckConstant<<<1, 1>>>(device_result);
  int host_result = 0;
  const cudaError_t copied = cudaMemcpy(&host_result, device_result,
                                        sizeof(int), cudaMemcpyDeviceToHost);
  cudaFree(device_result);
  const bool pass = copied == cudaSuccess && host_result == 1;
  std::printf("CUDA_RUNTIME status=%s devices=%d\n",
              pass ? "PASS" : "FAIL", device_count);
  return pass ? 0 : 1;
}
