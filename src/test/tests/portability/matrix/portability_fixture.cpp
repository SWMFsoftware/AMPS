/*
 * Dependency-light T12 fixture for real AMPS array storage.
 *
 * This is not a replacement for configured PIC timing.  It deliberately uses
 * array_1d from src/general so compiler, sanitizer, OpenMP, leak and bounded
 * throughput profiles all exercise the same production ownership/indexing
 * implementation.  The independent checksum is the closed-form sum
 * N(N-1)/2; timing is emitted separately and can never turn a wrong checksum
 * into a performance PASS.  All values are dimensionless fixture values.
 */

#include "array_1d.h"

#include <omp.h>

#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <stdexcept>
#include <string>
#include <utility>

/* The GNU-labelled matrix cases must not silently run through Clang's GNU
 * compatibility macros or an unrelated configured wrapper.  The runner
 * enables this guard only for those cases; the separate Clang case compiles
 * the same production fixture without it. */
#if defined(AMPS_T12_REQUIRE_GNU)
#if !defined(__GNUC__) || defined(__clang__)
#error "T12 GNU profile requires the GNU C++ compiler family, not Clang"
#endif
#endif

int ThisThread = 0;
int TotalThreadsNumber = 1;
unsigned int ExitErrorCode = 0;
MPI_Comm MPI_GLOBAL_COMMUNICATOR = MPI_COMM_WORLD;

long int nint(double value) {
  const long int lower = static_cast<long int>(value);
  return value - lower > 0.5 ? lower + 1 : lower;
}

void PrintErrorLog(const char*) {}
void PrintErrorLog(long int, const char*, const char*) {}

/* Production arrays route invalid extents through the AMPS exit hook.  This
 * fixture does not provoke that path, but a throwing hook makes any accidental
 * call a contained nonzero result instead of terminating the test runner. */
void exit(long int, const char*, const char* message) {
  throw std::runtime_error(message == NULL ? "AMPS trap exit" : message);
}

namespace {

double ExpectedSum(long int count) {
  return 0.5 * static_cast<double>(count) *
      static_cast<double>(count - 1);
}

bool FillAndCheck(long int count, int repetitions, double& elapsed_seconds) {
  const auto start = std::chrono::steady_clock::now();
  array_1d<double> values(static_cast<int>(count));
  double observed = 0.0;
  for (int repetition = 0; repetition < repetitions; ++repetition) {
    observed = 0.0;
    for (long int index = 0; index < count; ++index)
      values(static_cast<int>(index)) = static_cast<double>(index);
    for (long int index = 0; index < count; ++index)
      observed += values(static_cast<int>(index));
  }
  const auto stop = std::chrono::steady_clock::now();
  elapsed_seconds = std::chrono::duration<double>(stop - start).count();
  return observed == ExpectedSum(count);
}

int Correctness() {
  array_1d<double> original(1024);
  for (int index = 0; index < original.GetSize(); ++index)
    original(index) = static_cast<double>(index);
  array_1d<double> copy(original);
  copy(0) = -1.0;
  array_1d<double> moved(std::move(copy));
  const bool pass = original(0) == 0.0 && moved(0) == -1.0 &&
      copy.GetSize() == 0 && original.sum() == ExpectedSum(1024);
  std::printf("CORRECTNESS status=%s checksum=%.17g\n",
              pass ? "PASS" : "FAIL", original.sum());
  return pass ? 0 : 1;
}

int OpenMPInvariance(int threads) {
  const int count = 200000;
  array_1d<double> values(count);
  omp_set_dynamic(0);
  omp_set_num_threads(threads);
#pragma omp parallel for schedule(static)
  for (int index = 0; index < count; ++index)
    values(index) = static_cast<double>(index);
  double observed = 0.0;
#pragma omp parallel for reduction(+:observed) schedule(static)
  for (int index = 0; index < count; ++index) observed += values(index);
  const bool pass = observed == ExpectedSum(count);
  std::printf("OPENMP status=%s requested_threads=%d checksum=%.17g\n",
              pass ? "PASS" : "FAIL", threads, observed);
  return pass ? 0 : 1;
}

int ArrayPerformance() {
  double seconds = 0.0;
  const long int count = 1000000;
  const int repetitions = 8;
  const bool pass = FillAndCheck(count, repetitions, seconds);
  std::printf("PERF_ARRAY status=%s count=%ld repetitions=%d seconds=%.9g\n",
              pass ? "PASS" : "FAIL", count, repetitions, seconds);
  return pass ? 0 : 1;
}

int IOPerformance(const char* file_name) {
  const int count = 250000;
  array_1d<double> source(count), restored(count);
  for (int index = 0; index < count; ++index)
    source(index) = static_cast<double>(index);
  const auto start = std::chrono::steady_clock::now();
  {
    std::ofstream output(file_name, std::ios::binary | std::ios::trunc);
    output.write(reinterpret_cast<const char*>(source.get_data_ptr()),
                 static_cast<std::streamsize>(count * sizeof(double)));
    if (!output) return 2;
  }
  {
    std::ifstream input(file_name, std::ios::binary);
    input.read(reinterpret_cast<char*>(restored.get_data_ptr()),
               static_cast<std::streamsize>(count * sizeof(double)));
    if (!input) return 2;
  }
  const auto stop = std::chrono::steady_clock::now();
  double checksum = 0.0;
  for (int index = 0; index < count; ++index) checksum += restored(index);
  const bool pass = checksum == ExpectedSum(count);
  const double seconds = std::chrono::duration<double>(stop - start).count();
  std::printf("PERF_IO status=%s bytes=%zu seconds=%.9g checksum=%.17g\n",
              pass ? "PASS" : "FAIL", static_cast<size_t>(count) * sizeof(double),
              seconds, checksum);
  return pass ? 0 : 1;
}

}  // namespace

int main(int argc, char** argv) {
  try {
    if (argc == 2 && std::string(argv[1]) == "correctness")
      return Correctness();
    if (argc == 3 && std::string(argv[1]) == "openmp")
      return OpenMPInvariance(std::atoi(argv[2]));
    if (argc == 2 && std::string(argv[1]) == "performance-array")
      return ArrayPerformance();
    if (argc == 3 && std::string(argv[1]) == "performance-io")
      return IOPerformance(argv[2]);
    std::fprintf(stderr, "usage: %s correctness|openmp N|performance-array|performance-io FILE\n",
                 argv[0]);
    return 2;
  }
  catch (const std::exception& error) {
    std::fprintf(stderr, "fixture error: %s\n", error.what());
    return 2;
  }
}
