// Hybrid/OpenMP restart contract for the production AMPS RNG.
//
// Every OpenMP worker owns its uniform seed and thread_local Box--Muller spare.
// A restart must reproduce both layers independently of the number of normal
// samples drawn before reseeding.  This executable is separate from the
// dependency-light statistical fixture because it compiles the production
// rnd.cpp under the hybrid configuration and initializes MPI/OpenMP normally.

#include "specfunc.h"

#include <mpi.h>
#include <omp.h>

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <vector>

int ThisThread = 0;
int TotalThreadsNumber = 1;
int deviceThisThread = 0;
int deviceTotalThreadsNumber = 1;
int cudaThreadLimitMallocHeapSize = 0;
unsigned int ExitErrorCode = 0;
MPI_Comm MPI_GLOBAL_COMMUNICATOR = MPI_COMM_WORLD;

void exit(long int, const char*, const char*) {
  MPI_Abort(MPI_GLOBAL_COMMUNICATOR, 97);
  std::abort();
}

struct Snapshot {
  std::vector<double> uniform;
  std::vector<double> normal;
  bool oracle_match;
};

Snapshot Capture(int seed, int thread_count) {
  Snapshot result;
  result.uniform.resize(thread_count);
  result.normal.resize(thread_count);
  int oracle_match = 1;

#pragma omp parallel reduction(&:oracle_match)
  {
    const int thread = omp_get_thread_num();
    const unsigned long stream_seed = (thread == 0) ?
        static_cast<unsigned long>(seed) :
        static_cast<unsigned long>(std::abs(seed) + thread);
    cRndSeedContainer oracle;
    oracle.Seed = stream_seed;
    result.uniform[thread] = rnd();
    if (result.uniform[thread] != rnd(&oracle)) oracle_match = 0;
    // One draw deliberately leaves the paired normal cached on every worker.
    result.normal[thread] = Vector3D::Distribution::Normal();
  }
  result.oracle_match = oracle_match != 0;
  return result;
}

int main(int argc, char** argv) {
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_GLOBAL_COMMUNICATOR, &ThisThread);
  MPI_Comm_size(MPI_GLOBAL_COMMUNICATOR, &TotalThreadsNumber);
  omp_set_dynamic(0);
  omp_set_num_threads(4);

  const int seed = 17;
  rnd_seed(seed);
  const int thread_count = omp_get_max_threads();
  Snapshot first = Capture(seed, thread_count);
  rnd_seed(seed);
  Snapshot replay = Capture(seed, thread_count);

  bool pass = TotalThreadsNumber == 1 && thread_count == 4 &&
              first.oracle_match && replay.oracle_match;
  for (int thread = 0; thread < thread_count; ++thread) {
    pass = pass && first.uniform[thread] == replay.uniform[thread] &&
           first.normal[thread] == replay.normal[thread] &&
           std::isfinite(first.normal[thread]);
    for (int other = 0; other < thread; ++other)
      pass = pass && first.uniform[thread] != first.uniform[other];
  }

  if (ThisThread == 0) {
    std::printf("[GENERAL-RNG-STREAMS-001] %s "
                "four-worker hybrid uniform/normal reseed contract\n",
                pass ? "PASS" : "FAIL");
    std::fflush(stdout);
  }
  MPI_Finalize();
  return pass ? 0 : 1;
}
