//======================================================================================
// Mode3DParallel.cpp
//======================================================================================
//
// Shared parallel-backend selection for mesh-backed Mode3D backward products.
//
// Both cutoff rigidity and density/flux are embarrassingly parallel over observation
// locations once the magnetic field has been materialized on the AMR mesh.  The two
// solvers used to carry separate copies of the same backend-resolution code.  Keeping
// the logic here guarantees that standalone Mode3D and SWMF-coupled Mode3D interpret
// DENSITY_PARALLEL / DENSITY_THREADS exactly the same way for both products.
//======================================================================================

#include "Mode3DParallel.h"

#include "pic.h"

#ifdef _OPENMP
#include <omp.h>
#endif

#include <cstdlib>
#include <sstream>
#include <string>
#include <thread>

namespace Earth {
namespace Mode3D {
namespace {

static int ParsePositiveIntEnv_(const char* name) {
  const char* v = std::getenv(name);
  if (v == nullptr || *v == '\0') return 0;

  try {
    const int n = std::stoi(std::string(v));
    return (n > 0) ? n : 0;
  }
  catch (...) {
    return 0;
  }
}

} // namespace

const char* ParallelBackendName(ParallelBackend backend) {
  switch (backend) {
  case ParallelBackend::OPENMP:  return "OPENMP";
  case ParallelBackend::THREADS: return "THREADS";
  case ParallelBackend::SERIAL:  return "SERIAL";
  }
  return "UNKNOWN";
}

ParallelBackend ResolveParallelBackend(const EarthUtil::AmpsParam& prm,
                                       const char* diagnosticContext) {
  std::string token = prm.mode3d.densityParallelBackend;
  if (token.empty()) {
    const char* env = std::getenv("AMPS_MODE3D_DENSITY_PARALLEL");
    if (env != nullptr) token = env;
  }

  token = EarthUtil::ToUpper(token);
  if (token.empty() || token == "OPENMP" || token == "OMP") {
#ifdef _OPENMP
    return ParallelBackend::OPENMP;
#else
    return ParallelBackend::SERIAL;
#endif
  }

  if (token == "THREAD" || token == "THREADS" || token == "STD_THREAD" ||
      token == "STD_THREADS" || token == "PTHREAD" || token == "PTHREADS") {
    return ParallelBackend::THREADS;
  }

  if (token == "SERIAL" || token == "NONE" || token == "OFF") {
    return ParallelBackend::SERIAL;
  }

  std::ostringstream msg;
  msg << diagnosticContext
      << ": unknown DENSITY_PARALLEL/MODE3D_DENSITY_PARALLEL backend '"
      << token << "'. Valid values: OPENMP, THREADS, SERIAL.";
  exit(__LINE__,__FILE__,msg.str().c_str());
  return ParallelBackend::SERIAL;
}

int ResolveParallelThreadCount(const EarthUtil::AmpsParam& prm,
                               ParallelBackend backend) {
  int n = prm.mode3d.densityThreads;
  if (n <= 0) n = ParsePositiveIntEnv_("AMPS_MODE3D_DENSITY_THREADS");
  if (n <= 0) n = ParsePositiveIntEnv_("OMP_NUM_THREADS");

#ifdef _OPENMP
  if (n <= 0 && backend == ParallelBackend::OPENMP) n = omp_get_max_threads();
#endif

  if (n <= 0) {
    const unsigned hw = std::thread::hardware_concurrency();
    n = (hw > 0) ? static_cast<int>(hw) : 1;
  }

  return (n > 0) ? n : 1;
}

void ApplyWideAffinityForDirectThreadsOnce(ParallelBackend backend,
                                           int workerCount,
                                           const char* diagnosticContext) {
  (void)diagnosticContext;

  if (backend != ParallelBackend::THREADS || workerCount <= 1) return;

  static bool applied = false;
  if (applied) return;
  applied = true;

  PIC::Parallel::SetWideAffinityForScheduler();
}

} // namespace Mode3D
} // namespace Earth
