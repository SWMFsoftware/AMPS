#ifndef _SRC_EARTH_3D_MODE3DPARALLEL_H_
#define _SRC_EARTH_3D_MODE3DPARALLEL_H_

#include "../util/amps_param_parser.h"

namespace Earth {
namespace Mode3D {

// Shared intra-rank backend used by mesh-backed backward products in Mode3D.
// The same controls intentionally drive cutoff rigidity and density/flux so a
// single AMPS_PARAM.in deck selects the backend for all expensive backtracking
// products that operate on the already-materialized AMR magnetic-field snapshot.
enum class ParallelBackend {
  OPENMP,   // OpenMP team inside each MPI rank.
  THREADS,  // Direct std::thread worker pool inside each MPI rank.
  SERIAL    // No intra-rank shared-memory parallelism.
};

// Return a stable printable name for logs and diagnostics.
const char* ParallelBackendName(ParallelBackend backend);

// Resolve the backend from prm.mode3d.densityParallelBackend, falling back to
// AMPS_MODE3D_DENSITY_PARALLEL.  The historical DENSITY_* keyword names are kept
// for compatibility, but in Mode3D they now mean the generic backward-product
// backend and therefore apply to both cutoff and density/flux.
ParallelBackend ResolveParallelBackend(const EarthUtil::AmpsParam& prm,
                                       const char* diagnosticContext="Mode3D");

// Resolve the number of intra-rank workers from prm.mode3d.densityThreads,
// AMPS_MODE3D_DENSITY_THREADS, OMP_NUM_THREADS, OpenMP runtime defaults, and
// finally std::thread::hardware_concurrency().  Returns at least 1.
int ResolveParallelThreadCount(const EarthUtil::AmpsParam& prm,
                               ParallelBackend backend);

// For the direct std::thread backend, widen the current MPI rank's CPU-affinity
// mask before worker threads are created.  This is a no-op for OPENMP/SERIAL and
// for one-worker THREADS runs.  The operation is intentionally process-global and
// is applied only once per process because later worker pools inherit the mask.
void ApplyWideAffinityForDirectThreadsOnce(ParallelBackend backend,
                                           int workerCount,
                                           const char* diagnosticContext="Mode3D");

} // namespace Mode3D
} // namespace Earth

#endif // _SRC_EARTH_3D_MODE3DPARALLEL_H_
