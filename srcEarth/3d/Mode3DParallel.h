#ifndef _SRC_EARTH_3D_MODE3DPARALLEL_H_
#define _SRC_EARTH_3D_MODE3DPARALLEL_H_

#include "../util/amps_param_parser.h"

#include <mpi.h>

#include <string>

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


// Inter-rank work scheduler used by standalone Mode3D and gridless backward products.
//
// BLOCK_CYCLIC preserves the previous deterministic rank assignment:
//   rank r computes r, r+nRanks, r+2*nRanks, ...
//
// STATIC assigns one contiguous block to each rank and is kept mostly as a
// regression/debug option.  It is usually the least balanced scheduler for shell maps.
//
// DYNAMIC uses an MPI one-sided atomic counter.  Ranks repeatedly fetch the next
// chunk of global locations and process that chunk locally.  This is the recommended
// scheduler for cutoff/density runs where trajectory cost varies strongly with
// latitude, rigidity, or connectivity.
enum class MpiScheduler {
  BLOCK_CYCLIC,
  STATIC,
  DYNAMIC
};

// Return a stable printable name for logs and diagnostics.
const char* MpiSchedulerName(MpiScheduler scheduler);

// Resolve the inter-rank scheduler from input/environment.
// The input keywords are MODE3D_MPI_SCHEDULER / GRIDLESS_MPI_SCHEDULER / BACKTRACK_MPI_SCHEDULER.
// The environment fallback is AMPS_MODE3D_MPI_SCHEDULER.
MpiScheduler ResolveMpiScheduler(const EarthUtil::AmpsParam& prm,
                                 const char* diagnosticContext="Mode3D");

// Resolve the MPI dynamic chunk size.  A value <=0 means automatic.
// The automatic value is deliberately proportional to the per-rank worker count:
// large enough to amortize MPI_Fetch_and_op overhead, but small enough for good
// load balance when expensive trajectories are clustered geographically.
long long ResolveMpiDynamicChunk(const EarthUtil::AmpsParam& prm,
                                 int workerCount,
                                 long long nGlobalLocations);

// Small RAII wrapper around an MPI one-sided atomic global-location counter.
//
// Design contract:
//   * Constructed collectively on all ranks of `comm`.
//   * Rank 0 owns the exposed counter memory; all other ranks expose zero bytes.
//   * FetchNextChunkStart() is called independently by each rank, not by worker
//     threads.  This avoids requiring MPI_THREAD_MULTIPLE.
//   * The returned value is the first global location in the next chunk.  The caller
//     should stop when start >= nGlobalLocations and otherwise compute
//     [start, min(start+chunkSize, nGlobalLocations)).
class DynamicMpiLocationScheduler {
public:
  DynamicMpiLocationScheduler(MPI_Comm comm,
                              long long nGlobalLocations,
                              long long chunkSize,
                              const char* diagnosticContext="Mode3D");
  ~DynamicMpiLocationScheduler();

  DynamicMpiLocationScheduler(const DynamicMpiLocationScheduler&) = delete;
  DynamicMpiLocationScheduler& operator=(const DynamicMpiLocationScheduler&) = delete;

  long long FetchNextChunkStart();

  long long ChunkSize() const { return chunkSize_; }
  long long GlobalLocationCount() const { return nGlobalLocations_; }

private:
  MPI_Comm comm_{MPI_COMM_NULL};
  MPI_Win win_{MPI_WIN_NULL};
  int rank_{0};
  long long nGlobalLocations_{0};
  long long chunkSize_{1};
  long long exposedCounter_{0};
};

// Small MPI RMA completion counter used for live progress reporting in dynamic
// or deterministic collective schedulers.
//
// This class is intentionally separate from DynamicMpiLocationScheduler:
//   * DynamicMpiLocationScheduler counts work that has been ASSIGNED.
//   * DynamicMpiProgressCounter counts work that has actually COMPLETED.
//
// Keeping these counters separate prevents the progress bar from racing ahead when
// a rank fetches a large chunk but has not finished computing it yet.  All ranks
// can call Add(delta) from the rank/main thread after a chunk of tasks finishes;
// rank 0 can periodically call Get() to print an accurate global completion count.
// No worker thread should call this object, so MPI_THREAD_MULTIPLE is not needed.
class DynamicMpiProgressCounter {
public:
  DynamicMpiProgressCounter(MPI_Comm comm,
                            long long totalWork,
                            const char* diagnosticContext="Mode3D progress");
  ~DynamicMpiProgressCounter();

  DynamicMpiProgressCounter(const DynamicMpiProgressCounter&) = delete;
  DynamicMpiProgressCounter& operator=(const DynamicMpiProgressCounter&) = delete;

  // Add completed work units to the global counter.  A non-positive delta is ignored.
  void Add(long long delta);

  // Atomically read the current global completion count without modifying it.
  long long Get() const;

  long long TotalWork() const { return totalWork_; }

private:
  MPI_Comm comm_{MPI_COMM_NULL};
  MPI_Win win_{MPI_WIN_NULL};
  int rank_{0};
  long long totalWork_{0};
  mutable long long exposedCounter_{0};
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
