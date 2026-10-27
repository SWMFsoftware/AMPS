#ifndef _SRC_GEOSPACE_GPU_DRIVER_SUPPORT_H_
#define _SRC_GEOSPACE_GPU_DRIVER_SUPPORT_H_

//======================================================================================
// DriverSupport.h -- helpers used by the application drivers to run their existing worker
// loops as fibers that feed the batched trajectory engine.
//
// The drivers keep their task decomposition, MPI scheduling, result accumulation and
// output.  In batch mode, the per-rank worker pool (std::thread / OpenMP / serial)
// is replaced by `nWorkers` fibers on the rank's main thread; the requested worker
// thread count becomes the number of OpenMP threads of the engine's host executor (or
// is irrelevant on a GPU).
//======================================================================================

#include <functional>

namespace Geospace {
namespace GPU {

// True when drivers should use fibers + the batched engine (AMPS_EARTH_ENGINE!=LEGACY).
bool DriverBatchEnabled();

// Number of fibers for a pool that has `nTasks` independent tasks (0 = unknown) when
// the engine host executor uses `hostThreads` threads.
int DriverFiberCount(long long nTasks,int hostThreads);

// Run worker(w), w=0..nWorkers-1, as fibers.  service() is called on the rank's main
// thread between trajectory-pool advances (progress publication through MPI is safe).
// hostThreads sets the engine's host-executor thread count.  Rethrows the first
// exception raised by a worker after all workers have stopped.
void RunDriverFibers(int nWorkers,int hostThreads,
                     const std::function<void(int worker)>& worker,
                     const std::function<void()>& service,
                     const char* label);

// Parallel loop for code that already runs inside a driver fiber (nested work such as
// per-energy or per-direction loops).  Outside fibers it falls back to a plain loop.
void FiberLoop(long long begin,long long end,
               const std::function<void(long long index,int worker)>& body);

} // namespace GPU
} // namespace Geospace

#endif
