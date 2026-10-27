#ifndef _SRC_GEOSPACE_GPU_FIBER_H_
#define _SRC_GEOSPACE_GPU_FIBER_H_

//======================================================================================
// Fiber.h -- cooperative fibers that turn the sequential legacy search algorithms
// into producers of batched trajectory requests.
//
// WHY FIBERS
// ----------
// The physics drivers (UPPER_SCAN/BINARY/PENUMBRA_SCAN cutoff searches, the
// adaptive DIRECT_ACCESS refinement, Mode3D isotropic bracket shrinking, density and
// spectrum direction loops) are written as ordinary sequential code that calls a
// blocking "trace one trajectory" function many times.  A GPU needs tens of thousands
// of trajectories in flight.  Rewriting every search as an explicit state machine would
// duplicate several thousand lines of validated logic.  Instead each legacy worker
// runs as a fiber: when it requests a trajectory, the request is queued and the fiber is
// suspended; once every fiber is suspended the scheduler advances the whole trajectory
// pool (one OpenACC/OpenMP kernel over all active trajectories) and resumes the fibers
// whose trajectories finished.  The search code, retry policy and output writers are
// unchanged, and each search sees exactly the trajectory results it would see serially.
//
// RULES
// -----
//  * Fibers run on the calling OS thread (the MPI rank's main thread).  MPI calls made
//    from fiber code are therefore made by the main thread.
//  * InFiber() is false inside an OpenMP parallel region; trajectory requests made there
//    fall back to synchronous tracing on the calling thread.
//  * An exception that leaves a fiber body is captured and rethrown by ParallelFor().
//======================================================================================

#include <cstddef>
#include <exception>
#include <functional>

namespace Geospace {
namespace GPU {
namespace Fiber {

struct FiberRec;

// True when the caller is a fiber body that may block (and not inside OpenMP).
bool InFiber();

// Suspend the current fiber until Wake() is called for it.  Only valid when InFiber().
void Block();

// Handle of the running fiber (nullptr outside fibers).
FiberRec* Current();

// Make a blocked fiber runnable again.  Called by the batch tracer.
void Wake(FiberRec* fiber);

// Idle hook: called by the scheduler when no fiber is runnable.  It must make progress
// (advance trajectories and Wake() fibers) and return true, or return false when it has
// nothing to do (which with blocked fibers is a deadlock and raises an error).
void SetIdleHook(bool (*hook)());

// Optional periodic hook executed from the scheduler context (not from a fiber) between
// batch advances; used to publish progress through MPI from the rank's main thread.
void SetServiceHook(std::function<void()> hook);

// Run body(index,worker) for every index in [begin,end) on at most maxWorkers fibers.
// Called outside fibers, it creates the scheduler, runs all fibers to completion and
// returns.  Called from a fiber, it spawns child fibers into the running scheduler and
// blocks the caller until they finish, so nested loops remain fully batched.
void ParallelFor(long long begin,long long end,int maxWorkers,
                 const std::function<void(long long index,int worker)>& body);

// Maximum number of simultaneously live fibers: AMPS_EARTH_FIBERS if set, otherwise the
// default installed by SetDefaultMaxFibers() (4096 unless changed; the batch tracer
// raises it to the device-pool capacity on accelerators).
int MaxFibers();
void SetDefaultMaxFibers(int n);

// Statistics of the most recent top-level ParallelFor.
struct Stats {
  long long fibersCreated;
  long long switches;
  long long idleCalls;
  double fiberSeconds=0.0;     // time inside fibers (driver code)
  double idleSeconds=0.0;      // time in the idle hook (engine progress)
  double serviceSeconds=0.0;   // time in the service hook
};
Stats LastStats();

} // namespace Fiber
} // namespace GPU
} // namespace Geospace

#endif
