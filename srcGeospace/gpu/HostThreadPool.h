#ifndef _SRC_GEOSPACE_GPU_HOST_THREAD_POOL_H_
#define _SRC_GEOSPACE_GPU_HOST_THREAD_POOL_H_

//======================================================================================
// HostThreadPool.h -- persistent worker threads for the CPU kernel launches.
//
// WHY NOT ONLY OPENMP
// -------------------
// The batched engine alternates between a short parallel kernel (advance every active
// trajectory by K iterations) and a serial host phase on the rank's main thread
// (resume fibers, run search logic, submit new requests).  With libgomp's default
// OMP_WAIT_POLICY, idle OpenMP workers spin between parallel regions.  Test runners
// run e.g. 4 MPI ranks x 16 threads on 36 hardware threads; the spinning workers of all
// ranks then starve every rank's serial phase.  Measured on F1: 28 s with the default
// policy, 3.1 s with OMP_WAIT_POLICY=PASSIVE.  The policy can only be set in the
// environment before the program starts (libgomp reads it when it is loaded), and the
// test runners do not set it.
//
// This pool provides the same "parallel for, dynamic schedule" with workers that spin
// only briefly (AMPS_EARTH_SPIN_US, default 50 us) and then block on a condition
// variable.  It is the default CPU launcher; AMPS_EARTH_HOST_RUNTIME=OPENMP selects the
// OpenMP pragma launch instead.  The kernel body is the same code in every case.
//======================================================================================

#include <functional>

namespace Geospace {
namespace GPU {

// Run body(i) for i in [0,n) on `nThreads` threads (the caller participates), with
// dynamic scheduling of single iterations.  Returns after every iteration finished.
// Not reentrant: one launch at a time per process (the engine launches from the rank's
// main thread only).
void HostParallelFor(int n,int nThreads,const std::function<void(int)>& body);

// Asynchronous variant used by the device pool on CPU builds: body(i) runs on
// `nThreads` worker threads while the caller continues (e.g. resuming fibers); the
// body is copied.  HostParallelWait() blocks until the launch finished.  The async
// launches use a pool instance separate from HostParallelFor(), so a synchronous
// host-pool launch may run while a device-pool launch is in flight.
void HostParallelForAsync(int n,int nThreads,std::function<void(int)> body);
void HostParallelWait();

// True when the CPU kernels should use OpenMP instead of HostParallelFor.
bool HostRuntimeIsOpenMP();

} // namespace GPU
} // namespace Geospace

#endif
