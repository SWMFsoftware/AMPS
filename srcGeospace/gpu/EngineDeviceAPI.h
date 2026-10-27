#ifndef _SRC_GEOSPACE_GPU_ENGINE_DEVICE_API_H_
#define _SRC_GEOSPACE_GPU_ENGINE_DEVICE_API_H_

//======================================================================================
// EngineDeviceAPI.h -- C interface of the device trajectory pool.
//
// The device pool advances trajectories whose magnetic field can be evaluated in a
// kernel (NONE, analytic DIPOLE, Mode3D compact mesh).  It is implemented in
// EngineDevice.cpp, which is compiled
//   * by the AMPS C++ compiler with OpenMP for CPU builds (linked into libAMPS), or
//   * by nvc++ -acc=gpu into libGeospaceEngine.so for GPU builds.
// Only plain C/POD types cross this interface so the two compilers never exchange C++
// library objects.
//======================================================================================

#include "EngineTypes.h"

namespace Geospace { namespace GPU { struct EgMeshHost; } }

extern "C" {

struct EgDevicePool;

// Returns a short description of the compiled backend, e.g. "OpenACC(gpu)/double".
const char* eg_dev_backend_name();
// 1 when kernels execute on an accelerator.
int eg_dev_is_accelerator();

EgDevicePool* eg_dev_create(int capacity,int stepsPerLaunch,int turnCapacity);
void eg_dev_destroy(EgDevicePool* pool);

// Replace the context table (only when the pool has no active trajectory).
int eg_dev_set_contexts(EgDevicePool* pool,const Geospace::GPU::EgContext* ctx,int n);

// Upload a flattened Mode3D mesh (only when idle).  nullptr clears it.
int eg_dev_set_mesh(EgDevicePool* pool,const Geospace::GPU::EgMeshHost* mesh);

int eg_dev_capacity(EgDevicePool* pool);
int eg_dev_active(EgDevicePool* pool);
int eg_dev_free(EgDevicePool* pool);

// Queue up to n requests into free slots; returns the number accepted.
int eg_dev_submit(EgDevicePool* pool,const Geospace::GPU::EgRequest* reqs,int n);

// One synchronous kernel launch over all active slots.  Returns the number of finished
// trajectories that can be collected.
int eg_dev_advance(EgDevicePool* pool);

// Asynchronous launch over the currently active slots (returns their number, 0 when
// none, -1 when a launch is already pending).  Requests submitted while the launch is
// running are integrated by the next launch.  eg_dev_wait() completes the launch and
// returns the number of collectable results.  This overlaps the host's fiber work
// (search logic, new requests) with trajectory integration on the device/threads.
int eg_dev_launch(EgDevicePool* pool);
int eg_dev_wait(EgDevicePool* pool);
int eg_dev_pending_launch(EgDevicePool* pool);

// Copy finished results (at most maxOut) and free their slots.  When turns!=nullptr,
// turnStates of trajectory k are written to turns[k*turnStride ...] and their number to
// res[k].turnStateCount.  Returns the number of results written.
int eg_dev_collect(EgDevicePool* pool,Geospace::GPU::EgResult* res,int maxOut,
                   Geospace::GPU::EgTurnState* turns,int turnStride);

// Advance steps per launch (K).
void eg_dev_set_steps_per_launch(EgDevicePool* pool,int K);

// Cancel the request with this tag (speculative evaluation).  A request that is still
// integrating is finished at the next launch with termination kTermCancelled and is
// collected like any other result.  Unknown or already finished tags are ignored.
void eg_dev_cancel(EgDevicePool* pool,long long tag);

// Host threads used by the OpenMP-host backend (ignored by OpenACC GPU builds).
void eg_dev_set_threads(EgDevicePool* pool,int nThreads);

} // extern "C"

#endif
