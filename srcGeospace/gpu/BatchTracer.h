#ifndef _SRC_GEOSPACE_GPU_BATCH_TRACER_H_
#define _SRC_GEOSPACE_GPU_BATCH_TRACER_H_

//======================================================================================
// BatchTracer.h -- front end of the batched trajectory engine used by the application
// drivers.
//
// Trace() is a drop-in replacement for one call of the legacy per-trajectory tracer
// (gridless TraceTrajectoryImpl / Mode3D TraceTrajectory3D).  When it is called from a
// fiber (see Fiber.h) the request is queued and the fiber is suspended until a pool
// kernel finishes the trajectory; otherwise the trajectory is traced synchronously on
// the calling thread with the same engine code.
//
// Routing:
//   NONE / DIPOLE / Mode3D mesh   -> device pool (OpenACC GPU or OpenMP host build)
//   Geopack / Tsyganenko models   -> host pool (OpenMP threads on the CPU)
//
// Environment (all optional; see srcGeospace/doc/BUILD.md for details):
//   AMPS_EARTH_ENGINE=BATCH|LEGACY      BATCH (default) uses this engine; LEGACY keeps the
//                                       original per-trajectory code paths.
//   AMPS_EARTH_ENGINE_VERIFY=1          also run the legacy tracer for every trajectory
//                                       and compare bit for bit (EngineGlue.cpp)
//   AMPS_EARTH_ENGINE_VERBOSE=1         print pool/fiber statistics after each batch
//   AMPS_EARTH_POOL_CAPACITY=N          device-pool slots (CPU 64*threads, GPU 65536)
//   AMPS_EARTH_STEPS_PER_LAUNCH=K       loop iterations per slot and launch (CPU 256,
//                                       GPU 2048)
//   AMPS_EARTH_HOST_POOL_CAPACITY=N     host-pool slots (Geopack/Tsyganenko models)
//   AMPS_EARTH_HOST_STEPS_PER_LAUNCH=K  host-pool iterations per launch (64)
//   AMPS_EARTH_ASYNC=0|1                overlap device-pool launches with fiber work
//                                       (default 1)
//   AMPS_EARTH_HOST_RUNTIME=OPENMP      CPU kernels via OpenMP instead of the thread pool
//   AMPS_EARTH_SPIN_US=N                thread-pool spin before blocking (50 us)
//   AMPS_EARTH_FIBERS=N                 maximum live fibers (4096; GPU: pool capacity)
//   AMPS_EARTH_FIBER_STACK_KB=N         fiber stack size (1024 KB, reserved lazily)
//   AMPS_EARTH_FIBER_GUARD=1            guard pages for every fiber stack
//   AMPS_EARTH_TURN_CAPACITY=N          drift-turn diagnostic states per trajectory (512)
//======================================================================================

#include "EngineHost.h"
#include "EngineDeviceAPI.h"
#include "Fiber.h"

#include <cstdint>
#include <deque>
#include <mutex>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace Geospace {
namespace GPU {

// Thrown by Trace() in a fiber whose request was cancelled (Speculation.h).  Only the
// speculative helpers cancel requests, and they catch it.
struct TraceCancelled : std::exception {
  const char* what() const noexcept override { return "trajectory request cancelled"; }
};

struct TraceOutput {
  EgResult res;
  std::vector<EgTurnState> turns;
  std::string error;           // host field exception message (res.error==kErrHostField)
};

class BatchTracer {
 public:
  static BatchTracer& Get();

  // AMPS_EARTH_ENGINE != LEGACY
  static bool Enabled();

  // Register (or find) an immutable context; returns its index.  Thread-safe.
  int RegisterContext(const EgContext& c);

  // Flattened Mode3D mesh used by kFieldMesh requests.  Must be called while no
  // trajectory is in flight; `generation` identifies the published snapshot.  In CPU
  // double/SI builds the engine aliases the arrays of `m` (no copy), so they must stay
  // valid and unchanged until the next SetMesh()/ClearMesh(); float and GPU builds copy.
  void SetMesh(const EgMeshHost& m,unsigned long long generation);
  void ClearMesh();
  bool HaveMesh(unsigned long long generation) const;

  // Trace one request (blocking).  Turn states are returned only when the context
  // records them.
  void Trace(const EgRequest& rq,const HostFieldRef& host,TraceOutput& out);

  // Cancel the trajectory request of a fiber (speculative evaluation).  A request that
  // is still queued is dropped, one in the device pool is stopped at the next launch,
  // and the fiber resumes at once with TraceCancelled; a late result is discarded.  If
  // the fiber is not waiting for a trajectory, its next Trace() throws TraceCancelled
  // until ClearCancel(fiber).
  void Cancel(Fiber::FiberRec* fiber);
  void ClearCancel(Fiber::FiberRec* fiber);

  // Kernels execute on an accelerator.
  bool Accelerator() const { return accelerator_; }

  // Concurrency hints for the drivers.
  int PreferredFibers() const;
  int HostThreads() const { return hostThreads_; }
  void SetHostThreads(int n);

  // Fiber idle hook: advance the pools and wake fibers.  Returns false when there is
  // nothing in flight.
  bool Progress();

  std::string BackendName() const;
  void PrintStatistics(const char* label);

 private:
  BatchTracer();
  void EnsurePools();
  void Deliver(const EgResult& res,std::vector<EgTurnState>* turns,const std::string* err);

  struct Pending {
    EgRequest rq;
    HostFieldRef host;
    TraceOutput* out;
    Fiber::FiberRec* fiber;
    bool device;
    bool submitted=false;     // handed to a pool
    bool cancelled=false;
    bool resultArrived=false; // result delivered (or nothing will arrive)
    bool fiberGone=false;     // the fiber has left Trace()
  };

  mutable std::mutex mutex_;
  std::vector<EgContext> contexts_;
  int devUploaded_=0;
  int hostUploaded_=0;

  EgDevicePool* dev_=nullptr;
  std::unique_ptr<HostPool> host_;
  HostMesh hostMesh_;
  unsigned long long meshGeneration_=0;
  bool haveMesh_=false;

  std::deque<Pending*> queue_;
  std::unordered_map<long long,Pending*> inflight_;
  std::unordered_map<Fiber::FiberRec*,Pending*> byFiber_;
  std::unordered_set<Fiber::FiberRec*> cancelledFibers_;
  long long statCancelled_=0;
  long long nextTag_=1;

  int capacity_=0;
  int steps_=0;
  int turnCapacity_=0;
  int hostThreads_=1;
  bool accelerator_=false;
  bool asyncDevice_=true;     // AMPS_EARTH_ASYNC (default 1)

  // statistics
  long long statTraced_=0,statSync_=0,statLaunches_=0,statDelivered_=0;
  double statKernelSeconds_=0.0;
  double statSubmitSeconds_=0.0,statWaitSeconds_=0.0,statLaunchSeconds_=0.0,
         statCollectSeconds_=0.0,statDeliverSeconds_=0.0;
  long long statActiveSum_=0;
};

} // namespace GPU
} // namespace Geospace

#endif
