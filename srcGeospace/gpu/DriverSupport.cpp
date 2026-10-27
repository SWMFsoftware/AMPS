//======================================================================================
// DriverSupport.cpp -- see DriverSupport.h.
//======================================================================================

#include "DriverSupport.h"
#include "BatchTracer.h"
#include "Fiber.h"

#include <algorithm>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <string>

namespace Geospace {
namespace GPU {

bool DriverBatchEnabled() { return BatchTracer::Enabled(); }

int DriverFiberCount(long long nTasks,int hostThreads) {
  BatchTracer::Get().SetHostThreads(hostThreads);
  const long long preferred=std::max(1,BatchTracer::Get().PreferredFibers());
  if (nTasks<=0) return static_cast<int>(preferred);
  return static_cast<int>(std::max(1LL,std::min(nTasks,preferred)));
}

void RunDriverFibers(int nWorkers,int hostThreads,
                     const std::function<void(int)>& worker,
                     const std::function<void()>& service,
                     const char* label) {
  BatchTracer& tracer=BatchTracer::Get();
  tracer.SetHostThreads(hostThreads);
  Fiber::SetServiceHook(service);
  const auto t0=std::chrono::steady_clock::now();
  try {
    Fiber::ParallelFor(0,nWorkers,nWorkers,[&](long long w,int) { worker(static_cast<int>(w)); });
  }
  catch (...) {
    Fiber::SetServiceHook(nullptr);
    throw;
  }
  Fiber::SetServiceHook(nullptr);
  const char* verbose=std::getenv("AMPS_EARTH_ENGINE_VERBOSE");
  if (verbose!=nullptr && std::string(verbose)=="1") {
    const double dt=std::chrono::duration<double>(std::chrono::steady_clock::now()-t0).count();
    const Fiber::Stats st=Fiber::LastStats();
    std::printf("[earth-engine] %s: workers=%d hostThreads=%d wall=%.3fs fibers=%lld "
                "switches=%lld idle=%lld inFibers=%.3fs inProgress=%.3fs inService=%.3fs\n",
                label!=nullptr ? label : "",nWorkers,hostThreads,
                dt,st.fibersCreated,st.switches,st.idleCalls,st.fiberSeconds,st.idleSeconds,
                st.serviceSeconds);
    tracer.PrintStatistics(label);
  }
}

void FiberLoop(long long begin,long long end,
               const std::function<void(long long,int)>& body) {
  if (end<=begin) return;
  if (!Fiber::InFiber()) {
    for (long long i=begin;i<end;++i) body(i,0);
    return;
  }
  Fiber::ParallelFor(begin,end,static_cast<int>(std::min<long long>(end-begin,1<<20)),body);
}

} // namespace GPU
} // namespace Geospace
