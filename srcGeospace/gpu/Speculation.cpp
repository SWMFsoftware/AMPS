//======================================================================================
// Speculation.cpp -- see Speculation.h.
//======================================================================================

#include "Speculation.h"

#include <cstdlib>
#include <cstring>
#include <unordered_map>

namespace Geospace {
namespace GPU {

namespace {
// Deferred-action list of each fiber that currently evaluates a speculative candidate.
// All fibers of a rank run on its main thread, so no locking is needed.
std::unordered_map<Fiber::FiberRec*,std::vector<std::function<void()>>*> gDeferred;
} // namespace

bool SpeculationEnabled() {
  static int enabled=-1;
  if (enabled<0) {
    const char* e=std::getenv("AMPS_EARTH_SPECULATION");
    if (e!=nullptr && (std::strcmp(e,"0")==0 || std::strcmp(e,"1")==0)) enabled=(e[0]=='1') ? 1 : 0;
    else enabled=(DriverBatchEnabled() && BatchTracer::Get().Accelerator()) ? 1 : 0;
  }
  return enabled==1 && DriverBatchEnabled();
}

int SpeculativeDepth() {
  static int depth=0;
  if (depth==0) {
    depth=4;
    if (const char* e=std::getenv("AMPS_EARTH_SPECULATIVE_DEPTH")) {
      const int v=std::atoi(e);
      if (v>=1 && v<=8) depth=v;
    }
  }
  return depth;
}

namespace Detail {
void SetDeferredList(std::vector<std::function<void()>>* list) {
  Fiber::FiberRec* f=Fiber::Current();
  if (f==nullptr) return;
  if (list!=nullptr) gDeferred[f]=list;
  else gDeferred.erase(f);
}
} // namespace Detail

void DeferOrRun(std::function<void()> action) {
  Fiber::FiberRec* f=Fiber::Current();
  if (f!=nullptr) {
    auto it=gDeferred.find(f);
    if (it!=gDeferred.end()) {
      it->second->push_back(std::move(action));
      return;
    }
  }
  action();
}

} // namespace GPU
} // namespace Geospace
