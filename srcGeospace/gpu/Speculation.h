#ifndef _SRC_GEOSPACE_GPU_SPECULATION_H_
#define _SRC_GEOSPACE_GPU_SPECULATION_H_

//======================================================================================
// Speculation.h -- speculative parallel evaluation for sequential cutoff searches.
//
// The cutoff searches visit trial rigidities one at a time: UPPER_SCAN walks a grid from
// Rmax downward and stops at the first forbidden sample, the bisections evaluate one
// midpoint per step.  Each search therefore keeps at most one trajectory in flight,
// which leaves a GPU almost idle (a Mode3D deck with 168 locations has 168 trajectories
// in flight).  The helpers below evaluate the candidates that the sequential search
// may visit concurrently, as fibers, and then walk the results in the original order:
//
//   * SpeculativeEval(n, eval, onDone): eval(k) for every candidate k; onDone(k,...)
//     may cancel candidates that the search can no longer use.  A cancelled candidate
//     that is waiting for a trajectory is resumed at once (BatchTracer::Cancel) and its
//     trajectory is stopped at the next pool launch.
//   * SpeculativeBisect(lo, hi, tol, depth, allowed): the bisection
//       while (hi-lo>tol) { mid=0.5*(lo+hi); if (allowed(mid)) hi=mid; else lo=mid; }
//     with the midpoints of the next `depth` levels evaluated concurrently.
//
// The result is identical to the sequential search: every candidate it uses is
// evaluated with the same arguments (midpoints are formed with the same floating-point
// operations), and the walk applies the same decisions in the same order.  Exceptions
// raised by a candidate are stored and rethrown only when the walk reaches that
// candidate, as the sequential search would.  Side effects that must follow the
// sequential order (diagnostic accumulators) are deferred with DeferOrRun() and replayed
// for the used candidates in visiting order (RunDeferred()).
//
// Speculation is used only inside fibers and when SpeculationEnabled(): by default on
// accelerators (where idle hardware makes the extra trajectories free), never in the
// CPU build.  AMPS_EARTH_SPECULATION=1/0 forces it on/off (1 on the CPU reproduces the
// GPU scheduling for verification), AMPS_EARTH_SPECULATIVE_DEPTH=D sets the bisection
// depth (default 4: 15 midpoints per round, 4 levels resolved).
//======================================================================================

#include "BatchTracer.h"
#include "DriverSupport.h"
#include "Fiber.h"

#include <exception>
#include <functional>
#include <vector>

namespace Geospace {
namespace GPU {

bool SpeculationEnabled();
int SpeculativeDepth();

// Deferred side effects of a speculative candidate.  Outside a speculative candidate
// DeferOrRun(f) calls f() immediately.
void DeferOrRun(std::function<void()> action);

template <class T>
struct SpecOutcome {
  enum State { NotRun=0, Done=1, Cancelled=2, Failed=3 } state=NotRun;
  T value{};
  std::exception_ptr error;
  std::vector<std::function<void()>> deferred;   // side effects, in execution order

  // Replay the deferred side effects (once) and return the value; rethrow a failure.
  const T& Use() {
    if (state==Failed) std::rethrow_exception(error);
    for (auto& a:deferred) a();
    deferred.clear();
    return value;
  }
};

namespace Detail {
// Install `list` as the deferred-action list of the current fiber (nullptr removes it).
void SetDeferredList(std::vector<std::function<void()>>* list);
} // namespace Detail

// Evaluate eval(k), k=0..n-1, concurrently on fibers.  onDone(k,outcome,cancel) is called
// as each candidate finishes (in completion order); cancel(j) cancels candidate j if it
// has not finished.  Not inside a fiber, or with speculation disabled, the candidates
// run one after another in index order and a cancelled candidate is skipped.
template <class T,class Eval,class OnDone>
std::vector<SpecOutcome<T>> SpeculativeEval(int n,Eval eval,OnDone onDone) {
  std::vector<SpecOutcome<T>> out(static_cast<std::size_t>(std::max(0,n)));
  std::vector<Fiber::FiberRec*> running(out.size(),nullptr);
  std::vector<char> cancelled(out.size(),0);
  auto cancel=[&](int j) {
    if (j<0 || j>=n || cancelled[(std::size_t)j]) return;
    if (out[(std::size_t)j].state!=SpecOutcome<T>::NotRun) return;
    cancelled[(std::size_t)j]=1;
    if (running[(std::size_t)j]!=nullptr) BatchTracer::Get().Cancel(running[(std::size_t)j]);
  };
  FiberLoop(0,n,[&](long long kk,int) {
    const std::size_t k=static_cast<std::size_t>(kk);
    if (cancelled[k]) { out[k].state=SpecOutcome<T>::Cancelled; return; }
    Fiber::FiberRec* self=Fiber::Current();
    running[k]=self;
    Detail::SetDeferredList(&out[k].deferred);
    try {
      out[k].value=eval(static_cast<int>(k));
      out[k].state=SpecOutcome<T>::Done;
    }
    catch (const TraceCancelled&) {
      out[k].state=SpecOutcome<T>::Cancelled;
    }
    catch (...) {
      out[k].error=std::current_exception();
      out[k].state=SpecOutcome<T>::Failed;
    }
    Detail::SetDeferredList(nullptr);
    running[k]=nullptr;
    if (self!=nullptr) BatchTracer::Get().ClearCancel(self);
    if (cancelled[k] && out[k].state!=SpecOutcome<T>::Done) out[k].state=SpecOutcome<T>::Cancelled;
    if (out[k].state!=SpecOutcome<T>::Cancelled) onDone(static_cast<int>(k),out[k],cancel);
  });
  return out;
}

// Sequential-equivalent bisection with speculative evaluation of `depth` levels per
// round (see the header comment).  allowed(R) traces one trajectory.
template <class Allowed>
double SpeculativeBisect(double lo,double hi,double tol,int depth,Allowed allowed) {
  depth=std::max(1,depth);
  while ((hi-lo)>tol) {
    // Complete binary tree in BFS order: node q has children 2q+1 (allowed: [lo,mid])
    // and 2q+2 (forbidden: [mid,hi]).  A node whose interval is already within the
    // tolerance is never visited by the sequential loop and is not evaluated.
    const int nNodes=(1<<depth)-1;
    std::vector<double> nlo(nNodes),nhi(nNodes),nmid(nNodes);
    std::vector<char> need(nNodes,0);
    nlo[0]=lo; nhi[0]=hi;
    for (int q=0;q<nNodes;++q) {
      if (q>0) {
        const int parent=(q-1)/2;
        if (!need[parent]) continue;
        if (q==2*parent+1) { nlo[q]=nlo[parent]; nhi[q]=nmid[parent]; }
        else               { nlo[q]=nmid[parent]; nhi[q]=nhi[parent]; }
      }
      if ((nhi[q]-nlo[q])>tol) {
        need[q]=1;
        nmid[q]=0.5*(nlo[q]+nhi[q]);
      }
    }
    // Index map over the nodes that need a trajectory.
    std::vector<int> nodeOf;
    std::vector<int> slotOf(nNodes,-1);
    for (int q=0;q<nNodes;++q) if (need[q]) { slotOf[q]=(int)nodeOf.size(); nodeOf.push_back(q); }
    auto outcomes=SpeculativeEval<char>(
        static_cast<int>(nodeOf.size()),
        [&](int s) -> char { return allowed(nmid[(std::size_t)nodeOf[(std::size_t)s]]) ? 1 : 0; },
        [&](int s,SpecOutcome<char>& o,const std::function<void(int)>& cancel) {
          if (o.state!=SpecOutcome<char>::Done) return;
          // The subtree on the side the search will not take is no longer needed.
          const int q=nodeOf[(std::size_t)s];
          const int drop=o.value ? 2*q+2 : 2*q+1;
          std::vector<int> stack{drop};
          while (!stack.empty()) {
            const int r=stack.back(); stack.pop_back();
            if (r>=nNodes) continue;
            if (slotOf[r]>=0) cancel(slotOf[r]);
            stack.push_back(2*r+1); stack.push_back(2*r+2);
          }
        });
    // Walk the path exactly as the sequential loop would.
    int q=0;
    for (int level=0;level<depth && (hi-lo)>tol;++level) {
      // nmid[q]==0.5*(lo+hi): the tree was built from the same lo/hi sequence.
      SpecOutcome<char>& o=outcomes[(std::size_t)slotOf[q]];
      const bool a=o.Use()!=0;   // a node on the path is never cancelled
      if (a) { hi=nmid[q]; q=2*q+1; }
      else   { lo=nmid[q]; q=2*q+2; }
      if (q>=nNodes) break;
    }
  }
  return hi;
}

} // namespace GPU
} // namespace Geospace

#endif
