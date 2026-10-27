#ifndef _SRC_GEOSPACE_GPU_TASK_STREAM_H_
#define _SRC_GEOSPACE_GPU_TASK_STREAM_H_

//======================================================================================
// TaskStream.h -- one rank's task stream for the batch-mode drivers.
//
// The application drivers assign global task ids to MPI ranks with the DYNAMIC (RMA
// fetch-and-add chunks), BLOCK_CYCLIC or STATIC scheduler and then compute each chunk
// with a worker pool.  In batch mode the rank's fibers consume a single stream that
// yields exactly the ids the rank would have processed, fetching DYNAMIC chunks on
// demand.  RunStreamBatch() computes the tasks with fibers and then accumulates the
// results in increasing task id, which is the order used by the original schedulers,
// so order-dependent floating-point sums are unchanged.
//======================================================================================

#include "DriverSupport.h"
#include "../3d/Mode3DParallel.h"

#include <algorithm>
#include <atomic>
#include <functional>
#include <memory>
#include <utility>
#include <vector>

namespace Geospace {
namespace GPU {

class TaskStream {
 public:
  // `comm` must be the communicator the original driver scheduled on (the gridless
  // drivers use MPI_COMM_WORLD; Mode3D drivers use MPI_GLOBAL_COMMUNICATOR, which
  // differs in EPOCH_PARALLEL/SWMF runs).
  TaskStream(Geospace::Mode3D::MpiScheduler scheduler,long long totalTasks,long long chunk,
             int mpiRank,int mpiSize,const char* label,MPI_Comm comm=MPI_COMM_WORLD)
      : scheduler_(scheduler), total_(totalTasks), rank_(mpiRank), size_(mpiSize) {
    if (scheduler_==Geospace::Mode3D::MpiScheduler::DYNAMIC && total_>0)
      dynamic_.reset(new Geospace::Mode3D::DynamicMpiLocationScheduler(
          comm,total_,chunk,label));
    rankCount_=(total_<=rank_) ? 0LL : 1LL+(total_-1LL-(long long)rank_)/(long long)size_;
    staticBegin_=(total_*(long long)rank_)/(long long)size_;
    staticEnd_=(total_*(long long)(rank_+1))/(long long)size_;
  }

  // Next global task id of this rank, or -1 when the rank's share is exhausted.
  long long Next() {
    if (total_<=0) return -1;
    if (dynamic_) {
      if (next_>=end_) {
        if (exhausted_) return -1;
        const long long start=dynamic_->FetchNextChunkStart();
        if (start>=total_) { exhausted_=true; return -1; }
        next_=start;
        end_=std::min(start+dynamic_->ChunkSize(),total_);
      }
      return next_++;
    }
    if (scheduler_==Geospace::Mode3D::MpiScheduler::BLOCK_CYCLIC) {
      if (counter_>=rankCount_) return -1;
      return (long long)rank_+(counter_++)*(long long)size_;
    }
    if (staticBegin_+counter_>=staticEnd_) return -1;
    return staticBegin_+(counter_++);
  }

 private:
  Geospace::Mode3D::MpiScheduler scheduler_;
  long long total_;
  int rank_,size_;
  std::unique_ptr<Geospace::Mode3D::DynamicMpiLocationScheduler> dynamic_;
  long long next_=0,end_=0,counter_=0,rankCount_=0,staticBegin_=0,staticEnd_=0;
  bool exhausted_=false;
};

// Compute every task of `stream` with fibers, then call accumulate() in task-id order.
// progress(delta) is called from the scheduler with the number of newly completed tasks.
template <class Result>
void RunStreamBatch(TaskStream& stream,int hostThreads,
                    const std::function<Result(long long taskId,int fiber)>& compute,
                    const std::function<void(const Result&)>& accumulate,
                    const std::function<void(long long delta)>& progress,
                    const char* label) {
  const int nFibers=DriverFiberCount(0,hostThreads);
  std::vector<std::pair<long long,Result>> results;
  std::atomic<long long> completed(0);
  long long published=0;
  auto service=[&]() {
    const long long c=completed.load(std::memory_order_acquire);
    if (progress) progress(c-published);
    published=c;
  };
  RunDriverFibers(nFibers,hostThreads,[&](int fiber) {
    for (;;) {
      const long long id=stream.Next();
      if (id<0) break;
      Result r=compute(id,fiber);
      results.emplace_back(id,std::move(r));
      completed.fetch_add(1,std::memory_order_release);
    }
  },service,label);
  service();
  std::sort(results.begin(),results.end(),
            [](const std::pair<long long,Result>& a,const std::pair<long long,Result>& b) {
              return a.first<b.first;
            });
  for (const auto& r:results) accumulate(r.second);
}

} // namespace GPU
} // namespace Geospace

#endif
