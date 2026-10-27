//======================================================================================
// HostThreadPool.cpp -- see HostThreadPool.h.
//======================================================================================

#include "HostThreadPool.h"

#include <algorithm>
#include <atomic>
#include <chrono>
#include <condition_variable>
#include <cstdlib>
#include <cstring>
#include <mutex>
#include <thread>
#include <vector>

#if defined(__x86_64__) || defined(__i386__)
#include <immintrin.h>
#define EG_CPU_RELAX() _mm_pause()
#else
#define EG_CPU_RELAX() do {} while (0)
#endif

namespace Geospace {
namespace GPU {

namespace {

class Pool {
 public:
  static Pool& Get() {
    static Pool pool;
    return pool;
  }
  static Pool& GetAsync() {
    static Pool pool;
    return pool;
  }

  // Launch without the caller participating; Wait() completes the launch.
  void RunAsync(int n,int nThreads,std::function<void(int)> body) {
    asyncBody_=std::move(body);
    if (n<=0) { asyncPending_=false; return; }
    const int helpers=std::max(1,nThreads);
    EnsureWorkers(helpers);
    {
      std::lock_guard<std::mutex> lock(mutex_);
      body_=&asyncBody_;
      n_=n;
      next_.store(0,std::memory_order_relaxed);
      active_=helpers;
      remaining_.store(helpers,std::memory_order_relaxed);
      ++generation_;
      generationAtomic_.store(generation_,std::memory_order_release);
    }
    asyncPending_=true;
    cv_.notify_all();
  }

  // Spin briefly (back-to-back short launches), then block until the last worker of
  // the launch signals completion, so the waiting main thread does not occupy a core
  // during long launches.
  void Wait() {
    if (!asyncPending_) return;
    const auto t0=std::chrono::steady_clock::now();
    while (remaining_.load(std::memory_order_acquire)>0) {
      if (std::chrono::duration_cast<std::chrono::microseconds>(
              std::chrono::steady_clock::now()-t0).count()>spinUs_) {
        std::unique_lock<std::mutex> lock(doneMutex_);
        doneCv_.wait(lock,[&]() { return remaining_.load(std::memory_order_acquire)==0; });
        break;
      }
      EG_CPU_RELAX();
    }
    asyncPending_=false;
  }

  void Run(int n,int nThreads,const std::function<void(int)>& body) {
    if (n<=0) return;
    if (nThreads<=1 || n==1) {
      for (int i=0;i<n;++i) body(i);
      return;
    }
    EnsureWorkers(nThreads-1);
    const int helpers=nThreads-1;
    {
      std::lock_guard<std::mutex> lock(mutex_);
      body_=&body;
      n_=n;
      next_.store(0,std::memory_order_relaxed);
      active_=helpers;
      remaining_.store(helpers,std::memory_order_relaxed);
      ++generation_;
      generationAtomic_.store(generation_,std::memory_order_release);
    }
    cv_.notify_all();
    Work(body,n);
    // Wait for helpers that are still finishing an iteration.
    while (remaining_.load(std::memory_order_acquire)>0) EG_CPU_RELAX();
  }

 private:
  Pool() {
    if (const char* e=std::getenv("AMPS_EARTH_SPIN_US")) spinUs_=std::atoi(e);
  }
  ~Pool() {
    {
      std::lock_guard<std::mutex> lock(mutex_);
      stop_=true;
      ++generation_;
      generationAtomic_.store(generation_,std::memory_order_release);
    }
    cv_.notify_all();
    for (std::thread& t:threads_) t.join();
  }

  void Work(const std::function<void(int)>& body,int n) {
    for (;;) {
      const int i=next_.fetch_add(1,std::memory_order_relaxed);
      if (i>=n) break;
      body(i);
    }
  }

  void EnsureWorkers(int count) {
    while (static_cast<int>(threads_.size())<count) {
      const int id=static_cast<int>(threads_.size());
      threads_.emplace_back([this,id]() { Loop(id); });
    }
  }

  void Loop(int id) {
    unsigned long long seen=0;
    for (;;) {
      // Short spin for back-to-back launches, then block.
      const auto t0=std::chrono::steady_clock::now();
      while (generationAtomic_.load(std::memory_order_acquire)==seen) {
        if (std::chrono::duration_cast<std::chrono::microseconds>(
                std::chrono::steady_clock::now()-t0).count()>spinUs_) break;
        EG_CPU_RELAX();
      }
      const std::function<void(int)>* body;
      int n;
      bool participate;
      {
        std::unique_lock<std::mutex> lock(mutex_);
        cv_.wait(lock,[&]() { return generation_!=seen; });
        seen=generation_;
        if (stop_) return;
        body=body_;
        n=n_;
        participate=(id<active_);
      }
      if (!participate) continue;
      Work(*body,n);
      if (remaining_.fetch_sub(1,std::memory_order_acq_rel)==1) {
        // Last participant: wake a caller blocked in Wait().
        std::lock_guard<std::mutex> lock(doneMutex_);
        doneCv_.notify_all();
      }
    }
  }

  std::vector<std::thread> threads_;
  std::mutex mutex_;
  std::condition_variable cv_;
  unsigned long long generation_=0;
  std::atomic<unsigned long long> generationAtomic_{0};
  const std::function<void(int)>* body_=nullptr;
  int n_=0;
  int active_=0;
  std::atomic<int> next_{0};
  std::atomic<int> remaining_{0};
  bool stop_=false;
  long spinUs_=50;
  std::function<void(int)> asyncBody_;
  bool asyncPending_=false;
  std::mutex doneMutex_;
  std::condition_variable doneCv_;
};

} // namespace

void HostParallelFor(int n,int nThreads,const std::function<void(int)>& body) {
  Pool::Get().Run(n,nThreads,body);
}

void HostParallelForAsync(int n,int nThreads,std::function<void(int)> body) {
  Pool::GetAsync().RunAsync(n,nThreads,std::move(body));
}

void HostParallelWait() { Pool::GetAsync().Wait(); }

bool HostRuntimeIsOpenMP() {
  static int v=-1;
  if (v<0) {
    const char* e=std::getenv("AMPS_EARTH_HOST_RUNTIME");
    v=(e!=nullptr && (std::strcmp(e,"OPENMP")==0 || std::strcmp(e,"openmp")==0)) ? 1 : 0;
  }
  return v==1;
}

} // namespace GPU
} // namespace Geospace
