#include "test_framework.hpp"

#include <swcme_sep_interface.hpp>

#include <algorithm>
#include <array>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <new>
#include <string>
#include <vector>

namespace pst08_allocation_probe {

// The allocation probe is enabled only around already warmed production hot
// calls.  Global new/delete are overridden in this validation executable so the
// test observes allocations made anywhere below the public API, including an
// accidental allocation added by a future digest/status/helper refactor.
std::atomic<bool> enabled{false};
std::atomic<std::size_t> allocations{0};

void record() noexcept {
  if (enabled.load(std::memory_order_relaxed))
    allocations.fetch_add(1,std::memory_order_relaxed);
}

}  // namespace pst08_allocation_probe

// These replacements are validation instrumentation, not production memory
// management.  They preserve the normal malloc/free behavior and cover scalar,
// array, sized, and C++17 aligned forms so zero means the complete hot call was
// allocation-free rather than only one spelling of operator new.
void* operator new(std::size_t size) {
  pst08_allocation_probe::record();
  if (void* pointer=std::malloc(size==0 ? 1 : size)) return pointer;
  throw std::bad_alloc();
}

void* operator new[](std::size_t size) {
  return ::operator new(size);
}

void operator delete(void* pointer) noexcept { std::free(pointer); }
void operator delete[](void* pointer) noexcept { std::free(pointer); }
void operator delete(void* pointer,std::size_t) noexcept { std::free(pointer); }
void operator delete[](void* pointer,std::size_t) noexcept { std::free(pointer); }

void* operator new(std::size_t size,std::align_val_t alignment) {
  pst08_allocation_probe::record();
  void* pointer=nullptr;
  const std::size_t actual_size=size==0 ? 1 : size;
  if (posix_memalign(&pointer,static_cast<std::size_t>(alignment),actual_size)!=0)
    throw std::bad_alloc();
  return pointer;
}

void* operator new[](std::size_t size,std::align_val_t alignment) {
  return ::operator new(size,alignment);
}

void operator delete(void* pointer,std::align_val_t) noexcept {
  std::free(pointer);
}
void operator delete[](void* pointer,std::align_val_t) noexcept {
  std::free(pointer);
}
void operator delete(void* pointer,std::size_t,std::align_val_t) noexcept {
  std::free(pointer);
}
void operator delete[](void* pointer,std::size_t,std::align_val_t) noexcept {
  std::free(pointer);
}

namespace {

using Clock=std::chrono::steady_clock;

// A volatile checksum makes every benchmark result observable and prevents an
// optimizing compiler from deleting a successful query whose arrays are reused
// on the next repetition.  It is updated outside timed production calls where
// possible so it does not inflate the measured ownership cost.
volatile double benchmark_checksum=0.0;

struct TimingSummary {
  double median_ns=0.0;
  double p95_ns=0.0;
};

// Warm up code/data pages, collect independent aggregate samples, and report
// median plus nearest-rank 95th percentile as required by the roadmap.  Each
// callable controls its own inner loop, allowing cheap scalar calls to run long
// enough that timer resolution is not the dominant measurement.
template <typename Function>
TimingSummary benchmark(Function&& function,
                        std::size_t warmup_count=5,
                        std::size_t sample_count=21) {
  for (std::size_t i=0;i<warmup_count;++i) function();
  std::vector<double> samples;
  samples.reserve(sample_count);
  for (std::size_t i=0;i<sample_count;++i) {
    const auto begin=Clock::now();
    function();
    const auto end=Clock::now();
    samples.push_back(std::chrono::duration<double,std::nano>(end-begin).count());
  }
  std::sort(samples.begin(),samples.end());
  const std::size_t p95_index=static_cast<std::size_t>(
      std::ceil(0.95*static_cast<double>(samples.size())))-1;
  return {samples[samples.size()/2],samples[p95_index]};
}

template <typename Function>
std::size_t count_allocations(Function&& function) {
  pst08_allocation_probe::allocations.store(0,std::memory_order_relaxed);
  pst08_allocation_probe::enabled.store(true,std::memory_order_release);
  function();
  pst08_allocation_probe::enabled.store(false,std::memory_order_release);
  return pst08_allocation_probe::allocations.load(std::memory_order_relaxed);
}

void print_timing(const char* label,const TimingSummary& timing,
                  double divisor=1.0,const char* unit="ns/call") {
  std::cout << "  " << std::left << std::setw(38) << label << std::right
            << " median=" << std::setw(11) << std::fixed << std::setprecision(2)
            << timing.median_ns/divisor << ' ' << unit
            << " p95=" << std::setw(11) << timing.p95_ns/divisor << ' ' << unit
            << '\n';
}

}  // namespace

// PST08 measures the complete model-identity/configuration/integrity guard and
// corrected scalar, batch, and AMPS paths independently of prepare_step().  The
// emulated legacy topology repeats the same current validation the number of
// times removed by PST08, providing an on-host before/after comparison without
// exposing an unsafe unvalidated production API or maintaining old physics.
void test_pst08(swcme_test::Context& context) {
  std::cout << "PST08 state ownership performance\n";

  // Coverage instrumentation uses non-production optimization and inserts a
  // counter update around basic blocks.  Timing that binary against the
  // optimized PST08 budgets would test gcov overhead, not SWCME performance.
  // COV01 still compiles this complete translation unit and records the
  // explicit exclusion branch; the dedicated `pst08-performance` target and
  // ordinary optimized registry remain the authoritative performance gates.
  const char* coverage_child=std::getenv("SWCME_COV01_CHILD");
  if (coverage_child!=nullptr && coverage_child[0]=='1' &&
      coverage_child[1]=='\0') {
    context.expect_true(true,
        "PST08 timing excluded from unoptimized coverage execution");
    return;
  }

  constexpr std::size_t scalar_inner=256;
  constexpr std::size_t three_scalar_inner=8;
  constexpr std::size_t one_batch_size=16384;
  constexpr std::size_t three_batch_size=256;

  swcme1d::Params p1;
  swcme1d::Model model1(p1);
  const swcme1d::StepState step1=model1.prepare_step(3600.0);
  swcme::sep::Interface1D adapter1(p1);
  const auto adapter_step1=adapter1.prepare(3600.0);

  swcme3d::Params p3;
  p3.shape=swcme3d::ShockShape::Sphere;
  swcme3d::Model model3(p3);
  const swcme3d::StepState step3=model3.prepare_step(3600.0);
  swcme::sep::Interface3D adapter3(p3);
  const auto adapter_step3=adapter3.prepare(3600.0);

  // Preparation is timed separately and never folded into per-point query
  // throughput.  Each iteration uses a fresh owner because PST01 intentionally
  // freezes a model after its first successful prepared state.
  const TimingSummary prepare1=benchmark([&] {
    swcme1d::Model model(p1);
    const auto state=model.prepare_step(3600.0);
    benchmark_checksum+=state.r_sh_m*1.0e-20;
  });
  const TimingSummary prepare3=benchmark([&] {
    swcme3d::Model model(p3);
    const auto state=model.prepare_step(3600.0);
    benchmark_checksum+=state.r_sh_m*1.0e-20;
  });

  const TimingSummary validation1=benchmark([&] {
    unsigned code_sum=0;
    for (std::size_t i=0;i<scalar_inner;++i)
      code_sum+=static_cast<unsigned>(
          model1.validate_prepared_state(step1,"PST08 1-D validation").code);
    benchmark_checksum+=static_cast<double>(code_sum);
  });
  const TimingSummary validation3=benchmark([&] {
    unsigned code_sum=0;
    for (std::size_t i=0;i<scalar_inner;++i)
      code_sum+=static_cast<unsigned>(
          model3.validate_prepared_state(step3,"PST08 3-D validation").code);
    benchmark_checksum+=static_cast<double>(code_sum);
  });

  double r1=swcme::constants::AU_M,n1=0.0,v1=0.0,Br=0.0,Bphi=0.0;
  double Bmag=0.0,div1=0.0;
  const TimingSummary scalar1=benchmark([&] {
    unsigned code_sum=0;
    for (std::size_t i=0;i<scalar_inner;++i)
      code_sum+=static_cast<unsigned>(model1.evaluate_radii_with_B_div_checked(
          step1,&r1,&n1,&v1,&Br,&Bphi,&Bmag,&div1,1).code);
    benchmark_checksum+=n1*1.0e-12+static_cast<double>(code_sum);
  });

  // Emulate the pre-PST08 composite topology: the full evaluator authenticated
  // once and its nested fast evaluator authenticated the same state again.
  const TimingSummary legacy_scalar1=benchmark([&] {
    unsigned code_sum=0;
    for (std::size_t i=0;i<scalar_inner;++i) {
      code_sum+=static_cast<unsigned>(
          model1.validate_prepared_state(step1,"PST08 legacy nested 1-D").code);
      code_sum+=static_cast<unsigned>(model1.evaluate_radii_with_B_div_checked(
          step1,&r1,&n1,&v1,&Br,&Bphi,&Bmag,&div1,1).code);
    }
    benchmark_checksum+=n1*1.0e-12+static_cast<double>(code_sum);
  });

  swcme::sep::BackgroundState background1;
  const TimingSummary adapter_timing1=benchmark([&] {
    unsigned code_sum=0;
    for (std::size_t i=0;i<scalar_inner;++i)
      code_sum+=static_cast<unsigned>(
          adapter1.evaluate_background(adapter_step1,r1,background1).code);
    benchmark_checksum+=background1.density_m3*1.0e-12+
                        static_cast<double>(code_sum);
  });

  // Before PST08 the adapter, full evaluator, and nested fast evaluator each
  // performed the complete validation.  Two explicit guards plus the corrected
  // one-call adapter reproduce that removed cost using the exact current hash.
  const TimingSummary legacy_adapter1=benchmark([&] {
    unsigned code_sum=0;
    for (std::size_t i=0;i<scalar_inner;++i) {
      code_sum+=static_cast<unsigned>(adapter1.model().validate_prepared_state(
          adapter_step1,"PST08 legacy adapter outer").code);
      code_sum+=static_cast<unsigned>(adapter1.model().validate_prepared_state(
          adapter_step1,"PST08 legacy full nested").code);
      code_sum+=static_cast<unsigned>(
          adapter1.evaluate_background(adapter_step1,r1,background1).code);
    }
    benchmark_checksum+=background1.density_m3*1.0e-12+
                        static_cast<double>(code_sum);
  });

  std::vector<double> radii(one_batch_size),density1(one_batch_size),
      velocity1(one_batch_size);
  for (std::size_t i=0;i<one_batch_size;++i)
    radii[i]=(0.25+1.25*static_cast<double>(i)/
        static_cast<double>(one_batch_size-1))*swcme::constants::AU_M;
  const TimingSummary batch1=benchmark([&] {
    const swcme::ModelStatus status=model1.evaluate_radii_fast_checked(
        step1,radii.data(),density1.data(),velocity1.data(),one_batch_size);
    benchmark_checksum+=density1[one_batch_size/2]*1.0e-12+
                        static_cast<double>(status.code);
  });

  const double x=swcme::constants::AU_M,y=0.0,z=0.0;
  double n3=0.0,vx=0.0,vy=0.0,vz=0.0,bx=0.0,by=0.0,bz=0.0,div3=0.0;
  const TimingSummary scalar3=benchmark([&] {
    unsigned code_sum=0;
    for (std::size_t i=0;i<three_scalar_inner;++i)
      code_sum+=static_cast<unsigned>(model3.evaluate_cartesian_with_B_div_checked(
          step3,&x,&y,&z,&n3,&vx,&vy,&vz,&bx,&by,&bz,&div3,1).code);
    benchmark_checksum+=n3*1.0e-12+static_cast<double>(code_sum);
  });

  swcme::sep::BackgroundState background3;
  const std::array<double,3> position3{{x,y,z}};
  const TimingSummary adapter_timing3=benchmark([&] {
    unsigned code_sum=0;
    for (std::size_t i=0;i<three_scalar_inner;++i)
      code_sum+=static_cast<unsigned>(
          adapter3.evaluate_background(adapter_step3,position3,background3).code);
    benchmark_checksum+=background3.density_m3*1.0e-12+
                        static_cast<double>(code_sum);
  });

  std::vector<double> x3(three_batch_size),y3(three_batch_size),
      z3(three_batch_size),density3(three_batch_size),vx3(three_batch_size),
      vy3(three_batch_size),vz3(three_batch_size);
  for (std::size_t i=0;i<three_batch_size;++i) {
    const double angle=2.0*swcme::constants::PI*
        static_cast<double>(i)/static_cast<double>(three_batch_size);
    x3[i]=swcme::constants::AU_M*std::cos(angle);
    y3[i]=swcme::constants::AU_M*std::sin(angle);
    z3[i]=0.1*swcme::constants::AU_M*std::sin(2.0*angle);
  }
  const TimingSummary batch3=benchmark([&] {
    const swcme::ModelStatus status=model3.evaluate_cartesian_fast_checked(
        step3,x3.data(),y3.data(),z3.data(),density3.data(),vx3.data(),
        vy3.data(),vz3.data(),three_batch_size);
    benchmark_checksum+=density3[three_batch_size/2]*1.0e-12+
                        static_cast<double>(status.code);
  },3,15);

  // Reproduce the former per-point nested shock and geometry guards.  This is
  // a comparison metric only; correctness gates never expose the private
  // validated kernels or use this emulation to produce reference physics.
  const TimingSummary legacy_batch3=benchmark([&] {
    unsigned code_sum=0;
    for (std::size_t i=0;i<2*three_batch_size;++i)
      code_sum+=static_cast<unsigned>(model3.validate_prepared_state(
          step3,"PST08 legacy 3-D per-sample validation").code);
    const swcme::ModelStatus status=model3.evaluate_cartesian_fast_checked(
        step3,x3.data(),y3.data(),z3.data(),density3.data(),vx3.data(),
        vy3.data(),vz3.data(),three_batch_size);
    benchmark_checksum+=density3[three_batch_size/2]*1.0e-12+
                        static_cast<double>(code_sum)+
                        static_cast<double>(status.code);
  },3,15);

  std::cout << "  build contract: optimized C++17, steady_clock, warm-up=5, "
               "samples=21 (3-D batch samples=15)\n";
  print_timing("prepare_step 1-D",prepare1);
  print_timing("prepare_step 3-D",prepare3);
  print_timing("state validation 1-D",validation1,
               static_cast<double>(scalar_inner));
  print_timing("state validation 3-D",validation3,
               static_cast<double>(scalar_inner));
  print_timing("corrected scalar full 1-D",scalar1,
               static_cast<double>(scalar_inner));
  print_timing("emulated pre-PST08 scalar 1-D",legacy_scalar1,
               static_cast<double>(scalar_inner));
  print_timing("corrected AMPS scalar 1-D",adapter_timing1,
               static_cast<double>(scalar_inner));
  print_timing("emulated pre-PST08 AMPS 1-D",legacy_adapter1,
               static_cast<double>(scalar_inner));
  print_timing("corrected batch fast 1-D",batch1,
               static_cast<double>(one_batch_size),"ns/point");
  print_timing("corrected scalar full 3-D",scalar3,
               static_cast<double>(three_scalar_inner));
  print_timing("corrected AMPS scalar 3-D",adapter_timing3,
               static_cast<double>(three_scalar_inner));
  print_timing("corrected batch fast 3-D",batch3,
               static_cast<double>(three_batch_size),"ns/point");
  print_timing("emulated pre-PST08 batch 3-D",legacy_batch3,
               static_cast<double>(three_batch_size),"ns/point");

  // Absolute budgets are intentionally conservative enough for shared CI
  // hosts, yet several orders tighter than a production field-update cadence.
  // Relative checks independently ensure the removed validation topology gives
  // a material scalar improvement on the exact machine running the test.
  context.expect_true(validation1.p95_ns/scalar_inner<=10000.0 &&
                          validation3.p95_ns/scalar_inner<=10000.0,
                      "complete state validation stays below 10 us p95");
  context.expect_true(scalar1.p95_ns/scalar_inner<=20000.0 &&
                          adapter_timing1.p95_ns/scalar_inner<=25000.0,
                      "1-D scalar direct/AMPS queries meet p95 budgets");
  context.expect_true(scalar3.p95_ns/three_scalar_inner<=500000.0 &&
                          adapter_timing3.p95_ns/three_scalar_inner<=600000.0,
                      "3-D scalar direct/AMPS queries meet p95 budgets");
  context.expect_true(batch1.p95_ns<=5.0e6 && batch3.p95_ns<=1.0e8,
                      "representative corrected batch runtimes meet targets");
  // validation1/3 summaries contain scalar_inner calls per timing sample;
  // normalize to one public-call guard before comparing with one batch call.
  context.expect_true(
                          (validation1.median_ns/scalar_inner)/
                                  batch1.median_ns<=0.05 &&
                          (validation3.median_ns/scalar_inner)/
                                  batch3.median_ns<=0.05,
                      "one validation is at most 5% of representative batch time");
  context.expect_true(scalar1.median_ns<=0.80*legacy_scalar1.median_ns,
                      "1-D composite removes material duplicate validation cost");
  context.expect_true(adapter_timing1.median_ns<=0.60*legacy_adapter1.median_ns,
                      "AMPS hot path removes material duplicate validation cost");

  // Allocation checks run after all vectors, timing samples, and iostream
  // buffers are warm.  Every destination is preallocated and each public call
  // must complete without invoking any form of operator new.
  const std::size_t validation_allocations=count_allocations([&] {
    benchmark_checksum+=static_cast<double>(
        model1.validate_prepared_state(step1,"PST08 allocation validation").code);
  });
  const std::size_t one_allocations=count_allocations([&] {
    benchmark_checksum+=static_cast<double>(model1.evaluate_radii_fast_checked(
        step1,radii.data(),density1.data(),velocity1.data(),one_batch_size).code);
    benchmark_checksum+=static_cast<double>(
        adapter1.evaluate_background(adapter_step1,r1,background1).code);
  });
  const std::size_t three_allocations=count_allocations([&] {
    benchmark_checksum+=static_cast<double>(model3.evaluate_cartesian_fast_checked(
        step3,x3.data(),y3.data(),z3.data(),density3.data(),vx3.data(),
        vy3.data(),vz3.data(),three_batch_size).code);
    benchmark_checksum+=static_cast<double>(
        adapter3.evaluate_background(adapter_step3,position3,background3).code);
  });
  std::cout << "  hot-path allocations: validation=" << validation_allocations
            << " 1-D=" << one_allocations << " 3-D=" << three_allocations
            << '\n';
  context.expect_true(validation_allocations==0 && one_allocations==0 &&
                          three_allocations==0,
                      "ownership validation and direct/AMPS hot paths allocate zero memory");
  context.expect_true(std::isfinite(benchmark_checksum),
                      "benchmark checksum remains finite");
}
