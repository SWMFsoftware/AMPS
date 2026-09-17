#include "test_framework.hpp"

#include <swcme_sep_interface.hpp>

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <string>
#include <thread>
#include <vector>

namespace {

using Clock=std::chrono::steady_clock;
volatile double perf01_checksum=0.0;

struct Timing {
  double median_ms=0.0;
  double p95_ms=0.0;
};

// PERF01 uses aggregate operations long enough to dominate timer resolution,
// warms each path, and reports robust order statistics rather than one noisy
// wall-clock observation.  Nine samples keep the registered gate practical;
// the dedicated output retains both median and nearest-rank p95.
template <class Function>
Timing measure(Function&& function,std::size_t warmups=2,std::size_t samples=9) {
  for (std::size_t i=0;i<warmups;++i) function();
  std::vector<double> milliseconds;
  milliseconds.reserve(samples);
  for (std::size_t i=0;i<samples;++i) {
    const auto begin=Clock::now();
    function();
    const auto end=Clock::now();
    milliseconds.push_back(
        std::chrono::duration<double,std::milli>(end-begin).count());
  }
  std::sort(milliseconds.begin(),milliseconds.end());
  const std::size_t p95=static_cast<std::size_t>(
      std::ceil(0.95*static_cast<double>(samples)))-1;
  return {milliseconds[samples/2],milliseconds[p95]};
}

void report(const char* name,const Timing& timing) {
  std::cout << "  " << std::left << std::setw(34) << name << std::right
            << " median=" << std::setw(10) << std::fixed
            << std::setprecision(3) << timing.median_ms << " ms"
            << " p95=" << std::setw(10) << timing.p95_ms << " ms\n";
}

}  // namespace

void test_perf01(swcme_test::Context& context) {
  std::cout << "PERF01 performance and scaling guardrails\n";

  // Sanitizer and coverage executables deliberately perturb cost by large,
  // tool-dependent factors.  They compile this code but acknowledge the timing
  // exclusion; only the normal parent launches the pinned optimized child.
  if (std::getenv("SWCME_SAN01_CHILD")!=nullptr ||
      std::getenv("SWCME_COV01_CHILD")!=nullptr) {
    context.expect_true(true,"PERF01 excluded from instrumented execution");
    return;
  }
  const char* child=std::getenv("SWCME_PERF01_CHILD");
  if (child==nullptr || child[0]!='1' || child[1]!='\0') {
    context.expect_true(
        std::system("make --no-print-directory perf01-performance")==0,
        "pinned optimized PERF01 child meets all budgets");
    return;
  }

  constexpr double AU=swcme::constants::AU_M;
  swcme1d::Params p1;
  p1.kinematics_mode=swcme::kinematics::Mode::Ballistic;
  p1.r0_Rs=20.0;
  p1.V0_sh_kms=1450.0;
  swcme3d::Params p3;
  p3.shape=swcme3d::ShockShape::SSE;
  p3.half_width_rad=55.0*swcme::constants::PI/180.0;
  p3.kinematics_mode=p1.kinematics_mode;
  p3.r0_Rs=p1.r0_Rs;
  p3.V0_sh_kms=p1.V0_sh_kms;
  const swcme1d::Model model1(p1);
  const swcme3d::Model model3(p3);
  const auto step1=model1.prepare_step(4.0*3600.0);
  const auto step3=model3.prepare_step(4.0*3600.0);
  const swcme::sep::Interface3D adapter3(p3);
  const auto adapter_step3=adapter3.prepare(4.0*3600.0);

  const Timing prepare=measure([&] {
    swcme1d::Model fresh1(p1);
    swcme3d::Model fresh3(p3);
    const auto a=fresh1.prepare_step(4.0*3600.0);
    const auto b=fresh3.prepare_step(4.0*3600.0);
    perf01_checksum+=(a.r_sh_m+b.r_sh_m)*1.0e-20;
  });

  constexpr std::size_t one_count=8192;
  std::vector<double> radii(one_count),density(one_count),speed(one_count);
  for (std::size_t i=0;i<one_count;++i)
    radii[i]=(0.25+1.25*static_cast<double>(i)/(one_count-1))*AU;
  const Timing batch1=measure([&] {
    const auto status=model1.evaluate_radii_fast_checked(
        step1,radii.data(),density.data(),speed.data(),one_count);
    perf01_checksum+=density[one_count/2]*1.0e-12+
                     static_cast<double>(status.code);
  });

  const std::array<double,3> shock_direction{{1.0,0.0,0.0}};
  const Timing shock=measure([&] {
    for (int i=0;i<32;++i) {
      swcme3d::LocalShockState state;
      const auto status=model3.shock_state_direction_checked(
          step3,shock_direction.data(),state);
      perf01_checksum+=state.compression+static_cast<double>(status.code);
    }
  });

  const Timing mesh=measure([&] {
    const auto value=model3.build_shock_mesh(step3,18,36);
    swcme3d::TriMetrics metrics;
    model3.compute_triangle_metrics(value,metrics);
    const auto box=model3.default_apex_box(step3,0.15,12);
    perf01_checksum+=static_cast<double>(metrics.area.size()+box.Ni);
  });

  const Timing source=measure([&] {
    swcme::sep::SourceSurface surface;
    const auto status=adapter3.build_shock_surface_source(
        adapter_step3,12,24,surface);
    perf01_checksum+=static_cast<double>(surface.patch_count)+
                     static_cast<double>(status.code);
  },1,7);

  swcme3d::ConnectivityOptions connectivity_options;
  connectivity_options.scan_intervals=128;
  const std::array<double,3> observer{{AU,0.0,0.0}};
  const Timing connectivity=measure([&] {
    for (int i=0;i<8;++i) {
      const auto state=model3.observer_connectivity(
          step3,observer.data(),connectivity_options);
      perf01_checksum+=static_cast<double>(state.roots.size());
    }
  },1,7);

  // The scaling workload writes disjoint result slots and performs the final
  // checksum in canonical order.  Therefore measured differences reflect
  // scheduler throughput rather than an unsafe shared reduction.
  constexpr std::size_t point_count=192;
  std::vector<double> x(point_count),y(point_count),z(point_count);
  for (std::size_t i=0;i<point_count;++i) {
    const double angle=2.0*swcme::constants::PI*i/point_count;
    x[i]=1.15*AU*std::cos(angle);
    y[i]=1.15*AU*std::sin(angle);
    z[i]=0.08*AU*std::sin(3.0*angle);
  }
  const std::array<std::size_t,4> thread_counts{{1,2,4,8}};
  std::array<Timing,4> scaling{};
  for (std::size_t level=0;level<thread_counts.size();++level) {
    const std::size_t nthreads=thread_counts[level];
    scaling[level]=measure([&] {
      std::vector<double> n(point_count),vx(point_count),vy(point_count),
          vz(point_count);
      std::vector<std::thread> workers;
      for (std::size_t thread=0;thread<nthreads;++thread) {
        workers.emplace_back([&,thread] {
          const std::size_t begin=point_count*thread/nthreads;
          const std::size_t end=point_count*(thread+1)/nthreads;
          if (begin==end) return;
          model3.evaluate_cartesian_fast_checked(
              step3,x.data()+begin,y.data()+begin,z.data()+begin,n.data()+begin,
              vx.data()+begin,vy.data()+begin,vz.data()+begin,end-begin);
        });
      }
      for (auto& worker:workers) worker.join();
      for (double value:n) perf01_checksum+=value*1.0e-16;
    },1,7);
  }

  report("prepare 1-D + 3-D",prepare);
  report("1-D batch 8192",batch1);
  report("32 directional shock solves",shock);
  report("18x36 mesh + output preflight",mesh);
  report("12x24 integrated source",source);
  report("8 connectivity queries",connectivity);
  for (std::size_t i=0;i<thread_counts.size();++i) {
    const std::string name="192 fields / "+std::to_string(thread_counts[i])+
                           " thread(s)";
    report(name.c_str(),scaling[i]);
  }

  // These frozen budgets are deliberately broad enough for shared CI hosts,
  // but tight enough to catch an accidental O(N^2) loop, per-point allocation,
  // repeated state preparation, or complete loss of parallel work sharing.
  context.expect_true(prepare.p95_ms<100.0,"prepare-step p95 budget");
  context.expect_true(batch1.p95_ms<100.0,"1-D batch p95 budget");
  context.expect_true(shock.p95_ms<3000.0,"shock-solver p95 budget");
  context.expect_true(mesh.p95_ms<3000.0,"mesh/output-preflight p95 budget");
  context.expect_true(source.p95_ms<5000.0,"source integration p95 budget");
  context.expect_true(connectivity.p95_ms<5000.0,
                      "connectivity p95 budget");
  for (std::size_t i=1;i<scaling.size();++i)
    context.expect_true(scaling[i].median_ms<2.5*scaling[0].median_ms,
        std::to_string(thread_counts[i])+"-thread scaling avoids severe regression");
  context.expect_true(std::isfinite(perf01_checksum),
                      "all benchmark outputs remain observable and finite");

  std::cout << "  compiler=" << __VERSION__
            << " hardware_threads=" << std::thread::hardware_concurrency()
            << " checksum=" << std::scientific << perf01_checksum << '\n';
}
