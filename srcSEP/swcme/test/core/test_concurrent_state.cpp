#include "test_framework.hpp"

#include <swcme_sep_interface.hpp>

#include <algorithm>
#include <array>
#include <atomic>
#include <condition_variable>
#include <cstdint>
#include <cstring>
#include <exception>
#include <limits>
#include <mutex>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

namespace {

// A Snapshot is an ABI-independent, exact representation of one public API
// result.  We serialize each integer, flag, string byte, and IEEE-754 double
// explicitly instead of comparing C++ object memory, whose padding bytes are
// unspecified.  Exact word-vector equality therefore implements PST04's
// strongest reproducibility criterion without confusing padding differences
// with a physics or concurrency defect.
struct Snapshot {
  std::vector<std::uint64_t> words;

  void add_uint64(std::uint64_t value) { words.push_back(value); }
  void add_bool(bool value) { add_uint64(value ? 1U : 0U); }

  void add_double(double value) {
    static_assert(sizeof(value)==sizeof(std::uint64_t),
                  "PST04 requires a 64-bit IEEE-754 double representation");
    std::uint64_t bits=0;
    std::memcpy(&bits,&value,sizeof(bits));
    add_uint64(bits);
  }

  void add_string(const char* value) {
    const std::size_t length=value ? std::strlen(value) : 0;
    add_uint64(static_cast<std::uint64_t>(length));
    for (std::size_t i=0; i<length; ++i)
      add_uint64(static_cast<unsigned char>(value[i]));
  }
};

bool operator==(const Snapshot& lhs,const Snapshot& rhs) {
  return lhs.words==rhs.words;
}

// ModelStatus is returned by value, and every field is part of the concurrent
// contract.  In particular, serializing context text rather than its pointer
// checks the diagnostic presented to a caller while avoiding address-layout
// dependence between builds.
void append_status(Snapshot& out, const swcme::ModelStatus& status) {
  out.add_uint64(static_cast<std::uint64_t>(status.code));
  out.add_uint64(static_cast<std::uint64_t>(status.sample_index));
  out.add_string(status.context);
  out.add_double(status.offending_value);
  out.add_bool(status.has_offending_value);
  out.add_uint64(status.expected_model_identity);
  out.add_uint64(status.supplied_model_identity);
  out.add_bool(status.has_model_identities);
  out.add_uint64(status.expected_configuration_digest);
  out.add_uint64(status.supplied_configuration_digest);
  out.add_bool(status.has_configuration_digests);
  out.add_uint64(status.expected_state_integrity);
  out.add_uint64(status.computed_state_integrity);
  out.add_bool(status.has_state_integrity);
  out.add_uint64(static_cast<std::uint64_t>(status.io_byte_offset));
  out.add_bool(status.has_io_byte_offset);
}

template <std::size_t N>
void append_array(Snapshot& out, const std::array<double,N>& values) {
  for (double value : values) out.add_double(value);
}

void append_primitive(Snapshot& out,
                      const swcme::shock::PrimitiveState& state) {
  out.add_double(state.rho_kg_m3);
  out.add_double(state.pressure_Pa);
  append_array(out,state.velocity_m_s);
  append_array(out,state.magnetic_T);
}

// Include the complete local shock record because it is the common dependency
// of directional diagnostics, source conversion, and connectivity.  A hidden
// shared scratch value in any of those paths would change at least one word.
void append_shock(Snapshot& out, const swcme3d::LocalShockState& state) {
  out.add_bool(state.surface_exists);
  out.add_bool(state.has_shock);
  out.add_bool(state.solver_converged);
  append_status(out,state.status);
  out.add_double(state.Rdir_m);
  for (double value : state.normal) out.add_double(value);
  out.add_double(state.Vsh_n_m_s);
  out.add_double(state.theta_Bn_rad);
  out.add_double(state.fast_speed_m_s);
  out.add_double(state.fast_mach);
  out.add_double(state.compression);
  out.add_double(state.upstream_n_m3);
  out.add_double(state.downstream_n_m3);
  append_primitive(out,state.upstream);
  append_primitive(out,state.downstream);
  out.add_double(state.mass_residual);
  out.add_double(state.normal_B_residual);
  out.add_double(state.electric_residual);
  out.add_double(state.momentum_residual);
  out.add_double(state.energy_residual);
  out.add_double(state.entropy_ratio);
}

void append_spectrum(Snapshot& out,
                     const swcme::sep::SpectrumConfig& spectrum) {
  out.add_double(spectrum.particle_mass_kg);
  out.add_uint64(static_cast<std::uint64_t>(
      static_cast<std::int64_t>(spectrum.charge_number)));
  out.add_double(spectrum.kinetic_energy_min_MeV);
  out.add_double(spectrum.kinetic_energy_max_MeV);
  out.add_double(spectrum.reference_energy_MeV);
  out.add_uint64(static_cast<std::uint64_t>(spectrum.normalization));
  out.add_double(spectrum.reference_differential_intensity_SI);
}

// Serialize every transport-facing source field, including the embedded
// status and spectrum configuration.  This prevents a concurrency regression
// from hiding behind agreement of only the headline compression or Mach value.
void append_source(Snapshot& out,
                   const swcme::sep::SEPSourceState& source) {
  append_status(out,source.status);
  out.add_bool(source.active);
  out.add_bool(source.connection_evaluated);
  out.add_bool(source.connected);
  out.add_uint64(static_cast<std::uint64_t>(source.acceleration_mode));
  out.add_double(source.time_s);
  out.add_uint64(static_cast<std::uint64_t>(source.source_id));
  append_array(out,source.position_m);
  append_array(out,source.normal);
  out.add_double(source.patch_area_m2);
  out.add_double(source.active_surface_area_m2);
  out.add_double(source.area_fraction);
  out.add_double(source.relative_patch_weight);
  out.add_double(source.compression);
  out.add_double(source.theta_Bn_rad);
  out.add_double(source.fast_mach);
  out.add_double(source.normal_speed_m_s);
  out.add_double(source.upstream_density_m3);
  out.add_double(source.upstream_B_T);
  out.add_double(source.q_phase_space);
  out.add_double(source.momentum_intensity_index);
  out.add_double(source.nonrel_energy_intensity_index);
  append_spectrum(out,source.spectrum);
  out.add_double(source.relative_source_weight_per_area);
}

void append_background(Snapshot& out,
                       const swcme::sep::BackgroundState& background) {
  append_status(out,background.status);
  append_array(out,background.position_m);
  out.add_double(background.density_m3);
  append_array(out,background.velocity_m_s);
  append_array(out,background.magnetic_T);
  out.add_double(background.magnetic_magnitude_T);
  out.add_double(background.div_velocity_s_inv);
}

// Connectivity owns dynamically allocated local vectors, so it is especially
// useful for distinguishing thread-local workspaces from accidental model- or
// state-owned scratch storage.  All roots and their complete shock records are
// included rather than comparing only the selected cobpoint.
void append_connectivity(Snapshot& out,
                         const swcme3d::ConnectivityState& connectivity) {
  out.add_uint64(static_cast<std::uint64_t>(connectivity.status));
  out.add_bool(connectivity.connected);
  for (double value : connectivity.observer_position_m) out.add_double(value);
  out.add_double(connectivity.observer_radius_m);
  // CON09 diagnostics are observable result state and must remain bitwise
  // deterministic under concurrent reuse of one prepared state, just like the
  // physical roots that follow them.
  out.add_uint64(static_cast<std::uint64_t>(
      connectivity.requested_scan_intervals));
  out.add_uint64(static_cast<std::uint64_t>(
      connectivity.achieved_scan_intervals));
  out.add_uint64(static_cast<std::uint64_t>(
      connectivity.scan_interval_budget));
  out.add_uint64(static_cast<std::uint64_t>(connectivity.roots.size()));
  out.add_uint64(static_cast<std::uint64_t>(connectivity.selected_root));
  for (const swcme3d::ConnectivityRoot& root : connectivity.roots) {
    out.add_double(root.radius_m);
    for (double value : root.position_m) out.add_double(value);
    out.add_double(root.path_length_m);
    out.add_double(root.surface_residual_m);
    append_shock(out,root.shock);
  }
}

// Test code uses exceptions only to transport an unexpected fixture failure
// out of a worker.  The main validation thread catches and records it after all
// workers join, so the non-thread-safe Context object is never shared.
void require(bool condition, const char* message) {
  if (!condition) throw std::runtime_error(message);
}

enum class Operation : std::size_t {
  OneScalarBackground,
  OneBatchBackground,
  OneShockSource,
  ThreeScalarBackground,
  ThreeBatchBackground,
  ThreeDirectionalShock,
  ThreeDirectionalSource,
  ThreeConnectivity,
  IndependentFailureDiagnostics,
  Count
};

constexpr std::size_t kOperationCount=
    static_cast<std::size_t>(Operation::Count);

struct Fixture {
  swcme::sep::Interface1D one;
  swcme::sep::Interface3D three;
  swcme1d::StepState one_state;
  swcme3d::StepState three_state;

  // Both adapters and both states are constructed before worker creation.
  // No thread prepares or mutates configuration; PST04 deliberately measures
  // the documented steady-state AMPS usage pattern only.
  Fixture(const swcme1d::Params& p1, const swcme3d::Params& p3)
      : one(p1), three(p3), one_state(one.prepare(8.0*3600.0)),
        three_state(three.prepare(8.0*3600.0)) {}
};

// Execute one logical operation using the same immutable models and prepared
// records for serial and concurrent runs.  Fixed inputs span the scalar AMPS
// adapter, direct batch fields, shock/source diagnostics, and connectivity.
Snapshot evaluate_operation(const Fixture& fixture, Operation operation) {
  constexpr double AU=swcme::constants::AU_M;
  Snapshot out;

  switch (operation) {
    case Operation::OneScalarBackground: {
      swcme::sep::BackgroundState background;
      const swcme::ModelStatus status=fixture.one.evaluate_background(
          fixture.one_state,0.83*AU,background);
      require(status.ok(),"PST04 1-D scalar background fixture failed");
      append_status(out,status);
      append_background(out,background);
      break;
    }

    case Operation::OneBatchBackground: {
      constexpr std::size_t N=7;
      const double radius[N]={0.18*AU,0.31*AU,0.47*AU,0.66*AU,
                              0.91*AU,1.20*AU,1.55*AU};
      double density[N]={},velocity[N]={},Br[N]={},Bphi[N]={},Bmag[N]={},
             divergence[N]={};
      const swcme::ModelStatus status=
          fixture.one.model().evaluate_radii_with_B_div_checked(
              fixture.one_state,radius,density,velocity,Br,Bphi,Bmag,
              divergence,N);
      require(status.ok(),"PST04 1-D batch background fixture failed");
      append_status(out,status);
      for (const double* values :
           {density,velocity,Br,Bphi,Bmag,divergence})
        for (std::size_t i=0; i<N; ++i) out.add_double(values[i]);
      break;
    }

    case Operation::OneShockSource: {
      swcme::sep::SEPSourceState source;
      const swcme::ModelStatus status=
          fixture.one.source_at_shock(fixture.one_state,source);
      require(status.ok() && source.active,
              "PST04 1-D source fixture is not an active shock");
      append_status(out,status);
      append_source(out,source);
      break;
    }

    case Operation::ThreeScalarBackground: {
      const std::array<double,3> position{{0.72*AU,0.19*AU,-0.11*AU}};
      swcme::sep::BackgroundState background;
      const swcme::ModelStatus status=fixture.three.evaluate_background(
          fixture.three_state,position,background);
      require(status.ok(),"PST04 3-D scalar background fixture failed");
      append_status(out,status);
      append_background(out,background);
      break;
    }

    case Operation::ThreeBatchBackground: {
      constexpr std::size_t N=7;
      const double x[N]={0.20*AU,0.36*AU,0.51*AU,0.68*AU,
                         0.82*AU,1.07*AU,1.35*AU};
      const double y[N]={0.03*AU,-0.09*AU,0.12*AU,-0.16*AU,
                         0.21*AU,-0.25*AU,0.18*AU};
      const double z[N]={-0.02*AU,0.04*AU,0.08*AU,-0.06*AU,
                         0.13*AU,0.09*AU,-0.17*AU};
      double density[N]={},vx[N]={},vy[N]={},vz[N]={};
      double bx[N]={},by[N]={},bz[N]={},divergence[N]={};
      const swcme::ModelStatus status=
          fixture.three.model().evaluate_cartesian_with_B_div_checked(
              fixture.three_state,x,y,z,density,vx,vy,vz,bx,by,bz,
              divergence,N);
      require(status.ok(),"PST04 3-D batch background fixture failed");
      append_status(out,status);
      for (const double* values :
           {density,vx,vy,vz,bx,by,bz,divergence})
        for (std::size_t i=0; i<N; ++i) out.add_double(values[i]);
      break;
    }

    case Operation::ThreeDirectionalShock: {
      const double direction[3]={0.91,0.37,-0.18};
      swcme3d::LocalShockState shock;
      const swcme::ModelStatus status=
          fixture.three.model().shock_state_direction_checked(
              fixture.three_state,direction,shock);
      require(status.ok() && shock.surface_exists,
              "PST04 directional shock fixture failed");
      append_status(out,status);
      append_shock(out,shock);
      break;
    }

    case Operation::ThreeDirectionalSource: {
      const std::array<double,3> direction{{0.91,0.37,-0.18}};
      swcme::sep::SEPSourceState source;
      const swcme::ModelStatus status=fixture.three.source_at_direction(
          fixture.three_state,direction,source,17,2.75e18);
      require(status.ok() && source.active,
              "PST04 directional source fixture is not active");
      append_status(out,status);
      append_source(out,source);
      break;
    }

    case Operation::ThreeConnectivity: {
      const double observer[3]={AU,0.0,0.0};
      swcme3d::ConnectivityOptions options;
      // A modest explicit lower bound keeps PST04 fast; the production solver
      // still raises it automatically if Parker winding requires more samples.
      options.scan_intervals=128;
      const swcme3d::ConnectivityState connectivity=
          fixture.three.model().observer_connectivity(
              fixture.three_state,observer,options);
      require(connectivity.connected && !connectivity.roots.empty(),
              "PST04 connectivity fixture did not produce a cobpoint");
      append_connectivity(out,connectivity);
      break;
    }

    case Operation::IndependentFailureDiagnostics: {
      // Two deliberately different failures run beside successful calls.  If
      // status were stored in a shared model/global buffer, context, sample
      // index, or offending value could be overwritten by another thread.
      const double radius[3]={0.75*AU,0.5*swcme::constants::SOLAR_RADIUS_M,
                              1.10*AU};
      double density[3]={101.0,102.0,103.0};
      double velocity[3]={201.0,202.0,203.0};
      const swcme::ModelStatus one_status=
          fixture.one.model().evaluate_radii_fast_checked(
              fixture.one_state,radius,density,velocity,3);
      require(one_status.code==swcme::StatusCode::OutsideModelDomain &&
                  one_status.sample_index==1 &&
                  one_status.has_offending_value,
              "PST04 1-D failure diagnostic fixture changed");
      append_status(out,one_status);
      for (double value : density) out.add_double(value);
      for (double value : velocity) out.add_double(value);

      const double x[3]={0.62*AU,std::numeric_limits<double>::quiet_NaN(),AU};
      const double y[3]={0.07*AU,0.0,0.0};
      const double z[3]={-0.03*AU,0.0,0.0};
      double n[3]={301.0,302.0,303.0};
      double vx[3]={401.0,402.0,403.0};
      double vy[3]={501.0,502.0,503.0};
      double vz[3]={601.0,602.0,603.0};
      const swcme::ModelStatus three_status=
          fixture.three.model().evaluate_cartesian_fast_checked(
              fixture.three_state,x,y,z,n,vx,vy,vz,3);
      require(three_status.code==swcme::StatusCode::NonFiniteInput &&
                  three_status.sample_index==1 &&
                  !three_status.has_offending_value,
              "PST04 3-D failure diagnostic fixture changed");
      append_status(out,three_status);
      for (const double* values : {n,vx,vy,vz})
        for (std::size_t i=0; i<3; ++i) out.add_double(values[i]);
      break;
    }

    case Operation::Count:
      throw std::logic_error("PST04 invalid operation sentinel");
  }
  return out;
}

// C++17 has no standard barrier, so this small one-shot gate releases every
// worker only after the complete pool has arrived.  It increases real overlap
// between distinct API calls without adding synchronization to production.
class StartGate {
 public:
  explicit StartGate(std::size_t participants)
      : participants_(participants) {}

  void arrive_and_wait() {
    std::unique_lock<std::mutex> lock(mutex_);
    ++arrived_;
    if (arrived_==participants_) {
      open_=true;
      condition_.notify_all();
      return;
    }
    condition_.wait(lock,[this]{ return open_; });
  }

 private:
  const std::size_t participants_;
  std::size_t arrived_=0;
  bool open_=false;
  std::mutex mutex_;
  std::condition_variable condition_;
};

enum class Schedule { ForwardInterleaved, ReverseInterleaved, OperationGrouped };

// Map a scheduled position back to a canonical result slot.  The three fixed
// orders vary which operations overlap while preserving a deterministic oracle
// index for every repetition.
std::size_t canonical_job(std::size_t position, std::size_t repetitions,
                          Schedule schedule) {
  const std::size_t total=kOperationCount*repetitions;
  if (schedule==Schedule::ForwardInterleaved) return position;
  if (schedule==Schedule::ReverseInterleaved) return total-1-position;
  const std::size_t operation=position/repetitions;
  const std::size_t repetition=position%repetitions;
  return repetition*kOperationCount+operation;
}

struct ConcurrentRun {
  std::vector<Snapshot> results;
  std::vector<std::string> worker_errors;
};

// THR01 uses an atomic work distributor rather than the deterministic strided
// ownership used by PST04.  The returned ledger records how many times each
// canonical job was executed; this makes a scheduler bug visible even when a
// duplicated calculation happens to overwrite a result with identical bits.
struct DynamicRun {
  std::vector<Snapshot> results;
  std::vector<unsigned int> visits;
  std::vector<std::string> worker_errors;
};

DynamicRun run_dynamically(const Fixture& fixture,
                           std::size_t thread_count,
                           std::size_t repetitions,
                           std::size_t chunk_size,
                           bool reverse_within_chunk) {
  const std::size_t total=kOperationCount*repetitions;
  DynamicRun run;
  run.results.resize(total);
  run.visits.resize(total,0U);
  run.worker_errors.resize(thread_count);
  std::atomic<std::size_t> next{0};
  StartGate gate(thread_count);
  std::vector<std::thread> workers;
  workers.reserve(thread_count);

  for (std::size_t thread=0; thread<thread_count; ++thread) {
    workers.emplace_back([&,thread] {
      try {
        gate.arrive_and_wait();
        for (;;) {
          // fetch_add models a dynamic OpenMP/AMPS work queue.  Each claimed
          // half-open interval is disjoint; the explicit visit ledger below
          // verifies that this remains true for partial final chunks too.
          const std::size_t begin=next.fetch_add(chunk_size);
          if (begin>=total) break;
          const std::size_t end=std::min(total,begin+chunk_size);
          for (std::size_t offset=0; offset<end-begin; ++offset) {
            const std::size_t position=reverse_within_chunk
                ? end-1-offset : begin+offset;
            const std::size_t job=canonical_job(
                position,repetitions,Schedule::OperationGrouped);
            const Operation operation=static_cast<Operation>(
                job%kOperationCount);
            run.results[job]=evaluate_operation(fixture,operation);
            // Exactly one worker owns this slot by construction, so no atomic
            // is needed for the ledger itself.  A defective distributor would
            // be reported after join rather than hidden by a data race here.
            ++run.visits[job];
            if (((position+thread)&3U)==0U) std::this_thread::yield();
          }
        }
      } catch (const std::exception& error) {
        run.worker_errors[thread]=error.what();
      } catch (...) {
        run.worker_errors[thread]="unknown THR01 worker exception";
      }
    });
  }
  for (std::thread& worker : workers) worker.join();
  return run;
}

// Fold serialized records in canonical job order.  Scheduler completion order
// is deliberately excluded: consumers receive an index-addressed result set,
// so a reproducible reduction must use that stable order rather than racing
// floating-point additions in worker completion order.
std::uint64_t canonical_record_hash(const std::vector<Snapshot>& records) {
  std::uint64_t hash=1469598103934665603ULL;
  for (const Snapshot& record : records) {
    for (std::uint64_t word : record.words) {
      hash^=word;
      hash*=1099511628211ULL;
    }
  }
  return hash;
}

long double canonical_numeric_reduction(const std::vector<Snapshot>& records) {
  long double sum=0.0L;
  for (const Snapshot& record : records)
    for (std::uint64_t word : record.words)
      sum+=static_cast<long double>(word&0xffffU);
  return sum;
}

ConcurrentRun run_concurrently(const Fixture& fixture,
                               std::size_t thread_count,
                               std::size_t repetitions,
                               Schedule schedule) {
  const std::size_t total=kOperationCount*repetitions;
  ConcurrentRun run;
  run.results.resize(total);
  run.worker_errors.resize(thread_count);
  StartGate gate(thread_count);
  std::vector<std::thread> workers;
  workers.reserve(thread_count);

  for (std::size_t thread=0; thread<thread_count; ++thread) {
    workers.emplace_back([&,thread] {
      try {
        gate.arrive_and_wait();
        // Strided ownership assigns every canonical output slot to exactly one
        // worker.  Threads therefore share only const model/state inputs; the
        // result harness itself introduces no write/write race.
        for (std::size_t position=thread; position<total;
             position+=thread_count) {
          const std::size_t job=
              canonical_job(position,repetitions,schedule);
          const Operation operation=static_cast<Operation>(
              job%kOperationCount);
          run.results[job]=evaluate_operation(fixture,operation);
          // Yield at a deterministic subset of positions to make different
          // calls overlap at additional internal points on common schedulers.
          if ((position+thread)%3==0) std::this_thread::yield();
        }
      } catch (const std::exception& error) {
        run.worker_errors[thread]=error.what();
      } catch (...) {
        run.worker_errors[thread]="unknown PST04 worker exception";
      }
    });
  }

  for (std::thread& worker : workers) worker.join();
  return run;
}

const char* schedule_name(Schedule schedule) {
  switch (schedule) {
    case Schedule::ForwardInterleaved: return "forward-interleaved";
    case Schedule::ReverseInterleaved: return "reverse-interleaved";
    case Schedule::OperationGrouped: return "operation-grouped";
  }
  return "unknown";
}

}  // namespace

// PST04: one immutable 1-D state and one immutable 3-D state are reused by
// 1, 2, 4, and 8 worker threads.  Every concurrent result is compared bit for
// bit with a serial oracle across three deterministic scheduling orders.
void test_pst04(swcme_test::Context& context) {
  std::cout << "PST04 concurrent prepared-state evaluation\n";

  swcme1d::Params p1;
  p1.kinematics_mode=swcme::kinematics::Mode::Ballistic;
  p1.r0_Rs=20.0;
  p1.V0_sh_kms=1450.0;
  p1.V_sw_kms=400.0;
  p1.region_mode=swcme::regions::Mode::ShockOnly;
  p1.shock_acceleration_mode=swcme::acceleration::Mode::Source;

  swcme3d::Params p3;
  p3.shape=swcme3d::ShockShape::Sphere;
  p3.kinematics_mode=p1.kinematics_mode;
  p3.r0_Rs=p1.r0_Rs;
  p3.V0_sh_kms=p1.V0_sh_kms;
  p3.V_sw_kms=p1.V_sw_kms;
  p3.region_mode=p1.region_mode;
  p3.shock_acceleration_mode=p1.shock_acceleration_mode;
  p3.cme_dir[0]=1.0; p3.cme_dir[1]=0.0; p3.cme_dir[2]=0.0;

  const Fixture fixture(p1,p3);
  constexpr std::size_t repetitions=8;
  std::array<Snapshot,kOperationCount> oracle;
  bool oracle_ready=true;
  std::string oracle_error;
  try {
    for (std::size_t operation=0; operation<kOperationCount; ++operation)
      oracle[operation]=evaluate_operation(
          fixture,static_cast<Operation>(operation));
  } catch (const std::exception& error) {
    oracle_ready=false;
    oracle_error=error.what();
  }
  context.expect_true(oracle_ready,
                      "serial PST04 oracle is physically valid: "+oracle_error);
  if (!oracle_ready) return;

  const std::array<std::size_t,4> thread_counts{{1,2,4,8}};
  const std::array<Schedule,3> schedules{{
      Schedule::ForwardInterleaved,
      Schedule::ReverseInterleaved,
      Schedule::OperationGrouped}};

  for (std::size_t thread_count : thread_counts) {
    for (Schedule schedule : schedules) {
      const ConcurrentRun run=run_concurrently(
          fixture,thread_count,repetitions,schedule);
      bool no_worker_error=true;
      for (const std::string& error : run.worker_errors)
        if (!error.empty()) no_worker_error=false;

      const std::string label=std::to_string(thread_count)+" thread(s), "+
                              schedule_name(schedule);
      context.expect_true(no_worker_error,label+" completes without exception");

      bool exact=true;
      for (std::size_t job=0; job<run.results.size(); ++job) {
        const std::size_t operation=job%kOperationCount;
        if (run.results[job].words!=oracle[operation].words) {
          exact=false;
          break;
        }
      }
      context.expect_true(
          exact,label+" is bitwise identical to the serial oracle");
    }
  }

  // A final seal check confirms the stress test itself never modified either
  // shared prepared record while exercising public compatibility fields.
  context.expect_true(
      swcme1d::prepared_state_integrity(fixture.one_state)==
          fixture.one_state.integrity_digest(),
      "concurrent reuse preserves the 1-D prepared-state seal");
  context.expect_true(
      swcme3d::prepared_state_integrity(fixture.three_state)==
          fixture.three_state.integrity_digest(),
      "concurrent reuse preserves the 3-D prepared-state seal");
}

// THR01: exercise the same public prepared-state/background/shock/source/
// connectivity records through a dynamic scheduler.  PST04 proves concurrent
// API safety for fixed ownership; THR01 additionally qualifies the scheduler
// boundary, exact work accounting, and order-independent campaign reduction.
void test_thr01(swcme_test::Context& context) {
  std::cout << "THR01 thread and scheduler reproducibility\n";

  swcme1d::Params p1;
  p1.kinematics_mode=swcme::kinematics::Mode::Ballistic;
  p1.r0_Rs=20.0;
  p1.V0_sh_kms=1450.0;
  p1.V_sw_kms=400.0;
  p1.region_mode=swcme::regions::Mode::ShockOnly;
  p1.shock_acceleration_mode=swcme::acceleration::Mode::Source;

  swcme3d::Params p3;
  p3.shape=swcme3d::ShockShape::Sphere;
  p3.kinematics_mode=p1.kinematics_mode;
  p3.r0_Rs=p1.r0_Rs;
  p3.V0_sh_kms=p1.V0_sh_kms;
  p3.V_sw_kms=p1.V_sw_kms;
  p3.region_mode=p1.region_mode;
  p3.shock_acceleration_mode=p1.shock_acceleration_mode;
  p3.cme_dir[0]=1.0; p3.cme_dir[1]=0.0; p3.cme_dir[2]=0.0;

  const Fixture fixture(p1,p3);
  constexpr std::size_t repetitions=11;
  std::vector<Snapshot> oracle(kOperationCount*repetitions);
  for (std::size_t job=0; job<oracle.size(); ++job)
    oracle[job]=evaluate_operation(
        fixture,static_cast<Operation>(job%kOperationCount));
  const std::uint64_t oracle_hash=canonical_record_hash(oracle);
  const long double oracle_sum=canonical_numeric_reduction(oracle);

  const std::array<std::size_t,4> thread_counts{{1,2,4,8}};
  const std::array<std::size_t,4> chunk_sizes{{1,3,7,19}};
  for (std::size_t repeat=0; repeat<3; ++repeat) {
    for (std::size_t thread_count : thread_counts) {
      for (std::size_t chunk_size : chunk_sizes) {
        const bool reverse=((repeat+thread_count+chunk_size)&1U)!=0U;
        const DynamicRun run=run_dynamically(
            fixture,thread_count,repetitions,chunk_size,reverse);
        const std::string label="repeat="+std::to_string(repeat)+
            ", threads="+std::to_string(thread_count)+
            ", chunk="+std::to_string(chunk_size);

        bool workers_ok=true;
        for (const std::string& error : run.worker_errors)
          workers_ok=workers_ok && error.empty();
        context.expect_true(workers_ok,label+" completes without worker error");

        bool exactly_once=true;
        for (unsigned int visits : run.visits)
          exactly_once=exactly_once && visits==1U;
        context.expect_true(exactly_once,label+" visits every job exactly once");
        context.expect_true(run.results==oracle,
                            label+" preserves every serialized result bit");
        context.expect_true(canonical_record_hash(run.results)==oracle_hash,
                            label+" preserves canonical record hash");
        context.expect_true(canonical_numeric_reduction(run.results)==oracle_sum,
                            label+" preserves canonical ordered reduction");
      }
    }
  }

  context.expect_true(
      swcme1d::prepared_state_integrity(fixture.one_state)==
          fixture.one_state.integrity_digest(),
      "THR01 preserves the shared 1-D prepared-state seal");
  context.expect_true(
      swcme3d::prepared_state_integrity(fixture.three_state)==
          fixture.three_state.integrity_digest(),
      "THR01 preserves the shared 3-D prepared-state seal");

  std::cout << "  jobs=" << oracle.size()
            << " schedules=" << 3*thread_counts.size()*chunk_sizes.size()
            << " canonical_hash=" << oracle_hash << '\n';
}
