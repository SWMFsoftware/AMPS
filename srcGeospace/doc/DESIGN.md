# srcGeospace: design of the batched trajectory engine

This document describes the code in `srcGeospace/gpu/` and how the application drivers
use it.  It reflects the sources as of 2026-10-26 (section 13 lists the additions of
2026-10-22/23); items that are still being written are marked **in progress**.  Source names are relative to `srcGeospace/`.  "Legacy" denotes the original per-trajectory tracer inherited from srcEarth (selected with `AMPS_EARTH_ENGINE=LEGACY`), the bitwise reference of the double build.

## 0. Status at a glance

| part | status |
|---|---|
| Engine core: movers, time-step selectors, boundary events, Shue boundary, exit kinematics, trap detector, dipole and Mode3D mesh fields (`gpu/Engine*.h`) | implemented for both legacy variants (gridless and Mode3D) |
| Device pool `gpu/EngineDevice.cpp`, CPU build (persistent host threads or OpenMP, `gpu/HostThreadPool.*`) | linked into the default CPU `amps` |
| Device pool built by nvc++ `-acc=gpu` / `-acc=multicore` | runs on RTX A6000 (2026-10-22): double results bitwise identical to the CPU build (`-gpu=nofma`), float results identical to the CPU float build; `fmod` replaced by `efmod` (`BUILD.md` §5.1).  The merged build (2026-10-23) matches the CPU on every finished test except C3 3d MESH (`../GPU_STATUS.md` §3.4); GPU work parked 2026-10-27 |
| Host pool for Geopack/Tsyganenko fields (`gpu/EngineHost.*`) | implemented |
| Asynchronous device-pool launches (`eg_dev_launch/wait`, `AMPS_EARTH_ASYNC`, default on) | used on CPU and GPU |
| Performance additions of 2026-10-22/23 (§13): batched request upload, x86-64 fiber switch and persistent fiber pool, speculative cutoff searches (`gpu/Speculation.*`), interior trilinear fast path and stencil-cell cache of the mesh interpolation | integrated; outputs identical to the sequential/uncached code (CPU; GPU as of 2026-10-22, the merged GPU build differs on C3 3d MESH, `../GPU_STATUS.md` §3.4); Mode3D MESH tracing 2.5x (CPU) and 2.8-3.0x (GPU) faster from the cache alone, C19 70 -> 35 min (CPU) |
| Fibers, batch front end, driver helpers (`gpu/Fiber.*`, `gpu/BatchTracer.*`, `gpu/DriverSupport.*`, `gpu/TaskStream.h`) | implemented |
| Gridless cutoff driver (`gridless/CutoffRigidityGridless.cpp`) | integrated |
| Gridless density/spectrum drivers (`gridless/DensityGridless.cpp`) | integrated |
| Mode3D cutoff glue (`3d/CutoffRigidityMode3D.cpp`: `TraceTrajectory3DAttempt_`, batch scheduling section 14.9b; nested fiber loops in PENUMBRA_SCAN and the companion rigidity list) | integrated; bitwise identical to the legacy tracer (verify: 0 mismatches) |
| Mode3D mesh flattening and upload (`gpu/Mode3DMeshExport.*`, `GlobalMagneticField::GetCompactMagneticFieldView`) | integrated; self-check bitwise (double) on 26.4M cells; CPU double builds alias the compact B array (no copy) |
| Mode3D density (`3d/DensityMode3D.cpp`); `main.cpp` (no tracing loop, unchanged) | integrated; outputs byte-identical to the legacy path |
| Forward-tracing modes (`3d_forward*`) | not ported |
| Float build (`-DGEOSPACE_PRECISION_FLOAT`) | runs every CPU test deck, gridless and Mode3D; mesh field matches double to ~6e-7 relative; penumbral cutoffs move by up to ~24 %, about 2-3 times the scatter of a rounding-level change of the double build (`PRECISION.md` §5.6, `testroot/RESULTS.md`) |

## 1. Goals and constraints

1. **Same physics, same numbers.**  In the default double/SI build every trajectory
   must reproduce the legacy tracer bit for bit (same termination, step count, time, path,
   exit state, trap diagnostics), so all existing tests and reference outputs
   remain valid.
2. **No rewrite of the search algorithms.**  The cutoff searches (UPPER_SCAN, BINARY,
   PENUMBRA_SCAN, adaptive DIRECT_ACCESS refinement), the Mode3D bracket shrinking, the
   density/spectrum direction loops, retry policies, MPI scheduling and output writers
   stay as they are.
3. **Throughput on a GPU** needs 1e4-1e5 trajectories in flight, each with its own
   adaptive time step and very different lengths (from a few steps to the 300000-step
   limit).
4. **One source** for CPU (g++, host threads, linked into the normal `amps`) and GPU (nvc++
   + OpenACC) builds.

## 2. Architecture

```
 application driver code (unchanged search logic)          one MPI rank, main OS thread
 ───────────────────────────────────────────────
  worker loop ──> RunDriverFibers / TaskStream / FiberLoop      gpu/DriverSupport, TaskStream.h
                     │  N fibers, each runs the original per-task code
                     ▼
  TraceTrajectoryAttempt_()  ──(AMPS_EARTH_ENGINE=LEGACY)──> original TraceTrajectoryImpl()
                     │ EgRequest + context index + HostFieldRef
                     ▼
  BatchTracer::Trace()   in a fiber: queue request, Fiber::Block()   gpu/BatchTracer
                         outside fibers: TraceSync() on the calling thread
                     │
       (all fibers blocked) ──> idle hook ──> BatchTracer::Progress()
                     │                                   │
          fieldKind NONE/DIPOLE/MESH            fieldKind HOST (Geopack/Tsyganenko)
                     ▼                                   ▼
  device pool (C API eg_dev_*)                 HostPool (CPU threads, g++)  gpu/EngineHost
  gpu/EngineDevice.cpp:                          RunSlot<R,U,HostField>
   OpenACC kernel (nvc++) or                     field evaluated by the legacy
   CPU threads / OpenMP loop (g++)               cFieldEvaluator through a fn pointer
   RunSlot<R,U,DeviceField>
                     │ finished EgResult (+ drift-turn states)
                     ▼
  Deliver(tag) ──> Fiber::Wake() ──> fiber resumes in the driver code with the result
```

| file | role |
|---|---|
| `EngineConfig.h` | portability macros, `Real`, unit systems (`UnitsSI`, `UnitsGeo`), exact `emin/emax` |
| `EngineTypes.h` | POD types crossing the host/engine boundary: `EgRequest`, `EgContext`, `EgResult`, `EgTrapConfig`, `EgShue`, `EgTurnState`, enums |
| `EngineMath.h` | `Vec3<R>` helpers with legacy operand order |
| `EngineMovers.h` | port of `gridless/GridlessParticleMovers.cpp` |
| `EngineField.h` | `DeviceField<R,U>`: NONE, analytic dipole, Mode3D mesh |
| `EngineMesh.h` | flattened AMR tree + port of the AMPS row-stencil interpolation |
| `EngineTrap.h` | fixed-capacity port of `Geospace::TrajectoryTrap::Detector` |
| `EngineStep.h` | resumable tracer `Tracer<R,U>::Begin/Advance`, boundary events, Shue, exit kinematics |
| `EnginePool.h` | `RunSlot()` (one slot, up to K iterations), `EgMeshHost`, arena sizing |
| `EngineDeviceAPI.h` / `EngineDevice.cpp` | C API and implementation of the device pool |
| `EngineHost.h/.cpp` | `HostPool`, `HostMesh`, `TraceSync()` |
| `HostThreadPool.h/.cpp` | `HostParallelFor()`: persistent worker threads that launch the CPU kernels (default), `HostRuntimeIsOpenMP()` |
| `BatchTracer.h/.cpp` | request queue, context table, pool driving, routing, statistics |
| `Fiber.h/.cpp` | cooperative fibers and scheduler |
| `DriverSupport.h/.cpp`, `TaskStream.h` | helpers that turn driver worker loops into fibers |
| `EngineGlue.h/.cpp` | legacy <-> engine conversions, `AMPS_EARTH_ENGINE_VERIFY` comparison, warnings |

## 3. Fibers: sequential drivers as request producers

**Why.**  Every application driver calls a blocking "trace one trajectory" function from
deeply nested sequential code (bisection loops, retries, adaptive refinement).  Turning
each of them into an explicit state machine would duplicate thousands of lines of
validated logic.  Instead each unit of driver work runs in a fiber.  A trajectory
request suspends the fiber, and when no fiber can run, the scheduler advances the whole
trajectory pool once and resumes the fibers whose trajectories finished.  Each search
sees exactly the sequence of results it would see serially, so search decisions and
outputs are unchanged.

**Implementation (`Fiber.cpp`).**
* Context switch: on x86-64 a small assembly routine (`geospace_fiber_switch`: pushes
  the callee-saved registers, saves MXCSR and the x87 control word, swaps the stack
  pointer); other architectures, or `-DGEOSPACE_FIBER_UCONTEXT`, use `ucontext_t`
  (`makecontext`/`swapcontext`, which also saves the signal mask with a system call).
  One scheduler per OS thread (`thread_local Scheduler*`), normally the MPI rank's main
  thread, so MPI calls made in fiber code come from the main thread.
* Stacks: `mmap(MAP_PRIVATE|MAP_ANONYMOUS|MAP_NORESERVE)` of
  `AMPS_EARTH_FIBER_STACK_KB` (default 1024 KB, minimum 64) plus one `PROT_NONE` guard
  page at the low end, so an overflow faults instead of corrupting memory.  Finished
  fibers are recycled through a free list, and the stacks stay in a persistent
  per-process pool that later loops reuse (no unmapping and no page faults between
  loops).
* `ParallelFor(begin,end,maxWorkers,body)`: workers share an atomic index counter.
  Called outside fibers, it creates the scheduler and runs to completion.  Called
  inside a fiber (nested loops: energies, then directions), the calling fiber itself
  executes iterations and, before each iteration, spawns further child fibers whenever
  budget `MaxFibers()-live` is available (up to `maxWorkers` in total).  A nested loop
  that starts while the budget is exhausted therefore fans out as soon as other loops
  finish, instead of running all its iterations one at a time.  The parent then waits
  for its children; a finishing child wakes the parent only when the parent is blocked
  in that join (not while it is blocked on one of its own trajectories).
* `InFiber()` is false inside an OpenMP parallel region.  Trajectory requests from
  OpenMP threads therefore fall back to synchronous tracing.
* The idle hook is `BatchTracer::Progress()`.  If it reports no progress while fibers
  are blocked, the scheduler throws a deadlock error.  An optional service hook runs in
  scheduler context between pool advances (the drivers use it for MPI progress
  publication).
* The first exception that escapes a fiber body stops the group (`stop` flag) and is
  rethrown by `ParallelFor`.

**Driver helpers.**
* `DriverSupport.cpp`: `DriverBatchEnabled()` (= `AMPS_EARTH_ENGINE` is not
  `LEGACY`), `DriverFiberCount()`, `RunDriverFibers(nWorkers,hostThreads,worker,service,
  label)` (sets the host-executor thread count and the service hook, runs the workers,
  and prints fiber and pool statistics when `AMPS_EARTH_ENGINE_VERBOSE=1`), and
  `FiberLoop()` for nested loops (plain loop outside fibers).
* `TaskStream.h`: one task stream per rank that yields exactly the global task ids the
  original scheduler would give this rank: DYNAMIC chunks fetched on demand through
  `DynamicMpiLocationScheduler`, BLOCK_CYCLIC and STATIC index maps.
  `RunStreamBatch()` computes the tasks in fibers, then sorts the results by task id and
  accumulates them in that order, so order-dependent floating-point sums are unchanged.

**Integration examples.**
* Gridless cutoff (`CutoffRigidityGridless.cpp`): the per-rank task loop becomes
  `RunDriverFibers` over an inline task stream.  Each fiber runs the unchanged
  `ProcessTask()`, and results are sorted by task id before accumulation.
  `TraceTrajectoryAttempt_()` routes each attempt to `EngineTraceTrajectory_()`
  (engine) or `TraceTrajectoryImpl()` (legacy).  The per-thread field-evaluator cache
  (`GetCachedFieldEvaluator`) returns `shared_ptr`s, so a fiber keeps its snapshot alive
  even if a later request evicts it.  Outside fibers it holds one entry, as the original
  cache did.  Inside fibers it holds an LRU set of 64 snapshots.  Its key holds the same
  values as the original `%.17g` string key, but is compared field by field without
  allocation.
* Per-trajectory host overhead matters because every request of a rank passes through
  the rank's main thread.  The main session cached the snapshot fingerprint, the engine
  field kind and the context index per evaluator/thread, and skipped the per-trajectory
  `AmpsParam` copy when it would change nothing the tracer reads.  Callgrind then
  measured about 4200 serial host instructions per trajectory on the rank main thread.
* Gridless density (`DensityGridless.cpp`): energies run as fibers (`RunDriverFibers`),
  each fans out its direction block with `FiberLoop`, and the directions are recorded in
  the original order.  POINTS/SHELL task loops use `TaskStream` + `RunStreamBatch`.
  With legacy (non-reentrant) Tsyganenko wrappers and several epochs, the host executor
  is limited to one thread.

## 4. Batch front end (`BatchTracer`)

* Singleton.  `Enabled()` reads `AMPS_EARTH_ENGINE` once.
* **Contexts.**  Run-level data shared by many requests (variant, stopping surfaces,
  trap configuration, field kind, dipole moment, Shue parameters) live in an
  `EgContext` table.  `RegisterContext()` deduplicates by `memcmp` under a mutex.
  Contexts are uploaded to a pool only while that pool is idle
  (`eg_dev_set_contexts`, `HostPool::SetContexts`), because the trap arena stride
  depends on the largest configuration.
* **Requests.**  `Trace()` in a fiber assigns a unique tag, pushes the request onto a
  FIFO queue and blocks.  Outside a fiber it calls `TraceSync()`: the same engine code
  on the calling thread, with K = 2^30 iterations.
* **`Progress()`** (the fiber idle hook):
  1. upload new contexts to idle pools;
  2. submit queued requests in FIFO order.  A request whose context is not yet resident
     blocks further submissions to that pool until the pool drains;
  3. device pool, asynchronous mode (default, `AMPS_EARTH_ASYNC` not `0`): if a launch
     is pending, `eg_dev_wait` completes it and `eg_dev_launch` immediately starts the
     next one, over the still-active slots plus the requests submitted meanwhile.  Only
     then are the finished trajectories collected and `Deliver()`ed.  The resumed fibers
     run their search logic and submit new requests while the next launch integrates.
     If nothing is pending but uncollected active slots exist, a launch is started.
     Synchronous mode (`AMPS_EARTH_ASYNC=0`): one `eg_dev_advance` (launch and wait),
     then collect and deliver.  Host pool, if active: one `HostPool::Advance()`
     (synchronous), then collect and deliver;
  4. if nothing progressed but requests wait for a context, upload it now.
* **Routing.**  `fieldKind != kFieldHost` (NONE, DIPOLE, Mode3D MESH) goes to the device
  pool; `kFieldHost` (Geopack IGRF, T96/T01/T05/TA15/TA16 Fortran models, and the
  Mode3D forced-analytic empirical path) goes to the host pool.  Both pools advance in
  the same `Progress()` call.
* **Defaults** (environment overrides in `BUILD.md`).  Device pool capacity and K:
  65536 and 2048 on an accelerator, otherwise max(256, 64 x host threads) and 256.
  Host pool: max(64, 16 x host threads) slots, K = 64.  Drift-turn capacity 512.
* `PrintStatistics()` reports batched/synchronous requests, launches, the mean number
  of active slots per launch, kernel+collect seconds and the backend name.

## 5. Persistent trajectory pool and K-step kernels

**Resumable tracer (`EngineStep.h`).**  The legacy tracer integrates one trajectory inside one
call (`TraceTrajectoryImpl` gridless, `TraceTrajectory3D` Mode3D).  The engine splits
the same algorithm into `Tracer::Begin()` (everything before the integration loop:
result initialization, momentum from rigidity, trap-detector reset and initial sample)
and `Tracer::Advance(K)` (at most K iterations of the loop body including the loop exit
tests).  All loop state lives in a slot:

```
template <class R> struct EgSlot {
  Vec3<R> x,p;       // state (Real)
  double t, s;       // trace time and path length, internal units (always double)
  int n, done;       // step count; done: -1 new, 0 running, 1 finished
  int lastNode;      // Mode3D tree-search hint
  EgHybrid<R> hyb;   // HYBRID mover context (thread_local in the legacy tracer)
  EgTrapState trap;  // trap detector scalars (profiles live in the arena)
  EgResult res;
};
```

`RunSlot()` (`EnginePool.h`) builds the trap detector on the slot's arena, calls
`Begin` for a new slot and then `Advance(K)`, stores the mesh hint, merges dipole-error
statistics and converts fatal field errors into a finished result with `res.error` set.

**Per-trajectory adaptive dt.**  Every iteration selects its own dt exactly as the legacy tracer
does: `SelectLegacyCutoffTraceDt` (gridless, Re-based distances),
`SelectLegacyCutoffTraceDt3D` (Mode3D, metre-based), or `SelectStructuredTraceDt[3D]`
(gyro-angle limit 0.15 rad, upper bound only).  No global time step exists; a fixed one
would need 86-120x more steps on a penumbral workload (`PRECISION.md` §5.5).

**Device pool (`EngineDevice.cpp`).**
* `capacity` slots with request copies, per-slot trap arenas, an active list, a free
  list (host side) and a finished list.
* `eg_dev_submit` takes free slots, sets `done=-1`, `lastNode=-1` and copies that one
  slot and its request to the device.
* `eg_dev_advance` (synchronous) uploads the active list, launches the kernel once over
  all active slots, downloads one `done` flag per active slot and moves finished slots
  to the finished list.
* Asynchronous launches (added by the main session on 2026-10-21):
  `eg_dev_launch` moves the active slots to an `inflight` list, uploads it and starts
  the kernel without waiting.  It returns the number of slots, or -1 if a launch is
  already pending.  `eg_dev_wait` completes the launch, downloads the done flags and
  sorts `inflight` into finished and still-active slots.  `eg_dev_pending_launch`
  reports whether a launch is pending.  Slots submitted while a launch runs wait for
  the next launch; contexts and the mesh cannot change while a launch is pending.
  OpenACC: `async(1)` on the compute region, `acc_wait(1)` in `WaitKernel()`.  CPU
  build: `HostParallelForAsync` runs the body on `threads` worker threads of a separate
  pool instance (the calling thread does not participate) and `HostParallelWait()`
  completes it.  With `AMPS_EARTH_HOST_RUNTIME=OPENMP` the launch is synchronous and the
  wait is a no-op.
* `LaunchKernel` (formerly `AdvanceKernel`) distributes `AdvanceOne(a,...)`, the complete kernel body, over the
  active slots a = 0..n-1.  `AdvanceOne` builds an `EgMeshView`, two row stencils and a
  `DeviceField` on the stack and calls `RunSlot`.  Launchers:
  * OpenACC build: `#pragma acc parallel loop gang vector vector_length(64)` with
    `deviceptr(...)` for every buffer and `firstprivate(ms)` for the mesh scalars;
  * CPU build, default: `HostParallelFor(n,threads,body)` (see "CPU kernel launcher"
    below);
  * CPU build with `AMPS_EARTH_HOST_RUNTIME=OPENMP`: `#pragma omp parallel for
    schedule(dynamic,1) num_threads(threads)`.
* `eg_dev_collect` uploads the indices of the finished slots, `GatherKernel` packs their
  `EgResult`s (and drift-turn states) contiguously, one download returns them, and the
  slots go back to the free list.  `GatherKernel` is an `acc parallel loop` in the
  OpenACC build and a serial loop in the CPU build.
* Long and short trajectories share the hardware: a slot that finishes is refilled at
  the next `Progress()`, and K bounds the duration of one launch.

**Host pool (`EngineHost.cpp`).**  Same `RunSlot` with a `HostField` functor that calls
the legacy evaluator through `HostFieldRef{fn,obj}`.  Slots are advanced by
`HostParallelFor`, or by `omp parallel for schedule(dynamic,1)` with
`AMPS_EARTH_HOST_RUNTIME=OPENMP`.  The host pool is always g++ code, also in GPU
builds.  Exceptions from Geopack/Tsyganenko are caught per slot and returned as
`kErrHostField` with the message, which the gridless glue rethrows as
`std::runtime_error` in the requesting fiber.

**CPU kernel launcher (`HostThreadPool.cpp`).**  The engine alternates short parallel
launches (advance all active slots by K iterations) with a serial phase on the rank's
main thread (resume fibers, run the search logic, submit requests).  With libgomp's
default `OMP_WAIT_POLICY`, idle OpenMP workers spin between parallel regions.  The
tests run, for example, 4 MPI ranks x 16 threads on 36 hardware threads, so the
spinning workers of all ranks starve every rank's serial phase.  The main session
measured gridless test F1 at 28 s with the default policy, 3.1 s with
`OMP_WAIT_POLICY=PASSIVE` and 2.8 s with the thread pool, against 2.9-4.1 s for the
the legacy srcEarth executable, with bitwise identical outputs.  libgomp reads the policy when it
is loaded, and the test runners do not set it.  `HostParallelFor`
therefore provides "parallel for, dynamic schedule of single iterations" with:
* persistent `std::thread` workers created on first use (`nThreads-1` helpers; the
  calling thread participates);
* an atomic iteration counter;
* workers that spin for `AMPS_EARTH_SPIN_US` microseconds (default 50) after a launch,
  for back-to-back launches, then block on a condition variable;
* one launch at a time per process (not reentrant), which holds because the engine
  launches only from the rank's main thread.

The kernel body is the same code for every launcher, so the results do not depend on
the launcher.

**Memory per slot** (g++ x86-64; the nvc++ layout of these PODs is expected to match):

| item | double build | float build |
|---|---|---|
| `EgSlot<R>` | 816 B | 768 B |
| `EgRequest` | 136 B | 136 B |
| `EgResult` (inside the slot) | 368 B | 368 B |
| `EgContext` (per context, not per slot) | 440 B | 440 B |
| `EgRowStencil<R>` (2 per kernel thread, on the stack) | 1928 B | 1604 B |
| trap arena, drift detection on, defaults (24 bins, window 1, lag 1) | 5496 B | 5496 B |
| trap arena, 72 bins, window 2, lag 4 | 38372 B | 38372 B |
| drift-turn states (only with `AMPS_TRAP_TURN_DIAGNOSTICS=1`) | 136 B x 512 | same |

A 65536-slot GPU pool with drift detection at the defaults needs about 65536 x (816 +
136 + 5496) B = 0.42 GB, plus 65536 x 368 B for the packed results.

**Kernel resources** (nvc++ 24.5 -O2 `-gpu=cc86,nordc,nofma,ptxinfo`, the
`build_gpu_engine.sh` default): `AdvanceKernel` uses 255 registers and a 6704-byte stack
frame per thread in double (6720 bytes with FMA), 168 registers and 5296 bytes in float;
`GatherKernel` uses 28 registers.  At `vector_length(64)` and 255
registers, an SM holds at most 4 gangs (256 threads).

## 6. Exact-port strategy

The engine is a port, not a reimplementation.  The rules that make the double/SI
build bitwise identical to the legacy tracer:

* **Literal expressions.**  Every ported function keeps the legacy operations in the
  same order: `vadd/vsub/vmul/vdot/vcross/vnorm` mirror the helpers in
  `GridlessParticleMovers.h`; `emin(a,b)=(b<a)?b:a` and `emax(a,b)=(a<b)?b:a` reproduce
  `std::min/std::max` (including NaN behaviour); `std::min({...})` is a left fold;
  Butcher tableaus keep every zero coefficient in the stage sums.
* **Ported code** (file: source):
  * `EngineMovers.h`: `BorisStep`, `HigueraCaryMagneticMidpointStep`, `HC4StepUnchecked`
    (applied as 4 checked substeps per step), `RKStepGeneric<2,4,7>`, `GCStepGeneric<2,4,7>`,
    `EvaluateFieldGeometry`, `ComputeGuidingCenterAdiabaticity`, `GuidingCenterRhs`,
    `ReconstructGuidingCenterMomentum`, `SelectPerpendicularUnitVector`, the HYBRID
    eligibility/prepare/commit logic, and `StepParticleChecked`.
  * `EngineStep.h`: the loops of `TraceTrajectoryImpl` and `TraceTrajectory3D`,
    `SelectLegacyCutoffTraceDt[3D]`, `SelectStructuredTraceDt[3D]`,
    `TrajectoryBoundary::FindFirstEvent`, `PopulateExitKinematics` /
    `CompleteExitPitchAngle`, `SWMFCoupledAccess::ClassifyPosition/ClassifySegment`
    (Shue), and the Mode3D field-availability checks.
  * `EngineTrap.h`: `Geospace::TrajectoryTrap::Detector` (§9).
  * `EngineField.h`: `Dipole::GetB_Tesla` / `EvaluateConfiguredDipoleField_`,
    `cMode3DMeshFieldEval::GetB_T`, `RecordDipoleFieldError_`.
  * `EngineMesh.h`: the AMPS tree search and row-stencil interpolation (§8).
* **Variants.**  `EgContext::variant` selects gridless or Mode3D semantics where the
  two legacy tracers differ: the legacy dt formula, the path-length accumulation
  (`s += finite ? len : 0` vs `if finite s += len`), field-availability tests and the
  optional Shue boundary.
* **Semantics of legacy `exit()`.**  Conditions on which the legacy code terminates the
  process (mesh consistency, row-stencil overflow, negative weight, invalid strict
  dipole) become `EgResult::error`; the slot finishes as `NumericalFailure` and the host
  decides how to abort.
* **Unit conversions** go through `ToLength<U>()` etc., which are `if constexpr`
  identities for `UnitsSI`, so the SI build contains no extra arithmetic.
* **Same compiler and flags.**  In the CPU build the engine is compiled by the AMPS
  compiler with the AMPS flags (`mpicxx -O3`, no `-march`, so no FMA contraction on
  x86-64).

**Verification (`AMPS_EARTH_ENGINE_VERIFY=1`).**  `TraceTrajectoryAttempt_()` also runs
the original `TraceTrajectoryImpl()` on the calling thread for every engine trajectory,
and `VerifyCompare()` compares bit patterns (`memcmp`) of: termination, steps,
traceTime, traceDistance, exit valid flag, exit position and momentum (3 components
each), exit cos(alpha), exit time, mirror points, bounce cycles, drift revolutions,
drift angle, minimum outer margin, momentum spread, trap mechanism, drift matched-bin
fraction and the number of drift-turn states.  The first 20 mismatches per process are
printed, and a summary line `[earth-engine][verify] trajectories compared=N bitwise
mismatches=M` is printed at exit.  Expected outcome:

* CPU build (g++ device pool, either launcher), double: 0 mismatches by construction.
* nvc++ builds, double: `scripts/build_gpu_engine.sh` compiles with `-gpu=nofma` /
  `-Mnofma` by default, because nvc++ 24.5 contracts `a*b+c` into FMA at -O2 (4867
  vs 3363 DFMA instructions in the kernel SASS with and without `nofma`; 544 vs 0
  `vfmadd` in the multicore object).  Without FMA, `+ - * /` and `sqrt` are IEEE
  correctly rounded on both sides, so the movers, dipole and mesh interpolation are
  expected to match.  The trap detector (`atan2`, `hypot`) and the Shue boundary
  (`pow`) call math-library functions whose GPU implementations may differ from glibc
  in the last bit, so trajectories that use them can show rounding-level mismatches.
  Run on an RTX A6000 (2026-10-22, `-gpu=nofma`): 3,301,344 gridless trajectories (C1,
  F1, F11) and 17,842 Mode3D trajectories (C1 3d MESH and ANALYTIC), 0 bitwise
  mismatches; the mesh self-check matched AMPS on 58,793 samples.
* Float build: not bitwise comparable; validate statistically.

## 7. Fields and routing

| `EgFieldKind` | evaluator | where |
|---|---|---|
| `kFieldNone` | B = 0 | device pool |
| `kFieldDipole` | centered dipole with moment `(dipoleMx,My,Mz)` prepared on the host exactly like the legacy evaluator (gridless: M x m_hat from `Dipole::MakeParams`; Mode3D: moment x sin/cos tilt); `dipoleStrict` reproduces the Mode3D fatal check | device pool |
| `kFieldMesh` | Mode3D compact global field + AMPS row stencil (§8) | device pool |
| `kFieldHost` | legacy `cFieldEvaluator::GetB_T` via function pointer | host pool |

* `DeviceField::Available()` reproduces the sticky "field sample unavailable" status of
  `cMode3DMeshFieldEval`: a point inside the configured AMPS domain that falls in an AMR
  coverage hole, or an empty interpolation row, makes the sample unavailable and the
  Mode3D tracer stops with `InvalidField`.  Gridless fields are always available.
* `sampleDipoleError` (Mode3D mesh run of a DIPOLE case) accumulates |B_mesh - B_dip| /
  |B_dip| per trajectory (count, sum, max, location of max) into `EgResult` for the
  driver's diagnostic.

## 8. Mode3D mesh: flattening and the AMPS row-stencil port (`EngineMesh.h`)

The legacy Mode3D code evaluates B with `findTreeNode(x,lastNode)`, then
`PIC::InterpolationRoutines::CellCentered::Linear::InitStencil(x,node,cRowStencil)`
and `B = sum_s w_s GlobalMagneticField_[3(Temp_ID Nc + i + Nx(j + Ny k))]`.  The AMPS
tree is a pointer structure on the host, so it is flattened into index arrays
(`EgMeshHost`, SI double; the pool converts to `EgMeshView<R>` in engine units):

| array | AMPS source |
|---|---|
| `xmin[3n]`, `xmax[3n]` | node box |
| `ixMin[3n]`, `sizeIndex[n]` | integer node position/size at the finest refinement (`xMinGlobalIndex`, `NodeGeometricSizeIndex`) |
| `up[n]`, `down[8n]` | parent and children (-1 when absent) |
| `level`, `minNeibLevel`, `maxNeibLevel` | `RefinmentLevel`, `min/maxNeibRefinmentLevel` |
| `used[n]` | `IsUsedInCalculationFlag` |
| `tempId[n]` | `Temp_ID` = compact leaf index (or -1) |
| `B[3 interiorCells]`, `presence[interiorCells]` | `GlobalMagneticField_`, `GlobalCellPresence_` (an `int` vector in the legacy code; `unsigned char` here, converted by the glue) |
| scalars | root, `_BLOCK_CELLS_X/Y/Z_`, global box, `dxMaxRefinement`, `GlobalUsedLeafBlocks_`, `GlobalInteriorCellCount_` |

Ported AMPS routines: `cMeshAMRgeneric::findTreeNode` (double and integer forms),
`cTreeNodeAMR::neibNodeFace/Edge/Corner` in the integer-index form used with
`_MESH_GLOBAL_NODE_CONNECTION_INFO_MODE_=_OFF_AMR_MESH_`, `AddRowStencilCellByCoordinate`,
`AddRowStencilCell`, `BuildConstantRowStencil`, `CellCentered::Linear::InitStencil`
(AMPS routine, non-uniform AMR, non-periodic), `GetTriliniarInterpolationStencil` and
`GetTriliniarInterpolationMutiBlockStencil` (row parts), `cRowStencil::AddCell /
MultiplyScalar / Add / Normalize`, and `CheckedGlobalCellIndex_` + `InterpolateField_`
from `GlobalMagneticField.cpp`.  Details that matter for bitwise identity:

* `cRowStencil::AddCell` merges an identical (node,i,j,k) entry by adding weights,
  else appends.  The row capacity is 80 (`kMaxRowStencil`; AMPS allows 128): a
  multiblock stencil (8 support points x 8 fine cells) blended with one trilinear stencil
  has at most 72 distinct cells.  Overflow and negative weights are reported as errors,
  as AMPS would `exit()`.
* `(int)floor(v)` follows x86-64 semantics (`INT_MIN` for out-of-range or NaN,
  `FloorToIntX86`) so degenerate positions give the same result on the device.
* The tree search result does not depend on the start-node hint (the integer search
  climbs until the index box contains the point, then descends); `EgSlot::lastNode`
  keeps the hint for speed only.

**Export and upload** (`gpu/Mode3DMeshExport.*`).  `Mode3DMesh::EnsureUploaded(label)`
is called by the Mode3D drivers (`Mode3D::Run` after `BeginFrozenFieldBatch`,
`RunCutoffRigidity`, `RunDensityAndFlux`) before any fiber starts.  It flattens every
node reachable from `rootTree` (pre-order DFS; the numbering does not affect results),
reads the compact arrays through `GlobalMagneticField::GetCompactMagneticFieldView()`,
and uploads them when the snapshot generation or a topology/Temp_ID checksum changes.
The flattened arrays stay resident in the exporter.  `BatchTracer::SetMesh(m,generation)`
passes them to the device pool and to `HostMesh` (synchronous tracing); both refuse to
change the mesh while trajectories are in flight.  In CPU double/SI builds both **alias**
the exporter arrays and the compact `GlobalMagneticField_` B array instead of copying
(0.6 GB per rank saved for C1 at 0.02 Re); float builds convert to normalized float and
GPU builds copy to device memory.  A fiber never uploads: a trace with a stale mesh
inside a fiber stops with a clear message.  The self-check (`AMPS_EARTH_MESH_VERIFY=1`,
or with `AMPS_EARTH_ENGINE_VERIFY=1`) compares node search and interpolation with AMPS
(bitwise in double; C1 MESH: 139785 nodes, 26.4M cells, 0 mismatches).

**Float-only boundary handling.**  Two AMPS round-off tolerances are adapted in the
float build only (double keeps the literal legacy/AMPS values): the 1e-10 cell-index
tolerance of `AddRowStencilCellByCoordinate` (in float `6+1e-10` rounds to 6.0) is
scaled to the float resolution, and samples within a few float ulps of a domain face are
treated as outside the configured domain.  Without the latter, a sample 9 m outside the
55 Re face (float resolution there ~24 m) rounded onto the face, where the AMR search
reports "outside the mesh" while the configured-box test said "inside", which the legacy code
classifies as an AMR coverage hole (INVALID_FIELD).

## 9. Trap detector with a fixed-capacity arena (`EngineTrap.h`)

The legacy `Geospace::TrajectoryTrap::Detector` keeps drift-recurrence profiles in growing
`std::vector`s.  A kernel cannot allocate, so the state is split into scalars
(`EgTrapState`, inside the slot) and a per-slot arena whose size follows from the
configuration:

```
nb       = max(8, TRAP_DRIFT_PROFILE_BINS)
window   = max(1, TRAP_DRIFT_WINDOW_REVOLUTIONS)
maxLag   = max(max(1, TRAP_DRIFT_MAX_RECURRENCE_LAG), window)
retained = maxLag + window                       completed profiles kept by the legacy detector
doubles  = (1 + retained) * 9 * nb + retained    current + ring of profiles (9 quantities/bin) + mean radii
ints     = (1 + retained) * nb + maxLag + 1      per-bin counts + stable-comparison counters per lag
```

* The legacy `push_back` + `erase(begin())` above `retained` is a ring buffer with the
  same logical order.  `CombineProfiles()` is evaluated bin by bin with the same
  summation order, so the detector is bitwise identical in double.
* The arena stride is the maximum over all registered contexts (`TrapArenaStride`), and
  the arena is reallocated when contexts are uploaded.  It is allocated only when trap
  detection and drift detection are both enabled.
* The detector always runs in SI double, also in the float build (state converted on
  each update).
* Drift-turn diagnostics (`AMPS_TRAP_TURN_DIAGNOSTICS=1`) are stored in a fixed
  per-slot buffer of `AMPS_EARTH_TURN_CAPACITY` states (default 512).  When it is full
  the result carries `kErrTurnStateOverflow`; the gridless glue prints a one-time
  warning and the recorded list is truncated.

## 10. Units and precision

Compile-time `Real` (`double` by default, `float` with `-DGEOSPACE_PRECISION_FLOAT`).
State, field samples and mover arithmetic use `Real`.  Time, path length, the trap
detector, Shue boundary, exit kinematics and all API types stay double.  The float
build uses normalized geospace units (L = Re, v = c, m = amu, q = e, B = amu c/(e Re))
because SI float underflows |p|^2 and (m c)^2, and the AMPS PIC::Units normalization
overflows |x|^2.  Full analysis and data: `PRECISION.md`.

## 11. GPU data management

* `DevBuf<T>` (in `EngineDevice.cpp`) owns a host mirror (`std::vector<T>`) and a device
  pointer.  In OpenACC builds the device pointer comes from `acc_malloc` (freed with
  `acc_free`) and transfers are explicit `acc_memcpy_to_device/from_device` of index
  ranges.  In CPU builds the device pointer aliases the host vector, and the
  transfer calls are no-ops, so the same driver code runs without copies.
* Kernels receive raw device pointers through `deviceptr(...)`.  There are no `acc
  data` regions, no implicit copies of large arrays and no managed memory.  Scalars
  are passed by value (`firstprivate(ms)`; the compiler stages `ms` in CUDA shared
  memory).
* Per launch: upload the active list (n ints), download n done flags (after
  `acc_wait(1)` in asynchronous mode).  New requests are not copied one by one: they are
  recorded on the host and uploaded once per launch, and a scatter kernel on queue 1
  writes them and a reset slot into their slots (`FlushPending`, section 13).  Per
  collection: upload the finished indices, download packed results (and turn states).
  Contexts, trap arenas and the mesh are uploaded only while the pool is idle.
* No device code uses virtual functions, exceptions, `std::vector`, `std::string` or
  dynamic allocation.  Everything reached from a kernel is an inline or template
  function in a header, so nvc++ generates `acc routine seq` device versions
  implicitly; `EngineDevice.cpp` is the single translation unit, which allows
  `-gpu=nordc`.
* nvc++ 24.5 has no device `fmod`.  `EngineTrap.h` needs it once, and the proposed
  `efmod()` replaces it with an exact `remainder()`-based form (`BUILD.md` §5.1).
* `eg_dev_is_accelerator()` returns `acc_get_device_type() != acc_device_host`.  The
  backend name is `OpenACC/<real>`, `CPU-threads/<real>` (CPU build, default) or
  `CPU-OpenMP/<real>` (`AMPS_EARTH_HOST_RUNTIME=OPENMP`); `BatchTracer::BackendName()`
  adds the host pool's launcher (`CPU-threads/<real>` or `CPU-OpenMP/<real>`).

## 12. Limitations and future work

* **GPU speed.**  The GPU build gives the CPU results (except C3 3d MESH,
  `../GPU_STATUS.md` §3.4) but, for the test decks, is slower than the CPU build: they
  keep a few hundred trajectories in flight (e.g. 168 locations in C1), far too few for a
  GPU, and the Mode3D mesh interpolation is latency-bound (tree searches, up to 72-entry
  stencils in local memory).  The additions of §13 reduced the gap (C1 3d MESH 601 ->
  487 s at the default K = 2,048, 345 s with K = 256; C1 gridless 11.2 -> 4.5 s) but the
  CPU is still faster (80 s and 1.6 s).  A GPU pays off only for decks with many
  locations or directions.  Measurements of the full 1-GPU passes (2026-10-26/27):
  `../GPU_STATUS.md` §4.
* **Two OpenMP runtimes.**  The nvc++ OpenACC runtime (`libacchost`) requires
  `libnvomp.so`, which exports the full `omp_*`/`GOMP_*` ABI, while AMPS uses
  libgomp.  In the linked `amps`, libgomp is a direct dependency and libnvomp an
  indirect one, so AMPS's OpenMP calls bind to libgomp.  Whether both runtimes
  coexist cleanly at run time is not verified.
* **Occupancy.**  255 registers and 6.7 KB of stack per thread (double) limit an SM to
  256 resident threads.  Options: split the kernel (field evaluation vs. mover), move the
  row stencils out of the per-thread frame, `-gpu=maxregcount`, structure-of-arrays
  state.
* **Fiber budget vs. pool size.**  Every fiber has at most one trajectory in flight.
  `AMPS_EARTH_FIBERS` defaults to 4096 on CPU builds (pools of 64 x host threads) and to
  the pool capacity (65536) when the device pool runs on an accelerator.  Up to 16384
  live fibers each stack has its own mapping with a guard page (two kernel mappings per
  fiber; `vm.max_map_count` is 65530 here); above that, stacks are carved from slabs of
  256 stacks sharing one mapping and without guard pages (`AMPS_EARTH_FIBER_GUARD=1`
  forces guards).  Stacks are reserved with `MAP_NORESERVE`, so only touched pages use
  memory.
* **Context barrier.**  A new context (different stopping surfaces, trap settings or
  field) can only be uploaded to an idle pool, so it drains the pool first.
* **Host-only field models** (Geopack, Tsyganenko) always run on CPU threads.
* **Multi-GPU.**  The pool uses the default device of the process.  With several MPI
  ranks per node run `amps` through `scripts/amps_gpu_select.sh`, which sets
  `CUDA_VISIBLE_DEVICES` from the local rank (`GEOSPACE_GPUS=0,1`); there is no
  `acc_set_device_num` call in the code.
* **Float build**: validated statistically on the CPU test list (`PRECISION.md` §5.6,
  `VALIDATION_CPU_2026-10-21.md` §4).  It does not reproduce the double results bit for bit, and
  it fails C1 3d MESH (Rc vs Stormer tolerance) and C17 (one reflection pair over the
  limit), which the rounding-level baseline passes.  Absolute lengths given in metres
  need float-specific handling near the outer box (`PRECISION.md` §6.1; one case is
  still open).
* **Forward-tracing modes** are not ported.

## 13. Performance additions (2026-10-22/23)

All of these keep the results identical to the code they replace (double: bitwise;
float: identical between CPU and GPU float builds; verified on the GPU on 2026-10-22 --
the merged GPU build of 2026-10-23 differs from the CPU on C3 3d MESH,
`../GPU_STATUS.md` §3.4).

* **Batched request upload** (`EngineDevice.cpp`, `FlushPending`).  New requests are
  recorded on the host and written into their slots by one scatter kernel per launch
  (one host-to-device copy) instead of one synchronous copy per request.
* **Fiber switch and pool** (`Fiber.*`).  On x86-64 the fibers switch with a small
  assembly routine (callee-saved registers, MXCSR and x87 control word, stack pointer)
  instead of `swapcontext()`, which saved and restored the signal mask with a system
  call on every switch; stacks are kept in a persistent pool.  `-DGEOSPACE_FIBER_UCONTEXT`
  restores `swapcontext()`.
* **Speculative cutoff searches** (`Speculation.*`; `AMPS_EARTH_SPECULATION`, default on
  for accelerators).  The sequential searches (UPPER_SCAN from Rmax down to the first
  forbidden sample, the bisections, the PENUMBRA_SCAN grid, BINARY endpoints) keep at
  most one trajectory in flight per location.  The candidates they may visit are
  evaluated concurrently on fibers and the results are walked in the original order;
  candidates the search can no longer use are cancelled (`BatchTracer::Cancel`,
  `eg_dev_cancel`, applied at the next launch).  Exceptions are rethrown only when the
  walk reaches the candidate; order-dependent side effects (diagnostic accumulators) are
  deferred (`DeferOrRun`) and replayed in visiting order.
* **Interior trilinear fast path** (`EngineMesh.h` `InterpolateInteriorTrilinear`).  When
  the AMPS stencil selection gives the trilinear stencil of the containing leaf and all
  8 cells lie inside it, B is accumulated directly with the same weights and order.
* **Stencil-cell cache** (`EngineMesh.h` `BuildStencilCache`; `AMPS_EARTH_MESH_CACHE`).
  At every mesh upload, for every used leaf, the owner leaf and cell index of each of the
  296 ghost cells of its 6x6x6 stencil lattice, the tree search of the multiblock
  support points and the 26 face/edge/corner neighbours are computed with the uncached
  functions and stored (about 2.4 KB per leaf, 0.8 GB for 330k leaves); interior cells
  are verified to resolve to the leaf itself.  The interpolation then replaces these
  tree searches by table lookups.  The AMPS core has an equivalent per-thread cache for
  the legacy path (`src/pic/pic_interpolation_routines.cpp`).
* **Normalization rule** (2026-10-23, also in the AMPS core).  A row stencil is
  renormalized exactly when a contribution with a positive weight was not added
  (previously: when the trilinear row did not have 8 entries); `Trilinear()` and
  `MultiBlock()` must stay in step with the core so that the engine reproduces the
  legacy row stencil bitwise.
* **Rank/thread override** (`test/tools/run_srcearth_tests.sh --np/--nt`, done by
  `test/tools/sublist.py`).  Rewrites the `-np`/`-nt` of the commands of srcEarth's list,
  e.g. one rank per GPU with more host threads; tests matching `--np-nt-exclude`
  (default C13, C18) keep their own counts.  (Until 2026-10-27 this was an option of a
  srcGeospace copy of the test runner; the tests and the runner are srcEarth's now.)

