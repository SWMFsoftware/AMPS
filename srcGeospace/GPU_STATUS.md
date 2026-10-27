# srcGeospace GPU development: status

As of **2026-10-27 12:30 CDT**, when the GPU work was parked.  Repository
`/data/vtenishe/GEO/AMPS`; builds and results in `/data/vtenishe/GEO/AMPS_runs/srcGeospace`
(below `$RUNS`).  Source paths are relative to `srcGeospace/`.

This file is the starting point for resuming the GPU work.  The state of the whole
application and its CPU validation is in [STATUS.md](STATUS.md); the history and the
reasoning behind decisions are in [handoff.md](handoff.md); the engine design is in
[doc/DESIGN.md](doc/DESIGN.md); build internals and all environment variables are in
[doc/BUILD.md](doc/BUILD.md).

---------------------------------------------------------------------------------------

## 1. Summary

* **Decision (user, 2026-10-27): keep the GPU development for the future and stop the GPU
  runs.**  All GPU processes were stopped at 12:20
  (`$RUNS/results/2026-10-26_parallel/lock_events.log`, `[manual]` line of 12:21:07).
  Nothing runs on the GPUs.
* **The engine stays in use.**  `gpu/` (28 files, 7,499 lines) is committed with
  srcGeospace in `465b15b` (2026-10-27 02:20); the working tree of `gpu/` equals that
  commit.  The same engine is the CPU path of srcGeospace: every CPU `amps` runs its
  trajectories through `gpu/` on CPU threads.  What is parked is the nvc++/OpenACC
  device pool (`gpu/EngineDevice.cpp` built by `scripts/build_gpu_engine.sh`) and the
  GPU-specific tuning, not the `gpu/` directory.
* **Correct except C3 3d MESH.**  On every test that finished on one GPU (54 in double,
  55 in float), the outputs equal those of the CPU pass of the same precision, apart from
  input decks with other rank/thread counts, AMPS mesh cache files and differences of at
  most 1.1e-13 in F5/F11 sums.  The exception is **C3 3d MESH**, whose cutoffs differ
  from the CPU in double and float (section 3.4).  Pass/fail outcomes equal those of the
  CPU passes on every finished test.
* **Slow.**  One RTX A6000 (1 rank x 32 host threads) was slower than the CPU runs
  (4 or 8 ranks x 16 threads) on the Mode3D decks: whole tests 7.6 times (double) and
  5.3 times (float) over 53 tests, tracing phases of the mesh decks 5-83 times
  (section 4).  Gridless decks were about even in double (0.96) and 2.4 times slower in
  float; the C6 IGRF decks among them trace on CPU threads in the GPU build too.  C18
  (4 ranks x 1-2 threads by design) was the one test clearly faster in the GPU build
  (552 against 1,285 s; not analysed).
* **Why it is slow** (evidence in section 4.5): few trajectories in flight (110-2,600
  active per kernel launch measured; the design assumes 1e4-1e5), a heavy kernel
  (255 registers and a 6.7 KB stack per thread in double, at most 256 resident threads
  per SM), and latency-bound mesh interpolation.  In the profiled C1 3d MESH runs the host
  waited for the GPU 99.9 % of the time.  On C1 3d MESH with one rank, one A6000
  delivered about 1/13 of the trajectories per second of **one** CPU thread.
* **Validation is incomplete**: C19 did not complete and C17, C16, C15, C13, C12 and C8
  did not run, in either precision; C2 3d MESH did not finish in double; the 2-GPU passes,
  the rebuild of the double binary with the allocator fix and the C9 tail study never
  started (section 3.2).
* **The merged GPU engine was never checked trajectory by trajectory**
  (`AMPS_EARTH_ENGINE_VERIFY=1`) on a GPU.  The last such check was the pre-merge GPU
  build of 2026-10-22 (0 mismatches).
* **Single-rank runs need the AMR allocator fix** (`b06e7e7`).  `builds/merged_gpu_double`
  predates it (2026-10-23); its single-rank Mode3D setup is up to 24 times longer (C19:
  1 h 36 min against 4 min for the float build with the fix).  Rebuild before any further
  GPU run.
* **Several documents overstated the GPU results** (bitwise claims, "never run on a GPU");
  corrected on 2026-10-27, list in section 5, item 9.

## 2. What exists

### 2.1 Engine (details: `doc/DESIGN.md`)

| part | where | GPU relevance |
|---|---|---|
| fibers: each driver work unit (location, task, energy) runs in a fiber; a trajectory request suspends it | `gpu/Fiber.*`, `DriverSupport.*`, `TaskStream.h` | search logic unchanged; one trajectory in flight per fiber |
| batch front end: request queue, contexts, routing, statistics | `gpu/BatchTracer.*` | `Progress()` is the fiber idle hook |
| device pool: persistent slots, K-step kernel `LaunchKernel`/`AdvanceOne`, gather kernel, asynchronous launch/wait | `gpu/EngineDevice.cpp`, `EnginePool.h` | the only nvc++ translation unit; OpenACC `parallel loop gang vector vector_length(64) async(1)` |
| host pool for Geopack/Tsyganenko fields | `gpu/EngineHost.*`, `HostThreadPool.*` | always CPU threads, also in GPU builds |
| resumable tracer, per-trajectory adaptive dt, movers, boundaries, Shue, exit kinematics | `gpu/EngineStep.h`, `EngineMovers.h`, `EngineField.h` | exact port of the legacy code (`DESIGN.md` section 6) |
| trap detector with a fixed-capacity per-slot arena | `gpu/EngineTrap.h` | always SI double; uses `atan2`/`hypot` and `efmod` |
| Mode3D mesh flattening and upload; AMPS row-stencil port | `gpu/Mode3DMeshExport.*`, `EngineMesh.h` | tree search and two row stencils (capacity 80 entries, at most 72 used) per thread |
| speculative cutoff searches with cancellation | `gpu/Speculation.*`, `eg_dev_cancel` | default on GPUs only |
| performance additions of 2026-10-22/23 | `FlushPending` (batched request upload), assembly fiber switch and stack pool, interior trilinear fast path, stencil-cell cache (`BuildMeshCache`), normalization rule | `DESIGN.md` section 13 |
| driver glue | `gridless/CutoffRigidityGridless.cpp`, `gridless/DensityGridless.cpp`, `3d/CutoffRigidityMode3D.cpp`, `3d/DensityMode3D.cpp`, `3d/GlobalMagneticField.cpp`, `3d/Mode3D.cpp`, `makefile` | |

Device-only code paths (`EG_HAVE_OPENACC` in `gpu/EngineDevice.cpp`): the scatter kernel
of `FlushPending` (writes new requests and resets the slot from `slotTemplate`; the CPU
build instead sets `done`/`lastNode` on the host slot), `ApplyCancels` as an `async(1)`
kernel, mesh arrays and the stencil-cell cache copied to device memory (the CPU double
build aliases the host arrays), and the nvc++ compilation of everything reached from the
kernel.  These are the first places to look for GPU-only differences (section 3.4).

### 2.2 Toolchain and build

* Machine: Intel i9-10980XE (18 cores, 36 hardware threads), 251 GiB; 2 x RTX A6000
  (49,140 MiB, compute capability 8.6, driver 550.163.01).  NVIDIA HPC SDK 24.5 (nvc++,
  CUDA 12.4) in `/data/vtenishe/openacc` (`source /data/vtenishe/openacc/env.sh`; the
  build scripts source it themselves).
* `scripts/build_amps.sh --stage DIR --gpu [--float] -j 12` (clean build; wipes the stage):
  1. `scripts/build_gpu_engine.sh` compiles only `gpu/EngineDevice.cpp`:
     `nvc++ -std=c++17 -O2 -fPIC -acc=gpu -gpu=cc86,cuda12.4,nordc,nofma
     [-DGEOSPACE_PRECISION_FLOAT] -Minfo=accel` -> `<stage>/gpu_engine/libGeospaceEngine.so`
     (flags recorded in `<stage>/gpu_engine/build_info.txt`);
  2. the g++ AMPS build with `GEOSPACE_DEVICE_OBJ=` and `EXTRALINKEROPTIONS=<library and
     RUNPATHs>` on the make command line; no repository file is edited.
* `nofma` is the default so that the double GPU build can reproduce the legacy tracer
  bit for bit (nvc++ 24.5 contracts `a*b+c` at -O2); `--gpu-opts "--fma"` allows FMA.
  All GPU builds in `$RUNS/builds` use `nofma`.
* One OpenMP runtime: `amps` needs `libgomp.so.1`, which its RUNPATH resolves to the SDK's
  `compilers/lib/libgomp.so.1` (= `libnvomp`), shared by AMPS's OpenMP regions and the
  OpenACC runtime (`ldd builds/merged_gpu_double/amps`).  All GPU runs of 2026-10-22 to
  10-27 used this arrangement; no failure was traced to it.
* nvc++ 24.5 has no device `fmod`: `Geospace::GPU::efmod()` (`gpu/EngineConfig.h`, patch
  `doc/patches/0001-device-fmod.diff`, already applied), exact by construction and tested
  on 4.1e7 double and 4.1e7 float cases (`tools/device_math/efmod_exactness.cpp`).

### 2.3 Running on GPUs

* `scripts/amps_gpu_select.sh` replaces the executable: `GEOSPACE_AMPS=<gpu amps>`,
  `GEOSPACE_GPUS=0,1`; each local MPI rank gets `CUDA_VISIBLE_DEVICES` = entry
  (rank mod n).  The device pool uses the default device of the process.
* `test/tools/run_srcearth_tests.sh --np N --nt M` (until 2026-10-27: an option of
  srcGeospace's copy of the runner, `test/test_runner.py`) rewrites `-np/-nt` of the commands; C13
  and C18 keep their own counts (`--np-nt-exclude`).  The 2026-10-26/27 passes used
  `--np 1 --nt 32` (1 GPU); `--np 2 --nt 16` was planned for 2 GPUs.
* Defaults on an accelerator (`gpu/BatchTracer.cpp` `EnsurePools()`, `Speculation.cpp`):
  65,536 device-pool slots, K = 2,048 iterations per launch, fiber budget = pool capacity
  (`Fiber::SetDefaultMaxFibers`), speculation on (depth 4), asynchronous launches on,
  stencil-cell cache on; host pool max(64, 16 x host threads) slots with K = 64.  CPU
  builds: max(256, 64 x host threads) slots, K = 256, 4,096 fibers, speculation off.
* Environment variables: `doc/BUILD.md` section 7.1 (verified against `gpu/`).  Useful
  for GPU work: `AMPS_EARTH_ENGINE_VERBOSE=1` (fiber and pool statistics, including the
  "progress split" of wait/launch/collect time), `AMPS_EARTH_ENGINE_VERIFY=1`,
  `AMPS_EARTH_MESH_VERIFY=1`, `AMPS_EARTH_SPECULATION=0|1`, `AMPS_EARTH_MESH_CACHE=0`,
  `AMPS_EARTH_STEPS_PER_LAUNCH`, `AMPS_EARTH_POOL_CAPACITY`, `AMPS_EARTH_FIBERS`,
  `AMPS_EARTH_ASYNC=0`.  `AMPS_EARTH_ENGINE_TRACE_S` (periodic pool trace) exists only in
  the private copy `$RUNS/work/c9_gpu_tail/srcGeospace` (section 4.6).
* Host memory of single-rank runs (peak resident memory, `memlog_proc.txt`): C2 3d MESH
  99.4 GiB (float, with the fix), C19 78.4 GiB (double); most other Mode3D decks 22-30 GiB.
  Device memory seen with `nvidia-smi` during the runs (not logged): about 1.3 GiB for
  C9 FULL_SCAN (float) and 3.9 GiB for C19 (double), of 48 GiB.

### 2.4 Equivalence checks

* `AMPS_EARTH_ENGINE_VERIFY=1`: every engine trajectory is also traced by the legacy code
  on the host and compared bit for bit; summary line `[earth-engine][verify]
  trajectories compared=N bitwise mismatches=M`.  On a GPU this compares device results
  with the host legacy tracer.
* `AMPS_EARTH_MESH_VERIFY=1`: node search and interpolation of the flattened mesh against
  AMPS after each upload (host-side `EgMesh`, not the device copy).
* `test/UEngineEquivalence` (in the list) compiles engine headers with g++ against the
  legacy code; it does not run device code.
* `test/tools/compare_outputs.py`, `compare_builds.py`: run root against run root.

## 3. Validation

### 3.1 GPU runs so far

| date | build | what | result | results (`$RUNS/results/`) |
|---|---|---|---|---|
| 2026-10-22 | `gpu_double`, `gpu_float` (pre-merge engine; the only GPU builds at that time) | smoke tests; in double with `AMPS_EARTH_ENGINE_VERIFY=1` | every verify summary 0 bitwise mismatches (gridless and Mode3D); mesh self-check 58,793 interpolated samples, 0 mismatches | `2026-10-22_srcGeospace_gpu_smoke/` |
| 2026-10-22 | `gpu_double` | 1-GPU pass of `list_cpu` with the list's own rank/thread counts (C3 3d MESH: 4 x 16 on one GPU) | stopped after 35 of 62 tests; #020 C2 3d MESH killed (exit 137, 792 s); #022 C3 3d MESH 925 s, outputs equal to today's CPU | `2026-10-22_srcGeospace_gpu/` |
| 2026-10-22 | `exp_gpu_*` (experimental copy) | speculation, fast path, K, float on C1 3d MESH; profiling of F11 and C1 MESH | section 4.4 | `2026-10-22_speculation_check/`, `2026-10-22_gpu_profiling/` |
| 2026-10-22 | `exp_gpu_double`, K = 256, 4 x 8 | C19 | stopped by the user during snapshot 1 after about 50 min | `2026-10-22_c19_gpu/STOPPED.txt` |
| 2026-10-26/27 | `merged_gpu_double`, `merged_gpu_float` | 1-GPU passes, 1 rank x 32 threads | section 3.2 | `2026-10-26_merged_gpu/` |

### 3.2 The 1-GPU passes of 2026-10-26/27

List `list_cpu_locked`: the full `test/list` without SA01-SA03 (62 tests), each command
prefixed by `2026-10-26_parallel/lockrun.sh` because the CPU passes ran at the same time
(STATUS.md section 8).  Runner: `test_runner.py -j 1 --np 1 --nt 32`.  CPU references:
`2026-10-26_merged_cpu/{double,float}` (builds `cpu_{double,float}_2026-10-26`, without
the allocator fix, six tests at a time).

| run root | GPU | binary | time | tests | outcome |
|---|---|---|---|---|---|
| `gpu1_double` | 0 | `merged_gpu_double` (2026-10-23, no allocator fix) | 10-26 16:48 - 10-27 02:44 | #001-#052 | as listed, except #020 C2 3d MESH killed at 2,817 s in single-rank setup (exit 143, user decision).  Runner killed at 23:50 after #051 (incident, STATUS.md section 7); #052 C9 FULL_SCAN finished as an orphan, PASS by its own log (no exit code) |
| `gpu1_double_rest` | 0 | same | 10-27 02:44 - 12:20 | #053-#056 | #053 C9 DIRECT_ACCESS P; #054 C10 F (as on the CPU, identical scores); #055 C18 P; #056 C19 killed at snapshot 2 of 3, task 2,535 of 3,530 |
| `gpu1_float` | 1 | `merged_gpu_float` (2026-10-26 20:05, with the fix) | 10-26 23:50 - 10-27 12:20 | #001-#056 | F: #016 C1 3d MESH, #020 C2 3d MESH, #052 C9 FULL_SCAN, #053 C9 DIRECT_ACCESS, #054 C10 = exactly the CPU float failures among these tests.  #056 C19 waited 07:07-12:14 for the heavy-test lock, then was killed 2.5 min into tracing |

**Not run** (both precisions): #056 C19 (unfinished), #057 C17, #058 C16, #059 C15,
#060 C13, #061 C12, #062 C8; in double also #020 C2 3d MESH.  Never started: the 2-GPU
passes `gpu2_double` and `gpu2_float`, the rebuild of `merged_gpu_double` with the fix
(`2026-10-26_allocfix/rebuild_gpu_double.sh`), the C9 tail study.  SA01-SA03 were never
part of the GPU list.  The three run roots have **no `test_runner_report.*` files**
(runners killed): results are in `<run root>_runner_stdout.log`, `test_runner_logs/` and
`test_output/`; `stream.log` records the stream.

### 3.3 Outputs against the CPU

`test/tools/compare_outputs.py` (CPU pass of the same precision as reference):

| tests | identical | identical after paths | nondeterministic (AMPS mesh cache files) | numeric | different |
|---|---|---|---|---|---|
| double #001-#051 (`gpu1_double_vs_cpu_double_partial.txt`) | 575 | 8 | 44 | 104 | 1 |
| float #001-#051 (`gpu1_float_vs_cpu_float_partial.txt`) | 473 | 2 | 40 | 104 | 1 |
| double #052-#055 (C9 FULL_SCAN, C9 DIRECT_ACCESS, C10, C18) | 775 | 15 | 136 | 0 | 0 |
| float #052-#055 | 775 | 15 | 136 | 0 | 0 |

* The 104 numeric differences: 42 input decks (other `np`/`nt`), 49 files of F11 and 9 of
  F5 (largest absolute difference 1.1e-13; STATUS.md attributes them to the summation
  order with one rank instead of four), and 4 files of C3 3d MESH.  The one "different"
  file is `C3_3d_mesh/C3_summary.csv`.
* The rows for #052-#055 were made for this document with
  `compare_outputs.py <cpu root> <gpu root> --subdir test_output/<test dir>` (test dirs
  `C9_full_85`, `C9_direct_85`, `C10_direct`, `C18_parallel_determinism`); the partial
  files cover #001-#051 only.
* Identical to the CPU apart from input decks and mesh cache files, among others: C1 3d
  MESH, C2 3d MESH (float), the gridded C6/C7 decks, C9 (both), C10, C11, C14, C18.

### 3.4 Open issue: C3 3d MESH differs from the CPU

The test (list line 64, `run_C3.py --mode 3d --mode3d-field-eval MESH -np 4 -nt 16`):
one Mode3D shell at 9000 km (12 longitudes x 7 latitudes = 84 locations), dipole field
on the AMR mesh, BINARY and UPPER_SCAN searches.  It passes on the GPU, but the cutoffs
of `binary/` and `upper_scan/cutoff_3d_shells.dat` differ:

| comparison (84 points each) | BINARY: points that differ; median / max relative difference | UPPER_SCAN: same |
|---|---|---|
| GPU double vs CPU double | 38; 0.3 % / 287 % (a point collapses to the 0.0433 GV low boundary on one side only) | 54; 1.0 % / 17 % |
| GPU float vs CPU float | 34; 0.2 % / 290 % | 48; 0.4 % / 10 % |
| CPU float vs CPU double (for scale) | 35; 0.2 % / 282 % | 49; 0.7 % / 9.5 % |
| CPU double 10-26 vs CPU double 10-23 | 0 | 0 |

`C3_summary.csv` (the 16 points at +-60 degrees) differs in the same way, including its
`collapsed_to_low_boundary` flags.

What is known:
* The size of the difference is that of a rounding-level change: this deck amplifies any
  last-bit difference (CPU float against CPU double differs just as much).  But the double
  GPU build is meant to be bitwise identical to the CPU (`nofma`), and the 2026-10-22 GPU
  build was: its C3 3d MESH outputs equal the 2026-10-26 CPU outputs.
* CPU C3 is reproducible (10-22, 10-23, 10-26 identical).
* Excluded **on the CPU** (`$RUNS/work/c3_spec_check/`, binary `cpu_double_2026-10-26`):
  speculative evaluation (run `cpu_double_spec1`; the environment is not recorded there,
  handoff.md section 11 says `AMPS_EARTH_SPECULATION=1`) and the rank count
  (`cpu_double_np1`: `--np 1 --nt 32`).  Both equal the CPU reference.
* Not excluded: **two things changed between the matching 2026-10-22 GPU run and the
  differing 2026-10-26 runs**: the engine (merge of 2026-10-23) and the run configuration
  (4 ranks x 16 threads on one GPU then; 1 rank x 32 threads now, which also means a
  different pool load and speculation timing).  The merged GPU build at 4 ranks and the
  pre-merge GPU build at 1 rank were never run.  So STATUS.md's conclusion "a GPU-only
  part of the merged engine" is likely but not established.
* The engine objects compiled by g++ are byte-identical between `merged_gpu_double` and
  `cpu_double_2026-10-26` (`build/main/gpu/*.o`); only the nvc++ device library differs.
  Engine sources changed since the 2026-10-22 GPU build: `BatchTracer.*`,
  `DriverSupport.cpp`, `EngineDeviceAPI.h`, `EngineDevice.cpp`, `EngineMesh.h`,
  `EngineTypes.h`, `Fiber.*`, new `Speculation.*`.  Unchanged: `EngineTrap.h`,
  `EngineStep.h`, `EngineMovers.h`, `EngineField.h`.

Hypotheses (not tested): a device-only path of section 2.1 (slot reset from
`slotTemplate` in `FlushPending`, `ApplyCancels` on the device, cached or fast-path
interpolation compiled by nvc++), or a scheduling-dependent effect of speculation and
cancellation that appears only with the GPU pool size and K.  Isolation plan in
section 5, item 2.

## 4. Performance

All numbers of 2026-10-26/27 are net of lock waits (start/end lines of
`2026-10-26_parallel/lock_events.log`; per-test table made with the scratch script
`lock_times.py`, not kept in the repository).  **Not a clean comparison**: the CPU passes
ran six tests at a time, mostly 4 x 16 = 64 threads each, on 36 hardware threads; the GPU
tests used 1 rank x 32 threads and shared the machine with the CPU passes: GPU double
#001-#051 ran next to the CPU double pass, its C9 FULL_SCAN (19:34-02:44) next to the rest
of the CPU double pass and the whole CPU float pass, GPU float #001-#051 next to the CPU
float SA tests; the later tests ran next to the other GPU pass only.  The GPU double
binary lacks the allocator fix, the float binary has it.

### 4.1 Whole tests, by group (tests finished in all four passes)

| group (tests) | CPU double | GPU double | GPU d / CPU d | CPU float | GPU float | GPU f / CPU f | GPU f / GPU d |
|---|---|---|---|---|---|---|---|
| unit tests, C4, C5, F4 self-tests (22) | 152 s | 162 s | 1.07 | 118 s | 142 s | 1.20 | 0.88 |
| gridless cutoffs: C1-C3 gridless, C3 UPPER_SCAN, C6 x3 (7) | 301 s | 290 s | 0.96 | 179 s | 437 s | 2.44 | 1.51 |
| C1-C3 3d ANALYTIC, C11, C14, F1-F5, F11, F12, F15, F16 (14) | 149 s | 495 s | 3.3 | 142 s | 267 s | 1.9 | 0.54 |
| Mode3D MESH/GRIDDED: C1, C3, C14 MESH, C6 GRIDDED x3, C7 (7) | 439 s | 5,583 s | 12.7 | 416 s | 2,920 s | 7.0 | 0.52 |
| C9 FULL_SCAN, C9 DIRECT_ACCESS, C10 (3) | 3,882 s | 31,086 s | 8.0 | 3,741 s | 20,518 s | 5.5 | 0.66 |
| **total (53)** | **4,923 s** | **37,616 s** | **7.6** | **4,596 s** | **24,284 s** | **5.3** | **0.65** |
| C18 (4 ranks x 1-2 threads, own counts) | 1,285 s | 552 s | 0.43 | 1,262 s | 471 s | 0.37 | 0.85 |
| C2 3d MESH | 672 s | killed at 2,816 s | | 599 s | 266 s | 0.44 | |

* The gridless C6 subtests use IGRF, which runs in the host pool (CPU threads) also in the
  GPU build; their GPU-build times measure load and the 1 x 32 rank layout, not the GPU.
* C18 runs Mode3D MESH tracing on the GPU (OpenACC backend in its log) with 4 ranks of
  1-2 threads; it is the one test clearly faster in the GPU build.  Not analysed.  Some
  tests were faster overall only through a shorter setup or less load: C2 3d MESH and
  C9 DIRECT_ACCESS in float (section 4.2), C6 COMPLETE (gridless) in double.
* STATUS.md section 5.3 (6.6 and 3.2) used other test sets; the numbers above supersede it.

### 4.2 Setup and tracing of the Mode3D decks

Setup = input deck written -> `standalone_mode3d_manifest*.json` written (mesh, field,
mesh upload); tracing = manifest -> last `cutoff_3d*.dat` (file times, 1 s resolution;
"at least" takes a CPU tracing time of 1 s as 2 s).  CPU: 4 x 16 threads (C9, C10:
8 x 16).

| test | CPU double: setup / tracing | GPU double: setup / tracing | GPU float: setup / tracing (CPU float tracing) | tracing GPU d / CPU d | tracing GPU f / CPU f |
|---|---|---|---|---|---|
| C1 3d MESH | 73 / 2 s | 430 / 68 s | 72 / 49 s (2 s) | about 34 | about 25 |
| C2 3d MESH | 663 / 6 s | killed in setup | 227 / 33 s (6 s) | | about 5.5 |
| C3 3d MESH, BINARY / UPPER_SCAN | 15 / 1 s, 15 / 1 s | 35 / 113 s, 53 / 50 s | 17 / 88 s, 17 / 35 s (1 s) | at least 56 / 25 | at least 44 / 17 |
| C14 cross-solver Mode3D | 25 / 1 s | 124 / 65 s | 35 / 47 s (1 s) | at least 32 | at least 23 |
| C6 INITIAL GRIDDED | 68 / 30 s | 345 / 1,192 s | 62 / 797 s (33 s) | 40 | 24 |
| C6 COMPLETE GRIDDED | 54 / 14 s | 538 / 737 s | 88 / 539 s (20 s) | 53 | 27 |
| C6 MODERN GRIDDED | 49 / 16 s | 248 / 779 s | 58 / 565 s (17 s) | 49 | 33 |
| C7 IZMIRAN GRIDDED | 49 / 6 s | 268 / 498 s | 57 / 368 s (6 s) | 83 | 61 |
| C9 DIRECT_ACCESS, 7 intervals | 829 / 37 s | 1,493 / 565 s | 403 / 376 s (45 s) | 15 | 8.4 |
| C9 FULL_SCAN, 7 intervals | 481 / 1,322 s | 1,800 / 23,985 s | 467 / 17,904 s (1,342 s) | 18 | 13 |
| C10, 10 snapshots | 1,052 / 52 s | 2,156 / 1,028 s | 577 / 729 s (59 s) | 20 | 12 |

* The long GPU-double setups are mostly the single-rank allocator problem (STATUS.md
  section 3): with the fix the float setups are about those of the CPU.  C9 DIRECT_ACCESS
  in float was faster overall on the GPU (796 against 867 s) only because of the shorter
  setup; its tracing was 8.4 times slower.
* GPU float tracing takes 0.67-0.78 of GPU double tracing on these decks.

### 4.3 Long tests

* **C9 GRIDDED FULL_SCAN** (`test_output/C9_full_85/gridded/<interval>/`): tracing per
  interval GPU double 3,205-3,571 s, GPU float 2,355-2,712 s, CPU double 144-247 s, CPU
  float 144-252 s.  Whole test: GPU double 25,804 s (19:34-02:44), GPU float 18,389 s
  (01:17-06:23), CPU 1,818 / 1,781 s.  Outputs identical to the CPU in both precisions.
* **C19 SMOKE** (T05, GRIDDED, DIRECT_ACCESS, 3 snapshots, 3,530 tasks per snapshot;
  file times in `test_output/C19_goes_epead_ew/gridded/t05/batch_direct_access/`):

  | phase | CPU double (4 x 8) | GPU double (1 x 32, no fix) | GPU float (1 x 32, fix) |
  |---|---|---|---|
  | setup | 18:47:20-18:53:55 (395 s) | 04:21:42-05:57:55 (5,773 s) | 12:14:28-12:18:26 (238 s) |
  | snapshot 1 tracing | 416 s | 05:57:55-10:17:06 (15,551 s, 37 times the CPU) | killed after 147 s at task 0 |
  | snapshot 2 tracing | 583 s | 10:19:18-12:20:55, killed at task 2,535 of 3,530 | -- |
  | snapshot 3 tracing | 383 s | not reached | -- |
  | whole test | 2,059 s | -- | -- |

  A complete 1-GPU C19 would take on the order of 10 h (estimate from the two snapshots).
* **C2 3d MESH**: GPU float with the fix 266 s (227 s setup) against 100 s for the CPU
  with the fix on 4 ranks (STATUS.md section 3); outputs byte-identical to the CPU float run.
* **C10**: GPU double 3,210 s, float 1,333 s, CPU 1,161 / 1,093 s; identical outputs
  (the test fails on both, as in srcEarth).

### 4.4 Earlier measurements (2026-10-22, experimental and pre-merge builds)

* **One rank, one tracing thread, C1 3d MESH** (UPPER_SCAN, 168 locations;
  `$RUNS/work/simd_eval/cache_test/RESULTS.md`): CPU double 13.17 s without and 5.34 s
  with the stencil-cell cache (8,888 trajectories); GPU double 350.0 s and 125.7 s
  (16.3k trajectories, more because of speculation).  That is about 130 trajectories/s on
  the GPU against about 1,660/s on one CPU thread.  The cache made the GPU 2.8-3.0 times
  faster, outputs identical.
* **Launch size, C1 3d MESH, 4 ranks x 16 threads on one GPU**
  (`2026-10-22_speculation_check/gpu0_fastpath_double*`): default K = 2,048: 487.3 s,
  tracing 336-412 s per rank, 10-19 launches, 602-693 active slots per launch; run
  `gpu0_fastpath_double_k256` (57-100 launches, 182-202 active): 345.2 s, tracing
  230-266 s.  The smaller K was 1.4 times faster; the default was not changed.
  `DESIGN.md` section 12 quotes 345 s.  Float with K = 2,048: 282.8 s.
* **GPU-bound**: in those runs the "progress split" was wait = 335.5 s of 335.8 s;
  submit, launch, collect and deliver together took about 0.2 s.
* **C1 3d MESH float, 1 rank, 16 host threads, before the cache**
  (`2026-10-22_gpu_profiling/c1mesh/run_*_float.log`): Mode3D cutoff wall time CPU 1.50 s,
  GPU 170.8 s (117 launches, 109 active slots per launch).
* **F11 gridless density, dipole field, 1 rank x 16 threads**
  (`2026-10-22_gpu_profiling/f11/`): kernel+collect per batch CPU double 0.60 / 0.94 s,
  GPU double 3.58 / 3.64 s, GPU float 0.98 / 1.03 s; 1,444 / 2,601 active slots per
  launch.  Nsight Systems (`prof_gpu_float.sqlite`, `prof_exp.sqlite`): 60 kernel launches,
  15.7 ms each on average, on average 39 gangs x 64 threads, 168 registers (float);
  295,117 memory copies before the batched request upload, 218 after.
* **C19, 4 ranks x 8 threads on one GPU, K = 256**: snapshot 1 still running after about
  50 min, against 17.7 min on the CPU of that time (`2026-10-22_c19_gpu/STOPPED.txt`).

### 4.5 Why the GPU is slow

Measured:
1. **Too few trajectories in flight.**  Active slots per launch: 109 (C1 3d MESH float,
   1 rank, before the cache), 182-693 (C1 3d MESH, per rank of 4), 1,444-2,601 (F11);
   kernel grids of about 39 gangs x 64 threads (F11).  The design goal is 1e4-1e5
   trajectories in flight (`DESIGN.md` section 1).  The Mode3D drivers keep one fiber per
   location or task (C1: 168 locations, C3: 84, C9: 1,092, C19: 3,530 tasks per
   snapshot) and the cutoff searches are sequential per location; speculation adds
   candidates, not locations.
2. **Heavy kernel.**  255 registers and a 6,704-byte stack frame per thread in double
   (168 and 5,296 bytes in float); at most 4 gangs (256 threads) per SM
   (`DESIGN.md` section 5).  The mesh interpolation does tree searches and builds row
   stencils of up to 72 cells in per-thread local memory.
3. **The host waits for the GPU** (progress split above), so host-side work (fibers,
   search logic, about 4,200 serial instructions per trajectory, `DESIGN.md` section 3)
   is not the limit.
4. **Precision is not the limit on mesh decks.**  GPU float tracing takes 0.67-0.78 of
   double on the mesh decks, but 0.27 on the analytic-dipole F11 kernel (0.98 against
   3.58 s).  The A6000's FP64 rate is a small fraction of its FP32 rate
   (`PRECISION.md` section 6), so the mesh decks are bound by something else (hypothesis:
   memory latency of the interpolation and local-memory traffic).
5. **Launch granularity.**  K = 256 beat K = 2,048 by 1.4 times on C1 3d MESH.  Hypothesis:
   a launch lasts until every active slot has done K iterations, so finished trajectories
   and the next search steps wait for the longest slots.

Not measured: occupancy, register spills, local-memory traffic and divergence of the
kernel (no Nsight Compute run); throughput with a saturated pool.

### 4.6 C9 FULL_SCAN tail (prepared, never run)

Observation (handoff.md section 11, no timestamps, not verified): the progress lines of an
interval stay at 97.6 % (task 1,066 of 1,092) for many updates.  From the code
(`3d/CutoffRigidityMode3D.cpp` lines 2984-3004): PENUMBRA_SCAN traces the
160-point coarse rigidity grid of a location concurrently, but the transition refinement
(`RefineCutoffTransition`) and `IntegrateEffectiveCutoff` are sequential, one trajectory
in flight per location.  In the 2006-12-14 14:31 interval 151 of 1,092 locations have
transitions in 0.3-1.35 GV, 297 are all allowed and about 644 all forbidden.

Study set-up in `$RUNS/work/c9_gpu_tail/`: private copy `srcGeospace/` whose only engine
change is a periodic pool trace in `gpu/BatchTracer.{h,cpp}` (`AMPS_EARTH_ENGINE_TRACE_S=10`:
active slots, queue, in-flight requests, launches, mean active per launch, launch wait,
delivered trajectories, mean/max steps, count >= 1e5 steps, terminations); builds
`builds/c9tail_gpu_{double,float}` (with the allocator fix); `inputs/` = deck, driver
file and command of the 14:31 interval; `run_interval.sh` (section 6.3).  `runs/` does not
exist: the study never ran.  Candidate remedies once the tail is understood: speculative
evaluation of the refinement steps, or moving the last trajectories to the host when few
remain in flight.

## 5. Open issues and next steps (suggested order)

1. **Decide whether a GPU can pay off for these decks before more porting.**  Evidence:
   section 4.5.  Measure (a) the kernel with Nsight Compute (occupancy, local memory,
   divergence) on C1 3d MESH and C9, and (b) trajectories/s with a saturated pool (a deck
   with many more locations or directions, `AMPS_EARTH_FIBERS` and
   `AMPS_EARTH_POOL_CAPACITY` large), against 36 CPU threads.  Options listed in
   `DESIGN.md` section 12: split the kernel (field evaluation vs mover), move the row
   stencils out of the per-thread frame, `-gpu=maxregcount`, structure-of-arrays state;
   more concurrency per location (speculative refinement, item 7 (b) below).
2. **C3 3d MESH** (section 3.4).  Each GPU run takes a few minutes (254 s in double on
   2026-10-26).  In this order:
   * GPU with `AMPS_EARTH_ENGINE_VERIFY=1` (prints the first 20 mismatching trajectories);
   * GPU variants: `AMPS_EARTH_SPECULATION=0`, `AMPS_EARTH_MESH_CACHE=0`,
     `AMPS_EARTH_STEPS_PER_LAUNCH=256`, `AMPS_EARTH_ASYNC=0`, and without `--np/--nt`
     (4 ranks on one GPU, as on 2026-10-22);
   * CPU with GPU-like scheduling: `AMPS_EARTH_SPECULATION=1 AMPS_EARTH_POOL_CAPACITY=65536
     AMPS_EARTH_FIBERS=65536 AMPS_EARTH_STEPS_PER_LAUNCH=2048`, `--np 1 --nt 32`.
3. **Rebuild the GPU binaries from the current code** (HEAD includes `b06e7e7` and
   `465b15b`) and run the merged engine with `AMPS_EARTH_ENGINE_VERIFY=1` on the decks of
   the 2026-10-22 smoke tests (never done for the merged engine).
4. **Finish the validation**: C19, C17, C16, C15, C13, C12, C8 on 1 GPU in double and
   float; C2 3d MESH in double with the fixed binary; the 2-GPU passes; compare with the
   CPU passes (`compare_outputs.py`, `compare_builds.py`).  Plan for C19 (about 10 h on one
   GPU).
5. **C9 tail study** (section 4.6), then a fix proposal.
6. **Launch configuration**: K (2,048 default against 256), several ranks per GPU, pool
   and fiber sizes.
7. **Engine items not done** (handoff.md section 3.2): (b) direction blocks / fewer fiber
   switches; (c) multi-threaded host request generation (MPI is initialized
   `MPI_THREAD_FUNNELED`: only the main thread may call MPI); skip `BuildMeshCache` when
   only B changes (C10 uploads 10 snapshots).
8. **Two-GPU runs**: never run with the merged code (`amps_gpu_select.sh`, `--np 2 --nt 16`).
9. **Documents to correct** (done on 2026-10-27 13:00: each statement below was
   qualified in place, with a pointer to this file; kept here as the record):
   * `doc/DESIGN.md` section 0, rows "Device pool built by nvc++" ("double results bitwise
     identical to the CPU build") and "Performance additions" ("outputs identical ... CPU
     and GPU"): true for the 2026-10-22 build; the merged build differs on C3 3d MESH.
   * `BUILD.md` section 4: "Double GPU results are bitwise identical to the CPU double
     build" -- same.
   * `doc/DESIGN.md` section 12: "C1 3d MESH 601 -> 345 s" -- 345 s was measured with
     K = 256; the default K = 2,048 gave 487 s.
   * `doc/BUILD.md` section 5.2 ("not exercised at run time") and section 6.1 ("None of the
     GPU executables has been run"): outdated.
   * `doc/BUILD.md` section 7.1, practical notes: "By default only about 4096 trajectories
     are in flight" -- on a GPU the fiber budget defaults to the pool capacity (65,536).
   * `doc/PRECISION.md` section 6, item 5 ("Speed-ups have not been measured"): measured
     now (section 4 here).
   * `doc/README.md`: `gpu/` "6,576 lines"; `wc -l gpu/*` gives 7,499.
   * `STATUS.md` section 4.3: "UPPER_SCAN cutoffs differ by about 1 %" is the median (the
     maximum is 17 %), and the matching 2026-10-22 GPU run used 4 ranks (section 3.4 here).
10. **Cleanup** (when no longer needed): GPU builds `exp_gpu_*`, `gpu_double`, `gpu_float`
    (2026-10-22), `c9tail_gpu_*` (about 250 MB each); `/data/vtenishe/GEO/AMPS/srcEarth_gpu/`
    (old prototype with the old engine).  `/data/vtenishe/GEO/AMPS/srcEarth_cuda_prototype/`
    (untracked, 2026-10-20) is a separate standalone MPI/CUDA prototype (its own README and
    DESIGN.md; gridded dipole, no AMPS mesh); it was not part of this work and holds no
    recorded measurements.

## 6. Locations and how to resume

### 6.1 Builds (`$RUNS/builds/`)

| stage | date | contents |
|---|---|---|
| `merged_gpu_double` | 2026-10-23 12:21 | merged engine, double, `nofma`; **no allocator fix**; used by `gpu1_double`, `gpu1_double_rest` |
| `merged_gpu_float` | 2026-10-26 20:05 | merged engine, float, with the fix; used by `gpu1_float` |
| `c9tail_gpu_double`, `c9tail_gpu_float` | 2026-10-26 21:55/57 | instrumented copy `work/c9_gpu_tail/srcGeospace`, with the fix |
| `gpu_double`, `gpu_float` | 2026-10-22 10:54/56 | pre-merge engine (the C3-matching GPU run) |
| `exp_gpu_double`, `exp_gpu_float`, `exp_gpu_double_cache` | 2026-10-22 | experimental copy during the 10-22 tuning |

The engine sources compiled into `merged_gpu_double`, `merged_gpu_float` and
`cpu_double_2026-10-26` are identical (their stage copies `build/main/gpu/` match) and
equal the committed `gpu/` apart from AMPS configure processing.

### 6.2 Results and work directories (`$RUNS/`)

| path | contents |
|---|---|
| `results/2026-10-26_merged_gpu/` | `gpu1_double`, `gpu1_double_rest`, `gpu1_float` run roots; `<run root>_runner_stdout.log`; `stream.log`; `list_cpu`, `list_cpu_locked`, `list_cpu_locked_after179` (lines 182-228 = the 10 tests after #052); `gpu1_{double,float}_vs_cpu_*_partial.{csv,txt}` |
| `results/2026-10-26_merged_cpu/{double,float}` | CPU references |
| `results/2026-10-26_parallel/` | `gpu_stream.sh`, `gpu_stream2.sh` (the stream), `lockrun.sh` and `make_locked_list.py` (cross-stream lock), `lock_events.log`, `memmon.sh`, `memlog.txt`, `memlog_proc.txt`, `c19_overlap_helper.sh` |
| `results/2026-10-26_allocfix/` | allocator-fix checks; `rebuild_gpu_double.sh` (never ran) |
| `results/2026-10-22_*gpu*`, `2026-10-22_speculation_check/`, `2026-10-22_gpu_profiling/` | runs of section 3.1 and 4.4 (Nsight Systems reports in `gpu_profiling/f11/`) |
| `work/c9_gpu_tail/` | C9 tail study (section 4.6) |
| `work/c3_spec_check/` | CPU checks of C3 3d MESH |
| `work/gpu_opt/srcGeospace` | experimental copy, same engine as the merged code |
| `work/simd_eval/` | `REPORT.md` (AVX), `cache_test/RESULTS.md` (stencil-cell cache, CPU and GPU) |

### 6.3 Commands

From `/data/vtenishe/GEO/AMPS`.  Since 2026-10-27 srcGeospace's end-to-end tests are
srcEarth's tests, run with srcEarth's runner through `test/tools/run_srcearth_tests.sh`
(`--np/--nt` rewrite the commands' `-np/-nt`; C13 and C18 keep theirs).  The passes of
2026-10-26/27 used srcGeospace's former copy of the tests and runner.

```
RUNS=/data/vtenishe/GEO/AMPS_runs/srcGeospace
R=$RUNS/results/<date>_gpu

# clean builds of the current code (each call wipes its stage first)
srcGeospace/scripts/build_amps.sh --stage $RUNS/builds/gpu_double_<date> --gpu -j 12
srcGeospace/scripts/build_amps.sh --stage $RUNS/builds/gpu_float_<date> --gpu --float -j 12
srcGeospace/scripts/build_amps.sh --stage $RUNS/builds/cpu_double_<date> -j 12   # CPU reference

T=srcGeospace/test/tools

# one GPU (device 0), one rank with 32 host threads; srcEarth's list without SA01-SA03
GEOSPACE_AMPS=$RUNS/builds/gpu_double_<date>/amps GEOSPACE_GPUS=0 \
$T/run_srcearth_tests.sh --run-root $R/gpu1_double --amps srcGeospace/scripts/amps_gpu_select.sh \
    --exclude 'test/SA0[123]/' --np 1 --nt 32 -- -j 1

# two GPUs, one rank per GPU
GEOSPACE_AMPS=$RUNS/builds/gpu_double_<date>/amps GEOSPACE_GPUS=0,1 \
$T/run_srcearth_tests.sh --run-root $R/gpu2_double --amps srcGeospace/scripts/amps_gpu_select.sh \
    --exclude 'test/SA0[123]/' --np 2 --nt 16 -- -j 1

# srcGeospace's own unit test (engine components vs legacy code; host only)
$T/run_srcearth_tests.sh --run-root $R/unit --amps $RUNS/builds/cpu_double_<date>/amps \
    --list srcGeospace/test/list

# compare with a CPU run root of the same precision
python3 srcGeospace/test/tools/compare_outputs.py <cpu root> $R/gpu1_double --out $R/gpu1_double_vs_cpu.csv
srcGeospace/test/tools/compare_builds.py <cpu root> $R/gpu1_double -o $R/report_gpu1_double \
    --label-a CPU-double --label-b GPU-double

# C3 3d MESH isolation (section 5, item 2): one test, per-trajectory check on the GPU
AMPS_EARTH_ENGINE_VERIFY=1 GEOSPACE_AMPS=$RUNS/builds/gpu_double_<date>/amps GEOSPACE_GPUS=0 \
$T/run_srcearth_tests.sh --run-root $R/c3_gpu_verify --amps srcGeospace/scripts/amps_gpu_select.sh \
    --only 'C3/run_C3.py --mode 3d --mode3d-field-eval MESH' --np 1 --nt 32 -- -j 1
grep -r 'earth-engine\]\[verify\]' $R/c3_gpu_verify/test_output/C3_3d_mesh

# C9 tail study: one interval, double then float, on GPU 0 (about 1 h each)
TAIL_GPU=0 $RUNS/work/c9_gpu_tail/run_interval.sh gpu_double $RUNS/builds/c9tail_gpu_double/amps
TAIL_GPU=0 $RUNS/work/c9_gpu_tail/run_interval.sh gpu_float $RUNS/builds/c9tail_gpu_float/amps
```

`run_interval.sh` writes `runs/<name>/amps_ts.log` (time-stamped output with the
`[engine-trace]` lines), `gpu<N>_smi.log` and `run_info.txt`; it starts `amps` through
`results/2026-10-26_parallel/lockrun.sh light`, which does not wait on an idle machine.
If the engine changes, rebuild the copy with
`srcGeospace/scripts/build_amps.sh --app $RUNS/work/c9_gpu_tail/srcGeospace --stage ... --gpu`.

### 6.4 Lessons from running GPU streams here

* Never stop or kill the shell of a Claude Code background task whose children must keep
  running: the task manager then kills the whole process group (2026-10-26 23:50, the
  `gpu1_double` runner died; STATUS.md section 7).  Make the next step a no-op instead.
* Heavy single-rank tests need 75-100 GiB of host memory (section 2.3); when CPU and GPU
  streams share the machine, use the memory-based lock (`lockrun.sh`) so that heavy tests
  do not overlap.  Under that rule GPU 1 sat idle 07:07-12:14 on 2026-10-27 waiting for
  the other GPU's C19.
* Use binaries with the allocator fix for single-rank runs: without it the single-rank
  setup alone took more than 47 min (C2 3d MESH, never finished) and 1.6 h (C19).
