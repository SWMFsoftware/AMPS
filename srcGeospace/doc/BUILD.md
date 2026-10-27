# Building srcGeospace

All builds happen in a private AMPS stage created by `scripts/build_amps.sh`, by
default `$GEOSPACE_RUNS/builds/<cpu|gpu|multicore>_<double|float>` with
`GEOSPACE_RUNS=<repository>/../AMPS_runs/srcGeospace`.  The repository's own `build/`,
`amps`, `share/`, `Makefile*` and `.last_input_info` are never written, and nothing is
written into `srcGeospace/`.

## 1. Prerequisites

* The repository AMPS configuration (`Makefile.local`, `Makefile.def`, `Makefile.conf`,
  `.amps.conf`): `mpicxx` (g++ 9.4) for C++, `mpif90` as the linker (`INTERFACE=on`),
  OpenMP on, SPICE linked.  The stage copies or links these files unchanged.  The
  application input comes from `srcGeospace/input/` (`geospace.input` and the deck
  `srcGeospace.input`).
* GPU builds only: NVIDIA HPC SDK 24.5 (`nvc++`, CUDA 12.4), activated by
  `source /data/vtenishe/openacc/env.sh`.  The build scripts source it themselves when
  `nvc++` is not on `PATH` (override with `NVHPC_ENV=<file>` or `--nvhpc-env`).
  Target GPUs: 2 x RTX A6000 (compute capability 8.6).

## 2. CPU build (default)

```
srcGeospace/scripts/build_amps.sh -j 8                   # -> $GEOSPACE_RUNS/builds/cpu_double/amps
srcGeospace/scripts/build_amps.sh --stage DIR -j 8
```

Every recompilation is a clean build: `build_amps.sh` first **removes the whole stage
directory** (the previous AMPS `build/` tree, `amps`, `.last_input_info`, configuration
copies, logs and any GPU library) and recreates it.  For safety it only removes a
directory that is empty or carries the marker file `.srcGeospace_stage` written by an
earlier `build_amps.sh` run, and it refuses the repository root, the repository
`build/` and anything inside `srcGeospace/`.  It then links the read-only parts of the
repository into the stage, copies `share/`, `input/`, `Makefile.local` and the `.conf`
files, adds `srcGeospace/input/geospace.input` and `srcGeospace/input/srcGeospace.input`,
rewrites `PTDIR/MYDIR` in a private `Makefile.def`, replaces `ProjectSourceDirectory`
in the stage copy of the deck by the absolute path of the application (`ampsConfig.pl`
copies it with `cp -r`, which would copy a relative symbolic link as a link) and runs

```
make -j N amps EXTRACOMPILEROPTIONS="<--extra-cxx>"
```

The device pool `gpu/EngineDevice.cpp` is compiled by `mpicxx -fopenmp -O3` like every
other application file (`GEOSPACE_DEVICE_OBJ ?= gpu/EngineDevice.o` in
`makefile`).  Its kernels, and those of the host pool, launch through the
persistent thread pool `gpu/HostThreadPool.o` by default, or through `omp parallel for`
with `AMPS_EARTH_HOST_RUNTIME=OPENMP` (section 7.1).  Backend name at run time:
`CPU-threads/double`, or `CPU-OpenMP/double`.

## 3. Float build

```
srcGeospace/scripts/build_amps.sh -j 8 --float         # -> builds/cpu_float
# equivalent: --extra-cxx -DGEOSPACE_PRECISION_FLOAT
```

`-DGEOSPACE_PRECISION_FLOAT` reaches every AMPS compile through
`EXTRACOMPILEROPTIONS`.  It changes only the engine (`Real=float`, normalized geospace
units; see `PRECISION.md`).  Backend name: `CPU-threads/float` (or `CPU-OpenMP/float`).  The float build is not
bitwise comparable with the legacy tracer, so `AMPS_EARTH_ENGINE_VERIFY=1` reports differences
by design.

## 4. GPU build

### 4.1 One command

```
srcGeospace/scripts/build_amps.sh -j 8 --gpu                          # -> builds/gpu_double
srcGeospace/scripts/build_amps.sh -j 8 --gpu --float                  # -> builds/gpu_float
srcGeospace/scripts/build_amps.sh -j 8 --gpu --gpu-target multicore   # -> builds/multicore_double
```

Steps:

1. `scripts/build_gpu_engine.sh` compiles **only** `gpu/EngineDevice.cpp` with nvc++
   into `<stage>/gpu_engine/libGeospaceEngine.so`:
   ```
   nvc++ -std=c++17 -O2 -fPIC -acc=gpu -gpu=cc86,cuda12.4,nordc,nofma [-DGEOSPACE_PRECISION_FLOAT] -Minfo=accel \
         -c gpu/EngineDevice.cpp -o EngineDevice.o
   nvc++ -shared -fPIC -acc=gpu -gpu=cc86,cuda12.4,nordc,nofma -Wl,-soname,libGeospaceEngine.so \
         -o libGeospaceEngine.so EngineDevice.o
   ```
   (`-acc=multicore -Mnofma` for `--gpu-target multicore`.)  The SDK driver records a
   RUNPATH to its own `compilers/lib` and `cuda/12.4/lib64`, so the library finds
   `libacchost`, `libaccdevice`, `libaccdevaux`, `libcudadevice`, `libnvomp`,
   `libnvhpcatm`, `libnvcpumath` and `libnvc` without `LD_LIBRARY_PATH`.  The script
   checks that at least the 14 original `eg_dev_*` functions are exported (17 since the
   asynchronous `eg_dev_launch/wait/pending_launch` were added) and writes
   `link_flags.txt`.
2. The AMPS build runs with two extra make command-line variables:
   ```
   make -j N amps EXTRACOMPILEROPTIONS="..." GEOSPACE_DEVICE_OBJ= \
        EXTRALINKEROPTIONS="<stage>/gpu_engine/libGeospaceEngine.so -Wl,-rpath,<stage>/gpu_engine -Wl,-rpath,<SDK>/compilers/lib"
   ```
   * `GEOSPACE_DEVICE_OBJ=` (empty; command-line variables propagate to the sub-make
     in `build/main` through `MAKEFLAGS`, and `?=` keeps the empty value) drops the g++
     `gpu/EngineDevice.o` from `mainlib.a`.  Every other engine object stays g++ and is
     linked as usual.  `gpu/HostThreadPool.o` is still needed by `EngineHost.o` (host
     pool), but the nvc++ `EngineDevice.cpp` does not use it: under `_OPENACC` it
     includes no `HostThreadPool.h`.
   * `EXTRALINKEROPTIONS` is assigned empty in the repository `Makefile` and is used
     only at the end of the `amps_link` rule, after `libAMPS.a`.  A command-line value
     overrides the in-file assignment, so the library is linked **without editing any
     repository file**.  No other makefile in the tree uses this variable.  The
     resulting link line is
     ```
     mpif90 -o amps build/main/main.a build/libAMPS.a -lstdc++ -lmpi_cxx -fopenmp -lmpi++ \
            /home/vtenishe/SPICE/cspice/lib/cspice.a \
            <stage>/gpu_engine/libGeospaceEngine.so -Wl,-rpath,<stage>/gpu_engine -Wl,-rpath,<SDK>/compilers/lib
     ```
3. The script checks that `amps` does not define `eg_dev_create` itself (no g++ device
   pool linked) and prints its RUNPATH.

With `--float` the same define goes to both sides: the nvc++ library and every g++
compile.  The C API carries only SI-double POD types, so a mismatch would still link.
It would make the host pool and the synchronous `TraceSync()` path integrate in a
different precision from the device pool.  `eg_dev_backend_name()` reports the
library's precision (`OpenACC/double`, `OpenACC/float`).

### 4.2 Library only

```
srcGeospace/scripts/build_gpu_engine.sh [--app DIR] [--out DIR] [--precision double|float]
     [--target gpu|multicore] [--arch ccNN] [--cuda X.Y] [--rdc] [--fma] [-O N]
     [--nvhpc-env FILE] [--extra-nvcxx "FLAGS"]
```

Default output: `$GEOSPACE_RUNS/builds/gpu_engine/<target>-<precision>/`.  Any existing
stage can be relinked by hand:

```
cd <stage> && make amps_link EXTRALINKEROPTIONS="$(cat <out>/link_flags.txt)"
```

(only valid when the stage was built with `GEOSPACE_DEVICE_OBJ=`).  Replacing the
`.so` in place with the same path and SONAME needs no relink.

### 4.3 Options and why they are the defaults

* `-gpu=nordc`: `EngineDevice.cpp` is the only translation unit with device code and
  every device function is an inline template in a header, so no device link step is
  needed.  `--rdc` also works (tested): nvc++ then runs `nvlink` when it links the
  shared object.
* `-gpu=nofma` / `-Mnofma`: nvc++ 24.5 contracts `a*b+c` into FMA at -O2 on both
  targets, although `nvc++ -help` lists FMA as "default at -O3".  Measured on the
  kernel SASS (`cuobjdump -sass`): 4867 DFMA with default flags, 3363 with `nofma` (the
  remaining DFMA belong to division, sqrt and libdevice routines).  The multicore object
  has 544 `vfmadd` with default flags and 0 with `-Mnofma`.  The g++ AMPS build never
  contracts, because it targets baseline x86-64 without `-march`.  `nofma` keeps the
  double GPU build comparable with the legacy tracer (`DESIGN.md` §6); use `--fma` for production
  speed once validated.
* `-O2`: -O0 to -O3 all compile (section 6).

## 5. Known build issues

### 5.1 `fmod` has no device implementation in nvc++ 24.5 (fixed: `efmod`)

`EngineTrap.h` (`CurrentProfileBin`) called `std::fmod`.  Before the fix:

```
-gpu=nordc (compile):   ptxas fatal   : Unresolved extern function 'fmod'
                        NVC++-F-0155-Compiler failed to translate accelerator region (see -Minfo
                        messages): Device compiler exited with error status code
                        (.../gpu/EngineDevice.cpp: 224)
-gpu=rdc  (-c passes, link fails):
                        nvlink error   : Undefined reference to 'fmod' in '<object>.o'
                        pgacclnk: child process exit status 2: .../compilers/bin/tools/nvdd
```

So `nvc++ -acc=gpu -c EngineDevice.cpp` "compiles" with the default rdc, but no
GPU library or executable can be linked.  `std::fmod`, `::fmod`, `__builtin_fmod` and
`fmodf` all fail; `std::remainder`, `remquo`, `trunc` and `modf` work
(`tools/device_math/fmod_probe.cpp`).

Fix: [`patches/0001-device-fmod.diff`](patches/0001-device-fmod.diff) adds
`Geospace::GPU::efmod()` to `EngineConfig.h` and uses it in `EngineTrap.h`.  Under
`_OPENACC` it computes `r = remainder(x,y)` and corrects the sign
(`x >= 0 && r < 0 -> r + |y|`, `x < 0 && r > 0 -> r - |y|`).  Both operations are exact,
so the result is bitwise identical to `std::fmod`
(`tools/device_math/efmod_exactness.cpp`: 0 mismatches in 4.1e7 double and 4.1e7 float
cases, including edge cases).  g++ builds keep `std::fmod` unchanged.  The patch is
applied in `gpu/`.
`-acc=multicore` builds compile without the patch (no device code), but use the
same `efmod` path once it is applied.

### 5.2 OpenMP runtime of the GPU executable

The OpenACC host runtime `libacchost.so` references `__kmpc_*` and `__nvomp_*` and
therefore needs `libnvomp.so`.  `-nomp` removes the OpenACC runtime libraries
altogether, so it is not an option.  `libnvomp` exports the full `omp_*` (466) and
`GOMP_*` (268) ABI.  The g++-compiled AMPS objects record `NEEDED libgomp.so.1`, but
the `amps` RUNPATH contains the HPC SDK library directory, and in the SDK
`libgomp.so.1 -> libgomp.so -> libnvomp.so`.  `ldd $GEOSPACE_RUNS/builds/gpu_double/amps`
therefore shows `libgomp.so.1 => .../compilers/lib/libgomp.so.1`: the dynamic loader
maps the same file once (identical inode), so the process has a **single** OpenMP
runtime, NVIDIA's GOMP-compatible `libnvomp`, used by both AMPS's OpenMP regions and the
OpenACC runtime.  This is NVIDIA's intended way of running GCC-compiled OpenMP code in
an nvc++-linked program and avoids two OpenMP runtimes in one process.  Every GPU run of
2026-10-22 to 10-27 used this arrangement; no failure was traced to it
(`../GPU_STATUS.md` §2.2).  (Removing the SDK
directory from the `amps` RUNPATH would resolve `libgomp.so.1` to the system GNU
libgomp and load both runtimes; the engine library has its own RUNPATH for its
dependencies.)

### 5.3 `-acc=multicore` internal compiler error (not reproduced)

An earlier small probe was reported to hit an nvc++ 24.5 internal compiler error with
`-acc=multicore` at -O2.  The probe files that still exist
(`acc_both.cpp`, `off.cpp` in the design-phase scratch directory) compile at
-O1/-O2/-O3 now, and the engine
compiles with `-acc=multicore` at -O0, -O1, -O2 and -O3 in double and float (section 6).
No ICE was observed.

## 6. Recorded build outcomes

### 6.1 srcGeospace, 2026-10-22 (nvc++ 24.5-1, g++ 9.4), clean stages under `$GEOSPACE_RUNS/builds`

| command (`scripts/build_amps.sh -j 8 ...`) | stage | result | wall |
|---|---|---|---|
| (none) | `cpu_double` | OK, `amps` 6.5 MB | 117 s |
| `--float` | `cpu_float` | OK; `Tracer<float,...>` symbols present (host pool) | 117 s |
| `--gpu` | `gpu_double` | OK; `libGeospaceEngine.so` 3.5 MB, 17 `eg_dev_*` exports, `amps` NEEDED `libGeospaceEngine.so`, no g++ `EngineDevice.o` linked | 149 s |
| `--gpu --float` | `gpu_float` | OK; library 3.3 MB, compiled with `-DGEOSPACE_PRECISION_FLOAT` (`gpu_engine/build_info.txt`) | 147 s |

Two builds ran in parallel, so the wall times are upper bounds.  When these builds were
made none of the GPU executables had been run; the GPU runs since are listed in
`../GPU_STATUS.md` §3.1.

### 6.2 Earlier results (srcEarth_gpu prototype, 2026-10-21)

These were measured on the same engine sources under their earlier directory name; the
stages no longer exist.

* `nvc++ -c EngineDevice.cpp` compiles at -O2 and -O3, in double and float, for
  `-acc=gpu` (22-25 s) and for `-acc=multicore` at -O0 to -O3 (2-16 s).  Before the
  `efmod` fix (5.1), `-gpu=nordc` failed with `ptxas fatal : Unresolved extern function
  'fmod'` and `-gpu=rdc` failed at the device link.
* Shared libraries built in every combination: gpu double/float (nordc, nofma: 2.6 /
  2.5 MB), multicore double (`-Mnofma`, 0 `vfmadd`), gpu double with FMA allowed (nordc
  2.5 MB; rdc 3.3 MB with an nvlink step), gpu/multicore float with FMA allowed.
* Complete `amps` builds: `--float`, `--gpu`, `--gpu --float`, `--gpu --gpu-target
  multicore` all succeeded (119-199 s).  Adding the asynchronous launches grew the
  library from 14 to 17 exports (synchronous and `async(1)` kernel instances).

## 7. Running

```
srcGeospace/test/tools/run_srcearth_tests.sh --run-root <run root> --amps <stage>/amps -- -j 6
srcGeospace/scripts/make_run_root.sh <run root> <stage>/amps     # for one runner started by hand
```

See `../test/README.md` and `../test/tools/README.md` (CPU and GPU runs, cross-build
report).

* GPU selection: one device per process; use `ACC_DEVICE_NUM=0|1`.  With several MPI
  ranks per node, set it per rank (for example from `OMPI_COMM_WORLD_LOCAL_RANK` in a
  wrapper).  `NVCOMPILER_ACC_NOTIFY=1` prints every kernel launch, which shows that
  the kernels really run on the GPU.
* The GPU `amps` needs a visible NVIDIA GPU.  On a machine without one, use the CPU
  build.

### 7.1 Environment variables (verified against `gpu/`)

Integer variables are parsed with `atoi`; non-positive or invalid values select the
default.  "host threads" is the worker-thread count the driver passes to
`RunDriverFibers` / `DriverFiberCount` (e.g. the gridless thread count), taken when
the pools are first created.

| variable | default | effect | read in |
|---|---|---|---|
| `AMPS_EARTH_ENGINE` | engine on | `LEGACY` or `legacy`: use the legacy per-trajectory code paths (no fibers, no engine).  Any other value or unset: batched engine | `BatchTracer::Enabled()` (once) |
| `AMPS_EARTH_ENGINE_VERIFY` | off | `1`: also trace every engine trajectory with the legacy tracer and compare bitwise (`DESIGN.md` §6).  At exit prints `[earth-engine][verify] trajectories compared=N bitwise mismatches=M`; prints the first 20 mismatches | `EngineGlue.cpp` (once) |
| `AMPS_EARTH_ENGINE_VERBOSE` | off | `1`: after each driver fiber batch, print workers, host threads, wall time, fibers created, context switches, idle calls, and pool statistics (batched/sync requests, launches, mean active slots, kernel+collect time, backend) | `RunDriverFibers()` |
| `AMPS_EARTH_POOL_CAPACITY` | GPU: 65536; CPU: max(256, 64 x host threads) | device-pool slots; also the default number of driver fibers requested | `BatchTracer::EnsurePools()`, `PreferredFibers()` |
| `AMPS_EARTH_STEPS_PER_LAUNCH` | GPU: 2048; CPU: 256 | K = loop iterations per slot per device-pool launch | `EnsurePools()` |
| `AMPS_EARTH_FIBERS` | CPU: 4096; GPU: pool capacity | maximum number of simultaneously live fibers; when set, also the number of driver fibers requested.  Above 16384, stacks come from 256-stack slabs without guard pages | `Fiber::MaxFibers()` (once), `PreferredFibers()` |
| `AMPS_EARTH_FIBER_GUARD` | off | `1`: guard page for every fiber stack even above 16384 fibers (needs 2 kernel mappings per fiber, `vm.max_map_count`) | `Fiber.cpp` |
| `AMPS_EARTH_DEBUG_INVALID` | off | `1`: Mode3D prints, for every engine trajectory that ends with INVALID_FIELD, which check fired (`EgResult::invalidReason`), the final position and the first unavailable field sample (no AMR leaf inside the configured domain / empty interpolation row) | `CutoffRigidityMode3D.cpp` |
| `AMPS_EARTH_MESH_VERIFY` | off (also on with `AMPS_EARTH_ENGINE_VERIFY=1`) | `1`: after every Mode3D mesh upload compare `EgMesh::Find`/`Interpolate` with AMPS `findTreeNode`/`InterpolateMagneticField` on random and leaf-face points; double build: bitwise; float build: unavailable/extra samples and max relative `|dB|` | `Mode3DMeshExport.cpp` |
| `AMPS_EARTH_MESH_VERIFY_SAMPLES` | 60000 | number of mesh self-check samples | `Mode3DMeshExport.cpp` |
| `AMPS_MODE3D_CUTOFF_FIBERS` | heuristic | Mode3D cutoff batch: number of driver fibers (default min(fair share, max(MPI chunk, 4 x host threads)) on CPU executors with several ranks) | `CutoffRigidityMode3D.cpp` |
| `AMPS_EARTH_FIBER_STACK_KB` | 1024 | fiber stack size in KB (values < 64 ignored); stacks use `MAP_NORESERVE` | `Fiber.cpp` (once) |
| `AMPS_EARTH_TURN_CAPACITY` | 512 | drift-turn diagnostic states per trajectory (with `AMPS_TRAP_TURN_DIAGNOSTICS=1`); overflow truncates and warns once | `BatchTracer` constructor |
| `AMPS_EARTH_HOST_POOL_CAPACITY` | max(64, 16 x host threads) | host-pool slots (Geopack/Tsyganenko fields) | `EnsurePools()` |
| `AMPS_EARTH_HOST_STEPS_PER_LAUNCH` | 64 | K for the host pool | `EnsurePools()` |
| `AMPS_EARTH_ASYNC` | 1 | `0`: synchronous device-pool launches (`eg_dev_advance`).  Otherwise: asynchronous launch/wait, overlapping integration with fiber work (`DESIGN.md` §4).  In the CPU build the async launch runs on `host threads` workers while the main thread continues | `BatchTracer` constructor |
| `AMPS_EARTH_HOST_RUNTIME` | thread pool | `OPENMP` or `openmp`: CPU kernels (g++ device pool and host pool) launch with `#pragma omp parallel for schedule(dynamic,1)`.  Otherwise: persistent `std::thread` pool (`HostParallelFor`).  Results are identical; only the launcher differs | `HostRuntimeIsOpenMP()` (once) |
| `AMPS_EARTH_SPIN_US` | 50 | microseconds a thread-pool worker spins after a launch before blocking on a condition variable (`atoi`; 0 or negative values give practically no spinning) | thread-pool constructor (once) |
| `AMPS_TRAP_TURN_DIAGNOSTICS` | off | the legacy tracer's own switch; `1` records drift-turn states (engine and legacy) | `MakeTrapConfig()` |
| `AMPS_EARTH_SPECULATION` | GPU: on; CPU: off | `1`/`0`: speculative parallel evaluation of the cutoff searches (UPPER_SCAN candidates, bisection midpoints, PENUMBRA_SCAN grid) inside fibers (`gpu/Speculation.h`, `DESIGN.md` §13).  Results are identical to the sequential search; `1` on the CPU reproduces the GPU scheduling for verification | `SpeculationEnabled()` (once) |
| `AMPS_EARTH_SPECULATIVE_DEPTH` | 4 | bisection levels evaluated concurrently per round (2^D-1 midpoints) | `SpeculativeDepth()` (once) |
| `AMPS_EARTH_MESH_CACHE` | on | `0`: no stencil-cell cache in the device pool (uncached mesh interpolation; identical results).  Otherwise the cache of `EngineMesh.h` `BuildStencilCache` is built at every mesh upload (verbose: `[earth-engine] stencil-cell cache: ...`) | `eg_dev_set_mesh()` |
| `AMPS_INTERPOLATION_STENCIL_CACHE` | on | AMPS core switch: `0` disables the core stencil-cell cache and row-structure reuse (`src/pic/pic_interpolation_routines.cpp`); affects the legacy path and the field evaluation outside the engine | core (once) |

Compile-time option: `-DGEOSPACE_FIBER_UCONTEXT` (through `--extra-cxx`) makes the fibers
switch with `swapcontext()` instead of the x86-64 assembly switch of `gpu/Fiber.cpp`
(other architectures always use `swapcontext()`).

Practical notes:

* Fiber budget: on a GPU it defaults to the pool capacity (`AMPS_EARTH_POOL_CAPACITY`,
  65536; `BatchTracer` constructor), on a CPU to 4096; `AMPS_EARTH_FIBERS` overrides
  both.  In practice the test decks keep far fewer trajectories in flight than either
  (one fiber per location or task; `../GPU_STATUS.md` §4.5).  Each live fiber uses two
  memory mappings, so
  `vm.max_map_count` (65530 here) limits the fiber count to about 32000 unless the
  system setting is raised.
* The pool sizes are fixed at the first use in a process.
* Do not combine `AMPS_EARTH_HOST_RUNTIME=OPENMP` with libgomp's default active wait
  when several ranks share the cores.  Idle OpenMP workers spin between the short
  launches and starve the serial fiber phase.  The main session measured gridless F1 at
  28 s with the default `OMP_WAIT_POLICY`, 3.1 s with `OMP_WAIT_POLICY=PASSIVE` and
  2.8 s with the default thread pool (legacy srcEarth executable 2.9-4.1 s), with bitwise
  identical outputs.  This is why the thread pool is the default.
