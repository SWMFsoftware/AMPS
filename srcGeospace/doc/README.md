# srcGeospace: batched trajectory engine documentation

srcGeospace is the AMPS geospace application (cutoff rigidity, directional access,
density/flux/spectrum; gridless and Mode3D; static and time-dependent fields).  It was
created from `srcEarth` at commit `a29ec499` (`../SOURCE_REVISION`) and adds a batched
trajectory engine (`../gpu/`) that runs on CPU threads (persistent thread pool or
OpenMP) or on an NVIDIA GPU (OpenACC, nvc++), from one source.  The original
per-trajectory tracer is kept as the *legacy* path (`AMPS_EARTH_ENGINE=LEGACY`) and is
the bitwise reference of the double build.

| document | contents |
|---|---|
| [DESIGN.md](DESIGN.md) | Architecture: fibers, persistent trajectory pool with K-step kernels and per-trajectory adaptive dt, device and host pools, field routing, exact-port strategy and bitwise verification, Mode3D mesh flattening and the AMPS row-stencil port, fixed-capacity trap-detector arena, GPU data management, limitations, and the status of each part |
| [BUILD.md](BUILD.md) | Building `amps`: CPU (default), float, GPU (nvc++ device library + g++ AMPS link), `-acc=multicore`, run-time environment variables, recorded build outcomes and known build issues |
| [PRECISION.md](PRECISION.md) | Precision and unit policy (compile-time `Real`, double time/path/trap, normalized geospace units for float), the float-vs-double experiment and the in-AMPS float validation |
| [VALIDATION_CPU_2026-10-21.md](VALIDATION_CPU_2026-10-21.md) | CPU validation of 2026-10-21/22 (made while the code was called srcEarth_gpu): srcEarth reference vs double vs float, rounding baseline, timings |
| [INVENTORY.md](INVENTORY.md) | Every source file and header with what builds or uses it (`../scripts/inventory.py`), and the files removed as unused |
| [patches/0001-device-fmod.diff](patches/0001-device-fmod.diff) | Device `fmod` replacement (`efmod`) for nvc++, already applied |
| [../test/README.md](../test/README.md), [../test/tools/README.md](../test/tools/README.md) | srcGeospace's own unit test; running srcEarth's tests against a srcGeospace executable (`run_srcearth_tests.sh`); cross-build comparison and timing |

## Directory layout

The compiled sources have srcEarth's layout (`3d/`, `3d_forward/`, `3d_forward_swmf/`,
`boundary/`, `gridless/`, `util/` and the files at the top).  Not in srcGeospace, because
srcEarth provides them (2026-10-27, user decision): the tests except srcGeospace's own,
the studies, examples, release-validation and standalone-campaign tools, plotting and
post-processing scripts and the inherited application documentation (`srcEarth/README.md`
and the READMEs next to srcEarth's sources).  Additions:

| path | contents |
|---|---|
| `gpu/` | the engine (DESIGN.md section 2 lists every file) |
| `input/` | `geospace.input` and the deck `srcGeospace.input` (`ProjectSourceDirectory=srcGeospace`) |
| `scripts/build_amps.sh` | builds `amps` in a private AMPS stage under `$GEOSPACE_RUNS/builds` (CPU, `--float`, `--gpu`) |
| `scripts/build_gpu_engine.sh` | builds `libGeospaceEngine.so` (device pool) with nvc++ |
| `scripts/make_run_root.sh` | run root (links `srcEarth`, `srcGeospace`, `amps`) for starting one of srcEarth's test runners by hand |
| `scripts/inventory.py` | source inventory (INVENTORY.md) |
| `test/list`, `test/UEngineEquivalence/` | srcGeospace's own unit test: engine components against the legacy code, bit for bit |
| `test/tools/` | `run_srcearth_tests.sh` (srcEarth's tests with srcEarth's runner against a srcGeospace executable), list subsets with `-np/-nt` rewriting, output comparison, cross-build report with plots, timing |
| `tools/precision_experiment/`, `tools/device_math/` | standalone float-vs-double experiment; nvc++ device `fmod` probe |

Builds and test results are kept outside the tree, under `$GEOSPACE_RUNS` (default
`<repository>/../AMPS_runs/srcGeospace`): `builds/<cpu|gpu>_<double|float>/` and
`results/<run>/`.

## Quick start

```
srcGeospace/scripts/build_amps.sh -j 8                     # CPU, double -> $GEOSPACE_RUNS/builds/cpu_double
srcGeospace/scripts/build_amps.sh -j 8 --gpu               # GPU device pool -> builds/gpu_double
srcGeospace/test/tools/run_srcearth_tests.sh --run-root $GEOSPACE_RUNS/results/<run>/cpu_double \
    --amps $GEOSPACE_RUNS/builds/cpu_double/amps -- -j 6      # srcEarth's tests, srcGeospace amps
AMPS_EARTH_ENGINE_VERBOSE=1 ./amps ...                     # pool statistics
AMPS_EARTH_ENGINE_VERIFY=1  ./amps ...                     # bitwise check against the legacy tracer
AMPS_EARTH_ENGINE=LEGACY    ./amps ...                     # legacy per-trajectory paths
```

## Changes relative to srcEarth (a29ec499)

* **Renaming.**  Paths and names `srcEarth/...` became `srcGeospace/...` in documents
  and source comments.  The C++ namespaces
  `Earth` (with all nested namespaces, e.g. `Earth::GPU`, `Earth::Mode3D`),
  `EarthUtil` and `Earth3DForward` became `Geospace`, `GeospaceUtil` and
  `Geospace3DForward`, including the function names in the input-deck macros, the
  source-contract tests and the error messages.  Include guards `_SRC_EARTH_*` became
  `_SRC_GEOSPACE_*`.  The engine's build macros are `GEOSPACE_PRECISION_FLOAT`
  (Real=float build) and `GEOSPACE_DEVICE_OBJ` (makefile), and its device library is
  `libGeospaceEngine.so`.  Of the 84 files srcGeospace shares with srcEarth, 80 differ
  (2026-10-27, after the removal of the tests, studies and non-build copies), 65 of them
  only in these names.
  Unchanged: file names (`Earth.h`, `Earth.cpp`, `Earth.dfn`, `Earth_Sampling.cpp`)
  and the run-time environment variables (`AMPS_EARTH_*`).
* **Engine integration** (lines added/removed):

| file | +/- | change |
|---|---|---|
| `gridless/CutoffRigidityGridless.cpp` | +420/-52 | engine attempt (`EngineTraceTrajectory_`, `TraceTrajectoryAttempt_` with verify), fiber batch scheduling of the task stream, fiber-safe LRU field-evaluator cache with allocation-free keys, cached field kind/fingerprint/context, copy-free request budgets |
| `gridless/DensityGridless.cpp` | +114/-18 | fiber direction loops (results recorded in direction order), energy loops and MPI task streams in batch mode |
| `3d/CutoffRigidityMode3D.cpp/.h` | +768/-6 | Mode3D engine attempt, mesh upload, evaluator accessors, dipole-error merge, nested fiber loops (PENUMBRA_SCAN, companion list), batch scheduling (section 14.9b), INVALID_FIELD diagnostics |
| `3d/DensityMode3D.cpp/.h` | +440/-5 | fiber direction/energy loops and location/task streams in batch mode |
| `3d/GlobalMagneticField.cpp/.h` | +60/-3 | read-only compact-field view and snapshot generation |
| `3d/Mode3D.cpp` | +12/0 | mesh upload after the frozen-field batch starts |
| `makefile` | +14/-1 | engine objects; `GEOSPACE_DEVICE_OBJ` for GPU links |
| `gpu/` (new) | 28 files, 7,499 lines | the engine |

* **Tests** (2026-10-27, user decision).  srcGeospace's end-to-end tests are srcEarth's
  tests, run against a srcGeospace executable with srcEarth's runner
  (`test/tools/run_srcearth_tests.sh`); srcGeospace keeps only its own unit test
  (`test/UEngineEquivalence`, `test/list`) and its tools (`test/tools/`).  srcEarth's
  unit and source-text contract tests check srcEarth's sources, not srcGeospace's.
  Until then srcGeospace had a renamed copy of srcEarth's tests and studies, with a
  runner that had `--run-root/--amps/--np/--nt` (now `run_srcearth_tests.sh` and
  `sublist.py`) and an adapted source-text check of `3d/CutoffRigidityMode3D.cpp` in
  `studies/dec2006_ts05_cutoff_erosion/tests/test_pipeline.py` (both execution paths).
* **Removed** as unused (not compiled, included or referenced): 6 files, listed at the
  end of INVENTORY.md.
* **shieldSim removed** (2026-10-27, user decision): the shielding model `shieldSim/` is
  kept only in `srcEarth/shieldSim`; nothing in srcGeospace used it.
* **Not duplicated** (2026-10-27, user decision: srcGeospace keeps what it needs to
  compile and what srcEarth does not have): `test/` except srcGeospace's own unit test
  and tools, `studies/`, and 56 files that srcEarth has and the build does not use
  (`examples/`, `release_validation/`, `standalone_campaign/`, `gridless/test/`, `plot/`,
  `READ.md`, `TrajectoryFilter.pl`, `tranj.mcr`, `earth-continental-map.dat`,
  `development_and_validation_plan.docx` and the READMEs next to the sources); the
  inherited application documentation was dropped from `../README.md`.  List: end of
  INVENTORY.md.

Every legacy code path is kept and selected with `AMPS_EARTH_ENGINE=LEGACY`.
