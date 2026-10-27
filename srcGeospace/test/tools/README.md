# srcGeospace test tools

Helpers for running srcEarth's tests against srcGeospace executables with srcEarth's
runner `srcEarth/test/test_runner.py`, and for comparing and timing the runs.  Every tool
reads run roots and writes its output where it is told to; nothing is written into the
source tree.  Builds and results live under `$GEOSPACE_RUNS` (default
`<repository>/../AMPS_runs/srcGeospace`), see `../README.md` and `../../doc/BUILD.md`.

| file | purpose |
|---|---|
| `run_srcearth_tests.sh --run-root DIR --amps EXE [--list LIST] [--only/--exclude RE] [--np N --nt M] [-- RUNNER-OPTIONS]` | Runs a list (default `srcEarth/test/list`; `srcGeospace/test/list` for srcGeospace's own unit test) with srcEarth's runner in a run root with links `srcEarth`, `srcGeospace`, `amps` and `tmp/`.  Always passes `--commit-id` and a last-pass results file inside the run root; refuses `--update-last-pass`. |
| `sublist.py LIST [--only RE] [--exclude RE] [--np N --nt M] -o OUT` | Writes a subset of a list file (a command keeps its `! runner:` directives, `for` header and `last pass:` line).  `--np/--nt` replace the values of existing `-np/-nt` tokens (not in C13/C18, `--np-nt-exclude`).  Example: the GPU scope `--exclude 'test/SA0[123]/' --np 1 --nt 32`. |
| `compare_outputs.py A B [--out files.csv]` | Classifies every file under `A/test_output` against `B/test_output`: IDENTICAL (bytes), IDENTICAL_AFTER_PATHS (equal once run-root paths and the application directory name `srcEarth`/`srcGeospace` are replaced), NONDETERMINISTIC (AMPS mesh cache files, which differ between two runs of the same executable), NUMERIC (same tokens; max abs/rel difference), DIFFERENT, MISSING. Plots, logs and JSON are skipped by default. |
| `variable_stats.py A B [--out variables.csv]` | Per-variable statistics of all Tecplot `.dat` outputs: fraction identical as printed, median/p90/p99/max relative difference. |
| `compare_builds.py A B -o DIR --label-a CPU --label-b GPU` | Cross-build report: test status and run time of both roots, the two tables above, and plots (difference maps of the cutoff shell grids, B-vs-A scatter for other outputs) in `DIR/REPORT.md` and `DIR/plots/`. |
| `summarize.py --root L1=R1 --root L2=R2 ...` | Joins the `test_runner_report.json` of several roots into one Markdown table (status, elapsed time, ratio to the first root). |
| `run_timing.sh DIR L1=AMPS1 L2=AMPS2 ...` | Runs `list_timing` for several executables (srcGeospace or srcEarth builds) in two rounds of opposite order (run roots `DIR/timing_<label>_<round>/`). |
| `list_timing` | 12 representative tests of `srcEarth/test/list` for timing. |
| `list_rounding_subset` | The 24 tests of `srcEarth/test/list` whose float outputs differ from double; used with a double build compiled with `-mfma -ffp-contract=fast` as a rounding-level baseline (`../../doc/VALIDATION_CPU_2026-10-21.md` section 4.2). |

## CPU and GPU runs of srcEarth's tests

```
R=$GEOSPACE_RUNS/results/2026-10-30
T=srcGeospace/test/tools
$T/run_srcearth_tests.sh --run-root $R/cpu_double --amps $GEOSPACE_RUNS/builds/cpu_double/amps \
    --exclude 'test/SA0[123]/' -- -j 1
GEOSPACE_AMPS=$GEOSPACE_RUNS/builds/gpu_double/amps GEOSPACE_GPUS=0 \
$T/run_srcearth_tests.sh --run-root $R/gpu_double --amps srcGeospace/scripts/amps_gpu_select.sh \
    --exclude 'test/SA0[123]/' --np 1 --nt 32 -- -j 1
$T/compare_builds.py $R/cpu_double $R/gpu_double -o $R/report_gpu_vs_cpu \
    --label-a CPU-double --label-b GPU-double
```

Each test's own runner writes its usual comparison plots (model vs reference data or
analytic solution) under `<run root>/test_output/`, for the CPU and the GPU executable
alike; `compare_builds.py` adds the comparison of the two executables with each other.

## Equivalence checks of the engine

1. **In-situ, per trajectory**: `AMPS_EARTH_ENGINE_VERIFY=1` makes every engine
   trajectory also run through the legacy per-trajectory tracer and compares
   termination, step count, time, path length, exit state and trap diagnostics bit for
   bit (`[earth-engine][verify] trajectories compared=N bitwise mismatches=M` per rank).
   Bitwise agreement is expected for the double CPU build and for the double GPU build
   made with `nofma` (the default of `scripts/build_gpu_engine.sh`).
2. **Component test** `../UEngineEquivalence/run_test.sh` (in `../list`):
   `gpu/EngineTrap.h` vs `util/TrajectoryTrapDetector.h` (440 configurations) and
   `gpu/EngineMovers.h` vs `gridless/GridlessParticleMovers.cpp` (BORIS, HC4, RK2/4/6,
   GC2/4/6, HYBRID).
3. **End-to-end**: `compare_outputs.py` / `compare_builds.py` on two run roots.
