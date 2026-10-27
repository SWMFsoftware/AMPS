# srcGeospace tests

srcEarth is the main application and keeps the test suite; srcGeospace keeps only what
srcEarth does not have.

* **End-to-end tests of srcGeospace are srcEarth's tests**: `srcEarth/test/list` and the
  per-test runners next to it, run against a srcGeospace executable with srcEarth's runner
  `srcEarth/test/test_runner.py`.  `tools/run_srcearth_tests.sh` sets up the run root and
  starts that runner.
* **srcGeospace's own unit tests**, which srcEarth does not have, are listed in `list`:
  `UEngineEquivalence/` compares the batched engine's trap detector and particle movers
  bit for bit with the legacy implementations they replace (`gpu/EngineTrap.h` vs
  `util/TrajectoryTrapDetector.h`, `gpu/EngineMovers.h` vs
  `gridless/GridlessParticleMovers.cpp`).  It compiles srcGeospace's sources into `$TMPDIR`.
* **Tools** for runs, comparisons and timing: `tools/` (see `tools/README.md`).

Builds and results stay outside the source tree, in `$GEOSPACE_RUNS` (default
`<repository>/../AMPS_runs/srcGeospace`); see `../BUILD.md`.

## Running

From the repository root:

```
RUNS=$GEOSPACE_RUNS; R=$RUNS/results/<run>

# srcGeospace's own unit test
srcGeospace/test/tools/run_srcearth_tests.sh --run-root $R/unit \
    --amps $RUNS/builds/cpu_double/amps --list srcGeospace/test/list

# srcEarth's tests against a CPU build, six at a time
srcGeospace/test/tools/run_srcearth_tests.sh --run-root $R/cpu_double \
    --amps $RUNS/builds/cpu_double/amps -- -j 6

# one GPU: one rank with 32 threads; SA01-SA03 excluded as in the GPU passes so far
GEOSPACE_AMPS=$RUNS/builds/gpu_double/amps GEOSPACE_GPUS=0 \
srcGeospace/test/tools/run_srcearth_tests.sh --run-root $R/gpu1_double \
    --amps srcGeospace/scripts/amps_gpu_select.sh --exclude 'test/SA0[123]/' \
    --np 1 --nt 32 -- -j 1

# compare two run roots (e.g. GPU against CPU)
srcGeospace/test/tools/compare_builds.py $R/cpu_double $R/gpu1_double -o $R/report \
    --label-a CPU --label-b GPU
```

The run root holds symbolic links `srcEarth`, `srcGeospace` and `amps` and the directory
`tmp/` (`TMPDIR` of every test); tests write `test_output/`, logs and reports there.
`--np/--nt` replace the values of the `-np/-nt` options the commands already have (C13 and
C18 keep their own counts).  Options after `--` go to srcEarth's runner; the wrapper
always passes `--commit-id` and a `--last-pass-results-file` inside the run root, and it
refuses `--update-last-pass`, which would rewrite `srcEarth/test/list`.  The runner's
closing hint to apply the saved results with `--commit-last-pass` is not meant for
srcGeospace runs of srcEarth's list: its last-pass lines record srcEarth's own passes.

A single test runner can be started by hand in a run root made by
`../scripts/make_run_root.sh`, e.g.
`TMPDIR=$PWD/tmp PYTHONDONTWRITEBYTECODE=1 python3 srcEarth/test/C1/run_C1.py --mode gridless -np 4 -nt 16`.

## What srcEarth's tests check for srcGeospace

The end-to-end tests run srcGeospace's executable, so they check its results.  srcEarth's
unit tests and source-text contract tests compile or read srcEarth's sources, not
srcGeospace's; srcGeospace's own code is checked by the end-to-end tests,
`UEngineEquivalence` and the engine's in-run check `AMPS_EARTH_ENGINE_VERIFY=1` (every
engine trajectory also traced by the legacy code and compared bit for bit).
