# srcEarth validation tests

This directory contains executable regression/validation tests for the AMPS
Earth SEP/geospace backward products.  Each test is stored in its own directory
(`C1`, `C2`, ..., `F1`, `F2`, ...).  Test scripts are intended to be executed
from the directory containing the `amps` executable, not from inside the test
subdirectory.

Common convention:

```bash
python srcEarth/test/<TEST_ID>/run_<TEST_ID>.py -np 4 -nt 16
```

where `-np` is the number of MPI ranks passed to `mpirun` and `-nt` is the number
of threads per MPI rank.  Defaults are `-np 4` and `-nt 16`.

## C1 — Pure dipole vertical Størmer cutoff

C1 compares the numerical vertical cutoff rigidity with the analytical Størmer
formula in a centered aligned dipole.  It can be run in either Mode3D or gridless
mode:

```bash
python srcEarth/test/C1/run_C1.py --mode 3d --mode3d-field-eval ANALYTIC
python srcEarth/test/C1/run_C1.py --mode 3d --mode3d-field-eval MESH
python srcEarth/test/C1/run_C1.py --mode gridless
```

The Mode3D `--mode3d-field-eval ANALYTIC` option passes
`-mode3d-field-eval ANALYTIC` to AMPS.  It is useful for separating the pusher and
classifier from mesh interpolation error.  `--mode3d-field-eval MESH` exercises
the full mesh-backed Mode3D field path.


### Parser-compatible inputs

C1 input templates follow the CCMC/RoR-style `AMPS_PARAM_test.in` layout that is
known to pass the current parser.  The test harness passes newer validation
controls through the AMPS command line instead of placing them as active input
keywords.  This avoids parser failures in code checkouts where the CLI support is
newer than the input-file parser.

## C2 — Dipole longitude and north/south symmetry

C2 checks that a centered aligned dipole solution is independent of longitude and
symmetric between the northern and southern hemispheres.  It uses parser-compatible
shell syntax rather than validation-plan shorthand:

```text
OUTPUT_MODE            SHELLS
SHELL_COUNT            1
SHELL_ALTS_KM          9000.0
SHELL_RES_DEG          30
```

The Python harness extracts latitudes `-60, -30, 0, +30, +60` and longitudes
`0, 30, ..., 330` from the shell output and reports longitude spread and paired
north/south differences.

```bash
python srcEarth/test/C2/run_C2.py --mode 3d --mode3d-field-eval ANALYTIC
python srcEarth/test/C2/run_C2.py --mode 3d --mode3d-field-eval MESH
python srcEarth/test/C2/run_C2.py --mode gridless
```

C2 also accepts `--max-trace-time` and `--max-trace-distance`; these values are
written into the generated input file so the effect of finite-time or
finite-distance trajectory classification can be tested explicitly.

## C3 — Penumbra, Rc_upper, and UPPER_SCAN regression

C3 targets the high-latitude outer-shell dipole point where a simple endpoint
binary search can collapse to the lower rigidity bound.  It runs BINARY and
UPPER_SCAN by default and applies the PASS/FAIL gate to UPPER_SCAN.

```bash
python srcEarth/test/C3/run_C3.py --mode 3d --mode3d-field-eval ANALYTIC
python srcEarth/test/C3/run_C3.py --mode 3d --mode3d-field-eval MESH
python srcEarth/test/C3/run_C3.py --mode gridless
```

The target is `lat=±60 deg`, `alt=9000 km`, where the analytical vertical
Størmer cutoff is about `0.160 GV`.  Since `CUTOFF_EMIN=1 MeV/n` corresponds to
`Rmin≈0.043331 GV`, a numerical result near `Rmin` indicates lower-boundary
collapse or long-time leakage rather than the true upper cutoff.

C3 accepts `--max-trace-time`, `--max-trace-distance`, and `--dt-trace` so the
finite-time / finite-distance trajectory classification sensitivity can be
studied explicitly.

## C4 — Parser-safe trajectory integration convergence

C4 replaces the earlier non-parser-safe invariant diagnostic.  The current code
does not recognize `CUTOFF_DEBUG_EXIT_TRACE` or the related
`CUTOFF_DEBUG_EXIT_*` keywords, so C4 uses only production parser keywords and
checks ordinary cutoff shell output against the analytical vertical Størmer
solution.

The script runs a centered dipole shell at 9000 km for a sweep of `DT_TRACE`
values and checks selected longitudes/latitudes for convergence toward the
Størmer cutoff.

```bash
python srcEarth/test/C4/run_C4.py --mode 3d --mode3d-field-eval ANALYTIC
python srcEarth/test/C4/run_C4.py --mode 3d --mode3d-field-eval MESH --dt-sweep 1.0,0.5,0.25
python srcEarth/test/C4/run_C4.py --mode gridless --max-trace-distance 300
```

C4 records `MAX_TRACE_TIME` and `MAX_TRACE_DISTANCE`, because finite-time and
finite-distance caps can affect near-cutoff trajectory classification.


### ADAPTIVE_DT time-step control

`ADAPTIVE_DT` is a `#NUMERICAL` switch shared by Mode3D and gridless backward
cutoff/density tracing:

```text
ADAPTIVE_DT T   # default: DT_TRACE is the maximum allowed adaptive step
ADAPTIVE_DT F   # fixed-step regression mode: use DT_TRACE directly
```

The command-line override is:

```bash
-adaptive-dt T|F
```

Aliases `-fixed-dt`, `-no-adaptive-dt`, and `-use-adaptive-dt` are also accepted.
Use `ADAPTIVE_DT F` for pusher convergence tests where changing `DT_TRACE` must
change the actual integration step.  Production runs should normally use the
default adaptive mode.

## C6 — IGRF / realistic global cutoff morphology

C6 computes a global vertical cutoff-rigidity shell map at 500 km using a realistic internal geomagnetic-field geometry.  It verifies broad map morphology rather than a point-by-point analytical solution: low-latitude cutoffs should be high, polar cutoffs should be low, the equator-to-pole contrast should be large, and the map should show longitude structure.  When requested, the same map is compared across `DYNAMIC`, `BLOCK_CYCLIC`, and `STATIC` MPI schedulers.

Default run:

```bash
python srcEarth/test/C6/run_C6.py -np 4 -nt 16
```

Full scheduler-regression run:

```bash
python srcEarth/test/C6/run_C6.py --schedulers DYNAMIC,BLOCK_CYCLIC,STATIC -np 4 -nt 16
```

Fast smoke test:

```bash
python srcEarth/test/C6/run_C6.py --quick -np 2 -nt 8
```

The current standalone Mode3D field dispatcher may not expose pure `FIELD_MODEL IGRF`; therefore the default C6 input uses `FIELD_MODEL T96` with quiet external drivers, exercising the IGRF/Geopack setup and realistic internal-field coordinate transforms without requiring a separate IGRF-only selector.
