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

## C17 — Dipole charge-sign and velocity-reversal symmetry

C17 checks the exact Lorentz-force time-reversal symmetry in a static centered
dipole with E=0.  It runs a positive-charge and a negative-charge case with the
same mass and compares directional cutoff maps using the reversal relation
`(lon, lat) -> (lon + 180 deg, -lat)`.

```bash
python srcEarth/test/C17/run_C17.py -np 4 -nt 16
python srcEarth/test/C17/run_C17.py --dir-lon-res 60 --dir-lat-res 30 -np 2 -nt 8
python srcEarth/test/C17/run_C17.py --movers BORIS,RK4,RK6
```

The main outputs are `C17_summary.csv`,
`C17_pairwise_directional_residuals.csv`, `C17_result.json`, and
`reference_C17_symmetry.csv` under `test_output/C17_gridless`.

## F1 — Zero-field density/flux normalization

F1 is the first density/flux validation test from `validation.docx`. It isolates
the normalization and output-integration path by using a zero magnetic field and
no inner absorbing boundary:

```text
FIELD_MODEL            NONE
EFIELD_MODEL           NONE
R_INNER                0.0
DS_TRANSMISSION_MODE   DIRECT
SPECTRUM_TYPE          POWER_LAW
SPEC_J0/SPEC_E0/GAMMA  1.0 / 10.0 MeV / 3.5
SPEC_EMIN/SPEC_EMAX    1.0 / 1000.0 MeV
```

In this setup every backtraced direction from every point exits the domain
without magnetic bending or absorption, so the reference solution is
`T(E)=1`, `J_local(E)=J_boundary(E)`, total flux `4π∫J(E)dE`, and density
`4π∫J(E)/v(E)dE`. The test samples ten Cartesian points and requires zero
spatial variation within roundoff.

```bash
python srcEarth/test/F1/run_F1.py -np 4 -nt 16
python srcEarth/test/F1/run_F1.py -np 1 -nt 1
python srcEarth/test/F1/run_F1.py --dry-run
```

F1 writes `F1_summary.csv` and `F1_result.json` under
`test_output/F1_gridless`. The reference values are stored in
`srcEarth/test/F1/reference_F1_zero_field.csv`. The gridless field evaluator now
supports `FIELD_MODEL NONE` by returning `B=(0,0,0)` and skipping
Geopack/Tsyganenko initialization.

