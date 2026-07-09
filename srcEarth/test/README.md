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

## C4 — Mode3D trajectory-exit classifier and invariant diagnostic

C4 uses the Mode3D `CUTOFF_DEBUG_EXIT_LIST_FILE` diagnostic to trace many selected
vertical trajectories in one AMPS run and write one combined
`cutoff_3d_debug_exit_trace.dat` file.  It is not a shell-map cutoff-accuracy
test.  Instead, it checks whether individual trajectories terminate for the
right reason and conserve the quantities they should conserve in a static
centered dipole with `E=0`.

The harness generates `c4_debug_trajectories.dat` with all requested
`lon_deg lat_deg alt_km R_GV label` cases, launches AMPS once per mover/DT_TRACE
configuration, and then checks:

- high-rigidity cases above the analytical Størmer cutoff exit through `OUTER_BOX`;
- low-rigidity cases below the cutoff are not counted as allowed;
- `rel_dR` stays below the requested tolerance;
- `rel_dP_axis` is recorded and can be made a hard failure with `--fail-on-paxis`.

```bash
srcEarth/test/C4/run_C4.py --factors=0.5,2.0 --lats=-60,-30,0,30,60
srcEarth/test/C4/run_C4.py --adaptive-dt F --dt-sweep=1.0,0.5,0.25
srcEarth/test/C4/run_C4.py --fail-on-paxis --paxis-tol=1e-5
```

A single mover and single `DT_TRACE` value produces one AMPS run regardless of how
many latitudes, longitudes, or rigidity factors are requested.  Mover and
`DT_TRACE` sweeps still launch one run per configuration because those are global
AMPS settings.  The debug-exit file is written by rank 0 before the normal Mode3D
MPI scheduler starts, so it remains one file even with multiple MPI ranks and
threads.


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
