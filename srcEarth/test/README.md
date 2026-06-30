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
