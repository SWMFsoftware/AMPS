# C14 — Mode3D versus gridless cross-solver consistency

C14 compares the standalone Mode3D and gridless backward cutoff solvers on the
same centered aligned-dipole shell problem.  This is a cross-solver consistency
test: the two solver paths use different field/location representations, but they
should produce the same vertical cutoff map when supplied with the same physical
model and `UPPER_SCAN` cutoff search.

For the dipole case, C14 also compares both solvers with the independent
analytical vertical Størmer cutoff,

```text
Rc = R0 cos^4(lambda) / r_RE^2 .
```

Default grid:

```text
alt = 500, 9000 km
lon = 0, 30, ..., 330 deg
lat = -60, -30, 0, +30, +60 deg
FIELD_MODEL = DIPOLE
CUTOFF_SAMPLING = VERTICAL
```

## Running

Run from the directory containing the `amps` executable:

```bash
python srcEarth/test/C14/run_C14.py -np 4 -nt 16
python srcEarth/test/C14/run_C14.py --mode3d-field-eval MESH -np 4 -nt 16
python srcEarth/test/C14/run_C14.py --scheduler STATIC --dynamic-chunk 0
python srcEarth/test/C14/run_C14.py --lons 0,90,180,270 --lats -60,-30,0,30,60
python srcEarth/test/C14/run_C14.py --dry-run
```

The `--lats -60,...` form is supported directly; `--lats=-60,...` also works.

## What is checked

The runner checks:

1. Mode3D/gridless relative difference at every target point.
2. Mode3D versus analytical Størmer cutoff.
3. Gridless versus analytical Størmer cutoff.
4. Longitude invariance for each solver.
5. North/south symmetry for each solver.

Default cross-solver tolerances are `1%` for `|lat| <= 30 deg` and `5%` at
higher latitudes.  Analytical Størmer tolerances are looser at high latitude
because those points are more sensitive to penumbra and finite trajectory
classification.

## Files

Input templates:

```text
AMPS_PARAM_C14_mode3d.in
AMPS_PARAM_C14_gridless.in
```

Reference/check description:

```text
reference_C14_cross_solver.csv
```

Run artifacts under `test_output/C14_cross_solver`:

```text
mode3d/AMPS_PARAM_C14.in
gridless/AMPS_PARAM_C14.in
mode3d/C14_3d_amps.log
gridless/C14_gridless_amps.log
C14_summary.csv
C14_result.json
C14_mode3d_vs_gridless.png     # written when matplotlib is available
reference_C14_cross_solver_generated.csv
```
