# C1 — Pure dipole vertical Størmer cutoff

This test validates the backward cutoff calculation against the analytical
vertical Størmer cutoff for a centered, aligned dipole field:

```text
Rc = R0 cos^4(lambda) / r_RE^2
```

where `lambda` is magnetic latitude and `r_RE` is geocentric radius in Earth
radii.

The same Python harness can run the test in either standalone Mode3D or gridless
mode:

```bash
python srcEarth/test/C1/run_C1.py --mode 3d -np 4 -nt 16
python srcEarth/test/C1/run_C1.py --mode gridless -np 4 -nt 16
```

For Mode3D, the script exposes the field-evaluation backend:

```bash
python srcEarth/test/C1/run_C1.py --mode 3d --mode3d-field-eval ANALYTIC
python srcEarth/test/C1/run_C1.py --mode 3d --mode3d-field-eval MESH
```

`ANALYTIC` passes the production CLI key:

```bash
-mode3d-field-eval ANALYTIC
```

This bypasses mesh interpolation for analytic fields such as `DIPOLE` and the
Tsyganenko family, and isolates the pusher, cutoff search, and classifier.  The
`MESH` option keeps the Mode3D mesh-backed path and therefore validates field
materialization and interpolation as well.

The script runs AMPS with `mpirun`.  The requested MPI rank and thread counts are
controlled with:

```bash
-np 4
-nt 16
```

Defaults are `-np 4` and `-nt 16`.

The active AMPS input keywords in `AMPS_PARAM_C1_mode3d.in` are intentionally
kept close to the known parser-compatible `AMPS_PARAM_test.in` layout.  Current
validation controls are passed through the AMPS command line by the Python
harness, for example `-cutoff-search UPPER_SCAN`, `-mode3d-field-eval
ANALYTIC`, and the scheduler/thread options.  For older CLI checkouts,
`--no-cutoff-search-cli` suppresses the cutoff-search CLI options while keeping
the parser-compatible input.

The script writes:

```text
C1_amps.log
C1_summary.csv
C1_result.json
C1_stormer_comparison.png      # if matplotlib is available
```

It returns exit code `0` for pass and nonzero for an AMPS execution or validation
failure.
