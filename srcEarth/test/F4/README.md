# F4 — Transmission reconstruction consistency

F4 validates the internal consistency of the gridless SEP density/flux output. It
uses a centered aligned dipole, `DS_TRANSMISSION_MODE SCAN`, and a power-law
boundary spectrum at diagnostic points on the 9000 km shell.

This test does **not** compare to an external cutoff model. Instead, it checks
identities that must hold exactly for the stored transmission and spectrum:

```text
J_local(E) = T(E) J_boundary(E)
n          = 4π ∫ J_local(E) / v(E) dE
F          = 4π ∫ J_local(E) dE
F_channel  = 4π ∫_channel J_local(E) dE
```

The default points are at longitude 0 and latitudes
`-60, -30, 0, 30, 60 deg`; this includes low-latitude/high-shielding and
higher-latitude/more-open cases.

## Files

```text
srcEarth/test/F4/run_F4.py
srcEarth/test/F4/AMPS_PARAM_F4_gridless.in
srcEarth/test/F4/reference_F4_reconstruction.csv
```

The runner writes, by default:

```text
test_output/F4_gridless/AMPS_PARAM_F4.in
test_output/F4_gridless/gridless_points_density.dat
test_output/F4_gridless/gridless_points_flux.dat
test_output/F4_gridless/gridless_points_spectrum.dat
test_output/F4_gridless/F4_summary.csv
test_output/F4_gridless/F4_result.json
test_output/F4_gridless/reference_F4_reconstruction_used.csv
```

## Running

From the directory containing the AMPS executable:

```bash
python srcEarth/test/F4/run_F4.py
python srcEarth/test/F4/run_F4.py -np 4 -nt 16
python srcEarth/test/F4/run_F4.py --scan-n 160 --nintervals 240
python srcEarth/test/F4/run_F4.py --lats -60,-30,0,30,60
python srcEarth/test/F4/run_F4.py --dry-run --workdir test_output/F4_dryrun
python srcEarth/test/F4/run_F4.py --skip-run --workdir test_output/F4_gridless
```

The `--lats -60,...` form is supported directly even though the value starts
with a minus sign.

## Pass/fail checks

The expected value is zero for all rows in `F4_summary.csv` because the test
reports residuals/errors, not physical density or flux values. The main checks
are:

1. Saved `J_boundary(E)` matches the input power law.
2. Saved `J_local(E)` equals `T(E)J_boundary(E)`.
3. `T(E)` remains in `[0,1]`.
4. The density file matches direct reconstruction from `gridless_points_spectrum.dat`.
5. The total flux and every requested energy-channel flux match direct
   reconstruction from `gridless_points_spectrum.dat`.

Default tolerances are intentionally tight but allow the precision loss of ASCII
Tecplot-style output.
