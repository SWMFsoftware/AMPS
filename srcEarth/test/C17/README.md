# C17 — Dipole charge-sign and velocity-reversal symmetry

C17 checks the Lorentz-force sign convention independently of observational or
cross-code reference data.  In a static magnetic field with no electric field,
changing the particle charge sign and reversing the launch/arrival direction must
produce the time-reversed access relation

```text
T_q(x0, v, R) = T_-q(x0, -v, R)
```

For the current AMPS gridless cutoff output, the most direct available diagnostic
is the directional cutoff map.  Therefore the runner checks

```text
Rc_q(lon, lat) == Rc_-q(lon + 180 deg, -lat)
```

for every non-polar directional-map cell and every selected observation point.
It also derives a step-transmission proxy at the rigidity list used in the
validation plan:

```text
0.1, 0.2, 0.5, 1, 2, 5, 10 GV
```

A future AMPS output containing the full `T(R,Omega)` table can be compared with
the same mapping directly.

## Files

```text
srcEarth/test/C17/AMPS_PARAM_C17_gridless.in
srcEarth/test/C17/reference_C17_symmetry.csv
srcEarth/test/C17/run_C17.py
srcEarth/test/C17/README.md
```

## Default run

Run from the directory that contains the `amps` executable:

```bash
python srcEarth/test/C17/run_C17.py -np 4 -nt 16
```

The default calculation uses:

```text
FIELD_MODEL            DIPOLE
DIPOLE_TILT            0.0
E_FIELD                off/default
OUTPUT_MODE            POINTS
DIRECTIONAL_MAP        T
DIRMAP_LON_RES         30
DIRMAP_LAT_RES         30
observation lats       -60, -30, 0, 30, 60 deg
observation altitude   9000 km
mover                  BORIS
```

The two physical runs are:

```text
charge_plus:            SPECIES=PROTON,          CHARGE=+1, MASS_AMU=1.0073
charge_minus_reversed:  SPECIES=NEGATIVE_PROTON, CHARGE=-1, MASS_AMU=1.0073
```

The negative-charge run deliberately uses the same mass as the proton run.  This
isolates the charge-sign / velocity-reversal symmetry from kinetic-energy to
rigidity conversion differences that would appear if a true electron mass were
used.

## Useful variants

Faster directional-grid smoke test:

```bash
python srcEarth/test/C17/run_C17.py --dir-lon-res 60 --dir-lat-res 30 -np 2 -nt 8
```

Test several movers:

```bash
python srcEarth/test/C17/run_C17.py --movers BORIS,RK4,RK6 -np 4 -nt 16
```

Use fixed-step integration for a stricter pusher regression:

```bash
python srcEarth/test/C17/run_C17.py --adaptive-dt F --dt-trace 0.25 -np 4 -nt 16
```

Prepare inputs and commands without running AMPS:

```bash
python srcEarth/test/C17/run_C17.py --dry-run
```

Analyze an existing output directory:

```bash
python srcEarth/test/C17/run_C17.py --skip-run --workdir test_output/C17_gridless
```

## Output

The runner writes:

```text
test_output/C17_gridless/reference_C17_symmetry.csv
test_output/C17_gridless/C17_summary.csv
test_output/C17_gridless/C17_pairwise_directional_residuals.csv
test_output/C17_gridless/C17_result.json
test_output/C17_gridless/C17_residual_histogram.png   # if matplotlib is available
```

Each mover has two run directories:

```text
test_output/C17_gridless/<mover>/charge_plus/
test_output/C17_gridless/<mover>/charge_minus_reversed/
```

The AMPS directional-map files are expected to have names such as:

```text
cutoff_gridless_dir_map_point_0000.dat
cutoff_gridless_dir_map_point_0001.dat
...
```

## Pass/fail rule

For each observation point and non-polar directional-map cell,
`run_C17.py` finds the reversed cell in the negative-charge map:

```text
(lon, lat) -> ((lon + 180) mod 360, -lat)
```

The default hard checks are:

```text
abs(Rc_plus - Rc_minus_reversed) <= 1e-8 GV
or
abs(Rc_plus - Rc_minus_reversed)/max(Rc) <= 1e-6
```

and the derived step-transmission values at the selected rigidity list must also
match.  The ±90° directional-map cells are skipped by default because longitude
is degenerate at the poles.  Use `--include-poles` to include them.

## What a failure means

A persistent nonzero residual in C17 usually indicates one of the following:

```text
1. charge sign is not being passed correctly to the pusher;
2. the backtracing velocity/sign convention is inconsistent;
3. the directional-map index or coordinate mapping is wrong;
4. the output direction labels are not using the assumed lon+180, -lat reversal;
5. trajectory time/distance limits are clipping one charge/sign case differently.
```

This test does not prove that the cutoff values are physically correct.  It only
checks an exact internal symmetry.  C1/C12 should be used for analytical Størmer
values, and C8 remains the observational AMS-style East-West validation.
