# C12 — Particle-mover cross-validation against analytical Størmer cutoff

C12 validates the selectable gridless particle movers by running the same pure
centered-dipole vertical cutoff problem with each mover and comparing the result
against the analytical Størmer cutoff,

```text
Rc = R0 cos^4(lambda) / r_RE^2 .
```

This is a **numerical-integrator validation test**.  C1 already checks that the
production/default vertical cutoff path reproduces the analytical dipole result.
C12 asks a more specific question: if the particle mover is changed from BORIS to
RK2/RK4/RK6 or to a guiding-center mover, does the computed cutoff remain
consistent with the same closed-form reference?

## What C12 tests

C12 tests the following:

1. The selected mover can be chosen from the command line with `-mover`.
2. Each full-orbit mover (`BORIS`, `RK2`, `RK4`, `RK6`) reproduces the analytical
   Størmer vertical cutoff on the standard outer shell.
3. Each guiding-center mover (`GC2`, `GC4`, `GC6`) reproduces the analytical
   Størmer vertical cutoff on a high-altitude shell where the guiding-center
   approximation is more appropriate.
4. The high-altitude results from the high-accuracy movers are mutually
   consistent.  This cross-mover comparison is diagnostic by default and can be
   made a hard pass/fail criterion with `--enforce-cross-mover`.

C12 is deliberately a **gridless** test.  The mover dispatch layer used by the
cutoff and density/flux gridless solvers is the intended target of this test.
Mode3D mesh interpolation, AMR field materialization, and mesh-backed magnetic
field errors belong to C5, not C12.

## Reference solution

The reference is analytical:

```text
Rc(lambda, r) = R0 cos^4(lambda) / r_RE^2,
```

where `lambda` is magnetic latitude and `r_RE` is radial distance in Earth radii.
The coefficient used by the runner is the same convention used in C1:

```text
R0 = 0.299792458 * 0.25 * B_eq(Re) * Re,
B_eq(Re) = 3.12e-5 T,
Re = 6371.2 km.
```

The runner writes the reference table to:

```text
reference_C12_stormer_movers.csv
```

and copies it into the run directory.

## Why two shell altitudes are used

The full-orbit movers are checked on the standard C1/C3 outer shell:

```text
full_alt = 9000 km
```

The guiding-center movers are checked on a higher shell by default:

```text
gc_alt = 25000 km
```

The high shell is used because guiding-center motion is an approximation.  It is
not a useful validation target very close to the inner boundary where the particle
motion can be strongly nonadiabatic.  The runner may still record the GC result at
9000 km as a diagnostic, but it does not use that shell as a hard pass/fail gate
for GC movers.

The default `CUTOFF_EMIN` is `0.1 MeV/n`, lower than C1's usual `1 MeV/n`, so the
high-altitude, high-latitude Størmer cutoff remains above the lower scan bound.

## Default command

Run from the directory containing the `amps` executable:

```bash
python srcEarth/test/C12/run_C12.py -np 4 -nt 16
```

This runs:

```text
BORIS, RK2, RK4, RK6, GC2, GC4, GC6
```

with the gridless cutoff solver.

## Useful variants

Run only the full-orbit movers:

```bash
python srcEarth/test/C12/run_C12.py --movers BORIS,RK2,RK4,RK6 -np 4 -nt 16
```

Run only one mover:

```bash
python srcEarth/test/C12/run_C12.py --movers RK6 -np 2 -nt 8
```

Use a fixed time step for pusher convergence investigation:

```bash
python srcEarth/test/C12/run_C12.py --adaptive-dt F --dt-trace 0.25 --movers BORIS,RK4,RK6
```

Make the high-altitude cross-mover comparison a hard pass/fail criterion:

```bash
python srcEarth/test/C12/run_C12.py --enforce-cross-mover
```

Dry-run command generation without launching AMPS:

```bash
python srcEarth/test/C12/run_C12.py --dry-run
```

## Outputs

The default output directory is:

```text
test_output/C12_gridless/
```

Each mover gets its own subdirectory, for example:

```text
test_output/C12_gridless/boris/
test_output/C12_gridless/rk4/
test_output/C12_gridless/gc4/
```

Top-level summary files are:

```text
C12_summary.csv
C12_result.json
reference_C12_stormer_movers.csv
```

Each mover subdirectory also contains:

```text
AMPS_PARAM_C12.in
C12_<mover>_amps.log
C12_<mover>_summary.csv
cutoff_gridless_shells_dipole_compare.dat   # or cutoff_gridless_shells.dat
```

## Acceptance criteria

The default relative-error tolerances are:

| Mover | Default tolerance |
|---|---:|
| BORIS | 1e-3 |
| RK2 | 5e-3 |
| RK4 | 1e-3 |
| RK6 | 5e-4 |
| GC2 | 5e-3 |
| GC4 | 5e-3 |
| GC6 | 5e-3 |

Full-orbit movers are checked at both the standard and high shell.  GC movers are
hard-checked only at the high shell by default.

## Interpreting failures

A failure of only one mover usually indicates a mover-specific integration or
step-control problem.  A failure of all movers in the same direction usually
indicates a shared issue such as cutoff-search settings, boundary classification,
rigidity/energy conversion, or the analytical reference convention.

If high-latitude points fail while mid-latitude points pass, first check whether
`CUTOFF_EMIN` is below the analytical Størmer cutoff at that altitude and
latitude.  If the analytical cutoff is below the lower scan bound, the numerical
cutoff cannot converge to the reference.

If `--adaptive-dt T` gives a failure that disappears with `--adaptive-dt F` and a
smaller `--dt-trace`, the test is identifying a time-step control issue rather
than a field or cutoff-search issue.
