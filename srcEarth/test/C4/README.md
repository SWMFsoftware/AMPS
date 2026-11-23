# C4 — Mode3D trajectory-exit classifier and invariant diagnostic

C4 is a low-level trajectory-quality test for the Mode3D backward-tracing kernel.  It is not a global cutoff-map validation test.  C1, C2, and C3 compare cutoff values against the analytical Størmer result and exercise the penumbra-aware cutoff search.  C4 instead asks whether individual trajectories are integrated and classified in a physically credible way before the same mover is trusted by cutoff, density, spectrum, or flux products.

## What C4 tests

For each selected start location and rigidity, the test launches the normal AMPS Mode3D vertical backtrace and saves the detailed terminal state.  The Python harness then checks three things:

1. **Trajectory-exit classification.**  A trajectory with rigidity above the analytical vertical Størmer cutoff should escape through the Mode3D outer box and be reported as `allowed=1` with `reason=OUTER_BOX`.  A trajectory below the cutoff should not be counted as allowed; acceptable forbidden outcomes include inner-sphere loss, time limit, step limit, distance limit, or another non-outer-boundary terminal reason.

2. **Rigidity/energy conservation.**  In the C4 input file the field is a static centered dipole and the electric field is off.  The magnetic force should not change kinetic energy, so the diagnostic records `rel_dR`.  The default hard tolerance is `--dR-tol=1e-6`.

3. **Canonical dipole-axis momentum diagnostic.**  For an aligned centered dipole, the canonical angular momentum about the dipole axis should be approximately conserved.  The diagnostic records `rel_dP_axis`.  This is recorded by default and can be made a hard failure with `--fail-on-paxis --paxis-tol=1e-5`.

The purpose is to catch problems such as incorrect terminal-state classification, artificial energy drift, bad boundary crossing logic, or an inconsistent pusher before those issues are hidden inside a shell map or a density/flux integral.

## Single-run many-trajectory diagnostic

Older versions of the C4 harness restarted AMPS once for every trajectory because the debug-exit interface accepted only one `(lon, lat, alt, R)` case at a time.  The current implementation fixes that by using:

```text
CUTOFF_DEBUG_EXIT_TRACE      T
CUTOFF_DEBUG_EXIT_LIST_FILE  c4_debug_trajectories.dat
CUTOFF_DEBUG_EXIT_FILE       cutoff_3d_debug_exit_trace.dat
```

The list file contains all trajectories to trace in one AMPS run:

```text
# lon_deg lat_deg alt_km R_GV label
0.0 -60.0 9000.0 8.0000e-02 low_latm60
0.0 -60.0 9000.0 3.2000e-01 high_latm60
0.0   0.0 9000.0 1.2800e+00 low_lat0
0.0   0.0 9000.0 5.1200e+00 high_lat0
```

The output is one combined Tecplot-style ASCII file, `cutoff_3d_debug_exit_trace.dat`.  AMPS rank 0 writes this diagnostic before the normal Mode3D MPI location scheduler starts.  Therefore the diagnostic output remains a single file even when AMPS is launched with multiple MPI processes and multiple Mode3D worker threads.

## Running the test

Run from the directory containing the `amps` executable:

```bash
srcEarth/test/C4/run_C4.py -np 4 -nt 16
```

Recommended explicit command:

```bash
srcEarth/test/C4/run_C4.py \
  --factors=0.5,2.0 \
  --lats=-60,-30,0,30,60 \
  --lons=0 \
  --alt=9000 \
  -np 4 \
  -nt 16
```

The harness also accepts the form without `=` for negative latitude lists, for example:

```bash
srcEarth/test/C4/run_C4.py --factors 0.5,2.0 --lats -60,-30,0,30,60
```

The script normalizes that command before calling `argparse` so `-60,-30,...` is not mistaken for an option.

## Important options

```text
--lats             comma-separated start latitudes, default -60,-30,0,30,60
--lons             comma-separated start longitudes, default 0
--factors          rigidities as factors times Rc_Stormer, default 0.5,2.0
--dt-sweep         DT_TRACE values; each value is one AMPS run, default 0.25
--movers           mover list; each mover is one AMPS run, default BORIS; HC4 is accepted
--adaptive-dt      T or F, default T
--dR-tol           hard tolerance on |rel_dR|, default 1e-6
--fail-on-paxis    make |rel_dP_axis| a hard pass/fail quantity
--paxis-tol        tolerance used with --fail-on-paxis, default 1e-5
--min-paxis-convergence-order P
                   hard check (requires --adaptive-dt F and two or more distinct
                   --dt-sweep values): on every above-cutoff trajectory and for each
                   consecutive pair of steps, the observed order
                   ln(e_coarse/e_fine)/ln(dt_coarse/dt_fine) of e = |rel_dP_axis|
                   must be >= P; both runs must end on OUTER_BOX.  Writes
                   C4_convergence_order.csv.  Default: off
```


Useful HC4 C4 check:

```bash
srcEarth/test/C4/run_C4.py --movers=HC4,RK4,BORIS --adaptive-dt T
```

HC4 should be interpreted as the high-accuracy full-orbit Boris-family mover: it
uses a fourth-order Yoshida/Forest-Ruth composition of a symmetric
Higuera-Cary/Boris magnetic midpoint step.  In the C4 E=0 dipole diagnostic it
should conserve rigidity like Boris while reducing smooth-orbit phase error toward
RK4-like accuracy.

A command with one mover and one `DT_TRACE` value launches AMPS once and traces all requested trajectories inside that single run.  A mover sweep or a `DT_TRACE` sweep still launches one AMPS run per configuration because the mover and timestep are global run settings.

## Outputs

For each mover/timestep case, the run directory contains:

```text
AMPS_PARAM_C4.in                    rendered AMPS input file
c4_debug_trajectories.dat           many-trajectory input list for AMPS
cutoff_3d_debug_exit_trace.dat      one combined AMPS diagnostic output
C4_expected_cases.csv               expected cases generated by the harness
C4_amps.log                         AMPS stdout/stderr log
```

The top-level C4 work directory is `test_output/C4_exit` by default.
When the script is run with options and without `--workdir`, the directory
name gets a suffix `_<10 hex digits>`, the SHA-1 of the options (`--skip-run`
and `--keep` left out), and the options are recorded in
`workdir_configuration.json` inside it. Each entry of `srcEarth/test/list`
therefore has its own directory, and entries run in parallel by
`test_runner.py -j N` no longer delete each other's files.
It contains:

```text
C4_summary.csv       one row per trajectory with pass/fail checks
C4_case_metrics.csv  one row per mover/timestep configuration
C4_result.json       machine-readable test result
```

## Interpreting failures

- `expected OUTER_BOX but got ...`: the high-rigidity trajectory did not escape as expected.  This can indicate a pusher problem, a boundary configuration issue, too small a maximum trace time/distance, or a field-evaluation problem.

- `expected not OUTER_BOX but got OUTER_BOX`: the low-rigidity trajectory was classified as allowed.  This may indicate a true penumbral allowed island, but for the default factor `0.5` it is usually a useful warning that the chosen case is not a clean below-cutoff trajectory.  Try a smaller factor such as `--factors=0.25,2.0` if this occurs systematically.

- `|rel_dR| > tolerance`: the trajectory is gaining or losing rigidity in a static magnetic field.  Try a fixed-step sweep such as `--adaptive-dt F --dt-sweep=1e-3,5e-4,2.5e-4` to separate fixed-step pusher convergence from classifier effects.  The steps must resolve the gyration: fixed steps of 0.25-1 s are longer than the gyroperiods of the default cases and give misclassifications with every mover.

- `observed order ... < P` (with `--min-paxis-convergence-order`): the invariant drift does not fall as dt^P on an escaping orbit.  A value near 1 for BORIS means a first-order step, such as a field sampled at the wrong point of the step.

- large `rel_dP_axis`: the dipole-axis canonical invariant is drifting.  This is not a hard failure unless `--fail-on-paxis` is set, but it is useful for comparing movers and timestep settings.

## List entries

The five C4 entries of `srcEarth/test/list` (lines 67-77) are all expected to
pass (`P`).

| Line | Command options | Checks |
|---|---|---|
| 67 | `--factors=0.5,2.0 --lats=-60,-30,0,30,60` | BORIS, adaptive steps: classification, `\|rel_dR\| <= 1e-6` |
| 69 | `--adaptive-dt F --dt-sweep=1e-3,5e-4,2.5e-4 --min-paxis-convergence-order 1.9` | BORIS, fixed steps that resolve the gyration: classification, `rel_dR`, and observed second-order convergence of the dipole-axis invariant |
| 71 | `--fail-on-paxis --paxis-tol=1e-5 --movers HC4` | HC4, adaptive steps: classification, `rel_dR`, `\|rel_dP_axis\| <= 1e-5` |
| 74 | `--adaptive-dt T --movers HC4` | HC4, adaptive steps: classification, `rel_dR` |
| 76 | `--adaptive-dt T --movers BORIS` | BORIS, adaptive steps: classification, `rel_dR` |

**Revision of 2026-11-23** (user-approved).  Lines 69 and 71 were expected
failures (`F`), and a runner `F` entry is satisfied by any nonzero exit, so they
enforced nothing:

* Old line 69, `--adaptive-dt F --dt-sweep=1.0,0.5,0.25` (BORIS).  These
  fixed steps are longer than the gyroperiods of the C4 orbits, so no mover can
  resolve them.  The gyroperiods are estimated from the dipole field at the
  start points: about 0.017 s at ±60°, f=0.5, to 0.17 s at 0°, f=2.  In AMPS
  BORIS misclassified 3-5 of 10 orbits.  A standalone reproduction of the
  production formulas also gave 3 of 10 for HC4
  (`test_output/srcEarth/investigation_2026-11-23/c4_fixed_dt.out`).
* Old line 71, `--fail-on-paxis --paxis-tol=1e-5` (BORIS, adaptive steps of
  0.15 rad).  A second-order mover cannot reach 1e-5 there; the corrected BORIS
  gives 5.2e-3.

The replacements keep every tolerance.

The convergence check on line 69 exists because the classification and `rel_dR`
checks did not detect a real defect.  Until 2026-11-23, `BorisStep()` sampled B
at the start of the step instead of its midpoint
(`srcEarth/gridless/GridlessParticleMovers.cpp`).  That scheme is first order: in
AMPS its invariant drift halved when dt halved, and reached 136% on a
below-cutoff orbit, yet it passed lines 67 and 76 and the resolved fixed-step
sweep without the order check.

On the escaping orbits of line 69 (AMPS, fixed steps 1e-3 -> 5e-4 -> 2.5e-4 s):

| Build | `\|rel_dP_axis\|`, worst orbit | observed order |
|---|---|---|
| before the fix | 1.35e-2 -> 6.71e-3 -> 3.35e-3 | 0.98-1.01 (FAIL) |
| corrected BORIS | 8.47e-5 -> 2.12e-5 -> 5.30e-6 | 2.000 on every orbit and pair (PASS) |

The required order 1.9 is the method's order 2 with a 5% margin (chosen by the
user).  Below-cutoff orbits are excluded from this check: they end on
inner-sphere, distance or step limits whose integration span changes with dt.
At 2.5e-4 s the quasi-trapped orbits reach the 300000-step limit, so their
drifts are not comparable between steps.

Line 71 keeps the strict 1e-5 invariant tolerance for the fourth-order HC4
(6.4e-7).  The BORIS invariant is guarded by the convergence order on line 69.

Records: `test_output/srcEarth/investigation_2026-11-23/` and
`test_output/srcEarth/fix_verification_2026-11-23/` (not committed).

## What C4 does not test

C4 does not validate the global cutoff map, longitude symmetry, mesh interpolation convergence, or realistic IGRF/Tsyganenko morphology.  Those are covered by other validation tests.  C4 is intentionally narrower: it verifies that a single backtraced trajectory terminates for the right reason and conserves the quantities it should conserve in the clean dipole, `E=0` limit.
