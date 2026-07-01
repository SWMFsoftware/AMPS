# C4 — Parser-safe trace-control and mover convergence

C4 validates the trajectory integration path without using any non-implemented
parser keywords.  It can be run as a trace-control sweep, a mover
cross-comparison, or both.  The earlier design used a hypothetical detailed exit diagnostic
(`CUTOFF_DEBUG_EXIT_TRACE` and related `CUTOFF_DEBUG_EXIT_*` keywords).  Those
keywords are **not recognized by the current parser**, so this test was rewritten
as a parser-safe convergence test.

The test runs an ordinary centered-dipole vertical cutoff calculation on a shell
at 9000 km for a sweep of `DT_TRACE` values and, optionally, a sweep of particle
movers selected with `--movers`.  The output cutoff rigidities are
compared with the analytical vertical Størmer cutoff:

```text
Rc = 14.9 cos^4(lambda) / r_RE^2  GV
```

The checked locations are, by default:

```text
altitude:   9000 km
latitudes:  -60, -30, 0, +30, +60 deg
longitudes: 0, 90, 180, 270 deg
```

This is an indirect pusher/regression test.  As the trajectory step is reduced,
the cutoff values should remain stable and close to the analytical Størmer
reference, with no artificial longitude or north/south asymmetry.  Because
`ADAPTIVE_DT=T` can make `DT_TRACE` only an upper bound, mover-to-mover
comparisons are often more informative than a pure `DT_TRACE` sweep.

## Run

Run from the directory containing the `amps` executable:

```bash
python srcEarth/test/C4/run_C4.py -np 4 -nt 16
```

Defaults are `-np 4` and `-nt 16`.  AMPS is launched through `mpirun`.

Useful variants:

```bash
python srcEarth/test/C4/run_C4.py --mode 3d --mode3d-field-eval ANALYTIC
python srcEarth/test/C4/run_C4.py --mode 3d --mode3d-field-eval MESH --dt-sweep 1.0,0.5,0.25
python srcEarth/test/C4/run_C4.py --mode gridless --max-trace-distance 300
python srcEarth/test/C4/run_C4.py --dt-sweep 1.0,0.5,0.25,0.125
python srcEarth/test/C4/run_C4.py --movers BORIS,RK4,RK6 --adaptive-dt F
```

## Files

```text
run_C4.py                           Python harness
AMPS_PARAM_C4.in                    backward-compatible alias for Mode3D template
AMPS_PARAM_C4_mode3d.in             parser-compatible Mode3D input template
AMPS_PARAM_C4_gridless.in           parser-compatible gridless input template
reference_C4_stormer_convergence.csv analytical target values
```

The input template uses the same parser-compatible CCMC/RoR-style layout as the
working C1 input.  The Python harness edits only parser-supported keywords:

```text
FIELD_EVAL_METHOD
CUTOFF_MAX_TRAJ_TIME
SHELL_ALTS_KM
SHELL_RES_DEG
DT_TRACE
MAX_TRACE_TIME
MAX_TRACE_DISTANCE
```

Mover, scheduler, and thread controls are passed on the AMPS command line:

```text
-mover BORIS|RK2|RK4|RK6|GC2|GC4|GC6|HYBRID
-mode3d-field-eval ANALYTIC|MESH|GRID_3D
-mode3d-parallel THREADS
-mode3d-threads <nt>
-mode3d-mpi-scheduler DYNAMIC|BLOCK_CYCLIC|STATIC
-mode3d-mpi-dynamic-chunk <N>
-gridless-mpi-scheduler DYNAMIC|BLOCK_CYCLIC|STATIC
-gridless-mpi-dynamic-chunk <N>
-density-parallel THREADS
-density-threads <nt>
```

## Outputs

The default run directory is `test_output/C4_3d_analytic`.  Summary files are
written at the root of the run directory:

```text
reference_C4_stormer_convergence.csv
C4_summary.csv
C4_case_metrics.csv
C4_result.json
C4_convergence.png        # written when matplotlib is available
```

Each `(mover, DT_TRACE)` case has its own subdirectory containing the generated
input, AMPS log, and cutoff shell output.

## Notes

`MAX_TRACE_TIME` and `MAX_TRACE_DISTANCE` are part of the validation state.  If
the caps are too loose, low-rigidity quasi-trapped trajectories can leak to the
finite outer box and be misclassified as allowed.  If the caps are too tight,
truly allowed near-cutoff trajectories can be stopped too early.  C4 exposes
both controls so the convergence of the cutoff result can be tested explicitly.


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


### Run-banner diagnostics

Mode3D and gridless cutoff runs print the trace-control state in the AMPS banner:

```text
Particle mover
ADAPTIVE_DT
DT_TRACE and whether it is a fixed step or maximum allowed step
effective dt rule
MAX_TRACE_TIME
CUTOFF_MAX_TRAJ_TIME
effective cutoff trace-time cap
MAX_TRACE_DISTANCE
MAX_STEPS
adaptive limiter constants when ADAPTIVE_DT=T
```

This is important for interpreting C4.  If `ADAPTIVE_DT=T`, decreasing
`DT_TRACE` may not change the actual pusher step because the gyro-angle or
boundary-distance limiter can already be smaller.  In that case, the useful
convergence check is the mover/cap sensitivity, not the nominal `DT_TRACE`
sequence alone.
