# F3 — Dipole cutoff-filtered flux

F3 validates the first end-to-end magnetic-access flux calculation in the SEP
validation campaign.  It links the dipole cutoff benchmark to the density/flux
folding path by running gridless density/spectrum in a centered aligned dipole:

```text
FIELD_MODEL               DIPOLE
DS_TRANSMISSION_MODE      SCAN
DS_TRANSMISSION_SCAN_N    100
SPECTRUM_TYPE             POWER_LAW
SPEC_GAMMA                3.5
SPEC_EMIN/SPEC_EMAX       1 / 1000 MeV
```

The validation document describes F3 as a 9000 km shell latitude profile.  The
runner implements this with explicit `POINTS` rather than `SHELLS`, so it can use
the standard `gridless_points_density.dat`, `gridless_points_spectrum.dat`, and
`gridless_points_flux.dat` files.  The default point set is:

```text
alt = 9000 km
lon = 0, 90 deg
lat = -70, -60, -45, -30, 0, 30, 45, 60, 70 deg
```

## Reference solution

For each point, `run_F3.py` computes the vertical Størmer cutoff rigidity

```text
Rc(lambda,r) = 14.9 cos^4(lambda) / r_RE^2 GV
```

then converts that rigidity into a proton kinetic-energy cutoff, `Ecut`.  The
reference flux in each energy channel is

```text
F_ch = 4*pi*T_open*int_{max(E1,Ecut)}^{E2} J0*(E/E0)^(-gamma) dE
```

where `T_open` is the analytic straight-line open-sky fraction outside the inner
absorbing sphere.  This is an approximate external reference: the AMPS result is
computed from a full directional access map and contains penumbra structure,
while the Størmer reference collapses access to a vertical hard cutoff.

The runner also applies tighter internal checks that should be exact apart from
ASCII output precision:

- `J_local(E) = T(E) J_boundary(E)`
- density and flux reconstructed from the saved spectrum match the AMPS files
- centered-dipole longitude symmetry
- centered-dipole north/south symmetry
- total flux increases toward high absolute latitude

## Running

From the directory containing the `amps` executable:

```bash
python srcEarth/test/F3/run_F3.py -np 4
python srcEarth/test/F3/run_F3.py --fast -np 4
python srcEarth/test/F3/run_F3.py --mover BORIS --lons 0 --lats 45 --scan-n 40 --directions-per-energy 1152 -np 18
python srcEarth/test/F3/run_F3.py --mover RK4 --lons 0 --lats 45 --scan-n 40 --directions-per-energy 1152 -np 18
python srcEarth/test/F3/run_F3.py --lons 0,90,180,270
python srcEarth/test/F3/run_F3.py --analytic-flux-tol 0.75 --rc-tol 0.75
python srcEarth/test/F3/run_F3.py --runner-heartbeat-sec 15 --task-heartbeat-sec 15 --slow-task-sec 30
python srcEarth/test/F3/run_F3.py --dry-run --workdir test_output/F3_dryrun
python srcEarth/test/F3/run_F3.py --skip-run --workdir test_output/F3_gridless
```

## Fast regression profile

Use `--fast` for routine regression testing when the full directional-access
calculation is too expensive:

```bash
python srcEarth/test/F3/run_F3.py --fast -np 4
```

The fast profile reduces the dominant traced-trajectory workload by exactly a
factor of sixty:

```text
profile   longitudes   latitudes                         scan N   directions/E   trajectories
full      0,90         -70,-60,-45,-30,0,30,45,60,70     100      1152           2,073,600
fast      0,90         -70,-45,0,45,70                     12       288              34,560
```

The five fast-profile latitudes still provide:

- an equatorial point with the strongest shielding;
- an intermediate-cutoff north/south pair at ±45 degrees;
- a high-latitude north/south pair at ±70 degrees;
- two longitudes for the centered-dipole longitude-symmetry check;
- three absolute-latitude levels for the expected access trend.

The fast profile also uses the existing deterministic `DS_MAX_PARTICLES`
mechanism to retain 288 directions at each rigidity.  The native sky grid is
ordered as 24 zenith levels by 48 azimuths.  Selecting 288 entries therefore
keeps all 24 zenith levels and every fourth azimuth, giving 12 uniformly spaced
azimuths on every zenith ring.  This remains a full-sky directional test rather
than a vertical-only or small-direction shortcut.

The 12 log-rigidity samples and 288-direction sky are deliberately coarse
compared with the full profile, but they remain sufficient for the approximate
Størmer comparison and for detecting major regressions in cutoff filtering,
channel folding, symmetry, or latitude dependence.  Use the full profile for
release validation, angular/energy convergence studies, or investigating a
marginal fast-profile failure.

Explicit `--scan-n`, `--directions-per-energy`, `--lats`, or `--lons` values
override the corresponding `--fast` defaults.  The runner converts the angular
request to `DS_MAX_PARTICLES = scan_n * directions_per_energy`, prints the
resolved profile, trajectory count, and effective speedup before launching
AMPS, and stores the resolved values in the generated input comments and
`F3_result.json`.

F3 gridless tracing is intentionally MPI-only.  The generated input sets
`MODE3D_PARALLEL SERIAL` and one computational thread per MPI rank.  This avoids
concurrent calls into Geopack/Tsyganenko implementations that may use
thread-unsafe Fortran COMMON-block state.  The legacy `-nt` runner option is
accepted for command-line compatibility, but values other than one are ignored
with an explicit diagnostic.  Increase `-np` to allocate more parallel workers.

## Particle mover and adaptive time step

Use `--mover` to select the shared gridless backtracing integrator explicitly:

```text
BORIS  HC4  RK2  RK4  RK6  GC2  GC4  GC6  HYBRID
```

For example:

```bash
python srcEarth/test/F3/run_F3.py --mover BORIS --fast -np 18
python srcEarth/test/F3/run_F3.py --mover RK4 --lons 0 --lats 45 \
  --scan-n 40 --directions-per-energy 1152 -np 18
```

When `--mover` is omitted, the runner does not pass `-mover` and preserves the
AMPS executable's existing default-selection behavior.  The selected value is
printed before launch, written to the generated `AMPS_PARAM_F3.in` comments, and
stored in `F3_result.json`.

The adaptive selector treats `DT_TRACE`, the gyro-angle limit, and the remaining
trace time as **upper bounds**.  The former `100 km / v` minimum-displacement floor
incorrectly increased low-energy full-orbit steps and has been removed.  The former
`0.2*distance_to_boundary/v` upper bound has also been removed: it approached a
boundary geometrically and could prevent the endpoint from ever crossing it.

Every accepted mover step is now checked with an exact intersection of the numerical
trajectory chord against the inner absorbing sphere and the six faces of the outer
Cartesian box.  The first event along the chord determines the physical result.  A
valid small time step is never increased, and an invalid/non-positive step returns an
explicit `INVALID_TIME_STEP` termination instead of continuing with an arbitrary
fallback.  Gridless and Mode3D use the same event policy.


## Boundary events and unresolved trajectories

The production tracer reports one of the following explicit outcomes:

```text
OUTER_BOUNDARY_ALLOWED
INNER_BOUNDARY_FORBIDDEN
MAGNETICALLY_TRAPPED_FORBIDDEN
TIME_LIMIT
STEP_LIMIT
DISTANCE_LIMIT
INVALID_TIME_STEP
INVALID_FIELD
NUMERICAL_FAILURE
```

Outer escape, inner-sphere impact, and conservatively detected magnetic trapping
are resolved physical classifications.  Inner impact and trapping are both forbidden.
All numerical-limit outcomes remain unresolved.  The density/transmission solver defines

```text
N_forbidden = N_inner_forbidden + N_trapped
T_resolved = N_allowed / (N_allowed + N_forbidden)
unresolved_fraction = N_unresolved / N_sampled
```

and never counts an unresolved trajectory as physically forbidden.  With
`--retry-unresolved`, each unresolved direction is repeated once with half
`DT_TRACE`, twice the step limit, and twice the time limit.  The final result and
number of retried directions are saved in:

```text
gridless_termination_summary.dat
```

The file contains one row per point and energy with `N_sampled`, `N_retried`,
`N_inner_forbidden`, `N_trapped`, every numerical termination category, `T_resolved`,
`T_all`, and `unresolved_fraction`.  F3 requires
all count identities to close and fails when any point/energy unresolved fraction
exceeds `--unresolved-tol` (default 0.01).

Useful controls are:

```text
--max-steps N
--max-trace-time S
--boundary-refine-tol-m M
--trap-detection | --no-trap-detection
--trap-min-mirror-points N
--trap-min-bounces N
--trap-outer-margin-re R
--trap-radial-growth-tol-re R
--trap-energy-rel-tol F
--trap-parallel-deadband F
--unresolved-tol F
--retry-unresolved
--no-save-termination-summary
```

`--boundary-refine-tol-m` sets the geometric tolerance used by the exact chord
intersection and by the inward field-evaluation offset at an outer crossing.  The
`BOUNDARY_EVENT_MAX_ITER` input is retained as a reserved compatibility control;
the current line/sphere and line/box intersections are analytic and do not require
an iterative solve.


### Static-field trapped-orbit classification

F3 enables conservative trapping detection because its centered dipole is static.  A
trajectory is classified as `MAGNETICALLY_TRAPPED_FORBIDDEN` only after all of the
following are satisfied:

- repeated sign changes of `p_parallel` provide at least the requested mirror points;
- at least the requested number of complete bounce cycles is observed;
- the radial envelope is stable across the last two cycle comparisons;
- the entire observed orbit remains the requested distance inside the outer box; and
- relative momentum-magnitude variation stays below the configured tolerance.

Elapsed time by itself never proves trapping.  The feature is disabled by default in the
core parser and should remain disabled for time-dependent fields unless a separate
validated policy is provided.  The F3 runner enables it by default and provides
`--no-trap-detection` for convergence and diagnostic comparisons.

The lightweight geometry/status unit test can be run without AMPS:

```bash
python srcEarth/test/F3/run_boundary_unit_test.py
```

## Progress and hang diagnostics

F3 is substantially more expensive than F1/F2 because every energy/rigidity
sample requires a directional access calculation.  With the default 18 points,
100 scan samples, and the 24 x 48 direction grid, the run evaluates approximately
2,073,600 trajectories.  The runner prints this resolved workload estimate before
launching AMPS so an unintentionally large point or scan configuration is visible
immediately.

AMPS output is streamed live to the terminal while also being retained in
`F3_amps.log`.  If AMPS produces no output for the configured interval, the runner
prints a heartbeat containing elapsed time, process ID, log size, the last AMPS
line, and the status of the three expected Tecplot output files.

The gridless MPI driver also supports F3-specific per-rank diagnostics controlled
by the runner:

```text
--runner-heartbeat-sec N   Runner message after N seconds without AMPS output.
--task-heartbeat-sec N     Each MPI rank periodically identifies the task it is
                           about to process: point, energy, and direction block.
--slow-task-sec N          Report a completed MPI task whose wall time is at least N.
```

The default values are 30, 30, and 60 seconds, respectively.  Set any value to
zero to disable that diagnostic.  For example, when diagnosing an apparent stall:

```bash
python srcEarth/test/F3/run_F3.py -np 4 \
  --runner-heartbeat-sec 10 --task-heartbeat-sec 10 --slow-task-sec 20
```

A task-heartbeat line is written *before* the corresponding trajectory block is
started.  Therefore, if one block takes an unusually long time, the last line from
that MPI rank identifies the point index, energy index and value, direction-block
index, and exact direction-index range.  A later `slow-task` line confirms the wall
time if the block eventually completes.

## Files

Repository files:

```text
AMPS_PARAM_F3_gridless.in
reference_F3_dipole_cutoff.csv
run_F3.py
run_boundary_unit_test.py
test_trajectory_boundary.cpp
README.md
```

Run-directory outputs:

```text
AMPS_PARAM_F3.in
F3_amps.log
F3_summary.csv
F3_result.json
reference_F3_dipole_cutoff_used.csv
gridless_points_density.dat
gridless_points_spectrum.dat
gridless_points_flux.dat
gridless_termination_summary.dat
```
