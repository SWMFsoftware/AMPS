# srcEarth validation tests

This directory contains executable regression/validation tests for the AMPS
Earth SEP/geospace backward products.  Each test is stored in its own directory
(`C1`, `C2`, ..., `C6`, `C11`, `C14`, `F1`, `F2`, `F3`, `F4`, `F5`, `F11`, `F12`, `F15`, `F16`, ...).  Test scripts are intended to be executed
from the directory containing the `amps` executable, not from inside the test
subdirectory.

Common convention:

```bash
python srcEarth/test/<TEST_ID>/run_<TEST_ID>.py -np 4 -nt 16
```

where `-np` is the number of MPI ranks passed to `mpirun` and `-nt` is the number
of threads per MPI rank.  Defaults are `-np 4` and `-nt 16`.

## Main test-list runner and parameter loops

`srcEarth/test/test_runner.py` reads the `srcEarth/test/list` file and executes
each listed command as an independent pass/fail regression.  The original list
format remains unchanged:

```text
P srcEarth/test/C3/run_C3.py --algorithms UPPER_SCAN --cutoff-scan-n 200
last pass: 95469cd297ec8494614c0adefee04da69d6488e6
```

`P` means that exit code zero is expected. `F` means that a nonzero exit code is
expected. The optional `last pass:` value records a commit at which that exact
command exited zero. If no `for` declaration is present, parsing, execution,
reporting, and `--update-last-pass` behavior are identical to the historical
runner.

### Looping one command over movers or other CLI values

A `for` declaration expands the immediately following test command:

```text
for $m={RK4,HC4}
{P,P} srcEarth/test/C3/run_C3.py --algorithms UPPER_SCAN --cutoff-scan-n 200 --mover $m
last pass: {95469cd297ec8494614c0adefee04da69d6488e6,95469cd297ec8494614c0adefee04da69d6488e6}
```

This creates two independent scheduled tests:

```text
P srcEarth/test/C3/run_C3.py --algorithms UPPER_SCAN --cutoff-scan-n 200 --mover RK4
P srcEarth/test/C3/run_C3.py --algorithms UPPER_SCAN --cutoff-scan-n 200 --mover HC4
```

Each expansion receives its own log, exit status, expected result, elapsed time,
report row, and `last pass` value. A failure of the RK4 variant does not change
the HC4 result, and vice versa.

The declaration applies to exactly one following source command. A later test is
not part of the loop unless it has its own `for` line.

### Variable syntax

A variable name begins with a letter or underscore and may contain letters,
digits, and underscores. Both command-reference forms are accepted:

```text
$m
${m}
```

Values are comma-separated inside braces:

```text
for $m={BORIS,HC4,RK4,RK6}
```

Whitespace after commas is ignored. CSV-style quoting can protect commas inside
a value, although ordinary mover and CLI values normally do not need quotes.
Every declared variable must be referenced by the command template. This
prevents a misspelled or unused loop variable from silently launching duplicate
commands.

### Expected-result vectors

A loop may provide one expected result per expansion:

```text
for $m={RK4,HC4}
{P,F} command --mover $m
```

Here RK4 is expected to pass and HC4 is expected to fail. The vector length must
exactly match the number of expanded commands.

A scalar marker is shorthand for the same expected result for every expansion:

```text
for $m={RK4,HC4}
P command --mover $m
```

This is equivalent to `{P,P}`.

A vector marker without a preceding `for` declaration is rejected as a list-file
error rather than being interpreted ambiguously.

### Per-variant `last pass` metadata

For a loop, the recommended form is a vector in the same expansion order:

```text
last pass: {commit-for-RK4,commit-for-HC4}
```

An empty element means that no passing commit is currently known for that
variant:

```text
last pass: {known-rk4-commit,}
```

A scalar commit is accepted as shorthand for applying the same known-passing
commit to every expansion:

```text
last pass: shared-commit
```

When `--update-last-pass` is used, looped entries are rewritten in explicit
vector form. Only variants that actually exit zero are updated. Values for
failed variants are preserved. For example, if RK4 passes and HC4 fails:

```text
last pass: {old-rk4,old-hc4}
```

becomes:

```text
last pass: {new-current-commit,old-hc4}
```

This partial-update behavior prevents one mover's successful run from erasing
the provenance of another mover.

### Multiple loop variables

More than one assignment may appear on a `for` line:

```text
for $m={RK4,HC4} $policy={LEGACY,ACCURATE}
P command --mover $m --cutoff-trace-policy $policy
last pass: {c1,c2,c3,c4}
```

The runner forms a Cartesian product in declaration order, with the rightmost
variable changing fastest:

```text
m=RK4, policy=LEGACY
m=RK4, policy=ACCURATE
m=HC4, policy=LEGACY
m=HC4, policy=ACCURATE
```

Therefore the expected-status and `last pass` vectors in this example must each
contain four elements in that order.

### Comments, blank lines, and scope

Blank lines may appear between the `for` declaration, the command, and its
metadata. A comment line beginning with `!` may appear between a loop declaration
and its command. As in the historical parser, a comment after a command closes
the metadata association, so a later `last pass:` line will be rejected rather
than attached across a section boundary.

### Runner examples

Preview expansion without running commands:

```bash
python srcEarth/test/test_runner.py --dry-run --no-memory-gate srcEarth/test/list
```

Run up to six tests concurrently using the normal memory gate:

```bash
python srcEarth/test/test_runner.py -j 6 srcEarth/test/list
```

Run and update scalar or per-variant provenance to the current Git commit:

```bash
python srcEarth/test/test_runner.py -j 6 srcEarth/test/list --update-last-pass
```

The JSON and CSV reports include the expanded command, source command template,
loop-variable bindings, one-based variant index, and total variant count. Text
reports and terminal progress also identify the loop binding, for example
`variant 2/2 [m=HC4]`.

## F3 and legacy-cutoff compatibility architecture

The F3 work exposed an important difference between two uses of the same
backtracing physics:

1. A density/transmission calculation needs a **three-state result**: allowed,
   physically forbidden, or unresolved because a numerical safety limit was
   reached.
2. The historical cutoff search API is Boolean.  It must return either allowed
   or forbidden for every sampled rigidity so that C1, C2, C3, C11, and other
   pre-existing cutoff regressions can continue their search.

Using one interpretation for both callers caused the regression that originally
broke the C-series tests.  F3 correctly introduced explicit `TIME_LIMIT`,
`STEP_LIMIT`, and `DISTANCE_LIMIT` states, but the Boolean cutoff wrappers began
throwing when one of those states remained after a retry.  Low-rigidity trapped
or quasi-trapped trajectories routinely reach a configured cap, so one such
sample could abort an entire shell calculation.

The implementation now separates the **trajectory result** from the
**caller-specific interpretation**.  It does not weaken F3 and it does not add a
special case keyed to a test name.

### Public contracts

| API family | Primary users | Limit termination | Boundary policy | Intended result |
|---|---|---|---|---|
| `TraceTrajectoryShared()` / `TraceTrajectorySharedEx()` | Gridless density, F3, diagnostics | Remains unresolved | Exact first chord event | Full `TrajectoryResult` |
| `TraceTrajectoryMesh()` | Mode3D density/diagnostics | Remains unresolved | Exact first chord event | Full `TrajectoryResult` |
| `TraceAllowedShared()` / `TraceAllowedSharedEx()` | Gridless cutoff and legacy Boolean callers | Returns `false` | Legacy cutoff-compatible | Boolean |
| Internal `TraceAllowedImpl()` | Gridless cutoff search | Returns `false` | Legacy cutoff-compatible | Boolean |
| `TraceAllowedMesh()` / `TraceAllowedMeshEx()` and internal `TraceAllowed3D()` | Mode3D cutoff search | Returns `false` | Legacy cutoff-compatible | Boolean |

The structured APIs are the authoritative interfaces for any calculation that
must distinguish physical loss from numerical incompleteness.  The Boolean APIs
are compatibility interfaces for cutoff searches, where the historical meaning
of “did not escape before the configured cap” is forbidden within the requested
search configuration.

### Explicit termination taxonomy

`srcEarth/util/TrajectoryTermination.h` defines the shared outcomes:

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

The helper functions deliberately keep the categories separate:

```text
IsPhysicalForbiddenTermination = inner impact or validated trapping
IsTraceLimitTermination        = time, step, or distance cap
IsCutoffForbiddenTermination   = physical forbidden or trace limit
IsRetryableNumericalTermination= invalid step, invalid field, or numerical failure
```

`IsResolvedTermination()` was **not** changed to include trace limits.  Therefore
F3 still sees time/step/distance exits as unresolved and excludes them from its
physical transmission denominator.  Only the Boolean cutoff wrappers apply
`IsCutoffForbiddenTermination()`.

### Two integration policies in one tracing implementation

Both gridless and Mode3D tracing now select an explicit internal policy:

```text
StructuredAccurate
LegacyCutoffCompatible
```

The field model, particle mover, geometry, species, and input limits remain
shared.  The policy controls only the numerical details that must differ between
resolved density work and historical cutoff regression behavior.

#### `StructuredAccurate` — F3 and structured density/diagnostic callers

This policy implements the corrected F3 trajectory behavior:

- `DT_TRACE`, the gyro-angle restriction, and the remaining trace time are upper
  bounds.  A valid small step is never enlarged.
- The former `100 km / v` minimum-displacement floor is not applied.  That floor
  could override a smaller gyro-resolution step and change a resolved orbit.
- The former `0.2 * distance_to_boundary / v` limiter is not applied.  Near a
  boundary it could approach the surface geometrically and consume many steps
  without crossing it.
- After every accepted mover step, the numerical chord from the old position to
  the new position is intersected analytically with the inner sphere and all six
  outer-box faces.
- If the chord intersects more than one boundary, the first event along the chord
  determines the result.  This avoids classifying a large step solely from its
  endpoint.
- Allowed exits can retain a refined exit position, momentum, asymptotic
  direction, and field value for anisotropic boundary-spectrum calculations.
- Invalid/non-finite steps and fields are returned as explicit termination
  states rather than hidden by a fallback.
- Time, step, and distance caps remain unresolved.

This is the policy used by F3 transmission sampling and by the structured
Mode3D/gridless trajectory APIs.

#### `LegacyCutoffCompatible` — C-series Boolean cutoff callers

The pre-existing cutoff reference solutions were generated with the earlier
Boolean trajectory behavior.  To avoid changing their sampled penumbra while F3
uses the corrected tracer, cutoff searches retain the old numerical policy:

- The adaptive step includes the gyro-angle restriction.
- It also includes the historical distance-to-nearest-boundary upper limiter,
  `0.2 * distance / speed`.
- It retains the historical `100 km / speed` minimum-displacement floor.  The
  floor is intentionally isolated here; it is not used by F3.
- Inner/outer classification follows the legacy endpoint-oriented cutoff path.
  The mover's existing checked step still protects against an inner-sphere
  crossing, while outer escape is recognized by the cutoff-compatible endpoint
  logic.
- A configured time, step, or cumulative-distance limit returns Boolean
  `false`, exactly as the old cutoff implementation did.

Keeping this policy local to Boolean cutoff APIs is essential.  Reintroducing
its minimum-step floor into the structured tracer would invalidate F3's resolved
trajectory accounting; removing it from cutoff searches changes C1/C2/C3/C11
reference behavior.

### Limit and retry policy

A numerical safety limit is not retried in a Boolean cutoff calculation:

```text
TIME_LIMIT     -> false
STEP_LIMIT     -> false
DISTANCE_LIMIT -> false
```

Retrying such a trajectory merely increases cost and usually repeats the same
classification.  In particular, doubling the time and step limits does not
remove an unchanged `MAX_TRACE_DISTANCE` cap.

Only genuine numerical failures are eligible for one retry:

```text
INVALID_TIME_STEP
INVALID_FIELD
NUMERICAL_FAILURE
```

The retry uses half `DT_TRACE`, up to twice `MAX_STEPS`, and twice the effective
trace-time cap.  If it still ends in a genuine numerical failure, the wrapper
throws with the termination reason and trajectory diagnostics.  This preserves
fail-fast behavior for a broken field or integrator while preventing ordinary
trapped trajectories from aborting a cutoff map.

F3 has its own optional unresolved retry, controlled by the density workflow.
That retry does not change the definition of resolved transmission and is
reported in `gridless_termination_summary.dat`.

### `MAX_TRACE_DISTANCE` semantics

`MAX_TRACE_DISTANCE` is cumulative path length, expressed in Earth radii.  It is
not radial distance from Earth.  A particle can remain at a few Earth radii while
bouncing repeatedly and still accumulate hundreds of Earth radii of traveled
path.

```text
MAX_TRACE_DISTANCE 300.0  -> stop after 300 Re of accumulated trajectory length
MAX_TRACE_DISTANCE 0.0    -> disable only the cumulative-distance cap
```

Setting it to zero does not disable `MAX_TRACE_TIME` or `MAX_STEPS`.  Under the
structured F3 contract, reaching any enabled cap remains unresolved.  Under the
Boolean cutoff contract, reaching an enabled cap returns forbidden.

### Conservative trapped-orbit detection

The new trap detector is retained for F3 and is available to cutoff calculations,
but cutoff correctness does not depend on every trapped orbit being recognized.
The detector requires repeated mirror points, complete bounce cycles, a stable
radial envelope, clearance from the outer boundary, and bounded relative
momentum variation.  When all conditions are satisfied it returns
`MAGNETICALLY_TRAPPED_FORBIDDEN`, a resolved physical outcome.

Elapsed time alone is never considered proof of trapping.  The detector is most
appropriate for static fields such as the centered dipole used by F3.  It is
disabled by default in the core parser and should not be enabled for a
 time-dependent field without separate validation.

### Penumbra-safe and performance-safe `UPPER_SCAN`

The production cutoff search must find the final forbidden-to-allowed transition
below the continuously allowed high-rigidity branch.  A single endpoint binary
search is not sufficient when low-rigidity allowed pockets or forbidden islands
exist.

`UPPER_SCAN` builds the same logarithmic grid as before, but evaluates it from
`Rmax` downward:

1. Evaluate `Rmax`.  If it is forbidden, no allowed upper branch exists in the
   requested bracket and the routine returns the legacy negative sentinel.
2. Move downward until the first forbidden sample is found.
3. The sample immediately above it is known to be allowed, so this pair brackets
   the highest transition.
4. Refine only that bracket by local bisection and return the allowed side.
5. If every sample is allowed, return `Rmin`, preserving the historical finite
   convention.

This is equivalent to evaluating every grid point and selecting the highest
forbidden sample, but it avoids all lower-rigidity trajectories after the upper
transition has been bracketed.  Those lower samples are commonly the slowest,
especially with Mode3D `MESH` field evaluation, because they can remain trapped
until a safety limit.  The descending scan therefore prevents an apparent
0-percent stall in large C2 shell runs without changing the selected upper
cutoff.

### Source files implementing the split

```text
srcEarth/util/TrajectoryTermination.h
srcEarth/gridless/CutoffRigidityGridless.cpp
srcEarth/gridless/CutoffRigidityGridless.h
srcEarth/3d/CutoffRigidityMode3D.cpp
srcEarth/3d/CutoffRigidityMode3D.h
```

The policy is implemented at shared API boundaries rather than in the Python
tests.  Consequently, all cutoff products using the Boolean wrappers receive the
same compatibility behavior, and all structured density/diagnostic products
retain the F3 behavior.

### Validation expectations

The intended regression matrix is:

| Test family | Expected policy exercised |
|---|---|
| C1 | Boolean cutoff compatibility and absolute Størmer agreement |
| C2 | Boolean cutoff compatibility, MESH/ANALYTIC symmetry, finite-cap behavior |
| C3 | Descending `UPPER_SCAN`, upper-cutoff selection, finite-cap compatibility |
| C11 | BINARY-collapse diagnostic plus descending `UPPER_SCAN` recovery |
| F3 | Structured exact-boundary tracing, unresolved accounting, trapping, folding |

For release validation, run both field backends where supported and retain finite
caps in the C-series inputs.  F3 should additionally run its standalone boundary
unit test and verify termination-count closure and the unresolved-fraction gate.

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

## C4 — Parser-safe trajectory integration convergence

C4 replaces the earlier non-parser-safe invariant diagnostic.  The current code
does not recognize `CUTOFF_DEBUG_EXIT_TRACE` or the related
`CUTOFF_DEBUG_EXIT_*` keywords, so C4 uses only production parser keywords and
checks ordinary cutoff shell output against the analytical vertical Størmer
solution.

The script runs a centered dipole shell at 9000 km for a sweep of `DT_TRACE`
values and checks selected longitudes/latitudes for convergence toward the
Størmer cutoff.

```bash
python srcEarth/test/C4/run_C4.py --mode 3d --mode3d-field-eval ANALYTIC
python srcEarth/test/C4/run_C4.py --mode 3d --mode3d-field-eval MESH --dt-sweep 1.0,0.5,0.25
python srcEarth/test/C4/run_C4.py --mode gridless --max-trace-distance 300
```

C4 records `MAX_TRACE_TIME` and `MAX_TRACE_DISTANCE`, because finite-time and
finite-distance caps can affect near-cutoff trajectory classification.

## C6 — global IGRF effective vertical-cutoff validation

C6 is the first global external-reference cutoff test. It evaluates the internal
IGRF field at 20 km geodetic altitude and compares the calculated **effective**
vertical cutoff with three published world-grid families:

```text
INITIAL  Smart--Shea epoch-2000 printed 5 x 30 degree table
COMPLETE FAA CARI-7 tables for 1965, 1980, 1990, 1995, 2000, 2010
MODERN   Gerontidou et al. 2010/2015/2020 grid family
```

The repository contains the full 386,640-row CARI archive and the
machine-readable Smart--Shea and author-supplied Gerontidou-2015 grids. The
modern loader accepts authoritative 2010/2020 CSVs through
`--modern-reference`; it does not substitute CARI values for missing modern
data.

C6 adds and uses `FIELD_MODEL IGRF`, `CUTOFF_SCAN_SPACING LINEAR`, independent
shell longitude/latitude spacing, and `SHELL_GEOMETRY GEODETIC`. The
`PENUMBRA_SCAN` output now includes `Rc_effective_GV`, calculated from the total
allowed rigidity width between `Rc_lower_GV` and `Rc_upper_GV`. Numerical
safety-limit samples leave the effective cutoff unresolved.

```bash
python srcEarth/test/C6/run_C6.py --subtest INITIAL -np 8 -nt 16
python srcEarth/test/C6/run_C6.py --subtest COMPLETE --complete-epochs 2010 --complete-grid-step-deg 5
python srcEarth/test/C6/run_C6.py --subtest COMPLETE --complete-epochs 2000 --full-grid
python srcEarth/test/C6/run_C6.py --subtest MODERN --modern-epochs 2015
python srcEarth/test/C6/run_C6.py --validate-references
python srcEarth/test/C6/run_C6.py --subtest INITIAL --dry-run
```

The complete one-degree CARI reference is bundled, but the default execution
uses a 10-degree subset because all six epochs at full spatial and 0.01-GV
rigidity resolution require hundreds of millions of trajectories. See
`srcEarth/test/C6/README.md` for reference provenance, optional modern files,
comparison tolerances, outputs, and the full-grid workflow.

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

## C11 — Penumbra-pocket regression: BINARY Rmin-collapse defect

C11 targets the code-level cutoff-search defect that motivated `UPPER_SCAN`.  At
the high-latitude dipole-shell point `lon=0 deg, lat=-60 deg, alt=9000 km`, the
legacy endpoint `BINARY` search can return the lower search boundary when an
isolated low-rigidity allowed pocket is present.  The production `UPPER_SCAN`
search should instead return the final upper cutoff, close to the analytical
vertical Størmer value.

```bash
python srcEarth/test/C11/run_C11.py -np 4 -nt 16
python srcEarth/test/C11/run_C11.py --algorithms UPPER_SCAN --cutoff-scan-n 200
python srcEarth/test/C11/run_C11.py --mode3d-field-eval MESH -np 4 -nt 16
python srcEarth/test/C11/run_C11.py --dry-run
```

The test is Mode3D-only because it uses the `CUTOFF_DEBUG_RIGIDITY_SCAN`
diagnostic.  The primary acceptance check requires `BINARY` to reproduce the
known Rmin-collapse signature and requires `UPPER_SCAN` to recover the Størmer
upper cutoff.  The runner also checks repeat shell-output points at `lat=+60 deg`
and at `alt=500 km` as diagnostics.

C11 writes `C11_summary.csv`, `C11_result.json`, and the generated reference
table under `test_output/C11_mode3d`.  Each algorithm subdirectory contains the
rendered input, AMPS log, shell cutoff output, and a debug rigidity scan file.

## C14 — Mode3D/gridless cutoff-band consistency

C14 runs the same centered aligned-dipole vertical-cutoff shell problem through
Mode3D and gridless.  Unlike the historical scalar `UPPER_SCAN`, C14 requests
`PENUMBRA_SCAN`, evaluates one complete rigidity grid, and reports:

```text
Rc_lower
Rc_upper
PenumbraWidth = Rc_upper - Rc_lower
n_allowed_intervals
n_transitions
n_unresolved
```

The analytical vertical Størmer expression is compared with `Rc_lower`; the
upper endpoint is checked directly between Mode3D and gridless.  This avoids the
old false failure produced by comparing the Størmer lower cutoff with the upper
edge of the penumbra.  Time, step, and distance limits remain explicit
unresolved samples and cannot silently define a strict analytical cutoff.

The default `CROSS_SOLVER` profile checks both endpoints, penumbra width,
topology, north/south symmetry, and the 90-degree symmetry of the finite
Cartesian box.  `STRICT_STORMER` additionally requires analytical Mode3D field
evaluation, the accurate cutoff trace policy, and lower-cutoff agreement with
Størmer.

```bash
python srcEarth/test/C14/run_C14.py -np 4 -nt 16
python srcEarth/test/C14/run_C14.py --strict-stormer -np 4 -nt 16
python srcEarth/test/C14/run_C14.py --cutoff-definition LOWER -np 4 -nt 16
python srcEarth/test/C14/run_C14.py \
  --mode3d-field-eval MESH \
  --mode3d-mesh-res-earth-re 0.02 \
  --mode3d-mesh-res-boundary-re 3 \
  --mode3d-mesh-coarsening LINEAR -np 4 -nt 16
python srcEarth/test/C14/run_C14.py --dry-run
```

C14 writes dedicated `cutoff_3d_shells_penumbra.dat` and
`cutoff_gridless_shells_penumbra.dat` files, plus `C14_summary.csv`,
`C14_result.json`, a comparison plot, and a generated reference table.  See
`srcEarth/test/C14/README.md` for the complete algorithm description and a
detailed explanation of every CLI option.

## C17 — Dipole charge-sign and velocity-reversal symmetry

C17 checks the exact Lorentz-force time-reversal symmetry in a static centered
dipole with E=0.  It runs a positive-charge and a negative-charge case with the
same mass and compares directional cutoff maps using the reversal relation
`(lon, lat) -> (lon + 180 deg, -lat)`.

```bash
python srcEarth/test/C17/run_C17.py -np 4 -nt 16
python srcEarth/test/C17/run_C17.py --dir-lon-res 60 --dir-lat-res 30 -np 2 -nt 8
python srcEarth/test/C17/run_C17.py --movers BORIS,RK4,RK6
```

The main outputs are `C17_summary.csv`,
`C17_pairwise_directional_residuals.csv`, `C17_result.json`, and
`reference_C17_symmetry.csv` under `test_output/C17_gridless`.

## F1 — Zero-field density/flux normalization

F1 is the first density/flux validation test from `validation.docx`. It isolates
the normalization and output-integration path by using a zero magnetic field and
no inner absorbing boundary:

```text
FIELD_MODEL            NONE
EFIELD_MODEL           NONE
R_INNER                0.0
DS_TRANSMISSION_MODE   DIRECT
SPECTRUM_TYPE          POWER_LAW
SPEC_J0/SPEC_E0/GAMMA  1.0 / 10.0 MeV / 3.5
SPEC_EMIN/SPEC_EMAX    1.0 / 1000.0 MeV
```

In this setup every backtraced direction from every point exits the domain
without magnetic bending or absorption, so the reference solution is
`T(E)=1`, `J_local(E)=J_boundary(E)`, total flux `4π∫J(E)dE`, and density
`4π∫J(E)/v(E)dE`. The test samples ten Cartesian points and requires zero
spatial variation within roundoff.

```bash
python srcEarth/test/F1/run_F1.py -np 4 -nt 16
python srcEarth/test/F1/run_F1.py -np 1 -nt 1
python srcEarth/test/F1/run_F1.py --dry-run
```

F1 writes `F1_summary.csv` and `F1_result.json` under
`test_output/F1_gridless`. The reference values are stored in
`srcEarth/test/F1/reference_F1_zero_field.csv`. The gridless field evaluator now
supports `FIELD_MODEL NONE` by returning `B=(0,0,0)` and skipping
Geopack/Tsyganenko initialization.

## F2 — Power-law energy integration

F2 uses the same zero-field transport setup as F1 but focuses on the energy folding. It samples one Cartesian point with `T(E)=1` and checks the imposed power-law spectrum against analytical energy integrals.

The validation-plan energy bins `1, 3, 10, 30, 100, 300, 1000 MeV` are encoded as parser-compatible `#ENERGY_CHANNELS`, while `DS_NINTERVALS=960` supplies a fine log-spaced quadrature grid.

```bash
python srcEarth/test/F2/run_F2.py -np 4 -nt 16
python srcEarth/test/F2/run_F2.py --nintervals 1920 --integration-tol 5e-5
python srcEarth/test/F2/run_F2.py --dry-run
```

F2 writes `F2_summary.csv` and `F2_result.json` under `test_output/F2_gridless`. The physical target values are stored in `srcEarth/test/F2/reference_F2_power_law.csv`. The summary file uses `check_type` and `expected_value`, so zero expected values appear only for residual/error metrics such as `T(E)-1` or file-vs-spectrum reconstruction error.

## F3 — Dipole cutoff-filtered flux

F3 is the first magnetic-access flux test. It runs gridless density/spectrum in a centered aligned dipole using `DS_TRANSMISSION_MODE SCAN`, a power-law boundary spectrum, and explicit points on the 9000 km shell. The default points follow the validation-plan latitude profile `-70, -60, -45, -30, 0, 30, 45, 60, 70 deg` at longitudes `0` and `90 deg`.

The external reference is the vertical Størmer hard-cutoff approximation folded analytically with the power-law spectrum:

```text
Rc(lambda,r) = 14.9 cos^4(lambda) / r_RE^2 GV
F_ch = 4*pi*T_open*int_{max(E1,Ecut)}^{E2} J0*(E/E0)^(-gamma) dE
```

where `Ecut` is the proton kinetic energy corresponding to `Rc`, and `T_open` is the analytic straight-line open-sky fraction outside the inner absorbing sphere. The F3 analytical comparison is intentionally looser than F1/F2/F15 because the numerical solver evaluates a full directional access map with penumbra, while the reference collapses access to a vertical step function. The runner also applies tighter internal checks: `J_local(E)=T(E)J_boundary(E)`, file-vs-spectrum flux reconstruction, longitude symmetry, north/south symmetry, and the expected increase of flux toward high absolute latitude.

```bash
python srcEarth/test/F3/run_F3.py -np 4
python srcEarth/test/F3/run_F3.py --fast -np 4
python srcEarth/test/F3/run_F3.py --mover BORIS --lons 0 --lats 45 --scan-n 40 --directions-per-energy 1152 -np 18
python srcEarth/test/F3/run_F3.py --lons 0,90,180,270
python srcEarth/test/F3/run_F3.py --analytic-flux-tol 0.75 --rc-tol 0.75
python srcEarth/test/F3/run_F3.py --dry-run
```

For routine regression testing, `--fast` retains both longitudes and the
representative latitudes `-70,-45,0,45,70 deg`, while reducing the log-rigidity
scan from 100 to 12 samples and using a deterministic 288-direction full-sky
subset.  The subset preserves all 24 zenith levels and 12 evenly spaced azimuths
per level.  This lowers the trajectory count from 2,073,600 to 34,560—exactly
60x—while preserving the equatorial/intermediate/high-latitude access contrast,
north/south symmetry, longitude symmetry, channel folding, and latitude-trend
checks.  The full profile remains the release/accuracy-validation configuration.
The runner uses MPI ranks only and forces one computational thread per rank;
increase `-np` rather than `-nt` for additional parallelism.  `--mover` accepts
`BORIS`, `HC4`, `RK2`, `RK4`, `RK6`, `GC2`, `GC4`, `GC6`, or `HYBRID` and passes
the explicit selection to AMPS.  If it is omitted, the executable default is
preserved.  The shared gridless and Mode3D adaptive selectors no longer impose either the
former `100 km / v` minimum-displacement floor or the asymptotic
`0.2*distance_to_boundary/v` limiter.  Every accepted numerical trajectory chord is
intersected explicitly with the inner sphere and outer box.  The tracer returns
separate allowed, forbidden, time-limit, step-limit, distance-limit, invalid-step,
invalid-field, and numerical-failure states.  Density/transmission excludes unresolved
states from the physical denominator and writes `gridless_termination_summary.dat`.
For the static dipole, F3 also enables conservative trapped-orbit classification:
repeated mirror points, stable bounce envelopes, outer-boundary clearance, and momentum
conservation are all required before a trajectory is counted as trapped/forbidden.  F3
checks count closure and enforces `--unresolved-tol`.  Use `--no-trap-detection` for a
diagnostic comparison, `--retry-unresolved` for one stricter retry, and
`python srcEarth/test/F3/run_boundary_unit_test.py` for the standalone
boundary/termination/trapping unit test.

F3 writes `F3_summary.csv`, `F3_result.json`,
`gridless_termination_summary.dat`, and `reference_F3_dipole_cutoff_used.csv` under
`test_output/F3_gridless`. The default repository reference table is
`srcEarth/test/F3/reference_F3_dipole_cutoff.csv`.

## F4 — Transmission reconstruction consistency

F4 isolates the exact file/spectrum reconstruction identities in the gridless density/flux path. It runs a centered aligned dipole with `DS_TRANSMISSION_MODE SCAN` and power-law boundary flux at the validation-plan diagnostic latitudes `-60, -30, 0, 30, 60 deg` on the 9000 km shell.

The reference solution is internal and exact:

```text
J_local(E) = T(E) J_boundary(E)
n          = 4π ∫ J_local(E) / v(E) dE
F          = 4π ∫ J_local(E) dE
F_channel  = 4π ∫_channel J_local(E) dE
```

The summary file therefore reports residuals with `expected_value=0`; those zero entries mean zero reconstruction error, not zero physical density or flux.

```bash
python srcEarth/test/F4/run_F4.py -np 4 -nt 16
python srcEarth/test/F4/run_F4.py --scan-n 160 --nintervals 240
python srcEarth/test/F4/run_F4.py --lons 0,90 --lats -60,-30,0,30,60
python srcEarth/test/F4/run_F4.py --dry-run
```

F4 writes `F4_summary.csv`, `F4_result.json`, and `reference_F4_reconstruction_used.csv` under `test_output/F4_gridless`. The repository reference/check table is `srcEarth/test/F4/reference_F4_reconstruction.csv`.

## F5 — Directional anisotropy / PAD mapping

F5 validates the nontrivial pitch-angle-distribution mapping in the gridless anisotropic density/spectrum path. It runs three parser-supported PAD cases with the same dipole field, energy grid, and boundary spectrum:

```text
ISOTROPIC
COSALPHA_N, BA_PAD_EXPONENT=2
SINALPHA_N, BA_PAD_EXPONENT=2
```

The primary reference is exact: `sin^2(alpha) + cos^2(alpha) = 1`. Therefore the runner checks the complement identities

```text
T_sin(E) + T_cos(E) = T_iso(E)
J_local_sin(E) + J_local_cos(E) = J_local_iso(E)
n_sin + n_cos = n_iso
F_sin + F_cos = F_iso
```

at three diagnostic locations: 8 Re equator, 8 Re magnetic pole, and 6 Re equator. It also writes and checks a high-energy straight-line semi-analytic PAD reference for `density_aniso/density_isotropic`; those checks are intentionally looser because they neglect finite magnetic bending.

```bash
python srcEarth/test/F5/run_F5.py -np 4 -nt 16
python srcEarth/test/F5/run_F5.py --scan-n 128 --nintervals 96
python srcEarth/test/F5/run_F5.py --semi-n-theta 240 --semi-n-phi 480
python srcEarth/test/F5/run_F5.py --dry-run
```

F5 writes one case directory per PAD model under `test_output/F5_gridless`, plus `F5_summary.csv`, `F5_result.json`, and `reference_F5_pad_mapping_used.csv`. The repository reference table is `srcEarth/test/F5/reference_F5_pad_mapping.csv`.

## F11 — Anisotropic PAD model sum-check

F11 validates exact identities in the gridless anisotropic boundary-spectrum path. It uses `FIELD_MODEL DIPOLE`, `DS_BOUNDARY_MODE ANISOTROPIC`, `DS_TRANSMISSION_MODE SCAN`, and a parser-compatible `#BOUNDARY_ANISOTROPY` section.

The reference solution is pairwise and exact:

```text
BA_PAD_EXPONENT=0:  ISOTROPIC = SINALPHA_N = COSALPHA_N = BIDIRECTIONAL
COSALPHA_N(n):      COSALPHA_N = BIDIRECTIONAL for the same n
```

The runner compares total density, integral-flux channels, `T(E)`, `J_boundary(E)`, and `J_local(E)` for the identity pairs at one low-latitude and one mid/high-latitude point.

```bash
python srcEarth/test/F11/run_F11.py -np 4 -nt 16
python srcEarth/test/F11/run_F11.py --models ISOTROPIC,COSALPHA_N,BIDIRECTIONAL --exponents 0,2,4
python srcEarth/test/F11/run_F11.py --dry-run
```

F11 writes one case directory per PAD model/exponent under `test_output/F11_gridless`, plus `F11_summary.csv`, `F11_result.json`, and `reference_F11_pad_identities_used.csv`. The repository reference table is `srcEarth/test/F11/reference_F11_pad_identities.csv`.

## F12 — Day/night spatial boundary anisotropy identities

F12 validates exact identities in the `DAYSIDE_NIGHTSIDE` spatial boundary model used by the gridless anisotropic density/spectrum path. It uses `FIELD_MODEL NONE`, `R_INNER=0`, `BA_PAD_MODEL=ISOTROPIC`, and diagnostic points on the GSM Y-Z plane (`X=0`) so every sampled trajectory is allowed and the direction grid is paired under `x -> -x`.

The exact references are:

```text
DAYSIDE_NIGHTSIDE(1,1)   = UNIFORM
DAY_ONLY(1,0)            = 0.5 * UNIFORM
NIGHT_ONLY(0,1)          = 0.5 * UNIFORM
DAY_ONLY + NIGHT_ONLY    = UNIFORM
DAYSIDE_NIGHTSIDE(2,0.5) = 1.25 * UNIFORM
```

The runner checks density, total flux, channel fluxes, saved `T(E)`, `J_boundary(E)`, and `J_local(E)`.

```bash
python srcEarth/test/F12/run_F12.py -np 4 -nt 16
python srcEarth/test/F12/run_F12.py --ratio-tol 1e-12 --spectrum-tol 1e-12
python srcEarth/test/F12/run_F12.py --dry-run
```

F12 writes one case directory per spatial model/factor combination under `test_output/F12_gridless`, plus `F12_summary.csv`, `F12_result.json`, and `reference_F12_daynight_step_used.csv`. The repository reference table is `srcEarth/test/F12/reference_F12_daynight_step.csv`.

## F15 — Density normalization from differential flux

F15 checks the conversion from isotropic directional differential flux to density through the relativistic speed factor:

```text
n = 4π ∫ J(E)/v(E) dE
```

It uses `FIELD_MODEL NONE`, `R_INNER=0`, and four independent narrow TABLE top-hat spectra centered at `1, 10, 100, 1000 MeV`.  For each top-hat interval `[E1,E2]`, the references are closed form:

```text
F = 4π J0 (E2 - E1)
n = 4π J0/c * [sqrt(E2(E2+2m)) - sqrt(E1(E1+2m))]
```

The runner samples one point at the box center and one off-center point near the boundary; both must agree because the transport is all-open.

```bash
python srcEarth/test/F15/run_F15.py -np 4 -nt 16
python srcEarth/test/F15/run_F15.py --nintervals 256 --density-tol 1e-6
python srcEarth/test/F15/run_F15.py --dry-run
```

F15 writes one case directory per center energy under `test_output/F15_gridless`, plus `F15_summary.csv`, `F15_result.json`, and `reference_F15_top_hat_generated.csv` at the top level.  The default repository reference is stored in `srcEarth/test/F15/reference_F15_top_hat.csv`.
## F16 — Blocked-access zero-flux regression

F16 validates the exact zero-access limit in the gridless density/flux path. It uses `FIELD_MODEL NONE`, a nonzero power-law boundary spectrum, and places every diagnostic point inside the inner absorbing sphere. Since the trajectory classifier checks `R_INNER` before the first push, every sampled direction and energy must be forbidden:

```text
T(E)       = 0
J_local(E) = T(E) J_boundary(E) = 0
n          = 4π ∫ J_local(E)/v(E) dE = 0
F          = 4π ∫ J_local(E) dE = 0
```

The runner checks zero total density, zero total flux, zero channel fluxes, zero saved transmission, zero local spectrum, exact `J_local(E)=T(E)J_boundary(E)` closure, and confirms that the incident boundary spectrum is nonzero.

```bash
python srcEarth/test/F16/run_F16.py -np 4 -nt 16
python srcEarth/test/F16/run_F16.py --nintervals 160 --zero-flux-tol 1e-250
python srcEarth/test/F16/run_F16.py --dry-run
```

F16 writes `F16_summary.csv`, `F16_result.json`, and `reference_F16_blocked_access_used.csv` under `test_output/F16_gridless`. The repository reference/check table is `srcEarth/test/F16/reference_F16_blocked_access.csv`.

