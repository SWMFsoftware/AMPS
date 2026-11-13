# Running srcMoon diagnostic and demonstration runners

## Current-core integration regressions

This package compares every registered case in separate upstream/candidate
checkouts, then executes the Moon configuration-separated M0 and supporting
production-path checks. See [CORE_INTEGRATION.md](CORE_INTEGRATION.md) for the
full preparation, native build, application campaign and per-case audit
commands, including the authentic missing-fixture provenance and checkout-path
metadata preparation. Each campaign writes its own aggregate JSON/JUnit,
selection, owner logs and build/input/executable receipts. The audit writes
`result.json` and `comparison.csv`; its PASS means every upstream PASS was
preserved, with existing nonpasses still reported. It does not close M1/M2 or
approve a new I04/I13 reference. Use fresh output directories for every run.

## M2 supporting verification before M1 closure

### Description

This is the approved bounded development campaign. It clean-builds a dedicated
Na/Ar executable with chemistry OFF and individual weights ON, runs local
U08-U13, audits D01/D02/D03 provenance, and exercises the production LOLA
geometry/mover (I14 scope), surface distributions (I18 scope) and gas/surface
accounting (I19 scope). It then restores the canonical Na-only Moon/SPICE build.
Supporting results have separate counts. Formal I14-I20 remain SKIPPED until
same-campaign I13 is PASS; this command does not close M1 or M2.

### Commands

From the AMPS checkout, with no other build or model process using it:

```sh
python3 -B srcMoon/test/test_runner.py --phase m2-development \
  --mpi-ranks 8 \
  --output-dir test_output/srcMoon/M2-development-NEW
```

The runner uses `/home/vtenishe/SPICE/cspice` and
`/home/vtenishe/SPICE/Kernels`; override with `--spice-path` and
`--spice-kernels`. It removes `build/` before configuring and runs
`make -j 'OPT3=-O3 -g'`. Use a fresh output directory. Only when the current
generated tree already matches a retained dedicated M2 executable, reuse it:

```sh
python3 -B srcMoon/test/test_runner.py --phase m2-development \
  --mpi-ranks 8 --m2-amps /absolute/path/to/amps-moon-m2-na-ar \
  --output-dir test_output/srcMoon/M2-development-reuse-NEW
```

### Expected results and artifacts

`test_summary.json` records registered tests and separate `development_counts`.
The `m2-development/` child contains `development_summary.json`, build and
restore logs, `local/`, `data-audit/`, `supporting/I14/`, `supporting/I18/`,
`supporting/I19/` and genuine predecessor-blocked `formal/I14` through `I20`.
Each linked supporting run retains input, manifest, receipt, stdout/stderr and
resource history. I14/I18 write `comparison.csv/.png/.eps`; I19 writes
`m2_budget.csv`, `budget.png` and `budget.eps`. The I19 fixture is 20 ordinary
steps with declared incoming/outgoing particles and initial surface atoms;
these are verification initial conditions, never observations.

The November 12 retained example is
`test_output/srcMoon/m2-development-20261112T084925Z/campaign-5/`:
six local PASS, seven formal SKIPPED; supporting checks six PASS and one D02
SKIPPED, with zero FAIL/ERROR. D03 selected inputs are READY; D02 remains
PARTIAL pending authoritative illumination-integration epochs and operator
qualification. See STATUS.md section 5.7 for metrics and limitations.

## M1 event-budget and physical-duration diagnostic runner

This numerical diagnostic uses the approved spherical 480 km M1 configuration.
It checks real sampling storage and the configured Na chemistry/event ledger,
then advances the ordinary production model on eight MPI ranks under an
explicit wall/resource budget. It is separate from the demonstrations below
and cannot complete M1 or produce observational validation.

With a clean configured Moon/SPICE executable, run from the AMPS root:

```sh
python3 -B srcMoon/test/integration-tests/I11/run_duration_pilot.py \
  --amps ./amps --input srcMoon/examples/m1_memory200.in --mpi-ranks 8 \
  --wall-seconds 3600 --memory-limit-gib 200 \
  --output-dir test_output/srcMoon/I11-duration-pilot-NEW
```

The output root preserves prerequisite subdirectories, manifest and
`pilot_result.json`. `480km/` contains model stdout/stderr, input/receipt,
`particle_budget.csv`, `particle_sources.csv`, sampled `block_metrics.csv` and
resource histories. With completed samples, `480km/diagnostics/` contains
`metric_statistics.csv`, statistics/budget PNG and EPS, and provenance JSON.
The stage diagnostic is also CSV/PNG/EPS. A normal wall-budget stop is PARTIAL
with no duration candidate; COMPLETE means only that the unchanged coarsest
diagnostic gates passed and proposed durations still require review. There is
no particle-restart/resume claim. Use a new output directory for every run.

This document describes three independent demonstration runners:

1. `environment_demo/run.py` generates illumination, shadowing, and thermal
   state products from one production-kernel run.
2. `lola_surface_demo/run.py` generates and checks the triangulated LOLA
   surface without advancing the illumination or thermal models.
3. `sodium_exosphere_demo/run.py` executes the bounded production sodium
   exosphere for Tenishev et al. (2013) Case 1.

The runners have different purposes and output directories. A successful
execution of one runner does not imply that the other runner, or a formal
I-series validation test, has passed.

## Runner 3 — Tenishev-2013 Case-1 sodium exosphere

### Description

`srcMoon/test/sodium_exosphere_demo/run.py` runs the normal MPI production
executable with the exact accepted Case-1 physical receipt. It uses the fixed
uniform `2.65e16 atoms m^-2` Na reservoir, complete sticking/permanent loss,
the paper's `100+300 cos(chi)^0.25 K` surface temperature, MIV at
`1.67e22 atoms s^-1` and 6000 K, production PSD and TD, Combi radiation
pressure, and a 129600 s illuminated photoionization lifetime. Direct
sputtering is disabled. Ten 500-step sampling intervals are bounded by a
5000-step limit and particle restart writes remain disabled.

The paper associates Figs. 8--9 with 30 September 1991 but does not state the
observation time. The runner therefore has no UTC default: supply an exact UTC
on that date only after establishing it independently. A process-completion
PASS remains an `EXPERIMENTAL` demonstration. It is not I11 convergence or I32
validation, and no Flynn/Mendillo observation values are digitized from the
published plot.

### Commands

```sh
cd /data/vtenishe/Moon-test/AMPS

rm -rf build
./Config.pl -application=moon -amps-test=off \
  -spice-path=/home/vtenishe/SPICE/cspice \
  -spice-kernels=/home/vtenishe/SPICE/Kernels
make -j

RUN_OUTPUT=test_output/srcMoon/sodium-exosphere-demo/case1-19910930

python3 -B srcMoon/test/sodium_exosphere_demo/run.py \
  --simulation-start-utc 1991-09-30THH:MM:SSZ \
  --mpi-ranks 8 \
  --output-dir "$RUN_OUTPUT"

python3 -m json.tool "$RUN_OUTPUT/result.json"
```

Replace `HH:MM:SS` with an independently supported observation time. Use a new
output directory for each run. Alternatively add `--rebuild` to the Python
command to request the same clean configuration/build sequence from the
runner.

### Expected results and artifacts

The derived `application.in`, `stdout.log`, `stderr.log`, local U03/U04/U05/
U12/U14/U22 gate artifacts, `run_manifest.json`, and `result.json` are written
directly under `RUN_OUTPUT`. Untouched AMPS mesh, density, source-budget, limb,
column-density, and Na D-line brightness files are under
`RUN_OUTPUT/runtime/PT/plots/`. `result.json` records the required science
families and formal convergence/validation status separately. Completion
requires nonempty volumetric
Na density plus the per-stage `pic.LimbColumnDensity.*`,
`pic.ColumnDensityMap.spherical.*`, and
`pic.Moon.Anti-sunwardColumnIntegrals.*` products. Their headers contain Na
column density and both D-line brightness variables, and the limb product
contains tangent altitude in metres; an orbital ephemeris alone cannot satisfy
the runner.

The separate `pic.Moon.LimbIntegrals.*` and
`pic.Moon.AltitudeDistribution.*` tables are month-scale accumulated products:
their production writer is triggered only after an orbital-direction cycle.
They are not a valid completion requirement for this approximately 19-hour
bounded demonstration. To recover receipt generation after an older runner
rejected an otherwise completed preserved run, repeat the Python command with
`--finalize-existing`. This mode does not launch AMPS and fails closed unless
the exact UTC receipt, normal completion trailer, six local PASS artifacts,
and all four bounded-run product families are present.

The run manifest hashes the exact paper PDF, application receipt, executable,
production library, source/data roadmaps, and every model output. A missing or
changed paper, failed local gate, incompatible generated source configuration,
nonzero MPI exit, or absent science output is an `ERROR`, never a validation
PASS.

## Runner 1 — Illumination, shadowing, and thermal state

### Description

`srcMoon/test/environment_demo/run.py` produces a self-contained lunar-surface
demonstration from the normal `srcMoon` production kernels. It reconstructs a
triangulated lunar surface from the D01 `LDEM_4` LOLA digital elevation model,
evaluates SPICE Sun/Earth geometry, separates geometric night, terrain shadow,
and Earth occultation, and advances the production depth-resolved thermal
columns. The runner then writes complete CSV tables, CEA/Tecplot surface files,
PNG/EPS figures, logs, hashes, and a run manifest.

The Python runner does not contain an alternative illumination or thermal
model. A thin compiled adapter calls the same production geometry, absorbed
solar-flux, and thermal-column routines used by `srcMoon`. Before exporting a
case, it requires local tests U08, U09, and U11 to pass.

This is an `EXPERIMENTAL` meeting-demonstration runner, not a formal validation
test. A successful run does not establish I15, I17, or I20, and its model
output must not be represented as D02 illumination or D03 Diviner observation
data.

### Commands: one-hour demonstration

The single `environment_demo/run.py` invocation below generates all three
products together: illumination, geometric/terrain shadowing, and the lunar
thermal state. Separate physics runs are not required. Static PNG/EPS figures
are created automatically; the later visualization command is needed only for
an additional rotatable view.

Run these commands from a shell. The example uses a 100 km target surface
resolution and advances the thermal model by 3,600 seconds from the explicit
initial condition in `srcMoon/examples/lola_surface.in`.

```sh
cd /data/vtenishe/Moon-test/AMPS

test -r /home/vtenishe/SPICE/cspice/include/SpiceUsr.h
test -d /home/vtenishe/SPICE/Kernels
test -r /data/vtenishe/moon_validation_data/lola/raw/LDEM_4.IMG
test -r /data/vtenishe/moon_validation_data/lola/raw/LDEM_4.LBL

rm -rf build
./Config.pl -application=moon -amps-test=off \
  -spice-path=/home/vtenishe/SPICE/cspice \
  -spice-kernels=/home/vtenishe/SPICE/Kernels
make -j

RUN_OUTPUT=test_output/srcMoon/environment-demo/one-hour-100km

python3 -B srcMoon/test/environment_demo/run.py \
  --input srcMoon/examples/lola_surface.in \
  --build-dir build \
  --surface-resolution-m 100000 \
  --advance-seconds 3600 \
  --initialization input-initial \
  --output-dir "$RUN_OUTPUT"

python3 -m json.tool "$RUN_OUTPUT/result.json"
```

The output directory must either not exist or be empty. The runner deliberately
refuses to overwrite an earlier result. Choose a new `RUN_OUTPUT` name for each
case instead of deleting a scientifically meaningful run.

The one-hour case is an initial-condition transient. It is useful for checking
the complete data and visualization workflow, but it is not a converged lunar
day/night thermal solution. For an observation-independent periodic thermal
initialization, use the input-controlled production spin-up:

```sh
cd /data/vtenishe/Moon-test/AMPS

RUN_OUTPUT=test_output/srcMoon/environment-demo/periodic-spinup

python3 -B srcMoon/test/environment_demo/run.py \
  --input srcMoon/examples/lola_surface.in \
  --build-dir build \
  --surface-resolution-m 100000 \
  --initialization periodic-spinup \
  --snapshot-utc 2011-06-15T16:50:00.000Z \
  --output-dir "$RUN_OUTPUT"
```

Periodic spin-up can be substantially longer than the one-hour demonstration.
It uses the convergence criteria and cycle definition recorded in the input
file and has no arbitrary cycle cap. The requested snapshot must be at or
after the converged cycle end.

### Expected results and artifacts

For a successful run, the terminal prints a JSON object containing:

- `execution_status: PASS` for artifact generation and integrity checks;
- `local_prerequisites` with U08, U09, and U11 all equal to `PASS`;
- `surface_element_count: 6000`;
- `thermal_profile_row_count: 480000`, corresponding to 80 thermal layers for
  every logical surface cell;
- separate counts for terrain-visible and terrain-occluded cells;
- minimum and maximum modeled surface temperatures;
- I15, I17, and I20 recorded separately as `SKIPPED / NOT_VALIDATED`.

The exact shadow counts and temperatures depend on the selected epoch,
resolution, input parameters, and thermal initialization. They must not be
replaced with expected constants or adjusted to match a previous run.

With the `RUN_OUTPUT` value above, results are written under:

```text
test_output/srcMoon/environment-demo/one-hour-100km/
```

The principal files are:

| File | Contents |
| --- | --- |
| `surface_environment.csv` | One row per logical surface cell: location, terrain normal, local solar time, incidence, night/terrain/eclipse flags, direct-beam fraction, absorbed flux, and surface temperature. |
| `thermal_profiles.csv` | Every finite-volume thermal layer for every logical surface cell, including depth and temperature. |
| `lola_surface.cea` | LOLA-adjusted triangulated surface in an AMPS-loadable format. |
| `lola_surface.dat` | The same triangulated surface for Tecplot. |
| `moon_environment_illumination_shadowing.png` and `.eps` | Static illumination, terrain-shadow, eclipse, and three-dimensional direct-beam summary. |
| `moon_environment_thermal_state.png` and `.eps` | Static surface-temperature, local-time, distribution, and subsurface-profile summary. |
| `application.in` | Exact input snapshot used by the adapter. |
| `stdout.log`, `stderr.log` | Compiled-adapter output. A successful checked run should have an empty `stderr.log`. |
| `visualizer_stdout.log`, `visualizer_stderr.log` | Visualization receipt output and diagnostics. |
| `run_manifest.json` | Source revision, dirty paths, commands, build identity, compiler, roadmap hashes, and external-data hashes. |
| `result.json` | Machine-readable execution result, metrics, artifact hashes, and formal V&V status. |

### Visualization

Static figures are generated automatically. To open a rotatable 3-D view of
the completed case without replacing its original figures, write the custom
rendering into a separate directory:

```sh
cd /data/vtenishe/Moon-test/AMPS

RUN_OUTPUT=test_output/srcMoon/environment-demo/one-hour-100km
VIEW_OUTPUT="$RUN_OUTPUT/custom-temperature-view"

python3 -B srcMoon/tools/visualize_environment.py \
  "$RUN_OUTPUT/surface_environment.csv" \
  "$RUN_OUTPUT/thermal_profiles.csv" \
  --output-dir "$VIEW_OUTPUT" \
  --interactive-3d \
  --color-by surface_temperature_k \
  --selected-view-output "$VIEW_OUTPUT/selected_temperature_view.png"
```

Drag with the mouse to rotate the globe, press `S` to save the selected view,
and press `Q` to close it. Other supported interactive quantities are listed
by:

```sh
python3 -B srcMoon/tools/visualize_environment.py --help
```

### V&V status

All generated plots are explicitly labelled `EXPERIMENTAL MODEL OUTPUT - NOT
V&V EVIDENCE`. Formal comparison with illumination or Diviner observations
must follow the predecessor, convergence, data-readiness, and frozen-selection
requirements in `srcMoon/development_and_validation.md`.

## Runner 2 — Standalone LOLA surface map

### Description

The separate `srcMoon/test/lola_surface_demo/run.py` runner creates a
presentation-oriented surface package without running the illumination or
thermal model. It uses the production Moon input parser, native LOLA decoder,
icosphere generator, radial vertex adjustment, and CEA/Tecplot writers. It
also independently decodes the native D01 elevation selected for every mesh
vertex and records the comparison in `lola_dem_comparison.csv`.

### Commands: 100 km surface demonstration

Starting from the repository root, configure and build the SPICE-enabled Moon
application and run the 100 km surface case as follows:

```sh
cd /data/vtenishe/Moon-test/AMPS

test -r /home/vtenishe/SPICE/cspice/include/SpiceUsr.h
test -d /home/vtenishe/SPICE/Kernels
test -r /data/vtenishe/moon_validation_data/lola/raw/LDEM_4.IMG
test -r /data/vtenishe/moon_validation_data/lola/raw/LDEM_4.LBL

rm -rf build
./Config.pl -application=moon -amps-test=off \
  -spice-path=/home/vtenishe/SPICE/cspice \
  -spice-kernels=/home/vtenishe/SPICE/Kernels
make -j

LOLA_OUTPUT=test_output/srcMoon/lola-surface-demo/one-hundred-km

python3 -B srcMoon/test/lola_surface_demo/run.py \
  --build-dir build \
  --surface-resolution-m 100000 \
  --output-dir "$LOLA_OUTPUT"

python3 -m json.tool "$LOLA_OUTPUT/result.json"
```

The output directory must not contain an earlier run. The 100 km target is the
maximum requested pre-topography great-circle edge length, not an elevation
sampling interval. In the checked example it produced 10,242 vertices and
20,480 triangles; these counts are consequences of the icosphere refinement
and should be read from each run's `result.json`, not treated as acceptance
constants.

### Expected results and artifacts

The surface-runner output includes:

| File | Contents |
| --- | --- |
| `lola_surface.cea` | AMPS-loadable triangulated surface with radially adjusted vertex coordinates. |
| `lola_surface.dat` | Tecplot `FETRIANGLE` representation of the same surface. |
| `surface_vertices.csv` | Complete generated vertex table. |
| `lola_dem_comparison.csv` | Production radius/elevation and independently decoded D01 value for every vertex. |
| `lola_surface_demonstration.png` and `.eps` | Static surface maps, production-versus-D01 comparison, and residual distribution. |
| `application.in` | Exact derived input used for the run. |
| `verification/U08/result.json` | Required local LOLA geometry verification result. |
| `stdout.log`, `stderr.log`, `run_manifest.json`, `result.json` | Execution logs, provenance, hashes, metrics, and status. |

### Visualization

To inspect the generated geometry interactively, run:

```sh
cd /data/vtenishe/Moon-test/AMPS

LOLA_OUTPUT=test_output/srcMoon/lola-surface-demo/one-hundred-km

python3 -B srcMoon/tools/visualize_surface.py \
  "$LOLA_OUTPUT/lola_surface.cea" \
  --interactive-3d \
  --color-by elevation-m \
  --reference-radius-m 1737400 \
  --vertical-exaggeration 20 \
  --frame-label "MEAN EARTH/POLAR AXIS OF DE421" \
  --selected-view-output "$LOLA_OUTPUT/selected_view.png"
```

The reference radius and frame label above come from the native `LDEM_4.LBL`
metadata and must not be reused for another product without checking its
label. Drag with the mouse to rotate the surface, press `S` to save the chosen
camera view, and press `Q` to close it. The vertical exaggeration affects only
visualization; it does not alter the CEA/Tecplot surface coordinates.

### V&V status

This runner records artifact generation as `execution_status=PASS` only after
U08 and its independent per-vertex comparison pass. It remains an
`EXPERIMENTAL / MEETING_DEMONSTRATION` product: formal I14 is recorded as
`SKIPPED / NOT_VALIDATED` until its linked-production and predecessor gates
are satisfied.
