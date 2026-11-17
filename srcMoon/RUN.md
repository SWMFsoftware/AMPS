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

To select every registered core suite, independently of its routine campaign
label, run this from the active AMPS root:

```sh
/usr/bin/python3 -B src/test/run_tests.py \
  --config test/qualification-machine/passion.json --match 'CORE-*' \
  --jobs 8 --build-jobs 16 \
  --output-dir "test_output/core-full-$(date -u +%Y%m%dT%H%M%SZ)"
```

`--jobs` controls safe concurrent suites, not every test's MPI rank count.
The machine profile and registered case contracts determine MPI layouts;
native configuration changes remain serialized. The campaign writes
`amps-core-tests.json`, JUnit XML, per-suite `logs/` and retained `artifacts/`.
Replacing `--match 'CORE-*'` with `--all` also selects registered applications;
`--full-validation` additionally requests expensive qualification modes and
cannot be combined with either selector. `make core_test_nightly` selects the
routine nightly campaign rather than every registered core suite.

The historical current-core integration baseline had 1,085 core cases:
843 PASS, 240 SKIPPED, one strict-coverage FAIL and one CME-refinement timeout
ERROR. These are retained baseline dispositions, not results of the current
M2 changes. Later registries may contain additional cases; use `--list` and
the actual aggregate report for current membership/outcomes. The November 13
standalone runs also exposed Python dependency failures; those reports remain
preserved and must not be silently replaced by the older baseline.

### Python environment for this machine

The maintained `passion.json` profile selects `/usr/bin/python3`. Use that
same interpreter for development runners and check its dependencies explicitly:

```sh
/usr/bin/python3 -c 'import numpy, scipy, matplotlib, mpmath; print("Python numerical/plot dependencies available")'
```

The top-level control layer uses the standard library, while independent
numerical references, confidence intervals and PNG/EPS output need these
packages. A package installed in a shell-selected `python3` virtual environment
does not establish availability to `/usr/bin/python3`. On another machine use
a reviewed profile with actual interpreter/compiler/MPI/SPICE paths. These
dependency checks establish execution readiness only, not lunar validation.

## M2 supporting verification before M1 closure

### Description

This is the approved bounded development campaign. It requalifies the
configuration-separated M0 baseline, clean-builds dedicated neutral Na/Ar
(chemistry OFF, individual weights ON, explicit Mersenne Twister), runs U08-U13
and data QA, and exercises the original terrain/surface/budget fixtures. It
then runs 36 short cases: source-fed startup and source-free surface cohorts,
three timesteps, three independent particle weights and two declared seeds.
Finally it prepares a deterministic source-rate integral and restores the
canonical Na-only Moon/SPICE build. Numerical diagnostic PASS does not close
formal I14-I20: same-campaign I13 and all individual data/operator gates remain
required. No observation is generated, fitted, or scored by the short sweep.

The complete selected command has one 3,600-second monotonic wall budget,
including all required builds, tests, QA, plots and cleanup; 600 seconds are
reserved for termination/cleanup. Descendant MPI/build processes are monitored
across process sessions. Resource/deadline exhaustion is ERROR with preserved
logs and completed results. The guard does not make an incomplete campaign a
PASS. Long qualification is a separate, explicit cluster tier.

### Commands

Use a dedicated checkout with no concurrent configuration/build process:

```sh
/usr/bin/python3 -B srcMoon/test/test_runner.py --phase m2-development \
  --mpi-ranks 8 --stream \
  --m2-input srcMoon/examples/m2_development.in \
  --campaign-wall-seconds 3600 --cleanup-reserve-seconds 600 \
  --output-dir test_output/srcMoon/M2-development-NEW
```

Default SPICE paths are `/home/vtenishe/SPICE/cspice` and
`/home/vtenishe/SPICE/Kernels`; override with `--spice-path` and
`--spice-kernels`. Every generated tree is removed before configuration;
native builds use `make -j16 'OPT3=-O3 -g'`. Always select a new output path.
The development sweep is diagnostic, with 24/48/96 production steps and a
matched physical duration; it does not pretend to finish a 500-step science
sample. Sampling length remains 500; restart writes are suppressed.

For a focused rerun only, `--m2-amps /absolute/path/to/amps-moon-m2-na-ar`
requires matching binary/native-source/generated chemistry/layout/RNG hashes.
Reuse skips clean M0 requalification and the final normal-build restoration;
therefore it is not evidence that the full clean command meets the hour.

### Expected results and artifacts

`test_summary.json` records formal results and separate `development_counts`.
The child `m2-development/` contains `budget_result.json`, process-tree RSS and
allocation history in `budget_resources.csv`, complete `campaign.log`, and
`development_summary.json`. Actual phases are under its `campaign/` directory:

- `m0/`: all configuration-separated I01-I10 evidence and strict nightly diff;
- `build/`: dedicated executable, clean build log and source/config/RNG hashes;
- `local/` and `data-audit/`: local oracle results and package readiness;
- `supporting/I14`, `I18`, `I19`: the original positive-coverage fixtures;
- additional native parser/operator checks retain their own supporting logs
  and receipts; actual I15/I20 exporters are exercised without observations,
  while formal statuses remain separately gated;
- `numerical/<source-fed|cohort>-d<D>-w<W>-s<SEED_INDEX>/`: each case's
  `application.in`, `run_manifest.json`, receipt, streams, `result.json`,
  `m2_budget.csv`, completed `window_budget.csv`, and `comparison.csv/.png/.eps`;
- cohort cases additionally retain `m2_cohorts.csv`, `m2_hazards.csv`, and
  `cohorts.png/.eps`, with censored first residence, occupancy and tag closure;
- `source-normalization/`: source-only integral, preserved input/receipt and
  `result.json`, plus the requested actual full I19 initial-state CSV. The
  candidate records `qualified:false`; READY describes only a model-input
  normalization and supplies no initial-state approval or observations;
- `formal/I14` through `I20`: genuine predecessor-blocked results when I13 is
  absent; no synthetic comparison graphics are produced for these skips.

At present formal M1/M2 closure remains open and D02 is PARTIAL. Use STATUS.md
for the measured campaign outcome; code or preflight availability is not
qualification. See [M2_NUMERICAL_PROTOCOL.md](test/M2_NUMERICAL_PROTOCOL.md)
for numerical margins and detailed physics/implementation decisions.

### New production parser and environment-operator checks

The bounded development command now includes the new input-control/native
adapter tests. They exercise explicit finite-disc settings and rejected
incomplete controls, production fixed-site SPICE geometry, and actual
illumination/thermal exporters. They do not fabricate an I13 predecessor,
replace a D02/D03 reference, or assign an observational acceptance threshold.
The current native campaign determines their outcome; registered code or a
successful parser preflight is not an executed numerical PASS.

For a focused operator rerun with a dedicated Na/Ar executable and its matching
current generated build tree, build in an idle private checkout first if that
matching tree is unavailable:

```sh
/usr/bin/python3 -B srcMoon/test/build_m2_qualification.py \
  --output-dir test_output/srcMoon/M2-operator-build-NEW \
  --spice-path /home/vtenishe/SPICE/cspice \
  --spice-kernels /home/vtenishe/SPICE/Kernels

/usr/bin/python3 -B srcMoon/test/run_m2_operator_checks.py \
  --amps test_output/srcMoon/M2-operator-build-NEW/amps-moon-m2-na-ar \
  --input srcMoon/examples/m2_development.in --mpi-ranks 8 \
  --test-runtime-seconds 600 \
  --output-dir test_output/srcMoon/M2-operator-wiring-NEW
```

The 600 seconds is an explicit focused command budget, not a thermal or
observational tolerance or a claim that the separate build was timed inside
that focused command. The complete global development invocation covers all
required phases within its one-hour envelope. The exporter fixture uses declared geometry queries
and a deliberately millisecond forcing replay to check state/coordinate/depth
export wiring. It does not represent a synodic spin-up, qualified initial
lunar thermal state or annual polar observation. Expect
`operator_checks_summary.json`, separate native query/model CSVs, inputs,
manifests, stdout/stderr, receipts and resource logs. Formal I15/I20 remain
SKIPPED. The strict build-identity check refuses live-header/library mismatch;
see [BUILD_M2.md](test/BUILD_M2.md) before reusing a retained executable after
normal-profile restoration.

The new optional production finite-disc model samples terrain and opaque-Earth
visibility on the same rays. Thermal forcing integrates the local incidence
and directional absorptivity before entering the unchanged conduction solver.
PSD instead uses unabsorbed direct photon incidence: no thermal Bond albedo is
multiplied into its source. The serial SOURCE-BEGIN cache supplies one epoch
to ordinary per-element rates, total/CDF selection and finite-inventory hazards;
particle departure timing does not silently change that coefficient epoch.
The core callback defaults to NULL, preserving M0/M1 and other applications.
Details, units, API names and limits are in
[PSD_ILLUMINATION_CONTRACT.md](test/PSD_ILLUMINATION_CONTRACT.md).

### Deferred cluster tier

Formal preparation preserves the reviewed volume AMR, sources, thermal state
and epoch. It requires Ar already enabled, an evolving environment and the
declared zero inventory used by the numerical protocol. The short I18 helper
cannot silently replace those physical settings. Initial-state evidence must
identify the actual applied uniform/periodic method, native full-column bytes,
epoch, executable, physical input and actual furnished SPICE/LOLA bytes.
Missing evidence remains SKIPPED; no state becomes qualified by being exported.

`thermal_initialization_mode = periodic-spinup` invokes the existing solver
before source/weight queries; its explicit cycle end must be simulation t0.
`declared-uniform` retains the development transient and is not automatically
a lunar equilibrium. I19 can record actual initial columns through
`m2_verification_record_initial_thermal_state = true` within the complete
verification group. I20 records its own full starting columns and keeps annual
thermal and cohort epochs distinct. These model diagnostics are separate from
particle/inventory checkpoint writes.

Generate a **formal** source normalization only with the explicit
`prepare_m2_normalization.py --tier cluster` option and the reviewed full input;
the default development normalization uses the bounded fixture and cannot
qualify a different full mesh. Initialization compute/storage is part of the
measured cluster cost for all cases. See
[the initialization contract](test/M2_INITIAL_STATE_CONTRACT.md) and
[formal staging](test/FORMAL_STAGING_CONTRACT.md) for required evidence and
normalization commands. No actual qualifying state, interval or READY approval
has been created by these instructions.

Export the fixed 288-case plan without executing native cases:

```sh
python3 -B srcMoon/test/test_runner.py --phase m2-cluster --mpi-ranks 8 \
  --stream --output-dir test_output/srcMoon/M2-cluster-plan-NEW
```

Expect `m2-cluster/job_plan.json`, `preflight.json` and `cluster_summary.json`;
formal I19/I20 remain SKIPPED when execution is unselected or required evidence
is unavailable. The runner submits no remote jobs. All necessary execution,
proof and resource arguments are documented in
[M2_NUMERICAL_PROTOCOL.md](test/M2_NUMERICAL_PROTOCOL.md#cluster-qualification).
The 300 core-hour default is one aggregate planning allowance; a measured
forecast is required before launch and may establish that this design is too
expensive. No timestep/source/resolution physics is changed to force it to fit.

### D02 research and native-metadata audit

This read-only preflight is also registered automatically in the global M2
development campaign. It verifies frozen native bytes, labels, frame and
observable; it cannot supply missing averaging epochs or validate illumination.
To run just this audit from the AMPS root:

```sh
python3 -B srcMoon/test/data_processing/audit_d02_temporal_contract.py \
  --output-dir test_output/srcMoon/D02-native-audit-NEW
```

Expect `result.json` with SKIPPED and exit 77 while the exact-product temporal
contract is unresolved; corruption is ERROR/exit 2. No model is launched and
no comparison plot is produced. Use a new output directory outside the
canonical validation-data tree. For the combined package audit and its ten
automatically registered metadata protocol checks:

```sh
python3 -B srcMoon/test/data_processing/audit_m2_packages.py \
  --output-dir test_output/srcMoon/M2-data-audit-NEW
```

That command preserves individual package statuses in `result.json`, the
native audit in `d02_native_temporal_audit.json`, and infrastructure checks in
`d02_protocol_checks.log`. Its zero exit means no audit ERROR occurred; it
does not mean every dataset is READY. Read the counts and scopes. The current
D02 reference stays PARTIAL. See
[the research report](test/D02_ILLUMINATION_RESEARCH.md) and
[the approved I19/I20 numerical protocol](test/M2_NUMERICAL_PROTOCOL.md).
Package audit does not score observations or complete a formal gate. The
active protocol controls numerical diagnostics; existing conservation,
predecessor, observation-selection and hold-out gates remain binding.

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

## Initial M3 atomic rate verification

This runner verifies the approved Maxwellian/Voronov rate service. It runs the
unchanged clean M0 prerequisite (including the reviewed four-rank nightly
reference), registered U16 kernel/runtime controls, and an **eight-rank** native
callback fixture with all four neutral/ion pairs compiled. The native callback
uses real generated particle-record/species layouts; no particle reaction or
product trajectory is simulated. Formal I21-I24 remain explicitly SKIPPED.

From the selected checkout root:

```sh
OMP_NUM_THREADS=1 /usr/bin/python3 -B srcMoon/test/test_runner.py \
  --phase m3-development --mpi-ranks 8 \
  --campaign-wall-seconds 3600 --cleanup-reserve-seconds 600 \
  --spice-path /home/vtenishe/SPICE/cspice \
  --spice-kernels /home/vtenishe/SPICE/Kernels \
  --output-dir test_output/srcMoon/M3-development-NEW
```

Replace `NEW` for each execution; prior evidence is not overwritten. This
native campaign rebuilds the selected checkout several times and restores
ordinary Moon/SPICE afterward. Use an isolated checkout when the active tree
must retain another application/configuration. Configuration/build logs,
executable/source/generated-header identities and resource monitoring are
preserved within the one-hour campaign budget. As in M0, each configuration
starts by removing its generated `build` directory.

Expected files under `M3-development-NEW/m3-development/`:

- `m0/`: ten real linked baseline outcomes, complete U-series and zero-byte
  `I04/test_Moon.diff` when the reviewed nightly passes.
- `native-rate-callback/`: derived eight-species compile input identity,
  configure/build/probe logs, generated ABI headers, `application.in`,
  `run_manifest.json`, per-rank callback CSVs and `result.json`.
- `restore/`: normal Moon/SPICE clean build log and executable.
- `development_summary.json`, `execution_budget.json`, `campaign.log`: distinct
  local/supporting/formal outcomes, errors and process-tree resource evidence.

Run only the registered fast U16 case without reconfiguring AMPS:

```sh
/usr/bin/python3 -B srcMoon/test/run_all_test.py --skip-build --only U16 \
  --output-dir test_output/srcMoon/U16-Maxwellian-NEW
```

U16 independently compiles the maintained production rate and runtime-input
services. It preserves native/reference hashes, comparison CSV/PNG/EPS and
malformed-input/control outcomes. Its PASS does not qualify I23 survival,
product conservation or measured forcing. The complete `.in` example is
`examples/m3_atomic_rates.in`; its declared electron bath is an analytical
fixture, not a lunar observation. Charge-exchange and molecular-kinetics
proposals are in `test/M3_*_PROPOSAL.md`, and STATUS.md section 5.15 records
current evidence and unresolved decisions.
