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

## M3 supporting rates, packet transactions and force tracing

The full remaining work is specified in the
[M3 completion proposal](test/M3_COMPLETION_PROPOSAL.md) and its approved
[implementation contract](test/M3_IMPLEMENTATION_CONTRACT.md). The command
below runs local/supporting verification; it does not qualify I21–I24.

This runner verifies the production rate and event services. It runs the
unchanged clean M0 prerequisite (including the reviewed four-rank nightly
reference), registered U15/U16/U17 support, and an **eight-rank** native rate
callback with all four neutral/ion pairs compiled. With the explicit isolated
checkout acknowledgment it additionally runs actual PIC transaction controls,
both tracers with the actual Moon forces, and the normal-path constant-history
reaction controller on one/eight ranks. Those use
declared analytical storage/force fixtures, without assuming physical CEX or
photochemical recoil distributions. Formal I21-I24 remain explicitly SKIPPED.

From a disposable selected checkout root containing the current source changes
and the normal AMPS build prerequisites:

```sh
OMP_NUM_THREADS=1 /usr/bin/python3 -B srcMoon/test/test_runner.py \
  --phase m3-development --mpi-ranks 8 --isolated-checkout \
  --campaign-wall-seconds 3600 --cleanup-reserve-seconds 600 \
  --spice-path /home/vtenishe/SPICE/cspice \
  --spice-kernels /home/vtenishe/SPICE/Kernels \
  --output-dir test_output/srcMoon/M3-development-NEW
```

Replace `NEW` for each execution; prior evidence is not overwritten. The flag
acknowledges a disposable checkout; it does not create one or copy uncommitted
changes automatically. Omit it to skip the derived packet/force builds. This
native campaign rebuilds the selected checkout several times and restores
ordinary Moon/SPICE afterward. Use an isolated checkout when the active tree
must retain another application/configuration. Configuration/build logs,
executable/source/generated-header identities and resource monitoring are
preserved within the one-hour campaign budget. As in M0, each configuration
starts by removing its generated `build` directory.

The tracer examples select frozen prescribed LSO fields and
`m2_environment = disabled`. Their local particle clocks cannot be silently
combined with the current globally timed dynamic thermal/reservoir update;
that combination is refused by both parser and direct installation. The new
controller fixtures use explicit constant analytical hazards and providers,
without attaching the physical D04/D05/CEX registries automatically.

Expected files under `M3-development-NEW/m3-development/`:

- `m0/`: ten real linked baseline outcomes, complete U-series and zero-byte
  `I04/test_Moon.diff` when the reviewed nightly passes.
- `native-rate-callback/`: derived eight-species compile input identity,
  configure/build/probe logs, generated ABI headers, `application.in`,
  `run_manifest.json`, per-rank callback CSVs and `result.json`.
- `native-particle-backend/`: fresh all-species producer and generated ABI,
  compile/build logs, one/eight-rank real PIC allocation/weight/rollback controls
  and per-rank metrics. No artificial allocator replaces the production buffer.
- `native-moon-force-tracers/`: both actual Moon tracer orders, per-layout
  position/velocity/reference CSV and diagnostic PNG/EPS, production-force
  decomposition and independent integration-refinement evidence.
- `native-reaction-controller/`: both orders and layouts through normal Moon/PIC
  callbacks, competing constant analytical hazards, preserved clocks/event times,
  complete-packet continuation, recycled-source and capacity-refusal controls,
  event comparison CSV/PNG/EPS and the exact producer/ABI identities. These
  mechanics do not install a physical photon/electron/CEX recoil model.
- `formal/I21` through `I24`: the actual registered owners' predecessor and
  missing-capability receipts; no synthetic observation or placeholder PASS.
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

### M3 independent exposure, common-clock and protocol preparation

These runners verify declared mathematics and protocol/mechanical behavior;
they do not produce lunar observational validation or enable missing daughter
physics. Run from the AMPS repository root, with a new output directory each
invocation. A C++17 compiler is required. PNG/EPS diagnostics require a Python
installation with matplotlib (the development machine supplies this through
`/usr/bin/python3`).

```bash
/usr/bin/python3 -B srcMoon/test/reaction_execution/piecewise_hazard/run.py \
  --output-dir test_output/srcMoon/m3-piecewise-NEW
/usr/bin/python3 -B srcMoon/test/reaction_execution/exposure_bridge/run.py \
  --output-dir test_output/srcMoon/m3-exposure-bridge-NEW
/usr/bin/python3 -B srcMoon/test/run_m3_transport_clock.py --local-only \
  --output-dir test_output/srcMoon/m3-clock-local-NEW
/usr/bin/python3 -B srcMoon/test/test_m3_qualification_protocol.py \
  --output-dir test_output/srcMoon/m3-protocol-controls-NEW
/usr/bin/python3 -B srcMoon/test/m3_protocol.py \
  --output-dir test_output/srcMoon/m3-protocol-preparation-NEW \
  --plan-comparisons 48
```

Each runner writes its own `result.json` and logs beneath the selected output.
The numerical runners retain source/compiler/executable identity, independent
analytical CSV comparisons and PNG/EPS diagnostics where declared. Protocol
preparation currently exits 77/SKIPPED: structural identities can pass while
physical inputs and formal qualification are incomplete. No particle run is
launched by the preparation command. The comparison family size 48 is only a
planning example, not a frozen or approved physical campaign.

The global U-series invokes these services under U15 supporting children and
U16 orchestration controls. The existing extended `m3-development` phase in a
disposable checkout additionally runs the linked common-clock fixture at one
and eight MPI ranks against its freshly built packet backend, then restores
ordinary Moon configuration. It shares one hour across setup/build/tests/plots/
restoration/cleanup. See [the protocol](test/M3_QUALIFICATION_PROTOCOL.md),
[physical model](test/M3_PHYSICAL_MODEL_SPECIFICATION.md) and
[clock contract](test/M3_COMMON_TRANSPORT_CLOCK_CONTRACT.md) before interpreting
results; local or mechanical PASS does not close I21–I24 or prior milestones.

The extended campaign also runs full normal application startup after restoring
its canonical Na-only SPICE build: both orders, eight ranks, two actual
`amps_time_step` calls. This checks UTC-to-SPICE initialization, normal
sampling/source ordering and completed common endpoints. It reuses the restored
archive and does not perform another application build. Its result, generated
header snapshots and `normal_build_evidence.json` remain under the campaign's
restore/startup directories. With that authentic build still present, repeat
only this check in the same disposable checkout using a fresh output:

```bash
/usr/bin/python3 -B srcMoon/test/run_m3_transport_clock.py \
  --normal-startup-only \
  --build-evidence /absolute/path/to/campaign/restore/normal_build_evidence.json \
  --output-dir test_output/srcMoon/m3-normal-clock-startup-NEW \
  --wall-seconds 600
```

The runner verifies the actual archive, source, generated headers and successful
build log before linking the thin driver. Missing build provenance is SKIPPED;
changed bytes or contradictory proof are ERROR. Use only the remaining original
campaign budget when running this inside an existing development package.

## Approved M3 reduced-model components and ion boundary

### Source preparation

Run from the repository root. These commands create research/source evidence,
not observations or a complete physical state model. Use a Python interpreter
with Matplotlib for the component comparison plots (`/usr/bin/python3` on the
development machine). Use a new output directory each time:

```bash
MOON_COMPONENT_OUTPUT="test_output/srcMoon/m3-components-$(date -u +%Y%m%dT%H%M%SZ)"
/usr/bin/python3 -B srcMoon/test/physical_state_sources/prepare.py \
  --fetch --output-dir "$MOON_COMPONENT_OUTPUT/state-sources"
/usr/bin/python3 -B srcMoon/test/sodium_partial_capture/prepare_source.py \
  --source-html test_output/srcMoon/m3-atomic-online-databases-20261118T083400Z/collisiondb_Na_H1s_119990.html \
  --output-dir "$MOON_COMPONENT_OUTPUT/na-partial-source"
```

The first command downloads public NIST isotope/constants and H/Na ground-
threshold records. The second reproduces twelve measured partial capture
nodes from the already retained native IAEA response. Supply the actual
source response path when using another machine; a missing or changed source
is explicitly unavailable/corrupt, never replaced. Raw packet bytes and
provenance are retained. Neither command modifies canonical validation data.

Acquire the separate oxygen state/transition research packet, then replay it
offline into a fresh directory if needed:

```bash
/usr/bin/python3 -B srcMoon/test/state_cascades/prepare_sources.py \
  --download --output-dir "$MOON_COMPONENT_OUTPUT/oxygen-sources" \
  --wall-seconds 3600
/usr/bin/python3 -B srcMoon/test/state_cascades/prepare_sources.py \
  --source-packet "$MOON_COMPONENT_OUTPUT/oxygen-sources" \
  --output-dir "$MOON_COMPONENT_OUTPUT/oxygen-replay" \
  --wall-seconds 3600
```

This preserves every native level/line, missing A value and selected fine-
structure/multipole component. It does not normalize stored transitions to a
physical branching distribution or infer a total lifetime. The source QA
can pass while cascade physics remains NOT QUALIFIED. Offline replay verifies
the full acquisition manifest before processing; source packets are immutable.

### Individual component runners

```bash
/usr/bin/python3 -B srcMoon/test/test_m3_reduced_kinematics.py \
  --output-dir "$MOON_COMPONENT_OUTPUT/reduced-recoil"
/usr/bin/python3 -B srcMoon/test/test_m3_ion_collector.py \
  --output-dir "$MOON_COMPONENT_OUTPUT/ideal-collector"
/usr/bin/python3 -B srcMoon/test/test_m3_sodium_partial_capture.py \
  --reference-dir "$MOON_COMPONENT_OUTPUT/na-partial-source" \
  --mass-inventory "$MOON_COMPONENT_OUTPUT/state-sources/reference_mass_inventory.json" \
  --output-dir "$MOON_COMPONENT_OUTPUT/na-partial-capture"
/usr/bin/python3 -B srcMoon/test/test_m3_helium_silica.py \
  --source-pdf test_output/srcMoon/m3-ion-surface-cex-data-sources-20261117T225600Z/tellez2012_silica_neutralization.pdf \
  --output-dir "$MOON_COMPONENT_OUTPUT/he-silica"
/usr/bin/python3 -B srcMoon/test/test_m3_oxygen_sources.py \
  --source-packet "$MOON_COMPONENT_OUTPUT/oxygen-sources" \
  --output-dir "$MOON_COMPONENT_OUTPUT/oxygen-source-controls"
```

Expect structured result/manifest, source/compiler/executable identities,
compile/stdout/stderr logs and CSV/PNG/EPS comparisons in each requested
directory. Reduced photon and three-neutral models use explicit supplied
states and are not a complete moving-parent LSO chemistry provider. The Na
fixture uses native partial-channel beam nodes and an experimental forward
angle; native initial Na state metadata remain unresolved. The He/silica
fixture evaluates a conditional charge-survival law in its source coordinate;
conversion from arbitrary physical vectors remains Unavailable. None is a
lunar observational comparison or formal M3 PASS.
The oxygen runner has fifteen offline source/parser/corruption controls; it
builds no application and enables no physical cascade. The global U15 child
uses the retained research packet by default, makes no network requests and
reports SKIPPED when that packet is unavailable on another machine.

### Global/native execution and runtime boundary selection

The global runner includes these children through the existing U15/U17
supporting namespace, retaining each runner/log and actual status. To include
normal clean application builds, reviewed `make -j test_Moon`, actual PIC
packets, sphere collector dispatch and both force-tracing orders, run in a
**disposable, correctly installed checkout**:

```bash
/usr/bin/python3 -B srcMoon/test/test_runner.py \
  --phase m3-development --isolated-checkout --mpi-ranks 8 \
  --campaign-wall-seconds 3600 \
  --output-dir test_output/srcMoon/m3-native-NEW
```

`--isolated-checkout` acknowledges that the selected checkout owns its derived
build/configuration; it does not create a copy. Install/relocate that checkout's
site `PTDIR/MYDIR` before execution so it cannot compile against another tree's
generated headers. Include setup/cleanup in the one-hour package allowance.
Full U/I outcomes and supporting counts remain separate. Model-data/predecessor
gaps correctly leave formal chemistry SKIPPED even when components pass.

In a complete existing tracing `.in` receipt, `ion_surface_boundary =
unsupported` is the default. Explicit `ideal-collector` additionally requires
an installed state/formation-energy resolver and ledger matching actual compiled
species/masses/charges. Native analytical fixtures install that provider;
an ordinary CLI switch alone does not. The collector has no neutralization,
heat, implantation depth, sputtering, charging-feedback or release model.
Detailed assumptions, source gaps and exact approval limits are in
[the decision document](test/M3_DAUGHTER_AND_ION_SURFACE_DECISIONS.md).

### Source-state admission and conditional physical Na cohort

The new [cohort runner instructions](test/sodium_cohort/README.md) give the full
source-profile generation, deterministic native input projection and eight-rank
source-preflight commands. It uses an authentic current native producer from an
isolated normal build; it does not rebuild the active checkout or guess an
application CLI. The original partial profile correctly gives SKIPPED/exit77
for unavailable Na+. The separate source-derived atomic profile now supplies
an opt-in ground-ion/energy crosswalk. Its source inspection, physical execution,
numerical/statistical qualification and formal predecessor gates remain distinct.

Independent supporting checks, from this repository root and with NEW output
directories, are:

```bash
/usr/bin/python3 -B srcMoon/test/test_m3_physical_reaction_events.py \
  --output-dir test_output/srcMoon/photon-composition-NEW
/usr/bin/python3 -B srcMoon/test/test_m3_trajectory_exposure.py \
  --output-dir test_output/srcMoon/trajectory-exposure-NEW
/usr/bin/python3 -B srcMoon/test/test_m3_qualification.py \
  --output-dir test_output/srcMoon/qualification-controls-NEW
/usr/bin/python3 -B srcMoon/test/test_m3_sodium_ledger.py \
  --output-dir test_output/srcMoon/sodium-ledger-controls-NEW
```

These produce structured results, logs, hashes and comparison CSV/PNG/EPS.
Analytical state fixtures verify software and are labelled accordingly; they
are not physical sodium profiles or observations. The normal native campaign
also checks pre-buffer registration, full daughter genealogy/counters, actual
PIC initial-state readback, explicit SI charge versus legacy constants and
accepted light-to-shadow endpoint semantics on1/8MPI with both tracing orders.

Once authoritative source/state inputs and all predecessors exist, the frozen
per-cohort execution interface is:

```bash
MOON_NATIVE_CHECKOUT=/ABSOLUTE/PATH/TO/DISPOSABLE-AMPS
/usr/bin/python3 -B "$MOON_NATIVE_CHECKOUT/srcMoon/test/m3_qualification_executor.py" \
  --case-plan /ABSOLUTE/PATH/TO/FROZEN_COHORT_PLAN.json \
  --campaign-root /ABSOLUTE/PATH/TO/SAME_CAMPAIGN \
  --output-dir test_output/srcMoon/physical-cohort-NEW
```

The installed producer identity is `normal-source-bound-sodium-cohort-v1` for
I21. The plan's `native_adapter` schema is
`srcMoon-M3-normal-sodium-execution-v1`; it requires `producer_manifest`
`{path,sha256}`, `state_mapping` with parent/ion/electron IDs, and `runs` with
`mpi_ranks`, `physical_steps`, `input {path,sha256}` and the COMPLETE unique
`input_documents [{path,sha256},...]` including the root and application include
graph. Predeclare at least three identical step-count levels for1/8MPI. The
plan also contains all existing approved cohort/statistical fields. The controlled workflow below uses an authenticated real atomic profile;
an analytical state fixture cannot replace that physical source input.
Native77 remains SKIPPED, malformed/changed inputs are ERROR, evaluated numerical
or statistical disagreement is FAIL. A successful cohort is evidence toward the
whole frozen matrix, not a milestone PASS. All prior attempts remain immutable.


### P02-approved controlled sodium commissioning

This runner exercises a finite source-bound Na→Na+ cohort through normal AMPS
startup, both production species, the force-aware tracing dispatcher, chemistry
controller and atomic particle transaction. The supplied example selects the
second-order tracer and the explicitly compiled zero-force analytical reference.
It uses 128 independent original parents, equal represented weights, actual
frozen SPICE geometry, the quiet D04 endpoint and the approved P02 mechanical
limits. It is a commissioning experiment, not a complete lunar-exosphere run,
statistically adequate I21 campaign, observation comparison or M3 closure.
Read [the controlled protocol](test/M3_P02_CONTROLLED_NA_PROTOCOL.json) and
[the preparation contract](test/sodium_cohort/PREPARATION.md).

First prepare the real atomic profile and source-matched conservative producer
as documented in [the isolated build instructions](test/source_qualification_build/README.md).
The producer checkout must remain in that configuration. The active checkout
and its build do not need to be reconfigured. Supply actual absolute paths:

```sh
MOON_ACTIVE_ROOT=/absolute/path/to/AMPS
MOON_SOURCE_CHECKOUT=/absolute/path/to/disposable-source-matched-AMPS
MOON_ATOMIC_PROFILE=/absolute/path/to/atomic-profile/profile.json
MOON_ENERGY_CROSSWALK=/absolute/path/to/atomic-profile/na_source_energy_crosswalk.csv
MOON_PRODUCER=/absolute/path/to/source-build/normal_build_evidence.json
MOON_CASE_ROOT=/absolute/path/to/previously-absent-commissioning-package
mkdir "$MOON_CASE_ROOT"
/usr/bin/python3 -B "$MOON_ACTIVE_ROOT/srcMoon/test/sodium_cohort/prepare_physical_case.py" \
  --profile "$MOON_ATOMIC_PROFILE" --crosswalk "$MOON_ENERGY_CROSSWALK" \
  --application-input "$MOON_ACTIVE_ROOT/srcMoon/examples/m3_sodium_p02_commissioning.in" \
  --case-json "$MOON_ACTIVE_ROOT/srcMoon/examples/m3_sodium_p02_commissioning.json" \
  --output-dir "$MOON_CASE_ROOT/prepared"
cd "$MOON_SOURCE_CHECKOUT"
/usr/bin/python3 -B srcMoon/test/run_m3_sodium_cohort.py \
  --producer-manifest "$MOON_PRODUCER" \
  --input "$MOON_CASE_ROOT/prepared/native.in" \
  --mode run --mpi-ranks 8 --wall-seconds 600 \
  --output-dir "$MOON_CASE_ROOT/native-owner"
```

The JSON contains explicit native scientific controls; it is input to the case
producer, not an AMPS application input. The `.in` application is a complete
normal receipt. Site/data paths must exist. No source threshold, physical mass,
rate or missing model parameter is manufactured by the preparer. Changing the
application or controls creates a new experiment and requires new retained
input identities. All output directories must be new.

Expect logs/resources, executable/build/source/input identities, native input
receipt, every original parent, per-rank event and exact hexadecimal ledger CSVs
under `native-owner/`. Native PASS means the finite execution completed; the
independent checks are separate. Pin its original `run_manifest.json` and
`result.json` in the binding described in
[P02_README.md](test/sodium_cohort/P02_README.md), then audit from the active tree:

```sh
/usr/bin/python3 -B "$MOON_ACTIVE_ROOT/srcMoon/test/sodium_cohort/audit_commissioning.py" \
  --binding /absolute/path/to/pinned-owner-binding.json \
  --output-dir "$MOON_CASE_ROOT/independent-commissioning"
/usr/bin/python3 -B "$MOON_ACTIVE_ROOT/srcMoon/test/sodium_cohort/verify_p02.py" \
  --binding /absolute/path/to/pinned-owner-binding.json \
  --output-dir "$MOON_CASE_ROOT/independent-p02"
```

The commissioning audit retains conservation/rate/clock criteria and comparison
CSV/PNG/EPS, including per-event ledger plots. P02 retains per-state mechanical
criteria in CSV/JSON and an audit log. The survival curve is an analytical
source-law comparison, with no observational data. Refinement and statistical
qualification remain NOT_EVALUATED; formal M3 remains NOT QUALIFIED. Unsupported
physical source/state/bath/surface selections keep their ordinary data gates.

### I21 matched timestep and MPI development cohort

Use the same retained source-bound producer and a hash-pinned commissioning
geometry receipt to freeze and run a six-launch cohort. The exact preparation
and `run_all_test.py --only I21` commands are in the [I21 runner README](test/integration-tests/I21/README.md).
Look under `I21/homogeneous/` or `I21/trajectory/` in the chosen output directory for the native
run manifests, parent and event records, and CSV/PNG/EPS timing comparisons.
The outer result retains formal I21 status and reports any supporting test
failure. An ordinary aggregate run without a cohort plan stays lightweight.

### I23 stationary helium pilot

Use the profile and exact hash from the atomic-family preparation. From the
AMPS root the [I23 runner](test/integration-tests/I23/README.md) compiles the
helium reservoir provider and writes coefficient, conditional-moment and
admission-control comparisons under the chosen output directory:

```bash
/usr/bin/python3 -B srcMoon/test/integration-tests/I23/run_helium_pilot.py \
  --profile /path/to/profile.json --profile-sha256 PROFILE_SHA256 \
  --output-dir test_output/srcMoon/i23-helium-NEW
```

See `result.json`, `comparison.csv`, the exact profile snapshot and the
separate build/native logs. This controlled run has no measured plasma input.

### I23 prescribed electron-bath coupling

This runs the approved controlled 10 eV bath through the production owner.
The ensemble needs only the profile and D05. The PIC matrix additionally
needs a source-bound Na/He/He⁺ producer; build it as described in the
[I23 README](test/integration-tests/I23/README.md#prescribed-electron-bath-coupling-2026-11-25).

```bash
producer=test_output/srcMoon/i23-bath-uGAYf6qI/helium-producer/normal_build_evidence.json
out=$(mktemp -d -p "$PWD/test_output/srcMoon" i23-formal-XXXXXXXX)
/usr/bin/python3 -B srcMoon/test/run_all_test.py --only I23 --skip-build \
  --campaign-mode regression --mpi-ranks 8 --test-runtime-seconds 1800 \
  --i23-helium-profile test_output/srcMoon/m3-followup-MfThQfny/atomic-family/profile.json \
  --i23-pic-producer "$producer" --i23-pic-refinement --output-dir "$out"
```

`$out/I23/result.json` lists both supporting children and keeps the formal
status SKIPPED. Open `$out/I23/helium-electron-bath/comparison.png` for
survival against `exp(−n_e k t)`, and
`$out/I23/helium-pic-electron-bath/step_populations.png` for the represented
He and He⁺ after each step. The producer lives on temporary storage; after a
reboot, rebuild it and pass the new manifest.

### I24 sodium charge exchange at 4.0 keV/u

This development runner uses the pinned r14 partial cross section for a
monoenergetic proton beam incident on ground-state Na. The ensemble
estimates removal frequency and lifetime. The independent linked test
tracks both outgoing particles through each production force-cut-surface
mover. Prepare a three-species build by following the
[source-qualification instructions](test/source_qualification_build/README.md#build-the-i24-sodium-capture-case).

```bash
db=/data/vtenishe/CHEMICAL_DATABASE/snapshots/amps-chemistry-20261124-r14
producer=test_output/srcMoon/m3-production-aPrmWNe5/closure-capture-build/normal_build_evidence.json
out=$(mktemp -d -p "$PWD/test_output/srcMoon" i24-capture-XXXXXXXX)
/usr/bin/python3 -B srcMoon/test/run_all_test.py --only I24 --skip-build \
  --campaign-mode regression --mpi-ranks 8 --test-runtime-seconds 1800 \
  --i24-sodium-lock "$db/database-lock.json" --i24-pic-producer "$producer" \
  --output-dir "$out"
```

On another machine set `producer` to the new local producer manifest.
Supply `--i24-mass-inventory` when the NIST source packet is in another
directory. Add `--i24-pic-refinement` to the command to compare one, two and four
common time steps and three reference-weight scales on both MPI layouts.
The aggregate retains each child status in `$out/I24/result.json`.
Open `$out/I24/sodium-ground-capture/comparison.png` for the population
comparison and
`$out/I24/sodium-pic-ground-capture/second-8mpi/continuation.png` for the
two products' ballistic continuations. Corresponding CSV and EPS files are
available. Each run uses a fresh output directory. Formal I24 has separate
predecessor and physical-coverage requirements.
[Physics and output contract](test/integration-tests/I24/SODIUM_CASE.md)

### Run a frozen family of M3 development cases

The family runner sequences the existing illuminated, dark and moving Na
photoionization comparisons and the Na charge-exchange ensemble and
particle matrix. It uses one development-wide deadline and keeps each
case's original native logs, signed ledgers and independent figures.
Build or retain the current-source Na/Na⁺ and Na/Na⁺/H SPICE producers,
then prepare and execute the pinned JSON plan using the
[chemistry family commands](test/qualification/README.md#frozen-chemistry-family).
Read `result.json` and `comparison.csv` under the selected output directory.
The formal I21–I24 decisions and their predecessor receipts are under
`formal/`.
