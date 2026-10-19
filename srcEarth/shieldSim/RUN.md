# Building, configuring, testing, and running shieldSim

This is the operational guide for `shieldSim`. The scientific scope and
validation status are summarized in `README.md`; the normative staged plan is
`development_and_validatoin.md`; input examples are explained in
`examples/README.md`. Commands below assume the current directory is the
`shieldSim` source directory.

`shieldSim` is an environment-effects model for spacecraft-anomaly analysis and
preliminary mission assessment. It is not an instrument-design, device-
qualification, organ-dose, or crew-risk tool. A successful run is not evidence
that its requested quantity is physically qualified: always inspect the
quantity-level `validation_status`, `accuracy_status`, limitations, and
uncertainty fields in the result.

## 1. Software prerequisites

Required to compile and run the model:

- a C++17 compiler;
- CMake 3.16 or newer;
- a Geant4 installation with CMake package files, shared libraries, and the
  data sets required by the selected physics list;
- Python 3 for test drivers and response-weight construction; and
- Git if a real code revision and clean/dirty state are to be embedded in
  output provenance.

`matplotlib` is additionally required for validation figures. The core model
does not require it. No network access is needed for normal builds, tests, or
runs because validation inputs are checksum-pinned locally. The larger CRaTER
package is intentionally outside the source repository at:

```bash
export SHIELDSIM_REFERENCE_DATA_DIR=/data/vtenishe/shieldSim_reference_data
test -f "${SHIELDSIM_REFERENCE_DATA_DIR}/reference_crater/manifest.sha256"
```

The tested local Geant4 configuration can be selected explicitly as follows:

```bash
export Geant4_DIR=/data/vtenishe/G4/install/lib/cmake/Geant4
```

For another installation, set `Geant4_DIR` to the directory containing
`Geant4Config.cmake`. Do not set it to the executable or library directory.
If the Geant4 installation provides an environment script, source it before
configuring and again in every new shell used for direct model runs:

```bash
source /path/to/geant4/bin/geant4.sh
```

That script normally supplies the shared-library path and Geant4 data-set
variables. High-precision neutron configurations such as `FTFP_BERT_HP` need a
valid `G4NEUTRONHPDATA`; shieldSim tests intentionally fail if a required data
set is unavailable rather than silently reducing the physics configuration.

Useful checks are:

```bash
test -f "${Geant4_DIR}/Geant4Config.cmake"
cmake --version
c++ --version
python3 --version
python3 -c 'import matplotlib; print(matplotlib.__version__)'
```

## 2. Configure and compile

Build artifacts and generated tests belong below `test_output/`, not in the
source directories or `/tmp`. A release build with all registered tests is:

```bash
cmake -S . -B test_output/build \
  -DCMAKE_BUILD_TYPE=Release \
  -DBUILD_TESTING=ON \
  -DGeant4_DIR="${Geant4_DIR}" \
  -DSHIELDSIM_CRATER_REFERENCE_DIR="${SHIELDSIM_REFERENCE_DATA_DIR}/reference_crater"
cmake --build test_output/build -j 8
```

Choose a job count appropriate for the machine. To see the exact compiler and
link commands, rebuild with:

```bash
cmake --build test_output/build --verbose
```

After changing CMake options, rerun the configure command. After changing only
C++ sources, `cmake --build` is sufficient. CMake generates
`test_output/build/generated/BuildProvenance.hh`; results record the revision,
tree state, compiler, flags, platform, build type, Geant4 version, and relevant
data-library paths. Configure only after the intended source changes are
present so the recorded clean/dirty state is meaningful.

The build creates these executables:

- `shieldSim`: case validation and Geant4 transport;
- `shieldSimSourceTool`: strict conversion of pinned AMS-02, ACE/CRIS, and
  GOES-16 source products;
- `shieldSimResponseTool`: response-matrix assembly and folding;
- `shieldSimExposureTool`: Stage-7 time integration, spectrum-response
  quadrature, and event allocation; and
- `shieldSimReferenceTool`: controlled Geant4 reference-table probes used by
  validation.

Do not copy only an executable to another machine. It still depends on the
matching Geant4 shared libraries and data sets, and its provenance refers to
the build in which it was created.

## 3. Runtime environment and paths

Direct Geant4 runs must use the same runtime environment as the build. In a
new shell, source the Geant4 environment again and confirm that the executable
loads:

```bash
source /path/to/geant4/bin/geant4.sh
test_output/build/shieldSim --help
```

If the dynamic loader reports a missing `libG4*.so`, the Geant4 library
directory is absent from `LD_LIBRARY_PATH` (or the platform equivalent). Fix
the installation environment; do not copy individual Geant4 libraries into
the shieldSim tree.

Paths inside a versioned case are resolved relative to the case file where the
schema defines a file reference, while relative output prefixes are written
relative to the process working directory. The safest workflow is to run from
an owned output directory while referencing a case by an explicit path:

```bash
mkdir -p test_output/manual_stage4
cd test_output/manual_stage4
../build/shieldSim --validate-case=../../examples/stage4_tid_case.json
../build/shieldSim --case=../../examples/stage4_tid_case.json
cd ../..
```

The case's `output.prefix` is used for filenames. When a case contains several
absorbers, shieldSim runs each absorber in an independent Geant4 world and
appends `__<absorber_id>` to the prefix. An absorber therefore never becomes
unintended shielding for another absorber.

## 4. Preparing a versioned input case

Use versioned JSON cases for all new scientific work. Start from the example
closest to the requested observable:

- `examples/stage1_case.json`: explicit multilayer/slab and case-v1 contract;
- `examples/stage2_ams02_proton_validation_case.json`: downloaded AMS A=1
  proton component in its narrow verification-only admission;
- `examples/stage3_response_cell_case.json`: monoenergetic response cell;
- `examples/stage4_tid_case.json`: finite material-slab TID;
- `examples/stage5_ddd_neq_case.json`: natural-silicon DDD and neutron-only
  ASTM `n_eq`;
- `examples/stage6_let_hardness_case.json`: electronic LET and primary-proton
  `H100/10`; and
- `examples/stage7_time_resolved_exposure.json`: explicit time intervals,
  gaps, configuration changes, and incomplete uncertainty handling.

A transport case must explicitly define these sections:

1. `schema_version`, a unique `case_id`, and `purpose`;
2. material-library identity/path and any immutable case overrides;
3. world, source plane, ordered shield stack, and independent absorbers;
4. source mode and spectrum/response-node definition;
5. requested quantities;
6. accuracy scope, target, coverage, uncertainty budget, and sampling limits;
7. physics list, EM/ion choices, production cut, maximum step, and random seed;
8. output prefix.

The parser rejects unknown fields and units. This is deliberate: a misspelled
scientific option must not turn into a hidden default.

### 4.1 Geometry and shielding

All dimensions and placements belong in the input. A shield has one through
five ordered layers. Each layer may be supplied as geometric thickness or
areal density; the latter is converted using the resolved material density.
The resulting geometric thickness must be no more than 10 cm per layer.
Layer order is physical and changes the geometry/response hashes.

Case v2 surfaces use a translation plus scalar-first `(w,x,y,z)` unit
quaternion mapping local coordinates to world coordinates. For slabs, local
`+z` is downstream and the absorber transform origin is the upstream face.
The current absorber provider is `slab`; all three slab dimensions must be
explicit. The provider/scoring boundary is designed for a future arbitrary-
shape implementation, but no arbitrary shape is currently admitted for a
science-facing result.

### 4.2 Materials and overrides

The default versioned material library is
`data/material_library.v1.json`. A case names both its library ID/version and
path. Each material expands to explicit density/composition/provenance and a
content hash. Case overrides are allowed, but must provide the schema-required
physical properties and create a different resolved hash. A generic commercial
or screenshot label is not a unique material unless its properties have been
supplied.

### 4.3 Sources

Case v2 supports exact monoenergetic response nodes and normalized discrete
response nodes. The narrow Stage-2 AMS fixture is admitted only for
`purpose: verification_test`, the A=1 proton component, and its declared
isotropic boundary. It is not a full GCR model.

The legacy three-column spectrum interface requires its metadata sidecar and
is retained for compatibility. Do not route published AMS, ACE/CRIS, or GOES
data through it. GOES SGPS directions are ingestion/time-integration data only;
transport remains disabled until the sensor-to-simulation transform and
angular response are supplied.

### 4.4 Quantities and accuracy

Only request quantities supported by the selected case schema, species,
absorber material, response package, and validated energy domain. The default
requested relative accuracy is 0.10 (`k=1`). The request may be overridden per
quantity. It does not establish accuracy: an incomplete physical/systematic
budget yields `not_assessable` even when Monte Carlo precision is good.

The enabled domains are intentionally bounded:

- TID: finite Al, natural-Si, or water slab for proton/He-4 incidence;
- DDD: natural silicon with registered proton, He-4, or neutron response;
- `n_eq`: natural silicon, ASTM E722 neutron contribution only;
- LET: registered electronic-stopping/material/species domain; and
- `H100/10`: forward primary-proton shield-exit fluence ratio.

The result, not the request, is authoritative. Inspect `validation_status`,
`accuracy_status`, uncertainty components, intervals, rejected/out-of-domain
counts, and limitations for every output.

### 4.5 Validate before transport

Always parse and resolve a new or edited case before spending histories:

```bash
mkdir -p test_output/case_check
cd test_output/case_check
../build/shieldSim --validate-case=../../examples/stage6_let_hardness_case.json
cd ../..
```

Validation writes the resolved case and hashes but does not start transport.
Then run the identical file:

```bash
mkdir -p test_output/stage6_example
cd test_output/stage6_example
../build/shieldSim --case=../../examples/stage6_let_hardness_case.json
cd ../..
```

Do not combine `--case` or `--validate-case` with legacy physical CLI options.
To inspect compatibility options and their units:

```bash
test_output/build/shieldSim --help
test_output/build/shieldSim --list-materials
test_output/build/shieldSim --list-target-materials
test_output/build/shieldSim --list-quantities
```

## 5. Published-source adapters

The source tool validates strict metadata and writes a canonical versioned
product. Example commands using committed fixtures are:

```bash
mkdir -p test_output/source_adapters

test_output/build/shieldSimSourceTool ams02 \
  --info=tests/reference_spectra/gcr/ams02_20131118_hapi_info.json \
  --data=tests/reference_spectra/gcr/ams02_20131118_hapi_data.csv \
  --interpretation=tests/reference_spectra/gcr/ams02_20131118_interpretation.v1.json \
  --output=test_output/source_adapters/ams02.json

test_output/build/shieldSimSourceTool ams02-cycle \
  --interpretation=tests/reference_spectra/gcr/ams02_2025/interpretation.v1.json \
  --bartels-rotation=2581 \
  --output=test_output/source_adapters/ams02_cycle_br2581.json

test_output/build/shieldSimSourceTool ace-cris \
  --info=tests/reference_spectra/gcr/ace_cris/ac_h3_cris_hapi_info.json \
  --data=tests/reference_spectra/gcr/ace_cris/ac_h3_cris_20191201_all.csv \
  --interpretation=tests/reference_spectra/gcr/ace_cris/interpretation.v1.json \
  --output=test_output/source_adapters/ace_cris.json

test_output/build/shieldSimSourceTool goes16 \
  --csv=tests/reference_spectra/sep/goes16_sgps_gle72_timeseries.csv \
  --definition=tests/reference_spectra/sep/goes16_gle72_event_definition.v1.json \
  --output=test_output/source_adapters/goes16_gle72.json
```

These commands establish ingestion and unit/time semantics. They do not by
themselves validate transport or authorize a complete mission environment.

## 6. Response matrices

Response cells are produced by running case-v2
`monoenergetic_response_cell` cases. Every cell assembled into one matrix must
have the same geometry, materials, absorber, scoring contract, physics/build
identity, cuts, and response-package identity. Energy and direction identify
the cell and may differ. The response key enforces this boundary.

Assemble full cells with:

```bash
test_output/build/shieldSimResponseTool assemble \
  test_output/response/matrix.json \
  test_output/response/cell_001_response_cell.json \
  test_output/response/cell_002_response_cell.json
```

For a large qualification matrix, retain only explicitly registered primitive
scores:

```bash
test_output/build/shieldSimResponseTool assemble-selected \
  test_output/response/matrix_tid.json \
  absorber/silicon_slab/ionizing_deposited_energy_MeV,shield_exit/primary/proton/crossing_count \
  test_output/response/cell_001_response_cell.json \
  test_output/response/cell_002_response_cell.json
```

Selected assembly still fully validates every input cell before discarding
unrequested scores/covariances. It fails if any selected score is absent.

Construct explicit interpolation/quadrature weights, then fold:

```bash
python3 scripts/build_response_weights.py \
  test_output/response/matrix_tid.json \
  test_output/response/source_nodes.json \
  test_output/response/weights.json

test_output/build/shieldSimResponseTool fold \
  test_output/response/matrix_tid.json \
  test_output/response/weights.json \
  test_output/response/fold.json
```

The matrix hash must match the weight file. Extrapolation is prohibited.
Nonlinear outputs must fold their primitive numerator and denominator scores,
not a precomputed ratio.

## 7. Stage-7 exposure, quadrature, and event allocation

Integrate explicit UTC intervals with:

```bash
mkdir -p test_output/stage7_manual
test_output/build/shieldSimExposureTool integrate \
  examples/stage7_time_resolved_exposure.json \
  test_output/stage7_manual/exposure.json
```

The input must use strict `YYYY-MM-DDTHH:MM:SSZ` timestamps. Intervals are
half-open, ordered, and contiguous in the ledger. Missing time must appear as
an explicit `quality_state: gap` interval; the tool never fills or interpolates
a hidden gap. Configuration changes remain separate segments. Anomaly
timestamps are retained as markers and do not alter integration boundaries.
Covariance `rate_unit_squared` values use explicit parentheses, for example
`(Gy/s)^2`; `Gy/s^2` is rejected because it is not the square of a rate unit.

The same tool provides two deterministic supporting operations:

```bash
test_output/build/shieldSimExposureTool fold-spectrum \
  test_output/stage7_manual/quadrature_input.json \
  test_output/stage7_manual/quadrature_output.json

test_output/build/shieldSimExposureTool allocate \
  test_output/stage7_manual/allocation_input.json \
  test_output/stage7_manual/allocation_output.json
```

`fold-spectrum` uses response interpolation linear in log energy, adaptive
quadrature, and no extrapolation. `allocate` performs pilot-based stratified
Neyman history allocation with a required per-stratum minimum. It reallocates
unbiased histories; it does not enable Geant4 physics biasing. If a rare tail
is unresolved, the affected output remains `not_assessable`.

The committed AMS A=1 proton spectrum is sufficient for the Stage-7 component
test but remains an incomplete GCR environment. Published source covariance is
not available for the fold, so the physical source uncertainty is explicitly
incomplete. GOES transport remains disabled pending its boundary mapping.

## 8. Running tests

CTest is the authoritative entry point. With CMake 3.16, change into the build
directory rather than using newer `ctest --test-dir` syntax:

```bash
cd test_output/build
ctest -N
ctest --output-on-failure
cd ../..
```

Generated inputs, logs, results, matrices, evidence, temporary files, and
figures remain under `test_output/build/`. They are test artifacts, not pinned
reference data.

Run one stage or one test with a regular expression:

```bash
cd test_output/build
ctest --output-on-failure -R '^stage7_'
ctest --output-on-failure -R '^stage7_published_qualification$'
ctest --output-on-failure -R '^stage8_reference_inputs$'
ctest --output-on-failure -R '^stage4_(core|nist_reference|tid_transport)$'
cd ../..
```

### Parallel Stage-7 transport

One `shieldSim` case is currently serial. Do not launch the executable itself
with `mpiexec`: shieldSim has no MPI rank partitioning or rank-aware output,
so every rank would repeat the same histories and write the same files.

The Stage-7 response campaign instead supports bounded **case-level**
parallelism. Every worker launches one independent serial shieldSim process
with its registered seed and history count in a unique directory. Completion
order cannot alter pooling or evidence order. Set `SHIELDSIM_STAGE7_WORKERS`
to a positive integer no larger than the CPUs available to the process. The
default is `1`; a conservative starting value on a workstation is `4` or `8`.

This is a complete working example from the shieldSim source directory:

```bash
cmake -S . -B test_output/build \
  -DGeant4_DIR=/data/vtenishe/G4/install/lib/cmake/Geant4
cmake --build test_output/build --parallel 8
cd test_output/build
SHIELDSIM_STAGE7_WORKERS=8 \
  ctest --output-on-failure -R '^stage7_published_qualification$'
cd ../..
```

The build's Geant4 data environment must still be available as described in
Section 3. `ctest -j` does not accelerate this particular gate because it is a
single CTest test; the environment variable controls its internal queue of
independent cases. Avoid simultaneously running several transport-heavy CTest
gates with both `ctest -j` and a large Stage-7 worker count, since their product
can oversubscribe CPUs and memory.

An interrupted run is restartable. Rerun the identical command: the driver
hashes its executables, plan, source fixtures, material library, and controlling
script, then reuses only a case whose canonical input, requested history count,
result, and response cell all match. A changed campaign identity invalidates
normal reuse. Worker count is execution scheduling, not physics input, and may
change between restarts without changing a valid case result.

For response-only cases, restart retention is intentionally compact. The C++
response tool first parses the complete generated cell, including its full
covariance triangle. The driver then retains the two Stage-7 registered scores,
their three covariance entries, and a hash-checked completion marker, and
removes the redundant full result. Published-source and direct-quadrature
results remain complete. A missing, malformed, wrong-history, or hash-mismatched
compact artifact is rerun automatically.

Case-level parallelism deliberately does not use Geant4 event multithreading.
Moving one case to `G4MTRunManager` would require thread-local actions, tally
merging, deterministic random-stream policy, output serialization, and new
serial-versus-threaded validation. MPI would be useful later for distributing
the same independent-case queue across compute nodes, but is not required for
single-node execution.

The Stage-7 suite includes deterministic analytic time/quadrature/covariance
tests, strict rejection boundaries, five-seed published-component transport,
both Al/HDPE layer orders, natural-silicon and water slabs, direct-versus-
response-matrix comparisons, response-grid refinement, and exact 27-day
integration. Each direct row uses five 8,000-history seeds. Each response stack
uses independent 256-history pilots for a frozen 73,728-history base allocation
and a 56,320-history interleaved energy-refinement allocation, followed by a
112,640-history exact-node supplement, for each of five independent response
replicates. Every allocation has a 512-history floor per cell. The combined
45-energy-node campaign contains 1,213,440 production histories per stack. Its
reported response error is the larger of the propagated within-replicate and
empirical between-replicate standard errors. Its passing status is deliberately
labeled component/framework
qualification rather than full physical GCR/SEP qualification.

The equal-history direct seed ensembles use the same conservative
maximum-within/between standard-error rule. This is important for the rare TID
tail: at the exact 45-node grid, the reversed-stack empirical between-seed
error was 20.8% larger than the pooled event estimate. Neither comparison mean
nor the registered 10% and two-sided 95% acceptance limits is changed.

The final `validation_figures` test runs only after its numerical dependencies.
Figures are review aids; JSON evidence and numeric assertions are authoritative.
After a successful full suite, inspect:

```text
test_output/build/validation_evidence/
test_output/build/validation_figures/
test_output/build/stage7_published_runs/
```

The older layered shell scripts in `tests/` remain useful compatibility and
regression drivers. Use the CTest suite for staged qualification because it
also executes the strict analytic/reference gates and produces the unified
figure manifest.

## 9. Foreground-only axial CRaTER benchmark

This validation fixture is a separate executable, not a production shieldSim
geometry. Build it with the rest of the project or explicitly with:

```bash
cmake --build test_output/build --parallel 8 \
  --target shieldSimCraterBenchmark stage8_crater_core_tests
```

The example pins the first CRaTER observation interval and its time-matched
BON2020 H--Ni spectrum. From the repository root, a serial run is:

```bash
source /path/to/geant4/bin/geant4.sh

test_output/build/shieldSimCraterBenchmark \
  examples/stage8_crater_foreground_axial.json

python3 scripts/plot_crater_foreground_benchmark.py \
  test_output/crater_foreground_axial/crater_bon2020_2009259_2010065__foreground_only_axial_result.json \
  test_output/crater_foreground_axial/crater_bon2020_2009259_2010065__foreground_only_axial_comparison
```

The distributed example records the requested installation path
`/data/vtenishe/shieldSim_reference_data/reference_crater`. If that package is
moved again, update the case's `source.path` and configure
`SHIELDSIM_CRATER_REFERENCE_DIR` to the same directory. CTest passes the CMake
path explicitly and does not depend on the source package being inside Git.

The case's `histories` is the global history count and its accuracy field is
the requested relative Monte Carlo uncertainty (10% in the example). The
executable reports whether every nonzero bin from 0.5 to 10 keV/um reached that
target; it never changes the event count or hides under-resolved bins. The JSON
is authoritative. The CSV and PNG/SVG files are derived review artifacts.

### Local process parallelism

Use the supplied launcher on a workstation. This working four-process example
runs four independent Geant4 processes concurrently, checks their exit status,
merges their raw moments, and creates the comparison figure:

```bash
python3 scripts/run_crater_parallel.py \
  --benchmark test_output/build/shieldSimCraterBenchmark \
  --case examples/stage8_crater_foreground_axial.json \
  --workers 4
```

The launcher inherits the current shell environment, so load `geant4.sh`
first as shown in the serial example. An immediate
`libG4*.so: cannot open shared object file` error means the Geant4 runtime
environment was not loaded; no histories have run in that case. Use no more
workers than available physical CPU capacity and memory: every worker is a
complete Geant4 process with its own geometry, physics tables, and log.

The event histories are independent, so process sharding reduces the transport
part of elapsed time without changing the global requested history count. It
does not make initialization, final merging, or plotting parallel. On the
development host, the same 560-history smoke workload took 109.21 s serially
and 65.38--68.14 s with two workers (about 1.6--1.7x faster). These values are
an example, not a runtime guarantee; speedup depends on ion mix, secondary
production, storage, CPU count, and concurrent system load.

The launcher retains one log, partial JSON, and partial CSV per shard under
`test_output/crater_foreground_axial/`. It then writes the same unsuffixed
merged result/spectrum names used by the serial command. A rerun may overwrite
files for the same case ID, so archive a completed comparison before changing
the case or code.

### MPI launcher example

`shieldSimCraterBenchmark` does not link MPI. Instead it reads the standard
Open MPI `OMPI_COMM_WORLD_RANK/SIZE` variables or the MPICH/PMI
`PMI_RANK/SIZE` variables. This avoids an MPI build dependency while allowing
the same executable to be distributed across nodes:

```bash
mpiexec -n 4 test_output/build/shieldSimCraterBenchmark \
  examples/stage8_crater_foreground_axial.json

python3 scripts/merge_crater_shards.py \
  --output test_output/crater_foreground_axial/crater_bon2020_2009259_2010065__foreground_only_axial_result.json \
  test_output/crater_foreground_axial/crater_bon2020_2009259_2010065__foreground_only_axial__shard_*_partial_result.json

python3 scripts/plot_crater_foreground_benchmark.py \
  test_output/crater_foreground_axial/crater_bon2020_2009259_2010065__foreground_only_axial_result.json \
  test_output/crater_foreground_axial/crater_bon2020_2009259_2010065__foreground_only_axial_comparison
```

The merge tool refuses missing/duplicate ranks, different cases, source or
observation hashes, geometry, physics, Geant4/code identity, bin definitions,
or history totals. It also requires distinct deterministic rank seeds. It adds
the raw importance-weight sums and sums of squares first, then recomputes the
global rate and standard error. Averaging shard spectra or their error bars
would be statistically wrong and is not supported.

Explicit shard arguments are available for schedulers that do not export a
recognized MPI environment:

```bash
test_output/build/shieldSimCraterBenchmark \
  examples/stage8_crater_foreground_axial.json \
  --shard-index=0 --shard-count=4
```

Run the command once for each index `0` through `3`, then use the same merge
and plot commands. Histories are divided by quotient/remainder so their sum is
exactly the case's global count. Rank `r` uses seed `random_seed + r`; these are
distinct deterministic Geant4 pseudo-random streams. This is not a mathematical
claim that two finite pseudo-random streams can never share a value.

Do not run `mpiexec -n N` with an older executable lacking shard-aware output:
that launches duplicate full simulations and can collide on output files. The
current result must say `complete_parallel_merge` and
`complete_all_shards_present` before it is used. The routine CTest gate uses
two local workers and verifies this complete merge boundary.

Every result and figure is labeled `foreground_only_axial`. It includes the
published axial stack and D1--D6 acceptance but excludes housing, boards, dead
layers, thermal blanket, spacecraft, lunar albedo, and out-of-aperture/side
background. It is a diagnostic comparison, not full-instrument validation or
Stage-8 sign-off.

## 10. Stage-9 verification interface

Stage 9 is implemented as a **verification-only** interface while earlier
scientific gates remain open. The capability registry
`data/stage9_capabilities.v1.json` is authoritative for enabled profiles,
quantities, badges, and disabled reasons. A browser label cannot enable
spacecraft-anomaly or preliminary-mission purposes, combine response contracts
from different profiles, select a template material, introduce a non-slab
shape, or select the Stage-11 calculated-spectrum source. The latter is
available through direct versioned case JSON, but is not enabled in the
Stage-9 verification UI registry.

### 10.1 Request and preflight

Start with `examples/stage9_ui_request.json`. Its top-level fields are:

- `profile_id`: one bounded TID, silicon-damage, or LET/hardness verification
  source/response contract;
- `accuracy_target`: relative fraction; when omitted the compiler records the
  approved 0.10 default;
- `shield_configurations`: independent stacks, each with full transverse size,
  transform, and one through five source-to-absorber layers;
- `absorbers`: independent finite slabs with material, X/Y sizes, thickness,
  and transform; and
- `runtime_calibration`: either `not_provided`, or a measured
  `seconds_per_history` with provenance.

Layer thickness is either `geometric` in `mm`/`cm` or `areal_density` in
`g/cm2`. The compiler resolves density through the versioned material library,
reports its uncertainty metadata, calculates the other thickness form, and
enforces the inclusive 100 mm derived limit. A fully documented
`material_overrides` record can change density/composition in the request; the
C++ material parser performs final provenance, composition, role, and status
validation. Templates such as unresolved lunar/Mars regolith remain disabled.

Run a read-only preflight:

```bash
python3 scripts/stage9_interface.py preview \
  --request examples/stage9_ui_request.json
```

The summary separately reports:

```text
transport_case_count = number of shields * number of absorbers
requested_quantity_series_count = transport_case_count * number of quantities
```

Minimum/maximum histories are exact. Wall time is `not_assessable` without a
measured throughput; the compiler never invents one. Combined uncertainty is
also `not_assessable` when its systematic budget is incomplete. A declared
systematic coverage floor at or above the target is reported unachievable
before submission.

### 10.2 Compile, validate, and run without the browser

Create a new bundle only in an empty/nonexistent directory:

```bash
python3 scripts/stage9_interface.py compile \
  --request examples/stage9_ui_request.json \
  --output test_output/stage9_ui_jobs/example

for case_file in test_output/stage9_ui_jobs/example/cases/*.json; do
  test_output/build/shieldSim --validate-case="${case_file}"
done

python3 scripts/stage9_interface.py run \
  --bundle test_output/stage9_ui_jobs/example \
  --executable test_output/build/shieldSim
```

The bundle contains `request.canonical.json`, generated v2 cases,
`manifest.json`, `job_state.json`, and one `runs/<case_id>/` directory per
generated shield case. The manifest embeds the complete UI request, content
hashes, exact control-to-case JSON paths, every shield/absorber matrix mapping,
seed rule, material/layer resolution, cost/uncertainty preflight, validation
badges, and disabled reasons. Each run directory receives
`stage9_provenance.json` before transport, tying its outputs to the manifest,
case hash, matrix cells, and non-production validation badges. A subprocess
exit code of zero is not sufficient for job completion: the runner requires
every expected absorber-specific resolved-case/result pair and verifies their
case/absorber/hash linkage before checkpointing that generated case.

To request cancellation from a second terminal:

```bash
python3 scripts/stage9_interface.py cancel \
  --bundle test_output/stage9_ui_jobs/example
```

The runner terminates an active subprocess and records
`interrupted_restart_required`. It checkpoints only completed generated case
files; it does not claim to checkpoint Geant4 histories. Resume skips completed
cases and restarts an interrupted case from its declared seed:

```bash
python3 scripts/stage9_interface.py run \
  --bundle test_output/stage9_ui_jobs/example \
  --executable test_output/build/shieldSim \
  --resume
```

### 10.3 Local browser

The server uses only the Python standard library and binds to loopback by
default:

```bash
python3 scripts/stage9_ui_server.py \
  --host 127.0.0.1 --port 8765 \
  --output-root test_output/stage9_ui_jobs \
  --executable test_output/build/shieldSim
```

Open `http://127.0.0.1:8765/`. The layer editor displays density, density
uncertainty, mass thickness, derived geometric thickness, and order. The form
also exposes explicit shield/absorber dimensions and positions. Non-identity
unit quaternions and full material overrides can be entered in the advanced
versioned JSON view, then round-tripped back into the controls. `Preview`,
`Compile`, `Run`, `Cancel`, `Resume`, and `Refresh status` call the same core
compiler/runner used above. When `--executable` is supplied, `Compile` also
passes every generated case through the production C++ `--validate-case`
parser and writes `compile_validation.json`; `Run` is unavailable without that
server option. Do not expose this development server to an
untrusted network; it has no authentication or TLS.

### 10.4 Tests and qualification boundary

After configuring with `BUILD_TESTING=ON`:

```bash
cd test_output/build
ctest --output-on-failure -R '^stage9_interface$'
cd ../..
```

The test covers UI-control presence, JavaScript syntax when Node is available,
serialization round-trip, exact matrix cardinality, complete UI-to-case
traceability, deterministic hashes/seeds, material override and thickness
conversion, every major disabled guard, C++ validation of generated cases,
provenance badges, active cancellation, case-boundary cancellation, and
resume. `tests/reference/stage9_candidate_baseline.v1.json` is explicitly
`candidate_unreviewed_not_acceptance`; it protects V6 software behavior but is
not a reviewed regression acceptance baseline. Stage 9 cannot be marked
complete until its entry gates pass and that review exists.

## 11. Stage-10 analytic and tessellated-mesh absorbers

Stage 10 currently provides `analytic_csg.v1` and `tessellated_mesh.v1` for
verification work. Neither is enabled for spacecraft-anomaly or
preliminary-mission purposes because the independent nontrivial-shape
validation gate is open. Do not change a case's purpose merely to bypass this
boundary.

Build it with the normal executable; no optional library is required beyond
Geant4:

```bash
source /data/vtenishe/G4/install/bin/geant4.sh
cmake -S . -B test_output/build \
  -DGeant4_DIR=/data/vtenishe/G4/install/lib/cmake/Geant4 \
  -DCMAKE_BUILD_TYPE=Release \
  -DBUILD_TESTING=ON
cmake --build test_output/build --target shieldSim -j 8
```

Validate the complete Boolean-enclosure example without running transport:

```bash
test_output/build/shieldSim \
  --validate-case=examples/stage10_analytic_csg_case.json
```

Run it in a controlled output directory. Input-relative material paths resolve
relative to the case file, while result files are written under the process
working directory:

```bash
mkdir -p test_output/stage10_example
cd test_output/stage10_example
../build/shieldSim --case=../../examples/stage10_analytic_csg_case.json
cd ../..
```

The example is a synthetic 30×30×20 mm aluminum box shell made by subtracting
a concentric 26×26×16 mm void, with a separate 5×5×1 mm silicon sensitive box
in the cavity. It verifies region/material/Boolean plumbing only; it is not an
authoritative spacecraft or electronics assembly.

The directly runnable shape gallery contains three analytic and three mesh
cases:

```text
examples/stage10_shapes/analytic_box.json
examples/stage10_shapes/analytic_cylinder.json
examples/stage10_shapes/analytic_sphere.json
examples/stage10_shapes/mesh_box.json
examples/stage10_shapes/mesh_tetrahedron.json
examples/stage10_shapes/mesh_octahedron.json
```

The complementary selection gallery under `examples/configuration_cases/`
contains seven complete cases in which source species, absorber material,
shape, shield footprint, layer count/material/thickness, and calculation list
actually vary. Its detailed matrix and editable JSON fragments are in
`examples/configuration_cases/README.md`. Validate the entire set with:

```bash
for case_file in examples/configuration_cases/*.json; do
  test_output/build/shieldSim --validate-case="$case_file" || exit 1
done
```

For example, run the 1.05 MeV neutron/analytic-cylinder silicon DDD and `n_eq`
case entirely below `test_output/`:

```bash
mkdir -p test_output/configuration_cases/02_si_cylinder_neutron_damage
cd test_output/configuration_cases/02_si_cylinder_neutron_damage
../../build/shieldSim \
  --case=../../../examples/configuration_cases/02_si_cylinder_neutron_damage.json
cd ../../..
```

The proton cases selecting all five names correctly produce zero neutron-only
`n_eq`. Use the neutron files for an active `n_eq` score; neutron nodes cannot
be combined with TID because they are outside the qualified Stage-4 domain.
The executable regression for all seven configurations is:

```bash
cd test_output/build
ctest --output-on-failure -R '^configuration_examples$'
cd ../..
```

Validate all of them without transport:

```bash
for case_file in examples/stage10_shapes/*.json; do
  test_output/build/shieldSim --validate-case="$case_file" || exit 1
done
```

Run one mesh case and keep every product under `test_output/`:

```bash
mkdir -p test_output/stage10_shapes/mesh_tetrahedron
cd test_output/stage10_shapes/mesh_tetrahedron
../../build/shieldSim \
  --case=../../../examples/stage10_shapes/mesh_tetrahedron.json
cd ../../..
```

Each gallery case selects its shielding at `geometry.shield_stack.layers`, its
sensitive absorber material at
`geometry.absorbers[0].geometry.regions[*].material`, and its calculations in
the top-level `quantities` array. The committed default requests `TID`, `DDD`,
`n_eq`, `LET`, and `H100/10` for natural silicon behind Al/HDPE. It uses only
20 histories so it is a quick executable example, not a converged science
result. Increase `accuracy.sampling` and supply a defensible uncertainty budget
for a real verification study; the provider purpose guard still remains in
force.

Shield layers accept the screenshot material aliases `Al`, `Cu`, `W`, `Ta`,
`HDPE`, `BPE`, `Kevlar`, `Kapton`, `Water`, and `Si`; template composites and
regolith require an explicit project material override. Sensitive regions may
select `Skin`, `EyeLens`, `BFO`, `CNS`, `SoftTissue`, `Si`, `SiO2`, `SiC`,
`GaAs`, `InGaAs`, `Ge`, or `Water`, subject to the quantity-domain gates below.
The material records and user changes remain in input files, not source code.

The current validated calculation/material combinations are intentionally
narrow: TID and LET accept Al, natural Si, or water; DDD and `n_eq` require
natural Si; `H100/10` is a primary-proton shield-exit quantity and does not
depend on absorber material. All five can therefore be selected together only
with natural silicon. `DDD`/`n_eq` additionally require `damage_response`, and
LET requires `let_response`; the supplied files include both. The parser
rejects an unused response object, so remove `damage_response` when neither
damage quantity is selected and remove `let_response` when LET is not selected.
This explicit failure is preferable to silently applying a silicon response to
a tissue or compound semiconductor.

Additional downloaded reference data are staged for future expansion, not yet
accepted by the input validator. Ge and SiO2 have exact NIST proton/He-4
stopping tables; GaAs has candidate proton/He-4 DDD and a distinct ASTM E722
neutron-equivalence package; Ge and Ga0.47In0.53As have candidate proton/He-4
DDD only. Authenticate these files with:

```bash
cd test_output/build
ctest --output-on-failure -R '^material_reference_data$'
cd ../..
```

See `tests/reference_material_functions/README.md` before attempting to widen
a case domain. It documents the missing neutron responses, target-specific
conventions, and why tissue/water do not receive fabricated lattice-damage
outputs.

The separate Boolean-enclosure example
`examples/stage10_analytic_csg_case.json` uses these same three input fields.
One JSON case defines one ordered shielding stack; make distinct case files
with distinct `case_id` and `output.prefix` values when comparing stacks.

### 11.1 Required absorber fields

An analytic-CSG absorber contains no top-level `material`. Its `geometry`
object requires:

```json
{
  "type": "analytic_csg",
  "provider_version": "1.0.0",
  "coordinate_frame": {
    "name": "absorber_local_cartesian",
    "handedness": "right",
    "origin": "provider_origin",
    "length_unit": "mm"
  },
  "transform": {"translation": {}, "rotation_quaternion_local_to_world": {}},
  "source": {
    "kind": "inline",
    "source_id": "project-controlled ID",
    "source_version": "project-controlled version",
    "provenance": "where the dimensions came from"
  },
  "regions": []
}
```

Use full measure objects for every translation and dimension; bare numbers are
rejected. Quaternions are active local-to-parent rotations, scalar-first
`w,x,y,z`, finite, and unit-normalized. The absorber transform maps provider
coordinates into the case world; a region transform maps solid coordinates
into provider coordinates; a Boolean `right_transform` maps its right child
into its left-child frame.

Every region requires a unique `region_id`, `role` (`passive` or `sensitive`),
material library ID/alias, transform, solid, and:

```json
"mass_properties": {
  "volume": {"value": 0.025, "unit": "cm3"},
  "centroid": {
    "x": {"value": 0.0, "unit": "mm"},
    "y": {"value": 0.0, "unit": "mm"},
    "z": {"value": 0.0, "unit": "mm"}
  },
  "relative_tolerance": 1.0e-9,
  "provenance": "exact rectangular-prism formula"
}
```

Centroid coordinates are in the region solid's local frame. Exactly one region
must be sensitive. The sensitive solid must be a primitive box, sphere, or
cylinder in provider version 1.0.0; passive regions may be Boolean. This
restriction keeps dose mass and region meaning unambiguous with the current
one-sensitive-region quantity accumulators.

Solid forms are:

- box: `node_id`, `type: "box"`, and `size_x`, `size_y`, `size_z`;
- sphere: `node_id`, `type: "sphere"`, `inner_radius` (zero allowed), and
  `outer_radius`;
- cylinder: `node_id`, `type: "cylinder"`, radii, and `length`; local z is the
  cylinder axis;
- Boolean: `node_id`, `type`, `left`, `right`, `right_transform`, and
  `relationship`. Use `disjoint` for union, `right_fully_contained` for
  subtraction, and `general` for intersection.

All node IDs are unique across the absorber. The parser allows at most 128
nodes and 32 nesting levels. Exact analytic volume/centroid checks apply to
primitives, declared-disjoint unions, and declared-contained subtractions.
General intersections rely on the imported properties and Geant4 volume check
and remain outside the completed validation domain. Region siblings must not
overlap; detected Geant4 overlaps terminate the run.

The sensitive region also requires `entry_surface`. Choose deliberately:

- `local_negative_z_face` for a box or cylinder planar entry patch. This uses
  `1/(A |u·n|)` and reproduces slab-v1 upstream-face normalization.
- `entire_exterior` for a primitive closed boundary. This uses
  `1/(S |u·n|)` and reports one-sided scalar fluence averaged over all exterior
  surface area. It is not equivalent to a projected front face for a
  directional source.

The resolved case archives the canonical provider geometry,
`geometry_content_sha256`, material hashes, sensitive stable ID, provider
version, and full absorber geometry hash. The response-cell key also records
provider/version and the geometry hash. Never reuse a response after any of
those values changes.

The legacy `TargetThickness [mm]` table field reports the sensitive primitive's
local-z extent for CSG cases. It does not report the complete passive assembly
depth. Use resolved bounds and `absorber_mass` in JSON for physical analysis.

### 11.2 Tessellated-mesh input contract

A mesh absorber uses the same provider-level coordinate frame, transform, and
source-provenance fields as the analytic provider, but selects
`"type":"tessellated_mesh"`. Each region has an explicit transform, material,
role, mass properties, and this mesh object:

```json
"mesh": {
  "length_unit": "mm",
  "validation_tolerance": {"value": 1.0e-9, "unit": "mm"},
  "vertices": [[0,0,0], [10,0,0], [0,10,0], [0,0,10]],
  "facets": [
    {"facet_id":"z0_entry", "vertex_indices":[0,2,1]},
    {"facet_id":"y0",       "vertex_indices":[0,1,3]},
    {"facet_id":"x0",       "vertex_indices":[0,3,2]},
    {"facet_id":"slanted",  "vertex_indices":[1,2,3]}
  ]
}
```

Coordinates are governed only by the required `length_unit: "mm"`; indices
are zero-based and triangle winding must be outward. Provider v1 accepts a
single connected, closed, consistently oriented convex shell per region,
limited to 10,000 vertices and 20,000 triangles. Its pre-Geant4 checks include:

- finite unique vertices and unique facet IDs;
- three distinct in-range indices and nonzero area per facet;
- every undirected edge exactly twice, once in each direction;
- one connected shell with Euler characteristic `V-E+F=2`;
- every oriented triangle acting as an outward supporting plane for all
  vertices, which enforces the v1 convexity/orientation boundary;
- positive signed volume; and
- exact triangle-derived volume, centroid, and surface area, with declared
  volume/centroid agreement inside the explicit tolerances.

The required `mass_properties` record has the same form as Section 11.1 and is
in the mesh-local frame. Geant4's constructed volume and sensitive mass are
then checked independently against the exact oriented-triangle result.

A sensitive mesh requires one of:

```json
"entry_surface": {
  "type": "tagged_facets",
  "facet_ids": ["z0_entry"]
}
```

or:

```json
"entry_surface": {"type": "entire_exterior"}
```

Tagged facets must form a nonempty coplanar patch with a common outward normal.
The normal, summed triangle area, boundary point-in-triangle test, and
`1/(A |u·n|)` projected-cosine weight define entry. `entire_exterior` uses the
exact summed surface and `1/(S |u·n|)`. The entry definition is part of the
canonical geometry and response key.

No STL, OBJ, CAD, or GDML import occurs in v1. The model never guesses external
units or materials and never closes holes, flips facets, merges vertices,
repairs self-intersections, or convexifies a shape. General concave meshes,
multiple shells, non-manifold topology, and automated repair require a later
versioned provider and independent validation. `scripts/generate_stage10_shape_examples.py`
reproducibly generates the six committed gallery inputs; rerun it after an
intentional gallery change and review the resulting JSON diff.

### 11.3 Run the complete Stage-10 gallery

To execute all gallery cases outside CTest while isolating outputs:

```bash
for case_file in examples/stage10_shapes/*.json; do
  case_name=$(basename "$case_file" .json)
  output_dir="test_output/stage10_shapes/$case_name"
  mkdir -p "$output_dir"
  (
    cd "$output_dir"
    ../../build/shieldSim --case="../../../$case_file"
  ) || exit 1
done
```

Every result JSON reports the provider's stable sensitive ID, Geant4 mass,
geometry hash, and the explicit Stage-10 verification-only status. The resolved
case also archives the canonical vertices/facets, derived exact properties,
entry tags, material hashes, and provider/content hashes.

### 11.4 Stage-10 verification tests and figures

Run only the bounded Stage-10 gate:

```bash
cd test_output/build
ctest --output-on-failure -R '^stage10_(analytic_csg|tessellated_mesh)$'
cd ../..
```

The pair took about 26 seconds on the development workstation on 2026-10-18;
runtime depends on processor, Geant4 build, and data cache. The analytic test runs
legacy/CSG equivalent boxes, sphere and cylinder mass cases, a rotated box,
the Boolean enclosure, DDD/n_eq entry-surface equivalence, response-key checks,
disjoint union/general intersection construction, and rejection cases. Evidence
and figures are written to:

```text
test_output/build/validation_evidence/stage10_analytic_csg.json
test_output/build/validation_figures/stage10_analytic_csg_comparisons.png
test_output/build/validation_figures/stage10_analytic_csg_comparisons.svg
test_output/build/stage10_analytic_csg_runs/
test_output/build/validation_evidence/stage10_tessellated_mesh.json
test_output/build/validation_figures/stage10_tessellated_mesh_comparisons.png
test_output/build/validation_figures/stage10_tessellated_mesh_comparisons.svg
test_output/build/stage10_tessellated_mesh_runs/
```

The mesh test runs every committed gallery case, checks exact box/tetrahedron/
octahedron volume, centroid, surface area, and Geant4 mass, compares the
analytic and triangulated boxes through same-seed TID, DDD, and `n_eq` paths,
verifies cache separation, and exercises malformed-input rejection. The
DDD/`n_eq` pair specifically executes the tagged-facet area, normal,
point-in-triangle, and projected-cosine entry estimator.

The Stage-10 test is serial because each bounded case is short. Launching the
same test under `mpiexec` would duplicate it; `shieldSim` has no Stage-10 MPI
sharding/merge contract. Use CTest's independent-test parallelism only when
running a larger selected suite, for example `ctest -j 4 ...`, and avoid
concurrent tests that target the same evidence filenames.

Current fatal checks cover explicit units/frame/version, finite positive
dimensions, radius ordering, unit quaternions, unique IDs, node limits,
material mapping, exactly one nonempty sensitive region, primitive sensitive
mass, declared mass/centroid consistency, conservative world/shield bounds,
and sampled Geant4 overlaps, plus the mesh topology/orientation/convexity and
entry-patch checks in Section 11.2. No geometry is automatically repaired. Open
qualification work includes a separately specified general concave/CAD
provider, voxel/GDML providers, deterministic exhaustive navigation/reference
rays, curved-surface tessellation/resolution convergence, multiple sensitive
responses, and an independent authoritative nontrivial shape and
radiation-response comparison.

## Stage-11 calculated spectrum exchange

The complete producer and consumer contract is
[`SPECTUM_INPUT_FORMAT.md`](SPECTUM_INPUT_FORMAT.md). The current implementation
reads the versioned JSON and runs the selected component through the normal
source sampler, Geant4 transport, and requested proxy calculations. Code under
the parent `srcEarth/` tree has not been modified; its native writer still has
to emit this contract.

Configure and compile the reader, its standalone schema test, and shieldSim:

```bash
cmake -S . -B test_output/build -DBUILD_TESTING=ON \
  -DGeant4_DIR=/path/to/Geant4/lib/cmake/Geant4
cmake --build test_output/build -j 8 \
  --target stage11_spectrum_exchange_tests shieldSim
```

The file is selected directly in the case input:

```json
"source": {
  "mode": "isotropic",
  "spectrum": {
    "format": "spectrum_exchange_v1",
    "path": "../calculated_spectra/point_000007_proton.json",
    "component_id": "local_proton"
  },
  "energy_limits": {
    "proton": {
      "min": {"value": 10.0, "unit": "MeV"},
      "max": {"value": 10.0, "unit": "GeV"}
    },
    "alpha": {
      "min": {"value": 40.0, "unit": "MeV"},
      "max": {"value": 40.0, "unit": "GeV"}
    }
  }
}
```

`path` is relative to the case file unless absolute. `component_id` must match
one component exactly. Both energy-limit blocks are required by the current
case schema; they are total kinetic energy after any per-nucleon conversion.
Schema v1 requires `mode: "isotropic"` and admits only explicit proton or alpha
components to continuous transport. A science-facing purpose additionally
requires `science_status: "calculated_validation_qualified"`; the committed
fixture is verification-only.

Validate the committed example without transport, then run it:

```bash
test_output/build/shieldSim \
  --validate-case=examples/stage11_spectrum_exchange_case.json

mkdir -p test_output/stage11_example
cd test_output/stage11_example
../build/shieldSim \
  --case=../../examples/stage11_spectrum_exchange_case.json
cd ../..
```

The example input spectrum is
`tests/data/stage11_spectrum_exchange.v1.json`. It contains a point-sampled
proton component and a per-nucleon, bin-averaged alpha component. The example
case selects the proton component. To select alpha, change only `component_id`
to `local_alpha_bins`; the existing alpha total-energy limits already cover
the fixture's converted `[40,80] MeV` bin.

Run both Stage-11 gates from inside the build directory:

```bash
cd test_output/build
ctest -R '^stage11_spectrum_exchange(_transport)?$' --output-on-failure
cd ../..
```

The first gate checks schema and case admission without starting Geant4. The
second performs a 256-history Geant4/proxy run and compares source
normalization with an independent analytic power-law plus linear-zero
integral. On the development build the pair takes only a few seconds. Evidence
and full run artifacts are written under:

```text
test_output/build/validation_evidence/stage11_spectrum_exchange_core.json
test_output/build/validation_evidence/stage11_spectrum_exchange_transport.json
test_output/build/stage11_spectrum_exchange_runs/
test_output/build/validation_figures/stage11_spectrum_exchange_analytic.png
test_output/build/validation_figures/stage11_spectrum_exchange_analytic.svg
```

Do not feed `J_omni_perMeV` or `J_planar_perMeV` into `value`. The exchange
requires AMPS `J_local_perMeV`; shieldSim applies the isotropic `pi` plane
factor. It also converts `m^-2` to its internal `cm^-2` source normalization
and applies the exact per-nucleon Jacobian where required. A producer-rejected
quality status, invalid node, unknown schema/version/field, missing
interpolation, ambiguous species, non-isotropic v1 component, or unsupported
transport species is fatal.

## 12. Understanding outputs

Depending on the case, shieldSim writes a subset of:

- `*_resolved_case.json`: canonical resolved inputs and build/data provenance;
- `*_result.json`: requested outputs, uncertainty/status, and audit ledgers;
- `*_response_cell.json`: per-primary primitive responses and complete
  covariance triangle;
- `*_run_summary.dat` and diagnostic tables: human-readable/plotting aids.

JSON is authoritative. CSV and plotting products are derived review artifacts.
Preserve the case/result/response hashes together with the exact input and log.
Do not compare results solely by filename.

For each science quantity verify:

1. units and normalization;
2. incident species and energy range;
3. shield order, material hashes, geometry hash, absorber ID/provider;
4. physics list, EM/ion choices, cut, step, seed, Geant4/data versions;
5. statistical and systematic uncertainty components and correlations;
6. `accuracy_status` and confidence interval;
7. unsupported/out-of-range contributions and limitation text;
8. capability `validation_status`.

For anomaly work, retain exact event-onset, event-end, interval midtime/source
coordinates, gaps, attitude/configuration changes, and anomaly timestamps. A
radiation proxy is an environmental covariate, not proof of anomaly causation.

## 13. Reference data and controlled downloads

The compact committed validation data live under `tests/reference_spectra/`,
`tests/reference_tid/`, `tests/reference_ddd/`, and `tests/reference_let/`.
The larger CRaTER package lives at
`/data/vtenishe/shieldSim_reference_data/reference_crater/`. Both locations
retain provenance and checksums; normal tests use those immutable local files
and do not download mutable upstream content.

Scripts such as `scripts/download_gcr_validation_data.py` and
`scripts/download_nist_stage4.py` are maintenance tools for a deliberate data
update, not normal run prerequisites. A data refresh requires review of source
version/licensing/semantics, new immutable filenames or versions, regenerated
checksums/derived tables, and all affected validation gates. Never overwrite a
pinned file merely to make a checksum test pass.

The Stage-8 acquisition can be independently verified or restored with:

```bash
python3 scripts/download_crater_stage8_inputs.py
python3 tests/test_crater_reference_inputs.py \
  "${SHIELDSIM_REFERENCE_DATA_DIR}/reference_crater"
```

The downloader verifies exact upstream bytes and uses `test_output/` for
scratch files. It range-extracts the two small BON2020 spectra embedded in a
383 MB GLACE archive instead of downloading unrelated response kernels. It
does not download the multi-terabyte CRaTER event archive; the selected
published LET aggregate is already pinned. See
`${SHIELDSIM_REFERENCE_DATA_DIR}/reference_crater/acquisition.v1.json` for the files that remain
unavailable or require preregistration.

BON2020 is not used in production and is not required to run Stage 2 or Stage
7. Calculated AMPS spectra now have a tested shieldSim consumer through
`shieldSim.spectrum_exchange.v1`; the native `srcEarth` producer remains to be
implemented. The exchange path changes only source data/provenance and does not
change transport or proxy equations.

## 13. Troubleshooting

**CMake cannot find Geant4:** set `Geant4_DIR` to the directory containing
`Geant4Config.cmake`, remove only the failed dedicated build directory if a
fresh configuration is necessary, and configure again.

**`libG4*.so` cannot be opened:** source the Geant4 environment used at build
time or correct its shared-library path. Confirm with `ldd
test_output/build/shieldSim`.

**A Geant4 data-set variable is missing:** use the matching data installation
and environment script. Do not substitute a different data version without
rerunning the affected response/validation tests.

**The parser reports an unknown field:** correct the field name or schema;
unknown fields are never ignored. Compare with the matching example and
`examples/README.md`.

**A material/library hash changes:** archive the resolved material expansion,
invalidate affected response matrices, and rerun tests traced to that material.

**Matrix cells are incompatible:** compare their `response_key_document`
objects. Geometry, materials, absorber, scoring, physics, cuts, build, Geant4
data, and registered response package must match exactly.

**A fold requests a missing node:** enlarge and rerun the response grid. Do not
enable extrapolation or clamp the source to conceal missing support.

**Accuracy is `not_met`:** increase histories only if Monte Carlo uncertainty
is the limiting component and numerical convergence is already established.
Do not turn an incomplete systematic/source budget into a number by assumption.

**Accuracy is `not_assessable`:** supply the documented missing covariance,
physics-reference, material, numerical, or source evidence. More histories do
not repair missing physical information.

**A long test appears quiet:** transport drivers archive per-run logs and CTest
may not stream them. Inspect the relevant directory below
`test_output/build/`, or run CTest with `-V`. Do not redirect test work into
`/tmp`.

## 14. Reproducible run checklist

Before using a result in an anomaly analysis or preliminary assessment:

- use a reviewed source revision and record whether the tree is dirty;
- record compiler/CMake/Geant4 versions and data paths;
- validate the versioned case before transport;
- keep every material and source file named by the case;
- write into a dedicated `test_output/<run_id>/` directory;
- use recorded, distinct random seeds for an ensemble;
- verify numerical convergence for the observable and geometry;
- preserve full logs, resolved case, result, response artifacts, and hashes;
- confirm the requested quantity lies inside its validated domain;
- report incomplete uncertainty or environment coverage explicitly; and
- never relabel the AMS proton component as full GCR or GOES directional
  ingestion as an isotropic transport boundary.
