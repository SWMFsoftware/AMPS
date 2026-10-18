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
runs because validation inputs are checksum-pinned in the repository.

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
  -DGeant4_DIR="${Geant4_DIR}"
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

## 9. Understanding outputs

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

## 10. Reference data and controlled downloads

Committed validation data live under `tests/reference_spectra/`,
`tests/reference_tid/`, `tests/reference_ddd/`, and `tests/reference_let/` with
provenance and checksums. Normal tests use those immutable files and do not
download mutable upstream content.

Scripts such as `scripts/download_gcr_validation_data.py` and
`scripts/download_nist_stage4.py` are maintenance tools for a deliberate data
update, not normal run prerequisites. A data refresh requires review of source
version/licensing/semantics, new immutable filenames or versions, regenerated
checksums/derived tables, and all affected validation gates. Never overwrite a
pinned file merely to make a checksum test pass.

BON2020 is not used in production and is not required to run Stage 2 or Stage
7. AMPS is also not currently a runtime source. A later AMPS adapter must map a
versioned AMPS output into the existing canonical source contract without
changing transport or proxy equations.

## 11. Troubleshooting

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

## 12. Reproducible run checklist

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
