# shieldSim development and operations handoff

Prepared: 2026-10-19  
Repository: `/data/vtenishe/GEO/AMPS/srcEarth/shieldSim`  
Prepared at source revision: `195aad67d9696a4cba6bd3409bb00ebcba359966`  
Latest implemented development scope: Stage 11 calculated-spectrum exchange
consumer

## 1. Purpose of this document

This file is the working handoff for continuing development, running the model,
reviewing results, and understanding what the current implementation does and
does not establish. It consolidates operational information that is otherwise
spread across:

- `README.md`: capability summary and input/output reference;
- `RUN.md`: authoritative build, run, testing, and troubleshooting guide;
- `quick_guide.md`: concise CRaTER, Stage-9, Stage-10, Stage-11, and recorded
  high-statistics workflows;
- `development_and_validatoin.md`: normative development gates and validation
  boundaries; and
- `SPECTUM_INPUT_FORMAT.md`: normative Stage-11 calculated-spectrum exchange
  contract. The two filenames retain their originally requested spellings.

If this handoff conflicts with a strict schema check in the code, the code and
the matching versioned schema are authoritative. If it conflicts with a
scientific completion criterion, `development_and_validatoin.md` is
authoritative.

The model is intended for spacecraft-anomaly analysis and preliminary mission
environment/effect assessment. It is not an instrument-design model, a device-
qualification model, an anatomical phantom, or proof that radiation caused a
spacecraft anomaly.

## 2. Immediate repository state

At the start of this handoff, the Stage-11 work had been committed at
`195aad67`. The source tree was clean before this new `handoff.md` was created.
Always verify the current state rather than assuming that remains true:

```bash
git status --short
git log -1 --oneline
```

Stage 11 was implemented under a strict requirement not to modify code outside
`shieldSim/`. The parent `srcEarth` AMPS flux and spectrum writers were audited
read-only. Their native exchange-format writer remains future work.

Generated builds, tests, plots, and simulations belong under `test_output/`.
That directory is intentionally ignored by Git. It is persistent project
scratch/output, but it is not reference truth and must not be committed as if
it were validated source data. Python bytecode and `__pycache__` are also
ignored and should not be created in source directories; use
`PYTHONDONTWRITEBYTECODE=1` for ad hoc Python workflows when appropriate.

The current `test_output/build/validation_figures/` contains the Stage-11
analytic comparison PNG and SVG. Other stage figures are regenerated when
their numerical tests and the full `validation_figures` step are rerun.

## 3. Development and validation status

“Implemented” and “formal stage completion” are deliberately different. A
regression test can establish software behavior without establishing physical
validity over an unrestricted material, species, geometry, or environment
domain.

| Stage | Current status | Principal remaining boundary |
|---|---|---|
| 0 | Complete requirements and validation-protocol freeze | No code or physical validation was claimed by this stage |
| 1 | Complete and acceptance-tested core data/correctness boundary | Does not itself validate a radiation-effect quantity |
| 2 | Complete and acceptance-tested published-source adapters | AMS proton admission is a bounded verification fixture, not a complete GCR environment; GOES transport remains disabled |
| 3 | Geometry/response engine implemented and acceptance-tested | Required independent review has not been recorded |
| 4 | Finite-slab TID implemented and acceptance-tested | Required independent review; implemented/tested domain remains proton/He-4 incidence in Al, natural Si, or water slabs |
| 5 | Natural-Si DDD and neutron-only `n_eq` implemented and reference-tested | Project-selected accelerator/device validation and independent review |
| 6 | Electronic LET and primary-proton `H100/10` implemented and reference-tested | Required independent review |
| 7 | Time/exposure framework and downloaded AMS A=1 proton-component qualification implemented | Full isotope-resolved GCR, GOES boundary mapping, complete source covariance, and prerequisite reviews |
| 8 | Public CRaTER inputs and foreground-only axial benchmark implemented | Full engineering response/mass model, lunar/spacecraft/background terms, covariance, preregistration, blind test, and review |
| 9 | Verification UI/compiler/preflight/resume framework implemented and tested | Earlier science gates and a reviewed production regression baseline |
| 10 | Analytic CSG and bounded convex tessellated-mesh providers implemented and tested for verification | Independent nontrivial shape/response evidence, curved-mesh convergence, broader import providers, and review |
| 11 | Exchange contract, shieldSim reader, exact sampler, and Geant4 consumer tests implemented | Native AMPS writer and a frozen real native-export conservation/comparison test |

Consequences for users:

- Results must retain their emitted `validation_status`; do not shorten that
  status to “validated.”
- Verification-only inputs or geometry cannot be made science-facing by
  changing a label or purpose field. The parser rejects many such attempts,
  but a passing parser is still not a substitute for scientific qualification.
- More Monte Carlo histories can reduce statistical error. They cannot repair
  missing source covariance, material evidence, geometry fidelity, systematic
  uncertainty, or independent review.

## 4. Important source layout

The most important implementation boundaries are:

```text
shieldSim.cc                         main case/transport executable
sourceTool.cc                        published-source adapter CLI
responseTool.cc                      response-matrix assembly/folding CLI
exposureTool.cc                      Stage-7 time/quadrature/allocation CLI
referenceTool.cc                     controlled Geant4 reference-table probe
craterBenchmark.cc                   isolated CRaTER validation executable

include/CaseDefinition.hh            resolved, strict case contract
src/CaseDefinition.cc                JSON admission and capability guards
include/SpectrumDataset.hh           canonical spectrum/component container
src/SpectrumDataset.cc               legacy and Stage-11 spectrum readers
src/PrimaryGeneratorAction.cc        species/energy/direction sampling
include/AbsorberGeometry.hh          geometry-provider abstraction
src/SlabAbsorberGeometry.cc          science-facing slab provider
src/AnalyticCsgAbsorberGeometry.cc   verification-only analytic provider
src/TessellatedMeshAbsorberGeometry.cc verification-only convex mesh provider
src/DetectorConstruction.cc          shield and absorber construction
src/SteppingAction.cc                transport ledgers and scoring hooks
src/RunAction.cc                     proxy calculations and result writers
src/ResponseMatrix.cc                response assembly/folding/covariance
src/Stage6Physics.cc                 LET and H100/10 supporting algorithms
src/Stage7Qualification.cc           exposure and qualification calculations

data/material_library.v1.json        versioned material definitions
data/niel/                           authenticated damage responses
data/let/                            authenticated LET response inputs
examples/                            directly runnable versioned cases
tests/reference_*                    compact checksum-pinned reference inputs
scripts/                             download, plotting, UI, and run helpers
test_output/                         ignored build/run/evidence/figure products
```

Large Stage-8 CRaTER inputs are intentionally external to Git:

```text
/data/vtenishe/shieldSim_reference_data/reference_crater/
```

The external package has its own README, acquisition record, and checksum
manifest. The code must not silently substitute a different directory or an
unverified file with the same name.

## 5. Software and environment setup

Required for the core build and run:

- C++17 compiler;
- CMake 3.16 or newer;
- Geant4 with CMake configuration, matching shared libraries, and required
  physics data packages;
- Python 3 for tests and helper tools; and
- Git for meaningful build/source provenance.

Plotting additionally uses `matplotlib` and NumPy. The dose-sweep plot with an
exponential fit uses SciPy. These Python packages are not core transport
dependencies.

On the development host:

```bash
export SHIELDSIM_ROOT=/data/vtenishe/GEO/AMPS/srcEarth/shieldSim
export Geant4_DIR=/data/vtenishe/G4/install/lib/cmake/Geant4
export SHIELDSIM_REFERENCE_DATA_DIR=/data/vtenishe/shieldSim_reference_data

source /data/vtenishe/G4/install/bin/geant4.sh

test -f "${Geant4_DIR}/Geant4Config.cmake"
test -f "${SHIELDSIM_REFERENCE_DATA_DIR}/reference_crater/manifest.sha256"
cmake --version
c++ --version
python3 --version
python3 -c 'import matplotlib; print(matplotlib.__version__)'
```

For another system, set `Geant4_DIR` to the directory containing
`Geant4Config.cmake` and source that installation's environment script. The
environment must expose the same Geant4 shared libraries and data packages at
runtime. In particular, HP neutron configurations require a valid
`G4NEUTRONHPDATA`.

## 6. Configure and compile

From `SHIELDSIM_ROOT`:

```bash
cmake -S . -B test_output/build \
  -DCMAKE_BUILD_TYPE=Release \
  -DBUILD_TESTING=ON \
  -DGeant4_DIR="${Geant4_DIR}" \
  -DSHIELDSIM_CRATER_REFERENCE_DIR="${SHIELDSIM_REFERENCE_DATA_DIR}/reference_crater"

cmake --build test_output/build --parallel 8
```

Adjust the build parallelism to the machine. This parallelizes compilation,
not one normal shieldSim transport case.

The important executables are:

- `test_output/build/shieldSim`;
- `test_output/build/shieldSimSourceTool`;
- `test_output/build/shieldSimResponseTool`;
- `test_output/build/shieldSimExposureTool`;
- `test_output/build/shieldSimReferenceTool`; and
- `test_output/build/shieldSimCraterBenchmark`.

CMake generates `test_output/build/generated/BuildProvenance.hh`. Results
record the revision, clean/dirty state, compiler, flags, platform, build type,
Geant4 version, and relevant data paths. Configure after the intended source
changes exist so the embedded tree state is meaningful. A copied executable is
not a self-contained release because it still depends on the matching Geant4
libraries/data and its build provenance.

After a C++-only edit, a build is normally sufficient:

```bash
cmake --build test_output/build --parallel 8
```

After changing CMake options or target definitions, rerun the configure step.

## 7. Versioned case setup

Use versioned JSON cases for all new work. The case parser deliberately rejects
unknown fields, missing units, ambiguous material/species identities, and
unsupported quantity combinations.

A normal case specifies:

1. schema version, unique case ID, and purpose;
2. material library, version, path, and optional immutable overrides;
3. world and source-plane dimensions;
4. an ordered shield stack and absorber geometry;
5. source mode, spectrum or response nodes, and energy limits;
6. requested quantities;
7. requested accuracy, coverage, uncertainty budget, and sampling limits;
8. physics list, EM/ion configuration, production cut, maximum step, and seed;
9. output prefix.

The core shield contract is one through five ordered layers. Each layer may be
specified by geometric thickness or areal density. Areal density is converted
using the resolved material density, and the resulting geometric thickness may
not exceed 10 cm. Layer order is physical and is part of the geometry and
response hashes. Every geometry dimension must appear in the input; a material
name never implies size or placement.

For new case-v2 work, start from the closest example:

```text
examples/stage3_response_cell_case.json
examples/stage4_tid_case.json
examples/stage5_ddd_neq_case.json
examples/stage6_let_hardness_case.json
examples/stage7_time_resolved_exposure.json
examples/stage10_analytic_csg_case.json
examples/stage11_spectrum_exchange_case.json
```

The shape gallery is under `examples/stage10_shapes/`. Seven full combinations
that vary incident species, absorber material, geometry, layers, thickness
form, and quantity selection are under `examples/configuration_cases/`.

Always validate before spending histories:

```bash
"${SHIELDSIM_ROOT}/test_output/build/shieldSim" \
  --validate-case="${SHIELDSIM_ROOT}/examples/stage6_let_hardness_case.json"
```

Then run in a unique output directory. Input file paths are resolved relative
to the case where defined; output prefixes are relative to the process working
directory:

```bash
mkdir -p "${SHIELDSIM_ROOT}/test_output/runs/stage6_demo"
cd "${SHIELDSIM_ROOT}/test_output/runs/stage6_demo"

"${SHIELDSIM_ROOT}/test_output/build/shieldSim" \
  --case="${SHIELDSIM_ROOT}/examples/stage6_let_hardness_case.json" \
  > run.log 2>&1
```

Do not mix products from different cases in one directory. Use unique case IDs,
seeds, and output prefixes for ensembles.

## 8. Core physical and numerical algorithms

### 8.1 Source normalization and angular sampling

Tabulated source species are sampled with weight proportional to
`J(E) dE`. The integrated spectrum is retained as the physical source
normalization. Beam mode has angular factor 1. For an isotropic directional
intensity crossing an inward plane, shieldSim samples the cosine-weighted
hemisphere and applies

```text
integral(cos(theta) dOmega over inward hemisphere) = pi.
```

The source plane is finite. Very oblique tracks can encounter finite world or
side boundaries; source-plane and transverse-size convergence must therefore
be checked for a serious calculation.

Legacy three-column spectra require a versioned metadata sidecar. The reader
does not silently skip malformed rows, sort energies, take absolute values, or
invent an undeclared width. Published AMS/ACE/GOES data must use their specific
adapters rather than the legacy interface.

### 8.2 Stage-11 calculated-spectrum interpolation

`shieldSim.spectrum_exchange.v1` supports one explicit time/location context
per document and multiple separately identified species components. A case
selects exactly one component.

For point samples:

- two positive endpoints define a power law, linear in `log(E)` and `log(J)`;
- a segment with either zero endpoint is linear in `E` and `J`;
- zero nodes are retained;
- no extrapolation is permitted; and
- energy-limit clipping preserves the original segment interpolation.

Each retained segment is integrated analytically. Energy is then sampled from
the exact inverse cumulative distribution of that segment, not from an
approximate center value. Bin-average components use exact producer edges and
piecewise-constant sampling. Gaps remain gaps.

For a per-nucleon coordinate `x=E/A`, admission converts to total kinetic
energy using:

```text
E = A x
dJ/dE = (dJ/dx) / A
```

This conserves the integrated intensity. Exchange intensity uses SI area
`m^-2`; shieldSim converts once to its `cm^-2` source-plane convention with
the factor `1e-4`.

Schema v1 accepts isotropic local differential directional intensity only.
AMPS `J_local_perMeV` is the corresponding current product. Do not supply
`J_omni=4*pi*J_local` or `J_planar=pi*J_local` as the exchange value because
shieldSim applies the plane angular factor itself.

### 8.3 Geometry and shielding

Each shield layer is an explicit Geant4 solid. Each absorber is run in an
independently rebuilt world so one absorber cannot accidentally shield another.
Case-v2 transforms use explicit translations and scalar-first `(w,x,y,z)` unit
quaternions mapping local coordinates to world coordinates.

The provider boundary is `AbsorberGeometry`:

- `slab.v1`: current science-facing finite slab;
- `analytic_csg.v1`: verification-only boxes, spheres, z-axis cylinders, and
  bounded Boolean trees;
- `tessellated_mesh.v1`: verification-only closed, connected, outward-oriented
  convex triangle shells.

Analytic and mesh providers currently allow exactly one sensitive region but
may contain multiple passive material regions. Scoring uses stable logical
region IDs rather than provider-specific copy numbers. Geometry provider,
content hash, material hashes, transforms, mass properties, and entry-surface
definition are included in the response identity.

The mesh provider does not import or repair STL, OBJ, CAD, or GDML. It rejects
holes, non-manifold edges, reversed winding, multiple shells, concavity,
ambiguous units, and invalid entry facets. The analytic provider does not make
sampled Geant4 overlap checking a mathematical proof for arbitrary Boolean
trees.

### 8.4 Transport and secondaries

The default physics list is `FTFP_BERT`. Supported runtime selections include
`FTFP_BERT`, `FTFP_BERT_HP`, `Shielding`, and `QGSP_BIC_HP`, subject to the
case and validation gate. Reference-list EM is the Stage-4 baseline;
`G4EmStandardPhysics_option4` is used as an explicit model-form sensitivity.
Ion and neutron-HP choices are recorded in provenance.

Geant4 generates and transports all secondaries allowed by the configured
physics and cuts. The raw case-v2 ledgers retain particle identity, ancestry
generation, root primary, creator process, layer deposition, track length,
step count, and rear-face crossings. Proton, alpha, and neutron exit spectra
have dedicated compatibility outputs. Creator/process partitions are
diagnostic and are not a unique causal attribution of deposited energy.

The kinetic-energy ledger is not a complete nuclear rest-mass/Q-value
conservation equation. This limitation is stated in the result.

### 8.5 Quantity definitions

- **TID:** ionizing deposited energy in the complete finite absorber divided by
  its Geant4 mass. Total and non-ionizing deposition remain audit ledgers. The
  bounded implemented/tested domain is material-slab dose, not organ or device
  dose.
- **DDD:** a NIEL exposure proxy using authenticated response data. It is not
  direct device degradation. The accepted current target is natural silicon
  for registered proton, He-4, or neutron responses.
- **`n_eq`:** ASTM E722 one-MeV-neutron-equivalent fluence, enabled only for
  natural silicon and the registered neutron contribution.
- **LET:** step-fluence electronic LET in the absorber, with explicit bins,
  charged-particle component accounting, and track-mean/dose-mean summaries.
  The separate CRaTER deposited-energy/chord estimator is validation-only.
- **`H100/10`:** primary-proton shield-exit fluence above 100 MeV divided by
  that above 10 MeV, using total proton kinetic energy. Alphas and all
  secondaries are excluded. Case-v2 uses covariance-aware ratio moments and a
  Fieller interval.

All five quantity names can be requested together only when their material and
source-domain gates are simultaneously satisfied; natural silicon is the
current common absorber. A proton-only source can legitimately produce zero
for neutron-only `n_eq`. A neutron source activates `n_eq` but is outside the
qualified TID incidence domain. DDD/`n_eq` require `damage_response`; LET
requires `let_response`. Unused response objects are rejected.

### 8.6 Response matrices and time integration

Monoenergetic response cells store per-primary sums, squared sums, means,
sample variances, standard errors, and the complete primitive-score covariance
upper triangle. Response keys bind geometry, materials, provider, scoring,
physics, cuts, build/code, Geant4, and data identities. They deliberately do
not bind seed or event count.

`shieldSimResponseTool` assembles compatible cells and folds explicit source
weights. It refuses missing/duplicate nodes, response-key mismatch,
extrapolation, invalid covariance, and unit mismatch. Nonlinear ratios fold
their primitive numerator/denominator responses rather than a precomputed
ratio.

`shieldSimExposureTool` performs strict UTC interval integration,
spectrum-response quadrature, and pilot-based event allocation. GOES intervals
retain start, midpoint, end, direction, quality, rejected/missing ledgers, and
exact overlap with the approved half-open event window. Missing or rejected
samples are not interpolated.

### 8.7 Accuracy and uncertainty

The default requested relative combined standard uncertainty is 0.10 at `k=1`
with nominal 68% coverage. A two-sided 95% interval is also reported.
Per-quantity overrides are allowed. Adaptive transport evaluates cumulative
event statistics only between bounded batches.

The result reports `met`, `not_met`, or `not_assessable`. A low Monte Carlo
standard error does not imply a complete physical uncertainty budget.
Systematic floors, missing covariance, unsupported response contributions, and
incomplete source/environment coverage remain explicit.

## 9. Stage-11 input and run procedure

The complete format is in `SPECTUM_INPUT_FORMAT.md`. The case binding is:

```json
"source": {
  "mode": "isotropic",
  "spectrum": {
    "format": "spectrum_exchange_v1",
    "path": "../calculated_spectra/run42_point7.json",
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

`path` is resolved relative to the case unless absolute. `component_id` must
match exactly. Both energy-limit blocks are currently required and refer to
total kinetic energy after per-nucleon conversion.

Transport admission currently supports explicit proton `(Z=1,A=1)` or alpha
`(Z=2,A=4)` components. Other valid species can be preserved by the exchange
reader but are not yet admitted to the continuous source generator. Schema v1
is isotropic only. A science-facing case requires
`science_status: calculated_validation_qualified`; the committed fixture is
`verification_fixture`. Merely changing this string is not physical
qualification.

Run the committed example:

```bash
"${SHIELDSIM_ROOT}/test_output/build/shieldSim" \
  --validate-case="${SHIELDSIM_ROOT}/examples/stage11_spectrum_exchange_case.json"

mkdir -p "${SHIELDSIM_ROOT}/test_output/runs/stage11_example"
cd "${SHIELDSIM_ROOT}/test_output/runs/stage11_example"

"${SHIELDSIM_ROOT}/test_output/build/shieldSim" \
  --case="${SHIELDSIM_ROOT}/examples/stage11_spectrum_exchange_case.json" \
  > run.log 2>&1
```

The fixture contains proton point samples and an alpha per-nucleon bin. The
example selects `local_proton_points`; select `local_alpha_bins` to exercise
the alpha conversion.

## 10. Tests and validation evidence

CTest is the authoritative staged entry point. This host uses CMake/CTest 3.16;
change into the build directory because its CTest does not reliably support the
newer `ctest --test-dir` syntax:

```bash
cd "${SHIELDSIM_ROOT}/test_output/build"
ctest -N
ctest --output-on-failure
```

The full suite contains transport-heavy Stage-7 and Stage-8 work and may be
long. Use a targeted expression while developing, then run the applicable
upstream gates before claiming completion:

```bash
ctest --output-on-failure -R '^stage1_'
ctest --output-on-failure -R '^stage2_'
ctest --output-on-failure -R '^stage4_(core|nist_reference|tid_transport)$'
ctest --output-on-failure -R '^stage6_'
ctest --output-on-failure -R '^stage10_(analytic_csg|tessellated_mesh)$'
ctest --output-on-failure -R '^stage11_spectrum_exchange(_transport)?$'
```

Registered test groups are:

```text
stage1_core, stage1_case_boundaries, stage1_transport_smoke
stage2_source_adapters, stage2_reference_spectra, stage2_source_tool,
stage2_published_transport
stage3_core, stage3_response_weights, stage3_transport
stage4_core, stage4_nist_reference, stage4_tid_transport
stage5_core, stage5_transport
stage6_core, stage6_case_boundaries, stage6_reference,
material_reference_data, stage6_transport
stage7_core, stage7_boundaries, stage7_published_qualification
stage8_reference_inputs, stage8_crater_core, stage8_crater_transport
stage9_interface
stage10_analytic_csg, stage10_tessellated_mesh, configuration_examples
stage11_spectrum_exchange, stage11_spectrum_exchange_transport
validation_figures
```

Latest Stage-11 focused verification at handoff:

```text
stage11_spectrum_exchange              PASS
stage11_spectrum_exchange_transport    PASS
2/2 tests passed in approximately 3 seconds
```

Stage-1 and Stage-2 source regressions were also exercised during Stage-11
implementation and passed. The complete all-stage suite was not rerun as part
of the final Stage-11 check; do not report an all-stage pass without running it
again.

Numerical JSON evidence is authoritative. Figures are review aids. Normal CTest
artifacts are under:

```text
test_output/build/validation_evidence/
test_output/build/validation_figures/
test_output/build/stage*_runs/
```

## 11. Plotting and comparison figures

### 11.1 Unified validation figures

`tests/generate_validation_figures.py` reads already-passed numeric evidence
and emits both PNG and SVG plus a manifest. It does not introduce pixel-based
acceptance criteria. In a complete suite, `validation_figures` is ordered after
its numerical dependencies.

Run the complete applicable suite from the build directory to populate all
evidence, then run or inspect `validation_figures`. Selecting only that one test
does not manufacture missing prerequisite evidence.

The current Stage-11 comparison is:

```text
test_output/build/validation_figures/stage11_spectrum_exchange_analytic.png
test_output/build/validation_figures/stage11_spectrum_exchange_analytic.svg
```

It compares the source normalization used by the full Geant4/proxy run with an
independent analytic integral of a positive power-law segment plus a
linear-to-zero segment, and separately checks the isotropic `pi` factor and
25 cm2 source-plane rate. It validates the consumer/sampler, not AMPS physics
or a native AMPS export.

Stage-10 figure names are:

```text
stage10_analytic_csg_comparisons.{png,svg}
stage10_tessellated_mesh_comparisons.{png,svg}
```

### 11.2 Transport spectra

For a generated `*_spectra.dat` file:

```bash
python3 "${SHIELDSIM_ROOT}/scripts/plot_shield_spectra.py" \
  path/to/run__absorber_spectra.dat \
  --save --format png --dpi 240 --no-show

python3 "${SHIELDSIM_ROOT}/scripts/plot_shield_spectra.py" \
  path/to/run__absorber_spectra.dat \
  --save --format svg --no-show
```

The plot is written beside the input with the requested extension. The script
plots raw per-primary and source-normalized proton/alpha transport spectra. A
gap or separated line segment normally means zero/empty bins; do not draw
through those bins or replace them with a small positive number merely to make
the plot continuous.

### 11.3 LET spectrum

```bash
PYTHONDONTWRITEBYTECODE=1 python3 \
  "${SHIELDSIM_ROOT}/scripts/plot_let_spectrum.py" \
  path/to/run__absorber_let_spectrum.dat \
  --output-prefix path/to/run__absorber_let_spectrum \
  --formats pdf png svg --dpi 240
```

The plotter uses the authoritative `total` zone, displays Monte Carlo
one-standard-error bands, and exposes per-bin relative precision. It omits
zero-content bins on logarithmic axes rather than altering them. For the
recorded eight-process example, `--summary-json` adds pooled track- and
dose-mean reference lines; see the exact procedure in `quick_guide.md`.

### 11.4 Dose sweep

```bash
python3 "${SHIELDSIM_ROOT}/scripts/plot_dose_sweep.py" \
  path/to/shieldSim_dose_sweep.dat \
  --xaxis gcm2 --yunit mGy --fit --halving \
  --save --format png --dpi 240 --no-show
```

`--fit` requires SciPy. The exponential fit and halving thickness are derived
review products, not additional validated physics. A single-panel alternative
is `scripts/plot_dose_sweep_single_panel.py`.

### 11.5 CRaTER diagnostic comparison

```bash
python3 "${SHIELDSIM_ROOT}/scripts/plot_crater_foreground_benchmark.py" \
  path/to/crater_result.json \
  path/to/crater_comparison
```

This writes PNG and SVG. It is intentionally labeled foreground-only axial and
must not be presented as a full CRaTER validation.

### 11.6 Plot-review rules

- Numeric JSON assertions and evidence are authoritative, not visual agreement.
- Preserve zeros, invalid masks, gaps, and under-resolved bins.
- Always include units, normalization, source identity, and validation scope.
- Never compare files solely by filename; compare hashes and resolved case
  identities.
- Archive vector SVG or PDF alongside a high-resolution PNG for publication
  preparation, while retaining the numeric evidence used to draw it.

## 12. Parallel execution

### 12.1 Main shieldSim executable

One normal `shieldSim` case is serial. It does not have an MPI rank-partition
or canonical merge contract. Do not run:

```text
mpiexec -n N shieldSim --case=...
```

Every rank would repeat the full histories and could overwrite the same files.

Independent cases can run concurrently in separate directories with different
output prefixes and seeds. CTest can run independent tests with `-j`, but do
not allow concurrently selected tests to write the same evidence filenames.

The recorded 100,000-history, eight-process BON2020 proton example in
`quick_guide.md` creates eight independent case files, fixes 12,500 histories
per worker, uses distinct seeds and directories, checks source/geometry/physics
identity, and pools only review spectra and summary statistics. It explicitly
does not create a canonical merged result JSON. Reuse that procedure only with
the same integrity checks; simple averaging of ratios or standard errors is
wrong.

### 12.2 Stage-7 campaign workers

Stage-7 parallelism is case-level. One driver schedules independent serial
shieldSim processes:

```bash
cd "${SHIELDSIM_ROOT}/test_output/build"
SHIELDSIM_STAGE7_WORKERS=8 \
  ctest --output-on-failure -R '^stage7_published_qualification$'
```

The worker count is a scheduling choice, not a physical input. Avoid combining
a large worker count with aggressive `ctest -j`, which can oversubscribe CPUs
and memory.

### 12.3 CRaTER supported sharding

The separate `shieldSimCraterBenchmark` has a validated shard/merge boundary.
On one workstation:

```bash
python3 "${SHIELDSIM_ROOT}/scripts/run_crater_parallel.py" \
  --benchmark "${SHIELDSIM_ROOT}/test_output/build/shieldSimCraterBenchmark" \
  --case "${SHIELDSIM_ROOT}/examples/stage8_crater_foreground_axial.json" \
  --workers 4
```

For an MPI launcher, the executable reads Open MPI or PMI rank/size variables
without linking MPI:

```bash
cd "${SHIELDSIM_ROOT}"
mpiexec -n 4 test_output/build/shieldSimCraterBenchmark \
  examples/stage8_crater_foreground_axial.json

python3 scripts/merge_crater_shards.py \
  --output test_output/crater_foreground_axial/crater_merged_result.json \
  test_output/crater_foreground_axial/*_partial_result.json
```

The merge tool rejects missing/duplicate ranks, different case/source/
observation/geometry/physics identities, inconsistent bins/history totals, and
non-distinct seeds. It combines raw importance-weight sums and sums of squares
before recomputing global values and errors.

## 13. Output and run storage

Use a unique directory under `test_output/` for every run or ensemble:

```text
test_output/runs/<campaign>/<case_id>/<revision_or_timestamp>/
```

Recommended retained package:

- original case JSON and all referenced source/material/response files;
- `run.log` and exact command/environment record;
- canonical `*_resolved_case.json`;
- authoritative `*_result.json`;
- `*_response_cell.json` when generated;
- spectra, LET, quantities, exit histograms, and run summary;
- numeric validation/comparison evidence;
- generated figures and figure manifest;
- source revision and dirty-state record.

Important products:

- `*_resolved_case.json`: resolved inputs, materials, geometry/source hashes,
  build/Geant4/data provenance, seed, and capability state;
- `*_result.json`: quantity values, units, normalization, ledgers,
  uncertainty/intervals, `accuracy_status`, and `validation_status`;
- `*_response_cell.json`: primitive per-primary moments and full covariance;
- `*_spectra.dat`: raw per-primary and source-normalized transport spectra;
- `*_quantities.dat`: compatibility scalar table;
- `*_let_spectrum.dat`: total and component LET spectra;
- `*_dose_sweep.dat`: legacy sweep product; and
- `*ProtonExitE.csv`, `*AlphaExitE.csv`, `*NeutronExitE.csv`: Geant4 exit
  histograms.

JSON is authoritative. Text/CSV files and plots are derived review products.
An interrupted adaptive run rewrites the result after completed batches; the
final file represents all completed cumulative histories.

Do not store normal test or simulation output in `/tmp`. Do not commit
`test_output/` as reference data. Reference inputs belong in a checksum-
controlled reference directory with provenance and licensing information.

## 14. Reference data governance

Compact committed references are under:

```text
tests/reference_spectra/
tests/reference_tid/
tests/reference_ddd/
tests/reference_let/
tests/reference_material_functions/
```

The larger CRaTER package is external. Normal builds/tests are offline and use
the pinned local bytes. Download scripts are maintenance tools, not automatic
runtime dependencies. A data refresh requires:

1. review of source, version, licensing, units, and scientific semantics;
2. immutable new filenames or version identifiers;
3. regenerated checksum manifests and transformation records;
4. rerun of all affected tests and figures; and
5. explicit review of changed validation conclusions.

Never overwrite a pinned file merely to satisfy a checksum. BON2020 is not a
production runtime model and is not required for Stage 2 or Stage 7. It is used
in bounded CRaTER/reference examples with explicit provenance.

## 15. Known issues and hard limitations

### 15.1 Scientific qualification

- Stages 3, 4, and 6 still require independent review.
- Stage 5 lacks a project-selected accelerator/device validation data set and
  independent review.
- Stage 7 is not a complete GCR/SEP physical qualification.
- The CRaTER result is a foreground-only axial diagnostic, not a full
  observational validation.
- The Stage-9 UI is verification-only and its candidate regression snapshot is
  unreviewed.
- Stage-10 non-slab providers are verification-only.
- Stage 11 lacks a native AMPS producer and native-export validation.

### 15.2 Source and environment limitations

- The AMS A=1 proton fixture is not a full isotope-resolved GCR environment.
- AMS inclusive helium and heavier elemental data remain ingestion-only until
  isotope composition is defined for transport.
- GOES directions and quality/gap semantics are preserved, but transport is
  disabled until a sensor-to-boundary angular mapping is approved and tested.
- Stage-11 v1 is isotropic and transports one selected proton or alpha
  component per case.
- Current AMPS Tecplot spectrum products are not a robust direct interface.
  Flux files contain integrated totals and cannot reconstruct a unique
  differential spectrum.

### 15.3 Geometry limitations

- Science-facing qualified geometry remains the finite slab.
- Analytic CSG has bounded primitive/Boolean semantics and sampled overlap
  checks, not exhaustive navigation proof.
- Mesh v1 is convex, closed, connected, triangle-only, and inline. It does not
  support general concave CAD/STL repair/import.
- Only one sensitive region is supported per independently run absorber.
- Voxel, GDML, project-CAD adapters, curved-mesh refinement, and
  multiple-sensitive-region response contracts remain future work.

### 15.4 Quantity/material limitations

- TID does not validate photons, electrons, heavy ions, devices, organs, or
  crew dose outside its registered bounded domain.
- DDD is a NIEL proxy, not observed component degradation.
- `n_eq` is currently natural-silicon, neutron-only.
- Candidate GaAs, Ge, SiO2, and InGaAs reference data do not automatically
  widen the admitted calculation domain.
- `H100/10` excludes secondaries and non-proton primaries by definition.
- Unsupported response contributions result in `not_assessable`; they are not
  clamped or silently omitted.

### 15.5 Numerical and operational limitations

- Main shieldSim has no same-case MPI merge.
- Finite source/world geometry can affect oblique isotropic transport.
- Results are sensitive to cuts, maximum step, response-grid coverage,
  geometry size, physics list, and nuclear-data configuration; convergence
  evidence must accompany serious use.
- Sparse secondary spectra can look like disconnected plot segments because
  empty bins are correctly retained as gaps.
- CTest 3.16 should be run from the build directory; the newer `--test-dir`
  option is not reliable on this host.
- Full validation figures require the complete prerequisite evidence set. A
  plot-only command cannot replace missing numeric tests.

## 16. Stage-11 producer work still required

The parent AMPS audit found:

- `boundary/spectrum.h/.cpp` defines per-particle/per-nucleon energy, SI
  differential intensity, time policies, and boundary uncertainty;
- `gridless/DensityGridless.cpp` and `3d/DensityMode3D.cpp` write local energy
  nodes, conservative lower/upper curves, sampled/resolved counts, unresolved
  fractions, and derived `4*pi*J`/`pi*J` columns;
- Mode3D additionally carries allowed counts and field/product fingerprints;
  and
- current flux files contain integrated total/channel/planar values but are
  not invertible to a spectrum.

The next Stage-11 step requires authorization to modify parent `srcEarth` code:

1. add an atomic writer in gridless and/or Mode3D that emits
   `shieldSim.spectrum_exchange.v1` without removing the existing diagnostics;
2. emit one unambiguous location/epoch document, or an explicit manifest of
   such documents;
3. use `J_local` and preserve exact energy basis, interpolation, valid mask,
   bounds, counts, frame, UTC, source/snapshot/input fingerprints, and quality;
4. independently reconstruct a frozen native export and compare every field;
5. integrate the exchange curve independently and require agreement with the
   AMPS integrated flux products inside preregistered numeric tolerances;
6. run that frozen export through shieldSim and archive source, response,
   quantity, build, and data identities.

Do not parse current Tecplot files directly inside shieldSim, infer metadata
from filenames/zone order, treat unresolved bounds as Gaussian uncertainties,
bridge invalid/zero nodes, or relabel `J_omni`/`J_planar` as local intensity.

## 17. Recommended next actions

In priority order:

1. Commit this handoff after review.
2. Rerun the complete applicable CTest suite on a clean intended revision and
   archive the evidence/figure manifest.
3. Obtain permission and define the target gridless/Mode3D producer scope for
   the native Stage-11 writer.
4. Freeze one small real AMPS output with complete run and snapshot provenance
   for the Stage-11 producer/consumer comparison.
5. Arrange the outstanding independent reviews for Stages 3, 4, and 6.
6. Select a Stage-5 accelerator/device validation data set before widening
   damage claims.
7. Keep Stage-8 and Stage-10 claims bounded until authoritative geometry and
   independent response evidence exist.

## 18. Reproducible-analysis checklist

Before handing a result to an anomaly or mission study, confirm:

- [ ] source revision and dirty state recorded;
- [ ] compiler, CMake, Geant4, physics data, and reference package recorded;
- [ ] versioned case validated before transport;
- [ ] source, materials, responses, geometry, and hashes archived;
- [ ] unique run directory, output prefix, and random seed used;
- [ ] requested species/energy/material/geometry lie inside the admitted
      domain;
- [ ] numerical convergence checked for source plane, world, cut, step,
      histories, and response grid as applicable;
- [ ] statistical and systematic uncertainties kept separate;
- [ ] `accuracy_status`, `validation_status`, limitations, and rejected or
      out-of-domain contributions reviewed;
- [ ] JSON evidence retained with plots;
- [ ] no source component relabeled as a complete environment;
- [ ] radiation result described as an environmental proxy, not causal proof.
