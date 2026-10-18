# shieldSim

`shieldSim` is a Geant4 environment-effects model for spacecraft-anomaly
analysis and preliminary mission assessment. It transports energetic particles
through an ordered planar shielding stack into a finite absorber slab. It is
not an instrument-design model, a device-response model, or an anatomical
phantom.

## Implementation status

Stages 1 and 2 are complete. Stages 3 and 4 are implemented and acceptance-
tested but are **pending the required independent review**, so neither is
formally complete. Stage 5 is implemented and reference-tested for natural-
silicon DDD from protons, helium-4, and neutrons and for neutron-only ASTM E722
`n_eq`; formal Stage-5 completion is pending a project-selected accelerator or
device-degradation validation data set and independent review. Stage 6 is
implemented and reference-tested for generic electronic LET and primary-proton
`H100/10`; formal Stage-6 completion is pending independent review. Stage 7
implements strict time-resolved exposure integration, response/source
quadrature, pilot-based event allocation, and the approved downloaded AMS A=1
proton component representative-matrix qualification. It remains explicitly
incomplete as a physical GCR/SEP qualification because a full isotope-resolved
GCR scenario, the GOES sensor-to-simulation angular mapping, complete physical
uncertainty budgets, and required prior-stage reviews are unavailable. Stage 4 enables a bounded finite-material-slab TID
capability for proton and helium-4 incidence in aluminum, silicon, and water.
It does not validate photon TID, electrons, heavier ions, devices, organs, or
crew dose.
Stage 2 has strict
adapters and pinned data for the original AMS-02 proton/helium fixture, the
2025 AMS-02 He--O solar-cycle tables, ACE/CRIS B--Ni, and the GOES-16 GLE72
time series/event integration. Its downloaded AMS-02 proton fixture also
drives an end-to-end validation case through the primary generator, Geant4,
TID, and `H100/10`. BON2020 is not used in production and is not a Stage-2
dependency. The science-facing case path is a strict
`shieldSim.case.v1` JSON file and provides:

- one through five ordered shield layers, each limited to an inclusive 100 mm
  derived geometric thickness;
- either `g/cm2` areal density or geometric length per layer, with the other
  value derived from the resolved material density;
- an external, versioned material library with immutable case overrides,
  status gates, property provenance, canonical expansion, and SHA-256 hashes;
- explicit world, source-plane, shield, placement, and absorber dimensions—no
  hidden geometry defaults in a versioned case;
- one independent Geant4 world per shield × absorber calculation, so an
  absorber never shields another absorber;
- a geometry-neutral absorber-provider contract and a current `slab.v1`
  implementation;
- strict spectrum metadata, input-range rejection, and unknown-field
  rejection;
- separate total, ionizing, and non-ionizing deposition ledgers;
- the approved primary-proton H100/10 definition;
- a default 0.10 combined standard relative-uncertainty request, quantity
  overrides, adaptive batches, 95% intervals, and
  `met`/`not_met`/`not_assessable` states; and
- resolved case, geometry, material, source, Geant4, seed, and output
provenance.

The Stage-3/4/5/6 `shieldSim.case.v2` path adds explicit
local-to-world translation/unit-quaternion transforms for the source plane,
shield reference surface, and absorber reference surface; monoenergetic
proton, He-4, and explicit-ion response nodes; exact discrete-node direct
sampling; ancestry-aware raw transport ledgers; and deterministic response
matrix assembly/folding. It permits `transport_raw`, the Stage-4 `TID`, the
separately gated Stage-5 `DDD` and `n_eq`, and the Stage-6 `LET` and `H100/10`
quantities. TID is ionizing energy deposited in the complete finite slab divided
by that slab's Geant4 mass. It is labeled material-slab dose, never organ or
device dose. Crew, device, and instrument claims remain disabled on case v2.
Gamma nodes are accepted only for controlled attenuation
verification. The case-v2 parser rejects a `TID` request containing a gamma or
explicit heavy-ion response node, so those species cannot enter the validated
TID domain accidentally. It likewise rejects a TID absorber outside the
validated aluminum, silicon, and water material identities; this does not
restrict the independently declared shielding-stack materials.

Stage 4 adds an explicit EM-physics choice (`reference_list_default` for the
production baseline or `G4EmStandardPhysics_option4` for sensitivity),
particle/ancestry deposition partitions, separately labeled process-name
diagnostics, exact slab-mass and Gy/rad checks, and a structured uncertainty
budget. A complete case-v2 TID budget must name material/mass, source,
numerical, and physics-reference standard-uncertainty components, declare every
pairwise correlation, pass a positive-semidefinite covariance check, and select
the combined 95% coverage method. A missing budget remains runnable but yields
`not_assessable` rather than an unsupported accuracy claim.

Stage 5 adds a checksum-authenticated damage-response package. Its supported
matrix is intentionally narrow: natural silicon (`Si_natural_G4_11_4`) with
proton, helium-4, and neutron DDD, plus neutron-only ASTM E722-19 silicon
`n_eq`. Charged particles use SR-NIEL-7 v11.1 positive nodes with log-log
interpolation; neutrons use explicit ASTM groups through 20 MeV. Exact
absorber-entry crossing energies are authoritative. A nonzero unsupported or
out-of-domain contribution makes the affected quantity `not_assessable`; the
code never clamps, extrapolates, or silently drops it. Geant4 nonionizing
deposition per slab mass is retained as a separate diagnostic and is never
added to the table fold.

Stage 6 adds step-fluence electronic LET and covariance-aware spectral
hardness. LET uses unrestricted Geant4 electronic stopping at the midpoint of
each positive-length charged absorber step, weighted by track length divided
by absorber volume. It reports all charged species—including delta electrons—
as distinct charge/ancestry/generation components, plus track- and dose-mean
ratios in keV/um and MeV cm2/g. `H100/10` weights every forward primary-proton
rear-face crossing by `1/(A|n dot u|)` and reports a ratio of event means with
event covariance and a 95% Fieller interval. Zero or statistically unresolved
denominators are undefined. The validated reference matrix is p/He-4 in
natural Al, natural Si, and water; C-12 is limited to natural Si and requires
`G4EmStandardPhysics_option4`. Values outside declared LET bins are audited and
make the result not assessable; they are never clamped.

Stage 7 adds the Geant4-free `shieldSimExposureTool` for strict UTC interval
integration, explicit gaps and configuration segments, anomaly markers,
covariance and sampled-ensemble propagation, adaptive spectrum/response
quadrature, and deterministic pilot-based Neyman history allocation without
physics biasing. Its registered transport campaign uses five seeds, both
Al/HDPE layer orders, independent natural-Si and water slabs, full downloaded
AMS proton bins, and direct/quadrature/response-matrix paths. Log-energy
response interpolation cannot extrapolate; the final registered grid monitor
must change by no more than 2%, and direct/matrix TID must satisfy both the 10%
relative and combined two-sided 95% gates. Missing published source covariance
is reported as `not_assessable`, never replaced by an assumed diagonal model.
The two GOES directions remain separate and transport-disabled until their
boundary transform and angular response are supplied.

Stage 1 is a software/correctness gate, not physical validation of every
radiation proxy. Versioned case v1 enables its earlier `TID` and `H100/10`
path; case v2 adds the bounded Stage-4 TID, Stage-5 silicon damage, and Stage-6
LET/hardness capabilities described above. Every result retains its validation
and uncertainty status. The older
CLI can still produce its historical analytic DDD/LET/n_eq estimates for
regression compatibility; those results are prototype diagnostics, not
science-facing outputs.

The implementation and evidence needed to turn the current prototype into a
physically qualified radiation-proxy model are defined in
[`development_and_validatoin.md`](development_and_validatoin.md). Published
AMS-02 and ACE/CRIS GCR observations, GOES-16 SEP fixtures, and a CRaTER
observational LET product are documented under `tests/reference_spectra/`.
Dedicated Stage-2 adapters preserve their published binning,
species, counts/uncertainties, timing, and quality semantics. They remain
ingestion/event-integration products except for one narrow transport-test
admission: the explicit A=1 proton component of the original AMS fixture may
be used only by `verification_test` cases in isotropic mode. It is not a
complete mission GCR environment or a production default. Published data must
not be passed through the legacy three-column spectrum interface.

The package retains the old CLI for compatibility while adding versioned core
components:

- CLI-selectable Geant4 physics list (`FTFP_BERT`, `FTFP_BERT_HP`, `Shielding`, `QGSP_BIC_HP`),
- a strict case/material/spectrum schema with explicit units and provenance,
- a separate legacy catalog/CLI for older regression cases,
- source normalization using `J(E) dE`,
- `beam` and `isotropic` source modes controlled by CLI,
- local-coordinate shield rear-face scoring,
- separate total, ionizing, and non-ionizing deposit ledgers,
- versioned Stage-4 TID, Stage-5 silicon DDD/neutron-only `n_eq`, and Stage-6 LET/H100/10, plus explicitly labeled legacy diagnostics,
- detailed help text describing input/output units and radiation-effect approximations,
- Layer-1/2/3/4 scripts plus Stage-5/6/7 CTest gates with deterministic diagnostic and regression output,
- optional test/numerical controls: `--dump-run-summary`, `--production-cut`, and `--max-step`.

## Directory structure

```text
shieldSim/
  ├── CMakeLists.txt
  ├── README.md
  ├── RUN.md                         # build/run/test operations guide
  ├── shieldSim.cc
  ├── sourceTool.cc                  # published-source validation/conversion CLI
  ├── responseTool.cc                # response-matrix assembly/folding CLI
  ├── exposureTool.cc                # Stage-7 time/quadrature/allocation CLI
  ├── referenceTool.cc               # Stage-4 Geant4 EM table probe
├── include/
│   ├── AbsorberGeometry.hh       # future-shape provider contract + slab
│   ├── AccuracyController.hh     # running uncertainty/stopping state
│   ├── CaseDefinition.hh         # strict resolved case data
│   ├── DamageResponseLibrary.hh  # authenticated Stage-5 NIEL/E722 responses
│   ├── Stage6Physics.hh           # LET bins, analytic hardness, Fieller ratios
│   ├── Stage7Qualification.hh     # time, quadrature, allocation contracts
│   ├── GeometryTransform.hh       # canonical local-to-world transforms
│   ├── MaterialLibrary.hh        # versioned records and overrides
│   ├── ResponseMatrix.hh         # response cell/matrix/fold contracts
│   ├── SpectrumDataset.hh        # canonical spectrum contract
│   ├── StrictJson.hh             # dependency-free strict JSON
│   └── ...                       # Geant4 actions and legacy modules
├── src/
│   ├── CaseDefinition.cc
│   ├── DamageResponseLibrary.cc
│   ├── MaterialLibrary.cc
│   ├── MaterialFactory.cc
│   ├── SpectrumDataset.cc
│   ├── StrictJson.cc
│   └── ...
├── data/
│   ├── material_library.v1.json
│   └── niel/                    # Stage-5 manifest, tables, and provenance
├── examples/
│   ├── README.md
│   ├── stage1_case.json
│   ├── stage3_response_cell_case.json
│   ├── stage4_tid_case.json
│   ├── stage5_ddd_neq_case.json
│   ├── stage6_let_hardness_case.json
│   ├── stage7_time_resolved_exposure.json
│   ├── sep_spectrum.dat
│   └── sep_spectrum.metadata.v1.json
├── macros/
│   └── example_commands.sh
└── tests/
    ├── README.md
    ├── stage1_core_tests.cc
    ├── stage3_core_tests.cc
    ├── test_stage3_transport.py
    ├── stage4_core_tests.cc
    ├── test_stage4_nist_reference.py
    ├── test_stage4_tid_transport.py
    ├── stage5_core_tests.cc
    ├── test_stage5_transport.py
    ├── stage6_core_tests.cc
    ├── test_stage6_reference.py
    ├── test_stage6_transport.py
    ├── stage7_core_tests.cc
    ├── test_stage7_boundaries.py
    ├── test_stage7_published_qualification.py
    ├── reference_stage7/
    ├── reference_ddd/
    ├── reference_let/                # pinned GSI ATIMA heavy-ion responses
    ├── generate_validation_figures.py
    ├── generate_layer_validation_figures.py
    ├── test_stage1_case_boundaries.py
    ├── test_reference_spectra.py
    ├── run_layer1_tests.sh
    ├── run_layer2_tests.sh
    ├── run_layer3_tests.sh
    ├── run_layer4_tests.sh
    ├── regression_tools.py
    ├── data/
    │   ├── mono_spectrum_metadata.v1.json
    │   └── mono_*.dat
    ├── reference_spectra/
    ├── reference_tid/nist/          # raw + canonical pinned NIST tables
    ├── expected/
    │   └── layer4/
    └── logs/
```

## Automated tests

The package includes layered test scripts under `tests/`. Layer-1 covers build
and CLI smoke tests, Layer-2 covers geometry/source/scoring diagnostics, Layer-3
covers physics-output sanity checks, and Layer-4 performs baseline-driven
end-to-end regression tests. See `tests/README.md` for details.

Layer-4 baselines are created only after manual review:

```bash
tests/run_layer4_tests.sh --update-baseline
```

Subsequent checks compare against that accepted baseline:

```bash
tests/run_layer4_tests.sh --check
```

## Build

The complete operational instructions—including prerequisites, Geant4 paths,
case setup, every executable, response matrices, Stage-7 integration, tests,
outputs, and troubleshooting—are in [`RUN.md`](RUN.md).

Source your Geant4 environment first, then build with CMake:

```bash
cmake -S . -B test_output/build -DBUILD_TESTING=ON
cmake --build test_output/build -j
```

If CMake cannot find Geant4 automatically, provide `Geant4_DIR`, for example:

```bash
cmake -S . -B test_output/build \
  -DBUILD_TESTING=ON \
  -DGeant4_DIR=/path/to/geant4/lib/cmake/Geant4
```

Run all registered Stage-1 through Stage-7 tests from the build directory:

```bash
cd test_output/build
ctest --output-on-failure
cd ../..
```

`stage1_core` verifies canonical JSON/hashing, material resolution and
overrides, spectrum metadata, unit conversion, and accuracy mathematics.
`stage1_case_boundaries` exercises the production CLI parser at the one/five/
six-layer and 100 mm boundaries and checks required fields, units, ranges,
accuracy precedence, material-dependent thickness, and absorber-order
invariance. `stage2_reference_spectra`, `stage2_source_adapters`, and
`stage2_source_tool` validate pinned AMS/ACE/GOES bytes and transformations.
`stage2_published_transport` independently integrates the downloaded AMS-02
proton bins and verifies that they drive the normal case parser, exact-bin
source sampler, Geant4 transport, and TID/H100/10 outputs. Together they cover
strict schemas, dual event windows, sample
timing, quality rejection, coverage accounting, direction retention,
determinism, uncertainty-component retention, zero-count masking, and analytic
integration/statistical-error cases.

`stage3_core` checks canonical rigid transforms, exact multilayer interfaces
and mass/path formulas, order-sensitive hashes, response identity/zero/
linearity, covariance propagation, positive-semidefinite covariance rejection,
and the no-extrapolation boundary. `stage3_transport` then creates all of its
cases below `build/stage3_transport_runs/` and checks rotated/translated local
boundary scoring, independent absorber permutations, five-seed single/split
equivalence, W/HDPE and HDPE/W direct-vs-folded response, actual response-grid
refinement, cut/step/source-plane/transverse-size/event convergence, explicit
energy-ledger scope, proton/He-4/Fe-56 construction, FTFP_BERT_HP, the
G4IonQMDPhysics alternative, W/water transport, and an exact stopped-primary
zero. The registered test uses ≤2% Monte Carlo standard uncertainty, a 95%
confidence interval wholly inside ±5% for equivalence claims, and ≤2% for the
last numerical/grid refinement.
`stage3_response_weights` independently drives the production assembler,
weight builder, and folder with an analytic log-energy/`mu`/periodic-`phi`
tensor, and checks exact-node, extrapolation, and incomplete-grid rejection.

`stage4_core` checks exact MeV/kg-to-Gy and Gy-to-rad conversion, the case-v2
TID gate, the four-category all-pairs covariance calculation, rejection of an
indefinite correlation matrix and aggregate-only v2 claim, combined 95%
interval construction, gamma attenuation-node parsing, and continued rejection
of both gamma TID and DDD. `stage4_nist_reference` authenticates the pinned raw and canonical
NIST files, then makes 63 exact-node comparisons of Geant4 electronic/total
stopping power, numerically integrated CSDA range, and XCOM photon mass attenuation
for Al, Si, and water. It also records the option-4 EM sensitivity. X-ray
attenuation is never treated as finite-slab absorbed dose.

`stage4_tid_transport` writes all generated cases below
`build/stage4_tid_transport_runs/`. It runs the complete primary-generator,
multilayer Geant4, absorber-scoring, accuracy, and result path for proton and
He-4 incidence in Al, Si, and water. It proves analytic slab mass and dose
conversion for every accepted case, exact closure of particle/ancestry
partitions and diagnostic process views, ≤2% accepted-case Monte Carlo
standard uncertainty, three-level cut/step/lateral-escape convergence, option-4
TID sensitivity, and five-seed `1/sqrt(N)` behavior. Changing transverse slab
area changes its mass and therefore Gy per compact-plane primary; lateral-size
convergence is consequently tested on deposited energy while dose conversion
is checked separately at every fixed geometry.

The dependent `validation_figures` test renders every current analytic or
observational comparison after its numerical tests pass. It writes PNG, SVG,
and a checksum/caption manifest under the build tree:

```text
build/validation_figures/
├── stage1_core_analytic.{png,svg}
├── stage1_case_boundaries_analytic.{png,svg}
├── stage1_transport_smoke_analytic.{png,svg}
├── stage2_ams02_analytic_comparisons.{png,svg}
├── stage2_ams02_published_transport_analytic.{png,svg}
├── stage2_ams02_observation_adapter.{png,svg}
├── stage2_ams02_solar_cycle_observation_adapter.{png,svg}
├── stage2_ace_cris_observation_adapter.{png,svg}
├── stage2_ace_cris_observation_comparisons.{png,svg}
├── stage2_goes16_analytic_comparisons.{png,svg}
├── stage2_goes16_observation_event.{png,svg}
├── crater_observation_reference.{png,svg}
├── stage3_core_analytic.{png,svg}
├── stage3_response_weights_analytic.{png,svg}
├── stage3_transport_comparisons.{png,svg}
├── stage4_core_analytic.{png,svg}
├── stage4_nist_reference_comparisons.{png,svg}
├── stage4_tid_transport_comparisons.{png,svg}
└── manifest.json
```

The figures require Python `matplotlib`. They are human-review artifacts, not
pixel-based acceptance tests: the underlying numeric comparisons remain the
pass/fail authority. CRaTER panels are explicitly observation-only because the
instrument-response model required for a model-to-observation comparison has
not been implemented.

The standalone Layer-2 and Layer-3 runners also generate PNG/SVG review
figures after all of their prerequisite numerical checks pass. Layer 2 covers
the pencil beam, isotropic plane-crossing distribution, and near-vacuum
transmission limits. Layer 3 covers source-normalization scaling and H100/10
threshold limits; it also plots the existing controlled proton/alpha LET sanity
case, explicitly labeled as neither analytic nor observational validation.
These artifacts and their checksum/caption manifest are written to
`layer2_test_runs/validation_figures/` and
`layer3_test_runs/validation_figures/` by default. Set `FIGURE_DIR` to choose a
different destination.

## Stage-2 published-source tool

`shieldSimSourceTool` validates and canonicalizes published data without
starting Geant4. A successful conversion does not authorize transport.

```bash
./shieldSimSourceTool ams02 \
  --info=tests/reference_spectra/gcr/ams02_20131118_hapi_info.json \
  --data=tests/reference_spectra/gcr/ams02_20131118_hapi_data.csv \
  --interpretation=tests/reference_spectra/gcr/ams02_20131118_interpretation.v1.json \
  --output=ams02_canonical.json

./shieldSimSourceTool ams02-cycle \
  --interpretation=tests/reference_spectra/gcr/ams02_2025/interpretation.v1.json \
  --bartels-rotation=2581 \
  --output=ams02_br2581_canonical.json

./shieldSimSourceTool ace-cris \
  --info=tests/reference_spectra/gcr/ace_cris/ac_h3_cris_hapi_info.json \
  --data=tests/reference_spectra/gcr/ace_cris/ac_h3_cris_20191201_all.csv \
  --interpretation=tests/reference_spectra/gcr/ace_cris/interpretation.v1.json \
  --output=ace_cris_20191201_canonical.json

./shieldSimSourceTool goes16 \
  --csv=tests/reference_spectra/sep/goes16_sgps_gle72_timeseries.csv \
  --definition=tests/reference_spectra/sep/goes16_gle72_event_definition.v1.json \
  --output=goes16_gle72_fluence.json
```

The original AMS output keeps proton and inclusive mixed-isotope helium as
separate components and retains the energy-per-nucleon coordinate and total
uncertainties. The 2025 AMS output keeps He--O on the published rigidity grid,
retains all three uncertainty components, and assigns no isotope or Geant4 ion
to inclusive elemental measurements. The ACE/CRIS output retains B--Ni CDF
energy bounds, daily counts, count-derived statistical errors, and explicit
invalid masks for zero-count bins; systematic errors remain disabled pending
resolution of the V04/current-factors version mapping.
The GOES output contains one ledger for every event-window × sensor × channel,
including valid, rejected, and missing seconds. It does not combine published
uncertainties because their temporal covariance is unavailable.

The `ams02_hapi_component` case format is intentionally narrower than the
source tool. It accepts only `component_id: "proton"`, requires an isotropic
`verification_test`, and samples piecewise-constant intensity over the exact
published bin edges. Inclusive helium, AMS He--O, ACE B--Ni, and GOES are not
transport-admitted because their isotope or directional boundary information
is insufficient. See `examples/README.md` for the runnable validation case.

`data/bon2020_v2_3/README.md` records a retired intake proposal and explains
why the unverified BON2020-named legacy files remain prohibited. No BON2020
software or data are required for Stage 2 or planned as a production runtime
dependency.

The complete acquisition recipe for the new GCR observations is:

```bash
python3 scripts/download_gcr_validation_data.py
```

It verifies reviewed SHA-256 values and refuses changed upstream or local
bytes. CTest is offline and uses the committed immutable copies.

## Versioned cases

Validate and serialize the example without starting transport:

```bash
./shieldSim --validate-case=examples/stage1_case.json
```

Run its two independent absorbers:

```bash
./shieldSim --case=examples/stage1_case.json
```

Outputs use stable absorber IDs, for example
`stage1_example__water_proxy_*` and `stage1_example__silicon_proxy_*`. Every
absorber receives the same shield-exit environment in its own rebuilt world;
absorbers are never placed behind one another.

The complete annotated case is in `examples/stage1_case.json`. Its key
contracts are:

- Every physical number carries a unit. Accepted Stage-1 length units are
  `mm`, `cm`, and `m`; layer areal density is exactly `g/cm2`; spectrum energy
  is `MeV` total kinetic energy per particle.
- `material_library.library_id` and `library_version` must match the loaded
  file. A case override uses a new `material_id` and pins its parent as
  `<id>@<revision>`.
- `geometry.shield_stack.layers` is ordered upstream to downstream.
- Every absorber is currently `type: "slab"` and explicitly gives transverse
  size, thickness, and front-z position.
- Omitting only `accuracy.target` selects and archives `0.10`. A per-quantity
  override must replace target, scope, and coverage together.
- `min_histories`, `max_histories`, and `batch_size` are mandatory resource
  controls. An incomplete combined systematic budget produces
  `not_assessable`, never a false claim that the requested accuracy was met.
- A numeric spectrum is accepted only with a
  `shieldSim.spectrum_metadata.v1` sidecar. Rows must already be strictly
  increasing and finite; the program never sorts, clamps, or discards them.
- A downloaded AMS spectrum uses `format: "ams02_hapi_component"` with pinned
  info, data, and interpretation paths plus `component_id: "proton"`. This
  path is validation-only; the exact contract and runnable example are in
  `examples/stage2_ams02_proton_validation_case.json`.

### Stage-3 transformed response cases

`examples/stage3_response_cell_case.json` is the complete `shieldSim.case.v2`
example. Every transform is active local-to-world and contains a translation
with length units plus a scalar-first unit quaternion `(w,x,y,z)`. Local `+z`
is the surface normal and direction nodes use `mu=cos(theta)` in `(0,1]` plus
periodic `phi`. Quaternions must be normalized; `q` and `-q` serialize to one
canonical representation. The source, stack, and absorber transforms are
independent input fields, but validation requires the source upstream, the
absorber downstream, all oriented boxes inside the world, and no shield/
absorber positive-volume overlap.

Run one response cell:

```bash
./shieldSim --case=examples/stage3_response_cell_case.json
```

The response-cell source is exactly one declared species, energy, and angular
node per history. `particle` must be proton `(1,1)`, alpha/He-4 `(2,4)`, gamma
`(0,0)`, or `ion` with explicit positive `Z,A`. Gamma is verification-only and
requires total kinetic energy; gamma and explicit-ion nodes may request only
`transport_raw`, not the Stage-4 `TID` claim. Charged-particle energy is explicitly total `MeV`/`GeV` or per
nucleon `MeV/n`/`GeV/n`; Geant4 always receives the resolved total kinetic
energy. The case also declares the reference physics list, EM configuration,
either its standard ion configuration or `G4IonQMDPhysics`, range cut, step
limit, seed, and fixed or adaptive history plan.

For a grid, run one case per incident node and assemble only compatible cells:

```bash
./shieldSimResponseTool assemble response_matrix.json cell_*.json
python3 scripts/build_response_weights.py response_matrix.json source_nodes.json weights.json
./shieldSimResponseTool fold response_matrix.json weights.json folded.json
```

The JSON matrix manifest is authoritative and its long-form CSV derivative is
checksum-recorded. A response key includes the complete geometry/material
hashes, `slab.v1` provider, scoring contract, Geant4/code/build identity,
physics/ion settings, cuts/step, and Geant4 data-library paths. Seed, history
count, output path, and incident node are cell metadata and do not change that
key. Assembly rejects mismatched keys and forged key hashes. Folding requires
normalized explicit weights, propagates matrix and optional source covariance,
rejects a non-positive-semidefinite covariance, and never extrapolates.
Each cell carries the complete upper triangle of event-sample covariance for
its primitive score means. Folding propagates that matrix covariance and the
optional source-weight covariance into a complete folded primitive-score
covariance, so a later nonlinear quantity can retain numerator/denominator
correlation. `build_response_weights.py` supports exact nodes or a declared
tensor-grid
policy that is linear in response versus log-energy and bilinear in `mu` and
periodic `phi`; it rejects incomplete grids and out-of-range requests.

`discrete_response_nodes` is the independent direct-sampling cross-check. Its
node weights must sum to one and use the same exact species/energy/direction
definitions. A nonlinear derived result must never be folded directly: fold
its primitive numerator and denominator, including their covariance, and form
the ratio only in the later quantity implementation.

### Stage-4 finite-slab TID cases

`examples/stage4_tid_case.json` is the complete case-v2 TID example. Its
default 10% combined target is explicit. The example's zero material/source
terms apply only to its ideal analytic slab and exactly one-primary
verification normalization; they are not defaults for an as-built mission
case. Numerical and physics-reference values are linked to the Stage-4
acceptance evidence. Production inputs must replace every component with
case-specific evidence.

A complete structured budget supplies components in the `material_mass`,
`source`, `numerical`, and `physics_reference` categories. Every component has
a stable ID, relative standard uncertainty, and provenance. It also supplies
all unordered component pairs under `correlations`, including explicit zeros,
with a coefficient and rationale, and sets `pair_policy` to
`all_pairs_explicit`. The parser rejects missing or duplicate pairs, unknown
categories, coefficients outside `[-1,1]`, and a non-positive-semidefinite
correlation matrix. The legacy single aggregate percentage is rejected for a
complete case-v2 TID claim.

Combined v2 TID selects
`student_t_mc_plus_normal_systematics_quadrature`; both the MC-only and combined
95% intervals are written. An `incomplete` budget may instead report a known
floor and yields `not_assessable` for combined accuracy. See the example for a
fully expanded four-component/six-correlation record.

### Adding a future arbitrary absorber shape

`include/AbsorberGeometry.hh` is the extension boundary. Transport currently
constructs only `SlabAbsorberGeometry`, but scorers consume the provider's
stable region IDs and logical regions rather than a `G4Box` name or copy
number. A provider must implement `Validate`, `Build`, `SensitiveRegionIds`,
`MassProperties`, `BoundsAndReferenceFrame`, and
`CanonicalGeometryDescription`; its geometry hash is part of every response
key. This keeps source generation, ancestry/scoring, response matrices, and
quantity algorithms independent of whether a future implementation uses CSG,
CAD/mesh, voxels, GDML, or an application-specific spacecraft model.

A future provider must return stable sensitive-region IDs, logical volumes,
analytic or independently verified volumes/masses/centroids, material hashes,
bounds/reference frame, provider version, source-file hashes, and a canonical
geometry description. It must preserve explicit source units, scale and
placement transform, region-to-material mapping, sensitive-region mapping,
and any overburden. A geometry/material/region/placement change must invalidate
the response cache even when nominal mass is unchanged.

Before it can be enabled, such a provider must reject missing or ambiguous
units, non-finite vertices, non-watertight or self-intersecting solids,
inconsistent surface orientation, overlaps, unmapped regions, empty sensitive
regions, and mass/volume disagreement. Required tests include primitive
analytic volumes and chords, analytic/mesh/voxel convergence, rotation and
translation behavior, reference-ray boundary navigation, region ledger
closure, harmless tessellation invariance, and direct-vs-response agreement.
Automatic mesh repair is prohibited unless it is an explicit versioned
transformation with input/output hashes. See Section 3.3.4 and Stage 10 of
`development_and_validatoin.md` for the full gate.

## Legacy compatibility examples

The following CLI is retained for regression and diagnostic work. It uses the
historical single-shield/default-geometry model and is not the versioned
science-facing interface.

Single 2 mm Al shield with the built-in historical GCR-like approximation and
the default physics list (regression compatibility only; not a physical GCR
environment):

```bash
./shieldSim --source-mode=beam --shield=Al:2 --events=50000
```

Isotropic source over the upstream plane using high-precision neutron transport:

```bash
./shieldSim --physics-list=FTFP_BERT_HP --source-mode=isotropic \
            --shield=Al:2 --events=50000
```

Tabulated SEP-like source spectrum:

```bash
./shieldSim --source-mode=isotropic --spectrum=examples/sep_spectrum.dat \
            --legacy-spectrum-metadata=examples/sep_spectrum.metadata.v1.json \
            --shield=Al:2 --events=100000
```

Dose-vs-thickness sweep over HDPE:

```bash
./shieldSim --sweep --source-mode=isotropic --sweep-material=HDPE \
            --sweep-tmin=0.5 --sweep-tmax=30 --sweep-n=15 \
            --sweep-log --events=20000
```

## Material library and overrides

Versioned cases load `data/material_library.v1.json`. The runnable statuses are
`reference`, `trade_study`, and fully populated `project_specific` records.
`template` records are rejected before geometry construction. In particular,
generic CFRP, SiC composite, lunar regolith, and Mars regolith names are
disabled templates: their composition and density cannot be invented safely.

The initial library includes explicit Al, Cu, W, Ta, natural B, HDPE, 5 wt%
natural-boron polyethylene, dry Kevlar fiber, Kapton, water, silicon, skin,
eye-lens, soft-tissue/BFO/CNS composition slabs, SiO2, 6H-SiC, GaAs,
Ga0.47In0.53As, and Ge records. Each record identifies whether it is a
reference or preliminary trade recipe and states its applicability limits.
Unreported property uncertainties remain explicitly `not_provided`; the code
does not manufacture them.

A case-local override must use a new stable ID, pin its exact parent revision,
state `trade_study` or `project_specific`, and cite the changed property. The
parent stays unchanged. Density changes alter mass-thickness conversion and
geometry hashes; composition changes alter the material hash and invalidate
dependent response data. All transport-relevant properties in a complete
override—including density, state, temperature, composition, and provenance—
are input data rather than compiled constants.
Natural abundances use the pinned Geant4 element data; an override can instead
declare explicit isotope components such as `B-10`/`B-11`. Gas records require
pressure, and an explicitly overridden mean excitation energy is applied to the
Geant4 material. Unsupported thermal-scattering-law associations are rejected
rather than recorded without transport effect.

## Legacy shield-material selection

Use `--shield=<material>:<thickness_mm>` for a single run, or `--sweep-material=<material>` in sweep mode.  The material name can be either a normal Geant4/NIST material such as `G4_Al`, `G4_WATER`, or `G4_Si`, or one of the built-in keys below.

Print the full material table, including aliases, composition notes, and reference notes, with:

```bash
./shieldSim --list-materials
```

The same table is also printed inside:

```bash
./shieldSim --help
```

Built-in shielding keys:

```text
Structural Metals
  Al               Aluminum (Al), alias to G4_Al
  Cu               Copper (Cu), alias to G4_Cu
  W                Tungsten (W), alias to G4_W
  Ta               Tantalum (Ta), alias to G4_Ta

Hydrogenous Polymers
  HDPE             Custom polyethylene approximation, CH2, density 0.95 g/cm3
  BPE              Custom borated polyethylene, 95 wt% HDPE + 5 wt% B
  Kevlar           Custom aramid/Kevlar approximation, C14H10N2O2
  Kapton           Polyimide/Kapton, alias to G4_KAPTON

Composites
  CFRP             Custom 60 wt% carbon fiber + 40 wt% epoxy approximation
  SiCComposite     Custom 70 wt% SiC + 30 wt% epoxy/plastic approximation

Water / Regolith
  Water            Water, alias to G4_WATER
  LunarRegolith    Approximate bulk lunar regolith
  MarsRegolith     Approximate bulk Martian regolith
```

Examples:

```bash
./shieldSim --shield=HDPE:10 --events=50000
./shieldSim --shield=BPE:10 --physics-list=Shielding --events=50000
./shieldSim --sweep --sweep-material=CFRP --sweep-tmin=1 --sweep-tmax=20 --sweep-n=10
./shieldSim --shield=LunarRegolith:50 --source-mode=isotropic --events=100000
```

These names describe the compatibility-only compiled catalog. Its historical
generic composite/regolith recipes are not imported into the versioned
material library. Do not edit compiled source to represent a flight material;
create a cited `project_specific` case record instead.


## Legacy detector / absorber / target material selection

Use `--scoring=<material>:<thickness_mm>,...` to define downstream scoring slabs.  The alias `--target=<material>:<thickness_mm>,...` is equivalent and is often clearer when the slabs represent tissue, detector, or electronics media.

Print the full detector/target catalog with:

```bash
./shieldSim --list-target-materials
```

Built-in detector/target keys:

```text
Crewed Mission - Tissue Targets
  Skin             Skin/epidermis proxy; use e.g. Skin:0.07 for 0.07 mm
  EyeLens          Crystalline eye-lens tissue proxy
  BFO              Blood-forming-organ soft-tissue proxy; BFO:50 for 5 cm
  CNS              Brain/CNS tissue proxy
  SoftTissue       ICRU-style 4-element soft tissue / muscle proxy

Electronics - Silicon Family
  Si               Silicon, alias to G4_Si
  SiO2             Silicon dioxide / gate-oxide proxy
  SiC              Dense stoichiometric silicon carbide

Electronics - III-V and Compound Semiconductors
  GaAs             Gallium arsenide
  InGaAs           In0.53Ga0.47As proxy
  Ge               Germanium, alias to G4_Ge
  H2O              Water reference, alias to G4_WATER
```

Examples:

```bash
# Tissue and electronics scoring behind 2 mm aluminum
./shieldSim --shield=Al:2 --target=BFO:50,Si:1 --events=100000

# Shallow skin and eye-lens scoring
./shieldSim --shield=HDPE:5 --target=Skin:0.07,EyeLens:1 --source-mode=isotropic

# Electronics detector materials
./shieldSim --shield=Al:2 --target=Si:1,SiO2:0.01,GaAs:1,InGaAs:1,Ge:1
```

Organ/tissue names define material composition only. `shieldSim` does not apply
NASA limit checks, quality factors, LET/RBE weighting, organ weighting, or
device response functions internally.

## Stage-5 silicon DDD and neutron-equivalent fluence

Use the versioned case-v2 path for the implemented damage calculation:

```bash
./shieldSim --validate-case=../examples/stage5_ddd_neq_case.json
./shieldSim --case=../examples/stage5_ddd_neq_case.json
```

The case must request `DDD`, `n_eq`, or both and must provide a
`damage_response` object naming `shieldSim.silicon_displacement_damage` version
`1.0.0`. Every absorber in that case must resolve to the authenticated natural-
silicon material identity. The manifest and each table are hash-checked before
transport; the fully expanded manifest and its canonical hash are retained in
the resolved case/result provenance.

For slab entry area `A`, inward surface normal `n`, and crossing direction `u`,
the authoritative event scores are

```text
DDD  = sum_crossings S_NIEL(species, E_entry) / (A |n dot u|)
n_eq = sum_neutron_crossings k_E722(E_entry) / (A |n dot u|)
```

where re-entry through the declared upstream surface is another physical
fluence crossing. The projected-cosine factor converts crossings of the fixed
face into scalar fluence and is tested with an oblique analytic case. Side
entry and downstream-face backscatter are excluded by
the `slab.v1` surface definition. The surface classifier and normalization area
belong to the `AbsorberGeometry` provider, not to the quantity calculation. A
future arbitrary-shape provider must define and hash its entrance patch(es),
orientation, area/weighting rule, sensitive-region mapping, and boundary
tolerance, then pass the shape-validation ladder in the roadmap. It can reuse
the same transport hook and damage fold; syntax acceptance alone will not
enable a new shape.

The current table domains are 0.0002–10000 MeV for protons, 0.0001–4000 MeV
for helium-4, and 1e-10–20 MeV for neutrons. Energies are total kinetic energy
per particle. Charged-particle values are interpolated log-log between positive
nodes. Neutron values are group constants on explicit half-open ASTM intervals,
with the final 20 MeV endpoint included. No extrapolation is allowed. Any
nonzero missing contribution makes the corresponding output and its accuracy
`not_assessable`, while the supported partial value and rejected-crossing count
remain auditable.

`n_eq` includes neutron damage only and uses the ASTM E722-19 95 MeV-mb
natural-silicon reference. Proton and helium-4 DDD are not converted to
neutron-equivalent fluence. The result’s `diagnostic_Geant4_nonionizing_deposition`
is a separate finite-slab diagnostic and is never summed with DDD.

The usual 10% default accuracy request still applies. Statistical precision can
be assessed from event histories. A `combined` Stage-5 claim additionally
requires a case-supplied `physics_reference` standard-uncertainty component
with explicit provenance; the code does not invent an uncertainty for the
SR-NIEL/ASTM response. Without it, combined accuracy is `not_assessable`.
See [`data/niel/README.md`](data/niel/README.md) and
[`tests/reference_ddd/README.md`](tests/reference_ddd/README.md) for data and
test provenance. This is an environmental exposure proxy, not a device-
degradation or instrument-design calculation.

## Stage-6 electronic LET and spectral hardness

The case-v2 example exercises both implemented Stage-6 quantities:

```bash
./shieldSim --validate-case=../examples/stage6_let_hardness_case.json
./shieldSim --case=../examples/stage6_let_hardness_case.json
```

Every LET case must declare `let_response`. The response ID is
`shieldSim.electronic_let.v1`; bins are logarithmic in the primary unit
`keV/um`, and the only supported out-of-range policy is
`not_assessable_no_clamping`. The recommended example spans 0.001--10000
keV/um with 20 bins/decade. Resolved bin edges are included in the case hash.

For an absorber of volume `V`, every positive-length charged step contributes
`dx/V` to its LET bin. Geant4 unrestricted electronic stopping is evaluated at
the arithmetic midpoint of the transported pre/post kinetic energies. Track
and dose means are calculated from event-level sufficient statistics:

```text
L_track = sum(L dx)   / sum(dx)
L_dose  = sum(L^2 dx) / sum(L dx)
```

The total spectrum includes every transported charged particle. Species, PDG,
charge, Z/A, primary/secondary identity, root particle, and generation are
preserved in component spectra; delta electrons remain electron components.
Zero-length charged steps are counted but add no fluence. Ionizing deposited
energy is still the separate TID ledger and is not used as a substitute for
electronic stopping. Result JSON includes bin standard errors, component and
aggregate means, event moment/covariance records, under/overflow fluence, and
the density conversion `1 keV/um = 10/rho MeV cm2/g`.

`H100/10` is evaluated at the shield rear face from primary protons only:

```text
H100/10 = mean_event(sum crossings[E>100] 1/(A |n dot u|))
          / mean_event(sum crossings[E>10] 1/(A |n dot u|))
```

Thresholds are strict and use exact total kinetic energy. Every accepted
forward crossing contributes; repeated crossings remain correlated within the
same primary event. The result retains numerator/denominator sums, squares,
products, covariance of means, the ratio point estimate, and a two-sided 95%
Fieller interval. A zero denominator or a denominator whose Fieller confidence
set includes zero is emitted as JSON `null`, not zero or infinity.

The stopping-reference domain is proton/He-4 incidence in natural aluminum,
natural silicon, or water using the pinned NIST tables. The selected heavier-
ion screen is fully stripped C-12 in natural silicon at 10, 100, and 1000
MeV/u against GSI ATIMA v1.41. C-12 requires
`G4EmStandardPhysics_option4`: the reference-list default is rejected because
it misses the fixed 10 MeV/u 15% screen. That screen is not an experimental
uncertainty or device/instrument validation. A combined Stage-6 accuracy claim
requires a complete structured budget containing a provenance-bearing
`physics_reference` term; otherwise combined accuracy is `not_assessable`.
See [`tests/reference_let/README.md`](tests/reference_let/README.md) and the
Stage-6 roadmap section for the exact evidence boundary.


## Legacy computed-quantity selection

This section documents compatibility output only. The analytic DDD, non-silicon
`n_eq`, and transmitted-spectrum LET implementations do not satisfy the
production definitions in `development_and_validatoin.md` and are disabled in
`shieldSim.case.v1`.

The code now computes the quantities shown in the UI for every selected `shield x target` combination.  By default all computed quantities are enabled:

```bash
./shieldSim --shield=Al:2 --target=BFO:50,Si:1 --quantities=all --events=100000
```

You can restrict the output with a comma-separated list:

```bash
./shieldSim --shield=Al:2 --target=Si:1 --quantities=TID,DDD,n_eq,H100/10
./shieldSim --shield=Al:2 --target=Si:1 --quantities=TID,LET
./shieldSim --list-quantities
```

Implemented quantities:

```text
TID       Total Ionizing Dose
          D_ion = (E_dep,total - E_dep,nonionizing) / m in the finite target slab.
          The total and non-ionizing ledgers are retained for closure.
          Output units: Gy/primary, rad/primary, Gy/s, rad/s.

DDD       Displacement Damage Dose proxy
          D_d = integral Phi(E) NIEL(E) dE.
          Phi is the area-averaged transmitted fluence at the downstream shield face.
          NIEL is a documented analytic surrogate in src/ComputedQuantities.cc.
          Output units: MeV/g/primary and MeV/g/s, plus rad-equivalent columns.

n_eq      1-MeV Neutron Equivalent fluence proxy
          n_eq = DDD / NIEL_1MeV_neutron(material).
          Historical surrogate only; it is distinct from the case-v2 neutron-only ASTM E722 implementation.
          Output units: cm^-2/primary and cm^-2/s.

LET       LET spectrum
          dPhi/dLET = integral dPhi/dE delta(LET - LET(E,material)) dE.
          LET(E,material) is estimated with a Bethe-Bloch electronic mass-stopping-power approximation.
          Output units: area-averaged fluence per primary and fluence rate per LET bin.

H100/10   Spectral hardness index
          H100/10 = primary-proton shield-exit fluence(E>100 MeV)
                    / primary-proton shield-exit fluence(E>10 MeV).
          Exact crossing energies are used; secondaries, alphas, and neutrons are excluded.
          A zero denominator is undefined, not zero.
```

The historical equations are documented in `include/ComputedQuantities.hh` and
`src/ComputedQuantities.cc`; their presence is not validation evidence. They
are not used by the versioned Stage-5/6 calculations above. Formal sign-off
still requires the external evidence/review stated in the roadmap.

## Legacy compiled catalog maintenance

The material-property references used by the built-in catalogs are documented in `src/MaterialCatalog.cc` and summarized by `./shieldSim --list-materials` and `./shieldSim --list-target-materials`.  In brief:

- `Al`, `Cu`, `W`, `Ta`, `Kapton`, and `Water` are aliases to Geant4/NIST `G4_*` materials.
- `HDPE`, `BPE`, `Kevlar`, `CFRP`, `SiCComposite`, `LunarRegolith`, and `MarsRegolith` are custom ShieldSim trade-study approximations.
- Composite and regolith definitions are intentionally approximate. Replace their density and mass fractions with measured project-specific values before using them for final quantitative analysis.
- Tissue/organ target definitions are simplified material proxies. They are not anatomical phantoms and do not include biological weighting factors.
- Electronics target definitions are transport media only. The built-in DDD, n_eq, and LET outputs are documented engineering proxies; charge collection, dark-current response, threshold shifts, and device-specific damage response still require external post-processing.

To add a new material:

1. Edit `src/MaterialCatalog.cc` and add a `MaterialCatalogEntry` in `ShieldMaterialCatalog()` for a new shield/absorber material, or in `DetectorMaterialCatalog()` for a new detector/target scoring material.
2. If the material exists in Geant4/NIST, set `canonicalName` to the corresponding `G4_*` name and set `isCustom=false`.
3. If the material is custom, implement a `Build...()` function in `src/MaterialCatalog.cc`, use explicit Geant4 units such as `g/cm3`, register the builder in `BuildCustomMaterial()`, and document the density/composition/reference above the builder.
4. Add useful aliases. The parser ignores case and punctuation, so aliases can be user-friendly names.
5. Rebuild. The `--help`, `--list-materials`, and `--list-target-materials` output is generated from the catalogs automatically.

## Physics-list selection

Use `--physics-list=<name>` to select the Geant4 reference physics list at runtime. Supported values are:

```text
FTFP_BERT
FTFP_BERT_HP
Shielding
QGSP_BIC_HP
```

The default is `FTFP_BERT`. For shielding calculations where low-energy secondary neutrons can affect dose behind the absorber, compare against one of the high-precision neutron options, especially `FTFP_BERT_HP`, `QGSP_BIC_HP`, or `Shielding`.

Example:

```bash
./shieldSim --physics-list=Shielding --source-mode=isotropic \
            --shield=Al:5 --events=100000
```

## Numerical controls for verification

The default behavior is to use the selected Geant4 physics list with its normal
production cuts and stepping behavior.  For tests and numerical-convergence
studies, the CLI also provides:

```bash
--production-cut=<mm>
--max-step=<mm>
```

`--production-cut=<mm>` calls `SetDefaultCutValue()` on the selected Geant4
physics list.  This controls the range cut used for explicit secondary
production.  It can affect low-energy secondaries, TID, DDD, `n_eq`, and LET
tails, so production studies should record the selected value.

`--max-step=<mm>` attaches `G4UserLimits` to the shield and scoring slabs and
registers `G4StepLimiterPhysics`.  This is useful for thin-target dose checks
and LET-spectrum convergence studies.  Omit it to use the normal physics-list
stepping behavior.

These options are primarily intended for Layer-3 tests and convergence studies,
not as substitutes for physics validation.

## Legacy spectrum adapter

A tabulated spectrum file has three columns and now requires an explicit
`--legacy-spectrum-metadata=<file.json>` sidecar:

```text
E[MeV]   protonSpectrum   alphaSpectrum
```

`E` is total kinetic energy per particle. For alpha particles, `E` is total alpha kinetic energy, not MeV/nucleon.

The sidecar declares that coordinate, the energy and value units, value type,
angular model, bin representation/edge rule, science status, and provenance.
Single-row spectra must declare a bin width. Unknown fields, malformed rows,
negative/fill/NaN values, duplicate or descending energies, and an angular-mode
conflict are fatal. This adapter is for synthetic tests and explicitly
qualified user data. Published AMS-02 and GOES-16 products and future AMPS
products require their dedicated adapters; relabeling them—or an unverified
BON2020-named legacy file—as legacy columns is prohibited.

In `beam` mode, spectrum columns are interpreted as differential pencil-beam source rates, for example:

```text
particles / s / MeV
```

In `isotropic` mode, spectrum columns are interpreted as differential intensities, for example:

```text
particles / cm2 / s / sr / MeV
```

The code applies the `pi` angular factor for isotropic plane-crossing flux.

## Output files

### `<prefix>__<absorber>_resolved_case.json`

Canonical input/provenance manifest written before transport (and by
`--validate-case`). It contains the resolved 0.10 default when applicable,
both geometric and areal layer thicknesses, full resolved material records,
source metadata and hash, geometry and absorber hashes, compiler/build/source
revision, Geant4 version/data paths, seed, and capability status.

### `<prefix>__<absorber>_result.json`

Canonical result written after each adaptive batch, with
the final file representing all cumulative histories. It records quantity
definition, value/unit/normalization, deposition ledgers or H100/10 numerator
and denominator, statistical and systematic uncertainty state, 95% interval,
accuracy status/reason, and a capability-specific validation status. Versioned
Stage-6 H100/10 uses a covariance-aware Fieller 95% interval; the legacy case-v1
path retains its earlier compatibility behavior. JSON `null`, rather than zero
or infinity, represents an undefined value. A case-v2 run writes
`shieldSim.result.v2` with `transport_raw` moments, ancestry-qualified primary/
secondary rear-face crossings, per-layer deposit/track/step and secondary-
birth ledgers, absorber deposition moments, physics/EM/ion provenance, and an
explicit statement that nuclear rest mass/Q values are outside the kinetic-
energy ledger scope. For Stage-4 TID it additionally records slab mass,
Gy/primary, rad/primary, source-normalized rate, statistical and combined 95%
intervals, the complete systematic component/correlation records, and the
`stage4_finite_slab_tid_implemented_and_acceptance_tested_pending_independent_review`
status. Particle/ancestry deposition partitions are authoritative. Creator-
process and step-defining-process partitions are labeled diagnostic only; a
process name is not treated as a unique causal allocation of deposited energy.
For Stage 6 it additionally records LET response edges and out-of-range policy,
per-bin event uncertainty, charged-particle component metadata, track/dose mean
ratio moments and intervals, zero-length-step audit counts, and H100/10 weighted
fluence moments/covariance.

### `<prefix>__<absorber>_response_cell.json`

Written only for `monoenergetic_response_cell`. It stores every raw score as
the history count, sum, sum of squares, mean, sample variance, standard error,
and unit; zero-scoring histories are included in the denominator. The complete
score-covariance upper triangle stores sample covariance, covariance of the
mean, and sum of event products. It also stores the exact incident coordinate,
layer-resolved secondary ledger, and
hashed response-key document. Raw response-cell scores retain the Stage-3
acceptance status; requesting Stage-4 TID adds the separately gated TID output
to the result file and does not turn every raw score into a validated proxy.

### `shieldSim_spectra.dat`

Tecplot-format differential spectra. Energy is MeV total kinetic energy per particle.

- `*_MC` columns: raw Monte Carlo spectra in `counts/(MeV primary)`.
- `*_Norm` columns: source-normalized spectra.

For beam mode, normalized units are `particles/s/MeV` if the input spectrum is `particles/s/MeV`.

For isotropic mode, normalized units are `particles/(cm2 s MeV)` if the input spectrum is `particles/(cm2 s sr MeV)`.

### `shieldSim_quantities.dat`

Tecplot-format compatibility scalar table. A versioned case creates one file
per independent absorber. In legacy sweep mode a new zone is appended for each
shield thickness.

- TID columns: `Gy/primary`, `rad/primary`, `Gy/s`, `rad/s`.
- DDD columns: `MeV/g/primary`, `MeV/g/s`, and rad-equivalent conversions.
- `n_eq` columns: `cm^-2/primary` and `cm^-2/s`.
- `H100_10`: dimensionless hardness index.

### `shieldSim_let_spectrum.dat`

Tecplot-format LET spectra, one zone for each target material and shielding configuration when LET is enabled.  LET is in `MeV cm2/mg`.  Per-primary columns are area-averaged fluence per source primary per LET unit; rate columns multiply by the source-primary rate.

### `shieldSim_dose_sweep.dat`

Written only in sweep mode.

- `Dose_*_perPrimary` columns are `Gy/primary` when TID is enabled.
- `DoseRate_*` columns are `Gy/s` when the source spectrum has physical units.
- DDD, `n_eq`, and `H100_10` sweep columns are also included when those quantities are enabled.

### `shieldSimOutput*.csv`

Geant4 analysis histograms for exit energies of protons, alphas, and neutrons.

### Run-summary diagnostic file

Written only when requested with:

```bash
--dump-run-summary=<file>
```

This is a machine-readable diagnostic used by the Layer-3 tests. It repeats
integrated quantities such as source normalization, transmitted counts, TID,
DDD, `n_eq`, and H100/10 in simple keyword rows. It also includes H100/10
numerator/denominator counts and total/ionizing/non-ionizing deposit ledgers.
For a versioned run, the `_resolved_case.json` manifest is the authoritative
input/provenance record; output validation status must be considered before use.

## Notes and limitations

- The default physics list is `FTFP_BERT`, and the runtime option `--physics-list=` supports `FTFP_BERT`, `FTFP_BERT_HP`, `Shielding`, and `QGSP_BIC_HP`.
- The compiled built-in catalog exists for legacy regression only. Versioned
  cases use `data/material_library.v1.json` and immutable overrides.
- The isotropic source uses a finite upstream plane, not an infinite half-space source. Very oblique trajectories can interact with the finite side boundaries.
- The 5 cm × 5 cm detector is a legacy default only. Versioned cases require
  every transverse dimension. Dose per primary depends on finite absorber mass.
- Legacy DDD/n_eq use an analytic NIEL surrogate and legacy LET uses a
  transmitted-spectrum Bethe-Bloch estimate. They are disabled in versioned
  cases and must not be used for final electronics-damage analysis.

## Testing and Diagnostic Options

The package includes automated test scripts under `tests/`:

```bash
tests/run_layer1_tests.sh   # build and CLI smoke tests
tests/run_layer2_tests.sh   # geometry, source, and shield-exit scoring tests
tests/run_layer3_tests.sh   # physics-output and numerical-sanity tests
```

Layer-1 tests check build/CLI behavior.  Layer-2 tests check geometry, source
sampling, and shield rear-face scoring using short deterministic runs.  Layer-3
tests check source-normalization scaling, computed-quantity output, H100/10,
LET folding, physics-list selection, sweep output, and optional numerical
controls.

Diagnostic CLI options used by the tests are:

```bash
--random-seed=<integer>
--output-prefix=<name>
--dump-source-samples=<file>
--dump-exit-particles=<file>
--dump-run-summary=<file>
--diagnostic-max-rows=<n>
--production-cut=<mm>
--max-step=<mm>
```

`--dump-source-samples` writes one row per primary generated by
`PrimaryGeneratorAction`:

```text
row species E_MeV x_mm y_mm z_mm ux uy uz
```

`--dump-exit-particles` writes one row for each particle accepted by
`SteppingAction` as crossing the downstream shield face:

```text
row species E_MeV xg_mm yg_mm zg_mm xl_mm yl_mm zl_mm uxl uyl uzl
```

`--dump-run-summary` writes one machine-readable scalar block per run:

```text
scalar H100_10 <value>
count output_proton <value>
target index name thickness_mm TID_Gy_perPrimary TIDRate_Gy_s DDD_MeV_g_perPrimary DDDRate_MeV_g_s n_eq_cm2_perPrimary n_eq_rate_cm2_s
```

Run the tests from the top-level package directory:

```bash
chmod +x tests/run_layer1_tests.sh tests/run_layer2_tests.sh tests/run_layer3_tests.sh
tests/run_layer1_tests.sh
tests/run_layer2_tests.sh
tests/run_layer3_tests.sh
```

Each script accepts `-h`, `-help`, and `--help`:

```bash
tests/run_layer3_tests.sh --help
```

Layer-3 runtime can be increased for higher-statistics checks:

```bash
EVENTS_SHORT=10000 EVENTS_SMOKE=1000 tests/run_layer3_tests.sh
```

Generated build/test directories and logs are ignored by `.gitignore`.
