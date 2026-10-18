# shieldSim Test Suite

This directory contains automated and semi-automated test scripts for the
`shieldSim` Geant4 shielding application.

The tests are organized in layers. Lower-numbered layers test basic software
functionality, while higher-numbered layers test geometry, source sampling,
scoring, physics-output behavior, numerical controls, and regression stability.

## Stage-1 contract tests

The CMake/CTest suite is the acceptance entry point for the versioned Stage-1
core:

```bash
cmake -S . -B build -DGeant4_DIR=/path/to/geant4/lib/cmake/Geant4
cmake --build build -j
cd build
ctest --output-on-failure
```

`stage1_core_tests.cc` runs without transport and checks:

- the SHA-256 test vector and strict duplicate-key/NaN rejection;
- the explicit SiO2 and Ga0.47In0.53As corrections;
- material status gates, constituent-to-element conversion, immutable
  override precedence, and content-hash invalidation;
- mandatory legacy-spectrum metadata and declared one-bin width;
- areal-density-to-geometric-thickness conversion;
- the exact 0.10 default accuracy; and
- statistical `1/sqrt(N)` behavior and `not_assessable` combined status when
  the systematic budget is incomplete.

`test_stage1_case_boundaries.py` invokes the production executable with
`--validate-case` and checks:

- one and five layers accepted, six rejected;
- 100.000 mm accepted and 100.001 mm rejected per layer;
- missing dimensions, unknown units/fields, reversed energy ranges, invalid
  accuracy, and Stage-1-unsupported quantities rejected;
- global/per-quantity accuracy precedence in serialized provenance;
- absorber permutation preserves each absorber's independent geometry hash;
  and
- a density override changes derived thickness, material hash, and geometry
  hash.

`test_stage1_transport_smoke.py` performs short Geant4 runs in both absorber
orders. It verifies independent-world order invariance, resolved Geant4
material densities and analytic slab masses, and exact closure of
`total = ionizing + non-ionizing` deposition ledgers.

These tests validate schemas and exact analytic relationships. The short
Geant4 scripts below remain useful regression/sanity tests, but neither kind of
test substitutes for later comparison to analytic transport limits, published
stopping/NIEL data, or observations.

Published-data integrity tests are kept separate from the existing Geant4
layers:

```bash
python3 tests/test_reference_spectra.py
```

They verify checksum-pinned AMS-02 GCR, ACE/CRIS GCR, GOES-16 SGPS SEP, and
CRaTER LET fixtures under `tests/reference_spectra/`. These checks establish
data provenance and retained semantics. The separate
`stage2_published_transport` test establishes that one admitted published
proton component reaches the normal Geant4/proxy path. Neither result is a
claim of observational agreement or complete mission-environment coverage.
See `development_and_validatoin.md` for the staged validation gates.

The separate Stage-8 input-acquisition integrity gate is:

```bash
python3 tests/test_crater_reference_inputs.py
```

It checks the package under `tests/reference_crater/`: the CRaTER instrument
and method documents, exact published background cuts, PDS record/quality
metadata and processing source, and time-matched GLACE/BON2020 H--Ni spectra.
This test does not transport particles or compare shieldSim output with the
observation. Its acquisition record deliberately keeps the missing engineering
mass model, calibration plan, bin covariance, preregistration, and reviewer
sign-off visible; passing it must not be reported as Stage-8 completion.

The bounded Stage-8 benchmark gates are:

```bash
ctest -R '^stage8_(reference_inputs|crater_core|crater_transport)$' --output-on-failure
```

`stage8_crater_core` protects the publication-derived 19-entry axial stack,
132.02606 mm total length, six detector thicknesses, GLACE archive-derived
H--Ni isotope map, analytic log--log interpolation, D5/D6 splice, and exact
published cut logic. `stage8_crater_transport` divides one 560-history smoke
case between two concurrent processes, verifies every shard, merges raw
importance-weight sums and sums of squares, recomputes the global uncertainty,
and renders PNG/SVG observation comparisons. It asserts transport, accounting,
merging, and labeling—not observational agreement.

## Stage-2 published-source adapter tests

The CTest targets `stage2_reference_spectra`, `stage2_source_adapters`, and
`stage2_source_tool` exercise the production Stage-2 ingestion boundary
without Geant4 transport. `stage2_published_transport` then exercises the
end-to-end executable boundary. Together they cover:

- AMS primary-table rigidity-edge to kinetic-energy-per-nucleon conversion,
  the differential-flux Jacobian, area/energy unit conversion, independent
  proton and helium grids, and the prohibition against relabeling inclusive
  helium as pure alpha particles;
- the 2025 AMS-02 He--O solar-cycle tables: all 147 available rotations, all
  40 rigidity bins, separate statistical/time-dependent/time-independent
  errors, missing-rotation rejection, and no assumed isotope assignment;
- the optional preregistered AMS chronological development/holdout partition
  for a future source-environment comparison;
- ACE/CRIS B--Ni daily values, exact energy bounds from the pinned V04 CDF,
  counts, analytic `flux/sqrt(count)` statistical errors, and explicit
  zero-count masks;
- an enforced prohibition against combining the V04 CDAWeb values with the
  nonidentical bounds in the separately pinned current ACE factors table;
- strict hashes, versions, field sets, canonical round trips, deterministic
  output, and deliberate source-contract mutations that must fail;
- NOAA interval-start to midpoint/end conversion for every five-minute SGPS
  record and exact half-open overlap integration;
- separate `>100 MeV` primary and `>10 MeV` sensitivity windows, with the NMDB
  GLE72 onset retained only as an independent marker;
- both SGPS sensor directions and all 13 differential channels, without
  angular averaging;
- independent rejection of each record-count, valid-sample, quality-flag,
  missing-flux, and missing-uncertainty condition;
- closed valid/rejected/missing time ledgers and explicit incomplete status;
  and
- analytic constant/linear/partial-boundary cases plus an independently
  calculated published P7 event-fluence spot check;
- independent reconstruction of `sum(J_i Delta E_i)` from the downloaded AMS
  rigidity edges and proton values, plus the isotropic `pi` factor, source
  area, and incident-rate identity; and
- a deterministic 256-history multilayer Geant4 run producing TID and
  `H100/10`, with rejection tests for inclusive helium, beam relabeling, and
  mission-purpose use of the validation fixture.

`stage2_source_adapters` deliberately prints one ledger assertion per output
channel so a direction-, window-, or energy-specific failure is localized.
The tests do not claim that the near-Earth AMS spectrum is a complete
deep-space scenario or that directional GOES measurements are an isotropic
transport boundary. BON2020 is not used in production and is not required for
Stage 2.

## Stage-3 geometry, response, and convergence tests

`stage3_core_tests.cc` is the deterministic V1 boundary. It checks canonical
rigid transforms and inverse transforms, unit-quaternion rejection and `q/-q`
canonicalization, exact two-layer interface/areal-density accounting, physical
layer-order retention and hash invalidation, slab volume/mass/path formulas,
and response-matrix identity, exact zero, linearity, statistical/source-
covariance propagation, positive-semidefinite covariance enforcement, and
fatal no-extrapolation behavior.

`test_stage3_response_weights.py` builds an eight-corner synthetic tensor whose
response is affine in log-energy, `mu`, and periodic `phi`. It drives the
production assembler, weight builder, and folder, verifies the analytic
one-eighth corner weights and folded value, and requires exact-node,
extrapolation, and incomplete-grid requests to fail.

`test_stage3_transport.py` is the registered V3 boundary. CMake passes it a
dedicated `<build>/stage3_transport_runs` directory; the driver recreates only
that owned directory and writes no test case to `/tmp`. It archives every input
case, full Geant4 log, result, response cell, pooled five-seed cell, JSON/CSV
matrix, source-weight file, fold, and acceptance summary. It covers:

- direct-object translation/rotation and local rear-face scoring;
- independent absorber selection-order invariance;
- five independent seeds for single-vs-split Al and for direct-vs-folded
  W/HDPE and HDPE/W responses;
- per-layer step/deposit and ancestry-aware secondary-production ledgers;
- actual coarse/medium/fine response-grid folding through the production
  weight builder, with no extrapolation;
- three-level production-cut, maximum-step, common-transverse-size, and
  source-plane-size convergence, plus fourfold event-count scaling;
- proton, He-4, and Fe-56 construction; the `FTFP_BERT`, `FTFP_BERT_HP`,
  reference-ion, and `G4IonQMDPhysics` configurations; and representative Al,
  W, HDPE, Si, and water transport; and
- exact zero handling in a stopped-primary thick-W case.

The preregistered stochastic acceptance rules are a maximum 2% relative Monte
Carlo standard uncertainty for accepted nonzero comparisons, a two-sided 95%
difference interval wholly contained in ±5% for equivalence claims, and a
maximum 2% change at the final cut/step/size/grid refinement. Exact zero cases
use absolute equality. Those thresholds qualify only Stage-3 raw transport
infrastructure. Science quantities remain disabled and formal Stage-3
completion remains pending independent review.

## Stage-4 finite-slab TID tests

Stage 4 has three independent CTest gates. All generated inputs, logs, JSON,
and evidence are written below the CMake build tree; no Stage-4 test uses
`/tmp`.

`stage4_core_tests.cc` is deterministic and transport-free. It verifies exact
MeV/kg-to-Gy and Gy-to-rad conversion, the case-v2 TID
capability gate, DDD rejection when its separate Stage-5 response package is
absent, gamma `(0,0)` attenuation-node
admission, explicit rejection of gamma and heavy-ion TID, and the structured
uncertainty contract. It also rejects TID absorbers outside the registered
Al/Si/water material matrix without restricting shield-layer choices. A
complete budget must
contain the material/mass, source, numerical, and physics-reference categories
plus every unordered correlation pair. The test hand-calculates their
covariance, rejects an indefinite matrix and a legacy aggregate-only v2 claim,
and verifies the explicitly requested combined 95% interval.

`test_stage4_nist_reference.py` authenticates every entry in
`reference_tid/nist/manifest.sha256.json` before running. The raw official NIST
responses and their canonical CSVs are both pinned. The test invokes
`shieldSimReferenceTool` twice—`FTFP_BERT`'s reference EM constructor and
`G4EmStandardPhysics_option4`—and performs 63 exact-table-node comparisons:

- proton and helium-4 electronic and total stopping power in Al, Si, and water;
- CSDA mass range obtained by explicit log-energy quadrature of the same
  unrestricted Geant4 stopping model; and
- XCOM total photon mass attenuation with coherent scattering in Al, Si, and
  water.

Each comparison has its own registered relative limit in
`reference_tid/stage4_reference_plan.v1.json`; there is no universal tolerance
and interpolation uncertainty is exactly zero at the selected nodes. Photon
attenuation is a narrow-beam EM verification only. It is not converted into
finite-slab dose, and gamma TID is not validated. ESTAR/electrons and ions
heavier than helium remain outside Stage 4.

`test_stage4_tid_transport.py` is the full stochastic proxy-path gate. Its
registered plan is `reference_tid/stage4_transport_plan.v1.json`; artifacts go
to `<build>/stage4_tid_transport_runs`. It runs proton and helium-4 incidence
through the normal case parser, primary generator, multilayer stack, Geant4
transport, finite absorber, adaptive accuracy controller, and canonical result
writer for Al, Si, and water. For every accepted configuration it requires:

- exact analytic slab mass and `ionizing deposited energy / mass` dose;
- the exact rule that the numeric rad value is 100 times the numeric Gy value;
- independent closure of total/ionizing/nonionizing ledgers when partitioned
  by transported species, propagated ancestry, creator process, or step-
  defining process;
- process-name partitions explicitly labeled diagnostic rather than unique
  causal attribution;
- no more than 2% Monte Carlo relative standard uncertainty; and
- a complete four-category uncertainty record with the Stage-4 pending-review
  capability label.

Three levels each of production cut, maximum step, and transverse size must
change the final accepted observable by no more than 2%. Cut and step are
tested on TID at fixed geometry. Changing transverse slab area changes mass
and therefore Gy per one compact-plane primary, so lateral-escape convergence
is correctly tested on ionizing deposited energy; exact dose conversion is
still checked for every size. The test also limits option-4 TID sensitivity to
the registered 4% allocation and verifies `1/sqrt(N)` behavior with five
independent low/high-history seed pairs. Low-history points are convergence
diagnostics, not members of the ≤2% final acceptance population.

Together these gates qualify only proton/helium-4 finite material-slab TID in
Al, Si, and water. Formal Stage-4 completion still requires the independent-
review record. Passing Stage 4 does not enable the separately gated Stage-5
DDD/`n_eq` capability, LET, `H100/10`, device response, organ dose, or crew-
health claims.

## Stage-5 silicon DDD and neutron-equivalent tests

Stage 5 has two CTest gates and a downstream figure gate. Its response files
and registered campaign are documented in `reference_ddd/README.md`; generated
artifacts remain below the build tree and no Stage-5 test writes to `/tmp`.

`stage5_core_tests.cc` is deterministic and transport-free. It loads the same
strict `DamageResponseLibrary` used in production, authenticates each table,
and verifies:

- exact 100 MeV proton and helium-4 SR-NIEL nodes;
- the geometric-mean identity for log-log interpolation and a closed-form
  power-law interval fold over explicit source-bin edges;
- the exact `[1.0, 1.1)` MeV ASTM neutron-group fold and the independently
  converted 95 MeV-mb silicon reference;
- inclusive 20 MeV neutron support followed by immediate out-of-domain
  rejection, partial-bin rejection, and harmless zero fluence outside domain;
- explicit unsupported-electron response; and
- case-v2 rejection of a missing response package, a non-silicon absorber, or
  an unsupported incident species.

`test_stage5_transport.py` drives the full executable path. The preregistered
exact-node cases use a 0.001 mm upstream aluminum layer and require the
absorber-entry result to agree with the independently computed table value to
0.5%; crossing count must be at least 0.995 per primary. It checks proton and
helium-4 DDD, neutron DDD, and neutron-only `n_eq`, and closes the result to the
event-level raw score and response-cell primitive. It also verifies that the
response cache key contains the damage-manifest hash and that the Geant4
nonionizing-deposition diagnostic is not added to the authoritative fold.

A transported 25 MeV neutron must produce null DDD/`n_eq` values with
`not_assessable` accuracy and nonzero out-of-domain counters; it may not be
evaluated at 20 MeV. A separate 1 GeV proton/10 cm aluminum sensitivity case
must produce absorber-entry secondary neutrons. Its expected high-energy
neutrons exercise the same strict-domain policy while leaving supported partial
scores auditable. These are implementation and reference checks, not an
accelerator/device validation.

An additional 100 MeV proton case at `mu = 0.5` checks the analytic
`1/(A*mu)` scalar-fluence normalization while compact source placement prevents
finite-face losses from contaminating that identity.

Formal Stage-5 completion remains pending a project-selected accelerator or
device-degradation data set and independent review. No test assigns a
fabricated SR-NIEL uncertainty: a combined accuracy claim requires the case to
supply a provenance-bearing `physics_reference` component.

## Stage-6 LET and H100/10 tests

Stage 6 adds four CTest gates plus the downstream figure gate. Generated cases,
logs, result JSON, Tecplot spectra, reference-tool output, and evidence stay
under the configured build tree (`stage6_*_runs/` and
`validation_evidence/`); these tests do not use `/tmp`.

`stage6_core_tests.cc` is Geant4-free. It checks the exact 140-bin recommended
LET grid (0.001--10000 keV/um at 20 bins/decade), inclusive final edge,
unclamped under/overflow, exact constant and power-law partial-bin threshold
integration, a closed-form hardness fold, correlated ratio-of-means, and both
zero and Fieller-unresolved denominators.

`test_stage6_case_boundaries.py` drives the production parser. It rejects a
missing or unused `let_response`, wrong LET units, nonintegral log-bin spans,
gamma incidence, C-12 outside natural silicon, and C-12 with the reference-list
EM default. C-12 in silicon is accepted only with option 4 because that is the
configuration that passes the fixed ATIMA screen.

`test_stage6_reference.py` authenticates all reused NIST files and the pinned
GSI ATIMA responses under `reference_let/`. It compares unrestricted electronic
stopping at all 18 registered p/He-4 NIST nodes in Al/Si/water. Fully stripped
C-12 in silicon is checked at 10, 100, and 1000 MeV/u with a fixed 15% option-4
implementation screen. The screen is not an uncertainty assignment or an
experimental-validation claim; the rejected default-model 10 MeV/u result is
documented in the roadmap as the reason for the parser boundary.

`test_stage6_transport.py` exercises the full Geant4 path. Thin-silicon proton,
He-4, and C-12 primary component track means are compared with deterministic
midpoint stopping queries to 3%. The aggregate spectrum still contains delta
electrons as separately identified components. Every run closes binned plus
under/overflow track fluence to `sum(dx/V)`. Separate max-step, production-cut,
and absorber-thickness perturbations use the approved 10% Stage-6 convergence
limit. A direct 50/150 MeV proton mixture checks analytic `H100/10=2/3`,
including Fieller coverage, while an alpha-only run requires null,
`zero_denominator`, and `not_assessable`.

Stage-6 numerical evidence files are `stage6_core.json`,
`stage6_reference.json`, and `stage6_transport.json`. Formal Stage-6 completion
still requires the independent-review record; CRaTER remains a separate future
instrument-response validation and is not used as generic LET truth.

## Stage-7 published-component and time-integration tests

Stage 7 adds three CTest gates and participates in the downstream figure gate.
All cases, logs, response cells, selected-score matrices, weights, folds,
exposure ledgers, and evidence stay below the configured build tree. The test
temporary directory is also inside that tree; no Stage-7 test uses `/tmp`.

`stage7_core_tests.cc` is deterministic and Geant4-free. It verifies exact
log-energy interpolation, refusal to extrapolate, an analytic spectrum-response
integral with adaptive refinement, exact pilot-based Neyman allocation, strict
UTC/leap-day parsing, piecewise integration across partial coverage and an
explicit gap, full interval covariance propagation, and a sampled-source
ensemble. It also rejects an implicit time gap, overlap, indefinite covariance,
a false full-GCR label on the AMS proton component, and GOES transport without
the missing boundary mapping.

`test_stage7_boundaries.py` drives all three `shieldSimExposureTool` commands.
It checks deterministic output, strict unknown-field and timestamp handling,
exact gap-preserving integration, no-extrapolation quadrature, and an exact
history allocation with no physics biasing. The committed example has an
incomplete physical uncertainty block and must therefore remain
`not_assessable`.

`test_stage7_published_qualification.py` authenticates the registered Stage-7
plan and production source adapters. It retains both GOES directions and
confirms their transport-disabled status. All 45 downloaded AMS A=1 proton bins
then drive five-seed production runs for Al→HDPE and HDPE→Al stacks with
independent 1 mm natural-silicon and water slabs. Only Stage-4-qualified TID is
requested. For each silicon stack, an explicit energy/angle quadrature is run
both directly and through monoenergetic response matrices. Selected-score
assembly retains the TID primitive and primary-proton crossing monitor after
full validation of every response cell. Independent 256-history pilots drive a
frozen 73,728-history 12-node base allocation and a 56,320-history interleaved
energy-refinement allocation for each of five response replicates. The
23-node reversed-stack result passed the 10% difference gate but returned
`z=2.078`, above the unchanged `1.96` gate. The registered corrective
refinement therefore adds the remaining 22 published representatives with a
separate 112,640-history allocation per replicate. Every allocation has a
512-history-per-cell floor; the complete 45-node response campaign contains
1,213,440 production histories per stack. The
response uncertainty is the conservative
larger of propagated within-replicate and empirical between-replicate standard
error. This explicitly tests the rare TID tail instead of trusting one random
stream's variance estimate. Pilot moments are excluded and production seeds
are disjoint.

The same maximum-within/between rule is applied to the five equal-history
direct seeds. This is required symmetrically because the 45-node campaign
removed energy interpolation but the reversed-stack direct seed means still
ranged from 0.598 to 0.731 MeV/primary. Their empirical between-seed standard
error was 20.8% larger than the pooled event estimate. The direct mean, 10%
relative gate, and `z <= 1.96` gate are unchanged; the test also requires the
resulting direct relative standard error to remain at or below 10%.

Every full response cell is parsed by `shieldSimResponseTool` before retention.
The driver then keeps only the registered TID and primary-crossing scores and
their complete 2x2 covariance block, plus a content-hashed completion marker.
The redundant full result is removed for response-only cases; published and
direct-quadrature results remain unabridged. This compaction prevents the full
diagnostic covariance triangles from exhausting the test filesystem without
weakening cell validation or matrix inputs.

The driver can execute independent case directories concurrently with
`SHIELDSIM_STAGE7_WORKERS=N`. Each worker is a separate serial `shieldSim`
process; seeds, histories, case documents, output prefixes, result ordering,
and acceptance calculations do not depend on completion order. The default is
one worker, invalid/nonpositive values are rejected, and the requested count
may not exceed the CPUs visible to the process. This is task parallelism, not
MPI or Geant4 event multithreading. Exact-identity restarts reuse only complete
canonical cases, so an interrupted parallel campaign can be resumed safely.

The registered gates require five seeds, no more than 2% change in the final
response-grid monitor, and direct-versus-matrix TID agreement both within 10%
and within the combined two-sided 95% statistical interval. The normal
published sampler is also compared with the explicit quadrature. Each matrix
row is converted to an exact 27-day exposure. Because the published source
does not provide the covariance needed by this fold, uncertainty remains
incomplete and accuracy remains `not_assessable`; the test never creates a
diagonal covariance assumption.

Every direct row uses 8,000 histories for each of its five seeds. This
40,000-history total was fixed after a smaller development pilot exposed an
unstable heavy TID tail; no acceptance threshold was relaxed.

Passing these gates implements and qualifies the Stage-7 software framework
and the downloaded proton component. It does not complete physical GCR/SEP
qualification. Full completion still needs an approved isotope-resolved GCR
scenario, the GOES sensor-to-simulation transform/angular response, complete
physical uncertainty budgets, and the required earlier-stage reviews.

## Validation figures

CTest target `validation_figures` depends on all current analytic and
observational test targets. After their numerical assertions pass, it reads
their compact JSON evidence and the pinned observational fixtures, reruns the
production source tool, and writes review artifacts to:

```text
<build>/validation_evidence/*.json
<build>/validation_figures/*.png
<build>/validation_figures/*.svg
<build>/validation_figures/manifest.json
```

The twenty-six figure subjects are Stage-1 core formulas, Stage-1 case boundaries,
Stage-1 analytic slab masses/ledger closure, original AMS primary-table
analytic checks, the end-to-end AMS source-integral check, original AMS
observation/adapter agreement, AMS-02 2025
solar-cycle observation/adapter agreement, ACE/CRIS spectra and their
observation/statistical-error spot checks, GOES analytic/oracle integration,
GOES event observations/coverage, the CRaTER observational reference, Stage-3
geometry/response algebra, Stage-3 response-grid interpolation algebra,
Stage-3 transport/convergence comparisons, Stage-4 exact dose/covariance
algebra, all Stage-4 NIST comparisons, Stage-4 TID
transport/convergence/ledger comparisons, Stage-5 analytic/table folds, and
Stage-5 full-transport absorber-entry comparisons, Stage-6 exact LET/hardness
mathematics, Stage-6 NIST/ATIMA stopping references, Stage-6 full-transport
LET/hardness comparisons, Stage-7 exact time/covariance/quadrature algebra,
Stage-7 executable-boundary algebra, and Stage-7 downloaded proton-component
direct/matrix/grid comparisons. The
manifest records captions, byte sizes, and SHA-256 hashes.

The standalone Geant4 layer runners extend that coverage. They invoke
`generate_layer_validation_figures.py` only after their numerical assertions
pass and write to `${FIGURE_DIR}` (default:
`${RUN_DIR}/validation_figures`):

- Layer 2: exact pencil-beam source contract, isotropic plane-crossing
  `p(mu)=2*mu` distribution/moments, and near-vacuum rear-face limits.
- Layer 3: tenfold source-normalization scaling and the H100/10 below/above
  100 MeV threshold limits.
- Layer 3 also renders the existing controlled proton/alpha LET sanity case,
  but its caption and manifest do not represent it as an analytic or
  observational validation.

Every subject is saved in both PNG and SVG form. Each output directory has its
own `manifest.json` with captions, byte sizes, SHA-256 hashes, and
`physical_validation_claim: false`.

Figure generation requires Python `matplotlib` and uses its noninteractive
`Agg` backend. Figures are deliberately downstream of the numerical checks:
no test compares pixels, a plotting change cannot make an incorrect numeric
result pass, and observation-only figures contain explicit warnings against a
physical-validation interpretation.

## Test Philosophy

The purpose of the test suite is to separate several different questions:

1. Does the code build?
2. Does the command-line interface work?
3. Are invalid inputs rejected cleanly?
4. Are the geometry, source, and scoring algorithms behaving correctly?
5. Do the physics-facing outputs respond correctly to simple controlled cases?
6. Are Monte Carlo results statistically stable enough for the intended use?
7. Do future code changes preserve validated behavior?

Layer-1 tests should be quick and should be run frequently. Layer-2 tests are
also intended to be lightweight, but they run short Geant4 simulations and use
diagnostic output to check geometry/source/scoring behavior. Layer-3 tests run
more physics-facing simulations and are intended as sanity checks, not as a full
validation against external reference data. Layer-4 tests implement a baseline-driven regression workflow: generate an accepted baseline from a trusted version, then compare future runs against it with documented tolerances.

---

## Directory Layout

```text
tests/
├── README.md
├── stage1_core_tests.cc
├── stage2_source_adapter_tests.cc
├── stage3_core_tests.cc
├── stage8_crater_core_tests.cc
├── test_stage3_response_weights.py
├── test_stage3_transport.py
├── generate_validation_figures.py
├── generate_layer_validation_figures.py
├── test_stage1_case_boundaries.py
├── test_stage2_source_tool.py
├── test_stage1_transport_smoke.py
├── test_reference_spectra.py
├── test_crater_reference_inputs.py
├── test_stage8_crater_transport.py
├── reference_crater/
│   ├── README.md
│   ├── acquisition.v1.json
│   ├── manifest.sha256
│   ├── glace/
│   ├── pds/
│   ├── publications/
│   └── zenodo/
├── reference_spectra/
│   ├── README.md
│   ├── manifest.sha256
│   ├── gcr/
│   ├── sep/
│   └── observations/
├── run_layer1_tests.sh
├── run_layer2_tests.sh
├── run_layer3_tests.sh
├── run_layer4_tests.sh
├── regression_tools.py
├── data/
│   ├── mono_spectrum_metadata.v1.json
│   ├── mono_50MeV_proton.dat
│   ├── mono_100MeV_proton.dat
│   ├── mono_100MeV_proton_rate10.dat
│   ├── mono_100MeV_alpha.dat
│   ├── mono_100MeV_alpha_rate10.dat
│   └── mono_150MeV_proton.dat
├── expected/
│   ├── README.md
│   └── layer4/
│       ├── README.md
│       └── layer4_reference.json    # generated after accepting a baseline
└── logs/
    └── .gitignore
```

### `run_layer1_tests.sh`

Build and smoke tests.

This script verifies that the code builds, that basic CLI commands work, and
that invalid input is rejected cleanly.

### `run_layer2_tests.sh`

Geometry, source, and scoring tests.

This script verifies beam-source sampling, isotropic-source sampling,
near-vacuum transmission through the shield rear face, target ordering, and
basic material alias behavior.

### `run_layer3_tests.sh`

Physics-output and numerical-sanity tests.

This script verifies that computed-quantity outputs are written, that source
normalization scales rates correctly, that H100/10 behaves correctly for simple
below/above-threshold spectra, that alpha LET is larger than proton LET in a
controlled silicon case, that all supported physics lists run, and that sweep
output has monotonic areal density.


### `run_layer4_tests.sh`

End-to-end regression tests.

This script runs deterministic reference cases, collects compact signatures from
run-summary and output files, and compares the current results against an
accepted JSON baseline. It can also create or update that baseline after manual
review.

### `regression_tools.py`

Python helper used by `run_layer4_tests.sh` to collect run summaries into JSON
and compare actual results against the accepted baseline with configurable
numerical tolerances.

### `data/`

Small deterministic input spectra for tests. These files are intentionally
simple and should be version controlled. Rate-scaled spectra have the same shape
as their unit-rate counterpart, but a larger absolute source normalization. They
are used to test that per-primary Monte Carlo quantities are unchanged while
physical rates scale with the input spectrum intensity.

Every `mono_*.dat` use in the compatibility CLI passes
`mono_spectrum_metadata.v1.json`. A three-column file is intentionally
insufficient on its own: it does not define energy coordinate, units, angular
semantics, bin edges, or whether values are physical rates or test weights.
The shared sidecar labels these fixtures `synthetic_test_only` and declares the
one-bin width used by the sampler.

### `expected/`

Reference output summaries and accepted regression baselines. The Layer-4 baseline is stored under `tests/expected/layer4/layer4_reference.json` after it is generated with `tests/run_layer4_tests.sh --update-baseline`. Large raw
Geant4 output files should generally not be committed unless they are small and
essential.

### `logs/`

Generated test logs. This directory should not be version controlled except for
its `.gitignore` file.

---

## Diagnostic and Numerical CLI Options Used by Tests

The test suite uses several CLI options implemented in the main code:

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

### `--random-seed=<integer>`

Sets the CLHEP/Geant4 random seed before the run. This makes source-sampling and
short Monte Carlo tests reproducible.

### `--output-prefix=<name>`

Changes standard output file names. With the default prefix `shieldSim`, files
keep their historical names, such as:

```text
shieldSim_spectra.dat
shieldSim_quantities.dat
shieldSim_let_spectrum.dat
shieldSim_dose_sweep.dat
shieldSimOutput*.csv
```

With `--output-prefix=caseA`, they become:

```text
caseA_spectra.dat
caseA_quantities.dat
caseA_let_spectrum.dat
caseA_dose_sweep.dat
caseAOutput*.csv
```

### `--dump-source-samples=<file>`

Writes one row per generated primary:

```text
row species E_MeV x_mm y_mm z_mm ux uy uz
```

Positions are global coordinates in millimeters. The vector `(ux,uy,uz)` is the
unit momentum direction requested by `PrimaryGeneratorAction`.

This file is used to verify that:

- beam mode generates `x = y = 0` and direction `+z`,
- isotropic mode keeps positions on the upstream plane,
- isotropic mode samples the inward-hemisphere crossing distribution
  `p(mu) = 2 mu`, where `mu = uz`.

### `--dump-exit-particles=<file>`

Writes one row for each particle accepted as crossing the downstream shield
face:

```text
row species E_MeV xg_mm yg_mm zg_mm xl_mm yl_mm zl_mm uxl uyl uzl
```

The `g` coordinates are global coordinates. The `l` coordinates and directions
are in the local coordinate frame of the shield. This file is used to verify
that shield rear-face scoring is based on the correct local-coordinate crossing
condition.

### `--dump-run-summary=<file>`

Writes a machine-readable scalar summary for compatibility Layer-3 tests. It
gives stable keyword rows that are easier to parse in automated tests:

```text
scalar H100_10 <value>
count output_proton <value>
count H100_10_primary_proton_numerator <value>
count H100_10_primary_proton_denominator <value>
target index name thickness_mm TID_Gy_perPrimary TIDRate_Gy_s DDD_MeV_g_perPrimary DDDRate_MeV_g_s n_eq_cm2_perPrimary n_eq_rate_cm2_s
deposit_ledger index total_Gy_perPrimary ionizing_Gy_perPrimary nonionizing_Gy_perPrimary
```

In sweep mode, the file contains one `begin_run` / `end_run` block per sweep
thickness.

### `--production-cut=<mm>`

Sets the Geant4 default production range cut in millimeters. The Layer-3 script
uses this as a smoke/numerical-control test. Production runs should document the
chosen value because low-energy secondaries, TID, DDD, `n_eq`, and LET tails can
be cut-sensitive.

### `--max-step=<mm>`

Applies a maximum step length to the shield and scoring slabs and registers
`G4StepLimiterPhysics`. This is mainly for numerical convergence checks in thin
targets and LET-spectrum tests. Omit it to use the normal physics-list stepping
behavior.

### `--diagnostic-max-rows=<n>`

Limits the number of rows written to each diagnostic file. This prevents
accidentally generating huge diagnostic dumps. Use a value less than or equal to
zero to remove the explicit cap.

---

## Layer-1 Tests

Layer-1 tests are software-level build and CLI smoke tests.

They check:

- clean CMake configuration,
- clean build,
- executable discovery,
- `--help`,
- `--list-materials`,
- `--list-target-materials`,
- `--list-detector-materials`,
- `--list-quantities`,
- rejection of invalid shielding material,
- rejection of invalid physics list,
- rejection of invalid source mode,
- rejection of negative target thickness,
- rejection of invalid computed quantity,
- rejection of malformed shield specification,
- rejection of negative event count.

Run from the top-level package directory:

```bash
chmod +x tests/run_layer1_tests.sh
tests/run_layer1_tests.sh
```

To see script options:

```bash
tests/run_layer1_tests.sh --help
```

Optional environment variables:

```bash
BUILD_DIR=build_test JOBS=8 tests/run_layer1_tests.sh
```

Available variables:

```text
BUILD_DIR   Build directory used by the test script.
            Default: build_layer1_tests

EXE_NAME    Executable name.
            Default: shieldSim

LOG_DIR     Directory where test logs are written.
            Default: layer1_test_logs

JOBS        Number of parallel build jobs.
            Default: nproc, or 2 if nproc is unavailable.
```

Expected result:

```text
All Layer-1 tests passed.
```

---

## Layer-2 Tests

Layer-2 tests verify geometry, source sampling, and scoring logic using short
Geant4 runs and diagnostic files.

The implemented Layer-2 script checks:

1. The diagnostic options are present in `--help`.
2. Beam mode generates 100 MeV protons at `x = y = 0` with direction `+z`.
3. Isotropic mode samples the source plane and direction distribution correctly:
   - all directions point inward, `0 < mu <= 1`,
   - `<mu> ≈ 2/3`,
   - `<mu^2> ≈ 1/2`,
   - source positions stay within the finite upstream plane.
4. A near-vacuum shield transmits most 100 MeV protons through the downstream
   shield face.
5. Exit diagnostics use shield-local coordinates and positive downstream local
   direction cosine.
6. Detector/target material order is preserved in output metadata.
7. `Al` and `G4_Al` both run successfully as shield material specifications.

Run from the top-level package directory:

```bash
chmod +x tests/run_layer2_tests.sh
tests/run_layer2_tests.sh
```

To see script options:

```bash
tests/run_layer2_tests.sh --help
```

Optional environment variables:

```bash
BUILD_DIR=build_l2 JOBS=8 tests/run_layer2_tests.sh
```

Available variables:

```text
BUILD_DIR   Build directory used by the test script.
            Default: build_layer2_tests

EXE_NAME    Executable name.
            Default: shieldSim

LOG_DIR     Directory where test logs are written.
            Default: layer2_test_logs

RUN_DIR     Directory where each Layer-2 test case writes its run outputs.
            Default: layer2_test_runs

FIGURE_DIR  Directory for PNG/SVG figures and their manifest.
            Default: RUN_DIR/validation_figures

JOBS        Number of parallel build jobs.
            Default: nproc, or 2 if nproc is unavailable.

SKIP_BUILD  Set to 1 to skip CMake configure/build.
            Default: 0

EXE_PATH    Explicit executable path when using SKIP_BUILD=1.
```

Example using an already-built executable:

```bash
SKIP_BUILD=1 EXE_PATH=build/shieldSim tests/run_layer2_tests.sh
```

Expected result:

```text
All Layer-2 tests passed.
```

Layer-2 test outputs are written under:

```text
layer2_test_runs/
```

The three analytic figure subjects, each in PNG and SVG form, and their
manifest are written under:

```text
layer2_test_runs/validation_figures/
```

Logs are written under:

```text
layer2_test_logs/
```

These generated directories should not be committed.

---

## Layer-3 Tests

Layer-3 tests verify physics-output behavior and numerical-control plumbing with
short, deterministic Geant4 runs. They are sanity tests, not a substitute for
full validation against NIST stopping powers, SR-NIEL/device-response tables, or
experimental data.

The implemented Layer-3 script checks:

1. `--help` lists `--dump-run-summary`, `--production-cut`, and `--max-step`.
2. A monoenergetic 100 MeV proton run writes spectra, scalar quantity, LET, and
   summary outputs.
3. Source normalization is applied correctly: scaling the input spectrum by 10
   leaves per-primary TID unchanged but scales the source-normalized TID rate by
   10.
4. H100/10 behaves correctly for simple monoenergetic spectra:
   - 50 MeV protons give `H100/10 ≈ 0`,
   - 150 MeV protons give `H100/10 ≈ 1`.
5. The folded mean LET for 100 MeV total-energy alpha particles in silicon is
   larger than the folded mean LET for 100 MeV protons.
6. All supported physics lists run in a short smoke test:
   - `FTFP_BERT`,
   - `FTFP_BERT_HP`,
   - `Shielding`,
   - `QGSP_BIC_HP`.
7. A three-point Al sweep writes scalar outputs and increasing areal density.

Run from the top-level package directory:

```bash
chmod +x tests/run_layer3_tests.sh
tests/run_layer3_tests.sh
```

To see script options:

```bash
tests/run_layer3_tests.sh --help
```

Optional environment variables:

```bash
BUILD_DIR=build_l3 JOBS=8 tests/run_layer3_tests.sh
EVENTS_SHORT=10000 EVENTS_SMOKE=1000 tests/run_layer3_tests.sh
```

Available variables:

```text
BUILD_DIR       Build directory used by the test script.
                Default: build_layer3_tests

EXE_NAME        Executable name.
                Default: shieldSim

LOG_DIR         Directory where test logs are written.
                Default: layer3_test_logs

RUN_DIR         Directory where each Layer-3 test case writes its run outputs.
                Default: layer3_test_runs

FIGURE_DIR      Directory for PNG/SVG figures and their manifest.
                Default: RUN_DIR/validation_figures

JOBS            Number of parallel build jobs.
                Default: nproc, or 2 if nproc is unavailable.

SKIP_BUILD      Set to 1 to skip CMake configure/build.
                Default: 0

EXE_PATH        Explicit executable path when using SKIP_BUILD=1.

EVENTS_SHORT    Events for ordinary Layer-3 cases.
                Default: 2000

EVENTS_SMOKE    Events for quick physics-list/sweep smoke cases.
                Default: 300
```

Example using an already-built executable:

```bash
SKIP_BUILD=1 EXE_PATH=build/shieldSim tests/run_layer3_tests.sh
```

Expected result:

```text
All Layer-3 tests passed.
```

Layer-3 test outputs are written under:

```text
layer3_test_runs/
```

The two analytic figure subjects plus the explicitly non-validation controlled
LET figure, each in PNG and SVG form, and their manifest are written under:

```text
layer3_test_runs/validation_figures/
```

Logs are written under:

```text
layer3_test_logs/
```

These generated directories should not be committed.

---

## Layer-4 Regression Tests

Layer-4 tests verify regression stability of complete end-to-end cases. They are
intended to answer: did a code change alter previously accepted results beyond
the configured tolerance? They do **not** independently validate Geant4 physics.

The implemented Layer-4 workflow uses three modes:

```bash
tests/run_layer4_tests.sh --check            # default: compare against baseline
tests/run_layer4_tests.sh --update-baseline  # replace accepted baseline
tests/run_layer4_tests.sh --actual-only      # run and collect actual JSON only
```

The first trusted Layer-4 run should be reviewed manually, then accepted with:

```bash
chmod +x tests/run_layer4_tests.sh
tests/run_layer4_tests.sh --update-baseline
```

After a baseline exists, normal regression checks are run with:

```bash
tests/run_layer4_tests.sh
```

The accepted baseline is written to:

```text
tests/expected/layer4/layer4_reference.json
```

The current-run result is written by default to:

```text
layer4_test_runs/layer4_actual.json
```

The implemented Layer-4 script runs these deterministic cases:

```text
regression_001_beam_100MeV_p_Al2mm_Si1mm
regression_002_beam_100MeV_alpha_Al2mm_Si1mm
regression_003_isotropic_SEP_Al2mm_BFO50mm_Si1mm
regression_004_builtin_HDPE10mm_BFO50mm_Si1mm
regression_005_highZ_Shielding_W5mm_Si1mm
regression_006_sweep_Al_0p5_to_2mm_Si1mm
```

For each case, the script compares compact signatures derived from:

```text
run-summary metadata
input/output particle counts
source normalization
incident particle rate
shield thickness and areal density
H100/10
target TID, TID rate, DDD, DDD rate, n_eq, and n_eq rate
LET-spectrum integral and mean LET
compact numeric signatures of spectra and quantities output files
```

The comparison is performed by:

```text
tests/regression_tools.py
```

Default tolerances are intentionally moderate because these are short Monte
Carlo regression cases. The defaults can be overridden:

```bash
REL_TOL=0.05 ABS_TOL=1e-30 tests/run_layer4_tests.sh --check
EVENTS_REGRESSION=20000 EVENTS_SMOKE=5000 tests/run_layer4_tests.sh --check
```

Available variables:

```text
BUILD_DIR          Build directory used by the test script.
                   Default: build_layer4_tests

EXE_NAME           Executable name.
                   Default: shieldSim

LOG_DIR            Directory where test logs are written.
                   Default: layer4_test_logs

RUN_DIR            Directory where each Layer-4 test case writes run outputs.
                   Default: layer4_test_runs

JOBS               Number of parallel build jobs.
                   Default: nproc, or 2 if nproc is unavailable.

SKIP_BUILD         Set to 1 to skip CMake configure/build.
                   Default: 0

EXE_PATH           Explicit executable path when using SKIP_BUILD=1.

EVENTS_REGRESSION  Events for ordinary regression cases.
                   Default: 4000

EVENTS_SMOKE       Events for high-Z and sweep smoke cases.
                   Default: 1000

REL_TOL            Relative tolerance for numerical comparison.
                   Default: 0.15

ABS_TOL            Absolute tolerance for numerical comparison.
                   Default: 1e-30

BASELINE_FILE      Accepted baseline JSON.
                   Default: tests/expected/layer4/layer4_reference.json

ACTUAL_FILE        Actual JSON produced by the current run.
                   Default: layer4_test_runs/layer4_actual.json
```

Example using an already-built executable:

```bash
SKIP_BUILD=1 EXE_PATH=build/shieldSim tests/run_layer4_tests.sh --check
```

Layer-4 generated outputs are written under:

```text
layer4_test_runs/
```

Logs are written under:

```text
layer4_test_logs/
```

Do not update the Layer-4 baseline just to make tests pass. Update it only when
the changed behavior is expected and has been reviewed.

---

## What Should Be Committed

Commit:

```text
tests/*.sh
tests/README.md
tests/data/*.dat
tests/expected/*.md
tests/expected/*.json
tests/expected/layer4/*.md
tests/expected/layer4/*.json
tests/expected/*.txt
tests/logs/.gitignore
```

Do not commit generated build directories, logs, test-run directories, or large
raw output files.

Recommended `.gitignore` entries:

```text
build_layer1_tests/
build_layer2_tests/
build_layer3_tests/
build_layer4_tests/
layer1_test_logs/
layer2_test_logs/
layer3_test_logs/
layer4_test_logs/
layer2_test_runs/
layer3_test_runs/
layer4_test_runs/
tests/logs/
```

---

## Notes on Monte Carlo Testing

Monte Carlo tests are not exactly reproducible unless the random seed, Geant4
version, physics list, production cuts, and geometry are fixed. For reliable
regression testing, use:

```bash
--random-seed=<integer>
--output-prefix=<name>
--dump-run-summary=<file>
```

Useful future diagnostic or convenience options include:

```bash
--mono-proton=<MeV>
--mono-alpha=<MeV>
--disable-protons
--disable-alphas
```

The current package already implements source, exit, and run-summary diagnostic
dumps, which are enough for the implemented Layer-2, Layer-3, and Layer-4 tests.
