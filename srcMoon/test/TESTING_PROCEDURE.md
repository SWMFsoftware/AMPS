# `srcMoon` comprehensive testing procedure

## 1. Purpose and roadmap consistency

This procedure implements the layered campaign in
[`development_and_validation.md`](../development_and_validation.md):

```text
local production-kernel correctness (U01-U26)
→ linked production wiring (I01-I10, M0)
→ numerical convergence and controlled sphere (I11-I13, M1)
→ LOLA/illumination/thermal/surface physics (I14-I20, M2)
→ chemistry (I21-I24, M3)
→ drivers and sources (I25-I31, M4)
→ independent observations (I32-I40, M5)
→ cross-species freeze (I41, M6)
→ UQ/reproducibility/release (I42-I46, M7)
```

The requested architecture—one global runner, one runner and reference
contract per test, and visualizations of evaluated comparisons—is consistent
with the roadmap. It is also necessary for I42. It does **not** permit later
tests to substitute for failed earlier layers. A downloaded dataset, a
successful build, a source-code guard, or a demonstration plot is not a linked
physics PASS.

The authoritative global entry point is:

```bash
python3 -B srcMoon/test/test_runner.py --phase <scope>
```

The global runner delegates to the existing scientific runners; it contains no
lunar physics and does not recompute their metrics. `run_all_test.py` remains
the U/M1/M2 phase engine, and `run_integration_tests.py` remains the
configuration-separated M0 engine.

## 2. Status and milestone interpretation

Two questions must always be reported separately:

1. **How far does implementation reach?** Production code presently reaches
   into M2, with isolated experimental work beyond it.
2. **What is the latest qualified milestone?** M0 is complete. M1 remains the
   active frontier until I11-I13 all pass. Therefore M2 cannot be declared
   qualified even when individual M2 kernels or demonstrations run.

The maintained answer and change history are in [`../STATUS.md`](../STATUS.md).
Every work package that changes production capability, tests, data readiness,
thresholds, inputs, or campaign evidence must update that file.

## 3. Runner scopes

First audit the test architecture without compiling or running AMPS:

```bash
python3 -B srcMoon/test/test_runner.py --phase audit \
  --output-dir test_output/srcMoon/audit-YYYYMMDD
```

The orchestrator itself has a standard-library-only infrastructure self-test:

```bash
PYTHONDONTWRITEBYTECODE=1 python3 -B \
  srcMoon/test/test_runner_self_test.py
```

This checks exact implemented/planned registries, all per-ID contract files,
and fail-closed SKIPPED semantics.  It is not a U/I physics test and cannot
open a milestone gate.

List all frozen IDs. Implemented and planned IDs are labelled separately:

```bash
python3 -B srcMoon/test/test_runner.py --list
```

Run only the local production-kernel layer:

```bash
python3 -B srcMoon/test/test_runner.py --phase local \
  --output-dir test_output/srcMoon/local-YYYYMMDD
```

Run the configuration-separated M0 baseline. This performs the required clean
SPICE, no-SPICE, Na/Na+ and canonical regression builds:

```bash
python3 -B srcMoon/test/test_runner.py --phase m0 \
  --output-dir test_output/srcMoon/m0-YYYYMMDD
```

Run the current numerical/environment frontier only. This is useful for a
diagnostic rerun but does not prove that M0 passed in the same invocation:

```bash
python3 -B srcMoon/test/test_runner.py --phase frontier \
  --campaign-mode regression --mpi-ranks 8 \
  --output-dir test_output/srcMoon/frontier-YYYYMMDD
```

Run every implemented layer in roadmap order. The runner stops I11-I20 when
the M0 aggregate gate does not pass:

```bash
python3 -B srcMoon/test/test_runner.py --phase implemented \
  --campaign-mode regression --mpi-ranks 8 \
  --output-dir test_output/srcMoon/implemented-YYYYMMDD
```

Select the bounded M1 qualification campaign explicitly:

```bash
python3 -B srcMoon/test/test_runner.py --phase implemented \
  --campaign-mode qualification --mpi-ranks 8 \
  --output-dir test_output/srcMoon/qualification-YYYYMMDD
```

The release view adds I21-I46 as explicit `SKIPPED / not-run` records until
their individual production-linked runners exist. Release scope always rejects
skips:

```bash
python3 -B srcMoon/test/test_runner.py --phase release \
  --campaign-mode qualification --mpi-ranks 8 \
  --output-dir test_output/srcMoon/release-YYYYMMDD
```

The normal build sequence remains:

```bash
rm -rf build
./Config.pl -application=moon \
  -spice-path=/home/vtenishe/SPICE/cspice \
  -spice-kernels=/home/vtenishe/SPICE/Kernels
make -j "OPT3=-O3 -g"
```

Generated `Makefile` changes are local build products and are not committed.
Use `--skip-build` only when the caller deliberately accepts responsibility
for the current generated-tree identity. Qualification paths that require the
dedicated Na+Ar executable still perform their required build.

## 4. Mandatory individual-test layout

Every implemented permanent test has its own directory:

```text
UXX_descriptive_name/ or integration-tests/IXX/
├── README.md
├── test.py
└── reference/
    ├── acceptance.json
    └── independent reference files, when required
```

`test.py` must exercise production entry points or the normal linked
executable. It must write `result.json` containing its exact frozen ID and one
of `PASS`, `FAIL`, `ERROR`, or `SKIPPED`; exit codes are respectively 0, 1, 2,
and 77. The global/phase runners reject missing, malformed, ID-mismatched,
status-invalid, exit-inconsistent, or stale results.

Each `acceptance.json` must state:

```text
ID / name
purpose
production capability and entry points
configuration and inputs
independent oracle/reference
units and frame/time convention
deterministic/stochastic controls
procedure, metrics, and acceptance criteria
expected PASS and failure modes
required data and hashes
predecessor gates
artifacts
status semantics
```

Do not allocate a new ID without updating the roadmap, registry, procedure,
and status consistently. Do not create a generic test that returns PASS for an
unimplemented capability. Missing I21-I46 runners remain visible roadmap gaps.

## 5. Reference and comparison policy

Use the following oracle order:

1. closed-form analytical solution;
2. independently implemented numerical reference;
3. conservation or invariant identity;
4. authoritative external benchmark/observation;
5. reviewed regression baseline only when first-principles truth is absent.

Current production output cannot be frozen and called an independent physical
truth. A regression baseline must identify the reviewed originating run,
configuration, uncertainty, and hash. Stochastic references must record seed,
sample count/weight, confidence rule, and MPI/OpenMP layout.

When a numeric reference/model comparison is actually evaluated, preserve the
row-level values in `comparison.csv` and render those exact rows as both
`comparison.png` and `comparison.eps` with
`common/write_comparison_artifacts.py`. A skipped test produces none of these
three files: an empty or synthetic figure would make missing evidence look
like a comparison. Demonstration-only figures must say `EXPERIMENTAL` or
`VISUALIZATION ONLY` and remain outside formal gate evidence.

## 6. Component-by-component campaign

| Model component | Local kernel/oracle | Linked production gate | Required comparison/evidence | Present state |
|---|---|---|---|---|
| initialization, configuration, species, weights, time steps | U02, U25 | I01-I03 | manifests, generated inputs, exact species/weight/time-step receipt | M0 verified |
| lunar rotation, frames, Sun/Earth differential gravity, non-inertial terms | U04 analytical SPICE/transform checks | I06 trajectory through production mover | term vectors, kernel hashes, trajectory residuals | M0 verified |
| lunar gravity | U03 closed-form point gravity | I05 and I10 trajectory/convergence | state/energy/angular-momentum errors | M0 verified |
| Na radiation pressure and lunar/Earth shadow gates | U05 Combi-1997/reference geometry | I07 production mover | acceleration/trajectory and shadow classifications | M0 verified |
| Lorentz force and ion mover | U06 analytic gyro fixture | I08 dedicated Na/Na+ build | trajectory/invariant comparison | M0 verified fixture; normal build remains Na-only |
| particle injection and source registration | U22 Na kernels; U19-U21/U23 later sources | I26-I30; I29 is Na partition | injected weight/rate budgets and seeded distribution comparisons | Na partial/experimental; later species/drivers incomplete |
| surface impact and boundary classification | U01/U12 geometry and interaction kernels | I09, then I18-I20 | hit/miss/escape/domain counts, mass/event budgets | M0 classification verified; M2 bookkeeping not qualified |
| sticking, accommodation, residence, release | U12 analytical/table/Arrhenius laws | I18 stochastic production callback; I19 conservation | binomial/KS/rate comparisons, inventory closure | implemented experimentally; predecessor-gated |
| accumulated inventories, desorption, checkpoints | U12/U13 local kernels | I19, I20, later I44-I45 restart/reproducibility | global weighted conservation, checkpoint continuity | incomplete linked qualification |
| analytic spherical surface and angular refinement | U07 | I12, I13 | requested/realized resolution; area invariant; reviewed baseline | I12 implemented; M1 blocked by I11/I13 |
| LOLA LDEM ingest and triangulated surface generation | U08 native-label decode, radial displacement, topology/normals | I14 actual production particle boundary | vertex/radius/area/topology/intersection comparison; D01 hashes | implemented; formal I14 predecessor/data gates unresolved |
| ray/surface intersections and moving topographic boundary | U01, U08 | I09, I14 | analytic sphere rays and independent LOLA facet/BVH intersections | sphere baseline verified; LOLA linked qualification incomplete |
| illumination, terrain ray tracing, Earth eclipse, PSR | U09, U13 | I15, I20 | independent facet rays/D02 comparisons; preserve real gaps | local kernels implemented; linked M2 incomplete |
| simple, Diviner-backed, dynamic thermal state | U10, U11 | I16, I17 | interpolation/analytic thermal cases; D03 rows; dt/dt2/dt4; PNG/EPS | experimental; I17 numerical/hold-out gates incomplete |
| Na photoionization and survival | U14 | I21 | closed-form survival plus linked Na-to-Na+ budget | local Na kernel verified; I21 runner absent |
| multispecies photochemistry | U15, U24 | I22, I31 | D04 rates/branching and conserved atoms | unavailable/incomplete |
| electron impact and charge exchange | U16, U17 | I23, I24 | D05/distribution integrals and authoritative channel data | incomplete |
| plasma/magnetic driver | U18 | I25 | D06 native metadata, interpolation, frames and gaps | incomplete; D06 qualification required |
| He, Ne, Ar sources | U19-U21 | I26-I28 | D06/D10 source response and linked budgets | He/Ne incomplete; uniform Ar branch experimental only |
| meteoroid forcing and H2O/OH migration | U23, U24 | I30, I31 | D11/D15 timing, conversion, chemistry and budgets | incomplete |
| LOS, column-density and Na D-line brightness operators | U01, U26 | I32, I33 | independent geometry/operator fixtures before observations | generic output exists; full instrument qualification incomplete |
| Kaguya/LADEE/LAMP/LACE/PACE operators and observations | U25, U26 | I32-I40 | D07-D09/D12-D15 qualified packages, frozen selections and uncertainties | not yet validated |
| cross-species freeze, UQ, seeds, MPI/OpenMP, release | U25 plus all predecessors | I41-I46 | frozen parameters, ensembles, parallel equivalence, scorecard | not reached |

### 6.1 Particle sources

For each enabled source, test in this order:

1. deterministic rate/normalization and units at the production kernel;
2. source registration exactly once in the generated configuration;
3. seeded production injection through the normal callback;
4. position, direction, speed/energy and time distributions against an
   independent analytic CDF or authoritative driver;
5. injected statistical weight versus integrated requested source rate;
6. linked inventory/source/sink budget closure;
7. sensitivity to resolution/time step before observations.

Na PSD, thermal desorption, impact vaporization and sputtering must be reported
separately before their sum. He, Ne, Ar, meteoroid and H2O/OH sources remain
disabled or experimental until their named I26-I31 gates pass.

### 6.2 Boundary interaction, sticking and release

Use the exact production intersection and surface callback. Cover sphere and
LOLA facets, grazing rays, inside/outside origins, face/edge/vertex hits,
moving body-frame transforms, physical escape and numerical-domain exit.
For every species and surface-temperature regime, record impact weight,
sticking Bernoulli probability, accommodation law, adsorbed inventory,
residence/desorption probability, released speed/direction, cold-trap state,
and physical/numerical loss. I19 must demonstrate global weighted conservation
before any surface-controlled observational validation.

### 6.3 Lunar rotation and reference frames

Freeze UTC/TDB conversion, J2000 inertial, LSO rotating and
MOON_ME_DE421 body-fixed conventions plus exact SPICE hashes. Compare the
production transform derivative with an independent finite difference, then
propagate through the same mover used by the application. Geometry-only or
term-only success cannot replace the linked trajectory tests.

### 6.4 LOLA mesh and surface ray tracing

Before I14:

1. read the native LDEM label and verify dimensions, endian, scale, missing
   values, longitude direction/range, latitude definition and frame;
2. generate a near-uniform triangulated sphere at the input resolution;
3. radially displace every vertex by the decoded LOLA altitude—never merely
   attach altitude as a plotting scalar;
4. verify vertex radii independently from native pixels;
5. verify watertightness, manifold edges, positive area, outward orientation,
   face degeneracy, normals and total area;
6. compare BVH and brute-force intersections on deterministic rays;
7. save CEA for AMPS and Tecplot only when requested;
8. initialize the normal particle boundary from the saved/loaded result and
   prove the linked mover reaches it;
9. render geometry and residual maps as diagnostic artifacts, clearly
   separated from D01 qualification and I14 status.

### 6.5 Illumination, shadow and thermal state

Test geometric night, terrain obstruction, Earth occultation and PSR as
separate classifications. Use independent facet-ray fixtures locally and D02
only after package readiness/hashes pass. Apply Diviner native calibration and
quality flags before aggregation, preserve detector geometry and gaps, and
compare bolometric brightness temperature through the matching
`epsilon*sigma*T^4` operator. I17 requires the frozen June/December selection,
thermal time-step convergence, observation-free spin-up and declared
uncertainties before untouched hold-out scoring.

## 7. Data-dependent gate preflight

Before a Dxx test, read the package README, `provenance.json`, `qa_report.json`,
native labels and the exact processed product; verify every hash and READY
state. Directory presence is insufficient. The fixed dependency map is:

```text
D01 -> I14                  D09 -> I35
D02 -> I15, I20            D10 -> I28, supports I37/I38
D03 -> I16, I17, I20       D11 -> I30, supports I39
D04 -> I21, I22, I31       D12 -> I38
D05 -> I23                 D13 -> I40
D06 -> I25, I26, I35, I40 D14 -> I33
D07 -> I32                 D15 -> I39
D08 -> I34, I36-I39
```

Raw data are immutable. Every processed value must be reproduced by a
committed script. Never infer units, frame, time system, flags, uncertainties,
instrument response or observation geometry.

## 8. Required artifacts

Every linked run preserves, at minimum:

```text
stdout.log / stderr.log
result.json
run_manifest.json
application.in and configuration snapshot
source revision and executable/build hashes
external-data and processing-script hashes
random seed and MPI/OpenMP layout
science metrics and comparison rows
resource metrics for costly runs
```

Do not overwrite a scientifically meaningful output directory. The global
runner writes `test_summary.json`; delegated phase and individual artifacts
remain in their own directories and stay authoritative.

## 9. Adding or changing a test

1. Inventory the production implementation, callers, switches and data.
2. Identify an independent oracle and write/update `acceptance.json` first.
3. Extend an existing ID only if its scientific meaning is unchanged;
   otherwise record the gap in `TEST_GAP_PROPOSALS.md` until the roadmap is
   revised consistently.
4. Implement a thin adapter or linked runner using the production path.
5. Add deterministic/stochastic controls and failure fixtures.
6. Add row-level comparison and PNG/EPS rendering when a comparison runs.
7. Register the runner in its phase engine and the global architecture audit.
8. Run `test_runner.py --phase audit` and focused tests.
9. Update READMEs and `srcMoon/STATUS.md` with commands, hashes, outcomes and
   remaining blockers.

Never relax a threshold after failure, tune on hold-out data, edit an
observation, or report a demonstration/preflight result as verification or
validation.
