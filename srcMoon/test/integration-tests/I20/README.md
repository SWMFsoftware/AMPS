# I20 — Permanent-shadow and polar cold-trapping physics

Purpose: Collocate independent D02 PSR and D03 polar temperatures with production residence statistics, warming release, and inventory sensitivity.

Run this ID through the aggregate runner so its I13 predecessor is evaluated
in the same immutable campaign directory:

```bash
python3 srcMoon/test/run_all_test.py --campaign-mode qualification --only I11 --only I12 --only I13 --only I20 --mpi-ranks 8
```

The runner always writes `stdout.log`, `stderr.log`, `application.in`,
`run_manifest.json`, and `result.json`. If I13 is not `PASS`, this test is
`SKIPPED` with `execution=not-run`; a U-series pass is not promoted to linked
evidence. Required external evidence: D02 LPSR/AVGVISIB and D03 PRP. Raw archive files are immutable,
all units/frames come from native labels, and missing package provenance or QA
is `SKIPPED`, never `PASS`.

The runner requires both package-level statuses, not merely generated CSVs.
D02 remains `PARTIAL` for the exact-product temporal/operator gap documented in
[D02_ILLUMINATION_RESEARCH.md](../../D02_ILLUMINATION_RESEARCH.md). D03 now has
scoped READY metadata for its qualified GCP/PRP/frozen RDR inputs; each run must
verify the exact required product hashes and selection. I20 also waits
for linked I15, I17, I18, and I19 PASS evidence, so it remains
`SKIPPED / NOT VALIDATED` until every earlier layer is satisfied.

An executed comparison must additionally preserve the exact row-level
`comparison.csv` and render those unchanged rows as `comparison.png` and
`comparison.eps`. A skipped/non-executed test creates none of these three
files; empty or synthetic comparison graphics are prohibited.

## Approved numerical/cohort implementation before formal qualification

[M2_NUMERICAL_PROTOCOL.md](../../M2_NUMERICAL_PROTOCOL.md) specifies the
November 13 approved extension. The development runner records passive Na/Ar
cohort survival, censored first residence, current occupancy, actual competing
release hazards, physical clocks and tag-mass closure through production
injection/surface/mover paths. U12 supplies an independent imposed-warming
kernel oracle. These are supporting conditional kinetics and accounting checks;
they do not establish that the production thermal history is qualified or
that a selected logical element is an observed PSR.

The cluster planner preserves all 288 declared cases, native site metadata,
shared logical support and frozen normalizations, with explicit cost/resource
proof required before launch. Even a full numerical-equivalence PASS leaves
I20 SKIPPED until independent D02/D03 geometry/time/operator comparison and
all earlier linked gates pass. No frozen observation selection or December
hold-out residual is altered or examined by this work.

## Implemented formal comparison path — 2026-11-16

The missing code path is now implemented in
[`m2_polar.py`](../../m2_polar.py) and
[`MoonPolarDiagnostics.cpp`](../../../MoonPolarDiagnostics.cpp).
This is implementation status, not a declaration of I20 PASS. The actual D02
temporal/operator metadata are still unavailable, M1/I13 is open, and no real
polar residual is evaluated by the development checks.

`prepare_input(template, destination, campaign_root)` first validates package
readiness, all same-campaign formal I13–I19 predecessors, and a reviewed frozen
polar contract. It stages **geometry-only** queries and saves
`polar_reference_manifest.json` and `polar_input_binding.json`. It does not
require a checksum of a future model output: the existing template/reference
files are frozen first; the fresh native model gets a separate post-runtime
manifest bound to that contract. The preflight authenticates the explicitly
reviewed pilot's complete model initial state; neither this stage nor the
preflight reads reference/model temperature residuals.

The normal native I20 adapter then calls
`Moon::PolarDiagnostics::RunPolarDiagnostics`. Rank zero evaluates the same
production thermal columns used by the application; the other requested MPI
ranks remain initialized but do not duplicate this deterministic history. This
is recorded as one thermal worker, not claimed as eight-rank thermal speedup.
Spin-up is observation-free and uses the explicitly configured interval and
unchanged convergence criteria. Every query must occur at or after its end.

`finish_run(work, receipt)` joins the real thermal export, the same-campaign
I15 `illumination_model.csv`, I19's independently recomputed full numerical
score, and all sixteen finest-grid native cohort/hazard records. It produces
normalized model-only rows and a post-runtime manifest, then calls the strict
independent scorer. Direct import of a previously retained model manifest is
also supported, with the same checks; imported status labels alone do not
provide qualification.

### Geometry-only native query format

The input setting `m2_polar_queries_file` names a CSV with this exact header:

```text
query_id,site_id,logical_surface_element,epoch,time_system,depth_m,sample_weight,time_support_id
```

Each ID is explicit and unique. `time_system` is exactly `UTC` or `TDB`;
numeric TDB epochs are seconds. Depth is metres and the weight is a
nonnegative explicit sampling measure. The query carries no observed
temperature, native PSR flag or source abundance. Queries may share a model
element; that shared support remains visible and is not counted as independent
spatial resolution.

The exporter writes `receipt.txt.I20_thermal_model.csv` with actual body-fixed
coordinates, epoch, depth, weight and model temperature. At depth zero it also
emits `epsilon^(1/4) T` as bolometric brightness temperature. A subsurface
kinetic temperature has an **empty brightness field**. It is not treated as
radiation emerging directly from that depth.

The native receipt now declares
`surface_kinetics_export_schema=native-species-residence-retention-v1`.
Five additional CSV columns report surface physics independently of whether
inventory checkpoint files are enabled:

| Column | Meaning and units |
|---|---|
| `surface_retention_timescale_s` | Explicit input reporting timescale, seconds. |
| `na_surface_residence_time_s` | Current Na residence time from the production species kernel, seconds. |
| `na_surface_thermal_retention` | Native diagnostic `residence_time >= timescale`, exactly 0 or 1. |
| `ar_surface_residence_time_s` | Current Ar residence time from the production species kernel, seconds. |
| `ar_surface_thermal_retention` | The same comparison for Ar, exactly 0 or 1. |

Only depth-zero rows carry these fields; all five are empty for a subsurface
query. Surface residence uses the **kinetic surface temperature**, including
when the independent surface observable is brightness temperature. Python
never evaluates another Arrhenius law to manufacture these values. It checks
the native diagnostic relation and preserves a literal `inf` residence time
as the valid zero-release limit. Zero residence is the instantaneous-release
limit. NaN, negative residence, nonbinary labels and mismatched timescales
are rejected. An inconsistent evaluated native label is `FAIL`; a malformed
value or declared-but-incomplete new schema is `ERROR`; an older exporter
without this capability is `SKIPPED`.

The input `cold_trap_retention_timescale_s` must equal the explicitly approved
M2 protocol's `retention_reporting_timescale_s` for this formal design. It is
not a fitted species threshold. The species-specific effective temperature
boundary follows from its existing residence law and that timescale. The
diagnostic neither overrides continuous release nor creates a permanent
reservoir. A site may be geometrically dark and fail thermal retention, or be
intermittently illuminated and exceed the timescale at a particular epoch.

The production interpolation kernel `InterpolateThermalProfileK` evaluates
only the MODEL boundary/layer-centre state. It reproduces the boundary and
exact layer centres and interpolates linearly between them. No extrapolation
beyond the deepest centre is allowed. It does not interpolate, smooth or fill
any observational gap. Existing U11 adds independent affine-profile and invalid
depth/order controls for this actual production kernel.

### Site-bound cohort evidence and separate classifications

Every one of the sixteen predeclared finest I19 seed runs now supplies its
actual `run_manifest.json`, `receipt.txt`, `application.in`, `m2_budget.csv`,
`m2_cohorts.csv`, `m2_hazards.csv` and complete
`receipt.txt.I19_initial_thermal_state.csv` descriptors. Each exact path/hash
pair must occur in the same-campaign formal I19 aggregate's
`native_source_evidence` (the separate `native_evidence` family continues to
contain exactly 288 actual receipts).
A convenient PASS label, a copied unqualified CSV or a run at a different
site cannot replace that evidence. The scorer verifies the native case,
protocol, seed, executable/build identity, MPI/OpenMP layout, start UTC/TDB,
initial/final elapsed interval, synchronized species clocks and actual input
hashes. The sixteen seed runs must share the same physical initial epoch.

The cohorts' actual applied full thermal states are additionally verified against
I19's explicitly reviewed state proof and `applied_initial_state` per-case
hashes. [M2_INITIAL_STATE_CONTRACT.md](../../M2_INITIAL_STATE_CONTRACT.md)
defines that proof. An epoch-only READY claim or an advertised-but-unapplied
periodic state cannot qualify the cohorts. Site elements must exist in the
qualified complete grid. This authenticates their physical initial-value
problem; a model state is not promoted to observational truth.

The annual I20 history has a **separate explicit**
`initial_state_qualification` path/SHA-256 descriptor in its reviewed polar
contract. It authenticates a real native I20 pilot, or explicitly reuses an
I19 proof only when the actual I20 state, mode, TDB epoch and canonical
physical settings match. There is no automatic cohort-proof fallback and no
requirement that independently declared annual/cohort intervals share a
phase. Missing I20 initial-state qualification is SKIPPED.

The fresh I20 exporter retains its complete pre-query state in
`receipt.txt.I20_initial_thermal_state.csv`. The scorer checks it against that
chosen proof, including native initialization mode/epoch, physical input,
executable and full layer grid, before annual summaries or observational
residuals are scored. Its receipt uses `initial_thermal_state_file`,
`initial_thermal_surface_elements`, `initial_thermal_layer_count`,
`initial_thermal_state_tdb_s` and shared `thermal_initialization_*` metrics.
An explicit legacy I20 spin-up can establish its own later applied epoch;
that epoch must be reported and qualified rather than relabelled as a cohort
start or assumed identical to `simulation_start_utc`.

The ordered `m2_verification_surface_elements` list must exactly match the
frozen I20 native-site registry in its preserved order. Cohort index selects
that site, so every Na/Ar checkpoint and hazard row must carry its declared
logical element. Multiple native sites sharing an element retain separate
labels and cohorts; they are **correlated model support**, not independent
spatial resolution. Missing species, changed labels/elements, duplicate
checkpoints/hazard intervals and incomplete endpoints are `ERROR`.

Controlled initial physical stock equals the configured cohort packet scale
times the native **initialized, pre-sweep** particle weight. The production
fixture retains that stock while refining emitted packet weights. The
scorer checks this arithmetic identity and constant initial stock across
checkpoints. It does not change the initial-value problem with the numerical
weight divisor. The existing independent hazard integrator continues to
check survival, censored first residence, bounds and tagged inventory closure;
its evaluated failures remain `FAIL`.

The resulting model diagnostics are:

* `surface_retention.csv` and `surface_retention.json`: exact per-query,
  per-species residence/retention, kinetic temperature, element, depth-zero
  TDB epoch and thermal support; independently supplied terrain-only PSR and
  benchmark visibility are carried with their **own** illumination support.
* `site_kinetics.csv`: exact per-seed/site/species cohort survival, restricted
  first residence, current surface/gas occupancy and loss fractions, with
  absolute and elapsed clocks.
* `conditional_kinetics.json`: retained per-seed manifests/build identities,
  independent conditional checks and shared-support interpretation.

Geometric PSR is bound to I15's
`terrain_any_disc_occurrence_fraction`, never to its joint terrain/Earth
visibility. The existing `any_disc_occurrence_fraction` retains the
benchmark's declared Earth-occultation policy. A terrain-visible site whose
sampled sunlight is eclipsed by Earth is not silently reclassified as terrain
PSR. Missing terrain-only evidence is unavailable capability; it is not
inferred from the joint mask.

These are finite-support diagnostics. In particular,
`permanent_retention_established=false` is explicit in the reporting artifacts:
instantaneous residence exceeding a finite timescale is not proof of
indefinite capture. The cohort's approved numerical horizon is also not
equated with the PRP annual window or I15 illumination window. Formal annual
comparison still requires independently qualified temporal coverage and the
correct physical/brightness temperature/depth definition. No diagnostic file
is treated as an independent observed temperature, residence time or trap
population. Numerical/cohort failures stop scoring before reference residuals
are opened; evaluated diagnostic files remain available for investigation.

Local infrastructure coverage runs with:

```bash
/usr/bin/python3 -B srcMoon/test/test_m2_polar.py
```

The ordinary `integration-tests/I20/test.py` runner also executes both
`test_m2_polar.py` and `test_m2_operator_checks.py` automatically before its
shared scientific dispatcher. Thus the literal I20 entry in `run_all_test.py`
includes these local contracts even when I13 or D02 keeps formal science
SKIPPED. Commands, source hashes and individual stdout/stderr are retained
under `I20/protocol_checks/`. A subordinate failure creates an I20
infrastructure `ERROR` result and prevents native/scientific execution;
infrastructure success leaves every original scientific gate unchanged.

Its temporary metadata/arithmetic fixtures check the masks, zero/infinite
residence limits, site/species/depth/clock/seed/build binding and
PASS/FAIL/ERROR/SKIPPED control flow. They are not lunar observations or formal
I20 evidence, and they do not read December hold-out residuals.

### Frozen site selection and comparison contract

`m2_formal_contract_file` may name the direct
`srcMoon-I20-polar-contract-v1` document, or the shared reviewed
`srcMoon-M2-formal-contract-set-v1` dispatch bundle with a hashed `contracts.I20`
entry. The shared resolver records both selected-document and bundle identity.
An unreviewed bundle or child is not approval.

The direct document requires the following independently supplied evidence:

| Field | Required meaning |
|---|---|
| `status`, `selection_frozen`, `review` | Actual investigator approval before model evaluation, named reviewer and timezone-qualified approval time; code never creates approval. |
| `frame`, `coordinate_convention` | Native `MOON_ME_DE421`, planetocentric/east-positive geometry. |
| `application_template_sha256` | Hash of the already reviewed source `.in`; the fresh staged input also gets its own exact runtime hash. Formal staging preserves source/AMR/thermal settings and changes output ownership/diagnostic paths only. |
| `candidates`, `collocation_contract`, `site_registry` | Path/SHA-256 descriptors for independent native metadata, explicit reviewed collocation distance/rule, and the unchanged frozen selection. |
| `reference_rows` | Hashed independently processed native D02/D03 comparison rows with source product/record identities and original quality/missing-uncertainty states. |
| `thermal_queries`, `thermal_supports` | Hashed geometry-only schedule and explicit per-quantity support: interval, endpoint convention, query membership/weights, depth, temperature kind and independent sampling/coverage evidence. |
| `operator_qualification` | Real independent operator evidence, exact site/collocation hashes, current scorer hash, and the supporting evidence-file hashes. |
| `criteria` | Predeclared reviewed quantitative limits for each distinct observable. No polar thermal tolerance is silently borrowed from I16's non-PSR 15 K gate. |

Native metadata selection is reproducible with `prepare_site_registry`. It
applies the already approved SHA-256 ranking of compact JSON
`[product_id,line_1_based,sample_1_based]`, separately in each hemisphere,
absolute-latitude band and native PSR class. It preserves the eight-neighbour
edge exclusion, three available sites per stratum, missing strata and shared
logical elements. The resulting candidate has
`REVIEW_REQUIRED_NOT_FROZEN`; the script cannot approve its own selection.
No model residual/value is allowed in the selection inputs.

The following four observables remain distinct throughout native aggregation,
reference matching and plotting:

* `any_solar_disc_visibility_occurrence_fraction`, dimensionless, surface;
* `permanent_shadow_classification`, binary, surface;
* `annual_maximum_surface_temperature_K`, K, surface;
* `annual_average_temperature_at_2cm_K`, K, depth 0.02 m.

An instantaneous sample is not an annual statistic. Mean visible disc area
is not occurrence frequency. Geometric permanent shadow is not permanent
kinetic trapping. Native temperature kind must match the model observation
operator; surface brightness and subsurface kinetic temperatures cannot be
substituted for each other. The scorer recomputes every thermal aggregate
from the actual native samples and frozen weights/query membership, and
requires qualified temporal coverage. Missing support remains SKIPPED.

### Commands and expected results

Run the fast protocol/arithmetic checks independently of AMPS compilation:

```bash
python3 -B srcMoon/test/test_m2_polar.py
```

These temporary metadata/arithmetic fixtures exercise failure semantics and
do not produce an observational plot. The global M2 development runner also
includes them as supporting infrastructure checks.

The separate bounded native exporter checks use an actual retained Na/Ar
executable and the existing development receipt:

```bash
python3 -B srcMoon/test/run_m2_operator_checks.py \
  --amps /path/to/retained/amps-moon-m2-na-ar \
  --input srcMoon/examples/m2_development.in \
  --output-dir test_output/srcMoon/operator-wiring-NEW \
  --mpi-ranks 8
```

This supporting runner intentionally uses a **1 ms forcing replay**, not a
lunation. It preserves the 0.01 K/10000 J m^-2 spin-up convergence limits,
samples two real production columns at 0/0.02 m and two explicit epochs, and
checks identity/bounds/finite state and the grey surface operator. It also
executes eight predeclared I15 terrain/SPICE geometry queries. This validates
new API/production wiring at bounded cost; it establishes no lunar periodic
initial state, annual climatology, statistical convergence or observational
agreement. It writes separate support results with formal I15/I20 SKIPPED,
all native streams/input/query/receipt/build/data hashes and resource logs.
The enclosing development deadline still applies; no gate is relaxed if the
fixture fails or exceeds resources.

For a later fully qualified campaign, include the complete predecessor chain:

```bash
python3 -B srcMoon/test/run_all_test.py --campaign-mode qualification \
  --only I11 --only I12 --only I13 --only I14 --only I15 --only I16 \
  --only I17 --only I18 --only I19 --only I20 --mpi-ranks 8
```

That command cannot supply a missing frozen contract, D02 metadata, I17
review or long numerical qualification. With current prerequisites, formal
I20 remains SKIPPED and no observational plot is produced.

Once real gates exist, evaluated I20 rows retain `comparison.csv/.png/.eps`,
`polar_rows.json`, `conditional_kinetics.json`, `polar_preflight.json`,
`polar_result.json`, exact native model rows and all source/runtime bindings.
A failed conditional kinetic identity stops before observational residuals
are opened. Corrupt hashes/identities are ERROR, unavailable scientific input
SKIPPED, and an executed frozen-criterion miss FAIL.
