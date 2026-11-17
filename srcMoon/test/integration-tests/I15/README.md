# I15 — Topographic illumination and shadowing

Purpose: Compare the production terrain-horizon and Earth-eclipse service with independent D02 polar visibility/PSR evidence and published temporal sampling.

The individual `I15/test.py` runner first executes `test_m2_illumination.py`,
`data_processing/test_i15_spatial_support.py` and `test_m2_contracts.py`.
They use temporary metadata/arithmetic fixtures only, so their success is
infrastructure evidence and cannot qualify D02 or an illumination comparison.
The literal `run_all_test.py` registry includes these checks even while
missing timing metadata or I13 keeps scientific I15 `SKIPPED`. Each command,
stdout/stderr and infrastructure-only summary is retained under
`I15/protocol_checks/`. A subordinate failure instead writes per-ID
`result.json` as `ERROR` and prevents the shared scientific runner launching
a linked model; successful checks leave its science gates unchanged.

## Terrain-only PSR correction contract — 2026-11-16

This contract precedes the scorer correction. The native exporter already
reports terrain-only and jointly terrain/Earth-visible flags separately.
At inventory time the scorer derived geometric PSR from joint occurrence, so a site
whose terrain-visible Sun is Earth-eclipsed at every selected epoch can be
mislabelled as terrain-shadowed. Correct only that semantic join: preserve
the reviewed Earth inclusion policy for the benchmark visibility/history and
derive geometric PSR from terrain-only occurrence on the same declared epochs.

The existing I15 control/arithmetic test is extended with a predeclared site
that is terrain-visible but Earth-eclipsed throughout its finite history.
Its joint benchmark occurrence must be zero while terrain occurrence is one
and `geometric_psr` is zero. A separate direct epoch/timing case retains the
original 99%/native-area/timing gates. These are temporary Boolean/arithmetic
fixtures, not physical observations or a second geometry model. Require both
native flags to agree with their respective exported projected-area masks;
missing terrain flags must be ERROR, never inferred from joint visibility.

Inputs/units/frames remain the existing I15 contract: native mask fractions
are dimensionless, site/frame/epoch identities and hashes are exact, no RNG is
used by the join, and the native production path is unchanged. Preserve both
terrain and joint occurrence in `illumination_model.csv`, with explicit Earth
policy and finite-sampled-history scope. The corrected geometric field feeds
I20; it never means infinite-time darkness or permanent kinetic trapping.
Missing authoritative temporal support stays SKIPPED, corrupted mask/identity
data stay ERROR, and evaluated original criterion misses stay FAIL. No new ID,
new observational tolerance, changed selection or permanent-PSR assertion is
introduced.

The [D02 research](../../D02_ILLUMINATION_RESEARCH.md) distinguishes native
any-disc occurrence visibility from mean visible-disc area/incident flux.
Exact 201608 calculation epochs and sampling controls remain unresolved;
general publication values cannot substitute for them. Average maps cannot
provide sunrise/sunset or individual-epoch truth. D02 stays PARTIAL and the
production extended-disc/temporal operator remains unqualified. The package
audit automatically preserves native metadata and runs ten protocol checks;
these are supporting preflight, never linked I15 PASS.

Run this ID through the aggregate runner so its I13 predecessor is evaluated
in the same immutable campaign directory:

```bash
python3 srcMoon/test/run_all_test.py --campaign-mode qualification --only I11 --only I12 --only I13 --only I15 --mpi-ranks 8
```

The runner always writes `stdout.log`, `stderr.log`, `application.in`,
`run_manifest.json`, and `result.json`. If I13 is not `PASS`, this test is
`SKIPPED` with `execution=not-run`; a U-series pass is not promoted to linked
evidence. Required external evidence: D01 LOLA topography; D02 AVGVISIB/LPSR. Raw archive files are immutable,
all units/frames come from native labels, and missing package provenance or QA
is `SKIPPED`, never `PASS`.

An executed comparison must additionally preserve the exact row-level
`comparison.csv` and render those unchanged rows as `comparison.png` and
`comparison.eps`. A skipped/non-executed test creates none of these three
files; empty or synthetic comparison graphics are prohibited.

## Completed configurable production and formal comparison paths

The November 16 implementation adds `EvaluateSolarDisc` and equal-epoch
`TemporalVisibility` to the production `MoonSurfaceEnvironment` service.
`illumination_terrain_solar_model=uniform-disc` explicitly selects the new
terrain/opaque-Earth ray integration. Both `illumination_disc_radial_points`
and `illumination_disc_azimuth_points` must be supplied; the azimuth count is
at least four and divisible by four. These are numerical controls, not an
automatically qualified angular resolution. Legacy inputs retain point-Sun
forcing. `illumination_aberration_correction` explicitly selects the supported
SPICE convention (`NONE`, `LT`, `LT+S`, `CN`, `CN+S`); it must agree with the
authoritative reference contract.

The quadrature samples equal projected-disc areas. It retains terrain-only,
Earth-only and **joint** visible fractions rather than multiplying correlated
masks. Irradiance additionally integrates each ray's positive incidence and
the production directional Bond albedo, with the solid-angle Jacobian.
The thermal boundary receives that integrated absorbed flux through the same
conduction solver. Uniform radiance, opaque solid Earth and absent atmospheric
transmission remain explicit model limitations. A visible sky ray below an
elevated observer's local horizontal is not discarded as facet self-shadow;
I15 uses actual terrain intersections for that observer geometry.

The following observables are intentionally distinct:

* sampled any-disc occurrence at the explicitly listed epochs (D02 AVGVISIB);
* temporal mean visible-disc projected area;
* temporal mean projected incident irradiance.

Terrain-only sampled occurrence is retained separately from the benchmark's
Earth-policy joint occurrence. `geometric_psr` is exactly the zero-occurrence
classification of the terrain-only history on the declared sampled epochs.
A wholly Earth-eclipsed but terrain-visible history therefore has joint
occurrence zero, terrain occurrence nonzero and geometric PSR false. This
finite-history bit does not establish infinite-time permanent shadow or a
species' kinetic retention probability.

A finite ray grid can miss an arbitrarily thin visible sliver. Matching two
sparse grids alone is therefore insufficient for formal any-disc validation.
The formal driver requires reviewed independent angular/sliver and temporal
error evidence tied to real refinement runs before scoring. U09 analytical
full-disc/half-disc/mask/radiometry checks support implementation verification;
they do not establish the missing D02 product-specific timing or error proof.

### Reviewed reference schema

The `.in` setting `m2_formal_contract_file` names either a direct
`srcMoon-I15-reference-contract-v1` JSON document or the shared reviewed
`srcMoon-M2-formal-contract-set-v1` bundle whose `contracts.I15` descriptor
binds that document by path and SHA-256. Relative paths resolve against the
owning input/JSON document; quoted paths preserve spaces, case and literal
`!`. No READY contract is generated automatically.

The I15 document requires:

| Field | Required meaning |
|---|---|
| `schema`, `status`, `frozen_before_model_evaluation` | Schema above, `READY`, and literal `true`. |
| `review` | Explicit `approved: true`, named `reviewer`, timezone-qualified `date_utc`. Review must be supplied; code never approves evidence. |
| `selection_role` | Explicit `calibration` or `hold-out`, frozen before AMPS comparison. |
| `independent_source_ids` | Nonempty unique native benchmark/product IDs used by epoch/site/timing rows. |
| `operator` | `frame=MOON_ME_DE421`, `time_system=UTC` or `TDB`, positive per-site `sample_count`, exact `epoch_list_endpoint_convention`, `cadence_description`, `observer_height_m`, `height_direction=radial` or `facet-normal`, supported `aberration_correction`, `earth_occultation=excluded` or `opaque-sphere`, `solar_radius_m`, `solar_radiance_model=uniform`, and a frozen `time_support_id`. Included Earth also requires its explicit `opaque_earth_radius_m`. |
| `authoritative_documents` | Nonempty path/SHA-256 records for actual exact-product producer/native documentation. |
| `operator_source_locators` | For each required operator field, a valid `document_index` plus exact nonempty `locator`. General publication settings cannot silently become product-generation settings. |
| `geometry_queries`, `epoch_reference`, `site_reference`, `timing_reference`, `resolution_qualification`, `spatial_support` | Each a path/SHA-256 descriptor; no path-only evidence or inferred fallback. |

The geometry CSV fields are exactly `query_id,site_id,latitude_deg,longitude_east_deg,epoch,time_system,observer_height_m,height_direction,include_earth`.
IDs are explicit stable tokens; longitude is positive east in `[0,360)`.
All epochs of one site retain the same coordinates/height, every site has the
same complete schedule and actual epochs cannot be duplicated by alternate
numeric/fractional spellings. UTC uses explicit ISO calendar/`Z` spelling;
possible `23:59:60` spelling is retained for CSPICE's leap-second authority.
Python never substitutes POSIX arithmetic for UTC-to-TDB conversion.

The remaining independent tables require these columns:

* Epoch: `query_id,reference_visible,source_id`, with exact binary visibility.
* Site: `site_id,reference_average_occurrence,reference_psr,psr_boundary_ambiguous,pixel_area_m2,psr_resolution_allowance_m2,logical_surface_element,time_support_id,native_product_id,native_line,native_sample,source_id,latitude_deg,longitude_east_deg`.
* Timing: `event_id,before_query_id,after_query_id,reference_epoch_tdb_s,event_kind,maximum_timing_error_s,source_id`; `event_kind` is `sunrise` or `sunset`. Timing tolerances must also match independently documented, frozen support evidence.

`resolution_qualification` must bind the **actual** `source_input_sha256`,
`executable_sha256`, `production_source_signature_sha256`, canonical
`geometry_queries_sha256`, and selected `disc_radial_points` /
`disc_azimuth_points`. Its literal `angular_and_temporal_resolution_qualified`
is true and `status=PASS`. At least two distinct actual `native_runs` must
contain hashed `receipt`, `manifest` and `model_outputs` records. Native
test/layout, physical terrain/radii/aberration controls, geometry identity,
query count, ray counts and full output joins are checked. The selected grid
must be represented. `independent_resolution_assessment` identifies a hashed
JSON proof with `status=PASS`, explicit `any_disc_sliver_error_qualified` and
`temporal_event_brackets_qualified`, the same geometry hash, an independent
`method` (`closed-form-bound` or `independent-numerical-bound`), and hashed
supporting `evidence`. A boolean attached to an arbitrary file is insufficient.
These numerical review records cannot be inferred from AMPS agreement with
the observational hold-out.

### Spatial support and unchanged gates

Sparse site probes cannot establish PSR areas or boundaries. `spatial_support`
therefore requires `status=READY`, `frozen_before_model_evaluation=true`,
`scope=resolved-regional-native-cell-coverage`, and explicit
`area_measure=native-polar-stereographic-map-plane`, `native_pixel_width_m`,
`native_pixel_height_m`. Its `domains` identify product and integer
`line_min,line_max,sample_min,sample_max` plus the corresponding
`scoring_line_min,scoring_line_max,scoring_sample_min,scoring_sample_max`.
Every native cell is present exactly once and the scoring rectangle has a
complete one-cell neighbour padding ring. This states a predeclared regional
scope; it makes no unsupported whole-Moon area claim.

Eight-neighbour native PSR changes independently identify ambiguous boundary
cells. Their map-cell areas supply the resolution allowance once; interior
rows have zero allowance. Arbitrarily assigning a tolerance to every site
would permit all classifications to disagree and is rejected. Interior
mismatches fail the boundary criterion. `timing_tolerances_s` in the support
proof binds each event's previously documented acceptance tolerance; no fixed
error is invented. Native pixel area is quantization/resolution, not a claimed
measurement uncertainty. Real gaps and unknown time support remain missing.

The original visibility gate stays **at least 99%** agreement with independent
epoch classifications. PSR disagreement must satisfy its native-resolution
bound, and every independently bounded timing bracket must meet its frozen
criterion. Model event brackets are not substituted for unknown observed
sunrise times. Mean-occurrence comparisons are also exported as diagnostics;
no new agreement tolerance is invented for them. All conditions are conjunctive.

The native adapter emits per-rank geometry CSVs and echoes input epoch/time
system, site and coordinates. The scorer requires exact joins, source/grid
identity and consistent area/occurrence/mask values. Missing metadata or real
numerical proof is SKIPPED; corrupted hashes, duplicate/ragged rows, wrong
epochs/frames/coordinates/controls are ERROR; evaluated criterion misses are
FAIL. `comparison.csv/.png/.eps` and `illumination_model.csv` retain only
actually evaluated rows. No canonical observation file is modified.

`illumination_model.csv` preserves the existing
`any_disc_occurrence_fraction` benchmark quantity and its explicit
`earth_occultation_policy`. It additionally records
`joint_any_disc_occurrence_fraction`, `terrain_any_disc_occurrence_fraction`
and `mean_terrain_visible_area_fraction`; the existing
`mean_visible_area_fraction` remains the joint projected-area statistic.
`geometric_psr` is bound to the terrain-only occurrence, with
`geometric_psr_scope=terrain-dark-over-declared-sampled-epochs` and unchanged
`time_support_id`/`epoch_count`. I20 consumes this terrain-only classification
and keeps the benchmark visibility separate. There is no inference of a
missing terrain flag from joint darkness or area. Both native Booleans must
agree with their respective native area masks or the result is ERROR.

### Protocol tests and current limitation

```bash
python3 -B srcMoon/test/test_m2_illumination.py
```

The temporary metadata/arithmetic fixtures test missing exact-product
metadata, corrupt hashes, source/input/executable/grid binding, native epoch
and geometry joins, the distinction between occurrence and area, malformed
CSV/JSON, UTC/leap spelling, unchanged 99% and timing gates, complete regional
support and rejection of arbitrary PSR allowances. They never execute a fake
binary, invoke MPI, render synthetic observations, inspect canonical residuals
or qualify I15. U09 and linked production runs remain separate gates.
The terrain/occultation regression also exercises a terrain-visible,
all-epochs-Earth-eclipsed site through the actual scorer, plus refusal of a
missing/inconsistent terrain flag and unavailable temporal support. Its
temporary values are Boolean/arithmetic control inputs, never observations
or a new permanent-shadow benchmark.

The user confirmed on November 16 that no additional producer metadata are
available. The configurable production/export/scoring code is implemented;
formal D02 comparisons remain **SKIPPED / NOT VALIDATED** until authoritative
operator settings, independent timestamped truth and reviewed numerical/spatial
support exist, in addition to same-campaign I13 and the original predecessors.

The reproducible explicit-region preparation CLI and schema are documented in
[I15_SPATIAL_SUPPORT.md](../../data_processing/I15_SPATIAL_SUPPORT.md). It reuses
the frozen native decoder/projection and writes complete regional cells,
padding/edge flags, provenance and QA outside the canonical data tree. It
requires a reviewed region/area/time-support contract; no real region or
generation epoch was chosen by this implementation work.
