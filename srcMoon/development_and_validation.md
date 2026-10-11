# AMPS Lunar Exosphere — Development and Validation Roadmap

> **Purpose:** repository-level roadmap for evolving `srcMoon` from a legacy spherical lunar-exosphere application into a numerically verified, topography-aware, thermally realistic, multi-species, observation-validated AMPS application.

**Roadmap baseline:** U01–U26 stand-alone tests are the local verification layer. I01–I46 are the linked AMPS integration/validation layer. A source/preflight PASS is not a linked-runtime PASS.

> **Non-negotiable data rule:** no synthetic, mock, assumed, hand-tuned, or model-generated observational data may be used to claim validation. If an authoritative data product, required production capability, or independent comparison is unavailable, the corresponding result is `SKIPPED / NOT VALIDATED`, never `PASS`.

## 1. How to use this roadmap

Work through the phases in order. Do not begin calibration to observations until the linked baseline, numerical convergence, surface geometry, and mass/chemistry budgets are qualified. Each roadmap item is a release gate with explicit development work, data requirements, tests, and exit criteria.

A capability may be labeled:

- **VERIFIED** — analytical/unit and linked numerical gates pass.
- **VALIDATED** — verified **and** an independent observational comparison passes with frozen parameters.
- **EXPERIMENTAL** — implementation exists but one or more required gates/data are incomplete, skipped, or failed.
- **DISABLED** — not built/enabled in the production configuration.

Status semantics for the runner must remain distinct:

- `PASS`: acceptance criteria were evaluated and met.
- `FAIL`: the capability/data existed, the test ran, and the acceptance criterion was not met.
- `ERROR`: execution, parsing, build, or test infrastructure failed.
- `SKIPPED`: required data/capability/configuration was unavailable; this is not success.

## 2. Current entry point

| Component | Current roadmap status |
|---|---|
| U01–U26 | Implemented stand-alone verification layer; keep fast and independent of the full AMPS link. |
| I01–I46 | Integration/validation framework exists; actual linked numerical/science status must come from a fully configured Moon executable, not preflight alone. |
| Production surface | Spherical baseline remains essential; LOLA mode must become the actual particle-surface boundary before I14 can pass. |
| Production thermal state | Analytical/cosine baseline is retained for regression; Diviner-backed and dynamic thermal modes must be promoted separately. |
| Observation bundles | Must be acquired from authoritative archives, checksummed, processed by committed scripts, and frozen before scoring. |

## 3. Validation architecture and common engineering rules

- Tests must call the same production kernels, dispatch macros, geometry services, surface callbacks, and movers used by the science application. Do not maintain a second test-only physical formula.
- Every external raw product is immutable. Store it under a `raw/` directory with URL/DOI/product ID, access date, native label/metadata, and SHA-256 checksum.
- Every processed validation file must be reproducible from a committed script. Manual spreadsheet editing is not an acceptable processing pipeline.
- When a machine-readable archive does not exist and a paper figure must be digitized, record paper/DOI, page/figure, axis calibration, digitization tool, digitization uncertainty, and an independent cross-check.
- Declare calibration/training cases and hold-out cases before final scoring. Never retune parameters after inspecting hold-out residuals.
- Numerical, observational, and Monte Carlo uncertainties must all be smaller than or explicitly propagated into the scientific acceptance threshold.
- Every linked run must write a manifest containing source revision, build identity, input hashes, data-product hashes, species table, mesh/time-step settings, random seed, MPI/OpenMP layout, and enabled physics.
- Acceptance criteria are changed only through a reviewed roadmap/test revision, never simply because a run failed.

## 4. Milestones

| Milestone | Meaning | Exit condition |
|---|---|---|
| M0 | Linked baseline qualified | I01–I10 run in the actual linked AMPS Moon executable; no unresolved compile/API mismatches; provenance frozen. |
| M1 | Numerics qualified | I11–I13 pass; time-step and AMR error are below declared science tolerances; spherical reference is frozen. |
| M2 | Real lunar environment qualified | I14–I20 pass; LOLA production boundary, terrain illumination, Diviner/dynamic thermal state, adsorption/inventory, and cold trapping are operational. |
| M3 | Chemistry qualified | I21–I24 pass; photon, electron-impact, and documented charge-exchange channels are connected and budget-closing. |
| M4 | Drivers and sources qualified | I25–I31 pass; measured plasma/meteoroid forcing and Na/He/Ne/Ar/H2O/OH source modules are operational. |
| M5 | Species/phenomena validated | The independent observational tests required for each advertised species/phenomenon in I32–I40 pass. |
| M6 | Cross-species parameter set frozen | I41 passes with explicit training/hold-out data and no post-hoc retuning. |
| M7 | Release qualified | I42–I46 pass; UQ, seed convergence, parallel reproducibility, provenance, and the capability scorecard are complete. |

## 5. Recommended validation-data directory contract

```text
moon_validation_data/
├── baseline/
│   ├── README.md
│   ├── provenance.json
│   └── SHA256SUMS.txt
├── lola/
│   ├── raw/
│   ├── scripts/
│   ├── lon.txt
│   ├── lat.txt
│   ├── height_m.txt
│   ├── control_points.csv
│   ├── illumination_reference.csv
│   ├── provenance.json
│   └── SHA256SUMS.txt
├── diviner/
│   ├── raw/
│   ├── scripts/
│   ├── temperature_table.csv
│   └── eclipse_thermal_history.csv
├── photochemistry/
├── electron_impact/
├── drivers/
├── comparisons/
├── model/
└── campaign/
```

For every data subdirectory:

- `raw/` contains immutable downloads.
- `scripts/` contains only reproducible converters/collocators.
- `README.md` records source, product/version, rationale, and exact processing command.
- `provenance.json` records product identifiers, URLs/DOIs, access date, units, coordinate convention, transformations, and hashes.
- `SHA256SUMS.txt` hashes raw and frozen processed products.
- Processed comparison CSVs contain observation/product ID, UTC, geometry, observed value, uncertainty/quality, model value, and calibration/hold-out flag where applicable.

## 6. Codex instructions for preparing validation data and test inputs

Codex can be used after authoritative raw data have been downloaded manually. Its role is to make the preparation **reproducible, auditable, and directly consumable by the U/I test suite**. Codex must automate data handling; it must not manufacture scientific evidence.

OpenAI Codex CLI supports repository instructions through `AGENTS.md`. For this project, keep the science requirements in this roadmap and place the invariant operational rules below in the AMPS repository-root `AGENTS.md`.

### 6.1 Repository-level `AGENTS.md` rules for lunar validation

Recommended block:

```markdown
# srcMoon / lunar validation rules

When working on srcMoon or moon_validation_data:

1. Read srcMoon/development_and_validation.md and the README plus
   reference/acceptance.json for every affected test.
2. Treat moon_validation_data/**/raw/** as immutable.
3. Never create synthetic, mock, assumed, fabricated, or model-derived
   observations and use them as validation data.
4. Never invent uncertainties, quality flags, timestamps, instrument geometry,
   scale factors, coordinate metadata, reaction branches, or missing measurements.
5. If required authoritative metadata/data are missing, leave the test
   SKIPPED / NOT VALIDATED and explain what is missing.
6. Decode every scientific product from its native label/XML/CDF/FITS metadata.
   Native metadata overrides values copied into notes or old scripts.
7. Every processed product must be reproduced by a committed script.
8. Preserve raw units until a documented conversion step.
9. Record every unit, coordinate-frame, longitude, projection, time-system,
   quality-filter, interpolation, binning, and background-subtraction operation.
10. Preserve real gaps. Do not smooth, fill, extrapolate, or re-normalize data
    unless the target test specification explicitly requires that operation.
11. Never modify an acceptance threshold because a run failed.
12. Never convert SKIPPED, ERROR, or source-preflight success into a physics PASS.
13. Freeze calibration and hold-out sets before final scoring.
14. Do not tune model parameters on hold-out observations.
15. Linked tests must exercise the normal production AMPS path; do not create a
    second test-only physics implementation.
16. Add detailed comments to preprocessing code and update the relevant README.
17. At completion report raw files, product IDs/versions, checksums, output files,
    transformations, filtering statistics, unresolved issues, tests run, and
    PASS/FAIL/ERROR/SKIPPED status.
```

### 6.2 Standard Codex data-preparation workflow

Codex should follow these steps for **every** external dataset.

#### C1 — Read the target requirements

Before processing, read:

```text
srcMoon/development_and_validation.md
srcMoon/test/integration-tests/<TEST_ID>/README.md
srcMoon/test/integration-tests/<TEST_ID>/reference/acceptance.json
moon_validation_data/<dataset>/README.md          # if already present
```

Identify:

- the production capability being validated;
- the exact processed files expected by the test;
- required units, coordinates, timestamps, and uncertainties;
- whether the data are calibration, hold-out, or an independent reference;
- which U-tests must pass before the linked I-test is meaningful.

#### C2 — Inventory the immutable raw package

Create:

```text
moon_validation_data/<dataset>/raw_inventory.json
moon_validation_data/<dataset>/raw/SHA256SUMS.txt
```

For each raw file record:

- relative path;
- byte size;
- SHA-256;
- mission/instrument;
- product ID;
- product version/revision;
- calibration level;
- source archive/DOI;
- native file type;
- native units;
- coordinate/projection metadata;
- time system;
- documented fill values;
- quality flags;
- access date when known.

Unknown metadata must be written as `unknown`; Codex must decide whether that prevents scientific use rather than guessing a value.

#### C3 — Inspect native metadata before science values

Before decoding arrays or records, verify the native label/header/attributes:

- dimensions and record counts;
- integer/float representation;
- endianness;
- scale and offset;
- units;
- fill/special constants;
- map projection and coordinate frame;
- longitude direction/range;
- latitude definition;
- UTC/TT/TDB/ET or other time convention;
- calibration state;
- instrument mode;
- quality flags.

If raw size/records conflict with metadata, stop with `ERROR`. Do not silently reinterpret the product.

#### C4 — Write the processing script first

All transformations must live under:

```text
moon_validation_data/<dataset>/scripts/
```

A processing script must:

- accept explicit input/output paths;
- print product identity and raw hashes;
- fail on unexpected metadata;
- report raw/retained/rejected sample counts;
- preserve documented gaps and invalid samples;
- produce deterministic outputs;
- write provenance and QA information;
- contain detailed comments explaining physical and coordinate transformations.

Do **not** manually edit the generated CSV/TXT files afterward.

#### C5 — Categorize every transformation

Codex must label each processing operation as one of:

**Structural / normally allowed**
- decompression;
- binary/CDF/FITS/PDS decoding;
- axis reordering;
- documented unit conversion;
- documented timestamp conversion;
- documented coordinate/frame conversion;
- selecting named calibrated variables.

**Filtering / selection**
- documented quality flags;
- documented instrument mode;
- declared altitude/local-time/event selections.

**Information-changing / requires explicit justification**
- interpolation;
- regridding;
- smoothing;
- gap filling;
- extrapolation;
- background subtraction;
- temporal/spatial binning;
- figure digitization;
- deconvolution.

The third category must be explicitly justified by this roadmap, the product documentation, or the target-test README.

#### C6 — Write complete provenance

Generate:

```text
moon_validation_data/<dataset>/provenance.json
```

Minimum schema:

```json
{
  "dataset": "...",
  "mission": "...",
  "instrument": "...",
  "source_product_ids": [],
  "source_product_versions": {},
  "source_urls_or_dois": [],
  "raw_files": {},
  "processing_script": "...",
  "processing_script_sha256": "...",
  "transformations": [],
  "quality_filters": [],
  "interpolation": "none",
  "smoothing": "none",
  "manual_edits": false,
  "synthetic_data_used": false,
  "output_files": {},
  "notes": []
}
```

#### C7 — Produce a QA report before AMPS is run

Generate:

```text
moon_validation_data/<dataset>/qa_report.json
moon_validation_data/<dataset>/qa_report.txt
```

Check, where applicable:

- array/table dimensions;
- exact row counts;
- coordinate monotonicity;
- coordinate ranges;
- unit ranges;
- fill/NaN counts;
- duplicate times;
- time monotonicity;
- data gaps and maximum gap duration;
- quality-flag distributions;
- reaction-branch sums;
- conservation/normalization;
- product control points;
- native-record spot checks;
- coordinate transform round trips;
- observation geometry.

A malformed/corrupt data package is `ERROR`, not failure of lunar physics.

#### C8 — Freeze calibration and hold-out membership

Observation tests must have a committed file such as:

```text
moon_validation_data/<dataset>/calibration_holdout.json
```

It must be created **before final scoring** and contain product IDs/epochs/records assigned to calibration and hold-out sets. Once final scoring begins, Codex must not change membership after seeing residuals.

#### C9 — Run tests in increasing scope

Preferred order:

```text
data QA
→ relevant U-test
→ one target I-test
→ parent phase
→ full linked suite when appropriate
```

Example:

```bash
python3 test/run_tests.py --test I14 --amps ../amps --data-path /path/to/moon_validation_data --output-dir test_output/I14
```

If the data are valid but a production capability is absent, keep the result `SKIPPED` and identify the missing production code. Do not alter the data to make the test runnable.

#### C10 — Codex completion report

At the end of each processing task, Codex must report:

- authoritative source;
- product IDs and versions;
- raw file paths and hashes;
- processed file paths and hashes;
- exact processing command;
- transformation list;
- filtering counts;
- gap statistics;
- QA status;
- calibration/hold-out definition if relevant;
- exact U/I commands run;
- PASS/FAIL/ERROR/SKIPPED;
- missing production capability or unresolved data issue.

### 6.3 Common output structure Codex should create

```text
moon_validation_data/<dataset>/
├── raw/
│   ├── <original downloaded files>
│   └── SHA256SUMS.txt
├── scripts/
│   └── prepare_<dataset>.py
├── raw_inventory.json
├── provenance.json
├── qa_report.json
├── qa_report.txt
├── README.md
├── calibration_holdout.json      # observation tests only
├── <processed files required by tests>
└── SHA256SUMS.txt
```

### 6.4 Dataset-specific Codex processing recipes

These recipes supplement the D01–D15 source descriptions in the next section.

#### D01 — LOLA global topography → I14

**Raw input**

```text
lola/raw/LDEM_4.IMG
lola/raw/LDEM_4.LBL
```

Codex must read the PDS label and verify grid dimensions, sample type, bit depth, endianness, scaling factor, reference-radius/offset convention, map resolution, projection, longitude/latitude ranges, and special constants.

**Processing**

1. Decode the binary strictly from the label.
2. Convert DN to elevation in meters **relative to the documented lunar reference sphere**.
3. Reorder latitude/longitude only as required for a strictly increasing AMPS grid.
4. Normalize longitude only through a documented convention change.
5. Do not smooth or resample the first production DEM.
6. Extract geographically distributed **real native grid samples** as control points.
7. Never invent control-point heights.

**Outputs**

```text
lola/lon.txt
lola/lat.txt
lola/height_m.txt
lola/control_points.csv
```

**QA**

- raw byte count matches label;
- longitude and latitude arrays are strictly monotonic;
- row/column counts agree with label;
- fill/special constants handled explicitly;
- elevation range is physically plausible;
- decoded native control points agree exactly before interpolation;
- no undocumented smoothing/interpolation.

**Run:** `U08` then `I14`. If the production boundary remains spherical, I14 must remain `SKIPPED`; Codex should then implement the LOLA production boundary as the next roadmap task without changing the data.

#### D02 — LOLA/PGDA illumination and permanent shadow → I15/I20

Use GDRPSR and GDRVIS/PGDA as **independent references**, not as the model's own terrain input.

Codex must:

- preserve north/south polar products and native projection;
- decode scale/offset from labels;
- collocate reference pixels with AMPS surface points/facets;
- identify edge/ambiguous pixels with a predeclared rule;
- generate `lola/illumination_reference.csv` containing lat/lon, product value, classification, product ID, native resolution, and edge/uncertainty flag.

QA must verify pole orientation, projection, known geographic control locations, and a distributed set of unambiguous illuminated/shadowed points.

Run `U09`, `U13`, `I15`, and `I20`.

#### D03 — Diviner temperatures → I16/I17/I20

Codex should use:

- GCP/PCP/PRP for climatological/global/polar temperature constraints;
- RDR for time-resolved sunrise/sunset/eclipses.

It must inspect PDS4 XML labels, calibrated variable names, units, fill values, quality flags, longitude convention, local-solar-time convention, and time tags.

Generate separately:

```text
diviner/temperature_table.csv
diviner/eclipse_thermal_history.csv
```

Do not synthesize eclipse histories from climatological products.

Report raw and retained counts by quality flag, cadence, gap durations, and native-node interpolation checks.

Run `U10`, `U11`, `I16`, `I17`, `I20`.

#### D04 — PHIDRATES photochemistry → I21/I22/I31

Codex must build a versioned reaction table from real PHIDRATES/literature values.

For each channel preserve:

- parent;
- products;
- reaction type;
- 1-AU rate or wavelength-dependent source;
- branch fraction;
- solar spectrum/activity assumption;
- units;
- literature/source reference.

Never invent a branch to make probabilities sum to one. Unsupported channels remain explicitly unsupported.

QA: branch sums, unit conversion, 1-AU rate checks, heliocentric scaling, shadow gating.

Run `U14`, `U15`, `U24`, then `I21`, `I22`, `I31`.

#### D05 — Voronov/Verner electron-impact data → I23

Codex should preserve the original coefficient table and method description under `electron_impact/raw/`, parse the actual atomic-stage coefficients, and map them explicitly to AMPS He/Ne/Na/Ar species.

Generate:

```text
electron_impact/rate_vs_Te.csv
```

Implement the published rate formula once in shared production/unit code. QA must check coefficient transcription, temperature domain, units, and several reference temperatures.

Run `U16`, then `I23`.

#### D06 — ARTEMIS plasma forcing → I25/I26/I35/I40

Codex should process downloaded CDAWeb CDF/CSV products such as `THB_L2_MERGED`, `THC_L2_MERGED`, or the selected instrument-specific ESA/FGM products.

It must:

1. inspect CDF variable attributes;
2. extract time, density, bulk velocity, electron/ion temperature as required, B, spacecraft position, and quality flags;
3. convert to SI;
4. transform to the model coordinate frame through documented geometry;
5. preserve gaps;
6. enforce a declared maximum interpolation gap;
7. write `drivers/plasma_driver.csv`.

QA: duplicate times, monotonic times, gap histogram, quality counts, native-record spot checks, coordinate-transform round trips, interpolation exactness at native nodes.

**Do not assume a fixed alpha/proton ratio.** If alpha flux is needed for He validation and the downloaded product does not resolve it, Codex must flag the missing alpha measurement and obtain an appropriate authoritative source before validation.

Run `U18`, `I25`, then dependent I26/I35/I40 tests.

#### D07 — Kaguya UPI-TVIS Na Level 2A → I32

Codex must process the real image/label pairs, index/catalog files, documented dark/background data, and observation geometry.

Required processing:

- decode each image from its PDS label;
- preserve Level-2A Rayleigh units;
- apply only documented dark/background corrections;
- propagate quality/mask information;
- reconstruct spacecraft LOS/tangent geometry;
- create a comparison table with image ID, UTC, LOS/bin, tangent altitude/location, observed brightness, uncertainty/quality, and calibration/hold-out membership.

Do not normalize each image to the AMPS result.

Run `U26`, then `I32`.

#### D08 — LADEE NMS He/Ne/Ar and event data → I34/I36/I37/I38/I39

Codex should preserve PDS4 XML labels and collection inventories, and process both the derived noble-gas products and the calibrated/context products needed to audit them.

For each retained record include:

- product ID;
- UTC;
- altitude;
- latitude/longitude;
- local solar time;
- species;
- abundance/density;
- uncertainty where provided;
- quality/instrument mode;
- calibration/hold-out flag.

Keep ordinary noble-gas time series separate from the water-event product.

QA: record spot checks, quality-mode counts, continuity/gaps, geometry, units, and uncertainty provenance.

Run I34/I36/I37/I38/I39 as applicable; freeze the cross-species inputs before I41.

#### D09 — LRO/LAMP helium + ARTEMIS → I35

Codex must use the LAMP product class and epoch selection defined by the chosen published helium analysis.

It must reproduce the documented:

- calibrated spectral variable;
- background subtraction;
- wavelength/spectral window;
- LOS/integration geometry;
- spatial/time aggregation.

Then collocate the result with D06 ARTEMIS forcing using a documented matching/propagation rule.

QA must reproduce at least one published summary/retrieval value and report unmatched LAMP epochs caused by ARTEMIS gaps.

Run `I35`.

#### D10 — Kaguya GRS potassium map → I28

Codex must decode the real K intensity/count-rate map from its native label, convert coordinates into the AMPS lunar convention, and create a **non-negative spatial source-weight field**.

It may normalize the spatial weights to unit integral for use as a source shape.

It must **not** infer the absolute radiogenic 40Ar production rate from K count rate alone. Absolute normalization requires an independent cited Ar constraint/radiogenic model.

QA: grid dimensions, coordinate control points, background/invalid handling, and normalized spatial integral.

Run `U21`, `I28`, then observational I37/I38.

#### D11 — IAU meteoroid streams → I30/I39

Codex must preserve stream/orbit tables, classification/status, source version, and an independent background flux/yield reference.

Generate:

```text
drivers/meteoroid_driver.csv
```

with UTC, stream ID, status/confidence, encounter speed, radiant/impact geometry, relative forcing, and uncertainty.

The driver must be constructed independently of the LADEE water-event response; do not tune event amplitudes using I39 hold-out data.

Run `U23`, `I30`, then `I39`.

#### D12 — Apollo 17 LACE → I38

Preferred data hierarchy:

1. machine-readable public table;
2. tabulated values in NASA/publication;
3. figure digitization only when no numerical source exists.

If digitization is required, Codex must store:

- report/paper ID;
- page and figure;
- axis bounds/scales;
- digitization tool/method;
- extracted points;
- digitization uncertainty;
- independent spot check or second digitization.

Keep Apollo-site geometry and local-time/lunation conversion explicit.

Run `I38`.

#### D13 — Kaguya PACE pickup ions → I40

Codex must use the real PBF1 sensor products, labels, selected published event interval, spacecraft trajectory, and plasma/magnetic context.

It must reproduce the publication's energy/mass/angle selection rather than define cuts after examining AMPS.

Generate comparison bins matching the model synthetic observation.

QA: sensor ID/mode, energy-bin edges, time matching, viewing-angle/frame transformation, and absolute-flux calibration applicability.

Run `I40`.

#### D14 — Ground-based sodium tail → I33

Preferred source hierarchy:

1. author/institution machine-readable profile;
2. numerical table/supplement;
3. figure digitization.

Preserve observation UTC, aperture/slit, projected distance, Na brightness, uncertainty, seeing/spatial resolution, and geometry.

If digitizing a figure, use the D12 traceability rules.

Synthetic observations must reproduce aperture/FOV and spatial convolution. Never independently renormalize each observed profile to the model.

Run `I33`.

#### D15 — LADEE water events + meteoroid association → I39

Codex must combine:

- the published Benna et al. event list/supplement;
- matching LADEE NMS data from D08;
- independent meteoroid-stream context from D11.

Before running AMPS, freeze:

- event IDs;
- event windows;
- onset definition;
- detection metric;
- false-positive windows;
- stream association rule.

Create:

```text
comparisons/I39_ladee_water_events.csv
```

Keep the underlying NMS time series separate.

QA must cross-check every event against the publication and must not add/remove events after seeing AMPS output.

Run `I39`.

### 6.5 Codex handoff prompt for raw-data preparation

Use this template after manually downloading a dataset:

```text
Read srcMoon/development_and_validation.md, the README for test <TEST_ID>,
and <TEST_DIR>/reference/acceptance.json.

The authoritative raw files are already in:

    moon_validation_data/<DATASET>/raw/

Prepare the complete real-data package required by <TEST_ID>.

Rules:
- Use only the downloaded authoritative data and native metadata.
- Never generate synthetic, mock, assumed, fabricated, or model-derived
  observational values.
- Do not modify raw/.
- Inspect label/XML/CDF/FITS metadata before decoding science values.
- Verify dimensions, representation, scale/offset, units, fill values, time
  system, coordinates, calibration level, instrument mode, and quality flags.
- Create/update a deterministic script under
  moon_validation_data/<DATASET>/scripts/.
- Generate raw_inventory.json, provenance.json, qa_report.json, qa_report.txt,
  README.md, SHA256SUMS.txt, plus the processed files required by <TEST_ID>.
- Preserve real gaps.
- Do not interpolate, smooth, extrapolate, gap-fill, regrid, or renormalize
  unless the roadmap or test README explicitly requires it.
- Print raw/retained/rejected sample counts for every filter.
- Add detailed comments to all new processing code.
- Run the relevant U-test first when one exists.
- Then run <TEST_ID> if the required production AMPS capability exists.
- Never weaken an acceptance criterion.
- Missing capability/data remains SKIPPED/NOT VALIDATED, not PASS.
- Finish with a report of inputs, hashes, outputs, transformations, filtering,
  QA, exact commands, test status, and unresolved issues.
```

### 6.6 Codex handoff prompt when data are ready but production capability is missing

```text
Read srcMoon/development_and_validation.md and the README/acceptance criteria
for <TEST_ID>. The authoritative validation-data package is already prepared
and must not be modified.

Determine the exact production capability that keeps <TEST_ID> SKIPPED or
prevents it from reaching the acceptance check. Implement that capability in
the normal srcMoon/AMPS production path.

Constraints:
- do not create a test-only physics path;
- preserve the existing baseline mode;
- expose the new capability through an explicit production configuration;
- use normal AMPS mover/boundary/chemistry/source/sampling dispatch;
- add detailed code comments and README documentation;
- run relevant U-tests, then <TEST_ID>, then its parent phase;
- do not alter validation data or acceptance thresholds to obtain PASS;
- if the model disagrees with the real observation/reference, report FAIL and
  diagnose the physics rather than editing the evidence.
```

### 6.7 Codex handoff prompt for frozen observational scoring

```text
Use the frozen processed observation package and frozen
calibration_holdout.json for <TEST_ID>.

Do not alter hold-out membership and do not tune on hold-out data.

Reproduce the exact observation operator, including as applicable:
- time and spacecraft geometry;
- LOS/tangent geometry;
- altitude/local-time binning;
- field of view/aperture/slit;
- spectral response;
- spatial convolution;
- integration interval;
- temporal averaging.

Propagate observational uncertainty and model numerical/Monte-Carlo
uncertainty. Write point-by-point residuals and machine-readable summary
metrics. Report calibration and hold-out metrics separately.

If held-out observations fail the declared acceptance criteria, return FAIL and
diagnose the likely physical cause. Do not renormalize the observation or
retune parameters after inspecting the hold-out residuals.
```

### 6.8 Operations Codex must never perform automatically

Codex must never:

- substitute a model curve for unavailable observations;
- fabricate an uncertainty or quality column;
- infer undocumented scale factors because values “look reasonable”;
- erase or fill real observational gaps to obtain a continuous driver;
- extrapolate beyond the calibrated product domain without an explicit rule;
- smooth LOLA/Diviner data solely for numerical convenience;
- select only observations that agree with the model;
- renormalize each observational case independently;
- change calibration/hold-out membership after final residuals are known;
- silently drop low-quality points without counts and criteria;
- edit immutable raw files;
- hand-edit frozen processed tables;
- weaken test tolerances after failure;
- treat preflight or source inspection as linked numerical validation;
- turn `SKIPPED`, `ERROR`, or `FAIL` into `PASS`.

### 6.9 Codex completion checklist

Before declaring a data-preparation or validation work package complete:

- [ ] Target roadmap item, test README, and acceptance JSON were read.
- [ ] Raw files remain byte-identical.
- [ ] Raw and processed SHA-256 hashes are recorded.
- [ ] Product IDs, versions, source URL/DOI, and calibration level are recorded.
- [ ] Native metadata controlled decoding.
- [ ] Units/time/frame/coordinates are explicit.
- [ ] Filtering criteria and removed/retained counts are recorded.
- [ ] No undocumented smoothing/interpolation/extrapolation/renormalization occurred.
- [ ] No synthetic/mock/assumed observational values were introduced.
- [ ] Processing is deterministic from a committed script.
- [ ] `qa_report.json` passes.
- [ ] Calibration/hold-out split is frozen where applicable.
- [ ] Relevant U-test was run.
- [ ] Target linked I-test was run when the capability exists.
- [ ] Status is correctly classified as PASS/FAIL/ERROR/SKIPPED.
- [ ] Acceptance criteria were not weakened.
- [ ] Data/code README files were updated.
- [ ] Remaining science/model/data limitations are explicitly listed.

## 7. Authoritative data acquisition and processing

This section defines *what to download*, *where to get it*, and *how to turn it into a validation input*. Raw mission files should be retained even when only the processed tables are read by the test runner.

<!-- EXACT_DATA_ACQUISITION_V1 -->

### Exact manual-download matrix

The URLs below point to the **specific dataset, product directory, or direct file** to be acquired. They are intentionally more specific than archive home pages. When a provider does not expose a stable static file URL, the exact dataset application/article is given together with the exact product-selection procedure.

| ID | Required product | Exact data/product URL |
|---|---|---|
| D01 | LRO/LOLA `LDEM_4` global DEM | `https://imbrium.mit.edu/DATA/LOLA_GDR/CYLINDRICAL/IMG/` |
| D02 | LOLA average visibility + PSR, 65° polar, 240 m | `https://imbrium.mit.edu/EXTRAS/ILLUMINATION/IMG/` |
| D03 | Diviner GCP / PRP / time-resolved RDR | `https://pds-geosciences.wustl.edu/lro/urn-nasa-pds-lro_diviner_derived1/data_derived_gcp/` ; `https://pds-geosciences.wustl.edu/lro/urn-nasa-pds-lro_diviner_derived1/data_derived_prp/` ; `https://pds-geosciences.wustl.edu/lro/lro-l-dlre-4-rdr-v1/lrodlr_1001/data/` |
| D04 | PHIDRATES photochemistry | `https://phidrates.space.swri.edu/` |
| D05 | Voronov electron-impact coefficients | `https://www.pa.uky.edu/~verner/dima/col//cfit.dat` |
| D06 | ARTEMIS P1/P2 merged plasma+B+position | `https://cdaweb.gsfc.nasa.gov/cgi-bin/eval2.cgi?dataset=THB_L2_MERGED&index=sp_phys` ; `https://cdaweb.gsfc.nasa.gov/cgi-bin/eval2.cgi?dataset=THC_L2_MERGED&index=sp_phys` |
| D07 | Kaguya UPI-TVIS Na L2A + dark + main-orbiter trajectory | `https://data.darts.isas.jaxa.jp/pub/pds3/sln-e-tvis-5-level2a-v1.0/sln-e-tvis-5-na-level2a-v1.0/` ; `https://data.darts.isas.jaxa.jp/pub/pds3/sln-e-tvis-5-level2a-v1.0/sln-e-tvis-5-dark-level2a-v1.0/` ; `https://data.darts.isas.jaxa.jp/pub/pds3/sln-l-rise-5-traj-main-v1.0/` |
| D08 | LADEE NMS derived + calibrated data | `https://pds.nasa.gov/ds-view/pds/viewCollection.jsp?identifier=urn%3Anasa%3Apds%3Aladee_nms%3Adata_derived` and direct NMS tar files listed below |
| D09 | LRO/LAMP He campaign EDR list and EDR products | `https://academic.oup.com/mnras/article/501/3/4438/6041041` ; `https://ode.rsl.wustl.edu/moon/DataSetExplorer.aspx?datasetid=LRO-L-LAMP-2-EDR-V1.0&instrumenthost=LRO&instrumentid=LAMP` |
| D10 | Kaguya GRS potassium nuclide map v2.0 | `https://data.darts.isas.jaxa.jp/pub/pds3/sln-l-grs-5-nuclide-map-v2.0/` |
| D11 | Current IAU MDC established/full shower tables | `https://www.ta3.sk/IAUC22DB/MDC2022/Etc/streamestablisheddata2026.txt` ; `https://www.ta3.sk/IAUC22DB/MDC2022/Etc/streamfulldata2026.txt` |
| D12 | Apollo 17 LACE NASA-CR-150946 | `https://ntrs.nasa.gov/citations/19760025001` |
| D13 | Kaguya PACE PBF1 v3.0 | `https://data.darts.isas.jaxa.jp/pub/pds3/sln-l-pace-3-pbf1-v3.0/` |
| D14 | Published ground-based Na-tail observations | `https://sirius.bu.edu/aeronomy/Matta%20et%20al%202009.pdf` ; `https://www.sciencedirect.com/science/article/pii/S0019103512001303` |
| D15 | Published LADEE water-event data + supplement | `https://pmc.ncbi.nlm.nih.gov/articles/PMC7306913/` |

**Manual-download rule:** when a directory index is provided, download the native science file **and its matching label/metadata file**, plus the dataset index/catalog/documentation needed to interpret it. Never download only a browse PNG/JPEG when a numeric IMG/TAB/CDF/FITS product exists.

### D01 — LRO/LOLA global topography

#### Exact acquisition — download these files

Primary source: MIT LOLA Data Node, cylindrical numeric IMG products.

**Directory**
`https://imbrium.mit.edu/DATA/LOLA_GDR/CYLINDRICAL/IMG/`

**Download exactly**
- `https://imbrium.mit.edu/DATA/LOLA_GDR/CYLINDRICAL/IMG/LDEM_4.IMG`
- `https://imbrium.mit.edu/DATA/LOLA_GDR/CYLINDRICAL/IMG/LDEM_4.LBL`

Expected archived file sizes for this product:
- `LDEM_4.IMG`: **2,073,600 bytes**
- `LDEM_4.LBL`: about **5,107 bytes**

PDS mirror, if the MIT node is unavailable:
- `https://pds-geosciences.wustl.edu/lro/lro-l-lola-3-rdr-v1/lrolol_1xxx/data/lola_gdr/cylindrical/img/ldem_4.img`
- `https://pds-geosciences.wustl.edu/lro/lro-l-lola-3-rdr-v1/lrolol_1xxx/data/lola_gdr/cylindrical/img/ldem_4.lbl`

**Manual procedure**
1. Create `moon_validation_data/lola/raw/`.
2. Download the `.IMG` and matching `.LBL`; do not rename or edit the raw copies.
3. Check the IMG byte count before processing.
4. Save hashes:
   `sha256sum LDEM_4.IMG LDEM_4.LBL > SHA256SUMS.txt`.
5. Give the raw directory to Codex with the D01/I14 instructions from Section 6.
6. Codex must read the downloaded label rather than trusting hard-coded dimensions/scales.

**Quick local verification**
```bash
cd moon_validation_data/lola/raw
stat -c '%n %s bytes' LDEM_4.IMG LDEM_4.LBL
grep -E 'LINES|LINE_SAMPLES|SAMPLE_TYPE|SAMPLE_BITS|SCALING_FACTOR|OFFSET|MAP_RESOLUTION' LDEM_4.LBL
sha256sum LDEM_4.IMG LDEM_4.LBL
```


**Authoritative source:** LRO/LOLA GDR `LDEM_4` global gridded digital elevation model.

**Exact acquisition URL(s):** `https://imbrium.mit.edu/DATA/LOLA_GDR/CYLINDRICAL/IMG/LDEM_4.IMG` and matching `LDEM_4.LBL`.

**Recommended product:** LRO/LOLA `LDEM_4` global gridded DEM.

**Acquire**
1. Start at the PDS Geosciences LOLA landing page: <https://pds-geosciences.wustl.edu/missions/lro/lola.htm>.
2. Use the global/cylindrical gridded DEM products. For the initial production surface, download the `LDEM_4.IMG` binary and its matching `LDEM_4.LBL` label. The nominal grid is 1440 × 720 at 4 pixels/degree.
3. Keep the `.IMG` and `.LBL` together under `lola/raw/`; do not edit them.
4. Record the access date and compute `sha256sum` for both files.

**Process**
- Read the label first; do not hard-code endianness, dimensions, scale, offset, or map convention without checking the downloaded label.
- For the standard signed-integer `LDEM_4`, decode the image according to the PDS label, convert DN to elevation relative to the 1737.4-km reference sphere, and write strictly increasing longitude/latitude vectors plus `height_m.txt`.
- Reorder the latitude axis only if needed by the AMPS loader; do not resample or smooth the source DEM for the first I14 implementation.
- Extract real grid-cell control points into `control_points.csv`; every control-point value must come directly from the downloaded DEM.
- Preserve `provenance.json`, `SHA256SUMS.txt`, and the processing script.

**Validation use:** production geometry for I14; the DEM itself is *not* an independent validation source for illumination/PSR tests that are computed from the same terrain.

**Minimum roadmap processing requirement:** Use LDEM_4 IMG/LBL for the first global production surface (4 pixels/degree, 1440x720). Preserve labels/checksums; decode heights relative to the 1737.4-km reference sphere; reorder lat/lon only as required by the AMPS loader.

### D02 — LOLA illumination / permanent shadow

#### Exact acquisition — illumination and permanent-shadow reference products

Use the numeric LOLA illumination IMG products, not browse images.

**Directory**
`https://imbrium.mit.edu/EXTRAS/ILLUMINATION/IMG/`

For the baseline I15/I20 reference, download the **65°-to-pole, 240 m** north/south products:

**Average solar visibility**
- `https://imbrium.mit.edu/EXTRAS/ILLUMINATION/IMG/AVGVISIB_65N_240M_201608.IMG`
- `https://imbrium.mit.edu/EXTRAS/ILLUMINATION/IMG/AVGVISIB_65N_240M_201608.LBL`
- `https://imbrium.mit.edu/EXTRAS/ILLUMINATION/IMG/AVGVISIB_65S_240M_201608.IMG`
- `https://imbrium.mit.edu/EXTRAS/ILLUMINATION/IMG/AVGVISIB_65S_240M_201608.LBL`

**Permanent shadow**
- `https://imbrium.mit.edu/EXTRAS/ILLUMINATION/IMG/LPSR_65N_240M_201608.IMG`
- `https://imbrium.mit.edu/EXTRAS/ILLUMINATION/IMG/LPSR_65N_240M_201608.LBL`
- `https://imbrium.mit.edu/EXTRAS/ILLUMINATION/IMG/LPSR_65S_240M_201608.IMG`
- `https://imbrium.mit.edu/EXTRAS/ILLUMINATION/IMG/LPSR_65S_240M_201608.LBL`

Optional newer/high-resolution polar products may be added later, but must be versioned separately rather than silently replacing the baseline.

**Manual procedure**
1. Store these eight files under `moon_validation_data/lola/raw/illumination/`.
2. Hash every file.
3. Keep the north and south products separate.
4. Codex reads scale/offset/projection from each `.LBL` and creates the independent `illumination_reference.csv`.
5. Do **not** derive this reference from the same D01 terrain algorithm being tested; D02 is independent evidence for I15/I20.

**Quick verification**
```bash
cd moon_validation_data/lola/raw/illumination
ls -lh *.IMG *.LBL
sha256sum *.IMG *.LBL > SHA256SUMS.txt
grep -E 'LINES|LINE_SAMPLES|SAMPLE_TYPE|SAMPLE_BITS|SCALING_FACTOR|OFFSET|MAP_' *.LBL
```


**Authoritative source:** LOLA illumination products `AVGVISIB_65N/S_240M_201608` and permanent-shadow products `LPSR_65N/S_240M_201608`.

**Exact acquisition URL(s):** `https://imbrium.mit.edu/EXTRAS/ILLUMINATION/IMG/` — download the exact IMG/LBL pairs listed above.

**Recommended products:** LOLA GDRPSR (permanent shadow) and GDRVIS / NASA PGDA lunar polar illumination products.

**Acquire**
1. Use the ODE LOLA gridded-product documentation: <https://ode.rsl.wustl.edu/moon/pagehelp/Content/Missions_Instruments/LRO/LOLA/GDR/Intro.htm>.
2. For permanent-shadow maps use GDRPSR; for average solar visibility use GDRVIS.
3. NASA PGDA provides north/south polar illumination and permanent-shadow products at several resolutions: <https://pgda.gsfc.nasa.gov/products/69>.
4. Prefer PDS IMG+LBL for exact numeric reproduction or full-resolution GeoTIFF. If using GeoTIFF/COG, do not score overview pixels produced by pyramid averaging.
5. Store north/south, resolution, product version, projection metadata, labels, and checksums.

**Process**
- Apply `value = DN * SCALING_FACTOR + OFFSET` exactly as stated by the product label.
- Reproject only for collocation; retain the native polar stereographic files unchanged.
- Build an independent `illumination_reference.csv` containing latitude, longitude, product value, classification threshold, product ID, resolution, and uncertainty/edge flag.
- Keep D02 separate from the model-computed horizon solution to avoid circular validation.

**Validation use:** independent references for I15 and I20.

**Minimum roadmap processing requirement:** Use published illumination/PSR products only as independent classification references. Keep them separate from LOLA DEM used by the model, avoiding circular validation.

### D03 — LRO/Diviner temperature

#### Exact acquisition — Diviner products

##### A. Global temperature climatology for I16

**Exact GCP directory**
`https://pds-geosciences.wustl.edu/lro/urn-nasa-pds-lro_diviner_derived1/data_derived_gcp/`

The collection is split into 10° latitude bands. Download the `.tab` **and matching `.xml`** for every latitude band used by the validation. For a full global package, download all 18 band pairs:
`global_cumul_avg_cyl_00n10n_002`, ..., `global_cumul_avg_cyl_80n90n_002`, and the corresponding southern bands through `global_cumul_avg_cyl_90s80s_002`.

Also download:
- `https://pds-geosciences.wustl.edu/lro/urn-nasa-pds-lro_diviner_derived1/data_derived_gcp/collection_data_derived_gcp_inventory.csv`
- `https://pds-geosciences.wustl.edu/lro/urn-nasa-pds-lro_diviner_derived1/data_derived_gcp/collection_data_derived_gcp.xml`

Example exact pair:
- `https://pds-geosciences.wustl.edu/lro/urn-nasa-pds-lro_diviner_derived1/data_derived_gcp/global_cumul_avg_cyl_00n10n_002.tab`
- `https://pds-geosciences.wustl.edu/lro/urn-nasa-pds-lro_diviner_derived1/data_derived_gcp/global_cumul_avg_cyl_00n10n_002.xml`

##### B. Polar-resource / cold-region reference for I20

**Exact PRP directory**
`https://pds-geosciences.wustl.edu/lro/urn-nasa-pds-lro_diviner_derived1/data_derived_prp/`

Download:
- `https://pds-geosciences.wustl.edu/lro/urn-nasa-pds-lro_diviner_derived1/data_derived_prp/dlre_prp_north.tab`
- `https://pds-geosciences.wustl.edu/lro/urn-nasa-pds-lro_diviner_derived1/data_derived_prp/dlre_prp_north.xml`
- `https://pds-geosciences.wustl.edu/lro/urn-nasa-pds-lro_diviner_derived1/data_derived_prp/dlre_prp_south.tab`
- `https://pds-geosciences.wustl.edu/lro/urn-nasa-pds-lro_diviner_derived1/data_derived_prp/dlre_prp_south.xml`

##### C. Time-resolved temperature/radiance for I17

For 2009–2016 time-resolved Diviner RDRs use:
`https://pds-geosciences.wustl.edu/lro/lro-l-dlre-4-rdr-v1/lrodlr_1001/data/`

Direct year selectors include:
- `https://pds-geosciences.wustl.edu/lro/lro-l-dlre-4-rdr-v1/lrodlr_1001/data/2013/`
- `https://pds-geosciences.wustl.edu/lro/lro-l-dlre-4-rdr-v1/lrodlr_1001/data/2014/`

Inside a year choose the `YYYYMM/` directory, then the date/product covering the exact validation interval. Download every selected science table/product together with its detached label/metadata.

**Manual procedure**
1. Put GCP files under `diviner/raw/gcp/`, PRP under `diviner/raw/prp/`, and RDR files under `diviner/raw/rdr/YYYYMMDD/`.
2. Never use the GCP climatology as a substitute for a transient eclipse/sunrise history.
3. Hash the raw products before Codex processing.
4. Codex should prepare `temperature_table.csv` from GCP/PRP and `eclipse_thermal_history.csv` from RDR records.

**Verification**
```bash
find moon_validation_data/diviner/raw -type f -maxdepth 4 -print
find moon_validation_data/diviner/raw -type f -exec sha256sum {} \; > moon_validation_data/diviner/raw/SHA256SUMS.txt
```


**Authoritative source:** LRO/Diviner GCP (global climatology), PRP (polar resource), and Level-4 RDR (time-resolved) products.

**Exact acquisition URL(s):** GCP: `https://pds-geosciences.wustl.edu/lro/urn-nasa-pds-lro_diviner_derived1/data_derived_gcp/`; PRP: `https://pds-geosciences.wustl.edu/lro/urn-nasa-pds-lro_diviner_derived1/data_derived_prp/`; RDR 2009–2016: `https://pds-geosciences.wustl.edu/lro/lro-l-dlre-4-rdr-v1/lrodlr_1001/data/`.

**Recommended products:** Diviner GCP/PCP/PRP for climatology/poles and RDR for time-resolved thermal histories.

**Acquire**
1. Start at <https://pds-geosciences.wustl.edu/missions/lro/diviner.htm>.
2. For global temperature climatology download the Derived Data Products / Global Cumulative Products (GCP).
3. For polar/cold-trap work download Polar Cumulative Products (PCP) and, where relevant, Polar Resource Products (PRP).
4. For sunrise/sunset/eclipses or arbitrary epoch histories use Reduced Data Records (RDR), not a climatological map.
5. Preserve PDS4 XML labels and any quality/geometry metadata together with the data files.

**Process**
- Convert temperatures to kelvin only when required by the label; do not alter valid calibrated values.
- Apply documented quality flags and fill-value handling.
- Convert longitude/local-time conventions into one documented AMPS convention.
- Produce `temperature_table.csv` for I16 and `eclipse_thermal_history.csv` for I17, including time, lat/lon, local solar time, observed temperature, quality flag, and uncertainty/dispersion when available.
- Partition training/calibration and hold-out sites before tuning the thermal model.

**Validation use:** I16, I17, I20, and later Ar/H2O surface-physics validation.

**Minimum roadmap processing requirement:** Use Global Cumulative Products for local-time/latitude/longitude thermal climatology; use polar/resource products where appropriate. Preserve product XML/TAB labels; convert to SI K and documented coordinates; apply quality flags.

### D04 — Photoionization / photodissociation

#### Exact acquisition — PHIDRATES

PHIDRATES is itself the interactive data application; it does not expose a single stable static URL for every species/product combination.

**Exact application**
`https://phidrates.space.swri.edu/`

For this roadmap acquire the real PHIDRATES values for:
- atomic neutrals: **He, Na, Ne, Ar**;
- molecular species: **H2O** and **OH**;
- every photodissociation/photoionization branch actually enabled by the production model.

**Manual procedure**
1. Open the exact application URL above.
2. Select the species.
3. Select the appropriate photoionization/photodissociation cross-section or photorate view.
4. Select and record the Solar Radiation Field / activity assumption used for the rate.
5. Use the application's download/export facility for the numeric table when provided.
6. If the application presents a numeric text table without a separate download button, save that numeric table verbatim and record the species, branch, solar spectrum, and page state in `provenance.json`.
7. Save the method/reference information associated with each species.
8. Do not manually “repair” branch sums or fill missing channels.

Because PHIDRATES is an application rather than a conventional directory tree, Codex must record the exact species/channel/solar-spectrum selection in provenance for reproducibility.


**Authoritative source:** PHIDRATES photochemical cross sections/rates for the exact species/reaction branches enabled in the Moon model.

**Exact acquisition URL(s):** `https://phidrates.space.swri.edu/` — select the exact species, reaction branch, and solar-radiation field and export/save the numerical values.

**Source:** PHIDRATES, SwRI: <https://phidrates.space.swri.edu/>.

**Acquire**
1. Open the PHIDRATES species page for each parent species used by the Moon model.
2. Record the exact solar spectrum/activity case, wavelength/rate convention, branch definitions, source references, and access date.
3. Export/download machine-readable values where the site provides them. If values must be transcribed, preserve the source page/PDF/screenshot and double-check the transcription.
4. Never create an undocumented channel to make branching sum to unity.

**Process**
- Store raw/reference values by species in `photochemistry/raw/`.
- Generate a versioned table with parent, product(s), 1-AU rate, branch fraction, spectrum identifier, units, source reference, and validity notes.
- Verify branch sums and unit conversions independently in the processing script.
- Scale heliocentric rates in production only through the documented physical law (normally inverse-square when appropriate).

**Validation use:** U14/U15/U24 and linked I21/I22/I31.

**Minimum roadmap processing requirement:** Export/reference rates and branching for the required species and selected solar spectrum. Preserve retrieval date, species page, spectrum choice, and source references. Never interpolate across undocumented channels.

### D05 — Electron-impact ionization

#### Exact acquisition — Voronov electron-impact coefficients

Use the original Voronov/Verner tabulation and reference implementation.

**Download**
- coefficient table: `https://www.pa.uky.edu/~verner/dima/col//cfit.dat`
- table description: `https://www.pa.uky.edu/~verner/dima/col//cfit.txt`
- reference Fortran: `https://www.pa.uky.edu/~verner/dima/col//cfit.f`
- dataset description page: `https://www.pa.uky.edu/~verner/col.html`

**Manual procedure**
1. Save all three source files under `moon_validation_data/electron_impact/raw/`.
2. Hash them before parsing.
3. Codex must parse the actual neutral-atom entries for He, Ne, Na, and Ar from `cfit.dat`.
4. Use `cfit.txt` to interpret every coefficient and valid temperature range.
5. Use `cfit.f` as an independent implementation cross-check.
6. Do not transcribe coefficients into code by hand without an automated table comparison.

**Verification**
```bash
cd moon_validation_data/electron_impact/raw
wc -l cfit.dat cfit.txt cfit.f
sha256sum cfit.dat cfit.txt cfit.f > SHA256SUMS.txt
```


**Authoritative source:** Voronov (1997) electron-impact ionization fit coefficients as hosted by D. Verner.

**Exact acquisition URL(s):** `https://www.pa.uky.edu/~verner/dima/col//cfit.dat` with `cfit.txt` and `cfit.f` from the same directory.

**Source:** Voronov/Verner collisional-ionization tables: <https://www.pa.uky.edu/~verner/col.html>.

**Acquire**
1. From the page, download the ASCII table for the Voronov (1997) practical fit formula and the accompanying description.
2. Preserve the original ASCII coefficients and bibliography under `electron_impact/raw/`.
3. Record the atomic/ion stage used for He, Ne, Na, and Ar.

**Process**
- Implement the published rate formula exactly once in the shared production/unit kernel.
- Generate `rate_vs_Te.csv` over the electron-temperature range used by the lunar campaign.
- Unit-test the implementation at the published table nodes and several intermediate temperatures.
- If a later campaign uses non-Maxwellian electron distributions, archive the cross sections/distribution source separately and do not silently reuse the Maxwellian fit.

**Validation use:** U16 and linked I23.

**Minimum roadmap processing requirement:** Use published fit coefficients or tabulated rate coefficients. Generate species-specific rate tables versus electron temperature and compare production integration against those tables.

### D06 — Solar-wind / magnetospheric plasma

#### Exact acquisition — ARTEMIS P1/P2

Use the dataset-specific CDAWeb pages, not a generic CDAWeb search.

**ARTEMIS P1 / THEMIS-B merged L2**
`https://cdaweb.gsfc.nasa.gov/cgi-bin/eval2.cgi?dataset=THB_L2_MERGED&index=sp_phys`

**ARTEMIS P2 / THEMIS-C merged L2**
`https://cdaweb.gsfc.nasa.gov/cgi-bin/eval2.cgi?dataset=THC_L2_MERGED&index=sp_phys`

The merged products provide ESA plasma moments, FGM magnetic field, and GSE/SSE position at standard merged cadence.

**Manual CDAWeb procedure**
1. Open the P1 or P2 URL above.
2. Enter the exact **Start time** and **Stop time** for the validation interval.
3. For reproducible raw input, select **Download original files** / original CDF rather than relying only on a plotted/binned product.
4. If a compact human-readable check is useful, additionally select **List Data (ASCII/CSV)** for the same interval and variables.
5. Select at least:
   - FGM/FGS magnetic field vector;
   - ESA good ion density;
   - ion thermal velocity/temperature quantity needed by the driver;
   - plasma flow speed / ion velocity vector;
   - spacecraft position in SSE;
   - any data-quality/mode variables used in filtering.
6. Save the original CDF(s) under `drivers/raw/artemis/<probe>/<YYYYMMDD>/`.
7. Save the CDAWeb dataset metadata/master information with the files.
8. Hash before Codex processing.

**Important for He source validation:** `THB_L2_MERGED`/`THC_L2_MERGED` standard moments do not by themselves justify a fixed alpha/proton fraction. If I26/I35 requires alpha flux, download an alpha-resolved ARTEMIS ESA product or reproduce the alpha-flux extraction used by the cited He paper. Codex must not assume a constant He++ fraction.


**Authoritative source:** ARTEMIS P1/P2 Level-2 merged ESA plasma + FGM magnetic field + position products.

**Exact acquisition URL(s):** P1: `https://cdaweb.gsfc.nasa.gov/cgi-bin/eval2.cgi?dataset=THB_L2_MERGED&index=sp_phys`; P2: `https://cdaweb.gsfc.nasa.gov/cgi-bin/eval2.cgi?dataset=THC_L2_MERGED&index=sp_phys`.

**Source:** NASA SPDF/CDAWeb ARTEMIS.

**Acquire**
1. Open CDAWeb: <https://cdaweb.gsfc.nasa.gov/>.
2. For ARTEMIS P1 (THEMIS-B), use `THB_L2_MERGED` for a convenient combined product containing ESA plasma moments, FGM magnetic field, and geocentric/selenocentric position. For ARTEMIS P2 (THEMIS-C), use `THC_L2_MERGED`.
3. When higher cadence or instrument-specific handling is required, use `THB_L2_ESA` / `THC_L2_ESA` and `THB_L2_FGM` / `THC_L2_FGM` separately.
4. In CDAWeb select the event interval, choose `List Data (ASCII/CSV)` for compact validation inputs or download the original CDFs for a reproducible archive.
5. Always retain ESA data-quality flags and document whether FULL, REDUCED, or another mode is used.

**Process**
- Convert density, velocity, temperature/thermal velocity, magnetic field, and SSE position to SI units and a single epoch/time standard.
- Preserve gaps; do not interpolate across gaps longer than a declared maximum.
- Interpolate only at AMPS update times and test interpolation against exact input samples.
- If alpha-particle flux is needed, use an alpha-resolved archived/published product or reproduce a documented extraction from distribution data. Do **not** assume a fixed alpha/proton ratio.

**Validation use:** U18; linked I25; driver input for I26, I35, and I40.

**Minimum roadmap processing requirement:** Download time-tagged ion moments, velocities, magnetic field and selenocentric position; retain quality flags. If alpha flux is required, use an archived/published alpha-specific product or reproduce a published extraction from distributions—do not assume a constant alpha/proton ratio.

### D07 — Kaguya/SELENE sodium imaging

#### Exact acquisition — Kaguya UPI-TVIS sodium images

##### Na Level 2A

Dataset page:
`https://darts.isas.jaxa.jp/datasets/darts:sln-e-tvis-5-na-level2a-v1.0`

Direct distribution:
`https://data.darts.isas.jaxa.jp/pub/pds3/sln-e-tvis-5-level2a-v1.0/sln-e-tvis-5-na-level2a-v1.0/`

Products are organized by volume directory `NNNN/`. Each volume contains:
`data/tvis_{YYMMDDhhmmss}_open.img` plus matching `.lbl`, and `index/index.tab`.

##### Dark Level 2A

Dataset page:
`https://darts.isas.jaxa.jp/datasets/darts:sln-e-tvis-5-dark-level2a-v1.0`

Direct distribution:
`https://data.darts.isas.jaxa.jp/pub/pds3/sln-e-tvis-5-level2a-v1.0/sln-e-tvis-5-dark-level2a-v1.0/`

##### Main-orbiter trajectory

Dataset page:
`https://darts.isas.jaxa.jp/datasets/darts:sln-l-rise-5-traj-main-v1.0`

Direct distribution:
`https://data.darts.isas.jaxa.jp/pub/pds3/sln-l-rise-5-traj-main-v1.0/`

##### Product-format documentation

`https://darts.isas.jaxa.jp/app/pdap/selene/help/en/UPI_Format_en_V01.pdf`

**Manual procedure**
1. Choose candidate observation epochs **before** model comparison.
2. In the Na distribution, open each volume `index/index.tab` and locate products covering those epochs.
3. Download each selected Na `.img` and matching `.lbl`; also download that volume's `index.tab`, `index.lbl`, `aareadme.txt`, catalog, and relevant document files.
4. Download the corresponding dark product(s) required by the documented reduction.
5. Download the RISE main-orbiter trajectory file covering the same epoch.
6. Preserve the original volume/epoch directory structure locally.
7. Hash all raw files.
8. Codex then creates the calibrated/geometry-aware comparison table for I32.

Do not select “good-looking” Na images after examining the AMPS prediction; freeze the epoch list first.


**Authoritative source:** Kaguya/SELENE UPI-TVIS Na Level-2A + dark Level-2A + RISE main-orbiter trajectory.

**Exact acquisition URL(s):** Na: `https://data.darts.isas.jaxa.jp/pub/pds3/sln-e-tvis-5-level2a-v1.0/sln-e-tvis-5-na-level2a-v1.0/`; dark: `https://data.darts.isas.jaxa.jp/pub/pds3/sln-e-tvis-5-level2a-v1.0/sln-e-tvis-5-dark-level2a-v1.0/`; trajectory: `https://data.darts.isas.jaxa.jp/pub/pds3/sln-l-rise-5-traj-main-v1.0/`.

**Source:** DARTS Kaguya/SELENE UPI-TVIS Na Level 2A, dataset `darts:sln-e-tvis-5-na-level2a-v1.0`.

**Acquire**
1. Dataset page: <https://darts.isas.jaxa.jp/datasets/darts:sln-e-tvis-5-na-level2a-v1.0>.
2. The distribution tree contains per-volume directories with `data/*.img` and matching `*.lbl`, plus `index/index.tab`, catalog files, and documentation.
3. Download the Na image/label pairs for selected epochs, the index/catalog/documentation, and the matching TVIS dark-image dataset from the Kaguya archive.
4. Download/retain the spacecraft/geometry information needed to reconstruct the line of sight; use the archive product documentation rather than inferring the camera orientation.

**Process**
- Read every image according to its PDS3 label.
- Preserve the Level-2A calibrated Rayleigh units and image mask/quality information.
- Apply dark/background correction exactly as supported by the archive documentation/published TVIS processing.
- Generate an observation table with image ID, UTC, spacecraft state, pixel or binned LOS, tangent altitude/location, observed Rayleigh brightness, uncertainty/quality, and `calibration_or_holdout`.
- Freeze the calibration/hold-out image list before final I32 scoring.

**Validation use:** I32; geometry also exercises U26.

**Minimum roadmap processing requirement:** Download IMG/LBL plus dark/background products and relevant SPICE/geometry. Convert calibrated Rayleigh brightness and image geometry into line-of-sight comparison tables; retain original image products.

### D08 — LADEE NMS He/Ne/Ar and calibrated/derived data

#### Exact acquisition — LADEE NMS

##### Derived He/Ne/Ar collection

Exact PDS collection:
`https://pds.nasa.gov/ds-view/pds/viewCollection.jsp?identifier=urn%3Anasa%3Apds%3Aladee_nms%3Adata_derived`

Collection identifier:
`urn:nasa:pds:ladee_nms:data_derived::2.0`

Direct packaged download:
`https://atmos.nmsu.edu/PDS/data/PDS4/LADEE/nms_data_derived.tar.gz`

##### Calibrated NMS data

Direct packaged download:
`https://atmos.nmsu.edu/PDS/data/PDS4/LADEE/nms_data_calibrated.tar.gz`

##### Calibration and documentation packages

- `https://atmos.nmsu.edu/PDS/data/PDS4/LADEE/nms_calibration.tar.gz`
- `https://atmos.nmsu.edu/PDS/data/PDS4/LADEE/nms_document.tar.gz`

##### LADEE SPICE geometry

Meta-kernel directory:
`https://naif.jpl.nasa.gov/pub/naif/LADEE/kernels/mk/`

Use the latest archive meta-kernel compatible with the mission data; the directory includes `ladee_v03.tm`.

Kernel root:
`https://naif.jpl.nasa.gov/pub/naif/LADEE/kernels/`

**Manual procedure**
1. Download `nms_data_derived.tar.gz` for I34/I36/I37 noble-gas validation.
2. Download `nms_data_calibrated.tar.gz` plus calibration/document packages for audit and for I39 water-event reconstruction.
3. Download the LADEE meta-kernel and all kernel files it references, or download the complete required kernel subdirectories.
4. Preserve the tarballs unchanged and hash them **before extraction**.
5. Extract into versioned subdirectories while retaining all PDS4 XML labels and collection inventories.
6. Codex must derive geometry from SPICE and product metadata; do not infer spacecraft longitude/local time from plot axes.

**Verification**
```bash
sha256sum nms_data_derived.tar.gz nms_data_calibrated.tar.gz nms_calibration.tar.gz nms_document.tar.gz
tar -tzf nms_data_derived.tar.gz | head
tar -tzf nms_data_calibrated.tar.gz | head
```


**Authoritative source:** LADEE NMS PDS4 derived He/Ne/Ar collection plus calibrated NMS products and LADEE SPICE.

**Exact acquisition URL(s):** Derived collection: `https://pds.nasa.gov/ds-view/pds/viewCollection.jsp?identifier=urn%3Anasa%3Apds%3Aladee_nms%3Adata_derived`; packaged downloads are listed in the exact-acquisition block above.

**Source:** PDS4 collection `urn:nasa:pds:ladee_nms:data_derived::2.0` (DOI `10.17189/1408897`) plus the exact calibrated/calibration/document tarballs listed above.

**Acquire**
1. Download the exact derived/calibrated packages from the URLs above; use the PDS derived-collection page to verify collection identifier/version and inventories.
2. Download the NMS Derived Data Collection for He/Ne/Ar and the calibrated collection/geometry/context files needed to audit the derived values.
3. Preserve the collection inventory CSV/XML and product XML labels.
4. Use product versions from the archive; do not mix old local copies with a newer PDS release without recording the version.

**Process**
- Apply instrument-mode and quality filtering from the product documentation.
- Compute/collocate spacecraft altitude, latitude, longitude, local solar time, and model epoch.
- Keep observations in native physical units and preserve reported uncertainties/quality.
- Build separate calibration and hold-out intervals for each species.
- For water-event work, use archived NMS measurements together with the published D15 event list rather than manufacturing event windows.

**Validation use:** I34, I36, I37, I38, I39, and cross-species I41.

**Minimum roadmap processing requirement:** Download derived He/Ne/Ar plus calibrated products/geometry needed for independent checks. Filter by quality flags and instrument mode; compute spacecraft altitude, local solar time, latitude/longitude, and collocate AMPS without fitting the hold-out samples.

### D09 — LRO/LAMP helium

#### Exact acquisition — LRO/LAMP helium campaign

Use the Grava et al. 2021 campaign because its supplement explicitly identifies the LAMP EDR files and the intervals used for the He analysis.

**Exact paper / observation-selection source**
`https://academic.oup.com/mnras/article/501/3/4438/6041041`

DOI:
`https://doi.org/10.1093/mnras/staa3884`

On that article page, under **SUPPORTING INFORMATION**, download **`Tab_full.dat`**. This table lists the LAMP EDR filenames, UT interval of interest, LRO maneuver type, He g-factor, and integrated ARTEMIS alpha-particle flux used in the published analysis.

**Exact LAMP EDR dataset**
Dataset ID: `LRO-L-LAMP-2-EDR-V1.0`

Dataset browser:
`https://ode.rsl.wustl.edu/moon/DataSetExplorer.aspx?datasetid=LRO-L-LAMP-2-EDR-V1.0&instrumenthost=LRO&instrumentid=LAMP`

**Manual procedure**
1. Download `Tab_full.dat` first and freeze it as the campaign product list.
2. For every EDR filename in `Tab_full.dat`, open the exact LAMP EDR dataset browser above and retrieve the corresponding `.FIT` and detached `.LBL`.
3. Do not substitute a different LAMP orbit because it is easier to locate.
4. Save all EDRs under `lamp/raw/edr/` and the supplement under `lamp/raw/publication/`.
5. For independent driver reproduction, also download the corresponding ARTEMIS interval via D06; use the published integrated alpha flux in `Tab_full.dat` as a cross-check, not as a reason to assume alpha/proton ratio.
6. Hash all files before processing.
7. Codex should reproduce at least one published He column-density/source-rate result before I35 is considered data-ready.


**Authoritative source:** LRO/LAMP EDR products selected by Grava et al. (2021) `Tab_full.dat`, with ARTEMIS forcing.

**Exact acquisition URL(s):** Campaign/file list: `https://academic.oup.com/mnras/article/501/3/4438/6041041`; LAMP EDR dataset ID `LRO-L-LAMP-2-EDR-V1.0`: `https://ode.rsl.wustl.edu/moon/DataSetExplorer.aspx?datasetid=LRO-L-LAMP-2-EDR-V1.0&instrumenthost=LRO&instrumentid=LAMP`.

**Source:** LRO/LAMP EDR dataset `LRO-L-LAMP-2-EDR-V1.0`, with the exact campaign EDR filenames/UT windows frozen by Grava et al. (2021) `Tab_full.dat`.

**Acquire**
1. Download `Tab_full.dat` from the exact Grava et al. article page above and use its listed EDR filenames/UT windows as the frozen observation list.
2. Download the calibrated spectral/brightness products, labels, geometry files, and the relevant documentation.
3. Record the exact LAMP product IDs and the publication used to define the He retrieval/selection.
4. Obtain ARTEMIS forcing for the same epochs from D06.

**Process**
- Reproduce the published background subtraction, spectral window, and He retrieval as closely as the archive products permit.
- Convert the result into a time-tagged comparison table with geometry, brightness/column quantity, uncertainty, and quality.
- Match ARTEMIS forcing to LAMP epochs with an explicitly documented time window and propagation assumption.

**Validation use:** I35.

**Minimum roadmap processing requirement:** Use archived spectral/brightness products when available; document exact product IDs. Coordinate with ARTEMIS only after matching epochs and propagating timing/geometry uncertainties.

### D10 — Lunar potassium map for radiogenic Ar proxy

#### Exact acquisition — Kaguya GRS potassium abundance map

Use the **Nuclide Map v2.0** potassium product, not the older preliminary gamma-ray count-rate map.

Dataset page:
`https://darts.isas.jaxa.jp/datasets/darts:sln-l-grs-5-nuclide-map-v2.0`

Direct distribution:
`https://data.darts.isas.jaxa.jp/pub/pds3/sln-l-grs-5-nuclide-map-v2.0/`

Download the **unsmoothed potassium map** and label:
- `https://data.darts.isas.jaxa.jp/pub/pds3/sln-l-grs-5-nuclide-map-v2.0/data/GRS_NMAP_K_SPA_090210_090527.img`
- `https://data.darts.isas.jaxa.jp/pub/pds3/sln-l-grs-5-nuclide-map-v2.0/data/GRS_NMAP_K_SPA_090210_090527.lbl`

Also download:
- `aareadme.txt`;
- `index/index.tab` and `index/index.lbl`;
- relevant `document/` and `catalog/` files.

A smoothed K product also exists, but it must be treated as a separate sensitivity product and not silently replace the unsmoothed baseline.

**Manual procedure**
1. Store the unsmoothed `.img/.lbl` under `grs/raw/`.
2. Hash all files.
3. Codex decodes the map from the native label and creates a non-negative K spatial proxy.
4. Normalize only the *spatial source weights* when constructing the Ar source pattern.
5. Do not derive the absolute radiogenic Ar production rate from K abundance alone; absolute normalization remains an independent physical parameter/constraint.


**Authoritative source:** Kaguya/SELENE GRS Nuclide Map v2.0 potassium abundance product `darts:sln-l-grs-5-nuclide-map-v2.0`.

**Exact acquisition URL(s):** `https://data.darts.isas.jaxa.jp/pub/pds3/sln-l-grs-5-nuclide-map-v2.0/data/GRS_NMAP_K_SPA_090210_090527.img` and matching `.lbl`.

**Source:** DARTS Kaguya GRS Nuclide Map v2.0, dataset `darts:sln-l-grs-5-nuclide-map-v2.0`.

**Acquire**
1. Dataset page: <https://darts.isas.jaxa.jp/datasets/darts:sln-l-grs-5-nuclide-map-v2.0>.
2. Download the unsmoothed `GRS_NMAP_K_SPA_090210_090527.img` and matching `.lbl` from the direct distribution specified above.
3. Download the dataset index, catalog, and document files so the abundance units/projection are preserved.
4. Keep the smoothed map only as an explicitly labeled sensitivity product.

**Process**
- Decode the K abundance map according to its label.
- Convert grid coordinates to the lunar convention used by AMPS.
- Create a non-negative spatial source-weight map and normalize only the integrated spatial weight to 1.0.
- Do **not** infer the absolute 40Ar source rate from K abundance alone; absolute normalization must come from independent Ar constraints or a cited radiogenic production model.
- Keep the unsmoothed product as the baseline source shape.

**Validation use:** U21 and linked I28; later I37/I38.

**Minimum roadmap processing requirement:** Extract K-map layer, product geometry and uncertainties/limitations. Normalize the map to unit integrated source weight, then determine the absolute 40Ar production scale from independent Ar data—not from the K map itself.

### D11 — Meteor shower forcing

#### Exact acquisition — IAU Meteor Data Center

For current/future campaign forcing use the current machine-readable MDC text files.

**Established showers**
`https://www.ta3.sk/IAUC22DB/MDC2022/Etc/streamestablisheddata2026.txt`

**Full shower solutions**
`https://www.ta3.sk/IAUC22DB/MDC2022/Etc/streamfulldata2026.txt`

These files contain the format/header definitions in the file itself and are versioned by their internal “Last update” line.

**Manual procedure**
1. Download both text files unchanged to `meteoroids/raw/iau_mdc/`.
2. Record the internal `Last update` line in provenance.
3. Hash both files.
4. Codex parses the pipe-delimited/fixed-format fields and preserves shower status, solution number, radiant/orbit parameters, activity window, velocity, and references.
5. For I30, derive the lunar encounter from the published orbit/geometry, rather than treating an Earth peak time as automatically equal to the lunar peak.

**Important historical-reproduction rule for I39:** Benna et al. used an MDC snapshot dated **18 January 2018**. The current 2026 MDC files are not equivalent to that historical snapshot. To reproduce the published LADEE water-event association, use the event/stream assignments in Benna et al. **Supplementary Table 1 (D15)** as the frozen historical reference. Use the 2026 MDC tables only for new/current forward campaigns or sensitivity checks.


**Authoritative source:** IAU Meteor Data Center machine-readable shower-solution tables.

**Exact acquisition URL(s):** Established: `https://www.ta3.sk/IAUC22DB/MDC2022/Etc/streamestablisheddata2026.txt`; full solutions: `https://www.ta3.sk/IAUC22DB/MDC2022/Etc/streamfulldata2026.txt`.

**Source:** Current MDC machine-readable shower tables at the exact `.txt` URLs given above; use D15 Supplementary Table 1 for the historical 2018 shower assignment used by Benna et al.

**Acquire**
1. Download the exact `streamestablisheddata2026.txt` and `streamfulldata2026.txt` files given above for current/new campaign work.
2. Record whether each shower is established, working, or otherwise classified by the IAU MDC.
3. Preserve the source table/version/access date.
4. For each LADEE/Kaguya epoch, compute lunar encounter geometry from the published shower orbit rather than assuming the Earth encounter time applies unchanged at the Moon.

**Process**
- Build `meteoroid_driver.csv` with UTC, stream ID, confidence/status, encounter speed, radiant/impact geometry, relative forcing, and uncertainty.
- Combine shower geometry with an independently cited/background meteoroid flux model and impact-vaporization yield.
- Keep observed water-event labels separate from the driver calculation to avoid tuning the driver to the answer.

**Validation use:** U23, I30, and I39.

**Minimum roadmap processing requirement:** Download shower list/orbital parameters for relevant epochs. Compute lunar encounter geometry from published shower parameters; combine with independently specified/background meteoroid flux model. Record shower status (established/working) and uncertainty.

### D12 — Apollo 17 LACE Ar

#### Exact acquisition — Apollo 17 LACE

NASA report record:
`https://ntrs.nasa.gov/citations/19760025001`

Direct PDF endpoint:
`https://ntrs.nasa.gov/api/citations/19760025001/downloads/19760025001.pdf`

Document ID: `19760025001`  
Report: `NASA-CR-150946`, *Lunar atmospheric composition experiment*.

If the direct PDF endpoint is blocked by the browser/network, open the NTRS record above and click **19760025001.pdf** under **Available Downloads**.

**Manual procedure**
1. Save the PDF unchanged under `lace/raw/`.
2. Hash it.
3. Search the report and its cited LACE papers for tabulated Ar measurements before digitizing plots.
4. If no machine-readable/table value exists for the needed I38 quantity, digitize only the specific published figure and record report ID, page, figure, axis scale, digitization uncertainty, and independent spot checks.
5. Do not invent a global Ar map from this single Apollo-site measurement.


**Authoritative source:** Apollo 17 LACE report NASA-CR-150946 (`NTRS 19760025001`) and its cited LACE data publications.

**Exact acquisition URL(s):** `https://ntrs.nasa.gov/citations/19760025001`; direct PDF endpoint `https://ntrs.nasa.gov/api/citations/19760025001/downloads/19760025001.pdf`.

**Source:** Apollo 17 LACE public reports, especially NASA-CR-150946: <https://ntrs.nasa.gov/citations/19760025001>.

**Acquire**
1. Download the PDF from NTRS and preserve it unchanged.
2. Search the report and cited LACE publications for machine-readable/tabulated Ar and He measurements before digitizing plots.
3. If only figures are available, use a traceable digitization tool and record report ID, page, figure, axis limits/scales, and the digitized point set.

**Process**
- Digitize independently twice or perform an independent manual spot check.
- Include a digitization uncertainty term in addition to the measurement/model uncertainty.
- Convert mission time/lunation/local-time descriptions into UTC/local solar time with documented assumptions.
- Keep Apollo-site geometry explicit; this is a single-site surface measurement, not a global profile.

**Validation use:** I38.

**Minimum roadmap processing requirement:** Prefer machine-readable numerical tables if located. If only plots are available, digitize the original report/paper with a documented tool, record page/figure and axis calibration, digitize twice or cross-check manually, and carry digitization uncertainty.

### D13 — Kaguya PACE pickup ions

#### Exact acquisition — Kaguya PACE PBF1

Dataset page:
`https://darts.isas.jaxa.jp/en/datasets/sln-l-pace-3-pbf1-v3.0`

Direct distribution:
`https://data.darts.isas.jaxa.jp/pub/pds3/sln-l-pace-3-pbf1-v3.0/`

Each observation date has its own `YYYYMMDD/` directory. Download all required sensor files and labels for a frozen pickup-ion event date:

```text
IPACE_PBF1_YYMMDD_ESA1_V003.dat.gz + .lbl
IPACE_PBF1_YYMMDD_ESA2_V003.dat.gz + .lbl
IPACE_PBF1_YYMMDD_IEA_V003.dat.gz  + .lbl
IPACE_PBF1_YYMMDD_IMA_V003.dat.gz  + .lbl
```

Example real date directory:
`https://darts.isas.jaxa.jp/pub/pds3/sln-l-pace-3-pbf1-v3.0/20090426/data/`

Also download from the selected date directory:
- `index/index.tab` and `.lbl`;
- catalog/document files;
- `software/read_pbf_v2.c` and its header if present.

For spacecraft geometry use D07's RISE trajectory product covering the same interval; acquire magnetic/plasma context required by the selected published event before scoring.

**Manual procedure**
1. Choose the event interval from the published pickup-ion study **before** looking at AMPS.
2. Download the corresponding date directory's required PBF1 products and metadata.
3. Hash compressed raw files before decompression.
4. Codex reproduces the publication's species/energy/angle selection; it must not choose cuts after inspecting model output.


**Authoritative source:** Kaguya/SELENE PACE PBF1 v3.0 processed high-resolution electron/ion spectra.

**Exact acquisition URL(s):** `https://data.darts.isas.jaxa.jp/pub/pds3/sln-l-pace-3-pbf1-v3.0/`; select the frozen `YYYYMMDD/` event directory.

**Source:** DARTS Kaguya PACE PBF1, dataset `darts:sln-l-pace-3-pbf1-v3.0`.

**Acquire**
1. Dataset page: <https://darts.isas.jaxa.jp/en/datasets/sln-l-pace-3-pbf1-v3.0>.
2. Select the event date directory. Each date contains PBF1 products for ESA1, ESA2, IMA, and IEA as `*.dat.gz` plus matching labels.
3. Download the sensor products required by the cited pickup-ion paper/event, the labels, index/catalog/documentation, and any included reader software.
4. Obtain spacecraft trajectory and magnetic/plasma context for the same interval.

**Process**
- Reproduce the energy/mass/angle selection criteria used in the publication; do not create ad hoc species cuts after inspecting the AMPS result.
- Convert the observation into the same energy/angle/mass-bin convention used by the model observation operator.
- Preserve absolute flux only where the selected instrument product/calibration supports it.

**Validation use:** I40.

**Minimum roadmap processing requirement:** Download the relevant PBF1 intervals, trajectory/SPICE and magnetic/plasma context. Reproduce energy/mass-angle selections from the cited pickup-ion papers rather than creating ad hoc event cuts.

### D14 — Ground-based sodium tail

#### Exact acquisition — ground-based lunar sodium tail

There is no stable public machine-readable archive for the complete Matta et al. 2009 tail-spot campaign identified in this roadmap. Therefore use the actual observational publications as the source and digitize only when no numerical table is available.

**Matta et al. 2009, Icarus 204, 409–417**
Author-hosted PDF:
`https://sirius.bu.edu/aeronomy/Matta%20et%20al%202009.pdf`

DOI:
`https://doi.org/10.1016/j.icarus.2009.06.017`

**Additional spatial/velocity validation: Smith/Wilson et al. 2012 WHAM study**
Article:
`https://www.sciencedirect.com/science/article/pii/S0019103512001303`

DOI:
`https://doi.org/10.1016/j.icarus.2012.04.001`

**Manual procedure**
1. Save the source PDF/article materials locally and hash them.
2. First extract values from any numerical table (e.g. the 2012 study's tabulated brightest-beam quantities).
3. Where the needed 2-D brightness/profile quantity exists only in a figure, digitize that exact figure with traceability metadata.
4. Record observing UTC/night, aperture/FOV, projected coordinates, brightness, spectral/velocity quantity where applicable, seeing/resolution, page/figure/table, and digitization uncertainty.
5. Perform a second independent digitization/spot check.
6. Freeze the selected nights/profiles before I33 final scoring.
7. Do not renormalize each tail observation to AMPS.


**Authoritative source:** Published ground-based lunar sodium-tail measurements (Matta et al. 2009; complementary WHAM spatial/velocity observations in the 2012 Icarus study).

**Exact acquisition URL(s):** 2009 author PDF: `https://sirius.bu.edu/aeronomy/Matta%20et%20al%202009.pdf`; 2012 article: `https://www.sciencedirect.com/science/article/pii/S0019103512001303`.

**Source:** peer-reviewed lunar sodium-tail observations, including Wilson et al., *Icarus* 204 (2009), DOI `10.1016/j.icarus.2009.06.017`, plus the cited observing papers/data repositories.

**Acquire**
1. First search for author/institution-hosted reduced profiles or machine-readable tables.
2. If a numerical table is unavailable, digitize the published tail-axis/brightness profile from the original paper.
3. Preserve the PDF/figure citation, page/figure number, axis calibration, aperture/slit geometry, observing time, and seeing/resolution information.

**Process**
- Digitize at least twice or independently check a subset of points.
- Carry digitization uncertainty separately from observational uncertainty.
- Reproduce the observation aperture, Earth/Moon geometry, and spatial smoothing in the synthetic observation.
- Do not independently renormalize every held-out profile.

**Validation use:** I33.

**Minimum roadmap processing requirement:** Prefer reduced numerical profiles. If raw/reduced tables cannot be obtained, digitize published profiles with traceable figure/page metadata and measurement/digitization uncertainty.

### D15 — LADEE water release events

#### Exact acquisition — LADEE water-release events

**Open-access article containing the published event analysis**
`https://pmc.ncbi.nlm.nih.gov/articles/PMC7306913/`

Nature DOI:
`https://doi.org/10.1038/s41561-019-0345-3`

On the PMC article page, download the **Supplementary Information** file named:

`NIHMS1535822-supplement-Supplementary_Information.pdf`

Direct PMC file URL (when the browser permits direct binary access):
`https://pmc.ncbi.nlm.nih.gov/articles/PMC7306913/bin/NIHMS1535822-supplement-Supplementary_Information.pdf`

The supplement contains the event/stream information used in the published analysis, including Supplementary Tables referenced by the paper.

**Underlying NMS measurements**
Use the exact LADEE NMS calibrated package from D08:
`https://atmos.nmsu.edu/PDS/data/PDS4/LADEE/nms_data_calibrated.tar.gz`

Also acquire:
- `https://atmos.nmsu.edu/PDS/data/PDS4/LADEE/nms_calibration.tar.gz`
- `https://atmos.nmsu.edu/PDS/data/PDS4/LADEE/nms_document.tar.gz`
- LADEE SPICE kernels from `https://naif.jpl.nasa.gov/pub/naif/LADEE/kernels/`

**Manual procedure**
1. Download and hash the article supplement.
2. Extract Supplementary Tables 1/2 into machine-readable CSV using a committed script or carefully audited PDF-table extraction; preserve the PDF as immutable raw evidence.
3. Independently spot-check every extracted event time/stream assignment against the PDF.
4. Download the calibrated NMS tarball and supporting calibration/docs.
5. Reproduce the paper's NMS water selection from the raw/calibrated data: closed-source `m/z = 18`, first minutes following instrument turn-on, plus the documented background/temperature corrections from the Methods/Supplement.
6. Preserve all 743 candidate water measurements if reproducing the complete analysis; do not keep only detected events.
7. Freeze event windows, detection metric, false-positive windows, and stream associations before running AMPS.
8. For historical meteoroid associations, use the supplement's frozen 2018-era assignments rather than silently substituting the current 2026 MDC list.


**Authoritative source:** Benna et al. (2019) published LADEE water-event analysis/supplement plus calibrated LADEE NMS data.

**Exact acquisition URL(s):** Open-access article/supplement: `https://pmc.ncbi.nlm.nih.gov/articles/PMC7306913/`; calibrated NMS: `https://atmos.nmsu.edu/PDS/data/PDS4/LADEE/nms_data_calibrated.tar.gz`.

**Sources:** LADEE NMS PDS (D08) and Benna et al. (2019), DOI `10.1038/s41561-019-0345-3`.

**Acquire**
1. Download the paper and its supplementary information.
2. Extract the published event list/stream association and reported timing/uncertainty from the article/supplement.
3. Download the underlying LADEE NMS measurements for the same intervals from D08.
4. Match established meteoroid streams using D11; do not add unreported events to improve the model score.

**Process**
- Create `water_events.csv` with event ID, event time/range, stream association, published confidence/status, observed response metric, and citation location.
- Keep the NMS measurement series separately and derive the model-comparison metric with a committed script.
- Predeclare detection threshold, onset-time metric, and false-positive window before running I39.

**Validation use:** I39 and H2O/OH model qualification.

**Minimum roadmap processing requirement:** Use the published event list as event metadata and the archived NMS data as measurement data. Match meteor streams through D11. Do not invent unlisted events.

## 8. Sequential development and validation roadmap

Every roadmap item below should be tracked as an issue/work package. A step is complete only when its production change, data/provenance work, linked tests, and exit criteria are all satisfied.

## Phase 0 — Establish a Trusted Linked Baseline

### [ ] R0.1 — Freeze the current code, tests, and data contract

**Why this step exists**

All later scientific comparisons depend on knowing exactly which source revision, compile-time macros, input, and data versions produced a result.

**Capability/development work**

Create a release-candidate branch/tag; freeze the U01–U26 and I01–I46 IDs; require run_manifest.json for every I-series run; add SHA-256 for input, executable, external data, and processed data.

**Primary code areas**

srcMoon/test/run_tests.py; test/integration_common.py; MoonIntegrationTests.cpp; top-level Moon README.

**Data required**

No external science data required. Existing source/test tree is the reference.

**Where to obtain it**

Repository itself.

**How to prepare/process the data**

Generate a manifest schema and a small validator. Capture compiler, MPI/OpenMP, build ID, git/source revision, species table, random seed, rank/thread layout, physics flags, and data hashes.

**Tests to run at completion**

Run U01–U26 and I01 preflight; repeat identical manifest generation twice; deliberately change one control.

**Acceptance / exit criteria**

All literal U01–U26 IDs remain registered and execute with honest four-state
status.  Local gates directly required by the M0 I01–I10 baseline PASS;
stand-alone tests for production capabilities assigned to later milestones may
remain explicitly SKIPPED and do not become PASS through preflight.  I01–I10
must all PASS for M0.  Identical-control manifests are equivalent; a changed
control is recorded; no unknown/unhashed external file is allowed in a release
run.  This phase-appropriate rule prevents missing M2–M5 physics from being
silently waived while preserving the stated sequential milestone order.

**Expected result**

Every later result is reproducible and attributable to a specific code/data configuration.

**Dependencies / gate**

Do not begin observational scoring until this gate passes.

**Required completion artifacts**

- Updated production code with detailed comments explaining the physical assumption, units, coordinate frame, and failure behavior.
- Updated relevant `README.md` files describing configuration, data dependencies, and exact run command.
- Machine-readable `result.json`/comparison output from the linked test.
- Run manifest and SHA-256 hashes for all external data used.
- Test gates in this step:
  - `U01` — Line-of-sight geometry and output-record verification
  - `U26` — Observation-operator geometry and instrument-coordinate verification
  - `I01` — Reproducible build, configuration, and provenance capture

### [ ] R0.2 — Build and execute the linked I01–I10 baseline

**Why this step exists**

The unit layer can be correct while AMPS wiring, movers, species initialization, callbacks, or generated configuration are wrong. The previous compile error demonstrated that integration APIs must be verified in the actual AMPS tree.

**Capability/development work**

Build the Moon application through the normal AMPS Config.pl/input workflow; remove/recreate generated build directories after source/config changes; resolve any current-version API mismatch without reintroducing obsolete runtime setter APIs.

**Primary code areas**

Generated build/main plus srcMoon/main_lib.cpp, MoonIntegrationTests.cpp, Moon.h, exosphere configuration macros.

**Data required**

No mission data. Analytical references from U03–U06 and committed I01–I10 references.

**Where to obtain it**

Existing test/reference files.

**How to prepare/process the data**

Use the actual generated executable. Run one I-test at a time first, then --phase baseline. Preserve stdout/stderr and manifests. Keep the normal M0 baseline species list as neutral Na. I08 uses a separate committed Na/Na+ configuration with analytical uniform electric and magnetic fields to reach the production ion mover; that fixture does not enable Na+ in the normal science input. Register impact vaporization only through the named built-in process in the M0 input. Retain the disabled historical `MySource` definition for later conserved night-to-day surface-release development, but do not count its present impact-vaporization aliases as reservoir physics.

**Tests to run at completion**

I01–I10 linked runtime, including mover dispatch, particle weights/time steps, gravity, radiation pressure, Lorentz force, termination classification, and dt convergence.

**Acceptance / exit criteria**

No compile errors; no ERROR; all required baseline tests PASS. I10 refined trajectory error <=1% against analytical/reference state and specified convergence ratios.

**Expected result**

A trusted linked application exists before new physics is promoted.

**Dependencies / gate**

Milestone M0.

**Required completion artifacts**

- Updated production code with detailed comments explaining the physical assumption, units, coordinate frame, and failure behavior.
- Updated relevant `README.md` files describing configuration, data dependencies, and exact run command.
- Machine-readable `result.json`/comparison output from the linked test.
- Run manifest and SHA-256 hashes for all external data used.
- Test gates in this step:
  - `I01` — Reproducible build, configuration, and provenance capture
  - `I10` — Production particle-time-step convergence

## Phase 1 — Numerical Convergence and Controlled Baseline

### [ ] R1.1 — Parameterize and converge the AMR mesh

**Why this step exists**

Near-surface Na/Ne/Ar scale lengths can be far below the legacy coarse debug mesh. Scientific validation is meaningless until discretization error is smaller than observational/model uncertainty.

**Capability/development work**

Introduce named debug/production mesh profiles and runtime/configured near-surface target resolution. Preserve the existing AMR logic but remove hidden constants from campaign operation.

**Primary code areas**

main_lib.cpp mesh resolution callbacks; I11 diagnostics and runner.

**Data required**

Campaign-generated convergence outputs only; optional model/I11_mesh_science_convergence.csv.

**Where to obtain it**

Generated by AMPS, not downloaded.

**How to prepare/process the data**

Run the same frozen physics at roughly 40, 20, 10 and 5 km near-surface resolution where practical. Export altitude profiles, surface return flux, column density and synthetic brightness on a common grid.

**Tests to run at completion**

I11 plus I10. Compare successive mesh levels and record wall time/memory.

**Acceptance / exit criteria**

Integrated observables change <5%; pointwise quantities change <10% relative to the next finer mesh. If not met, refine further or restrict claimed spatial resolution.

**Expected result**

A justified production mesh profile and a numerical-error estimate for each key observable.

**Dependencies / gate**

Required before calibration to observations.

**Required completion artifacts**

- Updated production code with detailed comments explaining the physical assumption, units, coordinate frame, and failure behavior.
- Updated relevant `README.md` files describing configuration, data dependencies, and exact run command.
- Machine-readable `result.json`/comparison output from the linked test.
- Run manifest and SHA-256 hashes for all external data used.
- Test gates in this step:
  - `I11` — AMR spatial-resolution convergence
  - `I10` — Production particle-time-step convergence

### [ ] R1.2 — Verify angular surface refinement and terminator behavior

**Why this step exists**

The original source contained a forced SubsolarAngle=0 path; a correct angular refinement profile is essential around terminators, polar terrain and source gradients.

**Capability/development work**

Use the U07-tested selector in production; expose parameters controlling subsolar/terminator/polar refinement; document the profile.

**Primary code areas**

main_lib.cpp localSphericalSurfaceResolution / surface refinement service.

**Data required**

No external data.

**Where to obtain it**

Analytical reference only.

**How to prepare/process the data**

Evaluate requested resolution at SZA 0,30,60,90,120,180 degrees; then build the mesh and measure realized boundary-adjacent cells.

**Tests to run at completion**

I12 linked mesh diagnostic.

**Acceptance / exit criteria**

Requested formula agrees to machine precision; realized cells within +/-25% of requested target except documented AMR constraints.

**Expected result**

Refinement is physically placed rather than accidentally constant.

**Dependencies / gate**

Pass before LOLA production meshing.

**Required completion artifacts**

- Updated production code with detailed comments explaining the physical assumption, units, coordinate frame, and failure behavior.
- Updated relevant `README.md` files describing configuration, data dependencies, and exact run command.
- Machine-readable `result.json`/comparison output from the linked test.
- Run manifest and SHA-256 hashes for all external data used.
- Test gates in this step:
  - `I12` — Surface angular-refinement logic

### [ ] R1.3 — Preserve a controlled spherical-Moon baseline

**Why this step exists**

A simple sphere is needed for analytical verification and to isolate changes caused by topography, thermal physics, or new sources.

**Capability/development work**

Make SURFACE_GEOMETRY=SPHERE an explicit mode. Freeze a corrected-code baseline case and metrics after accepted bug fixes.

**Primary code areas**

main_lib.cpp surface registration; future Moon surface factory; I13.

**Data required**

baseline/I13_archived_baseline_metrics.csv, generated from the approved corrected spherical model—not observations.

**Where to obtain it**

Generate once from the tagged M0/M1 baseline, review, then freeze and checksum.

**How to prepare/process the data**

Store mesh count, species inventory, source totals, selected profiles/brightness, and exact configuration. Do not overwrite automatically.

**Tests to run at completion**

I13 plus U03/U04/U05/U06 regression.

**Acceptance / exit criteria**

Sphere area relative error <=1e-10; approved baseline metrics remain within declared statistical/numerical tolerances; any intentional change is documented.

**Expected result**

Topography can later be assessed by paired SPHERE vs LOLA runs.

**Dependencies / gate**

Milestone M1 when R1.1–R1.3 pass.

**Required completion artifacts**

- Updated production code with detailed comments explaining the physical assumption, units, coordinate frame, and failure behavior.
- Updated relevant `README.md` files describing configuration, data dependencies, and exact run command.
- Machine-readable `result.json`/comparison output from the linked test.
- Run manifest and SHA-256 hashes for all external data used.
- Test gates in this step:
  - `I13` — Spherical-surface baseline preservation
  - `U03` — Lunar gravity kernel and reference trajectories
  - `U04` — Third-body gravity and rotating-frame kernels
  - `U05` — Sodium radiation-pressure geometry and shadow kernels
  - `U06` — Lorentz-acceleration kernel and species-scaling verification

## Phase 2 — Real Lunar Surface, Illumination, Thermal State, and Surface Physics

### [ ] R2.1 — Promote LOLA topography into the production particle-surface boundary

**Why this step exists**

Topography controls impact altitude, local normal, horizon, illumination and cold trapping. A LOLA reader alone does not validate the production surface; the AMPS boundary must actually use it.

**Capability/development work**

Refactor MoonSurfaceData/Topography utilities into a production surface service; add SURFACE_GEOMETRY=LOLA; initially use a global coarse DEM; later allow regional/polar refinement. Build robust ray/triangle or cut-cell intersections and cache a checksum-tagged surface representation.

**Primary code areas**

MoonSurfaceData.h; Topography/lola.cpp; main_lib.cpp surface registration; new/updated surface factory and intersection code.

**Data required**

D01: LDEM_4.IMG + LDEM_4.LBL plus real control points extracted from the same product.

**Where to obtain it**

PDS Geosciences LOLA / LOLA Data Node (D01).

**How to prepare/process the data**

Preserve raw IMG/LBL. Verify 1440x720, 16-bit little-endian integer and label scale. Convert DN to elevation in meters relative to the 1737.4-km reference sphere; do not add the radius to height_m.txt. Normalize longitude to [0,360), sort latitude south-to-north, write lon.txt/lat.txt/height_m.txt, extract actual grid-cell control points, save SHA-256/provenance.

**Tests to run at completion**

U08 ingest/coordinates/normals; I14 production-boundary control points, min/max, crater profiles and vertical-ray intersections.

**Acceptance / exit criteria**

Control-point error <= max(one vertical quantization unit,50 m); 100% ray-hit consistency in designated regions; no inverted/self-intersecting surface elements; production manifest states LOLA product/hash.

**Expected result**

I14 changes from SKIPPED to PASS only when particles actually collide with the LOLA-backed boundary.

**Dependencies / gate**

Foundation for I15–I20 and H2O/OH. Milestone M2 cannot be reached without this.

**Required completion artifacts**

- Updated production code with detailed comments explaining the physical assumption, units, coordinate frame, and failure behavior.
- Updated relevant `README.md` files describing configuration, data dependencies, and exact run command.
- Machine-readable `result.json`/comparison output from the linked test.
- Run manifest and SHA-256 hashes for all external data used.
- Test gates in this step:
  - `U08` — LOLA data ingest, coordinate conversion, and topography geometry verification
  - `I14` — LOLA topography integration into production surface boundary

### [ ] R2.2 — Integrate terrain illumination, local horizons, Earth eclipse and PSR classification

**Why this step exists**

Illumination drives surface temperature, PSD and thermal desorption. Polar volatile science fails if terrain shadowing is simplified to cos(SZA).

**Capability/development work**

Promote Topography/shadow_calc logic into a reusable production illumination service; keep terrain blocking and Earth eclipse as separate terms; optionally precompute horizon azimuth tables for performance.

**Primary code areas**

Topography/shadow_calc.cpp; Moon surface illumination service; source callbacks; temperature boundary.

**Data required**

D01 DEM plus D02 independent LOLA/LRO illumination and PSR products.

**Where to obtain it**

LOLA/ODE GDRPSR or NASA polar illumination products (D02).

**How to prepare/process the data**

Transform reference products to the same lunar coordinate convention without using them to construct the model horizon. Build illumination_reference.csv with lon,lat,epoch/solar direction, observed/reference illuminated flag and product ID. Maintain independent provenance.

**Tests to run at completion**

U09 analytical terrain cases; I15 polar/selected-site classification and sunrise/sunset timing.

**Acceptance / exit criteria**

Sun-visibility agreement >=99% on reference samples; PSR boundaries/areas agree within the resolution of the independent product; flat/sphere limits are exact.

**Expected result**

A trusted illumination state is available to thermal/source modules.

**Dependencies / gate**

Requires R2.1.

**Required completion artifacts**

- Updated production code with detailed comments explaining the physical assumption, units, coordinate frame, and failure behavior.
- Updated relevant `README.md` files describing configuration, data dependencies, and exact run command.
- Machine-readable `result.json`/comparison output from the linked test.
- Run manifest and SHA-256 hashes for all external data used.
- Test gates in this step:
  - `U09` — Topographic illumination and horizon-shadow verification
  - `I15` — Topographic illumination and shadowing

### [ ] R2.3 — Implement a Diviner-backed surface-temperature boundary

**Why this step exists**

Residence times and volatile migration depend exponentially on temperature; the legacy instantaneous dayside + 100 K nightside prescription is inadequate for Ar/water science.

**Capability/development work**

Provide SIMPLE_COSINE (debug), DIVINER_MAP (data-driven), and THERMAL_MODEL (predictive) modes. DIVINER_MAP must interpolate vetted local-time/lat/lon temperature products and expose quality flags.

**Primary code areas**

MoonSurfaceTemperature.h; Moon.cpp GetSurfaceTemperature adapter; data loader/configuration.

**Data required**

D03: Diviner Global Cumulative Products and selected benchmark locations/epochs.

**Where to obtain it**

PDS Geosciences Diviner derived archive (D03).

**How to prepare/process the data**

Download TAB/XML pairs. Select bolometric surface temperature and quality fields; convert local time and coordinate conventions; reject fill/bad-quality samples. Create temperature_table.csv and benchmark_points.csv with product ID, source row/pixel and uncertainty/quality.

**Tests to run at completion**

U10 data ingest/interpolation; I16 production temperature sampling at benchmark points.

**Acceptance / exit criteria**

Non-PSR RMSE <=15 K for the selected validation set; no use of flagged/missing data; coordinate/local-time interpolation reproduced independently.

**Expected result**

Data-driven temperature mode is scientifically traceable and can be used for noble-gas validation.

**Dependencies / gate**

Requires R2.1 and preferably R2.2.

**Required completion artifacts**

- Updated production code with detailed comments explaining the physical assumption, units, coordinate frame, and failure behavior.
- Updated relevant `README.md` files describing configuration, data dependencies, and exact run command.
- Machine-readable `result.json`/comparison output from the linked test.
- Run manifest and SHA-256 hashes for all external data used.
- Test gates in this step:
  - `U10` — Diviner temperature-product ingest and interpolation verification
  - `I16` — Diviner-driven surface-temperature boundary

### [ ] R2.4 — Add time-dependent thermal inertia and eclipse response

**Why this step exists**

A map-only climatology cannot predict arbitrary epochs, transient eclipse cooling, or phase lag through the lunar day.

**Capability/development work**

Implement a 1-D conduction/thermal-inertia boundary model or an equivalent reduced model calibrated only on a training subset of Diviner. Couple illumination from R2.2; support terrain class/albedo/emissivity parameters with provenance.

**Primary code areas**

MoonThermalModel.h; MoonSurfaceTemperature.h; epoch update logic.

**Data required**

D03 time histories / Diviner RDR or time-resolved products; selected eclipse or day/night tracks.

**Where to obtain it**

PDS Diviner RDR/GDR (D03).

**How to prepare/process the data**

Extract temperature-versus-time histories for sites/latitudes used in validation; keep training and hold-out intervals separate. Store eclipse_thermal_history.csv with time, site, observed T, uncertainty, illumination state and quality.

**Tests to run at completion**

U11 solver response; I17 linked surface boundary through sunrise, sunset and eclipse.

**Acceptance / exit criteria**

No unphysical discontinuity; temperature jump at a forcing transition <50 K per integration update and converges with dt; phase/amplitude are within declared Diviner validation tolerance.

**Expected result**

A predictive thermal state suitable for residence-time evolution and water migration.

**Dependencies / gate**

Requires R2.2–R2.3.

**Required completion artifacts**

- Updated production code with detailed comments explaining the physical assumption, units, coordinate frame, and failure behavior.
- Updated relevant `README.md` files describing configuration, data dependencies, and exact run command.
- Machine-readable `result.json`/comparison output from the linked test.
- Run manifest and SHA-256 hashes for all external data used.
- Test gates in this step:
  - `U11` — Thermal-inertia and eclipse-response solver verification
  - `I17` — Time-dependent thermal inertia and eclipse response

### [ ] R2.5 — Promote species-specific adsorption, residence and re-emission physics

**Why this step exists**

A single sticking probability cannot represent Na, noble gases and water. Surface physics strongly controls diurnal structure and residence.

**Capability/development work**

Expose per-species sticking law, adsorption/binding energy or distribution, attempt frequency, accommodation and re-emission kernel. Use the U12-tested bounded interpolation and Arrhenius kernels. Record parameter citations and allowed ranges.

**Primary code areas**

MoonSurfacePhysics.h; Moon.cpp surface interaction; species configuration.

**Data required**

Laboratory/published parameter values; no observation curve is required for numerical verification. Parameter provenance must be recorded.

**Where to obtain it**

Peer-reviewed adsorption/sticking literature referenced in the model configuration; for Ar, include the adsorption framework used in lunar Ar studies; for H2O use experimental/regolith ranges with uncertainty.

**How to prepare/process the data**

Create a parameter table with units, source, temperature validity and uncertainty. Never choose a value solely to improve a held-out comparison.

**Tests to run at completion**

U12; I18 Monte Carlo re-emission distributions and linked surface interaction.

**Acceptance / exit criteria**

Probabilities remain [0,1]; table nodes interpolate exactly; stochastic moments/CDF lie within 95% confidence of analytic distributions; no out-of-range access.

**Expected result**

Every species has a physically explicit surface model instead of implicit defaults.

**Dependencies / gate**

Requires R2.3/R2.4 for temperature-dependent residence.

**Required completion artifacts**

- Updated production code with detailed comments explaining the physical assumption, units, coordinate frame, and failure behavior.
- Updated relevant `README.md` files describing configuration, data dependencies, and exact run command.
- Machine-readable `result.json`/comparison output from the linked test.
- Run manifest and SHA-256 hashes for all external data used.
- Test gates in this step:
  - `U12` — Species sticking, residence-time, and re-emission distribution verification
  - `I18` — Species-specific sticking, residence time, and re-emission distributions

### [ ] R2.6 — Close surface inventories and global mass budgets

**Why this step exists**

Source, adsorption, chemistry, escape and numerical deletion must balance. Without a budget, an apparently good density map can hide a conservation defect.

**Capability/development work**

Add particle/s and kg/s counters for injection, desorption, adsorption, interspecies conversion, escape, domain loss and inventory change; distinguish physical sinks from numerical losses.

**Primary code areas**

MoonIntegrationTests.cpp diagnostics; Moon surface reservoir/source modules; common budget recorder.

**Data required**

No external data.

**Where to obtain it**

Analytical/deterministic bookkeeping references plus campaign output.

**How to prepare/process the data**

For each process, record signed flux and uncertainty over synchronized windows. For conversions, conserve nuclei/mass across parent/product species. Write a machine-readable budget table.

**Tests to run at completion**

I19 across Na and at least one noble gas; later rerun for H2O/OH.

**Acceptance / exit criteria**

Deterministic closure to machine precision; full Monte Carlo integrated imbalance <0.5% of the dominant source after sufficient averaging, or a stricter/justified process-specific limit.

**Expected result**

Budget closure becomes a blocking release criterion.

**Dependencies / gate**

Requires R2.5; rerun after every new source/chemistry module.

**Required completion artifacts**

- Updated production code with detailed comments explaining the physical assumption, units, coordinate frame, and failure behavior.
- Updated relevant `README.md` files describing configuration, data dependencies, and exact run command.
- Machine-readable `result.json`/comparison output from the linked test.
- Run manifest and SHA-256 hashes for all external data used.
- Test gates in this step:
  - `I19` — Surface inventory and global mass-conservation accounting

### [ ] R2.7 — Implement PSR/cold-trap retention as a physical state

**Why this step exists**

A geometric PSR flag is not equivalent to permanent trapping. Species retention must depend on temperature and residence time relative to relevant timescales.

**Capability/development work**

Combine terrain/illumination, thermal state and species binding kinetics; distinguish temporary adsorption from cold-trapped inventory; allow re-release if a site warms.

**Primary code areas**

MoonColdTrap.h; MoonSurfacePhysics.h; surface inventory.

**Data required**

D02 PSR reference + D03 polar thermal products; optional published constraints on Ar/volatile trapping.

**Where to obtain it**

LOLA/ODE polar illumination and PDS Diviner polar/resource products (D02/D03).

**How to prepare/process the data**

Collocate PSR and Diviner polar temperatures on the production surface; create independent reference sites with illumination fraction and temperature statistics. Do not classify traps solely from a binary external PSR map.

**Tests to run at completion**

U13; I20 site residence statistics, trapped fraction, release on warming, global inventory sensitivity.

**Acceptance / exit criteria**

Trapping probabilities are bounded/converged; site residence follows the kinetic law; geometric PSR and thermal-retention classifications are both reported; no permanent trap is assigned above the species-specific thermal threshold.

**Expected result**

Topography+temperature+surface physics is ready for Ar and H2O/OH science.

**Dependencies / gate**

Completes Milestone M2 with R2.1–R2.6.

**Required completion artifacts**

- Updated production code with detailed comments explaining the physical assumption, units, coordinate frame, and failure behavior.
- Updated relevant `README.md` files describing configuration, data dependencies, and exact run command.
- Machine-readable `result.json`/comparison output from the linked test.
- Run manifest and SHA-256 hashes for all external data used.
- Test gates in this step:
  - `U13` — PSR/cold-trap classification and trapping-kernel verification
  - `I20` — Permanent-shadow and polar cold-trapping physics

## Phase 3 — Chemistry Promotion and Budget Closure

### [ ] R3.1 — Qualify Na photoionization and Na-to-Na+ conversion

**Why this step exists**

Na neutral loss and Na+ production affect exosphere/tail lifetime and ion access; deterministic rate and stochastic conversion must agree.

**Capability/development work**

Use MoonPhotochemistry shared kernel through the configured AMPS photolytic macros. Make the selected 1-AU rate/spectrum a documented configuration item and keep lunar/Earth shadow gating explicit.

**Primary code areas**

MoonPhotochemistry.h; Moon.h lifetime/conversion adapter; AMPS exosphere photolytic configuration.

**Data required**

D04 PHIDRATES or another explicitly selected peer-reviewed Na rate source; model-generated I21 survival table.

**Where to obtain it**

PHIDRATES (D04) and cited Na photoionization literature.

**How to prepare/process the data**

Record the selected solar spectrum/activity case and 1-AU rate. Generate deterministic survival curves from the production rate; separately run particle ensembles at multiple distances/shadow states.

**Tests to run at completion**

U14; I21 survival and conversion/budget test.

**Acceptance / exit criteria**

Ensemble lifetime relative error <=2%; r^-2 rate scaling error <=0.1%; Na->Na+ conversion closes the species budget.

**Expected result**

Na chemistry can be used in observation validation without hidden rate choices.

**Dependencies / gate**

R2.6 budget gate required.

**Required completion artifacts**

- Updated production code with detailed comments explaining the physical assumption, units, coordinate frame, and failure behavior.
- Updated relevant `README.md` files describing configuration, data dependencies, and exact run command.
- Machine-readable `result.json`/comparison output from the linked test.
- Run manifest and SHA-256 hashes for all external data used.
- Test gates in this step:
  - `U14` — Na photoionization-rate and survival-kernel verification
  - `I21` — Sodium photoionization survival and Na→Na+ conversion

### [ ] R3.2 — Add documented photo-process tables for He, Ne, Ar, H2O and OH

**Why this step exists**

Multi-species validation requires species-specific photon rates and branching rather than reusing Na behavior.

**Capability/development work**

Create a versioned photochemical table with species, parent, products, total rate, branching, solar-spectrum identifier and validity. Only enable channels with documented data.

**Primary code areas**

MoonPhotochemistry.h; product-species mapping; configuration/provenance.

**Data required**

D04 PHIDRATES species data.

**Where to obtain it**

PHIDRATES (D04).

**How to prepare/process the data**

Export or transcribe authoritative rates/branching with citations; store raw/reference values and a processing script; independently recompute total/branch sums. For H2O/OH preserve separate ionization and dissociation products.

**Tests to run at completion**

U15/U24; I22 linked ensemble rate/branching tests.

**Acceptance / exit criteria**

Ensemble rate error <=2%; branching fractions statistically consistent at 95% confidence; all branches sum to 1 within numerical precision after excluding explicitly unmodeled channels.

**Expected result**

Photon chemistry is complete enough for species/source validation.

**Dependencies / gate**

Requires species definitions/product masses/charges to be stable.

**Required completion artifacts**

- Updated production code with detailed comments explaining the physical assumption, units, coordinate frame, and failure behavior.
- Updated relevant `README.md` files describing configuration, data dependencies, and exact run command.
- Machine-readable `result.json`/comparison output from the linked test.
- Run manifest and SHA-256 hashes for all external data used.
- Test gates in this step:
  - `U15` — Multi-species photoionization table/rate verification
  - `U24` — H2O/OH surface-migration and photochemistry kernel verification
  - `I22` — Photoionization for He, Ne, Ar and future H2O/OH

### [ ] R3.3 — Replace the legacy electron-impact approximation with species-specific rate coefficients

**Why this step exists**

The old sigma times 400 km/s approximation is not an electron-impact rate model. Electron ionization depends on the electron energy distribution.

**Capability/development work**

Use species-specific Maxwellian rate coefficients or integrate cross sections over an imported electron distribution; share the U16 kernel with production; feed measured/model electron density and Te.

**Primary code areas**

MoonElectronImpact.h; Moon.cpp electron-impact adapter; plasma driver interface.

**Data required**

D05 published fit coefficients/tables; D06 electron moments for event forcing.

**Where to obtain it**

Voronov/Verner data and ARTEMIS when running real events.

**How to prepare/process the data**

Freeze coefficient table; make a temperature grid reference; when using ARTEMIS, filter on quality and interpolate n_e/Te to simulation time without filling large gaps.

**Tests to run at completion**

U16; I23 direct coefficient and ensemble lifetime tests.

**Acceptance / exit criteria**

Rate coefficient relative error <=2%; ensemble lifetime error <=3%; zero/invalid density/temperature handled explicitly.

**Expected result**

Electron-impact chemistry is physically interpretable and event-driven.

**Dependencies / gate**

Requires R4.1 for measured event forcing; kernel can be verified earlier.

**Required completion artifacts**

- Updated production code with detailed comments explaining the physical assumption, units, coordinate frame, and failure behavior.
- Updated relevant `README.md` files describing configuration, data dependencies, and exact run command.
- Machine-readable `result.json`/comparison output from the linked test.
- Run manifest and SHA-256 hashes for all external data used.
- Test gates in this step:
  - `U16` — Electron-impact ionization coefficient verification
  - `I23` — Electron-impact ionization rate coefficients

### [ ] R3.4 — Add only documented charge-exchange channels

**Why this step exists**

Charge exchange can control ion/neutral coupling, but invented cross sections are worse than an omitted channel.

**Capability/development work**

Build a reaction registry keyed by reactants/products and source reference; implement interpolation/fit only over valid ranges; expose unsupported channels as disabled.

**Primary code areas**

MoonChargeExchange.h; chemistry execution/budget integration.

**Data required**

Published charge-transfer rate/cross-section datasets for selected channels.

**Where to obtain it**

Atomic-data papers/databases cited in the reaction registry. The existing U17 documented channel can be the first enabled reaction.

**How to prepare/process the data**

Store raw coefficients/range/units and a small independent reference table. Do not extrapolate silently beyond validity.

**Tests to run at completion**

U17; I24 direct-rate and ensemble lifetime/product tests.

**Acceptance / exit criteria**

Direct rate error <=1%; ensemble lifetime error <=3%; product identities and charge/mass conservation exact; unsupported channels remain explicitly OFF.

**Expected result**

Chemistry scope is scientifically defensible.

**Dependencies / gate**

Milestone M3 when R3.1–R3.4 and budgets pass.

**Required completion artifacts**

- Updated production code with detailed comments explaining the physical assumption, units, coordinate frame, and failure behavior.
- Updated relevant `README.md` files describing configuration, data dependencies, and exact run command.
- Machine-readable `result.json`/comparison output from the linked test.
- Run manifest and SHA-256 hashes for all external data used.
- Test gates in this step:
  - `U17` — Charge-exchange reaction-kernel verification
  - `I24` — Charge-exchange reactions with solar-wind ions

## Phase 4 — Measured Drivers and Species Source Models

### [ ] R4.1 — Implement measured time-dependent plasma forcing

**Why this step exists**

He/Ne implantation, sputtering, electron impact and pickup-ion transport respond to the real solar wind/magnetospheric environment. Constant “typical” conditions cannot reproduce event variability.

**Capability/development work**

Connect MoonPlasmaDrivers to production; ingest time-tagged density, velocity, electron temperature and B (and E if available/derived with documented convention); include gap/quality handling and coordinate transforms.

**Primary code areas**

MoonPlasmaDrivers.h; coupler/background-field adapters; driver configuration.

**Data required**

D06 ARTEMIS P1/P2 data; optionally OMNI for upstream context, but lunar-local forcing should preferentially use ARTEMIS.

**Where to obtain it**

CDAWeb THB/THC L2 MERGED/MOM/FGM and ARTEMIS metadata (D06).

**How to prepare/process the data**

Download CDF/CSV for event windows with position. Keep quality flags. Convert coordinates to the AMPS lunar frame with tested rotations; interpolate only across permitted gaps; write plasma_driver.csv with source variables, uncertainty/quality and spacecraft distance to Moon.

**Tests to run at completion**

U18; I25 driver ingest, interpolation, gap behavior and linked field/source response.

**Acceptance / exit criteria**

At original time stamps, production values reproduce downloaded variables within floating/interpolation precision; bad/gap intervals are flagged rather than filled silently; coordinate transforms pass U18 references.

**Expected result**

A single authoritative event driver feeds sources, chemistry and ion dynamics.

**Dependencies / gate**

Required before observational He/ion validation.

**Required completion artifacts**

- Updated production code with detailed comments explaining the physical assumption, units, coordinate frame, and failure behavior.
- Updated relevant `README.md` files describing configuration, data dependencies, and exact run command.
- Machine-readable `result.json`/comparison output from the linked test.
- Run manifest and SHA-256 hashes for all external data used.
- Test gates in this step:
  - `U18` — Solar-wind/magnetospheric driver ingest and interpolation verification
  - `I25` — Measured/time-dependent solar-wind and magnetospheric forcing

### [ ] R4.2 — Promote the He implantation/reservoir source

**Why this step exists**

Helium provides a clean solar-wind-driven validation target and tests the surface reservoir response timescale.

**Capability/development work**

Use measured alpha-particle forcing when available; model implantation/neutralization efficiency, surface reservoir, thermal release/escape and optional endogenous component. Do not treat total ion density as alpha flux.

**Primary code areas**

MoonHeliumSource.h; Exosphere_Helium.cpp production adapter; surface inventory.

**Data required**

D06 alpha-specific forcing or published/archived alpha series; D08/D09 observations for later validation.

**Where to obtain it**

ARTEMIS distributions/published extraction; LADEE/LAMP observations.

**How to prepare/process the data**

If a directly archived alpha flux is unavailable, reproduce a published extraction from ARTEMIS distributions and document energy/species selection. Keep source forcing separate from validation observations.

**Tests to run at completion**

U19; I26 source/reservoir conservation and response tests; later I34/I35.

**Acceptance / exit criteria**

Reservoir ODE/analytic limits match U19; source/inventory closes I19; response is causal and positive; no calibration to I34/I35 hold-out cases.

**Expected result**

He source ready for two independent validation datasets.

**Dependencies / gate**

Requires R2.5/R2.6 and R4.1.

**Required completion artifacts**

- Updated production code with detailed comments explaining the physical assumption, units, coordinate frame, and failure behavior.
- Updated relevant `README.md` files describing configuration, data dependencies, and exact run command.
- Machine-readable `result.json`/comparison output from the linked test.
- Run manifest and SHA-256 hashes for all external data used.
- Test gates in this step:
  - `U19` — Helium source/reservoir response-kernel verification
  - `I26` — Helium source: solar-wind alpha implantation/re-emission plus endogenous component
  - `I34` — Helium local-time/altitude distribution — LADEE NMS
  - `I35` — Helium solar-wind response — LRO/LAMP coordinated with ARTEMIS

### [ ] R4.3 — Promote the Ne source and surface accommodation model

**Why this step exists**

Ne is largely solar-wind related but has different mass/accommodation behavior from He, providing an independent test of the source/surface framework.

**Capability/development work**

Drive Ne implantation with measured or documented solar-wind heavy-ion forcing; preserve immediate-release versus accommodated-reservoir fractions as explicit parameters.

**Primary code areas**

MoonNeonSource.h; Exosphere_Neon.cpp; surface reservoir.

**Data required**

D06 heavy-ion context; D08 LADEE Ne for later validation.

**Where to obtain it**

ARTEMIS or a documented solar-wind composition product plus LADEE NMS.

**How to prepare/process the data**

Use a documented Ne/solar-wind abundance source rather than an assumed constant if event-specific composition is required; record uncertainty and separate source input from validation data.

**Tests to run at completion**

U20; I27; later I36.

**Acceptance / exit criteria**

Correct limiting cases (zero illumination/source, zero accommodation, full accommodation); inventory closes; source scaling and temporal response are reproducible.

**Expected result**

Ne source ready for LADEE validation.

**Dependencies / gate**

Requires R4.1 and surface physics.

**Required completion artifacts**

- Updated production code with detailed comments explaining the physical assumption, units, coordinate frame, and failure behavior.
- Updated relevant `README.md` files describing configuration, data dependencies, and exact run command.
- Machine-readable `result.json`/comparison output from the linked test.
- Run manifest and SHA-256 hashes for all external data used.
- Test gates in this step:
  - `U20` — Neon source and accommodation-kernel verification
  - `I27` — Neon source and noncondensable surface accommodation
  - `I36` — Neon local-time distribution — LADEE NMS

### [ ] R4.4 — Implement radiogenic 40Ar source geography and transients

**Why this step exists**

Argon is fundamentally radiogenic rather than a solar-wind species; source geography and thermal adsorption drive observed variability.

**Capability/development work**

Add a baseline 40Ar production map tied to a documented K-abundance proxy and a separately parameterized transient/outgassing component. Absolute normalization must be calibrated on designated Ar training data only.

**Primary code areas**

MoonArgonSource.h; surface source distribution; configuration.

**Data required**

D10 Kaguya GRS K map; D08/D12 Ar observations for calibration/validation.

**Where to obtain it**

DARTS Kaguya GRS (D10).

**How to prepare/process the data**

Extract the K layer and map geometry; mask invalid cells; resample conservatively to the production surface; normalize to unit integrated weight. Keep raw count/intensity units and limitations. Fit only one absolute Ar production factor on a designated training subset, then freeze it.

**Tests to run at completion**

U21; I28 source-map normalization/transient gating; later I37/I38.

**Acceptance / exit criteria**

Spatial weights integrate to 1; no negative source; map transform is reproducible; transient source is additive and independently identifiable; held-out Ar validation is not used for retuning.

**Expected result**

Physically motivated Ar source ready for thermal/cold-trap validation.

**Dependencies / gate**

Requires R2.7 and frozen training/hold-out policy.

**Required completion artifacts**

- Updated production code with detailed comments explaining the physical assumption, units, coordinate frame, and failure behavior.
- Updated relevant `README.md` files describing configuration, data dependencies, and exact run command.
- Machine-readable `result.json`/comparison output from the linked test.
- Run manifest and SHA-256 hashes for all external data used.
- Test gates in this step:
  - `U21` — Radiogenic 40Ar source-map and transient-source verification
  - `I28` — Radiogenic 40Ar source map and transient outgassing option
  - `I37` — Argon diurnal and selenographic structure — LADEE NMS
  - `I38` — Argon sunrise pocket, lunation behavior, and polar sequestration — Apollo 17 LACE + LADEE synthesis

### [ ] R4.5 — Make Na source partition explicit and independently testable

**Why this step exists**

Na brightness/tail can be matched for the wrong reason if PSD, thermal desorption, sputtering and impact vaporization are freely traded against one another.

**Capability/development work**

Use separate source modules with independent parameters, source IDs and budget counters. Allow mechanisms to be enabled one at a time and in combinations.

**Primary code areas**

MoonNaSources.h; exosphere source callbacks; source budget output.

**Data required**

Laboratory/literature source parameters; D06 for sputtering forcing; D11 for impacts; D07/D14 for later validation.

**Where to obtain it**

Cited Na desorption/sputtering/impact literature plus measured drivers.

**How to prepare/process the data**

Record each parameter with unit/source/range. Build a source decomposition table over local time/latitude and verify normalization of every velocity distribution.

**Tests to run at completion**

U22; I29; budget I19; later I32/I33.

**Acceptance / exit criteria**

Each mechanism reproduces its unit/analytic normalization; summed production equals the sum of source counters; disabling a source makes its contribution exactly zero.

**Expected result**

Na comparisons can diagnose source physics rather than only total amplitude.

**Dependencies / gate**

Required before final Na observational validation.

**Required completion artifacts**

- Updated production code with detailed comments explaining the physical assumption, units, coordinate frame, and failure behavior.
- Updated relevant `README.md` files describing configuration, data dependencies, and exact run command.
- Machine-readable `result.json`/comparison output from the linked test.
- Run manifest and SHA-256 hashes for all external data used.
- Test gates in this step:
  - `U22` — Na source-process kernel verification
  - `I29` — Sodium source partition: PSD, thermal desorption, sputtering, and impact vaporization
  - `I19` — Surface inventory and global mass-conservation accounting
  - `I32` — Sodium limb/exosphere brightness — SELENE/Kaguya UPI-TVIS
  - `I33` — Sodium extended tail and Earth-focusing validation

### [ ] R4.6 — Add time-dependent meteoroid impact forcing

**Why this step exists**

Impact vaporization and water-release events are episodic; a constant source cannot validate meteor-associated observations.

**Capability/development work**

Ingest background and shower-dependent meteoroid mass flux/velocity/impact geometry; connect to species yields and impact vaporization; preserve event identity.

**Primary code areas**

MoonMeteoroidDriver.h; impact-vaporization source.

**Data required**

D11 IAU MDC plus published lunar meteoroid flux model; D15 water-event list for later validation.

**Where to obtain it**

IAU MDC (D11) and cited flux/yield literature.

**How to prepare/process the data**

Download shower parameters/status; compute Moon encounter epochs/radiants using SPICE and published orbital parameters; propagate uncertainty. Do not derive forcing amplitude from the same NMS water signal used for validation.

**Tests to run at completion**

U23; I30 driver/source response; later I39.

**Acceptance / exit criteria**

Driver interpolation and event timing exact at tabulated times; impact source scales monotonically with mass flux and velocity under the chosen yield law; event/background contributions are separately reported.

**Expected result**

A traceable transient source exists for Na and water studies.

**Dependencies / gate**

Requires source budget and SPICE geometry.

**Required completion artifacts**

- Updated production code with detailed comments explaining the physical assumption, units, coordinate frame, and failure behavior.
- Updated relevant `README.md` files describing configuration, data dependencies, and exact run command.
- Machine-readable `result.json`/comparison output from the linked test.
- Run manifest and SHA-256 hashes for all external data used.
- Test gates in this step:
  - `U23` — Meteoroid forcing and impact-vaporization driver verification
  - `I30` — Time-dependent meteoroid forcing and impact vaporization
  - `I39` — Meteoroid-driven exospheric water releases — LADEE NMS

### [ ] R4.7 — Promote H2O/OH surface migration and photochemistry

**Why this step exists**

Water/OH transport, hopping, photodissociation and PSR trapping are absent from the legacy application but are required for modern volatile science.

**Capability/development work**

Add H2O and OH species; adsorption/desorption/hopping, photodissociation branches, ballistic transport, cold-trap inventory and optional impact source. Use the already verified local kernels as the production implementation.

**Primary code areas**

MoonWater.h; MoonPhotochemistry.h; MoonColdTrap.h; surface/source modules.

**Data required**

D03 polar temperatures, D04 H2O/OH photochemistry, D11 meteoroids, D15 LADEE water events.

**Where to obtain it**

Diviner + PHIDRATES + IAU MDC + LADEE NMS.

**How to prepare/process the data**

Freeze photochemical branches/rates and binding-energy priors before looking at hold-out water events. Use real polar thermal/topographic state. Track H/O nuclei across reactions.

**Tests to run at completion**

U24; I31 local migration/budget; later I39 event validation.

**Acceptance / exit criteria**

Local rate/branch tests pass; global H/O nuclei budget closes; cold-trap inventory converges with mesh/time step/particle number; no negative/undefined residence times.

**Expected result**

Water/OH becomes EXPERIMENTAL/VERIFIED before observational validation.

**Dependencies / gate**

Milestone M4 when R4.1–R4.7 relevant enabled capabilities pass.

**Required completion artifacts**

- Updated production code with detailed comments explaining the physical assumption, units, coordinate frame, and failure behavior.
- Updated relevant `README.md` files describing configuration, data dependencies, and exact run command.
- Machine-readable `result.json`/comparison output from the linked test.
- Run manifest and SHA-256 hashes for all external data used.
- Test gates in this step:
  - `U24` — H2O/OH surface-migration and photochemistry kernel verification
  - `I31` — H2O/OH surface migration and photochemistry
  - `I39` — Meteoroid-driven exospheric water releases — LADEE NMS

## Phase 5 — Observation-Based Species and Phenomenon Validation

### [ ] R5.1 — Validate Na limb/exosphere brightness against Kaguya UPI-TVIS

**Why this step exists**

This tests the complete neutral Na chain: sources, ballistic transport, radiation pressure, photoionization, surface interaction, observation geometry and brightness conversion.

**Capability/development work**

Implement a Kaguya TVIS observation operator using actual spacecraft/line-of-sight geometry and instrument calibration; establish explicit calibration and hold-out image sets.

**Primary code areas**

tvis.Kaguya.cpp; MoonObservationGeometry.h; brightness sampler; I32 processor.

**Data required**

D07 TVIS Na Level 2A IMG/LBL, dark/background products, SPICE/trajectory, published calibration context.

**Where to obtain it**

DARTS dataset darts:sln-e-tvis-5-na-level2a-v1.0.

**How to prepare/process the data**

Download original image+label pairs. Apply only documented calibration/background treatment; convert to Rayleigh if not already calibrated; generate per-pixel/annular comparison tables with observed brightness, uncertainty, geometry and held_out flag. Store product IDs and hashes.

**Tests to run at completion**

I32 on multiple images/local times; use U26 geometry checks.

**Acceptance / exit criteria**

Reduced chi-square <=2 for calibration cases and <=3 for held-out cases; integrated brightness within 25%; profile correlation r>=0.9.

**Expected result**

Na near-Moon exosphere is observationally VALIDATED for the tested conditions.

**Dependencies / gate**

Requires M1–M4 Na-related gates.

**Required completion artifacts**

- Updated production code with detailed comments explaining the physical assumption, units, coordinate frame, and failure behavior.
- Updated relevant `README.md` files describing configuration, data dependencies, and exact run command.
- Machine-readable `result.json`/comparison output from the linked test.
- Run manifest and SHA-256 hashes for all external data used.
- Test gates in this step:
  - `I32` — Sodium limb/exosphere brightness — SELENE/Kaguya UPI-TVIS
  - `U26` — Observation-operator geometry and instrument-coordinate verification

### [ ] R5.2 — Validate the extended Na tail and Earth focusing

**Why this step exists**

Near-Moon agreement does not validate long-distance radiation-pressure and Earth-gravity transport.

**Capability/development work**

Generate Earth-view/tail-axis synthetic observations with the same apertures/resolution as selected ground-based observations; do not renormalize each held-out case independently.

**Primary code areas**

Moon_SampleVelocityDistribution.cpp; tail samplers; I33.

**Data required**

D14 reduced ground-based sodium-tail profiles or carefully digitized published figures.

**Where to obtain it**

Peer-reviewed sodium-tail datasets/papers; DOI 10.1016/j.icarus.2009.06.017 is a key reference.

**How to prepare/process the data**

Prefer author tables. If unavailable, digitize published axes/profiles with recorded figure/page, pixel-to-axis calibration, two-pass digitization and uncertainty. Preserve seeing/angular resolution and observation time.

**Tests to run at completion**

I33 axis, width, integrated brightness and velocity centroid where measured.

**Acceptance / exit criteria**

Tail axis within one resolution element; width within 20%; integrated brightness within 30% after the previously frozen source calibration; velocity centroid within measurement uncertainty.

**Expected result**

Na tail dynamics is validated independently of the near-Moon TVIS fit.

**Dependencies / gate**

Requires R5.1 parameter freeze.

**Required completion artifacts**

- Updated production code with detailed comments explaining the physical assumption, units, coordinate frame, and failure behavior.
- Updated relevant `README.md` files describing configuration, data dependencies, and exact run command.
- Machine-readable `result.json`/comparison output from the linked test.
- Run manifest and SHA-256 hashes for all external data used.
- Test gates in this step:
  - `I33` — Sodium extended tail and Earth-focusing validation

### [ ] R5.3 — Validate He altitude/local-time structure against LADEE NMS

**Why this step exists**

LADEE provides direct in-situ abundance profiles and tests source/reservoir/escape physics independently of remote sensing.

**Capability/development work**

Collocate AMPS He density to actual LADEE trajectory and observation times. Designate a calibration interval for at most a small set of He source parameters; hold out the rest.

**Primary code areas**

Observation collocator; I34.

**Data required**

D08 LADEE NMS derived He plus spacecraft geometry/quality.

**Where to obtain it**

NASA PDS LADEE NMS derived collection DOI 10.17189/1408897.

**How to prepare/process the data**

Download derived He and supporting calibrated/geometry products. Filter quality; compute altitude/local time/latitude; bin only with a documented rule. Keep individual points so binning can be audited.

**Tests to run at completion**

I34 time series, local-time and altitude profiles.

**Acceptance / exit criteria**

Correlation >=0.7; model/observation median ratio 0.67–1.5; peak local-solar-time error <=1 h for the selected comparison.

**Expected result**

He is VALIDATED against in-situ observations.

**Dependencies / gate**

Requires R4.2 and thermal/surface gates.

**Required completion artifacts**

- Updated production code with detailed comments explaining the physical assumption, units, coordinate frame, and failure behavior.
- Updated relevant `README.md` files describing configuration, data dependencies, and exact run command.
- Machine-readable `result.json`/comparison output from the linked test.
- Run manifest and SHA-256 hashes for all external data used.
- Test gates in this step:
  - `I34` — Helium local-time/altitude distribution — LADEE NMS

### [ ] R5.4 — Validate He response to solar-wind forcing using LAMP + ARTEMIS

**Why this step exists**

This tests whether temporal He variability is driven correctly rather than merely matching a static LADEE profile.

**Capability/development work**

Collocate LAMP He retrieval epochs with ARTEMIS solar-wind/alpha forcing and model reservoir response. Keep forcing data and validation data independent.

**Primary code areas**

He campaign post-processing; I35.

**Data required**

D09 LAMP helium observations + D06 ARTEMIS plasma/alpha forcing.

**Where to obtain it**

PDS Imaging LAMP archive; CDAWeb/ARTEMIS.

**How to prepare/process the data**

Download exact LAMP products used for the published He retrieval; reproduce retrieval/preprocessing as documented. Match to ARTEMIS epochs/positions. If an alpha-specific time series is not directly archived, reproduce a cited extraction—do not substitute a constant composition ratio.

**Tests to run at completion**

I35 temporal correlation and amplitude response.

**Acceptance / exit criteria**

Correlation >=0.7 and typical absolute relative error <=50% under the declared uncertainty model; response lag/timescale is physically consistent and not tuned per event.

**Expected result**

Independent remote-sensing validation of the He source/reservoir response.

**Dependencies / gate**

Requires R5.3 frozen He parameters.

**Required completion artifacts**

- Updated production code with detailed comments explaining the physical assumption, units, coordinate frame, and failure behavior.
- Updated relevant `README.md` files describing configuration, data dependencies, and exact run command.
- Machine-readable `result.json`/comparison output from the linked test.
- Run manifest and SHA-256 hashes for all external data used.
- Test gates in this step:
  - `I35` — Helium solar-wind response — LRO/LAMP coordinated with ARTEMIS

### [ ] R5.5 — Validate Ne against LADEE NMS

**Why this step exists**

Ne provides an independent noncondensable species test of solar-wind source and surface accommodation.

**Capability/development work**

Collocate the frozen Ne model with LADEE Ne observations using the same geometry/quality pipeline as He.

**Primary code areas**

I36.

**Data required**

D08 LADEE NMS derived Ne.

**Where to obtain it**

NASA PDS LADEE NMS.

**How to prepare/process the data**

Use the same reproducible PDS loader; do not refit common surface/thermal parameters already frozen by He unless the parameter is demonstrably species-specific and designated for training.

**Tests to run at completion**

I36 local-time/altitude distribution.

**Acceptance / exit criteria**

Correlation >=0.8; median model/observation ratio 0.5–2.0.

**Expected result**

Ne capability becomes observationally validated for the tested interval.

**Dependencies / gate**

Requires R4.3 and M2 thermal/surface physics.

**Required completion artifacts**

- Updated production code with detailed comments explaining the physical assumption, units, coordinate frame, and failure behavior.
- Updated relevant `README.md` files describing configuration, data dependencies, and exact run command.
- Machine-readable `result.json`/comparison output from the linked test.
- Run manifest and SHA-256 hashes for all external data used.
- Test gates in this step:
  - `I36` — Neon local-time distribution — LADEE NMS

### [ ] R5.6 — Validate Ar diurnal and selenographic structure against LADEE NMS

**Why this step exists**

Ar is the strongest integrated test of radiogenic geography, adsorption/desorption, thermal inertia and cold trapping.

**Capability/development work**

Run the frozen K-proxy source + surface/thermal model along the actual LADEE trajectory; allow only predeclared Ar-specific normalization/binding parameters to be calibrated on training data.

**Primary code areas**

I37.

**Data required**

D08 LADEE NMS Ar + D10 K map.

**Where to obtain it**

NASA PDS LADEE NMS; DARTS Kaguya GRS.

**How to prepare/process the data**

Prepare time/altitude/local-time/longitude observations with quality filters. Partition by time or longitude into calibration and hold-out sets before fitting.

**Tests to run at completion**

I37 abundance ratio, spatial/local-time shape and source-location diagnostics.

**Acceptance / exit criteria**

At least 80% of scored samples within factor 2; longitude/peak location error <=15 degrees for designated structures.

**Expected result**

Ar diurnal/geographic behavior is validated without using a free map fitted directly to LADEE.

**Dependencies / gate**

Requires R2.4–R2.7 and R4.4.

**Required completion artifacts**

- Updated production code with detailed comments explaining the physical assumption, units, coordinate frame, and failure behavior.
- Updated relevant `README.md` files describing configuration, data dependencies, and exact run command.
- Machine-readable `result.json`/comparison output from the linked test.
- Run manifest and SHA-256 hashes for all external data used.
- Test gates in this step:
  - `I37` — Argon diurnal and selenographic structure — LADEE NMS

### [ ] R5.7 — Cross-check Ar sunrise and long-timescale behavior with Apollo 17 LACE

**Why this step exists**

LADEE alone may leave degeneracies in adsorption/cold-trap parameters; the independent Apollo epoch adds a powerful temporal cross-check.

**Capability/development work**

Reproduce LACE observation geometry/site local time and compare the already frozen Ar model; avoid retuning to Apollo unless a new epoch-dependent source hypothesis is explicitly tested.

**Primary code areas**

I38.

**Data required**

D12 LACE report/data plus D08 LADEE.

**Where to obtain it**

NASA NTRS report 19760025001 and peer-reviewed Ar analyses.

**How to prepare/process the data**

Prefer numeric tables. If digitizing a plot, preserve raw image/PDF, page/figure, axis calibration, digitized points and estimated digitization error. Independently repeat a subset to quantify extraction reproducibility.

**Tests to run at completion**

I38 sunrise peak timing/amplitude and long-timescale background.

**Acceptance / exit criteria**

Sunrise peak timing error <=2 lunar hours; peak/background relative error <=30% for the declared comparison; polar-sequestration conclusions consistent across LADEE and LACE within uncertainties.

**Expected result**

Ar model gains cross-mission validation.

**Dependencies / gate**

Requires R5.6 parameters frozen.

**Required completion artifacts**

- Updated production code with detailed comments explaining the physical assumption, units, coordinate frame, and failure behavior.
- Updated relevant `README.md` files describing configuration, data dependencies, and exact run command.
- Machine-readable `result.json`/comparison output from the linked test.
- Run manifest and SHA-256 hashes for all external data used.
- Test gates in this step:
  - `I38` — Argon sunrise pocket, lunation behavior, and polar sequestration — Apollo 17 LACE + LADEE synthesis

### [ ] R5.8 — Validate meteoroid-driven water releases against LADEE NMS

**Why this step exists**

This tests the new H2O/OH model and time-dependent impact driver against discrete events rather than climatology.

**Capability/development work**

Run event windows using the frozen meteor driver and water physics; score detection, onset and amplitude without adding events after seeing NMS residuals.

**Primary code areas**

I39.

**Data required**

D15 published LADEE water-event list + D08 NMS + D11 meteor shower data.

**Where to obtain it**

Benna et al. 2019 DOI 10.1038/s41561-019-0345-3; LADEE PDS; IAU MDC.

**How to prepare/process the data**

Create a table of published event times/associations and uncertainties from the paper/supplement; obtain NMS H2O-related measurements from PDS; match established shower metadata from IAU MDC. Maintain an event/non-event evaluation window list fixed before model scoring.

**Tests to run at completion**

I39 event detection skill, onset timing and relative event response.

**Acceptance / exit criteria**

Detection skill >0.7; onset error <= one observation cadence; false positives and misses reported explicitly; no synthetic events.

**Expected result**

H2O/OH/meteoroid capability becomes validated for the tested event class.

**Dependencies / gate**

Requires R4.6–R4.7.

**Required completion artifacts**

- Updated production code with detailed comments explaining the physical assumption, units, coordinate frame, and failure behavior.
- Updated relevant `README.md` files describing configuration, data dependencies, and exact run command.
- Machine-readable `result.json`/comparison output from the linked test.
- Run manifest and SHA-256 hashes for all external data used.
- Test gates in this step:
  - `I39` — Meteoroid-driven exospheric water releases — LADEE NMS

### [ ] R5.9 — Validate pickup-ion transport where measurements support it

**Why this step exists**

Ion production/transport is a distinct capability; neutral validation does not validate pickup ions or Lorentz dynamics.

**Capability/development work**

Construct observation-space ion energy/mass/angle distributions from AMPS at Kaguya/PACE locations and compare only species/events with defensible instrument identification.

**Primary code areas**

MoonObservationGeometry.h; ion sampler; I40.

**Data required**

D13 Kaguya PACE PBF1 plus trajectory/field context and cited pickup-ion event papers.

**Where to obtain it**

DARTS PACE dataset darts:sln-l-pace-3-pbf1-v3.0.

**How to prepare/process the data**

Download event intervals, trajectory and magnetic/plasma context. Reproduce paper selection criteria for energy/mass/direction; do not invent species separation beyond published instrument capability. Convert AMPS particles through comparable energy/angle bins.

**Tests to run at completion**

I40 peak direction, energy, and flux where absolutely calibrated.

**Acceptance / exit criteria**

Peak angular bin error <=1 bin; peak energy relative error <=15%; flux within factor 2 when absolute calibration/source uncertainty supports such a test.

**Expected result**

Ion capability receives a separate evidence label from neutral species.

**Dependencies / gate**

Requires R3 chemistry + R4.1 plasma forcing.

**Required completion artifacts**

- Updated production code with detailed comments explaining the physical assumption, units, coordinate frame, and failure behavior.
- Updated relevant `README.md` files describing configuration, data dependencies, and exact run command.
- Machine-readable `result.json`/comparison output from the linked test.
- Run manifest and SHA-256 hashes for all external data used.
- Test gates in this step:
  - `I40` — Pickup-ion / ion escape comparison where observational constraints exist

### [ ] R5.10 — Freeze one cross-species parameter set and execute held-out validation

**Why this step exists**

A model that is re-tuned for every species/event is not a validated predictive framework.

**Capability/development work**

Write campaign/frozen_parameters.json identifying every free parameter, training data used, fixed literature values, priors and hold-out datasets. Lock it before final scoring.

**Primary code areas**

Campaign runner; I41.

**Data required**

Results/data from I32/I34/I36/I37 plus frozen parameter file; optionally I33/I35/I38–I40 for extended scorecard.

**Where to obtain it**

All preceding authoritative sources.

**How to prepare/process the data**

Generate frozen_parameters.json automatically from the exact run configuration and a human-reviewed calibration ledger. Hash it. Rerun held-out cases from clean output directories.

**Tests to run at completion**

I41 aggregated scorecard.

**Acceptance / exit criteria**

No post-hoc parameter changes after hold-out results are inspected. Any failed species remains failed/experimental rather than triggering silent retuning.

**Expected result**

A defensible cross-species validation statement is possible.

**Dependencies / gate**

Milestone M6.

**Required completion artifacts**

- Updated production code with detailed comments explaining the physical assumption, units, coordinate frame, and failure behavior.
- Updated relevant `README.md` files describing configuration, data dependencies, and exact run command.
- Machine-readable `result.json`/comparison output from the linked test.
- Run manifest and SHA-256 hashes for all external data used.
- Test gates in this step:
  - `I41` — Cross-species, cross-instrument validation without retuning

## Phase 6 — Uncertainty, Reproducibility, and Release Qualification

### [ ] R6.1 — Operationalize the one-command verification/validation runner

**Why this step exists**

A scientific validation campaign must be repeatable from a clean checkout and must not depend on manual test selection or interpretation.

**Capability/development work**

Finalize --list, --test, --phase, --all, --all-integration, --all-tests, --np, --nt, --launcher, --data-path and --strict-skips behavior; generate aggregate JSON/Markdown scorecards.

**Primary code areas**

test/run_tests.py; integration_common.py; integration_extended.py; READMEs.

**Data required**

No new external data.

**Where to obtain it**

Repository.

**How to prepare/process the data**

Run command matrix from a clean checkout; deliberately force one FAIL, one ERROR and one missing-data SKIPPED case to verify exit/status semantics.

**Tests to run at completion**

I42.

**Acceptance / exit criteria**

Statuses preserved exactly; strict-skips is nonzero when a required validation dataset/capability is missing; logs/manifests retained per test.

**Expected result**

Campaign can be run by another developer without tacit knowledge.

**Dependencies / gate**

Required for release automation.

**Required completion artifacts**

- Updated production code with detailed comments explaining the physical assumption, units, coordinate frame, and failure behavior.
- Updated relevant `README.md` files describing configuration, data dependencies, and exact run command.
- Machine-readable `result.json`/comparison output from the linked test.
- Run manifest and SHA-256 hashes for all external data used.
- Test gates in this step:
  - `I42` — Automated one-command Moon verification/validation runner

### [ ] R6.2 — Quantify parameter sensitivity, uncertainty and identifiability

**Why this step exists**

Agreement is not meaningful if equally plausible parameter combinations produce contradictory predictions or if numerical error exceeds physical effects.

**Capability/development work**

Define cited prior/range for uncertain source yields, binding energies, thermal parameters, rates and forcing; run local sensitivities then a reduced multivariate ensemble; propagate observation and Monte Carlo uncertainty into scores.

**Primary code areas**

Campaign parameter schema/post-processing; I43.

**Data required**

campaign/I43_sensitivity.csv generated from the model; observational uncertainty from I32–I40 data bundles.

**Where to obtain it**

Generated by the campaign, with priors cited to D03–D15/literature.

**How to prepare/process the data**

Use deterministic parameter grids/Latin hypercube with recorded seeds. Store every sample configuration and output metrics, not only summary plots.

**Tests to run at completion**

I43; rerun key observational metrics under parameter ensembles.

**Acceptance / exit criteria**

Dominant sensitivities physically interpretable; validation conclusions stable over credible parameter ranges; parameters claimed as constrained show sensitivity larger than relevant observation/noise floor.

**Expected result**

Quantitative uncertainty accompanies every validated observable.

**Dependencies / gate**

Requires frozen nominal model R5.10.

**Required completion artifacts**

- Updated production code with detailed comments explaining the physical assumption, units, coordinate frame, and failure behavior.
- Updated relevant `README.md` files describing configuration, data dependencies, and exact run command.
- Machine-readable `result.json`/comparison output from the linked test.
- Run manifest and SHA-256 hashes for all external data used.
- Test gates in this step:
  - `I43` — Sensitivity, uncertainty propagation, and identifiability

### [ ] R6.3 — Establish Monte Carlo convergence and seed robustness

**Why this step exists**

Sparse exosphere/tail/ion/PSR bins may look correct by chance in one particle realization.

**Capability/development work**

Expose/control random seed and particle statistics; run replicate ensembles at increasing particle count; compute standard errors and variance scaling.

**Primary code areas**

AMPS random initialization; runner; I44.

**Data required**

campaign/I44_seed_ensemble.csv generated from >=5 independent seeds per representative case.

**Where to obtain it**

Generated by AMPS.

**How to prepare/process the data**

Use identical physical input with independent seeds; record particle weights/counts. Compare integrated observables and bins above a predeclared effective-count threshold.

**Tests to run at completion**

I44 on representative Na, He, Ar and at least one low-density ion/PSR case.

**Acceptance / exit criteria**

>=5 seeds; key integrated observables <=5% 1-sigma Monte Carlo scatter; variance trends approximately N^-1/2 in well-behaved sampling regimes; validation classification stable.

**Expected result**

Statistical uncertainty is separated from model discrepancy.

**Dependencies / gate**

Required before final scores.

**Required completion artifacts**

- Updated production code with detailed comments explaining the physical assumption, units, coordinate frame, and failure behavior.
- Updated relevant `README.md` files describing configuration, data dependencies, and exact run command.
- Machine-readable `result.json`/comparison output from the linked test.
- Run manifest and SHA-256 hashes for all external data used.
- Test gates in this step:
  - `I44` — Monte Carlo statistical convergence and seed robustness

### [ ] R6.4 — Verify MPI/OpenMP reproducibility and production scaling

**Why this step exists**

Parallel scheduling can change stochastic streams/reductions or expose race conditions; high-resolution LOLA runs also need acceptable performance.

**Capability/development work**

Record rank/thread layout; use reproducible reductions where practical; otherwise require statistical equivalence. Add timing for mesh, mover, surface intersections, sources, chemistry and sampling.

**Primary code areas**

AMPS runtime + Moon diagnostics; I45.

**Data required**

campaign/I45_parallel_science_metrics.csv generated at 1x1, 4x1, 4x8 and a production layout.

**Where to obtain it**

Generated by AMPS.

**How to prepare/process the data**

Run the same tagged case/layout matrix on the same platform where possible. Compare science metrics with combined Monte Carlo uncertainty; preserve wall time and hardware/job metadata.

**Tests to run at completion**

I45 physical equivalence and scaling.

**Acceptance / exit criteria**

Integrated observables agree within combined 2-sigma statistical uncertainty; no budget/race errors; >20% performance regression versus approved baseline is flagged and explained.

**Expected result**

Parallel production runs are scientifically equivalent and operationally viable.

**Dependencies / gate**

Requires R6.3 uncertainty estimates.

**Required completion artifacts**

- Updated production code with detailed comments explaining the physical assumption, units, coordinate frame, and failure behavior.
- Updated relevant `README.md` files describing configuration, data dependencies, and exact run command.
- Machine-readable `result.json`/comparison output from the linked test.
- Run manifest and SHA-256 hashes for all external data used.
- Test gates in this step:
  - `I45` — MPI/OpenMP reproducibility and scaling

### [ ] R6.5 — Generate a release-readiness capability scorecard

**Why this step exists**

The application contains capabilities at different maturity levels; a single “validated model” label would overstate evidence.

**Capability/development work**

Map each capability to required U/I tests and evidence; automatically label VERIFIED, VALIDATED, EXPERIMENTAL or DISABLED; block release on P0 failures.

**Primary code areas**

README/release notes; I46 scorecard generator.

**Data required**

All test results + campaign/provenance_manifest.json.

**Where to obtain it**

Generated from the campaign and authoritative raw/processed data hashes.

**How to prepare/process the data**

Assemble a provenance manifest listing every raw product, processed product, script, DOI/URL, checksum, test result and capability mapping. Human-review the final scientific claims.

**Tests to run at completion**

I46 clean full-campaign scorecard.

**Acceptance / exit criteria**

Zero P0 failures; every VALIDATED label has all required independent observational tests PASS; missing/failed gates downgrade the capability rather than being waived silently; provenance complete.

**Expected result**

A release can state exactly what is verified/validated and under which conditions.

**Dependencies / gate**

Milestone M7 / definition of done.

**Required completion artifacts**

- Updated production code with detailed comments explaining the physical assumption, units, coordinate frame, and failure behavior.
- Updated relevant `README.md` files describing configuration, data dependencies, and exact run command.
- Machine-readable `result.json`/comparison output from the linked test.
- Run manifest and SHA-256 hashes for all external data used.
- Test gates in this step:
  - `I46` — Release-readiness science scorecard and capability levels

## 9. Full test map

### 8.1 Stand-alone lunar-module tests (no full AMPS link)

- `U01` — Line-of-sight geometry and output-record verification
- `U02` — Legacy-target isolation and source-level guard
- `U03` — Lunar gravity kernel and reference trajectories
- `U04` — Third-body gravity and rotating-frame kernels
- `U05` — Sodium radiation-pressure geometry and shadow kernels
- `U06` — Lorentz-acceleration kernel and species-scaling verification
- `U07` — Surface-refinement selector unit verification
- `U08` — LOLA data ingest, coordinate conversion, and topography geometry verification
- `U09` — Topographic illumination and horizon-shadow verification
- `U10` — Diviner temperature-product ingest and interpolation verification
- `U11` — Thermal-inertia and eclipse-response solver verification
- `U12` — Species sticking, residence-time, and re-emission distribution verification
- `U13` — PSR/cold-trap classification and trapping-kernel verification
- `U14` — Na photoionization-rate and survival-kernel verification
- `U15` — Multi-species photoionization table/rate verification
- `U16` — Electron-impact ionization coefficient verification
- `U17` — Charge-exchange reaction-kernel verification
- `U18` — Solar-wind/magnetospheric driver ingest and interpolation verification
- `U19` — Helium source/reservoir response-kernel verification
- `U20` — Neon source and accommodation-kernel verification
- `U21` — Radiogenic 40Ar source-map and transient-source verification
- `U22` — Na source-process kernel verification
- `U23` — Meteoroid forcing and impact-vaporization driver verification
- `U24` — H2O/OH surface-migration and photochemistry kernel verification
- `U25` — Observation/reference data provenance and immutable-package verification
- `U26` — Observation-operator geometry and instrument-coordinate verification

### 8.2 Linked AMPS integration and validation tests

- `I01` — Reproducible build, configuration, and provenance capture
- `I02` — Configured force, photochemistry, and surface dispatch audit
- `I03` — Species particle-weight and local-time-step initialization
- `I04` — Production sampler regression after LOS and record-format corrections
- `I05` — Lunar two-body gravity trajectory in the linked AMPS mover
- `I06` — Sun/Earth differential gravity and non-inertial terms
- `I07` — Sodium radiation pressure and shadow gating
- `I08` — Ion Lorentz-force mover in analytic uniform fields
- `I09` — Surface impact, re-emission path, physical escape and numerical-domain classification
- `I10` — Production particle-time-step convergence
- `I11` — AMR spatial-resolution convergence
- `I12` — Surface angular-refinement logic
- `I13` — Spherical-surface baseline preservation
- `I14` — LOLA topography integration into production surface boundary
- `I15` — Topographic illumination and shadowing
- `I16` — Diviner-driven surface-temperature boundary
- `I17` — Time-dependent thermal inertia and eclipse response
- `I18` — Species-specific sticking, residence time, and re-emission distributions
- `I19` — Surface inventory and global mass-conservation accounting
- `I20` — Permanent-shadow and polar cold-trapping physics
- `I21` — Sodium photoionization survival and Na→Na+ conversion
- `I22` — Photoionization for He, Ne, Ar and future H2O/OH
- `I23` — Electron-impact ionization rate coefficients
- `I24` — Charge-exchange reactions with solar-wind ions
- `I25` — Measured/time-dependent solar-wind and magnetospheric forcing
- `I26` — Helium source: solar-wind alpha implantation/re-emission plus endogenous component
- `I27` — Neon source and noncondensable surface accommodation
- `I28` — Radiogenic 40Ar source map and transient outgassing option
- `I29` — Sodium source partition: PSD, thermal desorption, sputtering, and impact vaporization
- `I30` — Time-dependent meteoroid forcing and impact vaporization
- `I31` — H2O/OH surface migration and photochemistry
- `I32` — Sodium limb/exosphere brightness — SELENE/Kaguya UPI-TVIS
- `I33` — Sodium extended tail and Earth-focusing validation
- `I34` — Helium local-time/altitude distribution — LADEE NMS
- `I35` — Helium solar-wind response — LRO/LAMP coordinated with ARTEMIS
- `I36` — Neon local-time distribution — LADEE NMS
- `I37` — Argon diurnal and selenographic structure — LADEE NMS
- `I38` — Argon sunrise pocket, lunation behavior, and polar sequestration — Apollo 17 LACE + LADEE synthesis
- `I39` — Meteoroid-driven exospheric water releases — LADEE NMS
- `I40` — Pickup-ion / ion escape comparison where observational constraints exist
- `I41` — Cross-species, cross-instrument validation without retuning
- `I42` — Automated one-command Moon verification/validation runner
- `I43` — Sensitivity, uncertainty propagation, and identifiability
- `I44` — Monte Carlo statistical convergence and seed robustness
- `I45` — MPI/OpenMP reproducibility and scaling
- `I46` — Release-readiness science scorecard and capability levels

### 8.3 Tests that require external validation/campaign data

- `I11` — campaign-generated mesh-convergence table
- `I13` — frozen spherical baseline metrics
- `I14` — LOLA LDEM production-surface package
- `I15` — independent LOLA/PGDA illumination/PSR reference
- `I16` — Diviner temperature products
- `I17` — time-resolved Diviner thermal histories
- `I21` — frozen model survival table plus authoritative Na photo-rate selection
- `I25` — ARTEMIS plasma-driver product
- `I30` — meteoroid-driver product
- `I32` — Kaguya TVIS Na images/geometry
- `I33` — ground-based sodium-tail observations
- `I34` — LADEE NMS He
- `I35` — LRO/LAMP He + ARTEMIS forcing
- `I36` — LADEE NMS Ne
- `I37` — LADEE NMS Ar + Kaguya GRS K proxy
- `I38` — Apollo 17 LACE + LADEE Ar
- `I39` — LADEE NMS water events + meteoroid metadata
- `I40` — Kaguya PACE pickup-ion intervals/context
- `I41` — frozen cross-species parameter file and preceding observational results
- `I43` — campaign sensitivity ensemble + observational uncertainties
- `I44` — multi-seed campaign ensemble
- `I45` — parallel-layout science/performance matrix
- `I46` — full provenance manifest and all prerequisite test results

## 10. Standard campaign commands

**List tests**

```bash
python3 test/run_tests.py --list
```

**Run all stand-alone tests**

```bash
python3 test/run_tests.py --all --output-dir test_output/U01-U26
```

**Run linked baseline**

```bash
python3 test/run_tests.py --phase baseline --amps ../amps --np 1 --output-dir test_output/baseline
```

**Run mesh/surface phase**

```bash
python3 test/run_tests.py --phase mesh-surface --amps ../amps --data-path /path/to/moon_validation_data --output-dir test_output/mesh-surface
```

**Run chemistry/source phase**

```bash
python3 test/run_tests.py --phase chemistry-sources --amps ../amps --data-path /path/to/moon_validation_data --output-dir test_output/chemistry-sources
```

**Run observation phase**

```bash
python3 test/run_tests.py --phase observations --amps ../amps --data-path /path/to/moon_validation_data --output-dir test_output/observations
```

**Run all linked I-tests**

```bash
python3 test/run_tests.py --all-integration --amps ../amps --data-path /path/to/moon_validation_data --np 4 --nt 8 --output-dir test_output/I01-I46
```

**Release gate — skipped validation is an error**

```bash
python3 test/run_tests.py --all-integration --amps ../amps --data-path /path/to/moon_validation_data --np 4 --nt 8 --strict-skips --output-dir test_output/release
```

**Run everything**

```bash
python3 test/run_tests.py --all-tests --amps ../amps --data-path /path/to/moon_validation_data --output-dir test_output/all
```

## 11. Release definition of done

- [ ] U01–U26 PASS from a clean checkout.
- [ ] I01–I46 are executed through the linked Moon executable for every capability claimed as validated; required gates are not replaced by preflight.
- [ ] Mesh, time-step, and particle-number errors are quantified and are smaller than the tolerances used for observational scoring.
- [ ] `SURFACE_GEOMETRY=LOLA` uses the actual LOLA-derived production particle-surface boundary; `SPHERE` remains available and regression-controlled.
- [ ] Terrain illumination, epoch geometry, and Diviner/dynamic thermal state refer to the same physical production surface.
- [ ] Surface inventories and global species/mass budgets close to the declared tolerance; numerical deletions are separated from physical sinks.
- [ ] Every enabled chemistry/source channel has an authoritative source/literature reference and both local and linked verification.
- [ ] Every capability labeled VALIDATED has at least one independent observational validation data set and an explicit calibration/hold-out split.
- [ ] Parameters are frozen before final hold-out scoring; a failed held-out case is reported rather than silently retuned.
- [ ] Numerical, observational, parameter, and Monte Carlo uncertainties are represented in the final scores.
- [ ] MPI/OpenMP layouts are physically equivalent within combined uncertainty and production performance is characterized.
- [ ] The release provenance manifest maps every capability to code revision, configuration, raw/processed data hashes, processing scripts, tests, and evidence.
- [ ] The final scorecard accurately distinguishes VERIFIED, VALIDATED, EXPERIMENTAL, and DISABLED capabilities.

## 12. Primary archive and literature references

- [LRO/LOLA — PDS Geosciences](https://pds-geosciences.wustl.edu/missions/lro/lola.htm)
- [LOLA Data Node](https://imbrium.mit.edu/)
- [NASA PGDA Lunar Polar Illumination](https://pgda.gsfc.nasa.gov/products/69)
- [LRO/Diviner — PDS Geosciences](https://pds-geosciences.wustl.edu/missions/lro/diviner.htm)
- [PHIDRATES](https://phidrates.space.swri.edu/)
- [Voronov/Verner collisional ionization](https://www.pa.uky.edu/~verner/col.html)
- [CDAWeb / ARTEMIS](https://cdaweb.gsfc.nasa.gov/)
- [Kaguya/SELENE DARTS mission archive](https://darts.isas.jaxa.jp/en/missions/kaguya)
- [Kaguya TVIS Na Level 2A](https://darts.isas.jaxa.jp/datasets/darts:sln-e-tvis-5-na-level2a-v1.0)
- [Kaguya GRS K/global intensity map](https://darts.isas.jaxa.jp/datasets/sln-l-grs-5-gamma-ray-map-v1.0)
- [Kaguya PACE PBF1](https://darts.isas.jaxa.jp/en/datasets/sln-l-pace-3-pbf1-v3.0)
- [LADEE NMS bundle](https://pds.nasa.gov/ds-view/pds/viewBundle.jsp?identifier=urn:nasa:pds:ladee_nms&version=2.0)
- [LRO/LAMP — PDS Imaging Node](https://pds-imaging.jpl.nasa.gov/)
- [IAU Meteor Data Center](https://iaumeteordatacenter.org/)
- [Apollo 17 LACE report NASA-CR-150946](https://ntrs.nasa.gov/citations/19760025001)
- [Wilson et al. lunar sodium tail, Icarus 204 (2009)](https://doi.org/10.1016/j.icarus.2009.06.017)
- [Benna et al. LADEE water events, Nature Geoscience 12 (2019)](https://doi.org/10.1038/s41561-019-0345-3)

- [OpenAI Codex prompting guide — `AGENTS.md` behavior](https://developers.openai.com/cookbook/examples/gpt-5/codex_prompting_guide)
- [OpenAI — Using `PLANS.md` for long-running Codex work](https://developers.openai.com/cookbook/articles/codex_exec_plans)

> **Archive-freshness rule:** before a campaign is frozen, open every product URL in this section and confirm that the product ID/version and filename still match this roadmap. Product-specific processing must always be driven by the metadata/label downloaded with the product. File dimensions, scale factors, quality flags, coordinate conventions, product versions, and archive paths may change; do not rely on a hard-coded value in this roadmap when the native label says otherwise.

## 13. Change-control rule for this roadmap

Any change to a physical assumption, observational data set, calibration/hold-out partition, or acceptance threshold must update this file, the affected test README, and the provenance/acceptance JSON together. Test gates must never be weakened silently to obtain a PASS.
