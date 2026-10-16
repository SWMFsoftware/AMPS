# M2 validation-data qualification helpers

These scripts inspect authoritative files under the canonical external data
root without changing anything under `raw/`. The raw-package qualifier updates
the LOLA external package's reproducible provenance/QA records; the D02/D03
processors write derived science tables only into the explicitly supplied
output directory. D02/GCP/PRP examples use campaign-owned directories under
`test_output/srcMoon/`; the I17 RDR processor is the deliberate exception and
is run explicitly into the canonical `diviner/processed/i17` directory so a
linked campaign can stage byte-identical evidence without reprocessing 11 GB
of raw data. A processor-level
`status=READY` means that transformation and its hashes completed; it never
promotes the external Dxx package. Linked tests separately require the package
`README.md`, `provenance.json`, and `qa_report.json` and an explicit, consistent
`READY` dataset status.

Run the LOLA D01/D02 qualification from the AMPS root:

```bash
python3 -B srcMoon/test/data_processing/qualify_lola_package.py
```

The script freezes and verifies every expected raw SHA-256 digest, derives
each binary byte-count contract from its own PDS detached label, and reports
dataset readiness separately. D01 can become `READY` for U08/I14 input use.
D02 remains `PARTIAL` after the raw-package check alone.  Create the independent
spatial table in a campaign-owned output directory with:

```bash
python3 -B srcMoon/test/data_processing/prepare_d02_reference.py \
  --output-dir test_output/srcMoon/M2-data/D02
```

The processor verifies the eight frozen raw hashes, decodes the PDS scale and
offset, applies the label-native polar-stereographic inverse, and selects a
deterministic 100 km comparison lattice.  It writes
`illumination_reference.csv` and `d02_processed_qa.json`.  The QA record keeps
PSR-edge rows identifiable and records the independent product's extended-Sun,
18.6-year, one-hour temporal operator.  Generating this table qualifies the
reference data; it does not establish I15/I20 PASS until the linked production
operator is scored.

Prepare the D03 global cumulative product used by I16 with:

```bash
python3 -B srcMoon/test/data_processing/prepare_d03_gcp.py \
  --output-dir test_output/srcMoon/M2-data/D03
```

This streamed conversion requires all 18 immutable GCP TAB/XML pairs.  It
selects only label-defined bolometric brightness temperature (`tbol`), rejects
the native missing constant, and uses a predeclared checkerboard split: coarse
nodes supply `temperature_table.csv`, while spatial midpoints at a fixed
six-hour cadence supply `benchmark_points.csv`.  The two-decimal local-time
centres alternate by 0.24/0.26 hour in text; the script validates them against
the exact quarter-hour bin centres instead of treating the rounded strings as
an exact uniform axis.  D03 supplies no per-bin GCP uncertainty, so the output
field remains empty and `d03_gcp_qa.json` records that fact.  Empty never means
zero uncertainty.

Prepare the independent D03 polar-resource rows for I20 with:

```bash
python3 -B srcMoon/test/data_processing/prepare_d03_prp.py \
  --output-dir test_output/srcMoon/M2-data/D03-PRP
```

The script retains the nearest real PRP triangle for every point on a
predeclared approximately 100 km polar lattice.  It does not regrid or average
temperature.  `polar_temperature_reference.csv` keeps annual-average
temperature at 2 cm depth distinct from annual-maximum surface temperature,
and its QA report records the absence of per-triangle uncertainty.

Prepare the frozen, time-resolved D03 RDR product used by I17 with:

```bash
python3 -B srcMoon/test/data_processing/prepare_d03_rdr_eclipse.py \
  --output-dir /data/vtenishe/moon_validation_data/diviner/processed/i17
```

This is an explicit data-maintenance operation and is never invoked implicitly
by `run_all_test.py`.  It applies the pre-hashed I17 selection contract to the
June training and December hold-out dates identically.  Each channel retains
its actual UTC, footprint centre, detector identities, viewing angle, and
native qca/qge/qmi values.  Calibrated radiance is averaged in fixed
600-second/0.5-degree bins before the native radiance-to-temperature inverse is
applied; gaps and incomplete C3--C9 bins are rejected rather than filled.

The published Paige C3--C9 sensitivity substitution and wavelength partitions
produce `bolometric_brightness_temperature_K`.  Calibration and within-bin
radiance errors are combined before conversion, then propagated through T_bol
with a branch-preserving local Jacobian.  Cross-channel covariance is not
available and is recorded as a limitation, not guessed.  The processor never
reads AMPS output, writes the science CSV under a temporary name, publishes it
atomically, and publishes the matching QA receipt last.

The 2026-10-13 reproducibility run produced 6,353 complete bins from 87
archives.  Its history, frozen-contract, and processor SHA-256 values are
`870577c82b430b098d816d210b528c714ec6cfa3a20ae9be3f6482a3639ae182`,
`88fc5196fee5863691937a505873d420990b98056db8a1e8167255cd736189b0`, and
`5a5d29d0f6dba88cd49c899e541a980eb84ce24717481b0cad471f8012459a65`.
This processor receipt does not make the D03 package or I17 observational
comparison READY; see the I17 README for predecessor and scoring gates.

The default root is `/data/vtenishe/moon_validation_data`; `--data-root` may
select an equivalent package for an isolated reproducibility check. The exact
authoritative source URLs and the native body-fixed frame are written to
`lola/provenance.json`.

The source-tree roadmap (SHA-256
`8a9a7f16cd929b4848da869ba9e41de0368e0633a6e444d3d4b38a96d2fb068a`)
contains the current implementation status and
controls code/test decisions.  The external-data roadmap (SHA-256
`d10eca1c539f9cf82940236d9e1aa36af665d43e82ba475925950d582fbf49c4`)
is the processing contract for immutable raw data.  They differ because the
source copy contains newer implementation-status entries; no scientific
requirement is silently reconciled or discarded.
