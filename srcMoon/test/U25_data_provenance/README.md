# U25 — Validation-data provenance and immutability

The executable contract is in [reference/acceptance.json](reference/acceptance.json):
purpose, procedure, acceptance criteria, failure modes, status semantics, and the
three decisions of 2026-11-23 (dataset selection, timestamp rule, C6 reporting).
The per-dataset gates and per-package recipes are in
[reference/packages.json](reference/packages.json). The verifier is
[package_verification.py](package_verification.py).

U25 judges validation-data **packages** only. A PASS means the selected package
is complete, unaltered, QA-ready, and reproducible from its immutable raw files
by its own committed processing script. It never compares AMPS with
observations, so it is never an observational validation.

## Running

From the AMPS repository root:

```bash
python3 srcMoon/test/U25_data_provenance/test.py                    # all declared datasets
python3 srcMoon/test/U25_data_provenance/test.py --dataset D04      # one dataset
python3 srcMoon/test/U25_data_provenance/test.py --for-test U15     # datasets declared by U15
```

Options: `--data-root` (default `$MOON_VALIDATION_DATA_ROOT` or
`/data/vtenishe/moon_validation_data`) and `--output-dir` (default
`test_output/srcMoon/unit`). `run_tests.py --all` and the nightly run the first
form. Exit codes: `0` PASS, `1` FAIL, `2` ERROR, `77` SKIPPED.

## Which datasets

The contract says the dataset list is "declared by the requesting test". The
requesting tests are the U contracts that list U25 among `predecessor_gates`.
The verifier reads their `required_data` at every run:

| Requesting test | Declared datasets |
|---|---|
| U15 | D04 |
| U16 | D05 |
| U18 | D06 |
| U21 | D10 |
| U23 | D11 |
| U26 | D07–D09, D12–D15 |

Without a target, U25 evaluates these 12 datasets. The U25 status is PASS only
if all 12 pass; otherwise it is the most serious of ERROR, FAIL and SKIPPED.
Each dataset's own verdict is in `U25/<Dxx>/result.json` and in the reason
line. A test that needs only its own package uses `--for-test` or `--dataset`.

## Gates, in order

Each gate stops the evaluation of a dataset; later gates are then reported as
not evaluated.

1. **Registry entry and package directory.** If absent: SKIPPED.
2. **QA readiness.** The registered key of `qa_report.json` must be `READY` or
   `PASS`. Any other native value: SKIPPED. A corrupt file or a missing key:
   ERROR.
3. **Frozen selection (observation datasets).** Roadmap C8 requires a
   calibration/hold-out split frozen before scoring. A selection counts as
   frozen only if its status is listed in that package's
   `frozen_status_values` in `packages.json`. No package had frozen its split
   on 2026-11-23, so the lists are empty and the gate fails closed: SKIPPED.
4. **Integrity.**
   - A missing required artifact gives SKIPPED.
   - These give ERROR:
     - a hash or size disagreeing with `SHA256SUMS.txt`, `raw_inventory.json` or `provenance.json`;
     - a file under `raw/` or `processed/` that no manifest lists;
     - a raw file not declared immutable, or one that differs from the inventory;
     - a provenance output or script list that does not cover its directory;
     - a symbolic link.
   - A provenance file declaring a prohibited method (manual edit, synthetic
     data, digitization, OCR) as used gives FAIL.
5. **Reproduction.**
   - The package's documented processing command runs on a private copy under
     the output directory. The canonical package is never written.
   - Network proxies point to a discard port, so an undeclared download fails.
   - These give ERROR: a nonzero exit, a timeout, or any change to the
     canonical package during the run (checked by re-hashing every file).
   - Any reproduced file that differs gives FAIL. The only exception is the
     timestamp rule below.

A dataset without a reviewed integrity/reproduction recipe (`verification: null`)
is SKIPPED if it ever reaches gate 4. It is never PASS.

## Timestamp rule

Rerunning a pipeline rewrites its provenance and QA files with a new generation
time. A reproduced file may differ from the package copy only on the lines of
fields that are both:

- frozen in the contract: `generated_utc`, `processing_utc`, or the text prefix
  `Generated UTC: `;
- registered for that file in `packages.json`.

Both values must be ISO-8601 UTC strings. A JSON file must also be equal to the
package copy after deleting exactly the registered key paths, so a field of the
same name elsewhere in the file cannot hide a change. The package's
`SHA256SUMS.txt` may then differ only in the entries of those files. Everything
else must be byte-identical.

## Metadata completeness

For each package, U25 reports which fields of the roadmap's C6 minimum
provenance schema are present, using the literal field names only. It does not
gate on them: no current package uses that schema literally. Byte-identical
reproduction is the direct evidence that no processed product was edited by
hand.

## Verifier controls

`test_package_verification.py` runs first on every U25 invocation and must
pass; otherwise U25 is ERROR. Its fixture packages, an integer list and its sum
in temporary directories, plant one defect each:

- a changed raw byte;
- an unlisted file;
- a hand edit with consistently updated hashes;
- a non-deterministic product;
- an unregistered time field;
- a declared manual edit;
- a failing or escaping processing script;
- a missing artifact;
- a symbolic link;
- PARTIAL or corrupt QA;
- an unfrozen or missing selection.

Each case checks the status that the contract requires. The controls also check
that every dataset declared by a requesting test is registered. The fixtures are
U25 test infrastructure, never validation data.

## Artifacts

Under `<output-dir>/U25/`:

- `result.json`: the aggregate.
- `<Dxx>/result.json`, `<Dxx>/verified_hashes.json`,
  `<Dxx>/reproduction_report.json`: the per-dataset results.
- `verified_hashes.json`, `reproduction_report.json`: indexes.
- `verifier-controls.log`.
- `packages/<UTC stamp>/<package>/`: the private copy and `reproduction.log`.

## State on 2026-11-23

- **D04 and D05 (photochemistry): PASS.**
  - 209 manifest entries, 178 raw files, 11 processed outputs and 6 scripts verified.
  - The rerun with the package's own `.venv/bin/python` (Python 3.8.10) took
    about 1 s. It reproduced 206 files byte for byte; `provenance.json`,
    `qa_report.json` and `qa_report.txt` were identical after the timestamp
    rule; `SHA256SUMS.txt` differed only in those three entries.
  - 14 of the 17 C6 fields are absent under their literal names.
- **D08: SKIPPED.** QA PASS, but `calibration_holdout.json` is
  `PROVISIONAL_NOT_FROZEN`.
- **D06, D07, D09–D15: SKIPPED** on their native QA status: ERROR, INCOMPLETE,
  `D10 DATA INCOMPLETE`, `D11 STREAM DATA PARTIAL`, or PARTIAL.

The untargeted U25 is therefore SKIPPED.

When one of these packages becomes ready, its integrity/reproduction recipe has
to be reviewed and added to `packages.json`. The package survey found several
obstacles to rerunning on a copy:

- `ladee_nms` refuses a separate output root.
- `lamp` and `lace` rewrite files under `raw/`.
- `pace` and `grs` write `.gz`/`.npz` outputs with embedded write times.
