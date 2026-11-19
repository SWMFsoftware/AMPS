# Native workbook and published-table source access

Preimplementation contract, 2026-11-19; source basis `8fbbd5e`.

The shared database currently accepts strict lexical CSV tables through explicit
source packets. Author-deposited Chang et al.2021 ZIPs contain five XLSX figure
workbooks. Numerical model/measurement columns and units are not uniformly
machine-labelled. A general source reader must therefore expose native cells
without inventing a column interpretation, normalization or missing observation.

`native_workbooks.py` will provide standard-library ZIP/XML decoding, an offline
`project_zip(source_zip, output_dir)` entry point and CLI. Every worksheet is
included in workbook order, including hidden worksheets. Cell coordinates,
explicit blanks, text `--`, raw numerical strings, shared/inline strings,
formula text and cached values are retained. Missing sparse cells remain absent;
rows are not padded. Formulae are never evaluated; dates remain serial numbers
plus source style identity. A long-cell CSV is source access, not a probability
sampler or a typed physical state/energy table. Native XLSX/ZIP bytes remain the
primary oracle and authoritative formatting/style metadata.

Errors include duplicate/escaping ZIP paths, external worksheet relationships,
unsupported namespace/type, malformed XML, duplicate coordinates or invalid
shared-string references. Existing output directories are never overwritten.
Archive expansion and XML sizes have explicit engineering limits; a refusal is
an infrastructure error, not a scientific discrepancy. No native source changes.

Protocol tests use clearly administrative ZIP/XML fixtures, never observations.
Real-source tests are marked `requires_authoritative_source=True` so the global
CORE software-only lane can omit data-dependent methods honestly. The original
archive/member XML is the independent lexical oracle. Tests cover all sheets,
blanks/sparse order, formulas without evaluation, shared/rich text, structural
corruption, traversal/external links, missing sources and deterministic output.
Actual Chang source tables will be independently checked against XML cells and
retained rights/producer checksums. JWST observed and model columns remain
separate; absent uncertainties remain absent. Poppe measured grid and model-rate
compilation remain distinct products. PASS denotes lexical/source integrity;
missing required source is SKIPPED; corruption/parser/IO failure is ERROR;
changed native values are FAIL. No U/I meaning or physical activation changes.

## API and replay

`decode_xlsx(raw_bytes, workbook_name)` returns `(rows, metadata)` without file
writes. `project_zip(source_zip, new_output_dir)` decodes every XLSX member and
emits one `workbook-NNN.cells.csv` per workbook plus `projection.json`.
`project_text(source, new_csv, skip_rows=..., source_width=...,
selected_columns=..., columns=...)` projects an **explicitly declared** native
whitespace-token table; it never guesses headers or skips unexpected blank rows.
This latter API separates published measurement and model columns by explicit
indices already documented in the original header.

The long-cell CSV fields are `workbook,sheet,sheet_state,cell,cell_type,
value_token,resolved_text,formula_token,formula_attributes_json,style_id`.
`value_token` retains original numeric strings, shared-string indices and formula
caches; `resolved_text` resolves display text without appending phonetic guides.
Coordinates and worksheet order locate every cell. Original ZIP/XML, styles,
merged-cell definitions and publisher captions remain available through the
source API; long cells do not imply a rectangular physical probability table.
The current supported namespace is the explicit transitional SpreadsheetML
namespace; unsupported namespaces/chartsheets are refused rather than omitted.

```sh
python3 -B utility/chemical_database/native_workbooks.py \
  --source-zip AUTHOR_DATA.zip --output-dir NEW_CELL_OUTPUT
python3 -B -m unittest discover -s test/chemical_database \
  -p test_native_workbooks.py -v
```

The 2026-11-19 retained source packets provide five Chang workbooks with 102,816
explicit cells; two JWST views each retain 2,851 native samples; the ARTEMIS
aggregate contains 1,600 rows and the separate Poppe model-rate compilation 25.
These counts are source-structure checks, not effective sample sizes. The JWST
source has no uncertainty columns. Native negative subtracted spectra, zeros,
missing cells and `--` are preserved. No selection, error bars, normalization,
physical rest-state reference, observation geometry or complete daughter PDF is
invented from those tables. Package preparation is reproducible by the retained
`prepare_packets.py` recipe, with the actual maintained reader copied/pinned in
each packet. Existing output directories require a new replay destination.
