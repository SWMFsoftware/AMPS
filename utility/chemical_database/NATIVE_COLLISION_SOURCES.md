# Native collision-source reader

`native_collision_sources.py` provides a strict, standard-library-only lexical
reader for native IAEA CollisionDB HTML fit-coefficient tables. It preserves
the published decimal spellings and row order. It does not evaluate a rate,
change a unit, choose a species state, infer a missing range, run the displayed
source Python, or construct a product distribution. Permission is checked
separately by the shared source-packet importer.

The audited 2026-11-19 acquisition contains 21 IAEA Na/He electron-ionization
pages and 172 printed coefficient rows. Native download/API endpoints returned
HTTP500 and are retained as failed acquisitions. The reader does not turn those
error responses into data. A source-page table is not a replacement for a
missing full-precision export. Current IAEA differential searches found three
records, all BeO, and none for the requested Na/He collisions. This is a report
of the captured queries, not an assertion that no differential data exist
elsewhere. No full joint collision-product distribution was obtained.

`parse_iaea_fit_page(html, expected_dataset_id)` returns exact source parameters,
rendered reaction, native columns, threshold, uncertainty, method, comment,
printed domain and fit-function links. Missing formula/domain remains explicit.
`parameters_csv(parsed)` emits only a deterministic lexical projection with
parser-assigned header `parameter_label,value_native`; dimensions remain
parameter-specific and must come from the original expression. `project(input,
sha256,output,dataset_id)` authenticates input bytes and refuses an existing
output. `NativeCollisionSourceError` means invalid/missing/inconsistent native
metadata. No fallback number is returned. Ordinary missing-file/output
collisions remain filesystem errors.

Reproduce one acquired table, using its actual packet SHA-256:

```sh
python3 -B utility/chemical_database/native_collision_sources.py \
  --input /path/to/packet/raw/iaea-D111150.html \
  --expected-sha256 <SHA256-from-source_packet.json> \
  --dataset-id D111150 --output /new/directory/D111150-fit-parameters.csv
```

The packet includes an exact copy/hash of this recipe and original native-parent
IDs. Its IAEA CC BY4.0 evidence controls managed storage. Neither a successful
parse nor an allowed licence qualifies a model. CCC data and OPEN-ADAS retain
unknown/pending managed-storage policies respectively and must not be copied
through an allowed derived-data record.

The native pages have unresolved source issues: Na pages can omit the rendered
formula/upper validity limit; coefficient tables show rounded precision;
JBORN1 Python/Fortran disagree about the minimum applicable energy; ALBELI's
displayed Python has an unassigned `ion` path and an inclusive array-index issue.
These are retained, not silently corrected or installed as physical kernels.
A source-referenced, independently verified evaluator needs a separate contract.

Tests are administrative/parser verification, not synthetic observations:

```sh
python3 -B -m unittest discover -s test/chemical_database \
  -p test_native_collision_sources.py -v
```

They exercise exact lexemes/order, declared missingness, identity mismatch,
malformed/duplicate/nonfinite rows, HTTP-error-shaped input, SHA refusal,
deterministic reproduction and refusal to overwrite source-derived outputs.
Units, frames, stochastic controls, independent observational reference and
formal U/I predecessors are inapplicable to these lexical fixtures. Every
source record retains its native physical conventions independently. Existing
U16/U17/I23/I24 meanings, thresholds and qualification order are unchanged.
