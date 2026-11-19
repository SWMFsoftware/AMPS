# Real physical-source/model release QA

`verify_physical_model_release.py` is an explicit real-data audit, outside unittest
discovery. It opens a supplied immutable database lock with production Python
`Database`, then uses the actual compiled C++ `Catalogue` probe for every allowed,
available CSV table. No rate, surface law, spectrum fit or model formula is
implemented by this audit.

```sh
python3 -B test/chemical_database/verify_physical_model_release.py \
  --database-lock /path/to/database-lock.json \
  --output-dir /path/to/NEW-QA-DIRECTORY \
  --native-probe /path/to/compiled/catalogue-probe \
  --source-packet /path/to/original/source_packet.json \
  --source-packet /path/to/another/original/source_packet.json \
  --wall-seconds 600
```

Use `test/chemical_database/probe.cpp` compiled with the actual maintained
DataIO/Runtime sources. Its new `denial-sweep` operation calls the production
permission gate for **every** held source and table, checking that refused output
arguments remain unchanged. It reads one authenticated Catalogue; this avoids
hundreds of redundant opens without sampling IDs or duplicating permission logic.
The audit compares its exact ID membership and typed result with the Python gates.

The independent lexical oracle is original source table bytes where available,
and otherwise the exact release-pinned CSV, explicitly distinguished in metrics.
Supply every admitted explicit packet as a repeated `--source-packet`; a missing
argument is SKIPPED and an unselected extra packet is ERROR. The admitted set is
recorded in release `provenance/import.json`. These arguments authenticate every declared member, native
parent and recipe against the admitted metadata. All archived original packets
and declared recipe bytes must match both originals and release-member hashes.
The model descriptor list, individual retrieval and published CSV JSON fields
must agree; metadata remains immutable and `physical_activation=false`.

All numbers remain strings in their native units. Negative observed spectra,
zero bins, missing cells and formula caches are neither repaired nor normalized.
Mixed observed/model roles are retained from each packet's source metadata.
The tool does not turn spectral or aggregate grid rows into independent samples,
invent uncertainties, establish reaction/Q/state mappings or interpret missing
coverage. Source integrity and software accessibility do not constitute U/I
physical verification or observational validation.

Outputs are `result.json`, `audit.log`, and individual command/stdout/stderr/hash
receipts under `native-queries/`. Existing output directories are refused. A
native executable is necessary for the declared scope: absence yields SKIPPED,
not a Python-only PASS. Missing originals/capabilities are SKIPPED; evaluated
lexical/access mismatches are FAIL; parser, hash, runtime, signal and budget
failures are ERROR. Exit codes are 0,1,2,77 respectively. Actual passing subchecks
remain visible if a later prerequisite is unavailable. The final step rehashes
originals, lock and executable to detect mutation during the audit.
