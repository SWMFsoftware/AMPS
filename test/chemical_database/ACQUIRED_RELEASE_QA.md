# Explicit real acquired-release audit

`verify_acquired_release.py` is an explicit real-data audit, excluded from
unittest discovery. It opens a supplied immutable release through the maintained
Python API. It compares actual stored data with original acquired source bytes,
not synthetic observations or AMPS output. A missing original is SKIPPED; an
evaluated preservation mismatch is FAIL; parser/hash/runtime corruption is ERROR.
The audit never edits an original, chooses a chemistry law or activates physics.

```sh
python3 -B test/chemical_database/verify_acquired_release.py \
  --database-lock /path/to/prepared/database-lock.json \
  --output-dir /new/audit/directory \
  --native-probe /path/to/already-compiled/catalogue-probe
```

The optional native probe is the existing `test/chemical_database/probe.cpp`
executable. No compilation/configuration occurs in this runner. With that option,
all45 new CSV tables are independently read through the actual C++ Catalogue,
compared against Python and the original oracle, and retained with exact argv,
executable hash, stdout/stderr and response hashes under `native-queries/`.
Without it, the result explicitly states that native comparisons were not requested.

The source packet declarations supply explicit original paths and SHA256 pins;
there is no latest-data discovery. All six explicit declarations must also be
archived byte-for-byte in the immutable release and listed in release.members.
Additional datasets such as ExoAtom are accepted with their original packet,
byte and quantity metadata checked; no inferred atomic state/model is selected.

The data-specific checks cover21 IAEA native fit tables/172 unchanged coefficient
lexemes and their source domains; NASA19 species/56 interval records with native
fixed-column source lines and Fortran-D coefficients; ExoMol7933x2 photo rows,
10606 OH states,10000 partition rows and every original negative/zero sample;
and18 original ZIP-member surface tables, with four measured tables/17 rows
separated from14 model tables/2033 rows. All held sources/tables must deny actual
payload operations. Labels, frames, units, source zeros, missing values and
uncertainty conventions remain native; no conversion, interpolation, clipping,
normalization, isotope/Q crosswalk or new statistical threshold is introduced.

Artifacts are `result.json`, `audit.log`, any `exception.log`, and per-query native
receipts. They include exact lock/export/release identities, script/source revision,
original input hashes, outcomes/counts, numerical source metrics and actual command.
Output directories must be new, preserving earlier failures or missingness.
A source-integrity PASS is not U/I model qualification or observational validation.
This extends shared database verification without assigning a new lunar test ID.
