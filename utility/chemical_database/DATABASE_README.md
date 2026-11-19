# AMPS shared chemical database

This installation is a portable, versioned chemical **source catalogue and
data-access library** shared by AMPS applications. It preserves original data,
native units, references, source conditions, missing fields and permission
holds. It supplies Python/CLI access and a dependency-light C++17 SDK.

The initial `catalogue-v1` profile supports metadata, source-object and native
table access. It does not select a physical reaction network, activate particle
chemistry, imply complete daughter distributions, or establish validation.
An integrity PASS means the requested bytes/schema were checked. Scientific
readiness and application qualification are separate.

The requested installation location is:

```text
/data/vtenishe/CHEMICAL_DATABASE
```

All internal runtime references are relative and hash-pinned. The directory
can be copied to another machine without the original AMPS checkout or
`moon_validation_data`. Use a Python 3 interpreter (tested with 3.8 and the
development interpreter) and a C++17 host compiler for native examples. Python
queries require only the standard library. Runtime uses no SQL server, network,
HDF5 library, PIC, MPI or SPICE. Original-format readers remain the caller's
responsibility for non-CSV sources.

## 1. Quick start

Run these commands from any working directory:

```bash
CHEM_DB_ROOT=/data/vtenishe/CHEMICAL_DATABASE
"$CHEM_DB_ROOT/bin/chemdb" --lock "$CHEM_DB_ROOT/database-lock.json" identity
"$CHEM_DB_ROOT/bin/chemdb" --lock "$CHEM_DB_ROOT/database-lock.json" verify
"$CHEM_DB_ROOT/bin/chemdb" --lock "$CHEM_DB_ROOT/database-lock.json" list --kind tables
"$CHEM_DB_ROOT/bin/chemdb" --lock "$CHEM_DB_ROOT/database-lock.json" query IAEA
```

`database-lock.json` selects exactly one release, profile and detached runtime
export by ID/path/SHA-256. It never means "latest". `identity` prints these
pins. `verify` checks stored source/member bytes and supported CSV schemas;
its output explicitly says that physical validation is NOT EVALUATED.

The first retained inventory contains 27 source groups, 622 source files and
35 table descriptors. Under the reviewed source policy, 31 source payloads
are importable; 591 source records are metadata-only. Seven table objects
are available, including four CSV representations and three native ASCII
tables. Exact installed counts are in `publication.json`; scientific source
issues and permission state are recorded per dataset, not hidden in totals.

To inspect chemical channel identities and source-labelled species:

```bash
"$CHEM_DB_ROOT/bin/chemdb" --lock "$CHEM_DB_ROOT/database-lock.json" list --kind species
"$CHEM_DB_ROOT/bin/chemdb" --lock "$CHEM_DB_ROOT/database-lock.json" list --kind reactions
"$CHEM_DB_ROOT/bin/chemdb" --lock "$CHEM_DB_ROOT/database-lock.json" describe phidrates-d04
```

Species entries are explicit source labels with unresolved state/isotope/mass
fields left missing. Reaction records preserve product occurrences, including
two separate H occurrences. Non-reactive absorption totals are labelled as
such. These metadata are not a compiled physical species/state mapping.

## 2. Read a real numerical table

The acquired IAEA D119990 partial-capture table has twelve original nodes:

```bash
CHEM_DB_ROOT=/data/vtenishe/CHEMICAL_DATABASE
"$CHEM_DB_ROOT/bin/chemdb" --lock "$CHEM_DB_ROOT/database-lock.json" \
  describe collisiondb_119990_embedded_arrays.csv
"$CHEM_DB_ROOT/bin/chemdb" --lock "$CHEM_DB_ROOT/database-lock.json" \
  table collisiondb_119990_embedded_arrays.csv --csv
```

The native header is:

```text
energy_eV_per_u,cross_section_cm2
```

The first source row is `4000.0,6.34e-15`. The API retains the original lexical
values and order. The energy unit is **eV per projectile atomic-mass unit**;
the native energy is in the **target-rest frame** (`Frame: Target`).
The cross-section unit is **cm²**, not m². The represented process is the
H(1s) ground-product partial capture channel, distinct from a total CEX law.
The source does not establish the initial Na electronic state. Reading these
nodes supplies neither interpolation, a daughter angle distribution nor an
ordinary solar-wind model. Other D107581 records include fit/plot/composite
representations and must not be combined with this measured partial as if
they were additional independent pure-capture measurements.

Python accepts a stable complete table ID or an **exact unique native table
name**. Ambiguous names fail; there is no silent first/latest-source choice.
C++ uses complete stable IDs. `describe` reports the actual `table_id`.

## 3. Python API

No pip installation is needed. Point `PYTHONPATH` at this installation:

```bash
CHEM_DB_ROOT=/data/vtenishe/CHEMICAL_DATABASE
export PYTHONPATH="$CHEM_DB_ROOT/python${PYTHONPATH:+:$PYTHONPATH}"
python3 -B - "$CHEM_DB_ROOT/database-lock.json" <<'PY'
import sys
from amps_chemical_database import Database, DatabaseError

database = Database.from_lock(sys.argv[1])
print("Pinned release:", database.identity["release_id"])
print("Source groups:", [item["id"] for item in database.datasets()])

table = database.table("collisiondb_119990_embedded_arrays.csv")
print("Columns:", table.columns)
print("Number of native rows:", table.row_count)
print("First lexical row:", table.rows[0])
energy = table.number(0, "energy_eV_per_u", expected_unit="eV/u")
cross_section = table.number(0, "cross_section_cm2", expected_unit="cm^2")
print("Exact Decimal energy and cross section:", energy, cross_section)

for row in table.dictionaries():
    print(row["energy_eV_per_u"], row["cross_section_cm2"])
PY
```

The public interface is:

| Operation | Result/contract |
|---|---|
| `Database.from_lock(path)` | Open the exact root-local release/profile/export; refuse malformed or conflicting identity |
| `Database.open(root, release_id, expected_manifest_sha256, profile_id, expected_profile_sha256, export_path, expected_export_sha256)` | Explicit alternative to a lock; every pin is required |
| `identity` | Read-only export identity, root-relative members and digests |
| `datasets()`, `sources()`, `tables()` | Immutable metadata views; tables include native names and documented units |
| `species()`, `reactions()` | Source identity/topology metadata, preserving unresolved physical fields |
| `describe(identifier)` | Dataset/source/table metadata, provenance, issues and permissions |
| `query(text)` | Case-insensitive metadata search; never ranks physical models |
| `require(capability)` | Check the selected profile's declared scope; no physical completeness inferred |
| `source_path(source_id)` | Resolve/authenticate a permitted source object inside this archive |
| `table(table_id_or_exact_name)` | Owned `Table` of exact native CSV strings, with shape/hash/permission checks |
| `verify(deep=True)` | Verify release members and source/table bytes; no observational score |

`Table.columns` is a tuple; `Table.rows` is a tuple of row tuples. Header
padding, signs, source precision and empty cells survive. `value(row,column)`
returns a string. `number(row,column,expected_unit)` returns a `Decimal` only
when that column's documented unit matches; missing/qualified/nonfinite cells
or unknown units do not become a number. It performs no unit conversion or
positivity repair. Physical rate evaluators must enforce their own domains.
`metadata` and records returned by query methods are defensive immutable views.

For JSON serialization, use the CLI or:

```python
from amps_chemical_database.api import plain
import json
print(json.dumps(plain(database.describe("phidrates-d04")), indent=2))
```

## 4. Native ASCII and original documents

ExoMol `.photo` files and CODATA fixed-width text retain their native format.
They are not automatically treated as a CSV or given invented headers.
Use their source IDs and your documented native-format decoder:

```bash
CHEM_DB_ROOT=/data/vtenishe/CHEMICAL_DATABASE
"$CHEM_DB_ROOT/bin/chemdb" --lock "$CHEM_DB_ROOT/database-lock.json" source-path \
  nist-codata-2022-unmarked:cd94f444f2c77cb3eec7603974f35b5b7319352f818fe0cfa455788e1c994e5a
"$CHEM_DB_ROOT/bin/chemdb" --lock "$CHEM_DB_ROOT/database-lock.json" source-path \
  molecular-exomol-research:66929c3806178e46a9a56b50397f0dd69d88397bc56970a5aa214ac78270159b
```

The second object is the acquired OH native spectrum. Its 19,173 actual rows
differ from the 19,172 declared in its metadata. The discrepancy remains an
issue; the archive does not remove a row to enforce agreement. A source-path
lookup is permission-cleared byte access, not a statement that absorption
data are usable daughter kinetics. Historical original paths remain only in
provenance; runtime never follows them to recover a missing archive object.

## 5. C++17 API and compilation

The SDK is under `sdk/src/models/chemical_database`. Public headers:

| Header | Interface |
|---|---|
| `Runtime.h` | `BundleRequest`, typed `Status`, immutable owning `Catalogue`, source/table records |
| `DataIO.h` | SHA-256, strict UTF-8/CSV and bounded-file identity helpers |
| `Photochemistry.h` | Source-authenticated D04 registry with explicit quiet/active cases |
| `ElectronImpact.h` | Native Voronov table and Maxwellian coefficient service |
| `ChargeExchange.h` | Existing named total-capture services with source/domain guards |
| `RateServices.h` | Catalogue-bound preparation factories requiring explicit table roles, permissions and scalar-rate capability |

The base reader has no application dependency. Compile its supplied inspection
example in a writable working directory (the archive itself remains immutable):

```bash
CHEM_DB_ROOT=/data/vtenishe/CHEMICAL_DATABASE
export PYTHONPATH="$CHEM_DB_ROOT/python${PYTHONPATH:+:$PYTHONPATH}"
c++ -std=c++17 -O2 -pthread \
  -I "$CHEM_DB_ROOT/sdk/src/models/chemical_database" \
  "$CHEM_DB_ROOT/sdk/src/models/chemical_database/DataIO.cpp" \
  "$CHEM_DB_ROOT/sdk/src/models/chemical_database/Runtime.cpp" \
  "$CHEM_DB_ROOT/sdk/src/models/chemical_database/example.cpp" \
  -o ./chemdb-inspect
python3 -B - "$CHEM_DB_ROOT/database-lock.json" ./chemdb-inspect <<'PY'
import sys
import subprocess
from pathlib import Path
from amps_chemical_database import Database

database = Database.from_lock(sys.argv[1])
pin = database.identity
subprocess.run([str(Path(sys.argv[2]).resolve()), str(database.root),
    pin["export_path"], pin["export_sha256"], pin["release_id"],
    pin["release_manifest_sha256"], pin["profile_id"], pin["profile_sha256"],
    "catalogue"], check=True)
PY
```

This executes the actual native reader. Expect `Evaluated` and the selected
release/profile/source/table counts. The sample takes eight explicit arguments;
the Python launcher extracts their values from the same pinned lock.

To use the API inside an application, prepare it once before transport:

```cpp
#include "Runtime.h"
namespace Chemical = AMPS::ChemicalDatabase;

Chemical::BundleRequest request;
// Supply these from the application's parsed, preserved configuration:
// root, exportPath, expectedExportSha256, expectedReleaseId,
// expectedReleaseSha256, expectedProfileId, expectedProfileSha256.
request.requiredCapabilities = {"table-access"};
Chemical::Catalogue catalogue;
Chemical::Status status = Chemical::OpenRuntimeBundle(request, &catalogue);
if (!status.ok()) {
  // Report CodeName(status.code) and status.message; do not substitute zero.
}
std::shared_ptr<const Chemical::ReactionData::Csv> rows;
// tableId is the complete stable ID returned by GetTable/tables metadata.
// status = catalogue.TableRows(tableId, "table-access", &rows);
```

`TableRows` contains data rows only; column names are in `TableRecord.columns`.
The native reader verifies and prepares eligible CSVs at open, then provides
immutable cached rows. A returned `shared_ptr` owns its rows beyond the lifetime
of a catalogue handle. Catalogue copies share const storage. Failed opens and
lookups preserve the caller's previous output. Concurrent const host queries
are supported; mutating/reloading the same handle concurrently is not.

For model evaluation link the additional `Photochemistry.cpp`,
`ElectronImpact.cpp`, `ChargeExchange.cpp` and `RateServices.cpp` from the SDK.
Their formulas, frozen-source identities, units and domains were moved from
the existing production services without changing the scientific selection.
They do not infer a spectrum, isotope mix, electron distribution or event law.
Factories require explicit table IDs, correct schema roles, `scalar-rates`
scope and cleared source dependencies. `catalogue-v1` does **not** grant that
scope. This prevents a lexical data archive from becoming an implicitly
selected chemistry model. Direct legacy source readers retain their original
explicit-source contracts; they do not authorize restricted managed storage.

`ResourceLimits.maximumFileBytes` and `maximumTotalBytes` are explicit optional
reader ceilings; zero means no artificial ceiling. Representability and
allocation failures still return a diagnostic. Applications map fatal required
failures to their normal termination handler, e.g.
`exit(__LINE__,__FILE__,message)`, with instructions to increase available
resources or reduce the explicitly selected data. No truncation occurs.

V1 is host code. No GPU chemistry/device-query qualification is claimed.
AMPS keeps the normal Config.pl/Make build; the shared sources have unique
`ChemDb*.o` archive members. Existing Moon namespace aliases call those same
implementations. This release does not add a Moon/Hartley2 `.in` switch that
silently replaces their current physical source selections. Other applications
can use this documented API through their own explicit configuration adapter.

## 6. Errors, missing data and permissions

`DatabaseError.code` and C++ `Code` distinguish:

```text
Evaluated, Unavailable, InvalidInput, OutsideDomain, CorruptSource,
UnsupportedSchema, PermissionRequired, NumericalFailure, ResourceLimit
```

The CLI exits 0 on successful access/integrity checks, 77 for an unavailable
or permission-held prerequisite, and 2 for invalid/corrupt/infrastructure
operations. There is no implicit fallback to another source, numerical zero,
ground state or repaired table. A wrong format request is an explicit unsupported
schema operation; use source-path plus its proper decoder for native ASCII.

Storage, evaluation, redistribution, source interpretation and application
qualification are separate fields. Read `RIGHTS.md` and each dataset record.
Current policy admits the reviewed IAEA, ExoMol, unmarked NIST constants/isotope
and selected licensed deposit sources with their stated conditions. ADAS needs
permission for insertion into a managed database. Several other acquired
sources have unresolved permission scope for this managed operation and are
indexed with source identities and issues, without copying their payloads.
No provider was contacted and no permission was manufactured by installation.

Existing PHIDRATES/Voronov originals and legacy application use remain in their
original archive. Holding a new managed copy does not disable those old modes.
It does mean that this shared release cannot claim numerical availability or
scalar/event qualification for those held catalogue entries. Readiness gaps in
states, thresholds, uncertainties, joint daughters, cascades/quenching and
material outcomes remain explicit even for a permission-cleared source.

## 7. Layout and provenance

```text
README.md, RIGHTS.md, publication.json, database-lock.json
bin/chemdb                            portable Python CLI
python/amps_chemical_database/        Python API, CLI and publisher source
sdk/src/models/chemical_database/     C++17 source/headers/example/schema
sdk/src/general/input/amps_input.h    exact shared AMPS grammar dependency
objects/sha256/<digest>/payload       exact permitted original bytes
releases/<release-id>/                canonical JSON catalogue/profile/provenance
exports/<release-sha>/<profile-sha>/csv-block-v1/
  export.in                          detached pinned runtime manifest
  source_index.csv, tables_index.csv  exact runtime record projections
SHA256SUMS                           published data/export identities
```

The catalogue records each original source location, URL, observed/expected
hash, native filename, permission evidence, original table descriptor and
scientific issue. Empty/missing fields and unsupported native schemas are
not silently repaired. Unrepresentable held CSV headers remain documentary
metadata with an opaque runtime descriptor. No source number is converted,
clipped, sorted, smoothed or renormalized by this publication process.

Exact acquisition receipts retain their recorded request/final URL, retrieval
time and HTTP outcome where present, together with the original manifest hash
and JSON locator. All matching receipts and discrepancies survive; a missing
receipt remains missing. A package creation time is not substituted for a
file retrieval time. Runtime uses no historical original path for recovery.

The four imported IAEA CSV representations are independently reproduced from
the hash-bound native typed arrays during publication. Their table metadata
records the source-object dependency, recipe hash, target frame, renderer role
and exact output digest. `Patch` geometry and the fitted composite are retained
with their actual meaning. Reproduce D119990 independently into a new output:

```bash
CHEM_DB_ROOT=/data/vtenishe/CHEMICAL_DATABASE
export PYTHONPATH="$CHEM_DB_ROOT/python${PYTHONPATH:+:$PYTHONPATH}"
python3 -B - "$CHEM_DB_ROOT/database-lock.json" /tmp/NEW_IAEA_REPRODUCTION <<'PY'
import sys
from pathlib import Path
from amps_chemical_database import Database
from amps_chemical_database.native_iaea import reproduce

database = Database.from_lock(sys.argv[1])
metadata = database.describe("collisiondb_119990_embedded_arrays.csv")
source = database.source_path(metadata["processing"]["source_id"])
records = reproduce(source, Path(sys.argv[2]), dataset_id=119990)
print(records)
PY
```

This performs source extraction only and refuses an existing output directory.
The publisher/reader/source-policy/recipe software is retained under the
release's provenance as well as in the installed package. Original native
arrays and CSVs remain byte-identical to their acquired versions.

The hash graph is acyclic: canonical members and profile precede release;
the detached export cites the finished release and its runtime indexes;
the installation lock pins the finished export. There is no manifest hashing
itself or cross-reference cycle. Import recipes/source revision are preserved.
An explicitly stored lock is a selection identity, not a scientific approval.

## 8. Relocation, rebuilding and maintenance

To relocate, copy the **whole installation** to a new empty directory, then
run its CLI against its own `database-lock.json`. No original checkout path
is used for runtime data recovery. Preserve attribution, native licence
conditions and the redistribution holds described in RIGHTS.md; local copy
permission does not automatically authorize public distribution of every item.

Do not edit objects, release files, indexes or lock hashes in place. A new
source/version needs a new publication and reproducible source-policy review.
The maintained publisher is available in AMPS:

```bash
python3 -B utility/chemical_database/build_database.py \
  --amps-root /path/to/AMPS \
  --validation-root /path/to/moon_validation_data \
  --output-root /path/to/NEW_EMPTY_PUBLICATION \
  --release-id NEW_RELEASE_ID
```

`output-root` must not already exist. The source inventory explicitly accounts
for unavailable packets. Original bytes remain read-only; partial publication
artifacts survive errors for diagnosis and are never mistaken for a complete
release. Rebuilding from another source checkout requires its corresponding
retained source packets; runtime access to an existing installation does not.

The shared verification suite is maintained at
`test/chemical_database/run_all_tests.py` in AMPS. It covers the actual Python/
C++ implementation, corruption/path/schema/permission refusal, immutable row
lifetime, exact native-source agreement and relocation. The source-independent
protocol is registered in the global core routine. Full application regression
results and resource/layout evidence remain in the dated AMPS test outputs.
No source/integrity/software PASS should be reported as physical validation.

The original source archives have been preserved. Backup scheduling and the
long-term maintenance/distribution owner are not silently invented by this
installation; retain the immutable snapshot and original archives until those
operational arrangements are established.
