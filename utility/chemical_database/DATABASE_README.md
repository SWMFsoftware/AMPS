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

The first November 19 extension is a separate immutable installation at
`/data/vtenishe/CHEMICAL_DATABASE/snapshots/amps-chemistry-20261119-r2`.
For all examples below, set `CHEM_DB_ROOT` to that snapshot to select r2, or
to the original root to retain r1. Each has its own explicit lock; publishing
r2 does not retarget the original installation or an application's input.

The subsequent physics/source extension is
`/data/vtenishe/CHEMICAL_DATABASE/snapshots/amps-chemistry-20261119-r3`.
Select **that snapshot's own lock** for the new sources and model definitions.
Neither r1 nor r2 is edited, and no application input is retargeted. This README
is shipped inside r3; earlier snapshots retain their original documentation.

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

The first (r1) retained inventory contains 27 source groups, 622 source files and
35 table descriptors. Under the reviewed source policy, 31 source payloads
are importable; 591 source records are metadata-only. Seven table objects
are available, including four CSV representations and three native ASCII
tables. Exact installed counts are in `publication.json`; scientific source
issues and permission state are recorded per dataset, not hidden in totals.

The r2 extension contains 33 source groups, 721 source records and 80 table
descriptors: 119 permitted source payloads, 602 metadata-only records and 52
available tables. The 99 added source records include 88 admitted payloads
and 11 metadata-only records. All 45 new CSV tables are available. These counts
include native sources, rights/metadata files and reproducible projections;
they are not counts of independent experiments or physically qualified laws.

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

### 2.1 New r2 tables and their physical scope

Use these complete, stable table IDs with either the Python API or CLI:

```bash
CHEM_DB_ROOT=/data/vtenishe/CHEMICAL_DATABASE/snapshots/amps-chemistry-20261119-r2
"$CHEM_DB_ROOT/bin/chemdb" --lock "$CHEM_DB_ROOT/database-lock.json" \
  table nasa-selected-species --csv
"$CHEM_DB_ROOT/bin/chemdb" --lock "$CHEM_DB_ROOT/database-lock.json" \
  table nasa-selected-coefficients --csv
"$CHEM_DB_ROOT/bin/chemdb" --lock "$CHEM_DB_ROOT/database-lock.json" \
  describe collision-D107513-fit-parameters
"$CHEM_DB_ROOT/bin/chemdb" --lock "$CHEM_DB_ROOT/database-lock.json" \
  table surface.tuwien.hjn8d-hyh84.Fig1a_Experiment.table --csv
"$CHEM_DB_ROOT/bin/chemdb" --lock "$CHEM_DB_ROOT/database-lock.json" \
  describe exomol-20261119:h2o-dtu-423-15K-photo:table
"$CHEM_DB_ROOT/bin/chemdb" --lock "$CHEM_DB_ROOT/database-lock.json" \
  describe exomol-20261119:oh-mythos-states-bz2:table
```

Native originals remain available by `source-path SOURCE_ID`. For example,
`exoatom-sodium:exoatom-Na__NIST.states` and
`exoatom-sodium:exoatom-Na_p__NIST.states` provide the 423 Na I and 186 Na II
source-relative state rows. Both source definitions have `Ionisation=null` and
the same listed mass; neither supplies a physical ion rest mass or the energy
offset between charge spectra. Their underlying NIST lineage is retained under
the separately published ExoAtom/ExoMol product, without clearing direct ASD
sources held in the original catalogue.

Important interpretation requirements are recorded in the dataset `files`
metadata and accompanying `NATIVE_*_SOURCES.md` notes:

* NASA's native September 2021 thermo.inp is pinned to the acquired CEA Git
  commit. The 19 selected gas-species rows and 56 interval rows retain source
  FORTRAN `D` exponents. Elemental molecular weights and thermochemical
  Hf(298.15 K)/H(298.15 K)-H(0) reference quantities are not state-resolved
  isotope masses or directly usable zero-temperature ionization energies.
  Nonincreasing condensed-phase native domains remain documented in the
  original full file; no interval is silently repaired.
* IAEA's 21 acquired Na/He source fit pages yield 172 printed coefficients.
  The source formulas, thresholds, fit domains and reported uncertainties
  remain native metadata. No source Python/Fortran expression is executed or
  repaired by the database. A coarse fit threshold is not a formation energy;
  integral ionization fits do not determine collision-conditioned daughters.
* ExoMol H2O absorption measurements have 7933 rows each at 423.15 and
  573.15 K. Their 679 and 1448 negative samples are preserved. Do not use
  these signed measurements as nonnegative probabilities or invent branching.
  OH MYTHOS has 10606 source states and a 10000-row partition table. Its actual
  maximum listed state energy differs from the definition's declared maximum;
  both remain recorded. Excitation coordinates are not formation energies,
  and no initial-state population or photodissociation law is selected.
* TU Wien supplies four measured sputter-yield files (17 rows) and fourteen
  simulation files (2033 rows). The sample-identifier discrepancy between
  article and methods remains unresolved. Missing experimental temperature,
  confidence convention, retention and outgoing charge distributions remain
  missing. Sputtering cannot substitute for an incident-ion sticking model.

Full acquisition/gap details are in `PHYSICAL_INPUT_STATUS_20261119.md`.
SpuBase, CCC, ADAS and other uncleared sources remain metadata-only. No r2
source-access PASS closes I21–I24 or activates a production reaction family.

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
  --release-id NEW_RELEASE_ID \
  --source-packet /path/to/explicit/source_packet.json
```

`output-root` must not already exist. The source inventory explicitly accounts
for unavailable packets. Original bytes remain read-only; partial publication
artifacts survive errors for diagnosis and are never mistaken for a complete
release. Rebuilding from another source checkout requires its corresponding
retained source packets; runtime access to an existing installation does not.
The optional `--source-packet` can be repeated; omitting it retains the original
inventory. Exact admitted manifests are preserved under the release's
`provenance/source_packets/`; native/derived objects and processing recipes
are separately hash-bound. Acquisition is a separate explicit online command,
never an automatic runtime fallback. See `SOURCE_PACKET_CONTRACT.md`.

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

## 9. Shared physical inputs and explicit model definitions (r3)

These interfaces are application-independent. The same immutable catalogue,
native data views and C++ helpers can be used by the Moon, comet, Earth or any
other application. There is no PIC, lunar-frame or lunar-species dependency.

The r3 additions preserve source roles rather than turning every downloaded
column into a reaction law:

| Input | Exposed representation | Scope and limits |
|---|---|---|
| CHIANTI O I | Native `.elvlc`/`.wgfa`, five low levels and nine links as CSV | Observed/theoretical energies separate; vacuum angstrom wavelengths; no supplied A-value uncertainty. Radiative data are one coherent source set. |
| Jones2020 Apollo 14163/10084 TPD | Native 14-file deposit and 13 reversible line/field-token tables | Keep measurements, background subtraction, fitted components, inverted site distributions and lunar-cycle predictions distinct. Malformed original headers, blanks, negative values and distinct temperature axes survive. |
| Chang2021 water photofragments | Original supplementary ZIP/XLSX and sparse long-cell tables | Preserve cell coordinates, types, cached values and unevaluated formulae. Narrowband spectra, fitted curves and state populations do not specify a full solar-spectrum joint daughter generator. Unknown units remain unknown. |
| JWST prompt OH | Author spectrum and separately labelled observed/model CSV views | The observed spectrum has no supplied uncertainty column. These astrophysical data are not lunar validation observations. |
| Poppe2022 ARTEMIS pickup products | Author aggregate observation grid, separate modeled rate table | Preserve source units and bins; do not infer Na-only composition, unreported error bars or the meaning of zero/coverage gaps. This is not a time-dependent plasma driver. |
| NASA-CR114270 proton–helium scattering/capture | Original NTRS19700032941 PDF and exact official citation/rights record | Tables8–12 have distinct angular/state scopes, units and errors. No OCR/hand-transcribed cross-section CSV is admitted; genuine blanks and ambiguous print remain source evidence. |

Each table's `describe` metadata records its exact native parents, projection
recipe, units/roles and known limitations. Native objects and original rights
notices remain accessible via `source_path`. Read the shipped
`NATIVE_CHIANTI_SOURCES.md`, `NATIVE_APOLLO_SOURCES.md`,
`NATIVE_WORKBOOK_SOURCES.md` and `RIGHTS.md` before interpreting a table.
`NATIVE_COLLISION_BENCHMARK_SOURCES.md` records the NASA table-specific scope
and narrow record-specific public-use determination; it is not a blanket
NASA public-domain or CEA licence assertion.

The r3 catalogue contains 41 source groups, 782 source records and 105 table
descriptors: 180 admitted source payloads, 602 metadata-only records and 77
available tables. Relative to r2, all 61 added source objects and 25 added CSV
tables are available, including the 12-row authored model-definition table.
These are access counts, not counts of independently measured physical laws.

The model registry includes mathematical mass balance, coherent radiative
branching, the Opal lower-energy electron marginal, conservative two-body final
states, an explicitly unresolved reservoir closure and first-order thermal
release. It also describes published water-site and H-reflection scenarios.
Scenario descriptions are not implemented providers, measured universal
parameters or automatic source selections. In particular, a descriptor with
`implementation=descriptor-only` has no numerical kernel.

```bash
CHEM_DB_ROOT=/data/vtenishe/CHEMICAL_DATABASE/snapshots/amps-chemistry-20261119-r3
"$CHEM_DB_ROOT/bin/chemdb" --lock "$CHEM_DB_ROOT/database-lock.json" list --kind models
"$CHEM_DB_ROOT/bin/chemdb" --lock "$CHEM_DB_ROOT/database-lock.json" describe opal-energy-marginal-v1
"$CHEM_DB_ROOT/bin/chemdb" --lock "$CHEM_DB_ROOT/database-lock.json" table chianti-oi-levels-five --csv
"$CHEM_DB_ROOT/bin/chemdb" --lock "$CHEM_DB_ROOT/database-lock.json" table amps-reduced-model-definitions-v1 --csv
PYTHONPATH="$CHEM_DB_ROOT/python" python3 -B - "$CHEM_DB_ROOT/database-lock.json" <<'PY'
import sys
from amps_chemical_database import Database
db = Database.from_lock(sys.argv[1])
print(db.model("polanyi-wigner-fixed-sites-v1"))
print(db.table("chianti-oi-radiative-lowest-five").rows)
PY
```

`models()` returns recursively immutable definitions selected by the lock;
`model(id)`/`describe(id)` return one. Old catalogues return an empty tuple for
`models()` and `Unavailable` for a new model ID. CLI `list --kind models` gives
the same definitions. C++ applications read the identical registry as the
ordinary `amps-reduced-model-definitions-v1` table through `Catalogue::TableRows`;
no new export codec or runtime ABI is required.

For calculations include `sdk/src/models/chemical_database/ReducedModels.h`
and use `AMPS::ChemicalDatabase::ReducedModels`. Its six header-only helpers
take **all** physical parameters explicitly, return the existing typed `Status`,
and preserve output objects when rejecting input. Read the installed
`sdk/src/models/chemical_database/REDUCED_MODELS.md` for signatures, units,
frame conventions, domain/resource behavior and caller responsibilities.
The helpers have no built-in species constants, surface parameters, source
ranking or reaction network. Kernel access does not bypass source permissions.

A standalone C++ caller example is shipped as `sdk/examples/reduced_models.cpp`:

```bash
c++ -std=c++17 -O2 -Wall -Wextra -pedantic \
  -I "$CHEM_DB_ROOT/sdk/src/models/chemical_database" \
  "$CHEM_DB_ROOT/sdk/examples/reduced_models.cpp" -o /tmp/amps-reduced-models
/tmp/amps-reduced-models "$NU_S" "$BINDING_J" "$KB_J_PER_K" \
  "$SURFACE_T_K" "$POPULATION" "$DURATION_S"
```

Set those six variables from the caller's explicit admitted parameter set;
none is assigned a species/material default by this example. Output reports
the rate, remaining/released inventory and representability flags. A refusal
returns exit2 with its reason. No chemistry-source download, parameter fit,
temperature-to-brightness conversion or application activation occurs.

The Opal marginal omits recoil and joint angles; no Na width is invented.
Two-body conservation requires the caller to supply a supported direction,
state-dependent Q, and total momentum including incident photons/plasma.
The parent-velocity approximation transfers exact signed momentum, energy and
charge differences to an unresolved bath; it predicts no electron spectrum or
heating partition. Fixed-site release assumes constant temperature and
parameters over the integration interval and integrates only that finite
inventory, without readsorption or competing channels. Source-state mappings,
observable-specific model-error bounds, numerical gates and application
admission remain required before physical activation.

To reproduce an r3-style publication, supply the original r2 packet selection
**and** each reviewed new packet, then add `--include-model-definitions` to the
builder command in section8. This flag is opt-in: omitting it preserves the
historical default selection. The release archives the builder/model registry,
exact admitted packet JSON, native parents, derived outputs and recipe hashes.
Locally authored definitions are labelled software under AMPS's Apache2.0
licence; they have no invented HTTP acquisition record. They do not copy or
clear held ASD, ADAS, PHIDRATES or other third-party data.

All source/integrity/API and analytical-helper PASS outcomes are software or
local kernel verification. They establish neither full daughter distributions
nor linked application verification or observational validation.

## 10. Additional physical inputs: explicit r4 selection

The r4 release adds molecular radiative state/transition arrays, seven IAEA
excited-Na composite records, Apollo15221 TPD/derived site distributions, six
soil-water workbooks and three OH(A)+H2 quenching source-data archives. Read
`PHYSICAL_INPUT_STATUS_20261119_R4.md` for the exact source roles, rights,
metadata discrepancies and remaining physical admission requirements. Original
root/r2/r3 locks and payloads are unchanged; an existing application keeps its
previous selection until its own explicit input is changed.

```bash
CHEM_DB_ROOT=/data/vtenishe/CHEMICAL_DATABASE/snapshots/amps-chemistry-20261119-r4
"$CHEM_DB_ROOT/bin/chemdb" --lock "$CHEM_DB_ROOT/database-lock.json" identity
"$CHEM_DB_ROOT/bin/chemdb" --lock "$CHEM_DB_ROOT/database-lock.json" verify
"$CHEM_DB_ROOT/bin/chemdb" --lock "$CHEM_DB_ROOT/database-lock.json" query ohplus
"$CHEM_DB_ROOT/bin/chemdb" --lock "$CHEM_DB_ROOT/database-lock.json" describe surface-apollo15221-water-tpd-csv
```

This example opens the actual new radiative matrix and measured TPD view
through the installed production API. Native energy/wavenumber units are
cm^-1; A is s^-1. Rows are exact strings, not converted floating-point values.

```bash
PYTHONPATH="$CHEM_DB_ROOT/python" python3 -B - "$CHEM_DB_ROOT/database-lock.json" <<'PY'
import sys
from amps_chemical_database import Database
db = Database.from_lock(sys.argv[1])
states = db.table("molecular-20261119-r4:ohplus-mollist:states-table")
transitions = db.table("molecular-20261119-r4:ohplus-mollist:transitions-table")
print("OH+ actual states/transitions:", len(states.rows), len(transitions.rows))
print(transitions.columns)
print(transitions.rows[:3])
tpd = db.table("surface-apollo15221-water-tpd-csv")
print(tpd.columns)
print(tpd.rows[:3])
PY
```

C++ callers use the same table IDs with the existing `Catalogue::TableRows`;
no Moon namespace, MPI, SPICE or new runtime ABI is required. Opaque native
originals remain available through `source-path` only when their source policy
permits storage/evaluation. Metadata-only records keep their IDs/URLs/hashes
and refuse payload access with `PermissionRequired`.

Each OH/H2 state view names only the four unambiguous mandatory fields; later
native quantum tokens remain positional where their metadata disagree. State
counts that disagree with source declarations are retained, without creating
missing rows. Internal spectroscopic energies and A values do not establish
absolute formation energies or complete photon/collision daughter laws.
IAEA capture-plus-ionization composites remain separate from partial CEX.
Source temperature, pressure, specimen, collision energy and calculated versus
measured labels accompany surface/quenching data; no lunar default is inferred.

The new release corrects an earlier MYTHOS JSON `max_energy` interpretation:
the matching native definition labels it a transition-completeness energy,
not the maximum state energy. `amend_exomol_metadata.py` reproduces this
metadata-only correction from pinned old/new evidence. Previous snapshots and
every retained native/CSV value remain unchanged. New source-specific readers
and exact packet recipes/dependencies are archived for offline reproduction.
Neither publication nor these queries activates physics or changes M3 status.

## 11. Targeted remaining-input sources: r5 selection

The new independent snapshot is
`/data/vtenishe/CHEMICAL_DATABASE/snapshots/amps-chemistry-20261120-r5`.
It adds the NIST H/O evaluation's exact table-cell views, Laffon's native
sticking deposit and primary/model-source metadata from the targeted follow-up.
The previous root/r2/r3/r4 selections remain unchanged. Read
`PHYSICAL_INPUT_STATUS_20261120.md` and `NATIVE_NIST_HXOY_SOURCES.md` before
interpreting the new entries: access to a source is not a state/Q profile,
physical reaction model or formal M3 qualification.

```bash
CHEM_DB_ROOT=/data/vtenishe/CHEMICAL_DATABASE/snapshots/amps-chemistry-20261120-r5
"$CHEM_DB_ROOT/bin/chemdb" --lock "$CHEM_DB_ROOT/database-lock.json" verify
"$CHEM_DB_ROOT/bin/chemdb" --lock "$CHEM_DB_ROOT/database-lock.json" describe nist-hxoy-2016-t13-cells
"$CHEM_DB_ROOT/bin/chemdb" --lock "$CHEM_DB_ROOT/database-lock.json" source-path laffon-sticking-13274840-v1-native
PYTHONPATH="$CHEM_DB_ROOT/python" python3 -B - "$CHEM_DB_ROOT/database-lock.json" <<'PY'
import sys
from amps_chemical_database import Database
db = Database.from_lock(sys.argv[1])
table = db.table("nist-hxoy-2016-t13-cells")
print(table.columns)
print(table.rows[:3])
raw = db.source_path("laffon-sticking-13274840-v1-native").read_bytes()
print("Original sticking table bytes:", len(raw))
PY
```

Table3 concerns OH formation enthalpies at298.15K and mixed historical methods;
Table13 mixes ion-energy quantities. Both preserve native cells, spans, units,
footnotes and source conventions without selecting a physical value or assuming
a zero-K isotope state. The sticking file retains its undeclared eight-bit
encoding and original layout. Its XPS-derived values use a gold S=1 calibration,
20K surface and300K gas; it supplies no energy-accommodation law or broad
lunar-temperature extrapolation. The normal C++ Catalogue reads the same two
CSV views; opaque originals use the ordinary source-access path.

## 12. Helium charge-exchange archive and sodium source cross-check

Select `/data/vtenishe/CHEMICAL_DATABASE/snapshots/amps-chemistry-20261124-r6`
with its own lock. The Swaczyna Zenodo deposit (DOI 10.5281/zenodo.4555716,
CC BY 4.0) includes the unmodified producer ZIP and two exact numerical
views for its calculated resonant `He0He1ct` branch. Both views support
the ordinary C++ Catalogue and Python CSV APIs. For a single-energy lookup,
use the smaller index to select the native ZIP member.

```bash
CHEM_DB_ROOT=/data/vtenishe/CHEMICAL_DATABASE/snapshots/amps-chemistry-20261124-r6
LOCK="$CHEM_DB_ROOT/database-lock.json"
"$CHEM_DB_ROOT/bin/chemdb" --lock "$LOCK" table swaczyna-4555716:he-ct-index --csv
"$CHEM_DB_ROOT/bin/chemdb" --lock "$LOCK" describe iaea-d119990-20261124:live-html
PYTHONPATH="$CHEM_DB_ROOT/python" python3 -B - "$LOCK" <<'PY'
from decimal import Decimal
import sys
from zipfile import ZipFile

from amps_chemical_database import Database

db = Database.from_lock(sys.argv[1])
index = db.table("swaczyna-4555716:he-ct-index")
row, = (row for row in index.dictionaries()
        if Decimal(row["energy_cm_eV_native"]) == Decimal("10"))
with ZipFile(db.source_path("swaczyna-4555716:native-zip")) as archive:
    member = archive.read(row["native_member"]).decode("ascii")
print(row)
print(member.splitlines()[18])
PY
```

The second table, `swaczyna-4555716:he-ct-angles`, joins to the index
via `log10_e_cm_x1000` and contains 507,648 native rows. Energies are
center-of-mass collision energies in eV. Angles are in radians and follow
the neutral direction in the [source paper's Appendix A](https://pmc.ncbi.nlm.nih.gov/articles/PMC8861972/).
The tabulated differential already includes `2π sin(θ)`. Integrate by
weighting that column with the source quadrature weights. The native
area unit is `2.8002852056e-17 cm²`. The publisher adds identifiers and
CSV delimiters. It preserves every numeric spelling and angle order.

The new live IAEA D119990 HTML confirms the existing 12-knot sodium
partial table. Groningen and Courbin evidence, the NORAD He photon
tables, Curtin electron-impact ionization and other inputs can be found
through `describe chemistry-archive-20261124` as source records. The
general shared snapshot holds only metadata for sources with restricted
or uncleared managed access. Describe the record before selecting a source
for a particular process. See the source/rights specifics in
`NATIVE_M3_ARCHIVE_SOURCES.md`.

## 13. Published laboratory observables and neutral ionization fits (r7)

Select this snapshot to access four small source tables through the
same CLI, Python API or C++ Catalogue. Example queries follow.

```bash
CHEM_DB_ROOT=/data/vtenishe/CHEMICAL_DATABASE/snapshots/amps-chemistry-20261124-r7
LOCK="$CHEM_DB_ROOT/database-lock.json"
"$CHEM_DB_ROOT/bin/chemdb" --lock "$LOCK" table cui2025-he-capture:table1 --csv
"$CHEM_DB_ROOT/bin/chemdb" --lock "$LOCK" table kim1992-neutral-ionization:table1 --csv
"$CHEM_DB_ROOT/bin/chemdb" --lock "$LOCK" table dibeler1966-water:observed-onsets --csv
"$CHEM_DB_ROOT/bin/chemdb" --lock "$LOCK" table dibeler1966-water:relative-ion-signals --csv
```

Cui et al. report separate measured and calculated helium capture
percentages. The missing excited-state measurement at 0.5 keV/u is blank.
These projectile energies are in keV/u. The r6 Swaczyna table uses
center-of-mass collision energy in eV. Agreement checks require an
explicit convention match. Kim's coefficients are read with his Eq. (2).
The historical water data are measured onsets and parent-relative ion
signals, with no detector mass-response correction. Source locators,
limitations, and source-specific rights appear in
`NATIVE_M3_LITERATURE_SOURCES.md` and the catalogue.

## 14. Subshell photoionization and helium-ion excitation (r8)

The r8 release provides the original 103-element Kalha et al. workbook,
the NASA Sutton report and the Bhatia paper. Use the existing read-only
API from any application. To retrieve a different element, use the original
workbook and its explanatory sheet.

```bash
CHEM_DB_ROOT=/data/vtenishe/CHEMICAL_DATABASE/snapshots/amps-chemistry-20261124-r8
LOCK="$CHEM_DB_ROOT/database-lock.json"
"$CHEM_DB_ROOT/bin/chemdb" --lock "$LOCK" table kalha-yeh-lindau-v2:he-ne-na-ar --csv
"$CHEM_DB_ROOT/bin/chemdb" --lock "$LOCK" table sutton1972-he:excitation-ratios --csv
"$CHEM_DB_ROOT/bin/chemdb" --lock "$LOCK" table bhatia2022-he:hybrid-table2 --csv
"$CHEM_DB_ROOT/bin/chemdb" --lock "$LOCK" table bhatia2022-he:excited-table8 --csv
PYTHONPATH="$CHEM_DB_ROOT/python" python3 -B - "$LOCK" <<'PY'
import sys
from amps_chemical_database import Database

db = Database.from_lock(sys.argv[1])
for row in db.table("kalha-yeh-lindau-v2:he-ne-na-ar").dictionaries():
    if row["neutral_atom"] == "Na" and row["subshell"] == "3s":
        print(row["photon_energy_eV_native"], row["photoionization_sigma_Mb_native"])
print(db.source_path("kalha-yeh-lindau-v2:workbook"))
PY
```

The four-element view preserves the exact spreadsheet strings and native
blanks. Its cross sections are theoretical, in Mb, at tabulated photon
energies in eV. Sutton reports gas-phase orbital cross-section ratios at
200 eV incident electron energy with instrument-width sensitivity and the
author's estimated errors. Bhatia's helium table uses outgoing electron
momentum in the source convention. His two Table 8 ratios describe
excited He⁺ relative to ground-state He⁺ and preserve the printed decade
jumps that need independent review. Check the source locator, exact units,
licence and domain details in `NATIVE_M3_PARTIAL_SOURCES.md`. The release
also records the embargoed Science Data Bank helium angular dataset and
the water-cation study as metadata.

## 15. Charge-exchange scattering and the D119990 review (r9)

The `isotropic-com-two-body-scattering-v1` entry describes uniform
solid-angle sampling in the incoming pair's center-of-mass frame.
It describes an experimental direction law for any two-body channel.
The caller supplies the rate and product states. The conservative
two-body model supplies the momentum magnitude and transformation with
the declared reaction energy.

```bash
CHEM_DB_ROOT=/data/vtenishe/CHEMICAL_DATABASE/snapshots/amps-chemistry-20261124-r9
LOCK="$CHEM_DB_ROOT/database-lock.json"
"$CHEM_DB_ROOT/bin/chemdb" --lock "$LOCK" describe isotropic-com-two-body-scattering-v1
"$CHEM_DB_ROOT/bin/chemdb" --lock "$LOCK" describe proton-sodium-ground-capture-evidence-v1
"$CHEM_DB_ROOT/bin/chemdb" --lock "$LOCK" describe collisiondb_119990_embedded_arrays.csv
PYTHONPATH="$CHEM_DB_ROOT/python" python3 -B - "$LOCK" <<'PY'
import sys
from amps_chemical_database import Database

db = Database.from_lock(sys.argv[1])
print(db.model("isotropic-com-two-body-scattering-v1")["equations"])
integral = db.describe("collisiondb_119990_embedded_arrays.csv")
print(db.model(integral["source_scope_review_id"])["limitations"])
PY
```

C++ applications can read the entries through `Catalogue::TableRows`
on `amps-reduced-model-definitions-v1` and the D119990 table note
through `Catalogue::GetTable`. The D119990 CSV retains the 12 native
integral cross-section rows. The table now links to the unresolved
ground-product rate/state comparison and the 4 keV/u angle gap.
The entries describe evidence and an experimental model. Both have
`physical_activation=false`. [Source review](NATIVE_CEX_SCATTERING_SOURCES.md)

## 16. OH surface loss and helium excitation-ionization (r10)

```bash
CHEM_DB_ROOT=/data/vtenishe/CHEMICAL_DATABASE/snapshots/amps-chemistry-20261124-r10
LOCK="$CHEM_DB_ROOT/database-lock.json"
"$CHEM_DB_ROOT/bin/chemdb" --lock "$LOCK" table wang2023-quartz-oh:table2 --csv
"$CHEM_DB_ROOT/bin/chemdb" --lock "$LOCK" table wang2023-quartz-oh:table3 --csv
"$CHEM_DB_ROOT/bin/chemdb" --lock "$LOCK" table sutton-kay1972-he:table3 --csv
PYTHONPATH="$CHEM_DB_ROOT/python" python3 -B - "$LOCK" <<'PY'
from amps_chemical_database import Database
import sys

db = Database.from_lock(sys.argv[1])
for row in db.table("wang2023-quartz-oh:table2").dictionaries():
    print(row["Y_mm_native"], row["plasma_feed"], row["loss_probability_beta_native"])
print(db.source_path("sutton-kay1972-he:preprint"))
PY
```

Wang and Bruggeman Table 2 gives atmospheric-pressure fused-silica OH
loss estimates with the source's two distinct flux scales. Native
missing values stay blank. Table 3 compares two spatial resolutions
at the same dry-helium location. Use the source conditions when
interpreting any loss probability. The Sutton and Kay file gives
200 eV gas-phase He⁺ excitation cross-section ratios relative to 4s.
It preserves the asymmetric and symmetric uncertainty estimates from
Table III. The earlier single-author NASA report remains separately
queryable. The ordinary shared C++ Catalogue uses the same CSV IDs.
[Source details and limits](NATIVE_M3_LAB_BENCHMARKS_SOURCES.md)

## 17. Reaction-specific laboratory and solar inputs in r11

Select the snapshot by its lock. The table interface returns the original
column tokens. The dataset description exposes experimental conditions,
theoretical assumptions, source locators, and individual measurements
from papers retained as metadata.

```bash
CHEM_DB_ROOT=/data/vtenishe/CHEMICAL_DATABASE/snapshots/amps-chemistry-20261124-r11a
LOCK="$CHEM_DB_ROOT/database-lock.json"
"$CHEM_DB_ROOT/bin/chemdb" --lock "$LOCK" table moore-mapaye2025:oh-water-fits --csv
"$CHEM_DB_ROOT/bin/chemdb" --lock "$LOCK" table laurell2025:argon-peak-centers --csv
"$CHEM_DB_ROOT/bin/chemdb" --lock "$LOCK" table iupac2010-oh-ice:uptake --csv
PYTHONPATH="$CHEM_DB_ROOT/python" python3 -B - "$LOCK" <<'PY'
from amps_chemical_database import Database
import sys

db = Database.from_lock(sys.argv[1])
for dataset in ("sodium-ground-capture-benchmark-20261124", "he-electron-ion-emission-20261124"):
    for source in db.describe(dataset)["files"]:
        print(source["id"], dict(source["quantity_metadata"]))
PY
```

Use the two sodium calculations as integral theory benchmarks with their
individual incident-energy conventions. The NIFS figure identification
disqualifies D119990 as a ground-product experimental rate. The helium
measurements are photon-emission cross sections for electron impact. The
flow-tube OH results describe irreversible uptake on ice. The solar fits
provide frequencies for individual products. The argon fitted centers
measure a fine-structure energy separation under narrow-band illumination.
The reports document parent states, apparatus, uncertainties and
applicability. [Sodium](NATIVE_M3_SODIUM_BENCHMARK_SOURCES.md)
[Helium](NATIVE_M3_HELIUM_EMISSION_SOURCES.md)
[OH on ice](NATIVE_M3_OH_ICE_SOURCES.md)
[Partial photorates and argon](NATIVE_M3_PHOTON_BENCHMARK_SOURCES.md)

For the monochromatic water experiment, query
`db.describe("kronebusch-water-photofragmentation")["files"]`. Kronebusch
and Berkowitz report fragment ions relative to total ionization at 21.2 eV.
The collection efficiency drops for ions above 1 eV. Use the printed
proceedings page when checking the weak O⁺ signal.
[Water source and replay](NATIVE_M3_WATER_BRANCH_SOURCES.md)
[ASMS primary report, printed page 114](https://www.asms.org/docs/default-source/proceedings-archive/1975_asms_23rd_conference.pdf#page=117)

## 18. Reproduce published exosphere model choices in r12

Poppe, Halekas and Harada's original 1 au rate file covers 25 species. The
CSV exposes five atomic rows with the native species tokens and all nine
frequencies per row. The article identifies the archive's `Ar` row as ⁴⁰Ar.
Each rate is per neutral per second. Use the minimum, median and maximum as
the authors' three environmental conditions. [Deposited rates](https://doi.org/10.5281/zenodo.6994648)

```bash
db=/data/vtenishe/CHEMICAL_DATABASE/snapshots/amps-chemistry-20261124-r12
"$db/bin/chemdb" --lock "$db/database-lock.json" table poppe2022:atomic-1au-rates --csv
PYTHONPATH="$db/python" python3 -B - "$db/database-lock.json" <<'PY'
from amps_chemical_database import Database
import sys

db = Database.from_lock(sys.argv[1])
print(db.source_path("poppe2022:ionization-rates"))
for dataset in ("water-photofragment-model-precedents-20261124",
                "lunar-oh-surface-model-precedents"):
    for source in db.describe(dataset)["files"]:
        print(source["id"], dict(source["quantity_metadata"]))
PY
```

The water records distinguish neutral fragment speeds, assumed charged
branch fractions, and the products returned by a published sampler. The OH
records distinguish incident retention, residence time and thermal launch.
Each record identifies its native source and the limits on its application.
The cross-model report links the individual studies, code and regression uses.
[Comparison and sources](NATIVE_M3_MODEL_PRACTICE_SOURCES.md)

## 19. Lunar neutral losses, water coverage and noble-gas product assumptions

The r13 snapshot adds three separately identifiable model inputs. The first
query exports the exact quiet-Sun frequencies from Peschel et al. Table 1.
The other two return the parameters for molecular water on lunar grains and
the single-channel He, Ne and Ar model used by EXOSPHID.

```bash
db=/data/vtenishe/CHEMICAL_DATABASE/snapshots/amps-chemistry-20261124-r13
"$db/bin/chemdb" --lock "$db/database-lock.json" table peschel2026:quiet-sun-neutral-loss --csv
PYTHONPATH="$db/python" python3 -B - "$db/database-lock.json" <<'PY'
from amps_chemical_database import Database
import sys

catalogue = Database.from_lock(sys.argv[1])
for dataset in ("schorghofer-lunar-water-models-v123", "exosphid-noble-gas-photoionization"):
    for source in catalogue.describe(dataset)["files"]:
        print(source["id"], dict(source["quantity_metadata"]))
PY
```

Compare the Ne frequency under the author's stated solar conditions and
the Ar threshold under the source's stated product convention. The water
desorption flux depends on coverage in molecules/m² and surface temperature
in K. Its species is adsorbed molecular H₂O. The source notes document
the differences from PHIDRATES, NIST and incident-OH sticking.
[Noble-gas loss rates](NATIVE_M3_LUNAR_NOBLE_GAS_SOURCES.md)
[Water surface model](NATIVE_M3_LUNAR_WATER_COVERAGE_SOURCES.md)
[EXOSPHID atomic channels](NATIVE_M3_EXOSPHID_ATOMS_SOURCES.md)

## Sodium pickup and binary charge exchange (r14)

Select a reference case by its model ID. A source-bound model appears only in
a release containing its reviewed source records. Python exposes the selectors
as `db.model(id)["source_bindings"]`. The CLI and C++ table reader expose the
same selectors through `amps-model-source-bindings-v1`, keyed by `model_id`.

| Model ID | Physical input and intended use |
| --- | --- |
| `poppe2022-na-pickup-median-v1` | Total Na⁺ pickup production at 1 au for the authors' median conditions. The Na row supplies photoionization, solar-wind proton charge exchange and electron impact. The published frequencies to three significant digits are 7.59 × 10⁻⁶, 8.40 × 10⁻⁷ and 5.00 × 10⁻⁷ s⁻¹ per neutral. The ion starts at zero velocity in the lunar frame. Multiply each frequency by neutral density where that process operates. The authors neglect gravity and use steady perpendicular wind and magnetic field for their pickup trajectories. [Poppe et al., Table 1 and §2.2](https://doi.org/10.1029/2022JE007422) · [Deposited data, CC BY 4.0](https://doi.org/10.5281/zenodo.6994648) |
| `hplus-na3s-ground-capture-4kev-u-isotropic-v1` | Monoenergetic H⁺ + Na(3s) → H(1s) + Na⁺(2p⁶). Shingal's A34 calculation supplies 0.22 Å² = 2.2 × 10⁻²¹ m² at 4.0 keV/u in the target-rest frame. AMPS supplies the angular assumption. Draw uniform cos θ on [−1,1] and uniform azimuth on [0,2π) in the incoming pair center of mass, and give the products opposing momenta. Use the actual reactant momentum, product masses and state-dependent reaction energy Q for recoil. The published point has no stated numerical uncertainty. Off-energy use is undefined. [Shingal et al., Table 4](https://doi.org/10.1088/0022-3700/19/3/014) |

The first case includes all contributing hydrogen final states. Adding the
second case to that total would count ground-state capture twice. To study the
effect of the birth law on Na⁺ energies and surface returns, hold the production
rate, source locations and fields fixed across the sensitivity runs. Such a
comparison measures dependence on the approximation.

Hartle et al. justify the zero birth velocity by comparing typical exospheric
and background plasma speeds. Shen et al. use the same assumption to trace
lunar pickup trajectories in their ARTEMIS analysis. These studies support the
bulk transport approximation. The Na frequencies come from the Poppe deposit.
[Hartle et al.](https://doi.org/10.1029/2011JA016859) · [Shen et al.](https://doi.org/10.3847/1538-4357/ad40a4)

From any directory, inspect both cases and obtain the values from the pinned sources:

```bash
db=/data/vtenishe/CHEMICAL_DATABASE/snapshots/amps-chemistry-20261124-r14
"$db/bin/chemdb" --lock "$db/database-lock.json" table amps-model-source-bindings-v1 --csv
PYTHONPATH="$db/python" python3 -B - "$db/database-lock.json" <<'PY'
from decimal import Decimal
import sys
from amps_chemical_database import Database

db = Database.from_lock(sys.argv[1])
pickup = db.model("poppe2022-na-pickup-median-v1")["source_bindings"]
na, = (row for row in db.table(pickup["table_id"]).dictionaries()
       if all(row[key] == value for key, value in pickup["row_match"].items()))
rates = {process: Decimal(na[column]) for process, column in pickup["columns"].items()}
print(rates, pickup["unit"], tuple(pickup["birth_velocity_m_s"]), pickup["velocity_frame"])
print("charge-exchange fraction when all three processes apply",
      rates["solar_wind_proton_charge_exchange"] / sum(rates.values()))

collision = db.model("hplus-na3s-ground-capture-4kev-u-isotropic-v1")["source_bindings"]
point, = (source for source in db.describe(collision["dataset_id"])["files"]
          if source["id"] == collision["source_id"])
for key in collision["metadata_path"]:
    point = point[key]
print(point["projectile_energy"], point["energy_unit"], point["energy_frame"],
      point[collision["area_field"]], collision["area_unit"])
print(collision["angular_model_id"], collision["recoil_model_id"])
PY
```

The original Poppe numerical file and its projection retain their attribution.
The Shingal selector resolves to the published fact preserved in bibliographic
metadata. Both sources retain their original access policies. Database tests
verify the selections independently of an application enabling transport.
