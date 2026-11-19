# Detached chemical runtime export: csv-block-v1

This is the implemented runtime projection, not a new canonical JSON schema.
The source catalogue/release/profile remain authoritative and are reproduced
by the maintained Python importer/exporter. C++ authenticates their bytes but
has no JSON parser. The detached export cites their final hashes; the final
release JSON must not in turn cite the hash of this same detached export.
An application pins release, profile **and export** identities independently.

One manifest contains exactly one block, without includes or child blocks:

```text
#begin chemical_database: runtime
schema = amps-chemical-runtime-v1
codec = csv-block-v1
release_id = <reviewed-release-id>
release_manifest_path = <root-relative-final-release-json>
release_manifest_sha256 = <actual-lowercase-sha256>
profile_id = <reviewed-profile-id>
profile_path = <root-relative-profile-json>
profile_sha256 = <actual-lowercase-sha256>
source_index_path = <root-relative-source-index-csv>
source_index_sha256 = <actual-lowercase-sha256>
table_index_path = <root-relative-table-index-csv>
table_index_sha256 = <actual-lowercase-sha256>
profile_capabilities = catalogue;source-access;table-access
topology_fingerprint = none
#end chemical_database
```

Placeholders deliberately make this example non-runnable; they are not science
evidence. All fourteen keys are mandatory. Unknown keys/blocks/schema/codec
are refused. Hashes have exactly64 lowercase hexadecimal digits. Identifiers
use ASCII letters/digits plus `_-. :` (without the displayed separating blank).
A topology fingerprint is `none` or a SHA-256; the application separately checks
actual compiled topology before activating a physical model. The generic
`AMPS::Input` grammar handles assignments, quoting, comments and source locations.
This owner forbids includes because every child file must be explicitly pinned.

All paths are forward-slash **archive-root-relative**. Absolute paths, `..`
components, backslashes and canonical/symlink escapes are refused. Relative
references are resolved against the explicit archive root, never the export
file directory or an application working directory. The input owner resolves
the archive root itself from its normal `.in` declaring-file convention.

## Source index

Exact ordered header:

```text
source_id,path,sha256,media_type,availability,storage_permission,evaluation_permission,capabilities,source_uri,notes
```

`source_id` is unique. `path/sha256` identify an optional retained payload;
`source_uri` retains its distinct publisher locator. An available object needs
both path and digest. For an unavailable/metadata-only/quarantined source, both
may be empty; one without the other is corruption. An observed historical digest
for an absent or permission-held object stays in canonical metadata rather than
pretending a retained runtime payload exists. `media_type` is explicit/nonempty.

Availability: `available`, `unavailable`, `metadata-only`, `quarantined`.
Permissions: `allowed`, `pending`, `prohibited`, `unknown`, `not-applicable`.
Only **allowed** authorizes payload access; not-applicable is representable but
is not an authorization. Metadata lookup remains available for all records.

## Table index

Exact ordered header:

```text
table_id,path,sha256,format,schema_id,source_ids,availability,evaluation_permission,capabilities,columns,row_count,notes
```

`table_id` is unique; `schema_id` names its versioned scientific/table role.
`source_ids` contains at least one existing source ID. Lists use semicolons,
without trimming or duplicate entries. Empty capability lists are allowed;
`none` is the manifest/list empty sentinel. `format` is exactly `csv`,
`native-ascii` or `opaque`; an unsupported source format can be retained as
an opaque source object without an invented evaluator.

For CSV, `columns` is the exact ordered header joined by semicolons, with
nonempty unique column labels. Spaces, units and punctuation are retained;
semicolon/newline column labels require another codec. `row_count` is an
explicit nonnegative integer and must match parsed rows. Every row must have
the declared number of fields. Quoted commas, doubled quotes and quoted
newlines are preserved; bare CR, NUL, malformed quotes, invalid UTF-8 and
trailing text after a closed quote are rejected. No source value is normalized,
converted, clipped, smoothed or filled by this parser.

For native-ascii/opaque, `columns` is empty and `row_count` is `unknown`.
Their bytes are authenticated; a separately selected source-specific parser
owns their scientific interpretation. `TableRows` refuses these formats.

A table access needs its explicit capability in BOTH profile and table, allowed
table evaluation permission, available table payload, and allowed storage AND
evaluation permissions plus available payloads for EVERY dependent source.
Dependent sources need not redundantly declare `table-access`; their own
capability list applies to direct `SourcePath` calls. A blocked source cannot be
laundered through an otherwise allowed derived table.

## Lifetime, resources and status

Objects/CSV rows are prepared before transport. Available, permitted entries
whose capabilities intersect the profile are authenticated; blocked metadata
is not opened. Immutable storage owns metadata and parsed rows. Copies share
that storage. `TableRows` returns an owning `shared_ptr<const Csv>` that outlives
the catalogue. Const concurrent lookups do not mutate tables or caches.
Legacy solver non-owning views must retain an owning catalogue/rows handle.

`RequireCapability` checks declared **selection scope**, not global scientific
completeness or an execution qualification. Requesting `scalar-rates` never
turns a missing product/state distribution into an available kinetic model.
`CodeName` spells outcomes exactly `Evaluated`, `Unavailable`, `InvalidInput`,
`OutsideDomain`, `CorruptSource`, `UnsupportedSchema`, `PermissionRequired`,
`NumericalFailure`, `ResourceLimit`. These are data-access outcomes, not test
PASS/FAIL or observational-validation labels.

File and total byte ceilings are explicit `BundleRequest::limits`; zero means
no artificial ceiling. Files are not truncated. Unrepresentable sizes, explicit
limit exhaustion and allocation failure refuse publication. Applications map
fatal production requirements to their detailed exit diagnostics; tests retain
ERROR versus genuine unavailable-prerequisite SKIPPED semantics. Opening a bad
bundle preserves the previous catalogue; failed lookups preserve outputs.
