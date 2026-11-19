# Shared chemical database implementation contract

2026-11-18. The user authorized development and installation at
`/data/vtenishe/CHEMICAL_DATABASE`, overriding the proposal's suggested local
archive path. Existing originals, scientific gaps and source permissions remain
binding. No commit, push, reference replacement, provider contact or automatic
physical chemistry activation is authorized by this work.

## Inventory, interfaces and scope

The source inventory and independent audits are recorded in
`doc/CHEMISTRY_DATABASE_PROPOSAL.md` and its dated evidence. Source revision is
`74974c2b2fc1b4cf09a73c320064a22370c0ca94` plus the recorded working tree.
Source/data roadmap hashes remain respectively
`712e1152bf50af6d9ee1dad8f8a585809ca8b92cf42ac6fc7a009f6c0253ee17`
and `d10eca1c539f9cf82940236d9e1aa36af665d43e82ba475925950d582fbf49c4`.
Their shared qualification rules/native metadata control science; approved
source-tree development/resource exceptions control work order.

Production services currently reside in Moon Photochemistry, ElectronImpact
and ChargeExchange; shared ReactionNetwork/CollisionDatabase containers and
Hartley2's startup executor have separate semantics. Extract the pure rate
services unchanged into `src/models/chemical_database`, retaining Moon aliases
and updating actual standalone build dependencies. Add a C++17 immutable
catalogue reader using existing AMPS block grammar, strict CSV and SHA-256.
No C++ SQL, JSON, PIC, MPI or SPICE dependency is introduced in that reader.

Python preparation/API/CLI lives in `utility/chemical_database`. It inventories
every retained provider, imports only permission-compatible payload, preserves
exact originals and source locators, and publishes immutable JSON/CSV releases
plus acyclic detached runtime exports. It never fetches data implicitly. An
explicit installation lock pins one release/profile/export, without searching
for latest data. Unknown/permission-restricted sources remain metadata-only.

The installed archive includes its Python package/CLI, portable C++ SDK,
examples, frozen data/indexes and detailed README. A prepared installation is
verified in a writable staging directory before the authorized external copy.
Existing unknown destination contents must not be removed or silently replaced.
The source-of-truth implementation remains in AMPS for normal review.

## Runtime export contract

One `#begin chemical_database: runtime` block, no includes/children, schema
`amps-chemical-runtime-v1`, codec `csv-block-v1`. It pins release/profile JSON
and source/table CSV indexes by root-relative path and SHA-256. Canonical
release does not include the digest of an export that cites it. Source/table
permission and availability are independent; absent payload has no invented
value. Native observed digests remain in catalogue provenance even where an
object cannot be imported. Production only reads; no automatic model selection,
clipping, source fusion, ground state or missing distribution is permitted.

## Test specification

| Field | Contract |
|---|---|
| ID/name | Shared chemical-database component and existing Moon supporting tests; no permanent new U/I IDs |
| Purpose | Actual loader/query/evaluator correctness, provenance, permission refusal, compatibility and relocation |
| Production capability/entry points | Shared C++ catalogue/rate services and Python inventory/builder/API/CLI; actual existing Moon wrappers |
| Configuration | Explicit root/release/profile/export pins, declared capabilities and original numerical model selection |
| Inputs/reference | Real native packets plus clearly analytical parser/invariant fixtures; independent source rows, Decimal/math formulas, exact integer identities and reviewed regression where applicable |
| Units/frame/time | Retain native typed units and spectrum/beam/state/material conventions; data access does not choose simulation geometry/epoch |
| Random controls | No RNG in loader or scalar services; existing application seeds/order remain unchanged |
| Procedure | Inventory/rights evidence, immutable build, native/independent comparisons, corruption/missingness refusals, relocation, normal application preservation and external installation verification |
| Metrics/acceptance | Exact byte/row/ID/hash preservation where declared; unchanged existing independent arithmetic/statistical/numerical criteria; no relaxed physical thresholds |
| PASS/failure modes | Only evaluated criterion met; wrong hashes/schema/path/rights/units/domain or silent missing-data zero fails explicitly; outputs unchanged on refusal |
| Required data/hashes | Frozen input manifests plus SHA-256 in each release/run; source omissions and unsupported material/state/kinetic data remain explicit |
| Predecessors | Source integrity then scalar/binding/wiring; formal M1/M2/M3 predecessors unchanged |
| Artifacts | result/manifest, source/build/data identities, exact commands, stdout/stderr, detailed comparisons, install lock and resource histories |
| Status | PASS evaluated/met; FAIL evaluated/missed; ERROR parser/build/runtime/corruption/infrastructure; SKIPPED unavailable prerequisite; source/data QA is not physical validation |

Tests must include actual deployed Python and C++ examples, relocation without
source-tree dependency, corrupted/escaped paths and changed bytes, permission/
capability refusal, immutable lifetime and compiled namespace compatibility.
Keep old neutral defaults and reviewed four-rank `make -j test_Moon` zero-byte
diff. Shared-core changes require preserving historical PASS membership; scope
and any deferred long gate must remain explicit. Native packages use the normal
remove-build/Config.pl/make path, eight MPI ranks/OMP1 plus explicit controls,
200GB process-tree ceiling and one original 3600-second deadline including
setup/tests/restoration. A source-only or analytical PASS never closes M3.

Every new implementation path needs detailed comments, explicit units and
failure behavior; application/test/module READMEs and STATUS record actual
scope/results. The remaining missing state/energy/kinetic/material inputs and
provider permissions cannot be supplied by implementing this database.
