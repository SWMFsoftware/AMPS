# Shared chemical database verification contract

Written before test/helper changes,2026-11-18, source revision
`74974c2b2fc1b4cf09a73c320064a22370c0ca94` plus existing uncommitted work.
The user authorized the shared database and installation at
`/data/vtenishe/CHEMICAL_DATABASE`. Tests stage disposable copies under their
owned output directory; this test work does not deploy or edit that database.

## Implementation inventory and production boundary

Moon currently exposes `species/Photochemistry`, `ElectronImpact` and
`ChargeExchange` and standalone compilers link those CPPs directly. These
public APIs are becoming namespace compatibility aliases for maintained
`src/models/chemical_database` services: Runtime, DataIO, Photochemistry,
ElectronImpact, ChargeExchange and RateServices. Existing numerical formulas,
frozen source identities, units, domains and scientific thresholds must stay
unchanged. No test may compile copied formulas or replace a real installer.

The standalone compile owners requiring shared dependencies are:
`test_m3_runtime_configuration.py`, `test_m3_runtime_rates.py`,
`run_m3_transport_clock.py`, `test_m3_ion_collector.py`,
`test_m3_helium_silica.py`, `test_m3_sodium_partial_capture.py`,
`U16_electron_impact/voronov_verification.py`, `common/m3_rates/test.py`, and
`common/m3_rates/test_photon_selection.py`. A shared source-list helper will
retain both compatibility-wrapper and actual shared-source/header identity.
The existing `production_source_identity` recursively hashes `src/` already;
the new production subtree is included without an added exclusion or shortcut.

The native runtime uses `export.in`: one strict `chemical_database: runtime`
block, `amps-chemical-runtime-v1` and `csv-block-v1`, explicit release/profile
identities and SHA-256-bound relative source/table indexes. Native labels and
the actual Runtime.h public schema control field interpretation. The Python
builder/API and C++ Catalogue must agree on availability, rights, capability,
source ordering and table data. A content hash authenticates bytes; it cannot
convert missing metadata, unknown rights or an unqualified model into science.

## Authority and scientific scope

Source roadmap SHA-256:
`712e1152bf50af6d9ee1dad8f8a585809ca8b92cf42ac6fc7a009f6c0253ee17`.
Canonical-data roadmap SHA-256:
`d10eca1c539f9cf82940236d9e1aa36af665d43e82ba475925950d582fbf49c4`.
Approved source development/resource amendments govern work order; native
data labels/external AGENTS govern real inputs and original formal gates
continue to govern qualification. No observation, uncertainty, cross section,
state mapping, permission or scientific capability is invented by a fixture.

Parser/metadata fixtures are explicitly administrative controls. Altered-byte
copies exercise rejection only and never replace canonical raw/processed
science data. Missing sources/permissions/capabilities remain unavailable;
directory existence or successful import is not readiness or validation.

## Test specification

| Required field | Contract |
|---|---|
| ID/name | Shared chemical-database software controls; descriptive component names only, no new permanent Moon U/I ID. |
| Purpose | Verify production shared C++/Python database access, integrity, honest missingness/rights, immutable ownership and source-equivalent numerical services. |
| Production entry points | `Catalogue`/`OpenRuntimeBundle` and shared rate APIs; actual `amps_chemical_database.Database` builder/API/CLI; existing Moon alias APIs. Exact method signatures come from maintained headers/package. |
| Configuration | Standalone C++17 and Python3, no AMPS application configuration or MPI build. Explicit source roots, lock/release/profile identities and caller-selected capabilities. |
| Inputs | Real source-bound staged release where available; clearly labelled copied metadata/parser rejection controls. Relative manifest/index/table paths and exact expected hashes. |
| Independent references | SHA-256 from Python hashlib; exact source/index ordering and native CSV rows; existing independent Decimal/analytical rate oracles and frozen source specifications; explicit availability/rights/capability rules. |
| Units/frame/time | Preserve native units/schema. Metadata/path tests have no physical frame or epoch; rate tests preserve existing SI/domain conventions and explicit quiet/active endpoints. |
| Deterministic controls | No new random selection; existing statistical controls retain their seeds/counts/criteria. No observational selection or residual tuning. |
| Procedure | Build/stage real database into new evidence root; compile actual shared services and thin native adapter; compare Python/native records; mutate independent copies for hash/path/schema/permission tests; move root and repeat; prove snapshots/lifetimes survive reload attempts. |
| Metrics | Exact status/order/schema/permission/capability/missing-value identities; hash equality; table bytes/rows; numeric source agreement at unchanged existing tolerances; zero partial publication on failed load/query. |
| Acceptance | All evaluated exact invariants match. Existing rate tolerances remain unchanged. Corrupt/path-escaping/self-inconsistent input refused before data use; missing/denied data remains unavailable rather than zero or fabricated READY. Native service output agrees with existing independent reference. |
| Expected PASS | Declared software/integrity/API criterion evaluated and met. Source/storage QA and model validation remain separate. |
| Expected failures | Missing hash/source, changed bytes, traversal/symlink escape, duplicate keys/IDs, unknown schema/codec, malformed rows/types, inconsistent counts, rights/capability denial, lost source order, stale mutable aliases, partial reload. |
| Required data/hashes | Source-bound staged release and its builder/lock/profile/export hashes; canonical data/root remains read-only. Missing real evidence produces SKIPPED for that dependent group; parser controls can still run separately. |
| Predecessors | Software tests do not close scientific M1/M2/M3 gates or authorize absent physical providers. |
| Artifacts | Global result/counts, individual test outcomes, exact runner/log paths for FAIL/ERROR, commands, compiler/executable/source identity, staged release/input copies, actual stdout/stderr and mutation fixtures. |
| Status semantics | PASS evaluated criterion; FAIL evaluated disagreement; ERROR build/runtime/parser/corrupt test infrastructure; SKIPPED real unavailable prerequisite. Expected rejection is PASS only for its declared specific guard. |

The single global runner must preserve every child outcome, list executed,
PASS/FAIL/ERROR/SKIPPED counts, and identify failing runner/log paths. It must
not classify scientific FAIL as ERROR or hide a child error behind SKIPPED.
Prior meaningful evidence is never overwritten. One bounded development
package includes setup, compilation, tests, summaries and cleanup; longer
qualification is not silently shortened. No full application build is run here.
