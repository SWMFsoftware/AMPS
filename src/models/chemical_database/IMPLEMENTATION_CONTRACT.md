# Shared chemical runtime implementation contract — 2026-11-18

Inventory before change: Moon Photochemistry, ElectronImpact and ChargeExchange
already provide dependency-free PHIDRATES/Voronov/total-capture production
services. Photochemistry owns the reviewed SHA-256 and strict CSV helpers;
source schemas/hashes, rate equations, native domains, named interpolation and
transactional refusal semantics are frozen. Existing shared ReactionNetwork and
collision containers are solver/execution views; no executor is migrated here.

This package moves the SAME services to AMPS::ChemicalDatabase, preserving
public class/function names and frozen historical schema strings. The parent
will replace Moon definitions with namespace aliases and normal-build wiring.
There is no PIC, MPI, SPICE, application header or mandatory third-party runtime
dependency. Existing generic AMPS block grammar is the sole runtime manifest
parser. No JSON parser is introduced: JSON catalogue/release bytes are hashed,
while their detached immutable block/CSV export supplies typed metadata.

Entry points: ReactionData file/byte SHA and strict CSV; original named rate
service APIs; Catalogue/OpenRuntimeBundle; capability/permission-gated source
and table lookup; factories that prepare the original services only from
permitted, source-authenticated selected table roles. No missing data returns
physical zero. No loader or service activates particle chemistry.

The runtime export contract is csv-block-v1, schema amps-chemical-runtime-v1,
with one chemical_database: runtime block. Source/table indexes and all selected
files are hash-bound; immutable catalogue storage owns every exposed view.
Paths resolve from a supplied root and reject absolute/traversal/symlink escape.
Unknown keys/schema, duplicates, bad hashes/columns/counts, dangling source IDs,
unpinned release/profile/export identities and unsupported capability requests
are refused. Permission- or science-unavailable entries remain visible as
metadata without exposing their numerical bytes. Their presence is not PASS.

Model formulas, coefficient units, optical/solar endpoint conventions and
published validity domains remain unchanged. Legacy source-authenticated APIs
remain available as compatibility services; the Catalogue factories add the
selected release's capability/rights guard before calling them. Metadata-only,
quarantined and incomplete datasets are never silently repaired or promoted.

Independent oracles: published SHA-256 vectors, native frozen source hashes,
Python hashlib/strict CSV/schema rejection fixtures, existing independent
PHIDRATES/Voronov/CEX reference tests and full parent-preservation campaigns.
This owned package only compiles/runs a standalone source smoke; it does not
configure AMPS, touch canonical sources or claim linked/observational PASS.
The parent owns exhaustive shared tests, old adapter aliases and core builds.

Failure contract: bool legacy calls preserve previous output/registry; typed
runtime outcomes distinguish evaluated, unavailable, invalid input, outside
domain, corrupt source, unsupported schema, permission required, numerical
failure and resource limit. Publication occurs only after validation. Explicit
caller resource ceilings never truncate tables/channels; zero ceiling means no
artificial cap and allocation failures return a detailed ResourceLimit status.
No network access, mutable latest release, RNG, runtime writes or hot reload is
performed. Source/provider selection and actual geometry/bath state remain
application-owned. Full activation/qualification order is unchanged.
