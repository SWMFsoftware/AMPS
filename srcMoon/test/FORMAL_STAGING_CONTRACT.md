# Formal I19 staging and initialization identity

This contract precedes the 2026-11-16 correction. It extends existing I19
without changing its numerical design, source laws, criteria or test ID.

The audit found that `m2_protocol.write_case` called the bounded I18
`prepare_na_ar_input` helper even for the formal cluster tier. That helper
enables Ar and replaces the template volume AMR. The native entry is
`MoonIntegrationTests::RunM2NumericalCase`; production initialization is
`Moon::SurfaceTemperature::InitializeProductionState`. The cluster adapter
previously authenticated a template rather than proving that its physics
survived staging. A separate audit identified unapplied thermal-state proof;
the native initialization adapter and reviewed evidence validator address
that distinct issue.

For cluster cases, preserve all supplied AMR, source, surface, thermal and
clock settings. Require explicitly enabled Ar, an evolving environment and
zero initial surface inventory; missing prerequisites are SKIPPED rather
than silently enabled. Only predeclared numerical/seed controls and output
ownership or suppression may differ. Resolve input paths relative to the
file containing the assignment, preserving case and the native no-escape
value semantics. The native grammar treats `!` as a comment delimiter even
inside quotes; staging does not replace that rule with shell parsing.
Development retains its approved bounded fixture transformation.

The independent oracle is a strict canonical configuration comparison before
launch: template and staged cluster physics must agree. Unknown assignments
remain part of identity, so a future field cannot silently disappear under a
prefix exclusion. Explicit numerical keys, seed and output-only settings are
excluded by a named list. Input paths contribute resolved meaning; external
file contents remain authenticated by the existing native data manifests.
The exact completion-only approval/diagnostic transport keys
`m2_formal_contract_file`, `m2_illumination_queries_file` and
`m2_polar_queries_file` are also excluded from the initial physical identity;
their files are independently hash-bound by the corresponding formal
contracts. This exception does not extend to arbitrary new file settings.
The initialization mode and spin-up settings remain physical identity fields.
The native declared-uniform default is represented explicitly in identity.

Parser/arithmetic tests use temporary input text only, never observations.
They cover distinct formal/development staging, literal path case, relocation,
source/AMR/thermal drift, unknown fields, duplicate/malformed input and missing
prerequisites. Expected PASS means the staging/identity invariant was actually
evaluated. Contradictory syntax/identity is ERROR; unavailable required formal
configuration is SKIPPED. No fixture result qualifies M1, the long I19 design,
thermal initialization, or lunar observations.

Required retained artifacts are the unchanged reviewed template, actual case
input, canonical physics identity/hash, native initial full-column state,
source/build/executable/data hashes and ordinary runner logs/receipts. The
existing same-campaign predecessor and reviewed thermal-state gates remain.

The same distinction applies to `prepare_m2_normalization.prepare`: its
default development staging keeps the existing bounded physical fixture and
full-duration deterministic source quadrature. Explicit `tier=cluster` /
`--tier cluster` preserves the reviewed physical template. READY references
must retain their staging tier, actual input and canonical physics identity,
native manifest/receipt hashes and provenance of the ordinary-initialization
weights used by the controlled initial stock. A development reference cannot
qualify formal normalization merely because its original-template hash
matches. The source integral is model-input evidence, never an observation or
transport/convergence PASS. If selected production initialization is periodic,
that initialization still runs before this source-only native adapter.

The source-only normalization additionally accepts `record_initial_state=True`
or `--record-initial-state`, a diagnostic request that changes no thermal
method, epoch, coefficient or numerical acceptance. Whenever the staged input
requests export (automatically in cluster tier), an actual native full-grid
CSV and its receipt counts are required. The same streaming state-grid and
actual furnished-SPICE/configured-DEM validators are used; missing/corrupt
requested evidence is ERROR. A model-input READY result may retain an
explicitly **unqualified** `candidate_initial_state`, with actual file/hash,
shape, method, TDB epoch and physical input/data identities. It never creates
an approved/READY initial-state qualification document or review.

This adapter launches only the retained native binary and therefore checks
its preserved generated-header identity without requiring the live generated
tree to be the same configuration. This permits a focused post-restore
execution without compiling a probe against unrelated headers. The existing
source/executable/header hash checks remain mandatory.
