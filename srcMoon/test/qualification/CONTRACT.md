# M3 qualification infrastructure and stream contract

Written before this implementation,2026-11-18, revision
`74974c2b2fc1b4cf09a73c320064a22370c0ca94` plus recorded working changes.
The original `M3_QUALIFICATION_PROTOCOL.json` and per-ID criteria stay immutable.
An explicit approved addendum records the user's newly approved numerical and
statistical choices; previous runs are never reinterpreted with the addendum.

## Inventory and production boundaries

`integration_m3.py` currently owns I21–I24 refusal/artifact handling, while
`m3_protocol.py` validates the original preparation-only specification. Existing
rate oracles and transactional/controller tests are supporting evidence only.
The production `MoonReactionController::OpenUniformProvider` requires strict
open uniforms and currently relies on a caller-supplied provider. No Philox or
Random123 implementation is available. SEP's `KeyedUniform01` can return zero
and is not silently substituted. Local/source/kernel and linked conservation,
numerical, statistical and observational gates remain distinct.

Source roadmap SHA256
`712e1152bf50af6d9ee1dad8f8a585809ca8b92cf42ac6fc7a009f6c0253ee17`;
data roadmap SHA256
`d10eca1c539f9cf82940236d9e1aa36af665d43e82ba475925950d582fbf49c4`.
Approved source development-order/resource amendments govern work order;
unchanged formal criteria and native data labels govern qualification.

## Approved design and unresolved science

The user approved10% of the unchanged ensemble error gates as numerical bias
budgets: I21/I22 relative0.002 and I23/I24 relative0.003. The primary homogeneous
I22 estimator is the endpoint-survival estimator `-log(S/N)/T`; its original2%
criterion applies to total parent removal. Branch frequencies and partial-rate
products are diagnostics plus separately frozen branch gates. Linked branches
use simultaneous95% exact Clopper–Pearson intervals over the entire declared
family, with10% relative interval half-width for every positive source branch.
U15's99% and U17's100000-event minimum stay unchanged. Complete positive/zero
channel lists are preserved; no channel is dropped/renormalized. Large cohorts
belong to a separately requested cluster tier rather than a shortened test.

The independent unit is an original parent, with equal initial represented
number within each binomial cohort. Siblings, splits and repeated samples are
not additional independent trials. Unequal product reference weights/fractional
corrections are allowed only when physical packets preserve that initial
cohort; genuinely unequal initial trial weights require another reviewed
statistical design. Homogeneous cohorts must not silently censor surface loss
or domain escape as photo/electron/CEX survival. Time-dependent transport needs
its own independent cumulative-exposure oracle and cannot be fit as one
constant-rate cohort. Product/state/source profiles, native executors and
same-campaign predecessor evidence remain mandatory and may still be absent.

## Random stream representation

New PIC-independent `QualificationRandom` exposes caller-owned semantic keys:
campaign seed, immutable original-parent ID, genealogy path, purpose, and draw
index. Genealogy appends explicit event/daughter ordinals for a NEW history;
MPI/ownership clones preserve the same key/counter. Rank, process ID, memory
address, allocator slot, timestep and number of rejected trajectory previews
are absent. A queried index is stateless; only the accepted consumer advances
its counter. No rollback/capacity retry may draw a second physical event.

HMAC-SHA256 uses the maintained production `DataIO::BytesSha256`, standard
64-byte inner/outer pads, an8-byte big-endian campaign seed, a versioned ASCII
domain separator, length-prefixed big-endian genealogy/purpose and64-bit
original-parent/draw indices. The first52 digest bits select the exact binary64
midpoint `(k+1/2)*2^-52`, strictly in(0,1). This is a deterministic finite-grid
PRF construction, not a mathematical claim of infinite independent real
uniforms. Its smallest value is2^-53; exponential-tail discretization and the
maximum planned draw family must be documented before physical qualification.
Independent Python `hmac/hashlib` vectors, exact midpoint endpoints, stream-key
distinctions, reordered/sharded reads, invalid transactional outputs and
predeclared distribution checks qualify the software construction only.

## Test specification

| Field | Contract |
|---|---|
| ID/name | Descriptive M3 qualification infrastructure; extends existing owners/support groups, no new permanent U/I ID. |
| Purpose | Correct approved planning/scoring, semantic RNG, fail-closed execution eligibility and retained evidence. |
| Production entry points | `QualificationRandom` using actual shared DataIO; actual controller/provider consumers remain parent-owned. Python independent analytical/statistical evaluators score native records. |
| Configuration | C++17/Python3; explicit frozen case/source/profile/protocol identities, seed, genealogy, full channel family, initial cohort, physical horizon and numerical levels. |
| Inputs | Actual native receipts when available; labelled analytical/binomial/arithmetic and infrastructure fixtures only for local controls. No manufactured observations or substituted source values. |
| Independent oracle | Python stdlib HMAC/SHA256; analytical exponential survival and source branching; exact binomial/Clopper–Pearson equations; high-precision/polynomial controls and conservation identities. |
| Units/frame/time | Rates s^-1, relative physical seconds; source/native event frame explicitly bound by the actual execution receipt. RNG fields dimensionless semantic integers/bytes. |
| Deterministic/stochastic controls | Approved source seeds731993,104729,130363; immutable original-parent trials; exact frozen family membership; deterministic counter keys; no post-output estimator/cohort selection. |
| Procedure | Verify addendum/legacy/source identities, freeze exact family and sample plan, validate actual execution/provenance/records, independently score and preserve CSV/PNG/EPS. Missing prerequisites refuse launch/qualification. |
| Metrics | Exact key/vector/ownership and count identities; endpoint rate/lifetime errors; conditional branch confidence and precision; numerical biases relative to explicit independent references; actual ledgers remain separately required. |
| Acceptance | Unchanged I21–I24 and local gates plus approved numerical/statistical addendum. No tolerance increases after failures. Invalid inputs leave pure-layer outputs unchanged. |
| Expected PASS | Declared supporting software criterion evaluated/met, or full formal criteria only after actual native/provider/predecessor evidence exists. |
| Failure modes | Stream/key reuse, false independent trials, dropped channels, nonfinite/corrupt records, accidental censoring, wrong estimator/horizon, unjustified finite-grid tail, missing budgets/providers/predecessors, contradictory source/status/artifacts. |
| Data/hashes | Original protocol/per-ID contracts and both roadmap hashes, exact consumed source/profile/build/executable/record hashes; no raw mutation. |
| Predecessors | Formal same-campaign I13–I20 and applicable U gates retained; early supporting development does not satisfy them. |
| Artifacts | Contract/addendum, source/build/data hashes, command/stdout/stderr, seed/stream keys, native records, result/manifest, comparisonCSV/PNG/EPS and resource histories where execution occurs. |
| Status | PASS evaluated/met; FAIL evaluated/missed; ERROR parser/build/runtime/corruption/infrastructure; SKIPPED genuine missing dependency/provider/approval/predecessor. No fake native receipt for not-run. |

One development package includes preparation, builds, tests, plots and cleanup
under3600s and200GB aggregate memory. No full AMPS configuration/build is run
by this subtask. Actual native campaigns require coordinated source freeze.

## Native evidence seam extension (before implementation)

There is currently no admitted source-complete physical cohort executable.
The receipt boundary therefore binds genuine optional artifacts without
guessing a native CLI or promoting an arbitrary CSV. A caller must supply an
explicit expected SHA256 for its producer-evidence manifest. The manifest
binds the exact plan/addendum/source profile, named producer, execution token,
common-time1/8MPI/OMP1 and zero native completion code, plus actual manifest,
parent/event records, numerical/budget evidence, executable, build manifest
and input snapshot. Every supplied existing file is hash-checked even while
another prerequisite is missing; corrupt evidence is ERROR. Missing actual
files/provider are SKIPPED, and metadata-only/all-green administrative fixtures
cannot install a physical producer. The normal run manifest must independently
match the same token/plan/case/profile. Hash binding is provenance, not scientific
validation or authority to invent a state/energy/ledger tolerance.
