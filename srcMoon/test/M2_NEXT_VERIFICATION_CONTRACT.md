# M2 research and next verification contract — 2026-11-12

The user requested D02 literature/archive research, a concrete I19/I20 proposal,
and continued verification/testing. This contract precedes changes to tests.
It does not approve new numerical tolerances, observation selections, or a
formal M2 predecessor bypass. No commit/push is authorized.

The follow-up request to prepare detailed suggestions for numerical-protocol
review authorizes proposal/documentation changes only. The v2 proposal must
identify existing versus proposed gates, expose the actual Ar timestep copy
and hard-coded fixture seed, define first-release/censored residence and
re-adsorption separately, and state the true two-profile simulation count.
Review uses existing source/receipts and primary statistical documentation;
no new model run, native-data processing, acceptance change or long campaign
is implied. Numerical examples are planning arithmetic, never observations.

## Inventory and scope

Source starts at `06eafde4802f0344694612c8f3c21f2b28e1a09f`. The user's
checkout has a separate configured application and generated Makefile changes;
all native verification therefore uses an independent checkout under `/tmp/`.
The active checkout's build and local configuration are preserved.

The production terrain service currently traces a point-Sun ray. The generic
formal runner explicitly skips I15/I19/I20 when their complete capabilities
are unavailable. Existing U08-U13 exercise production geometry, thermal,
surface and retention kernels. The supporting M2 driver runs the real I14,
I18 and finite-stock I19 fixture, records separate supporting outcomes, and
preserves formal I14-I20 as I13-blocked. No new permanent U/I ID is introduced.

The new read-only D02 metadata audit extends the existing M2 package audit.
It records the frozen raw products, native quantity, projection/frame/scaling,
and unresolved temporal operator. It cannot qualify the reference by guessing
calculation epochs from a file name, label time span, or another publication.
It is invoked by the existing global M2 development mode, not an alternative
validation runner. No model/observation residuals are computed here.

## Scientific/test contract

| Field | Contract |
|---|---|
| IDs/name | Existing D02 data-readiness support for I15/I20 and existing U01-U26/I01-I20 verification; no new permanent test ID. |
| Purpose | Resolve what the frozen illumination products actually measure, expose remaining metadata gaps, and verify the current linked/local production paths. |
| Production capability/entry | `Moon::SurfaceEnvironment` and existing Moon/PIC surface/thermal/mover/source entry points remain unchanged; package audit uses `audit_m2_packages.inspect` and the new read-only D02 inspector. |
| Configuration | Pinned source; existing Config.pl builds; clean generated trees; explicit production TESTMODE=off, separate canonical nightly ON; eight MPI ranks/OpenMP=1 for supporting Moon runs; existing dedicated Na/Ar chemistry-OFF/individual-weight build. |
| Inputs | Canonical D01/D02/D03 roots, frozen raw hash list, detached PDS labels, archived primary documentation, approved Moon examples and existing seeds. |
| Independent oracle | Native product labels and exact downloaded documentation for metadata; analytical/conservation/stochastic oracles already declared in U/M0/supporting contracts. A metadata-protocol fixture is not an observation. |
| Units | Native DN and scale/offset retained; visibility is dimensionless occurrence frequency; pixel scale m/pixel, latitude/longitude degrees; model SI. |
| Frame/time | Native DE421 mean-Earth/polar-axis body frame; no time-system inference. Label START_TIME/STOP_TIME are recorded without assigning them the unproven role of simulation integration bounds. |
| Deterministic controls | Frozen files/hashes and fixed parser checks. Existing linked test seeds/layouts are unchanged. |
| Procedure | Read native metadata and primary sources; freeze proposal before any new formal run; exercise metadata parser/corruption controls; run existing local/M0/supporting campaigns in an isolated checkout; preserve all attempts and source/config/build identities. |
| Metrics | Hash/byte contracts, native observable, documentation claims, missing operator fields; unchanged local/link/budget metrics; separately counted PASS/FAIL/ERROR/SKIPPED. |
| Acceptance | Metadata matches frozen native products; unresolved authoritative temporal metadata stays SKIPPED/PARTIAL. Missing required input is SKIPPED, altered/corrupt input ERROR. Existing scientific gates are unchanged. Proposed numerical budgets are not acceptance criteria until reviewed. |
| Expected PASS | Declared integrity/protocol or existing physical criteria actually evaluated and met. A subcheck PASS never becomes I15/I20 validation. |
| Failure modes | Corrupt bytes, missing labels, wrong dimensions/sample type/scaling/frame, confusing visibility with mean visible area, interpreting acquisition/product dates as simulation bounds, fabricating a temporal contract, build/runtime failure. |
| Required data | D01 LDEM_4; D02 AVGVISIB/LPSR 65N/65S 240 m 201608 IMG/LBL; D03 selected package and SPICE for existing linked tests. |
| Data hashes | The committed `qualify_lola_package.EXPECTED` list and generated audit/manifests; document snapshots have URL/access time/SHA-256. |
| Predecessors | Existing local/M0 requirements for supporting verification. Formal I14-I20 still require same-campaign I13 and all individual gates. No December residual inspection. |
| Artifacts | Research source snapshots/provenance, D02 audit result, proposal Markdown/JSON, native/local logs/results/manifests/config snapshots, exact commands, source/build/executable/data hashes, resource history and completion report. |
| Status semantics | PASS for evaluated declared scope; FAIL for evaluated scientific miss; ERROR for corruption/infrastructure; SKIPPED for unavailable data/capability. M1/M2 remain open unless their original full gates pass. |

## Roadmap authority

Source roadmap SHA-256:
`2761b447b46539edb21395550cd1b5d9ff6add4839b0db5a7e6d734cf2a4f02a`.
External roadmap SHA-256:
`d10eca1c539f9cf82940236d9e1aa36af665d43e82ba475925950d582fbf49c4`.
They differ as documented in ROADMAP_RECONCILIATION.md: source amendments govern
approved execution/development exceptions, common text governs scientific
gates, and native data labels/external AGENTS govern data interpretation.
Neither roadmap, calibration/hold-out membership, nor acceptance JSON is
changed by this work package.

## Discovered metadata defect and minimal correction

The existing D02 spatial processor hard-codes 18.6 years, one hour and 0.5 deg
into `illumination_operator` and returns READY even though no exact-product
averaging epochs or version-specific sampling contract is established. The
research distinguishes general PGDA claims from the 2011 methods and the
201608 native labels. Preserve all native CSV science rows byte-for-byte;
replace unsupported operator values with explicit unresolved fields, report
spatial processing READY separately from overall PARTIAL, and verify both
parser protocol and old/new CSV identity. This correction strengthens truthful
metadata; it changes no observation, model parameter or numerical threshold.


## User-directed execution-budget follow-up

The later user instruction requires an end-to-end development campaign within
one hour, including builds/plots/cleanup; longer qualification may be deferred
to a larger cluster for a few hundred CPU-core hours. Only this operational
policy is adopted in the source Phase-2 roadmap. V3 numerical criteria, short
workload settings and the concrete 300 core-hour planning target remain
proposals; no runner/physics implementation is changed by this documentation
work. Current per-test timeouts do not enforce the required aggregate hour.

The initial contract/historical results above retain their original roadmap
hash. The new source-roadmap hash is `bde5f001bdcbe072d7a8452d3f2afec93856f50cd1146560da3cd3027a96d2dc`;
external roadmap remains
`d10eca1c539f9cf82940236d9e1aa36af665d43e82ba475925950d582fbf49c4`.
ROADMAP_RECONCILIATION.md records the operational-only difference. No scientific
threshold, permanent ID, M1 500/5000 rule, observation or formal predecessor
is changed. This follow-up uses documentation/arithmetic preservation checks
only, with no new native build, model run or external raw-data processing.
