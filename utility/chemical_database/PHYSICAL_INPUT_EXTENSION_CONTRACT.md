# Additional physical-input release contract — 2026-11-19

This package implements the user's request to search online databases first,
then relevant primary publications, and integrate newly obtained physical
inputs into the shared chemical database. It changes source access and
reproducible source projections, not a Moon reaction rate, event model or
qualification threshold. Source HEAD at the start is
`26f1a65c58a8a70f1c6301b1a6c4d283a4558c94`. The pre-existing local Makefile
change is outside this package.

## Inventory and production entry points

Existing releases r1, r2 and r3 supply authenticated native objects, permission
metadata, exact CSV projections and twelve inert reduced-model definitions.
`source_inventory.inventory` and `source_packets.load_packets` are the
publication admission path. `amps_chemical_database.build.build_database`
constructs a new immutable release. The installed Python `Database`, CLI and
C++ `Catalogue` are the production query paths used by source-access checks.
The SDK's mathematical helpers and particle-chemistry providers remain separate
from source admission; adding a source does not enable them.

New molecular, collision and surface lanes work in independent source packets.
They preserve exact original bytes, retrieval status, source-native units,
states, domains, uncertainties, measured/model roles and pinned rights notices.
Only a reviewed existing source policy admits payload storage. A public URL,
article access, or successful HTTP request does not grant numerical reuse.
Failed requests stay in acquisition reports and cannot become source records.
No original under canonical validation-data `raw/` is changed.

## Source scope and independent oracles

* Molecular spectroscopy: native ExoMol definition/state/transition files
  supply source-labelled internal energies and Einstein A coefficients. Count
  mismatches remain explicit; optional quantum labels that disagree with native
  tokens remain positional. Spectroscopic excitation is not an absolute
  molecular formation energy or a complete PHIDRATES product-state mapping.
* Collision data: new IAEA native records and printed fit coefficients retain
  projectile/target states and the distinction between capture, ionization and
  a composite total. Missing displayed equations/domains are not reconstructed
  from neighboring datasets. Printed precision is not full native precision.
* Surface/relaxation: native author deposits retain specimen, preparation,
  incident conditions and the distinction between measured data and inverted
  or calculated quantities. An adsorption isotherm or TPD inversion is not an
  unconditional lunar sticking/retention law. Fixed-energy quantum scattering
  data are not Maxwellian thermal rates.
* Atomic publication search: native documents and their explicit permissions
  control admission. Missing isotope/state/uncertainty mappings and restrictive
  or unclear payload rights remain holds; no indexed search snippet or manually
  transcribed numerical value is used as managed science data.

Native source bytes and published schemas are the independent reference for
lexical projections. Mathematical identities and explicitly labelled
administrative fixtures test parsers, refusal behavior and the shared API;
they are never observations. No comparison with AMPS residuals selects data.

## Verification contract

| Field | Declared contract |
|---|---|
| ID/name | Supporting shared chemical-source release QA; no new lunar U/I ID |
| Purpose | Authenticate new source data and expose them through shared production APIs without altering previous releases |
| Configuration | Explicit packet list; explicit release ID; twelve existing inert model definitions retained; physical activation false |
| Inputs | Pinned successful native downloads, source policies, packet recipes and all previous r3 packets |
| Reference/oracle | Original byte hashes and source tokens, native metadata, independent cross-language API queries |
| Units/frame/time | Preserve native units, source energy zero, isotope/state labels, retrieval UTC and source conditions; no inferred common frame or epoch |
| Stochastic controls | None; deterministic source/lexical/software checks; no Monte Carlo or observational scoring |
| Procedure | Authenticate inputs; run source-parser controls and real-source tests; build new release; compare Python/C++ table cells and permissions; test relocated installed SDK; install into a new snapshot; recheck old identities |
| Metrics | Exact source hashes, tokens, row counts, declared/actual discrepancies, queried table counts, software outcomes and installed identities |
| Acceptance | All evaluated software/source checks pass; unchanged original lock and protected old members; missing physical metadata remains explicit |
| Expected PASS | Source/integrity/API operations meet their stated contract; no implied physical validation or milestone closure |
| Failure modes | Changed/truncated source, wrong units/schema, unresolved state ID, duplicate ID, unauthorized payload, changed recipe, unavailable original, build/query/install failure |
| Required data/hashes | Exact selected packet manifests and every member/recipe SHA-256 recorded in the publication plan and release provenance |
| Predecessors | Source rights/schema review and pinned originals; formal M3 predecessors are not replaced by these supporting checks |
| Artifacts | Acquisition receipts, packet/recipe originals, parser QA, software logs/results, source/release/build identities, install receipt and preservation report |
| Status | PASS/FAIL/ERROR/SKIPPED have the AGENTS.md meanings; unresolved physical inputs stay NOT QUALIFIED even when source QA passes |

## Preservation and corrections

Publish a new snapshot; do not replace r1/r2/r3, retarget the default lock,
alter source observations, activate a channel or change an acceptance gate.
The work package pins 494 pre-existing installed members before acquisition.
All are rechecked after installation. Corrections to an earlier interpretation
are documented in the new release with the new native evidence; published old
metadata and original bytes remain unchanged.

The source roadmap SHA-256 is
`712e1152bf50af6d9ee1dad8f8a585809ca8b92cf42ac6fc7a009f6c0253ee17`;
the canonical data roadmap SHA-256 is
`d10eca1c539f9cf82940236d9e1aa36af665d43e82ba475925950d582fbf49c4`.
Their differences remain explicit: approved source early-development/resource
amendments govern development order; canonical scientific/data and predecessor
contracts govern qualification. This package changes neither roadmap.
