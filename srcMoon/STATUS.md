# `srcMoon` development and validation status

**Status date:** 2026-11-25 UTC

**I23 electron-bath coupling.** The user-approved controlled 10 eV electron
bath is coupled to the normal reaction controller through a new startup
owner. Events use the production helium reservoir composer under a reviewed
narrow-scope admission. Supporting tests passed: 165 ensemble checks and 223
linked PIC checks. The rebuilt seven-case family passed 722 checks, and the
nightly diff is zero bytes. Formal I21–I24 remain 4 SKIPPED.
[Implementation, evidence and limits](#559-i23-prescribed-electron-bath-coupling--2026-11-25-utc)

**M3 campaign runner.** The pinned five-case sodium family ran with one and
eight MPI ranks. All 334 supporting checks passed. The four formal
chemistry gates remain open. [Campaign implementation and results](#558-frozen-m3-development-family--2026-11-24-utc)

**Current M3 chemistry development.** The 4.0 keV/u sodium capture
provider now produces both heavy fragments and accounts for the incoming
proton. Normal Moon startup accepts explicit neutral-H timestep and
statistical-weight controls for this build. The I21–I24 formal gates are open.
[Source, current executions, and remaining gates](#557-linked-m3-chemistry-development--2026-11-24-utc)

**Latest source update — shared chemistry r10.** The catalogue adds
fused-silica OH surface-loss data and the 200 eV Sutton–Kay He⁺
ratios. See section 5.52 for the experimental conditions and the
remaining M3 gates.

**Previous source update — shared chemistry r9.** The D119990 integral
table now points to a source review of the ground-product rate and
angular coverage. The shared catalogue also describes the isotropic
center-of-mass sensitivity model. See section 5.51.

**Previous source update — shared chemistry r8.** The release adds
theoretical He/Ne/Na/Ar subshell photoionization, two theoretical He
photoionization tables, and measured 200 eV He⁺ excitation ratios.
The [source report](../utility/chemical_database/NATIVE_M3_PARTIAL_SOURCES.md)
records the evidence gaps and table restrictions. See section 5.49
for the release and verification record.

**Latest development — I21 homogeneous native execution, ground-ion profiles,
and I23 fit slope.** The Moon aggregate accepts an explicit pinned I21 cohort
plan and preserves the formal predecessor decision. The scoped 8/16/32-step
Na run at 1 and 8 MPI ranks exercised real event conversion and independent
auditing. Its rate, chronology and ledger checks passed. Its strict timing
refinement requirement failed at the arithmetic floor. Native ⁴He, ²⁰Ne and
⁴⁰Ar profiles and the U16 fit-slope controls passed. See section 5.43 for the
open physics requirements and run evidence. Formal M3 remains unqualified.

**Latest completed package — measurement and community-model literature follow-up:**
The new [evidence review](test/M3_PHYSICS_LITERATURE_20261119.md) resolves several
acquisition questions: Na+ mass can be derived from measured inputs; actual
CHIANTI O I radiative arrays are available; original water-fragment spectra,
Apollo14163/10084 desorption data and reduced ARTEMIS ion observations were
acquired. Published heavy-ion/reservoir, electron-energy-sharing, forward-CEX
and heterogeneous-site models provide concrete reduced options with different
scopes. Source mapping, model selection/error assessment and linked qualification
remain required. New research does not activate a provider or alter installed
r1/r2. A native-source check corrected the earlier truncated Na threshold note;
the actual D04 threshold remains unchanged. See section5.34.

**Earlier completed package — online physical-input acquisition and database integration:**
The new immutable snapshot
`/data/vtenishe/CHEMICAL_DATABASE/snapshots/amps-chemistry-20261119-r2`
adds NASA CEA thermochemistry, ExoAtom Na/Na+ states, ExoMol H2O absorption/OH
rovibronic data, IAEA Na/He ionization fits and TU Wien regolith sputter tables.
It adds **99 source records and 45 available tables**, with native values,
provenance, permissions and unresolved metadata preserved. All **143 API/parser/
deployment checks PASS**; the registered CORE case PASS contains **122**
site-independent checks. Actual installed original-source QA passes **36/36**,
including all **45 new C++/Python table comparisons**. All **84 original r1
installation files** and the default root lock remain byte-identical.
The physical Na ion formation/mass/threshold crosswalk, complete daughter and
cascade inputs and applicable ion/regolith response remain unresolved. No model
was activated; M1/M2/M3 formal status is unchanged. See section5.33 and
[the physical-input matrix](../utility/chemical_database/PHYSICAL_INPUT_STATUS_20261119.md).

**Earlier completed package — upstream integration and corrective verification:**
Master integrates `60a156ddcc3deeaf48e2981b0beb14738245457e` with the preserved local source (eleven incoming commits). The completed corrective CORE registry has **844 PASS, 240 SKIP, 1 FAIL and 1 ERROR**. Exact historical/current-upstream PASS preservation audits pass; the known nonpasses remain visible. Moon M0 is **10/10 PASS**, with a zero-byte reviewed nightly diff, and **1533 targeted native checks PASS**. The U-series has 15 PASS/11 SKIPPED and zero FAIL/ERROR.
All local source/site configuration is preserved; no staging, commit or push.
The original parallel-supervisor attempt remains ERROR with its original clock;
a separate single-process corrective campaign supplies complete authentic CORE
rows. Pristine/Moon 822dec07 and delivered/corrective 60a156dd source bases remain distinct.
M1/M2/M3 qualification and missing physics/predecessors remain open. See
section5.32 and [CORE_INTEGRATION.md](CORE_INTEGRATION.md).

Earlier dated checkpoints retain their original source bases and counts. This
package supersedes their upstream-currency statements.

**Latest completed package — M3 independent implementation and verification:**
The source-bound profiles, bounded moving-parent Na composer, accepted-trajectory
exposure, exact compensated clocks, semantic daughter RNG, normal registration,
conservative native publication and independent cohort auditors are implemented.
The latest normal SPICE build passes M0 **10/10**, with a **zero-byte I04 diff**;
the new environment stage passes **154** native wiring checks. The full U15
supporting collection passes **20/20 groups**, while literal U15 remains SKIPPED
for missing physical prerequisites. Eight selected component suites comprise
**625 local PASS criteria**. The separate pre-stage campaign passes normal
registration20, legacy packet945, variable exposure56 and SI-semantic84 checks.
All844 protected CORE PASS cases remain preserved; unchanged shared-source
identities are recorded, without claiming a later redundant CORE rerun.

Three source-scoped native packages finished in3443.73s,3272.40s and3484.83s,
each within its ORIGINAL3600s deadline, including retained failures and required
restoration/reporting. Peak sampled package RSS remains below15.76GB, versus
200GB allowed. Earlier source bases and all oracle/fixture/resource ERRORs are
preserved explicitly. No root Makefile/site/configuration, reviewed reference,
scientific threshold or immutable validation/installed-database bytes changed.
No Git staging, commit or push was performed.

**M3 remains EXPERIMENTAL / NOT QUALIFIED.** Real Na source admission returns77
at1/8MPI for unavailable Na+ state/mass/threshold evidence. Molecular, electron
and complete CEX daughter physics and physical regolith response remain data
blocked; formal I21–I24 are **4 SKIPPED**,0PASS/FAIL/ERROR. M1/M2 predecessors
still control qualification. Downloaded rates, software checks and ideal
collectors do not supply missing physical evidence. See section5.31 and
`test_output/srcMoon/m3-physical-completion-20261118T223145Z/work_package_report.json`.
The final latest-code report is
`test_output/srcMoon/m3-coefficient-stage-20261119T003919Z/final_followup_report.md`.

**Initial database installation — regression preserved:**
`/data/vtenishe/CHEMICAL_DATABASE` now contains the immutable source catalogue,
Python API/CLI, portable C++17 SDK and detailed runnable examples. The installed
archive indexes 27 source groups/622 source records/35 table descriptors; 31
permitted payload records occupy 30 unique content objects, while 591 entries
remain metadata-only. The actual installed API/data/SDK suite passes **56/56**,
the nine migrated Moon scalar/adapter suites pass **558** checks, and exact
installed README examples pass. Linked I01–I10 pass **10/10**, the actual
reviewed four-rank nightly diff is zero bytes, and an additional eight-rank
production rate-callback probe passes all 32 rows. All 843 historical CORE PASS
cases remain PASS; the new shared case also passes. The reconciled 1086 cases are 844 PASS/240 SKIPPED/1 existing FAIL/
1 existing ERROR. Both isolated builds restore ordinary SPICE Moon successfully.
Native execution completed in 3097.47 s (51.6 min); final reporting finished
within 3360 s (56 min), with sampled aggregate peak RSS 9.79 GB. Details and remaining M3 work are in section 5.30.
No physical chemistry capability or formal milestone is activated by catalogue
integrity. Root configuration/Makefile and original scientific data are preserved.

**Historical planning update — shared chemical database:** the original
proposal and audit remain in section 5.29. The delivered initial implementation
and actual storage root supersede its proposal-only status. Broader physical
profile/application migration remains explicitly incomplete; see
[the decision register](../doc/CHEMISTRY_DATABASE_DECISIONS.json).

**Latest completed work — approved P01–P09 components:** reduced recoil/energy,
strict/collector ion-boundary, conditional He/silica and partial Na CEX
components are implemented. Their supporting controls pass (264,62,15,162
respectively); fifteen O I source controls also pass. The corrected full
U-series has **15 PASS /11 SKIPPED /0 FAIL/ERROR**. Actual linked I01–I10
all pass, with the reviewed nightly reference unchanged and a zero-byte diff.
Native packet/collector, force-tracer, controller, common-clock and full normal
startup checks pass (945/89/73/80/32). Setup, preserved earlier attempts,
builds, tests/plots and ordinary-Moon restoration took **2837.98s (47.3min)**;
peak measured aggregate process-tree RSS was **15,745,789,952B**. The original
bootstrap ERRORs and stale U02 source-guard FAIL remain visible; corrected
same-source evidence is explicitly reconciled, not retroactively relabelled.
No model/build/test process remains active. Full details are in section 5.28. Detailed
contract: [M3_APPROVED_MODELS_IMPLEMENTATION_CONTRACT.md](test/M3_APPROVED_MODELS_IMPLEMENTATION_CONTRACT.md).
Missing physical inputs and formal predecessors remain binding; no milestone
closure or physical activation follows from approval. Source preparation
retains native NIST/IAEA provenance, source uncertainties and physical-state
gaps; formal I21–I24 remain **4 SKIPPED /NOT QUALIFIED**.

**Earlier decision proposal:** the concrete daughter-physics/ion-surface
[P01–P09 recommendations](test/M3_DAUGHTER_AND_ION_SURFACE_DECISIONS.md)
and their JSON decision register retain the approved choices. They distinguish
source-supported inputs, explicitly experimental recoil/phase-space/forward-CEX
closures, a verification-only ideal collector and a separate He laboratory
charge-survival benchmark. Section5.27 retains the original proposal scope;
the user's subsequent approval is recorded in the decision register.

**Latest research update:** online chemical databases and the reachable Hartley2
H2O/OH implementation were audited. Available rates, new numerical sources,
source restrictions and the proposed reuse boundary are recorded in section5.26
and [the research report](test/M3_CHEMISTRY_DATABASE_RESEARCH.md). This adds no
formal physics qualification or changed source selection.

**Earlier completed package — independent M3 software/model preparation, 2026-11-18:**
Implemented declared exposure mathematics, its optical-clock/instantaneous-channel
bridge, opt-in common transport chronology/source wiring, concrete physical-model
and qualification specifications, and global-runner coverage. Actual local,
linked-hook, conservative-daughter and full normal startup controls passed.
The original thin-runner ERROR and its corrected continuation are preserved.
Formal M1/M2/M3 qualification remains open; no physical product law or observation
was invented or enabled. The user-authorized U21 directory repair preserved all
three file contents. Full details and exact counts are in section 5.25 and
`test_output/srcMoon/m3-independent-20261118T070304Z/completion_report.md`.
No model/build/test run remains active; root configuration and shared core are
unchanged by this package. No staged files, commit or push were created.

**Current state:** the approved M3 mechanics/rate-service implementation package
and its local/native regression work are complete within their declared scopes.
The current-source M0 baseline passes all ten linked gates with a zero-byte
nightly diff. The two force-aware cut-cell tracing orders, actual PIC packet
backend, constant-history controller and global status/resource handling have
passing supporting evidence. The final shared-core audit preserves all843
historical PASS IDs; its one existing coverage FAIL and one existing timeout
ERROR remain visible. See sections5.24–5.25 for exact final receipts/counts.

**Formal development stage:** M0 baseline preserved; M1/M2 and M3 qualification
remain open under the approved early-development exceptions. Full M3 physical
implementation is still incomplete: qualified recoil/joint kinetic providers,
variable exposure, ion surface fate and asynchronous science sampling/time
contracts remain unavailable. I21-I24 are gating/evidence owners, not completed
physical execution/scoring. No observation or milestone PASS is inferred from
mechanical fixture results. Concrete primary-source proposals are linked below.

The source basis is `74974c2b2fc1b4cf09a73c320064a22370c0ca94` plus uncommitted
working changes. Three agents and the parent implemented/reviewed the work.
All required builds/runs used disposable temporary checkouts; active root
Makefile/site configuration and staging remain unchanged. No commit, push,
reviewed-reference replacement or canonical data write occurred. The detailed
ledger below preserves failed attempts, source/build/data identities, original
deadlines and unresolved limits; earlier dated checkpoints are historical.

**Earlier proposal checkpoint — 2026-11-17, before implementation:**
[test/M3_COMPLETION_PROPOSAL.md](test/M3_COMPLETION_PROPOSAL.md) records the
remaining R3.1–R3.4 implementation, conservative event/mover design, actual
I21–I24 runner contracts, data/physics decisions and one-hour/cluster testing
tiers. Three parallel read-only audits inspected current source, data and
statistical/qualification requirements. All 59 checked data identities match;
this is an identity audit, not 59 science-test PASS results. No production
code, configuration, threshold, reference or roadmap gate changed. No build
or model test ran; formal M1/M2/M3 status remains unchanged. See section 5.17.
The source revision for this proposal is
`74974c2b2fc1b4cf09a73c320064a22370c0ca94`; earlier evidence retains its own basis.

**Previous completed work — upstream core integration, 2026-11-17:**
the working Moon code is integrated with upstream
`5ce92cb2feecadd11a1d60ace69153a7c7da4122`, preserving local work based on
`17d49e00d1d2ccf33b06d11a3e6e5ffef34ac810`. Native tests ran against
`4bff09213b1bfa3f59e113ba8cbb2d56f22d30ef`; the final two upstream commits
change unrelated comet Python/tooling and the Earth nightly launch command.
All **656 native/config source identities are unchanged** at the final head;
the latest nightly tooling separately passes 64/64 tests and its official check.

All **919 historical PASS cases are preserved**: 843 CORE, 62 Enceladus and
14 Moon local aliases. Pristine upstream and the candidate each have 843 PASS,
240 SKIP, one existing strict-coverage FAIL and one existing CME timeout ERROR;
there are no missing IDs, new core failures or status changes. M0 I01-I10 PASS
with unchanged reviewed I04 reference and zero-byte `test_Moon.diff`. Actual
one-/eight-rank Na event budgets and nine-window sampling checks PASS.

The integration retains upstream core implementations and the existing optional
Moon hooks. The M2 builder now derives both chemistry controls consistently,
checks effective compiler macros with actual build flags, and rejects stale
receipts. Its original ERROR is retained. The corrected eight-rank run has
**6 local PASS / 7 formal SKIPPED**, **56 supporting PASS / 2 expected SKIPPED**,
**0 FAIL/ERROR**, including all 36 short numerical cases. Diagnosis, original
attempt, correction, normal-build restoration and retained-ABI checks completed
in **3074.63 s**, within the original 3600-second envelope. Enceladus exporter
paths are reproducible across checkouts; only three provenance-comment lines
changed per example, with every other byte unchanged. Its E3 run completed
10,000 iterations with zero H2O/CO2 ledger residuals. Five additional upstream
mesh/surface targets PASS (34 reports, 643 PASS report lines including summaries).

The active branch remains `master`; no new commit or push was made. The local
SHARE linker setting and other site settings are retained. Generated
`Makefile.test` includes upstream registrations and eleven existing local
include overrides. Active application/build artifacts were not rebuilt during
source delivery; start the next application build by removing `build/` before
normal Config.pl and make -j. No model/test run remains active. Formal M1/M2/M3
qualification and unavailable D02 metadata remain open; this integration is
regression preservation and supporting verification, not observational validation.

Evidence and exact commands/hashes:
`test_output/srcMoon/upstream-core-integration-20261117T185033Z/completion_report.md`
and `.json`; source/adapters/build choices in [CORE_INTEGRATION.md](CORE_INTEGRATION.md).
See section 5.16. Earlier checkpoints retain their original revision/basis.

**Source revision inspected for the integration below:** `5ce92cb2feecadd11a1d60ace69153a7c7da4122` plus preserved uncommitted development and integration changes through section 5.16

**Previous completed work — initial M3 implementation/research, 2026-11-17:** the user
approved early M3 implementation/support before M1/M2 closure, explicitly
configured Maxwellian/Voronov atomic electron rates, and literature proposals
before charge-exchange/molecular daughter-kinetics choices are implemented.
Three agents completed the atomic kernel, literature investigations and full
core audit. The parent completed runtime configuration and native verification.
Analysis/contract: `test/M3_INITIAL_IMPLEMENTATION_CONTRACT.md`. The final
eight-rank development campaign completed in **980.26 s (16.3 min)**, with
**three supporting PASS / zero FAIL/ERROR**, one U16 PASS and four formal M3
SKIPPED results. All ten M0 gates PASS and the nightly diff is zero bytes.
Full CORE preserves **843/843** previous PASS cases in **40.8 min**. No M3
qualification PASS or linked product-budget result exists. M1/M2 formal
status and all original criteria remain unchanged. See section 5.15.
No model/build run remains active. Root Makefile/configuration are preserved;
no commit/push was made. Exact commands, hashes, plots and scope are in
`test_output/srcMoon/m3-initial-20261117T071134Z/completion_report.md` and `.json`.

**Purpose:** living implementation/evidence ledger, not a validation certificate.

**Current follow-up — remaining M2 code and combined verification:** the user
requested completion/testing of unfinished M2 code and explicitly authorized
parallel agents. Finite-disc production terrain forcing, I15 staging/scoring,
review-gated independent December execution, formal I19 dispatch and I20
production export/collocated comparison are implemented. Local protocol checks
pass. The final frozen eight-rank campaign completed in **2242.36 s (37.4 min)**
with **58 supporting PASS / 2 SKIPPED / 0 FAIL / 0 ERROR**, including native
verification of the actual-initial-state and formal-input-preservation fixes.
A subsequent 22.64-second native source-only check verified the requested I19
full initial-state export (486,000 rows), without creating a qualification.
All builds ran in an isolated checkout outside this source tree. The full
CORE audit preserves 843/843 previous PASS cases.
These supporting results do not close M1/M2 or constitute validation.
The user confirms no additional exact-product D02 metadata are available and
approves independent December spin-up with separately explicit UTC bounds.
M1/I13 remain open; D02 and unresolved review/physical-operator evidence keep
formal comparisons SKIPPED. See section 5.14 (appended below),
`test/M2_COMPLETION_CONTRACT.md` and
`test_output/srcMoon/m2-completion-20261116T195511Z/`.
No model/build run remains active; the active configuration is unchanged and
no commit or push was made. Detailed commands, identities, results and limits
are retained in that evidence directory's `completion_report.json` and `.md`.

**Previous November 13 follow-up — approved M2 numerical/execution package implemented and
verified:** the final clean eight-rank development command completed in
1849.678 s (30 min 50 s), including required clean builds, data QA, local and
linked checks, plots and final normal-build restoration. Results: six local
PASS, seven formal M2 SKIPPED; supporting 44 PASS, two SKIPPED; zero FAIL or
ERROR. All 36 short numerical cases PASS. M0 has ten linked PASS and zero-byte
nightly diff. All 123 selected shared-hook core regressions are preserved, and
the additional explicit dt/weight/two-tags-on-one-element check PASSes.
Formal M1/I13 and M2 remain open; no observational validation is claimed.
No model/build run remains active. The active checkout configuration is
unchanged; no commit/push was made. See section 5.13 and
`test_output/srcMoon/m2-implementation-20261113T060631Z/completion_report.json`.

**Previous current-core integration:** its audited upstream basis was
`09fabea55e367f07a530460b488fa7a981ff6176`. The integration work made no
commit or push; the user's later commit records its changes. All 1,193
registered core/application cases are accounted for: every one of 919 baseline
PASS results is preserved, with no missing cases or status changes. Both sides
have 919 PASS, 272 SKIPPED, one pre-existing strict-coverage FAIL and one
pre-existing CME timeout ERROR. This is regression preservation, not an
all-green registry or observational validation.

The user separately approved the exact current-core I04 reference after that
audit. Its strict checker now PASSes with a zero-byte diff, and M0 reconciles to
10 PASS, zero FAIL/ERROR/SKIPPED. Normal-mode one-/eight-rank Na event budgets
and sampling storage PASS; the eight-rank Na/Ar supporting campaign has six
local PASS, seven formal I13-gated SKIPPED, six supporting PASS and one D02
supporting SKIPPED, zero FAIL/ERROR. The final normal Na-only Moon/SPICE build
is restored with TESTMODE=off and -O3 -g. No model/test run remains active.

Durable evidence, source/configuration backups and complete attempt/retry
history: `test_output/srcMoon/core-integration-20261112T224309Z/`. Working
checkouts remain outside the tree under
`/tmp/amps-moon-core-integration-20261112T224309Z-ehx53mh8/`. See
[CORE_INTEGRATION.md](CORE_INTEGRATION.md) and section 5.9. Earlier checkpoints
below retain their original source/configuration basis.

**Previous M2 checkpoint, before core integration:** the user authorized M2 implementation/local and
supporting production-path verification before M1 closure. The bounded
eight-rank development campaign now PASSes U08-U13 and supporting terrain,
surface-distribution and coupled Na/Ar budget checks. Formal I14-I20 remain
SKIPPED without I13. D03 selected inputs are package-READY after 214 official
archive byte matches; D02 remains PARTIAL. Full M0 requalification PASSed after
the shared mover/source changes, as did the restored eight-rank Na chemistry/
event-ledger prerequisite. No model run remains active. See section 5.7.

**Previous chemistry checkpoint:** the November 12 actual-event audit found that the
configured generic Exosphere chemistry processor excludes Na and bypassed the
Moon lifetime callback. Production dispatch is corrected without changing the
existing 130500 s lifetime/shadow law. Linked survival/ledger fixtures PASS on
one/eight ranks; the changed nightly output is reproducible in two four-rank
runs. The user approved the replacement reference; the strict checker now
PASSes with a zero-byte diff. The eight-rank coarsest diagnostic pilot
completed PARTIAL with a closed budget, but no accepted statistical duration.
Full M0 requalification PASSes. See section 5.6; earlier results below
retain their original configuration and cannot qualify the corrected case.

**Resumed checkpoint:** work resumed after the requested safe stop. The
previous checkpoint is preserved in [`HANDOFF_M1.md`](HANDOFF_M1.md). The
`test_Moon` discrepancy was traced to a November 9 update in the shared
SWMF_data repository: this checkout's output exactly matches the October 11
I04 reference already pinned in its acceptance contract. The canonical checker
now uses a reproducible exact copy of that reviewed Git object and verifies
its checksum; `make -j test_Moon` PASSed with a zero-byte diff. Shared external
files and accepted numerical reference values were not changed.

## 1. Current development stage

Current work implements **experimental M2 capabilities**, while **M1 numerical
qualification remains open** under the explicit development-order amendment.
The current-core integration preserves all passing core/application regressions;
M0 requalification has ten linked PASS after the separately reviewed I04
replacement. Normal-mode supporting checks also pass. M1 remains incomplete,
formal I14-I20 remain predecessor-gated, and no later milestone is opened by
regression preservation or supporting evidence.

| Milestone | Formal state | Evidence and interpretation |
|---|---|---|
| M0 linked baseline | **CURRENT CORE REQUALIFIED** | Nine unchanged clean configuration-separated linked gates plus the user-approved, strictly checked I04 yield 10 PASS, zero FAIL/ERROR/SKIPPED. All required local prerequisites PASS. Three four-rank nightly outputs agree exactly; the installed reference is ed875cef… and test_Moon.diff is zero bytes. Original rejected comparisons remain retained. |
| M1 numerical qualification | **INCOMPLETE — active frontier** | I11-I13 runners exist. The approved bounded I11 campaign is 480/240/120/60 km, 500 steps/sample, 5000 steps maximum per level, restart writes off, with unchanged independence/equilibrium/CI/grid gates. No retained evidence establishes I11 PASS; I13 lacks a reviewed frozen passing baseline. |
| M2 lunar environment | **ACTIVE IMPLEMENTATION; EXPERIMENTAL / NOT QUALIFIED** | Supporting eight-rank geometry, surface distributions and Na/Ar gas+surface accounting PASS; U08-U13 PASS. D03 selected inputs READY, D02 PARTIAL. Formal I14-I20 are SKIPPED without same-campaign I13; development results do not complete M2. |
| M3 chemistry | **DEVELOPMENT / NOT QUALIFIED** | The atomic Maxwellian rate service and U16 are implemented. Narrow-scope particle conversions now run through the normal reaction controller for Na photoionization, 4.0 keV/u Na ground-state charge exchange and stationary-⁴He electron impact in a 10 eV bath; their supporting ensemble, PIC and family checks PASS (sections 5.57-5.59). Formal I21-I24 stay SKIPPED: they need same-campaign I13, and molecular kinetics, product kinematics, other bath conditions and measured forcing remain open. |
| M4 drivers/sources | **NOT REACHED** | I25-I31 runners do not exist. Some Na and uniform-Ar implementation is experimental; measured drivers, complete species sources, meteoroid forcing and H2O/OH remain incomplete. |
| M5 observations | **NOT REACHED / NOT VALIDATED** | I32-I40 runners do not exist; qualified operator/data/predecessor gates are incomplete. |
| M6 parameter freeze | **NOT REACHED** | I41 absent. |
| M7 UQ/release | **NOT REACHED** | I42-I46 absent; the new global driver is enabling infrastructure only, not I42 PASS. |

Implementation reach must not be confused with milestone qualification. The
code reaches well into M2 and contains isolated later experiments. The approved
exception permits M2 implementation/supporting tests now; formal qualification
still requires I11, I12 and I13 before I14-I20.

## 2. Roadmap reconciliation

The two required roadmaps are not byte-identical:

- source roadmap SHA-256:
  `eba1571968afff4f948505b77578ee782c4fb54d02ac786ab9ecaa7ea00438d2`
- validation-data roadmap SHA-256:
  `d10eca1c539f9cf82940236d9e1aa36af665d43e82ba475925950d582fbf49c4`

The source copy adds dated implementation status, the user-approved
480/240/120/60 km 200-GB campaign, the 500-by-5000 bounded execution contract,
the non-qualifying 480 km physical-duration diagnostic, diagnostic stage
plots, M2 implementation notes, partial Ar status, and the user-directed
one-hour development / deferred-cluster execution policy and approved
I19/I20 numerical-protocol amendment. The external copy
retains the common scientific roadmap. Therefore:

- the source roadmap controls current code, campaign selection and bounded
  execution;
- the shared/common text controls milestone meanings, physics requirements,
  thresholds, predecessor gates and U/I identities;
- `/data/vtenishe/moon_validation_data/AGENTS.md`, package metadata and native
  labels control data immutability/readiness/units/frames/times.

No difference authorizes a relaxed scientific gate.
The November 12 development-order and finite-release amendments are detailed
in `test/ROADMAP_RECONCILIATION.md`, including prior/current source hashes.

## 3. Test architecture status

| Layer | Implementation |
|---|---|
| U01-U26 | all 26 have individual `test.py`, README and `reference/acceptance.json` |
| I01-I10 | all 10 have individual linked runners/contracts; configuration-separated M0 orchestrator exists |
| I11-I20 | all 10 have individual runners/contracts; I14-I20 declare CSV/PNG/EPS comparison evidence |
| I21-I46 | frozen names are listed, but individual production-linked runners/contracts are not implemented |
| global runner | `srcMoon/test/test_runner.py` audits and delegates roadmap-ordered phases; it records planned gaps as SKIPPED in release scope |
| comparison renderer | `common/write_comparison_artifacts.py` strictly renders evaluated rows only |

Architecture audit executed on 2026-10-30:

```bash
PYTHONDONTWRITEBYTECODE=1 python3 -B srcMoon/test/test_runner.py \
  --phase audit --output-dir /tmp/srcmoon-audit-20261030-2
```

Result: `PASS`; implemented layout U01-U26/I01-I20 complete, required contract
fields present, comparison artifact contracts present for I14-I20, and planned
registry exactly I21-I46. This is an infrastructure audit, not a physics PASS.

## 4. Component status

| Component | Production state | Verification/validation state | Next required gate |
|---|---|---|---|
| parser/configuration/species/mover | operational | M0 verified | preserve through I11-I13 |
| lunar gravity/rotation/differential/non-inertial forces | operational | U03/U04 and I05/I06 verified at M0 | later parallel/long-run reproducibility |
| Na radiation pressure/shadow gating | operational | U05/I07 verified; Combi reference qualified locally | observational use waits for M1/operator gates |
| ion Lorentz mover | dedicated fixture operational | U06/I08 verified fixture | normal science species/driver qualification later |
| sphere boundary/ray tracing | operational | M0 classification and regression verified | I13 frozen controlled baseline |
| LOLA decode/mesh/radial displacement/BVH | implemented | local and demonstrations exist; I14 not formally qualified | M1 PASS, D01 READY, I14 |
| illumination/terrain shadow/Earth eclipse/PSR | implemented in M2 paths | local kernels/demos; I15/I20 incomplete | M1, D02 and linked gates |
| simple/Diviner/dynamic thermal state | implemented experimentally | U10/U11; I16/I17 incomplete | qualified D03, convergence, frozen hold-out scoring |
| impacts/sticking/accommodation/residence/re-emission | Na/Ar implementation exists | U12; I18 predecessor-gated | I13 then I18 |
| inventory/conservation/release/cold trapping | experimental | U13; I19/I20 incomplete | weighted global event/budget closure |
| Na sources | MIV/PSD/TD operational; sputtering configuration-specific | U22; Case-1 demonstration experimental | I29 after preceding gates/drivers |
| He/Ne/Ar/meteoroid/H2O-OH | partial or disabled | not linked-qualified | I25-I31 in roadmap order |
| chemistry | Na lifetime path operational; later networks incomplete | U14 local; I21-I24 absent | M3 after M2 |
| LOS/column/D-line output | generic production output exists | geometry/operator coverage incomplete | U26 extensions, I32/I33 after M4 |
| mission observation operators | incomplete/partially dormant | no advertised observational validation | I32-I40 with READY Dxx packages |

## 5. Latest retained execution evidence

### 5.1 I11 qualification evidence — 2026-10-31

The approved eight-rank 480/240/120/60 km qualification campaign completed
all four bounded production runs.  The first sandboxed launch is retained as
an infrastructure `ERROR` because OpenMPI could not see a network interface;
the unrestricted retry reused the hashed executable from the successful clean
Moon/SPICE build and completed normally:

```text
test_output/srcMoon/I11-qualification-20261031T001923CDT-mpi-retry/
```

Formal result: **I11 FAIL** — 0 PASS, 1 FAIL, 0 SKIPPED, 0 ERROR.  Every level
executed ten 500-step samples (5000 particle steps), kept the external state
exactly frozen, and emitted no particle restart files.  All four levels
reached the cap before the unchanged independence, equilibrium, and confidence
interval gates were met, so no cross-resolution means were accepted.

At the time, the sampled inventory trend was interpreted as physical startup.
The November 11 audit below invalidated that interpretation: stale sampling
windows also contributed, and the old campaign did not record an independent
particle-list inventory. The production timestep halves with each
mesh refinement: 27.4234612, 13.7117306, 6.85586530, and 3.42793265 s.  A fixed
5000-step cap therefore represents approximately 38.1, 19.0, 9.52, and 4.76
physical hours at 480, 240, 120, and 60 km. The refined cases receive
progressively less physical spin-up. None of the old raw means met the
statistical gates, so none is a valid spatial-convergence input. Corrected
sampling and direct particle-list diagnostics must be checked before choosing
a new duration or equilibrium-state contract; weakening thresholds is not an
acceptable remedy.

| resolution | wall time | peak process-tree RSS | retained disk | formal level |
|---|---:|---:|---:|---|
| 480 km | 1438.65 s | 9,735,204 KiB | 20 GB | FAIL |
| 240 km | 1674.51 s | 12,331,880 KiB | 25 GB | FAIL |
| 120 km | 2455.14 s | 17,585,100 KiB | 39 GB | FAIL |
| 60 km | 4191.39 s | 28,798,564 KiB | 71 GB | FAIL |

The aggregate retained directory is 153 GB.  Evidence identities are:

- I11 result SHA-256:
  `2841a30ae38ebb5829af217bea135992e7bc4a092dfa0388d00927d277c4c5ad`
- I11 campaign manifest SHA-256:
  `d8c1e542547d378358e9ab873dd082b1463f56478f268e550a1a64af857b8496`
- aggregate summary SHA-256:
  `b23873377daed811a1c10b0c7f5e5d800d8f74fb4efca2441a69c54b42abbe80`

The per-stage CSV/PNG/EPS figures are explicitly marked diagnostic-only and
`accepted_convergence_input=false`.  They are preserved failure evidence, not
validation plots.

### 5.2 I11 physical-duration pilot — ERROR on 2026-11-03

The formal 500-step-sample/5000-step-resolution I11 contract and all numerical
acceptance thresholds remain unchanged. A separate, explicitly non-qualifying
480 km pilot now exercises the same production mover, sources, surface
callbacks, completed production sampling buffers, frozen epoch, random seed,
and statistical scorer beyond 5000 steps. Its purpose is limited to finding a
candidate common physical duration for later investigator review; it cannot
report I11 PASS, create cross-grid convergence evidence, open I13, or open M2.

To prevent the earlier 20--71 GB per-level visualization output from obscuring
the numerical question, linked I11 now suppresses only completed-sample
volumetric and application visualization file writers. Per-step production
sampling, buffer exchange, I11 metric capture, and all simulated physics remain
enabled. The retained `analysis_diagnostics.csv` identifies the first failed
independence, equilibrium, or confidence-interval gate without changing its
threshold. The pilot runner has no scientific step cap, but terminates the
complete MPI process group at the approved 200 GiB process-tree RSS guard, an
8 GiB host-available floor, or a 100 GiB filesystem-free-space floor.

The clean SPICE Moon build and regression prerequisites passed before launch:

```sh
rm -rf build
./Config.pl -application=moon \
  -spice-path=/home/vtenishe/SPICE/cspice \
  -spice-kernels=/home/vtenishe/SPICE/Kernels
make -j "OPT3=-O3 -g"
make -j test_Moon
PYTHONDONTWRITEBYTECODE=1 python3 -B \
  srcMoon/test/integration-tests/I11/test_infrastructure.py
```

`test_Moon.diff` was zero bytes and all six focused Python infrastructure tests
passed. The eight-rank pilot was then launched with:

```sh
PYTHONDONTWRITEBYTECODE=1 python3 -B \
  srcMoon/test/integration-tests/I11/run_duration_pilot.py \
  --amps ./amps \
  --input srcMoon/examples/m1_memory200.in \
  --mpi-ranks 8 \
  --memory-limit-gib 200 \
  --minimum-host-available-gib 8 \
  --minimum-free-disk-gib 100 \
  --output-dir \
    test_output/srcMoon/I11-duration-pilot-20261101T0145CDT
```

The launch manifest correctly preserves the pre-documentation-edit source
roadmap hash
`808bcfc5858415a174b0685fa0521340c2ff3a26582ea876a2000efe09e1489c`.
The current hash in Section 2 differs only because the dated diagnostic
amendment and reconciliation text were added while this already-built process
was running; no compiled source, input, executable, or acceptance threshold was
changed after launch.

The pilot ended on 2026-11-03 at 01:44 CST with **ERROR**, exit 137: MPI rank
2 received SIGKILL. The logs do not establish who or what sent that signal.
The last printed iteration was 1,101,999 and 2,204 windows were retained.
Wall time was 176,296.93 s, peak aggregate RSS 9.35 GiB, minimum available RAM
206 GiB, and minimum free disk 6.0 TiB. None of the runner's resource guards
was reached. Formal I11 remains `NOT EVALUATED` by this pilot. Result hash:
`6e9f059032c17a1c63f8e2cee0ecfbdf805a9f6189efa93b911352101f17a2b9`.

On 2026-11-11 inspection found a cell-sampling lifecycle defect:
`PIC::Sampling::Sampling()` cleared the new collecting buffer only inside the
species loop guarded by `SaveOutputDataFile[s]`. I11 disabled all these flags,
so alternating stale windows accumulated. The apparent sampled inventory
growth to `4.488e30` therefore cannot determine a physical equilibration time.
This is a sampled-buffer defect; actual particle-list growth was not measured
by that pilot and must not be inferred from its density integral.

### 5.3 Sampling repair implementation and restart checkpoint — 2026-11-11

Evidence root: `test_output/srcMoon/m1-sampling-repair-20261112T021011Z/`.

- Common `src/pic/pic.cpp` now resets new collecting cell storage once after
  the buffer switch, independent of species file flags and output cadence.
  Completed storage remains readable by the existing callbacks.
- A documented I11 storage prerequisite exercises nine deterministic windows
  through the production sampler with all-off, first-on and all-on output
  layouts. Its one-rank SPICE run PASSed with exactly zero completed-marker,
  collecting-buffer and callback errors. The tested build predates the final
  particle-list diagnostic addition; the final eight-rank check remains pending.
- Direct instantaneous Na model-particle count and weighted inventory are
  written to `particle_diagnostics.csv` independently of sampled density.
  This is not a complete source/escape/numerical-loss budget or I19 PASS.
- Optional four-mesh scheduling, aggregate guards and shared I12/I13
  initialization receipts are implemented and documented. Nine isolated
  infrastructure tests PASS; four global-runner tests PASS; architecture audit
  PASS. Actual concurrent MPI scheduling, companion scoring and timing are
  pending and must not be advertised as tested qualification.
- `make -j test_Moon` completed but its scientific comparison FAILed: four
  speed/velocity values differ from the frozen four-rank reference and the
  resulting diff is 927 bytes. Repeating with the original PIC reset placement
  gave the identical microscopic output and diff. This isolates the mismatch
  from the reset relocation, but its underlying cause remains unresolved.
  The target ignores the DiffNum failure and returns zero, so exit status alone
  is not a valid regression PASS. The reference has not been changed.
- An attempted storage check against that no-SPICE nightly control executable
  returned ERROR because it emitted no integration receipt. It does not count
  as a successful regression demonstration against the original implementation.
- The corrected PIC source is restored and the clean Na-only Moon/SPICE build
  completed PASS (`OPT3=-O3 -g`). Executable SHA-256 is
  `c9dd1585f9dbb5e1e4a9e847421f13fd5ae8f8e29615a7ce053e695092da211e`.
  No AMPS process or temporary source rollback remains; final evidence and
  exact source identities are in HANDOFF_M1.md and the checkpoint receipt.

The one-hour campaign, eight-rank storage test, second-species build coverage,
and any required baseline correction are pending at the user's explicit stop.

### 5.4 Resumed verification — 2026-11-11

- The reviewed I04 reference hash remains
  `abc1d9c91f61a1c44e06e98bde8739b2de1650b7a6caef783bb97081f2eb80b4`.
  Its authoritative SWMF_data revision is
  `7a1c437369a600bd248bcf8846b66b286e7109c1`. The current shared file hash is
  `9a0bd0f9dca19cba820565990de846d6acc0c1efdaf3fea97c3d4a2f050cec3a`,
  changed by external commit `6d2644aa2994beaedc9fc76c97c10bba3abb99ac`
  for newer source numerics. It remains untouched.
- `MakefileTest/Table` now calls the checksum-aware I04 checker. The source
  copy is reproduced by `check_reference.py --import-from-repository ...`,
  using `git show` and the pre-existing frozen hash. Three reference-identity
  tests PASS. Canonical regression PASS is preserved in
  `build/test-Moon-pinned-reference.log`, `.diff` and `.json` under the work
  package evidence root. Root `test_Moon.diff` is zero bytes.
- The final eight-rank storage check initially FAILed because the test adapter
  used `GetCenterNode(0)`, a ghost slot rather than an interior cell. That
  adapter is corrected to the three-coordinate interior overload. The
  production reset itself retains the interior-cell loop and is unchanged.
- A separately recorded Na+Ar storage build uses the same documented species
  source derivation as the M2 builder with explicit `purpose=sampling`; it
  makes no I19 chemistry-contract claim. One-rank and eight-rank linked storage
  checks now PASS with two compiled species and nine windows, and exactly zero
  marker, reset and callback errors. Source-qualified fixtures use the existing
  documented positive Ar source to initialize its weight, spherical geometry,
  frozen epoch and simple-cosine temperature; no particles are advanced.
- Canonical Na-only SPICE restoration PASSed. The final executable SHA-256 is
  `b2a200158c540f91481aa75a2f7469b57f9aeb39df936525011e9516135e15cb`.
  Na-only storage checks on one and eight ranks both PASS with exact-zero
  errors. U03/U04/U05/U07/U14 are 5 PASS, zero FAIL/SKIPPED/ERROR. Linked I10
  PASSes with observed orders 1.9999167, 1.9999591 and 1.9999798 and finest
  scaled-state error `3.463247343503966e-9`.
- A real 32-rank deadline smoke test ran four isolated eight-rank meshes,
  stopped after its declared 20 s launch budget, and preserved four ERROR
  child records plus `execution_termination.json`. This is expected guard
  evidence, not physical qualification; MPI launches were reaped.
- Final infrastructure coverage now has ten PASSes, including species-builder
  collision preservation, plus five I04 reference/error-classification checks.
  The exact build-helper source is exempted from the broad `*build*` ignore
  rule so its reproducible species derivation is visible to Git. Existing
  colliding temporary inputs remain intact on builder refusal.
- The first three corrected coarse-grid density integrals match the old pilot
  for windows 1 and 2. Old window 3 minus corrected window 3 equals window 1
  within approximately `1e-12` relative roundoff, as expected for the stale
  alternating-buffer defect. The new direct particle-list inventory also
  increases during startup; that independent trend is real but does not, by
  itself, establish a source/sink budget or its eventual equilibrium time.
- The bounded M1 campaign completed in
  `test_output/srcMoon/m1-sampling-repair-20261112T021011Z/m1-parallel-qualification/`.
  It uses four meshes with eight ranks each, a 3600 s parallel-launch budget,
  the ordinary 500/5000 particle-step policy and the unchanged independent
  scorers. Aggregate wall time was **2228.11 s (37.14 min)**, with peak aggregate
  RSS **64.35 GiB**, minimum host-available memory **178.90 GiB**, and
  **1,281,568,768 bytes** retained at completion. No resource guard fired;
  every MPI launch exited normally and no AMPS process remains.

Final M1 campaign counts: **3 recorded IDs: 1 PASS, 1 FAIL, 1 SKIPPED, 0 ERROR**.

| ID / resolution | Formal result | Evidence / first failed gate |
|---|---|---|
| I11, 480 km | FAIL | 5000 steps, ten samples; 28/1175 metrics fail independence at the final evaluation. |
| I11, 240 km | FAIL | 5000 steps, ten samples; 16/1175 metrics fail independence. |
| I11, 120 km | FAIL | 5000 steps, ten samples; 14/1175 metrics fail independence. At nine samples independence passed but 3195 equilibrium checks failed. |
| I11, 60 km | FAIL | 5000 steps, ten samples; 17/1175 metrics fail independence. |
| I12 | PASS | Independent formula error at most `2.4253192047278088e-16`; 60/30 km targets realize 45.7058/22.8529 km cells. The 120-degree transition uses the predeclared neighbor-balance exception; its parent meets the unchanged 25% gate. |
| I13 | SKIPPED | Same-campaign I11 is not PASS; no reference was automatically frozen or scored. |

The campaign's frozen external-state drifts are exactly zero; all four stages
kept restart writes off and emitted scored-status diagnostic CSV/PNG/EPS/JSON
sets. No `qualified_metrics.csv` or accepted cross-grid convergence table was
created. Late coarse-grid instantaneous inventories still change by 3.60% and
4.28% between the last two samples; the window means change by 3.83% and
4.52%. These are startup diagnostics, not accepted steady means or a closed
physical budget. M1 is **INCOMPLETE** and M2 remains **EXPERIMENTAL**.

OpenMPI 4.0.3 emitted TCP connection warnings while all four processes kept
advancing; their cause is unresolved and stderr is retained. No MPI abort or
nonzero child exit occurred. This does not turn the evaluated scientific FAIL
into a runtime ERROR or a PASS.

The active campaign loaded the earlier scheduler before the progress-label
refinement. Its early coarse-stage plots with placeholder SKIPPED labels are
preserved under `provisional-before-scoring/`; its final four plots all say
FAIL. Future runs preserve separate `stage_diagnostic_pending.*` previews,
whose PENDING label is excluded from result counts. The final independent
particle-diagnostic hash indexing was also added after launch; the completion
audit records those hashes as supplemental evidence without rewriting the
completed campaign manifests. Compiled sources, executable, physics, inputs
and numerical/statistical thresholds were stable throughout the campaign.

The read-only closure script and supplemental hash index are retained as
`complete_work_package.py` and `completion_audit.json` in the work-package
root. The audit binds final campaign artifacts, particle diagnostics, eight
SPICE files, executable/configuration identity, source snapshots and the
latest focused test cohort. Its identity-audit PASS is infrastructure evidence
and does not change the I11 FAIL or the M1 milestone state.

Evidence SHA-256:

- aggregate summary:
  `84f18743314b0770d71fe0bf2680879a2b6034f20f165b98ce1025a53ddcc6c9`
- I11 result:
  `6b9f2bce24cc5cfb11020f8ef6fec991edbb9e2e0848befc257e21e48896293f`
- I11 campaign manifest:
  `82dff6671a1d3f4dc9c5798dd6d263e63b737e48373df1a9820df6a16022b43b`
- I12 result:
  `9c6113f4db1f1fb41f5733913c019a081297fc8ce07c31a1cf5ed85761143076`
- I13 result:
  `1cc03213ee7b7c61f8629997c416e020bea5815fe5fc0d414b16ee6ede230597`

The principal production changes are in `src/pic/pic.cpp` (collecting-buffer
lifecycle), `MoonIntegrationTests.{h,cpp}` and `main.cpp` (thin storage adapter,
direct particle diagnostics and companion receipts). The execution changes
are in `test/m1_execution.py`, `test/integration_common.py`,
`test/run_all_test.py`, `test/integration-tests/I11/run_sampling_check.py` and
`test/build_m2_qualification.py`. I04's checksum-aware comparator and negative
infrastructure tests are registered through `MakefileTest/Table`.
Corresponding PIC/Moon/test READMEs and contracts document units, scope and
failure behavior. The generated root Makefile is excluded from source changes
to commit. Existing unrelated working-tree changes are preserved.

Principal commands (from the AMPS root; build logs are under the evidence
root above):

```sh
# The nightly target configures its own reference case:
rm -rf /data/vtenishe/Moon-test/AMPS/build
make -j test_Moon
# Restore the canonical Moon/SPICE application before linked I11-I13:
rm -rf /data/vtenishe/Moon-test/AMPS/build
./Config.pl -application=moon \
  -spice-path=/home/vtenishe/SPICE/cspice \
  -spice-kernels=/home/vtenishe/SPICE/Kernels
make -j 'OPT3=-O3 -g'
# The final checker rerun uses the preserved nightly output without rebuilding:
make -j test_Moon_check

python3 -B srcMoon/test/run_all_test.py --skip-build \
  --only I11 --only I12 --only I13 --campaign-mode qualification \
  --m1-input srcMoon/examples/m1_memory200.in --mpi-ranks 8 \
  --i11-parallel-levels 4 --test-runtime-seconds 3600 \
  --output-dir test_output/srcMoon/m1-sampling-repair-20261112T021011Z/m1-parallel-qualification
```

The existing campaign path is evidence and must not be reused for another
execution; select a new path when reproducing. The clean final compiler was
g++ 9.4.0 with `-O3 -g`; DE421 CSPICE was enabled, OpenMP=1 and seed=731993.
Every child manifest records the eight SPICE file hashes, derived input hash,
configuration hashes, dirty source revision and executable identity. No Dxx
observational package is used by this numerical work package. Required
observational data and their readiness remain gates for later milestones.

### 5.5 Demonstration evidence

The 8-rank Tenishev-2013 Case-1 demonstration starting
`1991-09-30T00:00:00Z` completed 5000 iterations and ten 500-step samples.
Its bounded execution result is `PASS / EXPERIMENTAL`, not validation:

```text
test_output/srcMoon/sodium-exosphere-demo/
  case1-19910930T000000Z-20261017-retry1/
```

It preserves 40 required density/column/D-line artifacts, uses about 33 GB,
and records I11/I13 `NOT RUN` and I32 `SKIPPED`. The run advanced to
`1991-09-30T19:02:39`; midnight is its start epoch, not a fixed-time snapshot
or an independently established observation time.

LOLA surface and environment demonstration runners also exist under
`srcMoon/test/`. Their successful execution and visualizations remain
separate from I14-I20 qualification.

### 5.6 Actual-event audit and Na chemistry correction — 2026-11-12

Evidence root: `test_output/srcMoon/m1-budget-audit-20261112T064009Z/`.

- Wrote `test/integration-tests/I11/BUDGET_CONTRACT.md` before implementation.
  It declares production entry points, global Na/sphere scope, actual-event
  count/weighted oracles, units, timing, fixture controls, statistical criterion,
  resource guards, artifacts, data and honest status meanings. No test ID,
  observation selection, physical coefficient or formal 500/5000 limit changed.
- Added a null-by-default `PIC::ParticleEvents` observer. Accepted Exosphere
  allocations, adsorption/reflection, photochemical parent/product events,
  selected first/second-order mover escapes/known numerical removals, and
  signed prescribed-reservoir supply now have read-only diagnostic hooks.
- Added `MoonParticleBudget.{h,cpp}`. I11 independently sums allocated owned
  particle lists and compares with cumulative actual events, including
  individual weights. Count residual must be exactly zero; weighted residual
  uses the declared long-double summation bound. Reservoir arrays are
  replicated and counted once; their external replenishment is not a gas birth.
  Snapshots preserve before-step sampling and a final post-step closing row.
- Found the Na chemistry bypass: canonical input selected
  `Exosphere::ChemicalModel::PhotochemicalModelProcessor`, which uses its
  private `TotalLifeTime`/H-O species table instead of the configured Moon
  lifetime. That table excludes Na. Canonical input now selects
  `Moon::Chemistry::ProcessParticle`; non-Na neutral species retain generic
  dispatch. Na+ production weights/kinetics remain explicitly unsupported.
- The new actual configured-dispatch fixture compares 4096 sunlit particles
  with p=1/2 at dt=tau*ln(2) and exact shadow survival, using a predeclared 95%
  Hoeffding bound. Its generic-processor control retains all 8192 particles.
  On eight ranks corrected sunlight retains 2062/4096, shadow 4096/4096; on one
  rank sunlight retains 2050/4096 and shadow 4096/4096. Both chemistry and
  20-step normal transport ledgers independently close with count residual 0
  and no numerical removals. Emission, adsorption and escape are observed.
- The initial eight-rank fixture FAIL is preserved in `budget-8/`: the fixture
  called the chemistry walker before the first ordinary step initialized its
  block table and frozen ephemeris. Its setup was corrected; the successful
  rerun is `budget-8-initialized/`, with `budget-1-restored/` for one rank.
- Clean Moon/SPICE `make -j 'OPT3=-O3 -g'` restoration PASSed. Current executable
  SHA-256 is `e82a204cbcb2741f473c3474c4479611ba11fceef17f3b2cecd4cfc9c1ad6780`.
  Build logs preserve the initial parser/compile errors and successful fixes;
  every recompile first removed only the root `build`. Generated Makefile is
  not a source change to commit.
- U03/U04/U05/U07/U14/U22: **6 PASS, 0 FAIL/SKIPPED/ERROR**. Independent
  statistics analytical/parser checks: **5 PASS**. Ten I11 infrastructure
  checks also PASS. Four global-runner self-tests PASS; the final architecture
  audit PASSes its registry/contract inspection with zero physics test records
  executed. These infrastructure results are not observational validation.
  The independent processor reconstructs the previous
  coarse run's 28 failing independence metrics exactly; its plots are
  diagnostics, not accepted convergence input.
- `make -j test_Moon` completes its production run but its strict comparison
  FAILs with four numerical differences and a 1035-byte diff. Two corrected
  four-rank runs agree byte-for-byte at candidate SHA-256
  `4bb811b3e05d72fa708359f1bea0344efd1a9f7e89609d7e8a3ca7eade808cba`.
  At that failed comparison the reference was still
  `abc1d9c91f61a1c44e06e98bde8739b2de1650b7a6caef783bb97081f2eb80b4`.
  Candidate/review JSON is under `reference-review/`; raw former nightly
  output, strict diff, both corrected runs and configurations are preserved.
- The user then approved the exact corrected I04 candidate. The current
  bundled reference is `4bb811...08cba`; the pre-correction reference is
  archived unchanged. Its reviewed provenance is pinned at
  `2facd9d35ef2017e7d0e6c8ab7da19e1f055e70e42957f9bf85e530231f2bfb0`.
  Six reference-integrity/error checks and `make -j test_Moon_check` PASS;
  root `test_Moon.diff` is zero bytes. No tolerance or I13 baseline changed.
- The pilot launcher now executes both linked prerequisites, prints live
  progress, and requests a safe stop between full steps at an explicit wall
  budget (default 3600 s, 30 s grace). A real eight-rank one-second-budget
  smoke run stopped after two completed steps and wrote PARTIAL, a closing
  zero-residual ledger and receipt. It is operational evidence, not I11 PASS.
- The corrected 480 km eight-rank pilot completed **PARTIAL** after 3601.63 s,
  26022 completed steps, 52 samples, and 9.245 GiB peak process-tree RSS.
  All 54 ledger checkpoints independently PASS. Maximum count residual is 0;
  weighted residual/bound is at most 0.000724; numerical removals are 0.
  Total accepted model-particle events are 4701685 emissions, 3053793
  adsorptions, 591809 photolytic losses and 939778 escapes, leaving 116305
  particles (2.69511538e27 neutral Na atoms). This closes gas events only,
  not the full I19 evolving surface/multispecies system.
- At the final 52-sample evaluation, 414, 11 and 30 of 1175 observables fail
  independence at batch sizes 1, 2 and 4. The independent Python analysis
  exactly reconstructs all three reached production gate results. At batch
  size 2 the global inventory lag is 0.505041 versus bound 0.384388; the other
  10 failures are density samples. Independently computed later gates also
  fail (2414 equilibrium checks and 760 CI metrics at that batch size);
  those explanatory diagnostics do not override production short-circuiting.
  No all-zero science series was present. No candidate duration was produced.
- `duration-pilot/480km/diagnostics/` and `stage_diagnostic.*` contain the
  retained CSV/PNG/EPS figures and identities. A separate `diagnostics-rates/`
  reprocessing preserves the original figures while displaying weighted
  source/sink rates using actual elapsed checkpoint intervals, including the
  short final interval. Output occupies approximately 435 MiB before this
  additional figure set. Formal I11 remains NOT EVALUATED; I13 is unfrozen.
- Full M0 requalification under `m0-requalification/` completed **PASS**:
  I01-I10 **10 PASS, 0 FAIL, 0 SKIPPED, 0 ERROR** and all nine required local
  prerequisite IDs PASS. Initial neutral-SPICE U01-U26 records 14 PASS and
  12 SKIPPED; U06 then PASSes in the dedicated ion build. Across the appropriate
  configurations, 15 distinct U IDs PASS and 11 later-capability IDs remain
  SKIPPED, with no U FAIL/ERROR. Six reference-integrity tests also PASS.
  I10 observed orders are 1.999917/1.999959/1.999980; its finest scaled state
  error is 3.46325e-9. I08 successive error ratios are about 3.966/4.000/4.000.
- Every configuration starts from a removed `build/`, uses
  `make -j 'OPT3=-O3 -g'`, and retains its binary/logs. Canonical Moon/SPICE
  restoration PASSed; root `amps` again has SHA-256 `e82a204c...ad6780`,
  identical to the retained pilot executable. The new nightly microscopic
  output exactly matches approved reference `4bb811...08cba`, and root and
  retained I04 `test_Moon.diff` are zero bytes. The campaign's own nightly
  directory preserves the earlier root nightly evidence.

No Dxx observational package is used by these diagnostics. The same eight
SPICE files and both roadmap hashes are recorded by manifests. The complete
multispecies/evolving-geometry I19 budget and D04 chemistry qualification remain
later gates; the local legacy lifetime and its wiring are not a rate validation.

Reproduction commands for this work package (from the repository root, using
a new output directory for each run):

```bash
python3 -B srcMoon/test/integration-tests/I11/run_budget_check.py \
  --amps ./amps --input srcMoon/examples/m1_memory200.in --mpi-ranks 8 \
  --output-dir test_output/srcMoon/budget-check-NEW
python3 -B srcMoon/test/integration-tests/I11/run_duration_pilot.py \
  --amps ./amps --input srcMoon/examples/m1_memory200.in --mpi-ranks 8 \
  --wall-seconds 3600 --output-dir test_output/srcMoon/duration-pilot-NEW
python3 -B srcMoon/test/integration-tests/I11/analyze_diagnostics.py \
  --level-dir test_output/srcMoon/m1-budget-audit-20261112T064009Z/duration-pilot/480km \
  --output-dir test_output/srcMoon/independent-diagnostics-NEW
python3 -B srcMoon/test/run_integration_tests.py \
  --output-dir test_output/srcMoon/M0-requalification-NEW
```

The pilot executable is retained as `amps-budget-spice` within the evidence
root; its source/config identities are in `source-snapshot/` and linked
manifests. Compiler is g++ 9.4.0 with `-O3 -g`; MPI is Open MPI/OpenRTE 4.0.3,
pilot MPI layout is eight ranks with OpenMP=1 and seed 731993. No particle
restart is written. Full M0 additionally runs local U08 against READY D01:
`LDEM_4.IMG` SHA-256 is
`972d66d7f330abca62abc21aae270661516c00dd1a881ec102b94df9563656b5` and
`LDEM_4.LBL` is
`b5a197b4d417e3b5e0e84f8fe15f9a1d2203b8eb5d7213236c0ffbc86e940cc6`.
D02 remains PARTIAL. No new observational scoring or external-data edits occur.

Production files changed in this audit are `input/moon.input`,
`src/pic/{pic.h,pic.cpp,pic_mover.cpp}`,
`src/models/exosphere/{Exosphere.cpp,Exosphere_Chemistry.cpp,Exosphere_Parallel.cpp}`,
and `srcMoon/{Moon.cpp,Moon.h,MoonParticleBudget.cpp,MoonParticleBudget.h,`
`MoonIntegrationTests.cpp,MoonIntegrationTests.h,main.cpp,makefile}`.
Supporting changes include the I11 budget/statistics/pilot runners and
contracts, configured integration/global runner wiring, the I04 reviewed
reference/integrity/review tooling and canonical `MakefileTest/Table`, U14
scope documentation, and relevant README/RUN/HANDOFF/reconciliation/gap files.
Existing unrelated dirty source remains preserved; root `Makefile` is generated
and must not be committed. No roadmap hash, permanent ID or acceptance
threshold was changed by this audit.

### 5.7 Approved M2 development/supporting verification — 2026-11-12

Work root: `test_output/srcMoon/m2-development-20261112T084925Z/`.
Source HEAD remains `22ad748432ce032775189ba234be7a3d339a4d8d` plus preserved
dirty/untracked development. No commit/staging/reset was performed. The source
roadmap and AGENTS now explicitly permit M2 implementation/local and supporting
verification before M1 closure. Formal qualification, selection, thresholds
and I13 review remain unchanged. Source roadmap SHA-256 is
`2761b447b46539edb21395550cd1b5d9ff6add4839b0db5a7e6d734cf2a4f02a`;
external roadmap remains
`d10eca1c539f9cf82940236d9e1aa36af665d43e82ba475925950d582fbf49c4`.
The reconciliation document records both controlling scopes and hash history.

**Implemented production changes and scientific contract.** The user approved
continuous finite-inventory release. `FiniteSurfaceRelease.h` integrates
competing constant first-order thermal/PSD hazards analytically; the ordinary
production source generators supply their unchanged rates, velocity and angle
laws. Full and fractional weighted packets represent released atoms, with
conditional exponential departure times. Analytic remaining stock is stored
once on every replicated rank, without a second desorption-queue debit or
fictitious replenishment. Incoming adsorption stays pending for normal
exchange. Rank-zero draws are broadcast and exactly one final AMR owner
allocates. Individual weights ON, global species steps and OpenMP=1 are
explicit requirements; unsupported finite-stock sputtering fails with a
diagnostic instead of being converted into an assumed first-order law.
This is a continuous-reservoir quadrature, not atom shot noise. Temperature/
illumination are constant over the split release interval; temporal/weight
convergence and within-step deposition accuracy remain open linked gates.

The new `MoonM2Budget` observes actual accepted source, adsorption, surface
deposit/release and mover events, independently sums owned gas and replicated
surface stock counted once, includes pending fluxes, and records counts,
atoms and kilograms. Independent Decimal reconstruction checks gas, surface
and combined identities. Arithmetic error bounds account for actual double
surface updates and long-double accumulation; the 0.5% physical criterion is
unchanged. Physical escape uses the existing I09 instantaneous lunar two-body
energy after LSO-to-inertial velocity correction. Bound finite-box exits are
numerical loss, not physical escape; this is not ultimate three-body escape.

Supporting geometry tests exposed two production defects: strict BVH slab/
barycentric bounds missed shared vertices, and the configured force SecondOrder
mover bypassed the moving-surface callback. Arithmetic enclosures fix the first;
the actual selected entry now dispatches the boundary tracer. Collision search
uses the same accelerated second-order position extension as the mover and a
shortened midpoint velocity update. Source terrain placement and shortened
flight now share the sampled departure epoch and exact SPICE transform.
The I14 fixture initializes the ordinary orbital/force state through one
normal step before probing; it does not clamp Na radiation-table input or
replace the physical force.

**Retained latest evidence.** `campaign-5/test_summary.json` and its
`m2-development/development_summary.json` distinguish:

| Scope | PASS | FAIL | SKIPPED | ERROR |
|---|---:|---:|---:|---:|
| Registered U08-U13 and formal I14-I20 | 6 | 0 | 7 | 0 |
| Supporting data/production checks | 6 | 0 | 1 | 0 |
| Requalified I01-I10 | 10 | 0 | 0 | 0 |

The supporting PASSes are D01 readiness, D03 readiness, frozen I17-chain
integrity, and production I14/I18/I19 scopes. D02 is SKIPPED/PARTIAL. No
supporting PASS is counted as a formal I-series PASS.

- I14: 513/513 radial hits; maximum independent DEM control-point error
  `4.656612873077393e-10 m`, ray/vertex difference
  `9.313225746154785e-10 m`; zero degenerate/inverted faces; all seven ranks
  owning eligible facets among eight ranks report actual mover collisions.
- I18: 32768 Na trials and 32768 Ar samples; 17042 Na re-emissions against
  expected 17060.33 (frozen 95% half-width 177.25). Four angular/energy KS
  distances lie below unchanged limits; Ar source is `3.8e20 atoms/s`.
  Desorption relative errors are `3.9774e-15` (Na) and zero (Ar).
- I19: 20 synchronized ordinary steps, 42 checkpoint/species rows, exact count
  closure, maximum relative imbalance `1.066445872840518e-15`, no numerical
  removals. Final Na/Ar actual emission counts are 116/150, adsorption 39/65,
  surface release 78/98 and escape 1/1. The fixture explicitly seeds 32
  model-weight surface equivalents/species, 64 incoming particles/species and
  one energetic outgoing particle/species before the ledger initial row.
  These are verification initial conditions, never synthetic observations.
- U09 adds a closed-icosphere shared-vertex oracle, U11 exact isothermal
  radiative equilibrium at 1/60/3600 s, and U12 competing/stiff/semigroup
  finite-release cases. All U08-U13 PASS through production kernels.
- I19 infrastructure: five PASS (exact internal transfers, false receipt,
  unexplained gas atoms, duplicated surface stock, distinct build purposes).
  Global-runner infrastructure: five PASS, including subcases proving that
  supporting PASS stays separate and FAIL/ERROR controls aggregate exit status.
  I11 infrastructure's ten checks also PASS; these are not extra physics IDs.

**Build/baseline preservation.** Every reconfiguration removes root `build/`
first, uses normal Config.pl and `make -j 'OPT3=-O3 -g'`, and retains logs and
executables. The dedicated chemistry-OFF/individual-weights-ON Na/Ar executable
is `build-8/amps-moon-m2-na-ar`, SHA-256
`5d838d0898828d811a568c3889c2720b249e2cbbc43f73783fe8dd982d1c9a11`.
The configuration-separated M0 campaign includes SPICE, no-SPICE and Na/Na+
builds; all nine required local prerequisites PASS. Its first neutral U-series
has 14 PASS/12 SKIPPED, and the dedicated ion U06 adds the fifteenth distinct
passing U ID. I01-I10 are 10 PASS, zero FAIL/ERROR/SKIPPED.
`make -j test_Moon` (with campaign-owned nightly directory and debug compiler
override) produces a zero-byte root `test_Moon.diff`; reviewed reference/output
SHA-256 remains `4bb811b3e05d72fa708359f1bea0344efd1a9f7e89609d7e8a3ca7eade808cba`.
Canonical Moon/SPICE is restored; root executable SHA-256 is
`91619a05443336c42f738395f081b9ffea4584b4a7e89028b63ed8e9a13bb7f9`.
The post-restoration eight-rank I11 event prerequisite PASSes with 2062/4096
sunlit survivors, 4096/4096 shadow survivors, exact count residual zero,
3478 accepted emissions, 198 adsorptions, one escape and zero numerical loss.
This does not re-evaluate formal I11 convergence.

**Data readiness, immutable raw evidence and remaining facts.** D01 IMG/LBL
hashes remain `972d66d7f330abca62abc21aae270661516c00dd1a881ec102b94df9563656b5`
and `b5a197b4d417e3b5e0e84f8fe15f9a1d2203b8eb5d7213236c0ffbc86e940cc6`.
The D02 processor produced 366 spatial rows (SHA-256
`7c26d67be6a6b98e27c13f05a560607bea76ad03210d2820fb97caa39ed6af14`),
but actual 18.6-year calculation epochs and the matching extended-disc
temporal comparison remain unresolved. Native 2009-2013 product dates cannot
silently supply those epochs; the pending factual question includes code in
`question.md`. D02 remains PARTIAL.

Official calibrated and derived Diviner manifests match all 214 selected raw
files. `verify_diviner_archive.py` preserves the authoritative evidence;
`qualify_d03_package.py` reproduces scoped READY metadata from evaluated raw
and processor QA. Published D03 metadata and processed GCP/PRP products are
under `/data/vtenishe/moon_validation_data/diviner/`. No raw file, frozen I17
selection or existing I17 history was changed; December residuals were not
examined. Calibration release 3.1 is distinct from native product-label version
1.0. READY covers 18 GCP, two PRP and 87 frozen RDR pairs only. GCP/PRP
uncertainties remain absent where undocumented, GCP brightness temperature is
not assumed to equal kinetic temperature, and PRP 2-cm mean versus surface
maximum temperatures remain distinct. Processed science SHA-256:

- GCP table: `4a1c4eed1e1004c45848d6966bc58fd988b2da4236c67667f14ab3c0ab1c106d`.
- GCP independent points: `9e3d2635303cb6065152de2606115518506749c038c2621a3005f5563b10abba`.
- PRP reference: `804c748e21e4b4f651fbd578d44c9c7fffa157a3a2815ffec9f281cd57e66af3`.

**Reproduction and artifacts.** The ordinary fresh-build command is:

```sh
python3 -B srcMoon/test/test_runner.py --phase m2-development --mpi-ranks 8 \
  --output-dir test_output/srcMoon/M2-development-NEW
```

The recorded latest campaign explicitly reused the matching retained M2
build, before M0 reconfiguration:

```sh
python3 -B srcMoon/test/test_runner.py --phase m2-development --mpi-ranks 8 \
  --m2-amps test_output/srcMoon/m2-development-20261112T084925Z/build-8/amps-moon-m2-na-ar \
  --output-dir test_output/srcMoon/m2-development-20261112T084925Z/campaign-5
python3 -B srcMoon/test/run_integration_tests.py \
  --output-dir test_output/srcMoon/m2-development-20261112T084925Z/m0-requalification
python3 -B srcMoon/test/integration-tests/I11/run_budget_check.py \
  --amps ./amps --input srcMoon/examples/m1_memory200.in --mpi-ranks 8 \
  --output-dir test_output/srcMoon/m2-development-20261112T084925Z/canonical-budget-regression
```

Do not reuse those nonempty output paths. Manifests preserve exact inputs,
SPICE/data/source identities, seed 731993, eight ranks/OpenMP=1 and resource
histories. I14/I18 exact comparison CSV/PNG/EPS and I19 budget CSV/PNG/EPS are
under `campaign-5/m2-development/supporting/`. Operational monitoring enforces
3600 s per linked launch, 200 GiB process-tree RSS, 8 GiB host-memory floor and
100 GiB free-disk floor; these are execution guards, not scientific tolerances.
RUN.md and data_processing/README.md give complete reproducible commands.

Earlier failures remain under `campaign`, `campaign-2` through `campaign-4`,
`shared-vertex-before` and GDB logs. They retain evaluated ray misses/mover
bypass, fixture initialization ERROR, source-epoch ERROR and residual-receipt
precision ERROR; none was relabelled PASS. The packet-weight build requirement
and analytic-stock representation fix their identified mechanisms. The latest
run evaluates all supporting criteria afresh.

Production changes span generic PIC event/mover code, Exosphere finite release
and source/surface exchange, Moon terrain/source adapters, Ar ownership and
the new Moon budget/fixture. Test changes include production probes/contracts,
the dedicated builder, global development runner and data-qualification scripts.
Relevant README/RUN/HANDOFF/gap/question files are updated. `source-snapshot/`
and `completion_report.json` bind final source/build/data/result identities.
Root `Makefile` is generated and must not be committed. M1 remains INCOMPLETE;
formal M2 is EXPERIMENTAL / NOT QUALIFIED, with I14-I20 SKIPPED until I13 and
their individual production/operator/convergence gates pass.

### 5.8 Canonical test completion reporting — 2026-11-12

The user's 16:00 canonical run completed PASS with a zero-byte diff, but the
custom check recipe omitted the traditional completion line and the checker
self-tests printed expected FAIL/ERROR fixture messages into its terminal
output. This was an infrastructure reporting defect, not a new physics failure.

`integration-tests/I04/test_reference.py` now uses unittest output buffering:
passing invalid-input fixtures stay quiet; a genuinely failing test replays
its captured diagnostics and retains its assertion/nonzero status. The separate
real production checker remains unbuffered. `MakefileTest/Table` restores
`Test test_Moon has finished: ...` after successful comparison and the diff
listing, and silences the explanatory recipe comment. The generated
`Makefile.test` was rebuilt with the existing C++ table parser; its only changes
are those two reporting lines. Root `Makefile` was not edited.

Evidence: `test_output/srcMoon/I04-reporting-20261112T221625Z/` preserves the
previous recipe/diff/receipt, generation/build logs and exact generated diff.
Source HEAD is `22ad748432ce032775189ba234be7a3d339a4d8d` plus the recorded
working-tree changes. Application physics/build configuration, both roadmaps,
reference and provenance hashes are unchanged. Commands:

```sh
python3 -B srcMoon/test/integration-tests/I04/test_reference.py
make -j test_Moon_check
```

All six infrastructure self-tests PASS, with no misleading FAIL/ERROR lines.
The real production comparison PASSes against the previously completed
four-rank `run_test_Moon` output; `test_Moon.diff` remains zero bytes and the
final completion line is present. A separate deliberately failing assertion
probe confirms that buffering still displays failure diagnostics. This is a
reporting check and reuse of completed production output, not a new simulation,
M0 campaign, or milestone qualification. I04 README is updated; no merge,
commit, reference replacement or application rebuild was performed.

### 5.9 Current-core integration and regression preservation — 2026-11-12

The user authorized implementation without commits or pushes, preferring the
current core and retaining necessary shared extensions only when regressions
remain passing. Master fast-forwarded from `22ad748432ce032775189ba234be7a3d339a4d8d`
to existing upstream `09fabea55e367f07a530460b488fa7a981ff6176`. The initial 104
source/configuration files, tracked patch and local settings were backed up;
all local development remains unstaged and generated Makefile is not staged.
Temporary checkouts are outside the tree. The final evidence directory is
`test_output/srcMoon/core-integration-20261112T224309Z/`.

The integration keeps upstream Poisson waiting times, cumulative source
selection and corrected source-angle samplers. Necessary shared extensions
retain actual-event diagnostics, output-independent sampling reset, optional
moving departure epoch/force route, runtime source/stock policy and the already
approved finite-inventory injector. A NaN departure scope is now truly inert in
legacy parallel sources. No physical coefficient or acceptance tolerance was
changed. CORE_INTEGRATION.md inventories the source files/entry points and
contains exact build/registry/audit commands. Production C/C++/Fortran/header/
definition/library-makefile bytes agree across the tested candidate and the
core/application checkouts (653 compared files); focused Python changes and
metadata/source hashes are recorded separately.

| Registry namespace | PASS | FAIL | SKIPPED | ERROR | Result |
|---|---:|---:|---:|---:|---|
| Core, 1,085 cases | 843 | 1 | 240 | 1 | Exact baseline dispositions preserved. |
| Applications, 108 cases | 76 | 0 | 32 | 0 | Exact baseline dispositions preserved. |
| Combined, 1,193 cases | 919 | 1 | 272 | 1 | All 919 baseline PASS identities preserved; zero missing/new failing cases or status changes. |

The remaining nonpasses are unchanged upstream `T12-COVERAGE-STRICT-001`
(open coverage) and `CSWC-FULL-CSWC0628-MASS-TIME` (routine 180 s timeout).
The native Enceladus E3 test completed all 10,000 steps and 13 assertions.
Missing ignored native inputs, the shape manifest and reduced-CME task document
were recovered byte-for-byte from the existing local checkout, with explicit
hashes/source-tracking provenance; both sides use identical authentic fixtures.
Initial incomplete-fixture/SHARE setup attempts remain retained. One setup
attempt overlapped SHARE installation with native owners and was correctly
rejected; the completed same-criteria native retry has 296 PASS/35 SKIPPED.
Checkout-specific Enceladus exporter path comments are regenerated only after
verifying both complete non-path payloads remain byte-identical.

Two additional infrastructure defects were fixed without changing a physical
oracle. I12's qualification-only shared-initialization query no longer prevents
CUSTOM regression production wiring from executing; its eight-rank rerun
retains all artifacts and correctly SKIPs qualification. Normal builders and
the global Moon adapter explicitly set -amps-test=off, because application
selection alone retains nightly TESTMODE (default iteration/debug/sorting/
load-balancing policy). The initial inherited-mode restore remains preserved;
the final normal-mode Na/Ar supporting campaign and one-/eight-rank sampling/
Na event-budget rechecks PASS. The focused global Moon/harness retry has
27 PASS/22 SKIPPED, zero FAIL/ERROR. Twelve campaign protocol checks PASS.

The canonical nightly remains explicitly NO-SPICE/TESTMODE=on/four ranks.
Three completed outputs have identical SHA-256
`ed875cefc604e4ae3580fb89b8c528146929b6b9be099b41b61697faee7f0080`.
The user separately approved those bytes after regression preservation; the
previous 4bb811b3… reference and its provenance/acceptance remain archived.
Eight comparator/provenance protocol checks and the actual strict comparison
PASS. Root test_Moon.diff is zero bytes and the traditional completion line is
present. M0 reconciliation records nine original passing linked tests plus the
approved I04 PASS; original ERROR/FAIL attempts are never overwritten.

The final eight-rank supporting campaign has U08-U13 six PASS, formal I14-I20
seven SKIPPED, supporting six PASS/one D02 SKIPPED, zero FAIL/ERROR. I19's
42 rows close counts exactly; maximum relative Na/Ar imbalance is about
1.0664104644330185e-15. Normal Na budget tests have zero count residual and zero
numerical removals; eight-rank survival is 2062/4096 sunlit and 4096/4096 shadow
under the unchanged confidence criterion. Sampling's nine windows have exact
zero collecting/completed/callback errors on one/eight ranks.

Data raw bytes, frozen observation selections and both roadmaps are unchanged.
D01 LDEM_4.IMG SHA-256 is
`972d66d7f330abca62abc21aae270661516c00dd1a881ec102b94df9563656b5`;
its label is `b5a197b4d417e3b5e0e84f8fe15f9a1d2203b8eb5d7213236c0ffbc86e940cc6`.
D01 is READY, scoped D03 remains READY, D02 PARTIAL. Per-package and SPICE
hashes are retained in manifests/data audit. Source roadmap hash remains
2761b447… and data-root roadmap remains d10eca1c…; reconciliation still applies.
M1/I13 remains the formal frontier. No observation was scored, I13 baseline
frozen, or skipped/preflight result promoted to validation.

### 5.10 D02 research, I19/I20 proposal and verification follow-up — 2026-11-12

The user authorized the missing-metadata literature/archive research, a concrete
I19/I20 proposal and remaining approved verification. Source revision is
`06eafde4802f0344694612c8f3c21f2b28e1a09f`; production C++ physics is unchanged
in this work package. The analysis-before-change contract is
`test/M2_NEXT_VERIFICATION_CONTRACT.md`. No commit, push, staging, observation
edit or acceptance change was made. Existing user changes and the active
checkout's build configuration are preserved.

Durable evidence: `test_output/srcMoon/m2-next-20261113T042634Z/`.
Native builds/runs used the isolated checkout
`/tmp/amps-m2-next-20261113T042634Z-ihhn0q48/verification/`, never the active
checkout's generated tree. `completion_report.json` records exact commands,
configuration/executable/source/data hashes, counts, metrics and limitations;
the original configuration/patch and source snapshots remain retained.

**Research result.** `test/D02_ILLUMINATION_RESEARCH.md` records native labels,
NASA PGDA, the LOLA archive README, the 2011 publication's NASA record and the
author's LEAG abstract. Four primary-document snapshots retain URLs/access
times/hashes. General hourly sampling and older six-hour calculations do not
establish the exact 201608 simulation epochs/cadence/geometry. Native label
dates are retained without assigning an unproven role as averaging bounds.
AVGVISIB is occurrence of any visible solar-disc portion, distinct from mean
visible area or incident flux. Exact-product metadata and an independent
timestamped I15 reference remain missing. D02 stays PARTIAL / NOT VALIDATED.

The spatial processor previously placed unsupported timing/azimuth defaults
in QA. Its corrected QA v2 reports spatial processing READY separately from
overall PARTIAL and explicit unresolved operator fields. All **366** CSV rows
remain byte-identical before/after, SHA-256
`7c26d67be6a6b98e27c13f05a560607bea76ad03210d2820fb97caa39ed6af14`.
No canonical raw or processed data was rewritten. The new read-only D02 audit
and ten parser/corruption/status checks are registered automatically through
the existing global M2 data-audit entry, with no new permanent U/I ID.

**Concrete proposal.** `test/M2_I19_I20_PROTOCOL_PROPOSAL.md` and adjacent JSON
are **PROPOSED_NOT_APPROVED** and are not consumed by runners. They specify
synchronized physical-time accounting, common-global Na/Ar particle steps,
independent 3x3 dt/weight refinement, 16 fixed seeds, simultaneous paired-seed
equivalence intervals, reversible/censored residence accounting and frozen
polar-site selection. Proposed numerical margins are 1% of a declared
available-atom scale, one percentage point of retention/release fraction,
1% of the reporting duration for restricted mean residence, and 0.1% for
numerical-domain losses. Existing exact count/roundoff and 0.5% conservation
gates remain unchanged. These numerical budgets and the long common-duration
protocol require review before contract changes or formal execution. M1's
500-step samples/5000-step cap and all formal predecessors are unchanged.

| Executed scope | PASS | FAIL | SKIPPED | ERROR |
|---|---:|---:|---:|---:|
| Clean configuration-separated M0, I01-I10 | 10 | 0 | 0 | 0 |
| Neutral-SPICE U01-U26 invocation | 14 | 0 | 12 | 0 |
| Dedicated-ion U06 prerequisite | 1 | 0 | 0 | 0 |
| One-/eight-rank normal Na budget and sampling fixtures | 4 | 0 | 0 | 0 |
| M2 dedicated Na/Ar local U08-U13 invocation | 6 | 0 | 0 | 0 |
| M2 formal I14-I20 | 0 | 0 | 7 | 0 |
| M2 supporting scope: three linked checks plus four data/integrity checks | 6 | 0 | 1 | 0 |
| D02 metadata protocol checks, supporting infrastructure only | 10 | 0 | 0 | 0 |
| Global/I11/I17 protocol checks, supporting infrastructure only | 21 | 0 | 0 | 0 |

These rows are separate invocations/scopes, not unique test-ID counts: U06
passes in the proper ion build; U08-U13 are repeated in M2. The architecture
audit also PASSes with zero physical test records. Expected injected failures
inside protocol unit tests are not real campaign FAIL/ERROR outcomes.

Selected metrics:

- I04 strict approved-reference comparison produces a zero-byte diff against
  `ed875cefc604e4ae3580fb89b8c528146929b6b9be099b41b61697faee7f0080`.
  I10 finest scaled error is 3.46325e-9 and observed orders are approximately
  1.999917/1.999959/1.999980. No nightly/I13 reference was changed.
- One-/eight-rank Na gas ledgers have exact zero count residual and zero
  numerical removals. Sunlit/shadow survival is 2050/4096 and 4096/4096 on
  one rank, 2062/4096 and 4096/4096 on eight ranks, passing the unchanged
  analytical statistical criterion. Nine sampling windows have zero buffer/
  callback errors at both layouts.
- Supporting I14 has 513/513 radial hits, zero inverted/degenerate faces,
  maximum DEM control error 4.66e-10 m and normal mover collisions on all
  seven eligible facet-owner ranks. I18 passes the frozen residence-rate,
  sticking and velocity-distribution criteria, including 32768 Ar samples.
- Supporting I19 has 42 ledger rows, exact model-particle count closure and
  maximum Na/Ar relative imbalance 1.0664104644330185e-15. It remains a
  20-step scoped fixture, not full temporal/weight convergence or I19 PASS.
  Its sampled peak process-tree RSS is 22,297,004 KiB (about 21.26 GiB),
  elapsed resource trace 23.73 s. This does not estimate a month-long campaign.
- D01 and scoped D03 input/integrity checks PASS. D02 native integrity PASSes
  while its overall operator readiness is SKIPPED. No June/December residuals
  are scored and no observation is synthesized.

Reproduction commands, **inside an isolated checkout with installed SHARE and
the same frozen data/dependencies**, using new absolute output directories:

```sh
./Config.pl -compiler-option=-g
make -j4 share_lib
python3 -B srcMoon/test/run_integration_tests.py --output-dir /NEW/m0
python3 -B srcMoon/test/integration-tests/I11/run_budget_check.py \
  --amps ./amps --input srcMoon/examples/m1_memory200.in --mpi-ranks 8 \
  --output-dir /NEW/budget-8
python3 -B srcMoon/test/integration-tests/I11/run_sampling_check.py \
  --amps ./amps --input srcMoon/examples/m1_memory200.in --mpi-ranks 8 \
  --output-dir /NEW/sampling-8
python3 -B srcMoon/test/test_runner.py --phase m2-development --mpi-ranks 8 \
  --output-dir /NEW/m2-development
```

The two fixture commands also executed with `--mpi-ranks 1`. Exact actual
argv/log/exit receipts are `verification-commands.json`,
`supporting-prerequisite-commands.json`, `m2-command.json` and
`protocol-regression/result.json`. Native helpers remove `build/` before each
configuration, use explicit TESTMODE=off for normal builds and TESTMODE=on /
NO-SPICE for the four-rank nightly, and build with `make -j16 'OPT3=-O3 -g'`.
Compiler is g++ 9.4.0, MPI Open MPI 4.0.3; M2 uses eight ranks, OpenMP=1 and
seed 731993. The dedicated Na/Ar binary has chemistry OFF/individual weights
ON and SHA-256 `dda735ddacc4c28794b953e22dc136da0e239b21dd5dab3b25f8d40a4cf27fc1`.
Normal Na-only Moon/SPICE configuration is restored in the private checkout.
No model/build remains active.

D01 native IMG/LBL hashes remain `972d66d7...656b5` / `b5a197b4...940cc6`
(full values in section 5.9 and manifests). The D02 audit verifies all eight
frozen IMG/LBL hashes; D03 and eight SPICE hashes are in the linked manifests.
Both roadmap hashes are unchanged from section 2. Source amendments govern
the approved supporting-development exception; common scientific gates and
native archive metadata remain controlling. Required roadmap changes are
limited to reviewing/adopting a future numerical protocol consistently; none
was silently made by this work package.

Files changed are the D02 processor/new native audit/protocol tests, M2 package
audit, I15 missing-capability explanation, research/proposal/contract documents,
gap registry, question record, STATUS and relevant README/RUN files. Production
surface/mover/source/thermal/chemistry files were inspected through existing
contracts and executed without modification. Formal M1/M2 remain incomplete.

### 5.11 Detailed numerical-protocol recommendations — 2026-11-12

The user requested detailed suggestions for item 1, the I19/I20 numerical
protocol decision. The existing proposal Markdown/JSON is expanded to version
2, explicitly PROPOSED_NOT_APPROVED and unconsumed by runners. This is a
documentation-only work package: no production code, active acceptance JSON,
roadmap, frozen observation selection or build configuration changed. No
native build/model run was launched, and no commit/push/staging occurred.

The review inspects MoonIntegrationTests, MoonM2Budget, main_lib, finite-stock
injection, surface kinetics, the production RNG, I19/I20 contracts, existing
receipts and primary NIST statistical/survival documentation. It makes these
previously underspecified choices explicit:

- Common Na/Ar particle time after the final 0.1 Ar timestep copy; integer
  alignment of all physical checkpoints and independently pinned atom weights.
- Real seed dispatch rather than changing a manifest while the fixture still
  hard-codes 731993+rank; a documented production RNG/stream qualification.
- Separate source-fed zero-stock and controlled surface-cohort cases; first
  release, re-adsorption/current occupancy and censored restricted residence.
- Frozen 18-pair refinement comparisons, a conservative maximum confidence
  family of 8640 and the limitations of sixteen-seed Student-t intervals.
- Unchanged exact/roundoff/0.5% conservation requirements, explicit review of
  external-supply versus gross recycling normalization, and proposed new
  convergence/loss budgets that remain unapproved.
- Bounded capacity/precision planning before deciding on the full 288-run
  physical design. At the retained 0.6855865301040126 s bound, checkpoint
  alignment yields 3,719,800 coarse and 14,879,200 fine steps per case; the full
  design totals about 2.50 billion timesteps. This is planning arithmetic,
  not a measured runtime estimate or an assertion that the campaign fits
  routine regression resources.

Existing source revision remains `06eafde4802f0344694612c8f3c21f2b28e1a09f`.
Both roadmap hashes remain those in section 2; reconciliation still controls.
The active I19/I20 physics/numerical gates and D02/M1/I13 prerequisites remain
unresolved as before. Previous execution counts in section 5.10 are unchanged;
no fresh physical PASS is claimed. The proposal/JSON consistency and contract-
preservation checks, before/after proposal snapshots, source bibliography and
completion receipt are under
`test_output/srcMoon/m2-protocol-review-20261113T052112Z/`.

### 5.12 User-directed one-hour development / cluster qualification tiers — 2026-11-12

The user requires development tests to execute within an hour and permits long
qualification on a larger cluster for a few hundred CPU hours. The source
Phase-2 roadmap now records a single 3600-s end-to-end development requirement,
including required builds, QA, selected tests, figures and cleanup. The v3
proposal targets 3000 s work plus 600 s reserve and explicitly separates
short development evidence from deferred cluster qualification. This replaces
the v2 planning schedule that excluded builds/grace from its 50-minute slots.

The proposed short design uses two physical cases, a 3x3 dt/weight grid, two
development seeds and 24/48/96 steps at matched physical horizons (36 short
runs). These are proposed workload settings pending implementation/timing
qualification; they do not establish month-long convergence or sixteen-seed
confidence. Long cases are excluded from development defaults. The candidate
288-run design remains deferred, with a proposed 300 CPU-core-hour aggregate
planning target and mandatory measured feasibility estimate. Neither an
incomplete matrix nor a short PASS closes formal I19/I20.

The current run_m2_development.py gives I14/I18/I19 separate child timeouts;
it does not enforce the required whole-command hour. Aggregate deadline,
build/child budget propagation, safe shutdown and opt-in cluster accounting
remain implementation work. No runner or production physics changed in this
documentation amendment, no new model/build was launched, and no acceptance
threshold, source coefficient, frozen observation or I13 baseline changed.
M1's 500/5000 rule remains intact.

Current source-roadmap SHA-256: `bde5f001bdcbe072d7a8452d3f2afec93856f50cd1146560da3cd3027a96d2dc`.
Pre-amendment source hash: `2761b447b46539edb21395550cd1b5d9ff6add4839b0db5a7e6d734cf2a4f02a`.
External roadmap remains `d10eca1c539f9cf82940236d9e1aa36af665d43e82ba475925950d582fbf49c4`.
The source execution amendment controls runtime policy; common scientific
gates and native data rules still control qualification. Historical results
and hashes are retained. Scope/consistency/preservation checks and snapshots
are under `test_output/srcMoon/m2-execution-budget-20261113T054530Z/`.

## 6. Known blockers and risks

The current-core I04 review is resolved: the user approved the exact ed875cef…
replacement after all 919 upstream PASS identities were preserved. The former
reference/provenance/acceptance and rejected comparisons remain archived.
Current M0 has ten PASS. The pre-existing global strict-coverage FAIL and CME
timeout ERROR remain visible and are not repaired by this integration.

1. The corrected Na dispatch changes the physical configuration relative to
   previous duration studies. Its I04 replacement was explicitly reviewed and
   the strict checker and complete clean M0 requalification PASS.
   The retained 2026-10-31 I11 campaign FAILs the unchanged statistical gates.
   It and the former duration pilot have stale-window sampling contamination.
   The corrected sampler is now linked and its deterministic storage checks
   PASS. The corrected replacement campaign also FAILs all four stages and
   directly confirms startup in actual particle inventory.
   The fixed step cap still supplies unequal physical durations across grids.
2. I13 requires investigator review and freezing of a passing I11 candidate;
   the runner must never install its own baseline.
3. Formal M2 qualification cannot open until M1 passes in the required campaign
   context. M2 implementation/local verification is explicitly authorized now.
4. D03 selected-input package readiness is now established, but
   linked predecessor/operator/threshold/hold-out requirements must be checked
   before any new score. Directory presence is insufficient.
5. I15 lacks the qualified D02 temporal/extended-disc comparison. I19 now has
   supporting finite Na/Ar gas+surface accounting, but its complete formal
   campaign and predecessor gates remain open; these still block I20.
   Research additionally identifies the missing independent time-resolved
   I15 benchmark; average maps cannot supply individual epoch labels. The
   concrete I19/I20 numerical protocol remains a proposal awaiting review.
6. I21-I46 require real production capabilities, independent references and
   data packages; generic placeholder PASS runners are prohibited.
7. The checkout is dirty and contains user work. `Makefile` is generated and
   must not be committed. Do not use blanket restore/clean/reset operations.

## 7. Next work in roadmap order

Current-core integration, approved I04 comparison and the normal-mode restored
build are complete. Preserve the source/build identities and complete regression
evidence before continuing the unchanged M1/M2 qualification work below.

1. Review the completed corrected-chemistry PARTIAL pilot. Full M0
   requalification PASSes, and the actual-event audit closes the neutral-Na gas
   budget and identifies the former missing Na photolytic sink. Previous
   bounded I11 runs used different chemistry and cannot qualify the corrected
   case. Independence was their first final gate failure at every stage;
   transient inventory and unequal physical durations remain measured
   concerns. The corrected pilot supplies exact per-observable gate failures,
   but no accepted duration; preserve those facts when designing the next
   reviewed numerical protocol.
2. Review/freeze a justified spin-up and common physical-duration contract, or
   independently qualify an equilibrium initial state, before changing the
   5000-total-step policy. A corrected non-qualifying duration pilot is already
   available; its result cannot automatically change the qualification
   contract. Preserve I12 PASS, and defer I13 reference review until I11 PASS.
3. Continue approved M2 implementation/supporting verification, retaining
   separate results and data QA. Only after M1 passes, execute formal I14-I20
   in order with package preflights and comparison CSV/PNG/EPS.
4. Complete M2 budget/operator gaps before creating I21-I24.
5. Implement later I-runners only when the named production path and oracle or
   authoritative data exist.

### 7.1 Exact M1 completion checklist

M1 is complete only when all of the following evidence exists in a retained
qualification campaign; a QUICK or bounded-regression wiring result is not a
substitute:

- **I11 PASS:** execute the approved 480/240/120/60 km sequence with identical
  physics, epoch, seed, MPI layout, extents and sampling grids.  At every
  level require zero external-state drift, equilibrium within 2% for three
  consecutive blocks, at least five independent science blocks, 95% CI
  half-widths no larger than 1.67% integrated and 3.33% pointwise, and then
  successive-grid changes below 5% integrated and 10% pointwise.  Sampling is
  500 iterations and each level is capped at 5000 iterations; reaching the cap
  without satisfying the unchanged gates is FAIL, not completion.
- **I12 PASS:** run the production angular/SZA/terminator/polar selector and
  realized-mesh diagnostic with the approved `m1_memory200.in` memory
  profile (or explicit legacy FULL). Requested values must
  agree to `1e-14` relative error and realized
  boundary cells must satisfy the frozen 25%/documented-dyadic-AMR rule.  The
  earlier QUICK wiring result remains `SKIPPED / NOT QUALIFIED`.
- **I13 reviewed reference and PASS:** use the exact finest input from the
  passing I11 campaign, review the candidate science CSV and deterministic
  configuration JSON, freeze both with hashes, and rerun.  Sphere-area error
  must be at most `1e-10`, configuration identity must be exact, and every
  approved observable must remain within its predeclared I11 numerical plus
  statistical tolerance.  The runner is intentionally unable to approve or
  rewrite its own reference.
- Preserve the required logs, input snapshots, executable/source/SPICE hashes,
  MPI/OpenMP/seed identity, resource histories, block metrics, comparison
  tables and manifests.  Required predecessor regressions (including I10 and
  U03-U06 where specified) must remain PASS.

No additional physics choice is currently needed merely to attempt I11/I12.
After I11 passes, investigator review/approval of the proposed I13 reference
is a required scientific decision before the reference can be frozen.

### 7.2 M2 work remaining in addition to M1

Closing M1 opens formal M2 qualification; approved supporting development can
already proceed. M2 additionally requires
all I14-I20 linked gates to execute in order and PASS:

- **I14:** run the production LOLA boundary against READY D01.  D01 is READY
  in the current package; the linked comparison still needs a same-campaign
  I13 PASS and retained CSV/PNG/EPS evidence.
- **I15:** complete and verify the production extended-solar-disc temporal
  illumination operator. D02 spatial rows exist, but actual integration epochs
  and a matching temporal comparison remain unresolved; the package stays
  PARTIAL. Then meet the 99% visibility and PSR-resolution gates.
- **I16:** preserve the established scoped D03 package READY state and verify
  current hashes/selection before each run, then run
  the flagged-data-only Diviner GCP comparison and meet the frozen 15 K
  non-PSR RMSE gate without gap filling or assumed uncertainty.
- **I17:** complete qualified spin-up and time-step convergence, review and
  freeze the June training threshold before any December inspection, then run
  the untouched December hold-out.  The current experimental June plotting
  path and threshold candidate are not a PASS.  D03 package qualification is
   established only for selected inputs and must remain valid for each run.
- **I18:** build the dedicated Na+Ar executable and pass the linked stochastic
  sticking, residence, re-emission, desorption-rate and Ar-source gates using
  the already frozen literature parameter contracts.
- **I19:** extend the now-passing supporting Na/Ar finite-stock ledger to the
  complete formal campaign and required source/mover paths; verify temporal/
  weight convergence and preserve exact distinction between physical and
  numerical loss. Demonstrate deterministic closure to machine precision
  and imbalance below 0.5% of the dominant source. Chemistry remains OFF in
  this approved M2 case; later conversion networks need separate coverage.
- **I20:** after I15/I17/I18/I19 pass, qualify the combined D02/D03 polar
  comparison and demonstrate kinetic cold-trap retention/release rather than
  equating a geometric PSR flag with permanent trapping.

Every evaluated I14-I20 comparison must preserve the exact scored rows in CSV
and render those rows to PNG and EPS.  Missing capability or authoritative
metadata remains SKIPPED; an evaluated threshold miss remains FAIL.

## 8. Change log

### 2026-11-11 — sampling repair and concurrent M1 execution

- Repaired common production cell-sampling buffer reset independently of
  species visualization output and added deterministic linked storage checks.
- Verified Na/Na+Ar storage layouts on one/eight ranks with exact-zero errors;
  added instantaneous weighted particle diagnostics and shared mesh receipts.
- Traced the nightly mismatch to a shared-data reference update, restored the
  exact already-reviewed Git reference locally, and obtained a zero-byte diff.
- Added guarded concurrent mesh execution, explicit pending plot artifacts,
  active-runner announcements and reproducible species-build source tracking.
- Executed I11-I13 in 37.14 min: 1 PASS, 1 FAIL, 1 SKIPPED, zero ERROR. All
  I11 levels completed ten samples but failed independence at the final
  evaluation. Preserved every final stage plot, manifest and resource record.
- Kept M1 incomplete, I13 unfrozen, M2 experimental and all thresholds intact.
  Recorded the remaining statistical-support audit as an I11 extension
  proposal without assigning a new test ID.

### 2026-10-30 — comprehensive campaign architecture

- Added `srcMoon/test/TESTING_PROCEDURE.md` with the complete component and
  milestone test sequence, reference hierarchy, visualization rules, data
  preflight, artifacts and commands.
- Added `srcMoon/test/test_runner.py` as the single top-level audit/delegation/
  aggregation entry point for U01-U26/I01-I46.
- Kept I21-I46 explicit as planned `SKIPPED` gaps rather than fabricated tests.
- Added fail-closed architecture checks for individual runner/README/reference
  presence, mandatory acceptance fields and I14-I20 comparison artifacts.
- Added `test_runner_self_test.py`; its four infrastructure checks pass and do
  not claim any U/I physics result.
- Corrected `run_integration_tests.py` to pass explicit U01-U26 selections so
  its local-prerequisite step cannot accidentally launch I11-I20.
- Recorded the formal frontier as M1 and M2 as experimental/not qualified.
- Recorded exact M1 and M2 completion checklists, including the independent
  D02/D03, I15, I17 and I19 blockers that remain after M1.

### 2026-10-31 — bounded I11 qualification execution

- Performed a clean Moon/SPICE `-O3 -g` build and ran I11 with eight MPI ranks.
- Preserved the sandboxed OpenMPI pre-launch `ERROR` separately, then executed
  the unrestricted retry using the same hashed executable.
- Completed all four 5000-step levels and retained 153 GB of manifests, raw
  production output, metrics and diagnostic-only CSV/PNG/EPS figures.
- Classified I11 as FAIL at every level because statistical qualification was
  not reached; no accepted cross-grid comparison was produced.
- Identified monotonically increasing inventory and unequal physical duration
  across resolution-dependent timesteps as the first blocking mechanism.

### Earlier frozen evidence

- `aa160e0`: M0 production-linked verification completed.
- `92386a3`: configurable M1 mesh/observation framework.
- `fb2008e`: experimental M2 surface/thermal validation framework.
- `b337b8b`: bounded I11 diagnostics and LOLA surface visualization.
- `22ad748`: illumination/thermal-state demonstration runner.

Future entries must state files changed, commands/tests, exact status counts,
data/hashes, scientific decisions, unresolved blockers and roadmap impact.

### 5.13 Approved M2 execution/numerical implementation — 2026-11-13

Implementation and execution evidence are being assembled under
`test_output/srcMoon/m2-implementation-20261113T060631Z/`. This checkpoint records completed implementation and supporting verification.
It is not a formal M1/M2 qualification report.
The native verification checkout is
`/tmp/amps-m2-implementation-20261113T060631Z-t0joj5cp/verification`.
The campaign ran from the global entry with eight MPI ranks and a complete
clean-build deadline. Root configuration hashes are preserved in
`initial_state.json`; initial smoke source and receipts are versioned.

Initial protocol checks: six PASS (expanded final suite: thirteen PASS); global-runner infrastructure checks: five
PASS; original I19 budget infrastructure: five PASS. These checks exercise
protocol/control/process handling and do not independently verify lunar
physics. Final native results, clock/weight/tag metrics and full wall timing
are recorded below and in completion_report.json. Formal I14-I20 remain gated by I13
and all original environment/data/operator conditions.

#### 5.13.1 Preserved attempts and corrections

- `development-clean/`: clean M0 has ten PASS. Dedicated M2 configuration
  ERROR because its explicit RNG option was placed in `#main`; inspection of
  ampsConfig.pl identified `#General` as the native owner. Corrected the
  authoritative derived input generator; canonical input is unchanged. The
  normal restore completed; whole attempt 753.387 s.
- `development-clean-2/`: six local PASS, seven formal SKIPPED; supporting
  43 PASS, one FAIL, two SKIPPED, zero ERROR; 1905.334 s. All 36 numerical
  cases PASS, source-normalization READY and canonical restore completed.
  The failed original I18 adapter inherited the numerical template's base
  2000003, contradicting its frozen base 731993 and recorded manifest. The
  original FAIL/manifest are preserved with a separate provenance correction.
- `frozen-I18/`: one targeted eight-rank rerun at the already-frozen 731993
  PASSes all unchanged gates. Na angular KS 0.00376990 <0.01039501; Na energy
  KS 0.00667384 <0.01039501; Ar energy/cosine KS 0.00387311/0.00539714
  <0.00751301. No seed was selected from a pool, no threshold was widened,
  and the numerical sweep's 2000003/2001012 seeds remain unchanged.
- `core-check/`: five selected suites, 123 actual cases, all PASS in 280.316 s,
  peak aggregate process-tree RSS 4,748,115,968 bytes. The exact selected
  case audit confirms all 123 were previous baseline PASS and are preserved.
  This is not a fresh execution of the entire historical 1,193-case registry.
- `development-clean-3/`: final corrected clean command completed with zero
  FAIL/ERROR. It owns
  restoration from the first native configuration change, avoids a discarded
  intermediate M0 restore build, pins the original fixture seed, and includes
  the configured LOLA image/label hashes in linked manifests.

All attempts, original reports, source snapshots, inputs, executables, streams,
resources and comparison CSV/PNG/EPS are retained. Root configuration remains
unchanged. Formal M1/I13, D02 temporal/operator evidence, full stochastic/long
qualification and independent PSR/thermal comparisons remain unresolved.

#### 5.13.2 Final measured result and limits

The final command (from the isolated verification checkout) was:

```sh
python3 -B srcMoon/test/test_runner.py --phase m2-development \
  --mpi-ranks 8 --stream --m2-input srcMoon/examples/m2_development.in \
  --output-dir /data/vtenishe/Moon-test/AMPS/test_output/srcMoon/m2-implementation-20261113T060631Z/development-clean-3
```

Whole global wall time was **1849.678 s**, including clean SPICE/no-SPICE/ion/
nightly M0 variants, the dedicated M2 build, QA, all selected tests, comparisons
and final ordinary Na-only Moon/SPICE restoration. This meets the one-hour
requirement. The monotonic supervisor measured 1849.558 s, peak descendant RSS
**45,977,960,448 bytes** (45.98 GB), and approximately 11.59 GB of campaign
artifacts, below the 200/50 GB limits. Sampled local CPU is a 5.0378 core-hour
lower bound, not scheduler charging. The global result is SKIPPED/exit 0
because formal gates remain unavailable; it must not be described as M2 PASS.

- Local U08-U13: **6 PASS**. Formal I14-I20: **7 SKIPPED** without I13.
- Supporting namespace: **44 PASS, 2 SKIPPED, 0 FAIL, 0 ERROR**. The skips
  are D02 package readiness and unselected long cluster qualification.
- M0: **10 PASS, 0 FAIL/ERROR/SKIPPED**, approved reference unchanged;
  `m0/I04/test_Moon.diff` is **0 bytes**.
- Fixed short numerical design: **36 PASS**, two cases, nine dt/weight cells
  and two predetermined seeds on eight ranks. All use H=16.454076722496303 s
  at 24/48/96 steps with independent W, W/2 and W/4 controls. Exact gas count and
  independently evaluated signed atom/kg identities pass. Maximum relative
  ledger imbalance is **5.16450e-10**, within the frozen arithmetic bound;
  maximum relative passive-tag residual is **5.41589e-16**. These are short
  accounting/conditional kinetic results, not full statistical convergence.
- Independent warming release: finest restricted-residence error about
  **1.24e-4 s** for Na and Ar, versus the approved **0.6 s** allocation.
- Source-only normalization: READY; at H=2550240 s the native integral gives
  **4.3099056e28 Na atoms** and **9.690912e26 Ar atoms**. The independent
  constant-Ar Q*H gate passes. This does not simulate one month of transport.
- `explicit-controls-cohort/`: additional eight-rank production-path PASS
  for explicit coarse dt=0.5 s, dt divisor 2, weight divisor 4, actual dt=0.25 s,
  weights 2.5e20/2.5e19 Na/Ar atoms, H=12 s, 48 steps, and distinct tags
  sharing element 3000. Full CSV and unchanged arithmetic/hazard gates pass.
- Infrastructure/oracle checks: **13 protocol**, **5 global**, **5 original
  I19 budget**, and **10 D02 metadata** checks PASS. The final campaign itself
  ran the then-current ten protocol checks; three subsequent targeted checks
  cover explicit input preservation and cluster false-green/corruption cases.
- Cluster default planning and missing-proof execution both correctly retain
  **I19/I20 SKIPPED**, no native long run; the fixed plan has 288 cases.

Every completed numerical stage retains input/manifest/receipt/streams,
cumulative and window CSV, `comparison.csv/.png/.eps`, hazard/cohort data
where enabled, and resource logs. Manifests now include the actual configured
LOLA raw image/label hashes for every terrain-consuming ID. Earlier diagnostic
receipts have supplemental `external_data_manifest.json`; originals are
unchanged. All failed attempts and the Python-3.8 infrastructure-fixture error
are retained; no failed evidence or stochastic seed was discarded.

The native/science build was GCC/OpenMPI, C++17, `-O3 -g`, OpenMP=1. Dedicated
M2 explicitly selects Mersenne Twister, neutral Na/Ar, chemistry OFF and
individual weights ON. Normal M0/M1 RNG/input/reference settings are preserved.
Current source revision is 06eafde4802f0344694612c8f3c21f2b28e1a09f plus the
uncommitted changes captured by source/build hashes. The final source-roadmap
hash is f1a2f9265f97f36da9ef9f6a17b3c83783d9f4b51882e6d81081ef7d567757cd;
external roadmap remains d10eca1c539f9cf82940236d9e1aa36af665d43e82ba475925950d582fbf49c4.

**Remaining:** M1/I13 closure; D02 actual temporal/operator evidence; qualified
actual-epoch thermal state and independent PSR/thermal comparisons; frozen
native-site/collocation/statistical support; measured feasibility and explicit
long cluster qualification. No monthly precision, observational validation,
M2 completion, or M3 readiness follows from this short supporting package.

### 5.14 Remaining M2 code completion — 2026-11-16

The analysis-before-change inventory/test contract is `test/M2_COMPLETION_CONTRACT.md`.
Source HEAD is `17d49e00d1d2ccf33b06d11a3e6e5ffef34ac810` plus retained dirty changes.
The new evidence root is `test_output/srcMoon/m2-completion-20261116T195511Z/`;
native builds use `/tmp/amps-m2-completion-20261116T195511Z/verification`.
Root generated Makefile/configuration and pre-existing changes are preserved.

Implemented: production finite-disc correlated masks/incidence/albedo; shared
thermal solver flux entry; complete geometry-only I15 runtime/export with
reviewed independent scoring and exact-product/epoch/resolution gates; I17
selected-only release with compiled native authorization guard, independent
explicit December spin-up and reviewed metric groups; normal I19 dispatch
using real native-case receipts; I20 production temperature/depth exporter and
strict native-site/time/quantity comparison. Retained Na/Ar binaries now carry
actual generated ABI/RNG header snapshots so later Na-only builds cannot
invalidate authentic retained evidence, while probe reuse remains live-ABI strict.

The user confirms no additional D02 metadata exist. A geometry-only June
audit counted 1,580 sites: 1,126 have one bin, 449 have two, five have three.
These are sampling counts, not residuals; no local eclipse phase curve is
invented from them. Formal I17 also needs actual reviewed phase/amplitude
support and approximation-scope evidence, explicit December UTC bounds and
untouched hold-out execution. I20 additionally needs reviewed native annual
time/temperature/depth/collocation evidence. No real observational residual
was inspected during implementation, and no frozen threshold was auto-created.

Combined clean-build/global test results are pending; local protocol reports
are retained separately and are not linked model qualification. Existing
M1/I13 and all formal I14-I20 gates remain intact. Current source roadmap
SHA-256: `eba1571968afff4f948505b77578ee782c4fb54d02ac786ab9ecaa7ea00438d2`; external: `d10eca1c539f9cf82940236d9e1aa36af665d43e82ba475925950d582fbf49c4`.

Production review also identified the PSD source's stale radial-cosine
illumination/cache path. An optional NULL-default core geometry hook now
uses the same production terrain/Earth projected direct-photon integral,
actual point-to-Sun distance and explicit source-BEGIN TDB epoch. Its cache
is prepared serially before callbacks; the source CDF, total selector and
finite-inventory hazard use the same element kernel. Thermal albedo is not
part of PSD photon incidence. Disabling M2 removes only its own provider
before destroying the cache. U09 extends its analytical/default/lifecycle
checks; the native I15 support fixture now has nine queries (the original
eight plus one independently converted UTC startup query proving the actual
Moon provider's clock). This does not establish observational source strength.

Execution checkpoint: the isolated `development-clean-5` campaign is running
with eight MPI ranks and the unchanged 3600-second whole-command budget.
Its nine new infrastructure groups have completed **149 tests, all PASS**;
native M0/M2 build and production results are pending. Earlier attempts are
preserved as ERROR evidence: the shell's default Python lacks scientific
dependencies; the initial isolated clone omitted ignored installed machine
configuration; a process-supervisor orphan race was reproduced and corrected.
The unmodified orphan test passed ten host repetitions and the complete
protocol suite passed 14 tests. Machine configuration was copied from the
real checkout, not invented. `/usr/bin/python3` supplies the installed
NumPy/SciPy/Matplotlib/mpmath/psutil dependencies. The clean native build has
started successfully after copying `.general.conf` as well as `Makefile*`.
A separate, isolated full CORE regression retry is running with host MPI
permissions; the restricted-sandbox attempt is retained separately. No code
was committed/pushed and the active checkout was not configured or rebuilt.

The frozen `development-clean-5` attempt has now ended after **1467.38 s**.
Its M0 baseline is **10 PASS, zero FAIL/ERROR/SKIPPED**, including the reviewed
nightly reference and zero-byte `test_Moon.diff`; its U-series is **14 PASS,
12 SKIPPED, zero FAIL/ERROR** and its dedicated M2 U08-U13 are **6 PASS**.
Supporting linked I14/I18/I19 all PASS. Fourteen native input controls PASS,
but the negative PSD fatal-path adapter crashes after printing the intended
rejection. Actual host GDB identified `fprintf(NULL)` in the diagnostic writer:
the fixture lacked its configured `PT/amps_error.log` parent directory.
The fix is isolated diagnostic-directory/MPI setup in the thin test adapter,
not changing the shared production fatal handler or accepting a SIGSEGV.
Earlier results and GDB streams remain intact; the stopped attempt is ERROR,
not a successful combined campaign. The private ordinary build restored
successfully; active-root Makefile and `.amps.conf` hashes remain unchanged.

Independent review found three further implementation gaps: stale prescribed
source temperatures, checkpoint-dependent site retention reporting, and a
joint-Earth visibility quantity labelled geometric PSR. The user approved the
explicit `isothermal-grey` Diviner interpretation on November 16. Production
now has one prescribed boundary function and refreshed source-temperature
cache, a matching brightness-to-kinetic operator, a native per-species surface
residence/retention exporter independent of checkpoints, exact cohort/site
binding, and terrain-only PSR separate from the benchmark Earth policy.
Real map gaps remain unavailable and never inherit initialization temperatures.
Post-review infrastructure checks PASS (I15 19 tests, I20 57 tests, native
operator reporting 19 tests); **new native verification is pending** for this
later source state. These counts do not promote the frozen earlier run to
verification of changes made after its source snapshot.

The next frozen `development-clean-6` completed in **2772.32 s (46.2 min)**:
local results **6 PASS / 7 SKIPPED / 0 FAIL / 0 ERROR** and supporting checks
**58 PASS / 2 SKIPPED / 0 FAIL / 0 ERROR**. Supporting checks include the clean
M0 prerequisite (all ten linked gates PASS), eight-rank I14/I18/I19, the actual
I15/I20 exporters, corrected initialized fatal adapter, all 36 short numerical
cases, source quadrature and retained-build verification after restoration.
Its ten infrastructure groups contain 190 passing unit/control cases.
The two supporting skips are unavailable D02 metadata and unselected long
cluster qualification. Formal I14-I20 remain SKIPPED for their unchanged
predecessor/data/qualification gates. See `development-clean-6/` and `launch-6.json`.

The independent full CORE audit is complete: **843/843 prior PASS identities
preserved**, no missing/newly failing/status-changed case. Reconciled full
counts are **843 PASS / 1 FAIL / 240 SKIP / 1 ERROR**, not an all-green core
claim. The nonpasses remain strict coverage and the existing 180-second CME
routine ceiling. Three nested-Python dependency failures were explicitly
rerun with the installed system interpreter; all attempts remain retained.
Whole verification including setup/retries was **3286.84 s (54.8 min)**;
peak process-tree RSS **5,192,458,240 bytes**. Shared core source/configuration
identities match the tested snapshot and active checkout. Durable evidence:
`core-regressions/core_completion_report.json`,
`core-regressions/core-regression-comparison/`, and its preservation mapping.

The final formal-path audit found an unapplied-initial-state loophole and
bounded-I18 staging in long I19 cases. Current root source now has named
uniform/actual periodic initialization, actual full initial-column exports
before cohort seeding, exact applied-state/epoch/source/executable/physics
binding, preserved formal AMR/source/thermal inputs, and a separate cluster
normalization path. No qualifying state or interval has been invented or
auto-approved. I20 requires its own explicit initial-state qualification and
keeps annual thermal and cohort intervals distinct. Those later changes were
then frozen for the separate `development-clean-7` verification below; the
earlier 46.2-minute run remains preserved with its own source state. Contracts:
`test/M2_INITIAL_STATE_CONTRACT.md` and `test/FORMAL_STAGING_CONTRACT.md`.

Final frozen native verification, `development-clean-7`: **2242.36 s (37.4 min)**,
eight MPI ranks, OpenMP=1. The registry reports **6 PASS / 7 SKIPPED / 0 FAIL /
0 ERROR**; supporting results separately report **58 PASS / 2 SKIPPED / 0 FAIL /
0 ERROR**. All ten M0 gates PASS and `m0/I04/test_Moon.diff` is zero bytes.
U-series: 14 PASS / 12 SKIPPED, no FAIL/ERROR; dedicated U08-U13: six PASS.
Ten infrastructure groups contain 211 passing unit/control cases; five budget
and 26 protocol cases also pass. Eighteen native input/fatal controls pass.
All 36 short numerical cases and source normalization pass. The peak sampled
process-tree RSS is **45,979,766,784 bytes** (about 46 GB); sampled CPU is a
6.12 core-hour lower bound, not scheduler accounting. Supporting I19's maximum
relative atom imbalance is 9.5924e-16. Native I20 exports all 486,000 initial
boundary/layer rows; its millisecond periodic replay is wiring coverage,
explicitly not lunar spin-up or annual convergence qualification. The retained
Na/Ar executable remains authentic after ordinary Na-only restoration.

A final Python-only extension makes the global source normalizer explicitly
request the real I19 initial-state export. The focused eight-rank native
`initial-state-source-normalization/` check then completed in **22.64 s** after
build restoration, with the same retained executable and unchanged reviewed
bounded physical input. It validated all **486,000 rows / 6000 logical elements /
80 layers plus boundary**, 250 K declared-uniform state, actual TDB epoch,
ordered furnished SPICE pool and configured LOLA hashes. Its candidate carries
`qualified:false`; READY applies only to the deterministic model-input
normalization. All 27 updated protocol tests pass. This focused run exercises
the final Python extension; the frozen global campaign is not relabelled as
having executed that later extension.

The implementation/supporting work package is complete. **M2 is not qualified**:
formal I14-I20 remain SKIPPED, M1/I13 remains open, D02 lacks authoritative
timing/geometry conventions, I17 needs the explicit December interval and
reviewed frozen observational/operator evidence, and long I19/I20 need real
initial-state/statistical/cost evidence and execution. No hold-out residuals
were exposed, observations changed, thresholds weakened or approvals created.
No source/configuration build is active and no commit/push was made.

### 5.15 Initial M3 implementation and physics proposals — 2026-11-17

The user approved the recommendations to extend early development to M3,
use explicitly configured Maxwellian/Voronov electron rates, and research
concrete charge-exchange and molecular-product proposals before implementing
those unresolved choices. The approved exception is recorded in AGENTS.md
and the source roadmap; formal milestone order and thresholds are intact.

Current inventory: Na neutral-loss production dispatch is active but rejects
unqualified Na+ products. No D04 production reaction loader or lunar
charge-exchange channel exists. U16/U17 remain unimplemented at this initial
checkpoint. The dormant electron rate uses fixed cross sections/solar-wind
speed. Native D05 cfit.dat/cfit.txt/cfit.f and all 11 processed D04/D05 outputs
are present and match provenance hashes. D06 does not provide qualified
electron-density/temperature forcing. The shared electron executor also needs
separate product-transfer/cell-index review; a new rate kernel will not enable
that executor or claim particle-budget closure.

Evidence root: `test_output/srcMoon/m3-initial-20261117T071134Z/`.
The initial source/configuration snapshot and git status are preserved there.
Kernel and input integration are implemented. The registered U16 runner has
170 passing kernel/parser controls over 132 independent 70-digit rate points,
plus 17 passing compiled production runtime-parser/install/query controls.
Maximum coefficient/rate relative errors are 1.5637e-15/1.5135e-15, below the
unchanged 1e-12 criterion. Original native Fortran comparisons report their
single-precision/rounded-temperature discrepancy separately. These are local
rate results, not I23 ensemble/product/forcing qualification. Research
proposals and the new clean native/global campaign are being finalized.
Root Makefile/configuration are not being rebuilt or rewritten;
canonical validation data remain read-only. No commit/push is authorized.

Research proposals are now complete: `test/M3_CHARGE_EXCHANGE_PROPOSAL.md`
documents original He/Na total-capture references, supported energy ranges,
rate integration, physical uncertainties and missing joint product kinetics;
`test/M3_MOLECULAR_KINETICS_PROPOSAL.md` inventories all seven H2O/four
experimental-OH photon channels and the absent daughter/state/threshold data.
The molecular archive audit matches 16/16 hashes and makes no observations.
U17 remains SKIPPED. No proposed new reaction/kinematic closure is activated.

The first isolated native attempt ended after 44.45 seconds with a build
infrastructure ERROR: copied machine definitions retained the active tree's
absolute include paths. Its original coarse FAIL aggregate is preserved,
with a separate explicit ERROR audit. Only private machine paths were rebased;
root Makefile/.amps.conf hashes remain unchanged. The M3 supervisor now
inspects nested M0 ERROR records instead of treating every coarse gate FAIL
as scientific disagreement. All 25 orchestration protocol tests pass.
`global-clean-2` is running against frozen native/Python sources with the
unchanged one-hour whole-campaign limit. Native linked results remain pending.

Execution checkpoint: the corrected SPICE build succeeds. The complete
registered U-series is **15 PASS / 11 SKIPPED / 0 FAIL / 0 ERROR**; U16 now
also runs all 25 orchestration controls and writes mathematical-reference
CSV/PNG/EPS. Linked I01/I02/I03/I06/I07/I09/I10 PASS. The configuration-
separated no-SPICE, ion/nightly and eight-species callback stages still need
to finish before any combined package result is reported.

Subsequent native fixtures identified repository ion-name mismatches and a
pre-existing shared `He_PLUS` input typo: `mas =4.0*_AMU_` silently left its
mass zero in native preprocessing. The fixture now uses `NA_PLUS/HE_PLUS/
NE_PLUS/AR_PLUS`, and only the misspelled key in `input/species.input` was
corrected to `mass`, preserving the intended 4.0-AMU value. The new native
CSV checks actual generated masses/charges. Build identity now binds the
shared species include; its new retained-input regression passes (17 total
identity controls). Full CORE preservation is being executed independently
against all 843 previously passing cases; those original profiles do not
compile He+, so that audit alone cannot establish its corrected parameter.

The corrected actual eight-rank callback preflight PASSes **32 rows** across
all four neutral/ion mappings. Actual He+ mass is
`6.6421562663999999e-27 kg` and charge `1.602176565e-19 C`, as generated by the
existing AMPS constants/input model. It is rate/parameter wiring, not particle
kinematics or reaction-budget qualification. The preflight uses no SPICE
kernels at runtime; unused installed kernel inventory is explicitly labelled
build/configuration context. The final provenance helper hashes the actual
configured coefficient file and assignment-owning included input documents.
All 31 expanded orchestration controls pass; canonical validation data remain
read-only. Earlier global attempts and metadata errors are retained intact.

`global-clean-4` is the current final one-hour-budget campaign. Its ten M0
baseline gates PASS again after the shared-input correction, with zero-byte
nightly diff; callback and normal-build restoration are in progress. The
separate 51-suite/1085-case CORE campaign is in its final CME/portability
owners. Full aggregate completion and prior-PASS preservation remain pending.

Final `global-clean-4` campaign: **980.26 seconds (16.3 minutes)** whole-command
wall time. Registry: **1 PASS / 4 SKIPPED / 0 FAIL / 0 ERROR**; separate support:
**3 PASS / 0 SKIPPED / 0 FAIL / 0 ERROR** (complete M0 prerequisite, actual
eight-rank rate callback, ordinary Moon/SPICE restoration). The complete
U-series remains **15 PASS / 11 SKIPPED**, no FAIL/ERROR. Within U16, **170
kernel controls**, **17 native runtime controls** and **31 orchestration
controls** pass. Numerical maxima remain 1.5637e-15/1.5135e-15; mathematical
CSV/PNG/EPS figures are preserved. The native callback checks all **32
species/rank rows**, including actual generated ion masses and charges.
Peak sampled process-tree RSS is **5,354,135,552 bytes** (about 5.35 GB).

The separate full core audit completes in **2446.84 seconds (40.8 minutes)**
including setup/reporting, with peak **5,338,165,248 bytes**. All **843/843
previous PASS identities are preserved**, with no new failures, missing cases
or status changes. Full counts remain **843 PASS / 240 SKIP / 1 FAIL / 1 ERROR**;
the pre-existing strict-coverage FAIL and 180-second CME ceiling ERROR remain
explicit. Evidence: `core-species-regressions/core_completion_report.json`.
This is regression preservation, not an all-green core claim or independent
qualification of the newly enabled He+ table entry. The native M3 callback
supplies its separate parameter/rate evidence.

A final diagnostics-only Python extension creates stdout/stderr placeholders
when a runner fails before MPI and preserves the infrastructure reason. The
same 31 orchestration controls verify those paths, and all 17 retained-build
identity controls pass. Native source/physical-input bytes and computation
commands match the final campaign; the small post-campaign logging diff and
its controls are recorded separately rather than rewriting frozen results.

The approved initial M3 work package is complete; **M3 itself is not complete**.
The D04 photon-table production dispatch, conservative ion/product transfers,
competing reactions/ensemble budgets, measured electron forcing, molecular
state/joint daughter kinetics and charge-exchange activation remain to be
implemented or qualified. The completed literature proposals and `question.md`
give concrete next scopes for investigator review; no unsupported channel,
recoil distribution, bath assumption or source-range extrapolation was
silently selected. Formal I21-I24 and all M1/M2 predecessor requirements remain
intact. No observations/references were altered, and no commit/push occurred.

## 5.16 November 17 upstream integration and adapter repairs

The completed package preserves all 919 earlier PASS identities on the tested
core basis 4bff092. The final delivery target is 5ce92cb; Git-object and native
source-hash audits prove unchanged core/Moon/Enceladus/input/test/build sources.
Late changes were checked through the actual 64-case nightly tooling suite and
official registry check. Existing strict-coverage FAIL and CME 180-second ERROR
remain visible; no thresholds or physical references were replaced.

The linked I19 guard exposed a previously false-positive generated-mode check:
an early OFF definition was followed by an explicit application ON definition.
Moon's derived build input now changes both controls. Real compiler preprocessing
after the complete production include chain, four authenticated snapshots and
recipe/executable identity reject stale or contradictory retained builds. The
original failed attempt is retained; the corrective run meets the same original
budget and all 36 numerical criteria. Maximum relative atom imbalance in the
42-row I19 supporting ledger is 9.5924e-16; final Na/Ar count residuals, chemistry
event counts and numerical removals are zero. Formal I19 remains SKIPPED without
I13. The corrected retained ABI is accepted after ordinary Na restoration, while
strict reuse against that different live tree is correctly refused.

Enceladus's clone-dependent exact-byte exporter failure reproduced on pristine
upstream. Its owner now emits repository-relative provenance for owned inputs
and explicit absolute paths for external inputs. Existing byte comparisons stay
unchanged; two additional path controls run through the existing global owner.
All 62 historical Enceladus passes are preserved. No native Enceladus physics,
reference or science setting was changed.

M3's prior 32-row scalar-rate result retains its original build-adapter basis.
No common native source or canonical physics input changed afterward; only the
builder recipe was added to source identity. Old receipts were not backfilled,
and future retained reuse must satisfy the tightened identity. There is no new
M3 particle-chemistry qualification or December hold-out scoring.

Detailed manifests, exact commands, source/config snapshots, input/data hashes,
all original/retry outputs, restoration and delivery checks are retained under
`test_output/srcMoon/upstream-core-integration-20261117T185033Z/`. Active checkout
source was delivered only after completed tests and independent preservation
review. No commit, push or stash was made, and no unrelated local work was discarded.

## 5.17 M3 completion proposal — 2026-11-17

The user requested a detailed proposal, not activation of unresolved chemistry.
The new [M3_COMPLETION_PROPOSAL.md](test/M3_COMPLETION_PROPOSAL.md) inventories
production entry points and reusable upstream reaction/batching services,
defines staged deliveries and independent U14–U17/U24/I21–I24 verification,
and identifies missing data and decisions. Existing I21–I24 placeholders have
not been converted into executable runners by this documentation package.

The recommended engine uses one combined hazard owner and one correlated
parent/product transaction, individual weights and a common physical step.
It preserves legacy Na regression behavior. Normal Moon's force-aware
second-order tracer does not currently invoke the available first-order
trajectory-event hooks; the proposal specifies the smallest optional event/
sibling-continuation extension and full previous core PASS preservation if
that extension is needed. No core code was edited for this proposal.

Unresolved approvals/evidence include product recoil/electron/bath and mass/
energy conventions, charge-exchange channel/kinematic support, molecular
states/joint distributions/threshold tails, numerical and statistical design,
I22 total-versus-partial rate interpretation, U24's surface/chemistry milestone
overlap and measured-forcing scope. Ion fields/surface fate and the gas
solid-angle/origin operator are explicit contracts too; the M2 ion-impact
refusal is preserved. Atomic rates alone cannot close full M3.
Source and external roadmap hashes/differences are recorded without modifying
either roadmap. Formal M1/M2 predecessors and all existing criteria remain.

The read-only data audit finds 59/59 matching identities, including D04/D05
chemistry, charge-exchange research and the intact but empty/error D06 package.
No observations, processed science values or immutable raw files changed.
Retained audits, analytical sample-design examples and documentation checks:
`test_output/srcMoon/m3-completion-proposal-20261117/`. Sampling examples are
pre-output proposals, not model results or approved statistical thresholds.

Files changed: this status ledger, srcMoon/test and model READMEs, RUN.md,
the test-gap cross-reference and the new proposal. Source inspected:
`74974c2b2fc1b4cf09a73c320064a22370c0ca94`. No configure/build/model test ran;
there are zero new U/I evaluations and no qualification-state change. The
existing local Makefile/.gitignore edits are preserved. No commit/push/staging
was performed. The retained completion report records exact audit/document
checks, hashes and the unresolved physics/data/capability gates.

## 5.18 M3 implementation and local force-tracing checkpoint — 2026-11-17 23:31 UTC

Implementation is active against HEAD `74974c2` plus uncommitted source changes.
The user explicitly requires two particle-tracing orders with local physical
time, cut cells, lunar/Earth gravity, radiation pressure, Lorentz and Coriolis
forces. The new optional core movers preserve elapsed/remaining seconds,
epoch, event/face history and MPI continuation, with explicit numerical
controls and unchanged legacy defaults. Normal Moon dispatch selects the new
orders through `.in`; the initial application adapter explicitly requires a
frozen LSO environment and prescribed E/B. General core callbacks carry time;
an evolving Moon per-particle ephemeris/force provider is not silently assumed.

Implemented M3 services include the complete frozen D04 rate/channel inventory,
He/Na source-identified CEX rate kernels, conservative combined-hazard/event
plans and the actual PIC packet transaction adapter. No source-specific
daughter distribution is invented or enabled. Individual I21–I24 owners now
retain formal predecessor/capability refusals and exact input/manifest evidence;
they remain **gate-only owners**, not completed particle-execution/scoring
implementations. Their registration must not be called M3 qualification.

Completed supporting evidence at this checkpoint:

- Rate registries: 11 groups PASS; 103 independent CEX numerical evaluations
  have maximum relative implementation error `2.22e-16`. This is fit/software
  agreement, not the physical accuracy of the published He approximation.
- Photon selection: 12 families, 30 frozen branch claims and 1.2 million actual
  production selections PASS. The first statistical artifact exposed a SciPy
  degenerate-binomial quantile error; it is retained with an ERROR audit. The
  corrected oracle uses exact `[N,N]` bounds for single-channel families and
  unchanged seed/N/alpha for the final PASS. Observations were not edited.
- Reaction primitives: 56 controls PASS, including 100000 channel draws and
  local interval/weight/bath/state bookkeeping.
- Runtime input/installation: 30 native controls PASS, including table identity,
  transactional refusal and explicit first-/second-order tracing settings.
- Actual PIC storage/transactions: 828 rank/check evaluations PASS (92 on one
  rank and 736 across eight ranks). These are not 828 independent physical
  validation cases. Unequal dt/weights, rollback/deferral, fresh clocks and a
  real interleaved-ledger allocation fault are exercised.
- Actual Moon force/mover integration: 81 supporting checks PASS on one/eight
  ranks for both orders. All applicable force terms participate. Observed
  orders are approximately one/two; an RK4 comparison using the same production
  force verifies integration, while U03–U06 retain independent physical-term
  oracles. It is not new observational validation.

The last two results above belong to the **earlier 56-byte tracer record**.
The final record is 64 bytes and preserves last-face identity through zero-time
AMR/MPI continuation. A fresh normal producer build and native reruns are active;
old library/header bytes and results are retained, never relabeled as final ABI.

The six-body adversarial campaign found real failures beyond simple controls:
an AMR crossing coordinate `5.42e-20 m` from the exact shared plane caused a
zero-time owner loop; a rounded `D=v*h` classified a stored leaving velocity
as entering at a grazing face. Fixes project only solved AMR-face coordinates
within the declared geometry tolerance and apply a strict opt-in actual-velocity
filter; the original five-argument service/default behavior remains available.
The guard does not increase a limit or relax an acceptance tolerance. Accelerated
grazing contacts also require consistent incoming integrated velocity before any
state/clock mutation. Final fixed/adaptive, sphere/moving-plane and six-body
one/eight-rank checks continue. Earlier first-order six-body support completed
1350797 independent particle-step comparisons without mismatches; final ABI and
all requested variants still need their own results.

Native application build diagnosis also preserved infrastructure ERROR attempts:
the legacy typical B table is const (fixed through an immutable runtime field
provider, never const-cast); flat archive member naming replaced the Moon CEX
object (fixed with unique `MoonChargeExchange.o`); an unrelated configured
velocity sampler lay outside the analytical fixture mesh (disabled only in
that derived diagnostic recipe); and the thin MPI mesh probe lacked `MPI_ON`
(fixed to match the actual producer compile mode). These are not scientific
FAILs or discarded attempts.

Current core campaign has reached suite 48/51 with no new evaluated failures
reported so far. It retains its exact frozen source basis and original one-hour
deadline; later opt-in source corrections are recorded separately. A final
current-source preservation run is planned after source/native freeze, with its
own explicit <=3600-second setup/run/cleanup package. Previous known coverage
FAIL/CME timeout ERROR must remain visible; all 843 historical core PASS IDs
must be preserved. No unexecuted final result is claimed here.

Evidence directories:

- `test_output/srcMoon/m3-rate-registry-local-final2-20261117/`
- `test_output/srcMoon/u15-photon-selection-corrected-oracle-20261117T231600Z/`
- `test_output/srcMoon/m3-reaction-execution-20261117-verified/`
- `test_output/srcMoon/m3-force-cut-tracing-20261117T222159Z/`
- `test_output/srcMoon/m3-force-cut-native-20261117T225717Z/`
- `test_output/srcMoon/m3-full-core-20261117T225706Z/`

The user confirms no preferred missing kinetic references. Concrete primary
source/acquisition proposals are
[M3_KINEMATIC_REFERENCES_PROPOSAL.md](test/M3_KINEMATIC_REFERENCES_PROPOSAL.md)
and [M3_ION_SURFACE_AND_CEX_DATA_PROPOSAL.md](test/M3_ION_SURFACE_AND_CEX_DATA_PROPOSAL.md).
They distinguish measured tables, theoretical/simulation products, conditional
angular/energy support and author-request data. Missing joint kinetics, ion
surface fate and formal M1/M2 predecessors still prevent full M3 closure.
Canonical raw/processed observations, active root build/site configuration,
reviewed references and original numerical gates remain unchanged. No commit,
push or staging operation has occurred.


## 5.19 Final tracer package and active reaction-controller integration — 2026-11-17 23:58 UTC

The force-aware cut-cell work package has completed within its original
one-hour budget. Both production orders and the legacy six-body default pass
on one/eight MPI ranks. Each new order/layout performs 1350797 independent
particle-step comparisons, with zero mismatches and maximum position error
about 2e-12 m. Deterministic layout digests agree. The final opaque tracer
record is 64 bytes; the earlier 56-byte producer/results remain archived.

The committed core runner passes seven groups: ten reporting controls,
27 production numerical checks, 104 fixed and 114 adaptive checks on each
one/eight-rank layout (436 native criteria), and the expected fixed-grazing
refusal. Diagnostic rotation comparisons are retained as CSV/PNG/EPS.
Separate clean periodic, backward and self-consistent-field producers all
refuse through their specific production diagnostic. These PASS refusal
controls do not qualify those unsupported modes. Twelve final build/run
stages PASS, with no FAIL/ERROR/SKIPPED; earlier diagnostic attempts remain
preserved. Maximum recorded native-job RSS is approximately 4.14 GB.
See `test_output/srcMoon/m3-force-cut-native-20261117T225717Z/completion.json`
and the two core mover READMEs for exact commands, sources and build hashes.

The final 64-byte normal Moon producer has separately passed 828 PIC backend
rank/check evaluations and 81 actual-Moon-force tracer checks, with CSV/PNG/EPS
and source/generated-layout evidence under
`test_output/srcMoon/m3-force-cut-tracing-20261117T222159Z/particle-backend/`.
The completed combined backend/force package took 37.4 minutes, including its
preserved infrastructure attempts; it does not include later controller work.
These results verify the declared analytical packet/force integration scope,
not unknown ion, CEX or molecular daughter distributions.

The first complete CORE run preserved every one of the 843 prior PASS IDs:
1085 cases comprise 843 PASS, 240 SKIPPED, one existing strict-coverage FAIL
and one existing CME timeout ERROR. Its whole package took 2424.6 seconds.
That source snapshot predates the final optional-mode guards. The final
current-shared-source CORE run is active, with its original package deadline
2026-11-18 00:39:17 UTC; it must finish before claiming current-core preservation.
Later Moon-only controller edits are outside the core suite's compiled scope.
No current-core PASS is inferred from the earlier freeze.

The combined reaction controller is now being integrated at the normal
`Moon::MoveParticle` entry and the chemistry reinsertion callback. It requires
an explicit immutable constant-history hazard/topology configuration and
qualified kinetic provider, is disabled by default and rejects unimplemented
trajectory-dependent hazard integration. It retains optical depth through
refinement/migration and prevents duplicate legacy step-end photolysis when
enabled. Native controller tests and a new final Moon/M0/U09/global campaign
are still pending; this paragraph records implementation progress, not PASS.
Unknown physical providers remain unavailable rather than approximated.

Global U15/U17 now preserve actual supporting services and their separate
outcomes. I21-I24 owners remain gating/evidence owners, without full physical
execution/scoring. Protocol-only tests pass (36 orchestration, five formal-owner
and five registry controls) in
`test_output/srcMoon/m3-protocol-final-20261117T235715Z/`; those tests use
explicit infrastructure fixtures, not science data. Root Makefile/site config,
reviewed nightly references, raw validation data and qualification thresholds
remain unchanged. No staging, commit or push has been performed.


## 5.20 Frozen Moon source and final native campaign — 2026-11-18 00:19 UTC

The final Moon development package started at00:18:45UTC, with an absolute
01:18:45UTC deadline including clone/prerequisite setup, clean builds, native
runs, comparison figures, aggregation and restoration. The active checkout
remains unconfigured/unbuilt. Private source snapshot:
`/tmp/amps-m3-final-moon-20261118T001845Z/candidate`; durable evidence:
`test_output/srcMoon/m3-final-moon-20261118T001845Z/`. This package is ACTIVE,
not yet PASS. It requests the existing global M3 owner with eight MPI ranks,
OMP1 and the explicit disposable-checkout flag. Reviewed nightly remains four
ranks. The newly registered controller uses the same fresh all-species
producer as the real PIC packet and actual-Moon-force fixtures.

Independent review found and corrected additional production seams before
this snapshot:

- Moving-LOLA finder previews must not publish a collision frame or increment
  a physical collision counter. New-tracer queries are now pure; the actual
  selected impact callback commits its frame/count from its accepted normal
  and particle epoch. Legacy behavior remains unchanged.
- The explicit frozen tracer's source pose and material velocity must agree
  with the boundary pose. Source transforms and initial/updated global body
  matrices now share the frozen-only derivative-block policy for both sphere
  and LOLA. Actual birth offsets remain physical clock information. U09 now
  registers a real-SPICE production transaction regression; that new linked
  mode is pending the fresh build, not claimed passed from review.
- A gas volume point has no self facet. A finite ignored `sourceFace` is now
  refused, and U09 adds a genuine blocked gas-disc BVH case plus transactional
  refusal. This prevents a physical blocker from being accidentally omitted.
- Generic source boundary injection previously bypassed the selected Moon
  wrapper through the old boundary tracer. Both configured boundary macros
  now route to `Moon::MoveParticleBoundaryInjection`; legacy forwards the exact
  old call, while new tracing prepares clocks and applies its selected force
  order to the remaining birth interval. Source epochs use the scoped moving
  departure or the stationary source's already sampled remaining seconds.
- New local tracing explicitly refuses `m2_environment=enabled`: frozen global
  ephemeris gives zero thermal-step duration, and no shared evolving thermal/
  reservoir clock is defined for asynchronous local particle histories. No
  prescribed-stock subset or physical time mapping was assumed. The parser
  and programmatic installer both enforce this scope;34 native configuration
  controls PASS in `test_output/srcMoon/m3-runtime-temporal-guard-20261118T001600Z/`.

The normal controller is source-frozen, default DISABLED and constant-history
only. Its separate native package began00:17:56UTC with deadline01:17:56UTC in
`test_output/srcMoon/m3-reaction-controller-20261118T001756Z/`; it performs a
fresh producer build plus one/eight-rank, both-order normal-dispatch and
specific capacity-refusal checks. It has not completed. The final global
orchestration now actually invokes this runner after actual Moon-force checks;
38 protocol controls PASS in
`test_output/srcMoon/m3-controller-orchestration-20261118T000817Z/`.
Missing physical kinetic providers and formal M1/M2 predecessors remain open.
The final shared CORE run remains separate and source-frozen; no PIC source
changed during these Moon-only corrections. No stage, commit, push, reference
replacement or canonical data modification has occurred.


## 5.21 Completed shared-core audit; controller endpoint failure retained — 2026-11-18 00:27 UTC

The final shared-core audit is COMPLETE. Exact protected-case comparison
preserves843/843 prior PASS identities, with no missing case and no newly
failing case. Actual1085 outcomes remain843PASS,240SKIPPED, one known strict
coverage FAIL and one known CME0628 timeout ERROR. This is regression
preservation PASS, **not an all-green suite**. All12 shared freeze hashes
match both tested candidate and active source. Whole package2419.97s
(40.3min), guarded run2094.0s, peak owned RSS5.571GB; the original deadline
was never reset and the HOLD-interrupted earlier attempt remains an ERROR.
Receipt:
`test_output/srcMoon/m3-full-core-final-20261117T233917Z/guarded-attempt/current_shared_core_preservation.json`.
Raw report SHA256:
`94514896ea7a310b6acc4bb94092b6e23d83a5f2254ea7a0a908070c69337835`.
Protected exact-ID list SHA256:
`fc52fd33f3f81131ae5fa27710a454f27f5d6cf81a79db14f4d4c5fbdfd0f3a0`.
Later Moon/controller/doc changes are explicitly outside this compiled core
scope; no new shared PIC change is proposed.

The first fresh normal Moon build and U09 PASS, including real-SPICE
moving-surface transaction and the new gas-volume blocker/refusal controls.
The complete local series records15PASS,11SKIPPED,0FAIL/ERROR; partial U15/U17
retain actual component evidence and their full scientific SKIPs. The linked
baseline/nightly/remaining native global campaign is still active.

The separate fresh PIC backend passes828 checks. Its first controller run
completes with53PASS,8FAIL,0ERROR/SKIPPED: first order passes both layouts;
second order leaves the endpoint event's parent alive on both layouts. It
reports epoch1001s but residual interval1.1102230246251565e-16s and optical
depth1.1449174941446927e-16. The specific capacity-refusal control passes.
Original results, candidates, particles and plots remain under
`test_output/srcMoon/m3-reaction-controller-20261118T001756Z/controller/`.
This is an evaluated numerical FAIL, not an observational discrepancy or a
missing-data SKIP.

The first failing invariant is inconsistent interval partition arithmetic:
controller exposure accumulates accepted binary64 substeps, while the core
uses subtractive remaining time and publishes the original requested endpoint.
A Moon-only interval-anchored correction is being developed, using the original
immutable interval/optical clock and the core remaining-time coordinate.
It must not introduce an optical-depth tolerance or weaken endpoint criteria;
strictly future nextafter events must remain pending. Independent review and
new partition/endpoint controls precede a fresh native build. Old results
will not be relabeled, and any corrected global attempt must consume only
the remaining time before the original01:18:45UTC deadline, including setup,
rebuilds, comparisons and restoration.


## 5.22 Corrected native mechanics pass; final global/core campaigns active — 2026-11-18 00:56 UTC

The final controller package completes in1732.6s (28.9min), including setup,
three clean builds, the retained original failures, corrective checks, plots
and reporting. Its original01:17:56UTC deadline was not reset. The matched
final producer passes828 real PIC backend checks,89 actual Moon-force tracer
checks and73 normal-path controller checks on1/8MPI and both orders. Local
reaction primitives pass65 controls; resource-result classification passes
three controls. All664 non-test native source identities match the active
source and producer; six active configuration files remain byte-identical.
Report:
`test_output/srcMoon/m3-reaction-controller-20261118T001756Z/completion_report.json`
(SHA256 `b95ddd98ff15750a5a67ebfab88a611aa22701bcb6d93191ba6175a299273bae`).

The optical-clock correction anchors combined exposure to the immutable
interval-start clock and actual core remaining-time coordinate. A shared
production long-double exposure evaluator serves both primitive and controller
paths; compensated relative-time arithmetic avoids cancellation. There is no
optical-depth epsilon. Exact endpoint events fire, strictly future nextafter
events survive, tiny pending clocks continue through the next outer interval,
and queries/retries redraw nothing. Older53PASS/8FAIL outputs remain unchanged.

A separate new control then reproduced exact AMR-face/event/end coincidence
on the old core (69PASS/4FAIL, one failure per order/layout). The final optional
core dispatch processes an already due zero-horizon event after accepted
internal-face projection and before zero-duration publication. Current-owner
processing retains core ownership of parent deletion and all zero-time
children; Continue uses post-event species/velocity for administrative handoff.
No event time/position is nudged. Material-surface and external-boundary reaction
ties retain their previous priority and remain unqualified physical cases.
New shared tracer SHA256:
`33eefb1de600960bae98751e56b05e53e32e54093f767ad0ac6dc99e7c4c2f77`.
Public API and64-byte tracer record remain unchanged. The fresh core runner
passes516 fixed/adaptive criteria (124/134 on each1/8layout),27 math and10
reporting controls plus the expected grazing refusal. Both new orders also
pass all six-body criteria on1/8MPI:1350797 independent particle-steps per
run, zero mismatches, maxposition2e-12m. Legacy default runs are still active.

The Python audit also fixed real runner defects: owned descendants are reaped
on supporting-service timeout; one deadline covers all groups; nested failures,
counts and exit contradictions cannot hide behind an outer SKIPPED; retained
M0 ERROR prevents dependent admission; the supervisor leaves time for the
worker's ordinary-build restoration; resource refusal requires exact exit1,
not a signal plus matching text. Eleven new bounded protocol controls and54
combined orchestration controls pass. All are registered through existing U16;
no permanent ID or physics criterion changed.

The first whole global Moon attempt completed in1276.4s with all ten M0 gates
PASS, zero-byte reviewed nightly diff, backend828PASS and actual force89PASS;
its controller53PASS/8FAIL and aggregateFAIL remain preserved. The corrected
whole global attempt started about00:42UTC using only the remaining2205s
before the ORIGINAL01:18:45UTC deadline. Its local series again records
15PASS/11SKIPPED/0FAIL/ERROR; linked/native/restoration work continues. The
latest Python-only delta was independently run through the current ROOT global
entry with `--skip-build --only U15 --only U16 --only U17`:96.6s, U16PASS,
U15/U17SKIPPED, four supporting service groups PASS, no FAIL/ERROR. Its U16
actually runs all54 protocol controls; native helpers report65 primitive and
34 runtime-input PASS. This focused check performs no active AMPS configure or
build. It is retained under the original package's
`latest-wrapper-global-support/` and `latest-wrapper-supervisor/`.

The new shared-core change requires a fresh CORE audit; the earlier843-ID
preservation PASS belongs to the older source. The new package began00:42:59UTC,
with separate fixed01:42:59UTC deadline including setup, and is active under
`test_output/srcMoon/m3-core-amr-endpoint-20261118T004259Z/`. It binds13 shared/
fixture identities and protects all843 historical PASS IDs. No current-source
CORE result is claimed before its aggregate exists.

M3 remains EXPERIMENTAL / NOT QUALIFIED. Implemented mechanics do not supply
physical Na/electron recoil, CEX joint state/angle/projectile distributions,
molecular daughter/state cascades, variable geometry-dependent exposure or ion
surface fate. The initial Moon tracer still requires frozen external inputs
and `m2_environment=disabled`; no evolving thermal/reservoir clock mapping is
inferred. Formal M1/M2 predecessors and I21-I24 physical gates remain unchanged.


## 5.23 Global Moon package complete within original hour; separate CORE pending — 2026-11-18 01:19 UTC

The corrected global development owner completes with10 supporting PASS,
0 supporting FAIL/ERROR/SKIPPED. The registered result remains correctly
SKIPPED: U16PASS, U15/U17 and formalI21-I24 SKIPPED (onePASS/sixSKIPPED), not
M3 qualification. The current-source M0 linked baseline is10PASS/0FAIL/ERROR,
its reviewed `make -j test_Moon` diff has zero bytes, and the full local series
is15PASS/11SKIPPED/0FAIL/ERROR. Actual backend828, Moon-force89 and controller73
checks pass. Ordinary Moon/SPICE restoration passes.

Original setup, the failed first attempt and complete corrected native attempt
finish in2819.3s (47.0min), before the original01:18:45UTC deadline. Latest-root
Python wrapper controls, metadata audit and completion-report generation also
finish within the same original hour: report at3397.1s, final source checks at
3429.1s. No child received a fresh hour. The separate standalone controller
and mover packages have their own explicitly recorded bounds. Their results
must not be summed as independent physical validation cases.

Final global report:
`test_output/srcMoon/m3-final-moon-20261118T001845Z/completion_report.{json,md}`.
It records exact commands, modified code/docs,664 matching native source
identities, compiler/MPI/build/generated ABI, data-hash receipt links, all
retained failures and scope limits. `final_source_checks.json` records27 Python
syntax checks, eight JSON parses,76 new text-file whitespace checks, no
protected root changes and no staged paths. The later documentation-only
clarifications change no native source or numerical gate.

The fresh endpoint mover package is also complete:14 build/run stages PASS;
seven committed runner groups include516 fixed/adaptive,27 math,10 reporting
and expected grazing-refusal controls. Both new orders and the legacy default
pass all72 six-body criteria on1/8MPI,120steps each. All four new-tracer cases
compare1350797 independent particle-steps, zero mismatches, maxposition2e-12m;
layout digests match each other and the earlier package. Its preparation plus
execution/report take23.2min, peak4.14GB. Receipt:
`test_output/srcMoon/m3-force-cut-endpoint-20261118T003500Z/completion.json`.

A further read-only production audit documents the remaining sampling/time
scope precisely. Canonical Na-only **legacy** operation uses a spatially
global species step and keeps its existing snapshot meaning. Normal cell and
Moon LOS accumulators average occupancy by iteration, without particle-epoch
or residence-time normalization. The new frozen setup repeats external step
and within-step birth epochs while survivors advance locally, so it does not
establish a common advancing transient ensemble. General asynchronous
reactive density/LOS/brightness/flux science is therefore unqualified; this is
not a claim that every canonical snapshot is wrong. No estimator/time mapping
was invented. Production snippets and the scope analysis are retained in
`sampling_scope_audit.md`; README, capability register and gap proposals now
state the explicit limitation. Existing common-time M2 fixtures remain a
separate declared scope.

The current optional core hash remains33eefb1..., and its separate full CORE
package is still running within00:42:59–01:42:59UTC. Earlier843-ID preservation
receipts remain attached to their earlier source; this checkpoint does not
infer the final current-core result. The user-facing final status will be
updated when that actual aggregate exists. No commit, push, staged file,
active-checkout configure/build, reviewed-reference change or raw-data write
has occurred.


## 5.24 Final current-source verification and remaining M3 work — 2026-11-18 01:25 UTC

Final shared-core regression preservation is PASS on tracer
`33eefb1de600960bae98751e56b05e53e32e54093f767ad0ac6dc99e7c4c2f77`.
All843 protected PASS IDs remain PASS across four actual baselines; no case
is missing, newly failing or status-changed. Actual1085 outcomes are843PASS,
240SKIPPED, one existing strict-coverage FAIL and one existing CSWC0628
mass/time timeout ERROR. This is **preservation PASS, not an all-green routine**.
All13 frozen shared/fixture hashes match both candidate and active checkout.
No native production bytes changed after this source freeze; later changes
are documentation/Python-only and separately identified.

Receipt:
`test_output/srcMoon/m3-core-amr-endpoint-20261118T004259Z/current_shared_core_preservation.json`
(SHA256 `c13e95f9b0e34b8123e5d1ce8915c4cd9db557683adde49dcaca4f2f7dcc0357`).
Raw1085-case report SHA256:
`08ba768fc8448a4db0b3b141f5ad5386093a39a1e67a612e2a36fbd533ba4051`.
Full package including setup/aggregation is2458.0s (40m58s); native runner
2258s (37m38s), peak owned RSS5574389760bytes (5.57GB). Its original separate
00:42:59–01:42:59UTC budget was preserved. The earlier shared-core receipts
remain attached to their own freezes.

The separate Moon global package is complete as recorded in5.23: ten M0
PASS/zero-byte nightly diff;15 localPASS/11SKIPPED; ten supportingPASS; formal
I21-I24 SKIPPED. Matched native checks are828 backend,89 actual-force and73
controller PASS; local reaction/input/protocol controls are65/34/54PASS.
The mover package includes516 fixed/adaptive,27 math,10 reporting and72
six-body criteria PASS, plus specific expected refusal controls. Counts retain
their declared local/rank/group scopes and are not summed as independent
physical validation cases. Each package includes setup/builds/plots/cleanup
within its own explicit <=3600s limit; the global failed+corrected native
attempts finish in47.0min, and its reporting/source checks finish within the
same original hour. No physical horizon, confidence level or numerical gate
was shortened or relaxed to fit.

The global completion report is frozen with its truthful status at report
creation (when the separate CORE run was still active):
`test_output/srcMoon/m3-final-moon-20261118T001845Z/completion_report.{json,md}`.
The final separate CORE receipt above closes that pending regression work;
it does not alter the report's original deadline or imply formal M3 closure.
All build logs, executable/generated ABI/source identities, input snapshots,
external-data hash receipts, comparison CSV/PNG/EPS and original failed
attempts remain discoverable through those reports.

Remaining full-M3 implementation/physics work is explicit:

- Acquire/review actual spectral/state/angular recoil and bath-energy evidence
  for atomic ionization, molecular channels/cascades and selected CEX support;
  total rates and mean energies do not determine these distributions.
- Install and verify those physical providers through the conservative normal
  path; generic constant-history analytical dispatch is not their substitute.
- Define/implement the trajectory-dependent exposure contract and its independent
  geometry/time/conservation/convergence tests. Current controller refuses it.
- Qualify ion surface fate and physical surface/external-boundary reaction ties.
- Declare stationary residence-time versus common advancing transient sampling
  and source-clock semantics before producing asynchronous density/LOS/
  brightness/flux science. No estimator or temporal mapping has been assumed.
- Close the unchanged formal M1/M2 predecessors and each I21-I24 physics/data/
  numerical gate before claiming M3 VERIFIED/VALIDATED or advancing qualification.

The concrete acquisition routes, source validity and remaining missing arrays
are [the kinematics proposal](test/M3_KINEMATIC_REFERENCES_PROPOSAL.md) and
[the ion-surface/CEX proposal](test/M3_ION_SURFACE_AND_CEX_DATA_PROPOSAL.md).
Downloaded files remain research evidence; no unavailable joint distribution,
source metadata, observational uncertainty or hold-out result was invented.
The existing U01-U26/I01-I46 meanings, raw data, reviewed references and formal
acceptance criteria remain unchanged. No active-tree build/configuration,
staged file, commit or push was made.


## 5.25 Completed independent software/model preparation — 2026-11-18

This package starts from HEAD `74974c2b2fc1b4cf09a73c320064a22370c0ca94`
plus the previously recorded uncommitted development. The user's concurrence
covers independent software and concrete model/protocol specifications, with
common physical time as the initial qualification direction. It does not
supply missing recoil distributions, state populations, mass/energy conventions
or ion-surface laws. Those remain BLOCKED/null; formal M3 remains open.

The source/data roadmap identities remain respectively
`712e1152bf50af6d9ee1dad8f8a585809ca8b92cf42ac6fc7a009f6c0253ee17` and
`d10eca1c539f9cf82940236d9e1aa36af665d43e82ba475925950d582fbf49c4`.
Source development/resource exceptions control implementation order, native
labels control data interpretation, and unchanged formal criteria/predecessors
control qualification. Neither roadmap was edited. The implementation contract
is `test/M3_INDEPENDENT_SOFTWARE_CONTRACT.md`.

### Implemented capabilities and exact scope

- `species/PiecewiseHazard`: declared affine histories in relative seconds,
  analytic dimensionless integrals, exact rational supplied-binary crossing
  comparisons, explicit support gaps/jumps, first representable roots and
  approved left-limit endpoint rates. No observational interpolation or RNG.
- `species/ReactionExposure` and the shared `ReactionExecution` selector:
  immutable preview/retry, compensated optical debit and instantaneous event
  channel probabilities. A pending compensated depth that cannot be represented
  exactly by the scalar exposure interface is refused transactionally. The
  controller remains constant-history-only; no physical variable provider was
  installed by the pure service.
- `MoonTransportClock`, `MoonTracing`, `MoonInput`, production source epoch
  helpers and normal `amps_time_step`: optional `tracing_clock_mode=common-time`,
  explicit normal UTC/SPICE origin, compensated relative time, frozen external
  epoch, actual source departures and collective Begin/Commit. Every actual
  block/species dt must be identical and constant. MPI/clones retain common
  targets. Unequal steps, unfinished packets and unqualified restart are refused.
  Omitted mode preserves independent/legacy behavior.
- `M3_PHYSICAL_MODEL_SPECIFICATION.{md,json}` and
  `M3_PHYSICAL_DATA_REQUIREMENTS.csv`: complete D04/D05/CEX inventory, physical
  state/mass/energy/bath/kinetic evidence requirements and genuine gaps.
- `M3_QUALIFICATION_PROTOCOL.{md,json}` and `m3_protocol.py`: identity/contract
  preparation, frozen existing thresholds, separately proposed statistics and
  numerical budgets, and explicit full-molecular/U24/measured-I23 dependencies.
  The proposed numerical steps are T/8, T/16 and T/32 to one representable final
  horizon; the statistical half-survival calculation is a planning example.
  Preparation is structurally PASS and formally SKIPPED, never an executor.

The normal library compiles the same services used by thin probes. Global U15
includes seven supporting groups; U16 discovers the orchestration/protocol
controls. The extended development owner runs the same fresh PIC producer,
force/controller/common-clock fixtures and full normal startup after ordinary
restoration. No permanent U/I ID or acceptance threshold changed. Application,
species/test READMEs, RUN.md, INPUT_FILE.md and the gap registry document the
syntax, commands, scope and failure behavior.

### Actual verification results

| Evaluated scope | PASS | FAIL | ERROR | SKIPPED |
|---|---:|---:|---:|---:|
| Clean full U01–U26 inventory | 15 | 0 | 0 | 11 |
| Linked I01–I10 baseline | 10 | 0 | 0 | 0 |
| Actual PIC packet backend | 828 | 0 | 0 | 0 |
| Actual Moon force/tracer integration | 89 | 0 | 0 | 0 |
| Constant-history conservative controller | 73 | 0 | 0 | 0 |
| Extended common-clock/source/MPI/real-daughter composition | 80 | 0 | 0 | 0 |
| Full normal initializer + two real steps/order at eight ranks | 32 | 0 | 0 | 0 |
| Declared exposure mathematics, 84 fixtures | 634 | 0 | 0 | 0 |
| Exposure/clock bridge including precision refusal | 25 | 0 | 0 | 0 |
| Existing reaction primitives after shared-selector refactor | 65 | 0 | 0 | 0 |
| Local common clock + actual parser | 51 | 0 | 0 | 0 |
| Latest orchestration/protocol controls, including nested controls | 81 | 0 | 0 | 0 |
| Latest registered U15/U16/U17 top-level scope | 1 | 0 | 0 | 2 |
| Actual formal I21–I24 capability/predecessor owners | 0 | 0 | 0 | 4 |

Rows describe distinct retained runs and nested scopes; they must not be summed
as independent trials. U15/U17 remain formally incomplete despite eight actual
supporting groups PASS. The 81 orchestration controls include 17 qualification
protocol and five clock path/error controls; a further five formal-owner
infrastructure controls pass. Source-audit documentary checks are not science
qualification.

`make -j16 ... test_Moon` passed against the unchanged reviewed I04 reference
`ed875cefc604e4ae3580fb89b8c528146929b6b9be099b41b61697faee7f0080`;
`test_Moon.diff` is zero bytes. Normal clean Config.pl builds used SPICE,
separate no-SPICE/ion fixtures, and `-O3 -g`; eight MPI/OMP1 is the linked
baseline, explicit one-rank analytical controls are retained, nightly uses
its reviewed four-rank layout. Every generated build was removed before
configuration in the disposable checkout.

Full normal startup uses each complete common-time input, the actual normal
UTC-to-SPICE initializer and two actual `amps_time_step` calls. Both orders
use dt `6.85586530104012581432 s`, end at relative
`13.7117306020802516286 s`, keep external ET
`286027266.18458545208 s TDB`, and increment the sampling counter to 2.
Actual fixed-seed source emissions are 8 then 97, leaving 105 gas particles.
These are chronology/source-wiring diagnostics, not converged physical or
observational results. The observer chains the existing production ledger.

The additional common-time/controller case uses a declared nominal-mass,
zero-Q analytical packet only. Actual Na+/He daughters inherit the parent's
common target through the real backend clone, consume their unused interval,
retain fresh independently drawn optical clocks and finish at the common
endpoint with verified integer/bath/kinematic ledger. It supplies composition
coverage; no physical CEX kinematic law is inferred.

### Preserved failures, resources and provenance

Independent review found three preparatory issues: lossy projection of a
compensated optical depth, permissive/malformed JSON handling and admission of
an empty probe receipt. All were reproduced, fixed and independently checked;
no tolerance or scientific threshold was relaxed. Missing-compiler/plot and
strict-compiler development attempts remain preserved as their actual outcomes.

The first full native package is explicitly ERROR: its thin common-clock
runner passed JSON path strings into a Path-only hasher before native clock
execution. Its M0/backend/force/controller PASS receipts and ordinary restore
remain intact. The test-only path normalization and runner-boundary error
artifacts were corrected. A byte-identical saved backend/archive/header tree
was used for fresh clock controls, then the exact ordinary tree was restored.
The corrected continuation passed 76 hook controls and 32 full startup checks;
the later 80-check extension adds the actual daughter/common-target case.
The failed attempt is not relabeled PASS or overwritten. The original full
global invocation remains ERROR; corrected stages are authenticated same-deadline
continuations. Current orchestration also has passing protocol controls; a second
complete global invocation was not performed after those test-only fixes.

The original end-to-end native package starts 07:29:58 UTC and has the unchanged
08:29:58 UTC deadline, including setup, clean builds, failed attempt,
restoration, corrected checks, plots/report and cleanup. Actual test execution
finished before that deadline. Peak monitored owned-process RSS is
15,733,587,968 bytes, below the 200,000,000,000-byte limit. Each run retains
commands, stdout/stderr, result, source/input/build/executable identities,
actual furnished SPICE hashes, MPI/seed data and resource history.

Evidence roots:

- `test_output/srcMoon/m3-independent-20261118T070304Z/`: baseline, local/focused
  registry results, independent review and final completion/source-delta report.
- `test_output/srcMoon/m3-independent-native-20261118T072958Z/`: original raw
  package, successful M0/native receipts, preserved ERROR, corrected continuation,
  full normal startup and common-controller composition artifacts.
- `test_output/srcMoon/m3-piecewise-hazard-20261118T070200Z/`: independent exposure
  references/CSV/PNG/EPS and strict native mathematical checks.
- `test_output/srcMoon/m3-physical-model-spec-20261118T013000Z/`: reproducible
  read-only physical source audit and requirements documentary checks.

The pre-existing `U21_argon_source:Q!` directory contained exactly the three
HEAD files. The user approved restoring `U21_argon_source`; all file hashes
are unchanged and the global registry works. All seven protected active-root
configuration/Makefile files remain unchanged. All 13 shared-core identities
still match the earlier 843-PASS preservation audit; the full CORE suite was
not rerun for this Moon-only package. No native production source changed
after the clean-build snapshot. Later thin-test/protocol/documentation changes
are explicitly identified; src/INPUT_FILE.md changed only as documentation.
No staged file, commit, push, raw-data edit or reference replacement occurred.
No model/build/test process remains active.

### Remaining physical and formal blockers

The read-only audit matches 11 processed files and 54 provenance/raw identities,
but rates/mean energies do not supply joint daughter kinetics. Two positive
atomic samples below their listed thresholds remain source-interpretation
questions; no clipping or altered scalar rate was applied. Required state,
isotope, mass/formation/internal-energy and bath conventions, photon/electron/
CEX joint distributions, complete molecular thresholds/cascades and ion-surface
outcomes remain unavailable/unselected. D06 still lacks measured electron
forcing. The pure exposure service is not physical trajectory-provider wiring;
its compensated-depth storage refusal remains explicit. Asynchronous science
estimators, evolving M2 coupling, physical boundary reaction ties and complete
common-time checkpoint/restart are unqualified. Same-campaign M1/M2 and each
unchanged I21–I24 predecessor still control formal closure. The current stage
is **M3 implementation/supporting verification, EXPERIMENTAL**.

## 5.26 Chemical database and Hartley2 research — 2026-11-18

The user requested online chemical-rate sources and an audit of comet H2O
photochemistry. Three parallel read-only investigations covered the reachable
Hartley2 implementation, molecular databases and atomic/electron/capture data.
The consolidated report is `test/M3_CHEMISTRY_DATABASE_RESEARCH.md`; complete
source/download manifests are linked there. Source revision and both roadmap
hashes are unchanged from section5.25. This is research; no physics source,
compiled input, chosen rate model, acceptance criterion or canonical data was
modified, and no simulation/build/regression was run for this package.

Concrete findings:

- The existing D04 baseline supplies30 quiet/active rate/mean-energy rows for15
  reactive channels. A new read-only audit verified all11 processed-product
  hashes and30 exact catalog values. Seven active H2O rows were also checked
  against native TableB1 of the acquired Carberry Mogan2022 paper, including
  the stated5.2AU-to1AU normalization and printed rounding. They agree. This
  is a same-Huebner-source transcription check, not independent validation.
- Reachable `input/i3.conf`/`input/i3.in` provides17 channels, including9 neutral
  ionization sinks; the canonical hartley2.input has chemistryOFF. The empty
  source generated header is a sentinel. I3 already uses the Moore/Mapaye2025
  activity law. Its fixed water recoil energies4.04/3.94/.697eV correspond to
  an active source, while runtime changes rates only. Rubin2011 uses a
  separate1.7/3.8/.7eV choice. Their conditions must remain explicit.
- Existing two-body recoil math is a candidate reduced model after reviewing
  internal states, mean/translational energy, isotropy, photon momentum and
  validity. Three-body Gaussian/projection/rescaling is not generally uniform
  microcanonical sampling for unequal masses. The cometW/dt resampling,
  parent internal-byte inheritance, collapsed states, legacy mass tolerance
  and neutral-only ionization sinks cannot be imported unchanged into Moon.
- Leiden supplies a rate table, spectrum and partial cross sections; its
  documented normalized-Solar conversion is37700. Those are separate spectra
  from D04. Native waterHDF5 has negative values and nonclosing aggregate
  partials; no clipping or renormalization was applied. UMIST/KIDA ISRF rates
  and collapsed product networks are not solar1AU state-resolved sources.
- New acquired atomic/electron sources include NIST waterBEB andHe energy-
  differential tables, OPEN-ADAS neutral/metastable coefficients and an
  evaluated He capture table at0.3–500keV/u, Phelps Na3s/3p nodes and a measured
  IAEA partial Na capture channel. IAEA D107581 is composite capture+ionization
  and is not a pure-CEX replacement. Complete outgoing angular/state closure,
  domain/quality and redistribution terms remain explicit source constraints.

An existing Markdown transcription error in the physical specification was
corrected from Ar active9.07e-7 to7.74e-7s^-1, matching its unchanged frozen
JSON, processedCSV and original response. Previous document bytes and exact
before/after hashes are retained. No data/model/reference was changed.

Research evidence roots:

- `test_output/srcMoon/m3-database-hartley-research-20261118T083645Z/`
- `test_output/srcMoon/m3-hartley-h2o-code-audit-20261118T083318Z/`
- `test_output/srcMoon/m3-molecular-databases-20261118T083551Z/`
- `test_output/srcMoon/m3-atomic-online-databases-20261118T083400Z/`

No new U/I result or milestonePASS is assigned. New downloads remain
ACQUIRED_RESEARCH_NOT_READY; existing successful M0/native results retain
their original scope and were not rerun. The next development proposal is to
reuse the available rates and selected recoil mathematics through Moon's
conservative physical-time path after explicit physical-model review. M1/M2
formal predecessors and fullM3 product/surface/kinetic requirements remain.
No active-root build/configuration change, staging, commit or push occurred.

## 5.27 Concrete daughter-physics and ion-surface decisions — 2026-11-18

The user requested a concrete proposal identifying exact approval choices.
[M3_DAUGHTER_AND_ION_SURFACE_DECISIONS.md](test/M3_DAUGHTER_AND_ION_SURFACE_DECISIONS.md)
and its JSON register specify nine PENDING decisions; these labels are not
new U/I IDs. Source basis remains
`74974c2b2fc1b4cf09a73c320064a22370c0ca94` plus existing working changes.
Source/data roadmap hashes remain respectively
`712e1152bf50af6d9ee1dad8f8a585809ca8b92cf42ac6fc7a009f6c0253ee17`
and `d10eca1c539f9cf82940236d9e1aa36af665d43e82ba475925950d582fbf49c4`.
Approved source development/resource exceptions govern work order; common
formal gates and native labels retain their authority. No roadmap changed.

Recommendations are complete D04 endpoint/reference-population selection;
consistent additive isotope/electron reference masses and formation/excitation
energy; opt-in mean-energy/isotropic two-body recoil; separately opt-in
mass-weighted microcanonical three-neutral breakup; explicit product states
and supported cascades; prescribed electron-reservoir accounting while
unsupported impact products remain OFF; a bounded measured Na partial-CEX
beam experiment with an explicitly forward angular approximation; strict
ion-boundary behavior plus a verification-only ideal collector; and a
source-specific He+/silica conditional charge-survival benchmark. The
optional reduced models are proposed experimental closures, not acquired
distributions or qualified lunar physics. Missing channels remain in the
rate family; incomplete full-family transport requests are refused.

The proposal gives equations, source/validity limits, production mappings,
missing data, implementation order and independent verification oracles.
It distinguishes charge fractions conditional on scattering from incident
retention, and secondary sputtering yields from the projectile-fate
probability partition. Ion absorption into the proposed artificial collector
does not imply neutralization, implantation depth, heat or later release.
He laboratory results are not generalized to Na+/Ne+/Ar+ lunar mixtures.
Author-data requests are specified but no contact is authorized or sent.

Existing source/native evidence was inspected, with targeted primary-source
web checks of PHIDRATES methods, NIST electron/ASD data, versioned ATcT and
He/silicate publications. Some direct publisher/IAEA requests remain
unavailable; retained native packets supply the existing numerical evidence.
No data value was generated, filtered, repaired or promoted to READY.
Documentary schema/link/hash checks are separate from science results.

Evidence: `test_output/srcMoon/m3-daughter-ion-proposal-20261118T161952Z/`.
This package changes only proposal/README/STATUS documentation, including a
machine-readable proposal register. The frozen physical-specification and
qualification JSON, all acceptance/reference files, active-root build/site
configuration and existing production sources remain unchanged. Science
tests/builds executed:0; no new PASS/FAIL/ERROR/SKIPPED science count is
assigned. Previous M0/native results retain their original scope and were
not rerun. Formal M1/M2/M3 status is unchanged. No stage, commit, push,
canonical data write, runtime policy activation or observational scoring
occurred. Approval of these recommendations will authorize only the stated
scope; it cannot supply missing physical data or close M3 automatically.

## 5.28 Approved daughter/ion-boundary components and native verification — 2026-11-18

The user approved all P01–P09 recommendations. The implementation follows
`test/M3_APPROVED_MODELS_IMPLEMENTATION_CONTRACT.md`, written before code
changes, and records approval/completed scope in the decision JSON. The source
basis remains `74974c2b2fc1b4cf09a73c320064a22370c0ca94` plus existing working
changes. Three parallel agents handled recoil, collector and laboratory/source
subpackages; the parent reviewed numerical edge cases, integrated actual
callbacks and ran the isolated native package. Neither roadmap changed:
source SHA256 `712e1152bf50af6d9ee1dad8f8a585809ca8b92cf42ac6fc7a009f6c0253ee17`,
data SHA256 `d10eca1c539f9cf82940236d9e1aa36af665d43e82ba475925950d582fbf49c4`.
Source development/resource exceptions control work order; common formal gates
and native data metadata still control qualification and interpretation.

### Implemented and deliberately limited scope

- P01/P02: reproducible NIST isotope, CODATA2022 and H/Na threshold acquisition
  retains source values, uncertainty strings and unresolved isotope shifts.
  `species/ReducedKinematics` resolves immutable additive reference states;
  separate formation/excitation terms survive cancellation in Q. Explicit
  model-error limits and numerical/speed controls are required. Compiled PIC
  masses and the reviewed legacy reference are unchanged; no newer threshold
  is silently substituted into frozen D04 channels.
- P03/P04: the same production kernel implements general two-body closure,
  parent-rest inertial photon recoil with photon momentum, and mass-weighted
  three-neutral microcanonical breakup with both identical-H occurrences.
  Inputs supply all states, energies, reference zeros and random variates.
  Nonzero parent momentum at the parent-rest entry is refused. This is not a
  full moving-parent LSO physical event provider. Independent review found
  and fixed loss of small excitation during large-formation cancellation,
  admission of zero directions by a wide tolerance, and tiny-polar-angle
  loss of representable transverse recoil. Original reproductions survive.
- P05: NIST ASD 5.12 source acquisition retains 614 levels, 910 lines, 56 missing
  A values, eight selected O(1D)/O(1S) components and four fine-structure
  descendants. The stored graph has five levels and twelve distinct multipole
  components. All six processed products reproduce byte for byte. This
  establishes source-record processing only: completeness of physical exits,
  isotope applicability, alignment, quenching and omission bounds remain
  unresolved. Total physical lifetimes and normalized branches stay null;
  no physical cascade transport is enabled. OH populations/relaxation and
  the literal U24 surface/chemistry requirement remain unresolved.
- P06: explicit signed prescribed-reservoir accounting is retained. A scalar
  Maxwellian/Voronov rate does not supply collision-conditioned incoming
  electrons or joint outgoing products; unsupported impact kinetics stay OFF.
- P07: `SodiumPartialCapture` authenticates the twelve native IAEA D119990
  4–25keV/u partial-channel nodes and returns immutable density/node/constants
  receipts. The experimental forward-COM products use the shared two-body
  solver, explicit target-rest frame and source-derived Q. Native initial Na
  state is unspecified; Na3s is the approved reference cohort, not measured
  initial-state metadata. Physical state qualification remains false. No
  total-capture substitution, interpolation/extrapolation or physical PIC CEX
  provider was installed. Mass evidence is independently re-extracted from
  authenticated raw bytes; edited derived masses are rejected.
- P08: new `IonCollector` and `ion_surface_boundary=unsupported|ideal-collector`
  retain a strict default and require explicit physical-time tracing plus
  source/state/energy resolver installation. The normal sphere dispatcher and
  static/moving terrain callbacks dispatch ions before neutral sticking.
  The adapter uses actual local particle weight/correction and particle ET;
  geometry epoch, normal, surface velocity and inverse state transforms are
  retained separately. The collector records material/charge/momentum and
  kinetic/formation/excitation transfers; it neither deletes particles nor
  invents heat, neutralization, implantation, sputtering or release. The mover
  owns deletion. Neutral species retain the original callback/RNG path.
  Allocation is prepared before ledger publication; unsupported/missing state
  data terminate explicitly rather than losing an ion. Real moving-SPICE
  terrain ion qualification is still a separate requirement.
- P09: `HeliumSilicaChargeSurvival` preserves all seven native Tellez 2012
  Table I coefficients/uncertainties, regimes and gaps. It evaluates only the
  explicitly supplied published reciprocal-velocity coordinate; conversion
  from arbitrary physical velocity vectors returns Unavailable. Conditional
  scattered-ion survival is not incident retention or a lunar material law.

Production changes are confined to Moon modules: `MoonInput.{h,cpp}`,
`MoonSurfaceEnvironment.{h,cpp}`, `MoonSurfaceTemperature.cpp`, `Moon.cpp`,
`Moon.h`, `main_lib.cpp`, `makefile` and the four new species components.
Thin runners/fixtures/source tools occupy their own `test/reduced_kinematics`,
`ion_collector`, `sodium_partial_capture`, `helium_silica`,
`physical_state_sources` and `state_cascades` directories. The actual PIC
fixture adds thirteen normal-sphere collector controls per native layout.
Detailed comments describe units, frames, failure behavior and ownership.
The component inventory, application/species/test READMEs, RUN.md, input
syntax, gap registry and capability/approval registers are updated. No
permanent U/I ID or acceptance threshold changed.

### Actual results and preserved earlier attempts

Counts below describe different, sometimes nested scopes; do not sum reruns
or infer full owner PASS from a supporting component.

| Scope | PASS | FAIL | ERROR | SKIPPED |
|---|---:|---:|---:|---:|
| Reduced recoil controls |264|0|0|0|
| Collector kernel/actual parser controls |62|0|0|0|
| Na native partial-rate/forward-product controls |162|0|0|0|
| He/silica source/kernel controls |15|0|0|0|
| O I offline source controls |15|0|0|0|
| Corrected global U01–U26 |15|0|0|11|
| Actual linked I01–I10 |10|0|0|0|
| Eight-species production-rate callback rows |32|0|0|0|
| Actual PIC packet/normal-sphere collector |945|0|0|0|
| Actual Moon-force tracers |89|0|0|0|
| Constant-history reaction controller |73|0|0|0|
| Common transport-clock hooks |80|0|0|0|
| Full normal startup, both orders at 8 MPI ranks |32|0|0|0|
| Formal I21–I24 |0|0|0|4|

The kernel ensemble used 40,000 fixed-seed production draws and 17 statistics
with the frozen simultaneous 99% bound; its largest scaled deterministic
discrepancy is 1.19e-19 versus 1.39e-17. Na's twelve-node independent chemical-Q
closure diagnostic is checked on the Q scale, not hidden by beam energy;
maximum relative discrepancy is 3.40e-16. At the first 4000 eV/u control node,
speed is 878568.526738 m/s and Q is 1.355338514898e-18 J. These are source/model
controls, not measured lunar predictions. Thirty-five independent He law
comparisons have zero reported relative numerical discrepancy. CSV/PNG/EPS
comparisons, seeds, source/compiler/executable identities and logs survive.

The original native invocation and first continuation report ERROR: isolated
installation template timestamp and retained absolute PTDIR/MYDIR paths
respectively. Candidate-only timestamp preservation/site relocation corrected
them, with receipts. The next invocation really ran all ten linked M0 gates
successfully but reported overall FAIL because U02 still expected direct
legacy sphere registration. Its source guard now checks the actual ion-aware
registration **and** unchanged neutral delegation, together with all original
target/frame/force/source requirements. Original FAIL/output are retained;
the complete corrected U-series passes its available criteria. All 698 native/
configuration source identities match the actual ten-PASS M0 producer.
`m0_reconciliation.json` explicitly composes this same-source evidence; no
single original failed invocation was retroactively changed to PASS. No
physics criterion was waived. Actual `make -j16 test_Moon` produced a zero-byte
diff with the reviewed I04 hash unchanged:
`ed875cefc604e4ae3580fb89b8c528146929b6b9be099b41b61697faee7f0080`.
Earlier source/parser/interpreter errors and deliberate missing/corrupt-source
controls are retained separately from final successful scopes.
I01 runs its manifest check in process and I04 uses the parent nightly make
stage. Their inherited generic per-case stdout/stderr/receipt/command fields
are not independent execution receipts. `linked_log_index.json` identifies
the actual combined campaign/make output, numeric reference receipt and diff.
The command index labels declarations separately from execution histories;
no transcript or standalone per-case execution was fabricated.

### Build, resources, evidence and reproduction

The disposable candidate is `/tmp/amps-m3-approved-native-3l1e8dzm/candidate`;
the active source checkout's core-mesh configuration was never reconfigured.
Candidate support libraries are copies, so builds cannot modify active-tree
libraries. Each actual build removes `build` before ordinary Config.pl/SPICE
and make. Compiler: GCC 9.4.0 C++17, production `-O3 -g`; native default 8 MPI ranks,
OMP 1, explicit 1 MPI rank controls and reviewed 4 MPI rank nightly. Native rate-only cases
do not initialize SPICE; real force/startup/baseline cases retain actual
kernel hashes and frame/epoch metadata. Nominal packet/collector fixtures
make no physical isotope, thermal or moving-LSO qualification claim.

The original native deadline was 18:00:39 UTC. Setup started 17:00:39 UTC; final
tests and normal restoration finished 17:47:57 UTC, **2837.98s** including earlier
attempts. Peak measured process-tree RSS across supervision was
**15,745,789,952B**, below 200 GB. Remaining reporting/preservation stays within
that same deadline. Normal startup exercised both orders at 8 MPI ranks through real
`amps_init` and two actual time steps. Final supervision returned0/PASS; no
model/build/test process remains active. The isolated candidate is retained
for inspection; no temporary source directory was added to the active tree.
Eleven disposable Python bytecode files under srcMoon/test were removed;
maintained Python scripts and tracked source files were preserved.

Main evidence:
`test_output/srcMoon/m3-approved-products-20261118T163142Z/`:
`completion_report.md/.json`, `native-followup/result.json`,
`native_followup_execution.json`, the original native campaign directories,
`native-followup/u-series/run_all_summary.json`, per-run manifests/logs and
`native-followup/restore/normal_build_evidence.json`. Separate retained source/
component roots:
`m3-reduced-kinematics-20261118T163050Z`,
`m3-helium-silica-final-20261118`,
`m3-oxygen-cascade-sources-20261118T170845Z`.

Exact commands and reproduction alternatives are in RUN.md, per-run command
receipts and the completion report. Representative global invocation in the
correctly installed disposable checkout:

```bash
/usr/bin/python3 -B srcMoon/test/test_runner.py \
  --phase m3-development --isolated-checkout --mpi-ranks 8 \
  --campaign-wall-seconds 3600 --output-dir test_output/srcMoon/m3-native-NEW
```

`--isolated-checkout` acknowledges an existing disposable checkout; it does
not copy/relocate one. This package's retained launch records use only the
remaining original time budget after setup and each failed attempt. The final
artifact-only continuation delegates to the maintained runners and records
its exact calls, same-source M0 reconciliation and original absolute deadline.

Key actual source identities are retained in each manifest: NIST isotope v2
response `15de684d4be01faa70b0d72a19ddd5051cf330cf3c246bbcf36e423451a36983`;
CODATA `77fb90e66c40db3e6eb16630bc9c88e4c7c8beddbe5e71be406f2f26e3f67e67`;
IAEA twelve-node CSV `2a887ad43497a7f9e2b3d7b6173ae1b2597cc893d897926762c181946c92db93`;
He/silica PDF `bf554c667bb01d03952d151f86e4e348b706f20630c8fe896eedf15f6fa30ed2`;
OI level export `d331651674508628d040ca931019db9ec3de5774b19d37cc2a71c11788cab049`;
OI line export `43de51535539bc3de51795e337e79903156feb285c0fb75d8091a2c2c2b95b16`;
D05 native coefficient table
`470ea3078c64159eece6a469314dde2356e2d6ccaf618f4eaecb46649d61a5cd`.
Canonical validation-data/raw files were not modified or promoted to READY.

### Preservation and unfinished requirements

All seven protected active-root build/site/ignore files are unchanged from
package start, including the already locally modified Makefile. Both frozen
physical/qualification JSON identities are unchanged. All thirteen shared-core
identities match the retained 843-PASS preservation audit; full CORE was not
rerun in this Moon-only package. Its historical strict-coverage FAIL and CME
timeout ERROR remain visible. No staged files, commit, push, author contact,
reference replacement, hold-out inspection or observational scoring occurred.

Full M3 implementation/qualification remains **EXPERIMENTAL /NOT QUALIFIED**:
complete state/threshold/internal-energy/isotope crosswalks and reviewed model
error budgets; moving-parent photon and trajectory-dependent exposure/provider
wiring; collision-conditioned electron products; physical CEX families; complete
cascades/quenching and molecular/water-surface laws; physical ion-regolith
fate; and common-time science-estimator/environment/restart qualification
remain open. M1/M2 formal predecessors and unchanged I21–I24 statistical/
numerical requirements remain binding. Acquired rates, supported mathematical
closures, source QA and successful mechanical fixtures supply none of that
missing evidence automatically. The next package must address these specific
provider/data/qualification gaps, rather than running this fixture longer or
declaring an ideal collector a lunar material law.

## 5.29 Shared chemical database proposal — 2026-11-18

The user requested consolidation of all acquired chemical data into a single
database reusable by applications beyond Moon. The proposal is maintained in
`doc/CHEMISTRY_DATABASE_PROPOSAL.md` and its machine-readable decision register
`doc/CHEMISTRY_DATABASE_DECISIONS.json`. The eight SDB labels are proposed
design choices, not permanent test IDs or approved implementation. Source
basis remains `74974c2b2fc1b4cf09a73c320064a22370c0ca94` plus pre-existing
uncommitted work. Two parallel read-only audits covered data/permissions and
shared/comet architecture. The parent inspected native source terms, actual
build/input conventions, existing Moon services and primary-source software/
provider documentation before writing the concrete design.

The recommendation is one immutable, relocatable JSON/CSV catalogue with
hashed permitted source objects, named profiles and an application-independent
C++17 API. Detached runtime exports use typed CSV and the existing AMPS block
grammar, with acyclic release/profile/export hash bindings; no new native JSON
parser is required. Optional SQLite search is derived tooling, not a required runtime
dependency. A proposed persistent external root is
`/data/vtenishe/amps_chemistry_data`; it was not created. Code/tools/tests would
live in `src/models/chemical_database`, `utility/chemical_database` and
`test/chemical_database`. Existing source objects remain in place until a
reviewed permission-compatible copy/verification/migration completes. No
canonical/raw data is moved or edited by this proposal.

Every source must receive an explicit inventory disposition. Source integrity,
interpretation, model policy, permissions and application qualification remain
separate. Alternative rates/spectra, unresolved states, partial/composite
channels, non-reactive absorption totals, numerical anomalies, missing
uncertainties and incomplete decay graphs remain explicit. The database never
manufactures missing distributions or promotes an import to a science PASS.

The data audit found zero hash mismatches among 178 canonical raw files,
11 processed products and six processing scripts; also 153 atomic packet
files, 44 molecular retention entries, 30 ion/surface packet files and 16 O I
responses. Counts describe the respective manifests, not independent science
experiments. Canonical D04/D05 share one physical photochemistry package:
17 catalogue records (15 reactive), 30 solar endpoint rates, 2305 cross-section
rows, 1460 derived branch rows, four thresholds and 406 Voronov fits/four
neutral selections. Source/path/hash/count evidence is retained in the dated
proposal audit; no new readiness classification is assigned.
The reproducible audit records 449 hash checks/nine count checks across 448
unique observed inputs, zero mismatches and 0.826s execution. Its initial
root-discovery guard ERROR happened before source-data reads, was corrected
by requiring Config.pl plus srcMoon to establish the root, and is retained
separately. This is documentary evidence, not a science-test count.

Native OPEN-ADAS terms explicitly cover insertion into a managed database and
redistribution with modelling code. Current official terms confirm this. The
proposal therefore keeps such entries metadata/link-only until the necessary
written permission exists; local/private managed storage is not an automatic
exception. LXCat and publisher-specific distribution restrictions also require
source-specific treatment. Downloadability and webpage software copyrights
are not assumed to grant numerical-data redistribution. No permission request
or external contact was sent. Full inventory coverage does not imply every
source can already be packaged numerically.

Existing generic reaction/collision containers are retained as execution
views. Shared pure loaders/evaluators are reused behind compatibility façades.
The plan preserves Hartley2 startup rate-cache freeze, source/completion
precedence and deliberate neutral sinks, while Moon keeps its conservative
physical-time packets and actual particle-number accounting. Hartley2's
number-rate/normalized-time executor is not interchangeable with Moon's.
Charged identifiers require structured terms rather than the existing comet
parser's `+` splitting. Application geometry, thermal/plasma contexts, local
weights/time, MPI ownership, movers and samplers remain application-owned.
The reviewed design also retains fractional continuum mechanisms separately
from integer particle events, compile/runtime topology fingerprints, owning
immutable storage for legacy non-owning views, const host query concurrency
and an explicitly unqualified device-execution scope.

Six implementation packages cover full consumer inventory, permission-aware
source import, shared runtime/schema, Moon compatibility, Hartley2 second-
consumer verification and release/other-consumer migration. Proposed build
wiring uses the existing normal Config.pl/Make path and unique common object
names, retaining the user's uncommitted root Makefile constraint. Required
shared-core changes, if implemented, must preserve all historical PASS case
membership. The proposal gives common test contracts and twelve descriptive
coverage groups, plus the 3600-second monitored development profile, 8 MPI/
OMP1 and 200GB limits. Required expensive preservation/science checks may run
in separate cluster packages, retaining their gates until actual completion.
No new acceptance thresholds or automatic reference changes are proposed.

Source/data roadmaps remain respectively:
`712e1152bf50af6d9ee1dad8f8a585809ca8b92cf42ac6fc7a009f6c0253ee17`
and `d10eca1c539f9cf82940236d9e1aa36af665d43e82ba475925950d582fbf49c4`.
The source's approved development/resource amendments control work order;
unchanged common qualification gates and native metadata control scientific
interpretation. The differing roadmap bytes and explicit diff are preserved;
neither roadmap is silently rewritten. M1/M2/M3 formal status is unchanged.

Evidence:
`test_output/srcMoon/shared-chemistry-database-proposal-20261118T190424Z/`,
including `data-audit/`, `baseline.json`, inspected-file/command/source
identities and completion/preservation records. Updated documentation consists
of the proposal/register, STATUS and relevant shared/Moon/comet README links.
Production source/configuration, data, frozen physical/protocol specifications,
reviewed I04 reference, root site/Makefile files and Git staging are unchanged.
No application builds or science tests ran: science result counts are not
assigned. Documentary hash/schema/link checks are not validation. No stage,
commit, push, runtime source activation, hold-out scoring, data permission
request or external message occurred.


## 5.30 Shared chemical database implementation and M3 exit work — 2026-11-18

### Delivered scope and provenance

The user authorized implementation and installation of one shared chemical
archive at `/data/vtenishe/CHEMICAL_DATABASE`, with detailed API documentation
and examples. The source basis is
`74974c2b2fc1b4cf09a73c320064a22370c0ca94` plus the pre-existing uncommitted
work. Three agents worked on the shared C++ implementation/native preservation,
verification/registry integration, and source rights/provenance. The parent
implemented the Python API, inventory publisher, installer and documentation.
No Git stage, commit or push was requested or performed.

Maintained code lives in `src/models/chemical_database/`,
`utility/chemical_database/` and `test/chemical_database/`. The normal common
build compiles six uniquely named shared objects through `src/general/makefile`.
Moon's Photochemistry/ElectronImpact/ChargeExchange interfaces alias those
implementations; compatibility translation units contain no second copy of
source formulas. The nine thin verification clients compile the actual shared
production sources using `srcMoon/test/chemical_database_build.py`. The generic
reaction network and Hartley2 physical source/profile/executor remain unchanged;
second-consumer profile migration is further work.

The installation contains 84 files/4,469,969 bytes, 31 permitted payload records
and 30 deduplicated content objects, 27 source groups,622 source records,
35 table descriptors and7 available numerical/native-ASCII table products.
The remaining 591 records are metadata-only: 37 ADAS permission-pending and 554
with unresolved rights for this managed-database operation. Presence in the
catalogue never implies permission, interpreted numerical availability, a
selected transport law or physical validation. D04/D05 legacy archives remain
in their original locations; their existing research use is not reclassified
as prohibited by this publication. No canonical/raw source bytes are edited,
moved or silently replaced.

The native IAEA/ExoMol/unmarked-NIST entries use their source-specific rights.
NIST ASD copyright is not treated as a blanket public-domain .gov grant.
ADAS managed-database restrictions remain explicit; no author/provider was
contacted. `RIGHTS.md`, exact acquisition receipts and reproduction scripts
are installed. Four IAEA CSV representations reproduce their source arrays
bit-for-bit, retaining eV/u,cm² and target-rest frame semantics. No fitted or
synthetic source values are introduced. Unresolved isotope/state/mass fields
remain null rather than borrowing a neighbouring state.

Installed identity:

```text
release_id: amps-chemistry-20261118-r1
release_manifest_sha256: 5eac627528856de56a4e6b515c31f124a860060e9c02ecb5f3ec440e006876b7
profile_id: catalogue-v1
profile_sha256: 9431815d24f6954eb2fbf8a464650e432a052668b781814ec09980af620bcea6
export_manifest_sha256: 838ad22ee19037de52d9c92d87f5e5ed0acdfa7d4b7fbbb7cec3718a6e548ed5
lock_file_sha256: fdf5311274871d7abe884b3d8aec163bfdc52e8c3295adb09e4c5cd42813d945
```

The profile provides catalogue/source/table access. It deliberately lacks a
physical scalar-rate/particle-transport selection. Typed capability/domain/
permission failures prevent promoting held data through a factory. Runtime
queries are immutable and owning; no SQLite, PIC, MPI, SPICE or native JSON
parser is needed by the portable C++17 reader. GPU/device execution is not
qualified. The SDK, Python API and CLI relocate independently of the AMPS
checkout. The installed README gives exact source/table, units, compiler,
Python and native-reproduction examples. Published files are read-only; the
convenience lock is installed last after all content/software checks.

### Verification and preserved attempts

The main evidence directory is
`test_output/chemical_database/development-20261118T204811Z/`. Exact build,
installation, SDK/example and query argv, logs, identities and counts are in
its receipts/results. Independent scopes, which must not be double-counted:

| Scope | Actual final outcome | Evidence |
| --- | --- | --- |
| Installed shared native/Python/data/relocated SDK |56 PASS,0 FAIL/ERROR/SKIPPED |`installed-api-verification/result.json` |
| Final pre-install shared suite |56 PASS,0 FAIL/ERROR/SKIPPED; repeated identity check |`final-api-verification/result.json` |
| Nine migrated Moon scalar/adapter suites |558 PASS |`test_output/srcMoon/chemical-database-verification-20261118T204349Z/completion_report.json` |
| Inventory controls |13 PASS |`test_output/srcMoon/chemical-database-inventory-check-20261118T210741Z/result.json` |
| Acquisition controls |14 PASS |`test_output/srcMoon/chemical-database-acquisition-check-20261118T212452Z/result.json` |
| IAEA native reproduction controls |11 PASS |`test_output/srcMoon/chemical-database-native-iaea-20261118T213100Z/result.json` |
| Exact installed README examples |PASS; CLI,Python,SDK compile/query and native reproduction |`readme-examples/result.json` |
| Installed manifest/software integrity |PASS;86 nested data/index checks,27 software checks |`installation-receipt.json` |
| Linked M0 |I01–I10:10 PASS,0 FAIL/ERROR/SKIPPED |native package `moon-m0/baseline_summary.json` |
| Actual reviewed nightly |PASS;four MPI ranks,zero-byte diff |native package `moon-m0/I04/` |
| Four-neutral/eight-species coefficient dispatch |32 PASS;actual8MPI,OMP1 |native package `eight-species-callback/result.json` |
| Full CORE protected membership |PASS;843/843 historical PASS retained,plus1 new shared PASS;844 PASS/240 SKIP/1 existing FAIL/1 existing ERROR |native package `core-regression-final/result.json` |

Native evidence:
`test_output/srcMoon/chemical-database-core-preservation-20261118T210311Z/`.
All actual builds use isolated `/tmp` checkouts and the normal Config.pl/Make
path after removing `build`. Active root application/site/Makefile stay intact.
SPICE points to `/home/vtenishe/SPICE/cspice` and
`/home/vtenishe/SPICE/Kernels`. Compiler/MPI/configuration/executable/source and
external-data identities are retained per native manifest. The native package
has one original21:03:11–22:03:11UTC3600-second budget,8MPI/OMP1 for the
additional callback, the deliberate four-rank nightly baseline and200GB
aggregate RSS ceiling. This is the native verification package budget, not a
claim that the entire database design/implementation took less than an hour.

Initial real infrastructure failures remain visible: strict-compile diagnostics,
relative runner path, a malformed-lock raw exception, failed-subTest result
bookkeeping, and a padded native column-name issue were corrected and retested.
A held CSV with an empty native header is now an opaque metadata descriptor,
not a silently repaired active numerical table. The first native CORE discovery
ERROR identified missing registry coverage for the new component. The exact
registry/component/traceability entry was added. The resulting old harness
51-owner/1085-case pin failed; its exact expected totals were updated to
52/1086 and actual 13-case harness retry passed. Original reports remain intact;
final reconciliation must cite both the original full report and the actual
same-candidate retry. The first full CORE report also retained 97 reduced-SEP infrastructure ERROR
rows because the isolated snapshot lacked the untracked root task fixture.
The exact 15,578-byte existing root fixture was found byte-identical to the
protected prior snapshot, copied with its SHA-256 receipt, and the actual
97-row owner rerun passed 62 cases with 35 declared SKIPs in 42.4 s. No synthetic
sentinel, science code or reference was substituted. The final auditor cites
all three actual reports: full CORE, 13-case harness retry and 97-row reduced
retry. Its three self-controls pass; all 843 historical PASS IDs are present
and remain PASS, with no regression or newly failing case. Reconciled counts
are 844 PASS/240 SKIP/1 FAIL/1 ERROR. The existing strict-coverage failure and
full-volume CME 180-second refinement timeout remain open; this is a regression
preservation PASS, not an all-green CORE claim. Both private checkouts restore
ordinary SPICE Moon successfully. Setup, all failed/corrected attempts, builds,
tests and restores completed in 3097.47 s (51.6 min); final report-file writes
finished within 3360 s (56 min), below the original 3600-second limit. Sampled
simultaneous process-tree peak RSS was
9,787,752,448 bytes against 200,000,000,000 bytes. Active root configuration
and 632 native source files retain their declared identities.
No scientific acceptance threshold or baseline bytes
were altered. The reviewed I04 reference remains
`ed875cefc604e4ae3580fb89b8c528146929b6b9be099b41b61697faee7f0080`.

The new portable CORE owner is `CORE-CHEMICAL-DATABASE-001`, case
`CHEMICAL-DATABASE-PROTOCOL-001`. It wraps 52 site-independent protocol controls;
its single wrapper result is not 52 additional independent CORE cases.
Existing U01–U26/I01–I46 identities remain intact. Real-data tests add four
controls when a real lock is explicitly supplied. A source QA/API PASS remains
software verification and never substitutes for linked physics validation.

Useful exact commands from the checkout:

```bash
/data/vtenishe/CHEMICAL_DATABASE/bin/chemdb --lock /data/vtenishe/CHEMICAL_DATABASE/database-lock.json verify
python3 -B test/chemical_database/run_all_tests.py --output-dir /tmp/chemdb-verification-NEW --database-lock /data/vtenishe/CHEMICAL_DATABASE/database-lock.json
```

The output directory must be new. Native full CORE commands, temporary paths,
machine overlays, build/removal/restore commands and resource records are in
that package's `execution.json`, checkpoint and final completion report; these
commands preserve deliberate baseline exceptions rather than forcing every
owner to use the same MPI rank count.

### What is needed to complete the current milestone

The current implementation milestone remains **M3 chemistry and budget
closure**. M0 regression is preserved; formal M1/M2/M3 qualification remains
open. Database delivery closes the common storage/API package, not the physical
provider or predecessor gates. The following work remains:

1. Freeze usable, permission-compatible numerical application profiles and
   complete their species/isotope/state/threshold and formation/internal-energy
   crosswalks. The current catalogue-only profile cannot select numerical
   transport. Resolve molecular threshold/tail interpretation and O-state
   cascades/omission or quenching bounds, plus the U24 H2O/OH surface dependency.
2. Install production physical event providers for photon ionization/
   dissociation, electron impact and charge exchange. Reuse the approved
   P01–P09 reduced closures only within their declared source/state/energy and
   model-error domains. The existing parent-rest photon mathematics needs
   moving-parent/frame-aware coupling. A partial Na capture kernel is not a
   complete total-capture family; unsupported channels must not be discarded
   and renormalized.
   Approved reduced closures remain legitimate explicitly experimental choices;
   this does not reopen their approval or require a measured joint distribution
   for every reduced-model case. Their source/state/energy mappings, declared
   approximation error and actual production coupling still need qualification.
3. Couple actual trajectories to reaction exposure through the production
   illumination/eclipse/terrain/radius/plasma services. Declared affine
   exposure mathematics and the constant-history controller are supporting
   checks, not that physical provider. Preserve collision-conditioned bath
   selection, event ordering, unused daughter intervals and AMR/MPI handoff.
4. Qualify ion-surface behavior for surface-coupled cases. The strict boundary
   and approved ideal collector support bounded verification cases; neither
   supplies lunar neutralization/trapping/re-emission. The He/silica scattered
   charge-survival operator has a distinct limited laboratory scope. Measured
   plasma cases also require the missing D06 electron/ion forcing products;
   explicitly prescribed analytical baths are allowed only in declared tests.
5. Complete the linked I21–I24 executor and physical ledgers: represented
   number,nuclei,charge,reference mass,momentum,energy and signed bath/radiation
   exchanges. Run lifetime/branch checks and time-step/particle-weight/MPI
   controls at the same physical final time with both force-aware mover orders.
6. Freeze remaining qualification choices before execution. The current
   `M3_QUALIFICATION_PROTOCOL.json` explicitly still marks additional
   statistical/numerical design as proposed: numerical error budget is null,
   rare-channel precision/power design and original-parent stream independence
   need completion, and I22's total-removal versus partial-rate estimator needs
   an explicit definition. This file was not silently promoted by P01–P09 or
   by database approval.
7. Satisfy same-campaign formal predecessors I11–I13 and I14–I20 in the declared
   order. Approved early M3 implementation/local verification continues, but
   no database/core regression result bypasses M1/M2 closure. Missing
   authoritative data still produces SKIPPED / NOT QUALIFIED. Observational
   validation belongs to later layers; it is not manufactured as an M3 result.

A concrete next implementation slice is the normal-path Na photoionization
case with a complete Na+/electron/radiation ledger and actual trajectory
exposure, using the approved bounded verification boundary. That slice still
requires an explicit state/mass/energy mapping and reviewed qualification
budget. Its successful software fixture would then be extended to complete
families and linked formal cases; it would not close M3 alone. No missing
physical law or numerical acceptance value has been filled with an assumption.

## 5.31 M3 physical-provider infrastructure, strict source profiles and independent audits — 2026-11-18/19

**Disposition:** M3 implementation and supporting verification advanced; formal
M3 remains **SKIPPED / NOT QUALIFIED**. The approved early-development exception
does not close formal M1/M2 predecessors. No observational result, complete
physical Na run, full I21–I24 result or milestone completion is claimed here.
This section preserves the earlier source bases, failed attempts and explicit
scope of each numerical/software/native result.

### Physical source profiles and actual remaining data

`species/PhysicalStateProfile.{h,cpp}` authenticates exact CSV/native-source
hashes and resolves explicit isotope/nucleus/integer-charge/state entries through
the SAME production `ReducedKinematics` resolver. It retains independently
sourced physical rest masses separately from additive reference inertia; an
unknown physical mass is never replaced by reference mass. Formation and
excitation are independent sourced terms, not a fitted Q or an implied ground
state. A loaded profile is not a qualified reaction family or permission grant.

The real reproducible partial profile has **seven available reference states**:
`1H:ground`, `4He:ground`, `16O:ground`, `20Ne:ground`, `23Na:ground`,
`40Ar:ground`, and `electron:free`. It has **fourteen explicitly unavailable
states**: six singly ionized isotope entries, H2(16)O,16OH,1H2, O(1D),O(1S),
unlabelled O product, H2(16)O+ and16OH+. Unavailable entries retain reasons and
blank physical energies/masses, rather than fabricated zeros. The three source
roles are frozen NIST isotope masses, CODATA2022 constants and the approved P02
reference-coordinate definition; seven available states are not seven
complete physics datasets or seven qualified channels. The maintained source
preparation and its real-source83-control suite are under
`utility/chemical_database/state_profiles/`.

Additional actual public-source research resolved the previously incomplete
Na ASD uncertainty metadata. The retained official ASD response explicitly gives
5.13907696eV with supplied uncertainty0.00000025eV and ground configuration
[Ne].3s,2S1/2. Frozen D04 gives5.13907097104346eV. The arithmetic difference is
**5.98895654e-6eV**, or23.95582616 times the supplied ASD uncertainty; this is
not a sigma/confidence claim. The source threshold/state/energy-zero crosswalk
remains unresolved. No printed-digit rounding interval, isotope/hyperfine
mapping or uncertainty is invented; no source value or threshold is replaced,
clipped, widened or renormalized. The source/Q admission now refuses a physical
threshold mismatch independently of a caller's general numerical tolerance.
The older365-control composer result remains its historical basis; the later
additional source/Q refusal control yields366 controls in the parent's updated
suite and must be cited separately.

Fourteen exact retained research response SHA-256/MD5 pairs were checked;
seven ATcT requests returned native403 challenge responses. These are retained
HTTP outcomes, not acquired numerical thermochemical tables. No challenge was
bypassed and no provider contacted. The public NIST H response adds an explicit
uncertainty record, but no H/proton state profile is auto-selected. The NIST2025
H2(16)O rovibrational/partition-function paper is a documented lead, not an
acquired formation-energy/joint-daughter table. D04 still supplies15 reactive
channels but only4 atomic thresholds; all11 molecular thresholds remain absent.
Original positive He/Ar below-listed-threshold samples and OH_exp tail/state
ambiguities are preserved. Scalar electron/CEX rates or partial finite-domain
capture tables do not silently establish missing collision-conditioned complete
product/state/angular laws or unsupported tails.

Research, native responses, original URLs/outcomes, source audit, reproducible
extraction and limitations are retained at
`test_output/srcMoon/m3-physical-completion-20261118T223145Z/state-profiles/research/`.
The source-bound profile remains external original research/preparation output.
The installed `/data/vtenishe/CHEMICAL_DATABASE` catalogue release and its
immutable lock/profile are not changed. Held sources are not copied into a new
managed numerical archive or cleared by changing a readiness flag. Approved
reduced source-mean/isotropic or other bounded P01–P09 closures remain legitimate
EXPERIMENTAL choices within their approved source/state/energy/error domains;
a measured full joint PDF is not retroactively required for every approved
reduced case. The actual compatible source inputs and their domain evidence
are still required.

### Implemented independent software and production coupling

The bounded moving-parent Na composition service binds the complete D04 single
reactive channel to explicit parent/ion/electron states, reference/physical
masses, source energy convention, exact SI eV definition, frame/velocity offset,
photon direction, point visibility and caller numerical/model-error controls.
It produces a conservative tracked-ion plus signed electron-minus-photon bath
receipt at the actual event frame/epoch. Physical activation still refuses
unavailable Na+ data. Molecular, electron-impact and full CEX physical provider
families remain separately incomplete; the Na slice does not close them.

The trajectory exposure service now operates on the actual quadratic position
extension and production point-Sun geometry. Independent review found and
retained three geometry failures: merged nearby roots understated shadow error,
far-Sun subtraction cancelled silhouette normals, and binary64 midpoint
projection changed visibility outside ideal-arc root brackets. The current
implementation uses outward interval coefficients and Horner complement-sign
proof, the actual unresolved-bracket union, explicit geometry work bounds, and
full-cell uncertainty whenever the callback point cannot be proved in the same
strict sign cell as the exact midpoint. No root-count-times-epsilon estimate or
chosen position epsilon remains. The total declared error bound must hold or
calculation refuses. The37-control suite includes all counterexamples. Unknown
BVH/intersection/frame-transform errors remain separate geometry gates.

The variable clock preserves both compensated binary64 components through
arbitrary-size exact-rational location/comparison and a source-matching opaque
update receipt. Repeated7e-6s^-1 weak hazards, outer intervals, partitioning,
near-endpoint and full exponent-span cases no longer fail because of lossy
long-double projection. Existing constant/scalar APIs remain preserved;
`ExactPendingDepth` remains an explicit legacy projection check. Strong cache
identity compares the full clock. A rigorous reported low-component storage
rounding bound persists with the particle/accepted-event ledger. Supporting
counts are32 compensated controls,634 existing piecewise controls and25 bridge
controls; these nested/repeated scopes are not extra permanent U/I IDs.

Normal pre-buffer registration, semantic source/daughter genealogy/clock
storage, rollback-safe actual child initialization, explicit source-SI charge
policy and the conditional normal Na owner have been implemented. Its separate
pre-stage native package now passes20 normal-startup,945 legacy packet,56
variable-exposure and84 SI-semantic checks; restoration and final source-only
reporting retain their own receipts. Latest environment-stage production code
also has a fresh M0 and154 native wiring checks. These source bases remain
explicit; earlier results are not silently relabelled. The maintained
source compile recipe is now
`utility/chemical_database/state_profiles/compile_recipe/`: it authenticates
real available profile entries, stages an application-owned species include,
normal compiler/kernel companions with the selected basename, and independently
compiles the actual numeric charge-token expression. It preserves global
ElectronCharge/eV2J and canonical `input/species.input`. Missing Na+ refuses
before output. Real source constants and additive masses do not automatically
qualify a physical threshold crosswalk or model-error budget.

The independent event verifier parses authenticated source CSVs and native
hexadecimal moments directly to Fraction. It recomputes source additive masses,
formation Q, exact integer nuclei/charge, represented stoichiometry, signed
photon/electron baths, proper-frame transformations and actual committed initial
PIC child state/position/velocity/weight/epoch/interval values. It never accepts
a native `verified` Boolean. A valid empty event-file set has one explicit
empty-set identity criterion and zero evaluated physical events; the owning
audit must reconcile all original-parent histories, and unresolved source data
still skip before this identity.

The independent homogeneous timing auditor uses real frozen D04 endpoint
scalars, RFC2104 source-key HMAC and high-precision Decimal first-hit references.
It preserves actual original-parent membership, censoring, matched1/8 layouts,
three strict refinements and approved0.002/0.003 numerical allocations. Variable
trajectory cases cannot enter a constant-hazard oracle. Complete native include
graphs, source/profile/model identities, normal producer ABI, actual child logs,
failed attempts and original resource deadlines remain bound. Receipt schema,
rate-bound names, invalid-numeric classification, duplicate include graphs,
macro CALL validation and attempted-execution manifest propagation were reviewed
and corrected with retained regression controls. None of those controls is a
mocked successful physical cohort.

### Actual completed evidence and source-basis limits

| Scope | Actual evaluated outcome | Evidence |
|---|---|---|
| Initial source profile parser/source/arithmetic |83 PASS |`m3-physical-completion-20261118T223145Z/state-profiles/attempt3/result.json` |
| Corrected trajectory local geometry |37 PASS |same root `independent-trajectory-review/position-projection/maintained-fix/result.json` |
| Latest independent ledger software |21 PASS; retained labelled arithmetic reconstruction81 PASS |same root `final-independent-support-registration/sodium-ledger/` |
| Latest independent timing/provenance/executor software |37 PASS, including labelled generated arithmetic CSV/PNG/EPS |same root `independent-sodium-numerical/portable-clock-source-frozen/` |
| Maintained real-source compile recipe |19 PASS including opt-in homogeneous zero-acceleration compiler/parser binding; earlier12-control basis retained |same root `final-reference-recipe-registration/source-compile-recipe/result.json` |
| Registered combined supporting invocation |3 groups PASS/0 nonpasses,54 child criteria,28.70s |same root `final-independent-support-registration/result.json` |
| Later source-frozen registered invocation |8 groups PASS/0 nonpasses,618 child criteria,114.83s; profile-count protocol separately83 PASS |same root `final-source-frozen-supporting/result.json` and `profile-summary-protocol/result.json` |
| First isolated native package |M0 10/10; reviewed I04 diff0 bytes; backend/force/constant-controller/variable-controller/common-clock/startup945/89/73/52/80/32 PASS |`test_output/srcMoon/m3-physical-native-20261118T225229Z/completion_report.json` |
| First isolated full CORE preservation |844/844 protected PASS preserved; actual844 PASS/240 SKIP/1 preexisting FAIL/1 preexisting ERROR |same native report and exact case comparison |
| Pre-stage semantic/native follow-up |M0 ten/zero diff; normal20,legacy945,variable56,SI84 PASS; ordinary restore PASS; real source1/8MPI both SKIPPED77 |`test_output/srcMoon/m3-phase2-native-20261119T002618Z/completion_report.json` |
| Latest environment-stage extraction |Fresh M0 ten/zero diff;154 native frozen/evolving wiring checks PASS in111.34s |`test_output/srcMoon/m3-coefficient-stage-20261119T003919Z/completion.json` and `coefficient-stage/result.json` |

The first native package completed in3443.73s within its original one-hour
22:52:29–23:52:29UTC budget, with15.75GB sampled aggregate peak RSS under the
200GB ceiling. Isolated normal builds and original scientific references were
preserved. The two existing CORE nonpasses remain explicit: strict component
coverage FAIL and the180-second full-volume refinement timeout ERROR. This is
regression preservation, not an all-green full CORE claim. That package's M0
source basis preceded the final Moon-only projection fix; its actual fresh
all-species geometry/backend tests and37-control private check cover the later
geometry basis. Earlier M0 results were not relabelled as new-source results.

The separate phase-two native package finished in3272.40s within its original
00:26:18–01:26:18UTC budget, with15,261,839,360B sampled peak RSS. Its709
private production/build inputs remain unchanged, and611 shared CORE inputs
match the earlier844-PASS basis without claiming a redundant CORE rerun. The
two later root production deltas are exactly Moon.h/main_lib.cpp for the
separately tested environment-stage extraction. Actual source-preflight1/8MPI
returns77/SKIPPED for missing23Na+; this is not physical event execution. The
STATE-versus-position SPICE oracle ERROR, immutable-runtime test setup ERROR
and uncompiled-document freeze guard ERROR are retained with their mechanisms
and corrected source scopes. No physical acceptance criterion changed.
A further Moon-only environment-stage extraction package has original
00:39:19–01:39:19UTC budget and retains earlier attempts. The extracted complete
SPICE/orbit/thermal block is byte-identical, SHA-256
`e4665a6ee10b63972fc00afb97f757fb26dbccd27878a3ed9b1a8742cc40dc59`
(excluding the final newline). Independently reversing only that extraction
reconstructs both exact declared pre-stage mainlib/Moon-header hashes. Normal
order remains Begin→prepare environment→PIC::TimeStep→Commit; the finite owner
uses Begin→same helper→ExecutionTrackDefault→Commit. Frozen/evolving ET formulas,
SI conversion, frames, body transforms and thermal/inventory side effects are
unchanged by extraction. Source equivalence and runtime proof remain distinct:
the later actual fresh M0 passes all ten gates with zero-byte I04 diff;154
state/force-wiring checks pass on first/second-order frozen1/8MPI and ordinary
legacy evolving1/8MPI. Phase work completed in2065.17s with13.83GB sampled peak.
The original sandbox-launch ERROR and initial U15 support deadline ERROR remain
preserved. The final full U15 run takes301.915s and passes ALL20 supporting
groups with0FAIL/ERROR. Literal U15 remains SKIPPED for missing physical gates.
The named600s supporting allowance remains inside the SAME original01:39:19
deadline. The intermediate correction copied an earlier parser fixture; its
source-compile ERROR is retained, while the complete stable tool passes19
controls and the full final rerun passes. MPI admission/source-only follow-up
has its own retained output; current-source preflight1/8MPI both returns77.
This is not a positive physical Na or chemistry qualification.

Latest eight selected components comprise625 evaluated local criteria, all
PASS:83 source-profile,366 photon composition,37 trajectory geometry,30
approved protocol/statistics/RNG,32 exact-clock,21 independent ledger,37
timing/provenance and19 compile-recipe controls. These overlap full U15 and
earlier invocations and are not additional permanent U/I IDs. Final formal
I21–I24 are4SKIPPED,0PASS/FAIL/ERROR; their infrastructure controls pass6/6.
The consolidated machine report is
`test_output/srcMoon/m3-physical-completion-20261118T223145Z/work_package_report.json`.

### Governing contracts and reproducible execution

The two roadmaps remain different, with unchanged SHA-256 identities:

- source `srcMoon/development_and_validation.md`:
  `712e1152bf50af6d9ee1dad8f8a585809ca8b92cf42ac6fc7a009f6c0253ee17`;
- canonical `/data/vtenishe/moon_validation_data/development_and_validation.md`:
  `d10eca1c539f9cf82940236d9e1aa36af665d43e82ba475925950d582fbf49c4`.

Their actual difference is retained as
`test_output/srcMoon/m3-physical-completion-20261118T223145Z/roadmaps.diff`.
The source additions record approved memory/mesh/iteration/pilot resource
profiles, early M2/M3 implementation and partial Ar implementation notes. Those
approved source amendments govern development/resource scheduling. Original
scientific gates, literal IDs, native-source interpretation, calibration versus
hold-out rules and predecessor order continue to govern formal qualification.
The original model/protocol JSON remains unchanged; the separately approved
qualification addendum is pinned to
`551bb63b5de5be5766137bafb456aa140ae70fba5c38b4d9d8e9f70e74c15dbd`.

Exact commands, compiler/MPI/configuration/generated-header/executable identities,
source revisions, seeds, external file hashes, stdout/stderr, resource traces
and CSV/PNG/EPS are in each component and native manifest. Maintained commands
are in `RUN.md`, `test/README.md`, `test/sodium_cohort/README.md`,
`test/semantic_cohort/README.md` and `test/coefficient_stage/README.md`. Native
application builds use disposable checkouts and the normal remove-build,
Config.pl, make path. The active checkout is not reconfigured. A wrapper must
run from the same source checkout recorded by its actual compiled producer;
matching filenames or changing shell cwd does not establish that identity.

The active root Makefile/site/configuration files, original protocol/model-spec
and reviewed I04 reference still match all eleven protected package identities
in an independent current check. The I04 SHA remains
`ed875cefc604e4ae3580fb89b8c528146929b6b9be099b41b61697faee7f0080`.
Source revision remains74974c2b2fc1b4cf09a73c320064a22370c0ca94 plus recorded
working-tree changes. No canonical raw edit, automatic reference refresh,
threshold relaxation, commit or push is part of this package.

### Remaining completion conditions

Conditional executable software now exists, but actual source-admitted physical
cohorts still need compatible isotope/state/formation/threshold/physical-mass
crosswalks and independently frozen model/numerical error inputs. Current Na+
remains unavailable. Full molecular and other atomic photon families, electron
collision-conditioned products, complete CEX channels/state/angle support,
applicable state cascades and physical lunar regolith response remain separate
M3 work. Approved ideal collectors are bounded verification boundaries, not
lunar sticking/neutralization/implantation laws. U24 and applicable measured
R4.1/D06 dependencies remain explicit. Formal same-campaign I13–I20 predecessors,
complete I21–I24 family matrix, actual statistical/convergence/physical-budget
criteria and later observational validation are not waived by software PASS.

The later independent audit explicitly preserves the original I21 inverse-square
0.001 relative gate on BOTH actual source and accepted event rates, distinct
from0.002 timing allocation and caller ledger bounds. It authenticates consumed
bytes against ORIGINAL completed owner artifact maps and checks exact CSV/hex
relative-time projection plus explicit transport-origin absolute ET. The
chronology test's first31PASS/4ERROR is retained: systemPython3.8 has no
math.nextafter. Exact IEEE binary64 adjacent-value bit stepping resolved this
API gap without changing the preexisting two-adjacent-ULP epoch convention.
Independent signed-zero/subnormal/sign/binade/max-finite and2-versus3-ULP
controls pass in the final37-control run. The earlier combined54-child receipt
retains its earlier source basis; new/repeated counts are not added to claim
additional permanent U/I tests or a completed physical milestone.


## 5.32 Upstream integration and corrective CORE verification — November 19

Master integrates `60a156ddcc3deeaf48e2981b0beb14738245457e` with the preserved local source (eleven incoming commits). The completed corrective CORE registry has **844 PASS, 240 SKIP, 1 FAIL and 1 ERROR**. Exact historical/current-upstream PASS preservation audits pass; the known nonpasses remain visible. Moon M0 is **10/10 PASS**, with a zero-byte reviewed nightly diff, and **1533 targeted native checks PASS**. The U-series has 15 PASS/11 SKIPPED and zero FAIL/ERROR.

The original concurrent CORE attempt is **ERROR**: two CampaignBudget calls in
threads shared the same parent PID and adopted each other's children. Its
conservative outer reserve interrupted the merged lane before aggregate JSON
publication. Completed owner logs do not supply invented individual PASS rows.
The original failed deadline, logs and supporting short-tail retry are retained.
The completed pristine baseline is authentic: 843 PASS/240 SKIP/1 FAIL/1 ERROR.

A separately recorded corrective campaign uses one OS supervisor and one full
CORE candidate, with 24 build jobs on 36 logical cores. It retains its own ORIGINAL
3600-second deadline including setup/build/run/report/delivery, without resetting
that clock after failures. No owner timeout, fixture, reference or acceptance
threshold changes. Exact case-identity audits protect both historical 844 PASS
and current-upstream 843 PASS, and their union; incomplete/preflight evidence never
substitutes for actual runtime rows. Known incomplete-coverage FAIL and full-CME
180-second ERROR remain separate from preservation PASS.

The pristine CORE baseline and complete Moon runs retain actual 822dec07 source
revisions. Late upstream 60a156dd changes only an Earth source-contract Python
test; five focused checks pass. All native/CORE/Moon inputs and registry bytes
are unchanged by that addition. The corrective candidate executes on 60a156dd.
The 742 native-input hashes and the late-delta proof preserve these distinct bases.

Three-way merges retain upstream APIs and Moon's existing opt-in force-tracer
declarations in `src/pic/pic.h`, plus the local registration in
`MakefileTest/Table`. Its build uses `$(TESTCOMPILEJOBS)` as required by upstream.
The nightly planner now honors an explicit Table over an implicit cached
Makefile manifest, retaining explicit configured-tree priority and host checks.
Four new precedence controls, all 83 nightly unit tests, the official registry
check and six field-line helper cases pass. Three existing maintained chemistry
build helpers are explicitly unignored so future commits include required source.

Builds use fresh generated trees and normal Config.pl/GCC 9.4/MPI/SPICE. Native
manifests retain exact commands, input/configuration/header/executable hashes,
consumed data hashes, seeds, rank layout and stdout/stderr. Targeted Moon uses
one/eight MPI ranks; the reviewed nightly uses four MPI ranks/OMP 1 and produces the exact
reviewed data hash. I04's original generic scorer manifest is retained alongside
a SHA-bound actual-runtime supplement, rather than relabelling an unexecuted
command as the actual invocation. No observed random seed is invented.

Fast-forward delivery preserves unrelated local source, new files and all eight
site configuration files. No staging, user commit, push, reference replacement
or raw/installed-data edit occurs. Active generated build/executable and site
Makefile remain unchanged; remove build and configure/build normally before
using an active application executable after source delivery. Original source,
index and overlap backups remain available with hashes and recovery notes.

Formal M1/M2/M3 qualification remains open. Passing software fixtures and ideal
collectors supply no missing physical inputs or observational agreement. The
source/data roadmaps retain their different recorded hashes; approved source
resource/development amendments govern this work while original physical/data,
calibration/hold-out and predecessor gates remain binding.

Final evidence:
`test_output/srcMoon/upstream-core-merge-corrective-20261119T025503Z/completion_report.json` and `.md`.
The original error/native package is
`test_output/srcMoon/upstream-core-merge-20261119T015550Z/`;
its `latest-target/` holds late-delta and actual-application evidence. Actual
commands and hashes are in the individual native run manifests and package
contracts; preserve their source bases when reusing evidence.

## 5.33 Online physical inputs and immutable chemical snapshot — November 19

The user requested research of remaining physical inputs and their integration
into the existing shared database. Three research lanes covered atomic/thermal,
collision and surface inputs; the parent acquired molecular sources and owned
publication/integration. The source basis is
`8fbbd5e4b1a5dbaaa5329c30844c009a2b81d129`, plus the explicitly recorded
working changes. No commit, push or staging is performed by this package.

### Acquired source scope

* NASA CEA `thermo.inp`, pinned to official commit
  `2f79a647e85737742fb08a127da00d1ea1d069da`, contains 2111 native records.
  Exact projections contain 19 selected gas species and 56 coefficient intervals.
  Source molecular-weight, reference-enthalpy and Fortran-D conventions remain
  unchanged; eleven invalid condensed-phase domains in the full original remain
  explicit and are not repaired.
* ExoAtom supplies 423 Na I and 186 Na II source-relative states and definitions.
  Both definitions give the same mass and null ionization field. This supplies
  useful state identities without supplying a physical ion mass, absolute
  inter-charge formation-energy bridge or D04 threshold compatibility.
* ExoMol H2O DTU supplies two 7933-row absorption files at 423.15 and 573.15 K.
  Their 679/1448 negative and 575/574 zero samples are preserved. OH MYTHOS
  supplies 10606 rovibronic states and a 10000-row partition table. Actual and
  declared maximum state energy differ and remain separately documented.
* IAEA supplies 21 Na/He integral-ionization fit pages and 172 exact printed
  coefficient strings. Native labels, fit domains and uncertainties remain
  metadata. Failed full-precision export endpoints remain failures; no source
  formula bug or inconsistent domain is repaired or executed by this import.
* TU Wien's pinned regolith archive supplies four measured sputter-yield tables
  with 17 rows, separate from fourteen simulation tables with 2033 rows. Source
  errors, units and specimen-identifier disagreement remain native. Sputter
  yields do not supply an incoming-ion retention or charge-state kernel.

The maintained [physical-input disposition](../utility/chemical_database/PHYSICAL_INPUT_STATUS_20261119.md)
records source URLs, all six packet hashes and remaining requirements. NASA,
IAEA, ExoMol/ExoAtom and TU Wien have separate pinned rights evidence. CCC,
ADAS, SpuBase and other uncleared numerical sources remain metadata-only.
Specific HTTP403/404/500 failures are retained, not treated as source data or
proof that no authoritative source exists elsewhere. No canonical raw file is
written. Source signs, order, precision, gaps and missing fields are preserved.

### Publication, reproducibility and actual tests

The completed installation is
`/data/vtenishe/CHEMICAL_DATABASE/snapshots/amps-chemistry-20261119-r2`.
It contains **33 datasets, 721 source records, 80 table descriptors, 119
permitted payload records, 602 metadata-only records and 52 available tables**.
Relative to r1, that adds 99 records (88 admitted/11 held) and 45 available CSV
tables. These are catalogue counts, not independent measurements or qualified
reaction families. The installer verifies 195 files, 321 integrity members and
32 software/documentation identities. The exact lock pins are:

```text
release_id: amps-chemistry-20261119-r2
release_manifest_sha256: d03c1db383554b7f8aa16c08c3813790e595c04ace9838852871c2f1dc85bd6e
profile_id: catalogue-v1
profile_sha256: a99c3432d1f586dfee6923f0455ad2655b43293e3dfdafe5a562201054bc9fd4
export_manifest_sha256: 63e60b6ea0980308c1b335328c7c181694861da8a4ea05ae1b89640a5f374183
database-lock.json SHA256: 6c532f4ff5e37e4d766919c172b11848375f33d8296bb8bde5aa110b161ffd33
```

The maintained publisher now accepts repeatable explicit `--source-packet`
arguments; absence retains the original inventory. The offline strict importer
checks identities, paths, hashes, native-parent/recipe bindings, permissions,
rights evidence and false activation. Exact packet declarations are archived as
release members. The separate bounded HTTPS acquisition utility preserves
responses/statuses in a new packet and is never a runtime network fallback.
Native parsers and their documented units/failure behavior are shipped with
the portable Python package; SDK/runtime numerical formulas are unchanged.

The original root lock and all 84 original installed files are rehashed after
installation and remain byte-identical. The active local Makefile remains
byte-identical to the package-start copy. No active AMPS build/configuration or
source selection was changed. Only the standalone production C++17 data-access
service is compiled for these checks, with GCC9.4, optimization/debug/warnings
enabled; no application rebuild or MPI simulation is needed for this package.

Final evaluated counts are separate, overlapping software scopes:

| Scope | PASS | FAIL | ERROR | SKIPPED |
|---|---:|---:|---:|---:|
| Original r1 backward API/source/SDK checks, initial source basis | 140 | 0 | 0 | 0 |
| Final r2 API/parser/relocated SDK checks | 143 | 0 | 0 | 0 |
| Registered CORE case (122 source-independent child checks) | 1 | 0 | 0 | 0 |
| Original-source QA on actual installed r2 | 36 | 0 | 0 | 0 |

The last QA includes 45 actual C++ table queries, independent native-source
token/column/member checks, six archived packet pins, denial of all 602 held
source/26 held table operations, and 114 unchanged original inputs. A pre-install
check is reported as pre-install only; actual installed CLI/deep verification
and the post-install QA supply the deployment evidence. No new permanent U/I
ID is assigned. These repeated/nested counts are not summed into physical tests.

Earlier attempts remain preserved: the first real-release audit had one FAIL
because its expected simulation-role label differed from the native packet; the
test expectation was corrected to the source label, without data/model changes.
An intermediate candidate lacks archived packet manifests and remains SKIPPED.
The first global CORE attempt selected optional source-dependent classes and
therefore reported SKIPPED. Explicit pre-execution scope selection now keeps
those classes in the full runner and excludes them only from the declared
site-independent CORE scope; three administrative selection controls pass.
No missing original, failed physical result or acceptance threshold is hidden.

Exact publication argv, acquisition logs, compiler/executable/source hashes,
actual API/CORE/native-query stdout/stderr, original-preservation checks,
installation receipt and final report are under
`test_output/chemical_database/physical-input-integration-20261119T062501Z/`.
The maintained run commands are in
[the shared test README](../test/chemical_database/README.md),
[the actual-source audit instructions](../test/chemical_database/ACQUIRED_RELEASE_QA.md)
and the installed snapshot's README. To verify the installed archive:

```bash
/usr/bin/python3 -B /data/vtenishe/CHEMICAL_DATABASE/snapshots/amps-chemistry-20261119-r2/bin/chemdb \
  --lock /data/vtenishe/CHEMICAL_DATABASE/snapshots/amps-chemistry-20261119-r2/database-lock.json verify
```

### Physical and roadmap conditions that remain open

**M3 remains EXPERIMENTAL / NOT QUALIFIED.** No new physical model, reaction
network or source profile is activated. Compatible isotope/state/formation/
threshold/rest-mass mappings, supported molecular and collision-conditioned
products, complete applicable cascade information and physical lunar regolith
response remain incomplete. Approved reduced closures still require their real
source mappings and independent error bounds; this import adds no assumptions.
Formal same-campaign predecessors and I21–I24 qualification were not rerun or
closed by source-access tests. M1/M2 remain open as previously recorded.

The source roadmap SHA remains
`712e1152bf50af6d9ee1dad8f8a585809ca8b92cf42ac6fc7a009f6c0253ee17`;
the canonical data roadmap SHA remains
`d10eca1c539f9cf82940236d9e1aa36af665d43e82ba475925950d582fbf49c4`.
Approved source development/resource exceptions govern acquisition/development
order. Original scientific, data, calibration/hold-out and predecessor gates
continue to govern qualification. No roadmap is silently reconciled or relaxed.

## 5.34 Measurements and community approximation options — November 19

The user requested a new online search, prioritizing observations/measurements
and community-used approximations for remaining M3 processes. Three parallel
research lanes examined atomic/state/cascade, electron/CEX and material/surface
physics; the parent examined molecular products, published observation archives
and actual electron-forcing product metadata. This is research and source
inspection, not a new production implementation or a milestone closure.

The source basis remains `8fbbd5e4b1a5dbaaa5329c30844c009a2b81d129` plus recorded
working changes. The source and canonical data roadmaps retain their different
SHA256 values `712e1152bf50af6d9ee1dad8f8a585809ca8b92cf42ac6fc7a009f6c0253ee17`
and `d10eca1c539f9cf82940236d9e1aa36af665d43e82ba475925950d582fbf49c4`.
Approved source work-order/resource exceptions control development timing;
original scientific/data/predecessor/hold-out requirements remain binding.

### Findings that change the next development step

* **Na+ physical mass need not await another measurement.** The published
  CODATA atom/ion relation, NIST23Na mass, electron mass and ionization energy
  give nominal `3.8174499149007744e-26 kg`. The retained calculation propagates
  actual source uncertainties, with a first-order covariance-agnostic standard-
  uncertainty upper bound, rather than assuming independent inputs. An explicit
  state/constant/profile mapping and normal generated mass configuration remain
  implementation work.
* **Na energy conventions require a crosswalk.** Authenticated D04 supplies
  `2412.58 Å` / processed `5.13907097104346 eV`. Its rate HTML contains rate and
  mean excess, not a second independent `5.13907 eV` threshold. The earlier
  source-tree atomic note is corrected; no canonical data or immutable r2 file
  is edited. The proposed explicit reduced crosswalk keeps D04 hazard and source-
  effective photon energy while accounting for the `5.98895654e-6 eV` difference
  to the physical ionization energy. This is a model/source choice for review,
  not a tolerance relaxation or an already established absolute-rate error bound.
* **A concrete O I cascade source is acquired.** CHIANTI's actual native level
  and radiative arrays contain all six listed O(1D)/O(1S) exits and three ground
  fine-structure links in its lowest-five-level model. Its coherent A-set yields
  nominal lifetimes116.21666s /0.74695463s and green/UV photon ratio16.06733.
  Modern simultaneous Martian branching measurements give independent checks.
  The red-doublet observation versus matrix discrepancy and photon/energy/air-
  vacuum conventions remain explicit; no universal cascade-validation PASS is
  claimed and unrelated A-sets are not spliced.
* **Actual molecular product spectra are acquired.** The two Chang2021 publisher
  Source Data archives contain five XLSX originals, distinguishing measured
  energy/polarization spectra, fitted components, internal-state populations
  and mixed spectral inputs. These enable bounded benchmark cases; missing joint
  correlations, wavelength support and source units are not supplied by assumption.
  Their O(1S)+H2 channel is outside the frozen seven-channel D04 H2O family and is
  documented as a proposed versioned I22/source-registry extension.
* **Original surface measurements are acquired.** Jones2020's14 author objects
  (656794bytes, CC-BY4.0) contain Apollo14163/10084 TPD and separately inverted/
  simulated quantities. Exact specimen, coverage, heating, units, prefactor,
  sticking assumption and native error bounds must accompany a future U24
  benchmark. These are not14 independent measurements or a global OH sticking law.
* **There are useful observation archives.** JWST d203-506 Fig2 has2851 native
  spectral rows, with observations and synthetic alternatives in distinct columns.
  Poppe2022's acquired1600-row ARTEMIS grid distinguishes total/E-aligned/high-
  energy ion fluxes; its species-rate table is a model source. Neither product
  supplies missing uncertainty/composition/selection metadata automatically.
  Non-lunar OH observations are process evidence, not lunar validation cases.
* **Real electron products exist beyond merged D06.** Official THB_L2_ESA
  metadata supplies electron density/temperature, and THC_L2_MOM supplies onboard
  electron moments. Metadata are acquired. No new campaign/probe interval is
  selected, no CDF campaign is downloaded, and the canonical empty D06 package
  remains incomplete. A measured moment does not prove a Maxwellian distribution.

### Practical approximation choices, without hidden physical scope changes

The [full evidence/model review](test/M3_PHYSICS_LITERATURE_20261119.md) provides
primary links and exact next decisions. Existing fixed-mean/isotropic P03 remains
an approved named experimental closure. Published dipole/anisotropic alternatives
can assess the relevant observable; a full differential distribution is not
required to test neutral removal with terminal products and no feedback.

A new, separately named heavy-ion/reservoir variant is proposed for number and
specified pickup observables: create the ion at actual parent velocity/time and
transfer exact signed charge/momentum/energy differences to an unresolved bath.
That restricted model does not require outgoing-electron spectra. It still needs
source/state/Q and model-error/observable bounds; exact reservoir bookkeeping
does not prove recoil accuracy near escape, weak fields or surface interception.
It is a new model/specification choice, not automatic removal of current guards.

Opal/BED electron-energy marginals and angular closures have concrete published
precedents; their species-dependent widths, domains and inconsistencies remain
explicit. Proton–Na H(2s) measurements add0.5–2.5keV partial support. Actual NASA
proton–He differential/state-partial tables cover3–20keV/selected1–7degree lab
angles, with per-table uncertainties and genuine blank cells. OCR is an inspection
aid, not an accepted numerical source. Energy-corrected forward CEX is a possible
bounded model; copying unchanged nuclear velocities for nonzero Q is generally
not conservative. No H resonant law is silently transferred to Na or He.

CENA/ARTEMIS H-reflection measurements support empirical energy/angular models
for their specific H regime. They do not determine Na/He outcomes or permanent
retention as one minus reflected fraction. Published unity-sticking/fixed-energy
water-cycle and heterogeneous-site/efficient-diffusion models are explicitly
assumptions. The new Apollo data supply a material-specific alternative; a
prefactor from one source is not mixed with another source's inversion energy.

### Reproducibility, checks and unchanged qualification

All research originals, acquisition plans/receipts, source URLs/DOIs and SHA256,
successful and failed HTTP outcomes, scripts, lane reports and source ledgers are
under `test_output/srcMoon/m3-physics-literature-20261119T072531Z/`.
The maintained `utility/chemical_database/acquire_sources.py` executes the parent
plans. Native sources and analytical calculations are checked separately:

| Research scope | Evaluated result | Scientific meaning |
|---|---|---|
| Atomic source/independent arithmetic controls |13PASS,0FAIL/ERROR | Native-record and arithmetic QA |
| Collision retained-response identities |24PASS | Byte checks; actual failures/access pages remain recorded |
| Apollo author-file identities |14PASS,0FAIL/ERROR/SKIPPED | Producer size/MD5 and actual SHA checks |
| Apollo source-only packet admission |15allowed objects,0table descriptors | Importer admission; no numerical model or database installation |
| Parent molecular original checks |10PASS,1SKIPPED,0FAIL/ERROR | Failed author-version retrieval remains unavailable; original ZIP/XML/producer checks pass |
| Parent plasma metadata/author-file checks |5PASS,0FAIL/ERROR/SKIPPED | Native product definitions and author grids; no event forcing qualification |

These are overlapping acquisition/software scopes, not U/I physical PASS counts.
No AMPS build, MPI run, new permanent test ID, observation residual selection,
scientific threshold/reference change, canonical raw write, installed database
edit, staging, commit or push is performed. Newly acquired inputs are research
sources for a later reviewed immutable publication. Both r1/r2 remain unchanged.

Changed maintained documentation comprises this STATUS, srcMoon README,
`test/M3_PHYSICS_LITERATURE_20261119.md`, `test/TEST_GAP_PROPOSALS.md`, and the
shared database README/physical-status/atomic-source notes. Older findings retain
their date/source basis and the new follow-up explicitly supersedes acquisition
gaps where evidence now exists. Numerical implementation gates remain unchanged;
their1–3% tolerances are not a statement of physical source accuracy.

**M1/M2/M3 formal qualification remains open.** Actual source/profile admission,
reviewed new source/model choices, applicable production implementation, numerical/
statistical/budget checks and same-campaign predecessors remain to be completed.


## 5.35 — Shared physical-input/model integration, 2026-11-19

**Implemented, verified and installed as a shared source/software release.** The user approved
incorporating the section 5.34 sources and reduced-model specifications into the
application-independent chemical database. The work package is
`test_output/chemical_database/physical-model-integration-20261119T075920Z/`.
The preimplementation contract is
`utility/chemical_database/PHYSICAL_MODELS_INTEGRATION_CONTRACT.md`.

The release selects all six existing r2 packets plus seven new packets:
CHIANTI O I; Apollo sample TPD; two Chang water-product supplements; JWST OH;
Poppe ARTEMIS products; and the original NASA-CR114270 proton–helium benchmark.
Original tokens, source roles, gaps, units and record-specific permissions are
preserved. General explicit-input mathematical helpers, inert model descriptions
and shared Python/C++ access pass the declared software/source checks. The new
immutable r3 installation is published. No application provider, input or default lock is
retargeted; original r1/r2, canonical validation data and Makefile are protected.

No linked M3 or observational result is claimed. The full completion report,
exact commands, source/recipe hashes and retained attempts are in
`test_output/chemical_database/physical-model-integration-20261119T075920Z/WORK_PACKAGE_REPORT.md`.

### Delivered data, shared APIs and software

The installation is
`/data/vtenishe/CHEMICAL_DATABASE/snapshots/amps-chemistry-20261119-r3`, selected
through **its own** `database-lock.json`. Installed README section 9 gives
executable CLI/Python examples and shared C++ SDK guidance. The new services
have no PIC, MPI, SPICE, lunar-frame or lunar-species dependency.

The catalogue contains 41 source groups, 782 source records, 105 table descriptors,
180 admitted source payloads,602 metadata-only records and77 available tables.
Relative to r2, 61 source records and 25 CSV tables are added: CHIANTI five-level/
nine-link views;13 reversible Apollo TPD token views; five water-product
workbooks containing 102,816 explicit sparse cells; separate observed/model JWST
and ARTEMIS views; and the original NASA benchmark PDF/record. Native gaps,
negative values, malformed original headers, unknown units, cached values and
unevaluated formulae remain. No A-value or observational uncertainty is invented;
no OCR numerical cross-section table is admitted.

`model_definitions.py` supplies 12 explicitly selected, recursively immutable
model/scenario definitions through Python `models()`/`model(id)` and CLI
`list --kind models`. C++ reads the identical ordinary CSV registry through
`Catalogue::TableRows`; no runtime ABI/codec change is required. The publisher's
`--include-model-definitions` flag is opt-in. Old locks have no implicit models;
descriptions grant no physical capability or source permission.

`ReducedModels.h` adds six explicit-input mathematical helpers: singly ionized
mass balance, coherent radiative branching, Opal lower-energy marginal quantiles,
nonrelativistic two-body products, parent-velocity inheritance with signed bath
charge/momentum/energy, and first-order constant-temperature finite-inventory
thermal release. SI quantities, one declared event frame, source/state/Q and
all physical parameters are caller inputs. Typed refusals preserve outputs.
The installed standalone example takes all six thermal inputs explicitly.
Existing Runtime/DataIO and scalar-rate numerical implementations remain
byte-identical to HEAD. These helpers do not activate a Moon provider.

### Verification and identities

Source revision is `8fbbd5e4b1a5dbaaa5329c30844c009a2b81d129` plus recorded working
hashes. Standalone builds use GNU C++ 9.4.0, C++17, `-O2 -g -Wall -Wextra -Werror
-pedantic` and Python 3.8. No application/MPI/SPICE build is required here.

| Scope | PASS | FAIL | ERROR | SKIPPED | Evidence under this work package |
|---|---:|---:|---:|---:|---|
| Full shared suite |245|0|0|0| `shared-r3/result.json` |
| Registered CORE adapter |1|0|0|0| `core-chemistry/amps-core-tests.json` |
| CORE child methods, overlapping full suite |201|0|0|0| `core-chemistry/artifacts/CORE-CHEMICAL-DATABASE-001/shared-protocol/result.json` |
| Real-release source/API audit |19|0|0|0| `real-release-qa-candidate1/result.json` |
| Actual installed CLI/SDK commands |6|0|0|0| `installed-smoke/result.json` |

The audit compares all 74 available CSV tables and 1,841,608 exact cells through
original bytes, Python and C++; all 13 packet JSON and 10 declared recipes are
release members. All 602 held sources and 26 held tables remain denied in both
APIs with unchanged outputs. All 286 audited inputs remain unchanged. The full
suite also relocates a copy, runs its **own** CLI and compiles its **own** SDK.

Seven helper methods contain 38 native domain/transaction checks and 210
independent Decimal comparisons with the predeclared 1e-12 numerical budget.
Retained development failures exposed Opal conditioning and extreme two-body
mass/bulk arithmetic; the mechanisms are fixed without weakening a gate.
`test/chemical_database/REDUCED_MODEL_NUMERICS.md` preserves the binary-product
certificate proof and distinguishes mathematical error from physical accuracy.

```text
release 70a85232af7c9a24e669245e1f6d53fe70f26e37594cf75e4da08ea1a45f4f2e
profile 0f004933d7fdcb06b65c9250f4df51bf6680a505a539e6d0087356b1c3f8e729
export  5e93ac2132b25ba269b2336328e64e05c11886549ef9a7abf6019dfd44b89fd8
lock    52e382894da078b184559aacd73ecebc8dd8bc7f9fd61a5d08bfe33ee204c881
header  125dc21a5821bfbd4da07e91f955cbf82d7f06954ccbe312d3ce0294b27e3312
```

`publication-plan.json` records every packet path/SHA and exact builder argv;
`shared-suite-command.json` records the full-run argv/source selectors. The
full report gives repeatable commands with new output directories. Installer
preflight and actual copy pass; `installation.json` verifies 296 installed files
and the actual installed CLI. The preservation receipt checks 279 protected
r1/r2 files, including the original default lock. Makefile SHA remains
`0963a7bb50c03457e8a8c582649a1cf1eb4db8aca7d15574d7a0afa252a1f311`.
No canonical raw data, old snapshot, production input, Git index, commit or
remote is changed.

The different roadmaps remain unchanged: source SHA
`712e1152bf50af6d9ee1dad8f8a585809ca8b92cf42ac6fc7a009f6c0253ee17`,
data SHA `d10eca1c539f9cf82940236d9e1aa36af665d43e82ba475925950d582fbf49c4`.
Source development/resource amendments and the shared integration contract
govern additive code/publication; canonical scientific/data/predecessor rules
continue to govern qualification. No roadmap ID or threshold changes.

### Remaining physical/application gates

This release does not close M1/M2/M3. Direct ASD/ADAS and other managed-source
holds remain intact. A mass helper does not admit held Na threshold/mass data;
the source-effective photon/physical-state crosswalk and applicable source/rate
consistency require an admitted profile. Narrowband molecular measurements do
not determine every state/Q or solar-weighted joint daughter law. The reservoir
closure requires observable-specific error bounds and correct bath membership;
it predicts no electron spectrum/heating partition. CEX still requires supported
states, rate and differential/observable scope. Apollo material-specific TPD and
H-only CENA constraints do not define universal Na/He ion-regolith outcomes.
Descriptor-only scenarios stay non-executable, unsupported application channels
stay disabled, and linked numerical/statistical/budget/operator/predecessor
gates keep their previous status. No observational residuals are scored here.


## 5.36 — M3 source-bound implementation and verification, 2026-11-19

**Implementation advanced; formal M3 remains SKIPPED / NOT QUALIFIED.**
Work package: `test_output/srcMoon/m3-completion-20261119T091742Z/`.
Revision remains `8fbbd5e4b1a5dbaaa5329c30844c009a2b81d129`, with uncommitted
working-source changes; no commit, staging, push or reference update is made.
The full completion report and machine-readable results are retained in that
package. The existing P01–P09 approvals, immutable source hazards, U/I meanings,
numerical/statistical criteria and formal predecessors remain unchanged.

### Implementation and production paths

- The opt-in atomic profile producer authenticates the retained NIST isotope,
  CODATA, Na ground-state and ionization evidence. It derives physical Na+ mass
  and adds `23Na+:ground` without replacing the fourteen unavailable records or
  changing the old seven-state partial profile. Source-reported threshold
  uncertainty is not relabelled as a standard deviation. Its final profile SHA
  is `f739b4adeb50eae8feadc842c11e918ae7379341cc1395be7cd291fb5f730a26`;
  crosswalk SHA is `b45ce6f469a8244916522b8b933bd940031b41261c483f3916ef5b610b8c2217`.
- `ReducedKinematics` and `PhysicalReactionEvents` add an explicit
  `source-effective-photon` policy. The default source-equals-state guard is
  unchanged. `MoonSodiumCohort` authenticates the eight-column crosswalk and
  derives physical-minus-source threshold from actual state/source evidence.
  Frozen D04 rates and representative source photons are preserved. The
  physical formation gap, moving-parent transform, affine LSO moments and
  complete electron/photon ledger are accounted independently. The ledger
  verifier checks this policy without trusting the production receipt.
- The native owner supports explicit actual parent/ion block reference weights
  for the approved unequal-weight matrix. The new physical-case preparer
  requires every scientific control explicitly, checks the actual parser
  grammar, preserves source/include hashes and publishes no runnable input on
  missing configuration. It selects no numerical or model-error default.
- `build_m3_source_qualification.py` builds only in a disjoint disposable
  checkout, using normal Config/remove-build/make and actual producer receipts.
  The species recipe preserves inactive canonical symbols referenced by the
  current core and changes only selected physical masses/charges. Its explicit
  conservative-products option enables individual weight correction through
  normal input and proves the effective generated ABI after compilation.
  Default recipe layout is unchanged; no core compatibility patch is added.
- Shared `ReservoirProducts.h` supplies caller-evidenced collision-biased
  Maxwellian moments, inherited-velocity ionization and Q-corrected partial
  forward CEX accounting. It preserves separate isotope carriers, frame/energy
  conventions and complete signed bath moments. It supplies no electron PDF,
  bath identity, state evidence or physical-error allowance.
- `OxygenCascade` loads the selected real CHIANTI O I levels/links, computes
  lifetime/branch selection and finite-emitter photon recoil, with independent
  Cartesian dispersion and energy/momentum certification. Source isotope/state
  populations, quenching and omitted-network admission remain caller gates.

All new tests are globally registered under existing U15/U17 supporting owners;
new helpers hidden by the broad `*build*` ignore pattern have narrow source
exceptions. Heavy source builds remain explicit opt-in commands. Formal U/I
status is never promoted by a supporting PASS. Detailed contracts/commands
are in the atomic/compile-recipe, source_qualification_build, sodium_cohort,
physical_reaction_events, reservoir and oxygen_cascade READMEs.

### Evaluated results and limits

| Scope | PASS | FAIL / ERROR / SKIPPED | Evidence in work package |
|---|---:|---|---|
| Atomic profile/native registry | 78 | 0 / 0 / 0 | `registered-new-services/atomic-state-profile/` |
| Original partial-profile regression | 83 | 0 / 0 / 0 | `state_profiles/default-partial-regression/` |
| Strict/effective photon composer | 732 | 0 / 0 / 0 | `registered-new-services/physical-reaction-events/` |
| Independent sodium ledger controls | 24 | 0 / 0 / 0 | `registered-new-services/sodium-ledger/` |
| Reservoir products | 48 | 0 / 0 / 0 | `registered-new-services/reservoir-products/` |
| Oxygen cascade/recoil | 243 | 0 / 0 / 0 | `registered-new-services/oxygen-cascade/` |
| Source-build protocol | 10 | 0 / 0 / 0 | `registered-preparation-controls/source-build-protocol/` |
| Explicit physical-case preparation | 36 | 0 / 0 / 0 | `registered-preparation-controls/physical-case-preparation/` |
| Compile/species/conservative layout recipe | 27 | 0 / 0 / 0 | `state_profiles/compile-recipe-conservative-layout-20261119T1001/` |
| Registered shared CORE owner | 1 (249 child checks) | 0 / 0 / 0 | `core-chemical-database/amps-core-tests.json` |
| Clean source-bound conservative build | 1 | 0 / 0 / 0 | `source-build-conservative/` |
| Actual eight-rank source inspection | 1 | 0 / 0 / 0 | `native-source-eight-rank/` |

Overlapping local/shared counts are not unique tests and must not be summed.
The source-bound build takes 141.24 s, peak observed process-tree RSS
5,731,258,368 bytes; MPI source inspection takes 12.45 s with eight ranks and
OMP1. The latter inspects source states/masses and compiled bindings, launches
no particles and evaluates no transport/composer/numerical/statistical criterion.
Private source producer: `/tmp/amps-moon-m3-completion-20261119-gpbol3u9/AMPS`.
Root configuration/build and Makefile remain as found. All bounded campaigns
use at most 3600 s and a 200 GB ceiling. Ordinary nightly regression ran in a
second private checkout and PASSed: `nightly-default/result.json`, 280.33 s,
3,265,122,304 B peak observed RSS, zero-byte `test_Moon.diff`. The existing
reviewed four-rank I04 layout/reference was preserved. All thirteen selected
work-package scopes PASS; the full U/I and full CORE campaigns were not run.

Earlier build errors are preserved: missing normal generated `Makefile.test`
was fixed by normal private `Config.pl -install`; omitted inactive ENA species
symbols were fixed in the source-owned species recipe. No generated primary
fix, tolerance change or core patch was used. A subsequent canonical-layout
build passed, followed by the conservative-layout build above. Exploratory
kernel controls and corrected counterexamples retain their original outcomes
and independent oracle rationale in each lane report.

### Research, remaining implementation and required decision

Primary-source research now establishes isotope-specific 4He/20Ne/40Ar ground
thresholds with original uncertainty semantics. Their exact source ledgers,
D04 differences and remaining yield/source-state limitations are in
`atomic-family-research/`. These are research evidence, not admitted production
profiles. `cex-state-research/` preserves all twelve D119990 nodes (four QA
checks PASS), but recovers no initial Na state/isotope/citation/uncertainty.
Knoop's state-explicit alternative cannot be asserted to be D119990's source;
NASA's selected H(2s) angular measurements do not resolve residual He+ state.
No scans are silently digitized or missing product populations inferred.

The P02 numerical proposal was pending at the end of this earlier package.
The user subsequently approved relative inertia/acceleration limits 1e-8,
kinetic approximation limit 1e-4 and maximum speed 3e6 m/s; section 5.37 records
its implementation and actual controlled execution. The earlier package's
pending-review artifact remains historical evidence. The approval does not
certify source-mean/isotropic model accuracy or alter I21–I24 acceptance limits.

After this admission decision: execute and independently audit actual Na
cohorts, refinement/weight/layout/statistical matrices; complete admitted
multi-species/molecular state-energy mappings and production providers;
connect reservoir/cascade services only within supported bath/source domains;
close source-specific CEX/model-error and applicable surface gates. Formal
M1/M2 same-campaign predecessors remain required for I21–I24 qualification.
Neither this package nor the unchanged installed r1/r2/r3 database closes M3.

Source/data roadmaps remain different and unchanged: respectively
`712e1152bf50af6d9ee1dad8f8a585809ca8b92cf42ac6fc7a009f6c0253ee17` and
`d10eca1c539f9cf82940236d9e1aa36af665d43e82ba475925950d582fbf49c4`.
Source development/resource exceptions control implementation order; canonical
scientific/data/predecessor rules control qualification. No hold-out residual,
canonical raw data, installed release or reviewed reference is changed.


## 5.37 — Approved P02 bounds: actual Na commissioning and independent audits, 2026-11-19

**This controlled work package is complete. M3 remains EXPERIMENTAL / NOT
QUALIFIED.** The user approved relative reference-inertia and acceleration
errors 1e-8, relative nonrelativistic kinetic error 1e-4 and maximum inertial
particle speed 3e6 m/s for the controlled Na family. Additive authority:
`test/M3_P02_NUMERICAL_BUDGET.json`, with its companion explanatory Markdown.
The original profile-pinned decision JSON is unchanged. There is no pending
P02 numerical-budget question for this scope.

Package: `test_output/srcMoon/m3-p02-approved-20261119T175502Z/`; source revision
remains `8fbbd5e4b1a5dbaaa5329c30844c009a2b81d129` plus working edits.
`WORK_PACKAGE_REPORT.md`, `completion.json` and the package manifest give exact
commands, source/build/input/data identities, resources and historical failures.
No commit, staging, push, reviewed-reference update or canonical raw-data change.

### Applied inputs and minimal fixes

`examples/m3_sodium_p02_commissioning.{in,json}` and
`test/M3_P02_CONTROLLED_NA_PROTOCOL.json` freeze a quiet-source, 128-original-parent
commissioning experiment: frozen2009-01-24SPICE, actual LSO source at
(10000000,0,0)m, initial LSO velocity0, eight8192s intervals (horizon65536s),
equal reference weights7 and represented number2.625 per parent. ExactSI
149597870700m normalization and Earth occultation are explicit. The force-aware
second-order production dispatcher uses the existing compiled zero-acceleration
reference, not a duplicated force formula. The ideal collector is mathematical,
not lunar neutralization. This is source/kinematic commissioning, not a full
force-driven exosphere, observation or half-survival assertion.

1. The initial native adapter attempted to raise global dt through the core's
   minimum-reduction setter. Actual Na/Na+ steps stayed unequal; the normal
   common-clock guard correctly stopped it. The adapter now resets those
   reduction accumulators before applying the explicit reference dt and reads
   actual public table values back before creating particles. Core setters,
   normal application modes and common-clock criteria are unchanged.
2. The next run exposed a production rounded-frame inversion error. GDB's
   actual SPICE matrix and29.386km/s affine offset predict a spurious axisY
   momentum residual -3.660784e-37kgm/s; the actual refused residual is
   -3.660751e-37 against the unchanged3.137089e-37 bound. The composer now solves
   the supplied matrix instead of treating its transpose as an exact inverse,
   then checks the mass-weighted forward roundtrip. No SPICE coefficient, bath,
   source rate or acceptance tolerance is adjusted to fit a residual.
3. Both independent verifiers had applied a relative rotation comparison to a
   near-zero off-diagonal dot product. Their corrected exact COLUMN Gram
   `R^T R-I` and determinant checks use the original dimensionless UNIT scale,
   matching the declared production contract. The native0+1e-15 bound is
   unchanged. Initial FAIL artifacts remain; controls admit within-bound
   perturbations and refuse above-bound/exact-zero violations.

### Actual results

| Declared scope | Final result | Evidence relative to package |
|---|---|---|
| Production composer + independent inverse regression | 981 PASS, including all732 priorPASS | `rounded-affine-inverse-final/` |
| Approved/source mechanical-verifier controls | 30 PASS | `registered-final-controls/p02-independent-controls/` |
| Commissioning analytical/provenance controls | 20 PASS | `registered-final-controls/commissioning-audit-controls/` |
| Independent ledger controls | 27 PASS | `registered-final-controls/sodium-ledger/` |
| Explicit case preparation controls | 36 PASS | `registered-final-controls/physical-case-preparation/` |
| Fresh source-matched conservative build | PASS,248.884s | `source-build-affine-inverse/` |
| Actual eight-rank finite physical cohort | PASS,48.801s including thin-owner build | `commissioning-eight-rank-affine-fixed/` |
| Complete original-parent/source/rate/first-hit audit | 5 PASS | `independent-commissioning-unit-scale-corrected/` |
| Independent actual-event ledger reconstruction | 3655 PASS for43events | preceding audit's `ledger-00/` |
| Independent source/actual-event P02 cap audit | 1008 PASS for43events | `independent-p02-unit-scale-corrected/` |
| Clean canonical make-j test_Moon after production fix | PASS; zero-byte diff;309.802s | `nightly-after-affine/` |

All final rows have0FAIL/ERROR/SKIPPED. Counts overlap and are not unique tests.
The 128 original parents give85survivors and43Na→Na+ removals. Every original ID,
accepted event, actual initial committed ion, signed electron/photon ledger and
source clock is independently reconciled. The selected source's independent
rate is7.530630918928475e-6s^-1, with expected survival0.6104692; measured
survival85/128 is a diagnostic only, not a lifetime/statistical-acceptance PASS.
No probability, count, seed or source geometry was changed after this output.

`verified-data-repeat/` repeats the SAME128histories after prechecking exact
SPICE pins and D04 provenance/QA/output hashes; every external hash is unchanged
after the run. Parent CSV and all8exact-event files are byte-identical on this
host. The repeat is not another independent sample, cross-platform guarantee or
formal I45 result. Its separate independent audits are retained under
`repeat-independent-commissioning/` and `repeat-independent-p02/`.

Figures/CSV/EPS: `independent-commissioning-unit-scale-corrected/comparison.*`
and `ledger-00/comparison.*`; source-law timing and conservation diagnostics
contain no observations. P02 guards have scalar JSON/CSV/log artifacts. Run and
visualization/interpretation instructions are in `RUN.md` and the cohort
READMEs. New supporting runners are globally registered under U15 without a
new permanent ID. Whole-U/I statuses remain separate.

### Preservation, resources and remaining gates

All original failed runtime/verifier/control/debug attempts are retained with
mechanism analyses; corrected reruns use new directories. The final clean
producer is `/tmp/amps-moon-p02-20261119-_rpwgd3p/AMPS`; old producers remain
preserved. Core source, root Makefile/configuration/build and canonical sources
were not changed by this package. Builds remove their private build tree first,
use normal Config/make-j16 -O3/-g, and keep full logs. Physical runs use8MPI/OMP1
with200GB/600s phase guards inside the one-hour development package.

This closes the approved P02/first-controlled-Na commissioning work, not M3.
Full frozen timestep/order/weight/MPI/statistical matrices, final global
population/surface/escape budgets, admitted multi-species/molecular/CEX source
and product providers, and same-campaign M1/M2 predecessors remain. Statistical
qualification and refinement independence are NOT_EVALUATED by commissioning;
no observational agreement is claimed. Installed chemical-database releases
and all acceptance thresholds retain their previous status.

Roadmaps remain different and unchanged: source SHA
`712e1152bf50af6d9ee1dad8f8a585809ca8b92cf42ac6fc7a009f6c0253ee17`, data SHA
`d10eca1c539f9cf82940236d9e1aa36af665d43e82ba475925950d582fbf49c4`.
Source approved early-M3/resource amendments govern development order; unchanged
scientific/data/predecessor contracts govern qualification.


## 5.38 — Commit staging and upstream-core check, 2026-11-19

The user requested staging of the Moon/shared-chemistry implementation without
committing, plus a fresh upstream-core check. `git fetch --no-tags origin`
completed. Working HEAD is `8fbbd5e4b1a5dbaaa5329c30844c009a2b81d129`;
`origin/master` is `5941162c107411b8d9253d0b2243a61064b6d4b7`, four commits ahead.
The incoming core changes are:

- `src/general/specfunc.cpp`
- `src/meshAMR/meshAMRdef.h`
- `src/meshAMR/meshAMRinternalSurface_sphere.h`
- `src/models/exosphere/Exosphere.cpp`

These paths have no local edits and do not overlap the pending Moon or shared
chemical-database additions. The upstream core IS newer than this checkout.
The shared chemical-database component itself has no incoming changes. This
check does not silently merge/rebase/update core files or qualify the new base.
The preceding build/native/core-regression evidence belongs to the previous
HEAD and working-source identities; integration/regression on the newer core
is still required before calling that combination verified.

Source, tests, explicit examples, contracts and READMEs are selected for the
index. The local Makefile, generated root input decks/executables, run outputs,
archives, dependency trees and historical root handoff remain unstaged. Exact
selection/source hashes, upstream diff and index-preservation checks are in
`/tmp/amps-moon-staging-review-20261119-a555l7zo/`. No commit or push is made; M3
scientific status and all source/data/acceptance gates remain unchanged.


## 5.39 — Additional physical-input acquisition and shared database update

**Source acquisition/integration package COMPLETE; M3 remains NOT QUALIFIED.**
The requested database-first search, primary-publication follow-up and shared
database update are implemented. Source HEAD is
`26f1a65c58a8a70f1c6301b1a6c4d283a4558c94`, with exact working-source hashes in
the build/test manifests. Three agents handled independent molecular,
collision and surface/relaxation searches; the root handled atomic research,
shared admission/provenance, documentation, verification and installation.

Artifacts: `test_output/chemical_database/m3-physical-inputs-20261119T231417Z/`.
`work_package.json`, `publication-plan.json`, per-lane source ledgers/reports,
acquisition receipts, parser QA, full software logs and installation/preservation
receipts retain the complete source-to-release path. No source original,
acceptance criterion, application profile or observation was changed.

### New physical inputs and exact source roles

The new immutable installation is
`/data/vtenishe/CHEMICAL_DATABASE/snapshots/amps-chemistry-20261119-r4`.
Its own README section10 gives actual CLI/Python/C++ usage; select its lock
explicitly. The original root lock and r1/r2/r3 remain selected/available as
before. The new catalogue contains48datasets,885source records,267admitted
payloads,618metadata-only records,143table descriptors,115available tables
(112CSV) and the same12inert model definitions. Relative to r3, this adds
103source records and38accessible CSV tables; counts are source-access counts,
not independently validated laws.

| New inputs | Actually acquired numerical scope | Retained limitations |
|---|---|---|
| ExoMol OH+ MoLLIST, OH MoLLIST-OH, H2 ARLR and H2 RACPPK | Eight exact state/transition views:5046actual states and103947transitions, with native cm^-1/s^-1 units, isotopologue and dataset/version metadata | OH+,OH and RACPPK declared-state counts disagree with actual rows. Optional OH quantum labels are ambiguous and remain positional. No missing state, formation/Q crosswalk, full cascade or coefficient uncertainty is invented. |
| IAEA D107579/580/582–586 | Seven new native excited-Na AOCC composite pages and49exact printed coefficient rows | Capture+ionization is not pure CEX. Missing equation/upper domain and source-attributed10% model experience error remain explicit; no fit evaluator or daughter law is activated. |
| Zenodo5044539 Apollo15221 | Four native TPD/energy/site tables;9598rows | Experimental TPD remains separate from model components and inversions; signed values, specimen/preparation and distribution conventions are preserved. |
| Zenodo14588634 soil water | Six native XLSX files; eleven worksheets and19151represented cells | Simulated lunar soil, loess and montmorillonite remain distinct. Pressure/history, measured/calculated labels, units and unresolved figure numbering stay native. No prompt sticking/accommodation law is inferred. |
| Zhao2021 OH(A)+H2 | Three publisher ZIPs, thirteen long-cell views and35976retained records including138blank-row markers | Measured and quantum-calculated products, source energies, Ntot=0, duplicate headings and missing markers remain distinct. No Maxwellian total, missing elastic/inelastic outcome or complete joint product law is inferred. |

ExoMol's CC BY-SA4.0, IAEA's CC BY4.0, each Zenodo deposit's CC BY4.0 and
the quenching article's own CC BY4.0 evidence are pinned separately. Native
LAMDA/BASECOL/KIDA and atomic/publication records with unclear managed rights
remain metadata-only. ATcT403 and EuropePMC500 responses are retained acquisition
ERRORs, never successful source objects. CODATA2022's primary article is a2025
publication with its own rights-reserved notice. Public access/search snippets
are not substituted for a numerical source or permission. Exact URLs, hashes,
native roles and unresolved conditions are in
`utility/chemical_database/PHYSICAL_INPUT_STATUS_20261119_R4.md`.

### Metadata and provenance corrections

Fresh MYTHOS native definition SHA
`4da22ac7d223b6b92174f716e52c80ecf4fa77ca8387d21f56c457f455816f4c`
matches version20240526 and explicitly labels131898.260953cm^-1 as the
transition-completeness energy. Comparing it with the maximum state energy
132196.65537cm^-1 as a source defect was incorrect. A reproducible r4-only
metadata amendment corrects this interpretation; every old native/CSV member
and source ID remains byte-identical, and old release metadata is unchanged.

New source-specific readers preserve lexical fields and refuse changed hashes,
bad schemas, missing state references and overwritten outputs. Their exact
helpers/dependencies are archived beside their recipes. Independent review
found that a held data file could otherwise be declared a parser dependency;
the admission path now rejects native/derived/rights payload aliases by resolved
path or hash, authenticates software dependencies and checks Python syntax
without executing recipes. Regression controls cover the refusal. No held
payload was copied by the actual selected packets.

### Verification, commands and preserved identities

No AMPS application rebuild or linked U/I execution is required for this
source-catalogue/parser package. The portable shared C++17 API/SDK was compiled
with the normal shared test runner; compiler/source/build identities and all
stdout/stderr are retained. No PIC/MPI/SPICE configuration or root Makefile is
changed. The user's pre-existing Makefile edit remains untouched.

| Check | Actual outcome | Artifact |
|---|---|---|
| New publication and deep source integrity | PASS,14.095s;727objects checked | `build-result.json`, `staged-r4/publication.json` |
| Full shared software/API/source suite | 341PASS,0FAIL/ERROR/SKIPPED;69.757s | `software-qa/result.json` |
| Registered CORE protocol scope | 293PASS,0FAIL/ERROR/SKIPPED | `core-protocol-qa/result.json` |
| All selected packets/recipes and all available CSV Python/C++ comparisons | 26PASS,0FAIL/ERROR/SKIPPED;113.099s;112tables/2653685exact cells | `release-qa/result.json` |
| Original source/table preservation | 782old sources and105old table payloads unchanged | `release-delta.json` |
| Installed package and fresh installed CLI | PASS;412installed members | `installation.json` |
| Original installation and actual new installed API | 494protected old members unchanged;412members equal verified staging;9new table queries PASS | `installed-verification.json` |

The CORE293checks are a subset/scope repeat, not293additional independent
physics tests. Per-lane source/parser results are separate, overlapping source
QA. Full software includes actual relocated CLI/SDK checks. Release QA checks
all20selected packet declarations,16archived recipe/dependency members and
refusals for all618held sources/26held tables in both APIs. An initial helper
setup error used the wrong preservation-manifest key; its ERROR receipt is
retained. The corrected check uses the original pre-task494expected hashes,
not recalculated expected values. No installed data was changed by that error.

Exact publication and20packet arguments are in `publication-plan.json`.
`software-command.json` records all actual-source environment selectors. The
principal commands, from the repository root, are:

```sh
/usr/bin/python3 -B utility/chemical_database/build_database.py --amps-root "$PWD" --validation-root /data/vtenishe/moon_validation_data --output-root NEW_PUBLICATION --release-id NEW_RELEASE --include-model-definitions --source-packet PACKET_JSON
/usr/bin/python3 -B test/chemical_database/run_all_tests.py --database-lock DATABASE_LOCK --output-dir NEW_SOFTWARE_QA
/usr/bin/python3 -B test/chemical_database/run_all_tests.py --protocol-only --output-dir NEW_PROTOCOL_QA
/usr/bin/python3 -B test/chemical_database/verify_physical_model_release.py --database-lock DATABASE_LOCK --output-dir NEW_RELEASE_QA --native-probe CATALOGUE_PROBE --source-packet PACKET_JSON
/usr/bin/python3 -B utility/chemical_database/install_database.py --staged-root VERIFIED_PUBLICATION --destination NEW_SNAPSHOT --receipt NEW_RECEIPT
```

Repeat `--source-packet` with the full frozen family in the plan. Use new output
directories. Runtime selection and exact release hashes are:

```text
release: amps-chemistry-20261119-r4
release SHA256: 08153b0f7efb888260628adf87d9f4e80df4728ceb50d77303a7c26a846cdcba
profile SHA256: 8aafa28cbd73e12262c2a1ac4f4eafcbc414b21141ad7fc25665b2c47efdbfaa
export SHA256: a43c260082bdcdc0b2cd2bf2f64a3c1758e59e1af48d7d35b66148c30429517c
lock SHA256: db5d167351c31244a86df61d5eaf3fe0f7211dbd49cc5fb1e413ccc8e9e4b38a
```

### Remaining M3 requirements and roadmap control

New state/transition arrays reduce source gaps. Physical molecular admission
still needs explicit formation/threshold/asymptote/state and compatible-mass
crosswalks, populations, complete exits or omitted-exit bounds, and uncertainty/
model-error scope. The approved heavy-ion reservoir does not require a full
joint outgoing-electron PDF, but does require compatible state/rate/domain/
moment and ledger/error evidence. The new Na composites do not supply complete
CEX products, the missing D119990 provenance/uncertainty or a forward-angle
error bound. The approved ideal collector remains distinct from a real lunar
retention/neutralization/accommodation/release law.

Full linked numerical/statistical matrices, global budgets and same-campaign
M1/M2 predecessors remain required. Nothing in this source/software package
constitutes formal I21–I24 PASS or observational validation. No lunar test ID,
criterion, calibration/hold-out membership or roadmap is changed. No git stage,
commit or push is made.

Roadmap hashes remain source
`712e1152bf50af6d9ee1dad8f8a585809ca8b92cf42ac6fc7a009f6c0253ee17`
and canonical data
`d10eca1c539f9cf82940236d9e1aa36af665d43e82ba475925950d582fbf49c4`.
Approved source development/resource amendments govern implementation order;
canonical scientific/data/predecessor contracts govern qualification.


## 5.40 — Commit staging of the additional physical-input update

The user requested staging of the completed r4 source-input work without a
commit. The reviewed selection contains26files: shared chemical-database
readers/publication/provenance code, five new source/amendment test modules,
existing parser/release-audit tests and the corresponding contracts/READMEs/
STATUS. The source revision is
`26f1a65c58a8a70f1c6301b1a6c4d283a4558c94`. No AMPS core `src/` file is part
of this selection.

The local Makefile, generated inputs/executables, archives, dependency trees,
external source payloads and all `test_output/` artifacts remain outside the
index. Source/installed-data scientific scope and the341software/293protocol
subset/26release-audit PASS results in section5.39 are unchanged. All selected
Python files pass a fresh syntax check. Staged whitespace, exact selected
working/index bytes and preservation of excluded paths are checked after
staging. The operation makes no commit or push and changes no external database
selection. Its detailed before/after selection and hashes are retained in a
unique `/tmp/amps-chemdb-stage-*` review directory.


## 5.41 — Targeted remaining-physics search and r5 source integration

**Research/source-access package COMPLETE; formal M3 remains NOT QUALIFIED.**
The user requested online databases first and relevant primary publications
for inputs still missing from M3. Work began2026-11-19T23:54:59Z at revision
`ef73fd281f279e7538a7ead67620c02a629f0270`. Three agents independently covered
molecular energies/states, EI/CEX and surface/relaxation; the root covered
atomic sufficiency, common source policies, publication and verification.
Artifacts are in
`test_output/chemical_database/m3-targeted-inputs-20261119T235459Z/`.

### Findings that change the implementation plan

The [shared assessment](../utility/chemical_database/PHYSICAL_INPUT_STATUS_20261120.md)
and [M3 follow-up](test/M3_PHYSICAL_INPUT_FOLLOWUP_20261120.md) separate available
data, derivation/admission, a new model choice and genuinely unresolved
evidence. Detailed per-lane reports and source ledgers preserve actual native
HTTP/UTC/hash/rights/state/unit records. No original source, reaction rate,
acceptance criterion or model input is silently changed.

- Existing isotope-specific He/Ne/Ar ionization primitives are available;
  three retained-primary/hash/SI-conversion checks PASS. Profiles and energy
  crosswalks need implementation. Compatible ion masses are derivable, so
  separately measured ion masses are not a necessary missing input.
- Molecular binding/ionization/appearance energies have primary anchors.
  The report derives all eleven frozen H2O/OH reference thresholds symbolically
  using one energy zero and explicit states. H2+ is not a product of this
  network. Existing OH/OH+/H2 and O data are reusable rather than new-data gaps.
  Product/reference/spin selection, missing printed molecular-threshold
  convention and H2O04–06 dissociative-ionization partition remain explicit
  model choices. The approved neutral microcanonical law does not automatically
  authorize a charged/electron three-body partition.
- Bulk inherited-heavy-velocity electron-reservoir bookkeeping can proceed
  from selected scalar rates, states/Q and signed exchanges without a full
  outgoing joint PDF. Source-fit conditional moments are derivable but not
  measured moments or automatically justified microscopic samplers. The
  rounded Voronov Na I=5.1eV is below the physical threshold; no silent cutoff,
  renormalization or altered rate is authorized.
- Igenbergs2008's actual2007 preprint recovers source equations/domain/state
  interpretations for previously held Na EI/composite fits. Elow is distinct
  from physical ionization I. The coefficient audit has **104PASS/2FAIL**:
  two halfway printed-rounding differences remain unresolved and unchanged.
  An extra native composite A7 and exact keV versus eV/u bindings still require
  care. These optional source comparisons are separate from software QA and
  the unchanged selected Voronov production coefficients.
- Reference CEX Q and all-angle bounds are derivable. D119990's original
  experiment/preparation/uncertainty and prompt-versus-eventual H1s ambiguity
  remain unresolved; no small forward-angle accuracy bound was obtained.
- Source-defined OH first-order and heterogeneous-site water-release models
  exist. KIDA and the final Minissale review agree on the paired OH model;
  the older manuscript pair is not mixed with it. OH diffusion/recombination
  barriers are not desorption energies or impact sticking. The new water
  sticking observations retain their gold-relative calibration and20K domain.
  Matching OH impact/accommodation and extensions beyond measured conditions
  still need evidence or reviewed model scope. Explicit collisionless
  reference conditions have zero collisional quenching; nonzero baths require
  applicable state/rate/density bounds.

### New shared source release and preservation

Installed snapshot:
`/data/vtenishe/CHEMICAL_DATABASE/snapshots/amps-chemistry-20261120-r5`.
The original default/r2/r3/r4 locks remain unchanged. This release adds29source
records, including7allowed payloads and22metadata-only records, plus two CSV
views. Totals are52datasets,914sources,274admitted payloads,640metadata-only
records,145table descriptors,117available tables (114CSV), and the same12inert
model definitions. All885previous source records and143table payloads retain
their identities/bytes.

The NIST H/O evaluation's official PDF, native PMC HTML and notice are admitted
under a narrowly reviewed technical-publication policy. Its Table3 preserves
208cells of OH formation-enthalpy material at298.15K with historical methods;
Table13 preserves375cells of mixed ion-energy quantities. No source value,
confidence convention, isotope or zero-K state is selected by the reader.
The new maintained `native_nist_hxoy.py` and its automatically discovered test
retain spans, original HTML/text and footnotes. All five actual selected
source/parser controls PASS.

Figshare13274840.v1/file25563521 is a separate narrow CC BY4.0 source. Its
2136native bytes are preserved with producer MD5/SHA. It is XPS-derived
sticking relative to gold assigned S=1, with20±1K surface/300K gas. The
undeclared eight-bit encoding is not forced into UTF8/CSV; ordinary source
access is available. New policies check exact article/file/grant identity;
unrelated Figshare or NIST/SRD objects are not cleared by them. NIST requested
attribution, source notices and existing holds remain attached.

### Actual verification and retained failures

| Scope | Final result | Artifact |
|---|---|---|
| Source publication | PASS,13.990s | `build-v4-result.json`, `staged-r5-v4/publication.json` |
| Full shared software/source suite | 355PASS,0FAIL/ERROR/SKIPPED;68.945s | `software-qa-v4/result.json` |
| CORE protocol subset | 306PASS,0FAIL/ERROR/SKIPPED | `core-protocol-qa/result.json` |
| All selected packets/recipes/CSV Python–C++ audit | 30PASS,0FAIL/ERROR/SKIPPED;114.882s | `release-qa-v4/result.json` |
| Installed source/CLI integrity | PASS;752objects checked | `installation-r5.json` |
| Previous releases and installed equality | 906protected previous members unchanged;432installed files equal verified staging | `installed-verification-r5.json` |
| Optional Na source-coefficient correspondence | 104PASS,2FAIL retained | `collisions/coefficient_crosswalk.json` |

Software scopes overlap and are not independent physical tests. The release
audit compares all114available CSV tables and2658349exact cells through the
actual Python/C++ APIs, authenticates all24selected packets/17recipe members
and preserves held-source refusals. Actual installed Python/CLI queries also
verify both new table shapes and the original sticking-file producer checksum.

The first table-registration packet reused T13's ID for T3; publication
correctly ERRORed before writing a completed lock. Its next candidate used
a documentary `notes` field that conflicted with authenticated runtime notes;
the actual API correctly refused it (354PASS/1ERROR in that retained suite).
New packet versions correct identifiers/field names only; all source/CSV
values and recipe bytes are unchanged. Old manifests/staging/logs remain
preserved. The selected final packet is `molecular/source_packet_tables_v4.json`,
SHA`52d3f7055deb3a4636ac287f18c2a3ec016fb518d66c1afb5c0aee75ee6a1e03`.
These corrected infrastructure errors do not erase the separate two unresolved
source-coefficient FAIL comparisons. Acquisition403/404/API-error responses
also remain in their original manifests, never successful scientific sources.

`publication-plan-v4.json` records the complete24packet publication command;
`software-command-v4.json` records exact CLI/native source selectors;
`release-qa-command-v4.json` records the actual C++ probe/audit invocation.
Normal shared C++17/Python source tests were used. No AMPS application build,
MPI cohort or formal U/I execution was performed for this source package.
The existing application configuration, core files and local Makefile remain
as found. No git staging, commit or push is performed by this package.

```text
release: amps-chemistry-20261120-r5
release SHA256: 484268dc6f13cbdaeae519bdfc021dc166339474e619ffd8b0972b23c3abdcce
profile SHA256: eee2d54134db7b3fc1b6eb6a772d16e9ae3d61ca4355cbc69dc845bd9a39dc73
export SHA256: cd9cf3675b3ea0e2fba2541451f9e43e931fbc5e3166b2ba5d357c80a5ebe242
lock SHA256: 7b588795ba3dcbd17491d5afbb9a417f62824970afa94b10645e8eda835799a2
```

Formal M3 still requires actual production integration, frozen numerical/
statistical/budget qualification and applicable same-campaign M1/M2 gates.
The approved P01–P09 and numerical/statistical choices remain in force; no
renewed blanket approval or missing full-joint distribution is introduced.
Roadmaps remain different and unchanged: sourceSHA
`712e1152bf50af6d9ee1dad8f8a585809ca8b92cf42ac6fc7a009f6c0253ee17`,
canonical dataSHA
`d10eca1c539f9cf82940236d9e1aa36af665d43e82ba475925950d582fbf49c4`.
Source development/resource amendments govern work order; canonical scientific/
data/predecessor rules govern qualification.

## 5.42 Nightly chemical-database and Moon registration audit — 2026-11-23

**Audit complete; coverage limits remain explicit.** Source revision
`ef73fd281f279e7538a7ead67620c02a629f0270` plus the existing working changes.
No production/test runner, registry, application configuration or reference
was modified. Only the shared-test README, lunar-test README and this status
record were updated. Artifacts and exact commands/source hashes are under
`test_output/srcMoon/nightly-registration-audit-20261123T233024Z/`.

The source trace is `utility/nightly/runners.json`'s `core` entry
→ `src/test/run_tests.py --campaign nightly`
→ manifest owner `CORE-CHEMICAL-DATABASE-001`
→ `src/test/tests/models/chemical_database/run_test.py`
→ `test/chemical_database/run_all_tests.py --protocol-only`.
The owner is registered in smoke/nightly/extended, requires Python/C++17, and
uses the configured compiler. It reports one CORE case and preserves its
individual method results/counts in the nested shared-protocol artifact.

The read-only server SQLite query and original archived suite log confirm
actual scheduled execution, not registration alone. Latest observed run
`20261123T0201Z-passion-gnu` tested committed master
`b4980503e35c498d42a689188e0f6fc2e4913d94` and finished
2026-11-23T13:11:43Z. Its chemical owner PASS contains **293 PASS, 0 FAIL,
0 ERROR, 0 SKIPPED** in 22.821 seconds. Original log:
`/data/vtenishe/AMPSTEST/www/archive/2026/11/22/passion-gnu/logs/core/files/run/logs/CORE-CHEMICAL-DATABASE-001.log`.
The original bytes/hash are retained in this audit. The three preceding
archived chemical-owner entries also PASS.

A fresh execution of the same nightly adapter against this working tree used:

```bash
AMPS_TEST_TEMP_ROOT="$PWD/test_output/srcMoon/nightly-registration-audit-20261123T233024Z/work" \
/usr/bin/python3 -B src/test/run_tests.py \
  --config test/qualification-machine/passion.json \
  --campaign nightly --test CORE-CHEMICAL-DATABASE-001 \
  --jobs 1 --streams 1 \
  --output-dir "$PWD/test_output/srcMoon/nightly-registration-audit-20261123T233024Z/core-nightly" \
  --json "$PWD/test_output/srcMoon/nightly-registration-audit-20261123T233024Z/core-report.json"
```

It reports **1 CORE suite/case PASS**, containing **306 PASS, 0 FAIL,
0 ERROR, 0 SKIPPED**; child scope 25.133 seconds, outer command 25.436 seconds.
The 13 additional local methods (nine source-policy and four NIST parser
controls) all PASS. They remain uncommitted and therefore are not yet included
in the scheduled committed-master checkout. This is software/native/API scope,
not physical chemistry qualification. No AMPS build, MPI simulation, lunar U/I
execution, new chemical-data acquisition or observation processing was needed.

| Entry point | Actual selected/registered scope |
|---|---|
| Nightly `moon` → `srcMoon/test/run_tests.py --all` | 26 U cases, U01–U26 |
| Local `srcMoon/test/run_all_test.py --list` | 40 cases, U01–U26 plus I11–I24; registry validation finds no error |
| Global `APPLICATION-MOON-REGRESSION-001` | 36 cases, U01–U26 plus I11–I20; extended/application-smoke, not nightly |
| Global/local M0 owner | Separate I01–I10; not substituted by the U runner |

The latest archived `app-moon` family has exactly **26 records: 15 PASS,
11 SKIPPED, 0 FAILED**. The separate legacy `planets` entries Moon,
Moon-restart and Moon-tetra-cut-cells all PASS; they are additional legacy
regressions, not additional U IDs. Linked I11–I24 are not selected by the
nightly Moon entry, and I21–I24 are also missing from the global extended
adapter's case registry. These are documented coverage limits, not new tests
or silently expanded nightly work.

The current CORE protocol selection explicitly excludes **45 authoritative-
source methods**, recording their IDs in `outside_selected_scope`, and does
not select the **four real-release PythonProduction methods**. The independent
full acquired-release/source-packet audit is likewise not invoked nightly.
Those scopes require explicit pinned locks/original-file selectors and retain
their separate scientific/source-admission limitations. A CORE PASS does not
claim that the installed r5 release, all source bytes or linked chemistry were
qualified in the scheduled nightly.

The audit's own seven registration/history/execution checks PASS; they are
not permanent U/I identities or scientific acceptance tests. Full child
results are in `core-nightly/artifacts/CORE-CHEMICAL-DATABASE-001/shared-protocol/result.json`;
`historical-nightly-results.json` records the read-only archive query;
`result.json` records coverage, excluded method IDs and the missing global
I21–I24 registration. No data/reference hash, acceptance criterion, hold-out
membership or roadmap was changed. Source/data roadmap hashes remain
`712e1152bf50af6d9ee1dad8f8a585809ca8b92cf42ac6fc7a009f6c0253ee17` and
`d10eca1c539f9cf82940236d9e1aa36af665d43e82ba475925950d582fbf49c4`.
Their scoped development/resource amendments control execution order; the
scientific/data/predecessor contracts control qualification. Formal M1/M2/M3
status is unchanged. Nothing was staged, committed or pushed.

## 5.43 I21 homogeneous execution and I22/I23 preparation — 2026-11-24 UTC

Source revision `267fdc25fba6f9b510143425f390958801a64cd1` plus the recorded
changes. The work package is under
`test_output/srcMoon/m3-i21-execution-606xjgOZ/`. The two governing roadmap
SHA-256 values are recorded in the
[I21 execution contract](test/integration-tests/I21/EXECUTION_CONTRACT.md).
The source roadmap supplies the approved early M3 development/resource
exception. The common scientific roadmap supplies the frozen gates and data
semantics.

### Implemented scope

- `integration_m3.py` and the Moon aggregate accept one frozen I21 homogeneous
  cohort plan and its expected digest. The existing `m3_qualification_executor`
  selects the matching private source build and preserves the nested native
  artifacts. A qualification-candidate cohort requires U14 and same-campaign
  I13–I20 PASS. A development cohort retains those statuses and runs with
  diagnostic endpoint survival. The whole I21 gate remains SKIPPED.
- `I21/prepare_cohort.py` stages the existing sodium production owner at
  8/16/32 dyadic steps and both ranks. The same 65,536 s, 128 original parents,
  weight 2.625, seed 731993 and quiet fixed-point geometry from the commissioned
  case are used at all six levels. D04 source, geometry receipt, actual producer
  and the ordered set of eight SPICE kernels are pinned in the final plan.
- `prepare_atomic_family.py` builds independent ground-ion references for
  ⁴He, ²⁰Ne and ⁴⁰Ar and checks them with the actual production state resolver.
  The original isotope ledger and PDFs, tested base profile and D04 inputs are
  hash-bound. The derived threshold differences and original uncertainties are
  retained. I22 still needs the explicit photon product-state/energy policy.
- The shared Voronov fit exposes the analytic `d ln k / d ln Te` on its native
  domain. U16 independently checks 132 points using a high-precision difference
  of the direct rate expression. I23 still needs the source-to-state threshold,
  applicable bath and model-error scope, and production event wiring. The
  existing Maxwellian rates and frozen acceptance criteria are preserved.

### Retained execution and open gates

The first I21 matrix produced six completed native runs. Its independent
audit reports **18 PASS, 1 FAIL, 0 ERROR, 0 SKIPPED**. Source rates, chronology
and event ledgers passed. There were 85 survivors and 43 reactions in each
matched launch. These are the same parent histories across the six launches.
The observed endpoint fraction is diagnostic. The separate P02 audit reports
**6048 PASS** on those six realizations (258 repeated event records).

The single independent numerical check failed at both MPI layouts. The
maximum normalized errors are `7.3782e-17`, `7.3782e-17` and `9.0911e-17` at
8, 16 and 32 steps, respectively. These are well below `0.002` and do not
decrease monotonically. Inspecting the exact hexadecimal native event times
gives `4.9572e-17`, `4.9572e-17` and `4.9725e-17`. The plateau persists before
the parent CSV's binary64 projection. The frozen scorer remains strict. An
independent varying-hazard case is still needed to test the step-dependent
exposure error separately from this constant-hazard clock diagnostic.

The standalone results are **70 PASS** for the native atomic profiles,
**302 PASS** for U16, **37 PASS** for the timing/executor controls and
**306 PASS** for the shared core chemical-database protocol. The U16 rate
error stays below `1.6e-15` and the largest normalized slope error is
`2.14e-16`. The isolated final `make -j test_Moon TESTCOMPILEJOBS=16` passes
with a zero-byte `test_Moon.diff`. The no-plan I21–I24 aggregate is 4 SKIPPED.
The explicit I21 qualification request skips before launch when prerequisites
are missing. A deliberately wrong plan digest produces the expected ERROR.

The first final-source build attempt reported a source mismatch on the
updated U16 probe. Its test-control source had been omitted when copying the
private checkout. The refusal is retained under `final-native/`. The rerun in
`final-native-retry/` includes both shared source edits and the updated probe.
The final source-matched rerun reproduced all six native completions, the
same **18 PASS / 1 FAIL** audit, and **6048 PASS** in the separate P02 audit.
The final matrix took 388.10 seconds and peaked at 42,933,399,552 bytes of
sampled process-tree RAM. The zero-byte nightly diff and both independent
attempts remain available. `completion.json` indexes the final results and
commands, including the precise native plan and source-data hashes.
Build and run commands, exact input hashes, compiler/ABI provenance and
resource logs are retained alongside each result.

Formal M3 still requires I13–I20 and the complete linked I21–I24 cases. The
frozen distance/shadow, force, weight and statistical cases remain open for
I21. I22 still needs the molecular reference-state and charged-partition
choices plus the U24 dependencies. I23 still needs source-supported bounds
for the reduced reservoir. I24 still needs documented product-state and
angular support for its selected channel. No new U/I ID, scientific criterion,
chemical-database installation or reference baseline was introduced.

## 5.44 Remaining M3 work proposal — 2026-11-24 UTC

The [M3 completion follow-up](test/M3_COMPLETION_FOLLOWUP_20261124.md) proposes
an independent varying-hazard case for I21, atomic and state-defined molecular
I22 cases, a controlled ⁴He I23 pilot, and a ranked evidence-acquisition plan
for I24. It retains the I21 constant-hazard numerical FAIL and gives separate
source-admission, linked numerical and statistical gates. The proposal also
records the applicable M1/M2 predecessors, the one-hour development budgets
and the separate powered cluster campaign.

This is a documentation package against revision
`267fdc25fba6f9b510143425f390958801a64cd1` with the section 5.43 work in the
working tree. Proposed case and physical-model choices are pending review.
Formal I21–I24 remain SKIPPED. See the proposal for both governing roadmap
hashes and exact data lock. No build or campaign was run for the proposal.

## 5.45 M3 follow-up implementation and verification — 2026-11-24 UTC

The work follows the [approved completion proposal](test/M3_COMPLETION_FOLLOWUP_20261124.md)
on revision `267fdc25fba6f9b510143425f390958801a64cd1` plus the recorded
working-source hashes. Artifacts are under
`test_output/srcMoon/m3-followup-MfThQfny/`. Source roadmap SHA-256
`712e1152bf50af6d9ee1dad8f8a585809ca8b92cf42ac6fc7a009f6c0253ee17` supplies
the development and resource exception. The canonical roadmap SHA-256
`d10eca1c539f9cf82940236d9e1aa36af665d43e82ba475925950d582fbf49c4`
controls the shared scientific criteria and predecessor gates.

### I21 varying-hazard case

`prepare_cohort.py`, the existing aggregate owner and the independent sodium
auditor now support an explicit straight-line trajectory plan. The native
adapter retains launch velocity, zero-force evidence and spherical occultors.
The independent oracle proves full finite-ray visibility, brackets the
transverse Sun offset, locates parent-keyed first hits and checks event
position, source/event rates and the committed ledgers. The varying plan
uses integrated exposure; the homogeneous statistical scorer refuses it.

The prescribed radial fixture reports **44 PASS**. Its normalized first-hit
errors at 8, 16 and 32 steps are `7.70281e-4, 1.96110e-4, 4.96111e-5`.
The linked compact case reports **19 PASS** after six actual normal source-
bound launches at the same time levels and 1/8 MPI, OMP1. Exact native time
errors at both layouts are `9.95673e-11, 2.47282e-11, 6.14480e-12` normalized
by the 65,536 s horizon. The largest event-position discrepancy was
0.00168 m. The scope uses 128 initial histories and seed 731993. Each
matched launch retained 77 survivors and 51 reactions. The additional
independent mechanical audit reports **7,152 PASS**. The development
statistics are diagnostic. Re-auditing the original constant case gives
**18 PASS / 1 FAIL** with its original roundoff-plateau result.

The first attempt encountered restricted local MPI sockets. The second
attempt with MPI access stopped at the host available-memory reserve. The
retained retry increased near/far/angular/polar AMR target sizes eightfold
for this straight-line test only. It preserved the physical sphere,
ephemeris, forces, transport horizon and source-data identities. The AMR
tree dropped from 11,684 to 792 blocks. The final six launches and audit
took 190.39 s and peaked at 28,427,436,032 bytes of sampled process-tree
RSS. See the [case contract](test/integration-tests/I21/VARYING_EXPOSURE_CONTRACT.md)
and `trajectory_compact_command.json`, the exact native/build manifests
and CSV/PNG/EPS under the successful compact result.

### I22, I23 and I24 source and implementation scope

- The I22 ground-ion family re-ran with **70 PASS** and emitted six candidate
  effective-photon rows from exact D04 quiet/active data and validated
  state profiles. The rows record the source and physical thresholds
  separately. Atomic photon event admission still needs a supported
  mapping or bounded reduced closure. Full I22 also needs the three
  charged molecular breakup laws, populated states and U24 surface closure.
- The compiled I23 helium provider connects the actual D05 coefficient
  and analytic fit slope to the shared Maxwellian moment and signed
  electron-bath planner. It requires a stationary target and zero drift
  at 10 eV. Composition checks authenticated model-error and state
  evidence and uses the accepted interval and original represented
  weight. The probe reports **24 PASS** for independent rate/mean/threshold
  comparisons and admission guards. Its coefficient is
  `7.760630099173217e-16 m³/s`. Its model incident mean is
  45.41818 eV. No actual helium event was admitted in this pilot.
- The optional I23 source comparison extracted the neutral helium row
  from [Kim's NIST paper](https://doi.org/10.6028/jres.097.032), PDF SHA-256
  `b777277069851715e5810c17f36ad8fcbaf005fd6c0071c2b9baba28082be961`.
  Integrating the paper's fit over its stated finite range gives a
  coefficient 6.35% below the prescribed Voronov value and a conditional
  mean 1.20030 eV higher. This is a fit-to-fit sensitivity. It leaves
  the D05 coefficient reference, actual source-uncertainty bounds and
  residual-ion population evidence as distinct obligations.
- The [I24 source review](test/integration-tests/I24/SOURCE_ADMISSION_20261124.md)
  verified the already archived r5 D119990 typed and extracted bytes
  and found state-resolved Na/He leads. It produced zero new numerical
  collision records. The preferred Na channel still lacks primary
  provenance and angular support. Finite proton–He lab angle plots are
  candidates for source data recovery, with no claimed physical angle
  cap. No new chemical-database installation was made.

### Registry, provenance and remaining gates

The global extended Moon suite includes U01–U26 and I11–I24 (**40 cases**).
The nightly Moon U-series selection retains its 26 cases. The catalogue
generator and traceability validator pass. The focused chemical-database
core suite passes all **306** selected shared protocol controls. The
cohort qualification suite passes **31** checks. A no-plan aggregate
selecting all four M3 IDs gives **0 PASS / 0 FAIL / 0 ERROR / 4 SKIPPED**.
The global adapter's unavailable-SPICE control reported all four exact
I21–I24 IDs with typed SKIP.

The source-qualified SPICE build completed in 181.18 s in a disposable
checkout from a clean generated tree. Compiler, ABI, source, input,
D04/D05 and SPICE kernel hashes are in the corresponding manifests.
The clean isolated `make -j test_Moon TESTCOMPILEJOBS=16` run returned
zero in 185.18 s and produced the zero-byte `nightly-final/test_Moon.diff`.
The earlier attempt that reached the host memory reserve is retained
under `nightly-clean/`.
The optional NIST source and source-database r5 identities are retained
alongside the corresponding calculations. No original observation,
acceptance criterion, reference, or permanent U/I identity changed.
The package's `completion.json` indexes the evaluated outcomes and exact
command, manifest and plot files.

M3 remains open. Complete the rest of I21 (solar endpoints/distances,
shadow, both movers, weight and powered survival cases), the pending
physical source packages and I22–I24 production activation and
qualification. I13–I20 remain same-campaign prerequisites. The
development numerical PASS does not supersede the constant-clock FAIL.

## 5.46 M3 archive and publication search — 2026-11-24 UTC

The [source search report](test/M3_ARCHIVE_SOURCE_SEARCH_20261124.md)
records the datasets, source-domain limits, exact hashes and reproduction
commands against revision `e5f2978ac0576f42dc510e0da110dcbe0029d69a`.
The live IAEA D119990 catalogue reproduced the 12 frozen sodium partial
cross sections. Groningen MOTRIMS measurements identify the same
ground-state reaction and report differential extractions from 7 keV/u,
with published peak-normalized figures. Their numerical differential
bins, acceptance, uncertainty, normalization and D119990 lineage are
still needed for I24. The verified CC BY 4.0 helium resonant-exchange
archive contains 661 center-of-mass energy tables with 768 angles each.
Native quadrature reproduces all archive integrals. A separate helium
case would require a reviewed scope.

NORAD supplies ground-residual and total He photon cross sections, and
Curtin single-ionization cross sections give a finite-table incident
mean of 46.08014 eV for the 10 eV Maxwellian helium diagnostic. Broadband
state weighting, joint charged-water product kinetics, incident OH
surface response and the I23 physical error envelope remain open.
Source bytes and audit results are in
`test_output/srcMoon/m3-source-search-9wqcZyI0/`. This package ran
zero builds and zero linked cases. The formal I21–I24 status is
0 PASS, 0 FAIL, 0 ERROR and 4 SKIPPED.

## 5.47 Shared chemistry archive integration — 2026-11-24 UTC

The new snapshot is
`/data/vtenishe/CHEMICAL_DATABASE/snapshots/amps-chemistry-20261124-r6`.
Use its explicit `database-lock.json` and README section 12 to access
the tables from any application. This run used AMPS revision
`e5f2978ac0576f42dc510e0da110dcbe0029d69a` with the
[new source packet recipe and contract](../utility/chemical_database/NATIVE_M3_ARCHIVE_SOURCES.md).
Both roadmap hashes and source/admission scopes are recorded there.

The r6 packet adds 15 records, with five archived payloads and ten
metadata-only entries. The two new CSV views contain 661 native
center-of-mass He⁺–He charge-transfer energy rows and 507,648 angular
rows, covering 0.0001–10,000 eV. The original Zenodo ZIP is attributed
and stored under its deposit-specific CC BY 4.0 grant. All 661 native
quadrature sums agree with printed integrals within `4.60e-11` relative.
The live IAEA HTML is archived as corroborating evidence for the
existing 12-knot sodium partial. It reproduces the same numerical
SHA-256. The replay retained all r5 sources, table payloads, models,
reactions and species. Four old IAEA processing paths reflect the
release relocation, with their source and recipe hashes preserved.

The general snapshot retains metadata for NORAD's noncommercial He
cross sections, Curtin CCC and the other held sources. Molecular daughter
data, sodium differential bins and uncertainty, and physical error
bounds remain open. The snapshot selects no reaction and does not alter
the original lock or the supported production profiles.

Artifacts are in
`test_output/chemical_database/archive-integration-20261124-vUSDrSTa/`.
`build-command.json` pins the 24 previous packets and the new packet.
`packet-final.stdout.json` and `r5-preservation-verified.json` give
the native and preservation checks. `installation-r6.json` records
the verified external copy and its lock SHA-256
`c7c87fb0b00f98c1c402e552b297347bd1ba56e935b634b3de04c68a1ac07689`.
The full shared software suite
reported **337 PASS, 6 SKIPPED, 0 FAIL, 0 ERROR**. All three explicit
real-source methods passed. Missing optional inputs from other sources
account for the skips. The actual nightly CORE adapter passed its
one registered case with **307 passing protocol checks** and zero
nonpasses. Formal I21–I24 remain SKIPPED. No AMPS application build or
linked Moon run was required.

## 5.48 M3 publication follow-up and shared r7 release — 2026-11-24 UTC

The [source report](../utility/chemical_database/NATIVE_M3_LITERATURE_SOURCES.md)
records the archive search, native SHA-256 hashes, permissions, scientific
limitations and the offline extraction command. The release uses AMPS
revision `e5f2978ac0576f42dc510e0da110dcbe0029d69a` plus the recorded
working-source changes. The report identifies which roadmap amendments
control development order and which shared criteria control qualification.

r7 adds 11 source records and four queryable tables with 12 published rows.
The NBS water onsets and relative ion signals, Kim's neutral H/He total
ionization fits, and Cui's measured and calculated helium capture fractions
retain the exact native values. The Cui excited measurement at 0.5 keV/u
remains blank. Existing r6 data and selectors were replayed and compared.
The full r7 catalogue has 54 datasets, 940 source records and 123 queryable
tables. The explicit lock is in
`/data/vtenishe/CHEMICAL_DATABASE/snapshots/amps-chemistry-20261124-r7`.
Its SHA-256 is
`949a06675f1d58a261f35c9f4eb4dd4b06c096e18d121a5173013115df1ec020`.

Artifacts are in `test_output/chemical_database/m3-source-r7-xBR6k1oF/`.
`build-final-command.json` records the exact packet selection and builder
command. `r6-preservation.json` records source and numerical-table identity.
The installer verified **457 files**. The four installed CLI exports match
the packet projections byte for byte. The full shared test runner produced
**341 PASS, 6 SKIPPED, 0 FAIL, 0 ERROR**, including all three r7 source
checks. The six skips are optional sources outside this acquisition.
The nightly CORE case passed with **308 protocol checks** and zero nonpasses.
No AMPS application was built or run. Formal I21–I24 remain SKIPPED.

I22 continues to require joint single-photon charged-water product data or
an independently bounded reduced model, the broadband atomic state mapping,
and incident OH sticking and accommodation. I23 requires independent
gas-phase product-state and recoil bounds. The preferred I24 sodium
channel still needs native angular coverage and provenance. The Wang CDS
table acquisition was unsuccessful. The Cui results support a separate
helium comparison within their reported domain.

## 5.49 Atomic and molecular source follow-up for shared r8 — 2026-11-24 UTC

The shared snapshot
`/data/vtenishe/CHEMICAL_DATABASE/snapshots/amps-chemistry-20261124-r8`
adds nine source records and four queryable tables with 218 rows.
All 940 r7 source entries, 151 table entries, 21 species, 17 reactions
and 12 model definitions were preserved. Four existing IAEA table
locators were relocated and their referenced recipe and native file
hashes agree with r7. The installer checked 472 files. The complete r8
catalogue contains 55 datasets and 949 sources, with 127 queryable
tables and 28 unavailable or held table entries. Select its explicit
lock, SHA-256
`0d5cf5a26f2c2b6d154e7732f6690993629b7325fbabf6c4983c6b4a95549f13`.
The default root and r7 locks retain their previous hashes.

The original Kalha workbook includes 103 elements. The four-element
CSV preserves 192 source combinations with native energy strings,
cross sections and blanks. Sutton's 200 eV He⁺ 4p/4s, 4d/4s and 4f/4s
ratios, their reported uncertainties and instrument-width dependence
give an independent laboratory observable. Bhatia's two theoretical
tables keep their native momentum convention and the suspect Table 8
exponents as printed. The archive and publication evidence leaves the
Wang helium angular dataset and Suárez water surface study as
metadata records. The former is embargoed. The latter treats
prepared H₂O⁺ dissociation and gives no joint photon-event product law.

This work used AMPS revision
`e5f2978ac0576f42dc510e0da110dcbe0029d69a` with working-source changes
to packet admission, the offline projection, shared tests and
documentation. The [source report](../utility/chemical_database/NATIVE_M3_PARTIAL_SOURCES.md)
provides both roadmap hashes, native source and rights record hashes,
and the extraction command. Artifacts in
`test_output/chemical_database/m3-source-r8-tU34BqBU/` include the
published packet `packet-final-v2/source_packet.json`, exact
`build-final-command.json` using Python 3.14, byte-preservation audit
`r7-preservation.json`, installer receipt and four matching installed
CLI CSV exports. No AMPS application build or linked Moon run was needed.

The full shared suite reported **346 PASS, 6 SKIPPED, 0 FAIL, 0 ERROR**
in 90.3 seconds. All five explicit r8 checks passed. The six skips
concern optional inputs outside this acquisition. Reproduce the suite
with `AMPS_CHECK_M3_PARTIAL_RELEASE=1 AMPS_CHECK_M3_LITERATURE_RELEASE=1`,
`AMPS_M3_ARCHIVE_PACKET` pointing to the r6 packet and
`python3 -B test/chemical_database/run_all_tests.py --database-lock <r8-lock> --output-dir <new-directory>`.
The shared nightly protocol case `CORE-CHEMICAL-DATABASE-001` passed
308 checks with no skips, failures or errors under
`python3 -B src/test/run_tests.py --config test/qualification-machine/passion.json --campaign nightly --test CORE-CHEMICAL-DATABASE-001 --output-dir <new-directory>`.
I21–I24 remain formally SKIPPED. I22 needs the broadband state weighting, charge-resolved
water products and incident OH surface response. I23 needs an
energy-dependent gas-phase residual-state and recoil bound. I24
needs the native ground-capture proton–sodium differential angles
and uncertainty. The helium data cover a distinct reaction.

## 5.50 I24 ground-capture angular-source follow-up — 2026-11-24 UTC

The [source review](test/integration-tests/I24/ANGULAR_SOURCE_FOLLOWUP_20261124.md)
finds published ground-product laboratory angular evidence from
7 keV/u upward. The reported 4 keV/u curves cover other channels.
The same source shows that the D119990 integral at 4 keV/u is far
larger than the published prompt ground-product partial. Both
conditions require channel-specific provenance and source data for
physical use. The review records the uniform solid-angle
center-of-mass option for a controlled sensitivity case and an
independent scale calculation for the sodium recoil. It admitted
zero angular rows and generated no linked or software test outcomes.
Formal I24 remains SKIPPED. Exact source identities and reviewed
scope are in the source review.

## 5.51 I24 evidence and isotropic model in the shared database — 2026-11-24 UTC

The shared r9 catalogue publishes two entries through the existing model API.
`proton-sodium-ground-capture-evidence-v1` records the D119990
rate/state crosswalk discrepancy and the available laboratory angular
coverage. `isotropic-com-two-body-scattering-v1` defines the experimental
solid-angle distribution in the pair center-of-mass frame. Its parameters
and assumptions are explicit, and it supplies no rate or event provider.
The [shared database contract](../utility/chemical_database/NATIVE_CEX_SCATTERING_SOURCES.md)
records the permitted source scope, both roadmap hashes and native
evidence identities. Formal I24 remains SKIPPED.

The D119990 table metadata links to the evidence entry. The preserved
12-node integral has SHA-256
`2a887ad43497a7f9e2b3d7b6173ae1b2597cc893d897926762c181946c92db93`.
No angular data rows were admitted. All 27 source packets, 949 source
entries, 155 table entries, and previously published species and reaction
records were carried forward. The authored definition table grew from 12
to 14 entries. The IAEA table notice is readable in C++ and Python.

The publication uses AMPS revision
`e5f2978ac0576f42dc510e0da110dcbe0029d69a` plus the recorded shared-source
changes. The explicit installed lock is
`/data/vtenishe/CHEMICAL_DATABASE/snapshots/amps-chemistry-20261124-r9/database-lock.json`
with SHA-256
`4f49bfcea74ebc015948d8068fc9315f953086037bcb23a65ddf51b291340cd8`.
Artifacts are in `test_output/chemical_database/cex-scattering-r9-2GDzNL/`.
`build-final-command.json` records the exact replay command and all
selected packets. The preservation audit checked the root and r8 locks
and the D119990 bytes. No AMPS application build or linked lunar test
was required. The first stage logged an error on a string-valued
processing field in an unrelated source. The successful retry adds the
review link at the native IAEA reproduction step.

The installer verified 474 files. The installed CLI resolved the two
entries and the original table. Its D119990 CSV export matches r8
byte for byte. The installed C++ `Catalogue::GetTable` returned the
scope note and original table hash. The nightly protocol case
`CORE-CHEMICAL-DATABASE-001` reported **309 PASS, 0 SKIPPED, 0 FAIL,
0 ERROR** using
`python3 -B src/test/run_tests.py --config test/qualification-machine/passion.json --campaign nightly --test CORE-CHEMICAL-DATABASE-001 --output-dir <new-directory>`.

The full shared suite reported **349 PASS, 6 SKIPPED, 0 FAIL, 0 ERROR**
in 89.2 seconds. The six skips required optional sources outside this
publication. Both explicit r9 checks passed. Reproduce with
`AMPS_CHECK_CEX_SCATTERING_RELEASE=1`, the r8/r7 real-source selectors,
and `python3 -B test/chemical_database/run_all_tests.py --database-lock <r9-lock> --output-dir <new-directory>`.
The results and installed query commands are in `shared-tests/result.json`
and `installed-queries.json` under the artifact directory.

## 5.52 OH surface-loss and helium excitation sources in shared r10 — 2026-11-24 UTC

The [lab benchmark report](../utility/chemical_database/NATIVE_M3_LAB_BENCHMARKS_SOURCES.md)
records the database and literature search, both roadmap hashes, native
identities, permissions, and the exact table extraction. Wang and
Bruggeman Tables 2–3 give six dry/humid OH flux/loss conditions and
two spatial-resolution estimates on fused silica in an atmospheric
helium plasma. The dry β estimates are 0.010–0.011 in Table 2.
The 1.0 mm humid diffusive flux and β remain blank. The paired
10 and 50 µm results quantify the effect of measurement resolution.
Sutton and Kay give three 200 eV ratios relative to He⁺ 4s with
reported asymmetric or symmetric errors. The visible 4p central
ratio in the preprint is 16.00. The earlier single-author r8 source
gives 19 for the same experimental programme and is separately
retained. Four other leads are catalogued as metadata.

The new snapshot is
`/data/vtenishe/CHEMICAL_DATABASE/snapshots/amps-chemistry-20261124-r10`.
Its lock SHA-256 is
`b685b842cb3549655c4f09596e8c628397d01de9539af7515a2ea41bed63728a`.
It adds nine source records and three data tables with 11 published
rows. All 27 prior packets, 949 r9 source entries, 155 table entries,
21 species, 17 reactions and 14 model definitions were preserved.
Four older IAEA locator pairs refer to the same hashed data and
recipes after relocation. The default and r9 locks were verified.

This source package used AMPS revision
`e5f2978ac0576f42dc510e0da110dcbe0029d69a` with working-source changes.
The directory `test_output/chemical_database/m3-source-r10-CMOILD/`
contains authenticated retrieval records, the publication packet
`packet-reviewed/source_packet.json`, the exact
`build-reviewed-command.json`, and `r9-preservation.json`.
No AMPS application build or linked lunar run was needed.

The installer verified 486 files. The three installed CLI CSV exports
match the reviewed packet byte for byte. An offline replay reproduced
the packet and its projections byte for byte. The full shared suite gave
**353 PASS, 6 SKIPPED, 0 FAIL, 0 ERROR** in 92.8 seconds. All four explicit
r10 checks passed. The six skips required optional sources outside this
package. Reproduce with `AMPS_CHECK_M3_LAB_RELEASE=1`, the preceding
r9/r8/r7/r6 real-source selectors and
`python3 -B test/chemical_database/run_all_tests.py --database-lock <r10-lock> --output-dir <new-directory>`.
The nightly `CORE-CHEMICAL-DATABASE-001` gave **309 PASS, 0 SKIPPED,
0 FAIL, 0 ERROR** under
`python3 -B src/test/run_tests.py --config test/qualification-machine/passion.json --campaign nightly --test CORE-CHEMICAL-DATABASE-001 --output-dir <new-directory>`.
See `installation-r10.json`, `installed-checks.json`, and
`shared-tests/result.json` in the artifact directory.

I22 still needs a source-qualified OH sticking/accommodation law for
the lunar surface, joint charged-water breakup products and
broadband atomic daughter-state weighting. I23 still needs
energy-dependent excited-state fractions and recoil bounds for the
configured electron distribution. The I24 sodium rate/state discrepancy
and the physical angular-data gap remain open. Formal I21–I24 remain
SKIPPED.

## 5.53 Shared M3 physical inputs and the D119990 figure identification — 2026-11-24 UTC

The r11a shared snapshot adds five source packets, 16 source entries and
three tables with 814 rows. The explicit lock is
`/data/vtenishe/CHEMICAL_DATABASE/snapshots/amps-chemistry-20261124-r11a/database-lock.json`,
SHA-256 `39ae1187e765e2417fbee9d07afb3b18aac60eaffe7c51e5e6aa2bbd7ba1ef5b`.
The catalogue now has 61 datasets, 974 sources, 161 tables (133 available),
21 species, 17 reactions and 14 model definitions.

| Gate | Physical result | Scope |
|---|---|---|
| I22, H₂O and OH | Moore–Mapaye give eleven solar-cycle OH/H₂O partial-rate fits, including the three charged-water breakup channels. At 584 Å (21.2 eV), Kronebusch–Berkowitz give OH⁺, H⁺ and O⁺ fractions of 0.29, 0.05 and 0.0009 of total ionization. | Source fits use TIMED/SEE and PHIDRATES. The laboratory detector discriminates against ions above 1 eV. Joint correlated electron–ion–neutral products and continuum state weighting are still open. |
| I22, atomic states | Laurell's 800 paired Ar peak-center fits give a mean separation of 173.3706 meV, consistent with the published 173 ± 9 meV. | One narrow-band laser experiment. Broadband populated He⁺/Ne⁺/Ar⁺ states remain open. |
| U24 / I22 surface dependency | Two laboratory ice uptake estimates and the distinct IUPAC evaluation are available. The experiment reports approximately 0.1 initially and 0.03 ± 0.02 after conditioning at 205–230 K. | Uptake includes reactive loss. The incident-OH sticking and translational-energy accommodation coefficients on lunar material remain open. |
| I23 | Two 200 eV electron-impact emission measurements give (5.6 ± 1.4) and (6.2 ± 1.7) × 10⁻¹⁹ cm² for the 30.4 nm He II line. | Emission is sensitive to cascade feeding. Direct all-state He⁺ production and recoil over a 10 eV Maxwellian still require additional evidence. |
| I24 | IAEA release member `119990.txt` names Knoop et al. 2004. NIFS CHART 38917 explicitly names their Fig. 2. Its outer-shell capture signal sums over outgoing H states. Reject D119990 as experimental H(1s)-specific evidence. Shingal's 1986 calculation gives 2.2 × 10⁻¹⁷ cm² at 4.0 keV/u. The 1987 larger-basis calculation gives 4.073 × 10⁻¹⁷ cm² at 4.00 keV laboratory energy. | Both are theoretical and have distinct native energy conventions. No uncertainty on the 4 keV ground partial or source-qualified angular law was established. The isotropic option retains experimental status. |

The [photon](../utility/chemical_database/NATIVE_M3_PHOTON_BENCHMARK_SOURCES.md),
[water](../utility/chemical_database/NATIVE_M3_WATER_BRANCH_SOURCES.md),
[ice](../utility/chemical_database/NATIVE_M3_OH_ICE_SOURCES.md),
[helium](../utility/chemical_database/NATIVE_M3_HELIUM_EMISSION_SOURCES.md) and
[sodium](../utility/chemical_database/NATIVE_M3_SODIUM_BENCHMARK_SOURCES.md)
reports record the online database search, original files, article locators,
hashes, rights, scripts, physical interpretations and replay commands.
The text of the weak water O⁺ measurement was checked independently on the
original printed page. The historical IAEA data remain byte-identical.
The catalogue displays the resolved figure and its measurement scope beside
the retained integral table.

AMPS source revision `e5f2978ac0576f42dc510e0da110dcbe0029d69a` with
working-source changes. The source roadmap SHA-256 is
`712e1152bf50af6d9ee1dad8f8a585809ca8b92cf42ac6fc7a009f6c0253ee17`;
the validation-root roadmap SHA-256 is
`d10eca1c539f9cf82940236d9e1aa36af665d43e82ba475925950d582fbf49c4`.
The source roadmap permits early M3 verification and clarifies development
resource limits. The common science criteria govern qualification.
Formal I21–I24 remain SKIPPED. This change builds the database and tests
the C++/Python interfaces. It does not require an AMPS application rebuild
or scoring against lunar observations.

Artifact root `test_output/chemical_database/m3-physical-r11-fS2iFdak/`.
The preliminary r11 release carried an earlier packet hash in one source
report. Its scientific data are identical to r11a. The preliminary suite
reported one test-harness failure because its assertion compared the
immutable API energy-range tuple against a mutable list. The test now
checks the published 80–1000 eV range using the API's declared type.
The original logs are retained. The `final/` subdirectory holds the
r11a build command, checks of all previous scientific identities,
installation receipt, export checks and test evidence. The default
database lock remains pinned to its existing snapshot.

From the checkout, replay the complete installed-source test with a new
output directory. The retained shell file sets the source selectors for
all five new packets and the original CEA, CHIANTI, NIST and ExoMol tests.

```bash
work=test_output/chemical_database/m3-physical-r11-fS2iFdak/final
bash "$work/test-release.sh" "$work/shared-tests-repeat"
python3 -B src/test/run_tests.py --config test/qualification-machine/passion.json \
  --campaign nightly --test CORE-CHEMICAL-DATABASE-001 --output-dir "$work/nightly-repeat"
```

The installed CLI exports of the three CSVs matched their packet bytes.
The installer verified 523 files and 871 object references. Every previous
source/table digest passed its preservation check. The native sodium
cross-section labels and the original 12 numerical knots can be inspected
alongside the new figure cross-reference. The nightly case passed its
309 source-independent checks with zero failures, errors or skips.

Two additional archive paths help delimit the atomic evidence gap.
The CDS Bochsler catalogue tabulates total H, He, O and Ne photoionization
frequencies at Earth. Heimann et al. measured relative He and Ne satellite
yields near individual thresholds for outgoing electron energies of 0–1 eV.
Neither provides a populated-state distribution across the D04 solar spectrum.
They are recorded here as follow-up references. The r11a numerical addition
for noble-gas product energies is the Laurell argon workbook.
[CDS J/ApJS/210/12](https://cdsarc.cds.unistra.fr/viz-bin/cat/J/ApJS/210/12)
[Heimann et al.](https://doi.org/10.1103/PhysRevA.34.3782)

The full C++/Python suite on the installed r11a lock reported **386 PASS,
0 FAIL, 0 ERROR, 0 SKIPPED** in 107.2 seconds. All 11 new source-specific
checks passed. The results are in `final/shared-tests/result.json` and
`final/nightly-core/artifacts/CORE-CHEMICAL-DATABASE-001/shared-protocol/result.json`
under the artifact root above.

## 5.54 What other exosphere models use for M3 chemistry — 2026-11-24 UTC

The [model comparison](../utility/chemical_database/NATIVE_M3_MODEL_PRACTICE_SOURCES.md)
links the original lunar, cometary, heliospheric and Mercury studies and the
specific code paths behind the approximations. Three new source packets record
the Poppe–Halekas–Harada lunar 1 au frequencies, water breakup choices, and
OH return and release models. The original licensed Poppe file contains 25
species. Its five-row projection preserves 45 numerical inputs verbatim.
For Na, its median photon rate corresponds to 36.6 h, close to the 36 h
value adopted by the 2013 lunar AMPS paper. Both draw on the Huebner family
of compilations. [Poppe dataset](https://doi.org/10.5281/zenodo.6994648)
[Tenishev et al.](https://doi.org/10.1016/j.icarus.2013.08.021)

Hartle et al. also modeled lunar pickup from a cold stationary neutral gas,
combining photoionization, charge exchange and electron impact as scalar
sources. Shen et al. later used the zero-speed approximation with measured
fields for selected ARTEMIS events. The Mercury hybrid model of Paral et al.
instead launches Na⁺ from the cell's neutral velocity distribution. A comparison
of stationary and moving-parent initial conditions with otherwise identical
forcing can show how much this choice affects the pickup spectrum. [Hartle et al., section 2](https://doi.org/10.1029/2011JA016859)
[Shen et al.](https://doi.org/10.3847/1538-4357/ad40a4)
[Paral et al., paragraph 12](https://doi.org/10.1029/2010GL044413)

For I22, the water breakup models supply wavelength-dependent branch choices
or neutral H/OH speeds. The surface models separately prescribe reactive
sticking, residence or a launch distribution. For I23, the WarpX
electron-impact sampler supplies a conservation benchmark. For I24, an
environmental total Na ion-production rate and an isotropic sensitivity case
still leave the state-selected 4 keV/u recoil unqualified. Broadband populated
ion states, joint charged-water and electron-impact daughter distributions,
and the incident OH sticking/accommodation law on regolith remain open.
The source changes introduce no production defaults or test ID. Formal
I21–I24 remain **4 SKIPPED / NOT QUALIFIED** with I13 and M2 still open.

The [I23 runner notes](test/integration-tests/I23/README.md) additionally
identify the published Turner helium benchmark's equal electron energy
partition and the incident-energy coverage of the Peterson et al.
secondary-electron measurements.

The working code revision is `e5f2978ac0576f42dc510e0da110dcbe0029d69a`.
The roadmap hashes and scientific qualification order in §5.53 still apply.
This work uses the source-tree early-development exception for a shared
catalogue. The reviewed database lock SHA-256 is
`05a4942e7ba115f7c51743fc1cfa6bfa3de4c531d2b5cc7bda3e7d6b8b44c879`.
The build, checksums, source-specific transcripts and qualification results
are under `test_output/chemical_database/m3-model-practice-mKGmz860`.

```bash
w=test_output/chemical_database/m3-model-practice-mKGmz860
db=/data/vtenishe/CHEMICAL_DATABASE/snapshots/amps-chemistry-20261124-r12
replay=$(mktemp -d "$w/replay-XXXXXXXX")
cp "$w/reviewed/build_release.py" "$replay/"
python3 -B "$replay/build_release.py" \
  test_output/chemical_database/m3-pickup-precedent-LUM6SwNP/final/source_packet.json \
  test_output/chemical_database/lunar-oh-precedent-twDfB4/packet/source_packet.json \
  test_output/chemical_database/water-precedent-Pcq5D0ih/packet/source_packet.json
bash "$w/test-release.sh" "$w/installed-tests-repeat" "$db"
python3 -B src/test/run_tests.py --config test/qualification-machine/passion.json \
  --campaign nightly --test CORE-CHEMICAL-DATABASE-001 --output-dir "$w/nightly-repeat"
```

The build helper and its complete command in `reviewed/build-command.json`
record this publication. The release has 64 datasets, 986 source
entries, 162 tables with 134 available through the catalogue, and 36 explicit
source packets. All 974 previous source entries and 161 previous tables
passed the preservation comparison. Database software verification uses
Python 3.14 and g++ 9.4 with the normal native C++ catalogue probe.

The full suite on the installed r12 lock recorded **394 PASS**, including
all eight new real-source tests, in 106.9 seconds. The nightly core database
protocol recorded **309 PASS** in 23.4 seconds. Both runs had zero failures,
errors or skips. The installer checked 542 files and verified 886 object
references. A fresh installed CLI export of the five-row table reproduced
SHA-256 `76fb5cdb5181a5a4a6ceb0e9fff40da28fb3069affe107e440644863ed083b61`.
The default database lock still has SHA-256
`fdf5311274871d7abe884b3d8aec163bfdc52e8c3295adb09e4c5cd42813d945`.
[Installed suite](../test_output/chemical_database/m3-model-practice-mKGmz860/reviewed/installed-tests/result.json)
[Nightly protocol](../test_output/chemical_database/m3-model-practice-mKGmz860/nightly-core/artifacts/CORE-CHEMICAL-DATABASE-001/shared-protocol/result.json)
[Receipt](../test_output/chemical_database/m3-model-practice-mKGmz860/reviewed/installation.json)

## 5.55 Lunar model inputs and AMPS reference solutions — 2026-11-24 UTC

The [reference comparison](test/M3_LUNAR_MODEL_REFERENCE_CROSSWALK.md)
records how to use output from a reviewed AMPS run, the same model at a
different MPI layout, a sibling exosphere model, and an independent physical
calculation. Rechecking the retained four-rank run against the actual I04
baseline gave **PASS, zero diff bytes**. Native fates and times from one and
eight MPI ranks matched for **384/384** parent/refinement pairs in the I21
varying-exposure case. Both checks used previously executed runs. The existing
analytical event comparison provides a separate check of numerical accuracy.

Inspecting the Europa input, its shared transport callback and native Huebner
arrays established that its three charged water breakup channels match all
12 frozen D04 quiet/active rate and mean-energy values. The arrays describe
the available energy for each channel. The transport callback samples each
species separately. The historical Moon electron-impact estimate was found in
an archived checkout whose generated configuration disabled the process.
The lunar source papers and native Julia/Fortran programs supply independent
model assumptions about neutral loss and water retention. They leave the
state-resolved and correlated product distributions in §5.54 open.

The shared r13 release adds seven source entries in three packets. The
reviewed publication ties the Peschel PDF and three-row CSV to the article's
CC BY 4.0 notice. The Schörghofer and EXOSPHID packets preserve model
parameters with the version and operating conditions of their origin.
Read the [database examples](../utility/chemical_database/DATABASE_README.md#19-lunar-neutral-losses-water-coverage-and-noble-gas-product-assumptions)
for the lock and public queries. All 36 prior packet digests, 986 source
entries and 162 tables passed the preservation check. The new catalogue
has 67 datasets, 993 sources and 163 tables, including 135 available tables.
Lock SHA-256 `6890f21220d6065a19e10f19469210cad66702babcddc083854a6cff7617664a`.
D04's source CSV was checked against the canonical package manifest,
SHA-256 `3aa54915a08258371dbe975949e6afb2ab74618164fadca8bceebe876c926847`.

AMPS revision `e5f2978ac0576f42dc510e0da110dcbe0029d69a` with working-tree
changes. The two roadmap hashes and their respective scopes are recorded
in §5.53. The only build needed was the shared Python/C++ database suite.
The original three model packets, explicit release command and independently
repeatable comparison commands are recorded with the artifacts. I21–I24
retain their formal **4 SKIPPED** status.

```bash
w=test_output/chemical_database/lunar-source-comparison-pTDdeg8Y
db=/data/vtenishe/CHEMICAL_DATABASE/snapshots/amps-chemistry-20261124-r13
python3 -B "$w/compare_amps_mpi_outputs.py"
python3 -B "$w/check_europa_d04.py"
python3 -B srcMoon/test/integration-tests/I04/check_reference.py \
  --run-dir run_test_Moon --diff-file "$w/i04-recheck/test_Moon.diff" --diffnum share/Scripts/DiffNum.pl
bash "$w/test-release.sh" "$w/repeat-suite" "$db"
python3 -B src/test/run_tests.py --config test/qualification-machine/passion.json \
  --campaign nightly --test CORE-CHEMICAL-DATABASE-001 --output-dir "$w/repeat-nightly"
```

The full staged r13 suite reported **400 PASS** in 111.9 s, including
six new source-backed tests. The nightly database protocol reported
**309 PASS** in 23.9 s. Both reported zero failures, errors and skips.
Installation verified 560 files and 901 object references. Installed CLI
and API checks reproduced the new table digest and seven source records.
The default database lock remains the one recorded in §5.54. No new
AMPS particle generation or observational scoring was performed.
[Full suite](../test_output/chemical_database/lunar-source-comparison-pTDdeg8Y/stage-tests/result.json)
[Nightly](../test_output/chemical_database/lunar-source-comparison-pTDdeg8Y/nightly/artifacts/CORE-CHEMICAL-DATABASE-001/shared-protocol/result.json)
[Installed checks](../test_output/chemical_database/lunar-source-comparison-pTDdeg8Y/installed-checks.json)

## 5.56 Sodium pickup and ground-product collision reference models — 2026-11-24 UTC

The user approved publishing two distinct reference cases in the shared
database. The review used the existing Poppe deposit and §2.2 of its paper,
and the authenticated Shingal 1986 Table 4 benchmark. The r14 publication
replays all 39 r13 source packets from their SHA-256-pinned originals. It adds
two software model records and one small source-selector table. Fourteen
previous model definitions and the payload digests of 162 other tables were
compared to r13. This work was performed against AMPS revision
`e5f2978ac0576f42dc510e0da110dcbe0029d69a`.

* `poppe2022-na-pickup-median-v1` selects the original Na median photo,
  proton-exchange and electron-impact frequencies at 1 au. Charge exchange
  amounts to 9.4064945% of the three-process total under those conditions.
  The published pickup source assigns zero newborn ion velocity in the lunar frame.
* `hplus-na3s-ground-capture-4kev-u-isotropic-v1` selects the theory partial
  2.2 × 10⁻¹⁷ cm² for H⁺ + Na(3s) → H(1s) + Na⁺(2p⁶) at 4.0 keV/u,
  the approved isotropic direction in the incoming pair COM, and the shared
  conservative two-body recoil definition. Other projectile energies and the
  uncertainty on the theory point are unspecified. Actual Na⁺ energy spectra
  and return fractions require transport runs.

The source-tree roadmap controls early M3 software work. The common R3.4
criteria still control formal I24. This data-access work executed no linked
lunar transport or observational comparison. For the dated differences between
the roadmaps and their SHA-256s see [roadmap reconciliation](test/ROADMAP_RECONCILIATION.md).
[Source selections and API commands](../utility/chemical_database/DATABASE_README.md#sodium-pickup-and-binary-charge-exchange-r14)

The full test protocol detected an edit to the citation list during the first
development run. It recorded **405 PASS and one infrastructure ERROR** because
the source identity changed during execution. The final snapshot includes the
Hartle and Shen pickup references. The repeats use a new publication directory
and a new test result. [Initial run](../test_output/chemical_database/sodium-scenarios-nB2oDjV5/full/result.json)

The following commands reproduce the release checks. Run from the AMPS root and
use fresh output directories. The recorded build command includes all 39 packet
paths. The inherited suite script supplies the original-source paths used on
this machine.

```bash
w=test_output/chemical_database/sodium-scenarios-nB2oDjV5
db=/data/vtenishe/CHEMICAL_DATABASE/snapshots/amps-chemistry-20261124-r14
scratch=$(mktemp -d)
AMPS_CHECK_NA_REFERENCE_MODELS=1 bash test_output/chemical_database/lunar-source-comparison-pTDdeg8Y/test-release.sh "$scratch/full" "$db"
python3 -B src/test/run_tests.py --config test/qualification-machine/passion.json --campaign nightly --test CORE-CHEMICAL-DATABASE-001 --output-dir "$scratch/nightly"
PYTHONPATH="$db/python" python3 -B "$w/check_publication.py" "$db"
```

Final r14 results were **405 PASS** in the full shared database suite and
**310 PASS** in the core nightly database protocol. Each reported zero FAIL,
ERROR and SKIPPED. The installer verified 562 files and 905 object references.
Rerunning the published CLI and Python example from the installed snapshot
reproduced the native medians, theory point and model IDs. The r14 lock SHA-256
is `9fd6272f652fa14e33517f1354d9dc3ee6bf47948a07b3438b3b9f4c7b14f0d8`.
The canonical top-level lock still identifies the original release.
[Full suite](../test_output/chemical_database/sodium-scenarios-nB2oDjV5/release-tests/result.json) ·
[Nightly protocol](../test_output/chemical_database/sodium-scenarios-nB2oDjV5/release-nightly/artifacts/CORE-CHEMICAL-DATABASE-001/shared-protocol/result.json) ·
[Installed snapshot checks](../test_output/chemical_database/sodium-scenarios-nB2oDjV5/installed-checks.json) ·
[Installation receipt](../test_output/chemical_database/sodium-scenarios-nB2oDjV5/installation.json) ·
[Exact publication command](../test_output/chemical_database/sodium-scenarios-nB2oDjV5/release-command.json)

## 5.57 Linked M3 chemistry development — 2026-11-24 UTC

Source revision `e5f2978ac0576f42dc510e0da110dcbe0029d69a` plus the
working sources recorded by the clean producer. Source roadmap SHA-256
`712e1152bf50af6d9ee1dad8f8a585809ca8b92cf42ac6fc7a009f6c0253ee17`
governs the early-development and one-hour development exception. Data-root
roadmap SHA-256
`d10eca1c539f9cf82940236d9e1aa36af665d43e82ba475925950d582fbf49c4`
governs the common scientific criteria and formal predecessor order.
The global I21–I24 run reported **4 SKIPPED** with zero FAIL or ERROR.
Its manifests identify the missing same-campaign I13–I20 results and
capabilities for each test. [Formal aggregate](../test_output/srcMoon/m3-production-aPrmWNe5/formal/run_all_summary.json)

* **I21, Na photoionization.** The fixed illuminated cohort now compares
  the exact retained native time with an independent 80-digit first-hit
  calculation. Six original launches, with 43 conversions in each,
  passed 18 checks and failed the strict refinement check. Maximum
  normalized errors at 8, 16 and 32 steps were 4.957, 4.957 and 4.972
  × 10⁻¹⁷. The calculation separates the rounded initial optical depth
  from time integration and CSV projection. The transport contribution
  decreases on refinement. The existing varying-exposure case passed its
  per-cohort tests. The new zero-hazard case checks that the finite
  particle-to-Sun ray intersects the Moon independently of the production
  visibility flag. [Numerical method and results](test/sodium_cohort/NUMERICAL_README.md)
* **I22, photodissociation and photoionization.** The new parent-rest
  photon/three-body solver conserves momentum and energy for an explicit ion,
  electron and neutral with source-supplied states, energy partition and
  directions. Its controls passed 456 checks. Existing two-body photon
  controls passed 981 and D04 atomic-source preparation passed 70.
  Complete H₂O/OH source populations and correlated daughter distributions,
  normal molecular dispatch and U24 remain open. [Charged breakup results](../test_output/srcMoon/m3-i22-i23-E1QrJUHL/charged-photon/result.json)
* **I23, electron impact.** The stationary-He pilot passed 25 checks with
  the authenticated D05 coefficients and independent Kim fit. At 10 eV,
  the 100–2000 eV incident-energy interval in the Peterson measurement
  covers 0.5321% of the ionizing events predicted by that fit. It cannot
  establish the error in the bulk of the pilot distribution. The linked
  approved reduced He reservoir still needs a source-to-state threshold
  crosswalk, coverage and uncertainty of the incoming collision-weighted
  moment, a recoil/observable error bound, and normal bath-exposure coupling.
  Fully differential outgoing electron spectra are optional for this
  reduced closure. [Pilot result](../test_output/srcMoon/m3-i22-i23-E1QrJUHL/helium-final/result.json)
* **I24, charge exchange.** The production provider reads the r14 0.22 Å²
  partial cross section for H⁺ + Na(3s) → H(1s) + Na⁺(2p⁶) at 4.0 keV/u.
  It uses the approved isotropic pair-COM approximation and source-based
  mass and formation energies. In the prescribed homogeneous beam the direct
  rate was 4.8321268971 × 10⁻⁹ s⁻¹. Of 32,768 independent parents, 20,775
  reacted and 11,993 survived, with identical histories at one and eight
  ranks and a lifetime error of 0.51%. All 35 ensemble checks passed.
  The broader charge-exchange family still requires state-resolved products
  over the selected projectile-energy support and applicable bath
  conditions. [Ensemble results](../test_output/srcMoon/i24-sodium-audit-J0bii1U5/final/I24/sodium-ground-capture/result.json) ·
  [Linked run instructions](test/integration-tests/I24/README.md)

The first Na/Na⁺/H full-particle run stopped in the production
`localTimeStep` callback because it recognized no neutral H. The corrected
startup requires an explicit cell-crossing speed and H/Na reference-weight
ratio when H is compiled. The source-bound build uses the normal SPICE
configure/build path with optimization and debug information, individual
particle weights and a zero-acceleration test equation. The default
four-rank `make -j test_Moon` finished with a **0-byte test_Moon.diff**
against the reviewed I04 reference. [Original startup log](../test_output/srcMoon/m3-production-aPrmWNe5/pic/first-1mpi.log) ·
[Fresh build](../test_output/srcMoon/m3-production-aPrmWNe5/tracked-h-build-run/result.json) ·
[Nightly regression](../test_output/srcMoon/m3-production-aPrmWNe5/tracked-h-nightly.log)

The corrected H build then completed the full I24 PIC matrix. First- and
second-order movers at one and eight ranks all returned identical live
phase-space records. All **71 individual checks passed**. For a Na packet
representing 2.625 atoms, Na⁺ and H each carried 2.625 at reference weights
56 and 0.875. The bath lost 2.625 protons. With 0.0625 s of remaining flight,
the relative momentum and energy errors were 1.81 × 10⁻¹⁶ and 1.49 × 10⁻¹⁶.
The first scorer incorrectly interpreted the two integer charge changes
per *reaction* as packet-weighted moments. The repeated run checked both
stoichiometry and represented charge against the native packets and the
signed bath record. Formal I24 remains SKIPPED. [Native PIC matrix](../test_output/srcMoon/m3-production-aPrmWNe5/verified/I24/sodium-pic-ground-capture/result.json) ·
[Example displacement figure](../test_output/srcMoon/m3-production-aPrmWNe5/verified/I24/sodium-pic-ground-capture/second-8mpi/continuation.png) ·
[Ensemble and formal result](../test_output/srcMoon/m3-production-aPrmWNe5/verified/I24/result.json)

The maintained Moon configuration runner passed **54** checks, including
14 for H speed and packet weights. The source-bound build protocol passed
11 checks. The full shared database rerun with the real source selectors
and r14 lock passed **405** with zero FAIL, ERROR or SKIPPED. The release
lock SHA-256 is `9fd6272f652fa14e33517f1354d9dc3ee6bf47948a07b3438b3b9f4c7b14f0d8`.
The Na/H state inventory SHA-256 is `7c255c05398c9dbaf4d55a8b17550576884dd0c5710773d3f8bfb20cde6f7fa1`.
[Moon input controls](../test_output/srcMoon/i24-sodium-audit-J0bii1U5/hydrogen-input/result.json) ·
[Shared database results](../test_output/chemical_database/r14-all-source-selectors-4r6ikioh/suite/result.json)

Reproduce the conditional I24 case with the following commands after creating
the named three-species producer. The input, native command lines, build
identity, source hashes, SPICE kernels, root random key and per-rank logs
are retained in the individual manifests. The four-rank nightly uses its
separately reviewed regression configuration.

```sh
db=/data/vtenishe/CHEMICAL_DATABASE/snapshots/amps-chemistry-20261124-r14
pkg=test_output/srcMoon/m3-production-aPrmWNe5
/usr/bin/python3 -B srcMoon/test/run_all_test.py --only I24 --skip-build \
  --mpi-ranks 8 --campaign-mode regression --test-runtime-seconds 1600 \
  --i24-sodium-lock "$db/database-lock.json" \
  --i24-pic-refinement --i24-pic-producer "$pkg/closure-capture-build/normal_build_evidence.json" \
  --output-dir "$(mktemp -d)"
make -j test_Moon
```

The [I21 fixture](test/integration-tests/I21/README.md#reproduce-the-lunar-shadow-fixture)
and [I24 builder](test/source_qualification_build/README.md#build-the-i24-sodium-capture-case)
give the commands for rebuilding inputs on another working checkout.

A second, clean Na/Na⁺ production build of the current working sources
replayed the dark I21 fixture. Each of six runs (1/8 MPI ranks, 8/16/32
steps) retained all 128 originals and produced zero removals and zero bath
entries. The independent finite segment passed within 9.1431 × 10⁻¹⁰ m of
the Moon centre for a 1,737,100 m radius. All **25 checks passed**, including
exact zero numerical error in both rank layouts. The build took 153 s and
the matrix plus audit took 179 s. [Production build](../test_output/srcMoon/m3-production-aPrmWNe5/dark-current/native-build/result.json) ·
[Dark cohort audit](../test_output/srcMoon/m3-production-aPrmWNe5/dark-current/aggregate/I21/homogeneous/independent-audit/result.json) ·
[One shadow geometry receipt](../test_output/srcMoon/m3-production-aPrmWNe5/dark-current/aggregate/I21/homogeneous/independent-audit/lunar_shadow-05.json)

The registered core nightly selection `CORE-CHEMICAL-DATABASE-001` also
passed. Its child suite reported **310 PASS, 0 FAIL, 0 ERROR, 0 SKIPPED**.
[Core aggregate](../test_output/srcMoon/m3-production-aPrmWNe5/nightly-database/amps-core-tests.json) ·
[Database checks](../test_output/srcMoon/m3-production-aPrmWNe5/nightly-database/artifacts/CORE-CHEMICAL-DATABASE-001/shared-protocol/result.json)

The r14 review identified a reproducible conditional H₂O/OH residence-time
comparison against Laferrière's prescribed smooth-surface scenario. The
actual Arrhenius surface kernel passed nine numerical checks. Boundary
conditions and source data for a complete regolith sticking/accommodation
and migration calculation are still required for U24. [U24 scope and
results](../test_output/srcMoon/m3-i22-i23-E1QrJUHL/u24-r14-final/U24/result.json)


In the first I24 common-time refinement, the four-step case stopped at a
valid interior event. A constant positive endpoint rate had been evaluated
as the sum of two rounded interpolation terms. The independently
re-evaluated physical rate therefore differed in its last bit. The production
piecewise-hazard evaluator now returns the identical endpoint rate for an
exactly constant segment. The strict equality test remains in the event
controller. The focused regression reproduced 5.18750000000000000043
from an exact 5.1875 before the fix; the corrected local suite passed 644
checks. [Failed native launch](../test_output/srcMoon/i24-sodium-audit-J0bii1U5/pic-refinement/second-steps4-1mpi.log) ·
[Constant-rate unit result](../test_output/srcMoon/i24-sodium-audit-J0bii1U5/constant-rate-unit/result.json)

The I21 precision experiment resolved the illuminated cohort discrepancy with
an opt-in numerical change. The source draws exactly the same uniform
variates. The new mode evaluates −ln(U) in extended precision and preserves
that result in the clock's existing high and low components at birth,
migration and daughter initialization. The existing default and its reported
refinement failure are retained. The predeclared experiment produced
normalized maxima 1.3722 × 10⁻¹⁷, 6.8884 × 10⁻¹⁸ and 3.3138 × 10⁻¹⁸
at 8, 16 and 32 steps, identically on one and eight ranks. The strict
decrease and 0.002 fine-level threshold both passed. The keyed physical
cohort had 43 removals in 128 parents per launch. Select this mode explicitly in the sodium cohort input. A host without
extended precision cannot execute that selection. [Audited six-run experiment](../test_output/srcMoon/m3-production-aPrmWNe5/i21-opt-in-precision/native-matrix/independent-audit/result.json) ·
[Mode selection and scorer](test/sodium_cohort/NUMERICAL_README.md)

### Final integrated charge-exchange and nightly results

After the constant-rate and optical-clock changes, a clean SPICE build
completed in 178 s. The registered I24 development run finished in 211 s
and reported both supporting children PASS. The ensemble passed 35 checks.
The live PIC matrix passed **227 checks in 12 cases**, using both movers,
1 and 8 ranks, 1/2/4 common steps and 1/½/¼ reference-weight scales.
The physical beam, keyed original, step horizon and 2.625 represented atoms
were identical across the matrix. The maximum absolute event-epoch error
was 4.47 × 10⁻⁸ s (bound 9.54 × 10⁻⁷ s), the largest component position
error was 0.00305 m, and the largest relative momentum and energy errors
were 1.81 × 10⁻¹⁶ and 1.49 × 10⁻¹⁶. The registered owner retains formal I24
as SKIPPED. [Global I24 result](../test_output/srcMoon/m3-production-aPrmWNe5/final-i24/I24/result.json) ·
[12-case table](../test_output/srcMoon/m3-production-aPrmWNe5/final-i24/I24/sodium-pic-ground-capture/refinement.csv) ·
[Eight-rank four-step comparison](../test_output/srcMoon/m3-production-aPrmWNe5/final-i24/I24/sodium-pic-ground-capture/second-steps4-8mpi/continuation.png) ·
[Source-bound build](../test_output/srcMoon/m3-production-aPrmWNe5/closure-capture-build/normal_build_evidence.json)

The final default `make -j test_Moon` returned the reviewed I04 digest
`ed875cefc604e4ae3580fb89b8c528146929b6b9be099b41b61697faee7f0080` and
a zero-byte diff. The separate global qualification-mode check recorded
**U14 PASS, I21/I22/I23/I24 SKIPPED**. Its per-ID results list the required
predecessors and outstanding family or forcing evidence. [Nightly reference receipt](../test_output/srcMoon/m3-production-aPrmWNe5/closure-nightly-reference.json) ·
[Qualification-mode summary](../test_output/srcMoon/m3-production-aPrmWNe5/qualification/run_all_summary.json)

On the integrated current tree the Na/Na⁺ clean build took 177 s and the
six illuminated native launches plus audit took 250 s. All six launches and
all **19 independent checks passed** with the same three timing errors given
above. The historical default-mode comparison still records its actual
failure. Larger predeclared statistical families, physical source-model
qualification and the same-campaign predecessors remain open.
[Current-source audit](../test_output/srcMoon/m3-production-aPrmWNe5/i21-current-precision/native-matrix/independent-audit/result.json) ·
[Exact inputs, build and run evidence](../test_output/srcMoon/m3-production-aPrmWNe5/i21-current-precision/README.md)

U15's registered owner passed all 27 supporting groups when run from its
matched source checkout with the authentic dependencies. Formal U15 remains
SKIPPED for incomplete product mapping/dispatch and U25. U16 passed its
302 checks. U17 passed four supporting services and remains SKIPPED for
the broader charge-exchange product/branch contract. These results identify
which local gates still block the formal I-series runs. [U15](../test_output/srcMoon/m3-i22-i23-E1QrJUHL/u15-complete-packages/U15/result.json) ·
[U16](../test_output/srcMoon/m3-i22-i23-E1QrJUHL/u16/U16/result.json) ·
[U17](../test_output/srcMoon/m3-production-aPrmWNe5/local-u14-u17/U17/result.json)

Two compatibility campaigns then ran against the final same Na/Na⁺ producer.
The checked-in shadow fixture passed **25/25** across 1/8 MPI and 8/16/32
steps. The moving, illuminated case using the new precision mode passed
**19/19** over the same matrix. At +256 m/s along LSO x, all six launches
produced 51 removals and 77 survivors. Normalized event-time errors fell
from 9.96 × 10⁻¹¹ to 2.47 × 10⁻¹¹ to 6.14 × 10⁻¹². The maximum
position error fell from 1.67 mm to 0.103 mm. The latter campaign took
169 s and peaked at 28.4 GB. [Checked-in shadow replay](../test_output/srcMoon/m3-production-aPrmWNe5/final-shadow/aggregate/I21/homogeneous/independent-audit/result.json) ·
[Moving cohort summary](../test_output/srcMoon/m3-production-aPrmWNe5/i21-current-trajectory/summary.json) ·
[Timing comparison figure](../test_output/srcMoon/m3-production-aPrmWNe5/i21-current-trajectory/campaign/I21/trajectory/independent-audit/numerical_comparison.png)

## 5.58 Frozen M3 development family — 2026-11-24 UTC

The new `run_m3_campaign.py` pins an ordered set of source-bound cases,
provides one shared wall, RSS and artifact budget, and retains the existing
individual physics runners. It indexes case and independent check outcomes
in one CSV. Each I21 cohort additionally compares original parent fates and
channels across one and eight MPI ranks at every matched timestep. Native
first-hit and conservation scoring still uses independent source references.
The four formal results are preserved with their own prerequisite records.
[Campaign contract](test/qualification/CAMPAIGN_CONTRACT.md) ·
[Plan and execution commands](test/qualification/README.md#frozen-chemistry-family)

The active source revision is `74b15db8ac95dc72f355d065b83fe950e215ed6c`
plus the new Python runner and audit changes. The source roadmap identity
`712e1152bf50af6d9ee1dad8f8a585809ca8b92cf42ac6fc7a009f6c0253ee17`
controls the early-development exception and resource budget. The canonical
data roadmap identity
`d10eca1c539f9cf82940236d9e1aa36af665d43e82ba475925950d582fbf49c4`
controls the common science criteria and formal predecessor order.
No raw files or roadmap thresholds changed.

Two private source-matched SPICE builds passed, each at about 173 s, with
`make -j16`, optimization and debug information after normal clean
configuration. The standard nightly `make -j test_Moon` passed with a
zero-byte diff on the current source revision. The family-runner controls
passed eight tests, formal owner controls passed 11, the independent
numerical suite passed 42 and the statistical/qualification suite passed
31. [Build preparation and commands](../test_output/srcMoon/m3-family-t1soz40x/build.py) ·
[Nightly receipt](../test_output/srcMoon/m3-family-t1soz40x/nightly-reference.json) ·
[Family runner controls](../test_output/srcMoon/m3-family-t1soz40x/family-controls.log)

The family uses the previously declared Na/Na⁺ quiet-Sun fixtures, the
immutable D04 rate table SHA-256
`3aa54915a08258371dbe975949e6afb2ab74618164fadca8bceebe876c926847`,
the shared r14 lock SHA-256
`9fd6272f652fa14e33517f1354d9dc3ee6bf47948a07b3438b3b9f4c7b14f0d8`
and the NIST state inventory SHA-256
`7c255c05398c9dbaf4d55a8b17550576884dd0c5710773d3f8bfb20cde6f7fa1`.
The I21 source geometry, statistical keys and scientific case membership
were held fixed when the plans were bound to the new builds. The formal
I21–I24 owners still report SKIPPED for the recorded physics and I13–I20
prerequisites. [Frozen family](../test_output/srcMoon/m3-family-t1soz40x/campaign.json) ·
[Exact invocation](../test_output/srcMoon/m3-family-t1soz40x/campaign-command.json)

The complete frozen family reported **5 PASS, 0 FAIL, 0 ERROR, 0 SKIPPED**
and **334 independent checks passed** in **772.6 s**. The peak supervised
process-tree RSS was 42.9 GB, and the final family output size was 4.45 GB.
The separate formal readiness summary reported **4 SKIPPED**. Counts by
supporting case were 22 illuminated I21, 28 shadowed I21, 22 moving I21,
35 I24 ensemble and 227 I24 live particle/refinement checks. I21 matched all
128 original parent histories in each of nine 1/8-rank pairs, across the 8,
16 and 32-step cases. The maximum rank-to-rank serialized first-hit
difference was exactly zero. [Family result](../test_output/srcMoon/m3-family-t1soz40x/campaign/result.json) ·
[Comparison table](../test_output/srcMoon/m3-family-t1soz40x/campaign/comparison.csv) ·
[Run manifest](../test_output/srcMoon/m3-family-t1soz40x/campaign/run_manifest.json)

A separate mixed-input check reported SKIPPED for a missing I21 input and
ERROR for an intentionally changed SHA-256 pin, and then passed all 35 checks
of the real sodium capture ensemble. Neither failed input was launched.
Follow-up reporting changes make `progress.json` finish at its terminal
status and reject malformed or inconsistent independent child receipts.
All 334 retained check records replayed identically against the final report
schema, and all nine final family controls passed. [First mixed result](../test_output/srcMoon/m3-family-t1soz40x/mixed/result.json) ·
[Retained result replay](../test_output/srcMoon/m3-family-t1soz40x/final-check-schema-replay.json) ·
[Final controls](../test_output/srcMoon/m3-family-t1soz40x/family-controls-final.log)

The final reporting rerun exercised both real I24 children after the
intentionally unavailable and changed-pin cases. It retained **2 PASS,
1 SKIPPED, 1 ERROR, 0 FAIL** and all **262 I24 checks passed**. The terminal
`progress.json` and `result.json` are byte-identical. The global U16 owner
passes and now executes the nine new controls in its existing 101-test
orchestration battery. [Final reporting rerun](../test_output/srcMoon/m3-family-t1soz40x/final-reporting-replay/result.json) ·
[Global U16](../test_output/srcMoon/m3-family-t1soz40x/global-u16/run_all_summary.json) ·
[Orchestration log](../test_output/srcMoon/m3-family-t1soz40x/global-u16/U16/orchestration-controls.log)

## 5.59 I23 prescribed electron-bath coupling — 2026-11-25 UTC

Source revision `74b15db8ac95dc72f355d065b83fe950e215ed6c` plus the working
tree. The source roadmap SHA-256 `712e1152…` governs the approved early-M3
development, and the data-root roadmap `d10eca1c…` governs the common
criteria and predecessor order. No raw data, roadmap threshold, test ID or
reviewed reference changed. Evidence is under
`test_output/srcMoon/i23-bath-uGAYf6qI/`.

**User decisions.** Before implementation the user approved:

- the controlled bath: `k_B T_e = 10 eV`, `n_e ∈ {0, 1e6, 1e7} m⁻³`,
  isotropic and zero-drift, at rest in frozen LSO, uniform and constant,
  stationary ⁴He only (this resolves decision 3 of the 2026-11-24
  follow-up);
- production reservoir products with a narrow observable scope: heavy
  budgets and chronology only;
- the threshold/state crosswalk (all products He⁺(1s), excited yield
  unknown) and the Kim − Voronov slope spread (0.120) as the conditional-mean
  bound;
- conservative upper-end recoil tolerances (≤1095 m/s RMS, ≤0.0249 eV) and
  He P02 limits equal to the Na values;
- explicit He/He⁺ startup keys.

[Approval record](test/integration-tests/I23/HELIUM_ADMISSION_20261125.md)

**Production changes.**

- New `MoonHeliumImpactCohort.{h,cpp}`, linked through the makefile, is an
  explicit owner registered before particle allocation. It checks the
  stationary, force-free domain on every accepted segment, uses the
  production `HeliumElectronImpact::Compose`, and refuses invalid
  registrations before startup.
- `MoonInput` and `main_lib.cpp` gain the He⁺ timestep/weight pair and an
  optional He weight ratio. A He⁺ build previously stopped in
  `localTimeStep` ("species is not recognized"); a source-less He build
  stopped dividing by a zero He production rate. Legacy behaviour is
  unchanged when the keys are absent.

**Admission evidence.** `prepare_helium_admission.py` derives the crosswalk
(threshold difference 0.012611 eV plus a binary64/extended rounding
allowance of about 6 × 10⁻³⁴ J) and the model-error file (slope error
0.12003; upper RMS He⁺ speed bound 1094.7 m/s; kinetic bound 0.02486 eV).
Both files are frozen and hash-pinned, and every run regenerates them byte
for byte. The first ensemble attempt was ERROR: the stored binary64 24.6 eV
fit threshold exceeded an extended-only rounding allowance. The allowance
was corrected to cover the binary64 representation. The second attempt
FAILed on an incorrect audit expectation (per-event rather than
represented-population ledger totals); the audit now scores the totals, as
the I24 audit does. Production code was not changed for it. Both attempts
are retained.

| Evidence | Result |
| --- | --- |
| Ensemble (`ensemble-3`) | PASS, 165 checks in 73 s. Coefficient `7.7606300991732e-16 m³/s`, matching Decimal Voronov to `1.2e-16`; 2% gate met. At `1e7 m⁻³`, 20,628 of 32,768 parents ionized at one optical depth; endpoint lifetime error 0.71% (gate 3%). Event-time error ≤ `2.3e-16` of the horizon. Zero density keeps every parent. Every event matched the reservoir: +1 e, +mₑ, −e, −I, zero momentum, `E_in = 45.418 eV`, `E_out = 20.831 eV`. Identical 1/8-rank histories. The `1e6` case reuses the keyed draws and tests linear scaling only. |
| Linked PIC (`pic`) | PASS, 209 checks in 14 launches (176 s): both movers, 1/8 ranks, zero density, 1/2/4 steps, ½/¼ weights. The independent first event at `9.5855e7 s` was reproduced with zero ET error (bound `9.5e-7 s`). He⁺ held 2.625 represented atoms at rest. The bath recorded 2.625 incoming and 5.25 outgoing electrons; ledger closure met. Normal startup installed the He/He⁺ weights and He⁺ timestep from the new keys. |
| Rebuilt producers | Na/Na⁺ (165 s), Na/Na⁺/H (162 s) and Na/He/He⁺ (153 s) all PASS through the clean source-bound Config/make. |
| Family attempt 1 (`campaign`) | ERROR, preserved. Six cases PASS (499 checks). The I23 PIC probe failed to compile after the native refusal controls were added (C++ most-vexing-parse). |
| Family attempt 2 (`campaign-2`) | PASS: 7 cases, **722 independent checks**, 0 FAIL/ERROR/SKIPPED, 1239 s, peak 42.9 GB, 4.7 GB output. Na cases 22/28/22/35/227 checks; I23 ensemble 165 and PIC 223, including six native registration refusals. Formal I21–I24 4 SKIPPED. |
| Public aggregate (`formal-aggregate`) | `run_all_test.py --only I23` with the new flags: both supporting children PASS; I23 formal SKIPPED. |
| Nightly (disposable checkout) | `make -j test_Moon` zero-byte diff, reviewed I04 `ed875cef…`. |
| Controls | Orchestration battery 109 (7 new I23 admission controls), formal owners 14 (3 new), builder protocol 11, runtime configuration 74 (20 new He-key variants). The existing helium pilot passes 25 and reproduces the Kim comparison byte for byte. |

**Open limits.**

- Formal I23 needs same-campaign I13–I20, U16 and the frozen full family.
- The admission does not cover He⁺ velocity or energy observables, moving
  targets, drifting or structured baths, other temperatures, Ne/Ar/Na
  electron impact, excited-He⁺ yields, or measured R4.1/D06 forcing.
- The bath is not a lunar plasma model.
- Neutral He keeps its legacy hard-coded timestep ratio (0.1 × Na).
- A pre-existing `_NE_PLUS_SPEC_` weight-copy line in `main_lib.cpp` copies
  Ne to Ne instead of Ne⁺. No current build compiles Ne⁺; it was not changed.

[Contract](test/integration-tests/I23/ELECTRON_BATH_COUPLING_CONTRACT.md) ·
[Run commands](RUN.md#i23-prescribed-electron-bath-coupling)
