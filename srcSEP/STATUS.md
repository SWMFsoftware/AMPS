# srcSEP development status

Last updated: **2026-11-23**. This file tracks the active implementation of
both DOCX plans and their corrective proposal. **The full development and
validation program is not complete.** A passing software or isolated-physics
test does not qualify an unimplemented production path.

Update this file whenever a change is implemented, a verification gate changes,
a physical decision is clarified, or a blocker is discovered. Keep failed
attempts and their diagnosis visible until the corrected path has passed. The
primary requirement inventory is
[development_requirements.json](validation/development_requirements.json);
the original requirements remain in
[the development plan](improvement_and_validation_plan.docx) and
[the validation plan](srcSEP_Physics_Validation_Campaign_Implementation_Plan.docx).

## Stage 12 upstream physics decision, 2026-11-23

The user approved the prescribed outer-wind candidate for the new SWCME1D
`SHOCK_ONLY` material profile. It uses constant radial wind, inverse-square
mass density, an induction-compatible Parker field and explicit sustaining
force and work accounting. The approved equations and independent kinematic
checks are in [D10-1](validation/STAGE10_BACKEND_CAPABILITY_AUDIT.md#decision-records-and-frozen-acceptance).
The user also approved provider coverage of 0.30–1.20 AU and adiabatic expansion
with `gamma = 5/3` for the first controlled case. The 1 AU anchor uses 400 km/s
radial wind, 6 cm^-3 electron density, 120000 K electron and proton temperatures,
zero alpha abundance and 5 nT equatorial field magnitude. Test both polarities.
The independent reference values were calculated from the checked-in SWCME
constants. This closes the physics inputs for the first controlled upstream
case. `FULL_ICME` and the Stage 13
downstream volume need their separate closure decisions. The Stage 14 shared-core
background assembly can use its existing verification histories independently.

## Stage 12 upstream and Stage 14 core assembly, 2026-11-23

The authorized work covers the SWCME1D prescribed upstream and the shared
coronal ellipsoid/piston background assembly. The SWCME1D core now supplies its
analytic parcel deformation and Jacobian, charge-neutral e/p thermodynamics,
signed magnetic field, both Cartesian gradient tensors, sustaining force and
mechanical power. Its typed event uses the existing ballistic, DBM, and
data-driven front histories and local fast-shock solver. A separate detached
SEP adapter exports the actual canonical state and upstream/downstream surface
limits. Native field-line/source/restart admission still rejects that owner.
The legacy diagnostics and native AMPS timestep contract are retained.

The first strict adapter run passed MBC01–20. The global SWCME owner executed
all 14 bounded cases successfully, including ten independent upstream cases;
the original 127 dormant and 16 observational cases remained typed SKIP. The
same ten physics cases also passed strict warnings and ASan/UBSan. Nine
deliberate, independently compiled physical regressions each returned the
required FAIL. The Stage 12 AMPS build linked successfully, and the analytical
Parker case passed on two MPI ranks. The core qualification and exact process
receipts are in `/data/vtenishe/TEST/stage12-core-tests-20261123T172105Z`.

Initial Stage 14 tests exposed a misplaced smooth-history reference value,
an uncovered radial strip between quadrature centers, and caller-reusable
event identity. All three regressions are covered. The canonical ambient
binding passed separate PFSS/open/closed/Parker, EOS, signed sector, and
physical source-surface side checks. The solar cut-cell helper recovered a
cap whose cells were entirely removed by midpoint rejection and matched the
independent area at the predeclared 1e-9 relative gate. The full coronal registry
contains the original 222 cases plus ten background-assembly cases. Final
core checks cover checksummed complete history and piston assets, disabled
and enabled piston ownership, real expanding oblique shock patches in the
canonical Parker ambient, and failed candidate rollback. The original RH
relative bound is 1e-9. The controlled ill-conditioned high-beta case still
rejects a candidate whose RH residual cannot be resolved. It cannot publish
an epoch or serve as evidence of a physical fast shock.

The final coronal orchestrator run executed the complete 232-case inventory
with 229 PASS and three existing campaign-evidence SKIPs. All ten new cases
passed. The separate reduced-model compatibility run recorded 62 PASS and
35 existing native SKIPs. The final source-synchronized native build linked
successfully. Two-rank native admission rejection passed for all three movers.
The full maintained SEP Make target returned exit 0. It executed all 222 linked
cases with PASS, then ran all 14 bounded SWCME cases with PASS.

| Verification | Result | Retained evidence |
| --- | --- | --- |
| SWCME material upstream and original R1 | 14 PASS. All ten new cases passed strict warnings and ASan/UBSan. Nine deliberate physical regressions produced the required failures. | `/data/vtenishe/TEST/stage12-core-tests-20261123T172105Z/FINAL_QUALIFICATION.json` |
| Shared material interface and Stage 12 adapter | 20 strict PASS and 20 ASan/UBSan PASS | `material-delivery/summary.json`, `material-asan-sealed/summary.json` |
| Coronal assembly and canonical ambient binding | 10 strict PASS and 10 ASan/UBSan PASS | `stage14-final-owning/summary.json`, `stage14-sanitizer-final/owning-10/summary.json` |
| Core orchestrator | SWCME 14 PASS and 143 existing SKIP. Coronal 229 PASS and three existing SKIP. No FAIL or ERROR. | `global-swcme-final`, `global-coronal-delivery` |
| Native SEP Make target | Actual exit 0, 222 linked PASS, 14 embedded SWCME PASS | `bounded-final/process.json`, `bounded-final/native.json`, `bounded-final/native.xml` |
| Catalog, source ownership and traceability | 67 suites and 1328 cases. Existing SWCME and field-line audits passed with all 15 projection rejection controls. | `traceability-final.json`, `current-source-audits.log` |

The retained evidence root is `/data/vtenishe/TEST/sep-stage12-upstream-stage14-20261123T172952Z`.
The final native binary SHA256 is
`20a5f2782fd8130b0af1a5ff9aed0ca45a9f79cc859405fba52fb3f5fc8f6a21`.
Rebuilding the final piston-input guard and the matching source tree yielded
the byte-identical executable used by the full native regression. The exact
Make arguments and elapsed time are in the process receipt. ASan and UBSan
ran with leak detection disabled because this execution environment uses
ptrace and the default leak probe aborted.

The first bounded attempt retained an application-local forward declaration
and stale source tests copied from the earlier configured build. Moving the
declaration to the SWCME core and synchronizing the full tracked test/source
closure resolved the ownership and source-audit failures. The current
bounded run repeats the maintained checks. The first sanitizer attempts that
overlapped documentation edits were rejected by the source hash check. The
reported sanitizer results use the final frozen sources. A separate raw
PARKER02 probe lacked its material fixture and returned the same error on the
earlier and new executable. The maintained analytical PARKER01 probe passed
on two ranks with its declared fixture.

## Stages 10–11 complete for the approved scope, 2026-11-23

Stage 10's baseline and capability inventory are complete. Stage 11's common
contract, supported reduced-model physical campaigns and maintained regression
all pass. The reviewed fixtures and test entry points are integrated. Four
matching application builds pass their validation campaigns. The prototype
receipts preserve the before/after fixture evidence.
Earlier dated sections below are historical records, not current build receipts.

The repository is `/data/vtenishe/TEST/AMPS`; the recorded starting HEAD is
`80466e68e37fd776bae10436e66c04cbbc88c4a7`. The durable evidence directory is
`/data/vtenishe/TEST/sep-stages10-11-20261123T030052Z`, also recorded in
`/data/vtenishe/TEST/sep-stages10-11-current.txt`. Its `RESUME_STATE.md` contains
the detailed implementation handoff. Its baseline manifests, patch and
`baseline-repo/` preserve the original dirty checkout. No files have been staged
or committed in this batch. The pre-existing root `Makefile` edits and unrelated
untracked application assets must be preserved.

### Physics and implementation boundaries

- AMPS physical timesteps remain constant. Do not shorten a host step, change
  sampling assumptions, weaken numerical/statistical gates or regenerate
  reference solutions to obtain a PASS.
- The new common material-background interface declares SI/HCI conventions,
  composition, density basis, complete velocity/magnetic gradients, immutable
  owner/interval identity, geometry generation and coverage. Missing mandatory
  values are invalid, rather than invented physical zeros.
- The twelve explicit backend profiles distinguish snapshot exchange, material
  transport and physical downstream support. Only the already supported reduced
  SSE/Parker exterior is qualified for material transport. SWCME1D, SWCME3D,
  coronal-CME and full regional material transport remain explicitly rejected
  where their required physics is not implemented.
- A real defect in the reduced material adapter was fixed: its neutral `Move`
  interface accepted points outside the declared Parker region even though
  sampling rejected those points. It now checks descriptor/region validity at
  entry, every existing RK evaluation and the endpoint. The legacy trajectory,
  equations and integration controls are unchanged. The focused test first
  reported 19 PASS / 1 FAIL, then 20 PASS after the fix.
- Registered native geometry generation zero is valid when explicitly present;
  an absent generation or unregistered owner remains invalid. Public snapshot
  ABI/schema 5 and the canonical seven-member core archive are preserved. The
  new schema-1 codec carries material metadata only and is not a restart image.

### Verified application results before prototype integration

Matching builds are retained as `native-region`, `native-checkpoint-region`,
`native-thread2-region` and `native-thread4-region` under the evidence directory.
Each build completed successfully after recompiling the affected PIC consumers;
copied dependency files were not rewritten to impersonate a rebuild. Maintained
source audits pass. The ordinary executable SHA256 is
`1230472e7037b676e31942f4f17c1b03be1815c09677913513940a1b020b888b`.

| Gate | Actual result | Evidence directory |
| --- | --- | --- |
| Material contract, strict warnings | 20 PASS | `contract-region-strict` |
| Material contract, ASan/UBSan | 20 PASS; leak detection disabled, no leak-coverage claim | `contract-region-sanitizers` |
| Shared material registry | 20 PASS | `shared-material-region` |
| PCLI01–12, two ranks and repeated execution | 12 PASS | `pcli-region` |
| Native integration, 1/2/4/8 ranks | 48 successful jobs covering all 15 required IDs | `native-region-tests` |
| Unsupported backend rejection | 3 PASS with required native rejection exits | `backend-native-region` |
| Material transport campaign | 58 PASS; maximum budget residual 1.071e-15 | `material-region` |
| Checkpoint/threaded campaign | 44 PASS; maximum budget residual 1.467e-15; byte-exact images up to 21,425,968 bytes | `checkpoint-region` |
| Controlled CV01–12 / IV01–06 | 18 PASS | `controlled-region/results` |
| Full 400-step applications, three movers | All three MPI runs exit 0; three original verification gates and eleven corruption controls PASS | `application400-region` |
| Orchestrator / Python harness | 51 / 27 PASS; catalog 67 suites / 1308 cases and traceability PASS | `region-harness-compatibility.json` |

Each full application retains 400 physical steps of 300 s, 357 segments and
sample exposures from 5 through 400. The largest relative budget residual is
2.203e-15 against the unchanged 1e-10 limit. These are actual completed runs,
not extrapolation from a shorter campaign.

MPI is available. The earlier blanket statement that MPI initialization was
blocked was incorrect: a direct `mpiexec -n ... --bind-to none
--output-filename ...` invocation succeeds. Some outer wrappers selected a
restricted execution path and failed before MPI initialization. The retained
SWMF compatibility evidence also contains nine successful MPI phases and nine
successful numerical comparisons. Those historical results establish launcher
availability, not qualification of later source changes. Sequential Python
campaigns use the documented evidence-only `independent_runner.py` to remove
inherited finalized-rank environment state before independent native children.

### Original complete-regression failure and fixture repair

The initial `make -j1 -k test test-controlled-analytical` run exited 2. Source
prerequisites and the controlled comparisons pass, but bare native
`--all-tests` chooses one maximum initialization level for 222 heterogeneous
Routine cases. Its legacy field-line bootstrap reaches the required material
owner guard before callbacks, so the authoritative report records 222 ERROR.
These are not 222 independently evaluated numerical failures. The failed run
and diagnosis are retained in `bounded-region`; no ERROR has been waived.

The roster is 193 cases without that native initialization requirement plus
29 needing specific fixtures: sixteen storage/accessor cases, five existing
reduced CP05 cases, two zero-centerline FTE/Parker cases and six birth-lifecycle
cases. A diagnostic selection of the 193 gives 175 PASS and eighteen explicit
manifest-dependent SKIP; those eighteen separately pass the full controlled
CV/IV campaign. All 222 cases remain required in the proposed isolated driver.

Prototype source is `routine-prototype/AMPS`; its separate build is
`native-routine-prototype`. The seven reviewed C++ files now match the
application checkout and the four new configured builds. The maintained
Python, Make and orchestrator entry points now use the same reviewed fixtures.
The isolated driver preserves exact test IDs, thresholds,
seed derivation and repeat semantics, routes cases to their actual required
fixtures and rejects missing roles, missing cases and native failures.
Storage initialization does not grant production transport capability.

- The first complete prototype run reports **216 PASS / 6 ERROR**, actual exit
  2. The remaining errors are BIRTHN01–06. Results are retained in
  `routine-prototype/full222`.
- Focused prototype FTE/Parker, their seed/repeat checks, twelve storage cases
  and four PNW snapshot-only cases pass. PNW uses its actual SWCME snapshot
  diagnostic and does not admit unsupported moving transport.
- The matching v5 birth campaign completes successfully: **all six stable IDs
  PASS, all 24 native role jobs PASS**, two fresh repetitions on two MPI ranks,
  campaign seed 123. Every assertion count and Cartesian residual is zero.
  BIRTHN01–03 retain the exact 0.0025-s remaining-duration comparison and full
  parent-step negative; BIRTHN04 retains a positive late duration of
  1.1102230246251566e-18 s; BIRTHN05 retains actual interrank transfer and six
  representatives after the exact split. Evidence is
  `routine-prototype/birth-focused-v5-all6`.
- BIRTHN06's v4 ledger and foreign-list failures were caused by the fixture
  omitting the production post-mover `ExchangeFieldLineParticles` phase. With
  that actual migration phase restored in v5, both repetitions pass: sixteen
  genuine source survivors, exactly zero number/birth-energy/energy residuals,
  a successful completed checkpoint, and the intended active-window/pending
  negatives without mutation of particles or the existing destination. No
  production mover, source law, boundary policy or timestep changed.
- Independent review subsequently found a real reporting defect in v5:
  external replicate 1 used the correct distinct random streams, but its
  `test_case_seed` metadata was stamped after restoring replicate-0 context.
  The v6 correction extends the context lifetime through the complete registry
  call; random-stream laws and values are unchanged. Role-specific clock and
  duration wording and aggregate multi-seed reporting are also corrected.
  The v6 campaign now closes these focused provenance gates as well as the
  numerical assertions; the earlier defective reports remain archived.
- The matching v6 all-six/two-repetition/two-rank campaign completes with
  **all six IDs PASS, all 24 native roles PASS, actual launcher exit 0**, with
  strict seed-metadata checks. Its immutable executable SHA256 is
  `df20103f7ec60304422cae3b7d1d7eb706cbe76057d31c3cc0bd3ac42960d35f`;
  evidence is `routine-prototype/birth-focused-v6-all6`. Actual two-rank
  missing-reference and raw-repeat negatives also reject for the intended
  reasons. The same v6 executable then passes **all 222 routine cases** in
  the prototype with launcher exit 0 and no FAIL, SKIP or ERROR. Evidence is
  `routine-prototype/full222-v6`. The final maintained checkout qualification
  is recorded below.
- All five v4 admission controls now pass with actual launcher exit 0. The
  initial multiselection check expected the wrong exit status; the existing CLI
  correctly rejects earlier with exit 1 and its exact singleton-input reason.
  The corrected harness checks that precise exit/reason. Its first failed
  receipt is retained separately; native admission behavior was not weakened.
- Thirteen isolated-driver protocol controls pass, including manifest startup
  override and midrun-mutation rejection. Actual reduced input/event provenance
  and runner identity are explicitly bound; no arbitrary manifest argument may
  override the selected checkpoint fixture.

Birth-lifecycle comparisons use the documented operator splitting: the particle
operator consumes frozen departure geometry/fields, then the material map is
published after readers close. All comparison roles keep the physical parent
step at 0.01 s. The existing controlled diagnostic mover API supplies a 0.0025 s
suboperator reference; this does not change a physical AMPS host timestep or
introduce dense-time production evolution. Repeat role groups require fresh
processes with the existing derived seeds. Checkpoint negatives must fail for
their intended reason after actual geometry/source lifecycle updates.

Independent historical nonpasses are also preserved: the unsupported full
regional MASS-TIME campaign reaches its original 180 s ceiling, and a separate
packaging helper has an existing CLI mismatch. Neither is silently reclassified
as supported physics or a passing gate.

### Final checkout qualification and reboot recovery

The new `native-integrated`, `native-checkpoint-integrated`,
`native-thread2-integrated` and `native-thread4-integrated` builds all exit 0.
Each build recompiles 70 PIC consumers to regenerate dependencies for its own
configured tree. `v6-reviewed-cpp-applied.patch` and the four build receipts
bind the exact matching C++ bytes. Fresh MPI evidence is stored in
`native-integrated-tests`, `material-integrated`, `checkpoint-integrated`,
`backend-native-integrated`, `pcli-integrated` and `application400-integrated`.
The 48 native jobs, 12 CLI cases and three unsupported-backend negatives have
passed. The material contract passes all 20 cases in the strict warning,
ASan/UBSan and shared registry runs. The complete controlled analytical Make
target passes in `controlled-integrated` and the selected real orchestrator
steps pass in `orchestrator-integrated`. The standard bounded Make suite passes
in `bounded-integrated`. `maintained-wiring-applied-files.json` binds the
reviewed Python, Make and orchestrator patch. The current checkpoint collector
passes all 44 roles. An external logging command overflowed the process
argument limit after two threaded roles. Both raw logs are archived and the
matching roles were rerun successfully with bounded process receipts.

The matching requalification results are in the following directories.

| Evidence | Result |
| --- | --- |
| `material-integrated` | 58 PASS with maximum relative budget residual 1.071e-15 |
| `checkpoint-integrated` | 44 PASS with maximum relative budget residual 1.467e-15 and up to 21,425,968 exactly matched checkpoint bytes |
| `application400-integrated` | All three application exits and independent verifications pass. Each run completes 400 steps at 300 s, 357 sampling segments and exposures from 5 through 400. The largest budget residual is 2.203e-15 against 1e-10. All eleven corruption controls pass. |
| `orchestrator-integrated` | Ten PASS at two ranks, comprising all four PNW and six BIRTH native cases |
| `birth-repeat-integrated` | Six BIRTH cases and all 24 two-rank role jobs PASS with exact campaign-123 seed metadata across two fresh repetitions |
| `selected-cp05-integrated` | `make test-case CASE=CP05C01` passes and exports the actual JSON/JUnit records |
| `controlled-integrated` | `make test-controlled-analytical` exits 0. Every CV01–12 and IV01–06 comparison passes. |
| `bounded-integrated` | `make test` exits 0. All source prerequisites, exactly 222 SEP Routine cases, and all four embedded SWCME Routine cases pass. The SEP and SWCME reports each contain zero FAIL, SKIP or ERROR. |

No additional physical-model decision is needed to complete this bounded
stage-10/11 batch. Stage 10's exact audit, twelve-profile matrix and open
decision records are established. The final evidence audit in
`stage10-11-qualification.json` verifies all 222 maintained process receipts,
the full source and fixture identities, and the four configured executable
hashes. The full-regional backend remains unsupported and its original 180 s
timeout remains recorded. Physics development for the additional SWCME,
coronal and full regional profiles belongs to the later stages.

The completed agent handoffs are `field_line_scope_gate` (isolated driver and native
initialization), `material_contract_review` (birth fixtures), and
`threaded_checkpoint_gate` (independent review). Agent/tool process handles are
ephemeral. Their handoff files are
`routine-prototype/FIELD_AGENT_RESUME.md` and
`routine-prototype/BIRTH_AGENT_RESUME.md` in the evidence directory.

The installed Codex CLI supports resuming this exact session after reboot:

```sh
cd /data/vtenishe/TEST
codex resume 01a1fa84-5ac6-74c3-8297-995ef2a4c6db
```

After resuming, request: "Read `AMPS/srcSEP/STATUS.md` and
`sep-stages10-11-20261123T030052Z/RESUME_STATE.md` and inspect the current
checkout and completed stage 10–11 receipts. Preserve the fixed physical
timestep and all validation gates when continuing SEP work."

The session transcript is retained under
`/home/vtenishe/.codex/sessions/2026/11/20/` with this session UUID. Preserve that
home-directory session store as well as `/data/vtenishe/TEST`. A reboot stops
running MPI tests and agents; it does not turn partial output into a completed
PASS. Inspect existing receipts, verify source/binary provenance, and restart
unfinished runs in fresh evidence directories. If the CLI session record is
unavailable, start `codex` in `/data/vtenishe/TEST` and give the same handoff
instruction. Do not discard the dirty checkout or treat old process IDs as live.

## Rebase conflict resolution, 2026-11-17

The user-initiated rebase applies `7984c03ac5d4aaa6f9d6e7c7b4b189716effc798`
onto upstream `6bd40eacde14b810a199151c1c884513c7fe7cc9`. The two conflicts
in `src/meshAMR/cut_cell.hpp` and `src/meshAMR/meshAMRgeneric.h` contain only
the four whitespace cleanups recorded below. Both are resolved using the
exact upstream files. The automatically merged `MakefileTest/Table` also
duplicated four core campaign registrations already present upstream;
removing those exact duplicate stanzas restores the exact upstream table.
The SEP hooks in `src/pic/pic.h` and `srcSEP/BoundaryInjection.cpp` remain
byte-identical to the replayed commit. No numerical algorithm changes.

The existing table-parser suite passes **4/4**; the complete parsed registry
matches upstream (**162** entries, including existing host-selector variants).
Each affected core campaign occurs twice in the saved failed-merge table and
once after resolution. The broader native/core campaigns are not rerun for
these whitespace/registration changes. Checks, saved conflict stages and
the original index are in `/tmp/sep-rebase-conflicts.j5ge1br4`; see
[CORE_UPSTREAM_INTEGRATION.md](CORE_UPSTREAM_INTEGRATION.md) for the command.
The conflicts are marked resolved, but the rebase is left paused in accordance
with the user's instruction not to create commits. The following integration
and executable records describe the earlier `74974c2` integration, not a
fresh qualification of the newer rebase target.

## Upstream core integration completed

The user requested integration of the major upstream core update, with maximum
reuse of upstream implementations and preservation of previously passing core
tests. Local base: `c65494f5a6eff5988136db181abd3750b7cec8a7`.
The verified `origin/master` target is
`74974c2b2fc1b4cf09a73c320064a22370c0ca94`.

The dirty local source and proposal files are backed up before any integration
in `/tmp/srcsep-upstream-merge.1q31gpks` (264 files plus a binary Git diff).
The upstream source and clean three-way candidate were verified in isolated
checkouts. Complete unchanged-upstream and merged-source `CORE-*` campaigns each
select 51 suites / 1,085 cases with identical verified tools, fixtures and
independently built dependency archives. Initial sandbox MPI and
missing-private-SHARE attempts are archived as setup failures, not physical
regressions. The corrected campaigns ran outside the socket-restricted
sandbox. New upstream native core campaigns and the SEP build used separate
private checkouts. No build11 PASS is evidence for this new upstream core;
the complete core PASS-preservation gate now passes. Existing user changes
and unrelated application files are preserved.
The clean 65-path upstream candidate, one authenticated required external
test fixture, one retained-local compatibility amendment, source hashes and
acceptance gates are described in
[CORE_UPSTREAM_INTEGRATION.md](CORE_UPSTREAM_INTEGRATION.md).

Focused progress: the principal native PIC owner passes **95/95** cases in
both pristine and merged trees. The new upstream tetrahedral cut-cell matrix
passes **22/22** launches, **184** assertions and **11/11** exact one/four-rank
cell-file comparisons on the merged tree. All seven shell/specular surface-mover
launches pass, including both exact one/four-rank digest comparisons. All four
new upstream campaigns now pass **32/32** expected launch outcomes, **578/578**
assertions and **13/13** exact rank comparisons; the gas profile passes 310
checks and both missing-input modes reject with the required diagnostics.
All 327 copied source files remain unchanged in its final audit. The complete
core comparison subsequently passes after the compatibility amendment below.

The isolated merged-core SEP build succeeds, archived SHA256
`cc6ec6b7be9a23636eab85b94aded66a63e1c41465d4994f38e37f6568900f72`.
Its actual two-rank combined regression passes **121/121**, HOST01–03 **3/3**,
RADIUSN01–02 **2/2**, fresh-process CP05F01–03 **3/3**, and PNW01–04 **4/4**
at each of N=128 and N=129. Native CLI passes **9/9**; production CLI software
retains **18 PASS / 2 explicit SKIP**. Separate checkpoint process pairs have
identical final artifacts. Reports are under
`/tmp/srcsep-upstream-merge.1q31gpks/sep-evidence-merged`; its independent
audit passes. The two production/SWMF qualification SKIPs remain open.

The completed full-core comparison exposed a local shared-header compatibility
defect: the optional test-ID field broke an existing two-field aggregate under
`-Werror`. The original merged receipt retains its 97 reduced-owner ERRORs;
no physical assertion was relabelled. An explicit empty member initializer
preserves the upstream application unchanged. Strict fail-before/pass-after
compile and layout/default checks pass. Fresh retries of only the two affected
owners execute 240 cases (153 PASS / 87 SKIP / no FAIL or ERROR).

The reconciled complete core result is **843 PASS / 1 FAIL / 240 SKIP / 1 ERROR**,
exactly matching pristine upstream and the historical 843-PASS identity set.
There are no missing cases, new nonpasses or changed statuses. The two existing
nonpasses are strict incomplete coverage and the 180-second CME refinement
ceiling. Raw full receipts, scoped retry provenance and source/dependency
identity checks remain archived. The post-amendment SEP executable SHA256 is
`13dd767b86ed65174af628cba6f70a875a1b486c914bbc4802452ce89e6bf1e1`;
its generated configuration is unchanged and all nine native/CLI/restart/
spectral campaigns pass. An independent post-amendment audit also passes.
The accepted 67-path source payload is now published and the rebuilt workspace
executable passes the same gates, as recorded at the end of this file.

## Previous verified milestone (before upstream integration)

**Build11 passes the bounded native regression. The complete production and
physical-validation plan remains open.** These results use the same archived
executable. Component, host, radius and restart campaigns use two MPI ranks;
the CLI campaigns also include fresh-process and source checks:

| Gate | Latest result | Scope |
|---|---|---|
| Combined native regression | **121 PASS / 0 FAIL / 0 SKIP / 0 ERROR** | Controlled components and existing native seams; includes the five source-epoch cases and material numerical repairs. |
| Actual host/PIC driver | **HOST01–03: 3/3 PASS** | Two real steps and eight first-birth moves per mover case; fixed analytic diagnostic. |
| Translated provider radius | **RADIUSN01–02: 2/2 PASS** | Actual initialization and scalar/library coefficient paths; unknown origins remain unavailable. |
| Separate-process restart | **CP05F01–03: 3/3 PASS** | Separate write/resume MPI jobs, byte-identical final artifacts; fixed analytic geometry, same rank count/ownership. |
| Native test CLI | **9/9 PASS** | Applied seed/repeat/reporting and failure propagation. |
| Production CLI software | **18 PASS / 2 explicit SKIP** | Post-allocation production evolution and actual coupled SWMF execution remain unqualified. |

The largest remaining development items are actual native material-grid/ALE
evolution, positive-energy joint particle/wave publication with consistent
stochastic hazards, the Maxwell–Jüttner gyrotropic source bridge, live
four-spectrum restart and broader MPI/thread qualification. Equivalent-input
AMPS/M-FLAMPA experiments, full SWMF replay and instrument-aware observational
validation are not complete. Reflection/cascade remain disabled for the first
qualification case; their original requirements remain open.

The finite-amplitude study has **31 numerically resolved / 5 phase-unresolved
cases out of 36**. No physical accuracy tolerance or production closure limit
has been selected. Component PASS results do not resolve those physical inputs.

## Scientific decisions adopted

- Exactly three field-line movers: `parker`, `fte-dmumu`, and `fte-mfp`.
- Local wave storage is magnetic fluctuation energy. The user selected standard
  BATSRUS imports: raw total wave energy density converts by a factor of one
  half to the local magnetic-energy view. Spectral entries are joules
  integrated over each bin; divide by volume for J/m³. The conversion is
  `deltaB² = 2*mu0*W_B`.
- Local `+`/`-` propagation is along/opposite **B**, separately from geometric
  field-line orientation. The inspected standard-energy SWMF producer maps Wave1 along B and Wave2
  opposite B; arbitrary exports/representatives still need their exact mode and
  energy/composition contract. See the source audit linked below.
- Adopt the small-amplitude nondispersive MHD Alfvén closure, with equal
  magnetic/kinetic wave energies and an explicit prescribed-background energy
  account until plasma feedback exists. No numerical domain limits have been
  adopted.
- Both focused movers must use the same actual spectrum-derived joint
  pitch/momentum law. Preserve finite-band resonance gaps; report infinite or
  undefined derived mean free paths where appropriate.
- Reuse the existing initializer. Begin each branch with equal L/R energies,
  then evolve all four independently. Existing logarithmic k values are bin
  centers, with half-logarithmic outer edges.
- Continuously supply the initial **magnetic energy density per bin** on each
  incoming spatial characteristic, determined relative to grid motion. Let
  outgoing waves leave. At spectral edges supply zero incoming power and keep
  outgoing energy/action in a distinct external spectral account.
- Disable reflection and cascade for the first qualification case. Their full
  original-plan requirements remain open pending additional physical closures.
- Use a **material magnetic tube** for that case: vertices move with plasma U
  and particle coordinates use motion relative to the moving grid. Do not add
  the full U convection again in a material coordinate.
- Replace unqualified thermal injection with a drifting relativistic
  **Maxwell–Jüttner flux** distribution and verify its low-temperature Maxwellian
  limit. The full-particle material-plane kernel now exists; its native
  gyrotropic source reduction remains pending.
- Support corrected-model checkpoints only.
- The dispersion study uses explicitly declared parallel proton/electron
  Maxwellian cases. Report error curves; proposed 1%, 5%, and 10% contour labels
  do not select accepted physical tolerances.

## Implemented and verified portions

Rows retain evidence from their recorded builds. A new source revision does
not inherit a PASS merely because the preceding build passed a related test.

| Area | Implemented behavior | Evidence and limits |
|---|---|---|
| Persistent randomness | Versioned operator streams, cached normal state, source event identities, MFP residual hazard and split lineage | `RNG01`–`RNG10` source gates pass. Native source/population seams have separate tests; full decomposition/thread qualification remains open. |
| Provider derivatives | Explicit radial/tensor gradients; no temporal-density-only divergence or projected-velocity substitute for strain | Native `CP03G01`–`CP03G04`, `FTE01`, `PARKER01` pass. Missing SWMF gradients and FULL_ICME thermodynamics reject explicitly. |
| Wave normalization | Common `2*mu0*W_B` amplitude conversion, relativistic gyrofrequency, polarity handling; corrected one-helicity and balanced Jokipii prefactors | Independent `RWN01`/`RWN02` Lorentz/Kubo regressions pass. Two-rank normalization/value campaign: 14 PASS. No tolerance relaxation. |
| Population control | Actual PIC exact splitting preserves parent history; children get new identities/streams. Positive-weight resampling retains existing phase points and declared relativistic/pitch moments | Two-rank native `CRS01`–`CRS10` + `POPN01`–`POPN03`: 13 PASS. Distribution changes and unconstrained moments remain explicit. |
| Joint collision kernel | Conservative four-polarization, finite-band pitch/momentum operator; fixed-wave event history and independent PDE refinement | `RC01`–`RC06` source gates pass. Full PIC particle/wave transactions are not implemented or qualified. |
| Four-spectrum authority | Immutable value owner, exact bin integrals, dynamic native diagnostic shape, collective staging/readback/restore readiness | `PWS01`–`PWS08` source gates pass. `PNW01`–`PNW04`: 4 PASS at N=128 on one rank, 4 PASS at N=128 on two ranks, 4 PASS at N=129 on two ranks. Global production startup/evolution remains disabled. |
| Wave initialization | Collective validation/staging replaces malformed-B/geometry/profile fallbacks and partial writes | Native `PNW04` passes. Further whole-state geometry-cache transaction qualification remains open. |
| Analytic startup | Existing analytic current/previous density law now installs in shared startup before particle reads or test return | Fixed the actual CP05 zero-density failure. Current checkpoint continuation now passes on one and two ranks. |
| Physics checkpoint, current verified scope | Actual attached PIC particles, RNG/source history, fixed analytic geometry, native wave arrays/private caches and physics clock; staged atomic load | `CP05C01`–`CP05C04`: earlier one/two-rank gates pass; build 6 additionally passes all four transport-schema-4/hardened cases on two ranks, using same-process save/load. Inner checkpoint schema 4 explicitly encodes XYZ and validates decoded position/hazard consistency; build8 CP05C01–05 and fresh CP05F01–03 pass on two ranks. Build 7 (inner schema 3) adds actual flux-table persistence and passes CP05C01–05 plus CP05F01–03 in separate two-rank process pairs, with byte-identical uninterrupted/resumed final artifacts. Same rank count/ownership only; output continuation and live four-spectrum codec are unsupported. |
| Source boundary safeguards | Unsupported Cartesian/thermal paths reject before particle/source RNG mutation; generic Maxwellian geometry/initial owner fixes | Two-rank `SOURCE01`/`SOURCE02`: 2 PASS. This does **not** qualify the Maxwellian sampler, physical frame, birth timing or crossing behavior. |
| Test CLI | Applied seed/repeat selection, reports, per-rank evidence, failure propagation and bounded unsupported-repeat handling | Build 6: native MPI campaign 9/9 PASS; production software suite 18 PASS and 2 explicit qualification SKIPs. PCLI11 production advancement/PCLI12 actual SWMF execution remain open. |
| Frozen-wave transaction staging | Stable accepted event/phase/weight/provenance validation, deterministic all-bin budgets and representable wave transfers; no public commit | FCT01–08 strict/ASan/UBSan and two-rank build-6 callbacks: 8 PASS. Encoded particle-energy follow-up passes FCT01–09 source tests and two-rank native build-6c callbacks (9/9). Native coupled publication/rates remain unqualified. |
| Birth-time value layer | Explicit immutable driver window, positive remaining birth duration, late-time rounding safeguards, actual clock-policy preflight and atomic staged consumption | `BIRTH01`–`BIRTH08`: 8 PASS with ASan/UBSan; related RNG 10/10 and snapshot/clock source gates pass. Build 6b BIRTHN01–06: 6 PASS on one rank and 6 PASS on two ranks, including actual mover, clone/split/migration and unfinished-checkpoint rejection. |
| Warm-plasma reference | Independent complex parallel Maxwellian dispersion, precision/continuation/isolation checks and reproducible archived plotting | Mathematical/protocol tests: 19 PASS plus independent review. Controlled support scan: 231 resolved / 29 right-branch unresolved roots. No production cutoffs or amplitude qualification. |
| Controlled four-spectrum WKB | Staged signed reservoirs, work/ray transport, side/cap pressure-work ledger and exact binary64 external accounts; no native publication | WKB01–08 strict/ASan/UBSan: 8 PASS on the pinned source revision; independent exact-rational accumulator review: 2,007 vectors PASS. Build 7 also passes all eight callbacks on two MPI ranks. Actual native moving-tube publication remains pending. |
| Maxwell–Jüttner flux | Exact full-particle material-plane rate and gamma-mixture sampler; declared rest/simulation momenta, staged keyed randomness | Build-7 MJ01–06 passed before independent review exposed cold-boost encoded-crossing loss and high-Gamma conditioning defects. The expanded source gates and build8 MJ01–06 now pass with explicit resolution rejection and unchanged RNG/output; native injection remains disabled. |
| Numerical M-FLAMPA baseline | Pinned upstream source/dependency build and exact upstream `test_mflampa` experiment against already-local references | Two-rank solver and upstream DiffNum exit 0; 33 output blocks / 75 checks pass. This establishes the upstream regression baseline, not AMPS equivalence or observations. |

## Implementation history and outstanding integration

1. **Native birth bookkeeping passed:** schema 4 and BIRTHN01–06 pass on one
   and two ranks after the explicit later-epoch fixture correction. Keep the
   original build-6 duplicate-publication failure archived. An actual full
   PIC::TimeStep/host-driver smoke test subsequently passes build11 HOST01–03
   on two ranks. Physically qualified thermal/moving source paths remain open.
2. **Checkpoint persistence:** schema-4 decoded B/U/n/Alfvén and CRC-valid
   malformed/pending-state rejection pass all four native CP05 cases on two
   ranks. Geometric flux-registry persistence, opaque callback rejection,
   CP05C05, FLCP01–04 and paired CP05F01–03 pass build 7 on two ranks.
   Each fresh case uses separate write/resume MPI jobs and requires identical
   final checkpoint bytes. CRC-valid derived-position and strict hazard
   validation subsequently pass build8, including all three fresh-process pairs.
3. **Four-state authority guards:** integrated legacy wave writers now reject
   a live four-state owner in source. Current build8 PNW01–04 pass on two
   ranks at both N=128 and N=129; complete global-owner production integration
   remains pending.
   Full four-state checkpoint serialization remains explicitly unsupported.
4. **Controlled wave transport:** the adopted signed-characteristic reservoirs,
   spatial/spectral accounts and ideal-MHD WKB work/ray equations now pass their
   eight controlled source cases and all eight native-selectable callbacks
   in build 7 on two ranks.
   Actual native material geometry/provider motion and coupled publication are
   not implemented or inferred from those fixtures.
5. **Runner isolation corrected:** native `run_tests.py` used the source directory as its
   working directory and wrote `PT` diagnostics there. Native execution
   now uses artifact-owned working directories with preserved input-path
   semantics and verified isolation. This work's generated PT logs are archived in
   `/tmp/srcsep-development.AilCx1/native-runner-source-output`; the initially
   clean tracked log has been restored. Source-only Make execution keeps its
   source cwd.
   Focused tooling regressions now pass: Python runner 27/27, orchestrator
   28/28, paired-process harness 3/3, actual two-rank CLIEX01 seed/repeat execution and serial CV02 with
   its independent reference. Broader native suites must use the fixed runner.
6. **Frozen-wave transactions:** a prepare-only neutral component validates
   stable-keyed weighted proposals, preflights every
   bin budget and stages particle/wave changes with signed ledgers. It does not
   publish native state or establish an exact coupled competing-hazard process.
   Build-6 FCT01–08 pass on two ranks. The expanded FCT01–09 source revision
   now checks actual encoded particle energy with a stable difference formula,
   separate particle/wave roundoff and sub-resolution rejection; strict
   ASan/UBSan passes 9/9 with unchanged input hashes, and two-rank build-6c
   native callbacks also pass 9/9. The combined native regression passes 74/74.
7. **Independent position/thermal review:** independent review found that old mover commits left
   PB::X stale relative to the field-line coordinate. New BIRTHN01–04 controls
   independently check those views and invalid-geometry atomicity. MJ05 exposed an encoded cold-boost crossing that became outgoing after
   binary64 storage. The fixes now reject unresolved encoded crossings and
   ill-conditioned high-Gamma arithmetic; all expanded BIRTHN and MJ cases
   pass build8 on two ranks.
   Build 7a captures test-only changes before the production corrections;
   native fail-before is 7 PASS / 5 FAIL, exit 1, no ERROR. The failures are
   BIRTHN01–04 and MJ05. The analytic ARBN01 check also exposes missing
   previous B/U. ARBN02 uncovered a separate fatal geometry-builder/full-run
   validator mismatch; a nonfatal fixture now preserves that failure and a
   geometry-only API validation repair now passes build8 ARB01–06 and
   ARBN01–04, including rank-specific bad-input collective rollback and
   mover-independent background identity. Build 7c makes
   that fixture nonfatal and returns the expected 0 PASS / 2 FAIL, exit 1,
   no ERROR; build 7b retains the original fatal failure diagnosis.
8. **Material kinematics:** the neutral Cauchy/Nanson and global-affine
   cap/side GCL kernel passes six source cases, including exact curved
   toroidal reference integrals and explicitly measured quadrature error.
   Final hash-pinned strict/ASan-memory/UBSan checks pass 6/6; native factory
   execution also passes build8 MTG01–06. Subsequent independent extreme-shear
   review exposes loss of a small nonzero flux after binary64 narrowing; the expanded
   MTG06 regression fails before the fix (5 PASS / 1 FAIL) and all six now
   pass with certified signed-dot comparisons of reference, wide and stored
   B/area. An independent exact-rational review passes 2,000 cases, including
   334 accepted states and 66 accepted exact-zero fluxes; accepted nonzero
   flux errors stay below the declared numerical bound. Curved fixture mass
   now matches its explicit nonuniform-density integral (1.06e-20 kg).
   Material-only build9 links successfully, SHA256
   `ae18e39ca7169943734ea2be4f53f4c99b4f608e2bd3aac3b8bd4f3dd0b486b6`.
   It passes all **31/31** MTG/WKB/transaction/value callbacks on two ranks,
   including the expanded encoded-flux and rounding-mode controls. LeakSanitizer detection
   is disabled and is not claimed. No
   native material remap or particle/wave joint commit is inferred.
9. **Standard BATSRUS view:** user-selected total-to-magnetic conversion and
   along/opposite-B mapping are implemented as an immutable read-only
   adapter with explicit producer provenance. Source SBW01–05 and native
   SBWN01–02 now pass build8. No imported L/R or composition is
   invented; full provider execution remains unqualified.
10. Add detailed physics/units/ownership/numerical/failure comments to all new
   code and reconcile its relevant Markdown documentation. This is a standing
   user requirement, not a final cleanup substitute.

## Current verification record

**Build 8 linked successfully** from a pinned 162-file native-source/build-copy
manifest plus 23 canonical shared-source hashes. It includes corrected X and
checkpoint inner schema 4, expanded MJ resolution controls, ARB01–06 /
ARBN01–04 collective startup/identity checks, MTG01–06 and SBW01–05 /
SBWN01–02 read-only BATSRUS views. Source gates pass; native results are from the actual build8 execution,
not inherited from build7. Its SHA256 is `049a1959918dd154e5a81d75ed7ebb18a57764f29c09c3917c81b4068f9c58eb`.
The archived `amps.build8` binary passes the combined two-rank native
regression **116/116**, all three separate-process checkpoint pairs **3/3**,
native CLI **9/9**, and production CLI software **18 PASS / 2 explicit SKIP**.
All child jobs exit 0. The two SKIPs are post-allocation production evolution
and actual coupled SWMF execution; neither is counted as qualified. The
independent strict paired-process reader also passes all three cases after
rehashing the archived binary and reading all six native child reports and
actual continuation artifacts.

The combined two-rank build-6c regression passes **74/74**, exit 0, with no
FAIL, ERROR or SKIP. It includes birth windows, same-process checkpoints,
source safeguards, native population control, RNG, provider derivatives,
mover fixtures, wave normalization, four-state values, joint collision and
prepare-only transactions. Report:
`/tmp/srcsep-development.AilCx1/combined-mpi-build6c/srcsep-tests.json`.
This remains controlled/component evidence, not full material-tube production
evolution or observational qualification.

The earlier **test-only build 7c** has SHA256
`3e821fd1eb7275d4dab3d1d8a687928ad95ee59782458539790b37e3b1a1831f`.
It preserves the two clean analytic fail-before results. It does not contain
the subsequent fixes. Build **7** succeeds with `make -j4 amps`, with SHA256
`ecde38fd9dbb07228c36d74620d65ba7fe467ed79ecd98c2c9b4ca5c66cbff43`.
Its new-kernel/flux/checkpoint campaign passes **23/23** on two ranks; paired
fresh-process checkpoint continuation passes **3/3** across all three movers.
The first invocation used the Python alias as a native group and rejected the
unknown selector before initialization; corrected native `flux-checkpoint`
selection passes. Reports are in `new-kernels-mpi-build7-corrected-selector`
and `fresh-checkpoint-build7` under the evidence root. Build **7a** adds only
the independent position and thermal negative controls, SHA256
`0158318da83901e45157fb21a960dd2cfa8027c4b31ff5aedd1781ebf703c998`;
its fail-before gate returns the expected 7 PASS / 5 FAIL, exit 1, no ERROR. Build 6c retains SHA256
`4759f4aa39139591a0d70c3131d9d89e7a1c8ba4e58be17666630cd617d10ca4`.
Authoritative sources have subsequent changes; this binary is not evidence for
those unbuilt changes. Configured build copies preserve generated AMPS macros.
Build 6 uses a pinned 120-file source/build-copy manifest, including schema 4,
birth-native fixtures and initial frozen-wave staging. Its first birth runtime
gate returned 1 PASS / 5 ERROR on each rank layout.
Build 6b changes only the birth fixture's static epoch setup; corrected reruns
pass 6/6 on one rank and 6/6 on two ranks without changing snapshot invariants. Build 6c adds only the three FCT files for encoded particle energy/roundoff;
its nine transaction cases pass on two ranks. The combined two-rank regression
passes 74/74. Build 6b retains hash
`4c0472f3b2d83ea314d4329e45c3277e558b83c22a2cda9fe950410fb84be0bb`.
Build 6's original hash is
`4052d211343e4c714a5316c67e6f21c89748fdd30cd5aa083c2a8c44269ce0ac`.
The previous build-5 binary hash is
`ab9fb5608b437205adf0b208f0b29d60ecc23f685179f3d053c7eef84fc90916`;
those reports retain their original identity.

Current evidence is outside the source tree:

- `/tmp/srcsep-development.AilCx1/`: native source manifests, build logs,
  normalization, population, checkpoint, source and PNW reports.
- `/tmp/srcsep-validation-audit-y17wx199/`: MPI CLI results and independent
  warm-solver/provenance review.
- `/tmp/srcsep-warm-source-pinned/protocol.json`: reproducible study bundle,
  archived source/input hashes, resolved/unresolved coverage and plots.

These `/tmp` paths are session-local evidence, not a durable release archive.
Use the source-owned commands in [test/README.md](test/README.md) and
[the study README](validation/physics/README.md) to reproduce them. No release,
commit, push or observational qualification is claimed.

The bounded sanitizer aggregate passed its physics callbacks but then exposed
a source ownership assertion: a checkpoint test cleanup introduced a second
queue flush. Cleanup now asserts an empty queue; production retains the sole
flush owner. The focused WP31–WP41 rerun passes. The following documentation
gate then found generated Python bytecode; this work's cache was removed and
the documentation/warning/sanitizer gate now passes (optional STATIC01 skips
because this compiler does not advertise `-fanalyzer`). The same focused run
passes BIRTH 7/7 and warm mathematical/protocol 19/19. Keep these diagnoses distinct from
physical test failures and retain their original logs.

## Open scientific and original-plan requirements

- Select physical accuracy limits only after the requested numerical,
  finite-amplitude Lorentz and thermal-dispersion evidence supports them.
  High-k unresolved roots stay unresolved; no automatic support truncation.
- Complete native joint particle/wave evolution, finite-bin energy budgets,
  competing-hazard consistency and atomic positive-energy transactions.
- Establish the actual background/tube/ALE contract. Prescribed profiles do
  not automatically satisfy ideal-MHD wave work. The configured finite-source-
  radius Parker winding with purely radial U has a nonzero ideal-induction
  residual if B is held static (derived in the material staging document);
  the startup-history repair does not qualify it as an evolved material field.
- The analytic origin/speed, list-velocity ownership, static B/U/n histories,
  collective startup rollback and background-identity repairs pass source and
  build8 two-rank native gates.
  The subsequent provider-radius repair now passes build11 RADIUSN01–02 on
  two ranks, exercising both scalar and library paths with a translated origin.
  Unknown imported origins still return an explicit unavailable status.
- Complete the adopted Maxwell–Jüttner source bridge: proper density/field
  mapping, guiding-center cap flux, event counts and actual driver timing.
  The full-particle kernel is verified, but no truncated/resampled classical
  sampler or unproved gyroaveraged projection substitutes for that bridge.
- Extend actual PIC::TimeStep/source-driver qualification. Build11 HOST01–03
  pass the controlled two-step fixed-analytic cases through the real host
  driver; explicit source-epoch queries and completed shock/cache geometry
  correct the clean build10a failures. The compiled global time counter is ON;
  counter-OFF native execution, moving/thermal sources and additional ownership
  layouts remain unqualified. The core field-line manager's local-AMR/global-dt
  concern did not fail in these cases and has not received a speculative patch.
- Qualify live-four-state restart, supported MPI/thread layouts, migration,
  moving geometry, full manufactured evolution and the original IV06 coupling.
- The fixed-line API now passes separate-process continuation with explicit
  persisted flux caches. Qualify actual production resume entry points and
  decoded Cartesian-position/hazard consistency now passes build8; the
  live-four-state codec and output continuation remain unsupported.
- Remove inferred isotope metadata from adapters or require an explicit
  species source. The legacy nearest-mass/proton nucleon estimate is not a
  governing isotope contract; the parallel-diffusion adapter correctly keeps
  unavailable per-nucleon input unsupported.
- Complete equivalent-input numerical M-FLAMPA comparisons, actual SWMF replay,
  instrument-aware observations, untouched event cohorts and independent
  clean-checkout reproduction. The pinned upstream numerical M-FLAMPA regression now passes against
  those local references; it does not establish an equivalent AMPS experiment.

Consult [WAVE_EVOLUTION_CONTRACT_REVIEW.md](validation/physics/WAVE_EVOLUTION_CONTRACT_REVIEW.md)
for the governing wave-action/work equations, flux conventions and remaining
physical questions. No missing physical input is to be replaced by an assumption.
The source-owned [material staging design](validation/physics/MATERIAL_TUBE_WKB_STAGING.md)
and [Maxwell–Jüttner flux design](validation/physics/MAXWELL_JUTTNER_FLUX_DESIGN.md)
separate proposed native interfaces and unresolved reduction contracts from
implemented code.

The [SWMF producer audit](validation/physics/SWMF_FIELD_LINE_CONTRACT_AUDIT.md)
records the exact conditional branch mapping and distinguishes live PT Kelvin
input from archived MFLAMPA keV temperature. The user selected standard
BATSRUS total-energy exports and explicit division by two for the magnetic
view. That bounded read-only adapter passes build8; arbitrary producer and
composition contracts and full imported-provider qualification remain open. None of
these source findings enables unqualified import physics.

Build8 four-spectrum native ownership also passes PNW01–04 at N=128 and
N=129 on two ranks (4/4 each), with separate retained reports. These are
dynamic-storage/ownership gates, not global production evolution.

The earlier **build9** is a material-only three-file delta over
build8, with its own captured source changes and archived `amps.build9`.
The expanded material plus related native gate passes 31/31, exit 0, no SKIP.
Build8 remains independently archived and its 116-case/three-pair evidence
retains that exact source point.

The finite-amplitude closure study is archived as a separate prescribed-
circular-wave Lorentz experiment. It compares finite-time response with an
independent unperturbed-orbit reference, separates orbit/phase/k-mixture
refinement and preserves undefined relative error at zero response. It will
not label coherent-mode variance as Markov diffusion or select production
limits from these incomplete physical checks.

Build10 is a test-only six-file delta over build9: isolated HOST01–03 and
RADIUSN01 factories plus native registry/Make wiring and process exclusivity.
Its archived binary SHA256 is
`44206ff58475b7673d2ea9a7029436e41f8c6e898fd5d37d7cc00d0af86ecfbf`.
No source-epoch or radius-provider fix is included in that executable. Initial
harness failures are retained; build10a below provides clean fail-before results.

Build10a is the corrected **test-only** loader/trace delta, SHA256
`f44093105797048e8fc58db1f2ebfa8654d9e3ba8ceb29c9a76eaee3593714ac`.
It retains old physical producers. Clean native fail-before now executes
RADIUSN01/RADIUSN02 (both FAIL) and HOST01–03 (all FAIL), all on two ranks
with actual structured results and exit1. HOST cases each complete two real
PIC steps/eight first-birth moves, exposing only source X/coordinate mismatch
and lagged end-shock/cache endpoints. No core AMR/global-dt failure occurred;
no speculative core fix is made. Earlier input/trace harness failures remain
archived separately and are not counted as physical defect evidence.

Finite-amplitude response v2 now retains36controlled rows:31numerically
resolved and5phase-unresolved, using unchanged physical inputs, execution
code and numerical criteria with finer orbit/phase axes. All31,457,280updates
are archived with immutable parent lineage. Independent review passes the
lineage/coverage/error-account and Lorentz-oracle checks. These are conditional
convergence estimates: combined uncertainty/full response reaches6.13% in
one resolved case, and near-zero relative contrasts remain undefined. No
Markov scattering rate or production limit is inferred. The detailed record
is in `validation/studies/weak_alfven_domain.json`.

Build11 links from175captured native files plus23canonical source hashes,
SHA256 `82b824f27668b3cbb19161e6c6b1fe2c1b3d5ee17345cd3056a3459174994744`.
It includes source epoch/known-S fixes, ASE01–05, both provider-radius routes
and the tagged immutable startup loader. Actual **HOST01–03 pass 3/3** and
**RADIUSN01–02 pass 2/2**, all on two ranks against the archived `amps.build11`,
with native child exit 0 and retained structured reports. Each HOST case
executes two real host/PIC steps and records eight first-birth moves across
immutable per-rank traces. The radius cases exercise the actual translated
analytic initializer and both scalar and library coefficient paths. The
build10a failures remain archived as fail-before evidence. Reports:
`/tmp/srcsep-development.AilCx1/host-postfix-build11/summary.json`,
`radius-postfix-build11/summary.json`, and
`radius-library-postfix-build11/summary.json` under that same evidence root.
No core global-dt or counter-OFF native qualification is claimed. The combined
regression now passes **121/121**, exit 0, with no FAIL, SKIP or ERROR. All
three separate-process checkpoint pairs pass with byte-identical final
artifacts and all six native child jobs exit 0. Native CLI passes **9/9**;
production CLI software reports **18 PASS / 2 explicit SKIP**, exit 0. The
SKIPs remain post-allocation production evolution and actual coupled SWMF
execution. Their reports are in `combined-mpi-build11`,
`fresh-checkpoint-build11`, `native-cli-build11` and
`native-production-cli-build11` under the same evidence root. A post-gate
identity check rehashes all 175 native sources/build copies and 23 shared
sources (373 checks), with no mismatch. These records supersede the pending
build11 status without assigning production qualification to component tests.

## Published workspace build after core integration

All 67 reviewed delivery paths are now applied, including the required exact
external test fixture and the optional-context compatibility amendment.
`publication.json` under the integration evidence root binds every delivered
hash to the accepted candidate; `.git` index/refs remain unchanged. The source
backup is `pre-publication-core-backup`. The old generated build is archived as
`build.before-core-74974c2-20261117T223520Z` in this workspace, with settings and
the previous executable backed up in `workspace-pre-rebuild-backup` under the
evidence root.

The standard existing application selector was refreshed and `make -j16 amps`
passes. All existing `.amps.conf` settings, the original user Makefile/linker
edit and the compile-time profile are preserved. The actual workspace `amps`
and archived `sep-evidence-workspace/amps.workspace-core` both have SHA256
`dcad7ae35263e61cf123a8b7022235c337eb34231be30e4bd4180d0ff4b3f1c6`.
Its final nine native/CLI/restart/spectral campaigns all pass, with actual
zero child exits: combined **121/121**, HOST **3/3**, radius **2/2**, fresh
checkpoint **3/3** with byte-identical finals, native CLI **9/9**, production
CLI software **18 PASS / 2 explicit SKIP**, and PNW **4/4** at each of N=128
and N=129. The independent workspace audit verifies the raw child reports,
clocks, CLI controls, source/data/settings identities and all 67 published
hashes. Evidence: `sep-evidence-workspace/regression-summary.json` and
`sep-workspace-independent-audit/summary.json` under the integration root.
Core integration is complete; the full physical-development plan remains open.
