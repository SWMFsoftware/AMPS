# `srcMoon` development and validation status

**Status date:** 2026-11-12  
**Source revision inspected:** `09fabea55e367f07a530460b488fa7a981ff6176` plus preserved uncommitted Moon development  
**Purpose:** living implementation/evidence ledger, not a validation certificate.

**Current-core integration complete:** master and origin/master both point to
`09fabea55e367f07a530460b488fa7a981ff6176`; the preserved Moon development and
integration changes remain unstaged. No commit or push was made. All 1,193
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
| M3 chemistry | **NOT REACHED** | I21-I24 production-linked runners do not yet exist. Local Na photoionization exists; multispecies/electron-impact/charge-exchange qualification is incomplete. |
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
  `2761b447b46539edb21395550cd1b5d9ff6add4839b0db5a7e6d734cf2a4f02a`
- validation-data roadmap SHA-256:
  `d10eca1c539f9cf82940236d9e1aa36af665d43e82ba475925950d582fbf49c4`

The source copy adds dated implementation status, the user-approved
480/240/120/60 km 200-GB campaign, the 500-by-5000 bounded execution contract,
the non-qualifying 480 km physical-duration diagnostic, diagnostic stage
plots, M2 implementation notes, and partial Ar status. The
external copy retains the common scientific roadmap. Therefore:

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
  readiness must remain valid for the exact selected inputs.
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
