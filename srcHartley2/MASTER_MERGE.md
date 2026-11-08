# Generic-comet later-master integration plan

This document prepares a future integration of the generic-comet development
without authorizing one. It is a review contract: it defines what must be
identified, merged, rebuilt, tested, and presented before the user can decide
whether `master` should change.

## Absolute stop and authority boundary

No instruction in this file is permission to switch to `master`, create an
integration branch from it, rebase, cherry-pick, merge, push, reset a ref, or
otherwise change either a local or remote master branch. The current work must
remain on `1d_comet_dev`. Before the first operation involving a selected
master target, report the proposed target ref and full commit, the development
range, the dirty-worktree inventory, and the intended operation to the user;
then stop and obtain explicit approval.

Preparation may use read-only Git commands to inspect refs and diffs. It must
not silently turn a read-only review into an integration attempt. Even after a
disposable integration branch eventually passes all tests, present the final
evidence and stop again for explicit approval before moving or updating
`master`. `READY_FOR_USER_REVIEW` therefore means “the evidence package is
complete enough for a decision,” not “approved to merge.”

## Current snapshot: focused evidence only

The following snapshot was recorded while Stage 10 was prepared. It is
historical evidence, not a target selection and not a substitute for
re-recording the repository state immediately before a future authorized
integration:

```text
development branch:       1d_comet_dev
focused Stage-9 tip:      bff4c532f34515365fd910691595f261db96847f
local master:             cb7be2bba41f373352ab853abed23955ea69c8b9
origin/master:            635d1622f80b445cf83130e90e5539bfa872e453
local/origin divergence:  1 local-only commit, 17 origin-only commits
merge base with local:    cb7be2bba41f373352ab853abed23955ea69c8b9
merge base with origin:   f9c2c00da9c8348b9d4003e2a82c1fe139b300cc
local-master range:       20 development commits
local-master diff:        163 files, 67,547 insertions, 631 deletions
```

The `bff4c53` snapshot passed focused Stage 0--9 implementation tests: 38 of
38 tests in 10 of 10 groups, with no failures or skips. It did **not** complete
the balanced independent-seed Stage-9 science campaign and is not an accepted
Rubin comparison. Subsequent Stage-10 source, documentation, tests, commits,
ref movement, or remote updates make every count and commit above
non-authoritative. Historical commit counts, diff statistics, test totals, and
hashes elsewhere in the documentation must never be copied into a future
readiness report without rechecking them.

## Select an explicit target before integration

Local `master` and `origin/master` have diverged. A bare reference to
“master,” “latest master,” or “the merge base” is ambiguous. Before any
authorized integration work:

- the user must choose the exact target namespace, such as local `master`, the
  currently recorded `origin/master`, or another named ref;
- resolve the chosen ref to a 40-character commit and record both values;
- record when and how any remote-tracking ref was refreshed; a stale
  remote-tracking ref must not be described as the remote repository state;
- verify that the selected commit still resolves to the recorded value
  immediately before the disposable integration branch is created;
- calculate the merge base, commit range, name-status diff, and diffstat
  against that exact commit; and
- if the chosen ref moves at any later point, stop. Create a new disposable
  integration attempt from the new commit and repeat conflict resolution and
  the full validation matrix. Never force, reset, or conceal target movement.

The readiness manifest must contain separate fields for
`target_ref_requested`, `target_commit_resolved`,
`target_resolution_timestamp`, `development_commit`, and `merge_base`. It
must reject a symbolic ref without its full resolved commit.

## Source ownership and integration waves

Integrate semantically on a disposable clean worktree only after the initial
user approval. Do not merge the entire development branch as an opaque unit.
Each wave below has an owner, dependencies, and a review boundary. A file that
appears to overlap two concerns must be reviewed in the earliest wave that
defines its interface and again in the later wave that consumes it.

### Wave 0: immutable inventory and evidence contract

Owner: integration and validation tooling.

- Freeze the selected target commit, development commit, merge base, ordered
  development commits, source allowlist, excluded-artifact list, and dirty
  state.
- Preserve compact checksums and summaries for accepted focused evidence;
  keep raw runs outside the source commit.
- Confirm that the unrelated root `Makefile` modification is preserved in its
  original worktree but absent from the integration diff.
- Confirm that no branch or physics code changes occur in this wave.

Wave 0 must be repeated if either selected endpoint changes.

### Wave 1: reusable chemistry and cometary-physics libraries

Owner: reusable core model libraries; no Hartley/Halley geometry or source
identity is allowed here.

```text
src/models/chemical_solver/ChemistryDiagnostics.h
src/models/chemical_solver/ReactionNetwork.h
src/models/chemical_solver/SystematicResidualAccumulator.h
src/models/chemical_solver/test_reaction_network.cpp
src/models/chemical_solver/README.md
src/models/cometary_physics/CollisionDatabase.h
src/models/cometary_physics/RadiativeCooling.h
src/models/cometary_physics/test_collision_database.cpp
src/models/cometary_physics/test_radiative_cooling.cpp
src/models/cometary_physics/README.md
```

Merge pure data structures and unit tests before any application consumer.
Preserve unit conventions, reaction stoichiometry and mass checks, collision
provenance, and the default behavior of applications that do not use these
libraries.

### Wave 2: dimension-neutral mesh and PIC foundations

Owner: AMPS core. These files must not acquire comet species, Halley
parameters, radial source laws, or application callbacks with non-null default
behavior.

```text
src/meshAMR/meshAMRcutcell.cpp
src/meshAMR/meshAMRgeneric.h
src/pic/pic.h
src/pic/picGlobal.dfn
src/pic/picGlobalDefinitions.dfn
src/pic/pic_block_send_mask.cpp
src/pic/pic_field_line.cpp
src/pic/pic_mesh.cpp
src/pic/pic_mesh_search.cpp
src/pic/pic_mover.cpp
src/pic/pic_parallel.cpp
src/pic/pic_particle_collision_models.cpp
src/pic/pic_pbuffer.cpp
src/pic/pic_photolytic_reactions.cpp
src/pic/pic_restart.cpp
src/pic/pic_time_step.cpp
src/pic/parallel/init_send_halo.cpp
src/pic/README.md
```

Review this wave by feature, not by whole-file replacement: DIM-neutral child
counts and searches; local time step; block/species particle weights; mover
boundary continuation; conservative split/delete logic; restart layout and
ghost publication; collision and photochemistry handling of individual
carriers; buffer preflight; and the Stage-8 pre-move/null-default callback.
Run unrelated DIM3/global-time-step/global-weight applications after each
interface-bearing subset. If the selected target already has an equivalent
hook or repair, reconcile with that interface instead of installing a second
implementation.

### Wave 3: input databases, grammar, and generated contracts

Owner: input configuration. Parsers describe policies; they do not execute
particle transport or silently choose scientific thresholds.

```text
ampsConfig.pl
input/CometPhysicsDatabase.pm
input/CometReactionNetwork.pm
input/CometPopulationPolicy.pm
input/CometCorrectionWindowPolicy.pm
input/CometDaughterPacketPolicy.pm
input/CometChannelImportancePolicy.pm
input/CometTotalReactionForcingPolicy.pm
input/hartley2.pl
input/species.input
input/hartley2.input
input/hartley2.radial.input
input/halley1d.input
input/test_comet_physics_database.pl
input/test_comet_reaction_network.pl
input/test_comet_population_policy.pl
input/test_comet_correction_window_policy.pl
input/test_comet_daughter_packet_policy.pl
input/test_comet_channel_importance_policy.pl
input/test_comet_total_reaction_forcing_policy.pl
```

Generated application headers are reviewed outputs of this wave and inputs to
later application waves:

```text
srcHartley2/CometChemistry.generated.h
srcHartley2/CometPhysics.generated.h
srcHartley2/CometDaughterPacketPolicy.generated.h
srcHartley2/CometChannelImportancePolicy.generated.h
srcHartley2/CometTotalReactionForcingPolicy.generated.h
```

Regenerate them with the selected target toolchain, compare them to the
reviewed policy, and record their hashes. Do not copy similarly named root
parser/input artifacts into source. Omitted switches must preserve legacy
defaults; malformed or incomplete importance policies must fail closed; and
the Stage-9 policy template must retain literal unresolved decisions until the
user supplies scientific values.

### Wave 4: generic-comet application foundation

Owner: `srcHartley2`; comet-generic geometry, sources, physics, and runtime
selection, with named inputs carrying comet-specific values.

```text
srcHartley2/Comet.cpp
srcHartley2/Comet.h
srcHartley2/CometModel.cpp
srcHartley2/CometModel.h
srcHartley2/CometPhysics.cpp
srcHartley2/CometPhysics.h
srcHartley2/CometRadiativeCooling.cpp
srcHartley2/CometRadiativeCooling.h
srcHartley2/RadialMain.cpp
srcHartley2/RadialMesh.cpp
srcHartley2/RadialMesh.h
srcHartley2/RadialMeshProfile.cpp
srcHartley2/RadialSimulation.cpp
srcHartley2/RadialSimulation.h
srcHartley2/LocalTimeStepValidation.cpp
srcHartley2/LocalTimeStepValidation.h
srcHartley2/main.cpp
srcHartley2/makefile
srcHartley2/test_radial_mesh.cpp
```

Preserve the explicit distinction between Hartley 2 and Halley input files.
The one-dimensional radial backend may specialize traversal and boundary
geometry, but the chemistry, collision, cooling, species, and source policy
interfaces must remain usable by future DIM2/DIM3 comet backends.

### Wave 5: population and chemistry estimators

Owner: generic-comet application algorithms. These files consume the core and
generated contracts established by Waves 1--4:

```text
srcHartley2/CometChemistry.cpp
srcHartley2/CometChemistry.h
srcHartley2/CometDiagnostics.cpp
srcHartley2/CometDiagnostics.h
srcHartley2/ChannelImportancePlanner.h
srcHartley2/ComponentCarrierPlanner.h
srcHartley2/CorrectionWindowPlanner.h
srcHartley2/DaughterPacketPlanner.h
srcHartley2/TotalReactionForcingPlanner.h
```

Review Stage 0 diagnostics, Stage 1 policy selection, Stage 2 component
transactions, Stage 3 residual accumulation, Stage 4 individual corrections,
Stage 5 windows, Stage 6 packets, Stage 7 channel importance, Stage 8 total
reaction forcing, and Stage 9 estimator arms in dependency order. Do not
weaken exact carrier accounting, rollback, source locks, capacity preflight,
continuation records, recoil, sink handling, or correction composition to
resolve a target conflict.

### Wave 6: tests and validation tooling

Owner: verification. Tests must exercise production paths and must not embed a
second implementation of the algorithm they claim to verify.

- Integrate every tracked `srcHartley2/tests/run_stage*_tests.sh`,
  `run_stage*_integrated_tests.sh`, and `test_stage*` source/script.
- Integrate `srcHartley2/tests/run_all_stage_tests.sh` only after all child
  groups exist; verify its plan, archived/live mode, counts, skip policy, and
  failed-log reporting.
- Integrate the tracked `srcHartley2/validation/` runners, evaluators, policy
  templates, self-tests, comparison exporter, and validation README.
- Integrate `write_all_stage_result.py` with the aggregate runner and verify
  the schema-v1 receipt against every retained log and `summary.txt`. Then
  integrate `stage10_release_policy.json` and
  `stage10_merge_readiness.py`; retain exact target-ref/commit/base/dev/tree
  binding, strict canonical-policy equality, Stage-9 science recomputation,
  unknown-path rejection, and immutable `NOT_GRANTED` authorization.
- Replay `run_stage10_tests.sh`, including receipt/log tampering, wrong
  topology, forged-science, untracked-source, generated-artifact exclusion,
  overwrite, and read-only-Git mutation cases. Never weaken the policy merely
  because a target worktree contains additional local artifacts.
- Integrate the tracked `input/test_comet_*.pl` parsers with their Wave-3
  modules, even though the aggregate runner consumes them here.
- Keep validation archives and generated plots outside source. A compact
  checksum manifest may refer to them, but must not pretend that a missing
  external archive is committed data.

Test totals from an older wave are not acceptance criteria. The final report
must use the totals printed by the exact integrated runner and must identify
every failed or skipped test and its log.

### Wave 7: references, analysis tools, and documentation

Owner: scientific provenance and user documentation.

```text
srcHartley2/halley2011/COLLISION_DATA_REVIEW.md
srcHartley2/halley2011/PAPER_PARAMETERS.md
srcHartley2/halley2011/derive_halley1d_controls.py
srcHartley2/halley2011/digitize_rubin2011.py
srcHartley2/halley2011/reference/halley1d_stage0_baseline.json
srcHartley2/halley2011/reference/rubin2011_neutral_profiles.csv
srcHartley2/halley2011/reference/rubin2011_neutral_profiles.metadata.json
srcHartley2/comet_shape_properties.py
srcHartley2/shou2025_spherical_parameters.py
srcHartley2/visualization/
srcHartley2/CHEMISTRY.md
srcHartley2/GENERIC_COMET_DESIGN.md
srcHartley2/README.md
srcHartley2/STATUS.md
srcHartley2/MASTER_MERGE.md
handoff.md
```

Paper digitizations and their metadata are reviewed source references, not raw
run output. Preserve their checksums and provenance together. Reconcile
append-only status/handoff history after code conflicts; never replace current
documentation wholesale with an older branch version. Update the live handoff
snapshot so it supersedes, but does not falsify, historical records.

## Semantic conflict invariants

Every conflict resolution must preserve the following behavior. A textual
merge that compiles but violates one of these constraints is a failed merge.

- Core defaults remain global time step, global particle weight, no individual
  correction, and null application callbacks unless an input explicitly
  enables the new behavior.
- DIM1 fixes use `2^DIM` or geometry-neutral iteration where required and do
  not hard-code one child or silently damage DIM2/DIM3 behavior.
- A particle is advanced to the current block boundary before ownership
  changes. Its remaining local-step fraction is rescaled by the old/new time
  step ratio exactly once.
- Boundary split/delete probabilities compose particle-weight and time-step
  ratios without changing represented flux, mass, momentum, or energy.
- Owner blocks are authoritative. New local weights, time steps, and
  application block data are published to ghosts before a consumer uses them.
- Adaptation uses powers of two, capacity preflight precedes mutation, rejected
  plans roll back completely, empty cells remain empty, and neighbor-carrier
  constraints prevent catastrophic boundary cloning.
- Emergency adaptation is requested from sampled population evidence and is
  executed at the established safe pre-exchange point, before receive-slot
  allocation. The callback default remains null for unrelated applications.
- Local-weight restart restores owner data, publishes it to ghosts, and retains
  the existing global-weight file behavior. No unsupported old local format is
  silently accepted.
- Chemistry preserves stoichiometric carriers, likelihood multiplication,
  parent sinks, product recoil and excess energy, continuation state, reaction
  rate, and exact rollback on allocation failure.
- Stage 6 packet size, Stage 7 channel probability, and Stage 8 total-reaction
  forcing are independent likelihood factors. Their product is applied once;
  no conflict resolution may double-correct or omit a factor.
- `configured`, `analog`, and `enhanced` Stage-9 arms use one executable and
  identical physical tables, layout, callbacks, mesh, collision physics, and
  source model. Only the estimator selection and independent seed may differ
  as declared by the campaign.
- Stage-9 safe checkpoints require no unfinished local-step fraction, no
  pending reaction event, an exact live-list count, and carrier
  `q = W_block,species * c_particle`; ESS is `(sum q)^2/sum(q^2)`.
- Missing provenance, convergence, collision feedback, buffer headroom, ESS,
  flux closure, profile agreement, or cost evidence is a failure, never a
  default, interpolation, or skipped science metric.

Never resolve conflicts by copying an entire old `pic_mover.cpp`,
`pic_restart.cpp`, `CometChemistry.cpp`, `RadialSimulation.cpp`, or
`input/hartley2.pl` over the selected target. Compare current interfaces and
apply the smallest semantic change that retains these invariants.

## Source allowlist and artifact exclusions

The eligible source set is limited to the tracked paths explicitly assigned
to Waves 1--7, their exact tracked test files, and Stage-10 review tooling and
documentation. Before staging, compare the proposed index against that
allowlist. A clean `git status` is not an allowlist, and an untracked file is
not source merely because it is useful locally.

The Stage-10 collector enumerates untracked paths. Any untracked path under
`src/`, `input/`, or `srcHartley2/`, and any root Markdown document, blocks
readiness because tests could otherwise execute code absent from the commit
range. Only canonical policy-listed generated/runtime families may be
summarized as exclusions; all other untracked paths fail closed.

Always preserve and exclude the current unrelated root `Makefile`
modification. It belongs to the user and must not be restored, overwritten,
staged, or represented as comet development. Also exclude:

```text
build/ and every other build tree
amps and other executables or object/library/module files
PT, PT/, PT.*, and PT* run trees
validation_runs/ raw archives and attempts
root generated halley1d.input, hartley2.input, hartley2.radial.input
root generated hartley2.pl and species.input
table_parser and other generated parser copies
comet.radial.mesh.msh, AMR binary measures, and generated meshes
run/error logs, profiler output, core dumps, and transient summaries
temporary Stage-9 enabled inputs and resolved local policy copies
failed, partial, interrupted, or pre-acceptance validation archives
SWMF_data, share/, output/, data/, and util/ runtime copies
.amps.conf, .general.conf, .last_input_info, and machine-local settings
editor swap/backup files, including `.CHEMISTRY.md.swp`
Python caches and temporary plotting/viewer files
```

Generated headers under `srcHartley2/` are exceptions only when their exact
tracked names appear in Wave 3, they are deterministically regenerated from
reviewed tracked input, and their hashes are recorded. Paper reference CSV,
JSON, and baseline files listed in Wave 7 are also source exceptions with
explicit provenance. These narrow exceptions do not permit other generated
files.

Raw archives remain external evidence. Record absolute path, size, purpose,
manifest checksum, and required-file checksums in the readiness package; do
not stage the archive. If an external archive is unavailable or fails its
checksums, report it as unavailable or invalid rather than reconstructing an
apparently equivalent result.

## Clean build and test replay

After all waves are resolved on an authorized disposable integration branch:

1. Verify the worktree path, branch name, selected target commit, staged
   allowlist, and absence of the unrelated `Makefile` edit.
2. Remove exactly that integration worktree's `build/` directory before the
   recompile. Never use an unresolved variable, home directory, repository
   root, wildcard, or broad recursive target.
3. Configure from the tracked canonical input. Check `.amps.conf` and the
   generated root parser/input copies so they cannot override or lag the
   selected input.
4. Record compiler, MPI implementation, build command, source commit, target
   commit, generated definitions, and hashes of the executable, effective
   input, `picGlobal.dfn`, chemistry table, and Stage-6/7/8 generated headers.
5. Run reusable library unit tests and parser tests first.
6. Run DIM1 radial mesh, mover, local-time-step, local-weight, restart,
   collision, photochemistry, capacity, and ghost-publication regressions.
7. Run representative target-owned DIM3/global-time-step/global-weight and
   non-comet tests selected from the target's actual test inventory. Do not
   invent names or silently omit unavailable target tests; record each choice
   and result.
8. Run every Stage 0--10 focused group through
   `srcHartley2/tests/run_all_stage_tests.sh`, with a new durable log root.
   Require zero failures and zero unexplained skips. Record the runner plan,
   modes, exact test/group totals, and hashes of the summary and all logs.
9. Re-run live MPI matrices required by the changed core/application paths,
   including the Stage-9 transported-daughter production-general-NTC fixture
   at 1, 4, and 12 ranks. Archived replay alone cannot validate conflict
   resolution against a new target.
10. Measure build/runtime memory and performance against the documented
    baseline using equivalent commands, ranks, inputs, and hardware. Explain
    rather than hide any regression.

`git diff --check` must cover the selected target-to-integration range as well
as unstaged changes. At the `bff4c53` snapshot that range reported a blank
line at EOF in `src/models/chemical_solver/README.md` and line-ending/trailing-
whitespace findings throughout the Rubin reference CSV. Normalize them only
as an explicit reviewed source change and update every affected checksum, or
record a narrow reviewed exception. Never weaken the check globally.

## Current Stage-9 science blockers

Focused implementation validation is complete at the historical `bff4c53`
snapshot, but Stage-9 science acceptance remains blocked by missing decisions
and evidence:

- the minimum carrier-ESS fraction for every claimed daughter is unresolved;
- the minimum variance-per-wall-time gain is unresolved;
- the required fraction of resolved radial bins meeting that gain is
  unresolved;
- flux closure has not been selected as either an independent-seed confidence
  interval containing zero with a documented null tolerance or an explicit
  positive relative tolerance;
- the distributed policy still contains `DECISION_REQUIRED` and therefore is
  not a runnable scientific policy;
- no balanced campaign with at least four fresh independent seeds in each
  analog and enhanced arm has been accepted;
- convergence, complete post-burn-in exposure, buffer headroom, collision
  feedback, daughter ESS, profile invariance, flux closure, and
  cost-normalized gain have not passed together under an approved policy; and
- authenticated final Rubin PNG, EPS, PDF, CSV, and JSON products from an
  accepted ensemble do not exist.

The accepted transported-collision matrix and one-window analog/enhanced
checkpoint smoke are focused contract evidence. The checkpoint smoke is
intentionally nonconverged. Neither archive supplies production policy values
or establishes physical agreement with Rubin et al. (2011). Values used by
synthetic mutation tests or broad validation fixtures must not be promoted to
scientific thresholds or tuned to the Rubin curves.

## Evidence required for `READY_FOR_USER_REVIEW`

Stage 10 may produce a report before the science blockers are closed, but
that report must say `FOCUSED_ONLY_NOT_SCIENCE_READY` and enumerate them. The stronger
`READY_FOR_USER_REVIEW` label is permitted only when every item below is
present and authenticated:

- [ ] The user-selected target ref and full target commit are recorded, and
      the ref still resolves to that commit.
- [ ] The full development commit, merge base, ordered commit range,
      name-status diff, diffstat, and dirty-state inventory are frozen.
- [ ] The proposed source index matches the explicit allowlist; unrelated
      `Makefile` and every generated/raw artifact are listed as exclusions.
- [ ] Every integration wave has a recorded reviewer disposition and all
      semantic conflict invariants pass focused tests.
- [ ] A clean build from the exact integration source completes after removal
      of only the exact integration `build/` directory.
- [ ] Executable, effective input, generated definitions, chemistry and policy
      tables, compiler/MPI environment, and launch commands have checksums.
- [ ] Reusable unit/parser tests, core regressions, representative target
      applications, all Stage 0--10 focused groups, and required live MPI
      matrices pass with exact counts and log hashes.
- [ ] Default-off behavior is demonstrated for non-comet and canonical comet
      configurations.
- [ ] Performance and peak-memory comparisons use equivalent configurations
      and disclose every material regression.
- [ ] The Stage-9 policy contains explicit user-approved thresholds and no
      `DECISION_REQUIRED` values.
- [ ] A balanced analog/enhanced ensemble with at least four independent fresh
      seeds per arm passes provenance, convergence, exposure, collision,
      buffer, ESS, parent/profile invariance, flux-closure, and
      variance-per-wall-time gates.
- [ ] The accepted campaign manifest, compact summaries, evaluator result,
      source/executable hashes, and raw-archive checksum manifests are frozen.
- [ ] Rubin comparison PNG, EPS, PDF, CSV, and JSON products are regenerated
      from the accepted summaries, mutually consistent, and authenticated.
- [ ] Performance, memory, uncertainty, closure, convergence, scientific
      limitations, failed attempts, and all exclusions are disclosed in the
      final review report.
- [ ] No target/master ref has been changed by the preparation or integration
      validation work.

If a science item is missing, the report must identify the exact item and use
`FOCUSED_ONLY_NOT_SCIENCE_READY`; it must not reinterpret the item as optional. A tool or focused
test PASS is never a waiver for missing science evidence.

## 3I/ATLAS campaign overlay for a later integration

The 3I inputs are an application-data consumer of the generic interfaces and
must be merged only after their dependencies, never as a reason to introduce
case constants into a reusable layer.

Merge order within the existing waves:

1. Before application consumers, merge
   `src/models/resonance_fluorescence/` and the targeted root `Makefile`
   wiring. This is a reusable core-data commit: it contains published H/OH
   models and immutable reference data, but no Halley, Hartley, 3I, mesh,
   species-index, MPI, or PIC dependency. Preserve unrelated target Makefile
   changes while adding the compile switch, include path, model build, and
   `libAMPS.a` object glob.
2. In Wave 1, merge the generic brightness math/configuration types. This
   includes the exact inverse of the core DIM1 spherical mover's local-radial
   velocity rotation; keep its round-trip test beside the shell-chord and
   Rayleigh tests. Keep the fixed-scalar path for backward compatibility and
   use named core model IDs, not duplicated H/OH tables or case constants.
3. In Wave 3, merge the strict brightness grammar and generic collision and
   reaction input machinery. Re-run its invalid-input tests before accepting
   any concrete 3I file. Enabled brightness must require explicit
   heliocentric-radial/LOS/impact geometry, accept either scalar radial or
   full-vector comet velocity, and require `ABORT` or
   `EXCLUDE_AND_REPORT` velocity-domain policy. Named H still requires solar
   Ly-alpha flux and Lowell OH still requires an effective wavelength.
4. In Wave 4, merge the DIM1 particle sampler and the single sampling-callback
   registration in `RadialSimulation.cpp`. Preserve the core callback order:
   sample every core sampling iteration; on core output, MPI-reduce, divide by
   `PIC::LastSampleLength`, write, clear, and synchronize. Demonstrate
   brightness omission/default-off with `halley1d.input` so the Rubin file
   schema and comparison remain byte-compatible. Keep model IDs, provenance,
   exact reference hashes, estimator, vectors, and policy as Tecplot dataset
   metadata. Output names must remain
   `pic.<species>.brightness.s=<species-index>.out=<output-index>.dat`.
5. In Wave 6, merge `srcHartley2/i3/tests/test_prepare_i3_cases.py`, the core
   fluorescence tests, and the expanded brightness parser tests. Run them
   with the reusable chemistry, collision, injection-distribution,
   brightness, and all-stage tests.
6. In Wave 7, merge `input/i3.conf`, the single generated `input/i3.in`,
   `srcHartley2/i3/prepare_i3_cases.py`, the compact I3 case/reference CSVs,
   and the I3 README as one reviewed provenance unit. Regenerate `i3.in` and
   require no diff. The guarded fluorescence acquisition utility
   and its authoritative snapshot already belong to the earlier core module;
   an upstream network refresh is a separate scientific-data change and must
   not be mixed into conflict resolution.
7. After the case tables and visualization utility are present, merge
   `srcHartley2/visualization/references/i3_brightness_MAVEN/`,
   `extract_i3_brightness_maven.py`, `plot_i3_brightness_comparison.py`, and
   `test_i3_brightness_reference.py` as one reviewed visualization/provenance
   unit. The committed CSV is a raster digitization, not original instrument
   data. Preserve its pixel-anchor and metadata files, verify the archive and
   member hashes, and require byte-exact CSV/metadata regeneration when the
   external archive is available. Runtime plotting must use the committed
   numerical table and `--reference-dir`; it must not regain an undeclared
   dependency on untracked `Archive.zip`.

The compile/runtime split is intentional: `i3.conf` is the review surface for
fixed executable layout and capability, while `i3.in` is the exact scientific
receipt for every run. During conflict resolution, do not recreate the former
six compile-time files or hand-edit one generated case section. Change the
reviewed CSV or generator, regenerate `i3.in`, review the single-file diff,
and rerun the generator test. A changed solar
dataset, JPL ephemeris, reaction coefficient, energy, collision source, or g
factor requires an explicit provenance update; do not retain the old hash or
silently treat a new value as equivalent.

Pre-merge I3 validation must include:

- verification of the supplied archive and embedded PDF checksums documented
  in `srcHartley2/i3/README.md`;
- all 12 numerical MAVEN/IUVS subsets pass schema, ordering, units, positive
  values, radial-bin and paired-uncertainty checks; the 87-row CSV and metadata
  regenerate byte-for-byte from the authenticated archive and reviewed pixel
  anchors;
- exact regeneration of the one combined `input/i3.in`, including all six
  named cases, on the integration commit;
- validation of all 1,201 core Lowell OH records as non-extrapolated and positive,
  exact coverage of `[-60,60] km/s`, and retention of the energy g factor as
  the authoritative downloaded value rather than only its 3090-A photon
  conversion;
- parsing all nine official IKL-2020 H coefficients, evaluating the analytic
  model at all six case velocity/irradiance pairs across its strict
  `[-200,200] km/s` domain, and checking the documented force-balance
  conversion from IKL `mu` to photon g factor;
- one clean preprocess and compile from `input/i3.conf` after removing only
  that exact worktree's `build/` directory; the same executable must load
  cases 1 and 5 (the low- and highest-Q mesh/population extremes) from
  `input/i3.in` without recompilation;
- parser/runtime selection validation for all six cases plus negative tests
  for an omitted selector, an unknown selector, and an unknown runtime key;
- one-/multi-rank pilot output with identical particle-resolved brightness
  columns after normalization by the same sampling length, and exact Tecplot
  species/index/output filenames;
- no missing collision pair and the expected Combi/Rubin/CHEMKIN provenance
  for all 28 unordered pairs among seven species;
- source-injection tests at 175 K for H2O, CO, and CO2 against the exact
  crossing-flux speed and cosine-angle distributions;
- an isolated parser assertion of the exact generated brightness species,
  model, effective-wavelength, three geometry vectors, comet velocity,
  velocity-domain policy, case condition, tau mode, cross-section, and
  threshold arrays/declarations; checking only the
  enabled flag or metadata is insufficient because a Perl regex-escape error
  can otherwise leave the compiled arrays at zero;
- chemistry validation that accepts a physically valid equal-`W/dt` network
  such as the I3 H2O family, while retaining the existing Halley Stage-4 log
  requirement for a positive unequal-carrier count so split/roulette
  coverage is not lost;
- a converged production run or a clearly labeled nonconverged status for
  each case that is presented scientifically;
- line-center tau review for every H comparison point, and no claim that a
  flagged thin result is a solved resonant-transfer result;
- a focused check that each g factor uses
  `e_r dot (V_comet + R(n)*v_local)` for the individual model particle, with
  the same deterministic local basis as the core mover and no bulk-only H
  substitution;
- an out-of-domain particle fixture for both supported policies: fatal under
  `ABORT`, and excluded-count/column accounting under
  `EXCLUDE_AND_REPORT`, with no clamping or extrapolation;
- comparison only within the finite 1e9-m mesh unless a separately reviewed
  domain change or outer continuation is implemented; and
- a fresh Halley build/test and Rubin comparison proving that I3 additions did
  not change the pre-existing default-off workflow; preserve its explicit
  `BrightnessMode=off` record.

Do not stage `Archive.zip`, raw Horizons/NRC/LASP downloads, `PT` trees,
generated root copies, executables, or build products. The deliberately
retained core
`reference/g_factors/oh_ax00_schleicher_ahearn1988/lowell_oh00_raw.jsonl` is
an exception: it is the compact audit source for the normalized OH table and
its non-extrapolation test, not a runtime artifact. The reviewed CSV data,
official IKL coefficient file, Lowell snapshot, generated immutable include,
and manifest hashes are source artifacts. When applying the root `Makefile`
patch, preserve the pre-existing `EXTRALINKEROPTIONS_F=-lmpi++` modification;
it predates this module and must not be attributed to it or overwritten. As
with all other work in this document, these instructions do not authorize a
merge, branch switch, cherry-pick, or master update.

## Two-layer runtime-input integration wave (I3 redesign)

This section supersedes every earlier instruction that describes six
`i3.caseN.input` executables. The redesign is reusable infrastructure, not an
I3-only file shuffle, and must be transferred in dependency order:

1. Merge the generic grammar layer first:
   `src/general/input/amps_input.h`, its INP-05 regression in
   `test_amps_input.cpp`, `src/general/input/README.md`, and the normative
   `src/INPUT_FILE.md`. Preserve `ResolveCase` as an owner-invoked helper;
   never make the generic parser guess which block owns a case or read argv.
   The `src/INPUT_FILE.md` case semantics are part of the public core input
   contract and **must be merged into master**, not left as comet-only notes.
2. Merge the `Config.pl` selection change independently. Verify exact
   precedence: `input/<application>.input` first, then same-named `.conf`,
   otherwise fatal. Preserve target-side configuration changes and remove
   stale `APPLICATION`/`InputFileAMPS` selectors without deleting unrelated
   user options. Re-run `perl -c Config.pl`, configure one historical
   `.input` application, then configure `i3` and verify `i3.conf` was selected.
3. Merge the application runtime bridge:
   `CometRuntimeInput.*`, the `RadialMain.cpp` `-i/-case` wiring,
   `RuntimeInputMode` handling and `CometRuntimeInput.o` registration in
   `input/hartley2.pl`, plus the narrowly scoped pre-network chemistry and
   pre-sampler brightness configuration APIs. Keep
   `Comet::RuntimeInput::Required=false` in source; only `i3.conf` changes it
   to true. That default is the compatibility boundary for Halley and Hartley.
4. Merge `input/i3.conf`, `input/i3.in`, the generator, generator tests,
   updated surface-injection runner, I3 README, comet README, STATUS and
   handoff notes as the application-data wave. Do not stage the root copies
   `i3.conf`/`i3.in` made by `Config.pl`; `input/` owns the source files.
   Delete the tracked former `input/i3.input` in the same integration commit
   so the `.input`-first rule cannot silently select the retired paradigm.
5. Regenerate `input/i3.in` and require a byte-empty diff. Then remove only
   the integration worktree's exact `build/`, configure `-application=i3`,
   and compile once. With that one binary require: no `-i` fails; `-i i3.in`
   without `-case` fails; an unknown case fails and lists 1--6; a valid case
   reaches a production-model validator; all six case sections resolve.
6. Reconfigure and cleanly rebuild `halley1d.input`. Run the maintained
   Halley chemistry/collision/Rubin regressions and verify there is no runtime
   input requirement and no brightness-schema change. This is mandatory
   compatibility evidence, not an optional application check.

Conflict-resolution invariants:

- `.conf` contains fixed species, dimension, compiled module/table topology,
  block/ghost layout and allocation maxima; `.in` contains physical run and
  observation settings within those capabilities. Do not move an observation
  case, source rate, case reaction rate, mesh target, random seed, run length,
  or brightness geometry back into `.conf`.
- Case blocks remain generic direct named children. Common assignments are
  copied first and one explicitly selected case overlays them. Do not add a
  default case, environment-variable fallback, numerical special case, or an
  I3 branch to `amps_input.h`.
- Runtime parsing remains strict and collective: unknown/missing keys fail,
  each rank parses the same path, canonical hashes agree before model state is
  installed, and case provenance is consumed explicitly.
- Chemistry topology/storage stays compiled, but case rate coefficients must
  be installed before the first `GetNetwork()` call. Brightness metadata must
  be installed before sampler registration. A conflict that initializes
  either singleton early invalidates the one-executable/many-case contract.
- Never stage `build/`, `amps`, `PT/`, root-staged input copies, logs,
  `Archive.zip`, or validation-run trees. No instruction here authorizes a
  branch switch, commit, cherry-pick, merge, push, or modification of master.

## Final user decision

When, and only when, all readiness evidence is assembled, present the exact
target and integration commits, complete diff, wave dispositions, exclusions,
test and performance results, independent-seed statistics, convergence and
closure results, Rubin products, and known limitations to the user. Then
stop. Do not switch to, fast-forward, merge, push, or otherwise modify
`master` until the user gives a new explicit instruction approving that exact
target and exact integration result.

If approval later names a different target or the recorded target has moved,
the approval does not transfer. Repeat the disposable integration and full
validation process against the new full commit and request approval again.

## Brightness physical-cell index repair

When transferring the particle-resolved brightness sampler, preserve the
separate index spaces in `SampleParticleResolvedBrightness()`. AMPS center
nodes include ghost cells, so `_getCenterNodeLocalNumber()` is valid for
`GetCenterNode()` but must never index `FirstCellParticleTable`; the latter is
dense over physical cells only. Merge the `centerNodeNumber` /
`particleCellNumber` split, its explicit bound, detailed free-list-alias
comment, and the configured-source parser assertion as one unit.

This is correctness-critical rather than a performance-only optimization. In
the I3 DIM1 build, physical cell 6 mapped to center-node index 8 while the
particle table had valid indices 0--7. The adjacent zero was interpreted as
the head of the free-particle chain, causing all 20,000,000 unused buffer
slots to be sampled as H particles. Later integration validation must include
an empty first sampling iteration, first/last physical-cell coverage, and a
short multi-rank completed sampling window that writes finite H/OH Tecplot
brightness files. Do not resolve a conflict by restoring the former shared
`cellNumber` index.
