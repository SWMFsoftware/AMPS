# M1 implementation checkpoint — updated 2026-11-12

## Current-core integration complete — November 12

Master/origin/master are `09fabea55e367f07a530460b488fa7a981ff6176`; local
Moon development/integration remains unstaged. No commit/push was made.
All 1,193 core/application cases preserve the baseline dispositions:
919 PASS, 272 SKIPPED, one pre-existing strict-coverage FAIL, one pre-existing
CME timeout ERROR. All 919 baseline PASS identities remain PASS. The exact
ed875cef… I04 reference was separately approved/installed, older evidence
archived, and its strict comparison gives a zero-byte diff. M0 reconciles to
ten linked PASS. Normal-mode one-/eight-rank budget/sampling and eight-rank
Na/Ar supporting checks PASS; the final normal Moon/SPICE build explicitly
sets TESTMODE=off and -O3 -g. No run remains active.

Read STATUS.md section 5.9 and CORE_INTEGRATION.md before continuing. Durable
evidence/backups: `test_output/srcMoon/core-integration-20261112T224309Z/`.
Private workspaces remain under
`/tmp/amps-moon-core-integration-20261112T224309Z-ehx53mh8/`. Earlier checkpoints
below retain their source basis. M1 numerical qualification remains open;
formal M2/I14-I20 stay I13-gated. I04 review did not freeze an I13 reference or
authorize observational validation.

## M2 development continuation — November 12

The user authorized M2 implementation/local verification before M1 closure
and approved integrated finite-inventory release. Current work/evidence:
`test_output/srcMoon/m2-development-20261112T084925Z/`.

- `campaign-5/` records U08-U13 six PASS, formal I14-I20 seven SKIPPED,
  and supporting checks six PASS/one D02 SKIPPED, zero FAIL/ERROR.
- Eight-rank production I14 has 513/513 radial hits, 7/7 eligible-rank moving
  collisions and no invalid faces. I18 surface/source distributions PASS.
  I19 has exact count closure and maximum relative atom imbalance
  `1.066445872840518e-15` across 42 checkpoint/species rows.
- Fixed shared-vertex BVH misses, actual force-mover terrain routing and
  accelerated collision positioning; source placement and partial flight now
  share the actual departure epoch. Finite thermal/PSD release uses integrated
  competing first-order hazards and weighted packets, preserving analytic
  remaining stock. Operator-split temporal/weight convergence is still open.
- Dedicated binary/build identity is preserved under `build-8/`; do not use
  it with an unrelated generated configuration. New code/docs and the exact
  contracts are in STATUS.md section 5.7 and M2_IMPLEMENTATION_PLAN.md.
- D03 selected-input archive verification matches 214/214 official manifest
  entries. Reproducible qualification published scoped READY metadata and
  GCP/PRP processed products. Raw data and frozen June/December I17 selections
  were not changed. D02 remains PARTIAL: actual 18.6-year integration start/end
  epochs are unresolved; the question and relevant code are in question.md.
- Source roadmap is now
  `2761b447b46539edb21395550cd1b5d9ff6add4839b0db5a7e6d734cf2a4f02a`;
  external roadmap remains
  `d10eca1c539f9cf82940236d9e1aa36af665d43e82ba475925950d582fbf49c4`.
  Only the development-order exception changes work order; formal scientific
  gates are unchanged. The earlier statement below that no external-data file
  changed applies to the preceding Na-chemistry work package only.

Full shared-code M0 requalification and canonical restoration are recorded
under `m0-requalification/`; consult STATUS.md for their final result and the
post-restoration eight-rank Na ledger regression. Do not restart long M1
convergence runs automatically. M1 is incomplete; formal M2 is not qualified.

## November 12 continuation

The earlier stop instruction below is historical; the user resumed and
authorized this work. The current evidence root is
`test_output/srcMoon/m1-budget-audit-20261112T064009Z/`.

- The actual-event audit identified a production chemistry bypass: the
  configured generic Exosphere H/O lifetime table excludes Na and never calls
  the Moon lifetime callback. `input/moon.input` now dispatches
  `Moon::Chemistry::ProcessParticle`, with the same 130500 s lifetime and
  shadow law. One-/eight-rank linked sunlit/shadow survival and exact gas-event
  budget fixtures PASS. This does not qualify D04/I21 or Na+ daughter kinetics.
- The user explicitly approved the reproducible corrected I04 reference
  `4bb811b3e05d72fa708359f1bea0344efd1a9f7e89609d7e8a3ca7eade808cba`.
  It is installed with pinned provenance; the old reference is archived.
  `make -j test_Moon_check` PASSes with zero-byte `test_Moon.diff`.
- The corrected eight-rank 480 km diagnostic pilot finished PARTIAL after
  3601.63 s, 26022 steps, and 52 completed 500-step samples. Peak process-tree
  RSS was 9.245 GiB. Its 54 ledger rows independently PASS: exact zero count
  residual, weighted residual within the declared roundoff bound, no numerical
  removals. The final inventory is about 2.6951e27 neutral Na atoms.
- The pilot has no accepted duration candidate. Independence still fails for
  414/11/30 of 1175 observables at batch sizes 1/2/4. The independent Python
  reconstruction matches every reached production gate. Downstream
  equilibrium/CI diagnostics do not override the first failing gate.
- Figures and tables are under `duration-pilot/480km/diagnostics/`,
  `diagnostics-rates/`, and `stage_diagnostic.*`. The separately reprocessed
  rates figure divides actual interval totals by actual elapsed time so the
  short closing checkpoint does not appear as a sudden source/sink change.
  Earlier figures and their processor identity remain preserved.
- Full configuration-separated M0 requalification completed PASS in
  `m0-requalification/`, with outer log `m0-requalification.log`. It removes
  `build/` before each configuration, runs `make -j 'OPT3=-O3 -g'`, preserves
  each executable/log, runs the canonical nightly target under its own
  `nightly/`, and restores canonical Moon/SPICE. Its final summary records
  I01-I10 10 PASS, zero FAIL/SKIPPED/ERROR, with all nine required local
  prerequisites PASS. Across the appropriate neutral/ion builds, 15 distinct
  U IDs PASS and 11 later-capability IDs remain SKIPPED. The canonical nightly
  output exactly matches the approved reference; the diff is zero bytes.
  Restored root `amps` has the same `e82a204c...ad6780` SHA-256 as the pilot.
  No AMPS/MPI process remains. Formal I11 was not rerun by this baseline.
- `source-snapshot/` binds the pilot's dirty tracked/untracked source and
  compiled identities. `amps-budget-spice` preserves its executable. Do not
  confuse the later rebuilt root executable with the pilot identity.

Formal I11 still has 500-step samples, a 5000-total-step per-level cap and
disabled restart writes. The pilot is NOT EVALUATED for formal I11; I13
remains unfrozen and M1 incomplete. There are no observations in this work
package and no external validation-data file was modified. See STATUS.md
section 5.6 for current results and the unchanged roadmap hashes below.

## Earlier November 11 safe-stop checkpoint

**Resumed:** the user subsequently requested `continue`. The stop-state text
below is preserved as history. Current verification and campaign status are in
[`STATUS.md`, section 5.4](STATUS.md#54-resumed-verification--2026-11-11).
The shared-reference discrepancy and one-/eight-rank Na/Na+Ar storage checks
are now resolved. The bounded four-grid campaign finished in 37.14 minutes:
I11 FAIL, I12 PASS, I13 SKIPPED, zero ERROR; every mesh reached 5000 steps.
Peak aggregate RSS was 64.35 GiB. No AMPS/MPI process remains. Do not
interpret the historical nonzero diff or untested-scheduler notes below as
the resumed work package's current result.

The user requested: "stop when it is safe to restart with --approve-for-me".
Stop at this checkpoint; do not start a new long run until the user resumes.
The option text is the user's requested restart option, not a verified Codex
CLI recommendation. No processes need a scientific restart file.

## Checkout and scientific contract

- Checkout: `/data/vtenishe/Moon-test/AMPS`, branch `master`.
- HEAD: `22ad748432ce032775189ba234be7a3d339a4d8d`.
- Preserve the substantial existing dirty/untracked tree. No commit, stash,
  reset or blanket restore was performed. `Makefile` is generated and must
  not be committed. Existing root `handoff.md` is an older historical record;
  this file is the current M1 restart checkpoint.
- Read `AGENTS.md`, both roadmaps, the data-root `AGENTS.md`, STATUS.md,
  `test/ROADMAP_RECONCILIATION.md` and the affected I11/I12/I13 contracts.
- Source roadmap SHA-256:
  `08d53e71daa8e27bde364f4d46a82323dabeaf5b518440d9072491df1d93f18a`.
- Data roadmap SHA-256:
  `d10eca1c539f9cf82940236d9e1aa36af665d43e82ba475925950d582fbf49c4`.
- The source roadmap controls the approved 480/240/120/60 km campaign and
  500-step samples / 5000-step cap. Common scientific thresholds, reference
  review, milestone meanings and IDs remain binding. External raw data were
  not changed or used. SPICE hashes are in each linked manifest.
- Historical M0 is qualified. Current working-tree regression is FAIL and
  cannot inherit that qualification automatically. M1 remains incomplete;
  M2 remains EXPERIMENTAL / NOT QUALIFIED.

## Implemented during this work package

1. `src/pic/pic.cpp`: move collecting-cell buffer clearing immediately after
   the production buffer switch, outside the species file-output loop and
   output-cadence threshold. Clear once for all species; completed storage
   remains readable. Handle DIM=1/2/3 and null center nodes.
2. `srcMoon/MoonIntegrationTests.cpp`, `.h`, `main.cpp`: add
   `--i11-sampling-check` for I11, a nine-window constant-storage identity
   through the actual `PIC::Sampling::Sampling()` and completed callback.
   The fixture covers all-off / first-species-on / all-on flags, reports all
   compiled species and MPI-reduced exact errors. No particle advances or
   physical observations are supplied in this storage fixture.
3. Normal I11 now writes `particle_diagnostics.csv`: instantaneous Na particle
   count and statistical-weighted inventory from real owned particle lists,
   separate from the window-mean density integral. Full injected/escaped/
   numerical-loss event budgets are still missing and belong to the next
   diagnostic audit / I19. Do not call this a closed budget.
4. I11 writes I12 and I13 initialization receipts in the same allocated mesh.
   Python reuses the finest receipt only after matching executable/input/rank
   identities, then calls the ordinary independent scorer. I13 still requires
   passing I11 and investigator-approved frozen reference files.
5. `test/m1_execution.py`: concurrent isolated MPI launches, per-child logs
   and /proc histories, aggregate 200 GiB RSS guard, 8 GiB available-host RAM
   floor, 100 GiB disk floor and a default 3600 s parallel-launch budget.
   Stop/reap only this scheduler's process groups; preserve completed children.
6. `test/integration_common.py`: shared preparation, optional parallel I11
   execution, harvesting child errors, mandatory linked sampling prerequisite,
   shared initialization scoring. `run_all_test.py` forwards
   `--i11-parallel-levels` and `--test-runtime-seconds` to I11. Serial remains
   default. Scientific 500/5000 controls and scoring thresholds are unchanged.
7. Nine infrastructure tests now include child-failure harvesting, deadline
   cleanup and aggregate-RSS termination. Relevant READMEs, I11 acceptance,
   SAMPLING_CONTRACT.md, reconciliation and STATUS.md describe the changes.

## Preserved work-package evidence

```text
test_output/srcMoon/m1-sampling-repair-20261112T021011Z/
```

- `build/configure.log`, `make.log`: first clean SPICE/O3/debug build PASS.
- `sampling-1-rank/result.json`: PASS, nine windows, one species, one rank;
  completed-marker, collecting-buffer and callback errors are all exactly zero.
  Its executable predates the later instantaneous-inventory addition, which
  compiled successfully in the next clean build.
- `build/configure-final.log`, `make-final.log`: second clean SPICE/O3/debug
  build completed successfully before the nightly regression.
- `build/test-Moon.log`, `test-Moon-with-reset.diff`,
  `regression-with-reset/`: required four-rank `make -j test_Moon` result.
- `build/test-Moon-original-reset.log`, `test-Moon-original-reset.diff`,
  root `run_test_Moon/`: controlled rerun with original PIC reset placement.
  The corrected source was restored afterward. Both runs produce identical
  four microscopic-value differences and 927-byte diffs. The frozen reference
  was not changed. Make ignores DiffNum's failure and returns zero; check the
  diff itself. Current `test_Moon.diff` is therefore nonzero.
- `sampling-original-reset/result.json`: ERROR. The nightly no-SPICE control
  executable returned without the requested integration receipt, so this
  attempt does not prove the new fixture fails against the original code.
  Inspect its stdout/stderr and the `-application=Moon` versus `moon` build
  paths before using the nightly binary as a CLI integration adapter.
- `architecture-audit/test_summary.json`: PASS.
- Python test results: I11 infrastructure 9/9 PASS; global-driver self-tests
  4/4 PASS. These are infrastructure evidence, not M1 physics PASS.
- `build/configure-checkpoint.log`, `make-checkpoint.log`: final restoration
  of canonical Na-only Moon/SPICE/O3/debug build at user stop; completion and
  identity are recorded in `checkpoint.json` in the evidence root above.

Final checkpoint (2026-11-11 21:52 CST): clean restoration build PASS,
`APPLICATION=moon`, CSPICE `/home/vtenishe/SPICE/cspice`, kernels
`/home/vtenishe/SPICE/Kernels`, compiler override `OPT3=-O3 -g`. No AMPS
processes remain. Executable SHA-256:
`c9dd1585f9dbb5e1e4a9e847421f13fd5ae8f8e29615a7ce053e695092da211e`.
Source PIC reset file SHA-256:
`f9e6e3d72512b62b777e3435be834dc8eb3d03d57b8abad21a204c84bf166b67`.
The existing Python cache files under srcMoon/test predate this work and
were preserved. All Python checks in this work used `-B`.

Regression numbers (mean speed and LSO velocity components, m/s):

| quantity | current | frozen reference |
|---|---:|---:|
| mean speed | 2224.937 | 2229.331 |
| mean vx | -38.64898 | -21.92451 |
| mean vy | 11.29292 | 25.02429 |
| mean vz | 8.845359 | -32.30533 |

The reset relocation is isolated from this mismatch by the control run, but
the mismatch's underlying source/configuration cause is not diagnosed yet.
Do not silently refresh the reference to make the zero-byte requirement pass.

## Old duration pilot: final outcome

`test_output/srcMoon/I11-duration-pilot-20261101T0145CDT/pilot_result.json`
records ERROR, exit 137, MPI rank 2 SIGKILL on November 3. Approximately 49 h
wall time and 350 physical days were represented, with 2204 stored windows.
Its peak RSS was only 9.35 GiB; resource guards did not fire. The apparent
sampled inventory increase was contaminated by the now-identified stale cell
buffers and cannot establish a physical spin-up duration. There is no passing
I11 candidate from that pilot and no process to resume.

## Resume work in this order

1. Check the checkpoint receipt, source diff, clean build result and effective
   `.amps.conf`. No temporary rollback of `src/pic/pic.cpp` may remain.
2. Resolve the identical nightly reference differences in the control and
   corrected runs. Preserve both existing runs. Confirm a new
   `make -j test_Moon` has a zero-byte diff; restore canonical SPICE afterward
   using a clean generated tree.
3. Run the storage prerequisite on the final executable at one and eight MPI
   ranks, in new directories. Add/execute a separate Na+Ar compiled storage
   fixture to demonstrate multi-species layout coverage. Current fixture code
   covers every compiled species, but executed evidence so far is Na-only.
4. Inspect/test scheduler and companion-score reuse with real MPI launches;
   their unit tests alone do not establish a working concurrent campaign.
   In particular confirm I12 exact child-input derivation, manifest identities,
   per-stage plots and preservation of completed child results at termination.
5. Run the short corrected physical campaign and compare live particle-list
   inventory with the completed density mean. Audit source, photoionization,
   impacts/re-emission, escape and numerical losses if a real trend remains.
6. Only then run the one-hour M1 campaign below. If unchanged gates miss the
   5000-step cap, retain FAIL and request review of any proposed physical-
   duration contract. Do not weaken gates or automatically freeze I13.

Exact commands after resumption (use new output directory names):

```sh
cd /data/vtenishe/Moon-test/AMPS
python3 -B srcMoon/test/integration-tests/I11/run_sampling_check.py \
  --amps ./amps --input srcMoon/examples/m1_memory200.in --mpi-ranks 8 \
  --output-dir test_output/srcMoon/sampling-check-NEW

python3 -B srcMoon/test/run_all_test.py --skip-build \
  --amps ./amps --campaign-mode qualification \
  --m1-input srcMoon/examples/m1_memory200.in \
  --only I11 --only I12 --only I13 --mpi-ranks 8 \
  --i11-parallel-levels 4 --test-runtime-seconds 3600 \
  --output-dir test_output/srcMoon/m1-parallel-NEW
```

Hardware: 18 physical cores / 36 hardware threads, 251 GiB usable RAM. Four
eight-rank children use 32 hardware threads; do not describe these as 32
physical cores. Runtime has not been measured for this new concurrent path.
MPI execution needs permission outside the restricted network sandbox.
