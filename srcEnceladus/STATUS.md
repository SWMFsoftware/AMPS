# Enceladus forward-model status

Last updated: 2026-11-28 (America/Chicago)

## Merge with upstream `master` `6791b5a` — 2026-11-28 (second)

User authorization: pull upstream AMPS and merge if needed; nothing
committed. Fast-forward `5ce1335` → `6791b5a` (three SEP/SEP3D/nightly
commits; long-tier SEP tests on 16 ranks, routine `max_ranks` 16, new
`sep3d-event-spice` runner). Backup:
`test_output/upstream-merge-20261128b/backup/`. Only `MakefileTest/Table` and
`utility/nightly/runners.json` changed on both sides; both merged cleanly
three ways (the four Enceladus end-to-end blocks keep `Ranks=8` <= 16 and the
per-test `S_` runner note is kept). Upstream changed nothing under `src/`,
`srcEnceladus/` or `input/`, so the production executable is unchanged.
Checks: nightly tooling tests OK, Table check no problems, regenerated
`Makefile.test`/`.long` contain all four targets, test catalog 67/1342 PASS,
Enceladus application suite 117/117 PASS.

## End-to-end Table tests: E3 INMS and UVIS 2010 — 2026-11-28

User request: E3 and UVIS 2010 end-to-end runs in `MakefileTest/Table` that
load fitted plume parameters and plot the model against INMS and UVIS.
Decisions (user): `tenishev2014_hybrid_8` for both; nightly 5000 steps and
long tier 10000 steps; fitted sources in
`SWMF_data/PT/AMPS/data/input/Enceladus/fitted_plume` and observations in
`…/observation_files`; gate on an exact regression, plots and fit statistics
reported but never gating.

**Inputs.** Byte-identical copies (each directory with `README.md` and
`SHA256SUMS`) of the 2026-11-11 campaign's
`independent_each_observation` absolute H2O fits for E3 (3.142482e28 s^-1,
reduced chi-square 6.44) and UVIS 2010 E10 solar (4.576082e27 s^-1, shape
reduced chi-square 3.88), with manifests and the two geometry files their
headers cite, and of `plume_fit_v1` E3 H2O density and UVIS 2010 optical
depth plus ray table.

**Tests** (`srcEnceladus/test/end_to_end/README.md`): `Enceladus--E3-INMS`
and `Enceladus--UVIS-2010` (nightly, family `enceladus-e2e`, `${MPIRUN}`) and
`Long--Enceladus--E3-INMS` and `Long--Enceladus--UVIS-2010` (Tier=long,
`Ranks=8`, family `enceladus-long`; a long name may not begin with the nightly
test name `Enceladus`). Production build (`-amps-test=off`); decks copy the
production cases' physics, numerics and explicit mesh (E3 with its 20 km
tube); `check.sh` gates on completion, rank count, zone count, ledger closure
and `DiffNum.pl` equality of the observation and ledger products with
`output/test_Enceladus/<NAME>--*.ref_np8.gz`; `compare.sh` writes the
comparison PNGs and prints INMS chi-square/ratio and UVIS shape-fit/peak-scale
metrics. A perturbed one-particle ledger reference, a wrong zone count and a
wrong rank count each fail the gate.

**Reproducibility and `reproducible_parallel_run`.** Two 8-rank 40-step runs
were bit-identical, but two 8-rank 5000-step UVIS runs with the same seed were
not (resident 1,132 vs 1,135; photolytic losses 289 vs 295; columns up to
about 0.8% apart). Two load-dependent inputs were found: (1) the build enables
both core rebalancing paths (emergency, by measured latency, and at every
sampling break), which move AMR blocks by wall-clock time; (2) the unbuffered
particle exchange inserts received particles in MPI message-arrival order
(`MPI_Testany`), and particles draw random numbers in list order. With only
(1) removed, both E3 tests reproduced exactly but both UVIS tests still
differed after the first window (7 particles; columns up to 4%). The core
already sorts the cell lists by particle signature after each exchange, but
only in the nightly test mode (`-amps-test=on`), which these production-mode
tests do not use. Single-rank runs were reproducible. By user choice (keep 8
ranks; option 1), two core runtime switches were added, both defaulting to
the previous behavior for every application: `PIC::Parallel::DynamicLoadRebalancingEnabled`
(true; gates both rebalancing sites in `pic.cpp`) and
`PIC::Parallel::DeterministicParticleExchangeOrder` (false; also enables the
existing signature sort in the unbuffered and non-CUDA buffered exchange,
`pic_parallel.cpp`), declared and documented in `src/pic/pic.h`. The optional
Enceladus key `reproducible_parallel_run` (default false; parsed in
`ForwardModel.cpp`, allowed in schema-3 cases by `SourceCatalog.cpp`, applied
in `FinalizeNumericsAfterMeshAllocation`, shown in the receipt) sets both; the
four decks set it `true`. The runtime-parser cases U04/U04b also check default
false, `true`, and rejection of an invalid value.

**Verification** (evidence `test_output/enceladus-e2e-20261128/`; the four
tests built and staged through their own Table compile/rundir targets after
the native preflight, `stage.log`; runs through their Table run/check targets
with `MPIRUN`/`TESTMPILAUNCH` = `mpirun --bind-to none -np 8`, all four at once
on a shared node, `runs.log`):

- application suite (`srcEnceladus/test/run_tests.py`): **117/117 PASS**;
  nightly tooling tests **116 PASS**, Table check no problems;
- run A2 → references `output/test_Enceladus/<NAME>--{observation,ledger}.ref_np8.gz`
  (SWMF_data); run B2 against them: **4/4 empty diffs**, products
  byte-identical to A2; receipts `reproducible_parallel_run=on`, zero
  rebalancing events, 5 or 10 zones, ledger residual 0;
- run times at 8 ranks (compile ~1.7 min each): nightly E3 256-267 s, UVIS
  118-122 s; long E3 387-401 s, UVIS 187-191 s (Table `Time=` 15/8/20/12 min
  allow for a loaded node: up to 874 s was observed for nightly E3);
- comparison (not gating), last window: E3 INMS H2O chi2/sample 64.0 (5000
  steps) and 33.8 (10000), median model/observation 0.95 and 0.90 over the 16
  in-domain samples; UVIS 2010 shape chi2 9992 and 9922 over 1133 samples, model
  peak column 0.564 and 0.582 x the published 9.0e19 m^-2. Plots:
  `runB2.<NAME>/comparison/`.

**Notes.** The nightly publishes test logs, not run directories, so the
plots appear in the log only as metrics. `Makefile.test`/`Makefile.test.long`
(generated, git-ignored) were regenerated from the Table with
`utility/TestScripts/table_parser.cpp` to run the targets locally. The
legacy-test manual chapter (`doc/update_legacy_test_catalog.py --check`)
already fails on pristine upstream (`Mars-ions` lacks an editorial purpose)
and will need entries for these four tests when it is next updated. The
SWMF_data fixtures and references must be committed and pushed in SWMF_data
for the nightly, which tests a fresh SWMF_data clone.

## One receipt row per standalone test — 2026-11-28

User request: report one row per standalone plume-model test instead of one
row per module. `srcEnceladus/test/run_tests.py` now runs `standalone/tests`
through the new recorder `test/standalone_test_receipt.py` and writes one row
`S_<module>.<Class>.<method>` per test (the unittest selector without
`standalone.tests.test_`), with its own status, reason line, traceback and
duration. Module grouping and environments are unchanged: the four
exporter/publisher modules run together in the caller's environment and
still make up the registered aggregate `U00_exporter`; every other module runs
in its own interpreter with the data/SPICE roots. Rules: a unittest error is
`FAIL`; a failing subtest fails its test; a test not run because a class or
module fixture was skipped or failed inherits that outcome; a fixture failure
that explains no test (e.g. `tearDownClass`) is its own `FAIL` row; an
unimportable or empty module is one `FAIL` row `S_<module>`; a recorder crash
turns each of its modules into a `FAIL` row. All of these paths were exercised
on a throwaway fixture package (`--root`, not used by the runner).

The receipt grows from 55 to **117 rows** (14 `U`, 31 `R_`, 72 `S_`); the
registered 45 rows, their order and the global registry/catalog are
unchanged. `APPLICATION-ENCELADUS-FORWARD-001` accepts the dotted names and
reports `standalone_test_case_count`/`standalone_test_skipped`. The nightly
reader (`enceladus-receipt`) needed no change. Evidence under
`test_output/upstream-merge-20261128/`: `pertest-app-suite-roots` 117 PASS;
`pertest-app-suite-noroots` 107 PASS / 10 SKIP (the `S_spice_ground_tracks`,
`S_source_cases` and two `S_fitting` tests, which the per-module row had shown
as PASS); `pertest-global-forward` 45/45 PASS with 72 standalone rows; nightly
tooling unit tests 116 PASS. The nightly's history for the old ten
`S_<module>` rows ends when this change is committed; the new recorder file
must be committed with it.

## Merge with upstream `master` — 2026-11-28

User authorization: merge the local source with the current upstream core,
without committing. Local `master` again had no commits of its own (0 ahead /
15 behind), so `git merge --ff-only origin/master` fast-forwarded `74b15db` to
upstream `5ce1335` (2026-11-27); no merge commit exists and nothing is staged.
The same 35 tracked files remain modified, unstaged, on top of it.

**Preservation.** Every locally modified tracked file, the full local patch,
the before/target refs, the pre-merge `Makefile.def`, `AGENTS.md` and
`doc/AMPS.pdf` are saved under
`/home/vtenishe/Mars2/AMPS/test_output/upstream-merge-20261128/backup/`. Only
`Makefile` and `src/test/tests/harness/self/run_test.py` changed on both sides;
they were set to the committed base for the fast-forward and re-merged three
ways. No untracked file blocked the merge.

**Resolution.**

- `Makefile`: clean three-way merge (local `EXTRALINKEROPTIONS` and
  `FAST_WAVE_CXXFLAGS`; upstream `allclean` removes `Makefile.test.long`).
- `src/test/tests/harness/self/run_test.py`: upstream `cb91ebb` now pins the
  Moon cases itself (`[1220, 62, 50]` plus exact Moon identities); the local
  pin keeps upstream's form with `[1222, 70, 50]` for INP-06/07 and
  U09/U10/N06–N11. Its README gains a matching note; `src/test/README.md`
  names the new upstream base.
- The 2026-11-24 fix of `APPLICATION-ENCELADUS-FORWARD-001` (`S_` rows) is
  still required: upstream did not change that owner.
- Installed `Makefile.def` (git-ignored, written by `Config.pl -install`):
  upstream added `-include ${PTDIR}/Makefile.test.long` to `Makefile.def.amps`,
  after which the SWMF `Config.pl` refused every configuration ("Makefile.def.amps
  is newer than Makefile.def"). This was the first native-owner ERROR below.
  The existing site header (12 lines) was kept and the new template appended,
  which is exactly what install writes for that file; `-install` itself was not
  run because it would regenerate `Makefile.conf`. `Makefile.test.long` and the
  `table_parser` that writes it are produced only by a reinstall, so the
  Tier=long Table tests are not available from this checkout's `make`.
- Catalogs needed no regeneration (`--check` PASS, 67 suites / 1342 cases);
  `doc/AMPS.pdf` was rebuilt for upstream's `doc/ApplicationRegression.tex`
  change (PASS, snapshot 67/1342).

Upstream changed no `srcEnceladus/`, `src/general/`, `src/pic/` or registry
file; its core changes are SEP models, chemistry database, nightly tooling and
the Tier=long build hooks. No Enceladus model or physics change was needed.

**Validation** (evidence root
`/home/vtenishe/Mars2/AMPS/test_output/upstream-merge-20261128/`):

- Standalone discovery with data/SPICE roots: **72 PASS** (run twice).
- `srcEnceladus/test/run_tests.py` with data/SPICE roots: **55 PASS**
  (45 registered + ten `S_` modules), run twice.
- `global-focused/`: CORE-HARNESS-001 13, CORE-GENERAL-INPUT-001 7 and
  APPLICATION-ENCELADUS-FORWARD-001 45 PASS; the native owner ERROR (12) was
  the `Makefile.def` freshness refusal above, before any build step.
- `global-native/` (after the `Makefile.def` merge): native owner clean
  configure/hooks/`-j16` build and **12/12 PASS** (N00–N11).
- `global-forward-harness/`: **65/65 PASS, 3/3 suites**; the owner without
  data/SPICE roots records `S_source_cases` and `S_spice_ground_tracks` as SKIP.
- `production-build/`: preflight (limited to processes in this checkout;
  unrelated SEP/Earth jobs ran elsewhere on the node), `rm -rf -- build`,
  `Config.pl -application=enceladus
  -model-data-path=/data/vtenishe/cassini_validation_data/studies`, hooks,
  `env MAKEFLAGS=-j16 make amps` — PASS; `Makefile.conf` unchanged; the
  executable SHA-256 `1a6af0b7…` is identical to the 2026-11-24 build.
- `production/`: startup **6 PASS**; 1 READY / 27 MISSING_PUBLICATION;
  describe identical to 2026-11-24; local and four-rank validation PASS
  (fingerprint `1fdfc8af…`). The bounded E3 deck (`real-e3-short.in`) gives the
  same receipt as 2026-11-24: 40 injected, 39 resident, ledger residual 0.
- `APPLICATION-ENCELADUS-E3-001` (10,000-step extended campaign) was not run.

**Nightly coverage note.** The nightly's `enceladus` runner
(`utility/nightly/runners.json`) runs the committed upstream code, family
`app-enceladus` ("srcEnceladus (plume forward model)"): 53 rows, among them
`U00_exporter` (three exporter modules) and ten `S_<module>` rows. The local
`test_prepare_forward_cases.py` (untracked), the four data-root exporter tests,
U09 and U10 reach the nightly only after they are committed.

## Merge with upstream `master` — 2026-11-24

User authorization: fetch/pull upstream `origin/master` and merge it with the
local source, adjusting the Enceladus model if needed, without committing.
Local `master` had no commits of its own (0 ahead / 68 behind), so the pull was
a fast-forward from `d068476` to upstream `74b15db`; no merge commit exists.
All local work remains unstaged and uncommitted on top of it.

**Preservation.** Before the pull every locally modified tracked file, a
binary patch of all local changes (verified to reverse-apply), the 21 untracked
files the pull would overwrite and the before/after refs were saved under
`/home/vtenishe/Mars2/AMPS/test_output/upstream-merge-20261124/backup/`. The 13
files changed on both sides were set to the committed base for the
fast-forward and then re-merged three ways (base `d068476`, local backup,
upstream); 19 of the 21 untracked blockers were byte-identical to upstream.

**Conflict resolution.**

- `srcEnceladus/test/run_tests.py`: upstream `822dec0` (one `S_<module>` case
  per other standalone test module) is kept; `test_prepare_forward_cases`
  joins `EXPORTER_MODULES`, so it stays inside `U00_exporter` as documented.
- `src/test/tests/harness/self/run_test.py`: stream pin `[1222, 70, 50]` =
  upstream core 1220 + INP-06/07, Enceladus 62 + U09/U10/N06–N11, and Moon 50.
  Pristine upstream already failed `HARNESS-SELECTION-001` (observed
  `[1220, 62, 50]` against its pin `[1220, 62, 46]`) because `e5f2978`
  registered MOON-I21..I24 without updating the pin; verified on a
  `git archive` export of `74b15db`.
- `src/test/tests/applications/enceladus/run_test.py`: upstream's new `S_`
  rows made `APPLICATION-ENCELADUS-FORWARD-001` error on pristine upstream
  (exact-list check, no `SKIP` status). The owner now validates `S_` rows
  separately (unique, PASS/FAIL/SKIP, any FAIL raises a suite error), records
  skipped modules in its metrics, and still requires the 45 registered rows
  exactly and in order.
- `src/test/README.md` keeps both sections with merged totals;
  `doc/build_manual.py` is upstream's file, which already contains the local
  project Tectonic lookup; `TEST_CATALOG.md` and `doc/TestSuiteCatalog.tex`
  were regenerated (67 suites / 1342 cases, both `--check` PASS);
  `doc/AMPS.pdf` was rebuilt with `doc/build_manual.py --engine pdflatex`.
- Root `AGENTS.md`: the tracked upstream file (srcMoon V&V rules) is followed
  verbatim by the previously untracked local repository-wide file (CME program,
  read-only Git, native preflight), also verbatim.

Upstream's Enceladus changes (`74974c2`: portable repository-relative
provenance in exported source receipts and new exporter tests; `822dec0`:
per-module standalone cases and data/SPICE path options) needed no model or
physics change. No source law, reaction, transport, sampling or generic PIC
behavior was altered by the merge.

**Validation** (evidence root
`/home/vtenishe/Mars2/AMPS/test_output/upstream-merge-20261124/`):

- Standalone discovery with validation data and SPICE roots: **72 PASS**, no
  skips (66 earlier + upstream exporter tests + four data-root tests).
- `srcEnceladus/test/run_tests.py` with data/SPICE roots: **55 PASS**
  (45 registered + ten `S_` modules), no skips.
- Focused global campaign (`global-focused/`): native owner
  `APPLICATION-ENCELADUS-NATIVE-STARTUP-001` performed the clean
  configure/hooks/`-j16` build and passed all 12 cases (six legacy gates and
  N06–N11). The same run exposed the two pristine-upstream harness/owner
  regressions above (`HARNESS-SELECTION-001` FAIL, forward owner ERROR).
- After those fixes (`global-focused-rerun/`): CORE-HARNESS-001 13,
  CORE-GENERAL-INPUT-001 7 and APPLICATION-ENCELADUS-FORWARD-001 45 —
  **65/65 PASS, 3/3 suites**. Without data/SPICE roots in the owner
  environment `S_source_cases` and `S_spice_ground_tracks` are SKIP (recorded
  in metrics, not imported); both pass when the roots are supplied.
- Production-root rebuild (`production-build/`): preflight, `rm -rf -- build`,
  `Config.pl -application=enceladus
  -model-data-path=/data/vtenishe/cassini_validation_data/studies`, hooks,
  `env MAKEFLAGS=-j16 make amps` — PASS; `Makefile.conf` unchanged; the same
  four inherited warnings (dust format, `register`, two path-format sites);
  `srcEnceladus/test_output` held and restored.
- Production checks (`production/`): startup gates **6 PASS**; inventory
  exactly 1 READY / 27 MISSING_PUBLICATION; describe shows fitted
  `Q=3.142481780338366e28 s^-1`; local and four-rank collective validation
  PASS with no particle buffer. Bounded real E3 run (`real-e3-short.in`:
  4 steps, 2-step windows, 10 primaries/step, no tube): exit 0, 40 injected,
  39 resident, 1 escape, ledger residual 0, minimum edge 126050 m,
  `dt=210.0833 s` (CFL 0.5, 300 m/s). Execution/provenance only, not a
  converged plume.

The shared build was previously configured for `core-mesh-3d`; it is now the
production-root Enceladus build. The 2026-11-13 evidence directory
`test_output/enceladus-case-runtime-20261113/` cited below is no longer present
in this checkout; its paths are historical references.

## Implemented case runner and runtime AMR controls — 2026-11-13

This section supersedes the earlier proposal-only entry retained below. User
authorization: implement `CASE_RUNNER_PROPOSAL.md` and move domain/refinement
control into runtime `.in`. No fitting strategy is interpreted by Monte Carlo.
Git operations were read-only; source/data changes remain unstaged/uncommitted.

### Implementation

- New dependency-light `CaseSelection.{h,cpp}`, `SourceCatalog.{h,cpp}` and
  `MeshControls.cpp`; native application makefile includes all modules.
  The established core `ResolveCase` handles common/direct named overlays.
  CLI accepts `-i`/`-input`/`--input[=]`, `-case`/`--case[=]`, list, describe,
  local validation, and collective MPI validation. Unknown/duplicate/missing
  arguments and selectors are errors. Schema 2 explicit-file behavior remains.
- Schema 3 `source_set` resolves below `PIC::UserModelInputDataPath`, generated
  by existing `Config.pl -model-data-path=PATH`. Catalog/source/observation
  paths and hashes are checked without selecting an alternative fit.
  Absolute/traversing paths and escaping symlinks are rejected. All selected
  bytes are retained; native table readers use those verified snapshots.
  Actual ranks compare input/asset/executable identity before mesh allocation.
  MPI inspection explicitly verifies a null particle buffer and exits early.
- Runtime mesh keys define SI body-fixed box bounds, background cell maximum
  edge, optional hemisphere surface sizes, optional south-polar cone and
  repeatable named box/sphere/shell refinements. Requests combine by minimum.
  The exact block-intersection criterion catches small off-centre features;
  exact segment/AABB distance supports the existing flyby tube. Noncubic cell
  aspect ratio is handled explicitly. Explicit mesh law bypasses hidden
  TESTMODE scaling; absent new keys preserves legacy behavior. Compiled AMR
  limits are checked rather than clipping unattainable requests.
- Final numerical plan records actual minimum allocated edge,
  `dt=min(CFL*dx/v_characteristic)`, and `W=Q*dt/N` for each primary.
  Case inspection does not pretend these are known before mesh construction.
  Fresh case directories are required. Rank zero snapshots raw root/includes,
  canonical selected input, complete sources/assets and phase/iteration/rank/
  binary/source identity in `case_inputs/` and `run_manifest.json`; final
  receipts include realized cell size, timestep, limiting species and weights.
- Core input `Document::sourceBytes` appends exact bytes without changing the
  grammar/canonical form or existing leading fields. New dependency-free
  `amps_input_digest.h` supplies SHA-256, independently tested with known
  vectors. No CME physics, generic PIC particle algorithm, reaction constants
  or source emission law was altered by this task.
- `standalone/prepare_forward_cases.py` prepares a NEW immutable local bundle
  from explicit complete fit files and native observation blocks. It preserves
  nonuniform/zero strengths, validates full geometry/kinematics/mode coverage,
  merges disjoint species only when compatible, excludes nuisance backgrounds,
  and copies hash-bound observations/observer tables/trajectories. It refuses
  replacement of an existing bundle. “Publication” is local file packaging;
  there is no Git commit, network upload or public release.

### Real-data availability

The executable's final model-data root is
`/data/vtenishe/cassini_validation_data/studies`.
`input/enceladus_cases.in` declares all seven measurements × four bases.
Only one real prescription was already explicitly selected by the user:
`E3--tenishev2014_hybrid_8`, independent H2O, absolute
`Q=3.1424817803383656e28 molecules/s`. Its complete fitted export,
qualified `inms_e3_h2o` block/data and SPICE-derived E3 trajectory are copied
without changing original files into:

`/data/vtenishe/cassini_validation_data/studies/enceladus_sources_20261113`

The selection receipt is retained at
`test_output/enceladus-case-runtime-20261113/e3-source-selection.json`.
Other 27 real identities intentionally report `MISSING_PUBLICATION` until the
upstream selection is specified. No independent/pooled alternative was chosen
silently. Absolute E3/E5 physical CO2 remains deferred. A separately labelled
28-case synthetic bundle at
`.../studies/enceladus_case_validation_20261113_synthetic` exercises interfaces
only. Its geometry/rates are not inferred from Cassini and are not production
source prescriptions.

### Validation ledger

Evidence root (absolute):
`/home/vtenishe/Mars2/AMPS/test_output/enceladus-case-runtime-20261113`.
The native owner now serializes a clean configuration/hooks/top-level `-j16`
build, safely holds/restores the pre-existing 19 GiB application output during
source copying, and imports six unchanged legacy gates plus six new groups.
Machine tool/rank-flag settings remain configurable. MPI-unavailable groups
report SKIP; data/physics failures remain FAIL. New IDs are registered in the
manifest, case registry, traceability, generated `TEST_CATALOG.md` and manual
suite snapshot (56 suites / 1203 registered cases; registration is not PASS).

- Full standalone unittest discovery with
  `ENCELADUS_VALIDATION_DATA_ROOT=/data/vtenishe/cassini_validation_data`:
  **66 PASS**, no skips; `standalone-final.log`.
- Final dependency-light application campaign: **45 PASS**, no failures;
  `application-final/test_receipt.json`. U00 includes four new publisher tests;
  U09 includes nine C++ resolver-driven test methods; U10 checks independent
  refinement/distance/capability fixtures. Original 43 outcomes/references stay.
- Generic parser/digest tests: **7 PASS** (original five plus INP-06/07).
  Harness self-tests: **13 PASS**. See global run receipts below.
- Final combined report: **77 PASS / 0 FAIL / 0 SKIP / 0 ERROR**, four suites
  all PASS, 240.488 s. Global run:
  `test_output/enceladus-case-runtime-20261113/global-results/run-20261113T092040Z-2290774`.
  `amps-core-tests.json` and `.xml` are the single combined reports; per-stream
  core/application logs and exact per-test artifact paths are retained there.
  Native receipts contain all six legacy gates plus N06–N11. Nine transport
  runs cover all four bases and both primaries, including four ranks; two
  UVIS runs each output three rays in two windows. All 40-injection global
  ledgers close. Noncubic refinement gives minimum edge 63025 m and
  dt=7.878125 s for CFL 0.5 / characteristic speed 4000 m/s.
- Restored production-root native build: **PASS** through recorded root/idle,
  `rm -rf -- build`, Config.pl, ampsConfig.pl hooks and top-level `-j16` make;
  `production-build/{native-preflight.json,native-config.log,native-hooks.log,native-build.log}`.
  Site `Makefile.conf` remains byte-for-byte identical to preflight. The
  original application test output is restored. No `.amps.conf` was edited
  directly. The four existing build warnings (dust format, legacy register,
  two legacy path formatting sites) remain visible in the build transcript.
- Restored executable legacy startup: **6 PASS**; `production-startup/`.
  Production case inventory/describe: **PASS**, exactly one READY and 27
  MISSING_PUBLICATION, correct fitted Q; `production-list.log` and
  `production-describe.log`.
- Real E3 four-rank validation: **PASS**, null particle buffer;
  `production-mpi-validate.log`. Bounded four-step one-rank E3 H2O run:
  **PASS**, 40 injected primaries, closed final ledger, two sampling windows,
  complete source/input/observation snapshots, independently verified SHA-256,
  fitted Q=3.142481780338366e28 s^-1. Deliberately coarse test mesh minimum
  edge=126050 m; CFL 0.5 and 4000 m/s give dt=15.75625 s. Numerical overrides
  are confined to `real-e3-short.in` (4 steps, sample length 2, 10 primaries
  per step, no flyby tube). Original `input/enceladus_cases.in` remains the
  10,000-step/1000-window/100-primary configuration with the E3 tube.
  Output directory:
  `test_output/enceladus-case-runtime-20261113/production-e3-short-final`;
  log `production-e3-short-final.log`, receipt
  `production-e3-short-receipt.json`. Binary SHA-256:
  `54b776cea6298a128ff316b454d17ef09e5eb46d5421d2fef079a270dfcf94de`.
  This short run verifies execution/provenance, not a steady Cassini profile.

Earlier attempts remain visible rather than overwritten. The first native
legacy run exposed N02/N03 moving behind MPI initialization; restoring the
case-free pre-MPI installer fixed both unchanged gates. A first new UVIS test
looked for a nonexistent `instrument=UVIS` header; it was corrected to assert
actual documented UVIS title/event/schema, all ray rows and positive path
lengths. Global run `run-20261113T091105Z-2285726` recorded 74 PASS / 3 FAIL:
N08/N09/N10 checked newly added raw-input snapshots against generated sources
copied before that addition. A clean rebuild includes the final receipt code.
Run `run-20261113T091533Z-2287986` used a newly added context fixture that still
contained density column mappings with `comparison_species=none`; native
validation correctly rejected it. The fixture now removes those mappings,
preserving the original parser gate; the 45-result application rerun and
subsequent complete 77-check combined campaign pass.
The first private real-E3 short deck disabled the tube but retained its
inactive controls; unchanged validation correctly rejected that deck before
allocation. The final test deck removes the inactive controls;
`production-e3-short-input-rejected.log` and
`real-e3-short-disabled-refinement-attempt.in` retain the rejection evidence.
These are superseded development attempts, not omitted failure evidence.

### Scientific limits and operation

The bounded synthetic cases exercise all bases, both neutral primaries,
one/four-rank bookkeeping, both UVIS operators, output cadence and realized
runtime mesh controls. They do not establish observational agreement, steady
state, mesh/particle/seed/timestep convergence or physical CO2 calibration.
No new full 10,000-step fitting/validation campaign is claimed here. The
example retains the user's 300 m/s numerical characteristic speed and CFL
0.5; source gas speeds/temperatures remain the fitted file's values. These
controls require scientific convergence studies separately.

Complete key syntax, build/run commands, selection JSON and per-test purpose
are in `RUNTIME_INPUT.md`, `standalone/README.md` and `test/README.md`.
Automatic refinement from arbitrary UVIS/INMS paths remains future work;
explicit regions and the existing E3 tube are available now.

### Historical proposal record (superseded by implementation above)

### Proposed measurement/source-basis CLI and publication interface

`CASE_RUNNER_PROPOSAL.md` gives the detailed design for I3-style named cases
selected as `measurement--base`, including E3/UVIS examples, input/catalog
syntax, portable data roots, MPI startup, compatibility, implementation
stages and required tests. Inspection confirmed that `Config.pl
-model-data-path` already generates `PIC::UserModelInputDataPath`, and the
generic AMPS input parser already implements `ResolveCase`.

The forward model is to consume complete published source files, with no
fitting-strategy selector or interpretation. The standalone publisher must
map each case ID to its chosen source prescription and qualified observation
assets. Absolute files contain source rates; shape files retain explicit
user totals. The existing E3/E5 physical-CO2 deferral is preserved.

This change adds proposal/documentation only. No CLI/parser/application
implementation, native build, simulation or new validation result is claimed.
Existing explicit-file runtime inputs and all user-owned changes remain
untouched. The proposal is linked from `RUNTIME_INPUT.md`.

### Ready-to-run independent E3 H2O example

Added `input/enceladus_e3_h2o_fitted.in` for the independent hybrid-fit E3
H2O absolute source vector, with `Q=3.1424817803383656e28 molecules/s`.
It retains the user's 10,000/1,000-step, 100-particles/step, 300 m/s,
Courant-0.5 settings, active-Sun H2O parent loss and E3 refinement. It
selects H2O only and the E3 H2O observation; inactive CO2 controls and other
observations are removed. The default deck and fitted source files remain
unchanged. `./amps --validate-input --input
input/enceladus_e3_h2o_fitted.in` passes with the current native binary.
This was input validation only; the new long run was not started.
`RUNTIME_INPUT.md` gives the isolated-output run command and explains the
distinction between characteristic timestep speed and fitted gas velocities.

## Current outcome — automatic fitted sources and manuscript campaign

The complete requested manuscript fitting matrix has been repeated at:

`/data/vtenishe/cassini_validation_data/studies/fitting_compain_111126`

It contains 124 target fits: 48 base H2O/UVIS fits, 45 base CO2 fits and 31
complete 99-source candidate refits. All four source geometries are supported
by both automatic standalone export and the actual AMPS reader: hybrid 8+7,
Porco 98, Portyankina 41 and explicit Porco 99. The E17 search, window and
direction sensitivities, component ablations, cross predictions, quadrature
checks, geometry/candidate comparisons and orbital-proxy analyses were also
executed. The new E17 search selected the same coordinates as the historical
candidate; this remains a selected development hypothesis, not a vent discovery.

### Scientific decisions and export behavior

- Stripe vertices are projected radially to the 252100 m surface **before**
  fitting; the exporter serializes the same vertices and fitted kinematics.
  Library loading of old JSON is unchanged; the CLI's new projected default
  and explicit legacy/no-export options are documented in the standalone README.
- Pooled template training has separate nuisance backgrounds for every
  observation. Only physical source amplitudes are shared. Normalization and
  target-specific rate/background refits remain the manuscript's template
  procedure, not a newly introduced bilinear shared-shape inversion.
- `plume_fit.py`, `joint_fitting_study.py` and `co2_fitting_study.py` export
  automatically via `fit_source_export.py`. Each strategy/target receives its
  own manifest, complete canonical vector, resolved geometry, hash index and
  schema-2 source files. Exact zero components remain present; empirical
  backgrounds are excluded. No preferred strategy or new default forward
  source distribution is selected automatically.
- There are 124 conditional shape and 100 absolute files (224 total). The
  24 E3/E5 apparent-CO2 absolute exports are explicitly deferred as requested;
  no physical fragmentation factor in [2,5] was chosen. Apparent totals and
  conditional shape files remain available. Single-profile E3/E5 CO2 export
  follows the same deferral rule.
- The catalog's 360-degree azimuth endpoint is recorded as its identical
  0-degree direction in the half-open native syntax. Finite representable
  subnormal strengths are preserved, including the E18 candidate-fit tail
  near 3.8e-317. Overflow and underflow to zero remain fatal input errors.

### Numerical implementation and interpretation

Tilted catalog sources require negative shifted-Maxwellian speed tails.
The old analytic branch is retained for `a=s*cos(theta)>=-2`; for more
negative shifts, positive scaled Gauss-Legendre quadrature and inverse-CDF
bisection avoid recurrence cancellation. Independent moment/quantile tests
extend through shift -8 without changing the old numerical references.

The 98-source CO2 E3/E5-pooled TRF solve stagnated. Recovery solves the exact
same scaled NNLS problem with an independent active-set method, requires
objective nonincrease beyond roundoff and normalized dual/complementarity
KKT residuals <=1e-12, and retains the original failure and recovery losses.
No source prior, regularization, data reweighting or gate relaxation was added.
An exactly singular condition diagnostic is encoded as the strict-JSON string
`Infinity`; nonfinite rates and residuals are errors.

The report separates actual pooled training losses from post-normalization
target-refit scores. Its input-array audit also shows that independently
calibrated UVIS column proxies differ between 98/99 bases, so their pooled
objectives are not the same nested problem. One returned H2O INMS-only
99-source TRF result has normalized KKT residuals about 8.6e-10, above the
report's 1e-12 stationarity audit. That numerical limitation is reported, not
silently refitted or declared qualified. Large independent catalog totals
remain inverse-problem artifacts rather than physical production estimates.

### Validation evidence

| Check | Final result |
| --- | --- |
| Standalone mathematical/data/export tests | 62 PASS, no FAIL/SKIP |
| Global harness + Enceladus application | 56 cases / 2 suites PASS, no FAIL/SKIP/ERROR |
| Clean native Enceladus build | PASS, required root/idle/clean/configure/-j16 workflow |
| Native startup gates | 6 PASS |
| Every generated source through production validation | 224 PASS |
| Actual one-rank H2O/CO2 runs, four geometries | 8 PASS; 20 iterations, 200 injections, zero ledger residual each |
| Generated global catalog / manual snapshot | PASS; 56 suites / 1193 cases |

The final source-only catalog case additionally requires exact subnormal
preservation and rejection of range-losing literals. The global registration
adds `ENCELADUS-U03C-CATALOG-SOURCE-GEOMETRY`; the application-independent
receipt has 43 stable results. Native smoke runs validate source integration,
not statistically converged plume physics or observational agreement.

Local detailed evidence is under
`test_output/enceladus-fitting-111126/`. Final native build logs are in
`native-final-build/`; startup is `native-startup-final/test_receipt.json`;
all exports are in `native-all-exports-final/receipt.json`; evolving runs are
`native-catalog-final/<case>/<species>/receipt.json`. The final global receipt
is `global/run-20261111T073812Z-235/amps-core-tests.json`. Copies, original
failed attempts, commands, input products and software snapshots are retained
inside the external campaign's `validation/`, `execution/` and `provenance/`.
The old 19 GiB application test-output directory was moved only temporarily
during source-tree generation and restored intact. Git remained read-only;
unrelated user changes and site configuration were preserved.

### Delivered report and numerical products

The Word report, PDF preview, HTML, eight comparison PNG/PDF pairs and full
tables are in `<campaign>/report/`. The main document is
`Enceladus_fitting_campaign_111126.docx`. Tables include all 124 target rates,
7843 component rows including zeros, 20 actual pooled losses/KKT diagnostics,
224 export paths/hashes, and all 93 comparisons with rounded manuscript
Tables 2–5. Source files are in `<campaign>/{h2o,co2}/<case>/amps_sources/`.
`FITTING_CAMPAIGN_111126.md` and the standalone README give exact commands,
case definitions, software changes, provenance and reproduction instructions.

Known inherited limits are retained: steady collisionless/gravity-free
physics, digitized outbound E3/E5 water, common calibration systematics,
non-identifiable source partitions, approximate 2007 rays and 308 outward
2010 rays beginning at the tangent instead of the actual observer. This
campaign does not claim to repair those observation operators or qualify a
new long AMPS flyby simulation. The checked-in equal-rate startup sources and
default runtime deck are unchanged.

## Previous outcome — E7 density and E21 native-count observations

The runtime observation mechanism now distinguishes three contracts without
mixing units:

- response-qualified local number density (`comparison_species=H2O|CO2`);
- model-only geometry (`comparison_species=none`); and
- native instrument context (`comparison_species=none` plus the three
  `context_*` keys), which is displayed but never scored against density.

`data_file=input:PATH` resolves a repository-released table relative to the
runtime deck. Ordinary relative observation paths continue to resolve through
`observation_data_root`, so future external products require only a new named
`observation` block and no recompilation. Duplicate names, unknown keys,
missing columns, a mixed density/context contract, invalid uncertainties, and
nonfinite values fail before MPI. `RUNTIME_INPUT.md` contains the complete
syntax and units.

The checked-in `/home/vtenishe/Mars2/AMPS/input/enceladus.in` now requests 13
products:

- quantitative INMS H2O density for E3 and E5;
- quantitative INMS apparent-CO2 density for E3, E5, E7, E14, E17, and E18;
- model-only H2O/CO2 profiles on E13 geometry;
- E21 H2O and CO2 open-source neutral-beaming counts as separate native-count
  context products; and
- the 2007 zeta-Orionis and 2010 E10 UVIS occultation products.

### New observation releases and provenance

`standalone/export_forward_observations.py` reproducibly exports the added
tables from `/data/vtenishe/cassini_validation_data` and writes per-product
metadata plus `release_manifest.json`.

- E7: 43 mass-44 density rows from `Flybys2.mat`, matched count-for-count to
  official PDS Level-1A geometry. The optimized time alignment is 0.321 s and
  its maximum residual is 0.006 s. The point uncertainty propagates Poisson
  counting error through the local density/count factor; a separate 15%
  common calibration uncertainty is recorded but not folded into each point.
  CSV SHA-256:
  `1df70acf4549488aac187b0d783da01edd53c562be3181ec286b82cd0c461acb`.
- E21 H2O: 704 native count rows with interpolated supplied trajectory
  geometry. The 322 earlier count rows outside the supplied trajectory range
  are excluded rather than extrapolated. CSV SHA-256:
  `7b4507bf58b3015a23987815f68cc14b448cc4b57a3006f613ec97211ddbc292`.
- E21 CO2: 481 native count rows, all within the supplied trajectory range.
  CSV SHA-256:
  `7a8e7867bef97a16df6a98f969875b0fdf854334ac718a36531741847f44f2a5`.

The E21 values remain counts per integration with Poisson count uncertainty.
A quantitative AMPS comparison needs the model phase-space distribution
accumulated through the released velocity-, voltage-, and attitude-dependent
SIMION response. No scalar count-to-density factor was invented.

The exporter now forces LF-only records. The native validator caught the
initial CRLF final-header token before any run; the regenerated release passes
the strict C++ table reader on all three products.

### Output control

The new `volume_field_output=false` setting suppresses only the generic
full-volume active-species AMR dumps. Completed density buffers, window resets,
surface output, conservation ledgers, dense E3 flyby output, and all 13
observation operators remain active. This avoids approximately 14 GiB of
redundant H2O/CO2 field files across ten windows. A two-window native smoke run
proved that the collecting buffer is flushed between observation-only outputs;
memory also stayed constant at 2.889238 GiB throughout the long run.

### Native 10,000-iteration run

The current equal-per-source E3 distribution was transferred without
refitting. Its 15 vertically projected sources sum to
H2O=`2.060909950673409e28 s^-1` and apparent-CO2=`5.356e25 s^-1`.
The one-rank neutral run used 100 injected macroparticles per species per
iteration, the 26.260416666666668 s global timestep, H2O parent-loss-only
photolysis, no daughter particles, no plasma, no splitting/merging, 10,000
iterations, and ten 1,000-step sampling windows.

```text
cd /home/vtenishe/Mars2/AMPS/srcEnceladus
python3 test/run_e3_simulation.py \
  --executable ../amps \
  --input ../input/enceladus.in \
  --output test_output/srcEnceladus/e3_10000_extended_observations \
  --launcher direct --timeout-s 10800 \
  --expected-iterations 10000 --expected-sample-length 1000 \
  --expect-volume-fields off
```

The process completed in **2,521.6 s** and passed all **13/13** native gates.
Exactly 1,000,000 particles of each species were injected. Both ledgers close
with zero residual:

| Species | escaped | surface loss | photolytic loss | resident |
| --- | ---: | ---: | ---: | ---: |
| H2O | 990,189 | 3,234 | 368 | 6,209 |
| CO2 | 986,230 | 6,531 | 0 | 7,239 |

All 13 observation products contain ten zones. The dense E3 output contains
3,201 rows/zone and 186 in-domain rows/zone. Every in-domain E3 sample uses an
actual 15,756.25 m cell edge, below the requested 20 km maximum.

The complete 938 MiB evidence directory is:

```text
/home/vtenishe/Mars2/AMPS/srcEnceladus/test_output/srcEnceladus/e3_10000_extended_observations
```

Important retained artifacts and hashes are:

- `test_receipt.json`:
  `bef3784396a587b24bea65fa4ab50470c95a2c78f166479d413627c5dea5c2e3`;
- `enceladus.in`:
  `c48fa7c5ecd6b5e9cb312098a5e0bca3446c1aeb7bd730cc01121c084a7c423f`;
- `enceladus_e3_equal_h2o_co2_sources.in`:
  `eb30bcbcf5852b60ebdd072faf740e5d4f776f01561a8b80865dded11c70986d`;
- `forward_observation_release_manifest.json`:
  `0a15046d633381823f2e7519f1484c9e597193e665bcf03fc77cae0ff6377357`;
- `PT/plots/observation.<request>.dat` for each request;
- `PT/plots/comparisons/observation.<request>.comparison.png` for each
  request; and
- `PT/plots/comparisons/observation_comparison_summary.json`:
  `c7b1c21542f68dfbb34a8e4f1fd3b751474a4968ebfb3c466f859e4af2a43e32`.

### Comparison results

The plotter passed all **13/13** schema/zone/PNG gates. PASS means that a
consistent product was generated; it is not a claim that the unrefitted E3
source distribution agrees with every event. Final-window density metrics are:

| INMS product | in-domain rows | median model/observed | chi-square/row |
| --- | ---: | ---: | ---: |
| E3 H2O | 16 | 1.177 | 300.54 |
| E5 H2O | 13 | 0.499 | 89.97 |
| E3 apparent CO2 | 44 | 3.141 | 284.15 |
| E5 apparent CO2 | 44 | 0.390 | 226.83 |
| E7 apparent CO2 | 43 | 2.505 | 617.98 |
| E14 apparent CO2 | 371 | 1.665 | 721.00 |
| E17 apparent CO2 | 377 | 1.516 | 809.76 |
| E18 apparent CO2 | 366 | 1.266 | 588.71 |

E13 has 88 in-domain model-only rows. E21 H2O contains 704 count samples, of
which 530 lie in the current AMPS box; all 481 E21 CO2 count samples lie in the
box. E21 reports `display_only_no_density_or_count_response` and no goodness
of fit.

| UVIS event | samples | shape chi-square/sample | AMPS peak H2O column (m^-2) | published column (m^-2) | ratio |
| --- | ---: | ---: | ---: | ---: | ---: |
| 2007 zeta Orionis | 97 | 9.213 | 3.366e20 | 1.400e20 | 2.404 |
| 2010 E10 solar | 1,133 | 7.557 | 3.625e20 | 9.000e19 | 4.028 |

The mismatch is real evidence that one equal-source E3 distribution is not an
all-event fit. No observation was used to tune this rerun.

### Build, smoke, and regression evidence

The required native preflight was followed from
`/home/vtenishe/Mars2/AMPS`: the working directory and inactive process state
were checked, only `build/` was removed, `./Config.pl
-application=enceladus` was used without editing `.amps.conf`, production
hooks were regenerated with `./ampsConfig.pl -input enceladus.input
-no-compile`, and the top-level `env MAKEFLAGS=-j16 make amps` build passed.

The first hook-generation attempt exposed an existing packaging behavior:
the configured application source copy recursively included the prior 18 GiB
untracked `test_output` tree and filled the partial build. The process was
stopped, only the failed `build/` was removed, `test_output` was moved to a
sibling hold path during source-tree generation, and it was restored intact
before compilation. The repeated clean workflow passed.

- checked-in native input validation: **PASS**;
- `python3 test/run_tests.py`: **PASS, 42/42**;
- `python3 test/run_native_startup.py --executable ../amps`: **PASS, 6/6**;
- Python syntax checks for exporter, runner, and plotter: **PASS**;
- 200-step/two-window native smoke run: **PASS, 13/13**;
- smoke comparison generation: **PASS, 13/13**;
- native 10,000-step integration: **PASS, 13/13**;
- ten-zone comparison generation: **PASS, 13/13**;
- `git diff --check`: **PASS** after the final documentation audit.

### Current scientific limitations

- Only the E3 trajectory has a dedicated 20 km refinement tube. This run does
  not establish mesh convergence for E5/E7/E13/E14/E17/E18/E21.
- The finite box includes only 16/79 E3 and 13/79 E5 qualified H2O rows.
- E21 has no quantitative model-count response yet; a detector convolution
  must use phase space, voltage, attitude, and SIMION acceptance.
- E13 has no qualified H2O/CO2 measurement and remains model-only.
- UVIS columns include only the ray segment inside the finite AMPS domain.
  The 2007 ray is reconstructed from the tangent track, and 308 outward 2010
  branches lack a retained full Cassini vector.
- Published inventories identify additional 2005, 2011, 2016, and 2017 UVIS
  events, but the current local release lacks the qualified reduced spectra
  and explicit ray tables needed to add them without inventing data.
- The one-rank, 100-particle/species/step run has not passed particle-count,
  source-rate, mesh, timestep, box-size, or multi-rank convergence studies.

## Previous refined E3-only outcome

The requested refined E3 run is complete: 10,000 iterations in ten 1,000-step
sampling windows. All 13 run-level checks pass. A trajectory-following AMR
refinement capsule reduces the Cartesian cell edge at every in-domain E3
spacecraft sample from 126.05 km to 15.75625 km, an eightfold linear and
512-fold volume refinement. At each AMPS output the completed H2O/CO2
number-density buffer was evaluated along the frozen 3,201-point Cassini E3
SPICE trajectory and appended as one Tecplot zone. Model-only and INMS-overlay
PNGs were regenerated from the ten refined zones.

The coarse-mesh, effectively empty plateau from -14 through +10 s is removed:
the refined profile changes on one-to-two-second scales and no longer contains
the old 25-second, approximately 15-order-of-magnitude cellwise plateau. The
geometry is now resolved substantially better, but the density estimate is not
particle-converged. With the user-requested 100 injected
macroparticles/species/iteration and no splitting, 106 of the 186 in-domain
final-window H2O samples and 99 of 186 CO2 samples are exactly zero. In
particular H2O is empty from -10 through +7 s in the final window and CO2 is
empty from -8 through +4 s. These gaps are finite-particle occupancy, not a
resolved physical zero. Increasing spatial resolution alone therefore fixes
the cut-cell artifact but exposes Monte Carlo sparsity.

This remains a one-rank neutral-transport and H2O parent-loss photolysis run,
not a qualified observational validation or a mesh/timestep/particle-converged
production calculation. The direct data overlay has no instrument-response
convolution. The flyby sampler uses the compiled cell-centered-constant AMPS
stencil; it does not interpolate across neighboring cells.

The qualified evidence is:

- receipt: `test_output/srcEnceladus/e3_first_run_qualified/test_receipt.json`;
- native log: `test_output/srcEnceladus/e3_first_run_qualified/run.log`;
- AMPS fields and surface samples:
  `test_output/srcEnceladus/e3_first_run_qualified/PT/plots/`;
- application-independent receipt:
  `test_output/srcEnceladus/test_receipt.json`;
- native startup receipt:
  `test_output/srcEnceladus/native_startup/test_receipt.json`.

The retained superseded coarse long-run evidence is:

- receipt: `test_output/srcEnceladus/e3_10000_run/test_receipt.json`;
- native log: `test_output/srcEnceladus/e3_10000_run/run.log`;
- all AMPS fields and surface samples:
  `test_output/srcEnceladus/e3_10000_run/PT/plots/`;
- flyby Tecplot product:
  `test_output/srcEnceladus/e3_10000_run/PT/plots/flyby_density.E3.dat`;
- flyby PNG:
  `test_output/srcEnceladus/e3_10000_run/PT/plots/flyby_density.E3.png`;
- AMPS/INMS comparison PNG:
  `test_output/srcEnceladus/e3_10000_run/PT/plots/flyby_density.E3.inms_comparison.png`.

The current refined evidence is retained at the following absolute directory:

```text
/home/vtenishe/Mars2/AMPS/srcEnceladus/test_output/srcEnceladus/e3_10000_refined
```

- receipt: `test_output/srcEnceladus/e3_10000_refined/test_receipt.json`;
- native log: `test_output/srcEnceladus/e3_10000_refined/run.log`;
- exact runtime deck and equal-source startup receipt:
  `test_output/srcEnceladus/e3_10000_refined/enceladus.in` and
  `test_output/srcEnceladus/e3_10000_refined/enceladus_e3_equal_h2o_co2_sources.in`;
- refined mesh: `test_output/srcEnceladus/e3_10000_refined/mesh.dat` and
  `test_output/srcEnceladus/e3_10000_refined/mesh.msh`;
- flyby Tecplot product:
  `test_output/srcEnceladus/e3_10000_refined/PT/plots/flyby_density.E3.dat`;
- model-only PNG:
  `test_output/srcEnceladus/e3_10000_refined/PT/plots/flyby_density.E3.png`;
- AMPS/INMS comparison PNG:
  `test_output/srcEnceladus/e3_10000_refined/PT/plots/flyby_density.E3.inms_comparison.png`.

The complete 15 GiB AMPS field-output directory remains at
`/tmp/enceladus_e3_10000_refined_20261101`; the durable 171 MiB evidence set
above intentionally retains the mesh, log, receipt, inputs, Tecplot profile,
and PNGs rather than duplicating ten 742 MB volume files per active species in
the nearly full repository filesystem.

## Refined closest-approach E3 run (2026-11-01)

### Refinement and output implementation

Five fail-closed runtime controls were added to `input/enceladus.in`:

```text
flyby_mesh_refinement=true
flyby_refinement_target_cell_size_m=20000
flyby_refinement_radius_m=150000
flyby_refinement_time_min_s=-120
flyby_refinement_time_max_s=100
```

The mesh-resolution callback loads the selected segment of the same frozen E3
trajectory used by the observation operator and applies the target resolution
inside a 150 km-radius polyline capsule. The AMPS mesh builder tests a cell
characteristic diagonal, so the callback converts the requested maximum
Cartesian edge to its corresponding diagonal. Dyadic AMR then realizes an
actual 15,756.25 m edge throughout the in-domain track. The legacy south-polar
source refinement remains active outside this superposed capsule.

The runtime parser requires all four dependent refinement values when the mode
is on, rejects them when it is off, and validates positive lengths plus an
increasing time interval. The flyby Tecplot schema now records
`cell_max_edge_m` for each row. The integration runner requires every in-domain
value to be positive and at most the configured 20 km limit.

The refined mesh contains 3,165 leaf blocks and reaches level 5. Its minimum
edge is 15,756.25 m, compared with 126,050 m at level 2 in the superseded
coarse run. The global rule
`dt=min(Courant*minimum_cell_width/characteristic_speed)` consequently gives
`dt=26.260416666666668 s` for Courant 0.5 and characteristic speed 300 m/s.
Each 1,000-step sampling window spans 26,260.4167 s (7.29456 h), and the full
run spans 262,604.1667 s (3.03940 days) of model time.

The compiled species union remains unchanged, but generic volume-field output
is now enabled only for runtime-simulated species. This prevents inactive
neutral, ion, electron, and dust fields from consuming hundreds of MB apiece
on the refined mesh. The legacy zero-dust `mesh.data.dat` snapshot is also
suppressed only when the forward runtime is installed. Active H2O/CO2 mover,
sampling, chemistry, storage, and surface behavior are unchanged. Small
generic surface diagnostic files for the compiled species union remain.

### Build and staged run evidence

An initial clean-build invocation was not accepted because the configuration
generator was still running when compilation was started. No binary from that
overlap was qualified. The complete required preflight was repeated from
`/home/vtenishe/Mars2/AMPS`: the root and absence of active build/test processes
were confirmed, only `build/` was removed, and the application was regenerated
and built with the established workflow:

```text
cd /home/vtenishe/Mars2/AMPS
rm -rf -- build
./Config.pl -application=enceladus
./ampsConfig.pl -input enceladus.input -no-compile
env MAKEFLAGS="-j16" make amps
```

The final build passed. `.amps.conf` was not edited directly. Only inherited
warnings in legacy dust, Enceladus, and filename-formatting code appeared; the
modified forward-model sources compiled without a reported warning. Binary:

- `/home/vtenishe/Mars2/AMPS/amps`;
- 31,656,784 bytes, timestamp `2026-11-01 01:42:30.966129988 -0600`;
- SHA-256
  `1f3110e4447d77742f7563c683da17f67034cc6f8396308bb2f265e972499a02`.

Application-independent tests pass **41/41**. They now include parsed
refinement target/radius references and 17 typed negative configuration checks.
Receipt SHA-256:
`6c756838f7e607964d960d47d46678b2ae9f41a149e2a004d6327792eefbaa1b`.
Native startup tests pass **6/6**; refined-startup receipt SHA-256:
`84434d49d142258017a55df77e8596892e45d974219911230044dc3704163bcd`.

The first 100-step refinement pilot demonstrated the correct 15,756.25 m edge,
3,165-block mesh, and 26.2604167 s timestep, but it was deliberately stopped
during output after generic AMPS began writing the full refined field for each
compiled inactive species. That was an output-volume issue, not a solver
failure. After restricting generic volume output to active H2O/CO2, the fresh
pilot at `/tmp/enceladus_e3_refinement_smoke_active_20261101` passed all
**13/13** gates in 246.2 s. It injected 10,000 particles of each species,
closed both ledgers exactly, wrote only active volume fields, and verified all
186 in-domain trajectory points at 15,756.25 m.

The full command was:

```text
cd /home/vtenishe/Mars2/AMPS/srcEnceladus
python3 test/run_e3_simulation.py \
  --executable /home/vtenishe/Mars2/AMPS/amps \
  --launcher mpiexec \
  --expected-iterations 10000 \
  --expected-sample-length 1000 \
  --maximum-flyby-cell-edge-m 20000 \
  --timeout-s 14400 \
  --output /tmp/enceladus_e3_10000_refined_20261101
```

Result: **PASS, 13 checks, 0 failures**, elapsed wall time 5,235.4 s. The
process returned zero through ordinary MPI finalization, wrote outputs 0--9
with exactly 3,201 flyby rows and 186 in-domain positions in every zone, wrote
only H2O/CO2 generic volume fields, produced no plasma file, and verified
15,756.25 m as both the minimum and maximum in-domain flyby cell edge.

| species | injected | photolytic loss | surface absorption | outer escape | resident | ledger residual |
|---|---:|---:|---:|---:|---:|---:|
| H2O | 1,000,000 | 368 | 3,234 | 990,189 | 6,209 | 0 |
| CO2 | 1,000,000 | 0 | 6,531 | 986,230 | 7,239 | 0 |

The final native particle-buffer total was 13,448, exactly the sum of the two
independently traversed species lists. CO2 remained transport-only. H2O alone
underwent the selected Burger (2007) parent-loss channels, with daughters
discarded as configured.

### Refined flyby result and closest approach

The final 9,001--10,000 window has the following raw maxima inside the AMR
overlap:

| species | maximum density, m^-3 | E3 reference time, s | positive/186 cells |
|---|---:|---:|---:|
| H2O | 1.6602874525068922e14 | +27 | 80 |
| CO2 | 6.148649007991746e11 | +33 | 87 |

The nearest 1 s trajectory sample is at `t=0 s`, altitude 52.51995134 km; both
species have zero sampled particles in its containing 15.75625 km cell in the
final window. Selected points demonstrate the difference from the old coarse
plateau:

| E3 time, s | altitude, km | H2O, m^-3 | CO2, m^-3 |
|---:|---:|---:|---:|
| -14 | 111.497 | 0 | 1.438280e10 |
| -12 | 96.321 | 9.685010e12 | 0 |
| -11 | 89.421 | 1.245216e13 | 0 |
| -10 | 83.021 | 0 | 1.078710e10 |
| 0 | 52.520 | 0 | 0 |
| +5 | 62.937 | 0 | 3.595701e9 |
| +7 | 71.653 | 0 | 7.191402e9 |
| +8 | 76.890 | 5.534292e12 | 1.438280e10 |
| +11 | 95.762 | 2.628788e13 | 4.314841e10 |

The old values of order `1e-4 m^-3`/`1e-6 m^-3` from -14 through +10 s no
longer occur. Instead the finer mesh distinguishes plume-bearing and empty
cells at roughly one-second trajectory cadence. The remaining exact zeros must
not be interpreted physically: reducing the cell edge by eight reduces cell
volume by 512, while the fixed-particle timestep reduction increases the
resident macroparticle population only about eightfold. The expected local
occupancy is therefore much lower. A statistically resolved closest-approach
profile requires authorization to increase injected particles, enable a
qualified splitting scheme, lengthen/ensemble the sampling, or implement a
separately validated spatial observation estimator. None was assumed here.

The refined INMS overlay was regenerated with:

```text
python3 test/plot_flyby_density.py \
  test_output/srcEnceladus/e3_10000_refined/PT/plots/flyby_density.E3.dat \
  --output test_output/srcEnceladus/e3_10000_refined/PT/plots/flyby_density.E3.inms_comparison.png \
  --inms-h2o /data/vtenishe/cassini_validation_data/data/processed/plume_fit_v1/e3_h2o_sticking_corrected.csv \
  --inms-co2 /data/vtenishe/cassini_validation_data/data/processed/co2_fit_v1/e3_co2_density.csv \
  --closest-approach-offset-s -0.5
```

Result: **PASS**, ten zones, 186 in-domain positions, and 123 observational
rows loaded. The H2O observations remain outbound-only by product design, and
the mass-44 observable remains labeled apparent CO2. No source distribution
was fitted during this run: both species still use the equal 1/15 startup
fixture with E3 total rates inherited from the standalone observational
studies.

Artifact identities:

| artifact | bytes | SHA-256 |
|---|---:|---|
| run log | 2,254,848 | `cd986a979d8c3d9a26e19fa466ee8ee47699836d98e53cfb4000690bcf63a442` |
| run receipt | 4,775 | `81518157ed1bb0ad139b00325682c5e0d9be51cd4f2017fb09dcf8cd64a6f89b` |
| refined `mesh.dat` | 166,871,507 | `82b49cfa8693d0258d84f534c1cfacfd48a9331b08dbab4dbb67b51e537ff423` |
| refined `mesh.msh` | 716,945 | `c54a48b70ba5b4481341f7da8b34d58ce249cfb880e9ecb5550564711017637b` |
| runtime deck | 4,628 | `0b0e7c61579370ce6a081ef6cb5c49255b536b6d56bccee8d5be23d1dd84caca` |
| equal-source receipt | 15,874 | `eb30bcbcf5852b60ebdd072faf740e5d4f776f01561a8b80865dded11c70986d` |
| flyby Tecplot | 7,749,445 | `4a0770c56f6a441116f5d240efa9766584afc451e3a0d3029a3b93a6ad8f18e5` |
| model PNG, 2160x1512 | 389,007 | `a4a5753640769b53b9d48caf09fd91ebfd4d964ebb8237eaad7d115ac861069a` |
| INMS PNG, 2160x1512 | 465,977 | `104547f6743b6f51c4202971028351810526767369cedea40e2172e796d4cabb` |

## Superseded coarse 10,000-iteration E3 run (2026-11-01)

The checked-in `input/enceladus.in` now selects 10,000 total iterations,
1,000-step sampling, and the E3 flyby density observation operator. All other
physical/numerical controls are unchanged from the first qualified E3 run:
active Sun, neutral H2O+CO2, seven Burger (2007) H2O loss-only channels,
plasma off, 100 injected macroparticles/species/iteration, characteristic
speed 300 m/s, Courant 0.5, random seed 184467, and split/merge off.

Command:

```text
cd /home/vtenishe/Mars2/AMPS/srcEnceladus
python3 test/run_e3_simulation.py \
  --executable /home/vtenishe/Mars2/AMPS/amps \
  --launcher mpiexec \
  --expected-iterations 10000 \
  --expected-sample-length 1000 \
  --output test_output/srcEnceladus/e3_10000_run
```

Result: **PASS, 11 checks, 0 failures**, elapsed wall time 473.5 s. The
process returned zero through ordinary MPI finalization, produced no plasma
output, wrote outputs 0--9 with exactly 3,201 flyby rows per zone, and retained
exactly 186 in-domain positions in every zone.

| species | injected | photolytic loss | surface absorption | outer escape | resident | ledger residual |
|---|---:|---:|---:|---:|---:|---:|
| H2O | 1,000,000 | 425 | 3,738 | 994,987 | 850 | 0 |
| CO2 | 1,000,000 | 0 | 7,774 | 991,242 | 984 | 0 |

The particle-buffer resident total was 1,834, exactly the sum of the two
independently traversed species lists. The global timestep remained
210.08333333333334 s; each sampling window therefore spans 210,083.3333 s
(58.3565 h) of model time and the complete run spans 2,100,833.3333 s
(24.3152 days). This is model evolution time, not Cassini trajectory time;
every AMPS window is separately evaluated against the same E3 geometry.

### Flyby product and coverage

The trajectory contains E3 times -1600 through +1600 s at 1 s cadence. The
current Cartesian AMR root box covers only the 186 spacecraft locations from
-100 through +85 s. The Tecplot product retains all rows and marks coverage
with `in_domain`; outside-domain rows contain explicit zero placeholders and
are excluded from the PNG and physical analysis.

For the final 9,001--10,000 sampling window, the raw cell-sampled maxima
within that covered interval are:

| species | maximum density, m^-3 | E3 reference time, s |
|---|---:|---:|
| H2O | 2.3369410627079513e14 | +22 |
| CO2 | 9.020153492626351e11 | +22 |

These raw values are not convergence-qualified. The profile has abrupt cell
transitions and near-zero cell values around closest approach; with a
126.05 km minimum cell width and finite Monte Carlo population, the plot
diagnoses the present numerical resolution rather than validating an
instrument-resolved flyby profile. Mesh refinement, an expanded domain,
particle/seed studies, and an observational response operator remain required
before a quantitative physical validation against Cassini data.

### Closest-approach depression diagnosis (2026-11-01)

The large depression in both species is not a fitted physical feature. In the
final sampling window the decisive transitions are:

| E3 reference time, s | altitude, km | H2O, m^-3 | apparent CO2, m^-3 |
|---:|---:|---:|---:|
| -15 | 119.71 | 2.291543e12 | 6.966671e9 |
| -14 | 111.50 | 3.184311e-4 | 1.283412e-6 |
| 0 | 52.52 | 3.062876e-4 | 1.588486e-6 |
| +10 | 88.97 | 5.059817e-4 | 1.798881e-6 |
| +11 | 95.76 | 2.008343e13 | 4.651938e10 |

The effectively empty interval is -14 through +10 s. The entry at -14 s
coincides with the trajectory crossing the `x=-252.10 km` and `y=-252.10 km`
cell planes; the exit at +11 s coincides with crossing the `z=-252.10 km`
plane. The 16-order changes occur at cell boundaries, and the reported values
remain exactly constant for several consecutive spacecraft positions within a
cell.

The implementation comment in `ForwardModelAMPS.cpp` currently says "linear
center-node interpolation", but that statement is incorrect for this build.
`PIC::Mesh::GetCenterNodesInterpolationCoefficients` dispatches through the
configured coupler stencil, and `src/pic/picGlobal.dfn` selects
`_PIC_COUPLER__INTERPOLATION_MODE__CELL_CENTERED_CONSTANT_`. The selected
stencil assigns weight one to the containing cell. This code/comment and
documentation discrepancy was discovered after the retained run; no results
were regenerated or changed during the diagnosis.

The south-polar source geometry and solid-body occultation can produce a real
density reduction along this part of E3, but the current calculation does not
resolve its depth or timing. The temporary forward source file also contains
only the 15 surface components and omits the standalone H2O fit's constant and
radial background terms, so a shadowed coarse cell can become effectively
empty. An explicit spatial observation operator plus mesh/particle convergence
is required before assigning physical meaning to the depression.

### INMS comparison overlay (2026-11-01)

The model-only PNG and its run receipt remain unchanged. The new comparison
PNG was rendered from the same ten-zone Tecplot file with two locally released
E3 observational products:

| observable | released CSV | CSV SHA-256 | loaded | plotted in -100 to +85 s AMR overlap |
|---|---|---|---:|---:|
| sticking-corrected, INMS-derived H2O | `/data/vtenishe/cassini_validation_data/data/processed/plume_fit_v1/e3_h2o_sticking_corrected.csv` | `5233ac2f3bf6677839089fae2ac03edf5a9cf95bb34407192c848051f719fd05` | 79 | 16 |
| INMS mass 44, apparent CO2 | `/data/vtenishe/cassini_validation_data/data/processed/co2_fit_v1/e3_co2_density.csv` | `6c9ea9fd234d5ac5a169bf650bc36e1e6aa7f60f094894fdd505811a58955a8a` | 44 | 44 |

The corresponding metadata hashes are
`9f4bf56f92be2550e32852595f5ae24f0c69a7781eae1f9833cdca37ab0457aa`
for H2O and
`702dd6368ed03368fb9682270becb55832178c58389f3c966949aa25c7205e4c`
for apparent CO2. The 79-point H2O release spans +10 to +400 s from closest
approach; 16 points lie inside the present model coverage. All 44 CO2 rows lie
inside it: 37 have positive counts and seven are zero-count samples.

The absence of inbound H2O points in this PNG is therefore a data-product
selection, not an INMS non-detection. The selected fit-authoritative H2O
release is deliberately outbound-only because the raw mass-18 counts require
the time-dependent inlet/surface-sticking response; a scalar count-to-density
conversion was not substituted. By contrast, the selected closed-source mass
44 product covers approximately +/-50 s and consequently includes inbound
points. A separate digitized full-flyby published H2O trace is retained in
`standalone/figure_data/figure3_published_curves.csv` over approximately
-400 to +395 s. It is not merged with the sticking-corrected fit table because
it is a digitized published curve with different provenance and uncertainty.

The observational tables use time from closest approach. The plot applies
`t_reference = t_CA - 0.5 s`, because the independent 0.1 s SPICE search
resolved E3 closest approach at -0.5 s relative to reference UTC
`2008-03-12T19:06:12`. Thus the plotted H2O observations occupy +9.5 through
+84.5 s, and the mass-44 data occupy -49.975 through +48.743 s. This offset is
an explicit command-line argument and is not inferred from the data.

H2O error bars show the released independent digitization uncertainty; the
23% common density-scale uncertainty is stated below the plot but is not added
independently to every point. E3 mass 44 is labeled **apparent CO2**, not
composition-resolved CO2. Its error bars show Poisson counting uncertainty;
the 15% common calibration uncertainty and published factor-2-to-5 impact
fragmentation enhancement range are stated but not converted into pointwise
Gaussian bars. On the logarithmic plot, seven zero-count samples appear as
open downward triangles at the one-count sensitivity. Lower Poisson bars that
terminate at zero are omitted rather than extended to an arbitrary positive
density.

Exact rendering command:

```text
cd /home/vtenishe/Mars2/AMPS/srcEnceladus
python3 test/plot_flyby_density.py \
  test_output/srcEnceladus/e3_10000_run/PT/plots/flyby_density.E3.dat \
  --output test_output/srcEnceladus/e3_10000_run/PT/plots/flyby_density.E3.inms_comparison.png \
  --inms-h2o /data/vtenishe/cassini_validation_data/data/processed/plume_fit_v1/e3_h2o_sticking_corrected.csv \
  --inms-co2 /data/vtenishe/cassini_validation_data/data/processed/co2_fit_v1/e3_co2_density.csv \
  --closest-approach-offset-s -0.5
```

Result: **PASS**, ten model zones, 186 in-domain trajectory points, and both
observation series parsed and plotted. The comparison PNG is 2160x1512 RGBA,
292,394 bytes, SHA-256
`65702676de31f1d222886c4c22cf30de344e4c94d9e6bcb5551a0bdc7f995e43`.
The plotting utility is SHA-256
`fb9d0761e10abfc513f5105a94ed0526e262e16b8c87a6c84340232ef1c1db79`.
Python bytecode validation, the observational rendering, and a model-only
backward-compatibility rendering all passed. A negative CLI check omitted the
required closest-approach offset and correctly exited 2 before creating an
output. No AMPS rebuild or native rerun was needed because this change reads
the retained Tecplot and observational products without altering the
simulation.

This overlay is not an independent fit result. The forward file assigns the
same 1/15 relative strength to every one of the eight jets and seven stripe
components for both H2O and apparent CO2. Only the species totals came from
the E3 standalone studies: H2O from the independent E3 fit and apparent CO2
from the all-five-flyby shared-shape study. The strongly unequal fitted H2O
component fractions and the fitted CO2 shared shape were not transferred to
the forward run. Apparent visual normalization agreement, particularly for
outbound H2O, is therefore partly expected from reuse of the same observation
family and does not validate the equal-strength spatial source fixture.

Artifact identities and sizes:

| artifact | size | SHA-256 |
|---|---:|---|
| Tecplot flyby density | 6,981,187 bytes | `22199dfd18438e37442177737dcdd559670614490c59a935ea9ce01a6d50d876` |
| PNG (1800x1368 RGBA) | 183,956 bytes | `031959cba97ea9756991e2dde5a462f4a46fa99e6c7d3e9bf9a41d4ffa405d93` |
| INMS comparison PNG (2160x1512 RGBA) | 292,394 bytes | `65702676de31f1d222886c4c22cf30de344e4c94d9e6bcb5551a0bdc7f995e43` |
| native log | 2,256,569 bytes | `e8347db0fdd8c62fc329e9f865aa1c66e6fac89794f8f197ec9a2032481e694d` |
| run receipt | 4,488 bytes | `e770d3562ef917378f8e6227f5de6044a6c6a12625351b7428615bda4862a023` |

The complete run directory is 2.6 GiB because generic AMPS also writes field
and surface files for the compiled species union.

### SPICE provenance and independent reproduction

The suitable E3 kernels already present under `~/SPICE/Kernels` are:

- LSK `NAIF/naif0012.tls`, SHA-256
  `678e32bdb5a744117a467cd9601cd6b373f0e9bc9bbde1371d5eee39600a039b`;
- event PCK `CASSINI/cpck10Mar2008.tpc`, SHA-256
  `e1a4e70401b0c25abd1f9d7efa038c162596083f51c5fba8d736521b2169094d`;
- reconstructed SPK `CASSINI/080228AP_SCPSE_08050_08076.bsp`, SHA-256
  `b3eee97309e6ae59e0e00292c9d60aeb5b90c2040ba982b01157378475b1859e`.

The local CSPICE N0067 static library hash is
`53d74623fcfd22dcd5f9d01d1ebc4889f7b1034a5bd6738a681c1361b4d60ee1`.
`standalone/tools/spice_ground_tracks.cpp` was rebuilt and run with observer
CASSINI, target ENCELADUS, frame IAU_ENCELADUS, aberration `NONE`, reference
UTC `2008-03-12T19:06:12`, and -1600:1:+1600 s. The independently regenerated
CSV and checked-in trajectory both hash to
`e5707e456fb03586a36958524007394c5238b5e48fbfa2db3bf41030ace5e663`.
AMPS itself remains built without SPICE; it consumes this frozen, independently
reproducible trajectory and therefore does not silently select kernels at run
time.

### Implementation and verification changes

- Runtime keys `flyby_density_output`, `flyby_id`, and
  `flyby_trajectory_file` select the observation operator. Dependent controls,
  unsafe IDs, malformed headers/rows, mixed flyby IDs, non-finite values,
  negative altitude, and non-monotonic times fail before MPI.
- The trajectory file is parsed transactionally during runtime installation.
  At each ordinary AMPS output a registered external sampling callback reads
  the completed number-density buffer, uses the build-selected native
  cell-centered constant stencil, enforces one AMR owner per point, reduces
  across MPI, and appends a Tecplot zone with SI values and coverage flags.
- `test/plot_flyby_density.py` strictly parses the zone schema and renders only
  in-domain values. Its optional INMS readers validate E3 identity, finite and
  strictly increasing values, the mass-44/apparent-CO2 interpretation, and
  count/density consistency before plotting. `test/run_e3_simulation.py`
  checks requested iteration, sample, and particle totals plus all Tecplot
  zones and the model-only PNG artifact.
- Application-independent tests now pass **41/41**; receipt SHA-256 is
  `6c756838f7e607964d960d47d46678b2ae9f41a149e2a004d6327792eefbaa1b`.
- Native startup tests pass **6/6**; refined-startup receipt SHA-256 is
  `84434d49d142258017a55df77e8596892e45d974219911230044dc3704163bcd`.

The first clean rebuild exposed and stopped on an incorrect namespace for the
generic interpolation-stencil capacity. After correcting that compile-time
API name, the full mandatory preflight was repeated from the repository root:
process check, `rm -rf -- build`, `./Config.pl -application=enceladus`,
`./ampsConfig.pl -input enceladus.input -no-compile`, and
`env MAKEFLAGS="-j16" make amps`. The second clean build passed. `.amps.conf`
was never edited directly. The resulting executable is 31,591,216 bytes,
SHA-256 `5b7643cebb05129058cec10906373c00bd493309e1bc9cb01841a5f6cd8cafc3`;
the long-run deck hash is
`733f6fbaa0e307e15c19453365fe4ef37263db958486d75c137afe431c7da0bc`.

## Qualified E3 run configuration

The following section records the retained, earlier 800-step qualification
run and should not be confused with the current 10,000-step result.

The run deck used for that retained qualification was complete rather than a
fail-closed template. It used:

| control | value |
|---|---|
| source case | `tenishev2014_hybrid_8`, E3 absolute receipt |
| primary species | neutral H2O and neutral CO2 |
| chemistry | seven Burger (2007) H2O photo channels |
| reaction products | discarded; parent loss only |
| solar activity | active |
| plasma | `none`; allocated coupled-background arrays explicitly zeroed |
| characteristic speed | 300 m/s for H2O and CO2 |
| Courant number | 0.5 for H2O and CO2 |
| primary particles | exactly 100/species/iteration globally |
| sample length | 800 steps |
| total length | 800 steps, one sampling window |
| random seed | 184467, folded with MPI rank for native RNG |
| split/merge | off |

After AMR allocation, the global rule
`dt=min(Courant*minimum_cell_width/characteristic_speed)` produced
`dt=210.08333333333334 s`, limited by H2O (the two configured limits are
equal). The implied minimum allocated cell width is 126,050 m. The 800-step
window represents 168,066.6667 s (46.6852 h).

The fixed primary statistical weights printed by the executable are:

| species | absolute receipt rate, molecules/s | macroparticle weight |
|---|---:|---:|
| H2O | 2.0609099506734099e28 | 4.3296283213730557e28 |
| apparent CO2 | 5.356000000000001e25 | 1.1252063333333338e26 |

The CO2 rate is the fitted mass-44 apparent rate. The retained factor-2--5
impact-fragmentation interpretation brackets physical E3 CO2 at
`1.0712e25--2.678e25 molecules/s`; no value inside that interval was selected.
Consequently this startup run must not be cited as a physical-CO2 inference.

## Run results and conservation

Command:

```text
cd /home/vtenishe/Mars2/AMPS/srcEnceladus
python3 test/run_e3_simulation.py \
  --executable /home/vtenishe/Mars2/AMPS/amps \
  --launcher mpiexec \
  --output test_output/srcEnceladus/e3_first_run_qualified
```

Result: **PASS, 9 checks, 0 failures**, elapsed wall time 40.0 s. The native
loop completed all 800 steps, every rank reported 800 injection calls, and the
executable terminated through the ordinary MPI-finalization path.

| species | injected | photolytic loss | surface absorption | outer escape | resident | ledger residual |
|---|---:|---:|---:|---:|---:|---:|
| H2O | 80,000 | 28 | 288 | 78,896 | 788 | 0 |
| CO2 | 80,000 | 0 | 568 | 78,513 | 919 | 0 |

For each row, `injected = photolytic + surface + escape + resident` exactly.
CO2 is transport-only and recorded zero photolytic losses. No
`ices.data.dat` was produced. The combined resident count of 1,707 agrees
with the particle-buffer count and the sum of per-species linked-list counts.

The run-level receipt SHA-256 is
`4f9570d54c460a4ad25b0b5290e03ae79c9a192fb31e36e31c218d133adb53eb`;
the native log SHA-256 is
`690c4aad886fe06eecb1bfd90a5a941bcf2346751e43a57cc824fcb6d718dc3d`.

## Implemented forward-model behavior

### Source handoff

- `tenishev2014_hybrid_8` is the only runtime source distribution currently
  implemented. The schema is extensible through the source-case identifier.
- The canonical source set has eight point jets and seven tiger-stripe
  polylines. Every published stripe direction is projected radially to the
  252.100 km Enceladus sphere. All 102 generated stripe nodes pass a `1e-6 m`
  compiled-surface check.
- The E3 and E5 startup receipts distribute each species total equally over
  the 15 sources only to make the forward model runnable. This is not a claim
  that the physical sources have equal activity. The standalone fit will
  replace these fixtures with a shape receipt or individually fitted rates.
- E3 totals are H2O `2.060909950673409e28` and apparent CO2 `5.356e25`
  molecules/s. E5 totals are H2O `4.520015885781454e28` and apparent CO2
  `1.013e26` molecules/s.
- `SOURCE_DISTRIBUTION_FILE.md` is the normative standalone/forward handoff
  syntax for both absolute component rates and normalized shape files.

### Runtime parser and physics adapter

- `input/enceladus.input` remains the infrequently changed build-generator
  deck. `input/enceladus.in` is parsed for every run. Unknown, duplicate,
  inactive, missing, or inconsistent controls fail before the simulation.
- `--validate-input` validates the complete per-run/source pair without
  initializing MPI or allocating the AMR mesh.
- Runtime controls select simulated/injected/chemistry species, reactions,
  solar activity, plasma authority, particle counts, speeds, Courant numbers,
  sampling and total lengths, seed, and generic split/merge mode.
- `plasma_background_mode=none` is a typed mode. The adapter zeroes every
  allocated DATAFILE background variable because the compiled Boris mover
  still interpolates E/B even when the active particles are neutral.
- The production injector creates the exact requested global count on each
  step. A fixed species weight is `Q*dt/N`; event time is uniform within the
  step; deterministic event replay and native AMR ownership choose one owner.
- Source velocities use the production-normalized Tenishev flux law. Gravity,
  the internal spherical boundary, and outer-domain deletion remain active.
- H2O loss uses the exact per-step constant-rate probability
  `1-exp[-dt*sum(nu_r)]`. No daughter, charge, recoil, or plasma feedback is
  created. CO2 has no invented Burger-2007 reaction and is transport-only.
- The final receipt counts injection, photolysis, surface absorption, true
  outer-domain escape, and resident particles separately. A thread-local
  synchronous marker distinguishes internal-sphere deletion from Boris's
  shared `_PARTICLE_LEFT_THE_DOMAIN_` return; per-species conservation is
  then checked against independently traversed native particle lists.

## Build and validation evidence

The final rebuild followed the required repository workflow. `.amps.conf` was
not edited directly:

```text
cd /home/vtenishe/Mars2/AMPS
# confirmed this root and no active build/test process
rm -rf -- build
./Config.pl -application=enceladus
./ampsConfig.pl -input enceladus.input -no-compile
env MAKEFLAGS="-j16" make amps
```

Build result: **PASS**. Only inherited warnings in legacy dust, legacy
Enceladus, and legacy filename-formatting code appeared; `ForwardModel.cpp`
and `ForwardModelAMPS.cpp` compiled without a reported warning or error.

- executable: `/home/vtenishe/Mars2/AMPS/amps`;
- size: 31,656,784 bytes;
- timestamp: 2026-11-01 01:42:30.966129988 -0600;
- SHA-256:
  `1f3110e4447d77742f7563c683da17f67034cc6f8396308bb2f265e972499a02`;
- run deck SHA-256:
  `0b0e7c61579370ce6a081ef6cb5c49255b536b6d56bccee8d5be23d1dd84caca`;
- `.amps.conf` SHA-256:
  `affef8bbc6d754159a42a4298faa26c43ae8de78c23a5cd89511c6d5ca0b029d`.

Application-independent suite:

```text
python3 test/run_tests.py
```

Result: **PASS, 41 checks, 0 failures**. It includes independent quadrature,
source normalization, parser/reference, reaction constants, timestep/weight,
projected-stripe, flyby-refinement configuration references, and 17 typed
negative checks. Receipt SHA-256:
`6c756838f7e607964d960d47d46678b2ae9f41a149e2a004d6327792eefbaa1b`.
The adjacent `reference_comparison.svg` visualizes calculated/reference error.

Native startup suite:

```text
python3 test/run_native_startup.py \
  --executable /home/vtenishe/Mars2/AMPS/amps \
  --output test_output/srcEnceladus/native_startup_refined_active_output
```

Result: **PASS, 6 checks, 0 failures**. Receipt SHA-256:
`84434d49d142258017a55df77e8596892e45d974219911230044dc3704163bcd`.

### Retained superseded attempts

Two first attempts are retained under `e3_first_run/` and
`e3_first_run_singleton/`. Both failed before `MPI_Init` because the restricted
execution sandbox exposed no network interface for Open MPI. They are
environmental failures, not model runs. The successful command was therefore
run with host-level local-MPI socket access.

An earlier successful native run is retained under `e3_first_run_native/`.
Its dynamics and total particle count completed, but its diagnostic wrapper
classified sphere absorptions both as surface losses and domain escapes
because Boris shares one return code for those deletions. That receipt is
superseded by `e3_first_run_qualified/`, whose separated sink ledger closes
exactly.

## Limits and next validation gates

- The current result is one seed, one rank, one refined mesh, one timestep and
  one particle count over ten 1,000-step windows. It checks execution, repeated
  output, requested trajectory resolution and exact discrete bookkeeping; it
  is not a mesh/timestep/particle convergence study.
- Fresh particles receive chemistry beginning with the following full step;
  their initial partial-birth-step exposure is omitted. This operator-split
  error is `O(nu*dt)` and requires a timestep study or an exact birth-to-end
  loss treatment before chemistry accuracy is qualified.
- Only H2O parent decay is modeled. No daughters, CO2 photolysis, charged
  species, charge exchange, H3O+, plasma, or back-coupling is active.
- The apparent CO2 rate and equal 15-way fixture are not final physical source
  products. A fitting-model receipt with justified CO2 interpretation remains
  required.
- Four-rank ownership/rank-equivalence, multiple seeds, expanded domain,
  mesh/timestep/particle convergence, an explicit spatial observation
  operator, analytic profile comparison, and observation-operator comparison
  remain open. The 15.75625 km cell-centered-constant flyby sampler resolves
  the old cut-cell plateau but remains too sparsely populated at 100 injected
  particles/species/iteration to interpret exact closest-approach zeros.
- Generic refined volume fields are now limited to runtime-active H2O and CO2.
  Small surface files for compiled inactive species remain generic diagnostics;
  they do not mean those species were simulated or produced.
- `PT/amps_error.log` contains the generic pre-allocation warning that no mesh
  blocks yet existed and that connectivity-list mode was adjusted. The same
  run subsequently allocated and advanced 3,165 blocks and returned process
  code zero; the warning is retained rather than suppressed.

No observational agreement, production readiness, or unique source recovery
is claimed by the current E3 run.


## Case runner validation commands (2026-11-13 final receipt)

All commands below run from `/home/vtenishe/Mars2/AMPS`; native rebuilds use
only the required root `build` cleanup and established configuration workflow.
The machine overlay changes output/work locations without changing physics.

```bash
env ENCELADUS_VALIDATION_DATA_ROOT=/data/vtenishe/cassini_validation_data \
  python3 -B -m unittest discover -s srcEnceladus/standalone/tests -v
python3 -B src/test/run_tests.py \
  --config test_output/enceladus-case-runtime-20261113/machine.json \
  --test CORE-HARNESS-001 --test CORE-GENERAL-INPUT-001 \
  --test APPLICATION-ENCELADUS-FORWARD-001 \
  --test APPLICATION-ENCELADUS-NATIVE-STARTUP-001 --jobs 1 --streams 1
python3 -B test_output/enceladus-case-runtime-20261113/production-build/build.py
python3 -B srcEnceladus/test/run_native_startup.py --executable ./amps \
  --output test_output/enceladus-case-runtime-20261113/production-startup
./amps -i input/enceladus_cases.in --list-cases
./amps -i input/enceladus_cases.in -case E3--tenishev2014_hybrid_8 --describe-case
mpiexec -n 4 ./amps -i input/enceladus_cases.in \
  -case E3--tenishev2014_hybrid_8 --validate-runtime-input
make -C src/general/input test
python3 -B src/test/generate_test_catalog.py --check
python3 -B doc/update_test_catalog.py --check
```

The standalone/publisher/native test output directories are retained and
creation-only. Choose NEW output directories when repeating the scripts;
do not overwrite qualified evidence. The production build driver likewise
uses creation-only logs; its recorded command sequence is the reproducible
build recipe, not a command to overwrite the existing receipt. MPI runs on
this host require normal host local-socket access; sandbox-only socket errors
are an execution-environment failure, not a model result.
