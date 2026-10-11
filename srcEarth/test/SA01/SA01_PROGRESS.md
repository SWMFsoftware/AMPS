# SA01 progress and handoff record

## Current state

- Status 2026-10-11T19:29Z: SA01 now has two independent repository entry points.
  `srcEarth/test/list` contains the full QUICK campaign as an exclusive `P` test;
  `srcEarth/studies/list` contains the full STANDARD campaign as a study. Each entry
  executes both spacecraft, evaluates the unchanged comparison gates, and requests
  all comparison figures in PNG and EPS. Each uses a separate output parent and the
  runner creates a new UTC/PID campaign directory, so repeated or concurrent QUICK
  and STANDARD invocations cannot reuse or overwrite one another's evidence.
- The list commands use the explicit `WORKFLOW` exit-code policy. This policy returns
  success only after physical execution, evaluated comparison, report generation,
  and verified plot-manifest creation all complete. It does not convert a scientific
  FAIL into a scientific PASS: `comparison_summary.json` and
  `qualification_summary.json` retain the actual unchanged result, and the new
  `workflow_completion.json` records that result plus hashes of the evidence and
  figures. Direct CLI use retains the stricter `QUALIFICATION` default.
- The QUICK list entry is preceded by `! runner: exclusive`. The repository runner
  therefore waits for its process pool to drain and admits no other test until QUICK
  exits. The command itself has no list-only CLI prerequisite and is independently
  executable exactly as written; standalone resource isolation remains the
  operator's responsibility.
- Verification at current source: 38/38 SA01 unit/contract tests pass; the complete
  Step-12 release suite passes all eight validation and five source-contract tests;
  the normal test runner dry-run parses 62 tests and recognizes QUICK as the second
  exclusive test. Independent real-data
  preparation passes for QUICK (59 RBSPA / 66 RBSPB epochs) and STANDARD
  (251 / 261). Prior completed QUICK and STANDARD plot manifests each identify 12
  PNG and 12 EPS files, all present. No new physical SA01 campaign was launched.
- Status 2026-10-11T19:04Z: the earlier QUICK and STANDARD physical campaigns are
  complete. QUICK completed 59/59 RBSPA and 66/66 RBSPB snapshots; STANDARD
  completed 251/251 and 261/261. Both native branches exited 0. Both campaign
  qualification summaries remain `FAIL`: 10/34 observational gates pass, 24/34
  fail, and the unchanged 2% common-population convergence gate fails. These are
  retained scientific results, not infrastructure failures and not candidates for
  gate relaxation.
- The current change adds an explicit two-level snapshot scheduling policy without
  altering SA01 physics, sampling, trajectory tolerances, or acceptance gates.
  `COLLECTIVE` is the safe default and keeps one distributed AMR mesh while all MPI
  ranks process epochs sequentially. `EPOCH_PARALLEL` is an explicit opt-in for
  independent epoch teams; it retains the private-domain requirement and now refuses
  startup when the measured first-domain resident-memory increment projects that the
  remaining replicas would consume more than the configured fraction of the node's
  available memory. SA01 preparation records `EPOCH_PARALLEL` and the default 0.80
  fraction explicitly so its intended cross-epoch parallelism is never inferred.
- A clean build at current HEAD `af26f41725eb42f15f4d895ee8c916234f177425`
  passed. Executable SHA-256 is
  `31ebdbb25327a81a95e324e183fcb0b6180ceeda1af94adb5dfde1d73fabb769`;
  build log: `test_output/C19_snapshot_policy_20261011/build.log`.
- C19 was run through the repository test runner after the build using explicit
  `COLLECTIVE`. It passed the execution/regression contract in 4744.6 s with zero
  mismatches and zero timeouts. The run produced 12/12 model rows, completed all
  three field epochs with status 0, passed trajectory-resolution and plot-consistency
  checks, and generated 27 PNG plus 27 EPS files. C19's provisional observational
  thresholds were not enforced by the SMOKE profile and its scientific result remains
  `observational_passed=false`; that distinction is preserved in `C19_result.json`.
- No physical SA01 run was started after the scheduling/list-integration changes.
  SA01's current 38 implementation/contract tests were executed; all passed.
- Standing instructions: `AGENTS.md`.
- Execution specification: `SA01_TASK.md`.
- Detailed initial handoff: `/data/vtenishe/GEO/AMPS/handoff.md` (historical state).
- Target repository: `/data/vtenishe/GEO/AMPS`.
- Target data package: `/data/vtenishe/GEO/sep-in-geospace-sa01-sa04-standalone-20260926`.
- Source dirty state includes this policy implementation and documentation. The
  pre-existing unrelated `Makefile` modification and generated/untracked files remain
  preserved and were not incorporated into this work.
- C19 evidence root:
  `test_output/C19_snapshot_policy_20261011/C19_goes_epead_ew`; runner report:
  `test_output/C19_snapshot_policy_20261011/test_runner_report.json`.
- Next action: execute QUICK through `srcEarth/test/test_runner.py` when a full QUICK
  validation run is intended, or execute the STANDARD command in
  `srcEarth/studies/list` for the study. Inspect each new campaign's
  `workflow_completion.json`, `qualification_summary.json`, and `figures/`; do not
  report workflow completion as scientific qualification.

## Evidence carried from this conversation

The supplied archive was previously inspected in a separate environment. Its
four-case data-integrity validator and ten pipeline reference tests passed there.
This is historical package/interface evidence, not evidence that current target-host
AMPS runs, scientific comparisons, or convergence have passed. Verify local files
and revision before relying on prior results.

The archive contains preprocessed SA01 observations and prepared RBSPA/RBSPB inputs.
Its package status states no physical AMPS runs were included. No physical SA01
execution is claimed by the preparation of these instructions.

## Stage tracking

| Stage | State | Evidence / notes |
| --- | --- | --- |
| Target repository and instructions inspection | Complete | Revision, dirty state, executable/build, process state and available reference implementation inspected 2026-10-08 |
| Local SA01 data validation | Complete | Package auditor PASS: 14 raw files, 13 processed hashes, 983 complete boundary spectra; report under `test_output/SA01/inspection` |
| Runner implementation | Complete | SA01-specific prepare/run/compare/plot/common modules, template, README and requirements created |
| Numerical/interface verification | Complete for implementation | 38/38 current SA01 tests pass; prior physical repeatability controls pass exactly; snapshot and list-integration policies are covered by parser/CLI and source-contract tests |
| Build / executable verification | Complete | Mandatory clean `make -j` passed; current MPI+OpenMP executable SHA-256 `31ebdbb...`; retained C19-policy build log |
| QUICK physical run | Complete, unqualified baseline | Corrected HC4 completed RBSPA 59/59 and RBSPB 66/66; native runs PASS, campaign qualification FAIL |
| STANDARD RBSPA run | Complete, unqualified baseline | Corrected HC4 completed 251/251 snapshots; native run PASS |
| STANDARD RBSPB run | Complete, unqualified baseline | Corrected HC4 completed 261/261 snapshots; native run PASS |
| Observational comparisons | Complete, FAIL | QUICK and STANDARD each pass 10/34 unchanged gates and fail 24/34 |
| Visualization | Complete for prior campaigns | Reports and PNG/EPS artifacts exist under both completed campaign roots |
| Numerical convergence | FAIL | Fixed-field timestep and full common-population refinement fail the unchanged 2% criterion |
| Snapshot-policy validation | C19 PASS; SA01 prepared | C19 collective run passed and generated 54 plots; both SA01 profiles independently pass real-data preparation under explicit `EPOCH_PARALLEL`; no new physical SA01 run was started |
| Repository list integration | Complete | QUICK is an exclusive master-list test; STANDARD is an independent study; both require full execution, comparison, reporting, and PNG/EPS plotting before workflow completion |
| Final report and documentation | Complete for this change | Runtime and workflow/qualification policies are documented in test and SA01 READMEs; this journal records exact evidence and scientific limitations |

## Session history

### Instruction package prepared — 2026-10-08 UTC

- Created standing instructions, detailed task specification and this initial journal.
- Preserved user-specified paths, rebuild procedure, scientific scope and fixed gates.
- Did not access or change the user's target-host AMPS repository.
- Did not compile or execute physical AMPS in preparing this instruction package.
- Next action: begin target-host inspection and data validation.

### 2026-10-08T05:11:43Z — target-host repository inspection

- Goal: establish the live repository, build, executable, process, and reference-code state before any mutation or run.
- Repository revision and dirty state: HEAD `6104004ca26d77881ad9d05c6a6a21dd9a9577cc`; modified `Makefile` changes `EXTRALINKEROPTIONS_F` to `-lmpi++`; untracked `SWMF_data`, `amps`, `gitinfo.txt`, and `srcEarth/test/SA01/` were already present.
- Files changed (and why): this progress journal only, to record the required handoff evidence.
- Commands executed, working directories, and return codes: read all three SA01 instruction/task/progress files; `git status`, `git diff`, revision inspection, process listing, executable/build inspection, and data/reference-tree discovery from repository root; all substantive inspection commands returned 0 (one bounded `sort | head` inventory emitted benign SIGPIPE diagnostics).
- Build/executable identity: `amps` SHA-256 `a59dbb15b6da3ee64d3400d750392aaa929fdaddd3a7d2dfc535af0e209d7eee`; ELF x86-64, not stripped; source/configuration match not yet proven.
- Dataset/input/output identities: required processed SA01 files and packaged RBSPA/RBSPB baseline decks are present under `/data/vtenishe/GEO/sep-in-geospace-sa01-sa04-standalone-20260926`; no run output is claimed.
- Result and status: PASS for initial inspection only; no scientific or execution gate evaluated.
- Evidence paths: repository state is summarized above; package reference implementation is under the extracted package's `srcEarth/studies/rbsp_sep_validation` and `srcEarth/test/URBSPSEPValidation` trees, while those trees are absent from the live repository.
- Scientific/numerical decisions and justification: none yet.
- Remaining issue or blocker: validate package hashes/schema and determine what reference routines can be reviewed/reused without copying packaged source over the live repository.
- Active process state: no matching `amps`, `mpiexec`, `mpirun`, or `make` process observed.
- Exact next action / resume command: inspect packaged provenance/manifests and execute its read-only SA01 validator against the extracted data root.

### 2026-10-08T05:11:43Z — immutable data and reference-pipeline validation

- Goal: verify the local SA01 data package and establish which historical routines are scientifically safe to reuse.
- Repository revision and dirty state: unchanged from the preceding entry.
- Files changed (and why): this progress journal only.
- Commands executed, working directories, and return codes: package `validate_data_package.py --cases SA01` returned 0; packaged `URBSPSEPValidation/test_pipeline.py` ran 10 tests and returned 0; dependency/build/parser source inspection returned 0.
- Build/executable identity: existing executable remains SHA-256 `a59dbb15b6da3ee64d3400d750392aaa929fdaddd3a7d2dfc535af0e209d7eee`; `gitinfo.txt` identifies AMPS `6104004`; `Makefile.local` configures `/home/vtenishe/SPICE/cspice` and kernels, and current source supports `OUTPUT_MODE TRAJECTORY` in gridless density/spectrum mode.
- Dataset/input/output identities: SA01 auditor checked 14 raw manifest files, all 14 provenance raw hashes (including the manifest), 13 processed hashes, 983 positive ten-channel boundary rows, 787 RBSPA and 814 RBSPB eligible samples, and six original-data figure signatures.
- Result and status: PASS for immutable package integrity and historical synthetic/reference tests. This is not a physical-run or observation-agreement PASS.
- Metrics and comparison counts: 10/10 historical pipeline tests passed. No AMPS comparison rows exist yet.
- Evidence paths: `test_output/SA01/inspection/sa01_data_validation.json`; immutable package provenance at `validation_data/processed/SA01/provenance.json` under the extracted package.
- Scientific/numerical decisions and justification: retain the package's exact-boundary-row trajectory policy and rectangular REPT folding. Do not reuse its `bin_average_loglinear` unchanged because it filters zero model nodes and can interpolate across physical-zero support, contrary to the task contract; the SA01 implementation will preserve zeros and use explicit linear transition handling.
- Remaining issue or blocker: implement and validate the live SA01 tools before launching expensive work.
- Active process state: no physical run or build active.
- Exact next action / resume command: create `sa01_common.py`, the four entry points, tests, and README under `srcEarth/test/SA01`.

### 2026-10-08T05:11:43Z — SA01 runner implementation and unit verification

- Goal: implement the complete SA01-specific workflow without mutating observations or importing the older packaged tree into the live source tree.
- Repository revision and dirty state: HEAD remains `6104004ca26d77881ad9d05c6a6a21dd9a9577cc`; pre-existing dirty state preserved; new files are confined to `srcEarth/test/SA01` and generated inspection evidence under `test_output/SA01`.
- Files changed (and why): added `AMPS_PARAM.in.template`, `sa01_common.py`, `prepare_SA01.py`, `compare_SA01.py`, `plot_SA01.py`, `run_SA01.py`, `README.md`, `requirements.txt`, and `tests/test_sa01.py`; updated this journal.
- Commands executed, working directories, and return codes: Python byte compilation returned 0; CLI help checks returned 0; initial 10-test run returned 1 because this host's `Path.write_text` lacks a `newline` parameter; the compatibility defect was fixed and the repeated 10-test run returned 0; `git diff --check` returned 0.
- Build/executable identity: no rebuild and no physical run in this stage.
- Dataset/input/output identities: immutable data validator remains PASS; no prepared physical inputs yet.
- Result and status: PASS for implementation unit/failure tests; physical execution and observational status remain unevaluated.
- Metrics and comparison counts: 10/10 SA01 tests pass, covering units/relativistic density, analytic integration, zero preservation, nonfinite/support failure, duplicate/malformed Tecplot input, exact boundary availability, absolute rendering, gate denominator, near-zero convergence, and rebuild-target safety.
- Evidence paths: `srcEarth/test/SA01/tests/test_sa01.py`; prior data report `test_output/SA01/inspection/sa01_data_validation.json`.
- Scientific/numerical decisions and justification: positive-positive model segments use exact log-linear integration; a segment touching physical zero uses an explicit local linear transition, never node deletion/bridging. Density reconstruction integrates `J/v_rel` continuously, while the REPT reference retains its documented center-speed bin approximation.
- Remaining issue or blocker: exercise the rendered deck against the live executable and fix any parser/runtime-interface mismatch.
- Active process state: no build or physical run active.
- Exact next action / resume command: `python3 srcEarth/test/SA01/run_SA01.py --prepare-only --data-root /data/vtenishe/GEO/sep-in-geospace-sa01-sa04-standalone-20260926 --output-dir test_output/SA01/quick_prepare_20261008 --profile QUICK --mpi-ranks 1 --threads 1`.

### 2026-10-08T05:20:00Z — real-data preparation and live parser smoke

- Goal: validate generated QUICK selections/decks against immutable SA01 data and the live executable before a full physical run.
- Repository revision and dirty state: unchanged HEAD and pre-existing dirty state.
- Files changed (and why): `prepare_SA01.py` replaced Python-3.9-only `str.removesuffix` after this host's Python 3.8.10 exposed the compatibility failure; journal updated.
- Commands executed, working directories, and return codes: first real-data prepare returned 1 at `str.removesuffix` and left a clearly incomplete preserved directory `test_output/SA01/quick_prepare_20261008`; after the fix, tests returned 0 and preparation into `_v2` returned 0. An in-sandbox MPI smoke returned 1 because Open MPI could not open a local network interface; the approved unsandboxed repeat parsed/initialized successfully and was intentionally stopped by `timeout 20s` (exit 124) after entering the first cutoff task.
- Build/executable identity: no rebuild; live executable SHA-256 remains `a59dbb15b6da3ee64d3400d750392aaa929fdaddd3a7d2dfc535af0e209d7eee`.
- Dataset/input/output identities: `quick_prepare_20261008_v2` contains absolute decks; RBSPA kept 59 epochs (726 cadence exclusions, 2 missing-boundary exclusions), RBSPB kept 66 (745 cadence exclusions, 3 missing-boundary exclusions).
- Result and status: PASS for preparation and parser/runtime startup. The smoke timeout is diagnostic only and not a completed physical run.
- Metrics and comparison counts: no observation comparisons yet.
- Evidence paths: `test_output/SA01/quick_prepare_20261008_v2/{preparation_summary.json,RBSPA/prepared_run.json,RBSPB/prepared_run.json}`; isolated smoke manifest/log-side products under `test_output/SA01/parser_smoke_20261008`.
- Scientific/numerical decisions and justification: the live parser accepted T05, zero electric field, trajectory mode, exact boundary table, strict gap policy, seven channels, and requested QUICK controls. The driver reader warned that the table header lacks a column-map annotation but loaded all 1,044 positional records and proceeded through validated standalone initialization; no data were altered.
- Remaining issue or blocker: full QUICK cost/runtime and possible unresolved-trajectory failure remain unknown.
- Active process state: the time-bounded smoke process exited; no detached process remains.
- Exact next action / resume command: run the unified CLI with `--profile QUICK --mpi-ranks 4 --threads 16` into `test_output/SA01/quick_physical_20261008`.

### 2026-10-08T05:37:00Z — QUICK physical interface defect and correction

- Goal: execute physical QUICK RBSPA/RBSPB and proceed automatically to comparison/plots.
- Repository revision and dirty state: HEAD unchanged; no rebuild; new SA01 code/journal plus generated evidence only.
- Files changed (and why): corrected the density template/preparer so `DS_TRANSMISSION_SCAN_N=DS_NINTERVALS+1` and `DS_MAX_PARTICLES=energy_nodes*angular_samples`; enhanced the launcher to record and terminate cleanly on interruption; updated tests and this journal.
- Commands executed, working directories, and return codes: full unified QUICK launch began in `test_output/SA01/quick_physical_20261008`; RBSPA cutoff completed all 59 locations with four-rank accounting exactly 59/59. Density startup then printed the decisive effective configuration. The run was intentionally interrupted with SIGINT before invalid science work continued; the Python wrapper exited by `KeyboardInterrupt` before the subsequent interruption-handling patch.
- Build/executable identity: unchanged executable SHA-256 `a59dbb15b6da3ee64d3400d750392aaa929fdaddd3a7d2dfc535af0e209d7eee`.
- Dataset/input/output identities: invalid attempt used immutable SA01 inputs and produced a complete diagnostic cutoff file, but no complete density/spectrum/flux products. None of its partial outputs may be compared.
- Result and status: ERROR (diagnostic interface failure), not a scientific FAIL. Current-source semantics are `DS_TRANSMISSION_SCAN_N = energy-node count` and `DS_MAX_PARTICLES = total trajectories per point`; the historical packaged template incorrectly assigned the profile's 96 angular samples to both.
- Metrics and comparison counts: cutoff task accounting 59/59; density revealed 96 energy nodes, one direction per energy, and 4/5664 tasks before interruption; zero observation comparisons.
- Evidence paths: `test_output/SA01/quick_physical_20261008/RBSPA/amps.log`, `cutoff_gridless_points.dat`, and `interrupted_run.json`.
- Scientific/numerical decisions and justification: profile columns are now interpreted literally: QUICK has 24 energy intervals = 25 nodes and 96 angular directions at every energy; its cap is therefore 2,400 trajectories per point. A one-direction result cannot represent the specified isotropic angular integration and was rejected rather than mislabeled QUICK.
- Remaining issue or blocker: correct QUICK is intrinsically expensive under multi-epoch T05 because live AMPS safely falls back from requested 16 threads to one intra-rank worker while preserving MPI parallelism. Runtime must be measured and recorded; difficulty alone is not treated as a blocker.
- Active process state: SIGINT terminated the Python launcher and MPI children; no physical run remains active.
- Exact next action / resume command: prepare `test_output/SA01/quick_corrected_20261008`, inspect its effective controls, and start the corrected physical run.

### 2026-10-08T05:49:20Z — corrected QUICK launched and plotting compatibility checked

- Goal: launch the corrected literal QUICK workload and finish independent implementation verification while it runs.
- Repository revision and dirty state: HEAD unchanged; no rebuild; SA01-only changes plus generated evidence.
- Files changed (and why): plotting now uses the actual `soho_erne_5min.csv` flux column and converts timezone-aware pandas timestamps to timezone-naive UTC arrays for compatibility with installed Matplotlib 3.1/Pandas 2.0; added this active-run record/journal entry.
- Commands executed, working directories, and return codes: corrected prepare returned 0 and rendered `DS_NINTERVALS=24`, `DS_TRANSMISSION_SCAN_N=25`, `DS_MAX_PARTICLES=2400` in both branches; 10/10 tests returned 0. A plot-context smoke first exposed two real compatibility errors (wrong ERNE column name, then old-Matplotlib timestamp handling); both were fixed and the repeat generated a valid 1857x1358 PNG. A `nohup` attempt returned a PID but was immediately reaped by the execution supervisor before Python/log output; it is not active. The attached corrected run then launched successfully.
- Build/executable identity: current-revision executable SHA-256 unchanged.
- Dataset/input/output identities: corrected root `test_output/SA01/quick_corrected_20261008`; RBSPA/RBSPB prepared hashes are in their `prepared_run.json` files; active RBSPA log is `RBSPA/amps.log`.
- Result and status: RUNNING. At journal update, corrected RBSPA cutoff is 12/59 (20.3%); no error observed.
- Metrics and comparison counts: none yet. Plot smoke is interface evidence only, not an observational figure or PASS.
- Evidence paths: `test_output/SA01/quick_corrected_20261008/active_run.json`, `RBSPA/amps.log`, and prepared records; temporary plot smoke `/tmp/sa01_plot_smoke/sa01_orbit_upstream_storm_context.png`.
- Scientific/numerical decisions and justification: literal profile mapping retained. Requested 16 threads are expected to fall back to one worker per rank for changing T05 snapshots; four-rank MPI scheduling remains effective and this requested/effective distinction must be reported.
- Remaining issue or blocker: long physical runtime; no scientific status can be assigned until complete products exist.
- Active process state: managed exec session 15967; PIDs and exact command are recorded in `active_run.json`. Do not rebuild or launch a competing SA01 run.
- Exact next action / resume command: poll managed session 15967 (or read `RBSPA/amps.log` and verify PIDs 269077/269144/269149--269152), then append the next meaningful execution stage.

### 2026-10-08T06:14:15Z — corrected RBSPA cutoff complete; density workload verified

- Goal: verify the corrected profile from live solver output and hand off the long active calculation without losing process identity.
- Repository revision and dirty state: unchanged; no rebuild during the active run.
- Files changed (and why): `active_run.json` and this append-only journal were updated with the live stage/resource evidence.
- Commands executed, working directories, and return codes: continued polling managed exec session 15967; read-only process/log snapshot returned 0.
- Build/executable identity: unchanged current-revision executable.
- Dataset/input/output identities: active output root `test_output/SA01/quick_corrected_20261008`; complete corrected RBSPA cutoff diagnostic is `RBSPA/cutoff_gridless_points.dat`; primary density/flux products do not exist yet.
- Result and status: corrected RBSPA cutoff completed 59/59 with rank task sum 59. Density/spectrum is RUNNING.
- Metrics and comparison counts: live solver reports `Npoints=25`, `Directions grid=96`, `DS_MAX_PARTICLES=2400`, exactly matching QUICK's 24 intervals and 96 angular samples. No observation comparison exists.
- Evidence paths: `RBSPA/amps.log`, `RBSPA/cutoff_gridless_points.dat`, and root `active_run.json` under the active output directory.
- Scientific/numerical decisions and justification: retain the exact corrected workload despite cost. Requested 16 threads safely fell back to one per MPI rank for multi-epoch T05; four ranks remain continuously utilized, so recorded effective layout is 4 MPI x 1 intra-rank worker.
- Remaining issue or blocker: calculation is long; no completed density block/progress estimate was emitted by 2026-10-08T06:14:15Z. This is active computation, not a blocker or failure.
- Active process state: runner PID 269077, MPI PID 269144, AMPS PIDs 269149--269152; all four AMPS processes in `R` state at approximately 99.8% CPU; managed session 15967.
- Exact next action / resume command: poll session 15967 or tail `test_output/SA01/quick_corrected_20261008/RBSPA/amps.log`; do not rebuild or launch another physical run while these PIDs are active.

### 2026-10-08T06:25:46Z — active-run status check

- Goal: answer whether the recorded physical run remains active and capture its first density progress estimate.
- Commands executed, working directories, and return codes: read-only process, artifact, run-status, and log inspection returned 0.
- Result and status: RUNNING. Runner and MPI parent remain alive; all four AMPS ranks are in `R` state at 99.8--99.9% CPU.
- Metrics and comparison counts: RBSPA density/spectrum is 12/4425 MPI tasks (0.3%); solver ETA is 107:47:06. No complete native product or observational comparison exists yet.
- Evidence paths: active output `test_output/SA01/quick_corrected_20261008`, especially `active_run.json` and `RBSPA/amps.log`.
- Active process state: managed session 15967; runner 269077, MPI 269144, ranks 269149--269152.
- Exact next action / resume command: continue polling session 15967 or the RBSPA log; do not start another run or rebuild.

### 2026-10-08T06:33:37Z — multi-epoch shared-memory parallelism diagnosis

- Goal: determine why requested intra-rank threads are disabled for a trajectory file with many epochs and identify the changes required to use the reentrant T05/GEOPACK implementation safely.
- Repository revision and dirty state: unchanged; no source, executable, input, or observation data were modified during this diagnostic pass.
- Files changed (and why): this progress journal only, to preserve the architectural finding and active-run handoff.
- Commands executed, working directories, and return codes: read-only inspection of the gridless cutoff/density schedulers, T05 wrapper/Fortran source, GEOPACK wrapper/Fortran source, and active process/log state from the repository root returned 0.
- Build/executable identity: no rebuild; the active run continues to use executable SHA-256 `a59dbb15b6da3ee64d3400d750392aaa929fdaddd3a7d2dfc535af0e209d7eee`.
- Result and status: DIAGNOSTIC ONLY. `T05.for` is explicitly reentrant and `geopack2008.for` already provides immutable multi-epoch IGRF context routines, but the gridless C++ evaluator still routes through namespace-global `T05::PS`/`T05::PARMOD` and legacy common-block `Geopack::IGRF::GetMagneticField`. Cutoff therefore imposes the `oneFrozenFieldSnapshot` fallback. The multi-rank density path separately processes each MPI chunk serially and calls `ResolveMpiDynamicChunk(...,1,...)`, so removing the cutoff guard alone would not thread density.
- Metrics and comparison counts: active RBSPA density reached 20/4425 tasks (0.5%, ETA 96:35:29); no observational comparison exists.
- Evidence paths: `srcEarth/gridless/CutoffRigidityGridless.cpp`, `srcEarth/gridless/DensityGridless.cpp`, `src/interface/T05Interface.cpp`, `src/interface/T05.for`, `src/interface/GeopackInterface.cpp`, `src/interface/geopack2008.for`, and `test_output/SA01/quick_corrected_20261008/RBSPA/amps.log`.
- Scientific/numerical decisions and justification: enable multi-epoch threading by preparing one immutable field context per unique epoch serially on the rank/main thread, including epoch-specific IGRF coefficients/rotation, dipole tilt, and T05 `PARMOD`; worker field calls must pass that context and parameters explicitly to reentrant T05/IGRF entry points. Do not merely remove the guard or call legacy `Geopack::Init` from workers. Add a density intra-rank batch worker pool with main-thread-only MPI and serial post-join accumulation (or worker-private reductions) to avoid same `(point,energy)` accumulator races.
- Remaining issue or blocker: implementation and rebuild must wait until the active physical run completes or the user explicitly chooses to stop it; the mandatory rebuild procedure forbids rebuilding while it is active.
- Active process state: host process inspection confirms runner 269077, MPI parent 269144, and AMPS ranks 269149--269152 remain active in managed session 15967.
- Exact next action / resume command: keep polling session 15967; after it ends, implement/expose the GEOPACK context API, refactor gridless field snapshots, add the density worker pool, then rebuild and validate one-thread versus many-thread multi-epoch equivalence.

### 2026-10-08T15:54:12Z — gridded density-path capability diagnosis

- Goal: determine whether density/flux can use the same gridded background magnetic-field representation as cutoff rigidity.
- Repository revision and dirty state: unchanged; no production source, executable, input, or observation data were modified.
- Files changed (and why): this progress journal only.
- Commands executed, working directories, and return codes: read-only tracing of the gridless and Mode3D dispatch, density trajectory kernels, Mode3D field preparation, temporal snapshot batching, documentation, active log, and host process state returned 0.
- Build/executable identity: no rebuild; active executable unchanged.
- Result and status: DIAGNOSTIC ONLY. Mode3D density/flux is mesh capable and calls `TraceTrajectoryMesh`, using the same prepared compact read-only AMR field snapshot as Mode3D cutoff. Gridless density always calls the shared gridless/direct evaluator and has no switch to consume the Mode3D mesh. The active SA01 command is `-mode gridless`, so both its cutoff and density stages are direct-field calculations.
- Metrics and comparison counts: active RBSPA density reached 316/4425 tasks (7.1%); no observational comparison exists.
- Evidence paths: `srcEarth/main.cpp`, `srcEarth/3d/Mode3D.cpp`, `srcEarth/3d/DensityMode3D.cpp`, `srcEarth/3d/CutoffRigidityMode3D.cpp`, `srcEarth/gridless/DensityGridless.cpp`, and `test_output/SA01/quick_corrected_20261008/RBSPA/amps.log`.
- Scientific/numerical decisions and justification: a gridded SA01 trajectory run should use `-mode 3d -mode3d-field-eval INTERPOLATION` and explicit `TEMPORAL_MODE SNAPSHOT_LIST` batching. Mode3D prepares one epoch-specific mesh snapshot at a time, selects the matching trajectory rows, runs both requested products against that frozen snapshot with MPI and intra-rank threads, then advances to the next epoch. This is not simultaneous multi-epoch mesh residency and requires mesh-resolution convergence plus SA01 runner/output aggregation support before it can replace the existing gridless baseline.
- Remaining issue or blocker: no Mode3D SA01 run has been prepared or validated, and the active gridless run must not be disturbed or the repository rebuilt without explicit direction.
- Active process state: runner 269077, MPI parent 269144, and AMPS ranks 269149--269152 remain active in managed session 15967 after about 10 hours elapsed.
- Exact next action / resume command: continue monitoring the active gridless baseline; after it ends, choose between implementing reentrant multi-epoch gridless threading or adding a separately labeled Mode3D mesh-convergence branch to the SA01 runner.

### 2026-10-08T15:57:14Z — comparison and figure artifact check

- Goal: determine whether the SA01 execution has produced any valid observational comparisons or plots.
- Files changed (and why): `active_run.json` and this progress journal were refreshed with the latest recorded progress and artifact state; no source or scientific input changed.
- Commands executed, working directories, and return codes: recursive comparison/metrics/report/PNG/EPS inventory, complete active-output inventory, log tail, and run-record inspection returned 0.
- Result and status: no valid SA01 comparison table, metrics report, PNG, or EPS has been produced. The corrected run currently has only the complete RBSPA cutoff table plus incomplete density execution evidence. RBSPB has prepared inputs but has not started. The earlier `/tmp/sa01_plot_smoke/sa01_orbit_upstream_storm_context.png` was only a plotting-interface smoke check and is not a model/reference result.
- Metrics and comparison counts: RBSPA density 316/4425 tasks (7.1%); zero completed model/reference comparisons.
- Evidence paths: `test_output/SA01/quick_corrected_20261008/RBSPA/cutoff_gridless_points.dat`, `RBSPA/amps.log`, and root `active_run.json`; no comparison/figure files exist under `test_output/SA01`.
- Active process state: unchanged from the 2026-10-08T15:54:12Z host check; runner 269077, MPI parent 269144, and ranks 269149--269152 were active.
- Exact next action / resume command: allow RBSPA and RBSPB physical products to complete; only then run the saved-output comparison and plotting stages and report their actual gates/artifact paths.

### 2026-10-08T15:58:46Z — active-run progress check

- Goal: report current completion of the corrected QUICK physical run.
- Files changed (and why): refreshed `active_run.json` and this journal; no scientific source or input changed.
- Result and status: RUNNING after approximately 10 h 13 min. RBSPA cutoff remains complete at 59/59; RBSPA density/spectrum is 316/4425 tasks (7.1%), leaving 4109 RBSPA density tasks. RBSPB has not started.
- Metrics and comparison counts: four AMPS ranks are in `R` state at 99.9% CPU each. The most recent printed solver ETA is 92:27:26, but the ETA has varied substantially with trajectory cost and is not a reliable completion prediction. With 66 prepared RBSPB epochs and the same 25-energy/3-direction-block layout, 316 completed tasks represent approximately 3.4% of the inferred 9375 density tasks across both spacecraft; this combined percentage excludes cutoff overhead.
- Evidence paths: `test_output/SA01/quick_corrected_20261008/RBSPA/amps.log` and root `active_run.json`.
- Active process state: runner 269077, MPI parent 269144, AMPS ranks 269149--269152.
- Exact next action / resume command: continue monitoring without rebuilding or starting a competing run.

### 2026-10-08T16:02:38Z — slow direct-field run terminated for optimization

- Goal: stop the prohibitively slow direct-field baseline at the user's request, preserve evidence, and clear the repository for code improvement and rebuild.
- Files changed (and why): the runner wrote `RBSPA/run_status.json`; `active_run.json` and this journal were updated to mark the run interrupted. No partial science output was relabeled complete.
- Commands executed, working directories, and return codes: sent Ctrl-C to managed session 15967; the runner terminated its MPI child, returned the expected interruption error, and host process inspection confirmed all recorded SA01 PIDs exited.
- Result and status: INTERRUPTED BY REQUEST, not a scientific result. RBSPA cutoff completed 59/59; density stopped at 316/4425 tasks (7.1%); required spectrum/flux/density/termination products are absent; RBSPB never started.
- Evidence paths: `test_output/SA01/quick_corrected_20261008/RBSPA/{amps.log,run_status.json,cutoff_gridless_points.dat}` and root `active_run.json`.
- Scientific/numerical decisions and justification: preserve the partial run only as performance and interface evidence. It is excluded from comparisons because the density product transaction did not complete.
- Active process state: none; former runner 269077, MPI parent 269144, and ranks 269149--269152 are absent.
- Exact next action / resume command: implement the reentrant per-epoch T05/IGRF snapshot path and threaded MPI-density batches, rebuild according to `AGENTS.md`, validate equivalence, and launch a new corrected QUICK run.

### 2026-10-08T16:10:00Z — multi-epoch T05 threading implementation

- Timestamp correction recorded 2026-10-08T16:05:36Z: the heading time above was entered prospectively in error; this stage was recorded at 16:05:36Z.
- Goal: remove the direct-field serialization bottleneck without changing trajectory physics, task definitions, or scientific acceptance thresholds.
- Repository revision and dirty state: HEAD remains `6104004ca26d77881ad9d05c6a6a21dd9a9577cc`; the pre-existing user `Makefile` change and untracked items remain preserved.
- Files changed (and why): added an immutable reentrant GEOPACK context API in `GeopackInterface.{h,cpp}`; extended `geopack2008.for` context preparation to export epoch-specific dipole tilt; changed the gridless T05 evaluator to call reentrant `t04_s_` plus context IGRF with explicit per-epoch state; narrowed the cutoff multi-epoch fallback to models still using legacy globals; added an MPI-density intra-rank THREADS/OPENMP batch worker pool, main-thread driver preprocessing, actual-worker dynamic chunk sizing, private task results, and deterministic serial accumulation.
- Commands executed, working directories, and return codes: source inspection, patch application, process-safety inspection, deletion-target verification, `git diff --check`, and diff summary returned 0. No AMPS/build process is active.
- Result and status: IMPLEMENTED, NOT YET BUILT OR PHYSICALLY VALIDATED. The worker path contains no MPI calls; multiple direction blocks targeting one `(point,energy)` cannot race because workers write private result slots.
- Scientific/numerical decisions and justification: T05 and IGRF use immutable contexts; other empirical models retain the conservative one-snapshot fallback. GEOPACK RECALC/context preparation is serialized and cached once per unique epoch, while the expensive T05/IGRF field evaluations are concurrent. The exported GEOPACK tilt replaces the gridless evaluator's historical fixed `0.170481` T05 tilt, which requires explicit legacy/context and serial/thread regression checks before physical use.
- Remaining issue or blocker: mandatory clean rebuild and validation are pending.
- Active process state: no SA01, AMPS, MPI, or build process found during the pre-build check.
- Exact next action / resume command: remove the verified non-symlink `/data/vtenishe/GEO/AMPS/build`, run `make -j`, then execute focused field/thread equivalence tests before relaunching SA01.

### 2026-10-08T16:16:17Z — reentrant T05 threading build and equivalence validation complete

- Goal: build and validate the multi-epoch T05/IGRF threading changes before restarting the corrected SA01 QUICK calculation.
- Repository revision and dirty state: HEAD remains `6104004ca26d77881ad9d05c6a6a21dd9a9577cc`; intended changes are in `src/interface/GeopackInterface.{h,cpp}`, `src/interface/geopack2008.for`, `srcEarth/gridless/CutoffRigidityGridless.cpp`, and `srcEarth/gridless/DensityGridless.cpp`; the pre-existing user `Makefile` change and untracked items remain preserved.
- Files changed (and why): completed the immutable context/reentrant field path and worker scheduler recorded in the preceding stage; follow-up fixes added the missing C++ string header, froze driver-table interpolation before density workers enter CSPICE-sensitive code, and serialized only the short GEOPACK/CSPICE context-preparation section while leaving field and trajectory integration concurrent.
- Commands executed, working directories, and return codes: performed four mandatory clean rebuild attempts under the repository-root build policy; attempts 1--3 exposed and corrected compile/CSPICE setup issues, and attempt 4 completed successfully. Final one-rank SERIAL and two-rank/four-thread THREADS cutoff and density smokes all returned 0. `python3 -m unittest discover -s srcEarth/test/SA01/tests -v` passed 10/10. `git diff --check` returned 0.
- Build/executable identity: final executable `/data/vtenishe/GEO/AMPS/amps`, SHA-256 `c84264bed08e70226202b7efebf938c639bd3160e23ad8bf1a08f5a41b6e5bfe`; successful final build log `test_output/SA01/threading_rebuild_20261008/build_attempt4.log`.
- Dataset/input/output identities: focused validation uses two distinct epochs from the SA01 T05 driver table and trajectory fixture under `test_output/SA01/threading_validation_20261008`; it is a bounded engineering regression, not a science comparison.
- Result and status: PASS. Multi-epoch T05 cutoff and density both report `shared backend: THREADS`, `workers/rank: 4`, and complete without the former fallback or CSPICE/GEOPACK race. SERIAL and THREADS output tables are byte-identical.
- Metrics and comparison counts: cutoff 2/2 flattened tasks; density 4/4 tasks and 8/8 trajectories; five density `.dat` files and the cutoff `.dat` file compared byte-for-byte identical; SA01 unit tests 10/10 passed.
- Evidence paths: `test_output/SA01/threading_rebuild_20261008/build_attempt{1,2,3,4}.log` and `test_output/SA01/threading_validation_20261008/{serial_v3,threads_v3,cutoff_serial_v3,cutoff_threads_v3}/amps.log` plus their output tables.
- Scientific/numerical decisions and justification: serial and threaded executions now share the same epoch-specific immutable field context and actual epoch dipole tilt. MPI remains main-thread-only; density accumulation remains deterministic after worker joins. The forthcoming run will use four MPI ranks and eight threads per rank (32 workers on the available 36-CPU set), avoiding the earlier requested 64-way oversubscription.
- Remaining issue or blocker: no full SA01 physical result or observational comparison exists yet; the corrected QUICK run must be prepared and launched into a new output root.
- Active process state: no SA01 physical calculation is active at this checkpoint.
- Exact next action / resume command: prepare `test_output/SA01/quick_threaded_20261008` with the corrected runner, then start it with `--mpi-ranks 4 --threads 8` and verify live logs show THREADS without fallback.

### 2026-10-08T16:16:41Z — fresh threaded QUICK run prepared

- Goal: create immutable inputs for a new full SA01 run using the validated threaded executable.
- Files changed (and why): created new output root `test_output/SA01/quick_threaded_20261008` through the SA01 preparation workflow; no production source changed in this stage.
- Commands executed, working directories, and return codes: `python3 srcEarth/test/SA01/run_SA01.py --prepare-only ... --profile QUICK --mpi-ranks 4 --threads 8 --spacecraft RBSPA RBSPB --formats png eps` returned 0 and reported PASS.
- Build/executable identity: `/data/vtenishe/GEO/AMPS/amps`, SHA-256 `c84264bed08e70226202b7efebf938c639bd3160e23ad8bf1a08f5a41b6e5bfe`.
- Dataset/input/output identities: data root `/data/vtenishe/GEO/sep-in-geospace-sa01-sa04-standalone-20260926`; output root `test_output/SA01/quick_threaded_20261008`; RBSPA input SHA-256 `8f9028e84c8e62137e2b261a8892ead81a1e832736fd2954f1f62442f57fbd76`; T05 driver SHA-256 `79a8bfb9199b4db39be67e20b2da6e5ed55e48fd404574bd66546cf668862be4`; boundary spectrum SHA-256 `dda64f9fa1fa0d3b54a3fdb8112e96846f08079f470ff1322d141fa4a62aa8a1`.
- Result and status: PREPARED/PASS. Both spacecraft branches use `FIELD_MODEL T05`, `FIELD_EVAL_METHOD GRIDLESS`, combined cutoff+density targets, dynamic MPI scheduling, THREADS, and eight workers per rank.
- Metrics and comparison counts: QUICK parameters remain 24 cutoff nodes, 25 spectrum nodes, 96 directions per energy, and 2400 trajectory cap per point. No physical result or comparison yet.
- Evidence paths: `test_output/SA01/quick_threaded_20261008/{preparation_summary.json,data_validation.json,RBSPA/prepared_run.json,RBSPB/prepared_run.json}` and each branch input deck.
- Scientific/numerical decisions and justification: retain direct gridless T05 as the requested physical model but amortize it across 32 worker threads total; do not substitute a gridded interpolation model because that would be a distinct mesh-convergence branch.
- Remaining issue or blocker: physical run has not yet started.
- Active process state: none.
- Exact next action / resume command: invoke the same runner with `--run-only --mpi-ranks 4 --threads 8`, then confirm the solver reports THREADS and eight workers per rank without fallback.

### 2026-10-08T16:21:44Z — corrected threaded QUICK physical run started

- Goal: restart the full two-spacecraft SA01 QUICK calculation after validating the reentrant multi-epoch T05 path.
- Files changed (and why): runner-created RBSPA log/output state plus root `active_run.json`; no source or prepared input changed after launch.
- Commands executed, working directories, and return codes: invoked `run_SA01.py --run-only` from the repository root with four MPI ranks and eight threads per rank; the long command remains active in managed session 5600. Process/log checks returned 0. A stale two-rank process rooted in the earlier failed `threading_validation_20261008/threads` smoke was positively identified by PID, process group, command, and working directory and terminated with SIGTERM; it is no longer present.
- Build/executable identity: `/data/vtenishe/GEO/AMPS/amps`, SHA-256 `c84264bed08e70226202b7efebf938c639bd3160e23ad8bf1a08f5a41b6e5bfe`.
- Dataset/input/output identities: unchanged from the preceding preparation stage; active output root `test_output/SA01/quick_threaded_20261008`.
- Result and status: RUNNING. RBSPA cutoff initialized successfully. Live solver banner reports direct field evaluation, `shared backend: THREADS`, `workers/rank: 8`, DYNAMIC MPI scheduling, and no multi-epoch fallback.
- Metrics and comparison counts: 59 RBSPA cutoff tasks total; automatic chunk is 32, so the two available chunks currently engage two rank-local eight-thread teams while the other ranks await work. Density has 4425 finer tasks and is expected to engage all four rank-local teams. At the checkpoint, active worker ranks showed approximately 491--609% CPU and MPI-waiting ranks approximately 100% CPU.
- Evidence paths: `test_output/SA01/quick_threaded_20261008/{active_run.json,RBSPA/amps.log}`.
- Scientific/numerical decisions and justification: this remains direct/gridless T05, not a gridded-field approximation. Four ranks x eight workers avoids oversubscribing the available 36-CPU affinity set. Preserve the previously interrupted output separately; it is not reused as science data.
- Remaining issue or blocker: the physical calculation and all observation comparisons/plots are pending completion; no blocker has appeared.
- Active process state: managed session 5600; runner PID 522992; MPI PID 523040; AMPS PIDs 523044--523047.
- Exact next action / resume command: poll session 5600 or tail `test_output/SA01/quick_threaded_20261008/RBSPA/amps.log`; confirm cutoff progress/completion and that density also reports THREADS with eight workers per rank.

### 2026-10-08T17:47:00Z — threaded run progress: RBSPA cutoff complete

- Goal: report how far the restarted SA01 calculation has proceeded and preserve its first full threaded density checkpoint.
- Files changed (and why): refreshed the active-run record and this append-only journal; no source, executable, or scientific input changed.
- Commands executed, working directories, and return codes: read-only host process, live log, status, and artifact inspection returned 0.
- Build/executable identity: unchanged threaded executable, SHA-256 `c84264bed08e70226202b7efebf938c639bd3160e23ad8bf1a08f5a41b6e5bfe`.
- Dataset/input/output identities: active output root remains `test_output/SA01/quick_threaded_20261008`; completed RBSPA cutoff table is `RBSPA/cutoff_gridless_points.dat` (7122 bytes).
- Result and status: RUNNING. RBSPA cutoff completed 59/59 and the solver entered density/spectrum with THREADS and eight workers per rank. RBSPB has not started.
- Metrics and comparison counts: latest density progress is 192/4425 MPI direction-block tasks (4.3%), leaving 4233 RBSPA tasks. Latest printed solver ETA is `07:21:18`, which is workload-dependent. Total run elapsed time at the process check was 01:29:58. All four AMPS ranks were active at approximately 676--719% CPU with 9--11 live threads each.
- Evidence paths: `test_output/SA01/quick_threaded_20261008/{active_run.json,RBSPA/amps.log,RBSPA/cutoff_gridless_points.dat}`.
- Scientific/numerical decisions and justification: no change; continue the validated direct/gridless T05 calculation. The completed cutoff is a native model product, but no model/reference comparison or plot exists until density and both spacecraft branches finish.
- Remaining issue or blocker: RBSPA density, the entire RBSPB branch, comparison gates, and plots remain pending. No execution error or blocker is present.
- Active process state: managed session 5600; runner PID 522992; MPI PID 523040; AMPS PIDs 523044--523047.
- Exact next action / resume command: poll session 5600 or tail `test_output/SA01/quick_threaded_20261008/RBSPA/amps.log`; next meaningful milestone is RBSPA density completion and RBSPB launch.

### 2026-10-08T19:29:32Z — active density check: 192/4425 remains latest committed count

- Goal: report current SA01 progress and distinguish a long in-flight density batch from a stopped process.
- Files changed (and why): refreshed `active_run.json` and this progress journal only.
- Commands executed, working directories, and return codes: read-only process, log, artifact, modification-time, and accumulated-CPU inspection returned 0.
- Build/executable identity: unchanged threaded executable, SHA-256 `c84264bed08e70226202b7efebf938c639bd3160e23ad8bf1a08f5a41b6e5bfe`.
- Dataset/input/output identities: active output root remains `test_output/SA01/quick_threaded_20261008`; no new completed product appeared after the RBSPA cutoff table.
- Result and status: RUNNING/COMPUTE-ACTIVE. The latest committed RBSPA density count remains 192/4425 (4.3%); RBSPB has not started. No error or termination record exists.
- Metrics and comparison counts: elapsed wall time 03:12:36. All four AMPS ranks are alive and average 705--749% CPU, with approximately 22:39--24:04 accumulated CPU-hours per rank. The log last changed at 2026-10-08 11:45:37 America/Chicago because progress is published only after a full dynamic batch joins; no later 32-task batch has committed.
- Evidence paths: `test_output/SA01/quick_threaded_20261008/{active_run.json,RBSPA/amps.log,RBSPA/cutoff_gridless_points.dat}` and live PIDs 522992/523040/523044--523047.
- Scientific/numerical decisions and justification: do not treat an unchanged batch counter as completion or failure while workers continue accumulating CPU. Continue the exact direct/gridless T05 calculation without changing inputs mid-run.
- Remaining issue or blocker: current batches contain costly trajectory stragglers, so the earlier `07:21:18` ETA is stale and should not be trusted as a current forecast. RBSPA density, RBSPB, comparison, and plotting remain pending.
- Active process state: managed session 5600; runner PID 522992; MPI PID 523040; all four AMPS ranks live and compute-active.
- Exact next action / resume command: continue monitoring session 5600 or `RBSPA/amps.log`; the next reliable progress change occurs when another rank-local 32-task batch finishes and commits.

### 2026-10-08T19:36:02Z — density/flux acceleration options diagnosed

- Goal: identify concrete ways to accelerate density/flux calculation without silently changing the SA01 science contract.
- Repository revision and dirty state: unchanged during this read-only diagnostic; the active executable and run were not modified.
- Files changed (and why): this append-only journal only.
- Commands executed, working directories, and return codes: inspected density task mapping, MPI chunk resolution, rank-local worker loop, trajectory retry policy, SA01 numerical/sampling controls, and build optimization flags; all read-only commands returned 0.
- Build/executable identity: active executable remains SHA-256 `c84264bed08e70226202b7efebf938c639bd3160e23ad8bf1a08f5a41b6e5bfe`, compiled at `-O3`; `AVXMODE` is currently off and the runtime uses the custom THREADS backend.
- Result and status: DIAGNOSTIC ONLY. The main exposed load-balancing bottleneck is two-level batching: `kTrajBatch=32` makes one density work item contain up to 32 serial trajectories, while automatic MPI chunking gives each rank 32 work items (up to 1024 trajectories) before the rank can fetch more work. `ProcessMappedTaskRange` joins the whole chunk, so completed workers cannot take global work while one costly block remains.
- Metrics and comparison counts: RBSPA contains 59 points x 25 energies x 96 directions = 141600 primary trajectories, represented as 4425 tasks (three 32-direction blocks per point/energy). The active layout is 4 MPI x 8 workers and automatic MPI chunk 32. The same sampling structure gives 158400 primary RBSPB trajectories for 66 points, before unresolved retries.
- Evidence paths: `srcEarth/gridless/DensityGridless.cpp`, `srcEarth/3d/Mode3DParallel.cpp`, `srcEarth/test/SA01/AMPS_PARAM.in.template`, `srcEarth/gridless/READ.ME`, and the active `RBSPA/amps.log`.
- Scientific/numerical decisions and justification: highest-priority science-preserving work is to make direction/trajectory the scheduling unit (or make the current block size configurable and reduce it), reduce the MPI fetch chunk to approximately one local worker team, and eventually use a persistent worker queue. A separately labeled Mode3D epoch-snapshot/interpolated-field branch could remove most repeated T05 evaluations but requires mesh-resolution convergence. Reducing directions, energies, cadence, trace caps, or retry policy changes uncertainty/support and must be treated as an explicitly converged profile rather than a transparent optimization.
- Remaining issue or blocker: no profiling/benchmark matrix has yet quantified each option. The active run is not checkpointable within density, so rebuilding or restarting now would discard its current RBSPA density work; no such action was authorized or taken.
- Active process state: unchanged active run in managed session 5600; all four AMPS ranks remained compute-active at the preceding checkpoint.
- Exact next action / resume command: allow the current run to continue unless explicitly stopped; for the next engineering cycle, benchmark identical serial/thread results using direction block sizes 32/8/4/1 and MPI chunks 32/8, then implement persistent trajectory-level dispatch and compare a Mode3D mesh-resolution series against direct T05.

### 2026-10-08T19:38:22Z — threaded physical run terminated by request

- Goal: stop the current SA01 execution at the user's request and preserve an unambiguous interrupted-run record.
- Files changed (and why): the runner wrote `RBSPA/run_status.json`; `active_run.json` and this journal were updated to record termination. No source, executable, prepared input, or completed cutoff output was changed.
- Commands executed, working directories, and return codes: sent Ctrl-C to managed session 5600; the runner terminated the MPI child and returned the expected interruption status. Host process inspection confirmed runner PID 522992, MPI PID 523040, and AMPS PIDs 523044--523047 are absent.
- Build/executable identity: interrupted run used executable SHA-256 `c84264bed08e70226202b7efebf938c639bd3160e23ad8bf1a08f5a41b6e5bfe`.
- Dataset/input/output identities: output root `test_output/SA01/quick_threaded_20261008`; completed RBSPA cutoff is retained as `RBSPA/cutoff_gridless_points.dat`; final density, flux, spectrum, and termination-summary files do not exist.
- Result and status: INTERRUPTED BY REQUEST, not a completed science result. RBSPA cutoff completed 59/59; latest committed density scheduler count was 192/4425 (4.3%); RBSPB never started.
- Metrics and comparison counts: runner record finished at 2026-10-08T19:38:16Z with return code 1 and `interrupted=true`; four required RBSPA density products are explicitly listed as missing. No comparison or plot was produced.
- Evidence paths: `test_output/SA01/quick_threaded_20261008/{active_run.json,RBSPA/run_status.json,RBSPA/amps.log,RBSPA/cutoff_gridless_points.dat}`.
- Scientific/numerical decisions and justification: preserve the complete cutoff as diagnostic/native output, but exclude the interrupted density work from observational comparison because accumulation/output is transactional and did not complete.
- Remaining issue or blocker: no physical process remains. Finer trajectory scheduling and validation are the next proposed engineering work but were not requested in this termination action.
- Active process state: none; all recorded SA01 processes have exited.
- Exact next action / resume command: before another full run, implement and benchmark smaller direction work units and MPI chunks as recorded in the preceding diagnostic stage; prepare a new output root rather than overwriting this interrupted evidence.

### 2026-10-08T21:54:22Z — repository Earth-test execution inventory reviewed

- Goal: incorporate the repository's established Earth-test inventory and execution conventions before further SA01 scheduler work.
- Repository revision and dirty state: unchanged during this read-only review.
- Files changed (and why): this append-only journal only; `srcEarth/test/list` was not modified.
- Commands executed, working directories, and return codes: read all 224 lines of `srcEarth/test/list` and searched its setup, F3, and SA01 entries; commands returned 0.
- Result and status: DIAGNOSTIC/PASS. The established suite is compiled with `make test_Earth-cutoff-rigidity-3d-analytic_compile -j` and orchestrated by `./srcEarth/test/test_runner.py -j 6 ./srcEarth/test/list`. Independent unit/contract entries and physical C/F test commands carry per-entry last-pass provenance.
- Metrics and comparison counts: the active F3 density/flux validation command uses HC4, 24 scan nodes, 288 directions per energy, 600-s trace cap, 1000000 steps, unresolved retry/tolerance controls, `--dynamic-chunk 1`, and 18 MPI ranks. The list currently contains no SA01 entry.
- Evidence paths: `srcEarth/test/list`, especially lines 1--3 and 117--118.
- Scientific/numerical decisions and justification: future density scheduling changes must be exercised through the existing F3/list workflow, not only through ad-hoc SA01 smokes. The recorded F3 layout also demonstrates the established strategy for highly variable density trajectories: fine MPI chunking (`1`) and more MPI ranks, which should inform the next scheduler benchmark.
- Remaining issue or blocker: no tests were executed and no code was changed in this review-only stage. An SA01 list entry should be added only when its bounded command/runtime and pass criteria are defined.
- Active process state: no SA01 run is active.
- Exact next action / resume command: before implementing scheduler changes, inspect the F3 runner defaults/contracts and the test-runner entry semantics; after implementation, run the focused unit tests and F3 command through the established test infrastructure, then add a bounded SA01 regression entry if appropriate.

### 2026-10-08T22:05:20Z — Earth regression failures and missing-output behavior diagnosed

- Goal: investigate the reported test-list regressions and output-tree differences without weakening gates or changing production/test behavior.
- Repository revision and dirty state: no production source, runner, list, executable, or test gate was modified; this journal entry is the only write.
- Commands executed, working directories, and return codes: inspected all three last-pass regression logs, the complete runner reports, C10's additional actual failure, current/prior lists and output roots, runner subprocess I/O behavior, Step-12 source contract/history, runner defaults, executable identities, artifact timestamps, and live process state. Read-only commands returned 0; no test process remains active.
- Build/executable identity: correct current Earth executable is `/data/vtenishe/GEO/AMPS/amps`, SHA-256 `c84264bed08e70226202b7efebf938c639bd3160e23ad8bf1a08f5a41b6e5bfe`, built 2026-10-08. The C9/C10 list entries instead select `/home/vtenishe/T11/AMPS/amps`, SHA-256 `dd058d678fbe44d2e23431817add2ab22aa1804362751e8a1b25d7a739289498`, rebuilt 2026-09-28 as a SEP executable after the 2026-09-26 last-pass commit.
- Result and status: DIAGNOSTIC COMPLETE; NO FIX IMPLEMENTED. Of 61 entries, 52 actually passed and 9 actually failed; five of the failures are expected-F gates, leaving four mismatches. `UStep12Release` fails one source-contract assertion; both C9 entries and C10 fail before physics because the stale absolute executable rejects Earth option `-mode`.
- Metrics and comparison counts: C9 FULL_SCAN and DIRECT_ACCESS exit after 3.192 s and 2.406 s; C10 exits after 0.727 s. Current and prior top-level output directories are identical except prior-only `C10_direct` and current-only `SA01`. C9 current roots contain reference/driver/deck/error-log setup only. C10 explicitly writes into the prior tree, where a 2026-09-26 result coexists with the new 2026-10-08 failed-launch log.
- Evidence paths: `test_runner_report{,_to_address,_actual_failed,_last_pass_failed}.{txt,csv,json}`, `test_runner_logs/test_{011,051,052,053}_*.log`, `srcEarth/test/list`, `srcEarth/test/UStep12Release/test_step12_source_contract.py`, `srcEarth/test/test_runner.py`, current `test_output/C9_{full,direct}_85`, and prior `/home/vtenishe/T11/AMPS/test_output/C10_direct`.
- Scientific/numerical decisions and justification: no physical gate failed because C9/C10 never entered Earth physics. Repair should keep every existing numerical threshold and command role: use the repository-local Earth executable and repository-local output roots, add an Earth-executable preflight, and retain fail-closed exit behavior. The Step-12 contract is self-invalidating after the runner updates mutable `last pass:` provenance; repair should normalize/validate provenance separately while retaining an immutable structural baseline for commands/P-F markers and explicitly reviewed list maintenance.
- Remaining issue or blocker: two user-intent points need confirmation before editing: whether `/data/vtenishe/GEO/AMPS`/`./amps` is the authoritative active test installation, and whether “stopped producing output” refers to missing artifacts or lack of live terminal streaming. The runner intentionally sends child stdout/stderr only to per-test logs and prints only start/completion records to the terminal.
- Active process state: no Earth test runner, MPI test, or SA01 process is active.
- Exact next action / resume command: after clarification, patch portable C9/C10 paths, add executable/output provenance preflights and transactional failure status, repair the Step-12 structural contract without changing scientific gates, run focused unit/source-contract tests, then rerun C9/C10 through `test_runner.py` with their original gate arguments.

<!-- Append one section after each meaningful stage or before a handoff.
Use actual timestamps and observed evidence. Remove no previous session entry.

### YYYY-MM-DDTHH:MM:SSZ — stage/session title

- Goal:
- Repository revision and dirty state:
- Files changed (and why):
- Commands executed, working directories, and return codes:
- Build/executable identity if relevant:
- Dataset/input/output identities if relevant:
- Result and status: PASS / FAIL / ERROR / BLOCKED / diagnostic only.
- Metrics and comparison counts (actual values only):
- Evidence paths: logs, CSV/JSON, figures, report.
- Scientific/numerical decisions and justification:
- Remaining issue or blocker:
- Active process state (PID/job identity and command, if any):
- Exact next action / resume command:
-->

### 2026-10-08T22:22:09Z — portable Earth-test launch and provenance repairs implemented

- Goal: implement the diagnosed C9/C10/list/Step-12 repairs while preserving concise test-runner start/completion reporting and every existing scientific acceptance gate.
- Repository revision and dirty state: the pre-existing SA01 scheduler/Geopack/Makefile work remains dirty and was not reverted; this stage adds focused Earth-test runner, list, contract, unit-test, and README changes.
- Files changed (and why): added `srcEarth/test/earth_test_runtime.py` and `UTestRunner/test_earth_test_runtime.py`; updated C9/C10 runners for Earth-application preflight, executable SHA-256 provenance, atomic run status, and recoverable prior-result archival; updated `srcEarth/test/list` to use checkout-local `./amps` and output roots; updated the Step-12 list contract; updated C9, C10, UTestRunner, and UStep12Release READMEs and the UTestRunner entry point. The production `test_runner.py` output policy was intentionally not changed.
- Commands executed, working directories, and return codes: Python syntax compilation returned 0; the new UTR-F06 tests returned 0 (2/2); the focused Step-12 source contract returned 0 (5/5); probing the real checkout executable returned 0. A combined UTestRunner/Step-12 command was interrupted after the existing deferred-last-pass test produced no completion for more than two minutes; a bounded 45-second rerun of that pre-existing test also timed out, so it is being investigated separately rather than reported as a pass.
- Build/executable identity: `/data/vtenishe/GEO/AMPS/amps` was positively identified from its `-h` contract as `AMPS Earth: Gridless Energetic Particle Solver`; SHA-256 `c84264bed08e70226202b7efebf938c639bd3160e23ad8bf1a08f5a41b6e5bfe`.
- Dataset/input/output identities: no physical C9, C10, or SA01 solver was launched in this stage; no scientific output or previous result tree was overwritten.
- Result and status: IMPLEMENTATION COMPLETE, FOCUSED TESTS PASS, broader focused-suite investigation still in progress. C9/C10 now fail before MPI for an incorrect AMPS application, and each output root has an authoritative `C9_run_status.json` or `C10_run_status.json` whose `result_current` is false until a result is finalized.
- Metrics and comparison counts: Python compile 5 files passed; shared runtime tests 2/2 passed; Step-12 source-contract tests 5/5 passed; the first 9 exclusive-runner tests in the combined UTestRunner invocation passed before the unrelated deferred test stalled.
- Evidence paths: modified source/tests/READMEs above; focused command output in this session; status schema and comments in `srcEarth/test/earth_test_runtime.py`.
- Scientific/numerical decisions and justification: no tolerance, observation, reference, grid, trajectory policy, P/F expectation, or comparison gate changed. The Step-12 digest now removes only syntactically valid runner-owned commit values; all command text, activation state, ordering, expectations, and comments remain hash-frozen. Live child output remains captured in per-test logs, with terminal reporting left at start/completion granularity as requested.
- Remaining issue or blocker: exercise C9/C10 dry-run and wrong-binary integration paths, determine why the pre-existing deferred-last-pass unit test stalls in this execution environment, and run the complete Step-12 validation suite.
- Active process state: none; both bounded test commands have exited and no physical solver is running.
- Exact next action / resume command: run C9/C10 smoke dry-runs into temporary roots, inspect their status JSON, then run `test_release_validation.py` and diagnose the deferred-runner timeout without launching the expensive physical list entries.

### 2026-10-08T23:35:40Z — runner completion fix and focused regression verification complete

- Goal: finish the non-scientific execution/provenance repairs, retain concise start/completion terminal reporting, and verify them without launching the expensive C9/C10 physics calculations.
- Repository revision and dirty state: existing SA01 numerical/threading changes remain intact; the final focused change set adds the shared runtime module, C9/C10 wiring/docs, portable list entries, Step-12 structural contract/docs, test-runner completion heartbeat/docs, and UTestRunner coverage.
- Files changed (and why): `srcEarth/test/test_runner.py` now wakes its asyncio process wait every 0.25 s so an already-queued child-exit notification cannot stall an unlimited run; child output remains log-only. `UTestRunner/test_exclusive_runner.py` now verifies heartbeat-compatible timeout cleanup. `earth_test_runtime.py`, `UTestRunner/earth_test_runtime_test.py`, C9/C10 runners, list, source contract, and relevant READMEs retain the executable/provenance repairs recorded in the prior stage.
- Commands executed, working directories, and return codes: C9 and C10 SMOKE dry-runs into `/tmp/amps-c{9,10}-provenance-smoke` returned 0; C9 and C10 invocations selecting the known SEP executable returned the required code 2 and the enclosing assertion command returned 0; `UTestRunner/run_test.sh` returned 0; `UStep12Release/run_test.sh` returned 0; `git diff --check` returned 0.
- Build/executable identity: correct Earth binary remains `/data/vtenishe/GEO/AMPS/amps`, SHA-256 `c84264bed08e70226202b7efebf938c639bd3160e23ad8bf1a08f5a41b6e5bfe`. Wrong-application integration used `/home/vtenishe/T11/AMPS/amps` and was rejected from missing all three required Earth help tokens before MPI.
- Dataset/input/output identities: dry-run evidence is under `/tmp/amps-c9-provenance-smoke` and `/tmp/amps-c10-provenance-smoke`; rejected-launch evidence is under `/tmp/amps-c9-wrong-app-smoke` and `/tmp/amps-c10-wrong-app-smoke`. No repository scientific output root was overwritten and no physical solver ran.
- Result and status: PASS for all focused tests. Dry-run status files end `DRY_RUN_COMPLETE`, code 0, `result_current=false`; wrong-application status files end `ERROR`, code 2, `result_current=false`, and contain no solver artifacts. The previously stalled deferred-last-pass suite now completes in about 2.9 s while retaining log-only child output.
- Metrics and comparison counts: UTestRunner passed 10 scheduler/environment/timeout tests, 3 deferred-provenance tests, and 3 Earth runtime/wiring tests (16/16 total). Step-12 passed 8 release-validation and 5 source-contract tests (13/13 total). Both C9/C10 dry runs generated two SMOKE snapshot commands each; both wrong-application checks rejected 1/1 invocations.
- Evidence paths: `/tmp/amps-c9-provenance-smoke/C9_run_status.json`, `/tmp/amps-c10-provenance-smoke/C10_run_status.json`, `/tmp/amps-c9-wrong-app-smoke/C9_run_status.json`, `/tmp/amps-c10-wrong-app-smoke/C10_run_status.json`, and the modified focused test files.
- Scientific/numerical decisions and justification: all C9 PAMELA and C10 POES gates, reference products, shell/mesh settings, rigidity sets, trajectory controls, expected P/F markers, and tolerances are unchanged. The runner heartbeat changes only process-completion observation; it does not read or print child output, impose a new timeout, or change process-group timeout cleanup.
- Remaining issue or blocker: the repaired active C9/C10 entries still require their normal expensive physical reruns before their last-pass provenance can legitimately advance; those calculations were intentionally not started in this focused repair stage.
- Active process state: none; no test runner, MPI launcher, AMPS solver, or SA01 process remains active.
- Exact next action / resume command: execute the normal list workflow when compute time is available; the repaired runner will show `[START]` and `[OK]`/`[MISMATCH]` while detailed child output remains in `test_runner_logs`, and C9/C10 status files will distinguish in-progress, error, and current results.

### 2026-10-09T06:51:38Z — Mode3D single-location density parallelization plan completed

- Goal: define a science-preserving implementation plan that allows one gridded
  Mode3D density location/snapshot to use multiple MPI ranks and multiple thread
  workers without sharing mutable field-evaluator or accumulation state.
- Repository revision and dirty state: HEAD
  `648f06f341d6e2f2cc9ff18d4ef11633a1142e69`; observed dirty state is the modified
  `Makefile` plus untracked `SWMF_data`, `amps`, `gitinfo.txt`, `srcEarth/test/SA01/`,
  and `srcEarth/test/__pycache__/`. No production source was changed in this stage.
- Files changed (and why): this progress journal only, to preserve the requested
  design and exact next action.
- Commands executed, working directories, and return codes: reread `AGENTS.md`,
  `SA01_TASK.md`, and the complete progress record; inspected current process/repository
  state; traced `DensityMode3D.cpp`, `CutoffRigidityMode3D.cpp`, `Mode3D.cpp`,
  `Mode3DParallel.{h,cpp}`, `FluxNumerics.h`, the Earth test list, and relevant READMEs.
  All read-only commands returned 0. No solver, build, or test was launched.
- Build/executable identity: current untracked `amps` is 6,008,896 bytes, SHA-256
  `fea7c6dc390c7c2d5b273657b1a33335a92fca9e787e5edcdc73c7e5624bda05`, mtime
  2026-10-08 18:59:22 -0500; it was not used and current source correspondence was
  not asserted in this planning stage.
- Result and status: PLAN COMPLETE; IMPLEMENTATION/VALIDATION PENDING. The current
  Mode3D THREADS backend schedules whole locations through a rank-local atomic counter.
  OpenMP can enter energy/direction loops, but those loops are suppressed inside direct
  thread workers. With one trajectory sample selected for an epoch, one MPI rank and one
  direct worker therefore own the useful density work.
- Scientific/numerical decisions and justification: the recommended inter-rank work
  item is `(location,energy)`, not a whole location and not a cross-rank partial direction
  block. The owning MPI rank dynamically distributes that energy's directions among its
  worker threads, records one result per direction, and reduces those results in canonical
  direction order. One rank owns each complete `AccessAccumulator`, so existing
  `ResolveAccess` semantics and bitwise-stable per-energy aggregation can be retained;
  the global post-scheduler reduction combines disjoint `(location,energy)` slots. Only
  after all energy slots are present are spectrum, density, flux, channel, and detector
  products evaluated. Snapshot construction/advancement remains outside the worker
  region and the compact field arrays remain immutable.
- Validation requirements: add pure task encode/decode and ownership tests; require
  exact SERIAL versus THREADS output/count equivalence for a one-location Mode3D case;
  require one-rank versus multi-rank and DYNAMIC/STATIC/BLOCK_CYCLIC equivalence; cover
  unresolved retry/termination accounting and anisotropic weights; prove more than one
  worker and, when `nE>1`, more than one MPI rank executed useful work; run existing
  Step 11/12 product contracts and the complete relevant Earth regression entries without
  relaxing any gate. F3 remains required scientific density evidence but is gridless and
  cannot substitute for a new bounded Mode3D density scheduler regression.
- Evidence paths: production call path begins at `srcEarth/3d/DensityMode3D.cpp`;
  snapshot isolation is in `srcEarth/3d/Mode3D.cpp`; shared scheduler support is in
  `srcEarth/3d/Mode3DParallel.{h,cpp}`; accumulator semantics are in
  `srcEarth/util/FluxNumerics.h`.
- Remaining issue or blocker: no code has been implemented, no test ID has been assigned
  for the new bounded regression, and no performance result exists. A persistent worker
  pool and finer cross-rank direction blocks are deferred until the correctness-first
  `(location,energy)` design is measured; they are not prerequisites for useful
  single-location parallelism.
- Active process state: none matching SA01, AMPS, MPI launchers, the Earth test runner,
  or a build was observed.
- Exact next action / resume command: implement the scheduler/refactor and tests, update
  the Mode3D parallel documentation, then follow the mandatory verified-build-directory
  removal plus repository-root `make -j` procedure before executing the validation matrix.

### 2026-10-09T07:01:44Z — GRIDDED QUICK repeat prepared and clean rebuild completed

- Goal: repeat the interrupted SA01 QUICK physical calculation with the mesh-backed
  Mode3D/GRIDDED field representation requested by the user, without overwriting the
  interrupted GRIDLESS evidence.
- Repository revision and dirty state: HEAD
  `648f06f341d6e2f2cc9ff18d4ef11633a1142e69`; pre-existing modified/untracked items
  were preserved. SA01 runner/template/tests/README and this journal were extended for
  the new solver branch.
- Files changed (and why): `AMPS_PARAM.in.template`, `prepare_SA01.py`, `run_SA01.py`,
  `compare_SA01.py`, `tests/test_sa01.py`, and `README.md` now support an explicit
  GRIDDED solver, SNAPSHOT_LIST batching, explicit mesh controls, log-only child output,
  per-snapshot completion checks, and comparison discovery of suffixed Mode3D products.
- Commands executed, working directories, and return codes: Python syntax compilation
  returned 0; focused SA01 tests passed 11/11; GRIDDED QUICK preparation returned 0;
  `git diff --check` returned 0. The verified non-symlink live `build` directory was
  removed and repository-root `make -j` returned 0 with output captured in the build log.
- Build/executable identity: `/data/vtenishe/GEO/AMPS/amps`, 6,008,768 bytes, SHA-256
  `c84264bed08e70226202b7efebf938c639bd3160e23ad8bf1a08f5a41b6e5bfe`, mtime
  2026-10-09 02:01:33 -0500; build log
  `test_output/SA01/threading_rebuild_20261009/build.log`.
- Dataset/input/output identities: new output root
  `test_output/SA01/quick_gridded_20261009`; RBSPA has 59 selected epochs and RBSPB
  66. Both decks specify GRID_3D, SNAPSHOT_LIST, `0.05 Re` near-Earth resolution,
  `2.0 Re` boundary resolution, LINEAR coarsening, dynamic MPI chunk 1, and OPENMP
  with eight workers per rank.
- Result and status: PREPARATION/BUILD PASS; PHYSICAL RUN NOT YET STARTED. The selected
  mesh is an initial diagnostic resolution and does not satisfy the required mesh-
  convergence gate by declaration; convergence must be measured on common epochs.
- Scientific/numerical decisions and justification: every trajectory row is evaluated
  only in its matching frozen T05 snapshot. OPENMP is used because the current direct
  THREADS Mode3D density scheduler owns a whole location, whereas OpenMP enters the
  energy/direction loops for the one active location. No physical sampling count,
  unresolved rule, observation selection, or acceptance threshold changed.
- Evidence paths: preparation JSON/decks under the new output root; build log and the
  modified SA01 files listed above.
- Remaining issue or blocker: first-snapshot runtime, mesh size, and physical product
  behavior have not yet been observed. No comparison or plot exists.
- Active process state: no SA01/AMPS/MPI/build process was active at this checkpoint.
- Exact next action / resume command: run the prepared branches sequentially with
  `run_SA01.py --run-only ... --mpi-ranks 4 --threads 8`, inspect the first RBSPA
  snapshot startup and retain the run only if it is truly INTERPOLATION/GRIDDED.

### 2026-10-09T07:05:03Z — First GRIDDED snapshot diagnosed; serial fallback stopped

- Goal: verify that the prepared GRIDDED repeat actually uses mesh interpolation and
  threaded trajectory work before committing to the complete 125-epoch calculation.
- Commands executed and result: launched the prepared RBSPA branch with four MPI ranks
  and eight requested density threads. Snapshot 1/59 constructed a T05 GSM compact
  field with 27,840 leaf blocks and 6,013,440 stored cells; cutoff completed and wrote
  its first epoch product. The executable banner then reported that it had been built
  without OpenMP, and both cutoff and density selected SERIAL. The diagnostic was
  deliberately terminated before serial density could consume the long production
  runtime. No matching SA01/AMPS/MPI process remained afterward.
- Evidence paths: `test_output/SA01/quick_gridded_20261009/RBSPA/amps.log`,
  `active_run.json`, `RBSPA/run_status.json`, and the first cutoff/manifest products.
  No density product was completed, so this interrupted root cannot be used for the
  requested comparison.
- Build correction stage: changed repository `Makefile` from `OPENMP=off` to
  `OPENMP=on`, preserving the pre-existing `EXTRALINKEROPTIONS_F=-lmpi++` modification.
  After verified removal of the exact non-symlink `build` directory, repository-root
  `make -j` failed at `src/pic/pic_mover.cpp`: `_OPENMP` enabled calls to
  `omp_get_thread_num()` and `omp_get_num_threads()`, while `pic.h` only included
  `omp.h` for AMPS HYBRID mode. This is a genuine supported-build defect rather than a
  test-gate issue. The guarded use-site include has now been added.
- Build evidence: failed build log
  `test_output/SA01/threading_rebuild_20261009/build_openmp.log`.
- Status and exact next action: DIAGNOSTIC COMPLETE; OPENMP REBUILD FIX IN PROGRESS.
  Perform another mandatory clean repository-root build, prepare a fresh output root,
  and relaunch only after the runtime banner confirms INTERPOLATION plus OPENMP backends.

### 2026-10-09T07:11:51Z — Threaded GRIDDED repeat launched and runtime-verified

- Goal: rebuild the executable with working OpenMP support, start the SA01 QUICK
  GRIDDED repeat in a clean evidence root, and verify its actual runtime backend before
  allowing the long calculation to continue.
- Files changed (and why): `src/pic/pic_mover.cpp` now includes `omp.h` at the guarded
  `_OPENMP` use site because `OPENMP=on` does not imply AMPS HYBRID mode;
  `srcEarth/gridless/DensityGridless.cpp` and `srcEarth/3d/DensityMode3D.cpp` snapshot
  immutable spectrum-unit metadata outside their `default(none)` worker regions. This
  both fixes strict OpenMP compilation and keeps mutable spectrum epoch selection out of
  threaded work. Detailed comments record both invariants. `Makefile` remains
  `OPENMP=on`; the pre-existing Fortran linker option remains preserved.
- Commands and validation: a second clean build exposed the undeclared global spectrum
  dependency in three density OpenMP loops and failed (log `build_openmp_fixed.log`).
  After the explicit unit snapshots were added, `git diff --check` passed. The exact
  non-symlink `build` directory was again verified and removed, and repository-root
  `make -j` completed successfully. The final link contains `-fopenmp`, and `ldd amps`
  resolves `libgomp.so.1`.
- Build identity: `/data/vtenishe/GEO/AMPS/amps`, 6,029,056 bytes, SHA-256
  `157c584db4d39a5ff6eeebc6f554aa32fd8c6582b9e41d34af14594d3ef6ebf1`, mtime
  2026-10-09 02:10:02 -0500. Successful build log:
  `test_output/SA01/threading_rebuild_20261009/build_openmp_fixed2.log`.
- Fresh evidence root: `test_output/SA01/quick_gridded_openmp_20261009`; preparation
  passed for RBSPA (59 epochs) and RBSPB (66 epochs). The decks retain `GRID_3D`,
  `SNAPSHOT_LIST`, 0.05/2.0 Re LINEAR mesh controls, OPENMP with eight requested
  workers, and dynamic MPI chunk one.
- Runtime verification: RBSPA snapshot 1/59 prepared a T05 GSM compact field with
  27,840 leaf blocks and 6,013,440 field cells. The runtime banner reports
  `OpenMP threads : 8 per rank`, `Cutoff backend : OPENMP`, and
  `Density backend : OPENMP, threads/MPI rank=8`. The first cutoff product completed;
  density is actively processing 2,400 trajectory tasks for the first epoch. This is
  mesh interpolation, not direct/gridless T05 tracing.
- Active process state: runner PID 1130510 and launcher PID 1130560 are recorded in
  `test_output/SA01/quick_gridded_openmp_20261009/active_run.json`; the managed execution
  session is 99469. Start time is 2026-10-09T07:10:35Z. RBSPA is running; RBSPB will
  start sequentially after RBSPA completes.
- Status and exact next action: RUNNING AND BACKEND-VERIFIED. Continue monitoring
  snapshot completion/product counts from `RBSPA/amps.log` without streaming child
  output. Do not compare or plot until all 59 RBSPA and 66 RBSPB epoch products pass
  the runner's exact-count checks. The 0.05/2.0 Re mesh remains diagnostic and does not
  satisfy mesh-convergence merely by completing.

### 2026-10-09T07:13:26Z — GRIDDED baseline stopped by strict unresolved gate

- Goal: determine whether the backend-verified GRIDDED repeat can complete the first
  physical epoch under the unchanged packaged 300-s trace limit and 0.01 unresolved
  acceptance threshold.
- Result and status: PHYSICAL/NUMERICAL FAIL, not an interruption and not PASS. RBSPA
  snapshot 1/59 completed all 2,400 primary density trajectory tasks. The termination
  transaction measured a maximum unresolved fraction of
  `0.22916666666666666 > 0.01` and invoked the configured fatal gate before writing
  spectrum, flux, or density products. RBSPB never started. No process remains active.
- Termination evidence: at 18.5 MeV, 22 of 96 directions remained `TIME_LIMIT` after
  the configured retry; 29 reached the outer allowed boundary, 43 reached the inner
  forbidden boundary, and 2 were classified magnetically trapped. Other low-energy
  rows also contain time-limit terminations. Numerical/field failures and step-limit
  terminations are zero, so this is specifically insufficient physical trace-time
  closure, not an OpenMP crash or invalid interpolated field.
- Gate handling: `DS_UNRESOLVED_TOL 0.01` and `DS_FAIL_ON_UNRESOLVED T` were preserved.
  TIME_LIMIT was not converted to forbidden, no rows/channels were omitted, and the
  observation selection was not changed. The runner additionally rejected AMPS's
  nominal return code zero because exact artifact counts were incomplete (1/59 cutoff,
  1/59 termination, 1/59 manifest, and 0/59 spectrum/flux/density).
- Exit-status defect: the fatal gate called the repository's generic `exit(...)`, whose
  current packed error code was `0.0`; that propagated to `MPI_Abort(..., 0)` and made
  `mpiexec` report success despite the fatal message. The SA01 artifact gate still
  caught the failure. This misleading zero status is a separate production error-path
  defect and must be fixed without weakening the scientific gate.
- Backend/performance evidence: field construction was INTERPOLATION/GRIDDED with
  27,840 leaf blocks, 6,013,440 compact cells, four MPI ranks, and parallel field
  initialization. Cutoff and density both reported OPENMP with eight workers. The
  complete first density workload reached 2,400/2,400, proving that the new OpenMP
  build and within-location threaded density path executed before the scientific gate.
- Evidence paths: `test_output/SA01/quick_gridded_openmp_20261009/RBSPA/amps.log`,
  `RBSPA/run_status.json`, `RBSPA/mode3d_termination_summary_snapshot_000000_2014_01_06T09_00_00Z.dat`,
  and root `active_run.json`.
- Remaining decision / exact next action: the package baseline explicitly uses 300 s,
  so it must remain recorded rather than silently replaced. A completion attempt now
  requires either (a) a separately labelled trace-limit convergence run with longer
  `DS_MAX_TRAJ_TIME`/`MAX_TRACE_TIME` and proportional `MAX_STEPS`, retaining the 0.01
  gate, or (b) a focused investigation of why long-lived low-energy trajectories are
  not closed by the existing trap classifier. Obtain the user's choice before changing
  that scientific numerical control.

### 2026-10-09T07:17:37Z — Non-weakening resolution options identified

- Goal: distinguish scientifically valid remedies for the first-epoch TIME_LIMIT
  failure from changes that would merely weaken or bypass SA01 validation.
- Evidence reviewed: the Mode3D density retry currently reruns every unresolved
  direction from its initial state with half `DT_TRACE`, twice `MAX_STEPS`, and twice
  the 300-s time budget. The termination summary contains TIME_LIMIT outcomes only for
  the unresolved population; it contains no step, distance, invalid-field, invalid-step,
  or numerical-failure outcomes. The SA01 deck enables the mirror/bounce trap detector
  but leaves the existing drift-shell recurrence detector disabled by default.
- Valid remedies, in recommended order: (1) run a separately labelled trace-time
  convergence sweep while retaining `DS_UNRESOLVED_TOL=0.01`; (2) enable the existing
  C19-validated drift-shell recurrence settings and require positive recurrence before
  classifying a trajectory as trapped; (3) improve density retry into deterministic
  unresolved-only continuation/extension with full initial/final limit provenance;
  (4) investigate mesh-resolution sensitivity against direct T05 on identical
  trajectories; and (5) make the fatal path return a nonzero MPI/process status while
  keeping the artifact-count gate.
- Explicitly rejected shortcuts: lowering the unresolved gate, setting
  `DS_FAIL_ON_UNRESOLVED=F` for a validation PASS, mapping every TIME_LIMIT to forbidden,
  dropping low-energy directions/channels, or publishing nominal values without their
  unresolved bounds. Increasing only `MAX_STEPS` is also ineffective for the observed
  TIME_LIMIT-only failure.
- Status / next action: DESIGN OPTIONS RECORDED; no new run or parameter change was
  made. Choose between trace-time convergence first (lowest semantic risk) and a
  drift-detector qualification run (potentially much cheaper) before preparing a new,
  separately fingerprinted output root.

### 2026-10-09T07:22:24Z — Non-weakening density-closure fix implemented

- Goal: resolve the observed Mode3D TIME_LIMIT population without lowering the 0.01
  unresolved gate, remapping generic timeouts to forbidden, or changing sampling and
  observation selections.
- Production changes: `amps_param_parser.{h,cpp}` now parses and validates bounded
  `DS_UNRESOLVED_EXTENSION_PASSES` and `DS_UNRESOLVED_EXTENSION_FACTOR` controls.
  Gridless and Mode3D density separate genuine numerical-failure retry from physical
  time convergence: invalid integration may receive one smaller-step retry, while only
  TIME_LIMIT/STEP_LIMIT directions restart from the identical seed at explicitly larger
  budgets with a scaled step allowance. DISTANCE_LIMIT is never relaxed. Structured
  result fields retain primary termination/time and initial/final extension limits.
- Failure-path correction: all three density unresolved-gate fatal paths now set a
  nonzero packed AMPS function code immediately before `exit(...)`, preventing
  `MPI_Abort(...,0)` from reporting a rejected validation transaction as success.
- SA01 controls: the 300-s packaged primary trace remains unchanged. The rendered deck
  enables three factor-two unresolved-only extensions (600, 1200, 2400 s) and enables
  the existing C19-qualified drift-shell recurrence classifier with its three-turn,
  profile-coverage, phase-space, and secular-radius checks. Positive recurrence is
  required; arbitrary timeout is still unresolved. Prepared provenance now records the
  closure policy explicitly.
- Documentation/tests: detailed implementation comments were added to all modified
  paths; `srcEarth/README.md`, `3d/README.md`, `gridless/READ.ME`, and the SA01 README
  document semantics and rejected shortcuts. SA01 render tests now require the strict
  gate, extension controls, and drift detector.
- Checks executed: `git diff --check` passed; all 11 SA01 Python tests passed; the
  dependency-free UTrajectoryCore U-F10..U-F13 policy suite passed.
- Status / exact next action: IMPLEMENTED; CLEAN REBUILD PENDING. Verify that no AMPS,
  MPI, SA01, or build process is active; verify the exact non-symlink live `build`
  directory; remove only it; run repository-root `make -j` with a new captured log.

### 2026-10-09T07:24:28Z — Density-closure build and fresh GRIDDED preparation passed

- Goal: compile the new closure policy with OpenMP and prepare immutable inputs for a
  clean physical restart without overwriting either failed GRIDDED baseline.
- Build procedure/result: process inspection found no active AMPS/MPI/SA01/build work;
  `/data/vtenishe/GEO/AMPS/build` resolved exactly to a non-symlink directory. It alone
  was removed, then repository-root `make -j` returned 0 with no compiler/linker errors.
  The final link includes `-fopenmp`, and `ldd` resolves `libgomp.so.1`.
- Build identity: `/data/vtenishe/GEO/AMPS/amps`, 6,029,456 bytes, SHA-256
  `d1b771bda423c4da144e8045cfd447a370362a2c26346071bd7e5921b556cee9`, mtime
  2026-10-09 02:24:08 -0500. Build log:
  `test_output/SA01/threading_rebuild_20261009/build_density_closure.log`.
- Preparation result: PASS in fresh root
  `test_output/SA01/quick_gridded_closure_20261009`. Both branches render GRID_3D,
  SNAPSHOT_LIST, OPENMP/8 workers, the unchanged 0.01 fatal gate, three factor-two
  extensions, and the full C19 drift recurrence controls. `prepared_run.json` records
  the primary 300-s and maximum 2400-s budgets plus recurrence provenance.
- Status / exact next action: BUILD/PREPARATION PASS; PHYSICAL RESTART PENDING. Launch
  `run_SA01.py --run-only` with four MPI ranks and eight threads, inspect snapshot 1
  for parser acceptance/OpenMP/INTERPOLATION, then require the previously failing first
  density epoch to produce zero-or-acceptable unresolved support before continuing.

### 2026-10-09T07:27:30Z — Repaired GRIDDED restart cleared the failed first epoch

- Goal: prove that the closure fix resolves the exact epoch that failed at 0.2291667
  unresolved fraction before allowing the multi-epoch calculation to continue.
- Result: PASS for RBSPA snapshot 1/59. Field evaluation is INTERPOLATION/GRIDDED;
  cutoff and density report OPENMP with eight workers. Density completed 2,400/2,400
  primary direction/energy tasks and wrote cutoff, density, spectrum, flux, termination,
  and manifest products. Maximum final unresolved fraction is exactly zero against the
  unchanged 0.01 threshold.
- Termination comparison: the original 18.5-MeV row had 22 TIME_LIMIT outcomes after
  its single legacy retry. The repaired row has 96/96 resolved: 34 allowed outer exits,
  60 inner losses, 2 mirror/bounce trapped, and zero TIME_LIMIT/STEP_LIMIT/numerical
  failures. `N_retried=41` records directions that required additional integration;
  no generic timeout was counted as trapped.
- Work distribution: each snapshot contains exactly one active trajectory location.
  RBSPA has 59 locations/epochs and RBSPB 66 (125 total), but they are processed as
  one-location snapshots to preserve immutable epoch-specific grids. Consequently one
  MPI rank owns density for a snapshot and uses eight OpenMP workers (~800% CPU); the
  other ranks participate in distributed field construction but have no second density
  location to claim.
- Active state: snapshot 2/59 density is running in managed session 14675. Runner PID
  1146912, launcher PID 1146966, and AMPS PIDs 1146976--1146979 are recorded in
  `test_output/SA01/quick_gridded_closure_20261009/active_run.json`.
- Exact next action: keep the run active, inspect each termination summary for the
  unchanged gate, and report the next meaningful milestone rather than streaming child
  output. RBSPB starts only after all 59 RBSPA snapshots pass exact artifact counts.

### 2026-10-09T07:30:39Z — Epoch-level MPI concurrency design clarified

- Goal: explain why 59 RBSPA epochs currently engage only one rank for density and
  whether dynamic one-epoch-per-rank scheduling is technically valid.
- Finding: the design is valid because MPI ranks have independent address spaces; each
  can own one epoch-specific field while its OpenMP threads share that rank-local frozen
  field. The current implementation cannot diverge epochs because `Mode3D::Run` advances
  the snapshot loop collectively. Field construction assembles compact arrays with
  global reductions, while cutoff/density use many `MPI_GLOBAL_COMMUNICATOR` reductions,
  broadcasts, barriers, and rank-0 writers. Each snapshot work parameter therefore has
  one location, and the inner dynamic location scheduler has only one item to award.
- Recommended implementation boundary: dynamically schedule epoch/location indices on
  the world communicator, but execute each awarded epoch inside an independent
  rank-local (`MPI_COMM_SELF`) or worker-group communicator with a complete mesh/compact
  field and uniquely suffixed output. No product or field collective may remain on the
  world communicator inside that worker transaction. Gather only final status and
  artifact inventory globally. The lower-risk near-term equivalent is an external pool
  of independent `mpiexec -n 1` jobs, each with eight OpenMP workers and one epoch deck.
- Resource/validation constraints: four concurrent rank-local jobs would use roughly
  the same aggregate memory already observed (~4.2--4.6 GiB each) and 32 trajectory
  threads; field initialization must avoid oversubscribing the 36-CPU allocation.
  Require exact artifact counts, per-epoch snapshot fingerprints, unchanged unresolved
  gates, and serial-versus-epoch-parallel product equivalence.
- Active run state: the existing collective restart was not altered. RBSPA snapshots
  1 and 2 have completed; snapshot 3/59 density is running in session 14675.

### 2026-10-09T07:33:24Z — Scheduling policy corrected by output mode

- Scientific/runtime decision: epoch/location ownership by one MPI rank applies to
  timestamped `OUTPUT_MODE=TRAJECTORY`, where different samples can require different
  immutable background-field snapshots. Each rank should dynamically claim one such
  epoch/location transaction, build and retain its private gridded field, and use its
  rank-local threads to divide that location's trajectory work.
- `OUTPUT_MODE=POINTS` and the locations within one `OUTPUT_MODE=SHELLS` snapshot share
  one epoch and one read-only background field. They must retain finer global load
  balancing instead of reserving a complete location for one MPI rank. If SHELLS is
  evaluated with a snapshot list, this statement applies independently inside each
  epoch: change the field between snapshots, then distribute that snapshot's complete
  location/task set across ranks.
- Existing-code audit: Mode3D cutoff already flattens independent cutoff,
  directional-map, or rigidity-list products into global task IDs and distributes them
  through the MPI RMA counter. Mode3D density still gives the MPI scheduler whole
  locations; OpenMP then divides the energy work and its direction traces within the
  owning rank. The density path therefore needs a separate shared-epoch task scheduler
  for POINTS/SHELLS, while the new outer epoch scheduler is specific to TRAJECTORY.
- Required validation: preserve unchanged physical/unresolved gates and prove
  one-rank versus multi-rank equivalence for each mode, including POINTS/SHELLS task
  reduction, TRAJECTORY epoch-product identity, unique output ownership, and dynamic
  scheduling determinism. The active SA01 run was not changed or interrupted.

### 2026-10-09T07:31:09Z — Corrected closure run stopped at snapshot 3/59

- Correction to the previously recorded active state: the RBSPA GRIDDED run stopped
  normally through the runner with AMPS return code 1 at snapshot 3/59; no process is
  active. Snapshots 1 and 2 completed, while snapshot 3 wrote cutoff, termination, and
  manifest artifacts but no density/spectrum/flux products.
- The unchanged fatal gate worked as intended. At `2014-01-06T11:00:00Z`, one of 96
  directions remained unresolved after the configured bounded extensions, producing
  `0.010416666666666666 > DS_UNRESOLVED_TOL=0.01`. AMPS reported nonzero error code
  `1794.1`, MPI aborted with code 1, and the runner recorded `status=ERROR`.
- Evidence: `test_output/SA01/quick_gridded_closure_20261009/RBSPA/amps.log`,
  `run_status.json`, and root `active_run.json`. This failure is retained and is not a
  reason to weaken the gate. The next implementation stage is the mode-aware scheduler;
  any later closure refinement must remain separately justified and validated.

### 2026-10-09T07:50:00Z — Mode-aware MPI/thread scheduler implemented in source

- `Mode3D::Run()` now selects independent-domain AMPS initialization specifically for
  `TRAJECTORY + SNAPSHOT_LIST`. Every world rank owns a complete mesh and singleton
  physics communicator; an `MPI_COMM_WORLD` RMA queue assigns one epoch/location per
  atomic fetch. Rank-local field construction, cutoff, density, writers, and frozen
  snapshot checks cannot mix epochs. Unique existing snapshot suffixes make output
  ownership independent of completion order.
- Legacy initialization files are protected by one-time world-rank serialization before
  epoch concurrency begins. A world barrier prevents the RMA window from being freed
  while another rank might still fetch; final status is reduced only after collective
  scheduler teardown.
- `POINTS` and per-snapshot `SHELLS` retain the shared distributed field. Mode3D density
  now schedules flattened `(location,energy)` work IDs across MPI ranks and threads,
  reduces disjoint transmissivity/bounds/termination rows, and only then integrates
  density, total/planar flux, channels, and detector responses from complete curves.
  Cutoff retains its already finer flattened task scheduler.
- Detailed design/safety comments were added to `Mode3D.cpp`, `DensityMode3D.cpp`, and
  `Mode3DParallel.{h,cpp}`. User documentation was updated in `srcEarth/3d/README.md`,
  `srcEarth/README.md`, and `srcEarth/test/SA01/README.md`. No physical grids, trace
  limits, classification rules, unresolved gates, or acceptance thresholds changed.
- Status / exact next action: SOURCE/DOCUMENTATION IMPLEMENTED; UNBUILT. Add source
  contract regression checks, run the focused Python tests, then perform the mandated
  clean repository build because no AMPS/SA01/build process is active.
- Timestamp correction recorded at `2026-10-09T07:40:15Z`: the `07:50:00Z` heading
  above was a transcription error; this stage was recorded at `07:40:15Z`.

### 2026-10-09T07:42:09Z — First clean scheduler build failed; scoped edit corrected

- The mandatory clean build reached `srcEarth/3d/Mode3D.cpp` and failed because the
  new world-rank log filter had been applied to `ConfigureMeshResolutionProfile()`
  instead of the later snapshot-count banner. The compiler correctly reported that
  `trajectoryEpochParallel` and `worldRank` were out of scope. No executable was linked.
- The condition was restored to `PIC::ThisThread==0` in the mesh-profile helper and
  applied to the intended banner inside `Mode3D::Run()`. Full failed-build evidence is
  retained in `test_output/SA01/mode_aware_scheduler_20261009/build.log`.
- Status / exact next action: rebuild fix applied; repeat the required verified removal
  of only `/data/vtenishe/GEO/AMPS/build` followed by repository-root `make -j`.

### 2026-10-09T07:43:55Z — Clean scheduler rebuild passed

- Reverified that no AMPS, SA01, or build process was active and that
  `/data/vtenishe/GEO/AMPS/build` was the exact, non-symlink build directory. Removed
  only that directory and repeated the mandatory repository-root `make -j`.
- The complete compile, archive, and link finished with exit code 0. The retained log is
  `test_output/SA01/mode_aware_scheduler_20261009/build_retry.log`; emitted diagnostics
  are pre-existing compiler path/format warnings and Open MPI interface-probe warnings,
  not build failures.
- Executable evidence: `/data/vtenishe/GEO/AMPS/amps`, 6,051,264 bytes, SHA-256
  `92ab356f6e74a7506fbb021d67f5452b9c45ebfb3d7218c9fe16d30f1ba3b8d1`.
  Dynamic linkage includes both `libmpi`/`libmpi_cxx` and `libgomp`, confirming the
  rebuilt binary contains the required MPI and OpenMP runtimes.
- Status / exact next action: BUILD PASSED. Re-run focused source/contract tests, then
  exercise the new multi-rank TRAJECTORY epoch queue and POINTS/SHELLS shared-epoch
  density task path with bounded integration fixtures before considering a long SA01
  restart. No physical or validation gate has been changed.

### 2026-10-09T07:45:00Z — Focused SA01 scheduler contracts passed after rebuild

- Executed `python3 -m unittest discover -s srcEarth/test/SA01/tests -v` against the
  rebuilt source tree. All 13 tests passed.
- The suite includes explicit contracts for (1) singleton physics communicators plus a
  world-communicator epoch queue in timestamped TRAJECTORY mode and (2) flattened
  `(location,energy)` density work only in shared-epoch POINTS/SHELLS modes, followed by
  post-reduction integrated-product evaluation.
- These are source/runner and numerical failure-handling checks, not a substitute for
  multi-rank physical execution. Status / exact next action: prepare bounded copies of
  existing physical fixtures and exercise both scheduling branches with the unchanged
  scientific gates.

### 2026-10-09T07:45:30Z — First two-epoch integration launch blocked by sandbox networking

- Prepared `test_output/SA01/scheduler_smoke_20261009/trajectory` from the real GRIDDED
  SA01 RBSPA deck, restricted only to the already-resolved 09:00 and 10:00 UTC epochs.
  The physics grids, trace limits, direction count, and `DS_UNRESOLVED_TOL=0.01` remain
  unchanged; this is a scheduler integration check, not a replacement SA01 study.
- The first `mpiexec -n 2` attempt exited before AMPS startup because the restricted
  execution environment denied Open MPI access to every out-of-band network interface.
  Evidence is retained in the fixture's `amps.log`; no physics ran and no code defect is
  inferred from this launcher failure.
- Status / exact next action: relaunch the identical command with the already-authorized
  unsandboxed MPI execution path, keeping child output in the log.

### 2026-10-09T09:26:22Z — Two-epoch fixture coverage mismatch rejected before physics

- The authorized two-rank launch entered independent-domain initialization and built a
  complete cached mesh for each rank, proving that both processes can initialize the
  new layout. It then stopped at `Mode3DValidateSnapshotListCoverage` before the epoch
  queue or any product calculation.
- Cause: the smoke fixture reduced `snapshot_epochs.txt` to 09:00 and 10:00 UTC but had
  copied the original 59-row trajectory unchanged. The strict validator correctly
  rejected trajectory row 3 (11:00 UTC), which matched zero requested epochs. This is a
  fixture error and does not weaken or implicate the physical/scheduler implementation.
- Evidence is retained in
  `test_output/SA01/scheduler_smoke_20261009/trajectory/amps.log`. The shell observed an
  MPI return of zero because the legacy AMPS fatal path calls `MPI_Abort(...,0)`; the
  explicit fatal log and missing products determine failure here.
- Status / exact next action: trim only the copied smoke trajectory to the same two
  epochs, retain the exact SA01 physics controls and unresolved gate, and repeat using
  the generated mesh cache.

### 2026-10-09T09:30:40Z — Epoch queue exercised; strict physical gate stopped the batch

- The corrected two-rank GRIDDED fixture entered the new world queue. Rank 1 claimed
  09:00 UTC while rank 0 claimed 10:00 UTC, and both prepared separate complete compact
  fields concurrently. Both cutoff calculations completed with status zero; their
  SHA-256 hashes exactly match the corresponding earlier collective-run artifacts
  (`525f302a...` and `7c32fe80...`), confirming epoch identity and unique output ownership.
- Each rank then entered one-location density with `OPENMP`, eight requested workers,
  25 energy points, and 96 directions per energy. The 10:00 transaction completed all
  2,400 traces but retained one unresolved direction at each of its first two energies:
  `1/96 = 0.010416666666666666`, just above the unchanged `0.01` gate. AMPS therefore
  aborted with nonzero code 1 before 09:00 could finish. No tolerance, trace limit,
  direction grid, or acceptance rule was relaxed.
- This result proves concurrent epoch assignment/field ownership and thread entry, but
  it is intentionally recorded as a failed physical batch, not a completed integration
  PASS. The density termination table differs from the earlier executable because the
  current bounded unresolved-only extension implementation replaced the old generic
  retry; the identical cutoff hashes rule out an epoch/field-selection mix-up.
- Evidence: `test_output/SA01/scheduler_smoke_20261009/trajectory/amps.log` (SHA-256
  `ce77d851...`) and the 10:00 termination summary (SHA-256 `0f9579c9...`). The prior
  coverage-failure log was preserved as `amps_coverage_failure.log`.
- Status / exact next action: use a separately labelled, inexpensive near-boundary
  scheduler fixture with the strict unresolved gate intact to verify clean queue
  teardown, then validate shared-epoch POINTS density by comparing one-rank and
  two-rank products. Do not present either smoke fixture as SA01 observational evidence.

### 2026-10-09T09:33:21Z — TRAJECTORY epoch queue and threaded density integration passed

- Ran the separately labelled bounded fixture at
  `test_output/SA01/scheduler_smoke_20261009/trajectory_near_boundary` with two GRIDDED
  T05 epochs, two MPI ranks, and eight OpenMP workers per rank. The fixture changes only
  workload size/location to keep the integration test short; it retains the fatal
  `DS_UNRESOLVED_TOL=0.01` rule and is not SA01 observational evidence.
- The world queue assigned 09:00 UTC to rank 1 and 10:00 UTC to rank 0 concurrently.
  Each transaction reported `MPI work unit: rank-owned epoch/location; threads divide
  energy/directions`, eight energy nodes, eight directions per node, and eight OpenMP
  workers. Both termination summaries reported maximum unresolved fraction zero.
- Both ranks completed density/flux with status zero, emitted their epoch-complete
  records, safely left the RMA fetch loop, collectively destroyed the scheduler, and
  returned process status zero. Exactly two uniquely suffixed density, spectrum, flux,
  termination, and manifest artifacts exist. Log SHA-256 is
  `05ae7be3714d4ae8bf9f474ef08dac95f3ad1a7ed58b92c6ee170e38e57f0f46`.
- Status / exact next action: TRAJECTORY scheduler integration PASS. Run the prepared
  shared-epoch POINTS density fixture once with one rank and once with two ranks, then
  compare numeric artifacts and confirm the two-rank log used flattened
  `(location,energy)` MPI work.

### 2026-10-09T09:36:09Z — Shared-epoch fine density scheduling passed exact equivalence

- Ran the same bounded GRIDDED T05 `OUTPUT_MODE=POINTS` fixture first with one MPI rank
  and then with two, using `DENSITY_PARALLEL=OPENMP`, eight requested workers, and the
  dynamic MPI scheduler. Both runs retained `DS_FAIL_ON_UNRESOLVED=T` and
  `DS_UNRESOLVED_TOL=0.01`; both reported maximum unresolved fraction zero and exited 0.
- The one-rank reference used the historical whole-location path. The two-rank log
  explicitly reported `flattened (location,energy) shared-epoch task`, dynamic chunk 1,
  two locations, four energy nodes, and all 32 direction traces completed.
- The density, total/channel flux, spectrum, and termination files are byte-identical
  between one and two ranks. Their respective SHA-256 hashes are `12082638...`,
  `f026d0ea...`, `edec4122...`, and `bdb68f99...`. This proves that post-reduction
  integration reconstructs the same complete access curves as the single-rank path.
- Evidence:
  `test_output/SA01/scheduler_smoke_20261009/{points_rank1,points_rank2}`.

### 2026-10-09T09:36:33Z — Final scheduler verification and handoff

- Re-ran `UTrajectoryCore` (U-F10 through U-F13), `UBoundaryProducts`, and
  `UFluxNumerics`; all passed. Re-ran all 13 focused SA01 tests; all passed.
  `git diff --check` also passed.
- No AMPS, MPI launcher, or build process remains active. The current-state section and
  stage table above were updated to replace stale process/build information.
- Implementation result: PASS for build, unit contracts, concurrent TRAJECTORY epoch
  ownership/thread subdivision, safe RMA teardown, and shared-epoch POINTS exact
  one-rank/two-rank equivalence. Full SA01 observational validation remains incomplete;
  its unchanged unresolved gate still rejects current difficult directions, and no
  observational comparison or figure is claimed from the scheduler smoke fixtures.
- Exact next action: diagnose the current-executable 10:00 and 11:00 unresolved
  directions with convergence evidence, without altering the fixed `0.01` gate, before
  launching another full GRIDDED QUICK output root.

### 2026-10-09T09:41:27Z — Detailed next-session handout added

- Added `srcEarth/handout.md` as the next-session entry point. It documents the current
  build/executable identity, implemented communicator and work-unit design, precise
  thread subdivision, successful TRAJECTORY and POINTS integration evidence, preserved
  build/fixture/MPI/physical failures, memory and initialization costs, absence of a
  completed observational comparison or plot set, mandatory rebuild procedure, and a
  gate-preserving sequence for diagnosing the remaining 10:00/11:00 unresolved rows.
- The handout links the authoritative progress journal and all relevant SA01, Mode3D,
  gridless, boundary-spectrum, C19, and repository-test documentation. It explicitly
  states that it does not replace this append-only progress record.
- No source behavior, validation threshold, executable, generated physical product, or
  process state changed in this documentation-only stage. No AMPS/build process is
  active. Exact next action remains targeted convergence diagnosis before any new full
  GRIDDED QUICK launch.

### 2026-10-09T10:13Z — Unresolved-direction investigation started

- Verified HEAD and executable SHA-256 against the handoff; preserved the dirty tree.
- The 10:00 termination file has one TIME_LIMIT at each of 18.5 and
  19.996638408877409 MeV, with zero STEP_LIMIT or invalid-field/numerical outcomes.
  Runtime logs identify RK4. No cause is inferred from the aggregate counts alone.
- Added opt-in `AMPS_DENSITY_TRACE_DIAGNOSTICS=1` records for primary unresolved
  directions and every extension: exact seed, index, rigidity, budgets, termination,
  trace time/steps, mirror/bounce/drift counts, and momentum spread. This changes
  diagnostics only. Mandatory clean build and bounded reproduction are next.

### 2026-10-09T10:22Z — RK4 invariant drift measured; controlled alternatives running

- Mandatory verified remove-build then `make -j` succeeded; log:
  `test_output/SA01/closure_diagnosis_20261009/build.log`. Previous executable
  preserved as `amps.before`. UTrajectoryCore and all 13 SA01 tests passed.
- `rk4_baseline` reproduced the 11:00 failure with 1/96 unresolved. Direction 24,
  rigidity 0.18723866390774777 GV, ran 2400 s/492468 steps, 1218 bounce cycles,
  and 94 drift revolutions, with momentum spread 0.028276386470777262. This
  exceeds the unchanged trap momentum tolerance 1e-4. It is a numerical-invariant
  violation, not evidence supporting a relaxed trap gate. The 10:00 epoch resolved
  in this repeat; its older failure remains preserved and no exact reproduction
  is claimed for that row.
- Running matched Boris and RK4 DT_TRACE=0.001 experiments with the same epochs,
  grids, field mesh, trace/extension limits and strict unresolved gate. Each has a
  diagnostic_status.json containing launcher PID, command, times and input hashes.
- SA01 workflow review also identified missing channel-integral convergence checks,
  missing prepared-epoch/hash validation in comparisons, and Matplotlib 3.1 plotting
  incompatibilities. Focused implementation repairs are in progress while physics runs.

### 2026-10-09T10:31Z — Boris closure and workflow completion repairs

- Matched `boris_reference` completed both 10:00 and 11:00 with all products and
  zero unresolved; exit 0. This is bounded diagnostic evidence, not observational
  validation. A matched Boris DT_TRACE=0.001 refinement is active.
- RK4 DT_TRACE=0.001 still violates the unchanged 1e-4 invariant tolerance:
  the original 11:00 direction 24 has momentum spread 0.003238605 after 852.85 s
  and exits the outer boundary. Reduced dt alone does not establish a trustworthy
  classification; no agreement or convergence claim is made from it.
- Added explicit recorded mover/timestep/cadence controls, input/output/epoch-bound
  comparison provenance, seven-channel convergence and physical-input identity
  checks. Fixed one-spacecraft and Matplotlib 3.1 plots, gap breaks, physical-zero
  annotations, and report separation of observational/numerical qualification.
- Found the production 1152-direction hard ceiling incompatible with PUBLICATION's
  requested 1536. Added positive parser-backed parent-grid dimensions in both
  production density paths, retaining QUICK/STANDARD defaults and archived default
  product fingerprints; non-default grids extend the control fingerprint.
- All 18 focused SA01 tests passed. USWMFCoupledProducts numerical, comparator and
  source contracts passed, including new angular-control identity/rejection tests.
  Source additions still require the mandatory rebuild after diagnostics finish.

### 2026-10-09T10:38Z — Convergence prevents a false numerical success

- RK4 DT_TRACE=0.001 exited 1 at the unchanged strict gate: 10:00 direction 95
  at 19.996638 MeV ended at 2400 s/2452456 steps with momentum spread 0.0060653.
- Boris DT_TRACE=0.001 completed with zero unresolved and ~1e-14 momentum drift.
  However comparison to Boris DT_TRACE=0.1 FAILS: total flux changes 24.36%/45.36%,
  density 27.37%/49.41%, maximum channel change 79.92%. Evidence:
  `closure_diagnosis_20261009/boris_timestep_comparison.json`. This two-epoch test
  cannot substitute for the required >=12-epoch campaign convergence population.
- Code inspection explains why momentum preservation alone is insufficient:
  historical BORIS samples B before the spatial half drift. Existing HC4 uses a
  symmetric midpoint magnetic map and fourth-order composition. No mover kernel
  was changed. Matched HC4 dt=0.1/0.001 experiments are now active.
- Verified no AMPS/MPI/build process before mandatory clean rebuild; build succeeded
  10:35:25Z. New executable hash and log are in current state above.
- Actual Mode3D products revealed UTC in SNAPSHOT_EPOCH_UTC metadata, not TimeUTC
  columns/zone names. Added validated single-location metadata handling, conflicting
  epoch and ambiguous-row rejection, and point-table declared-row count checks.
  All 20 SA01 tests passed. Live diagnostic products now parse successfully.
- Added source/build-file fingerprints, prepared-input rechecks, streamed physical
  logging, and report refresh for standalone compare/plot modes. Observation data
  and all acceptance requirements remain unchanged.

### 2026-10-09T10:46Z — HC4 baseline closure and complete plotting smoke

- `hc4_reference` completed all 4800 primary direction traces, both epochs and all
  products, exit 0, maximum unresolved 0. The original 11:00/18.5-MeV direction 24
  classified MAGNETICALLY_TRAPPED_FORBIDDEN at 48.1301 s/8828 steps, 22 bounce
  cycles, momentum spread 1.4435e-14. The matched dt=0.001 experiment remains
  active; no HC4 convergence PASS is claimed before its products exist.
- All 21 SA01 regression tests passed, including all eleven plotting families.
  Persistent synthetic PNG/EPS smoke evidence is under
  `closure_diagnosis_20261009/plot_smoke/synthetic_output/figures` (22 files,
  largest 2,328,266 bytes). These are manufactured fixtures, NOT observation plots.
- Visually inspected representative spectra and channel time series; repaired
  title/legend overlap, crowded energy tick labels, isolated samples across gaps,
  zero annotations and explicit EPS opaque rendering. Added units to scatter axes.
- Comparison now verifies declared intensity units, per-particle energy basis,
  static-magnetic mapping, effective profile grid sizes and sampled/resolved count
  closure. Runner streams output and records dependencies/source hashes before
  physical launch, including failed campaigns. A failed spacecraft branch does
  not prevent the other branch from executing; explicit interruption still stops.

### 2026-10-09T10:52Z — HC4 timestep dependence remains visible

- Both HC4 dt=0.1 and dt=0.001 completed with zero unresolved. The integrated
  comparison FAILS: total-flux changes are 1.736%/3.978%, density 1.829%/4.363%,
  and worst channel 17.475% (11:00 REPT_P2). Do not call this a qualified fix.
- Added `analyze_diagnostics.py` for reproducible bounded comparisons: all nine
  integrated quantities, product hashes, original run controls, and final direction
  classification transitions. Evidence: `hc4_timestep_comparison.json` in the
  diagnosis root. All 4800 exact direction seeds are represented at both timesteps;
  classification differences are not hidden by the zero-unresolved result.
- HC4 dt=0.002 and dt=0.0005 are running to check the smaller-step neighborhood
  and nonmonotone classification behavior, using the same source deck/population.
- All 22 SA01 tests pass, including native-product-versus-reconstruction and
  declared-unit rejection. C19 portable checks initially failed because Tk could
  not connect to DISPLAY=:1; rerunning with MPLBACKEND=Agg, no C19 source change.
- Asked the user asynchronously about a wall-clock budget for the multi-hour full
  campaign; bounded numerical investigation continues without waiting for a reply.

### 2026-10-09T18:00Z — Reproducibility correction and field audits

- Correction to the earlier timestep interpretations: `hc4_reproducibility` and
  `hc4_reference` used identical settings and executable, yet their maximum channel
  difference was 20.3376%. All earlier timestep differences are therefore confounded.
  Evidence: `closure_diagnosis_20261009/hc4_reproducibility_comparison.json`.
- HC4 dt=0.002 and dt=0.0005 both completed with zero unresolved. No job was killed.
  Zero unresolved does not establish either repeatability or numerical convergence.
- Low-level compiled T05+IGRF serial/concurrent probe: 60000 exact comparisons,
  zero differences. This does not prove the production mesh/wrapper path is safe.
- C19 portable checks passed with MPLBACKEND=Agg; SA01 22 checks and the trajectory,
  flux, boundary-product, and SWMF-coupled-product checks passed. These are regression
  results, not observational qualification.
- Mandatory clean build completed 17:55:51Z, exit 0. Added read-only sampled owner
  versus serial field audit, concurrent frozen interpolation audit, and full packed
  field content hash. These probes diagnose discrepancies; they do not repair values.
- Verified no AMPS/MPI/build job before starting `hc4_field_audit_a` and `_b`, the
  exact two-epoch repeats listed above. No data, trajectory selection, or gate changed.
- User authorized continued full-campaign execution and requested detailed comments
  in all new code and relevant READMEs. Full runs follow evidence-backed reproducibility
  repair and controlled refinement; scientific failure remains a valid test outcome.

### 2026-10-09T18:08Z — Field divergence localized; alias regression implemented

- Both instrumented repeats completed all products, exit 0. Exact comparison FAIL:
  `closure_diagnosis_20261009/pre_fix_field_repeatability.json`, maximum integrated
  channel difference 0.119893. Both epochs have differing packed-field fingerprints
  despite identical driver snapshot IDs and geometry. Each interpolation audit
  compared 12288 evaluations without a discrepancy. Owner-vs-serial audits had
  85750–93521 differing samples out of 222720; maximum relative error ~6.25e-8.
- Code evidence: meshAMRgeneric.h shares same-level center-node pointers between
  neighbor interior/ghost slots. The initializer partitions ghost-inclusive logical
  slots and writes their shared buffers concurrently; local coordinate arithmetic
  is not guaranteed bit-identical across aliases. This is a demonstrated write race,
  independent of whether it explains every trajectory discrepancy.
- Added production UniqueFieldWriterPlan: select one physical-buffer writer before
  threading, prefer owner-interior coordinates, retain deterministic ghost-only
  initialization, reject duplicate interior ownership. Serial/threaded paths consume
  the same immutable plan. No physics/gates changed; build awaits serial control exit.
- Permanent C++ regression reproduces old schedule-dependent alias writes, verifies
  fixed owner/ghost/null handling under reversed completion, rejects duplicate
  interior ownership. SA01 suite now 26 tests PASS. Detailed comments and SA01/3d
  README describe ownership, audit limits, and exact-repeat versus convergence gates.
- Serial field control is still running using the preserved pre-fix executable;
  it varies construction only, retaining two ranks/eight trajectory workers.

### 2026-10-09T18:13Z — Permanent validation integration and evidence hardening

- SA01 now has 28 passing regression tests. `SA01/list` registers the suite with
  the repository test runner without changing the hash-frozen `test/list` or its
  approved digest. Actual runner execution PASS, 15.3 s, no mismatches/timeouts;
  reports/logs: `closure_diagnosis_20261009/routine_regressions/`. Deferred provenance
  is saved there; no list last-pass entries or git commits were written.
- Re-executed UTrajectoryCore, UFluxNumerics, UBoundaryProducts,
  USWMFCoupledProducts and C19 portable checks (MPLBACKEND=Agg): all PASS.
- Exact-repeat evidence now requires independent direction coverage from each
  native spectrum's sampled counts. Energy identity alone permits eight machine
  epsilons because the logger converts rigidity back to energy; verified actual
  differences were a few ulps. Products/trajectory records remain exact. Negative
  tests reject truncated logs and genuinely changed energy. Current evidence:
  `pre_fix_field_repeatability_v3.json` (all 4800 directions; 3008 changed records).
- New diagnostics record referenced-input/output/log hashes, reject nonfinite dt
  and duplicate trajectory rows, and isolate their launcher process group for
  interruption cleanup. Reports reject stale convergence by revalidating both
  populations and the status hashes recorded in the convergence comparison.
- Detailed docstrings/comments added to new preparation, comparison, reporting,
  plotting, audit and ownership code. Relevant SA01, 3d, application and test
  READMEs updated. `git diff --check` passes.
- Preserved pre-fix executable as `closure_diagnosis_20261009/amps.field_alias_before`.
  Serial control remains active; no rebuild while it uses repository products.

### 2026-10-09T18:49Z — Serial control complete; clean fix build underway

- Serial field control completed 18:24:12Z, exit 0, both epochs/all products.
  Exact comparison with pre-fix parallel audit A FAILS: maximum channel difference
  28.8276%, different field hashes. Evidence: `pre_fix_serial_field_comparison.json`.
  Its owner audit has 89308/97840 discrepant sampled cells; concurrent interpolation
  again has zero discrepancies at both epochs. These are pre-fix observations.
- Process inspection resumed after host-tool approval and confirmed no AMPS/MPI
  job. Mandatory verified remove-build then make -j started, session 31822.
- UStep12Release numeric and source-contract suites PASS, including unchanged
  normalized test/list digest. No release gate/list checksum was relaxed.

### 2026-10-09T18:52Z — Corrected executable built; physical repeat matrix active

- Mandatory clean rebuild passed 18:51:11Z, log and executable SHA-256 in current
  state. Source changes are compiled; no acceptance setting was modified.
- Started two exact HC4 dt=0.1 repeats, two ranks/eight workers, unchanged original
  10:00/11:00 epochs and complete 25x96 populations. Also started a control changing
  only field initialization to SERIAL. New diagnostics retain all input/output/log
  hashes. Field agreement, final trajectory agreement and strict closure remain
  separate checks; none is assumed from successful compilation.

### 2026-10-09T18:54Z — Corrected fields match; trajectories still running

- Corrected initializer resolves 27840000 logical slots to 7174592 unique physical
  buffers per rank. Both repeats have zero owner-audit discrepancies out of 222720
  samples per epoch, and zero interpolation discrepancies out of 12288 comparisons.
- Full compact field fingerprints match exactly across repeats:
  10:00 `standalone-field-v1-9189195e40171081`,
  11:00 `standalone-field-v1-b3ef6ca975a07722`;
  geometry `mesh-geometry-v1-abfb6f6f56316fe8` unchanged from pre-fix runs.
  This is field-stage evidence, not yet a complete trajectory-repeatability PASS.
- SA01 suite now 29 PASS. Full runner/report share qualification logic: missing
  convergence, partial spacecraft coverage or ERROR cannot produce full PASS/exit 0.
  Explicit execution-only/comparison-only modes remain separately labelled.

### 2026-10-09T18:56Z — Exact physical repeatability PASS after ownership fix

- Both corrected repeats completed exit 0 with zero unresolved. Strict checker PASS:
  `closure_diagnosis_20261009/post_fix_exact_repeatability.json`.
  All 4800 final per-direction records match exactly (zero changed terminations,
  steps, invariant diagnostics or initial states); all 18 integrated comparisons
  (nine quantities x two epochs) have fractional difference exactly zero.
  Full packed-field identities and both sampled audits agree/are clean.
- This is bounded physical repeatability evidence, NOT >=12-epoch numerical
  convergence or observational qualification. The original nonrepeatability is
  reproduced before the fix and absent after it, without tolerance changes.
- Started same-settings one-rank/eight-worker control, two-rank/eight-worker HC4
  dt=0.001 refinement and two-rank/two-worker original RK4 retest. Corrected serial
  field control is still active. Counts/layouts are recorded in each command;
  the serial-field source path omits the pthread initializer's affinity widening,
  but its effective CPU mask has not been measured. Do not infer numerical validity
  or resource utilization from requested thread counts alone.

### 2026-10-09T19:03Z — MPI-layout parity PASS; fixed-field refinement active

- One-rank/eight-worker run completed exit 0 and matches the two-rank/eight-worker
  run exactly: 4800 final direction records, all nine integrated products at both
  epochs, full field identities and sampled audits. Evidence:
  `closure_diagnosis_20261009/post_fix_mpi_repeatability.json` (PASS, max change 0).
- Started two-rank/one-worker check, explicitly recorded `--mpi-binding NONE`.
  Host ps showed serial-field control ranks using ~105–113% CPU while concurrent
  initializing ranks averaged ~28–32%; small-job core contention is an observed
  resource issue. Existing active runs were not altered. New optional launcher
  binding control and tests/docs make future layout experiments explicit.
- HC4 dt=0.001 uses the exact corrected baseline field fingerprints and clean
  audits. Thus this is genuinely a fixed-field integration comparison, unlike the
  pre-fix experiments. Its physical products remain pending. The diagnostic checker
  now requires field evidence for controlled numerical PASS; old unaudited runs
  are INCOMPLETE with their raw numerical-difference gate still separately visible.
- Full QUICK prepared at `test_output/SA01/quick_unique_hc4_20261009`, both spacecraft,
  GRIDDED/HC4/dt=0.1, four ranks/eight threads. Preparation PASS only; not launched.
- Updated 29-test suite passed through repository test runner, 15.8 s;
  evidence `closure_diagnosis_20261009/routine_regressions_29/`.

### 2026-10-09T19:09Z — Serial field-initialization parity PASS

- Corrected SERIAL field control completed exit 0, zero unresolved, and matches
  corrected parallel initialization exactly in all 4800 final direction records,
  both field fingerprints, sampled audits, and all 18 integrated comparisons.
  Evidence: `closure_diagnosis_20261009/post_fix_serial_field_repeatability.json`.
- Full prepared QUICK populations are RBSPA 59 and RBSPB 66 epochs; no favorable
  selection or failure omission. Remaining active bounded cases are the one-worker
  layout, HC4 dt=0.001 refinement and original RK4 control. RK4 still shows momentum
  drift on identical corrected fields; wait for its actual final strict-gate result.

### 2026-10-09T19:44Z — Repeatability matrix complete; full QUICK launched

- One-worker control completed exit 0 and exact comparison PASS:
  `post_fix_worker_repeatability.json`. All 4800 final direction records and all
  integrated quantities match baseline. The tested repeatability matrix now covers
  identical settings, one/two MPI ranks, one/eight workers and serial/parallel field
  initialization. Every comparison retains exact equality, not a relaxed tolerance.
- Corrected HC4 dt=0.001 completed with zero unresolved but genuine fixed-field
  timestep comparison FAILS, maximum 0.18159821057233794 (11:00 REPT_P1). Total flux
  differences 1.6767%/2.0126%, density 1.9502%/1.4802%. Evidence:
  `post_fix_timestep_comparison.json`. This is not >=12-epoch campaign convergence.
- Original RK4 control exited 1 at strict gate 1/96=0.010416666666666666 >0.01.
  10:00 E18.5 direction35 remained TIME_LIMIT at 2400 s/552671 steps, momentum
  spread 0.0297866270. Correct field construction does not fix RK4 invariant drift.
- After host approval, full corrected QUICK started 19:43:47Z, four ranks/six
  workers, prepared HC4 dt=0.1, unchanged gates/population. Requested six workers
  are an explicit execution override of the prepared eight-worker request; actual
  command/environment/status record the override. It is an UNQUALIFIED numerical
  baseline because the timestep gate failed, not a claimed production PASS.
- Half-step HC4 dt=0.0005 diagnosis started alongside QUICK, two ranks/two workers,
  explicit MPI binding NONE. It compares against dt=0.001 on identical corrected
  fields; no scientific setting other than timestep is changed.
- All 30 SA01 tests PASS. The preserved source list and release gate digest remain
  untouched. No active run was killed or changed during the host-approval wait.

### 2026-10-09T19:49Z — Full QUICK passes the formerly blocking epochs

- Full RBSPA QUICK completed epochs 09:00, 10:00 and 11:00, all status 0, and is
  progressing to later epochs. This includes the formerly blocking physical rows;
  it is still a partial campaign, not full observational/numerical qualification.
- UFieldProvider (10 strict contracts), UDirectionalAccess (U-F14–18 including
  production CLI/C19 reader/C8 auditor), and UStandaloneProducts numeric/source
  suites all PASS. Existing tests/gates were not modified to achieve these results.

### 2026-10-09T19:59Z — Post-processing provenance and zero/gap presentation

- Comparison and convergence summaries now fingerprint their analysis/common source.
  Plot generation records source, input comparison/context, and every figure hash in
  `plot_manifest.json`. Time-energy plots mark physical model zeros separately from
  unmarked missing cells; driver/upstream/orbit curves break across genuine gaps.
  These changes affect presentation/provenance, not selection or numerical gates.
- All 30 SA01 tests PASS after these changes (15.761 s), including actual generation
  of all 11 plot families and verification of their manifest. `git diff --check`
  passed. No production C++ changed or rebuild occurred during active physical runs.
- Full QUICK has completed 12/59 RBSPA epochs, all status 0; the half-step diagnostic
  remains active. No full campaign observational or numerical PASS is claimed.

### 2026-10-09T20:03Z — Persistent failure evidence and full STANDARD continuation

- Added explicit `--diagnostic-evidence` registration, bound to prepared campaign
  hashes. Reports revalidate comparison and native diagnostic bytes and retain
  registered checks on later invocations. A bounded FAIL overrides another axis's
  PASS; a bounded PASS never supplies full-population convergence coverage. Registered
  the actual 18.1598% timestep FAIL for current QUICK; STANDARD will register the same
  explicitly relevant dt=0.1 baseline failure during preparation.
- Added recorded sequential waits with process-start identity checks and terminal
  artifact requirements. Missing terminal evidence after process exit is an ERROR.
  Continuation does not signal/restart other jobs. Added common-population refinement
  and automatic coarse/fine report/figure generation after STANDARD, including
  explicit refinement errors and machine-readable qualification status.
- All 34 regression tests PASS (14.698 s), including persisted failures, stale
  diagnostics, process disappearance and successful atomic wait completion.
  Detailed code comments and SA01 README explain these boundaries and usage.
- Full STANDARD continuation launched at 20:02:55Z, runner PID 1385739/session 21878.
  Exact command: `test_output/SA01/standard_unique_hc4_20261009.continuation.json`.
  It is WAITING for QUICK and `hc4_unique_dt00005`; no additional physical resources
  are consumed until both finish. Then four ranks/eight workers execute RBSPA/B
  sequentially with unchanged gates and the full STANDARD profile.
- This deliberate sampling-refinement baseline holds HC4 dt=0.1 and mesh fixed;
  it does not claim timestep or mesh convergence. Known numerical failure remains
  explicit even if sampling or observational metrics later pass. No C++ rebuild.

### 2026-10-09T20:05Z — Routine runner and updated plot smoke evidence

- The normal repository runner executed the dedicated SA01 list: all 34 tests PASS,
  no reference mismatches, 19.8 s. Logs/reports and unapplied deferred last-pass
  records are under `closure_diagnosis_20261009/routine_regressions_34/`. The frozen
  master list remains unchanged; no last-pass metadata or commits were applied.
- Regenerated all 11 PNG/EPS plot families (22 files, 3.7 MiB combined) under
  `closure_diagnosis_20261009/plot_smoke_v2/synthetic_output/figures/`; visually
  inspected the time-energy map's distinct zero markers and labels. This is a
  manufactured fixture, not an observational model/reference result.
- An automatic half-step comparison observer is waiting in session 93066. It will
  validate completed `hc4_unique_dt0001` and `hc4_unique_dt00005`, then write
  `closure_diagnosis_20261009/post_fix_halfstep_comparison.json`, preserving FAIL
  and exit 2 if the fixed 2% test fails. Its first sandbox attempt could not see
  the host runner PID and exited before writing results; the identical host-approved
  command now observes runner 1379911 successfully. No physical job was altered.
- QUICK RBSPA reached 23/59 completed status-zero epochs. STANDARD continuation
  remains WAITING. Updated the top-level handoff with a dated current-state notice
  while preserving its historical scheduler evidence verbatim below the notice.

### 2026-10-09T20:06Z — Live continuation verification

- Host process inspection confirms QUICK runner 1379833/launcher 1379888, diagnostic
  runner 1379911/launcher 1379955, STANDARD waiting runner 1385739, and half-step
  comparison observer 1386624 are active. Four QUICK ranks use about 510–518% CPU
  each; two half-step ranks about 194–195% each. No build is running.
- QUICK RBSPA has completed 25/59 epochs, all status 0. Executable hash remains
  `745fcb14ef5e54ac6e14cea6db9bf454aa0bde3411f49bd42a40458002183986`;
  source whitespace check PASS. Full observational and common-population results
  remain pending; the continuation records how they will be produced automatically.
- Resume by inspecting current artifacts/process identities, not by rerunning the
  launch commands. In particular, do not duplicate STANDARD while its waiting or
  active runner exists. Inspect the actual half-step JSON when generated and both
  campaigns' `qualification_summary.json`/reports/PNG/EPS after STANDARD finishes.
  Preserve numerical FAIL; further convergence studies require explicit new roots.

### 2026-10-11T19:04Z — Explicit snapshot scheduling policy and C19 validation

- Goal: remove the unsafe implicit rule that made every trajectory snapshot-list run
  allocate a private full AMR domain per MPI rank, while preserving C19 and SA01 as
  integral validation workloads with unchanged physics, sampling, and scientific
  gates.
- Repository revision and dirty state: HEAD
  `af26f41725eb42f15f4d895ee8c916234f177425`; the pre-existing unrelated `Makefile`
  modification and generated/untracked files were preserved. Policy changes affect
  Mode3D configuration/dispatch, CLI and input parsing, C19/SA01 launch provenance,
  contract tests, and the corresponding READMEs.
- Policy implemented: `MODE3D_TRAJECTORY_SNAPSHOT_SCHEDULER=COLLECTIVE` is the safe
  default. It creates one domain distributed across the MPI communicator and processes
  snapshot epochs sequentially with all ranks. `EPOCH_PARALLEL` must be selected
  explicitly, keeps the existing private-domain execution contract, and is rejected
  for incompatible Mode3D workflows. The inner trajectory scheduler
  (`MODE3D_MPI_SCHEDULER`) remains independent of this outer snapshot policy.
- Memory safety implemented for `EPOCH_PARALLEL`: after the first private domain is
  built, each shared-memory node measures the rank-zero resident-memory increment and
  projects the remaining replicas. On Linux, startup is refused collectively if the
  projected allocation exceeds `MODE3D_EPOCH_PARALLEL_MEMORY_FRACTION` times
  `/proc/meminfo` `MemAvailable` (default 0.80). If `/proc` telemetry is unavailable,
  the code reports that the preflight is unavailable and makes no false safety claim.
- Quality boundary: C19 explicitly requests `COLLECTIVE`; SA01 preparation and launch
  explicitly request `EPOCH_PARALLEL` with fraction 0.80. No grid resolution,
  trajectory integrator, timestep, trace limit, energy/angular sampling, comparison
  threshold, or PASS/FAIL interpretation was changed.
- Verification commands and results: C19 self-tests PASS; SA01 unit/contract tests
  34/34 PASS; `UDirectionalAccess/run_test.sh` PASS including U-F14--U-F18 and C8;
  source whitespace checks PASS. The mandatory clean `make -j` completed with exit 0.
  Executable SHA-256 is
  `31ebdbb25327a81a95e324e183fcb0b6180ceeda1af94adb5dfde1d73fabb769`;
  build log is `test_output/C19_snapshot_policy_20261011/build.log`.
- Physical C19 result: repository test runner exit 0 after 4744.6 s, one expected PASS
  and one actual PASS, zero mismatches/timeouts. The collective mesh contained 330856
  blocks and completed the three T05 epochs with status 0. C19 produced 12 model rows,
  12 reference rows, zero run failures, `trajectory_resolution_passed=true`, and
  consistent plotted/direct row counts. It generated 27 PNG and 27 EPS files.
- Scientific interpretation: C19's SMOKE profile leaves provisional observational
  thresholds unenforced. `C19_result.json` therefore records execution complete and
  numerically complete but `observational_passed=false` and `passed=false`; four of
  the six reported aggregate metric groups pass their provisional thresholds. The
  repository runner PASS verifies the C19 execution/regression contract and restored
  plot production, not publication-grade observational agreement.
- Evidence: `test_output/C19_snapshot_policy_20261011/test_runner_report.json`,
  `test_output/C19_snapshot_policy_20261011/runner_logs/`, and
  `test_output/C19_snapshot_policy_20261011/C19_goes_epead_ew/C19_result.json`.
  Comparison plot:
  `test_output/C19_snapshot_policy_20261011/C19_goes_epead_ew/C19_comparison_gridded_t05.png`.
- SA01 execution state: no physical SA01 campaign was run in this stage, per the user
  instruction that SA01 requires a small additional tweak first. The earlier completed
  QUICK/STANDARD results remain preserved and scientifically unqualified.
- Exact next action: implement and verify the requested SA01 tweak, perform the
  mandated clean rebuild, then run SA01 from a new evidence root using the explicit
  `EPOCH_PARALLEL` policy and retain the memory-preflight record.

### 2026-10-11T19:29Z — independent QUICK test and STANDARD study integration

- Goal: make both SA01 profiles independently executable from their appropriate
  repository procedures, with comparison plots required from each completed run and
  without weakening any scientific comparison or convergence gate.
- Repository revision and dirty state: HEAD
  `af26f41725eb42f15f4d895ee8c916234f177425`; the pre-existing unrelated `Makefile`
  modification, snapshot-policy implementation, generated evidence, and untracked
  files were preserved. No commit was created.
- Entry points: appended the full QUICK command to `srcEarth/test/list` as an
  exclusive expected-pass test and appended the full STANDARD command to
  `srcEarth/studies/list`. Both commands use the authoritative external data root,
  explicit GRIDDED/HC4 controls, both spacecraft, PNG and EPS output, and distinct
  output parents. QUICK requests four MPI ranks and six workers; STANDARD requests
  four ranks and eight workers.
- Output isolation: added `--output-parent` to reserve a new empty child named from
  profile, UTC timestamp, and process ID. Reservation uses atomic directory creation
  and a deterministic collision suffix. Existing `--output-dir` remains available
  for controlled one-off runs, while stage-only modes deliberately reject
  `--output-parent` because they operate on an already named campaign.
- Completion semantics: added `--exit-code-policy`, whose default remains
  `QUALIFICATION`. List entries use `WORKFLOW`, which can return zero only after the
  two physical branches pass execution, comparisons are evaluated, the report is
  written, and every plot listed in the manifest exists. The emitted
  `workflow_completion.json` records evidence/figure hashes and copies the actual
  scientific qualification. It never changes comparison thresholds, convergence
  thresholds, comparison counts, or `qualification_summary.json`.
- Exclusivity: the QUICK list entry includes `--require-runner-exclusive`. The runner
  checks `AMPS_TEST_RUNNER_EXCLUSIVE=1` before reserving output or starting MPI, so a
  future parser/runner regression cannot silently run this high-memory test in a
  shared slot. Standalone QUICK and the separate STANDARD study do not depend on
  that list-specific handshake.
- Code and documentation: added detailed comments/docstrings for output reservation,
  the exclusivity handshake, and workflow evidence; documented exact commands and
  the workflow-versus-qualification boundary in `srcEarth/test/SA01/README.md` and
  `srcEarth/test/README.md`; updated the focused SA01 list comment. Advanced the
  unchanged Step-12 structural digest to the reviewed master-list content and added
  exact assertions for the new command instead of weakening the frozen-list guard.
- Verification: `python3 -B -m unittest discover -s srcEarth/test/SA01/tests -v`
  passed 38/38; `srcEarth/test/UStep12Release/run_test.sh` passed all eight release
  validation and five source-contract tests; runner `--dry-run srcEarth/test/list`
  parsed 62 tests and classified SA01 QUICK as exclusive; `run_SA01.py --help`
  passed; `git diff --check` passed.
- Real-data preparation: independent QUICK and STANDARD `--prepare-only` executions
  both passed against the authoritative dataset. QUICK retained 59 RBSPA and 66
  RBSPB epochs; STANDARD retained 251 and 261. These temporary preparation checks
  are interface evidence only and do not claim physical execution or qualification.
- Plot evidence: the prior completed QUICK and STANDARD campaign manifests under
  `test_output/SA01/{quick,standard}_unique_hc4_20261009` each contain hashes for 12
  PNG and 12 EPS comparison figures, and all 48 referenced files are present. The
  new list paths invoke the same plot stage and require its manifest/files before
  workflow completion. No new expensive physical campaign was launched in this
  implementation stage.
- Result: PASS for repository integration and non-physical verification. Scientific
  status remains unchanged: both prior campaigns fail 24/34 observational gates and
  the common-population convergence requirement. A future completed list/study run
  must expose that result separately even when its workflow completes successfully.

### 2026-10-11T19:36Z — standalone QUICK invocation clarification

- A direct QUICK command copied with the master-list-only
  `--require-runner-exclusive` option correctly refused to start because a standalone
  shell does not receive `AMPS_TEST_RUNNER_EXCLUSIVE=1` from `test_runner.py`. This
  was an operator-facing documentation defect, not a physical runner failure.
- Preserved the safety boundary: `srcEarth/test/list` still includes both the
  exclusive directive and handshake option, so the repository procedure cannot
  silently run QUICK concurrently under an obsolete or incompatible runner.
- Corrected the documented independently executable QUICK command to omit the
  private option. Expanded both READMEs to state that standalone resource exclusivity
  is the operator's responsibility, while STANDARD remains directly executable as
  written. Expanded the runtime error with the exact remediation: omit
  `--require-runner-exclusive` for a direct campaign or use the unchanged list entry
  through the current test runner.
- The user explicitly made `run_SA01.py` executable; that mode change was preserved.
  The shebang invocation and `python3` invocation remain equivalent.
- Verification after the correction: 38/38 SA01 implementation/contract tests pass;
  the complete Step-12 suite passes eight release-validation plus five source-contract
  tests; `git diff --check` remains required at final handoff. No physical campaign
  was launched or modified during this clarification.

### 2026-10-11T19:38Z — exact registered QUICK command made standalone-executable

- The preceding documentation-only correction was insufficient because copying the
  exact payload from `srcEarth/test/list` still copied its private handshake option.
  Removed that SA01 option from the registered command and from `run_SA01.py` rather
  than requiring an operator to edit a reviewed test command or forge runner state.
- Test-runner resource quality is unchanged: `! runner: exclusive` remains directly
  before QUICK, and `test_runner.py` still drains the process pool before launch and
  blocks all other launches until QUICK exits. The direct command uses the identical
  physical, sampling, comparison, reporting, and plotting configuration; only the
  redundant environment prerequisite was removed.
- Added a regression that tokenizes the exact master-list payload and parses it
  unchanged as a standalone SA01 CLI. Retained independent assertions for QUICK and
  STANDARD ownership, unique output parents, GRIDDED/HC4 controls, both plot formats,
  and WORKFLOW completion. The Step-12 source contract now freezes adjacency of the
  exclusive directive as well as the exact SA01 command and advances its reviewed
  normalized-list digest.
- Updated both READMEs to describe exclusivity as a repository scheduler property,
  not a command-line prerequisite. Historical journal entries above describe the
  intermediate handshake implementation and are intentionally retained; this entry
  supersedes that implementation.
- Verification: 38/38 SA01 tests pass, including the exact-command regression; all
  eight Step-12 release-validation and five source-contract tests pass; the full
  master-list dry run parses 62 tests and reports SA01 QUICK as exclusive with the
  independently executable command. No physical campaign was launched in this fix.
