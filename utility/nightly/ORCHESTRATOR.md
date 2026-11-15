# AMPS nightly test orchestrator — design

Status: phase 1 (sequential) implemented 2026-11-14: `amps-nightly run`, `amps_nightly/{orchestrator,workspace,
procs,readers,makefile_runner}.py`, `runners.json`, and `prune` (run by `update`). Not yet: the concurrent
scheduler (section 7).
Context: [DESIGN.md](DESIGN.md) (record format, results store, report, machine files).

## 1. What it is

The orchestrator is `amps-nightly run`. It runs on each test machine every night, and does this for
each enabled configuration (compiler) in the machine file:

```
check out the night's commit
for each runner listed in runners.json, in order:          (phase 1: one at a time)
    make a private copy of the checkout and configure it   (Config.pl calls from the machine file)
    run the runner's own build steps, then its command, under a time limit
    read the runner's own report into the common result rows
    save the runner's logs into the night's archive folder
write DONE                                                 (script 1: nightly_run.sh)
tar the run folder, send the tarball to passion            (script 2: nightly_pack.sh)
passion unpacks it, updates results.db, archive and pages  (script 3: nightly_update.sh)
```

A **runner** is an existing test program that already works on its own:

- `src/test/run_tests.py`
- `tools/sep_test_orchestrator.py` (srcSEP, srcSEP3D)
- `srcEarth/test/test_runner.py`
- `srcMoon/test/run_tests.py`
- `srcEnceladus/test/run_tests.py`
- the `Makefile.test` tests

The orchestrator never re-implements a test. It only runs the commands listed for each runner and
interprets the report that runner already writes. Adding a runner means adding an entry to
`runners.json`. If the runner's report format is new, it also means a reader function of about 50 lines.

## 2. The runners file: `utility/nightly/runners.json`

This file is the list of test-runner commands that the orchestrator executes. It replaces the `sources`
list of `routine.json`, and that file is renamed `runners.json`. The `disabled` list stays in it,
unchanged.

```json
{
  "schema_version": 1,
  "checkout": {
    "repository": "git@github.com:SWMFsoftware/AMPS.git", "branch": "master",
    "companions": [{"repository": "git@github.com:SWMFsoftware/BATL.git", "path": "BATL"},
                   {"repository": "git@github.com:SWMFsoftware/SWMF_data.git", "path": "SWMF_data", "shared": true}]
  },
  "runners": [
    {
      "id": "core",
      "title": "src/test core suites",
      "workspace": "configured",
      "prepare": [],
      "command": ["{python}", "src/test/run_tests.py", "--config", "{core_profile}", "--campaign", "nightly",
                  "--jobs", "{jobs}", "--output-dir", "{out}/run", "--json", "{out}/report.json"],
      "env": {"AMPS_TEST_TEMP_ROOT": "{tmp}"},
      "time_limit_min": 120,
      "results": {"reader": "src-test-json", "file": "{out}/report.json"},
      "families": "by-prefix",
      "logs": ["{out}/run/logs/*.log", "{out}/run/streams/*.log"],
      "requires": [], "requires_vars": ["core_profile"],
      "resources": {"cores": 16, "memory_gb": 16, "exclusive": false}
    }
  ],
  "disabled": []
}
```

| field | meaning |
|---|---|
| `id`, `title`, `enabled`, `note` | as for sources today |
| `workspace` | `configured` (default): a private worktree of the night's commit, configured with the configuration's `Config.pl` calls. `plain`: a private worktree without them, for runners that never build AMPS (Enceladus). |
| `prepare` | commands run before `command` (application build: `Config.pl -application=...`, `make amps`); their failure is stage `build` |
| `command` | the runner command (argv list, run from the workspace root) |
| `env` | extra environment variables |
| `time_limit_min` | wall-clock limit for `prepare` + `command` together |
| `results.reader`, `results.file` | how to read the runner's report (see §4) |
| `families` | `"by-prefix"` (test IDs mapped by `suite_prefixes` in families.json) or one family id |
| `jobs` | upper limit for `{jobs}` (default: the machine's `cores`) |
| `expect.src_test_campaign`, `expect.match` | the src/test suites the runner should report; listed by `plan`, and a missing one becomes FAILED "no result" |
| `logs` | globs of the runner's own log files to archive (see §6); the full console output is always archived |
| `requires`, `requires_vars` | machine capabilities and machine `vars` the runner needs; if one is missing the runner's tests are reported as SKIPPED with that reason |
| `resources`, `locks`, `after` | used only by the concurrent scheduler of phase 2 (see §7); ignored in phase 1 |

**Placeholders** in `prepare`, `command`, `env`, `results.file` and `logs`:

| placeholder | value |
|---|---|
| `{python}` | the machine's Python (`/usr/bin/python3` on passion) |
| `{ws}` | the workspace root |
| `{out}` | the runner's output directory; outside the workspace, so a runner can never delete it |
| `{tmp}` | a private temp directory, also exported as `TMPDIR` |
| `{jobs}` | cores granted to the runner; in phase 1 all of the machine's cores |
| `{mpiexec}` | the machine's MPI launcher |
| `{table_host}` | the machine's Table host name |
| `{<name>}` | any entry of the machine file's new `vars` object, e.g. `{spice_path}`, `{core_profile}` |

The machine file gains three fields: `python`, `vars` (machine-specific paths such as SPICE kernels and
data roots) and `max_parallel` (1 = phase 1). Runner entries therefore contain no machine-specific paths.

## 3. Workspaces: every runner gets its own copy of the checkout

The study of the runners showed that most of them reconfigure and rebuild AMPS in the repository root,
and none of them restores it:

- `src/test` native suites run `rm -rf build`, `Config.pl -application=test/...` and `make amps` in the root.
- The SEP orchestrator builds with `Config.pl -application=test/sep_parker_spiral__field_line`.
- srcSEP3D's build test needs the tree configured with `-application=sep3d`.
- Earth needs `test/geospace/geospace`.
- Moon needs `-application=moon`.
- The legacy tests run `Config.pl -application=<test>` for every test.

Config.pl settings also leak from one runner into the next. The Moon adapter already passes
`-amps-test=off` because a legacy run left `TESTMODE` on.

Two more constraints:

- `src/test` native suites refuse to run while any `make`, `Config.pl` or `amps*` process has its
  current directory inside the checkout.
- SEP, Earth and Moon each need their own `./amps`.

So, already in phase 1:

- The night's commit is taken once from a persistent mirror (`{workdir}/mirror/AMPS.git`, updated with
  `git fetch`), so every runner and configuration tests the same commit.
- Each runner with `workspace: configured` gets `git worktree add --detach <commit>`, the companions (BATL
  cloned; SWMF_data is one persistent shared clone, linked in), and the configuration's `Config.pl`
  calls run with `AMPS_TABLE_HOST={table_host}`.
- A configured tree cannot simply be copied, because `Makefile.def` holds absolute paths
  (`PTDIR`, `MYDIR`). Phase 1 measures how long the per-workspace `Config.pl -install` takes.
- The orchestrator's own working directory is `{workdir}`, never inside a workspace. Runners are started
  directly, not through `make`.
- Clean environment for each runner:
  - `PATH=/usr/bin:/bin` plus the configuration's `setup`; make recipes call bare `python3`, which must
    be the system Python with numpy, matplotlib and aacgmv2.
  - `HOME`, `USER` and `LANG` are kept.
  - Every inherited `AMPS_TEST_*` variable is removed.
  - `TMPDIR={tmp}`, because several runners use fixed temp-file names.
  - A marker `AMPS_NIGHTLY_TAG=<run-id>/<runner>` (see §5).
- Workspaces of runners whose tests all passed are deleted once their logs are archived. Failed ones are
  kept until the next night for debugging (`keep_workspaces: "failures"`). The checkout is about 140 MB;
  a `src/test` run directory is about 700 MB, mostly work directories that are not archived.

This isolation is also what makes phase 2 a scheduling change only (§7).

## 4. Readers: from each runner's report to the common result rows

Each reader turns one report format into rows `{family, test, status, stage, reason, duration_s,
cases_total, cases_failed, log}`.

| reader | runner, report | mapping |
|---|---|---|
| `src-test-json` | run_tests.py `--json` | one row per suite (`results[]`), case counts from `cases[]`. FAIL → FAILED; ERROR → FAILED with "harness: " before the message; SKIP → SKIPPED. Stage from the message: "did not build" → build; "timed out" → run; else check. |
| `sep-orchestrator-json` | tools/sep_test_orchestrator.py `run_all_test.json` | one row per record; `kind: check` → stage build; per-test log from `log`/`line` |
| `earth-report` | test_runner.py `<prefix>.json` | one row per list entry; name from the command and list line; reason from the last FAIL line of its log. Entries the list expects to fail stay FAILED, with reason "expected to fail (srcEarth/test/list)". |
| `moon-summary` | srcMoon `summary.json` (+ `U*/result.json`) | one row per U-test; SKIPPED keeps its reason |
| `enceladus-receipt` | srcEnceladus `test_receipt.json` | one row per result (PASS/FAIL only) |
| `junit` | any JUnit XML | one row per testcase |
| `exit-code` | anything | one row: PASS if the exit code is 0 |
| `amps-results` | a runner that writes our `results.jsonl` itself | taken as is; the recommended format for new runners |
| `makefile` | built in | the Makefile.test tests (§5.2) |

**Rules that apply to every runner**, so that a failure can never look green:

| situation | rows written |
|---|---|
| report missing or unreadable | one FAILED row `<runner>` with the last `ERROR:` line of the runner's output, or "exited with N, no report" |
| exit code non-zero but no failing row | one FAILED row `<runner>: exit status N` |
| time limit exceeded | the rows of a partial report if there is one, plus one FAILED row `<runner>: time limit of N min exceeded` |
| a test listed by `plan` has no row | that test FAILED with "no result (runner stopped early)" |
| prepare (build) failed | every planned test of the runner FAILED with stage `build` and the build log linked; one row if the list is unknown |

## 5. Running a runner (phase 1: sequential)

### 5.1 Common procedure

1. Take a lock (`flock {workdir}/lock`); if last night's run is still going, log it and exit.
2. Resolve the plan (`amps-nightly plan`). Runners that are switched off or lack a capability or `var`
   produce their SKIPPED rows without running.
3. For each runner, in file order:
   1. create the workspace;
   2. run `prepare`, then `command`, with stdout and stderr to `logs/<runner>/runner.log`;
   3. apply the reader and the rules above;
   4. archive the logs;
   5. append the rows to `results.jsonl` and the runner's line to `runners.json` right away, so a crash
      still leaves everything finished so far;
   6. clean up stray processes.

   A failing runner never stops the procedure.
4. Write `run.json` (with `complete: true` if every runner finished), then `DONE`. `run.partial.json`,
   written when the run starts, and `progress.jsonl` are removed.
5. `amps-nightly pack` (a separate step, `transfer.py`):
   - tars each finished run folder that is not `PACKED` into `{workdir}/outbox/<run_id>.tar.gz`, `DONE` last;
   - sends it to the machine file's `upload`: `rsync` to `passion:AMPSTEST/nightly/incoming`, or a local
     copy on passion; the file shows up under its final name only once it is complete;
   - an unsent tarball stays in the outbox and is retried by the next `pack`; sent ones are deleted after
     14 days;
   - a run folder without `DONE` from a night that is no longer running (the orchestrator died) is
     completed first: rows from `progress.jsonl`, run record from `run.partial.json`, a FAILED row for
     the runner that was running (found in `orchestrator.log`), and `complete: false` with a note;
   - `run --no-upload` marks its run folders `PACKED`, so trials are never sent;
   - the next `run` deletes old nights' workspaces but keeps run folders that are not packed yet.
6. `amps-nightly update` on passion unpacks each `incoming/*.tar.gz` into a hidden staging folder after
   checking every member (regular files and folders under one top folder: no absolute paths, no `..`, no
   links, at most 8 GB), then ingests it like a run folder. A bad tarball goes to `rejected/` with
   `REJECTED.txt`.

**Starting the night (`bin/amps_nightly_test.sh`).** Cron starts the copy of this script in `~/bin`
(installed with `amps-nightly install`; README "Installing the nightly on a machine"). It downloads AMPS
(`git clone --depth 1` of the branch), updates the nightly's own SWMF_data copy (reset and cleaned every night:
Makefile.test checks edit reference files in place), and runs `nightly_run.sh --commit <sha>` from the
download, so the tool, `runners.json`, the machine files and the tested code are all the same commit; then
`nightly_pack.sh` and, on the server, `nightly_update.sh`. `--commit` must be on the branch. The commits of
share, util, BATL and SWMF_data go into the run's notes, and so do the script's warnings (`--note`), such as
"the `~/bin` copy differs from the branch's version"; the night's summary page shows both.

**Process control.**

- Each command is started in its own session, with the marker variable `AMPS_NIGHTLY_TAG` set.
- At the time limit, the orchestrator sends SIGINT to the runner's process group. run_tests.py handles
  SIGINT and kills its suites, but it ignores SIGTERM.
- After 60 s it sends SIGKILL to every process of this user whose environment carries the runner's
  marker. This also reaches grandchildren that started their own sessions (run_tests adapters,
  Earth tests, mpiexec ranks).
- The same sweep runs after every runner, so no orphaned job can disturb the next one.

**Measurements kept for phase 2.** For each runner, `runners.json` records start and end times, exit
code, whether the time limit hit, and CPU time. It also records peak memory, from the summed resident
memory of the marker processes sampled every 30 s.

### 5.2 The Makefile.test runner (built in)

The test list comes from `make -s test_manifest`. For each test, the stages
`make test_<X>_compile`, `_rundir`, `_run` and `_check` are run, each with the machine's time limit for
that stage. The output goes to one log per test: `logs/legacy/tests/<X>.log`. That replaces the 96 MB
single `test_amps.log`, which is mostly about 12,000 compiler command lines.

Status and reason are decided from the `.diff` exactly as the legacy importer does now; `.diff` contents
are still not stored. Compiling changes the tree, so compiles run one after another. The run and check
stages work in their own `run_test_<X>` directories, so up to `lanes` of them run while the next test
compiles. That is the old scheduler's split, now with time limits and stage logs.

### 5.3 Proposed runner list for passion

| order | runner | workspace / prepare | command (abridged) | reader | estimated time |
|---|---|---|---|---|---|
| 1 | `core` | configured | `src/test/run_tests.py --campaign nightly` | src-test-json | ~27 min |
| 2 | `enceladus` | plain (prepare copies the git-ignored data file) | `srcEnceladus/test/run_tests.py --output {out}` | enceladus-receipt | < 1 min |
| 3 | `sep` | configured (the SEP orchestrator builds) | `tools/sep_test_orchestrator.py --app srcsep --output-dir {out}` | sep-orchestrator-json | ~15 min |
| 4 | `sep3d` | configured | `tools/sep_test_orchestrator.py --app srcsep3d --output-dir {out}` | sep-orchestrator-json | ~15 min |
| 5 | `moon` | configured; prepare: `Config.pl -application=moon …`, `make` | `srcMoon/test/run_tests.py --all --output-dir {out}` | moon-summary | ~10 min |
| 6 | `earth` | configured; prepare: `make test_Earth-cutoff-rigidity-3d-analytic_compile`, then check `./amps -h` prints the Earth banner | `srcEarth/test/test_runner.py srcEarth/test/list -j {jobs} --log-dir {out}/logs --report-prefix {out}/report --last-pass-results-file {out}/last_pass.json --timeout …` | earth-report | ~1.5–3 h |
| 7 | `legacy` | configured | built-in makefile runner | makefile | 6–7.5 h, the legacy figure with 10 lanes |

Runner order is shortest first, so most results exist even if the night runs long.

The Moon regression (4 h, 8 ranks) and the Enceladus native/E3 runs stay off. They are long
qualification runs, available through the `src/test` campaign `application-smoke`/`extended`.

**Time.** Run one after another, the total is about 9–11 h, so it does not fit 22:30–08:00. Phase 1
starts at 20:00 with an overall limit of 14 h. Phase 2 (§7) is what brings it back inside the night.

## 6. Logs and the archive

The archive has the same shape as the legacy `AMPS_TEST_RESULTS/YYYY/MM/DD/<folder>/`, one folder per
column and night:

```
~/AMPSTEST/nightly/www/archive/2026/11/14/passion-gnu/
    index.html            summary of the night: commit, runners with time, exit code and status, links to logs
    run.json  results.jsonl  runners.json
    logs/
        orchestrator.log          checkout, workspaces, Config.pl calls, timings, cleanup
        <runner>/runner.log       complete console output of prepare + command
        <runner>/tests/<test>.log per-test logs the runner wrote (copied), or cut from a shared log
        <runner>/report.json      the runner's own summary report
    DONE
```

- **Folder date:** the local date on which the night started. A second run of the same column on the
  same night goes to `<column>.2`. This replaces the current `www/runs/YYYY/MM/DD/<run-id>/`; the 628
  imported nights are moved, and their database paths updated.
- **What is archived:**
  - each runner's console output;
  - the files matched by its `logs` globs:
    - `src/test`: `logs/<ID>.log` and `streams/*.log`;
    - Earth: one log per entry;
    - SEP: the step logs and `failed-tests/**`;
    - Moon: `U*/result.json`, `stdout.log` and `stderr.log`;
    - legacy: the per-test make logs;
  - the runner's summary report.
- **Not archived:** work directories, executables, plots, artifacts, and legacy `.diff` contents.
- **Size:**
  - A legacy night's `test_amps.log` was 96 MB, which is how the legacy archive reached 47 GB; gzip
    shrinks it to 8 MB.
  - Every step of a log (each `### ` header that procs.run writes starts one: compile, rundir, run,
    check, prepare, command) keeps its first 256 KB and last 768 KB, with a line saying how much was left
    out. The full run of 2026-11-14 had 104 MB, 87 MB of it in eight Makefile.test logs: those test
    programs print 80,000–240,000 lines each (`magnetic-reconnection`: "need to inject npart" 110,000
    times; `cdap_test`: about 95,000 lines of "matvecs ... NaN"). Trimmed it is 38 MB, and its tarball 2.1 MB.
  - A log copied from a runner's own files that is over 50 MB is first cut to its first 10 MB and last 40 MB.
  - Logs stay plain text for 14 nights, so the pages link to them directly (Firefox shows `.gz` only as
    a download). After that, `amps-nightly prune` gzips them.
  - Logs are deleted after 400 days; `run.json`, `results.jsonl` and `runners.json` are kept forever.
  - Expected size: about 40 MB per column-night before compression, a few GB per column-year after.
- **Report:**
  - A failed test's "log" link opens its own log, or else its runner's log.
  - The column header links to the night's `index.html`, which plays the role of the legacy `test_amps.html`.
- **Legacy archive:** `AMPS_TEST_RESULTS` stays where it is as the read-only history before the switch.

## 7. Phase 2: running runners concurrently

Phase 1 already gives each runner a private workspace, temp directory, output directory and process
marker. Phase 2 is therefore only a scheduler in place of the loop of §5.1, switched on with
`"max_parallel"` > 1 in the machine file.

- **Resources:**
  - The machine declares `resources: {cores, memory_gb}`.
  - Each runner declares `resources: {cores, memory_gb, exclusive}`, set from the phase-1 measurements.
  - Its `{jobs}` becomes the granted cores, which it passes to its own parallelism: run_tests `--jobs`,
    Earth `-j`, `MAKEFLAGS=-j{jobs}`, the legacy lanes.
- **Ordering:**
  - Runners are started in file order whenever their cores and memory are free, their `after` runners
    have finished, and their named `locks` are free.
  - The machine's free memory is also checked before each start, as Earth's runner does.
  - An `exclusive` runner waiting at the head of the queue blocks new starts until it has run, so it
    cannot starve.
- **Locks for shared paths outside the workspaces:**
  - srcSEP3D's HARN04 writes fixed `/tmp/sep3d_harn04_test.json`, so two srcSEP3D runs (two
    configurations) must not overlap: lock `sep3d-tmp`.
  - Timing-sensitive suites, such as the PIC performance suite or Earth's exclusive C19 entry, need
    their runner to be `exclusive`.
- **Configurations:** all configurations of one machine feed the same queue, since each runner already
  has its own workspace.
- **Legacy runner:** with more than one lane, compiles also run in parallel in extra worktrees, one per lane.
- **Batch machines (pleiades):** the orchestrator runs inside one batch job; that job's allocation is
  the `resources`.

## 8. Problems found in the runners (to fix in their own code, independent of the orchestrator)

- **srcEnceladus:** `test/data/shape_manifest.json` is excluded by the `.gitignore` pattern `data`, so
  U01 and nine later cases fail in a fresh checkout. A copy exists in `~/Mars2/AMPS`.
- **srcEarth:**
  - `make test_Earth-…_compile` ignores a failed `make amps` (leading `-`); the orchestrator checks
    `./amps -h` instead.
  - SA01–SA03 pass `--amps-root /data/vtenishe/GEO/AMPS`, so their provenance names another checkout.
  - C10 asks for matplotlib ≥ 3.7, but the system has 3.1.2.
  - The list starts up to 64 threads on 36 cores.
  - No timeouts by default.
- **srcSEP:**
  - `--routine` passes the executable's exit code through and writes no JSON on a crash.
  - A plotting failure exits 2 even when the tests passed; nightly uses `--plot none`.
  - With `--timeout`, a hung MPI grandchild can block the runner.
- **srcSEP3D:**
  - The 300 s per-command default also covers cold builds.
  - A build timeout ends in a traceback (exit 1, no JSON).
  - BLDL3D01 fails unless the tree is configured for sep3d.
  - HARN04 uses a fixed `/tmp` path.
  - The tracked `test/stage1` binary is rewritten.
- **src/test:**
  - SIGTERM is not handled, so killing the process orphans its suites.
  - Running it through `make -C src/test nightly` probably makes every native suite report "native
    build/run already active", because the parent `make` has its working directory in the checkout.
    Reasoned from the code, not tested.
- **srcMoon:**
  - No timeouts anywhere.
  - U02 needs the `moon.input` that `Config.pl -application=moon` copies in.
  - `--strict-skips` always fails today, because 10 tests are not implemented.
- **Legacy:**
  - The old scheduler ran 10 `make` lanes in one tree.
  - Config.pl state carries over between tests and runners.

## 9. Implementation steps

**Phase 1:**

1. `routine.json` becomes `runners.json` with the fields of §2. `check` validates readers, placeholders
   and `vars`, and `plan` shows each runner's command with the placeholders filled in.
2. Workspace manager: mirror, worktrees, companions, configuration.
3. Process control: sessions, marker variable, time limit, SIGINT then SIGKILL sweep, memory sampling.
4. Readers of §4, and the built-in makefile runner.
5. Archive writer: layout, size rules, per-night `index.html`. Ingest and report move to `archive/`
   with log links, and `prune` is added.
6. `amps-nightly run`, with:
   - `run --only <runner> --keep` for trying one runner entry before committing it;
   - `run --dry-run`, which prints the resolved commands;
   - the cron line.
7. Trial on passion: each runner alone, then a full night; compare the legacy results with the old scripts.

**Phase 2:** the scheduler of §7, with `resources` set from the measurements of phase 1.
