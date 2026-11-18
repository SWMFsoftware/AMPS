# AMPS nightly testing — design

Status: step 1 of 4 implemented (record format, results store, report pages, legacy importers);
machine files and the routine file are defined and checked (`check`, `machines`, `plan`), the runner is next.
Started 2026-11-13. Replaces the per-machine csh/perl scripts in `utility/TestScripts`.

## Goals

- One procedure that runs every kind of AMPS test on any number of machines:
  the legacy application tests in `Makefile.test` (generated from `MakefileTest/Table`),
  the core tests of `src/test/run_tests.py`, and the applications' own test runners.
- Adding a machine means adding a configuration file, not copying scripts.
- Every result is stored with the commit it tested, so the report can say when a test last
  passed and at which commit.
- A web page: test families × columns, then one page per family with its tests.

## Decisions

| Decision | Choice (2026-11-14) |
|---|---|
| Grouping of legacy tests | a `Family=<id>` line in every `MakefileTest/Table` block |
| Report columns | one per **machine + compiler configuration** (e.g. `passion · GNU`) |
| Server and web host | **passion** (`/data/vtenishe/AMPSTEST`; tests run in `~/Tmp_AMPS_test`) |
| Implementation language | Python 3, standard library only, 3.8 compatible (cron on passion runs `/usr/bin/python3` 3.8) |
| Test statuses | `PASS`, `FAILED`, `SKIPPED`; harness errors and timeouts are `FAILED` with a reason |
| Machines | one file per machine, `machines/<machine>.json`; its enabled configurations are the report columns |
| What runs | `routine.json` lists the test sources; tests are switched off there (`disabled`, shown as SKIPPED with the reason), not by editing the sources |

## Concepts

- **Column** — one machine + compiler configuration, id `<machine>-<config>` (`passion-gnu`).
- **Run** — the results of one column for one night. Each run records the commit it tested.
- **Test** — the smallest item with a status: a `Makefile.test` test, a `src/test` suite
  (with its case counts), or an application test.
- **Family** — a named group of tests, defined in `families.json`; families are grouped by source
  (`legacy`, `core`, `app`, `other`).
- **Family status** in a column: `FAILED` if any test failed; `SKIPPED` if every test was skipped;
  otherwise `PASS`. Counts are of individual tests; a `src/test` suite counts its cases.

## Data flow

```
 each test machine (cron, or a PBS/Slurm job)          server: passion
 ───────────────────────────────────────────           ──────────────────────────────────────────
 1 nightly_run.sh  (amps-nightly run)                  /data/vtenishe/AMPSTEST/
   checkout → run tests → run directory                   incoming/<run-id>.tar.gz
     run.json, results.jsonl, logs/, DONE    2 pack      3 nightly_update.sh (amps-nightly update):
 2 nightly_pack.sh (amps-nightly pack)     ───────────►   ingest  → results.db (SQLite)
     outbox/<run-id>.tar.gz, rsync/copy                           → www/archive/YYYY/MM/DD/<column>/
                                                          amps-nightly report  → www/index.html
                                                                               → www/family/<family>/<column>.html
```

## Run directory format (`amps_nightly/records.py`)

`run.json`:

| field | meaning |
|---|---|
| `schema_version` | 1 |
| `run_id` | unique id, `YYYYMMDDTHHMMZ-<machine>-<config>` |
| `machine`, `config` | column; `[A-Za-z0-9][A-Za-z0-9._+-]*` |
| `compiler`, `toolchain` | display text, e.g. `GNU`, `GCC 13.2 · OpenMPI 4.1` |
| `branch`, `commit` | tested branch and full SHA |
| `commit_inferred` | true when the commit was reconstructed (legacy archive), shown as `≈` |
| `started`, `finished` | ISO-8601 with offset, stored as UTC |
| `complete` | false when the run was cut short |
| `producer`, `notes` | who wrote the record |

`results.jsonl`, one object per line: `family`, `test`, `status`, and optionally `stage`
(`compile`, `rundir`, `run`, `check`, …), `reason`, `source`, `duration_s`, `cases_total`,
`cases_failed`, `log` (path inside the run directory). `(family, test)` is unique within a run.
Unknown fields are dropped.

The contents of the legacy `.diff` files are **not stored** (decision 2026-11-14): they are read only
to decide the status, the failing stage and a one-line reason.

`DONE` is written last. The server never reads a directory without it, so an interrupted upload
is never half-ingested.

## Test families (`families.json`, `MakefileTest/Table`)

- Legacy tests: each Table block has `Family=<id>` right after `Name=`. `amps-nightly check-table`
  verifies that every block has exactly one registered legacy family (run it after editing the Table).
- `table_parser.cpp` writes `TEST<name>FAMILY=<id>` into `Makefile.test` and a `test_manifest` target:
  `make -s test_manifest` prints `name|family|ref|exception-code|type` for every test selected for this
  host and compiler. The runner uses it instead of re-implementing the `@@(machines)` and `{compiler}` rules.
- `src/test` suites and application campaigns are assigned by `suite_prefixes` (first match wins).
- Results whose family is not registered are shown in an "Other" group, never dropped.

## Machines and the routine (`machines/`, `routine.json`)

The aim is that adding or removing a machine or a test touches one file, and `amps-nightly check` /
`amps-nightly plan` show the effect before anything runs.

**Machine file** `machines/<machine>.json` (`amps_nightly/machines.py`; keys starting with `_` are comments):

| field | meaning |
|---|---|
| `machine` | name, equal to the file name; first half of the column id |
| `enabled` | `false` pauses the machine: nothing runs, its columns leave the report (history is kept) |
| `hostnames` | `hostname` values of the machine; the runner and `plan` find their machine file by them |
| `table_host` | name used by the `@@(...)` selectors of `MakefileTest/Table` (default: `machine`); the runner passes it to `table_parser` as `AMPS_TABLE_HOST`, so a new machine needs no edit of the host aliases in `table_parser.cpp`, and can reuse another machine's Table settings |
| `workdir`, `upload`, `start` | work area, upload target (local directory or `[user@]host:dir` for rsync), nightly start time |
| `core_profile` | the machine's `test/qualification-machine/*.json` for `src/test/run_tests.py` (that tool rejects unknown keys, so its settings stay in its own file); `null` reports core suites as skipped |
| `capabilities` | what the machine offers (`spice`, `swmf-data`, `gpu`, …); sources that `require` a missing one are reported as skipped |
| `lanes`, `timeouts_s` | parallel test lanes; per-stage time limits (`compile`, `rundir`, `run`, `check`) |
| `compile_lanes`, `compile_jobs` | Makefile.test tests compiled at once (each in its own copy of the workspace), and `make -j` of each compile (`TESTCOMPILEJOBS`; default `cores / (compile_lanes x max_parallel)`; a `configs[]` entry may set its own `compile_jobs`) |
| `configs[]` | one entry per compiler configuration = one report column `<machine>-<config>`: `compiler`, `toolchain` (display), `setup` (module loads), `config_pl` (Config.pl calls), `after_install`, optional `sources` (only these sources), `enabled` |

**Routine file** `routine.json` (`amps_nightly/routine.py`):

- `sources[]` — what the runner executes, each through one adapter (`kind`):
  `makefile` (the Table tests: `make -s test_manifest` of the configured tree),
  `core` (`src/test/run_tests.py`, a `campaign` and optional `match` globs on suite IDs),
  `command` (any other test program: `command` run in `cwd`, results from JUnit XML or the exit code;
  its results go to `family`). A source has `enabled`, `requires` (capabilities) and a `note`.
  The application runners are declared but switched off until their adapters are verified (step 3).
- `disabled[]` — tests switched off without being deleted: `test` (glob), `reason` (required), optional
  `source`, `machines`, `configs`, `until` (last day; afterwards the test runs again and `check` warns).

**Resolution** for each enabled configuration (`amps_nightly/plan.py`; `amps-nightly plan` prints it, and the
runner will use the same code): every enabled source that the configuration selects lists its tests; a test is
reported as SKIPPED with a reason when the source needs a capability the machine lacks, a core source runs on a
machine without `core_profile`, the Table marks it `ExeptionCode=SKIP` for the compiler, or a `disabled` entry
matches; otherwise it runs. Two sources producing the same (family, test) is an error. Without a configured
tree, `plan` evaluates the Table's `@@(...)` rule itself (same substring rule as `table_parser.cpp`; checked
against `make test_manifest`: 100 tests for passion, 95 for an unknown host).

**Report columns**: every enabled configuration of every enabled machine file, in file-name order. A declared
column without runs is shown as "No results yet"; stored columns that no enabled file declares are hidden and
named in a note under the table. Without machine files the report falls back to "columns with a run in the
last 60 days".

`amps-nightly check` validates all four files together: Table families, families of `command` sources, core
campaigns and `match` globs against `src/test/manifest.json`, `sources` and `core_profile` of each machine,
duplicate host names, `disabled` entries that expired, name unknown machines or match no test, and machines
whose `table_host` appears in no Table selector.

## Server

- `amps-nightly ingest` validates each complete upload, moves it to `www/archive/YYYY/MM/DD/<machine>-<config>/`
  (the local date the night started; same layout as the legacy `AMPS_TEST_RESULTS`), writes its `index.html`
  and inserts it into `results.db`. Invalid or duplicate uploads go to `rejected/` with `REJECTED.txt`.
- `amps-nightly report` rebuilds the pages from the store:
  - the columns declared by the machine files (above); a column whose newest run is older than 36 h is flagged;
  - "Changes since the previous run" per column: PASS → FAILED (new failure) and FAILED → PASS (fixed);
  - for each test: current status, failing stage and one-line reason, the date of the newest
    passing run in that column and that run's commit;
  - pages are written to a temporary directory and swapped in, so readers never see a partial page.
- `amps-nightly update` = ingest + report (for cron). Pages are static and also work from `file://`.

## Legacy importers (step 1, `amps_nightly/legacy.py`)

They let the new pages show real data before the new runner exists.

- `import-legacy-archive` converts the server archive `AMPS_TEST_RESULTS/YYYY/MM/DD/<folder>/test_amps.res`.
  Each archived night's commit is inferred from the checkout time printed at the top of
  `test_amps.log` and the first-parent history of `master` (`--history`, a full or blob-less clone).
  Nights that repeat the previous night's checkout (stale files re-sent by the old scp) are skipped.
- `import-legacy-run` converts the `*.diff` files left in a legacy work tree after a night
  (exact commit from `git rev-parse HEAD`).
- Result evaluation follows the legacy rules: an empty `.diff` is a pass; otherwise the stage markers
  and the `DiffNum.pl` output give the failing stage and reason. Extra result files a test writes
  (`test_<name>--DUST.diff`, `test_<name>1.diff`, files named in its Table block) are folded into it.
  Only the status, stage and reason are kept; the `.diff` text itself is discarded.

Limits of imported history: only the first 100 lines of each `.diff` were archived (difference counts
are then lower bounds, shown as "at least N"), the archive only knows when results reached the server
("Results collected"), and the archived column folder `amps-gpu_gnu` was produced on passion
(checkout times are in Central time), so it is imported as `passion-gnu`.

## Next steps

The orchestrator (runners file, workspaces, readers, log archive, sequential then concurrent
execution) is designed in [ORCHESTRATOR.md](ORCHESTRATOR.md); it supersedes the runner sketch below.

2. **Runner** (`amps-nightly run [--machine <name>]`): reads the machine file and the plan above, then per
   configuration: checkout, `Config.pl` calls with `AMPS_TABLE_HOST=<table_host>`, the sources, a run
   directory with `SKIPPED` rows for everything the plan skips, upload; adapters:
   - `makefile`: per test `make test_<X>_compile/_rundir/_run/_check` with timeouts, using
     `make test_manifest`; independent lanes in separate `git worktree`s so compiles run in parallel;
   - `core`: `src/test/run_tests.py --config <core_profile> --campaign <campaign> --json`, suites mapped by
     `suite_prefixes`;
   - `command`: the source's command with its time limit; JUnit cases or the exit code become results;
   - upload with rsync, `DONE` last.
3. Application adapters (srcSEP, srcSEP3D, srcEarth, srcMoon, srcEnceladus) and more machines (pleiades).
4. Retire `utility/TestScripts`.
