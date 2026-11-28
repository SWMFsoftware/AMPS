# utility/nightly

New AMPS nightly testing: result records, results store and report pages. Design: [DESIGN.md](DESIGN.md);
test orchestrator: [ORCHESTRATOR.md](ORCHESTRATOR.md).

Everything is one tool, `amps-nightly` (Python 3.8+, standard library only):

```sh
utility/nightly/amps-nightly check            # Table, families.json, runners*.json and machines/ agree
utility/nightly/amps-nightly machines         # machine files, their report columns and last uploads
utility/nightly/amps-nightly plan [--list]    # what the routine runs and skips on this machine, and why
utility/nightly/amps-nightly run [--dry-run]  # 1: the nightly test procedure on this machine (runners.json)
utility/nightly/amps-nightly run --routine long   # the same for the long (periodic) routine (runners.long.json)
utility/nightly/amps-nightly pack             # 2: tar each finished run with its logs, send it to the server
utility/nightly/amps-nightly update           # 3: (server) take the tarballs, update results.db, archive, pages
utility/nightly/amps-nightly --help
```

For planning, an explicit `--table /path/to/Table` controls test selection even
if the default checkout contains a previously generated `Makefile.test`.
An explicit `--amps-dir /path/to/configured/AMPS` instead requests that tree's
exact compiler-selected manifest and checks its configured host; it retains
priority when both options are supplied. With neither override, an existing
manifest for the matching host is used as before. This precedence prevents
unrelated local configuration from changing a deliberately selected Table.

The three nightly steps also exist as scripts, like `run_test_amps.sh`, `scp_test_amps.sh` and
`process_amps_tests.sh` of the old routine: `nightly_run.sh`, `nightly_pack.sh` and `nightly_update.sh`
(each passes its arguments on to the command above).

What runs where is decided by these files:

| file | decides |
|---|---|
| `machines/<machine>.json` | the test machines; each configuration in it (compiler) is one report column; `routines`: the routines other than the nightly that run there, and when |
| `runners.json` | the test runners the nightly executes, in order, with their commands; tests switched off for now (`disabled`) |
| `runners.<routine>.json` | the same for another routine, e.g. `runners.long.json`: the long tests, run periodically (section "Long tests") |
| `MakefileTest/Table` | the Makefile.test tests, which machines run them (`@@(...)`), their family (`Family=`), their tier (`Tier=long`, `Ranks=`) |
| `families.json` | the rows of the report |

## Adding and removing machines

**Add a machine** (no other file needs to change):

1. `cp machines/passion.json machines/<name>.json` and edit it: `machine` (= file name), `hostnames`
   (`hostname` output there), `workdir`, `upload` (`passion.ndc.nasa.gov:/data/vtenishe/AMPSTEST/incoming`
   from any other machine), `start`, `lanes`, `compile_lanes`, `timeouts_s`, `capabilities`, and one entry in
   `configs` per compiler (`config_pl` = the `Config.pl` calls, `setup` = module loads). `compile_jobs` (machine,
   or one `configs` entry) sets `make -j` of the Makefile.test compiles; by default the compile lanes share the
   `cores`.
2. src/test suites need that machine's profile: `cp test/qualification-machine/template.json
   test/qualification-machine/<host>.json`, edit it, and set `core_profile` to it (or `null`: core suites
   are then reported as skipped there).
3. `table_host` is the name the Table's `@@(0+...)` selectors use for this machine. Use a new name and add it
   to the selectors of the tests that need machine-specific `Ref=`/`Keys=` lines, or reuse the name of a
   machine with the same paths (e.g. `"table_host": "passion"`).
4. `amps-nightly check` and `amps-nightly plan --machine <name> --list`, then commit and push.
   The report shows the new columns at once, marked "No results yet" until the first upload.
5. On the new machine: clone AMPS anywhere and run `utility/nightly/amps-nightly install` (section
   "Installing the nightly on a machine").

**Pause a machine**: `"enabled": false` in its file (or in one entry of `configs` to pause one compiler).
Nothing runs there, and its columns leave the report; their history stays in the results store and comes
back when it is switched on again. **Remove a machine**: `amps-nightly uninstall` there, then delete its file.

## Adding and removing tests

| to … | Makefile.test tests | src/test suites | another test program |
|---|---|---|---|
| add | a block in `MakefileTest/Table` with `Name=`, `Family=` (+ reference files) | add `nightly` to the suite's tiers in `src/test/manifest.json` | a `"kind": "command"` source in `runners.json` and its family in `families.json` |
| run only on some machines | `Name=<test>@@(0+passion+pleiades)` | `requires.capabilities` of the suite (run_tests.py checks them against the machine's profile) | `requires` a capability, or list the source in the machines' `sources` |
| switch off for now (shown as SKIPPED with the reason) | an entry in `runners.json` `disabled` | same | same |
| remove | delete the block, or `Name=<test>@@(0+)` | remove `nightly` from its tiers | delete the source, or `"enabled": false` |

A `disabled` entry needs a `reason`; `until` (a date) makes it lapse by itself, and `check` warns once it has:

```json
"disabled": [
  {"test": "SEP-CME--*", "reason": "reference files are being regenerated", "until": "2026-12-15",
   "machines": ["passion"], "configs": ["gnu"], "source": "legacy"}
]
```

`test` is a shell glob on the test name (Table name or src/test suite ID); `machines`, `configs` and
`source` are optional filters. A removed test disappears from the report; its history stays in the store.

The Makefile.test tests are compiled one after another in one workspace. Before each one the runner puts
the files Config.pl rewrites (`Makefile`, `Makefile.local`, `Makefile.conf`, `.amps.conf`, `.general.conf`)
back to their state after the install, and the test's log says so when an earlier test had changed them.
So a test's compile rule must set everything it needs (`-application`, data paths, `-openmp`, link and
compiler options) and cannot rely on what an earlier test configured.

## Long tests (the periodic "long" routine)

Tests too long for a night (the 2013-04-11 SEP real-event runs: a full event on 32 MPI ranks) run in a second
routine, `long`, once a week. It is generic: any application's tests can join it.

**A Makefile.test test** becomes a long test with two lines in the keyword part of its `MakefileTest/Table` block
(before `Compile=`/`Rundir=`/`Run=`/`Check=`):

```
#>
Name=SEP--event-2013-04-11-long--Parker
Family=sep-long
Tier=long
Ranks=32
...
Run=
>>>cd $(TESTSEP--event-2013-04-11-long--ParkerDIR) && $(TESTMPILAUNCH) $(TESTSEP--event-2013-04-11-long--ParkerRANKS) ./amps ...
```

- `Tier=long` puts the block into `Makefile.test.long` (its list: `make -s test_long_manifest`) instead of
  `Makefile.test`, so the nightly never runs it. `make test_<name>` works from the AMPS root as for any test.
- `Ranks=<N>` (required with `Tier=long`, not allowed without) is the test's fixed number of MPI ranks: its run
  rule launches `$(TESTMPILAUNCH) $(TEST<name>RANKS) ./amps` (`TESTMPILAUNCH` defaults to `mpirun -np`), and its
  references hold for that count.
- `Family=` names a family of the group `long` in `families.json` (e.g. `sep-long`); nightly tests use the
  group `legacy`.
- The name must not start with the name of a nightly test (`X--long` would be counted as a result file of `X`).

**Another application's test program** joins the long routine as a runner entry in `runners.long.json` (same
fields as in `runners.json`, ORCHESTRATOR.md section 2; runner ids must differ from the nightly's), with its
results in a family of the group `long`. `runners.long.json` holds one runner now, `legacy-long`: the built-in
Makefile.test runner with `"manifest_target": "test_long_manifest"`.

**Where and when it runs** is the `routines` object of a machine file; a machine without it does not run the
long routine. On passion:

```json
"routines": {"long": {"enabled": true, "start": "Sat 10:00", "time_limit_h": 32, "lanes": 1, "compile_lanes": 1,
                      "timeouts_s": {"compile": 3600, "rundir": 600, "run": 86400, "check": 1800},
                      "make_vars": {"TESTMPILAUNCH": "mpirun -np"}, "max_ranks": 16, "stale_after_h": 192}}
```

While the routine runs, `time_limit_h` replaces `night_limit_h`, `lanes` and `compile_lanes` replace the
machine's, and `timeouts_s` and `make_vars` are merged into the machine's. `start` is `<day> HH:MM` (weekly) or
`HH:MM` (daily). `max_ranks` (default: `cores`) is the most ranks a long test may use there (`check` reports a
larger `Ranks=`). `stale_after_h` (default 192) is when the report flags the routine's results as old.

**By hand** (results stay on the machine; `pack` does not send them):

```sh
utility/nightly/amps-nightly plan --routine long --list
utility/nightly/amps-nightly run --routine long --dry-run
utility/nightly/amps-nightly run --routine long --only legacy-long --test 'SEP--event-2013-04-11-long--Parker' --no-upload
```

`run --routine long` works in `{workdir}/nights/<YYYYMMDD-HHMM>-long/`; its run ids end in `-long`
(`20261128T1600Z-passion-gnu-long`) and its runs are recorded with `"routine": "long"`. The report shows the
nightly table as before (the newest **nightly** run of each column), and below it a section "Long validation
runs (periodic)" with each long test's newest result, the date and commit of the run it comes from, and when
it last passed.

## Running the nightly procedure

Three steps, like the three scripts of the old routine:

1. **`nightly_run.sh`** (`amps-nightly run`) checks out the head of `master` (local mirrors in
   `{workdir}/mirrors`), then for each enabled configuration runs the runners of `runners.json` one after
   another. Each runner gets its own installed copy, in one folder per compiler as in the old routine:

   ```
   ~/Tmp_AMPS_test/                       the machine file's "workdir" (TESTDIR of the ~/bin script)
       AMPS/  AMPS.prev/                  tonight's download (the orchestrator, runners.json, machines/)
       mirrors/                           AMPS.git share.git util.git BATL.git
       GNU/<runner>/AMPS                  installed copy of each runner (Intel/<runner>/AMPS, ... for more
       GNU/<runner>/out  GNU/<runner>/tmp   compilers); replaced the next night, deleted early if all passed
       nights/<YYYYMMDD-HHMM>/            orchestrator.log, runs/<column>/ (results and logs until pack sends them)
       outbox/                            the tarballs
   ```

   No runner, and no Makefile.test test, starts with less than `min_free_gb` (default 5 GB) free there;
   it is reported SKIPPED with the reason.
   Each step of a log keeps its first 256 KB and last 768 KB; the middle of very long output is left out.
2. **`nightly_pack.sh`** (`amps-nightly pack`) writes every finished run not sent yet into
   `{workdir}/outbox/<run_id>.tar.gz` (run.json, results.jsonl, runners.json, logs/: everything the pages
   show) and copies it to the machine file's `upload` destination (a local directory, or
   `[user@]host:directory` through rsync). A tarball that could not be sent stays in the outbox and goes
   next time; sent ones are deleted after 14 days. A run whose `run` died (reboot, killed) is completed from
   the progress saved after each runner and sent marked incomplete.
3. **`nightly_update.sh`** (`amps-nightly update`, on the server) unpacks every tarball in `incoming/` (or
   the tarballs given as arguments), checks it, adds it to `results.db`, moves it to
   `www/archive/YYYY/MM/DD/<column>/`, rebuilds the pages and gzips logs older than 14 nights.

```sh
utility/nightly/nightly_run.sh --dry-run                      # the commands each runner would execute here
utility/nightly/nightly_run.sh --only enceladus --no-upload   # try one runner (pack will not send it)
utility/nightly/nightly_run.sh --only legacy --test 'fast-wave*' --no-upload   # a few Makefile.test tests
utility/nightly/nightly_run.sh                                # the whole night
utility/nightly/nightly_pack.sh --no-send                     # only write the tarballs
utility/nightly/nightly_update.sh ~/somewhere/<run_id>.tar.gz # a tarball brought by hand
```

Workspaces of runners with failures are kept until the next night starts (`keep_workspaces`).

## Installing the nightly on a machine (cron)

As in the old routine (`Config.pl -set-test`, `~/bin/run_test_amps.sh`), cron starts one script in `~/bin`.
Its version-controlled original is `bin/amps_nightly_test.sh`; the copy in `~/bin` is installed with

```sh
<AMPS>/utility/nightly/amps-nightly install --dry-run       # show what would change
<AMPS>/utility/nightly/amps-nightly install                 # make the changes (asks first)
```

which writes `~/bin/amps_nightly_test.sh`, its settings `~/bin/amps_nightly_test.conf` (`AMPSTEST`,
`BRANCH`, `UPDATE_SITE`, `PYTHON`), and one crontab line per routine below a marker comment: the nightly's at the
machine file's `start`, and one for each routine its `routines` enables, at that routine's `start`, passing
`--routine <name>`. The crontab is saved to `$AMPSTEST/backup/` first. `amps-nightly uninstall` removes them
again. `install --time HH:MM` changes the nightly's time, `install --routine long --time 'Sun 06:00'` the long
routine's. On passion:

```
# amps-nightly: AMPS nightly testing (added by 'amps-nightly install'; remove with 'amps-nightly uninstall')
00 20 * * *  /home/vtenishe/bin/amps_nightly_test.sh >> /data/vtenishe/AMPSTEST/logs/cron.log 2>&1
00 10 * * 6  /home/vtenishe/bin/amps_nightly_test.sh --routine long >> /data/vtenishe/AMPSTEST/logs/cron.log 2>&1
```

Each night the script

1. takes a lock (one night at a time, of any routine: two routines never run at once on a machine). The
   nightly gives up when the lock is taken; another routine (`--routine long`) waits for it up to
   `LOCK_WAIT_H` hours (default 6), then gives up; both say so in the log;
2. downloads AMPS: `git clone --depth 1` of `BRANCH` into `$TESTDIR/AMPS` (`~/Tmp_AMPS_test/AMPS`; the
   night before is kept as `AMPS.prev` and used if the download fails);
3. pulls SWMF_data: the nightly's own copy at the machine file's `swmf_data`, reset to the newest commit
   and cleaned every night, because some Makefile.test checks edit reference files in it
   (`RemoveColumn.pl`, `perl -i.tmp`) and later runs would see the edits;
4. writes the night's plan to the log (where each runner is installed and with which `Config.pl` calls:
   `run --dry-run`), runs `$TESTDIR/AMPS/utility/nightly/nightly_run.sh --commit <downloaded commit>`, then
   `nightly_pack.sh`,
   and on the server (`UPDATE_SITE=yes`) `nightly_update.sh`;
5. logs everything to `$AMPSTEST/logs/nightly-YYYYMMDD-HHMM.log` (kept 60 days); another routine to
   `$AMPSTEST/logs/<routine>-YYYYMMDD-HHMM.log` (`long-20261128-1000.log`).

So the orchestrator, `runners.json` and the machine files come from the branch on GitHub, and a change to them
takes effect the night after it is pushed; a working checkout does not affect the nightly. The tests run on
the same commit as the tool. The `~/bin` script itself changes only by `amps-nightly install`: when the
branch has a different version of it, the night's log and summary page show a warning with the command to run.
Arguments of the script go to `nightly_run.sh`, so `~/bin/amps_nightly_test.sh --only moon --no-upload` is a trial
of the whole chain.

Another test machine installs with `"upload": "passion.ndc.nasa.gov:/data/vtenishe/AMPSTEST/incoming"` in its
machine file (so `UPDATE_SITE=no`), and passion picks its tarballs up with one more crontab line,
`30 10 * * * <AMPS>/utility/nightly/nightly_update.sh` (or with its own night's update).

**Adding a runner** = one entry in `runners.json` (fields: ORCHESTRATOR.md section 2). If its report has a new
format, add a reader function to `amps_nightly/readers.py`. Try it with `run --only <id> --no-upload`.

## Server data

Server data lives in `--root` (default `$AMPS_NIGHTLY_ROOT`, else the `server_root` of this host's machine
file, else `~/AMPSTEST/nightly`). On passion that is `/data/vtenishe/AMPSTEST`, which also holds the nightly's
`logs/` and `SWMF_data/` (`~/AMPSTEST/nightly` is a link to it, for the nights before the move):

```
incoming/   <run_id>.tar.gz from `pack` (hidden names are still arriving); run directories with DONE also work
rejected/   uploads that failed validation, with REJECTED.txt
results.db  SQLite store of every run and result
www/        index.html, family/, assets/, archive/   (open www/index.html in a browser)
www/archive/index.html   every stored night of every column, newest first, linked from "All nights"
www/archive/YYYY/MM/DD/<machine>-<config>/   one night of one column, like the legacy AMPS_TEST_RESULTS:
            index.html (night summary, with links to the previous and next night of the column),
            nav.json (those links), run.json, results.jsonl, runners.json, logs/<runner>/...
www/archive/YYYY/MM/DD/<machine>-<config>-<routine>/   one run of another routine (passion-gnu-long), linked to
            the previous and next run of that routine
```

On passion, Firefox runs under an AppArmor profile that lets it read files under `$HOME` and `/tmp`
but not under `/data`: a page there just keeps loading (a link from `$HOME` does not help). Open the report
with `google-chrome /data/vtenishe/AMPSTEST/www/index.html`, or serve it and use any browser:
`python3 -m http.server 8000 --directory /data/vtenishe/AMPSTEST/www`, then `http://localhost:8000/`.

## Feeding it from the legacy nightly (until the new runner exists)

One-time backfill from the legacy archive on passion (commits are inferred from checkout times,
which needs a clone with full history, e.g. `git clone --bare --filter=blob:none git@github.com:SWMFsoftware/AMPS.git`):

```sh
utility/nightly/amps-nightly import-legacy-archive \
    --archive ~/AMPSTEST/Sites/AMPS_TEST_RESULTS \
    --map amps-gpu_gnu=passion/gnu --map amps-gpu_cuda=passion/cuda \
    --history /path/to/AMPS-history.git
utility/nightly/amps-nightly update
```

After each legacy night (suggested crontab on passion, after the old `run_test_amps.sh` has finished):

```
# m  h  dom mon dow  command
15 09 * * *  A=/path/to/AMPS/utility/nightly; /usr/bin/python3 $A/amps-nightly import-legacy-run --amps-dir $HOME/Tmp_AMPS_test/GNU/AMPS --machine passion --config gnu && /usr/bin/python3 $A/amps-nightly update
```

`import-legacy-run` takes the families from `<amps-dir>/MakefileTest/Table`; until the `Family=` lines are
committed and the nightly checks them out, add `--table /path/to/AMPS/MakefileTest/Table`.

## Tests

```sh
cd utility/nightly && /usr/bin/python3 -m unittest discover -s tests -t .
```
