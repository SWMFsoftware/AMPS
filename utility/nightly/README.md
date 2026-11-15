# utility/nightly

New AMPS nightly testing: result records, results store and report pages. Design: [DESIGN.md](DESIGN.md);
test orchestrator: [ORCHESTRATOR.md](ORCHESTRATOR.md).
Look-and-feel mockup: [mockup/index.html](mockup/index.html).

Everything is one tool, `amps-nightly` (Python 3.8+, standard library only):

```sh
utility/nightly/amps-nightly check            # Table, families.json, runners.json and machines/ agree
utility/nightly/amps-nightly machines         # machine files, their report columns and last uploads
utility/nightly/amps-nightly plan [--list]    # what the routine runs and skips on this machine, and why
utility/nightly/amps-nightly run [--dry-run]  # 1: the nightly test procedure on this machine (runners.json)
utility/nightly/amps-nightly pack             # 2: tar each finished run with its logs, send it to the server
utility/nightly/amps-nightly update           # 3: (server) take the tarballs, update results.db, archive, pages
utility/nightly/amps-nightly --help
```

The three nightly steps also exist as scripts, like `run_test_amps.sh`, `scp_test_amps.sh` and
`process_amps_tests.sh` of the old routine: `nightly_run.sh`, `nightly_pack.sh` and `nightly_update.sh`
(each passes its arguments on to the command above).

What runs where is decided by four files:

| file | decides |
|---|---|
| `machines/<machine>.json` | the test machines; each configuration in it (compiler) is one report column |
| `runners.json` | the test runners the nightly executes, in order, with their commands; tests switched off for now (`disabled`) |
| `MakefileTest/Table` | the Makefile.test tests, which machines run them (`@@(...)`), their family (`Family=`) |
| `families.json` | the rows of the report |

## Adding and removing machines

**Add a machine** (no other file needs to change):

1. `cp machines/passion.json machines/<name>.json` and edit it: `machine` (= file name), `hostnames`
   (`hostname` output there), `workdir`, `upload` (`passion.ndc.nasa.gov:AMPSTEST/nightly/incoming`
   from any other machine), `start`, `lanes`, `timeouts_s`, `capabilities`, and one entry in `configs` per
   compiler (`config_pl` = the `Config.pl` calls, `setup` = module loads).
2. src/test suites need that machine's profile: `cp test/qualification-machine/template.json
   test/qualification-machine/<host>.json`, edit it, and set `core_profile` to it (or `null`: core suites
   are then reported as skipped there).
3. `table_host` is the name the Table's `@@(0+...)` selectors use for this machine. Use a new name and add it
   to the selectors of the tests that need machine-specific `Ref=`/`Keys=` lines, or reuse the name of a
   machine with the same paths (e.g. `"table_host": "passion"`).
4. `amps-nightly check` and `amps-nightly plan --machine <name> --list`, then commit.
   The report shows the new columns at once, marked "No results yet" until the first upload.

**Pause a machine**: `"enabled": false` in its file (or in one entry of `configs` to pause one compiler).
Nothing runs there, and its columns leave the report; their history stays in the results store and comes
back when it is switched on again. **Remove a machine**: delete its file (and its cron line).

## Adding and removing tests

| to … | Makefile.test tests | src/test suites | another test program |
|---|---|---|---|
| add | a block in `MakefileTest/Table` with `Name=`, `Family=` (+ reference files) | add `nightly` to the suite's tiers in `src/test/manifest.json` | a `"kind": "command"` source in `routine.json` and its family in `families.json` |
| run only on some machines | `Name=<test>@@(0+passion+pleiades)` | `requires.capabilities` of the suite (run_tests.py checks them against the machine's profile) | `requires` a capability, or list the source in the machines' `sources` |
| switch off for now (shown as SKIPPED with the reason) | an entry in `routine.json` `disabled` | same | same |
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

## Running the nightly procedure

Three steps, like the three scripts of the old routine:

1. **`nightly_run.sh`** (`amps-nightly run`) checks out the head of `master` (local mirrors in
   `{workdir}/mirrors`), then for each enabled configuration runs the runners of `runners.json` one after
   another, each in its own worktree. Results and logs go to `{workdir}/nights/<YYYYMMDD-HHMM>/runs/<column>/`.
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

Crontab on passion, which is test machine and server (the night starts at 20:00 and stops starting runners
after `night_limit_h` = 14 h):

```
# m  h  dom mon dow  command
00 20 * * *  D=/data/vtenishe/TEST/AMPS/utility/nightly; $D/nightly_run.sh; $D/nightly_pack.sh; $D/nightly_update.sh
```

Another test machine runs `nightly_run.sh; nightly_pack.sh` with `"upload": "passion.ndc.nasa.gov:AMPSTEST/nightly/incoming"`,
and passion adds `30 10 * * * $D/nightly_update.sh` to pick its tarballs up in the morning.

The tool, `runners.json` and the machine files are read from the checkout that holds `amps-nightly`; the
tests run on the commit at the head of `master`. Workspaces of runners with failures are kept until the next
night starts (`keep_workspaces`).

**Adding a runner** = one entry in `runners.json` (fields: ORCHESTRATOR.md section 2). If its report has a new
format, add a reader function to `amps_nightly/readers.py`. Try it with `run --only <id> --no-upload`.

## Server data

Server data lives in `--root` (default `$AMPS_NIGHTLY_ROOT`, else `~/AMPSTEST/nightly`):

```
incoming/   <run_id>.tar.gz from `pack` (hidden names are still arriving); run directories with DONE also work
rejected/   uploads that failed validation, with REJECTED.txt
results.db  SQLite store of every run and result
www/        index.html, family/, assets/, archive/   (open www/index.html in a browser)
www/archive/YYYY/MM/DD/<machine>-<config>/   one night of one column, like the legacy AMPS_TEST_RESULTS:
            index.html (night summary), run.json, results.jsonl, runners.json, logs/<runner>/...
```

On passion, Firefox runs under an AppArmor profile that lets it read files under `$HOME` and `/tmp`
but not under `/data`: a page there just keeps loading. Keep `--root` under `$HOME` (the default), and to
view `mockup/index.html` from a `/data` checkout use `google-chrome` or a copy under `$HOME`.

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
