# srcSEP3D Testing Procedure

`test/run_tests.py` is the single user-facing test interface. Its selectors
match `srcSEP/test/run_tests.py` so the two applications can use the same
automation habits even though srcSEP3D currently has only the R0 evidence set.

## Quick commands

```bash
python3 test/run_tests.py --list
python3 test/run_tests.py --routine --amps-source /path/to/AMPS
python3 test/run_tests.py --test LAY01
python3 test/run_tests.py --group HARN --group BLDL3D \
  --amps-source /path/to/AMPS
python3 test/run_tests.py --all --amps-source /path/to/AMPS \
  --output-dir test_output/all
```

The runner does not accept an implicit mode. Choose exactly one of `--list`,
`--test`/`--group`, `--routine`, `--all`, or `--suite`. This prevents an empty
or misspelled selection from exiting successfully.

## CLI correspondence with srcSEP

| srcSEP selector | srcSEP3D behavior |
|---|---|
| `--list` | prints stable IDs, groups, evidence kind, and suites |
| `--test ID` | runs one ID; repeatable and case-insensitive |
| `--group GROUP` | expands a registered group; repeatable |
| `--routine` | bounded normal-development set; excludes shell exit-code probes |
| `--all` | runs every public test separately and continues after failures |
| `--suite NAME` | runs one dependency/evidence class; repeatable |
| `--output-dir DIR` | writes the complete report bundle below `DIR` |
| `--amps PATH` | reserves the linked-executable path for later linked phases |
| `--timeout SEC` | applies a per-command timeout |

Additional R0 setup options are `--amps-source`, `--make-config`,
`--sep-common-dir`, and `--sep-common-archive`.

## Evidence classes

### Standalone C++ registry

`test/stage1` is compiled with no AMPS include path and no MPI library. It links
only the retained L0 source, test callbacks, and the shared `sep_common.a` test
registry. The runner adds `-Wall -Wextra -Wpedantic -Werror`.

| Group | IDs | Purpose |
|---|---|---|
| `HARN` | `HARN01`–`HARN04` | registry selection, PASS/FAIL/SKIP/ERROR exits, JSON, and JUnit |
| `LAY` | `LAY01`, `LAY02` | L0/L1 dependency rule and a negative control proving the guard fires |
| `BLD` | `BLD01` | `nm -u` confirms the standalone binary has no AMPS/MPI symbols |
| `UTIL` | `UTIL02` | byte-exact shared-kernel reference record |
| `RUNNER` | `RUN3D01` | Python selector, de-duplication, usage-error, JSON, and JUnit contract |

`HARN02-EXITCODE` and `HARN03-EXITCODE` are shell-level probes. They launch the
internal fail/skip beacons and verify the outer process status. The beacons are
not public `--all` tests because one is intentionally failing.

### Phase R0 source and ABI gates

#### BLDL3D01 — configured enclosing AMPS build

The runner locates or accepts a configured `Makefile.conf`, then executes:

```bash
make -f makefile strict-production \
  AMPS_ROOT=/path/to/AMPS \
  AMPS_CONFIG=/path/to/AMPS/Makefile.conf
```

That target invokes the top-level `make amps`, which owns generated headers,
include paths, compile definitions, libraries, and the final executable link.
It then audits `AMPS/build/main/mainlib.a` and `main.a` with `nm`. A direct
`make lib` from `AMPS/srcSEP3D` is not production evidence because it does not
inherit the enclosing include configuration and cannot reliably locate
`build/pic/pic.h`. If the real configuration is absent, the result is SKIP. A
source-only compile or mock header is never reported as a production pass.

#### BLDL3D02 — retired-source exclusion

This check reads the live source tree and active makefile lines. It requires:

- no `SEP3D.cpp` file or makefile object;
- no retired mover, acceleration, sampler, or field-line injection symbol in
  production headers/sources;
- no former wedge bounds;
- no `PrepopulateDomain`, placeholder mesh output, or mesh-file write in
  `main_lib.cpp`;
- only the R0 production manifest described in the root README.

The check deliberately scans production code, not documentation, because the
migration record must be allowed to name what was removed.

#### BLDL3D03 — AMPS mover status mapping

This check verifies both sides of the boundary:

1. `amps/amps_mover_status.h` contains the required `static_assert` and mapping
   branches.
2. The supplied AMPS `src/pic/pic.h` defines deleted-on-face `0`, left-domain
   `2`, and motion-finished `3`.

Use `--amps-source /path/to/AMPS` when the source is not in a discoverable
parent directory. Missing `pic.h` produces SKIP, not PASS.

#### BLDL3D04 — AMPS macro hygiene

AMPS `general/constants.h` defines `Pi` as a preprocessor macro. Namespaces do
not prevent macro substitution, so the former `Const::Pi` declaration broke
unrelated translation units after `pic.h` included the application header.
This compile-only regression defines the same macro before including
`core/sep3d_types.h` and requires the collision-safe `Const::kPi` API.

#### BLDL3D05 — copied-build makefile paths

AMPS copies the application into `AMPS/build/main` before compiling it. This
test constructs both `AMPS/srcSEP3D/makefile` and copied
`AMPS/build/main/makefile` fixtures, then invokes each with `make -f` from an
unrelated working directory. Both must resolve the same absolute `AMPS_ROOT`,
`AMPS_CONFIG`, and `SEP_COMMON_DIR`; the active makefile directory itself must
match its source or copied location. This directly protects against the
`build/main: ../Makefile.conf: No such file or directory` failure.
The same fixture also supplies a synthetic enclosing `amps` target and verifies
that `strict-production` delegates to it and audits `build/main` archives rather
than attempting a bare compile in `srcSEP3D`.

## Named suites

| Suite | Contents |
|---|---|
| `standalone` | C++ registry, shell exit-code probes, and `RUN3D01` |
| `r0` | R0 source/ABI/production gates plus RUN3D01, LAY01, and BLD01 |
| `production` | BLDL3D01–05 |

Suites can be repeated. Overlapping IDs are de-duplicated in stable order.

## Shared utility discovery

In the normal AMPS layout, the runner discovers `src/models/sep_common` from
the AMPS root and builds `sep_common.a` when needed. For a detached tree:

```bash
python3 test/run_tests.py --routine \
  --sep-common-dir /path/to/AMPS/src/models/sep_common \
  --sep-common-archive /path/to/sep_common/sep_common.a \
  --amps-source /path/to/AMPS
```

The archive must be the shared build used by the applications. The runner does
not compile private copies of shared kernel sources.

## Reports and exit codes

Every execution writes:

- `srcsep3d-tests.json` — commands, messages, elapsed times, totals, and final
  exit code;
- `srcsep3d-tests.xml` — JUnit representation for CI;
- one native `<ID>.json` for each C++ registry test.

Exit codes are:

| Code | Meaning |
|---:|---|
| 0 | all selected tests passed or were explicitly skipped |
| 1 | at least one test executed and failed an acceptance criterion |
| 2 | setup/usage error or at least one test could not be evaluated |

SKIP is an evidence limitation. It is acceptable in source-only development
but cannot close a release gate that requires the skipped test.

## Adding a test

For a dependency-free C++ test:

1. add a callback file under `test/individual-test/`;
2. declare its registration function in `core/sep3d_test_registry.h`;
3. register it in `test/stage1.cpp`;
4. add its source to the makefile and runner build command;
5. add one `TestDefinition` to `test/run_tests.py`;
6. document setup, independent reference, tolerance, and failure meaning here;
7. run `--list`, the individual ID, its group, `--routine`, and `--all`.

For an AMPS-linked test, do not add a mock replacement to the standalone
binary. Register it in the production executable and make the Python runner
invoke that exact linked callback, following the srcSEP pattern.

## Troubleshooting

| Symptom | Interpretation and action |
|---|---|
| `cannot locate srcSEP utility headers` | pass `--sep-common-dir` |
| `cannot locate sep_common.a` | build the shared archive, then pass `--sep-common-archive` |
| `BLDL3D01 SKIP` | run in a configured AMPS checkout or pass `--make-config` |
| `BLDL3D01` reports `pic.h: No such file or directory` from `srcSEP3D` | update the makefile; the gate must delegate to top-level `make amps`, not compile `main_lib.cpp` directly in the source directory |
| `BLDL3D02` reports stale `SEP3D.cpp` | remove the obsolete file from the installed `srcSEP3D`; overlay extraction does not delete files left by an older version |
| `BLDL3D03 SKIP` | pass `--amps-source` pointing to a tree containing `src/pic/pic.h` |
| `build/main` cannot find `../Makefile.conf` | run `make print-layout-paths`; `AMPS_CONFIG` must resolve to the absolute `AMPS/Makefile.conf` path |
| standalone compile failure | rerun with `--rebuild --verbose` |
| unknown test/group | use `--list`; unknown selectors are usage errors |
| report missing after a C++ test | treat as ERROR; inspect verbose subprocess output |
