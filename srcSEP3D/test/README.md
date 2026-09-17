# srcSEP3D Testing Procedure

`test/run_tests.py` is the single user-facing test interface. Its selectors
match `srcSEP/test/run_tests.py` so the two applications can use the same
automation habits. The runner combines R0/R1/R2 foundation evidence with
Phase-M mesh, Phase-B background, and Phase-T turbulence/coefficient evidence
without requiring AMPS for dependency-free tests.

## Quick commands

```bash
python3 test/run_tests.py --list
python3 test/run_tests.py --routine --amps-source /path/to/AMPS
python3 test/run_tests.py --test LAY01
python3 test/run_tests.py --suite r1 --suite r2
python3 test/run_tests.py --suite phase-m --suite phase-b --suite phase-t
python3 test/run_tests.py --group HARN --group BLDL3D \
  --amps-source /path/to/AMPS
python3 test/run_tests.py --all --amps-source /path/to/AMPS \
  --output-dir test_output/all
```

The runner does not accept an implicit mode. Choose exactly one of `--list`,
`--test`/`--group`, `--routine`, `--all`, or `--suite`. This prevents an empty
or misspelled selection from exiting successfully.

### Parallel production compilation

The `BLDL3D01` gate delegates compilation to the enclosing AMPS GNU Make
build. Set `MAKEFLAGS` for this runner invocation to allow that build and its
recursive make operations to compile in parallel. The portable `env` form
works with `tcsh`, `csh`, `bash`, and `zsh`:

```bash
env MAKEFLAGS="-j16" test/run_tests.py --all --amps-source .. \
  --make-config ../Makefile.conf --output-dir test_output/all --rebuild
```

Adjust `16` to the CPU and memory available on the build host. Parallelism
applies to GNU Make compilation performed by the runner, including independent
standalone object compilation and the enclosing AMPS build. The tests
themselves remain sequential. The shorter `MAKEFLAGS="-j16" command`
assignment-prefix form is valid
in `bash` and `zsh`, but not in `csh` or `tcsh`; use the documented `env` form
for a shell-independent command.

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

Additional setup options are `--amps-source`, `--make-config`,
`--sep-common-dir`, and `--sep-common-archive`.

## Evidence classes

### Standalone C++ registry

`test/stage1` is compiled with no AMPS include path and no MPI library. It links
the Phase-M mesh model, Phase-B providers/snapshots, Phase-T turbulence bridge,
R2 Runtime/adapters, test callbacks, and shared `sep_common.a`. The runner adds
`-Wall -Wextra -Wpedantic -Werror`.

| Group | IDs | Purpose |
|---|---|---|
| `HARN` | `HARN01`–`HARN04` | registry selection, PASS/FAIL/SKIP/ERROR exits, JSON, and JUnit |
| `LAY` | `LAY01`, `LAY02` | core/background/runtime/mesh/turbulence dependency rule and a negative control proving the guard fires |
| `BLD` | `BLD01` | `nm -u` confirms the standalone binary has no AMPS/MPI symbols |
| `UTIL` | `UTIL02` | byte-exact shared-kernel reference record |
| `LIFE3D` | `LIFE3D01`–`LIFE3D04` | immutable configuration, state machine, frozen layout, counters, adapter parity, and no-parser boundary |
| `MSH3D` | `MSH3D01`–`MSH3D09` | resolution, tube geometry, balance, octree budget/ownership, presets, gradients |
| `BGP3D` | `BGP3D01`–`BGP3D06` | analytic Parker field/plasma identities and polar limits |
| `SNAP3D` | `SNAP3D01`–`SNAP3D08` | snapshot completeness, coupling conversion, atomicity, interpolation, batch/frame policy |
| `TUR3D` | `TUR3D01`–`TUR3D04` | spectrum, AWSoM convention, resonance, missing-data policy |
| `COEF3D` | `COEF3D01`, `COEF3D02` | conversion round trips and shared-kernel identity |
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
It then audits `AMPS/build/main/mainlib.a` and `main.a` with `nm`/`ar`, including
exactly one member for each shared kernel and SWCME. A direct
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
- only the current R0–R2/M/B/T production manifest described in the root README.

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
`AMPS_CONFIG`, `SEP_COMMON_DIR`, and `SWCME_DIR`; the active makefile directory
itself must match its source or copied location. This directly protects against the
`build/main: ../Makefile.conf: No such file or directory` failure.
The same fixture also supplies a synthetic enclosing `amps` target and verifies
that `strict-production` delegates to it and audits `build/main` archives rather
than attempting a bare compile in `srcSEP3D`.

#### BLDL3D06 — production turbulence-header boundary

AMPS compiles `build/main/main_lib.cpp` through a generic `Makefile.conf` rule.
Some deployed rules do not consume application additions to `CPPFLAGS`,
`CXXFLAGS`, or `INCLUDE`. This regression therefore compiles
`turbulence_models.h` with only the srcSEP3D include root and requires it to be
independent of `sep_common`. It then compiles the opt-in
`coefficient_bridge.h` with the canonical `SEP_COMMON_DIR` include path. The
test reproduces the former `sep_coefficient_physics.h: No such file or
directory` failure without requiring a configured AMPS build.

### Phase R1 shared-library gates

`ARCH3D02` runs both canonical archives' `verify` targets, compares exact `ar`
membership, and checks that the srcSEP3D makefile consumes the canonical object
lists. It deliberately does not inspect or require the independent `srcSEP`
application tree. `SWCME3D01` invokes the relocated SWCME runner and requires
`R1CFG01`, `R1D01`, `R3D01`, and `R13D01` to pass.

```bash
../src/models/swcme/test/run_tests.py --routine \
  --output-dir test_output/swcme-r1 --rebuild
```

### Phase R2 lifecycle gates

| ID | Acceptance contract |
|---|---|
| `LIFE3D01` | canonical legal path reaches `Finalized`; step/output/checkpoint counters remain inside `Runtime` |
| `LIFE3D02` | all 80 operation/state pairs match the transition table; rejected calls preserve state, counters, and snapshot generation |
| `LIFE3D03` | pre-mesh layout is deterministic; output-only settings do not change physics identity; mismatched binding is atomic |
| `LIFE3D04` | standalone Parker and SWMF adapters reach `SnapshotReady` through `Runtime`; coupled sources contain no process/file/environment parser |

```bash
python3 test/run_tests.py --suite r2 --rebuild \
  --output-dir test_output/r2
```

### Phase M mesh/storage gates

| ID | Acceptance contract |
|---|---|
| `MSH3D01` | one million deterministic points stay within the configured cell-size bounds |
| `MSH3D02` | radial surface, transition, and octave values match closed forms |
| `MSH3D03` | generated Parker centrelines have negligible tube distance for both polarities |
| `MSH3D04` | fast tube-distance approximation converges above second order |
| `MSH3D05` | shoulder mesh is 2:1 balanced and a deliberately illegal jump is detected |
| `MSH3D06` | co-rotation preserves resolution |
| `MSH3D07` | five octrees reproduce exact leaf/memory counts and reject non-owner writes |
| `MSH3D08` | Earth and Mars domains enclose exact declared outer spheres |
| `MSH3D09` | mixed coarse/fine gradients are linear-exact and rank-deficient stencils fail |

```bash
python3 test/run_tests.py --suite phase-m --rebuild \
  --output-dir test_output/phase-m
```

### Phase B background/snapshot gates

| IDs | Acceptance contract |
|---|---|
| `BGP3D01–03` | analytic Parker divergence, components, and tangency |
| `BGP3D04–06` | focusing derivative, radial-wind derivatives, and finite polar limits |
| `SNAP3D01–02` | missing/non-finite required fields reject the candidate |
| `SNAP3D03–05` | coupling units, epoch consistency, and atomic failure |
| `SNAP3D06` | compatible snapshots interpolate only inside their time bracket |
| `SNAP3D07–08` | per-sample batch status and explicit coordinate-frame policy |

```bash
python3 test/run_tests.py --suite phase-b --rebuild \
  --output-dir test_output/phase-b
```

### Phase T turbulence/coefficient gates

| ID | Acceptance contract |
|---|---|
| `TUR3D01` | log-space quadrature recovers prescribed magnetic variance |
| `TUR3D02` | AWSoM energies use `deltaB²=mu0*w` and directions follow field polarity |
| `TUR3D03` | proton/electron resonances below, inside, and above the band apply the selected policy exactly |
| `TUR3D04` | incomplete waves fail unless ballistic mode is explicit and typed |
| `COEF3D01` | Dmumu/mean-free-path/kappa conversions round-trip below 1e-12 over six decades |
| `COEF3D02` | srcSEP3D bridge and direct `sep_common` Jokipii calls are bitwise identical |

```bash
python3 test/run_tests.py --suite phase-t --rebuild \
  --output-dir test_output/phase-t
```

## Named suites

| Suite | Contents |
|---|---|
| `standalone` | C++ registry, shell exit-code probes, and `RUN3D01` |
| `r0` | R0 source/ABI/production gates plus RUN3D01, LAY01, and BLD01 |
| `r1` | canonical shared-archive audit, relocated SWCME suite, and frozen common kernels |
| `r2` | LIFE3D01–LIFE3D04 immutable configuration and lifecycle gates |
| `phase-m` | MSH3D01–MSH3D09 mesh/storage gates |
| `phase-b` | BGP3D01–06 and SNAP3D01–08 background/snapshot gates |
| `phase-t` | TUR3D01–04 and COEF3D01–02 turbulence/coefficient gates |
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
| `cannot locate sep_common directory` | pass `--sep-common-dir` |
| `cannot locate sep_common.a` | build the shared archive, then pass `--sep-common-archive` |
| `BLDL3D01 SKIP` | run in a configured AMPS checkout or pass `--make-config` |
| `BLDL3D01` reports `pic.h: No such file or directory` from `srcSEP3D` | update the makefile; the gate must delegate to top-level `make amps`, not compile `main_lib.cpp` directly in the source directory |
| `BLDL3D02` reports stale `SEP3D.cpp` | remove the obsolete file from the installed `srcSEP3D`; overlay extraction does not delete files left by an older version |
| `BLDL3D03 SKIP` | pass `--amps-source` pointing to a tree containing `src/pic/pic.h` |
| `build/main` cannot find `../Makefile.conf` | run `make print-layout-paths`; `AMPS_CONFIG` must resolve to the absolute `AMPS/Makefile.conf` path |
| standalone compile failure | rerun with `--rebuild --verbose` |
| unknown test/group | use `--list`; unknown selectors are usage errors |
| report missing after a C++ test | treat as ERROR; inspect verbose subprocess output |
