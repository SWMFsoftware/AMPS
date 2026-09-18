# srcSEP3D Testing Procedure

`test/run_tests.py` is the single user-facing test interface. Its selectors
match `srcSEP/test/run_tests.py` so the two applications can use the same
automation habits. The runner combines R0/R1/R2 foundation evidence with
Phase-M mesh, Phase-B background, Phase-T turbulence/coefficient, Phase-P
transport, Phase-A adapter, and Phase-O sampling/restart evidence without
requiring AMPS for dependency-free tests. Phase V adds controlled integration
and physics checks plus explicit linked/cross-model/observational evidence
gates. R01–R07 production-runtime improvements are native registry tests, not
an external checklist.

## Quick commands

```bash
python3 test/run_tests.py --list
python3 test/run_tests.py --routine --amps-source /path/to/AMPS
python3 test/run_tests.py --test LAY01
python3 test/run_tests.py --suite r1 --suite r2
python3 test/run_tests.py --suite improvements-c --rebuild
python3 test/run_tests.py --suite improvements-r --rebuild
python3 test/run_tests.py --suite improvements-v --rebuild
python3 test/run_tests.py --test V2D01 --sep1d-source /path/to/AMPS/srcSEP
python3 test/run_tests.py --suite phase-m --suite phase-b --suite phase-t
python3 test/run_tests.py --suite phase-p --suite phase-a --suite phase-o
python3 test/run_tests.py --suite phase-v --rebuild
python3 test/run_tests.py --suite phase-v --amps /path/to/amps \
  --validation-data /path/to/evidence \
  --validation-launch-prefix "mpiexec -n 8"
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
| `--amps PATH` | supplies the configured linked executable for Phase-V native cases |
| `--timeout SEC` | applies a per-command timeout |

Additional setup options are `--amps-source`, `--make-config`,
`--sep-common-dir`, `--sep-common-archive`, `--validation-data`, and
`--validation-launch-prefix`.

## Evidence classes

### Standalone C++ registry

`test/stage1` is compiled with no AMPS include path and no MPI library. It links
the Phase-M mesh model, Phase-B providers/snapshots, Phase-T turbulence bridge,
Phase-P transport cores, Phase-A neutral adapters, Phase-O output/restart,
Phase-V metrics/audits, R2 Runtime, test callbacks, and shared `sep_common.a`.
The runner adds
`-Wall -Wextra -Wpedantic -Werror`.

| Group | IDs | Purpose |
|---|---|---|
| `HARN` | `HARN01`–`HARN04` | registry selection, PASS/FAIL/SKIP/ERROR exits, JSON, and JUnit |
| `LAY` | `LAY01`, `LAY02` | AMPS-free core/background/runtime/mesh/turbulence/transport/adapters/output/validation rule and negative control |
| `BLD` | `BLD01` | `nm -u` confirms the standalone binary has no AMPS/MPI symbols |
| `UTIL` | `UTIL02` | byte-exact shared-kernel reference record |
| `LIFE3D` | `LIFE3D01`–`LIFE3D04` | immutable configuration, state machine, frozen layout, counters, adapter parity, and no-parser boundary |
| `R3D` | `R3D01`–`R3D07` | mover hook, requested-time loop, snapshot transaction, tick/events, source, observers, restart |
| `CFG3D` | `CFG3D01`–`CFG3D05` | C01-C05 schema/CLI, typed contracts, domains, Parker geometry, and mesh/memory preflight |
| `MSH3D` | `MSH3D01`–`MSH3D09` | resolution, tube geometry, balance, octree budget/ownership, presets, gradients |
| `BGP3D` | `BGP3D01`–`BGP3D06` | analytic Parker field/plasma identities and polar limits |
| `SNAP3D` | `SNAP3D01`–`SNAP3D08` | snapshot completeness, coupling conversion, atomicity, interpolation, batch/frame policy |
| `TUR3D` | `TUR3D01`–`TUR3D04` | spectrum, AWSoM convention, resonance, missing-data policy |
| `COEF3D` | `COEF3D01`–`COEF3D05` | conversion/shared identity plus tensor assembly, Itô drift, and rejection |
| `PRK3D` | `PRK3D01`–`PRK3D08` | Parker moments, characteristics, PDE/first passage, and named limits |
| `FTE3D` | `FTE3D01`–`FTE3D07` | focused streaming, focusing, pitch scattering/boundaries, momentum, strong-scattering limit |
| `RNG3D` | `RNG3D01`–`RNG3D03` | worker/order independence and random-purpose isolation |
| `V1D` | `V1D01`–`V1D05` | V01 tensor, cross-field moments, drift, focused invariance, and timestep |
| `V2D` | `V2D01` | distinct compiled srcSEP/srcSEP3D production-core parity |
| `V5D` | `V5D01` | native profiles, deferred-R8 campaign block, and release governance |
| `ADP3D` | `ADP3D01` | exact two-core production registry and one validating dispatch |
| `NAT3D` | `NAT3D04`–`NAT3D08` | boundary outcomes, ledger closure, sampling isolation, output schema, shock crossing |
| `SHK3D` | `SHK3D01`–`SHK3D04` | common source identity, moving-sphere geometry, guards, and normalization |
| `RST3D` | `RST3D01`–`RST3D03` | full round trip, transactional rejection, and snapshot policy |
| `INT3D` | `INT3D01`–`INT3D03` | stable-ID rank merge, global conservation, and resource budgets |
| `VFY3D` | `VFY3D01`–`VFY3D05` | comparison metrics and analytical Parker/focused/SWCME validation |
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
- only the current R0–R2/M/B/T/P/A/O production manifest described in the root README.

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
than attempting a bare compile in `srcSEP3D`. The synthetic build uses valid
ELF fixture objects, creates the complete application/shared member manifest,
and exports the normalized-domain and C05 mesh-preflight ABI names required by
the production audit. It therefore exercises the strict archive contract while
leaving real AMPS compilation coverage to BLDL3D01.

#### BLDL3D06 — production transitive-header boundary

Every generic AMPS translation unit reaches `SEP3D.h` through generated
`pic.h`. Those compiler commands do not necessarily contain the canonical
`sep_common` or SWCME include paths. This regression rejects concrete
source/restart/particle-adapter includes from the umbrella header and rejects a
concrete `source_runtime.h` include from the AMPS adapter declaration. The
corresponding types must remain forward declarations until an implementation
file that is compiled with model include paths.

The same gate compiles `turbulence_models.h` with only the srcSEP3D include
root and requires it to be independent of `sep_common`. It then compiles the opt-in
`coefficient_bridge.h` with the canonical `SEP_COMMON_DIR` include path.

Finally, the gate creates a mock installed `Makefile.conf` whose generic
`main_lib.o` recipe intentionally ignores `CPPFLAGS`, `CXXFLAGS`, and
`INCLUDE`. The command has the same fixed shape seen in production AMPS logs.
It must nevertheless compile a translation unit that includes both
`sep_injection_spectrum.h` and `swcme_sep_source.hpp`, proving that the
makefile's target-scoped `CPLUS_INCLUDE_PATH` reaches `mpicxx`. The test thus
reproduces both former `sep_coefficient_physics.h` and
`sep_injection_spectrum.h: No such file or directory` failures without a
configured AMPS build.

#### BLDL3D07 — application-object ABI freshness

Reproducible source packages normalize file timestamps. If such a package is
overlaid on an existing configured tree, timestamp-only dependency checks can
mistake an older `build/main/mesh/mesh_model.o` for a current object. The
symptom is a final-link failure for the new normalized-domain or C05 preflight
symbols even though other, unchanged mesh symbols resolve.

The production makefile therefore forces every srcSEP3D-owned object to be
recompiled whenever its `lib` or `amps` target runs. It recreates indexed
archives, verifies exactly one copy of every application/shared member, and
uses `nm -C` to require the current `MakeDomain(RunConfiguration3DOptions)` and
`BuildRefinementPreflight` definitions before AMPS reaches `mpif90`. BLDL3D07
protects that makefile contract without requiring a configured host.

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

### R01–R07 production-runtime gates

| ID | Acceptance contract |
|---|---|
| `R3D01` | generated `picGlobal.dfn` hook names the exact mover signature and strict production depends on hook installation |
| `R3D02` | one AMPS request is consumed by multiple accepted substeps with fresh local resolution and exact final time |
| `R3D03` | failed fill preserves the active generation; a collectively accepted staged pair commits atomically |
| `R3D04` | integer clock/event schedule agrees with host time and restores without cadence drift |
| `R3D05` | physical source number survives rounding/cap, policies are counted, and successive ticks use distinct identities |
| `R3D06` | moving observer geometry, acceptance/uncertainty metadata, and commit-only accumulator reset are enforced |
| `R3D07` | schema-2 restart round-trips identities/layout, clocks/events, providers, shock, RNG tuple, ledgers, and sampling state |

```bash
python3 test/run_tests.py --suite improvements-r --rebuild \
  --output-dir test_output/improvements-r
```

These tests are AMPS-independent and always participate in `--all`.
`BLDL3D01/03/05` remain the configured-host evidence for the generated mover
ABI and copied `build/main` production layout. On a configured tree the build
sequence is:

```bash
make -C srcSEP3D prepare-production
env MAKEFLAGS="-j16" srcSEP3D/test/run_tests.py --all \
  --amps-source . --make-config Makefile.conf \
  --output-dir srcSEP3D/test_output/all --rebuild
```

### C01-C05 configuration and preflight gates

| ID | Acceptance contract |
|---|---|
| `CFG3D01` | complete versioned input parses; file and typed construction have one fingerprint; documented CLI, early-error, and dry-run contracts hold |
| `CFG3D02` | output-only changes preserve physics identity; physical changes alter it; invalid shock/source intent fails; SWMF uses the same factory |
| `CFG3D03` | solar, one-AU, Mars, and explicit radii normalize exactly; invalid observers fail; boundary status respects crossing direction |
| `CFG3D04` | mesh centerline/tangent and analytic field use one Parker geometry; polarity reverses `B` without moving the tube |
| `CFG3D05` | composite profiles are monotone, tube width scales from its reference, all memory categories/levels report, and an impossible level cap fails |

```bash
python3 test/run_tests.py --suite improvements-c --rebuild \
  --output-dir test_output/improvements-c
```

### Phase M mesh/storage gates

| ID | Acceptance contract |
|---|---|
| `MSH3D01` | one million deterministic points stay within the configured cell-size bounds |
| `MSH3D02` | linear radial surface, midpoint, and transition values match closed forms |
| `MSH3D03` | generated polarity-independent Parker centreline has negligible tube distance |
| `MSH3D04` | fast tube-distance approximation converges above second order |
| `MSH3D05` | composite tube mesh is 2:1 balanced and a deliberately illegal jump is detected |
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

`phase-t` deliberately contains `COEF3D01–02`; the later tensor/drift
coefficient tests belong to Phase P even though they share the `COEF3D` group.

### Phase P transport gates

| IDs | Acceptance contract |
|---|---|
| `COEF3D03–05` | rank-one tensor assembly, complete numerical/analytic Itô divergence, invalid/reserved coefficient rejection |
| `PRK3D01–08` | diffusion moments, advection/rotation, nonuniform equilibrium, cooling, first passage, radial PDE, all named limits |
| `FTE3D01–07` | ballistic characteristic, focusing, eigenmodes, reflecting boundaries, momentum, Parker limit, zero-perpendicular guard |
| `RNG3D01–03` | worker partition, list order, and future-purpose changes cannot alter keyed histories |

```bash
python3 test/run_tests.py --suite phase-p --rebuild \
  --output-dir test_output/phase-p
```

### Phase A mover/source adapter gates

| ID | Acceptance contract |
|---|---|
| `ADP3D01` | registry contains only Parker tensor and focused split; both pass one validator |
| `NAT3D04` | inner absorption, outer escape, and invalid background remain distinct |
| `NAT3D05` | integer particle ledger closes exactly and mismatch is transactional |
| `NAT3D08` | first expanding-sphere root is recorded once per generation |
| `SHK3D01` | identical common SWCME records produce bitwise-identical dimensional source samples |
| `SHK3D02–04` | analytic moving shock, ownership/generation guards, and event-weight sum |

```bash
python3 test/run_tests.py --suite phase-a --rebuild \
  --output-dir test_output/phase-a
```

The standalone gate proves the AMPS-independent conversion and conservation
logic. `BLDL3D01` remains required to compile `amps_particle_adapter.cpp`
against the real particle-buffer/list ABI and the generated mover macro.

### Phase O sampling/output/restart gates

| ID | Acceptance contract |
|---|---|
| `NAT3D06` | stable-ID order and repeated sampling are bitwise identical; input remains untouched |
| `NAT3D07` | atomic output bundle has SI headers, complete identity manifest, and verified hashes |
| `RST3D01` | all state round-trips and future keyed normal draws are identical |
| `RST3D02` | fingerprint/checksum errors leave destination state unchanged |
| `RST3D03` | unavailable snapshot is rejected or awaited only through explicit bounded policy |

```bash
python3 test/run_tests.py --suite phase-o --rebuild \
  --output-dir test_output/phase-o
```

### Phase V integration/scientific-validation gates

The Phase-V suite combines immediately executable prerequisites with external
release evidence; their statuses must be interpreted separately.

| IDs | Acceptance contract |
|---|---|
| `INT3D01–03` | rank/order-independent stable-ID gather, exact global integer conservation, and explicit load/wall/memory budgets |
| `VFY3D01–02` | normalization, log-space metrics, coverage rejection, and malformed-evidence classification |
| `VFY3D03` | 3-D Parker projections reproduce the independent one-dimensional Gaussian Green function |
| `VFY3D04` | focused transport converges at second order to the exact focusing characteristic |
| `VFY3D05` | sampled SWCME/DSA momentum CDF and total represented event weight agree with their declared laws |
| `NAT3D01–03/09–12` | configured AMPS mesh, storage, gradients, balance, budgets, coupled cadence, and output grammar |
| `MPI3D01–02` | multi-rank sampling and restart continuation are decomposition independent |
| `XM3D01–06` | checksum-owned cross-model profiles, exact source identity, convergence, and independent PDE evidence |
| `OV3D01–04` | reviewed event comparisons, with release-gate and diagnostic roles retained in reports |
| `VALRUN3D01` | runner lists all classes, preserves SKIP, verifies checksums, and evaluates convergence bundles |

```bash
# Runs controlled prerequisites now and records unavailable external work as SKIP.
python3 test/run_tests.py --suite phase-v --rebuild \
  --output-dir test_output/phase-v

# Adds independently exported cross-model/observational evidence.
python3 test/run_tests.py --suite phase-v \
  --validation-data /path/to/evidence \
  --output-dir test_output/phase-v-evidence

# Adds native tests from a configured linked application.
python3 test/run_tests.py --suite phase-v --amps /path/to/amps \
  --validation-launch-prefix "mpiexec -n 8" \
  --output-dir test_output/phase-v-linked
```

The launch prefix is parsed into process arguments and is never passed to a
shell. A linked binary must advertise the requested ID through `--list-tests`;
a stale binary is `ERROR`. Scientific evidence must use the templates under
`validation/templates/`, remain inside `EVIDENCE_ROOT/CASE_ID`, and match every
declared SHA-256. Missing evidence is `SKIP`; checksum/schema/provenance failure
is `ERROR`; a valid metric outside tolerance is `FAIL`.

`XM3D03`, `OV3D03`, and `OV3D04` remain diagnostic. V01 now supplies controlled
perpendicular diffusion and guiding-centre drift, but these cases still need
reviewed evidence before cross-field transport can be called scientifically
validated. See `INTEGRATION_SCIENTIFIC_VALIDATION.md` for the full
equations, algorithms, case roles, and evidence schemas.

## Named suites

| Suite | Contents |
|---|---|
| `standalone` | C++ registry, shell exit-code probes, `RUN3D01`, and `VALRUN3D01` |
| `r0` | R0 source/ABI/production gates plus RUN3D01, LAY01, and BLD01 |
| `r1` | canonical shared-archive audit, relocated SWCME suite, and frozen common kernels |
| `r2` | LIFE3D01–LIFE3D04 immutable configuration and lifecycle gates |
| `improvements-c` | CFG3D01–CFG3D05 production configuration and preflight gates |
| `improvements-r` | R3D01–R3D07 production runtime integration gates |
| `improvements-v` | V1D01–05 controlled physics, V2D01 true parity, and V5D01 governance |
| `phase-m` | MSH3D01–MSH3D09 mesh/storage gates |
| `phase-b` | BGP3D01–06 and SNAP3D01–08 background/snapshot gates |
| `phase-t` | TUR3D01–04 and COEF3D01–02 turbulence/coefficient gates |
| `phase-p` | COEF3D03–05, PRK3D01–08, FTE3D01–07, RNG3D01–03 |
| `phase-a` | ADP3D01, NAT3D04–05/08, SHK3D01–04 |
| `phase-o` | NAT3D06–07 and RST3D01–03 |
| `phase-v` | INT3D/VFY3D prerequisites, external NAT3D/MPI3D/XM3D/OV3D cases, and VALRUN3D01 |
| `production` | BLDL3D01–07 |

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
| `main_lib.cpp` reports `sep_injection_spectrum.h: No such file or directory` | install the updated srcSEP3D makefile in the source tree and refresh the copied `build/main`; `BLDL3D06` verifies that the fixed generic recipe receives the target-scoped canonical model search path |
| final link reports undefined `Mesh::MakeDomain(RunConfiguration3DOptions)` or `BuildRefinementPreflight` | stale pre-C03/C05 `mesh_model.o`; install the updated makefile, run `make clean`, and rebuild. BLDL3D07 prevents recurrence |
| standalone compile failure | rerun with `--rebuild --verbose` |
| Phase-V linked case `SKIP` | supply `--amps`; use `--validation-launch-prefix` when MPI launch arguments are required |
| XM3D/OV3D case `SKIP` | supply `--validation-data` containing `CASE_ID/manifest.json` and its declared artifacts |
| Phase-V checksum `ERROR` | regenerate the SHA-256 only after reviewing the changed evidence; never edit a hash merely to silence the gate |
| linked executable does not advertise a case | rebuild the configured AMPS application from this source tree and verify its production test registry |
| unknown test/group | use `--list`; unknown selectors are usage errors |
| report missing after a C++ test | treat as ERROR; inspect verbose subprocess output |
