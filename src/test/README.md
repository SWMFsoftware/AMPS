# AMPS core unit and module tests

`src/test` is the global test-control layer for production code under `src`
and for explicitly registered application-local campaigns. It does not
replace application validation or model-specific qualification. It gives
those tests one portable registry, one result vocabulary, and one
nightly-compatible evidence format while preserving each application's own
oracles and acceptance boundaries.

The initial migration is intentionally honest: every core component is in
`component_inventory.json`, but only entries marked `covered` have complete
registered coverage. `partial` and `planned` are visible technical work, not
implicit passes. The strict coverage target remains red until that work is
actually complete.

The frozen T1 qualification registry contains 11 executable suites and 98
individual cases, including 31 cases in the configured native suite.  Its
retained normal, reverse, and seed-1729 shuffled campaigns each report 93
PASS, five typed provider SKIPs, and no FAIL/ERROR; the retained native
executable has SHA-256
`17b98e5ecd22915fc181cce51b8d71bd6abef517802a456ff2efba659c581919`.
Optional interface providers remain unqualified where their external
libraries and independent fixtures are not configured.

At the frozen T2 checkpoint the registry had 16 suites and 152 cases, without
rewriting the T1 receipt. Six exclusive native owners formed the 85-case `profile`
campaign: the 43-case packed/fixed baseline, 13-case packed individual-weight
split/merge, 7-case aligned-16 fixed layout, 13-case aligned-16 individual
layout/split-merge, 5-case species-local weight/time, and 4-case single-global
weight/time profiles. `/tmp/amps-t2-native-full-normal` passes 43/43; both
`/tmp/amps-t2-profile-reverse` and
`/tmp/amps-t2-profile-shuffle-1729` pass 85/85. These receipts qualify the
declared one-rank host contracts for record layout/lifecycle, finite-positive
weights/time steps, local rescaling, and conservative split/merge. CUDA/device
diagnostics, semantic averaging of opaque application fields, unequal-weight
collision physics, injection/sampling normalization, restart, and distributed
ownership remain open under the T2 contract.

T3 extends the live registry, without rewriting the T1 or T2 receipts, to 17
suites and 210 cases. `CORE-MESH-NATIVE-001` owns 58 one-rank cases: 14
topology/adaptation cases for each of `DIM=1`, `DIM=2`, and `DIM=3`, eight 3-D
interpolation cases, and eight 3-D cut-cell/NASTRAN cases. The retained normal,
reverse, and seed-1729 shuffled receipts at `/tmp/amps-t3-native-normal`,
`/tmp/amps-t3-native-reverse`, and `/tmp/amps-t3-native-shuffle-1729` each pass
58/58, but an independent audit subsequently corrected a lower-dimensional
bounds defect. Those receipts remain superseded bring-up evidence. The final
corrected-tree receipts `/tmp/amps-t3-qualified-normal`,
`/tmp/amps-t3-qualified-reverse`, and
`/tmp/amps-t3-qualified-shuffle-1729` each pass 58/58 and qualify the declared
Cartesian, one-rank contracts. Restart, MPI ownership/ghosts,
non-Cartesian coordinates, deep imbalance, and general production-surface
topology remain open. The legacy `OutputLayerNodeMap` path remains hard-coded
for 3-D and is not exercised or qualified in the 1-D/2-D profiles.

The implementation includes lower-dimensional loop/type/neighbor hardening,
dimension-bounded global initialization and `buildMesh`, and a private
production-`buildMesh` topology oracle in every dimension; split/coarsen
restriction and prolongation; complete coarsest-neighbor
selection for production interpolation, including the coarse-surrounded-by-
fine `<=` selector and corrected edge-1/edge-3 predicates; open-interior
cut-surface contact and `O(A*h)` integration-depth control; translated/scaled
sphere volume recursion; and sparse-ID NASTRAN read/write/read validation.
Details and independent oracles are in `test/srcCoreMeshTests/README.md`.

T4/T5 extend the live registry to 22 suites and 261 cases. The baseline native
owner now contains 87 cases, with 51 selected T4/T5 mover, boundary, injection,
sampling, collision, and IDF cases passing in normal, reverse, and seed-1729
shuffled order. Five new compile-profile owners exercise periodic boundaries,
VHS, VSS, and `Zrot=2/10`; the packed individual-weight owner adds two
selector-specific collision cases. The six-profile variant campaign passes
7/7 selected cases. These receipts establish the bounded contracts documented
in `test/srcCoreUnitTests/README.md`; they do not qualify the remaining
coarse/fine and MPI mover paths, every production surface type, source
ownership/rate normalization, histograms, Yinsi individual weights, an
unequal cross-species Zrot law, or VT/VV relaxation.

The final cross-stage receipt `/tmp/amps-t2-t5-profile-qualified` passes all
194 cases in all 12 configured T2-T5 suites. It supersedes the first
`/tmp/amps-t2-t5-profile-final` diagnostic, which exposed and retained one
stale seven-species oracle after the charged profile added H+ and electrons.

T8/T9 extend the live registry to 41 suites and 359 cases. T8 owns 34 cases:
five native manufactured ECSIM operator/moment/identity checks and 29
dependency-light boundary, deposition, field, field-line, geometry,
background and callback/output checks. T9 owns 13 versioned persistence/data/
output cases plus one native AMPS sampling/particle restart case selected from
the baseline owner. The native restart restores AMR cell placement, particle
state and extension bytes, global weights/time steps, completed sampling data
and the scalar MPI RNG seed. `STATUS.md` records the exact receipts and the
explicit remaining legacy-format, hybrid-RNG and distributed boundaries.

## Architecture and ownership

The execution path is deliberately one-way:

```text
manifest + machine profile -> global runner -> per-suite context
                           <- atomic result <- individual runner
                           -> JSON/JUnit aggregate
```

`manifest.json` owns what may run. `component_inventory.json` owns the honest
coverage state. The global runner owns selection, prerequisite probes,
timeouts, work/artifact paths, protocol validation, and aggregate reports. An
individual runner owns only the mechanics and assertions for one registered
suite. This boundary prevents a suite from selecting extra tests, changing the
global output location, or declaring itself green while returning a failed
case.

[`TEST_CATALOG.md`](TEST_CATALOG.md) is the generated human-readable projection
of the expanded manifest: all suites and cases with concise what/how/why/where
descriptions. `generate_test_catalog.py --check`, structural traceability, and
the harness reject catalog drift. The catalog documents registered intent; the
JSON receipt from a particular run remains authoritative for what actually
executed and how it ended.

Every suite has two writable locations. `work_dir` is disposable compiler and
execution state and may be placed on node-local scratch. `artifact_dir` is
retained evidence such as a comparison table or executable transcript. The
suite publishes `result.json` atomically in its work directory; the global
runner validates it before copying its meaning into the aggregate. A killed
process, truncated JSON file, mismatched ID, duplicate case, hidden failure,
or status/exit-code disagreement is therefore reported as `ERROR`.

Three governance ledgers refine that ownership without changing execution:

- `source_traceability.json` assigns every current production translation unit
  and public `.h`, `.hpp`, `.dfn`, or `.inc` contract to exactly one area with
  `tested`, `partial`, `planned`, or `exempt` status, stable owning cases or a
  model-local registry, a reason, and explicit remaining gaps.
- `configuration_contracts.json` records machine-profile, runtime-input,
  `Config.pl`, generated-file, compile-profile, and runtime-mode contracts.
  A successful build is never by itself evidence that a mode is tested.
- `t1_public_algorithm_inventory.json` enumerates the public general,
  species, units, generated-species, ABI, Geopack, and provider families in
  the T1 scope.  It distinguishes `tested`, `partial`, `unavailable`, and
  `planned` behavior so an aggregate case count cannot hide an untested
  public family.

Validate and optionally expand the rule-based ledger from the repository root:

```sh
python3 src/test/validate_traceability.py \
  --output /tmp/amps-core-traceability.json
```

The dependency-free validator rejects an unclassified production file,
overlapping ownership, a missing authority/local registry, stale suite or case
ownership, component-ledger disagreement, and invalid status/gap combinations.
The expanded output lists every discovered file and its reviewed status. It is
a structural receipt, not a physics-test result.

## Run

From the AMPS repository root:

```sh
PROFILE=test/qualification-machine/passion.json
python3 src/test/run_tests.py --config "$PROFILE" --campaign smoke
python3 src/test/run_tests.py --config "$PROFILE" --campaign nightly
python3 src/test/run_tests.py --config "$PROFILE" --campaign profile
python3 src/test/run_tests.py --config "$PROFILE" --campaign nightly --jobs 4
python3 src/test/run_tests.py --config "$PROFILE" --all --streams 3 --jobs 2
python3 src/test/run_tests.py --config "$PROFILE" --list --streams 3
python3 src/test/run_tests.py --config "$PROFILE" --full-validation --list
python3 src/test/run_tests.py --config "$PROFILE" --coverage
make -C src/test smoke
make -C src/test STREAMS=3 JOBS=2 nightly
make core_test_nightly
```

`--config` is mandatory, including for `--list` and `--coverage`. Omitting it
terminates before manifest selection and points to
`test/qualification-machine/`, which contains the maintained
[`passion.json`](../../test/qualification-machine/passion.json), a
[`template.json`](../../test/qualification-machine/template.json), and profile-creation
instructions. The make wrappers select
`../../test/qualification-machine/$(hostname -s).json` by default; override with
`make -C src/test CONFIG=/path/to/machine.json nightly` when needed.

`--streams N` is the maximum number of independent execution streams and
defaults to `1`. The complete registry has three streams: `core`,
`application-srcEnceladus`, and `application-srcMoon`. `--jobs N` separately
limits concurrent nonexclusive suites inside each active stream. For example,
`--streams 3 --jobs 2` permits all three streams and up to two safe suites per
stream to make progress; it does not permit native builds to race. The make
wrapper exposes the same controls as `STREAMS=N JOBS=N`.

A bare `--list` is a complete registry inventory. It first summarizes every
registered application, then prints the core and each selected application in
separate `STREAM` sections with exact suite/case counts. Explicit selectors
such as `--campaign smoke --list`, `--component srcMoon --list`, or
`--test ID --list` still narrow the selected stream plan while retaining the
registered-application summary.

`--tier` remains a backward-compatible alias for `--campaign`. The top-level
`test_all` target depends on `core_test_nightly`. The live registry has 54
suites and 1,135 cases; use `--campaign NAME --list` for an exact campaign
membership or bare `--list` for the complete registry. The
routine `nightly` campaign intentionally remains the established 42-suite,
375-case deterministic core selection, while `extended` owns the T11 model
imports and bounded T12 portability matrix. The `distributed` campaign selects
the five-case two-rank channel/RNG suite and the 16-case T10 collective
1/2/4-rank matrix (two suites and 21 cases total).
`--all` ignores campaign membership and selects the complete registry.
`--full-validation` also selects the complete registry, but additionally marks
the execution as an explicit request for expensive and qualification modes:
the Enceladus E3 owner attempts its 10,000-step run and the Moon aggregate uses
`--campaign-mode qualification`. The flag is intentionally incompatible with
`--campaign`, `--tier`, `--all`, `--test`, `--match`, `--component`, and
`--case`; a full-validation receipt can therefore never describe a narrowed
selection. Use a machine overlay to provide actual SPICE paths, compilers,
launchers, external datasets, and other site prerequisites. The explicit flag
authorizes costly execution but does not fabricate missing resources or turn a
typed FAIL/SKIP/ERROR into PASS.

Run the complete execution only on an appropriately provisioned worker:

```sh
python3 src/test/run_tests.py \
  --config test/qualification-machine/passion.json \
  --full-validation --streams 3 --jobs 2
```

Both `selection.json` and the final `amps-core-tests.json` record
`full_validation: true` and the union of included campaign names. The Moon M0
owner remains an explicit SKIP until its historical multi-build procedure is
converted to the mandatory native-build preflight; full validation exposes
that open gate rather than bypassing it.
Site/nightly orchestration must preserve the native suite's exclusive barrier
because it selects a repository-wide application.

The `application-smoke` campaign exposes both application trees requested for
global orchestration:

```sh
PROFILE=test/qualification-machine/passion.json
python3 src/test/run_tests.py --config "$PROFILE" --campaign application-smoke
python3 src/test/run_tests.py --config "$PROFILE" --test APPLICATION-ENCELADUS-FORWARD-001
python3 src/test/run_tests.py --config "$PROFILE" --test APPLICATION-ENCELADUS-NATIVE-STARTUP-001
python3 src/test/run_tests.py --config "$PROFILE" --test APPLICATION-ENCELADUS-E3-001
python3 src/test/run_tests.py --config "$PROFILE" --component srcMoon --list
```

`APPLICATION-ENCELADUS-FORWARD-001` executes 42 dependency-light checks.
`APPLICATION-ENCELADUS-NATIVE-STARTUP-001` performs the mandatory clean native
build and imports all six pre-MPI startup/configuration gates. The separate
13-case `APPLICATION-ENCELADUS-E3-001` owner is in `extended`, not the routine
application-smoke campaign: it reports every identity as `SKIP` unless a
machine overlay enables `enceladus-e3` or the operator explicitly requests
`--full-validation`. A capability skip is visibility, not validation.
`APPLICATION-MOON-REGRESSION-001` owns all U01-U26/I11-I20 outcomes but is a
prerequisite SKIP unless the machine profile declares `moon-spice` and supplies
both SPICE paths. `APPLICATION-MOON-M0-001` lists I01-I10 as typed SKIPs until
the local multi-build orchestrator satisfies the current mandatory native
preflight. These statuses are intentional and cannot be read as application
qualification.

Run one suite or component:

```sh
PROFILE=test/qualification-machine/passion.json
python3 src/test/run_tests.py --config "$PROFILE" --test CORE-GENERAL-INPUT-001
python3 src/test/run_tests.py --config "$PROFILE" --test CORE-GENERAL-MPI-001
python3 src/test/run_tests.py --config "$PROFILE" --test CORE-INTERFACE-ABI-001
python3 src/test/run_tests.py --config "$PROFILE" --component src/general
python3 src/test/run_tests.py --config "$PROFILE" --match 'CORE-GENERAL-*'
```

Run exact cases or deterministic isolation orders without broadening the
selection:

```sh
PROFILE=test/qualification-machine/passion.json
python3 src/test/run_tests.py --config "$PROFILE" --case INP-03
python3 src/test/run_tests.py --config "$PROFILE" --case PHYS-IDF-RT-001 --order reverse
python3 src/test/run_tests.py --config "$PROFILE" --all --order shuffle --shuffle-seed 1729
```

`--case` may be repeated. A requested ID must have exactly one manifest owner;
the runner forwards only those selected IDs to the owning adapter and validates
their ordered result list. Reverse/shuffle support is an execution mechanism,
not proof of isolation by itself; the retained normal/reverse/shuffle and
one-case evidence is recorded in `STATUS.md`.

Each configured AMPS application test is exclusive and is selected by its
declared campaigns. The baseline profile can be run directly with:

```sh
python3 src/test/run_tests.py \
  --config test/qualification-machine/passion.json \
  --test CORE-PIC-NATIVE-001
```

Run the complete T3 dimensional matrix and its isolation orders with:

```sh
python3 src/test/run_tests.py --config test/qualification-machine/passion.json \
  --test CORE-MESH-NATIVE-001 --order normal \
  --output-dir /tmp/amps-t3-qualified-normal
python3 src/test/run_tests.py --config test/qualification-machine/passion.json \
  --test CORE-MESH-NATIVE-001 --order reverse \
  --output-dir /tmp/amps-t3-qualified-reverse
python3 src/test/run_tests.py --config test/qualification-machine/passion.json \
  --test CORE-MESH-NATIVE-001 --order shuffle \
  --shuffle-seed 1729 --output-dir /tmp/amps-t3-qualified-shuffle-1729
```

The retained executable SHA-256 values are order-independent: 1-D
`54e79ddf1c34e2910a059a23e38232bea6322ab7a1ee6d31e35204acd8dbf4e0`, 2-D
`8d5311b99bfba3812599336b2cdbb4fbcd53ab16a3e625122b6db44a30b4aad9`, and
3-D `96f13a4be39ed221b7a252bb8313db3371ed49b2fd464ea1b5d09a545bf42e34`.
The retained T2+T3 checkpoint command
`python3 src/test/run_tests.py --campaign profile --order normal --output-dir
/tmp/amps-t2-t3-profile-qualified` passes 143/143 cases in all seven suites
(the six T2 owners plus T3), with no FAIL, SKIP, or ERROR. That receipt was
created when `profile` contained only those seven owners; the same command now
selects the current 17-suite/211-case profile, so use explicit suite IDs to
reproduce only the historical T2+T3 subset.

Run the bounded T4/T5 baseline and compile variants with:

```sh
python3 src/test/run_tests.py --config test/qualification-machine/passion.json \
  --test CORE-PIC-MOVER-CHARGED-001 \
  --test CORE-PIC-MOVER-RELATIVISTIC-001 --test CORE-PIC-GYRO-001 \
  --test CORE-PIC-BOUNDARY-001 --test CORE-PIC-INJECTION-001 \
  --test CORE-PIC-SAMPLING-001 --test CORE-PIC-COLLISION-ALGORITHMS-001 \
  --test CORE-PIC-IDF-001

python3 src/test/run_tests.py --config test/qualification-machine/passion.json \
  --test CORE-PIC-PERIODIC-001 \
  --test CORE-PIC-COLLISION-VHS-001 --test CORE-PIC-COLLISION-VSS-001 \
  --test CORE-PIC-IDF-ROTZ-Z2-001 --test CORE-PIC-IDF-ROTZ-Z10-001 \
  --test CORE-PIC-COLLISION-INDIVIDUAL-WEIGHT-001
```

Run the bounded T6 reaction and chemistry gate with:

```sh
python3 src/test/run_tests.py --config test/qualification-machine/passion.json \
  --test CORE-CHEMICAL-SOLVER-001 \
  --test CORE-PHOTOCHEM-EXECUTION-001 \
  --test CORE-PHOTOCHEM-DAUGHTER-001 \
  --test CORE-CHARGE-EXCHANGE-001 \
  --test CORE-ELECTRON-IMPACT-001
```

The charge-exchange and electron-impact IDs are selection aliases for one
exclusive configured reaction-profile owner. The photochemistry execution ID
is the baseline native alias; the daughter ID owns a separate
individual-weight build. The solver is dependency-light and may run
concurrently with non-native suites, but all configured owners remain
repository-wide scheduler barriers.

Run the implemented T7 surface, exosphere, sputtering and dust slice with:

```sh
python3 src/test/run_tests.py --config test/qualification-machine/passion.json \
  --test CORE-SURFACE-001 \
  --test CORE-EXOSPHERE-SOURCES-001 \
  --test CORE-EXOSPHERE-SPHERE-NATIVE-001 \
  --test CORE-SPUTTERING-001 --test CORE-DUST-001 \
  --test CORE-DUST-NATIVE-001 \
  --test CORE-DUST-INJECTION-NATIVE-001
```

`CORE-SURFACE-001` and `CORE-EXOSPHERE-SOURCES-001` are selection aliases for
one dependency-light owner and therefore compile once. The three configured
T7 owners are exclusive: they perform the mandatory clean
`test/core-exosphere-sphere`, `test/model-dust`, and
`test/model-dust-injection` rebuilds and run one MPI rank. The command selects
29 cases in six canonical suites; the retained campaign passes all 29 with no
FAIL/SKIP/ERROR. `STATUS.md` records the exact aggregate and
normal/reverse/shuffled dependency-kernel receipts plus the remaining MPI,
nonspherical, target-specific source, native sputtered-product insertion, and
restart boundaries.

Run the bounded T8 and T9 owners with:

```sh
python3 src/test/run_tests.py --config test/qualification-machine/passion.json \
  --test CORE-ECSIM-OPERATORS-001
python3 src/test/run_tests.py --config test/qualification-machine/passion.json \
  --test CORE-ECSIM-BOUNDARIES-001 \
  --test CORE-PIC-DEPOSITION-001 --test CORE-PIC-FIELD-SOLVER-001 \
  --test CORE-PIC-AUXILIARY-001 --test CORE-PIC-FIELD-LINE-001 \
  --test CORE-PIC-CALLBACK-OUTPUT-001 --jobs 4
python3 src/test/run_tests.py --config test/qualification-machine/passion.json \
  --test CORE-PIC-RESTART-001 \
  --test CORE-PIC-DATAFILE-001 --test CORE-PIC-OUTPUT-001 --jobs 3
python3 src/test/run_tests.py --config test/qualification-machine/passion.json \
  --test CORE-PIC-RESTART-NATIVE-001
```

The first and last commands are exclusive configured builds and obey the
mandatory clean-build preflight. `CORE-PIC-RESTART-NATIVE-001` is a selection
alias for one case in `CORE-PIC-NATIVE-001`, so it reuses the established
native application rather than adding a duplicate owner. The dependency-light
T9 fixtures validate canonical continuation and corruption/compatibility; the
native case separately proves the real sampling and particle restart paths.

The historical particle-buffer, packed-fixed layout, physics, species, and
species-global weight/time suite IDs are manifest aliases for selection
compatibility; they are not separate executable suites and must not trigger
two rebuilds. `--campaign profile` currently serializes 18 exclusive owners
and 217 cases. The six T2 owners select mutually exclusive record/weight
profiles; the T3 owner itself performs three clean builds for its 1-D, 2-D,
and 3-D applications. Later configured reaction and T7 owners add their own
clean application selections.

It performs the required native rebuild sequence: validate the repository and
concurrent-process state, remove only `build`, select
`test/core-unit-tests` through `Config.pl`, and run `make amps` with
`MAKEFLAGS=-j16`. It never edits `.amps.conf` directly.

The adapter builds once, retains that executable and its generated
configuration hashes, then launches every selected executable case as a fresh
one-case MPI process. This is the isolation boundary for persistent PIC mesh,
spherical-boundary, callback, and random-generator singletons. The C++
fixtures additionally use the guards in `test/srcCoreUnitTests/core_test_state.h`
to restore the particle buffer, intrusive-list heads, lifetimes, cell sampling
bytes, weights, time steps, and scoped scalar state. Direct multi-case C++
execution remains useful for diagnosis, but nightly qualification does not
depend on every process-global singleton being safely reset in one process.

`--streams N` controls concurrent core/application streams; `--jobs N`
controls suite concurrency inside each stream. Only consecutive nonexclusive
suites within a stream may overlap. An `exclusive: true` suite is also guarded
by a writer-priority checkout gate shared by every stream: it starts only after
all nonexclusive readers finish, runs without another stream inspecting or
mutating generated configuration, and blocks new readers while waiting. This
serializes `Config.pl` and root-build transitions without unnecessarily
serializing independent dependency-light work. Stream and suite limits are
independent, so operators must choose values compatible with machine CPU and
memory capacity. `--fail-fast` lets already-running work publish its receipt
but prevents subsequent batches from starting.

The baseline native suite owns three generated-configuration adapter cases,
three particle-buffer cases, seven packed fixed-layout cases, 19 original
physics cases, 19 T4 cases, 25 T5 cases, four compiled/generated species
cases, two production-linked global-unit cases, and five species-global
particle weight/time-step cases, plus one native sampling/particle restart
case. Ten other native owners cover packed/aligned
individual weights and split/merge, aligned fixed records, species-local and
single-global weight/time, periodic boundaries, VHS/VSS, and two Zrot
profiles. Across the matrix, T2 checks independently computed alignment and
stride, clone/reuse/list/capacity behavior, isolated invalid/exhaustion paths,
finite-positive host values, local physical-distribution rescaling, and
split/merge conservation of count, physical weight, species, weighted
position, momentum, translational and configured internal energies. Opaque
extension bytes use an explicit clone/survivor policy. The physics cases add
analytical/invariant checks for Huebner
rates, branching, products, and exponential parent decay; Maxwellian moments;
hard-sphere scattering; populated-cell NTC collision rates; Larsen-Borgnakke
rotational/vibrational laws and RT equilibration; neutral Simple/Boris motion;
and analytic spherical-surface reflection. Its real three-dimensional AMR
fixture allocates particle populations in owner cells but remains fixed at one
MPI rank. T10 supplies its separate collective ownership/migration matrix;
that evidence does not retroactively make the baseline physics fixture
distributed. Exact equations, seeds, sample counts, tolerances,
limitations, and case IDs are documented in
`test/srcCoreUnitTests/README.md`; the legacy GoogleTest mapping is in
`SRC_UNIT_TESTS_MIGRATION.md`.

The remaining T1 suites exercise the following bounded contracts:

- general arrays through rank five, managed stacks, algebraic helpers, and
  independent quadrature oracles;
- deterministic uniform/normal/gamma streams and isotropic directions, plus
  two-rank independent-stream evidence;
- two-rank scalar/vector/structure/zero-length `CMPI_channel` transfers;
- thread-safe diagnostic formatting, rotation, durable flush, and fatal
  subprocess behavior;
- SI/CGS/normalized unit factors and scalar/vector round trips;
- He/Na retained table-regression points, finite clamping, distance scaling,
  acceleration units, and typed invalid/unrepresentable domains;
- synthetic C and optional Fortran compiler/toolchain ABI echoes plus direct
  production Geopack layout, coordinate, unit, and no-SPICE behavior.

The ABI echoes do not qualify an AMPS production C/Fortran wrapper. Likewise,
the retained He/Na table literals exercise computational behavior but are not
a separately versioned digitization of the cited papers; that scientific
provenance remains visible in the component and source-traceability ledgers.

`CORE-INTERFACE-ABI-001` also owns explicit ARMS, BATSRUS, CCMC, Kameleon,
and SWMF provider cases.  On a machine without the required provider library,
configured coupling state, and versioned independent data fixture, each case
must report `SKIP`; a synthetic fixture must not be reported as provider
qualification.  Legacy/ISO-C Fortran cases similarly skip only when no
Fortran compiler is available.  Available ABI and Geopack cases still execute
and report their own states.

The native runner takes an advisory lock on the checkout's read-only
`Config.pl` inode and, on Linux, also scans `/proc` for legacy AMPS/
configuration processes that do not use that lock. This makes cooperating
jobs on cluster nodes share a lock whenever they share the checkout, even if
their `TMPDIR` values differ; independent checkouts use different inodes. On a
non-Linux machine the lock still serializes this runner, but the site scheduler
must keep unrelated legacy builds/runs from using the same checkout
concurrently.

Each run writes:

- `effective-config.json`: resolved machine settings and their source files;
- `selection.json`: immutable selected registry entries and stream plan;
- `streams/<stream>.log`: concise RUN/DONE history for one independent stream;
- `streams/<stream>.json`: that stream's selection, results, counts, and
  not-executed inventory;
- `logs/ID.log`: complete output of each individual runner;
- `work/ID/`: isolated transient state (moved below the configured machine
  scratch root when `temporary_root` is set);
- `artifacts/ID/`: runner-owned evidence;
- `amps-core-tests.json`: structured aggregate result referencing every stream
  receipt;
- `amps-core-tests.xml`: JUnit for CI/nightly systems.

After execution, the terminal prints `FINAL INDIVIDUAL TEST RESULTS` with one
row per reported case (or one explicit suite-level row when a suite could not
produce case evidence), followed by separate individual and suite summaries.
The JSON report exposes the same normalized view as `individual_results` and
`individual_counts`; `counts` remains the suite-level compatibility summary.
Thus `--all` ends with the complete case/state inventory even when earlier
case-start output was interleaved by multiple streams or suite workers.

Statuses are `PASS`, `FAIL`, `SKIP`, and `ERROR`. A missing prerequisite is a
typed `SKIP`; malformed/missing runner evidence, timeout, or harness failure is
`ERROR`; a tested scientific/software assertion is `FAIL`. A runner cannot
hide a failing case behind a passing suite status or a zero exit code.

During execution the terminal shows both suite and individual-case progress:

```text
[stream core] [2/3] RUN  CORE-GENERAL-INPUT-001 (src/general/input)
    [stream core] [CASE] RUN  INP-01 nested/named blocks, comments, continuations, and locations
```

Each adapter streams the `[CASE] RUN` record immediately before entering the
case body. The global runner relays those records live while retaining all
compiler, MPI, diagnostic, and result output in `logs/ID.log`. This keeps long
native builds readable without delaying the identity of a running case until
the entire suite finishes.

## Portable machine settings

Test implementations contain no site installation paths. The global CLI now
requires an explicit machine profile. Configuration precedence is:

1. committed `config/default.json`;
2. JSON files named by `AMPS_TEST_CONFIG` (colon-separated on Unix);
3. each repeated `--config FILE` from left to right;
4. scalar `AMPS_TEST_*` environment variables;
5. command-line overrides.

Copy the maintained template to a short-hostname profile and adjust it for a
workstation or cluster:

```sh
cp test/qualification-machine/template.json \
  test/qualification-machine/$(hostname -s).json
python3 src/test/run_tests.py \
  --config test/qualification-machine/$(hostname -s).json --tier nightly
```

The explicit `--config` option is required even if `AMPS_TEST_CONFIG` is set;
environment-named overlays remain an additional lower-precedence layer rather
than an implicit choice of execution machine. Repeated CLI profiles are still
supported and merge from left to right.

Relative paths resolve from the repository root. `${ENVIRONMENT_VARIABLE}` and
`~` are expanded. Commands are argument lists and are executed without a
shell, so wrapper commands and paths containing blanks remain unambiguous.
JSON arrays are recommended for tools that need wrapper arguments, for
example `"mpiexec": ["srun", "--mpi=pmix"]`. A string is accepted and parsed
with shell-like quoting only while loading the profile; it is still executed
later as an argument vector, not through a shell.

Keep site profiles outside the repository when they contain scheduler paths or
filesystem layout. A profile is an overlay: it needs to name only settings
that differ from `default.json`. Lists, including command vectors, replace the
default list as a whole; JSON objects merge recursively. Unknown top-level
keys, missing required tools, invalid positive integers, and invalid
`keep_work` values stop before a test runs.

Supported environment overrides are:

| Variable | Setting |
|---|---|
| `AMPS_TEST_CONFIG` | Additional JSON overlay path(s) |
| `AMPS_TEST_OUTPUT_ROOT` | Result root |
| `AMPS_TEST_TEMP_ROOT` | Machine scratch root reserved for runners |
| `AMPS_TEST_DATA_ROOT` | Optional root for large machine-local reference data |
| `AMPS_TEST_SPICE_PATH` | Optional CSPICE installation for registered application suites |
| `AMPS_TEST_SPICE_KERNELS` | Optional SPICE kernel tree for registered application suites |
| `AMPS_TEST_PYTHON` | Python command |
| `AMPS_TEST_CXX` | Dependency-light C++ compiler command |
| `AMPS_TEST_MPI_CXX` | MPI C++ compiler command |
| `AMPS_TEST_MAKE` | Make command |
| `AMPS_TEST_MPIEXEC` | MPI launcher command |
| `AMPS_TEST_MPI_RANK_FLAG` | Launcher rank option (`-n`, `-np`, etc.) |
| `AMPS_TEST_BUILD_JOBS` | Compile concurrency for non-native adapters |
| `AMPS_TEST_MPI_RANKS` | MPI ranks |
| `AMPS_TEST_TIMEOUT` | Default per-suite timeout in seconds |

The native AMPS rebuild remains fixed at `-j16` by repository policy even if a
general adapter uses a different configured build concurrency.

### Generated-configuration boundary

`CORE-PIC-NATIVE-001` owns the baseline `test/core-unit-tests` compile profile;
five additional exclusive owners select the individual-weight, aligned-record,
species-local, and single-global T2 variants. `CORE-MESH-NATIVE-001` separately
owns the `test/core-mesh-{1d,2d,3d}` profiles and verifies the effective
generated `DIM` before executing profile-owned cases. The baseline configuration cases are intended
to prove application selection through `Config.pl`, the final effective value
of each selected mutually exclusive macro in this profile, an invalid-deck
rejection, and retained SHA-256 identities. They do not qualify every
`Config.pl`/`ampsConfig.pl` option, species permutation, invalid application,
or every generated profile. Those exact partial and tested contracts
are enumerated in `configuration_contracts.json`.

Future generator tests must run in a validated isolated checkout/workspace or
use a non-mutating generator interface. They must never edit `.amps.conf`
directly. Required negative cases include unknown applications/options,
unknown species, missing values, and incompatible mode pairs; duplicate
top-level `SpeciesList` rejection is now a registered native case, but other
species permutations and generator modes remain open;
required positive cases compare generated macro values and species tables with
independent expectations, record a configuration identity, and verify clean-
build reproducibility. Individual-weight/split-merge and alternate
weight/time/layout profiles now have separate T2 contracts; alternate
collision, field-solver, and distributed profiles remain open because they
cannot be inferred from the qualified neutral one-rank matrix.

The mesh input decks match the repository's existing `*.input` ignore rule.
Codex does not modify Git control files; a human preparing a commit must use an
explicit force-add for `input/test/core-mesh-{1d,2d,3d}.input` so the three
configured T3 profiles are not omitted.

### T10 collective MPI matrix

`CORE-MESH-MPI-OWNERSHIP-001` is the exclusive owner of 16 collective cases.
Its aliases `CORE-PIC-MPI-HALO-001`, `CORE-PIC-MPI-MIGRATION-001`, and
`CORE-PIC-MPI-RANK-INVARIANCE-001` select bounded subsets without duplicate
builds. The owner retains nonperiodic and periodic executables, then launches
every selected case independently at 1, 2, and 4 ranks. A matrix element is
accepted only when every rank publishes its own PASS record and rank zero
publishes exactly one collective summary.

```sh
python3 src/test/run_tests.py --config test/qualification-machine/passion.json \
  --test CORE-MESH-MPI-OWNERSHIP-001 \
  --output-dir /tmp/amps-t10
```

The two input decks are ignored by the repository's broad `*.input` rule; a
human preparing a commit must explicitly force-add
`input/test/core-mpi.input` and `input/test/core-mpi-periodic.input`. The
fixture contract, independent equations, exact case inventory, qualified
boundary, and evidence files are documented in
`test/srcCoreMPITests/README.md`.

### T11 local-registry adapters

Six extended suites import 637 exact identities without duplicating their
model-local physics oracles. Large identity lists are stored in
`case_registries.json`; `lib/manifest.py` resolves them through the same
validation path used by selection, reporting and traceability. Every adapter
validates the complete ordered central/local registry before applying an exact
case selection, so deleting an expensive or failing case cannot create a
smaller green campaign.

The current T11 evidence is 402 PASS, 233 SKIP, one FAIL and one ERROR. The
full-volume `CMBGU05` leakage failure and `CSWC0628-MASS-TIME` refinement
timeout remain visible. SEP-common's 52 long Monte-Carlo acceptance cases,
SWCME's unavailable extended/observational registries, SCCM's three missing
campaign-evidence cases, and reduced-front's 35 native cases remain typed
SKIPs. A reduced-front PASS never qualifies BG3D-4/G1, and the legacy audit
proves only explicit directory/owner dispositions rather than complete legacy
behavioral migration. Exact commands, evidence boundaries and limitations are
documented beside each adapter under `tests/models/` and `tests/legacy/`.

### T12 bounded portability matrix

Run the registered matrix with:

```sh
python3 src/test/run_tests.py --config test/qualification-machine/passion.json \
  --test CORE-T12-PORTABILITY-001 \
  --output-dir /tmp/amps-t12-portability --jobs 1
```

The matrix uses a private general-header mirror and writes all binaries below
the suite work directory; it never changes `.amps.conf`, invokes `Config.pl`,
or rebuilds native AMPS. It covers GNU debug/optimized, Fortran ABI, ASan,
UBSan, selected four-thread TSan, Valgrind, OpenMP invariance, explicit
OpenMPI/MPICH pairs, CUDA compile/runtime, two broad performance guardrails,
and strict coverage as 16 separate cases. Optional tool absence is SKIP, a
present tool that violates its contract is FAIL, and the strict-coverage case
reports the real inventory exit.

Machine files may override optional argv vectors using `fortran`,
`clang_cxx`, `valgrind`, `nvcc`, `openmpi_cxx`, `openmpi_exec`, `mpich_cxx`
and `mpich_exec`; see `test/qualification-machine/template.json`. This is the
mechanism for site wrappers and paths—test code contains no site-specific
absolute path.
On the current host the matrix is intentionally red at 10 PASS, five SKIP and
one FAIL because strict coverage still has 2 covered and 22 incomplete
components. See `tests/portability/matrix/README.md` for the exact boundary;
whole-AMPS sanitizers and production mesh/push/collision/halo/restart
performance budgets remain open.

## Add tests to the global procedure

This section is the contributor procedure for registering a test, writing its
adapter, producing independent evidence, and verifying that the global runner
cannot silently omit it. Architectural rationale and coverage policy remain in
[DESIGN.md](DESIGN.md); the staged coverage roadmap remains in
[`TEST.md`](TEST.md).

> **Mandatory rule for agents and contributors:** Any change that adds,
> removes, renames, reorders, or changes ownership/description of a suite or
> individual case **must** regenerate [`TEST_CATALOG.md`](TEST_CATALOG.md) with
> `python3 src/test/generate_test_catalog.py` in the same change. The work is
> incomplete until `python3 src/test/generate_test_catalog.py --check` passes.
> Do not hand-edit the generated catalog and do not defer its update to a later
> agent.

### 1. Decide whether to add a case or a suite

Prefer a new **case in an existing suite** when it uses the same production
component, compiled configuration, prerequisites, timeout class, isolation
boundary, and runner protocol. Add a **new suite** when any of those differ—for
example, a new generated AMPS configuration, MPI rank layout, accelerator,
external dataset, long acceptance campaign, or repository-wide build owner.

Do not create a new suite merely to avoid updating an existing exact case
registry. Conversely, do not hide mutually incompatible compile-time modes in
one runner. One suite is one schedulable isolation and evidence owner.

Use these levels consistently:

| Level | Intended boundary |
| --- | --- |
| `unit` | One dependency-light function, class, parser, or exact algebraic contract |
| `module` | Several production translation units or one configured subsystem |
| `integration` | Cross-module, native application, MPI, restart, or external-provider boundary |
| `system` | Complete application workflow with an explicitly bounded scientific claim |

Campaigns express cost and environment, not importance:

- `smoke`: fast, deterministic bring-up checks;
- `nightly`: deterministic routine CPU coverage;
- `profile`: mutually exclusive generated configurations;
- `extended`: long convergence, stochastic, compiler, or sanitizer work;
- `distributed`: multi-rank ownership/communication/restart work;
- `application-smoke`: application plumbing without scientific qualification.

A suite may belong to multiple campaigns. Do not add a test to `smoke` unless
its prerequisites and runtime are appropriate for every supported smoke host.

### 2. Choose stable identities

Suite, alias, registry, and case IDs must match `[A-Z][A-Z0-9-]*`. Case IDs are
globally unique across the complete manifest. Use a component/family prefix,
the behavior being asserted, and a numeric suffix, for example
`GENERAL-ARRAY-RESIZE-001`. Never encode a temporary hostname, output path,
random seed, date, or current implementation detail in an ID.

IDs are evidence keys and must not be renamed casually. If behavior is
replaced, preserve the prior result/evidence boundary or document a migration;
do not reuse its ID for a different assertion.

For a short registry, declare cases inline in `manifest.json`. Large imported
or staged registries belong in `case_registries.json` and are referenced by
`expected_cases_ref`. The loader validates every stored registry—even one not
yet referenced—and rejects empty, duplicate, malformed, or globally duplicate
case IDs. A suite must specify exactly one of `expected_cases` and
`expected_cases_ref`.

### 3. Register the suite in `manifest.json`

A minimal dependency-light entry is:

```json
{
  "id": "CORE-EXAMPLE-001",
  "component": "src/models/example",
  "level": "unit",
  "tiers": ["nightly"],
  "runner": "tests/models/example/run_test.py",
  "description": "Independent example-model identities and invalid-domain behavior",
  "requires": {
    "tools": ["cxx"],
    "paths": ["src/models/example/example.h"],
    "capabilities": []
  },
  "expected_cases": [
    "EXAMPLE-IDENTITY-001",
    "EXAMPLE-INVALID-DOMAIN-001"
  ],
  "timeout_seconds": 120
}
```

Manifest fields have these contracts:

| Field | Requirement |
| --- | --- |
| `id` | Globally unique stable suite ID; aliases also share this namespace |
| `component` | Existing production path below `src/`, or an explicitly approved application root such as `srcMoon`/`srcEnceladus`; never a test directory |
| `level` | One of `unit`, `module`, `integration`, or `system` |
| `tiers` | Nonempty unique subset of the declared campaigns |
| `runner` | Existing file below `src/test`, written relative to `src/test` |
| `description` | Specific tested behavior and boundary, not “run unit tests” |
| `requires.tools` | Configured tool-vector names, not executable paths |
| `requires.paths` | Repository-relative prerequisite files/directories |
| `requires.capabilities` | Named optional/required machine capabilities |
| `expected_cases` | Ordered, unique, nonempty exact case list |
| `expected_cases_ref` | Out-of-line registry ID used instead of `expected_cases` |
| `timeout_seconds` | Positive suite deadline including build and execution |
| `mpi_ranks` | Positive rank request when the suite has one declared layout |
| `exclusive` | `true` when checkout-wide generated state can be changed |

`aliases` and `alias_cases` are compatibility/selection mechanisms, not a way
to give one result several identities. Each alias subset must be a nonempty
subset of the owning suite's expected cases. Add aliases only when an existing
stable selection name must continue to work.

Set `exclusive: true` if the runner or any transitive command invokes
`Config.pl`/`ampsConfig.pl`, removes or rebuilds root `build`, changes generated
headers or `.general.conf`, writes source-adjacent objects, or otherwise changes
checkout-wide state. Unknown isolation is exclusive. A private source copy
does not become nonexclusive until all transitive writes are proven to remain
below its unique work/artifact roots. See
[`../../TEST_IMPROVEMENT_PLAN.md`](../../TEST_IMPROVEMENT_PLAN.md) for the
isolation audit and future resource model.

### 4. Create the runner

Put the adapter at
`src/test/tests/<component>/<suite>/run_test.py`. A fully commented,
fail-closed starting point is
[`templates/run_test.py.template`](templates/run_test.py.template). After
copying it:

1. replace every `TEMPLATE` marker;
2. keep its complete ordered registry equal to the manifest owner registry;
3. implement independent case assertions;
4. retain useful transcripts/references below `artifact_dir`; and
5. run the validation sequence in Section 11.

The global runner invokes an adapter only as:

```text
python3 path/to/run_test.py --context resolved-context.json
```

Use `load_context()` rather than parsing additional suite-specific CLI flags.
The schema-1 context provides:

- `repo_root`, `work_dir`, `artifact_dir`, and the controlled `result_file`;
- the resolved suite manifest entry, including `expected_cases` and
  `selected_cases`;
- portable tool argv vectors and resolved data/output paths; and
- execution settings such as build jobs, MPI ranks, and timeout policy.

Use `command(context, "cxx")` (or another configured tool name) to obtain an
argv vector. Use `run_command()` to execute it without a shell, stream output,
enforce a timeout, and clean up the owned process group. Do not concatenate a
shell string, interpolate an untrusted path, or discover a site compiler from
a hard-coded absolute path.

Before executing a selected subset, compare the adapter's complete ordered
registry with `context["test"]["expected_cases"]`. Then honor
`selected_cases` in its supplied order. This proves that `--case`, reverse,
and shuffle selection cannot conceal a deleted manifest or local case.

Emit and flush this record immediately before each assertion:

```text
[CASE] RUN  <stable-case-id> <short description>
```

The line is live progress, not evidence. The atomic schema-1 result is
authoritative.

### 5. Publish the exact result protocol

Use `case()` and `write_result()` from `lib/runner_support.py`. A runner result
contains:

```json
{
  "schema_version": 1,
  "id": "CORE-EXAMPLE-001",
  "status": "PASS",
  "message": "two example cases passed",
  "cases": [
    {"id": "EXAMPLE-IDENTITY-001", "status": "PASS", "message": "residual <= tolerance"},
    {"id": "EXAMPLE-INVALID-DOMAIN-001", "status": "PASS", "message": "typed rejection observed"}
  ],
  "evidence": ["/absolute/artifact/path/example-output.txt"],
  "metrics": {"maximum_residual": 1.2e-13}
}
```

The global runner requires exact selected case content **and order**, unique
IDs, a suite status at least as severe as every child, and the corresponding
process exit code:

| Status | Meaning | Exit |
| --- | --- | ---: |
| `PASS` | The requested assertion ran and met its independent acceptance rule | 0 |
| `SKIP` | A genuinely optional capability/reference is unavailable or the declared profile does not execute it | 0 |
| `FAIL` | The tested behavior ran and violated its assertion, residual, conservation law, or required build contract | 1 |
| `ERROR` | No trustworthy tested conclusion exists because setup, parsing, timeout, filesystem, protocol, or control infrastructure failed | 2 |

Severity is `ERROR > FAIL > SKIP > PASS`. A mixed PASS/SKIP suite is SKIP; a
suite cannot claim PASS while hiding a failed child. Missing required
capability policy is handled by the global prerequisite check and may be ERROR
rather than SKIP.

Catch common setup/control exceptions inside `main()` and return
`error_result(context, error)`. It attributes the infrastructure ERROR to all
selected case IDs, preventing a secondary empty-registry diagnostic from
hiding the original failure. `main_guard()` is only the last boundary for a
failure that occurs before a trustworthy context/result can be written.

### 6. Build independent scientific and numerical oracles

Do not calculate an expected value with the production routine under test.
Use one or more of:

- an exact identity or conservation law;
- an independently evaluated analytical solution;
- a manufactured solution with known convergence order;
- an independently implemented high-accuracy quadrature/reference algorithm;
- a versioned external fixture with bibliographic/file provenance; or
- cross-decomposition/rank/thread invariants where exact ordering permits.

Every numerical reference must document:

- physical quantity, SI or declared native units;
- coordinate frame, orientation, and sign convention;
- equations and assumptions;
- fixture source/version and any transformation performed;
- absolute/relative tolerance and its truncation, conditioning, or uncertainty
  rationale; and
- valid domain and typed invalid-domain behavior.

Test limits, zero states, invalid values, unsupported modes, and conservation
budgets—not only one nominal value. Convergence claims require multiple
resolutions and an independently specified expected order. Stochastic tests
require deterministic seeds for reproducibility plus acceptance bounds derived
from the distribution/sample count; one recorded production output is not a
statistical oracle.

Keep small references beside the runner. Large or machine-provided references
must be addressed through the configured `data_root` and a committed manifest
containing provenance and integrity information.

### 7. Add C++ production-linked cases

Dependency-light C++ tests may use `native/test_harness.h`; it has no
GoogleTest, PIC, MPI, or application dependency. Compile objects and
executables below `work_dir`, not next to sources. Parse the C++ harness's exact
`[ID] PASS|FAIL message` rows and reject missing, duplicate, foreign, or
out-of-order records even when the process exits zero.

For a configured AMPS application, use `test/srcCoreUnitTests` and keep the
suite exclusive. AMPS links application objects through a static archive, so a
translation unit containing only static `AMPS_TEST_CASE` registrations may be
discarded by the linker. Export a small link-anchor function from that
translation unit, call it from `main.cpp`, and add every registered case ID to
the owning adapter. The exact-ID registry must fail if a source file or link
anchor is omitted even when all remaining cases pass.

Before every native rebuild, follow the repository preflight exactly:

1. confirm `/home/vtenishe/Mars2/AMPS` and that no build/test process is active;
2. remove only root `build` with `rm -rf -- build` from that confirmed root;
3. regenerate the selected application and production hooks through the
   established configuration workflow; and
4. compile through the top-level workflow with `-j16`.

Never edit `.amps.conf` directly. A direct source-only compile is useful for a
dependency-light fixture but is not native AMPS build evidence.

### 8. Add MPI, threaded, or accelerator tests

Declare MPI tools/capabilities and ranks in the manifest/machine profile. Use
the configured launcher argv and rank flag; never assume `mpiexec`, `-n`, a
hostfile, scheduler command, or network interface. Every rank must publish its
local assertion or contribute to a checked exact reduction. Rank-zero-only
success is insufficient for ownership, migration, halo, or restart claims.

Use distinct identities for serial, MPI implementations, thread layouts,
compile-only accelerator evidence, and runtime accelerator evidence. Missing
optional hardware can be SKIP; a compiler PASS cannot qualify a device runtime.
Use unique output/restart paths and account for every allocated rank/thread or
device. Tests that share launch fabrics, devices, generated configurations, or
fixed services remain exclusive until the resource scheduler described in
`TEST_IMPROVEMENT_PLAN.md` is implemented and qualified.

### 9. Keep paths and machine configuration portable

All writes belong below `work_dir` or `artifact_dir`:

- `work_dir`: disposable objects, executables, scratch files, and temporary
  copied source trees;
- `artifact_dir`: retained transcripts, residual tables, reference comparison,
  plots, and suite-owned evidence;
- `result_file`: controlled by `write_result()`; never choose or overwrite it;
- global logs: written by the scheduler below `<output>/logs/<suite-id>.log`.

Do not write into production sources, root `build` from a nonexclusive suite,
a shared application output directory, `$HOME`, or a fixed `/tmp/name`. Do not
embed `/home/...`, compiler installations, SPICE/data roots, MPI launchers, or
site scratch paths in test code. Add a documented argv/path/capability field to
a machine overlay when a site needs one; see
`test/qualification-machine/template.json`.

`run_command()` receives argv elements, not shell syntax. Do not use pipes,
redirections, command substitution, or environment-variable expansion as
arguments and expect a shell to interpret them. Perform small file operations
with Python APIs and retain the resulting paths as evidence.

### 10. Update coverage and traceability ledgers

Every new suite must be named in the owning entry of
`component_inventory.json`. Keep the entry `partial` unless the documented
remaining behaviors are genuinely closed; one green case never automatically
makes a component `covered`. An `exempt` state needs an explicit reviewed
rationale.

Also update the relevant ledger when the change alters its ownership:

- `source_traceability.json` for a new production source/public-contract
  ownership rule or changed suite mapping;
- `configuration_contracts.json` for a new generated/configured mode;
- `t1_public_algorithm_inventory.json` for a newly tested or reclassified
  public algorithm family;
- `case_registries.json` for a large out-of-line exact registry; and
- component/model README and `src/test/STATUS.md` for commands, evidence,
  limitations, and qualification state.

Run `validate_traceability.py` after any registry or ownership change. Its PASS
proves structural classification integrity; it does not substitute for a
behavioral test.

### 11. Validate the contribution

From the repository root, replace the example IDs and output paths below:

```sh
# Confirm manifest loading, ownership, and exact registered case count.
PROFILE=test/qualification-machine/passion.json
python3 src/test/run_tests.py --config "$PROFILE" --test CORE-EXAMPLE-001 --list

# Run one case, then the complete owning suite.
python3 src/test/run_tests.py --config "$PROFILE" --case EXAMPLE-IDENTITY-001 \
  --output-dir /tmp/amps-example-one --jobs 1
python3 src/test/run_tests.py --config "$PROFILE" --test CORE-EXAMPLE-001 \
  --output-dir /tmp/amps-example-suite --jobs 1

# Prove ordering does not change case meaning or final identity accounting.
python3 src/test/run_tests.py --config "$PROFILE" --test CORE-EXAMPLE-001 --order reverse \
  --output-dir /tmp/amps-example-reverse --jobs 1
python3 src/test/run_tests.py --config "$PROFILE" --test CORE-EXAMPLE-001 --order shuffle \
  --shuffle-seed 1729 --output-dir /tmp/amps-example-shuffle --jobs 1

# Validate the control plane and structural ledgers.
python3 src/test/run_tests.py --config "$PROFILE" --test CORE-HARNESS-001 \
  --output-dir /tmp/amps-example-harness --jobs 1
python3 src/test/generate_test_catalog.py --check
python3 src/test/validate_traceability.py \
  --output /tmp/amps-example-traceability.json

# Run every campaign to which the suite was added, then the full selection
# when its cost and prerequisites make that appropriate.
python3 src/test/run_tests.py --config "$PROFILE" --campaign nightly \
  --output-dir /tmp/amps-example-nightly --jobs 8
```

Inspect all of:

- terminal individual results and counts;
- `<output>/amps-core-tests.json`;
- `<output>/amps-core-tests.xml`;
- `<output>/logs/<suite-id>.log`;
- the suite artifact directory and every path named by `evidence`; and
- `not_executed_suites` if fail-fast was used.

Confirm that JSON/JUnit contain every selected identity exactly once and that
FAIL/ERROR/SKIP messages retain their real reasons. Run `git diff --check` as a
read-only formatting check. The repository-wide strict-coverage gate remains
red while inventory components are partial; do not misrepresent that known
global state as a failure of an otherwise valid new test, and do not weaken it
to obtain a green contribution.

Native suites additionally require their clean-build preflight evidence and
the applicable one-/multi-rank execution receipts. Optional-provider tests
require both the available behavior and the correctly typed unavailable path.

### 12. Contributor checklist

- [ ] The case belongs in this suite, or the new suite has a distinct justified isolation/configuration owner.
- [ ] Suite and case IDs are stable, valid, globally unique, and exactly ordered.
- [ ] Manifest component, level, campaigns, requirements, timeout, ranks, and exclusivity are accurate.
- [ ] The runner accepts only `--context`, validates the complete registry, and honors exact selection order.
- [ ] Every case emits `[CASE] RUN` and publishes one structured result.
- [ ] Suite severity and process exit agree with all child statuses.
- [ ] Common infrastructure failures use `error_result()` and retain selected identities.
- [ ] Expected values are independent and units, frame, provenance, domains, and tolerances are documented.
- [ ] Invalid/limit states and conservation or convergence behavior are tested where applicable.
- [ ] All generated files stay below declared work/artifact roots.
- [ ] Tool and data locations come from the context/machine profile; no site path is hard-coded.
- [ ] Configured native builds are exclusive and follow the mandatory clean preflight.
- [ ] MPI/thread/device evidence accounts for every participant and does not broaden another backend's claim.
- [ ] Component inventory and every affected traceability/configuration ledger are updated honestly.
- [ ] `TEST_CATALOG.md` is regenerated in this change and its `--check` gate passes.
- [ ] Exact-case, complete-suite, reverse, shuffle, harness, traceability, and applicable campaign runs are checked.
- [ ] JSON, JUnit, logs, artifacts, counts, and typed non-PASS outcomes reconcile.

For complete worked examples, compare a dependency-light adapter such as
`tests/species/radiation/run_test.py`, a native owner such as
`tests/pic/native_core/run_test.py`, an MPI owner such as
`tests/pic/mpi_distributed/run_test.py`, and an imported-registry adapter under
`tests/models/cme/`.
