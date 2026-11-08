# AMPS core unit and module tests

`src/test` is the global test-control layer for production code under `src`.
It does not replace application validation or model-specific qualification.
It gives those tests one portable registry, one result vocabulary, and one
nightly-compatible evidence format.

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
python3 src/test/run_tests.py --campaign smoke
python3 src/test/run_tests.py --campaign nightly
python3 src/test/run_tests.py --campaign profile
python3 src/test/run_tests.py --campaign nightly --jobs 4
python3 src/test/run_tests.py --all --list
python3 src/test/run_tests.py --coverage
make -C src/test smoke
make core_test_nightly
```

`--tier` remains a backward-compatible alias for `--campaign`. The top-level
`test_all` target depends on `core_test_nightly`. The live T10 registry has
42 suites and 375 cases; use `--list` for exact campaign memberships because
later roadmap stages extend those sets. The `distributed` campaign selects
the five-case two-rank channel/RNG suite and the 16-case T10 collective
1/2/4-rank matrix (two suites and 21 cases total).
`--all` ignores campaign membership and selects the complete registry.
Site/nightly orchestration must preserve the native suite's exclusive barrier
because it selects a repository-wide application.

Run one suite or component:

```sh
python3 src/test/run_tests.py --test CORE-GENERAL-INPUT-001
python3 src/test/run_tests.py --test CORE-GENERAL-MPI-001
python3 src/test/run_tests.py --test CORE-INTERFACE-ABI-001
python3 src/test/run_tests.py --component src/general
python3 src/test/run_tests.py --match 'CORE-GENERAL-*'
```

Run exact cases or deterministic isolation orders without broadening the
selection:

```sh
python3 src/test/run_tests.py --case INP-03
python3 src/test/run_tests.py --case PHYS-IDF-RT-001 --order reverse
python3 src/test/run_tests.py --all --order shuffle --shuffle-seed 1729
```

`--case` may be repeated. A requested ID must have exactly one manifest owner;
the runner forwards only those selected IDs to the owning adapter and validates
their ordered result list. Reverse/shuffle support is an execution mechanism,
not proof of isolation by itself; the retained normal/reverse/shuffle and
one-case evidence is recorded in `STATUS.md`.

Each configured AMPS application test is exclusive and is selected by its
declared campaigns. The baseline profile can be run directly with:

```sh
python3 src/test/run_tests.py --test CORE-PIC-NATIVE-001
```

Run the complete T3 dimensional matrix and its isolation orders with:

```sh
python3 src/test/run_tests.py --test CORE-MESH-NATIVE-001 --order normal \
  --output-dir /tmp/amps-t3-qualified-normal
python3 src/test/run_tests.py --test CORE-MESH-NATIVE-001 --order reverse \
  --output-dir /tmp/amps-t3-qualified-reverse
python3 src/test/run_tests.py --test CORE-MESH-NATIVE-001 --order shuffle \
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
python3 src/test/run_tests.py --test CORE-PIC-MOVER-CHARGED-001 \
  --test CORE-PIC-MOVER-RELATIVISTIC-001 --test CORE-PIC-GYRO-001 \
  --test CORE-PIC-BOUNDARY-001 --test CORE-PIC-INJECTION-001 \
  --test CORE-PIC-SAMPLING-001 --test CORE-PIC-COLLISION-ALGORITHMS-001 \
  --test CORE-PIC-IDF-001

python3 src/test/run_tests.py --test CORE-PIC-PERIODIC-001 \
  --test CORE-PIC-COLLISION-VHS-001 --test CORE-PIC-COLLISION-VSS-001 \
  --test CORE-PIC-IDF-ROTZ-Z2-001 --test CORE-PIC-IDF-ROTZ-Z10-001 \
  --test CORE-PIC-COLLISION-INDIVIDUAL-WEIGHT-001
```

Run the bounded T6 reaction and chemistry gate with:

```sh
python3 src/test/run_tests.py --test CORE-CHEMICAL-SOLVER-001 \
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
python3 src/test/run_tests.py --test CORE-SURFACE-001 \
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
python3 src/test/run_tests.py --test CORE-ECSIM-OPERATORS-001
python3 src/test/run_tests.py --test CORE-ECSIM-BOUNDARIES-001 \
  --test CORE-PIC-DEPOSITION-001 --test CORE-PIC-FIELD-SOLVER-001 \
  --test CORE-PIC-AUXILIARY-001 --test CORE-PIC-FIELD-LINE-001 \
  --test CORE-PIC-CALLBACK-OUTPUT-001 --jobs 4
python3 src/test/run_tests.py --test CORE-PIC-RESTART-001 \
  --test CORE-PIC-DATAFILE-001 --test CORE-PIC-OUTPUT-001 --jobs 3
python3 src/test/run_tests.py --test CORE-PIC-RESTART-NATIVE-001
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

`--jobs N` controls suite concurrency and defaults to `1`. Only consecutive
nonexclusive suites may overlap. An `exclusive: true` suite forms a barrier:
all earlier work completes before it starts, it runs alone, and later work is
not submitted until it finishes. This is required because a native suite
selects a checkout-wide application and rebuilds `build`. With `--jobs > 1`,
streamed case-start lines include the suite ID because independent output can
interleave. `--fail-fast` finishes an already-submitted parallel batch and
then prevents later batches from starting.

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
- `selection.json`: immutable selected registry entries;
- `logs/ID.log`: complete output of each individual runner;
- `work/ID/`: isolated transient state (moved below the configured machine
  scratch root when `temporary_root` is set);
- `artifacts/ID/`: runner-owned evidence;
- `amps-core-tests.json`: structured aggregate result;
- `amps-core-tests.xml`: JUnit for CI/nightly systems.

After execution, the terminal prints `FINAL INDIVIDUAL TEST RESULTS` with one
row per reported case (or one explicit suite-level row when a suite could not
produce case evidence), followed by separate individual and suite summaries.
The JSON report exposes the same normalized view as `individual_results` and
`individual_counts`; `counts` remains the suite-level compatibility summary.
Thus `--all` ends with the complete case/state inventory even when earlier
case-start output was interleaved by `--jobs`.

Statuses are `PASS`, `FAIL`, `SKIP`, and `ERROR`. A missing prerequisite is a
typed `SKIP`; malformed/missing runner evidence, timeout, or harness failure is
`ERROR`; a tested scientific/software assertion is `FAIL`. A runner cannot
hide a failing case behind a passing suite status or a zero exit code.

During execution the terminal shows both suite and individual-case progress:

```text
[2/3] RUN  CORE-GENERAL-INPUT-001 (src/general/input)
    [CASE] RUN  INP-01 nested/named blocks, comments, continuations, and locations
```

Each adapter streams the `[CASE] RUN` record immediately before entering the
case body. The global runner relays those records live while retaining all
compiler, MPI, diagnostic, and result output in `logs/ID.log`. This keeps long
native builds readable without delaying the identity of a running case until
the entire suite finishes.

## Portable machine settings

No test contains an absolute installation path. Configuration precedence is:

1. committed `config/default.json`;
2. JSON files named by `AMPS_TEST_CONFIG` (colon-separated on Unix);
3. each repeated `--config FILE` from left to right;
4. scalar `AMPS_TEST_*` environment variables;
5. command-line overrides.

Copy `config/machine.example.json` outside the repository (or to the ignored
`config/local.json`) and adjust it for a workstation or cluster:

```sh
python3 src/test/run_tests.py \
  --config src/test/config/local.json --tier nightly
```

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
python3 src/test/run_tests.py --test CORE-MESH-MPI-OWNERSHIP-001 \
  --output-dir /tmp/amps-t10
```

The two input decks are ignored by the repository's broad `*.input` rule; a
human preparing a commit must explicitly force-add
`input/test/core-mpi.input` and `input/test/core-mpi-periodic.input`. The
fixture contract, independent equations, exact case inventory, qualified
boundary, and evidence files are documented in
`test/srcCoreMPITests/README.md`.

## Add an individual test

1. Give the suite a stable ID and add it to `manifest.json`. Declare the exact
   production component, level, tiers, requirements, timeout, and whether it
   mutates repository-wide configuration (`exclusive`).
2. Put its runner below `tests/<component>/<suite>/run_test.py`. It receives
   only `--context FILE`; use `lib/runner_support.py` to load tools/paths, run
   commands, and atomically write the versioned result.
3. Put hand-derived fixtures and references beside that runner. Never generate
   an expected solution with the production routine under test.
4. Add the suite ID to `component_inventory.json`, but mark a component
   `covered` only when its required unit, module, MPI, restart, invalid-input,
   and numerical-convergence behavior is genuinely exercised.
5. Run the harness self-test, the new suite alone, smoke, and the appropriate
   nightly/extended tier. Inspect both JSON and JUnit receipts.

Individual C++ tests may use `native/test_harness.h`; it has no GoogleTest,
PIC, MPI, or application dependency. Tests that need configured PIC use the
application in `test/srcCoreUnitTests` and remain serialized because AMPS
application selection is repository-wide.

For a configured application, remember that AMPS links application objects
through a static archive. A translation unit containing only static
`AMPS_TEST_CASE` registrations may be discarded by the linker. Export a small
link-anchor function from that translation unit and call it from `main.cpp`,
then add every expected case ID to the suite adapter. The exact-ID registry
will fail if a source file is omitted even when all remaining cases pass.

See [DESIGN.md](DESIGN.md) for test levels, component plans, acceptance rules,
and the staged route to complete nightly coverage.
