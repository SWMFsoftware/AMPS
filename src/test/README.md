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

T2 extends the live registry, without rewriting that T1 receipt, to 11 suites
and 103 cases.  The native suite now owns 36 cases: five new cases qualify
particle statistical-weight and time-step round trips, species-distribution
copies, and typed invalid-value/index rejection in the host species-global
profile.  The retained normal T2 run in
`/tmp/amps-t2-all-qualified-normal` reports 98 PASS, five provider SKIPs, and
no FAIL/ERROR; focused evidence is in `/tmp/amps-t2-weight-time-initial`, with
the collision regression recheck in
`/tmp/amps-t2-weight-time-collision-recheck`.  All 36 native cases also pass
in reverse order and deterministic seed-1729 shuffled order under
`/tmp/amps-t2-native-qualified-reverse` and
`/tmp/amps-t2-native-qualified-shuffle-1729`; the corrected normal, reverse,
and shuffled builds have identical executable SHA-256
`6329674e470d9cadda436aca2d410f1ef5ecad27f4918108d56bd3078a1f4925`.
Species-local and single-global profiles, individual weight correction and
split/merge, unequal-weight/time collision physics, and CUDA device
diagnostics remain open.

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
python3 src/test/run_tests.py --campaign nightly --jobs 4
python3 src/test/run_tests.py --all --list
python3 src/test/run_tests.py --coverage
make -C src/test smoke
make core_test_nightly
```

`--tier` remains a backward-compatible alias for `--campaign`. The top-level
`test_all` target depends on `core_test_nightly`. In the current manifest,
`smoke` selects eight dependency-light suites (51 cases), `nightly` selects
all 11 suites (103 cases), `profile` selects the 36-case configured native
suite, `extended` selects the 11-case interface suite plus the native suite,
and `distributed` selects the five-case two-rank MPI channel/RNG suite.
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

The configured AMPS application test is exclusive and is selected by the
`nightly`, `profile`, and `extended` campaigns:

```sh
python3 src/test/run_tests.py --test CORE-PIC-NATIVE-001
```

The historical particle-buffer and physics suite IDs are manifest aliases for
selection compatibility; they are not separate executable suites and must not
trigger two rebuilds.

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

The combined native suite owns three generated-configuration adapter cases,
three particle-buffer cases, 19 physics cases, four compiled/generated
species cases, two production-linked global-unit cases, and five particle
weight/time-step cases.  The T2 cases exercise the configured host
species-global minimum-reduction and distribution-copy behavior and require
typed rejection of nonpositive, nonfinite, and out-of-range inputs.  They do
not yet qualify species-local or single-global profiles, device-side error
behavior, individual correction, split/merge, or unequal-weight/time collision
physics.  The physics cases add analytical/invariant checks for Huebner
rates, branching, products, and exponential parent decay; Maxwellian moments;
hard-sphere scattering; populated-cell NTC collision rates; Larsen-Borgnakke
rotational/vibrational laws and RT equilibration; neutral Simple/Boris motion;
and analytic spherical-surface reflection. Its real three-dimensional AMR
fixture allocates particle populations in owner cells but remains fixed at one
MPI rank. Distributed ownership/migration evidence belongs in later
collective suites. Exact equations, seeds, sample counts, tolerances,
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

`CORE-PIC-NATIVE-001` is the single registered owner of the current
`test/core-unit-tests` compile profile. Its configuration cases are intended
to prove application selection through `Config.pl`, the final effective value
of each selected mutually exclusive macro in this profile, an invalid-deck
rejection, and retained SHA-256 identities. They do not qualify every
`Config.pl`/`ampsConfig.pl` option, species permutation, invalid application,
or alternative generated profile. Those exact partial and planned contracts
are enumerated in `configuration_contracts.json`.

Future generator tests must run in a validated isolated checkout/workspace or
use a non-mutating generator interface. They must never edit `.amps.conf`
directly. Required negative cases include unknown applications/options,
unknown species, missing values, and incompatible mode pairs; duplicate
top-level `SpeciesList` rejection is now a registered native case, but other
species permutations and generator modes remain open;
required positive cases compare generated macro values and species tables with
independent expectations, record a configuration identity, and verify clean-
build reproducibility. Individual-weight/split-merge, alternate collision,
field-solver, and distributed profiles remain separate contracts because they
cannot be inferred from the active neutral one-rank profile.

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
