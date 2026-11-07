# AMPS core-test status

Last updated: 2026-11-07

## Implemented

- Portable global runner with validated registry, isolated work/artifact
  directories, timeouts, typed PASS/FAIL/SKIP/ERROR states, JSON, and JUnit.
- The global runner now ends every execution with a stable
  `FINAL INDIVIDUAL TEST RESULTS` table.  It prints one PASS/FAIL/SKIP/ERROR
  row for each executed case, then separate individual-case and suite
  summaries.  The JSON receipt preserves the same normalized rows in
  `individual_results` and their totals in `individual_counts`; suites which
  cannot publish case evidence are explicitly labeled as suite-level rows.
- Layered machine configuration with external JSON profiles, environment
  overrides, relative-path resolution, and no hard-coded installation paths.
- Explicit component inventory with structural drift audit and a separate
  strict-completeness gate. `source_traceability.json` now classifies every
  discovered production translation unit and public compile contract, while
  `configuration_contracts.json` records the machine, Config.pl, generated-
  file, compile-profile, and runtime-mode boundaries.
- The dependency-free `validate_traceability.py` expands those rules, rejects
  unclassified/overlapping files and stale suite/case/local-registry owners,
  and cross-checks the 24-area component inventory. This is structural
  governance evidence, not behavioral physics coverage.
- Eleven registered harness self-tests cover the original configuration,
  registry, selection, protocol, scheduler, and report contracts plus schema,
  hostile-subprocess failures, real timeout and SIGINT process-tree cleanup,
  bounded concurrency/fail-fast, and the full source/configuration
  traceability contract.
- Generic runtime-input parser adapter (four existing independent fixtures).
- General utility tests for arrays, bit arrays, stack allocation/reuse, linked
  containers, string helpers, and sparse 4-D indexing.
- New dependency-free C++ harness shared by standalone and configured tests.
- One combined `CORE-PIC-NATIVE-001` registered owner for the configured
  `test/srcCoreUnitTests` application. It replaces two executable native suite
  entries while preserving their old IDs as selection aliases, so compatible
  particle-buffer and physics cases require one configured build.
- Expanded 19-case native physics suite:
  - four photochemistry kernel cases for Huebner H2O rates/yield/radial
    scaling, exponential decay probability/event time, product momentum/
    energy and two-body angular isotropy, and the seven-channel branching
    distribution plus neutral-only compiled-product availability;
  - two hard-sphere collision-kernel cases for unequal-mass momentum/energy/
    relative-speed conservation and spherical isotropy;
  - one drifting-Maxwellian distribution case for component means, analytic
    mean thermal speed, and recovered kinetic temperature;
  - six Larsen-Borgnakke cases for optional particle storage, rotational
    equipartition, vibrational harmonic-oscillator means, beta energy
    partition, two-stage RT redistribution/conservation, and the analytic
    `1/Z_rot=1/5` relaxation occurrence rate;
  - three mover cases for neutral ballistic Simple/Boris limits, destination-
    cell insertion, and an analytic specular reflection at a registered sphere;
  - three populated 3-D-cell cases: NTC/LB H2O rotational equilibration with
    momentum/total-energy conservation, the independent Maxwellian H2O/O
    collision rate `n0*n1*sigma*<g>`, and 2,400-particle photolytic parent
    survival with no daughter production.
- Native physics statistics use fixed case-local seeds, documented sample
  counts, analytical expectations, and explicit confidence/tolerance bounds.
- The photoproduct angular test exposed and fixed uniform-polar-angle sampling
  in `GetProductVelocity_2`; the implementation now samples `cos(theta)`
  uniformly for the correct solid-angle measure.
- The neutral-only branching test exposed and fixed a second production defect:
  product IDs were compacted around unavailable charged daughters while their
  velocity vectors retained original channel slots. Both ordinary and
  given-product generators now compact IDs and velocities together, with a
  NaN-sentinel regression for a surviving slot-1 daughter.
- Native-run safety includes a repository-specific advisory lock, an active
  build/run preflight, exact validation of the generated `build` target,
  application selection through `Config.pl`, fixed `MAKEFLAGS=-j16`, and a
  retained configure/build/MPI transcript.
- The native lock now uses the shared checkout's read-only `Config.pl` inode,
  so cooperating jobs on different nodes or with different `TMPDIR` values
  cannot reconfigure the same checkout concurrently. Explicit output reuse
  invalidates the old per-suite receipt before launch; an adapter must publish
  fresh evidence. PASS receipts must contain at least one executed case.
- Exact order-independent case registries detect a missing, renamed,
  duplicated, or unlinked C++ test even when the remaining executable returns
  success, without relying on unspecified cross-translation-unit static
  initialization order.
- Exact `--case` selection and deterministic `--order normal|reverse|shuffle`
  with a recorded shuffle seed are implemented. The harness checks selection,
  aliasing, and deterministic order construction. Final-tree retained
  normal/reverse/shuffle campaigns each pass all 45 cases. All 21
  dependency-light cases also pass through one-case global selections; the 22
  executable native cases are launched as exact one-case MPI processes from
  one clean build, while the two configuration cases execute as individually
  identified adapter cases.
- The global runner accepts `--jobs N` (default `1`) and executes only
  consecutive nonexclusive suites concurrently. Every `exclusive: true`
  native suite has a leading and trailing barrier, so repository-wide
  configuration/build state cannot overlap another suite. Results are restored
  to manifest order; concurrent case-start lines carry their suite ID.
- The runner now streams `[CASE] RUN <ID> <description>` immediately before
  every Python, dependency-light C++, parser, and configured native case. The
  terminal shows case-level progress while full compiler/MPI output remains in
  the suite log; streaming retains timeout enforcement and kills the complete
  POSIX child process group on timeout or user interrupt. The harness exercises
  both paths with a real descendant process.
- The legacy `test/srcUnitTests` entry point now returns GoogleTest's result
  instead of reporting process success after failed assertions.
- `SRC_UNIT_TESTS_MIGRATION.md` maps every legacy source/case to its registered
  replacement or an explicit open item. It records the legacy out-of-bounds RT
  assertion, allocator-mismatched fixture teardown, non-reset statistical
  accumulators, and still-unported individual-weight split/merge profile.
- `array_1d` now has explicit deep-copy/resizing-copy and ownership-transferring
  move semantics. The container tests exposed that the previous implicit copy
  constructor aliased owned storage and could double-delete when copy elision
  was unavailable.
- CRC diagnostic formatting no longer uses overlapping source/destination
  `sprintf`; bounded append formatting preserves termination on long paths.
- The manifest distinguishes smoke, nightly, profile, extended, distributed,
  and application-smoke campaigns. Smoke selects the three dependency-light
  suites, nightly selects all four current suites, and profile/extended select
  the configured native suite. The harness pins these distinct sets. Empty
  campaign selections do not imply implemented coverage; the exclusive native
  suite remains a scheduler barrier wherever selected.
- Framework, suite-adapter, native-fixture, and core-change comments document
  protocol/failure semantics, path containment, portability, concurrency,
  memory ownership, MPI reduction, test oracles, and validity limits. The
  framework/design and configured-application READMEs carry the matching
  maintainer and machine-configuration contracts.

## Not yet complete

Only `src/general/input` is currently marked `covered`. The test harness is
`partial`, and all `partial` and `planned` entries in
`component_inventory.json` remain open.
In particular, the legacy individual-weight split/merge cases, generic optional
particle-layout alignment, full-buffer fatal path, and four non-NTC collision
selectors have not yet been migrated. Existing model/application runners also
remain to be adapted; mesh/MPI/restart/rank-invariance coverage is open.

The PIC and photolytic-reaction entries remain `partial`. Production NTC cell
selection, a cross-species analytic rate, and H2O RT equilibration now have
one-rank populated-cell evidence, but the other collision selectors, unequal
weights/time steps, collision-limiting branch, injection, general sampling,
and split/merge remain open. Movers still need charged uniform-E/B convergence,
outer/refinement crossings, absorption/moving surfaces, and rank migration.
Photochemistry now has sampled branching and 3-D no-daughter parent decay, but
still needs three-body angular qualification and daughter insertion/statistical-
weight bookkeeping. VT/VV relaxation is not claimed because the current VT
implementation explicitly exits as unverified.

The RT occurrence case deliberately qualifies the current LB closure's
hard-coded `Zrot=5`. The input generator emits and initialization validates
`RotationZnumber`, but `LB::RedistributeEnergy` does not currently use that
array. Selecting a cross-species combination rule is a production-physics
decision, so this test campaign records the discrepancy instead of inventing
one. A later correction must set a non-5 input and verify the configured rate.

`HARNESS-TRACEABILITY-001` imports and runs the same exhaustive validation used
by the standalone `validate_traceability.py` command. It expands and audits
every current production translation unit/public contract and validates the
configuration ledger, so structural drift is part of smoke/nightly execution.
The explicit standalone command remains useful for retaining the expanded
per-file JSON receipt.

`make -C src/test strict-coverage` is expected to fail until those items are
implemented and qualified. A passing smoke run must not be described as
complete AMPS coverage.

## Current validation evidence

The source tree registers four executable suites and 45 individual cases: 11
harness, four generic-input, six general-container, and 24 combined native
cases. T0's final-tree normal, reverse, and seed-1729 shuffled `--all`
campaigns each pass 45/45 with matching terminal, JSON, and 45-test JUnit
receipts. Their retained native executables have the identical SHA-256
`446e920d95fd2f174c7cfce8fadede469ed2af7905116a37fb608c3b9529484a`.
The 38-case all-run and split receipts below are historical progression, not
the current qualification evidence.

| Date | Command/evidence | Outcome | Scope |
|---|---|---|---|
| 2026-11-07 | `python3 src/test/validate_traceability.py --output /tmp/amps-core-traceability.json` | PASS | Structural classification: 24 areas, 209 production translation units, 181 public contracts, 12 configuration contracts, 4 suites, 45 registered cases. No behavioral case execution. |
| 2026-11-07 | `python3 src/test/run_tests.py --all --order normal --output-dir /tmp/amps-t0-all-final-normal` | PASS 45/45, 4/4 suites | Final-tree normal order; JSON/JUnit reconcile; 22 native executable cases each ran in a fresh exact-case MPI process. |
| 2026-11-07 | `python3 src/test/run_tests.py --all --order reverse --output-dir /tmp/amps-t0-all-final-reverse` | PASS 45/45, 4/4 suites | Final-tree reverse suite/case order; exact reverse receipt and identical native executable hash. |
| 2026-11-07 | `python3 src/test/run_tests.py --all --order shuffle --shuffle-seed 1729 --output-dir /tmp/amps-t0-all-final-complete-shuffle-1729` | PASS 45/45, 4/4 suites | Final-tree deterministic shuffled suite/case order; includes real timeout/SIGINT cleanup, terminal/JSON/JUnit reconciliation, distinct campaign assertions, effective generated-mode checks, and an explicit retained native preflight transcript. |
| 2026-11-07 | `/tmp/amps-t0-harness-final-complete` | PASS 11/11, 1/1 suite | Final harness-only receipt after campaign-set, interrupt, traceability, and JSON/JUnit reconciliation hardening. |
| 2026-11-07 | `/tmp/amps-t0-individual-<CASE-ID>` for all 21 dependency-light IDs | PASS 21/21 | Each dependency-light case selected and executed alone through the global `--case` interface. Native executable cases use exact one-case MPI subprocesses in every all-order campaign. |
| 2026-11-07 | `/tmp/amps-t0-individual-CONFIG-GENERATOR-001`, `/tmp/amps-t0-individual-CONFIG-INVALID-DECK-001` | PASS 1/1 each | Both adapter-owned configuration cases pass as exact one-case global selections after independent required clean rebuilds. Together with the 21 dependency-light probes and 22 one-case native MPI processes, every registered case has individual-selection evidence. |
| 2026-11-07 | `/tmp/amps-t0-smoke-final-normal` | PASS 21/21, 3/3 suites | Final dependency-light registry in normal order, including exact campaign-set assertions. |
| 2026-11-07 | `/tmp/amps-t0-smoke-final-reverse` | PASS 21/21, 3/3 suites | Same dependency-light registry in reverse suite/case order. |
| 2026-11-07 | `/tmp/amps-t0-smoke-final-shuffle-1729` | PASS 21/21, 3/3 suites | Same dependency-light registry with deterministic shuffle seed 1729. |
| 2026-11-07 | `/tmp/amps-core-native-t0-b` | PASS 24/24, 1/1 suite | Combined native normal order from one clean configured build; includes both Config.pl/generated-deck cases, three particle-buffer cases, and 19 physics cases. |
| 2026-11-07 | `/tmp/amps-core-native-t0-b-reverse` | PASS 24/24, 1/1 suite | Combined native reverse order from one clean configured build; JSON case selection/results reconcile exactly. |
| 2026-11-07 | `test_output/core-tests/run-20261107T072523Z-3412060/` | PASS 38/38, 5/5 suites | Last complete pre-consolidation execution; historical registry only. |
| 2026-11-07 | `/tmp/amps-core-report-final` | PASS 16/16, 3/3 suites | Last retained dependency-light receipt before the 11-case harness expansion. |
| 2026-11-07 | `/tmp/amps-core-physics-audit-final` | PASS 19/19 | Last retained native-physics family receipt before suite consolidation. |
| 2026-11-07 | `/tmp/amps-core-pbuffer-expanded` | PASS 3/3 | Last retained particle-buffer family receipt before suite consolidation. |

The following records are retained as historical bring-up/progression. Their
different case counts are expected for their then-current registries; none
supersedes the consolidated final-tree 45-case evidence above.

## Historical validation evidence

- `/tmp/amps-t0-final-report-reconciliation` correctly failed
  `HARNESS-REPORT-001` while adding direct JUnit coverage because the synthetic
  self-test receipt omitted the real writer's required `seconds` field. The
  fixture was corrected without changing the writer; the immediate rerun at
  `/tmp/amps-t0-final-report-reconciliation-2` and the final 11-case harness
  receipt both pass.

- `python3 src/test/run_tests.py --all` — PASS on 2026-11-07 after both
  required clean configured native builds: 38/38 individual tests and 5/5
  suites passed with 0 FAIL/SKIP/ERROR.  The command printed the final 38-row
  individual result table and the two summaries.  JSON/JUnit receipts and all
  suite transcripts are retained in
  `test_output/core-tests/run-20261107T072523Z-3412060/`; the JSON report
  independently records `individual_counts={TOTAL:38,PASS:38,FAIL:0,SKIP:0,ERROR:0}`
  and suite `counts={TOTAL:5,PASS:5,FAIL:0,SKIP:0,ERROR:0}`.
- `python3 src/test/run_tests.py --tier smoke --jobs 3 --output-dir
  /tmp/amps-core-report-final` — PASS after strengthening the formatter
  regression: `HARNESS-REPORT-001` directly captures and checks the heading,
  per-case states, and typed suite-level fallback; all 16 dependency-light
  cases were printed and recorded, with 3/3 suites and 0 FAIL/SKIP/ERROR.

- `python3 src/test/run_tests.py --tier smoke --jobs 3 --output-dir
  /tmp/amps-core-concurrent-smoke` — PASS: the three nonexclusive smoke suites
  executed concurrently, all 15 individual cases passed, suite-qualified case
  starts were streamed, and aggregate results remained in manifest order.
- `python3 src/test/run_tests.py --test CORE-PIC-PHYSICS-001 --output-dir
  /tmp/amps-core-physics-expanded-4` — PASS after the required clean configured
  build and one-rank MPI run: the exact 19-case registry passed with 0
  FAIL/SKIP/ERROR in 54.629 seconds. JSON/JUnit reports and the full native
  transcript are retained in that directory. This was the expanded
  native-physics evidence for that registry state.
- `python3 src/test/run_tests.py --test CORE-PIC-PARTICLE-BUFFER-001
  --output-dir /tmp/amps-core-pbuffer-expanded` — PASS after the required
  clean rebuild with RT and collision-frequency storage enabled: exact 3-case
  registry, 0 FAIL/SKIP/ERROR in 37.360 seconds.
- `python3 src/test/run_tests.py --tier smoke --jobs 3 --output-dir
  /tmp/amps-core-concurrent-final` — PASS after interruption/scratch/fail-fast
  hardening: 3 suites and 15 cases, 0 FAIL/SKIP/ERROR in 3.007 seconds.
- `python3 src/test/run_tests.py --test CORE-PIC-PHYSICS-001 --output-dir
  /tmp/amps-core-physics-expanded-final` — PASS after the independent audit
  fixes and required clean rebuild: exact 19-case registry, 0
  FAIL/SKIP/ERROR in 54.818 seconds. This supersedes
  `/tmp/amps-core-physics-expanded-4` as the then-current native evidence.
- `python3 src/test/run_tests.py --tier smoke --jobs 3 --output-dir
  /tmp/amps-core-concurrent-final-2` — PASS after the final scheduler timeout,
  interrupt-cleanup, fail-fast, and scratch-identity changes: 3 suites, 15
  cases, 0 FAIL/SKIP/ERROR in 3.015 seconds. This supersedes the earlier
  concurrent smoke receipts.
- `python3 src/test/run_tests.py --test CORE-PIC-PHYSICS-001 --output-dir
  /tmp/amps-core-physics-audit-final` — PASS after all final independent-audit
  fixes and the required clean rebuild: exact 19-case registry, 0
  FAIL/SKIP/ERROR in 55.407 seconds. This was the final native-family evidence
  before consolidation.
- `python3 src/test/run_tests.py --tier smoke --jobs 3 --output-dir
  /tmp/amps-core-smoke-audit-final` — PASS: final scheduler/protocol tree, 3
  suites and 15 cases, 0 FAIL/SKIP/ERROR in 2.987 seconds. This was a prior
  dependency-light receipt.
- `mpiexec -n 1 ./amps --filter PIC-PBUFFER-` — PASS on the executable from
  the final clean native build: 3/3 particle-buffer cases after hybrid-thread
  fixture restoration was added.

- `python3 src/test/run_tests.py --tier smoke --output-dir /tmp/amps-core-test-final-2`
  — PASS: 3 suites, 14 individual cases, 0 FAIL/SKIP/ERROR. JSON and JUnit
  receipts are in that output directory.
- After the detailed code/documentation audit,
  `python3 src/test/run_tests.py --tier smoke --output-dir
  /tmp/amps-core-test-final-3` — PASS: the same 3 suites and 14 individual
  cases, 0 FAIL/SKIP/ERROR. This was the final dependency-light receipt for
  that earlier registry state.
- `AMPS_TEST_OUTPUT_ROOT=/tmp/amps-core-top-level make core_test_nightly`
  — PASS: validated the top-level Make/nightly entry point with the same 14
  cases.
- `SCRATCH=/tmp python3 src/test/run_tests.py --config
  src/test/config/machine.example.json --test CORE-HARNESS-001 --output-dir
  /tmp/amps-core-machine-profile-2` — PASS and verified that transient work is
  routed through the configured scratch root.
- `python3 src/test/run_tests.py --coverage --require-complete-coverage` —
  expected exit 1: 2 complete and 22 incomplete inventory components.
- `python3 src/test/run_tests.py --test CORE-PIC-PARTICLE-BUFFER-001
  --output-dir /tmp/amps-core-native-final-3` — PASS after a clean configured
  native build and one-rank MPI run: 1 suite, 3 cases, 0 FAIL/SKIP/ERROR. The
  JSON/JUnit reports, full build log, and particle-buffer transcript are in
  `/tmp/amps-core-native-final-3`.
- `python3 src/test/run_tests.py --test CORE-PIC-PHYSICS-001 --output-dir
  /tmp/amps-core-physics-5` — PASS after a clean configured build and one-rank
  MPI run: 1 suite, 11 cases, 0 FAIL/SKIP/ERROR. JSON/JUnit reports and the
  complete native transcript are retained in that directory.
- `python3 src/test/run_tests.py --test CORE-PIC-PARTICLE-BUFFER-001
  --output-dir /tmp/amps-core-pbuffer-physics-final` — PASS after a second
  clean configured build: the pre-existing exact 3-case registry remains green
  with the expanded species and optional particle-data layout.
- `python3 src/test/run_tests.py --tier smoke --output-dir
  /tmp/amps-core-smoke-physics-final` — PASS: 3 dependency-light suites, 14
  cases, 0 FAIL/SKIP/ERROR after registering the physics suite.
- `python3 src/test/run_tests.py --all --list` and `python3
  src/test/run_tests.py --coverage` — PASS registry/inventory audit: 5 suites,
  24 components, 2 covered, and 22 honestly partial/planned.
- `python3 src/test/run_tests.py --tier smoke --output-dir
  /tmp/amps-core-case-progress-final` — PASS: all 14 cases were announced by
  stable ID and description before execution, with unchanged JSON/JUnit
  results.
- `python3 src/test/run_tests.py --test CORE-PIC-PARTICLE-BUFFER-001
  --output-dir /tmp/amps-core-case-progress-native` — PASS after the required
  clean native build: all 3 MPI-launched cases were announced live and the
  exact registry remained green.
- `git diff --check` — PASS after the implementation and documentation audit.
- `PYTHONPYCACHEPREFIX=/tmp/amps-core-pycache python3 -m compileall -q
  src/test` — PASS for every Python runner/framework module without writing
  bytecode into the source tree.

## Native bring-up record

The first clean native attempt (`/tmp/amps-core-native-final`) failed while
compiling `meshAMR`: the new input deck had not invoked `exosphere.pl`, which
is the established generator for AMPS' default internal-boundary user-data
classes. The fix adds a target-free, SPICE-free parser block. It supplies type
definitions only; exosphere physics and source modules remain excluded.

The second clean attempt (`/tmp/amps-core-native-final-2`) compiled and linked
successfully but reported no cases. AMPS stores application objects in the
static archive `main.a`, and the linker did not load `particle_buffer.o` when
it contained only static registration constructors. The fixture now exports
`LinkParticleBufferTests()`, called by `main.cpp`, to make archive membership
deterministic. The third clean attempt passed the exact three-case registry.
Both failed attempts remain useful diagnostic evidence and were not described
as test passes.

## Native physics bring-up record

- `/tmp/amps-core-physics-1` failed during input generation because the
  collision table used a single backslash instead of AMPS' double-backslash
  continuation token.
- `/tmp/amps-core-physics-2` reached the IDF parser and rejected unsupported
  `RotZnum all` shorthand. The unused rotation-collision and pair-index values
  were removed because every relaxation switch is off in this kernel fixture.
- `/tmp/amps-core-physics-3` built and passed all nine non-mover cases. The
  mover fixture incorrectly selected an unallocated AMR internal/root node;
  the second case inherited it and triggered the expected signal handler.
- `/tmp/amps-core-physics-4` resolved the owner leaf and produced correct
  ballistic states, but inspected the starting leaf's temporary list for an
  endpoint exactly on a leaf face. The fixture now avoids an ambiguous face
  endpoint and resolves the destination owner explicitly.
- `/tmp/amps-core-physics-5` passed the exact 11-case registry. None of the
  failed bring-up runs is counted as validation evidence.
- `/tmp/amps-core-physics-expanded-1` executed all 19 expanded cases and
  reported 16 PASS / 3 FAIL. The failures were test-fixture/oracle defects:
  the neutral-only build intentionally omits unavailable ionic daughters, the
  RT oracle initially used one molecule's rotational exponent rather than the
  pair's six rotational degrees, and the Boris free-flight trajectory crossed
  the newly registered sphere. All three definitions were corrected without
  relaxing a production tolerance.
- `/tmp/amps-core-physics-expanded-2` was manually interrupted because the
  240,000-draw branching check also generated full two/three-body kinematics
  on every draw and was unsuitable for routine execution. Its orphaned adapter
  and MPI processes were explicitly cleaned up. The final test separates the
  high-statistics production channel selector from a smaller product-contract
  sample; product kinematics retain their dedicated conservation case.
- `/tmp/amps-core-physics-expanded-3` produced a typed `ERROR` because the
  interrupted adapter still held the repository-specific native lock. No test
  result from that invocation is counted as scientific validation.
- `/tmp/amps-core-physics-expanded-4` passed the final exact 19-case registry
  in 54.629 seconds. It supersedes the three diagnostic expanded bring-up runs;
  those failures/interruption remain recorded and are not described as passes.
