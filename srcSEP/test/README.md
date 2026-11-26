# srcSEP standalone component tests

## CME shock evolution regression

Legacy case `SEP--CME-shock-evolution` checks the Python fitter's rejection
of the archived 29 November 2020 CME at the spherical model's width bound.
It runs both the observation preview and packing scripts with their own
offline references before the native test.
It also runs the linked `./amps --test CME` on the manufactured event. The
checker uses an independent trajectory oracle for the initial shock map
and fifteen 20-minute outputs, with a frozen five-hour numerical reference.
See the [case guide](cme_shock_evolution/README.md) and
[observation package](cme_shock_evolution/observed_2020_11_29/README.md).

## 2013-04-11 real-event tests

Legacy cases `SEP--event-2013-04-11--Parker`, `--FTE-Dmumu` and `--FTE-MFP`
transport protons along the field line through Earth for the SEP event of
11 April 2013 with each production mover. The shock, solar wind and field
come from real observations (Lario et al. 2014 shock heights, the Wind
shock, OMNI plasma, the WSO magnetogram); the comparison data are SOHO/ERNE,
ACE/EPAM and GOES-13. The nightly cases run the first hour on a coarse grid
and compare four products with DiffNum against AMPS references; the
full-event validation tier is proposed separately. See the
[case guide](event_2013_04_11/README.md).

## Stage 13 spherical SWCME3D

`make -C srcSEP test-stage13-cases SEP_EXECUTABLE=/absolute/build/amps`
checks the canonical spherical core, native `.in` selection, boundary values
and observation package integrity. `make -C srcSEP test` invokes it in the
routine sequence and executes linked adapter case `S13AD01`. `S13SC11` and
`S13SC12` cover spatial gradients and the fixed-step material map. The linked
case checks both profile identities and the actual local species states. The [Stage 13
campaign guide](../validation/STAGE13_SPHERICAL_VALIDATION.md) gives the
commands for a single observation or a full campaign.

## Stage-8 backend admission and stage-9 full application

The [support matrix](../validation/STANDALONE_BACKEND_SUPPORT.md) and
[test contract](../validation/STAGES8_9_TEST_CONTRACT.md) distinguish diagnostic
exchange/rejection from actual supported reduced-profile propagation. Run:

```sh
make -C srcSEP test-backend-admission-unit test-standalone-application-runner-unit
python3 -B srcSEP/test/run_backend_admission_tests.py --native --amps /absolute/build/amps --output-dir /fresh/backend-native
python3 -B srcSEP/test/run_standalone_application.py --build-root /absolute/build --output-dir /fresh/application400
```

The latter requires a configured proton reduced-front field-line build whose
generated source, compiler dependencies, objects/archive and standalone main
match this checkout. It copies an immutable executable and runs all three
movers for 400 complete 300-s physical steps. Physical input, finite release
distance, prescribed coefficients, source interpretation, absorbing return,
1800-s per-run deadline and 1e-10 budget criterion are the existing gate.
Missing build/input/output/MPI evidence is ERROR, never a short-run PASS.

Where MPI must be launched directly by the host, use `--prepare-only`, execute
the exact jobs in `launch-plan.json`, preserve each real launcher exit in its
`process-result.json` with matching `argv` and `cwd`, then use `--collect`.
Preparation performs no physical test. The native backend runner supports the
same distinction, using `prepared.json` and its named per-job receipt paths.
Never infer or supply a successful return code without actual execution.

Shared suites are `CORE-SEP-BACKEND-ADMISSION-001` (SBACK01–05),
`CORE-SEP-BACKEND-NATIVE-001` (SBACK06–08) and
`CORE-SEP-FULL-APPLICATION-001` (SAPP01–03). Their individual readers require
the complete catalog and preserve original FAIL/ERROR/SKIP statuses. Existing
source-suite aliases expose `backend-admission`, `backend-admission-native`,
`standalone-application-runner` and `standalone-application-native`.
The orchestrator's `standalone-application-runner` step requires all eleven
software controls individually, including preserved existing evidence and
complete malformed-plan ERROR reports. Expected child diagnostics are captured
and asserted as fixture data, so they cannot be misread as unittest verdicts.
Case identities are preserved across Python's class-only and fully qualified
method output formats.

```sh
python3 -B tools/sep_test_orchestrator.py --app srcsep --skip-build --only backend-admission --only backend-admission-native --only standalone-application-runner --only standalone-full-application --output-dir /fresh/sep-gates
```

This orchestrator command runs from the matching configured checkout. The
full application is a long step; `--quick` intentionally omits it and cannot
satisfy stage-9 full-output qualification. Existing native/material/thread
restart, core suites and four AMPS–BATS-R-US profiles remain separate required
campaign roles. No regression reference is regenerated by these procedures.

The correctness amendment adds source sanitizer gates selectable through
`make -C srcSEP test-rng-lifecycle-unit`,
`make -C srcSEP test-background-derivatives-unit`, and
`make -C srcSEP test-wave-normalization-unit`. The exact descriptor callbacks
are also linked into the native catalog under `rng-lifecycle`,
`background-derivatives`, and `wave-normalization`; their pure-kernel evidence
does not certify native particle storage or restart. `POPN01` and `POPN02`,
group `population-native`, exercise actual PIC field-line splitting and
byte-exact particle cloning. All four groups are reachable from the shared
orchestrator. Source and native execution records retain their evidence level.

Current requirements and unresolved qualification gates are maintained in
[`../validation/development_requirements.json`](../validation/development_requirements.json).
External reference locators and their limitations are in
[`../validation/reference_sources.json`](../validation/reference_sources.json).

## Build and run everything: `run_all_test.sh`

`test/run_all_test.sh` performs a clean native srcSEP build (AGENTS.md
preflight, `rm -rf -- build`, `./Config.pl
-application=test/sep_parker_spiral__field_line`, `make -j`) and then runs every
srcSEP test procedure recorded in `srcSEP3D/LOG.md`: `test/run_tests.py --all`
against the new executable (its per-test logs are cut from `test-run.log`),
`test/run_parallel_diffusion_mover_tests.py` (per-suite logs), the generic
input-reader test (`amps-input`: `make -C src/general/input test`), and the
orchestrator self-test. It replaces whatever application `build/` held before.
The script is a thin wrapper around the shared orchestrator
`tools/sep_test_orchestrator.py --app srcsep`.

Each step's complete terminal output goes to `OUTPUT_DIR/NN-step.log`
(`--stream` also shows it live), and everything printed on the terminal is
appended to the run log `OUTPUT_DIR/run_all_test.log` (or `--log FILE`).

**Summary.** After each step, a *result reader* in
`tools/sep_test_orchestrator.py` interprets the report that step's runner
already writes (JSON where available, otherwise the runner's stable text
format) and the script prints, per step and in total, how many tests ran,
passed, failed, were skipped, and errored. "Test executions" is the plain sum;
"unique tests" counts a test that several steps run (same namespace and ID)
once, with its worst status. The build and initialization previews are
counted separately as checks.

**Failed tests.** Every failed or errored test is listed with its name,
status, message, the runner (script and step) that ran it, and its test log:
the runner's own per-test log where one exists, otherwise the test's section
copied from the runner's shared log (or its captured output) into
`OUTPUT_DIR/failed-tests/<step>/<test>.log`, otherwise the step log and line
number (several native tests in one MPI process cannot be separated). The same
list is written to `OUTPUT_DIR/failed_tests.txt`; all records and counts go to
`OUTPUT_DIR/run_all_test.json`, and `--junit FILE` adds JUnit XML. A step
whose report is missing counts as one ERROR (`<step>:results-unavailable`), and
a step that exits nonzero without a failing test record adds an ERROR
(`<step>:exit-status`), so a failure can never look green. The script exits 0
only when every check passed and no test failed or errored.

`tools/test_sep_test_orchestrator.py` (step `orchestrator-selftest`) tests every
result reader, the log extraction, the step guards, and the summary.

```sh
srcSEP/test/run_all_test.sh --help        # options and a description of each step
srcSEP/test/run_all_test.sh               # build + all tests
srcSEP/test/run_all_test.sh --list        # steps and exact commands
srcSEP/test/run_all_test.sh --stream --log my-run.log   # live test output, custom run log
srcSEP/test/run_all_test.sh --skip-build --only pd-mover --junit pd.xml
```

> **Required for new tests:** when a new srcSEP test, test target, or runner
> is developed, add it as a step to the srcSEP step table in
> `tools/sep_test_orchestrator.py` in the same change, with a result reader
> for its report and the location of its per-test logs (and a fixture in
> `tools/test_sep_test_orchestrator.py` for a new report format); the runner's
> output must follow the contract in the next subsection. A test not
> reachable from `run_all_test.sh` is not part of the srcSEP regression series;
> a step without a reader can only count as a single check.

### What a new test runner must output for `run_all_test.sh`

`run_all_test.sh` does not re-implement any test: it runs each step's command
from the AMPS root, captures its stdout and stderr in the step log
`OUTPUT_DIR/NN-<step>.log`, and then a *result reader* in
`tools/sep_test_orchestrator.py` turns the runner's output into one record per
test. A new runner (or native test target) is understood only if its output
follows the contract below. Any output that does not fit is counted as ERROR,
never as PASS.

**1. One record per test, with a stable unique ID and a status.**

- *Test ID*: a stable, unique name such as `CLI07` or `PDB06`. It must not
  change between runs, because de-duplication and the failed-test list are
  keyed on it. Use only `A-Z a-z 0-9 _ . - : /`.
- *Status*: one of `PASS`, `FAIL`, `SKIP`, `ERROR`. The spellings `passed`,
  `ok`, `failed`, `failure` and `skipped` are also accepted (case does not
  matter). Any other word becomes ERROR.
  - `FAIL`: the test ran and an assertion or physics criterion was violated.
  - `ERROR`: the test could not be evaluated, for example because the program
    crashed, timed out or wrote no output.
  - `SKIP`: a declared precondition is missing, such as a data file, a tool,
    MPI or an unbuilt target.
  A missing input must never be reported as `PASS`.
- *Message*: a one-line reason. It is required for `FAIL`, `ERROR` and `SKIP`;
  for `SKIP` it should name the missing precondition. Longer messages are
  collapsed to one line of at most 300 characters.

**2. Preferred format: a JSON report at a path given on the command line.**

The runner takes an output location (`--json FILE` or `--output-dir DIR`) and
writes the report there on every run. It writes the report even when tests
fail, and it writes it after the last test so that the report is complete.
When the report path is fixed rather than given on the command line, the
reader must reject a report older than the step; `reader_coronal_cme` shows
how. The shape below is read unchanged by the existing
`reader_reduced_front(DIR)` (it reads `DIR/summary.json`; despite its name it
is format-generic), so a runner that writes it needs no new reader code:

```json
{
  "results": [
    {"scope": "", "test_id": "CLI07", "status": "PASS", "message": "",
     "log": "/abs/path/OUTPUT/logs/CLI07.log", "report": ""},
    {"scope": "", "test_id": "CLI08", "status": "FAIL",
     "message": "expected mover parker, got fte-dmumu",
     "log": "/abs/path/OUTPUT/logs/CLI08.log",
     "report": "/abs/path/OUTPUT/CLI08.json"}
  ]
}
```

- `scope` (optional) prefixes the ID as `scope/test_id`. Use it when the same
  ID is legitimately run in several configurations, such as 1 rank and
  4 ranks.
- `log` is the test's own log: its complete output and the verdict, written
  for every test and not only for failures. It may be an absolute path or a
  path relative to the AMPS root (steps run from the root); absolute paths are
  safest. It is printed in the failed-test list. When it is empty, the failed
  list falls back to the step log.
- `report` (optional) is a per-test evidence file, listed as evidence.
- A per-suite log shared by several tests is acceptable if each result also
  gives the 1-based `line` of its result line, as
  `run_parallel_diffusion_mover_tests.py` does (`reader_pd_mover`).

**3. Accepted text formats (when a JSON report is impractical).** The reader
then works on the step log, so each test's record points to its line in that
log instead of a separate test log.

- One line per test on stdout, at the start of the line:
  `[ID] PASS|FAIL|SKIP|ERROR message`, or `ID PASS|FAIL|SKIP|ERROR: message`
  (`read_bracket_text`). If the same ID is printed more than once, the worst
  status is kept.
- Several test binaries run by one make recipe: run them with `make -i` and
  use `reader_make_recipe(r"^build/(test_\w+)$")`. Each binary is then counted,
  and a binary that exits nonzero without a `FAIL`/`ERROR` line becomes an
  ERROR.
- Python `unittest` run with verbosity 2 (`read_unittest`). Each failure
  block is copied to that test's own log.
- A native AMPS registry test run as `./amps --test ID --test-json FILE`
  (`reader_native_json`).

**4. Exit status.** The runner exits nonzero when any test FAILed or ERRORed,
and exits 0 otherwise; SKIPs alone do not make the exit nonzero. If the runner
exits nonzero but its report contains no failing test, `run_all_test.sh` adds
the ERROR `<step>:exit-status`. If the report is missing, unreadable or empty,
it records the ERROR `<step>:results-unavailable`. A crashed runner can
therefore never look green, but it also cannot say which test failed.

**5. Registration.** Add a `Step(...)` to this application's step table in
`tools/sep_test_orchestrator.py`, together with any directories the runner
needs (`make_dirs`). For example:

```python
Step("my-runner",
     "one-line description shown by --help and --list",
     [["python3", "srcSEP/test/run_my_tests.py", "--output-dir", str(out / "my-runner")]],
     reader=reader_reduced_front(out / "my-runner"),   # JSON report of section 2
     runner="srcSEP/test/run_my_tests.py",
     long=False),   # True if it takes many minutes; such steps are skipped by --quick
```

Give `namespace=` only if the step re-runs tests that another step also runs,
so that they are counted once as unique tests. A new report format needs its
own reader plus a fixture test in `tools/test_sep_test_orchestrator.py` that
uses a hand-written report rather than one produced by the runner. Then check:

```sh
python3 tools/test_sep_test_orchestrator.py
srcSEP/test/run_all_test.sh --list
srcSEP/test/run_all_test.sh --skip-build --only my-runner
```

The new step's test count must equal the runner's own total. A deliberately
failing test must appear in the failed-test list with its runner and log.

## Run the complete test series

From the `srcSEP` directory, run every test advertised by the linked AMPS
application with:

```sh
test/run_tests.py --amps ../amps --all --output-dir test_output/all
```

When registry source has changed, rebuild the configured enclosing application
before discovery with the single command:

```sh
env MAKEFLAGS="-j16" test/run_tests.py --amps ../amps --all \
  --amps-source .. --make-config ../Makefile.conf --rebuild \
  --output-dir test_output/all
```

`MAKEFLAGS` affects the enclosing clean/build only. Registry callbacks still
run one isolated process at a time so process-global AMPS/SWCME state and each
test's log/report remain deterministic.

The runner discovers the complete native registry, executes every test in an
isolated process, and automatically uses the registered input/reference
workflow for CV, IV, and XM validation cases. A failed, errored, or crashed
test does not prevent later tests from running. The final output reports the
PASS/FAIL/SKIP/ERROR totals and lists every failed or errored test with its
diagnostic. Detailed JSON, JUnit, logs, plots, and per-test artifacts are saved
under `test_output/all`.

Improvements D01 and D02 and the bounded D03 native preflight are now public
C++ registry entries named `D01`, `D02`, and `D03PRE`. Consequently the command
above discovers and executes them automatically; selecting a Make suite is no
longer required to make their results visible. They are marked `extended`
because they reconfigure process-owned SWCME state, so `--routine`/native
`--all-tests` excludes them while the isolated Python `--all` includes them.
`D03PRE` verifies only the linked mover/refresh evidence hooks. The configured
serial/MPI/restart D03 campaign remains the separate external gate documented
below and must pass before a release can claim D03 completion.

Step 1 provides one catalog and result contract for standalone component tests.
The catalog lives in `component_tests.cpp`; generic deterministic selection,
metadata validation, timing, status handling, and output formatting live in
the canonical `src/models/sep_common/sep_test_registry.*`.  The production CLI
remains application-owned in `util/sep_cli.*`.

## Running tests

Use the linked executable directly or the Make wrappers:

```sh
# Dependency-light parser/registry contract test; does not need AMPS or MPI.
make test-cli-unit

# Dependency-light immutable-background and single-clock contract tests.
make test-state-unit

# Dependency-light SI flux-tube geometry and source-normalization tests.
make test-geometry-source-unit

# Dependency-light three-mover registry and CLI contract tests.
make test-mover-api-unit

# Source-scope checks for the field-line-only Step 5 boundary.
make test-field-line-scope-unit

# Step 6 shared state, coordinate, cooling, limits, and RNG kernels.
make test-transport-common-unit

# Step 7 canonical Itô Parker solver.
make test-parker-unit

# Step 8 coefficient-driven focused-transport Dmumu solver.
make test-fte-dmumu-unit

# Step 9 event-driven focused-transport mean-free-path solver.
make test-fte-mfp-unit

# Step 10 coefficient registries, validation, and SI conversions.
make test-coefficients-unit

# WP11-WP20 bounded diffusion, error controls, repaired coefficients,
# species-source normalization, named scales, and ballistic compatibility.
make test-wp11-wp20-unit

# WP21-WP30 shock/source/flux/sampling/output/configuration contracts.
make test-wp21-wp30-unit

# WP31-WP41 stiffness, population, system, statistical, robustness,
# observation, performance, and evidence-governance contracts.
make test-wp31-wp41-unit

# B05 experimental WP42-WP64 mathematical/component contracts. These
# implementations are intentionally not linked into the production archive.
make test-wp42-wp64-experimental

# Native/external WP59, WP60, WP63, and WP64 gates. Without an approved site
# command this produces a first-class SKIP record, never a synthetic PASS.
make test-wp59-wp64-native

# Print the exact 90-row production compatibility matrix from its registry.
make print-configuration-matrix

# Step 11 authoritative turbulence state, ledger, remap, and restart.
make test-turbulence-core-unit

# Step 12 deterministic worker/rank reduction and keyed RNG evidence.
make test-reproducibility-unit

# Step 13 registered background/cross-mover fixtures and JSON/JUnit evidence.
make test-acceptance-unit

# Step 14 docs, public-symbol, archive-hygiene, warning, and analyzer gates.
make test-documentation-unit

# Step 15 numerical, cross-mover/model, SWCME replay, and release-report gates.
make test-scientific-validation

# Canonical SWCME ownership plus srcSEP/build-main path equivalence.
make test-swcme-relocation-unit

# D01 strict SWCME query/preparation status, recovery policy, and counters.
make test-swcme-fail-closed-unit

# D02 canonical presets, units, layers, source/event config, and fingerprint.
make test-swcme-configuration-unit

# D03 manifest/coverage self-test. The configured native command is below.
make test-d03-runner-unit

# Canonical seven-member SEP-common archive, one-definition, consumer-link,
# srcSEP/build-main path equivalence, sibling-AMPS public-header fallback, and
# AMPS-root-only production-source include coverage (including field_line.cpp).
make test-sep-common-ownership-unit

# Dependency-light analytical component catalog. CV01 runs additionally when
# SEP_EXECUTABLE identifies the linked application.
make test-controlled-analytical

# Test 01 / CV01 strict compile gate plus linked end-to-end execution.
make test-cv01-unit SEP_EXECUTABLE=/path/to/amps

# CV02-CV05 strict source/registry gate plus linked end-to-end execution.
make test-cv02-cv05-unit SEP_EXECUTABLE=/path/to/amps

# CV06-CV12 strict source/registry gate plus optional linked execution.
make test-cv06-cv12-unit SEP_EXECUTABLE=/path/to/amps

# IV01-IV06 integrated manufactured source gate plus linked execution.
make test-iv01-iv06-unit SEP_EXECUTABLE=/path/to/amps

# XM01-XM03 cross-model/observation source gate plus linked execution.
make test-xm01-xm03-unit SEP_EXECUTABLE=/path/to/amps

# OV01-OV05 observational registry/fixed-input gate. OV01-OV02 are release
# gates; OV03-OV05 remain diagnostic. Linked execution is optional here.
make test-ov01-ov05-unit SEP_EXECUTABLE=/path/to/amps

# EV01-EV02 training/holdout pilot registration, dispatch, provenance, and
# source-history checks. These cases intentionally remain pilot-incomplete.
make test-ev01-ev02-unit

# Cross-file Stage 3 manifests, retired files, OV/EV roles, and srcSEP3D
# complete compiled-species ownership. This gate needs neither AMPS nor MPI.
make test-stage3-contracts

# Strict four-component source-tree audit and negative controls.
make test-package-hygiene-unit

# All dependency-light ASan/UBSan suites from Steps 3 and 6–13.
make test-sanitizer

# Repeat keyed campaigns and compare physics evidence after timing normalization.
make test-stochastic-repeat

# Discover the catalog without initializing the model.
make test-list

# Run one stable ID, a group, or all bounded routine tests.
make test-case CASE=DXX01
make test-group GROUP=parker
make test-turbulence
make test-stress
make test-mpi MPI_NP=4
make -j test
```

The makefile expects the linked executable at `../amps`.  If the enclosing AMPS
application writes it elsewhere, add
`SEP_EXECUTABLE=/absolute/path/to/executable`.  `CASE` and `GROUP` are required
for their respective targets, and any child failure is returned by Make.

Equivalent CLI examples are:

```sh
../amps --list-tests
../amps --test DXX01
../amps --test=DXX01 --test=TURB01
../amps --test-group parker
../amps --test-group fte-dmumu
../amps --test-group fte-mfp
../amps --test PARK07 --test FTEM08
../amps --test TURB21 --test TURB22 --test TURB23
../amps --test D01
../amps --test D02
../amps --test D03PRE
../amps --test-group swcme-background
../amps --test-group swcme-configuration
../amps --test-group native-integration
../amps --all-tests
../amps --all-tests --test-json results.json --test-junit results.xml
```

`make -j test` intentionally sequences shared-state component execution even
when Make is given `-j`, then runs the canonical
`src/models/swcme/test/run_tests.py --routine` suite. SWCME is not embedded in
srcSEP. Expensive extended tests are not part of `--all-tests`.

The dependency-light D01 gate can also be selected through the Python runner:

```sh
test/run_tests.py --suite d01-background --output-dir test_output/d01
```

This source-only gate compiles the production private adapter and the same
`SwcmeImprovementDescriptors()` factory linked into the application against
the canonical `src/models/swcme` implementation. It covers unprepared state,
boundary and sub-boundary radii, NaN/infinity, negative density, zero speed,
transactional failed preparation, strict/clamp/fallback policy semantics,
recovery counters, complete rank/epoch/location/state diagnostics, and the
unchanged fast/slow preset trajectories. It does not substitute for the D03
native MPI gate, which verifies cross-rank agreement in a configured AMPS run.

The D02 gate is selectable with:

```sh
test/run_tests.py --suite d02-configuration --output-dir test_output/d02
```

It verifies fast/slow equivalence, supported unit conversion, input-to-CLI
precedence, duplicate/unknown/unsupported fields, incompatible shock modes,
data-driven tables, deterministic normalized fingerprints, parser-free host
construction, validity intervals, and explicit rejection of an absolute source
normalization not yet implemented by srcSEP.

The three native entries can be checked directly after rebuilding the linked
application:

```sh
../amps --test D01 --test-json d01.json --test-junit d01.xml
../amps --test D02 --test-json d02.json --test-junit d02.xml
../amps --test D03PRE --test-json d03pre.json --test-junit d03pre.xml
```

`D01` and `D02` run the identical callbacks used by the corresponding focused
Make targets. `D03PRE` checks the exact three canonical movers, their
coefficient contracts, a nonempty canonical SWCME fingerprint, and two
monotonic prepared-state generations. A PASS means the linked executable is
ready to emit the evidence consumed by the outer campaign; it does not claim
that MPI decomposition or checkpoint/restart equivalence was executed.

The D03 gate uses the same runner entry point but requires an enclosing AMPS
configuration and, for execution rather than a development SKIP, a reviewed
site manifest:

```sh
# Development: missing native prerequisites are retained as a structured SKIP.
env MAKEFLAGS="-j16" test/run_tests.py \
  --suite d03-native-integration \
  --amps-source .. --make-config ../Makefile.conf \
  --output-dir test_output/d03

# Release: clean build and complete site campaign are mandatory.
env MAKEFLAGS="-j16" test/run_tests.py \
  --suite d03-native-integration \
  --amps-source .. --make-config ../Makefile.conf \
  --native-manifest test/native_integration_manifest.site.json \
  --rebuild --release --output-dir test_output/d03-release
```

`test/native_integration_manifest.example.json` is a coverage template. Copy it
to a site file, remove `"template_only": true`, and replace all `--site-input`,
`--checkpoint`, and artifact placeholders with the actual AMPS post-compile
input/restart convention. It is intentionally not runnable unchanged. The
validator requires all three movers,
serial plus two distinct MPI decompositions for each, a source state ID of at
least two, at least one completed particle dispatch per case, one exact output
group spanning each mover's decompositions, and uninterrupted/checkpoint/resume
coverage with a shared exact final-output group.

The native executable, not the Python orchestrator, supplies the physics
evidence. At shutdown it compares SWCME state identity and epoch across ranks,
reduces background diagnostics and completed mover dispatches, and prints an
explicit consensus marker. The gate additionally compares configuration
fingerprints and SHA-256 artifacts. Full logs, commands, compiler/MPI versions,
checksums, timings, case results, and equivalence groups are saved in
`native_integration_evidence.json`. See
[`../NATIVE_INTEGRATION.md`](../NATIVE_INTEGRATION.md) for the complete evidence
and claim contract.

```sh
test/run_tests.py --suite swcme-relocation \
  --output-dir test_output/swcme-relocation
```

In the normal AMPS layout it discovers `../src/models/swcme` automatically. A
detached source checkout may set `SWCME_DIR=/absolute/path/to/swcme`; this is an
explicit development override, never a fallback to a directory within srcSEP.

The Step 6–15, WP11–WP20, WP21–WP30, and WP31–WP41 focused targets do not
require the linked executable, AMPS, PIC, MPI, SWMF, or field-line host classes.
Each compiles the exact production numerical cores with C++11,
`-Wall -Wextra -Werror -pedantic`, AddressSanitizer, and
UndefinedBehaviorSanitizer. The B05 WP42–WP64 focused target uses the same
compiler discipline but compiles a deliberately experimental overlay that is
absent from `MAINLIBOBJ`; its result must not be described as production-path
coverage. LeakSanitizer alone is disabled because the managed test environment
does not expose the required `/proc` task data.

## Python runner and analytical-comparison plots

`run_tests.py` provides one Python 3 entry point for native registry tests and
dependency-light Make suites. It does not contain test callbacks, movers, or
analytical physics: stable IDs/groups are discovered from the linked executable
and execution status comes from the registry JSON. Matplotlib is the only
additional runtime dependency for figures.

Native discovery, selected groups, and each `--all` case now run in unique
`native-work/<label>-<unique>/` directories below the requested evidence root.
This includes linked CV/IV/XM/OV/EV realizations and the orchestrator's PNW
group. PIC `PT/thread=*.log`, mesh caches, and incidental files belong to that
directory; concurrent campaigns do not share source-tree scratch names.
Source-only Make suites keep their existing `srcSEP` working directory.

Python-owned `--amps`, `--output-dir`, and `--case-input` paths remain relative
to the caller. Forwarded native path arguments retain their historic `srcSEP`
base and become absolute before launch; no existence-based fallback search is
performed. The reviewed flat native CV/controlled argument manifests and their
`--initial-particles` / `--source-history-csv` leaf inputs are copied byte for
byte into the owned run, then only those declared leaf references are rewritten
to the copied absolute paths. `native-execution.json` records requested and
resolved arguments, the original base, retained input paths, and SHA-256 hashes.
Unknown nested manifest schemas fail explicitly and require a reviewed adapter.
Versioned production `--input` retains its own parser's nested-path contract;
the component-test CLI currently rejects combining it with a test selection.
Bare `--list` uses a disposable temporary working directory. Campaign discovery
retains `discovery.log` and its working directory with the other evidence.

The runner isolation regression uses an actual child that writes PIC-style
scratch before both discovery and execution. It verifies that a source-tree
sentinel and the original relative manifest/CSV remain unchanged, all three
native working directories differ, and the linked portfolio launcher also
owns its discovery and execution scratch. These checks validate infrastructure,
without implying a physical model qualification:

```sh
python3 -B srcSEP/test/test_python_test_runner.py
python3 -B tools/test_sep_test_orchestrator.py
```

Run `python3 test/run_tests.py --help` for annotated, copy-and-paste examples of
every selection mode, MPI execution, report replotting, plot controls, and
forwarding model-specific arguments to the linked executable.

```sh
# Run the first numbered validation case through the linked application.
python3 test/run_tests.py --amps /path/to/amps --validation-case CV01 \
  --output-dir /evidence/srcsep/CV01

# Repeat --validation-case to run any subset of the linked portfolio.
python3 test/run_tests.py --amps /path/to/amps \
  --validation-case CV02 --validation-case CV03 \
  --validation-case CV04 --validation-case CV05 \
  --output-dir /evidence/srcsep/CV02-CV05

# Run CV06-CV12 through the same linked executable and evidence contract.
python3 test/run_tests.py --amps /path/to/amps \
  --validation-case CV06 --validation-case CV07 \
  --validation-case CV08 --validation-case CV09 \
  --validation-case CV10 --validation-case CV11 \
  --validation-case CV12 --output-dir /evidence/srcsep/CV06-CV12

# Run integrated manufactured geometry/operator/coupling cases.
python3 test/run_tests.py --amps /path/to/amps \
  --validation-case IV01 --validation-case IV02 \
  --validation-case IV03 --validation-case IV04 \
  --validation-case IV05 --validation-case IV06 \
  --output-dir /evidence/srcsep/IV01-IV06

python3 test/run_tests.py --amps /path/to/amps \
  --validation-case XM01 --validation-case XM02 --validation-case XM03 \
  --output-dir /evidence/srcsep/XM01-XM03

# Registered observational comparisons. Their publication-owned defaults are
# immutable, so --case-input is intentionally rejected for these IDs.
python3 test/run_tests.py --amps /path/to/amps \
  --validation-case OV01 --validation-case OV02 \
  --validation-case OV03 --validation-case OV04 --validation-case OV05 \
  --output-dir /evidence/srcsep/OV01-OV05

# Calibration/held-out ensemble pilots. These emit campaign evidence but do
# not acquire release-gate status when they pass.
python3 test/run_tests.py --amps /path/to/amps \
  --validation-case EV01 --validation-case EV02 \
  --output-dir /evidence/srcsep/EV01-EV02

# List all native registry IDs without initializing AMPS.
python3 test/run_tests.py --amps /path/to/amps --list

# Select arbitrary IDs and/or groups. Repeat either selector as needed.
python3 test/run_tests.py --amps /path/to/amps \
  --test PARK01 --test FTED08 --group fte-mfp \
  --output-dir /evidence/srcsep/selected

# Bounded routine registry versus every discoverable routine/extended case.
python3 test/run_tests.py --amps /path/to/amps --routine \
  --output-dir /evidence/srcsep/routine
python3 test/run_tests.py --amps /path/to/amps --all \
  --output-dir /evidence/srcsep/all

# Rebuild the configured binary before discovering the native catalog.
env MAKEFLAGS="-j16" python3 test/run_tests.py --amps /path/to/amps --all \
  --amps-source /path/to/AMPS --make-config /path/to/AMPS/Makefile.conf \
  --rebuild --output-dir /evidence/srcsep/all-rebuilt

# Dependency-light focused suites are repeatable and use the Make build rules.
python3 test/run_tests.py --suite parker --suite fte-dmumu \
  --output-dir /evidence/srcsep/focused
python3 test/run_tests.py --suite controlled-analytical \
  --output-dir /evidence/srcsep/controlled

# Plot a previously archived component-test report without rerunning tests.
python3 test/run_tests.py --from-json /evidence/results.json \
  --output-dir /evidence/srcsep/replot --formats png,eps
```

`--validation-case ID` selects a registered end-to-end case from
`validation/case_registry.json`; `--validation-all` runs the complete currently
implemented portfolio. Both use `--amps PATH` (or `SEP_EXECUTABLE`) and refuse
to run if the linked binary is missing, rejects `--list-tests`, or does not
advertise the selected case. `--case-input PATH` can override one CV/IV case,
but is rejected for XM02/XM03 because those tests have one registered
publication-derived input each. To use sanitizers, build the linked application with the desired
sanitizer flags and pass that executable; the runner never compiles a substitute
driver. These cases retain richer native/model/reference artifacts while using
the same aggregate report and plotting contract as other runner modes.
See [../validation/cases/README.md](../validation/cases/README.md) for the
required structure and the CV01-CV12, IV01-IV06, and XM01-XM03 subdirectory READMEs for their equations,
inputs, gates, outputs, and failure interpretation.

For CV01 the Python runner executes six commands of the form `amps --test CV01
--test-input <generated-native-args> --test-output-dir <isolated-directory>`.
Each command also requests native JSON and JUnit and must return a registry
`PASS` plus a nonempty model CSV before the independent reference runs. The
generated native argument file is an internal, one-token-per-line protocol;
`resolved_input.json` remains the authoritative unit-bearing configuration.
CV02-CV12 use the same protocol once per case. Their C++ callback writes raw
model evidence only; Python then runs an independent analytical or finite-
volume reference and generates both case-specific and common PNG/EPS overlays.
The advanced cases additionally retain modal leakage, event/front statistics,
censored arrival histories, spectral fits, per-bin wave state, and energy ledgers.

### Application runner and build ownership

`srcSEP/test/run_tests.py` is the field-line application runner. Installing the
3-D application's identically named script here produces the reported
`No rule to make target 'test/stage1'` or, with `--rebuild`,
`No rule to make target 'clean-standalone'`. Those targets are intentionally
absent from srcSEP. Restore this application's runner; do not add no-op targets
or use `--no-build` to bypass the error in a mismatched harness.

The supported native command is:

```sh
test/run_tests.py --amps ../amps --all --output-dir test_output/all
```

Without `--rebuild`, it launches `../amps --list-tests` and then one native
process per discovered ID; it does not invoke `make`. With an explicit
`--rebuild`, it validates the enclosing AMPS root and `Makefile.conf`, calls
`make -C <AMPS-root> clean`, and invokes this application's `strict-production`
target with absolute `AMPS_ROOT`/`AMPS_CONFIG`. That target enters the enclosing
`make amps` build and audits the shared model archive before registry discovery.
Use this only when the configured AMPS checkout builds srcSEP; an executable
built for srcSEP3D is not interchangeable with the field-line application.

Source-only work uses `--suite python-runner`, `--suite cli` or the documented
focused Make targets. It has no monolithic `test/stage1` executable. The
subprocess regression in `test_python_test_runner.py` runs the public native
CLI with a literal protocol fixture and checks both build ownership and report
retention. This is an orchestration test, not evidence of native AMPS physics.

`--routine` forwards the native `--all-tests` policy and therefore excludes
extended cases. `--all` first calls `--list-tests`, accepts only identifier
tokens containing a digit (so the printed `ID | ...` heading can never become
a selection), and runs every discovered ID in a separate process. This includes
the extended `D01`, `D02`, and `D03PRE` registry entries. Ordinary
component tests use an explicit command of the form `amps --test PARK01
--test-json ... --test-junit ...`. Registered CV/IV/XM IDs instead use one
isolated `validation/run_case.py --case ID` command so their reviewed input,
native argument manifest, and independent reference are honored. That case
runner prints the exact nested `amps --test ID --test-input ...` command before
starting the linked application.

Isolation means an assertion failure, nonzero exit, timeout, abort, or
segmentation fault produces a `FAIL`/`ERROR` record for only that ID. If the
child could not write JSON, the Python runner creates an explicit synthetic
`ERROR` record, prints it immediately, and proceeds with the next test. At the
end it merges all per-test records into top-level JSON/JUnit and prints
`Overall test summary: TOTAL=... PASS=... FAIL=... SKIP=... ERROR=...`.
After the counts, `Failed tests (N)` and `Error tests (N)` list the affected
IDs and their diagnostics in stable ID order. Both headings are always shown;
a clean category contains `none`. This actionable block is deliberately
printed after the results-directory and plotting messages, so it remains at
the bottom of even a long terminal transcript.
The process exit remains 0 when all non-skipped tests pass, 1 when at least one
test fails, and 2 when at least one test errors. Consequently `--all` always
continues; `--keep-going` controls only repeated source-suite targets.
`--mpi-np N` and `--mpiexec PATH` launch ordinary native selections under MPI;
registered validation cases remain serial. Arguments after a literal `--` are
passed unchanged to direct AMPS invocations. `--timeout` applies per command.

The output directory contains `test-run.log`, `run_manifest.json`,
`srcsep-tests.json`/`srcsep-tests.xml`, `analytical_plot_manifest.json`, and a
`plots/` directory. Isolated `--all` runs additionally retain each native or
validation report beneath `individual/<ID>/`; their run manifest records the
ordered command list rather than a single monolithic command. The manifest
also hashes both the aggregate report and native executable. JSON is
authoritative; images are review aids.

### B04 build and orchestration result records

Every `--suite` selection now creates an outer result even when no C++ test
executable is produced. Its stable ID is `SUITE-<NORMALIZED-SELECTION>` and its
record contains `status`, monotonic elapsed seconds, UTC start/end timestamps,
the exact normalized Make command, raw return code, and a bounded diagnostic
excerpt. The complete interleaved output remains in `test-run.log`; the JSON
excerpt retains only the final 16 KiB so a compiler-error cascade cannot
exhaust runner memory.

The outer record and any child registry records copied through
`SRCSEP_REPORT_DIR` are merged before the terminal summary,
`srcsep-tests.json`, and `srcsep-tests.xml` are written. Thus a compilation
failure cannot disappear merely because the test executable never created its
own report. Reusing an output directory removes only stale focused JSON/JUnit
files before launch, preventing an earlier PASS from masking a current build
failure.

External status classification is explicit:

- return 0 is `PASS`;
- direct return 77 is an intentional `SKIP`;
- the exact output line `SRCSEP_SUITE_RESULT=SKIP` plus return 0 is also
  `SKIP`; this is the supported path through Make, which converts recipe 77
  into its own failure status;
- an ordinary nonzero Make/shell return is `FAIL`;
- launcher/permission errors, timeout, and signal termination are `ERROR`;
- without `--keep-going`, later selected suites receive explicit `SKIP`
  records naming the failed/error suite that blocked them.

Overall status is calculated only from merged records: any `ERROR` gives exit
2; otherwise any `FAIL` gives exit 1; otherwise `PASS` plus optional `SKIP`
gives exit 0. `run_manifest.json` stores that same code, and JUnit recomputes
its counts from the same records rather than trusting copied summary fields.
The Python self-test exercises a real invalid-C++ compile plus synthetic
success, timeout, missing executable, signal, skip, fail-fast blocking, and
cross-format count/identifier checks.

There are two deliberately distinct image types:

- **Solution-series comparison:** a declared artifact CSV has at least two
  finite rows and recognized coordinate (`x`, `time`, `s`, `mu`, and documented
  unit-bearing variants), numerical (`numerical`, `model`, `simulated`), and
  analytical (`analytic`, `analytical`, `exact`, `reference`) columns. The
  figure overlays those two series.
- **Metric-level comparison:** no eligible CSV exists, so reported numerical
  errors, moments, orders, probabilities, or conservation residuals are shown
  against their analytical reference/acceptance tolerance. The title labels
  this explicitly; it is not a fabricated pointwise solution.

Only controlled IDs explicitly classified in `ANALYTICAL_IDS`, or results that
declare a recognizable analytical-series CSV, are plotted by default.
`--plot all` also plots other tests with comparable finite metrics; `--plot
none` disables Matplotlib entirely. `--formats png`, `--formats eps`, or
`--formats png,eps` selects outputs. A failed test retains its figures and the
runner returns the originating status (`1` for FAIL, `2` for ERROR); plots can
never convert a failed registry result to PASS.

The focused Parker, Dmumu, MFP, and turbulence scripts normally delete all
temporary evidence. When invoked through this runner they honor the internal
`SRCSEP_REPORT_DIR` contract and copy only validated JSON/JUnit reports before
cleanup; sanitizer executables and objects are never retained. Source suites
that do not emit child reports still produce their B04 outer result in console,
JSON, and JUnit. They have no physical metrics to plot, but their build status
is part of campaign evidence. Verify orchestration and both figure paths with:

```sh
make test-python-runner-unit
```

## Focused WP21–WP30 tests

`test/run_wp21_wp30_tests.sh` compiles the production dependency-free kernels
and executes one controlled assertion per work package. Coverage includes exact
shock knots/restart partitioning, oriented multi-crossings, normalized spectra,
position-sensitive source keys, relative-normal shock-wave energy and branch
closure, conserved flux/refinement/restart, relativistic Larmor and invalid
exclusion, no-clamp bins plus physical products/uncertainty, transactional
checksum/duplicate/corruption behavior, and run-configuration
precedence/fingerprint/restart mismatch. A final source gate verifies each core
is reached from its production adapter and that sampling output contains no
shell invocation or unchecked `sprintf`.

The target is not a substitute for `test-native-amps-validation`,
`test-swmf-validation`, or `test-observational-validation`.

## Focused WP31–WP41 tests

`test/run_wp31_wp41_tests.sh` compiles the production dependency-light
contracts under C++11, strict warnings, ASan, and UBSan. It checks:

- `WP31`: a shared stiff reflection limit and automatically fitted first-order
  cascade refinement (`p=0.969891` in the delivered environment);
- `WP32`: typed correction/rejection accounting and nonfinite-state failure;
- `WP33`: exact split/merge moments plus stable lineage identifiers;
- `WP34`: refusal to promote a source double to native evidence and validation
  of the complete native-observation record;
- `WP35`: all 90 mover/source/ownership/coupling rows classified, preflighted,
  and rendered from one registry;
- `WP36`: global number/charge/energy/momentum closure and checksummed restart;
- `WP37`: versioned domain-separated seed panels and explicit mean statistics;
- `WP38`: IEEE boundary generation, reproducible counterexamples, and fault hits;
- `WP39`: analytical response-folded instrument counts and uncertainty;
- `WP40`: exact normalized-work and environment-specific timing decisions;
- `WP41`: claim/evidence validation and exactly one queue-flush owner.

WP34 native mover execution, WP36 full PIC seam activation, WP37 scheduled
ensembles, WP38 native adapter fuzzing, WP39 held-out events, and WP40 MPI/OpenMP
scaling require the enclosing application or external evidence. Their absence
is BLOCKED, never converted into a source-only PASS.

## B05 WP42–WP64 disposition and focused tests

`WP42_WP64_DISPOSITION.json` is the authoritative inventory for all 23 work
packages. The checker requires the exact ordered WP42–WP64 set, one approved
status and reason for every entry, a named implementation symbol or evidence
gate, and a concrete promotion condition. It also proves that the four
experimental extension implementations are absent from `MAINLIBOBJ`, that the
accepted supporting contracts are reached by active production code, and that
both runner selections remain exposed. Its negative control mutates a copy of
the manifest and requires the audit to reject that copy.

`test/run_wp42_wp64_tests.sh` compiles the experimental component surface with
C++11, strict warnings, ASan, and UBSan, then executes one stable WP42–WP64
check per identifier. The checks cover persistent turbulence/checkpoint
concepts, transaction and rollback semantics, geometry/measure contracts,
adaptive integration and convergence, validation/evidence records, and
performance accounting. This is deliberately an algorithm/component suite:
`sep_runtime_contracts.cpp`, `sep_physics_extensions.cpp`,
`sep_numerical_extensions.cpp`, and `sep_validation_extensions.cpp` are not
production archive members and have no production selector.

Two narrow supporting changes are production-active and are audited separately:

- canonical `sep_transport_common.h` defines `ParkerMeasure`, making the
  distribution measure explicit instead of relying on an undocumented factor;
- the event-driven MFP core records `resonantBranch` and the actual
  `preWaveMomentumKgMPerS`/`postWaveMomentumKgMPerS` at a completed scattering
  event, and `focused_transport_mfp.cpp` deposits those values only after the
  event/no-event transaction succeeds.

The native/external target covers WP59, WP60, WP63, and WP64. With no reviewed
site command, `run_wp59_wp64_native_gates.sh` emits
`SRCSEP_SUITE_RESULT=SKIP`; B04 therefore records
`SUITE-WP59-WP64-NATIVE` as `SKIP` in terminal, JSON, and JUnit output. A site
may run an approved command explicitly:

```sh
make test-wp59-wp64-native \
  SRCSEP_NATIVE_GATE='/absolute/path/to/reviewed-native-gate --options'
```

That command is responsible for retaining its native AMPS/SWMF, external-data,
or scaling evidence. A successful external gate does not by itself insert the
experimental objects into `mainlib.a`; production promotion additionally
requires a selector/call path, restart compatibility, native observations, and
an explicit reviewed update to the disposition and source manifest.

## Focused WP11--WP20 tests

`test/run_wp11_wp20_tests.sh` compiles the exact dependency-light production
kernels and executes one controlled case per work package:

- `WP11`: stationary isotropic bounded Milstein diffusion, endpoint symmetry,
  and proof that no coefficient call leaves `|mu|<=1`;
- `WP12`: tolerance validation, full-step/two-half-step error estimate,
  accepted-step accounting, and named limiter histogram;
- `WP13`: authoritative source-field perturbation and irrelevant-field
  isolation through the pure coefficient view;
- `WP14`: constant Dmumu value/derivative and invalid-value rejection;
- `WP15`: Jokipii analytic derivative versus finite differences and finite
  endpoint limits;
- `WP16`: Florinskiy branch mirror symmetry, output assignment, and bounded
  derivative behavior;
- `WP17`: adaptive recovery of `kappa=v^2/(6D0)` and explicit resonance-gap
  reject/ballistic states;
- `WP18`: proton/alpha/electron scaling, exact species-abundance closure,
  energy-per-nucleon conversion, and incomplete-definition rejection;
- `WP19`: validation of named turbulence/spectrum scales and removal of the
  sampling hard-coded field;
- `WP20`: typed ballistic lambda, event-mover compatibility, and Parker
  preflight rejection.

The runner adds source gates for production adapter delegation and uses strict
warnings plus ASan/UBSan. It does not link the PIC adapter; native AMPS, real
SWMF, and observational evidence remain separate validation gates.

## Focused Step 15 scientific-validation tests

`test/run_step15_tests.sh` runs four source-distribution cases and assembles an
auditable campaign report:

- `VAL01`: analytical Parker displacement moments and adiabatic cooling;
- `VAL02`: matched `fte-dmumu`/`fte-mfp` pitch moments, diffusion, detector
  onset, peak time, and fluence;
- `VAL03`: production Dmumu particles versus an independently coded
  conservative finite-volume pitch-angle solver;
- `VAL04`: actual `swcme::sep::Interface1D` background/source records consumed
  by immutable srcSEP snapshots and the production Parker core.

The runner also verifies that missing real SWMF and held-out spacecraft inputs
remain `INCOMPLETE`, that `--release` returns nonzero, and that manufactured
labels cannot acquire observational status. Set `STEP15_REPORT_DIR` to an
absolute directory to retain JSON, JUnit, and campaign reports; by default all
outputs are disposable. Full formats and external-manifest requirements are in
[../validation/README.md](../validation/README.md).

## Focused Step 13 acceptance tests

`test/run_step13_tests.sh` compiles the exact descriptors linked into the
production registry and verifies:

- `BG01`: standalone analytic and SWCME snapshots preserve distinct provider,
  epoch, ownership, validity, generation, and configuration identity;
- `BG02`: a mock SWMF snapshot is imported read-only and can become locally
  evolved only through an explicit handoff copy;
- `CROSS01`: `fte-dmumu` at `Dmumu=0` and `fte-mfp` at
  `lambda_parallel=+infinity` produce the same ballistic state;
- `CROSS02`: Parker `kappa`, focused-transport `Dmumu`, and event-driven
  `lambda_parallel` round-trip under their declared isotropic closure;
- `HIDDEN`: a callback cannot return PASS with a retained assertion failure;
- `REPORT`: JSON and JUnit retain every result's diagnostic evidence.

The fixtures are deterministic and small enough for routine acceptance. Large
Monte Carlo checks remain in the explicitly selected Parker/FTE stress cases.

## Focused Step 14 cleanup tests

`test/run_step14_tests.sh` verifies:

- `DOC01`: every documented focused command and structured-report option maps
  to a real Make target or parser option;
- `DOC02`: public headers, dispatch, and build inputs contain only the three
  canonical movers and no retired mutable function-pointer API;
- `DOC03`: the source handoff contains no objects, libraries, coverage data,
  generated test evidence, or nested archive;
- `WARN-FULL-SOURCE`: the root-only sampling writer uses structured control
  flow rather than a jump across initialized C++ objects, and field-line list
  sizes use a `size_t`-compatible output format;
- `WARN01`: every dependency-light Step 13/14 translation unit compiles with
  C++11, `-Wall -Wextra -Wpedantic -Werror`;
- `STATIC01`: GCC's `-fanalyzer` checks the new infrastructure when supported;
- `SAN01`: the exact Step 13 production callbacks pass ASan and UBSan.

The full native `SAN01`/MPI acceptance command is `make test` in an enclosing
AMPS checkout. A source archive cannot emulate generated PIC types or MPI rank
execution and reports that boundary explicitly.

## Focused Step 6 common-transport tests

`test/run_step6_tests.sh` verifies:

- `CORE01`: relativistic SI speed/momentum round trips and rejection of
  luminal input;
- `CORE02`: a physical arc-length advance does not leak host coordinate units;
- `CORE03`: absorbing and reflecting boundary policies are explicit;
- `CORE04`: focusing length is exactly `-1/(d ln|B|/ds)`, including the uniform
  field limit;
- `CORE05`: exact plasma-frame adiabatic momentum evolution;
- `CORE06`: named timestep limits, limiter diagnostics, and explicit underflow;
- `CORE07`: keyed stochastic streams reproduce independently of particle order;
- `CORE-SOURCE`: all three production shells use common loading, advancement,
  commit, and field-line attachment code.

## Focused Step 7 Parker tests

`test/run_step7_tests.sh` verifies:

- `PARK01`: zero-diffusion convection;
- `PARK02`: Gaussian displacement mean and variance `2*kappa*dt` over 120,000
  deterministically keyed samples;
- `PARK03`: the Itô variable-diffusion drift has the `+d(kappa)/ds` sign;
- `PARK04`: the manufactured adiabatic momentum solution;
- `PARK05`: explicit absorbing boundary crossing;
- `PARK06`: refinement convergence for a variable-kappa manufactured drift;
- `PARK07`: absorbing first-passage probability and mean exit time for
  drift-free Brownian transport on a finite interval;
- `PARK-SOURCE`: the public registry selects the canonical provider-driven
  implementation and its source contains no hidden/global random draw.

## Focused Step 8 coefficient-driven FTE tests

`test/run_step8_tests.sh` verifies:

- `FTED01`: constant-`Dmumu` Itô ensemble moments;
- `FTED02`: inclusion and sign of `dDmumu/dmu` drift;
- `FTED03`: the shared magnetic-focusing convention;
- `FTED04`: combined focusing, velocity gradient, streaming, and cooling;
- `FTED05`: regular `mu=0` behavior and arbitrary reflective overshoot;
- `FTED06`: QLT normalization/units plus matching coefficient and turbulence
  identity for wave deposition;
- `FTED07`: refinement convergence of the symmetric split;
- `FTED08`: decay of a Legendre `P2` eigenmode according to
  `<P2>(t)=(a/5) exp(-6 D0 t)`;
- `FTED-SOURCE`: canonical dispatch, coefficient-provider isolation, explicit
  randomness, and deferred wave coupling.

These stable IDs are emitted by the focused test binary. They complement the
linked standalone CLI catalog because they are designed to run from a
source-only handoff; native adapter and coupled-execution evidence remains a
separate gate.

## Focused Step 9 event-driven FTE tests

`test/run_step9_tests.sh` verifies:

- `FTEM01`: exponential waiting-time mean `1/nu`;
- `FTEM02`: Poisson event-count mean `nu*dt`;
- `FTEM03`: infinite-lambda ballistic transport;
- `FTEM04`: wave-frame speed conservation under Lorentz scattering;
- `FTEM05`: magnetic-focusing refinement;
- `FTEM06`: exact adiabatic cooling between events;
- `FTEM07`: combined focus/cooling/streaming event-split refinement;
- `FTEM08`: finite-time persistent-flight mean-square displacement and its
  long-time `kappa_parallel=v*lambda_parallel/3` diffusion limit;
- `FTEM-SOURCE`: canonical dispatch and removal of the legacy MFP object from
  the production archive.

## Focused Step 10 coefficient tests

`test/run_step10_tests.sh` verifies:

- `COEF01`: canonical registry names, SI units, and parameter schemas;
- `COEF02`: invalid parameter/source combinations and conversion cycles;
- `COEF03`: lambda/kappa and isotropic-Dmumu conversion round trips;
- `COEF04`: analytic/imported source identity with identical SI conversion;
- `COEF05`: structured invalid-domain and ownership status;
- `COEF06`: parsed constant Dmumu reaches the same pure provider value;
- `COEF-SOURCE`: all three canonical movers and CLI use the shared registry.

## Focused Step 11 turbulence tests

`test/run_step11_tests.sh` executes registered `TURB02`–`TURB23` plus
`TURBOWN01`: source ownership,
particle-wave closure, units/volume, initialization, integrated/spectral
projection, constant/variable-speed advection, boundary policies, reflection,
cascade/dissipation, growth, resonance/coefficient closure, shock injection,
operator order, positivity limiting, conservative remap, restart,
reproducibility, mover/source separation, and standalone/coupled driver
equivalence. `TURB21` compares nonuniform advection with an exact translated
sine profile at two resolutions; `TURB22` compares every time level of a
linear-in-time wave-energy source with its quadratic integral. `TURB22`
isolates the source-application operator: calculating the source from a
particle distribution remains covered by the native
`growth_rate_validation_test.cpp` path, which requires the full AMPS/PIC/MPI
environment. The focused runner also validates each Step 11 CLI family and
source contract. `TURB23` converts kinetic-energy changes from controlled
outward/minus-branch and inward/plus-branch wave-frame scattering events into
opposite signed wave increments and checks total particle-plus-wave energy and
ledger closure without invoking an external turbulence data source.

## Analytical refinement evidence

`PARK06`, `FTED07`, `FTEM05`, `FTEM07`, and `TURB21` call the common
`EstimateRefinementOrder` helper. Given positive errors `e_c`, `e_f` and
resolutions `h_c>h_f`, it records
`p=log(e_c/e_f)/log(h_c/h_f)` plus both input ratios. Invalid, zero-error, or
reversed-resolution inputs fail instead of producing a manufactured order.
The Step 1 `REFINE` fixture independently checks an exact second-order example
and the invalid-input paths.

## Separate validation gates

Use `make test-controlled-analytical` for dependency-light equation-level
tests. Native AMPS, real SWMF, and held-out spacecraft evidence are separate:

```sh
make test-native-amps-validation SEP_EXECUTABLE=/path/to/amps
make test-swmf-validation SWMF_MANIFEST=/evidence/swmf-replay.json
make test-observational-validation \
  OBSERVATIONAL_MANIFESTS="/evidence/event-1.json /evidence/event-2.json"
```

The latter two commands reopen and hash every manifest input. The
observational command also requires complete onset, anisotropy, spectra,
fluence, decay, and multi-spacecraft-longitude coverage. Missing arguments are
errors, not `SKIP`, and no gate inherits PASS from another evidence class.

## Focused Step 12 reproducibility tests

`test/run_step12_tests.sh` executes:

- `PAR01`: worker-local accumulation followed by a single canonical writer;
- `PAR02`: identical evidence across worker/scheduler layouts;
- `PAR03`: identical evidence across synthetic MPI partitions;
- `PAR04`: mover RNG independence from a diagnostic purpose stream;
- `PAR05`: stable evidence hashes and resettable atomic integer counters.

The source gate confirms that rank/thread IDs are absent from physical keys and
that MPI policy is gather followed by the same canonical reduction.

## Focused Step 3 geometry and source tests

`test/run_step3_tests.sh` builds the production
`util/sep_flux_tube_geometry_core.cpp` with C++11, strict warnings,
AddressSanitizer, and UndefinedBehaviorSanitizer. It requires neither AMPS nor
MPI and verifies:

- `GEOA01`: magnetic-flux conservation, `A|B| = constant`;
- `GEOA02`: fourth-order convergence of segment-volume integration;
- `SRC01`: swept-volume dimensions and numerical value;
- `SRC02`: equal injected physical weight for provider-equivalent shock states;
- `SRC03`: spectral normalization and the MeV-to-joule API boundary.

The PIC-facing adapter is exercised by the native regression gate because its
vertex and segment types are supplied by the enclosing AMPS checkout. See
[../FLUX_TUBE_GEOMETRY.md](../FLUX_TUBE_GEOMETRY.md).

## Focused Step 4 production-mover tests

`test/run_step4_tests.sh` compiles the exact production registry and CLI parser
with C++11 and strict warnings. It verifies:

- `MOVCLI01`: canonical `parker`, `fte-dmumu`, and `fte-mfp` parsing and a
  three-entry help/discovery surface;
- `MOVCLI02`: hard rejection of every retired transition alias;
- `MOVCLI03`: coefficient/state capability reporting and capability-based
  main-loop policy;
- `MOVCLI04`: rejection of ambiguous, direct-wave, legacy, and 3-D movers.

The linked native gate remains responsible for passing a real particle through
each PIC adapter mapping. See
[../PRODUCTION_MOVER_API.md](../PRODUCTION_MOVER_API.md).

## Focused Step 5 field-line-scope tests

`test/run_step5_tests.sh` verifies the source and build boundary without PIC or
MPI:

- `SCOPE01`: transferred Parker3D, 2019 He, Kartavykh, Borovikov, drift, Boris,
  and spatial-neighborhood sampling symbols and object files are absent from
  production sources;
- `SCOPE02`: production mover sources cannot write a Cartesian particle
  position, call a Boris pusher, apply cross-field diffusion, or attach to an
  AMR-node particle list; the adapter enforces segment attachment;
- `SCOPE03`: Cartesian field-line embedding, vector magnetic-field access,
  SI flux-tube geometry, and field-line observer sampling remain present.

The script then runs all `MOVCLI` tests, ensuring the removal leaves exactly the
three canonical public movers. Native srcSEP and receiving-application builds
remain required on a complete AMPS checkout. See
[../STEP5_CHANGE_MANIFEST.md](../STEP5_CHANGE_MANIFEST.md).

## Focused Step 2 state and clock tests

`test/run_step2_tests.sh` builds the exact production
`util/sep_background_snapshot.cpp` implementation in a disposable directory
with C++11, `-Wall -Wextra -Werror`, and pthread support. It verifies:

- `BKG01`: required metadata, validity domains, read-only SWMF ownership, and
  compile-time snapshot non-assignability;
- `BKG02`: one snapshot per particle read phase, mover acquisition, stale-time
  rejection, and no publication while particles are active;
- `BKG03`: cross-provider overwrite rejection and an explicit, provenance-bearing
  SWMF-to-local handoff with a new field-line generation;
- `BKG04`: concurrent scheduler threads acquire the identical const generation;
- `BKG05`: deterministic configuration fingerprints change when background
  configuration changes;
- `CLK01`: production code contains no standalone elapsed-time/launch clock and
  reads `PIC::SimulationTime::Get()` only through the runtime clock adapter.

These checks need neither AMPS nor MPI. The linked native executable is still
required to prove the snapshot boundary around the real `PIC::TimeStep()` in a
standalone and SWMF-coupled run. See
[../BACKGROUND_STATE.md](../BACKGROUND_STATE.md) for the full runtime contract.

## Stage 1–2 integration-boundary gates

`make test-stage1-stage2-contracts` runs
`test/check_stage1_stage2_contracts.py` without AMPS or MPI. The audit checks
properties that equation-level tests cannot infer from numerical output:

- `main_lib.cpp` contains exactly one live turbulence `Advance` call and
  standalone `main.cpp` contains none; the call must belong to common
  `amps_time_step()`, and four negative controls reject duplicate, missing or
  foreign owners;
- `main.cpp` does not include the retired application-local coefficient
  registry after `sep.h` has loaded the canonical `sep_common` registry;
- the legacy unconditional initializer is not executable and local wave
  initialization names only the two sources owned from startup;
- analytic density and SWCME preparation are guarded by the configured
  provider;
- `FTE01`, `PARKER01`, and `PARKER02` acquire `BeginParticleRead` before their
  first production mover call; and
- the CV01–CV12 and IV01–IV06 detached compiles include the canonical
  repository root.

`make test-controlled-analytical` runs this audit before its numerical suites.
`make test-cv01-unit` then strictly compiles the native CV01 model in a
disposable directory. If `SEP_EXECUTABLE` is set to a linked application, the
same target also runs the complete linked comparison; otherwise that portion
is explicitly `SKIP`.

## Stage 3 release, validation, and species gates

`make test-stage3-contracts` runs `test/check_stage3_contracts.py`. This is a
source-structure test for invariants that cannot be inferred from a successful
numerical trajectory. It requires all four `SOURCE_MANIFEST.json` files, proves
that `srcSEP/mover.cpp` is absent and retired, checks the unified runner/Make
and native OV/EV registration, verifies each evidence role and fixed-input
policy, and follows the `srcSEP3D` species index from configuration through
observer selection and injection.

The case-specific targets provide complementary executable evidence:

- `test-ov01-ov05-unit` compiles the production validation callback, verifies
  immutable publication inputs/references and registry metadata, and runs the
  linked cases when `SEP_EXECUTABLE` is available;
- `test-ev01-ev02-unit` verifies production dispatch, the nine-event CCMC
  source/observation provenance, strictly increasing source histories,
  plotting attribution, and campaign-runner syntax; and
- `test-package-hygiene-unit` executes the four-component allowlist audit plus
  positive and negative self-tests.

Evidence labels remain binding. `OV01`/`OV02` are release gates,
`OV03`–`OV05` are diagnostic only, and `EV01`/`EV02` are incomplete pilots.
A source-only linked-case `SKIP` is not converted into a release pass.

## Registered tests

| ID | Group | Class | Initialization | What is asserted | State/artifacts |
|---|---|---|---|---|---|
| `CV01` | `controlled-analytical` | routine | linked srcSEP/AMPS registry | A localized 10 MeV proton packet follows exact ballistic characteristics for four pitch angles, two boundary policies, and three timesteps; packet moments, crossing times, momentum, and active/escaped weight close against an independent solver. | Linked executable hash; case-registry input; seed 10101; native JSON/JUnit and raw model/reference CSV; aggregate JSON/JUnit; PNG/EPS overlay, residual, and four-panel evidence. |
| `CV02` | `controlled-analytical` | routine | linked srcSEP/AMPS registry | Constant spatial diffusion matches Gaussian moments, exact bin probabilities, and fitted kappa across 10 seeds, three particle counts, and three timesteps. | Seed 20202; raw ensemble moments/profiles; exact-bin reference; negative-control ratio; JSON/JUnit, hashes, PNG/EPS. |
| `CV03` | `controlled-analytical` | routine | linked srcSEP/AMPS registry | Sinusoidal diffusion preserves uniform equilibrium and matches an independent conservative finite-volume transient with analytic/numerical derivative paths. | Seed 30303; four seeds; reversed-drift control; refinement/closure metrics; JSON/JUnit, hashes, PNG/EPS. |
| `CV04` | `controlled-analytical` | routine | linked srcSEP/AMPS registry | Constant-divergence and spherical-wind adiabatic momentum changes match exact proton/alpha relativistic characteristics and designed order. | Seed 40404; three energies/species and timesteps; momentum/energy histories; JSON/JUnit, hashes, PNG/EPS. |
| `CV05` | `controlled-analytical` | routine | linked srcSEP/AMPS registry | Zero-scattering magnetic focusing matches exact position/pitch characteristics for both gradient signs while preserving bounds, momentum, angular moments, and invariant. | Seed 50505; endpoint/near-endpoint pitch angles; three timesteps; JSON/JUnit, hashes, PNG/EPS. |
| `CV06` | `controlled-analytical` | extended | linked srcSEP/AMPS registry | Production pitch diffusion evolves Legendre modes 1–6 with the exact eigenvalue decay and bounded leakage. | Seed 60606; three seeds/timesteps; raw modal matrix and boundary counts; JSON/JUnit, hashes, PNG/EPS. |
| `CV07` | `controlled-analytical` | extended | linked srcSEP/AMPS registry | Persistent random flights reproduce telegraph fronts, causal support, MSD, events, and late diffusion. | Seed 70707; three rates, four regimes, ten seeds; profiles/moments; JSON/JUnit, hashes, PNG/EPS. |
| `CV08` | `controlled-analytical` | extended | linked srcSEP/AMPS registry | Absorbing Parker trajectories match the inverse-Gaussian first-passage CDF and overshoot refinement. | Seed 80808; two drifts, three timesteps, ten seeds; censored histories and negative control. |
| `CV09` | `controlled-analytical` | extended | linked srcSEP/AMPS registry | Controlled shock cycles reproduce planar DSA indices and acceleration times for r=2,3,4. | Seed 90909; 50k particles/ratio; spectrum, timing, accounting, diffusion-length evidence. |
| `CV10` | `controlled-analytical` | routine | linked srcSEP/AMPS registry | Both spectral wave branches advect conservatively across fixed, expanding-area, and remapped grids. | Three resolutions; per-bin profiles, invariant/non-negativity/order metrics; PNG/EPS. |
| `CV11` | `controlled-analytical` | routine | linked srcSEP/AMPS registry | One-hot wave energy matches constant and time-dependent exponential growth/damping histories. | Three timesteps; active/inactive bins, cancellation, positivity, second-order rate integration. |
| `CV12` | `controlled-analytical` | routine | linked srcSEP/AMPS registry | Coupled wave-frame scattering closes total energy while uncoupled controls remain fixed and suppressed deposition fails. | Three counts/timesteps; branch exchange and per-step ledger residuals; PNG/EPS. |
| `IV01` | `integrated-manufactured` | routine | linked srcSEP/AMPS registry | Parker-spiral focusing and flight time agree with independent characteristics and are invariant to vertex order. | Two wind speeds; ballistic/weak scattering; three timesteps; per-particle CSV and PNG/EPS. |
| `IV02` | `integrated-manufactured` | routine | linked srcSEP/AMPS registry | Coupled scattering/focusing converges from isotropic and beam states to the zero-flux exponential PDF. | Three focusing ratios; 40-bin PDFs; normalization and first-moment metrics. |
| `IV03` | `integrated-manufactured` | routine | linked srcSEP/AMPS registry | The full smooth manufactured residual converges at second order without loss of positivity. | Exact/operator residual rows at three levels; L2/order evidence. |
| `IV04` | `integrated-manufactured` | routine | linked srcSEP/AMPS registry | Uniform physical state survives conservative remap on translated-shape grids at roundoff. | Four node motions, three resolutions, node/segment/integral/ledger diagnostics. |
| `IV05` | `integrated-manufactured` | routine | linked srcSEP/AMPS registry | Moving and stationary shock frames give the same crossings, momentum, and DSA spectrum. | Two shock speeds, three timesteps, 5k histories, exact-node flag. |
| `IV06` | `integrated-manufactured` | routine | linked srcSEP/AMPS registry | Resonant wave growth strengthens scattering and self-limits streaming while total energy closes. | Frozen/one-way/two-way timelines, resonant bin, Dmumu, streaming and ledger evidence. |
| `XM01` | `cross-model` | extended | linked srcSEP/AMPS registry | Production focused-transport samples agree with an independent conservative PDE solver across isolated and combined operators. | Two sampling refinements; full `(s,mu)` probability, intensity, anisotropy, momentum, JSON/JUnit, and PNG/EPS. |
| `XM02` | `cross-model` | extended | linked srcSEP/AMPS registry | A controlled production-core first-passage ensemble uses the reported 0.05/0.3/1.0 au MFPs and compares unit-peak profiles with Zhao et al. Figure 7. | Fixed publication input, native model CSV, digitized reference/provenance, metrics, and PNG/EPS; no external CSV. |
| `XM03` | `cross-model` | extended | linked srcSEP/AMPS registry | Event-informed one-field-line Parker transport is compared with ACE/EPAM, GOES-13/EPEAD, and SOHO/ERNE Earth spectra in Liu et al. Figure 12. | Fixed paper-derived input/source trace, 80 vector-extracted observations, one global amplitude, metrics, and individual 4/12/36 h publication PNG/EPS figures with in-panel legends; no external CSV. |
| `BG01` | `background` | routine | none | Standalone analytic and SWCME snapshots preserve provider, epoch, ownership, validity, generation, and distinct configuration identity. | Stack-owned immutable snapshots; no external provider or artifact. |
| `BG02` | `background` | routine | none | A mock SWMF import is read-only and becomes locally evolved only through an explicit handoff copy. | Resets the snapshot store before/after; no external SWMF process. |
| `CROSS01` | `cross-mover` | routine | none | `fte-dmumu` and `fte-mfp` agree in the matched ballistic limit. | Keyed seed 1301; stack-owned state; no artifact. |
| `CROSS02` | `cross-mover` | routine | none | Parker `kappa`, FTE `Dmumu`, and event-driven `lambda` round-trip under the declared isotropic closure. | Pure SI conversion; no RNG or artifact. |
| `D01` | `swcme-background` | extended | none | The production SWCME adapter rejects invalid/unprepared state, preserves transactional preparation, and exposes only explicit, counted clamp/fallback recovery. | Deterministic, no RNG; restores fast/strict state; identical callback used by `test-swcme-fail-closed-unit`. |
| `D02` | `swcme-configuration` | extended | none | The canonical resolver enforces presets, units, layer precedence, provenance, validity, supported source modes, and order-independent fingerprints. | Deterministic, no RNG; restores fast/strict state; identical callback used by `test-swcme-configuration-unit`. |
| `D03PRE` | `native-integration` | extended | none | The linked catalog contains the three field-line movers with the expected coefficient contracts, and SWCME refresh exposes a monotonic state ID/epoch plus fingerprint. | Preflight only; full D03 requires `--suite d03-native-integration` with reviewed site inputs. |
| `DXX01` | `diffusion` | routine | field-line model | `GetDxx` agrees with the constant-coefficient analytical result and the existing million-panel independent quadrature at relative tolerance `1e-5`. | Temporarily replaces and restores the pitch-angle diffusion function pointer; no artifact. |
| `FTE01` | `transport` | routine | field-line model | The focused-transport mover follows the expected field-line displacement while preserving velocity in a static-plasma fixture, using the existing `1e-2`/`1e-5` checks. | Uses fixed registry seed 1002; opens one immutable background read phase around all mover calls; then restores vertex data and diffusion pointer, deletes its particle, and clears test lists. |
| `PARKER01` | `parker` | routine | field-line model | Parker convection keeps the line coordinate stationary in the fixture and matches the analytical density-driven momentum update within `1e-5`. | Uses fixed registry seed 1001; opens one immutable background read phase around all mover calls; then restores vertex data and diffusion pointer, deletes its particle, and clears lists. |
| `PARKER02` | `parker` | extended | field-line model | Four million legacy stochastic trials produce at least one in-range displacement sample and rank 0 successfully writes the histogram.  This is an execution/output assertion, not yet a Gaussian-shape validation. | Uses fixed registry seed 1003; holds one immutable generation for the Monte Carlo loop; rank 0 writes `dxParker.dat`; particle/lists are cleaned. |
| `SCAT01` | `scattering` | extended | field-line model | Currently reports `SKIP`: the legacy return-probability diagnostic has no approved reference/tolerance and contains a singular zero-energy case. | Registry path performs no mutation. Historical TestManager may write `rmax-E=...` and `time-E=...` files. |
| `TURB01` | `turbulence` | routine | none | The production 1-AU helper equals the independently evaluated magnetic-pressure closure `delta_B^2/(2 mu_0)` within 32 machine epsilons. | Pure deterministic calculation; no RNG, model state, or artifact. |
| `PARK01`–`PARK07` | `parker` | routine except `PARK02`/`PARK07` extended | none | Controlled convection, Gaussian diffusion, Itô drift, cooling, boundary, refinement, and first-passage references. | Injected providers; fixed seeds where stochastic. |
| `FTED01`–`FTED08` | `fte-dmumu` | routine except `FTED01`/`FTED08` extended | none | Controlled Itô, focusing, cooling, reflection, QLT, refinement, and Legendre-mode references. | Injected providers; fixed seeds where stochastic. |
| `FTEM01`–`FTEM08` | `fte-mfp` | routine except `FTEM01`/`FTEM02`/`FTEM08` extended | none | Waiting-time, Poisson, ballistic, wave-frame, cooling, refinement, persistent-flight, and diffusion-limit references. | Injected providers; fixed seeds where stochastic. |
| `TURB02`–`TURB23`, `TURBOWN01` | `turbulence` | routine | none | Controlled ownership, ledger, operator, restart, translated-profile, time-dependent-growth, and particle-wave total-energy checks. | Stack-owned state; fixed seed only for `TURB23`; JSON/JUnit in disposable focused builds. |
| `VAL01` | `validation` | routine | none | Parker displacement moments and cooling match independent analytical solutions. | Fixed seed 1501001; report only when requested. |
| `VAL02` | `validation` | extended | none | Matched Dmumu/MFP ensembles agree in pitch, diffusion, onset, peak, and fluence. | Fixed seed 1502001; report only when requested. |
| `VAL03` | `validation` | extended | none | Dmumu ensemble agrees with an independent finite-volume solver. | Fixed seed 1503001; report only when requested. |
| `VAL04` | `validation` | routine | none | Actual SWCME SI background/source records drive srcSEP Parker transport. | Fixed API seed 1504001; focused C++17 registry; report only when requested. |

The list printed by native `--list-tests` is authoritative for C++ component
callbacks and also includes supported build modes, seed policy, and
state/isolation notes. End-to-end portfolio cases CV01-CV12, IV01-IV06, and XM01-XM03 are listed by
`python3 validation/run_case.py --list`; both paths are selected through the
common `test/run_tests.py` orchestration interface. Entries are sorted by ID
regardless of construction or registration order.

## WP01--WP10 focused contract gate

Run `./test/run_wp01_wp10_tests.sh` from the srcSEP directory. The disposable
C++11 strict-warning executable checks physical snapshot limits and generations,
physical background epochs, full gyrotropic coefficients, location-aware
midpoint sampling, carried nonhomogeneous optical depth, zero-probability empty
branches, and duplicate physical-key rejection. The same script reuses the
Stage 1–2 ownership audit: the common `amps_time_step()` owns exactly one
turbulence PIC-adapter call and the standalone driver owns none. Four in-memory
negative controls reject duplicate driver/library calls, a missing call and a
foreign owner. Requiring both entry points to advance waves independently would
consume the fixed physical interval twice. All three movers must still publish
self-contained coupling records without a stored particle pointer.

This is a source-only controlled gate. It does not claim a native AMPS build,
MPI decomposition, SWMF coupling, restart continuity, or observational PASS;
those remain separate targets and must run in their configured environments.

## Results, exit codes, and MPI

Each test returns one of `PASS`, `FAIL`, `SKIP`, or `ERROR`, plus a message,
duration, optional seed, metrics/tolerances, and artifact paths.  Exit status is:

- `0`: every requested result is `PASS` or `SKIP`;
- `1`: one or more scientific/test assertions are `FAIL`;
- `2`: one or more tests encounter `ERROR`.

In MPI execution every required rank calls the adapter. Status severity and
duration are reduced deterministically; the most severe rank status and maximum
duration are reported. Before reduction, root gathers each rank's complete
status/message/seed/configuration/metric/artifact evidence for JSON and JUnit.
Only rank 0 prints the registry summary. Shared output
from `PARKER02` is also root-only.  Legacy test bodies may still print detailed
per-rank diagnostics; their removal requires a later test refactor.

Selectors are case-insensitive and de-duplicated.  Unknown IDs/groups, missing
values, `--list-tests` combined with execution, or `--all-tests` combined with
explicit selectors fail before model initialization.  `--help` and
`--list-tests` are pre-initialization success paths.  A new execution selector
always exits before the production loop.

## Focused Step 1 contract tests

`test/run_step1_tests.sh` builds in a temporary directory using `-Wall -Wextra
-Werror`.  It compiles the exact production CLI source with
`SEP_CLI_PARSE_ONLY`, which excludes only the function that writes AMPS globals,
and exercises:

- `CLI01`: help/list mode selection, and acceptance of `-input FILE` as an
  exact synonym of `--input FILE`/`--input=FILE`;
- `CLI02`: both `--test` forms and exactly-once execution;
- `CLI03`: deterministic ordering, de-duplication, groups, routine/extended;
- `CLI04`: malformed/conflicting/unknown selection rejection;
- `CLI05`: result exit propagation and frozen legacy TestManager parsing;
- `CLI06`: mover precedence (`--particle-mover` > schema-6
  `application.sep.transport` > default) and deferred mover-dependent
  validation when an input file (`-input`/`--input`) is given without
  `--particle-mover`;
- `HIDDEN`: positive `assertion_failures` cannot be masked by callback PASS;
- `REPORT`: JSON/JUnit schema, outcome, metric, and artifact preservation;
- registry completeness: required metadata, callbacks, and unique IDs.

`HIDDEN01` is an intentional inner negative control: its fixture falsely marks
itself `PASS` while reporting one failed assertion. The registry must print
that inner result as `FAIL`; the outer `HIDDEN PASS` then confirms that the
contradiction was detected. The script's final summary labels this expected
internal failure explicitly and reports only unexpected outer failures in its
`FAIL` count.

The linked-host targets are still required for final evidence that early exits
precede real AMPS initialization and that test-only execution never enters the
production loop.  A source-only archive cannot substitute static scans for that
runtime evidence.

## Adding a component test

1. Write an adapter returning `SEP::Testing::Result`; never call `exit()` or
   infer PASS from printed output.
2. Add a descriptor in `ComponentTestRegistry()` with a unique stable ID,
   nonempty name/group/description, exact initialization level, supported build
   modes, runtime class, seed policy, and state-isolation contract.
3. Return metrics with comparison operators, tolerances, and units.  List every
   created artifact and make its shared writer MPI-root-only.
4. Restore production globals, function pointers, field-line data, particle
   lists, output settings, and deterministic random state with RAII where
   available.  Otherwise isolate the test and document the limitation.
5. Add focused selection/status tests, run `make test-cli-unit`, then run the
   linked individual/group target and the complete bounded `make -j test` suite.

Self-consistent Alfvén turbulence remains a production provider/evolution
subsystem, not a particle mover.  Its registry group is only a component-test
entry point and does not replace the later dedicated turbulence verification and
validation campaign.

## Initialization input and mesh gates

### Parallel-diffusion / particle-mover runner

`test/run_parallel_diffusion_mover_tests.py` runs, from any directory and
without an AMPS build, the suites that cover the shared parallel-diffusion
library binding and the schema-6 mover key `application.sep.transport`: `library` (library `make verify`,
including PD13-\*), `binding` (PDB01–PDB05), `cli` (CLI01–CLI06), and
`coefficients` (COEF01–COEF06). It parses every per-test result line, prints
per-suite and total PASS/FAIL/SKIP/ERROR/MISSING counts, and lists each failed
test with its suite and message. An expected test that reports nothing is
`MISSING`; a nonzero suite exit not explained by a test failure is recorded as
an `ERROR`; the step-1 suite's intentional `[HIDDEN01] FAIL` negative control
is ignored. Exit status is 0 only when nothing failed. With `--json FILE` (or `--log-dir DIR`) each suite's
complete output is saved to `<json name>-logs/<suite>.log` (or `DIR/<suite>.log`),
and every result in the JSON carries `log` and the 1-based `line` of its result
line, so a failed test can be located directly.

```sh
python3 srcSEP/test/run_parallel_diffusion_mover_tests.py --help   # describes every test
python3 srcSEP/test/run_parallel_diffusion_mover_tests.py          # all suites
python3 srcSEP/test/run_parallel_diffusion_mover_tests.py --suite binding --suite cli
python3 srcSEP/test/run_parallel_diffusion_mover_tests.py --verbose --json report.json
```

`make test-parallel-diffusion-binding-unit` compiles the production
`util/sep_initialization.cpp`, `adapters/parallel_diffusion_adapter.cpp`, the
canonical coefficient registry, and the shared parallel-diffusion library
without AMPS or MPI; its inputs are shared SEP schema-6 files
(`../../src/INPUT_FILE.md` section 5). `PDB01` checks that the schema-6
reader stores the `physics: parallel_diffusion` assignments as `key = value`
records (comments removed by the generic AMPS input reader) with their line
numbers without altering case-sensitive library keys, that the block is
optional, and that a retired schema number (4) and duplicate blocks are
rejected. `PDB02` checks every rejected
selection (block without provider, provider without an input file or block, focused
mover, models needing unavailable host inputs, malformed values with the deck
line number) and that nothing is installed. `PDB03` checks the registry
provider and its Parker-only compatibility. `PDB04` installs models through
the adapter and compares `kappa`/`lambda` returned through
`ActiveParallelDiffusion` with independent closed forms (`v*lambda/3`, the
rigidity/radial power law). `PDB05` checks the schema-6
`application.sep.transport` key: `parker`, `focused-diffusion`, and
`focused-scattering` map to the movers `parker`, `fte-dmumu`, and `fte-mfp`
(expected values written in the test); srcSEP mover names, unknown and empty
values, the retired `application.sep.particle_mover` key ("unknown
configuration key"), a duplicate key, and a missing key are rejected; the
value is part of the startup fingerprint. The PIC sampling wrapper is not
exercised here; it requires a native build.

`make test-initialization-unit` (`test/run_initialization_tests.sh`) compiles
the production `util/sep_initialization.cpp` and the exact native callbacks
without AMPS or MPI, with `-std=c++17` and `-I src/general/input` for the
header-only generic AMPS input reader (`amps_input.h`).  `INIT01` parses a complete schema-6 srcSEP input and checks the installed
values and fingerprint (16 hex digits, independent of SWCME assignment
order); the D4 ownership (Parker winding wind `4.0e5` m/s from
`ambient.wind_speed`, rotation `2.86533e-6` rad/s, SWCME
`source.energy_min = 1.602176634e-14 J` synthesized from `source: shock`);
`application.sep.transport` and the turbulence block (128 spectral bins,
shock efficiency 0.02, diagnostic cadence 100). Eleven negative controls must
each fail with their own diagnosis: a retired schema (`retired`); a
srcSEP3D-only block (`mesh.memory`), key (`run.campaign_seed`), provider value
(`analytic-parker`), and observer kind (`fixed-cartesian`), each `not
supported by srcSEP`; an SWCME `source.*` key in `physics: swcme` (`written in
source: shock`); a `physics: background` / `parker` record; a 3-D SWCME
geometry key; a 3-D SWCME rotation rate different from SWCME1D; the retired
`particle_mover` key (`unknown configuration key`); and an incomplete
turbulence block (missing `spectral_bins`). The SWCME1D rotation rate
written as the 3-D key must be accepted with an identical fingerprint, and
missing fields, unknown blocks, the retired INI syntax, and an incomplete
SWCME model setting are rejected.  `INIT02` checks requested Parker point
count, configured arc length, the surface resolution identity, and the tube
cross-section identity at its reference radius. `INIT03` writes the production
finite-line Tecplot zone, verifies its unit-qualified schema and configured
point count, and removes the temporary product.

The same callbacks are appended by `ComponentTestRegistry()`.  Consequently a
rebuilt linked executable advertises both IDs through `--list-tests`, and the
normal Python `--all` procedure executes them individually.  Do not remove
them from the native registry or replace their negative controls with a source
scan; the focused gate is fast feedback, while the linked registry proves the
packaged production objects contain the feature.

## Conservative population control and native lifecycle gates

`sh test/run_conservation_resampling_tests.sh` builds the same `CRS01`–`CRS11`
callbacks used by the native registry, with strict C++11 warnings, ASan and
UBSan. Set `SRCSEP_REPORT_DIR` to retain its JSON/JUnit reports. The group is
`conservative-resampling`; these callbacks require no PIC initialization.

| ID | Independent gate |
| --- | --- |
| CRS01 | Both positive weight endpoints preserve number, charge, relativistic kinetic energy, parallel momentum, first/second pitch moments and local spatial moment. |
| CRS02 | Pivoted Householder QR matches an analytical sixth finite-difference null vector and its endpoint probabilities. |
| CRS03 | Full-rank sparse populations return explicit no-reduction; cold/identical constraints retain their actual rank. |
| CRS04 | Every phase point has its original conditional expected weight; endpoint selection obeys a fixed five-sigma frequency gate. |
| CRS05 | SI-scale negative controls reject 50% losses at approximately `1e-14 J` and `1e-19 kg m/s`, including signed-momentum cancellation. |
| CRS06 | Invalid physical/numerical inputs and an unmet zero residual budget reject without substituting another model. |
| CRS07 | Current and stress changes are reported separately; their ensemble expectations remain unchanged. |
| CRS08 | Repeated conservative deletion retains original phase points and stops at the constraint rank. |
| CRS09 | A bounded 1000-to-800 particle stress run retains original SI moments; measured duration is reported without a portable speed claim. |
| CRS10 | A predictable stop retains expected weights; a selective whole-pass rollback negative control demonstrates its bias. |
| CRS11 | Immutable source birth energy is conserved independently of current kinetic energy; an omitted-row control changes the birth measure, and missing/mixed/nonfinite histories reject. |

The weight-only algorithm follows [Gonoskov, CPC 271 (2022) 108200](https://doi.org/10.1016/j.cpc.2021.108200)
([primary preprint](https://arxiv.org/html/1607.03755v2)). Its two nonnegative
endpoints lie along a null direction of the moment matrix. Sampling them with
the derived probabilities preserves each input weight in expectation. Existing
phase points survive or receive zero weight; the algorithm creates no new
positions, velocities, gyrophases, identities, RNG histories or hazard clocks.

The fixed same-species, same-segment basis is
`[1, K_rel, p_parallel, mu, mu^2, local_fraction]`. Charge follows number for a
fixed charge state. Each residual uses its own `sum(abs(weight*quantity))` SI
scale, with no dimensionless `max(1,...)` floor. Current and pressure/stress
moments need additional independent constraints to be exact when energies
vary; preserving `mu^2` alone does not conserve the pressure tensor. Production
reports their signed changes explicitly in `C m/s` and `J`, before any sampling
volume normalization.

The production adapter retains configured spatial/parallel-velocity/normal-
velocity bins, using local segment fractions and ordered tuple keys. It does
not coarsen bins to reach a target. Sparse or rank-limited bins report
`sparse-no-reduction`. Both possible endpoints are checked against the original
whole-segment residual budget before drawing. A numerical stop retains previous
validated updates; the failed attempt changes no weights. This predictable,
bounded stopping preserves unbiasedness. A redundant final-gate inconsistency
terminates explicitly. Population randomness uses the persisted source
scheduler's separate `PopulationControl` channel and physical line/species/
segment keys; it never consumes a retained particle's mover stream.

After a native rebuild, `POPN01`–`POPN03` (group `population-native`) exercise
actual PIC storage and the existing core controllers:

| ID | Native gate |
| --- | --- |
| POPN01 | Actual exact splits from 10 to 12 to 14 records retain the true parent's ID, cached RNG, event/hazard history and persisted split counter; children receive their own keyed future. |
| POPN02 | Ordinary `CloneParticle` retains every requested byte and destination linkage, including split serial and cached-normal continuation. |
| POPN03 | Actual conservative merging on a nonzero segment preserves relativistic energy and every survivor byte except weight/list links; a sparse pass retains the exact population without legacy fallback. |

These fixtures restore original lists, temporary lists, allocations, callback
pointers and source counters. The private population pass cache is regenerated
by the next required `BeginConservativePopulationPass`. Native execution still
requires real serial/MPI reports; these tests alone do not qualify checkpoint
I/O, decomposition invariance or a long physical campaign. Distribution/bin/
particle-count and multi-seed refinement remain required for unpreserved
observables, along with moving/imported geometry once its stable physical-line
identity contract is supplied.

### Native test campaign CLI (DEV01)

The existing parser also accepts a selector union through `--tests`: canonical
IDs, group names, explicit tags, and `id:`, `group:`, or `tag:` prefixes. Commas
or quoted whitespace separate selectors; matches are case insensitive,
deduplicated and sorted by canonical ID. `geometry` and `mover:parker` are
explicit catalog tags. The original plan's complete `all` selection is
`--tests all`, including extended cases. The existing `--all-tests` remains the
bounded routine set for compatibility. A complete catalog selection does not
supply reviewed physical inputs or constitute a successful full validation
campaign.

```sh
./amps --tests 'RNG01,RNG04,group:conservative-resampling' --test-seed 123 --test-repeat 2 --test-output /tmp/sep-cli-seeds
./amps --tests 'tag:geometry,mover:parker' --test-output /tmp/sep-selected-tests
./amps --test CLIEX01 --test-repeat 2 --test-format junit --keep-test-artifacts --test-output /tmp/sep-cli-probe
python3 srcSEP/test/run_native_cli_campaign_tests.py --amps /absolute/path/to/amps --output /tmp/sep-native-cli-proof
```

`--test-seed` is a strict unsigned 64-bit decimal value. `--test-repeat` is a
strict positive 32-bit count. Each converted stochastic callback uses a
versioned case/replicate key independent of selection order, MPI rank and
thread scheduling; `TestSeed(default,purpose)` maps its actual stream seed and
retains that mapping. With no seed/repeat request, historical stream defaults
are unchanged. The catalog records `context`, `deterministic`, or
`fixed-or-unconverted` capability plus explicit repeat support. Unsupported
seed/repeat controls produce `ERROR` without executing that callback. Reviewed
CV/IV/XM/OV/EV campaign inputs retain their declared seeds and cannot be
overridden by generic test controls.

Converted families include RNG lifecycle, conservative resampling, pure
mover/neutral resonant tests, and native population fixtures. RNG07/RNG08 and
deterministic endpoint/mover cases consume no random variates; a supplied
campaign seed is therefore inapplicable to their numerical result. All MFP
fixtures draw a persistent optical-depth clock, including ballistic cases,
and are identified as stochastic. The legacy FTE01/PARKER01/PARKER02 fixtures
map active keys to a positive odd 31-bit seed accepted by `rnd_seed`; the
actual int seed is reported. Their override is qualified only for MPI with
CUDA disabled; hybrid/device builds explicitly reject that unsupported RNG
restoration path. The declared `positive-odd31-v1` bridge avoids even
fixed/short-period states of the legacy multiplicative generator, including
the `2^30` seed whose uniforms are all `0.5`. Historical defaults remain
unchanged. The native SOURCE02 fixture has a private actual-core negative
control for this failure mode.

The default evidence root is a timestamp/process-specific directory below
`test_output`, or use `--test-output`. JSON and JUnit are always retained;
`--test-format text|json|junit` chooses the displayed summary. Legacy
initialization diagnostics may precede that summary on stdout, so use the
committed report file for machine parsing. `--test-json` and `--test-junit`
remain exact report path overrides. Output directories are checked before
model initialization. Root alone commits shared summaries; callback outputs
live under `ID/replicate-N`, with a separate `rank-N` directory for non-root
ranks. Existing reviewed one-case `--test-output-dir` remains exact when no new
campaign controls are requested.

JSON retains one stable aggregate ID with complete `invocations` records;
replicate metrics and seeds are never collapsed into the first run. JUnit
names executed cases `ID[replicate=N]`. `--fail-fast` waits for MPI outcome
reduction, stops every rank on the first `FAIL` or `ERROR`, and records compact
unexecuted ranges with their invocation count. Status exits remain 0 for
PASS/SKIP, 1 for FAIL and 2 for ERROR. Worker joins are required before fixture
execution. Native fixture scopes restore supported complete RNG/cache,
source-scheduler, campaign and background state; failed restoration blocks
dependent callbacks and reports `ERROR`.

Declared artifacts and per-rank `native-log.txt` transcripts remain available.
Harness-owned successful scratch directories are removed by default;
`--keep-test-artifacts` retains them, and failed/error scratch is retained.
This flag does not change a callback's private `tmpfile` lifecycle. The native
CLI subprocess proof checks real file retention rather than a metadata flag.

CLIEX01 (group `cli-infrastructure`, initialization `none`) is an
**infrastructure probe, not physics validation**. It performs actual keyed
uniform/normal draws and transactional per-rank artifact writes. Controlled
fault testing requires explicit `--test-input` with schema
`srcsep-cli-fault-v1` and exactly `fault_rank=N` and `fault_replicate=N`; unknown,
duplicate, malformed or out-of-range input is rejected. The exact input is
archived. The small reviewed fixture
`test/fixtures/cli_fault_rank1_rep0.txt` requests rank 1 to fail on replicate 0
of a two-rank campaign, proving non-root failures drive collective fail-fast.
Normal catalog invocation supplies no fault and passes. Native two-rank
evidence must be produced by the subprocess proof after rebuilding; parser
tests alone do not qualify MPI behavior.

The dependency-light RNG and conservative-resampling sanitizer launchers also
run fixed campaigns with seeds 123 and 456, two replicates, and reversed
execution order. They assert changes in **actual stochastic observables**, not
only reported keys. `SRCSEP_REPORT_DIR` retains all default and seeded JSON /
JUnit reports. Mover/neutral resonant families have a separate actual-draw
trace check. Statistical tolerances and physical invariants are unchanged.

Unsupported fixed/unconverted repeats produce one contract `ERROR` and a
compact range for all remaining requested invocations, including a
`UINT32_MAX` request without fail-fast. Later supported descriptors still run
when fail-fast is disabled. The parser test checks this bound and its exact
not-run count.

`run_native_cli_campaign_tests.py` writes transactional `summary.json` with
stable infrastructure IDs NCLI01–NCLI09, each command, expected/actual process
exit, individual postchecks and process-log path. An expected child failure is
an infrastructure PASS only after its provenance and collective stopping
checks pass. Launch/timeout/postcheck failures are recorded before propagation,
and blocked jobs are SKIP. The `native-cli-campaign` orchestrator step reads
these individual records. Its reader rejects a false PASS containing a failed
postcheck. Independent subprocess negative controls test both launch and exit
failures; they do not count as native model execution.

The `production-cli` native group contains PCLI01–PCLI12. It calls the same
production parser and `ApplyTurbulenceOptions`, then checks the actual coupling,
cascade/reflection, injection count, mover registry, wave representation,
coefficient provider, tolerances and diagnostic cadence. It restores every
changed production variable and requires zero prior mover-dispatch evidence
before mutation, since selecting a mover resets that counter. The tests cover
defaults, help, frozen Boolean aliases, case/value syntax, numeric boundaries,
errors, actual application, input-layer precedence, configuration warnings,
legacy TestManager control mapping, 192 applied option combinations, and parser
repetition.

PCLI09 checks warning output against the preserved advertised capability table.
All three selectable movers currently advertise streaming accumulation, so
the historical missing-capability warning branch is unreachable and remains
explicitly unexercised. That metadata is not proof of physically qualified
four-spectrum streaming; the separate production closure gate remains active.

```sh
python3 srcSEP/test/run_native_production_cli_tests.py --amps /absolute/path/to/amps --output /tmp/sep-production-cli
```

This companion runner executes both MPI ranks twice, retains individual native
PCLI results, checks fresh executable help/errors with a hostile `amps.in`, and
records source hashes for the limited input-order/coupled-ownership callsite
audit. The orchestrator step is `production-cli-compatibility`. **These software
checks do not complete the original PCLI11 production-evolution or PCLI12 full
SWMF-run gates.** Those missing gates remain explicit SKIP records. The heavy
legacy TestManager is not silently executed by a registry test. Its parser
control mapping is covered; executing that historical diagnostic remains a
separate gate. Current canonical mover names follow the documented Step14
migration, while the original misspelled cascade alias remains accepted.

An enabled legacy TestManager and component execution selector now reject as
conflicting modes. Previously, the test-only return silently bypassed that
enabled legacy request. Unmixed legacy diagnostics retain their behavior;
`--run-test-manager --no-test-manager --test PCLI01` remains valid under the
existing last-occurrence-wins rule. Parser and linked PCLI10 checks verify the
conflict without invoking the heavy diagnostic.

### Independent physical normalization and polarized wave values

`make -C srcSEP test-real-wave-normalization-unit` runs RWN01/RWN02, the
independent real-field normalization gates. RWN01 integrates full Lorentz
orbits in declared real circular fields and compares finite-time Kubo variance
and time-step refinement with the one-helicity resonant coefficient. RWN02
uses explicitly balanced independent real L/R fields and positive-wavenumber
spectrum integrals to check the Jokipii and static Florinskiy normalization
against total magnetic variance. These are deterministic prescribed-wave
experiments, not evidence that the native particle-wave adapter is fully
qualified. They are extended native cases in group `resonant-collision`; the
orchestrator source step is `real-wave-normalization`.

`make -C srcSEP test-polarized-wave-state-unit` runs the private immutable-value
tests in group `polarized-wave-state`; the orchestrator step has the same name.

| ID | Value contract |
| --- | --- |
| PWS01 | Missing directional/helicity allocation rejects; polarization must be supplied explicitly. |
| PWS02 | Declared logarithmic centers and width give geometric internal edges and outer half bins. |
| PWS03 | Exact finite-bin integrals preserve initial branch energy with signed rounding residuals. |
| PWS04 | The resonant kernel consumes the authoritative four directional/helicity arrays. |
| PWS05 | Invalid initialization leaves the previously published immutable value unchanged. |
| PWS06 | Native spectral count/grid mismatches reject before source reads or destination mutation, including counts above allocated NK. |
| PWS07 | Restore owns an exact immutable copy of all four arrays and signed ledgers; changed bins, grid, revision or ledger reject. |

These value tests do not prove production ownership across a native coupled
step. In particular, restoring evolved helicities must use the stored arrays;
an initial equal-L/R prescription is not a replacement for evolved state.

### Narrow native Maxwellian birth geometry evidence

SOURCE02, in `source-native`, exercises the existing generic PIC field-line
Maxwellian routine with an explicitly short diagnostic interval and global
time-step/statistical-weight mode. It forces the segment owner to differ from
the AMR endpoint owner on multiple ranks, checks two actual initialized births,
and verifies physical vertex coordinates, attachment, fresh extension state and
post-birth stable identities. The fixture restores geometry, owner, numerics,
full legacy RNG state and allocations. A missing fixed line, nonempty process,
or unsupported local numerical mode is an explicit SKIP.

This is **geometry/initial-owner software evidence only**. It excludes the
physical plasma-versus-lab injection frame, within-step birth-time law,
segment-crossing routing, purpose-keyed Maxwellian sampler, decomposition-
independent sampling and production source normalization. The short interval
prevents a segment crossing; it does not validate crossing behavior. The
production Maxwellian source remains unqualified while those contracts are
unresolved.

### Explicit birth windows and native four-spectrum gates

`bash srcSEP/test/run_birth_window_tests.sh` runs BIRTH01–BIRTH11 with strict warnings,
ASan and UBSan; `SRCSEP_REPORT_DIR` retains its JSON/JUnit. The `birth-window`
orchestrator step reads each stable ID, rather than counting the executable as
one test. These are time-lifecycle value tests, with seconds as the time unit;
they do not qualify an untested native source adapter or choose its physical
birth-time law.

| ID | Time contract |
| --- | --- |
| BIRTH01 | The driver supplies a finite, representable complete host-step window. |
| BIRTH02 | A late birth retains its exact positive remaining duration even when its displayed absolute epoch rounds to the step end. |
| BIRTH03 | First-movement consumption is staged; the original birth history stays unchanged until commit. |
| BIRTH04 | Unknown history, invalid phases, mismatched windows and unequal host durations reject without output mutation. |
| BIRTH05 | Immutable snapshot coverage includes the complete parent window and the birth-relative movement interval. |
| BIRTH06 | An explicit initial condition is distinct from a scheduled source birth. |
| BIRTH07 | A prescribed positive terminal remainder is completed; an accuracy-requested subdivision below the numerical floor still rejects. |
| BIRTH08 | An externally owned clock without an explicit window remains unsupported; standalone clock capability is unchanged. |
| BIRTH09 | A bound external window requires the physical departure, full nominal duration, and the audited nonadvancing PIC clock policy. |
| BIRTH10 | Complete fixed-duration steps follow nondivisible host limits without counting skipped calls as sampling exposure. |
| BIRTH11 | Timestep mutation, premature/duplicate phase transitions, invalid limits and unrepresentable endpoints reject. |

BIRTH09–11 exercise the same value/phase contract used by the SWMF wrapper;
their private exposure counter is a scheduling oracle, not qualification of the
native sampler. The original eight-rank coupled run and unchanged energy and
pitch-angle comparisons remain mandatory. Its configured timestep must not be
changed to fit coupling intervals. A skipped call contributes no physical step
or accumulated particle sample. Real evolving nondivisible backgrounds and full
PT restart continuation require their own explicit provider/restore contracts.

`make -C srcSEP test-maxwell-juttner-volume-unit` runs MJV01–03 with strict
warnings, ASan and UBSan. The same callbacks are linked into the native registry
as group `maxwell-juttner-volume`. MJV01 compares cold and relativistic energy
moments against independent quadrature and the ideal-gas pressure identity;
MJV02 checks the boosted simulation-volume number and stress-energy measure;
MJV03 checks invalid-input rejection without RNG/output mutation. Campaign-seed
controls verify that stochastic callbacks actually consume their declared seed.
These checks qualify the momentum kernel, not a complete native SWMF population.

CP03G05–07, run by `bash srcSEP/test/run_background_derivative_tests.sh`, test
the production receive-time reconstruction kernel on full Cartesian affine
fields, coplanar negative controls and a symmetric cubic convergence sequence.
Its condition-number gate bounds roundoff amplification, not spatial truncation
error on arbitrary imported data. Native packet units, epoch consistency,
initialization, transport, sampling and original reference comparisons remain
separate integration obligations.

The tiny-remainder negative control uses `nextafter(1,0)` at epoch `1e12 s`.
It preserves the positive duration without clipping or redrawing the supplied
phase. Permitting that prescribed remainder does not relax adaptive error
controls.

PNW01–PNW04, group `polarized-wave-native`, require spectral storage to be
registered **before** field-line initialization. The default integrated model
cannot exercise these native buffer contracts. Run the explicitly configured
two-rank gate, also registered as orchestrator step `polarized-wave-native`:

```sh
mpiexec -n 2 ./amps --test-group polarized-wave-native --turbulence-model wave-number-resolved --turbulence-source self-consistent-spectral --coupling off --cascade off --reflection off --test-output /tmp/sep-polarized-native
```

| ID | Native contract |
| --- | --- |
| PNW01 | Actual four-spectrum values reproduce finite-bin and branch moments, initial signed rounding ledgers, and MPI diagnostic readback. |
| PNW02 | Malformed collective initialization rejects before replacing the owner or changing host diagnostic buffers. |
| PNW03 | Staged restore retains the exact evolved four arrays and ledgers; it does not reconstruct initial equal helicities. |
| PNW04 | Invalid amplitudes, terminal magnetic field and radial endpoints reject before host writes; explicit zero amplitude remains zero. |

These are native initialization/restore contracts, not qualification of coupled
physical evolution. The orchestrator expects all four cases in this configured
step. A missing case or prerequisite SKIP becomes a step ERROR while the child
report preserves its original status and explanation; an accidentally omitted
spectral flag cannot produce a green native gate.

### Validation-only warm parallel Maxwellian study

`make -C srcSEP test-warm-parallel-study-unit` runs the Python mathematical
verification suite. The `warm-parallel-study` orchestrator step uses the actual
verbose unittest reader, preserving each test name, status and failure block.
It covers independent Landau-contour integration, signed cold-plasma limits,
the first warm MHD correction, precision/continuation refinement, extreme
neutral cancellation, serialization round trips, winding aliasing, singular
predictor rejection and mandatory equilibrium inputs.

The separate study requires explicit proton/electron densities, temperatures,
rest masses, charge magnitude, magnetic field, zero drifts, physical constants,
distribution, propagation and provenance. The small checked-in manifest is a
declared controlled verification case, not an inferred production equilibrium:

```sh
python3 srcSEP/validation/physics/warm_parallel_maxwellian.py --case srcSEP/validation/physics/controlled_ep_maxwellian_baseline.json --kdi .001 .01 .1 1 --helicities Right Left --digits 40 --verification-digits 70 --output /tmp/sep-warm-controlled.json
```

The example requests study points; it selects no accepted production domain.
The solver retains displacement current, both species' inertia/thermal
responses and the entire Landau continuation. It records convergence-based
frequency uncertainty, output rounding, raw requested points, input hash and a
content-addressed implementation archive. A serialized-frequency regression
checks the emitted values against their reported uncertainty at higher
precision. Winding uses noncommensurate grids plus an independently averaged
analytic logarithmic derivative; a manufactured 65-root negative control
detects the former doubled-grid alias.

Root isolation remains numerical, with no interval-arithmetic proof. Damping
below estimated numerical uncertainty is unresolved, not exact zero. The
Hermitian wave-energy ratio is a weak-damping diagnostic; no damping acceptance
threshold is supplied. The declared model is homogeneous, isotropic,
nonrelativistic and strictly parallel; it does not establish oblique,
inhomogeneous, nonlinear or observational validity. Study artifacts keep
`physical_production_qualified=false` and `domain_limits_selected=false`.

## Frozen accepted-proposal transaction: FCT01–FCT09

`SEP::Transport::FrozenCollisionTransaction::Prepare` is a neutral bookkeeping
component for an explicitly declared frozen-wave operator batch. It accepts
already weighted resonant-collision proposals, pins their snapshot, epoch,
geometry, background, physical cell, parent state, lifecycle token and stable
particle/event key, and calls the canonical `RC::ValidateProposal` for the
elastic orbit and branch/helicity/bin partition. It consumes no random numbers.
The enclosing driver supplies the accepted-event phase identity: this helper
does not reprove stochastic acceptance or neighbor membership in an unspecified
generator, and it does not certify the owner's opaque RNG/hazard token.
Exact encoded parent momenta must match: an approximate inverse Lorentz
transform cannot silently replace the parent's actual stored phase point.
Identity zero remains the native/source missing-identity sentinel; event zero
is the valid first accepted event.

Events reduce in stable particle/event order into canonical physical bins.
Statistical weight is applied once by the collision proposal and checked at
the transaction seam. All net magnetic bin budgets are checked before any
candidate is exposed. A negative budget or an unresolvable nonzero stored
transfer returns `UnresolvedCoefficient`; no clipping, accepted-event removal,
redraw or partial update occurs. The caller supplies a dimensionless numerical
tolerance. Each bin's stored-change roundoff must satisfy that tolerance
relative to its requested net transfer, with no SI floor or reservoir-sized
normalization. Ordinary bounded roundoff has a signed ledger and is not heat.
Total nondispersive Alfvén wave energy is twice magnetic energy. The canonical
proposal's requested invariant increment is kept separate from the actual
energy difference of the encoded before/after particle momenta. That difference
uses a stable relativistic difference-of-squares calculation without subtracting
rest-mass-dominated energies. Every event's particle rounding is signed and must
resolve its requested transfer at the same explicit tolerance. A nonzero
requested increment represented as zero, the wrong sign, or unresolved encoded
energy returns `UnresolvedCoefficient` before any wave debit. The total actual
residual uses actual encoded particle energy plus twice the actual stored
magnetic change and must also pass the explicit numerical tolerance.

The helper returns staged particle records plus an immutable wave candidate
from `Polarized::PrepareEvolutionCandidate`, which alone owns canonical wave
hashing/revision. It has no public commit operation and owns no production
state. A private test-only bundle demonstrates an atomic pointer swap after
precondition checks. Future native publication must coordinate actual PIC
particles, lists, RNG/hazard state and wave ownership in one outer transaction.
Frozen-wave batching does not recompute hazards after intermediate deposits
and is not an exact coupled CTMC or a production-qualified collision driver.

| ID | Actual focused contract |
| --- | --- |
| FCT01 | Weighted SI particle energy agrees with an independent relativistic reference; both propagation branches close against twice-magnetic wave transfer. |
| FCT02 | Ordinary storage roundoff is bounded and signed; sub-ULP debit fails unresolved with no particle or wave candidate. |
| FCT03 | A whole batch containing an over-budget event is rejected after private staging, with unchanged inputs. |
| FCT04 | Stale generation/epoch/identity, mismatched phase/weight/lifecycle, one-ULP wrong parent, and tampered orbit or branch/helicity/bin exchanges fail. |
| FCT05 | Exact chained parent states and consecutive events succeed; duplicate/missing/overflow events and reserved particle identity zero fail. |
| FCT06 | Reversed particle/event input gives identical staged values, signed reductions and canonical snapshot identity. |
| FCT07 | Missing frozen phase, stale snapshot/background contracts and invalid explicit numerical tolerances fail. |
| FCT08 | A private fixture rejects stale wave/particle publication preconditions and swaps one complete immutable bundle. |
| FCT09 | Actual encoded particle energy and signed roundoff are checked; an RC-valid roundoff-scale wire with resolvable wave debit but zero encoded particle transfer fails unresolved. |

```sh
make -C srcSEP test-frozen-collision-transaction-unit
python3 -B srcSEP/test/run_tests.py --suite frozen-collision-transaction --output-dir /tmp/sep-fct-source --plot none
python3 -B srcSEP/test/run_tests.py --amps ./amps --group frozen-collision-transaction --output-dir /tmp/sep-fct-native --plot none
python3 -B tools/sep_test_orchestrator.py --app srcsep --skip-build --only frozen-collision-transaction --output-dir /tmp/sep-fct-orchestrator
```

The source command compiles the same native-selectable deterministic callbacks
with C++11 strict warnings, ASan and UBSan in `/tmp`. Native commands require an
executable rebuilt with this factory. The dedicated orchestrator reader retains
individual FCT PASS/FAIL/ERROR/SKIP records; a source PASS makes no PIC-commit claim.
The focused source runner records hashes of all source-owned translation units,
local headers and its own script before compilation and after compilation and
execution. A changed input discards the entire run with nonzero status.

## Native source birth lifecycle: BIRTHN01–BIRTHN06

The `birth-window-native` group retains its six stable required IDs. BIRTHN01,
BIRTHN02 and BIRTHN03 exercise the actual Parker, Dmumu and MFP adapters;
BIRTHN04 retains invalid-window/history/Cartesian atomic rejection and the
strictly positive late interval at epoch `1e12`; BIRTHN05 retains generic
clone, actual MPI exchange and exact split history checks. BIRTHN06 tests the
real completed-state checkpoint boundary. These are lifecycle/operator
regressions, not a qualification of a new source sampler or a different time
integration method.

Use `run_birth_operator_stage_tests.py --amps /absolute/amps --mpi-ranks 2
--output /absolute/new-directory`. Each operator experiment is a fresh process
initialized with the same physical AMPS timestep `D=0.01 s`; no callback writes
`GlobalTimeStep`. For BIRTHN01–03 the runner executes independent `direct`,
`full`, and `actual` roles. All three consume the exact same frozen departure
realization R0: geometry/owner metadata, provider and coefficient identities,
all consumed vertex/segment bytes, and initial particle/RNG/hazard/identity
bytes must match. The actual role freezes the full parent window and creates
a phase-0.75 birth with the unchanged exact remaining duration `h=0.0025 s`.
The direct control calls the existing native diagnostic dispatcher for `h`,
without a physical source window, counter increment or sampling. Its complete
particle image must match the actual pending-birth operator, excluding only
the original birth-record/list-link exclusions. The independent full-`D`
control must differ. Missing controls, changed inputs, native SKIP, contradictory
exit/status/metrics, or unequal byte images fail the gate.

This comparison grades the maintained departure-snapshot/operator-split
method. It does not claim that a fresh physical simulation starting at the
birth epoch uses the same fields. The diagnostic owns the actual declared
nonzero radial wind, and after readers join it publishes the same complete
material map `M_D` through `PrepareGeometryStep`/`CommitGeometryStep` in every
role. Existing vertices, including both endpoints, move; inner extension
traces along the pushed magnetic field, new segments start empty, and normal
outer retirement applies. The initial curved Parker field is pushed using
`B=F B0/J` and number density using `n=n0/J`; full radial strain, independent
Cauchy-field/continuity checks, owner epoch/generation, and exact
`PB::X == GetCartesian(S)` are checked. The Cartesian assertion runs once on
the frozen operator chart and again on the completed material chart. Thermal
and wave evolution are outside this controlled fixture. At epoch `1e12`, the
map captures the explicit physical `D`, retaining the positive birth remainder
separately; rounded endpoint subtraction never supplies either duration.

BIRTHN06 uses the existing explicit reduced-material checkpoint protocol2
startup and its unchanged initialized timestep. It injects genuine source
records, completes their actual mover, absorbing-boundary, source-ledger and
material-geometry transactions, and requires a nonempty surviving source
cohort. A positive completed save establishes ordinary prerequisites. Fault
injection then changes only the pending-first-move bit of those actual birth
records: source identity, birth energy and cohort bytes remain unchanged.
The save must reject for the original unconsumed-movement reason and preserve
both the full particle image and an existing destination. Active-window save
and load retain their separate rejection/nonmutation checks. No source history
or material-checkpoint capability is fabricated.

`--repeat-count N` runs complete independent role groups in fresh processes,
preserving the registry's existing campaign/ID/replicate seed schedule. The
native descriptor retains repeat capability with `repeat:process-isolated`;
a raw in-process repeat request returns an explicit ERROR directing callers to
the runner. Per-role requests/logs/reports and binary images remain under the
owned output, together with detailed `summary.json`, standard
`srcsep-tests.json` and JUnit output. Every requested role and repetition must
pass; no required ID is dropped or converted to SKIP.

## Local numerical MFLAMPA comparison: verified feasibility and required contracts

The local `SWMF_data/SP/MFLAMPA/data` files are archived numerical regression
inputs and outputs. They are distinct from the digitized MFLAMPA publication
curves used by XM02/XM03 and from measured detector observations. A source/build
audit found a clean local MFLAMPA checkout at official commit
`38e9defc9c525bd3d49669bc2e7d8014987de488`; an isolated standalone build from
that commit and pinned local `share`/`util` sources compiled and linked in
`/tmp/srcsep-mflampa-audit-kbou2cbd`. The exact named upstream `test_mflampa`
regression subsequently passed on two MPI ranks with its unchanged parameter
file, full grid/defaults, and all available input states through 600 s. The
actual solver exited zero in 9.79 s; 33 selected output files/blocks passed
75 setup/header/shape/name/time checks, and the upstream
`DiffNum -BLESS=NO -t -r=1e-6 -a=1e-6` comparison exited zero with an empty diff.
The initial wall-time bound was 120 s; no physics or source reductions were
made. Protected source/binary/input/reference hashes remained unchanged.
This is upstream software-regression evidence, with no AMPS cross-model score,
observational validation or adopted polarized-gap qualification.
`validation/reference_sources.json` records the source commits,
executable checksum, binary-input coverage and detailed feasibility evidence.
This new executable is not identified as the historical reference producer.

The [official MFLAMPA Makefile](https://github.com/SWMFsoftware/MFLAMPA/blob/38e9defc9c525bd3d49669bc2e7d8014987de488/Makefile)
maps the local references as follows:

| Upstream target | Parameter realization | Archived numerical output |
| --- | --- | --- |
| `test_mflampa` | `Param/PARAM.in.test`; Parker upwind | `MH_data.ref.gz` |
| `test_poisson` | `Param/PARAM.in.test.poisson`; Parker Poisson brackets | `MH_poisson_data.ref.gz` |
| `test_steady` | steady start and restart parameter files | `MH_steady_data.ref.gz` |
| `test_spectra` | `Param/PARAM.in.test.spectra` | `MH_spectra_data.ref.gz`, `Distr_spectra_data.ref.gz`, `Distraj_spectra_data.ref.gz` |

The existing AMPS `Makefile.test` target
`test_SEP--focused_transport-IH-mflampa-fl` configures a coupled SC/IH/PT build
and compares AMPS energy/pitch histograms against `output/test_SEP3D` references.
It does not establish agreement with these standalone SP numerical references.
Its configuration/build/run recipes modify their build trees, so they must not
be launched in the user's current checkout as an isolated comparison.

The checked local input archive contains 112 numerical files: 16 field lines
at seven actual times, `0,100,...,600 s`. All 112 files passed a source-guided
Fortran-record/header/shape audit: little-endian four-byte record markers,
binary64 arrays, one Lagrangian coordinate, 13 state variables, six scalar
parameters, and 3942–4200 vertices. The control header and tag list declare
13 epochs through 1200 s, but the files after 600 s are absent. The official
two-session stop at 600 s fits the available coverage; the tag list cannot
justify extrapolation or a longer run.

Source-resolved settings also matter. The official test uses
`nVertexMax=20000`, `nMomentum=100`, `nPitchAngle=1`, initial integral flux
`0.01 pfu`, finite kinetic-energy range `10 keV`–`10^7 keV`, default parallel
`sokolov2004` diffusion, and disabled perpendicular diffusion/drift. Its second
session's `MeanFreePath0InAu=0.3` selects an upstream normalization:
`lambda_parallel = 0.3 AU (r/AU) (p c / 1 GeV)^(1/3)` where
`r > 1.1 r(iShock)`. It is not a constant 0.3 AU at every location and energy.
The near-Sun/momentum boundary policies, thermal injection, nonshock injection
floor, solver switches and session transition must be preserved explicitly.
These are verified defaults of the inspected source, not proof of the exact
compiler/dependency tuple that generated each archived reference.
[The upstream diffusion implementation](https://github.com/SWMFsoftware/MFLAMPA/blob/38e9defc9c525bd3d49669bc2e7d8014987de488/src/ModDiffusion.f90)
also consumes `deltaB^2 = mu0 (Wave1 + Wave2)`; units `J/m3` alone do not
identify magnetic versus total wave energy or provide circular polarization.

An equivalent native replay needs a reviewed input contract before scoring:

1. Pin producer/dependency/compiler versions and physical constants. Validate
   every field's units and full Lagrangian vertex/line identity across actual
   snapshots, including offsets, new nodes and three-dimensional geometry.
2. Resolve density/composition and velocity-frame semantics from producer code
   and provenance. The current coupler divides SI mass density by proton mass;
   the archive header's `amu/m3` label cannot independently define the import.
   The current MFLAMPA file reader explicitly treats stored temperature as
   thermal energy in keV. The PT coupling default expects Kelvin, so those
   input arrays cannot be fed directly to that interface without conversion.
   Use explicit UTC/calendar, Julian date and simulation-time provenance;
   `StartTime` is seconds from the compiled 1965 epoch, not Unix seconds.
3. Implement a declared, read-only offline background adapter with immutable
   epochs. Match the initial distribution, both spatial and momentum boundaries,
   represented volume/flux measure, source normalization and observer projection.
4. Approve the common scattering experiment. A prescribed finite-lambda Parker
   diagnostic tests transport numerics; it does not validate the adopted
   four-spectrum finite-band polarized collision law. MFLAMPA's optional FTE
   implementation floors its near-zero-pitch Kolmogorov scattering factor;
   suppressing native gaps or inferring four helicities would change the model.
5. Predeclare matched time/energy/channel bins, sampling uncertainty and
   space/time/momentum refinement gates. Start with separately identified
   upstream regression reproduction, then equivalent-input comparison; neither
   is observational validation. Do not fit an arbitrary global amplitude.

The comparison wrapper must capture the actual solver exit separately from
`tee`, enumerate/sort output files explicitly, and validate headers, shapes and
provenance before numerical tolerances. The upstream aggregate Make target
ignores some child errors, its checks use brace globs under `/bin/sh`, and
`DiffNum -t` ignores text between numbers. A zero Make exit alone is insufficient
evidence. Copies of immutable references must remain outside any `BLESS=YES`
path; the isolated build's `data` directory is deliberately private and empty.

## Flux-owner, controlled WKB and relativistic material-plane source suites

The Python catalog exposes the following source suites. Each invokes the
existing Make target and source runner; none cleans, configures or relinks the
enclosing native AMPS executable. The `run_all_test` orchestrator exposes the
same three step names and retains each callback's bracketed verdict. Its
readers require every listed ID: missing callbacks and SKIP become ERROR for
these explicitly requested source-family gates. Compiler, sanitizer and child
process failures also propagate through the outer step exit status.

| Python suite / orchestrator step | Source runner | Required callback IDs and evidence scope |
| --- | --- | --- |
| `flux-tube-state` | `test/run_flux_tube_state_tests.sh` | FLCP01–04: exact explicit magnetic-flux owner values/provenance, physical-state rejection, transactional codec rejection, and `A=Phi/B` SI area normalization. These are private value/codec tests, without native owner mutation. |
| `polarized-wkb` | `test/run_polarized_wkb_tests.sh` | WKB01–08: immutable incoming four-spectrum reservoirs, controlled ideal-MHD work/ray escape, independent spatial/spectral convergence, ALE geometric conservation, rejection, material side pressure and exact multiscale external accounts. Native material-grid publication remains separate. |
| `maxwell-juttner-flux` | `test/run_maxwell_juttner_flux_tests.sh` | MJ01–06: proper-density normalization, material-plane Lorentz measure, radial/angular sampling, cold limits, atomic numerical rejection and keyed six-draw lifecycle. This is a full-particle mathematical kernel; the gyroaveraged native source bridge remains gated. |

For example, run one source gate or select all three orchestrator steps:

```sh
python3 -B srcSEP/test/run_tests.py --suite flux-tube-state --output-dir /tmp/sep-flux-state --plot none
python3 -B srcSEP/test/run_tests.py --suite polarized-wkb --output-dir /tmp/sep-wkb --plot none
python3 -B srcSEP/test/run_tests.py --suite maxwell-juttner-flux --output-dir /tmp/sep-mj --plot none
python3 -B tools/sep_test_orchestrator.py --app srcsep --skip-build --only flux-tube-state --only polarized-wkb --only maxwell-juttner-flux --output-dir /tmp/sep-new-source-gates
```

The FLCP native descriptor group is currently `flux-checkpoint`; the Python
source-suite alias and orchestrator step are `flux-tube-state`. WKB and MJ use
the same group names as their Python aliases. A passing catalog/reader self-test
proves discovery and failure accounting, not the mathematical callback verdicts.
In particular, independent review found an extreme cold-boost encoded-momentum
counterexample now included in MJ05. The kernel must reject a thermal crossing
that binary64 momentum cannot resolve before committing output or RNG. Track
that correction and its actual re-execution in
[`MAXWELL_JUTTNER_FLUX_IMPLEMENTATION.md`](../validation/physics/MAXWELL_JUTTNER_FLUX_IMPLEMENTATION.md)
and [`STATUS.md`](../STATUS.md); historical six-case PASS evidence does not
cover the added negative control.

## Paired fresh-process native checkpoint gate

The `fresh-process-checkpoint` orchestrator step runs
`test/run_fresh_process_checkpoint_tests.py` on two MPI ranks. CP05F01, CP05F02
and CP05F03 cover the moving prescribed reduced-front Parker, Dmumu and MFP API
continuations respectively. Each case launches a write job and then a separate
resume job. Each role owns a distinct artifact directory and private PIC cwd;
only the explicitly named checkpoint/final artifacts cross the process boundary.
An ordinary single `--test CP05F01` invocation cannot establish this paired
protocol. A missing material fixture now produces structured ERROR before
allocation; it cannot qualify through an absent role or a frozen nonzero-U line.

```sh
python3 -B srcSEP/test/run_fresh_process_checkpoint_tests.py --amps ./amps --mpi-ranks 2 --output /tmp/sep-fresh-checkpoint
python3 -B tools/sep_test_orchestrator.py --app srcsep --skip-build --only fresh-process-checkpoint --output-dir /tmp/sep-fresh-checkpoint-orchestrator
```

Choose new or empty output directories. The runner retains
`summary.json` with schema `srcsep-fresh-process-checkpoint-v1`, the executable's
initial SHA-256, role exit/status/duration, child JSON reports, complete native
logs and owned working directories. Its own per-role postchecks require an
unchanged executable, zero exit and the requested native PASS. Final
`physics.expected` and `physics.resumed` files must exist and match byte for
byte. The summary is published by temporary-file rename after all selected
cases have been recorded.

The orchestrator independently reads that evidence. It requires all three case
IDs exactly once, both write/resume roles in order, integer zero exits, both
role statuses PASS, exactly the requested PASS in each child native report,
the current executable matching the initial hash, the runner's equality check,
and actual final-file byte equality. Missing/unreadable reports or final files,
SKIP, duplicate/missing cases, changed binaries and contradictory PASS labels
cannot satisfy the gate. Keep the binary and retained files available while
reading the result; a removed/replaced binary is explicitly unavailable evidence.
Harness tests damage each prerequisite independently while leaving the summary's
PASS labels intact. These synthetic fixtures test software accounting only.

An independent application of this strict reader to the retained build-8
three-case campaign also passed all three cases, including all six actual
write/resume child reports and final-byte comparisons. Its executable was pinned
outside the mutable build path with SHA-256
`049a1959918dd154e5a81d75ed7ebb18a57764f29c09c3917c81b4068f9c58eb`.
This evidence establishes that campaign's declared diagnostic API scope only.

This gate establishes the declared fresh-process diagnostic API continuation.
It does not establish a production resume CLI, output-file continuation,
MPI repartitioning, material-grid evolution, adopted four-spectrum collision
physics or native relativistic source qualification. Those require their own
input contracts and independently executed gates.

## Fixed analytic radial startup and geometry-only Parker construction

`analytic-radial` is a source suite/orchestrator step using
`test/run_analytic_radial_background_tests.sh`; its strict C++17, ASan and UBSan
gate compiles the actual shared radial helper and initialization builder. It
hashes 23 listed source inputs before compilation, after compilation and after
execution, and discards verdicts if those inputs change. Set
`SRCSEP_REPORT_DIR` to retain its JSON/JUnit reports and input-hash record.
The native `analytic-radial-native` step uses two MPI ranks through the isolated
Python runner and requires all four native cases; missing or skipped native
cases cannot satisfy this explicit gate.

| ID | Scope |
| --- | --- |
| ARB01 | A translated 3-4-5 radial direction with explicit 725 km/s speed, including invariance under a joint point/origin translation. |
| ARB02 | Independent Cartesian projection-matrix tensor, `div U=2v/r`, radial zero strain and transverse `v/r` strain. |
| ARB03 | Central differences of velocity independently recover the full tensor, with a declared displacement and dimensionless truncation/roundoff gate. |
| ARB04 | Invalid speed/origin and undefined radial positions reject without changing output; a legacy-origin/wind mismatch is rejected rather than supplying another provider. |
| ARB05 | Geometry-only Parker construction needs no invented run/source/output fields, while complete configuration validation and rejected-output atomicity remain strict. |
| ARB06 | The existing canonical background hash distinguishes adjacent representable origin/speed values and preserves physical signed-zero equivalence. |
| ARBN01 | Read-only inspection of actual generated analytic startup U and exact current/previous B/U/n agreement with the existing inverse-square density. |
| ARBN02 | Actual translated-origin/nondefault-speed Parker generator and PIC registration, followed by exact retention of arbitrary imported list U; appended fixture lines are deleted and the original registry count restored. |
| ARBN03 | A finite early history sentinel on every rank and bad later current U on only rank one (rank zero in serial) require collective whole-transaction rejection with byte-identical B/U/n arrays everywhere; every modified datum is restored afterward. |
| ARBN04 | Actual frozen run mover values and provenance are varied/restored while the background identity stays fixed; runtime mover selection and dispatch diagnostics remain untouched. |

```sh
make -C srcSEP test-analytic-radial-background-unit
python3 -B srcSEP/test/run_tests.py --suite analytic-radial --output-dir /tmp/sep-radial-source --plot none
python3 -B srcSEP/test/run_tests.py --amps ./amps --group analytic-radial-native --mpi-np 2 --output-dir /tmp/sep-radial-native --plot none
python3 -B tools/sep_test_orchestrator.py --app srcsep --skip-build --only analytic-radial --only analytic-radial-native --output-dir /tmp/sep-radial-gates
```

The recorded before-fix native run reached a clean ARBN01/ARBN02 FAIL, exit one:
startup previous B/U did not match current storage, and the configured generator
was blocked by unrelated full-runtime validation. An earlier diagnostic run
also retained that generator's fatal PIC exit/signal-handler failure; it is
failure evidence, not a clean component verdict. The corrected builder uses one
geometry validator shared by full validation and geometry construction. Source
ARB05 failed before that split and passed afterward, while INIT01–03 and the
shipped input still pass their complete configuration gates. Native after-fix
results belong to their actual executable/source manifest in
[`STATUS.md`](../STATUS.md); source passes alone do not establish native PASS.
The build-8 two-rank report now records ARB01–06 and ARBN01–04 PASS, with both
rank dispositions retained. ARBN01 inspected 4001 actual startup vertices with
zero reported velocity error; ARBN03 passed its one-rank invalid-endpoint
collective rejection, and ARBN04 preserved background identity across frozen
mover/provenance changes. The retained build-8 binary hash is given in the paired
checkpoint section above. These are the bounded fixed-background gates, not an
extension to imported/moving backgrounds or production closure.

The helper's explicit wind declaration supplies `U=v rhat` and
`G=(v/r)(I-rhat*rhat)` about the same origin. It creates no mutable background
authority or thermal profile. At native fixed-analytic startup, the whole
realization is validated before the existing `n(1 AU)=5e6 m^-3`, inverse-square
density and identical initial previous B/U/n are committed collectively. The
application background identity schema is version three: a private identity-only
run projection excludes just mover value/source, without changing the complete
run/checkpoint fingerprint or executing its canonical sentinel. Registration retains
imported U but does not infer its derivatives or previous physical epoch. These
tests do not qualify imported/moving material backgrounds, a source frame,
temperature closure, full wave physics or production collision closure. The
existing schema-6 SWCME provider selection is unchanged.

## Controlled material-tube and standard BATSRUS wave-view suites

The source suite `material-tube-kinematics` maps to
`test-material-tube-kinematics-unit` and
`test/run_material_tube_kinematics_tests.sh`. Its shared descriptor group is
`material-tube`, with MTG01–06. These checks cover pointwise Cauchy induction and
Nanson area, exact material dilation/mass/cap-side geometric conservation,
shear/translation/tilted caps, independently integrated curved-reference
geometry, the explicit nonaffine cell limitation, and atomic invalid-geometry
rejection. The `material-tube-native` orchestrator step executes those same six
value callbacks in the linked native registry; that proves packaging/execution
and does not create a PIC material-grid publication path.

The source suite `standard-batsrus-waves` maps to
`test-standard-batsrus-wave-unit` and
`test/run_standard_batsrus_wave_tests.sh`. SBW01–05 use the user-selected
standard BATSRUS producer contract: Wave1/Wave2 contain total Alfvén-wave energy
densities in the along-B/opposite-B branches, with magnetic fraction one half.
The typed view checks conversion, orientation/order, complete endpoint identity,
interpolation numerical limits and atomic rejection. It does not derive four
circular polarizations from two integrated values. SBWN01/02, group
`standard-batsrus-waves-native`, exercise controlled imported native wave buffers
and read-only ownership/rejection. They do not demonstrate a live SWMF producer,
complete imported velocity-gradient/composition contracts or full production
closure. The background fingerprint uses the adapter's single contract identity;
the Python tooling neither duplicates nor selects a different normalization.

```sh
python3 -B srcSEP/test/run_tests.py --suite material-tube-kinematics --output-dir /tmp/sep-material-source --plot none
python3 -B srcSEP/test/run_tests.py --amps ./amps --group material-tube --output-dir /tmp/sep-material-native --plot none
python3 -B srcSEP/test/run_tests.py --suite standard-batsrus-waves --output-dir /tmp/sep-batsrus-source --plot none
python3 -B srcSEP/test/run_tests.py --amps ./amps --group standard-batsrus-waves-native --mpi-np 2 --output-dir /tmp/sep-batsrus-native --plot none
python3 -B tools/sep_test_orchestrator.py --app srcsep --skip-build --only material-tube-kinematics --only material-tube-native --only standard-batsrus-waves --only standard-batsrus-waves-native --output-dir /tmp/sep-material-batsrus-gates
```

All four orchestrator steps require complete callback coverage. A source/native
SKIP or missing MTG/SBW/SBWN case becomes ERROR for its explicitly requested
gate, retaining the child's original report. Native steps use two MPI ranks and
owned cwd/output directories. Negative reader tests leave an initial SKIP and
omit the remaining cases to verify that partial or unexercised runs cannot
appear green. These tooling fixtures prove failure accounting; the actual
mathematical/native outcomes are recorded against their tested build separately.
The build-8 six-case material baseline predates the later extreme-shear
flux-precision negative control. That uncovered condition is being corrected
and revalidated within MTG06; the earlier six PASS verdicts do not qualify it.

## Actual host-step birth and analytic source-epoch gates

`analytic-shock-epoch` selects ASE01–05 and the strict/sanitized source runner
`test/run_analytic_shock_epoch_tests.sh`. Independent constant-speed and
linear-speed exponential references check radius and speed at an explicit
departure [s], arbitrary query order, atomic invalid inputs, and mixed-view
provenance rejection. ASE05 rejects unprepared/foreign/stale owners, including
different future knot tables with exactly coincident departure radius and
speed. The helper evaluates the existing trajectory law; it
neither reads PIC's mutable clock nor advances a native owner. All 13 compiled
source/header/launcher inputs are hashed before and after compilation and
execution. A source change discards the verdicts. These five value tests do
not establish that a real PIC timestep dispatches a source or moves a newborn.

HOST01–03 (`host-step-native`) address that separate integration seam, once per
production mover, using the actual `amps_time_step` → `PIC::TimeStep` → source
and mover callbacks. Each ID requires a fresh native process and exactly one
explicit input protocol with repeat one. The ordinary catalog invocation
without the protocol returns SKIP; that SKIP cannot qualify a requested host
gate. The dedicated runner supplies
`test/fixtures/host_step_smoke.json` through the versioned
`SEP_NATIVE_HOST_STEP 1` protocol. Its numerical/source controls are a declared
software diagnostic, not a selected physical thermal, guiding-center, source
frame, or observational model. The fixture uses the retained zero-origin
analytic shock law and explicitly constant source energy/pitch controls,
prescribed frozen waves and disabled wave coupling/evolution.

The native callbacks observe two actual host steps without setting or manually
advancing PIC's clock. They require pending step-start birth records, their
full first movement exactly once, matching initial Cartesian X and field-line
coordinate S, one actual clock increment per step, and completed analytic shock
geometry/ledger endpoints. Public fixture controls are restored; PIC private
first-call/statistics state is contained by process termination. The runner
also proves the existing multiple-input-ID rejection and the additive
exclusive-process repeat rejection before mesh initialization.

Each rank publishes immutable cumulative `host-step-trace-000000.json`,
`host-step-trace-000001.json`, … snapshots and normal transactional manifests.
Trace schema `srcsep-native-host-step-trace-v2` carries the sequence, compiled
counter policy, requested/completed step counts, a final-complete marker, and
the actual source/mover event records. The reader requires contiguous names,
unchanged cumulative prefixes and a complete last snapshot; early files remain
available after a fatal mover exit. It verifies every rank and first-birth
identity coverage. A real native assertion FAIL with exit one remains FAIL;
missing traces, SKIP, changed executable/input, contradictory exits and harness
errors become ERROR. Both outcomes return nonzero from the campaign.
Both the JSON fixture and the actual generated native `fixture.args` are
hashed; the strict reader also requires the latter to be the exact versioned
serialization of the former. A changed native input cannot hide behind an
unchanged higher-level fixture.

```sh
make -C srcSEP test-analytic-shock-epoch-unit
python3 -B srcSEP/test/run_tests.py --suite analytic-shock-epoch --output-dir /tmp/sep-source-epoch --plot none
python3 -B srcSEP/test/run_native_host_step_tests.py --amps ./amps --fixture srcSEP/test/fixtures/host_step_smoke.json --mpi-ranks 2 --output /tmp/sep-real-host-steps
python3 -B tools/sep_test_orchestrator.py --app srcsep --skip-build --only analytic-shock-epoch --only host-step-native --output-dir /tmp/sep-host-gates
```

The archived build-10b test-only native run reached all three real host
callbacks, completed two steps and eight first-birth moves per case, and failed
the initial X/S and completed shock-epoch assertions. It also exposed the first
fixture's repeated immutable-output name; the sequenced trace revision corrects
that harness defect without changing the physical assertions. No remote-AMR
block/timestep failure occurred in this run, so no core timestep repair is
claimed or applied. Current native results are recorded against their exact
staged binary in [`STATUS.md`](../STATUS.md); an ASE or harness PASS does not
replace the native gate. The tested native counter policy is ON; an OFF build
requires its own evidence. No full production/source-frame/thermal/material/
polarized closure or source-driven observation comparison is qualified here.

The `coefficient-radius-native` orchestrator step separately invokes
`test/run_native_coefficient_radius_tests.py`, default RADIUSN01. This uses the
tagged schema-6 native startup protocol and the reviewed complete input
template, with declared translated geometry/nondefault wind edits retained in
the manifest. It checks the actual read-only QLT1 provider against an independent
translated-radius law, distinguishes the former global-zero radius, and rejects
an imported radius with no declared origin. No host timestep/source or inferred
SWMF coordinate frame is included. Its strict reader requires native zero exit,
the exact requested PASS with zero assertions, and unchanged executable,
template and owned input hashes. `--prepare-only` is PREPARED and never PASS.
The separate `coefficient-radius-library-native` step requests RADIUSN02 through
`--case RADIUSN02` and the reviewed parallel-diffusion input template. It checks
the actual configured radius/speed-dependent library route with its declared
power-law inputs, using the provider's translated heliocentric position. Its
reader pins that selected ID and its own template hash; a RADIUSN01 result
cannot satisfy the library gate. Neither test invents an origin for imported
data whose coordinate frame is unknown.

```sh
python3 -B srcSEP/test/run_native_coefficient_radius_tests.py --amps ./amps --mpi-ranks 2 --output /tmp/sep-radius-native
python3 -B srcSEP/test/run_native_coefficient_radius_tests.py --amps ./amps --case RADIUSN02 --mpi-ranks 2 --output /tmp/sep-radius-library-native
python3 -B tools/sep_test_orchestrator.py --app srcsep --skip-build --only coefficient-radius-native --only coefficient-radius-library-native --output-dir /tmp/sep-radius-gate
```

## Bounded finite-amplitude response study checks

The source-suite alias and orchestrator step `finite-amplitude-study` invoke
`test-finite-amplitude-study-unit` and
`test/run_finite_amplitude_study_tests.sh`. This launches the 16 verbose Python
unit tests in `test_finite_amplitude_response.py`, with before/after hashes of
the listed solver/protocol/plotter/input/test closure. It does not launch the
broad study or choose physical limits. The checks cover independent exact gyro
motion/order, analytic finite-time Lorentz forcing and helicity signs, weak/full
orbit comparison, DOP853 and energy/momentum invariants, separate phase/mode
quadratures, finite-gap semantics, numerical underflow/undefined pitch, explicit
resource failures, raw-coordinate/provenance binding and portable pinned plots.

The result reader uses the actual verbose-unittest class/method IDs and
retains individual failure blocks. All 16 named checks are required; missing or
SKIP prerequisites become ERROR. A source-only success is bounded mathematical
and reproducibility evidence, not a production error budget, kinetic dispersion
qualification, adopted gap bridge or full coupled/native convergence. The
experiment's exact physical/numerical inputs and remaining domain limits are
documented in
[`FINITE_AMPLITUDE_RESPONSE_STUDY.md`](../validation/physics/FINITE_AMPLITUDE_RESPONSE_STUDY.md).

```sh
make -C srcSEP test-finite-amplitude-study-unit
python3 -B srcSEP/test/run_tests.py --suite finite-amplitude-study --output-dir /tmp/sep-finite-bounded --plot none
python3 -B tools/sep_test_orchestrator.py --app srcsep --skip-build --only finite-amplitude-study --output-dir /tmp/sep-finite-gate
```

CP03G08–10 extend `bash srcSEP/test/run_background_derivative_tests.sh` with
strict warnings, ASan and UBSan: both focused kernels and Parker follow exact
plasma-comoving ballistic characteristics for every actual provider declaration
(analytic, SWCME, reduced-front, SWMF and local evolution); material-label changes preserve the
fraction under translation/stretching; unknown frames reject; and physical
isotropic compression remains active in all three kernels. Fixed-grid controls
retain bulk convection only in the shared helper outside the srcSEP provider
contract. Unknown provider values reject. These are production-kernel manufactured tests; the
native coupled run and restart/coupling regressions remain separate gates.

## Standalone core integration, stages 0–4

`test-standalone-integration-unit` is a prerequisite of `test` and
`test-controlled-analytical`. It runs SINT01–07, SINT07R and SINT10–16 through
production core/adapter/geometry/transport code in an owned temporary source
snapshot. It emits individual JSON/JUnit results, exact commands and input
hashes; strict warnings are errors. `--sanitizers` adds ASan/UBSan to the new
adapter and numerical paths. Fixed statistical seeds, sample counts and
acceptance bounds are in the source and retained logs. LeakSanitizer is not
claimed in a ptrace-managed environment.

`test-standalone-integration-native` requires a clean matching configured
field-line archive selected with `AMPS_FIELD_LINE_NATIVE_BUILD`. Its actual
MPI 1-/4-rank fixture exercises SINT08/09, SINT13N and SINT14N–16N. Every
requested rank count/mode must finish. Startup errors, missing cases, stale
objects and compile failures produce ERROR. No MPI dependency substitute is
accepted. Native testing uses manufactured initial ions and a conditioned
MFP no-event path; existing FTED/FTEM suites independently test scattering,
hazard distributions and conservation. It introduces no production source.

```sh
python3 srcSEP/test/run_tests.py --suite standalone-integration \
  --output-dir /tmp/sep-pure-FRESH --plot none
python3 srcSEP/test/run_tests.py --suite standalone-integration-native \
  --amps-source /path/to/matching/configured/AMPS \
  --output-dir /tmp/sep-native-FRESH --plot none
```

Both suites retain per-case compatibility reports for the existing source
runner, while command/hash files live in owned child directories. The shared
registry names are `CORE-SEP-STANDALONE-INTEGRATION-001` and
`CORE-SEP-STANDALONE-NATIVE-001`; the latter requires the declared
`native-field-line` capability. The application orchestrator steps are
`standalone-core-integration` and `standalone-core-native`, with readers that
require the complete stable ID list. A missing/partial child report cannot
become PASS. Optional capability SKIPs are never native qualification.

The bounded evidence and known full-volume prerequisite timeout are recorded
in [STANDALONE_INTEGRATION_STATUS.md](../validation/STANDALONE_INTEGRATION_STATUS.md).
Moving restart and AMPS–BATS-R-US coupling remain their separate later-stage
gates; historical runs are not attributed to this source snapshot.
## Moving material checkpoint verification (roadmap stage 7)

The prescribed reduced-core consumer uses checkpoint schema 5. The tests
exercise production material ownership, particle histories, source ledgers,
collective output and segment exposure. They do not promote a frozen-line
analytical fixture or a live-wave subsystem to production qualification.

| ID | Required native check |
| --- | --- |
| MCP01–03 | Parker/Dmumu/MFP: independent continuous, write and fresh-resume jobs on MPI 1, 2, 4 and 8 |
| MCP04–06 | The same movers with MPI ownership changed 1→8 and 8→1 |
| MCP07 | Corrupt-file collective rejection on MPI 1 and 8 with no publication/counter mutation |
| MCP08 | Actual production population-controller reduction of a mixed cooled/fresh source, preserving immutable birth energy as a seventh moment |
| MCP09–11 | Each mover resumes a later save with changed segment count, then advances through the next output |
| MCP12 | Last-rank rejection after a valid nonempty detached replacement is prepared; no live pointer/state/counter changes, then an ordinary retry succeeds |

MCP08 is in `run_native_standalone_integration.py`; the other IDs form
`run_material_checkpoint_tests.py`'s complete 58-job campaign. Every resume
requires nonempty, unique particle identities. Same-decomposition final images
are byte-exact. Changed-MPI comparisons retain every history byte and exact
exposure; only six-decimal native text serialization has its justified 2e-6
comparison bound. SI number, birth-energy and signed current-energy balances
retain their separate 1e-10 bound. No numerical tolerance is enlarged after a
failed run. Unmeasured segments and measured vacuum remain distinct.

Use a freshly configured, clean `test/sep_reduced_front__field_line` native
build, which has the core source's physical proton mass and registered mover
dispatch. The runner audits all PIC-header consumers; selected movers being
fresh is insufficient when another object can supply an old inline pool method.

```sh
make -C srcSEP test-material-checkpoint-unit
make -C srcSEP test-material-checkpoint-native \
  NATIVE_BUILD_ROOT=/path/to/configured/AMPS \
  MATERIAL_CHECKPOINT_OUTPUT=/path/to/FRESH/material-evidence
```

For a site that launches MPI externally, preserve the prepared commands and
record **actual** return codes in each `process-result.json`; missing/crashed
jobs cannot pass. Preparation alone is not execution evidence:

```sh
python3 srcSEP/test/run_native_standalone_integration.py \
  --build-root /path/to/configured/AMPS --ranks 1 2 4 8 \
  --output-dir /path/to/FRESH/native --prepare-only
python3 srcSEP/test/run_material_checkpoint_tests.py \
  --native-prepared /path/to/FRESH/native/prepared.json \
  --output /path/to/FRESH/material --prepare-only
# Execute prepared jobs in dependency order; prepare the isolated CRC control
# after its writer succeeds, using --prepare-corruption. Then collect:
python3 srcSEP/test/run_material_checkpoint_tests.py \
  --output /path/to/FRESH/material --collect
```

The shared test registry exposes `CORE-SEP-MATERIAL-CHECKPOINT-001`; the
application orchestrator exposes `material-checkpoint-native`. Both require
the complete mover/role/decomposition catalog and import individual failures.
`test-material-checkpoint-unit` is a routine `test` prerequisite. Its isolated
negative controls detect changed birth/work accounts, missing/duplicate
histories, wrong/global sampling denominators, missing/deadline/failed process
receipts and changed source identities; four dependency controls detect stale
headers and cross-build PIC authority. These controls are software checks,
not substitute native evidence.

Full-driver CLI restart is additionally checked after a complete physical
step: save at 45 and resume to the same intended 64 steps, with fixed dt 300 s.
The 64-step two-/four-thread worker checks do not complete the separate
400-step production output gate. Current executed results, source snapshots,
and the repaired legacy-fixture initialization are recorded in
[the status report](../validation/STANDALONE_INTEGRATION_STATUS.md).

### Complete material, operator and threaded continuation gate

`test-stage7-checkpoint-native` runs a fixed catalog of **44 real MPI jobs**.
CP05C01–05 are five two-rank material routines; CP05F01–03 have independent
two-rank write/resume roles. Both use the real host step, explicit reduced
event/startup deck, fixed dt 300 s, a nonempty accelerated source population,
and prescribed scattering. MFP continuation must observe a physical scattering
event. The existing malformed-checkpoint, hazard, Cartesian and flux-owner
assertions remain mandatory.

CP05O01–05 and CP05OF01–03 retain the old wave/codec assertions in a separate
native executable. They advance isolated operators and serialization records,
never particles or the host on a frozen nonzero-U line. They qualify codec and
operator state only. No production material-owner guard is bypassed. MCP13
requires a bare material case to exit 2 with a native ERROR report rather than
abort without output.

Material checkpoint records are serialized in immutable-particle-ID order.
Worker arrival and native allocation order carry no physical information in
this keyed model. Every complete normalized history record remains exact;
geometry, samples and ledgers are still compared in full. CP05C01–03 reverse
actual segment lists and require an unchanged complete image, restoring both
link directions before subsequent host steps. MFP requires at least one
multi-record segment. The fixed-line shared-stream operator codec retains its
semantic list order. This is canonical serialization, not a tolerance or a
partial-record comparison.

TCR01–03 run each production mover continuously and through a separate-process
save at step 45/resume to step 64, with two and four actual worker threads.
Dmumu additionally runs on eight MPI ranks with four workers per rank. Final
checkpoints must match byte for byte and every fixed-step number, birth-energy
and signed current-energy ledger must satisfy the unchanged 1e-10 bound.
These short restart tests supplement the separate 400-step output campaign.

```sh
AMPS_FIELD_LINE_NATIVE_BUILD=/path/to/base \
AMPS_FIELD_LINE_THREAD2_BUILD=/path/to/thread2 \
AMPS_FIELD_LINE_THREAD4_BUILD=/path/to/thread4 \
make -C srcSEP test-stage7-checkpoint-native \
  STAGE7_CHECKPOINT_OUTPUT=/path/to/FRESH/stage7
```

The `physics-checkpoint-native` Python **suite** alias selects this maintained
gate. The obsolete bare-group orchestrator launch is replaced by
`stage7-checkpoint-native`, preserving every CP05C/F assertion. The shared
registry exposes all 20 logical IDs through `CORE-SEP-STAGE7-CHECKPOINT-001`.
The collector requires the complete role/thread/rank catalog, actual process
receipts, matching source/helper/archive/executable hashes and all native
reports. Missing evidence is ERROR; an evaluated numerical mismatch is FAIL.
`--prepare-only` and `--collect` support externally launched MPI, using each
job's explicit `argv`, `cwd` and dependency list. Preparation is not a pass.

Material fixture requests use the existing test-input channel, with a strict
`SEP_CHECKPOINT_FRESH_PROCESS 2` header followed by quoted `role`, `prefix`,
`startup_input` and `event` values in that order. Paths are absolute. Routine
cases use role `routine`; fresh cases use `write`/`resume`. Operator-only roles
retain protocol 1. The helper supplies the reviewed verification clock and
source/coefficient settings before normal allocation; production defaults do
not change. The material request requires actual reduced ownership, and its
results report material epochs and physical dt.

MCP12's callback receives only staged object counts. It may reject preparation
but cannot publish or modify a replacement. It is rejected in production run
purpose, has no CLI/environment hook, and must not call MPI because another
rank can fail before reaching it. The normal collective agreement discards all
detached replacements before returning the failure. The test checks complete
state equality, caller counters and actual live vertex/segment pointers, then
loads the same valid file normally.

# Reduced source and return-boundary verification

The maintained standalone core suite includes SINT17–SINT24 for source
configuration, the DSA number Jacobian, actual-epoch offset intersections,
proton/tube flux, outward pitch support and analytic first-passage laws.
`make test-standalone-integration-unit` retains strict warning and structured
failure propagation. SINT23 includes independent stationary and moving-front
first-passage CDFs. SINT26 checks the exact rotation-axis background shortcut
against off-axis dipole/wind identities while preserving every primitive
validity check. SINT27 grades the core's exact parallel RH bracket selection
against independent gas-dynamic identities, including transverse-flow controls
that retain the general scan. No reference regeneration or relaxed tolerance is used.
SINT28 adds exact closed-endpoint MFP clock checks over 4096 fixed-key numerical
durations/early scatters, plus absorbed-clock and no-later-noise controls.
SINT29 checks 64 constant-flow trajectories with short final numerical
intervals against independent endpoints at the unchanged 1-mm tolerance.
It is a positive stress check, not a reproduction of the long MFP failure.
SINT30 retains the exact iteration-84 one-/four-rank trajectory failure and
compares the corrected endpoint against an independent, converged midpoint
oracle. The original position tolerance, maximum numerical step and work
bound remain fixed; only the internal endpoint partition changes.

SINT25N–SINT32N extend the linked native fixture to actual net-source births,
all three movers, signed energy/number ledgers, labelled manufactured monoenergy,
zero efficiency, atomic rejection and finite-front disconnection. The fixture
also rejects a low-energy band with no outward escape: dropping those proposals
would silently change the specified net rate. It
checks actual output metadata against the core apex at a non-cadence epoch,
and checks serialized nanotesla against the native tesla-valued record. It
joins every native destination list after direct dispatch and performs the
production MPI exchange before closing the source window. A separate complete application run
is required to qualify the PIC::TimeStep driver and output sampling sequence.

The native runner supports `--prepare-only` to compile/hash and list launch
jobs when an environment restricts MPI inside Python subprocesses. Preparation
produces no test PASS. Execute every listed real MPI job, retain each process's
actual return code in `process-result.json`, then invoke the same runner/output/
rank selection with `--collect`. Collection checks unchanged source/binary
hashes and complete native results. Missing logs, return codes or case IDs are
ERROR. Ordinary environments should use the default direct runner workflow.

The complete-driver stage-6 gate also has a maintained evidence checker:

```sh
make -C srcSEP test-standalone-integration-application EVIDENCE_DIR=/path/to/application-campaign
```

The campaign directory must contain `launch-plan.json`, `build-sha256.json`,
the immutable executable, and one `parker/`, `fte-dmumu/`, `fte-mfp/` directory
with actual `process-result.json` receipts and `application.rank0.log` output.
SAPP01–03 require all 400 physical 300-s steps, the complete sampling output,
segment lifetimes reconstructed independently from topology commits, active
background metadata, and number/birth/current-energy balances recomputed from
printed measures at the original 1e-10 relative bound. It requires identical
source histories across movers. Missing output, crashes or modified executable/
inputs remain ERROR/FAIL; partial campaigns cannot satisfy the gate. The checker
writes `verified-results.json` and `junit.xml` beside the retained launch evidence.
The same Make target runs eleven negative controls against isolated copies of the
completed output. Wrong physical clocks/cadence, a lost step, broken birth
energy, altered segment dt, inactive shock metadata, nonfinite active fields,
changed executable/source identities and missing/failed process receipts must all be
rejected. Those deliberate corruptions are audit tests, not new simulations.
The gate verifies every recorded source/header/archive hash as well as the
executable and physical inputs; old output cannot certify a changed checkout.

These tests quantify software/numerical behavior for the declared synthetic
event. A manufactured release distance is not a production SEP calibration,
and passing an unsupported-mode rejection is not downstream qualification.


### Shared material contract: stages 10–11

MBC01–20 test exact profile admission, immutable SI/HCI identity and coverage,
explicit species/EOS and all nine derivative entries, canonical samples and
trajectories, local shock limits, descriptor round-trip/rejection and failed
publication rollback. The independent radial trajectory check integrates the
actual accelerating Parker wind; it does not assume constant solar-wind speed.
The metadata controls distinguish an absent geometry generation from the valid
initial value zero; the actual native fixture also checks geometry registration
ownership before permitting publication, without renumbering PIC geometry.
No additional backend, downstream volume or partial physical step is admitted.

The neutral material trajectory has the same Parker-exterior region contract
as a material state query. MBC18 rejects actual open/closed PFSS points, a
zero-duration PFSS query, a verified reverse-time crossing from the Parker
exterior into the PFSS corona, and a foreign adapter descriptor. The wrapper
checks every canonical primitive used by the existing RK controller and its
accepted endpoint; it does not infer that endpoint membership proves the
whole path qualified. Its positive trajectory remains bitwise equal to the
legacy diagnostic mover and satisfies the unchanged independent flight-time
oracle. The generic diagnostic `MoveVertex` API remains available under its
existing broader core-query contract.

MBC02 also exercises the detached Stage 12 SWCME1D interval against its actual
core event. It checks the thermodynamic basis, Cartesian derivatives, parcel
map, local fast-shock limits, radial/time coverage, identity, and rejection of
native binding or a sheath-volume query. The independent physics oracle lives
with the SWCME owner. Its ten cases are individually selectable as
`MAT1D01`–`MAT1D10` and are included alongside the four R1 cases in the existing
bounded SWCME routine. The global extended registry uses IDs
`SWCME-MAT1D01`–`SWCME-MAT1D10`.

```sh
python3 src/models/swcme/test/run_tests.py --group MATERIAL1D --output-dir /tmp/sep-stage12-upstream
```

```sh
make -C srcSEP test-material-background-unit
python3 -B srcSEP/test/run_tests.py --suite material-background --plot none --output-dir /tmp/sep-material-suite
python3 -B srcSEP/test/run_material_background_tests.py --sanitizers --output-dir /tmp/sep-material-asan
python3 -B tools/sep_test_orchestrator.py --app srcsep --skip-build --only material-background --output-dir /tmp/sep-material-orchestrated
python3 -B src/test/run_tests.py --config test/qualification-machine/passion.json --test CORE-SEP-MATERIAL-BACKGROUND-001 --output-dir /tmp/sep-material-shared
```

Use fresh output directories. The maintained runner builds actual canonical
libraries in isolation, enforces `-Werror`, hashes compiled inputs and rejects
changed inputs. It records actual assertion counts, per-case durations, reference residuals,
unchanged tolerances and source configuration in JSON/JUnit; missing build/output is
ERROR for every required case and exits nonzero. Existing evidence is immutable.
The orchestrator self-tests also exercise malformed/missing catalogs, retained
FAIL/ERROR statuses, actual compiler failure and output-directory reuse.
This suite is a prerequisite of routine `test`/`test-controlled-analytical` and
belongs to the shared **extended** campaign. It is not selected by smoke alone.
Native MPI, full output windows and restart/thread gates remain separately
required; source assertions never substitute for them. See the
[qualification record](../validation/STAGE11_MATERIAL_BACKGROUND_CONTRACT.md).

CV01 provenance hashes the canonical shared transport implementation at
`src/models/sep_common/sep_transport_common.cpp`, alongside the actual focused
transport kernel. A missing source remains an error; the former application-local
path cannot silently stand in for that relocated production implementation.
The relocation gate likewise requires the complete six-root canonical include
contract in both source and copied build layouts. Its expectation includes the
shared transport, SWCME, both CME libraries, shared front and generic input
owners; extra or omitted roots fail the exact comparison.

SCOPE02 distinguishes the authoritative one-dimensional field-line coordinate
from PIC's derived Cartesian position cache. `CommitAndAttach` must project the
exact committed coordinate on the particle's own immutable line and commit that
SI position to the same particle; it must not add a separate solar-wind
displacement or perform an independent Cartesian advance. Motion of the line's
vertices already enters this view through the geometry projection, so adding
`U*dt` again would double-count that motion. `check_field_line_position_contract.py` permits
exactly that one projection/commit across the seven audited advance units and
rejects additional position accesses, pointer escapes, changed coordinates or
owners, extra writers and geometry mutation. Its fifteen in-memory negative
controls run as part of `make -C srcSEP test-field-line-scope-unit`; the existing
Boris, perpendicular-diffusion and AMR-list bans remain enforced. These are source
contract controls, not numerical physics evidence. Native SINT08 independently
checks the actual stored Cartesian view against `GetCartesian` after movement
and remains required on the declared MPI decompositions.

### Isolated Routine initialization and material fixture ownership

`make test`, `make test-native-amps-validation`, `make test-mpi`, and the
maintained Python `--routine` entry use `routine_fixtures.py`. The native
`amps --all-tests` CLI retains its existing selection and startup semantics.
The maintained dispatcher reads the actual compiled catalog and requires all
223 reviewed Routine IDs, including their initialization, seed/repeat,
isolation and tag metadata. Each case receives its declared startup protocol;
CV01–CV12 and IV01–IV06 execute their complete independent comparisons. No
Routine ID is removed to avoid a prerequisite. Missing fixtures, changed
metadata, duplicate/foreign reports, contradictory child exits and incomplete
role groups are ERROR. Native SKIP remains visible and makes this required
aggregate unsuccessful. Extended cases remain separate explicit campaigns.
The existing `--all` orchestrator surface reuses these Routine results and
continues through **every Extended descriptor** in its existing independent
workflow; it is never replaced with `--routine`. Selected CP05C01–CP05C05 reuse
the same reviewed hashed source/event protocol through the existing Python
`--test` and Make `test-case` interfaces. An explicit native role manifest still
selects the raw native path.

Examples (the linked executable is reused):

```sh
python3 test/run_tests.py --amps ../amps --routine --plot none
make test JSON=/absolute/results.json JUNIT=/absolute/results.xml
make test-case CASE=BIRTHN01
python3 test/run_tests.py --amps ../amps --test BIRTHN01 --plot none
python3 test/run_tests.py --amps ../amps --group birth-window-native --mpi-np 2 --plot none
```

`make test-case CASE=BIRTHN01` and the maintained single-ID selector execute
that ID's complete reviewed operator/control roles. An explicitly forwarded
native `--test-input` retains the raw role-protocol path. Birth roles default
to two MPI ranks, including the real BIRTHN05 migration check. Other native
cases retain the requested/default rank count. Every repeated Birth experiment
uses fresh complete role processes and the existing shared case/replicate seed
law. The context remains bound through registry seed-report stamping, so actual
streams and reported case seeds agree. Physical startup dt is unchanged within
every process. A raw repeated role or multiselection is an explicit error.

JSON/JUNIT Make destinations and Python `--report-json`/`--report-junit` export
the same authoritative aggregate. Repeating an entry retains fresh nested raw
case directories instead of replacing earlier native evidence. The parsed
fixture declaration is saved and hash-bound throughout execution; a changed
manifest blocks dependent callbacks and invalidates every required aggregate
record while preserving raw child reports. Fixture arguments cannot override
the selected callback, hashed input or report destinations. Software boundary
and negative protocol tests exercise these failure paths without pretending
their labelled fake reports are model calculations. They are prerequisites of
`test-python-runner-unit` and therefore of `make test`.

`NativeStorage` preserves the existing initialization enum values and appends a
separate allocation capability. It means native buffers exist; it does not
claim a production particle snapshot, a material map, or the absence of
particle allocation inside a diagnostic. SOURCE02, for example, initializes a
real thermal population under its own tested contract. PNW's explicit
snapshot-only tag requests actual SWCME primitive export and spectral storage.
The orchestrator selects `--shock-model swcme1d` explicitly. Wave inertia uses
actual exported composition/mass density and magnetic/plasma data, rather
than the unconsumed energetic SEP source mass/charge. Production and all other
source hosts retain their original species checks and tolerances. This is
snapshot/export evidence, not admission of SWCME1D material transport.

FTE01 and PARKER01 explicitly select analytical startup and install their
existing manufactured fields before their own production particle preflight.
All three components of U vanish at every centerline vertex, so the exact
centerline material map is the identity. PARKER01 still supplies its full
nonzero transverse velocity-gradient tensor: zero centerline U does not imply
zero divergence or zero adiabatic momentum change. Scoped detached native
lines preserve and restore the exact original owner, pointers, buffers,
geometry and sampling vectors. Negative controls reject nonzero centerline
flow and a foreign registered owner. Provider names, scalar U-dot-B and an
ownership flag cannot substitute for this physical declaration.

BIRTHN01–05 are process-owned frozen-departure operator experiments. Every
role initializes the same fixed physical timestep of .01 s. Actual scheduled
births consume the existing .0025-s first-move remainder; the direct control
calls the existing diagnostic operator with that duration, and the full-step
negative control uses .01 s. These are native operator calls, not shorter
physical host steps or fractional sampling. All roles start with the same
curved Parker realization and publish the same completed full-step material
map. BIRTHN06 instead uses the qualified reduced source/checkpoint host with
its original fixed timestep, genuine injection/ledgers and preserved sampling
contract. The diagnostic role's operator duration and physical clock are
reported separately and truthfully; these tests do not introduce new source
physics or dense-time material transport equations.

After real material publication, restoring an obsolete fixed R0 source
namespace or snapshot over R1 would be invalid. The single-process-owned
protocol validates its final source, clock, geometry and immutable snapshot
before it exits. A fresh empty material source namespace is prepared and
independently validated before any diagnostic birth; existing counters or live
particles cannot be relabeled. Other native callbacks retain exact restoration
and failure propagation. The unsupported generic analytical checkpoint guard
remains intact; BIRTHN06 uses the actual supported reduced source owner. No
production guard or numerical acceptance limit is relaxed by these fixture
boundaries.
