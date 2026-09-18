# srcSEP standalone component tests

## Run the complete test series

From the `srcSEP` directory, run every test advertised by the linked AMPS
application with:

```sh
test/run_tests.py --amps ../amps --all --output-dir test_output/all
```

The runner discovers the complete native registry, executes every test in an
isolated process, and automatically uses the registered input/reference
workflow for CV, IV, and XM validation cases. A failed, errored, or crashed
test does not prevent later tests from running. The final output reports the
PASS/FAIL/SKIP/ERROR totals and lists every failed or errored test with its
diagnostic. Detailed JSON, JUnit, logs, plots, and per-test artifacts are saved
under `test_output/all`.

Step 1 provides one catalog and result contract for standalone component tests.
The catalog lives in `component_tests.cpp`; generic deterministic selection,
metadata validation, timing, status handling, and output formatting live in
`util/sep_test_registry.*`.  The production CLI remains `util/sep_cli.*`.

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
../amps --all-tests
../amps --all-tests --test-json results.json --test-junit results.xml
```

`make -j test` intentionally sequences shared-state component execution even
when Make is given `-j`, then runs the canonical
`src/models/swcme/test/run_tests.py --routine` suite. SWCME is not embedded in
srcSEP. Expensive extended tests are not part of `--all-tests`.

The relocation gate can also be selected through the Python runner:

```sh
test/run_tests.py --suite swcme-relocation \
  --output-dir test_output/swcme-relocation
```

In the normal AMPS layout it discovers `../src/models/swcme` automatically. A
detached source checkout may set `SWCME_DIR=/absolute/path/to/swcme`; this is an
explicit development override, never a fallback to a directory within srcSEP.

The Step 6–15, WP11–WP20, WP21–WP30, and WP31–WP41 focused targets do not require the linked executable, AMPS, PIC,
MPI, SWMF, or field-line host classes. Each compiles the exact production
numerical cores with C++11, `-Wall -Wextra -Werror -pedantic`, AddressSanitizer,
and UndefinedBehaviorSanitizer. LeakSanitizer alone is disabled because the
managed test environment does not expose the required `/proc` task data.

## Python runner and analytical-comparison plots

`run_tests.py` provides one Python 3 entry point for native registry tests and
dependency-light Make suites. It does not contain test callbacks, movers, or
analytical physics: stable IDs/groups are discovered from the linked executable
and execution status comes from the registry JSON. Matplotlib is the only
additional runtime dependency for figures.

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

`--routine` forwards the native `--all-tests` policy and therefore excludes
extended cases. `--all` first calls `--list-tests`, accepts only identifier
tokens containing a digit (so the printed `ID | ...` heading can never become
a selection), and runs every discovered ID in a separate process. Ordinary
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
cleanup; sanitizer executables and objects are never retained. Other source
suites that do not emit structured reports still run normally and produce a
log/manifest but cannot generate analytical plots. Verify the orchestration and
both figure paths with:

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
| `DXX01` | `diffusion` | routine | field-line model | `GetDxx` agrees with the constant-coefficient analytical result and the existing million-panel independent quadrature at relative tolerance `1e-5`. | Temporarily replaces and restores the pitch-angle diffusion function pointer; no artifact. |
| `FTE01` | `transport` | routine | field-line model | The focused-transport mover follows the expected field-line displacement while preserving velocity in a static-plasma fixture, using the existing `1e-2`/`1e-5` checks. | Uses fixed registry seed 1002, restores vertex data and diffusion pointer, deletes its particle, and clears test lists. |
| `PARKER01` | `parker` | routine | field-line model | Parker convection keeps the line coordinate stationary in the fixture and matches the analytical density-driven momentum update within `1e-5`. | Uses fixed registry seed 1001, restores vertex data and diffusion pointer, deletes its particle, and clears lists. |
| `PARKER02` | `parker` | extended | field-line model | Four million legacy stochastic trials produce at least one in-range displacement sample and rank 0 successfully writes the histogram.  This is an execution/output assertion, not yet a Gaussian-shape validation. | Uses fixed registry seed 1003; rank 0 writes `dxParker.dat`; particle/lists are cleaned. |
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
branches, and duplicate physical-key rejection. The same script checks that
both production drivers call the single turbulence PIC adapter and that all
three movers publish self-contained coupling records without a stored particle
pointer.

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

- `CLI01`: help/list mode selection;
- `CLI02`: both `--test` forms and exactly-once execution;
- `CLI03`: deterministic ordering, de-duplication, groups, routine/extended;
- `CLI04`: malformed/conflicting/unknown selection rejection;
- `CLI05`: result exit propagation and frozen legacy TestManager parsing;
- `HIDDEN`: positive `assertion_failures` cannot be masked by callback PASS;
- `REPORT`: JSON/JUnit schema, outcome, metric, and artifact preservation;
- registry completeness: required metadata, callbacks, and unique IDs.

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
