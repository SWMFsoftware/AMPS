# srcSEP field-line solar energetic particle transport

`srcSEP` models solar energetic particle (SEP) transport along magnetic field
lines embedded in three-dimensional space.  The field lines carry solar-wind,
interplanetary magnetic-field, shock, and turbulence state from analytic
providers, SWCME, or an SWMF coupling.  The source also contains the
self-consistent Alfvén-turbulence subsystem, including integrated and
wave-number-resolved representations and particle-wave coupling.

This source includes the Step 1 selectable standalone component-test registry,
the Step 2 immutable background/clock boundary, the Step 3 common SI flux-tube
geometry/source normalization, the Step 4 production mover API, the Step 5
field-line-only scope boundary, the Step 6 common transport numerics, the Step
7 canonical Parker solver, the Step 8 coefficient-driven focused-transport
solver, the Step 9 event-driven mean-free-path mover, the Step 10 shared
coefficient registries, the Step 11 authoritative/restartable turbulence
subsystem, and the Step 12 deterministic parallel-reduction contract.
Step 13 adds reportable, fail-closed acceptance fixtures, Step 14 removes
the retired mover implementations/aliases while enforcing a clean source
handoff, and Step 15 adds a claim-separated scientific-validation campaign.
The numbered campaign now includes linked-application cases `CV01`–`CV05`.
They cover ballistic streaming, constant and nonuniform spatial diffusion,
adiabatic cooling, and magnetic focusing through the same native-registry,
reviewed-input, independent-reference, metric, visualization, and provenance
contract.
See [STEP13_CHANGE_MANIFEST.md](STEP13_CHANGE_MANIFEST.md),
[STEP14_CHANGE_MANIFEST.md](STEP14_CHANGE_MANIFEST.md), and
[STEP15_CHANGE_MANIFEST.md](STEP15_CHANGE_MANIFEST.md). Retired names remain
mapped in [MIGRATION_MANIFEST.md](MIGRATION_MANIFEST.md).
Self-consistent Alfvén turbulence remains a production subsystem and uses the
same physical area and volume as injection and particle sampling.

## Field-line-only production scope

`srcSEP` advances particles only by a scalar coordinate on a magnetic field
line, plus parallel/normal momentum state. The line itself remains embedded in
three-dimensional space and continues to carry three-component IMF, plasma,
electric-field, and SWMF/SWCME state. This distinction is enforced at compile
time: the enclosing AMPS configuration must enable field-line mode and attach
particles to field-line segments.

Step 5 removed the Parker3D, 2019 He, Kartavykh, Borovikov, drift, and default
Boris dispatch implementations; their macros, Cartesian particle offsets,
cell-data derivative support, and mesh-cell particle-list paths are also gone.
The public mover set remains exactly `parker`, `fte-dmumu`, and `fte-mfp`.

At the user's direction, the spatial-neighborhood sampling module was removed
now. Field-line observer sampling remains available and reports density, flux,
return flux, pitch angle, energy, Larmor radius, and mean free path through
`SEP::Sampling::InitSingleFieldLineSampling`. Outputs that require collecting
particles in a three-dimensional volume around an arbitrary Cartesian point
are temporarily unavailable in srcSEP and should be hosted by the separate 3-D
application. See [STEP5_CHANGE_MANIFEST.md](STEP5_CHANGE_MANIFEST.md) for the
complete transfer boundary.

## Flux-tube geometry and source units

`SEP::FieldLine::FluxTubeGeometry` returns area in m² and volume in m³. With
valid magnetic data it enforces `A|B| = constant`; otherwise the application
must configure an explicit area profile. Segment population, solar-wind and
shock sources, wave-energy density, growth rates, turbulence transport, and
sampling all consume this interface.

Shock-source efficiency is configured once through `InjectionEfficiency` for
analytic, SWCME, and SWMF providers. The historical hard-coded particle-weight
override has been removed, and MeV input energies are converted to joules before
SI momentum routines. See [FLUX_TUBE_GEOMETRY.md](FLUX_TUBE_GEOMETRY.md) for the
API, fallback policy, normalization contract, and focused-test definitions.

## Production mover API

The supported public mover set is now exactly `parker`, `fte-dmumu`, and
`fte-mfp`. `--list-movers` reports their field-line representation, turbulence,
and coefficient capabilities without initializing AMPS. A single adapter
validates particle and field-line attachment state before dispatch. Main-loop
physics queries mover capabilities and no longer compares function addresses.

Step 14 closes the legacy-name transition period: only the three names printed
by `--list-movers` are accepted. Ambiguous movers, direct-wave experimental
movers, former aliases, and 3-D trajectory names all fail before model
initialization. Startup metadata prints the canonical choice, coefficient
authority, and active `Dxx`, `Dmumu`, and mean-free-path providers. See
[PRODUCTION_MOVER_API.md](PRODUCTION_MOVER_API.md) for the runtime contract and
[MIGRATION_MANIFEST.md](MIGRATION_MANIFEST.md) for every retired replacement.

## Common and canonical transport numerics

All public mover shells now share one validated particle representation, one
physical-distance field-line advance, one attachment path, one
`d ln|B|/ds` focusing convention, exact plasma-frame adiabatic cooling, and
explicitly keyed random streams. All three canonical paths use strict
relativistic conversion and named composable timestep limits;
invalid input and a stability-limit underflow are errors rather than silent
clamps.

The canonical Parker mover advances the Itô process
`ds=(U_parallel+d(kappa)/ds)dt+sqrt(2*kappa*dt)dW`, with `kappa` and its SI
gradient supplied together by `SpatialDiffusionProvider`. The canonical
`fte-dmumu` mover obtains `Dmumu`, `dDmumu/dmu`, provenance, and turbulence-state
identity from `PitchAngleDiffusionProvider`; it uses symmetric deterministic /
stochastic / deterministic splitting, reflective pitch-angle boundaries, and
midpoint streaming. The canonical `fte-mfp` mover samples exact exponential
waiting times with event rate `nu=v/lambda`, treats `lambda=+infinity` as the
ballistic limit, applies focusing/cooling between events, and redistributes
pitch angle isotropically in the selected Alfvén-wave frame with exact Lorentz
velocity transforms. The QLT Kolmogorov spectrum is normalized so its integral
is `deltaB^2` and `Dmumu` is returned in s^-1.

Neither mover mutates shared wave arrays from a particle worker. Coupling
records are self-contained (stable particle ID, species, statistical weight,
snapshot generation, pre/post momentum, interval, event, branch, and boundary
metadata) and are published once each in-domain interval completes. Boundary
exit records are clipped to the exact endpoint before particle deletion. The
post-`PIC::TimeStep()` update rejects duplicate physical keys and sorts without
using particle-buffer addresses, thread IDs, or ranks. See
[TRANSPORT_NUMERICS.md](TRANSPORT_NUMERICS.md) and the Step 6–10 change manifests
for equations, contracts, tests, and limitations.

WP01–WP10 additionally make `SEP::Turbulence::PICAdapter::Advance` the sole
reachable production turbulence mutation entry for standalone and coupled
library stepping. The adapter converts legacy shock/coupling sources inside one
transaction, maps authoritative PIC segment state into the common core, runs
the configured operator order, exports the authoritative representation, and
derives output fields. See
[WP01_WP10_IMPLEMENTATION.md](WP01_WP10_IMPLEMENTATION.md) for the work-package
mapping and validation boundary.

WP11–WP20 add bounded Milstein pitch-angle diffusion, named local error
budgets, source-bound coefficient inputs, repaired constant/Jokipii/Florinskiy
providers, adaptive spatial-diffusion quadrature, species-aware transport and
injection, named turbulence scales, and typed ballistic compatibility. The
implementation, equations, configuration behavior, per-package verification,
and native-review boundary are recorded in
[WP11_WP20_IMPLEMENTATION.md](WP11_WP20_IMPLEMENTATION.md).
The executed source-only evidence and explicit native/external-data boundary are
summarized in
[WP11_WP20_VALIDATION_REPORT.md](WP11_WP20_VALIDATION_REPORT.md).

WP21–WP30 repair analytical shock evolution/restart, general field-line shock
intersections, injection-spectrum normalization and semantic RNG keys, shock
turbulence ledgering, per-line magnetic flux, strict particle diagnostics,
unit-bearing sampling products, transactional output, and immutable run control.
See [WP21_WP30_IMPLEMENTATION.md](WP21_WP30_IMPLEMENTATION.md) for the formulas,
configuration/provenance contracts, production wiring, and limitations, and
[WP21_WP30_VALIDATION_REPORT.md](WP21_WP30_VALIDATION_REPORT.md) for executed
evidence versus blocked native/external gates.

WP31–WP41 add shared turbulence stiffness planning and typed exceptional paths,
configured particle-population invariants and lineage, a fail-closed native
adapter observation contract, a generated 90-case compatibility matrix, global
system checkpoint ledgers, versioned seed panels, bounded property/fault tests,
an instrument forward model, deterministic performance counters, and evidence-
level governance. See
[WP31_WP41_IMPLEMENTATION.md](WP31_WP41_IMPLEMENTATION.md) for each package's
why/what/how contract and
[WP31_WP41_VALIDATION_REPORT.md](WP31_WP41_VALIDATION_REPORT.md) for executed
versus blocked evidence.

All three movers obtain coefficients through one registry. Canonical CLI names
select coefficient authority (`prescribed`, `self-consistent`, `swmf`), spatial
closure (`from-dmumu`, `from-mfp`), pitch-angle provider (`configured`,
`constant`, `jokipii-1966`, `florinskiy`), MFP
model (`qlt`, `qlt1`, `tenishev-2005`, `chen-2024`, `from-spatial`), and invalid
value policy (`fail`, `ballistic`). Conversion cycles and source/ownership
mismatches are rejected rather than inferred.

The WP11–WP20 registry also exposes resonance-gap and turbulence-amplitude
policies and makes spectrum, correlation, quadrature, and mover-error scales
named configuration. Coupled sources require a source-bound provider and never
fall back to the legacy configured callback. Parker configurations that could
produce infinite spatial diffusion fail during preflight, while `fte-mfp`
retains the exact ballistic zero-event-rate state.

The controlled mover cases `PARK01`–`PARK07`, `FTED01`–`FTED08`, and
`FTEM01`–`FTEM08` are now descriptors in the same selectable component-test
registry used by the linked CLI. The source-only sanitizer targets execute
those exact callbacks. Analytical references include convection, diffusion
moments, cooling, Itô drift, waiting/event statistics, ballistic transport,
first passage, persistent-random-flight transport, automatically calculated
refinement order, and Legendre-mode pitch-angle decay. See
[CONTROLLED_COMPONENT_VALIDATION.md](CONTROLLED_COMPONENT_VALIDATION.md) for
equations, tolerances, and commands.

## Authoritative turbulence and reproducibility

Turbulence source is explicit and independent of the mover: `prescribed`,
`self-consistent-integrated`, `self-consistent-spectral`, `swmf-read-only`, or
`swmf-initial-then-local`. Exactly one integrated or spectral state is
authoritative. A fixed driver order applies sources, particle exchange,
CFL-subcycled advection, reflection, cascade/dissipation, synchronization,
diagnostics, and checkpoint accounting. Prescribed and SWMF-read-only states
are immutable; the SWMF-to-local path requires a one-time epoch/checksum
handoff.

Every update produces a signed energy ledger. Boundary policies, remap,
spectral grid, pending coupling, RNG metadata, operator phase, and accumulated
ledger are restartable. Worker contributions use physical keys and one
canonical reduction order, while random streams use campaign, particle, step,
and purpose keys. See [TURBULENCE_MODEL.md](TURBULENCE_MODEL.md),
[STEP11_CHANGE_MANIFEST.md](STEP11_CHANGE_MANIFEST.md), and
[STEP12_CHANGE_MANIFEST.md](STEP12_CHANGE_MANIFEST.md). The exact source-only
results and native-build boundary are recorded in
[STEPS10_12_VALIDATION_REPORT.md](STEPS10_12_VALIDATION_REPORT.md).

WP31 declares the production time integrator as first-order Lie splitting and
chooses one stage count from advection, source, reflection, and cascade limits.
`--turbulence-operator-safety`, `--turbulence-max-source-fraction`,
`--turbulence-max-cascade-fraction`, `--turbulence-min-substep`, and
`--turbulence-max-substeps` configure and fingerprint this policy. Corrections,
rejected source energy, physical zeros, and per-operator work are explicit.

The bounded WP01--WP10 contract gate is available as
`make test-wp01-wp10-unit` (or `./test/run_wp01_wp10_tests.sh`). It compiles
only dependency-free numerical contracts; the production PIC adapter still
requires the enclosing AMPS configuration.

The corresponding WP11--WP20 gate is available as
`make test-wp11-wp20-unit` (or `./test/run_wp11_wp20_tests.sh`). It executes ten
focused physics/numerics tests with ASan/UBSan and verifies production source
wiring without claiming a native AMPS result.

`TURB02`–`TURB23` and `TURBOWN01` are also registered with the common CLI.
The new `TURB21` test advects a nonuniform periodic sine profile and measures
L1 convergence and energy conservation against exact translated cell averages.
`TURB22` compares the complete wave-energy growth history with an analytically
integrated linear-in-time source. It validates application and accounting of a
known growth source; it deliberately does not claim to validate the separate
QLT calculation that derives that source from particle distributions. These
cases strengthen the earlier uniform-state and one-step invariants without
replacing the AMPS/PIC/MPI-dependent coupled growth-rate validation.
`TURB23` independently calculates the kinetic-energy change of controlled
macro-particles scattered in both resonant Alfvén-wave branches, applies the
opposite signed wave-energy increments, and requires particle-plus-wave total
energy and the turbulence ledger to close to `2e-12 J`.

## Background state and time

PIC owns the single authoritative simulation clock. The standalone global time
and separate SWCME launch counter have been removed; SWCME state is now prepared
at the exact PIC epoch before the particle step instead of one iteration late.
Shock motion and output timestamps read the same clock through
`SEP::Background::SimulationTimeSeconds()`.

Every call to a configured particle mover occurs inside one immutable background
read phase. Its const snapshot declares the provider, ownership, epoch, validity
interval, field-line generation, configuration fingerprint, and provenance.
Analytic/SWCME updates are model-owned, SWMF imports are read-only, and a
provider change is rejected unless SWMF data have first been copied and an
explicit local-evolution handoff is published. Background publication is
rejected while `PIC::TimeStep()` is moving particles.

See [BACKGROUND_STATE.md](BACKGROUND_STATE.md) for the data contract, update
order, API, ownership rules, units, tests, and current limitations.

## Step 15 scientific validation

Step 15 adds four stable validation cases after the numerical architecture has
stabilized. `VAL01` compares Parker diffusion/cooling with independent
analytical solutions. `VAL02` compares the two focused movers under a matched
mean-free-path closure using pitch moments, diffusion, onset, peak, and
fluence. `VAL03` compares the Dmumu ensemble with a separately implemented
finite-volume solver under an identical coefficient history. `VAL04-SWCME`
passes real SWCME adapter states through the immutable background boundary and
production Parker kernel and also requires an active SWCME shock source.

The release report keeps numerical, cross-mover, cross-model, coupled, and
observational evidence separate. A complete coupled result additionally
requires a checksum-verified real SWMF replay. Observational readiness requires
held-out spacecraft products, uncertainty metadata, instrument forward
operators, and onset, anisotropy, spectra, fluence, decay, and
multi-spacecraft-longitude metrics. Missing external evidence is
`INCOMPLETE`; manufactured data cannot be labelled observational. See
[validation/README.md](validation/README.md) and
[STEP15_VALIDATION_REPORT.md](STEP15_VALIDATION_REPORT.md).

The evidence classes have intentionally separate commands:

```sh
make test-controlled-analytical
make test-native-amps-validation SEP_EXECUTABLE=/path/to/amps
make test-swmf-validation SWMF_MANIFEST=/evidence/swmf-replay.json
make test-observational-validation \
  OBSERVATIONAL_MANIFESTS="/evidence/event-1.json /evidence/event-2.json"
```

A successful controlled or native command cannot satisfy either external-data
gate. SWMF and observational targets require checksum-verified `PASS`
manifests; the observational target also requires the union of supplied
held-out events to cover every declared metric family.

## Standalone component-test CLI

The registry extends the existing parser in `util/sep_cli.cpp`; there is no
second parser.  These commands are available only through the standalone
executable:

```sh
../amps --list-tests
../amps --test DXX01
../amps --test=DXX01 --test=TURB01
../amps --test-group parker
../amps --test-group=turbulence
../amps --test-group fte-dmumu
../amps --test-group fte-mfp
../amps --test PARK07 --test FTED08 --test FTEM08
../amps --test TURB21 --test TURB22 --test TURB23
../amps --all-tests
../amps --all-tests --test-json results.json --test-junit results.xml
```

IDs and group names are case-insensitive.  Repeated/overlapping selectors run a
test once, in stable ID order.  `--all-tests` includes only tests classified as
`routine`; expensive `extended` tests remain visible and can be selected by ID
or group.  `--list-tests` exits before post-compile input, AMPS, field lines,
turbulence, SWCME, or output initialization.

Every execution selector is test-only: the driver resolves it before expensive
initialization, initializes the maximum prerequisite declared by the selected
tests, prints one MPI-root summary, and exits before `TestManager()` and the
production timestep loop.  `FAIL` returns status 1, `ERROR` returns status 2,
and selections containing only `PASS`/`SKIP` return 0.  A skipped test remains
labelled `SKIP`; it is never presented as a successful scientific assertion.
Any positive/non-finite `assertion_failures` metric converts a nominal callback
PASS into FAIL. Requested JSON/JUnit output is acceptance evidence: inability to
create either report is ERROR/status 2, and both formats retain IDs, seeds,
configuration, metrics, messages, durations, and artifact paths.

### Python campaign runner and analytical figures

`test/run_tests.py` is the single orchestration interface for selecting one,
several, a group, the bounded routine set, or every registered test. It calls
the production registry rather than duplicating test selection or physics in
Python. The default figure formats are PNG and EPS:

```sh
# Run Test 01 / CV01 in the linked application, then compare its output with
# the independent characteristic and create PNG/EPS figures.
python3 test/run_tests.py --amps /path/to/amps --validation-case CV01 \
  --output-dir /absolute/path/to/evidence/CV01

# Run CV02-CV05 together; each still invokes the selected linked application.
python3 test/run_tests.py --amps /path/to/amps \
  --validation-case CV02 --validation-case CV03 \
  --validation-case CV04 --validation-case CV05 \
  --output-dir /absolute/path/to/evidence/CV02-CV05

# Discover stable IDs from the linked executable.
python3 test/run_tests.py --amps /path/to/amps --list

# Run any chosen cases; --test and --group may each be repeated.
python3 test/run_tests.py --amps /path/to/amps \
  --test PARK07 --test FTED08 --group turbulence \
  --output-dir /absolute/path/to/evidence/selected

# Run the bounded native set, or literally every discoverable test including
# extended Monte Carlo cases. The latter can be substantially more expensive.
python3 test/run_tests.py --amps /path/to/amps --routine \
  --output-dir /absolute/path/to/evidence/routine
python3 test/run_tests.py --amps /path/to/amps --all \
  --output-dir /absolute/path/to/evidence/all

# Run the dependency-light analytical mover/turbulence suite without AMPS.
python3 test/run_tests.py --suite controlled-analytical \
  --output-dir /absolute/path/to/evidence/controlled

# Re-render a previously retained registry report without rerunning physics.
python3 test/run_tests.py --from-json results.json \
  --output-dir /absolute/path/to/evidence/figures
```

For an analytical case, a reported CSV artifact with conventional coordinate,
`numerical`/`model`, and `analytical`/`exact` columns becomes a pointwise
solution overlay. When a C++ test reports only an error, moment, convergence
order, or conservation residual, the runner instead plots that numerical
metric against its analytical reference or acceptance limit and labels the
figure as metric-level evidence. It never reconstructs an expected solution
with production model code. `run_manifest.json`, the authoritative registry
JSON/JUnit reports, `analytical_plot_manifest.json`, the command log, and
`plots/*.{png,eps}` remain together in the selected output directory.
CV01 first verifies that the selected linked executable advertises the case,
then runs all six numerical realizations through `amps --test CV01`. It writes
native JSON/JUnit, raw model and independent reference states for periodic/open
boundaries at all three timesteps, particle and packet-moment tables, a
four-panel diagnostic, the resolved input, and checksummed executable
provenance. A missing or stale binary is an error, never a Python/standalone
fallback. See
[validation/cases/CV01/README.md](validation/cases/CV01/README.md).
CV02-CV05 use the same fail-closed application lifecycle and add exact
Green-function, independent finite-volume, adiabatic-characteristic, and
focusing-characteristic references. See the case README files under
[validation/cases](validation/cases/README.md).

The source-only Step 13 fixtures registered in the same production catalog are
`BG01` (analytic/SWCME provider epochs), `BG02` (mock read-only SWMF import and
handoff), `CROSS01` (matched ballistic focused movers), and `CROSS02` (the
isotropic `kappa`/`Dmumu`/`lambda` closure). Bounded fixtures run under
`--all-tests`; long Monte Carlo campaigns remain explicitly selectable.

The historical `--test-manager`, `--testmanager`, `--run-test-manager`, and
`--no-test-manager` switches retain their earlier behavior: when enabled they
run the legacy diagnostics after field-line initialization and then continue
into production.  They are deliberately not aliases for the new test-only
interface.

## Production and coupled behavior

A run with no arguments retains the existing defaults and enters the normal
production loop.  Existing option aliases, Boolean spellings, case-insensitive
values, `--option value`/`--option=value` forms, and last-occurrence-wins
semantics are preserved.  CLI values continue to be applied after the optional
post-compile input file.

The registry is called only from `main.cpp`.  SWMF-facing library entry points
do not parse process arguments, list tests, or launch standalone tests.  This is
important for a coupled execution in which SWMF owns initialization and run
control.

## Building and testing

The enclosing AMPS checkout supplies `Makefile.conf`, PIC/field-line headers,
MPI, and the final linked executable.  From `srcSEP`, the default expected path
is `../amps`; override it when necessary:

```sh
make test-cli-unit
make test-state-unit
make test-geometry-source-unit
make test-mover-api-unit
make test-field-line-scope-unit
make test-transport-common-unit
make test-parker-unit
make test-fte-dmumu-unit
make test-fte-mfp-unit
make test-coefficients-unit
make test-turbulence-core-unit
make test-reproducibility-unit
make test-wp21-wp30-unit
make test-wp31-wp41-unit
make test-python-runner-unit
make print-configuration-matrix
make test-acceptance-unit
make test-documentation-unit
make test-scientific-validation
make test-controlled-analytical
make test-sanitizer
make test-stochastic-repeat
../amps --list-movers
make test-list SEP_EXECUTABLE=/path/to/amps
make test-case CASE=DXX01 SEP_EXECUTABLE=/path/to/amps
make test-group GROUP=turbulence SEP_EXECUTABLE=/path/to/amps
make test-turbulence SEP_EXECUTABLE=/path/to/amps
make test-stress SEP_EXECUTABLE=/path/to/amps
make test-mpi MPI_NP=4 SEP_EXECUTABLE=/path/to/amps
make -j test SEP_EXECUTABLE=/path/to/amps JSON=results.json JUNIT=results.xml
```

The focused unit targets are dependency-light checks of the parser,
immutable background core, SI geometry/source core, production mover registry,
field-line-only source boundary, common transport kernels, Parker solver, and
coefficient-driven focused-transport solver, event-driven MFP solver, and
coefficient registry, authoritative turbulence driver, reproducible reduction,
reportable acceptance fixtures, and the cleaned public/source inventory. The
Step 6–14 and WP31–WP41 numerical targets use strict C++11 warnings plus AddressSanitizer and
UndefinedBehaviorSanitizer. The Step 15 numerical runner uses the same C++11
checks; its real SWCME replay uses C++17 because that is SWCME's public API
baseline. The remaining targets intentionally invoke the linked production CLI
and propagate its status.

The Python plotting test requires Python 3 and Matplotlib. It uses only a
temporary synthetic structured report, verifies PNG and EPS output paths, and
does not require the linked executable. See [test/README.md](test/README.md)
for all runner selectors, output semantics, MPI launch support, and exit codes.

For the native Step 5 completion gate, build this application in a field-line
AMPS configuration, then build the receiving Cartesian-transport application
from its own clean tree. That second application is not included in this
archive. A source-only run can establish `SCOPE01`–`SCOPE03`, but cannot claim
either linked build as complete. See [test/README.md](test/README.md) for the
catalog, prerequisites, isolation contracts, and extension procedure.

## Current limitations

The registry and focused numerical suites are not a claim of complete model or
observational validation. Step 15 records their evidence as separate classes
and intentionally reports the delivered campaign `INCOMPLETE`. A source-only archive cannot
exercise the enclosing AMPS/PIC adapter, MPI rank decomposition, coupled SWMF
epochs, or long campaign conservation behavior; those remain native integration
gates. Boundary-exit particle-to-wave flux is intentionally not deposited
because the legacy coupling callback still requires a live particle record.
The dependency-light TURB01–TURB20 and PAR01–PAR05 suites cover the new core and
synthetic decomposition invariance, but do not replace native AMPS/OpenMP/MPI,
SWMF handoff, long-campaign conservation, or observational validation.
`SCAT01` is explicitly `SKIP` because the historical stochastic diagnostic has
no quantitative acceptance rule and includes a singular zero-energy sample.
The historical routine remains callable through `--run-test-manager` until its
physics reference and valid input domain are approved in a later step.

The archive intentionally contains no build products, coverage data, runtime
output, test reports, or nested archives. `mover.cpp`, `fte_mover.cpp`, and
`fte_mover_dmumu.cpp` no longer exist; the three canonical adapters and
`mover_state.cpp` are the complete mover source layout. A source-only archive
still cannot claim the native AMPS strict-warning, MPI, or coupled-SWMF gates;
commands for those checks are recorded in the migration manifest.

The WP31–WP41 source gate establishes only analytical-core and source-integration
contracts. Native mover traversal, full seam conservation, long multi-seed
ensembles, native fuzzing, real instrument comparisons, and MPI/OpenMP scaling
remain BLOCKED until their named executable, data, or hardware baselines are
supplied. Documentation uses these evidence levels: `analytical-core`,
`source-integration`, `native-amps`, `swmf-replay`, and
`observational-validation`.
