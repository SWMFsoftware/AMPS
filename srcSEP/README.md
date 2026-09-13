# srcSEP field-line solar energetic particle transport

`srcSEP` models solar energetic particle (SEP) transport along magnetic field
lines embedded in three-dimensional space.  The field lines carry solar-wind,
interplanetary magnetic-field, shock, and turbulence state from analytic
providers, SWCME, or an SWMF coupling.  The source also contains the
self-consistent Alfvén-turbulence subsystem, including integrated and
wave-number-resolved representations and particle-wave coupling.

This source includes the Step 1 selectable standalone component-test registry,
the Step 2 immutable background/clock boundary, and the Step 3 common SI
flux-tube geometry and source normalization. Step 3 does not remove the
self-consistent Alfvén-turbulence subsystem; it makes that subsystem use the
same physical area and volume as injection and particle sampling.

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

Strictly equivalent legacy CLI names are accepted for one transition release
with warnings; ambiguous movers, direct-wave experimental movers, and 3-D
trajectory names are rejected. Startup metadata prints the canonical choice and
active `Dxx`, `Dmumu`, or mean-free-path provider. See
[PRODUCTION_MOVER_API.md](PRODUCTION_MOVER_API.md) for mappings, aliases,
rejections, capabilities, and examples.

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
../amps --all-tests
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
../amps --list-movers
make test-list SEP_EXECUTABLE=/path/to/amps
make test-case CASE=DXX01 SEP_EXECUTABLE=/path/to/amps
make test-group GROUP=turbulence SEP_EXECUTABLE=/path/to/amps
make test-turbulence SEP_EXECUTABLE=/path/to/amps
make -j test SEP_EXECUTABLE=/path/to/amps
```

The first three unit targets are dependency-light strict-warning checks of the
production parser, immutable background core, and SI geometry/source core. The
remaining targets intentionally invoke the
linked production CLI and propagate its status.  See [test/README.md](test/README.md)
for the catalog, prerequisites, isolation contracts, and extension procedure.

## Current limitations

The registry is test infrastructure, not a claim of complete model validation.
`TURB01` checks one deterministic wave-energy closure; it does not replace the
planned Alfvén-turbulence initialization, advection, cascade, reflection,
particle-coupling, conservation-ledger, MPI, restart, and observational tests.
`SCAT01` is explicitly `SKIP` because the historical stochastic diagnostic has
no quantitative acceptance rule and includes a singular zero-energy sample.
The historical routine remains callable through `--run-test-manager` until its
physics reference and valid input domain are approved in a later step.
