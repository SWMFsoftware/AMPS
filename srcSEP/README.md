# srcSEP field-line solar energetic particle transport

`srcSEP` models solar energetic particle (SEP) transport along magnetic field
lines embedded in three-dimensional space.  The field lines carry solar-wind,
interplanetary magnetic-field, shock, and turbulence state from analytic
providers, SWCME, or an SWMF coupling.  The source also contains the
self-consistent Alfvén-turbulence subsystem, including integrated and
wave-number-resolved representations and particle-wave coupling.

## B01 source-distribution baseline

`srcSEP` is distributed as source, documentation, reviewed validation inputs,
and test scripts only. Its authoritative classification is
[`SOURCE_MANIFEST.json`](SOURCE_MANIFEST.json). The manifest is enforced by the
AMPS-level package gate rather than by either SEP application, preserving the
rule that `srcSEP` and `srcSEP3D` are independent applications.

Run the focused gate with:

```sh
python3 test/run_tests.py --suite package-hygiene \
  --output-dir test_output/package-hygiene
```

The gate rejects retired mover/SWCME sources, native binaries, object and
archive files, dependency files, Python caches, test-output trees, and AMPS
runtime binary data. It also runs an isolated negative control that inserts a
deliberate stale object and verifies that the checker fails. Local ignore rules
name only generated products; validation JSON/CSV inputs remain visible.

The removed monolithic mover, full-three-dimensional drift/sample sources, and
application-local SWCME examples were not production objects. Their active
responsibilities are owned respectively by the three canonical field-line
movers, field-line sampling/output, `srcSEP3D`, and `src/models/swcme`. See
`MIGRATION_MANIFEST.md` for symbol-level mappings.

## B02 canonical SEP-common dependency

The seven provider-neutral kernels have exactly one source and binary owner:
`AMPS/src/models/sep_common`.  No corresponding `.h` or `.cpp` file remains
under `srcSEP/util`. Public srcSEP headers route the canonical names through
`util/sep_common_header_path.h`. In an installed AMPS tree it prefers
`src/models/sep_common/<header>` through the standard `-I$(AMPS_ROOT)`, which
also prevents a stale short-name header from winning by include order. A
detached component build falls back to the public short name supplied by
`-I$(SEP_COMMON_DIR)`. This is a path resolver, not a forwarding copy; both
branches open the same canonical file and the canonical archive remains the
only binary owner.

The fallback is necessary because `sep.h` is included transitively through
`pic.h` while AMPS compiles libraries such as meshAMR. Variables appended by
`build/main/makefile` apply to that child build only and cannot modify a sibling
submake's already-generated compile command. The global AMPS-root include is
part of those commands, so the resolver makes source `srcSEP/sep.h` and copied
`build/main/sep.h` behave identically without modifying `Makefile.conf` or
installing duplicate shared headers.

The same constraint applies to production translation units compiled by the
generic AMPS object rule. Files such as `field_line.cpp`, `diffusion.cpp`, the
private SWCME adapter, and the linked validation cases therefore include
`sep_common_header_path.h` and request each canonical header through
`SRCSEP_SEP_COMMON_HEADER(...)`. A bare `#include "sep_species_source.h"` can
work in a focused local build yet fail after the application is copied to
`build/main`, because the earlier child make cannot add its private include
variables to that later command. Test-only translation units retain short
includes intentionally and receive `-I$(SEP_COMMON_DIR)` explicitly.

`srcSEP/makefile` resolves `SEP_COMMON_DIR` from `AMPS_ROOT`, builds
`sep_common.a`, and inserts the seven canonical objects into `mainlib.a`.
Inserting objects preserves AMPS's existing application-archive link contract;
the archive audit then requires each member exactly once and rejects duplicate
strong symbols.  The same active-makefile algorithm works both at
`AMPS/srcSEP/makefile` and after AMPS copies the application to
`AMPS/build/main/makefile`.  An override is supported only for a deliberately
detached component checkout.

The focused ownership gate is:

```sh
make test-sep-common-ownership-unit
```

It rejects local copies and ad hoc path-qualified shared includes, verifies all
canonical sources, rebuilds and audits exact archive membership, links a real
consumer through that archive, and reproduces both installed makefile layouts
from an unrelated working directory. It also rejects bare shared-header
includes in production `.cpp` files, compiles public and production-source
header surfaces with only `-I$(AMPS_ROOT)`, and repeats the public-header probe
from a synthetic copied `build/main` layout. This directly covers both the
meshAMR/PIC include context and the `field_line.cpp` failure mode.
Dependency-light physics tests compile the same canonical sources with their
sanitizer/debug flags; they do not keep test-only copies.

## Canonical SWCME dependency

SWCME is not stored or built beneath this application. Its single source and
binary owner is `AMPS/src/models/swcme`, shared independently by `srcSEP` and
`srcSEP3D`. The application makefile follows the proven srcSEP3D pattern:

- it derives `AMPS_ROOT` from the active makefile rather than the shell's
  current directory, so the same file works in `AMPS/srcSEP` and after AMPS
  copies it to `AMPS/build/main`;
- it defines `SWCME_DIR=$(AMPS_ROOT)/src/models/swcme` and gives the private
  adapter and validation-registry objects that canonical include root;
- `sep.h` exposes no SWCME header or type. Only
  `adapters/swcme1d_adapter.cpp` includes `swcme1d.hpp`, while its application
  API remains provider-neutral;
- `util/sep_swcme_validation.cpp` privately includes `swcme1d_input.hpp`
  because the native D01-D03 callbacks validate the canonical production
  configuration rather than a copied test schema;
- both SWCME-facing objects have explicit make rules. Older AMPS
  `Makefile.conf` generic recipes do not consistently consume flags appended
  to `CPPFLAGS`, `CXXFLAGS`, or `INCLUDE`, so the explicit recipes pass
  `-I$(SWCME_DIR)` in both the source tree and copied `build/main` tree;
- it builds the canonical `swcme.a` and inserts its one production member,
  `swcme3d.o`, exactly once into `mainlib.a`; and
- source-only validation uses the same `SWCME_DIR` and never searches for an
  application-local fallback.

`SWCME_DIR` can be overridden for a deliberately detached checkout. No
override is needed in the installed AMPS layout. The former
`swcme----moved-out` tree and root-level SWCME demo are absent; the maintained
provider sources, demonstrations, tests, and documentation all live under
`src/models/swcme`.

The focused relocation gate is:

```sh
make test-swcme-relocation-unit
```

It reproduces both source and copied-build layouts, compiles both the
provider-neutral adapter API and the canonical-header implementation, and
rejects application-local SWCME sources, wrappers, public-header copies,
implementation namespaces, or path-qualified includes. It also rebuilds and
audits the one-member canonical archive before compiling the private adapter.
`make strict-production` delegates to the enclosing `make amps` build and then
audits `AMPS/build/main/mainlib.a` to require exactly one member for every
SEP-common kernel and exactly one `swcme3d.o`, with no duplicate strong symbol.

### D01 fail-closed SWCME background handling

The 1-D adapter no longer clips every sub-domain radius or converts malformed
density, speed, and divergence to zero. `QueryAtRadius` returns a typed,
transactional result containing one immutable SI sample. Production defaults to
`--swcme-failure-policy strict`; `clamp-radius` and `diagnostic-fallback` are
explicit, fingerprinted experimental policies. Every recovery is counted and
reported at shutdown. Shock injection treats an unrecovered query as a fatal
background-consistency error and prints rank, epoch, radius, failed field, and
source-state ID instead of silently injecting no particles.

The fallback sample can be set with `--swcme-fallback-density` `[m^-3]`,
`--swcme-fallback-speed` `[m/s]`, and `--swcme-fallback-divergence` `[s^-1]`.
All values are validated and stored in the frozen run/restart fingerprint even
when strict mode leaves them dormant. Detailed lifecycle and physics behavior
is documented in [BACKGROUND_STATE.md](BACKGROUND_STATE.md). Run the focused
gate with `make test-swcme-fail-closed-unit`. The identical callback is linked
into the native C++ registry as extended test `D01`, so
`test/run_tests.py --amps ../amps --all ...` executes it automatically in an
isolated process.

### D02 general SWCME and shock configuration

`src/models/swcme/swcme1d_input.hpp` now owns the one textual configuration
schema. It expands `fast`/`slow`, applies input-file and CLI/programmatic layers,
converts explicitly named units into canonical public units, calls the existing
model/source validators, and emits a deterministic normalized manifest plus
fingerprint. srcSEP transports assignments through its private adapter and does
not duplicate canonical model fields.

Standalone input uses a `SWCME1D on` block; see
[`PARAM.SWCME1D.example`](PARAM.SWCME1D.example). Repeatable
`--swcme-override 'key=value unit'` values have the final authority. The schema
covers ambient wind, Parker settings, ballistic/DBM/data-driven CME kinematics,
shock acceleration representation, region geometry/smoothing, launch and
validity epochs, source species, energy bounds, normalization, and injection
efficiency. The resolved energy interval and efficiency drive the existing 1-D
injector; relative shock source weight multiplies its swept-volume source. AMPS
species mass/charge are checked against the resolved source after AMPS setup.

Invalid or duplicate input is rejected before mesh construction with key,
authority layer, and source line. Equivalent effective configurations produce
the same field-order-canonical fingerprint regardless of assignment order. A
coupled host uses the parser-free `ConfigurationRequest` API. Focused coverage
is `make test-swcme-configuration-unit`.

The identical configuration callback is linked into the native C++ registry as
extended test `D02`. Both the focused Make target and `amps --test D02` use
`util/sep_swcme_validation.cpp`; there is no second expected-value
implementation that can drift away from the production-registry gate.

### D03 configured AMPS, MPI, refresh, and restart integration

D03 closes the source-only evidence boundary with a configured enclosing AMPS
build and short native campaigns. Every production mover must run in serial and
in at least two distinct MPI decompositions, refresh the SWCME state, complete
at least one validated particle dispatch, and produce a canonical artifact that
is byte-identical across decompositions. A checkpointed/resumed trajectory must
also match its uninterrupted trajectory exactly.

The driver reduces the process-local mover-dispatch and D01 query/recovery
counters after the final timestep. It independently checks that all ranks have
the same SWCME state ID and epoch before printing `mpi_consensus=pass`. The
Python gate archives the build and launch commands, tool versions, executable,
configuration and log checksums, run/SWCME fingerprints, artifacts, and exact
comparison groups in one JSON evidence record.

A missing native tree or reviewed campaign is an explicit `SKIP` during
development and a nonzero `INCOMPLETE` result in `--release` mode. Release mode
requires `--rebuild` and cannot reuse an existing binary. See
[NATIVE_INTEGRATION.md](NATIVE_INTEGRATION.md) for the manifest schema,
site-placeholder boundary, algorithms, evidence interpretation, and complete
commands.

Extended registry test `D03PRE` is a bounded prerequisite check executed by
Python `--all`: it verifies the exact three linked mover names/capabilities and
two monotonic SWCME refresh generations with a canonical fingerprint. It is
deliberately named as a preflight and cannot satisfy the D03 release gate; only
the external campaign can launch multiple decompositions, compare canonical
artifacts, and establish checkpoint/resume equivalence.

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
The numbered campaign now includes linked-application cases `CV01`–`CV12` and
integrated manufactured cases `IV01`–`IV06`.
Cross-model cases `XM01`–`XM03` add an independent focused-transport PDE
comparison, a publication-informed M-FLAMPA transport benchmark, and a linked
Earth-observation comparison for the 2013 April 11 event. XM03 now generates
its Parker-transport spectrum inside the selected executable and compares it
with ACE/EPAM, GOES-13/EPEAD, and SOHO/ERNE data from Liu et al. Figure 12; it
does not require or accept an external model export.
They cover ballistic streaming; constant/nonuniform spatial diffusion;
adiabatic cooling; magnetic focusing; Legendre pitch diffusion; telegraph and
first-passage transport; planar DSA; turbulence advection and growth; and
closed particle-wave exchange. IV01-IV06 then couple Parker geometry,
scattering/focusing, manufactured multi-coordinate transport, moving grids,
moving shocks, and nonlinear wave feedback through the same native-registry,
reviewed-input, independent-reference, metric, visualization, and provenance
contract.
See [STEP13_CHANGE_MANIFEST.md](STEP13_CHANGE_MANIFEST.md),
[STEP14_CHANGE_MANIFEST.md](STEP14_CHANGE_MANIFEST.md), and
[STEP15_CHANGE_MANIFEST.md](STEP15_CHANGE_MANIFEST.md). The CV06–CV12 design,
scope, evidence contract, and commands are summarized in
[CV06_CV12_VALIDATION_IMPLEMENTATION.md](CV06_CV12_VALIDATION_IMPLEMENTATION.md).
The integrated geometry/operator, moving-grid/shock, and nonlinear-feedback
cases are documented in
[IV01_IV06_VALIDATION_IMPLEMENTATION.md](IV01_IV06_VALIDATION_IMPLEMENTATION.md).
XM equations, publication sources, reference extraction, inputs, and evidence
limits are documented in
[XM01_XM03_VALIDATION_IMPLEMENTATION.md](XM01_XM03_VALIDATION_IMPLEMENTATION.md).
Retired names remain
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

B05 classifies every WP42–WP64 item before it can affect the production build.
The authoritative, machine-readable decision record is
[WP42_WP64_DISPOSITION.json](WP42_WP64_DISPOSITION.json): WP42–WP58 and
WP61–WP62 are dependency-light experimental component APIs, while WP59, WP60,
WP63, and WP64 are native/external evidence gates. The four experimental
implementation files compile and run under strict warnings and sanitizers, but
are deliberately absent from `MAINLIBOBJ`; therefore their focused tests are
not evidence that a production mover calls them. Two supporting contracts were
accepted into active code because existing production paths consume them:
`ParkerMeasure` names the phase-space measure used by Parker transport, and an
event-driven wave contribution now carries the actual resonant branch and
pre/post-event particle momenta rather than reconstructing them from limiter
shell endpoints. See
[WP42_WP64_IMPLEMENTATION.md](WP42_WP64_IMPLEMENTATION.md) for algorithms and
promotion criteria and
[WP42_WP64_VALIDATION_REPORT.md](WP42_WP64_VALIDATION_REPORT.md) for the exact
evidence boundary.

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

# Run the advanced controlled portfolio. CV06-CV09 are extended stochastic
# cases; CV10-CV12 are bounded deterministic/conservation cases.
python3 test/run_tests.py --amps /path/to/amps \
  --validation-case CV06 --validation-case CV07 \
  --validation-case CV08 --validation-case CV09 \
  --validation-case CV10 --validation-case CV11 \
  --validation-case CV12 \
  --output-dir /absolute/path/to/evidence/CV06-CV12

# Run all required-nightly integrated manufactured cases.
python3 test/run_tests.py --amps /path/to/amps \
  --validation-case IV01 --validation-case IV02 \
  --validation-case IV03 --validation-case IV04 \
  --validation-case IV05 --validation-case IV06 \
  --output-dir /absolute/path/to/evidence/IV01-IV06

python3 test/run_tests.py --amps /absolute/path/to/amps \
  --validation-case XM01 --validation-case XM02 --validation-case XM03 \
  --output-dir /absolute/path/to/evidence/XM01-XM03

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

# After changing C++ registry sources, rebuild the enclosing application before
# discovery. MAKEFLAGS parallelizes compilation, not registry execution.
env MAKEFLAGS="-j16" python3 test/run_tests.py --amps /path/to/amps --all \
  --amps-source /path/to/AMPS --make-config /path/to/AMPS/Makefile.conf \
  --rebuild --output-dir /absolute/path/to/evidence/all-rebuilt

# Run the dependency-light analytical mover/turbulence suite without AMPS.
python3 test/run_tests.py --suite controlled-analytical \
  --output-dir /absolute/path/to/evidence/controlled

# Re-render a previously retained registry report without rerunning physics.
python3 test/run_tests.py --from-json results.json \
  --output-dir /absolute/path/to/evidence/figures
```

The Python `--all` mode is fault-isolated. It filters the human-readable
`ID | ...` table heading from `--list-tests`, then launches one selected test
per process instead of passing every ID to one AMPS invocation. Each exact
command and normalized `RESULT ID: PASS|FAIL|SKIP|ERROR` line is printed as it
runs. The discovered set includes extended IDs `D01`, `D02`, and `D03PRE`.
A crash, timeout, or missing child report becomes a retained `ERROR` and
does not stop later IDs. CV/IV/XM IDs are run through their registered
validation entrypoints so required case inputs and independent references are
constructed before the linked `amps --test ID --test-input ...` call. The
aggregate JSON/JUnit and final TOTAL/PASS/FAIL/SKIP/ERROR summary cover the
whole discovered portfolio; per-ID evidence is retained under
`individual/<ID>/`.

Make-backed `--suite` selections use the same evidence model. Each receives a
stable `SUITE-*` record containing the exact command, start/end time, duration,
return code, and bounded output excerpt even if compilation fails before a C++
registry can start. Child registry reports are merged with those outer records.
Timeouts, missing launchers, and signals are `ERROR`; ordinary nonzero suite
exits are `FAIL`. A suite may request `SKIP` either by exiting 77 when invoked
directly or by emitting the exact `SRCSEP_SUITE_RESULT=SKIP` marker and exiting
successfully through Make, which does not preserve recipe status 77.
Fail-fast-blocked later suites are also `SKIP`. JSON, JUnit, terminal counts,
`run_manifest.json`, and the process exit code are all derived from the merged records with precedence
`ERROR (2) > FAIL (1) > PASS/SKIP (0)`.

The terminal summary ends with separate `Failed tests (N)` and `Error tests
(N)` lists. Each entry contains the test ID and its normalized diagnostic, so
the cases requiring attention remain visible without searching the complete
run log. A category with no affected tests explicitly prints `none`.

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
CV02-CV12 use the same fail-closed application lifecycle and add exact
Green-function, independent finite-volume, adiabatic-characteristic, and
focusing-characteristic references plus modal, telegraph, inverse-Gaussian,
DSA, manufactured-wave, exponential-rate, and energy-ledger checks. See the case README files under
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
make test-wp42-wp64-experimental
make test-wp59-wp64-native
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
Step 6–14, WP31–WP41, and the explicitly experimental WP42–WP64 component
targets use strict C++11 warnings plus AddressSanitizer and
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

WP42–WP64 remain governed by the B05 disposition record. Passing
`test-wp42-wp64-experimental` proves the isolated algorithms and contracts
compile and behave as specified; it does not promote the four extension source
files into `mainlib.a` and does not establish native traversal. WP59, WP60,
WP63, and WP64 report `SKIP` until a reviewed native AMPS/SWMF, external-data,
or scaling command is supplied through `SRCSEP_NATIVE_GATE`. Promotion
requires an explicit production selector/call path, restart compatibility,
native observations, and a disposition update reviewed with the build manifest.
