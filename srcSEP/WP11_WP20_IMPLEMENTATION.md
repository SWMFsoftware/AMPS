# WP11--WP20 implementation and validation record

This document maps Work Packages 11 through 20 from the physics, numerics, and
testing improvement plan to the delivered source. The packages were implemented
in sequence because later coefficient and mover policies depend on the bounded
stochastic state, error controls, and source-authority contracts established by
the earlier packages.

The dependency-light evidence command is:

```sh
make test-wp11-wp20-unit
```

It compiles the production numerical kernels with C++11, strict warnings,
AddressSanitizer, and UndefinedBehaviorSanitizer, executes controlled analytical
tests, and checks production adapter wiring. It does not replace a linked native
AMPS/PIC/MPI run or validation against real SWMF and observational products.

## WP11: bounded pitch-angle diffusion

### Why

The physical pitch-angle cosine satisfies `-1 <= mu <= 1`. A stochastic update
that evaluates a coefficient beyond this interval, or that clips overshoots,
changes the zero-flux boundary condition and can create artificial endpoint
probability. Euler--Maruyama also omits the multiplicative-noise correction
needed for useful weak accuracy when `Dmumu` varies with `mu`.

### What and how

`util/sep_focused_transport_core.*` now exposes a named bounded scheme. The
default `ReflectingMilstein` increment is

```text
delta_mu = D'(mu) dt + sqrt(2 D(mu)) dW
           + 0.5 D'(mu) (dW^2 - dt),
```

followed by exact repeated reflection into the closed physical interval. The
coefficient provider is called only with the pre-step bounded value. A
`ReflectingEulerMaruyama` enumerator remains explicit for controlled comparison;
there is no hidden clipping mode.

### Verification and limits

The WP11 ensemble starts from the exact uniform stationary distribution for
`D=D0(1-mu^2)`, includes particles at both endpoints, and checks the mean,
Legendre `P2` moment, endpoint population, symmetry, and provider domain-call
count. A native long-duration production ensemble remains part of the AMPS gate.

## WP12: configurable local error controls

### Why

Anonymous fractions such as `0.25` or `0.05` mix geometry, drift, stochastic,
cooling, focusing, shock, and minimum-step decisions. They cannot be audited or
refined independently and do not expose which constraint controls cost.

### What and how

`util/sep_transport_common.*` defines one validated `NumericalTolerances`
record containing geometry fraction, deterministic relative tolerance,
stochastic pitch RMS, cooling logarithmic change, focusing pitch change, shock
fraction, and minimum step. The three production movers read this active record.
`EstimateStepDoublingError` compares one full result with two half-step results,
and `StepDiagnostics` records accepted/rejected/error-control counts plus a
named limiter histogram.

The CLI exposes each tolerance with a dimensionally named option. Configuration
is validated transactionally before the active state is changed.

### Verification and limits

WP12 tests valid/invalid budgets, full-versus-two-half acceptance, minimum-step
selection, accepted-step accounting, and the limiter histogram. The present
production shells use all named physical limits and acceptance accounting; a
future package may extend automatic retry to every coupled operator without
changing this public error-budget contract.

## WP13: coefficient data bound to one authoritative source

### Why

Mixing magnetic field, turbulence variance, Alfvén speed, or spectrum bounds
from different generations can yield a plausible coefficient with no coherent
physical state. A fallback from self-consistent or SWMF data to a prescribed
global value also makes requested coupling silently ineffective.

### What and how

`util/sep_coefficient_physics.*` defines an immutable `LocalInputView` carrying
source, representation, generation, checksum, radius, magnetic field, total and
branch turbulence variances, and Alfvén speed. `coefficient_providers.cpp`
constructs that view through explicit prescribed, self-consistent, or SWMF
adapters, then passes it to pure coefficient kernels. Samples return the same
identity metadata. Incomplete selected data are errors; the adapter does not
consult an unrelated source as a fallback.

### Verification and limits

WP13 perturbs a selected turbulence variance and observes the expected
coefficient response, then perturbs an unused field and verifies invariance.
Source-level gates require the production provider to call the source-bound
builder. Native adapter tests must still verify host-specific datum offsets and
SWMF epoch publication.

## WP14: constant pitch-angle diffusion provider

### Why

The historical constant provider could fail to assign the configured value,
making a supposedly controlled scattering case behave ballistically.

### What and how

`EvaluateConstantDmumu` now returns the configured non-negative SI rate in
`s^-1` and an exact zero derivative. Negative, NaN, infinite, or out-of-domain
inputs return structured errors. The legacy callback delegates to this pure
kernel. `--constant-dmumu` is an explicit override; without that option a
post-compile input value is preserved and validated rather than reset.

### Verification

WP14 checks zero and positive values, the exact derivative, invalid values, and
production delegation.

## WP15: consistent Jokipii coefficient and derivative

### Why

Computing `Dmumu` and `dDmumu/dmu` from different expressions corrupts the Itô
drift. Finite spectral bands add genuine piecewise transitions that must not be
smoothed accidentally or confused with endpoint singularities.

### What and how

`EvaluateJokipiiSlab` normalizes a finite power-law slab spectrum to the selected
`deltaB^2`, computes the resonant wavenumber from the explicit species and local
field, and evaluates the value and derivative from the same smooth branch.
Outside-band resonance returns zero. Exact band transitions are typed as
`Nondifferentiable`; `mu=+/-1` uses finite one-sided derivative limits.

### Verification and limits

WP15 compares analytic derivatives with centered finite-difference oracles at
positive and negative pitch angles and checks both endpoints. Published-model
parameter choices remain run configuration, not hard-coded validation truth.

## WP16: repaired Florinskiy provider

### Why

The prior path omitted a final `D` assignment in one branch, read an
uninitialized derivative displacement, and used the plus-branch resonance
denominator in the minus contribution. These are direct defects rather than
model-choice differences.

### What and how

Both the pure and compatibility kernels assign `Dmumu`. The outward and inward
spectral contributions use `|v mu-vA|` and `|v mu+vA|`, respectively. Their
zero-denominator limits are evaluated analytically for spectral index greater
than one. Derivatives use bounded centered or one-sided stencils and never
sample beyond `|mu|=1`; the compatibility helper initializes its displacement
before reading it.

### Verification and limits

WP16 checks non-negativity, endpoint finiteness, and plus/minus mirror symmetry.
This repair establishes internal branch consistency. Equation-level provenance
for a selected Florinskiy publication and parameterization should additionally
be approved as part of scientific model review.

## WP17: error-controlled spatial-diffusion quadrature

### Why

`kappa_parallel=(v^2/8) integral (1-mu^2)^2/Dmumu dmu` may contain narrow
features, band transitions, or a resonance gap. A fixed rule can miss those
features, while a numerical floor silently turns a physical infinite/undefined
case into finite diffusion.

### What and how

`IntegrateSpatialDiffusion` uses adaptive Simpson subdivision with absolute and
relative tolerances, splits the domain at zero, counts evaluations, and returns
an error estimate. Interior zero `Dmumu` is handled only by the selected
`Reject` or `Ballistic` resonance-gap policy. Results distinguish finite,
ballistic, unavailable, and invalid states. The production spatial derivative
uses nested stencils in physical arc length and a refinement consistency check.

### Verification and limits

WP17 recovers the analytical result `v^2/(6D0)` for
`D=D0(1-mu^2)` and verifies that an interior gap is either rejected or returned
as typed infinity. No production kernel divides silently by zero.

## WP18: species-aware coefficients and sources

### Why

Gyrofrequency, resonance, rigidity, Larmor radius, energy per nucleon, and
injection abundance depend on mass and charge. A proton constant inside a
species-taking API can produce numerically smooth but physically incorrect ion
or electron transport.

### What and how

`SpeciesProperties` carries the PIC species identifier, name, signed charge,
rest mass, and nucleon count. Pure gyrofrequency, rigidity, Larmor-radius,
Jokipii, Florinskiy, and correlation-length MFP kernels require this record.
Absolute charge controls magnitude while the signed value remains available for
polarity-aware extensions.

`util/sep_species_source.*` adds a validated per-species source table with
normalized abundance, charge/mass metadata, injection efficiency, spectral
index, and an explicit total-energy or energy-per-nucleon convention.
`field_line.cpp` resolves one record per injection event, uses its abundance and
efficiency in source normalization, and converts energy per nucleon to total
joules before relativistic momentum conversion. Installing a partial explicit
table is an error; an empty table preserves the documented legacy defaults.

### Verification and limits

WP18 checks proton/alpha/electron gyrofrequency, Larmor-radius, and rigidity
scaling; exact two-species abundance closure; alpha energy-per-nucleon
conversion; and rejection of an electron per-nucleon definition. PIC hosts that
provide authoritative isotope metadata should populate nucleon count directly;
the compatibility adapter estimates it from mass only when that host metadata
is unavailable.

## WP19: named turbulence scales and amplitude policy

### Why

Embedded values for `B`, `deltaB/B`, correlation length, or spectral cutoffs
make diagnostics disagree with the selected background and prevent meaningful
configuration fingerprints.

### What and how

The coefficient configuration now owns prescribed `deltaB/B`, correlation
length at 1 AU, spectral reference radius and bounds, independent radial
exponents, and the Florinskiy correlation length. Production adapters use the
selected source's local magnetic field and turbulence energy. The sampling
Larmor-radius diagnostic uses the actual segment field rather than a QLT1
constant. If turbulence amplitude exceeds the mean-field policy boundary, the
run either rejects the state or applies the explicitly selected limiter and
increments a regularization ledger counter.

### Verification

WP19 rejects negative amplitude, invalid correlation length, and unordered
spectrum bounds. Source gates reject a reintroduced hard-coded QLT1 field.
Startup metadata prints the active named scales and policies.

## WP20: explicit ballistic mean-free-path compatibility

### Why

`lambda=+infinity` is a valid zero-event-rate limit for the event-driven mover,
but `kappa=v lambda/3=+infinity` is not a finite coefficient accepted by the
Parker diffusion operator. Treating every infinity alike causes either an
unphysical clamp or a numerical failure inside the mover.

### What and how

Mean-free-path and spatial results carry `ValueState`; ballistic is distinct
from unavailable or invalid. Preflight compatibility rejects Parker
configurations that can feed ballistic MFP or resonance-gap states into its
finite diffusion operator. `fte-mfp` accepts typed ballistic state and maps it
to zero scattering rate. Scalar legacy conversion helpers deliberately reject
infinity so unsupported code cannot inherit it accidentally.

### Verification

WP20 checks zero-turbulence ballistic MFP, event-mover acceptance, Parker
preflight rejection, and rejection by the finite `lambda`-to-`kappa` helper.

## Files and review boundary

| Area | Principal files |
|---|---|
| Bounded stochastic integration and error controls | `util/sep_focused_transport_core.*`, `util/sep_transport_common.*`, three canonical mover shells |
| Pure coefficient physics and source authority | `util/sep_coefficient_physics.*`, `util/sep_coefficient_registry.*`, `coefficient_providers.*` |
| Legacy callback repairs | `diffusion.cpp`, `diffusion_dxx.cpp`, `sep.h` |
| Species source and diagnostics | `util/sep_species_source.*`, `field_line.cpp`, `sampling.cpp` |
| Configuration and metadata | `util/sep_cli.*`, `production_mover_runtime.cpp`, `main.cpp` |
| Focused evidence | `test/test_wp11_wp20.cpp`, `test/run_wp11_wp20_tests.sh` |

Run `make test-wp11-wp20-unit` in a source-only checkout. In a configured AMPS
application, additionally build under strict warnings and run
`make test-native-amps-validation SEP_EXECUTABLE=/path/to/amps`. Real SWMF and
held-out observational evidence remain the separate gates described in
`validation/README.md`.
