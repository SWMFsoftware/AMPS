# WP21–WP30 physics, numerics, sampling, and run-control implementation

This document is the implementation contract for work packages WP21 through
WP30. All dimensional interfaces use SI units. Dependency-free kernels live in
`util/`; PIC-facing files translate field-line state into those records and do
not duplicate their formulas. The work packages were implemented in numerical
dependency order: shock state, geometry, spectrum/RNG, turbulence source,
magnetic flux, diagnostics/products, output transaction, then run configuration.

## WP21 — analytical shock kinematics and restart state

Why: the former 0.3–0.5 AU interpolation divided by `0.5-0.1`, producing a
speed discontinuity, and advanced radius by forward Euler from an unserialized
function-static time. A nonzero start or different timestep partition therefore
changed the trajectory.

What/how: `util/sep_shock_source_core.*` defines ordered `(radius [m], speed
[m/s])` knots and analytically integrates each linear `v(r)` interval. The
constant-speed limit uses `(r1-r0)/v`; nonzero slope uses the exact logarithmic
travel time and `expm1` inverse. Every epoch is re-evaluated from the serialized
launch datum, giving bitwise partition invariance. The versioned trajectory
checkpoint includes launch epoch/radius, current epoch/radius, and all knots.
`shock_analytical_model.cpp` updates this state before any radius-derived query.

## WP22 — robust shock/field-line intersections

Why: the old search required the first field-line point to be inside the shock,
accepted only inside-to-outside crossings, fabricated a segment-start result for
a negative discriminant, and divided by zero on degenerate segments.

What/how: `IntersectSphere` is a pure 3-D polyline/sphere kernel. It validates
segment length, applies a physical-distance tolerance to discriminants and
endpoints, coalesces a shared-vertex duplicate, and records segment fraction,
arc length, Cartesian point, and inward/outward/tangent orientation. Callers
choose `First`, `FirstOutward`, or `All`; no intersection, ambiguity, and invalid
geometry are typed results. The analytical injection adapter uses
`FirstOutward`, including when the polyline begins outside or crosses repeatedly.

## WP23 — defined and normalized injection spectra

Why: the Sokolov normalizer raised `log(p_min)` to a fractional power and could
produce NaN. The Tenishev path sampled uniformly in log momentum but did not
state whether its exponent described `f(p)`, `dN/dp`, or `dN/dlog(p)`.

What/how: `util/sep_injection_spectrum.*` names the sampled measure and evaluates
normalized power-law densities/inverse CDFs with explicit index-one logarithmic
limits. In the Tenishev production route, diffusive-shock `f(p)∝p^-q` becomes
the number spectrum `dN/dp∝p^(2-q)` before direct inverse-CDF sampling. The
Sokolov index is explicitly `dN/dp`. Both now have unit statistical correction;
the optional log-uniform importance ratio remains available and normalized.
Bounds and non-finite results return structured errors.

## WP24 — explicit source configuration and reproducible source RNG

Why: source draws used the process-global `rnd()` stream, and the earlier keyed
stream XORed independent field hashes. XOR is commutative, so swapping two key
fields could yield the same sequence; scheduling and unrelated diagnostic draws
could also perturb injection.

What/how: injection configuration declares campaign seed, macroparticles/event,
efficiency, coordinate measure/bounds/index, and angular law. `RandomKey`
contains campaign, physical event, field line, species, macro index, and a
semantic purpose (`count`, `spectrum`, `pitch`, `gyrophase`, or `position`).
Tagged sequential hash combination is position sensitive. The principal
field-line shock source now uses separate keyed streams for event rounding,
spectrum, pitch, and gyrophase, with no buffer address/rank/thread key.

## WP25 — shock turbulence through the common ledger

Why: the former function silently returned on bad input, converted radial
offsets into field-line arc coordinates, assumed absolute radial shock speed,
hard-coded source parameters, and mutated the legacy energy datum directly.

What/how: each segment computes its actual overlap with the swept spherical
shell, shared WP26 partial volume, upstream density, radial normal, and upstream
normal flow. `ComputeTurbulenceSource` applies
`E=eta rho V (Vshock,n-Upstream,n)^2/2` and the configured branch split. The
adapter queues typed `(line, segment, E+, E-, provenance)` records. Import maps
them to `pendingShockPlusJ/MinusJ`; only `Turbulence::Advance` applies them,
limits negative requests, and records the applied source in its signed ledger.
Zero relative speed is an exact zero source. Invalid geometry/physics are counted.

## WP26 — physical magnetic flux per field line

Why: all areas used one implicit global reference area even though each field
line represents a distinct seed-surface element.

What/how: `MagneticFluxRecord` owns positive finite `Phi_i [Wb]`, generation,
and provenance. The PIC adapter converts the configured seed area and the first
valid field value into a record once; all later areas use `A=Phi_i/|B|`.
Conservative refinement takes positive fractions summing to one and assigns the
floating-point residual to the final child. Versioned table serialization
preserves flux, generation, and provenance across restart. The historical
`pi m2` seed remains only as a named, fingerprinted compatibility default and
is overridable through `--field-line-seed-area`.

## WP27 — invalid particles and Larmor diagnostics

Why: diagnostics capped a luminal particle to `(1-epsilon)c`, thereby creating
a finite but fabricated event. Legacy Larmor calculations also risked using a
nonlocal field or nonrelativistic `m v_perp`.

What/how: `EvaluateKinematics` accepts only finite, nonzero, strictly subluminal
states with valid mass, charge, and local field. It returns relativistic
`p=gamma m v`, kinetic energy, `mu`, and `rL=p_perp/(|q|B_local)`. Invalid states
are excluded (or may be promoted to fatal by policy) and contribute no bin
weight. Production sampling now rejects rather than caps `v>=c`.

## WP28 — physical sampling products and overflow semantics

Why: underflow/overflow values were clamped into edge bins and output arrays
were first normalized by their integral and then again by their maximum, making
units and conservation unknowable.

What/how: `WeightedBins` keeps in-range `sum(w)`, `sum(w^2)`, and separate
underflow/overflow/invalid/missing counters. Bins are half-open `[edge_i,
edge_i+1)`; the final upper edge is overflow. `BuildProduct` separately creates
number density per energy, omnidirectional differential intensity, directional
crossing flux, or unit-integral shape, with the required volume/area/time/speed
normalization and propagated weighted-count standard error. Effective sample
size is `(sum w)^2/sum(w^2)`. Legacy shape outputs no longer receive a second
peak normalization and identify schema/product/units in their header.

## WP29 — safe transactional writes

Why: shell `mkdir -p`, fixed-size `sprintf`, unchecked `fopen`, and direct final
writes exposed command injection, truncation, partial artifacts, and corrupt
restart reads.

What/how: `util/sep_transactional_output.*` validates relative paths, creates
directories with POSIX APIs, refuses duplicate final names, writes a temporary
file with checked short-write handling, `fsync`s it, and atomically renames it.
A completion manifest records schema, record count, bytes, configuration
fingerprint, checksum, and `complete=true`; readers reject missing, partial, or
mismatched data. Energy, Larmor, MFP, and pitch-angle sampling outputs now stage
then commit through this API. No shell or unchecked `sprintf` remains in
`sampling_output.cpp`.

## WP30 — one frozen RunConfiguration

Why: shock model/scenario, iteration count, seed area, turbulence source,
source parameters, and merge/split thresholds were assigned in different driver
locations, sometimes after parsing. A static SWCME initialization flag also
made reinitialization order-dependent.

What/how: `util/sep_run_configuration.*` collects mover, shock/scenario,
iterations, per-line flux seed, shock-wave efficiency/split, merge limits,
sampling policy, coefficient/tolerance, turbulence, and injection records.
Layered values enforce `default < input-file < command-line`. Validation occurs
before `FrozenConfiguration::Create`; consumers receive const views. A canonical
startup fingerprint and versioned snapshot support exact roundtrip and restart
compatibility rejection. The standalone driver freezes before model/mesh init,
prints the fingerprint and value provenance, sets the shock/scenario and flux
seed from the frozen view, and the coupled library reads merge/split limits from
the same active record. The SWCME helper is deterministic and has no static
one-time flag.

New CLI controls are `--total-iterations`, `--shock-model`, `--cme-scenario`,
`--field-line-seed-area`, `--shock-turbulence-efficiency`,
`--shock-turbulence-plus-fraction`, `--merge-minimum`, and `--merge-maximum`.

## Verification and evidence boundary

Run `make test-wp21-wp30-unit`. The runner compiles the exact pure kernels with
C++11, strict warnings, ASan, and UBSan, then checks production wiring and the
absence of shell/`sprintf` use in sampling output. It covers one named assertion
per work package, including trajectory/restart bitwise equality, multi-crossing
orientation, spectrum normalization, key order, source-energy closure, flux
restart/refinement, relativistic diagnostics, bin/product normalization,
transaction corruption/duplicate rejection, and run-config precedence/roundtrip.

This source-only evidence does not establish a linked AMPS build, a real MPI
decomposition, SWMF replay, or observational validity. Those remain the separate
native, SWMF, and observational gates documented under `validation/`.
