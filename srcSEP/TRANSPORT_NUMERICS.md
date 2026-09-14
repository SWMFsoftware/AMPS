# Field-line transport and turbulence numerical contract

Steps 6–10 define one numerical vocabulary for every public srcSEP mover and one
canonical implementation each for the Parker and coefficient-driven focused
transport equations. The kernels under `util/` intentionally depend on neither
PIC nor MPI, while `transport_common.*`, `coefficient_providers.*`, and the
three canonical mover adapters connect those kernels to the host application.

## State, frames, and units

The common particle state contains species, field-line identifier, host line
coordinate, particle mass [kg], and local parallel/normal velocities [m/s]. The
host coordinate is not assumed to be arc length. Every core returns a physical
displacement [m], and only `PICAdapter::AdvanceAlongFieldLine` passes that
distance to `FieldLine::move`; this is the single metric conversion.

Momentum is relativistic SI momentum [kg m/s]. `MomentumFromSpeed` requires
`0 <= v < c`; `SpeedFromMomentum` accepts finite nonnegative momentum. Neither
routine clips an invalid state. Cooling is evaluated in the local plasma frame:

`dp/dt = -(p/3) div(U)`, therefore `p(t+dt)=p(t) exp[-div(U) dt/3]`.

The immutable background read phase supplies `U_parallel`, `dU_parallel/ds`,
`div(U)`, `d ln|B|/ds`, Alfvén speed, and snapshot validity from one generation.
The focusing length convention
is `L=-1/(d ln|B|/ds)`; a uniform magnetic field has infinite `L` and zero
focusing rather than an error.

## Status, boundaries, and substeps

Numerical boundaries return `Status`/`StatusCode`. Invalid particle state,
coefficient, argument, domain exit, and timestep underflow remain distinct.
Production adapters translate an unexpected non-OK result to the host fatal
error path; focused tests inspect it directly.

Physical line-coordinate boundaries declare either absorption or reflection.
Pitch angle always uses repeated reflection into `[-1,1]`; the modulo-folding
implementation handles a stochastic increment that crosses multiple boundaries
and keeps `mu=0` regular.

Each mover assembles named limits for its active operators and selects their
minimum with `SelectSubstep`. The selector records the limiting name and count.
If that minimum falls below the declared threshold it reports `StepUnderflow`;
it never increases the interval with an undocumented floor.

This complete composable-limit contract is used by all three canonical movers.
The event-driven mover additionally splits at the current shock and immutable
snapshot-validity boundaries.

## Deterministic stochastic streams

`KeyedRandomStream` is constructed from campaign seed, particle identity,
operator identity, and event/substep index. MPI rank and worker number are not
keys. The open-interval uniform generator supports Box–Muller normals without a
singular logarithm, and an identical key reproduces an identical stream
independently of scheduling.

WP24 additionally makes key composition ordered and purpose tagged. The source
tuple is `(campaign,event,field-line,species,macroparticle,purpose)`; swapping
two fields changes the stream, and count/spectrum/pitch/gyrophase draws cannot
perturb one another. The primary field-line shock source no longer consumes the
process-global RNG.

## Shock/source and diagnostic numerical contracts

WP21 integrates a piecewise-linear `v_shock(r)` exactly from a serialized launch
datum, and WP22 intersects its spherical surface with the complete 3-D
field-line polyline using oriented/tolerance-aware typed results. WP25 evaluates
the swept shell volume and upstream normal-relative kinetic energy, queues the
result, and lets the authoritative turbulence core apply and ledger it.

WP27 diagnostics reject invalid or luminal particle state rather than clipping
it. Relativistic Larmor radius is `p_perp/(|q|B_local)`. WP28 uses half-open bins
with separate underflow/overflow/invalid counters, preserves `sum(w^2)` for
effective sample size/error, and constructs number density, omnidirectional
differential intensity, directional crossing flux, and normalized shape as
separate unit-bearing products.

`CampaignRandomSeed` is the production adapter boundary. A campaign that needs
a nonzero seed must set it before particle motion; leaving it at zero is still
deterministic and explicit.

## Parker equation

`SpatialDiffusionProvider` returns `kappa_parallel` [m²/s],
`d(kappa_parallel)/ds` [m/s], status, and provenance together. The Itô update is

`ds = [U_parallel + d(kappa_parallel)/ds] dt
      + sqrt(2 kappa_parallel dt) dW`.

The production provider delegates to the configured srcSEP `GetDxx` contract;
a missing configured diffusion callback is the explicit zero-diffusion case.
Convection, variable-coefficient drift, Gaussian noise, and exact cooling are
composed in `AdvanceParker`. Segment-scale drift and diffusion constraints are
applied outside the pure core by the PIC adapter.

## Coefficient-driven focused transport

`PitchAngleDiffusionProvider` returns `Dmumu` [s^-1], `dDmumu/dmu` [s^-1],
status, provenance, and the immutable turbulence-state identity used to obtain
the coefficient. A configured analytical derivative is used directly; the
numerical mode evaluates a symmetric stencil except for one-sided endpoints.

One substep uses a symmetric deterministic/stochastic/deterministic split:

1. half of the focusing, parallel-flow-gradient, and cooling operators;
2. one Itô kick `dmu = (dDmumu/dmu)dt + sqrt(2 Dmumu dt)dW`;
3. the remaining deterministic half; and
4. streaming with the midpoint pitch angle and local plasma advection.

The deterministic pitch-angle half-steps use explicit midpoint integration.
The production adapter limits streaming, focusing, diffusion amplitude, and
diffusion drift separately.

## Event-driven mean-free-path focused transport

`MeanFreePathProvider` returns `lambda_parallel` [m], branch rates `nu+` and
`nu-` [s^-1], status, provenance, optional momentum validity bounds, and the
turbulence-state identity. Positive infinite lambda is the explicit zero-rate
ballistic limit. A particle stores the remaining unit-exponential optical depth
`tau=-log(xi)`. Each deterministic interval integrates the nonhomogeneous
hazard `integral(nu+ + nu-) dt`; a monotone root solve locates an event inside
the interval, and only a completed event draws the next optical depth. Segment,
shock, snapshot, and outer-timestep partitions therefore do not resample the
physical process.

Between events, the mover applies midpoint focusing and parallel-flow-gradient
drift, exact plasma-frame cooling, and midpoint streaming. At an event the
branch is sampled with probabilities `nu+/(nu+ + nu-)` and
`nu-/(nu+ + nu-)`; an empty branch is never selected. The particle direction
is redistributed isotropically in that signed Alfvén wave frame. Exact Lorentz
velocity transforms preserve wave-frame speed and keep the returned
plasma-frame state subluminal.

## Unified coefficient registry

`util/sep_coefficient_registry.*` publishes canonical names, SI units, and
parameter schemas for spatial, pitch-angle, and mean-free-path coefficients.
The PIC-facing providers add background ownership validation and provenance.
The only supported conversions are centralized as `kappa=v*lambda/3` and the
explicit isotropic `Dmumu` closure. A spatial-from-MFP/MFP-from-spatial cycle is
rejected. Self-consistent and SWMF MFP transport must use `from-spatial`; an
SWMF source is valid only with a read-only SWMF snapshot.

Invalid lambda output follows an explicit policy. `fail` returns a structured
error; `ballistic` substitutes positive infinity and increments both invalid-
sample and ballistic-substitution counters. The canonical mover has no silent
minimum/maximum lambda clamp.

## QLT coefficient normalization

For `P(k)=C k^(-5/3)` on `[kmin,kmax]`, the normalization is

`C=(2/3) deltaB^2/(kmin^(-2/3)-kmax^(-2/3))`,

so `integral P(k) dk = deltaB^2`. With
`k_res=Omega/(v |mu|)`, the slab expression is

`Dmumu=(pi/2) Omega^2 (1-mu^2) P(k_res)/(B^2 v |mu|)` [s^-1].

An invalid input or a resonance outside the represented band returns zero from
the historical QLT API. Provider adapters subsequently reject any non-finite or
negative coefficient before advancing a particle.

## Full focused equation and authoritative background view

Production `fte-dmumu` and `fte-mfp` use the local plasma-frame gyrotropic SDE

`dmu/dt = (1-mu^2)/2[-v dln|B|/ds + mu(divU - 3 bb:gradU)] + dDmumu/dmu`,

`dln(p)/dt = -1/2[(1-mu^2)divU + (3mu^2-1)bb:gradU]`.

The stochastic pitch increment is `sqrt(2 Dmumu dt) dW`. The explicit
`ReducedFieldAligned1D` mode preserves the earlier one-dimensional closure for
documented comparisons, but production adapters select `FullGyrotropic`.
`LocalBackgroundView` is the single SI record supplying magnetic field/vector,
unit field direction, line tangent, plasma velocity, number/mass density,
`divU`, `bb:gradU`, Alfvén speed, provider provenance, generation, and validity.
Providers evaluate requested physical arc-length displacements and the Dmumu
core samples its predicted midpoint rather than a hard-coded location.

Density-derived divergence uses the interval between physical background
epochs stored by `BackgroundSnapshot`, never a particle timestep. All three
movers compose the same snapshot validity/generation limit before every shell;
a stale generation or an attempt beyond the validity endpoint is an error.

## Wave feedback and concurrency

The numerical core deposits only a signed streaming record labelled with the
same turbulence identity carried by its coefficient. PIC-facing movers never
write shared G+/G- wave arrays from worker threads. They append a typed,
self-contained record after every completed deterministic interval/event.
Species and statistical weight are copied while the particle is alive; the
post-step queue never retains or dereferences a PIC handle. If an absorbing
boundary is crossed, the attempted trajectory is clipped to the exact endpoint,
interval time is scaled to the in-domain fraction, the boundary record is
enqueued, and only then is the particle deleted.

Step 12 formalizes the next reduction boundary in
`util/sep_reproducible_reduction.*`. Thread-local records carry a versioned
schema plus source, field line, segment, branch, spectral bin, species, stable
particle, step, event, interval, and purpose keys. Duplicate complete keys are
errors. After workers join (and after rank-local records are gathered), a
complete-key sort establishes one fixed
accumulation order with long-double intermediates. Non-additive authoritative
wave state uses owner/gather synchronization, never an `MPI_Allreduce`
overwrite. Canonical evidence hashes are bitwise invariant to the tested worker
and synthetic rank partitions.

The retained turbulence driver, ownership, boundary, ledger, remap, and restart
contract is documented in [TURBULENCE_MODEL.md](TURBULENCE_MODEL.md).

## WP33 population control and WP36 system accounting

`RunConfiguration::populationControl` owns the actual spatial/momentum/pitch
bin counts, minimum/maximum population, deterministic-selection requirement,
invariant tolerance, and lineage schema passed to PIC merging and splitting.
`sep_population_control.*` provides the independent SI moment observer used by
controlled and native fixtures. Split children inherit phase-space state and
partition weight with a final-child residual; stable child IDs hash parent,
generation, ordinal, and schema without rank or thread.

`sep_system_ledger.*` closes number, charge, energy, and parallel momentum over
named signed seam transactions. Its checksum-protected checkpoint carries
field-line generation, stable IDs, lineage, RNG counters, residual scattering
hazards, turbulence-state identity, and output-manifest identity. This source
contract is ready for native adapters; complete injection/escape/coupling/remap/
restart equality still requires a linked AMPS campaign.

## Verification boundary

The dependency-light suites are `make test-transport-common-unit`,
`make test-parker-unit`, `make test-fte-dmumu-unit`,
`make test-fte-mfp-unit`, `make test-coefficients-unit`,
`make test-turbulence-core-unit`, and `make test-reproducibility-unit`. They use strict C++11
warnings, AddressSanitizer, and UndefinedBehaviorSanitizer and cover `CORE01`–
`CORE07`, `PARK01`–`PARK07`, `FTED01`–`FTED08`, `FTEM01`–`FTEM08`, and
`COEF01`–`COEF05`, `TURB02`–`TURB23`, `TURBOWN01`, and `PAR01`–`PAR05`.
The mover and turbulence cases are descriptors in the same registry consumed
by `--test` and `--test-group`; their source-only runners no longer carry
parallel copies of the callbacks. `FTED08` adds the analytical Legendre-mode
decay reference, while `TURB21` and `TURB22` add translated-profile advection
and time-dependent wave-energy-source references. `TURB22` verifies how a
known source is integrated and accounted; it does not replace the native
particle-distribution-to-growth-rate test. Source checks additionally
verify public dispatch and the common adapter boundary.

Refinement evidence is evaluated uniformly as
`p=log(e_coarse/e_fine)/log(h_coarse/h_fine)`. The common helper rejects zero,
negative, nonfinite, or reversed-resolution inputs so an empty/error-free
special case cannot masquerade as convergence. `PARK07` adds exact absorbing
first-passage statistics, and `FTEM08` adds the finite-time persistent-flight
MSD and its `kappa_parallel=v*lambda_parallel/3` diffusion asymptote.

A complete AMPS checkout is still required to compile the PIC-facing adapters,
run a particle through every registered shell, compare serial/OpenMP/MPI
histograms, exercise SWMF coupling epochs, and execute long-run conservation
and observational validation. Those native and scientific gates cannot be
claimed from this source-only archive.
