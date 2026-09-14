# WP01--WP10 implementation and validation record

This document maps the physics/numerics/testing improvement plan to the source
changes in this delivery. All internal quantities are SI. “Validated” below
means the bounded source-only controlled gate passed; native AMPS, MPI, SWMF,
restart, performance, and observational gates remain separate evidence classes.

## WP01 — production turbulence-core integration

Why: standalone `main.cpp` previously owned a long direct sequence while the
common `Turbulence::Advance` driver was exercised primarily by controlled tests.
That allowed production ordering and tested ordering to diverge.

What/how: `turbulence_production_adapter.*` provides the sole reachable
production mutation entry, called by standalone and coupled-library stepping.
It flushes completed particle records, contains legacy shock/resonant-growth
logic as source adapters, maps PIC field-line geometry/background/energy into
the common state, calls `Turbulence::Advance`, exports the authoritative
integrated or branch-major spectral state, and derives W+/W-/cross helicity.
The old inline block is constant-false migration evidence and is not executable.

Gate: source checks assert both drivers call the adapter. Native field-line
mapping, MPI ownership, and long-run energy equivalence require a configured
AMPS executable and are not claimed here.

## WP02 — reproducible production reductions

Why: pointer and worker ordering are allocation/scheduling artifacts and cannot
define floating-point summation or keyed stochastic identity.

What/how: the reduction key now has a restart-visible schema plus source,
field line, segment, branch, spectral bin, species, stable particle, step,
event, interval, and purpose. Duplicate complete keys are rejected. Canonical
sorting and long-double accumulation remain common to worker and synthetic-rank
partitions. Production flush constructs these common contributions before
path-resolved deposition and no longer sorts by a particle-buffer address.

Gate: the focused test proves duplicate rejection and the existing PAR01--PAR05
suite proves scheduler/partition invariance. A real multi-rank gather must still
be exercised by the native MPI gate.

## WP03 — event-resolved, escape-safe coupling

Why: deferring a live particle handle until after an absorbing exit discarded
valid in-domain work and could dereference deleted storage.

What/how: `CouplingRecord` copies stable identity, species, statistical weight,
snapshot generation, pre/post momentum, midpoint velocity, event/interval,
path, branch, and boundary metadata. The coupling implementation has payload
overloads that need no live particle. Parker and both focused movers enqueue
completed intervals. Exits are clipped to the first physical endpoint; record
time is scaled to the in-domain path and publication precedes deletion. The MFP
core emits one wave record per deterministic/event interval.

Gate: static checks reject a stored `particlePointer` and require all movers to
use `QueueWaveContribution`. Native energy ledgers across actual exits remain a
linked integration test.

## WP04 — all-mover snapshot validity

Why: using a newly published geometry/background in the middle of a particle
step violates the frozen-coefficient SDE and can mix generations in feedback.

What/how: `ComposeSnapshotValidityLimit` compares the acquired and sampled
generations and converts the immutable validity endpoint into remaining seconds.
Parker, Dmumu, and MFP add it to their named limiter set before every shell.
Exact exhaustion is an error requiring the driver to end the read phase.

Gate: the focused test covers a valid residual interval and stale-generation
rejection. Native tests must shorten an actual provider validity interval for
each mover.

## WP05 — physical background epochs

Why: dividing density change by a particle timestep makes inferred `div(U)`
change when only numerical particle refinement changes.

What/how: immutable snapshots carry previous/current physical provider epochs.
Runtime publication advances that pair from analytic/SWCME or SWMF coupling
times. `EvaluateLocalBackgroundAt` computes density-derived divergence only
from that interval. A single static epoch is permitted only when current and
previous density agree; otherwise the state is rejected.

Gate: constructor/interval behavior is controlled-tested. Native expansion,
compression, refinement, coupling, and restart histories remain required.

## WP06 — authoritative local background view

Why: movers previously reconstructed partially overlapping B, flow, density,
Alfvén speed, and derivative quantities through heterogeneous fallbacks.

What/how: `LocalBackgroundView` publishes B vector/magnitude, field unit vector,
line tangent, plasma velocity, number/mass density, `dln|B|/ds`, projected flow,
parallel gradient, `div(U)`, `bb:grad(U)`, Alfvén speed, requested coordinate,
provider provenance, generation, identity, and validity. All values are produced
by one status-returning SI adapter; production movers consume this view.

Gate: full coefficient calculations use independently supplied divergence and
field-aligned strain. SWMF interpolation-order/provenance evidence remains a
native coupling gate.

## WP07 — focused transport equation contract

Why: the older reduced closure mixed focusing, parallel velocity gradients,
and isotropic Parker cooling without freezing the general gyrotropic equation.

What/how: public helpers implement

`dmu/dt=(1-mu²)/2[-v dln|B|/ds + mu(divU-3 bb:gradU)]`

and

`dlnp/dt=-[(1-mu²)divU+(3mu²-1)bb:gradU]/2`.

Production focused movers explicitly select `FullGyrotropic`. The earlier
one-dimensional equation remains only as `ReducedFieldAligned1D` for controlled
compatibility comparisons. Momentum is advanced exponentially for a frozen
midpoint pitch angle, retaining positivity.

Gate: the focused test compares both public coefficients with independent
formulas. Physics-owner approval and independent finite-volume combined-profile
validation remain release requirements.

## WP08 — location-aware variable coefficients

Why: accepting an `s` argument but evaluating the captured starting segment
locks coefficients to a grid location and destroys variable-profile convergence.

What/how: PIC coefficient providers translate every requested physical
arc-length displacement through the field-line metric before evaluating B,
turbulence, lambda, kappa, or Dmumu. The Dmumu split samples at a predicted
physical midpoint. MFP samples both ends of each limited deterministic interval.
Out-of-domain samples return explicit status rather than extrapolating.

Gate: a recording provider proves the core requests a moved midpoint. Native
manufactured smooth/discontinuous profiles and grid-lock counts remain needed.

## WP09 — nonhomogeneous scattering hazard

Why: resampling an exponential waiting time whenever geometry or snapshots split
a step changes a spatially varying Poisson process and makes answers depend on
the limiter partition.

What/how: MFP state carries residual unit-exponential optical depth and event
index. Each interval uses endpoint rates, integrates the linear hazard, locates
an interior event by stable monotone bisection, subtracts hazard when no event
occurs, and draws a new depth only after an event. Zero total rate is ballistic.

Gate: a constant-rate history is bitwise-consistent between one call and a
0.3+0.7 partition using the same keyed stream. Linear/exponential/piecewise
statistical campaigns should be added to the extended native suite.

## WP10 — branch-resolved event rates

Why: choosing a resonant branch from sign(mu), with a 50/50 special case at
zero, ignores imbalanced wave populations and can select a physically empty
branch.

What/how: `MeanFreePathSample` publishes `nuPlusPerS`, `nuMinusPerS`, and an
explicit availability flag. Total hazard is their sum; event-time conditional
selection uses each branch fraction. The production PIC adapter obtains the
total rate from its lambda closure and partitions it with the authoritative
local integrated E+/E- populations, so an empty wave population has zero rate.
Providers without active turbulence map explicitly to balanced rates. A future
spectral closure may refine these integrated proportions at the resonant bin.

Gate: the focused all-plus test produces plus events and exactly zero minus
events. Balanced, imbalanced, polarity, and spectral-provider campaigns remain
native/extended work.

## Commands and evidence boundary

Run:

```sh
./test/run_wp01_wp10_tests.sh
./test/run_step6_tests.sh
./test/run_step8_tests.sh
./test/run_step9_tests.sh
./test/run_step11_tests.sh
./test/run_step12_tests.sh
```

The archive does not contain the enclosing AMPS executable, so `make -j`,
native MPI/OpenMP, SWMF, restart, sanitizer coverage of PIC adapters, long-run
energy, performance, and observational gates must be reported as unavailable,
not PASS.
