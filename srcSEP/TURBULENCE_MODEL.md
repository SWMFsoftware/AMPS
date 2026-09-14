# Alfvén-turbulence ownership and evolution contract

Steps 11 and 12 retain self-consistent Alfvén turbulence as a physics
subsystem; it is not a fourth particle mover. The exactly three production
movers remain `parker`, `fte-dmumu`, and `fte-mfp`.

## Source and authoritative state

`util/sep_turbulence_core.*` defines five explicit source modes:

| Source | Mutable here? | Authoritative representation |
|---|---:|---|
| `prescribed` | No | selected prescribed state |
| `self-consistent-integrated` | Yes | segment-integrated `E+`, `E-` [J] |
| `self-consistent-spectral` | Yes | branch-major `E±(k)` [J per log bin] |
| `swmf-read-only` | No | immutable imported state |
| `swmf-initial-then-local` | After one handoff | selected local representation |

Integrated and spectral arrays cannot be simultaneous mutable owners. In
spectral mode, branch-integrated values are derived by summing bins. In
integrated mode, spectral storage must be absent. An SWMF-to-local handoff is a
one-time operation and records its epoch and source checksum.

Wave density is derived as `W±=E±/V` [J/m³], amplitude as
`deltaB²=2*mu0*(W+ + W-)`, and cross helicity as
`sigma_c=(E+ - E-)/(E+ + E-)`. Coefficient views derive resonance, branch,
`Dmumu`, parallel mean free path, and parallel spatial diffusion from this same
state; they do not own a second copy.

## Driver order and stability

`SEP::Turbulence::Advance` is the context-neutral driver for standalone and
coupled state adapters. Its fixed order is:

1. shock and other external wave sources;
2. particle-wave energy exchange;
3. finite-volume advection with explicit CFL subcycling;
4. branch reflection;
5. cascade and physical dissipation;
6. authoritative/derived representation synchronization; and
7. diagnostics, ledger accumulation, and checkpoint phase reset.

`SEP::Turbulence::PICAdapter::Advance` is the only reachable production
mutation entry in both `main.cpp` and `main_lib.cpp`. It flushes typed particle
records, runs the retained shock and resonant-growth implementations strictly
as source adapters inside one transaction, imports each field line into the
common `State`, invokes the core for advection/reflection/cascade, exports the
authoritative integrated or spectral representation, and refreshes derived
energy density/cross helicity. The prior inline driver block remains guarded by
constant false only as migration-review evidence and cannot be executed.
`LastStepLedger()` exposes measured shock/particle source changes plus summed
core boundary exchange, dissipation, limiter correction, closure residual, and
subcycle count; source-adapter mutations are therefore visible rather than
being hidden before the common-core ledger.

The standalone loop installs the same active source, operator, and coefficient
configuration. Prescribed and `swmf-read-only` sources bypass every mutating
operator but remain available for synchronization and diagnostics. An invalid
or non-positive production CFL limit is fatal; it is never replaced by the
remaining particle timestep.

Reflection is an internal signed branch transfer and preserves total energy.
Particle and shock withdrawals are limited by available branch energy, with
the requested/applied difference reported as a limiter correction. The
integrated reference cascade uses an implicit positive turnover sink.

### WP31 shared operator limits

The declared production scheme is `LieFirstOrder`. Before mutation,
`PlanAdvance` evaluates the finite-volume CFL bound, pending-source change,
reflection rate, and nonlinear cascade rate in every segment. The smallest
limit determines a common stage count. A plan below `minimumSubstepS` or above
`maximumSubsteps` is rejected with `StepUnderflow`; no operator substitutes the
global step as a fallback. Reflection uses its exact exponential branch exchange
at each stage and cascade uses a positive backward-Euler update. The focused
combined-operator case automatically measures first-order temporal convergence.

| Field | Meaning | Unit/domain |
| --- | --- | --- |
| `operatorAccuracySafety` | reflection/local-rate accuracy fraction | dimensionless `(0,1]` |
| `maximumSourceFraction` | pending source relative to current energy | dimensionless `(0,1]` |
| `maximumCascadeFraction` | cascade transfer per stage | dimensionless `(0,1]` |
| `minimumSubstepS` | smallest accepted common stage | seconds, positive |
| `maximumSubsteps` | work/underflow guard | positive count |

### WP32 statuses and WP40 work counters

Every correction creates an `OperatorEvent` with requested/applied joules and a
typed disposition. The signed unapplied amount enters `rejectedSourceJ`; exact
zero sources and empty physical wave states are counted as physical zeros.
Nonfinite rates, invalid density/volume, and zero-pitch resonance fail explicitly.
Diagnostics also count advection faces, reflection/cascade cells, spectral-bin
updates, limiters, rejections, and per-operator stages. The PIC adapter exports
these deterministic cost drivers with `LastStepLedger()`.

## Boundary and conservation ledger

Each field-line end owns one restartable boundary object:
`specified-incoming-energy`, `specified-incoming-flux`,
`transparent-outflow`, or `fixed-reservoir`. Energy values use joules and flux
values use inward watts. Non-transparent negative or non-finite values are
rejected.

Every update reports signed initial/final energy, inner and outer boundary
exchange, shock source, reflection transfer, cascade transfer, physical
dissipation, particle exchange, limiter correction, remap correction, and
closure residual. Reflection and internal cascade transfer cancel from the
total. Material positivity corrections increment the limiter diagnostic.

## Remap and restart

`RemapConservatively` overlaps old and new cells in physical arc length. It
conserves each integrated branch and, in spectral mode, every individual bin on
the covered domain. It increments the field-line generation and reports any
domain or roundoff correction.

The deterministic versioned checkpoint schema 2 contains source, representation,
operators, boundary objects, spectral grid, phase, epoch, field-line
generation, campaign seed, completed step, pending shock/particle exchange,
handoff epoch/checksum, provenance, accumulated ledger, geometry, and all
authoritative energies. Deserialization validates ownership and dimensions.

## Deterministic concurrency policy

`util/sep_reproducible_reduction.*` supplies worker-local append-only buffers.
Contributions are keyed by schema, source, field line, segment, branch,
spectral bin, species, stable particle, timestep, event, interval, and
purpose—never MPI rank, thread, queue order, or a pointer. Duplicate complete
keys are rejected. After workers join, the reducer sorts by the complete key
and accumulates with long-double intermediates. MPI rank blocks
must be gathered and passed through the same canonical reducer. `MPI_Allreduce`
is reserved for genuinely additive accumulators; non-additive authoritative
state uses owner/gather synchronization.

Warnings and limiter counts are integer atomics. Random streams are keyed by
campaign seed, particle ID, timestep, and purpose so optional diagnostics cannot
advance a mover stream. `EvidenceHashHex` hashes the canonical reduced result;
PAR01–PAR05 require the same hash across synthetic worker and rank layouts.

## Verification

```sh
make test-turbulence-core-unit
make test-reproducibility-unit
```

The first command executes registered `TURB02`–`TURB23` and `TURBOWN01` plus
CLI/source gates. `TURB21` compares a nonuniform periodic profile with exact
translated cell averages at two resolutions; `TURB22` compares the complete
time history with the analytical integral of a prescribed linear-in-time
wave-energy source. This tests source application and ledger accounting, not
the QLT calculation of that source from particles; the latter remains in the
native `growth_rate_validation_test.cpp` integration test. The
`TURB23` case then closes the weighted kinetic-energy change from controlled
outward and inward wave-frame scattering against opposite signed increments
in the resonant minus and plus wave branches. It requires both the total
particle-plus-wave energy and `particleExchangeJ` ledger to agree within
`2e-12 J`, with no positivity limiting. The
production catalog's established `TURB01` remains the independent 1-AU
magnetic-pressure closure. The second command executes PAR01–PAR05 plus
source-policy gates. Both compile the exact production
cores as C++11 with strict warnings and sanitizers. A complete AMPS checkout is
still required for native OpenMP/MPI evidence, SWMF handoff integration, and
long-run scientific conservation/observational campaigns.

The same controlled turbulence descriptors can be selected in a linked build:

```sh
../amps --test TURB21 --test TURB22 --test TURB23
../amps --test-group turbulence --test-json turbulence-results.json
```

Explicit group selection includes every registered turbulence case. Routine
`--all-tests` selection remains governed by each descriptor's runtime class.
# WP25 shock-source transaction

Shock injection no longer mutates `CellIntegratedWaveEnergy` before the common
driver. The shock adapter computes typed, provenance-bearing segment records
from swept spherical-shell overlap, per-line flux-tube volume, upstream mass
density, and the normal-relative speed. `PICAdapter::QueueShockContribution`
holds them until import maps them to `pendingShockPlusJ/MinusJ`. The common
`Turbulence::Advance` routine is the only writer and records the applied amount,
limiter correction, and closure residual in the signed energy ledger. Efficiency
and branch split come from the frozen WP30 `RunConfiguration`.
