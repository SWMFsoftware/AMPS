# Alfvén-turbulence ownership and evolution contract

Steps 11 and 12 retain self-consistent Alfvén turbulence as a physics
subsystem; it is not a fourth particle mover. The exactly three production
movers remain `parker`, `fte-dmumu`, and `fte-mfp`.

The adopted 2026-11-16 physical contract stores magnetic energy alone, with
E+ along B and E− opposite B. Linear small-amplitude MHD gives an equal kinetic
wave energy, so total wave energy is twice the stored magnetic energy. A
particle energy gain therefore debits half that amount from magnetic storage.
The selected resonance-resolved implementation must evolve separate left/right
polarizations in each branch. The user selected equal initial L/R energies in
each branch, using the existing wave-energy initializer; all four spectra then
evolve independently. Existing logarithmic k values are bin centers, with edges
half a logarithmic spacing to either side. The finite-band zero-scattering gaps
are retained and lambda must derive from the same law.
The replacement coupling is under development; retained legacy growth/source
adapters and scalar ledger tests do not qualify it. Domain limits remain unset
pending `validation/studies/weak_alfven_domain.json`.

The controlled real-wave regression `RWN01` integrates the relativistic Lorentz
equations independently of the collision operator. With the explicitly real
circular field `dBx=a*cos(phi), dBy=-a*sin(phi)`, the positive-k spectrum of
that helicity satisfies `integral P_h dk = <dB_h^2>`. Its one-helicity coefficient
is `D_mumu = (pi/2)*Omega^2*(1-mu^2)*(1-u*mu/v)^2*P_h /
(B0^2*abs(v*mu-u))`, where `Omega` is the relativistic gyrofrequency and
`u=sigma*v_A`. This convention removes the extra factor of two formerly in
the neutral collision kernel. `RWN02` independently checks the balanced real
L/R limit: the prescribed Jokipii coefficient multiplying a *total*-variance
spectrum needs `pi/4`. That balanced prescribed limit does not authorize
replacing four independently evolved spectra with an imposed equal split.

The four-spectrum value owner and collision kernel are being connected to PIC.
Until the joint update, domain study and background ledger are qualified,
production preflight rejects dependent scattering/local wave evolution.
Component-test processes explicitly select a controlled diagnostic purpose;
legacy two-branch operator tests describe software behavior, not qualification
of the adopted physical closure.

The next finite-wave-energy collision design is recorded in
[the CP05/Q04 derivation](validation/physics/FINITE_WAVE_ENERGY_COLLISION_DESIGN.md).
It distinguishes the fixed-wave QLT drift from a total derivative along a
finite-energy constraint: the latter adds an unadopted source at zero wave
energy. A budget-feasible interacting jump approximation and an independent
conservative continuum reference are proposals requiring depletion/interface,
event-clock and weight-convergence proofs. Existing frozen-batch budget
rejection remains a prepare-only bookkeeping gate. Canonical wave
pseudomomentum and physical prescribed-background impulse require separate
labels; this design adds no production update or qualification claim.

The user selected on 2026-11-28 (design section 10): weighted native particles
with the constrained interacting jump approximation, the phase-population
finite-volume system only as its independent reference, the continuum QLT
target with per-particle residual-hazard clocks, and the Stage 12 SWCME1D
`SHOCK_ONLY` material profile as the first coupled-case background. The
retained two-branch model below is a regression baseline only: its tests
describe software behavior and are not run as qualification of self-excited
waves. The next step is the independent stopped benchmarks of design section 8.

## Four-spectrum WKB boundary and material-tube contract

The first qualification case uses a **material magnetic tube**: vertices move
with the full plasma velocity U. Particle coordinate convection must therefore
use U minus grid velocity, which vanishes for this motion. Moving the vertices
while retaining an additional U convection would count the same motion twice.
This selection does not establish that a retained static Parker/provider state
obeys ideal induction, continuity, or the required material geometry.

On every incoming signed characteristic, the spatial reservoirs continuously
supply the **initial magnetic energy density [J/m³ per bin]**, separately for
all four components. Outgoing characteristics leave freely. Classification
uses `(U + sigma*v_A*b - v_grid) dot n_out`; no branch is assigned permanently
to the outer boundary. Initial reservoirs are captured from revision1 and must
be retained across evolution/restart rather than reconstructed from the
evolved interior. A versioned content identity detects accidental mutation of
the public reservoir value; trusted initial capture and immutable outer
ownership still establish its historical origin. Reflection and cascade are
disabled for this first case.

The modeled k edges remain fixed. Exterior incoming power is zero. Energy and
total-wave action leaving either spectral edge are recorded at the crossing
time in an external account; they are not heat, reflection, or a numerical
floor. Finite-band resonant gaps therefore remain physical gaps of this model.

`util/sep_polarized_wkb_transport.*` is a **prepare-only controlled component**.
It transports every bin using the same flat-per-linear-k `Ebin/Delta k`
reconstruction used by the collision view. Ideal-MHD ray/work coefficients are
`d ln k/dt = -b.grad(U).b - sigma*d_s(v_A)` and magnetic packet work
`d ln E_B/dt = -div(U)/2`. All velocity, deformation, cap flux, mass, lateral
sweep, and geometric inputs are explicit; discrete GCL, magnetic flux, and
source-free mass continuity are checked before a candidate is returned. The
current stage requires fixed B direction within a resolved interval. Affine
manufactured maps supply explicit background kinematics; the component does
not solve background momentum/thermal equations or qualify a provider's
forcing law. It neither advances native vertices nor publishes a native wave
state.

The mechanical ledger includes physical `U*P_wave` pressure work on **all**
surfaces, including material side walls even when their relative advective flux
is zero. `P_wave=W_B` and total wave energy is `2*E_B`. Pressure is frozen at the
donor stage; that quadrature, spatial/source splitting, and action projection
need convergence studies. The explicit action residual is numerical error,
never assigned to background heating. External nonnegative spectral accounts
use exact binary64 integer summation with a declared finite-total bound;
displayed nearest/remainder values carry their own rounding bound.

The shared callback family `PolarizedWkbDescriptors()` and dependency-light
runner `test/run_polarized_wkb_tests.sh` exercise signed boundaries, ideal
expansion/strain, finite support, spatial/ray convergence, arbitrary controlled
ALE GCL, rejected inputs, material-side pressure work, and multiscale external
accounting. These are synthetic equation/software verification, not native or
observational qualification. Native material staging, coupled transactions,
all-three-mover integration, full checkpoint state, and physical domain studies
remain required. See
[the material/WKB design](validation/physics/MATERIAL_TUBE_WKB_STAGING.md) and
[the governing-contract review](validation/physics/WAVE_EVOLUTION_CONTRACT_REVIEW.md).

The separate neutral material kinematic component
`util/sep_material_tube_kinematics.*` checks local Cauchy induction/mass/Nanson
area mapping and global-affine cell GCL with actual curved reference-surface
moments, explicit time quadrature and numerical budgets. `MTG01`–`MTG06` use
independent affine/curved/parabolic references and atomic rejection controls.
It publishes no native geometry; arbitrary nonaffine cell interpolation and
the full particle/wave/background material transaction remain unavailable.
See [its bounded contract](validation/physics/MATERIAL_KINEMATIC_CORE.md).

The finite-amplitude study now compares independent relativistic Lorentz
orbits with a derived finite-time unperturbed response in prescribed circular
waves. Its16 bounded mathematical/protocol checks pass. The initial unchanged
pilot retains9 numerically resolved and27 unresolved responses; it selects no
physical amplitude/k cutoff. A separate versioned time/phase refinement using
the exact parent execution code and unchanged acceptance criteria retains31
resolved and5 phase-unresolved responses; the original data remain unchanged.
The five remaining phase-variance triples and numerical allowances are
reported before any further refinement. Coherent single-wave or separate-mode ensemble
variance is not presented as Markov diffusion. See
[the controlled response study](validation/physics/FINITE_AMPLITUDE_RESPONSE_STUDY.md)
for inputs, separate numerical estimates, limited results and unresolved scope.

## Retained two-branch source and authoritative state

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

## Retained two-branch driver order and stability

`SEP::Turbulence::Advance` is the context-neutral driver for standalone and
coupled state adapters. Its fixed order is:

1. shock and other external wave sources;
2. particle-wave energy exchange;
3. finite-volume advection with explicit CFL subcycling;
4. branch reflection;
5. cascade and physical dissipation;
6. authoritative/derived representation synchronization; and
7. diagnostics, ledger accumulation, and checkpoint phase reset.

`SEP::Turbulence::PICAdapter::Advance` is called once by the common
`amps_time_step()` in `main_lib.cpp`; standalone `main.cpp` calls that shared
driver. It flushes typed particle
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

The following operators describe the retained software model; they do not
qualify reflection, cascade, or limited particle exchange for the adopted
four-spectrum closure. Reflection is an internal signed branch transfer and
preserves total energy.
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

## Retained two-branch boundary and conservation ledger

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

The pure-core deterministic checkpoint schema 3 contains source, representation,
operators, boundary objects, spectral grid, phase, epoch, field-line
generation, campaign seed, completed step, pending shock/particle exchange,
handoff epoch/checksum, provenance, accumulated ledger, geometry, and all
authoritative energies, including magnetic polarity relative to the line.
Deserialization validates ownership and dimensions. This codec does not save
the actual attached native PIC particle lists. The separate application physics
checkpoint is under native qualification and must capture live arrays, private
boundary/source caches, particle/source random state and snapshot/clock state.

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
