# Finite-wave-energy collision design: CP05 / Q04

This is a derivation and implementation proposal, not an implemented collision
algorithm or a production qualification. The current four-spectrum owner,
resonant-collision (RC) kernel and prepare-only frozen-collision transaction
(FCT) retain their existing behavior. In particular, an FCT batch with an
insufficient final bin budget fails atomically. That useful bookkeeping gate
does not supply a conservative finite-reservoir stochastic evolution law.

The recommended target remains the adopted **stimulated, continuum quasilinear
system**. A finite number of weighted numerical particles is an approximation
to that target, not an independently specified physical finite-particle theory.
A candidate interacting jump approximation is derived below. Its depletion,
spectral-interface and weight limits must be established before implementation
in a native mover. No physical amplitude, dispersion or numerical acceptance
cutoff is selected here.

## 1. Fixed contracts and source boundaries

The stored state is magnetic energy `W_j >= 0` [J] in each physical cell,
propagation branch and circular helicity. Total wave energy at the adopted
small-amplitude nondispersive MHD order is `2 W_j`. The four components evolve
independently after their explicitly selected equal-L/R initialization. The
positive-k bin edges are authoritative. The collision reconstruction is

```text
P_h(k) = 2 mu0 W_j / (V Delta k_j),     k in bin j,
integral P_h dk = <delta B_h^2>.
```

There is no power outside the modeled support, no numerical bridge through the
ninety-degree resonance gap, and no additive spontaneous wave source in the
adopted collision law. A zero-energy bin therefore stays zero under stimulated
collision exchange alone. Explicit boundary reservoirs or another separately
specified physical operator can supply power; the collision algorithm cannot
manufacture that source.

The inspected implementation boundaries are:

| Source | Existing responsibility | What it does not establish |
|---|---|---|
| [`sep_resonant_collision.h`](../../util/sep_resonant_collision.h), `BuildGenerator` | Fixed-spectrum reversible FV pitch generator, true finite support, wave-frame invariant | A coupled generator when the same event changes its spectrum |
| [`sep_resonant_collision.cpp`](../../util/sep_resonant_collision.cpp), `ElasticProposal` / `ValidateProposal` | Weighted elastic energy/momentum increments and exact branch/helicity/bin path partition | Positivity of the updated finite wave budget or the event's stochastic selection |
| [`sep_frozen_collision_transaction.h`](../../util/sep_frozen_collision_transaction.h), `Prepare` | Deterministic private staging, encoded-particle energy audit, net bin budgets, canonical candidate | Exact coupled hazards, a native commit, or complete physical momentum accounting |
| [`sep_polarized_wave_state.h`](../../util/sep_polarized_wave_state.h), `PrepareEvolutionCandidate` | One immutable four-array authority, canonical identity/revision and structural validation | Physical qualification of a new update law |
| [`sep_polarized_wkb_transport.h`](../../util/sep_polarized_wkb_transport.h) | Separate controlled ray/work/boundary staging with external ledgers | An event clock coupled continuously to changing wave power |

The existing RC conductance is valid with the wave spectrum held fixed. Reusing
its equilibrium measure unchanged after eliminating wave energy along a
particle-energy constraint is generally incorrect.

## 2. Independent one-branch derivation

First freeze one cell's background, volume, magnetic direction and `v_A` during
a collision-only interval. This is an explicitly bounded local problem, not a
claim about an evolving material tube. Use plasma-frame time throughout. For
branch speed `u = sigma v_A`, define

```text
Gamma = 1 / sqrt(1 - u^2/c^2),
zeta = p_parallel' / p',
C = E - u p_parallel = E'/Gamma.
```

The primed total particle energy `E'` [J] and momentum magnitude `p'` [kg m/s]
are constant on this branch's elastic orbit. Rest energy is included in `E`
but cancels in every kinetic-energy increment. The exact Lorentz map gives

```text
E(zeta)          = Gamma (E' + u p' zeta),
p_parallel(zeta)= Gamma (p' zeta + u E'/c^2),
J(zeta)          = Gamma E(zeta) p'/c^2.
```

Thus the phase measure is `J dC dzeta`, and, for statistical weight `w` [number
of physical particles represented],

```text
A = Gamma u p',                         [J per unit zeta per particle]
Delta K_particle = w A Delta zeta,
Delta p_parallel,particle = w Gamma p' Delta zeta,
Delta K_particle = u Delta p_parallel,particle.
```

Within one resonant bin and away from its physical edges, the already adopted
RC normalization has the form

```text
D(zeta,W) = W d(zeta),
J D = pi mu0 q^2 W (1-zeta^2) / (V Delta k |zeta|).
```

Here `D` has units s^-1 and `d` has units (J s)^-1. This identity follows from
the documented real-one-helicity variance convention; no further factor of two
is introduced for kinetic wave energy. The fixed-wave backward generator is

```text
L_W g = W d g'' + b_fixed g',
b_fixed = (1/J) partial_zeta|W (J W d)
        = W [d' + d J'/J].
```

The fixed-W partial derivative is a necessary part of the physical definition.
It is not interchangeable with a derivative along a finite-energy constraint.

### 2.1 An energy-conserving stimulated numerical lift

For this local one-bin orbit let

```text
a = w A/2,                              [magnetic J per unit zeta]
H = w K(zeta) + 2 W,
W(zeta) = W0 - a (zeta-zeta0).
```

One candidate interacting-particle approximation is the Itô system

```text
dzeta = b_fixed(zeta,W) dt + sqrt(2 W d(zeta)) dB,
dW    = -a dzeta.
```

It conserves `H` pathwise in exact arithmetic because `K(zeta)` is affine.
Both drift and diffusion vanish at `W=0`; this boundary is absorbing for this
stimulated operator. In joint coordinates the same generator is

```text
V = partial_zeta - a partial_W,
L F = W d V^2 F + b_fixed V F.
```

For observables depending only on the particle, it recovers `L_W` at the
current joint state. This is a candidate numerical approximation, not a proof
that finite macroparticles have these physical fluctuation statistics.
Other branch collisions change `C` for this branch, so the orbit and its
coefficients must subsequently be reconstructed from the actual updated phase
point. One cannot keep either branch's old orbit frozen after the other acts.

### 2.2 A tempting but different generator

If a constrained implementation substitutes `W(zeta)` into the reversible RC
form and differentiates the product totally, it instead obtains

```text
b_total = (1/J) d/dzeta [J W(zeta) d]
        = b_fixed - a d.
```

At `W=0`, this additional drift gives `dW/dt = a^2 d`, even though there is no
stimulated power. It introduces a finite-weight wave source that was not
adopted. Detailed balance on a constrained energy surface is therefore not a
license to use the original fixed-wave stationary measure.

For the stimulated lift, the forward density `f` relative to `J dzeta` satisfies

```text
partial_t(J f) = partial_zeta [J W d partial_zeta f + W' J d f].
```

Its zero-current *interior* density is proportional to `1/W`, with speed
measure proportional to `J/W`; an absorbing boundary also requires its boundary
probability mass. This interior expression need not be normalizable at
depletion and is not an equilibrium prescription for the full problem.

These distinctions persist even when each proposed jump satisfies exact
energy conservation. Conservation alone does not identify the right drift.

## 3. Finite weights: numerical fluctuations and the mean-field target

For independent numerical particle noises, the variance rate of a shared bin
contains `2 sum_i a_i^2 W d_i`. It scales with squared particle weights. These
are macroparticle sampling fluctuations; their size is changed by splitting a
numerical particle into otherwise identical lower-weight particles. The
adopted mean-field QLT equations do not specify them as physical fluctuations.

The required limit holds represented physical density fixed while
`sum_i w_i^2 -> 0`, resolving the phase distribution and every active bin. It
must be studied separately from pitch mesh refinement and time splitting.
Increasing the wave reservoir, choosing a favorable seed, or only decreasing
an outer timestep is not a substitute for this weight convergence.

For smooth phase densities and no-flux endpoint terms, the expected energy
transfer of the lift is `-a integral b_fixed J f dzeta`. Integration by parts
turns this into `a integral J W d partial_zeta f dzeta`, the stimulated
wave-growth structure proportional to W. In a finite phase-population
discretization with bounded growth coefficients, a positive mean-field bin
therefore remains positive over finite time. Finite-weight sample paths may
instead hit zero. Convergence must distinguish this sampling-induced
extinction from deterministic physical wave damping, and account for finite-N
particle/wave correlations rather than assuming an exact factorization.

If a finite-number physical fluctuation theory, spontaneous emission or
thermal wave noise is desired, it needs its own governing source/correlation
law and normalization. It cannot be inferred from macroparticle weights. None
is proposed for the first stimulated qualification case.

## 4. Candidate budget-feasible jump approximation

The smooth one-bin problem admits a concrete local construction. At the
current physical state choose neighboring pitches `zeta+h_plus` and
`zeta-h_minus`, with both step lengths positive, **before constructing hazards
or drawing an event**. Choose targets within the current physically connected
support and such that every affected wave bin remains nonnegative. In the
one-bin case the draining step is at most `W/|a|`; a boundary node at exactly
zero wave energy can be included in the numerical state space.

Match the first two conditional infinitesimal moments:

```text
q_plus h_plus - q_minus h_minus = b_fixed,
q_plus h_plus^2 + q_minus h_minus^2 = 2 D,

q_plus  = (2 D + b_fixed h_minus) / [h_plus (h_plus+h_minus)],
q_minus = (2 D - b_fixed h_plus ) / [h_minus(h_plus+h_minus)].
```

Local refinement must make both rates nonnegative and reduce the consistency
remainder. This is an approximation to a diffusion: its mesh jump frequency is
not a unique physical collision-event frequency. At `W=0`, use the specified
absorbing state, not an undefined division by a zero draining step.

Reducing a *prospective numerical step* and recomputing moment-matched rates
is fundamentally different from clipping an already sampled physical transfer
or suppressing its rate after a failed budget check. No event is drawn from an
infeasible generator. A representability failure returns Unresolved before
publication; it is not repaired by moving the pitch, discarding an event,
changing its weight after selection, or drawing again.

This construction is only an interior candidate. The following work is still
required before describing it as a usable collision solver:

1. Prove weak convergence at the accessible depletion boundary, including any
   limiting boundary probability mass and mesh-rate growth.
2. Derive the weak transmission law at every mapped k/helicity edge. `D` is
   piecewise and its divergence-form interface drift cannot be replaced by a
   smooth derivative that ignores the jump. An arbitrarily tiny step that
   remains on one side is not a derivation of the correct interface crossing.
3. Preserve the exact RC partition when a path encounters more than one bin.
   Each traversed bin has its own constraint; a total branch budget cannot
   authorize borrowing from the wrong helicity or k bin. A positive-width
   zero-power interval is disconnected and cannot be crossed.
4. Establish consistency of both branch generators in a common physical phase
   space, including changed orbits and conditional rates after each event.
5. Verify conservation using the actual encoded momentum records and stored
   wave values, not only the ideal orbit increment. Depletion targets may be
   numerically unresolved at a given precision/weight; a typed failure is
   preferable to an unreported energy correction.

The present reversible RC conductances and J-weighted within-cell target draw
cannot simply be reused for these tasks. They solve the fixed-spectrum
problem, not the constrained problem derived in section 2.

## 5. Independent conservative continuum discretization

A useful reference and possible alternative representation evolves nonnegative
phase-cell populations `n_i` [physical-particle number] and wave bins together.
For an elementary phase face entirely inside one resonant bin, let `m_i,m_j`
be its exact phase measures and let conductance `G = kappa W`. Define

```text
F = kappa W (n_i/m_i - n_j/m_j),
dot n_i = -F,                  dot n_j = F,
dot W = -(K_j-K_i) F/2.
```

Mass and `K_i n_i + K_j n_j + 2 W` are exactly conserved by these semidiscrete
equations. The parallel momentum balance follows if the two phase states lie
on the same elastic orbit, `K_j-K_i = u(p_parallel,j-p_parallel,i)`. The vector
field points into the nonnegative population domain and vanishes at `W=0`.
For this one-face subsystem, `n_i+n_j=N` and energy eliminate W, leaving a
scalar product-of-linear-factors ODE for `n_j`. Its exact flow is a candidate
positive conservative elementary update. This is a derivation of the local
subproblem, not a qualified multi-face timestep integrator.

A composition of such elementary flows needs timestep/order convergence,
interface analysis and a common physical representation for both branches.
Independently exponentiating a wave growth rate while advancing a particle
solver does not retain the coupled energy invariant automatically.

This approach advances population mass rather than indivisible, fixed-weight
macroparticle jumps. A particle realization needs explicit cohort splitting,
weights, identity/lifecycle and conservative phase remapping. It must not be
presented as an exact sampled combined-branch CTMC. Its main immediate value
is an independently formulated continuum reference for the interacting jump
approximation. The conservative Galerkin work cited below supports using a
common weak energy-transfer form, but it does not supply a positivity proof
for either proposed algorithm here.

## 6. Event time, shared-bin ordering and outer operators

For a chosen semidiscrete joint generator, the exact collision-only clock is
the sum of all current outgoing rates. With per-particle clocks, each rate is
the sum over **both branches and all its feasible targets**. Expire integrated
total hazard first, then choose the actual branch/target using their current
rate fractions. Never choose a future branch at the preceding event.

An event changes the spectra seen by every particle sharing affected bins.
Their rates must be updated from that same newly accepted wave state. One
valid schedule retains every particle's residual optical depth, subtracts its
integrated old hazard over the elapsed common interval, then recomputes the
affected rates. A zero rate pauses its clock. Restart and outer-step cuts must
retain this residual, and simultaneous-event/tie handling must be explicit.
Alternatively a single cell-wide total-hazard clock can be defined, but that
is a different persisted lifecycle schema, not an implicit reuse of the
current particle clock. Stable physical cell/particle keys must determine
ordering and RNG domains rather than MPI rank or thread arrival.

This statement is exact for the selected piecewise-constant semidiscrete
collision generator. Material transport, WKB work, ray motion, boundary supply,
particle entry/exit and changing background generally change rates. Their
event-time coupling or controlled time splitting needs a declared chronology
and an outer-splitting convergence axis. Frozen-wave FCT proposals cannot be
retroactively interpreted as having these coupled hazards.

A batch whose final net budget is positive can still have an impossible
intermediate chronological deficit. Conversely, reordering credits before
debits changes a coupled stochastic path. Deterministic reduction order alone
does not establish physical event-time consistency.

## 7. Energy and momentum: distinguish the conserved quantities

For fixed u, the adopted transfer relation gives the algebraic invariant

```text
sum_i w_i p_parallel,i + sum_j 2 W_j/u = constant.
```

The wave term here is **canonical wave pseudomomentum** `k E_wave/omega`, not
automatically electromagnetic or material momentum. It provides a valuable
independent numerical audit, but must be labeled with that meaning.

In SI, an ideal parallel traveling perturbation has
`delta E = -u b cross delta B`. Its electromagnetic momentum is
`epsilon0 integral delta E cross delta B dV = (u/c^2) 2 W b`, at the retained
linear-wave order. This differs from `2 W b/u` by `u^2/c^2`. In nonrelativistic
MHD the electromagnetic momentum is neglected, and the background receives
the corresponding particle recoil. Keeping electromagnetic inertia instead
requires a consistently extended wave/background model; it cannot be added
to the adopted nondispersive closure by relabeling the canonical ledger.

Thomas & Pfrommer distinguish these wave and background momentum contributions
and discuss the order of the MHD approximation. Their result supports tracking
a prescribed-background impulse alongside energy/work, rather than asserting
that particle plus electromagnetic wave momentum alone is closed. See their
sections 5.1 and 6.2–6.3, equations 74 and 120–131.
([Primary paper](https://arxiv.org/html/1805.11092))

At the already adopted MHD order, a first implementation should therefore
record actual weighted particle momentum transfer, canonical wave transfer,
and the signed prescribed-background impulse separately. It must state the
retained order and excluded electromagnetic/relativistic corrections. Exact
arithmetic closure of the adopted energy ledger is not proof of an all-orders
relativistic physical stress-energy conservation law. The earlier domain study
still has to delimit this approximation; no new default v_A/c limit is chosen.

For the frozen local problem at the retained nonrelativistic MHD order, the
parallel background receipt is `Delta P_background =
-sum_i w_i Delta p_parallel,i b`, because electromagnetic wave momentum is
excluded at that order. This is a signed reaction ledger for the prescribed
background; it does not advance its velocity or infer a new background mass.
If electromagnetic momentum is retained, even this receipt changes by
`-Delta P_EM`, and energy/dispersion must be extended consistently. Neither
receipt should be added to the canonical pseudomomentum audit as if both wave
momenta represented the same physical conserved quantity.

For the full material case, momentum direction, boundary momentum flux and
external work also depend on the same staged geometry/background. A local
collision ledger cannot invent their global native mapping.

The outer transaction must report the signed identity
`Delta E_particles(actual) + 2 Delta W_B(actual) + Delta E_external +
Delta E_background_receipt = 0` to its declared numerical resolution.
Here an external/background receipt is positive for energy received by that
account; supplied energy is negative. For the bounded collision-only plasma
rest-frame problem, the adopted transfer has no additional background energy
receipt. That local statement does not eliminate the separate material
pressure-work, boundary or frame-transformation ledgers. Momentum receipts
and truncation terms must be reported alongside, not hidden in an energy
roundoff account.

## 8. Independent mathematical benchmarks and negative controls

No implementation gate is claimed in this document. The following predeclared
benchmarks are designed to discriminate the physical target and the wrong
drifts, not only to check a candidate against its own generator.

| Benchmark | Independent expected result | Defect it exposes |
|---|---|---|
| Stopped square-root normal form | Exact exit probability and mean time below | Missing absorbing mass, hidden reflection/spontaneous drift |
| Actual one-bin stopped orbit | Independent scale/speed quadrature | Fixed-W versus total derivative; phase-measure errors |
| Fixed reservoir | Continuum FP moments/measure computed independently | Wrong normalization, J factor, branch chart or time frame |
| Two nonempty bins | Weak flux/transmission reference from the original divergence-form PDE | Smoothing, edge trapping, arithmetic averaging or wrong bin partition |
| Empty band and 90-degree gap | No probability/energy flux through disconnected support | Artificial crossing or hidden floor |
| Two branches sharing waves | Exact finite-state master-equation exponential for a separately specified small graph | Wrong combined clock, premature mark, stale shared-bin rates |
| Finite weight sequence | Fixed represented distribution; weight-noise reduction and continuum moments | A macroparticle-dependent result called physical convergence |
| Encoded transfer | Independent stable relativistic energy and actual momentum differences | Ideal-proposal conservation masking stored-state loss |
| Chronological budget | Every accepted intermediate state feasible | Positive final net batch hiding an impossible path |

### 8.1 Exact stopped normal form and wrong-drift negative control

This is a dimensionally scaled mathematical normal form, not a complete
physical orbit or a claimed Lorentz diffusivity. On `0 < W < M`, stop the
process at either endpoint and take

```text
L phi = q W phi'',          q > 0 [J/s].
Pr(hit 0 before M | W0) = 1 - W0/M,
E[t_exit | W0] = (W0/q) ln(M/W0).
```

Both expressions follow directly by solving `L h=0` and `L T=-1` with their
endpoint values. For `M=1 J`, `q=1 J/s`, `W0=0.25 J`, the independently exact
targets are probability `3/4` and mean time `ln(4)/4 s`. This stopped domain
tests the numerical process without prescribing a production upper boundary.

For locally constant d and J the stimulated constrained lift reduces to this
form with `q=a^2 d`. The incorrect total derivative adds `+q phi'`; at zero it
produces positive wave drift, and `h=1-W/M` then has residual `-q/M`, not zero.
Its zero boundary has a different accessibility classification. A conservation
test alone would accept both algorithms, so this is an essential negative
control.

### 8.2 Actual source coefficients without a fabricated constant diffusion

On a stopped interval `[l,r]` lying inside one actual resonant bin, choose an
explicit W(zeta) positive in the interior and vanishing at the draining end.
For the correct stimulated lift its scale derivative is

```text
s'(zeta) = exp[-integral b_fixed/D dzeta] = constant / [J(zeta) d(zeta)].
```

The probability of hitting r before l is
`[s(zeta)-s(l)]/[s(r)-s(l)]`, computed by independent high-precision quadrature.
Expected exit time follows from the scale/speed Green function, equivalently
the independently solved equation `D T''+b_fixed T'=-1`, `T(l)=T(r)=0`.
Endpoint integrability/accessibility must be checked explicitly for each
declared case. The wrong total derivative introduces an extra `1/W` into the
scale derivative. No fitted reference coefficient or collision-kernel call
should generate these expected values.

There is also an analytic check using the actual source coefficient. Let
`A0=pi mu0 q^2/(V Delta k)`, so on a positive-pitch one-bin interval
`J d=A0(1-zeta^2)/zeta`. Then

```text
s(zeta) = -ln(1-zeta^2)/(2 A0) + constant.
```

For the explicitly controlled stopped coordinates `l=1/4`, `zeta0=1/2`,
`r=3/4`, choose the along-B branch with `a>0`, `W0=a(r-zeta0)` and one bin
whose support contains this entire interval. These are mathematical test
coordinates, not production parameters or new physical limits. The exact
hit-r probability is

```text
ln(5/4) / ln(15/7) = 0.2927854935781266... .
```

For example, bin edges `[k_star,8 k_star]`, where
`k_star=|q| B Gamma/p'`, place both stopped endpoints strictly inside that
bin's resonance interval. This probability is independent of the coefficient
normalization and initial weight; the corresponding exit-time reference must
also be checked to test normalization and timing. The formula is derived from
the scale equation, not a call to RC quadrature.

## 9. Convergence and implementation sequence

Record each axis separately; agreement on one does not qualify the others:

1. One-bin pitch/joint-state mesh, depletion probability and exit-time error.
2. Sharp physical bin-interface treatment, with all real edges retained.
3. Particle weight/distribution resolution at fixed physical normalization.
4. Coupled event chronology and hazard integration, including time partition,
   paused clocks, both branches and affected shared-cell particles.
5. Collision/transport/WKB/boundary splitting and material geometry resolution.
6. Floating-point transfer resolution, precision and deterministic reductions.
7. Spectral discretization inside declared fixed support; physical support is
   never enlarged as a numerical gap-bridging device.

Ensemble uncertainty, deterministic discretization error and physical-closure
error need separate reports. Seeds, numerical controls, represented weights,
initial four spectra, physical inputs, unresolved cases and source identities
must be retained. Error budgets for acceptance require review; none is fitted
to obtain a favorable result.

Recommended bounded work after review is: (a) independent stopped benchmarks;
(b) one-bin constrained numerical kernel with no native publication; (c) sharp
two-bin weak-interface derivation/reference; (d) two-branch/shared-bin event
generator and weight convergence; (e) outer material/WKB chronology; and only
then (f) a single native transaction staging particle phase, all RNG/clock
state, wave owner, geometry identities and external ledgers before an atomic
publication. Existing owner/checkpoint/source lifecycle gates remain required.

## 10. Decisions and unresolved scope

**User decisions, 2026-11-28.** These select the design; they implement
nothing and qualify nothing.

1. **Coupling algorithm.** Production keeps weighted native particles and
   uses the constrained interacting jump approximation of section 4. The
   phase-population finite-volume system of section 5 is built only as its
   independent conservative continuum reference, not as a production
   representation. The five proofs listed in section 4 remain required
   before the jump approximation may serve as a collision solver.
2. **Finite budgets (Q04).** The target is the stimulated continuum QLT
   system. A bin emptied by finite macroparticle weights is a numerical
   sampling effect, studied with `sum_i w_i^2 -> 0` at fixed represented
   physical density (section 3). The generator chooses budget-feasible
   targets before any hazard is built, and `W=0` is absorbing (section 4).
   Event time uses per-particle clocks that keep their residual optical
   depth; every particle sharing an affected bin is re-rated from the
   accepted wave state (section 6). No finite-population fluctuation or
   spontaneous-emission law is adopted. Still open: what the run does after
   a typed Unresolved result, either stop or split the particle before a
   new attempt.
3. **Background of the first coupled case.** The approved Stage 12 SWCME1D
   `SHOCK_ONLY` material profile; see
   [the contract review](WAVE_EVOLUTION_CONTRACT_REVIEW.md).
4. **Retained two-branch model.** A regression baseline only (CV10–CV12
   audited and run routinely, IV06 rewritten with a computed resonance or
   retired). It is not a production closure for self-excited waves.

The first implementation step is section 9 (a): the independent stopped
benchmarks and negative controls of section 8.

Already fixed and **not** questions to reopen: magnetic-only storage and
`E_wave=2 W`, four independent spectra, the real-wave QLT normalization,
stimulated zero-power behavior, finite support/gaps, initial reservoirs,
material-tube kinematics, and no reflection/cascade for the first case.

The recommended numerical interpretation is continuum QLT with controlled
weight and mesh limits. The source equations do not identify a unique physical
finite-number jump process. If exact finite-physical-population statistics are
requested instead, the missing input is a governing finite-population wave
fluctuation/spontaneous-emission theory, not a depleted-bin rejection policy.

The algorithmic choice for review is whether to retain weighted native
particles with a convergent constrained interacting jump approximation, or
allow a phase-population/cohort FV representation. The first fits the existing
particle architecture but needs the boundary/interface and weight proofs
above; the second offers a direct conservative mean-field reference but needs
new representation/remap/lifecycle contracts. This is an explicit numerical
design choice, not permission to select new physical sources.

Before a claim of physical momentum conservation, adopt the clearly labeled
background impulse ledger at the already selected MHD order or supply a
consistently extended stress-energy model. No total-derivative drift,
microcanonical detailed balance, wave quantum size, seed power or finite-N
thermal noise should be inferred to fill these gaps.

## Primary references and limits of their use

- [Thomas & Pfrommer, 2019, MNRAS 485, 2977](https://arxiv.org/html/1805.11092):
  adopted small-amplitude Alfvén energy/force/background framework. The one-bin
  constrained drift and jump proposal here are independent derivations, not an
  algorithm claimed to appear in that paper.
- [Huang, Abdelmalik & Gamba, *A Conservative Galerkin Solver for Quasilinear
  Diffusion Equation with Arbitrary Wave Dispersion Relation*](https://web.ma.utexas.edu/users/gamba/papers/HAG-2022.pdf):
  equations 2.1, 2.10–2.11 exhibit the stimulated coupled structure, canonical
  parallel momentum and common weak energy-transfer form. Its dispersion and
  particle problem are not substituted for the adopted SEP coefficients, and
  no positivity/depletion result for this proposed algorithm is attributed to
  it.
- [Afanasiev & Vainio, 2013, ApJS 207, 29](https://doi.org/10.1088/0067-0049/207/2/29):
  existing adopted transport reference. This investigation did not obtain the
  full primary text from that endpoint, so it does not attribute a particular
  finite-budget or finite-population rule to that paper.
