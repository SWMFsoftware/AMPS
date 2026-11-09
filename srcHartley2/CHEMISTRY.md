# Minor-daughter chemistry statistics roadmap

## Purpose and scope

This roadmap describes how to improve the statistical quality of chemically
produced minor species in the generic AMPS comet application without
increasing the number of directly injected parent particles. C2 is the
immediate diagnostic species, but every design and input control must operate
on arbitrary reaction networks and species lists.

The work belongs on `1d_comet_dev`. Nothing in this roadmap authorizes a
switch to, commit on, or merge into `master`. Each stage must be implemented,
tested, documented, and reviewed on the development branch. Stop and present
the complete validation result to the user before any operation that changes
`master`.

This is a staged plan, not a statement that the described features already
exist. In particular, `ParticleWeightCorrectionMode=on` cannot simply be
selected for the current radial model: the local-timestep mover and radial
adaptation validator explicitly reject that combination.

## Measured baseline

The reference archive for this roadmap is:

```text
validation_runs/halley1d_16x_particles_40000iter_sampling10000_20261104_18mpi
```

Its 18-rank production run reached 40,000 iterations and exited normally, but
did not converge. The final convergence record was a maximum relative L1
change of `3.733287e-01`, limited by CH, with 0/3 consecutive passing windows
at a tolerance of `5.0e-2`. Therefore its profiles are diagnostic evidence,
not an accepted steady solution.

The final 10,000-step sample reports the following mean instantaneous model
particle occupancies:

| Species | Mean particles/cell | Maximum | Fraction of zero cells |
| --- | ---: | ---: | ---: |
| C2 | 0.140 | 1.559 | 71.6% |
| C2H | 0.073 | 0.594 | 56.1% |
| CH | 0.0115 | 0.119 | 47.9% |
| CN | 0.108 | 0.569 | 45.4% |
| C2H5 | 0.0716 | 0.394 | 40.7% |
| HCO | 0.155 | 0.793 | 17.1% |
| H2O | 231.8 | 1,036 | 0% |

For C2, the first nonzero final-window cell occurs at 6,384 km. At 10,000 km
the AMPS density is `0.681 cm^-3`, compared with `81.45 cm^-3` in the
digitized Rubin profile. The AMPS values at that radius in the last three
windows are `0.987`, `1.633`, and `0.681 cm^-3`. The variation and the zero
region show that the present run cannot distinguish a physical C2 deficit
from transient and Monte Carlo error.

The relevant final carrier rates are:

```text
C2H2 direct parent: W/dt = 6.5625e25 molecules/s
C2 daughter:        W/dt = 5.46875e27 molecules/s
C2H daughter:       W/dt = 5.46875e27 molecules/s
```

Thus one selected C2H2-to-C2 channel has expected fixed-weight C2
multiplicity

```text
mu = (W_parent/dt_parent)/(W_C2/dt_C2) = 0.012.
```

Approximately one in 83 selected channel events creates a model C2 particle.
Residual sampling is unbiased over an ensemble, but this carrier ratio makes
the variance unacceptable for a minor-species profile.

## Representation and invariants

Use the following notation throughout implementation and review:

- `W_s,b`: base statistical weight of species `s` in block `b`;
- `dt_s,b`: local timestep of species `s` in block `b`;
- `B_s,b = W_s,b/dt_s,b`: base flux carrier;
- `c_p`: optional individual correction on particle `p`; and
- `C_p = B_s,b*c_p`: physical flux carrier represented by `p`.

Every stage must preserve these invariants:

1. Direct physical source rates, reaction rates, branching fractions, cross
   sections, and excess energies are not altered to obtain a smoother plot.
2. Number, mass, momentum, kinetic energy, and internal-energy production are
   conserved exactly where possible and otherwise in a documented unbiased
   expectation with a measurable variance.
3. A particle crossing a block face carries the same expected physical flux
   into the destination as it carried out of the source.
4. A non-analog probability is always accompanied by its complete likelihood
   ratio. Biasing a reaction must not shorten the parent lifetime.
5. Empty cells and blocks are not seeded merely to satisfy a population
   target.
6. No reaction transaction, face transfer, adaptation, or MPI exchange may
   exceed the configured particle-buffer safety limits.
7. Direct-source parent populations can be locked. Improving daughters must
   not implicitly increase their injection counts.
8. Results created immediately by cloning are correlated. A burn-in or a
   measured decorrelation interval is required before those particles count
   toward a science sample.
9. Default behavior for applications that do not select a new control remains
   unchanged.
10. A Rubin curve is validation data, never a target used to tune a reaction
    probability or statistical weight.

## Stage 0: freeze the baseline and add observability

### Objective

Make every later statistical claim reproducible and expose where minor
species are lost: channel selection, residual roulette, transport, chemistry,
or escape.

### Implementation

1. Preserve the exact baseline input, executable hash, run log, comparison
   manifest, and final output-window index in the validation archive.
2. Add a read-only minor-species analysis tool under
   `srcHartley2/visualization/` that writes CSV and JSON summaries containing:
   - mean, median, quantiles, and maximum sampled particle occupancy;
   - zero-cell fraction and first/last nonzero radii;
   - base `W`, local `dt`, and `W/dt` ranges;
   - window-to-window density changes at requested radii;
   - integrated density changes using the same cell measures as convergence;
   - comparison values and ratios only where both AMPS and reference data
     exist; and
   - explicit convergence and sampling-window metadata.
3. Add chemistry counters, reduced collectively at each complete sampling
   break, for each reaction and product:
   - eligible parent visits;
   - physical reaction selections;
   - selected channels;
   - product multiplicity zero/one/more-than-one counts;
   - allocated and reused slots;
   - represented product number, mass, momentum, and energy;
   - destruction and neutral sink counts; and
   - allocation deferrals and emergency requests.
4. Add species production, destruction, inner-boundary, and outer-boundary
   flux tallies so steady-state closure can be evaluated independently of
   density convergence.
5. Record all diagnostic definitions and units in `srcHartley2/README.md`.

### Tests

- Unit-test counter accumulation with reactions having one, two, repeated,
  and zero products.
- Verify MPI reductions at 1, 4, and 12 ranks using identical physical event
  totals, allowing stochastic differences only where expected.
- Re-analyze the preserved final window and reproduce the baseline C2 values
  above within output rounding.
- Confirm that diagnostics disabled at input or compile time add no fields and
  do not change model results.

### Exit criteria

- A machine-readable baseline report is archived.
- Every physical C2 creation path can be reconciled from selected parent
  channel through allocated product, later loss, and escape.
- Diagnostic overhead is measured and reported; it is not guessed.

### Stage 0 closure implementation and evidence

The focused test entry point is:

```sh
srcHartley2/tests/run_stage0_tests.sh
```

It provides both focused tests and an optional production-report gate:

- `test_stage0_chemistry_diagnostics.cpp` exercises one-product, two-product,
  repeated-product-term, and zero-product/sink reactions. It separately
  checks eligible visits, selected reactions/channels, destruction/sink,
  zero/one/many carrier multiplicity, slot allocation/reuse, deferral,
  emergency, and represented number/mass/momentum/energy.
- The same deterministic physical event set is partitioned over 1, 4, and 12
  MPI ranks. Each run must reduce to the same count and physical-moment
  signature; no stochastic tolerance is hidden in this accounting test.
- Disabled diagnostics must retain no reaction table and serialize no count or
  moment fields. This is the unit-level acceptance condition for the future
  input/compile-time off path.
- `test_stage0_minor_species.py` verifies exact cell-measure-weighted L1
  changes, per-radius changes, conditional comparison ratios, convergence-log
  metadata, and explicit rejection of invented measures for legacy output.
- The external high-statistics Halley archive is checksum-verified and
  re-analyzed. C2 occupancy, populated-radius range, `W/dt`, density at
  10,000 km, and the Rubin ratio must reproduce the committed machine-readable
  baseline in
  `srcHartley2/halley2011/reference/halley1d_stage0_baseline.json`.
- `CometDiagnostics=ON|OFF` is parsed by the application and defaults to
  `OFF`. With it enabled, `CometChemistry.cpp` records a parent visit once per
  positive-rate channel and records expected and represented product number,
  mass, momentum, energy, and represented-number variance without adding or
  reordering random draws.
- `CometDiagnostics.cpp` collectively emits one CSV and JSON record per full
  sampling window. Species fluxes are represented physical rate contributions
  accumulated as `W/dt` and divided by the sampling-window length. Gross inner
  incident and re-emitted rates are separate even though the current
  thermalizing boundary makes their net contribution zero.
- `RadialSimulation.cpp` writes the exact AMR cell traversal and `Measure()`
  values to `PT/diagnostics/cell_measures.csv`; no shell-volume reconstruction
  is used.
- `test_stage0_production_report.py` independently checks CSV alignment,
  reaction/product closure, flux identities, the two C2 creation channels,
  the C2 destruction channel, exact cell measures, and the approved
  represented-production acceptance
  `abs(represented-expected) <= max(5 sigma, 0.01 abs(expected))`.

Stage 0 is complete on `1d_comet_dev`. The clean 14-rank run at
`validation_runs/stage0_closure_20261105/diagnostics_on_10000iter_14mpi/production`
generated 20 reports, 2,984 exact cell measures, and 540 model output files.
It observed both selected C2-producing channels and a represented C2 carrier;
the exposed C2 destruction and outer-escape terms were physically zero during
this finite transient. Every output file was byte-identical to the preserved
diagnostics-disabled 10,000-iteration baseline. A controlled 4-rank paired
smoke test measured 8.58 s disabled versus 8.66 s enabled (+0.93%) and
1,376,232 KiB versus 1,376,572 KiB maximum resident memory (+340 KiB). The
closure run itself did not converge (0/3 required windows), so it validates
observability and non-perturbation, not a steady Halley solution.

## Stage 1: input-selectable population policies

### Objective

Allow direct parents and chemical-only daughters to have different numerical
population policies while keeping the application generic.

### Implemented selector model and input spelling

The user approved supporting both explicit-species and derived-group
selectors. Precedence is fixed as explicit species, then derived group, then
global default. Stage 1 implements the following exact spelling:

```text
RadialParticleNumberPerCellRange=800 | 3200
RadialParticleNumberPerCellRange(CHEMICAL_ONLY)=50 | 200
RadialParticleWeightLock=DIRECT_SOURCE
RadialParticleWeightAdaptationInterval=1000
```

Selector precedence must be deterministic:

1. explicit species;
2. reaction/source-derived group; and
3. legacy global range.

`DIRECT_SOURCE` means a species with a positive configured source rate.
`CHEMICAL_ONLY` means a species reachable through enabled reactions but with
no direct source. Group membership is derived only after the generated source
and reaction tables are available. An explicit species override must be
accepted for mixed cases.

The existing global record remains backward compatible. If no group or
species records are supplied, generated C++ and runtime behavior must match
the current case.

### Implementation

- Extend `input/hartley2.pl` to validate selectors, ranges, precedence,
  duplicate records, unknown species, and contradictory lock/range records.
- Generate species-indexed policy tables rather than embedding comet names.
- Separate the adaptation cadence from `RadialSamplingLength`. Adapt only at
  a completed mover/chemistry boundary, never while continuation particles
  are in flight.
- Report the resolved policy for every species during startup.

### Tests

- Extend the Perl parser/generator suites with valid precedence cases and
  failures for duplicates, unknown selectors, invalid ranges, and locked
  species with an active range.
- Hash generated output for a legacy input before and after the parser change.
- Verify that changing a daughter range does not modify direct-source
  injection weights or rates.

### Exit criteria

- Policies are fully input-driven and comet-independent.
- Legacy inputs are unchanged.
- A startup table proves which species are locked and which target applies.

### Completion evidence (2026-11-05)

Stage 1 is complete on `1d_comet_dev`. `CometPopulationPolicy.pm` resolves
global, `DIRECT_SOURCE`, `CHEMICAL_ONLY`, and explicit-species records after
the source and arbitrary reaction tables are available. The parser rejects
duplicate selectors, unknown species/groups, invalid ranges, and a range and
lock assigned at the same precedence level. Reaction reachability is the
transitive closure of enabled, positive-rate channels; therefore a secondary
or later-generation daughter is classified without embedding a Halley species
list.

The runtime now consumes generated species-indexed ranges, locks, policy
sources, and derived classes. A positive
`RadialParticleWeightAdaptationInterval` schedules adaptation at completed
iteration boundaries independently of sampling. Omitting it selects the
legacy complete-sampling-break cadence. Regular population and
chemistry-allocation planning cannot modify a locked species. Emergency
common-factor coarsening remains the explicitly reported buffer-safety
override because refusing it could prevent a rank from reaching the next safe
MPI exchange.

The validation archive is
`validation_runs/stage1_population_policy_20261105`. Its canonical four-rank
validator resolves all 27 species and passes owner/ghost, face-transfer,
photochemical-allocation, and emergency checks. A second build changes only
the chemical-only range from 50--200 to 100--400 and proves every direct
source rate and `W/dt` carrier unchanged. A cadence run records adaptation at
iterations 0, 50, and 100. Finally, an exact archived legacy input produces
54 output `.dat` files byte-identical to the pre-Stage-1 baseline. Run all
available gates with the environment-variable command documented in
`srcHartley2/README.md`.

## Stage 2: component-wide, correction-free daughter carriers

### Objective

Make the existing fixed-weight method useful for minor daughters before
introducing individual corrections.

### Current defect to replace

The present controller selects a population-driven weight independently in
each block. `RadialMaxFaceTransferMultiplicity=1` subsequently selects the
largest carrier in each face-connected component and raises every smaller
carrier. An empty or coarse block can therefore undo daughter splitting over
the complete radial domain.

### Algorithm

For every species and face-connected component:

1. Gather authoritative owner counts, cell counts, carriers, direct-source
   locks, and buffer headroom.
2. Compute the component population statistic required by the resolved input
   policy. Empty cells contribute zero; they are never seeded.
3. Select one dyadic carrier multiplier for the entire component. Apply the
   same `W/dt` multiplier to populated and empty blocks.
4. Split or roulette only particles already present. Changing the carrier of
   an empty block changes no population.
5. Solve the reaction-network allocation constraints from locked parents
   through all reachable daughters. Raise daughter carriers only as needed to
   satisfy `RadialMaxPhotochemicalProductsPerReaction`.
6. Re-run the component equality constraint after chemistry planning and
   iterate to a fixed point.
7. Precompute the exact worst-case number of clones and reaction allocations
   on every rank. Commit the transaction only if all ranks pass capacity
   checks.
8. Synchronize owner values to ghosts and verify both senses of every face
   before transport resumes.

Applying one carrier to the whole connected component retains face
multiplicity one while allowing the component as a whole to lower a daughter
weight. This removes the present conflict between local population planning
and face safety.

### Provisional sizing study

Using the completed final sample and an illustrative target of 50 mean
daughter particles per cell produced the following offline estimate:

| Species | Carrier reduction | Predicted mean occupancy |
| --- | ---: | ---: |
| C2 | 512 | 71.7 |
| C2H | 1,024 | 75.2 |
| CN | 512 | 55.4 |
| C2H5 | 1,024 | 73.3 |
| HCO | 512 | 79.4 |

All estimated reaction allocations remained below the existing hard limit of
20; the maximum was 17. The mean live-population proxy increased from about
10.3 million to 12.7 million, approximately 24%, with direct parents locked.
These values size a pilot only. They are not defaults or validation results.

### Tests

- Pure planner tests on a one-block component, a long radial chain,
  disconnected components, and DIM2/DIM3 synthetic face graphs.
- Test empty components, one populated block among empty blocks, and
  populations on both ends of a long chain.
- Validate power-of-two selection at exact lower and upper boundaries.
- Validate every reachable reaction, including repeated products and sinks,
  against the hard allocation limit.
- Repeat the integrated adaptation suite at 1, 4, and 12 ranks.
- Verify identical owner/ghost carriers after synchronization.
- Verify exact split conservation and ensemble roulette conservation for
  number, mass, momentum, kinetic energy, and internal energy.
- Force a capacity-preflight failure and prove that no partial carrier or list
  change is committed.

### Pilot acceptance criteria

- Direct-source injection counts and carriers are bitwise unchanged.
- Every face reports multiplicity at most one within tolerance.
- Every reaction reports worst-case allocation at most the configured limit.
- C2, C2H, CN, C2H5, and HCO achieve their configured component target after
  burn-in, or the report identifies the exact chemistry constraint preventing
  it.
- No emergency adaptation, buffer exhaustion, or partial transaction occurs.

### Implementation and transaction-test evidence (2026-11-05)

The Stage-2 component transaction is implemented on `1d_comet_dev`. The
runtime builds a deterministic face graph from the replicated AMR leaf tree,
collects exactly one authoritative count/carrier row from each block owner,
and calls the geometry-independent planner once per species. Population
selection uses total particles divided by total component cells; an empty
block therefore affects the requested mean but is never seeded. All targets
remain on each species' power-of-two carrier lattice.

Chemistry is solved on the component scalars before mutation. Reachability is
closed transitively from represented and configured direct-source parents,
repeated product terms are aggregated, irreducible stoichiometric floors are
checked, and unlocked products are coarsened only as much as required by the
event cap. If locked products are limiting, an unlocked parent may be split.
The final audit revisits every reachable reaction. The transaction computes
the exact number of required split clones and compares it with the per-rank
buffer ceiling after reserving `RadialAllocationReserveFraction`; any rank's
failure rejects the collective plan before a weight or list is changed.

The old two-phase path (independent block decision followed by component
maximum propagation) has been removed. Commit now applies the already planned
split/roulette or empty-block weight change, synchronizes owners to boundary
copies once, verifies exact locked-source weight bits, and rechecks chemistry.
Transport resumes only after the common component carrier is visible on all
copies.

Pure planner tests are in
`srcHartley2/tests/test_stage2_component_carrier_planner.cpp`; their runner
also exercises synthetic DIM2 and DIM3 graphs without an AMPS mesh. The
integrated validator uses the canonical Halley mesh and one frozen executable
at 1, 4, and 12 MPI ranks. It checks sparse and two-ended chains, empty blocks,
exact split conservation and statistically bounded roulette conservation for
number/mass/momentum/kinetic/internal energy, bitwise source locks, the hard
reaction limit, owner/ghost synchronization, and a forced zero-capacity
rollback. The strict parser reports:

```text
PASS Stage-2 integrated component carriers 1rank-owners=1-ghosts=0 4rank-owners=4-ghosts=28 12rank-owners=12-ghosts=84
```

All runs report maximum face multiplicity 1, maximum reaction allocation 20
for a limit of 20, and `partial-mutations=0`. Evidence is archived under
`validation_runs/stage2_component_carrier_20261105`; the frozen executable
SHA-256 is
`c40c6bb44871209672cd6c7ecc6e16cc838ca5ea3e3e0184970696d273909939`.

This closes the Stage-2 implementation and transaction test matrix. The
science-pilot criteria above still require a new burn-in/production study of
C2, C2H, CN, C2H5, and HCO before the component method is claimed to improve a
converged Rubin comparison. Stage 3 must not use the transaction tests as a
substitute for that science evidence.

## Stage 3: correction-free variance reduction and decision gate

### Objective

Reduce residual-rounding noise without individual weights, then decide from
measurements whether the higher-risk individual-weight stages are necessary.

### Randomized residual accumulators

The current `floor(mu)+Bernoulli(frac(mu))` draw is unbiased independently for
every event but has binomial accumulated variance. Add an optional
random-phase systematic accumulator for each owner block, reaction, and
concrete product:

```text
accumulator += mu
multiplicity = floor(accumulator)
accumulator -= multiplicity
```

Initialize each accumulator with a uniform random phase in `[0,1)`. Over many
events the allocated total then differs from the exact cumulative expectation
by less than one particle per accumulator, while the random phase retains an
unbiased ensemble.

Accumulator state must:

- live in owner block data and migrate with dynamic load balancing;
- be checkpointed if restart support is used;
- never be read or updated by ghosts;
- preserve the all-products transaction and capacity preflight; and
- be flushed or transformed explicitly if a carrier changes.

For products sharing the same carrier, use a common event-packet phase where
possible so stoichiometric partners are created together. When carriers
differ, document which quantities are conserved per packet and which are
conserved only in cumulative expectation.

### Statistical output improvement

- Add optional adaptive logarithmic output bins that merge radial cells until
  a requested particle-time exposure is reached.
- Preserve the original cell output; the merged profile is an additional
  diagnostic and must record its radial support and exposure.
- Add subwindow or independent-seed uncertainty estimates rather than hiding
  zero and negative samples.

### Tests and decision gate

- Compare Bernoulli and accumulator ensembles against analytic mean and
  variance for `mu<1`, integer `mu`, and `mu>1`.
- Validate restart and MPI block migration with nonzero accumulators.
- Run at least four independent pilot seeds after a decorrelation burn-in.
- Compare zero-cell fraction, coefficient of variation, reaction-source
  closure, and window-to-window profile change.

Proceed to Stage 4 only if at least one important minor species remains below
its input target or has unacceptable measured uncertainty after Stages 2 and
3. If C2 and the other requested daughters are resolved, retain the simpler
correction-free model and skip individual weights.

### Stage 3 implementation status (2026-11-05)

The correction-free implementation and its focused test matrix are complete
on `1d_comet_dev`. `ChemistryResidualSampling=SYSTEMATIC` selects a
random-phase systematic accumulator; omission retains the legacy independent
Bernoulli path, allocates no accumulator block storage, and makes no new RNG
call. Systematic mode requires an explicit `CometRandomSeed`. Physical-event
RNG streams use `seed+rank`, while phases use a separate counter-based hash of
the seed, stable AMR node ID, carrier epoch, reaction, and concrete product.

State is registered with explicit double alignment in the application-owned
block buffer after core local weight/timestep storage and before block
allocation. Each owner block retains one phase for every concrete
stoichiometric product. A 53-bit layout signature covers the Stage-3 state
ABI and reaction/product order so a same-size but semantically different
layout is rejected. Exact-equal product carriers in one reaction receive a
common initial phase, but retain distinct
slots so later carrier changes remain unambiguous. `ProcessParticle()` plans
every multiplicity without changing state, performs the Stage-2 hard-cap and
particle-capacity checks, then commits phases only with the complete reaction.
An all-zero daughter event commits before deleting its parent; a capacity
deferral leaves the parent and phases bitwise unchanged. First-use phase
construction is provisional: an unaccepted event never marks or fills the
block state.

An accepted regular carrier transaction closes the current bounded-error
epoch and rephases lazily under the new carrier. Emergency common-factor
coarsening retains phases because it multiplies every block/species weight by
the same factor and leaves every chemistry source/product `W/dt` ratio
unchanged. Dynamic load balancing carries the state with its owner block.
Core local-weight restart now saves every owner block, including empty blocks,
and restores the full application-owned region. After owner records are read,
the core block-data exchange republishes weights, timesteps, and application
bytes to ghosts before the mover can use them. Populated restart records and
particle checksums are validated. The same core change replaces three
hard-coded eight-child restart traversals with `2^DIM`, which is needed for
valid DIM1/DIM2 recursion.

The Bernoulli analytic variance is never presented as a systematic variance.
Systematic diagnostic JSON reports `represented_number_sigma_s-1: null` and a
conservative fixed-carrier-epoch discrepancy bound. CSV leaves the Bernoulli
sigma cell blank and appends the bound plus its uncertainty method. Original
AMPS plot output is unchanged; complete diagnostic windows add authoritative
center-cell `stage3_cells.*.csv` tables because Tecplot rows are
node-interpolated.
`visualization/adaptive_radial_bins.py` adds inner-to-outer exposure-controlled
bins, retains zero and final-underexposed bins, and, when
given independent-seed replicas on a common mesh, reports sample standard
deviation, standard error, coefficient of variation, and nonpositive-seed
fraction on common radial supports. Embedded species/output identity, exact
geometry, sampling length, seed validity, residual mode, and MPI count are
checked; duplicate paths/seeds and mismatched provenance are rejected.
Boundaries remain exact cell faces on the logarithmically expanding radial
mesh; the tool does not interpolate data onto invented logarithmic edges.

Pure analytic tests cover `mu<1`, integer `mu`, `mu>1`, strict prefix error
below one, plan/retry, packet coupling, serialization, and the floating-point
`LONG_MAX` boundary. Eight isolated parser cases cover explicit/default and
invalid/duplicate inputs. Integrated tests
with one frozen input/seed at 1, 4, and 12 MPI ranks cover 373 owner blocks,
36,554 accumulators, common phases, non-mutating plan retries, first-use
rollback, carrier epochs, owner-only access, AMR migration, populated
local-weight restart continuation, particle checksums, and restored ghost
publication. Nine nonzero
states change owner in the four-rank test and 33 in the twelve-rank test; all
373 restart comparisons pass at every decomposition. They also require the
authoritative center-cell writer to complete; each run produces 2,984 C2 cell
records that pass the adaptive-bin command-line smoke test.

The Stage-3 validator does not mislabel its pure plan/retry check as a forced
`ProcessParticle()` capacity failure. The already passing Stage-2 integrated
transaction suite supplies the real capacity-rejection, particle-list, lock,
and no-mutation rollback test required by Stage 3.

The scientific decision policy is now approved and encoded in
`validation/stage3_science_policy.json`, with non-relaxable values checked by
`validation/stage3_gate_common.py`. One frozen executable is used for four
fixed Bernoulli seeds and the same four Systematic seeds at 12 MPI ranks and
one OpenMP thread. The first 10,000-iteration window is burn-in/reference;
the 20,000, 30,000, and 40,000 iteration outputs are the three accepted
10,000-iteration windows. Each must satisfy the existing 5% profile-change,
chemistry/source-closure, and zero-emergency gates.

The statistical comparison uses only exact cells fully contained from 10 to
100,000 km altitude above the 6 km surface. Common adaptive-bin support is
defined by the minimum exposure over all eight runs, not by a favorable
single mode or seed, and requires 100,000 particle-iterations. For C2, C2H,
CN, C2H5, and HCO, at least 90% of Systematic bins must have CV at most 25%,
none may exceed 50%, no bin may be empty, the raw zero-cell fraction may not
exceed 5%, and at least 90% of populated complete blocks must remain inside
the configured 50--200 `CHEMICAL_ONLY` target. Bernoulli and Systematic means
must agree within the larger of 5% of their scale and two combined standard
errors. Reaction/product production agreement, all-species source agreement,
and Rubin log-density regression are independent gates, not substitutes for
minor-species precision.

`validation/run_stage3_science_gate.py` freezes the executable, canonical
input, generated `picGlobal.dfn`, and policy; alternates the eight cases;
records atomic per-attempt status; validates and hashes every required
artifact; and resumes only by starting an incomplete case again at iteration
zero. `validation/evaluate_stage3_science_gate.py` is the read-only scientific
evaluator. It reports `PASS_SKIP_STAGE4`, `FAIL_CONSIDER_STAGE4`,
`INCONCLUSIVE_RERUN`, or `INVALID_EVIDENCE` and writes compact JSON, CSV, and
Markdown. Its first valid use writes an unaccepted Rubin candidate baseline;
human scientific review must explicitly accept that baseline before a final
pass is possible. Full commands, recovery rules, thresholds, and exit codes
are in `validation/README.md`.

The clean current build and 1/4/12-rank focused tests pass, as do 14 synthetic
science-gate boundary/failure cases. The long eight-run science matrix has not
yet been launched, so Stage 3 implementation is complete but its scientific
go/no-go outcome remains pending. Do not begin Stage 4 merely from the focused
tests.

The generic network executor is currently valid only with the radial
batch/capacity/emergency controller. Input preprocessing and runtime
validation reject chemistry in 2-D/3-D rather than inventing policy values.
Generic geometry and source support are not removed; a future geometry-neutral
batch policy is required before this chemistry implementation is enabled and
validated outside the radial backend.

## Stage 4: individual-weight compatibility foundation

### Objective

Make individual statistical corrections compatible with the radial
local-timestep mover, adaptive base weights, general NTC collisions, chemistry,
MPI exchange, and sampling before using them for science.

### Implemented support and remaining exclusions

- The particle buffer stores, clones, exchanges, migrates, and restarts an
  individual correction. Macroscopic sampling uses the physical carrier
  `W*c`.
- The local-timestep face mover now implements the fixed-correction residual
  policy below. It validates positive finite `c`, preserves `c` in every
  survivor/clone, and measures number, mass, momentum, kinetic energy, and
  internal-energy flux with `W*c/dt`.
- Adaptive and emergency block-weight changes preserve `c`; their conservation
  diagnostics use corrected carriers. Capacity failure remains transactional.
- General NTC now selects a particle in proportion to its physical carrier,
  rejects same-particle pairs without conditioning a redraw, and validates
  every unequal-carrier update probability. A disabled-by-default diagnostic
  exposes focused-test counts and extrema without changing production output.
- Comet chemistry creates exactly one corrected daughter for each concrete
  stoichiometric product in correction-on builds. This made the existing
  carrier-conversion formula reachable at the Stage-4 boundary; the
  separately tested Stage-5 windows are documented below. Channel/reaction
  importance sampling remains unimplemented.
- The validated collision path is general NTC only. Simplified NTC terminates
  explicitly at entry, and the radial setup rejects MF or another collision
  model when corrections and collisions are both enabled.
- MPI-only radial DIM1 validation is complete. Hybrid diagnostic reduction,
  CUDA continuation, field-line movers, periodic topology, and 2-D/3-D
  correction-on chemistry have not been validated and are not claimed.

### Face-transfer design

For a particle with correction `c`, fixed-correction residual transfer uses

```text
mu = (W_old*c/dt_old)/(W_new*c/dt_new)
   = (W_old/dt_old)/(W_new/dt_new).
```

The current base multiplicity remains mathematically valid if every survivor
and clone keeps `c`. The mover diagnostics must multiply incoming and outgoing
number, mass, momentum, and energy carriers by `c`.

Implement and validate one explicit policy rather than mixing two policies:

1. **Residual face transfer:** retain `c` and split/roulette by the base
   carrier ratio; or
2. **Correction transformation:** retain one particle and transform `c` by
   the old/new carrier ratio.

The first policy matches the current correction-free mover and is preferred
for initial implementation because it preserves the existing allocation
diagnostics. Do not retain the legacy correction-transform path as an
unreviewed second route.

### Adaptive base-weight interaction

- Halving a base weight and cloning each particle once preserves represented
  physical weight exactly when both copies retain `c`.
- Doubling a base weight and retaining each particle with probability one
  half preserves represented weight in expectation when survivors retain
  `c`.
- Component-wide carrier changes and emergency adaptation must include
  corrections in diagnostic totals but need not modify `c` under this policy.
- Capacity preflight must count all clones independently of correction.

### Tests

- Extend the face-flux validator with corrections spanning the intended
  operational range and with unequal base weights and timesteps.
- Verify number, mass, momentum, kinetic energy, internal energy, and remainder
  fraction through multiple local and remote faces.
- Exercise MPI continuation, ghost exchange, dynamic load balancing, inner
  re-emission, and outer deletion.
- Add general NTC relaxation and collision-rate tests with unequal individual
  corrections, species weights, timesteps, and internal degrees of freedom.
- Assert that all correction-derived collision update probabilities remain in
  `[0,1]` within tolerance.
- Confirm that simplified NTC and unsupported collision paths fail at setup
  rather than silently ignoring corrections.
- Repeat adaptation and emergency tests at 1, 4, and 12 ranks.

### Exit criteria

- The compile-time prohibition is removed only for the validated combination.
- Every new feature remains disabled by default.
- All conservation and MPI tests pass with corrections both equal to one and
  widely unequal.
- A correction-on run with importance disabled reproduces the
  correction-off reference within independent Monte Carlo uncertainty.

### Completion record (2026-11-05)

Stage 4 is complete on `1d_comet_dev` for the supported radial configuration.
The user explicitly authorized this implementation after the Stage-3 focused
gate was completed; this does not retroactively claim that the pending long
Stage-3 science matrix selected Stage 4.

The validated compile-time/runtime combination is:

```text
DIM1 radial mesh
SpeciesLocalTimeStep
SpeciesLocalParticleWeight
ParticleWeightCorrectionMode=on
local-timestep face-flux resampling=on
general NTC collisions
MPI-only focused diagnostics
internal degrees of freedom enabled by the configured species data
```

The correction-on archive is
`validation_runs/stage4-correction-on-full-20261105`. It exercises 200,000
mover particles per case and 200,000 chemistry events at 1, 4, and 12 MPI
ranks. Each decomposition covers five exact correction values from `1/16` to
`16`, all-correction-one identity runs, local and remote faces, multiple-face
continuation, adaptation, emergency adaptation, sampling, general NTC
relaxation, boundary deletion/re-emission, restart, dynamic migration, ghost
publication, and chemistry. The archive also proves that simplified NTC fails
explicitly. A separately configured MF setup was rejected before mesh setup in
`validation_runs/stage4-unsupported-mf-20261105b`.

The correction-off identity reference is
`validation_runs/stage4-correction-off-reference-20261105`. The strict checker
requires identical path-coverage counts at every MPI decomposition and
requires both correction-on/all-`c=1` and correction-off arms independently to
pass the same number, mass, momentum, kinetic/internal-energy, remainder, and
post-move flux tests. The maximum physical-flux checksum error was
`5.48e-12` with corrections on and `1.39e-11` with corrections off, both below
the frozen `1e-9` limit. Run the combined audit as documented in
`srcHartley2/README.md`.

Correction-off Stage-2 adaptation and Stage-3 systematic-residual regressions
also pass at 1, 4, and 12 ranks in
`validation_runs/stage4-correction-off-regression-20261105`. Stage-0, Stage-1,
Stage-2 planner, Stage-3 accumulator/parser/science-gate synthetic tests, and
adaptive-bin/diagnostic-semantics tests pass. These are compatibility tests;
they are not a replacement for a long converged comet-science comparison.

The repository's canonical Halley input remains correction-off with
Systematic chemistry residuals. Individual corrections and all Stage-4/5
diagnostics remain opt-in. Stage 5 is now implemented and focused-tested, but
its numerical window values are a validation policy, not a calibrated Rubin
science policy. The normal Halley input therefore intentionally contains no
correction-window records.

The complete generic-input, I3/brightness, and Stage-0--10 regression stack is
orchestrated by `tests/run_all_stage_tests.sh`. Its default archived mode reads
the checksummed compact evidence under `tests/data`, so a clean clone does not
depend on machine-local validation archives or a configured AMPS build. It
rejects missing data instead of accepting child-runner skip behavior. `--plan`
performs a read-only preflight; optional fresh modes use separately supplied
executables and never reconfigure the canonical workspace. Compact protocol
fixtures validate implementation contracts and do not replace complete
scientific archives or a new converged model run. See `tests/README.md` for the
operator contract and fixture boundary.

## Stage 5: exact individually weighted daughters and weight windows (complete)

### Objective

Create each stoichiometric product for every selected analog reaction event,
eliminating zero-daughter roulette caused solely by a coarse product carrier.

### Chemistry algorithm

Keep total reaction occurrence and channel selection analog. For a corrected
parent with physical carrier

```text
C_parent = W_parent*c_parent/dt_parent,
```

create one particle for each concrete stoichiometric product and assign

```text
c_daughter = C_parent/(W_daughter/dt_daughter).
```

This formula was made reachable by Stage 4. Stage 5 leaves that chemistry
transaction unchanged and adds a post-step correction controller around it.

With one individual particle per concrete product, the hard transaction count
is the stoichiometric product count rather than
`ceil(C_parent/C_product)`. Keep a separate input hard limit and preflight the
complete transaction before modifying the parent.

### Input policy

Stage 5 adds these records to an input that also enables
`ParticleWeightCorrectionMode`:

```text
ParticleWeightCorrectionWindowMode=on
ParticleWeightCorrectionWindow=LOWER | TARGET | UPPER
ParticleWeightCorrectionWindow(DIRECT_SOURCE)=LOWER | TARGET | UPPER
ParticleWeightCorrectionWindow(CHEMICAL_ONLY)=LOWER | TARGET | UPPER
ParticleWeightCorrectionWindow(SPECIES)=LOWER | TARGET | UPPER
ParticleWeightCorrectionWindowMaxSplits=N
```

The unqualified record is mandatory when window mode is on. Resolution is
deterministic: an explicit species overrides its derived source/chemical
group, which overrides the global record. All values must be positive finite
numbers with `LOWER < TARGET < UPPER`; a directly injected species must have a
window containing `c=1`. The parser also rejects a window inconsistent with
the fixed split rule and rejects every record/limit when the mode is off or
omitted. Omission therefore preserves legacy behavior exactly.

The focused acceptance build used global `0.25|1|4`, direct-source
`0.5|1|2`, chemical-only `0.0625|0.25|4`, and a per-particle maximum of 32.
Those values were staged only for the frozen Stage-5 executable. They are not
present in canonical `input/halley1d.input` and must not be interpreted as a
scientifically selected importance policy.

### Weight-window transaction

Individual corrections must not drift over an unbounded range. Add a
species/group input policy defining lower, target, and upper correction or
physical-carrier limits.

The inclusive band is inert: `LOWER <= c <= UPPER` leaves the particle
unchanged. For `c<LOWER`, roulette to `TARGET` with survival probability
`c/TARGET`. For `c>UPPER`, use `n=ceil(c/TARGET)`, limited by the input maximum,
and create `n` particles. The first `n-1` receive `c/n`; the final particle
receives the floating-point remainder, and the production-order sum must
equal the original `c` exactly. Thus the committed split preserves the
carrier without silently discarding roundoff.

The controller runs collectively after every completed `PIC::TimeStep()`,
after a coincident base-weight adaptation and before the next sample. It
traverses authoritative owner blocks only and rejects any nonzero mover
continuation remainder. Preflight performs no random draw or mutation, counts
the exact `sum(n-1)` clone demand without crediting possible roulette
deletions, retains the configured allocation reserve, and obtains an all-rank
accept decision. One rank's capacity/configuration failure rolls back the
entire batch before RNG, list, tracker, or diagnostic state changes.

The commit clones the complete particle record, so position, velocity,
internal modes, species, chemistry stage, and application state are retained;
only `c` and a clone's tracker identity change. Roulette deletion repairs the
cell list through the core particle-buffer API. Correlated split copies do not
qualify as independent samples until a measured transport/collision burn-in.

The pure implementation is in `CorrectionWindowPlanner.h`. It has no AMPS,
MPI, mesh, or particle-buffer dependency. `RadialSimulation.cpp` provides the
DIM1 owner-list executor. This separation is intentional preparation for a
future geometry-neutral DIM2/DIM3 executor; Stage 5 does not claim that
support yet.

### Tests and evidence

- `tests/run_stage5_correction_window_tests.sh` checks parser precedence,
  malformed/duplicate/off-mode cases, direct-source identity, keep bounds,
  exact split sums, roulette expectation, hard split/capacity failures,
  multiple species, repeated products, and empty input.
- `--validate-correction-windows` uses an initialized mesh and production
  executor to check inclusive bounds, full-state splitting, collective
  zero-capacity rollback, correction/restart moments, 65,536 roulette
  candidates per owner block, direct-source `c=1`, and an empty list.
- `tests/run_stage5_integrated_tests.sh` runs that validator and the exact
  chemistry ensemble at 1, 4, and 12 MPI ranks. The strict Python parser
  requires every named path and identical resolved policy/chemistry coverage.
- The completed 200,000-event archive is
  `validation_runs/stage5-correction-window-full-20261105`; it reports
  `PASS Stage-5 correction-window integrated matrix 1rank=pass 4rank=pass
  12rank=pass`.

### Completion boundary and remaining scientific work

Stage 5 is closed as an implementation/compatibility stage: exact concrete
products, bounded windows, transactional capacity behavior, restart, and MPI
decomposition invariants have passed. It is not a long-run C2 science pass.
Selecting species-specific windows, measuring clone decorrelation burn-in,
checking parent profiles and C2 flux closure over independent seeds, and
comparing a converged run with Rubin remain Stage 9 work. Packet splitting and
channel importance sampling remain Stages 6 and 7; the Stage-5 controller must
not be described as either feature.

## Stage 6: daughter packet splitting

Stage 6 is implemented and its focused compatibility gates are closed for the
same correction-on DIM1/general-NTC configuration qualified by Stages 4 and
5. It creates several independent recoil realizations for one selected
chemical event without changing the represented event rate.

### Input policy and precedence

The optional records are:

```text
ChemistryDaughterPacketMode=on
ChemistryDaughterPacketMaximum=8
ChemistryDaughterPacketCount(GLOBAL)=1
ChemistryDaughterPacketCount(CHEMICAL_ONLY)=2
ChemistryDaughterPacketCount(SPECIES:C2)=8
ChemistryDaughterPacketCount(REACTION:H2O_OH_H)=4
ChemistryDaughterPacketCount(REACTION:C2_2C)=8
```

The mode and maximum are required when packet splitting is enabled. A
reaction-ID record has highest precedence. Otherwise, the resolver evaluates
every concrete product and selects the maximum applicable product count; for
each product, explicit species overrides `CHEMICAL_ONLY`, which overrides the
global count. This rule handles repeated products and mixed-product reactions
without relying on generated reaction-array positions.

Product-free sinks always remain `N=1`. Explicitly targeting a sink is an
input error even when the requested count is one, because such a selector
cannot produce a packet. Enabled records that only match zero-rate channels,
duplicate selectors, unknown species/reactions, nonpositive/nonintegral
counts, counts above the declared maximum, inactive-mode packet records, and
any `K*N` transaction above
`RadialMaxPhotochemicalProductsPerReaction` stop preprocessing. Omitting the
mode and all packet records is the strict legacy default: generated counts
are one and production uses the historical one-recoil path.

### Production transaction

For a selected reaction with `K` concrete stoichiometric products and packet
count `N`, chemistry reserves and creates exactly `K*N` daughters. Every
packet calls the existing recoil generator independently and contains all `K`
products from that one recoil solution. Momentum and excess energy are
therefore conserved inside each packet; merely cloning one velocity is never
used as a variance-reduction substitute.

For product `k`, the Stage-4 unsplit correction is

```text
c_k = (W_parent*c_parent/dt_parent)/(W_k/dt_k).
```

`DaughterPacketPlanner.h` partitions it into `N` positive corrections. The
first `N-1` receive `c_k/N`, and the last receives the arithmetic remainder in
production order, so the stored corrections reconstruct `c_k` exactly. The
planner checks integer overflow, declared packet maximum, `K*N` hard limit,
and available particle capacity before mutation. Production reuses the parent
slot once, allocates the remaining `K*N-1` records, preserves position and
application state, assigns product species/velocity/correction, initializes
new tracker identities, and marks every committed daughter as already
processed for the current chemistry stage.

The batch capacity mirrors in `RadialSimulation.cpp` use the same `K*N`
cardinality for normal chemistry, emergency forecasts, and block-weight
adaptation. A failed preflight leaves the parent record byte-identical and
defers the event before RNG-dependent commit. Stage-5 correction windows run
normally after the fixed packet transaction; packet daughters receive no
special exemption and can subsequently be kept, rouletted, or split. This is
why scientific packet policies must be evaluated together with their window
policy.

### Focused validation

The pure suite is:

```sh
srcHartley2/tests/run_stage6_daughter_packet_tests.sh
```

It covers default/off identity, all precedence levels, maximum-across-product
resolution, repeated products, sink and zero-rate rejection, malformed and
duplicate records, generated-header integrity, `N=1,2,4,8`, exact correction
sums, overflow, hard-limit, and capacity rollback.

An initialized Stage-6 build exposes:

```text
--validate-daughter-packets --chemistry-validation-events N
```

For every one of the 46 generated neutral reactions, the validator compares
`N=1,2,4,8` using deterministic global event IDs. It checks all concrete
products, exact correction reconstruction, packet momentum/energy,
unchanged weighted velocity/density means, and the expected `1/N` estimator
variance. Centered online moments are accumulated in extended precision and
merged with the parallel `(count,mean,M2)` identity, preventing MPI grouping
roundoff from masquerading as a different event ensemble. A separate real
`C2_2C -> C+C` transaction checks 16 committed daughters for `N=8`, parent
slot reuse, list links, processing flags, represented moments, independent
recoils, and byte-exact capacity rollback.

The accepted full matrix is:

```text
validation_runs/stage6-daughter-packet-full-20261105-r4
PASS Stage-6 daughter-packet integrated matrix 1rank=pass 4rank=pass 12rank=pass reactions=46 sinks=23 packet-cases=4 events=200000
```

The validation policy used global 1, `CHEMICAL_ONLY` 2, species `C2` 8,
reaction `H2O_OH_H` 4, and reaction `C2_2C` 8, with packet maximum 8 and a
64-product transaction limit. Broad Stage-5 windows were used only so the
focused `N=8` trajectories survived to the validator. These are test values,
not calibrated Halley settings.

### Completion boundary

Stage 6 is closed as an implementation and focused-compatibility stage. It
does not yet select a production packet count from independent-seed
variance/cost, prove primary-profile invariance in a transported converged
run, establish minor-species flux closure, measure buffer headroom, or update
the Rubin comparison. Those scientific gates remain Stage 9. DIM2/DIM3
corrected chemistry also remains unsupported even though the resolver and
pure planner are geometry-neutral.

## Stage 7: reaction-channel importance sampling

Stage 7 implements stochastic conditional-channel biasing. It does not force
the parent to react: the parent total rate, exponential event probability,
reaction position, and destruction count remain analog. Only the mutually
exclusive channel selected after an ordinary total-reaction event is changed.

### Accepted input contract

The exact input grammar is:

```text
ChemistryChannelImportanceMode=ON
ChemistryChannelBiasedProbability(REACTION:C2_2C)=0.5
ChemistryChannelBiasedProbability(REACTION:C2_ION)=0.5
```

`REACTION:` is syntax, while the following reaction ID is stable and
case-sensitive. Generated reaction-array indices, species selectors, group
selectors, global defaults, and implicit remainders are intentionally not an
input interface. A parent family remains analog unless at least one of its
channels is named. Once one positive-rate unit-unimolecular channel is named,
every positive-rate sibling for the same parent must be named, including
product-free neutral sinks/ionization losses. Each requested `q*` must be
finite and in `(0,1]`, and the complete family must sum to one within the
strict parser tolerance. The runtime assigns the last positive interval the
arithmetic remainder, so the actual selection measure closes at exactly one.

The resolver rejects:

- records paired with omitted or explicit `OFF` mode;
- `ON` mode without any probability record;
- duplicate, empty, unknown, or case-mismatched reaction IDs;
- zero-rate and non-unit-unimolecular targets;
- partial positive-rate sibling families;
- zero, negative, non-finite, greater-than-one, or nonnumeric probabilities;
- family sums outside the normalization tolerance; and
- malformed reaction descriptors or generated arrays.

An untouched family has generated `Configured=0` and
`BiasedProbability=0`. With the mode omitted, every family follows the exact
historical raw-rate selector. A configured one-channel family with `q*=1` is
valid as an identity case but provides no variance reduction. A zero physical
rate can never receive positive statistical support.

Stage 7 requires comet chemistry, DIM1, individual particle corrections,
Bernoulli residual mode, Stage-5 correction windows, and comet diagnostics.
The mode is independently default-off. It does not require Stage-6 packet
splitting, but it composes with it when both switches are enabled.

### Runtime estimator and conservation contract

For parent `s`, the existing cache retains each physical rate `k_j` and

```text
k_total = sum_j k_j,
q_j     = k_j/k_total.
```

The total reaction test remains exactly
`1-exp(-k_total*dt_parent)`. Unconfigured families then retain the historical
`rnd()*k_total` and left-to-right raw-rate subtraction, including its random
draw count and floating-point channel boundaries. Configured families use one
uniform draw against the complete positive `q_j*` table. For the selected
channel,

```text
L_j = q_j/q_j*.
```

The parent is destroyed exactly once with its original analog carrier. This
destruction must not be multiplied by `L_j`, because every channel in the
family destroys the same already-reacted parent. A selected product-free sink
therefore creates no daughter but records its raw selection and likelihood so
its weighted channel yield remains auditable.

For product species `d`, Stage 4's unsplit corrected carrier becomes

```text
c_daughter =
    [W_parent*c_parent/dt_parent] /
    [W_daughter/dt_daughter] * L_j.
```

If Stage 6 requests `N` independent complete daughter packets, its exact
final-remainder planner partitions that value:

```text
c_daughter,i =
    [W_parent*c_parent/dt_parent] /
    [W_daughter/dt_daughter] * L_j/N,
```

with the last packet carrying the arithmetic remainder. Thus Stage 7 changes
no `K*N` allocation cardinality, recoil sampling, parent-slot reuse, tracker
state, or transaction-capacity bound. Capacity rejection remains mutation- and
commit-diagnostic-free. After the chemistry transaction, Stage-5 windows run
normally on the likelihood-weighted corrections; they may keep, roulette, or
split them without a Stage-7 exemption.

For any channel observable `f_j`, the estimator is unbiased because

```text
sum_j q_j* L_j f_j = sum_j q_j f_j.
```

The application repeats all generated-family, positive-runtime-rate,
normalization, and likelihood checks when it builds the rate cache. This
catches stale or hand-edited generated headers before particle processing.

### Diagnostics and effective sample size

When Stage 7 is enabled, each committed selection records `L`, `L^2`, and the
likelihood-weighted channel rate. Separate files preserve the legacy Stage-0
schema and off-mode output identity:

```text
PT/diagnostics/comet_stage7_importance.out=<index>.csv
PT/diagnostics/comet_stage7_importance.out=<index>.json
```

The Stage-7 schema reports stable reaction ID, parent, configured/sink flags,
physical `q`, actual biased `q*`, `L`, raw selections, likelihood sums,
weighted channel rate, and family totals. The sampled family effective sample
size is

```text
ESS = (sum L)^2/sum(L^2).
```

This is selection-weight ESS for one parent-channel family. It is not a claim
of independent transported particle samples after collisions, packet cloning,
or Stage-5 resampling. Production diagnostics are mandatory when importance
mode is on so sink yields and ESS cannot disappear silently.

### Focused tests and evidence status

The pure runner is:

```sh
srcHartley2/tests/run_stage7_channel_importance_tests.sh
```

Its policy suite covers default/off identity, stable case-sensitive IDs,
complete two-/many-/single-channel families, sinks, tiny/dominant `q*`, input
order, malformed/duplicate/unknown records, partial families, zero-rate and
non-unit targets, normalization, generated-header sentinels, and defensive
reaction descriptors. `ChannelImportancePlanner.h` covers exact analog
rate-subtraction identity, selection boundaries, likelihoods, weighted means,
variance, disabled channels, and analytic/sample ESS.

The initialized runtime entry point is:

```text
--validate-channel-importance --chemistry-validation-events N
```

The runtime test must compare deterministic global-event analog and biased
ensembles for every configured parent family, including sinks. It must require
unchanged total parent events; physical weighted channel yields; number, mass,
momentum, and energy agreement; the predicted variance change; valid ESS; and
identical global statistics at 1, 4, and 12 MPI ranks. The Stage-6 real
`K*N` transaction validator must also be replayed with Stage 7 enabled so
stored daughter corrections prove the combined `L/N` rule. Stage-5 window and
Stage-4 exact-chemistry dependencies must continue to pass with the same
frozen executable.

The completed focused evidence is:

```text
Comet channel-importance policy test passed
PASS Stage-7 channel-importance planner analog-identity=1 two-channel=1 many-channel=1 dominant-tiny=1 disabled=1 boundaries=1 weighted-mean=1 variance=1 ess=1 errors=1
PASS Stage-5 correction-window integrated matrix 1rank=pass 4rank=pass 12rank=pass
PASS Stage-6 daughter-packet integrated matrix 1rank=pass 4rank=pass 12rank=pass reactions=46 sinks=23 packet-cases=4 events=200000
PASS Stage-7 channel-importance integrated matrix 1rank=pass 4rank=pass 12rank=pass families=27 reactions=69 sidecar-rows=69 events=200000
accepted archive: validation_runs/stage7-channel-importance-full-20261105-r2
PASS Hartley2 all-stage suite tests=25 test-groups=8 stage4-mode=archived stage5-mode=archived stage6-mode=archived stage7-mode=archived
```

The archive freezes the Stage-7-on executable, validation input, generated
core configuration, generated importance policy, and generated packet policy
under `provenance/SHA256SUMS`. The first unsuffixed archive is development
evidence: its AMPS cases passed, but the original parser tolerance was too
tight for ordinary binary64 reduction-order roundoff. The accepted `-r2`
archive completed the end-to-end runner after limiting an `8*N*epsilon` bound
to accumulated diagnostic fields; integer counts, metadata, and `q/q*/L`
remain strict. A focused validator closes implementation compatibility only;
transported independent-seed science validation remains Stage 9.

### Scope limitation and later master merge

Deterministic all-channel branching is deferred. It is not an alternate
Stage-7 mode in the current implementation. Such a mode would have to retain
or explicitly represent the sink/loss branch, allocate several channel
transactions for one analog event, introduce a separate hard capacity bound,
and pass its own parent-survival and conservation proof. Do not infer this
behavior from stochastic `q*` selection.

For a later authorized master integration, preserve this order after Stages
0--6 have passed on a disposable clean integration branch:

1. add the defensive reaction-family descriptor, input policy module, and
   pure Perl tests;
2. add the pure planner and C++ unit test, proving analog identity before any
   application runtime change;
3. add the self-contained generated-header fallback and `hartley2.pl` parser,
   then inspect one omitted/off input and one complete C2 family;
4. add rate-cache policy validation and conditional selection while retaining
   the byte-for-byte analog selection branch;
5. apply `L` once to the Stage-4 unsplit correction before Stage-6 divides by
   `N`; preserve Stage-5 scheduling and all transaction rollback invariants;
6. add separate diagnostics/ESS output, CLI validation, strict log parser, and
   1/4/12-rank runner;
7. freeze and hash a clean Stage-7-on executable/input/generated policy, run
   all prerequisite matrices into a new immutable archive, then restore and
   clean-build the canonical Stage-7-off input;
8. rerun the complete Stage 0--7 suite and review off-mode RNG/output identity,
   parent lifetime, sink weighting, channel yields, conservation, ESS, and MPI
   decomposition invariance; and
9. present the exact diff, exclusions, hashes, validation records, and Stage-9
   science limitations to the user, then stop before any `master` operation.

Generated binaries, build copies, raw output, validation archives, and the
pre-existing unrelated `Makefile` change must not enter the source commit.
Documentation of a merge sequence is not permission to merge.

## Stage 8: optional total-reaction forcing

### Objective

Address a species that remains unresolved because the total parent reaction
probability—not merely the channel or product carrier—is too small.

This stage is optional. It was implemented after Stages 2 through 7 exposed
total reaction occurrence, rather than only conditional channel selection or
daughter carrier size, as a remaining source of variance.

### Implemented deterministic two-branch estimator

There is deliberately no biased occurrence probability `P*` in the input.
For a configured parent with total constant rate `lambda` during one logical
local-timestep segment `dt`, the planner evaluates

```text
tau = lambda*dt
P = -expm1(-tau)
S = 1-P.
```

An incoming physical carrier `C=(W_parent/dt_parent)*c_parent` is represented
by both outcomes: an unchanged parent carrier `C*S` and a reacted carrier
`C*P`. The reacted value is adjusted by at most a few binary64 ulps and the
survivor is its arithmetic remainder, so the committed pair adds to `C`
exactly in production precision. This is deterministic implicit branching,
not a stochastic `P*` estimator. The formerly proposed `P/P*` reacted and
`(1-P)/(1-P*)` survivor likelihoods remain a rejected alternative; omitting
the latter would bias the parent lifetime.

The reacted branch retains the parent species until its physical reaction
point. With `u` uniform on `[0,1)`, the conditional event time is

```text
t = -log1p(-u*P)/lambda,       0 <= t <= dt.
```

Stage 7 then selects a conditional channel. Its likelihood `L=q/q*` is
applied once; Stage 6 packet count `N` divides each concrete product carrier.
For every concrete product, the final represented carrier is therefore
`C*P*L/N` per packet. A sink has one pending marker and no daughter. A
productive reaction with `K` concrete products allocates `K*N` pending
markers plus the survivor. All markers initially remain parent particles so
ordinary mover geometry, face crossings, MPI exchange, and local-timestep
continuation determine the event location before conversion.

Each marker carries a compact 16-byte application record containing the
reaction/product code, deterministic recoil key/leader flag, and remaining
normalized parent interval. The first-order force-free traced mover queries
the event without side effects, advances the record by the elapsed segment,
and calls the application at the event point. Conversion changes species,
restores all products in a packet from the same deterministic conservative
recoil key, scales the correction by the parent/product `W/dt` ratio, and
rescales the unspent interval by `dt_product/dt_parent`. Only the packet
leader performs once-per-reaction destruction and selection accounting.

The transaction preflights the maximum `K*N` family allocation before the
first Stage-8 RNG draw or mutation. Failure leaves the parent record, path,
RNG state, and particle buffer unchanged. Stage-5 windows run at the normal
completed-cycle boundary; pending in-flight records must have completed
before then. The particle record is requested only when the mode is enabled,
and omitted/off inputs leave the three core callbacks null.

### Input and supported scope

The exact input spelling is:

```text
ChemistryTotalReactionForcingMode=on
ChemistryTotalReactionForcing(SPECIES:C2)=on
```

Only explicit, case-sensitive compiled species are accepted. Numeric values,
`P*`, duplicate/unknown species, zero-rate or non-unit-unimolecular parents,
records in off mode, and an enabled mode without a species are rejected.
The enabled mode requires DIM1 spherical radial chemistry, individual
corrections, Bernoulli residual sampling, Stage-5 correction windows, an
explicit photochemical product cap, and diagnostics. The generated table is
one immutable flag per species.

The production hook is currently qualified only for the legacy first-order,
force-free traced mover used by the DIM1 radial application. The generic core
interface is geometry-neutral, but no second-order, field-line, periodic,
CUDA, DIM2, or DIM3 event integration is claimed. The runtime currently
assumes a constant total reaction rate over one logical local-timestep
segment. The pure planner and analytic oracle cover piecewise-constant
hazards, but production does not yet accumulate a multi-block hazard plan.

### Tests

- `run_stage8_total_reaction_forcing_tests.sh` covers the strict Perl policy,
  parser syntax, exact carrier closure, homogeneous and piecewise hazards,
  conditional positions, capacity arithmetic, and invalid inputs. Its pure
  planner also advances an analytic discrete-cycle `A -> B -> C` network,
  proves `A+B+C` carrier closure at every cycle, compares 64 keyed analog
  ensembles with the closed form, and demonstrates that deterministic total-
  occurrence forcing removes occurrence-count variance.
- `--validate-total-reaction-forcing` exercises homogeneous survival,
  constant/piecewise positions, multistep exponential decay, the real
  `C2_2C` Stage-7 likelihood plus Stage-6 `N=8` composition, a production
  transaction, mover-event conversion/continuation, and byte/RNG/path
  rollback on forced capacity failure. The mover fixture invokes
  `ExecutionTrackDefault`, crosses a local face at one rank and an MPI-owned
  face otherwise, places sixteen productive events after the face, and puts a
  product-free sink bitwise at the face time. It requires face precedence,
  advance-before-copy, destination-owner zero-time processing, mover-owned
  deletion, completed product endpoints, and no residual marker.
- `run_stage8_integrated_tests.sh` first replays the Stage-4 local-timestep,
  unit-correction, and weight-adaptation validators and then the accepted
  Stage-5, Stage-6, Stage-7, and Stage-8 validators at exactly 1, 4, and 12
  ranks with 200,000 mover particles and 200,000 chemistry events. Independent
  strict parsers own every lower-stage and Stage-8 schema.
- The complete orchestrator adds five Stage-8 jobs and reports 30 tests in
  nine groups. Its archived mode consumes the immutable accepted archive;
  fresh mode requires a caller-supplied enabled executable and a new output
  root.

### Exit criteria

- Parent survival, reaction-position distribution, and daughter production
  pass analytic tests.
- Variance reduction justifies the additional branching and collision
  complexity.
- The feature remains opt-in and is not required by ordinary comet cases.

### Closure and remaining scientific gate

Stage 8 is closed as an implementation and focused-compatibility stage on
`1d_comet_dev`. The accepted archive is
`validation_runs/stage8-total-reaction-forcing-full-20261105-r4`; its
`provenance/SHA256SUMS` freezes the enabled executable, validation input, core
configuration, and generated Stage-6/7/8 policies. The `-r1` sibling is
retained as failed development evidence: every MPI case ran, but the first
validator reported product markers rather than one logical reacted branch in
its summary counter. The `-r2` retry fixed that counter but exposed an invalid
checker requirement that per-rank transaction roundoff maxima be identical.
Every residue was already below the strict physical tolerance. The `-r3`
archive then passed the direct callback algebra, but audit showed that it did
not run pending markers through the production mover and did not replay the
Stage-4 mover matrix; it is therefore pre-acceptance evidence. The accepted
`-r4` checker retains exact rank identity for deterministic analytic fixtures,
requires each independent transaction residue to pass its bound, and requires
the real mover/MPI/sink record plus the complete Stage-4 dependency replay.

This closure is not a claim of transported science convergence or useful
end-to-end variance reduction. The pure `A -> B -> C` gate establishes the
estimator's occurrence-variance mechanism; independent-seed transported
analog/forced comparisons, Stage-8 product collision feedback, primary-
profile invariance, daughter flux closure, buffer/runtime cost, effective
sample size, and a converged Rubin comparison remain Stage 9.

## Stage 9: long-run science validation

### Objective

Demonstrate that improved statistics produce a converged, physically
interpretable Halley result rather than merely a smoother curve.

### Validation sequence

1. Clean the exact `build/` directory and rebuild from the archived input.
2. Run all parser, reaction-network, chemistry-ensemble, collision, mover,
   adaptation, weight-window, MPI, and whitespace tests required by the
   enabled stage.
3. Run a short one-rank functional case and a multi-rank capacity case.
4. Run at least four independent statistical pilots using identical physical
   inputs and different documented random seeds.
5. Select numerical population and importance controls from measured
   uncertainty and cost. Do not tune them to the Rubin curve.
6. Start a fresh long run from iteration zero. Do not reuse the nonconverged
   final state as a science restart unless a separately validated restart
   procedure transforms all carrier/correction state.
7. Discard the measured clone/weight-window decorrelation burn-in.
8. Require the configured density-convergence criterion for every species and
   independent source/loss/escape flux closure for minor daughters.
9. Accumulate multiple complete steady windows for uncertainty estimates.
10. Regenerate Rubin comparison PNG, EPS, PDF, CSV, and JSON products with
    explicit statistical method, seed, effective exposure, and convergence
    metadata.

### Acceptance metrics

- No particle-buffer exhaustion or unplanned emergency coarsening.
- No unfinished local-timestep trajectory at a chemistry, collision,
  adaptation, sampling, or load-balance boundary.
- Direct-parent source rates and statistically resolved profiles agree with
  the correction-free reference within combined uncertainty.
- Weighted reaction production agrees with reaction counters and flux closure.
- Minor-species zero-cell fractions and effective sample sizes meet their
  input requirements over the comparison range.
- Independent-seed confidence intervals are reported for C2 and every other
  claimed minor-species result.
- The run reaches the configured convergence criterion; a hard-limit result
  remains labeled nonconverged regardless of visual smoothness.

### Stage-9 implementation status (2026-11-05)

The measurement and comparison infrastructure is now implemented on
`1d_comet_dev`; the long independent-seed science decision has **not** yet
been run or accepted.  One Stage-6/7/8-enabled executable can select either
estimator at startup with

```text
--chemistry-estimator-mode analog
--chemistry-estimator-mode enhanced
```

`analog` uses one daughter packet, physical conditional channel probabilities,
unit likelihood, and ordinary stochastic total reaction occurrence.
`enhanced` executes the compiled Stage-6 packet, Stage-7 channel-importance,
and Stage-8 total-occurrence policies.  Both choices retain the same particle
layout, generated tables, callbacks, physical reaction network, collision
database, mesh, and executable.  `configured` remains the default and keeps
the input-selected production behavior.  Runtime selection is rejected after
particle storage is reserved, and the compiled policies are still validated
in the analog arm so that an invalid enhanced configuration cannot hide behind
the comparison mode.

The opt-in input record

```text
Stage9ScienceDiagnostics=ON
```

requires DIM1 radial chemistry, comet diagnostics, individual particle
corrections, Bernoulli residual sampling, and enabled Stage-5/6/7/8 policies.
At each complete sampling window it writes

```text
PT/diagnostics/comet_stage9_checkpoint.out=<N>.csv
PT/diagnostics/comet_stage9_checkpoint.out=<N>.json
PT/diagnostics/stage9_carriers.<SPECIES>.out=<N>.csv
```

For a model particle the carrier used in the effective-sample-size audit is
`q=W_block,species*c_particle`, not the correction alone.  The per-cell Kish
ESS is `(sum q)^2/sum(q^2)`.  Checkpoints also record physical inventory,
correction extrema, occupied cells, the minimum positive cell ESS, sampling-
window wall time, rank capacity, peak buffer utilization, estimator mode,
seed, and MPI size.  Safe-boundary scans at load balance, regular adaptation,
and checkpoint time fail immediately if the linked-list live count disagrees
with `NAllPart`, a local-timestep remainder is pending, or a Stage-8 reaction
event remains armed.

The focused `--validate-transported-daughter-collision` mode extends the
accepted Stage-8 production-mover fixture.  Its actual transported C daughter
records from `C2 -> 2C` enter the production general-NTC kernel with H2O.  The
strict 1/4/12-rank parser requires pointer/owner identity, velocity feedback
in both populations, unchanged species/correction/event/remainder records,
zero unfinished trajectories, valid unit probabilities, represented-rate
closure, and weighted momentum closure.  This is a collision-wiring and
correction-compatibility test; it is not an independent-seed convergence or
Rubin-validation result.

The checksum-protected science gate consumes at least four balanced seeds per
arm and rejects missing or ambiguous data.  It checks fresh starts, identical
physics/mesh/network hashes, convergence and burn-in, steady windows,
buffer/emergency state, pending trajectories/events, parent invariance,
daughter confidence intervals and profiles, carrier ESS, zero-cell fraction,
source/production/loss/escape identities, transported daughter collisions,
and variance per wall-clock cost.  Its synthetic mutation suite validates
the evaluator only; a synthetic PASS is never comet science evidence.

The campaign driver exposes separate `plan`, `run`, `resume`, and read-only
`status` operations. Planning freezes one executable/input/core/generated-
policy provenance set before launch. Every case begins at iteration zero;
resume preserves accepted cases but makes a new attempt for interrupted or
failed cases and never reuses a partial PT tree. Enhanced cases run the
focused transported-daughter collision gate before their long run. The
distributed policy template contains literal `DECISION_REQUIRED` entries,
and planning rejects it until the scientific choices below are supplied.

The following policy values remain scientific decisions and therefore are
not inferred from Stage-6/7/8 fixtures or from visual agreement with Rubin:

- minimum carrier-ESS fraction for each claimed daughter;
- minimum variance-per-wall-time gain and the required fraction of resolved
  radial bins that must achieve it; and
- the daughter flux-closure rule: preferably a seed confidence interval that
  contains zero imbalance, or otherwise an explicitly supplied relative
  tolerance.

Until those values are approved and the complete balanced ensemble passes,
Stage 9 remains open as a science stage even when all focused implementation
tests pass.

## Stage 10: documentation and master-merge preparation

### Objective and non-authority boundary

Stage 10 packages the already implemented stages into a reproducible review
bundle. It does **not** merge, switch to, update, rebase onto, or push
`master`; it does not turn focused implementation tests into a Halley science
result. Every audit names the target, merge-base, and development commits by
full SHA-1. A branch name alone is deliberately insufficient because local
and remote `master` may point at different histories.

The readiness report has two independent evidence gates:

1. the complete Stage 0--10 focused implementation suite must be
   authenticated by its machine-readable receipt; and
2. Stage-9 science acceptance requires a separate, checksum-protected real
   campaign manifest and accepted Rubin-comparison manifest.

Passing only the first gate is reported as
`FOCUSED_ONLY_NOT_SCIENCE_READY`. Passing both permits
`READY_FOR_USER_REVIEW`, which still records `merge_authorization` as
`NOT_GRANTED`. No tool in this stage runs a mutating Git command or treats a
readiness result as user approval.

### Reproducible readiness bundle

`validation/stage10_merge_readiness.py` implements two fail-closed actions:

- `collect` authenticates the exact commit topology, categorized integration
  paths, dirty-worktree exclusions, the all-stage receipt, and any optional
  Stage-9 science evidence, then creates a brand-new checksummed report
  directory; and
- `verify` recomputes the bundle checksums and validates its focused or
  science/review state without modifying the repository.

The associated `stage10_release_policy.json` is the versioned classification
policy. Source, tests, tracked inputs, and documentation are reviewable;
build trees, executables, `PT*`, transient logs, root generated input/parser
copies, meshes, and raw validation archives are excluded. Unknown paths fail
closed for human classification instead of being silently staged. Untracked
source, tests, tracked-input candidates, and root documentation also block;
only the canonical policy's explicit generated/runtime artifact families may
be summarized as untracked exclusions.

`tests/run_all_stage_tests.sh --log-root NEW_DIRECTORY` writes an
`amps.comet.all-stage-test-result.v1` receipt containing the exact HEAD commit
and tree, branch/detached state, every test and group result, log paths and
SHA-256 values, aggregate counts, archived/live execution modes, and an
explicit interpretation of focused versus Stage-9 science evidence. The
human summary remains available, but Stage 10 consumes the JSON receipt so a
copied PASS line cannot stand in for test evidence.

The tool and receipt writers are covered by mutation tests for altered logs,
hashes, commits, classifications, result counts, missing evidence, output
reuse, and attempted readiness escalation. The tracked `MASTER_MERGE.md`
describes semantic integration by subsystem and is a review guide, not an
executable merge script.

### Documentation required after every stage

- Update `srcHartley2/STATUS.md` with exact modified files, rationale, tests,
  hashes, run paths, and unresolved limitations.
- Update `srcHartley2/README.md` and `GENERIC_COMET_DESIGN.md` for user-visible
  input and algorithms.
- Update `handoff.md` with current branch, executable/input hashes, active or
  completed run state, recovery commands, and next validation action.
- Update this file when implementation evidence changes a later stage.
- Add detailed comments beside every new conservation, likelihood,
  resampling, and capacity calculation.

### Commit structure

Keep stages independently reviewable. A suggested order is:

1. diagnostics and analysis tooling;
2. input policy tables;
3. component-wide correction-free carrier planner;
4. residual accumulator;
5. core individual-weight/local-timestep compatibility;
6. exact individually weighted chemistry products and windows;
7. packet splitting;
8. channel importance sampling;
9. optional total-reaction forcing; and
10. validation artifacts and documentation.

Do not mix core mover changes with input parsing or application chemistry in
one commit. Preserve default-off behavior in every core commit. Generated
binaries, large raw PT output, transient logs, and machine-local build edits
must not be staged merely because they are present in the worktree.

### Merge order after successful validation

Prepare, but do not execute, a master integration in this dependency order:

1. reusable reaction/statistics interfaces and unit tests;
2. input parser and generated policy tables;
3. application diagnostics and correction-free carrier planner;
4. narrowly scoped core mover support, if Stage 4 was needed;
5. individual chemistry and weight windows;
6. optional importance-sampling layers;
7. application inputs and documentation; and
8. regenerated validation manifests and comparison products.

Before any master operation, present the exact development commit range,
integration diff, excluded files, test matrix, performance/memory comparison,
independent-seed statistics, convergence record, and Rubin comparison
manifest. Wait for explicit user approval.

## Recommended immediate next action

Stage 9 is implemented and focused-tested for its explicitly supported radial
configuration. The canonical correction-off, window-off, packet-off,
importance-off, total-forcing-off, diagnostics-off build has been restored.
The complete archived Stage 0--9 orchestrator passes 38 tests in ten groups,
including the 1/4/12 transported-daughter/general-NTC matrix and paired
same-binary analog/enhanced checkpoint validation. These tests close the
implementation stage, not the scientific gate. Stage 10 adds one release-
evidence group, so the final focused worktree suite passes 39 tests in 11
groups; its machine receipt remains development evidence until the Stage-10
changes are committed and retested at that exact commit. Do not merge to
`master` without a separate user decision and completed science validation.

The next scientific action is an independent-seed Stage-9 pilot.
Compare analog chemistry with candidate Stage-6 packets, complete Stage-7
conditional policies, and Stage-8 forcing using identical physical inputs and
documented independent seeds. Measure uncertainty per daughter, family ESS,
flux closure, buffer use, and runtime; verify unchanged parent lifetime and
primary profiles; and only then select production controls. The equal C2 and
H2O probabilities and broad windows in the accepted archives are validation
fixtures, not calibrated Halley/Rubin parameters. Never select packets,
windows, `q*`, or forcing species by tuning to the Rubin curve.

Before planning that ensemble, obtain explicit values for the minimum
per-daughter carrier-ESS fraction, minimum variance-per-walltime gain and
required resolved-bin fraction, plus either the recommended independent-seed
flux-imbalance confidence interval containing zero or an approved positive
relative tolerance. The shipped template intentionally contains
`DECISION_REQUIRED` and refuses to launch without those choices.
