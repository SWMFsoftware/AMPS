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

### Existing support and blockers

- The particle buffer already stores and clones an individual correction.
- Macroscopic sampling multiplies the base weight by that correction.
- The general NTC collision implementation includes individual corrections in
  collision majorants and unequal-weight update probabilities.
- Comet chemistry already reads a corrected parent and can assign one
  corrected daughter per concrete product.
- The local-timestep flux-resampling mover has a compile-time prohibition.
- `RadialSimulation::Validate()` rejects adaptive block weights with
  individual correction enabled.
- Simplified NTC assumes corrections are off, and unsupported collision
  implementations must continue to fail explicitly.

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

## Stage 5: exact individually weighted daughters and weight windows

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

This is the formula already present in `CometChemistry.cpp`. The work in this
stage must make it reachable under the validated radial configuration and
remove fixed-weight constraints that no longer describe actual allocation.

With one individual particle per concrete product, the hard transaction count
is the stoichiometric product count rather than
`ceil(C_parent/C_product)`. Keep a separate input hard limit and preflight the
complete transaction before modifying the parent.

### Weight windows

Individual corrections must not drift over an unbounded range. Add a
species/group input policy defining lower, target, and upper correction or
physical-carrier limits.

For `c<c_target`, roulette to the target with survival probability
`c/c_target`. For `c>c_target`, split into a bounded number of particles whose
corrections sum exactly to `c`. Apply windows only at safe completed stages,
not during an unfinished mover continuation.

Every window operation must preserve position, velocity, internal state, and
chemistry stage. Split copies receive unique tracker identities. Correlated
copies do not qualify as independent samples until a measured burn-in.

### Tests

- Extend the chemistry ensemble to span correction ratios below and above one
  over every reaction in the generated network.
- Verify weighted number, mass, momentum, and excess-energy production for
  each channel and the complete ensemble.
- Verify repeated stoichiometric products and neutral sinks.
- Test weight-window mean/variance, split-sum exactness, roulette expectation,
  capacity failure, and restart.
- Confirm that direct-source particles retain correction one unless another
  explicitly selected model changes it.

### Pilot acceptance criteria

- Every selected C2 channel creates a model C2 product.
- The weighted C2 production rate agrees with the analog physical rate within
  its predicted uncertainty.
- Parent reaction probability, lifetime, density, velocity, and injection
  statistics are unchanged within measured uncertainty.
- Correction distributions remain inside the configured windows after the
  controller runs.

## Stage 6: daughter packet splitting

### Objective

Create several statistically useful trajectories from one minor-product
cohort while preserving its total individual weight.

### Algorithm

For a daughter correction `c_d` and requested packet count `N`, create `N`
particles with corrections summing exactly to `c_d`. Use equal corrections
unless a documented variance-minimizing allocation is introduced later.

For a multi-product dissociation, generate `N` complete reaction packets. In
each packet, sample one recoil solution and create every stoichiometric
product from that same solution. Different packets receive independent recoil
directions. This preserves momentum and excess energy within each packet and
provides independent velocity information; cloning one velocity `N` times
does not.

Packet selection must be input-driven by species, reaction ID, or a generic
chemical-only group. Enforce a hard packet/transaction limit and include it in
the chemistry batch capacity calculation.

### Tests

- Verify per-packet momentum and energy conservation.
- Verify that the sum of packet corrections equals the unsplit correction.
- Compare velocity and density variance for `N=1,2,4,8` over independent
  seeds.
- Confirm no change in weighted means as `N` changes.
- Profile memory and runtime; report the marginal variance reduction per
  added daughter particle.

### Exit criteria

- Packet count is chosen from measured variance/cost, not from agreement with
  the paper curve.
- Primary populations remain unchanged.
- C2 and other selected minor species meet their input statistical target with
  acceptable buffer headroom.

## Stage 7: reaction-channel importance sampling

### Objective

Reduce variance from rare channel selection while leaving the total parent
reaction probability analog.

### Algorithm

If physical channel probability is `q_j` and the biased selection probability
is `q_j*`, multiply every product packet from that channel by

```text
L_j = q_j/q_j*.
```

For `N` equal packets, the correction is

```text
c_daughter,i =
    [C_parent/(W_daughter/dt_daughter)] * L_j/N.
```

All enabled channel probabilities must remain normalized and positive where
the physical rate is positive. A channel with zero physical rate cannot be
made nonzero statistically. The input must refer to stable reaction IDs, not
generated array indices.

An optional deterministic branching mode may instantiate every neutral
channel for one analog parent reaction, with each branch weighted by its
physical probability. Treat sink/ionization channels as weighted loss. This
mode has a larger transaction size and requires its own capacity bound.

### Tests

- Validate likelihood ratios and normalization for two-channel, many-channel,
  dominant-channel, and very small-channel networks.
- Compare analog and biased ensembles for product number, mass, momentum,
  energy, and parent loss.
- Prove that changing `q_j*` changes variance but not weighted mean.
- Test reaction-ID parsing, duplicates, missing IDs, and biasing a disabled
  channel.

### Exit criteria

- Parent total reaction rate and lifetime remain analog and unchanged.
- Weighted channel yields agree with catalog branching fractions.
- Bias factors and resulting effective sample sizes are reported in output
  metadata.

## Stage 8: optional total-reaction forcing

### Objective

Address a species that remains unresolved because the total parent reaction
probability—not merely the channel or product carrier—is too small.

This stage is optional and must not begin unless Stages 2 through 7 are
insufficient and diagnostics identify rare total reaction occurrence as the
limiting variance.

### Required physics

If physical reaction probability is `P` and a biased probability `P*` is
used, a reacted branch carries `P/P*`, while a sampled survivor must carry

```text
(1-P)/(1-P*).
```

Removing parents at `P*` without the survivor correction is prohibited because
it biases parent lifetime and every downstream density.

Prefer forced-reaction segments or implicit branching over per-step
uncontrolled splitting:

- retain an unreacted branch with the correct survival weight;
- create a reacted branch with the complementary physical weight;
- sample reaction position from the conditional optical-depth distribution;
- distribute reacted weight among channels; and
- immediately apply weight windows and capacity limits.

### Tests

- Pure exponential-decay tests with known analytic survival and product
  distributions.
- Spatially varying rate tests with conditional reaction positions.
- Multiple-generation networks with biased and unbiased parents.
- Collision/transport tests proving weighted parents retain the correct
  lifetime and feedback.
- Direct comparisons among analog, forced, and analytic solutions over many
  independent seeds.

### Exit criteria

- Parent survival, reaction-position distribution, and daughter production
  pass analytic tests.
- Variance reduction justifies the additional branching and collision
  complexity.
- The feature remains opt-in and is not required by ordinary comet cases.

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

## Stage 10: documentation and master-merge preparation

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

Stages 0--2, the Stage 3 implementation, and the permanent scientific gate are
implemented and focused-tested. Review the read-only `plan` output, then run
the complete eight-case matrix with `validation/run_stage3_science_gate.py`.
Evaluate the first complete archive, review and explicitly accept its Rubin
candidate baseline, and rerun the evaluator with that accepted baseline. Use
the resulting outcome to skip Stage 4, consider Stage 4, or repeat an
inconclusive experiment. Do not merge to `master` or start Stage 4 before that
decision is recorded.
