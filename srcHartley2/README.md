# Generic AMPS comet application

`srcHartley2` retains its historical directory name so existing AMPS
application selection continues to work. Its new model, chemistry, collision,
and radial-mesh layers are generic. Hartley 2 and Halley values belong in
separate files under `input/`.

The application currently supports:

- a compiled and executed DIM1 spherically symmetric radial backend;
- the existing DIM3 triangulated-nucleus Hartley 2 backend; and
- a reserved DIM2 axisymmetric model type whose execution backend is not yet
  implemented.

`GENERIC_COMET_DESIGN.md` describes the interfaces and dimensional split.
`STATUS.md` records build and runtime verification. The Hartley 2 radial case
is `input/hartley2.radial.input`; the historical triangulated case remains
`input/hartley2.input`. The Rubin et al. case is explicitly named for comet
Halley in `input/halley1d.input`.

## Code organization

| File | Responsibility |
| --- | --- |
| `CometModel.*` | Geometry-neutral heliocentric conditions and per-species source descriptions |
| `CometChemistry.*` | Validated discrete execution of generated unimolecular reaction events |
| `CometPhysics.*` | Adapter from generated neutral collision data to the AMPS collision callback |
| `CometRadiativeCooling.*` | Species/model selection and AMPS internal-energy update for published cooling laws |
| `RadialMesh*` | DIM1 spherical domain and logarithmically expanding target resolution |
| `RadialSimulation.*` | DIM1 source injection, timestep selection, inner re-emission, outer loss, and time stepping |
| `halley2011/derive_halley1d_controls.py` | Reproducible Halley surface-resolution, injection-count, and daughter-speed derivation |
| `halley2011/digitize_rubin2011.py` and `halley2011/reference/` | Reproducible Rubin vector-curve extraction and committed neutral-profile reference data |
| `visualization/plot_1d.py` | CLI plotter for named quantities in one or more DIM1 output files |
| `visualization/compare_rubin2011.py` | Paper-ordered neutral-density and H2O-velocity comparison artifacts for the Halley run |
| `Comet.cpp` and legacy helpers | Existing DIM3 Hartley 2 implementation |

Reusable data containers live outside the application in
`src/models/chemical_solver` and `src/models/cometary_physics`. The Perl input
modules in `input/` generate species-indexed C++ tables after AMPS fixes the
`SpeciesList` ordering.

Collision data follow the same separation. `PhysicsDatabase=Rubin2011` loads
a named paper catalog, while repeated `NeutralCollision` records hold
case-specific or supplemental cited pairs. Hartley 2's former hard-coded
H2O/CO2 values are now explicit records in `input/hartley2.input`.

The Halley input separately selects
`NeutralCollisionCompletion=COMBI1996_CHEMKIN_LORENTZ_HS`. The generator
preserves the 66 Rubin parent rows, adds the 18 nonoverlapping Combi rows, and
uses CHEMKIN diameters with Lorentz mixing for only the remaining 294 pairs.
Those final rows are labeled analytical estimates in their generated source
strings; the policy is not a generic runtime fallback.

The core model-particle collision adapter must also select the user-defined
cross-section route. That route passes the current block and center cell to
`CometPhysics`, allowing the H2O-H2O Crifo law to use temperature. A generated
Halley build should contain the final override
`_PIC__PARTICLE_COLLISION_MODEL__CROSS_SECTION__USER_DEFINED_` in
`build/pic/picGlobal.dfn`.

Rubin Table 3 lists all 66 unordered pairs among the 11 parents but says that
the table is only a subset; it does not print the daughter-pair values used by
the paper model. Without an explicit completion keyword, the application
continues to abort on those missing pairs.

For the temperature-dependent H2O self-collision law, `CometPhysics` reads
`cDataCenterNode::GetTranslationalTemperature(H2O)`. That core accessor reads
the completed sample, never the sample being collected. Input must provide a
positive `CollisionTemperatureFloor`; the Halley case evaluates the law at
`max(T_sampled,180 K)`, so the first sampling window and colder later samples
use 180 K. The current DIM3 Hartley input explicitly selects the same rule.
Each temperature-law record also carries its cited validity interval. Sampled
temperatures above its maximum are evaluated at that maximum, so the Crifo
law is capped at 300 K without extrapolation. Cases can additionally set
`CollisionTemperatureMinimumSampleParticleNumber`; a cell whose completed-
window mean `Particle Number` is below that input cutoff uses the temperature
floor instead of its sampled kinetic temperature. Halley sets the cutoff to
`1.0` model particle per cell on average.

## Radial mesh

The DIM1 coordinate is cometocentric radius. The inner coordinate equals
`NucleusRadius`. Exactly one of `DomainRadialExtent` (height above the surface)
or `DomainOuterRadius` (cometocentric outer radius) must be supplied.

For requested surface width `dr0`, ideal cell count `N`, and radial extent
`L`, the mesh profile solves

```text
dr_i = dr0 q^i
L = dr0 (q^N - 1)/(q - 1),  q >= 1.
```

AMPS represents this target with dyadic AMR blocks. The allocated cell count
can therefore exceed `N`; startup diagnostics print both counts and the actual
minimum cell width.

## Radial particle backend

When `RadialSimulationMode=ON`, every compiled species requires an explicit
positive `RadialCharacteristicSpeed`. The local step is
`CourantNumber*cell_size/characteristic_speed`. This keeps energetic light
photoproducts from inheriting a parent-gas transport step.

Each direct parent source is injected once through the spherical inner face
at its total rate. Injection times are sampled as a Poisson process and
velocities come from the input temperature and radial bulk speed. Particles
leaving the outer face are deleted. Gas returning to the inner boundary is
fully thermalized and re-emitted from a half-Maxwellian at
`RadialInnerBoundaryTemperature`.

The chemistry executor samples one mutually exclusive channel per decay
event, creates all stoichiometric products from that event, and assigns
velocities that conserve event momentum and excess kinetic energy. With
individual weight correction enabled it uses one corrected model daughter per
stoichiometric product. With fixed global or species-local block weights and
correction disabled, it reads the parent and product weights from the current
block and converts the parent particle rate `W_parent/dt_parent` to each
product rate `W_product/dt_product` with unbiased residual split/roulette.
Number, mass, momentum, and energy production rates are then conserved
statistically, while event-by-event product correlations acquire Monte Carlo
variance.

### Local-timestep face transport

`LocalTimeStepFluxResampling=on` enables the core implementation used by the
Halley case. It is deliberately limited to `SpeciesLocalTimeStep`, either a
supported global weight or `SpeciesLocalParticleWeight`, correction mode off,
CPU execution, a spatial non-periodic AMR mesh, and the trace mover.
Unsupported CUDA, field-line, periodic, individual-correction, and unsupported
single-local-weight configurations fail at compile time. At every AMR block
face, the trace mover stops at the face instead of first calculating a final
position and subsequently looking up its block.

Let a particle cross from an old block with timestep `dt_old` and weight
`W_old` to a new block with `dt_new` and `W_new`. Its expected destination
multiplicity is

```text
mu = (W_old/dt_old) / (W_new/dt_new)
   = W_old*dt_new/(W_new*dt_old).
```

The mover creates `floor(mu)` particles plus one with probability
`mu-floor(mu)`. This is cloning for `mu>1` and Russian roulette for `mu<1`, and
therefore matches the outgoing and incoming number, mass, momentum, and
kinetic-energy fluxes in expectation without changing velocity. The current
Halley initializer sets each species' block weight proportional to its local
timestep. Thus `W_new/W_old=dt_new/dt_old` and ordinary same-species radial
transport has `mu=1`; the full expression remains necessary for arbitrary
neighboring weights and is tested independently of the timestep ratio.

If `tau_old` is the unused portion of the old block step, continuation uses

```text
tau_new = tau_old * dt_new/dt_old,
tau_new/dt_new = tau_old/dt_old.
```

The dimensionless unfinished fraction is stored in an optional particle-data
field. AMPS commits local lists and first exchanges particles stopped by the
main mover. Each continuation sweep then advances only unfinished particles,
unconditionally exchanges all particles deposited in remote boundary blocks,
and only then performs the global pending-work reduction. The post-move
exchange must not be conditional on a nonzero remainder: a particle may finish
its continuation in a remote boundary block and still need ownership transfer.
Collisions, chemistry, sampling, and generic particle splitting run only after
all trajectory segments are complete. This ordering prevents an interface
crossing from receiving more than one physical collision/chemistry pass in an
outer iteration.

The algorithm is a steady-state local-pseudo-time acceleration, not a
synchronized transient integrator: different blocks advance different
physical intervals during one outer iteration. Block-local weights remove the
global-weight population factor across the logarithmic domain: `W` grows with
`dt`, keeping the expected steady-flux model population approximately uniform
per radial cell. This does not turn the algorithm into synchronized physical
time, and a long population/convergence run remains an acceptance requirement.

Direct-parent surface weights retain the requested injection statistics and
are scaled through the mesh at constant `W/dt` for that species. Species with
no direct source use the largest direct-parent `W/dt` carrier rate to avoid
large direct-parent-to-daughter clone bursts. Reactions between species with
different carrier rates still perform
`(W_parent/dt_parent)/(W_product/dt_product)` split/roulette explicitly.

### Sampling-break block-weight adaptation

The radial input can independently enable population control for the
species-local block weights:

```text
RadialParticleWeightAdaptation=on
RadialParticleNumberPerCellRange=50 | 200
RadialEmergencyParticleBufferUtilization=0.75
RadialAllocationLowWatermark=0.50
RadialAllocationReserveFraction=0.10
RadialMaxFaceTransferMultiplicity=1
RadialMaxPhotochemicalProductsPerReaction=20
```

Both bounds are inclusive positive integers. The maximum must be at least
twice the minimum; that condition guarantees that an integral power-of-two
change can place any nonzero population in the requested band. The switch is
mandatory for a radial case so that adaptation is never enabled implicitly;
the range is required only when the switch is on. When chemistry and
adaptation are both enabled, the positive-integer photochemical-product limit
is also required. The emergency high watermark, allocation low watermark,
reserved fraction, and face-transfer limit are required whenever adaptation
is on. The low watermark must be below the high watermark, and the high
watermark cannot enter the reserved tail of the particle buffer. The currently
validated correction-free face controller requires a transfer limit of one.
These are numerical controls, not Rubin et al. model parameters.

Adaptation runs only after a complete sampling window has been printed and
used for the convergence comparison. It therefore cannot alter the sample
that triggered it. For each owner block and each species separately, let `N`
be the instantaneous number of model particles and `Nc` the number of cells
in that block. The controller uses the whole-block mean

```text
nbar = N/Nc.
```

Cells with no particles contribute zero to this mean. A block/species pair
with `N=0` is skipped by the population-range controller, and all population
resampling stays in the particle's current cell, so an empty cell is never
populated by adaptation.

The target weight is restricted to `W_new=W_old*2^k`. The controller chooses
the integer `k` that moves the expected post-adaptation mean into the input
range. Each factor-of-two step is applied as follows:

- `k<0`: set `W_new=W_old/2` and clone every current particle of that species
  exactly once. Repeating this operation implements larger negative powers.
- `k>0`: set `W_new=2*W_old` and retain each current particle independently
  with probability `W_old/W_new=1/2`. Repeating this operation implements
  larger positive powers.

Splitting conserves represented number, mass, momentum, energy, and internal
state exactly. Roulette conserves those quantities in expectation and adds
the corresponding Monte Carlo variance. A split improves later transport and
collision sampling only after the siblings decorrelate; copying a particle
does not create independent information at the adaptation instant. Clones
carry the complete particle record but receive a new tracker identity when
tracking is enabled. Before
changing any lists, the rank computes the exact worst-case clone demand and
aborts if the configured particle buffer cannot hold it; this prevents a
partially applied adaptation.

Chemistry allocation requires a hard bound, including products that are
already present. After planning ordinary population changes, the controller
starts from every species present in a block and every configured direct
surface source. It closes the enabled unit-unimolecular reaction graph
transitively, because a daughter can be created and react again before the
next sampling break. Including direct sources lets the same proof run during
initialization, before the first particle is injected.

For one selected reaction channel, a concrete daughter's residual-resampling
mean is

```text
mu = (W_parent/dt_parent)/(W_daughter/dt_daughter).
```

The actual allocation may be `ceil(mu)` because the residual Bernoulli draw
can round upward. For every reachable channel the controller therefore raises
product weights by powers of two until

```text
sum_products coefficient * ceil(C_parent/C_product) <= M,
```

where `C=W/dt` and `M` is
`RadialMaxPhotochemicalProductsPerReaction`. The reachability/carrier loop
continues to a fixed point, so the bound covers later-generation daughters and
all products of a selected channel, whether present or absent. No weight is
lowered and no empty cell is seeded. After transport constraints have been
applied, the hard chemistry bound is evaluated again with the final carriers;
a violation aborts before chemistry resumes.

Transport has a separate hard allocation bound. Its residual multiplicity at
a block face is

```text
mu_face = (W_source/dt_source)/(W_destination/dt_destination).
```

This ratio can be large for either a populated or empty destination. With
`RadialMaxFaceTransferMultiplicity=1`, both directed senses of every face must
satisfy `mu_face<=1`, so every block in a face-connected component must have a
common carrier for that species. After population and chemistry plans are
applied, the controller:

1. gathers authoritative owner counts and carriers for every calculation leaf;
2. constructs face adjacency using all `2*DIM` faces and the active AMR
   subface slots, without introducing edge or corner coupling;
3. finds each face-connected component and selects a common carrier reachable
   from every member by power-of-two increases;
4. applies unbiased factor-two roulette to populated blocks being raised,
   while empty blocks receive only a weight change and remain empty; and
5. verifies both directed senses of every face and aborts unless the maximum
   transfer multiplicity is at most one to numerical tolerance.

Because all adaptive changes are dyadic, no unrelated statistical scale is
introduced. The standard block-data exchange then distributes owner values to
ghost blocks before transport resumes. Diagnostics report changed
block/species pairs, component count, maximum graph distance, and maximum face
clone mean before and after enforcement.

For a future multidimensional comet calculation, the definition must remain
face connected and all active coarse/fine subfaces must be tested. The
implementation has the required DIM-independent enumeration, but DIM2/DIM3
refinement-interface tests are required before claiming validation outside
DIM1. On the connected radial mesh, a limit of one makes `W/dt`
component-global for each species; local weights still vary with local
timestep. This deliberate variance/capacity tradeoff is stronger than the
earlier empty-neighbor-only repair.

Only owner blocks are resampled. After any change, the normal AMPS block-data
exchange propagates the new weights to MPI boundary copies, and the legacy
global-minimum weight cache is recomputed because adaptation can increase as
well as decrease every local weight. Subsequent block-face transport and
species-changing chemistry already use the current source and destination
block weights, so their existing `(W/dt)` balance consumes the synchronized
carrier plan without a second correction mechanism.

### Emergency pre-exchange buffer protection

Sampling-break adaptation alone cannot protect a rank from population growth
inside a 500-iteration window. The 12-rank validation ending at iteration
5,690 demonstrated this: the global population grew to 29,159,492 after the
iteration-5,500 adaptation and rank 7 exhausted its 20,000,000-slot local
buffer before the next sampling break.

The core therefore provides an optional application callback immediately
before `PIC::Parallel::ExchangeParticleData()` allocates receive-side particle
slots. Generic AMPS cases leave the callback null. The radial comet
initializer installs it only when block-weight adaptation is enabled. Every
ordinary exchange, including local-timestep continuation exchanges,
collectively evaluates the rank-buffer trigger

```text
u_max = max_rank(N_allocated/MaxNPart).
```

Completed sampling breaks provide a second trigger without adding another
particle-list traversal. After counting every species for the ordinary
block-weight controller, each owner block evaluates

```text
mean_total(block) = sum_species N(block,species) / N_cells(block)
sampling_limit = 10 * nTotalSpecies * MaximumParticlesPerCell.
```

The maximum is reduced across MPI ranks. A strict exceedance records the same
pending request on every rank; it does not run roulette recursively from the
sampling routine. The next ordinary particle exchange consumes that request
in the established pre-exchange callback. The count is taken before regular
per-species sampling-break adaptation, so a transient population excursion
cannot be hidden by the ordinary controller before it is recorded. With the
Halley settings this limit is `10*27*200 = 54,000` total model particles per
cell as a block mean.

If neither the rank-buffer threshold nor a pending sampling-population request
is active, no state is changed. Either trigger causes one emergency
factor-two step:

1. retains every allocated model particle independently with probability
   `1/2`, including particles already staged in non-owner MPI boundary-block
   copies;
2. doubles every species weight in every allocated owner or boundary block;
3. reconstructs the core's global-minimum weight cache; and
4. consumes a pending sampling-population request; and
5. rechecks the maximum rank utilization, repeating collectively only if
   stochastic survival did not create the requested buffer headroom.

Thus a sampling-population trigger guarantees one factor-two step. The 75%
rank-buffer trigger may require additional steps until every rank is below its
configured threshold. Diagnostics distinguish `buffer`,
`sampling-population`, and combined triggers and retain the sampling iteration,
observed block mean, and threshold.

Applying one common factor to all species and blocks is deliberate. It
preserves every adjacent-block weight ratio and every inter-species `W/dt`
carrier ratio, so face-flux conversion and chemical residual sampling are
unchanged. Roulette conserves represented number, mass, momentum, translational
energy, and internal energy in expectation. It never splits particles or
populates an empty cell, and it is safe for a survivor with an unfinished
local-timestep continuation because the complete particle record is retained.
The price is an emergency increase in Monte Carlo variance; the ordinary
per-block 50--200 controller may later restore resolution at a completed
sampling break.

The ordinary threshold protects MPI receive allocation. Allocation-producing
chemistry adds a stronger path: it stops at collective batch checkpoints when
the 0.75 high watermark is reached, requests the same common-factor roulette,
and repeats factor-two coarsening until every rank is at or below the 0.50 low
watermark. The last 10% of each particle buffer is never offered to a chemistry
transaction. Transport cannot create a face clone when the face limit is one,
so it cannot produce the earlier first-arrival allocation burst.

### Transactional photochemistry

The generic AMPS photolysis driver exposes a nullable five-callback batch
interface. Applications that leave all callbacks null retain the historical
single pass. The comet application installs the complete set and reserves one
byte in each particle record before the particle buffer is initialized. At the
start of each chemistry stage that byte is cleared on existing particles.
Completed parents and all daughters created during the stage are marked, so a
new daughter cannot react a second time in the same outer AMPS iteration.

Each batch conservatively charges every attempted parent the configured hard
maximum of 20 daughter slots and stops before crossing the 0.75 high watermark.
For an event that reacts, all residual multiplicities are drawn before any list
mutation. If the complete event would cross the 90% transaction ceiling, the
unchanged parent is reinserted and retried after the collective checkpoint.
Otherwise the existing parent slot is reused for the first daughter and only
the remaining daughter slots are allocated. An actual event above the hard cap
aborts because truncation would bias chemistry. This callback mode is currently
MPI-only; a hybrid/OpenMP build terminates with an explicit unsupported-mode
error instead of entering unsafe shared counters.

The focused regression is selected instead of a science run:

```sh
mpirun -np 1 ./amps --validate-weight-adaptation
mpirun -np 4 ./amps --validate-weight-adaptation
mpirun -np 12 ./amps --validate-weight-adaptation
```

It exercises a multi-step exact split, factor-two binomial roulette, owner
blocks on every rank, and preservation of initially empty cells. It creates a
two-sided unsafe face ratio and verifies component-wide carrier equalization
without clone allocation. It independently exercises the buffer,
allocation-burst, and sampling-population triggers. For the sampling trigger
it constructs an exact synthetic
`10*nTotalSpecies*maximum*block-cell-count` threshold, exceeds it, verifies
that the sampling routine records but does not execute the request, and then
invokes the regular core particle-exchange API. The test verifies that the
emergency callback runs before migration, consumes the request, applies a
factor-two/binomial step, lowers maximum-rank occupancy for the independent
buffer test, and keeps owner/boundary weights synchronized. The allocation
test requires two consecutive weight doublings, verifies the low-water target,
request consumption, and statistical survival. With chemistry enabled it
forces one channel above the actual hard allocation cap, verifies transitive
reachability and direct-source initialization, confirms the bound falls to at
most 20 without seeding particles, and checks MPI boundary-copy
synchronization. It deletes its artificial population and must not be invoked
during a science calculation.

## Internal energy and radiative cooling

The generic input record pairs a compiled species with a named physical law:

```text
RadiativeCooling=H2O | GOMBOSI1986_H2O
```

`GOMBOSI1986_H2O` is restricted to H2O during preprocessing. The reusable
evaluator in `src/models/cometary_physics/RadiativeCooling.h` implements the
Crovisier/Huebner water rotational-cooling approximation printed by Gombosi
et al. (1986), equation (48b). `CometRadiativeCooling` traverses local AMPS
cells after each `PIC::TimeStep()`, evaluates the selected rate from the most
recent completed H2O density and rotational-temperature sample, and subtracts
the corresponding energy from each H2O particle's rotational state. It is
dimension-independent and is called by both DIM1 and DIM3 mains.

The Halley input enables the core Larsen-Borgnakke model, H2O's three
rotational degrees of freedom, `Zrot=5`, RT relaxation, and a hard-sphere
temperature index of 0.5 for every H2O collision partner. Surface injection
and thermalized return initialize internal energy through the core IDF API.
VT and VV relaxation are disabled because the present core VT routine exits
as untested and unverified; the case does not claim vibrational relaxation.

## Build and verification

AMPS preprocessing uses the application input and parser copies in the
repository root. After selecting the desired input in `.amps.conf`, ensure
those generated root copies reflect `input/`, then remove `build` before every
compile:

```sh
rm -rf build
make -j8
mpirun -np 1 ./amps
```

The DIM1 executable also contains a transport-only validation mode. It uses
the configured radial AMR mesh but does not start the comet source, chemistry,
collision, cooling, or sampling loop:

```sh
mpirun -np 1 ./amps --validate-local-timestep \
  --validation-particles 200000 --validation-mover both \
  --validation-skip-multiface

mpirun -np 4 ./amps --validate-local-timestep \
  --validation-particles 200000 --validation-mover both
```

`--validation-mover` accepts `first`, `second`, or `both`. For local-weight
builds, the single-face matrix tests `(dt ratio, weight ratio)` pairs
`(0.25,1)`, `(0.5,2)`, `(1,0.4)`, `(2,0.5)`, `(2.5,2.5)`, `(0.5,0.25)`, and
`(2,4)` in both radial directions. These give `mu` values from 0.25 through 4
and include the production cancellation case `dt ratio = weight ratio`.
Selected independent-ratio, clone, roulette, and cancellation cases are
repeated across an MPI ownership face. Four velocity groups include radial
and transverse components. Every case checks particle-list membership, final
block, finite state, unchanged speed, analytic radial position, zero remainder,
continuation lifecycle counts, and number/mass/three-momentum/kinetic-energy
flux. A stochastic moment is accepted inside the larger of a 2% relative band
and its five-standard-deviation Monte Carlo band, so a near-zero transverse
moment is not judged by an ill-conditioned relative denominator.

The multi-face case follows 32 particles across at least three faces and, on
three or more ranks, at least two ownership boundaries. It requires the number
of continuation sweeps to equal the number of faces, no surviving remainder,
the expected transfer count, and conservation of every recorded flux moment.
This mode is an algorithm regression; the full comet acceptance gate still
requires a long production run through sampling, load redistribution, and a
converged Rubin comparison.

Launch the Halley executable from the repository root. Its configured output
directory is `PT/plots`, relative to that working directory, so AMPS creates
both `/data/vtenishe/COMET1D/AMPS/PT` and its `plots` subdirectory during
runtime initialization. Compilation alone does not create run directories.

Parallel DIM1 runs require the dimension-aware neighbor-mask bounds in
`src/pic/pic_block_send_mask.cpp`. The AMR exchange visits two faces with one
slot per face in DIM1, four faces with two slots per face in DIM2, and six
faces with four slots per face in DIM3. This prevents lower-dimensional runs
from indexing the three-dimensional neighbor layout during MPI exchange.

After decomposition, only the rank that owns the nucleus-adjacent block can
read its surface weight and timestep. `RadialSimulation.cpp` verifies that
there is exactly one owner and distributes those values with `MPI_Allreduce`
before all ranks initialize their local blocks. Sampling-break redistribution
uses dimension-aware block packing in `src/pic/pic_mesh.cpp`: DIM1 serializes
only its positive radial face, while DIM2 and DIM3 retain their two and three
positive faces. These bounds are required when particle load balancing moves
radial blocks between MPI ranks.

The checked-in Hartley 2 radial input leaves `RadialSimulationMode=off`; its
runtime validates configuration and builds the mesh without inventing missing
science-run controls. The Halley file records the supplied 2,000 cells,
Courant number 0.25, 400 m/s parent baseline, disabled methanol branch, and
controls derived by `halley2011/derive_halley1d_controls.py`. The 500-step
sampling window and 10,000-iteration hard ceiling are explicit. The selected
0.36 eV OH branch uses the lower published rate, `3.0e-6 s^-1` at 1 AU. The
run compares the volume-weighted relative L1 change of every species' complete
radial density profile and terminates after all species remain within 5% for
three consecutive window-to-window comparisons. The first complete window is
the reference, so convergence needs at least four windows. Multi-rank tests
have crossed output and block-migration boundaries, verifying runtime plumbing.

`RadialLoadBalanceInterval` is a separate numerical control and does not
change the 500-step science window. After a complete outer iteration, the
radial driver may redistribute whole AMPS blocks through the core's standard
`CreateNewParallelDistributionLists()` path, reset the core rebalancing
accounting, and rebuild the rank-local block table. Redistribution is skipped
when the same iteration already triggered the core sampling-break path.
Halley uses an interval of 50 iterations: an uninterrupted 18-rank validation
with sampling-only redistribution reached about 12 million global particles,
but its outward-moving front filled rank 17's five-million-slot local buffer
at iteration 2344. Frequent redistribution changes only MPI ownership;
particle state, statistical weight, local timestep, accumulated samples,
chemistry, cooling, and internal energy move with each block. The buffer stays
at five million particles per rank so the corrected run tests the ownership
correction instead of masking the imbalance with additional memory. That
18-rank run passed the old failure point and reached iteration 2699. It wrote
five complete windows through iteration 2500, but none of the four profile
comparisons met the 5% convergence tolerance. Population grew from 49,983 at
the iteration-500 output to 26,220,686 at iteration 2500 and 55,330,093 in the
iteration-2696 core diagnostic. The attempted iteration-2700 redistribution
then filled rank 10 while migrating its new whole-block set. Scheduled
balancing therefore fixes premature ownership localization, but does not make
the species-global-weight population stable over the full logarithmic Halley
domain. The last census was 99.14% direct-parent particles, so this growth is
not a chemistry clone cascade. The detailed population history, 6,165.38 s
wall time, and transient Rubin comparison ratios are recorded in `STATUS.md`.
The corrected long run used the selected complete collision hierarchy. A
later audit of the earlier smoke/test binaries cited below found that they had
used AMPS's zero placeholder cross section because the Rubin callback was not
connected. That dispatch is now fixed. A new clean four-rank diagnostic used
the 180 K first-window rule,
completed the first 500-iteration sample, applied cooling from that completed
sample during subsequent iterations, and was deliberately stopped at
iteration 525. Its output is transient validation data, not a converged
physical Halley result and predates the selected complete collision hierarchy.
The earlier 0.5-Courant convergence run reached iteration 3500 but a fast
collision product then traversed beyond the mover's one-block MPI boundary
allocation. The input now uses 0.25; its clean four-rank restart crossed the
first two output/load-balance boundaries. Exact evidence, live output paths,
and remaining comparison limits are recorded in `STATUS.md`.

That global-weight run was explicitly stopped at iteration 2729 after
reaching 62,087,575 particles. The replacement configuration uses block-local
weights proportional to local timestep. Its clean focused regression passed
all first-/second-order serial and four-rank MPI face cases, including
independent weight/timestep ratios and a 95-face path through two ownership
transitions. The four-rank chemistry ensemble also passed all 46 neutral
reactions. A production smoke run completed outputs at iterations 500 and
1000 with 40,841 and 78,754 particles, respectively, versus 51,576 and 172,215
for the superseded global-weight run. In the H2O output, block weights span a
factor of `3.355443e7` while `W/dt` is constant to `5.3e-7` relative at printed
precision. Exact logs and paths are recorded in `STATUS.md`.

The subsequent four-rank full-domain run completed all 10,000 iterations in
approximately 3 h 18 min 41 s with 628,490 final particles and no buffer,
mover, chemistry, or MPI failure. It produced all 20 sample windows and the
complete Rubin PNG/EPS/PDF/CSV/manifest comparison set from output 19.
Population stability therefore passes for block-local weights. The run is
not converged, however: the final all-species relative L1 maximum is
`1.719734` for sparse CH with `0/3` required consecutive windows. The
comparison must remain labeled transient; paths, hashes, per-window metrics,
and selected primary-species ratios are in `STATUS.md`.

These focused results still do **not** authorize merging the full branch into
`master`. `STATUS.md` records the gate-by-gate disposition, independently
reviewable subsets, file dependency order, conflict checks, clean-commit
requirements, and the integration-branch procedure. No merge should begin
until a documented, statistically defensible convergence policy passes and
the resulting steady Rubin comparison is reviewed. The block-local-weight
run has demonstrated bounded population and no buffer exhaustion, but it did
not satisfy the current three-consecutive-window 5% criterion.

## Hot-loop indexes and allocation reuse

The Halley production profile identified two generic database scans and
per-cell/per-particle allocation as avoidable costs. The application and core
now use three coordinated optimizations:

- `cCollisionDatabase` constructs a symmetric dense pair index, and
  `CometPhysics` evaluates the pointer returned by the first lookup instead of
  searching the 378-record table a second time;
- `cReactionNetwork` constructs a per-parent unit-unimolecular index, while
  `CometChemistry` caches ordered channel rates and their total once per
  heliocentric-distance state; and
- the general core NTC routine builds reusable per-thread species buckets in
  one cell-list traversal instead of allocating scratch arrays and rescanning
  the linked list for each species.

These changes preserve input reaction order, mutually exclusive channel
selection, collision pair ordering, the NTC 100-probe majorant, random-number
call order, and all block-local weight/timestep corrections. The chemistry
cache is refreshed during validation and at a batched-stage boundary if the
heliocentric distance changes. It contains only positive unit-unimolecular
rates; arbitrary non-unimolecular records remain in the generic reaction
network for other solvers.

The clean 14-rank regression under
`validation_runs/halley1d_priority_optimization_preflight_20261104/benchmark14`
completed the full 10,000 iterations. All 540 AMPS `.dat` files are
byte-for-byte identical to the pre-optimization baseline. The timestamp span
from the first output file to the last fell from 3,052 s to 470 s, a 6.49x
end-to-end speedup for this host/run configuration. The final population
(391,034), convergence metric (1.946230, limited by CH3), and 0/3 passing
windows are unchanged. This is a performance/equivalence result, not a claim
that the science case converged.

A post-change Callgrind sample is retained in
`validation_runs/halley1d_priority_optimization_preflight_20261104/callgrind_optimized`.
It contains no linear collision-table scan. Temperature-law `pow` evaluation
is about 2.2% of sampled instructions, so no cell/state cross-section cache was
added; such a cache would need an explicit invalidation contract and is not
justified by the residual cost.

## Visualizing DIM1 output

`visualization/plot_1d.py` reads column names directly from AMPS Tecplot
output. The radial X axis is logarithmic; the CLI requires a physical quantity
and an explicit Y scale. Multiple input files are drawn on the same axes, and
each legend entry contains its source filename. For example:

```sh
python3 srcHartley2/visualization/plot_1d.py \
  PT/plots/pic.H2O.s=13.out=0.dat PT/plots/pic.CO2.s=26.out=0.dat \
  -q dens --yscale log --radius-unit km --rrange 6 100 \
  --output PT/plots/parent-density
```

Run with `--list-quantities` to print the exact physical quantities present in
a file. `--help` lists the short quantity aliases, including `dens`, `temp`,
`speed`, `vel`, and all species-specific collision frequencies, followed by
copyable usage examples. Each saved plot is written as both PNG and EPS and
opens in the desktop image viewer; use `--no-show` for batch execution.
`--rmax` selects a maximum cometocentric distance and `--rrange MIN MAX`
selects an inclusive interval in the chosen radius unit. Non-positive points
are skipped, with curve gaps, when logarithmic Y scaling is selected.
Additional examples are documented in `visualization/README.md`.

For the Halley paper case, generate the latest complete set of neutral plots,
the combined PDF, and the numerical comparison table with:

```sh
python3 srcHartley2/visualization/compare_rubin2011.py \
  PT/plots --run-log halley1d.run.log
```

This overlays AMPS results with vector-extracted Rubin neutral DSMC curves for
Figs. 1--17 and the modeled H2O radial-velocity curve from Fig. 18. The default
reference table and provenance metadata are in `halley2011/reference/`; the
reproducible extractor is `halley2011/digitize_rubin2011.py`. The paper's
Giotto NMS mass-channel curves require an instrument-response model that is
not part of this application and are not mixed with the neutral-density
comparison.
