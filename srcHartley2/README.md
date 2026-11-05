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
Halley in `input/halley1d.input`. `srcHartley2/CHEMISTRY.md` provides
the staged implementation and validation roadmap for improving minor
photochemical-daughter statistics without increasing direct-parent injection.

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
| `validation/` | Frozen Stage-3 science policy, durable eight-run Bernoulli/Systematic runner, read-only decision evaluator, and operating instructions |
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

The Stage-2 validator places sparse populations in one owner block per rank
and at both ends of the 373-block radial chain. It verifies that one component
carrier is selected from the mean over all 2,984 cells, including empty cells,
while empty cells remain empty. Exact splitting conserves represented number,
mass, three momentum components, kinetic energy, and internal energy. The
roulette case checks the same moments against the exact weighted Bernoulli
variance, including unequal local block weights.

The validator also forces a zero-capacity preflight rejection and compares
particle counts, moments, and exact carrier bits before and after planning;
the rejected transaction must make no partial mutation. Direct-source carrier
bits remain locked. With chemistry enabled, every reachable component channel
is solved before commit and the focused unsafe channel must fall to at most 20
products per event. Owner values are then checked in MPI boundary copies. The
earlier buffer, allocation-burst, and sampling-population emergency tests
remain mandatory regressions. This mode deletes its artificial population and
must never be invoked during a science calculation.

The geometry-independent planner suite covers one block, long chains,
disconnected components, empty and two-ended populations, exact range
boundaries, locks, capacity rejection, and synthetic DIM2/DIM3 face graphs:

```sh
srcHartley2/tests/run_stage2_planner_tests.sh
```

After producing the 1/4/12-rank runtime logs, run the strict cross-rank gate:

```sh
HARTLEY2_STAGE2_LOG_1=validation_runs/stage2_component_carrier_20261105/mpi1/run.log \
HARTLEY2_STAGE2_LOG_4=validation_runs/stage2_component_carrier_20261105/mpi4/run.log \
HARTLEY2_STAGE2_LOG_12=validation_runs/stage2_component_carrier_20261105/mpi12/run.log \
  srcHartley2/tests/run_stage2_integrated_tests.sh
```

The archived 2026-11-05 result passes with maximum face multiplicity 1,
maximum reaction allocation 20, zero partial mutations, and respectively
0/28/84 non-owner ghost checks at 1/4/12 MPI ranks. These transaction tests do
not by themselves establish convergence of the Rubin science comparison.

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

## Stage 0 chemistry-observability tests

The minor-daughter roadmap is in `srcHartley2/CHEMISTRY.md`. Stage 0 is
complete on `1d_comet_dev`; run its focused and archived tests with:

```sh
srcHartley2/tests/run_stage0_tests.sh
```

The suite compiles the geometry-independent accumulator in
`src/models/chemical_solver/ChemistryDiagnostics.h`, verifies one-, two-,
repeated-, and zero-product reaction accounting, and reduces an identical
physical event set at 1, 4, and 12 MPI ranks. The reduced count and physical
moment signatures must be identical at every rank count. The disabled path
must retain no records and emit no serialized fields.

Production instrumentation is selected in the comet input with:

```text
CometDiagnostics=ON
```

The default and canonical Halley setting is `OFF`. The off path returns before
allocating diagnostic tables, calling MPI, writing files, changing model
state, or drawing random numbers. With it enabled, each completed sampling
window writes `PT/diagnostics/comet_stage0.out=N.csv` and `.json` plus one
run-level `PT/diagnostics/cell_measures.csv`. The latter contains the exact
ordered AMR cell-center radius and core `Measure()` value used by integrated
convergence; it is not reconstructed from radial coordinates.

An eligible parent visit means one visit to one enabled, positive-rate
reaction channel for the current parent particle. Reaction rows report
selection, destruction, sink, deferral, and emergency counts. Product rows
retain each stoichiometric term separately and report zero/one/many residual
multiplicity, allocated/reused/model carriers, and expected and represented
number, mass, momentum, and energy rates. Represented-number variance is used
to report `sigma`; expected/represented production passes when
`abs(represented-expected) <= max(5*sigma, 1%*abs(expected))`.

Species flux columns are physical number rates in `s^-1`. Each event adds its
represented `W/dt`; the collective sum is divided by the completed sampling
length. Surface injection, expected/represented chemical production,
chemical destruction, gross inner-boundary incident and thermal re-emitted
rates, and outer escape are distinct. The two inner terms are intentionally
separate even though the current thermalizing boundary gives zero net number
flux. Momentum columns are newtons, mass rates are `kg s^-1`, and energy rates
are watts.

`visualization/analyze_minor_species.py` creates JSON and CSV reports from one
or more saved species windows. The committed acceptance record
`halley2011/reference/halley1d_stage0_baseline.json` freezes the executable,
input, run-log, C2 output, and Rubin-comparison checksums plus the final C2
occupancy and comparison values. The test re-analyzes that archive when it is
available. Exact volume-weighted L1 diagnostics require an AMPS cell-measure
column or an explicitly supplied aligned table; the tool marks legacy output
as unavailable instead of estimating cell volume from plot coordinates.

Validate an instrumented production archive as part of the same runner with:

```sh
HARTLEY2_STAGE0_PRODUCTION_ROOT=/absolute/path/to/production \
  srcHartley2/tests/run_stage0_tests.sh
```

The accepted closure run contains 20 reports, 2,984 exact cell measures, an
observed represented C2 carrier, and byte-identical model output relative to
the diagnostics-off baseline. The controlled four-rank enabled/disabled smoke
pair measured +0.93% wall time and +340 KiB maximum resident memory. These are
diagnostic-implementation results; the 10,000-iteration run remained
nonconverged and is not a science-quality Halley solution.

## Stage 1 species population policies

Population control can now be selected globally, by a derived comet-physics
group, or by one explicit species. The canonical Halley input uses:

```text
RadialParticleWeightAdaptation=on
RadialParticleWeightAdaptationInterval=10000
RadialParticleNumberPerCellRange=800 | 3200
RadialParticleNumberPerCellRange(CHEMICAL_ONLY)=50 | 200
RadialParticleWeightLock=DIRECT_SOURCE
```

`DIRECT_SOURCE` contains species whose configured surface or extended source
rate is positive. `CHEMICAL_ONLY` contains species transitively reachable from
those sources through enabled, positive-rate reactions but having no direct
source of their own. The input generator derives both groups after reading the
source and reaction databases; no comet-specific species list is compiled
into the controller. An explicit species name can replace a group selector.
Resolution is deterministic: explicit species overrides a derived group,
which overrides the unqualified global record.

A selector may have one range policy or one lock policy, not both at the same
precedence level. Ranges contain positive integers and require
`maximum >= 2*minimum`, matching the power-of-two controller. Duplicate
records, unknown selectors, invalid bounds, and range/lock contradictions are
fatal input errors. A locked species keeps its source-derived weight during
regular population adaptation and chemistry allocation protection. The
collective emergency common-factor coarsener may still change it as a
last-resort buffer-safety operation; that override is always logged.

`RadialParticleWeightAdaptationInterval=N` decouples numerical adaptation
from the science sampling length. Adaptation occurs only after a complete
mover/chemistry iteration and never while a continuation particle is in
flight. Omitting the interval preserves legacy behavior: adapt only at a full
sampling break. Startup prints one resolved row per species with its class,
winning selector, lock state, range, physical source rate, and actual surface
carrier `W/dt`.

Focused parser tests are always available:

```sh
srcHartley2/tests/run_stage1_tests.sh
```

The complete 2026-11-05 evidence, when the external validation archive is
present, is reproduced with:

```sh
HARTLEY2_STAGE1_POLICY_LOG=validation_runs/stage1_population_policy_20261105/canonical_policy/run.log \
HARTLEY2_STAGE1_ALTERNATE_LOG=validation_runs/stage1_population_policy_20261105/alternate_policy/run.log \
HARTLEY2_STAGE1_CADENCE_LOG=validation_runs/stage1_population_policy_20261105/cadence/run.log \
HARTLEY2_STAGE1_LEGACY_BASELINE=validation_runs/stage0_closure_20261105/diagnostics_off/production/PT/plots \
HARTLEY2_STAGE1_LEGACY_CANDIDATE=validation_runs/stage1_population_policy_20261105/legacy_candidate/PT/plots \
  srcHartley2/tests/run_stage1_tests.sh
```

The accepted result is 27 resolved species, unchanged locked direct-source
rates/carriers under the alternate daughter target, cadence records at
iterations 0/50/100, and 54 byte-identical legacy output files. These tests
close Stage 1 numerical-policy behavior; they do not change the existing
nonconverged status of the long Rubin science comparison.

## Stage 3 systematic chemistry residuals

The optional input records are:

```text
CometRandomSeed=104729
ChemistryResidualSampling=systematic
CometDiagnostics=on
```

`SYSTEMATIC` maintains a randomized cumulative residual for each owner block,
reaction, and concrete product. It reduces cumulative integer-allocation
error from independent binomial noise to strictly less than one model
particle per fixed-carrier epoch without enabling individual statistical
weights. A reaction is still one transaction: accumulator phases advance only
after the hard product cap and particle-buffer preflight accept all products.
Regular carrier changes start a recorded epoch; emergency common scaling
preserves phases because it does not change any relevant carrier ratio.

The source default is `BERNOULLI`. Inputs that omit both new records reserve
no residual block data and do not call `rnd_seed`; this preserves the radial
application's historical inherited core RNG state. Systematic mode requires
an explicit seed. `--random-seed N` overrides it for independent pilot runs.
The input parser rejects duplicate/invalid seeds and modes, systematic mode
without a seed, and a residual mode with chemistry disabled.

Residual state is explicitly aligned in owner-block application storage and
contains a 53-bit ABI/network-layout signature. The signature prevents a
same-size reordered chemistry layout from being silently read as valid state.
Local-weight restart saves empty and populated owner blocks, particle records,
and the entire application region, then publishes restored owner data to
ghost blocks before particle movement.

Run the pure tests with:

```sh
srcHartley2/tests/run_stage3_accumulator_tests.sh
PYTHONDONTWRITEBYTECODE=1 python3 srcHartley2/tests/test_stage3_adaptive_bins.py
PYTHONDONTWRITEBYTECODE=1 python3 srcHartley2/tests/test_stage3_diagnostic_semantics.py
srcHartley2/tests/run_stage3_parser_tests.sh
```

Validate archived 1/4/12-rank migration/restart logs with:

```sh
HARTLEY2_STAGE3_LOG_1=validation_runs/stage3_systematic_20261105/mpi1/stage3.mpi1.log \
HARTLEY2_STAGE3_LOG_4=validation_runs/stage3_systematic_20261105/mpi4/stage3.mpi4.log \
HARTLEY2_STAGE3_LOG_12=validation_runs/stage3_systematic_20261105/mpi12/stage3.mpi12.log \
  srcHartley2/tests/run_stage3_integrated_tests.sh
```

The focused executable mode is
`--validate-systematic-residual [--random-seed N]`. It forces a real AMR
redistribution and performs a local-weight restart round trip with one
controlled particle per MPI rank. The test record is strict: owner state and
particle checksums must survive, restored state must reach boundary ghosts,
and no partial transaction mutation is allowed. It also invokes the exact
center-cell writer and uncertainty report; the PASS record requires both.

The focused record calls its local check `plan-retry-checks`; it does not
claim to force `ProcessParticle()` through a full-buffer rejection. The
Stage-2 1/4/12 integrated transaction archive supplies the real capacity
rejection and complete particle/list/state rollback prerequisite.

For additional, non-destructive statistical output:

```sh
python3 srcHartley2/visualization/adaptive_radial_bins.py \
  seed1/PT/diagnostics/stage3_cells.C2.out=3.csv \
  --replicate seed2/PT/diagnostics/stage3_cells.C2.out=3.csv \
  --replicate seed3/PT/diagnostics/stage3_cells.C2.out=3.csv \
  --replicate seed4/PT/diagnostics/stage3_cells.C2.out=3.csv \
  --cell-measures seed1/PT/diagnostics/cell_measures.csv \
  --sampling-length 10000 --minimum-exposure 100000 \
  --output seed1/PT/plots/C2.stage3-adaptive
```

The tool selects common bin boundaries from model-particle exposure and then
estimates uncertainty across independent seeds. It verifies embedded
species/output identity, exact geometry, sampling length, explicit distinct
seeds, residual mode, and MPI metadata before combining replicas. Fewer than
four inputs are retained but marked not decision-ready. The scientific Stage-4
decision policy is now frozen in
`validation/stage3_science_policy.json`; the required
Bernoulli/Systematic matrix has not yet been run.

The production decision workflow is documented in
`validation/README.md`. In short, one cleanly built executable supplies both
runtime residual modes, four fixed seeds are run in both modes at 12 MPI
ranks, the first 10,000 iterations are burn-in, and three subsequent complete
10,000-iteration windows are evaluated through 40,000 iterations. The
evaluator uses all-eight-run minimum exposure on exact fully contained cells
from 10--100,000 km altitude and checks CV, zero cells, 50--200 population
coverage, mode-mean agreement, chemistry/source closure, emergencies,
convergence, and an explicitly accepted Rubin regression baseline.

Preflight without launching the long matrix:

```sh
python3 srcHartley2/validation/run_stage3_science_gate.py plan \
  validation_runs/stage3_science_gate_YYYYMMDD \
  --executable amps --input input/halley1d.input \
  --pic-global build/pic/picGlobal.dfn
srcHartley2/tests/run_stage3_science_gate_tests.sh
```

The operational runner supports `run`, `resume`, and `status`, always
re-hashes completed evidence, and never resumes a partial physical state. The
read-only evaluator returns `PASS_SKIP_STAGE4`, `FAIL_CONSIDER_STAGE4`,
`INCONCLUSIVE_RERUN`, or `INVALID_EVIDENCE`. Operational completion alone is
not a scientific pass.

The generic chemistry executor currently depends on the radial
batch/capacity/emergency policy. Chemistry is therefore rejected for the
2-D/3-D backend until a geometry-neutral policy with explicit input values is
implemented; this does not remove generic 3-D geometry/source support.

## Stage 4 individual-correction compatibility

Stage 4 is complete for the MPI DIM1 radial configuration using
species-local timesteps, species-local block weights, local-timestep face
resampling, and the general NTC collision model. It is an opt-in compatibility
foundation, not a change to the production Halley setup. The canonical
`input/halley1d.input` remains correction-off with Systematic chemistry
residuals; the repository-root `halley1d.input` is a symbolic link to that
file and must not be edited as though it were a disposable generated copy.

The correction-aware mover uses fixed-`c` residual transfer. At a face,

```text
mu = (W_old*c/dt_old)/(W_new*c/dt_new)
   = (W_old/dt_old)/(W_new/dt_new).
```

Every survivor and clone keeps `c`; diagnostics and physical moments use
`W*c/dt`. Block-weight adaptation likewise changes only `W`. General NTC
selects particles in proportion to `W*c`, uses the carrier ratio in each
post-collision update probability, and rejects a same-particle draw instead
of conditionally redrawing it. Simplified NTC explicitly terminates when
individual corrections are compiled on. The radial setup also rejects MF or
another unvalidated collision model before mesh construction.

The focused modes are:

```text
--validate-local-timestep [--validation-unit-corrections]
--validate-weight-adaptation
--validate-chemistry
--validate-unsupported-simplified-ntc
```

The integrated runner executes the wide-correction and all-`c=1` mover,
adaptation/core, and chemistry cases at 1, 4, and 12 MPI ranks. It also checks
the expected simplified-NTC failure:

```sh
srcHartley2/tests/run_stage4_integrated_tests.sh \
  --executable /absolute/path/to/correction-on/amps \
  --output-root /new/absolute/path/to/stage4-archive \
  --validation-particles 200000 \
  --chemistry-validation-events 200000
```

The completed archive is
`validation_runs/stage4-correction-on-full-20261105`. Its frozen executable,
input, generated `picGlobal.dfn`, and `SHA256SUMS` are under `provenance/`.
Rebuild the same configuration without modifying the canonical input by using
the archived correction-on input as a command-line make override. As required
for this application, remove exactly the root build directory first:

```sh
cd /data/vtenishe/COMET1D/AMPS
rm -rf /data/vtenishe/COMET1D/AMPS/build
make -j8 EXTRACOMPILEROPTIONS=-O3 \
  InputFileAMPS=validation_runs/stage4-correction-on-full-20261105/provenance/halley1d.correction-on.input
rg 'INDIVIDUAL_PARTICLE_WEIGHT_MODE|PARTICLE_COLLISION_MODEL_ ' \
  build/pic/picGlobal.dfn
```

The combined identity audit compares the all-`c=1` correction-on arm with the
independently built correction-off reference. It requires the same path
coverage and requires both arms to pass the same physical-flux tests:

```sh
stage4_on=validation_runs/stage4-correction-on-full-20261105
stage4_off=validation_runs/stage4-correction-off-reference-20261105
stage4_args=
for stage4_rank in 1 4 12; do
  stage4_args="$stage4_args \
    --local-log $stage4_rank $stage4_on/mpi$stage4_rank/local-timestep/run.log \
    --unit-local-log $stage4_rank $stage4_on/mpi$stage4_rank/local-timestep-unit/run.log \
    --off-unit-local-log $stage4_rank $stage4_off/mpi$stage4_rank/run.log \
    --adaptation-log $stage4_rank $stage4_on/mpi$stage4_rank/weight-adaptation/run.log \
    --chemistry-log $stage4_rank $stage4_on/mpi$stage4_rank/chemistry/run.log"
done
# shellcheck disable=SC2086 # the variable deliberately contains CLI words
PYTHONDONTWRITEBYTECODE=1 python3 \
  srcHartley2/tests/test_stage4_integrated_logs.py $stage4_args
```

That audit reports both `PASS Stage-4 all-c=1 correction identity` and
`PASS Stage-4 correction-on integrated matrix`. Correction-off Stage-2 and
Stage-3 regressions are retained in
`validation_runs/stage4-correction-off-regression-20261105`. To restore the
normal development executable after any correction-on rebuild, again remove
only the exact `build` directory and build with
`InputFileAMPS=input/halley1d.input`; confirm that the generated macro is
`_INDIVIDUAL_PARTICLE_WEIGHT_OFF_`.

The Stage-4 archive does not validate hybrid diagnostic counters, CUDA
continuation, field-line/periodic movers, or correction-on 2-D/3-D chemistry.
It predates the separately tested Stage-5 windows documented below and still
does not provide importance sampling. Do not merge to `master` from this
workflow; the dependency order and conflict invariants are recorded in
`STATUS.md`.

## Complete Stage 0--4 test runner

`tests/run_all_stage_tests.sh` is the single entry point for every focused
test and archived runtime-evidence checker used to close Stages 0 through 4.
Its default archived mode is non-destructive: it does not compile AMPS, edit
an input, remove `build`, launch a production simulation, or overwrite an
archive. It does run the Stage-0 standalone diagnostic at 1, 4, and 12 MPI
ranks and compiles the temporary standalone Stage-0/2/3 unit executables.

The runner makes all optional archive inputs of the individual stage runners
mandatory. Missing evidence fails preflight; it cannot turn into a successful
all-stage result containing `SKIP` records. It also verifies that the current
generated build is correction-off before using it as the read-only Stage-3
parser fixture. The runner and that parser test use standard `grep`; `rg` is
not a runtime prerequisite.

Inspect the resolved tools, evidence paths, configuration, and coverage
without executing tests:

```sh
srcHartley2/tests/run_all_stage_tests.sh --plan
```

Run the complete archived suite:

```sh
srcHartley2/tests/run_all_stage_tests.sh
```

The validated 2026-11-05 archive bundle completes 17 independently logged
test jobs organized into five stage groups and ends with:

```text
Hartley2 all-stage test summary status=PASS stage4-mode=archived
  tests: run=17 passed=17 skipped=0 failed=0
  test-groups: run=5 passed=5 skipped=0 failed=0
PASS Hartley2 all-stage suite tests=17 test-groups=5 stage4-mode=archived
```

Here a **test** is one independently logged orchestrated job and a **test
group** is one roadmap stage (`stage0` through `stage4`). The summary reports
each stage's contained pass/skip/fail counts. The runner continues after an
independent job fails so the final result lists every observed failure rather
than only the first one.

Every job writes a separate log while the suite is running. Successful default
runs remove those temporary logs. If any test fails or emits `SKIP`, the logs
are retained automatically and the final `failed-tests`/`skipped-tests` lists
give the test name, stage, exit code, and absolute log path. The summary itself
is retained as `summary.txt` in the reported log root. To retain logs even for
a successful run, select a new directory explicitly:

```sh
srcHartley2/tests/run_all_stage_tests.sh \
  --log-root validation_runs/all-stage-logs-YYYYMMDD
```

An existing log root is rejected so evidence from separate attempts cannot be
mixed. Any failed or skipped test makes the complete runner exit nonzero.

The groups include Stage-0 MPI/analyzer/baseline/production checks; all
Stage-1 parser, runtime, cadence, and legacy checks; Stage-2 pure and 1/4/12
integrated checks; every Stage-3 pure/parser/visualization/semantics/science-
gate/integrated check; the final correction-off Stage-2/3 regressions; frozen
Stage-4 executable provenance; MF and simplified-NTC rejection; and the
correction-on wide/all-`c=1` 1/4/12 matrix.

To repeat the Stage-4 MPI matrix rather than replay its archived compact
records, first provide an independently clean-built correction-on/general-NTC
executable. Compilation remains outside the runner: remove exactly the root
`build` directory before that build. Then use a new output path:

```sh
srcHartley2/tests/run_all_stage_tests.sh \
  --stage4-mode fresh \
  --stage4-executable /absolute/path/to/amps.correction-on \
  --stage4-output-root /new/absolute/path/to/stage4-validation
```

Fresh mode invokes `run_stage4_integrated_tests.sh`, then applies the combined
all-`c=1` identity checker against the frozen correction-off reference. It
never installs the supplied executable as the repository's `amps`, so a
failure cannot leave the workspace correction-on. The mover particle and
chemistry-event counts default to 200,000 and can be changed with the explicit
CLI options printed by `--help`. MPI rank counts remain the validated 1, 4,
and 12 and therefore remain below the approved 18-process limit.

`--evidence-root` may point to a relocated archive bundle, but the complete
dated subdirectory structure printed by `--help` is required. This is
intentional: automatically choosing a newest-looking run could combine
different executables or incomplete evidence and produce a false pass.

## Stage 5 correction windows

Stage 5 is implemented for the same correction-on DIM1 radial configuration
qualified in Stage 4. Exact chemistry still creates one carrier for every
concrete stoichiometric product with

```text
c_daughter = (W_parent*c_parent/dt_parent)/(W_daughter/dt_daughter).
```

The new controller bounds this dimensionless correction after each completed
time step. The configured interval is inclusive. A correction below its lower
edge is rouletted to the target with probability `c/target`; a correction
above its upper edge is split into `ceil(c/target)` complete-record copies,
subject to a per-particle hard maximum. Split corrections reconstruct the
input correction exactly in production summation order. The controller runs
after any coincident block-base-weight adaptation and refuses an unfinished
local-timestep continuation.

The transaction has a read-only preflight and an all-rank decision. It counts
all required clones, retains `RadialAllocationReserveFraction`, and does not
credit possible roulette deletions. Rejection therefore precedes RNG draws,
particle-list edits, tracker changes, or diagnostics on every rank. Accepted
splits preserve position, velocity, rotational/vibrational state, chemistry
stage, and all other particle fields; each new tracked copy obtains a new
tracker identity. Only authoritative owner-block particle lists are visited.

### Input records

Correction-window records are legal only when both the core correction mode
and the application window switch are on:

```text
ParticleWeightCorrectionMode=on
ParticleWeightCorrectionWindowMode=on
ParticleWeightCorrectionWindow=LOWER | TARGET | UPPER
ParticleWeightCorrectionWindow(DIRECT_SOURCE)=LOWER | TARGET | UPPER
ParticleWeightCorrectionWindow(CHEMICAL_ONLY)=LOWER | TARGET | UPPER
ParticleWeightCorrectionWindow(H2O)=LOWER | TARGET | UPPER
ParticleWeightCorrectionWindowMaxSplits=32
```

The unqualified global record and the split limit are required. An explicit
species record takes precedence over `DIRECT_SOURCE`/`CHEMICAL_ONLY`, which
takes precedence over global. Direct-source windows must contain `c=1`.
Malformed, duplicate, nonfinite, unordered, contradictory, or off-mode
records stop configuration. Omitting the mode and all records is the legacy
default-off behavior.

The acceptance build used `0.25|1|4` globally, `0.5|1|2` for direct sources,
`0.0625|0.25|4` for chemical-only species, and 32 maximum split particles.
Those are test values, not a calibrated Halley/Rubin policy. Canonical
`input/halley1d.input` remains correction-off and contains no window records.

### Focused and MPI tests

Run the parser and pure C++ planner tests without configuring AMPS:

```sh
srcHartley2/tests/run_stage5_correction_window_tests.sh
```

An initialized correction-on build exposes:

```text
--validate-correction-windows
```

The validator covers both inclusive boundaries, exact full-state splitting,
collective capacity rollback, restart, low-side roulette expectation and
physical moments, direct-source identity, empty lists, and owner-only MPI
traversal. The integrated runner pairs it with the exact chemistry ensemble at
1, 4, and 12 ranks:

```sh
srcHartley2/tests/run_stage5_integrated_tests.sh \
  --executable /absolute/path/to/correction-on/amps \
  --output-root /new/absolute/path/to/stage5-validation \
  --chemistry-validation-events 200000
```

The completed archive is
`validation_runs/stage5-correction-window-full-20261105`. Its strict checker
reported:

```text
PASS Stage-5 correction-window integrated matrix 1rank=pass 4rank=pass 12rank=pass
```

This closes the implementation stage, not the scientific C2 pilot. Window
selection, independent-seed uncertainty, clone decorrelation burn-in, parent
profile invariance, minor-species flux closure, and the converged Rubin
comparison remain required before enabling the feature in a production input.
Correction-on DIM2/DIM3 chemistry is still unsupported; the pure planner is
geometry-neutral, but the current executor intentionally remains radial.

## Stage 6 daughter-packet splitting

Stage 6 adds input-selected independent recoil packets for statistically weak
chemical daughters. For `K` concrete products and packet count `N`, one
selected physical reaction creates exactly `K*N` model particles. Every
packet contains every product from one recoil solution; separate packets use
separate recoil calls. Product `k` keeps its Stage-4 represented rate because
its unsplit correction is divided into `N` positive parts with an exact final
remainder.

Packet selection is resolved during preprocessing:

```text
ChemistryDaughterPacketMode=on
ChemistryDaughterPacketMaximum=8
ChemistryDaughterPacketCount(GLOBAL)=1
ChemistryDaughterPacketCount(CHEMICAL_ONLY)=2
ChemistryDaughterPacketCount(SPECIES:C2)=8
ChemistryDaughterPacketCount(REACTION:C2_2C)=8
```

Reaction ID has highest precedence. Otherwise, the maximum applicable count
among all concrete products is used; explicit species overrides
`CHEMICAL_ONLY`, which overrides global. Explicit sink selectors are illegal,
and the parser rejects any policy whose `K*N` exceeds the configured chemistry
transaction limit. With the mode and records omitted, all generated counts
are one and the historical one-recoil RNG path is retained.

Run the pure policy/planner tests with:

```sh
srcHartley2/tests/run_stage6_daughter_packet_tests.sh
```

The initialized-mesh validator is selected with
`--validate-daughter-packets`. The integrated runner exercises it together
with the Stage-5 correction-window and exact-chemistry validators at 1, 4,
and 12 ranks:

```sh
srcHartley2/tests/run_stage6_integrated_tests.sh \
  --executable /absolute/path/to/stage6/amps \
  --output-root /new/absolute/path/to/stage6-validation \
  --chemistry-validation-events 200000
```

The accepted archive is
`validation_runs/stage6-daughter-packet-full-20261105-r4`. It covers 46
neutral reactions and 23 sinks at `N=1,2,4,8`, plus a real repeated-product
`C2_2C` allocation/rollback transaction. The fixed test packet counts and
broad validation windows are not production Halley parameters. A transported
independent-seed variance/cost and flux-closure study remains necessary before
enabling packets in a science input.

## Stage 7 stochastic reaction-channel importance sampling

Stage 7 changes only the conditional channel choice after the ordinary total
reaction test succeeds.  A configured positive-rate, unit-unimolecular parent
family must list every sibling, including product-free ionization sinks, with
strictly positive target probabilities `q*` that sum to one:

```text
ChemistryChannelImportanceMode=on
ChemistryChannelBiasedProbability(REACTION:C2_2C)=0.5
ChemistryChannelBiasedProbability(REACTION:C2_ION)=0.5
```

Stable reaction IDs are case-sensitive.  Partial families, duplicates,
unknown/zero-rate/non-unit targets, disabled-mode records, and non-normalized
families stop preprocessing.  Unconfigured families retain the historical
`rnd()*totalRate` raw-rate subtraction and RNG count.  For a selected channel,
the product correction receives `L=q/q*` before Stage 6 divides it among `N`
packets, so each packet carries `L/N`.  Parent destruction is not weighted;
sinks create no daughter but remain present in diagnostics.

Run the pure input-policy and planner tests with:

```sh
srcHartley2/tests/run_stage7_channel_importance_tests.sh
```

The initialized validator uses `--validate-channel-importance` and writes a
separate Stage-7 CSV/JSON pair containing raw selections, `q`, `q*`, `L`,
likelihood sums, weighted channel rate, and parent-family ESS.  The integrated
runner replays Stages 4--6 and validates 1/4/12-rank deterministic identity:

```sh
srcHartley2/tests/run_stage7_integrated_tests.sh \
  --executable /absolute/path/to/stage7/amps \
  --output-root /new/absolute/path/to/stage7-validation \
  --chemistry-validation-events 200000
```

The accepted archive is
`validation_runs/stage7-channel-importance-full-20261105-r2`.  Its validation
policy uses equal target probabilities for the complete C2 and H2O families;
those values test sinks, rare branches, variance changes, and Stage-6 `L/N`
composition and are not production Halley bias parameters.  The accepted
matrix covers 27 parent families and all 69 reactions at 200,000 events per
family.  Deterministic all-channel branching remains deliberately deferred.

## Complete Stage 0--7 test runner

The general runner now has eight groups and 25 tests. Archived mode checks the
Stage-7 pure policy/planner job and frozen importance matrix, including its
Stage-4/5/6 dependency replay and CSV/JSON provenance, in addition to every
Stage 0--6 job:

```sh
srcHartley2/tests/run_all_stage_tests.sh --plan
srcHartley2/tests/run_all_stage_tests.sh
```

Fresh Stage-5 mode requires an independently clean-built, correction-enabled
binary and a new output directory:

```sh
srcHartley2/tests/run_all_stage_tests.sh \
  --stage5-mode fresh \
  --stage5-executable /absolute/path/to/correction-on/amps \
  --stage5-output-root /new/absolute/path/to/stage5-validation
```

Fresh Stage-6 mode likewise requires a frozen packet-enabled binary and a new
immutable output root:

```sh
srcHartley2/tests/run_all_stage_tests.sh \
  --stage6-mode fresh \
  --stage6-executable /absolute/path/to/stage6/amps \
  --stage6-output-root /new/absolute/path/to/stage6-validation
```

Fresh Stage-7 mode requires a frozen importance-enabled binary and a new
immutable output root:

```sh
srcHartley2/tests/run_all_stage_tests.sh \
  --stage7-mode fresh \
  --stage7-executable /absolute/path/to/stage7/amps \
  --stage7-output-root /new/absolute/path/to/stage7-validation
```

The current build used by archived mode must remain correction-off because it
is also the read-only Stage-3 parser fixture. The Stage-5/6/7 logs are replayed
from frozen archives; the runner never mutates input or `build/`. As before,
it counts every test and group, retains separate logs for failures/skips, and
prints every failed test name and absolute log path.

The final canonical-off replay completed with:

```text
Hartley2 all-stage test summary status=PASS stage4-mode=archived stage5-mode=archived stage6-mode=archived stage7-mode=archived
  tests: run=25 passed=25 skipped=0 failed=0
  test-groups: run=8 passed=8 skipped=0 failed=0
PASS Hartley2 all-stage suite tests=25 test-groups=8 stage4-mode=archived stage5-mode=archived stage6-mode=archived stage7-mode=archived
```

## Stage 8 deterministic total-reaction forcing

Stage 8 addresses rare total reaction occurrence without changing a physical
reaction rate or introducing a tunable occurrence probability. For each
configured unit-unimolecular parent and one logical local-timestep segment,
the application retains both exact physical outcomes:

```text
P = -expm1(-lambda*dt)
survivor carrier = C*(1-P)
reacted carrier  = C*P
```

The reacted branch remains a parent particle while it moves to a reaction
time sampled from the physical exponential distribution conditional on an
event inside the segment. It then converts to a product, or is removed for a
sink, and completes its unspent normalized interval using the product local
timestep. Stage-7 channel likelihood `L` and Stage-6 packet count `N` compose
as `C*P*L/N` for each concrete product packet. Products from one packet use a
common deterministic conservative recoil solution.

Enable it only with an explicit per-species policy:

```text
ChemistryTotalReactionForcingMode=on
ChemistryTotalReactionForcing(SPECIES:C2)=on
```

Species symbols are case-sensitive and the value must be `ON`; there is no
numeric `P*` record. Enabled forcing requires DIM1 radial chemistry,
individual particle corrections, Bernoulli residual sampling, Stage-5
windows, an explicit `RadialMaxPhotochemicalProductsPerReaction`, and
diagnostics. The parser rejects incomplete prerequisites, duplicates,
unknown species, and species without a positive-rate unit-unimolecular
family. Omitting the mode and all records preserves the default-off path.

The implementation stores a 16-byte pending-event record only in enabled
builds. Before any mutation or Stage-8 RNG draw it preflights the complete
`K*N` daughter transaction and retained survivor; an insufficient buffer
leaves bytes, path, and RNG state unchanged. The generic core event callbacks
are all null by default and must be installed as a complete set. A geometric
face at exactly the event time wins, so ordinary block/MPI continuation is
completed before conversion.

Current qualification is deliberately narrow: the production hook is the
first-order force-free traced mover in the DIM1 spherical radial application,
with a constant total reaction rate during a logical segment. The pure
planner tests piecewise-constant hazards, but production multi-block hazard
integration, second-order/field-line/periodic/CUDA movers, and DIM2/DIM3 are
not yet supported claims.

Run the pure policy/planner checks without changing a build:

```sh
srcHartley2/tests/run_stage8_total_reaction_forcing_tests.sh
```

An enabled executable exposes `--validate-total-reaction-forcing`. Its mover
fixture uses the production execution track: sixteen C2 products react after
a local/MPI face and one product-free sink is due exactly at the face. The
test requires face-first ordering, advance-before-copy, destination-owner
processing, core-owned deletion, complete remaining flight, and exact carrier
closure. The full runner replays the Stage-4 mover/adaptation matrix and
Stages 5--8 at 1, 4, and 12 ranks:

```sh
srcHartley2/tests/run_stage8_integrated_tests.sh \
  --executable /absolute/path/to/stage8/amps \
  --output-root /new/absolute/path/to/stage8-validation \
  --validation-particles 200000 \
  --chemistry-validation-events 200000
```

The output root must not exist. The accepted immutable evidence is
`validation_runs/stage8-total-reaction-forcing-full-20261105-r4`; the `-r1`
and `-r2` siblings are failed parser-development evidence. The `-r3` sibling
passed direct callback checks but did not exercise events through the real
mover or replay Stage 4, so it is pre-acceptance evidence. None may be
substituted for `-r4`.
Validation fixtures such as broad windows, `N=8`, and equal C2/H2O `q*` are
not science parameters.

## Complete Stage 0--8 test runner

The general runner now has nine groups and 30 independently logged tests. In
addition to every Stage 0--7 check it runs the Stage-8 policy/planner and
analytic-oracle job, verifies frozen provenance, replays every lower-stage
compatibility record from the Stage-8 archive, and strictly checks the forced
segment records at 1, 4, and 12 ranks:

```sh
srcHartley2/tests/run_all_stage_tests.sh --plan
srcHartley2/tests/run_all_stage_tests.sh
```

Fresh Stage-8 mode requires a separately clean-built enabled binary and a new
immutable directory:

```sh
srcHartley2/tests/run_all_stage_tests.sh \
  --stage8-mode fresh \
  --stage8-executable /absolute/path/to/stage8/amps \
  --stage8-output-root /new/absolute/path/to/stage8-validation
```

The archived runner remains read-only and the current repository build must
be canonical correction-off for the Stage-3 parser fixture. Successful
closure reports:

```text
tests: run=30 passed=30 skipped=0 failed=0
test-groups: run=9 passed=9 skipped=0 failed=0
PASS Hartley2 all-stage suite tests=30 test-groups=9 stage4-mode=archived stage5-mode=archived stage6-mode=archived stage7-mode=archived stage8-mode=archived
```

This is implementation compatibility evidence, not the Stage-9 transported
science result. Production use still requires independent-seed analog/forced
variance and cost comparisons, parent-profile invariance, collision feedback,
minor-species flux closure, buffer headroom, convergence, and regenerated
Rubin comparison products.
