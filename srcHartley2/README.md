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
stoichiometric product. With `SpeciesGlobalParticleWeight` and correction
disabled, it converts the parent particle rate `W_parent/dt_parent` to each
product rate `W_product/dt_product` with unbiased residual split/roulette.
Number, mass, momentum, and energy production rates are then conserved
statistically, while event-by-event product correlations acquire Monte Carlo
variance.

### Local-timestep face transport

`LocalTimeStepFluxResampling=on` enables the core implementation used by the
Halley case. It is deliberately limited to `SpeciesLocalTimeStep`, a particle
weight that is spatially constant for each species, correction mode off, CPU
execution, a spatial non-periodic AMR mesh, and the trace mover. Unsupported
CUDA, field-line, and periodic configurations fail at compile time. At every
AMR block face, the trace mover stops at the face instead of first calculating
a final position and subsequently looking up its block.

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
kinetic-energy fluxes in expectation without changing velocity. In the
configured species-global mode `W_old=W_new`, so `mu=dt_new/dt_old`.

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
physical intervals during one outer iteration. A spatially global weight also
means the model population needed for a steady flux grows in proportion to
the local timestep. Across the Halley logarithmic domain that ratio is very
large, so a long convergence run must monitor the particle buffer. A future
production improvement should combine block-local weights with a conservative
face reservoir or weighted-particle representation, or introduce asynchronous
block clocks; simply capping the clone count would bias the flux.

The application uses one block-independent weight per species rather than one
weight shared by every species. Direct-parent weights retain the requested
surface statistics. Species with no direct source use the largest parent
weight so a reaction cannot create an enormous instantaneous clone burst;
chemistry performs
`(W_parent/dt_parent)/(W_product/dt_product)` split/roulette explicitly.

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

`--validation-mover` accepts `first`, `second`, or `both`. The single-face
matrix tests timestep ratios 0.25, 0.5, 1, 2, and 2.5 in both radial directions
on a same-owner face; ratios 0.5, 1, and 2.5 are repeated across an MPI
ownership face when available. Four velocity groups include radial and
transverse components. Every case checks particle-list membership, final
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

The completed 18-rank validation does **not** authorize merging the full
branch into `master`: conservative face and chemistry tests pass, but
species-global-weight population stability and converged Rubin comparison
fail. `STATUS.md` records the gate-by-gate disposition, the required
block-local-weight/reservoir follow-up, independently reviewable subsets, file
dependency order, conflict checks, clean-commit requirements, and the
integration-branch procedure. No merge should begin until the population fix
has been validated on a clean, identified commit and the Halley run satisfies
three consecutive 5% convergence windows without buffer exhaustion.

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
