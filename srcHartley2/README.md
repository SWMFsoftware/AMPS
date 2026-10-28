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
| `RadialMesh*` | DIM1 spherical domain and logarithmically expanding target resolution |
| `RadialSimulation.*` | DIM1 source injection, timestep selection, inner re-emission, outer loss, and time stepping |
| `halley2011/derive_halley1d_controls.py` | Reproducible Halley surface-resolution, injection-count, and daughter-speed derivation |
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

The core model-particle collision adapter must also select the user-defined
cross-section route. That route passes the current block and center cell to
`CometPhysics`, allowing the H2O-H2O Crifo law to use temperature. A generated
Halley build should contain the final override
`_PIC__PARTICLE_COLLISION_MODEL__CROSS_SECTION__USER_DEFINED_` in
`build/pic/picGlobal.dfn`.

Rubin Table 3 lists all 66 unordered pairs among the 11 parents but says that
the table is only a subset; it does not print the daughter-pair values used by
the paper model. The application does not assign zero or invent a mixing rule
for those missing pairs.

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
event, creates all products from that event together, conserves event momentum
and excess kinetic energy, and applies exact individual statistical-weight
corrections. Chemistry therefore requires
`ParticleWeightCorrectionMode=on`.

For a wide logarithmic domain, `SpeciesLocalParticleWeight` scales each
species' base weight with its local timestep. A steady radial flux then uses
approximately the same number of model particles per cell instead of growing
with cell width. Chemistry retains exact event weights through the individual
correction stored on every particle.

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
Courant number 0.5, 400 m/s parent baseline, disabled methanol branch, and
controls derived by `halley2011/derive_halley1d_controls.py`. The 500-step
sampling window and 10,000-iteration hard ceiling are explicit. The selected
0.36 eV OH branch uses the lower published rate, `3.0e-6 s^-1` at 1 AU. The
run compares the volume-weighted relative L1 change of every species' complete
radial density profile and terminates after all species remain within 5% for
three consecutive window-to-window comparisons. The first complete window is
the reference, so convergence needs at least four windows. Multi-rank tests
have crossed output and block-migration boundaries, verifying runtime plumbing.
A later audit found that those test binaries used AMPS's zero placeholder
cross section because the Rubin callback was not connected. That dispatch is
now fixed, but the repaired diagnostic stops at iteration 1 because no sampled
cell temperature exists yet for the H2O law. No output currently qualifies as
a converged physical Halley solution. Exact evidence and unresolved physics
are recorded in `STATUS.md`.

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

This reproduces the neutral-species groupings and units of Rubin et al. (2011)
Figs. 1--17 and the modeled H2O radial velocity corresponding to Fig. 18. The
paper's Giotto NMS mass-channel curves require an instrument-response model
that is not part of this application and are not included.
