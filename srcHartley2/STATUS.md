# Generic comet / Hartley 2 status

Date: 2026-10-27
Branch: `1d_comet_dev`

## Current result

| Check | Result |
| --- | --- |
| Generic model input preprocessing | Pass for DIM1 and DIM3 |
| Existing DIM3 clean C++17/MPI compile and link | Pass with legacy warnings |
| Spherical DIM1 clean C++17/MPI compile and link | Pass |
| DIM1 radial mesh runtime | Pass on one MPI process |
| DIM1 source/boundary/time-step backend compile and link | Pass |
| Reaction database C++ unit test | Pass |
| Reaction input/generator Perl test | Pass |
| Collision database C++ unit test | Pass |
| Rubin chemistry/collision generator Perl test | Pass; 66 parent pairs and Table 4 channels |
| Shape-property extraction | Pass; closed 32,040-triangle surface |
| Full DIM3 runtime | Not run; required Hartley 2 mission SPICE kernels are unavailable locally |
| Hartley 2 DIM1 particle simulation | Backend present; input remains mesh-only pending explicit science-run controls |
| Halley DIM1 one-step plumbing test | Pass with documented temporary numerical choices |
| Halley 2,000-cell/local-weight clean build and one-step run | Pass; 373 blocks and 2,984 cells |
| Halley DIM1 four-rank runtime | Pass through one complete 500-step output window and the following block migration |
| Halley DIM1 convergence run | Blocked before a physical result; the first attempted run used an unwired zero cross section, and the repaired build stops at iteration 1 because no initial H2O cell temperature is available |
| Core collision cross-section dispatch | Fixed and clean-built; NTC/MF now route the selected user law with block/cell state instead of always using the zero HS placeholder |
| Rubin collision-data coverage | Incomplete by source: Table 3 supplies 66 parent-parent pairs but omits daughter pairs |
| DIM1 result visualization | Pass for short aliases, single-/multi-file PNG/EPS output, and viewer control |
| Rubin paper-comparison artifact generator | Pass on completed transient output; no converged physical Halley output exists yet |
| DIM2 axisymmetric application | Interface reserved; execution backend not implemented |

The `build` directory is removed before each compilation. Hartley 2 and Halley
case data are kept in different input files; `srcHartley2` is the shared
generic application despite its historical directory name.
The repository build selectors now name `halley1d.input`. An earlier clean
DIM1 executable used the lower OH endpoint and one iteration to smoke-test the
execution path; the final requested controls are now in the tracked input.
The Halley runtime output path is `PT/plots`, relative to the launch directory;
the intended launch directory is `/data/vtenishe/COMET1D/AMPS`. A compile does
not create this runtime directory.

### Collision-dispatch audit and invalidated run

The first optimized 18-rank convergence attempt reached iteration 2,533 and
wrote five complete 500-step windows, but it is **not a physical Halley
result**. Source inspection after the run showed that every NTC cross-section
evaluation still called
`PIC::MolecularData::MolecularModels::GetTotalCrossSection`. Its generated
initializer selected `HS::GetTotalCrossSection`, whose unchecked fallback is
a `1 x 1` zero table. The application-level Rubin callback and all 66
generated parent-pair entries were linked but unreachable. The run was stopped
and preserved as `PT.invalid-zero-cross-section` with
`halley1d.invalid-zero-cross-section.log`; those files must not be used for a
paper comparison.

The core collision algorithms now call one adapter selected by
`_PIC__PARTICLE_COLLISION_MODEL__CROSS_SECTION_`. The default mode retains the
existing HS/VHS/VSS function pointer. The user-defined mode passes velocities,
species, block, and center cell to the application callback, so a collision
law may depend on cell state without mutable global variables. All ten
cross-section evaluations in the simplified/full NTC and MF implementations
use the adapter. A clean `-O3` Halley build confirmed that the generated input
selects the user-defined mode and links the Rubin callback.

A four-rank diagnostic of that repaired binary reached the collision callback
at iteration 1 and then stopped with `non-positive translational temperature
for comet collision species 13`. The H2O-H2O Crifo law requires local gas
temperature, but the AMPS center-cell sampling accumulator has no completed
sample when the first collision step runs. No fallback temperature, temporal
lag, or relative-speed conversion has been added because Rubin et al. do not
specify that numerical initialization rule. This is a separate blocker from
the missing daughter collision matrix.

The first four-rank DIM1 runtime test exposed a core MPI neighbor-mask defect
in `src/pic/pic_block_send_mask.cpp`: the exchange code scanned the six faces
and four subface slots of a three-dimensional node even though a DIM1 AMR node
stores two faces and one slot per face. The resulting out-of-bounds neighbor
lookup caused `SIGSEGV` in `BlockElementSendMask::InitLayerBlock`. The core now
limits both the face count and subface indices to the compiled dimension; the
basic fallback mask likewise visits only the positive faces that exist.

After that core defect was removed, the four-rank test exposed a separate
application ownership error in local particle-weight initialization. Only one
MPI rank owns the nucleus-adjacent block after decomposition, but every rank
previously required a local copy. Initialization now verifies exactly one
global owner and uses `MPI_Allreduce` to distribute the owner's surface weight
and timestep before each rank scales the weights of its own radial blocks.

The first completed 500-step four-rank window wrote all species outputs, then
exposed another core DIM1 defect during sampling-break load redistribution.
The four symmetric block pack/unpack paths in `src/pic/pic_mesh.cpp` visited
positive x, y, and z faces unconditionally. They now serialize only the `DIM`
positive faces, keeping migrated-block size calculation and unpack ordering
valid for DIM1, DIM2, and DIM3.

The dimension and MPI corrections were tested in two stages on four MPI ranks. A
three-step diagnostic forced a sampling break on every step and exited
normally after an actual ownership change from `98/97/94/84` blocks to
`1/372/0/0`. The final binary then completed the configured 500-step sampling
window, wrote all 27 species files in `PT/plots`, and reported 87,485 model
particles. Its load balancer changed ownership from `98/97/94/84` blocks to
`1/2/5/365`; the run continued normally through iteration 529 and was then
stopped deliberately. Each output file has 3,011 lines and spans the full
radial mesh. This verifies initialization, particle transport, chemistry,
sampling output, and post-output block migration, but the later dispatch audit
showed that this test was collisionless. It therefore remains valid evidence
for MPI/runtime plumbing only. It does not demonstrate physical profile
convergence.

`visualization/plot_1d.py` provides a dimension-one result plotter that reads
the physical quantity names from each AMPS `VARIABLES` record. It always uses
a logarithmic radius axis, requires the Y scale on the command line, accepts
one or more output files, and labels every curve with its source filename.
Logarithmic Y plots explicitly mask and report non-positive samples. The tool
also lists available quantities and accepts documented short names including
`dens`, `temp`, `speed`, `vel`, velocity/temperature components, timestep,
weight, mesh metadata, and every `coll-SPECIES` column. Each output request
writes PNG and EPS files, then opens the PNG in the system desktop viewer
unless `--no-show` is selected or no graphical display is available. A
separately requested PDF or SVG is retained in addition to the mandatory PNG
and EPS outputs.

The visualization CLI accepts either a positive maximum cometocentric distance
with `--max-distance`/`--rmax` or a positive inclusive interval with
`--distance-range`/`--rrange`. Values use the selected radius unit and the two
forms are mutually exclusive. For logarithmic Y scaling, zero and negative
ordinates are replaced by gaps rather than passed to matplotlib or connected
across after row removal.

Multi-file logarithmic plots skip an entire file when no positive samples
remain after distance filtering, allowing populated species to be plotted when
an input glob also contains absent species. The script reports each omitted
file and returns an error only if every requested curve is empty.

## Merge inventory relative to `master`

This section is a file-level map for merging `1d_comet_dev` back into
`master`. The comparison base used while writing it is local `master` commit
`cb7be2b` (`Fix Hartley2 compatibility with current AMPS core`). It covers the
committed DIM1 mesh work and the later generic-comet, Halley, chemistry,
collision, runtime, and visualization work in the working tree. Recompute
`git diff --name-status master` immediately before the merge in case either
branch has moved.

No source file is intentionally deleted. Files marked **modified** must be
merged with their current `master` versions; files marked **new** must be
added to the branch before it is merged. The tables describe source-tree
files only. Generated build trees, run output, executable links, and parser
copies are listed separately below.

### Build and input files

| File | State | Why it changed and merge dependency |
| --- | --- | --- |
| `Makefile` | modified | Adds `-lmpi++` to `EXTRALINKEROPTIONS_F`. This is a toolchain/link setting, not comet physics. Its necessity on the target `master` toolchain must be checked with a clean link; keep it separate from the scientific changes if the master toolchain no longer needs it. |
| `input/hartley2.pl` | modified | Extends the historical Hartley parser into the generic comet preprocessor. It parses dimension/geometry, mesh units, source records, radial controls, chemistry, and neutral-collision data; selects DIM1 objects and spherical core options; validates required values; writes species-indexed arrays; and generates the chemistry and collision headers in the working source tree. This file depends on both new Perl modules and on the exact variable names in `CometModel.cpp`, `CometChemistry.cpp`, `CometPhysics.cpp`, `RadialMesh.cpp`, and `RadialSimulation.cpp`. Preserve the legacy Hartley keywords when resolving conflicts. |
| `input/hartley2.input` | modified | Adds the generic DIM3 case declarations, configurable NASTRAN mesh units, heliocentric distance, disabled chemistry, and the three neutral collision laws that were formerly hard-coded in `Comet.cpp`. The remaining historical Hartley 2 settings are intentionally retained. |
| `input/hartley2.radial.input` | new | Defines the mesh-only, gas-only, spherically symmetric Hartley 2/Shou 2025 DIM1 case. It contains the area-equivalent nucleus radius, 10,000 km altitude domain, 1,000-cell ideal logarithmic profile, and 1 m surface resolution estimate. It deliberately leaves `RadialSimulationMode=off`. |
| `input/halley1d.input` | new | Defines the gas-only Rubin 2011 Halley simulation: 27 species, 11 parent sources, 6 km inner radius, `10^6 km` outer radius, 2,000 ideal cells, derived surface resolution, local weights/timesteps, generated chemistry, collisions, and convergence controls. Dust, extended sources, and the optional methanol-to-formaldehyde channel are disabled. |
| `input/species.input` | modified | Adds the neutral species required by the Rubin network (`CH`, `CH2`, `N`, `CH3`, `NH`, `NH2`, `NH3`, `C2`, `C2H`, `C2H2`, `CN`, `HCN`, `C2H4`, `CO`, `C2H5`, `HCO`, `C2H6`, `H2CO`, and `CH3OH`). Merge by component name and avoid duplicate definitions if `master` has independently added any species. |
| `input/CometReactionNetwork.pm` | new | Parses arbitrary stoichiometric input reactions, resolves species names to compiled indices, validates rate-law fields, and emits `CometChemistry.generated.h`. It is the input-side producer for `ReactionNetwork.h` and `CometChemistry.cpp`. |
| `input/CometPhysicsDatabase.pm` | new | Parses explicit collision records and supplies the named `RUBIN2011` data set. It writes 66 unordered parent-pair collision entries and adds the selected Rubin Table 4 chemistry channels to the reaction-network generator. It emits `CometPhysics.generated.h`. |
| `input/test_comet_reaction_network.pl` | new | Tests reaction parsing, all supported rate-law syntax, loss channels, generated output, and rejection of unknown species. |
| `input/test_comet_physics_database.pl` | new | Tests the complete Rubin parent collision count, generated chemistry/collision data, the required OH-rate choice, methanol sensitivity handling, and custom collision records. |

### Reusable core and model-library files

| File | State | Why it changed and merge dependency |
| --- | --- | --- |
| `src/models/chemical_solver/ReactionNetwork.h` | new | Adds a small, header-only, validated reaction-data container with stoichiometric terms, constant/inverse-square/Arrhenius rate laws, mass-action evaluation, source/loss evaluation, and unit-unimolecular detection. It is independent of the incomplete legacy CHEMEQ2 driver. |
| `src/models/chemical_solver/test_reaction_network.cpp` | new | Unit-tests validation, rate evaluation, stoichiometry, source/loss construction, and unimolecular classification. |
| `src/models/chemical_solver/README.md` | new | Documents the supported network representation and the current limitation that general multi-reactant cell execution is not yet implemented. |
| `src/models/cometary_physics/CollisionDatabase.h` | new | Adds a header-only unordered species-pair database for constant and temperature-power collision cross sections, including explicit temperature validity ranges and distinct missing-pair/out-of-range results. |
| `src/models/cometary_physics/test_collision_database.cpp` | new | Unit-tests symmetric pair lookup, both collision laws, validation, missing pairs, and range enforcement. |
| `src/models/cometary_physics/README.md` | new | Documents database semantics, units, validation, and the rule against silent extrapolation or undocumented zero cross sections. |
| `src/meshAMR/meshAMRgeneric.h` | modified | Contains the largest core change. It separates the `findNeibNodeFace()` lookup method from the `neibNodeFace` storage member; implements DIM1 face lookup, spherical face area, neighbor refinement queries, dimension-specific ghost ranges, and symmetry checks; limits block allocation to `2*DIM` faces; fixes the maximum-refinement integer-width test; bounds maximum-refinement spacing scans by the compiled mesh dimension; and makes mesh diagnostics dimension-aware. All callers of the renamed method must be reconciled together. |
| `src/meshAMR/meshAMRcutcell.cpp` | modified | Bounds global-domain extent reads by `_MESH_DIMENSION_` while retaining three-component physical ray vectors. This removes optimized lower-dimensional out-of-bounds reads without changing DIM3 geometry. |
| `src/pic/parallel/init_send_halo.cpp` | modified | Updates the DIM1/DIM2/DIM3 halo constructors to call `findNeibNodeFace()`. This is an API-following change and must be merged with `meshAMRgeneric.h`. |
| `src/pic/pic.h` | modified | Uses the dimension-neutral `cAmpsMesh` alias in the particle tracker instead of the hard-coded `cMeshAMR3d` type, allowing the same code to compile in DIM1. |
| `src/pic/picGlobalDefinitions.dfn` | modified | Defines explicit built-in and user-defined model-particle cross-section modes. These constants are required by the collision adapter and must land with `picGlobal.dfn` and `pic_particle_collision_models.cpp`. |
| `src/pic/picGlobal.dfn` | modified | Selects the built-in molecular cross-section pointer by default. Application input can override it with the user-defined mode; Halley does so in `#General`. |
| `src/pic/pic_particle_collision_models.cpp` | modified | Routes every full/simplified NTC and MF cross-section evaluation through a common adapter. Built-in HS/VHS/VSS behavior is retained; user-defined laws receive velocity, species, block, and cell. This fixes the previously compiled-but-unreachable Rubin callback. |
| `src/pic/pic_parallel.cpp` | modified | Uses the common `cAMRnodeID` alias in buffered and unbuffered particle exchange, initializes three-component work vectors, and reads node coordinates/visits faces only in active dimensions. |
| `src/pic/pic_block_send_mask.cpp` | modified | Restricts send-mask neighbor scans to the active `2*DIM` faces and the number of subface slots represented in that dimension. This removes out-of-bounds DIM1/DIM2 neighbor reads and is required for multi-rank radial runs. |
| `src/pic/pic_mesh.cpp` | modified | Makes all four block pack/unpack and boundary send/receive paths visit positive faces only for active coordinates. Sender and receiver loop order must remain identical; this fix is required for load redistribution after a DIM1 sampling window. |
| `src/pic/pic_mesh_search.cpp` | modified | Traverses only child axes represented by the compiled mesh, copies only active node-bound components, and fixes inactive hash indices at zero. This removes optimized DIM1/DIM2 reads past node coordinate arrays. |
| `src/pic/pic_mover.cpp` | modified | Moves rotation-body declarations outside a dimension-only branch, corrects three duplicated z-coordinate typos, repairs variables used by the spherical second-order boundary-intersection calculation, and bounds node-coordinate scale calculations by the compiled dimension. The radial input still selects the tested first-order spherical mover; do not infer that the second-order spherical path is fully complete. |
| `src/pic/pic_field_line.cpp` | modified | Keeps physical three-component position, magnetic-field, and Sun vectors at length three even when the mesh is compiled with fewer dimensions. This is a lower-dimensional compile fix for code that still performs three-vector algebra. |

### Detailed core-source change record

This subsection records the core changes at symbol and dependency level so a
future merge does not have to infer their purpose from a large branch diff.
The comparison used here is local `master` at `cb7be2b` against branch HEAD
`5c01ecf` plus the current working tree. At that comparison point, 13 existing
files below `src/meshAMR` and `src/pic` differ from `master` (259 insertions and
120 deletions). Eight of those files contain additional uncommitted work from
the latest optimized-run audit. Recompute both the base and file counts before
merging because they cease to be authoritative if either branch moves.

#### New reusable core libraries

The following are additions rather than edits to an existing AMPS execution
path:

- `src/models/chemical_solver/ReactionNetwork.h` defines `RateLaw`,
  `cStoichiometricTerm`, `cReaction`, and `cReactionNetwork`. Validation checks
  indices, stoichiometry, rate-law parameters, and environmental requirements.
  `GetRateCoefficient()`, `GetReactionRate()`, and `EvaluateSourceLoss()`
  provide a data representation and rate evaluator; `IsUnitUnimolecular()`
  identifies the subset executed by the current particle chemistry adapter.
  It does not replace or repair the legacy CHEMEQ2 integration path.
- `src/models/chemical_solver/test_reaction_network.cpp` exercises validation,
  constant/inverse-square/Arrhenius rates, reaction rates, source/loss terms,
  and unit-unimolecular classification. Its adjacent `README.md` documents the
  unsupported general cell-solver case.
- `src/models/cometary_physics/CollisionDatabase.h` defines `CollisionLaw`,
  `cCollision`, and `cCollisionDatabase`. `Validate()` rejects malformed or
  duplicate unordered pairs, `Find()` performs symmetric lookup, and
  `GetCrossSection()` distinguishes missing data from a temperature outside a
  law's validity interval. The container does not choose a mixing rule or
  silently turn an absent pair into zero.
- `src/models/cometary_physics/test_collision_database.cpp` tests pair
  symmetry, constant and temperature-power laws, validation, missing-pair
  reporting, and temperature bounds. Its adjacent `README.md` records units
  and database semantics.

These headers have no AMR or PIC runtime dependency. Merge them before the
generic comet adapters so their unit tests can be reviewed independently.

#### Existing AMR and PIC files

- `src/meshAMR/meshAMRgeneric.h`: the AMR node member named
  `neibNodeFace` remains the neighbor storage, while face lookup is renamed to
  `findNeibNodeFace()` to remove the member/function collision that prevented
  current-core compilation. Dimension-specific implementations cover DIM1,
  DIM2, and DIM3. DIM1 gains radial neighbor lookup and spherical face area;
  neighbor refinement, ghost ranges, block allocation, MPI neighbor scans,
  and diagnostics use the face/child counts represented by the compiled mesh.
  The refinement-level guard now uses the usable bit width of `int`, and the
  latest audit restricts global-spacing initialization to
  `_MESH_DIMENSION_` entries instead of reading three values from a DIM1/DIM2
  array.
- `src/pic/parallel/init_send_halo.cpp`: all three halo constructors call
  `findNeibNodeFace()`. This is an API-following change and cannot be omitted
  when merging the method rename in `meshAMRgeneric.h`.
- `src/pic/pic.h`: the particle tracker uses the common `cAmpsMesh` type rather
  than a hard-coded `cMeshAMR3d` specialization. This is a type correction; it
  does not change the particle data layout.
- `src/pic/pic_block_send_mask.cpp`: `InitLayerBlockBasic()` and
  `InitLayerBlock()` scan `2*DIM` faces and only the valid subface count for
  that dimension. Previously DIM1 accessed faces and subfaces that do not
  exist, causing the observed four-rank `SIGSEGV` during exchange-mask setup.
- `src/pic/pic_mesh.cpp`: `PackBlockData_Internal()`,
  `UnpackBlockData_Internal()`, `sendBoundaryLayerBlockData()`, and
  `recvBoundaryLayerBlockData()` visit positive faces only along active axes.
  The sender and receiver order remains identical. This fixed corruption when
  a sampling break migrated DIM1 blocks between ranks.
- `src/pic/pic_parallel.cpp`: buffered and unbuffered particle exchange use
  the dimension-neutral `cAMRnodeID`. The latest changes zero inactive
  components of three-vector work arrays, compute cell sizes only from node
  coordinates that exist, and process only `2*_MESH_DIMENSION_` faces in
  `ProcessCenterBlockBoundaryNodes()`. Physical vectors stay length three;
  AMR bounds do not.
- `src/pic/pic_mesh_search.cpp`: `ScanTree()` descends only through child axes
  present in the AMR tree. `Init()` copies only active root bounds and gives
  inactive hash coordinates safe fixed values. `FindBlock()` fixes inactive
  hash indices at zero. These changes remove optimized DIM1/DIM2 reads past
  node-bound and child-pointer arrays.
- `src/meshAMR/meshAMRcutcell.cpp`:
  `CutCell::CheckPointInsideDomain_default()` computes the escape-ray length
  from only the domain extents stored by the mesh. Its ray and surface
  geometry remain three-dimensional. This is an out-of-bounds correction, not
  a conversion of cut-cell geometry to one dimension.
- `src/pic/pic_field_line.cpp`: work arrays used for physical position,
  magnetic field, and Sun direction remain three-component arrays even in a
  reduced-dimensional mesh. The change prevents vector algebra from writing
  beyond `DIM`-sized temporary arrays.
- `src/pic/pic_mover.cpp`: body-of-rotation declarations are available to the
  dimension-neutral boundary switch; three spherical rotations now use the
  actual z component instead of duplicating y; the second-order spherical
  block-intersection calculation uses `xInit` and `vMiddle` and declares its
  selected face; and the latest cut-surface block-scale calculation reads only
  active node-bound coordinates. The production radial path uses the
  first-order spherical mover. The second-order path compiles but still needs
  a dedicated runtime test before it is advertised as supported.

#### Collision-dispatch core files from the latest audit

These three files form one indivisible change:

- `src/pic/picGlobalDefinitions.dfn` adds symbolic values for the existing
  built-in molecular cross-section route and the cell-aware user-defined
  route.
- `src/pic/picGlobal.dfn` selects the built-in route by default, preserving
  existing applications unless their generated input explicitly overrides the
  mode.
- `src/pic/pic_particle_collision_models.cpp` adds
  `EvaluateCollisionCrossSection()` and routes all ten cross-section
  evaluations in simplified/full NTC, MF, and MF-Yinsi through it. The
  built-in branch continues to call
  `MolecularModels::GetTotalCrossSection()`. The user branch calls
  `ParticleCollisionModel::UserDefined::GetTotalCrossSection()` with both
  velocities, both species, the block, and the center cell.

The cell-aware branch is needed for laws such as Rubin's temperature-dependent
H2O-H2O cross section. Before this change, the comet callback was compiled and
linked but unreachable; NTC used the default `1 x 1` zero HS placeholder. The
Halley input override, `CometPhysics.cpp`, and the DIM3 callback in `Comet.cpp`
are application-side dependencies of this core change. A merge that retains
only the core adapter or only the application callback is incomplete.

#### Latest uncommitted core files

At the time of this record, the following core files have working-tree changes
after branch HEAD `5c01ecf`:

```text
src/meshAMR/meshAMRcutcell.cpp
src/meshAMR/meshAMRgeneric.h
src/pic/picGlobal.dfn
src/pic/picGlobalDefinitions.dfn
src/pic/pic_mesh_search.cpp
src/pic/pic_mover.cpp
src/pic/pic_parallel.cpp
src/pic/pic_particle_collision_models.cpp
```

Do not expect a later merge of `1d_comet_dev` to include these changes until
they have been committed on that branch.

### Core merge procedure

Before starting the merge, make the development branch clean and record the
exact bases:

```sh
git switch 1d_comet_dev
git status --short
git rev-parse master
git rev-parse HEAD
git diff --name-status master...HEAD -- src/meshAMR src/pic src/models
git diff --check
```

Commit the latest working-tree core fixes before switching branches. Keep the
collision-dispatch trio in the same commit, and do not mix generated `build/`
copies or run output into that commit. Once the branch is clean, merge it into
an integration branch based on the then-current `master`; use that integration
branch for conflict resolution and testing before updating `master`.

Resolve and review the core in these dependency batches:

1. **Reusable libraries:** add `ReactionNetwork.h`, `CollisionDatabase.h`,
   their tests, and their READMEs. Run their standalone tests before any
   application code is considered.
2. **AMR API and geometry:** merge `meshAMRgeneric.h` together with
   `parallel/init_send_halo.cpp` and `pic.h`. Search every caller after
   resolving the `findNeibNodeFace()` rename. Preserve newer upstream AMR
   behavior around the changed blocks rather than replacing the complete file
   with the development-branch version.
3. **MPI serialization and ownership:** merge `pic_block_send_mask.cpp`,
   `pic_mesh.cpp`, and `pic_parallel.cpp`. Inspect sender and receiver loops as
   a pair. A conflict resolution that changes only pack or only unpack can
   compile and still corrupt a migrated block.
4. **Reduced-dimensional bounds:** merge `meshAMRcutcell.cpp`,
   `pic_mesh_search.cpp`, `pic_field_line.cpp`, and the dimension-bound parts
   of `pic_mover.cpp`. Keep three-component loops for physical vectors and use
   `_MESH_DIMENSION_` only for AMR coordinate/bound arrays. Replacing every
   literal `3` mechanically would be incorrect.
5. **Collision dispatch:** merge `picGlobalDefinitions.dfn`, `picGlobal.dfn`,
   and `pic_particle_collision_models.cpp` together, followed by the comet
   callback adapters and the Halley input override. Confirm that the default
   remains the built-in route for applications that do not request a custom
   law.

Useful semantic checks after resolving conflicts are:

```sh
# Locate the old lookup spelling. The current tree has one commented example
# in init_send_halo.cpp; no compiled call may remain after conflict resolution.
rg '\bneibNodeFace\(' src

# Sender and receiver paths must retain the same active-face ordering.
rg -n 'PackBlockData_Internal|UnpackBlockData_Internal|sendBoundaryLayerBlockData|recvBoundaryLayerBlockData' \
  src/pic/pic_mesh.cpp

# All model-particle collision algorithms must use the common adapter.
rg -n 'MolecularModels::GetTotalCrossSection|EvaluateCollisionCrossSection' \
  src/pic/pic_particle_collision_models.cpp

# Generated Halley configuration must end with the user-defined override.
rg -n 'PARTICLE_COLLISION_MODEL__CROSS_SECTION' \
  build/pic/picGlobal.dfn build/pic/picGlobalDefinitions.dfn
```

After conflict resolution, remove `build` before every configuration compile,
as required for this application. Run the standalone reaction, collision, and
radial-profile tests; then perform clean DIM1 and DIM3 builds. The current
working tree has a clean optimized DIM1 build, but the latest uncommitted core
round has not yet received a clean DIM3 regression build. Do not merge to
`master` until that DIM3 build passes.

The DIM1 MPI acceptance run must cross an output boundary and actual block
redistribution, because a one-rank or pre-output run cannot exercise
`pic_block_send_mask.cpp` and the pack/unpack changes. Separately verify that a
user-defined collision case reaches its application callback. The present
Halley diagnostic does reach the Rubin callback but then stops because the
first collision step has no sampled H2O temperature; that expected physical
input gap is not a collision-dispatch regression.

Do not merge or commit `build/`, `PT*/`, `amps`, run logs, generated mesh and
signature files, root-level parser/input copies, or generated application
headers copied out of `build/main`. Regenerate those artifacts from the
tracked input and source files after the merge.

### Generic comet application files

| File | State | Why it changed and merge dependency |
| --- | --- | --- |
| `srcHartley2/makefile` | modified | Adds the generic model, chemistry, physics, radial mesh, and radial runtime objects to the application library. During DIM1 preprocessing, `hartley2.pl` replaces the DIM3 object list and main object with the smaller radial set and `RadialMain.o`; keep that substitution synchronized with this file. |
| `srcHartley2/main.cpp` | modified | Replaces the hard-coded `hartley2.shape.nas` kilometre reader with the mesh filename and length unit from `CometModel`. This preserves DIM3 behavior while allowing future comet meshes and units to be input data. |
| `srcHartley2/Comet.cpp` | modified | Removes the hard-coded H2O/CO2 collision table and delegates the legacy DIM3 AMPS collision callback to `Comet::Physics::GetNeutralCrossSection()`. DIM1 does not link this file, so `CometPhysics.cpp` supplies the callback only under `#if DIM == 1`; preserve that split to avoid either an undefined or duplicate symbol. |
| `srcHartley2/CometModel.h` | new | Declares geometry modes and per-species surface/extended source descriptors shared by all comet geometries. It is the stable interface between input-generated data and geometry-specific execution. |
| `srcHartley2/CometModel.cpp` | new | Stores generated model data and validates geometry, heliocentric distance, mesh units, source rates, source modes, temperatures, and bulk speeds before a run. Its initialized declarations are substitution targets for `hartley2.pl`; retain their spelling and one-line form. |
| `srcHartley2/CometChemistry.h` | new | Declares the generated-network adapter and the AMPS particle-reaction callback. It includes the reusable reaction header from either a generated or direct source-tree build. |
| `srcHartley2/CometChemistry.cpp` | new | Validates the generated network and executes supported unit-unimolecular events. It selects mutually exclusive channels, expands integer product multiplicity, handles explicit sinks, preserves event weight, and assigns isotropic product velocities with center-of-mass momentum and excess-energy conservation. |
| `srcHartley2/CometChemistry.generated.h` | new/generated template | Provides an empty source-tree fallback so direct compilation is possible. The preprocessor writes the case-specific version into `build/main`; never replace this template with a generated Halley table during conflict resolution. |
| `srcHartley2/CometPhysics.h` | new | Declares the collision-database adapter and the neutral-collision callback used by both DIM1 and DIM3. |
| `srcHartley2/CometPhysics.cpp` | new | Builds and validates the generated collision database, enforces missing-pair and temperature-range failures, evaluates temperature laws, and owns the AMPS callback only in DIM1. |
| `srcHartley2/CometPhysics.generated.h` | new/generated template | Provides the empty `NONE` database for direct source-tree builds. As with the chemistry template, case data is generated into `build/main` and must not be copied back from the build directory. |
| `srcHartley2/RadialMesh.h` | new | Declares explicit radial domain controls and the ideal logarithmic resolution profile used by the DIM1 mesh callback. Exactly one of radial extent above the surface and absolute outer radius is allowed. |
| `srcHartley2/RadialMeshProfile.cpp` | new | Solves the geometric cell-width ratio, maps cometocentric radius to ideal cell width, exposes ideal boundaries, and rejects invalid or decreasing-resolution requests. It contains no AMPS runtime dependency and is unit-testable by itself. |
| `srcHartley2/RadialMesh.cpp` | new | Connects parsed controls to AMPS, translates the ideal profile into a local-resolution callback, counts allocated dyadic AMR cells, and prints mesh diagnostics. |
| `srcHartley2/test_radial_mesh.cpp` | new | Tests expanding, uniform, one-cell, and invalid radial-profile configurations. |
| `srcHartley2/RadialMain.cpp` | new | Supplies the DIM1 executable entry point. It initializes AMPS, validates model/chemistry/physics/runtime data, builds or reads the distributed radial mesh, creates load lists, and starts the particle simulation only when enabled. |
| `srcHartley2/RadialSimulation.h` | new | Declares all required numerical controls and the AMPS callbacks for local timestep, injection, and domain-boundary processing. |
| `srcHartley2/RadialSimulation.cpp` | new | Implements spherical parent injection, species- and radius-local timesteps/weights, thermal inner re-emission, outer loss, daughter weight initialization, chemistry execution, MPI ownership synchronization, full radial-profile assembly, relative-L1 convergence, sampling-window output, load redistribution, and the hard iteration ceiling. |

### Documentation and reproducibility tools

| File | State | Why it belongs in the merge |
| --- | --- | --- |
| `srcHartley2/README.md` | new | User-facing build/run guide for generic DIM1 and existing DIM3 cases, including generated-input behavior and validation limits. |
| `srcHartley2/GENERIC_COMET_DESIGN.md` | new | Records the intended separation among geometry, source data, chemistry, collision physics, and case inputs, including the reserved DIM2 axisymmetric interface. |
| `srcHartley2/STATUS.md` | modified | Records implemented behavior, test evidence, known physical omissions, paper-derived parameters, and this merge map. |
| `srcHartley2/comet_shape_properties.py` | new | Reproducibly reads the Hartley NASTRAN surface and derives area, volume, topology checks, and equivalent radii. |
| `srcHartley2/shou2025_spherical_parameters.py` | new | Reproduces the Hartley spherical source temperatures, densities, mean-free-path estimate, and mesh-related scalar parameters from the documented inputs. |
| `srcHartley2/halley2011/PAPER_PARAMETERS.md` | new | Separates values stated by Rubin et al. from AMPS numerical choices and records equations and unresolved physics needed for a full paper reproduction. |
| `srcHartley2/halley2011/derive_halley1d_controls.py` | new | Reproduces the Halley surface-resolution estimate, particle-count choice, and daughter characteristic speeds from the selected reactions and excess energies. |
| `srcHartley2/visualization/plot_1d.py` | new | Plots named quantities from one or many DIM1 Tecplot files with aliases, radial filtering, explicit Y scaling, log-Y masking, filename legends, PNG/EPS output, and optional viewer launch. |
| `srcHartley2/visualization/compare_rubin2011.py` | new | Selects a complete Halley output index and creates paper-ordered neutral-density panels, H2O radial velocity, PNG/EPS files, a combined PDF, a radius-sampled CSV table, and a JSON manifest. It does not synthesize the paper's NMS instrument response or unpublished curves. |
| `srcHartley2/visualization/README.md` | new | Documents plotting aliases, CLI examples, output behavior, absent-species handling on logarithmic plots, and the Rubin comparison artifact set. |

No Python bytecode or temporary Python source is currently present under
`srcHartley2`. Keep `__pycache__/` and `*.pyc` excluded from commits.

### Merge dependencies and conflict handling

The safest integration order is:

1. Merge the reusable DIM1 core corrections in `src/meshAMR` and `src/pic`.
   Compile at least DIM1 and DIM3 before adding application code. In
   particular, the `findNeibNodeFace()` rename and every caller must land in
   the same commit.
2. Add `ReactionNetwork.h`, `CollisionDatabase.h`, their tests, and their
   READMEs. These components have no Hartley application dependency and can be
   reviewed independently.
3. Add the Perl generators and tests, then the generic `CometModel`,
   `CometChemistry`, and `CometPhysics` adapters. Verify that preprocessing
   writes headers only into `build/main`.
4. Merge the radial mesh and runtime files, followed by the
   dimension-selecting parts of `hartley2.pl` and `srcHartley2/makefile`.
5. Merge the Hartley DIM3, Hartley radial, and Halley inputs with the expanded
   species catalog. Add documentation and analysis/visualization scripts last.

The following conflict points require semantic review rather than accepting
one complete side of a conflict:

- In `meshAMRgeneric.h`, retain new upstream AMR behavior while applying the
  dimension-aware loop bounds, DIM1 spherical geometry, integer-width check,
  and neighbor-lookup naming consistently. After resolution, search for stale
  calls with `rg '\bneibNodeFace\(' src` and inspect every result.
- In `pic_mesh.cpp`, pack, unpack, boundary-send, and boundary-receive paths
  must use the same `DIM`-bounded face order. A one-sided resolution can
  compile and still corrupt load-migration messages.
- In `hartley2.pl`, keep upstream legacy keywords and finalization code while
  preserving generic-record parsing, validation flags, source-array writes,
  generated-header writes, and DIM1 object selection. The variable names used
  by `ChangeValueOfVariable` and `ChangeValueOfArray` must continue to match
  the initialized declarations in the application `.cpp` files exactly.
- In `Comet.cpp` and `CometPhysics.cpp`, retain a single definition of the
  AMPS user collision callback for each dimension: `Comet.cpp` owns DIM3 and
  delegates to the adapter; `CometPhysics.cpp` owns DIM1 under its preprocessor
  guard.
- Merge `picGlobalDefinitions.dfn`, `picGlobal.dfn`, and
  `pic_particle_collision_models.cpp` as one collision-dispatch change. If the
  mode constants or generated override are lost, a build can succeed while
  silently reverting to the zero placeholder table. Verify the final generated
  `build/pic/picGlobal.dfn` and run a callback-reaching diagnostic.
- In `species.input`, compare names before merging blocks. Species ordering is
  an input-case choice, but generated chemistry indices always follow the
  selected `SpeciesList`; regenerate after any ordering change.
- Treat `CometChemistry.generated.h` and `CometPhysics.generated.h` in
  `srcHartley2` as empty templates. The populated copies in `build/main` are
  build products and are never merge inputs.
- Review the root `Makefile` linker flag independently. There is no physical
  model dependency on `-lmpi++`, and a different MPI installation may not
  provide that compatibility library.

Before committing or merging, exclude `build/`, `PT/`, `output/`, `amps`,
`comet.radial.mesh.msh`, `error.log`, `gitinfo.txt`, root-level parser copies
such as `halley1d.input`, `hartley2.input`, `hartley2.pl`, and
`species.input`, and generated convenience links/directories such as
`SWMF_data`, `data`, `share`, `util`, and `table_parser`. These files are
created or exposed by preprocessing, compilation, or execution and are not
part of the source change set.

After resolving the merge, rerun the complete commands in **Verification**.
The minimum acceptance evidence is a clean DIM1 build, a clean DIM3 build,
all five standalone parser/database/profile tests, a multi-rank DIM1 run that
crosses a sampling/output/load-redistribution boundary, and
`git diff --check`. A successful compile alone does not exercise the two MPI
defects fixed in `pic_block_send_mask.cpp` and `pic_mesh.cpp`.

## Changes in this work

### Generic comet configuration

`CometModel.h/.cpp` now holds geometry, heliocentric distance, and per-species
surface/extended source descriptors. `hartley2.pl` accepts:

- `Dimension` and `Geometry`;
- `HeliocentricDistanceAU`;
- `SurfaceMeshFile` plus `M` or `KM` units;
- `SurfaceSource` records by species;
- `ExtendedSource` records by species;
- `ChemistryMode` and repeated `ChemistryReaction` records;
- `PhysicsDatabase` and its database-specific sensitivity choices; and
- `RadialSimulationMode` plus explicit DIM1 numerical controls.

Input species names are checked against `SpeciesList`. Rates, temperatures,
speeds, dimensions, mesh parameters, distributions, and chemistry species are
validated. The former hard-coded `hartley2.shape.nas` read in `main.cpp` now
uses the configured mesh name and length conversion.

`CometModel`, `CometChemistry`, and `CometPhysics` are compiled into both
dimension targets. The DIM1 main validates model, chemistry, collision, and
radial-run data before building the mesh.

### Chemistry network

The new core header `src/models/chemical_solver/ReactionNetwork.h` provides a
validated flat reaction database, constant/inverse-square/Arrhenius rates,
mass-action rate evaluation, source/loss evaluation, and identification of
unit unimolecular reactions.

`input/CometReactionNetwork.pm` parses arbitrary reaction stoichiometry from
the comet block and generates `CometChemistry.generated.h`. The checked-in
generated header is an empty fallback for direct source-tree builds. The
network is data-driven, but it does not make the incomplete legacy CHEMEQ2
code operational. Multi-reactant execution still requires a cell chemistry
solver.

`CometChemistry.cpp` implements the unit-unimolecular PIC executor. It samples
one mutually exclusive channel for a parent event, expands integer
stoichiometry into concrete daughters, removes neutral parents for explicit
ionization sinks, and generates isotropic daughter velocities with exact
center-of-mass momentum and tabulated excess kinetic energy. Product weights
preserve the physical parent-event weight. Startup validation rejects
noninteger product multiplicities, nonconserving neutral channels, and
multi-reactant networks.

### Generic collision physics

`src/models/cometary_physics/CollisionDatabase.h` provides validated,
unordered species-pair lookup for constant and temperature-power-law cross
sections. Each temperature law contains its cited validity interval. A missing
pair and an out-of-range temperature are separate error states; the
application adapter aborts on either rather than assigning zero or
extrapolating.

`input/CometPhysicsDatabase.pm` adds the `RUBIN2011` named database. It
transcribes the 66 parent-parent collision pairs in Table 3 and the neutral
photodissociation/photoionization channels in Table 4. The Crifo water
self-collision law is restricted to its published 50--300 K interval. Table 3
does not provide daughter-pair values, and those entries remain absent.

### Radial particle execution

`RadialSimulation.cpp` completes the generic DIM1 execution path when enabled:

- uniform parent rates are injected once through the spherical inner face;
- a Poisson process samples fractional-step injection times;
- species-local timesteps use an explicit characteristic speed and Courant
  number;
- species-local particle weights grow with the local timestep, keeping the
  expected model population approximately uniform per logarithmic cell;
- returning gas is fully thermalized and re-emitted at the inner boundary;
- particles crossing the outer boundary are deleted; and
- daughter species receive initialized base weights for exact event-weight
  correction;
- complete radial density profiles are assembled across MPI ranks after each
  sampling window using the core's spherical shell volumes; and
- the run terminates after every species passes the configured relative L1
  profile test for the required number of consecutive windows, or when the
  hard iteration ceiling is reached.

The backend requires its numerical controls in input. It supplies no hidden
defaults for speeds, timestep fraction, boundary temperature, iteration count,
sampling length, tolerance, or converged-window count. Block redistribution at
an output boundary does not change the comparison: every rank traverses the
same static AMR tree order and `MPI_Allreduce` reconstructs the global profile.

### Halley / Rubin et al. (2011)

The Rubin paper models Comet 1P/Halley. Its files are named accordingly:

```text
input/halley1d.input
srcHartley2/halley2011/PAPER_PARAMETERS.md
```

The Hartley 2 inputs remain `input/hartley2.input` and
`input/hartley2.radial.input`.

The Halley input contains the paper-derived 6 km nucleus radius, 0.90 AU
heliocentric distance, `10^6 km` cometocentric outer radius, 180 K zero-drift
uniform parent boundary, all 27 neutral species, and all 11 parent rates. The
rates sum to `8.54e29 s^-1`. Dust and extended sources are not enabled in this
gas-only input.

The paper describes the surface mean free path only as a fraction of a metre.
Using its parent rates, 6 km radius, 180 K half-Maxwellian boundary, and Table
3 H2O-parent cross sections gives a documented conventional estimate of
`0.329 m` (`0.466 m` without the equilibrium `sqrt(2)` factor). No exact mesh
width is attributed to the paper.

The requested 2,000 ideal cells and the `0.3293692849882 m` surface estimate
give a geometric width ratio of `1.008573952117` and require maximum AMR
refinement level 29. `ReferenceInjectionParticleNumber=16` gives an estimated
minimum of 52.454 direct-parent particles in the surface cell. The Courant
number is 0.5. Direct parents start from the supplied 400 m/s characteristic
speed; daughter recoil bounds are propagated through the reaction graph from
Table 4 excess energies by `halley2011/derive_halley1d_controls.py`.

The methanol-to-formaldehyde branch is disabled. The particle buffer begins at
five million slots, consistent with the paper's statement that the full model
uses several million representative particles, and must be monitored in a
convergence run. Species-local weights proportional to timestep prevent
particle count from increasing with radial cell width.

The clean one-process smoke run allocated 373 blocks and 2,984 radial cells;
its smallest cell is `0.2328292 m`. Summing the derived no-loss parent
occupancies over those cells gives approximately `2.21e6` direct-parent model
particles, leaving about `2.79e6` buffer slots for daughters and transients.

Rubin Table 4 gives the 0.36 eV OH branch rate as `[3e-6,6e-6] s^-1`; the
requested input selects its lower endpoint, `3.0e-6 s^-1`. Its other neutral
OH branches sum to `2.2e-6 s^-1`, and OH ionization is `1.2e-5 s^-1`. The
sampling window is 500 iterations and the hard ceiling is 10,000 iterations.
The run stops earlier when every species' shell-volume-weighted relative L1
radial density-profile change is at most 5% for three consecutive completed
window comparisons.

The paper includes Crovisier/Huebner radiative cooling and collisional
relaxation. The legacy DIM3 source contains a cooling routine, but the generic
DIM1 path has not yet been connected and validated. This is recorded as a
physical blocker rather than silently omitted from a claimed reproduction.

### Hartley 2 paper configuration

`input/hartley2.radial.input` now compiles H2O and CO2 and records values from
Shou et al. (2025):

```text
heliocentric distance       1.06 AU
H2O direct source           2.0e27 s^-1
CO2 direct source           3.5e26 s^-1
H2O effective temperature   188.350171377 K
CO2 effective temperature   182.852991011 K
surface bulk drift          0 m/s
surface mean free path      1.004209133 m (core-derived estimate)
```

Dust and both extended gas sources are omitted from this case. `DustMode=off`
sets the AMPS core switch explicitly, and the lack of `ExtendedSource` records
leaves every extended-source rate and mode disabled.

The direct sources are uniform because the requested DIM1 geometry is
spherically symmetric. `shou2025_spherical_parameters.py` reproduces their
scalar temperatures from the 1.29-AU solar-angle profile cited by the paper and
stored in `CG-BC.h`. For species night-side factor `a`, it evaluates

```text
T_eff = integral(T(theta) max(cos(theta),a) sin(theta) dtheta)
        / integral(max(cos(theta),a) sin(theta) dtheta),  0 <= theta <= pi.
```

The paper gives `a=0.1` for H2O and `a=0.3` for CO2. This number-weighted mean
also preserves the emitted thermal energy of a zero-drift half-Maxwellian. The
existing Hartley2 DSMC injection code passes a zero bulk vector to that
distribution, so the input now uses `0.0 m/s`. The paper's `600 m/s` is the
constant speed in its analytic inversion flux relation and Haser profiles; it
is not the drift passed to the DSMC surface sampler.

`comet_shape_properties.py` reads plain or gzip-compressed NASTRAN meshes and
reports topology, bounds, surface area, volume, and equivalent radii. For the
tracked Hartley 2 mesh:

```text
surface area                 5,242,876.570262 m^2
area-equivalent radius             645.921714366 m
volume-equivalent radius           578.189522337 m
```

The DIM1 input uses the requested area-equivalent radius, a supplied radial
extent of 10,000 km above the surface, and 1,000 ideal radial cells. The paper
does not specify the nucleus-surface resolution or a collision cross section.
The configured 1 m surface target comes from the mean-free-path estimate below.

### Surface mean-free-path estimate

The closest reusable collision prescription in the AMPS core is
`input/comet_sphere.input`, which assigns a constant `3.4e-19 m^2` cross
section to H2O-H2O, H2O-CO2, and CO2-CO2 collisions. Shou et al. do not state
a cross section or surface mean free path.

For an upward, zero-drift half-Maxwellian, the outgoing density of species `i`
immediately above the surface is estimated from its uniform spherical flux as

```text
Gamma_i = Q_i / A
u_normal,i = sqrt(2 k T_i / (pi m_i))
n_i = Gamma_i / u_normal,i.
```

The configured sources give

```text
H2O outgoing density        1.620900853e18 m^-3
CO2 outgoing density        4.501078225e17 m^-3
total outgoing density      2.071008675e18 m^-3
```

Using the conventional hard-sphere estimate
`lambda = 1/(sqrt(2) n sigma)` gives `lambda = 1.004209133 m`. The 1 m target
therefore resolves one estimated surface mean free path. The `sqrt(2)` form is
exact for an equilibrium single-species Maxwellian; its use for this outward,
two-temperature mixture is a documented resolution estimate. A
species-dependent collision model would be needed for a more exact value.

### Radial mesh

For surface width `dx0`, requested ideal count `N`, extent `L`, and expansion
ratio `q`, `RadialMeshProfile.cpp` solves

```text
dx_i = dx0 q^i
L = dx0 (q^N - 1) / (q - 1).
```

AMPS uses dyadic block refinement and approximates the ideal geometric
profile. The requested science-domain parameters are:

```text
nucleus radius                  645.921714366 m
radial extent above surface  10000000 m
outer radius                 10000645.921714366 m
requested surface width             1 m
requested ideal cells             1000
maximum AMR refinement level        21
geometric width ratio       1.011741915188
ideal outer-cell width       116057.4 m
allocated AMR blocks              189
allocated AMR cells              1512
smallest allocated cell width       0.5960464 m
```

With eight cells in the root block, level 20 has a minimum cell width of
`1.192 m` and cannot enforce the 1 m request. Level 21 permits `0.596 m` cells
and is the minimum sufficient AMR limit for this domain. The clean DIM1
runtime completed mesh construction with the values above.

The mesh driver now writes `comet.radial.mesh.msh`, reflecting its generic
role.

### Existing core dimension fixes

The earlier `1d_comet_dev` commit remains required. It repaired dimension
neutral AMR/PIC paths in:

- `src/meshAMR/meshAMRgeneric.h`;
- `src/pic/pic.h` and `src/pic/pic_parallel.cpp`;
- `src/pic/parallel/init_send_halo.cpp`;
- `src/pic/pic_mover.cpp`; and
- `src/pic/pic_field_line.cpp`.

Those changes implement DIM1 face-neighbor handling, dimension-aware node
indices and diagnostics, spherical face areas, current AMR node IDs, and
several mover/vector fixes. The radial input still selects the implemented
first-order spherical mover because the second-order spherical branch is not
complete.

The 10,000 km domain with a 1 m surface target needs refinement level 21. A
core validation check in `meshAMRgeneric.h` formerly compared that level with
`1 << sizeof(int)`, using the byte count as though it were a bit count and
incorrectly rejecting all levels 16 and above on a four-byte `int`. The check
now uses the number of usable signed-integer bits. Level 21 is representable;
levels that would make the core's `1 << level` expressions overflow remain
rejected.

The DIM1 particle smoke test also reached `AllocateTreeBlocks`, which exposed
a second core dimension defect in `cMeshAMRgeneric::AllocateBlock`: its common
node-reuse path always iterated six faces. A DIM1 tree node stores only its two
radial face neighbors, so faces 2--5 read beyond `neibNodeFace` and produced an
invalid pointer. The loop now uses `2*_MESH_DIMENSION_`; DIM2 and DIM3 retain
their four and six physical faces respectively.

After that core correction, a temporary Halley plumbing case completed one
iteration on one MPI process. The test used 100 ideal cells, a 1 m surface
target, refinement level 27, reference injection count 100, particle-buffer
capacity 100,000, Courant number 0.5, 5 km/s characteristic speeds, the lower
OH rate endpoint, and the methanol branch enabled. Those values verify program
flow only; they are not retained in the paper input and are not claimed as a
physical or converged Halley solution. The run built 28 blocks (224 AMR cells),
initialized all 11 parent weights, injected/moved particles, and returned
success after iteration 0.

## Verification

```sh
perl -c input/hartley2.pl
perl -c input/CometReactionNetwork.pm
perl -c input/CometPhysicsDatabase.pm
perl input/test_comet_reaction_network.pl
perl input/test_comet_physics_database.pl

g++ -std=c++17 -Wall -Wextra -pedantic \
  -Isrc/models/chemical_solver \
  src/models/chemical_solver/test_reaction_network.cpp \
  -o /tmp/test_reaction_network
/tmp/test_reaction_network

g++ -std=c++17 -Wall -Wextra -pedantic \
  -Isrc/models/cometary_physics \
  src/models/cometary_physics/test_collision_database.cpp \
  -o /tmp/test_collision_database
/tmp/test_collision_database

g++ -std=c++17 -Wall -Wextra -pedantic \
  -IsrcHartley2 srcHartley2/RadialMeshProfile.cpp \
  srcHartley2/test_radial_mesh.cpp -o /tmp/test_radial_mesh
/tmp/test_radial_mesh

python3 srcHartley2/comet_shape_properties.py \
  input/hartley2.shape.nas.gz --input-unit km
python3 srcHartley2/shou2025_spherical_parameters.py

rm -rf build
make -j8                         # DIM1 input selected
mpirun -np 1 ./amps              # mesh-only Hartley 2 radial input

rm -rf build
make -j8 EXTRACOMPILEROPTIONS=-O3  # final halley1d.input: 500/10000, 5%, 3 windows

rg 'PARTICLE_COLLISION_MODEL__CROSS_SECTION' \
  build/pic/picGlobal.dfn build/pic/picGlobalDefinitions.dfn

git diff --check
```

The final clean optimized Halley build generated `SamplingLength=500`,
`MaxIterations=10000`, `ConvergenceTolerance=0.05`,
`RequiredConvergedWindows=3`, and the 0.36 eV OH rate `3.0e-6 s^-1`; it linked
successfully. Its generated cross-section macro selects the user-defined
Rubin callback. A four-rank callback diagnostic then stopped at iteration 1
because the first collision step has no sampled H2O translational temperature.
The earlier multi-window output is invalid for physical comparison because it
used the zero placeholder cross section. The MPI compiler reports local socket
warnings in the restricted environment; they do not prevent compilation or
runtime. DIM3 retains legacy warnings for `register`, format widths, narrowing
conversions, and fixed-buffer `sprintf` calls.

## Remaining decisions for DIM1 science simulations

For Hartley 2, the paper, existing implementation, and supplied domain controls
determine the direct source rates, spherical effective temperatures, zero
surface drift, species, heliocentric distance, area-equivalent radius, domain
extent, cell count, and the core-derived surface resolution. Remaining controls
for a science run are:

- particle weights, timestep policy, sampling duration and convergence rule;
- whether to retain the provisional constant H2O/CO2 collision cross sections,
  and the outer-boundary behavior; and
- a chemistry network if chemistry is to be enabled.

Chemistry is off because Shou et al. do not publish a reaction network or rate
table for this run.

For Halley, the paper fixes the parent sources, domain endpoints, and parent
physics but leaves the input choices recorded in
`halley2011/PAPER_PARAMETERS.md`. The physically blocking items are the
daughter collision matrix, an explicit initialization/evaluation rule for the
temperature-dependent H2O collision law before cell sampling is available,
and a validated DIM1 radiative cooling/internal-energy path. Paper-curve comparison also needs
traceable source data or digitization with uncertainty.
