# Generic comet / Hartley 2 status

Date: 2026-11-01
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
| Rubin chemistry/collision generator Perl test | Pass; full hierarchy contains 378 pairs partitioned as 66 Rubin, 18 additional Combi, and 294 CHEMKIN/Lorentz estimates |
| Shape-property extraction | Pass; closed 32,040-triangle surface |
| Full DIM3 runtime | Not run; required Hartley 2 mission SPICE kernels are unavailable locally |
| Hartley 2 DIM1 particle simulation | Backend present; input remains mesh-only pending explicit science-run controls |
| Halley DIM1 one-step plumbing test | Pass with documented temporary numerical choices |
| Halley 2,000-cell/local-weight clean build and one-step run | Pass; 373 blocks and 2,984 cells |
| Halley DIM1 four-rank runtime | Reached iteration 3500 and seven complete 500-step outputs at Courant 0.5; then a fast particle traversed beyond the mover's one-block MPI boundary allocation |
| Completed-sample kinetic temperature | Pass; the Crifo law reads `GetTranslationalTemperature()` from the completed buffer and uses `max(Tsample,180 K)` |
| H2O internal degrees of freedom | Pass; core Larsen-Borgnakke rotation/RT relaxation enabled and initialized at injection/re-emission |
| Generic radiative cooling | Pass; named H2O Gombosi 1986 law unit-tested, clean-built, and executed after the first completed sample |
| Halley DIM1 convergence run | Courant-0.25 restart active; clean `-O3` four-rank run crossed four output/load-balance boundaries through iteration 2,000, currently transient at 0/3 convergence passes |
| Core collision cross-section dispatch | Fixed and clean-built; NTC/MF now route the selected user law with block/cell state instead of always using the zero HS placeholder |
| Halley collision-data coverage | Complete by explicit policy: 66 Rubin rows, 18 additional Combi rows, and 294 labeled CHEMKIN/Lorentz estimates |
| Completed-hierarchy clean build and smoke test | Pass; clean `-O3` DIM1 compile/link and isolated one-rank execution through iteration 124 |
| DIM1 result visualization | Pass for short aliases, single-/multi-file PNG/EPS output, and viewer control |
| Rubin paper-comparison artifact generator | Pass: every panel overlays AMPS with vector-extracted Rubin DSMC data; the available AMPS window is transient, not converged |
| Local-timestep face-flux transport | Pass for clean `-O2` DIM1 compile/link and a final two-rank smoke test through the start of iteration 153. The generated mode quartet is species-local timestep, species-global weight, correction off, and face resampling on. No runtime error occurred; timeout ended the test. A 500-step/MPI-ownership-boundary/long-population test remains required. |
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

### Local-timestep implementation added on 2026-11-01

The new transport path is opt-in through
`LocalTimeStepFluxResampling=on`. The Halley case combines species-local block
timesteps with one spatially global weight per species and disables individual
particle-weight corrections. A trajectory is integrated only to the next
block face. If `tau_old` is its unused old-block interval, the destination
interval is `tau_new=tau_old*dt_new/dt_old`; thus the unfinished fraction of a
local step is unchanged.

At the face, residual resampling draws an integer population with mean

```text
mu = (W_old/dt_old)/(W_new/dt_new)
   = W_old*dt_new/(W_new*dt_old).
```

The result is `floor(mu)` plus one Bernoulli draw for the fractional part.
The operation is unbiased for number, mass, momentum, and kinetic-energy flux
because every sibling retains the face position and velocity. It reduces to
`mu=dt_new/dt_old` for the configured per-species global weights. A deleted
particle terminates at roulette; surviving particles and clones carry their
normalized remainder in a dynamically reserved particle-data field.

AMPS then executes exchange/continuation sweeps until an `MPI_Allreduce`
reports no unfinished particles anywhere. Collisions, comet chemistry,
sampling, and the generic particle splitter remain outside that loop. This is
essential: a multi-block trajectory must receive one, not one-per-segment,
collision/chemistry/sampling pass during an outer iteration. MPI particle
serialization and restart records already copy the complete dynamically sized
particle record, so the new remainder field needs no separate message format.

The following files are modified by this implementation and must be reviewed
as one dependency set:

| File | Detailed reason for the change |
| --- | --- |
| `ampsConfig.pl` | Parses `LocalTimeStepFluxResampling=on/off` in `#main` and writes the selected core macro. The default remains off for every existing AMPS application. |
| `src/pic/picGlobal.dfn` | Defines the opt-in mode as off in the source core and documents its steady-state local-pseudo-time interpretation and configuration constraints. |
| `src/pic/pic.h` | Declares the continuation-field offset/accessors, unfinished-particle mover sweep, pending-particle counter, opt-in flux/lifecycle diagnostics, and validation-only mover override. These declarations are required by the allocator, mover, timestep driver, and focused application regression. |
| `src/pic/pic_mover.cpp` | Enforces compatible compile-time modes; reserves the optional particle field; implements `W/dt` residual resampling and normalized remainder scaling; stops first- and second-order trace movers at each face; normalizes DIM1 spherical position/velocity before cloning; resumes unfinished particles; and records flux moments and continuation lifecycle when diagnostics are enabled. A shared positive spherical-root solver fixes inward crossings, and the second-order DIM1 path now performs midpoint and no-intersection radial lookups. Production mover selection is unchanged unless validation explicitly installs an override. |
| `src/pic/pic_pbuffer.cpp` | Clears the continuation marker whenever a deleted particle slot is reused. Without this initialization, a newly injected or chemistry-created particle could inherit an old unfinished trajectory. |
| `src/pic/pic_time_step.cpp` | Exchanges main-mover face stops, then repeats continuation movement followed by an unconditional ownership exchange and global pending reduction. The post-move exchange retains completed particles deposited in remote boundary blocks. It records lifecycle diagnostics, accounts for timing, and delays splitting/list validation until no remainder remains. The legacy path remains under the default-off branch. |
| `input/halley1d.input` | Selects `SpeciesLocalTimeStep`, `SpeciesGlobalParticleWeight`, correction off, and the opt-in face algorithm. These settings are part of one physical/statistical policy and should not be merged independently. |
| `srcHartley2/RadialSimulation.cpp` | Initializes direct-source global weights from their requested injection statistics and gives non-source daughters the largest direct-parent weight. This bounds instantaneous chemistry clone bursts without changing expected production. The older local-weight and single-global modes remain compile-time alternatives. |
| `srcHartley2/CometChemistry.cpp` | Converts every stoichiometric product from parent rate `W_parent/dt_parent` to product rate `W_product/dt_product`. Correction-off builds use residual split/roulette and conserve in expectation; correction-on builds use one daughter with the corresponding exact correction. Differing species timesteps prevent a literal synchronized event interpretation. |
| `srcHartley2/LocalTimeStepValidation.h`, `srcHartley2/LocalTimeStepValidation.cpp` | Add the focused serial/MPI statistical regression: timestep-ratio and direction matrices for both movers, same-owner and ownership faces, particle state/list/lifecycle checks, flux moments with combined relative/statistical tolerance, and a multi-face/two-ownership-boundary trajectory. |
| `srcHartley2/RadialMain.cpp` | Adds validation CLI parsing and dispatch after normal mesh/core initialization but before the comet production loop. |
| `input/hartley2.pl`, `srcHartley2/makefile` | Include the validation object in the generated DIM1 application object list and direct application build. |
| `srcHartley2/README.md` | Documents the equations, sweep ordering, supported modes, steady-state meaning, limitations, and improvement paths. |
| `srcHartley2/STATUS.md` | Records implementation, verification evidence, known limits, and the master-integration procedure. |

The implementation intentionally does not claim synchronized transient time.
An outer iteration advances each block by its own local pseudo-time. It also
does not solve the population-growth cost of a spatially global weight: a
steady outward flux requires model population proportional to local timestep,
which becomes expensive across the full logarithmic Halley domain. Acceptable
future improvements are conservative face-population reservoirs, true local
weights with auditable weight transfer, or asynchronous block clocks. A hard
clone cap is not acceptable because it would bias flux.

The first clean two-rank smoke build used one weight shared by all species. It
reached iteration 29 without a mover or exchange failure, but by iteration 10
held 341,675 particles and injected about 34,307 particles per iteration,
including 281,261 H2O particles. That interpretation is computationally
unusable with a five-million-slot buffer and is not retained. The test ran in
`/tmp/amps-local-timestep-test.G72Idf` and did not alter repository results.

The retained one-weight-per-species configuration was then rebuilt cleanly
from a removed `build` directory. A two-rank 180-second transport run reached
the start of iteration 267 without an AMPS failure; that binary preceded the
final species-changing `W/dt` chemistry correction and is transport evidence
only. After inspecting the core exosphere chemistry convention, product
resampling was corrected from `W_parent/W_product` to
`(W_parent/dt_parent)/(W_product/dt_product)`. `build` was removed again and
the final source cleanly compiled and linked with `-O2`. Its two-rank test in
`/tmp/amps-local-timestep-final.bbNFsN` reached the start of iteration 153 in
30 seconds and ended only because of the external timeout. At iteration 10 it
held 887 particles and injected 89.4 particles per iteration. Neither retained
test reached the 500-step output/load-balance boundary or the MPI ownership
boundary, so those checks remain in the pre-merge acceptance gate.

### Focused local-timestep validation completed on 2026-11-01

The transport regression is implemented in
`LocalTimeStepValidation.{h,cpp}` and selected by CLI in `RadialMain.cpp`.
The final executable was produced only after removing `build` and cleanly
running `make -j8` with `halley1d.input`. Validation ran from
`/tmp/amps-lts-validation.0L6Pni`, not the repository run directory, so it did
not modify `PT` or any paper-comparison output.

The final commands were:

```sh
mpirun -np 1 /data/vtenishe/COMET1D/AMPS/amps \
  --validate-local-timestep --validation-particles 200000 \
  --validation-mover both --validation-skip-multiface

mpirun -np 4 /data/vtenishe/COMET1D/AMPS/amps \
  --validate-local-timestep --validation-particles 200000 \
  --validation-mover both
```

The one-rank command passed all 20 same-owner cases: five destination/source
timestep ratios (`0.25`, `0.5`, `1`, `2`, `2.5`), both radial directions, and
both first- and second-order trace movers. The four-rank command passed those
20 cases plus 12 MPI-face cases (`0.5`, `1`, `2.5`, both directions and both
movers), for 32 single-face cases total. Each case used 200,000 input
particles split among four distinct radial/transverse velocities and checked:

- expected residual split/roulette population and list membership;
- destination block, finite state, unchanged speed, and analytic position;
- one continuation pass, no continuation deletion, and zero final remainder;
- no loss between the first face stop and the first ownership exchange; and
- number, mass, all three momentum components, and kinetic-energy flux within
  the larger of 2% relative error and the five-sigma Monte Carlo uncertainty.

The four-rank long-trajectory test passed with 32 particles crossing 95 block
faces, including two MPI ownership transitions, in exactly 95 continuation
sweeps. It recorded 95 transfers per input particle, retained all particles,
and ended with no positive remainder marker.

The focused harness exposed and drove correction of four independent defects:

1. The old inline spherical-root calculation mutated coefficients between
   inner- and outer-face tests and could repeatedly select the same face for
   inward motion. Both movers now use a shared stable first-positive-root
   helper with the full three-component physical state.
2. The second-order DIM1 spherical path lacked its no-intersection radial node
   lookup and did not force face tests when the legacy face-check macro was
   off. Its midpoint/final radial lookups and opt-in face check are implemented.
3. DIM1 stores radial position only. Cloning a face state before rotating the
   velocity into that stored frame made some MPI clones inconsistent with
   their position. The state is normalized before split/roulette and copying.
4. A continuation could finish with marker zero in a remote boundary-layer
   block. Counting only positive markers before another exchange stranded that
   completed particle outside its owner list. Every continuation mover pass
   now performs an unconditional ownership exchange before the pending count.
   The formerly failing inward MPI ratio-2.5 case now reports all 500,230
   generated destination particles and zero pending work.

This completes the focused face-transport item of the acceptance plan. It does
not replace the production-run requirements: the Halley case must still run
through sampling and load redistribution, demonstrate manageable population
and convergence, exercise unequal-weight chemistry statistically, and produce
a converged Rubin comparison before any merge to `master` is considered.

### Courant-0.5 mover failure and restart policy

The repaired, capped-temperature, four-rank hierarchy run completed seven
500-iteration windows through iteration 3500 and wrote all 27 species for
output indices 0--6. It held 588,275 model particles at the last completed
window. The last available convergence records were `0/3` consecutive passes,
with maximum relative L1 change `0.8358229` for N at iteration 3000; the data
are therefore transient and must not be represented as a converged Halley
solution.

Immediately after the iteration-3500 output and load redistribution, the
legacy trace mover stopped in `pic_mover.cpp` because its current node had no
locally allocated block. The radial case uses eight cells per block and the
local timestep is `C*cell_width/characteristic_speed`. At `C=0.5`, a particle
moving faster than 16 times its species characteristic speed can traverse a
whole block in one move. Hard-sphere transfer from the roughly 57 km/s H
daughter tail can put a heavy collision partner near 6 km/s, close to 15 times
the 400 m/s parent baseline, leaving no stochastic margin. The failure is
therefore consistent with a collision-product tail crossing the one-block MPI
boundary allocation rather than with a mesh, chemistry, cooling, or collision
temperature error.

The restart uses `RadialCourantNumber=0.25`. A one-block traversal then
requires more than 32 times the configured characteristic speed, retaining
margin above the expected H-to-heavy transfer while preserving the explicitly
derived species characteristic speeds. The reference injection count is
reduced from 16 to 8 with the timestep so the requested roughly 50 parent
particles per surface cell are unchanged. The 500-iteration sampling length,
convergence test, and 10,000-iteration ceiling are unchanged.
The failed run and its Rubin comparison are archived rather than overwritten.
The exact paths are `PT.failed-mover-timestep-courant-0.5-20261028`,
`halley1d.failed-mover-timestep-courant-0.5-20261028.log`, and
`error.log.failed-mover-timestep-courant-0.5-20261028`.

After removing `build`, the Courant-0.25 case was rebuilt cleanly with `-O3`
and restarted on four MPI ranks from the repository root on 2026-10-31. It
crossed the complete output/load-redistribution boundaries through iteration
2,000 without either the earlier Crifo-temperature stop or mover stop. At
iteration 2,000 the maximum relative L1 change was `1.008269` for H2, with
`0/3` consecutive passes, so the new data remain
transient. The active log is
`halley1d.courant-0.25.run.log`; runtime data and the current Rubin-overlaid
comparison are under `PT/plots` and `PT/plots/rubin2011_comparison`.

### Full-hierarchy run stopped at the Crifo temperature limit

The clean optimized binary was launched from the repository root on four MPI
ranks with the generated 378-pair Rubin--Combi--CHEMKIN/Lorentz collision
database. It completed the sampling windows at iterations 500 and 1,000,
wrote all 27 species files for output indices 0 and 1, and reached 152,787
model particles. The first window established the convergence reference, so
two windows are insufficient for the configured three-consecutive-window
termination test.

Immediately after the second output and block redistribution, the collision
adapter encountered a completed-sample H2O temperature of
`395.80384407760528 K`. The H2O--H2O cross section is the Crifo relation quoted
by Rubin et al.; Rubin explicitly limits its cold-water-vapor use to
`50--300 K`. `CometPhysics.cpp` therefore stopped instead of silently
extrapolating. This is the intended failure mode for an undefined physical
policy, rather than a collision-database, MPI, memory, or output failure.

The second H2O output contains 41 sampled cells above 300 K, concentrated at
roughly 31--37 km cometocentric radius. These cells are statistically sparse:
the highest reported temperature is 487.525 K at 32.4643 km with only 0.282
sampled model particles in that cell. At the abort temperature, direct
extrapolation would give `1.405700619e-19 m^2`, while limiting the input
temperature to 300 K would give `1.66e-19 m^2`; neither behavior is stated by
Rubin et al. and must not be selected implicitly. The run is preserved as
`PT.failed-crifo-temperature-range-20261028`, with its complete log in
`halley1d.failed-crifo-temperature-range-20261028.log`. The selected repair
caps the collision-law input at the record's published upper endpoint and
uses the 180 K floor when the completed-window mean model-particle occupancy
is below `CollisionTemperatureMinimumSampleParticleNumber`. The numerical
Halley input sets the occupancy cutoff to `1.0`, equivalent to roughly 500
particle observations during a 500-iteration window before accounting for
temporal correlation.

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

A four-rank diagnostic of the first repaired binary reached the collision
callback at iteration 1 and stopped with `non-positive translational
temperature for comet collision species 13`. Core inspection established
that `cDataCenterNode::GetTranslationalTemperature()` reads only
`completedCellSampleDataPointerOffset` and returns zero while
`PIC::LastSampleLength` is zero. The requested rule is now explicit case data:
`CollisionTemperatureFloor=180.0`. `CometPhysics` evaluates the Crifo law at
the completed sampled temperature when it exceeds 180 K and at 180 K
otherwise. A clean four-rank run subsequently crossed the first 500-step
sampling boundary and continued through iteration 525 without that failure.

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
| `ampsConfig.pl` | modified | Adds the opt-in `LocalTimeStepFluxResampling` main-input keyword and rewrites `_PIC_LOCAL_TIME_STEP_FLUX_RESAMPLING_MODE_` in the generated core. Merge it with the core default/mover/timestep changes; without the parser change the application cannot enable the compiled path. |
| `input/hartley2.pl` | modified | Extends the historical Hartley parser into the generic comet preprocessor. It parses dimension/geometry, mesh units, source records, radial controls, chemistry, neutral collisions and completion policy, completed-sample collision-temperature controls, and per-species cooling; selects DIM1 spherical core options; includes `LocalTimeStepValidation.o` in the radial object set; validates required values; writes species arrays; and generates chemistry/collision headers. It depends on the new Perl modules and exact initialized declarations in the model/runtime sources. Preserve legacy Hartley keywords and the validation object when resolving conflicts. |
| `input/hartley2.input` | modified | Adds the generic DIM3 case declarations, configurable NASTRAN mesh units, heliocentric distance, disabled chemistry, the three neutral collision laws that were formerly hard-coded in `Comet.cpp`, and the explicit 180 K completed-sample collision floor. The remaining historical Hartley 2 settings are intentionally retained. |
| `input/hartley2.radial.input` | new | Defines the mesh-only, gas-only, spherically symmetric Hartley 2/Shou 2025 DIM1 case. It contains the area-equivalent nucleus radius, 10,000 km altitude domain, 1,000-cell ideal logarithmic profile, and 1 m surface resolution estimate. It deliberately leaves `RadialSimulationMode=off`. |
| `input/halley1d.input` | new | Defines the gas-only Rubin 2011 Halley simulation: 27 species, 11 parent sources, 6 km inner radius, `10^6 km` outer radius, 2,000 ideal cells, derived surface resolution, species-local timesteps with species-global weights and face-flux resampling, the explicit Rubin--Combi--CHEMKIN/Lorentz 378-pair collision hierarchy, generated chemistry, 180 K collision-temperature floor, H2O Larsen-Borgnakke rotation/RT relaxation, H2O Gombosi cooling, and convergence controls. Dust, extended sources, vibrational relaxation, and the optional methanol-to-formaldehyde channel are disabled. |
| `input/species.input` | modified | Adds the neutral species required by the Rubin network (`CH`, `CH2`, `N`, `CH3`, `NH`, `NH2`, `NH3`, `C2`, `C2H`, `C2H2`, `CN`, `HCN`, `C2H4`, `CO`, `C2H5`, `HCO`, `C2H6`, `H2CO`, and `CH3OH`). Carbon uses exact `12*_AMU_`, consistent with the nominal integer masses of the network, so compile-time reaction mass conservation is not defeated by a rounded kilogram literal. Merge by component name and avoid duplicate definitions if `master` has independently added any species. |
| `input/CometReactionNetwork.pm` | new | Parses arbitrary stoichiometric input reactions, resolves species names to compiled indices, validates rate-law fields, and emits `CometChemistry.generated.h`. It is the input-side producer for `ReactionNetwork.h` and `CometChemistry.cpp`. |
| `input/CometPhysicsDatabase.pm` | new | Parses explicit collision records and supplies the named `RUBIN2011` data set. It writes 66 unordered parent-pair collision entries and adds the selected Rubin Table 4 chemistry channels to the reaction-network generator. Its input-selected completion preserves those rows, adds 18 nonoverlapping Combi values, and generates 294 labeled geometric hard-sphere estimates from CHEMKIN diameters and Lorentz mixing. It rejects duplicates and emits `CometPhysics.generated.h`. |
| `input/test_comet_reaction_network.pl` | new | Tests reaction parsing, all supported rate-law syntax, loss channels, generated output, and rejection of unknown species. |
| `input/test_comet_physics_database.pl` | new | Tests the Rubin parent count, full 378-pair hierarchy and its 66/18/294 provenance partition, duplicate-completion rejection, generated chemistry/collision data, the required OH-rate choice, methanol sensitivity handling, and custom collision records. |

### Reusable core and model-library files

| File | State | Why it changed and merge dependency |
| --- | --- | --- |
| `src/models/chemical_solver/ReactionNetwork.h` | new | Adds a small, header-only, validated reaction-data container with stoichiometric terms, constant/inverse-square/Arrhenius rate laws, mass-action evaluation, source/loss evaluation, and unit-unimolecular detection. It is independent of the incomplete legacy CHEMEQ2 driver. |
| `src/models/chemical_solver/test_reaction_network.cpp` | new | Unit-tests validation, rate evaluation, stoichiometry, source/loss construction, and unimolecular classification. |
| `src/models/chemical_solver/README.md` | new | Documents the supported network representation and the current limitation that general multi-reactant cell execution is not yet implemented. |
| `src/models/cometary_physics/CollisionDatabase.h` | new | Adds a header-only unordered species-pair database for constant and temperature-power collision cross sections, including explicit temperature validity ranges and distinct missing-pair/out-of-range results. |
| `src/models/cometary_physics/test_collision_database.cpp` | new | Unit-tests symmetric pair lookup, both collision laws, validation, missing pairs, and range enforcement. |
| `src/models/cometary_physics/RadiativeCooling.h` | new | Adds a geometry-independent, unit-documented evaluator for the Gombosi et al. (1986) analytic H2O rotational-cooling law. It returns a positive per-molecule loss in J/s and contains no AMPS particle or mesh code. |
| `src/models/cometary_physics/test_radiative_cooling.cpp` | new | Unit-tests both temperature branches, optical-depth attenuation, SI/CGS conversion, and zero/invalid state handling. |
| `src/models/cometary_physics/README.md` | new | Documents collision and cooling semantics, units, validation, completed-sample temperature behavior, the published daughter-pair subset, Lennard-Jones analytical-estimate option, and the rule against silent extrapolation or undocumented zero cross sections. |
| `src/meshAMR/meshAMRgeneric.h` | modified | Contains the largest core change. It separates the `findNeibNodeFace()` lookup method from the `neibNodeFace` storage member; implements DIM1 face lookup, spherical face area, neighbor refinement queries, dimension-specific ghost ranges, and symmetry checks; limits block allocation to `2*DIM` faces; fixes the maximum-refinement integer-width test; bounds maximum-refinement spacing scans by the compiled mesh dimension; and makes mesh diagnostics dimension-aware. All callers of the renamed method must be reconciled together. |
| `src/meshAMR/meshAMRcutcell.cpp` | modified | Bounds global-domain extent reads by `_MESH_DIMENSION_` while retaining three-component physical ray vectors. This removes optimized lower-dimensional out-of-bounds reads without changing DIM3 geometry. |
| `src/pic/parallel/init_send_halo.cpp` | modified | Updates the DIM1/DIM2/DIM3 halo constructors to call `findNeibNodeFace()`. This is an API-following change and must be merged with `meshAMRgeneric.h`. |
| `src/pic/pic.h` | modified | Uses the dimension-neutral `cAmpsMesh` alias in the tracker and declares the local-timestep remainder field/accessors, continuation-sweep API, opt-in flux/lifecycle diagnostic record, reset function, and validation mover override. These declarations must land with `pic_mover.cpp`, `pic_pbuffer.cpp`, `pic_time_step.cpp`, and the validation harness. |
| `src/pic/picGlobalDefinitions.dfn` | modified | Defines explicit built-in and user-defined model-particle cross-section modes. These constants are required by the collision adapter and must land with `picGlobal.dfn` and `pic_particle_collision_models.cpp`. |
| `src/pic/picGlobal.dfn` | modified | Selects the built-in molecular cross-section pointer by default and defines local-timestep flux resampling off by default. Application input can independently override the collision and timestep modes; the default-off transport switch preserves existing applications. |
| `src/pic/pic_particle_collision_models.cpp` | modified | Routes every full/simplified NTC and MF cross-section evaluation through a common adapter. Built-in HS/VHS/VSS behavior is retained; user-defined laws receive velocity, species, block, and cell. This fixes the previously compiled-but-unreachable Rubin callback. |
| `src/pic/pic_parallel.cpp` | modified | Uses the common `cAMRnodeID` alias in buffered and unbuffered particle exchange, initializes three-component work vectors, and reads node coordinates/visits faces only in active dimensions. |
| `src/pic/pic_block_send_mask.cpp` | modified | Restricts send-mask neighbor scans to the active `2*DIM` faces and the number of subface slots represented in that dimension. This removes out-of-bounds DIM1/DIM2 neighbor reads and is required for multi-rank radial runs. |
| `src/pic/pic_mesh.cpp` | modified | Makes all four block pack/unpack and boundary send/receive paths visit positive faces only for active coordinates. Sender and receiver loop order must remain identical; this fix is required for load redistribution after a DIM1 sampling window. |
| `src/pic/pic_mesh_search.cpp` | modified | Traverses only child axes represented by the compiled mesh, copies only active node-bound components, and fixes inactive hash indices at zero. This removes optimized DIM1/DIM2 reads past node coordinate arrays. |
| `src/pic/pic_mover.cpp` | modified | Retains the reduced-dimension fixes and adds face-resolved local-timestep transport: compile-time mode checks, dynamic remainder storage, `W/dt` residual resampling, old-to-new timestep scaling, first-/second-order face stopping, local/MPI continuation, and pending-work counting. It uses a stable positive spherical-root solver for inward/outward crossings, normalizes the DIM1 spherical state before cloning/MPI transfer, implements second-order radial midpoint/no-intersection lookups, and exposes diagnostics plus a validation-only mover override. Both movers now have serial and four-rank focused runtime coverage. |
| `src/pic/pic_pbuffer.cpp` | modified | Clears the optional continuation fraction on every reused particle-buffer slot so injection, chemistry, and normal allocation cannot inherit stale unfinished work. This is a required correctness dependency of the dynamic field in `pic_mover.cpp`. |
| `src/pic/pic_time_step.cpp` | modified | Transfers main-mover face stops, then repeats unfinished movement, unconditional post-move ownership exchange, and global pending counting until completion before collisions/chemistry/sampling proceed. The post-move exchange prevents zero-remainder particles from being stranded in remote boundary layers. It retains legacy ordering when the opt-in macro is off and must merge with the mover continuation API. |
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
  than a hard-coded `cMeshAMR3d` specialization. The mover namespace also
  declares the dynamically allocated remainder offset, typed accessors,
  continuation sweep, and pending-particle counter. The optional field changes
  particle layout only when the generated opt-in macro is on.
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
  selected face; and the cut-surface block-scale calculation reads only
  active node-bound coordinates. The new opt-in path reserves one `double` in
  the dynamic particle record, stops both trace movers exactly at block faces,
  applies residual split/roulette with mean
  `W_old*dt_new/(W_new*dt_old)`, stores the normalized unused interval, and
  resumes only marked particles after list commit or MPI exchange. It assigns
  new tracker IDs to clones. Before cloning, the DIM1 path rotates the physical
  position and velocity into the stored radial frame so discarded transverse
  coordinates cannot make a clone's velocity inconsistent after MPI transfer.
  A common stable positive-root helper is used for inward and outward spherical
  face intersections; it fixes the historical inward-root coefficient mutation.
  The second-order DIM1 midpoint and no-intersection node lookups are now
  implemented and have the same focused runtime matrix as the first-order path.
  Flux/lifecycle instrumentation and the mover override are inert unless the
  application validation mode explicitly enables them.
- `src/pic/pic_pbuffer.cpp`: all three particle-slot allocation routes clear
  the optional remainder field immediately after marking a recycled slot
  allocated. Deliberate clone/MPI copy follows that initialization and may
  restore a nonzero marker. This prevents stale bytes from scheduling an
  unrelated new particle for continuation.
- `src/pic/pic_time_step.cpp`: the opt-in execution path first exchanges face
  stops from the main mover. Each sweep then advances marked remainders,
  unconditionally exchanges the resulting boundary-layer particles, and only
  afterward globally counts pending work. The unconditional second exchange
  fixes completed continuations that have a zero marker but still reside on a
  remote boundary block. It includes a defensive one-million-sweep guard.
  Generic splitting and list validation occur only after the loop; the
  default-off branch keeps the historical splitting/check/exchange order.

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

At the time of this record, branch HEAD is `905133b`. The following core and
core-configuration files have additional working-tree changes for the
local-timestep implementation:

```text
ampsConfig.pl
src/pic/pic.h
src/pic/picGlobal.dfn
src/pic/pic_mover.cpp
src/pic/pic_pbuffer.cpp
src/pic/pic_time_step.cpp
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
6. **Local-timestep transport:** merge `ampsConfig.pl`, the continuation
   declarations in `pic.h`, the default-off macro in `picGlobal.dfn`, and
   `pic_mover.cpp`, `pic_pbuffer.cpp`, and `pic_time_step.cpp` as one batch.
   Then merge the Halley input, `RadialSimulation.cpp`, and
   `CometChemistry.cpp`. Do not enable the input switch if any core member of
   this batch is absent. Conversely, leaving the source-core default off makes
   the merged machinery inert for unrelated applications.

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

# The final Halley build must select the complete compatible mode quartet.
rg -n 'LOCAL_TIME_STEP_FLUX|SIMULATION_TIME_STEP_MODE|SIMULATION_PARTICLE_WEIGHT_MODE|INDIVIDUAL_PARTICLE_WEIGHT_MODE' \
  build/pic/picGlobal.dfn
```

After conflict resolution, remove `build` before every configuration compile,
as required for this application. Run the standalone reaction, collision, and
radial-profile tests; then perform clean DIM1 and DIM3 builds. The current
working tree passed clean optimized DIM1 and DIM3 compile/link regressions
after the sampled-temperature and radiative-cooling changes. A full DIM3
Hartley runtime still requires the mission SPICE kernels listed above.

The DIM1 MPI acceptance run must cross an output boundary and actual block
redistribution, because a one-rank or pre-output run cannot exercise
`pic_block_send_mask.cpp` and the pack/unpack changes. Separately verify that a
user-defined collision case reaches its application callback. The present
Halley diagnostic reaches the Rubin callback, uses the explicit 180 K floor
before a completed sample exists, and crosses the first 500-step output and
redistribution boundary.

Do not merge or commit `build/`, `PT*/`, `amps`, run logs, generated mesh and
signature files, root-level parser/input copies, or generated application
headers copied out of `build/main`. Regenerate those artifacts from the
tracked input and source files after the merge.

### Generic comet application files

| File | State | Why it changed and merge dependency |
| --- | --- | --- |
| `srcHartley2/makefile` | modified | Adds the generic model, chemistry, collision physics, radiative cooling, radial mesh/runtime, and local-timestep validation objects. During DIM1 preprocessing, `hartley2.pl` replaces the DIM3 object list and main with the radial set and `RadialMain.o`; keep that generated object list synchronized, including `LocalTimeStepValidation.o`. |
| `srcHartley2/main.cpp` | modified | Replaces the hard-coded `hartley2.shape.nas` kilometre reader with the mesh filename and length unit from `CometModel`; validates generated chemistry, collision physics, and the selected cooling/IDF combination; and calls the generic cooling adapter between DIM3 iterations. |
| `srcHartley2/Comet.h` | modified | Removes declarations of the superseded Hartley-only Crovisier cooling routine and time-step wrapper. Cooling now has one geometry-neutral application interface. |
| `srcHartley2/Comet.cpp` | modified | Removes the hard-coded H2O/CO2 collision table and the duplicate Hartley-only water-cooling implementation. It delegates the legacy DIM3 AMPS collision callback to `Comet::Physics::GetNeutralCrossSection()`. DIM1 does not link this file, so `CometPhysics.cpp` supplies the callback only under `#if DIM == 1`; preserve that split to avoid either an undefined or duplicate symbol. |
| `srcHartley2/CometModel.h` | new | Declares geometry modes and per-species surface/extended source descriptors shared by all comet geometries. It is the stable interface between input-generated data and geometry-specific execution. |
| `srcHartley2/CometModel.cpp` | new | Stores generated model data and validates geometry, heliocentric distance, mesh units, source rates, source modes, temperatures, and bulk speeds before a run. Its initialized declarations are substitution targets for `hartley2.pl`; retain their spelling and one-line form. |
| `srcHartley2/CometChemistry.h` | new | Declares the generated-network adapter and the AMPS particle-reaction callback. It includes the reusable reaction header from either a generated or direct source-tree build. |
| `srcHartley2/CometChemistry.cpp` | new | Validates the generated network and executes supported unit-unimolecular events. It selects mutually exclusive channels, expands integer product multiplicity, handles explicit sinks, and assigns isotropic product velocities with center-of-mass momentum and excess-energy conservation. Correction-on builds carry the exact parent-to-product `W/dt` rate in one corrected daughter; correction-off species-global builds residual-resample by `(W_parent/dt_parent)/(W_product/dt_product)`, preserving source rates in expectation without hidden corrections. |
| `srcHartley2/CometChemistry.generated.h` | new/generated template | Provides an empty source-tree fallback so direct compilation is possible. The preprocessor writes the case-specific version into `build/main`; never replace this template with a generated Halley table during conflict resolution. |
| `srcHartley2/CometPhysics.h` | new | Declares the collision-database adapter, input-defined collision-temperature floor and minimum completed-sample occupancy, and neutral-collision callback used by both DIM1 and DIM3. |
| `srcHartley2/CometPhysics.cpp` | new | Builds and validates the generated collision database, enforces missing-pair failures, evaluates temperature laws from the completed AMPS translational-temperature sample, sends input-defined undersampled cells to the floor, caps high temperatures at each law's cited maximum, and owns the AMPS callback only in DIM1. Both initialized collision-temperature control declarations are input-substitution targets. |
| `srcHartley2/CometPhysics.generated.h` | new/generated template | Provides the empty `NONE` database for direct source-tree builds. As with the chemistry template, case data is generated into `build/main` and must not be copied back from the build directory. |
| `srcHartley2/CometRadiativeCooling.h` | new | Declares the per-species input model array, startup validation, and geometry-neutral between-iteration cooling hook. |
| `srcHartley2/CometRadiativeCooling.cpp` | new | Validates that cooling uses the AMPS Larsen-Borgnakke IDF model, reads completed density/rotational-temperature samples, evaluates the reusable selected law once per cell/species, and removes bounded rotational energy from locally owned particles with their species-local timestep. The initialized `SpeciesModel` array is an input-substitution target. |
| `srcHartley2/RadialMesh.h` | new | Declares explicit radial domain controls and the ideal logarithmic resolution profile used by the DIM1 mesh callback. Exactly one of radial extent above the surface and absolute outer radius is allowed. |
| `srcHartley2/RadialMeshProfile.cpp` | new | Solves the geometric cell-width ratio, maps cometocentric radius to ideal cell width, exposes ideal boundaries, and rejects invalid or decreasing-resolution requests. It contains no AMPS runtime dependency and is unit-testable by itself. |
| `srcHartley2/RadialMesh.cpp` | new | Connects parsed controls to AMPS, translates the ideal profile into a local-resolution callback, counts allocated dyadic AMR cells, and prints mesh diagnostics. |
| `srcHartley2/test_radial_mesh.cpp` | new | Tests expanding, uniform, one-cell, and invalid radial-profile configurations. |
| `srcHartley2/RadialMain.cpp` | new | Supplies the DIM1 executable entry point. It initializes AMPS, validates model/chemistry/collision/cooling/runtime data, builds or reads the distributed radial mesh, creates load lists, and starts the particle simulation only when enabled. It also parses `--validate-local-timestep`, particle-count, mover-selection, and multi-face controls and dispatches the focused test after normal mesh initialization. |
| `srcHartley2/LocalTimeStepValidation.h` | new | Defines the small validation-option interface kept separate from the production simulation API. |
| `srcHartley2/LocalTimeStepValidation.cpp` | new | Injects controlled ensembles at same-owner and MPI block faces; runs first-/second-order timestep-ratio matrices in both directions; checks list ownership, kinematics, remainder/lifecycle state, and number/mass/momentum/energy flux; and follows a long first-order trajectory through multiple faces and ownership transfers. |
| `srcHartley2/RadialSimulation.h` | new | Declares all required numerical controls and the AMPS callbacks for local timestep, injection, and domain-boundary processing. |
| `srcHartley2/RadialSimulation.cpp` | new | Implements spherical parent injection, AMPS internal-energy initialization at injection/re-emission, species- and radius-local timesteps, supported local/global weight initialization, thermal inner re-emission, outer loss, bounded daughter global-weight initialization, chemistry execution, between-iteration cooling, MPI ownership synchronization, full radial-profile assembly, relative-L1 convergence, sampling-window output, load redistribution, and the hard iteration ceiling. |

### Documentation and reproducibility tools

| File | State | Why it belongs in the merge |
| --- | --- | --- |
| `srcHartley2/README.md` | new | User-facing build/run guide for generic DIM1 and existing DIM3 cases, including generated-input behavior and validation limits. |
| `srcHartley2/GENERIC_COMET_DESIGN.md` | new | Records the intended separation among geometry, source data, chemistry, collision physics, and case inputs, including the reserved DIM2 axisymmetric interface. |
| `srcHartley2/STATUS.md` | modified | Records implemented behavior, test evidence, known physical omissions, paper-derived parameters, and this merge map. |
| `srcHartley2/comet_shape_properties.py` | new | Reproducibly reads the Hartley NASTRAN surface and derives area, volume, topology checks, and equivalent radii. |
| `srcHartley2/shou2025_spherical_parameters.py` | new | Reproduces the Hartley spherical source temperatures, densities, mean-free-path estimate, and mesh-related scalar parameters from the documented inputs. |
| `srcHartley2/halley2011/PAPER_PARAMETERS.md` | new | Separates values stated by Rubin et al. from AMPS numerical choices and records equations and unresolved physics needed for a full paper reproduction. |
| `srcHartley2/halley2011/COLLISION_DATA_REVIEW.md` | new | Records the online collision-data search for all 27 Halley species, the 84 directly published cometary pair values, CHEMKIN Lennard-Jones parameters covering every species, derived geometric estimates, provenance limits, and the selected input-level completion of the remaining 294 pairs. |
| `srcHartley2/halley2011/derive_halley1d_controls.py` | new | Reproduces the Halley surface-resolution estimate, particle-count choice, daughter characteristic speeds from the selected reactions and excess energies, and the post-failure Courant-0.25 timestep control. |
| `srcHartley2/halley2011/digitize_rubin2011.py` | new | Reproducibly extracts the neutral DSMC vector paths in Rubin Figs. 1--17 and modeled H2O velocity in Fig. 18, including the paper's split logarithmic distance axes. It verifies the exact input PDF digest before extraction. |
| `srcHartley2/halley2011/reference/rubin2011_neutral_profiles.csv` | new | Stores the extracted Rubin neutral density and velocity coordinates used by the overlay. This is reference plot data, not an original paper simulation table. |
| `srcHartley2/halley2011/reference/rubin2011_neutral_profiles.metadata.json` | new | Records the source-PDF digest, pages, axes, colors, variants, and vector-extraction provenance needed to audit or regenerate the CSV. |
| `srcHartley2/visualization/plot_1d.py` | new | Plots named quantities from one or many DIM1 Tecplot files with aliases, radial filtering, explicit Y scaling, log-Y masking, filename legends, PNG/EPS output, and optional viewer launch. |
| `srcHartley2/visualization/compare_rubin2011.py` | new | Selects a complete Halley output index and overlays AMPS profiles with the committed Rubin neutral DSMC vector digitization. It creates PNG/EPS panels, H2O radial velocity, a combined PDF, paired radius-sampled CSV values, and a provenance manifest. It does not synthesize the paper's NMS instrument response. |
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
   `CometChemistry`, `CometPhysics`, and `CometRadiativeCooling` adapters.
   Verify that preprocessing writes generated headers only into `build/main`.
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
  cooling-model writes, generated-header writes, and DIM1 object selection.
  The variable names used by `ChangeValueOfVariable` and
  `ChangeValueOfArray` must continue to match the initialized declarations in
  the application `.cpp` files exactly.
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
- Keep `RadiativeCooling.h` in `src/models/cometary_physics` independent of
  AMPS. Particle traversal and IDF access belong in
  `CometRadiativeCooling.cpp`. Preserve calls after `PIC::TimeStep()` in both
  mains so cooling runs between iterations and reads a completed sample.
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
The minimum acceptance evidence is a clean DIM1 build, the one- and four-rank
focused local-timestep harness, a clean DIM3 build, all standalone parser/
database/profile/cooling tests, a multi-rank DIM1 production run that crosses
a sampling/output/load-redistribution boundary, and `git diff --check`. A
successful compile or focused face test alone does not exercise the two MPI
redistribution defects fixed in `pic_block_send_mask.cpp` and `pic_mesh.cpp`.

## Post-test merge into `master`

Do not merge the comet branch directly into `master` merely because it
compiles or passes the focused transport harness. The local-timestep model
must first satisfy the complete acceptance gate below on `1d_comet_dev`. The
focused tests do not demonstrate bounded particle population over the full
Halley domain or a converged Rubin comparison.

### Merge disposition after the 2026-11-01 validation

**The full `1d_comet_dev` change set is not eligible for merge into `master`
yet.** The completed long run validates the face transport, chemistry
statistics, and scheduled redistribution implementation, but it fails the two
science-level gates that motivated the run: bounded population and a converged
Rubin comparison. At iteration 2696 the population was still increasing and
had reached 55,330,093 particles; the attempted iteration-2700 redistribution
filled rank 10's five-million-slot buffer. The last complete profile change at
iteration 2500 was 1.000366 rather than the required 0.05, with 0/3 consecutive
converged windows. Treating that run as acceptance would turn a documented
scaling failure into a `master` regression.

The validation state of the gate below is:

| Gate | Result | Evidence and consequence |
| --- | --- | --- |
| Clean optimized DIM1 configuration/build | **PASS** | A removed `build` tree was regenerated from `halley1d.input`; the executable linked with species-local timesteps, species-global weights, individual correction off, and face resampling on. |
| Same-owner/MPI face transport | **PASS** | The retained four-rank log passes all 32 single-face cases and the 95-face/two-ownership-transition continuation trajectory for both movers where applicable. |
| Unequal-weight/timestep chemistry | **PASS** | The retained 18-rank ensemble log passes all 46 neutral product reactions at 200,000 events per reaction; 23 loss channels are classified as sinks. |
| Scheduled block redistribution | **PASS for the original defect; insufficient for final stability** | The run passed the old iteration-2344/rank-17 failure and reached iteration 2699. Fifty-step redistribution repeatedly restored balanced measured work, but whole-block migration ultimately could not fit the growing population into equal fixed rank buffers. |
| Bounded production population | **FAIL** | Complete-output population grew from 49,983 at iteration 500 to 26,220,686 at iteration 2500, then to 55,330,093 in the iteration-2696 core diagnostic. Direct parents were 99.14% of the last census, so this is not a chemistry clone cascade. |
| Convergence termination | **FAIL** | Window changes were reference, 0.8129400, 1.000000, 1.005519, and 1.000366; no comparison window met 5%, so the required count remained 0/3. |
| Rubin comparison | **FAIL as acceptance; PASS as transient diagnostic tooling** | Output 4 produced the complete PNG/EPS/PDF/CSV/manifest artifact set, but the outward-increasing density deficit and manifest convergence record correctly identify it as non-converged. |
| DIM3/legacy regression | **PARTIAL** | Clean DIM3 compile/link passed earlier. A complete mission runtime still depends on the listed SPICE data, and an unrelated AMPS application with the new mode off must be rerun on the eventual integration commit. |

The local branch topology also matters. At the time of this note, local
`master` and the merge base are
`cb7be2bba41f373352ab853abed23955ea69c8b9`, while branch `HEAD` is
`905133b9ddcbfe43e987ec69712751c75d15051b`. The validated executable includes
working-tree changes newer than that `HEAD`; therefore the run cannot be
identified by the branch commit alone. Before any future integration, commit
the intended source on `1d_comet_dev`, make the tree clean, rebuild from a
removed `build` directory, and repeat the decisive gates on that exact commit.
Record `git rev-parse HEAD` in the new run archive.

### What must be completed before the merge is reconsidered

1. Replace the species-global statistical representation over the logarithmic
   domain with block-local weights, a conservative face reservoir, or an
   equivalent weighted-particle method. Preserve number, mass, all momentum
   components, and kinetic-energy flux in expectation across both local and
   MPI faces. A larger particle buffer or a hard clone cap is not an accepted
   fix.
2. Extend the focused harness to cover the selected new-weight algorithm in
   both radial directions, timestep and weight ratios below/equal/above one,
   both movers, a same-owner face, an MPI face, and the existing long
   multi-owner trajectory. Retain a complete PASS log.
3. Repeat the 46-reaction ensemble with the final block-local weights and
   timesteps. This is required because chemistry currently resamples the
   parent/product `W/dt` ratio and must follow the same weight policy as the
   mover.
4. Run Halley long enough to demonstrate a bounded or statistically stationary
   population without buffer exhaustion. Record global and maximum-rank
   populations at every sampling window and redistribution. The maximum-rank
   occupancy must retain an explicit safety margin; successful migration at
   one instant is not enough.
5. Reach three consecutive complete windows with every species' radial-density
   relative-L1 change at or below 5%, or define and document a different
   scientifically justified criterion before the run. The 10,000-iteration
   ceiling remains a safety limit, not evidence of convergence.
6. Regenerate the Rubin overlays from the final converged output. Preserve the
   input, exact commit, compiler/configuration record, MPI rank count, wall
   time, log, convergence history, population history, comparison CSV,
   manifest, and PDF together outside the source commits.
7. From that same clean commit, repeat the DIM1 and DIM3 clean builds, the
   parser/database/mesh/cooling tests, a DIM3 runtime when its kernels are
   available, and an unrelated AMPS application with
   `LocalTimeStepFluxResampling=off`.

Only after all seven items pass should the full comet/local-timestep feature
move to an integration branch. Until then, keep `1d_comet_dev`, its validation
logs, and its comparison artifacts intact and do not merge, rebase, squash, or
force-update it in a way that loses the failure provenance.

### Optional separation of independently reviewable work

If reusable work is needed in `master` before the population algorithm is
finished, prepare narrowly scoped integration commits rather than merging the
branch wholesale:

1. `ReactionNetwork.h`, `CollisionDatabase.h`, `RadiativeCooling.h`, their
   standalone tests, Perl database modules/tests, and documentation can be
   reviewed as libraries because they do not depend on the production
   local-timestep switch.
2. Reduced-dimensional AMR/PIC corrections can be reviewed separately, but
   only after DIM1/DIM2/DIM3 compile tests and MPI pack/unpack tests pass on
   the target `master`. Do not mix these fixes with the still-blocked Halley
   numerical policy.
3. Collision callback dispatch may be integrated as its own default-preserving
   batch only with a callback-reaching test; keep the built-in default for
   unrelated applications.
4. Do **not** integrate the global-weight face-resampling production setting,
   `halley1d.input` enablement, or claim a validated Rubin case until the
   population blocker above is closed. Off-by-default infrastructure could be
   proposed separately only if its scope, tests, and lack of production
   enablement are explicit.

Each optional subset needs its own clean branch from the current `master` and
must be selected by explicit paths or commits. Do not cherry-pick the large
working-tree diff, use `git add .`, or use the successful focused tests to
imply that the excluded full model passed.

### Acceptance gate before integration

1. Remove `build`, preprocess `input/halley1d.input`, and complete a clean
   optimized DIM1 link. Verify the generated timestep, final statistical-
   weight, individual-correction, and face-resampling modes as one compatible
   quartet. The already-tested species-global/correction-off quartet is a
   recorded failure over the full Halley domain and must not be presented as
   the final production selection.
2. Run one and multiple MPI ranks through at least one 500-iteration sampling
   boundary and the subsequent load redistribution. Record wall time, rank
   count, maximum particle-buffer occupancy, continuation-sweep count, and
   whether any remainder marker survives an outer timestep.
3. Run the implemented focused statistical harness with timestep ratios below,
   equal to, and above one. Over many crossings, measured incoming
   `N*W/dt`, mass, momentum, and kinetic-energy flux must agree with outgoing
   values within the documented 2%/five-sigma tolerance. Exercise crossings
   in both radial directions, both trace movers, an MPI ownership face, and a
   path spanning at least two MPI ownership boundaries. Preserve the command
   and complete PASS output with the eventual run archive.
4. Exercise every supported reaction using unequal parent/product effective
   `W/dt` ratios from the final local-weight policy. Check stoichiometric
   number, mass, momentum, and excess-energy sources over an ensemble. Exact
   eventwise equality is not expected when corrections are off; unbiased
   ensemble agreement is required.
5. Monitor the full logarithmic-domain population. The failed global-weight
   run establishes that growth to buffer capacity is not acceptable. The
   replacement must demonstrate bounded/statistically stationary global and
   maximum-rank populations with migration headroom. Never hide the problem
   with a clone cap.
6. Produce a converged Halley run and regenerate the Rubin overlay/table.
   Preserve input, executable git hash, log, rank count, wall time,
   convergence and population histories, and comparison artifacts together.
7. Rebuild the DIM3 Hartley input from a removed `build` directory and run its
   available initialization/runtime regression. Build at least one unrelated
   AMPS application with the new mode left off to prove legacy execution is
   unchanged.
8. Run every standalone chemistry, collision, cooling, mesh, parser, and
   visualization check listed under **Verification**, followed by
   `git diff --check`.

### Commit the tested development branch

Keep generated and unrelated files out of the commits. The existing local
`Makefile` linker edit is not part of the local-timestep implementation and
must be reviewed/staged separately. A reviewable commit split is:

```sh
git switch 1d_comet_dev

# Core configuration and transport are one dependency batch.
git add ampsConfig.pl \
  src/pic/picGlobal.dfn src/pic/pic.h src/pic/pic_mover.cpp \
  src/pic/pic_pbuffer.cpp src/pic/pic_time_step.cpp
git diff --cached --check
git diff --cached --stat
git commit -m "feat(pic): conserve flux across local-timestep block faces"

# The case and chemistry policy depend on that core batch.
git add input/halley1d.input input/hartley2.pl \
  srcHartley2/makefile srcHartley2/RadialMain.cpp \
  srcHartley2/LocalTimeStepValidation.h \
  srcHartley2/LocalTimeStepValidation.cpp \
  srcHartley2/RadialSimulation.h srcHartley2/RadialSimulation.cpp \
  srcHartley2/CometChemistry.h srcHartley2/CometChemistry.cpp
git diff --cached --check
git commit -m "feat(comet): add local-timestep validation and radial balancing"

# Documentation can record the final test paths and numbers last.
git add srcHartley2/README.md srcHartley2/STATUS.md
git diff --cached --check
git commit -m "docs(comet): document local-timestep validation and merge plan"

git status --short
git log --oneline master..1d_comet_dev
```

Before each commit, inspect `git diff --cached --name-status`; do not stage
`build/`, `amps`, root parser copies, mesh files, `PT*`, logs, `output/`,
`gitinfo.txt`, or generated links/directories. Do not use a broad `git add .`
in this worktree.

### Integrate without risking `master`

After the development branch is clean and the acceptance evidence is recorded,
resolve conflicts on a disposable integration branch created from the latest
`master`:

The commands below are future instructions only. In this workspace, stop and
report the completed gate evidence to the user before the first branch switch,
integration-branch creation, cherry-pick, or merge. Do not execute any of
those operations without the user's explicit approval.

```sh
git switch master
# If this clone tracks a remote and network policy permits it:
git pull --ff-only
git switch -c integrate/comet-local-timestep
git merge --no-ff 1d_comet_dev
```

Resolve conflicts by dependency batch, not by choosing an entire version of a
core file. Follow **Core merge procedure** and the conflict notes above. Then
remove `build` and rerun the entire acceptance matrix on the integration
branch. The generated root files must come from preprocessing on that branch;
never resolve a conflict by copying old generated files into source control.

Only after the integration branch passes should `master` move:

```sh
git switch master
git merge --ff-only integrate/comet-local-timestep
```

The fast-forward is intentional: it proves `master` has not moved since the
integration branch was based. If it fails, do not force or reset `master`;
create a new integration branch from the new tip, merge again, resolve, and
repeat testing. Keep `1d_comet_dev` and the validated run archive until the
post-merge DIM1 and DIM3 smoke tests pass. The feature can be disabled for an
emergency regression by setting `LocalTimeStepFluxResampling=off` and cleanly
rebuilding, because the source-core default is also off.

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
- `PhysicsDatabase` and its database-specific sensitivity choices;
- `CollisionTemperatureFloor` and per-species `RadiativeCooling` records; and
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
center-of-mass momentum and tabulated excess kinetic energy. Product
representation preserves the parent-to-product physical rate by converting
`W_parent/dt_parent` to `W_product/dt_product`; the correction-off Halley mode
uses residual split/roulette. Startup validation rejects noninteger product
multiplicities, nonconserving neutral channels, and multi-reactant networks.

### Generic collision physics

`src/models/cometary_physics/CollisionDatabase.h` provides validated,
unordered species-pair lookup for constant and temperature-power-law cross
sections. Each temperature law contains its cited validity interval. A missing
pair and an out-of-range temperature are separate low-level results. The
application adapter aborts on a missing pair. For temperature laws it applies
the input floor, sends cells below the configured completed-window mean model-
particle occupancy to that floor, and caps a hotter value at the law record's
maximum instead of extrapolating.

`input/CometPhysicsDatabase.pm` adds the `RUBIN2011` named database. It
transcribes the 66 parent-parent collision pairs in Table 3 and the neutral
photodissociation/photoionization channels in Table 4. The Crifo water
self-collision law is restricted to its published 50--300 K interval. Table 3
does not provide daughter-pair values; the named database alone leaves them
absent.

`CometPhysics.cpp` retrieves temperature-dependent collision state through
`cDataCenterNode::GetTranslationalTemperature()`. Core inspection confirms
that this accessor averages the completed sampling buffer and returns zero
when `PIC::LastSampleLength==0`. The collision adapter uses the input-defined
floor whenever the completed value is not warmer. Halley supplies 180 K; no
generic default is embedded in C++.

The collision-data search found Combi's complete H2O/OH/H2/H/O/CO subset,
also partially reproduced by Huang et al. (2016). Rubin and Combi together
cover 84 of the 378 unordered pairs in the 27-species network. The CHEMKIN
transport database supplies Lennard-Jones diameter and well-depth parameters
for every compiled species, including C2. The Halley input now explicitly
selects `COMBI1996_CHEMKIN_LORENTZ_HS`. The generator preserves 66 Rubin rows,
adds 18 nonoverlapping Combi rows, and applies a Lorentz mixed diameter with
`pi*d_ij^2` to only the remaining 294 pairs. Every estimated row is labeled in
its source string. `halley2011/COLLISION_DATA_REVIEW.md` records the values,
source links, provenance limits, and temperature-dependent alternatives.

### Internal energy and radiative cooling

`src/models/cometary_physics/RadiativeCooling.h` implements the
Gombosi et al. (1986) equation (48b) approximation to Crovisier/Huebner H2O
rotational cooling. The evaluator accepts K, m^-3, and m and returns a
positive per-molecule loss in J/s. Its unit test covers the two branches at
52 K, infrared attenuation, conversion, and invalid state.

`CometRadiativeCooling.cpp` maps the input-selected model for each species to
that reusable evaluator. It requires the AMPS Larsen-Borgnakke internal-energy
model, reads completed H2O number-density and rotational-temperature samples,
and removes `rate*local_timestep` from each locally owned H2O particle after
every particle iteration, clamped at zero. The legacy duplicate routine in
`Comet.cpp` was removed. Both DIM1 and DIM3 mains call the new adapter.

`halley1d.input` configures H2O with three rotational degrees of freedom,
`Zrot=5`, RT relaxation, and hard-sphere temperature index 0.5 against every
compiled partner. Radial surface injection and thermal inner re-emission call
the core IDF initializers. VT and VV remain off because the current core VT
transition explicitly exits as untested and unverified.

### Radial particle execution

`RadialSimulation.cpp` completes the generic DIM1 execution path when enabled:

- uniform parent rates are injected once through the spherical inner face;
- a Poisson process samples fractional-step injection times;
- species-local timesteps use an explicit characteristic speed and Courant
  number;
- the current Halley case uses one block-independent weight per species and
  performs unbiased `W/dt` split/roulette at each local-timestep block face;
- returning gas is fully thermalized and re-emitted at the inner boundary;
- particles crossing the outer boundary are deleted; and
- daughter species receive bounded global weights and chemistry converts
  parent to product `W/dt` rate statistically when individual corrections are
  off;
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
refinement level 29. `ReferenceInjectionParticleNumber=8` gives an estimated
minimum of 52.454 direct-parent particles in the surface cell. The Courant
number is 0.25 after the documented mover failure at the initially supplied
0.5 value. Direct parents start from the supplied 400 m/s characteristic
speed; daughter recoil bounds are propagated through the reaction graph from
Table 4 excess energies by `halley2011/derive_halley1d_controls.py`.

The methanol-to-formaldehyde branch is disabled. The particle buffer begins at
five million slots, consistent with the paper's statement that the full model
uses several million representative particles, and must be monitored in a
convergence run. The earlier species-local-weight configuration kept particle
count approximately uniform per radial cell. The newly requested
species-global-weight implementation instead conserves `W/dt` at faces by
cloning/roulette, so its steady-state population grows with radial timestep;
the former 2.21-million no-loss estimate does not apply to this new policy.

The clean one-process local-weight smoke run allocated 373 blocks and 2,984
radial cells; its smallest cell is `0.2328292 m`. Its approximate `2.21e6`
direct-parent estimate is retained only as historical evidence for that old
weight policy and must not be used to size the new global-weight run.

Rubin Table 4 gives the 0.36 eV OH branch rate as `[3e-6,6e-6] s^-1`; the
requested input selects its lower endpoint, `3.0e-6 s^-1`. Its other neutral
OH branches sum to `2.2e-6 s^-1`, and OH ionization is `1.2e-5 s^-1`. The
sampling window is 500 iterations and the hard ceiling is 10,000 iterations.
The run stops earlier when every species' shell-volume-weighted relative L1
radial density-profile change is at most 5% for three consecutive completed
window comparisons.

The paper includes Crovisier/Huebner radiative cooling and collisional
relaxation. Both are now connected through the generic application path. The
collision matrix is complete under the explicitly selected analytical
completion, whose 294 estimated rows must be identified in comparisons.

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
several mover/vector fixes. The production radial input still selects the
first-order spherical mover as its established policy. The second-order DIM1
trace path is now implemented sufficiently for the same focused local-
timestep face, continuation, and MPI regression matrix; changing the
production case to second order would remain a separate input/model decision.

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
  -Isrc/models/cometary_physics \
  src/models/cometary_physics/test_radiative_cooling.cpp \
  -o /tmp/test_radiative_cooling
/tmp/test_radiative_cooling

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
successfully. Generated source contains the 180 K collision-temperature floor,
assigns cooling model 1 only to H2O, enables Larsen-Borgnakke IDF and RT
relaxation, gives H2O three rotational modes and `Zrot=5`, and selects the
user-defined collision callback. After the collision-completion change, a
fresh `-O3` build generated database
`RUBIN2011+COMBI1996+CHEMKIN_LORENTZ_HS` with 378 rows: 66 source strings from
Rubin Table 3, 18 additional Combi rows, and 294 CHEMKIN/Lorentz estimates.
The optimized executable was launched from
`/tmp/amps_halley_collision_smoke_o3`, built the 373-block/2,984-cell mesh, and
ran through iteration 124 without a missing-pair or invalid-cross-section
abort before being deliberately stopped. This isolated smoke test did not
modify the repository `PT` directory and did not reach a complete 500-step
sampling window.

A separate clean optimized DIM3
Hartley compile and link also passed with the same generic collision and
cooling adapters before the tree was removed and rebuilt from
`halley1d.input`. The new four-rank diagnostic completed
the first 500-iteration sample with 77,482 particles, wrote all 27 outputs,
redistributed blocks, executed 26 subsequent between-iteration cooling calls
with the completed sample available, and was deliberately stopped at
iteration 525. That output predates the selected collision completion and
remains transient, so it cannot validate the newly completed model. The
earlier multi-window output remains invalid for physical comparison because it
used the zero placeholder cross section. The MPI compiler reports local socket
warnings in the restricted environment; they do not prevent compilation or
runtime. The optimized core also retains warnings in inactive/general-purpose
paths, including a `pic_fluid.cpp` aggressive-loop diagnostic. DIM3 retains
legacy warnings for `register`, format widths, narrowing conversions, and
fixed-buffer `sprintf` calls.

### 2026-11-01 long-run population limit and scheduled redistribution

The first corrected 18-rank production run is preserved under
`validation_runs/halley1d_lts_longrun_20261101`. It used the optimized DIM1
Halley executable, a 500-iteration sampling window, 10,000-iteration ceiling,
5%/three-window convergence criterion, global species weights, and
species-dependent block-local timesteps. It completed four full output sets
(iterations 500, 1000, 1500, and 2000). Sampled populations were 49,542,
164,020, 593,164, and 3,225,917. The corresponding maximum relative-L1
changes were reference, 0.8083251 (NH3), 0.7710021 (NH3), and 1.0 (H first
appeared as a sampled chemistry daughter). The production-helper ensemble
test immediately before the run passed all 46 neutral reaction channels at
200,000 events/channel on 18 ranks; all 46 exercised unequal weight/timestep
resampling, while 23 ionization-loss channels were classified as sinks.

The production run aborted at iteration 2344 when MPI rank 17 filled its local
5,000,000-particle buffer. The last complete core diagnostic before the abort
reported 11,636,990 global particles at iteration 2333, so aggregate capacity
was not the problem. Core review found that the sampling-break path moved
blocks only every 500 iterations. The emergency path uses accumulated
synchronization latency; that value stayed near zero because each outer
collision/chemistry step was already synchronized, even while the radial
population front accumulated on one owner. Increasing the buffer alone would
postpone, not correct, this ownership imbalance.

`RadialSimulation.{h,cpp}` now exposes and validates
`RadialLoadBalanceInterval`. After transport continuation, particle exchange,
collisions, chemistry, sampling, and between-iteration radiative cooling are
complete, the driver invokes the same core redistribution sequence used at
sampling breaks, resets the core rebalancing counters, and refreshes
`DomainBlockDecomposition::BlockTable`. It skips a duplicate call when an
output iteration already redistributed the mesh. Both parser copies require
the new positive integer for an enabled radial run, and `halley1d.input` sets
50 iterations. The particle buffer remains 5,000,000 per rank so the next
clean run directly tests whether redistribution prevents the observed
failure. No master-branch merge has been attempted.

### 2026-11-01 scheduled-balance long-run result

The corrected run is preserved under
`validation_runs/halley1d_lts_lb50_longrun_20261101`. It used the clean `-O3`
build generated from `input/halley1d.input`, 18 MPI ranks, the unchanged
5,000,000-particle buffer on each rank, 500-iteration sampling windows, a
10,000-iteration ceiling, and the 5%/three-consecutive-window convergence
rule. `/usr/bin/time -p` measured 6,165.38 s wall time (1 h 42 min 45.38 s).
The complete stdout/stderr record is `halley1d.long-run.log` in that
directory.

The focused production-path gates passed immediately before the run:

- the four-rank local-timestep test passed all 32 single-face cases and the
  95-face/two-ownership-transition continuation case with 200,000 particles
  per single-face case; and
- the 18-rank chemistry ensemble test passed all 46 neutral product reactions
  at 200,000 events per reaction. All 46 exercised unequal species weights
  and timesteps; 23 loss-only channels were classified separately as sinks.

The retained gate logs are
`validation_runs/local_timestep_lb50_20261101/local-timestep.validation.log`
and
`validation_runs/chemistry_ensemble_lb50_20261101/chemistry.validation.log`.

Scheduled redistribution corrected the failure that motivated it. The run
crossed iteration 2344, where the sampling-only run had filled rank 17 with
only about 12 million particles globally, and continued through iteration
2699. Redistribution every 50 iterations repeatedly reduced measured rank
work from spreads as large as 0.80--1.34 times the mean to approximately
0.91--1.08 times the mean. At iteration 2696 the largest pre-redistribution
rank population was 4,218,772 rather than five million. The attempted
iteration-2700 redistribution then aborted while rank 10 was receiving its
new block set: `pic_pbuffer.cpp:615`, `The particle buffer is full`. The core
printed the newly balanced work distribution before particle migration
exhausted the receiving rank. Thus the original ownership-localization defect
is fixed, but whole-block redistribution cannot fit this growing global
population into fixed equal-capacity rank buffers indefinitely.

Population did not approach a steady plateau:

| completed iteration | population at complete output | maximum relative-L1 profile change | limiting species |
|---:|---:|---:|:---|
| 500 | 49,983 | reference | -- |
| 1000 | 167,800 | 0.8129400 | NH3 |
| 1500 | 637,080 | 1.000000 | NH2 |
| 2000 | 3,678,905 | 1.005519 | NH2 |
| 2500 | 26,220,686 | 1.000366 | H2 |

The independent core diagnostics continued from 26,479,150 particles at
iteration 2501 to 55,330,093 at iteration 2696. The last five complete
convergence records were therefore reference, 0/3, 0/3, 0/3, and 0/3; the run
did not satisfy even one 5% window. The last diagnostic species census shows
that direct parents account for 54,856,785 particles (99.14%) and chemistry
daughters for 473,308 (0.86%). Representative daughter counts were H=15,225,
H2=1,748, OH=130,387, NH2=25,284, and HCO=2,820. Chemistry was active and
statistically well behaved, but it was not the cause of the population
growth. The growth is the expected numerical cost of representing the still
advancing spherical parent flux with species-global weights while local
timesteps, and therefore face-resampling multiplicities, grow across the
logarithmic domain.

The last complete output (`out=4`, nominal iteration 2500) was processed with
`visualization/compare_rubin2011.py`. Comparison artifacts are under
`PT/plots/rubin2011_comparison` inside the run directory: 18 PNG/EPS overlays,
`rubin2011_neutral_profiles.pdf`, `comparison_profiles.csv`, and
`manifest.json`. The manifest records the exact output index and the failed
convergence record. The comparison is scientifically useful as a transient
front diagnostic, not as a paper validation. For example, AMPS/Rubin HCN is
0.0328 at 100 km, 0.00275 at 1,000 km, and 0.000164 at 10,000 km; H2O is
0.000767 at 10,000 km and 0.0000231 at 100,000 km. The outward-increasing
deficit agrees with the order-unity convergence change and demonstrates that
the population front had not filled the paper domain.

Acceptance is consequently split: conservative local-timestep transport,
multi-rank continuation, chemistry kinematics/statistics, and frequent core
block redistribution pass their focused validations; this global-weight
Halley configuration fails population stability, convergence, and a
steady-state Rubin comparison. Increasing the per-rank buffer would postpone
the same scaling limit and is not a physical fix. The next algorithmic stage
must use block-local statistical weights (or an equivalent weighted-particle
or face-reservoir formulation) while preserving conservative interface flux.
No master-branch merge or commit was attempted.

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
`halley2011/PAPER_PARAMETERS.md`. Traceable vector-extracted Rubin DSMC curves
and their metadata are now committed under `halley2011/reference/`, and the
comparison utility overlays them with each completed AMPS window. A physically
meaningful comparison still requires a converged restarted output and must
state that 294 pair cross sections are analytical estimates outside Rubin's
printed table. The NMS observables remain out of scope until an instrument
response is implemented.
