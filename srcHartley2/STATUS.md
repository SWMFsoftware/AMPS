# Generic comet / Hartley 2 status

Date: 2026-11-03
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
| Halley DIM1 convergence run | Completed the block-local-weight run through the 10,000-iteration hard limit on four MPI ranks. Population remained bounded at 628,490 particles and no buffer or mover failure occurred, but the all-species convergence gate failed: final maximum relative L1 change `1.719734` (CH), `0/3` consecutive windows. |
| Core collision cross-section dispatch | Fixed and clean-built; NTC/MF now route the selected user law with block/cell state instead of always using the zero HS placeholder |
| Halley collision-data coverage | Complete by explicit policy: 66 Rubin rows, 18 additional Combi rows, and 294 labeled CHEMKIN/Lorentz estimates |
| Completed-hierarchy clean build and smoke test | Pass; clean `-O3` DIM1 compile/link and isolated one-rank execution through iteration 124 |
| DIM1 result visualization | Pass for short aliases, single-/multi-file PNG/EPS output, and viewer control |
| Rubin paper-comparison artifact generator | Pass as tooling and transient output: all 18 panels overlay output 19 with vector-extracted Rubin DSMC data and include PNG/EPS plus a combined PDF/CSV/manifest. Acceptance fails because the source run did not converge. |
| Local-timestep/local-weight face-flux transport | Pass for a clean `-O3` DIM1 compile/link, 200,000-particle statistical face tests for both movers on local and MPI faces, a 95-face/two-ownership-crossing trajectory, all 46 neutral chemistry reactions, an integrated production smoke, and bounded full-domain population through iteration 10,000. The generated mode quartet is species-local timestep, species-local block weight, correction off, and face resampling on. The strict all-species convergence gate still fails. |
| Sampling-break block-weight adaptation | Implemented for each species using the instantaneous whole-block mean and requested `50--200` particles/cell range. The selected `W/dt` carrier constraint is also propagated by powers of two through face-connected empty neighbors so a first arrival cannot clone against a stale low weight. Clean build and focused serial/MPI validation passed; evidence is recorded in the dated sections below. |
| Emergency pre-exchange weight adaptation | Implemented with input threshold `0.75`; clean optimized build and integrated one-/four-/twelve-rank callback tests pass. A new long convergence run is still required. |
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
2,000 without either the earlier Crifo-temperature stop or mover stop. During
the 2026-11-01 process audit this legacy process was found still running from
an already replaced/deleted executable and was stopped as required before the
block-local-weight work. It ended at iteration 9,305 with a latest observed
census of 750,220 particles. Its last complete convergence output was
iteration 9,000: maximum relative L1 change `0.9631707` for CH and `0/3`
consecutive passes. It therefore did not converge and is not evidence for the
current implementation. The retained log is
`halley1d.courant-0.25.run.log`; its runtime data and Rubin-overlaid comparison
remain under `PT/plots` and `PT/plots/rubin2011_comparison`.

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
| `input/halley1d.input` | new | Defines the gas-only Rubin 2011 Halley simulation: 27 species, 11 parent sources, 6 km inner radius, `10^6 km` outer radius, 2,000 ideal cells, derived surface resolution, species-local timesteps with block-local species weights and combined timestep/weight face-flux resampling, the explicit Rubin--Combi--CHEMKIN/Lorentz 378-pair collision hierarchy, generated chemistry, 180 K collision-temperature floor, H2O Larsen-Borgnakke rotation/RT relaxation, H2O Gombosi cooling, and convergence controls. Dust, extended sources, vibrational relaxation, and the optional methanol-to-formaldehyde channel are disabled. |
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
- the current Halley case uses one weight per species and block, proportional
  to that block's species-local timestep, and performs unbiased full-ratio
  `W/dt` split/roulette at each block face;
- returning gas is fully thermalized and re-emitted at the inner boundary;
- particles crossing the outer boundary are deleted; and
- daughter species receive a bounded local `W/dt` carrier rate and chemistry
  converts parent to product `W/dt` statistically when individual corrections
  are off;
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

### 2026-11-01 enlarged-buffer convergence rerun

At the user's request, `input/halley1d.input` now reserves 20,000,000 particle
slots per MPI rank, four times the capacity of the iteration-2700 failure.
The value is based on the observed 4,218,772-particle maximum rank population
immediately before the failed redistribution and provides room for both the
resident population and temporary whole-block migration imbalance.  With 18
ranks this exposes up to 360 million model-particle slots; the host had about
233 GiB available before the clean rebuild and launch.

This change is deliberately isolated to `RequestedParticleBufferLength`.
Source rates, statistical weights, species-local timesteps, 50-iteration load
balancing, the 500-iteration sampling window, the 5%/three-window convergence
test, and the 10,000-iteration safety ceiling remain unchanged.  Therefore a
different convergence history can be attributed to stochastic evolution and
additional runtime rather than a changed physical model.  The rerun must be
built only after removing `build`, launched from a new validation directory,
and retained with its complete log and comparison manifest.  If population
growth fills the enlarged buffer before convergence, that outcome must be
reported as a population-stability failure; increasing capacity again is not
by itself a validated remedy.

The clean rebuild completed successfully after deleting `build`; generated
`build/pic/pic.cpp` was checked for the 20,000,000-slot buffer, 10,000 maximum
iterations, 500-iteration sampling length, 5% tolerance, and three required
consecutive windows.  The 18-rank production rerun was then started in
`validation_runs/halley1d_lts_lb50_buf20m_convergence_20261101`; its live and
retained log is `halley1d.convergence.log`, and its sampled data are written to
`PT/plots` below that directory.  The first four completed windows contained
51,576, 172,215, 660,541, and 3,874,714 particles.  Their convergence records
were reference, 0.8129400 (NH3), 1.000000 (NH2), and 1.005519 (NH2), so no
window through iteration 2000 met the 5% threshold and the run correctly
continued.  This paragraph is a launch/progress record, not a claim of
convergence; the final termination reason and replacement policy are recorded
in the following block-local-weight section.

For a later master merge, keep the enlarged buffer as a case-level resource
choice rather than a core default.  Do not merge this input edit as evidence
that species-global weights solve the logarithmic-domain scaling problem, and
do not change `src/pic` buffer defaults.  The local-weight validation retains
the value as temporary headroom; production population evidence must determine
the final case capacity.

### 2026-11-01 block-local particle-weight transport

The 20-million-slot global-weight rerun was stopped at the user's request
after completing iteration 2729.  The last core census before the interrupt
contained 62,087,575 particles; its last complete convergence output remained
iteration 2500 with maximum relative L1 change `1.000366` (H2) and `0/3`
converged windows.  Exit status 130 records the requested interrupt, not a
physics, MPI, or particle-buffer failure.  Its retained directory remains
`validation_runs/halley1d_lts_lb50_buf20m_convergence_20261101`.

The same final host audit also found the older four-rank Courant-0.25 process
described above.  It was not the newly built block-local-weight executable:
Linux reported its executable as deleted, and its five-million-slot startup
record predates the current input and build.  That process group was stopped
at iteration 9,305.  A follow-up process scan found no remaining Halley AMPS
launcher or worker, so no production simulation is currently active.

The replacement Halley configuration now combines
`SpeciesLocalTimeStep`, `SpeciesLocalParticleWeight`, correction off, and
`LocalTimeStepFluxResampling=on`.  For a species `s`, the radial initializer
sets the weight in block `b` to

```text
W[s,b] = (W[s,surface]/dt[s,surface]) * dt[s,b].
```

Thus each species retains its requested surface injection statistics while
`W/dt` is constant with radius and expected steady-flux occupancy is roughly
constant per logarithmic radial cell.  Non-source daughters use the largest
direct-parent `W/dt` carrier rate.  This reduces direct-parent-to-daughter
clone bursts without changing any physical source rate; reactions into a
direct species with a smaller carrier rate can still clone and are handled by
the same unbiased chemistry conversion.

At every crossed block face the mover reads weights and timesteps from both
blocks and draws the destination multiplicity with mean

```text
mu = (W_old/dt_old)/(W_new/dt_new)
   = W_old*dt_new/(W_new*dt_old).
```

The draw is `floor(mu)` plus a Bernoulli residual.  When `mu<1`, the deletion
probability is `1-mu`; when `mu>1`, the original survives and the additional
destination multiplicity is cloned.  Every sibling retains the same face
position and velocity, so number, mass, momentum, and kinetic-energy flux are
conserved in expectation.  The unused trajectory interval remains a separate
kinematic invariant, `tau_new=tau_old*dt_new/dt_old`; particle weight changes
the population probability, not the remaining flight fraction.

#### Files changed for block-local weights

| File | Modification and merge reason |
| --- | --- |
| `input/halley1d.input` | Selects `SpeciesLocalParticleWeight`, explains the complete `W/dt` face conversion, and retains the 20-million-slot validation buffer.  This case edit is required to activate the new path but must not be treated as a core default. |
| `src/pic/picGlobal.dfn` | Expands the default-off mode documentation to permit species-local block weights and records that both weight and timestep ratios enter multiplicity.  The macro remains off for other applications. |
| `src/pic/pic_mover.cpp` | Removes the former compile-time rejection of `_SPECIES_DEPENDENT_LOCAL_PARTICLE_WEIGHT_`.  The already implemented full `W_old*dt_new/(W_new*dt_old)` expression is now active for local weights; all positive/finite checks, diagnostics, continuation storage, and unsupported-mode guards remain in force. |
| `srcHartley2/RadialSimulation.cpp` | Replaces the surface-pair helper with explicit constant-`W/dt` initialization on every allocated block.  Direct parents preserve source-derived surface weights; daughters use the maximum direct-parent carrier rate.  Invalid local timesteps or derived weights abort at initialization. |
| `srcHartley2/CometChemistry.cpp` | Allows correction-off chemistry with species-local block weights, reads parent and product weights from the current block, and retains the general source/product `W/dt` residual resampling.  Validation comments now describe the block-local production configuration. |
| `srcHartley2/LocalTimeStepValidation.cpp`, `LocalTimeStepValidation.h` | Extend the focused regression to set source/destination weight ratios independently of timestep ratios.  The matrix covers roulette, cloning, reinforcement, exact integer multiplicity, and the production cancellation `W_new/W_old=dt_new/dt_old` on same-rank and MPI faces for both trace movers. |
| `src/pic/README.md` | Documents the core local-weight face contract, deletion/cloning probability, continuation invariant, block-data ownership, and unsupported configurations. |
| `srcHartley2/README.md` | Documents the Halley constant-`W/dt` policy, revised validation matrix, smoke evidence, and remaining production acceptance gate. |
| `srcHartley2/halley2011/PAPER_PARAMETERS.md` | Reconciles the case table with the 20-million-slot validation capacity and correction-off local-weight policy; neither value is attributed to Rubin et al. |
| `srcHartley2/STATUS.md` | Records the stopped global-weight run, implementation, evidence, limitations, and integration procedure in this section. |

No new core source file was needed.  In particular, the existing
`cDataBlockAMR::SetLocalParticleWeight()`/`GetLocalParticleWeight()` storage in
`src/pic/pic.cpp`, block-associated-data migration, and complete particle
record MPI packing were reused without modification this round.

#### Verification evidence

`build` was deleted before configuration.  The optimized clean build
completed successfully, and generated `build/pic/picGlobal.dfn` contains
species-dependent local timestep, species-dependent local particle weight,
face resampling on, and individual correction off.  Existing aggressive-loop
warnings in unrelated legacy DIM1 paths were unchanged; there was no new
compile or link warning attributable to this implementation.

The retained validation root is
`validation_runs/halley1d_local_weight_face_validation_20261101`:

- `serial/local-weight.serial.log`: one rank, 200,000 source particles per
  case, both trace movers, all 28 local-face direction/ratio cases passed;
- `mpi4/local-weight.mpi4.log`: four ranks, all local and MPI-face cases for
  both movers passed, followed by a 95-face path crossing two MPI ownership
  transitions with 32/32 final particles and no pending remainder;
- `chemistry_mpi4/local-weight.chemistry.mpi4.log`: all 46 neutral reactions
  passed 200,000-event number/mass/momentum/energy ensemble checks; 25 had
  unequal source/product carrier rates and 23 sink reactions were identified;
- `production_smoke/local-weight.production-smoke.log`: the full collision,
  chemistry, cooling, sampling, and 50-iteration load-balance path completed
  outputs at iterations 500 and 1000 before the validation smoke was
  deliberately interrupted after iteration 1046.  Populations were 40,841
  and 78,754, compared with 51,576
  and 172,215 in the superseded global-weight run.  The expected transient
  convergence record at iteration 1000 was `0.8107134`, NH3, `0/3`.

The production H2O `out=1` file directly confirms spatially local weight:
printed block weights span `1.273285e25` to `4.272435e32`, a factor of
`3.355443e7`, while `W/dt` spans only `5.273e-7` relatively at output
precision.  This is the intended cancellation, not evidence of a globally
constant weight.

#### Full-domain 10,000-iteration result

The clean production run is retained at
`validation_runs/halley1d_local_weight_convergence_rubin_20261101`. Before
compilation, `build` was deleted. The frozen executable, input, and generated
macro file are `amps.local-weight`, `halley1d.input`, and `picGlobal.dfn` in
that directory. Their SHA-256 digests are, respectively,
`e679021a0900b5a12e55d27e86814dee3d6f8973c5a8a4b5f440049430017c11`,
`9c40df1cd1c4235fac25be9efac03fffee84b28dc1a03f43e54c550b5c50d757`,
and `f4719d752c828cfa0303c8b14c1959bd73cc16aa710fcc4f0a68e77cdcea6c72`.
The generated configuration confirms species-local timestep,
species-local block weight, individual correction off, and face-flux
resampling on.

The command `mpirun -np 4 ./amps.local-weight` ran from 18:17:06 to 21:35:47
on 2026-11-01, approximately 3 h 18 min 41 s wall time. It completed all 20
500-iteration samples and wrote 540 species files (`27 species x 20
outputs`) under `PT/plots`. The final population was 628,490 particles; no
particle-buffer exhaustion, invalid face conversion, chemistry failure, or
MPI error occurred. This is a successful full-domain population-stability
test for block-local weights, not a successful steady-state test.

The maximum relative L1 profile changes at iterations 500 through 10,000
were: reference, `0.8107134`, `0.7817780`, `1.0`, `1.0`, `1.006929`,
`1.212851`, `1.458767`, `1.312500`, `1.824477`, `1.834117`, `1.421037`,
`1.527195`, `1.748311`, `1.776113`, `1.716863`, `1.904352`, `1.754535`,
`1.796966`, and `1.719734`. The final limiting species was CH and no window
met the required three consecutive all-species changes at or below 5%.
Late windows are controlled by statistically sparse CH/CH2/CH3 daughter
profiles. Therefore output 19 must remain labeled transient/non-converged,
even though the direct-parent population no longer exhibits the catastrophic
growth of the global-weight calculation.

The Rubin overlay generated from output 19 is in
`PT/plots/rubin2011_comparison` below the retained run. The combined file is
`rubin2011_neutral_profiles.pdf`; the same directory contains 18 PNG/EPS
panels, `comparison_profiles.csv`, and `manifest.json`. At the five exported
comparison radii, H2O radial velocity differs from the digitized Rubin curve
by approximately +8.8%, +3.7%, -0.07%, -5.6%, and -7.2% at 10, 100, 1,000,
10,000, and 100,000 km. Where both H2O density curves are available, AMPS is
0.766 times Rubin at 10,000 km and 0.911 times Rubin at 100,000 km. These
numbers are diagnostics only: the failed convergence record is embedded in
the manifest and prevents treating the apparent agreement as model
validation. CO remains strongly underpredicted in the shared digitized
range, which also requires physics/source investigation after a statistically
appropriate convergence criterion is established.

#### Remaining acceptance and later master merge

The focused algorithm, MPI continuation, chemistry conversion, integrated
run, full-domain population stability, and comparison tooling pass.
Three-window convergence and a scientifically accepted Rubin overlay do not;
do not merge to `master` or claim paper validation yet. The next validation
step is to separate direct-parent convergence from the noise floor of sparse
daughter species, while retaining explicit per-species diagnostics and a
strict physical acceptance rule. Any replacement criterion must be
documented and regression-tested before another production run; it must not
silently ignore unresolved chemically important daughters.

When that gate passes, integrate on a temporary branch based on current
`master`, never directly on `master`, and stop for user review before the
merge.  Bring the dependency set in this order: (1) core default/documentation
and mover support (`src/pic/picGlobal.dfn`, `src/pic/pic_mover.cpp`,
`src/pic/README.md`); (2) application initialization, chemistry, and focused
tests; (3) Halley input and application/paper/status documentation.  Preserve
the core macro default as off, preserve existing master behavior for global
weights, and resolve conflicts in `PrepareLocalTimeStepFaceTransfer()` by
retaining the complete `W_old*dt_new/(W_new*dt_old)` expression and all
continuation/MPI fixes.  After conflict resolution, delete `build`, rebuild
the exact Halley input, verify the four generated macros, rerun serial and
four-rank focused regressions plus chemistry, then repeat the converged
production comparison.  Do not copy generated `build/`, `amps`, `PT/`, mesh,
log, or validation-run artifacts into the source commit.

#### Master-integration runbook for the block-local-weight round

This runbook is for future use; it does not authorize a merge. The current
10,000-iteration run failed the convergence gate, so the procedure stops at
the development branch until a documented convergence policy and a steady
Rubin comparison pass. In particular, do not switch to `master`, create an
integration branch, cherry-pick, rebase, or run `git merge` without first
reporting the completed evidence and receiving explicit user approval.

First make the development branch reviewable. Commit source changes from
`1d_comet_dev` with explicit paths; never use `git add .` in this worktree.
For the block-local-weight round, the intended path set is:

```sh
git switch 1d_comet_dev
git add input/halley1d.input \
  src/pic/picGlobal.dfn src/pic/pic_mover.cpp src/pic/README.md \
  srcHartley2/RadialSimulation.cpp srcHartley2/CometChemistry.cpp \
  srcHartley2/LocalTimeStepValidation.cpp \
  srcHartley2/LocalTimeStepValidation.h \
  srcHartley2/README.md srcHartley2/STATUS.md \
  srcHartley2/halley2011/PAPER_PARAMETERS.md
git diff --cached --name-status
git diff --cached --check
git commit -m "feat(comet): add conservative block-local particle weights"
```

The local `Makefile` edit adding `-lmpi++` is not part of this feature and
must be justified and tested separately or left unstaged. Also exclude
`build/`, `amps`, `PT*`, `validation_runs/`, meshes, logs, generated root
parser/input copies, `output/`, `gitinfo.txt`, and temporary files. Before
integration, require a clean source status or a written explanation for every
remaining path, and audit the complete branch delta:

```sh
git status --short
git log --reverse --oneline master..1d_comet_dev
git diff --name-status master...1d_comet_dev
git diff --check master...1d_comet_dev
```

**Stop here and report the commit list, file list, validation results, and
known limitations to the user.** The following integration commands are
permitted only after explicit approval. Fetch the current upstream state when
available, create a disposable branch from that exact `master`, and perform a
non-committing merge so the result can be inspected before recording it:

```sh
git fetch origin
git switch --create integrate/comet-local-weight origin/master
git merge --no-ff --no-commit 1d_comet_dev
```

If no remote is available, create the integration branch from the verified
local `master` instead. Resolve conflicts by semantic dependency batch, never
by accepting an entire side of a core source file. Preserve these invariants:

1. Core defaults remain off for applications that do not select comet local
   transport; global-timestep and global-weight behavior must be unchanged.
2. `PrepareLocalTimeStepFaceTransfer()` retains
   `W_old*dt_new/(W_new*dt_old)`, finite/positive checks, stochastic residual
   multiplicity, continuation state, and local/MPI ownership handling.
3. Sender and receiver packing retain identical active-face and block-data
   ordering, including local timestep and local particle-weight storage.
4. Halley direct-parent weights retain their source-derived surface rates;
   block weights remain proportional to local timestep, and daughters use the
   documented conservative `W/dt` carrier policy.
5. Chemistry uses the current block's parent and product weights and performs
   unbiased `W/dt` resampling with individual correction disabled.
6. Do not resolve generated files. Regenerate them from the selected inputs
   after source conflicts are resolved.

Before committing the integration merge, remove `build`, compile the exact
Halley input, and verify the generated mode quartet: species-local timestep,
species-local block weight, individual correction off, and face-flux
resampling on. Run the serial and four-rank 200,000-particle mover suites for
both movers, the 46-reaction ensemble, parser/database/unit tests, a DIM3
compile, and at least one unrelated AMPS build with the feature left off.
Then perform the accepted long Halley run and regenerate the Rubin comparison.
The integration branch must reproduce bounded population, have no buffer or
continuation failure, and satisfy the approved convergence criterion. Finish
with `git diff --check` and inspect the final source-only file list.

Once those checks pass, record the pending merge on the integration branch,
but do not update `master` yet:

```sh
git status --short
git diff --cached --check
git commit
git log --oneline --decorate master..integrate/comet-local-weight
git diff --stat master...integrate/comet-local-weight
```

**Stop a second time and request user review.** Only after a separate explicit
approval may `master` be fast-forwarded:

```sh
git switch master
git merge --ff-only integrate/comet-local-weight
```

If the fast-forward fails because `master` moved, do not force, reset, or
merge directly into it. Create a new disposable integration branch from the
new `master`, repeat conflict resolution and the full validation matrix, and
request review again. Keep `1d_comet_dev`, the integration branch, and the
validated run archive until post-merge DIM1/DIM3 smoke tests are complete.

## Species-local block-weight adaptation added on 2026-11-02

This round adds an application-level controller for the already implemented
species-local block weights. It does not change a core source file. The Halley
case explicitly selects

```text
RadialParticleWeightAdaptation=on
RadialParticleNumberPerCellRange=50 | 200
```

The range is interpreted separately for each species. At a complete sampling
break, the instantaneous count of that species in an owner block is divided
by all eight DIM1 cells in the block. Thus an empty cell contributes zero to
the block mean. The population-range controller skips a block/species pair
with no particles, and every clone remains in its source cell; adaptation
never seeds an empty cell. The later 2026-11-03 enhancement below may raise
the statistical weight of an empty photochemical daughter without creating a
particle, preventing its first reaction event from producing a clone burst.

The selected block weight is restricted to an integral power of two relative
to its old value. If the mean is below 50, each factor-two step halves the
weight and clones every existing particle of that species exactly once. If
the mean is above 200, each step doubles the weight and retains each particle
with probability `W_old/W_new=1/2`. The first operation preserves represented
moments exactly; the second preserves them in expectation. The input parser
requires `maximum >= 2*minimum`, which guarantees that some power-of-two
image of every positive population lies in the inclusive requested interval.

Adaptation occurs only after a full sampling output has been collected and
used by the convergence calculation. The completed sample is therefore not
retroactively rescaled. It is not run when convergence has already been
accepted and the simulation is returning. Before any split, the code computes
the exact total clone demand for the rank and verifies particle-buffer
capacity, preventing a half-completed transaction. It also verifies that no
particle has a pending local-timestep continuation.

Only owner blocks are modified. After at least one change, the standard AMPS
block-data exchange synchronizes local weights to MPI boundary copies. The
legacy per-species global-minimum weight cache is then recomputed because
adaptation can raise every local block weight, whereas the existing setter
only lowers that cache. Future block crossings and chemical species changes
are connected through the previously validated full `(W/dt)` source-to-target
resampling; no individual particle-weight correction is introduced.

### Files changed in this adaptation round

| File | Detailed reason and merge significance |
| --- | --- |
| `input/hartley2.pl` | Adds the explicit ON/OFF keyword and the `MINIMUM | MAXIMUM` range, generates the three radial application variables, validates positive integral bounds and the factor-two reachability condition, and requires an explicit switch for every enabled radial simulation. Preserve the existing generic comet/chemistry/physics handlers when resolving parser conflicts. |
| `input/halley1d.input` | Enables adaptation and records the supplied 50--200 per-cell range, species/block interpretation, factor-two mechanics, and empty-cell rule. These are numerical AMPS controls, not values attributed to Rubin et al. |
| `srcHartley2/RadialSimulation.h` | Declares the generated controls and the focused validator. The public validator is intentionally a replacement execution mode because it destroys artificial test particles. |
| `srcHartley2/RadialSimulation.cpp` | Implements per-species owner-block counting, power selection, exact splitting, unbiased roulette, clone-capacity preflight, continuation-state checks, MPI weight synchronization, global-minimum-cache repair, diagnostics, sampling-break invocation, and integrated split/roulette validation. It also initializes the core owner-block table before focused pre-timestep validators; normal `PIC::TimeStep()` already performs the same update. |
| `srcHartley2/RadialMain.cpp` | Adds the mutually exclusive `--validate-weight-adaptation` CLI path after ordinary radial initialization. |
| `srcHartley2/README.md` | Documents input syntax, equations, timing, conservation properties, empty-cell behavior, MPI synchronization, limitations, and validator commands. |
| `srcHartley2/halley2011/PAPER_PARAMETERS.md` | Separates the 50--200 adaptation policy from published Halley parameters and states exactly how its mean is evaluated. |
| `srcHartley2/STATUS.md` | Records this implementation, reproducible evidence, limitations, and future master-integration procedure. |

No file under `src/pic/` is newly modified for this round. The controller
deliberately reuses the core's local-weight accessors, owner-block table,
particle-list API, block-data exchange, and previously developed local-face
`W_old*dt_new/(W_new*dt_old)` conversion. Those earlier core changes remain a
prerequisite and must not be dropped during a later master merge.

### Focused validation evidence

Before compilation, `/data/vtenishe/COMET1D/AMPS/build` was deleted. The exact
Halley input then cleanly preprocessed, compiled, and linked with
`make -j8 EXTRACOMPILEROPTIONS=-O3`. The known DIM1 legacy aggressive-loop
warnings and sandbox-only OpenMPI interface notices remained, but there was no
new compile or link error. The generated source contained adaptation enabled
with bounds 50 and 200; the generated core configuration retained species-
local timestep, species-local block weight, individual correction off, and
local-face resampling on.

The frozen executable, input, macro file, and logs are under
`validation_runs/halley1d_weight_adaptation_20261102`. Their executable,
input, and macro SHA-256 digests are respectively
`3f1aac2c14f0f4f766c33bd5e6157a4e2e5e9df186c778e7efd8caf967c203b3`,
`47d54396c80b17f0dc52024f6c01bb23b2252f440c7d762b787e09c52ffd2ef3`,
and `36011722e4d1827ba0c73464dc27fb5a8b0e2efe87a3f6844beec68a9d816183`.

`mpirun -np 1 ./amps.weight-adaptation --validate-weight-adaptation`
passed. The split case changed 50 particles in one populated cell to 400,
moving the eight-cell block mean from 6.25 to 50 by three exact steps. The
roulette case changed 3,200 particles to 1,586, an observed mean of 198.25
within the 50--200 band and the validator's six-sigma binomial bound. All
other cells remained empty.

`mpirun -np 4 ./amps.weight-adaptation --validate-weight-adaptation` also
passed with one independently adapted owner block on each rank. The aggregate
split was 200 to 1,600 particles; the aggregate roulette was 12,800 to 6,383;
and the final report was
`PASS ... owner-blocks=4 ghost-weight-checks=8`. The eight direct checks are
non-owner copies of selected inter-rank boundary blocks after both split and
roulette synchronization; the decomposition contained six MPI boundary
blocks. `perl -c input/hartley2.pl`,
`input/test_comet_reaction_network.pl`, and
`input/test_comet_physics_database.pl` all passed. `git diff --check` is part
of the final source audit.

These tests validate controller mechanics and MPI synchronization. They do
not replace a new 10,000-iteration Halley convergence/population run with the
controller enabled. The earlier non-adaptive production run remains the last
long science result, so the feature is not yet accepted for a paper claim or
for merging to `master`.

### Merge notes specific to this round

Do not merge yet. Continue development and long validation on
`1d_comet_dev`, and report to the user before any switch, rebase, cherry-pick,
or merge involving `master`. When the long adaptive run is accepted, add this
round to the existing temporary integration-branch runbook in three semantic
batches: (1) parser and Halley controls; (2) radial application implementation
and focused CLI test; (3) README, paper-parameter, and status documentation.
The implementation is application-local, but batch 2 depends on the earlier
core local-weight/face-transfer batch and must follow it.

Resolve conflicts semantically. In `RadialSimulation.cpp`, retain both the
source-derived initial `W/dt` profile and the new sampling-break adaptation;
do not replace one with the other. Retain owner-only resampling, exact
factor-two steps, `W_old/W_new` roulette, preflight capacity checking,
`ParallelBlockDataExchange()`, global-minimum recomputation, and the rule that
converged output returns before another adaptation. In the parser, preserve
the explicit switch requirement and do not silently enable this behavior for
other radial comet inputs. Keep the range in the input rather than hard-coding
50 and 200 in C++.

For a source-only commit, explicitly review and add only the eight files in
the table above. Do not add generated root `hartley2.pl`, `build/`, `amps`,
`PT*`, meshes, AMR measure caches, logs, or `validation_runs/`. The unrelated
local `Makefile` linker edit remains outside this feature. After future
conflict resolution, delete `build`, rebuild `halley1d.input`, confirm the
generated controls and core mode quartet, rerun both focused adaptation tests,
the existing serial/four-rank face-mover suite, chemistry and parser tests,
then complete the approved long convergence and Rubin-comparison run. Stop
again for user review before recording or applying any merge to `master`.

## Photochemical first-daughter burst protection added on 2026-11-03

The adaptive-weight Halley validation under
`validation_runs/halley1d_adaptive_weight_convergence_rubin_20261102`
terminated at iteration 5516 after approximately 8 h 31 min on four MPI
ranks. Rank 2 reported `pic_pbuffer.cpp` error 225.1, "The particle buffer is
full." The last complete break at iteration 5500 had globally reduced the
population from 14,133,960 to 3,736,445, yet the failure followed only 16
iterations later. The preceding controller adapted only species already
present in each block. An absent daughter therefore retained an obsolete low
carrier rate after its parent weight increased; its first later chemical event
could request a very large residual-resampling multiplicity even though the
global post-adaptation population was modest. The last convergence metric was
1.669883, limited by CN, so this failed run is neither converged nor suitable
for a Rubin comparison claim.

The input now explicitly supplies

```text
RadialMaxPhotochemicalProductsPerReaction=20
```

The parser accepts only a positive integer and requires this record whenever
DIM1 radial particle-weight adaptation and chemistry are both enabled. The
value is a numerical AMPS population-control choice, not a parameter inferred
from Rubin et al. The generated application also validates the value before
initialization.

At each complete sampling break, ordinary 50--200 population plans are built
first. Starting with every species having a nonzero population in an owner
block, the controller closes the enabled unit-unimolecular reaction graph
transitively. Thus, a daughter that is currently absent but can be created and
react again before the next sampling break is protected as a later-generation
parent. For each reachable reaction, if its absent products' stoichiometric
coefficients sum to `Cabs`, the minimum carrier rate assigned to every absent
product is

```text
(W_product/dt_product)_min =
    (W_parent/dt_parent) * Cabs/20.
```

The parent carrier in this expression includes the power-of-two population
change planned in the same transaction. Reachability and carrier constraints
are iterated to a fixed point. Therefore, for every reachable generation, the
sum of the absent products' expected residual-resampling multiplicities for
one individual selected reaction event is at most 20. When several parents or
channels reach the same absent species, the largest required carrier is
retained. The daughter weight is only raised, by an integral power of two.
Because there is no particle list to resample, the operation creates no
particles and leaves every empty cell empty. The first later chemistry event
uses the new block weight through the existing unbiased
`(W_parent/dt_parent)/(W_product/dt_product)` residual sampling. Products
already present remain governed by the independent 50--200 mean-population
controller.

Every adaptation now reports `empty-photochemical-daughter` and the global
maximum expected absent-product multiplicity before and after the transaction.
The code checks the post-adaptation maximum against the input limit and aborts
before chemistry resumes if it is violated. Updated owner weights continue to
use the standard block-data exchange and global-minimum-cache reconstruction,
so MPI boundary copies see the same protected daughter weight.

### Files changed in this protection round

| File | Detailed reason and later-merge significance |
| --- | --- |
| `input/hartley2.pl` | Adds and validates `RadialMaxPhotochemicalProductsPerReaction`, generates the C++ value, and requires it only for the combination of enabled radial adaptation and enabled chemistry. Preserve the existing generic reaction-network generator and order-independent input handling during conflict resolution. |
| `input/halley1d.input` | Supplies the user-selected value 20 and documents that the limit is the total expected model multiplicity of all currently absent products in one selected channel. It is a numerical setting, not a Rubin parameter. |
| `srcHartley2/RadialSimulation.h` | Declares the generated cap and its block-local, absent-daughter semantics. |
| `srcHartley2/RadialSimulation.cpp` | Reads the generic generated reaction network, transitively closes all reactions reachable from present species, iterates absent-daughter carrier constraints to a fixed point using post-adaptation parent carriers, applies weight-only power-of-two increases, verifies the cap, extends diagnostics, and adds an integrated serial/MPI regression. No particle is created by the new planner. |
| `srcHartley2/README.md` | Documents the new required input, derivation, conservation interpretation, diagnostics, scope, and focused validator behavior. |
| `srcHartley2/halley2011/PAPER_PARAMETERS.md` | Separates the value 20 from published Halley inputs and records it as an AMPS numerical control. |
| `srcHartley2/STATUS.md` | Records the failed-run evidence, implementation, validation, limitations, and merge procedure in this section. |

The root `hartley2.pl` used by the AMPS preprocessing command was refreshed
mechanically from `input/hartley2.pl`; it remains an untracked generated
working copy and must not be included in a source commit. No file under
`src/pic/` or any other core directory was modified in this round. The change
uses the existing generic reaction-network API, local block weight/timestep
accessors, particle lists, and MPI block-data exchange.

### Build and focused validation evidence

Before the accepted compilation, `build` was removed and the exact
`halley1d.input` configuration was regenerated. The optimized command
`make -j8 EXTRACOMPILEROPTIONS=-O3` completed and linked successfully. The
generated `build/main/RadialSimulation.cpp` contains adaptation enabled,
the 50--200 range, and
`MaximumPhotochemicalProductsPerReaction=20`. The known DIM1 aggressive-loop
and sandbox OpenMPI interface warnings remained; there was no final compile or
link error. `perl -c input/hartley2.pl`, the reaction-network parser test, and
the comet physics-database test all passed.

The focused validator deliberately set an enabled channel's absent products
to an unsafe total expectation of 160 model particles per individual event.
It evaluates the full reaction reachability closure when checking this bound.
The serial run changed two empty daughter weights, created no daughters,
reduced the expectation from 160 to exactly 20, and reported PASS. The
four-rank run changed eight empty daughter block/species weights across four
owner blocks, again reduced 160 to 20, and passed four photochemical checks
plus 16 direct ghost-weight synchronization checks. Its final line was:

```text
PASS particle-weight adaptation validation range=50:200 max-photochemical-products=20 owner-blocks=4 photochemical-checks=4 ghost-weight-checks=16
```

These focused results validate the cap calculation, power-of-two weight-only
operation, empty-cell preservation, and MPI synchronization. They do not prove
long-run population stability or convergence. A fresh 10,000-iteration Halley
run must pass the former iteration-5516 failure point without a buffer error,
then satisfy the existing convergence criterion and regenerate the Rubin
comparison before this enhancement is accepted.

### Merge notes for this protection round

Do not merge yet. Continue work on `1d_comet_dev` and tell the user before any
operation involving `master`. During a later temporary integration-branch
merge, place this round after the existing generic chemistry and species-local
weight batches. Resolve `RadialSimulation.cpp` semantically: retain both the
50--200 present-species controller and the absent-daughter cap; calculate
daughter requirements from the parent's *planned* carrier; preserve the
transitive reachability/fixed-point iteration and maximum requirement over all
reachable channels; never seed an empty particle list; retain the postcondition
check, owner/ghost synchronization, and global-minimum repair. Do not replace
the generated generic network with a Halley-specific reaction table or
hard-code 20 in C++.

After conflict resolution, remove `build`, regenerate the exact Halley input,
confirm the generated cap and local timestep/weight mode quartet, and rerun
the parser/database tests plus serial and four-rank adaptation validators.
Then perform the full convergence/population/Rubin-comparison run and verify
that diagnostics never report a post-adaptation absent-product expectation
above 20. Include only the seven maintained source/documentation files in the
table above; exclude root generated copies, `build`, `amps`, `PT*`, meshes,
logs, caches, and `validation_runs`. Stop for user review before recording an
integration commit or changing `master`.

## Emergency pre-exchange particle-buffer protection added on 2026-11-03

The 12-rank run in
`validation_runs/halley1d_photochemical_cap_convergence_rubin_20261103_12mpi`
terminated at iteration 5,690 after approximately 5 h 15 min. Rank 7 aborted
in `pic_pbuffer.cpp:476` with error 225.1, "The particle buffer is full." The
iteration-5,500 sampling adaptation had reduced the global population from
12,209,935 to 3,706,593, but the population subsequently grew to 29,159,492
before the next sampling break. NH2 alone accounted for 16,021,072 particles
and NH3 for 6,381,797. The last complete convergence metric was 1.466266,
limited by CH with `0/3` accepted windows. Thus the run is useful failure
evidence but is not a converged Rubin validation result.

The Halley input now explicitly supplies:

```text
RadialEmergencyParticleBufferUtilization=0.75
```

The parser accepts only a finite fraction in `(0,1)` and requires the keyword
when radial particle-weight adaptation is enabled. The value is an AMPS
numerical population-control setting, not a quantity inferred from Rubin et
al. The configured 0.75 means that the emergency path is entered when the
largest local particle count among MPI ranks reaches 75% of that rank's
`MaxNPart` capacity.

### Algorithm and conservation

The generic core particle-exchange wrapper now contains a nullable callback
executed collectively immediately before either buffered or unbuffered normal
particle migration. Other applications retain the historical path because the
default pointer is null. `RadialSimulation::Initialize()` installs the comet
callback only when its existing weight-adaptation switch is on.

The callback first uses `MPI_Allreduce` to obtain the maximum rank utilization.
If it is below the input threshold, it returns without changing any state. If
it reaches the threshold, one emergency step visits every allocated owner and
MPI boundary-block copy on every rank, retains each particle independently
with probability one half, and doubles every block-local weight for every
species. It then rebuilds the core global-minimum weight cache and measures
utilization again. The factor-two step repeats only if stochastic survival
left a rank at or above the threshold.

The common factor across all blocks and species is essential at this point in
the iteration. The mover may already have placed completed or continuation
particles in a non-owner boundary-block copy while waiting for the exchange.
Changing only an owner block would give those staged particles inconsistent
weights. Applying the same factor everywhere preserves every block-to-block
weight ratio and every inter-species `W/dt` carrier ratio. Consequently the
existing face-flux and chemistry resampling probabilities are unchanged.
Surviving particles retain position, velocity, internal state, and any local-
timestep remainder. Roulette conserves represented moments in expectation,
does not split particles, and never seeds an empty cell.

This is deliberately an emergency variance-for-capacity trade. It is not a
replacement for the ordinary species/block 50--200 controller at completed
sampling breaks. It also cannot guarantee protection against an arbitrarily
large allocation burst that occurs before the next pre-exchange check; the
25% headroom and the separate per-reaction product cap must be assessed in a
new long run.

### Files changed and later-merge significance

| File | Detailed reason and merge guidance |
| --- | --- |
| `src/pic/pic.h` | Adds the generic `BeforeParticleDataExchange` callback declaration and invokes it immediately before the exchange-mode dispatch. During a master merge, retain the callback before both buffered and unbuffered implementations and do not move it after receive-slot allocation. If master already has an equivalent hook, use that hook instead of creating a duplicate. |
| `src/pic/pic_parallel.cpp` | Defines the new callback pointer as null so all unrelated AMPS applications remain behaviorally unchanged. Keep this definition paired with the header declaration. |
| `input/hartley2.pl` | Adds strict parsing/generation for `RadialEmergencyParticleBufferUtilization` and requires the explicit value when radial weight adaptation is enabled. Preserve the existing order-independent generic comet parser and do not introduce a hidden default. |
| `input/halley1d.input` | Supplies 0.75 and documents trigger timing, factor-two roulette, staged particles, and the reason a common weight factor preserves `W/dt` ratios. This value is numerical, not from Rubin et al. |
| `srcHartley2/RadialSimulation.h` | Declares the generated threshold and the callback whose address is stored by the core. |
| `srcHartley2/RadialSimulation.cpp` | Implements collective occupancy checks, all-copy/all-species factor-two roulette, cache reconstruction, diagnostics, callback installation, input validation, and an integrated test through the real exchange API. Preserve the distinction between this global emergency step and the owner-only per-block sampling adaptation. |
| `srcHartley2/README.md` | Documents syntax, trigger ordering, conservation, variance, interaction with staged continuation particles, and limitations. |
| `srcHartley2/halley2011/PAPER_PARAMETERS.md` | Records 0.75 as a supplied AMPS control rather than a published Halley parameter. |
| `srcHartley2/STATUS.md` | Records the failure that motivated the change, implementation, focused evidence, limitations, and this merge procedure. |

The root `hartley2.pl` is only the untracked mechanical preprocessing copy of
`input/hartley2.pl`; it must not be committed. Generated `build/`, `amps`,
`PT*`, mesh/cache files, logs, and `validation_runs/` also remain outside a
source commit.

### Clean build and focused evidence

`build` was explicitly deleted before preprocessing and compilation. The exact
tracked `halley1d.input` then generated
`EmergencyParticleBufferUtilization=0.75`, species-local timesteps,
species-local block weights, individual weight correction off, and local-face
flux resampling on. `make -j8 EXTRACOMPILEROPTIONS=-O3` compiled and linked
successfully. Only the previously documented aggressive-loop and sandbox
OpenMPI interface warnings appeared. Both parser syntax checks and the comet
reaction/physics database tests passed.

The integrated validator temporarily sets a small synthetic trigger equal to
75% of its initial test utilization, then calls the ordinary core
`ExchangeParticleData()` API. This exercises callback placement without
allocating 15 million particles per rank merely to reach 75% of the production
20-million-slot buffer; the production value is restored and printed in the
final result.

The one-rank test reduced 4,096 particles to 1,957 in one factor-two step and
reported maximum utilization `2.048e-4 -> 9.785e-5` before the regular MPI
exchange. It passed one emergency and one photochemical check. The four-rank
test reduced 16,384 particles to 8,133, with the maximum rank changing from
4,096 to 2,074 and utilization `2.048e-4 -> 1.037e-4`. It passed four
emergency checks, four photochemical checks, and 20 direct ghost-weight
checks. The matching 12-rank test reduced 49,152 particles to 24,449 and
passed 12 emergency checks, 12 photochemical checks, and 60 direct ghost-
weight checks. Its final line was:

```text
PASS particle-weight adaptation validation range=50:200 emergency-buffer-utilization=0.750000 max-photochemical-products=20 owner-blocks=12 emergency-checks=12 photochemical-checks=12 ghost-weight-checks=60
```

These focused results validate trigger collectives, factor-two statistics,
core-hook timing, owner/boundary weight agreement, and all previous adaptation
checks. They do not establish long-run population stability or physical
convergence. Acceptance requires a fresh 12-rank run to pass iteration 5,690,
avoid all buffer failures through termination, satisfy the existing three-
window convergence rule, and regenerate Rubin comparisons from the accepted
final output.

### Merge instructions for this round

Do not merge yet. Remain on `1d_comet_dev` until the new long validation is
accepted, and tell the user before any switch, rebase, cherry-pick, or merge
involving `master`. In a later disposable integration branch, apply the two
generic core-hook files as a small semantic batch after the existing local-
timestep/local-weight core prerequisites. Apply the parser/application batch
next, followed by documentation. The application batch depends on the core
hook and must not be cherry-picked alone.

During conflict resolution, confirm that the callback remains null by default,
is called collectively before receive allocation, and cannot recurse through
`ExchangeParticleData()`. In the radial code, retain the common all-species,
all-block factor for emergency operation; do not replace it with the ordinary
owner-only per-species planner while particles may be staged in boundary
copies. Retain the post-step utilization check, global-minimum-cache rebuild,
and detailed diagnostic. Preserve the independent 50--200 sampling-break
controller and absent-daughter cap.

After conflict resolution, delete `build`, regenerate the Halley configuration,
verify the generated 0.75 value and mode quartet, rerun the parser/database
tests and serial/four-/twelve-rank adaptation validators, then perform the full 12-rank
convergence/population/Rubin study. Stop for user review before recording an
integration commit or changing `master`.

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

## Fourteen-rank post-emergency long validation stopped on 2026-11-03

The requested repeat run was archived in:

```text
/data/vtenishe/COMET1D/AMPS/validation_runs/halley1d_emergency_weight_convergence_rubin_20261103_14mpi
```

It was launched at 10:23:21 CST with:

```text
mpiexec -n 14 ./amps.emergency-weight
```

The directory contains immutable snapshots of the executable, maintained
Halley input, Rubin comparison utility, generic one-dimensional plotting
utility, and digitized Rubin reference data. The executable SHA-256 is
`1dbbd0bf0dce67ab518747fdb90ae5cdd8afb1739a92b42e9ccb78614453ff34`;
the input SHA-256 is
`6f9a7de309e3d0f1d3274d0aa4786ae63e6b7659a150ac9c7867f1591e280cd7`.
The snapshot requests 20,000,000 particle slots per rank, block-local species
weights and timesteps, regular 50--200 mean-particles-per-cell adaptation,
the per-event product cap of 20, and emergency pre-exchange adaptation at 75%
rank-local buffer utilization. Termination remains three consecutive complete
500-iteration windows below 5% relative change, with 10,000 iterations as the
hard ceiling.

The job reported exactly 14 MPI ranks and completed ten sampling windows. At
iteration 5,000, the convergence metric was `1.946078`, limited by CH, with
zero of three required converged windows. Regular adaptation changed the
population from 10,667,398 to 3,317,073 and reduced the maximum expected
absent-product multiplicity from 65,537 to 16. The run did not exercise the
75% emergency-buffer trigger. It was stopped at iteration 5,492 at the user's
request before reaching the former iteration-5,690 failure point. No particle-
buffer error or MPI abort occurred before the requested stop. This is useful
population-control evidence but is neither a convergence result nor an
accepted Rubin comparison. Do not merge it to `master`.

## Sampling-population emergency request added on 2026-11-03

The rank-local 75% buffer trigger does not detect a severe population
concentration until that concentration occupies most of one rank's entire
20-million-slot buffer. A second trigger now reuses the particle counts that
the regular sampling-break controller already calculates. For every owner
block, it sums all compiled species and forms the whole-block mean

```text
mean_total = sum_species(block particle count) / block cell count.
```

The threshold is the explicitly requested

```text
10 * number of compiled species * maximum requested particles/species/cell.
```

For the current 27-species, 50--200 Halley case, this is 54,000 total model
particles per cell as a block mean. The implementation compares integer block
counts against the equivalent threshold multiplied by the block cell count,
with checked integer multiplication, so a value exactly at the boundary does
not trigger and there is no floating-point boundary ambiguity.

The comparison uses the pre-adaptation sampling snapshot. If the global
maximum owner-block mean strictly exceeds the threshold, all MPI ranks record
the same pending request and a diagnostic containing the sampling iteration,
observed mean, threshold, multiplier, and species count. The sampling routine
does not call particle exchange or global emergency roulette recursively. The
next normal `PIC::Parallel::ExchangeParticleData()` call consumes the request
inside the existing `BeforeParticleDataExchange` callback. It applies at least
one common all-species/all-allocated-block factor-two weight increase and
unbiased one-half roulette before normal MPI particle migration. If the
independent rank-buffer trigger is also active, roulette continues until all
ranks are below its configured utilization threshold. Diagnostics distinguish
`sampling-population`, `buffer`, and combined trigger reasons.

This timing preserves the established safety properties. A particle already
staged in a non-owner boundary block is resampled together with ordinary owner-
block particles; every owner and ghost weight changes by the same factor; and
all block-to-block and inter-species `W/dt` ratios remain unchanged. The
sampling request is cleared only after roulette and the global-minimum weight
cache rebuild complete successfully. A later sampling window can raise a new
one-shot request.

### Files changed in this enhancement

| File | Change and later-merge significance |
| --- | --- |
| `srcHartley2/RadialSimulation.cpp` | Adds checked threshold construction, all-species block-count accumulation, the collective pending request, trigger-aware pre-exchange diagnostics, one-shot request consumption, and a focused integrated regression. Keep the count before ordinary per-species resampling and keep execution in the pre-exchange callback when resolving master conflicts. |
| `srcHartley2/RadialSimulation.h` | Extends the emergency-control documentation to describe both independent trigger conditions. No new public variable or input keyword is introduced. |
| `input/halley1d.input` | Documents that the existing maximum of 200 also defines the fixed `10*nSpecies*maximum` sampling alarm. No Rubin parameter is inferred. |
| `srcHartley2/README.md` | Documents units, strict comparison, collective request timing, conservation, diagnostics, and interaction with the 75% trigger. |
| `srcHartley2/halley2011/PAPER_PARAMETERS.md` | Marks the new multiplier as a supplied numerical control rather than a paper value. |
| `srcHartley2/STATUS.md` | Records the stopped 14-rank run, algorithm, validation requirements, and merge guidance. |

No additional core source modification is required in this round. The
application continues to use the already-added generic nullable
`BeforeParticleDataExchange` hook in `src/pic/pic.h` and
`src/pic/pic_parallel.cpp`. The core hook must remain earlier than both normal
exchange implementations and receive-side allocation.

### Validation and merge requirements

The integrated weight-adaptation validator now constructs a small valid
synthetic 5--20 range, exceeds its exact
`10*nSpecies*maximum*block-cell-count` threshold, calls the real sampling
adaptation, verifies that a request was recorded but not executed there, then
calls the ordinary core particle-exchange API. It checks factor-two/binomial
statistics, request consumption, empty-cell preservation, and owner/ghost
weight agreement. The existing synthetic 75% buffer test remains independent.

Before accepting this change, delete `build`, compile the exact Halley input,
run parser/database checks, and pass the serial, four-rank, and twelve-rank
focused validators. Then start a fresh long convergence run and confirm both
trigger diagnostics when their respective conditions occur, population
stability past iteration 5,690, no buffer failures through termination, and a
new Rubin comparison from the accepted final window.

The clean-build and focused portion of that plan passed on 2026-11-03. The
workspace `build/` directory was removed before compiling the exact
`input/halley1d.input` configuration with `-O3`. The input parser and the
reaction-network and physical-database checks passed. The resulting executable
then passed `--validate-weight-adaptation` with 1, 4, and 12 MPI processes. The
sampling-population test deliberately observed a block mean of 10,800
particles per cell against its synthetic 5,400 threshold and produced a
`trigger=sampling-population` event only in the next pre-exchange callback.
The final MPI evidence was:

```text
1 rank:  owner-blocks=1  emergency-checks=1  sampling-population-emergency-checks=1  ghost-weight-checks=0
4 ranks: owner-blocks=4  emergency-checks=4  sampling-population-emergency-checks=4  ghost-weight-checks=24
12 ranks: owner-blocks=12 emergency-checks=12 sampling-population-emergency-checks=12 ghost-weight-checks=72
```

All three runs reported `PASS`. The remaining acceptance work is the new long
convergence/Rubin-comparison run described above; it has not been started by
this enhancement.

For a later master integration, do not cherry-pick this application change
without the earlier two-file core hook. Resolve `RadialSimulation.cpp`
semantically so that the ordinary per-species 50--200 controller, absent-
daughter product cap, 75% rank-buffer trigger, and new total-block-mean trigger
all remain independent. Do not move global roulette into the sampling routine,
do not change the strict `>` comparison to `>=`, and do not calculate the
sampling alarm after regular adaptation. Exclude generated `build/`, root
parser copies, executables, meshes, logs, `PT*`, and `validation_runs/` from a
source commit. Remain on `1d_comet_dev` and stop for user review before any
operation involving `master`.

## Fourteen-rank sampling-emergency validation failed at iteration 6,019

After the clean build and 1/4/12-rank focused regressions passed, a fresh long
Halley/Rubin validation was started from iteration zero with 14 MPI processes:

```text
/data/vtenishe/COMET1D/AMPS/validation_runs/halley1d_sampling_emergency_convergence_rubin_20261103_14mpi
mpiexec -n 14 ./amps.sampling-emergency
```

The previous stopped run directory and its iteration-5,492 evidence were left
unchanged. The new run uses the newly built executable whose SHA-256 is
`9fcbd3c0818767519dc1afc28aaa373580b39f3896f713bd20e33b55ba65486a`;
the staged `halley1d.input` is byte-identical to `input/halley1d.input`. Its
`PT/plots` tree was created normally, and eleven sampling windows had completed
through iteration 5,500 at the latest status check. The convergence metric was
1.592657 at iteration 5,500 (C2H2-limited, 0/3 consecutive windows), so this is
not yet a convergence result. The new sampling-population path received its
first long-run exercise at iteration 2,500: the maximum owner-block mean of
139,595.1 particles/cell exceeded the 54,000 threshold and raised a collective
request. At the next pre-exchange callback, one global factor-two step reduced
the represented model-particle population from 1,547,234 to 772,725 and the
maximum rank population from 229,100 to 114,351. The aggregate buffer use was
only 1.1455%, confirming that this event was caused by the new local-
concentration alarm rather than the independent 75% rank-buffer alarm. The run
triggered a second sampling-population event at iteration 3,500, where the
maximum block mean of 99,652.25 particles/cell again exceeded 54,000. The
pre-exchange step reduced the model population from 2,302,385 to 1,151,350 and
the maximum-rank population from 762,910 to 381,673; rank-buffer utilization
was only 3.81455%, again well below the independent 75% threshold. No third
emergency event occurred at iteration 4,000; regular adaptation changed the
population from 2,292,318 to 1,470,464. A third sampling-population event
occurred at iteration 4,500, with a block mean of 96,618.12 particles/cell.
The pre-exchange step reduced 2,646,904 model particles to 1,322,913 and the
maximum-rank population from 510,896 to 255,554. A fourth event at iteration
5,000 reduced 2,112,903 particles to 1,055,870 after a block mean of 68,605.5
particles/cell was detected. A fifth event at iteration 5,500 caught a much
larger 262,430.6-particle/cell block mean and reduced 1,880,894 particles to
940,366; the maximum-rank population fell from 376,750 to 188,505. The run
subsequently passed the former iteration-5,690 particle-buffer failure point
and advanced through iteration 5,950 without an AMPS, MPI, or particle-buffer
error. This is direct long-run stability evidence for the new sampling alarm,
although the model remained at 0/3 converged sampling windows.

The run did not complete. At iteration 6,000, its convergence metric was
1.635755 (H2-limited, 0/3). A sixth sampling alarm detected an extreme block
mean of 524,326.8 particles/cell. Regular 50--200 adaptation first reduced the
population from 6,432,566 to 1,675,482, and the deferred global emergency step
then reduced 1,678,520 particles to 839,312; maximum-rank occupancy fell from
339,908 to 169,649. Nevertheless, the population regrew extremely rapidly. At
iteration 6,008 the global population was 2,242,648 and rank 6 already held
1,416,169 particles. During iteration 6,019, before another emergency callback
could create headroom, rank 6 exhausted its 20-million-slot particle buffer:

```text
AMPS: exit: line=476, file=pic_pbuffer.cpp,
message=The particle buffer is full (error code=225.1)
MPI_ABORT was invoked on rank 6
```

The MPI launcher exited with status 1 and no simulation processes remain. This
failure demonstrates that sampling-window detection plus a pre-exchange 75%
occupancy test cannot protect against a very large intra-iteration allocation
burst. A subsequent fix must provide allocation-site or chemistry-stage
headroom control rather than merely lowering the sampling threshold. The run
is useful population-stability evidence, but it is not an accepted convergence
or Rubin-comparison result and must not be merged to `master` as validation of
the current algorithm.

## Empty-neighbor carrier propagation added on 2026-11-03

The iteration-6,019 failure above need not have been caused by chemistry alone.
The block-face mover uses the residual-resampling mean

```text
mu = (W_source/dt_source)/(W_destination/dt_destination).
```

Regular adaptation previously skipped a zero-population block/species pair.
Such a block could therefore retain an old, extremely small local statistical
weight. The first particle transported into it could have `mu >> 1` and create
a large clone burst before either the next 500-iteration sampling break or the
pre-exchange emergency callback. The chemistry daughter cap cannot control
this transport path.

`srcHartley2/RadialSimulation.cpp` now prepares empty transport regions in the
same sampling-break adaptation transaction. It gathers the authoritative
owner count and the planned carrier rate `W/dt` for every calculation leaf and
species onto every MPI rank. For each species, it constructs face adjacency,
finds every connected component whose blocks are all empty for that species,
and activates only components touching a currently populated block. The target
carrier of an active component is the maximum of all current empty-member
carriers and all populated boundary carriers. Every empty member below that
target is raised through exact factor-two weight steps. No particle is
allocated, removed, or cloned by propagation.

An explicit work list and multi-source breadth-first traversal carry the
constraint through a neighbor, its neighbor, and so on, with no recursive
call-stack limit. A dyadic fixed-point pass preserves the existing invariant
that adaptive weights differ by powers of two. After planning, the code checks
every directed face entering an active empty block and aborts unless its clone
mean is at most one within numerical tolerance. The established owner-to-ghost
block-data exchange then distributes the new weights before particle transport
resumes. Diagnostics add `empty-transport-neighbor`, `empty-components`,
`max-empty-hops`, and `max-empty-face-clone-mean=before->after`.

The adjacency builder is intentionally structured for later DIM2/DIM3 use: it
enumerates all `2*DIM` faces and the dimension-dependent face subslots. The
definition remains face connected because a corner contact is not a mover
crossing. Before accepting this path for higher dimensions, add dedicated
DIM2/DIM3 tests containing coarse/fine AMR faces, verify all active subfaces,
deduplicate repeated neighbors, and enforce the same directed-entry
postcondition. The present evidence validates DIM1 only.

### Files changed in this round

| File | Detailed reason and later-merge significance |
| --- | --- |
| `srcHartley2/RadialSimulation.cpp` | Adds global leaf/count/carrier collection, dimension-aware face adjacency, empty-component discovery, multi-hop work-list propagation, dyadic carrier raising, postcondition checks, diagnostics, and focused regression assertions. Preserve the `W/dt` calculation; substituting raw `W` breaks local-timestep flux balance. |
| `srcHartley2/RadialSimulation.h` | Documents that ordinary adaptation includes weight-only propagation into face-connected empty neighbors. It adds no public control or ABI object. |
| `input/halley1d.input` | Documents the automatic transport-consistency step next to the existing adaptation controls. No new keyword or Rubin-derived value is introduced. |
| `srcHartley2/README.md` | Describes the mover equation, component algorithm, conservation/empty-cell behavior, diagnostics, and requirements for a future higher-dimensional validation. |
| `srcHartley2/halley2011/PAPER_PARAMETERS.md` | Separates this numerical population safeguard from Rubin et al. physics and published parameters. |
| `srcHartley2/STATUS.md` | Records motivation, implementation, validation evidence, long-run disposition, and master-integration guidance. |

No core file under `src/pic/` was added or modified for this enhancement. It
uses the previously implemented local-timestep/local-weight face conversion,
local weight accessors, owner block table, and owner/ghost block-data exchange.
Those earlier core prerequisites remain necessary but should be reviewed as a
separate semantic batch during a future integration.

### Clean build and focused validation

Before compilation, the exact workspace `build/` directory was removed. The
Halley configuration then compiled and linked from scratch with
`make -j8 EXTRACOMPILEROPTIONS=-O3`. The reaction-network parser and comet
physics-database tests passed, and `git diff --check` reported no whitespace
error.

The focused `--validate-weight-adaptation` test passed with 1, 4, and 12 MPI
processes. Its roulette case deliberately creates a populated-to-empty carrier
ratio of two and requires propagation through successive empty blocks. Its
sampling-population stress case creates a clone mean of 1,024. In all cases the
post-propagation maximum was exactly one, no empty cell was populated, and the
multi-rank owner/ghost checks passed. The final evidence was:

```text
1 rank:  empty-transport-neighbor=372 max-empty-hops=372 max-empty-face-clone-mean=1024->1 empty-transport-propagation-checks=1  ghost-weight-checks=0
4 ranks: empty-transport-neighbor=369 max-empty-hops=97  max-empty-face-clone-mean=1024->1 empty-transport-propagation-checks=4  ghost-weight-checks=24
12 ranks: empty-transport-neighbor=361 max-empty-hops=32 max-empty-face-clone-mean=1024->1 empty-transport-propagation-checks=12 ghost-weight-checks=72
```

These tests establish deterministic multi-hop DIM1 planning, collective MPI
agreement, and ghost synchronization. They do not establish long-run
population stability, convergence, or agreement with Rubin et al.; those
remain requirements of the restarted validation run.

### Merge instructions for the empty-neighbor enhancement

Do not merge this work yet. Keep development and long validation on
`1d_comet_dev`, and stop for user review before any switch, rebase,
cherry-pick, or merge involving `master`. In a later disposable integration
branch based on the then-current master, first integrate and validate the
earlier generic core local-time/local-weight face machinery and nullable
pre-exchange hook. Then integrate this application change with the complete
sampling-break controller; do not extract only the empty-neighbor function
without its power-of-two plans and owner/ghost synchronization.

Resolve `RadialSimulation.cpp` semantically. Retain the order: count owner
particles, create ordinary population plans, constrain absent chemistry
daughters, propagate planned `W/dt` through empty face components, preflight
capacity, apply particle resampling and weight-only plans, synchronize ghost
weights, and rebuild the global-minimum cache. Preserve the postcondition and
diagnostics. Never replace `W/dt` with `W`, seed empty blocks, use corner
adjacency, or lower an empty carrier during propagation.

After conflicts are resolved, delete `build`, regenerate from
`input/halley1d.input`, repeat the parser/database tests and the 1/4/12-rank
focused validator, and run the full 14-rank convergence/population/Rubin
comparison. Exclude generated executables, `build/`, `PT*`, root parser
copies, meshes, logs, and `validation_runs/` from source commits. A successful
long run and regenerated comparison plots are prerequisites for review, not
permission to merge automatically.

### Restarted fourteen-rank long validation

A fresh run was launched from iteration zero at 18:46:50 CST on 2026-11-03:

```text
/data/vtenishe/COMET1D/AMPS/validation_runs/halley1d_empty_neighbor_weight_convergence_rubin_20261103_14mpi
mpiexec -n 14 ./amps.empty-neighbor-weight
```

The failed sampling-emergency run and all older archives remain unchanged.
The new directory contains the executable, exact Halley input, generic DIM1
plotter, Rubin comparison program, and digitized Rubin reference data. The
staged executable SHA-256 is
`2fa06c798ec8627cea87cb3eca4ec773dff25bef83742ff82b11de1ad75ef05d`;
the staged input SHA-256 is
`d2199c0af07ff346c934d7934d085bd4fc51ab2edc8cb5d09420adee91e19b96`.
Both staged files were byte-compared with the tested workspace artifacts.

Initialization reported exactly 14 MPI processes and one OpenMP thread per
process, and `PT/plots` was created normally. The run completed two sampling
windows and advanced beyond iteration 1,000 at this status update. Iteration
500 provided a production exercise of the new transport safeguard:

```text
radial particle-weight adaptation iteration=500 range=50:200
particles=39730->168757 changed-block-species=351 split=279 roulette=4
empty-photochemical-daughter=68 empty-transport-neighbor=671
empty-components=18 max-empty-hops=371
max-empty-face-clone-mean=6.400000e+01->1.000000e+00
max-absent-products-per-reaction=3.343333e+02->6.208333e+00 limit=20
```

Thus the real Halley population exposed a prospective factor-64 first-arrival
clone amplification, and the recursive empty-neighbor preparation reduced it
to one before transport resumed. At iteration 1,000 the maximum relative L1
change was `0.8181504` (NH3-limited, `0/3` consecutive windows). The second
adaptation reduced 1,520,239 particles to 484,013; no empty weight needed a
further increase and the maximum empty-entry clone mean remained one. No
particle-buffer, MPI, or postcondition error was present at this checkpoint.
The job remains active. It is not yet a convergence result or an accepted
Rubin comparison; acceptance still requires stable completion or convergence,
particular scrutiny past iteration 6,019, and plots regenerated from the
accepted final output.

## Allocation-burst control implemented and focused-tested on 2026-11-04

The preceding “job remains active” statement records the state at that earlier
checkpoint; the run later failed and is no longer active. The 14-rank
empty-neighbor run reached iteration 8,136, after completing its iteration-
8,000 sampling window with a maximum relative L1 change of `1.810334` for CN
and `0/3` converged windows. Rank 5 then exhausted its 20-million-slot particle
buffer inside a sampling window:

```text
AMPS: (11/3 23:36:43), Iteration: 8136
AMPS: exit: line=476, file=pic_pbuffer.cpp,
message=The particle buffer is full (error code=225.1)
MPI_ABORT was invoked on rank 5
```

The pre-exchange high-water check could not stop an allocation burst that
occurred before the next exchange. The empty-neighbor repair also left two
unbounded cases: a particle could enter a populated block with a much smaller
carrier, and residual chemistry could round an expected daughter count upward.
The current implementation controls the allocation sites rather than relying
only on a later population alarm.

### Implemented invariants and execution order

The complete controller has four coupled stages:

1. **Sampling-break and initialization planning.** Ordinary 50--200
   block/species population plans remain dyadic. The same planning pass closes
   the reaction graph from represented particles and configured surface-source
   species. For every reachable selected channel it enforces the hard bound
   `sum coefficient*ceil(C_parent/C_product) <= 20`, where `C=W/dt`. This
   includes products that are present, products that are absent, residual
   rounding, and later-generation daughters. Initialization performs the
   weight-only closure before the first source injection.
2. **Face allocation elimination.** For each species, every calculation leaf
   is placed in a graph using real mesh faces and all dimension-dependent
   subface slots. Because the requested maximum multiplicity is one in both
   directions, every block in a face-connected component is raised to a common
   dyadically reachable carrier. Populated blocks use unbiased factor-two
   roulette; empty blocks change weight only. Every directed face is checked
   after application, and owner values are distributed to ghost copies before
   transport resumes.
3. **Transactional chemistry.** The core photolysis driver makes bounded
   passes through eligible parents. The comet processor draws the complete
   event multiplicity before mutating the parent, reuses the parent slot for
   the first daughter, and commits only if every remaining daughter slot fits
   below the 90% transaction ceiling. Otherwise it reinserts the unchanged
   parent for the next pass. A per-particle stage byte prevents completed
   parents and new daughters from reacting again in the same outer iteration.
4. **Collective allocation emergency.** A chemistry checkpoint at the 75%
   high watermark requests common factor-two roulette and repeats it until the
   maximum MPI-rank population is at or below the 50% low watermark. The
   existing buffer and sampling triggers retain their 75% target behavior.
   Ten percent of each buffer remains outside application chemistry
   transactions. The callback is collective and never recursively invokes
   particle exchange.

The common-factor emergency conserves represented quantities in expectation
and preserves all carrier ratios. Ordinary block-specific roulette also
conserves in expectation and increases Monte Carlo variance. With a connected
DIM1 mesh and face multiplicity one, each species has a component-global
carrier after an adaptation pass; local weights still differ in proportion to
local timesteps. This is a deliberate zero-allocation transport policy, but it
is stronger than independently optimized block-local statistical weights.

### Files changed for the allocation-burst implementation

| File | Detailed modification and later-merge significance |
| --- | --- |
| `src/pic/pic.h` | Declares the nullable pre-exchange hook and the five nullable transactional-photochemistry callbacks. The complete batch set is all-or-none; unrelated applications keep null defaults. Merge declarations with their definitions and driver changes, not as an isolated header edit. |
| `src/pic/pic_parallel.cpp` | Defines the null pre-exchange callback used by the radial emergency controller. The invocation remains in `PIC::Parallel::ExchangeParticleData()` before buffered/unbuffered dispatch and receive-slot allocation. If current master has an equivalent hook, reconcile semantically rather than installing a duplicate. |
| `src/pic/pic_photolytic_reactions.cpp` | Defines null batch callbacks and extends the historical photolysis pass with bounded repeated passes, eligibility filtering, a collective application checkpoint, a global remaining-parent count, partial-set rejection, and a loop guard. The null path retains historical behavior. Installed batch mode explicitly rejects hybrid/OpenMP compilation until counters and cell-list reconstruction are made thread-safe. |
| `srcHartley2/CometChemistry.h` | Documents and exports the particle-byte reservation, batch callbacks, eligibility protocol, and boolean transactional processor contract. |
| `srcHartley2/CometChemistry.cpp` | Reserves/resets the stage flag, calculates a conservative parent budget from the high watermark and hard event cap, performs collective pressure checks, and makes each reaction atomic. It never truncates products: an event above the hard cap aborts, while a capacity-deferred event leaves the parent unchanged. |
| `srcHartley2/RadialSimulation.h` | Declares and documents high/low/reserve/face controls plus the collective allocation-emergency entry point. The photochemical input is documented as a hard actual allocation cap rather than an expected absent-product cap. |
| `srcHartley2/RadialSimulation.cpp` | Implements direct-source reaction reachability, hard `ceil`-based chemistry planning, all-face component carrier equalization, populated-block roulette, final post-face chemistry checks, low-water allocation emergencies, initialization-time planning, input validation, diagnostics, and expanded focused regressions. The face graph enumerates DIM1/2/3 topology, but only DIM1 has been validated. |
| `input/hartley2.pl` | Strictly parses and requires `RadialAllocationLowWatermark`, `RadialAllocationReserveFraction`, and `RadialMaxFaceTransferMultiplicity` when adaptation is enabled. It generates the corresponding application variables and retains explicit input rather than hidden defaults. |
| `input/halley1d.input` | Sets the agreed controls to `0.50`, `0.10`, and `1`; documents transactional chemistry and the hard 20-daughter allocation formula. This is the canonical tracked Halley input. |
| `srcHartley2/README.md` | Replaces the superseded expected/absent and empty-neighbor descriptions with the hard chemistry, all-face, transactional, and low-water algorithms and their limitations. |
| `src/pic/README.md` | Documents both generic nullable core interfaces, their no-op defaults, callback contracts, MPI-only validation status, and core merge grouping. |
| `srcHartley2/halley2011/PAPER_PARAMETERS.md` | Separates the 50%, 75%, 10%, face-one, and daughter-20 controls from Rubin et al. physics and states the exact hard bounds used in comparisons. |
| `srcHartley2/STATUS.md` | Records the motivating failure, implementation, evidence, limitations, and this master-integration procedure. |

The root `hartley2.pl` and `halley1d.input` currently mirror their canonical
tracked files because AMPS preprocessing stages inputs in the repository root.
They are untracked build inputs, not additional source files, and must not be
included in a later source commit or merge.

### Clean build and focused evidence

The workspace `build/` directory was removed before the final compile. The
exact Halley input then passed a clean optimized build with
`make -j8 EXTRACOMPILEROPTIONS=-O3`. Only the previously known DIM1
aggressive-loop warnings remained. The following standalone checks passed:

- `perl -c input/hartley2.pl`;
- `perl input/test_comet_reaction_network.pl`;
- `perl input/test_comet_physics_database.pl`;
- the C++ reaction-network test; and
- the C++ collision-database test.

The focused validation archive is
`validation_runs/halley1d_allocation_burst_control_20261104`. The exact
`--validate-weight-adaptation` executable passed with 1, 4, and 12 MPI ranks:

```text
1 rank:  owner=1  emergency=1  allocation-burst=1  sampling=1  face=1  photochemical=1  ghost=0
4 ranks: owner=4  emergency=4  allocation-burst=4  sampling=4  face=4  photochemical=4  ghost=28
12 ranks: owner=12 emergency=12 allocation-burst=12 sampling=12 face=12 photochemical=12 ghost=84
```

Initialization changed no population and reduced the maximum hard channel
allocation from 335 to 18 (`particles=0->0`, 746 chemistry-constrained plans).
The synthetic chemistry case forced 3,202 possible allocated products and
reduced that bound to 18. The all-face test reduced a two-sided clone mean from
2 to 1. The allocation emergency deliberately used a compact synthetic low
target and required two common-factor doublings; on one rank utilization fell
from `2.048e-4` to `5.195e-5`, proving that the request is consumed only after
the measured target is reached. Equivalent collective and ghost-copy checks
passed at 4 and 12 ranks.

A final four-rank production smoke is in
`validation_runs/halley1d_allocation_burst_control_20261104/production_smoke_final`.
It was deliberately bounded by an external 60-second timeout (exit 124), not
an AMPS failure. It reached iteration 1,074 and wrote 54 plot files: 27 species
at iterations 500 and 1,000. The production initialization again reduced
335 to 18 without creating particles. The iteration-500 face maximum fell from
512 to 1; at iteration 1,000 the adaptation changed `56,355` particles to
`54,692`, reduced the face maximum from 512 to 1, and retained a hard chemistry
maximum of 12 against the limit of 20. The iteration-1,000 convergence metric
was `0.8150373` for NH3 with `0/3` windows. No buffer-full, MPI, chemistry,
postcondition, or abort message is present. This smoke proves integration
through two output/adaptation boundaries; it is not a convergence or Rubin
agreement result.

### Known limitations and remaining validation

- Transactional batch chemistry is MPI-only. A hybrid/OpenMP run with the
  callbacks installed intentionally aborts. Thread-safe counters/list handling
  and hybrid tests are required before removing that guard.
- DIM2/DIM3 face and subface enumeration is implemented but not validated.
  Add refined coarse/fine interfaces, disconnected components, bidirectional
  crossings, and owner/ghost checks before claiming generic multidimensional
  allocation control.
- The current reachability preparation is conservative across a connected
  component. A future disconnected multidimensional mesh can reduce variance
  by computing chemistry reachability and caps per transport component, but it
  must preserve the hard final postcondition.
- Factor-two roulette is conservative only in expectation and changes sampling
  variance. The new population history and Rubin comparison must be assessed,
  not assumed equivalent to the failed run.
- The final smoke ended at iteration 1,074. A fresh 14-rank run must pass the
  previous iteration-8,136 failure, reach the 10,000 ceiling or the configured
  convergence criterion, and generate a new accepted Rubin comparison.

### Required master-merge sequence after validation

Do **not** merge yet. Remain on `1d_comet_dev`, and tell the user before any
switch, rebase, cherry-pick, integration-branch creation, or merge involving
`master`. After a long run and its comparison are accepted, perform the work
on a disposable integration branch based on the then-current master:

1. Record `git status`, `git diff --name-status master`, the branch tip, exact
   input hash, executable hash, MPI count, and accepted validation paths. Make
   source commits clean before integration; retain run archives outside those
   commits.
2. Merge the generic pre-exchange primitive as one nullable core unit:
   `pic.h` plus `pic_parallel.cpp`. Verify the callback is before both exchange
   implementations and does not recurse. Compile and smoke-test a representative
   non-comet application to prove the null default is inert.
3. Merge the generic photochemistry batch primitive as a second core unit:
   `pic.h` plus `pic_photolytic_reactions.cpp`. Preserve the historical null
   path, all-or-none validation, global remaining-parent reduction, loop guard,
   and explicit hybrid rejection. Compile both MPI and hybrid configurations;
   the MPI null path must run normally and the installed hybrid path must fail
   with the documented message until thread-safe support is implemented.
4. Merge the comet parser/input and application implementation together:
   `input/hartley2.pl`, canonical comet input, `CometChemistry.*`, and
   `RadialSimulation.*`. Do not cherry-pick only the input keywords, only the
   driver callbacks, or only the face planner. Preserve the execution order:
   population plans, hard chemistry plans, face-component enforcement, final
   hard chemistry verification, owner/ghost synchronization, and cache rebuild.
5. Resolve `RadialSimulation.cpp` semantically rather than selecting an entire
   side of a conflict. Preserve `W/dt`, power-of-two-only increases, actual
   `ceil` multiplicities, all reaction products, direct-source roots, both
   directions of every face, empty-cell preservation, and the 90/75/50
   transaction thresholds. Do not replace these with raw `W`, expected values,
   absent-products-only logic, or empty-neighbor-only propagation.
6. Remove `build`, stage the canonical parser and `input/halley1d.input` through
   the normal AMPS configuration, and rebuild from scratch. Repeat parser and
   database checks, focused 1/4/12-rank validation, chemistry ensemble,
   local-timestep face tests, and a non-comet regression.
7. Repeat the full 14-rank Halley run from iteration zero. Require stable
   passage beyond 8,136, either convergence or the 10,000 ceiling, consistent
   population diagnostics, and regenerated Rubin PNG/EPS/PDF/CSV/manifest
   artifacts from the accepted final window. Add DIM3 and refined-face tests
   before advertising the controller as validated for higher dimensions.
8. Stop for user review of the integration diff and validation evidence before
   recording or applying any merge to `master`.

Do not commit or merge `build/`, `amps`, root-staged parser/input copies,
`PT*`, `validation_runs/`, logs, meshes, `SWMF_data`, `share/`, generated
tables, or other run artifacts. The high-conflict core locations are the
parallel declarations/exchange wrapper in `pic.h`, the photolytic namespace in
`pic.h`, the driver loop in `pic_photolytic_reactions.cpp`, and the pre-exchange
definition in `pic_parallel.cpp`. Application conflicts will concentrate in
the adaptation order and initialization callbacks in `RadialSimulation.cpp`.
No commit, branch switch, or merge was performed in this implementation round.

## 2026-11-04: fresh 14-rank allocation-burst validation started

A new validation run was staged from iteration zero after removing the exact
workspace `build/` directory and compiling the Halley input with optimized
settings (`make -j8 EXTRACOMPILEROPTIONS=-O3`). The parser, reaction-network
database, physics database, and `git diff --check` gates passed. The resulting
executable also passed the allocation-burst validation mode on 1, 4, and 12
MPI ranks before the long launch. The focused logs are under
`validation_runs/halley1d_allocation_burst_preflight_20261104_clean/`.

The long-run archive is:

```text
/data/vtenishe/COMET1D/AMPS/validation_runs/halley1d_allocation_burst_convergence_rubin_20261104_14mpi_fresh
```

It was launched at 2026-11-04 12:25:32 CST with 14 MPI processes and one
OpenMP thread per process. The exact staged executable and input SHA-256 values
are recorded in `RUN_MANIFEST.md`. An initial sandboxed MPI attempt failed
before application initialization because local Open MPI networking was
denied; that diagnostic is retained as `error.sandbox-launch.log`. The same
executable then started successfully with host MPI access.

At the 12:27:05 snapshot the run had reached iteration 1,819 and had written
three complete 500-iteration windows. Initialization changed no particles and
reduced the worst chemistry allocation bound from 335 to 18. The adaptations
at iterations 500, 1,000, and 1,500 ended with 29,169, 55,138, and 81,182
particles, respectively. At iteration 1,500 the hard reaction multiplicity
was 12 against the configured limit of 20, and the maximum expected face clone
mean was reduced from 512 to 1. The latest convergence metric was 0.7885894,
limited by NH3, with 0/3 consecutive passing windows; this is an active run,
not a convergence claim. No AMPS error, NaN, buffer exhaustion, hard
postcondition failure, or MPI abort was present. The lone `No protocol
specified` message in `error.log` is from the optional GUI plot viewer in the
headless launch environment and is non-fatal.

No commit, branch switch, master operation, or merge was performed for this
launch. The run subsequently reached the documented 10,000-iteration ceiling,
and the next section records the diagnostic comparison generated from its
final complete window.

## 2026-11-04: completed-run Rubin comparison generated

The fresh 14-rank run finished normally at the 10,000-iteration hard ceiling
with exit status zero. It passed the historical iteration-8,136 buffer failure
and ended with 391,034 particles after the final adaptation. Allocation
stability passed, but convergence did not: window 20 reported maximum relative
L1 1.946230 for CH3 and 0/3 consecutive windows at tolerance 0.05.

The final complete `out=19` profiles were compared with the digitized Rubin
et al. (2011) neutral DSMC curves. The generated directory is:

```text
/data/vtenishe/COMET1D/AMPS/validation_runs/halley1d_allocation_burst_convergence_rubin_20261104_14mpi_fresh/PT/plots/rubin2011_comparison
```

It contains 18 PNG figures, 18 EPS figures, an 18-page combined PDF,
`comparison_profiles.csv`, and `manifest.json`. Both calculated and
published curves are present. The plot generator was updated so every detached
figure visibly carries `NONCONVERGED (0/3 windows)`; the manifest now also
stores `scientific_status: nonconverged` and an explicit diagnostic-only
limitation. This prevents the completed allocation-stability run from being
mistaken for an accepted steady-state science result.

```text
a445b9f23d920836f0f9bdc641bc34e92e44d826f2a5ae6bfe34d7e6aa92a8bf  manifest.json
223880ba79fe2a5d63295c05b1f1d10618da82ac3d520add93179676fd90f5a9  comparison_profiles.csv
ff6c7677a65f5c10a28ec511187590cfaa28e844914063726d5f6b45310dd005  rubin2011_neutral_profiles.pdf
```

The overlay is diagnostic, not a Rubin-agreement acceptance result. The next
step is to review the primary-species density and velocity discrepancies and
separate physical evolution from noise in sparsely sampled daughter species
before changing convergence controls or launching another long run. No commit,
branch switch, master operation, or merge was performed.

## 2026-11-04: master-integration appendix updated after `c76c488`

The allocation-burst development is now committed at
`c76c48820b4d3ee4ce28c2bb816ae44a2cafd641`, which matches
`origin/1d_comet_dev`. The current local `master` and merge base are
`cb7be2bba41f373352ab853abed23955ea69c8b9`. The comet branch contains seven
dependent commits after that base and changes 72 files (+29,224/-437).

Do not merge yet. Before any switch, fetch, rebase, cherry-pick, integration
branch, or merge involving `master`, tell the user. The latest commit cannot
be cherry-picked alone: it depends on the radial mesh, generic parser and
databases, collision/cooling/chemistry wiring, face-stopping local-timestep
mover, block-local weights, and validation framework introduced by the six
earlier commits.

The expanded whole-branch integration procedure is in `handoff.md`, under
"Whole-branch master-integration appendix after `c76c488`." It records:

- the exact seven-commit dependency chain and current history boundary;
- exclusions for the local root `Makefile` change and generated run data;
- semantic phases for parser/databases, dimension-neutral core repairs,
  physics, face-stopping transport, adaptive weights, core callbacks, and
  visualization/provenance;
- high-conflict core/application files and invariants that must survive
  conflict resolution;
- clean build, non-comet/DIM3, parser/database, local-timestep, chemistry,
  1/4/12-rank weight, and fresh 14-rank validation gates; and
- the final user-review stop before any commit or operation changing
  `master`.

The most important integration invariants remain `W/dt` face carriers,
face-first continuation, power-of-two adaptive weights, hard upward-rounded
complete reaction multiplicities, all-face bidirectional enforcement,
empty-cell preservation, owner/ghost synchronization, null generic-core
defaults, and explicit nonconvergence labels in Rubin products.

The completed 14-rank run validates allocation stability through iteration
10,000 but remains nonconverged (relative L1 1.946230 for CH3, 0/3 windows).
This result is evidence for allocation safety, not authority to merge or a
claim of Rubin agreement. No branch switch, integration-branch creation, or
merge was performed while adding these instructions.

## 2026-11-04: priority hot-loop acceleration implemented and validated

The first acceleration phase is implemented on `1d_comet_dev`. No branch
switch, commit, integration-branch operation, or merge was performed. Before
compilation, the exact workspace `build/` directory was removed. The canonical
`halley1d.input` was then preprocessed and the application rebuilt with
`make -j8 EXTRACOMPILEROPTIONS=-O3`. The executable SHA-256 is
`8ffab40d0bc09fc5e12c5a0288bd78efb390780dfd2009295dea94179e55e42d`;
the input SHA-256 remains
`5652c655e61e3d4ecbbdb5c477a60c9f89b493576f77440fa52edd079c6e346e`.

### Modified source files and purpose

- `src/models/cometary_physics/CollisionDatabase.h` now owns a dense symmetric
  species-pair pointer index. `Find()` is O(1), and new `Evaluate()` applies a
  law to the record already found. Construction tolerates malformed records
  only long enough for `Validate()` to return its detailed error; invalid or
  duplicate input is still rejected.
- `src/models/cometary_physics/test_collision_database.cpp` checks direct and
  reversed indexed lookup, invalid indices, and `Evaluate()` equivalence.
- `srcHartley2/CometPhysics.cpp` retains the record returned by `Find()` and
  evaluates that record, eliminating the old second database lookup. The
  temperature floor, 300 K cap, undersampling rule, and fatal missing-pair
  behavior are unchanged.
- `src/models/chemical_solver/ReactionNetwork.h` now builds flat offsets and
  reaction-index arrays grouped by valid unit-unimolecular parent species.
  Channels preserve input order. General binary/higher-order records remain
  available in the main network.
- `src/models/chemical_solver/test_reaction_network.cpp` checks counts, order,
  binary-reaction exclusion, and invalid accessor arguments.
- `srcHartley2/CometChemistry.cpp` builds a per-parent ordered channel/rate
  cache from the generic index. The cache is initialized after validation and
  refreshed at a serial batched-stage boundary if heliocentric distance
  changes. `ProcessParticle()` no longer scans all 69 reactions or constructs
  two temporary vectors. Its occurrence and channel-selection random draws,
  left-to-right rate sum, channel order, energy partition, and transactional
  behavior are unchanged.
- `src/pic/pic_particle_collision_models.cpp` changes only the general NTC
  cell routine used by the local-weight/local-timestep Halley run. Reusable
  `thread_local` species vectors replace per-cell arrays. A single linked-list
  traversal fills all buckets and weight statistics. The fixed 100-probe
  majorant, species/pair order, RNG order, acceptance rules, flag lifecycle,
  and writeback points remain unchanged.

The root `Makefile` modification shown by `git status` predates this work
(`EXTRALINKEROPTIONS_F=-lmpi++`) and was neither created nor altered in this
acceleration phase. Generated executables, meshes, `PT` output, and validation
archives are evidence only and are not source changes for a future merge.

### Validation evidence

The standalone reaction-network and collision-database C++ tests pass. Both
Perl generator/parser suites pass. The application chemistry ensemble passes
20,000 events with all 46 neutral reactions, including 25 unequal
weight/timestep and 23 sink cases. Weight-adaptation validation passes on 1,
4, and 12 MPI ranks; the 12-rank result includes 84 ghost-weight checks plus
all emergency, allocation-burst, sampling-population, face-carrier, and
photochemical checks. Focused artifacts are under:

```text
/data/vtenishe/COMET1D/AMPS/validation_runs/halley1d_priority_optimization_preflight_20261104
```

The `benchmark14` subdirectory is a fresh 14-rank run from iteration zero. It
completed all 10,000 iterations and produced 20 sampling windows. Comparing
every file matching `PT/plots/*.dat` against the prior fresh baseline found
`total=540`, `differing=0`, and `missing_baseline=0`. The first-to-last output
timestamp span is 470 s versus 3,052 s for
`halley1d_allocation_burst_convergence_rubin_20261104_14mpi_fresh`, a measured
6.4936x speedup. Final population is 391,034 in both runs; final relative L1 is
1.946230 for CH3 with 0/3 passing windows in both runs. The byte identity is
the strongest available check that reaction/collision order and RNG
consumption were preserved.

The pre-change clean Callgrind window recorded 3,476,809,796 instructions;
the linear collision lookup alone accounted for about 2.09 billion. The
post-change `callgrind_optimized/callgrind.out.2359084.1` window records
1,791,797,371 instructions for 25 late-run timesteps. Collision processing is
still the largest kernel because NTC requests many cross sections, but the
linear table scan is absent. `pow` for the H2O temperature law accounts for
about 40.1 million instructions (2.2%), so temperature caching was deliberately
not added without a state-invalidation design.

### Master-integration notes for this phase

Do not merge yet, and tell the user before any attempt involving `master`.
When the comet branch is eventually integrated, treat these as three semantic
units after the generic database/core prerequisites are present:

1. Merge `CollisionDatabase.h`, its test, and `CometPhysics.cpp` together.
   Preserve constructor-built symmetric indexing, `Validate()` diagnostics,
   fatal missing-pair policy, and the direct `Find()`/`Evaluate()` sequence.
2. Merge `ReactionNetwork.h`, its test, and `CometChemistry.cpp` together.
   Preserve input order, the general network table, cache refresh at a serial
   stage boundary, and the exact particle RNG sequence.
3. Merge the NTC bucket change independently. Confirm the target master path
   still uses the same `cParticleDataList` layout and writeback rules. Keep the
   storage `thread_local`; do not copy the optimization into specialized NTC
   implementations without their own tests.

After conflict resolution, remove only the integration worktree's `build/`,
rebuild Halley, run both standalone/generator tests, the chemistry ensemble,
and 1/4/12-rank weight matrix, then run a fresh deterministic 14-rank case.
Require all 540 output data files to match the accepted development reference
or document and review every intended physical difference. Reprofile rather
than assuming the same hotspot distribution on newer master code. Keep the
existing whole-branch master stop point and nonconvergence disclosure.

## 2026-11-04: optimized 14-rank validation repeated

The optimized production validation was repeated independently from iteration
zero with 14 MPI ranks and one OpenMP thread per rank. The run used the same
staged executable and input hashes documented above and completed normally at
the 10,000-iteration ceiling with exit status zero. Its archive is:

```text
/data/vtenishe/COMET1D/AMPS/validation_runs/halley1d_priority_optimization_repeat_20261104_14mpi
```

All 20 sampling windows and 540 `.dat` files were produced. Every data file is
byte-for-byte identical to both the first optimized run and the pre-optimization
fresh baseline: zero different and zero missing files in each comparison. The
repeat output timestamp span is 482 s. Final adaptation again gives
391,194 -> 391,034 particles; final relative L1 is again 1.946230 for CH3 with
0/3 passing windows. The only stderr text is the known headless-viewer
`No protocol specified` message. There is no AMPS, MPI, buffer, mover,
chemistry, or allocation failure.

This repetition strengthens deterministic and operational equivalence but
does not change the scientific disposition: the run is nonconverged. No
commit, branch switch, or merge was performed.

### Rubin comparison generated from the repeated optimized run

The final complete output (`out=19`, iteration 10,000) was compared with the
digitized Rubin et al. neutral DSMC profiles. The generated directory is:

```text
/data/vtenishe/COMET1D/AMPS/validation_runs/halley1d_priority_optimization_repeat_20261104_14mpi/PT/plots/rubin2011_comparison
```

It contains 18 PNG files, 18 EPS files, an 18-page PDF,
`comparison_profiles.csv`, and `manifest.json`. A representative H2O velocity
overlay was inspected and contains both AMPS and Rubin DSMC curves. Every
figure visibly states `NONCONVERGED (0/3 windows)`, and the manifest records
`scientific_status: nonconverged` and the iteration-10,000 convergence line.
The PDF SHA-256 is
`d6d7a9b25d4dd62a39bdf46bfd0cd50f1ae8af097da26a191698317cb3e1f35b`.
These products are diagnostic and do not change the failed convergence gate.

### Master-integration documentation extended

`handoff.md` now contains "Master-integration addendum for the optimized
model." It adds the exact uncommitted source/document set and exclusions,
orders the collision index, reaction/rate cache, and NTC bucket changes after
their existing comet prerequisites, lists conflict-resolution invariants, and
adds deterministic 540-file, performance-profile, and regenerated Rubin-plot
gates. The addendum preserves the explicit final stop before any operation
that changes `master`. No merge-related operation was performed while adding
it.

## 2026-11-04: four-times-particle run restarted with 1,000-step sampling

The requested higher-statistics Halley validation configuration now uses
`ReferenceInjectionParticleNumber=32`, a requested block-mean range of
200--800 particles/species/cell, a 20,000-iteration ceiling, and
`RadialSamplingLength=1000`.  The earlier attempt with a 500-step sampling
length was intentionally stopped at iteration 4,752 and preserved as a
superseded run; it did not terminate because of an AMPS failure.

The exact `build/` directory was removed before rebuilding with `-O3`.
Generated-source inspection confirmed 20,000 iterations, sampling length
1,000, range 200--800, and reference injection 32.  The rebuilt executable
SHA-256 is
`cf939ed9b2b8d4a27786ce359f9a64b26383d5ddbb1b42cd558560dcfedb2f79`;
the input SHA-256 is
`26e16d7a7f26b33626d17b6f5a3d212c14c8426862555a27a691829f1d6a9b8f`.

Before restart, the chemistry ensemble passed all 46 neutral reactions with
20,000 events per reaction, and particle-weight adaptation passed at 1, 4,
and 12 MPI ranks for the configured 200--800 range.  The new 14-rank run was
launched at 2026-11-04 15:26 CST from:

```text
/data/vtenishe/COMET1D/AMPS/validation_runs/halley1d_4x_particles_20000iter_sampling1000_20261104_14mpi/production
```

Runtime output explicitly reports `current sample length:1000`.  The run
completed its first sampling window and created all 27 species files before
this status update.  Full provenance and preflight logs are recorded in the
parent archive's `RUN_MANIFEST.md`.  No branch switch, commit, or merge was
performed.

## 2026-11-04: higher-statistics run completed and compared with Rubin

The 14-rank, 1,000-step-sampling run documented immediately above completed
normally at its 20,000-iteration ceiling with exit status zero.  It produced
20 complete windows and all 540 expected `.dat` files.  The final sampled
population was 2,021,443, and final weight adaptation produced 2,024,519
particles, about 10.1% of the 20-million-slot capacity.  No emergency
adaptation or runtime failure occurred.

The run remains scientifically nonconverged: the final maximum relative L1 is
1.521183 for CH with 0/3 passing consecutive windows.  It must not be treated
as an accepted steady Rubin solution.

The final `out=19` data were overlaid with the digitized Rubin DSMC profiles
at:

```text
/data/vtenishe/COMET1D/AMPS/validation_runs/halley1d_4x_particles_20000iter_sampling1000_20261104_14mpi/production/PT/plots/rubin2011_comparison
```

The validated artifact set contains 18 PNG files, 18 EPS files, an 18-page
PDF, `comparison_profiles.csv`, and `manifest.json`.  Both a density panel and
the H2O-velocity panel were visually checked for AMPS/published curve pairs and
the explicit `NONCONVERGED (0/3 windows)` label.  The PDF SHA-256 is
`19d12367eee0a22ae7281c1f10fbef6c671285a296dffeecd6521d05da09e381`.
No commit, branch switch, or merge was performed.

## 2026-11-04: second higher-statistics 18-rank run started

At the user's request, the Halley controls were increased by another factor of
four relative to the completed 20,000-iteration run.  The injection reference
is now 128 and the adaptive block/species target is 800--3,200 particles per
cell.  This is four times the immediately preceding run and sixteen times the
original eight-event/50--200 numerical baseline; physical source rates are
unchanged.  Sampling length is 10,000 and the hard ceiling is 40,000
iterations.

The exact `build/` directory was removed before an `-O3` clean rebuild.
Generated C++ inspection confirmed all four requested controls.  Both Perl
generator suites, the 20,000-event chemistry ensemble for all 46 reactions,
and the 1/4/12-rank weight-adaptation matrix passed.  The weight tests report
the configured 800--3,200 range and include emergency, allocation-burst,
sampling-population, face-carrier, photochemical, and ghost-weight checks.

The fresh run uses 18 MPI processes and was launched at 2026-11-04 16:47 CST
from:

```text
/data/vtenishe/COMET1D/AMPS/validation_runs/halley1d_16x_particles_40000iter_sampling10000_20261104_18mpi/production
```

The executable SHA-256 is
`bf32933f8155745d9623b5b2836f608234daaf32946284916506c649775113fd`
and input SHA-256 is
`3e9d2d70a2b19ffdedc2b4a58f441ef25b4dd9e031a9229d4665696856a0810e`.
Startup output explicitly reports `current sample length:10000`, range
800--3,200, and no error.  Because window 1 establishes the reference, only
iterations 20,000, 30,000, and 40,000 can provide the three comparisons
needed for an earliest possible 3/3 convergence result; every one must pass
for convergence before the ceiling.  Full provenance is in the archive's
`RUN_MANIFEST.md`.  No commit, branch switch, or merge was performed.

## 2026-11-05: 40,000-iteration run completed and Rubin comparison generated

The 18-rank, 10,000-step-sampling validation run completed all 40,000
iterations and exited with status zero after approximately 7 hours 39 minutes.
It produced four complete output windows and all 108 expected species `.dat`
files. The final sampled population was 10,733,211; final regular adaptation
changed 10,733,262 particles to 10,737,306. There was no particle-buffer
exhaustion, emergency adaptation, or AMPS abort.

The run did **not** satisfy the scientific convergence test. Its final record
was `max-relative-L1=3.733287e-01`, limiting species CH, and 0/3 consecutive
passing windows at iteration 40,000, versus a tolerance of 5.0e-2. The final
plots are therefore diagnostic comparisons and are explicitly marked
`NONCONVERGED (0/3 windows)`; they must not be presented as an accepted steady
Rubin solution.

The final `out=3` data were overlaid with the digitized Rubin profiles at:

```text
/data/vtenishe/COMET1D/AMPS/validation_runs/halley1d_16x_particles_40000iter_sampling10000_20261104_18mpi/production/PT/plots/rubin2011_comparison
```

The validated set contains 18 PNG files, 18 EPS files, an 18-page combined
PDF, `comparison_profiles.csv`, and `manifest.json`. The density and
H2O-velocity panels were visually checked for both AMPS and Rubin curves and
the nonconverged qualification. The combined PDF SHA-256 is
`68d66ac57b082ffde4832f0159c20deb19ef2fd594c389931ca58ca871f04080`.

The MPI error log retained TCP transport warnings around 22:38 CST, including
`Connection reset by peer` and an invalid-header report. All 18 ranks
continued, all expected outputs were written, and the process exited zero;
the warnings remain documented for subsequent operational investigation.
Complete run, convergence, artifact, and checksum provenance is in the
archive's `RUN_MANIFEST.md`. No commit, branch switch, or merge was performed.
