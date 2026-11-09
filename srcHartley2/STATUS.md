# Generic comet / Hartley 2 status

Date: 2026-11-09
Branch: `1d_comet_dev`

All paths in this file are relative to the AMPS repository root (the
directory that contains `Config.pl`, `input/`, and `srcHartley2/`), except
explicitly absolute scratch locations such as `/tmp/...`.

## 3I/ATLAS six-case campaign (updated 2026-11-08)

### Decision: six independent input files will not be used

The original layout, in which the generator template `input/i3.input`
produced six separately compiled inputs `input/i3.case1.input` through
`input/i3.case6.input` (one clean build and one executable per observation),
is retired and will not be used. Those files have been removed from the
source tree. Do not regenerate, restore, or build from them, and do not start
new runs from executables that were compiled from them.

All six observations now share one compile-time configuration and one runtime
input:

| File | Role |
| --- | --- |
| `input/i3.conf` | Compile-time capability and layout: seven species, DIM1 spherical backend, block/ghost layout, particle storage and buffer maximum, maximum AMR refinement, collision/reaction topology, and brightness linkage. Contains no observation case. |
| `input/i3.in` | Runtime scientific receipt: common physical and numerical settings once, plus six named sections `case: 1` through `case: 6` holding the observation-dependent values (production rates, surface resolution, solar-corrected reaction rates, geometry, radial velocity, Ly-alpha flux). |

Build once and select the case at launch:

```sh
./Config.pl -application=i3
rm -rf -- build
make -j8 EXTRACOMPILEROPTIONS=-O3
mpiexec -n 12 ./amps -i input/i3.in -case N        # N = 1..6
```

`Config.pl` selects `input/i3.conf` because no `input/i3.input` exists.
Changing the case or any runtime value requires no rebuild. `-i` and `-case`
are mandatory for the I3 executable, and an unknown case is fatal. Run each
case in its own directory. The implementation and validation record is the
section "I3 two-layer input and generic case sections (2026-11-08)" at the
end of this file.

Earlier sections that mention `i3.input`, `i3.caseN.input`, or a per-case
rebuild are historical. Their run paths and hashes still record what was
actually executed, but they do not describe the current workflow.

### Physics and data content (unchanged by the input redesign)

The generic DIM1 application has six deterministic, observation-driven
3I/ATLAS cases. `srcHartley2/i3/prepare_i3_cases.py` combines reviewed CSV
catalogs into the single runtime file `input/i3.in`. No 3I production rate,
date, solar index, nucleus radius, or mesh value was added to C++. Published
H/OH fluorescence data are in a reusable core module, not in the application.

The implementation uses:

- the table Q(H2O), Q(CO), and Q(CO2) values as molecules/s; raster-fit H2O
  values are retained only as comparison metadata;
- compiled species `H,H2,O,OH,H2O,CO,CO2`, with CO and CO2 collision-only;
- Moore--Mapaye per-channel `A0*F10.7P^alpha` H2O/OH/H2 chemistry plus neutral
  ionization sinks, scaled once by the existing heliocentric inverse-square
  law;
- state-distinct OH reaction records even though O(3P), O(1D), and O(1S)
  products share the compiled O species;
- a 2.5-km, 175-K uniform spherical source, cometocentric outer radius 1e9 m,
  2,000 radial cells, local timestep/local block weight, and convergence and
  emergency-population controls inherited from the validated Halley policy;
- direct-parent target range 800--3,200 particles/cell, chemical-only range
  50--200, direct-source weight locking, and a derived injection reference of
  114; and
- input-driven H Ly-alpha and OH(0,0) brightness on 200 logarithmic impact
  parameters, with the H thin result guarded by a line-center tau diagnostic.

`reference/solar_activity_audit.csv` records the approved JPL/Carrington
solar-longitude correction. The final exposure-averaged F10.7P values are
`136.194566908`, `142.874791643`, `139.711261811`, `143.568741662`,
`150.444624050`, and `138.071807780` sfu; provisional unshifted values are not
used. `reference/photochemistry_channels.csv` owns all rate coefficients,
solar exponents, energies, and citations. The H2O/H2 excess energies follow
Huebner--Mukherjee through Carberry Mogan et al.; the state-labeled OH values
are 1.27 eV for O(3P), 7.73 eV for O(1D), and 10.0 eV for O(1S), with the
specific published source recorded per row.

The surface cell width is calculated rather than chosen. For each parent,
the generator forms the residence-weighted outgoing half-Maxwellian density
and evaluates `lambda_i=1/(sqrt(2)*sum_j n_j sigma_ij)` with the exact selected
175-K parent collision laws. The minimum path is the CO2 value in all cases:
1.75261509304, 1.29785287642, 0.531312387076, 0.131442336866,
0.122734766819, and 0.131442336866 m.

Focused deterministic evidence recorded on 2026-11-07 for the retired
six-file layout (the current generator test runs 6 tests and regenerates only
`input/i3.in`; see the final section):

```text
python3 -m unittest srcHartley2/i3/tests/test_prepare_i3_cases.py
Ran 5 tests: OK

make -C src/models/resonance_fluorescence test
Ran 3 reference-data tests: OK
PASS resonance-fluorescence core evaluator

srcHartley2/tests/run_comet_brightness_tests.sh
PASS comet-brightness focused test group tests=7
```

The tests check all six corrected solar values, 17 reactions, minimum mixture
paths, direct injection count, recoil-derived characteristic speeds,
brightness records, and byte-for-byte regeneration of every input that was
tracked at that time.
The core fluorescence tests check all nine official IKL-2020 H coefficients,
all 1,201 non-extrapolated Lowell OH records over +/-60 km/s, authoritative
and generated-file SHA-256 hashes, strict model domains, inverse-square
scaling, and agreement with all six existing rounded case scalars. A clean
case-1 application compile also passes. These checks do not replace a
convergence run. Before any later recompile, remove exactly
`build`; preserve the pre-existing root
`Makefile` linker edit and all run archives.

Files added for this work and their reason are:

| File | Purpose |
| --- | --- |
| `input/i3.conf` | Compile-time capability/layout shared by all six cases (replaces the retired template). |
| `input/i3.in` | Single loadable runtime receipt with common settings and six named case sections (replaces the six per-case inputs). |
| `input/i3.input`, `input/i3.case1.input` ... `input/i3.case6.input` | **Retired; will not be used.** Former generator template and six per-case compile-time inputs. Deleted from the source tree; do not regenerate or restore. |
| `srcHartley2/i3/prepare_i3_cases.py` | Network-free calculation of rates, mean free paths, injection statistics, recoil speed bounds, and brightness records; writes only `input/i3.in`. |
| `srcHartley2/CometRuntimeInput.h/.cpp`, `src/general/input/amps_input.h` (`ResolveCase`), `Config.pl` | Strict runtime loading of `i3.in`, generic named-case overlay, MPI-consistent receipt check, and `.input`-then-`.conf` application selection (see the final section). |
| `srcHartley2/i3/reference/i3_cases.csv` | Archive transcription, solar-corrected indices, and fluorescence metadata. |
| `srcHartley2/i3/reference/solar_activity_audit.csv` | Endpoint/exposure longitude shifts and unshifted controls. |
| `srcHartley2/i3/reference/photochemistry_channels.csv` | Geometry-independent chemical network parameters and provenance. |
| `srcHartley2/i3/tests/test_prepare_i3_cases.py` | Determinism and physics-control regression. |
| `srcHartley2/i3/README.md` | Full derivation, citations, limitations, operation, and validation boundary. |
| `src/models/resonance_fluorescence/resonance_fluorescence.h/.cpp` | Application-independent SI API and strict IKL H / Lowell OH evaluators; no MPI, PIC, mesh, species index, or campaign dependency. |
| `src/models/resonance_fluorescence/reference/g_factors/` | Authoritative IKL coefficient source, exact Lowell REST snapshot, normalized OH energy table, and machine-readable provenance/hashes. |
| `src/models/resonance_fluorescence/generated/reference_tables.inc` | Deterministically generated immutable runtime representation; avoids file I/O in model loops. |
| `src/models/resonance_fluorescence/tools/` | Guarded upstream acquisition and deterministic reference-table generator with stale-file check. |
| `src/models/resonance_fluorescence/tests/` | Hash/grid/no-extrapolation audit plus six-case scientific evaluator regression. |
| `src/models/resonance_fluorescence/makefile`, root `Makefile` wiring | Builds the reusable object into `libAMPS.a`; preserves the pre-existing root linker-option change. |
| `srcHartley2/CometBrightnessMath.h` | PIC-independent exact shell chord, opaque-nucleus occultation, Rayleigh conversion, non-attenuating tau diagnostic, and inverse of the core DIM1 local-radial velocity rotation. |
| `srcHartley2/CometBrightness.h/.cpp` | Default-off particle-resolved external sampler, named core-model resolution, velocity-dependent per-particle g factors, MPI reduction, synchronized buffer lifecycle, and per-species Tecplot output. Fixed scalar inputs remain supported. |
| `input/CometBrightnessConfig.pm` | Strict fixed-or-named model grammar; requires explicit observation geometry, supports scalar-radial or full-vector comet velocity, enforces an out-of-domain policy, and retains H solar-flux/OH wavelength requirements. |
| `input/hartley2.pl` | Generates model selectors, geometry/vector/policy data, SI case conditions, fixed/model arrays, and the DIM1 link dependency. Its exact declarations and arrays are regression-tested. |
| `srcHartley2/RadialMain.cpp` | Validates brightness configuration, parses `-i`/`-case` and `--validate-runtime-input`, and provides the production surface-injection distribution regression. |
| `srcHartley2/RadialSimulation.cpp` | Registers the brightness sampling and output callbacks after PIC initialization; core output controls reduction, writing, and buffer clearing at the ordinary sampling boundary. |
| `srcHartley2/tests/run_comet_brightness_tests.sh` and brightness tests | Velocity-rotation/geometry/unit/policy/parser regression, including exact generated declarations and array values. |
| `srcHartley2/tests/run_i3_surface_injection_test.sh` and `I3_SURFACE_INJECTION.md` | Reproducible H2O/CO/CO2 flux-distribution runtime test and its mathematical contract. |
| `input/halley1d.input` | Explicitly selects `BrightnessMode=off`; no Rubin comparison physics is changed. |
| `srcHartley2/CometChemistry.cpp` | Treats equal carrier rates as a valid network while continuing to report unequal-carrier coverage for dedicated fixtures. |

### Clean case-1 build and focused runtime evidence (historical, retired layout)

This subsection records the 2026-11-07 build made from the retired
compile-time file `input/i3.case1.input`. It is kept as evidence only. The
current equivalent is the clean `i3.conf` build and the
`-i input/i3.in -case N --validate-runtime-input` checks recorded in the final
section.

The latest case-1 source was configured only after removing exactly
`build`, then compiled and linked with C++17/MPI.
The first attempted configure exposed that `.amps.conf` overrides
the Make command's `InputFileAMPS`; it was stopped while still compiling the
old Halley selection. The workspace selector was set to
`input/i3.case1.input`, `build/` was removed again, and every accepted build
reported that exact input. This configuration produced 376 AMR blocks and
3,008 allocated dyadic cells; the requested ideal radial profile remains
2,000 cells, and the finest allocated cell is 0.9313202 m versus the
1.75261509304-m case-1 surface request.

An earlier brightness-enabled configure also found a real generator defect:
the five `ChangeValueOfArray` regular expressions in `input/hartley2.pl`
used one Perl escape where two were required. Metadata and `Enabled=true`
were correct, but the generated numerical arrays remained zero. The parser
was fixed and `test_comet_brightness_parser.py` now asserts the exact species,
g-factor/model, wavelength, tau-mode, cross-section, and threshold array
contents. The latest generated case-1 selectors are OH index 3 with model 2,
effective wavelength 308.5 nm, and H index 0 with model 1. The generated
bulk velocity is -39.7 km/s and the solar Ly-alpha flux is
`5.1040518742e11 photon cm^-2 s^-1` at 1 AU. Runtime startup resolves the
case-1 photon rates to about `1.800241179e-4 s^-1` for OH and
`9.658586e-4 s^-1` for H; the small differences from archived rounded scalars
are covered by the core regression. H keeps `sigma_0=9.5e-18 m^2` and tau
threshold 0.1.

Focused runtime results against the final executable are:

```text
PASS I3 surface-injection validation mpi-ranks=1
     species=H2O,CO,CO2 temperature-K=175 particles-per-species=200000
PASS chemistry reaction=... for all 8 neutral-product channels
Chemistry ensemble validation result: PASS
     neutral-reactions=8 unequal-weight/timestep-reactions=0 sink-reactions=9
```

The zero unequal-carrier count is intentional for this input: H2O is the
largest direct-source `W/dt` carrier and every reacting H2O descendant
inherits it, so residual multiplicity is exactly one. `CometChemistry.cpp`
no longer declares an otherwise valid chemistry configuration failed merely
because its configured carriers are equal. It continues to report the count;
the independent Halley Stage-4 log checker still requires a positive unequal
count and therefore retains coverage of split/roulette behavior.

Final local build evidence (build products are not source and must not be
merged) is:

```text
9fd8a455d3d33d8858ab25d038bcf3c9742c6ebbb7857eeabf2a693f11674c9b  amps
3a8b9df77c9f15654b95bb8b5d33021b5c5e1722cee6f2cf6c5f7e52c4f2dfe7  input/i3.case1.input
d3820a51eec9077c99b78180c6581ac1a4e026f749226442b43f2099eaead104  build/main/CometBrightness.cpp
1fb36028af013311e14d9a7140fbdb08fba854c77efb11ec873a804869d6089d  src/models/resonance_fluorescence/generated/reference_tables.inc
068380d305aa16e276f16e36780c8b8204c988c99eb2bc5d0f7b1f946edb4662  build/main/CometChemistry.generated.h
2b45b41de73d6a98e1231c3d109af8adef348db99101fed92653e575f4587dd4  build/main/CometPhysics.generated.h
9ddabb43c544acaebb08f70732e417c25ecba4a47144e52702d98beb8c56282e  build/pic/picGlobal.dfn
```

### Current I3 status and remaining work (2026-11-08)

The I3 work is not yet a scientific validation result.

- Configuration: one clean optimized executable built from `input/i3.conf`
  loads all six cases from `input/i3.in`. Cases 1--6 each pass
  `--validate-runtime-input` on one rank, and case 1 also passes on two ranks.
  The surface-injection regression passes through `-i/-case` with 200,000
  samples per species for H2O, CO, and CO2.
- Production: all six cases ran with the `i3.conf` executable and
  `-i i3.in -case N` on 128 ranks, under
  `validation_runs/i3_campaign_20261108/caseN_np128/production`. Every run
  ended cleanly at the 30,000-iteration cap without converging (window-3
  change 28--37% against a 5% tolerance). See "I3 six-case production results
  (2026-11-08)" at the end of this file.
- Done in this round: 1/4/12-rank brightness equivalence, the `ABORT` and
  `EXCLUDE_AND_REPORT` out-of-domain fixture, sampler performance and peak
  memory, the H tau-flag review, and nonconverged MAVEN/IUVS comparisons for
  outputs 0--2.
- Superseded: the earlier 12-rank case-1 run in
  `validation_runs/i3_case1_rays200_sampling10000_max30000_indexfix_20261107_12mpi/production`
  was compiled from the retired `i3.case1.input` and is diagnostic only.

Remaining before scientific use:

- converged runs. Campaign 2 raised `simulation.max_iterations` to 150,000,
  and cases 1, 2, and 5 converged at iteration 70,000. See "I3 campaign 2
  interim results" at the end of this file;
- the H points beyond, or close to, the domain boundary in cases 4--6.
  Campaign 2 moved the boundary to 3e9 m and corrects the remaining 13--21%
  column loss with the collisionless reference;
- tracked executable-level I3 tests. The portable orchestrator
  `srcHartley2/tests/run_all_stage_tests.sh` (2026-11-09) covers the I3
  generator, MAVEN/IUVS regeneration, and the focused brightness tests. The
  checks that need the I3 executable still exist only as the untracked
  `validation_runs/i3_campaign_20261108/scripts/focused_tests.sbatch` and
  `validation_runs/i3_campaign2_20261108/scripts/smoke.sbatch`. These cover
  runtime input per case, the local-time-step, weight-adaptation, and
  chemistry validators, surface injection, and 1/4-rank brightness
  equivalence;
- comparisons that disclose that the observations are raster digitizations
  without the original numerical uncertainties.

All six cases select the named core H/OH models, evaluated with per-particle
Doppler velocities; the legacy fixed-scalar input form remains
regression-tested. The explicit 1e9-m boundary truncates archive H points
beyond 10^6 km and must be disclosed rather than filled by an undocumented
tail.

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
the intended launch directory is the AMPS repository root. A compile does
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

The final commands, run from that scratch directory with `AMPS_ROOT` set to
the AMPS repository root, were:

```sh
mpirun -np 1 "$AMPS_ROOT/amps" \
  --validate-local-timestep --validation-particles 200000 \
  --validation-mover both --validation-skip-multiface

mpirun -np 4 "$AMPS_ROOT/amps" \
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

Before compilation, `build` was deleted. The exact
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
validation_runs/halley1d_emergency_weight_convergence_rubin_20261103_14mpi
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
validation_runs/halley1d_sampling_emergency_convergence_rubin_20261103_14mpi
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
validation_runs/halley1d_empty_neighbor_weight_convergence_rubin_20261103_14mpi
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
validation_runs/halley1d_allocation_burst_convergence_rubin_20261104_14mpi_fresh
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
validation_runs/halley1d_allocation_burst_convergence_rubin_20261104_14mpi_fresh/PT/plots/rubin2011_comparison
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
validation_runs/halley1d_priority_optimization_preflight_20261104
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
validation_runs/halley1d_priority_optimization_repeat_20261104_14mpi
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
validation_runs/halley1d_priority_optimization_repeat_20261104_14mpi/PT/plots/rubin2011_comparison
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
validation_runs/halley1d_4x_particles_20000iter_sampling1000_20261104_14mpi/production
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
validation_runs/halley1d_4x_particles_20000iter_sampling1000_20261104_14mpi/production/PT/plots/rubin2011_comparison
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
validation_runs/halley1d_16x_particles_40000iter_sampling10000_20261104_18mpi/production
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
validation_runs/halley1d_16x_particles_40000iter_sampling10000_20261104_18mpi/production/PT/plots/rubin2011_comparison
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

## 2026-11-05: minor-daughter chemistry roadmap documented

`srcHartley2/CHEMISTRY.md` now records the measured C2 and other
minor-daughter sampling deficiencies from the completed 40,000-iteration run
and a decision-gated implementation roadmap. The stages cover baseline
diagnostics, input-selectable parent/daughter policies, component-wide
correction-free carriers, residual accumulators, individual-weight support for
the local-timestep mover and NTC collisions, exact weighted daughters, packet
splitting, channel importance sampling, optional total-reaction forcing, and
long-run science validation.

Each stage includes implementation tasks, conservation and capacity
invariants, tests, acceptance criteria, documentation requirements, and a
separate master-merge preparation sequence. The first recommended development
path retains correction mode off and locks direct-source parent carriers;
individual weights are a later decision gate only if the lower-risk method is
insufficient. `srcHartley2/README.md` links to the roadmap. No model source,
input value, branch, commit, or master state was changed while preparing it.

## 2026-11-05: Stage 0 chemistry closure tests implemented

The focused tests required by `srcHartley2/CHEMISTRY.md` Stage 0 are now implemented on
`1d_comet_dev`. No branch switch, merge, or master modification was made.

New and modified files are:

| File | Purpose and later-merge dependency |
| --- | --- |
| `src/models/chemical_solver/ChemistryDiagnostics.h` | New geometry-, PIC-, and MPI-independent event accumulator. It retains reaction-wide eligible/selection/destruction/sink/deferral/emergency counts and product-term zero/one/many multiplicity, allocated/reused carrier, number, mass, momentum, and energy values. Repeated product terms remain distinct. Disabled configuration retains no table and serializes empty vectors. Merge this generic core header before application instrumentation. |
| `src/models/chemical_solver/README.md` | Documents the new generic accumulator boundary, product-term indexing, reduction contract, disabled behavior, and Stage 0 test command. Merge with the core header. |
| `srcHartley2/tests/test_stage0_chemistry_diagnostics.cpp` | New deterministic unit/MPI test. It exercises one, two, repeated, and zero products and distributes the same 240 physical events across 1, 4, and 12 ranks. |
| `srcHartley2/tests/test_stage0_minor_species.py` | New analyzer regression. It checks the exact core L1 formula with explicit measures, missing-measure behavior, comparison ratios, convergence metadata, artifact hashes, and final archived C2 values. |
| `srcHartley2/tests/run_stage0_tests.sh` | New single entry point that compiles the C++ test, compares all three MPI signatures, and runs Python tests. |
| `srcHartley2/visualization/analyze_minor_species.py` | New read-only CSV/JSON analyzer for occupancy, populated radius, `W`, `dt`, `W/dt`, adjacent-window changes, convergence state, and available paper comparisons. |
| `srcHartley2/halley2011/reference/halley1d_stage0_baseline.json` | New immutable machine-readable final-window values and SHA-256 provenance for the external 18-rank archive. It does not claim that the run converged. |
| `srcHartley2/CHEMISTRY.md` | Records the implemented tests and explicitly lists the production instrumentation still required before Stage 0 can close. |
| `srcHartley2/README.md`, `srcHartley2/visualization/README.md` | Add commands, definitions, exact-measure policy, archive behavior, and closure limitations. |

Validation performed in the working tree:

- the counter executable compiled with the runner's C++11 warning flags and
  `-DOMPI_SKIP_MPICXX=1` (the test uses the MPI C API);
- 1-, 4-, and 12-rank executions all printed exactly
  `counts=55900 moments=5.14068525000000000e+07 events=240`;
- all three Python tests passed, including verification of every baseline
  artifact checksum and re-analysis of `pic.C2.s=14.out=3.dat`;
- the C2 regression reproduced 3,010 rows, mean occupancy
  `0.14003317275747507`, 2,154 zero rows, first positive radius
  `6,384,136 m`, density `681143.6827560823 m^-3` at 10,000 km, and
  AMPS/Rubin ratio `0.008362717017089849`;
- the baseline remains explicitly `nonconverged` with 0/3 passing windows.

The legacy Tecplot files do not contain the cell `Measure` used by
`CollectIntegratedDensityProfile()`. The analyzer therefore reports the exact
integrated L1 metric as unavailable for that archive. Tests supply explicit
measures and verify the same `n*V` formula as `RadialSimulation.cpp`; no radius
reconstruction or inferred shell volume is accepted.

### Master-merge notes for this test increment

After comet-model validation is complete and the user explicitly authorizes a
master integration attempt, merge this increment in dependency order:

1. bring in `ChemistryDiagnostics.h` as an isolated generic-core commit and
   run its standalone 1/4/12-rank test before connecting any application;
2. merge the analyzer and its Python test with the committed baseline JSON;
3. merge the runner and documentation after adjusting only repository-relative
   paths, never the frozen hashes or acceptance values;
4. connect production instrumentation in a later commit, with diagnostics off
   by default, then prove byte-identical model output for the disabled path;
5. add exact cell-measure and species flux reports and repeat C2 reconciliation
   before marking Stage 0 complete.

The external `validation_runs/...` directory is evidence, not source to be
merged into master. If it is absent in another checkout, set
`HARTLEY2_STAGE0_BASELINE_ROOT`; otherwise the test runner reports the archive
test as skipped. Do not replace that skip with fabricated fixtures or inferred
cell measures. The production hooks, boundary/source fluxes, full C2 path
reconciliation, and measured real-loop diagnostic overhead remain open.

## 2026-11-05: Stage 0 production diagnostics closed

The previously open production-observability work is now implemented and
validated on `1d_comet_dev`. This section supersedes the open items in the
preceding incremental test note. No branch switch, commit, merge, or master
operation was performed.

### Files added or modified

| File | Change, reason, and merge dependency |
| --- | --- |
| `src/models/chemical_solver/ChemistryDiagnostics.h` | Extended the generic, geometry-free accumulator to retain expected and represented number/mass/momentum/energy plus represented-number variance. Added stable flattened count/moment views and sampling-window reset. Merge this core header first and run its standalone MPI test before application hooks. |
| `src/models/chemical_solver/README.md` | Documents the accumulator layout, variance, reset, and disabled contract. Merge with the header. |
| `srcHartley2/CometDiagnostics.h`, `srcHartley2/CometDiagnostics.cpp` | New comet adapter. It builds counters from the generated arbitrary reaction network, tallies per-species physical rates, reduces them at complete sampling breaks, and writes aligned CSV/JSON. The default is off. It does not draw random numbers or change particles. Merge only after the core header and reaction-network API exist. |
| `srcHartley2/CometChemistry.cpp` | Adds one eligible visit per positive-rate channel; transactional selection, sink, deferral, and emergency counts; expected/represented product moments and variance; and species production/destruction tallies. Existing selection order and RNG order are preserved. Conflict resolution must retain the optimized parent-channel index and transactional allocation logic. |
| `srcHartley2/RadialSimulation.cpp` | Adds `W/dt` surface and boundary tallies, distinct gross inner incident/re-emission and outer escape, exact AMR `Measure()` audit output, and the diagnostic report call after each full sampling window. Keep the same traversal used by `CollectIntegratedDensityProfile()` when merging. |
| `input/hartley2.pl`, repository-root `hartley2.pl` symlink target | Adds validated `CometDiagnostics=ON|OFF` generation and the `CometDiagnostics.o` link dependency for the 1-D application. Default is off. Regenerate, do not hand-edit, generated headers. |
| `srcHartley2/makefile` | Links the new diagnostics adapter. Merge with parser-generated object selection so it is included only where the comet application uses it. |
| `input/halley1d.input` | Explicitly keeps `CometDiagnostics=off` in the canonical 128-injection, 800--3,200 population, 10,000-sample, 40,000-iteration case. The temporary closure run used a staged copy with diagnostics on and smaller controls; do not merge those temporary values into this input. |
| `srcHartley2/tests/test_stage0_chemistry_diagnostics.cpp` | Updated deterministic accumulator test for expected/represented moments and variance. Current 1/4/12-rank signature is `counts=55900 moments=2.17583523750000000e+08 events=240`. |
| `srcHartley2/tests/test_stage0_production_report.py` | New independent archive checker. It validates exact measures, CSV/JSON alignment, C2 channel lineage and actual carrier creation, product slots, `max(5 sigma,1%)` acceptance, flux identities, and explicit zero loss/escape observations. |
| `srcHartley2/tests/run_stage0_tests.sh` | Runs the focused suite and optionally the production archive when `HARTLEY2_STAGE0_PRODUCTION_ROOT` is supplied. |
| `srcHartley2/CHEMISTRY.md` | Moved the roadmap under the application, records Stage 0 closure, and freezes the approved Stage 1 selector decision: support explicit species and derived groups with precedence `explicit species > derived group > global`. |
| `srcHartley2/README.md`, `srcHartley2/visualization/README.md`, `srcHartley2/STATUS.md`, `handoff.md` | Document definitions, units, commands, evidence, limitations, and master-integration order. Reconcile these after code conflicts rather than taking one branch wholesale. |

### Definitions frozen by validation

- An eligible visit is one parent visit to one enabled channel whose evaluated
  rate is positive. Zero-rate channels are not eligible.
- Production agreement is
  `abs(represented-expected) <= max(5*sigma, 0.01*abs(expected))`.
- All species source/loss/boundary fields are sampling-window average physical
  number rates: accumulate represented `W/dt` for each event, reduce across
  ranks, and divide by the number of iterations in the completed window.
- Inner incident and thermal re-emission are separate gross terms; their
  equality for the current boundary is checked instead of silently canceling
  them before output.
- Emergency requests are global per sampling window. A collective request is
  counted once on rank zero before reduction, not once per MPI process.

### Closure evidence

The production archive is:

```text
validation_runs/stage0_closure_20261105/diagnostics_on_10000iter_14mpi/production
```

It used 14 MPI ranks, diagnostics enabled, 500-iteration sampling windows,
and a 10,000-iteration hard limit. The temporary input SHA-256 is
`06369c95d71f4c98b61dba5bce4912f886c2d9023d2d53c6e99cb247ce336d43` and
the executable SHA-256 is
`03099b09674e3077dcdd9415937223477d5a890f48acb1def8d51fe1ff0939ce`.
The run exited zero in 8:21.93 wall time, produced 20 diagnostic windows and
2,984 exact cell measures, observed selected C2 production and a represented
C2 carrier, and triggered no emergency adaptations. C2 destruction and outer
escape were zero during the finite transient; their explicit paths and
nonnegative zero fluxes reconcile rather than being inferred.

After restoring the canonical high-statistics input and removing the complete
workspace `build/` directory, a clean `-O3` build succeeded. The resulting
diagnostics-off executable SHA-256 is
`ce077706cb4aff32af92b339aae712252643b26c5407c8b9efe88d0f4f290ef6`; the
canonical input SHA-256 is
`cbaf3bb629199dc0f3779c7a0e0b5a4b88f5fd2d5677885ad028666aff76d37e`.

All 540 model `.dat` files are byte-identical to the preserved diagnostics-off
baseline in
`validation_runs/halley1d_priority_optimization_repeat_20261104_14mpi/PT/plots`.
The only directory-level difference is that the older baseline also contains
its `rubin2011_comparison` subdirectory. The closure run remained physically
nonconverged (`0/3`, limiting CH3 relative L1 `1.946230`), so this evidence
closes observability and non-perturbation only.

The controlled four-rank 200-iteration pair under
`validation_runs/stage0_closure_20261105` measured:

| Configuration | Wall time | Maximum RSS |
| --- | ---: | ---: |
| diagnostics off | 8.58 s | 1,376,232 KiB |
| diagnostics on | 8.66 s | 1,376,572 KiB |
| difference | +0.93% | +340 KiB (+0.025%) |

### Master-merge instructions for the closed Stage 0 increment

Do not merge yet. After the comet model is validated and the user explicitly
authorizes a master integration attempt:

1. integrate `ChemistryDiagnostics.h` and its README/test alone; require the
   identical 1/4/12-rank signature;
2. integrate the parser switch and diagnostics adapter with the existing
   reaction-network API, keeping `CometDiagnostics=off` as the generated and
   canonical default;
3. integrate chemistry hooks while auditing parent-channel order, occurrence
   and channel random draws, transactional allocation, and emergency request
   placement line by line;
4. integrate radial flux and exact-measure hooks while preserving local
   timestep/weight transport, owner/ghost synchronization, convergence order,
   and the exact AMR traversal;
5. run both Perl generator suites, the focused Stage 0 runner, and an
   instrumented production report check;
6. perform paired diagnostics-off/on runs from the same clean executable and
   input controls; require zero differing `.dat` files with diagnostics off
   versus the pre-integration reference and no unexplained enabled-path model
   differences;
7. report the code commit range, generated-file policy, full diff, tests,
   memory/time overhead, C2 path result, and nonconvergence status, then stop
   for user approval before any master commit or merge.

Never merge `validation_runs`, `PT`, generated executables, generated input
symlinks, meshes, or profiler output. Preserve the archive as external
evidence. Do not turn an absent archive into fabricated fixtures or inferred
cell measures.

## 2026-11-05: Stage 1 input-selectable population policies complete

Stage 1 of `CHEMISTRY.md` is implemented and its focused gates pass on
`1d_comet_dev`. No branch switch, commit, rebase, cherry-pick, integration
branch, or operation on `master` was performed. This stage separates numerical
population policy from comet identity; it does not claim that the existing
40,000-iteration Rubin comparison converged.

### Source files added or modified

| File | Change and reason | Later master-merge dependency |
| --- | --- | --- |
| `input/CometPopulationPolicy.pm` | **New.** Parses global, `DIRECT_SOURCE`, `CHEMICAL_ONLY`, and explicit-species range/lock records; validates conflicts and resolves `species > group > global` deterministically. Policy source/class values are numeric generated-data contracts, not Halley-specific tables. | Add before `hartley2.pl`; preserve exact numeric enum agreement with `RadialSimulation.h`. |
| `input/test_comet_population_policy.pl` | **New.** Exercises precedence, transitive reaction classes, explicit overrides, locks, and rejection of duplicate, unknown, contradictory, invalid, and global-lock records. | Merge with the module so parser behavior cannot land untested. |
| `input/CometReactionNetwork.pm` | Adds `reachable_species_from()`, a generic closure over enabled positive-rate reactions. This derives later-generation chemical-only species from arbitrary networks. | Reconcile with any master reaction-parser changes first; do not replace the closure with a fixed comet species list. |
| `input/hartley2.pl` | Delegates population records to the policy module, reads the independent adaptation interval, derives direct-source and reaction-reachable groups after all databases are loaded, and generates species-indexed ranges, locks, policy sources, and classes. | Merge after both Perl modules and before C++ consumers. Regenerate build output; do not hand-merge generated headers. |
| `input/halley1d.input` | Canonical Stage 1 policy: global 800--3,200, chemical-only 50--200, direct-source lock, and 10,000-iteration adaptation interval. | Merge only after parser/runtime support. Keep Hartley inputs separate; these values are Halley numerical controls, not universal defaults. |
| `srcHartley2/RadialSimulation.h` | Declares generated species policy arrays and stable policy source/class codes; documents legacy cadence and emergency-lock override semantics. | Must remain numerically consistent with `CometPopulationPolicy.pm`. |
| `srcHartley2/RadialSimulation.cpp` | Uses per-species targets, skips locked species in regular and chemistry planning, enforces locks transactionally through the face solver, schedules adaptation independently of sampling, reports all resolved policies, and bases sampling emergency thresholds on resolved maxima. The chemistry allocator now ignores irreducible one-carrier products and may refine an unlocked parent when locked products cannot be coarsened. | Merge after generated arrays. Audit the chemistry constraint and face-carrier transaction together; separating them can reintroduce allocation bursts or silently alter locked source weights. Emergency common-factor coarsening intentionally remains a safety override. |
| `srcHartley2/tests/test_stage1_population_policy_runtime.py` | Parses startup policy/cadence records, checks derived classes, verifies the alternate daughter target leaves every direct-source rate/carrier unchanged, and checks adaptation iterations 0/50/100. | Merge with runtime logging; field names are a tested diagnostic interface. |
| `srcHartley2/tests/test_stage1_legacy_equivalence.py` | Hash-compares the complete set of top-level output `.dat` files from an archived legacy input and a candidate run. | Keep archive paths external and configurable; do not add generated outputs to Git. |
| `srcHartley2/tests/run_stage1_tests.sh` | Runs mandatory Perl/syntax tests and optional archived runtime, cadence, and legacy gates selected by environment variables. It disables Python bytecode output. | Merge after all Stage 1 tests; an absent external archive is an explicit skip, never fabricated evidence. |
| `srcHartley2/CHEMISTRY.md`, `srcHartley2/README.md`, `srcHartley2/STATUS.md`, `handoff.md` | Record the implemented grammar, algorithms, validation evidence, limitations, recovery steps, and merge order. | Reconcile prose after code conflicts; do not take an older whole-file version that reopens Stage 0/1 or loses the master stop. |

The repository-root `hartley2.pl`, `amps`, `build/`, validation executables,
run logs, and `PT` output are generated or external evidence. They are not
Stage 1 source files and must not be staged for the later master integration.
The unrelated dirty top-level `Makefile` remains user-owned and was not
modified for Stage 1.

### Validation evidence

The archive root is:

```text
validation_runs/stage1_population_policy_20261105
```

- `canonical_policy/run.log`: four-rank runtime validator, 27 resolved
  species, fixed cadence 10,000, and a final `PASS` for all 4 owner blocks,
  emergency/allocation/sampling/face/photochemical checks, and 28 ghost-weight
  checks.
- `alternate_policy/run.log`: chemical-only target 100--400. Every locked
  direct-source physical rate and surface carrier matches the canonical
  50--200 build within the test's 5e-14 relative tolerance.
- `cadence/run.log`: short four-rank run with interval 50 and sampling length
  100; adaptation records occur at iterations 0, 50, and 100.
- `legacy_candidate/PT/plots`: exact pre-Stage-1 input semantics. All 54
  `.dat` files are byte-identical to
  `stage0_closure_20261105/diagnostics_off/production/PT/plots`.

The final canonical clean `-O3` executable SHA-256 is
`3661a7cdfb4b97a0ff011b23fc463e78d8ad7555302a49ba0b28db0b30871199`;
the canonical input SHA-256 is
`50cffc75488d26a0c3ccc1fb1289118b8729f756f08f2ade4e196a81cedc91c5`.
The entire `build/` directory was removed before each of the canonical,
alternate, cadence, legacy, and final canonical compilations.

### Stage 1 merge order and conflict rules

Do not merge yet. After the full comet model is tested and the user explicitly
authorizes master-integration work, create a disposable integration branch
from the then-current `master` and apply Stage 1 in this dependency order:

1. merge `CometReactionNetwork.pm::reachable_species_from()` and its existing
   reaction-network tests;
2. add `CometPopulationPolicy.pm` together with its Perl test;
3. merge `input/hartley2.pl`, regenerate from representative legacy and
   selector inputs, and inspect the generated species arrays;
4. merge the matching declarations and enum values in `RadialSimulation.h`;
5. merge runtime selection, cadence, reporting, lock preflight, and chemistry
   constraint changes in `RadialSimulation.cpp` as one reviewable application
   unit;
6. merge the Stage 1 runtime/legacy tests and runner;
7. merge `halley1d.input` and documentation last.

During conflict resolution preserve these invariants:

- group membership is database-derived and order-independent;
- precedence is explicit species, derived group, then global;
- the absent interval retains complete-sampling-break legacy behavior;
- adaptation runs only at a completed iteration, with no continuation
  particle in flight;
- regular and chemistry-driven adaptation cannot alter locked species;
- the face/chemistry plan is preflighted before any block mutation;
- emergency common-factor adaptation can override locks only for buffer
  safety and must report the event;
- the new parser/controller adds no random draw when no Stage 1 policy changes
  a legacy input.

Rerun the focused runner, a clean compile for the target master toolchain, the
four-rank weight validator, the alternate-source invariance gate, the cadence
run, and the 54-file legacy comparison on the integration branch. Then rerun
the existing DIM1 transport/chemistry suite and one unrelated DIM3 smoke test.
Present the exact integration diff and results to the user and stop before any
operation that changes `master`.

## 2026-11-05: Stage 2 component-wide carrier transaction implemented

Stage 2 of `CHEMISTRY.md` is implemented and its pure and integrated
transaction gates pass on `1d_comet_dev`. No commit, branch switch, rebase,
cherry-pick, integration branch, or operation on `master` was performed in
this round. The result replaces the Stage-1 two-phase controller that could
split a minor daughter locally and then undo that split when the face solver
raised the complete component to its largest pre-existing carrier.

### Source files added or modified

| File | Detailed change and reason | Later master-merge dependency |
| --- | --- | --- |
| `srcHartley2/ComponentCarrierPlanner.h` | **New.** Header-only, AMPS/MPI/geometry-independent graph planner. It discovers face-connected components, validates the positive dyadic carrier lattice, computes total-particle/total-cell means, selects one carrier per species/component, preserves locks, and calculates exact split capacity before returning a plan. Keeping this layer pure makes DIM2/DIM3 topology testable before those comet drivers use the controller. | Add before `RadialSimulation.cpp`. Preserve the integer power-of-two and empty-node contracts if the planner is moved into core later. |
| `srcHartley2/RadialSimulation.cpp` | Replaces independent per-block selection plus post-commit maximum propagation with one collective plan/preflight/commit transaction. It gathers authoritative owner counts/carriers, constructs dimension-independent face/subface adjacency, solves population and transitive chemistry constraints on component scalars, aggregates repeated products, audits every reachable reaction, preflights exact clone capacity with the configured reserve, commits only accepted plans, synchronizes once, verifies bitwise source locks, and rechecks committed carriers/reaction allocation. Empty blocks can be refined without creating particles. The superseded Stage-1 planning/face functions were removed. The integrated validator now measures all conserved moments, sparse/empty/two-ended behavior, capacity rollback, locks, reaction limits, and ghosts. | Merge only after Stage-1 generated policy arrays/locks. Treat the planner, chemistry fixed point, capacity preflight, commit, and validator as one logical unit; taking only the new commit loop or only the graph planner reopens the allocation defect. Reconcile any master changes in AMR face access, owner/ghost exchange, or particle-buffer reservation line by line. |
| `srcHartley2/RadialSimulation.h` | Updates the public numerical-control contract from block-local adaptation to component-wide planning and documents fixed-point chemistry, transactional capacity rejection, face multiplicity one, and the expanded validation mode. | Merge with `RadialSimulation.cpp`; stale Stage-1 comments would misdescribe runtime semantics. |
| `srcHartley2/tests/test_stage2_component_carrier_planner.cpp` | **New.** Pure C++ tests for a single block, exact lower/upper boundaries, long chains, empty and two-ended components, disconnected components, synthetic 2D/3D graphs, unequal dyadic inputs, source locks, invalid non-dyadic inputs, and exact capacity rejection. | No AMPS build dependency; run first during conflict resolution. Extend this test when higher-dimensional neighbor rules change. |
| `srcHartley2/tests/run_stage2_planner_tests.sh` | **New.** Builds the pure planner test in a temporary directory and leaves no repository binary. | Merge with its C++ test. |
| `srcHartley2/tests/test_stage2_integrated_logs.py` | **New.** Strictly parses exactly one legacy and one Stage-2 PASS record from each 1/4/12-rank log; rejects missing per-owner checks, absent MPI ghost checks, face/reaction limit violations, partial mutations, and decomposition-dependent chemistry results. | Diagnostic field names are a tested interface; update producer and parser together. |
| `srcHartley2/tests/run_stage2_integrated_tests.sh` | **New.** Routes externally archived 1/4/12 logs through the strict parser and explicitly skips when evidence paths are not supplied. Python bytecode output is disabled. | Keep validation archives external; never fabricate runtime evidence in Git. |
| `srcHartley2/CHEMISTRY.md` | Records the implemented algorithm, exact gates, archive/hash, and distinction between closed transaction tests and the still-required minor-daughter science pilot. | Reconcile after code and tests so later stages retain the true completion boundary. |
| `srcHartley2/README.md` | Documents the component transaction, conserved-moment/capacity tests, commands, and archived results. | Merge after runtime/test paths are final. |
| `srcHartley2/STATUS.md`, `handoff.md` | Add recovery state, validation evidence, and detailed master integration order/conflict invariants. | Documentation only; do not use it to authorize a master operation. |

The top-level `Makefile` was already dirty and remains user-owned. Root
`hartley2.pl`, `amps`, `build/`, root input links, meshes, `PT`, logs, and
`validation_runs` are generated or external evidence and must not be staged
as Stage-2 source.

### Validation evidence

The exact archive root is:

```text
validation_runs/stage2_component_carrier_20261105
```

The build directory was removed before each recompilation. A clean
`make -j8 EXTRACOMPILEROPTIONS=-O3` completed with only the known legacy DIM1
aggressive-loop warnings and sandbox OpenMPI interface notices. The frozen
executable/input/`picGlobal.dfn` SHA-256 values are:

```text
amps.stage2    c40c6bb44871209672cd6c7ecc6e16cc838ca5ea3e3e0184970696d273909939
halley1d.input 50cffc75488d26a0c3ccc1fb1289118b8729f756f08f2ade4e196a81cedc91c5
picGlobal.dfn  c3d2bca73f068cff357ae6d1f5002c97dc738b6a3190c712b2787890c92e2726
```

`run_stage2_planner_tests.sh` passes. The same executable then passes
`--validate-weight-adaptation` at 1, 4, and 12 MPI ranks. Each run reports one
owner check per rank for component population, empty blocks, the two-ended
chain, source locks, reaction allocation, split conservation, roulette
conservation, and capacity rollback. Ghost checks are 0/28/84 as expected;
maximum face multiplicity is exactly 1, maximum reaction allocation is 20 for
a limit of 20, and partial mutations are zero. The strict cross-rank parser
reports:

```text
PASS Stage-2 integrated component carriers 1rank-owners=1-ghosts=0 4rank-owners=4-ghosts=28 12rank-owners=12-ghosts=84
```

The existing Stage-0 suite also passes after the implementation. These gates
validate the transaction and decomposition behavior; they do not establish a
converged Rubin science result. A production minor-daughter pilot remains the
next acceptance step before Stage 3 or master integration.

### Stage 2 master-merge order and conflict rules

Do not merge yet. After the minor-daughter pilot and full comet validation are
complete, and only after explicit user authorization, create a disposable
integration branch from the then-current `master` and apply Stage 2 in this
order:

1. confirm that the full Stage-1 policy/parser/runtime commit range is already
   present and rerun the Stage-1 selector/legacy gates;
2. add `ComponentCarrierPlanner.h` and its pure C++ test/runner, then run that
   test before touching AMPS runtime code;
3. merge the authoritative owner gather and dimension-independent face graph
   construction in `RadialSimulation.cpp`;
4. merge component population selection and chemistry fixed-point planning as
   one unit, preserving repeated-product aggregation and locked-parent/product
   handling;
5. merge exact clone/reserve preflight and the plan-only rejection path before
   merging mutation code;
6. merge commit, one-time owner/ghost synchronization, exact lock verification,
   face/reaction postchecks, and removal of the old two-phase functions;
7. merge the integrated validator and log parser/runner;
8. update header/README/CHEMISTRY/STATUS documentation last, then rebuild from
   a deleted `build/` directory and rerun 1/4/12 MPI gates.

During conflict resolution preserve all of these invariants:

- every calculation leaf has exactly one authoritative owner contribution;
- component population is `sum(owner particles)/sum(component cells)`, not a
  mean of block means and not a populated-cell-only statistic;
- an empty block may change carrier but may not acquire a model particle;
- carrier changes are integral powers of two on the existing species lattice;
- direct-source regular weights remain bitwise unchanged;
- every reachable positive-rate unit-unimolecular reaction is audited, with
  repeated product terms summed before applying the hard cap;
- clone capacity is computed from the final population+chemistry plan and all
  ranks accept before the first mutation; stochastic future roulette may not
  be counted as available capacity;
- owners synchronize to boundary copies only after a successful commit, and
  both senses of every face remain at multiplicity at most one;
- emergency common-factor coarsening remains a separately logged buffer-safety
  override and must not be conflated with the regular component transaction;
- a rejected transaction leaves list links, represented moments, and carrier
  bits unchanged.

Do not merge generated build copies or validation output. On the disposable
integration branch, run the pure planner, Stage 0 and Stage 1 suites, the
1/4/12 Stage-2 validator, one unrelated DIM3 mesh smoke test, and the eventual
minor-daughter science pilot. Show the exact diff, source commit range,
excluded files, test logs, performance/memory impact, and science result to
the user, then stop before any operation that changes `master`.

## 2026-11-05: Stage 3 systematic residual implementation and validation

Stage 3 correction-free residual variance reduction is implemented on
`1d_comet_dev`. No commit, branch switch, rebase, cherry-pick, integration
branch, or master operation was performed. Three read-only sub-agent audits
were used to separate Stage-3 requirements from later individual-weight work,
identify the owner-block migration APIs, and define the validation matrix.

### Detailed source-file map

| File | Change and reason | Master-integration note |
| --- | --- | --- |
| `src/models/chemical_solver/SystematicResidualAccumulator.h` | **New.** Pure plan/commit accumulator and counter-based phase hash. `Plan()` returns multiplicity and next phase without mutation, allowing particle capacity to be accepted before state changes. | Add before comet runtime code. It has no AMPS mesh/MPI dependency and should remain reusable by other chemistry applications. |
| `srcHartley2/CometChemistry.h/.cpp` | Adds optional Bernoulli/systematic mode, separate phase seed, aligned owner-block state, stable concrete-product indexing, common initial phases for exact-equal carriers, lazy epoch initialization, regular-carrier rephasing, a 53-bit layout/ABI signature, diagnostic bounds, and transaction-safe integration in `ProcessParticle()`. Adds the migration/restart validator. | Merge as one logical transaction. Preserve all-zero commit and first-use/capacity-deferral no-mutation ordering. Do not replace concrete-product slots with one slot per species or term. Increment the layout ABI version whenever field meaning or ordering changes. |
| `srcHartley2/CometDiagnostics.h/.cpp` | Marks Bernoulli sigma unavailable in systematic mode and reports a conservative fixed-carrier-epoch error bound. Adds the selected residual mode and uncertainty method to JSON/CSV. | Keep CSV compatibility and JSON `null`/blank-CSV semantics; never reuse independent-Bernoulli variance for correlated systematic draws. |
| `srcHartley2/RadialSimulation.cpp` | Registers chemistry block data before block allocation; closes/rephases accumulator epochs before an accepted regular carrier transaction; explicitly preserves phases under ratio-preserving emergency common scaling; extends cell-measure audit with exact inner/outer radii; writes authoritative center-cell density/occupancy CSVs for complete diagnostic windows. | Registration order is core local W/dt fields, sampling layout, application residual data, then block allocation. Do not call the epoch hook for a rejected transaction. Do not substitute node-interpolated Tecplot rows for center-cell exposure. |
| `srcHartley2/CometModel.h/.cpp` | Adds the generic comet random seed and an explicit `RandomSeedConfigured` bit. | The flag is essential: an input that omits the seed must not receive a new `rnd_seed()` call. |
| `srcHartley2/RadialMain.cpp` | Adds `--validate-systematic-residual` and `--random-seed`; selects a physical seed+rank stream only when configured/overridden and retains a decomposition-independent keyed phase stream. | Preserve the no-call legacy path and the one-focused-validator rule. |
| `srcHartley2/main.cpp` | Replaces the historical literal seed 100 with the generic model value and registers systematic block storage in the legacy geometry backend. | Default remains 100, so the legacy 3-D source semantics do not change. Registration is a no-op in Bernoulli mode. |
| `input/hartley2.pl` | Adds validated, duplicate-rejecting `CometRandomSeed` and `ChemistryResidualSampling=BERNOULLI|SYSTEMATIC`. Systematic mode requires enabled chemistry and an explicit seed. | The residual-mode C++ definition must remain on one physical line because `ChangeValueOfVariable()` replaces complete definition lines. Merge parser and source definitions together. |
| `input/halley1d.input` | Selects seed 104729, systematic residuals, and diagnostics for the Stage-3 Halley case. | Case policy, not a generic default. Do not copy these records into Hartley input unless that case is intentionally converted. |
| `src/pic/pic_restart.cpp` | Extends previously unsupported species-local-weight restart: saves the global minimum-weight cache, layout length, and full application-owned data for every owner block including empty blocks; restores only the current owner, then publishes restored weights/timesteps/application bytes to ghosts before movement. Changes three restart/checksum tree recursions from hard-coded eight children to `2^DIM`. The global-weight file layout is unchanged. | This is the only existing PIC core source modified; the reusable accumulator header is a new core-model file. Split restart review into (1) dimension fix, (2) local-weight format, and (3) post-read ghost publication. Old local-weight restart aborted, so there is no supported old local file format to preserve. Still run DIM1 local and unrelated DIM3 global restart tests on master. |
| `srcHartley2/visualization/adaptive_radial_bins.py` | **New.** Non-destructive exposure-controlled radial aggregation of `stage3_cells.*.csv`, with exact support/volume, retained zero/underexposed bins, and independent-seed uncertainty on common supports. It explicitly rejects node-interpolated `pic.*.dat` as an exact exposure source. | Keep it read-only with respect to AMPS outputs. Replica meshes must match; never interpolate mismatched meshes silently. |
| `srcHartley2/tests/test_stage3_residual_accumulator.cpp`, `run_stage3_accumulator_tests.sh` | **New.** Analytic mean/variance, strict prefix error, plan/retry, common packet, raw round-trip, and invalid-input tests. | Run before AMPS build during conflict resolution. |
| `srcHartley2/tests/test_stage3_adaptive_bins.py`, `test_stage3_diagnostic_semantics.py` | **New.** Nine adaptive tests cover support selection, conservation, zeros, invalid exposure, exact geometry joins, provenance/mismatched replicas, distinct explicit seeds, sampling-length consistency, and complete CLI output. Four semantic tests prevent a Bernoulli sigma from being used in systematic mode while retaining readability of pre-Stage-3 Bernoulli reports. | Run with bytecode disabled in clean CI if repository cleanliness matters. |
| `srcHartley2/tests/run_stage3_parser_tests.sh` | **New.** Isolated generated-tree test for explicit and omitted defaults plus invalid/duplicate seed, invalid/duplicate mode, missing systematic seed, and chemistry-off rejection. | Restore source defaults before each fixture; a previously configured build must not supply omitted values. |
| `srcHartley2/tests/test_stage3_integrated_logs.py`, `run_stage3_integrated_tests.sh` | **New.** Strict protocol for exactly one PASS record at 1/4/12 ranks, including owner counts, packet/plan-retry/first-use/epoch checks, nonzero MPI migration, populated restart continuation/checksum, ghost publication, owner-only mutation, prefix error below one, and successful center/report output. | Producer field names are a tested interface; update validator and parser together. Logs remain validation artifacts, not source. Real capacity-rejection rollback remains covered by the inherited Stage-2 integrated transaction test. |
| `srcHartley2/CHEMISTRY.md`, `README.md`, `visualization/README.md`, `STATUS.md`, `handoff.md` | Documents algorithm, limitations, evidence, recovery, and merge sequence. | Merge after code/test paths stabilize. Documentation does not authorize a master operation. |

The top-level `Makefile` was already dirty and remains user-owned. Root
`hartley2.pl`, `halley1d.input`, `amps`, `build/`, meshes, `PT`, logs, and
`validation_runs/` are generated links/artifacts/evidence and must not be
staged as Stage-3 source.

### Focused validation evidence

The evidence root is:

```text
validation_runs/stage3_systematic_20261105
```

The final clean executable, canonical Stage-3 input, generated PIC definition,
exact commands, durations, and source/evidence fingerprints are recorded in
`RUN_MANIFEST.md`. The executable and input SHA-256 values are:

```text
amps.stage3    572abc22185a8ef6df16a22669491046e55fff6ff6e31fe8cfb7e50a3137f0e2
halley1d.input 0d89b536e6daa2dea7251bb57920fb6f1dd2546efd7c454d25b6dcccc26eccf5
```

The build directory was removed before every compilation. Pure accumulator
tests pass for means 0.25, 1.0, and 2.75; the systematic cumulative variance
matches `g(1-g)` and every fixed-phase prefix error is below one. The tests
also cover the double/`LONG_MAX` conversion boundary. All nine adaptive-bin
tests, four uncertainty-semantic tests, eight parser cases, and the real
2,984-cell command-line smoke test pass. Stage-0 counters at 1/4/12 ranks are
identical, and the archived 10,000-iteration report passes the observed-C2
selection/carrier gate. The complete archived Stage-1 runtime/cadence and
54-file legacy equivalence checks pass. Stage-2 pure planning and real 1/4/12
transaction rollback checks pass.

The integrated Stage-3 validator passes at 1, 4, and 12 MPI ranks. Every run
checks 373 owner blocks and 36,554 accumulator records, 15,666 exact-carrier
packet pairs, 373 non-mutating plan/retry cases, 373 first-use no-mutation
cases, 373 carrier-epoch resets, and 373 restart continuations. The forced
redistribution changes ownership of 0/9/33 nonzero block states at 1/4/12
ranks. Each restart includes 1/4/12 controlled particles, preserves the
particle checksum, and republishes restored state to 0/6/22 boundary ghosts.
Owner-only state checks pass and partial mutations are zero; maximum tested
prefix error is `4.99479983965357754e-01`. Every rank count also writes 2,984
exact center-cell records for C2 plus a complete uncertainty report. The
strict parser reports:

```text
PASS Stage-3 integrated systematic residual 1rank-owners=373-migrated=0 4rank-owners=373-migrated=9 12rank-owners=373-migrated=33
```

The focused Stage-3 plan/retry check is not described as a forced
`ProcessParticle()` capacity rejection. Actual particle/list/state rollback
under rejection is supplied by the passing Stage-2 integrated transaction
archive, which is a prerequisite of this implementation. Together these
tests close the intended implementation behavior. They do not close the
scientific decision gate: the required four post-burn-in Bernoulli and four
systematic pilots have not been run, because accepted burn-in duration,
comparison radius range, maximum CV, and maximum zero-cell fraction are
unspecified.
Stage 4 is neither approved nor rejected by these focused tests.

### Stage 3 master-merge sequence and conflict rules

Do not merge yet. After the comet science pilot is accepted and the user
explicitly authorizes integration, create a disposable branch from the
then-current `master` and proceed in this order:

1. rerun and verify the complete Stage-0/1/2 dependency stack on master;
2. merge the three `2^DIM` restart-recursion fixes independently and run DIM1,
   DIM2 (if available), and DIM3 restart/checksum smoke tests;
3. merge the species-local block-state restart format, then test empty and
   populated owner blocks across changed MPI decompositions and verify the
   post-read owner-to-ghost exchange;
4. add the pure systematic accumulator and run its no-AMPS test;
5. merge generic seed declarations plus parser controls, and explicitly test
   both an omitted legacy seed and an invalid/duplicate/systematic-without-seed
   input;
6. merge block-state registration/layout/phase initialization and verify one
   stable slot per concrete product, explicit alignment, and layout signature;
7. merge the `ProcessParticle()` plan/preflight/commit integration together
   with regular/emergency carrier epoch handling;
8. merge systematic diagnostic bounds and exact radial support output;
9. merge focused runtime validators, strict log parser, and adaptive-bin tool;
10. merge Halley case selection and documentation last, delete the exact build
    directory, rebuild, and rerun 1/4/12 tests with one frozen executable.

During conflict resolution preserve these invariants:

- omitted residual mode is Bernoulli, reserves no block bytes, consumes no
  additional random draw, and an omitted seed does not call `rnd_seed()`;
- the phase stream is separate from physical-event randomness and keyed by
  stable IDs rather than owner rank or first-use order;
- only owners initialize/read/update phase state; ghosts are synchronization
  copies, never stochastic actors;
- all product phases are planned before the capacity test and committed only
  with the complete physical reaction; all-zero committed events still
  advance their residual;
- exact-equal carrier siblings share initial packet phase but keep separate
  concrete-product records;
- a regular accepted carrier change starts a new bounded-error epoch, whereas
  emergency common scaling preserves phases only while all relevant ratios
  are exactly unchanged;
- systematic diagnostics report a bound and `sigma=null`, never a Bernoulli
  variance relabeled as systematic uncertainty;
- local restart contains empty owner blocks and application data; global
  restart layout remains unchanged, restored owner data is published to
  ghosts, and the application signature rejects same-size layout drift;
- tree recursion always uses `1<<DIM`, not a 3-D constant;
- original AMPS cell output remains unchanged and adaptive output never drops
  zero/underexposed samples.

On the disposable integration branch, also run one unrelated DIM3 global-
weight restart case and compare a legacy Bernoulli run byte-for-byte. Present
the exact diff, core/application split, executable hashes, test logs, and
science-gate report to the user, then stop before changing `master`.

The generic chemistry network executor remains intentionally radial-only.
Non-radial chemistry is rejected until master has a geometry-neutral
batch/capacity/emergency policy with explicit input values; this limitation
must not be weakened merely because generic 3-D geometry and sources exist.

## 2026-11-05: permanent Stage 3 scientific decision gate

The user-approved Stage-3 scientific decision is now encoded and
focused-tested on branch `1d_comet_dev` at starting HEAD `d40a6a39ccda`. No
merge, commit, or operation on `master` was performed. The long eight-case
science matrix has not been launched; this section records implementation
readiness, not a Stage-4 go/no-go result.

### Files added or modified in this gate round

| File | Change and reason | Later master-merge note |
| --- | --- | --- |
| `srcHartley2/RadialMain.cpp` | Adds `--chemistry-residual-sampling bernoulli|systematic`. The override is applied before `RadialSimulation::Initialize()` requests application block storage, allowing one frozen executable to run both arms without reserving residual bytes in Bernoulli mode. Systematic override still requires an explicit input or CLI seed. | Merge after the Stage-3 parser/model enums. Do not move the override below `PIC::Init_BeforeParser`/radial initialization or the block ABI can be selected for the wrong mode. |
| `srcHartley2/CometDiagnostics.h`, `CometDiagnostics.cpp` | Adds structured executed-emergency counters by buffer, sampling-population, and allocation-burst trigger; records seed/mode/MPI provenance; retains collective request counts; and emits these values in every complete JSON window. Existing systematic fixed-epoch bounds remain distinct from Bernoulli sigma. CSV rows were re-audited to exactly 50 columns. | Merge producer and strict Stage-0/Stage-3 consumers together. Preserve rank-zero event ownership before MPI reduction and do not reinterpret a request as an executed emergency. |
| `srcHartley2/RadialSimulation.cpp` | Records sampling emergency requests when first raised and executed trigger categories when emergency roulette actually commits. The Stage-3 center-cell CSV now carries block/local indices, block cell count, resolved per-cell target, lock, derived class, and policy source. | Merge only after Stage-1 population-policy arrays/classes exist. Keep writer order synchronized with `adaptive_radial_bins.py`; policy fields are science evidence, not inferred by the evaluator. |
| `srcHartley2/visualization/adaptive_radial_bins.py` | Reads and validates optional authoritative block-policy metadata, preserves legacy tables that contain none of those fields, rejects partial/corrupt metadata, and carries policy identity through replicate merging. | Backward compatibility is intentional for exploratory archives. The science evaluator separately requires the new complete metadata. |
| `srcHartley2/validation/stage3_science_policy.json` | **New.** Freezes four seeds, two modes, 12 MPI ranks, 10k burn-in, three accepted 10k windows through 40k, 10--100,000 km altitude, five important daughters, exposure/CV/zero/population/mode-mean/closure/emergency/Rubin gates. | Treat as reviewed science policy, not a generic comet default. Any threshold change needs explicit scientific review and matching code/tests. |
| `srcHartley2/validation/stage3_gate_common.py` | **New.** Supplies strict policy validation, canonical hashes, atomic JSON, exact eight-run order, and canonical-input checks. Approved numerical boundaries are pinned in code so replacing the JSON/hash cannot silently relax them. | Merge before runner/evaluator. Do not replace exact policy checks with permissive defaults. |
| `srcHartley2/validation/run_stage3_science_gate.py` | **New.** Durable operational runner with `plan/run/resume/status`, frozen executable/input/`picGlobal.dfn`/policy hashes, alternating sequential eight-case order, atomic attempt state, advisory locking, signal forwarding, fresh-at-zero retry, exact artifact validation, and per-artifact checksum verification. | Keep independent from the evaluator. It must never claim a scientific pass and must never reuse a partial `PT` tree as restart evidence. |
| `srcHartley2/validation/evaluate_stage3_science_gate.py` | **New.** Read-only archive evaluator. It independently re-hashes frozen configuration and every protected run artifact, then requires exact provenance and complete windows, all-eight-run minimum-exposure support, authoritative population policy, convergence/closure/emergency gates, Systematic statistical quality, Bernoulli/Systematic mean agreement, and explicit Rubin-baseline acceptance. It writes JSON/CSV/Markdown and stable outcome exit codes. | Merge after diagnostic/table producers. Preserve exact-cell support and do not interpolate missing/zero samples away. Baseline promotion remains a human decision. |
| `srcHartley2/validation/README.md` | **New.** Defines the experiment, thresholds, build/preflight, run/resume/reboot recovery, baseline review, evaluation outcomes, and focused test commands. | Merge with validation tools so operators do not launch an unproven ad-hoc matrix. |
| `srcHartley2/tests/test_stage3_science_gate.py`, `run_stage3_science_gate_tests.sh` | **New.** Fourteen synthetic full-archive cases cover exact inclusive thresholds, common support, valid/missing/mismatched Rubin baselines, statistical failures, emergencies, nonconvergence, malformed tables, bad policy metadata/manifests, and post-run artifact mutation. | Keep outcome and exit-code expectations synchronized with the evaluator. These tests intentionally do not launch AMPS. |
| `srcHartley2/tests/test_stage3_adaptive_bins.py` | Adds optional policy-metadata preservation and corrupt block-index rejection; the suite now has eleven cases. | Merge with the metadata reader/writer changes. |
| `srcHartley2/CHEMISTRY.md`, `README.md`, `visualization/README.md`, `STATUS.md`, `handoff.md` | Records the approved decision gate, its true completion boundary, recovery steps, and later merge dependencies. | Documentation does not authorize a master merge or count as science evidence. |

No new core source was added or modified in this gate round. The earlier
Stage-3 core restart changes in `src/pic/pic_restart.cpp` remain a separate
dependency and retain their previously documented DIM-aware traversal,
application-block restart, particle-record, and ghost-publication tests.

### Focused verification

Before compiling, exactly `build` was removed.
`make -j8 EXTRACOMPILEROPTIONS=-O3` passed. The generated artifacts are:

```text
amps                  3d5c0c3111b57a9c2e4b2d939d2252f19faa43e3ab89c08bbd6f0b2e6c6c8e70
input/halley1d.input  0d89b536e6daa2dea7251bb57920fb6f1dd2546efd7c454d25b6dcccc26eccf5
build/pic/picGlobal.dfn c3d2bca73f068cff357ae6d1f5002c97dc738b6a3190c712b2787890c92e2726
```

The same executable passed 10,000-event chemistry validation with both CLI
residual modes. The Systematic focused runtime passed at 1, 4, and 12 MPI
ranks, including 373 owner blocks, 36,554 accumulators, migration, restart,
0/6/22 restored ghost records, and the new center-table schema. The current
one-rank table has 2,984 exact rows and the adaptive-bin smoke test accepts
its authoritative 50--200 `CHEMICAL_ONLY` metadata. Additional results:

```text
Stage-0 suite and archived production report       PASS
Stage-1 parser/runtime/legacy-equivalence suite    PASS
Stage-2 pure planner and archived 1/4/12 parser    PASS
Stage-3 pure accumulator/parser                    PASS
Stage-3 archived 1/4/12 integrated parser          PASS
Stage-3 adaptive bins                              11/11 PASS
Stage-3 diagnostic semantics                       4/4 PASS
Stage-3 synthetic science gate                     14/14 PASS
Stage-3 operational plan preflight                 PASS
```

The operational preflight selected `/usr/bin/orterun`, 12 ranks, OMP=1, and
the approved alternating matrix. It did not create or launch a production
science archive.

### Merge into master only after scientific validation

Do not merge yet. After the complete eight-run archive has been evaluated,
the Rubin candidate has been explicitly reviewed/accepted, and the user has
authorized integration, prepare a disposable branch from the then-current
`master` and merge in this order:

1. merge the already validated Stage-0 diagnostic foundation and Stage-1
   population-policy tables/classes;
2. merge the Stage-2 component planner and transactional capacity behavior;
3. merge the Stage-3 accumulator/storage/restart implementation and all pure
   plus 1/4/12 focused tests, including the narrowly scoped core restart
   commit;
4. merge the runtime residual-mode CLI before the science runner so one binary
   retains the early storage-selection invariant;
5. merge structured emergency counters and the authoritative center-cell
   policy schema as one producer-interface change;
6. merge adaptive-bin metadata validation, then common policy, science policy,
   runner, evaluator, and synthetic tests in dependency order;
7. add the explicitly accepted Rubin baseline and compact evaluation products
   only after scientific review; keep raw `PT`, binaries, build trees, and
   machine-local logs outside Git;
8. merge application documentation/input last, delete the exact integration
   `build/` directory, rebuild, and repeat every focused gate plus one complete
   science evaluation against the retained archive.

During conflict resolution preserve the exact policy values, fresh-at-zero
attempt rule, all-eight-run minimum exposure, fully contained exact cells,
unmodified zeros, distinction between requests and executed emergencies,
Systematic-bound versus Bernoulli-sigma semantics, source/chemistry closure,
and explicit baseline acceptance. Present the commit range, integration diff,
excluded generated files, artifact hashes, 1/4/12 logs, eight-run wall/memory
cost, final outcome, and Rubin comparison to the user, then stop for explicit
approval before changing `master`.

## 2026-11-05: Stage 4 individual-correction compatibility complete

Stage 4 is implemented and focused-tested on branch `1d_comet_dev` at current
HEAD `c0dc90e`. No commit, merge, or operation on `master` was performed. The
user explicitly authorized this stage after the Stage-3 implementation; this
record does not claim that the still-pending long Stage-3 scientific matrix
made the go/no-go decision.

The supported configuration is the DIM1 radial backend with species-local
timesteps, species-local block weights, local-timestep face resampling,
individual particle corrections, the general NTC collision model, MPI-only
focused diagnostics, and the configured internal degrees of freedom.
Simplified NTC and MF are explicitly rejected. CUDA continuation, hybrid NTC
diagnostic reduction, field-line/periodic movers, and correction-on 2-D/3-D
chemistry are not qualified.

### Core files modified

| File | Detailed modification and reason | Master-merge dependency/invariant |
| --- | --- | --- |
| `src/pic/pic.h` | Extends the opt-in local-timestep flux diagnostic with separate internal-energy input/output/variance fields. Adds a disabled-by-default general-NTC validation record and reset/snapshot API for trials, accepted pairs, self-pair rejection, particle updates, invalid/clamped probabilities, and extrema. | Merge each declaration with its implementation. These counters are process-local validation observability, not production science samples; do not enable them by default or use them in hybrid mode without a thread reduction. |
| `src/pic/pic_mover.cpp` | Removes the correction-on compile-time prohibition only for the existing local-timestep face-resampling path. Validates positive finite `c`; keeps `c` fixed in survivors/clones; retains the base `W_old*dt_new/(W_new*dt_old)` multiplicity because `c` cancels; changes all diagnostic carriers to `W*c/dt`; and audits rotational/vibrational internal-energy flux. | Merge with the `pic.h` flux fields. Preserve full-record cloning and fixed-`c` policy. Do not revive the legacy correction-transform route, broaden to CUDA, or omit correction from diagnostics. |
| `src/pic/pic_particle_collision_models.cpp` | Makes general NTC carrier-correct: caches `c` and prefix sums in reusable species buckets, selects trials in proportion to `W*c`, uses independently drawn same-species pairs and discards self pairs, fixes the limiter to `ancoll/(N0*N1)`, validates every update probability against `[0,1]`, and exposes opt-in diagnostics. It also prevents the one-particle same-species loop. Simplified NTC now terminates explicitly in correction-on builds with a nonzero AMPS exit code. | Merge with the `pic.h` NTC API. Preserve weighted selection, independent-draw normalization, self-pair discard, the corrected parentheses, finite checks, and explicit simplified-kernel rejection. Re-run both correction-off and correction-on collision tests; do not fold MF/specialized kernels into this change. |

No new core source file was added. The three core modifications above are the
complete Stage-4 core delta. Existing Stage-2/3 core work remains separately
documented and must not be conflated with this review unit.

### Comet application and test files modified or added

| File | Detailed modification and reason | Master-merge dependency/invariant |
| --- | --- | --- |
| `srcHartley2/LocalTimeStepValidation.h`, `.cpp` | Adds exact correction histograms for `1/16,1/4,1,4,16`, an all-`c=1` identity mode, deterministic rotational/vibrational energies, corrected post-move flux checksums, and compact coverage/audit records for local, MPI, and multiple-face first/second-order mover paths. Correction-off builds report the unit pattern through the same protocol. | Merge after the core mover diagnostic. Preserve exact powers-of-two corrections, separate correction/long-sum tolerances, and rank-sensitive MPI coverage. |
| `srcHartley2/RadialMain.cpp` | Adds `--validation-unit-corrections` and a narrow `--validate-unsupported-simplified-ntc` expected-failure mode. | Merge with the validators that consume these options. The unsupported mode is a guard test, not a runtime collision selector. |
| `srcHartley2/RadialSimulation.cpp` | Allows corrections only with local face resampling and general NTC; keeps MF/other models rejected. Corrects chemistry capacity accounting to one daughter per concrete product in correction-on builds. Makes regular/emergency adaptation audits use `W*c` while preserving `c`. Adds focused sampling, general-NTC relaxation/probability, boundary deletion/re-emission, restart, migration, ghost publication, and rollback checks with compact records. | Merge after Stage-2 adaptive weights, transactional capacity hooks, restart support, and the Stage-4 core primitives. Preserve the radial-only guard, concrete stoichiometric transaction count, callback restoration after focused NTC, and no-partial-mutation requirements. |
| `srcHartley2/CometChemistry.cpp` | Reaches the exact corrected-daughter path in correction-on builds: one daughter per concrete product with `c_daughter=(W_parent*c_parent/dt_parent)/(W_daughter/dt_daughter)`. Extends the ensemble through the five-value correction range and validates parent/daughter physical carriers. | Merge after the Stage-4 radial setup and before any Stage-5 windows. This is carrier compatibility only; do not label Stage 5 complete or remove the Systematic/correction-on incompatibility. |
| `srcHartley2/tests/run_stage4_integrated_tests.sh` | **New.** Runs wide and unit correction mover, adaptation/core, chemistry, and expected simplified-NTC failure at exactly 1, 4, and 12 ranks without configuring or compiling AMPS or overwriting an archive. | Keep executable/input provenance outside the runner. Output roots must remain new and correction-on configuration must be verified before launch. |
| `srcHartley2/tests/test_stage4_integrated_logs.py` | **New.** Strictly parses compact evidence instead of trusting exit status. Requires wide corrections, exact unit identity, flux and probability bounds, NTC work/relaxation, adaptation/rollback/ghost paths, sampling/boundary/restart/migration, and chemistry. Optional correction-off logs enforce identical 1/4/12 path coverage and independent physical-test success for all-`c=1`. | Merge with the C++ record schema. A field rename is an interface change and requires synchronized updates; do not make missing evidence optional. |
| `srcHartley2/CHEMISTRY.md`, `README.md`; `src/pic/README.md`; `STATUS.md`; `handoff.md` | Records the algorithm, supported matrix, test evidence, recovery, exclusions, and later merge sequence. | Documentation is not authorization to merge or to activate corrections in canonical inputs. |

The canonical `input/halley1d.input` was restored exactly: correction mode is
off, the collision model is the default general NTC, and chemistry residuals
are Systematic. The root `halley1d.input` is a symlink to that tracked file.
The pre-existing dirty root `Makefile`, generated executable/parser links,
historical PT/log directories, and validation archives are not Stage-4 source
changes and must not be staged accidentally.

### Validation evidence

Every compile in this round first removed exactly
`build`. The full correction-on matrix passed with
200,000 mover particles per case and 200,000 chemistry events per case:

```text
validation_runs/stage4-correction-on-full-20261105
PASS Stage-4 all-c=1 correction identity 1rank=pass 4rank=pass 12rank=pass
PASS Stage-4 correction-on integrated matrix 1rank=pass 4rank=pass 12rank=pass
```

Frozen correction-on provenance is:

```text
amps.correction-on
  9e9f79398f5bf335f3ac0a283504c8782773d7a72f4bd6c51ce5fa2152fb5732
halley1d.correction-on.input
  81fe26cb4aaafff2af65ed33d3461da41da51e940678c7dea99b941f9bed84fc
picGlobal.correction-on.dfn
  eb968bf7d4b95f43564070c8f64fecf772a989aafa63da8477980483a726c971
```

The correction-off identity archive is
`validation_runs/stage4-correction-off-reference-20261105`. At all three MPI
sizes it exercised the exact same local/MPI/multiple-face, clone, correction-
preservation, and internal-energy coverage as the correction-on all-`c=1`
arm. Both independently passed the physical-moment tests. Maximum flux
checksum errors were about `5.47e-12` on and `1.39e-11` off, below `1e-9`.

`validation_runs/stage4-correction-off-regression-20261105` retains passing
1/4/12-rank Stage-2 adaptation and Stage-3 Systematic restart/migration logs.
The correction-on MF setup failed before mesh setup as required in
`validation_runs/stage4-unsupported-mf-20261105b`; its input and generated
macro file are retained with the rejection log. The full Stage-0 runner,
Stage-1 parser/reaction checks, Stage-2 planner, Stage-3 accumulator/parser,
14-case synthetic science gate, 11 adaptive-bin tests, and 4 diagnostic-
semantics tests also passed. No AMPS/MPI validation process is active.

The final clean executable is the canonical correction-off build:

```text
amps                   773d9606b7282b0105da59d3d996df135223522e6162b723284903ff5b73b6d5
input/halley1d.input   0d89b536e6daa2dea7251bb57920fb6f1dd2546efd7c454d25b6dcccc26eccf5
build/pic/picGlobal.dfn c3d2bca73f068cff357ae6d1f5002c97dc738b6a3190c712b2787890c92e2726
```

### Later merge to master: detailed sequence

Do not merge yet. After comet validation is accepted and the user explicitly
authorizes master integration:

1. create a disposable integration branch from the then-current `master` and
   record both base SHAs; do not switch the working validation branch with
   uncommitted artifacts present;
2. exclude the dirty root `Makefile`, root symlinks/generated parser,
   executable, `build`, meshes, `PT*`, logs, and raw validation archives from
   every source commit;
3. first integrate prerequisite Stage-2 adaptive block carriers/capacity hooks
   and Stage-3 restart/application-storage work, with their original tests;
4. integrate the mover unit: the relevant `pic.h` internal-flux fields plus
   `pic_mover.cpp`, then build correction-off and run the local-timestep tests;
5. integrate the collision unit: the relevant `pic.h` NTC diagnostic plus
   `pic_particle_collision_models.cpp`, then test correction-off and the
   correction-on simplified-kernel rejection;
6. integrate `LocalTimeStepValidation.*` and the CLI options, followed by the
   radial setup/adaptation/core validator and `CometChemistry.cpp` exact
   daughter path;
7. add the Stage-4 shell/Python tests and documentation, clean-regenerate a
   correction-on executable, run the complete 1/4/12 matrix, then clean-
   regenerate the canonical correction-off executable and repeat Stage-0--3
   regressions plus the combined identity audit;
8. review the integration diff for unintended changes to input defaults,
   collision selection, random streams in correction-off builds, particle
   record layout, or unsupported geometry; retain correction-off/Systematic
   in Halley and keep Hartley input separate;
9. present commit range, exact diff, excluded files, hashes, MPI logs,
   memory/wall cost, Stage-3 science outcome, and Rubin comparison to the
   user, then stop before changing `master`.

During conflicts, preserve: fixed correction at faces; `W*c/dt` diagnostics;
full-record cloning; weighted general-NTC selection; independent same-species
draws with self-pair discard; `ancoll/(N0*N1)`; finite probability checks;
explicit simplified/MF rejection; one corrected daughter per concrete
product; transactional capacity rollback; ghost publication after adaptation;
and default-off behavior. Stage-5 correction windows and importance sampling
must remain a later, separately approved change.

## 2026-11-05: unified Stage 0--4 test orchestrator

`srcHartley2/tests/run_all_stage_tests.sh` is now the authoritative single
entry point for the complete previously executed focused/archive test stack.
No commit, branch switch, or master operation was performed. The orchestrator
does not compile AMPS or mutate `build`; this separates testing from the
correction-off/correction-on build transition and prevents a failed test from
leaving the canonical workspace in the wrong configuration.

### Files changed

| File | Change and reason | Later master-merge note |
| --- | --- | --- |
| `srcHartley2/tests/run_all_stage_tests.sh` | **New.** Strictly preflights compilers/interpreters/MPI, the correction-off Stage-3 parser fixture, and every required external archive. Runs all Stage-0/1 runners, both Stage-2 runners, all Stage-3 runners plus the two standalone Python suites and real CLI smoke test, final correction-off Stage-2/3 regressions, Stage-4 provenance/rejection checks, and the combined correction-on/all-`c=1` matrix. Optional fresh mode delegates to the existing Stage-4 MPI runner using a caller-supplied executable/new archive. | Merge after the individual runners and their documented archive layouts. Preserve failure on missing evidence and the read-only build policy. Do not replace dated explicit paths with newest-directory discovery. |
| `srcHartley2/tests/run_stage0_tests.sh` | Adds `HARTLEY2_MPIEXEC` as an optional launcher override, resolved before compilation. Direct use still defaults to `mpiexec`. This lets the unified runner use the same explicit launcher for Stage 0 and fresh Stage 4. | Backward compatible. Keep the standalone 1/4/12 ranks and `--oversubscribe`; do not reinterpret the override as permission to change rank coverage. |
| `srcHartley2/tests/run_stage3_parser_tests.sh` | Replaces five `rg` checks with portable `grep -E/-F` checks so the all-stage suite does not depend on a developer-specific ripgrep installation. Test patterns and semantics are unchanged. | Preserve `--` before patterns and keep fixed-string matching for generated C++ definitions. |
| `srcHartley2/README.md` | Documents archived/fresh operation, build separation, evidence requirements, coverage, and commands. | Merge with the runner; documentation does not embed validation artifacts in Git. |
| `srcHartley2/STATUS.md`, `handoff.md` | Records evidence, recovery, and integration behavior. | Documentation only; it does not authorize a master merge. |

### Orchestration contract

Archived mode is the default and requires all previously retained evidence.
The child runners' standalone optional `SKIP` behavior is intentionally not
allowed: the wrapper supplies every Stage-0 production/baseline path, every
Stage-1 runtime/cadence/legacy path, and every Stage-2/3 1/4/12 log. A missing
file or directory fails preflight before any test starts.

The current `build` is read only and is used only by the Stage-3 parser test.
The wrapper requires the assembled block/species inputs, generated sources,
and `_INDIVIDUAL_PARTICLE_WEIGHT_OFF_`; it never regenerates them. Stage-4
archived mode verifies the three frozen provenance hashes before parsing
compact records. It independently checks the MF and simplified-NTC rejection
messages and applies the all-`c=1` correction-on/off identity comparison.

Fresh Stage-4 mode requires `--stage4-executable` and a nonexistent
`--stage4-output-root`. It delegates the 1/4/12 MPI runs to the existing
Stage-4 runner and then repeats the combined identity test against the frozen
correction-off archive. It does not copy or install the executable. A caller
that builds it must first remove exactly `build`,
and must separately restore a clean canonical correction-off build after any
configuration work.

### Validation result

`--plan` passed strict tool/artifact/configuration preflight. The first full
execution correctly caught an invalid hard-coded Stage-3 sampling length in
the new CLI smoke wiring. The runner was changed to read and validate the
authoritative `sampling_length` column rather than assuming the production
value. The second complete archived execution passed:

```text
Hartley2 all-stage test summary status=PASS stage4-mode=archived
  tests: run=17 passed=17 skipped=0 failed=0
  test-groups: run=5 passed=5 skipped=0 failed=0
PASS Hartley2 all-stage suite tests=17 test-groups=5 stage4-mode=archived
```

The complete suite was repeated with `PATH=/usr/bin:/bin`, which excludes the
workspace's `rg`; the Stage-3 parser and all 17 groups passed. Therefore `rg`
is not an undeclared runtime dependency.

The result collector was also tested with two deliberate compiler failures.
It continued through all remaining independent jobs, returned nonzero, and
reported `run=17 passed=15 skipped=0 failed=2` plus five group results
(`passed=3 failed=2`). Both `stage2.planner` and `stage3.accumulator` appeared
in `failed-tests` with their exit codes and distinct absolute log paths; both
logs and `summary.txt` existed and were nonempty. The injected archive was
removed after verification and is not validation evidence.

Default successful runs remove temporary logs. Failed/skipped runs preserve
their automatically allocated `/tmp/hartley2-all-stages.*/logs` directory.
`--log-root NEW_PATH` retains all per-test logs and `summary.txt` even on
success and refuses an existing path. A child `SKIP` is counted separately in
both test and stage-group summaries but still makes this complete suite exit
nonzero.

That result includes the Stage-0 standalone MPI signature at 1/4/12 ranks,
three Stage-0 Python checks and 20 production reports; complete Stage-1
runtime/cadence and 54-file legacy equivalence; Stage-2 pure and integrated
checks; Stage-3 accumulator, 11 adaptive-bin, 4 diagnostic-semantics, 8 parser
case, 14 science-gate, 2,984-cell CLI, and integrated checks; correction-off
Stage-2/3 regression replays; Stage-4 provenance and both unsupported-model
rejections; and the correction-on matrix/identity check. Temporary unit
binaries and CLI products were written under `mktemp` directories and removed
by traps. No Python bytecode directory remains.

For later master integration, add this orchestrator only after all child
runners and record schemas are stable. On the integration branch, assemble a
complete evidence bundle or pass an explicit relocated `--evidence-root`,
run `--plan`, then archived mode. Build a separate correction-on executable
with the exact clean-build procedure and run fresh mode into a new archive.
Finally restore/rebuild correction-off and rerun archived mode. Present both
terminal PASS records and hashes to the user before any operation on master.

## 2026-11-05: Stage 5 exact daughters and correction windows complete

Stage 5 is implemented and focused-tested on `1d_comet_dev`, based on commit
`5da7fe664f986ac74191dcffb010322cb6dcda53`. No commit, branch switch, merge,
or operation on `master` was performed. Three sub-agents independently audited
the existing chemistry/runtime algorithm, designed the input policy, and
designed the pure/integrated validation coverage. Their recommendations were
reconciled here; they did not merge or commit code.

The supported runtime remains deliberately narrow: DIM1 radial mesh,
species-local block weights and timesteps, local-timestep face resampling,
individual corrections, general NTC collisions, and Bernoulli chemistry
residual allocation. Correction windows with Systematic residuals, simplified
NTC, MF, CUDA, hybrid batched chemistry, or DIM2/DIM3 chemistry are not
qualified. The planner is geometry-neutral, but the particle-list executor is
not yet qualified outside the radial application.

### Algorithm and conservation contract

The Stage-4 exact chemistry path remains authoritative. Every selected analog
reaction creates one particle for every concrete stoichiometric product and
sets

```text
c_daughter = (W_parent*c_parent/dt_parent)/(W_daughter/dt_daughter).
```

Stage 5 does not alter reaction probability, channel selection, recoil, or
parent destruction. It applies species-resolved correction windows only after
a complete `PIC::TimeStep()`. The inclusive interval is inert. For `c<lower`,
roulette survival probability is `c/target` and a survivor receives `target`.
For `c>upper`, the fixed representation is `ceil(c/target)` particles, bounded
by the input hard maximum. The first children receive `c/n`; the final child
receives the arithmetic remainder. The planner requires the production-order
sum to equal the input correction exactly.

The runtime uses two owner-only traversals rather than retaining a plan vector
for every live particle. Preflight validates finite positive corrections,
zero mover remainder, species policy, exact clone count, per-particle split
limit, and reserve-aware local capacity. Roulette deletions are not credited
toward capacity. `MPI_Allreduce(MPI_MIN)` makes acceptance collective: one
rank's rejection occurs before all RNG draws, list changes, tracker changes,
and window diagnostics on every rank. Commit snapshots each original next
link, so head-inserted children are not processed recursively. Full-record
cloning preserves position, velocity, internal modes, species, chemistry
stage, and application bytes; only the child correction and tracker identity
are changed. Base-weight adaptation, when due on the same boundary, commits
first. No ghost block weight publication is required because `c` travels in
the full particle record.

### Files added or modified

| File | Detailed change and reason | Later master-merge dependency/invariant |
| --- | --- | --- |
| `input/CometCorrectionWindowPolicy.pm` | **New.** Geometry- and network-independent input resolver for global, `DIRECT_SOURCE`, `CHEMICAL_ONLY`, and explicit-species windows. Validates finite ordered triples, fixed-split compatibility, duplicate/unknown selectors, direct-source `c=1`, required global/limit records, and strict default-off behavior. | Merge before `hartley2.pl`. Preserve precedence `species > group > global` and reject inactive records instead of silently ignoring them. Do not combine this dimensionless-correction policy with block population ranges. |
| `input/test_comet_correction_window_policy.pl` | **New.** Parser/resolver regression covering precedence and invalid/off-mode configurations. | Merge with the module; keep it independent of a generated AMPS tree. |
| `input/hartley2.pl` | Parses the Stage-5 switch, window records, and split limit; derives direct/reaction-reachable groups after all databases load; generates one lower/target/upper/source row per compiled species; and requires radial DIM1 plus core correction mode when enabled. | Merge after the policy module and before C++ consumers. Keep generated arrays on one physical source line for `ChangeValueOfArray`; regenerate, never hand-merge, `build/main/RadialSimulation.cpp`. |
| `srcHartley2/CorrectionWindowPlanner.h` | **New.** Pure AMPS/MPI/mesh-independent keep/roulette/split planner and aggregate capacity preflight. Implements inclusive bounds, fixed ceil splitting, exact final remainder, hard split limit, overflow checks, and no-mutation batch rejection. | Merge before the runtime executor. This header is the future DIM2/DIM3 reuse boundary; do not add radial mesh or RNG calls to it. |
| `srcHartley2/RadialSimulation.h` | Declares generated policy arrays, split limit, switch, and the focused initialized-mesh validator. Detailed comments distinguish the correction window from the chemistry transaction bound. | Merge in lockstep with generated definitions and parser source codes. |
| `srcHartley2/RadialSimulation.cpp` | Adds policy reporting/validation, the two-pass collective owner-list executor, reserve-aware capacity transaction, full-state split/tracker commit, roulette deletion, safe post-step scheduling in every run branch, and the production-path validator. The validator covers keep boundaries, split moments/state, global zero-capacity rollback, restart, analytic roulette statistics, direct `c=1`, and empty populations. | Merge after Stage-2/3 capacity/restart and all Stage-4 mover/collision/chemistry prerequisites. Preserve owner-only traversal, zero-remainder gate, no roulette capacity credit, all-rank preflight, base-adaptation-before-window ordering, and full-record clone semantics. |
| `srcHartley2/RadialMain.cpp` | Adds mutually exclusive `--validate-correction-windows` CLI routing and usage text. | Merge with the validator; it is a focused test mode, not a production policy selector. |
| `srcHartley2/tests/test_stage5_correction_window_planner.cpp` | **New.** Pure C++ cases for inclusive policy, expected roulette, exact production-order split sum, hard-limit/capacity rollback, multiple species, repeated products, and empty batches. | Merge with the planner. Exact equality is intentional, not a tolerance typo. |
| `srcHartley2/tests/run_stage5_correction_window_tests.sh` | **New.** Compiles/runs the pure planner and executes the Perl policy suite in a temporary directory. | Merge after both pure tests; it must not configure AMPS or leave a binary. |
| `srcHartley2/tests/run_stage5_integrated_tests.sh` | **New.** Runs correction-window and exact chemistry validators at exactly 1, 4, and 12 ranks with Bernoulli residuals, immutable new output roots, separate commands/stdout/stderr, and a strict final checker. | Merge after CLI/runtime and parser. Preserve explicit caller-supplied executable and refusal to overwrite evidence. |
| `srcHartley2/tests/test_stage5_integrated_logs.py` | **New.** Requires one passing window record per rank, exact owner/rank and path flags, consistent policy, analytic seven-sigma roulette bounds, and Stage-4 exact chemistry coverage at every decomposition. | Merge with the C++ output schema. A renamed field is an interface change; zero exit status alone is not acceptance. |
| `srcHartley2/tests/run_all_stage_tests.sh` | Extends the general runner through Stage 5, adds archived/fresh Stage-5 modes, strict evidence preflight, the pure job and matrix job, a sixth group, and Stage-5 mode in the final summary. | Merge after the Stage-5 runner/checker. Preserve Stage-3's requirement that the current build be correction-off; archived Stage-5 evidence comes from a separate frozen binary. |
| `srcHartley2/CHEMISTRY.md`, `srcHartley2/README.md` | Records the implemented algorithm, grammar, validation values/evidence, unsupported matrix, science boundary, and operator commands. | Merge after code/test schemas stabilize. Do not treat the validation window values as a production Rubin policy. |
| `srcHartley2/STATUS.md`, `handoff.md` | Records complete file rationale, recovery state, and master-integration sequence. | Documentation is not permission to merge. Reconcile append-only history rather than replacing it with an older whole-file version. |

No core file under `src/pic`, `src/meshAMR`, or `src/models` was added or
modified for Stage 5. Stage 5 consumes the already documented Stage-4 core
individual-correction mover/collision interfaces and the Stage-2/3
application capacity/restart machinery. Keep this distinction when reviewing
a future integration diff.

The canonical `input/halley1d.input` was only staged temporarily to build the
validator. The accepted build used global `0.25|1|4`, direct-source
`0.5|1|2`, chemical-only `0.0625|0.25|4`, and a maximum split count of 32.
After the archive was produced, the tracked input was restored byte-for-byte
to its original correction-off/Systematic form (SHA-256
`0d89b536e6daa2dea7251bb57920fb6f1dd2546efd7c454d25b6dcccc26eccf5`).
Those focused values have not been scientifically calibrated.

### Validation evidence

Before each compile, exactly `build` was removed.
The correction-on build completed without compiler/linker errors. The pure
suite printed:

```text
Comet correction-window policy test passed
PASS Stage-5 correction-window planner policy=1 roulette=1 split=1 capacity=1 multi-species=1 repeated-products=1 empty=1
```

The full runtime archive is:

```text
validation_runs/stage5-correction-window-full-20261105
PASS Stage-5 correction-window integrated matrix 1rank=pass 4rank=pass 12rank=pass
```

At each rank count, the initialized-mesh controller passed keep, split,
roulette, capacity rollback, restart, direct-source, and empty-list checks.
It used 65,536 roulette candidates per owner block. The paired chemistry run
used 200,000 events per generated reaction and retained the complete exact
daughter number/mass/momentum/energy checks from Stage 4. The validation
binary hash was
`6f640dd3db0fad77f96b392b19e20fe9dcde65e698207db24a1dc9014a58d89f`;
its generated `build/pic/picGlobal.dfn` hash was
`eb968bf7d4b95f43564070c8f64fecf772a989aafa63da8477980483a726c971`.
The first sandboxed launch was preserved separately as
`stage5-correction-window-full-20261105.sandbox-launch-failed`; Open MPI had
no visible network interface and AMPS never started. The unrestricted rerun
is the acceptance archive.

This evidence closes implementation Stage 5, not the long scientific pilot.
No production convergence, independent-seed C2 uncertainty, clone
decorrelation time, parent-profile invariance, minor-species flux closure, or
new Rubin comparison is claimed here.

### Detailed later merge sequence

Do not merge yet. Tell the user before any switch, fetch/rebase, cherry-pick,
integration branch, commit on `master`, or merge involving `master`. After the
comet science validation is accepted and explicit authorization is given:

1. record the current `1d_comet_dev` tip, then-current `master` tip, merge
   base, dirty/untracked exclusions, and all Stage 0--5 archive hashes;
2. use a disposable integration branch created from that recorded master;
   never reuse this validation worktree with raw archives or generated files;
3. integrate prerequisite Stage 0--3 diagnostics, policies, component
   planner, transactional capacity, accumulator, and restart changes in their
   documented order and pass their correction-off suites;
4. integrate the Stage-4 core mover and general-NTC changes as separate,
   reviewable core commits, including correction-off identity and unsupported
   collision-model gates; do not mix them with Stage-5 parsing;
5. integrate the Stage-4 exact corrected-daughter application path and prove
   the 1/4/12 all-`c=1` identity before adding any window;
6. add `CometCorrectionWindowPolicy.pm` and its Perl test, then modify
   `input/hartley2.pl`; regenerate and inspect species arrays from at least one
   legacy off input and the focused on input;
7. add the pure planner/test, then declarations, runtime executor, scheduling,
   CLI, and focused validator. Resolve conflicts by preserving conservation
   and rollback invariants, not by accepting whichever branch compiles;
8. add the Stage-5 shell/Python runners and general-runner wiring. Build the
   correction-on validation binary only after deleting the exact build tree,
   rerun the 1/4/12 matrix into a new archive, and parse it strictly;
9. restore canonical correction-off `halley1d.input`, delete the exact build
   tree again, rebuild, and run all archived Stage 0--5 groups. Verify omission
   of window records is behavior-identical to the old input;
10. review for unintended changes to random streams in off mode, generated
    policy order, source classification, collision selection, chemistry event
    count, tracker state, capacity reserve, and unsupported geometry;
11. exclude `build/`, `amps`, root generated parser/input copies, meshes,
    `PT*`, raw validation archives, sandbox failure evidence, and the
    pre-existing dirty `Makefile` from source commits;
12. present the exact integration diff, commit range, exclusions, build/input
    hashes, 1/4/12 records, memory/wall comparison, science-pilot statistics,
    convergence/closure result, and Rubin comparison to the user, then stop
    before changing `master`.

During conflict resolution, preserve: correction mode default off; selector
precedence; direct-source `c=1`; one exact daughter per concrete product;
inclusive windows; fixed ceil splitting; exact correction sum; full-record
cloning; unique child tracker setup; no roulette capacity credit; collective
preflight/rollback; owner-only traversal; base adaptation before windows; and
the separation between the chemistry event-product limit and window split
limit. Stage 6 packet splitting and Stage 7 channel importance sampling remain
future, separately approved development.

### Final restored workspace and complete regression

After the correction-on acceptance archive, the canonical input was restored,
exactly `build` was removed again, and AMPS was
cleanly rebuilt from `input/halley1d.input`. The generated configuration has
`_INDIVIDUAL_PARTICLE_WEIGHT_OFF_`, general NTC, window mode false, and split
limit `-1`. Final hashes are:

```text
amps                         162707be2849ff3da60318389219c2c54bbdb4b5f6c9987f375f8cb25a4cf40d
input/halley1d.input         0d89b536e6daa2dea7251bb57920fb6f1dd2546efd7c454d25b6dcccc26eccf5
build/pic/picGlobal.dfn      c3d2bca73f068cff357ae6d1f5002c97dc738b6a3190c712b2787890c92e2726
build/main/RadialSimulation.cpp d54a0ce8bdf1d115324b3e2d126b0b0fe351d5657d5370f43c6b1c423eb70cae
```

`run_all_stage_tests.sh --plan` passed against that correction-off build. The
complete archived suite then passed all Stage 0--5 work:

```text
Hartley2 all-stage test summary status=PASS stage4-mode=archived stage5-mode=archived
  tests: run=19 passed=19 skipped=0 failed=0
  test-groups: run=6 passed=6 skipped=0 failed=0
PASS Hartley2 all-stage suite tests=19 test-groups=6 stage4-mode=archived stage5-mode=archived
```

There is no active AMPS/MPI validation process. The pre-existing tracked
`Makefile` difference and historical/generated untracked files are unrelated
to Stage 5 and were preserved. `git diff --check` passes. No commit was made.

## Stage 6 daughter-packet implementation and validation (2026-11-05)

This section supersedes the preceding Stage-5-era statement that Stage 6 was
future work. Stage 6 is now implemented and focused-tested on branch
`1d_comet_dev`, starting from commit
`9eeb06041893cbbd5fe8a92bef7a2ec5d285b61b`. No commit, branch switch, merge,
or operation involving `master` was performed. The user must be told and must
explicitly authorize the operation before any attempt to integrate this work
with `master`.

### Approved behavior and implementation invariants

The input policy selects a fixed daughter-packet count for each generated
reaction. An explicit reaction ID has highest precedence. Without an explicit
reaction rule, each concrete product resolves `species > CHEMICAL_ONLY >
global`, and the reaction uses the maximum resolved count across all of its
products. Repeated products participate normally. A product-free sink is
always `N=1`; explicitly targeting one is rejected even if the requested count
is one. Selectors that match only zero-rate reactions cannot create an
effective split, and enabled packet mode must contain at least one positive-rate
reaction with `N>1`.

For a reaction with `K` concrete stoichiometric products and selected packet
count `N`, the production transaction creates exactly `K*N` particles. Each
packet obtains an independent recoil solution containing all `K` products.
For each product, the original Stage-4 correction is divided across its `N`
representatives; the last representative receives the arithmetic remainder so
the corrections reconstruct the original value in production order. The
parent slot is reused once, the remaining records are allocated, all committed
daughters are marked processed for the current chemistry stage, and every new
particle receives a valid independent tracker identity.

Preflight proves `K*N` against integer overflow, the configured packet maximum,
the radial hard product limit, and available particle capacity before mutation.
Capacity rejection leaves the parent byte-identical and consumes no
commit-dependent random numbers. Normal chemistry, emergency population
forecasts, and block-weight adaptation all use the same `K*N` cardinality.
Stage-5 correction windows run normally after packet creation: packet daughters
have no exemption from keep, roulette, or split decisions. The packet and
window policies therefore must be qualified together in later science runs.

With all packet records omitted, the generated reaction counts are one and the
runtime retains the historical one-recoil random-number path. This strict-off
identity is a merge requirement. The implementation is currently supported
only by the qualified DIM1, individual-correction, general-NTC, Bernoulli
chemistry path. The policy resolver and pure planner are geometry independent,
but this work does not claim DIM2/DIM3 corrected-chemistry support.

### Complete Stage-6 source inventory

| File | Detailed Stage-6 change and reason | Later master-merge dependency/invariant |
| --- | --- | --- |
| `input/CometDaughterPacketPolicy.pm` | **New.** Parses mode, declared maximum, global/group/species/reaction selectors; resolves the approved precedence and maximum-across-products rule; rejects duplicates, unknown targets, sinks, zero-rate-only effective splits, invalid counts, inactive-mode records, and `K*N` transactions above the hard limit; emits the generated C++ policy header. | Merge before `input/hartley2.pl`. Preserve stable reaction IDs and strict default-off behavior. Do not resolve from generated array position. |
| `input/test_comet_daughter_packet_policy.pl` | **New.** Pure Perl regression for off identity, all precedence levels, mixed and repeated products, sinks, zero-rate cases, malformed/duplicate records, maximum enforcement, transaction limits, and generated-header integrity. | Merge with the policy module and keep independent of an AMPS build tree. |
| `input/CometReactionNetwork.pm` | Adds `packet_policy_reactions()`, a defensive description containing stable ID, nominal rate, and concrete product stoichiometry after database expansion. | Merge before policy resolution. Preserve existing reaction order and keep private network state immutable. |
| `input/hartley2.pl` | Instantiates the packet policy, parses the three Stage-6 record types, retains the hard product limit for cross-validation, enforces chemistry/DIM1/individual-correction/Bernoulli/Stage-5 prerequisites, resolves after the full network and chemical-only set exist, and writes the generated header. | Merge after both input modules. Regenerate outputs; do not hand-merge `build/main` copies. Omission must leave all packet counts at one. |
| `srcHartley2/CometDaughterPacketPolicy.generated.h` | **New.** Safe source-tree fallback with zero generated reactions and count-one sentinel arrays. The application generator replaces it only in the generated build tree. | Merge before compiling `CometChemistry.cpp`. Keep `<cstddef>` and the zero-reaction sentinel so direct source-tree compilation is well formed. |
| `srcHartley2/DaughterPacketPlanner.h` | **New.** AMPS/MPI/mesh-independent fixed-count correction partition and capacity planner for `N=1,2,4,8` and general valid counts. | Merge before runtime chemistry. Preserve exact final-remainder reconstruction and mutation-free rejection. |
| `srcHartley2/CometChemistry.h` | Declares packet switch/maximum, generated-policy accessors, focused validator, and comments defining the transaction contract. | Merge in lockstep with `CometChemistry.cpp` and generated header. |
| `srcHartley2/CometChemistry.cpp` | Loads and validates the generated policy; retains the off/`N=1` RNG path; implements independent recoil packets, exact correction partitions, `K*N` capacity transaction, parent reuse, tracker/processed-state commit, and policy diagnostics. Adds deterministic ensemble and real `C2_2C -> C+C` production-transaction validation. Validator-only centered extended-precision `(count,mean,M2)` moments are merged with the parallel identity so rank regrouping is numerically stable. | Merge after Stage-4 exact corrected daughters and Stage-5 windows. Preserve per-packet recoil, full transaction rollback, DIM-only coordinate comparisons, high-precision diagnostics, and centered-moment merge. |
| `srcHartley2/RadialSimulation.cpp` | Changes all chemistry capacity mirrors, emergency forecasts, and adaptive block-weight projections from `K` to the runtime `K*N` worst-case transaction size. | Merge with chemistry runtime. Missing any one mirror can admit a partial allocation burst even if the chemistry executor itself is correct. |
| `srcHartley2/RadialMain.cpp` | Adds mutually exclusive `--validate-daughter-packets` routing and help text, using the shared chemistry validation event count. | Merge with the validator; this is a test mode, not a production input selector. |
| `srcHartley2/tests/test_stage6_daughter_packet_planner.cpp` | **New.** Pure C++ tests for identity, exact sums, `N=1,2,4,8`, repeated products, overflow/hard-limit/capacity rejection, and sink cardinality. | Merge with the pure planner. Exact production-order equality is intentional. |
| `srcHartley2/tests/run_stage6_daughter_packet_tests.sh` | **New.** Runs the Perl policy suite and compiles/runs the pure C++ planner in a temporary directory. | Merge after both pure tests; it must not configure AMPS or retain a test binary. |
| `srcHartley2/tests/run_stage6_integrated_tests.sh` | **New.** Runs Stage-5 windows, exact chemistry, and the packet validator at exactly 1, 4, and 12 MPI ranks into a caller-supplied immutable root. | Merge after runtime/CLI. Preserve refusal to overwrite evidence and separate command/stdout/stderr records. |
| `srcHartley2/tests/test_stage6_daughter_packet_logs.py` | **New.** Strictly parses all 46 reactions, 23 sinks, `N=1,2,4,8`, policy sources, represented moments, variance scaling, conservation, MPI identity, production transaction, links, processed flags, tracker state, and rollback. | Merge with the C++ output schema; field or precision changes are validation-interface changes. |
| `srcHartley2/tests/run_all_stage_tests.sh` | Extends the general runner to seven groups and 22 checks; adds archived/fresh Stage-6 modes, immutable evidence preflight, packet pure/matrix jobs, Stage-5 dependency replay, and SHA-256 provenance verification. | Merge after Stage-6 runners. Keep the current build correction-off; archived Stage-6 evidence uses its separate frozen correction-on binary. |
| `srcHartley2/CHEMISTRY.md` | Records policy grammar, precedence, transaction math, Stage-5 interaction, validation boundary, and remaining science work. | Merge after behavior is stable; do not reinterpret validation counts as Halley science parameters. |
| `srcHartley2/README.md` | Adds operator commands for pure, focused, archived, fresh, and all-stage Stage-6 validation. | Keep paths and CLI synchronized with the runners. |
| `srcHartley2/STATUS.md`, `handoff.md` | Records implementation rationale, hashes, evidence, recovery procedure, exclusions, and ordered master integration. | Append/reconcile these histories; documentation never grants permission to merge. |

No core file under `src/pic`, `src/meshAMR`, or `src/models` was added or
modified by Stage 6. The implementation consumes the previously developed
Stage-4 individual-correction mover/collision interfaces and the Stage-2/3
application capacity/restart machinery. `RadialSimulation.cpp` is an
application source despite its infrastructure role. The pre-existing tracked
root `Makefile` difference is unrelated and must not be attributed to or
staged with Stage 6.

### Focused and complete validation evidence

Before every Stage-6 compilation, exactly
`build` was removed. The accepted packet-enabled
build used packet maximum 8; global 1; `CHEMICAL_ONLY` 2; species `C2` 8;
reaction `H2O_OH_H` 4; reaction `C2_2C` 8; a 64-product transaction limit;
and deliberately broad Stage-5 windows. These are validation settings, not a
scientifically selected Halley policy.

The pure test result is:

```text
Comet daughter-packet policy test passed
PASS Stage-6 daughter-packet planner identity=1 exact-sum=1 packet-counts=1,2,4,8 repeated-products=1 capacity=1 sink=1
```

The immutable accepted runtime archive is:

```text
validation_runs/stage6-daughter-packet-full-20261105-r4
PASS Stage-5 correction-window integrated matrix 1rank=pass 4rank=pass 12rank=pass
PASS Stage-6 daughter-packet integrated matrix 1rank=pass 4rank=pass 12rank=pass reactions=46 sinks=23 packet-cases=4 events=200000
```

Its `provenance/SHA256SUMS` verifies:

```text
amps.stage6-correction-on                 c70e51c7a5e36f8626d3966126aab15d6e82b7bd22081ee93ac7b64c18558170
halley1d.stage6-validation.input          6b9a4479e1f57eecb1c03d59a594518912e381fa1b7caf84f9368b20c243abc2
picGlobal.stage6-validation.dfn           eb968bf7d4b95f43564070c8f64fecf772a989aafa63da8477980483a726c971
CometDaughterPacketPolicy.generated.h     274d5170bbaae9b74ea41a81a81034ba8c79ce03e7fb014b848eb454f7264b87
```

The sibling archives without `-r4` and with suffixes `-r2` and `-r3` are
development diagnostics, not acceptance evidence. They exposed, in order,
an invalid DIM1 comparison of inactive coordinates, insufficient numeric
print precision for strict rank comparison, and cancellation in raw
sum/sum-square variance reconstruction. The accepted implementation compares
only active dimensions, prints deterministic statistics at full precision,
and uses centered extended-precision parallel moments. Retain the failed
archives only for debugging history; never combine them with `-r4`.

After acceptance, canonical `input/halley1d.input` was restored byte-for-byte,
the exact build directory was removed again, and a clean correction-off build
completed. Its final state is:

```text
amps                                      4ed4147d6c640bc473a4a87caf89808c0438d3a3cbfe2d99dd8853f6f943218e
input/halley1d.input                      0d89b536e6daa2dea7251bb57920fb6f1dd2546efd7c454d25b6dcccc26eccf5
build/pic/picGlobal.dfn                   c3d2bca73f068cff357ae6d1f5002c97dc738b6a3190c712b2787890c92e2726
build/main/RadialSimulation.cpp           b120edd6fe20c0965a7bc71491b75439cb787d731d9d86a7e2d409a6ff6fc2b0
srcHartley2/CometDaughterPacketPolicy.generated.h ffa6c0387ec0b3b5909eecf71ee3d0167642a168c12c889f7e4fcd08fc4a0fa6
```

The canonical generated configuration is individual-correction off, general
NTC, Stage-5 windows off, Stage-6 packet mode false, and every generated packet
count one. The final all-stage preflight and complete archived replay passed:

```text
Hartley2 all-stage test summary status=PASS stage4-mode=archived stage5-mode=archived stage6-mode=archived
  tests: run=22 passed=22 skipped=0 failed=0
  test-groups: run=7 passed=7 skipped=0 failed=0
PASS Hartley2 all-stage suite tests=22 test-groups=7 stage4-mode=archived stage5-mode=archived stage6-mode=archived
```

This closes Stage 6 as an implementation and focused-compatibility stage. It
does not close the scientific validation gate: production packet counts still
require independent-seed variance/cost selection, transported convergence,
primary-profile invariance, minor-species flux closure, buffer-headroom
measurement, and a new Rubin comparison. Stage 7 channel importance sampling
also remains separate future work.

### Detailed Stage-6 additions to the later master-integration procedure

Do not merge yet. After the comet model and the scientific packet pilot are
accepted and the user explicitly authorizes master integration, follow the
existing Stage 0--5 sequence above, then add these Stage-6-specific steps:

1. record the exact `1d_comet_dev` tip, target `master` tip, merge base, dirty
   tracked files, untracked/generated exclusions, and SHA-256 manifests for
   every accepted Stage 0--6 archive;
2. create a disposable integration branch from that recorded master tip in a
   clean worktree; do not switch this validation worktree to `master` and do
   not copy raw `PT`, mesh, build, executable, or archive contents into source
   commits;
3. finish and validate Stages 0--5 in their documented order before applying
   any packet patch; Stage 6 depends on Stage-4 exact corrected products and
   Stage-5 post-cycle correction windows;
4. merge `CometReactionNetwork.pm`, then add
   `CometDaughterPacketPolicy.pm` and its Perl test. Run the pure policy test
   before changing the application parser;
5. merge the `hartley2.pl` records and generated-header step. Configure one
   legacy input with all Stage-6 records omitted and one focused on-input;
   inspect stable reaction IDs, policy sources, counts, sinks, zero-rate rows,
   and `K*N` bounds rather than trusting a successful configure alone;
6. add the fallback generated header and pure planner/test. Prove `N=1`
   identity and exact production-order correction sums before runtime wiring;
7. merge `CometChemistry.h/.cpp`, then the three `RadialSimulation.cpp`
   capacity mirrors and `RadialMain.cpp` validator. During conflicts preserve
   independent recoil per packet, parent byte-identity on rejection, tracker
   initialization, processed flags, active-dimension checks, and centered
   parallel moments;
8. add the Stage-6 shell/Python runners and general-runner provenance job.
   Freeze a newly clean-built packet-enabled executable/input/generated header
   and run 1/4/12 ranks into a new immutable archive; never overwrite or
   relabel `stage6-daughter-packet-full-20261105-r4`;
9. restore the canonical no-packet input, remove exactly the build directory,
   rebuild, and prove the correction-off/off-mode identity plus the complete
   22-test, seven-group archived suite;
10. review random-stream identity when `N=1`, selector precedence, sink and
    zero-rate rejection, transaction cardinality in every allocation forecast,
    Stage-5 scheduling order, conservation, rank invariance, provenance, and
    unsupported geometry/model gates;
11. exclude the pre-existing `Makefile` change and all generated root inputs,
    `build/`, `amps`, meshes, `PT*`, logs, raw validation archives, and failed
    diagnostic archives from source commits unless the user separately asks
    to version a curated evidence manifest;
12. present the exact integration diff and commit range, exclusions, all
    hashes, 1/4/12 records, science-pilot uncertainty/cost/closure results, and
    Rubin comparison to the user, then stop before any operation on `master`.

No AMPS or MPI process is active at this handoff. No commit or merge was made.

## Stage 7 stochastic channel importance sampling completed (2026-11-05)

Stage 7 is implemented and accepted as a focused implementation/compatibility
stage on branch `1d_comet_dev`. The pre-Stage-7 branch tip remains
`83d24f9222c941deeb8e009a9e0ff89d4e56ce46`; no commit, merge, branch switch,
or operation on `master` was performed. The unrelated pre-existing tracked
`Makefile` edit remains outside this work.

The implementation biases only the channel selected *after* an ordinary
physical reaction event has occurred. For one configured unit-unimolecular
parent family, the input must supply a strictly positive target probability
`q*` for every positive-rate sibling, including product-free sinks, and the
family must sum to one. The total physical reaction rate, parent lifetime, and
parent deletion are unchanged. A selected productive channel receives
likelihood `L=q/q*`; Stage 4 applies `L` to the unsplit product correction and
Stage 6 then creates `N` packets carrying `L/N` each. Sinks produce no model
particle but remain in selection and diagnostics. Unconfigured families use
the former raw-rate subtraction selector and preserve its RNG path.

### Stage-7 file inventory and rationale

New source and test files:

- `input/CometChannelImportancePolicy.pm` parses stable, case-sensitive
  reaction selectors, resolves complete families, rejects partial/duplicate/
  unknown/zero-rate/non-unit policies, checks normalization, and emits the
  generated C++ arrays.
- `input/test_comet_channel_importance_policy.pl` covers default-off identity,
  two- and many-channel families, sinks, tiny/dominant probabilities, record
  order, malformed policies, and generated sentinels.
- `srcHartley2/ChannelImportancePlanner.h` is the geometry-neutral pure
  selector/likelihood planner used to test analog identity, boundary choices,
  weighted moments, variance, and ESS without initializing AMPS.
- `srcHartley2/CometChannelImportancePolicy.generated.h` is a safe source-tree
  fallback with no configured reaction. A configured build generates the
  69-row application-specific copy in `build/main/`; do not hand-edit either
  generated table.
- `srcHartley2/tests/test_stage7_channel_importance_planner.cpp` validates the
  pure planner for analog, two-/many-channel, dominant/tiny, disabled, boundary,
  moment, variance, ESS, and error cases.
- `srcHartley2/tests/run_stage7_channel_importance_tests.sh` runs the Perl and
  pure C++ tests without modifying a configured build.
- `srcHartley2/tests/test_stage7_channel_importance_logs.py` strictly checks
  1/4/12-rank runtime logs and CSV/JSON sidecars, including reaction coverage,
  family coverage, likelihood identities, uncertainty, ESS, weighted rates,
  provenance, and rank-decomposition equivalence. Its floating tolerance is
  limited to accumulated binary64 fields; identifiers, counts, metadata,
  `q`, `q*`, and `L` remain strict.
- `srcHartley2/tests/run_stage7_integrated_tests.sh` runs an immutable frozen
  executable at 1, 4, and 12 MPI ranks and replays the Stage-5 window,
  Stage-6 packet, and Stage-7 channel validators into a new output root.

Modified source and test files:

- `input/CometReactionNetwork.pm` now exposes a defensive ordered reaction
  view for policy construction without changing the physical database.
- `input/hartley2.pl` accepts `ChemistryChannelImportanceMode` and
  `ChemistryChannelBiasedProbability(REACTION:<ID>)`, validates prerequisites,
  and generates the policy header. The root `hartley2.pl` was refreshed only
  as the normal mechanical configure copy and must not be staged as source.
- `srcHartley2/CometChemistry.h/.cpp` cache and validate complete configured
  families, preserve the exact analog path for unconfigured families, sample
  configured `q*`, apply `L` exactly once before packet division, retain
  unweighted parent destruction, record sinks, and provide the deterministic
  initialized validator. Capacity forecasts use the combined Stage-6
  transaction cardinality, not an after-allocation repair.
- `srcHartley2/CometDiagnostics.h/.cpp` add a separate Stage-7 CSV/JSON schema
  with per-reaction selection, likelihood, likelihood-square, and weighted-rate
  sums plus parent-family ESS. Keeping a separate schema avoids changing the
  Stage-0 production diagnostic contract.
- `srcHartley2/RadialMain.cpp` adds
  `--validate-channel-importance --chemistry-validation-events N`.
- `srcHartley2/tests/run_all_stage_tests.sh` adds Stage-7 pure, provenance,
  and matrix jobs; archived/fresh Stage-7 selection; sidecar preflight; and
  the new 25-test/eight-group summary.
- `srcHartley2/CHEMISTRY.md`, `README.md`, and
  `GENERIC_COMET_DESIGN.md` document the algorithm, input, prerequisites,
  diagnostics, tests, generic-comet boundary, limitations, and merge order.
- This file and root `handoff.md` record validation and later integration
  instructions.

No AMPS core source under `src/`, `srcInterface/`, `srcGeneral/`, `meshAMR/`,
or the shared model libraries was added or modified for Stage 7. The feature
uses the application-side Stage-4 correction carrier, Stage-5 windows, and
Stage-6 transaction machinery already developed on this branch. The clean
build copied application sources into `build/`; those copies are disposable
artifacts, not core edits.

### Accepted Stage-7 evidence

The validation-only input used complete equal-probability C2 and H2O families:

```text
ChemistryChannelImportanceMode=on
C2_2C=0.5                 C2_ION=0.5
H2O_OH_H=0.25             H2O_O_H2=0.25
H2O_O_2H=0.25             H2O_ION=0.25
```

These are test controls, not a recommended Halley policy. The accepted,
immutable archive is:

```text
validation_runs/stage7-channel-importance-full-20261105-r2
```

Its end-to-end runner passed:

```text
PASS Stage-5 correction-window integrated matrix 1rank=pass 4rank=pass 12rank=pass
PASS Stage-6 daughter-packet integrated matrix 1rank=pass 4rank=pass 12rank=pass reactions=46 sinks=23 packet-cases=4 events=200000
PASS Stage-7 channel-importance integrated matrix 1rank=pass 4rank=pass 12rank=pass families=27 reactions=69 sidecar-rows=69 events=200000
```

The archive's `provenance/SHA256SUMS` verifies:

```text
62d3a1e53eb2fa873c2a86f07a5000cbcb5813817ebf8267ba07986bd9ef92a8  amps.stage7-channel-importance
41ca2daa708a5cfdeaf2e7293579ecaa1d91b8fb9b479a05dd16c0fdeada4891  halley1d.stage7-validation.input
eb968bf7d4b95f43564070c8f64fecf772a989aafa63da8477980483a726c971  picGlobal.stage7-validation.dfn
e4f958f113d41f7d71b63ed7e8b093543ac978b065355e97a8f1974943fd1cd3  CometChannelImportancePolicy.generated.h
274d5170bbaae9b74ea41a81a81034ba8c79ce03e7fb014b848eb454f7264b87  CometDaughterPacketPolicy.generated.h
```

The unsuffixed Stage-7 archive is development evidence only: every AMPS case
passed, but the first post-run parser used an unrealistically tight repeated-
addition tolerance. The accepted `-r2` parser applies a documented
`8*N*epsilon` bound only to accumulated floating fields. Never relabel,
combine, or substitute the unsuffixed archive for `-r2`.

After acceptance, `.amps.conf` and canonical `input/halley1d.input` were
restored. The exact `build` directory was removed
before a clean canonical compile. The resulting build is individual-particle-
correction off, Stage-6 packet counts all one, and all 69 Stage-7 policy rows
unconfigured. Current hashes are:

```text
a99f39b913f2d59109c87e724b1179549dbd8115ad13761c55a76a439ea94ed4  amps
0d89b536e6daa2dea7251bb57920fb6f1dd2546efd7c454d25b6dcccc26eccf5  input/halley1d.input
c3d2bca73f068cff357ae6d1f5002c97dc738b6a3190c712b2787890c92e2726  build/pic/picGlobal.dfn
b120edd6fe20c0965a7bc71491b75439cb787d731d9d86a7e2d409a6ff6fc2b0  build/main/RadialSimulation.cpp
7f4e0be8f6a6f09e0b39a28642200839a1d32a7d195953befd043609632a9fbf  build/main/CometChannelImportancePolicy.generated.h
c1330629ecb9e0048d716085dafba93ea098959f1bc9cb842ea8616d933cf5ae  srcHartley2/CometChannelImportancePolicy.generated.h
```

The final complete archived replay passed:

```text
Hartley2 all-stage test summary status=PASS stage4-mode=archived stage5-mode=archived stage6-mode=archived stage7-mode=archived
  tests: run=25 passed=25 skipped=0 failed=0
  test-groups: run=8 passed=8 skipped=0 failed=0
PASS Hartley2 all-stage suite tests=25 test-groups=8 stage4-mode=archived stage5-mode=archived stage6-mode=archived stage7-mode=archived
```

This closes Stage 7 as an implementation/focused-compatibility stage only.
It does not choose scientifically optimal `q*`, prove transported steady-state
variance reduction, or close the Rubin comparison. Deterministic all-channel
branching and total-reaction forcing are intentionally absent. The runtime
executor remains validated only for DIM1; the policy/parser/planner are
geometry-neutral but a future DIM2/DIM3 executor must satisfy the same
correction, capacity, conservation, and diagnostics contracts.

### Detailed Stage-7 additions to the later master-integration procedure

Do not merge yet. After Stage-9 science validation and only after explicit
user authorization, extend the prior Stage 0--6 merge sequence as follows:

1. Record the exact development tip, target `master` tip, merge base, dirty
   tracked files, untracked/generated exclusions, and SHA-256 manifests for
   all accepted archives. Work on a disposable clean integration branch; do
   not switch this validation worktree to `master`.
2. Integrate and validate Stages 0--6 first. Stage 7 depends on Stage-4
   per-particle correction, Stage-5 post-transaction windows, and Stage-6
   packet accounting; conflict resolution must not bypass those contracts.
3. Merge the defensive reaction descriptor in `CometReactionNetwork.pm`, then
   `CometChannelImportancePolicy.pm` and its Perl test. Prove strict complete-
   family validation before changing runtime code.
4. Merge `hartley2.pl` generation and the fallback header. Configure one
   omitted/off input and one complete-family on-input; inspect all stable IDs,
   sinks, zero-rate exclusions, normalization, and all-zero off arrays. The
   root parser copy is mechanical and build/root generated artifacts are not
   source files.
5. Add `ChannelImportancePlanner.h` and its pure test. Require exact analog
   selector identity, analytic weighted means/variance/ESS, boundary behavior,
   and all rejection cases before runtime integration.
6. Merge `CometChemistry.h/.cpp`. Preserve the raw-rate/RNG branch for
   unconfigured families, physical total event probability, unweighted parent
   deletion, sink selection, one application of `L`, Stage-6 `L/N`, allocation
   forecasts, and rollback invariants.
7. Merge the separate diagnostics schema and CLI validator, followed by the
   strict parser/runners and general-runner jobs. Do not fold Stage-7 fields
   into the stable Stage-0 production schema during conflict resolution.
8. Remove the exact build directory, make a clean Stage-7-on build, freeze its
   executable/input/generated headers, and run Stage-5/6/7 validators at 1, 4,
   and 12 ranks into a new immutable archive. Never overwrite the accepted
   `stage7-channel-importance-full-20261105-r2` evidence.
9. Restore canonical `input/halley1d.input`, remove the exact build directory
   again, clean-build the correction/packet/importance-off application, and
   rerun all 25 tests in eight groups. Review off-mode RNG identity, parent
   lifetime, sinks, yields, conservation, ESS, capacity, and MPI invariance.
10. Run independent-seed transported Stage-9 pilots and a new Rubin comparison
    before proposing production `q*`. Do not promote the equal C2/H2O test
    policy to a science input.
11. Exclude the unrelated `Makefile`, root generated inputs/parsers, `build/`,
    `amps`, meshes, `PT*`, logs, raw validation archives, and failed diagnostic
    archives unless the user separately requests a curated manifest.
12. Present the exact integration diff and commit range, exclusions, hashes,
    test matrix, performance/memory effects, science uncertainty/closure, and
    Rubin products to the user, then stop before any `master` operation.

No AMPS/MPI validation process is active after the final replay. No commit or
merge was made.

## Stage 8 deterministic total-reaction forcing (2026-11-05)

Stage 8 is implemented and focused-tested on `1d_comet_dev`. The work started
from commit `c2da2d752b0f09bb39abd75b1b4a1c7432adfc50`. No commit, branch switch,
merge, or operation on `master` was performed. The tracked `Makefile` edit
predates Stage 8 and remains user-owned; do not include it in a Stage-8 commit.

### Accepted algorithm and invariants

The approved estimator is deterministic physical two-branch splitting, not a
biased Bernoulli occurrence probability. For total rate `lambda`, logical
local step `dt`, and incoming physical carrier `C`, production evaluates
`P=-expm1(-lambda*dt)` and retains survivor/reacted carriers `C*(1-P)` and
`C*P`. The reacted value and arithmetic survivor remainder close exactly in
binary64. No `P*` exists in the input or generated policy.

The reacted branch stays a parent particle until the mover reaches a reaction
time sampled from the physical exponential law conditional on an event in the
segment. A configured Stage-7 channel contributes `L=q/q*`; a Stage-6 packet
count contributes `/N`. Each concrete product packet therefore represents
`C*P*L/N`. Products in one packet reconstruct a common conservative recoil
solution from a stored counter-based key. A sink remains a pending parent
marker until the same event point and is then deleted.

Before the first channel, event-position, or recoil-key random draw, the
application preflights the largest `K*N` transaction in the complete parent
family plus the retained survivor. A rejection leaves particle bytes, linked-
list path, RNG state, physical selection/product diagnostics, and source
carrier unchanged. Eligible-visit and allocation-deferral counters still
record the attempted work and capacity pressure. Pending state
is 16 application-owned bytes copied by ordinary block/MPI particle transfer.
The event callback converts species at the event, applies the parent/product
`W/dt` ratio once, and rescales the unused interval by
`dt_product/dt_parent`. A zero-length re-query proves the marker was consumed
exactly once.

### New and modified files

Core files:

- `src/pic/pic.h` adds the generic all-or-none `QueryTrajectoryEvent`,
  `AdvanceTrajectoryEvent`, and `ProcessTrajectoryEvent` callback types,
  status/result enums, and ownership contract. The interface has no comet or
  chemistry dependency and defaults to three null pointers.
- `src/pic/pic_mover.cpp` defines the null callbacks, rejects partial
  installation during initialization and at mover entry, and integrates an
  event only into the legacy first-order traced mover. Query is side-effect
  free; a strictly earlier event wins over geometry, while exact face equality
  gives the face precedence. Event clocks advance before a face clone/MPI
  exchange. The core owns deletion and destination-list insertion and checks
  allocation, species, finite phase space/remaining time, positive correction,
  and marker clearing after application processing.
- `src/pic/README.md` documents callback timing, face-equality semantics,
  particle/list ownership, local-timestep normalization, and unsupported
  movers. With all callbacks null, no event call or RNG draw occurs.

Input and pure-policy files:

- `input/CometTotalReactionForcingPolicy.pm` is the strict species-policy
  resolver/generator. It accepts only `ON`, case-sensitive known species, and
  positive-rate unit-unimolecular parent families; it rejects duplicates,
  off-mode records, zero-rate/noneligible species, and an empty enabled mode.
- `input/test_comet_total_reaction_forcing_policy.pl` covers valid generation,
  default-off identity, and malformed/incomplete policies.
- `input/hartley2.pl` recognizes the exact records
  `ChemistryTotalReactionForcingMode=ON` and
  `ChemistryTotalReactionForcing(SPECIES:name)=ON`, writes the generated table,
  and requires chemistry, DIM1 radial mode, individual corrections, Bernoulli
  residuals, Stage-5 windows, diagnostics, and an explicit product cap.
- `srcHartley2/CometTotalReactionForcingPolicy.generated.h` is a safe source-
  tree fallback with all species unconfigured. The build copy is always
  regenerated from the selected input.
- `srcHartley2/TotalReactionForcingPlanner.h` is AMPS/MPI/RNG independent. It
  implements homogeneous and piecewise optical-depth plans, stable `expm1`/
  `log1p` inversion, exact carrier closure, multigeneration decay, and checked
  capacity arithmetic.

Application runtime files:

- `srcHartley2/CometChemistry.h` declares the opt-in mode, resolved-species
  query, mover installer, and initialized validator.
- `srcHartley2/CometChemistry.cpp` contains the generated-policy validation,
  conditional split transaction, compact pending record, deterministic recoil,
  three mover callbacks, diagnostics composition, capacity/RNG rollback, and
  Stage-8 runtime fixtures. It reserves pending bytes only in enabled builds.
  The pre-existing Stage-6 validator temporarily disables Stage 8 around its
  direct endpoint replay so the lower-stage gate remains unchanged. The final
  mover fixture creates one real C2 transaction with sixteen C product markers
  plus a survivor and adds one explicit C2 ionization sink. Product events are
  scheduled after the selected face; the sink is scheduled at the exact stable
  q-form face-intersection double. The fixture runs
  `PIC::TimeStepInternal::ExecutionTrackDefault`, audits callback owner/radius/
  remaining time from inside `ProcessTotalReactionEvent`, and scans final
  particles for species, marker drain, correction sum, and predicted endpoint.
  Particle-list insertion deliberately converts `FindCellIndex`'s ghost-
  offset center-node result to the compact physical-cell index; the initial
  smoke test exposed and fixed this distinction before acceptance.
- `srcHartley2/RadialSimulation.cpp` installs the complete callback set after
  particle-data offsets exist and before `PIC::Mover::Init()` validates it.
- `srcHartley2/RadialMain.cpp` adds the mutually exclusive
  `--validate-total-reaction-forcing` command and dispatch.

Test files:

- `srcHartley2/tests/test_stage8_total_reaction_forcing_planner.cpp` covers
  split closure, constant/piecewise hazards, conditional placement,
  multigeneration, capacity limits, extreme inputs, and rejection paths. It
  now also advances a closed-form discrete-cycle `A -> B -> C` network,
  verifies `A+B+C` carrier closure after every cycle, and compares 64 keyed
  analog ensembles against the deterministic-forcing result and analytic
  occurrence variance.
- `srcHartley2/tests/test_stage8_forced_segment_reference.py` is an independent
  analytic oracle for survival, conditional CDFs, and piecewise positions.
- `srcHartley2/tests/test_stage8_total_reaction_forcing_logs.py` strictly
  parses and recomputes the runtime records. Global analytic fixtures must be
  rank-identical; each rank-local transaction residue must independently meet
  its physical bound because `MPI_MAX` can select different valid roundoff
  residues at different decompositions. It also requires one production-
  mover record with 18 face transfers/continuations, 16 product events, one
  face-tied sink, destination ownership, 17 completed particles, no pending
  marker/remainder, and bounded carrier/event/endpoint errors.
- `srcHartley2/tests/test_stage8_total_reaction_forcing_parser.py` mutation-
  tests missing/wrong survival, event bounds, `P*L/N`, rollback,
  double-chemistry, transaction roundoff, missing local/MPI face crossing,
  uncleared markers, unfinished remainders, bad product endpoints, and bad
  mover correction. All 13 mutations are independently rejected.
- `srcHartley2/tests/run_stage8_total_reaction_forcing_tests.sh` runs the pure
  Perl/C++ gate without configuring or compiling AMPS.
- `srcHartley2/tests/run_stage8_integrated_tests.sh` replays the Stage-4 local-
  timestep, unit-correction, and weight-adaptation validators before running
  Stage-5/6/7/8 at exactly 1, 4, and 12 ranks into a nonexisting immutable
  root. It exposes separate mover-particle and chemistry-event counts and uses
  the existing strict Stage-4 parser rather than duplicating its schema.
- `srcHartley2/tests/run_all_stage_tests.sh` adds archived/fresh Stage-8
  modes, five jobs, strict provenance, and a total of 30 tests in nine groups;
  it still requires no `rg` executable.

Documentation files:

- `srcHartley2/CHEMISTRY.md` records the equations, policy, mover transaction,
  validation, limits, and Stage-9 science gate.
- `srcHartley2/README.md` gives operator input/CLI/test commands and archive
  rules.
- `srcHartley2/GENERIC_COMET_DESIGN.md` records the reusable core/application
  boundary and future geometry/mover work.
- `srcHartley2/STATUS.md` and root `handoff.md` contain recovery state and the
  ordered later-master integration procedure.

The temporary validation edits to `input/halley1d.input` are not source
changes. They are frozen only as
`validation_runs/stage8-total-reaction-forcing-full-20261105-r4/provenance/halley1d.stage8-validation.input`.
The tracked canonical input is restored to SHA-256
`0d89b536e6daa2dea7251bb57920fb6f1dd2546efd7c454d25b6dcccc26eccf5`.
Root `hartley2.pl`, root input symlinks/copies, `build/`, `amps`, PT output,
meshes, logs, and validation archives are generated/local artifacts and must
not be staged with the implementation.

### Validation evidence

Pure policy, parser syntax, planner, analytic reference, 13 mutation tests,
shell syntax, and `git diff --check` pass. The Stage-8-enabled application was
built with `-O3` only after removing exactly the repository `build/` directory.
Both enabled and final canonical builds emitted the pre-existing DIM1
aggressive-loop-optimization warnings in unrelated generic 3-coordinate core
paths; neither build emitted a Stage-8 compile or link error.
The accepted runtime archive is:

```text
validation_runs/stage8-total-reaction-forcing-full-20261105-r4
PASS Stage-4 correction-on integrated matrix 1rank=pass 4rank=pass 12rank=pass
PASS Stage-5 correction-window integrated matrix 1rank=pass 4rank=pass 12rank=pass
PASS Stage-6 daughter-packet integrated matrix 1rank=pass 4rank=pass 12rank=pass reactions=46 sinks=23 packet-cases=4 events=200000
PASS Stage-7 channel-importance integrated matrix 1rank=pass 4rank=pass 12rank=pass families=27 reactions=69 sidecar-rows=69 events=200000
PASS Stage-8 total-reaction forcing integrated matrix 1rank=pass 4rank=pass 12rank=pass events=200000 homogeneous=1 piecewise=1 decay=1 composition=1 rollback=1 mover=1
```

The archive manifest freezes the acceptance executable, validation input,
core configuration, generated chemistry table, and Stage-6/7/8 policy tables.
The enabled executable SHA-256 is
`822905e735deef438b9c0cbe31bb2ef150ee23b60d912def5d33d8f522f09ce6`;
the generated Stage-8 table is
`f3a788fdb0e65457b8fea06986d172d1f6cac84fecd79bc488cf73ee439db8fa`.

The `-r1` archive is failed development evidence because its summary counted
`K*N` product markers as reacted branches. The `-r2` archive fixed that
definition but exposed an invalid parser demand for identical rank-local
roundoff maxima; all physical bounds passed. The `-r3` archive passed the
direct Query/Advance/Process algebra but an audit found that it never passed
pending markers through `ExecutionTrackDefault` and did not replay Stage 4.
It is pre-acceptance evidence. Do not use, merge, rename, or substitute any of
`-r1`, `-r2`, or `-r3` for accepted `-r4`.

### Supported scope and unresolved work

Stage 8 is complete as implementation/focused-compatibility evidence only for
DIM1 spherical radial, the first-order force-free traced mover, and a constant
total rate during one logical local-timestep segment. The pure planner tests
piecewise hazards, but production does not yet integrate a varying hazard
across several blocks. Second-order, Boris, guiding-center, field-line,
periodic, CUDA, hybrid, DIM2, and DIM3 paths are not supported claims.

Stage 9 still requires multiple independent transported seeds comparing
analog and forced occurrence, measured variance/cost, collision feedback,
unchanged parent lifetime and primary profiles, daughter production/loss/
escape closure, buffer headroom, effective sample sizes, convergence, and a
new Rubin comparison. Validation values such as broad correction windows,
`N=8`, and equal target channel probabilities are fixtures, not Halley input
recommendations. The accepted matrix does replay Stage-4 corrected NTC
collision invariants with callbacks installed, but it does not yet collide a
newly transported Stage-8 daughter; that specific feedback remains Stage 9.

### Later merge to master: ordered procedure, not authorization

Do not merge yet. After Stage-9 science acceptance and explicit user approval,
perform integration on a disposable branch made from a recorded current
`master` tip. Never switch the present dirty development worktree to master.

1. Record the exact `1d_comet_dev` commit range, merge base, working-tree
   diff, excluded artifacts, and current master tip. Preserve unrelated master
   changes and the user-owned `Makefile` edit.
2. Integrate and test Stages 0--7 first. Stage 8 depends on individual
   corrections, correction windows, packet counts, channel likelihoods,
   diagnostics, and the established `K*N` capacity contract.
3. Merge the three generic core callback declarations/definitions and
   first-order mover hook as a separate reviewable change. Retain all-null
   defaults, all-or-none initialization validation, face-equality precedence,
   core-only deletion/list ownership, post-callback validation, and zero-event
   re-query. Resolve mover conflicts semantically; do not copy a whole old
   `pic_mover.cpp` over master.
4. Clean-build representative non-comet/core applications with null callbacks
   and run their mover tests before adding comet registration. Confirm no RNG
   draw and no behavior/configuration change in off mode.
5. Add the strict input policy, pure planner, generated fallback, and pure
   tests. Run omitted/off, valid C2, duplicate, unknown, zero-rate, bad value,
   and missing-prerequisite parser cases before runtime code.
6. Add application declarations and chemistry implementation, then the
   initialization registration and CLI. Preserve transaction preflight before
   RNG/mutation, 16-byte memcpy access, one `P`, one `L`, one `/N`, common
   packet recoil, `W/dt` conversion, remaining-time rescaling, sink handling,
   leader-only diagnostics, and rollback invariants.
7. Add strict analytic/mutation/runtime tests and runner wiring. The runtime
   gate must use the production execution track, not direct callback calls:
   require local and MPI faces, an exact face/event tie, advance-before-copy,
   destination-owner processing, mover-owned sink deletion, marker/remainder
   drain, and predicted product endpoints. Preserve the compact particle-list
   index when seeding fixtures; `FindCellIndex` returns a ghost-offset center-
   node index. Do not weaken exact analytic rank identity; only rank-local
   transaction maxima use independently bounded roundoff.
8. Create an independent validation input outside the canonical file, remove
   exactly the integration worktree's `build/`, build Stage 8 on, freeze all
   provenance, and rerun the complete 1/4/12 Stage-4/5/6/7/8 matrix into a new
   nonexisting archive. Stage 4 must be replayed with the Stage-8 callback set
   installed so legacy face, MPI continuation, correction, collision, and
   adaptation paths cannot regress silently.
9. Restore the canonical off input, remove exactly `build/` again, clean-build
   and verify correction/packet/importance/forcing off. Run the complete
   30-test, nine-group orchestrator and relevant core/application regressions.
10. Run the Stage-9 transported science plan; do not promote validation
    windows, packets, `q*`, or forcing species into a production input based
    on visual agreement with Rubin.
11. Exclude `Makefile`, `build/`, `amps`, root generated parsers/inputs,
    meshes, `PT*`, logs, raw validation archives, failed `-r1`/`-r2`, and pre-
    acceptance `-r3` evidence from source commits unless the user explicitly
    requests a curated small manifest.
12. Present the exact integration diff/commits, exclusions, hashes, complete
    test matrix, performance/memory/uncertainty/closure/convergence results,
    and Rubin products to the user. Stop and obtain explicit approval before
    any merge, push, or operation affecting `master`.

### Final canonical workspace after Stage-8 validation

After freezing the enabled executable, `input/halley1d.input` was restored,
exactly `build` was removed, and the application
was clean-built with `make -j8 EXTRACOMPILEROPTIONS=-O3`. The generated core
macro is `_INDIVIDUAL_PARTICLE_WEIGHT_OFF_`; all 27 total-forcing flags are
zero; all 69 packet counts are one; and all 69 channel-importance flags are
zero. Final hashes are:

```text
7752f0967512c4de136c5ef9f284a682bbfc986aab6a89b8e2459262c6a96317  amps
0d89b536e6daa2dea7251bb57920fb6f1dd2546efd7c454d25b6dcccc26eccf5  input/halley1d.input
c3d2bca73f068cff357ae6d1f5002c97dc738b6a3190c712b2787890c92e2726  build/pic/picGlobal.dfn
ce925a583f1b1cffa641ac8336fb59eba850605ecb6730df35e09fa39507910b  build/main/RadialSimulation.cpp
311c78deeb0fcc2140f21d995d86b5024ac2c9c9b3091a46631e5caef055413a  build/main/CometTotalReactionForcingPolicy.generated.h
419fda57d23876a1837ab9b9f4867231331c05b05ba5fabd631d1c1184299fc1  build/main/CometDaughterPacketPolicy.generated.h
7f4e0be8f6a6f09e0b39a28642200839a1d32a7d195953befd043609632a9fbf  build/main/CometChannelImportancePolicy.generated.h
```

The final archived replay passed without skip or failure:

```text
tests: run=30 passed=30 skipped=0 failed=0
test-groups: run=9 passed=9 skipped=0 failed=0
PASS Hartley2 all-stage suite tests=30 test-groups=9 stage4-mode=archived stage5-mode=archived stage6-mode=archived stage7-mode=archived stage8-mode=archived
```

No AMPS or MPI process remains active. No commit or merge was made.

## Stage 9 transported-science implementation (2026-11-05, focused validation complete)

Stage 9 now has the application instrumentation and strict evidence readers
needed to compare the analog and fully enhanced estimators from one frozen
executable. This section records implementation state, not science closure.
The balanced long ensemble and the regenerated Rubin comparison have not yet
been accepted. Development remains on `1d_comet_dev`, based on Stage-8 commit
`c87bdcd0c7d2d06b0795bd7d5ec9e1a422f9e22b`. No branch switch, commit,
merge, push, or operation on `master` is part of this work.

### Modified application and input files

- `input/hartley2.pl` parses `Stage9ScienceDiagnostics=ON|OFF`, writes the
  generated `Stage9ScienceDiagnosticsEnabled` value, and fails closed unless
  an enabled case is DIM1 radial chemistry with comet diagnostics, individual
  corrections, Bernoulli residuals, correction windows, daughter packets,
  channel importance, and total-reaction forcing. The default is off.
- `srcHartley2/CometChemistry.h` declares the runtime estimator modes and
  accessors `ConfigureEstimatorMode`, `GetEstimatorModeName`,
  `IsDaughterPacketSplittingActive`, `IsChannelImportanceSamplingActive`, and
  `IsTotalReactionForcingActive`. It also declares the read-only
  `HasPendingTotalReactionEvent` safe-boundary audit and the optional
  transported-collision extension to `ValidateTotalReactionForcing`.
- `srcHartley2/CometChemistry.cpp` implements `configured`, `analog`, and
  `enhanced` modes without modifying the compiled enable flags. Analog uses
  `N=1`, physical channel selection, likelihood one, and ordinary stochastic
  total occurrence; enhanced uses the compiled Stage-6/7/8 policies. Both
  retain identical particle bytes and callbacks. The file also extends the
  live Stage-8 mover fixture so actual transported C daughters from `C2_2C`
  collide with H2O through production `ModelCellCollisions_ntc`. It audits
  exact transported records, ownership, feedback to both velocity
  populations, collision probabilities, correction/rate retention, completed
  trajectories, and represented momentum.
- `srcHartley2/RadialMain.cpp` adds
  `--chemistry-estimator-mode configured|analog|enhanced` and
  `--validate-transported-daughter-collision`. Estimator selection is applied
  after compiled input initialization and before validation or particle
  storage; a leader-only provenance record reports the active and compiled
  policies.
- `srcHartley2/RadialSimulation.h/.cpp` add the opt-in Stage-9 checkpoint
  writer. A per-iteration local peak update avoids an MPI collective in the
  hot loop. Complete sampling windows reduce carrier sums, squared carrier
  sums, counts, correction extrema, capacity, utilization, wall time, and
  safe-boundary audit counts. Carrier ESS uses
  `q=W_block,species*c_particle`, so local block-weight adaptation is included.
  Safe-boundary scans abort on an unfinished local-timestep fraction, pending
  Stage-8 event, or linked-list/live-buffer disagreement before output,
  adaptation, or load balancing.

No Stage-9 change is made to `src/pic`. The transported-collision test calls
the existing production general-NTC implementation, and the same-binary arms
use the previously accepted Stage-8 null-default trajectory-event contract.
The unrelated tracked `Makefile` linker edit remains user-owned and must not
be folded into this stage.

### New test and evidence files

- `srcHartley2/tests/run_stage9_integrated_tests.sh` runs the immutable
  transported-daughter collision fixture at exactly 1, 4, and 12 ranks into
  a new nonexisting root. It does not build, edit input, or overwrite evidence.
- `srcHartley2/tests/test_stage9_transported_collision_logs.py` strictly
  parses the focused collision records and compares only deterministic
  chemical/carrier invariants across rank counts; RNG-dependent trial counts
  and owner rank are bounded but not required to be identical.
- `srcHartley2/tests/test_stage9_transported_collision_logs_selftest.py`
  mutation-tests that parser.
- `srcHartley2/tests/run_stage9_science_gate_tests.sh` runs the dependency-free
  gate self-tests and optionally evaluates a real manifest.
- `srcHartley2/tests/test_stage9_science_gate.py` reads checksum-protected
  `amps.comet.stage9-science-manifest.v1` and
  `amps.comet.stage9-run-summary.v1` records. It requires balanced analog and
  forced arms, at least four seeds per arm, fresh starts, identical physical
  provenance, convergence/burn-in/windows, capacity and pending-state
  evidence, collision feedback, primary invariance, daughter uncertainty and
  ESS, zero-cell fractions, flux identities, and cost-normalized variance.
- `srcHartley2/tests/test_stage9_science_gate_selftest.py` supplies valid and
  deliberately corrupted synthetic archives. Those fixture thresholds are
  parser tests and are not Halley science choices.
- `srcHartley2/validation/stage9_run_summary.py` is the production summary
  verifier. It reads completed Stage-0 chemistry
  reports, Stage-9 checkpoint/carrier files, and unambiguous AMPS profiles; it
  refuses to repair, interpolate, or infer missing science data.
- `srcHartley2/validation/run_stage9_science_campaign.py` provides immutable
  `plan`, `run`, `resume`, and read-only `status` operations for four fixed
  seeds in balanced analog/enhanced arms. Every case starts from iteration
  zero; resume creates a new attempt rather than reusing a partial PT tree.
  Enhanced cases first run the focused transported-daughter collision gate in
  an isolated directory, avoiding a new observer in the NTC hot path.
- `srcHartley2/validation/stage9_science_policy.template.json` freezes the
  12-rank, one-thread, 10,000-iteration burn-in, three 10,000-step accepted
  windows, 40,000 maximum iterations, and four documented seeds. Its ESS,
  flux-closure, variance/cost gain, and bin-coverage entries remain literal
  `DECISION_REQUIRED`; campaign planning fails until they are resolved.
- `srcHartley2/validation/stage9_tooling_selftest.py` and
  `run_stage9_tooling_tests.sh` mutation-test summary construction, complete
  provenance freezing, balanced case ordering, interruption recovery, both
  flux-closure policies, and publication of an accepted manifest.
- `srcHartley2/validation/generate_stage9_rubin_comparison.py`,
  `stage9_comparison_selftest.py`, and `run_stage9_comparison_tests.sh`
  authenticate an accepted balanced manifest before producing non-overwriting
  PNG, EPS, PDF, CSV, and JSON ensemble products with Rubin overlays.
- `srcHartley2/tests/run_stage9_parser_tests.sh` and
  `test_stage9_science_diagnostics_parser.py` cover default/off/on parsing,
  invalid tokens, and every required Stage-5/6/7/8 prerequisite.
- `srcHartley2/tests/test_stage9_checkpoint_outputs.py` independently
  recomputes every per-cell carrier count, `sum(q)`, `sum(q^2)`, Kish ESS,
  correction range, and species aggregate from paired analog/enhanced smoke
  outputs. `test_stage9_checkpoint_outputs_selftest.py` mutation-tests that
  fail-closed contract.
- `srcHartley2/tests/run_all_stage_tests.sh` now includes all Stage-9 parser,
  science-gate, campaign, comparison, checkpoint, checksum, and collision
  groups without adding an `rg` runtime dependency.

### Output schemas and safe interpretation

An enabled science input writes

```text
PT/diagnostics/comet_stage9_checkpoint.out=<N>.csv
PT/diagnostics/comet_stage9_checkpoint.out=<N>.json
PT/diagnostics/stage9_carriers.<SPECIES>.out=<N>.csv
```

The checkpoint schema is `amps.comet.stage9-checkpoint.v1`. It records the
completed iteration, output index, sample length, mode, seed, ranks, maximum
sampling-window wall time, live population/capacity, peak rank utilization,
safe-boundary audits, zero pending remainders/events, and per-species physical
carrier inventory/ESS/correction range. Per-cell tables preserve cell index,
carrier count, `sum(q)`, `sum(q^2)`, ESS, and provenance. An empty cell has
zero ESS and blank correction extrema; the writer does not manufacture a
carrier for it.

The focused collision log is a compatibility proof, not a measured physical
collision rate. The controlled fixture chooses its carrier/time-step scale so
the two update probabilities are analytically one and can expose correction-
aware NTC wiring errors without statistical ambiguity. Long-run science must
use the case's unmodified physical weights, timesteps, and collision data.

### Accepted focused validation and exact evidence

The Stage-9-enabled build was created only after removing the exact
`build` directory. The accepted transported-
daughter archive is:

```text
validation_runs/stage9-transported-collision-full-20261105-r3
```

Its frozen executable SHA-256 is
`c10fca52f5be1c8bb6da79d1ac1daab6167422aeb5f5b440134fedb794d5efaf`
and its validation-input SHA-256 is
`8e3599d62b57b712075043ec284eca1a5d489aa1fb37b079e8ad56cb00c917d6`;
the generated-header hashes are in `provenance/SHA256SUMS`. The strict
1/4/12-rank matrix passed with 16 transported C daughters, 64 H2O partners,
production general NTC, rate closure, momentum closure, and zero unfinished
state. Earlier `-r1` and `-r2` directories are failed evidence and must not be
substituted for `-r3`.

The paired one-window same-binary smoke archive is:

```text
validation_runs/stage9-checkpoint-smoke-20261105-r2
```

Both arms used executable SHA-256
`75cedeb999287868dd4cc7a5a273de8967de7c848cc9b3d31abf47954136ff47`,
input SHA-256
`ba8c812616a6f47a64652e5ff35eeceb60625bfd07221ea268a62231a48054da`,
seed 104729, one rank, one completed iteration/window, 27 species, and 2,984
radial cells. The independent output checker passed with zero pending
remainders/events and identical analog/enhanced layout. The smoke is
deliberately nonconverged and is only an output-contract test.

Focused mutation/parser results were: collision parser 9/9, science gate
15/15, campaign/summary tooling 9/9, comparison tooling 5/5, input parser
13/13, and checkpoint validator 10/10. After restoring `halley1d.input`, the
exact build directory was removed again and the canonical configuration was
rebuilt with individual corrections, Stage-5 windows, Stage-6 packets,
Stage-7 channel importance, Stage-8 forcing, and Stage-9 diagnostics all off.
The complete aggregate result was:

```text
tests: run=38 passed=38 skipped=0 failed=0
test-groups: run=10 passed=10 skipped=0 failed=0
PASS Hartley2 all-stage suite tests=38 test-groups=10 stage4-mode=archived stage5-mode=archived stage6-mode=archived stage7-mode=archived stage8-mode=archived stage9-mode=archived
```

Logs are preserved under
`validation_runs/stage9-all-stage-tests-20261105`.
This closes Stage-9 implementation/focused validation, not its long-run
scientific gate. Do not describe a synthetic tooling PASS, the one-window
smoke, or the deterministic collision fixture as Halley/Rubin acceptance.

Two build-provenance hazards were found and must be checked during any replay:
`.amps.conf` can override a command-line input selection, and the root
generated `hartley2.pl` can be stale relative to tracked `input/hartley2.pl`.
Verify the effective `InputFileAMPS`, synchronize/regenerate the root parser,
and hash the executable/input/generated headers before launching evidence.

### Scientific decisions still required

The evaluator intentionally does not silently choose:

1. a minimum per-daughter carrier-ESS fraction;
2. a minimum variance-per-wall-time gain and required fraction of resolved
   radial bins achieving that gain; or
3. a relative flux-closure tolerance. The preferred alternative for item 3
   is to require the independent-seed confidence interval of the net
   source-plus-production-minus-loss-minus-escape imbalance to contain zero.

Existing Stage-3 profile/zero-cell thresholds may be retained only when their
meaning and comparison interval are unchanged. Stage-6 packet counts, Stage-7
target probabilities, Stage-8 forcing choices, and synthetic self-test values
are not substitutes for these decisions.

### Later merge to master: Stage-9 additions, not authorization

Do not merge yet. After the full science gate passes and the user explicitly
authorizes integration, append these steps to the Stage-8 master procedure:

1. Freeze the `1d_comet_dev` commit range, merge base, Stage-9 executable,
   validation input, generated policy headers, mesh/network/physics hashes,
   focused collision archive, balanced seed manifest, summaries, evaluator
   result, and Rubin comparison products. Record excluded local artifacts and
   preserve the unrelated `Makefile` edit.
2. Integrate the already reviewed Stages 0--8 first. Re-run the null-callback
   core mover regressions before Stage-9 application code; Stage 9 adds no new
   core API and must not be used to conceal a Stage-8 core conflict.
3. Merge estimator-mode declarations/implementation and CLI as one
   application change. Verify the default `configured` path first, then prove
   analog and enhanced reserve identical particle bytes, install identical
   callbacks, load identical generated tables, and reject selection after
   storage reservation.
4. Merge the input switch and checkpoint writer together. Replay parser cases
   for omitted/off/on, every missing prerequisite, and invalid values. Verify
   off mode adds neither files nor scans; verify an enabled short case writes
   schema-v1 JSON/CSV with exact seed/rank/mode provenance and zero pending
   states.
5. Merge the transported-collision fixture and strict parser. Run 1/4/12 MPI
   from one frozen executable, preserving production general NTC and the real
   Stage-8 mover route. Do not weaken pointer, owner, probability, carrier,
   momentum, or drain checks to accommodate a decomposition.
6. Merge the summary builder, manifest gate, mutation tests, campaign tooling,
   comparison exporter, checkpoint validator, and aggregate-runner wiring.
   Verify the staged file list explicitly and keep raw PT trees and large
   archives outside source.
7. Clean-build the canonical correction-off input after removing exactly the
   integration worktree's `build/`. Run every Stage 0--9 focused group plus
   representative master/core applications and confirm default-off behavior.
8. Reproduce the balanced analog/enhanced seed matrix from iteration zero.
   The manifest must authenticate identical physics and numerical controls
   except estimator mode and seed. Require convergence, accepted burn-in and
   windows, collision feedback, buffer headroom, ESS/profile/flux gates, and
   measured cost-normalized gain under the approved policy.
9. Regenerate Rubin PNG, EPS, PDF, CSV, and JSON products from accepted
   summaries with seed confidence intervals and convergence/cost metadata.
10. Present exact diffs, commits, hashes, test counts, failed/skipped logs,
    performance, uncertainty, closure, convergence, Rubin products, and all
    excluded artifacts. Stop for explicit approval before any merge, push, or
    operation affecting `master`.

Never copy an entire older `CometChemistry.cpp`, `RadialSimulation.cpp`, or
`hartley2.pl` over master. Resolve conflicts semantically around current
master interfaces. Never stage `build/`, `amps`, root generated inputs or
parser copies, meshes, `PT*`, transient logs, temporary validation inputs, or
raw validation archives merely because they exist in the worktree.

## Stage 10 merge-readiness evidence (2026-11-05, development branch only)

Stage 10 prepares a reproducible review package for the generic-comet work.
It does not change comet physics, chemistry, particle transport, weights,
timesteps, collisions, source rates, or output quantities. It also does not
switch branches, commit, merge, rebase, cherry-pick, push, or update
`master`. Development remains on `1d_comet_dev`, whose committed starting tip
for this stage is
`bff4c532f34515365fd910691595f261db96847f`.

### Repository topology and authority boundary

The read-only audit found two different possible master targets:

```text
local master:       cb7be2bba41f373352ab853abed23955ea69c8b9
origin/master:      635d1622f80b445cf83130e90e5539bfa872e453
divergence:         1 local-only commit, 17 origin-only commits
local merge base:   cb7be2bba41f373352ab853abed23955ea69c8b9
origin merge base:  f9c2c00da9c8348b9d4003e2a82c1fe139b300cc
```

Consequently no Stage-10 command infers what “master” means. Collection
requires an explicit allowed target ref, its full commit, the full merge-base
commit, and the full development commit. The report always serializes
`merge_authorization=NOT_GRANTED`. Even a future
`READY_FOR_USER_REVIEW` result is evidence for a user decision, not authority
to perform a Git operation.

The current scientific state is intentionally separate. The committed Stage-9
implementation has 38/38 focused tests in 10/10 groups, but the balanced
independent-seed Halley campaign and accepted Rubin products do not exist.
The minimum daughter carrier-ESS fraction, variance-per-walltime gain,
required resolved-bin fraction, and flux-closure policy remain explicit
decisions. Therefore the only valid current readiness interpretation is
`FOCUSED_ONLY_NOT_SCIENCE_READY`.

### New Stage-10 implementation files

- `srcHartley2/tests/write_all_stage_result.py` converts the shell runner's
  exact test records into schema `amps.comet.all-stage-test-result.v1`. It
  validates full Git object IDs, Stage-4--9 execution modes, unique test IDs,
  result/exit consistency, log containment, group/test totals, optional
  science-test coupling, and SHA-256 values. Publication is atomic and refuses
  overwrite.
- `srcHartley2/tests/test_stage10_all_stage_result.py` mutation-tests JSON
  escaping, focused/science separation, duplicate IDs, escaped log paths,
  inconsistent exits/results, missing science records, and overwrite
  refusal.
- `srcHartley2/tests/run_stage10_tests.sh` is the one-command dependency-free
  Stage-10 test group. It runs both the receipt and readiness mutation suites
  and labels its PASS as tooling-only.
- `srcHartley2/validation/stage10_release_policy.json` versions the allowed
  target refs, sole tracked-dirty exclusion, exact focused-suite contract,
  science schemas/minima, path classifications, and immutable
  `NOT_GRANTED` authorization state.
- `srcHartley2/validation/stage10_merge_readiness.py` uses only an explicit
  allowlist of read-only Git queries. `collect` verifies commits, branch,
  merge base, target ref, tracked dirty exclusions, the complete JSON test
  receipt and its logs, exact integration diffs, and any optional Stage-9/
  Rubin evidence before creating a new checksum-protected bundle. `verify`
  checks every member, canonical Markdown rendering, policy binding, status,
  blockers, and authorization without consulting or modifying Git.
- `srcHartley2/validation/stage10_merge_readiness_selftest.py` creates
  disposable repositories and evidence, then mutation-tests topology,
  receipt/log checksums, classifications, dirty exclusions, output
  immutability, forged/missing science, bundle tampering, and authorization.
- `srcHartley2/validation/run_stage10_merge_readiness_tests.sh` provides the
  dependency/tool preflight and invokes the self-test without writing Python
  bytecode into the source tree.
- `srcHartley2/MASTER_MERGE.md` is the tracked manual semantic-integration
  guide. It records the authority stops, explicit target selection, ownership
  and dependency waves, invariants, allowlist/exclusions, clean-build/test
  replay, science blockers, and the final evidence that must be shown to the
  user. It is not an executable merge script.

### Modified files and rationale

- `srcHartley2/tests/run_all_stage_tests.sh` retains its human summary and
  now writes `all_stage_result.json` only for an explicit new `--log-root`.
  The receipt binds the test list, log hashes, counts, modes, branch, commit,
  tree, and separate focused/science meanings. A supplied real Stage-9
  manifest is evaluated as its own test instead of being conflated with the
  synthetic science-gate mutation suite. Stage 10 is one additional group,
  making the focused no-manifest contract 39 tests in 11 groups.
- `srcHartley2/CHEMISTRY.md`, `README.md`,
  `GENERIC_COMET_DESIGN.md`, and `validation/README.md` document the evidence
  states, commands, schema ownership, default science blocker, and strict
  non-authority boundary.
- `handoff.md` now begins with the current branch/commit/topology, Stage-9
  commit state, Stage-10 recovery state, unresolved science decisions, and
  the prohibition on any master operation.
- `src/models/chemical_solver/README.md` removes one EOF whitespace finding
  exposed by the complete integration-range audit. No technical text changed.
- `srcHartley2/halley2011/reference/rubin2011_neutral_profiles.csv` is
  normalized from CRLF to LF because the committed range otherwise reports
  every reference row as trailing whitespace. The data fields are unchanged;
  the normalized SHA-256 is
  `a75b28ac9141255b2ee9ce4f5f1e1abbb8c1332b7d9e16d79fcae2eab0d19ace`.

The unrelated tracked root `Makefile` linker edit remains user-owned and is
not part of Stage 10. The release policy permits it only as an explicitly
hashed dirty exclusion. Untracked `build/`, `amps`, `PT*`, meshes, root input
and parser copies, logs, generated tables, and `validation_runs/` remain
outside merge scope.

### Later integration procedure

After Stage 9 science is accepted and the user separately authorizes work
against one exact target commit, follow `MASTER_MERGE.md` rather than copying
entire development files. Integrate reusable chemistry/cometary-physics
interfaces first, core particle/mover changes with null/default-off behavior
second, application input/generated policies next, then DIM1 application
consumers, tests/evidence tooling, reference inputs, and documentation. After
each semantic wave, replay its focused tests. Remove exactly the disposable
integration worktree's `build/` before every recompile; run representative
target-owned non-comet/DIM3 regressions as well as all Stage 0--10 groups and
the live MPI matrices. Present exact diffs, exclusions, hashes, performance,
memory, independent-seed statistics, convergence/closure, and Rubin products,
then stop again for approval. No current result authorizes these operations.

### Stage-10 focused validation result

The dedicated Stage-10 runner passed both dependency-free mutation suites:

```text
PASS Stage-10 all-stage result writer cases=8
PASS Stage-10 merge-readiness selftest
  cases=target-ref,branch,receipt,untracked-source,artifact-exclusion,
        dirty-exclusion,extra-dirt,overwrite,tamper,forged-science,
        unknown-path,read-only-git
```

The complete archived Stage 0--10 suite then passed:

```text
tests:       run=39 passed=39 skipped=0 failed=0
test-groups: run=11 passed=11 skipped=0 failed=0
Stage 4--9 modes: archived
Stage-9 science:  NOT_EVALUATED
```

Durable evidence is under:

```text
validation_runs/stage10-all-stage-tests-20261105-r3
```

The human summary SHA-256 is
`549fd51c44a18512e1dd18e3007e57280e10d761fbe742aecf68fc9f8e90baab`.
The schema-v1 receipt SHA-256 is
`ca55abc4e319b464a5600ed10d8a33be9b8fb7a042f169b3bb129a5190d4a7d5`.
The canonical Stage-10 release-policy SHA-256 is
`8fd7c4b76b45c03d00422bffd9ca0ef0eb628bc4cc4c4fe3ca3fcabb737570c2`.
It records committed HEAD
`bff4c532f34515365fd910691595f261db96847f` and tree
`7ede9f010babbb3146dfcc072405caa84c8a5291`.

This archive proves the current Stage-10 worktree passed focused tests. It is
not a final merge-readiness bundle because the Stage-10 source/documentation
is not yet committed, so its development commit cannot describe those
changes. After a later user-requested Stage-10 commit, rerun the suite into a
new directory and run `collect` against that new exact commit. The collector
must reject dirty or untracked reviewable source and the current missing
Stage-9 science/Rubin evidence; do not weaken those checks to reuse this
development receipt.

A real-workspace `collect` preflight was also run deliberately before commit.
It exited 2 and named every non-allowlisted tracked Stage-10 modification
rather than producing a bundle; the requested output directory was not
created. This is the expected fail-closed behavior, not a test failure.

No C++ or generated AMPS configuration changed in Stage 10, so no compile was
required and the existing `build/` directory was not touched. The standing
rule remains: remove exactly `build` before any
future recompile.

## Input-driven DIM1 brightness development

The generic radial application now has an optional, default-off resonance-
fluorescence output layer. `CometBrightnessMath.h` owns the PIC-independent
exact spherical-shell chord, opaque-nucleus occultation, energy-to-photon and
Rayleigh conversions, and line-center optical-depth diagnostic.
`CometBrightness.cpp/.h` own the
input-generated line configuration, owner-cell MPI reduction and separate
completed-window CSV writer. `RadialSimulation.cpp` invokes that writer only
after the full-sampling-window guard and before population adaptation;
`RadialMain.cpp` validates generated configuration before mesh construction.

`input/CometBrightnessConfig.pm` strictly parses the grid and emission
records. Species and transition, a positive fixed photon g factor or named
core model, provenance, the explicit `THIN` mode, and the `ISOTROPIC_4PI`
phase convention are required. Core-model records require bulk radial
velocity; H requires case solar irradiance; OH energy-to-photon conversion
requires an explicit wavelength. Atomic-H records must also supply a
positive line-center cross
section, tau flag threshold, and independent cross-section provenance. The calculated
thin brightness is never attenuated; the output states
`attenuation_applied=false` and reports tau and the threshold flag separately.
No case H or OH g factor is hardwired. The reusable core module owns the
published reference models and data. `input/hartley2.pl` emits selectors,
conditions, arrays, and metadata and links `CometBrightness.o` only in the
DIM1 object set.
Existing inputs that omit `BrightnessMode` remain off, so the Rubin/Halley
comparison data contract is unchanged.

Each CSV row explicitly records `radial_domain_outer_radius_m`, states that
the column stops at that boundary, and labels the outer continuation `NONE`.
The writer never silently extrapolates a density tail beyond the DIM1 mesh.

Focused evidence currently passes:

```text
PASS comet-brightness exact-chord=1 occultation=1 units=1 tau-flag=1 invalid=1
PASS resonance-fluorescence core evaluator
Ran 3 reference-data tests: OK
PASS comet-brightness input policy fixed=1 modeled=1 default-off=1 invalid=11
Ran 3 parser integration tests: OK
PASS comet-brightness focused test group tests=7
```

The test runner is `srcHartley2/tests/run_comet_brightness_tests.sh`; it does
not touch `build/` or launch AMPS. A clean final case-1 configure/compile now
passes and the generated numerical arrays were inspected as described in the
I3 section above. A short run reaching a complete sampling window and writing
the CSV remains pending. Before any further recompile, remove exactly
`build/` as required by the standing clean-build policy.

For later integration with `master`, merge the reusable core model, reference
manifest/data, tests, and targeted root build wiring first. Then merge the
pure brightness math header/tests, parser module/application adapter, and the
two small radial call-site changes. Preserve default-off and fixed-scalar
behavior and the separate CSV schema. Do not replace the DIM1 shell sum with
the core Cartesian column
integrator: that core path currently exits for `DIM!=3`. Do not merge a case
campaign g factor, solar flux, OH photon wavelength, or attenuation
approximation into the generic C++ source; those remain explicit
input/provenance decisions.

## Particle-resolved fluorescence sampling revision (2026-11-07)

This section supersedes every earlier statement in this file that describes
brightness as a completed-window shell-density CSV or a bulk-velocity-only
calculation. Brightness is now evaluated inside the AMPS sampling lifecycle
from individual model particles and written in Tecplot format.

The DIM1 mover representation was verified before implementing the inverse
mapping. `src/pic/pic_mover.cpp` advances the particle with all three Cartesian
velocity components, detects the radial block/sphere intersection, and then
rotates both position and velocity so the stored position is `(r,0,0)` and
the stored velocity is expressed in the local radial basis. The new
`RotateLocalVelocityToCartesian()` in `CometBrightnessMath.h` is the exact
transpose/inverse of that normalization. Its regression uses nonzero radial
and both transverse components and checks a round trip to floating-point
tolerance.

For an impact parameter `b` and LOS coordinate `s`, the sampler uses
`x=b*IMPACT+s*LOS`. It intersects the ray with every occupied owner-cell
spherical shell, removes the occulted far segment when `b<R_n`, and evaluates
each remaining segment with two-point Gauss--Legendre quadrature. At each
quadrature point, every relevant model particle in the shell contributes

```text
(W_block * individual_correction / shell_measure) * ds
```

to the column. Its local velocity is rotated to `x/|x|`, combined with the
full comet heliocentric velocity, and projected as

```text
v_helio,r = e_r dot (V_comet + R(x/|x|) v_local).
```

Named H/OH models receive that particle-specific Doppler velocity. The fixed
g-factor path uses the same particle column estimator but skips the model
evaluation. Static ray data (`b`, `b^2`, and the sky-plane offset) are
precomputed once; AMR shell intersections remain owner-cell calculations so
they stay valid after load redistribution.

`V_comet` defaults to exactly `(0,0,0)` and can be changed through
`SetCometHeliocentricVelocity()`. Enabled input requires explicit normalized
heliocentric radial, comet-to-observer LOS, and sky-plane impact directions.
The scalar legacy velocity is shorthand for `V_comet=Vr*e_r`; a full
three-component input is also accepted. Out-of-domain g-factor evaluations
are never clamped or extrapolated. `ABORT` reports the line, model, velocity,
and status before terminating; `EXCLUDE_AND_REPORT` records excluded column
and evaluation count.

The sampler is registered once after PIC/radial initialization. The core
invokes it every sampling iteration and calls its output callback at the same
boundary as ordinary AMPS sampling output. The callback MPI-reduces the
buffers, divides by `PIC::LastSampleLength`, writes rank-zero Tecplot files,
clears all ranks' external buffers, and synchronizes. There is one file per
configured species; multiple transitions are columns:

```text
PT/plots/pic.<species>.brightness.s=<species-index>.out=<output-index>.dat
```

The variables include impact parameter, column density, effective photon g
factor, R/kR brightness, line-center tau/flag, and excluded-domain fraction
and count. Dataset auxiliary metadata records estimator, geometry, comet
velocity, policy, model/provenance, and exact reference-data hash.

Files changed specifically for this revision:

| File | Revision and later-merge purpose |
| --- | --- |
| `srcHartley2/CometBrightnessMath.h` | Adds the PIC-independent inverse DIM1 velocity rotation; merge with its round-trip unit test before the runtime sampler. |
| `srcHartley2/CometBrightness.h/.cpp` | Replaces the old density CSV writer with the particle external sampler, public geometry/velocity API, strict domain policy, static ray geometry, MPI lifecycle, and Tecplot writer. Treat this as one semantic unit. |
| `srcHartley2/RadialSimulation.cpp` | Registers the sampler after PIC initialization and removes the old manual completed-window writer call. Merge only after the target core callback signatures are checked. |
| `input/CometBrightnessConfig.pm` | Adds vector/scalar velocity resolution, nine-component geometry, and velocity-domain policy with duplicate and orthogonality checks. |
| `input/hartley2.pl` | Parses/generates the new arrays and policy. The policy declaration must stay on one physical C++ source line because `ChangeValueOfVariable()` replaces the complete declaration line. |
| `srcHartley2/i3/reference/i3_cases.csv` | Stores per-case `e_r`, comet-to-Mars LOS, and projected antisolar impact vectors; scalar table radial velocities remain authoritative. |
| `srcHartley2/i3/prepare_i3_cases.py`, `input/i3.input`, `input/i3.case*.input` | Emit explicit geometry and `ABORT` policy without hardwiring a campaign case into C++. Regenerate all six together; do not hand-edit one concrete case. **Superseded 2026-11-08:** the template and six case inputs are retired and will not be used; the generator now writes the geometry and policy into the case sections of `input/i3.in`. |
| `srcHartley2/tests/test_comet_brightness_math.cpp` | Tests the mover/sampler velocity-transform inverse in addition to chord, occultation, units, and tau policy. |
| `srcHartley2/tests/test_comet_brightness_parser.py`, `input/test_comet_brightness_config.pl` | Assert exact generated vectors and one-line policy declaration and reject missing/nonorthogonal geometry and invalid policy. |
| `srcHartley2/README.md`, `srcHartley2/i3/README.md`, `srcHartley2/MASTER_MERGE.md`, `handoff.md` | Replace obsolete bulk/CSV descriptions and document input, estimator, Tecplot schema, limitations, recovery, and integration order. |

Focused evidence after this revision is:

```text
PASS comet-brightness rotation=1 exact-chord=1 occultation=1 units=1 tau-flag=1 invalid=1
PASS resonance-fluorescence core evaluator
Ran 3 reference-data tests: OK
PASS comet-brightness input policy fixed=1 modeled=1 default-off=1 invalid=13
Ran 3 parser integration tests: OK
PASS comet-brightness focused test group tests=7
```

A clean case-1 preprocess/compile was performed after deleting exactly
`build`; the generated policy-declaration defect
found by the first attempt was fixed and covered by the parser regression.
Runtime initialization subsequently passed brightness configuration
validation and built the 376-block/3,008-cell I3 mesh. The selected
one-rank weight-adaptation validator later failed its own multi-rank
ghost/boundary coverage checks; that is not a brightness failure and is not
claimed as a passing adaptation test.

Still required before scientific use are a short completed-window runtime
that produces and parses the new Tecplot files, one-/multi-rank estimator
equivalence, an explicit runtime out-of-domain policy fixture, performance
and peak-memory measurement for the 200-ray campaign, and converged I3
production comparisons. The direct reference estimator intentionally favors
traceability over an unvalidated velocity-histogram shortcut.

For a later master integration, merge the reusable core fluorescence module
first, then the pure math/API and tests, then parser/generator, then the
runtime sampler plus the one radial registration change, and only then the
I3 campaign files. Preserve the Halley `BrightnessMode=off` input and rerun
the Rubin workflow. Do not retain the removed manual writer call, do not add
the particle sampler alongside the old density writer, and do not change the
filename/metadata contract during conflict resolution. No merge, commit,
branch switch, or master modification was performed here.

## I3 case-1 12-rank production start (2026-11-07)

*Historical: this run and the regenerated per-case inputs belong to the
retired six-file layout. New runs use `-i input/i3.in -case N`; see the top of
this file.*

At the user's direction, the maintained I3 brightness grid was reduced from
2,000 to 200 logarithmic rays. The common six-case template now uses a
10,000-iteration sampling length and a 30,000-iteration maximum. All six
concrete inputs were regenerated, the five generator tests and seven focused
brightness test groups passed, and case 1 was cleanly rebuilt only after
removing exactly `build`.

The frozen production directory is:

```text
validation_runs/i3_case1_rays200_sampling10000_max30000_20261107_12mpi/production
```

It contains private copies of the exact executable and input used:

```text
1c1966cdad3313462fb6d773b7513de3d85c5b4d27463dd6588d0005d8f335b6  amps
063251ec7efb312f82799ccddf328bad8c4d80b4aca809a8f09f22491ee8029d  i3.case1.input
```

The launch command was `mpiexec -n 12 ./amps`, issued from that production
directory at 2026-11-07 21:23:11 CST. Initialization reported all 12 ranks,
distributed 376 blocks as 29--33 blocks/rank, allocated the 20,000,000-slot
particle buffer, resolved the direct/chemical population policies, and
entered iteration 0 with 10,000 iterations to the first output. Output will
be isolated under that directory's `PT/`; the pre-existing repository-root
`PT` tree is untouched. No case 2--6 run has been started.

## I3 case-1 iteration-0 debugger diagnosis (2026-11-07)

The 12-rank production process above was stopped at the user's request; its
partial directory and frozen production executable were retained. There is
no I3 production process running after this investigation. Before rebuilding,
exactly `build` was removed. The case-1 target was
then compiled with the explicit diagnostic flags
`-O0 -g3 -fno-omit-frame-pointer`; `file amps` confirmed that the resulting ELF
contained debug information and was not stripped. A private copy (SHA-256
`d49c9746990638e5e8244026f03ee453d95ee631af8161710cbebc04f7c9e5a5`) was
run under GDB from:

```text
validation_runs/i3_case1_debug_20261107/production
```

The single-rank debug run reproduced the apparent iteration-0 stall. It
initialized the same 376-block/3,008-cell mesh, reported zero active model
particles, and entered the external brightness sampling callback. A GDB
interrupt found the main thread in this call chain:

```text
AMPS::ResonanceFluorescence::Evaluate()
EvaluatePhotonGFactor()
SampleQuadraturePoint()
SampleSegment()
Comet::Brightness::SampleParticleResolvedBrightness()
PIC::Sampling::Sampling()
PIC::TimeStep()
```

This is an indexing defect, not an MPI deadlock and not legitimate work on
active particles. `CometBrightness.cpp` uses
`_getCenterNodeLocalNumber(localCell,0,0)` both to retrieve the center node and
to index `FirstCellParticleTable`. In this DIM1 configuration the center-node
macro expands to `localCell+2` because the center-node array includes two ghost
cells. `FirstCellParticleTable`, however, is a dense eight-entry table for the
eight physical cells and must be indexed by the physical particle-cell number
(`localCell` in this loop).

The decisive GDB values were:

```text
localCell = 6
cellNumber = 8
FirstCellParticleTable[0..7] = {-1,-1,-1,-1,-1,-1,-1,-1}
FirstCellParticleTable[8] = 0                 # out-of-bounds read
PIC::ParticleBuffer::GetAllPartNum() = 0
PIC::ParticleBuffer::GetMaxNPart() = 20000000
particles.size() = 20000000
```

Index 8 reads the adjacent value zero. Zero is the head of the unused
particle-buffer free list, so the brightness gather loop follows that list and
mistakes all 20,000,000 free slots for H particles. It then performs the
Doppler-sensitive g-factor calculation for all of them at each relevant ray
quadrature point. This explains both the CPU activity and the lack of
iteration progress in the production run.

No source fix was made during this diagnostic-only step. The required repair
is to keep the ghost-aware center-node index for `GetCenterNode()` but use a
separate dense physical-cell index for `FirstCellParticleTable`. Regression
coverage must include (1) an empty-buffer sampling call that gathers zero
particles, (2) particles in physical cells 0 and 7, (3) bounds validation of
every particle-list index, and (4) a short one-/multi-rank completed-window
brightness run. Do not restart production until that repair and those tests
are complete. The GDB inferior was killed and GDB exited cleanly; no merge,
commit, branch switch, or master operation was performed.

## I3 case-1 particle-list repair and 12-rank retry (2026-11-07)

*Historical input layout: this retry's executable was compiled from the
retired `i3.case1.input`. The brightness index repair itself is current and
is part of every `i3.conf` build.*

The debugger-confirmed brightness defect is repaired in
`CometBrightness.cpp`. The loop now uses a ghost-aware `centerNodeNumber` only
for `GetCenterNode()` and a separately named, explicitly bounded
`particleCellNumber=localCell` for the DIM1 `FirstCellParticleTable`. Detailed
source comments record why these index spaces must not be combined and how an
out-of-bounds zero aliases the free-particle chain. The configured-source
parser regression asserts the separate physical-cell index and rejects the
former `FirstCellParticleTable[cellNumber]` expression.

The focused brightness suite passed all seven groups after this change. A
clean case-1 rebuild was then made only after deleting exactly
`build`; the generated source was inspected to
confirm the fixed index at the runtime call site. The frozen retry is:

```text
validation_runs/i3_case1_rays200_sampling10000_max30000_indexfix_20261107_12mpi/production
```

Its immutable launch inputs are:

```text
264cf69f8c01abcf8ee250a88b7b755e2b04c2841b7f207ecca0f7d9031d08dc  amps
063251ec7efb312f82799ccddf328bad8c4d80b4aca809a8f09f22491ee8029d  i3.case1.input
```

The active command is `mpiexec -n 12 ./amps`, started at 2026-11-07 21:53:55
CST with terminal output also saved in `case1.run.log`. Initialization entered
iteration 0 at 21:54:00. Unlike the defective run, it immediately advanced:
iteration 10 reported 3,458 real model particles and approximately 5.1 ms per
iteration, and the run passed iteration 700 within the first half minute.
This proves that an empty sampler no longer walks the 20,000,000-slot free
list and that the 12-rank case progresses. It does not yet constitute a
completed-window or scientific validation result. Expected outputs remain
under this retry directory's `PT/plots`; cases 2--6 are not running.

The retry must remain under observation through its first 10,000-iteration
sampling boundary. Verify that both ordinary radial outputs and the H/OH
Tecplot brightness files are written, parse their finite values and metadata,
and compare one-/multi-rank estimators before declaring the runtime sampler
validated. No merge, commit, branch switch, or master operation was performed.

## I3 case-1 first brightness comparison (2026-11-08)

The repaired 12-rank retry crossed iteration 10,000 and wrote its first
complete radial and particle-resolved brightness window.  Both H and OH
brightness tables contain 200 Tecplot points with finite values.  A reusable
comparison utility was added as
`srcHartley2/visualization/plot_i3_brightness_comparison.py`.

The supplied observations exist only as PNG plots inside `Archive.zip`; the
utility therefore does not invent precise numerical points by digitizing the
raster.  Instead it detects the long black plot spines, extracts the interior
observation layer, maps every raster pixel to the panel's documented
logarithmic distance/brightness coordinates, and overlays the AMPS Tecplot
curve.  It removes the archive's red legacy-model pixels by default, including
their anti-aliased fringe; `--keep-archive-model` retains and labels that curve
for a three-way comparison.  The implementation works explicitly in log10
coordinates because placing a linear raster image directly on matplotlib log
axes would distort its pixel-to-data mapping.

The plot automatically expands below the archive's 0.09-kR lower edge when
the calculated profile is fainter, and visibly marks that original edge.  For
case 1 output 0 the chosen range is `1e-4--20 kR`.  This is necessary: over
the archived comparison range the first-window AMPS curves lie below the
observations rather than merely being hidden by the original plot bounds.
This visual discrepancy is diagnostic only; output 0 is not a converged
campaign result.

Generated artifacts are:

```text
validation_runs/i3_case1_rays200_sampling10000_max30000_indexfix_20261107_12mpi/production/PT/plots/i3_case1_H_brightness_comparison.png
validation_runs/i3_case1_rays200_sampling10000_max30000_indexfix_20261107_12mpi/production/PT/plots/i3_case1_OH_brightness_comparison.png
```

Validation exercised the plot-box detector against all 12 supplied H/OH
panels, successfully generated both case-1 PNGs, confirmed their PNG/RGBA
format, and passed `git diff --check`.  `visualization/README.md` documents the
CLI, provenance limitation, examples, lower-range behavior, and optional
legacy-model overlay.  This development changes no solver, core source, input,
or active run state and requires no special master-merge ordering beyond
copying the new script with its README documentation.  No merge, commit,
branch switch, or master operation was performed.

The comparison CLI subsequently made interactive Matplotlib display its
default behavior: it saves the PNG first, opens the live figure, and returns
to the shell when that window closes.  `--no-show` selects the Agg backend for
batch use.  A missing X11/Wayland display also falls back to saved-PNG-only
operation with an explicit diagnostic rather than a backend exception.

The implicit output destination was then aligned with the ordinary AMPS user
workflow: it is `PT/plots/i3_caseN_SPECIES_brightness_comparison.png` relative
to the command's current working directory.  This prevents an input archived
under `validation_runs` from silently placing the comparison away from the
repository-level `PT/plots` a user is inspecting.  Explicit `--output` paths
continue to take precedence for frozen validation artifacts.

The top of `srcHartley2/README.md` now contains a complete operator workflow
for both campaign families: compile-time input selection, exact build-tree
cleanup, per-case I3 regeneration/tests/build/freeze/MPI launch, independent
Halley build/freeze/MPI launch, output-directory semantics, Rubin comparison
with the current 10,000-iteration sampling length, and interactive/headless I3
brightness plotting.  It explicitly records that `Archive.zip` remains an
external regeneration dependency: the plotting script reads observation PNG
members directly without extraction or source embedding, while each generated
comparison PNG is self-contained for later viewing.  This is a documentation
change only and does not alter the active executable or production run.

## I3 numerical MAVEN/IUVS reference extraction (2026-11-08)

The earlier direct-raster comparison described above has been superseded for
ordinary plotting.  All 12 supplied H Ly-alpha and OH(0,0) panels were
digitized into the maintained reference unit:

```text
srcHartley2/visualization/references/i3_brightness_MAVEN/
```

`maven_iuvs_brightness.csv` contains 87 observations: H has respectively
11, 14, 14, 7, 8 and 9 points in cases 1--6; OH has 4, 8, 5, 3, 2 and 2.
Radial centers are exact integer multiples of each case's independently
transcribed `resolution_1e5_km`, and horizontal bounds are one half resolution
on each side, matching the contiguous source bins.  Brightness centers and
visibly separable OH vertical bars are mapped from the reviewed rows in
`digitization_pixels.csv` through each detected plot box on common log10 axes.
Blank vertical bounds mean the raster does not resolve a separate vertical
bar; they do not assert zero physical uncertainty.

`metadata.json` authenticates the source archive, case table, reviewed anchor
table, output CSV and every member PNG, and records dimensions, plot boxes,
units, limits and extraction policy.  The deterministic extractor is
`visualization/extract_i3_brightness_maven.py`.  It validates black observation
ink at every reviewed center and reproduces both the CSV and metadata byte for
byte.  The red legacy-model curve is not extracted.

`plot_i3_brightness_comparison.py` now reads this numerical reference by
default, draws the radial and available vertical bounds as error bars, and no
longer needs `Archive.zip` during ordinary plotting.  `--reference-dir DIR`
selects an alternate directory containing the same
`maven_iuvs_brightness.csv` schema.  The archive is retained only for audit or
regeneration.  Case-1 H/OH PNGs in the repository-root `PT/plots` were
regenerated from the new table and visually inspected.

The focused `test_i3_brightness_reference.py` runs four tests: all 12 subsets
and exact point counts, metadata hashes, alternate-directory loading, and
byte-exact archive regeneration when the external archive is present.  The
top-level and visualization READMEs document extraction limitations, the new
CLI, and detailed clean compile/run/output/comparison workflows for each I3
input and the canonical Halley input.  No solver, core source, simulation
input, active executable, or active production process was changed.

## I3 two-layer input and generic case sections (2026-11-08)

Status: implemented and focused validation passed on branch `1d_comet_dev`.
Nothing was committed or merged, and no master ref was changed. The six
independent per-case input files will not be used again (see the decision at
the top of this file).

The former one-executable-per-observation layout has been replaced by one
compile-time capability file, `input/i3.conf`, and one runtime scientific
receipt, `input/i3.in`. `./Config.pl -application=i3` selects `i3.conf` because
the former tracked `input/i3.input` has been deleted. The optimized executable
then runs any observation with `-i input/i3.in -case N`, where `N` is 1--6.
Case selection is explicit and fail-closed: a file containing cases cannot be
loaded without `-case`, and an unknown name lists the available cases.

### Files changed and why

| File | Role in this redesign |
| --- | --- |
| `Config.pl` | Resolves `input/<application>.input` first and falls back to same-named `.conf`; refreshes and de-duplicates persisted application/input selectors. Historical `.input` precedence was retained. |
| `src/general/input/amps_input.h` | Adds application-independent `ResolveCase`: common values form a base, one direct named case overlays matching assignments/child blocks, nested/unnamed cases and invalid selectors fail. Source locations survive overlay. |
| `src/general/input/test_amps_input.cpp` | Adds INP-05 coverage for overlay replacement/addition, recursive child overlay, required/unknown/spurious selector diagnostics, and nested-case rejection. |
| `src/INPUT_FILE.md`, `src/general/input/README.md` | Make `.input`/`.conf` precedence and named-case behavior part of the normative core input-file contract. These files must be merged to master with the parser, not left in comet-only documentation. |
| `input/hartley2.pl` | Adds `RuntimeInputMode=ON/OFF`, generates the `Comet::RuntimeInput::Required` compatibility flag, and links `CometRuntimeInput.o` in DIM1 when enabled. Legacy inputs default to OFF. |
| `srcHartley2/CometRuntimeInput.h/.cpp` | New strict runtime owner. It parses `application: comet`, invokes the generic case resolver, checks canonical configuration hashes across MPI ranks, rejects missing/unknown keys, and installs sources, heliocentric conditions, mesh targets, simulation controls, population policies, reaction rates, cooling and brightness before allocation. |
| `srcHartley2/RadialMain.cpp` | Adds `-i/--input`, `-case/--case`, mandatory runtime-input enforcement for I3, selected receipt output, and `--validate-runtime-input`. CLI chemistry/random overrides are applied after runtime loading. |
| `srcHartley2/CometChemistry.h/.cpp` | Keeps compiled reaction topology/storage but permits each case's `A` coefficient to be installed by stable reaction ID before the singleton network is constructed. Late/unknown/invalid updates fail. Legacy generated values are unchanged. |
| `srcHartley2/CometBrightness.h/.cpp` | Adds a pre-registration reset and runtime emission-metadata setter so one executable can install case geometry, solar conditions and line metadata before sampling. Mutation after registration fails. |
| `input/i3.conf` | New compile-time capability/layout configuration: seven species, DIM1 spherical backend, block/ghost layout, particle storage/buffer maximum, maximum AMR refinement, IDF layout, collision/reaction topology, and brightness linkage. Its placeholder run values cannot be executed because runtime input is required. |
| `input/i3.in` | New loadable scientific receipt: common physical/numerical settings once plus six `case: N` overlays containing observation provenance, production rates, surface resolution, solar-corrected reaction rates, geometry, radial velocity and Ly-alpha flux. |
| `input/i3.input` | Deleted so `.input` precedence cannot silently select the retired generator-template/one-case executable paradigm. Obsolete `i3.case1.input.Assembled*` artifacts were also removed from the working tree. |
| `srcHartley2/i3/prepare_i3_cases.py` | Generates exactly one `i3.in`; obsolete template substitution and six concrete compile-time renderers were removed. |
| `srcHartley2/i3/tests/test_prepare_i3_cases.py` | Verifies byte-exact combined-file generation, six case blocks, 102 case reaction rates, common settings and compile/runtime separation. |
| `srcHartley2/tests/run_i3_surface_injection_test.sh`, `I3_SURFACE_INJECTION.md` | Runner now supplies explicit runtime input/case (`I3_RUNTIME_INPUT`, `I3_CASE`) so it tests production wiring rather than compiled case constants. |
| `srcHartley2/README.md`, `srcHartley2/i3/README.md` | Document one clean I3 build, runtime case selection, isolated run directories, validator use, and corrected legacy Halley/Hartley `Config.pl` selection. |
| `srcHartley2/MASTER_MERGE.md` | Adds a dependency-ordered two-layer integration wave, compatibility invariants, exclusions, exact tests, and the explicit requirement to merge the generic input description into master. |
| `handoff.md` | Records recovery commands, validation evidence and the final I3-selected workspace state. |

The existing root `Makefile` resonance-fluorescence changes shown by the dirty
tree predate this input redesign. Preserve and review them in their documented
brightness/core-data merge wave; do not attribute the unrelated
`EXTRALINKEROPTIONS_F=-lmpi++` line to this work.

### Validation evidence

- `perl -c Config.pl` and `perl -c input/hartley2.pl`: pass.
- `make -C src/general/input test`: 5/5 groups pass, including INP-05.
- `python3 -m unittest srcHartley2/i3/tests/test_prepare_i3_cases.py`: 6/6 pass.
- `srcHartley2/tests/run_comet_brightness_tests.sh`: all seven focused groups pass, including legacy default-off generation.
- `Config.pl` precedence: `-application=hartley2` selected `hartley2.input`; `-application=i3` selected `i3.conf`; final selectors are `APPLICATION=i3`, `InputFileAMPS=i3.conf`.
- After deleting exactly `build`, `make -j8 EXTRACOMPILEROPTIONS=-O3` completed from `i3.conf`.
- Missing `-i` fails; `-i input/i3.in` without `-case` fails; `-case 9` fails and lists cases 1--6.
- Cases 1--6 each pass `mpiexec -n 1 ./amps -i input/i3.in -case N --validate-runtime-input`; case 1 also passes with two ranks, exercising collective canonical-hash agreement.
- The updated runner passes 200,000 production boundary samples per species for H2O, CO and CO2 through `-i/-case`; all energy/cosine/moment tests pass.
- A separate clean optimized build from `halley1d.input` passes `mpiexec -n 1 ./amps --validate-runtime-input` with `file=<compiled> case=<none>`, proving the Rubin input does not acquire a runtime-input requirement.
- The workspace was then restored with another exact clean build from `i3.conf`; the final two-rank case-1 runtime validation passes.
- `git diff --check`: pass before this status append; rerun before staging.

These checks establish configuration, build, backward-compatibility and
focused physical-source wiring. They are not a converged I3 transport result
and do not establish observational brightness agreement.

## I3 pre-production testing and six-case launch (2026-11-08)

Status: all six I3 production runs submitted to Slurm from the AMPS
repository root on cluster host `voyager1`. Nothing was
committed or merged, and no master ref was changed. The campaign directory
is `validation_runs/i3_campaign_20261108/`; its `scripts/` subdirectory holds
the build, focused-test, pilot, and production batch scripts used below.

### Environment setup on this host

The checkout had no build, `.amps.conf`, or run archive. `Config.pl` cloned
`share/` and `util/`, then
`./Config.pl -install=s -compiler=gfortran,gcc_mpicc` (GCC 11.5, OpenMPI
4.1.7rc1) and `./Config.pl -application=i3` selected `input/i3.conf`. The
login node exposes one core, so every build runs as a Slurm job
(`scripts/build_i3.sbatch`: `rm -rf -- build`, then
`make -j32 EXTRACOMPILEROPTIONS=-O3`). A clean build takes about 40-55 s.
Compute-node `/tmp` is mounted `noexec`, so the focused-test job points
`TMPDIR` at its own test directory. `/usr/bin/time` is not installed on the
compute nodes; peak RSS is read from `/proc/<pid>/status` (VmHWM).

### New code in this round

| File | Change and reason |
| --- | --- |
| `srcHartley2/RadialMain.cpp` | Adds `--validate-brightness [--brightness-out-of-domain exclude\|abort]`. It seeds a deterministic, decomposition-independent H/OH particle set (4,000 log-uniform particles per species, plus one in each of the first and last physical cells), checks an empty sampling window, writes independent column densities from `Math::VisibleShellChordLength()`, times the sampler, and optionally injects one 500 km/s H particle at r=1e7 m. Production paths are unchanged. |
| `srcHartley2/tests/test_i3_brightness_validation.py` | Checker for that mode. `equivalence` compares 1/4/12-rank tables, the zero empty window, and the independent columns; `exclude` checks EXCLUDE_AND_REPORT accounting against a base run. |
| `srcHartley2/CometBrightness.cpp` | Sampler load accounting. Under particle-number balancing (forced by `ForceRepeatableSimulationPath=on`), each owner block now adds 0.04 particle-move equivalents per g-factor evaluation to `ParallelLoadMeasure`; under execution-time balancing it adds its measured time instead. This changes MPI block ownership only; no sampled value changes. |
| `srcHartley2/README.md` | Documents the brightness fixture and checker commands. |

### Focused tests on the production executable

The final executable (SHA-256 `c33c69daf31dbda56d3967dd4df4f691035469146c0201cee332145e0335cbea`)
passed 33/33 tests in `validation_runs/i3_campaign_20261108/tests_122856/`:

- byte-exact `input/i3.in` regeneration, generator unit tests (6), the
  brightness group (7), the generic input parser (INP-01..05), and the
  MAVEN/IUVS reference tests;
- `--validate-runtime-input` for cases 1-6 (one rank) and case 1 (two ranks);
  missing `-i`, missing `-case`, and `-case 9` fail with their diagnostics;
- surface injection at 175 K (200,000 samples per parent), chemistry
  ensembles (cases 1 and 5, 1 and 4 ranks; 8 neutral channels plus 9 sinks),
  local-timestep face transport (4 ranks, both movers, 95-face path), and
  weight adaptation (4 ranks);
- the new brightness fixture for cases 1 and 5 at 1, 4, and 12 ranks.
  Maximum cross-rank relative difference 1.7e-14, maximum difference from
  the independent column 1.4e-14, 200/200 positive rays per species, and a
  zero empty window. EXCLUDE_AND_REPORT recorded 1,470 excluded evaluations,
  only on rays crossing the injected shell. ABORT terminated with
  `brightness g-factor evaluation failed ... velocity=-4.24e5 m/s`.

### Performance pilots and the load-balance defect

Pilots ran the production executable and `i3.in` with only the sampling
length and iteration cap changed (`scripts/pilot.sbatch`). Peak RSS is about
1.4 GB per rank, roughly 180 GB per 128-rank node.

The first 12,000-iteration pilots showed time per iteration rising much
faster than population (0.37 s/iteration near iteration 5,500 with about
1.5 M particles), and 128 ranks were no faster than 64. GDB stack sampling
of all 64 ranks found 63 waiting in the chemistry `MPI_Allreduce`
(`FinishBatchedPass`) while one rank evaluated g factors in
`SampleParticleResolvedBrightness()`. The AMPS "Photo Chemistry Time"
column was absorbing that wait. Particle-number balancing ignored the
sampler, whose cost scales with emitting particles times crossing rays, so
the rank owning the outer H/OH front did about 86% of the 200-ray work.

After the load-accounting change, the same pilots ran at 0.058-0.086
s/iteration over iterations 5,000-6,000 (4-6x faster), and 128 ranks beat
64. A repeated stack sample showed 37 of 128 ranks working in the sampler at
once. Later timings for the 128-rank case-5 pilot:

| Iterations | s/iteration | Particles at end |
| --- | ---: | ---: |
| 8,000-9,000 | 0.233 | 2.42 M |
| 9,000-10,000 | 0.315 | 2.69 M |
| 10,000-10,705 | 0.457 | 3.36 M |

The iteration-10,000 weight adaptation in the pilots was split-only
(2.72 M -> 3.24 M and 2.86 M -> 3.32 M), with at most 18-20 products per
reaction against a limit of 20 and no emergency coarsening. The sampler
remains the dominant cost.

### Production launch

Each case runs the unmodified `input/i3.in` with 128 MPI ranks on one
exclusive 128-core node, in the `superlong128cores` partition with a 4-day
limit (the runs cannot be restarted). Each run directory freezes `amps`,
`i3.in`, `i3.conf`, and `picGlobal.dfn` with `SHA256SUMS` and
`provenance.txt`, and writes per-rank peak RSS to `rss.*.txt`.

| Case | Slurm job | Run directory |
| --- | --- | --- |
| 1 | 122873 | `validation_runs/i3_campaign_20261108/case1_np128/production` |
| 2 | 122874 | `validation_runs/i3_campaign_20261108/case2_np128/production` |
| 3 | 122875 | `validation_runs/i3_campaign_20261108/case3_np128/production` |
| 4 | 122876 | `validation_runs/i3_campaign_20261108/case4_np128/production` |
| 5 | 122877 | `validation_runs/i3_campaign_20261108/case5_np128/production` |
| 6 | 122878 | `validation_runs/i3_campaign_20261108/case6_np128/production` |

All six started on nodes node032 and node009-node013. Each log reports
`comet runtime input=i3.in case=N` and the tested executable hash. Output
goes to each run directory's `PT/plots`; brightness tables follow
`pic.H.brightness.s=0.out=K.dat` and `pic.OH.brightness.s=3.out=K.dat` after
each 10,000-iteration window. The runs stop at convergence or at 30,000
iterations. From the pilots, about 6-17 h per case is expected if the late
cost levels off at 1-3 s/iteration; this is an estimate, not a measurement.

The post-run items listed at launch (convergence status, H tau review,
nonconverged comparison plots, and domain-truncation disclosure) are
reported in the next section. The H g-factor evaluation cost could also be
reduced (for example with a tabulated g(v)), but that would change numerical
values and needs its own validation.

## I3 six-case production results (2026-11-08)

All six jobs (122873--122878) exited with code 0. No log contains an error,
abort, or velocity-domain exclusion. Each run reached the 30,000-iteration
cap on 128 ranks. Peak RSS was 1.42 GB per rank, about 180 GB per node. The
late cost was 0.42--0.60 s/iteration with 3.9--4.3 M particles, well below
the 1--3 s/iteration pilot extrapolation.

### Convergence

No case converged. The criterion is the largest volume-weighted relative
L1 change between successive 10,000-iteration windows, below 0.05 for three
consecutive windows. The first window is only the reference, so it needs at
least four windows; 30,000 iterations provide two comparisons.

| Case | Elapsed | Window 2 (iteration 20,000) | Window 3 (iteration 30,000) | Limiting species, window 3 |
| --- | --- | ---: | ---: | --- |
| 1 | 2:40:17 | 0.959 | 0.367 | O |
| 2 | 3:18:25 | 0.953 | 0.350 | O |
| 3 | 3:20:24 | 0.924 | 0.298 | H |
| 4 | 3:37:29 | 0.903 | 0.288 | H |
| 5 | 3:36:26 | 0.905 | 0.276 | H |
| 6 | 3:36:28 | 0.906 | 0.280 | H |

Every brightness result below is therefore nonconverged.

### Brightness against MAVEN/IUVS

The table gives the median model/observed ratio over the observation points
inside the 1e9-m domain. AMPS brightness is log-interpolated to each point.

| Case | Species | Points in domain | Output 0 | Output 1 | Output 2 | Output-2 range |
| --- | --- | ---: | ---: | ---: | ---: | --- |
| 1 | H | 11/11 | 0.003 | 0.072 | 0.110 | 0.095--0.140 |
| 1 | OH | 4/4 | 0.026 | 0.226 | 0.250 | 0.207--0.325 |
| 2 | H | 14/14 | 0.005 | 0.083 | 0.120 | 0.107--0.160 |
| 2 | OH | 8/8 | 0.025 | 0.199 | 0.216 | 0.189--0.221 |
| 3 | H | 14/14 | 0.011 | 0.152 | 0.221 | 0.190--0.240 |
| 3 | OH | 5/5 | 0.057 | 0.313 | 0.332 | 0.321--0.356 |
| 4 | H | 6/7 | 0.019 | 0.236 | 0.339 | 0.078--0.636 |
| 4 | OH | 3/3 | 0.097 | 0.480 | 0.532 | 0.464--0.737 |
| 5 | H | 5/8 | 0.022 | 0.277 | 0.389 | 0.162--0.721 |
| 5 | OH | 2/2 | 0.159 | 0.737 | 0.792 | 0.720--0.864 |
| 6 | H | 4/9 | 0.026 | 0.315 | 0.442 | 0.232--0.702 |
| 6 | OH | 2/2 | 0.125 | 0.610 | 0.654 | 0.577--0.731 |

- OH grew 6--11% between outputs 1 and 2, so it is close to steady. It is
  still 1.3--4.6 times below the observations.
- H grew 1.4--1.5 times per window and is still rising. It is 2.3--9 times
  below the observations.
- The model remains below every observation point.

### H line-center optical depth

The sampler flags rays with H tau above 0.1 but never attenuates them.

| Case | Flagged rays | Maximum tau | tau > 0.1 out to | tau > 1 out to | Maximum tau at an observation point |
| --- | ---: | ---: | ---: | ---: | ---: |
| 1 | 58/200 | 0.22 | 101 km | none | 0.017 |
| 2 | 95/200 | 0.43 | 1,110 km | none | 0.036 |
| 3 | 153/200 | 1.72 | 47,500 km | 88 km | 0.096 |
| 4 | 182/200 | 12.1 | 311,000 km | 12,200 km | 0.206 |
| 5 | 182/200 | 12.7 | 311,000 km | 13,900 km | 0.196 |
| 6 | 181/200 | 11.0 | 292,000 km | 11,400 km | 0.156 |

The optically thick region lies well inside the innermost observation point
in every case. Only the innermost H point in each of cases 4--6 exceeds tau
0.1, with tau at most 0.21. Attenuating these points would lower them by at
most about 20% and would not change the comparison.

### Domain truncation

The H column is accumulated only inside r = 1e9 m (10^6 km), and the modeled
H brightness falls to zero at that impact parameter. H points beyond the
boundary are excluded: 1 of 7 in case 4, 3 of 8 in case 5, and 5 of 9 in
case 6. In-domain H points near the boundary are also biased low by the
missing outer column; for example, the case-4 minimum ratio of 0.078 is at
its outermost in-domain point. All OH points lie inside the domain.

### Output files

Plots are in `validation_runs/i3_campaign_20261108/plots_out0`, `plots_out1`,
and `plots_out2`, as `i3_caseN_{H,OH}_outK.png` plus six-case montages. The
ratio table is `validation_runs/i3_campaign_20261108/plots_out2/ratios_out0_out1_out2.txt`.

### Next step

Converged results need longer runs. At the late cost of about 0.5
s/iteration, each 10,000-iteration window takes about 1.4 h per case. The
runs cannot be restarted, so a longer run starts again from iteration 0 with
a larger `simulation.max_iterations` in `input/i3.in`; no rebuild is needed.

## I3 campaign 2 launch (2026-11-08)

Six new production runs started on 2026-11-08 at about 23:35 CST (Slurm jobs
123164–123169, one 128-rank node each) in
`validation_runs/i3_campaign2_20261108/caseN_np128/production`. The
executable is the campaign-1 binary (commit `d61307f`, SHA-256 `c33c69da...`).
Only `input/i3.in` changed; it was regenerated with these settings:

- Production rates: H2O, CO, and CO2 multiplied by 3.95, 4.81, 2.97, 1.74,
  1.25, and 1.53 for cases 1–6, from the OH match in
  `validation_runs/i3_campaign_20261108/analysis_q_scaling/README.md`. The
  factors live in the new `srcHartley2/i3/reference/i3_rate_scaling.csv`, so
  the archive transcription in `i3_cases.csv` is unchanged.
- Domain: `mesh.outer_radius_m` raised from 1e9 to 3e9 m, and the brightness
  impact-parameter grid extended to 3e9 m. The 30-level AMR limit makes the
  finest cell 0.349 m, so the generator floors the surface request at
  0.353 m. Cases 2–6 are therefore 1.3–4.6 times coarser than the surface
  mean free path (chosen over 1e10 m, which would have been 2.6–15 times
  coarser, and over keeping 1e9 m).
- Time step: Courant number 0.8 on conservative speed bounds, with the
  parent bound raised from 400 to 2,500 m/s. The H time step is 3.1 times
  longer than in campaign 1 and the parent time step 0.51 times as long.
- Iteration cap: 150,000.

Before launch, `smoke_123163` passed 25/25 checks (see the campaign README).
The generator changes (`prepare_i3_cases.py`, its unit tests, and
`i3_rate_scaling.csv`), the regenerated `input/i3.in`, and the
`srcHartley2/i3/README.md` update are committed together with the interim
results below.

## 2026-11-09: portable single-entry generic-comet test suite

Status: implemented and exercised on branch `1d_comet_dev`; no commit, branch
switch, merge, or master operation was performed. The authoritative operator
entry point is now:

```sh
srcHartley2/tests/run_all_stage_tests.sh
```

The default command no longer depends on the repository-root
`validation_runs/`, a configured/generated `build/`, root-level assembled
input files, or an AMPS science executable. It verifies and consumes the
source-controlled 8.9-MiB fixture bundle under `srcHartley2/tests/data/`,
builds only focused temporary test programs, and runs the Stage-0 diagnostic
at 1, 4, and 12 MPI ranks. `--evidence-root` and the existing Stage-4--9
`fresh` modes remain available for auditing full external archives or a
caller-supplied already-built executable.

### Modified and added files

| File or directory | Change and reason | Later master-integration note |
| --- | --- | --- |
| `srcHartley2/tests/run_all_stage_tests.sh` | Makes the compact tracked evidence tree the default; authenticates it before parsing; adds the generic runtime-input grammar, I3 six-case generator, MAVEN/IUVS regeneration, and brightness/resonance-fluorescence groups to the Stage-0--10 stack. The runner retains strict preflight, per-test logs, aggregate test/group counts, failed-log reporting, and optional durable JSON output. | Integrate after all child runners and `tests/data`. Preserve the distinction between portable protocol replay and fresh science validation. Do not restore a generated-build dependency. |
| `srcHartley2/tests/data/` | New 8.9-MiB source-controlled fixture bundle. It contains only checker-consumed logs, diagnostic sidecars, compact numerical tables, generated-policy provenance, a two-species/two-cell Stage-9 checkpoint pair, and the reviewed I3 source archive. `SHA256SUMS` covers 179 files, including nested provenance manifests. | Add the complete directory, including files named `run.log`, CSV/JSON tables, nested manifests, and `i3/Archive.zip`. Do not substitute the multi-gigabyte production archives or omit files because their names resemble generated output. Review any future numerical fixture change scientifically and regenerate the root manifest. |
| `srcHartley2/tests/data/README.md` | Documents provenance, reduction boundaries, the difference between regression fixtures and science evidence, and the exact deterministic manifest-regeneration command. | Keep beside the data. A compact fixture PASS cannot be cited as a new converged Rubin or I3 run. |
| `srcHartley2/tests/README.md` | Defines the one public command, portable requirements, log retention, external-archive audit mode, and scientific limitations. | Make this the test operator entry point; individual child scripts remain implementation components. |
| `srcHartley2/tests/test_stage3_parser.py`, `run_stage3_parser_tests.sh` | Replaces the prior copy of a developer's generated `build/` and root assembled files with source-only temporary parser fixtures. Covers valid, omitted-default, duplicate, invalid, missing-seed, and chemistry-off behavior. | Integrate the Python test with the small shell wrapper. It imports the maintained source-only `ParserFixture` from the Stage-9 parser tests, so merge those together. |
| `srcHartley2/tests/test_stage0_minor_species.py` | Accepts an explicit `HARTLEY2_STAGE0_BASELINE_MANIFEST` for the compact baseline while preserving the full immutable science manifest as the direct-invocation default. | Preserve both paths; the portable manifest must not silently replace the scientific provenance contract for external archives. |
| `srcHartley2/tests/test_i3_brightness_reference.py` | Makes archive-to-table regeneration mandatory using `tests/data/i3/Archive.zip`; a missing source artifact is now a failure rather than a skipped test. | Integrate with the bundled archive and updated reference metadata in the same change. |
| `srcHartley2/visualization/references/i3_brightness_MAVEN/metadata.json` and related README files | Records the canonical tracked archive path and updates regeneration commands. | Do not stage the unrelated root `Archive.zip`; do stage the checksummed test artifact. |
| `srcHartley2/README.md`, `srcHartley2/i3/README.md`, `srcHartley2/visualization/README.md`, `srcHartley2/CHEMISTRY.md` | Documents the single runner, portable data, archive location, dependencies, scope, and full-archive/fresh modes. | Reconcile these active instructions after code and test paths; older dated validation sections remain historical evidence. |
| `srcHartley2/MASTER_MERGE.md` | Adds the portable fixture/test wave and explicitly distinguishes the tracked I3 archive from the unrelated root copy and other excluded run artifacts. | Follow the documented order; no text in that file authorizes a merge. |

### Validation result

The final command was:

```sh
srcHartley2/tests/run_all_stage_tests.sh \
  --log-root /tmp/hartley2-portable-suite-final-20261109-r4
```

It passed all 44 independently logged jobs in all 14 groups, with zero skips
and zero failures. Coverage included the data manifest, generic case-overlay
grammar, I3 generation, byte-exact MAVEN/IUVS extraction, focused brightness,
Stage-0 MPI at 1/4/12 ranks, and all Stage-1--10 checks. The durable files are:

```text
/tmp/hartley2-portable-suite-final-20261109-r4/summary.txt
/tmp/hartley2-portable-suite-final-20261109-r4/all_stage_result.json
```

The receipt records committed HEAD
`d61307f929ab54dbc2f517be372e28d50311999d`; this implementation was still an
uncommitted working-tree change during the run. After the files are reviewed
and committed, rerun the same command with a new log directory to bind release
evidence to the resulting commit. The current portable PASS validates the
implementation and fixture contracts; it is not a new converged comet run.

## I3 campaign 2 interim results (2026-11-09 11:10 CST)

Cases 1, 2, and 5 converged at iteration 70,000. Their windows 5–7 passed the
5% criterion, and the runs ended normally after 6.8, 9.8, and 9.4 h. Cases 4
and 6 have passed two of the three windows and should end at iteration
80,000. Case 3 had passed none at iteration 45,600. No log contains an error.
Details, plots, and scripts are in
`validation_runs/i3_campaign2_20261108/README.md`.

### Brightness against MAVEN/IUVS

The table gives the median model/observed ratio at the latest output.
"Edge-corrected" multiplies the AMPS H by the ratio of the untruncated to the
3e9-m-truncated column from the collisionless test-particle reference.

| Case | Output | State | OH | H | H edge-corrected | H AMPS/collisionless |
| ---: | ---: | --- | ---: | ---: | ---: | ---: |
| 1 | 6 | converged | 0.95 | 0.55 | 0.60 | 1.01 |
| 2 | 6 | converged | 0.99 | 0.72 | 0.75 | 1.02 |
| 3 | 3 | running | 0.97 | 0.94 | 1.02 | not converged |
| 4 | 6 | running | 1.04 | 1.00 | 1.14 | 1.04 |
| 5 | 6 | converged | 1.00 | 0.73 | 0.89 | 1.07 |
| 6 | 6 | running | 1.03 | 0.89 | 1.16 | 1.02 |

- OH: the rate scaling works. OH matches within 5% in every case.
- Collisional H enhancement: 1–7%. This is smaller than the campaign-1
  estimate of 3–16%, so collisions with the added CO/CO2 do not explain the H
  deficit.
- H with OH matched:
  - too low in cases 1 and 2, by factors of 1.7 and 1.3, with a flatter
    radial slope than observed;
  - about 10% low in case 5;
  - 14–16% high in cases 4 and 6.
  A single H2O rate cannot fit both species in cases 1 and 2 with the
  current product energies and OH lifetime.
- Domain edge: the 3e9 m edge still removes 13–21% of the H column at the
  observation points in cases 4–6, and 4–5% in cases 1–2. The collisionless
  reference supplies the correction, so the domain need not grow further.

### Performance regression

The cost per iteration rose from about 0.4 s to 1.6–3 s per iteration in all
running cases, with no growth in particle count (3.2–3.7 M). Cases 2, 4, and
6 recovered briefly to about 0.7 s at about 08:10 CST. The recovery
coincided with neither an output nor a common iteration. The per-rank timing
table puts the time in the photochemistry phase and in particle exchange.
The time is equal on every rank, including ranks without particles, so most
of it is waiting at a synchronization point. The nodes are exclusive. The
cause has not been found yet.
