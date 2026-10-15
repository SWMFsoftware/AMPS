# Enceladus plume, plasma, and charged-dust model: development and validation plan

Revision: 2026-10-13. Status: proposed implementation roadmap, grounded in a read-only source inspection and the local Cassini inventory. No application was built or executed during preparation of this plan; no development gate below is yet qualified.

## 1. Objective and ordered milestones

Develop a reproducible Cassini-constrained model in four successive stages. Preserve the fitted neutral source and the observation operators across stages so that changes in plasma or dust physics can be assessed without silently changing the underlying plume.

| Stage | AMPS responsibility | BATS-R-US responsibility | Exit deliverable |
|---|---|---|---|
| 0 — neutral-source inference | Fit the multi-plume neutral model to INMS; verify kinetic transport and gravity | None | Versioned plume parameter ensemble, qualified neutral solver, independent validation results |
| 1 — coupled plasma interaction | Kinetic neutrals and pickup ions; reaction and exchange accounting | Ambient plasma and self-consistent electromagnetic/MHD response | Conservatively coupled AMPS–BATS-R-US flyby simulations |
| 2 — plasma inside AMPS | Neutrals, pickup ions, ambient plasma, and self-consistent fields | Reference/comparison only | Qualified AMPS-only plasma model and common-physics comparison to Stage 1 |
| 3 — charged dust | Stage 2 plus dust injection, motion, charging, collection, and feedback where required | None required | Joint gas–plasma–dust validation and quantified limits of applicability |

Sequence: data qualification → neutral fit → neutral kinetic verification → prescribed-field pickup tests → two-way plasma coupling → AMPS plasma benchmarks → AMPS-only flybys → passive dust → charging → coupled dusty plasma.

Development of isolated later-stage components can proceed once their interfaces are defined, but their science milestones remain open until predecessor gates pass. A compilation, initialization, or short smoke run cannot substitute for an evolving, converged physical validation case.

## 2. Evidence from the current implementation

Source root: `/home/vtenishe/Mars2/AMPS`. The observations below describe inspected source, not proven runtime behavior. Configuration-generated code must be inspected before assigning production capabilities.

| Evidence | Current implementation or concern | Required action |
|---|---|---|
| `srcEnceladus--multi-plume/EnceladusMultiPlume.h`; `input/enceladus--multi-plume.pl` | Source templates contain one individual plume and one stripe. The parser rewrites counts and source tables in generated files. | Preserve that workflow initially; test that the requested eight-vent configuration and stripe alternatives are actually generated. Do not label the template counts a physical limitation. |
| `EnceladusMultiPlume_SourceModel.cpp` | Aggregates per-species source rates and chooses point versus stripe sources by their rates. The selection ratio has no local zero-total guard. | Add validation of rates and explicit zero-source behavior; test inactive and empty source sets. |
| `EnceladusMultiPlume_SourceIndividualPlumes.cpp` | Rate-weighted source selection, Maxwellian injection, and local MPI ownership checks exist. Positions and velocities are copied between IAU and SO buffers. | Verify normalization, outward support, MPI weights, and frame assumptions. Copies are valid only when the frames coincide. |
| `EnceladusMultiPlume_SourceLocation.cpp` | Converts west longitude to east longitude, calls `srfrec_c`, and builds source direction from local north/azimuth/tilt. | Check signed southern latitude, planetographic versus planetocentric latitude, surface normal, polar singularity, and km-to-m conversion. |
| `EnceladusMultiPlume.cpp` | Both named frames are `IAU_ENCELADUS`; the application acceleration routine contains central gravity. `GetTAA` returns zero. | Define a consistent dynamics frame and prove which acceleration callback is compiled. A rotating frame requires the appropriate inertial terms or a justified short-time approximation. |
| `main.cpp` in the multi-plume application | Non-dust local timestep branch recognizes H₂O only; a file-background path and fixed vector rotation are present. | Add species-aware charged-particle timesteps, synchronized physical time, validated field transforms, and an explicit background-provider interface. |
| `EnceladusMultiPlume_ColumnIntegral.cpp` | Column-density variable list returns zero and its integrand is empty. | Implement and independently test line-of-sight sampling before claiming UVIS reproduction. |
| `EnceladusMultiPlume.cpp::GenerateInitialGrainParameters` | Returns `false`; the former body is commented out. | Stage 3 must restore/replace and verify grain injection. Dust declarations do not demonstrate an operational dust source. |
| `srcEnceladus/Enceladus.h`, `plume.h`, `plume.cpp`, `Enceladus.cpp` | Legacy eight-plume analytical machinery and a grain source exist. Fitting helpers use bounded temperature/speed ranges and mutable distribution-function state. | Reuse only after normalization, reentrancy, units, and parameter provenance tests. Treat legacy bounds as assumptions, not measured uncertainties. |
| `input/enceladus.input`; `input/enceladus--multi-plume.input` | Legacy deck is dust-focused; multi-plume deck includes H₂O and dust, enables file coupling, and contains old kernels/epochs. Dust acceleration is selected in the deck. | Create explicit stage-specific inputs; inspect generated mover dispatch so neutral particles cannot accidentally enter dust-only callbacks. |
| `src/models/dust/Dust.cpp`, `Dust_charging.cpp` | Shared dust/charging infrastructure exists; an inspected steady-state charging branch uses proton mass in its ion thermal-current expression. | Audit branch selection, water-group ion masses, drift, electron populations, current signs, and charging validity before reuse. |
| `src/pic/pic_field_solver_ecsim.cpp` | An electromagnetic ECSIM solver exists in the shared core. | Evaluate reuse in Stage 2; its presence does not establish an Enceladus configuration, physical fidelity, or feasible production cost. |
| `Makefile.test`, `test_Enceladus_compile` | Legacy test target exists and configures `nospice`; its run-directory recipe removes/recreates a directory. | Use it only as a reviewed regression reference. Add a non-destructive, SPICE-enabled science test workflow with unique output directories. |

The inspection did not establish a working Enceladus-specific pickup-ion chemistry network, conservative two-way BATS-R-US coupling, or self-consistent AMPS plasma configuration. These are development tasks, not assumed existing capabilities.

Repository rules in `/home/vtenishe/Mars2/AMPS/AGENTS.md` concern an active CME campaign and preservation of other applications. This document plans Enceladus work; it does not execute that campaign or authorize unrelated changes. During future implementation, preserve dirty-tree edits and site configuration, use the repository's required native build preflight, regenerate the selected application, and build through the established top-level workflow. Do not replace the current CME configuration in a shared running build. Record the exact configuration/build commands after confirming the application's supported selector. No new command proposed below is claimed to exist today.

## 3. Data qualification and observation contract

### 3.1 What is actually available locally

The acquisition inventory reports 383 downloaded files, approximately 5.60 GB, plus eight generated SPICE metakernels. Download verification establishes reproducible received bytes; it does not establish calibration accuracy or usable science coverage. The measurement-volume estimate is provisional until record selection and response processing are complete.

| Product | Use | Qualification needed before comparison |
|---|---|---|
| INMS L1A V2 counts for E3, E5, E7, E13, E14, E17, E18, E21 | Neutral source fitting; selected ion measurements when modes permit | Mode/mass-channel selection, exposure, sensitivity, background, saturation/dead time, ram direction, spacecraft velocity, chamber memory/sticking and possible grain contamination |
| CAPS calibrated sensor products and derived moments/potential | Pickup spectra, plasma flow, composition, temperatures, electron constraints | Actual measurement coverage, angular sampling, spacecraft potential, energy/mass response, moment validity in anisotropic distributions |
| MAG KRTP vectors | Magnetic perturbations, draping, interaction geometry | Quality flags, upstream/background estimation, correct spherical-basis conversion and vector-frame transformation |
| RPWS upper-hybrid density products | Electron-density profiles and depletion | Valid resonance detections, uncertainty and method assumptions; missing detections are missing data |
| CDA volumes and HRD processed bins | Dust encounter profiles, impact rates and selected size/speed constraints | Record-level coverage, live time, attitude, detection efficiency, calibration, dead time and background subtraction |
| SPICE | Trajectory, orientation, instrument look directions, source geometry | Actual kernel coverage and precedence, state/frame evaluation throughout selected intervals, reconstructed versus gap-filled attitude |

The manifest/inventory contains CAPS derived products for six encounters: E3, E5, E7, E13, E17, E18. **No E14 CAPS products are recorded.** The earlier `metadata/coverage_notes.md` assertion that E14 derived products were included is inconsistent with the inventory and must be corrected by the qualification task. E21 has no CAPS products. Partial CAPS blocks at other encounters require a coverage mask, not an assumption of full ±2 h observations.

UVIS E2 and the exact historical density/sticking-corrected data underlying Tenishev et al. (2010) are not established as available in the local package. The inspected INMS label references `L1A_STRUCT_06.FMT`, which was not found among the downloaded files. Acquire required formats, calibration/response documentation, and historical comparison products before claiming a complete reduction or exact paper reproduction.

The E3 CDA tarball disagrees with a published MD5 file that predates it. Repacking is a plausible explanation, not verified provenance. All eight gzip checks passed, but gzip integrity and locally calculated SHA-256 do not resolve publisher authenticity. Check extracted members against published member checksums or an authoritative corrected archive; mark affected products `provenance_pending` until resolved.

### 3.2 Required processing outputs

Create immutable raw inputs and versioned derived products under the data workspace:

```text
data/processed/<release>/<instrument>/<flyby>/
metadata/<release>/coverage.csv
metadata/<release>/quality_rules.yaml
metadata/<release>/calibration_manifest.json
metadata/<release>/validation_split.yaml
validation/<campaign>/<case_id>/
```

Each accepted record must carry UTC and ET, time bounds/exposure, spacecraft position and velocity, named frames and units, measurement channel/mode, observable and uncertainty, quality/exclusion reason, response/calibration version, and raw-file/checksum lineage. Preserve native sampling and missing-data masks. Store counts separately from any inferred density or flux. Products derived from the same underlying measurements must not be treated as independent likelihood terms.

Default context is closest approach ±2 h; determine exact closest approach from the adopted SPICE solution rather than minute-rounded metadata. Define narrower plume/core intervals and upstream/background intervals separately, without selecting them after seeing model residuals. Expand windows only when settling, upstream identification, or background estimation requires it.

Use geometric spacecraft states for local in-situ comparisons and document aberration conventions separately for remote sensing. Transform position/velocity with a full state transformation when frames rotate; vector rotation alone is insufficient for velocity. NAIF documents this distinction in its [state transformation interface](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/FORTRAN/spicelib/sxform.html).

### 3.3 Freeze the inference/validation split

1. Historical reproduction: E3/E5 INMS, and E2 UVIS only after acquisition, reproduce the paper's joint-fit setup. These are calibration observations in that experiment.
2. Modern Stage 0 calibration: E3/E5, with explicit calibration nuisance parameters and instrument memory treatment. Use whole-flyby cross-validation to assess transferability.
3. Development validation: E7/E13, used to choose physically motivated extensions. Once used for model selection, these cease to be blind tests.
4. Final transfer tests: E14/E17/E18/E21 where qualified data exist. Freeze predictions before opening scored intervals. E14/E21 cannot validate CAPS spectra using this package.
5. Separate campaigns can use different splits if registered before fitting. Every data release lists which observations fitted plume sources, upstream plasma, dust parameters, and instrument nuisance parameters.

A time-variable source is permissible, but a fit with newly inferred per-flyby source strengths is a conditional reconstruction, not a blind prediction of that flyby. Report both experiments separately. Identify upstream boundary data independently of the plume/wake intervals scored as validation.

## 4. Stage 0 — fit and qualify the neutral multi-plume model

### 4.1 Scientific baseline from Tenishev et al. (2010)

Reference: V. Tenishev, M. R. Combi, B. D. Teolis, and J. H. Waite, *An approach to numerical simulation of the gas distribution in the atmosphere of Enceladus*, JGR 115, A09302, DOI [10.1029/2009JA015223](https://doi.org/10.1029/2009JA015223). Inspected local copy: `/home/vtenishe/Downloads/Tenishev-2010-JGR.pdf`.

Reproduce the eight point sources plus spherical/background component before extending to distributed tiger stripes. The paper uses a convected Maxwellian source, a gravity-free semi-analytical calculation, and a gravitational test-particle comparison. Its two fits use original versus sticking-processed INMS E3/E5 data, together with UVIS E2. The historical INMS data already include an upward factor-of-three correction relative to the initial report; never apply that correction automatically to modern archive reductions.

Historical checks, not new hard priors: Table 2 gives total plume rates of 4.2×10²⁸ and 2.6×10²⁸ molecules/s, and spherical rates of 3.6×10²⁵ and 4.7×10²⁵ molecules/s for Cases 1 and 2. Only plumes I, II, III, VII are strongly constrained in that geometry; other sources have limits rather than uniquely recovered parameters. Table 1 uses planetographic latitude and west longitude. Preserve the southern hemisphere sign when transferring values to code.

Use the paper's opening-angle definition, α = atan(sqrt(3 kB T/m)/u), when comparing tables; this convention is not a fitted Gaussian angular standard deviation. Its empirical background n(r)=c₀+c₁/r² must be distinguished from an explicitly injected spherical source: c₀ represents ambient/background density and is not a surface production rate.

### 4.2 Implementation work packages

| ID | Work and concrete output | Dependencies |
|---|---|---|
| S0.1 | Audit source/configuration paths; capture revision, dirty diff, toolchain, generated species/source/mover configuration; create a neutral-only input with dust and plasma coupling disabled | Source inspection |
| S0.2 | Qualify INMS/geometry; acquire missing format and calibration resources; build counts and density comparison paths with uncertainty and memory treatment; create coverage and split manifests | Section 3 |
| S0.3 | Implement a versioned source schema: source ID, location convention, direction, Q in molecules/s, T in K, u in m/s, species fractions, activity, time law and prior bounds; export consistently to analytic and particle models | S0.1 |
| S0.4 | Independently verify analytical source normalization, angular support and inverse-square behavior; implement spherical/background component and finite exposure/trajectory sampling | S0.3 |
| S0.5 | Reconstruct historical cases and implement a deterministic fitting driver, bounded optimization, uncertainty estimation and source-identifiability analysis | S0.2–S0.4 |
| S0.6 | Transfer source ensembles to AMPS; qualify gravity, boundary losses/reimpact, MPI source ownership, steady/transient sampling, mesh and timestep convergence | S0.5 |
| S0.7 | Evaluate collisionless applicability through local mean free path/Knudsen number and flight/loss times; add near-vent DSMC or a matched source surface if required for scored intervals | S0.6 |
| S0.8 | Compare point-source and stripe alternatives with held-out data; freeze the simplest supported neutral model, source posterior and limitation report | S0.6–S0.7 |

Fit nonnegative source strengths first at fixed geometry/velocity shape. Then release speed, temperature, tilt and azimuth in bounded groups supported by the data. Explicitly test Q/u and temperature/opening-angle degeneracies. Do not fit every source property independently for every flyby without evidence of identifiability. Use profile likelihoods or posterior ensembles, multi-start optimization, and leave-one-flyby-out tests. Preserve correlations between parameters in downstream ensembles.

For counts, predict λ = background + exposure-integrated instrument response applied to the modeled velocity distribution; use an appropriate count likelihood. For published density products, use their response assumptions and covariance. Include uncertainty from sensitivity, geometry, sticking/memory and Monte Carlo noise. A flexible instrument nuisance model must not absorb arbitrary plume structure. For memory, predict both chamber accumulation and release or adopt a documented reduction; avoid applying both corrections to the same signal.

Do not use particle noise as the optimizer's principal gradient signal. Fit the verified analytical surrogate first, test its error against AMPS, and use sufficiently converged AMPS refinements or a separately validated emulator when gravity/collisions affect the inference. Full UVIS reproduction remains open until the observation data, ray geometry, and integrator exist; the INMS-only milestone can be completed and labeled separately.

### 4.3 Stage 0 tests and acceptance

| Test | Independent reference / acceptance |
|---|---|
| N01 source geometry | Hand-constructed cardinal locations/directions and independent SPICE transforms; correct south sign and west/east convention; no polar NaNs; reconstructed direction error <10⁻⁸ rad in non-singular synthetic cases |
| N02 normalization and selection | Numerical flux integral equals prescribed Q within 10⁻⁴ relative error after quadrature convergence; sampled source fractions and velocity moments agree within 3 standard errors |
| N03 zero/inactive sources | Zero injection, finite outputs, no divide-by-zero or unbounded rejection loop; malformed inputs fail explicitly |
| N04 analytical free expansion | Fixed-direction density follows r⁻² outside the source and visibility mask; independent reference profile agrees within 0.1% |
| N05 gravity and boundaries | Ballistic escape/return and Kepler-orbit checks; relative orbit energy error <10⁻³ over the test, converging with timestep; surface and escape losses reconcile with the particle ledger |
| N06 analytic–AMPS comparison | Gravity off and identical source law: resolved density profiles agree within combined 3σ Monte Carlo uncertainty plus 3% discretization target; test domain truncation separately |
| N07 gravity sensitivity | Run gravity on/off for E3/E5 and a northern/slow-particle diagnostic; quantify where the paper's approximation holds instead of imposing global agreement |
| N08 instrument response | Synthetic counts recover a known source, including exposure and chamber-memory cases; zero signal/background, saturated channels and invalid attitudes handled explicitly |
| N09 inverse recovery | Noise-free identifiable synthetic parameters recovered to optimization tolerance; noisy ensemble coverage tested; degenerate sources yield broad limits rather than false precision |
| N10 numerical convergence | Three mesh/timestep levels, particle-count and independent-seed study, expanded domain; important profile metrics change <5% with Monte Carlo error below that change or below observational uncertainty |
| N11 fit and transfer | Publish per-flyby residuals, peak timing/amplitude, integrated signal, uncertainty coverage, and whole-flyby holdout scores; use Section 8 rules |
| N12 reproducibility | One/four-rank and checkpoint/restart cases preserve budgets and statistically equivalent profiles; configuration generation reproduces all requested source parameters |

**Gate G0:** N01–N12 pass for the declared neutral regime, qualified INMS inputs and response exist, the source ensemble and validation split are frozen, and observational limitations are explicit. Historical figure reproduction is a separate flag if exact historical INMS/UVIS products remain unavailable. Do not claim unique source recovery where the data constrain only combinations or upper limits.

## 5. Stage 1 — AMPS neutrals/pickup ions coupled to BATS-R-US plasma

### 5.1 Physics and ownership contract

Begin with neutral H₂O and singly charged H₂O⁺ pickup, then add OH/O/H chemistry and further water-group ions only with documented rates and relevant observational constraints. Photoionization and electron-impact ionization create ion/electron pairs; charge exchange replaces the appropriate charged and neutral states and transfers momentum/energy. Charge exchange must not be implemented as an unconditional increase in ion number. Use reaction rates with units, temperature/energy validity ranges, branching fractions, source citations, and uncertainty.

AMPS advances pickup ions using dv/dt = (q/m)(E + v×B) + gravity and applicable collision/reaction terms. Sample birth velocities from the parent neutral distribution, with physically justified reaction recoil. Evaluate electron-impact rates from the declared electron distribution; temperature alone may be inadequate for a suprathermal population.

BATS-R-US owns the ambient fluid state and field evolution. AMPS owns the kinetic pickup distribution. **There must be exactly one authoritative representation of pickup mass, momentum, energy and current.** Plan the production system as a conservative kinetic–fluid coupling: derive how AMPS pickup moments/current enter the BATS-R-US field/fluid equations and how the equal-and-opposite exchanges enter AMPS. Retain pickup stress anisotropy when required, or quantify a closure that discards it. Evolving a BATS-R-US pickup fluid in parallel with a fully weighted AMPS pickup population is prohibited unless the latter is explicitly diagnostic and excluded from all feedback budgets.

Prescribed-field/one-way operation is the first verification milestone. It is useful for pickup trajectories but does not close the Stage 1 self-consistent mass-loading gate. If the existing BATS-R-US interface cannot support the selected kinetic feedback formulation, its development is an explicit Stage 1 task; do not silently replace physical pickup particles with duplicate fluid mass. Audit the actual equation set and SWMF interfaces at S1.1. The official [BATS-R-US repository](https://github.com/SWMFsoftware/BATSRUS) and [SWMF configuration workflow](https://github.com/SWMFsoftware/SWMF/blob/master/Config.pl) are implementation references, not evidence that this Enceladus coupling already works.

### 5.2 Coupling data and synchronization

| Direction | Required exchange |
|---|---|
| BATS-R-US → AMPS | E, B; per-fluid density, bulk velocity and thermal information; electron closure/temperature/distribution inputs; epoch, frame, units, cell geometry, interpolation validity |
| AMPS → BATS-R-US | Reaction mass/momentum/energy exchanges assigned by owner; kinetic pickup density, momentum, stress and current required by the chosen equations; sampling interval, statistical errors, spatial support |
| Shared | Species/charge mapping, reaction ID, coupling generation, ownership mask, cell volume, units, physical time, external reservoirs and boundary fluxes |

Define a common physical coupling interval, with documented particle/field subcycling and time-centered exchanges. Existing species-local steady-state timesteps cannot be assumed valid for transient feedback. Map extensive quantities conservatively across distinct meshes and convert to densities only using receiving cell volumes. Fields need a suitable interpolation and divergence/error analysis; they are not deposited like particle mass. Invalidate cached maps after regridding and checkpoint coupling histories together.

Specify incoming ambient plasma and field boundaries from qualified upstream intervals, open downstream boundaries, and a surface model for absorption/current closure. Start with a documented simple absorbing body; quantify sensitivity to conductivity/surface charging rather than letting a numerical boundary mimic inferred plume mass loading. Apply domain-size tests along the flow and magnetic field, including the interaction/wake extent.

### 5.3 Implementation work packages

| ID | Implementation | Completion evidence |
|---|---|---|
| S1.1 | Inspect BATS-R-US installation/equations, native SWMF adapters and mapping utilities; write the kinetic–fluid equations and ownership ledger | Reviewed equations, units, interface schema, no duplicate pickup accounting |
| S1.2 | Species metadata, ion birth/charge-exchange/loss network, consistent neutral sinks, Lorentz mover and timestep constraints | I01–I03 |
| S1.3 | Uniform/prescribed-field provider, followed by versioned BATS-R-US snapshots | I04; reproducible time/space interpolation |
| S1.4 | Conservative reaction exchange and kinetic moments; one-cell and two-mesh closed-box coupling | C01–C03 |
| S1.5 | Live two-way coupling, subcycling, regridding, restart, surface and inflow/outflow boundaries | C04–C06 |
| S1.6 | Controlled mass-loading/obstacle benchmarks; progressively stronger plume loading | C07–C08 |
| S1.7 | E3/E5 demonstration, E7/E13 development comparisons, frozen final flyby campaign | C09; INMS/MAG/CAPS/RPWS scores and budget report |

### 5.4 Tests and gate

- **I01 gyromotion:** analytical uniform-B orbit for each ion mass/charge; gyrofrequency, gyroradius and energy error <0.1% at the qualified timestep over ten gyroperiods. Demonstrate temporal convergence.
- **I02 crossed fields:** recover E×B drift and pickup ring/cycloid in a uniform-flow setup. Test the sign of charge and E = −u×B in the selected frame.
- **I03 reactions:** zero-rate, constant-rate exponential survival, branching ratios, and isolated charge-exchange events. Species/charge are exact per event where represented discretely; weighted mass/momentum/energy close to arithmetic tolerance including tracked electron/photon reservoirs.
- **I04 background provider:** independent constant/linear field fixtures, bracketing time interpolation, boundary cells and rotating-frame cases. Missing/stale/out-of-domain fields produce explicit failures or defined masks.
- **C01 exchange mapping:** constant and structured source maps on mismatched meshes conserve integrated transferred quantities to 10⁻¹⁰ relative error in deterministic fixtures; test rank and refinement boundaries.
- **C02 double-count prevention:** one neutral ionization event produces exactly one physical pickup ion and its electron counterpart, with one neutral loss. Charge exchange changes identities without false net ion creation. Diagnostic tracers contribute zero feedback.
- **C03 closed-box budgets:** include reaction energy, field work and external reservoirs; deterministic bookkeeping residual <10⁻⁸ relative to a nonzero reference scale. Separate physical closure error from arithmetic residual.
- **C04 synchronized coupling:** halve coupling interval twice and quantify lag; target <5% change in resolved mass-loading/field/profile metrics on the finest pair.
- **C05 restart/MPI/regrid:** one/four-rank runs, restart during a nonzero-source interval, and mesh changes preserve accumulated exchanges with no missing/double events.
- **C06 zero-source limit:** removing plume/ionization returns the chosen ambient/obstacle reference; disabling feedback returns the separately verified prescribed-field particle solution.
- **C07 mass loading:** homogeneous and one-dimensional reference problems conserve fluxes and reproduce independently derived source-driven deceleration/energy transfer under matched assumptions.
- **C08 full numerical verification:** mesh, timestep, particle count, smoothing, domain and boundary studies; open-system mass and charge residuals <0.1%, momentum/energy residuals <1% of integrated throughput or a documented nonzero reference scale, with convergence.
- **C09 observations:** trajectory-sampled MAG vectors, CAPS response-folded spectra/moments, RPWS density and INMS consistency using the same source ensemble. Record failures as scientific constraints on closures, not as reasons to retune held-out data.

**Gate G1:** two-way feedback and ownership tests pass; coupled solutions converge and meet preregistered observation criteria; neutral constraints remain satisfied within propagated uncertainties. One-way pickup calculations alone close only S1.3, not G1.

## 6. Stage 2 — neutrals, pickup ions and plasma entirely within AMPS

### 6.1 Resolve the plasma formulation before production implementation

The baseline candidate is AMPS electromagnetic PIC using the existing ECSIM infrastructure, with kinetic ambient ions, kinetic pickup ions, and kinetic electrons. This directly implements an AMPS-only plasma model. First establish whether realistic Enceladus scale separation can be resolved at affordable cost. Implicit numerical stability does not remove physical accuracy requirements.

S2.1 must tabulate electron/ion inertial lengths, Debye length, gyroradii, gyro/plasma frequencies, neutral flight times, reaction times, and anticipated interaction scales for upstream, plume and wake conditions. Estimate cells, particles/species, memory, timesteps and wall time from measured small-run costs. Audit field boundary conditions, charge/current deposition, AMR support, chemistry sources and solver tolerances in the actual ECSIM build.

A hybrid alternative—kinetic ambient/pickup ions plus an electron fluid and AMPS-owned fields—can still satisfy AMPS ownership, but it changes the electron physics. If selected for global feasibility, record that choice, derive its generalized Ohm law and electron energy closure, and validate it with local kinetic-electron cases. Kinetic-electron observables must then be labeled closure-dependent. Do not call a prescribed BATS-R-US field replay an AMPS-only self-consistent solution. Neither global full PIC nor an implemented hybrid solver is assumed to be ready from the source inspection.

Reduced mass ratio, reduced light speed, enlarged Debye length, or other scaling changes require their own controlled sensitivity study and an explicit list of dimensionless ratios preserved or altered. A scaled run cannot be relabeled as physical-parameter validation.

### 6.2 Implementation and tests

| ID | Work | Required tests |
|---|---|---|
| S2.1 | Select kinetic-electron or documented hybrid formulation using scale/cost analysis; define governing equations and closures | P01 feasibility and assumptions report |
| S2.2 | Integrate ambient ion/electron initialization and boundary injection, species weights, field solve and current deposition | P02 uniform plasma; P03 oscillations/waves |
| S2.3 | Port reaction ledger from Stage 1; couple neutral ionization, charge exchange and recombination to particle/field/electron evolution | P04 reaction conservation and charge continuity |
| S2.4 | Qualify surface currents, open boundaries, refinement, particle splitting/merging if used, subcycling and restart | P05–P07 |
| S2.5 | Controlled common-physics AMPS versus BATS-R-US benchmark; then restore resolved kinetic physics | P08 |
| S2.6 | AMPS-only flyby campaign with frozen neutral sources and held-out tests | P09; revised validity/cost report |

- **P01 feasibility:** measured memory/time for a pilot and extrapolated full-campaign cost with safety margin. If cost is prohibitive, resolve the formulation before committing to full simulations; this is an open gate, not a passed physics test.
- **P02 equilibrium:** periodic homogeneous plasma shows no secular bulk acceleration or artificial heating beyond the declared budget; total energy residual target <1% over at least ten relevant ion gyroperiods and decreases with refinement/solver tolerance.
- **P03 linear response:** recover plasma oscillations and electromagnetic modes for full PIC, or the supported Alfvén/ion modes for hybrid. Compare dispersion/damping with an independent analytical or linear-solver reference in its applicable regime; resolved phase-speed/frequency errors <5%.
- **P04 charge continuity:** ionization and recombination add/remove consistent electron and ion contributions; charge exchange conserves charge; the discrete continuity residual converges. Monitor Gauss-law residual for full PIC, or quasineutrality/current consistency for hybrid, including boundary currents and sources.
- **P05 boundaries/obstacle:** independent inflow/outflow and absorbing-body tests; expand the domain and test reflected-wave contamination. Demonstrate closure of surface-plus-volume current accounting.
- **P06 AMR:** verify constant/linear field transfers, divergence control, current/charge consistency, and energy changes during refine/coarsen. Quantify effects of particle resampling on moments and noise.
- **P07 repeatability:** one/four-rank and restart statistics; global physical time and field/particle time slots remain synchronized through chemistry and regridding.
- **P08 cross-model:** equal boundaries, neutrals, chemistry, species masses and pressure closure in a fluid-limit benchmark should agree within combined numerical uncertainty, with a provisional 5% target for common moments/fields. In kinetic regimes, publish expected physical differences; agreement with BATS-R-US alone is not observational validation.
- **P09 observations:** same INMS, MAG, CAPS and RPWS operators/splits as Stage 1; test whether resolved kinetic features improve held-out predictions. Hold boundary/source nuisance treatment consistent between stages.

**Gate G2:** AMPS evolves its own plasma and fields for the declared formulation, P01–P09 pass, conserved quantities and field constraints are monitored, and the physical/scaled regime is explicit. Stage 1 remains a retained independent comparison configuration.

## 7. Stage 3 — charged dust added to the AMPS simulation

### 7.1 Add dust in qualified increments

1. **Uncharged passive grains:** restore injection; use source-dependent size and launch distributions, neutral drag and gravity. Verify rate normalization, grain weights and return/escape before charging.
2. **Charged test grains:** fixed charge first, then time-dependent charging in a frozen qualified Stage 2 plasma. This isolates force and charging errors.
3. **Coupled dust:** include collection losses, dust charge/current, momentum transfer and energy exchanges in the plasma/neutral equations where loading warrants feedback. Quantify when passive grains suffice; electron-depletion science requires the appropriate electron collection and charge accounting.

Use each grain's physical radius a, material density, mass m=4πρa³/3, physical charge q, velocity and source identity. Keep statistical weight distinct from physical grain mass/charge. Define whether a power law is dN/da, dN/dln(a), or a mass distribution; normalize to a source number rate or mass rate and expose both. The existing hard-coded 180 kg/s and radius limits are legacy configuration values, not observational constraints.

The dynamical model includes neutral drag, Enceladus gravity and q(E+v×B). Add Saturn tides/rotating-frame terms and radiation pressure where domain and residence-time estimates show relevance. Avoid counting near-vent drag acceleration twice if a prescribed launch law already represents an accelerated outflow. Resolve the small-grain/slow-grain part of the distribution that controls charge and reimpact even if it contributes little mass.

Integrate dq/dt = I_e + ΣI_i + I_photo + I_secondary + other justified currents. Use species-dependent mass, charge, temperature and relative drift. Define illumination/shadowing and the electron distribution. Equilibrium charging is permitted only where charging time is short relative to trajectory/plasma-variation times. Assess discrete charge states for grains with |q|/e of order unity. Document when isolated-sphere capacitance q=4πɛ₀aφ and OML-type currents are valid; do not extend them uncritically to dense, screened, non-Maxwellian or nanograin regimes.

In coupled mode, collection transfers plasma charge to grains: the plasma decrement and grain increment must balance. Account for emitted electrons and photo/secondary energy sources. Include dust in the field/current or quasineutral closure consistently; for negative dust, electron density need not equal total positive-ion charge density. Local and global charge/current budgets must include boundaries and the absorbing surface.

### 7.2 Implementation and tests

| ID | Work | Test and exit evidence |
|---|---|---|
| S3.1 | Size/source schema, injection implementation, grain weighting and radius bins | D01 integrated number/mass moments and sampling; log-bin convention test; zero-dust limit |
| S3.2 | Gravity, drag and fixed-charge motion with appropriate subcycling | D02 analytic ballistic and constant-drag relaxation; D03 charged gyro/drift tests |
| S3.3 | Time-dependent multicomponent charging, equilibrium option, shadowing and small-charge regime | D04 equilibrium-current root; D05 charging transient and sign changes |
| S3.4 | Electron/ion collection, dust deposition/current, momentum and energy backreaction | D06 closed-box and open-system budgets; D07 weak/strong loading |
| S3.5 | CDA/HRD forward operator and E-ring/background treatment | D08 synthetic detection/recovery and actual live-time/attitude masks |
| S3.6 | Joint neutral/plasma/dust inference and held-out flyby campaign | D09 joint scores, uncertainty propagation and identifiability |

- **D01:** independently integrate the size distribution and recover number/mass/charge-weighted moments; deterministic normalization error <10⁻⁴, sampling within 3σ. Confirm rank-independent physical mass injection and proper weights after any velocity-group changes.
- **D02:** neutral-drag-only relaxation against an independent ODE/reference solution with <1% velocity error on the qualified timestep; gravitational escape/return reproduces ballistic references. Gas/dust momentum exchange balances when feedback is active.
- **D03:** fixed-q gyromotion/drift recovers the same charged-particle invariants as I01/I02, using grain q/m; q=0 exactly reduces to the uncharged-grain case.
- **D04:** charging equilibrium agrees with an independent current-balance root to <0.1% potential where the root is nonzero; bound absolute potential/current error near zero. Check positive and negative potential branches and multiple ion masses.
- **D05:** constant-current q(t) and independent transient charging ODE tests; illumination transitions and charge-sign crossings; equilibrium and dynamic solutions agree only in the short-charging-time limit. No negative radius, NaN currents or timestep-dependent sign oscillations.
- **D06:** collection/emission conserve combined particle/plasma/surface charge to the discrete bookkeeping tolerance; apply Section 8 open-system energy/momentum budgets and explicitly track radiative reservoirs.
- **D07:** dust number density →0 recovers Stage 2. Increase loading in controlled tests; verify electron depletion, positivity and closure validity. For hybrid, failure of the quasineutral closure must trigger a limitation or local kinetic treatment, not arbitrary density clipping.
- **D08:** predicted impacts integrate relative grain flux over detector area, look direction, efficiency, speed/size response and live time. Recover known synthetic size populations and test dead time/background. Density alone is not a detector count prediction; HRD counts alone do not uniquely measure grain charge.
- **D09:** qualify size-bin, timestep, charging subcycle, particle-count and domain convergence; score CDA/HRD alongside INMS/MAG/CAPS/RPWS without independently retuning each instrument. Report which dust properties remain prior-dependent.

**Gate G3:** grain injection and charging are operational and verified, the selected feedback regime is justified by loading/charge budgets, detector-level dust predictions meet registered criteria, and previously qualified gas/plasma observables remain consistent or have documented physically explained changes.

## 8. Common numerical and scientific acceptance policy

The numerical values in this plan are proposed engineering targets, not measured accuracy or literature guarantees. Final case-specific observational thresholds must be frozen after data qualification and before fitting/scoring held-out data. Changing a threshold requires a versioned rationale and reevaluation, not silent relaxation after a failure.

### 8.1 Numerical verification

- Use independent analytical solutions, quadrature or reference integrators; do not generate expected values with the production routine under test.
- Separate truncation error, Monte Carlo variance, interpolation error, field-solver error, and coupling lag. Use at least three mesh/timestep levels and independent seeds; increase particle statistics enough to resolve trends.
- For stochastic metrics, record effective independent sample count/autocorrelation and confidence intervals. One/four-rank results need statistical agreement; bitwise identity is required only for deliberately deterministic fixtures.
- Normalize budget residuals using integrated input/output throughput or an explicit nonzero scale, never an almost-zero net flux. Include particles in the domain, escapes, absorption, injection, reaction reservoirs, field energy/work and surface exchanges as applicable.
- Require both local source consistency and global conservation. Opposite mapping errors in different cells can hide in a global sum.
- Avoid using species-local steady-state stepping as transient physical time. Demonstrate synchronized clocks, physical sampling windows and coupling histories for Stages 1–3.
- Production convergence target: <5% change in primary observables on the finest tested pair, with numerical/Monte Carlo uncertainty preferably <one third of observational uncertainty in scored intervals. If that is infeasible, report the achievable numerical uncertainty in the likelihood and do not claim unresolved detail.

### 8.2 Observation metrics

For each instrument and flyby publish absolute/residual profiles with masks and uncertainty, not only a global goodness-of-fit number. Use Poisson/dead-time-aware statistics for counts, covariance-aware residuals for derived continuous measurements, and vector residuals for MAG. Account for shared calibration errors and time correlation.

Primary metrics are peak time/location, peak magnitude, integrated excess above background, profile width/asymmetry, upstream-to-wake changes, magnetic vector perturbations, and CAPS spectral shape/flow where identifiable. Report a goodness-of-fit or posterior predictive score, credible/predictive interval coverage, and residual structure. Compare against a registered simple baseline (e.g. single-source neutral model or no-plume plasma interaction) and account for extra fit parameters.

Suggested provisional science targets: peak timing consistent with exposure/time-registration and propagated source uncertainty; integrated/peak signals consistent with combined measurement and source uncertainty; residuals show no unexplained systematic structure larger than that uncertainty across the primary validation intervals. Freeze numerical score thresholds per product in `quality_rules.yaml` after calibration audit. A factor-of-two neutral-density envelope can be reported for comparison with the 2010 paper, but is not a universal criterion for modern instrument validation.

If held-out data reject a model, retain the failed prediction, identify the missing physics or unconstrained nuisance parameter, and define a new model/data-split version. A physically informative failure is a result; it does not close a passed-validation gate.

## 9. Architecture, artifacts and execution order

### 9.1 Proposed implementation organization

Keep application setup thin and consolidate reusable tested physics where practical. Proposed locations are design targets, not existing APIs:

```text
AMPS/src/models/enceladus/
  sources/       # source schema, geometry, analytical source, sampling
  chemistry/     # reaction tables, rates and conservative event ledger
  coupling/      # provider contracts and conservative exchange adapters
  observers/     # trajectory/LOS sampling; physical distributions and moments
  dust/          # Enceladus-specific injection and closure adapters
AMPS/srcEnceladus--multi-plume/
  adapters/      # application wiring to shared physics
  test/          # staged unit, native, MPI and regression cases
AMPS/input/      # versioned Stage 0/1/2/3 decks through existing parser
cassini_validation_data/
  data/processed/  # immutable qualified observation releases
  validation/      # response processing, scores and evidence bundles
```

Retain `srcEnceladus` as a legacy reference until a replacement reproduces its relevant behavior. Avoid duplicate production source models with diverging parameters. Keep instrument calibration/likelihood processing separable from the PIC solver so synthetic tests and observational provenance are independently reviewable.

Every test/campaign record contains case ID, predecessor gates, source revision/dirty diff, input and generated-code hashes, compiler/MPI settings, executable hash, data release, kernel/calibration hashes, species/closure/ownership configuration, RNG seeds, ranks/threads, exact command, start/stop times, resolution/timesteps/particle counts, metrics, tolerances and PASS/FAIL/SKIP/ERROR outcome. SKIP and unavailable data never count as PASS.

### 9.2 Evidence required at each release

| Release | Required reviewable artifacts |
|---|---|
| R0 data/neutral | Coverage and calibration audit; historical/modern data distinction; fit inputs and posterior; analytic/kinetic tests; source export; INMS holdout profiles; UVIS status |
| R1 coupled | Equation/ownership contract; interface schema; exchange ledgers; native two-way runs; mapping/synchronization/restart evidence; plasma observations and neutral consistency |
| R2 AMPS plasma | Formulation decision; scale/cost assessment; field/particle benchmarks; physical versus scaled parameter statement; cross-model comparison; AMPS-only validation |
| R3 dusty plasma | Grain source/size/charge definitions; charging validity; collection/feedback budgets; detector response; joint holdout results and remaining degeneracies |

Each release includes runnable small examples, exact supported commands, machine-readable metrics, plots and a limitations document. Publication plots should retain data uncertainty, coverage gaps, Monte Carlo error and calibration/model version labels.

### 9.3 Priority backlog and relative effort

Effort classes are relative complexity, not delivery promises: S = bounded adapter/test; M = several interacting routines plus native verification; L = a new physical workflow; XL = solver/coupling research and qualification. Calendar estimates require the first native build, data audit and performance pilots.

| Priority | Task | Effort | Blocking dependency |
|---|---|---|---|
| 1 | Fix coverage assertions; audit INMS formats/calibration, historical products, CDA provenance and SPICE interval coverage | M | G0 observational qualification |
| 2 | Reproduce generated neutral configuration and qualified source geometry/normalization | M | Native neutral runs |
| 3 | Fit E3/E5 and recover synthetic inputs; quantify source degeneracy and response uncertainty | L | Frozen neutral source ensemble |
| 4 | AMPS neutral convergence/gravity/applicability and holdout campaign | L | G0 |
| 5 | Pickup mover/chemistry and prescribed-field reference cases | M–L | S1.3 |
| 6 | Derive and implement conservative kinetic–fluid feedback and live BATS-R-US coupling | XL | G1 |
| 7 | AMPS plasma scale/cost decision, solver integration and common-physics benchmarks | XL | G2 |
| 8 | Restore passive dust, then dynamic multicomponent charging | L | S3.3 |
| 9 | Dust feedback and joint instrument-level validation | L–XL | G3 |

The immediate executable campaign is S0.1–S0.6. Its first deliverable is a qualified neutral input plus an INMS reduction/response fixture, followed by a reproducible E3/E5 fit. No BATS-R-US coupling, AMPS plasma, or dust release should rely on an unqualified neutral source normalization or an unverified instrument conversion.

## 10. Main risks and explicit decision points

| Risk | Diagnostic | Required response |
|---|---|---|
| INMS counts mistaken for density; memory/grain effects mistaken for plume structure | Response audit and synthetic instrument tests | Use documented response/nuisance model; exclude unsupported modes |
| Historical correction applied twice | Calibration lineage comparison with paper Cases 1/2 | Keep historical reproduction separate from modern reduction |
| Source degeneracy or temporal variation | Posterior correlations, source ablation, whole-flyby tests | Report limits/conditional reconstructions; restrict fitted degrees of freedom |
| Frame/mover mismatch in legacy deck | Generated callback inspection and rotation/force tests | Explicit species dispatch, frame and dynamics contract |
| Collisionless model invalid near vents | Local Knudsen/collision-time maps and solver comparison | DSMC or matched injection surface with uncertainty |
| Pickup mass/current double counted | One-event and closed-box ownership tests | Single physical owner; derive feedback equations before live coupling |
| Electron closure controls ionization/charging | Electron-distribution sensitivity and CAPS comparison | Include closure uncertainty or improve electron treatment |
| Global kinetic plasma unaffordable or artificially scaled | S2.1 scale table and measured pilot cost | Document hybrid/local-kinetic alternative and its validity before production |
| Dust current closure invalid at high depletion | Charge balance, positivity and screening diagnostics | Improve closure or limit regime; do not silently clip electrons |
| Broad archive coverage mistaken for encounter detections | Native-record coverage/quality audit | Mark unavailable observables; obtain additional products only where needed |

Source evidence and data references: the two inspected Enceladus application directories and input/parser files listed in Section 2; the local Tenishev PDF; [download manifest](metadata/download_manifest.jsonl), [inventory](metadata/inventory.csv), [verification report](metadata/verification_summary.json), and [volume estimate](metadata/usable_volume_estimate.md). These acquisition records support availability and byte counts; the scientific qualification work remains part of this roadmap.
