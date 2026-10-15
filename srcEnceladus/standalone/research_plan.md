# Research plan: multi-epoch dynamics of the Enceladus neutral plume

Revision: 2026-10-14  
Status: proposed science and implementation plan  
Primary model: AMPS Enceladus multi-plume application and its Stage-0 standalone fitting tools  
Primary observations: Cassini INMS and UVIS, with SPICE geometry; CDA, VIMS, AACS, CAPS, MAG, and RPWS as optional extensions or downstream validation

## 1. Executive summary

The central research objective is to determine which properties of Enceladus's plume are persistent and which change between observations. The analysis will jointly use local INMS traversals and remote UVIS occultations to distinguish:

1. a broad, comparatively slow and diffuse component emitted along the tiger stripes;
2. narrow, directed, higher-Mach jets associated with discrete source locations;
3. an empirical large-scale background where required by the data;
4. species-dependent transport of H2O and CO2; and
5. flyby-to-flyby changes in source rate, collimation, and composition.

The principal scientific hypothesis is that the total H2O output is relatively stable while the division of that output among narrow jets changes with epoch and orbital phase. The alternative hypotheses are that the entire plume varies coherently, that apparent variability is caused mainly by observing geometry and instrument response, or that composition and kinematics change together.

The investigation will not convert every observable to a nominal H2O local density. Each instrument will be compared in its native observable space:

- sticking-corrected local H2O density for the published E3/E5 products;
- prompt mass-44 CO2 density or counts for E7/E14/E17/E18;
- line-of-sight H2O transmission, optical depth, or column for UVIS;
- open-source INMS counts for E21 after applying the angular and velocity response.

The main deliverable will be a cross-validated ensemble of three-dimensional, multi-species plume solutions containing shared source geometry and kinematics together with epoch-dependent source strengths. The ensemble, rather than a single best fit, will be propagated into the later pickup-ion, plasma, and charged-dust stages.

## 2. Scientific questions and testable hypotheses

### 2.1 Primary questions

1. What fractions of the escaping neutral gas are supplied by diffuse tiger-stripe emission, discrete jets, and any required background component?
2. What bulk velocities, thermal spreads, angular divergences, and effective Mach numbers are supported for these components?
3. Which source regions are consistently required by independent INMS trajectories and UVIS lines of sight?
4. Do total production, individual jet strengths, or jet collimation change between epochs?
5. Are the changes correlated with orbital mean anomaly or published tidal-stress predictors?
6. Is CO2 emitted in a constant proportion to H2O, or is the CO2/H2O abundance spatially or temporally variable?
7. What fractions of the emitted gas escape Enceladus, reimpact the surface, or remain temporarily bound?
8. How much uncertainty in the neutral solution must be carried into plasma and dust calculations?

### 2.2 Competing hypotheses

| ID | Hypothesis | Observable prediction |
|---|---|---|
| H0 | Stationary plume | A common source-rate vector and common kinematics explain all observations after geometry and instrument response are applied. |
| H1 | Stable diffuse component plus variable jets | UVIS bulk columns remain comparatively stable while high-spatial-resolution INMS peak amplitudes and selected narrow UVIS features vary. |
| H2 | Coherent total-output variability | Diffuse and jet components scale together with epoch or orbital phase. |
| H3 | Kinematic variability | Peak widths and contrasts require changes in velocity, temperature, or angular divergence even after source amplitudes are allowed to vary. |
| H4 | Composition variability | A common CO2/H2O ratio fails across E3/E5 and later observations, while epoch- or source-dependent ratios improve held-out predictions. |
| H5 | Geometry/instrument explanation | Apparent differences disappear after exact SPICE geometry, exposure integration, inlet response, and calibration uncertainty are modeled. |

H1 is the leading hypothesis, but it must not be imposed by construction. H0 is the required baseline, and more flexible hypotheses must demonstrate improved prediction of observations that were not used to estimate their additional parameters.

### 2.3 Meaning of “dynamics” in this study

Two forms of dynamics must remain separate:

- **Transport dynamics:** motion after emission, including bulk outflow, thermal spreading, gravity, collisions where important, escape, and surface reimpact.
- **Source dynamics:** changes in emission rate, direction, temperature, composition, or collimation between observing epochs.

During one INMS flyby, elapsed time is primarily an along-track coordinate. A feature in a time series normally represents spatial plume structure passing the spacecraft, not demonstrated second-scale source variability. Source variability within a flyby will only be introduced if an independent simultaneous observation or an internally inconsistent spatial model requires it.

## 3. Observational foundation

### 3.1 Core dataset and role of each observation

| Observation | Native constraint | Principal research role | Important limitation |
|---|---|---|---|
| INMS E3/E5 | Published sticking-corrected H2O density; raw mass channels including 44 u | Absolute H2O normalization, diffuse plume width, initial CO2/H2O calibration | H2O inlet memory and grain contamination; relatively coarse source discrimination |
| INMS E7 | Processed mass-28/44 profiles | Intermediate-epoch structure and transfer test | About 12 km along-track resolution; not direct H2O density |
| INMS E14/E17/E18 | High-cadence mass-28/44 profiles | Primary discrete-jet and flyby-to-flyby variability constraints | CO2 is a separate species; source attribution is non-unique |
| INMS E21 | Open-source H2, H2O, and CO2 counts with sensitivity functions | Independent prediction in instrument-count space | No unique count-to-density inversion |
| UVIS 2007 stellar | Resolved H2O attenuation/column profile | Three-dimensional extent and independent H2O constraint | Line-of-sight integral, not local density |
| UVIS 2010 solar, E10 | Highest-resolution H2O occultation profile | Jet widths, projected positions, and total column | Requires exact solar-ray geometry and spectral/count calibration |
| UVIS 2005, 2011, 2016, 2017 | Raw spectra and published H2O column results | Multi-epoch H2O variability and additional viewing geometries | A reproducible reduction pipeline must be completed |
| UVIS non-detections | Upper limits | Control of broad/background plume and false-positive structure | Requires explicit detection thresholds |

E13 is retained for background and geometry studies but is not a primary south-polar plume-structure constraint. Its closest approach was in the northern hemisphere and published plume-structure studies exclude it for that reason.

The local validation archive already contains INMS Level-1A products for E3, E5, E7, E13, E14, E17, E18, and E21, together with SPICE and supporting plasma/dust instruments. The processed SwRI modeling package and raw UVIS occultation products must be placed under versioned, checksummed local storage before analysis. The current literature inventory and exact product identifiers are in [metadata/inms_uvis_literature_search.md](metadata/inms_uvis_literature_search.md).

### 3.2 Required derived data products

Every processed observation will contain:

- UTC, spacecraft clock where available, and SPICE ephemeris time;
- exposure start, midpoint, stop, and integration duration;
- spacecraft state in a named inertial frame and body-fixed IAU_ENCELADUS coordinates;
- instrument look direction, ram angle, and quality flags;
- native measured quantity, units, statistical uncertainty, and systematic-error metadata;
- background estimate and the samples used to obtain it;
- calibration/response version and source-file SHA-256 lineage;
- masks for off-ram observations, saturation, missing attitude, grain spikes, and other exclusions;
- closest-approach-relative time only as a convenience coordinate, never as the sole time reference.

Counts, calibrated species densities, inferred total densities, optical depths, and columns will be stored as distinct fields. A derived density will never overwrite the count data from which it was obtained.

### 3.3 Data products to acquire or construct

1. Archive the approximately 19 MB SwRI modeling release containing `Flybys2.mat`, E7/E14/E17/E18 trajectories, E21 count series and sensitivity functions, and the processed 2007/2010 UVIS series.
2. Extract E3/E5 mass-44 data from the local Level-1A records and synchronize them with the existing digitized sticking-corrected H2O profiles.
3. Acquire the identified UVIS PDS products for 2005, 2007, 2010, 2011, 2016, and 2017.
4. Reproduce published representative H2O columns before accepting the UVIS reduction pipeline.
5. Construct SPICE ray geometry for every UVIS exposure and spacecraft state/attitude for every retained INMS exposure.
6. Record the non-detection chords and their quantitative detection limits.

## 4. Physical model

### 4.1 Source architecture

The common physical source map will contain:

1. discrete jets at published source locations and nominal directions;
2. distributed emission along the tiger stripes;
3. an optional broad or spherical component;
4. an empirical constant background only if observations outside the modeled plume require it.

The empirical background must not be interpreted as a surface production rate. Model comparison will determine whether the data require both a physical broad source and an empirical offset.

For source element `j`, species `s`, and epoch `f`, define

\[
\Theta_{j,s,f}=\{Q_{j,s,f}, u_{j,s,f}, T_{j,s,f},
\hat{\mathbf d}_{j,f}, \Delta\theta_{j,f}\},
\]

where `Q` is production rate, `u` is bulk velocity, `T` is translational temperature, `d` is the nominal emission direction, and `Delta theta` describes nonthermal angular divergence when required.

The baseline velocity distribution is the normalized, outward-directed convected Maxwellian used in the Stage-0 model. Point, line/stripe, and surface-flux emission laws will retain distinct normalization conventions and will be tested independently.

### 4.2 Shared and epoch-dependent parameters

Parameters expected to remain common unless the data reject that constraint:

- surface source position;
- tiger-stripe geometry;
- nominal jet direction within its imaging uncertainty;
- species molecular mass;
- source-law definition and surface boundary condition.

Candidate shared population parameters:

- characteristic jet temperature and velocity distribution;
- characteristic diffuse-source temperature and velocity;
- distribution of angular widths within each source class;
- mean CO2/H2O abundance.

Parameters allowed to vary by epoch in the hierarchical model:

- jet and stripe source rates;
- total plume normalization;
- deviations in collimation or velocity if amplitude variation is insufficient;
- CO2/H2O abundance about a shared population mean;
- documented calibration nuisance factors.

The first fit will not assign every property independently to every candidate jet. Sources will initially be grouped by stripe and by geometric visibility. Individual parameters will only be released when multiple independent observations constrain them.

### 4.3 Species-dependent transport

H2O and CO2 will be transported as separate species. Even at a common source temperature and bulk speed, their thermal speeds differ:

\[
v_{\mathrm{th},s}\propto\sqrt{\frac{T_s}{m_s}}.
\]

Therefore a CO2 density profile cannot be converted point by point to H2O merely by multiplying by an abundance factor. The joint model will use

\[
Q_{j,\mathrm{CO2},f}=r_{j,f}Q_{j,\mathrm{H2O},f}
\]

and will propagate each species with its own molecular mass and permitted thermodynamic parameters. The initial models for `r` will be, in order:

1. one global ratio;
2. one ratio per epoch;
3. one population ratio per source class;
4. source-specific ratios only if demanded by held-out data.

E3 and E5 will constrain the prior distribution of `r`, not impose an exact universal conversion. Published composition estimates and independent density comparisons imply uncertainty substantially larger than counting statistics.

### 4.4 Transport hierarchy

The following hierarchy will prevent unnecessary physics from compensating for source errors:

1. collisionless straight-line analytical transport;
2. collisionless transport with Enceladus gravity and surface loss;
3. AMPS test-particle transport with identical sources;
4. collisional DSMC transport in the near-vent region if the calculated Knudsen number or data require it;
5. matched near-field/far-field model if fully collisional global calculation is unnecessary.

For each epoch and accepted posterior sample, calculate:

- total and species source rates;
- velocity and angular moments;
- local Knudsen number and collision frequency;
- escape, reimpact, and temporarily bound fractions;
- flight-time distribution to each observation;
- surface deposition map;
- three-dimensional density and bulk-velocity fields with uncertainty.

## 5. Observation operators

### 5.1 E3/E5 H2O

The historical baseline will compare the model with the published Teolis sticking-corrected H2O densities using their stated units and digitization uncertainty. A second, independent path will forward-model the raw mass-18 response when sufficient inlet-memory and calibration information is available.

The comparison must account for:

- finite exposure and spacecraft movement during an integration;
- ram enhancement and pointing;
- H2O conversion and wall adsorption/desorption;
- grain-impact contamination;
- correlated uncertainty introduced by deconvolution or smoothing.

The digitized corrected profile and the raw-count likelihood are alternative representations of related data and must not be counted as independent evidence.

### 5.2 Closed-source mass 44 for E3/E5/E7/E14/E17/E18

The preferred observable is the published or reconstructed CO2 density/count series. The operator will:

1. sample the modeled CO2 velocity distribution along the SPICE trajectory;
2. integrate over the exposure;
3. apply ram/pointing enhancement and closed-source sensitivity;
4. add a background term constrained outside the plume;
5. apply quality masks and a count likelihood where counts are retained.

Mass-44 grain spikes will be modeled as contamination or excluded using a rule fixed before plume fitting. They will not be interpreted as narrow gas jets solely because they are sharp.

### 5.3 UVIS

For each UVIS ray at time `t`, compute

\[
N_{\mathrm{H2O}}(t)=\int_0^\infty
n_{\mathrm{H2O}}\left(\mathbf x_{\mathrm{sc}}(t)+
\ell\hat{\mathbf k}(t)\right)d\ell .
\]

For spectral comparisons, predict transmission using

\[
T(\lambda,t)=
\left[\exp\left(-N_{\mathrm{H2O}}(t)
\sigma_{\mathrm{H2O}}(\lambda,T)\right)\right]
\otimes L_{\mathrm{UVIS}},
\]

where the cross-section source, temperature convention, spectral line-spread function, wavelength solution, and detector background are versioned inputs. For HSP or band-integrated products, the same spectral model will be integrated through the relevant response rather than replaced by an undocumented scalar conversion.

The UVIS pipeline will be accepted only after reproducing published columns within their reported uncertainties for at least two geometrically distinct events, including the 2010 solar occultation and one stellar occultation.

### 5.4 E21 open-source INMS

E21 will remain a count-space validation. The predicted signal is schematically

\[
\lambda_{s}(t)=b_s(t)+\Delta t
\int f_s(\mathbf x_{\mathrm{sc}},\mathbf v,t)
\,|\mathbf v-\mathbf v_{\mathrm{sc}}|\,
R_s(\mathbf v-\mathbf v_{\mathrm{sc}},\mathbf A,t)\,d^3v,
\]

where `R` is the species-dependent open-source response, including the SIMION-derived sensitivity and spacecraft attitude `A`. No single scalar count-to-density conversion will be used.

## 6. Inference design

### 6.1 Model sequence

Models will be introduced in the following order:

| Model | Source variability | Kinematic variability | Composition variability | Purpose |
|---|---|---|---|---|
| M0 | None | None | One global ratio | Stationary-plume null model |
| M1 | Per-epoch total scale only | None | One global ratio | Test coherent total-output variation |
| M2 | Stable diffuse source; per-epoch jet amplitudes | None | One global ratio | Leading physical hypothesis |
| M3 | M2 | Source-class or epoch kinematics | One global ratio | Test changing velocity/collimation |
| M4 | M2 or M3 | As supported | Hierarchical epoch/source-class ratios | Test composition variability |
| M5 | Supported model plus orbital/tidal law | As supported | As supported | Test a predictive orbital-phase relationship |

Each additional level must improve held-out predictive skill, not merely reduce residuals on fitted data.

### 6.2 Parameter estimation

At fixed geometry and kinematics, nonnegative source strengths will be solved first using a deterministic linear or constrained optimization step. Nonlinear parameters will then be released in groups:

1. diffuse and jet source-class velocities;
2. source-class temperatures/angular spreads;
3. tightly bounded direction corrections based on imaging uncertainties;
4. composition hierarchy;
5. epoch dependence.

The primary inference will retain full parameter correlations using either posterior sampling or a validated ensemble/profile-likelihood approximation. Multi-start optimization is required because different jet combinations can create similar trajectory peaks.

Useful regularization includes:

- nonnegative source rates;
- hierarchical shrinkage toward a shared source-class distribution;
- sparse deviations from persistent activity rather than sparse total emission;
- smooth distributed stripe emission at the chosen spatial resolution;
- direction priors based on imaging geometry;
- calibration priors derived independently of plume residuals.

Regularization strength will be selected by held-out prediction or synthetic recovery, not by appearance of the fitted curve.

### 6.3 Likelihoods and uncertainty

- Raw INMS counts: Poisson likelihood, or an overdispersed count model only after demonstrating excess variance.
- Published densities: Gaussian or lognormal likelihood according to their processing and error behavior, with covariance or an estimated correlation length for smoothed/deconvolved data.
- UVIS spectra/counts: detector-appropriate count likelihood including reference-spectrum uncertainty and background.
- Published UVIS columns: likelihood using reported statistical uncertainty plus separately represented calibration/systematic terms.
- Non-detections: censored likelihood or explicit detection-limit probability, not a zero-valued measurement.

Nuisance parameters will include instrument sensitivity, background, pointing/time offsets within documented uncertainty, E3/E5 memory correction, UVIS reference intensity, and selected composition factors. Nuisance flexibility will be bounded so that it cannot absorb plume peaks or broad physical structure.

Monte Carlo particle noise will be measured independently and included in the likelihood or made negligible relative to observational errors. It must not become the optimizer's apparent data noise.

## 7. Calibration and validation strategy

### 7.1 Registered data roles

The following is the initial split; a machine-readable manifest will freeze the exact intervals before fitting:

- **Historical reproduction:** E3 and E5 corrected H2O.
- **Composition calibration:** simultaneous E3/E5 H2O and mass 44.
- **Development structure set:** E7 and one of E14/E17 selected before detailed residual inspection.
- **Primary held-out structure tests:** the remaining E14/E17/E18 flybys, rotated in leave-one-flyby-out experiments.
- **UVIS geometry tests:** fit a subset of 2007/2010 and predict the other; later repeat with 2005/2011/2016/2017.
- **Final instrument transfer test:** E21 open-source counts.
- **Upper-limit tests:** UVIS non-detections and qualified off-plume intervals.

Once an observation is used to choose sources, priors, or a model family, it is no longer labeled blind validation. Results will clearly distinguish conditional reconstruction, cross-validation, and true frozen prediction.

### 7.2 Cross-validation experiments

1. Fit E3, predict E5; fit E5, predict E3.
2. Derive a CO2/H2O population prior from one early flyby and test it on the other.
3. Fit E14+E17, predict E18; rotate all three choices.
4. Fit INMS observations and predict UVIS profiles with no UVIS amplitude adjustment.
5. Fit UVIS and E3/E5 H2O, predict the amplitude and shape of later CO2 profiles through the inferred composition distribution.
6. Fit all closed-source data and predict E21 count rates through the open-source response.
7. Fit detected UVIS chords and verify that non-detection chords remain below their thresholds.

### 7.3 Evaluation metrics

Report metrics in physical and statistical form:

- count deviance and predictive log likelihood;
- weighted residual RMS and median absolute residual;
- integrated density/column and its uncertainty;
- peak time/location, amplitude, FWHM, and prominence;
- cross-correlation lag, with any fitted time shift reported separately;
- number of resolved peaks correctly predicted within spatial/time uncertainty;
- coverage of 50%, 90%, and 95% predictive intervals;
- source-rate, mass, and energy bookkeeping;
- information criteria or cross-validated score penalizing unnecessary flexibility.

A model will not be accepted solely because it reproduces integrated production while missing resolved structure, or because it reproduces peaks while violating UVIS columns and upper limits.

## 8. Work packages and decision gates

### WP1 — data qualification and common geometry

Tasks:

1. archive the SwRI package and UVIS products with checksums and licenses;
2. construct common time, frame, and unit conventions;
3. generate INMS spacecraft states and UVIS rays from SPICE;
4. create quality masks and background intervals;
5. publish data dictionaries and provenance manifests.

Gate G1: every accepted record is traceable to immutable raw data; trajectory and line-of-sight geometry pass independent checks; no count product is mislabeled as density.

### WP2 — observation operators

Tasks:

1. implement E3/E5 corrected-density comparison;
2. implement closed-source mass-44 forward sampling;
3. implement UVIS spectral/column reduction and ray integration;
4. implement E21 open-source count prediction;
5. construct synthetic observations for every operator.

Gate G2: synthetic truth is recovered without bias at the expected numerical tolerance, and the UVIS pipeline reproduces published benchmark columns.

### WP3 — stationary multi-species plume

Tasks:

1. reproduce the current E3/E5 H2O fits;
2. introduce CO2 with physically separate thermal spreading;
3. fit M0 to all calibration observations;
4. document structured residuals and parameter degeneracies;
5. compare point jets, stripes, and combined source architecture.

Gate G3: source normalization and transport tests pass; M0 results and failures are reproducible; the selected source architecture improves held-out predictions rather than only fitted residuals.

### WP4 — variable jets and composition

Tasks:

1. fit M1–M4 in the registered sequence;
2. perform E14/E17/E18 leave-one-flyby-out tests;
3. quantify persistent versus epoch-dependent jet activity;
4. estimate population and epoch-specific CO2/H2O ratios;
5. test whether changes in kinematics are actually identifiable.

Gate G4: variability is claimed only if it improves held-out predictive scores and exceeds combined observational, calibration, geometry, and model uncertainty.

### WP5 — orbital/tidal dependence

Tasks:

1. compute mean anomaly and preregister selected tidal-stress predictors;
2. regress total, diffuse, and jet emission separately against those predictors;
3. include observation year to test long-term drift or aperiodic variability;
4. use posterior predictive checks to distinguish a periodic relation from a few influential epochs.

Gate G5: report an orbital relation only if it predicts excluded epochs better than an exchangeable epoch model. Otherwise report an upper limit on supported modulation.

### WP6 — neutral fate and downstream products

Tasks:

1. propagate accepted source ensembles with gravity and surface interaction;
2. calculate escape, reimpact, residence-time, and deposition distributions;
3. determine where collisions affect the inferred observables;
4. export AMPS-ready source ensembles and three-dimensional neutral fields;
5. define uncertainty samples for plasma and charged-dust simulations.

Gate G6: derived fate fractions are numerically converged and stable across the accepted source ensemble; unsupported individual-jet detail is not passed downstream as certainty.

## 9. Verification and test matrix

| Test | Purpose | Acceptance evidence |
|---|---|---|
| T01 source normalization | Verify point, stripe, and surface emission laws | Numerical surface/velocity integral reproduces prescribed `Q`; convergence demonstrated |
| T02 ballistic inverse-square limit | Verify free expansion | Analytical far-field solution recovered for a controlled source |
| T03 gravity/energy | Verify trajectories | Conserved mechanical energy and known radial turning point reproduced at timestep convergence |
| T04 surface fate | Verify escape/reimpact classification | Analytical escape-speed cases and boundary crossings classified correctly |
| T05 species thermal width | Verify H2O/CO2 transport distinction | Sample moments recover the specified mass-dependent Maxwellians |
| T06 SPICE state | Verify geometry | Independent SPICE evaluation agrees for selected epochs and frames |
| T07 UVIS ray integral | Verify column operator | Uniform slab, spherical `1/r^2`, and compact-source integrals match analytical or high-accuracy quadrature results |
| T08 Beer–Lambert retrieval | Verify UVIS reduction | Synthetic columns recovered after spectral convolution and added noise |
| T09 INMS exposure integration | Verify trajectory sampling | Constant and linear fields reproduce exact exposure averages |
| T10 inlet memory | Verify E3/E5 response | Impulse, constant-input, and zero-input behavior conserve the documented response and show no double correction |
| T11 E21 response | Verify open-source operator | Isotropic beam and monoenergetic directional-beam fixtures reproduce independently calculated responses |
| T12 synthetic source recovery | Identify degeneracies | Known source groups and uncertainty coverage recovered across multiple viewing geometries |
| T13 MPI/random reproducibility | Verify AMPS sampling | Statistical agreement across rank counts and seeds; deterministic metadata/output ordering where applicable |
| T14 numerical convergence | Bound discretization error | Trajectory, ray integration, mesh, timestep, and particle-count refinements change scored observables less than the declared numerical budget |
| T15 blind prediction | Demonstrate science validity | Frozen held-out flyby/occultation predictions and metrics recorded before refitting |

Exact numerical tolerances will be strict for deterministic unit tests and based on convergence for stochastic simulations. Science acceptance thresholds will be frozen after the measurement uncertainties are qualified, not chosen after viewing final residuals.

## 10. Identifiability limits and safeguards

The following parameter combinations are expected to be correlated:

- source rate versus bulk speed;
- temperature versus nonthermal angular width;
- source direction versus source strength;
- several nearby jets intersecting the same trajectory segment;
- diffuse stripe emission versus many weak jets;
- CO2/H2O ratio versus flyby calibration scale;
- real time variability versus species-dependent instrumental response.

Safeguards:

1. use independent imaging locations/directions as uncertain geometry, not perfectly known truth;
2. reduce the source basis until synthetic tests demonstrate recoverability;
3. report source groups when individual jets are unresolved;
4. preserve posterior/profile correlations rather than quoting independent one-parameter errors;
5. compare profile shapes and absolute amplitudes separately;
6. publish null-space or singular-value diagnostics for the linear source-strength problem;
7. require multiple geometries before assigning a persistent individual jet;
8. do not interpret a prior-dominated parameter as measured.

The study may robustly resolve a stable diffuse component and groups of variable jets without uniquely determining all published candidate jets. That outcome is scientifically useful and preferable to an over-parameterized source map.

## 11. Expected scientific products

### 11.1 Primary results

1. A quantified decomposition of H2O and CO2 emission among background, distributed stripes, and discrete jets.
2. Credible ranges for source-class bulk velocity, temperature, angular spread, and effective Mach number.
3. A catalog of source regions classified as persistent, variable, unconstrained, or unnecessary.
4. Epoch-by-epoch total and component production rates with correlated uncertainties.
5. Constraints on spatial and temporal variability of CO2/H2O.
6. A test of whether total H2O output is stable while narrow jets vary.
7. A controlled test of orbital/tidal modulation, including a null result or upper limit if warranted.
8. Escape, reimpact, residence-time, and deposition distributions.
9. A neutral-source ensemble suitable for later plasma and dust calculations.

### 11.2 Reproducible technical products

- versioned processed INMS and UVIS datasets;
- a SPICE geometry generator for in-situ trajectories and occultation rays;
- closed-source and open-source INMS observation operators;
- a validated UVIS column/spectral operator;
- machine-readable source schemas and posterior ensembles;
- scripts for every figure and table;
- automated unit, regression, and synthetic-recovery tests;
- a frozen validation manifest and scored prediction reports.

### 11.3 Proposed paper structure

1. Data and native observables.
2. Multi-species plume and instrument-forward model.
3. Stationary-plume baseline.
4. Persistent and variable source components.
5. Composition and species-dependent dynamics.
6. Orbital-phase and long-term variability tests.
7. Escape, fallback, and consequences for the plasma interaction.
8. Identifiability, limitations, and predictions for future missions.

## 12. Risks and mitigation

| Risk | Consequence | Mitigation |
|---|---|---|
| UVIS calibration cannot be reproduced from archive documentation | Later occultations remain scalar constraints | Validate 2007/2010 first; use published columns with explicit provenance while resolving calibration; never invent a missing response |
| E3/E5 CO2/H2O estimates conflict | Weak absolute H2O inference from later mass 44 | Use a hierarchical abundance distribution and report sensitivity across published ranges |
| Too many candidate jets are degenerate | Unstable source map | Group sources, shrink epoch deviations, report resolution kernels/null space, use UVIS geometry |
| Particle noise disrupts fitting | Biased or inefficient inference | Fit an analytical surrogate first; validate it against converged AMPS samples |
| Near-vent collisions invalidate a collisionless source | Biased velocity/width inference | Calculate Knudsen maps; introduce DSMC or a matched exobase/source surface only where required |
| Instrument nuisance terms absorb physics | False fit quality | Bound nuisances independently and evaluate structure in held-out residuals |
| Apparent temporal variability is geometry | Incorrect tidal conclusion | Apply complete SPICE/response operators and test the stationary null model first |
| Absolute production remains uncertain | Weak downstream normalization | Propagate the normalization ensemble and use UVIS/AACS or other independent constraints without double counting |

## 13. Reproducibility and research integrity

All raw observations will remain immutable. Every derived product will record software revision, configuration, input checksums, units, coordinate frames, calibration sources, and exclusion rules. Random seeds and numerical convergence settings will be archived. Figures will be regenerated from tabular products rather than manually edited data.

The following labels will be used consistently:

- **measured:** native instrument observable;
- **calibrated:** observable transformed using an instrument calibration;
- **derived:** physical quantity obtained with additional assumptions;
- **digitized:** values recovered from a publication figure;
- **modeled:** forward-model output;
- **conditional reconstruction:** fit using the observation being shown;
- **cross-validation:** observation excluded during its parameter fit but possibly used in model-family development;
- **blind prediction:** model and scoring rule frozen before the observation is examined.

Negative results, failed transfer tests, excluded records, and parameter non-identifiability will be retained in the research record.

## 14. Principal references

- Tenishev et al. (2010), neutral plume model: <https://doi.org/10.1029/2009JA015223>
- Teolis et al. (2010), E3/E5 H2O inlet response and plume profiles: <https://doi.org/10.1029/2009JA015192>
- Perry et al. (2015), E3/E5/E7/E14/E17/E18 INMS plume data: <https://doi.org/10.1016/j.icarus.2015.04.037>
- Teolis et al. (2017), joint INMS, UVIS, and imaging plume structure: <https://doi.org/10.1089/ast.2017.1647>
- Hansen et al. (2011), 2010 solar occultation: <https://doi.org/10.1029/2011GL047415>
- Hansen et al. (2017), UVIS gas variability: <https://doi.org/10.1002/2016GL071853>
- Hansen et al. (2020), complete UVIS occultation analysis: <https://doi.org/10.1016/j.icarus.2019.113461>
- Portyankina et al. (2022), multi-epoch UVIS jet modeling: <https://doi.org/10.1016/j.icarus.2022.114918>
- Fu et al. (2021), high-resolution CO2 plume structure: <https://doi.org/10.1093/mnras/stab1265>
- Waite et al. (2017), E21 molecular hydrogen: <https://doi.org/10.1126/science.aai8703>
- SwRI INMS/UVIS modeling release: <https://inms-support.space.swri.edu/pages/modeling>

## 15. Immediate next actions

1. Archive and checksum the SwRI modeling package in this workspace.
2. Add the UVIS raw products and their labels to the acquisition manifest.
3. Build an E3/E5 synchronized H2O–mass-44 dataset with documented masks.
4. Implement and test a species-independent SPICE trajectory/ray geometry layer.
5. Implement the UVIS line integral before attempting a joint fit.
6. Extend the standalone Stage-0 model to propagate CO2 separately from H2O.
7. Fit the stationary M0 model and freeze its residual diagnostics.
8. Register the E14/E17/E18 leave-one-flyby-out experiments before fitting variable-jet models.
