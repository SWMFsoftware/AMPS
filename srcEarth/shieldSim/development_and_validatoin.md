# shieldSim development and validation roadmap

Status: **Stage 0 complete; Stage 1 complete; Stage 2 complete; Stages 3 and 4
implemented and acceptance-tested, pending independent review; Stage 5
implemented and reference-tested, pending accelerator/device validation and
independent review; Stage 6 implemented and reference-tested, pending
independent review**
Prepared: 2026-10-16; Stage-6 implementation/evidence updated 2026-10-17
Requirements decisions approved: 2026-10-16 and 2026-10-17
Scope: environment-effect radiation proxies for spacecraft-anomaly analysis
and preliminary mission assessment; ordered multilayer shielding; slab
absorbers initially; published GCR/SEP inputs now; arbitrary absorber geometry
and AMPS input through later adapters  
Release status: **the Stage-4 TID, Stage-5 silicon damage, and Stage-6 generic
LET/hardness implementations have passed their registered automated gates.
Stages 4 and 6 lack independent review; Stage 5 also lacks the project-selected
accelerator/device validation data set. None may be represented as a formally
qualified radiation-effects model.**

The filename retains the spelling requested in the task. This document uses
the normative words **must**, **should**, and **may** deliberately. A stage is
complete only when every stated exit criterion is satisfied and its evidence
is archived. Passing a regression test alone is never evidence of physical
validity.

## 1. Intended result and hard boundaries

### 1.1 Recorded project decisions

The following scope decisions were confirmed on 2026-10-16 and are no longer
open assumptions:

1. The model supports **spacecraft-anomaly analysis** and **preliminary mission
   environment/effect assessment**. It is not an instrument-design model.
2. A shield configuration may be an **ordered multilayer stack**, not only a
   homogeneous slab. Layer order and every interface are part of the physical
   configuration.
3. Absorbers are **finite material slabs in the initial implementation**. The
   data model and scoring interfaces must allow a later geometry provider to
   supply an arbitrary absorber shape without changing the source, transport
   response, or quantity APIs.
4. Requested relative accuracy is an **input field**, with a default of `0.10`
   (10%). Per-quantity overrides are permitted. The result must report whether
   the request was met; the input is not a license to claim accuracy that has
   not been demonstrated.
5. An ordered shield stack contains **one through five layers**. Each layer has
   a maximum geometric thickness of **10 cm**, inclusive. A mass-thickness input
   is converted with the recorded material density and rejected if its derived
   geometric thickness exceeds 10 cm.
6. **All geometry dimensions are supplied in the versioned input file.** No
   shield footprint, layer thickness, absorber dimension, placement, spacing,
   or overburden is inferred from a material name.
7. The default accuracy request is **10% combined standard relative
   uncertainty (`k=1`, nominal 68% coverage)**. A two-sided 95% coverage
   interval must also be reported. Monte Carlo statistical uncertainty remains
   a separately reported component and cannot be presented as the total
   physical accuracy.
8. **TID** is ionizing energy deposited per unit mass in the finite absorber.
   Total and non-ionizing deposited energy are retained as audit diagnostics.
9. **LET** is a step-scored electronic-LET distribution in the absorber, with
   explicitly defined track-mean and dose-mean summaries. A CRaTER deposited-
   energy/chord estimator is a separate, validation-only observable.
10. **`H100/10`** is the ratio of primary-proton shield-exit fluence above
    100 MeV to that above 10 MeV, using total proton kinetic energy. Alpha
    particles and all secondaries are excluded.
11. DDD follows ECSS-E-ST-10-12C and ASTM E3084-17(2022)e1 terminology. It is
    reported as a NIEL exposure proxy, not as device degradation. ASTM
    E722-19(2026) 1-MeV-neutron-equivalent fluence is enabled only for silicon.
12. Stage 2 uses checksum-pinned, downloaded **published spectra as validation
    inputs**. BON2020 is neither a Stage-2 dependency nor a planned production
    runtime component. The downloaded AMS-02 proton component must drive the
    normal case parser, primary generator, full Geant4 transport, and proxy
    writers in an end-to-end test. This is an executable-boundary validation,
    not permission to call one near-Earth Bartels rotation a complete mission
    GCR environment. Inclusive elemental populations remain ingestion-only
    until their isotope composition is supplied. The reference SEP event is
    the GOES-16 **September 10-11, 2017 GLE**, retained both as a time-resolved
    directional series and as a provenance-preserving integrated event
    fluence. A boundary angular transform still requires explicit physics
    justification; the adapter must not silently make the data isotropic.
13. CRaTER is an **instrument-faithful validation harness only**. Its geometry,
    selections, and response do not become the generic absorber model or expand
    the product into an instrument-design tool.
14. Screenshot material names resolve through a versioned material library.
    Every resolved material is explicitly labeled, fully serialized, and may
    be replaced or overridden in the user input. Product-, layup-, porosity-,
    site-, or process-dependent materials remain named trade-study references
    or disabled templates until the required properties are supplied; a generic
    label is never treated as a unique physical material.
15. Requirements are owned by the project owner. Stage 2 source-adapter
    completion is accepted by the project owner from its pinned-data,
    transformation, and test evidence and **does not require an independent
    reviewer**. Later physical capability and observational-validation gates
    retain their stated independent-review requirements.
16. GOES-16 event integration keeps the two SGPS directions separate, does not
    interpolate missing or rejected samples, integrates only valid five-minute
    sample intervals, reports time coverage for every sensor/channel, and marks
    an incomplete fluence as incomplete rather than silently filling gaps. The
    exact published event-onset timestamp is the inclusive integration start and
    the exact published event-end timestamp is the exclusive integration stop;
    both and their source/definition must be recorded before integration. NOAA
    defines the SGPS record coordinate as the start of its five-minute averaging
    interval. The adapter retains that source coordinate and derives the
    midpoint and end, giving half-open support `[start,start + 300 s)` or,
    equivalently, `[midpoint - 150 s, midpoint + 150 s)`. A boundary sample
    contributes only for its exact overlap with the event window. The event-name
    shorthand is not itself a timestamp rule.
17. The Stage-3 implementation decisions approved on 2026-10-16 are:
    (a) Stage 3 validates geometry/transport isolation, raw response matrices,
    particle ancestry, and raw ledgers only; science claims for TID, DDD,
    `n_eq`, LET, and `H100/10` remain disabled until their later gates;
    (b) `shieldSim.case.v2` uses explicit translation and scalar-first unit-
    quaternion local-to-world transforms for the source, shield, and absorber,
    while case v1 remains supported;
    (c) the response library is an authoritative deterministic JSON manifest
    plus a checksum-recorded long CSV, with response and covariance per primary
    and no extrapolation;
    (d) the validation species are proton, He-4, and Fe-56, all transported
    secondaries are retained, and Al/W/HDPE/Si/water are representative test
    materials rather than a claim of complete GCR/material coverage;
    (e) response grids are configurable, interpolation and quadrature are
    declared, grid refinement is tested, and nonlinear ratios fold primitive
    numerator/denominator responses rather than a precomputed ratio;
    (f) `FTFP_BERT` is the reference base, `FTFP_BERT_HP` exercises the neutron-
    HP configuration, and reference-list ion physics is compared as a
    configuration sensitivity with `G4IonQMDPhysics`;
    (g) statistical comparisons pool five independent seeds, require at most
    2% Monte Carlo relative standard uncertainty, require two-sided 95%
    equivalence intervals within +/-5%, require no more than 2% change at the
    last cut/step/size/grid refinement, and treat exact rare/zero outcomes with
    an absolute rule; and
    (h) implementation and tests may be finished without a reviewer, but Stage
    3 is not formally complete until its required independent review is
    recorded.
18. The Stage-4 implementation decisions approved on 2026-10-17 are:
    (a) the validated observable is finite-material-slab TID for proton and
    helium-4 incidence in aluminum, silicon, and water; photon transport is
    limited to X-ray mass-attenuation verification, while electrons, heavier
    ions, organ/device dose, and crew-health outputs remain unvalidated and
    disabled;
    (b) `FTFP_BERT` with its reference-list EM constructor is the production
    baseline, `G4EmStandardPhysics_option4` is an explicit model-form
    sensitivity, and HP/QMD configurations remain diagnostics rather than the
    production TID baseline;
    (c) official NIST PSTAR, ASTAR, XCOM, and X-Ray Mass Attenuation bytes and
    canonical tables are checksum-pinned. Comparisons occur only inside the
    tabulated domains and never extrapolate. Stopping/range and photon
    attenuation are distinct reference checks; attenuation is not relabeled as
    dose;
    (d) accepted transport cases require at most 2% Monte Carlo relative
    standard uncertainty and no more than 2% change at the last cut, step, and
    lateral-escape refinement. Physical reference comparisons use registered
    row-specific limits rather than a universal percentage;
    (e) the uncertainty report separates material/mass, source, numerical, and
    physics-reference components, requires every pairwise correlation and a
    positive-semidefinite covariance matrix, and reports an incomplete budget
    as `not_assessable`;
    (f) deposition by transported particle and propagated ancestry is the
    authoritative breakdown. Creator-process and step-defining-process views
    are diagnostic only and are not unique causal allocations of energy;
    (g) Stage-4 TID is enabled only through case v2. Case v1 remains a
    compatibility/earlier-stage path and does not inherit the Stage-4 claim;
    and
    (h) implementation and acceptance testing may finish without a reviewer,
    but formal Stage-4 completion requires the independent-review record in
    Section 10.
19. The Stage-5 implementation decisions approved on 2026-10-17 are:
    (a) the enabled target is natural silicon only. DDD includes proton,
    helium-4, and neutron contributions; electron, heavier-ion, and non-silicon
    damage responses remain unavailable;
    (b) charged-particle responses are SR-NIEL-7 v11.1 natural-Si tables with
    `Ed = 21 eV`, Coulomb plus hadronic contributions, and no hadronic-threshold
    scaling. No alternative scaling is a hidden default;
    (c) the neutron response is ASTM E722-19, reaffirmed 2026, with the
    prescribed 95 MeV-mb 1-MeV-neutron reference. The enabled neutron domain
    ends at 20 MeV;
    (d) E722 `n_eq` contains the neutron contribution only. Proton and helium-4
    DDD are never divided by the E722 reference or reported as neutron-
    equivalent fluence;
    (e) charged positive nodes use log-log interpolation; neutron values are
    groupwise on explicit boundaries; source folds require explicit bin edges.
    Any nonzero unsupported or out-of-domain contribution makes the affected
    result `not_assessable`; no clamping, extrapolation, or silent omission is
    permitted;
    (f) absorber-entry fluence folding is authoritative. Geant4 nonionizing
    deposition divided by slab mass is retained as a separate diagnostic and
    is never added to the table-folded result;
    (g) shieldSim does not invent an SR-NIEL/ASTM physics-reference
    uncertainty. A combined accuracy request requires a case-supplied
    `physics_reference` uncertainty component with provenance; otherwise it is
    `not_assessable`; and
    (h) implementation and reference testing may finish without external
    device evidence or a reviewer, but formal Stage-5 completion requires a
    project-selected accelerator/device-degradation data set and the independent-
    review record in Section 10.
20. The Stage-6 implementation decisions approved on 2026-10-17 are:
    (a) linear LET is Geant4 unrestricted electronic stopping evaluated at the
    arithmetic midpoint energy of every positive-length charged-particle step
    in the absorber. Its histogram weight is step track length divided by
    absorber volume. Local ionizing deposition remains a distinct TID ledger;
    it is not substituted for stopping power;
    (b) track mean is `sum(L dx)/sum(dx)` and dose mean is
    `sum(L^2 dx)/sum(L dx)`. The primary unit is `keV/um`; the equivalent
    `MeV cm2/g` value uses the resolved absorber density. Both ratio means must
    meet the requested statistical accuracy;
    (c) every charged primary and secondary is scored and identified by
    species, charge, mass, primary/secondary state, root ancestry, and
    generation. Delta electrons are their own component. Positive-length
    boundary/stopping steps are included; zero-length charged steps are counted
    but contribute no track fluence;
    (d) every LET case declares `let_response`. The recommended range is
    0.001--10000 keV/um at 20 logarithmic bins per decade. Underflow and
    overflow are audited and make the histogram result `not_assessable`; values
    are never clamped;
    (e) the p/He-4 stopping-reference domain is natural Al, natural Si, and
    water using the pinned NIST package. The selected heavier ion is fully
    stripped C-12 in natural Si at 10, 100, and 1000 MeV/u against pinned GSI
    ATIMA v1.41 responses with a fixed 15% implementation screen. The initial
    reference-list default missed the 10 MeV/u row by 21.6%; option 4 differed
    by 2.9%, so C-12 LET cases explicitly require
    `G4EmStandardPhysics_option4`. This screen is neither an uncertainty value
    nor experimental validation;
    (f) `H100/10` uses scalar-fluence weight `1/(A |n dot u|)` for every
    accepted forward primary-proton rear-face crossing. Its estimator is the
    ratio of event numerator and denominator means; event variance and
    covariance are retained and the two-sided 95% interval uses Fieller's
    construction. A zero or statistically unresolved denominator is undefined.
    Constant and power-law partial source bins are integrated analytically at
    the strict 10 and 100 MeV thresholds; and
    (g) a complete combined LET or hardness accuracy claim requires a
    provenance-bearing `physics_reference` uncertainty component. Otherwise it
    is `not_assessable`. Stage 6 may be implemented and reference-tested without
    a reviewer, but formal completion requires the independent-review record.

For anomaly analysis, a calculated proxy is an environmental covariate, not by
itself proof that radiation caused an anomaly. Results must preserve time
window, location, source uncertainty, shielding configuration, and the proxy's
physical definition so they can be compared with anomaly timing and other
telemetry. Preliminary mission results are screening evidence, not final
qualification, certification, crew-limit compliance, or detailed instrument
design.

The target product evaluates a matrix of independent cases:

```text
incident environment
  -> source adapter and unit validation
  -> shield transport response
  -> independently defined absorber/instrument response
  -> radiation proxy calculation
  -> uncertainty and provenance report
```

For every selected shield configuration and every selected absorber, the
product is intended to report only the requested quantities. The screenshot
catalogs imply the following coverage:

| Category | Requested choices |
|---|---|
| Shields | Al, Cu, W, Ta, HDPE, BPE, Kevlar, Kapton, CFRP, SiC composite, water, lunar/Mars regolith, silicon |
| Absorbers | skin, eye lens, BFO, CNS, soft tissue, Si, SiO2, SiC, GaAs, InGaAs, Ge, H2O |
| Quantities | TID, DDD, 1-MeV-neutron-equivalent fluence (`n_eq`), LET, `H100/10` |
| Shield input | ordered layers; each layer has material and mass or geometric thickness, with mass thickness in g/cm2 as the primary screenshot-facing form |
| Execution model | one result for each shield-configuration x absorber case |

The following boundaries are mandatory:

1. **AMPS is outside the current implementation scope.** The physics core must
   expose a versioned spectrum contract. A later adapter will translate AMPS
   output into that contract without changing transport or proxy algorithms.
2. **Each shield-configuration x absorber result is an independent
   experiment.** A shield configuration can contain many ordered layers, but
   selecting BFO and Si must not silently place Si downstream of BFO. Absorber
   selection never adds a shield layer.
3. **A material name is not an organ or device model.** Composition, density,
   dimensions, sensitive volume, overburden, and response/cuts must be supplied
   or cited. “BFO,” for example, cannot acquire a physical depth by convention
   without an approved requirement.
4. **Mass thickness and geometric thickness are not interchangeable.** For a
   homogeneous material, `x [g/cm2] = density [g/cm3] * t [cm]`. Both values,
   the density source, and the conversion must be recorded.
5. **Published spectra are not automatically mission environments.** Epoch,
   location, geomagnetic access, angular distribution, species coverage,
   differential variable, bin definition, and whether the values are flux,
   intensity, or fluence must remain attached to the data.
6. **A label is not an implementation or validation claim.** The five quantity
   definitions are frozen in Section 4, but each remains science-disabled until
   its applicable transport, reference-data, uncertainty, and validation gates
   pass.
7. **No causal anomaly attribution, detailed instrument design, medical,
   crew-limit, part-qualification, or flight-safety claim is in scope.** Such a
   claim requires an independently approved application model, uncertainty
   margin, applicable standard, corroborating evidence, and review authority.
8. **No validation parameter may be tuned against its final validation data.**
   Calibration cases and blind validation cases must be declared in advance.

## 2. Baseline audit that motivated the staged work

This section records the legacy behavior found before Stages 1-5. It is retained
for traceability, not as a description of the current versioned case-v2 path.
The stage-status sections below state which gaps have since been closed and
which formal validation gates remain open.

### 2.1 Useful foundation already present

The baseline repository provided a compact Geant4 planar-slab application with:

- CLI-selected Geant4 reference physics lists;
- a material catalog spanning the requested material families;
- beam and cosine-weighted plane-crossing source modes;
- energy-deposition and downstream-face particle scoring;
- deterministic seeds and diagnostic dumps;
- software, geometry, physics-sanity, and regression test layers;
- output paths for TID, DDD, `n_eq`, LET, and `H100/10`.

These are useful engineering foundations. They are not a completed physics
validation.

### 2.2 Gaps that block scientific use

| Area | Baseline behavior | Why it blocked a physical claim |
|---|---|---|
| Spectrum schema | Three columns: total kinetic energy, proton value, alpha value | Cannot preserve bin edges, energy/nucleon or rigidity, direction, epoch, uncertainty, covariance, or additional ions |
| Built-in GCR | Analytic GCR-like proton/alpha approximation | Not a traceable mission environment and lacks heavy GCR ions |
| Retired BON2020-named files | Unverified legacy artifacts with incompatible documented area units | They remain prohibited inputs; Stage 2 uses pinned published observations and does not depend on BON2020 |
| Example SEP | Explicitly synthetic | Suitable for software tests only, not radiation predictions |
| Multiple absorbers | Constructed as an ordered physical stack | Results become selection-order dependent and do not represent independent screenshot cases |
| TID | Uses Geant4 total energy deposition divided by slab mass | TID must separate ionizing from non-ionizing deposition; the current scorer does not record that separation |
| DDD | Folds a transparent analytic NIEL surrogate with particles crossing the shield rear face | It is not an authoritative, species- and material-specific NIEL calculation and does not model energy deposition inside each device |
| `n_eq` | Divides the DDD proxy by a material-dependent placeholder | The standard 1-MeV-neutron-equivalent convention is principally a silicon damage convention; extension to other media requires an explicit standard/reference |
| LET | Maps transmitted proton/alpha energy to a simplified Bethe stopping-power estimate | It is neither a step-scored absorber LET distribution nor a CRaTER-equivalent deposited-energy/track-length response; corrections and heavy ions are absent |
| `H100/10` | Uses transmitted p + alpha + neutron counts at total-particle-energy thresholds | Species, energy basis, flux/fluence basis, and measurement location are not specified by the label |
| Response location | DDD/LET use the shield rear-face spectrum | A detector/organ response generally depends on transport through and within the absorber, not only the upstream boundary spectrum |
| Material catalog | Several composites, regoliths, and organ proxies are documented approximations | Mission-quality composition, density, moisture/porosity, layup, and uncertainty are not defined |
| Uncertainty | Monte Carlo outputs and rates are written without a complete uncertainty budget | Statistical precision cannot stand in for source, material, cross-section, geometry, and response uncertainty |
| Validation | Regression and qualitative sanity tests exist | A regression preserves behavior, including incorrect behavior; no accepted analytic, tabulated, or observational validation baseline exists |

The first development release had to disable or prominently watermark
science-facing DDD, `n_eq`, LET, and hardness outputs until the corresponding
stage gates below are passed. The approximate functions may remain
available only under names containing `prototype` or `surrogate`.

## 3. Required architecture

### 3.1 Canonical incident-environment contract

All source readers, including the future AMPS reader, must create the same
validated `SpectrumDataset`. At minimum it must contain:

- dataset identifier, version, DOI/URL, retrieval time, license, and checksum;
- observation/model interval and spatial frame/location;
- particle identity: PDG/Geant4 name plus `Z`, `A`, charge state if relevant;
- independent variable type: total kinetic energy, kinetic energy per nucleon,
  momentum, or rigidity;
- bin lower/upper edges and, when supplied, representative/effective energy;
- differential variable and exact units;
- value meaning: differential intensity, directional flux, plane-crossing
  flux, omnidirectional flux, or fluence;
- angular representation and coordinate frame: beam, isotropic hemisphere,
  directional bins, pitch-angle bins, or an explicit angular map;
- statistical and systematic uncertainties kept separately when available;
- covariance/correlation representation when available;
- data-quality flags and fill-value rules;
- transformations already applied, including their equations and version;
- source metadata that could not be interpreted, retained verbatim.

The core must reject incomplete or dimensionally inconsistent data. It must not
guess units from magnitude, infer upper edges from a field named “low,” convert
a directional measurement to isotropic intensity, or replace missing species
with zero without an explicit scenario definition.

Canonical internal units should be chosen once and enforced through typed
quantities or dimension-checked conversion functions. Human-readable output
must repeat the original and canonical units.

### 3.2 Source adapters

Adapters must be thin, independently tested translators:

- `LegacyThreeColumnAdapter`: supports existing synthetic/monoenergetic tests,
  but requires a sidecar or CLI declaration for energy basis, value meaning,
  units, and angular model;
- `AMS02HapiAdapter`: imports the pinned measured GCR fixture after confirming
  the bin-coordinate semantics against its pinned primary supplemental table;
- `GOESRSGPSAdapter`: imports channel bounds, two sensor directions,
  uncertainties, samples, quality flags, and yaw state without angular
  averaging, then integrates valid records with explicit coverage ledgers;
- `PublishedTransportAdmission`: admits only a specifically selected component
  whose particle identity, energy coordinate, value units, angular meaning,
  bin edges, and validity are sufficient for transport. The initial admission
  is the A=1 AMS-02 proton component in validation-only isotropic cases;
- `AMPSSpectrumAdapter` (later): validates AMPS schema/version, grid,
  normalization, species, energy variable, frame, time interval, and angular
  representation before constructing `SpectrumDataset`.

The AMPS adapter must be a separate library/module. No AMPS file-format logic
may enter the transport, response-matrix, or quantity implementations.

### 3.3 Case definition and geometry

A `CaseDefinition` must identify:

1. incident dataset and angular boundary condition;
2. one shield configuration;
3. one absorber or instrument configuration;
4. requested quantities;
5. physics list, Geant4/data-library versions, cuts, step controls, and random
   seed policy;
6. global requested relative accuracy, optional per-quantity overrides, and
   event/convergence policy;
7. a machine-readable statement of whether the accuracy request refers to
   Monte Carlo precision or the combined reported uncertainty.

The screenshot selection matrix is expanded into independent case definitions.
No material in one case can affect another case.

#### 3.3.1 Ordered multilayer shield contract

`ShieldConfiguration` must contain a `ShieldStack` with a stable identifier,
schema version, declared incident side, reference surface/normal, transverse
extent, and an ordered list of one or more `ShieldLayer` objects. The list order
is the physical source-to-absorber order and must never be sorted by material or
thickness.

The stack must contain between one and five layers. The geometric thickness of
every layer must be greater than zero and no greater than 10 cm. If a layer is
specified by mass thickness, the parser must derive its geometric thickness
from the versioned nominal density before geometry construction and apply the
same 10 cm limit. Five layers at the maximum imply a possible 50 cm total
geometric stack, but no total mass-thickness limit is inferred because it is
material dependent.

Every layer must contain:

- a stable layer identifier and zero-based order index;
- a versioned material definition, including composition, density, state, and
  their provenance/uncertainties;
- exactly one user-specified thickness representation: mass thickness in
  `g/cm2` or geometric thickness with an explicit length unit;
- the derived thickness representation, calculated from the recorded density;
- transverse footprint and alignment, initially inherited from the planar
  stack unless explicitly overridden by a later geometry type;
- optional manufacturing metadata such as lot/configuration identifiers, which
  do not alter physics unless mapped to a material property;
- a content hash included in the full stack/response hash.

Layer boundaries are contiguous in the initial planar implementation. Adhesive,
coating, fabric/resin, void, stand-off, or vacuum regions that are physically
important must be explicit layers with their own material and thickness; the
model must not insert or infer them. Zero or negative thickness, missing
density, overlapping placement, and unresolved material keys are fatal input
errors. Repeated adjacent materials are allowed because they can represent
separately documented components.

For normal incidence the reported total mass thickness is the sum of layer
mass thicknesses. For oblique particles, Geant4 geometry determines the path
through each layer; the stack must not be replaced by a total equivalent areal
density or homogenized material. Layer order is therefore part of the response
matrix key even when total areal density is unchanged.

The initial shield geometry is a co-axial planar stack. Directional areal
thickness maps, curved vehicle shells, and regolith half-spaces are separate
geometry providers and require their own validation. A regolith half-space is
not represented by silently creating a thick finite slab.

Multilayer-specific verification must include:

- exact interface positions, layer volumes/masses, total physical thickness,
  and total normal-incidence mass thickness from a hand-calculated stack;
- serialization round-trip without layer reordering or unit loss;
- distinct configuration hashes for `A|B` and `B|A`;
- a one-layer stack statistically reproducing the former single-slab case;
- splitting one homogeneous layer into two adjacent layers of the same
  material statistically reproducing the unsplit case, subject to identical
  cuts and boundary handling;
- explicit gap-layer behavior and rejection of implicit gaps/overlaps;
- interface-crossing diagnostics in local and global coordinates;
- production-cut and step convergence near every material boundary;
- direct-spectrum and response-matrix agreement for at least one stack with a
  high-Z layer and one hydrogenous layer in each order;
- separate component scoring so secondary production in each layer can be
  traced to downstream proxies.

No expected physical ordering effect may be asserted without an analytic,
tabulated, experimental, or independently simulated reference. The software
test requires that order is represented and traceable; the physics-validation
test determines whether its predicted effect is accurate.

#### 3.3.2 Initial slab absorber

The initially enabled `SlabAbsorberGeometry` must define:

- material and material-version identifier;
- geometric thickness and transverse dimensions with units;
- density and computed volume/mass;
- position, orientation, and distance from the last shield layer;
- the sensitive region, initially the entire slab unless explicitly partitioned;
- stable scoring-region identifiers;
- overburden, if any, as explicit geometry rather than a name-based convention;
- scoring policies for deposited energy, non-ionizing deposit, step/track
  length, particle identity, ancestry, entry/exit state, and uncertainty.

The program must require slab dimensions; it must not derive organ or device
dimensions from the material label. Tissue endpoints remain material-slab
proxies, not anatomical models. Electronic materials remain material-slab
proxies, not devices. Each selected absorber produces a separate case/world so
one absorber never shields another.

The input file must provide the transverse dimensions, source plane, world
boundary, layer footprint, shield-to-absorber gap, absorber offset, and all slab
dimensions. The accepted input ranges must be confirmed by convergence tests.
There are no hidden geometry defaults.

#### 3.3.3 Versioned case-file input schema

The versioned case file is the source of truth. CLI/UI controls may generate or
modify it, but the exact serialized case must be archived with the result. The
schema must express at least the following fields; this is a semantic sketch,
not a commitment to YAML rather than JSON or another reviewed serialization:

```text
schema_version
accuracy.default_relative_target = 0.10
accuracy.per_quantity[TID|DDD|n_eq|LET|H100_10] = optional override
accuracy.scope = combined                         # approved default
accuracy.standard_relative_coverage = 0.68
accuracy.report_interval_coverage = 0.95

material_library.id_and_version
materials[] = optional immutable per-case definitions/overrides

geometry.coordinate_frame
geometry.world_dimensions
geometry.source_plane_dimensions_and_transform

shield_stack.incident_side
shield_stack.reference_surface_and_normal
shield_stack.footprint
shield_stack.layers[1..5]:
  layer_id, order, material_id
  exactly one of mass_thickness or geometric_thickness
  density_source/version/uncertainty
  derived thickness and material/content hash

absorber.geometry_provider = slab
absorber.material_id
absorber.slab_dimensions
absorber.transform
absorber.shield_gap
absorber.sensitive_region
absorber.explicit_overburden
```

Parsing must produce a fully resolved `CaseDefinition` or fail. Every library
entry and per-case material override must be expanded into its canonical
resolved record and hash before any mass-thickness conversion. A missing
dimension, unit, coordinate frame, transform, density needed for conversion,
or required provenance field is not filled from a material label or a previous
run. The approved accuracy defaults are 10% combined standard relative
uncertainty (`k=1`, nominal 68% coverage) and a reported 95% interval; all
other physical inputs remain explicit.

#### 3.3.4 Geometry-neutral absorber extension point

Transport and scoring code must depend on an `AbsorberGeometry` contract, not
on `G4Box` or slab copy-number assumptions. A geometry provider must expose the
following logical operations regardless of its internal representation:

```text
Validate()                         -> structured validation report
Build(mother_volume)               -> Geant4 placements and region map
SensitiveRegions()                 -> stable region IDs and scoring policies
MassProperties(region_id)          -> volume, mass, centroid, uncertainty
BoundsAndReferenceFrame()          -> bounds plus frame/transform metadata
CanonicalGeometryDescription()     -> deterministic serialization
GeometryContentHash()              -> response-cache/provenance key
```

`SlabAbsorberGeometry` is the only provider enabled initially. The schema must
reserve versioned provider types for later additions without pretending they
are equivalent:

- analytic/constructive solid geometry primitives and Boolean combinations;
- tessellated CAD/mesh geometry;
- voxelized material/region maps;
- GDML or another reviewed Geant4-native geometry description;
- an application-specific observer/spacecraft geometry adapter.

Every future provider must supply explicit length units, coordinate frame,
origin, orientation, scale, region-to-material mapping, sensitive-region
mapping, source file version/checksum, and mass-property provenance. Scorers
must attach results to stable logical region IDs returned by the provider,
never to incidental Geant4 volume names or copy numbers.

A response-cache key must include the canonical geometry content hash,
placement transform, material hashes, sensitive-region definitions, scoring
policy, and geometry-provider version. A shape or placement change invalidates
the response even if the material and nominal mass are unchanged.

Arbitrary-shape import must fail before transport when any applicable check
fails: missing/ambiguous units, non-finite coordinates, non-watertight or
self-intersecting solid, inconsistent surface orientation, duplicate/overlapping
regions, unmapped material, empty sensitive region, invalid transform, Geant4
overlap, or a mass/volume inconsistency outside an approved tolerance. Repairing
a mesh automatically is prohibited unless the repair is an explicit,
versioned preprocessing operation whose input/output hashes and geometric
changes are reported.

The future shape-validation ladder must include:

1. exact volume, surface intersection, chord-length, and mass tests for spheres,
   cylinders, boxes, and simple Boolean solids;
2. the same primitive represented analytically, as a mesh, and as successively
   refined voxels, with convergence of mass and radiation proxies;
3. translation/rotation invariance tests for an isotropic environment and
   predicted orientation dependence for a directional environment;
4. reference-ray navigation tests through every material/sensitive-region
   boundary;
5. energy-deposition and track-length ledger closure by logical region;
6. scoring invariance under harmless tessellation changes;
7. convergence against mesh/voxel resolution and Geant4 geometry tolerance;
8. a direct-vs-response-matrix comparison using the identical geometry hash;
9. an independent reference for at least one nontrivial shape before that
   provider can produce science-facing results.

An “equivalent slab” may be run as a separately labeled comparison, but it may
not stand in for an arbitrary shape or be used to validate that shape.

#### 3.3.5 Versioned, user-overridable material library

Material definitions are physics inputs, not display metadata. The production
catalog must be a versioned data file loaded before geometry construction; the
hard-coded builders in `src/MaterialCatalog.cc` are migration inputs, not the
future source of truth. Every material instance must have a stable library ID,
library version, status, applicability statement, property-level provenance,
and a canonical content hash.

The library uses the following statuses:

- `reference`: a reproducible reference material with enough sourced
  properties to instantiate it for transport. This status does not claim that
  it matches a particular spacecraft part or biological structure.
- `trade_study`: a fully explicit reference recipe usable for sensitivity and
  preliminary comparisons, but product, process, site, or condition differences
  prevent it from being treated as mission-specific.
- `template`: required fields are intentionally unresolved. It cannot be used
  in transport until a user supplies and cites them.
- `project_specific`: a user or project definition with measured/supplier
  properties, uncertainties, and approval record. Only this status may be used
  to claim representation of a named flight article or material lot.

The resolved material record must contain the following properties when they
are applicable:

| Property | Required representation and reason |
|---|---|
| Identity | stable ID, revision, explicit label, permitted shield/absorber roles, status, and applicability limits |
| Bulk state | density with unit and uncertainty, solid/liquid/gas state, reference temperature, and pressure where relevant; density is required for geometry, dose mass, and mass-thickness conversion |
| Elemental composition | exactly one declared basis: elemental mass fractions or integer/real atom counts; normalization rule and residual/unreported fraction must be explicit |
| Isotopes | natural abundance tied to the installed Geant4/NIST data version, or explicit isotope fractions; boron enrichment and similar changes may not hide behind an elemental label |
| Compound/composite recipe | constituent IDs and mass/volume basis, resin/filler or fiber/matrix fraction, void/porosity fraction, moisture/volatile fraction, and the deterministic conversion to final elemental mass fractions |
| Electromagnetic metadata | mean excitation energy only when explicitly overridden from the transport library, with source; otherwise record the value resolved by the pinned Geant4 material implementation |
| Low-energy neutron metadata | temperature and thermal-scattering-law association where supported and material-appropriate; absence must be recorded for neutron-sensitive results |
| Damage-response metadata | a separate NIEL dataset/version, target identity, displacement-threshold convention, incident species and energy domain, interpolation rule, and out-of-range policy; it is not inferred from elemental composition |
| Uncertainty/provenance | source URL/DOI/report, edition/version/page or table, retrieval date, property uncertainty/covariance when available, license, and any explicit transformation |

Mechanical strength, modulus, weave orientation, or grain size is not directly
consumed by homogeneous Geant4 particle transport. It must nevertheless be
retained when it determines effective density, porosity, constituent fraction,
or whether homogenization is physically acceptable. A weave, honeycomb, void,
or macroscopic grain structure that is not well represented by a homogeneous
mixture must be expressed as geometry in a later provider, not hidden in an
effective material name.

The case file may extend a library material or define a new one. Overrides must
be explicit data, never a mutation of the shared catalog entry. The required
semantic form is:

```text
material_library.id_and_version
materials[]:
  material_id                         # new immutable ID for this case
  extends                             # optional exact library ID@version
  status                              # normally project_specific or trade_study
  label
  density: {value, unit, uncertainty, reference_condition, source}
  state: {phase, temperature, pressure}
  composition:
    basis = elemental_mass_fraction | atom_count | constituent_mass_fraction |
            constituent_volume_fraction
    components[] = {element/isotope/material_id, value, uncertainty}
    residual_policy
  optional mean_excitation_energy and thermal_scattering_law
  optional damage_response_links[]
  provenance: {source, version, retrieval_date, transformation, approval}
  applicability_and_limitations
```

All transport-relevant properties are therefore changeable in the user input,
including density, composition, isotope abundance, state, temperature, and a
referenced damage-response association. An override produces a new resolved ID
and content hash; the output preserves both the parent record and exact changes.
Changing composition invalidates any inherited mean-excitation-energy override,
thermal-scattering association, NIEL table, response matrix, and validation
claim unless compatibility is re-established. Changing density invalidates the
geometry derived from mass thickness and every affected response matrix.

User inputs must be rejected for negative/non-finite values, unknown elements
or isotopes, duplicate components, inconsistent fraction bases, mass fractions
outside the declared normalization tolerance, unresolved constituent cycles,
unsupported state/temperature combinations, missing density, or missing
property provenance. Renormalizing an incomplete published analysis is allowed
only as a named transformation that retains the original sum and unreported
fraction; it may not happen silently.

##### Initial catalog disposition

The following is the approved starting disposition for the screenshot labels.
Numerical values identify a reference recipe; they are not hidden defaults for
a flight material. Standard entries follow the pinned [Geant4 11.4 material
database](https://geant4.web.cern.ch/documentation/dev/bfad_html/ForApplicationDevelopers/Appendix/materialNames.html).

| UI label | Initial explicit library record | Status and required use |
|---|---|---|
| Al | natural Al, `2.699 g/cm3`, solid | `reference`; project alloy/temper and measured density may override it |
| Cu | natural Cu, `8.96 g/cm3`, solid | `reference` |
| W | natural W, `19.3 g/cm3`, solid | `reference` |
| Ta | natural Ta, `16.654 g/cm3`, solid | `reference` |
| HDPE | `HDPE_CH2_0p950_trade`, CH2 repeat-unit composition, `0.950 g/cm3` | `trade_study`; the density is tied to the named [DOW 12450N HDPE](https://www.dow.com/content/dc/br/pt/pdp.dow-12450n-high-density-polyethylene-resin.47054z.html) reference, while the library retains `G4_POLYETHYLENE` at `0.94 g/cm3` separately; mission use requires supplier/grade/lot density, comonomers/additives, and uncertainty |
| BPE | `BPE_5wtpct_naturalB_trade`, 95 mass % `HDPE_CH2_0p950_trade` plus 5 mass % natural B, reference density `1.04 g/cm3` | `trade_study`; explicitly not a universal BPE. The density is consistent with a named 5% product, but filler chemistry and isotopic enrichment must be supplied for a flight product ([JCS JCPE5 product data](https://johncaunt.com/products/jcpe5/)) |
| Kevlar | para-aramid repeat unit C14H10N2O2, `1.44 g/cm3` | `reference` for dry fiber only, consistent with `G4_KEVLAR` and the [DuPont Kevlar technical guide](https://www.dupont.com/content/dam/aramids/amer/us/en/safety/public/documents/en/Kevlar_Technical_Guide_0319.pdf); cloth/resin/moisture/void content requires a separately named assembly |
| Kapton | C22H10N2O5, `1.42 g/cm3` (`G4_KAPTON`) | `reference`; a coated film or adhesive stack is not this material |
| CFRP | `CFRP_template` | disabled `template`; requires fiber identity, resin elemental composition, mass/volume fraction basis, cured density, void/moisture content, and uncertainty. [Hexcel HexPly 8552](https://www.hexcel.com/wp-content/uploads/2025/12/HexPly_8552_eu_DataSheet1.pdf) supplies product-specific fiber/resin/laminate densities and fractions, but not enough resin chemistry to invent a transport composition |
| SiC composite plastic | `SiC_composite_template` | disabled `template`; requires SiC polytype/filler, matrix chemistry, loading basis, density, porosity, additives, and uncertainties |
| Water | H2O, `1.0 g/cm3`, liquid (`G4_WATER`) | `reference` at the catalog reference condition; temperature/density and dissolved constituents may be overridden |
| Lunar regolith | `LunarRegolith_template` | disabled generic `template`; must resolve to a named Apollo/site composition or named simulant plus packing state. NASA reports material/site dependence and depth-dependent bulk density; e.g. the top 15 cm reference is `1.50 +/- 0.05 g/cm3`, not a universal lunar density ([NASA/TP-20220018746](https://ntrs.nasa.gov/api/citations/20220018746/downloads/TP-20220018746.pdf)) |
| Mars regolith | `MarsRegolith_template` | disabled generic `template`; must resolve to a named site/APXS analysis or named simulant with packing/moisture state. The [NASA PDS APXS archive](https://pds.nasa.gov/ds-view/pds/viewDataset.jsp?dsid=MSL-M-APXS-4%2F5-RDR-V1.0) is the preferred provenance path for in-situ chemistry |
| Silicon shield | natural Si, `2.33 g/cm3`, solid (`G4_Si`) | `reference` |

The absorber labels resolve as follows. Tissue labels describe only slab
composition; they do not encode anatomy, depth, biological weighting, or a
human phantom.

| UI label | Initial explicit library record | Status and required use |
|---|---|---|
| Skin | `G4_SKIN_ICRP`, `1.09 g/cm3`, published nine-element mass fractions | `reference` tissue-composition slab |
| Eye lens | `G4_EYE_LENS_ICRP`, `1.07 g/cm3`, published eight-element mass fractions | `reference` tissue-composition slab |
| BFO | `BFO_soft_tissue_proxy` extending `G4_TISSUE_SOFT_ICRU-4`, `1.00 g/cm3` | explicitly labeled `trade_study` material proxy; BFO is not a chemical material and geometry/overburden remain user inputs |
| CNS | `G4_BRAIN_ICRP`, `1.04 g/cm3`, published nine-element mass fractions | `reference` composition slab, not an anatomical CNS model |
| Soft tissue | `G4_TISSUE_SOFT_ICRU-4`, `1.00 g/cm3` | `reference` four-element tissue slab |
| Si | natural Si, `2.33 g/cm3` (`G4_Si`) | `reference`; silicon is the only initially enabled ASTM E722 `n_eq` target |
| SiO2 | stoichiometric `G4_SILICON_DIOXIDE`, `2.32 g/cm3` | `reference` for that tabulation; amorphous/fused or process oxide at another density must use a different/overridden record |
| SiC | stoichiometric 6H-SiC, `3.214 +/- 0.001 g/cm3` at 20 C | `reference`, sourced to the [NIST SiC review](https://srd.nist.gov/JPCRD/jpcrd529.pdf); another polytype, porosity, or device film requires an override |
| GaAs | stoichiometric `G4_GALLIUM_ARSENIDE`, `5.31 g/cm3` | `reference` |
| InGaAs | `Ga0.47In0.53As_300K`, atom ratios 0.47:0.53:1 and `5.50 g/cm3` | `reference` for this explicitly named alloy only, from the [Ioffe NSM compilation](https://www.ioffe.ru/SVA/NSM/Semicond/GaInAs/basic.html); the alloy fraction and density are user-overridable |
| Ge | natural Ge, `5.323 g/cm3` (`G4_Ge`) | `reference` |
| H2O | H2O, `1.0 g/cm3`, liquid (`G4_WATER`) | `reference` |

This audit exposes two required Stage 1 corrections: the current custom SiO2
record (`2.20 g/cm3`) must either be relabeled as a particular amorphous/fused
silica or replaced by the explicit `G4_SILICON_DIOXIDE` reference; and the
current `In0.53Ga0.47As` record uses `5.68 g/cm3`, whereas the cited 300 K
composition is approximately `5.50 g/cm3`. The present CFRP, SiC-composite,
lunar-regolith, and Mars-regolith hard-coded mixtures must not migrate as
`reference` records.

Material-library V0/V1 tests must verify schema validity, fraction and isotope
normalization, constituent-to-element conversion, property-level citations,
status gating, exact override precedence, round-trip serialization, and content
hash changes for every transport property. A Geant4 material dump must match
the resolved record. Boundary tests must prove that density overrides change
mass-thickness conversion and absorber mass, that composition changes
invalidate response/NIEL links, and that templates or generic aliases cannot
run. Reference-material values must be spot-checked against the cited source
and installed Geant4 version; these are independent data tests, not regression
snapshots generated by shieldSim.

### 3.4 Response-matrix transport model

The preferred architecture separates expensive Geant4 transport from changing
incident spectra. For each shield/absorber/physics configuration, compute a
monoenergetic response matrix

```text
R(q, output bin | species, input energy bin, input angular bin, case)
```

and fold it with an incident `SpectrumDataset`. The response must retain
primary and secondary species contributions. A response library records all
geometry, material, physics, cut, Geant4, nuclear-data, and statistical
metadata.

This structure is important for AMPS: a later AMPS dataset changes the source
vector, not the validated response implementation. It also enables exact
linearity and unit tests.

At selected points, a direct spectrum-sampled Geant4 simulation must agree with
the response-matrix fold within combined Monte Carlo uncertainty. Importance
sampling or biasing may be added only after unbiased equivalence tests pass.

### 3.5 Species and physics coverage

Protons and helium are insufficient for a quantitatively accurate GCR LET or
dose-equivalent environment. A production GCR scenario must retain the heavy
ion species required by its stated domain and approved accuracy target; the
Stage-2 proton-only transport case tests software integration, not that
coverage. SEP species coverage must be decided from the selected event product;
the current GOES fixture is proton-only.

Transport must create and score all relevant secondary charged particles,
neutrons, photons, electrons/positrons, ions, and nuclear fragments. Output may
aggregate them only after preserving component results.

At least one high-precision neutron configuration must be evaluated because
hydrogenous shielding and regolith can make secondary neutrons important.
Electromagnetic, hadronic, ion, neutron, decay, production-cut, and step-limit
choices must be stated in every result. Reference physics lists are starting
points, not validation evidence.

### 3.6 Result contract

Every scalar, spectrum, or matrix result must include:

- quantity name and unambiguous definition/version;
- units and normalization (per primary, per incident fluence, rate, or
  time-integrated result);
- source and response identifiers;
- shield and absorber definitions;
- particle-component breakdown;
- statistical uncertainty and event/effective-sample count;
- propagated source uncertainty where available;
- identified systematic uncertainty contributions or an explicit statement
  that they remain unevaluated;
- code revision, build flags, Geant4 and data-library versions, random seeds,
  platform, and run manifest;
- validation status of that exact quantity/domain.

CSV/JSON or HDF5 should be the canonical machine-readable output. Plotting
formats may be derived from it, not treated as the data model.

### 3.7 Time-resolved anomaly-analysis contract

Anomaly analysis requires a time-resolved exposure product, not only a mission
average. A `TimeResolvedExposureCase` must bind:

- monotonic interval start/end times in an explicit time system, normally UTC
  for exchange plus the original source time coordinate;
- spacecraft location/trajectory and coordinate frame when the environment is
  spatially dependent;
- the spectrum and angular distribution valid in each interval;
- shield-stack and absorber-geometry identifiers valid in each interval;
- spacecraft attitude or a documented isotropic simplification;
- data-quality state, coverage fraction, gaps, and interpolation policy;
- anomaly/telemetry timestamps as a separate input that is never used to alter
  the radiation calculation silently.

For each source interval, fold the spectrum with the applicable response and
report instantaneous/interval-average proxy rates. Compute accumulated proxies
only as an explicit time integral over documented interval widths. Preserve
both the time series and the integrated value. If the shield configuration,
attitude, source product, or response changes, split the interval; do not
average unlike configurations before response folding.

No gap may be filled, cadence changed, background subtracted, or clock shifted
without a named, versioned transformation and a parallel quality flag.
Correlation, lag, superposed-epoch, threshold-exceedance, or causal analysis of
spacecraft anomalies belongs in a downstream analysis layer. The transport
model supplies traceable environmental-effect proxies and their uncertainty;
it does not label an anomaly radiation-caused.

Time-series verification must include constant-rate integration, piecewise
rate/configuration changes, irregular intervals, gaps, duplicate/out-of-order
timestamps, time-zone/time-system rejection, cadence rebinning with conserved
fluence, and equality of a time-integrated response fold with the sum of its
component intervals.

## 4. Quantity definitions and decision gates

### 4.1 TID and absorbed dose

For target mass `m`, score separately:

```text
D_abs = E_deposited,total / m
D_ion = E_deposited,ionizing / m
```

The approved UI label **TID** maps to `D_ion` in the finite physical absorber,
not automatically to total deposited energy and not to a fluence-to-dose
coefficient. Geant4 stepping must accumulate total and non-ionizing energy
deposit separately so the relationship can be audited. Report gray and rad
with an exact unit conversion. Charged-particle equilibrium, finite-volume
escape, and target-size dependence must be documented rather than corrected
silently.

### 4.2 Displacement damage dose

For particle species `i` in a specified target material:

```text
DDD = sum_i integral Phi_i(E) * S_NIEL,i,target(E) dE
```

where `Phi` is fluence and `S_NIEL` is the non-ionizing energy-loss mass
stopping power in `MeV cm2/g`, giving DDD in `MeV/g`. Bin integration must use
edges and an explicit within-bin representation; multiplying a bin-center
value by an inferred width is not acceptable when the source defines another
quantity.

Two methods are retained with deliberately different status:

1. transport-scored non-ionizing deposition in the target, as a separate
   diagnostic;
2. absorber-entry fluence folded with the approved NIEL table, as the
   authoritative DDD proxy.

Their definitions and expected differences must be documented. The existing
analytic surrogate must not be used in a production result. DDD terminology
and reporting follow ECSS-E-ST-10-12C and ASTM E3084-17(2022)e1. The data implementation
must use a version-pinned SR-NIEL or equivalent approved table that records its
displacement-threshold convention, particle/material coverage, interpolation,
and behavior outside table bounds. DDD is an exposure proxy: ASTM E3084 warns
that correspondence between NIEL and a material/device property change must
not be assumed without a demonstrated correlation. Extrapolation is forbidden
unless explicitly approved and uncertainty-bounded.

### 4.3 One-MeV-neutron-equivalent fluence

For silicon, the conventional form is:

```text
Phi_eq,1MeV,n,Si = DDD_neutron,Si / S_NIEL,1MeV,n,Si
```

The numerator is the neutron contribution only; charged-particle DDD is not
converted to `n_eq`. The NIEL numerator and reference denominator must come
from the same ASTM E722-19(2026) silicon convention. The implementation must store the reference value, edition,
and citation with the result. For SiO2, SiC, GaAs, InGaAs, Ge, water, and tissue,
the program returns “not defined for selected convention.” It must not
manufacture a material-specific `n_eq` by analogy. A later non-silicon
equivalence quantity would require a new name, definition, standard, and gate.

### 4.4 LET, lineal energy, and spectra

At least three quantities must not be conflated:

- unrestricted electronic LET/collision stopping power in a material;
- step- or track-averaged LET inside a finite sensitive volume;
- lineal energy or a detector estimator based on deposited energy divided by a
  chord/track length.

The production model must score per-step electronic energy loss and path length
inside the absorber, retaining species, charge, and primary/secondary ancestry.
The approved output is the electronic-LET distribution together with explicitly
defined track-mean and dose-mean summaries. `keV/um` is the primary linear-LET
unit; `MeV cm2/g` may also be reported only through the resolved material
density, with that conversion and uncertainty recorded. Zero-length steps,
boundary steps, stopped tracks, delta-ray energy accounting, and secondary
handling must be specified and tested. Tabulated NIST PSTAR/ASTAR values can
verify proton/helium stopping-power calculations in their covered domains.
Heavy-ion stopping power requires a separately selected reference and charge
state treatment.

CRaTER validation requires an instrument-specific deposited-energy and track
selection estimator. It is a separate validation-only output and the generic
absorber LET result must not be compared to CRaTER merely because both axes are
labeled LET.

### 4.5 `H100/10`

The approved definition is

```text
H100_10 = Phi_exit,primary-proton(E_kinetic > 100 MeV)
          / Phi_exit,primary-proton(E_kinetic > 10 MeV)
```

`Phi` is shield-rear-face plane-crossing fluence accumulated with the declared
source angular distribution and finite plane acceptance. Energy is total proton
kinetic energy, not energy per nucleon. Alpha particles, neutrons, other species,
and all secondary protons are excluded. Both numerator and denominator and
their covariance must be reported with the dimensionless ratio. The calculation
must integrate partial edge bins rather than classify bin centers; a zero or
statistically unresolved denominator returns an undefined/upper-limit state,
not zero or infinity.

### 4.6 Uncertainty

Each result must distinguish:

1. Monte Carlo statistical uncertainty;
2. incident-spectrum statistical and systematic uncertainty, including known
   covariance;
3. response-matrix statistical uncertainty;
4. material density/composition/thickness uncertainty;
5. geometry and angular-model uncertainty;
6. electromagnetic/hadronic/nuclear-model and data-library uncertainty;
7. numerical discretization, production-cut, step, and finite-area effects;
8. detector calibration/response uncertainty for observational comparisons.

Independent contributions may be combined in quadrature only when independence
is justified. Correlated source errors must retain their covariance or be
propagated through documented ensembles. Outputs with an incomplete systematic
budget must say so.

#### 4.6.1 Requested accuracy input and stopping behavior

The case file contains a dimensionless relative-accuracy target. If omitted,
the parser records the approved default `0.10` for combined standard relative
uncertainty (`k=1`, nominal 68% coverage); it must not leave an implicit
code-path default. A two-sided 95% coverage interval is additionally mandatory
in the result. A global target applies to all requested quantities unless a
quantity supplies an override. Overrides must state the same semantics or
explicitly replace target, scope, and coverage together.

“Requested accuracy” and “achieved Monte Carlo precision” are not synonyms:

- a `statistical` scope is an allowed explicit diagnostic override that asks
  the sampler to control Monte Carlo uncertainty only;
- the default `combined` scope asks whether the reported statistical plus
  quantified systematic uncertainty satisfies the target;
- neither scope proves freedom from unmodeled bias or validity outside the
  completed validation domain.

The result must repeat the requested target and scope, then report each
uncertainty component, achieved statistical relative uncertainty, combined
standard relative uncertainty (`k=1`, nominal 68% coverage), the method and
interval at 95% coverage, and one of `met`, `not_met`, or `not_assessable`.
Normality must not be assumed merely to compute the 95% interval; use propagated
ensembles or another justified coverage method when the distribution is
asymmetric or bounded. An incomplete systematic budget makes a combined request
`not_assessable`, not `met`. If a known systematic floor already exceeds the
target, the program must report that the requested combined accuracy is
unachievable with the selected model rather than increase the event count
indefinitely.

Monte Carlo execution should proceed in independently seeded batches after a
required minimum history count. It may stop only when all applicable requested
outputs meet the statistical condition at the declared confidence/coverage, or
when an input resource limit is reached. The input file must include or resolve
the minimum histories, maximum histories or wall-time budget, batch size, and
confidence/coverage convention. Reaching a resource limit returns `not_met`;
it must not silently relax the 10% request.

For a nonzero scalar `q`, the ordinary check uses its estimated relative
uncertainty `u(q)/abs(q)`. Relative uncertainty is undefined for a zero or
statistically unresolved result. Spectrum and response-matrix requests must
declare whether the target applies to each reported bin, a selected range, or
an integral; empty/rare-event bins require an absolute criterion or upper-limit
procedure. The program must not hide a poorly resolved high-LET or neutron tail
behind the precision of a large integrated quantity.

Accuracy-control verification must include:

- omitted global target resolves to exactly `0.10` in the archived case;
- global and per-quantity override precedence and serialization round trips;
- invalid nonpositive/non-finite targets are rejected;
- a synthetic estimator with known variance shows the expected
  `1/sqrt(N)` behavior and batch stopping point;
- all requested quantities, not merely the easiest one, gate stopping;
- a maximum-history/resource stop returns `not_met` reproducibly;
- a combined target below a supplied systematic floor is reported
  `not_met`/unachievable;
- zero, signed, upper-limit, correlated, and rare-bin cases never divide by
  zero or report a false success;
- rerunning with tighter targets yields statistically consistent central
  values and nonincreasing reported statistical uncertainty.

## 5. Published data placed under test now

The repository now includes a documented, checksum-pinned fixture set under
`tests/reference_spectra/` and a dependency-free test:

```bash
python3 tests/test_reference_spectra.py
```

### 5.1 Measured GCR fixture

The NASA CDAWeb product `ISS_27DAY-AVERAGES_AMS-02` supplies a calibrated
AMS-02 Bartels-rotation average beginning 2013-11-18: 45 proton and 40 helium
flux values plus total uncertainties. Units in the HAPI metadata are
`1/(m2 sr s GeV/n)`.

The Stage-2 adapter uses it for immutable-data and ingestion tests. In
addition, its explicit A=1 proton component is admitted to validation-only
cases and drives the full Geant4/proxy path. The HAPI
coordinate objects are named `BinProLow`/`BinHeLow` but expose an array called
`centers`; the pinned AMS supplemental Table SM XXXV establishes that they are
transformed rigidity-bin lower edges. The adapter preserves independent grids,
converts rigidity edges to kinetic energy per nucleon, applies the differential
unit conversion, and retains total uncertainties. It is not a complete GCR
test environment because it lacks heavier ions, represents a specific near-
Earth interval, has limited energy coverage, and publishes inclusive helium
rather than a pure alpha population. The direct transport test uses only the
proton component: for A=1, kinetic energy per nucleon equals total proton
kinetic energy. It samples the published differential intensity as a
piecewise-constant density over the exact transformed bin edges, applies the
isotropic plane-crossing factor, and records the finite source-plane incident
rate. This narrow mapping is not generalized to inclusive helium.

The observational source-validation suite now also pins the AMS-02 2025
solar-cycle He, Li, Be, B, C, N, and O tables. They contain 147 available
Bartels rotations between May 2011 and November 2022 on 40 bins from 1.92 to
60.3 GV. `AMS02SolarCycleAdapter` keeps the independent variable as rigidity
and retains statistical, time-dependent systematic, and time-independent
systematic errors separately. It does not assign an isotope/mass number or
convert to energy per nucleon: the tables publish inclusive elements, and an
isotope mixture has not been supplied.

The optional AMS partition record freezes the first 118 published rotations through
BR2552 for development and the final 29 rotations beginning BR2553 as a
chronological holdout for a future source-environment comparison. It is not a
Stage-2 transport requirement and is not activated by the current direct-use
test, which consumes the separate BR2460 proton fixture. Missing rotations are
never interpolated. Per-bin quadrature error bars may be drawn for review, but
the time-independent term is not treated as independent between rotations.

ACE/CRIS now supplies a second pinned near-1-AU observation: the 2019-12-01
daily `AC_H3_CRIS` HAPI record contains seven intensity bins and seven counts
for each element B through Ni. `ACECRISHapiAdapter` takes exact energy bounds
from the label variables embedded in the pinned CDAWeb V04 CDF. For positive-
count bins it calculates only the analytic counting-statistic term
`flux/sqrt(count)`; zero-count bins are retained and explicitly invalid for
quantitative comparison.

The current ACE factors file and release notes are pinned as documentation,
but their revised energy bounds differ from the V04 CDF labels. The adapter
therefore prohibits combining those factors with this fixture. It also leaves
systematic uncertainty unassigned because the applicable processing-version
mapping has not been established. The single ACE day is adapter evidence, not
a model holdout; a defensible time partition requires a declared comparison
objective and an approved SEP/instrument-quality exclusion rule.

Production GCR cases will use the future AMPS adapter or another explicitly
selected, mission-applicable environment with the required species, epoch,
location, and angular boundary. The downloaded observational fixtures are
test assets, not production defaults. Existing repository files carrying a
BON2020 name remain prohibited because their generation, version, parameters,
species interpretation, and checksums have not been recovered; no BON2020
package or adapter is required to complete Stage 2.

### 5.2 Measured SEP fixture

The NOAA GOES-16 SGPS fixtures include the published September 2017 five-minute
record beginning 2017-09-11 07:30 UTC and the full selected GLE72 interval.
They contain both directional sensors, 13 proton channel bounds from about 1
to 404 MeV, uncertainties, sample counts, quality flags, and yaw state. The
original netCDF flux units were
`protons/(cm2 sr keV s)`; the extraction converts them to per MeV and records
the transformation.

The Stage-2 adapter retains the NOAA interval start, derives midpoint and end,
and produces separate time-resolved inputs and coverage-aware partial fluences
for the primary `>100 MeV` and sensitivity `>10 MeV` event windows. It must not
be fed to an isotropic source. Transport still requires a justified mapping
from geosynchronous directional measurements to the simulation boundary.

### 5.3 CRaTER observational fixture

The repository contains CRaTER GCR LET spectra version 1.3 from Zenodo. The
product contains 34 quiet-time, approximately half-year spectra from 2009 to
2026 and is described as measurements through the CRaTER telescope behind
9.9 g/cm2 of mostly tissue-equivalent plastic.

This file is a future observational benchmark, not a current generic-slab
expected answer. The published methods show that instrument geometry,
coincidence/track selections, acceptance, lunar viewing, and the deposited
energy-to-path estimator are integral to the product. The spectral ordinate's
exact normalization and uncertainty representation must be recovered from the
record/method before a statistical comparison is designed; the CSV header alone
does not define them.

## 6. Verification and validation hierarchy

The test hierarchy is cumulative. A higher layer cannot waive a failure below
it.

### V0: data, provenance, and schema tests

These tests require no Geant4 and run on every change:

- SHA-256 matches each reviewed fixture;
- source URL/DOI, version, license, retrieval record, and original units exist;
- arrays are rectangular and match published/declared dimensions;
- bins are positive and ordered when the source defines their semantics;
- values and uncertainties are finite and valid/fill values are rejected;
- directional and temporal metadata cannot be dropped by serialization;
- unit conversions round-trip and are dimensionally correct;
- ambiguous metadata causes an explicit failure, not a warning-and-guess.

`tests/test_reference_spectra.py` implements these checks for both AMS-02
products, ACE/CRIS, GOES-16 SGPS, and CRaTER.

### V1: exact software and analytic tests

These tests isolate algorithms from stochastic transport:

1. `g/cm2 <-> cm` conversion for known density, including dimensional-error
   rejection.
2. Exact conversion among `m^-2` and `cm^-2`, GeV/n and MeV/n, and total energy
   and energy/nucleon using explicit `A`.
3. Closed-form integration of constant and power-law spectra over known bin
   edges, including partial 10 and 100 MeV edge bins.
4. Isotropic plane crossing: `p(mu)=2mu`, CDF `mu^2`, moments
   `E[mu^k]=2/(k+2)`, azimuthal uniformity, and the `pi` angular factor.
5. Non-interacting slab path length `t/mu`, subject to explicitly tested finite
   transverse geometry.
6. Dose conversion from a prescribed deposited-energy ledger and target mass:
   joule/kg, Gy, and rad.
7. Separate total, ionizing, and non-ionizing ledgers, including
   `E_ion = E_total_deposit - E_nonionizing_deposit` if that is the approved
   Geant4 scoring convention.
8. Response-matrix identity, zero response, linearity, superposition, and exact
   convolution of a small hand-calculated matrix.
9. Conservation of integrated fluence when rebinned, with explicit treatment
   of partial overlap.
10. LET histogram conservation: the integral over LET bins equals the accepted
    charged-particle fluence for in-range entries.
11. Analytic uncertainty propagation for diagonal covariance and a hand-built
    correlated covariance matrix.
12. `n_eq` consistency when DDD and the 1-MeV neutron reference are prescribed;
    unsupported material/convention pairs must fail.
13. Exact multilayer placement/mass accounting, order-preserving round trips,
    one-layer compatibility, and same-material split-layer equivalence.
14. Slab-absorber analytic volume, mass, intercept, and chord/path calculations,
    plus deterministic geometry serialization and hashing.
15. Accuracy-default/override resolution, adaptive-batch stopping, resource-limit
    failure, systematic-floor rejection, and zero/rare-bin handling from
    Section 4.6.1.

Pure deterministic calculations should agree to a tolerance derived from
floating-point conditioning. Stochastic source-distribution tests should use a
predeclared false-rejection rate and sample size, not an arbitrary percent
band.

### V2: authoritative tabulated-reference tests

Use independent reference data rather than values generated by shieldSim:

- NIST PSTAR proton electronic/nuclear stopping power and CSDA range;
- NIST ASTAR helium stopping power and range;
- NIST ESTAR when electron transport becomes a reported source component;
- NIST XCOM photon mass attenuation and the separately published NIST X-Ray
  Mass Attenuation tables' mass energy-absorption coefficients for narrow-beam
  electromagnetic checks. XCOM attenuation must not be described as an
  absorbed-dose calculation;
- approved NIEL/SR-NIEL tables and the selected 1-MeV-neutron-equivalent
  standard;
- versioned Geant4 material compositions/densities for `G4_*` entries;
- project/supplier/laboratory data for custom composites and regoliths.

Reference tables must be pinned with license/provenance and compared only over
their stated domains. The acceptance rule for each table must be registered
before the results are examined, using the table uncertainty, interpolation
error, and Monte Carlo uncertainty. NIST notes, for example, that PSTAR/ASTAR
uncertainty depends strongly on material and energy; one universal percentage
is not defensible.

Suggested controlled configurations include:

- proton/helium range in homogeneous Al, water, and Si;
- thin-target electronic stopping power where energy loss is small;
- photon narrow-beam attenuation `I/I0 = exp[-(mu/rho)x]` with scattering
  rejection matching the XCOM definition;
- CSDA energy-loss curves with straggling reported separately;
- selected silicon NIEL points for every supported incident species.

### V3: Geant4 transport and numerical-convergence tests

For monoenergetic inputs and every material family:

- compare independent absorbers in reversed selection order; outputs must be
  statistically indistinguishable and geometry hashes identical per case;
- compare one-layer and same-material split-layer stacks, then exercise
  high-Z/hydrogenous interfaces in both physical orders without asserting an
  unreferenced expected magnitude;
- verify local boundary scoring in rotated/translated test geometry;
- expand transverse size/source plane until the reported quantity changes by
  less than its approved numerical-error allocation;
- vary production cuts and maximum step length over at least three resolved
  values and demonstrate a plateau;
- compare `FTFP_BERT`, an HP neutron list, and an approved ion configuration;
- check primary plus secondary particle/energy accounting with a ledger whose
  scope is explicitly defined for nuclear rest-mass/Q-value effects;
- verify mass-thickness equivalence only in cases for which geometry and
  material composition make it physically applicable;
- demonstrate statistical `1/sqrt(N)` behavior with independent seeds;
- show response-matrix and direct spectrum-sampling agreement;
- test extreme thin/thick cases, stopped primaries, backscatter, lateral
  escape, and zero-transmission bins.

Convergence is quantity-specific. A setting accepted for transmitted protons
does not automatically qualify TID, neutron DDD, or the high-LET tail.

### V4: published-spectrum end-to-end tests

After adapters and definitions pass V0-V3:

- ingest an approved full GCR scenario and the GOES-16 reference record without
  losing metadata;
- preserve and report source components through transport and folding;
- compare direct-sampling and response-matrix results;
- propagate spectrum uncertainties through at least linear covariance and
  sampled-ensemble paths;
- demonstrate stable results against input bin refinement/coarsening where the
  source resolution supports it;
- assert that an attempt to label the two-direction GOES snapshot “isotropic”
  without an approved transform fails;
- assert that incomplete AMS-02 proton/helium coverage cannot be labeled a
  complete GCR environment;
- retain the current simple monoenergetic and synthetic spectra only as
  software regressions.

### V5: observational validation

Observational validation must reproduce the measurement process, not compare
dissimilar quantities because their labels match.

#### V5a: CRaTER GCR LET benchmark

Build a versioned CRaTER instrument case from primary instrument/method
documentation. This is a validation harness for the environment-effect physics,
not a generally supported instrument-design workflow:

1. silicon detector and tissue-equivalent-plastic dimensions, composition,
   density, dead layers, housing/spacecraft overburden, alignment, and
   separations;
2. telescope angular acceptance and lunar-pointing/viewing selection;
3. the exact detector coincidence, energy-deposit, track/chord, quality, and
   quiet-time selection used by the published product;
4. electronics threshold, saturation, resolution, calibration, dead time, and
   efficiency required by the product;
5. lunar solid-angle obstruction, surface albedo/secondaries, and orbital
   geometry appropriate to the interval;
6. a time-matched incident GCR model including the species needed by the LET
   range;
7. the same reported LET/lineal-energy estimator and binning.

Use at least one interval for implementation/calibration and reserve later
intervals as blind tests. Compare absolute normalization and shape with a
likelihood or chi-square that uses the observational covariance when available;
otherwise report the missing correlations as a limitation. Also compare
integrated rates over preregistered, fully efficient LET intervals. Acceptance
limits must be set before unblinding from combined observational, source,
transport, and Monte Carlo uncertainty.

The 2013 CRaTER/Geant4 study, *The radiation environment near the lunar
surface: CRaTER observations and Geant4 simulations*, is a useful reproducibility
target because it explicitly combines the GCR environment, lunar albedo,
shielding, and sensor response. Reproduce one published configuration before
claiming independent agreement.

#### V5b: absorbed-dose cross-checks

A 2024 lunar-surface study reports a December 2009 CRaTER D1/D2 silicon dose
comparison of approximately `9.86 +/- 0.02 uGy/h` after its stated scaling,
versus its model value of about `9.3 +/- 0.1 uGy/h`. This is a candidate
cross-check, not a ready expected value: the comparison depends on detector
geometry/material, viewing, lunar shielding/albedo, neutral-particle response,
time selection, and the study's conversion. The primary CRaTER calibrated data
and exact analysis procedure must be recovered before it becomes an acceptance
test. The secondary table alone must not be used as ground truth.

LND/RADOM or other lunar measurements may provide complementary checks, but
each requires its own instrument-response model. Agreement with one detector
does not validate all absorber materials or quantities.

### V6: regression tests

Only after V1-V5 results have been reviewed should their compact outputs become
regression baselines. Store:

- approved mean and statistical uncertainty, not only one random realization;
- permitted numerical/physics drift and its rationale;
- random seed ensemble and effective event count;
- exact code, Geant4, data, and compiler versions;
- validation domain and known exclusions.

Regression tests answer “did validated behavior change?” They must never be
used to answer “was this behavior correct?” Baselines are updated only through
a review record explaining the physical or numerical reason.

## 7. Development stages and gates

### Stage 0 - requirements and validation protocol freeze

**Status: complete as a technical requirements freeze on 2026-10-16.** This
does not validate the current code or enable any science output. The project
owner approved the decisions in Section 1.1. Stage 2 is explicitly exempt from
independent review; the later physical capability and observational-validation
gates retain their stated review requirements.

Entry: this roadmap and the screenshots.  
Work:

- record the resolved decisions in Section 1.1 and the disabled/deferred inputs
  in Section 10;
- define project-owner acceptance for Stage 2 and independent-reviewer roles
  for the later physical transport, dosimetry, device-damage, and observational
  validation gates;
- write versioned definitions for every enabled output;
- record the confirmed anomaly-analysis/preliminary-assessment purpose, then
  declare the numerical environment domain, accuracy goals, and excluded uses;
- select the pinned CRaTER product as the observational target and require its
  exact calibration/blind interval split to be preregistered before Stage 8
  implementation sees the blind values;
- register acceptance rules before model/data comparisons;
- create a requirement-to-test traceability file.

Tests: V0 fixture tests and documentation/checksum tests. The committed
reference-fixture test passed five checks on 2026-10-16.  
Exit evidence: this approved decision record, output definitions, material
library contract/disposition, selected standards, pinned source/observational
fixtures, traceability table, and unresolved inputs reduced to explicitly
disabled capabilities. Named independent reviewers and their signed findings
are evidence only for later gates that retain that requirement, not for Stage 2
and not as a claim attached to Stage 0.  
Forbidden carry-over: no output with an unresolved definition may be called
validated.

### Stage 1 - correctness repairs and strongly defined core data

**Status: implemented and acceptance-tested on 2026-10-16.** This status means
the Stage-1 software/data boundary is complete; it does not promote TID or
H100/10 to a physically validated capability and it does not advance any
Stage-2-or-later source/response capability.

Implementation evidence:

- `CaseDefinition`, `SpectrumDataset`, `MaterialLibrary`, `MaterialFactory`,
  `AccuracyController`, strict JSON, canonical SHA-256 hashing,
  `AbsorberGeometry`, and `SlabAbsorberGeometry` are implemented in
  `include/` and `src/`;
- `data/material_library.v1.json` contains the approved explicit material
  dispositions. Generic CFRP, SiC composite, lunar regolith, and Mars regolith
  are disabled templates; SiO2 and Ga0.47In0.53As contain the approved density
  corrections. Case records may extend an exact parent or define a fully
  explicit new material, including explicit isotope recipes;
- `examples/stage1_case.json` is a complete two-layer/two-absorber example. It
  explicitly supplies world, source-plane, shield, and absorber size and
  placement fields and deliberately omits the global target to demonstrate
  resolved `0.10` default provenance;
- each absorber is constructed and scored in a separately rebuilt Geant4
  world. Stable absorber IDs control file naming, and the random seed/tallies
  are reset at the independent-case boundary;
- step scoring now closes separate total, ionizing, and non-ionizing deposited-
  energy ledgers. H100/10 uses primary protons and exact strict-threshold rear-
  face crossing energies;
- cumulative event-level variance is evaluated between bounded BeamOn batches.
  Complete/incomplete systematic budgets and known systematic floors produce
  explicit `met`, `not_met`, or `not_assessable` states without silently
  relaxing the requested target;
- each absorber writes a canonical resolved-case manifest and
  `shieldSim.result.v1` result. These carry material/source/geometry hashes,
  build/source/Geant4 provenance, seed, definitions, units, normalization,
  ledger/count evidence, uncertainty state, and the explicit
  `implemented_not_physically_validated_stage1` capability label; and
- legacy three-column spectra now require a versioned metadata sidecar; rows
  are never silently skipped, sorted, made positive, or assigned an undeclared
  single-bin width. Unknown CLI options and built-in GCR out-of-domain requests
  are fatal.

Acceptance evidence is registered in CTest:

```text
stage1_core
stage1_case_boundaries
stage1_transport_smoke
```

The tests cover canonical round trips and hashing, fraction/constituent/isotope
resolution, immutable overrides, template gates, G4 density and analytic slab-
mass checks, one/five/six layers, 100.000/100.001 mm, missing dimensions and
units, layer order, single-vs-split analytic layer ledgers, accuracy
default/override/`1/sqrt(N)`/systematic-floor behavior, independent-case order
invariance, and deposited-energy ledger closure. The pre-existing Layer-2/3/4
legacy spectrum calls were also changed to supply explicit metadata; those
scripts remain compatibility regression tests and are not substituted for the
versioned acceptance suite.

Entry: Stage 0 definitions needed for core semantics.  
Work:

- implement `SpectrumDataset`, `CaseDefinition`, `ShieldStack`,
  `SlabAbsorberGeometry`, the geometry-neutral `AbsorberGeometry` contract,
  typed unit conversions, and result provenance;
- implement the external versioned material library, per-case immutable
  material overrides, status gating, provenance, canonical expansion/hashing,
  and the material V0/V1 tests in Section 3.3.5;
- replace/relabel the current ambiguous SiO2 definition, correct the explicit
  `Ga0.47In0.53As` reference density, and migrate generic CFRP, SiC composite,
  and regolith entries only as disabled templates or clearly named trade-study
  records;
- make each shield x absorber case independent;
- make ordered multilayer stacks first-class inputs, with mass thickness the
  primary screenshot-facing value per layer while preserving an explicit
  geometric-thickness mode;
- enforce one through five layers and the inclusive 10 cm derived geometric
  thickness ceiling on every layer;
- require every shield, source-plane, world, spacing, placement, and slab
  dimension in the versioned input file;
- implement the recorded 10% default accuracy, per-quantity overrides, declared
  scope/coverage, adaptive batch evaluation, and `met/not_met/not_assessable`
  result state;
- separate total, ionizing, and non-ionizing deposited energy;
- reject unknown CLI fields, missing units, invalid ranges, NaN/fill values,
  and unsupported quantity/material combinations;
- preserve the legacy parser only behind explicit metadata declarations;
- fix any built-in energy-range clamping/reordering behavior so invalid ranges
  fail instead of silently changing the request.

Tests: V1 unit, schema, dose-ledger, independent-case, layer-order,
one/five/six-layer boundaries, 10 cm/over-limit boundaries, required-dimension
failures, accuracy default/override/stopping behavior, single-vs-split-layer
equivalence, geometry serialization/hash, and failure tests; all existing
Layer-1/2 tests revised to the new contracts.  
Exit evidence: zero order dependence, complete serialized provenance, exact
unit tests, and no silent defaulting at a scientific boundary.  
Forbidden carry-over: current `/m2` files cannot enter the `/cm2` legacy path.

### Stage 2 - published-source adapters

**Status: complete and acceptance-tested on 2026-10-16.** The two AMS-02
adapters, ACE/CRIS adapter, GOES-16 adapter, pinned source bytes,
transformations, uncertainty boundaries, event windows, and tests are
implemented. A downloaded AMS-02 proton spectrum also drives the full
production case/Geant4/proxy path. BON2020 is not part of production operation
and is not a Stage-2 completion dependency. No independent reviewer is
required for this stage.

Implemented evidence:

- `AMS02HapiAdapter` verifies the HAPI metadata/data and primary supplement
  hashes, resolves the misleading generic `centers` label with Table SM XXXV,
  transforms rigidity-bin edges to kinetic energy per nucleon, converts
  differential units, retains uncertainties, and keeps inclusive helium from
  being silently labeled as alpha particles;
- `AMS02SolarCycleAdapter` verifies all seven primary CSV hashes and layouts,
  common rotation/grid coverage, separate uncertainty components, explicit
  rotation selection, and missing-rotation rejection while preserving the
  inclusive elements on their published rigidity coordinate;
- `ACECRISHapiAdapter` verifies HAPI/CDF/documentation hashes, B--Ni coverage,
  exact V04 CDF-label energy bounds, daily values/counts, count-derived
  statistical errors, and zero-count masks. It rejects use of the nonidentical
  current factor-table bounds and does not invent a systematic error;
- `GOESRSGPSAdapter` verifies the event definition and extracted-data hash,
  retains every interval start/midpoint/end, both directions and all channels,
  enforces the approved quality predicate, and closes separate valid/rejected/
  missing-time ledgers for both event windows;
- the primary `>100 MeV` window is `[2017-09-10T16:25:00Z,
  2017-09-12T22:40:00Z)` and the sensitivity `>10 MeV` window is
  `[2017-09-10T16:45:00Z, 2017-09-14T18:50:00Z)`. The NMDB
  `2017-09-10T16:08:00Z` onset is retained as an independent marker only;
- `shieldSimSourceTool` emits deterministic canonical JSON without invoking
  transport; and
- `ams02_hapi_component` is a strict validation-only case source. It admits
  only the AMS A=1 proton component, requires isotropic source mode, samples
  exact transformed energy-bin edges, and rejects mission-purpose cases,
  inclusive helium, unsupported units, invalid masks, and non-overlapping
  energy limits;
- `stage2_published_transport` independently reconstructs the published
  proton source integral, verifies the `pi` plane-crossing factor and finite
  source-plane rate, transports 256 proton histories through the normal
  multilayer Geant4 geometry, and verifies TID and `H100/10` production
  outputs and source/hash provenance; and
- the official NOAA netCDF, complete extracted event series, AMS primary
  supplement, transformation records, and all compact fixtures are checksum-
  pinned in `tests/reference_spectra/manifest.sha256`.

The registered CTest evidence is `stage2_reference_spectra`,
`stage2_source_adapters`, `stage2_source_tool`,
`stage2_published_transport`, and the dependent
`validation_figures` review-artifact generator. The latter archives PNG/SVG
views and a checksum/caption manifest for every current analytic or
observational comparison; numerical assertions remain authoritative.
The standalone Layer-2 and Layer-3 runners likewise generate PNG/SVG and
manifest review artifacts for their analytic source, geometry, normalization,
and threshold-limit checks. Their controlled LET plot is explicitly marked as
a legacy sanity comparison rather than physical validation.

The retired BON2020 intake directory records why its unverified legacy files
remain prohibited; it is not an active completion gate.

Entry: validated canonical schema and unit system.  
Work:

- maintain the implemented lossless readers for the pinned AMS-02, ACE/CRIS,
  and GOES-16 fixtures and their strict canonical boundaries;
- retain the resolved AMS-02 energy-coordinate semantics and primary-table
  Jacobian evidence;
- preserve GOES sensor direction, time support, and quality metadata;
- keep the downloaded AMS proton path wired through the same production case,
  source sampling, Geant4 transport, and proxy writers used by later AMPS
  inputs; do not introduce a test-only transport implementation;
- sample documented bins using their exact edges and independently verify the
  integrated intensity, angular factor, plane area, and incident rate;
- maintain the implemented GOES-16 September 10-11, 2017 GLE directional time
  series and its provenance-preserving event-fluence integration. Keep both sensor
  directions separate, do not interpolate missing/rejected samples, and report
  per-channel coverage and incomplete status. Record the published event onset
  as the inclusive UTC start and the published event end as the exclusive UTC
  stop. Preserve every NOAA interval-start coordinate, derive the midpoint and
  end, and duration-weight only the exact overlap with the event window. Reject
  transport until the sensor
  look-direction mapping to the boundary is explicit and tested;
- deprecate the untraceable built-in GCR approximation for physical runs.

Tests: V0 plus round-trip serialization, analytic unit/Jacobian tests, rejected
ambiguity and contract-mutation tests, published-row spot checks, exact
start/midpoint/end conversion, constant/linear/partial-boundary integration,
independent quality-predicate mutations, and missing/flagged-bin coverage
tests. The source tests also cover AMS uncertainty-component retention and
chronological partition boundaries, ACE `flux/sqrt(count)` statistics,
zero-count masking, the V04/current-factor-table mismatch, and direct
downloaded-spectrum Geant4/proxy execution. The end-to-end test also rejects
inclusive-helium transport, beam relabeling, and mission-purpose use of the
near-Earth fixture. These implemented tests pass.
Exit evidence: byte-pinned raw inputs, transformation manifests, independently
calculated spot checks, an end-to-end Geant4/proxy result with source hashes,
explicit test-only labels, and project-owner acceptance. Stage 2 has no
independent-reviewer sign-off requirement.
Forbidden carry-over: the AMS proton fixture is a transport test input, not a
complete deep-space environment or production default. AMS inclusive helium,
AMS He--O, and ACE B--Ni remain ingestion-only because isotope composition is
not defined. GOES remains ingestion/integration-only until its directional
boundary mapping is approved. None may be relabeled as another environment
class.

### Stage 3 - independent geometry and response engine

**Status: implemented and acceptance-tested on 2026-10-16; pending independent
review and therefore not formally complete.** This status applies only to
geometry/transport isolation, raw response infrastructure, ancestry, and raw
ledgers. It does not enable a science-facing radiation proxy.

Implementation evidence:

- `shieldSim.case.v2` adds explicit source/shield/absorber local-to-world
  transforms. Unit quaternions are validated and canonicalized; oriented-box
  containment and overlap checks execute before Geant4. Geant4 placements use
  its direct-object transform API so the source and solids share the same
  rotation convention;
- `SlabAbsorberGeometry` implements validation, stable sensitive-region IDs,
  analytic mass properties, bounds/reference frame, canonical serialization,
  content hashing, and construction. Scorers depend on the returned logical
  regions rather than slab copy numbers;
- all ordered shield layers are separate Geant4 solids. Raw tallies retain
  layer ID, total/ionizing/nonionizing deposition, track length, step count,
  secondary particle identity, creator process, root primary, and ancestry
  generation. The energy-ledger output explicitly excludes a global nuclear
  rest-mass/Q-value conservation claim;
- monoenergetic response nodes accept proton, He-4, or an explicit `(Z,A)` ion,
  total or per-nucleon kinetic energy, and exact `(mu,phi)`. The independent
  `discrete_response_nodes` source samples an exact normalized node mixture;
- each response cell stores event sums, squared sums, means, sample variances,
  standard errors, and the complete primitive-score covariance upper triangle,
  plus normalization, incident coordinates, and a response key. Folding
  propagates both response-matrix and source-weight cross-covariance so later
  numerator/denominator algorithms need not discard correlation.
  That key covers geometry/provider/material/scoring, Geant4/code/build/data,
  physics/ion, cut, and step state and deliberately excludes seed, event count,
  output path, and incident grid coordinate;
- `shieldSimResponseTool` strictly assembles compatible cells into canonical
  `shieldSim.response_matrix.v1` JSON plus checksum-recorded CSV and folds
  explicit weights into primitive scores. It rejects key mismatch, key-hash
  mismatch, duplicate/unknown nodes, invalid normalization, unsupported
  interpolation, extrapolation, non-finite/asymmetric/non-positive-semidefinite
  source covariance, and score-unit mismatch;
- `scripts/build_response_weights.py` writes reviewable explicit weights for
  exact nodes or response interpolation linear in log energy and bilinear in
  `mu`/periodic `phi`, rejecting out-of-domain and incomplete tensor grids; and
- `stage3_core`, `stage3_response_weights`, `stage3_transport`, and
  `validation_figures` are registered in CTest. All generated Stage-3 transport artifacts are kept under
  `test_output/stage3_transport_runs` (the corresponding CMake build tree), not
  `/tmp`.

Acceptance evidence:

- deterministic tests cover transform inversion and `q/-q` identity, exact
  layer interfaces/mass thickness, A|B versus B|A hashes, slab volume/mass/path,
  response zero/identity/linearity/superposition, analytic covariance, strict
  response boundaries, and exact zero handling;
- the rotated/translated geometry, alternative-absorber permutations, and
  local rear-face scorer pass; single versus split Al satisfies the five-seed
  95% equivalence gate;
- independent direct sampling and matrix folding satisfy the same gate for
  both W-to-HDPE and HDPE-to-W stacks without asserting an unreferenced expected
  physical ordering; both layers retain secondary-production diagnostics;
- production-cut, maximum-step, common transverse-size, source-plane-size,
  response-grid, and event-count drivers pass their registered criteria;
- proton, He-4, Fe-56, `FTFP_BERT`, `FTFP_BERT_HP`, reference-list ions,
  `G4IonQMDPhysics`, and representative Al/W/HDPE/Si/water configurations run
  with explicit provenance; and
- a thick-W stopped-primary case supplies an exact zero-transmission boundary.
  PNG/SVG review figures and a hash/caption manifest are downstream of the
  numerical assertions.

Entry: independent case model and approved shield/absorber definitions.  
Work:

- build one shield/absorber world per case or demonstrate equivalent isolated
  scoring with no cross-case transport;
- build every ordered shield layer explicitly, validate interfaces/overlaps,
  and retain layer-resolved secondary-production diagnostics;
- construct current absorbers only through `SlabAbsorberGeometry`, with scoring
  attached to stable logical region IDs supplied by the geometry interface;
- implement versioned monoenergetic/species/angular response matrices;
- retain particle ancestry and primary/secondary component scoring;
- implement direct-spectrum mode as an independent cross-check;
- add finite-size, source-plane, cut, step, and event convergence drivers;
- version physics and nuclear-data configuration in response keys.

Tests: V1 matrix tests, all multilayer tests in Section 3.3.1, slab analytic
mass/path tests, and V3 monoenergetic/convergence/direct-fold tests.  
Exit evidence: statistical equivalence of absorber-case permutations,
single-vs-split-layer equivalence, exact layer/interface metadata,
direct-vs-folded results, and convergence envelopes for each quantity class.  
Formal sign-off additionally requires the independent-review record named in
Section 10; until it exists, the correct stage label is “implemented and
acceptance-tested, pending independent review,” not “complete.”
Forbidden carry-over: no response matrix is reused when any hashed geometry,
material, physics, scoring, or Geant4 field changes.

### Stage 4 - physically defined TID

**Status: implemented and acceptance-tested on 2026-10-17; pending independent
review and therefore not formally complete.** The accepted domain is proton
and helium-4 finite material-slab TID in aluminum, silicon, and water. It does
not include electron or heavier-ion incidence, gamma TID, arbitrary shapes,
device or organ response, or a crew-health quantity.

Implementation evidence:

- case v2 now admits `TID` alongside `transport_raw`; DDD, `n_eq`, LET, and
  `H100/10` remain rejected. A gamma `(Z=0,A=0)` total-energy node is admitted
  only for the photon attenuation harness. The parser rejects `TID` whenever a
  monoenergetic or discrete response source contains gamma or an explicit ion
  heavier than He-4. It also rejects absorbers outside the registered
  aluminum/silicon/water material identities; shielding-layer materials remain
  input-defined and are not conflated with the absorber-response domain;
- TID is the event-summed Geant4 ionizing-deposition ledger divided by the
  exact mass of the active finite slab. Total and nonionizing deposition remain
  separate and close by `total = ionizing + nonionizing`. Results label this as
  material-slab dose and write Gy/primary, the numeric rad value (100 times the
  numeric Gy value), and source-normalized rate without claiming device,
  organ, or tissue dose;
- every depositing step is classified authoritatively by transported species,
  primary/secondary identity, root primary, and ancestry generation. Separate
  creator-process and step-defining-process partitions are explicitly marked
  diagnostic because a process name is not a unique causal allocation of
  deposited energy. Each independent view closes to the aggregate ledger;
- `transport.em_physics` is mandatory in case v2 and participates in resolved
  case and response hashes. `reference_list_default` is the production
  baseline and `G4EmStandardPhysics_option4` is the registered sensitivity;
- a complete case-v2 TID systematic budget contains named material/mass,
  source, numerical, and physics-reference standard-uncertainty components,
  provenance, and every unordered correlation pair. The parser rejects
  missing/duplicate pairs, invalid coefficients, absent categories, legacy
  aggregate-only claims, and non-positive-semidefinite correlation matrices.
  It computes the correlated combined standard uncertainty. The explicitly
  selected coverage model combines a finite-sample Student-t MC interval with
  normal expanded systematic components in quadrature and writes both MC-only
  and combined 95% intervals. Missing budgets remain `not_assessable`;
- `shieldSimReferenceTool` exposes deterministic Geant4 electronic/total
  stopping, numerical CSDA integration, and gamma mass attenuation in the NIST
  mass-normalized units. CSDA is the log-energy quadrature of Geant4's
  unrestricted total stopping model from the 1 keV STAR lower boundary; it is
  not a stochastic slab-range estimator; and
- `scripts/download_nist_stage4.py` retrieves official PSTAR/ASTAR, XCOM, and
  NIST X-Ray Mass Attenuation pages. Tests use only the committed raw bytes,
  canonical CSVs, and SHA-256 manifest under `tests/reference_tid/nist/`; they
  do not depend on a live service. The X-ray tables retain both `mu/rho` and
  `mu_en/rho`, but Stage 4 compares `mu/rho` only and never substitutes either
  coefficient for finite-slab dose.

Acceptance evidence:

- `stage4_core` passes exact MeV/kg-to-Gy, numeric Gy-to-rad, and hand-
  calculated correlated-uncertainty comparisons. It also passes all
  capability/rejection boundaries—including gamma-node admission for raw
  attenuation checks, gamma/heavy-ion TID rejection, and rejection of an
  absorber outside the Al/Si/water matrix—and combined-interval checks.
  Analytic slab-mass checks reside in the full transport gate where Geant4's
  constructed mass is available;
- `stage4_nist_reference` authenticates all 24 pinned raw/canonical files and
  passes 63 exact-node comparisons across the three materials. The maximum
  reference relative difference is 0.8286%; the maximum default-vs-option-4
  deterministic table sensitivity is 1.1409%. Each comparison uses the
  row-specific registered limit in `stage4_reference_plan.v1.json` and has zero
  interpolation uncertainty at the selected node;
- `stage4_tid_transport` executes the full case/parser/generator/multilayer/
  Geant4/scorer/accuracy/result path for both particles and all three absorber
  materials. The maximum material/species MC relative standard uncertainty is
  0.8480%. All 44 exact mass, dose, and partition-closure comparisons pass;
- final production-cut, maximum-step, and lateral-escape refinements change
  their accepted observables by 1.4593%, 1.1620%, and 0.2868%, respectively.
  Cut and step use TID at fixed geometry. Lateral size uses deposited energy
  because changing slab area changes its mass and therefore Gy per one compact-
  plane primary by definition; exact TID conversion is still checked at every
  size;
- option-4 finite-slab TID differs from the production baseline by 0.8771%,
  inside its registered 4% physics-reference allocation; and
- the median high/low-history relative-standard-error ratio is 0.5872 for a
  fourfold history increase, inside the registered `[0.35,0.65]` interval
  around the `1/sqrt(N)=0.5` expectation. Low-history samples are convergence
  diagnostics, not accepted production cases.

The CTest evidence is written under `test_output/validation_evidence`, all
transport inputs/logs/results are under `test_output/stage4_tid_transport_runs`,
and three Stage-4 PNG/SVG subjects plus manifest hashes are under
`test_output/validation_figures`. Numeric JSON evidence, not the figures, is
the pass/fail authority.

Entry: approved TID observable and absorber geometry.  
Work:

- score ionizing/non-ionizing deposition separately;
- implement uncertainty accumulation and particle/process breakdown;
- validate target mass and finite-volume scoring;
- label material-slab dose separately from organ/device dose.

Tests: V1 exact energy/mass conversions, V2 PSTAR/ASTAR/XCOM controlled tests,
and V3 cut/step/size/statistical convergence.  
Exit evidence: all TID verification cases meet preregistered combined
uncertainty criteria.  
Formal sign-off additionally requires the Section-10 independent-review
record. Until it exists the correct label is “implemented and acceptance-
tested, pending independent review,” not “complete.”
Forbidden carry-over: passing TID does not enable DDD, `n_eq`, LET, or a crew
health claim.

### Stage 5 - physically defined DDD and `n_eq`

**Status: implemented and reference-tested on 2026-10-17; pending a project-
selected accelerator/device-degradation validation data set and independent
review, and therefore not formally complete.** The implemented domain is
natural-silicon DDD for proton, helium-4, and neutron absorber-entry fluence and
ASTM E722-19 neutron-only `n_eq`. It excludes electrons, ions heavier than
helium-4, non-silicon targets, neutron energies above 20 MeV, arbitrary absorber
shapes, and every claim of device degradation.

Implementation evidence:

- `DamageResponseLibrary` loads
  `data/niel/damage_response_library.v1.json`, verifies every canonical table
  SHA-256, and strictly validates identity, version, target, units, total-
  kinetic-energy coordinates, source provenance, positive/monotonic nodes,
  contiguous neutron groups, group hardness ratios, and the independent
  conversion of 95 MeV-mb to `0.002037006185398159 MeV cm2/g` using the
  manifest atomic weight and Avogadro constant;
- the package pins SR-NIEL-7 v11.1 proton and helium-4 responses for natural Si,
  `Ed = 21 eV`, Coulomb plus hadronic terms, and no hadronic scaling. Its
  neutron file pins 640 ASTM E722-19 groups with explicit boundaries from the
  SAND2019-15194 E192 appendix and a maximum energy of 20 MeV. Raw calculator-
  response and edge-source hashes, URLs, versions, retrieval times, and
  parameter choices travel in the runtime manifest;
- charged values use log-log interpolation only between positive nodes. A
  piecewise-constant differential source bin is integrated analytically across
  every intersected power-law segment. Neutron bins are intersected exactly
  with the explicit constant-response groups. A nonzero partially or wholly
  out-of-domain bin is rejected rather than clipped; a zero out-of-domain bin
  cannot invalidate an otherwise valid fold;
- the full transport path scores exact kinetic energy whenever a track enters
  the absorber's declared upstream response surface. For the current slab,
  side entries and downstream-face backscatter are excluded and upstream-face
  re-entry is counted. Each crossing contributes
  `S_NIEL(E)/(A_entry*|n dot u|)` so oblique crossings have the correct scalar-
  fluence normalization. The entry
  classifier and normalization area are supplied by `AbsorberGeometry`, so a
  future shape provider can define an entrance patch without changing the
  damage algorithm; that provider remains disabled until Stage 10 validation;
- DDD sums only supported proton, helium-4, and neutron crossings. `n_eq` sums
  neutron E722 hardness only; proton/helium-4 damage is explicitly excluded.
  Any nonzero unsupported DDD crossing or out-of-domain crossing makes the
  affected value, interval, rate, and accuracy `not_assessable`, while the
  partial supported sum and rejected-crossing counts remain in raw evidence;
- Geant4 nonionizing energy deposited in the finite slab divided by slab mass
  is emitted as `diagnostic_Geant4_nonionizing_deposition`. It is neither a
  replacement for nor an additive correction to absorber-entry DDD;
- case v2 requires an exact damage-library ID/version and the authenticated
  natural-silicon absorber identity. The expanded manifest/hash enters the
  resolved case, result provenance, and response-cache key, preventing reuse
  across a damage-table change. Response cells label the Stage-5 scoring
  contract separately from raw Stage-3 transport; and
- event-level DDD/`n_eq` values enter the normal accuracy controller. Statistical
  precision can be assessed, but a complete combined Stage-5 budget must
  include a case-supplied `physics_reference` component with provenance and the
  approved covariance/coverage method. No default SR-NIEL uncertainty exists.

Automated acceptance evidence:

- `stage5_core` checks authenticated loading, exact 100 MeV proton and helium-4
  nodes, the log-log geometric-mean identity, a closed-form charged interval
  fold, an exact `[1.0,1.1)` MeV ASTM fold, the 95 MeV-mb ratio, inclusive
  20 MeV endpoint, immediate out-of-domain rejection, zero-fluence semantics,
  unsupported species, and strict case capability boundaries;
- `stage5_transport` runs the normal parser, primary generator, Geant4
  multilayer geometry, provider-defined entry surface, event statistics, and
  output writers. Exact-node proton, helium-4, neutron DDD, and neutron `n_eq`
  agree with independent values within the preregistered 0.5% transport limit
  and close exactly to the raw entry/response-cell ledgers. A separate oblique
  proton case verifies the analytic `1/(A*mu)` projected-area factor;
- the transport gate requires a 25 MeV neutron to yield null/not-assessable
  outputs with explicit rejected crossings and demonstrates secondary-neutron
  sensitivity behind the maximum 10 cm aluminum layer; and
- `validation_figures` renders Stage-5 core and transport comparisons in PNG
  and SVG only after the numeric gates pass. JSON evidence is authoritative;
  figures are review aids.

Entry: selected NIEL source/standard and supported material/species matrix.  
Work:

- import pinned NIEL tables without extrapolation;
- implement edge-aware fluence integration and uncertainty propagation;
- retain target-scored non-ionizing deposition as a separate diagnostic and
  compare it with table folding only under a specifically registered benchmark;
- implement 1-MeV-neutron-equivalent fluence only for approved materials and
  convention;
- expose unsupported combinations as unavailable.

Tests: V1 analytic folds and V2 NIEL point/integral comparisons; V3 secondary
and neutron sensitivity.  
Exit evidence: every enabled particle/material pair has provenance, domain,
interpolation contract/reference tests, a justified physics-reference
uncertainty, applicable accelerator/device validation, and review. The
implementation/reference portions are complete; external physical evidence and
review remain open.
Forbidden carry-over: no analytic-surrogate value in a production column.

### Stage 6 - physically defined LET and hardness

Entry: approved LET and `H100/10` definitions.  
Implementation status: **implemented and all registered Stage-6 gates pass;
formal completion remains pending independent review.**

Implemented calculation boundary:

- case-v2 admits `LET` and `H100/10`. `LET` requires the explicit versioned
  `let_response` object. Unknown fields, wrong units, a non-integral logarithmic
  bin span, an unused response, and unsupported incident/material/model
  combinations are rejected before Geant4 initialization;
- `shieldSim.electronic_let.v1` evaluates `G4EmCalculator` unrestricted
  electronic stopping at each actual charged step's midpoint energy. The
  fluence contribution is `dx/V` in `cm^-2`; bins store event sums and squares,
  while track/dose means retain numerator/denominator squares and products;
- the total spectrum and every species/charge/ancestry/generation component are
  emitted. Component track/dose means allow a primary stopping comparison
  without silently mixing transported delta-electron fluence into that
  reference check. The all-charged spectrum deliberately retains those
  electrons as separate components;
- histogram underflow and overflow, positive- and zero-length step counts, all
  bin uncertainties, ratio sufficient statistics, Fieller intervals, and both
  LET unit systems are present in result JSON. Any nonzero out-of-range fluence
  prevents an assessable histogram result;
- the Tecplot LET file uses the same step-fluence scores and response edges as
  JSON for case-v2. The old transmitted-spectrum/Bethe file path remains only
  for the legacy CLI and is explicitly non-production;
- every forward primary-proton rear crossing contributes independently to the
  strict `>10 MeV` and `>100 MeV` event fluence sums with `1/(A |n dot u|)`.
  Multiple crossings are not de-duplicated; their within-history correlation
  is represented by the event product sum. Result JSON carries the point
  estimate, raw moments, covariance of means, Fieller status/interval, and
  crossing-count diagnostics;
- response-cell cache identity includes the complete LET-response definition
  and changes its scoring-contract identifier, preventing reuse across bin,
  policy, or scoring-definition changes; and
- the Stage-6 validation label is
  `stage6_LET_and_H100_10_implemented_and_reference_tested_pending_independent_review`.

Registered gates and their precise boundaries:

1. `stage6_core` checks the recommended 0.001--10000 keV/um grid at 20 bins per
   decade, exact endpoint/underflow/overflow behavior, constant and power-law
   partial-bin integrals, an exact `H100/10` fold, correlated ratio-of-means,
   Fieller coverage, and zero/unresolved denominators.
2. `stage6_case_boundaries` accepts the valid p/He-4 and hardness contracts and
   rejects missing/unused/malformed LET responses, gamma LET, C-12 outside
   natural Si, and C-12 with the EM reference-list default. It accepts the same
   C-12/Si case with option 4.
3. `stage6_reference` authenticates the pinned NIST and ATIMA bytes. It checks
   18 exact p/He-4 NIST electronic-stopping nodes across Al/Si/water using the
   reference-list baseline and C-12/Si at 10, 100, and 1000 MeV/u using option
   4 against ATIMA. The fixed C-12 limit is 15%; no value from this test is
   promoted into a systematic-uncertainty claim.
4. `stage6_transport` drives the case parser, primary generator, Geant4 step
   scorer, event moments, result writer, and Tecplot writer. Thin-silicon p,
   He-4, and C-12 primary-component track means must agree with the same
   deterministic midpoint stopping query within 3%. Bin plus under/overflow
   fluence must close exactly to `sum(dx/V)`. Separate maximum-step,
   production-cut, and absorber-thickness perturbations must remain within the
   approved 10% model accuracy. A weighted 50/150 MeV proton mixture must agree
   with analytic `H100/10 = 2/3` within 0.04 and its Fieller interval must cover
   `2/3`; an alpha-only case must emit null with `zero_denominator`.
5. `validation_figures` reads `stage6_core.json`, `stage6_reference.json`, and
   `stage6_transport.json` only after the numeric tests pass and renders PNG and
   SVG comparison figures plus hashes/captions in the figure manifest. Figures
   are review aids and cannot affect acceptance.

Reference provenance is under `tests/reference_let/`; all generated requests,
logs, results, evidence, and figures remain under the configured build/test
output tree rather than `/tmp`.

Exit evidence currently satisfied: each production output name maps one-to-one
to the versioned equation, scoring surface/volume, units, response range,
species/material/model domain, stochastic sufficient statistics, and pinned
reference source. Independent review is the only remaining Stage-6 sign-off
item.
Forbidden carry-over: a CRaTER estimator cannot be substituted for generic LET
or vice versa.

### Stage 7 - published-spectrum end-to-end qualification

Entry: Stages 2-6 for the requested outputs.  
Work:

- run approved GCR and SEP cases across a small representative shield/absorber
  matrix;
- generate time-resolved rates and explicitly integrated exposures for anomaly
  windows, retaining data gaps and configuration changes;
- validate response interpolation and source uncertainty propagation;
- establish event allocation/variance reduction for rare high-LET and neutron
  contributions;
- produce machine-readable reports and human review plots.

Tests: V4, the time-series tests in Section 3.7, direct-vs-matrix comparisons,
and independent-seed ensembles.  
Exit evidence: complete uncertainty/provenance reports and acceptance across
the preregistered representative matrix.  
Forbidden carry-over: results outside validated energy, species, areal-density,
material, and angular domains must be marked extrapolations or refused.

### Stage 8 - CRaTER observational validation

Entry: validated generic transport/scoring plus approved CRaTER implementation
scope and primary instrument documentation.  
Work: implement the V5a instrument benchmark, reproduce a published
CRaTER/Geant4 case, calibrate on declared intervals, then evaluate blind
intervals; optionally perform the V5b dose cross-check. Keep the CRaTER geometry
isolated as a validation fixture rather than exposing it as an instrument-design
feature.  
Tests: V5 with preregistered metrics and uncertainty model.  
Exit evidence: reproducible comparison package, residual diagnostics, passed or
failed criteria, and bounded conclusions specifying which quantities and
domains the observation supports.  
Forbidden carry-over: model parameters cannot be adjusted after blind-data
inspection without restarting the validation protocol.

### Stage 9 - production UI and regression qualification

Entry: each exposed quantity has passed its applicable prior gates.  
Work:

- map screenshot controls to versioned case definitions;
- provide an ordered layer editor that displays material, density, mass
  thickness, derived geometric thickness, uncertainty, and source-to-absorber
  order for every layer;
- show density and both mass/geometric thickness;
- identify the current absorber provider as `slab` and require all slab
  dimensions rather than filling them from the material name;
- show disabled states and reasons for unsupported combinations;
- correct case-count summaries to `number of shields x number of absorbers`;
- estimate run cost and uncertainty before submission;
- attach provenance/validation badges to outputs;
- establish reviewed regression baselines.

Tests: full UI-to-case traceability, matrix cardinality, cancellation/resume,
serialization round-trip, reproducibility, and V6.  
Exit evidence: every UI request is recoverable exactly from its run manifest,
and no UI state can bypass a scientific guard.  
Forbidden carry-over: approximate prototype labels may not appear as validated
production quantities.

### Stage 10 - arbitrary absorber-geometry extension

Entry: the slab provider, geometry-neutral interface, response hashing, and
slab physics outputs have passed their applicable gates; a concrete non-slab
use case and authoritative shape source have been supplied.  
Work:

- select one provider type at a time (analytic CSG, mesh, voxel, GDML, or a
  project adapter) and version its import contract;
- import units, reference frame, transform, material regions, sensitive
  regions, source checksum, and mass properties without inferred defaults;
- implement the validation and fatal-error checks in Section 3.3.4;
- keep geometry construction behind `AbsorberGeometry` so source adapters,
  response folding, and quantity algorithms remain unchanged;
- invalidate/rebuild response matrices by the full geometry content hash;
- document representational limitations such as mesh resolution, voxel partial
  volumes, navigation tolerance, and unresolved internal structures.

Tests: the full shape-validation ladder in Section 3.3.4, provider-specific
round trips, Geant4 overlap/navigation checks, independent reference geometry,
and all quantity convergence tests repeated for the new provider.  
Exit evidence: identical primitive cases agree across slab/analytic/mesh/voxel
representations within preregistered discretization plus Monte Carlo
uncertainty, and a nontrivial shape passes an independent reference comparison.  
Forbidden carry-over: accepting a file syntactically does not validate its
geometry or radiation response; a new provider remains disabled in
science-facing runs until its own evidence passes.

### Stage 11 - AMPS adapter

Entry: stable `SpectrumDataset` schema and a versioned AMPS output contract
provided by the AMPS team.  
Work:

- map AMPS species, grid, energy variable, normalization, frame, direction,
  location, time window, and uncertainties;
- reject unsupported schema versions and conservation/unit inconsistencies;
- compare a frozen AMPS export with an independently prepared canonical
  fixture;
- fold the same canonical dataset from AMPS and a test-file adapter and require
  identical radiation-proxy results.

Tests: adapter contract/round-trip tests, hand-calculated unit transforms,
schema-version rejection, missing-species/direction tests, and identical-input
end-to-end comparison.  
Exit evidence: AMPS changes only the source provenance/values; response and
quantity code hashes are unchanged.  
Forbidden carry-over: do not design the adapter around an inferred AMPS format.

## 8. Requirement-to-test traceability starter

| Requirement | Primary verification | Validation evidence | Earliest enabling stage |
|---|---|---|---|
| Independent shield x absorber cases | order/permutation and geometry-hash test | representative end-to-end matrix | 3 |
| Ordered multilayer shields | interface/mass checks, order-preserving round trip, split-layer equivalence | direct-vs-matrix high-Z/hydrogenous stacks in both orders | 1/3/7 |
| Five-layer/10 cm limits | exact one/five/six-layer and boundary/over-limit parser plus geometry tests | validation restricted to accepted case domain | 1/3 |
| Screenshot mass thickness | exact per-layer density conversion and unit rejection | material/thickness uncertainty study | 1/7 |
| Explicit, overridable materials | schema/fraction/isotope/source/status/hash and override tests; compare Geant4 dumps | supplier/project uncertainty studies and material-specific reference transport | 1/3-8 |
| Slab absorbers | analytic mass/path and independent-case tests | representative published-spectrum slab matrix | 1/3/7 |
| Input-defined dimensions | missing-unit/dimension rejection and serialized-case reproduction | convergence across the accepted geometry range | 1/3/7 |
| Default 10% accuracy request | default/override, adaptive-stop, systematic-floor, and resource-limit tests | combined uncertainty only within each quantity's validated domain | 1/4-8 |
| Future arbitrary absorber shape | provider validation, geometry hash, analytic/mesh/voxel equivalence | nontrivial independent shape reference | 10 |
| Time-resolved anomaly proxy | exact interval integration, gap/cadence/configuration tests | comparison with time-matched environmental observations; no causal claim | 7 |
| Published GCR ingestion | checksum/schema/unit/Jacobian tests | CRaTER GCR interval comparison | 2/8 |
| Published SEP ingestion | checksum, bounds, directions, quality flags | selected SEP observation/device benchmark, if approved | 2/7+ |
| TID | exact ledger plus NIST-controlled transport tests | CRaTER silicon dose candidate | 4/8 |
| DDD | analytic fold plus NIEL table comparisons | device/beam data selected by project | 5 |
| Silicon `n_eq` | exact standard-consistent DDD ratio | device/beam data selected by project | 5 |
| LET | histogram conservation plus stopping-power tables | instrument-matched CRaTER LET | 6/8 |
| `H100/10` | closed-form primary-proton exit-fluence threshold integrals, partial bins, covariance, zero denominator | direct-vs-response comparison; no separate observation presently identified | 6 |
| Secondary neutrons | HP-list and cut/model sensitivity | CRaTER/lunar neutral-sensitive case only if response supports it | 3-8 |
| AMPS input | adapter contract and canonical equivalence | re-run previously qualified canonical scenario | 11 |

DDD and `n_eq` need an observational or accelerator/device validation source
selected by the project. CRaTER alone does not validate semiconductor
displacement damage.

## 9. CI, data governance, and review artifacts

### 9.1 Test tiers

- **Per commit:** V0/V1, build/CLI tests, small deterministic geometry tests;
  no network access.
- **Nightly:** selected response rows, independent-seed ensembles, cut/step
  monitors, direct-vs-folded comparisons.
- **Release candidate:** full V2/V3 reference grid, published GCR/SEP matrix,
  physics-list/systematic ensembles, and observational holdouts.
- **On data/Geant4/material change:** invalidate affected response cache and run
  all tests traced to the changed item.

### 9.2 Data rules

- CI uses committed, checksum-pinned data; it never downloads mutable “latest”
  URLs.
- Raw and transformed data are stored separately. Every transformation has a
  tested script, input/output hashes, formula, and tool version.
- New upstream versions receive new filenames and review; checksums are not
  casually replaced.
- Licensing and attribution travel with the data.
- Large raw products may live in a versioned external artifact store, but the
  repository records immutable identifiers and retains compact test slices.

### 9.3 Mandatory review package per gate

Each gate archives:

- requirements and definition versions;
- source/response/result manifests;
- test code and logs;
- expected and observed uncertainty;
- residual/convergence plots and machine-readable values;
- failures and waived items, with reviewer and rationale;
- exact reproduction command/environment;
- a concise statement of what the passed gate does and does not validate.

## 10. Resolved decisions and deferred external inputs

The project decisions needed for the Stage 0 technical freeze are resolved in
Section 1.1. In particular, accuracy semantics, finite-absorber TID,
electronic-LET outputs, primary-proton `H100/10`, ASTM/ECSS device-damage
terminology, downloaded published spectra for Stage-2 tests, the GOES-16 event representations,
the validation-only CRaTER harness, explicit user-overridable materials, and
review authority are no longer open questions.

Explicit gap, adhesive, coating, or overburden regions use layer entries and
therefore count toward the five-layer limit. No later implementation stage may
silently reinterpret an approved definition; a change requires a new decision
record, schema version, and affected-test review.

The following are external/project-specific inputs, not gaps to fill with an
assumption. The named capability remains disabled until its input is available:

1. **Project materials:** A mission run using BPE, a Kevlar assembly, CFRP, an
   SiC composite, or regolith must supply the specific property record described
   in Section 3.3.5. Until then, only the explicitly labeled reference/trade-
   study records may run, and generic templates fail input validation.
2. **GOES boundary mapping:** The selected source event is the September 10-11,
   2017 GLE. The measurements retain both sensor directions, and the approved
   event-onset/event-end, source-start/derived-midpoint, and missing/quality-
   data rules are implemented and recorded in the event manifest. Transport
   requires source metadata
   sufficient to map the look directions/acceptance into the simulation frame;
   otherwise only data ingestion and event integration are enabled.
3. **Stage-5 external validation and review:** The SR-NIEL-7 v11.1 and ASTM
   E722-19 tables, exact natural-Si conventions, hashes, algorithms, strict
   domains, and automated reference/transport evidence are implemented. Formal
   completion still requires the project to select an accelerator or device-
   degradation data set with reconstructable particle/energy/fluence, silicon
   target/device, irradiation conditions, measured response, uncertainty, and
   calibration/validation partition. That evidence must support a
   provenance-bearing physics-reference uncertainty; shieldSim does not infer
   one. An independent reviewer must then examine data extraction, NIEL/E722
   conventions, domain/error policy, absorber-entry surface normalization,
   secondary sensitivity, uncertainty/coverage, generated figures, and the
   prohibition against interpreting DDD as device degradation. The review
   record must give name, role/affiliation, date, findings, deviations/waivers,
   and disposition.
4. **CRaTER preregistration:** Before Stage 8 implementation uses observational
   values, the independent reviewer must archive the instrument reconstruction,
   interval-quality rules, calibration interval(s), chronologically held-out
   blind interval(s), statistics, and acceptance criteria. The pinned product
   is selected; the split awaits recovery of its full method/quality metadata.
5. **Named reviewers:** Stage 2 is exempt from an independent-reviewer sign-off
   and closes on project-owner acceptance of its pinned-source and adapter-test
   evidence. Later physical capability and observational-validation gates that
   require independent review record the reviewer by name, affiliation/role,
   date, and finding; the roadmap does not invent personnel identities.
6. **Stage-3 independent review:** The implementation and automated acceptance
   package exist, but no reviewer has been designated. Stage 3 therefore
   remains pending. The reviewer must examine the transform convention,
   response-key invalidation boundary, scoring/ancestry and energy-ledger scope,
   five-seed equivalence/convergence evidence, physics/data provenance, and
   generated figures, then record name, role/affiliation, date, findings, and
   disposition. No code invents or substitutes this sign-off.
7. **Stage-4 independent review:** The implementation, pinned NIST package,
   63 reference comparisons, material/species transport matrix, structured
   uncertainty/correlation model, deposition partitions, numerical convergence,
   option-4 sensitivity, and review figures exist, but no reviewer has been
   designated. The reviewer must verify the registered domains and row-specific
   reference limits, NIST parsing and units, CSDA quadrature, finite-slab mass
   and TID definition, statistical/numerical gates, uncertainty provenance and
   coverage method, and all forbidden carry-over claims. The record must give
   reviewer name, role/affiliation, date, findings, deviations/waivers, and
   disposition. Until then Stage 4 remains pending, even though its automated
   acceptance suite passes.
8. **AMPS contract:** When available, the AMPS owner must provide the
   authoritative schema, units, coordinate frame, uncertainty semantics,
   version, and frozen example. No adapter implementation begins before then.

## 11. Primary references and source records

- NIST PSTAR/ASTAR stopping-power and range database: [NIST Physical
  Measurement Laboratory](https://physics.nist.gov/PhysRefData/Star/Text/intro.html).
  Stage-4 fixture material numbers are 013 (Al), 014 (Si), and 276 (water).
- NIST XCOM Photon Cross Sections Database: [XCOM home and
  documentation](https://physics.nist.gov/PhysRefData/Xcom/Text/XCOM.html).
  Stage 4 pins the exact 0.1, 1, and 10 MeV total-attenuation output including
  coherent scattering for Al, Si, and H2O.
- NIST X-Ray Mass Attenuation Coefficients: [tables and
  introduction](https://physics.nist.gov/PhysRefData/XrayMassCoef/cover.html).
  The pinned Al, Si, and water tables retain both mass attenuation and mass
  energy-absorption coefficients; Stage 4 compares attenuation only and does
  not infer absorbed dose from it.

- Badhwar-O'Neill 2020 GCR model: [Slaba et al., *Space Weather*,
  2020](https://doi.org/10.1029/2020SW002456). The publication documents the
  updated model, ion coverage, and validation against spacecraft measurements.
- AMS-02 time-resolved proton/helium data: [AMS Collaboration,
  2018](https://doi.org/10.1103/PhysRevLett.121.051101); [NASA SPASE data
  record](https://spase-metadata.org/NASA/NumericalData/ISS/AMS-02/P27D).
- AMS-02 solar-cycle He--O spectra: [AMS Collaboration,
  2025](https://doi.org/10.1103/PhysRevLett.134.051001); [official tables and
  publication record](https://ams02.space/publications/202403).
- ACE/CRIS B--Ni daily Level-2 spectra: [NASA SPASE data
  record](https://spase-metadata.org/NASA/NumericalData/ACE/CRIS/L2/P1D);
  [ACE Science Center Level-2 data](https://izw1.caltech.edu/ACE/ASC/level2/lvl2DATA_CRIS.html).
- GOES-R SGPS September 2017 data/calibration: [NOAA/NCEI event
  directory](https://www.ngdc.noaa.gov/stp/space-weather/satellite-data/satellite-systems/goesr/solar_proton_events/sgps_sep2017_event_data/);
  [Sandberg et al., 2022](https://doi.org/10.1029/2021SW002750).
- NIST charged-particle reference databases: [ESTAR/PSTAR/ASTAR
  introduction](https://physics.nist.gov/PhysRefData/Star/Text/intro.html) and
  [PSTAR/ASTAR quantities and methods](https://physics.nist.gov/PhysRefData/Star/Text/programs.html).
- NIST photon attenuation/energy absorption compilation:
  [introduction](https://physics.nist.gov/PhysRefData/XrayMassCoef/intro.html).
- Geant4 physics implementation reference: [Geant4 Physics Reference
  Manual](https://geant4.web.cern.ch/documentation/dev/prm_pdf/PhysicsReferenceManual.pdf).
- Geant4 material definitions: [Geant4 11.4 material
  database](https://geant4.web.cern.ch/documentation/dev/bfad_html/ForApplicationDevelopers/Appendix/materialNames.html)
  and [material construction/override guidance](https://geant4.web.cern.ch/documentation/dev/bfad_html/ForApplicationDevelopers/GettingStarted/materialDef.html).
- Radiation-effects calculation and device-damage conventions:
  [ECSS-E-ST-10-12C](https://ecss.nl/standard/ecss-e-st-10-12c-methods-for-the-calculation-of-radiation-received-and-its-effects-and-a-policy-for-design-margins/),
  [ASTM E3084](https://store.astm.org/e3084-17.html), and
  [ASTM E722](https://store.astm.org/standards/e722).
- NIEL implementation candidate: [SR-NIEL calculator and data project,
  version history](https://www.sr-niel.asi.it/index.php/history). Stage 5 pins
  the exact SR-NIEL-7 v11.1 natural-silicon charged tables and the ASTM E722-19
  neutron result package under `data/niel/`.
- ASTM/E192 group-edge provenance used for the Stage-5 neutron package:
  [SAND2019-15194](https://www.osti.gov/servlets/purl/1592863), Appendix B.
- CRaTER instrument data archive: [NASA PDS LRO CRaTER bundle,
  DOI 10.17189/jr9s-6678](https://pds.nasa.gov/ds-view/pds/viewBundle.jsp?identifier=urn%3Anasa%3Apds%3Alro-crater&version=2.4).
- CRaTER GCR LET method: [Looper et al., 2020](https://doi.org/10.1029/2020SW002543).
- Pinned CRaTER LET product: [Looper, version 1.3,
  2026](https://doi.org/10.5281/zenodo.22883077), CC-BY-4.0.
- CRaTER observation/Geant4 response precedent: [Spence et al.,
  2013](https://doi.org/10.1002/swe.20034).
- Candidate lunar-dose cross-check and its limitations: [Guo et al.,
  2024](https://doi.org/10.1029/2024SW004095).
