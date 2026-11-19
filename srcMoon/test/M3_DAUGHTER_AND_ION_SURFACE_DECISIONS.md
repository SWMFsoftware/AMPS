# M3 daughter physics and ion-surface models: decision proposal

Prepared 2026-11-18 from revision
`74974c2b2fc1b4cf09a73c320064a22370c0ca94` and the current uncommitted work.
**APPROVED SCOPE, 2026-11-18:** the user approved P01–P09 with `approve`.
Data/validity gates below remain mandatory. These are decision
labels, not new U/I test IDs. No physics is activated by this document.
The [decision register](M3_DAUGHTER_AND_ION_SURFACE_DECISIONS.json) records
each recommendation, approval scope and remaining evidence requirement.

The recommendation is to develop an explicitly named reduced event model
alongside the source-resolved model. The reduced model can establish
conservative execution and sensitivity to declared assumptions. Its mean
energies and angular approximations cannot establish accurate daughter escape,
pickup-ion distributions, OH brightness, or lunar ion recycling without
independent qualification. Approval of a model assumption is not evidence of
its accuracy. Missing scientific inputs cannot be supplied by approval.

The existing [physical specification](M3_PHYSICAL_MODEL_SPECIFICATION.md),
[qualification protocol](M3_QUALIFICATION_PROTOCOL.md), and
[database/Hartley2 audit](M3_CHEMISTRY_DATABASE_RESEARCH.md) remain authoritative
inventories of requirements, gates and acquired sources. Their frozen JSON
files, acceptance files, rates and roadmap are unchanged. This proposal
supersedes no approved contract and closes no milestone.

## 1. Approved decisions and limits

| Label | Recommendation for approval | What approval would authorize | What it would not establish |
|---|---|---|---|
| P01 | Retain complete D04 quiet/active endpoint families; use explicit reference cohorts and a point-Sun photon geometry for the first reduced tests | Source selection, declared initial populations and a bounded illumination model | An observed solar spectrum, warm-parent kinetics or an actual lunar isotope mixture |
| P02 | Use additive reference masses with separately tabulated formation/excitation energies and explicit isotope/state identities | A consistent nonrelativistic event ledger and its independent error analysis | Conservation of physical material rest mass, or permission to invent bath mass |
| P03 | Offer fixed source-mean excess energy and isotropic center-of-mass recoil as an opt-in reduced photon model | Two-body neutral breakup and atomic photoionization with complete products/bath accounting, once the required state/threshold crosswalk exists | A measured energy/angle distribution or a validated escape/brightness prediction |
| P04 | Use a correctly mass-weighted microcanonical sampler as the optional reduced three-neutral-body closure | A conservative, explicitly statistical model for H2O -> O + H + H after its O state and threshold are resolved | A measured dissociation mechanism; the same closure for dissociative ionization |
| P05 | Store product states explicitly and propagate source-supported radiative cascades; do not collapse excited products immediately | State tags, decay competition and acquisition/QA of complete transition data | Ground-state assignment where the source is silent, or omission of OH relaxation/quenching without a bound |
| P06 | Keep electrons in a prescribed reservoir, with signed per-event accounting; retain Voronov as the selected scalar rate model | A heavy-particle exosphere without PIC electrons, and development of a separate source-consistent electron collision provider | Derivation of outgoing electrons from a Maxwellian rate, or measured forcing from unavailable D06 data |
| P07 | For the first reduced CEX transport control, use the measured Na ground-product partial channel at its native beam-energy nodes, with an explicitly forward center-of-mass scattering approximation and physical Q | A bounded partial-channel transport experiment; retain existing He/Na total-rate services separately | A complete total-CEX law, all-state closure, or an ordinary solar-wind prediction |
| P08 | Keep a strict unsupported-ion-surface default; additionally offer an ideal absorbing collector only for verification and first-impact diagnostics | An explicit mathematical boundary with incident material, charge, momentum and energy recorded | Neutralization, implantation depth, thermalization, sputtering, residence time or lunar re-emission physics |
| P09 | Develop a separate He+/silica conditional charge-survival benchmark, then a CaSiO3 benchmark if numerical data are obtained | Source-specific laboratory operators and a concrete data-acquisition specification | A universal lunar sticking/retention coefficient or permission to contact authors |

The user approved all nine recommendations.
P03/P04/P07/P08 contain the substantive reduced-model approximations; they
must remain separately selectable. This approval does not enable a
channel whose required threshold, product state or source data are missing.
Author contact, if desired, needs a separate explicit instruction identifying
the recipients; none is requested as a prerequisite for preparing this work.

## 2. Starting implementation and the proposed production boundary

| Existing production component | Present scope | Proposed extension after review |
|---|---|---|
| `species/Photochemistry::Registry` | Complete frozen D04 rate/channel registry; daughter kinetics unavailable | Retain rates and multiplicities; bind each kinetic channel to a reviewed state/energy policy |
| `species/ReactionExecution` and `Moon::Chemistry::ParticleTransaction` in `MoonReactionAdapter` | Number-preserving all-or-nothing packets and explicit bath moments; nominal analytical fixtures | Implement P02 meanings before committing a physical event; keep retries free of new random draws |
| `Moon::Chemistry::Controller::{Install,Query,Advance,Process}` in `MoonReactionController` | Opt-in constant-history controller and continuation through actual PIC movers | Add actual trajectory exposure and an identified product-policy provider, without duplicating rate laws |
| `Moon::MoveParticle` and core `ForceCutSurfaceFirstOrder` / `ForceCutSurfaceSecondOrder` | Both force-aware tracing orders and local physical intervals | Continue every daughter through the unused original interval; propagate its state and fresh clock across MPI |
| `MoonTracing` common transport chronology | Actual Begin/Commit and source births verified for equal, constant species/block steps; external environment frozen | Use this bounded first ensemble scope; evolving environmental coupling and asynchronous science sampling need separate qualification |
| `MoonSurfaceEnvironment` moving collision service and `MoonSurfaceTemperature::ProductionStickingProbability` | Actual geometry/impact context; supported neutral Na/Ar laws; ions refused | Dispatch P08/P09 before any neutral sticking callback; preserve the actual collision epoch, local normal and surface velocity |

The names above identify source files and current namespace responsibilities;
new CLI options and provider entry points are not installed yet. The shared
core should continue to supply the tracer/transaction mechanics. Species and
surface physical policies belong in `srcMoon`, using the existing services.

## 3. P01 — source, population and illumination contract

Retain all 15 reactive D04 channels: four atomic, seven H2O and four OH_exp.
The two total-photoabsorption records are not additional deaths. Both quiet
and active endpoints are separate test cases; no interpolation between them,
new F10.7 proxy, or solar spectrum substitution is selected. Preserve
`lambda_j = lambda_j(1 AU) * visibility / r_AU^2` with each factor explicit.

For reference cohorts, propose explicitly declared ground-state parents and
the isotopes 1H, 4He, 16O, 20Ne, 23Na and 40Ar. Molecular reference parents
are H2-16O and 16OH in a declared initial electronic/rovibrational state,
initially their lowest state where a matching source is established. These
are test populations, not inferred lunar abundances. Natural-average masses
or a standard ATcT species label must not silently stand for those isotopes.
An unsupported source/parent-state mapping remains unavailable.

In particular, `OH_exp` identifies an experimental compilation, not an
electronic-state population. Review its native interpretation before mapping
it to OH(X). Initial alignment/orientation is proposed isotropic for the
reduced tests, not assumed measured. Actual warm source populations require
their own input and state-dependent cross-section coverage.

Use a point Sun in the first reduced photon tests: visibility is exactly the
production ray-classification value 0 or 1, with actual terrain and Earth
occultation where selected. The incoming photon direction is the same ray
that establishes visibility. A thermal facet cosine is not a gas-photon
visibility factor. A later finite-disc UV model must integrate and sample
the same surviving angular distribution; a fractional visibility with an
unrelated central-ray recoil direction is not the proposed model.

The current frozen-environment tracing mode is the first supported scope.
Prescribed reference epochs and common-time populations must be stated in
the input; this proposal does not infer an observation date. Evolving orbit,
thermal state and plasma require actual per-particle-time coupling before
an event-run claim.

## 4. P02 — masses, states and event energy

Recommend the reference-mass/formation-energy representation already listed
as option2 in the physical specification. For a species/state `s`, define:

```text
m_ref(s) = sum_a N_a(s) * m_ref(neutral isotope a) - z(s) * m_e
E_model(s) = |p|^2 / (2*m_ref(s)) + epsilon_formation(s)
            + epsilon_excitation(s)
```

Here `z` is integer charge in units of positive elementary charge; an electron
has `m_ref=m_e`. Nuclear identities and charge determine additive reference
mass. Physical rest mass differs by binding/internal energy divided by c².
The two are not interchangeable. The implementation must document and bound
the effect of using reference rather than physical inertial masses in both
the mover and event solver, without loosening a conservation tolerance.

The actual adapter/controller require exact agreement between catalog
`modelMassKg` and `PIC::MolecularData::GetMass`, and between catalog and PIC
charge. Preserve those checks. A new explicit physical configuration must
generate matching mass/charge tables through normal configuration; do not
use one mass in the ledger and another in acceleration. Preserve the current
legacy baseline configuration and reviewed nightly reference. Generated build
copies are not the place to implement the mapping.

Select one versioned isotope-mass inventory, NIST ASD atomic ionization/level
data, and isotope/state-appropriate 0-K thermochemistry. Recommend
[ATcT version1.220](https://atct.anl.gov/Thermochemical%20Data/version%201.220/index.php)
as the first thermochemical source to assess; it supplies 0-K and 298.15-K
values and provenance. Use gas-phase and matching isotope entries, and
transform every formation value to the declared separated-atom/electron
reference. A standard mixture entry is not an isotope-resolved measurement.
Missing entries or conflicting channel thresholds remain a data task.

For each event, define `Q = sum epsilon_in - sum epsilon_out`. Include all
photons, reservoir particles and material in incoming/outgoing accounting.
Evaluate Q from source thresholds/energies with consistent zeros; do not
subtract rounded heavy rest energies or subtract an excitation twice.
Never set ion and neutral physical masses equal merely to pass a check.

The proposed solver is nonrelativistic. Before physical activation, compute
the maximum speed and binding/reference-mass error from the selected
source/state domain, compare with an explicitly reviewed model-error budget,
and freeze a runtime speed guard. No numerical error allowance or universal
speed cutoff is selected in this proposal. A broad database energy range
does not by itself qualify the classical tracer.

Frames are explicit: perform event kinematics in a local inertial frame,
then transform to the actual LSO tracing frame at SPICE ET/TDB. At a material
boundary use the facet rest frame. Velocity transforms include the derivative
term of the SPICE state transform. During an instantaneous event, trajectory
gravity/Lorentz/non-inertial work is not a chemical energy defect.

## 5. P03 — first reduced photon-product model

**Approved development scope:** offer `source-mean-isotropic` as an opt-in experimental
policy, while retaining `source-resolved` as the accuracy-target model.
These are proposed policy names, not currently accepted input syntax.

For each channel, use its matching quiet/active mean excess energy; do not
take active energies with quiet rates. For example, H2O -> H+OH has source
mean excess3.41/4.04eV for quiet/active. Preserve the rate, energy source and
state/threshold crosswalk separately. The Huebner/Mukherjee source defines
mean excess through the photon energy above a channel threshold; its average
does not specify the event distribution.
[Source](https://doi.org/10.1016/j.pss.2014.11.022).

Define a representative incoming photon in the parent rest frame by
`E_gamma = E_threshold,channel + E_mean,channel`. Include its momentum
`(E_gamma/c) * n_hat`. With all incoming momenta/energies and explicit final
internal states, compute available product center-of-mass translation:

```text
P_total = sum incoming material momenta + photon momentum
K_available = E_including_photon - sum epsilon_out
              - |P_total|^2 / (2*M_products)
```

All terms use the same P02 energy zero. Photon recoil means K_available is
not identically the tabulated mean excess. Refuse a negative K_available,
unresolved threshold, or incompatible state mapping; do not truncate to zero.
Where a source channel includes unresolved internal excitation, assigning
all excess to translation is not automatic: obtain the missing partition
or leave that channel unactivated. An assumed ground state is not the fix.

For two massive products, `mu=m1*m2/(m1+m2)`, select a unit direction
uniformly on the sphere, and set `p_rel=sqrt(2*mu*K_available) * n_hat`.
The outgoing momenta are `m_i*P_total/M_products +/- p_rel`. This is the
part of the existing comet two-body mathematics worth reusing after audit.
It gives exact classical energy/momentum closure for the declared model.

For atomic photoionization the two products are the ion and electron. Trace
the ion; record the sampled electron in the prescribed reservoir. Include
the incident photon as a reservoir loss. An isotropic electron direction
is a reduced-model assumption, not a Na dipole law. It must not be described
as the measured solar-weighted electron angular distribution.

For a neutral two-body breakup, trace both products at the same event
position and time. Each product occurrence represents the original physical
number, has the correct species weight correction, receives a fresh reaction
clock, and consumes only the unused parent interval. No comet W/dt factor
or parent internal-state byte copying is allowed.

**Validity limit:** this model fixes only a representative energy and an
angular closure. It does not reproduce the true variance/tails or correlations.
It is suitable for declared-model conservation and sensitivity work, not an
accuracy claim for escape fractions, daughter spectra or state-resolved
brightness. Such quantities depend nonlinearly on the missing distributions.
No agreement with observations will be used to choose the approximation.

## 6. P04/P05 — molecular multiplicity, states and relaxation

For the three-neutral-body H2O -> O+H+H channel, propose a uniform classical
microcanonical distribution at the available translational energy as a
separately selected reduced model. This is a model choice, not measured water
breakup physics. In the outgoing center-of-mass frame, construct:

```text
g_i: independent standard-normal three-vectors
p_i = sqrt(m_i) * g_i
p_i <- p_i - (m_i/M) * sum_j p_j
K_draw = sum_i |p_i|^2/(2*m_i)
p_i <- p_i * sqrt(K_available/K_draw)
```

The sqrt(m_i) factor is essential. Drawing equal-variance velocities and
rescaling, as in the generic comet many-body construction, gives a different
mass-dependent distribution. Handle K_available=0 deterministically and
refuse invalid draws; an allocation retry cannot redraw kinematics. Preserve
the two H occurrences as equal-status daughters, with no doubled event rate.
An independent Jacobi-coordinate oracle has energy fraction
`x ~ Beta(3/2,3/2)` for the two relative modes. This verifies the selected
sampler, not its physical accuracy.

Do not automatically extend this statistical recipe to the three
dissociative-ionization channels. Their correlated electron/ion/neutral
partition needs its own physical closure. Recommend keeping their transport
disabled until a complete source-bound law is reviewed. Retain their original
hazards in the registry; a request for an incomplete full molecular kinetic
family is rejected before execution. No dropped-channel renormalization is
permitted. A separately labelled channel-isolated supporting experiment
cannot be reported as I22 completion.

P05 recommends explicit state tags and their independent production rules:

| Product family | Recommended state treatment | Activation condition |
|---|---|---|
| O(3P), O(1D), O(1S) | Distinct identities; resolved J where required by the adopted transition model | Source branch/state mapping and complete applicable NIST OI transitions frozen |
| OH products | Retain electronic and relevant rovibrational populations | Acquire source-supported OH state probabilities/cascade data; `OH_exp` is not a population |
| H2 and H2O+/OH+/O+ products | Retain states represented by source thresholds and future rates | State/threshold/formation-energy crosswalk exists; no automatic lowest-state assignment |
| Unspecified O in O+H+H | Keep unresolved in inventory | Resolve from native source before physical transport |

Use `tau_u=1/sum_l A_ul` and branch probability `A_ul/sum_l A_ul`, with an
independent exponential decay clock competing with chemistry/impact. Include
the emitted photon in the ledger. A decay changes state, not atom number.
Acquire a complete reviewed [NIST ASD level](https://physics.nist.gov/PhysRefData/ASD/levels_form.html)
and [transition export](https://physics.nist.gov/PhysRefData/ASD/Html/lineshelp.html);
an arbitrarily truncated wavelength query cannot establish a complete lifetime.

Do not impose unit O(1S)->O(1D) branching or instant OH ground-state relaxation.
Neglect of collisions/quenching, alignment or fine-structure effects requires
a documented domain and omission bound. The state-resolved work of
[Tabone etal2021](https://doi.org/10.1051/0004-6361/202039549) is a concrete OH
population/cascade acquisition route; its restricted wavelength/parent-state
coverage does not complete the solar channel family.

H2O/OH surface migration and binding remain a separate physical task in
literal U24/I31. The gas-product proposal supplies no water surface sticking,
site distribution, residence time or cold-trap release law.

## 7. P06 — electron reservoir and electron-impact products

Recommend prescribed electrons rather than transported PIC electrons for
the first exosphere scope. Preserve the already approved explicit scalar
Maxwellian/Voronov mode and its source domain. This decision selects the
environment representation, not a full electron collision law.

The signed reservoir change for a photon event is the gained photoelectron
minus the incident photon. For `e_in + A -> A+ + e_1 + e_2`, it is both
outgoing electrons minus the incoming electron. Net reservoir charge is -e,
balancing the new tracked +e. Every momentum/energy term comes from the
declared sampled states before closure is checked; it is not computed as
minus an unexplained tracked residual. Prescribed fields/reservoirs do not
solve sheath formation or self-consistent plasma feedback.

**Concrete recommendation for the missing impact sampler:** keep physical
electron-impact daughter transport OFF while acquiring a source-consistent
incident/outgoing provider. A Voronov coefficient cannot determine the
collision-conditioned electron distribution. Required incident sampling is
proportional to `f_e(v) * sigma(|v-v_A|) * |v-v_A|`; use the corresponding
joint outgoing law, exact threshold and ion state from the same model.

Start with a bounded He electron-beam subproblem, where the existing NIST
100-eV singly differential data can check the emitted-energy marginal. It
still lacks the complete correlated angular law. [NIST SRD107](https://physics.nist.gov/PhysRefData/Ionization/intro.html)
provides BEB total cross sections and selected singly differential electron
energy data; those are distinct from a fully differential distribution.
The Na CCC source is the corresponding acquisition priority for sodium.
Do not pair NIST/CCC products with a Voronov rate by an undocumented rescaling.

If a new differential source implies a different Maxwellian total, expose it
as a separately named rate-and-product model, retain the existing Voronov
mode and compare both totals without calibration. That source choice requires
its own numerical review before activation. Do not truncate a Maxwellian's
unbounded tails to fit a finite table. A bounded declared beam is a different
physical experiment and must retain that label.

D06 contains no qualified measured electron forcing for the intended event.
Prescribed reference n_e/Te are not observations. R4.1/I23 measured-event
requirements stay open; approval of a reservoir does not remove them.

## 8. P07 — bounded first CEX daughter experiment

The acquired IAEA D119990 data are the measured partial channel
`H+ + Na(3s) -> H(1s) + Na+(2p6)` at12 native energies4–25keV/u.
Recommend using those nodes without interpolation in the first supporting
beam experiment, beginning at the native4keV/u point. The specified
projectile density is a reference input, not a measured solar-wind density.
The existing Barnett-He and Morgan-Na total-rate registries remain unchanged.
[Native source](https://db-amdis.org/collisiondb/datasets/119990/).

Implementation audit clarification: that native label specifies the outgoing
states but leaves the initial Na state unspecified. Na(3s) is the approved
reference cohort, not recovered experimental metadata. Its physical
source-state crosswalk and bibliography/uncertainty qualification remain open;
the reduced component must preserve that limit rather than claim it measured
an initial 3s population.

Propose center-of-mass forward scattering as an explicitly named experimental
closure for this partial channel, with the physical energy defect retained:

```text
Q = I_H(1s) - I_Na(3s)
K_relative,out = K_relative,in + Q
p_relative,out = sqrt(2*mu_out*K_relative,out) * n_relative,in
```

Restore the full center-of-mass momentum with the P02 product masses. This
is an energy-corrected forward model, not an unconditional velocity swap.
Reject impossible energies, zero/undefined incident directions or out-of-range
sources. Transport both physical heavy products, with the removed bath proton
recorded if the incoming wind is prescribed rather than stored in PIC.
Each product occurrence represents the same reacted target number. No PIC
proton can be removed twice, and a reservoir proton cannot be counted as a
pre-existing tracked particle.

The forward angle is proposed for mechanism/software sensitivity only. The
measurement supplies the partial cross section and identities, not that angle.
A generic classical elastic scattering model is not an independently measured
nonresonant capture angular law. Do not combine this ground-state partial
sampler with the all-state Morgan total or infer the missing capture states.
IAEA D107581 is composite capture+impact-ionization and is unsuitable as its
replacement. The4keV/u lower endpoint also does not cover ordinary~1keV/u
wind; there is no low-energy extrapolation.

For a later He full-state model, the acquired OPEN-ADAS low-energy table at
0.3–500keV/u is a promising separate source. Its assessed state completion,
angle correlations, redistribution terms and cascades still require review.
Independent state/angle marginals do not establish their joint distribution.
Full U17/I24 stays open while those complete requested families are absent.

## 9. P08 — ion impacts: strict mode and an ideal collector

**Recommended default:** retain explicit refusal of an unsupported physical
ion-regolith outcome. Preflight should reject a request whose specified
species/material/energy domain has no complete boundary policy. If an
unexpected unsupported impact occurs, retain its incident record and report
the missing capability; a failed physical assumption must not be hidden by
falling through to neutral-Na/Ar sticking.

**Additional mode approved for development:** a perfectly absorbing collector
for verification, distinct from a lunar retention model. This is a prescribed
mathematical boundary, independent of regolith temperature:

1. At the actual collision point/epoch, transform incident velocity to the
   moving facet rest frame and record species, charge, state, represented
   number, energy, direction, facet and normal.
2. End that tracked history and transfer its full material, charge, momentum
   and model energy to a named collector reservoir. Preserve separate kinetic,
   internal/formation and charge entries; do not label their sum as heat.
3. Produce no reflected particle, neutralization event, sputtered target atom
   or delayed release. Preserve charge in the collector ledger; no electron
   is created/consumed implicitly. Fields remain explicitly prescribed, so
   this mode has no self-consistent surface charging response.
4. Record the boundary impulse and the corresponding rotating-frame/inertial
   bookkeeping consistently. The transferred momentum is the actual incident
   momentum under this boundary condition, not a fitted budget residual.

There is no implantation depth, neutral fraction, sticking coefficient,
residence time or accommodation factor in this mode. Its incident tally is
the transfer across an artificial boundary, not a prediction of lunar retained
abundance. First-impact flux/energy/angle diagnostics are useful; the gas
population after absorption is a prediction only of this ideal collector.
It must not appear as validated lunar recycling, surface heating or trapping.

An alternative diagnostic may terminate a trajectory at its first unsupported
impact and store it as a censored/unresolved history without choosing any
outcome. Such histories cannot supply post-impact gas densities or release
fluxes. Recommend retaining strict behavior for science and using the ideal
collector only when explicitly requested, rather than silently censoring a
normal science run.

### What a physical ion-regolith kernel must eventually contain

Condition on incoming species/isotope/charge/state, E, mu, material, surface
temperature and relevant irradiation/electrostatic state. Partition the
incident projectile's fate among all supported reflected charge states and
retained states, with normalized probabilities. **Sputtering and electron
emission are conditional secondary yields that can accompany those fates;
they are not extra mutually exclusive fates added to that probability sum.**

Reflection requires joint outgoing energy/angle/charge/state. Retention
requires depth/site and subsequent diffusion/trapping/release laws. Target
sputtering debits the material inventory; its secondary products carry their
own charge, energy, momentum and state. Surface neutralization must account
for transferred electrons and chemical energy. Do not replace these inputs
with a single neutral sticking or Arrhenius constant.

No complete numerical Na+/Ne+/Ar+ lunar-mixture outcome kernel was qualified
by the research. A constant albedo or all-ions-neutralized assumption would
be an additional physical hypothesis, not a value supplied by the database.
I do not recommend either as a default lunar model.

## 10. P09 — first source-supported ion-surface benchmark

Recommend **He+ on the specified isotope-exchanged silica** before a generic
lunar material mixture. The acquired [Tellez etal2012 paper](https://pure.tue.nl/ws/portalfiles/portal/3503988/583864115491353.pdf)
studies0.6–7keV, normal incidence,145-degree detection. Implement its
conditional charge-survival relation and TableI regimes only after exact
numerical transcription/QA, retaining its mechanism-specific velocity
definition. Do not bridge the oxygen transition regime with an invented fit.
The output is `P(charge=+1 | detected scattered He, experimental class)`.
It is not `P(retained | incident He+)`, nor a full outgoing-angle sampler.

The separate [Brotzner etal2026 CaSiO3 study](https://journals.aps.org/prb/abstract/10.1103/3ks5-1c83)
at1–3keV offers a more directly silicate benchmark, with different scattering
classes. Its native numerical spectra remain request-only in the acquired
paper. The BCA-derived fractions must remain distinguishable from measured
spectra. Neither experiment supplies Na+/Ne+/Ar+ material outcomes.

The concrete author-data request would ask for: incident species/isotope/
charge, energy/angle/dose and surface preparation/composition; calibrated
spectra and detector acceptance/response; uncertainties/covariance and class
definitions; charge-resolved outgoing energy/angle distributions; absolute
incident normalization and retention measurements if any; BCA inputs/version;
and data reuse terms. If unavailable, label each requested quantity unavailable
instead of reconstructing it from a charge fraction. No email is sent by this
proposal or by approval of its software scope.

## 11. Inputs, failure behavior and delivery order

All model selections must be explicit in the runtime .in file using the
repository grammar, not compiled defaults. Implement proposed policy names
only after review. Each run must preserve source endpoint, input population,
isotopes/states, energy convention, product model, photon geometry, bath,
ion boundary, data/policy hashes and reduced-model label. No runnable example
for these new policies is claimed here.

Bind model review, data identities and supported claim scope to the provider
and run manifest. The existing `KinematicPolicy::qualified` Boolean is not
permission to declare a reduced model physically accurate. Do not bypass its
current evidence/admission checks by setting that flag. Any required change
to distinguish reviewed declared-model execution from full physical
qualification must be specified and tested explicitly. The formal owners
must continue to refuse unsupported accuracy claims and predecessor gaps.

Missing required state/threshold/source or formal predecessor -> SKIPPED /
NOT QUALIFIED. An explicitly requested unsupported normal science
configuration is refused with an actionable diagnostic. Corrupt data,
nonfinite arithmetic, broken dispatch, resource failure or nonclosing ledger
-> ERROR; an evaluated physical/numerical acceptance miss -> FAIL. Neither
an unknown boundary nor a difficult channel becomes zero probability.

After approval, implement in this order:

1. Freeze the source/state/threshold and mass/energy crosswalk, with a
   reproducible acquisition/QA script and independently derived reaction Q.
   Unsupported entries stay explicit; never manufacture a complete table.
2. Implement P02 and P03 in the existing Moon packet provider, initially the
   atomic/two-body subset whose crosswalk is complete. Preserve all existing
   rates and complete-family admission checks. Add P04 and P05 where their
   real data exist. Supporting subset results are not full I22 PASS.
3. Wire actual trajectory exposure and event-conditioned channel selection
   through both production movers/common-time transport. Use the actual gas
   photon ray, force path and event epoch; retain fresh daughter clocks.
4. Implement P08 as a separate explicit boundary, then P09's conditional
   laboratory operator; continue the complete physical material-data work.
5. Implement P07's bounded partial-channel fixture. For P06, implement full
   electron transport products only after the chosen incident/outgoing source
   and rate consistency are reviewed. Keep scalar-only controls available.
6. Execute the fixed local/linked protocol and preserve all predecessor gates.
   A remaining data or accuracy gap keeps the corresponding formal ID open.

## 12. Verification contract and artifacts

Extend existing owners within their meanings; assign no new U/I ID. Update
the detailed executable test contract before code changes. Preserve every
frozen tolerance and the protocol's separately pending statistical choices.

| Existing owner/scope | Production check | Independent oracle and required metrics |
|---|---|---|
| U15/U14; I21/I22 | Complete rates, survival, channel selection and photon products | Exponential/matrix-exponential solution; independent source transcription; per-channel represented nuclei/charge, event timing, mass/energy convention |
| U24 supporting gas scope | Two-/three-body molecular products and cascades | Two-body reduced-mass identities; independent Jacobi Beta distribution; state-count/radiative-rate identities; H multiplicity and H/O budgets; full U24 surface scope remains required |
| U16/I23 | Selected scalar electron rates; future collision-conditioned joint products | Native fit evaluator, independent collision quadrature and supplied differential distribution; reservoir electron count and energy/charge balance |
| U17/I24 | Source-bounded CEX partial/full families as actually configured | Native node/fit evaluation, Q and COM identities; no ground-partial/total substitution; native domain and product counters |
| Existing tracer/boundary support; proposed gaps recorded if no frozen owner applies | P08 collector and P09 laboratory conditional operator | Analytic crossing plus actual owned gas/collector sums; independent TableI evaluation; conditional output normalization, frame/epoch/angle |

For each new executable test, record all standard specification fields: purpose,
production entry points/configuration, input/domain, oracle, units/frames/time,
seed and independent parent count, procedure/metrics/criteria, expected result
and failures, data identities, predecessors, artifacts and status semantics.
No new surface-physics gate is assigned to an unrelated chemistry ID.

Use both mover orders, one/eight MPI ranks, equal/unequal reference weights,
fractional corrections, allocation retries and source births. Include a
moving facet test and force/chemistry/boundary event ordering. Common-time
steps must meet the current equal-constant-step admission condition. A
supporting field fixture is not a measured solar wind. Continue to require
the source/geometry convergence and formal I13–I20 predecessors.

Development verification must share the already approved3600-second envelope
and200GB aggregate RAM limit, including setup, clean build, plots, restoration
and cleanup. Default8MPI/OMP1, with explicit1-rank controls. Do not promise
complete rare-branch statistics within this envelope; preserve SKIPPED where
the frozen statistical design needs the separate cluster tier. No new time,
RAM, seed, sample-size or confidence approval is requested by this proposal.

Preserve logs, result.json, run_manifest.json, exact input, source/build and
executable hashes, source data/policy hashes, RNG/parent-stream layout,
MPI/OMP layout, event/collector ledgers and resource history. Produce CSV,
PNG and EPS comparisons against the declared independent oracle. Label
reduced-model plots as experimental; do not call model outputs observations.

## 13. Authority, source identities and work-package record

Source roadmap SHA256:
`712e1152bf50af6d9ee1dad8f8a585809ca8b92cf42ac6fc7a009f6c0253ee17`.
Data-root roadmap SHA256:
`d10eca1c539f9cf82940236d9e1aa36af665d43e82ba475925950d582fbf49c4`.
The source adds approved development-order/resource/status amendments; those
control work order. Common U/I meanings and criteria still control formal
qualification; native labels/data-root instructions control interpretation.
Both roadmaps are unchanged and their difference is not silently reconciled.

Reused source identities:

| Actual local source | SHA256 |
|---|---|
| D04 `processed/phidrates/photo_rates_1au.csv` | `3aa54915a08258371dbe975949e6afb2ab74618164fadca8bceebe876c926847` |
| D04 `processed/phidrates/thresholds.csv` | `a3335b3b37c3d275a4f72e92ae8e5b1f6ecfcbc4e2d84735cb20694dc93b48f1` |
| Tellez2012 publisher PDF in the retained ion-surface packet | `bf554c667bb01d03952d151f86e4e348b706f20630c8fe896eedf15f6fa30ed2` |
| Brotzner2026 publisher PDF in that packet | `11a94ab8d81f99329fa6423f3adf41e892bcce9f82fb98f684b406ac90518bb6` |
| IAEA D119990 extracted native12-row CSV in the atomic packet | `2a887ad43497a7f9e2b3d7b6173ae1b2597cc893d897926762c181946c92db93` |

The research report links the immutable canonical source provenance and all
research manifests. The12 IAEA rows are published measurements extracted
from native numeric arrays, not digitized pixels. The proposed forward-angle
law is an approximation; it is not an extracted measurement. New NIST/ATcT
state/threshold/cascade products are still to be acquired/frozen and have no
invented hashes. Some publisher/IAEA fetches failed during this follow-up;
the already retained native files and manifests remain the numerical basis.

This work changes only proposal/status/README documents and adds a decision
register. It does not select a runtime model, change any source data, run a
simulation/build, inspect observational residuals, change acceptance/reference
files, stage, commit or push. Documentary link/schema/source-identity checks
are distinct from scientific PASS counts. Formal M3 remains SKIPPED / NOT
QUALIFIED, with M1/M2 predecessors open. Approval of P01–P09 authorizes only
their explicit scope and does not by itself complete M3.
