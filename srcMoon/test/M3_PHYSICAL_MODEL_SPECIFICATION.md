# M3 physical model specification and unresolved data contract

Inspection date: 2026-11-18. Source revision:
`74974c2b2fc1b4cf09a73c320064a22370c0ca94`, with concurrent approved work.

**Physical activation remains disabled pending qualified evidence.** This
document specifies what a complete physical event model must contain. It does
not supply missing distributions, approve a recoil approximation, set new
acceptance tolerances, or declare I21–I24 PASS. The companion
[versioned JSON](M3_PHYSICAL_MODEL_SPECIFICATION.json) is a machine-readable
inventory. Every current physical process is `BLOCKED_MISSING_DATA`; missing
physical quantities are `null`, not zero. The separate numerical/qualification
protocol owns test matrices, sample sizes, confidence limits and execution
budgets. Nothing here replaces that protocol.

The [per-process data requirements](M3_PHYSICAL_DATA_REQUIREMENTS.csv) distinguish
actual acquired arrays, published papers, request-only data, and unresolved
source interpretation. The earlier [kinematic source proposal](M3_KINEMATIC_REFERENCES_PROPOSAL.md),
[molecular inventory](M3_MOLECULAR_KINETICS_PROPOSAL.md), and
[ion-surface/CEX source proposal](M3_ION_SURFACE_AND_CEX_DATA_PROPOSAL.md)
provide acquisition detail. A research download is not a READY validation
package or a qualified production model.

## 1. Authority, approvals and scope

The source roadmap SHA-256 is
`712e1152bf50af6d9ee1dad8f8a585809ca8b92cf42ac6fc7a009f6c0253ee17`.
The canonical-data roadmap SHA-256 is
`d10eca1c539f9cf82940236d9e1aa36af665d43e82ba475925950d582fbf49c4`.
They differ: the source copy records approved development-order, resource and
implementation-status amendments absent from the external copy. Those user
approvals permit implementation/supporting verification before predecessor
closure. The unchanged U/I meanings, formal predecessor order, physical
requirements and native canonical data metadata still control qualification.
Neither document was edited or silently reconciled here.

| Decision | Authority/status | What it establishes |
|---|---|---|
| First physical qualification uses a common advancing ensemble time and qualified common initial state | APPROVED design | Sources, bath, force, geometry and sampling must refer to the same physical time. No initial epoch/population or qualification artifact is invented here. |
| Scalar prescribed Maxwellian electron mode with Voronov coefficients | APPROVED rate closure | Explicit electron density and temperature determine the approved scalar rate. It does not specify incident/outgoing event distributions or authorize PIC electrons. |
| Conservative complete product packets and qualified force/cut-surface tracing | APPROVED design | Preserve represented number, exact nuclei/charge, product state, remaining physical interval and fresh daughter clocks. Resource shortage cannot discard a daughter or redraw an event. |
| Source-supported total rates and named Na source-knot log-log CEX interpolation | APPROVED bounded rate services | Preserve complete rate families and source domains. These services do not qualify product transport. |
| Unsupported physical channels remain OFF | APPROVED requirement | Refuse an incomplete requested family rather than delete its difficult channels and redistribute their probability. |
| Physical mass/energy zero, electron-reservoir representation, joint kinetic laws, initial-state populations and surface outcome laws | PROPOSED / BLOCKED | Actual source data, complete domain definitions, independent checks and review are required before selection. No implicit ground state, isotropy, zero recoil or zero Q is selected. |

The JSON policy IDs identify these responsibilities. An approved design is
not a qualified instance. Its actual source files, hashes, production mapping,
initial-state proof and review are still required. Setting `qualified=true`
or writing an evidence-name string cannot satisfy those requirements.

## 2. Audited production boundary

The read-only source inventory is
[`source_audit.json`](../../test_output/srcMoon/m3-physical-model-spec-20261118T013000Z/source_audit.json).
It binds the inspected files; later concurrent code must supply its own build
identity. The table describes the audited baseline, not a claim that future
implementation is unchanged.

| Production entry point | Implemented support | Physical limitation |
|---|---|---|
| `Photochemistry::Registry::LoadFrozenD04`, `SelectChannel`, `TotalDestructionRateS` | Transactional native CSV/hash parser, all 17 catalogue rows, 15 reactive channels, quiet/active rates, states supplied by source, multiplicities and explicit photon visibility | `daughterKineticsAvailable=false`; scalar mean energy is not a daughter PDF. |
| `ElectronImpact::MaxwellianVoronovTable` and the Moon electron adapter | Complete native 406-row parse, four neutral-atom fits, explicit prescribed scalar electron density/temperature | No differential collision-conditioned electron sampler; no molecular electron-impact fit is supplied. |
| `ChargeExchange::Registry` | Source-bound He fit, named Na log-log law, strict source domains, explicit monoenergetic `n sigma g` | No state/Q/angular daughter policy, no automatic thermal-tail prescription and no production particle CEX activation. |
| `ReactionExecution::ValidateReaction`, `MakeEventPlan`, `PendingEvent` | Integer nuclei/charge and per-product occurrence accounting; explicit bath/classical moments; complete packet reserve/stage/commit | A caller must supply a genuinely qualified kinetic/state/mass-energy policy. The current additive model-mass check is not a proof of physical rest-mass/Q accounting. |
| `MoonReactionAdapter::ParticleTransaction` | Actual PIC species/weight/mass/charge checks, complete packet allocation and daughter continuation | Matching the generated mass table proves configuration consistency, not the physical accuracy of that table. |
| `MoonReactionController::Install` and actual Moon mover hooks | Actual conservative production-path integration for declared constant-history hazards in the inspected baseline | A constant-hazard mechanical control is not a physical, changing solar/plasma exposure model. New common-clock/exposure work must independently qualify its physical providers. |
| M2 surface callbacks | Reviewed Na/Ar neutral kinetics; ion requests fail closed in the M2 adapter | Neither an ion neutralization/trapping law nor a water/OH migration law follows from neutral Na/Ar coefficients. |

The native controller probe explicitly requires nominal neutral/ion mass
equality and labels its energy law zero-Q mechanics. Its `.25/.75 s^-1`
analytical hazards and elastic control are numerical fixtures, not lunar
photoionization/CEX physics. Their successful tests establish storage,
bookkeeping, RNG, event ordering and mover integration only.

## 3. Frozen available rate data

Canonical D04 and D05 are both physically stored under
`moon_validation_data/photochemistry/`; there is no separate local
`electron_impact/` root. The package README, provenance, native headers and QA
were read. The current archived extraction QA is PASS. This audit verified
**11 processed products, 21 provenance descriptors and 33 actually-used raw
descriptors, with zero hash mismatches**. Overlapping provenance descriptors
are counted as checks, not as additional independent science files.

| Product, relative to `photochemistry/` | SHA-256 |
|---|---|
| `processed/phidrates/reaction_channels.csv` | `2d9fa1c8eb72d3bbba0b91ab4c9a912e6ebce2abc80b251096df2f41a4354061` |
| `processed/phidrates/photo_rates_1au.csv` | `3aa54915a08258371dbe975949e6afb2ab74618164fadca8bceebe876c926847` |
| `processed/phidrates/thresholds.csv` | `a3335b3b37c3d275a4f72e92ae8e5b1f6ecfcbc4e2d84735cb20694dc93b48f1` |
| `processed/phidrates/cross_sections.csv` | `e900eb5514b616dcbbde7e3ca124b5c06d4d22a28a3cd429ed58a14c1933c860` |
| `processed/phidrates/branching_ratios.csv` | `430bb68558a43973f0622094049c87004dd522a6288a48935027856f48229ba7` |
| `processed/electron_impact/lunar_neutral_coefficients.csv` | `da7391979a2136c072bd8dc02932af7989a8030ac208fa5982ac975aa1394126` |
| `processed/electron_impact/voronov_coefficients.csv` | `8c1e33a9ed64d95de884b01548b7e2757339859b7c8090e996f842c4f6ca9225` |

The JSON/audit additionally bind all 11 processed products, the native
`cfit.dat/.txt/.f`, raw PHIDRATES responses, scripts, QA and provenance.
No canonical file was written, filtered, repaired, or reprocessed.

D04 supplies two separate solar endpoints: activity 0.0 and 1.0, optical
depth 0.0, reference distance 1 AU. Native wavelength Å becomes nm by 0.1;
native cross-section cm² becomes m² by 10^-4. Rates are s^-1 and the native
mean-excess-energy field is eV. These conversions are archived processing
operations, not a new source model. Endpoint interpolation, an intervening
solar spectrum, attenuation, or a new state population is not inferred.

### 3.1 Complete photon family inventory

Every row below is an original reactive channel. The photon is an incoming
environment term even though it is not printed as a material CSV product.
`e` has charge -1; each `+` ion has charge +1; neutral products have charge 0.
Every listed event preserves exact H/O/atomic nuclear identities and net
charge when electrons are included. For `O + H + H`, two separate H occurrences
must each represent the parent weight. The unexplained state of `O` stays
unknown; it must not be renamed O(3P).

| Parent | Original products | Quiet rate, s^-1 | Active rate, s^-1 | Source state information |
|---|---|---:|---:|---|
| He | He+ + e | 5.64e-8 | 1.68e-7 | Parent/residual detailed populations not supplied in these rows |
| Na | Na+ + e | 7.26e-6 | 7.91e-6 | Same gap; no automatic 3s/ground-state substitution |
| Ne | Ne+ + e | 1.84e-7 | 5.87e-7 | Same gap |
| Ar | Ar+ + e | 3.42e-7 | 7.74e-7 | Same gap |
| H2O | H + OH | 1.03e-5 | 1.76e-5 | OH electronic/rovibrational populations unspecified |
| H2O | H2 + O(1D) | 5.97e-7 | 1.48e-6 | O term explicit; H2 internal population unspecified |
| H2O | O + H + H | 7.54e-7 | 1.91e-6 | O term unspecified; correlated three-body momenta missing |
| H2O | OH+ + H + e | 5.54e-8 | 1.51e-7 | Joint ion/neutral/electron states and momenta missing |
| H2O | O+ + H2 + e | 5.85e-9 | 2.21e-8 | Same type of gap |
| H2O | H+ + OH + e | 1.31e-8 | 4.07e-8 | Same type of gap |
| H2O | H2O+ + e | 3.31e-7 | 8.28e-7 | Ion internal state and electron/recoil distribution missing |
| OH_exp | O(3P) + H | 1.20e-5 | 1.37e-5 | Term explicit; fine-structure populations unspecified |
| OH_exp | O(1S) + H | 8.33e-7 | 2.11e-6 | Excited O must remain state-resolved |
| OH_exp | OH+ + e | 2.43e-7 | 6.43e-7 | Ion/electron joint distribution missing |
| OH_exp | O(1D) + H | 7.01e-6 | 1.76e-5 | Excited O must remain state-resolved |

The exact source text, channel IDs, mean energies, partial domains, source
references and unknown fields are retained per row in JSON. The two additional
catalogue rows are H2O/OH_exp `total_photoabsorption`: **metadata, not extra
destructive events**. Do not add them to the 15 reactive hazards or assume
absorption implies destruction. Conditional channel sampling uses
`k_i / sum(all declared reactive k_i)`, leaving each original hazard intact.
An unsupported outgoing branch blocks requested family transport; it does
not authorize dropping or renormalizing that branch.

`OH_exp` preserves the source's experimental `OHe` compilation. Its physical
parent is OH. It is neither an isotope nor a defined OH(X)/OH(A), vibrational,
rotational or thermal population. Product OH may use this rate compilation
only after its physical state/population has been shown to match the source
scope or after a separately reviewed state-resolved law is supplied.

### 3.2 Source interpretation gaps that block spectral kinetics

The 2305 archived cross-section rows include duplicate samples and zeros.
Their minima/maxima describe archived samples; they are not permission to
interpolate across gaps, extrapolate, repair a threshold or infer a continuous
physical support. No such operation is performed by this specification.

Two **positive atomic samples fall below their own listed threshold**:

| Species | Wavelength, nm | Photon energy, eV | Cross section, m² | Listed threshold, eV | Raw locator |
|---|---:|---:|---:|---:|---|
| He | 51.75 | 23.958299214145 | 7.47e-22 | 24.5873554184746 | `brnout_001.raw`, line197/column2 |
| Ar | 82.49 | 15.0302095324524 | 2.98e-21 | 15.7596347408481 | `brnout_001.raw`, line125/column2 |

These are unchanged source bytes, not a failed hash or an AMPS residual.
Determine the native state/endpoint/threshold interpretation from authoritative
source documentation before an energy-resolved physical sampler consumes them.
Do not clip the samples, change thresholds, infer an excited parent population,
or alter the published integrated rates to remove the inconsistency.

The existing molecular audit also found tiny positive OH_exp O(1S)/H entries
extending to 6132 nm, about 0.202 eV. Their physical meaning is unresolved;
the molecular threshold table is absent. State-resolved molecular thresholds,
source floors/tails and source printing precision must be interpreted before
spectral transport. An arithmetic mismatch between rounded partial and total
cross sections is not authority to renormalize source numbers or invent errors.

### 3.3 Electron and charge-exchange rate scope

Voronov's native fit is
`U=dE/T_e[eV]`,
`k=A(1+P sqrt(U))/(X+U) U^K exp(-U)`.
The four selected neutral atoms each have a source temperature domain
**1–20000 eV**, where temperature means `k_B T/e`, not incident electron
kinetic energy. Native A in cm³/s converts by 10^-6 to m³/s. The complete
406-row table remains preserved. There are no H2O/OH electron coefficients in
this package. Its ionization-energy fit parameter is not automatically an
exact physical state threshold: for example, native Na `dE=5.1 eV` differs
from the separate D04 threshold `5.13907097104346 eV`. Preserve both source
roles; do not silently substitute or force equality.

| CEX rate family | Supported scalar source | Native domain | Missing full kinetics |
|---|---|---|---|
| H+ + He → H(unresolved) + He+(unresolved) | Barnett1990 ORNL-6086/V1 A-32, Appendix I Chebyshev convention | 99–11000000 eV/u, target-rest projectile energy per u | All captured/residual states, state-specific Q, joint recoil/angle law |
| H+ + Na(3s) → H(unresolved) + Na+(unresolved) | Morgan1985 Table2, named approved source-knot log-log law | 1000–200000 eV/u, same source energy convention | Same gaps; Na3s target scope cannot be transferred to another state |

Publication SHA-256 values are respectively
`3bb3fdb4be8f0f806d378083be8798cf3ce347714c3d06d58f90ee483b10c030`
and `7c412cb93f827b3412b6da000e674c63237cacb0e201733a6f5a4f431359952e`.
Na interpolation is a named model between actual source knots, not newly
measured data. There is no extrapolation, clamping, missing-range zero or
universal Ne/Ar/alpha law. `ColdBeamRateS` takes energy, relative speed and
density explicitly; a physical policy must establish their consistent mass/
energy mapping. The upper source range does not establish a classical-mover
error bound.

## 4. Physical state, mass and energy convention: required review

Exact material invariants are integer nuclear/isotope identities and electric
charge. **The sum of material rest masses is not a conserved invariant in an
ionization, dissociation or state-changing reaction.** Its change exchanges
rest/binding/internal energy with kinetic energy and the radiation/environment.
Therefore an artificial bath mass chosen to close a numerical residual cannot
establish physical conservation.

Two consistent representations can be reviewed; neither is selected here:

1. Use state-specific physical rest masses and total four-momentum. For a
   massive particle `E=sqrt(m^2 c^4+p^2 c^2)`; for a photon
   `p_gamma=(E_gamma/c) n_hat`. Total incoming and outgoing four-momenta,
   including the bath/surface/radiation, must agree. Electronic/rovibrational
   excitation already included in a state mass is not added again as a
   separate excitation energy.
2. Use an explicitly documented reference-mass and formation/excitation-energy
   decomposition, with a proved nonrelativistic error bound over the selected
   domain. Each species/state has one documented energy zero; physical
   threshold/Q values and excitation offsets must correspond to that same
   decomposition. Any difference between bookkeeping mass and physical
   inertial mass must remain explicit and enter the error analysis.

Under the second convention, define the signed state-dependent chemical
energy release `Q = sum(epsilon_reactants) - sum(epsilon_products)`. In an
isolated event, kinetic-energy gain plus emitted-photon energy is Q, after
including incoming photon/projectile/electron energy and every outgoing
environment term. A threshold already defined for O(1D), for example, must
not have the O(1D) excitation subtracted twice. Use stable source-threshold/Q
evaluators rather than subtracting rounded heavy rest energies and calling
the cancellation physical accuracy.

The existing production planner currently checks:

```cpp
// ReactionExecution.cpp: audited numerical model contract
afterMass - beforeMass + bathMass == 0;  // within explicit closure bounds
energy = 0.5 * modelMassKg * speed * speed + internalEnergyJ;
```

The actual implementation uses guarded long-double accumulation and explicit
bounds; the excerpt explains its scope, not exact executable syntax. Before
physical activation, map the approved physical convention into these fields
and checks, or amend their physical meaning/implementation with independent
tests. Merely broadening a mass tolerance, using equal neutral/ion masses,
setting internal energy to zero, or compensating with a fictitious bath term
does not implement a physical energy defect. Existing nominal controls remain
useful numerical tests and must retain their limited label.

Required input files are an authoritative isotope/species/state mass inventory,
state energies and uncertainties, a threshold/Q crosswalk for every channel,
and the reviewed convention/validity proof. None is supplied here by inference.

## 5. Signed bath and represented-event accounting

For each **one physical parent event**, define bath deltas as final minus
initial untracked environment. Then require separately:

```text
Delta N_nucleus,tracked + Delta N_nucleus,bath/surface = 0
Delta q_tracked         + Delta q_bath/surface         = 0
Delta p_tracked         + Delta p_bath/surface/radiation = 0
Delta E_tracked         + Delta E_bath/surface/radiation = 0
```

The energy convention in section4 determines E. Gravity/Lorentz work during
flights and work of a rotating/non-inertial frame belong in the trajectory/
field ledger, not in an invented instantaneous chemical Q. Bath moments must
be computed from declared sampled incoming/outgoing physical states **before**
checking conservation. They must not be defined afterward as minus the
measured tracked residual.

The parent represents `W_parent * correction_parent` physical particles.
Every unit product occurrence represents the same number; multiplicity two
means two such occurrences or an independently qualified equivalent packet.
Each daughter's correction uses its actual reference weight, without a
species-birth-timestep factor. These numerical requirements are approved; they
do not determine physical daughter velocities or internal states.

All state/chemical/environment transfers must share the actual event epoch
and immutable original interval. A repeated event-plan capacity retry reuses
the same sampled event and clocks. A daughter receives a fresh stochastic
reaction clock and only the parent's remaining physical time. The common-time
ensemble and initial condition require their own actual qualification proof;
asynchronous local mechanical success cannot substitute for it.

## 6. Photon operator and atomic product model

For a specified parent state, the physical channel hazard is the integral of
incident spectral photon flux times the corresponding conditional cross
section over the declared energies, directions and polarization. Event
selection samples that same integrand. Any use of the D04 integrated endpoint
rate with another spectral/differential source needs an explicitly reviewed
normalization/consistency contract; a total rate cannot be inverted uniquely
into a photon spectrum.

The source spectrum must identify energy versus wavelength density, photon
versus energy units, 1-AU normalization, exact quiet/active endpoint, epoch
support and angular/polarization convention. Doppler or recoil approximations
need a supported domain and independent bound; no numeric bound is selected
here. The gas photon-visibility factor is evaluated at the actual particle
position/epoch. It is not a facet projected-incidence cosine. Terrain, finite
solar disc, Earth occultation and any assumed spectral-disc weighting require
explicit production geometry and source conventions. Thermal finite-disc
settings do not automatically define UV photon intensity/polarization.

For `A + photon -> A+ + e`, sample the photon and residual-ion state, then
the conditional electron energy/direction and correlated ion recoil. The
incident photon momentum belongs in closure. Isotropic emission, a universal
dipole parameter, a fixed electron energy equal to the source mean, and
zero ion recoil are all physical approximations requiring separate evidence
and review; none is enabled by this specification.

**Frame/time contract:** production trajectories use the explicitly configured
LSO frame and SPICE TDB/ET seconds; surface material uses the documented Moon
body frame. Source data may use a laboratory, target-rest or center-of-mass
frame. Specify the event-inertial frame and every epoch-dependent transform.
For rotating frames, transforming a velocity requires the state-transform
derivative term, not just a three-vector rotation. Carry photon axes and
polarization bases as well as particle vectors. Store the source convention
and transformation in each physical policy's evidence.

## 7. Electron-impact event model and measured forcing

The approved scalar Maxwellian rate is `lambda = n_e k_Voronov(T_e)`.
Its positive-temperature domain and invalid/zero-density behavior remain the
production coefficient contract. This scalar closure does not specify
electron bulk velocity, collisions of a moving target with a drifting bath,
the conditional incident energy, or either outgoing electron.

For a supplied differential cross section, the collision-conditioned incident
electron distribution is proportional to
`f_e(v_e) sigma(|v_e-v_parent|) |v_e-v_parent|`, followed by its **joint**
two-electron/residual-ion distribution. A scalar Maxwellian coefficient does
not determine those factors. If that differential source integrates to a
different total from the approved fit or needs drift corrections, preserve
the discrepancy and resolve the model/source consistency before activation.
Do not silently rescale a differential table to reproduce a convenient total.

**Proposed representation for review:** retain an explicitly prescribed
electron reservoir for the first exosphere physical scope, with each incident
and outgoing electron recorded in signed charge, momentum and energy ledgers.
For `e + A -> A+ + 2e`, the bath loses one incident electron and gains two
outgoing electrons: net one bath electron, balancing the new tracked positive
charge. For photon ionization the bath gains the emitted electron and loses
the incident photon. Those moments still require real joint data. An
untracked electron is not an unaccounted electron.

This is a proposal, not selection of full kinetics. Explicitly transported
PIC electrons would additionally require their actual masses, time scales,
forces, boundaries, sampling and a consistent electromagnetic closure. The
current force tracer rejects self-consistent PIC/deposition modes. Scalar
prescribed Maxwellian approval does not authorize enabling that unsupported
solver or running both generic electron reactions and the Moon controller.

For a **measured event**, R4.1/D06 remains a real dependency. The current D06
QA is ERROR, with zero original campaign CDFs, zero selected intervals and
zero driver rows. Its selected THB/THC merged variables are bulk ion moments,
GSE magnetic/velocity vectors and SSE position; they do not provide the
needed electron moments or a resolved alpha abundance. Do not infer `n_e`
from bulk ion density or `T_e` from ion thermal speed.

Select an actual independent event/probe/electron product before examining
model residuals, acquire native quality/time/frame metadata, preserve gaps,
and qualify its physical mapping to the modeled environment. A prescribed
bath fixture can verify its declared subproblem earlier. It cannot turn the
missing measured-forcing branch into a formal PASS or remove its dependency
administratively. The parent formal protocol must preserve the original I23
scope and predecessor evidence.

## 8. Molecular kinetics and excited-state propagation

Full M3 includes all seven H2O and four OH_exp reactive channels, not merely
an atomic subset. For two-body dissociation, specified incoming photon/parent
and outgoing internal states determine the available relative translation;
the angular law remains an independent input. Three-body channels require
correlated energy/momentum distributions, including the two indistinguishable
H occurrences where applicable. Dissociative ionization also needs correlated
photoelectron/ion/neutral kinematics. One mean excess energy does not fix any
of those joint distributions or a solar-weighted escape/brightness prediction.

Preserve O(3P), O(1D), O(1S), OH electronic/rovibrational states and relevant
H2/H2O+ states explicitly. Exact state storage may use tags or distinct PIC
species, but the physical populations, mass/energy mapping, rates and future
relaxation must remain bound to those identities. No unspecified product is
made ground-state by a compile-time species index.

Acquire a versioned complete OI energy/transition export, its underlying
references, units and uncertainties. [NIST ASD levels](https://physics.nist.gov/PhysRefData/ASD/levels_form.html)
exposes critically evaluated level data; its
[line-data documentation](https://physics.nist.gov/PhysRefData/ASD/Html/lineshelp.html)
describes transition probabilities and state classifications. These official
interfaces were checked on 2026-11-18; a complete reviewed OI kinetic export
has **not** been acquired/frozen by this task.

For complete independent radiative contributions, `tau_u=1/sum_l A_ul` and
branch probabilities are `A_ul/sum_l A_ul`. Do not force O(1S) to O(1D) with
unit probability or collapse a state instantly because a quoted lifetime is
short. Decay preserves the atom, changes internal state, and accounts for
the emitted photon's energy/momentum on the actual trajectory/time. Quenching,
OH relaxation/predissociation and secondary chemistry need their own source
laws, or independently demonstrated omission bounds reviewed before execution.
No numerical lifetime, quenching constant or omission bound is selected here.

The literal **U24** also covers sticking/residence, thermal hops, migration and
closed H/O surface/gas budgets. A successful gas-only D04 table/packet test
does not complete it. Water/OH-specific state/material/temperature-dependent
binding/site, sticking, residence, relaxation and release laws plus the actual
production callbacks remain required. Na/Ar neutral laws cannot supply these
missing inputs. Local gas-only checks may be reported under their existing
meaning; U24 and dependent full I22/I31 scope stay SKIPPED until their real
requirements are evaluated. No new ID or PASS waiver is introduced.

## 9. CEX kinetic requirements

For a prescribed projectile population the rate is the collision integral
`n_projectile integral f(v) sigma(E_rel) g d^3v`; a monoenergetic beam reduces
to the existing explicit `n sigma(E) g`. The physical policy must state
projectile isotope/charge/state, target initial population, relative-frame
energy per u and its mapping to physical speed/mass. Source-valid cross
section ranges do not authorize an arbitrary energy-speed conversion.

Outgoing projectile/capture state, residual target-ion state and correlated
scattering/recoil determine Q and the joint momenta. Decide explicitly which
neutral/ion products are transported and which are recorded as bath terms.
Replacing projectile H+ by a bath neutral H changes its charge, momentum and
energy even if it is not stored in PIC. Neither a forward velocity swap,
all-ground-state capture nor factorization of separately measured state and
angle marginals is selected here.

The acquired Rice He angular table is a finite-angle **5-keV, state-summed**
measurement. Its first angle 0.002 degrees conflicts with the primary article
abstract's 0.02-degree limit; resolve that source discrepancy before using an
integration domain. Keep duplicate-angle/uncertainty records and missing tails.
The acquired OPEN-ADAS ADF01 He/H table contains state-resolved integral
cross sections at 70–15000 keV/u, not a complete ordinary-keV angular law.
Its n/nl/nlm hierarchical sums must not be added as independent branches.
Neither file supplies the missing joint ordinary-keV state-angle distribution.

Na state-differential theory at 0.5–5 keV is a concrete acquisition lead, not
an acquired numerical joint law. A Maxwellian projectile distribution has
unbounded energy support; no tail clipping/renormalization, zero extension
or extrapolation outside the rate/differential source domain is approved.
Acquire support or propose and independently bound a limited physical model
for review. A source-supported bounded beam can be considered separately;
it does not become solar-wind validation merely because its rate is known.

## 10. Ion-surface boundary outcome requirements

A physical ion reaching the lunar boundary requires an outcome law conditioned
on incident species/isotope/charge/state, energy and angle in the actual
surface-rest frame, material/composition/structure, temperature and relevant
irradiation/charge state. Use the actual terrain normal and moving/rotating
surface velocity. A laboratory 145-degree detector is not a generic facet
reemission distribution.

The joint outcome model must distinguish reflected ions, reflected neutrals,
implantation/trapping (including depth/state), charge transfer/electron emission,
sputtered target material and later diffusion/desorption/release. Its complete
probabilities must account for every incident event within its supported
domain. A measured charge fraction conditional on a detected scattered atom
is not the unconditional implantation or sticking probability.

Track incoming/outgoing nuclear inventory, surface and emitted charge,
particle/radiation momentum, thermal deposition, electronic/chemical storage
and outgoing excitation/kinetic energy under the same approved convention.
Do not declare all ions neutralized/absorbed, reuse neutral-Na sticking, or
discard incident charge/energy without a separately approved bounded model.
Gas-only collision fixtures that never reach a material boundary may verify
that limited subproblem; they do not qualify the coupled lunar ion boundary.

Available evidence includes He+ on specified silica at 0.6–7 keV and normal
incidence with 145-degree detection (Téllez2012), and CaSiO3 at 1–3 keV with
the same nominal detector angle (Brötzner2026). The latter numerical spectra
are request-only, not acquired. Neither gives a universal lunar-regolith
joint reflection/implantation/release law. Apollo/enstatite EELS/microscopy
constrain named reservoirs/materials; accompanying SDTrimSP output remains a
simulation rather than an observed event yield. Actual research-file hashes,
MD5 checks and acquisition failures are preserved in the linked source
proposal; no author was contacted.

## 11. Concrete evidence required before physical READY

For each process, the JSON `missing_inputs` and `required_policy_ids` are
mandatory, source-bound responsibilities. A future review package must contain:

1. Actual immutable physical data files and metadata, hashes, processing scripts,
   source domains/units/frames, uncertainties and independent numerical checks.
2. Initial/product state and isotope populations, consistent mass/threshold/
   internal-energy and environment representations, plus signed per-event
   energy/momentum/charge/nucleus accounting.
3. The complete conditional sampling law and its support, normalization,
   interpolation and any independently bounded approximation or omission,
   reviewed before production results. Missing joint data remain a blocker.
4. Production provider/dispatch mapping, source/build identities, common-time
   initial-state proof, and all predecessor/operator/numerical evidence required
   by the separate formal qualification protocol.
5. A material outcome law for a surface-coupled claim, literal U24 closure for
   the full molecular scope, and actual R4.1/D06 evidence for measured forcing.

Missing scientific inputs or supported capability cause formal
`SKIPPED / NOT QUALIFIED`. Corrupt hashes/schema/metadata or an inconsistent
implementation cause ERROR. An available supported criterion that is executed
and not met causes FAIL. Source extraction/hash PASS, a scalar coefficient,
a nominal zero-Q mechanics test or a `qualified` flag is never a physical or
observational PASS. No new acceptance threshold is introduced here.

## 12. Work-package record

New files are this specification, its JSON inventory and the CSV requirements
table. The audit inspected 33 actually-used raw descriptors and all 11 existing
processed products without changing any canonical data; no numerical record
was generated as an observation, filtered, repaired or discarded. Source
identity checks passed with zero mismatches. The two atomic threshold conflicts
and pre-existing molecular/state/differential/forcing gaps remain explicit.

No C++/runner/configuration change, application build, MPI simulation,
observational scoring, commit or push was performed by this specification
package. Documentary schema/identity/arithmetic checks are supporting checks
only. Formal I21–I24 remain BLOCKED/SKIPPED; no physical kinetic policy was
manufactured or activated.
