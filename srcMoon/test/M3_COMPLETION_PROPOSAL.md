# Proposal for completing M3: chemistry promotion and budget closure

Prepared 2026-11-17 against source revision
`74974c2b2fc1b4cf09a73c320064a22370c0ca94`.
**PROPOSAL: new designs and physics choices below require review before
activation. This document implements no chemistry and changes no acceptance
criterion, reference, roadmap dependency or qualification result.**

## 1. Recommended completion target

Complete the existing R3.1–R3.4 requirements through a common, conservative
production reaction service, with real parent/product transfer and integration
with the normal Moon mover. Keep the reviewed legacy Na mode available.
Deliver the atomic photon/electron components first, while resolving the
charge-exchange and molecular data contracts in parallel. An atomic subset is
a useful intermediate delivery; it is not completion of the unchanged M3.

Formal completion requires the declared I21–I24 gates, their local gates and
the applicable M1/M2 predecessors. Early implementation/supporting verification
is already authorized; early formal qualification is not. In particular, the
current scalar electron coefficient checks do not establish particle conversion,
product kinematics, measured plasma forcing or event-level conservation.

The recommended immediate work package is sections 5–7, stages A–C: freeze the
registry and test specifications, implement conservative transactions, then
verify their integration with the actual mover. Product transport remains OFF
where its physics contract is incomplete. This produces reviewable code without
inventing the missing daughter distributions.

## 2. Authority and present development state

The required source/data roadmaps and external `AGENTS.md` were read. The
roadmaps differ; neither was edited or silently reconciled for this proposal.

| Roadmap | SHA-256 | Controlling contract |
|---|---|---|
| [Source roadmap](../development_and_validation.md) | `712e1152bf50af6d9ee1dad8f8a585809ca8b92cf42ac6fc7a009f6c0253ee17` | Literal U/I meanings, production requirements, numerical exit criteria and approved implementation-order/resource amendments. |
| `/data/vtenishe/moon_validation_data/development_and_validation.md` | `d10eca1c539f9cf82940236d9e1aa36af665d43e82ba475925950d582fbf49c4` | Original qualification order and data dependencies. Native labels and external `AGENTS.md` control data interpretation. |

The source copy records the approved M2 development exception, one-hour
development/explicit cluster policy and November 17 early-M3 authorization.
The external copy retains the original ordering. Both preserve the R3.1–R3.4
numerical requirements. User approvals control development order; they do not
turn supporting checks into formal I-series PASS results.

| Area | Actual state at inspection | Consequence |
|---|---|---|
| M0 | I01–I10 PASS; reviewed I04 reference retained and zero-byte nightly diff. | Preserve legacy behavior and regression evidence during M3 changes. |
| M1 | I11–I13 closure remains open; a passing reviewed I13 convergence reference is unavailable. | No claim that M3 formal predecessors have passed. |
| M2 | Substantial surface/environment implementation and supporting verification exist; formal I14–I20 remain unqualified. D02 metadata and other operator/physical gates remain unresolved. | Implement M3 under the approved exception; keep formal qualification blocked. |
| M3 atomic electron rates | Production `species/ElectronImpact` and explicit Maxwellian/Voronov callback verified locally and by a native eight-species/eight-rank scalar fixture. | Reuse the service. The particle executor remains disabled. |
| M3 photon chemistry | Reviewed neutral-Na loss path exists; no qualified D04 multi-species loader or Na→Na+ packet transfer. | Implement a separate opt-in mode, preserving the old mode. |
| M3 charge exchange / molecules | Research proposals and source tables exist; complete channel/product kinetics are not installed. | Missing physics/data remain explicit blockers. |
| I21–I24 ownership | Frozen IDs exist, but individual integration directories/contracts/runners are absent; the current development runner emits four SKIPPED placeholders. | Implement actual runners and register them; placeholders are not executable integration tests. |

The previous core integration preserves all 919 historical PASS cases, including
843 CORE cases. Its known strict-coverage FAIL and CME timeout ERROR remain
documented; “preserved” does not mean the complete registry is all green.
These are existing results, not tests repeated for this proposal.

## 3. Current production paths and reusable upstream services

This inventory is the starting point for implementation, not a suggestion to
patch generated `build/` copies.

| Source / entry point | Audited behavior and proposed treatment |
|---|---|
| `input/moon.input`, `Moon::Chemistry::ProcessParticle` in `Moon.cpp` | The configured photolytic callback records/deletes neutral Na and explicitly refuses compiled Na+. Preserve it as the legacy branch. An unused inline conversion helper is not proof of working ion transfer. |
| `Moon::ExospherePhotoionizationLifeTime` in `Moon.h`; `species/SodiumCase1.cpp` | Legacy lifetime is 130500 s except explicitly configured paper-case behavior. No inverse-square scaling; orbit-on lunar/Earth gates and an interior-Earth short-lifetime shortcut exist. Preserve regression semantics in legacy mode; do not import that shortcut as a physical photon rate. |
| `species/ElectronImpact.h/.cpp`; `Moon::ElectronImpactIonizationRate` | Existing immutable D05 table, explicit density/temperature, native domain refusal and neutral/ion mapping. Reuse without loosening its verified tolerances. |
| `src/models/exosphere/Exosphere_Chemistry.cpp` | Generic product births and parent loss are independently sampled, with product dt used for births. The old cached photochemical table is not the new D04 registry. Do not use it as a conservative transient reaction executor. |
| `src/pic/pic_electron_impact_ionization.cpp` | Independent births/loss and product-dt weighting remain. A product-loop `i` shadows the spatial cell index. Keep this executor disabled for the new Moon mode; enabling the rate callback does not fix these problems. |
| `src/models/chemical_solver/ReactionNetwork.h` | Reuse validated reaction/stoichiometric containers and parent indexing where applicable. This class does not supply physical units, scattering, rates for every process or a particle integrator. |
| `src/models/chemical_solver/ChemistryDiagnostics.h` | Reuse per-reaction and per-product-term accounting, including repeated identical products. Add explicit Moon MPI reductions and bath/element/state metadata. |
| `src/pic/pic_photolytic_reactions.cpp`, `PIC::ChemicalReactions::PhotolyticReactions` | Reuse the existing five-callback batched walker and collective-safe deferral checkpoints. Install all callbacks or none; preserve NULL/default behavior for other applications. |
| `srcHartley2/CometChemistry.cpp`, `srcHartley2/RadialSimulation.cpp` | Existing batching/installation examples. Their stationary `W/dt` carriers and population controls are not automatic transient lunar inventory physics. |
| `PIC::ParticleBuffer`, `PIC::Mover::DispatchParticleMover` | Reuse record allocation, individual corrections, list ownership and actual movers. Register persistent chemistry fields before particle-record layout is frozen. Avoid moving/inserting a record twice. |
| `MoonSurfaceEnvironment.h/.cpp` | Reuse production BVH/Sun/Earth ray masks; add verified solid-angle gas exposure and volume-point origin support. The existing area fraction/surface offset is not automatically that operator. |
| `Moon::TotalParticleAcceleration` in `Moon.h` | Charged products use Lorentz force with typical compiled E/B when the coupler is OFF, otherwise interpolated fields. New-mode inputs must declare field values/provider and frame rather than inherit a nominal solar wind silently. Retain U06's actual Lorentz-kernel check. |
| `Exosphere::SurfaceInteraction::StickingProbability` in `Moon.cpp`; `ProductionStickingProbability` in `MoonSurfaceTemperature.cpp` | M2 deliberately refuses Na+/Ar+ impacts: ion neutralization/trapping laws are unconfigured. Do not fall back silently to legacy neutral laws for new ions. |
| `MoonM2Budget.cpp` | Existing two-species Na/Ar budget has an intentional limited contract. Extend through an explicit chemistry ledger; do not relabel it as an all-species/charge/energy budget. |

Important upstream limitation: the existing Query/Advance/Process trajectory
event hooks execute in the first-order tracer, whereas normal Moon uses the
force-aware second-order mover. The latter does not currently invoke those
hooks. The batched walker solves allocation/dispatch ownership, not event
location along an accelerating trajectory or multi-daughter continuation.
This is the one likely shared-core extension; section 6.4 limits its scope.

## 4. Data that are available, and data that are missing

Read-only package inspection covered README, provenance, QA, processed/native
labels, selection state and identities. The retained audit has **59/59 matching
file identities**: 11 processed chemistry products, 33 native inputs, five
charge-exchange research files, three driver configuration/QA files and seven
driver products. Matching identities do not establish physical completeness.
The complete audit and independent reviews are retained under
[`test_output/srcMoon/m3-completion-proposal-20261117/`](../../test_output/srcMoon/m3-completion-proposal-20261117/).

### 4.1 D04 photon data

Root: `/data/vtenishe/moon_validation_data/photochemistry`.
Extraction QA passes; 17 catalog entries include two absorption totals and
**15 reactive channels: four atomic, seven H2O and four experimental OH**.
Do not add absorption totals as extra destruction channels.

| Atomic parent | Quiet rate at 1 AU, s^-1 | Active rate at 1 AU, s^-1 |
|---|---:|---:|
| He | 5.64e-8 | 1.68e-7 |
| Na | 7.26e-6 | 7.91e-6 |
| Ne | 1.84e-7 | 5.87e-7 |
| Ar | 3.42e-7 | 7.74e-7 |

Both endpoints can be verified now. No intermediate solar-activity or
time-varying spectrum law follows from two endpoint rates. A requested run
must explicitly select its endpoint or supply a separately documented spectrum.

Molecular hazards/products are inventoried in
[M3_MOLECULAR_KINETICS_PROPOSAL.md](M3_MOLECULAR_KINETICS_PROPOSAL.md).
Preserve all seven/four branches, multiplicities and native identifiers:
`OH_exp` comes from PHIDRATES `OHe` and is a physical OH parent with experimental
provenance, not a new isotope or a declared rovibrational population.
Conditional probabilities `k_c / sum(k_c)` over the complete declared reactive
family do not alter the individual archived hazards. Omitting difficult
channels and renormalizing the retained ones is prohibited.

The archive does not supply joint daughter translation/internal-state/angular
distributions. Scalar mean excess energies are not such distributions. Some
channel-specific citations and numerical molecular thresholds are missing;
rounded wavelength partial sums differ from printed totals; tiny experimental
OH O(1S) tails extend to photon energies incompatible with an assumed ground
parent. Preserve these values and seek authoritative interpretation. Do not
repair, clip or renormalize them to satisfy a gate.

### 4.2 D05 electron-impact coefficients

The native Voronov table contains 406 atomic/ionic fits. The selected He/Na/Ne/Ar
neutral fits support inclusive `1 <= k_B T_e[eV] <= 20000`. Native coefficient
units cm^3/s convert once to m^3/s. This is an already Maxwellian-averaged rate
coefficient, not a differential cross section or a truncation of the electron
energy distribution. No additional Maxwellian-tail rule is needed to evaluate
this fit inside its domain.

It does not supply secondary-electron energies/angles, unique ion recoil or
H2O/OH electron-impact coefficients. Keep those unsupported processes OFF.
The original Fortran and independent high-precision reference remain useful
coefficient oracles; neither is a product-kinematics oracle.

### 4.3 D06 measured plasma

`drivers/qa_report.json` currently reports **ERROR**: zero raw CDFs and zero
processed driver rows; the CSV files are header-only. The selected variables
also do not supply the required electron moments, species-resolved proton
distribution or alpha fraction. Do not infer `n_e=n_i`, `T_e=T_i`, proton
temperature from an undocumented thermal-speed convention or a fixed alpha
fraction.

This blocks measured-event qualification. It does not block explicitly
prescribed Maxwellian-electron or monoenergetic-proton analytical verification.
Those are model inputs, never downloaded observations. R4.1/D06 remain required
when measured event forcing is claimed.

### 4.4 Charge-exchange candidates

The existing [charge-exchange proposal](M3_CHARGE_EXCHANGE_PROPOSAL.md) retains
primary ORNL/NIST sources and their extraction caveats.

| Candidate | Acquired support | Remaining limitation |
|---|---|---|
| H+ + He → H* + He+ | ORNL A-32 total capture fit/table, 99–11000000 eV/u. | Need reviewed product-state/Q/differential kinetics; high-energy velocity conversion needs an appropriate relativistic convention or bounded classical scope. |
| H+ + Na(3s) → H* + Na+ | NIST/Morgan total capture knots populated over 1–200 keV/u. | No acquired lower-energy table. Proposed log-log interpolation is only between positive populated knots and still requires review. Total capture does not select H states or angles. |
| Ne/Ar targets; alpha capture | Literature leads, incomplete numerical/product/forcing contracts. | Explicitly OFF; no borrowed/scaled cross section. |

Source energy is target-rest projectile energy per atomic mass unit, with
documented unit conversion. A 400 km/s proton example lies below the acquired
Na table's lower bound. Do not choose a wind speed to avoid missing data.

For normalized projectile distribution `f(u)`, the hazard is
`nu(v) = n_p integral f(u) sigma(E(|u-v|)) |u-v| d^3u`.
It is generally not `n_p sigma(mean g) mean g`. Accepted collision partners
must follow the collision-biased distribution proportional to `f(u) sigma(g) g`.
A finite cross-section interval does not define an unbounded Maxwellian's
complete integral; acquire the missing support or approve a physically justified
distribution/tail bound. This limitation is different from D05's complete
Maxwellian fit.

## 5. Decisions to freeze before the affected feature is enabled

These are concrete recommendations, not approvals inferred from this document.
Unresolved rows do not prevent inventory/parser/transaction development; they
do prevent activation and the corresponding full qualification claim.

| Decision | Recommended choice | Evidence/approval still required |
|---|---|---|
| M3 scope | Retain full R3.1–R3.4, including molecular R3.2; deliver atomic components first. | An atomic-only milestone would require an explicit roadmap amendment. |
| New transient execution | One common physical step across reactive species, chosen for the most restrictive mover; individual particle weights ON; OpenMP=1 initially. | Approve this initial execution scope. Refuse asynchronous species clocks until separately designed/verified. Different reference weights remain supported. |
| Photon input | Explicit quiet/active endpoint, table/hash, Sun-distance convention, solar angular-radiance assumption and lunar/terrain/Earth occultation policy. Test both endpoints. | Freeze each actual run's choice and justify any uniform-disc/point-source approximation; no invented intermediate spectrum/default. |
| Atomic products | Explicit tracked-ion/electron-or-bath policy, recoil/energy/angular model and uncertainty bound. | D04/D05 total rates do not define these. Parent-velocity inheritance is an approximation requiring a quantitative justification and approval. |
| Ion forces and surface fate | Explicit prescribed E/B or qualified field provider/frame; retain fail-closed unsupported ion impacts until a documented neutralization/trapping/re-emission policy is approved. | Existing U06 verifies force arithmetic, not arbitrary ion trajectories or surface physics. Chemistry-only microfixtures may avoid impacts but cannot qualify surface-coupled transport. |
| Mass and energy accounting | Versioned species/state mass convention, exact nuclei/charge accounting and explicit photon/electron/projectile bath transfers; binding/internal energy accounted consistently. | Current neutral/ion pairs have equal nominal masses (Na/Na+ both 3.82e-26 kg). Decide model-mass versus electron-inclusive mass-energy scope before claiming exact conservation; no silent shared mass-table edits. |
| First charge exchange | Begin with source-supported He and Na rate kernels; keep transport OFF until products/kinematics are complete. | Approve exact channels, valid energies, interpolation/fit and bath distribution; obtain state/differential data or review a named approximation with a bounded error. Rate-only/neutral-loss-only is not full I24. |
| Molecular transport | Complete seven/four hazards, explicit H/O/charge products and O(3P/1D/1S) states; acquire joint daughter distributions and state evolution. | Parent state populations, ambiguous O label, thresholds/tails, energy sharing, radiative lifetimes and any quenching policy remain to be resolved. A mean-energy isotropic recipe is not automatically acceptable. |
| Statistical protocol | Freeze independent trial unit, cases, seeds, coverage, family accounting, power/precision and counts before AMPS outputs. | Review the concrete planning method in section 8; confidence level alone does not specify adequate sample size. |
| Numerical protocol | At least three predeclared dt levels at the same physical final time, independently controlled hazard/mover errors and weighted-layout checks. | Freeze the numerical error allocation before results; the approved I17 thermal tolerance does not automatically apply to chemistry. |
| I22 rate metric | Propose 2% for total parent removal rate, source-precision deterministic partial-rate checks, and the separately frozen branch-frequency tests. | Review whether the existing 2% wording instead requires each partial ensemble rate. Rare-channel precision would need a different, potentially much larger cohort; do not silently choose the cheaper interpretation. |
| U24 / M4 overlap | Explicitly name the chemistry sub-contract for I22 while preserving full U24's surface-migration meaning. | Requires reviewed dependency clarification; otherwise implement required U24 surface capability early or leave full M3 open. Never emit full U24 PASS for gas chemistry alone. |
| I23 / R4.1 overlap | Qualify prescribed-bath rate/executor components separately from measured-event forcing. | Clarify the milestone's event-driven claim; if literal M3 closure requires measured forcing, complete R4.1 first. No D06 proxy. |

Also correct, through a reviewed roadmap reconciliation, the statement that
an existing documented U17 lunar channel can be enabled first. Current U17 is
not implemented and has no approved lunar channel. A core hydrogen model is
not that missing channel. This proposal does not edit either roadmap.

## 6. Proposed production architecture

### 6.1 Versioned reaction registry and input contract

Keep species-dependent services under `srcMoon/species/`. Proposed files are
`Photochemistry.h/.cpp`, `ChargeExchange.h/.cpp`, and a shared Moon chemistry
registry/executor in `srcMoon/`; these names are prospective, not current APIs.
Reuse existing `ElectronImpact` and upstream stoichiometry/diagnostics.

Each enabled reaction record must contain stable ID/version, process family,
reactants/products/multiplicities, electronic/internal states, compiled AMPS
species IDs, rate-provider ID, units/domain, solar or plasma input identity,
source files/hashes, named kinematic policy and declared bath transfers.
Validate before particle updates: masses/charges, complete product support,
duplicate channels, rate signs/finiteness, domains and all required data.
Load immutable data transactionally; callbacks perform no file I/O.
Require nonnegative integer stoichiometric multiplicities for discrete event
channels, with positive coefficients for listed reactants/products. The shared
continuum registry accepts noninteger coefficients; its validation alone does
not prove physical particle-event stoichiometry.

Keep full binary reactant stoichiometry for electron impact and charge exchange
even when one reactant is an external bath. The shared unit-unimolecular index
cannot directly represent their complete chemistry; provide an explicit
tracked-parent execution index linked to the bath. Shared `ExcessEnergyJ`
validation rejects negative values, so signed endothermic Q/threshold metadata
needs a separate reviewed representation, not a value silently replaced by zero.

Add opt-in configuration through the existing CLI `.in` parser using
`#begin application: moon` / `#end application`, consistent with
`src/INPUT_FILE.md`. Preserve current keys and legacy defaults. Proposed new
configuration must explicitly select execution mode, photon table/activity,
kinematics/bath policy, enabled channels, time-step policy and protocol file.
Exact new key names should be frozen with parser tests before documentation
shows them as executable syntax. Missing mandatory new-mode inputs are errors;
absence of the new mode retains legacy behavior.

The protocol file binds the selected case matrix, immutable data identities,
species/state catalog, oracle versions, seeds, numerical/statistical design,
predecessors and resources. A proposed protocol is not a reviewed/frozen one.

### 6.2 One physical event and one conservative transaction

Represented parent number is `N_p = W_parent(node) * individual_correction`.
For constant competing hazards, use `Lambda=sum(lambda_c)`,
`P(event)=1-exp(-Lambda*dt)`, an exponential clock and conditional channel
probability `lambda_c/Lambda`. Use stable `expm1`/`log1p` forms and a documented
open-interval uniform convention. Zero total hazard is an exact no-op.

For varying rates, consume the optical-depth clock `H=-log(U)` along actual
trajectory segments; locate the event where the integrated hazard reaches H.
Carry remaining H, physical epoch and continuation interval through MPI
transfers/clones. Bracket visibility/driver discontinuities explicitly.
Here continuation/cloning means the same physical trajectory/packet before its
event. After a real reaction, each product's next Markov clock must be newly
initialized according to the approved joint stochastic model; never copy a
consumed/zero parent clock to daughters and trigger spurious immediate cascades.

One accepted event removes/transforms the parent and creates the complete
correlated product packet. Product term j carries `nu_j*N_p`, with actual
reference weight converted to an individual correction. No product-dt/parent-dt
factor belongs in an instantaneous transient number transfer. Repeated H
products may require distinct correlated velocities; a single double-weight
H record does not automatically represent the same distribution.

Use plan → validate physics/capacity → reserve → commit. All parent/product
changes and diagnostics commit together. Capacity deferral leaves physical
state and budgets unchanged, preserves the pending event/random decisions,
and retries at a collective-safe checkpoint. Never redraw a channel because
its product packet was too large. Resource failure is an explicit ERROR,
not lost products, hidden merging or a changed reaction probability.

Install all upstream batching callbacks together:
`BeginBatchedModel`, `IsBatchedParticleEligible`, `GetBatchedParticleLimit`,
`BatchedReactionProcessor`, `FinishBatchedPass`. The old-list traversal skips
newly inserted records in the current pass, but later batch passes consult the
application eligibility callback. Moon must initialize daughter/processed fields
so those passes do not duplicate already consumed exposure. The Moon
continuation queue must handle permitted daughter cascades during their
remaining physical interval, or explicitly expose and converge a reviewed
splitting approximation. Do not call a delayed daughter reaction an exact
event solution.

### 6.3 Joint dispatcher, geometry and bookkeeping

One new-mode owner combines photon, electron and charge-exchange hazards.
The independent generic electron executor stays OFF, preventing duplicate
consumption of parents. Legacy mode keeps its existing selected dispatcher.
Tests must prove actual make-expanded/preprocessed macro selection, not merely
find a favorable token in generated text.

Gas photon exposure reuses the joint terrain/Earth ray masks with
`applyLocalSelfHorizon=false` and rays originating at the particle position.
Its scalar absorption visibility is `integral_visible I dOmega / integral_full
I dOmega`, under the explicitly selected angular-radiance/spectrum convention.
The existing `jointlyVisibleAreaFraction` is a projected-disc area fraction;
it is not that solid-angle ratio exactly. For the existing equal projected-area
grid and uniform radiance, use the `1/cos(rho)` solid-angle Jacobian per ray and
normalize by the full-disc weights. Expose a verified geometry accumulator;
do not reconstruct masks through a duplicate test-only ray service.
The present API also displaces the origin by a strictly positive surface
offset even with self-horizon OFF. Add an explicit volume-point origin policy
or verify an approved bounded-offset error; do not claim the current API
already launches exactly from an arbitrary gas particle. Preserve existing
surface self-intersection behavior and test near-boundary limits separately.

No facet cosine, projected-incidence factor or surface albedo belongs in the
gas ionization hazard. Bind the unobscured normalization to the selected D04
1-AU/r^-2 rate convention, resolve actual Sun distance/epoch and occultations,
and independently test distance and angular-quadrature errors. Any point-source
distance or uniform-disc approximation needs its named scope/error bound;
an exact finite-disc angular integral must not be mislabeled exact r^-2
normalization without that distinction. A particle inside solid Earth/Moon is
a boundary/geometry issue, not a reason to accelerate the photon rate artificially.

Electron and ion baths have independent geometry/driver rules; photon shadow
must not suppress an electron reaction unless the declared plasma input itself
has zero density. Transform parent, projectile and products into one documented
event frame. Use full SPICE state transformations when frame velocity terms
matter; a position rotation alone may not transform velocities correctly.

For ion trajectories, declare E/B and their frame/provider explicitly and
reuse the real Lorentz branch; add U06 and appropriate ion-mover prerequisites
to the affected linked wiring contract. A prescribed-bath coefficient does
not imply a consistent electric/magnetic field. Full M2 surface coupling also
requires an approved ion-impact/neutralization/trapping law: the current Na+/Ar+
adapter intentionally refuses these impacts. Preserve that guard; impact-free
chemistry fixtures have a limited scope and do not qualify surface transport.
An incoming thermal electron's reacting-energy distribution likewise cannot
be recovered uniquely from a total Maxwellian fit; explicit collision-energy/
angle sampling needs the corresponding collision-weighted differential data.

Record accepted events, parent loss, each product term, nuclei, charge,
represented model mass, momentum, kinetic/internal energy and signed bath
transfers. Keep surface loss, escape and numerical deletion distinct. Momentum
and energy across a reaction must include photons/projectiles/electrons under
the approved convention; external-force work over the step is a separate term.
Do not assert literal electron-inclusive rest-mass equality using equal nominal
neutral/ion masses. Reuse diagnostics without stealing existing observer slots.

### 6.4 Normal-mover integration and the limited potential core change

First prove transactions with the existing batched service in homogeneous
fixtures. Then connect event timing to the normal force-aware second-order
Moon tracer. A chemistry-only walker fixture cannot close this second gate.

Prefer extending the existing optional trajectory-event interface into that
tracer, preserving its default NULL behavior, geometry precedence and MPI/list
ownership. If its current single-particle result cannot return siblings, design
the smallest explicit pending-sibling/continuation contract. Do not duplicate
the force law or replace normal Moon with the first-order tracer merely because
that tracer already has hooks.

Every product originates at the physical event position/time, uses its own
post-reaction force/mover and advances exactly the remaining interval. Include
source birth offsets, surface/face/event ties, escape, remote MPI ownership and
daughter cascades. A species change cannot reset particle age or give a newborn
an unrelated full step. Unsupported asynchronous/local species clocks must be
refused in the initial common-step profile.

Reject an incompatible selected mover during initialization, even with zero
particles. The current all-or-none callback check does not prove that the
selected mover invokes the callbacks. No in-flight mode may advertise support
until its actual selected tracer and sibling continuation are qualified.

Before any necessary core edit, capture the passing baseline on the same
upstream/source/configuration; run targeted mover/particle/batching tests and
the full core routine afterward. Preserve every previous PASS and classify
known nonpasses separately. Do not change reference bytes automatically to
hide a changed trajectory. Keep application-specific rate/kinematic laws in
Moon. Core tests are a separate bounded work package if the combined one-hour
Moon development envelope cannot include them honestly.

## 7. Implementation sequence and concrete deliverables

These letters are work packages, not new permanent test IDs.

| Stage | Work and affected areas | Required reviewable result before moving on |
|---|---|---|
| A — freeze contracts | Registry/state inventory; .in/protocol specification; actual I21–I24 acceptance/README files; resolve section 5 dependency wording and per-channel evidence. | Complete independent oracles, case/seed matrix, hashes and failure semantics. Missing physics explicitly disables only its affected feature. |
| B — rates and transactions | D04 loader, existing D05 adapter, reviewed CEX rate kernels, upstream reaction-network/diagnostic reuse, correlated packet planner/committer. | Local U14–U17/U24 subchecks; exact nucleus/charge and approved mass/bath bookkeeping; no partial commits. Rate-only records do not activate transport. |
| C — mover wiring | Common physical clock, combined hazard owner, real force-aware event hook and sibling continuation if necessary. | Native tests through actual PIC stepping/mover, correct time/position/remaining interval, MPI ownership and source offsets; all prior core PASS retained if core changes. |
| D — atomic photon/electron transport | Opt-in D04 Na/He/Ne/Ar and approved ion/bath/recoil policy; actual Voronov particle conversion. | I21 and atomic I22/I23 components with independent survival/branch/budget oracles. Legacy M0/nightly still passes unchanged. |
| E — charge exchange | Approved channels, domains, conditional projectile sampler, state/Q/angle law and bath accounting. | U17/I24 direct rates, lifetime, actual products and kinematics; unsupported channels explicitly OFF. Missing differential data cannot be replaced by a total cross section. |
| F — molecular transport | Complete H2O/OH hazards, states, joint daughter distributions, cascades, state evolution and reviewed U24 chemistry dependency. | Full molecular I22 components; no omitted hazards, collapsed states or disappearing H/O/charge. Whole U24 retains its surface capability requirements. |
| G — convergence and closure | Three dt levels, weights/layouts/seeds, combined hazards, predecessor evidence, full registry aggregation and reproducible artifacts. | Applicable same-campaign M1/M2 gates and every declared I21–I24 criterion pass; source/recipe/data identities frozen and STATUS/READMEs reconcile. |

Rate loaders and event bookkeeping can progress while literature/data work for
E/F proceeds. D/E/F transport cannot be activated before its kinematic contract
is complete. Do not estimate a completion date until data decisions and a
measured cost pilot establish feasibility.

## 8. Verification contracts and statistical design

### 8.1 Shared fields inherited by each proposed I21–I24 specification

The following fields, together with each individual entry below, constitute
the proposed test specification. The implementation must copy them into each
runner's maintained acceptance contract rather than rely on this prose alone.

| Required field | Proposed common contract |
|---|---|
| Purpose / production capability | Distinguish scalar rate, event executor, force-aware wiring, convergence, conservation and any measured forcing; report separate subresults. |
| Production entry points | New Moon registry/dispatcher/committer, configured photolytic batch service, actual `PIC::TimeStep` and selected mover for wiring cases; existing `ElectronImpact` where applicable. Exact installed symbols recorded after implementation. |
| Configuration / inputs | Explicit opt-in .in; species/state/mass/charge/weight/common-dt configuration; process toggles; bath/spectrum/geometry; versioned pre-output protocol. Never a second test-only production law. |
| Independent reference/oracle | Closed-form survival/hazard/CDF and stoichiometry; independent high-precision coefficient/rate evaluator; independent reaction-chain matrix solution; reviewed product distributions and conservation identities. Regression output alone is not physical truth. |
| Units | Seconds, metres, m/s, m^-3, m^3/s, m^2, kg, coulomb, joule; explicitly marked source eV/keV/u/1-AU conversions. |
| Frame/time convention | Named inertial event frame and actual mover frame; explicit UTC input/SPICE ET conversion, physical epoch/step/birth offsets. Record each state transform and occultation policy. |
| Deterministic/stochastic controls | Frozen independent-parent trial unit, global cohort size, seed/stream policy, case membership, confidence/power and statistical method; OMP=1 and eight MPI by default, explicit one-rank control. |
| Procedure | Validate input/data/build and predecessors; run local oracles; execute real production path; preserve raw records; independently score; make comparison plots; aggregate without replacing child statuses. |
| Metrics / acceptance | Per-ID frozen numeric gates below plus reviewed state/budget/numerical contracts. Keep local stricter criteria. No tolerance chosen from observed residuals. |
| Expected PASS result | All required metrics actually evaluated and satisfied for the declared scope; no missing product/kinematic/predecessor component. |
| Expected failure modes | Wrong rate/unit/frame/branch, double executor, incorrect dt/weight, missing product/cascade, bad cell/ownership, nonclosing ledger, biased sampling, changed legacy regression. |
| Required data / hashes | Applicable D04/D05 and approved CEX/kinematics data; section 10 anchors and complete per-file manifest. D06 only where actual measured forcing is consumed/claimed. |
| Predecessor gates | Literal local U dependencies and formal M1/M2 qualification order; I21 explicitly requires R2.6. Each individual contract must enumerate its actual same-campaign gates. Supporting execution under the approved exception records missing predecessors and retains formal SKIPPED. |
| Artifacts | stdout/stderr, result/manifest/input/config/build identities, external hashes, seed/layout, independent references, raw science/budget records, CSV/PNG/EPS comparisons, resources and status reasons. |
| Status semantics | PASS evaluated/met; FAIL evaluated/not met; ERROR build/runtime/parser/corruption/infrastructure; SKIPPED missing capability/data/configuration/predecessor. Supporting PASS never overwrites formal SKIPPED. |

### 8.2 I21 — sodium photoionization and Na→Na+ conversion

Purpose: qualify the actual configured Na process, not just neutral loss.
Use D04 quiet/active rates, selected Sun distances including 1 and 2 AU,
explicit shadow/visibility cases and compiled Na/Na+ with the approved bath
and recoil policy. The independent homogeneous oracle is
`N(t)=N0 exp(-lambda*t)`; piecewise exposure uses the independently known
hazard integral. Zero illumination has exactly zero photon hazard.

Run the homogeneous production executor, then separate actual-mover cases
with known event clocks, acceleration/geometry, source birth offsets, face
crossing and remaining daughter interval. Controlled clock/bath fixtures are
verification inputs, not observations. A synthetic engine clock must not be
presented as a measured or physical D04 Na rate.

Metrics: **lifetime relative error <=2%; inverse-square rate error <=0.1%;
Na→Na+ species budget closes** under the frozen mass/charge/bath convention.
U14's legacy 130500-s test and `1e-14` criterion remain unchanged. Report new
mode U14 cases separately. Expected failure signatures include lost ions,
product-dt bias, shadow exposure after impact, duplicated chemistry and an
event at the wrong physical time. Keep survival/event-time, exposure and
parent/ion/bath ledger comparisons in CSV/PNG/EPS.

### 8.3 I22 — multispecies photon rates and branching

Purpose: extend the same production executor to He/Ne/Ar and the complete
approved H2O/OH reaction families. Include quiet/active endpoints, full product
states/multiplicities and subsequent OH chemistry during the correct remaining
interval. Independent rate/channel transcription comes from D04; independent
multinomial and reaction-chain/Bateman or matrix-exponential solutions test
frequencies and sequential products.

Metrics: **ensemble rate error <=2%; linked branch consistency at 95%
confidence; the complete declared branch family sums to one within numerical
precision**. U15 separately retains its **predeclared 99% intervals** and exact
transcription within storage precision. Arithmetic normalized conditional
probabilities over the complete rate family must not be confused with editing
rounded wavelength data to force partial/total equality.

Check all seven H2O/four OH channels, repeated H products, H/O/charge/state
bookkeeping and cascade timing. Missing daughter physics keeps its transport
subresult SKIPPED and full I22 unqualified. Atomic-only success is reported as
an atomic component. U24 chemistry subchecks cannot give aggregate U24 PASS
while its surface migration/sticking/residence requirements remain absent.
Preserve rate/branch, survival-chain and state-resolved ledger plots.

### 8.4 I23 — electron-impact coefficients and particle conversion

Purpose: qualify the real conservative electron reaction path using the
approved Maxwellian bath and existing four atomic Voronov fits. Input ne/Te
explicitly, including zero density and native temperature boundaries.
Independent high-precision/native-Fortran references test coefficients;
constant `lambda=ne*C(Te)` gives an analytical survival oracle.

Metrics: **linked coefficient relative error <=2%; ensemble lifetime error
<=3%; explicit zero/invalid density/temperature handling**. Existing U16
`<=1e-12` formula checks and domain refusals remain unchanged. Invalid values
must fail before state mutation; zero density must not create particles.

Check actual neutral/ion/bath transfer, non-unit individual corrections,
nonzero spatial cell indices, MPI movement, absence of photon-shadow suppression
and absence of the generic executor. A combined photon/electron fixture checks
sum-hazard survival and conditional process counts. State-specific electron
and ion recoil remains subject to the approved product model. Prescribed-bath
success is not R4.1/D06 measured-event qualification. Preserve coefficient,
survival, competing-process and bath-ledger comparisons.

### 8.5 I24 — documented charge exchange

Purpose: qualify an approved lunar neutral/projectile channel through the
actual executor, products and mover. Explicitly distinguish a monoenergetic
beam fixture from a thermal ion distribution or measured solar wind.
Use an independent source-fit/interpolation evaluator and independent
quadrature of `f*sigma*g`; sample collision partners with the corresponding
conditional density. Test valid endpoints, below/above-domain refusal and
species/state identity.

Metrics: **direct rate error <=1%; ensemble lifetime error <=3%; correct
products and exact charge/mass conservation under the reviewed convention;
unsupported channels explicitly OFF**. Add independently justified
momentum/energy/state/distribution checks before transport qualification;
their channel-specific limits must be frozen before activation. U17 already
requires fixed seeds and **at least 100000 events for probabilistic branches**;
retain that rule and freeze the authoritative per-channel criteria.

Total cross-section agreement or neutral disappearance alone does not pass
full I24. Missing H-state/Q/differential evidence remains SKIPPED for product
transport. Preserve rate/support, projectile-distribution, survival, product
phase-space and event-ledger comparisons. Do not hide out-of-range encounters
by clipping or renormalizing the projectile distribution.

### 8.6 Cross-cutting numerical, statistical and failure coverage

Use at least three dt levels at the same physical final time and a predeclared
numerical budget. Distinguish exact constant-rate event sampling from error
in trajectory/hazard integration. Test equal reference weights first, then
proposed stress ratios 1:8 and 8:1 plus fractional corrections, without changing
the represented physical population. These are numerical stress cases, not
assumed lunar particle weights. Unsupported unequal physical clocks get a
deliberate refusal test; they are not silently declared verified.

Eight-rank and explicit one-rank runs use the same global cohort size. No
bitwise cross-layout equality is assumed. Split records, sibling products,
repeated samples of one parent and copied RNG streams are not independent
additional trials. Freeze streams and inference for weighted/layout comparisons.

Exercise event/cell/surface ties, source injection offsets, remote exchange,
capacity deferral without random redraw, cascade queue exhaustion and malformed
data. Test the old generic cell-index failure mechanism by actual new-mode
product insertion at nonzero cell indices, not by enabling that old executor.
In a moving open domain, independently account for/censor escape and surface
events; they are not photoionization deaths in a survival estimator.

The existing percentage gates do not specify statistical power, branch
precision or a unique N. The recommended planning approach is to fix the
family and independent trials before outputs, then calculate exact binomial
tail risk and multinomial/branch intervals. Proposed family planning-risk
budget: 0.01, allocated over the frozen K stochastic comparisons. This is
**a recommendation requiring review, not a newly approved acceptance rule**.

For illustration, at fixed exposure `a=lambda*T=log(2)`, a constant-rate
equal-weight cohort has survival probability 1/2. For a lifetime point-error
gate delta, passing counts lie in
`ceil(N exp[-a/(1-delta)]) ... floor(N exp[-a/(1+delta)])`.
Exact binomial planning at risk `0.01/K` gives the following example global
parent counts, divisible by eight:

| Existing lifetime gate | K=1 example | K=48 example |
|---|---:|---:|
| 2% (I21) | 34496 | 71624 |
| 3% (I23/I24) | 15376 | 31880 |

These are analytical planning examples, not executed AMPS results, approved
case membership, branch precision guarantees or counts sufficient for a
different confidence-containment criterion. I22's rate estimator needs its
own corresponding design. Recompute rounding/tails at adequate precision
when freezing the actual protocol. The retained `sampling_design.py/.json`
documents the calculation.

Rare-branch detection also needs a pre-output event-count/power calculation;
observing a branch once is not precision. Do not select a wider interval after
failure. Keep U15's 99% and I22's 95% rules distinct, state pointwise versus
family coverage, and fix precision targets for small branches. A missing
statistical contract blocks qualification; an executed inadequate or failed
criterion cannot be relabeled PASS because runtime was expensive.

## 9. Runners, artifacts and runtime/resource policy

### 9.1 Global and individual ownership

Implement the existing I21–I24 as individual directories under
`srcMoon/test/integration-tests/`, each with README, versioned acceptance,
independent reference generator/data provenance, runner, scorer and plotting.
Use the repository's existing per-test result schema. Extend `run_all_test.py`
and the global `test_runner.py` registration/phase dispatch consistently, with
tests preventing omissions/duplicate IDs. Remove each planned placeholder
only when its real owner is installed; unavailable physics still reports
SKIPPED with the actual reason.

The global summary must preserve executed/PASS/FAIL/SKIPPED/ERROR counts;
for every FAIL/ERROR retain literal ID, individual runner and full log path.
Preserve linked children/artifacts and distinguish local/supporting/formal
outcomes. A protocol/preflight PASS must never become a linked PASS.

The current executable command remains **initial rate-only development**:

```sh
OMP_NUM_THREADS=1 /usr/bin/python3 -B srcMoon/test/test_runner.py \
  --phase m3-development --mpi-ranks 8 \
  --campaign-wall-seconds 3600 --cleanup-reserve-seconds 600 \
  --spice-path /home/vtenishe/SPICE/cspice \
  --spice-kernels /home/vtenishe/SPICE/Kernels \
  --output-dir test_output/srcMoon/M3-development-NEW
```

Replace `NEW` for each run and use an isolated checkout when preserving a
separately configured active tree. This command is documented in
[RUN.md](../RUN.md#initial-m3-atomic-rate-verification); it does not currently
execute completed I21–I24. Future full-chemistry/cluster CLI modes must be
implemented and documented before being presented as runnable commands.

### 9.2 One-hour development tier

Use one monotonic **3600-second** budget including input/data QA, removal of
generated `build/`, normal Config.pl/SPICE configuration, `make -j`, tests,
plots, cleanup and normal-configuration restoration. Keep the existing
600-second cleanup reservation and measured process-tree resource accounting.
Default linked Moon fixtures use eight MPI ranks/OMP=1; the reviewed nightly
continues to use its own four-rank contract.

First pilot the complete species recipe and fixture costs. Minimize redundant
builds by sharing an identity-verified M3 fixture among compatible cases.
Retain separate normal and chemistry-disabled M2 recipes; never reuse an M2
OFF binary for an M3 claim. The builder must preserve immutable generated
snapshots, executable/source/input/config hashes and actual compiler-macro
proof before/after build, including both application and PIC chemistry controls.

The previous initial M3 package took 16.3 minutes, corrected M2 package 51.2
minutes including diagnosis/restoration, and full CORE about 40.8 minutes in
separate campaigns. They cannot simply be summed and called a one-hour test.
Choose separately bounded work packages, measure the new cost and defer full
high-statistics qualification when it does not fit. Do not lower physical
rates, observation requirements, confidence or dt criteria to fit the clock.

### 9.3 Explicit cluster qualification tier

Only after reviewed physics/data/statistics and a measured pilot, prepare a
versioned execution contract with total CPU-core-hour/wall/allocation budget,
per-case costs, rank layout, memory/disk limits and retry/cleanup policy.
The user allows longer cases of a few hundred CPU hours; the exact M3 case
matrix/allocation remains to be frozen from measurements. A larger allocation
cannot repair missing kinematic data or failed predecessors.

Long execution requires explicit cluster opt-in; development commands must
not launch it automatically. Preserve the physical horizon, cohort/statistical
design and numerical refinement. Where affordable, homogeneous exact-rate
ensembles avoid unnecessary transport cost; retain separate real-mover gates
so this optimization does not bypass production wiring.

### 9.4 Resources and reproducibility

For local runs retain the user's **200 GB total RAM ceiling** on the 256-GB
machine; account for all ranks/build children and declare GB/GiB units. Set
disk quota/free-space floor and allocation limits explicitly. Monitor event
queues, particle capacity and actual process memory; on exhaustion terminate
with production `exit(__LINE__,__FILE__,msg)` explaining resource, measured
usage, reason and concrete ways to reduce resources. Preserve ERROR artifacts;
do not silently truncate products, time, trials or iterations.

Mesh/restart output stays opt-in for chemistry verification. Emit compact raw
science/ledger tables, not redundant full-domain dumps. Keep scientifically
meaningful prior results. Each linked case preserves stdout/stderr,
`result.json`, `run_manifest.json`, input/config snapshots, source revision,
compiler/MPI/OMP/build/executable identities, generated macros, seed/streams,
external hashes, science metrics and resource history. Comparison products
are CSV, PNG and EPS, explicitly labelled verification or qualification scope.
These are analytical/production comparisons, not observational validation.

## 10. Required data identity anchors

Paths below are relative to
`/data/vtenishe/moon_validation_data/photochemistry/`. The retained
`data_audit.json` records all 59 inspected identities, expected/actual hashes,
native labels, inventory and D06 gaps. Every future run must recheck its exact
consumed files and package readiness rather than trust this dated snapshot.

| Input | SHA-256 |
|---|---|
| `processed/combined/lunar_photochemistry_catalog.csv` | `25166581488316b5299fe40384b52c848ea96162d8579f8125c6d917cbd99d2b` |
| `processed/phidrates/photo_rates_1au.csv` | `3aa54915a08258371dbe975949e6afb2ab74618164fadca8bceebe876c926847` |
| `processed/phidrates/reaction_channels.csv` | `2d9fa1c8eb72d3bbba0b91ab4c9a912e6ebce2abc80b251096df2f41a4354061` |
| `processed/phidrates/cross_sections.csv` | `e900eb5514b616dcbbde7e3ca124b5c06d4d22a28a3cd429ed58a14c1933c860` |
| `processed/phidrates/branching_ratios.csv` | `430bb68558a43973f0622094049c87004dd522a6288a48935027856f48229ba7` |
| `processed/phidrates/thresholds.csv` | `a3335b3b37c3d275a4f72e92ae8e5b1f6ecfcbc4e2d84735cb20694dc93b48f1` |
| `raw/electron_impact/voronov/cfit.dat` | `470ea3078c64159eece6a469314dde2356e2d6ccaf618f4eaecb46649d61a5cd` |
| `raw/electron_impact/voronov/cfit.txt` | `5b12dc3b4b546ef863216d57df7d64c81a553cc5eb32a8cf2d3d275d4158f9ef` |
| `raw/electron_impact/voronov/cfit.f` | `7fe7d78c25166084041bea58bb55c2fcc6ad724fba8432c167cf044e6eb1bcc1` |

CEX publication identities: ORNL PDF
`3bb3fdb4be8f0f806d378083be8798cf3ce347714c3d06d58f90ee483b10c030`;
NIST PDF `7c412cb93f827b3412b6da000e674c63237cacb0e201733a6f5a4f431359952e`.
These identify publications, not a completed processed kinetic registry.
D06 `drivers/qa_report.json` SHA-256 is
`589342d3054c080359567290591e786adfe3e464bb1105bb3b838c6a75a8a2aa`;
its intact ERROR report does not establish readiness.

## 11. Definition of done and proposal work-package report

M3 may be declared complete only when:

- Reviewed scope/dependency, species/state/mass/bath, kinematic, statistical
  and numerical contracts are explicit; both roadmap copies are reconciled
  through an authorized update where needed.
- The new opt-in registry and conservative combined executor run through normal
  initialization, sources, geometry, chemistry and force-aware mover; legacy
  modes and historical core PASS cases remain preserved.
- All selected channels have authoritative rates, complete products and
  qualified kinetic models over the actually exercised domain; unsupported
  channels remain OFF and no original hazard is silently removed.
- Applicable local gates and same-campaign M1/M2 predecessors pass; I21–I24
  actually evaluate every declared criterion. U24's unresolved full-scope
  requirement and measured-forcing claims are resolved honestly.
- Counts, conservation, event timing, dt/weight/layout/statistical qualification
  and per-channel physical oracles pass under frozen criteria; failures and
  absent evidence are retained with their proper statuses.
- Individual/global runners, references, CSV/PNG/EPS comparisons, manifests,
  resource contracts, detailed code comments, READMEs, RUN.md and STATUS.md
  are complete and reproducible. No observation is fabricated or tuned.

For this proposal only: source/roadmaps, current production dispatch, upstream
reaction/mover APIs, test registries and D04/D05/D06/CEX evidence were inspected
with three parallel read-only reviews. Documentation and retained planning
artifacts are the only changes. **No configuration/build/model test was run;
there are zero new U/I PASS/FAIL/ERROR/SKIPPED evaluations.** Existing formal
SKIPPED results remain unchanged. Identity verification is 59 matches and zero
mismatches, not 59 science tests. Documentation consistency and protected-file
checks are recorded in the retained proposal completion report.

Build discipline for subsequent implementation remains: remove generated
`build/` before each required normal configuration/rebuild, use Config.pl and
`make -j`, preserve build logs and require `make -j test_Moon` to succeed with
a zero-byte `test_Moon.diff`. No new regression reference is frozen by this
proposal. No commit, push or staging operation is part of this work package.
