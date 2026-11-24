# P02–P04 reduced kinematics: analysis and test contract

Written before source implementation,2026-11-18, against revision
`74974c2b2fc1b4cf09a73c320064a22370c0ca94` and the user's explicit approval of
P01–P09. The decision documents still carried their older proposed status
when this work began; the current user approval controls implementation, not
automatic physical activation. These are existing-meaning M3 supporting
controls, not new permanent U/I identifiers or formal chemistry qualification.

## Inventory, authorities and production boundary

`species/ReactionExecution` provides catalog/state identities, exact material
stoichiometry, physical intervals and immutable conservative plans; it does
not supply source-dependent product kinematics. Existing Hartley2 recoil
sampling removes COM velocity and rescales equal-variance velocities. That
is a conservative isotropic two-body closure but generally not uniform
mass-weighted microcanonical sampling for unequal three-body masses.
`ReducedKinematics` will contain the independently reusable numerical service;
it will not change that application, PIC, a parser or the existing planner.

Source/external roadmap hashes are respectively
`712e1152bf50af6d9ee1dad8f8a585809ca8b92cf42ac6fc7a009f6c0253ee17` and
`d10eca1c539f9cf82940236d9e1aa36af665d43e82ba475925950d582fbf49c4`.
The source adds approved resource/development-order amendments; external
native metadata governs real data and common formal predecessors remain.
No isotope mass, threshold, internal population, excitation or source-state
mapping is inferred or acquired by this numerical work.

## Explicit bookkeeping and frame contract

Reference mass is the supplied neutral-isotope sum minus integer charge
times supplied electron mass. Every isotope/state, formation/excitation
energy and channel energy carries an explicit nonempty source identity.
Formation Q is incoming minus outgoing formation PLUS excitation energy,
with explicit occurrences; photons and other reservoirs enter separately.
No rounded heavy physical rest-energy subtraction or invented mass defect
is used. An independent helper reports reference/physical inertial-mass,
force-acceleration and relativistic-kinetic errors against explicit limits.

The initial photon service is **parent-rest inertial frame only**. Nonzero
parent momentum is refused, even if a caller labels that frame differently.
The parent approved this bounded API. A moving-parent/LSO provider still needs
consistent particle/photon state transformation; this service supplies none.
Representative photon energy is explicit threshold plus matching source mean
excess. Threshold must agree with the declared state-energy difference under
the caller's numerical crosswalk bound. Its momentum is `E_gamma/c * n`.
Available COM translation subtracts the product-state energy and photon
recoil COM energy; a negative value is refused, never clipped to zero.

A generic two-body final-state helper additionally accepts explicit total
momentum and available COM kinetic energy, product states and a unit relative
direction. This can serve a separately reviewed CEX provider; it chooses no
CEX source, Q, rate, direction or bath law. Photon two-body direction comes
from two supplied open uniforms, not an internal RNG.

Three-neutral sampling takes exactly three resolved states and nine supplied
finite standard normals. In mass-weighted momentum coordinates it projects
out total momentum and rescales to supplied COM energy. It preserves all
three occurrences, including identical H labels. Zero COM energy gives a
deterministic common velocity; a degenerate draw with positive energy is
refused without resampling. Numerical projection uses the heaviest product
as the dependent momentum to avoid cancellation for disparate masses; it is
the same constrained subspace, not a physical energy-partition adjustment.

## Test specification

| Field | Contract |
|---|---|
| ID/name | Supporting P02/P03/P04 reduced-model kernels under existing M3 meanings; no new ID. |
| Purpose/capability | Source-explicit reference-energy bookkeeping and conservative two-/three-body final states. |
| Production entries | `ReducedKinematics::{ResolveState,FormationQ,EvaluateModelErrors,GeneralTwoBodyFinalState,PhotonTwoBody,ThreeNeutralMicrocanonical}`. |
| Configuration | Pure C++17; explicit numerical mass/energy/momentum/unit-vector bounds and maximum classical speed; explicit independent model-error limits. No default physical tolerance, RNG or PIC allocation. |
| Inputs | Caller-supplied isotope/electron masses and source strings; state formation/excitation values; source threshold/mean; complete product occurrences; frame, normals/uniforms, total momentum/COM energy. |
| Independent oracle | Closed-form reduced-mass two-body identities and independent high-precision Decimal arithmetic; three-body Jacobi-coordinate kinetic identity and Beta(3/2,3/2) first/second moments/distribution; exact nuclei/charge/reference-mass identities. |
| Units | kg, J, kg m/s, m/s; dimensionless supplied random variates and relative bounds. |
| Frame/time | Explicit caller inertial frame; photon is parent-rest. No implicit SPICE epoch, LSO boost or moving-surface transformation. Instantaneous event only. |
| Deterministic/stochastic controls | No RNG inside production kernels. Deterministic analytical cases and a fixed-seed external normal ensemble with recorded count; no fabricated observations. |
| Procedure | Compile maintained production code with thin native adapter, compare actual results with independent references, exercise negative/refusal/unchanged-output cases, retain complete build/source/fixture/metrics evidence. |
| Metrics | Exact source/state/occurrence identities; mass, charge, Q, photon momentum, product COM/relative energy and momentum residuals; error-bound reports; Jacobi energy fraction and moments. |
| Acceptance | Explicit source/identity/refusal decisions exact; deterministic arithmetic reference errors within128 native long-double eps relative when conditioned; production closure passes ONLY caller-declared numerical bounds. Numerical tests use predeclared rounding bounds, never a substituted physical error budget. Microcanonical mean/second moment use frozen99% simultaneous concentration bounds, and symmetry/isotropy metrics have separately predeclared statistical limits. No threshold changes after observing outcomes. |
| Expected PASS | Declared local mathematical/transactional criteria met. Outputs/provenance stable; no `KinematicPolicy::qualified` flag set or guard bypass. |
| Failure modes | Unknown/duplicate isotope, absent source/state/energy, invalid charge/mass/threshold, nonunit direction, nonopen uniform, negative COM energy, invalid/degenerate normals, nonclassical speed, unrepresentable result, incompatible closure or model-error bounds. Failure preserves output. |
| Required data/hashes | No observational or canonical data used. Analytical fixture values explicitly labelled. Source, compiler, executable, numerical controls and source-roadmap identities retained. |
| Predecessors | Local work allowed by approved development exception; actual physical catalog/source-state crosswalk, native adapter/controller and unchanged formal gates remain separate. |
| Artifacts | `result.json`, `run_manifest.json`, compile/stdout logs, fixture/metrics CSV and mathematical-reference PNG/EPS; complete initial failure artifacts retained. |
| Status | PASS local criterion evaluated; FAIL evaluated mismatch; ERROR build/runtime/corrupt artifact; SKIPPED missing compiler/optional dependency. Reduced physical models remain EXPERIMENTAL. |

The fixed external ensemble has40000 events, seed731993 and17 bounded
statistics: Jacobi energy mean/second moment, three single-product energy
shares, three squared directional cosines, and nine CDF values at0.1,...,0.9.
Use the99% simultaneous Hoeffding bound
`sqrt(log(2*17/0.01)/(2*40000))` for each statistic. The analytical targets
are Beta(3/2,3/2), `E[x]=1/2`, `E[x^2]=5/16`, energy shares
`(1-m_i/M)/2`, and directional cosine squares1/3. These tests verify the
selected phase-space measure, not its accuracy for molecular breakup.

The development package is bounded to3600s including implementation, compile,
tests, figures and reporting. No active application build, core edit, canonical
data mutation, parser activation, commit or push is part of this work.

## Independent review corrections before source freeze

The review found three numerical seams, reproduced separately before fixing:
an excessively wide caller norm bound could admit a zero photon direction;
collapsing formation plus excitation before signed subtraction lost a1J Q
beside1e20J; and reconstructing polar sine from a rounded cosine erased the
representable transverse component at open uniform1e-25. Refuse zero vectors
unconditionally, normalize accepted near-unit encodings, sum signed energy
components through a floating expansion while retaining original occurrence
records, and use `2*sqrt(u)*sqrt(1-u)` for sine. The original criteria remain;
new explicit source/Decimal regressions exercise each mechanism.

### I22 explicit charged-breakup construction

Use one supplied fraction of the post-state, post-COM kinetic energy and two
supplied directions in the same inertial axes. The first direction gives the
outgoing electron relative to the heavy-pair center of mass. The second gives
the ion relative to its neutral partner. The event has one free electron, one
singly charged ion and one neutral with explicit nuclear composition. Each
state uses the same resolved mass and formation-energy reference. The caller
supplies the available translation after accounting for all chosen internal
energies, the incoming photon and common COM motion. Each fraction endpoint
is valid, and all nonzero relative modes require a supplied unit direction.
This algebra alone assigns no population, energy or angular distribution.

The predeclared deterministic controls use two noncollinear directions,
unequal masses, zero and nonzero bulk velocities, both partition endpoints,
zero available energy, invalid partitions, missing source provenance and
wrong product classes. Reconstruct both relative-mode energies and directions,
each product momentum, the total momentum and the total kinetic energy.
Compare against the independent Decimal two-body boost identities under
the existing 128-machine-epsilon arithmetic rule. The existing supporting
reduced-kinematics runner records its CSV, manifest and individual statuses.
No new observational data or stochastic physical law enters this control.

The parent-rest photon wrapper adds a resolved parent and independently sourced
threshold and mean excess energy. Verify the photon energy and momentum,
formation Q and all three material occurrences using the same additive mass
reference. Require exact isotope/charge conservation, refuse a moving parent
in a parent-rest call and reject an inconsistent source threshold under the
strict source-equals-state policy.
