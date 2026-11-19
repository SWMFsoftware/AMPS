# Source-explicit physical event composition — implementation contract

Written before implementation, 2026-11-18, revision
`74974c2b2fc1b4cf09a73c320064a22370c0ca94` plus the approved working tree.
The user approved P01–P09 and explicitly approved the bounded experimental
moving-parent Na extension: consistent incoming photon transformation and
complete electron/radiation ledgers, with source/state crosswalk and model-error
bounds required. Approval is not missing scientific data or M3 qualification.

## Inventory and ownership

`Photochemistry::Registry` supplies the immutable complete D04 source family,
channel multiplicities, quiet/active rates and matching mean excess energies.
`ReducedKinematics` supplies resolved additive-reference states, formation Q,
error analysis and the conservative two-body solver. `ReactionExecution`
supplies exact integer topology, source/product state identities, represented
number, immutable remaining intervals, fresh clocks and conservative plans.
The existing controller/PIC adapter owns event timing, RNG, allocation and
publication. This new PIC-independent module composes those same services;
it introduces no alternative rate law, mover, event scheduler or surface law.

The first scope is the single complete Na photoionization channel, one tracked
Na+ and an outgoing reservoir electron, with the incoming photon explicitly
removed from the radiation reservoir. Other source families are refused, not
partially sampled or renormalized. Sodium charge exchange is not substituted
for photoionization. The existing ion collector may receive a state converted
from this same explicit energy/catalog binding; it remains an artificial
absorbing boundary, never neutralization, heat, implantation or release.

## Inputs, frame and energy contract

Every state, isotope/electron mass, formation/excitation energy, physical rest
mass, source threshold, unit conversion, numerical bound and model-error limit
is explicitly supplied and source identified. Catalog masses must equal the
binary64 representation of the resolved additive-reference masses; the actual
PIC adapter still checks the generated mass/charge table independently.

The representative photon energy in the parent rest frame is the selected
D04 threshold plus the SAME quiet/active source mean. An actual inertial ray
and parent velocity determine the inertial photon energy through the exact
photon Doppler transformation. All material kinematics remain the approved
classical reference-mass model. The complete inertial input momentum/energy
and source-derived Q determine the available product COM translation. The
same production two-body solver receives an isotropic direction from two
caller-supplied open uniforms. Negative energy, unresolved metadata and
unrepresentable values are refused without clipping or resampling.

The supplied affine instantaneous velocity transform is
`v_event = R * v_inertial + offset`. It has a separately declared physical
event epoch and geometry/frame basis epoch, preserving the frozen-environment
meaning. R must be a proper orthogonal rotation within explicit numerical
bounds. All material particles share the same event position and transform.
The signed reservoir moments are calculated from the outgoing electron and
incoming photon independently of any residual. Their event-coordinate moments
follow the same classical affine transformation, including electron mass,
velocity-offset work and kinetic-energy terms. The receipt retains physical
inertial photon/electron components separately; coordinate energy is not a
claim about a photon frequency in a rotating frame.

## Admission and tests

A planner authorization template is an explicit caller input. This module
never grants physical qualification from a boolean or evidence string. The
runtime owner must first verify a reviewed source-bound state mapping and
select the experimental model; no application mode is enabled here. Missing
Na ion/threshold mapping is currently a real physical blocker: the D04
threshold and NIST ASD gap differ and no authoritative interpretation has yet
been established. Widening a numerical threshold tolerance is forbidden.

The source/state admission check now enforces this independently of every
caller tolerance. It compares the exact signed sum of all native state
formation/excitation components with the native source eV value times the
exact SI rational `1602176634/10^28 J/eV`, using the existing arbitrary-size
integer engine. Only two native long-double spacings at the small source
threshold cover coefficient/product/state representation. This is neither a
published source uncertainty nor a configurable acceptance threshold. A large
shared energy zero cannot enlarge that envelope. A regression deliberately
sets the caller threshold bound to1J and verifies that a1e-24J physical gap
disagreement is still refused without replacing the prepared provider.

Required local controls use the actual production composer and an independent
closed-form/high-precision photon/two-body/affine-moment oracle. Include zero
and moving parents, nonzero frame offset, both source endpoints, exact nuclei
and signed charge/mass, complete photon/electron moments, unequal statistical
weights, original remaining interval, fresh clocks, degenerate/nonopen draws,
unsupported families, wrong source/state/threshold, missing physical masses or
error budgets, unqualified authorization, invalid rotations and unchanged
outputs after every refusal. Analytical input states are explicitly numerical
fixtures, never lunar states, observations or evidence that the real source
crosswalk is ready. Preserve source/compiler/executable identities and logs.
Local PASS is numerical/software verification; formal I21–I24 and all preceding
qualification gates remain unchanged. No new permanent U/I ID is assigned.

The source/canonical roadmap hashes remain respectively
`712e1152bf50af6d9ee1dad8f8a585809ca8b92cf42ac6fc7a009f6c0253ee17` and
`d10eca1c539f9cf82940236d9e1aa36af665d43e82ba475925950d582fbf49c4`.
Approved development amendments control implementation order; unchanged formal
gates and native source metadata govern qualification and interpretation.

## Additive explicit source-effective policy — 2026-11-19

The earlier exact source/state admission rule above remains the **default**.
The separately named experimental policy and authenticated eight-column mapping
are specified in [M3_SOURCE_EFFECTIVE_PHOTON_CONTRACT.md](../M3_SOURCE_EFFECTIVE_PHOTON_CONTRACT.md).
They account for the resolved physical formation gap while preserving the source
photon threshold/mean and hazard. This explicit model convention does not widen
the default representation allowance or any numerical acceptance criterion.
Missing/wrong state-source identities, unsupported mappings and negative energy
remain refusals. The retained real-source Na profile resolves the opt-in atomic
mapping; physical model-error admission and linked/predecessor gates remain.
