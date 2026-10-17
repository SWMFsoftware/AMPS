# Versioned case examples

`stage1_case.json` is the smallest complete science-facing Stage-1 input. It
defines every world, source-plane, shield, placement, and absorber dimension;
uses two ordered shield layers; and requests two **independent** slab absorber
cases. Running it creates separate output prefixes ending in `water_proxy` and
`silicon_proxy`. The absorbers are never placed behind one another.

Transverse centers are explicit even when they are zero. Omitting a
`center_x`, `center_y`, front-z, size, or thickness field is an error; the
versioned parser will not infer a centered placement or legacy spacing.

Validate and archive the resolved case without transport:

```bash
./shieldSim --validate-case=../examples/stage1_case.json
```

Run the case:

```bash
./shieldSim --case=../examples/stage1_case.json
```

The omitted global accuracy target in this example resolves to exactly `0.10`.
That resolution is visible in each `_resolved_case.json`. Because Stage 1 does
not yet have a complete systematic-uncertainty budget, its default `combined`
accuracy result is `not_assessable`; the adaptive sampler still tests and
reports Monte Carlo precision and stops at the resource limit if necessary.

`sep_spectrum.dat` remains a synthetic format example, not a published SEP
environment. Its mandatory metadata is `sep_spectrum.metadata.v1.json`.

## Downloaded-spectrum Stage-2 validation case

`stage2_ams02_proton_validation_case.json` exercises the same case parser,
primary generator, Geant4 transport, and TID/H100/10 writers that a later AMPS
adapter will use:

```bash
./shieldSim --case=../examples/stage2_ams02_proton_validation_case.json
```

The case selects only the checksum-pinned AMS-02 BR2460 A=1 proton component.
It uses `purpose: verification_test` and an isotropic boundary because those
are enforced admission requirements. Published differential intensities are
sampled as piecewise-constant values over their exact transformed energy-bin
edges. For a proton, MeV/n is numerically equal to total kinetic MeV.

This file is a reproducible executable test, not a preliminary-mission source:
one near-Earth Bartels rotation is not a complete deep-space GCR environment.
Changing the purpose to a mission assessment, choosing inclusive helium, or
changing the source to `beam` is rejected. The alpha energy-limit block remains
required by the v1 case schema but is unused because the admitted component is
proton-only.

## Stage-3 monoenergetic response-cell case

`stage3_response_cell_case.json` demonstrates the verification-only
`shieldSim.case.v2` contract:

```bash
./shieldSim --validate-case=../examples/stage3_response_cell_case.json
./shieldSim --case=../examples/stage3_response_cell_case.json
```

The file uses explicit source-plane, shield-reference-surface, and
absorber-reference-surface transforms. Translations carry units. The
`rotation_quaternion_local_to_world` object uses scalar-first `(w,x,y,z)`, is
unit-normalized, and actively rotates local axes into the world frame. For a
slab, local `+z` points downstream; the absorber transform origin is the center
of its upstream face rather than its volume center.

The source creates one 100 MeV proton at `mu=1`, `phi=0` per history. The
two-layer Al/HDPE stack and silicon slab produce a
`*_response_cell.json` containing response moments per generated primary,
primary/secondary ancestry, per-layer production diagnostics, and a complete
response-cache key. This example requests only `transport_raw`; case v2 also
permits the separately gated Stage-4 finite-slab `TID` output. Raw tallies are
not automatically validated TID, DDD, `n_eq`, LET, hardness, device, organ, or
instrument outputs.

To create additional grid cells, copy the case and change only the declared
`incident_id`, species, energy, and/or direction. Do not change geometry,
materials, physics, cuts, step limits, build, or Geant4 data when assembling
the cells into one matrix; the response tool will reject any such key change.
Protons use `name: proton, Z: 1, A: 1`, He-4 uses
`name: alpha, Z: 2, A: 4`, and another isotope uses `name: ion` plus explicit
`Z,A`. Both total kinetic energy and kinetic energy per nucleon are supported
when their unit is declared.

Case v2 also accepts `name: gamma, Z: 0, A: 0` with total kinetic energy for
the Stage-4 photon attenuation checks. Gamma TID is not in the validated
domain. Every case-v2 transport block explicitly declares `em_physics` as
`reference_list_default` or `G4EmStandardPhysics_option4`; this choice changes
the case and response hashes.

A direct cross-check replaces the monoenergetic spectrum object with
`format: discrete_response_nodes` and a normalized list of exact node
definitions plus weights. The Stage-3 CTest driver generates examples of this
form under `test_output/stage3_transport_runs/`; generated test cases do not
belong in this source directory.

## Stage-4 finite-slab TID case

`stage4_tid_case.json` exercises the physically defined Stage-4 slab observable:

```bash
./shieldSim --validate-case=../examples/stage4_tid_case.json
./shieldSim --case=../examples/stage4_tid_case.json
```

TID is the Geant4 ionizing-deposition ledger divided by the full finite silicon
slab mass. Total and nonionizing ledgers remain audit fields. The result also
contains the exact rule that the numeric rad value is 100 times the numeric Gy
value, particle and propagated-
ancestry deposition partitions, and separately labeled creator/step-process
diagnostics. It is material-slab dose per incident primary, not device, organ,
or crew dose.

The Stage-4 gate accepts only proton or He-4 incidence and absorber material
identities for aluminum, silicon, or water. Gamma and explicit heavier-ion
nodes remain available for `transport_raw` verification, and other approved
materials remain available in the shielding stack, but the parser rejects a
TID claim outside the validated incident/absorber matrix.

The example demonstrates the complete structured uncertainty input. Four
named categories—material/mass, source, numerical, and physics reference—carry
provenance. All six unordered correlations are explicit and the parser checks
that their matrix is positive semidefinite before calculating the combined
standard uncertainty. The zero mass/source values are valid only for this
ideal verification fixture; mission cases must provide actual fabrication,
environment, numerical, and physics-reference uncertainty evidence. Omitting a
complete budget does not prevent a diagnostic run, but combined accuracy is
then `not_assessable`.

The example status is “implemented and acceptance-tested, pending independent
review.” Passing it does not itself enable the separately gated Stage-5 damage
calculation, LET, `H100/10`, or any health or device-response claim.

## Stage-5 silicon DDD and neutron-only `n_eq` case

`stage5_ddd_neq_case.json` exercises the authenticated absorber-entry damage
path:

```bash
./shieldSim --validate-case=../examples/stage5_ddd_neq_case.json
./shieldSim --case=../examples/stage5_ddd_neq_case.json
```

The source is a 1.05 MeV neutron directed normally through a negligible but
nonzero aluminum layer into a 50 mm × 50 mm natural-silicon slab. This exact
energy lies in the ASTM E722 `[1.0, 1.1)` MeV group, so it is useful as a
transparent reference example. Each accepted entry contributes the group NIEL
divided by 25 cm² to DDD and the group hardness factor divided by 25 cm² to
`n_eq`. The output records the exact table package/hash and states that charged-
particle DDD is excluded from `n_eq`.

The `damage_response` block is mandatory for a DDD/`n_eq` case and identifies
the strict manifest in `data/niel/`. Changing the absorber to a non-silicon
material, changing the declared library identity/version, or using an
unsupported incident species is rejected at case load. A transported secondary
outside a response domain instead makes the affected runtime result
`not_assessable`, with the rejected crossings and supported partial sum
retained.

This example requests statistical accuracy only. For a combined-accuracy
request, the case must supply a `physics_reference` uncertainty component and
its provenance; shieldSim does not assign an unsupported uncertainty to the
published damage function. The Stage-5 capability is implemented and
reference-tested, but formal completion remains pending project-selected
accelerator/device evidence and independent review. The result is a silicon
exposure proxy, not a device-degradation prediction.

## Stage-6 LET and H100/10 case

`stage6_let_hardness_case.json` demonstrates both Stage-6 outputs with an exact
two-node proton source:

```bash
./shieldSim --validate-case=../examples/stage6_let_hardness_case.json
./shieldSim --case=../examples/stage6_let_hardness_case.json
```

The 50 and 150 MeV nodes have equal probability and normal incidence. In the
negligible-loss limit both exceed 10 MeV and only the second exceeds 100 MeV,
so the source-level analytic hardness is 0.5. The actual result remains a
transport score at the shield rear face and reports event covariance plus a
95% Fieller interval.

The required `let_response` selects `shieldSim.electronic_let.v1`, explicit
0.001--10000 keV/um logarithmic edges at 20 bins/decade, and the mandatory
`not_assessable_no_clamping` policy. The resolved response and every edge enter
the case hash. The result reports aggregate and component track-fluence bins,
track/dose means, both linear and mass stopping units, delta-electron components,
and under/overflow and zero-length-step audits.

This proton case uses the reference-list EM baseline. A C-12 LET case must use
natural silicon and `G4EmStandardPhysics_option4`; the parser rejects the
reference-list default because it fails the fixed 10 MeV/u ATIMA screen. The
example requests statistical accuracy only. A combined claim requires a
complete structured uncertainty budget with a provenance-bearing
`physics_reference` component. Stage 6 is implemented and reference-tested but
remains pending independent review; this example does not validate CRaTER,
instrument response, device response, organ dose, or crew risk.
