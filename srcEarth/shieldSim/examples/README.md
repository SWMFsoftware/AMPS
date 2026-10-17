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
response-cache key. The only permitted requested quantity is `transport_raw`.
These are verification tallies and are not validated TID, DDD, `n_eq`, LET,
hardness, device, organ, or instrument outputs.

To create additional grid cells, copy the case and change only the declared
`incident_id`, species, energy, and/or direction. Do not change geometry,
materials, physics, cuts, step limits, build, or Geant4 data when assembling
the cells into one matrix; the response tool will reject any such key change.
Protons use `name: proton, Z: 1, A: 1`, He-4 uses
`name: alpha, Z: 2, A: 4`, and another isotope uses `name: ion` plus explicit
`Z,A`. Both total kinetic energy and kinetic energy per nucleon are supported
when their unit is declared.

A direct cross-check replaces the monoenergetic spectrum object with
`format: discrete_response_nodes` and a normalized list of exact node
definitions plus weights. The Stage-3 CTest driver generates examples of this
form under `test_output/stage3_transport_runs/`; generated test cases do not
belong in this source directory.
