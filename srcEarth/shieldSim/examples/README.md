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
