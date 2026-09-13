# Step 5 change manifest: field-line-only srcSEP

Step 5 removes full three-dimensional particle transport from `srcSEP` while
retaining magnetic field lines embedded in three-dimensional space, imported
three-component background vectors, focusing geometry, and field-line observer
sampling. The separate Cartesian-transport application was not supplied in
this source package, so this manifest records the transfer boundary explicitly.

## Removed production implementations and state

- Removed the Parker3D mean-free-path, 2019 He, Kartavykh, Borovikov, default
  Boris dispatcher, and drift mover implementations and declarations.
- Removed `drift.cpp` and its magnetic-gradient/curl cell-data offsets,
  derivative setup, and static-cell sampling hook.
- Removed Cartesian trajectory selectors, mover-selection macros, the
  Parker3D Runge–Kutta selector, perpendicular-displacement state, and the
  Cartesian momentum extensions used only by transferred movers.
- Removed all AMR-node particle-list alternatives from srcSEP movers and
  sampling. `production_mover_runtime.cpp` now rejects any build not configured
  for field-line mode and field-line-segment particle attachment.
- Removed the spatial-neighborhood `sample3d` module and its three object files.
  This intentionally removes volume-neighborhood population, flux, return-flux,
  and pitch-angle output. Field-line observer sampling remains available via
  `SEP::Sampling::InitSingleFieldLineSampling` and `cSamplingBuffer`; recreating
  the removed volume diagnostic belongs in the separate 3-D application.
- Removed the AMR-cell drift diagnostic in `output.cpp`. The retained sampling
  reports field-line density, flux, return flux, pitch angle, energy, Larmor
  radius, and mean free path.

## Retained three-dimensional model information

- Field-line vertices and segments retain three Cartesian coordinates.
- Magnetic field, plasma velocity, electric field, and SWMF/SWCME background
  vectors remain three-component data on the field line.
- `FluxTubeGeometry` continues to derive SI area and segment volume from the
  local magnetic-field magnitude and the embedded line geometry.
- Parker, `fte-dmumu`, and `fte-mfp` continue to advance the scalar field-line
  coordinate with parallel/normal momentum state and turbulence coupling.

## Scope verification

Run:

```sh
make test-field-line-scope-unit
```

The target implements `SCOPE01` (forbidden transferred symbols and objects),
`SCOPE02` (no Cartesian mover advance or AMR-cell attachment), and `SCOPE03`
(3-D field-line geometry and observer sampling retained). It also reruns the
Step 4 registry test to prove the three canonical mover choices remain intact.

The native completion gate additionally requires independent builds of srcSEP
and the receiving 3-D application. This source-only package lacks the enclosing
AMPS dependency tree and contains no receiving application, so those two linked
build results cannot be manufactured here; see the README for the exact gate.

## Files removed

- `drift.cpp`
- `output.cpp`
- `sample3d.h`
- `sample3d_init.cpp`
- `sample3d_sampling.cpp`
- `sample3d_output.cpp`

Generated objects, archives, sanitizer files, coverage files, and test output
must not be committed.
