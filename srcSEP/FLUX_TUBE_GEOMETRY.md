# Flux-tube geometry and source normalization

Step 3 replaces the historical `MagneticTubeRadius` convention with one
dimensionally explicit control-volume interface. A flux tube now exposes area
in square metres and segment volume in cubic metres. No caller squares an
unlabelled value or independently reconstructs a tube volume.

## Geometry contract

`SEP::FieldLine::FluxTubeGeometry` is the PIC-facing interface:

- `AreaAtVertexM2(vertex, line)` returns cross-sectional area in m²;
- `AreaAtSegmentFractionM2(segment, line, fraction)` linearly interpolates the
  conservative endpoint area;
- `SegmentVolumeM3(segment, line)` integrates flux-conserving area over segment
  length with begin/mid/end Simpson quadrature;
- `PartialSegmentVolumeM3(...)` integrates a subinterval; and
- `SweptVolumeM3(...)` returns `A v dt` in m³.

When finite nonzero magnetic fields are stored at both the first and requested
vertices, the default closure is magnetic-flux conservation:

```text
A(s) |B(s)| = A(reference) |B(reference)|.
```

The default reference area is π m², which makes the former implicit one-metre
reference radius explicit without preserving its erroneous radial scaling.
Applications with a measured normalization call `SetReferenceAreaM2()`.

If either magnetic magnitude is absent, zero, NaN, or infinite, geometry throws
unless the application has installed an `ExplicitAreaProfile`. The callback
accepts Cartesian coordinates in metres and returns area in square metres.
There is intentionally no automatic heliocentric-distance fallback: an `r²`
expansion is a physical assumption, not missing-data recovery.

## Shared normalization

The geometry is now used by particle population, solar-wind injection, shock
injection, wave-energy density, turbulence advection, particle-wave growth, and
field-line sampling. The dependency-light
`util/sep_flux_tube_geometry_core.*` also owns the typed operations for flux
conservation, volume integration, swept volume, injected physical count,
macroparticle weight, and spectral number density.

Shock injection uses one efficiency,
`SEP::FieldLine::InjectionParameters::InjectionEfficiency`, for analytic,
SWCME, and SWMF states. The separate hard-coded density/speed/efficiency weight
override in `main_lib.cpp` has been removed. AMPS retains the configured base
particle weight; the injector represents the calculated physical source using
the particle's statistical correction.

Input `emin`, `emax`, and `ConstEnergyInjectionValue` remain in MeV for backward
configuration compatibility. Every route converts with `EnergyFromMeV()`
before calling an SI relativistic momentum function. New dimensional variables
use suffixes such as `_J`, `_m2`, `_m3`, `_m_s`, and `_s`, and shared APIs use
the strong `SEP::Units` wrappers.

## Focused verification

Run:

```sh
make test-geometry-source-unit
```

The dependency-light, strict-warning ASan/UBSan test covers:

- `GEOA01`: `A|B|` invariance;
- `GEOA02`: fourth-order segment-volume convergence;
- `SRC01`: dimensional `A v dt` swept volume;
- `SRC02`: equal physical source for provider-equivalent states despite
  different macroparticle counts; and
- `SRC03`: spectral bins integrate to the source number density, including an
  exact MeV-to-joule boundary check.

The final native completion gate must also run the enclosing AMPS executable
with analytic, SWCME, and SWMF fixtures, plus sanitizer-enabled injection
branches. Those coupled prerequisites are not contained in a source-only
`srcSEP` archive.
