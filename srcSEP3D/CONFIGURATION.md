# Configuration, Geometry, and Preflight (C01-C05)

## V01 perpendicular diffusion and guiding-centre drift

Under `[transport]`, `perpendicular_diffusion` accepts `none`, `constant`, or
`constant-ratio`. The selected closure requires its positive SI coefficient:
`constant_kappa_perpendicular_m2_per_s` or
`kappa_perpendicular_to_parallel_ratio`. `drifts` accepts `none`,
`gradient-b`, `curvature`, or `gradient-curvature`. These choices enter the
physics fingerprint. Drift uses signed `[species] charge_c` and forces
magnetic-gradient storage before the mesh layout freezes. Current-sheet drift
is unsupported because its geometry is unspecified.

This document describes the production input boundary introduced by
improvements C01-C05. The implementation is split between
`runtime/configuration_io.{h,cpp}`, `runtime/run_configuration.{h,cpp}`,
`core/parker_geometry.{h,cpp}`, and `mesh/mesh_model.{h,cpp}`. None of these
files depends on AMPS or MPI.

## C01: input file and standalone CLI

The standalone executable accepts one versioned INI-style file:

```bash
./amps --input srcSEP3D/examples/sep3d_analytic_parker.in --dry-run
./amps --input run.in --output-dir products
./amps --input run.in --restart restart/sep3d.chk --log-level verbose
```

`--output-dir` and `--restart` override only their corresponding file fields.
`--dry-run` parses, normalizes, validates, freezes, fingerprints, samples the
resolution law, and estimates memory; it does not allocate the AMPS mesh.
Unknown options, conflicting test selectors, or a missing `--input` are usage
errors and return code 2 before AMPS initialization.

The public test interface remains `test/run_tests.py`; it provides `--list`,
`--test`, `--group`, `--routine`, `--all`, and `--suite`. The production
executable recognizes test-selector spellings only to reject an accidental
simulation launch with a diagnostic directing the user to the test runner.

Every dimensional file key includes its SI unit, for example
`inner_radius_m`, `time_step_s`, `radial_field_at_reference_t`, and
`minimum_energy_j`. A key without the declared suffix is unknown. The parser
does not guess whether a number is AU, solar radii, km/s, nT, or MeV.

The required groups are `run`, `domain`, `mesh`, `mesh.solar`, `mesh.tube`,
`memory`, `background`, `background.parker`, `turbulence`, `transport`,
`shock`, `source`, `species`, `storage`, `output`, `restart`, and at least one
`observer.ID`. A disabled feature still has a section, making the decision
visible in code review. Duplicate keys, unknown keys, invalid enumerations,
missing groups, non-finite values, and incompatible choices fail in the parser
or immutable factory.

The complete commented example is
[`examples/sep3d_analytic_parker.in`](examples/sep3d_analytic_parker.in).

## C02: typed run contract and fingerprints

The parser produces `RunConfiguration3DOptions`, the same SI-only record a
coupled SWMF host constructs directly. SWMF coupling does not parse a file and
does not maintain a second configuration representation. The successful
`RunConfiguration3D::Create` call returns an immutable normalized object.

The typed contract includes:

- run intent, transport core, time step, maximum steps, random campaign, and
  integer background/injection/sampling/checkpoint cadences;
- explicit domain, boundary, coordinate-frame, mesh, and storage choices;
- complete analytic Parker parameters;
- turbulence spectrum and out-of-range/missing-data policies;
- shock interval, radial extent, speed, and compression;
- source efficiency, physical particle rate, energy interval, spectrum, and
  maximum samples per injection event;
- species name, mass, signed charge, and macroparticle weight;
- named observer geometry/trajectory, collection or shell radius, cadence,
  species/pitch acceptance, energy range/bin counts, normalization, and products;
- output/restart paths and cadence.

The run intent prevents ambiguous source behavior. `transport-only` rejects an
enabled source. `shock-injection` requires both SWCME shock authority and an
enabled source. SWMF turbulence requires SWMF background authority, and an
`imported-coverage` outer boundary likewise requires SWMF background.

The physics fingerprint includes all trajectory, source, observer, mesh,
layout, and resource-model inputs. Output directory, output prefix, restart
paths, and output cadence are included in the resolved manifest but excluded
from the physics fingerprint. Consequently a cosmetic output relocation can
resume the same physics, while a species weight, observer, mesh, shock, or
transport change cannot masquerade as the same run. Restart loading compares
the frozen physics/layout identities transactionally.

## C03: domain and boundary contract

`outer_radius_mode=preset` resolves a named radius during normalization:

| Preset | Resolved outer radius |
|---|---:|
| `solar` | 0.30 AU |
| `one-au` (`earth` input alias) | 1.00 AU |
| `mars` | 1.666 AU |

`outer_radius_mode=explicit` uses `outer_radius_m`; no numeric sentinel means
“automatic.” The current geometry is heliocentric, so the coordinate origin
must be `(0,0,0)` and the named frame must agree with the Parker or supplied
SWMF frame. Observers, shock initial radius, shock maximum radius, Parker
reference radius, and tube reference radius are validated against the
normalized domain before mesh initialization.

The inner sphere is absorbing. An outward segment that crosses the outer
sphere returns `DomainExit`; the message distinguishes ordinary escape from
leaving imported SWMF coverage. Direction matters: a segment entering the
outer sphere is not an escape, and a segment moving outward across the inner
sphere is not solar absorption.

## C04: one Parker geometry

`ParkerSpiralGeometry` is the single geometric definition used by both the
mesh tube and analytic background. For source position \(\mathbf{x}_0\),
rotation axis \(\hat{a}\), source radius \(r_0\), rotation rate \(\Omega\), and
wind speed \(V\), the curve rotates the source direction by

\[
\Delta\phi(r)=-\Omega\max(0,r-r_0)/V.
\]

Its local outward tangent is

\[
\hat{t}=\operatorname{normalize}
\left[\hat{r}-\frac{\Omega\max(0,r-r_0)}{V}
(\hat{a}\times\hat{r})\right].
\]

Magnetic polarity is deliberately absent from this geometry. The analytic
provider applies polarity only after computing the unsigned tangent, reversing
`B` and pitch orientation without moving the mesh tube. Tube distance is the
spherical transverse arc
`r*atan2(|r_hat cross t_hat|, r_hat dot t_hat)` at the point's radius.

## C05: composite refinement and whole-run memory

The global target is `global_cell_size_m`. Near the Sun, normalized coordinate

\[
s=(r-r_{in})/(r_{transition}-r_{in})
\]

is clamped to `[0,1]` and mapped by `linear`, `power-law`, or `smoothstep`.
The result interpolates from `surface_cell_size_m` to the global target.

The tube has a declared physical radius at a declared heliocentric reference
distance. `physical-constant` keeps that radius fixed;
`constant-angular-width` scales it as `radius_at_reference_m*r/reference_r`.
The transverse profile interpolates from `center_cell_size_m` on the Parker
centerline to the global target at the tube boundary. Where near-Sun and tube
regions overlap, the finer request wins. The final value is clamped between
the declared minimum and global sizes.

The dry-run preflight reports sampled minimum/maximum requests and locations,
estimated blocks by AMR level, tube radius at the reference distance, and a
conservative whole-run memory model. Its categories are:

- base cell objects;
- frozen associated background/turbulence bytes;
- sampling bytes;
- mesh node objects;
- block structures and application block overhead;
- resident particles;
- communication buffers and halo fraction;
- explicit safety margin.

All byte/object coefficients and both fractions are input parameters because
the exact native sizes depend on the AMPS build. The estimate is a planning
gate; native integration tests compare it with actual AMPS counts and peak
memory. Configuration also proves that `maximum_level` can realize the finest
requested cell size. An impossible level cap or an estimate above the memory
budget fails before full allocation.

## R04–R06 runtime controls

The production clock and recurring services are configured with integer
cadences:

| Key | Meaning |
|---|---|
| `run.background_cadence_steps` | request the next background/turbulence generation |
| `run.injection_cadence_steps` | evaluate shock patches and inject a source event |
| `output.cadence_steps` | gather and publish observer products |
| `output.checkpoint_cadence_steps` | write a complete R07 restart; zero disables periodic checkpoints |
| `transport.maximum_substeps` | hard accepted-substep limit for one requested AMPS interval |
| `source.physical_particle_rate_per_s` | physical source rate normalized over the injection interval |

Observer sections additionally accept `kind`, Cartesian velocity, collection
radius, shell radius, minimum/maximum energy, minimum/maximum pitch cosine,
comma-separated species IDs, and `normalization`. Observer cadence in seconds
must be an exact integer multiple of `run.time_step_s`; normalization therefore
cannot acquire a drifting fractional-tick window.

All trajectory, cadence, source, and observer acceptance fields participate in
the physics fingerprint. Checkpoint and output paths remain relocation-only
manifest fields. `R3D04–06` cover clock/event restoration, source cadence
identity/conservation, and resolved observer geometry plus commit-only reset.

## Acceptance evidence

| ID | Contract |
|---|---|
| `CFG3D01` | schema, CLI, early errors, typed/file fingerprint parity, and dry-run |
| `CFG3D02` | typed groups, field classification, compatibility, and SWMF/analytic factory parity |
| `CFG3D03` | presets, explicit override, containment, and directional boundary status |
| `CFG3D04` | shared tangent and polarity-independent tube geometry |
| `CFG3D05` | monotone composite profiles, tube scaling, AMR levels, memory categories, and level rejection |
| `CFG3D06` | complete finite Parker-line input and fail-closed source consistency |

## Version 2 finite Parker-line section

Schema version 2 requires these eight keys:

```ini
[parker_spiral]
origin_x_m = 0
origin_y_m = 0
origin_z_m = 0
initial_x_m = 1.3914e10
initial_y_m = 0
initial_z_m = 0
length_m = 2.0e11
point_count = 4001
```

`point_count` includes both endpoints and must be at least two.  `length_m` is
arc length, not final heliocentric radius.  The origin must equal the domain
origin, and the initial point relative to that origin must lie on
`domain.inner_radius_m` in the direction declared by the tube source longitude
and colatitude.  Solar wind speed and rotation remain owned by
`[background.parker]`; the line section does not duplicate them.

The parser accepts schema version 1 for restart/campaign compatibility.  Its
finite line is normalized from the existing source direction, domain radii,
and a deterministic count before fingerprinting.  Schema version 2 never
falls back to those values: omitting any of the eight keys is an error.

Run only these gates with:

```bash
test/run_tests.py --suite improvements-c --rebuild \
  --output-dir test_output/improvements-c
```
