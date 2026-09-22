# Configuration, Geometry, Preflight, and Compiled Species Binding

## Compile-time species selection

The AMPS application deck and the srcSEP3D runtime deck have different roles.
`input/sep3d.input` is processed while AMPS is configured and compiled. Its
`SpeciesList` is the only authority for species count, order, chemical symbol,
mass, and signed charge. The file passed to `./amps --input FILE` is parsed only
after compilation and cannot select, relabel, add, remove, or mutate an AMPS
species. Its `[species]` section therefore contains only the common numerical
`macroparticle_weight`.

An executable generated with `SpeciesList=ELECTRON` is valid and injects an
electron. A mixed `SpeciesList=H_PLUS ELECTRON` build initializes and injects
both entries in generated index order. Runtime never searches for a proton
macro and never overwrites either molecular-table entry.

For an AMPS checkout already configured as `sep3d`, changing the source deck
does not necessarily refresh the root-level working copy. Use
`cp input/sep3d.input sep3d.input`, run
`./ampsConfig.pl -input sep3d.input -no-compile`, and then perform the normal
clean site build. Before execution, inspect the generated `nTotalSpecies` and
`ChemTable` declarations in `build/pic/pic.h`; startup prints the same table
with mass and signed charge after binding.

## V01 perpendicular diffusion and guiding-centre drift

Under `[transport]`, `perpendicular_diffusion` accepts `none`, `constant`, or
`constant-ratio`. The selected closure requires its positive SI coefficient:
`constant_kappa_perpendicular_m2_per_s` or
`kappa_perpendicular_to_parallel_ratio`. `drifts` accepts `none`,
`gradient-b`, `curvature`, or `gradient-curvature`. These choices enter the
physics fingerprint. Drift uses the signed charge from the compiled AMPS table and forces
magnetic-gradient storage before the mesh layout freezes. Current-sheet drift
is unsupported because its geometry is unspecified.

This document describes the production input boundary introduced by the
configuration/preflight improvements and extended by Stage 3 species
ownership. The implementation is split between
`runtime/configuration_io.{h,cpp}`, `runtime/run_configuration.{h,cpp}`,
`core/parker_geometry.{h,cpp}`, and `mesh/mesh_model.{h,cpp}`. None of these
files depends on AMPS or MPI.

## C01: input file and standalone CLI

The standalone executable accepts one versioned INI-style file:

```bash
./amps --input srcSEP3D/examples/sep3d_analytic_parker.in --dry-run
mpiexec -n 4 ./amps --input srcSEP3D/examples/sep3d_analytic_parker.in \
  --initialization-only --initialization-output-dir mesh-preview
./amps --input run.in --output-dir products
./amps --input run.in --restart restart/sep3d.chk --log-level verbose
```

`--output-dir` and `--restart` override only their corresponding file fields.
`--dry-run` parses, normalizes, validates, freezes, fingerprints, samples the
resolution law, and estimates memory; it does not allocate the AMPS mesh.
`--initialization-only` instead follows the production path through AMPS mesh
construction and complete application initialization, writes all three declared
Tecplot products, synchronizes all MPI ranks, and exits before the first
particle step. `--initialization-output-dir` requires that mode and changes
only the parent directory of the three initialization filenames.
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

Every assignment must follow a section header. A flat legacy assignment such
as `scattering = ...` fails with a section-aware diagnostic rather than being
guessed into `[turbulence]` or `[transport]`. This is intentional: those two
groups represent different physical contracts. The current example begins
with `[run]` and is exercised directly by `CFG3D01`.

The complete commented example is
[`examples/sep3d_analytic_parker.in`](examples/sep3d_analytic_parker.in).

### Schema 3 initialization contract

Schema 3 adds the required `[parker_spiral]` section from schema 2, a complete
canonical `[swcme]` section, explicit initialization Tecplot paths, complete
observer acceptance fields, and an exact particle source. It is intended for
standalone SWCME shock-injection runs and therefore requires:

- `run.intent = shock-injection`, `shock.authority = swcme`, and
  `source.enabled = true`;
- `run.injection_cadence_steps = 1`, so
  `source.samples_per_step` means exactly that many computational particles for
  every compiled AMPS species on every active simulation step;
- `background.provider = analytic-parker` and prescribed turbulence selected as
  `kolmogorov`, `kraichnan`, or explicit `power-law`. The recognized
  `python-interpolator` value is reserved and fails before allocation; SWMF/AWSoM
  coupling remains available through the parser-free typed host interface;
- a spherical canonical shock with `shock_only` region behavior,
  source-mode acceleration, and relative-only source normalization. The
  current AMPS crossing operator is spherical, so accepting ellipsoid/SSE here
  would be a physically false approximation;
- a complete `[swcme]` assignment layer: the 41 common SWCME1D fields plus
  shape, axis ratios, half width, CME direction, solar rotation axis/rate, and
  surface theta/phi resolution. Data-driven kinematics additionally owns its
  two knot lists;
- exact agreement between application and canonical descriptions of Parker
  wind/rotation/source radius/reference latitude, +Z rotation axis, magnetic
  normalization/polarity, density, temperature, energy range, and injection
  efficiency. After this check, the runtime Parker state is populated from the
  canonical SWCME result, including adiabatic index, composition, electron and
  alpha temperatures, and thermodynamic closure. Species identity remains
  exclusively in AMPS' compiled table;
  and
- nonempty `output.initialization_mesh_tecplot_file` and
  `output.initialization_parker_line_tecplot_file`, and
  `output.initialization_data_tecplot_file`.

The magnetic comparison respects the two public conventions. SWCME supplies
total `|B|` at one AU and its reference latitude; the analytic provider
supplies radial `Br` at its declared reference radius. The parser removes the
one-AU Parker winding and then applies `Br proportional to r^-2`. Consequently
a half-AU analytic reference requires four times the corresponding one-AU
radial component; it is not compared directly with the total SWCME magnitude.

Density uses the separate, unambiguous key
`background.parker.number_density_at_one_au_m3`. It must match SWCME
`ambient.density_1au` and remains a one-AU normalization when the magnetic
reference radius changes. The legacy `number_density_at_reference_m3` spelling
is accepted only by schemas 1–2 and rejected by schema 3.

Complete prescribed turbulence input additionally declares `model`,
`normalized_cross_helicity`, `reference_radius_m`, independent radial exponents
for both wave-number bounds, a correlation-length radial exponent, and provider
validity cadence. Kolmogorov requires `q=5/3`, Kraichnan requires `q=3/2`, and
the general power-law choice requires an explicit finite `q>1`; a contradictory
named model/index pair fails rather than being normalized silently.

All application fields and every field of every `[observer.ID]` are explicit,
including values inactive under the chosen mode. The `[swcme]` resolver is the
only code that interprets its unit-bearing values. Its normalized manifest and
fingerprint are frozen into the application identity; raw spelling is not used
as a substitute for resolved physics.

Before AMPS allocation, the standalone provider evaluates the canonical model
at `event.valid_from`, builds and validates the complete shock surface, solves
the MHD jump/source state on every patch, and verifies that
the per-species `source.samples_per_step` is at least the active patch count.
This distinguishes
a legitimately delayed event (inactive before `valid_from`) from an invalid or
empty source.

The required particle-weight identity is

```text
species.macroparticle_weight =
  source.physical_particle_rate_per_s
  * source.injection_efficiency
  * swcme.shock.relative_source_weight_per_area
  * run.time_step_s
  / source.samples_per_step
```

Each species' count is apportioned with a deterministic largest-remainder
allocation after reserving one representative per active physical patch.
Per-patch individual weight corrections preserve represented physical number
exactly. The global and every local AMPS timestep/weight are then set for every
compiled species from the frozen configuration.

The full field-by-field tables and `[swcme]` key list are in the top-level
[`README.md`](README.md); the annotated file is executable acceptance input,
not pseudocode.

## C02: typed run contract and fingerprints

The parser produces `RunConfiguration3DOptions`, the same SI-only record a
coupled SWMF host constructs directly. SWMF coupling does not parse a file and
does not maintain a second configuration representation. The successful
`RunConfiguration3D::Create` call returns an immutable normalized object.

The typed contract includes:

- run intent, transport core, time step, maximum steps, random campaign, and
  integer background/injection/sampling/checkpoint cadences;
- explicit domain, boundary, coordinate-frame, mesh, and storage choices;
- complete SWCME-backed analytic Parker/solar-wind parameters and a typed
  reserved Python-interpolator authority;
- turbulence model, amplitude, cross helicity, radial scalings, spectrum,
  cadence, and out-of-range/missing-data policies;
- legacy shock interval/radial/speed/compression for schemas 1–2, or the
  canonical SWCME3D manifest/fingerprint for schema 3;
- source efficiency, physical particle rate, energy interval, spectrum, and
  maximum samples per injection event;
- species AMPS index, name, mass, signed charge, and macroparticle weight;
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

### Complete compiled AMPS binding

Immediately after `PIC::Init_BeforeParser()`, production enumerates
`0..PIC::nTotalSpecies-1` and copies `GetChemSymbol`, `GetMass`, and
`GetElectricCharge` into AMPS-independent `CompiledSpeciesRecord` values.
`ValidateCompiledSpeciesBinding` requires the copied count to equal the
generated count, indices to be contiguous, symbols to be unique and nonempty,
masses to be finite and positive, and charges to be finite and nonzero. A
neutral entry fails because the selected focused SEP scattering operator
requires charge; treating it as charged would not be physically correct.

The adapter does not call the AMPS molecular-data setters. It then installs the
explicit runtime timestep and common base statistical weight in every compiled
global slot and every allocated local block. Observer comma-separated numeric
indices are checked against the generated count during this binding, when that
count is authoritative and available. `species = all` selects the entire
compiled table without embedding a build-specific count in the runtime file;
use numeric indices only for an intentional observer subset.

Source controls have per-compiled-species semantics. Every species receives
`source.samples_per_step` computational particles and the declared per-species
physical rate. The kinetic-energy bounds are converted to momentum separately
with that species' AMPS mass; SWCME shock geometry and compression remain
common. This is complete multi-species initialization/injection without a
second runtime species-definition authority.

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
`species = all` or comma-separated subset IDs, and `normalization`. Observer cadence in seconds
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
| `CFG3D07` | complete generated table, mixed signed charges, count/index/symbol/mass/charge failures, observer bounds, and fingerprint identity |
| `CFG3D08` | complete schema-3 SWCME input plus missing-field, weight, and per-step-cadence rejection |
| `MSH3D11` | deterministic unit-labeled initialization Parker Tecplot output |
| `R3D08` | canonical source-surface preflight and exact per-species, per-step particle allocation |

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

The parser accepts schema version 1 for restart/campaign compatibility. Its
finite line is normalized from the existing source direction, domain radii,
and a deterministic count before fingerprinting.  Schema version 2 never
falls back to those values: omitting any of the eight keys is an error. Schema
version 3 includes the same explicit line and adds the complete initialization
contract above.

Run only these gates with:

```bash
test/run_tests.py --suite improvements-c --rebuild \
  --output-dir test_output/improvements-c
```
