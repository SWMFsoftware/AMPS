# AMPS Runtime Input File

This document is the normative specification of the **generic AMPS runtime
input format**: the one file format used to set up a run of any AMPS
application (planets, moons, comets, spacecraft re-entry, SEP and GCR
transport).  The reader that implements the grammar is the header-only
`src/general/input/amps_input.h`; its regression test is
`src/general/input/test_amps_input.cpp` (`make -C src/general/input test`).

## 1. Two input layers

| Layer | File | Decides | Changes require |
|---|---|---|---|
| Compile time | `input/<application>.input`, read by `ampsConfig.pl` | which species, modules and physics libraries are compiled in | a rebuild |
| Run time | `amps.in` (or `-input FILE`) — this document | every value: mesh, model selections and parameters, sources, observers, output | an edit |

The runtime file can choose only among what was compiled in; it cannot add or
relabel species.

## 2. File selection

### 2.1 `srcMoon` application block

The lunar executable selects the same generic runtime file with `-input`:

```bash
./amps -input srcMoon/examples/lola_surface.in
```

Omitting `-input` deliberately retains the historical analytic-sphere
regression. The `application: moon` block has no implicit values; all keys below are
required even when `surface_geometry=sphere`, so a run receipt remains
complete and switching geometry cannot expose stale undeclared paths:

```text
#begin application: moon
spice_path = /home/vtenishe/SPICE
surface_geometry = lola
surface_mesh_resolution_m = 100000
lola_product_id = LDEM_4
lola_image_file = /data/vtenishe/moon_validation_data/lola/raw/LDEM_4.IMG
lola_label_file = /data/vtenishe/moon_validation_data/lola/raw/LDEM_4.LBL
surface_cea_file = test_output/srcMoon/surface/lola_surface.cea
surface_tecplot_file = test_output/srcMoon/surface/lola_surface.dat
mesh_profile = full
mesh_near_surface_cell_size_m = 10000
mesh_near_surface_extent_m = 200000
mesh_far_field_cell_size_m = 40000
mesh_far_field_extent_m = 10000000
mesh_outer_cell_size_m = 17371000
mesh_output_tecplot = false
mesh_angular_resolution_m_by_sza = 0:10000,30:10000,60:5000,90:5000,120:10000,180:10000
mesh_terminator_half_width_deg = 15
mesh_polar_cap_start_latitude_deg = 75
mesh_polar_resolution_m = 5000
simulation_start_utc = 2009-01-24T00:00:00
simulation_environment_mode = frozen
random_seed = 731993
profile_altitude_near_maximum_m = 200000
profile_altitude_near_step_m = 5000
profile_altitude_middle_maximum_m = 1000000
profile_altitude_middle_step_m = 25000
profile_altitude_far_maximum_m = 10000000
profile_altitude_far_points_per_decade = 8
profile_solar_zenith_angles_deg = 0,30,60,90,120,180
profile_include_body_fixed_poles = true
profile_include_shell_average = true
profile_shell_angular_weighting = solid-angle
limb_observer = earth
limb_plane = sun-moon-earth
limb_sides = both
limb_branch_reference = projected-sun
limb_tangent_maximum_altitude_m = 1000000
tail_center_direction = anti-lunar
tail_map_coordinate_system = j2000-ra-dec
tail_map_display_orientation = west-right
tail_map_half_width_deg = 7.5
tail_beam_diameter_deg = 1.0
tail_exposure_duration_s = 120
tail_beam_response = unresolved
tail_temporal_operator = exposure-average
#end application
```

Relative paths are resolved against the file containing each assignment. The
SPICE root must have `cspice/include/SpiceUsr.h`, `cspice/lib/cspice.a`, and
`Kernels/`. `surface_mesh_resolution_m` is the maximum requested great-circle
edge on the pre-topography reference sphere, not a latitude/longitude step or
an AMR cell size. The CEA and Tecplot output paths must differ. Runtime path
selection does not change the build-time SPICE mode.

`mesh_profile` is a required provenance label (`quick`, `full`, or `custom`);
it does not hide a second table of compiled defaults. `full` identifies the
frozen M1-resolution receipt, `quick` identifies the committed approximately
20-minute wiring profile, and `custom` records any investigator-selected
resolution set. All radial
controls are metres. The near-surface target applies from the lunar reference
radius through `mesh_near_surface_extent_m`; the selector then transitions in
dimensionless log-radius to `mesh_far_field_cell_size_m` at
`mesh_far_field_extent_m` altitude. Outside that declared science region it
coarsens monotonically to `mesh_outer_cell_size_m` at the compiled outer
boundary. The committed M1 values therefore request 40 km at 10,000 km but do
not impose that resolution throughout the 400-lunar-radius tail box.

The committed QUICK receipt uses three times every FULL spatial target while
retaining the same physical extents and non-mesh controls. This factor is
empirical, not hidden profile logic: the earlier factor-two receipt produced
1,074,788 blocks and its eight-rank I12 run was stopped after 28.7 minutes at
223.7 GiB RSS. Factor three crosses the next coarser dyadic AMR level. Its cold
eight-rank I12 wiring run completed in 1,409.9 seconds (23.50 minutes), with
191.9 GiB peak RSS, on `passion` on 2026-10-12. This timing is host/cache-
dependent operational guidance, not a validation criterion. QUICK cannot
substitute for FULL I11-I13 convergence. To choose
another resolution, edit the explicit target fields and set
`mesh_profile = custom`; no resolution is selected implicitly by the profile
name.

The aggregate srcMoon developer regression uses the separate committed
`srcMoon/examples/m1_regression.in` receipt.  It preserves the FULL physical
extents, epoch, sources, science grid, and observation geometry, but sets all
spatial targets to twelve times their FULL values and declares
`mesh_profile=custom`.  It exists to bound routine production-wiring runs;
even when its evaluated wiring criteria succeed, its M1 result is
`SKIPPED / NOT QUALIFIED`.  Formal I11--I13 qualification still requires the
FULL receipt and unchanged thresholds.

`mesh_output_tecplot` is a required Boolean. `false` suppresses the optional
full-volume `mesh.dat` and `tetra_mesh*.dat` visualization exports without
changing the AMR tree, particle calculation, samplers, or the internal
`mesh.msh` file used by the existing MPI initialization path. Set it to
`true` only when the output filesystem has been deliberately provisioned: a
fine three-dimensional mesh is written as text and can require hundreds of
gigabytes. This switch does not control `surface_tecplot_file`, which is the
separate triangulated lunar surface product.

`mesh_angular_resolution_m_by_sza` is a comma-separated `SZA:resolution`
table in degrees and metres. It must start at 0 degrees, contain an exact
90-degree terminator knot, end at 180 degrees, and have strictly increasing
angles and positive resolutions. Piecewise-linear interpolation defines the
requested boundary resolution. Within the declared terminator half-width the
90-degree value is an upper bound; poleward of the absolute body-fixed
latitude threshold, `mesh_polar_resolution_m` is an upper bound. Solar zenith
angle is evaluated in the Moon-centred LSO frame, while latitude is evaluated
in `MOON_ME_DE421` at `simulation_start_utc`. AMPS realizes only its discrete
power-of-two AMR levels, so the requested and realized sizes are both retained
by I12.

`simulation_start_utc` is the single run-wide UTC passed to CSPICE after the
configured kernels are furnished. It controls dynamics, illumination, and
diagnostic geometry together; it is not a separate observation-only clock.
`simulation_environment_mode=frozen` holds those ephemeris-driven quantities
at that UTC while particle transport continues, providing the stationary
problem required for I11 spatial convergence. `evolving` advances the SPICE
epoch normally and is intended for transient and observational campaigns.
`random_seed` is a positive integer installed before any production source,
surface interaction, or mover draws from AMPS's random stream. The frozen M1
campaign uses `731993`; the MPI/OpenMP layout is recorded separately in each
linked-run manifest.
The profile grid uses metres above the AMPS lunar reference sphere
(`_RADIUS_(_TARGET_)`, currently 1737.10 km), not above the 1737.4 km LOLA
datum and not above local terrain. Its near and middle endpoints must be
integral numbers of their respective linear steps.
Above the middle endpoint, points are logarithmically spaced at
`profile_altitude_far_points_per_decade`; the far endpoint is always included.
For the shown values the grid contains 81 points: 0–200 km by 5 km, 200–1000
km by 25 km, and 1000–10000 km at eight points per decade.

`profile_solar_zenith_angles_deg` is a strictly increasing comma-separated
list in `[0,180]` for profiles in the LSO equatorial plane. The pole flag adds
the `+Z` and `-Z` directions of `MOON_ME_DE421`; no longitude is assigned at a
pole. The shell-average flag requests a distinct angular average rather than
using any one radial profile as a surrogate. `solid-angle` is required so an
implementation cannot substitute an AMR-cell-count or unweighted-facet mean.
The approved limb geometry is an Earth observer and the Sun-Moon-Earth plane.
`limb_sides` can request `both`, `toward-projected-sun`, or
`away-from-projected-sun`; branch names refer to the projected Moon-to-Sun
direction in the observer sky plane. Its tangent maximum must be one of the
common altitude-grid points.

The tail map is centred on the Earth-sky anti-lunar direction—opposite the
instantaneous Earth-to-Moon vector—not the anti-solar direction. Tail-map
storage may be selected as `j2000-ra-dec` (the recommended native WHAM table
coordinates) or `eclipj2000-lon-lat`; both store the horizontal coordinate
positive east and the vertical coordinate positive north. `west-right` and
`east-right` affect display layout only and never change stored coordinates.
The beam diameter is a circular field-of-view diameter, not map-pixel spacing.
The exposure duration is seconds and `exposure-average` requires temporal
integration over that interval. Because D14 does not supply an authoritative
radial beam response, `tail_beam_response` must remain `unresolved`; the code
does not silently substitute a Gaussian or uniform disk.

The parser records and validates these controls. The existing subsolar-limb sampler uses
the exact altitude nodes through the declared tangent maximum, but it does not
yet implement the complete two-limb/pole/shell operators or circular-aperture
and exposure-time convolution; therefore configuration acceptance and partial
grid wiring are not U26/I33 verification.

#### M2 lunar environment fields

`m2_environment=disabled` is the exact M0/M1 path; every M2-only assignment is
then rejected so a misspelled or ignored physics receipt cannot look active.
`enabled` requires the complete group shown in
`srcMoon/examples/lola_surface.in`:

- the executable must be built with SPICE. Every M2 temperature mode needs an
  epoch-dependent LSO/MOON_ME_DE421 transform, and the predictive mode also
  needs Sun/Earth ephemerides. A no-SPICE build rejects M2 rather than using an
  identity transform or frozen illumination;

- `surface_temperature_mode` is `simple-cosine`, `diviner-map`, or
  `thermal-model`. `diviner-map` additionally requires an existing processed
  D03 CSV; `none` is valid only for the other modes.
- `diviner_eclipse_history_file` and
  `diviner_eclipse_selection_contract_file` are the paired I17 observation
  receipt. Both may be literal `none` in a non-I17 run. Otherwise both paths
  must exist. The history is a reproducibly processed PDS4 RDR table; the JSON
  contract freezes the June-training/December-holdout membership, actual
  footprint/time binning, native flags, Paige channel-sensitivity rule, and
  uncertainty propagation before hold-out residuals are evaluated. These
  files are observation inputs only: they never initialize, force, or tune the
  production thermal state.
- `illumination_ray_origin_offset_m` offsets a terrain ray from its source
  facet. `surface_collision_maximum_segment_length_m`,
  `surface_collision_time_tolerance_s`,
  `surface_collision_distance_tolerance_m`, and
  `surface_collision_maximum_iterations` control moving-boundary bracketing
  and bisection. All are positive and must be converged by I14.
- `thermal_*` fields explicitly provide solar flux, emissivity, albedo law,
  density/conductivity depth laws, the five heat-capacity polynomial
  coefficients, geothermal flux, depth/layers, update size, initial
  temperature, and nonlinear tolerances. Dimensional suffixes are normative.
- Eclipse forcing is also explicit. The currently implemented experimental
  contract requires
  `thermal_eclipse_direct_irradiance_model=apparent-disc-overlap`,
  `thermal_solar_disc_radiance_model=uniform`, and
  `thermal_earth_atmosphere_transmission_model=disabled`. At each surface
  point and SPICE epoch, the code computes Sun/Earth centre separation and
  angular radii, evaluates the analytical overlap of their circular apparent
  discs, and multiplies direct solar flux by the unocculted solar-disc area.
  `thermal_solar_disc_radius_m` and `thermal_opaque_earth_radius_m` state the
  physical radii used in that geometry. Terrain horizon blocking remains an
  independent all-or-none ray test to the solar centre. The unocculted area
  fraction scales the direct flux evaluated with the solar-centre incidence
  cosine; neither incidence nor terrain blocking is integrated across the
  finite solar disc. This first model assumes uniform projected solar
  radiance, a spherical opaque solid Earth, and zero refracted,
  scattered, or transmitted terrestrial-atmosphere irradiance. It therefore
  produces `EXPERIMENTAL / NOT VALIDATED` eclipse comparisons; the parser
  rejects names implying unimplemented limb darkening or atmosphere physics.
- `thermal_spinup_cycle_start_utc` and `thermal_spinup_cycle_end_utc` define
  the observation-free forcing interval replayed before an I17 comparison.
  No Diviner temperature initializes the model. Replay continues without an
  arbitrary cycle-count limit until the largest change at any surface or
  subsurface thermal node is no greater than
  `thermal_spinup_temperature_tolerance_k` and the largest whole-column
  sensible-enthalpy change is no greater than
  `thermal_spinup_energy_tolerance_j_m2`. Both tolerances are positive input
  values. CSPICE controls UTC/TDB conversion; a non-positive interval is a
  runtime error. Allocation or job-resource exhaustion is an `ERROR`, not
  evidence of convergence.
- Na and Ar each declare a `*_sticking_model`, a strictly increasing
  `temperature_k:probability` table, an `*_accommodation_speed_model`, a
  bounded accommodation coefficient, and a non-empty literature-reference
  receipt. `constant` sticking requires identical probability values at all
  table nodes; `temperature-table` uses bounded linear interpolation and
  clamps only outside the declared endpoint range. Accommodation may be
  `linear-speed`,
  `vf=a*vt+(1-a)*vi`, or `kinetic-energy`,
  `vf=sqrt(a*vt^2+(1-a)*vi^2)`. Here `vi` and `vf` are impact/re-emission
  speeds relative to the moving body-fixed surface and `vt` is the production
  Maxwellian draw at the local surface temperature. The M2 Na example selects
  the kinetic-energy definition and `a=0.2` from Tenishev et al. (2013),
  Eq. 14. Ar additionally uses `argon_desorption_energy_j_mol`; Na uses
  `sodium_activation_energy_ev` and `sodium_attempt_frequency_hz`.
  `ion_surface_interaction_mode=unconfigured` is presently the only accepted
  value and makes Na+/Ar+ impacts fail rather than inherit neutral physics.
- The five `argon_source_*` fields form a separate radiogenic-source receipt.
  `argon_source_enabled` is Boolean; the maintained Na-only example sets it
  to `false` but still records the candidate parameters. The only accepted
  geography and time law are `uniform-global` and `continuous`.
  `argon_source_rate_atoms_s` is the total lunar rate in atoms per second,
  not a surface flux; it must be positive and finite. The selected example
  uses the `3.8e20 atoms s^-1` lower end of the Killen (2002) diffusion range.
  `argon_source_map_file` must be literal `none`, because local D10 covers
  only the South Pole-Aitken region and is QA-incomplete; the parser rejects
  a potassium-weighted map rather than extrapolating it. The non-empty
  `argon_source_reference` is a provenance receipt. Enabling injection also
  requires an Ar-enabled build and the dedicated `RadiogenicArgon` Config.pl
  source block; a missing compile-time Ar species is a fatal configuration
  error, never an ignored source.
- `cold_trap_retention_timescale_s` is a diagnostic comparison timescale; it
  does not create a permanent sink. `surface_inventory_initial_state` is
  `zero` or `checkpoint`. The checkpoint setting is a base name: restart
  requires both `<base>.Na` and `<base>.Ar`. Literal `none` disables periodic
  writes and is normalized to an empty internal path.

The sodium surface receipt cites Tenishev et al. (2013),
doi:10.1016/j.icarus.2013.08.021. Its Fig. 3 supplies the temperature-dependent
sticking curve (adapted there from Yakshinskiy and Madey, 2005), while Eq. 14
supplies the energy-mixing definition and the paper evaluates `a=0.2`.
Selecting another model or coefficient is permitted only through the explicit
per-species fields above; AMPS does not silently substitute the legacy
linear-speed interpretation.

### 2.2 General selection rules

```bash
./amps -input path/to/run.in    # the AMPS-wide spelling
./amps --input path/to/run.in   # exact synonym, kept for existing scripts
./amps                          # reads ./amps.in
```

Without `-input`/`--input` both applications read `./amps.in` from the working
directory.  srcSEP3D requires it; srcSEP reads it when it exists and otherwise
keeps its historical no-argument run (the legacy built-in mesh).  Help, list
and native/component-test modes never pick up a working-directory file.

A coupled SWMF run has no command line; AMPS then reads `./amps.in` in the run
directory.  For srcSEP3D this happens in `amps_init_mesh` when no host has
installed a configuration; a coupled srcSEP run does not read a file yet.

## 3. Grammar

The file is a tree of named blocks.

```text
! a comment runs from '!' to the end of the line
#begin mesh                      ! opens a block of kind "mesh"
  maximum_level = 7
  #begin refinement: solar       ! a named instance of kind "refinement"
    model = radial-smoothstep
    surface_cell_size_m = \      ! a trailing '\' continues the line
      1.495978707e9
  #end refinement                ! closes the innermost block; kind must match
#end mesh
#include "parts/physics.in"      ! textual inclusion
```

Rules (every violation is fatal and reported as `file:line: reason`, followed
by the offending line):

* **Comments.** `!` starts a comment that runs to the end of the physical
  line.  `#` is reserved for directives; there is no `#` comment.
* **Continuation.** A `\` that is the last non-blank character before a
  comment joins the next physical line (one space is inserted).  A
  continuation at the end of a file is an error; diagnostics name the first
  physical line of the statement.
* **Blocks.** `#begin KIND` or `#begin KIND: NAME` opens a block and
  `#end KIND` closes the innermost open block, whose kind must match.
  `KIND` is a lower-case letter followed by lower-case letters, digits and
  `_` (directives and kinds are case-insensitive and stored in lower case).
  `NAME` contains letters, digits, `_` and `-` (case-sensitive; no dots).
  Blocks nest to any depth.  Two sibling blocks may not have the same kind and
  name.
* **Assignments.** `KEY = VALUE`, one per logical line, always inside a block.
  `KEY` starts with a letter and contains letters, digits, `_`, `.` and `-`;
  keys are case-sensitive.  `VALUE` is the trimmed rest of the line and must
  not be empty; it may contain blanks (`400 km/s`, `0, 0, 1`).  A key may
  appear only once in a block.
* **Includes.** `#include PATH`, `#include "PATH"` or `#include <PATH>`
  inserts a file at that point.  A relative path is resolved from the
  directory of the file that contains the directive, not from the working
  directory.  Includes may be nested (at most 64 deep); cycles are errors.
  Every block opened in a file must be closed in the same file, so an included
  file is always block-balanced; it may still be included inside an open block
  to supply that block's assignments.
* **Unknown directives** (`#section`, `#subsection`, ...) and the retired INI
  syntax (`[section]`) are errors.

## 4. Ownership: block kinds and owners

Every block is handed to exactly **one owner**, and the owner validates its
keys, values and units.  The owner is the code directory that parses the
block, so adding a model never touches a central parser.  By AMPS policy a
block whose owner is not part of the executable is a **fatal error**, never
ignored (`AMPS::Input::RequireOnlyBlocks`).

| Block | Owner | Contents |
|---|---|---|
| `run` | AMPS core / application | time stepping, stop condition, random seed; nested `restart` |
| `mesh` | AMPS core mesh / application | AMR levels and cell sizes; nested `domain`, `field_line`, `refinement: NAME`, `active_region`, `memory` |
| `species` | AMPS core | runtime settings of compiled species (e.g. statistical weight) |
| `body: NAME` | core or the planet/moon model | internal bodies (radius, surface model, boundary condition) |
| `physics: NAME` | the shared library in `src/models/...` | background plasma/field, collisions, chemistry, photolysis, turbulence, parallel diffusion, CME/shock models |
| `source: NAME` | core injection plus the supplying model | surface, boundary, shock or plume injection |
| `observer: NAME`, `output` | core sampling | observers/detectors, output cadence and files |
| `application: NAME` | `srcSEP/`, `srcSEP3D/`, ... | only what one application knows (mover numerics, population control) |

**Model blocks.** A `physics:` (or `source:`) block selects its model and holds
only that model's parameters.  A model library keeps its own key vocabulary:
SWCME keeps units in the value (`ambient.wind_speed = 400 km/s`), the shock
front and the CME event models keep dotted keys (`ambient.reference_radius_m`).
Inside a model block a nested block is the same key space as the dotted
prefix: `#begin ambient` + `x = 1` equals `ambient.x = 1`
(`AMPS::Input::Flatten`); writing the same key both ways is a duplicate.

**Units.** Application keys carry their SI unit in the key name
(`inner_radius_m`, `time_step_s`, `minimum_energy_j`); a bare number is never
interpreted in AU, R_sun, km/s or MeV.  Model blocks follow their library's
documented rule (SWCME: unit in the value).

**Qualified names.** A key's qualified name is the dot-joined path of the
blocks that contain it, where `#begin kind: name` contributes `kind.name`:
`mesh.domain.inner_radius_m`, `mesh.refinement.solar.profile`,
`physics.background.parker.temperature_k`, `observer.earth.cadence_s`.
Diagnostics, documentation and tests use these names.

## 5. Shared SEP layout (srcSEP and srcSEP3D, content schema 6)

srcSEP (field-line transport) and srcSEP3D (three-dimensional transport) read
**one** vocabulary.  Its complete table of blocks and keys, each marked with
the application(s) that implement it, is
`src/models/sep_common/sep_input_layout.h`; both parsers walk their documents
against that table, so the two applications cannot spell the same setting
differently.  Owners: `srcSEP3D/runtime/configuration_io.cpp` and
`srcSEP/util/sep_initialization.cpp`.  The full srcSEP3D key reference is
`srcSEP3D/CONFIGURATION.md`; the srcSEP keys are documented in
`srcSEP/README.md`.

### 5.1 Decisions

| | Rule |
|---|---|
| D1 | One application block, `#begin application: sep`, and one content schema, `schema_version = 6`.  The earlier schemas (srcSEP 1-4, srcSEP3D 1-5) are retired: a file that names one stops with a "retired" diagnosis and the converter command. |
| D2 | The turbulence description is `#begin turbulence` inside `application: sep`.  It is the **only** block whose content differs: srcSEP3D describes a prescribed spectrum, srcSEP its self-consistent Alfvén-wave model.  srcSEP's turbulence command-line options remain and override the file. |
| D3 | The vocabulary is identical; an application rejects any setting it cannot honour with a typed `... is a srcSEP3D setting that is not supported by srcSEP` (or the reverse) error, never ignoring it.  This holds for whole blocks (`mesh.active_region` in srcSEP), single keys (`observer.ID.heliocentric_radius_m` in srcSEP3D) and values (`physics.background.provider = analytic-parker` in srcSEP). |
| D4 | Every physical quantity is written once.  With `physics: swcme` the SWCME ambient owns the solar wind, so the Parker keys `reference_radius_m`, `radial_field_at_reference_t`, `solar_rotation_rate_rad_per_s`, `solar_wind_speed_m_per_s`, `magnetic_polarity`, `number_density_at_one_au_m3` and `temperature_k` are illegal in `physics: background` / `parker` (srcSEP3D derives them from the resolved SWCME model; without SWCME they are required).  `source: shock` owns the particle source, so the SWCME `source.*` keys are illegal in `physics: swcme`; both applications synthesize them from `source: shock`. |
| D5 | One observer vocabulary built on `kind`: srcSEP implements `field-line-radius` (`heliocentric_radius_m` on the field line), srcSEP3D the spatial kinds (`fixed-cartesian`, `fixed-heliographic`, `moving-cartesian`, `spherical-shell`, `field-connected`). |

### 5.2 Blocks

| Block path | srcSEP | srcSEP3D | Contents |
|---|---|---|---|
| `run` | ✓ | ✓ | `time_step_s`, `maximum_time_steps` (both); `time_step_mode`, `maximum_particle_speed_m_s`, `time_step_margin_factor`, `stop_shock_radius_m`, `campaign_seed` (srcSEP3D) |
| `run` / `restart` | – | ✓ | `input_path`, `output_path` |
| `application: sep` | ✓ | ✓ | `schema_version = 6`, `transport` (both); `intent`, `background_cadence_steps`, `injection_cadence_steps` (srcSEP3D) |
| `application: sep` / `turbulence` | ✓ | ✓ | the application-specific turbulence model (5.4) |
| `application: sep` / `transport`, `population_control`, `storage` | – | ✓ | 3-D mover numerics and coefficients, split/merge limits, storage |
| `mesh` | ✓ | ✓ | `global_cell_size_m`, `minimum_cell_size_m`, `maximum_level` (both); `cells_per_block_edge`, `memory_budget_bytes`, `block_overhead_bytes` (srcSEP3D) |
| `mesh` / `domain` | ✓ | ✓ | `inner_radius_m`, `outer_radius_m` (both); preset, boundaries, frame, origin, corner box (srcSEP3D) |
| `mesh` / `field_line` | ✓ | ✓ | `origin_*_m`, `initial_*_m`, `length_m`, `point_count` (both); `start_mode`, `end_radius_m` (srcSEP3D) |
| `mesh` / `refinement: solar`, `refinement: tube` | ✓ | ✓ | refinement profiles (both); `anchor`, `source_longitude_rad`, `source_colatitude_rad` (srcSEP3D) |
| `mesh` / `active_region`, `memory` | – | ✓ | active corridor and memory planning |
| `physics: background` | ✓ | ✓ | `provider` (srcSEP: `swcme`; srcSEP3D: `analytic-parker`, `swcme`, `runtime-model`, ...); `model_id` (srcSEP3D) |
| `physics: background` / `parker` | – | ✓ | `validity_cadence_s`; the seven ambient keys only without `physics: swcme` (D4) |
| `physics: shock` | ✓ | ✓ | `authority` (srcSEP: `swcme`; srcSEP3D: `none`, `swcme`) |
| `physics: swcme` | ✓ | rejected | srcSEP3D: **the SWCME model is not ready yet to be used** -- any `physics: swcme` block stops the parser at its `#begin` line with that diagnosis (user decision 2026-10-12; the six SWCME example decks are kept as rejected fixtures).  srcSEP: canonical SWCME model (library vocabulary, units in the value); srcSEP rejects the 3-D keys listed in `SwcmeThreeDimensionalKeys()` (CME shape, width, direction and rotation axis `geometry.*`, `surface.*`; the 1-AU thicknesses `geometry.sheath_thickness_1au`/`ejecta_thickness_1au` are shared) and accepts `parker.solar_rotation_rate_rad_per_s` only at its fixed SWCME1D rate |
| `physics: parallel_diffusion` | ✓ | ✓ | shared parallel-diffusion library |
| `physics: shock_front` | – | ✓ | reduced shock-front event, usually `#include`d from its event file |
| `source: shock` | ✓ | ✓ | srcSEP3D: `model` (required, explicit particle-source model: `reduced-shock-front`, the reduced front of `physics: shock_front`, is the only accepted value; `model = swcme` stops with the not-ready diagnosis; an enabled source needs the `physics: shock_front` block).  Both: `enabled`, `samples_per_step`, `minimum_energy_j`, `maximum_energy_j`, `injection_efficiency`, the reference spectrum `spectrum_reference_energy_j`, `spectrum_particle_mass_kg`, `spectrum_charge_number`, `spectrum_normalization`, `spectrum_reference_intensity_per_m2_s_sr_j` (both; both applications currently implement only `spectrum_normalization = relative-only`); rate, spectrum-model and derived-normalization keys (srcSEP3D) |
| `species` | ✓ | ✓ | `macroparticle_weight` |
| `observer: ID` | ✓ | ✓ | `kind`, `minimum_energy_j`, `maximum_energy_j`, `energy_bins`, `energy_spacing`, `pitch_angle_bins` (both); `heliocentric_radius_m` (srcSEP); position, velocity, cadence, normalization, μ range, species, products (srcSEP3D) |
| `output` | ✓ | ✓ | `initialization_mesh_tecplot_file`, `initialization_field_line_tecplot_file`, `initialization_data_tecplot_file` (both); `directory`, `prefix`, `cadence_steps`, `checkpoint_cadence_steps` (srcSEP3D) |
| `sampling: points` (optional) | – | ✓ | `products`, `point_file`; sub-blocks `species: SYMBOL` (`minimum_energy_mev`, `maximum_energy_mev`, `energy_intervals`, `energy_scale`) and `point: NAME` (`position_x_m`, `position_y_m`, `position_z_m`); see 5.7 |

Completeness: every block an application implements is required (srcSEP3D:
all except the optional model blocks and `species` in derived-numerics mode;
srcSEP: `run`, `application: sep` with `turbulence`, the four mesh blocks,
`physics: background`, `physics: shock`, `physics: swcme`, `source: shock`,
`species`, `output` and at least one observer), and every key of a required
block is required except the mode-dependent keys listed in the parsers
(srcSEP: only `run.maximum_time_steps` is optional).

### 5.3 What differs between an srcSEP and an srcSEP3D file

For the same physical event the two files share `run.time_step_s`,
`application.sep.schema_version`/`transport`, the mesh blocks' shared keys,
`physics: background`/`shock`/`swcme`/`parallel_diffusion`, `source: shock`'s
shared keys, `species` and the observers' energy/pitch grids.  They differ in

* the content of `application: sep` / `turbulence` (D2);
* the settings only one application implements (D3): the srcSEP3D numerics
  blocks and keys, and the observer position (`heliocentric_radius_m` versus
  a spatial `kind`);
* `physics: swcme` 3-D geometry keys, which srcSEP's spherically symmetric
  SWCME1D cannot honour.

### 5.4 The turbulence block

srcSEP3D (prescribed spectrum; `srcSEP3D/CONFIGURATION.md`): `authority`,
`model`, `amplitude_model`, `delta_b_over_b`,
`wave_energy_density_at_reference_j_per_m3`,
`wave_energy_density_radial_exponent`, `normalized_cross_helicity`,
`reference_radius_m`, `k_min_per_m`, `k_max_per_m`, `k_min_radial_exponent`,
`k_max_radial_exponent`, `spectral_index`, `correlation_length_m`,
`correlation_length_radial_exponent`, `validity_cadence_s`, `missing_data`,
`resonance_range`, `self_consistent_3d`.

srcSEP (self-consistent Alfvén waves; `srcSEP/README.md`).  Each key has the
command-line option in parentheses, which overrides it:

| Key | Values / unit | Option |
|---|---|---|
| `source` | `prescribed`, `self-consistent-integrated`, `self-consistent-spectral`, `swmf-read-only`, `swmf-initial-then-local` | `--turbulence-source` |
| `representation` | `integrated`, `spectral` (must agree with a self-consistent source) | `--turbulence-model` |
| `coupling` | `disabled`, `streaming-energy-exchange` | `--turbulence-coupling-policy`, `--coupling` |
| `advection`, `reflection`, `cascade`, `shock_injection` | `true`/`false` | `--turbulence-advection`, `--reflection`, `--cascade`, `--shock-injection` |
| `inner_boundary`, `outer_boundary` | `specified-incoming-energy`, `specified-incoming-flux`, `transparent-outflow`, `fixed-reservoir` | `--turbulence-inner-boundary`, `--turbulence-outer-boundary` |
| `inner_boundary_value_si`, `outer_boundary_value_si` | J (energy, reservoir) or W (flux) | `--turbulence-inner-value`, `--turbulence-outer-value` |
| `reflection_coefficient`, `cascade_coefficient` | dimensionless | `--reflection-coefficient`, `--cascade-coefficient` |
| `perpendicular_correlation_length_m` | m | `--turbulence-correlation-length` |
| `electron_heating_fraction` | dimensionless | (file only) |
| `spectral_k_min_per_m`, `spectral_k_max_per_m`, `spectral_bins` | m⁻¹, count | `--spectral-k-min`, `--spectral-k-max`, `--spectral-bins` |
| `cfl_safety`, `operator_accuracy_safety`, `maximum_source_fraction`, `maximum_cascade_fraction` | dimensionless | `--turbulence-cfl`, `--turbulence-operator-safety`, `--turbulence-max-source-fraction`, `--turbulence-max-cascade-fraction` |
| `minimum_substep_s`, `maximum_substeps` | s, count | `--turbulence-min-substep`, `--turbulence-max-substeps` |
| `conservation_relative_tolerance` | dimensionless | `--turbulence-conservation-tolerance` |
| `shock_source_efficiency`, `shock_source_plus_fraction` | [0, 1] | `--shock-turbulence-efficiency`, `--shock-turbulence-plus-fraction` |
| `diagnostic_cadence_steps` | main-loop steps (0 = off) | `--spectral-output-interval` |

Precedence in srcSEP: built-in defaults < input file < command line.  main.cpp
installs the file block as the base layer of the option record and re-parses
the command line on top of it, so every option keeps its usual semantics.

### 5.5 Derived particle numerics (srcSEP3D reduced-front injection)

The reduced-front particle source derives the global step and the statistical
weights instead of taking them from the file:

```text
#begin run
  time_step_mode = derived-from-mesh   ! dt = f * h_min / v_max after mesh allocation
  maximum_particle_speed_m_s = 2.0e8   ! v_max, in (0, c]
  time_step_margin_factor = 0.30       ! f, in (0, 1]
  ...                                  ! run.time_step_s must be absent
#end run
#begin source: shock
  enabled = true
  normalization = accepted-shock-incident-flux
  normalization_radius_m = 1.3914e10   ! apex radius where Ndot_s is evaluated; no default
  samples_per_step = 1000              ! N, expected model particles per step and species
  weight_model = constant-statistical-weight
  maximum_events_per_species_per_step = 1000000   ! fatal Poisson runaway guard
  ...                                  ! physical_particle_rate_per_s, injection_efficiency
#end source                            ! and species.macroparticle_weight must be absent
```

After AMPS allocates the mesh, all owner blocks reduce the smallest
characteristic cell length `h_min` and the one step shared by all species is

```text
dt = time_step_margin_factor * h_min / maximum_particle_speed_m_s .
```

At `normalization_radius_m` a side-effect-free epoch of the reduced front is
evaluated and, for every compiled species `s`,

```text
Ndot_s = sum over accepted curved faces of n_s [V_n - U_1 . n]_+ A_face
W_s    = Ndot_s * dt / samples_per_step .
```

`A_face` is the exact curved SSE area.  Sub-fast, non-forward and
below-support faces contribute zero; numerical-unknown area fails.  A
compiled species absent from the upstream composition fails rather than
borrowing another species' rate.  At a live epoch the constant-weight
macro-event rate is `lambda_s(t) = Ndot_s(t) / W_s`; waiting times are
`-log(U)/lambda_s`, faces are chosen with probability `Ndot_face/Ndot_total`,
points are uniform on the planar triangle and momenta are anti-sunward.  Only
the owning MPI rank allocates.  `log-uniform-momentum-importance` is a typed
reserved mode that stops startup.  Observer cadences are moved up to the first
exact multiple of `dt`.  The incident flux is a normalization convention, not
an acceleration efficiency.  The two derivation keys must be selected
together, and this mode requires `intent = shock-injection`,
`physics.shock.authority = none` and the reduced shock-front background (so
no `physics: swcme` and no reference-spectrum keys).

### 5.6 srcSEP command-line layering

srcSEP keeps its command line (srcSEP/AGENTS.md).  Each option overrides the
file value it shadows: `--particle-mover` overrides `application.sep.transport`,
`--total-iterations` overrides `run.maximum_time_steps`,
`--injection-particles-per-iteration` overrides `source.shock.samples_per_step`,
`--swcme-override KEY=VALUE` overrides one SWCME assignment (an override of
the ambient wind is rejected because the field-line mesh was built from the
file value), and the turbulence options override 5.4.  `--shock-model` and
`--reduced-shock-event` override `physics.background.provider = swcme` /
`physics.shock.authority = swcme`; a post-compile PARAM file that selects the
analytical shock without a command-line option contradicts the input file and
is rejected.

### 5.7 Point sampling (srcSEP3D, `#begin sampling: points`)

Fixed sampling points that record density, flux, differential intensity and
fluence per energy channel **in the AMPS output cycle**: the sampler is
registered with the AMPS user-defined sampling hook
(`PIC::Sampling::ExternalSamplingLocalVariables::RegisterSamplingRoutine`),
samples at every AMPS sampling iteration and writes one Tecplot file per AMPS
output with the same index and the same sampling length as the AMPS mesh
files (`<OutputDirectory>/<output.prefix>.point-sampling.out=N.dat` next to
`pic.<species>.s=<s>.out=N.dat`).  The block is optional; srcSEP rejects it
with the D3 diagnostic (its field-line sampling is a separate, deferred
procedure).

```text
#begin sampling: points
  products   = density, flux, spectrum, fluence   ! any non-empty subset
  point_file = points.txt     ! optional; relative to the file holding this block
  #begin species: H_PLUS      ! one block per sampled AMPS species (symbol)
    minimum_energy_mev = 0.01 ! kinetic energy per particle [MeV]
    maximum_energy_mev = 100
    energy_intervals   = 8
    energy_scale       = log  ! log | linear
  #end species
  #begin point: earth         ! optional named points (block points come first)
    position_x_m = 1.495978707e11   ! [m], AMPS heliocentric Cartesian frame
    position_y_m = 0
    position_z_m = 0
  #end point
#end sampling
```

The point file holds one point per line, `x y z [name]` in metres; `!` and
`#` start comments; a point without a name is named by its zero-based index
among the points of the file.  Every key of the block and of each sub-block
is required (`point_file` is optional); at least one species and one point
are required; species symbols and point names are unique.

* **Units exception.** Channel energies are written in MeV (key suffix
  `_mev`): the one deliberate exception to the SI suffix rule of section 4.
* **Sampling volume.** Only the macroparticles of the AMPS cell containing
  the point are sampled; the cell volume is the sampling volume.  The cell
  centre and edge lengths are written into the product.
* **Quantities** (both the local plasma frame -- the stored transport
  momentum -- and the inertial heliocentric frame obtained by a relativistic
  boost with the cell's bulk velocity): density [cm⁻³], omnidirectional flux
  and field-aligned flux [cm⁻² s⁻¹], differential intensity
  F/(4π ΔE) [cm⁻² s⁻¹ sr⁻¹ MeV⁻¹], interval and cumulative fluence [cm⁻²] and
  the fluence intensity [cm⁻² sr⁻¹ MeV⁻¹]; the cumulative fluence starts at
  the first sample of the run.  The equations are in
  `srcSEP3D/output/point_sampling.h`.
* **Fail-closed.** A species that is not compiled, or a point outside the
  allocated AMPS mesh (or the radial domain, or an inactive region) stops the
  run at start-up.

## 6. Model files

| Model | File | Block | Reader |
|---|---|---|---|
| Reduced shock front | `*.event` (e.g. `srcSEP3D/examples/shock-front/positive_1au.event`) | `physics: shock_front` | `ShockFront::ResolveConfiguration` (`src/models/sep_corona_swcme/shock_front`) |
| CME background event (BG3D) | `event.conf` (`src/models/sep_corona_swcme/examples/*/`) | `physics: cme_background` | `CoronaSwcme::ResolveEventConfiguration` |
| Coronal CME model (schema 5) | `src/models/sep_coronal_cme/test/data/schema5_stage0.in` | `physics: sep_coronal_cme` with one child block per schema section (`#begin observer: earth` for identified sections) | `CoronalCME::ParseConfiguration` |

Checksummed data assets referenced by these models (harmonic CSV, history
CSV, ambient/ejecta/contact `key=value` assets) are byte-identified data whose
SHA-256 enters the physics identity; they are not runtime input files and keep
their frozen data formats.

## 7. Diagnostics and provenance

Input is transactional: nothing is committed unless the complete expanded
input and the immutable configuration are valid.  Errors terminate before
mesh creation and report the actual (included) file, the physical line, the
cause and the offending line, for example

```text
/case/parts/sep3d.in:4: unknown configuration key 'source.shock.particlez'
  line: particlez = 10
```

Fingerprints are computed from normalized values (srcSEP3D/srcSEP option
records, model manifests), not from file bytes, so comments, layout and
include structure do not change a run's identity.  `AMPS::Input::ToCanonicalText`
gives the deterministic expanded text of a block (used to pass an included
model block to its library).

## 8. Migration from the retired formats

All retired SEP inputs -- the srcSEP3D INI decks (schemas 1-5), the srcSEP
INI decks (schemas 1-4), their generic-format successors with
`#begin application: sep3d` / schema 1-5 and `#begin application: sep` /
schema 1-4, the one-dash `#section begin: sep3d` container, the `key=value`
event records and the sep_coronal_cme INI deck -- were migrated; the old
readers were removed.  `tools/convert_amps_input.py` converts old files,
keeping every comment except those attached to a removed key:

```bash
python3 tools/convert_amps_input.py sep6-from-srcsep3d old.in -o new.in   # generic-format srcSEP3D 1-5
python3 tools/convert_amps_input.py sep6-from-srcsep   old.in -o new.in   # generic-format srcSEP 1-4
python3 tools/convert_amps_input.py srcsep3d old.ini -o new.in            # INI srcSEP3D (chains into schema 6)
python3 tools/convert_amps_input.py srcsep   old.ini -o new.in            # INI srcSEP (chains into schema 6)
python3 tools/convert_amps_input.py event    old.event -o new.event
python3 tools/convert_amps_input.py cme-event old.conf -o new.conf
```

The schema-6 rewrite renames `application: sep3d` to `application: sep`,
moves `physics: turbulence` into `application: sep` / `turbulence` (srcSEP
receives the complete block with its former command-line defaults), drops the
srcSEP3D `particle_mover` (the family follows from `transport`) and the
deprecated `external_script`, maps srcSEP's `particle_mover` to `transport`,
removes the SWCME-owned Parker keys and moves the SWCME `source.*` keys into
`source: shock` in SI (MeV → J and `mp` → kg in exact decimal arithmetic),
and renames `initialization_parker_line_tecplot_file`, srcSEP's
`mesh_/field_line_/data_tecplot_file`, `source: field_line` /
`macroparticles_per_step`, `particle_weight` and `energy_channels` to the
shared names.  Values are unchanged except for the unit conversions, which
are exact.

## 9. Current limitations

* Only srcSEP, srcSEP3D and the SEP model libraries use this format so far;
  the planetary applications still configure runs at compile time.
* A coupled SWMF srcSEP run does not read `./amps.in` (only the standalone
  srcSEP executable falls back to it).
* srcSEP and srcSEP3D files for the same event still differ where one
  application has no implementation (5.3); those settings are typed errors,
  not silent omissions.
* The SWCME, parallel-diffusion and model vocabularies are owned and
  documented by their libraries (`src/models/swcme/README.md`,
  `src/models/sep_common/parallel_diffusion/README.md`,
  `src/models/sep_corona_swcme/shock_front/README.md`).
