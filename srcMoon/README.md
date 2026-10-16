# AMPS lunar exosphere application (`srcMoon`)

For the 200 GB memory-budget volume mesh, use
[`examples/m1_memory200.in`](examples/m1_memory200.in) for the spherical
allocation check; [`examples/lola_surface.in`](examples/lola_surface.in)
uses the same coarsened volume targets without changing its terrain/thermal
resolution. See the [memory-profile contract](examples/README.md#200-gb-memory-budget-mesh-2026-10-14).
The spherical memory profile is now approved for the 480/240/120/60 km
qualification campaign; convergence must still be demonstrated, not assumed.

## Scope and evidence

This document describes the code that is present in this checkout, not a
desired lunar model inferred from file names or downloaded data.  The audit was
last rechecked at source revision `92386a37737eab45cc8a855587a881d16bcefcd1`.
Generated files under `build/` were inspected to determine the effective
configuration, but they are build products and must never be patched as the
primary implementation.

The controlling verification contract is
[`development_and_validation.md`](development_and_validation.md).  Its SHA-256
is
`8a9a7f16cd929b4848da869ba9e41de0368e0633a6e444d3d4b38a96d2fb068a`.
The external roadmap at
`/data/vtenishe/moon_validation_data/development_and_validation.md` has the
SHA-256
`d10eca1c539f9cf82940236d9e1aa36af665d43e82ba475925950d582fbf49c4`.
The source copy adds dated implementation-status paragraphs but changes no
test ID, threshold, dataset map, or milestone definition; the controlling
split is recorded in
[`test/ROADMAP_RECONCILIATION.md`](test/ROADMAP_RECONCILIATION.md).

Capability words in this document use the roadmap meanings.  In particular,
source presence, downloaded data, and preflight checks do not establish
observational validation.

## What the current executable is

The default generated `moon.input` requests one species, neutral Na, and builds
`srcMoon/main.cpp`, `srcMoon/main_lib.cpp`, `srcMoon/Moon.cpp`, and the generic
AMPS exosphere model.  With no application `-input` argument, `amps_init()`
retains the regression's analytic sphere.  A `moon` application-input section
can instead select a generated LOLA triangulation.  The production input now
enables SPICE/orbit evolution with the frozen DE421 M0 kernel set.  The
historical `make test_Moon` target deliberately creates a separate no-SPICE
nightly-test configuration; its zero-diff regression therefore does not
exercise the orbit-on branch.

| Item | Effective setting | Evidence and consequence |
|---|---|---|
| target | `_MOON_` | `moon.input` generated definition |
| object/inertial/body/solar-orbital names | `Moon`, `J2000`, `MOON_ME_DE421`, `LSO` | `Moon.h`, `Moon.cpp` |
| species | Na only (`_NA_SPEC_ == 0`) | `SpeciesList=Na`; ion/He/Ne/Ar branches are compiled out or unreachable |
| timestep | species-dependent global | `moon.input`; surface inventory exchange requires this mode |
| particle weight | species-dependent global; no individual correction | `moon.input` |
| reproducible path | on | `ForceRepeatableSimulationPath=on`; this is not by itself a recorded seed contract |
| surface | analytic sphere by default; input-selectable LDEM_4 triangulation | `amps_init()`, `MoonInput.cpp`, `LunarSurface.cpp` |
| production orbit/SPICE | on | `input/moon.input`; requires the exact kernel set and toolkit described below |
| surface temperature | analytic cosine baseline; input-selected processed-Diviner map or depth-resolved thermal model in M2 | the Diviner table and thermal columns are separate production modes; neither local availability nor downloaded D03 files establish I16/I17 |
| Na sticking | generated Yakshinskiy-2005 digitized table | parser rewrites the source default |
| accommodation | per-species input model and factor | M2 Na selects Tenishev et al. (2013) Eq. 14 with `omega=0.2`; legacy no-input runs retain their generated law |
| surface content | legacy prescribed `2.3e16 m^-2` for Na; M2 starts from explicit zero or a versioned checkpoint | the M2 flux-balanced path evolves population; conservation remains an I19 gate |
| active sources | impact vaporization, thermal desorption, and PSD | sputtering off; historical user-defined `MySource` is retained but disabled |
| chemistry | photolytic reactions on; legacy constant Na lifetime | no qualified multispecies network |
| restart | particle restart output requested every 20 iterations | no audited restoration of surface/thermal/sampler/RNG state |

The M0 baseline registers impact vaporization only through the named built-in
process, so `1.69e22 s^-1` is included once.  The historical user-defined
`MySource` block is retained in `input/moon.input` but configured `off`.  The
reported intent of that source is to represent Na that sticks on the nightside
and is later released on the dayside.  Its current callbacks merely alias the
impact-vaporization rate and generator: they do not implement or conserve such
a reservoir.  `MySource` must therefore remain disabled until the inventory,
transport/release law, restart state, and I19 conservation gate are defined and
verified.

## Execution and data flow

1. `main.cpp` initializes MPI, calls `amps_init()`, and repeatedly calls
   `amps_time_step()`.
2. `amps_init()` initializes Moon and PIC services, registers samplers, selects
   either the analytic sphere or generated/loaded LOLA boundary, allocates the
   logical spherical surface state, initializes global timesteps and weights,
   and initializes particle storage and boundary injection.
3. The logical surface dispatches source rates and particle creation through
   `Exosphere::SourceProcesses::totalProductionRate()` and
   `InjectionBoundaryModel()`.
4. The generic injector draws a Poisson arrival process, selects an enabled
   source using its rate, calls that source's production generator, records the
   source ID and origin element, creates the production particle, and moves it
   through the remaining randomized fraction of the timestep.
5. The mover calls `Moon::TotalParticleAcceleration()` through the configured
   production macro.  Lunar point gravity is always active; the SPICE-enabled
   production configuration also activates Na radiation pressure, Sun/Earth
   differential gravity, and LSO centrifugal/Coriolis terms.
6. A sphere hit calls the legacy sphere callback; a terrain-facet hit calls the
   same accommodation kernel with the facet normal.
   A sticking particle is deleted and its statistical weight is added to the
   body-fixed adsorption counter.  A non-sticking particle is re-emitted.
7. `Moon::ExchangeSurfaceAreaDensity()` aliases the generic
   `Exosphere::ExchangeSurfaceAreaDensity()`: MPI ranks reduce adsorption minus
   desorption, update the population, recompute surface source rates, then zero
   the flux counters.  Under the current user-defined surface-content mode the
   population is then overwritten by `area × 2.3e16 m^-2` on every exchange.
8. `amps_time_step()` updates SPICE/frame state only when orbit mode is compiled
   on, then calls `PIC::TimeStep()`.

## Frames, time, and units

- AMPS positions are metres, velocities are metres per second, accelerations
  are metres per second squared, temperatures are kelvin, number densities are
  per cubic metre, surface densities are per square metre, and source rates are
  particles per second unless the called API states otherwise.
- `J2000` is the inertial frame used for ephemeris states. `SO_FRAME` is the
  Moon-centred `LSO` rotating frame whose +x direction points toward the Sun
  and is used by analytic surface and shadow logic. The shared Exosphere
  variable named `IAU_FRAME` contains `MOON_ME_DE421`, the DE421 mean-Earth /
  polar-axis body-fixed frame used by the selected LOLA product.
- Six-by-six SPICE transformations carry both position and velocity.  Surface
  impacts transform SO to IAU before body-fixed binning and re-emission, then
  transform back to SO.
- The orbit-on branch advances ephemeris time by the smallest global species
  timestep. External states and the LSO angular velocity are frozen at the
  step midpoint; body-frame transformations are stored at the new endpoint.
  Force states use geometric aberration correction `NONE`. The former
  Mercury-specific `MSGR_HCI` dependency has been removed from this path.

### Frozen M0 SPICE contract

The M0 frame test uses `2009-01-24T00:00:00 UTC`. Configure the checkout with:

```sh
rm -rf build
./Config.pl -application=moon \
  -spice-path=/home/vtenishe/SPICE/cspice \
  -spice-kernels=/home/vtenishe/SPICE/Kernels
make -j
```

`input/moon.input.spicekernels` furnishes exactly `de421.bsp`,
`moon_pa_de421_1900-2050.bpc`, `moon_080317.tf`, `moon_assoc_me.tf`,
`naif0012.tls`, `pck00011.tpc`, and `Moon.LSO.tf`. U04 freezes and verifies the
SHA-256 of every file before it calls CSPICE. A missing file is `SKIPPED`; a
hash mismatch is `ERROR`; neither condition is a physics PASS.

## M2 production surface/environment state

An M2-enabled receipt promotes the generated LOLA triangulation into the
ordinary particle mover. For every particle segment, the boundary callback
evaluates the MOON_ME_DE421-to-LSO state transform at the collision epoch,
brackets the first outside-to-inside terrain crossing using the configured
maximum segment length, bisects it to the stated time/distance tolerances, and
returns both the LSO facet normal and moving-facet velocity. Every MPI rank
builds the same read-only triangulation/BVH; rank zero alone writes requested
CEA/Tecplot surface products.

The terrain BVH also supplies the production horizon test. Terrain blocking
and the Earth-eclipse cylinder are separate classifications. D02
AVGVISIB/LPSR products never construct the model horizon and remain independent
I15/I20 evidence.

`surface_temperature_mode` selects `simple-cosine`, `diviner-map`, or
`thermal-model`. The Diviner loader accepts only its documented processed CSV
schema and rejects duplicates, missing/non-GOOD interpolation corners, and
out-of-support latitude. The predictive mode advances one finite-volume 1-D
column per logical surface element using the cited Hayne et al. density,
conductivity, heat-capacity, albedo, emissivity, and geothermal parameters. It
has nonlinear convergence guards and an exact enthalpy-form energy audit. A
climatology is never substituted for a transient eclipse history.

Na and Ar adsorption update a continuous flux-balanced surface inventory.
Thermal source rates are stored population divided by the current residence
time. Cold-trap status is a diagnostic comparison with the configured
timescale, not an irreversible sink; warming increases release continuously.
Initial inventory is exactly zero unless a checkpoint is selected. Checkpoints
are independent `<base>.Na` and `<base>.Ar` files, and restart fails if either
is missing. Neutral laws are not applied to Na+/Ar+ while ion branching is
unconfigured.

Immediate non-sticking Na impacts use the explicitly selected kinetic-energy
law from Tenishev et al. (2013), Eq. 14,
`vf=sqrt(omega*vt^2+(1-omega)*vi^2)`, with the paper's evaluated
`omega=0.2` and isotropic outgoing direction. Sticking is an independent
Bernoulli event using the input-resident Fig. 3 table. Both the sticking model
and accommodation definition/factor are named in the `.in` receipt; the M2
adapter never inherits the legacy linear-speed blend implicitly. This resolves
the former Na-law configuration blocker. I18 now exercises both laws through
production callbacks in a separately built Na+Ar executable, but a
qualification result still requires I13 PASS in the same campaign. See
[`SURFACE_INTERACTION_LITERATURE.md`](SURFACE_INTERACTION_LITERATURE.md).

Species-specific new code is kept under `srcMoon/species/`. The current
`species/Argon.cpp` module owns the optional radiogenic `40Ar` source receipt,
reference-sphere area partition, source-rate callback, and local-temperature
Maxwell-flux particle generator. The normal Na-only Config.pl source block is
off, preserving the M0 baseline. An experimental Ar build must explicitly add
Ar to `SpeciesList`, enable only `RadiogenicArgon`, and set
`argon_source_enabled=true` in its runtime receipt. The maintained receipt
selects the Killen (2002) lower bound `3.8e20 atoms s^-1`; it does not use the
Kegerreis fitted localized rates as independent source physics. D10 remains a
regional, QA-incomplete spatial proxy, so potassium-weighted geography is
rejected and U21 remains `SKIPPED / NOT VALIDATED` despite its local uniform-
source kernel check.

## Mesh and surface geometry

### Production surface selection

`srcMoon` accepts an optional application file:

```sh
./amps -input srcMoon/examples/lola_surface.in
```

The file uses the repository-wide grammar implemented by
`src/general/input/amps_input.h`. It must contain exactly one flat
`#begin application: moon` / `#end application` block. Keys are
case-sensitive, assignments are `name = value`, `!` begins a comment, and
duplicate, unknown, or missing settings are errors. Unowned root blocks,
nested Moon blocks, and the retired `#section` spelling are rejected rather
than silently ignored. Relative paths are resolved against the file containing
the assignment, including when the generic `#include` directive is used. The
required settings are:

```text
spice_path
surface_geometry = sphere | lola
surface_mesh_resolution_m
lola_product_id
lola_image_file
lola_label_file
surface_cea_file
surface_tecplot_file
mesh_profile = quick | full | custom
mesh_near_surface_cell_size_m
mesh_near_surface_extent_m
mesh_far_field_cell_size_m
mesh_far_field_extent_m
mesh_outer_cell_size_m
mesh_output_tecplot = true | false
mesh_angular_resolution_m_by_sza
mesh_terminator_half_width_deg
mesh_polar_cap_start_latitude_deg
mesh_polar_resolution_m
simulation_start_utc
simulation_environment_mode = frozen | evolving
random_seed = 731993
profile_altitude_near_maximum_m
profile_altitude_near_step_m
profile_altitude_middle_maximum_m
profile_altitude_middle_step_m
profile_altitude_far_maximum_m
profile_altitude_far_points_per_decade
profile_solar_zenith_angles_deg
profile_include_body_fixed_poles
profile_include_shell_average
profile_shell_angular_weighting = solid-angle
limb_observer = earth
limb_plane = sun-moon-earth
limb_sides = both | toward-projected-sun | away-from-projected-sun
limb_branch_reference = projected-sun
limb_tangent_maximum_altitude_m
tail_center_direction = anti-lunar
tail_map_coordinate_system = j2000-ra-dec | eclipj2000-lon-lat
tail_map_display_orientation = west-right | east-right
tail_map_half_width_deg
tail_beam_diameter_deg
tail_exposure_duration_s
tail_beam_response = unresolved
tail_temporal_operator = exposure-average
```

`spice_path` must contain `cspice/include/SpiceUsr.h`,
`cspice/lib/cspice.a`, and `Kernels/`, matching `/home/vtenishe/SPICE`. It is
recorded and checked at runtime; it does not override the build-time SPICE
selection or kernel root supplied to `Config.pl`.

The FULL M1 production example (`srcMoon/examples/m1_sphere.in`) requests
10 km near-surface volume cells through
200 km altitude, transitions to 40 km at 10,000 km altitude, then coarsens to
17,371 km at the existing 400-lunar-radius outer boundary. The outer scale is
explicit because applying 40 km cells throughout that three-dimensional tail
box is computationally infeasible. M1 claims are limited to the common
0--10,000 km science grid; later tail work requires a separate convergence
study. Its
boundary table is
`0:10000,30:10000,60:5000,90:5000,120:10000,180:10000`; a 15 degree band
around the terminator and body-fixed latitudes at or above 75 degrees are
limited to 5 km. These are explicit input values, not compiled defaults. SZA
is an LSO angle. Polar membership is transformed to `MOON_ME_DE421` at the
configured epoch. The selector returns requested metres; linked I12 separately
records the power-of-two AMR size actually realized.

The QUICK profile (`srcMoon/examples/m1_quick.in`) keeps those physical
extents and every non-mesh control unchanged but uses three times each FULL
spatial target: 30 km near-surface, 120 km at 10,000 km altitude, 52,113 km at
the outer boundary, and 15 km minimum angular/polar resolution. The original
factor-two QUICK candidate built 1,074,788 blocks; its eight-rank I12 run was
stopped after 28.7 minutes at 223.7 GiB RSS with only 5.7 GiB host memory left.
Factor three crosses the next coarser dyadic AMR level. Its cold eight-rank
I12 run completed the wiring criteria in 1,409.9 seconds (23.50 minutes), with
191.9 GiB peak RSS and 39.8 GiB minimum host-available memory, on `passion` on
2026-10-12. The required `CenterNodeMeasure.bin` cache was 6,128,421,120 bytes;
optional `mesh.dat`/`tetra_mesh` files were absent. Runtime is host/cache-
dependent operational guidance, not an acceptance criterion. QUICK is for
production-wiring checks only; M1 numerical qualification continues to require
the FULL I11-I13 campaign.

The still-coarser `srcMoon/examples/m1_regression.in` receipt uses twelve
times the FULL spatial targets while preserving the physical domain, epoch,
sources, profile grid, and observation geometry. It is intentionally marked
`custom` and is the default bounded eight-rank linked run in
`srcMoon/test/run_all_test.py`; successful wiring is therefore retained as
`SKIPPED / NOT QUALIFIED`, never as I11--I13 evidence.

All target sizes are ordinary `.in` values rather than compiled defaults.
Copy a committed receipt, edit the radial, angular, polar, and (when relevant)
surface target sizes, and label the result `mesh_profile = custom`. Changing a
resolution target never inherits an undocumented change to a physical extent.

`mesh_output_tecplot` controls only the optional full-volume `mesh.dat` and
`tetra_mesh*.dat` visualization products. It is required in every configured
receipt and is `false` in the committed science examples because the text dump
can require hundreds of gigabytes without changing any modeled physics. The
AMR tree, particle mover, and samplers are identical in either setting. The
internal `mesh.msh` exchange used by the current MPI initialization and the
separate triangulated `surface_tecplot_file` are not controlled by this flag.

### M1 comparison and observation geometry controls

The runtime receipt also defines one run-wide UTC and the common grid on which
future I11–I13 convergence products must be compared. Altitude is measured
from the AMPS `_RADIUS_(_TARGET_)` reference sphere (1737.10 km), not from the
1737.4 km LOLA datum or local relief. The committed example
uses 0–200 km in 5 km steps, 200–1000 km in 25 km steps, and 1000–10000 km at
eight logarithmic points per decade. The parser constructs that exact grid in
metres and rejects non-integral linear endpoints, unordered/duplicate SZA
values, and a limb maximum that is not a grid node. Values above 10000 km are
reserved for integrated LOS/tail outputs; this configuration does not invent
pointwise observations there.

`simulation_start_utc` replaces the generated `SimulationStartTimeString` for
the whole CLI-configured run before `Exosphere::Init_SPICE()` converts UTC to
ET. This preserves a single physical clock for forces, rotation, illumination,
and diagnostic geometry. Omitting `-input` retains the compile-time epoch for
the historical regression.

`simulation_environment_mode=frozen` is the approved I11 steady-state
contract. Particle transport and source/surface processes continue with their
normal SI time steps, while SPICE states, illumination, body transforms, and
rotating-frame coefficients remain at `simulation_start_utc`. This isolates
spatial discretization from mesh-dependent epoch drift. `evolving` retains
normal midpoint ephemeris evaluation and advances ET once per particle step;
it is not used to score I11 unless matched physical time windows are defined.

`random_seed` initializes the production AMPS random stream before sources,
surface interactions, and transport begin. I11 requires `731993` in every
derived 40/20/10/5 km receipt, records one positive MPI layout for the complete
campaign, and fixes one OpenMP thread per rank.

The linked I11 sampler reads completed production density and surface buffers.
Its vector contains the six LSO SZA profiles, transformed north/south
`MOON_ME_DE421` pole profiles, solid-angle shell averages, six fixed 30-degree
SZA surface-return-flux bands, global return rate, finite-distance Earth limb column
density and both Na D brightnesses on both projected-Sun branches, and global
inventory. Batch duration doubles until all lag-one correlations are within
`1.96/sqrt(N)`. Four independent means establish three consecutive <=2%
equilibrium changes; five following disjoint means supply Student-t 95%
intervals. Runs have no wall-time or memory cap. Resource histories are
archived; actual allocation failure terminates through the AMPS fatal path
with resource-reduction guidance and cannot produce a PASS.

For I11 resolution derivation, the surface, near-surface, SZA/angular, and
polar targets scale through 40/20/10/5 km, while far-field and outer-domain
targets remain fixed. This isolates near-Moon convergence without multiplying
remote tail cells. After load balancing and before block allocation, srcMoon
also sums a conservative unavoidable block-object memory lower bound across
MPI ranks sharing each host. It terminates only if that lower bound already
equals or exceeds Linux `MemAvailable`; no empirical cap is imposed.

The agreed profile directions are SZA 0, 30, 60, 90, 120 and 180 degrees in
the LSO equatorial plane, the north/south `MOON_ME_DE421` pole directions, and
a separate solid-angle-weighted shell average. Solid-angle weighting is
physical; counting AMR cells or unweighted nonuniform facets is not an
equivalent shell average. Near-Moon LOS geometry uses an Earth observer, the
Sun-Moon-Earth plane, and the common grid through 1000 km. `limb_sides` may
select both branches or either branch named relative to the projected
Moon-to-Sun direction. A zero projected-Sun direction is a geometry
degeneracy, not permission to choose a branch silently.

Tail geometry uses a ±7.5 degree Earth-view map centred on anti-lunar with a
1 degree circular field and 120 s exposure. Anti-lunar is explicitly not
treated as anti-solar. The recommended storage system is native J2000 RA/Dec,
with east-positive RA offsets and north-positive Dec offsets. ECLIPJ2000
longitude/latitude is an adjustable alternative because D14 contains those
independently transformed coordinates. `west-right` reproduces the usual sky
display convention but does not negate or reorder stored coordinates.

D14 identifies the aperture diameter, shape, and exposure duration, but it
does not provide the WHAM radial beam-response function. The required input is
therefore `tail_beam_response=unresolved`: no Gaussian or uniform-disk response
is assumed. `tail_temporal_operator=exposure-average` records the physically
required time integration. The full raster spans hours and must ultimately be
evaluated at native pointing times rather than at the 2009 M0 baseline epoch.

These fields are now strict, adjustable input/provenance controls. The existing
subsolar-limb production sampler consumes the exact common altitude nodes only
through `limb_tangent_maximum_altitude_m`; the no-input regression retains its
historical geometric grid. The legacy samplers still do not produce both limb
branches, consume the pole/shell selections, or perform finite-beam and
exposure-time convolution, so parsing and one-limb grid wiring are not an
observation-operator PASS. U26 and the applicable I-series test remain
`SKIPPED / NOT VERIFIED` until the complete production operators and
independent geometry oracles are implemented.

For `surface_geometry=sphere`, the registered analytic sphere is unchanged.
For `surface_geometry=lola`, rank zero reads the native detached PDS IMG/LBL
pair, constructs an icosphere, refines it until the maximum great-circle edge
on the reference sphere is no larger than `surface_mesh_resolution_m`, and
bilinearly samples the native DEM at every vertex. An icosphere avoids the
polar edge collapse of a latitude/longitude mesh and gives approximately
uniform global edge lengths. No smoothing, gap-filling, or rescaling is
applied.

The triangulation is atomically written in CEA long format for the normal AMPS
loader and Tecplot `FETRIANGLE` format for inspection. Every MPI rank then
reloads the CEA file through `ReadCEASurfaceMeshLongFormat()`. That
triangulation is the sole geometric collision boundary. The existing 60×100
spherical grid remains only as the source and surface-inventory index. Source
locations are projected radially onto a terrain facet and emission velocities
are rotated from the radial normal to the facet normal. Impact
sticking/re-emission uses the facet normal but retains the established
body-fixed spherical inventory bins.

### Surface implementation invariants

The following constraints are part of the implementation contract and should
remain visible in code review:

- Parsing is transactional and occurs before MPI initialization. A failed
  input never installs a partial configuration.
- Relative paths are anchored to the application input file. They must not
  depend on the directory from which `mpiexec` is launched.
- The detached PDS label controls product ID, byte order, dimensions, scale,
  datum, projection, longitude direction, and coordinate-system name. A new
  product must extend the metadata checks rather than reuse LDEM_4 assumptions.
- `surface_mesh_resolution_m` limits great-circle edge length on the reference
  sphere before relief is applied. It is not an AMR cell size, DEM pixel size,
  triangle chord, or post-relief geodesic tolerance.
- Shared midpoint caching makes the icosphere watertight. Every face must have
  outward winding before it reaches the AMPS CEA loader.
- Only MPI rank zero reads the raster and writes surface products. All ranks
  load the completed CEA file after a collective success check and barrier.
- The registered triangulation is the mover/collision boundary. The
  unregistered 60×100 sphere is bookkeeping state for existing source and
  surface-population kernels; it must never be cast from the first registered
  boundary in LOLA mode.
- Terrain injection and interaction adapt geometry only. Rates, source
  distributions, sticking probabilities, accommodation, and inventory updates
  continue to use the shared production exosphere kernels.
- Candidate injection ownership must be recomputed after radial projection,
  because the final terrain point may lie in another AMR block or MPI rank.
- A re-emitted or newly injected particle is displaced outward by the AMR
  geometric tolerance. A point exactly on a facet has ambiguous parity and
  must not be handed back to the mover unchanged.
- Realistic terrain remains static. Orbit-on builds are rejected until an
  explicit, tested transformation/moving-boundary design exists.

Realistic mode is restricted to the static orbit-off build because the AMR
triangulation is stored in the body-fixed LDEM axes. Moving a cut-cell surface
during a run is not implemented. Terrain shadowing is also not implemented by
this change.

### LOLA frame report

The local `LDEM_4.LBL` controls interpretation. It declares `BODY-FIXED
ROTATING`, east-positive longitude, simple cylindrical projection, and
coordinate-system name `MEAN EARTH/POLAR AXIS OF DE421`. It contains 720×1440
signed 16-bit little-endian samples at four pixels per degree, a 0.5-metre
scale, and elevations relative to a 1,737,400 m sphere. Pixel centres are used;
longitude wraps at 0/360 degrees and latitude clamps at the centres of the
polar rows.

The corresponding NAIF DE421 mean-Earth lunar frame is `MOON_ME_DE421`.
`LunarSurface.cpp` does not silently substitute `IAU_MOON`: it consumes the
label-native longitude/latitude axes directly. A future orbit-on surface must
declare and apply the transform between this frozen DE421 mean-Earth frame and
the run's selected lunar orientation frame before it can be enabled.

Authoritative references are the [LOLA RDR Software Interface
Specification](https://imbrium.mit.edu/DOCUMENT/rdrsis.htm) and the [NAIF lunar
frames tutorial](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/Tutorials/pdf/individual_docs/23_lunar-earth_pck-fk.pdf).

`main_lib.cpp::localResolution()` supplies the volume callback.  The current
surface callback computes a subsolar angle and then forcibly sets it to zero;
therefore `localSphericalSurfaceResolution()` is constant at
`0.16 R_Moon = 277936 m` in this build.  At the lunar-radius U07 fixture,
`localResolution()` returns `8 R_Moon = 13896800 m`.  These values document the
compiled baseline; they are not claims of spatial adequacy.

`srcMoon/Topography/lola.cpp` is the superseded conversion prototype. It reads
hard-coded longitude, latitude, altitude, and sphere-mesh filenames, uses a
nearest grid cell, converts altitude from kilometres to metres, radially moves
mesh vertices, and writes a CEA mesh.  `srcMoon/Topography/shadow_calc.cpp` is
also a stand-alone program: it reads that mesh, loads hard-coded SPICE paths and
kernels, ray-traces sunlight, and writes shadowed meshes. The production code
does not call either prototype. `shadow_calc.cpp` remains outside production,
so its hard-coded inputs are not a reproducible D02 pipeline.

The production interaction normal is radial in sphere mode and the loaded
facet normal in LOLA mode. Terrain horizons, self-shadowing, and permanent
shadow regions still do not affect the model.

## Forces and lunar rotation

`Moon::TotalParticleAcceleration()` contains these branches:

- Lunar gravity is always evaluated as
  `a = -G M_Moon x / |x|^3`.
- Na radiation pressure is present only in an orbit-on build and is suppressed
  in the cylindrical lunar shadow and by `EarthShadowCheck()`.  The kernel uses
  heliocentric radial velocity and inverse-square heliocentric distance.
- Charged species receive `q(E + v × B)/m`.  With the coupler off, the fields
  come from the configured typical solar wind; otherwise they come from the
  coupler interpolation stencil.  No charged species is active now.
- Orbit-on builds add differential solar gravity, differential terrestrial
  gravity, centrifugal acceleration, and Coriolis acceleration.

The orbit-on update obtains lunar, solar, and terrestrial geometric SPICE
states in the declared frames. It gets the instantaneous LSO angular velocity
from the derivative in the `LSO`→`J2000` state transformation and fills
LSO↔`MOON_ME_DE421` transforms. This removes the former small-angle
`acos`/`sin(angle)` singularity. U04 independently reconstructs
`W=(dR/dt)R^T`, checks rotation orthogonality, and checks each differential,
centrifugal, and Coriolis component. This verifies the local production
kernels and frame convention; the full linked trajectory gates I06/I10 remain
required before M0 can be called complete.

`EarthShadowCheck()` implements a cylinder of Earth radius extending
anti-sunward from Earth; it is not a conical umbra/penumbra model.

The original machine-readable digitization behind the embedded Combi (1997)
Na radiation-pressure array is unavailable. U05 therefore compares the
compiled production function with a separately generated two-resolution
digitization of the published Figure 7 curve. Agreement is assessed at 14
unobscured velocities using a frozen `0.5 cm s^-2` graphical uncertainty. This
qualifies the local curve and unit conversion at publication resolution; it
does not establish exact historical table provenance or the linked mover and
shadow dispatch required by I07.

## Sources and particle injection

The configured Na source values are:

- impact vaporization: `1.69e22 s^-1`, 6000 K, zero heliocentric power index;
- thermal desorption: binding energy `1.85 eV`, vibrational frequency
  `1.0e13 s^-1`, multiplied by the surface population and the Boltzmann factor;
- photon-stimulated desorption: photon flux `2.0e18 m^-2 s^-1`, cross section
  `3.0e-25 m^2`, configured injection-speed interval 10–10000 m/s;
- solar-wind sputtering: configured off.

The historical `MySource` user-defined registration is also configured off.
Its retained callback settings currently alias impact vaporization and must not
be interpreted as verified night-to-day surface transport or delayed release.

The injector uses production particle weights and global timesteps.  Source
events follow exponentially distributed inter-arrival times.  Source process,
surface element, position within an element, energy/temperature distribution,
direction, and the fraction of a timestep already moved are stochastic.  The
source ID, origin element, injected weight/rate, and velocity statistics are
sampled by production code.

`Exosphere_Helium.cpp` contains a dayside projected-flux generator and samples
a Maxwellian at the analytic surface temperature.  Its header fixes the alpha
particle fraction to 5 percent of solar-wind proton flux.  The source is not
enabled, He is absent from the species table, and no complete reservoir balance
or authoritative driver is established.  `Exosphere_Neon.cpp` contains only a
dayside rate calculation in this tree; Ne is absent.  No production code for
radiogenic Ar geography/transients, D11 meteoroid forcing, or H2O/OH
migration/chemistry was found.

## Surface interaction, sticking, and release

At impact, the production callback:

1. transforms the state from LSO to `MOON_ME_DE421` (the shared arrays retain
   their historical `IAU` names);
2. computes analytic surface temperature;
3. samples return flux and impact speed;
4. draws a Bernoulli sticking decision;
5. if sticking, adds `reemission_fraction × particle_weight` to the body-fixed
   adsorption counter and deletes the particle;
6. otherwise draws a thermal speed and random outward direction, applies the
   per-species input-selected linear-speed or kinetic-energy accommodation
   law, converts back to LSO, and retains the particle.

For M2 Na, the input selects bounded linear interpolation through the frozen
Tenishev et al. Fig. 3/Yakshinskiy-Madey digitization (100–495 K, with endpoint
clamping); the 100 K control point is 0.9983. M2 Ar selects an explicit
constant-unity table. The legacy no-input build retains its generated Na/Ar
functions, including the older Ar piecewise base-10 regression controls at 88,
110, and 158 K. He and Ne surface models remain later-roadmap capabilities and
are not silently assigned Na or Ar parameters.

The retained analytic baseline is 100 K on the nightside or in terrestrial
shadow and `100 + 280 cos(theta)^(1/4) K` otherwise.  M2 can instead select a
strict processed-Diviner interpolation table or a depth-resolved nonlinear
thermal column advanced with production terrain and eclipse illumination.
Those modes are not interchangeable: the former is an observational boundary,
whereas the latter predicts physical kinetic temperature.  I17 converts the
thermal prediction through `epsilon sigma T^4` before comparison with Diviner
bolometric brightness temperature; it never compares those two temperatures
as though emissivity were unity.

The predictive eclipse path no longer switches the Sun discontinuously with
the legacy infinite-cylinder predicate. It evaluates the analytical overlap of
two circular apparent discs at each body-fixed surface point using SPICE
Sun/Earth vectors and the input-declared radii, then supplies the unocculted
fraction to the same nonlinear thermal boundary. Terrain occlusion is still a
separate all-or-none ray to the solar centre. The same centre direction
supplies the incidence cosine; neither incidence nor terrain penumbra is
integrated across the finite solar disc. The approved first-comparison receipt assumes a
uniform-radiance solar disc, an opaque spherical solid Earth, and no refracted,
scattered, or transmitted terrestrial-atmosphere light. It uses the IAU nominal
solar radius (695,700 km) and the AMPS mean solid-Earth radius (6,371 km), both
visible in the `.in` file rather than hidden constants. These limitations are
material during totality, so every resulting plot is labeled
`EXPERIMENTAL / NOT VALIDATED`.

I17 initializes the complete 60x100 logical-surface set of 80-layer production
columns without observations. It replays the input-declared pre-eclipse
forcing interval until every boundary/layer temperature and each column's
integrated sensible enthalpy meet the two declared periodic tolerances. The
maintained receipt uses a 2,551,442.976 s mean-synodic interval ending at the
first June training bin, 0.01 K maximum node change, and 10,000 J m^-2 maximum
column-enthalpy change. Repeating one ephemeris interval assumes periodic lunar
forcing even though real successive lunations differ; it avoids using Diviner
temperatures as an initial condition but remains an explicit approximation.
There is no cycle-count cutoff. Allocation failure or externally exhausted job
resources terminate as an error with the preserved resource history.

The generic exosphere library has flux-balanced adsorption/desorption machinery
and thermal/PSD source depletion counters.  In the active configuration,
however, the user-defined fixed surface density overwrites the evolved Na
population after applying net flux.  Thus it does **not** implement a conserved
accumulating reservoir, a physical residence-time distribution, PSR trapping,
or delayed release.  Those capabilities remain experimental/missing and must
not be inferred from the allocation of surface arrays.

## Chemistry and plasma

The active photochemical lifetime function accepts Na only and, with orbit mode
off, returns the constant
`3600 × 5.8 / 0.4^2 = 130500 s` with reactions allowed everywhere.  In an
orbit-on build it additionally gates reactions using lunar and terrestrial
shadow geometry.  The generic photochemical processor is selected.  Whether a
daughter ion can be represented depends on the species table; Na+ is absent in
the active table, so a local lifetime PASS is not linked neutral-to-ion
conversion verification.

`Moon::ElectronImpactIonizationRate()` contains fixed `2.0e-20 m^2` cross
sections and a fixed 400 km/s speed for Na, Ne, and Ar, multiplied by the
typical solar-wind density.  It is not wired as a qualified D05 rate model.
No configured charge-exchange reaction path was found.  A dormant BATSRUS data
file/coupler branch exists in `main_lib.cpp`, but it is not a qualified ARTEMIS
D06 plasma driver and the current D06 QA status is ERROR.

## Sampling and observation operators

The code registers:

- cell density/source-resolved sampling from the generic exosphere model;
- body-surface source, return, sticking, speed, and population diagnostics;
- LOS column density and density-weighted mean-speed integrands;
- Na D-line brightness using embedded g-factor functions and shadow gating;
- subsolar-limb and anti-solar column maps;
- velocity distributions at configured points;
- Kaguya TVIS geometry tables and output routines.

SPICE-dependent observer geometry is compiled into the production orbit-on
configuration, but it has not passed an independent observation-operator
geometry test. The embedded Kaguya tables and code presence are not operator
verification and are not D07 validation. No qualified production operators
for LADEE NMS, LAMP, LACE, PACE, or D15 water events were found.  Before any
observational score, each operator needs an independent geometry/units test and
all predecessor gates required by the roadmap.

## Component inventory and verification map

The table records production locations/callers, controls and I/O, present
coverage, and the first uncovered failure mode.  External package labels refer
to the roadmap; package readiness must still be checked from its README,
provenance, QA report, and hashes before use.

| Component | Production source / caller | Switches, inputs → outputs; units/frame | Data / stochastic state | Coverage and first gap |
|---|---|---|---|---|
| initialization/configuration | `main.cpp`, `main_lib.cpp::amps_init`, `Moon::Init_AfterParser` | `moon.input` → callbacks, mesh, weights, samplers | parser-generated tree | U02 source guard; linked baseline I01-I10 pending |
| species/mass/charge/weights/timestep | `species.input`, PIC molecular tables, `localTimeStep`, PIC weight initialization | species; cell size → seconds and statistical weights | Na only; RNG enters weight/injection use | no complete local contract; multi-species paths disabled |
| volume AMR/refinement | `localResolution`, mesh init/build; optional `outputMeshTECPLOT` | metres in LSO → target cell size; explicit `mesh_output_tecplot` gates only full-volume visualization | deterministic targets; I11 uses seed 731993 and MPI=1/OpenMP=1; text visualization may require hundreds of GB | U07 callback and linked I11 sampler/I12 realization gate implemented; I11 execution evidence and I13 reviewed freeze remain pending |
| spherical surface | `amps_init`, `cInternalSphericalData` | 60×100, lunar radius → areas/normals/elements | deterministic geometry | U07 arithmetic; intersections/area invariants pending U01/I14 |
| LOLA triangulated geometry | `MoonInput.cpp`, `LunarSurface.cpp`, `MoonSurfaceEnvironment.cpp`, `amps_init` | native IMG/LBL and requested max edge m → body-fixed terrain/BVH plus optional CEA/Tecplot mesh | D01 package READY with frozen raw hashes; moving geometry requires SPICE | U08 PASS; orbit-enabled moving intersection is wired; linked I14 remains predecessor-gated |
| intersections/normals | `MoonSurfaceEnvironment.cpp` called by the second-order mover | LSO segment and TDB epoch → hit time/position, LSO and body normals, facet velocity | deterministic BVH and SPICE state transforms | U09 analytic ray fixtures; linked I14 collision/convergence gate pending |
| illumination/Earth eclipse/PSR | production terrain BVH and separate Earth-eclipse test; thermal/source callers | body-fixed facet/Sun and SPICE epoch → separate terrain/eclipsed flags | D02 raw acquisition complete but package PARTIAL | U09 PASS; independent I15/I20 comparison SKIPPED |
| simple temperature | `Exosphere::GetSurfaceTemperature` | cosine, LSO position → K | deterministic | U12 linked probe PASS-capable |
| Diviner/dynamic thermal state and observation | `MoonSurfaceTemperature.cpp`, `MoonThermalModel.cpp`, `MoonDivinerObservation.cpp`, ordinary time-step update | processed D03 or forcing/thermal parameters → kinetic K, depth-column state, finite-disc direct irradiance, epsilon-sigma-T^4 and bolometric brightness K | GCP/PRP and frozen June/December RDR files exist; observation-free periodic spin-up and a June-only experimental comparison path are implemented, while D03 package metadata and a reviewed/frozen pre-holdout threshold remain incomplete | U09--U11 local gates cover overlap/solver/operator mathematics; linked I17 can emit June CSV/PNG/EPS only after I13+D03 readiness and remains SKIPPED / NOT VALIDATED until threshold freeze and untouched hold-out scoring |
| lunar gravity | `Moon::TotalParticleAcceleration`; mover macro | m → m s^-2, LSO | deterministic | U03 linked analytic probe; trajectory I05 pending |
| Sun/Earth differential and rotating forces | `Moon::OrbitalDynamics` helpers; `amps_time_step` populates state | DE421 at fixed UTC; J2000/LSO, m, m/s, rad/s → m/s² | exact seven-file kernel set and hashes frozen in U04 | U04 linked-kernel probe PASS-capable; full executable I06/I10 pending |
| Na radiation pressure | same kernel plus `Na.h` table | heliocentric radial speed/distance → m/s² | embedded table; Combi Figure 7 digitization and hashes frozen in U05 | U05 local curve/scaling/shadow qualified to graphical uncertainty; linked mover gating I07 pending |
| Lorentz/ion mover | same kernel; typical or coupled E/B | q, kg, V/m, T, m/s → m/s² | no active ion; D06 ERROR | U06/U18 SKIPPED |
| surface impact classification | generic surface-interaction callback plus `species/Argon.cpp` source generator | species/state/weight → delete or re-emit | Bernoulli RNG and Maxwell/Maxwell-flux draws | U12 PASS; I18 linked Na+Ar sampler implemented and awaiting the same-campaign I13 gate |
| sticking/accommodation/re-emission | `Moon.cpp` plus generic callback | K/probability; incident velocity → outward m/s | RNG for decision/speed/direction | U12 control points; distribution/conservation pending |
| surface inventory/release | runtime Exosphere flux-balance callback plus `MoonSurfaceTemperature.cpp` checkpoints | weighted adsorption/desorption per element → continuous Na/Ar population and source rates | MPI reduction; zero/checkpoint initial state | U12/U13 PASS; complete linked counters and I19 budget closure pending |
| cold trapping | `MoonSurfacePhysics.cpp`, production temperature/inventory service | residence time and explicit comparison timescale → diagnostic classification; warming releases continuously | D02/D03 packages not fully qualified | U13 PASS; I20 SKIPPED |
| photoionization/dissociation | `ExospherePhotoionizationLifeTime`, generic processor | species/position → lifetime s/reaction | reaction RNG in linked model; D04 QA PASS | U14 legacy lifetime only; U15 and I21/I22 pending |
| electron impact | `Moon::ElectronImpactIonizationRate` | species, fixed density/speed/cross section → s^-1 | D05 not qualified | U16/I23 SKIPPED |
| charge exchange | generic library available; no lunar configuration found | undeclared | required authoritative cross sections absent | U17/I24 SKIPPED |
| plasma/magnetic driver | typical fields or dormant coupler/datafile path | driver fields → E/B/plasma state | D06 QA ERROR | U18/I25/I26/I35/I40 SKIPPED |
| He source/reservoir | `Exosphere_Helium.cpp`; no production registration | projected alpha flux → rate/Maxwellian particle | RNG; no reservoir contract | U19/I27/I34/I35 SKIPPED |
| Ne source/accommodation | `Exosphere_Neon.cpp`; no production registration | projected flux → rate | incomplete generator path | U20/I27/I36 SKIPPED |
| radiogenic Ar | sticking function only | no geography/transient source | D10 supports later comparison but is incomplete | U21/I28/I37-I39 SKIPPED |
| Na impact/PSD/thermal/sputtering | generic exosphere called by sphere injection | configured rates/energies/population → particles/s | multiple RNG draws; built-in impact path active once, historical `MySource` disabled | U22 kernels; I02/I03 linked dispatch/rate audit and I29 budget/statistics pending |
| meteoroid forcing | no production path | undeclared | D11 required | U23/I30 SKIPPED |
| H2O/OH migration/chemistry | no production path/species | undeclared | D04/D15 required | U24/I31/I39 SKIPPED |
| LOS/brightness/column operators | `Moon.cpp`, subsolar/velocity samplers | density, velocity, ray, g factor → m^-2, mean m/s, rayleigh | mesh sampling; observer geometry | U01/U26 SKIPPED; operator gates required |
| mission operators | Kaguya tables/code; no qualified LADEE/LAMP/LACE/PACE adapters | epochs/attitude/FOV/mass-energy-angle cuts → observables | D07-D09/D12-D15 mixed readiness | U26 and I32-I40 pending |
| diagnostics/budgets | generic source, loss, return, sticking, surface samplers | statistical weights/time → rates/densities | MPI reductions and sampling window | no closed conservation gate yet |
| checkpoint/restart | PIC particle restart configured in `moon.input` | particle state → restart file | every 20 iterations | surface/RNG/sampler continuity unverified |
| MPI/OpenMP/stochastic behavior | generic PIC/exosphere plus input `random_seed` | layout and seed → statistical solution | I11 freezes seed 731993 and MPI=1/OpenMP=1; later layout sensitivity remains open | I11 controls implemented; I42-I46 pending |
| provenance/run manifests | linked Python runner in `test/integration_common.py` | revision/build/config/data hashes/layout/seed → manifest | canonical data root external | I01 and I11-I13 manifests implemented; release gates I42-I46 pending |

### M0 production mover and linked diagnostics

`input/moon.input` explicitly selects the production second-order force mover.
This is a physical wiring requirement, not a numerical preference: AMPS's
otherwise active default is named
`UniformWeight_UniformTimeStep_noForce_TraceTrajectory_SecondOrder` and does
not call the Moon acceleration callback. The selected force mover evaluates
the existing `Moon::TotalParticleAcceleration` at the midpoint and retains
the existing sphere/domain callbacks. Consequently lunar and differential
gravity, rotating-frame terms, illuminated Na radiation pressure, and Lorentz
acceleration for configured ions are applied in ordinary production dispatch.

The selected second-order mover now passes its located midpoint AMR leaf—not
the not-yet-assigned final-position node—to the midpoint acceleration callback.
It returns the existing domain-loss status before force evaluation if that
midpoint is outside the finite mesh. This is required for charged-species
callbacks, which legitimately inspect node/block data, and is covered by the
dedicated Na/Na+ U06 and I08 configuration.

When that force mover hands a sphere-crossing segment to AMPS's production
surface tracer, the tracer now resolves its final cell from the current
`startNode`. Surface re-emission can update `startNode` without assigning the
temporary block-face neighbor `newNode`; using the latter caused null/stale
cell lookups in the four-rank canonical Moon run. The final particle-list
insertion already uses `startNode->block`, so the lookup and insertion now
share one authoritative leaf.

The same surface-tracer fallback now releases particle-buffer storage whenever
the registered surface callback returns `_PARTICLE_DELETED_ON_THE_FACE_`, in
both its inside-sphere recovery and ordinary intersection branches. The
callback declares sticking/absorption; the mover owns buffer deletion, as in
AMPS's second-order and Boris movers. Without that deletion, absorbed
particles were allocated but absent from every cell list, eventually tripping
the four-rank `CheckParticleList()` conservation invariant.

The optional executable arguments
`--moon-integration-test ID --integration-output FILE` are limited to bounded
M0 diagnostics after the complete `amps_init()` path. They record raw SI
states and callback results; independent acceptance calculations remain in
`test/integration_common.py`. Normal no-argument application behavior is
unchanged. Diagnostic particles follow the compiled production data layout:
the unity individual-weight correction is initialized only when
`_INDIVIDUAL_PARTICLE_WEIGHT_MODE_` is `ON`; in the normal Moon `OFF` build no
such field exists and the block/species production weight applies directly.
See [`test/README.md`](test/README.md) for the clean-build campaign,
per-test configurations, preserved artifacts, and exact commands.

## Capability status

- **VERIFIED locally:** only the narrow kernels that pass implemented U tests:
  lunar point gravity; Earth-shadow classifications and radiation-pressure
  inverse-square invariant; current resolution callbacks; analytic temperature
  and sticking control points; legacy Na lifetime; impact normalization and PSD
  energy-density formula.  See test results rather than assuming PASS.
- **EXPERIMENTAL:** the assembled Na exosphere, input-selected LOLA surface,
  source mixture, surface
  interaction and fixed surface-density behavior, chemistry, and sampling.
  Production code exists, but required linked, convergence, conservation, and
  operator gates are incomplete.
- **DISABLED or missing:** charged species/Lorentz trajectories, terrain
  shadowing, Diviner/dynamic
  thermal state, cold trapping, qualified plasma/chemistry drivers, He/Ne/Ar
  campaigns, meteoroids, H2O/OH, and most mission operators.
- **VALIDATED:** none established by this audit.  No observational comparison
  may proceed until its predecessor gates and package readiness are satisfied.

## Verification tests

The detailed U-series contracts, per-test runners, explicit reference values,
and status handling are in [`test/README.md`](test/README.md).  The basic use is:

```sh
python3 srcMoon/test/run_tests.py --list
python3 srcMoon/test/run_tests.py --all \
  --output-dir test_output/srcMoon/unit
```

The common C++ adapter links `build/libAMPS.a` and calls compiled production
kernels.  It is intentionally thin; analytical calculations and control values
are the reference, not frozen AMPS output.  Missing capabilities return
`SKIPPED`, not PASS.

As a production-wiring smoke check, the 500 km U08 surface completed the
configured 100-step executable run with one and two MPI ranks on 2026-10-10.
This exercised CEA reload, cut-cell construction, injection, surface exchange,
facet impact/re-emission, and restart output. It is not an I14 result: the
coarse smoke configuration has no convergence claim and D01 lacks the required
README, provenance, and QA records, so the capability remains EXPERIMENTAL and
NOT VALIDATED.

The normal build is Make-based; no CMake project exists in this checkout.  When
recompilation is required, remove `build/` first:

```sh
rm -rf build
make -j
```

SPICE kernel paths and Moon/exosphere diagnostic output paths are constructed
with `std::string`, not fixed-size `sprintf` buffers. This preserves the full
configured path and removes the compiler's format-overflow diagnostics from
both SPICE-enabled and canonical no-SPICE builds.

For the canonical regression, use a clean build and require a zero-byte diff:

```sh
rm -rf build
make -j test_Moon
ls -l test_Moon.diff
```

The four-rank `test_Moon.ref_np4` baseline was deliberately reviewed and
updated on 2026-10-11 after the production input began selecting the
force-aware second-order mover. The superseded reference
(`9386ea2d353db638ce84bf1f202d55cb7276443f800938c6365807be8a8deef3`)
was generated with the former implicit no-force mover and is not retained as
an alternative runtime configuration: a no-force production mode would omit
the configured lunar and non-inertial physics. Before approval, two complete
four-rank runs reproduced the replacement output byte-for-byte. The reviewed
replacement reference has SHA-256
`abc1d9c91f61a1c44e06e98bde8739b2de1650b7a6caef783bb97081f2eb80b4`.
This is a whole-application regression baseline, not an independent physics
oracle; future changes still require an explicit review and must never update
the reference automatically.

### Debug builds and GDB

`DebuggerMode=on` in `moon.input` enables AMPS runtime checks, but it does not
add source-level debug information. For a GDB investigation, use a separate
clean build with symbols and reduced optimization:

```sh
rm -rf build
make -j DEBUGC='-g3 -Og -fno-omit-frame-pointer' test_Moon
gdb --args run_test_Moon/amps -input srcMoon/examples/lola_surface.in
```

For a noninteractive single-process backtrace:

```sh
gdb --batch -ex run -ex 'thread apply all bt full' \
  --args run_test_Moon/amps -input srcMoon/examples/lola_surface.in
```

Debug-mode output is diagnostic evidence, not the canonical optimized
regression. Recreate `build/` and rerun the normal `make -j test_Moon` before
recording a release or regression result.

Column-integral variable names are assembled in a bounded `std::string` and
copied once into the AMPS `_MAX_STRING_LENGTH_PIC_` output buffer. Exceeding
that contract is a fatal configuration error. This avoids the former
self-referential `sprintf` operation, whose overlapping source/destination was
undefined and produced `-Wformat-overflow` diagnostics in optimized builds.
The shared exosphere ground-observation initializer likewise validates each
generated UTC string against its destination before copying it; it fails
instead of truncating an epoch that would alter observation geometry.

## External validation data status

The canonical root is `/data/vtenishe/moon_validation_data`.  This audit found
D04 chemistry QA marked PASS, D06 driver QA marked ERROR, D10 incomplete, and
several observational packages partial, incomplete, or provisional. The D01
`lola/raw` directory contains `LDEM_4.IMG` and `LDEM_4.LBL`, but the package
lacks `README.md`, `provenance.json`, and `qa_report.json`; it is therefore not
qualified as READY. D02 now has immutable LOLA illumination/PSR products but
remains `PARTIAL`. D03 now contains the GCP, PRP, and both complete frozen RDR
date acquisitions, plus reproducibly processed tables. In particular, the I17
RDR transformation has a `READY` processor receipt with 6,353 complete bins,
but the parent `diviner/` package still lacks its package-level `README.md`,
`provenance.json`, and `qa_report.json`; D03 is therefore not READY and no
dataset-dependent test was scored. These are preflight facts, not validation.
A future linked test must re-read each selected package's README,
`provenance.json`, and `qa_report.json`, verify all referenced hashes, and
record the exact processed file in its run manifest.

## Clarification required before implementation proceeds

The source does not answer the following questions, so this audit does not
choose values or behavior:

1. Should surface abundance remain prescribed at `2.3e16 m^-2`, or is the
   target model a conserved reservoir?  If conserved, what residence-time,
   diffusion/migration, trapping, and delayed-release equations and parameters
   are authoritative?
2. What is the intended production species table beyond the confirmed M0
   separation (Na-only normal baseline and Na/Na+ I08 fixture), including He,
   He+, Ne, Ne+,
   Ar, Ar+, H2O, OH), including masses, charges, reactions, and particle-weight
   relationships?
3. Which publications or archived tables are authoritative for the absolute Na
   impact, PSD, thermal-desorption, radiation-pressure, and sticking parameters?
   Existing comments contain abbreviated or uncertain citations and cannot
   serve as provenance.

Until these are resolved, the corresponding tests remain explicitly skipped
and the related capabilities remain experimental or disabled.

The earlier LOLA-selection question is resolved for this implementation by the
input file plus native label: `LDEM_4`, label-defined datum/frame and
east-positive longitude, centre-registered bilinear interpolation with wrapped
longitude and clamped polar rows, icosphere resolution expressed as maximum
great-circle edge metres, CEA for AMPS, and `FETRIANGLE` for Tecplot.
