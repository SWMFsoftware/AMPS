# AMPS lunar exosphere application (`srcMoon`)

## Scope and evidence

This document describes the code that is present in this checkout, not a
desired lunar model inferred from file names or downloaded data.  The audit was
performed at source revision `17e86acbafa264a43f2f06b3d2113b29466b5a20`.
Generated files under `build/` were inspected to determine the effective
configuration, but they are build products and must never be patched as the
primary implementation.

The controlling verification contract is
[`development_and_validation.md`](development_and_validation.md).  Its SHA-256
is
`ca30ecc80e9eebc7741ea79eabde58852a7eae8e4e7d9bae5bdc6892257c76de`.
The external roadmap at
`/data/vtenishe/moon_validation_data/development_and_validation.md` has the
same hash, so no contract conflict was found during this audit.

Capability words in this document use the roadmap meanings.  In particular,
source presence, downloaded data, and preflight checks do not establish
observational validation.

## What the current executable is

The active `moon.input` requests one species, neutral Na, and builds
`srcMoon/main.cpp`, `srcMoon/main_lib.cpp`, `srcMoon/Moon.cpp`, and the generic
AMPS exosphere model.  With no application `-input` argument, `amps_init()`
retains the regression's analytic sphere.  A `moon` application-input section
can instead select a generated LOLA triangulation.  SPICE/orbit evolution is
still disabled in the generated build.

| Item | Effective setting | Evidence and consequence |
|---|---|---|
| target | `_MOON_` | `moon.input` generated definition |
| object/body/solar-orbital names | `Moon`, `IAU_MOON`, `LSO` | `Moon.cpp` |
| species | Na only (`_NA_SPEC_ == 0`) | `SpeciesList=Na`; ion/He/Ne/Ar branches are compiled out or unreachable |
| timestep | species-dependent global | `moon.input`; surface inventory exchange requires this mode |
| particle weight | species-dependent global; no individual correction | `moon.input` |
| reproducible path | on | `ForceRepeatableSimulationPath=on`; this is not by itself a recorded seed contract |
| surface | analytic sphere by default; input-selectable LDEM_4 triangulation | `amps_init()`, `MoonInput.cpp`, `LunarSurface.cpp` |
| orbit/SPICE | off | `_EXOSPHERE__ORBIT_CALCUALTION__MODE_ == _PIC_MODE_OFF_` |
| surface temperature | analytic cosine law | no Diviner or time-dependent thermal state |
| Na sticking | generated Yakshinskiy-2005 digitized table | parser rewrites the source default |
| accommodation | constant 0.2 for Na | applies only to a non-sticking collision |
| surface content | prescribed `2.3e16 m^-2` for Na | user-defined mode overwrites accumulated inventory each exchange |
| active sources | impact vaporization, thermal desorption, PSD, plus a user-defined source mapped to impact vaporization | sputtering off; duplicate impact registration is discussed below |
| chemistry | photolytic reactions on; legacy constant Na lifetime | no qualified multispecies network |
| restart | particle restart output requested every 20 iterations | no audited restoration of surface/thermal/sampler/RNG state |

The generated exosphere code includes the configured impact process twice:
once as the built-in impact-vaporization source and once as user source
`MySource`, which calls the same rate and generator.  Consequently
`totalProductionRate()` adds the same configured impact rate twice and the
source selector exposes both IDs.  This is a verified wiring fact, not an
assertion that doubling was intended.  It must be resolved before source-budget
or observational work.

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
   production macro.  In the present Na-only, orbit-off build the only active
   force is lunar point-mass gravity.
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
- `SO_FRAME` is the Moon-centred `LSO` frame whose x direction is used as the
  Sun direction by analytic surface and shadow logic.  `IAU_FRAME` is
  `IAU_MOON` and is used for body-fixed surface state.
- Six-by-six SPICE transformations carry both position and velocity.  Surface
  impacts transform SO to IAU before body-fixed binning and re-emission, then
  transform back to SO.
- The orbit-on branch advances ephemeris time by the smallest global species
  timestep.  Its use of the name `MSGR_HCI` is a retained Mercury-era dependency
  and has not been qualified for a lunar run.

## Mesh and surface geometry

### Production surface selection

`srcMoon` accepts an optional application file:

```sh
./amps -input srcMoon/examples/lola_surface.in
```

Only the unique section delimited by `#section begin: moon` and
`#section end` is interpreted. Names are case-insensitive, assignments are
`name = value`, `!` begins a comment, and duplicate, unknown, or missing
settings are errors. Relative paths are resolved against the input file. The
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
```

`spice_path` must contain `cspice/include/SpiceUsr.h`,
`cspice/lib/cspice.a`, and `Kernels/`, matching `/home/vtenishe/SPICE`. It is
recorded and checked at runtime; it does not override the build-time
`SPICE=off` configuration.

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

The orbit-on update obtains lunar, solar, and terrestrial SPICE states, derives
the rotation between consecutive `LSO` frames, and fills SO↔IAU transforms.
Because orbit mode is off and the retained inertial frame is unresolved, lunar
rotation and the orbit-dependent forces are currently **DISABLED**, not
verified by compilation of the dormant branch.

`EarthShadowCheck()` implements a cylinder of Earth radius extending
anti-sunward from Earth; it is not a conical umbra/penumbra model.

## Sources and particle injection

The configured Na source values are:

- impact vaporization: `1.69e22 s^-1`, 6000 K, zero heliocentric power index;
- thermal desorption: binding energy `1.85 eV`, vibrational frequency
  `1.0e13 s^-1`, multiplied by the surface population and the Boltzmann factor;
- photon-stimulated desorption: photon flux `2.0e18 m^-2 s^-1`, cross section
  `3.0e-25 m^2`, configured injection-speed interval 10–10000 m/s;
- solar-wind sputtering: configured off.

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

1. transforms the state from LSO to IAU_MOON;
2. computes analytic surface temperature;
3. samples return flux and impact speed;
4. draws a Bernoulli sticking decision;
5. if sticking, adds `reemission_fraction × particle_weight` to the body-fixed
   adsorption counter and deletes the particle;
6. otherwise draws a thermal speed and random outward direction, blends its
   magnitude with incident speed using accommodation coefficient 0.2, converts
   back to LSO, and retains the particle.

For Na, the generated build uses linear interpolation through the embedded
Yakshinskiy table (100–495 K, with endpoint clamping); the 100 K control point
is 0.9983.  Ar uses a piecewise base-10 law with control values tested at 88,
110, and 158 K.  He and Ne return zero sticking.  The re-emission fraction
returned by the Na and Ar kernels is 1.0.

The analytic temperature is 100 K on the nightside or in terrestrial shadow,
and `100 + 280 cos(theta)^(1/4) K` otherwise.  There is no Diviner-backed or
dynamic thermal state.

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

SPICE-dependent observer geometry is compiled out in the present orbit-off
build.  The embedded Kaguya tables and code presence are not operator
verification and are not D07 validation.  No qualified production operators
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
| volume AMR/refinement | `localResolution`, mesh init/build | metres in LSO → target cell size | deterministic | U07 callback only; I11-I13 convergence pending |
| spherical surface | `amps_init`, `cInternalSphericalData` | 60×100, lunar radius → areas/normals/elements | deterministic geometry | U07 arithmetic; intersections/area invariants pending U01/I14 |
| LOLA triangulated geometry | `MoonInput.cpp`, `LunarSurface.cpp`, `amps_init` | native IMG/LBL and requested max edge m → CEA/Tecplot mesh | raw D01 exists but package QA/provenance are absent; static orbit-off mode only | U08 kernel subcheck; campaign SKIPPED / NOT VALIDATED; I14 pending |
| intersections/normals | generic sphere/ray services called by mover/LOS | metres/direction → hit, element, radial normal | stochastic position for injection | U01 SKIPPED; no isolated production fixture |
| illumination/Earth eclipse/PSR | cosine and `EarthShadowCheck`; callers in temperature, force, brightness | LSO position → cosine/boolean | D02 absent | U05 cylinder invariant; terrain/PSR U09/U13 SKIPPED |
| simple temperature | `Exosphere::GetSurfaceTemperature` | cosine, LSO position → K | deterministic | U12 linked probe PASS-capable |
| Diviner/dynamic thermal state | no production implementation | required local time/history → K | D03 absent | U10/U11 SKIPPED |
| lunar gravity | `Moon::TotalParticleAcceleration`; mover macro | m → m s^-2, LSO | deterministic | U03 linked analytic probe; trajectory I05 pending |
| Sun/Earth differential and rotating forces | same kernel; `amps_time_step` populates state | SPICE state, m, m/s → m/s² | kernel/epoch not frozen | U04 SKIPPED; orbit mode off/frame unresolved |
| Na radiation pressure | same kernel plus `Na.h` table | heliocentric radial speed/distance → m/s² | embedded table; deterministic after state | U05 scaling/shadow only; absolute source qualification pending |
| Lorentz/ion mover | same kernel; typical or coupled E/B | q, kg, V/m, T, m/s → m/s² | no active ion; D06 ERROR | U06/U18 SKIPPED |
| surface impact classification | generic surface-interaction callback | species/state/weight → delete or re-emit | Bernoulli RNG | U12 deterministic kernels only; linked statistics I18 pending |
| sticking/accommodation/re-emission | `Moon.cpp` plus generic callback | K/probability; incident velocity → outward m/s | RNG for decision/speed/direction | U12 control points; distribution/conservation pending |
| surface inventory/release | `Exosphere_Parallel.cpp::ExchangeSurfaceAreaDensity`, source depletion | weighted particles per element → population/source rates | MPI reduction; source RNG | fixed-density overwrite prevents reservoir verification; I19 pending |
| cold trapping | no production PSR/trap model | not defined | D02/D03 absent | U13/I20 SKIPPED |
| photoionization/dissociation | `ExospherePhotoionizationLifeTime`, generic processor | species/position → lifetime s/reaction | reaction RNG in linked model; D04 QA PASS | U14 legacy lifetime only; U15 and I21/I22 pending |
| electron impact | `Moon::ElectronImpactIonizationRate` | species, fixed density/speed/cross section → s^-1 | D05 not qualified | U16/I23 SKIPPED |
| charge exchange | generic library available; no lunar configuration found | undeclared | required authoritative cross sections absent | U17/I24 SKIPPED |
| plasma/magnetic driver | typical fields or dormant coupler/datafile path | driver fields → E/B/plasma state | D06 QA ERROR | U18/I25/I26/I35/I40 SKIPPED |
| He source/reservoir | `Exosphere_Helium.cpp`; no production registration | projected alpha flux → rate/Maxwellian particle | RNG; no reservoir contract | U19/I27/I34/I35 SKIPPED |
| Ne source/accommodation | `Exosphere_Neon.cpp`; no production registration | projected flux → rate | incomplete generator path | U20/I27/I36 SKIPPED |
| radiogenic Ar | sticking function only | no geography/transient source | D10 supports later comparison but is incomplete | U21/I28/I37-I39 SKIPPED |
| Na impact/PSD/thermal/sputtering | generic exosphere called by sphere injection | configured rates/energies/population → particles/s | multiple RNG draws; duplicate impact registration | U22 kernels; I29 budget/statistics pending |
| meteoroid forcing | no production path | undeclared | D11 required | U23/I30 SKIPPED |
| H2O/OH migration/chemistry | no production path/species | undeclared | D04/D15 required | U24/I31/I39 SKIPPED |
| LOS/brightness/column operators | `Moon.cpp`, subsolar/velocity samplers | density, velocity, ray, g factor → m^-2, mean m/s, rayleigh | mesh sampling; observer geometry | U01/U26 SKIPPED; operator gates required |
| mission operators | Kaguya tables/code; no qualified LADEE/LAMP/LACE/PACE adapters | epochs/attitude/FOV/mass-energy-angle cuts → observables | D07-D09/D12-D15 mixed readiness | U26 and I32-I40 pending |
| diagnostics/budgets | generic source, loss, return, sticking, surface samplers | statistical weights/time → rates/densities | MPI reductions and sampling window | no closed conservation gate yet |
| checkpoint/restart | PIC particle restart configured in `moon.input` | particle state → restart file | every 20 iterations | surface/RNG/sampler continuity unverified |
| MPI/OpenMP/stochastic behavior | generic PIC/exosphere | layout and seed → statistical solution | random stream/layout not recorded by a Moon manifest | I42-I46 pending |
| provenance/run manifests | no Moon-specific manifest writer found | revision/build/config/data hashes/layout/seed → manifest | canonical data root external | U25/I42-I46 SKIPPED |

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
- **DISABLED or missing:** orbit/SPICE evolution, charged species/Lorentz
  trajectories, terrain shadowing, Diviner/dynamic
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

For the canonical regression, use a clean build and require a zero-byte diff:

```sh
rm -rf build
make -j test_Moon
ls -l test_Moon.diff
```

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

## External validation data status

The canonical root is `/data/vtenishe/moon_validation_data`.  This audit found
D04 chemistry QA marked PASS, D06 driver QA marked ERROR, D10 incomplete, and
several observational packages partial, incomplete, or provisional. The D01
`lola/raw` directory contains `LDEM_4.IMG` and `LDEM_4.LBL`, but the package
lacks `README.md`, `provenance.json`, and `qa_report.json`; it is therefore not
qualified as READY. D02-D03 were not present. Those observations are
only preflight facts; no dataset-dependent test was scored.  A future linked
test must re-read each selected package's README, `provenance.json`, and
`qa_report.json`, verify all referenced hashes, and record the exact processed
file in its run manifest.

## Clarification required before implementation proceeds

The source does not answer the following questions, so this audit does not
choose values or behavior:

1. Should the user-defined `MySource` impact mapping be removed, or is doubling
   the built-in impact-vaporization contribution intentional?
2. What is the authoritative lunar inertial frame and SPICE kernel/epoch set
   for the orbit-on branch, and should every `MSGR_HCI` occurrence be replaced?
3. Should surface abundance remain prescribed at `2.3e16 m^-2`, or is the
   target model a conserved reservoir?  If conserved, what residence-time,
   diffusion/migration, trapping, and delayed-release equations and parameters
   are authoritative?
4. What is the intended production species table (Na, Na+, He, He+, Ne, Ne+,
   Ar, Ar+, H2O, OH), including masses, charges, reactions, and particle-weight
   relationships?
5. Which publications or archived tables are authoritative for the absolute Na
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
