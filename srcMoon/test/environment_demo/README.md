# Illumination, terrain-shadow, and thermal-state demonstration

This runner exports presentation-ready **model output** from the production
LOLA terrain, horizon-ray, finite-disc Earth-eclipse, and depth-resolved thermal
kernels. It is deliberately outside the frozen U01--U26/I01--I46 registry and
does not assign a new test ID.

## Physics and implementation contract

The runner first requires U08, U09, and U11 to return `PASS`. A thin C++ adapter
then uses the production Moon input parser, native D01 decoder, icosphere and
radial topography adjustment, terrain BVH, SPICE `MOON_ME_DE421` Sun/Earth
geometry, uniform-solar-disc/opaque-Earth overlap, directional Bond albedo,
and `Moon::Thermal::Column` solver. The exported illumination calls the same
factored production routine used by `AdvanceProductionState`; neither Python
nor the adapter implements a second shadow or heating law.

The logical thermal/inventory surface is the same 60-by-100 equal-solid-angle
grid used by normal `srcMoon`. Each radial logical cell is intersected with the
LOLA triangulation, so exported positions and normals follow the adjusted
terrain rather than the reference sphere. `surface_environment.csv` keeps the
geometric nightside, topographic occlusion, and Earth occultation as separate
columns. `thermal_profiles.csv` contains all 80 configured finite-volume
layers for all 6,000 cells; no attractive sites are selected and no values are
smoothed, interpolated, clipped, or gap-filled.

The default initialization is intentionally literal: all columns start at the
input `thermal_initial_temperature_k` at `simulation_start_utc`, and the
default zero-second advance exports that declared initial condition. A
positive `--advance-seconds` evolves the exact production solver but remains
an initial-condition transient, not a periodic lunar thermal solution. Use
`--initialization periodic-spinup` only when the input's complete synodic cycle
and convergence tolerances are intended; it can be computationally long and
has no arbitrary cycle cap. In periodic mode the requested snapshot UTC must
not precede the converged cycle end.

The thermal model assumptions are those recorded in
`srcMoon/examples/lola_surface.in`: Hayne-type material parameters, grey
emissivity, direct solar heating, terrain-centre ray shadow, uniform apparent
solar radiance, opaque solid Earth, no solar limb darkening, and no refracted,
scattered, or transmitted terrestrial-atmosphere light. These omissions are
not fitted. A model image is not a D02/D03 observation.

## Running

After the normal clean SPICE Moon build, create an instantaneous forcing map
and the explicit 250 K initial thermal state:

```sh
python3 -B srcMoon/test/environment_demo/run.py
```

Advance the production columns for one explicit terrestrial day from the
input initial condition:

```sh
python3 -B srcMoon/test/environment_demo/run.py \
  --advance-seconds 86400 \
  --output-dir test_output/srcMoon/environment-demo/one-day-transient
```

Run the potentially long observation-independent periodic initialization and
sample at or after its configured cycle end:

```sh
python3 -B srcMoon/test/environment_demo/run.py \
  --initialization periodic-spinup \
  --snapshot-utc 2011-06-15T16:50:00.000Z \
  --output-dir test_output/srcMoon/environment-demo/periodic
```

The runner refuses to overwrite a non-empty output directory and disables
surface-inventory checkpoints in its preserved `application.in`.

## Artifacts and status

Outputs include `surface_environment.csv`, `thermal_profiles.csv`, CEA and
Tecplot LOLA surfaces, illumination/shadow and thermal PNG/EPS figures,
`application.in`, logs, visualization receipt, `run_manifest.json`, and
`result.json`. The manifest records source/build identities, dirty paths,
compiler, commands, both roadmap hashes, and D01 package/raw hashes.

`execution_status=PASS` means only that the demonstration artifacts were
generated and passed schema/range integrity checks. The result separately
records I15, I17, and I20 as `execution=NOT_RUN`, `status=SKIPPED`, and
`validation_status=NOT_VALIDATED`. It does not compare D02 illumination or D03
Diviner observations, run same-campaign I13, establish thermal time-step or
spatial convergence, or validate PSR/cold-trapping behavior.

For custom static or rotatable views, use
`srcMoon/tools/visualize_environment.py`; its complete CLI is available with
`--help`.
