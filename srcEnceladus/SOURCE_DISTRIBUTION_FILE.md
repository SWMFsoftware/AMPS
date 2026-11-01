# Enceladus fitted-source distribution file, schema 2

This file is the handoff between the standalone fit and the AMPS forward
model. It carries the complete `tenishev2014_hybrid_8` geometry, gas
kinematics, and either normalized source shapes or absolute source rates. The
forward model does not choose a fit strategy or silently substitute a compiled
source map. A production run must name an exported file explicitly.

## Complete example

```text
#begin source_distribution: e5_shared_shape
schema_version=2
source_case=tenishev2014_hybrid_8
frame=IAU_ENCELADUS
mode=shape
source_law=tenishev2010_point
stripe_surface_mapping=local_vertical_to_reference_sphere
stripe_reference_radius_m=252100

#begin source: I_baghdad
geometry=point_jet
latitude_deg=-81.5
west_longitude_deg=31.2
tilt_deg=9.7
azimuth_deg=228.7
bulk_speed_m_s=800.0
temperature_K=190.0
strength.H2O=0.12
strength.CO2=0.08
#end source

#begin source: alexandria
geometry=polyline_stripe
bulk_speed_m_s=350.0
temperature_K=180.0
strength.H2O=0.04
strength.CO2=0.06
#begin vertex: 0000
x_m=-79215.712740378852
y_m=-15727.617431889514
z_m=-238813.57353582291
#end vertex
#begin vertex: 0001
x_m=-74912.991729420974
y_m=-19258.989233198092
z_m=-239940.71143485271
#end vertex
#end source

...the other six jets and six stripes...
#end source_distribution
```

The numeric strengths above illustrate syntax only and are not an
observational result.

## Top-level contract

The file uses the common AMPS `#begin`/`#end` grammar described by
`src/INPUT_FILE.md`. Keys are case-sensitive; duplicate and unknown keys are
fatal.

| Key | Required value |
| --- | --- |
| `schema_version` | Integer `2`. |
| `source_case` | `tenishev2014_hybrid_8` in this implementation. Later source cases receive new explicit values and validators. |
| `frame` | `IAU_ENCELADUS`. Cartesian coordinates are metres in that body-fixed frame. |
| `mode` | `shape` or `absolute`, described below. |
| `source_law` | `tenishev2010_point`. Emission events have the production-flux measure described below; this is not the alternative `surface_flux_maxwellian` angular law. |
| `stripe_surface_mapping` | `local_vertical_to_reference_sphere`. Each published stripe direction is retained and its radius is replaced by the declared reference radius. |
| `stripe_reference_radius_m` | Finite positive spherical source radius in metres. It must agree with the compiled Enceladus internal boundary within `1e-6 m`. |

The distribution must contain exactly the canonical eight point jets and
seven stripe components:

```text
I_baghdad, II_damascus, III_damascus, IV_alexandria,
V_cairo, VI_baghdad, VII_baghdad, VIII_cairo,
alexandria, baghdad, cairo, damascus,
baghdad_extra, cairo_extra, damascus_extra
```

This strict receipt prevents a missing fitted coefficient from being treated
as a physical zero. Every source must have a `strength.SPECIES` key for every
species represented in the file. The first eight canonical IDs are required
to use `geometry=point_jet`; the seven named fissures are required to use
`geometry=polyline_stripe`. Swapping IDs between otherwise valid geometry
blocks is a startup error.

## Strength modes

In `mode=shape`, `strength.H2O`, `strength.CO2`, or a future
`strength.SPECIES` is a nonnegative dimensionless fraction. Each species'
fractions must sum to one within an absolute tolerance of `1e-12`. The total
physical rate is then mandatory in the main run file:

```text
total_source_rate.H2O_s-1=1.0e28
total_source_rate.CO2_s-1=2.0e26
```

In `mode=absolute`, each strength is the component's nonnegative physical
production in molecules/s. The forward model sums the components. A
`total_source_rate.SPECIES_s-1` key in the main run file is then an error, so
two competing normalization authorities cannot coexist.

The mode applies to the entire file. Mixing fractions and absolute rates in
one receipt is forbidden. H2O and CO2 may use different spatial shapes in the
same file.

## Point jets

`geometry=point_jet` requires:

- `latitude_deg`: signed planetocentric latitude in `[-90,90]`; the published
  south-latitude magnitudes are negative here.
- `west_longitude_deg`: west-positive longitude in `[0,360)`.
- `tilt_deg`: angle from the outward radial normal in `[0,90]`.
- `azimuth_deg`: clockwise from local north in `[0,360)`.
- positive `temperature_K` and nonnegative `bulk_speed_m_s`.

The position is the radial intersection with the configured Enceladus sphere.
The axis convention is identical to `standalone/plume_model.py`.

For an emitted direction \(\widehat{\mathbf v}\) and speed
\(t=v/\sqrt{2k_BT/m}\), the forward sampler uses the same production measure
that normalizes the standalone fit:

```text
p(t,v_hat) dt dOmega proportional to
    t^3 exp[-t^2 + 2 s c t - s^2] dt dOmega,
s = u/sqrt(2 k_B T/m),  c = v_hat dot source_axis,
v_hat dot surface_normal > 0.
```

Thus directions are weighted by the standalone radial moment `K3(c,s)`, and
the conditional speed density is proportional to the integrand above. It is
not an ordinary shifted Cartesian Maxwellian truncated at the surface. The
current analytic inverse-CDF recurrence is qualified for `s*c >= -2`; all
outward directions of the present `tenishev2014_hybrid_8` H2O/CO2 sources
satisfy that bound, and startup checks the complete hemisphere.

## Tiger stripes

`geometry=polyline_stripe` requires at least two ordered `vertex` children.
Each vertex has finite `x_m`, `y_m`, and `z_m`. Vertex names are unique labels;
zero-padded increasing ordinals are recommended because file order is the
physical order. A zero-length segment is an error.

The hybrid model uses precisely the standalone model's trapezoidal arc-length
node discretization: half of each adjacent segment length is assigned to its
endpoint, followed by normalization over the polyline. It does not replace
the fitted line by a continuous source or an equal-per-node distribution.

The source geometry file retains the published stripe nodes at radii
244.116--245.695 km. Schema 2 applies the selected vertical projection during
export: for each raw node \(\mathbf{x}\),

```text
x_surface = stripe_reference_radius_m * x / |x|.
```

On a spherical body the local vertical is radial, so this preserves the
planetocentric latitude and longitude of every fitted node and changes only
its radius to 252.100 km. This is not a constant-altitude translation. The
raw geometry file and its SHA-256 remain in the receipt comments, while the
runtime vertices are the projected coordinates. The adapter independently
requires both the declared reference radius and every serialized node radius
to agree with the compiled sphere within `1e-6 m`; a corrupted, unprojected,
inside, or outside receipt is rejected before MPI.

## Current E3 and E5 equal-rate startup receipts

The two flyby receipts contain both neutral H2O and CO2. They use
`mode=absolute`, and every species total is divided equally among all 15
canonical components only to provide deterministic startup inputs:

| flyby | H2O total (s^-1) | H2O per component (s^-1) | apparent-CO2 total (s^-1) | apparent-CO2 per component (s^-1) |
|---|---:|---:|---:|---:|
| E3 | `2.060909950673409e28` | `1.3739399671156061e27` | `5.356e25` | `3.570666666666667e24` |
| E5 | `4.520015885781454e28` | `3.013343923854303e27` | `1.013e26` | `6.753333333333333e24` |

The H2O authorities are the checked-in independent hybrid-8 fits
`standalone/results/e3_fit.json` and `e5_fit.json`. The CO2 totals are the
primary all-five shared-shape hybrid-8 results in `standalone/CO2.md`, section
5.1. For E3 and E5 the mass-44 observable is explicitly **apparent CO2**:
impact fragmentation can enhance it by a factor of 2--5. Therefore the
corresponding physical-CO2 brackets are `1.0712e25--2.678e25 s^-1` for E3 and
`2.026e25--5.065e25 s^-1` for E5. The generated receipt deliberately retains
the reported apparent total and records the bracket in its manifest; it must
not be presented as a uniquely determined physical CO2 production rate.

The manifests are
`standalone/results/e3_equal_h2o_co2_startup_manifest.json` and
`e5_equal_h2o_co2_startup_manifest.json`; the generated receipts are
`input/enceladus_e3_equal_h2o_co2_sources.in` and
`enceladus_e5_equal_h2o_co2_sources.in`. The default run deck selects E3; E5
is selected by changing only `source_distribution_file`. These equal
partitions are initialization fixtures, not fitted component estimates, and
must be replaced by the later standalone-fit handoff.

## Standalone exporter contract

`standalone/export_forward_sources.py` is the reference writer. The fitting
workflow supplies a small JSON manifest:

```json
{
  "schema_version": 1,
  "source_case": "tenishev2014_hybrid_8",
  "distribution_name": "E5_shared_shape",
  "mode": "shape",
  "species": {
    "H2O": {
      "jet:I_baghdad": 0.12,
      "stripe:alexandria": 0.04
    },
    "CO2": {
      "jet:I_baghdad": 0.08,
      "stripe:alexandria": 0.06
    }
  },
  "fit_provenance": {
    "study_result": "/absolute/path/to/study_results.json",
    "strategy": "shared_shape_inms_only",
    "target": "E5"
  }
}
```

All 15 labels are mandatory for each species, even when their value is zero.
Distribution names must match `[A-Za-z0-9][A-Za-z0-9_.-]*`, and species names
must match `[A-Za-z][A-Za-z0-9_]*`, so generated block names and
`strength.SPECIES` keys cannot alter the AMPS grammar.
The exporter reads the canonical geometry and kinematics from
`standalone/config/tenishev2014_stage0.json` and the stripe nodes from its
`stripe_geometry_file`; it validates normalization/rates, then projects the
raw stripe directions to that configuration's `reference_radius_m` before
creating the AMPS file. `fit_provenance` is emitted as comments so it remains
visible but does not create an unowned runtime key.

Example command:

```bash
python3 standalone/export_forward_sources.py \
  --manifest fitted_sources.json \
  --output enceladus_sources.in
```

The fitting code must convert its internal labels to the exact `jet:ID` and
`stripe:ID` names above. It must not export nuisance-background coefficients,
renormalize absolute mode, discard zero components, or select a strategy on
behalf of the forward run. The manifest, geometry configuration path and
SHA-256 hashes are recorded in comments for provenance.

## Failure behavior and extension rule

Nonfinite/negative values, incomplete species maps, shape sums outside the
tolerance, wrong geometry, unknown keys, or an unsupported source case stop
startup before MPI/PIC initialization. To add another spatial distribution,
add a new `source_case`, its geometry validator and independent reference
tests. Do not reinterpret `tenishev2014_hybrid_8` or make its required sources
optional.
