# LiouvilleCOPS fitting: quick guide

New results are written under **`<AMPS>/cops_output/`** (`../../cops_output`
from LiouvilleCOPS). For example, `--out output_fit_new` writes
`../../cops_output/output_fit_new/`; the same rule applies to C++
`OutputDirectory` and plot `--save`. Absolute output paths are honored only
outside LiouvilleCOPS; paths resolving inside it are rejected.
`export COPS_OUTPUT_ROOT=/absolute/path/outside/LiouvilleCOPS` overrides the common root.
Run/table inputs prefer this root; legacy local fallback remains supported,
but the old local scientific outputs have been removed.
Use an absolute or explicit `./` input path to force the old location.
Former local results have been removed. Rebuild the executable after source changes.
See [RUN.md](RUN.md#output-location) for complete path conventions.

External inputs default to `/data/vtenishe/rosetta_data_end_of_mission_cops`:
raw BG data are in `data_COPS_L2/mtp35_cops_bg`, and L4/image sources are in
`sources`. Set `COPS_DATA_ROOT` to relocate both, or use `COPS_L2_ROOT` and
`COPS_SOURCES_ROOT` separately.

For the all-source NG/RG comparison after output cleanup, regenerate every
derived prerequisite before invoking the collector:

```bash
python3 tools/cops_hk_table.py
python3 tools/digitize_cops_figures.py
env OMP_NUM_THREADS=4 ./liouville_cops fit_base.in
env OMP_NUM_THREADS=4 ./liouville_cops m004.in
python3 tools/collect_cops_sources.py --show
```

The first two commands are preprocessing steps. The next two create the SHAP5
fitting base and M004 altitude table under `<AMPS>/cops_output`. The M004 input
uses `GeometryOnly = on`: its exact SPICE positions and M004 nearest-surface
altitudes are usable, but its model-pressure columns are zero and it is not a
fitting base. The collector lists all missing prerequisites
without a Python traceback.

## Purpose of the fitting procedure

The purpose is to infer a surface gas-emission model that reproduces the
measured nude-gauge (NG) and ram-gauge (RG) pressures along Rosetta's trajectory.

The forward calculation predicts pressures from the nucleus shape, surface
activity, temperature, spacecraft geometry, and instrument response. Fitting
adjusts:

- H2O surface activity: up to 25 spatial coefficients.
- Surface temperature: `Tmax`, `T75`, and `Tmin`.
- CO2 activity: one multiplier on the base spatial distribution.

It minimizes logarithmic differences between modeled and observed pressures.
Optional regularization penalizes departure of the H2O activity map from the
base map.

Using both gauges provides complementary constraints: in the static model,
NG density depends approximately on emission flux divided by the square root
of temperature, whereas the main RG flux contribution depends on emission
flux. Their combination helps distinguish source strength from temperature.

The result is a fitted physical source model, not simply a smoothed pressure
curve or an automatic reconstruction of missing measurements. It can produce
predictions where measurements are absent, but those predictions remain
model-dependent. A good pressure fit does not establish a unique global
activity map or validate the collisionless approximation.

## Implemented surface boundary conditions

There are two outgoing angular emission laws, selected by `SurfaceBC`.
Let `F` be the total emitted molecular flux per unit surface area, and
`theta` the angle from the local outward surface normal.

| Setting | Emitted flux per unit surface area and solid angle | Interpretation |
|---|---|---|
| `SurfaceBC = flat` | `F cos(theta) / pi` | Lambertian emission: strongest along the normal, decreasing toward grazing directions. |
| `SurfaceBC = rough` | `F / (2 pi)` | Uniform emission per solid angle across the outward hemisphere, giving relatively more grazing emission. |

Both apply only to the outward hemisphere and integrate over solid angle
to the same total flux `F`.

`rough` does not explicitly resolve microscopic roughness. It substitutes a
different angular emission law; it does not introduce additional facets or
simulate multiple reflections.

Both use the temperature-dependent emitted speed distribution:

```text
h(v) = 2 beta^4 v^3 exp(-beta^2 v^2)
beta = sqrt(m / (2 k T))
```

Here `v` is speed relative to the emitting surface, `m` is molecular mass,
`k` is Boltzmann's constant, and `T` is surface temperature. The distribution
`h(v)` integrates to 1 for speeds from zero to infinity.

The nucleus blocks obstructed trajectories. There is no implemented
reflected/re-emitted return flux from such blocked paths. The emission laws
are implemented in [src/Liouville.cpp](src/Liouville.cpp).

### Related settings and fitting implications

Surface temperature can be prescribed using `TemperatureModel = analytical`,
`constant`, or `bjorn`; these are separate choices from the angular boundary
condition. The current fitter supports the analytical temperature law.

Both supplied observation inputs use:

```ini
SurfaceBC = flat
FlightTimeRotation = none
```

The fitting procedure does not choose or fit the boundary-condition type.
To compare `flat` and `rough`:

1. Generate separate forward bases with the respective setting and distinct
   output directories.
2. Fit each using the same observations and fitting options.
3. Compare residuals, fitted parameters, and convergence.

Changing `SurfaceBC` requires regenerating the cached kernels; editing the
input after generating a base does not change its boundary condition.

### Rebuild the base automatically for a different boundary condition

Use the new driver with the original observation input. First preview it:

```bash
python3 tools/rebuild_cops_boundary.py --input fit_base.in \
  --surface-bc rough --out output_fit_base_rough --threads 4 --dry-run
```

Then remove `--dry-run` to regenerate the base:

```bash
python3 tools/rebuild_cops_boundary.py --input fit_base.in \
  --surface-bc rough --out output_fit_base_rough --threads 4
```

`--surface-bc` accepts `flat` or `rough`; `--out` must name a directory that
does not already exist. `--input` defaults to `fit_base.in`; `--threads` is
optional, and `--executable FILE` selects a solver other than the default
application-root `liouville_cops`. Relative `--out` resolves under the output
root. Other relative source paths use the application root, with generated
input/HK/activity lookup as described in RUN.md; relative `Kernel` entries
use `KernelDirectory`.
The tool needs Python 3.8+, a compiled solver, and the configured external
data. Dry runs create nothing and do not check external-file availability.

The script preserves temperatures, activity, and data-selection settings,
forces kernel output on, and reruns the forward calculation. It requires
static transport, analytical surface temperature, and temperature correction
factor 1. It rejects incompatible inputs and does not perform fitting.

The new directory contains the recomputed model table, surface elements,
point kernels, `boundary.in`, `source.in`, `rebuild.log`, and `rebuild.json`.
Do not use partial output from a failed run. Keep the generated `boundary.in`
with the base; the fitter reads it.

After success, follow the fitting steps below, replacing `output_fit_base`
with `output_fit_base_rough` and choosing fresh selection/fit directories.
Use the same data-selection and fitting options to compare the two laws.

**Temperature and recalculation:** changing the boundary law requires new
`Kn`, `KPhi`, and `KV` angular kernels and new density/flux/pressure predictions.
Those static kernels contain no surface temperature. Surface temperature
enters through molecular speeds: NG density and the RG spacecraft-velocity
term contain `1/sqrt(T_surface)`; the main RG flux term does not. The fitter
can change analytical surface temperatures using the new kernels without
rerunning SPICE. Instrument NG/RG temperatures enter pressure conversion
separately and remain fixed in the fit. The rebuild does not automatically
adjust any temperatures: refit to obtain best-fit temperatures and source
parameters for the new boundary condition. With fixed source parameters,
total production is unchanged by switching the angular law.

Exact-rotation response kernels are temperature-dependent and are not
supported by this fitting-base driver or the current fitter. See
[README.md, section 6.6](README.md#66-change-the-surface-boundary-condition-and-rebuild-the-fitting-base)
for the detailed recalculation rules and [RUN.md](RUN.md) for all options.

## Nucleus rotation: what is included and how to control it

The current fit includes the nucleus's changing orientation throughout the
observation interval. SPICE updates the spacecraft/Sun directions and gauge
pointing at every observation time. Each surface facet then has a changing
solar zenith angle and, with shadow tracing enabled, changing sunlight
obstruction. The activity map stays attached to the nucleus.

The resulting source flux varies with illumination, with nonzero prescribed
night-side emission. Surface temperature varies through the analytical
solar-angle law. The activity coefficients and `Tmax`/`T75`/`Tmin` are shared
fit parameters across the whole interval, not separate parameters per time.
The fitter retains these temporal changes through each observation's cached
geometry, solar angle, and shadow flags; it does not rerun SPICE on every trial.

| Setting | Included effect / control |
|---|---|
| `FlightTimeRotation = none` | Current fitting mode: includes rotation between observations, but neglects rotation during each molecule's flight. It does **not** freeze the nucleus for the whole dataset. |
| `BodyFrame` and SPICE `Kernel` settings | Determine orientation as a function of time. No spin period, phase, or pole is fitted; changing the geometry requires a new base. |
| `CastShadows = on` | Includes topographic sunlight obstruction. Turning it off still leaves day/night and rotational changes; rebuild the base after changing it. |
| `SZAMode = continuous` | Smooth cosine source-rate dependence; `srccg_table` uses integer-degree source angles. This is not a rotation switch. |
| `TemperatureModel = analytical` | The supported fit temperature law; responds instantaneously to illumination. |
| `--fix-temperature` | Holds the three law parameters fixed, **not** each facet's temperature in time. To test a constant surface temperature, create an analytical base with all three parameters equal and use this option. |
| `--ordered-temperature`, `--tbounds LO HI` | Control temperature parameter constraints, not temporal lag. Ordered-mode bounds apply to Tmin, not an upper cap on Tmax/T75. |
| `FlightTimeRotation = exact` | Additional C++ forward calculation accounting for rotation during molecular flight; unsupported by the current fitter and boundary-rebuild driver. |

There is no thermal inertia, subsurface conduction, delayed heating/cooling,
or time-dependent activity map in the fit. Source flux and temperature are
separate illumination prescriptions, not a coupled sublimation/energy-balance
solution. `SpacecraftVelocity` is not a switch to disable nucleus rotation.
Changing `FitInterval` changes data selection, not the orientation model.

To assess molecular-flight rotation after a fit, copy the generated
`cops_fit.in` to a separate input, keep its fitted parameters, and edit:

```ini
FlightTimeRotation = exact
OutputDirectory = output_fit_rotation_check
```

Run that input directly with `./liouville_cops YOUR_COPY.in`, after choosing
an unused output directory. Compare with the static C++ verification run.
`SpeedNodes`, `SpeedRangeMax`, and `PathTolerance` control numerical accuracy
of this calculation. Despite its name, `exact` still assumes constant angular
velocity during each flight and uses reception-time cast-shadow flags. It
does not add thermal lag, and this forward comparison is not an exact-rotation
refit. Do not feed its version-2 kernels to the current fitter.

See [README.md, section 4.7](README.md#47-rotation-driven-source-variability-and-controls-in-the-fit)
for equations, cache dependencies, and the full control reference.

## Running the fitting procedure

Run these commands from the `LiouvilleCOPS` directory. Fitting uses cached
forward-model results, so SPICE does not need to be rerun for each optimization.
For compilation, external paths, and the complete script inventory, see
[RUN.md](RUN.md). For the physics, data limitations, and fitting equations,
see [README.md](README.md).

## 1. Check the fitter and selected data

The working fitter is included with shared-output path handling:

```bash
python3 tools/fit_cops.py --help
```

The former RCS archive has been removed; use this working file as the current
version. Python needs NumPy, SciPy, and Matplotlib.

First generate `output_fit_base/` with `./liouville_cops fit_base.in`; it will
be created under `<AMPS>/cops_output`. A usable base requires `COPS_model.dat`, `elements.bin`, and all
corresponding `kernels/point=*.bin` files, as well as the original input
recorded in the table header. Keep these consistent and preserve that input.
The supported fitting base uses:

```ini
WriteKernels = on
FlightTimeRotation = none
TemperatureModel = analytical
SurfaceTemperatureCorrectionFactor = 1
```

First export the observations that will be fitted:

```bash
python3 tools/fit_cops.py output_fit_base \
  --apply-exclude-ng \
  --exclude-ng-saturated \
  --data-only \
  --out output_selection_new
```

Inspect `output_selection_new/fit_data_used.png` and `fit_data_used.txt`.
The command also writes `fit_profiles_used.png`. Even `--data-only` loads
and checks the cached model; it is not a table-only operation.

These switches preserve the input's NG exclusion intervals and remove
saturated NG points. Without them, and without a replacement `--data` table,
the fitter normally re-enables otherwise valid NG points in those exclusion
intervals and retains saturated points. The housekeeping table is needed
even when saturation exclusion is not enabled.

## 2. Run a fit

```bash
python3 tools/fit_cops.py output_fit_base \
  --apply-exclude-ng \
  --exclude-ng-saturated \
  --degree 4 \
  --lambda 3e-3 \
  --ordered-temperature \
  --max-nfev 1000 \
  --out output_fit_new
```

Here `3e-3` is an example regularization weight, not an established optimal
value. Use a fresh output directory to avoid overwriting previous results.

### What is fitted?

By default, there are 29 free parameters:

| Quantity | Parameters | Control |
|---|---|---|
| H2O surface activity | 25 coefficients through nominal degree 4 | `--degree 0..4` controls which harmonic groups vary. Higher-degree coefficients retain their base values. |
| Surface temperature | `Tmax`, `T75`, `Tmin` | `--fix-temperature` keeps the base temperatures. |
| CO2 activity | One nonnegative multiplier, `kappa` | `--fix-kappa` holds it at 1. |

The fit minimizes logarithmic NG/RG pressure discrepancies, with equal
aggregate weight for each gauge, plus the optional H2O activity-prior
penalty. It does not independently fit a pressure at every observation.

H2O coefficients are unbounded, but the resulting activity map is clipped
at zero; the CO2 multiplier is nonnegative. By default, each temperature
parameter is independently bounded to 30–400 K through `--tbounds`.
`--ordered-temperature` enforces `Tmax >= T75 >= Tmin`; in that mode only
Tmin receives the stated upper bound, so Tmax/T75 can exceed it.

The number of free parameters is `(degree+1)^2`, plus 3 if temperature is
fitted, plus 1 if kappa is fitted. For degree 4:

| Options | Free parameters |
|---|---:|
| Neither freezing switch | 29 |
| `--fix-temperature` | 26 |
| `--fix-kappa` | 28 |
| `--fix-temperature --fix-kappa` | 25 |

Higher-degree coefficients retained by `--degree N` are not zeroed. There is
no dedicated switch to freeze all H2O coefficients: even degree 0 leaves
one free coefficient. Nominal parameter counts do not guarantee that all
parameters are independently constrained by the data.

### Which physics is held fixed?

"Fixed" means not adjusted by the optimizer, not necessarily constant in time.

| Quantity / assumption | Fixed part |
|---|---|
| Nucleus shape | Mesh, facet geometry, areas, and spatial basis. |
| Rotation and trajectory | SPICE orientation, spacecraft motion, observation times, gauge position and pointing. Spin period, pole, phase, and timing offsets are not fitted. |
| Angular emission | `SurfaceBC = flat` or `rough` chosen for the base, not selected by the optimizer. |
| Illumination | Cached solar angles/shadows and the chosen `CastShadows`/`SZAMode` prescriptions. |
| Source-rate law | Heliocentric-distance laws and night/day ratios: 0.2 for H2O and 0.1 for CO2. |
| CO2 map shape | Spatial pattern stays fixed; only its overall kappa multiplier can vary. |
| Temperature-law structure | Analytical solar-angle law, its 75-degree transition and floor; only the three temperature parameters can vary. |
| Gauge response | NG/RG acceptance, shielding, head-on weighting, offsets, pressure-conversion model, and species sensitivities. |
| Gauge temperatures | `TemperatureNG` and saved constant/HK-derived `T_RG`; distinct from fitted surface temperature. |
| Transport | Collisionless model, chosen spacecraft-velocity correction, and `FlightTimeRotation = none` for the fitting base. |

The nucleus still rotates, and illumination, shadows, source flux, facet
temperature, visibility, and HK sensor temperature can vary across the
observations. The fit adjusts one activity map and one set of temperature-law
parameters against that whole time series. In particular,
`--fix-temperature` does not freeze each facet's temperature history.

Changing the cached physical assumptions requires a compatible new forward
base and a new fit. Use `tools/rebuild_cops_boundary.py` for a boundary-law
change. There is no need to regenerate geometric kernels for the activity
or analytical temperature changes evaluated internally by the optimizer.

The fit does not estimate thermal inertia/lag, subsurface conduction,
time-dependent activity coefficients, calibration offsets, or missing
collision physics. Those effects are absent, not hidden adjustable parameters.
CO contributes no modeled pressure. Total production, facet fluxes, and
reconstructed pressures are derived outputs, not independent fitting parameters.

Finally, `--lambda`, data selection, bounds, and optimizer tolerances are
user-controlled fitting settings, not physical parameters learned by the
optimizer. Both gauges have equal aggregate misfit weight; no noise or
calibration parameters are fitted. See
[README.md, section 5.1](README.md#51-free-and-fixed-quantities) for the full
parameter-count table, constraints, and fixed-physics details.

### Parameters you can control

| Option | Default | Meaning |
|---|---|---|
| Positional `run` | Required | Compatible base directory, such as `output_fit_base`. |
| `--lambda L [L ...]` | `0` | Strength of the penalty for departing from the base H2O activity map; zero means no regularization. Multiple values produce a regularization scan. |
| `--degree N` | `4` | Highest harmonic degree allowed to vary, intended range 0–4. |
| `--fix-temperature` | off | Do not fit surface temperatures. |
| `--ordered-temperature` | off | Enforce `Tmax >= T75 >= Tmin`. |
| `--tbounds LO HI` | `30 400` | Temperature bounds in kelvin; see the ordered-mode caveat below. |
| `--fix-kappa` | off | Do not fit the CO2 multiplier; keep it at 1. |
| `--max-nfev N` | `1000` | Maximum function evaluations. |
| `--ftol`, `--xtol`, `--gtol` | `1e-10` each | Objective-change, step, and optimality stopping tolerances. |
| `--apply-exclude-ng` | off | Apply the base input's NG exclusion intervals. |
| `--exclude-ng-saturated` | off | Remove saturation-flagged NG points. |
| `--data FILE` | Base observations | Use a replacement pressure/selection table, such as `data/COPS_fit_data.txt`. |
| `--hk-table FILE` | `data/COPS_HK_20160930.txt`, output-root first | Housekeeping table for saturation checking, not replacement cached gauge temperatures. |
| `--out DIR` | `<run-name>_fit` | Destination under the output root. Existing same-named files are overwritten. |
| `--jobs N` | Minimum of number of lambdas and 12 | Concurrent fits when multiple lambda values are supplied. |
| `--data-only` | off | Export the selection without optimization, after model checks. |
| `--check-jacobian` | off | Compare analytic derivatives with finite differences at the initial parameters, then stop. Inspect the printed differences; there is no separate derivative pass/fail threshold. |
| `--check-tolerance VALUE` | `1e-6` | Allowed initial Python/C++ pressure discrepancy, not a Jacobian tolerance. |

Without temperature ordering, all three temperatures are independently
bounded. With `--ordered-temperature`, only Tmin receives the specified
bounds; Tmax and T75 can exceed the upper bound. The base temperatures must
already be ordered for this mode.

Shape, trajectory, instrument response, gauge temperatures, shadows, and
transport settings are fixed by the base run. Changing those requires
generating a compatible new base, not merely changing fitting options.
The current fitter cannot use exact-flight-time-rotation kernels.

### Use an exported pressure/selection table

To use the table produced by `tools/collect_cops_sources.py`:

```bash
python3 tools/fit_cops.py output_fit_base \
  --data data/COPS_fit_data.txt \
  --lambda 3e-3 --ordered-temperature \
  --out output_fit_table_new
```

The table's flags control selection; NG exclusions are not automatically
re-enabled in this mode. Add `--apply-exclude-ng` and/or
`--exclude-ng-saturated` to remove additional NG points explicitly.
The table must match the base observation IDs, exact UTC strings, and saved
RG temperatures. Inspect the resulting `fit_data_used.txt` rather than
assuming that the selected points match another run.

## 3. Examine and verify the result

Results appear under the directory selected by `--out`, for example
`output_fit_new/`:

| File | Contents |
|---|---|
| `fit_result.txt` | Parameters, residuals, convergence, and other fit diagnostics. |
| `ActivityFit.txt` | Activity coefficients: H2O fitted, CO2 scaled, CO row unchanged. |
| `fit_data_used.txt`, `fit_data_used.png`, `fit_profiles_used.png` | Actual observations selected for fitting. |
| `COPS_model.dat` | Fitted total NG/RG pressure predictions. |
| `cops_fit.in` | Input for full C++ verification. |

Reaching the evaluation limit is not convergence. Inspect `fit_result.txt`
for termination status, temperature bounds, residuals, and conditioning.

The Python-fitted table is not a complete forward-model recomputation:
species pressure/density/flux columns are set to -1, while fields such as
production rates are inherited from the base. Use C++ verification for
recomputed species results and production rates.

Review the generated input paths and output destination, then run:

```bash
env OMP_NUM_THREADS=4 ./liouville_cops output_fit_new/cops_fit.in

python3 tools/plot_cops_model.py \
  output_fit_base output_fit_new output_fit_new_verify \
  --ref output_fit_base --save fit_comparison.png
```

The generated input normally writes `output_fit_new_verify/`; reusing that
directory overwrites its existing products. An explicit `Threads` setting
in the input overrides `OMP_NUM_THREADS`.

The C++ rerun uses the original observation file and copied input exclusions;
it does not automatically reproduce a replacement `--data` table or Python
selection switches. `fit_data_used.txt` records the actual fit selection.

## 4. Optional regularization scan

To compare several weights using the same explicit selection:

```bash
env OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  python3 tools/fit_cops.py output_fit_base \
  --apply-exclude-ng --exclude-ng-saturated --ordered-temperature \
  --lambda 1e-4 1e-3 3e-3 1e-2 1e-1 \
  --jobs 2 --out output_lcurve_new

python3 tools/plot_lcurve.py output_lcurve_new \
  --save output_lcurve_new/lcurve.png
```

Each lambda starts an independent fit from the base parameters. Results are
written into `output_lcurve_new/lambda=<value>/`; `lcurve.txt` and selection
exports are in the parent directory. `--jobs` controls Python processes,
not C++ OpenMP threads. Use fewer jobs if memory is constrained.

The maximum-curvature point is a diagnostic, not a guarantee of convergence
or scientific validity. See [RUN.md](RUN.md) for the full workflow and
[README.md](README.md) for interpretation and limitations.
