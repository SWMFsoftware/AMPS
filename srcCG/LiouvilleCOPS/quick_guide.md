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
`sources`. The standard SHAP5 mesh is
`/data/vtenishe/rosetta_data_end_of_mission_cops/SHAP5_Volume.bdf`, selected
by `ShapeModel` in `cops.in` and `fit_base.in`. Set `COPS_DATA_ROOT` to
relocate the Python raw/L4/image inputs, or use `COPS_L2_ROOT` and
`COPS_SOURCES_ROOT` separately; edit `ShapeModel` if the mesh moves.

New C++ runs and fits require new output directories. C++ writes
`input.snapshot.in` plus `run_manifest.json`; the fitter writes
`fit_manifest.json` after validating every table, element, and static kernel.
The manifests keep process completion, optimizer convergence, cache validity,
C++ fitted-pressure agreement, and scientific review as separate states.
Passing validation never sets `scientific_review`.

```bash
# Read-only fitting-base validation
python3 tools/validate_cops_run.py output_fit_base

# After ./liouville_cops output_fit_new/cops_fit.in
python3 tools/validate_cops_run.py \
  --fit output_fit_new --verification output_fit_new_verify
```

Pre-W01 bases without `run_manifest.json` are refused by the fitter unless
`--allow-legacy-base` is supplied explicitly. Regenerating a manifested base
is preferred because missing input identities cannot be reconstructed safely.

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

### Named data profiles

Use the profile runner for repeatable comparisons of the three requested data
masks:

| Profile | NG/RG points | Meaning |
|---|---:|---|
| `legacy_base_mask` | 839 / 1,146 | Base `ExcludeNG`/`ExcludeRG` intervals applied. |
| `all_available` | 930 / 1,146 | All positive outside-shape data, including saturated and late NG. |
| `transition_cut_100m` | 869 / 1,146 | All RG; NG excludes `10:24:42`--`10:37:10` inclusive and resumes at `10:37:14` (99.48 m in the audited base). |

Inspect all selections without fitting:

```bash
python3 tools/run_fit_profile.py \
  --run output_fit_full_default_base --all-profiles \
  --data-only --jobs 14 --out-prefix a1_selection
```

Run all literal-default fits:

```bash
python3 tools/run_fit_profile.py \
  --run output_fit_full_default_base --all-profiles \
  --degree 4 --lambda 0 --jobs 14 --threads-per-process 1 \
  --tbounds 30 400 --max-nfev 1000 \
  --ftol 1e-10 --xtol 1e-10 --gtol 1e-10 \
  --check-tolerance 1e-6 --out-prefix a1_fit
```

Profiles are run sequentially. `--jobs 14` limits the fitter's pool of
independent lambda fits; the single default `lambda=0` has only one task, so it
cannot occupy 14 workers. Do not invent 14 lambda values merely to occupy all
cores. Each output stores `selection_profile.json`, selection/reason hashes in
`fit_manifest.json`, and reason codes in `fit_data_used.txt`. A count mismatch
stops before optimization.

### Exact end-to-end commands used for the full literal-default fit

The following is the fully explicit equivalent of the command sequence used
for the full 10-hour observation fit named `output_fit_full_default`. Options
that took their CLI defaults during the executed run are spelled out here so
the numerical setup is visible. Run it from the
`srcCG/LiouvilleCOPS` directory. It starts by regenerating a manifested static
kernel base, performs the Python optimization, reruns the fitted model with
the independent C++ forward solver, registers the C++/Python agreement, and
writes comparison figures versus time and altitude.

The commands assume that `liouville_cops` has already been compiled and that
the external shape, SPICE, observation, and housekeeping inputs configured in
`fit_base.in` exist. This definition keeps the output path portable relative
to the current AMPS checkout:

```bash
export COPS_OUTPUT_ROOT="$(cd ../.. && pwd)/cops_output"
```

Each `--out` destination must be new. The names below are the names of the
completed run and therefore already exist on the machine where that run was
made; use consistently renamed base, fit, and verification directories when
repeating it.

1. Regenerate the full SHAP5 static fitting base. This runs the C++ forward
   model for all observation points, uses the flat/Lambertian surface boundary
   condition, writes the per-point kernels required by the fitter, and records
   a run manifest. `--threads 16` is the C++ OpenMP thread count.

```bash
python3 tools/rebuild_cops_boundary.py \
  --input fit_base.in \
  --surface-bc flat \
  --out output_fit_full_default_base \
  --threads 16 \
  --executable ./liouville_cops
```

2. Fit the 25 degree-0-through-4 H2O activity coefficients, the three
   analytical surface temperatures, and the CO2 scale. The numerical defaults
   are written explicitly here. Limiting BLAS and OpenMP to one thread avoids
   hidden nested threading; this single-lambda optimizer itself uses one CPU
   core.

```bash
env OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  python3 tools/fit_cops.py output_fit_full_default_base \
  --degree 4 \
  --lambda 0 \
  --jobs 1 \
  --tbounds 30 400 \
  --max-nfev 1000 \
  --ftol 1e-10 \
  --xtol 1e-10 \
  --gtol 1e-10 \
  --check-tolerance 1e-6 \
  --hk-table "$COPS_OUTPUT_ROOT/data/COPS_HK_20160930.txt" \
  --out output_fit_full_default
```

The omissions in this command are deliberate and define the requested
literal-default selection and constraints:

- no `--ordered-temperature`: `Tmax`, `T75`, and `Tmin` are independently
  bounded, with no enforced ordering;
- no `--exclude-ng-saturated`: all 13 saturation-flagged NG samples are used;
- no `--apply-exclude-ng`: otherwise valid late NG data in the base input's
  `ExcludeNG` interval are re-enabled (91 samples for this base);
- no `--data`: the fit uses `refNG` and `refRG` from the base run;
- no `--fix-temperature` or `--fix-kappa`: temperature and the CO2 multiplier
  remain free;
- no `--allow-legacy-base`: the regenerated manifested base is required.

For this run, step 2 took about 2 hours 4 minutes. It reached
`--max-nfev 1000` rather than satisfying the strict convergence tolerances, so
its manifest correctly says `optimizer=not_converged`; completing the command
is not the same as optimizer convergence.

3. Run the generated fitted input through the C++ forward model. This produces
   the complete source/species diagnostics that the cached Python prediction
   does not calculate. The generated input contains `Threads = 16` and writes
   `output_fit_full_default_verify`; an input-file `Threads` value takes
   precedence over `OMP_NUM_THREADS`.

```bash
env OMP_NUM_THREADS=16 ./liouville_cops \
  "$COPS_OUTPUT_ROOT/output_fit_full_default/cops_fit.in"
```

4. Compare the fitted Python aggregate pressures with the independent C++
   pressures and register the result in `fit_manifest.json`. The explicit
   tolerance allows at most a `1e-6` relative NG or RG difference; this check
   does not mark the result scientifically reviewed.

```bash
python3 tools/validate_cops_run.py \
  --fit output_fit_full_default \
  --verification output_fit_full_default_verify \
  --tolerance 1e-6
```

5. Plot the base, cached fitted prediction, independently verified C++ run,
   and reference NG/RG measurements versus UTC. The figure is written inside
   the fit result directory under the shared output root.

```bash
python3 tools/plot_cops_model.py \
  output_fit_full_default_base \
  output_fit_full_default \
  output_fit_full_default_verify \
  --ref output_fit_full_default \
  --x time \
  --theme light \
  --save output_fit_full_default/comparison_time.png
```

6. Make the same comparison versus shape-model altitude.

```bash
python3 tools/plot_cops_model.py \
  output_fit_full_default_base \
  output_fit_full_default \
  output_fit_full_default_verify \
  --ref output_fit_full_default \
  --x altitude \
  --theme light \
  --save output_fit_full_default/comparison_altitude.png
```

The principal numerical report is
`$COPS_OUTPUT_ROOT/output_fit_full_default/fit_result.txt`; provenance and
authoritative process/optimizer/verification states are in
`fit_manifest.json`. The two comparison figures are `comparison_time.png` and
`comparison_altitude.png` in that directory. The complete C++ table is
`$COPS_OUTPUT_ROOT/output_fit_full_default_verify/COPS_model.dat`.

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

### Which pressure curves actually enter the fit?

The default fit uses two discrete observation series, not all curves shown by
the comparison plot:

| Fitted observation | Saved column | Immediate source used by Python |
|---|---|---|
| COPS nude-gauge pressure | `refNG` | `output_fit_base/COPS_model.dat` |
| COPS ram-gauge pressure | `refRG` | `output_fit_base/COPS_model.dat` |

The C++ base run copies these reference values from the file named by
`ObservationFile` in `fit_base.in`, currently
`/home/vtenishe/ROS-END-OF-MISSION/AMPS/srcCG/RosinaMeasurements.h`, and
converts them from mbar to Pa because `ReferenceDataUnit = mbar`. The inspected
SHAP5 reference series has 1,150 timestamps from
`2016-09-30T0:00:04` through `2016-09-30T10:39:26`. Nonpositive values are
missing-data markers and do not enter the objective.

For each selected sample, the residual is
`log(p_model) - log(p_reference)`. NG residuals are divided by the square root
of the selected NG count and RG residuals by the square root of the selected
RG count, so the two gauges have equal aggregate weight. The inspected
dataset gives these selection counts:

| Selection | NG samples | RG samples |
|---|---:|---:|
| Python default, no replacement `--data` table | 930 | 1,146 |
| Default plus `--apply-exclude-ng` | 839 | 1,146 |
| Default plus `--exclude-ng-saturated` | 917 | 1,146 |
| Both exclusion switches, as in the inspection command above | 826 | 1,146 |

These counts describe the inspected data, shape, and intervals; they are not
constants for another base. The 13 saturation-flagged NG samples from
`10:23:50` through `10:24:38` are included unless
`--exclude-ng-saturated` is supplied. Without `--apply-exclude-ng`, the fitter
also re-enables otherwise valid NG samples in the base input's `ExcludeNG`
interval.

The L4 density-derived series, digitized ESA-blog curve, digitized EGU NG/RG
curves, M004 comparison, and other plotted simulations are comparison sources;
they are not silently merged into the fit. The collector's default
`COPS_fit_data.txt` also takes its pressures from the base `refNG/refRG`
columns. Only an explicit `--data FILE` replaces the pressures and selection
flags, and that table must match the base point IDs, UTC strings, and saved RG
temperatures.

The saved reference series is closely related to the locally processed COPS
L2 housekeeping pressures, but the audit documented in `README.md` found small
numerical differences, two unmatched timestamps, and one NG missing-value
difference. Its conversion history is not established here, so do not claim
that it is byte-for-byte identical to the raw or processed L2 archive.

`fit_data_used.txt`, `fit_data_used.png`, and `fit_profiles_used.png` produced
by each run are the authoritative record and visualization of the curves and
masks that actually entered that fit. In the text table, `-1` means that the
corresponding gauge sample was not selected.

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
value. A fresh output directory is required; the fitter refuses an existing path.

### What is fitted?

By default, there are 29 free parameters:

| Quantity | Parameters | Control |
|---|---|---|
| H2O surface activity | 25 coefficients through nominal degree 4 | `--degree 0..4` controls which harmonic groups vary. Higher-degree coefficients retain their base values. |
| Surface temperature | `Tmax`, `T75`, `Tmin` | `--fix-temperature` keeps the base temperatures. |
| CO2 activity | One nonnegative multiplier, `kappa` | `--fix-kappa` holds it at 1. |

### Optional higher-degree H2O residual fit and C++ verification

Stages A2--A5 provide matching Python/C++ real, orthonormal harmonics, an
optional H2O forward-model residual, and residual fitting at degrees 5--8. The definition,
normalization, phase, coordinate frame, and coefficient order are fixed in
`HARMONIC_BASIS.md`. Validate them with:

```bash
python3 -m unittest tests.test_a2_activity_basis -v
python3 -m unittest tests.test_a3_activity_enhancement -v
python3 -m unittest tests.test_a4_enhancement_fit -v
python3 -m unittest tests.test_a5_reporting -v
python3 tools/benchmark_activity_basis.py \
  output_fit_full_default_base/elements.bin --repeat 5
```

The C++ source model can apply a compatible coefficient file. In a copy of an
observation input, choose a fresh output and set:

```ini
ActivityEnhancementModel = real_spherical_harmonics_v1
ActivityEnhancementDegree = 5
ActivityEnhancementFile = examples/activity_enhancement_zero_l5.txt
ActivityEnhancementScale = 1.0
```

Run it normally, for example
`OMP_NUM_THREADS=14 ./liouville_cops my_enhanced_run.in`. The example is an
all-zero regression file for the embedded SHAP5 baseline. Enabled runs write
`ActivityEnhancementFacets.dat`; `run_manifest.json` records the basis, frame,
degree, scale, clipping rule, legacy identity, and coefficient-file hash. The
default remains `ActivityEnhancementModel = none` with degree 4.

For a controlled residual fit, all legacy activity, temperature, and CO2
parameters must remain fixed. The complete one-fit command is:

```bash
python3 tools/run_fit_profile.py \
  --run output_fit_full_default_base \
  --profile all_available \
  --enhancement-degree 8 \
  --fix-legacy-activity --fix-temperature --fix-kappa \
  --lambda 1e-4 \
  --jobs 14 --threads-per-process 1 \
  --max-nfev 1000 --ftol 1e-10 --xtol 1e-10 --gtol 1e-10 \
  --check-tolerance 1e-6 \
  --out fit_enh_l8_all_available
```

The `b_lm` coefficients start at zero and are unbounded. `A_scale` is derived,
not fitted: it is the area-weighted RMS of the unclipped legacy linear map.
The `--lambda` penalty is the area-weighted squared unclipped enhancement,
relative to the legacy linear-map power, and pulls the residual to zero. Use
several lambda values to run independent fits in parallel; one lambda provides
only one optimizer task even with `--jobs 14`.

Run the independent C++ reconstruction and compare it with the Python result:

```bash
./liouville_cops ../../cops_output/fit_enh_l8_all_available/cops_fit.in
python3 tools/validate_cops_run.py \
  --fit fit_enh_l8_all_available \
  --verification fit_enh_l8_all_available_verify
python3 tools/plot_cops_model.py \
  fit_enh_l8_all_available fit_enh_l8_all_available_verify \
  --ref fit_enh_l8_all_available
```

Relative product paths resolve through `<AMPS>/cops_output`; replace the first
argument with the printed absolute path if required. Enhanced outputs add
`ActivityEnhancementFit.txt` and `ActivityEnhancementFacets.dat`, while
`ActivityFit.txt` remains the legacy three-by-25 format. Degree 4 with the
three fix switches is the zero-parameter control. Degrees 5--8 contain 11,
24, 39, and 56 free coefficients.

### Stage-A5 result package

Every completed fit or lambda subdirectory automatically contains:

- `fit_report.json`, whose `ReportIdentity` binds the profile/selection,
  legacy and fitted activity, basis/degree, regularizer, and lambda;
- `activity_facets.dat`, `activity_degree_power.dat`, and
  `activity_region_diagnostics.dat` for facet fields, power by degree, and
  all/seen/unseen area-weighted diagnostics;
- `activity_maps.png` for separate baseline, residual, and clipped total;
- `pressure_residuals.dat`, `pressure_comparison_time.png`,
  `pressure_comparison_altitude.png`, and `pressure_residuals.png`.

Selected observations are filled, excluded observations hollow, saturated NG
crossed, and the instrument-transition interval shaded. The log-residual table
uses `nan` for observations that did not enter the objective. “Seen” means a
facet contributes to at least one selected-data kernel; it does not mean the
facet is uniquely identifiable. Check the data/augmented Jacobian ranks and
singular values in `fit_result.txt` and `fit_manifest.json` before interpreting
individual modes.

For an additional multi-run comparison, use `plot_cops_model.py` as in the
command above (add `--save FILE` for a PNG). The reader verifies the fitted
table and report hashes plus profile, selection, report, and basis identities;
it rejects incompatible schema-v2 fits and visibly labels unmanifested legacy
tables. Its labels include profile, degree, lambda, and short report identity.

The fitter takes its initial physical values from the selected base run; it
does not keep a separate, universal starting vector. With the supplied
`fit_base.in`, which uses `ActivityModel = shap5` and does not override the
analytical-temperature defaults, the starting point and bounds are:

| Free quantity | Default initial value | Default optimizer bounds | Where/how to change it |
|---|---|---|---|
| `a00` through `a24` | The 25 H2O `Activity_SHAP5` coefficients in `src/SrcCGTables.h` | Unbounded | Select fewer free groups with `--degree N`. To change the initial map, create a three-row, 25-coefficient activity file, set `ActivityModel = file`, `ActivityFile = FILE`, and `ActivityFileHeliocentricModel = shap5` or `shap5_1` in a new base input, then regenerate the base. |
| `Tmax` | 182.09922 K | 30–400 K | Change the base input's `AnalyticalTmax`; change optimizer bounds with `--tbounds LO HI`; freeze it together with the other temperatures using `--fix-temperature`. |
| `T75` | 136.11149 K | 30–400 K | Change `AnalyticalT75` in the base input; use the same temperature controls. |
| `Tmin` | 133.0 K | 30–400 K | Change `AnalyticalTmin` in the base input; use the same temperature controls. |
| `kappa` | 1.0 relative to the base CO2 map | 0 to infinity | Freeze it with `--fix-kappa`. There is no command-line option for a different initial `kappa`; change the CO2 row of an `ActivityModel = file` base if a different effective starting CO2 source is required. |

The coefficient variables are internally scaled by the largest absolute base
H2O coefficient for numerical conditioning; `ActivityFit.txt` and reported
fit values are converted back to physical coefficient values. The scaling is
not an additional fitted parameter.

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
| `--enhancement-degree N` | `4` | Residual maximum degree, 4 for the zero-length control or 5–8 for degrees 5 through N. |
| `--fix-legacy-activity` | off | Keep all 25 legacy H2O coefficients fixed. Required, with fixed temperature and `kappa`, for degrees 5–8. |
| `--enhancement-scale VALUE` | derived | Optional assertion of the area-weighted legacy-linear RMS; it cannot override that definition. |
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
| `--out DIR` | `<run-name>_fit` | New destination under the output root; any existing path is refused. |
| `--jobs N` | Minimum of number of lambdas and 12 | Concurrent fits when multiple lambda values are supplied. |
| `--data-only` | off | Export the selection without optimization, after model checks. |
| `--allow-legacy-base` | off | Explicitly accept a base without `run_manifest.json` after strict cache validation; prefer regenerating the base. |
| `--check-jacobian` | off | Compare analytic derivatives with finite differences at the initial parameters, then stop. Inspect the printed differences; there is no separate derivative pass/fail threshold. |
| `--check-tolerance VALUE` | `1e-6` | Allowed initial Python/C++ pressure discrepancy, not a Jacobian tolerance. |

These are Python fitter options and must appear after the base-run argument,
for example:

```bash
python3 tools/fit_cops.py output_fit_base \
  --degree 2 \
  --lambda 3e-3 \
  --ordered-temperature \
  --tbounds 80 300 \
  --fix-kappa \
  --apply-exclude-ng \
  --exclude-ng-saturated \
  --out output_fit_degree2_new
```

Options alter that fit only. Initial activity and temperature values are
properties of the C++ base and must be changed before regenerating it. Make a
new input (and a new output name) rather than overwriting the existing base:

```text
OutputDirectory = output_fit_custom_base
ActivityModel = file
ActivityFile = /path/to/activity_3x25.txt
ActivityFileHeliocentricModel = shap5
AnalyticalTmax = 190.0
AnalyticalT75 = 145.0
AnalyticalTmin = 130.0
```

Then generate and fit the new base:

```bash
env OMP_NUM_THREADS=4 ./liouville_cops fit_custom.in
python3 tools/fit_cops.py output_fit_custom_base \
  --ordered-temperature --lambda 3e-3 \
  --out output_fit_custom_new
```

The numerical values in this example are illustrative configuration syntax,
not recommended physical values. Choose them only from an identified physical
or calibration basis.

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
| `fit_manifest.json` | Input/cache/selection hashes and separate process, optimizer, C++-verification, and scientific-review states. |
| `fit_report.json` | Self-identifying index for this exact profile/activity/basis/lambda result. |
| `activity_facets.dat`, `activity_degree_power.dat`, `activity_region_diagnostics.dat`, `activity_maps.png` | Facet, degree, region, and map diagnostics for baseline/residual/total activity. |
| `pressure_residuals.dat`, `pressure_comparison_time.png`, `pressure_comparison_altitude.png`, `pressure_residuals.png` | Selected/excluded/saturation/transition-aware pressure diagnostics. |

Reaching the evaluation limit is not convergence. Inspect `fit_result.txt`
for termination status, temperature bounds, residuals, and conditioning.

The Python-fitted table is not a complete forward-model recomputation:
aggregate NG/RG pressures are recomputed, while production, species
pressure/density/flux, and below-deck fields are explicitly `nan`. Use C++
verification for complete physical diagnostics.

Review the generated input paths and output destination, then run:

```bash
env OMP_NUM_THREADS=4 ./liouville_cops output_fit_new/cops_fit.in

python3 tools/validate_cops_run.py \
  --fit output_fit_new --verification output_fit_new_verify

python3 tools/plot_cops_model.py \
  output_fit_base output_fit_new output_fit_new_verify \
  --ref output_fit_new --save fit_comparison.png
```

Use the Python fit as `--ref` when the figure must distinguish observations
used by that fit. The fit prediction table carries the fitter's final
`useNG`/`useRG` masks; the base and C++ verification tables retain the base
input's exclusion flags and can therefore make re-enabled fit samples appear
hollow even though they entered the objective.

The generated input normally writes `output_fit_new_verify/`; the solver
refuses to reuse an existing directory. An explicit `Threads` setting
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
