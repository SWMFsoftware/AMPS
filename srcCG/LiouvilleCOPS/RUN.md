# LiouvilleCOPS: build and run guide

This guide covers only this standalone application and its subdirectories.
Sections 1–2 describe the multi-source pressure plot; sections 3 onward cover
compilation, paths, the forward model, fitting, verification, and every Python
entry point. See [README.md](README.md) for the model equations, data audit,
file formats, and scientific limitations.

## Output location

Run commands from LiouvilleCOPS. New scientific products default to
**`../../cops_output/`** (`<AMPS>/cops_output`), not this directory. Relative
output arguments and C++ `OutputDirectory` are placed under that root.
Absolute output paths override it only when outside LiouvilleCOPS; output paths
that resolve inside the application are rejected. Parent directories are created as needed.
Examples below use output-root-relative names:

```bash
./liouville_cops fit_base.in
# Results: ../../cops_output/output_fit_base/
python3 tools/fit_cops.py output_fit_base --out output_fit_new --lambda 3e-3
# Results: ../../cops_output/output_fit_new/
python3 tools/plot_cops_model.py output_fit_new --save plots/fit.png
# Image: ../../cops_output/plots/fit.png
```

To change the root for all tools and the rebuilt executable:

```bash
export COPS_OUTPUT_ROOT=/absolute/path/to/cops_output
```

Generated-product inputs prefer that root, with read-only fallback to existing
working/application-local products. An absolute input path or explicit `./` or
`../` selects that location directly. The former local results have been removed. Raw
data and source files are read from the external `COPS_DATA_ROOT` described
below. Configuration files stay local. Build files and the executable remain in LiouvilleCOPS. Rebuild with
`make -j4` after this source update; keep the executable in the application
directory. Shell `cp` and `>` redirection require actual filesystem paths,
not output-root-relative shorthand.

## COPS source comparison

`tools/collect_cops_sources.py` combines the available ROSINA/COPS nude-gauge
(NG) and ram-gauge (RG) measurements, converts them to common units, plots the
sources together, and marks the data selected for the LiouvilleCOPS fit.

`tools/collect_cops_sources.py` is included as a normal working file.

Run all commands below from the `LiouvilleCOPS` directory.

## 1. Build the COPS housekeeping table

```bash
python3 tools/cops_hk_table.py \
  /data/vtenishe/rosetta_data_end_of_mission_cops/data_COPS_L2/mtp35_cops_bg \
  data/COPS_HK_20160930.txt
```

`cops_hk_table.py` has two optional positional parameters:

| Parameter | Meaning |
|---|---|
| `/data/vtenishe/rosetta_data_end_of_mission_cops/data_COPS_L2/mtp35_cops_bg` | Directory containing raw `bg_*.tab` files; this is the default. |
| `data/COPS_HK_20160930.txt` | Optional output filename (also the default), under the output root. This is the table read by `collect_cops_sources.py`. |

## 2. Collect and plot all sources

The collector does not calculate spacecraft geometry itself. After a clean
checkout or output cleanup, generate its four prerequisite product groups:

```bash
# Archive housekeeping table used for RG temperature and L2 comparisons
python3 tools/cops_hk_table.py

# Five CSV curves digitized from the three external source images
python3 tools/digitize_cops_figures.py

# SHAP5 forward/fitting base and M004 altitude geometry table
env OMP_NUM_THREADS=4 ./liouville_cops fit_base.in
env OMP_NUM_THREADS=4 ./liouville_cops m004.in
```

The supplied `m004.in` differs from the fitting base in its ESA M004 DSK,
`OutputDirectory = output_M004`, `GeometryOnly = on`, and
`WriteKernels = off`. It calculates exact SPICE positions, inside flags, and
nearest-surface distances using the M004 mesh but skips unused source and gas
transport calculations. Its pressure columns are zero by design, so it is not
a forward-pressure result or fitting base. Review the external shape, observation, SPICE, and
HK-table paths in both input files before starting these potentially long
full-resolution runs; loading/building the four-million-plate M004 mesh still
takes time, while its 1,150 geometry evaluations are fast. The removed historical output directories are not
bundled inputs and are not recreated automatically.

Then use the collector defaults:

```bash
python3 tools/collect_cops_sources.py
```

Open the plots in interactive windows as well as writing the PNG files:

```bash
python3 tools/collect_cops_sources.py --show
```

Available options:

| Option | Default | Meaning |
|---|---|---|
| `--shap5-run DIR` | `output_fit_base` | LiouvilleCOPS SHAP5 run used for reference pressures and SHAP5 geometry. The directory must contain `COPS_model.dat`. |
| `--m004-run DIR` | `output_M004` | LiouvilleCOPS M004 run used for M004 altitude and cometocentric distance. |
| `--fit-run DIR` | `output_fit_base` | Base fit run whose `useNG` and `useRG` flags define the simulation/fit data selection. |
| `--fit-data FILE` | `data/COPS_fit_data.txt` | Output table containing the pressures and selection flags for `fit_cops.py --data`. Its PNG plot uses the same filename stem. |
| `--fit-used FILE` | `output_fit_base_selection/fit_data_used.txt` | Selection exported by the fitter, including a `--data-only` run, used to make the separate fit-profile plot. The plot is skipped if this file is absent. |
| `--apply-exclude-ng` | off | Apply the `ExcludeNG` intervals from the fit run instead of retaining otherwise valid NG points in those intervals. |
| `--exclude-ng-saturated` | off | Set `useNG=0` for the 13 saturated NG points near the end of the descent. |
| `--out FILE.json` | `data/COPS_EoM_sources.json` | Combined machine-readable output. |
| `--plots PREFIX` | `data/COPS_EoM_sources` | Prefix for the `_time.png` and `_altitude.png` plots. |
| `--show` | off | Open interactive Matplotlib windows. PNG files are written with or without this option. |

Example with explicit run directories and saturated NG points excluded:

```bash
python3 tools/collect_cops_sources.py \
  --shap5-run output_fit_base \
  --m004-run output_M004 \
  --fit-run output_fit_base \
  --exclude-ng-saturated \
  --show
```

The raw COPS directory is supplied only to `cops_hk_table.py`.
`collect_cops_sources.py` reads `data/COPS_HK_20160930.txt` from the output
root. Generate it first with `cops_hk_table.py`; no bundled copy remains.
The collector has no raw-data-directory option. Before creating any result,
it checks all required products and reports every missing pathname together
with the regeneration commands.

### Source-comparison results

The default run writes:

| Result | Contents |
|---|---|
| `data/COPS_EoM_sources.json` | All collected datasets, units, source information, geometry, and event metadata. |
| `data/COPS_EoM_sources_time.png` | NG and RG pressures from all sources versus UTC; red markers identify points selected for the fit. |
| `data/COPS_EoM_sources_altitude.png` | The same sources versus altitude. |
| `data/COPS_fit_data.txt` | Pressure values and `useNG`/`useRG` selection flags that can be passed to `fit_cops.py --data`. |
| `data/COPS_fit_data.png` | Selected and excluded fit data versus time and altitude. |
| `data/COPS_fit_profiles.png` | NG and RG profiles actually used by the fit, when `--fit-used` exists. |

The red candidate-fit markers are controlled by `--fit-run` and the selection
switches, not by `--fit-used`. The latter supplies only the separate profile
plot. Without `--apply-exclude-ng`, otherwise valid NG points in the base
run's `ExcludeNG` intervals are re-enabled. Saturated NG points are retained
unless `--exclude-ng-saturated` is supplied. Choose these switches explicitly
when comparing selections. A missing `--fit-used` file does not delete an old
`COPS_fit_profiles.png`, so check for stale plots.

## 3. Complete runnable inventory

All examples assume the current directory is the application root:

```bash
cd /home/vtenishe/ROS-END-OF-MISSION/AMPS/srcCG/LiouvilleCOPS
```

Replace that location on another machine. Shell examples use Bash-compatible
syntax. The `env OMP_NUM_THREADS=4 command ...` form also avoids needing a
shell-specific environment-variable assignment.

| Entry point | Purpose | How to invoke / result |
|---|---|---|
| `make` | Build all local `src/*.cpp` into one executable. | Section 4; produces `build/*.o` and `liouville_cops`. |
| `make clean` | Remove build products only. | Deletes `build/` and `liouville_cops`; does not remove run/data directories. |
| `liouville_cops selftest.in` | Internal analytic and geometry/transport checks. | Section 6; results on the terminal. |
| `liouville_cops cops.in` | Forward NG/RG pressure calculation. | Section 7; supplied destination is `output/`. |
| `liouville_cops fit_base.in` | Forward calculation with the cached geometry needed for fitting. | Section 7; supplied destination is `output_fit_base/`. |
| `liouville_cops FILE.in` | Any configured forward variant or fitted-model verification. | Sections 7 and 10; destination is the input's `OutputDirectory`. |
| `tools/cops_hk_table.py` | Convert raw L2 BG files to a housekeeping table. | Section 1; optional source directory and output filename, both with defaults. |
| `tools/rebuild_cops_boundary.py` | Regenerate a static fitting base for flat/rough emission without overwriting existing runs. | Section 7, "Change the surface boundary condition"; new table, elements, kernels, input snapshot, log, and metadata. |
| `tools/collect_cops_sources.py` | Compare pressure sources and export candidate fit data. | Section 2; JSON, PNGs, and selection table. |
| `tools/digitize_cops_figures.py` | Regenerate approximate pressure traces from the three supplied images. | Section 8; five CSV files and optional check plots. |
| `tools/fit_cops.py` | Select data, check the fit model/Jacobian, fit once, or scan regularization. | Section 9; included working file with shared output paths. |
| `tools/plot_cops_model.py` | Compare modeled and observed NG/RG pressures versus time or altitude. | Section 11; interactive figure or `--save` image. |
| `tools/plot_cops_model_altitude.py` | Same plotter, with altitude as the default horizontal axis. | Section 11; same arguments. |
| `tools/plot_lcurve.py` | Plot the regularization scan and diagnostics. | Sections 9 and 11; interactive figure or `--save` image. |

There is one C++ executable, not separate executables for different shapes,
gauges, or transport modes. The Makefile has no `test`, `install`, MPI, or
Python-fitting target. `__pycache__` files are interpreter caches, not
additional supported launch commands. Saved
`cops_fit.in` files under result directories are inputs to the same executable.

No scientific result directories are currently retained inside LiouvilleCOPS.
Output filenames are overwritten when an external destination is reused; use
new input filenames and output directories to preserve separate runs.
The commands below are examples, not a claim that the existing scientific
results have been regenerated.

## 4. Dependencies and compilation

### C++ executable: CSPICE and OpenMP

Required: `make`, a C++17 compiler with OpenMP support, and a CSPICE toolkit
with `include/SpiceUsr.h` and `lib/cspice.a`. The supplied Linux build uses:

| Make variable | Supplied default / role |
|---|---|
| `CSPICE` | `/home/vtenishe/SPICE/cspice`; toolkit root, not the mission-kernel directory. |
| `CXX` | `g++` unless overridden by make/environment settings. |
| `CXXFLAGS` | `-O3 -march=native -std=c++17 -fopenmp -Wall -Wextra -Wno-unused-parameter`. |
| `INCLUDES` | `-I$(CSPICE)/include`. |
| `LIBS` | `$(CSPICE)/lib/cspice.a -lm -fopenmp`. |

CSPICE is linked as a static archive. OpenMP is enabled at compilation and
linking; a GCC build normally also needs its OpenMP runtime (`libgomp`) on
the execution machine. This application does not require the surrounding
AMPS code to build, nor does it use MPI to distribute this calculation.

Build with the supplied toolkit path:

```bash
make -j4
```

Or set the toolkit and compiler explicitly:

```bash
make -j4 CSPICE=/path/to/cspice CXX=g++
```

Replace `/path/to/cspice` with a real toolkit directory. Changing `CSPICE`
does not change any input-file paths or tell SPICE where mission data live.
No download or installation of CSPICE is performed by this Makefile.

Make does not track changes to compiler flags or toolkit paths as file
dependencies. Force a rebuild after changing them:

```bash
make -B -j4 CSPICE=/path/to/cspice CXX=g++
```

The default `-march=native` optimizes for the build machine's CPU. To build
without that CPU-specific option, for example before copying the executable
to a different machine:

```bash
make -B -j4 CSPICE=/path/to/cspice CXX=g++ \
  CXXFLAGS='-O3 -std=c++17 -fopenmp -Wall -Wextra -Wno-unused-parameter'
```

Keep C++17 and OpenMP enabled when overriding the flags. There is no serial
build switch in this Makefile; run with one OpenMP thread if needed. A debug
compile can similarly use `-O0 -g -std=c++17 -fopenmp` in `CXXFLAGS`.

### Python tools

`cops_hk_table.py` uses the Python standard library.
`rebuild_cops_boundary.py` uses only the standard library with Python 3.8+,
plus the already compiled C++ executable. The collector and plotters
use NumPy and Matplotlib; fitting and image digitization also require SciPy.
There is no pinned Python environment supplied here. Use an existing
environment with these packages, or optionally create one:

```bash
python3 -m venv .venv
.venv/bin/python -m pip install numpy scipy matplotlib
```

For that environment, replace `python3` in the examples with
`.venv/bin/python`. Installing packages requires access to a package source.
The former RCS archives have been removed after verifying that every archived
file had a working counterpart. RCS tools are not required.

For servers without a display, use `--save` with the model/L-curve plotters
and omit `--show` from the collector. The digitizer writes files without
opening windows.

## 5. Paths and runtime input setup

### Four independent kinds of path

| Resource | Where to configure it | Current/example location |
|---|---|---|
| CSPICE software headers/library | `make CSPICE=...` | `/home/vtenishe/SPICE/cspice` |
| SPICE mission kernels | `KernelDirectory` and repeated `Kernel` in the `.in` file | `/home/vtenishe/SPICE/Kernels/ROSETTA/esa_2022` |
| External COPS data root | `COPS_DATA_ROOT` | `/data/vtenishe/rosetta_data_end_of_mission_cops` |
| Raw COPS L2 BG files | First argument or `COPS_L2_ROOT` | `<COPS_DATA_ROOT>/data_COPS_L2/mtp35_cops_bg` |
| L4 table and source images | collector `--sources-dir` or `COPS_SOURCES_ROOT` | `<COPS_DATA_ROOT>/sources` |
| Processed HK table | Second argument of the HK script; then `TemperatureRGTable` in the `.in` file | `data/COPS_HK_20160930.txt` |
| Observation times and reference pressures | `ObservationFile` in the `.in` file | `/home/vtenishe/ROS-END-OF-MISSION/AMPS/srcCG/RosinaMeasurements.h` |
| Surface mesh | `ShapeModel` and `ShapeModelFormat` in the `.in` file | `/home/vtenishe/ROS-END-OF-MISSION/AMPS/SHAP5_Volume.bdf`, `nastran` |
| Forward results | `OutputDirectory` in the `.in` file | `output` or `output_fit_base` in the supplied inputs |
| Fit results | Fitter `--out DIR` | Default is `<run-name>_fit` under the output root |

The observation file is read at runtime as a data table; its current location
outside LiouvilleCOPS is not a build dependency on another application.
The C++ observation reader expects the supported `RosinaMeasurements.h`
record format, not arbitrary CSV or the raw BG files. Moving the raw archive
does not automatically replace `ObservationFile` or change the observations.

The collector looks for HK and digitized products under the output root first.
It reads the L4 table from `--sources-dir`, whose default is
`<COPS_DATA_ROOT>/sources`. The digitizer reads its images from the same
external directory by default.

The current external defaults are:

```bash
export COPS_DATA_ROOT=/data/vtenishe/rosetta_data_end_of_mission_cops
# Optional finer-grained overrides:
# export COPS_L2_ROOT=/another/location/data_COPS_L2
# export COPS_SOURCES_ROOT=/another/location/sources
```
The fitter has a separate `--hk-table FILE` option for saturation screening;
this does not replace the gauge temperatures already saved in the base run.

To use the new raw archive location, run the command in section 1. Keep the
following entries in the forward input, or change the table path if you
deliberately wrote the processed table elsewhere:

```ini
TemperatureRGMode = hk_table
TemperatureRGTable = data/COPS_HK_20160930.txt
TemperatureRGTableColumn = 6
```

Column numbers are 1-based; column 6 is sensor temperature in kelvin. Raw
pressure fields are in mbar, and the supplied observation file also requires
`ReferenceDataUnit = mbar`. Forward result pressures are in Pa. Do not change
the unit setting merely because an output plot uses a different unit.

### Input-file syntax and path resolution

Files use `Key = Value`; `#` starts a comment. Edit actual entries in a copy
of `cops.in` or `fit_base.in`. The snippets here are modifications, not
complete standalone observation inputs.

- Paths are literal: do not put shell quotes, `$VARIABLE`, or `~` in an
  input-file path expecting shell expansion. Prefer paths without spaces
  for compatibility with the Python tools and table headers.
- Ordinary C++ source-data paths resolve against the launch directory, not
  the input-file directory. Exceptions: `OutputDirectory` uses the output root;
  `TemperatureRGTable` and `ActivityFile` use generated-product lookup above.
  An input filename not found as written also gets generated-product lookup,
  allowing `./liouville_cops output_fit_new/cops_fit.in`.
- Relative `Kernel` entries are joined to `KernelDirectory`; absolute ones
  are used directly. Merely setting the directory does not load its contents.
- A repeated scalar key uses its last value. Repeated `Kernel`, `ExcludeNG`,
  and `ExcludeRG` entries are meaningful lists. At most one `FitInterval`
  is supported.
- Unknown or unused keys are rejected in observation mode. For example,
  when switching to `TemperatureRGMode = constant`, remove/comment
  `TemperatureRGTable` and `TemperatureRGTableColumn`.
- Use a fresh `OutputDirectory`. Nested parent directories are created
  recursively; existing same-named products can still be overwritten.

Keep the original input used to generate a fitting base. The fitter reads its
name from the base table header and, for relative names, resolves it under
the application root. Do not move/delete that input or edit its exclusion
intervals after generating the base if you need a reproducible fit.

### SPICE runtime kernels

The supplied observation inputs explicitly load the following list. Retain
the list when moving the same kernel set, and change `KernelDirectory`:

```ini
KernelDirectory = /home/vtenishe/SPICE/Kernels/ROSETTA/esa_2022
Kernel = lsk/NAIF0011.TLS
Kernel = pck/PCK00010.TPC
Kernel = pck/ROS_CG_RAD_V10.TPC
Kernel = fk/ROS_V38.TF
Kernel = fk/RSSD0002.TF
Kernel = ik/ROS_ROSINA_V12.TI
Kernel = sclk/ROS_160929_STEP.TSC
Kernel = spk/DE405.BSP
Kernel = spk/CORB_DV_257_03___T19_00345.BSP
Kernel = spk/RORB_DV_257_03___T19_00345.BSP
Kernel = spk/ROS_STRUCT_V8.BSP
Kernel = ck/CATT_DV_257_03_______00344.BC
Kernel = ck/RATT_DV_257_02_04____00344.BC
Kernel = ck/ROS_SA_2016_V0042.BC
```

These files must exist and cover the observation times. Kernel availability
and time/frame coverage are separate requirements from successful linking.
Preserve/review `BodyFrame`, body/spacecraft names, `SpacecraftFrame`, and
the NG/RG instrument IDs when changing kernel sets. The checked-in input
notes an attitude-kernel revision difference from the cited ESA meta-kernel;
do not silently substitute a different revision and call it the same run.

The DSK alternatives are external files under the same SPICE root:

```ini
ShapeModelFormat = dsk
ShapeModel = /home/vtenishe/SPICE/Kernels/ROSETTA/esa_2022/dsk/ROS_CG_M004_OSPGDLR_N_V1.BDS
```

Use `ROS_CG_M001_OSPGDLR_N_V1.BDS` there for M001 instead. The DSK reader
expects its supported single type-2 segment with the requested body frame;
it converts DSK kilometers to meters. NASTRAN coordinates are assumed to be
meters. Keep `RequireClosedSurface = on` unless deliberately investigating
mesh topology; changing it is not a substitute for fixing an invalid mesh.

## 6. Smoke test and parallel execution

After building:

```bash
env OMP_NUM_THREADS=4 ./liouville_cops selftest.in
```

This uses a generated sphere and internal analytic/geometry/transport checks;
it does not need the external shape, observation table, raw COPS archive, or
mission SPICE kernels at runtime. The executable still has to be built with
CSPICE and OpenMP. Inspect the terminal summary and exit status; passing is
a numerical implementation check, not validation of the fitted science model.

The supplied self-test settings are `SelfTestRadius = 1000`,
`SelfTestRefinement = 6`, `SelfTestTolerance = 0.01`,
`SubdivisionAccuracy = 0.05`, and `MaxSubdivisionLevel = 8`.
Optional `SelfTestRays` defaults to 2000. Make a separate input to change them.

Runtime thread control is independent of compilation parallelism:

| Control | Meaning |
|---|---|
| `make -j4` | Up to four build jobs; no control over the running model. |
| `env OMP_NUM_THREADS=4 ./liouville_cops FILE.in` | OpenMP thread count for the C++ run when `Threads` is absent. |
| `Threads = 4` in `FILE.in` | Explicit OpenMP thread count; overrides the environment setting. |
| Fitter `--jobs 2` | Up to two independent regularization fits; not C++ OpenMP threads. |

The executable prints its OpenMP maximum thread count on startup. Use
`OMP_NUM_THREADS=1` for a one-thread run; no `mpiexec` wrapper is needed.

## 7. Run the forward model and generate a fitting base

The complete C++ command-line interface is:

```bash
./liouville_cops INPUT.in
```

It takes exactly one filename. Settings, including mode, output destination,
and paths, go inside that file; there are no `--output`, `--threads`, or
`--cspice` command-line options.

After reviewing the paths, the supplied full-run commands are:

```bash
env OMP_NUM_THREADS=4 ./liouville_cops cops.in
env OMP_NUM_THREADS=4 ./liouville_cops fit_base.in
```

These are alternative runs, not two required stages. They write under
`<AMPS>/cops_output/output/` and `output_fit_base/`, respectively. To create a
separately named base:

```bash
cp -i fit_base.in fit_base_new.in
```

Edit `fit_base_new.in`: set paths as in section 5, change
`OutputDirectory = output_fit_base_new`, and retain these fitting settings:

```ini
Mode = observations
WriteKernels = on
FlightTimeRotation = none
TemperatureModel = analytical
SurfaceTemperatureCorrectionFactor = 1
```

Then run, optionally retaining the terminal log:

```bash
mkdir -p "${COPS_OUTPUT_ROOT:-../../cops_output}"
env OMP_NUM_THREADS=4 ./liouville_cops fit_base_new.in \
  > "${COPS_OUTPUT_ROOT:-../../cops_output}/fit_base_new.log" 2>&1
```

Check the exit status and log before fitting. Preserve `fit_base_new.in`.
For a short trial, use a different input/output and a small range of
`FirstPoint`, `LastPoint`, and `PointStep`. The indices are zero-based and
the last index is inclusive; supplied full inputs use 0 through 1149, step 1.
A reduced run is not a replacement for a full-data fitting base.

### Runtime controls and forward variants

| Input keys | Purpose / choices |
|---|---|
| `GeometryOnly` | `off` normally. `on` skips source/transport evaluation and writes only usable geometry/reference/selection fields; pressure/source fields are zero. It requires `WriteKernels = off` and the result cannot be fitted. Used by `m004.in` for the collector. |
| `FitInterval`, `ExcludeNG`, `ExcludeRG` | UTC intervals defining candidate fitting flags, not a request to skip the forward calculation. Endpoints are inclusive. Python NG overrides are discussed in section 9. |
| `ActivityModel` | `shap5`, `shap5_1`, or `file`; the first two use local compiled coefficient tables. |
| `ActivityFile`, `ActivityFileHeliocentricModel` | For `ActivityModel = file`: three rows of 25 coefficients (H2O, CO2, CO), and the `shap5`/`shap5_1` heliocentric scaling to use. |
| `SZAMode`, `CastShadows` | `continuous` or `srccg_table` solar-angle treatment; shadow ray tracing on/off. |
| `TemperatureModel` | `analytical`, `constant`, or `bjorn`. Only the analytical base is supported by the current fitter. |
| `AnalyticalTmax`, `AnalyticalT75`, `AnalyticalTmin` | Three analytical surface-temperature parameters in kelvin. |
| `ConstantTemperature`, `BjornTemperatureDistanceIndex` | Temperature or table index for the corresponding alternative temperature model. |
| `SurfaceTemperatureCorrectionFactor` | Temperature multiplier; retain 1 for Python/C++ fitting consistency. |
| `SurfaceBC` | `flat` (Lambertian) or `rough` (uniform hemispherical emission). |
| `SubdivisionAccuracy`, `MaxSubdivisionLevel` | Adaptive surface-integration accuracy and maximum refinement level. |
| `SpacecraftVelocity` | Enable/disable the spacecraft-velocity contribution to RG flux. |
| `FlightTimeRotation` | `none` for the static calculation and fit base; `exact` for numerical flight-time rotation. |
| `SpeedNodes`, `SpeedRangeMax`, `PathTolerance` | Exact-rotation velocity quadrature and path approximation; supplied values 16, 5.5, and 0.5 m. |
| `NGAcceptance`, `NGHeadOnWeight` | NG acceptance: `srccg_shield`, `hemisphere`, `ik_fov`, or `full`; optional front-incidence weighting. |
| `NGShieldRmax`, `NGShieldHmax`, `NGShieldRtube` | Shield geometry ratios; supplied values 10.3, 33, 20.5. |
| `NGShieldTableBins`, `NGShieldTableSamples` | Optional shield lookup resolution, defaults 2001 and 20000; only used in shield mode. |
| `RGAcceptance`, `RGPressureModel` | Acceptance `hemisphere`/`ik_fov`; pressure conversion `equilibrium_sphere`, `calibration_srccg`, or `dsmc_srccg`. Supplied pressure conversion is `equilibrium_sphere`. |
| `TemperatureNG`, `TemperatureRG`, `TemperatureRGMode` | NG gauge temperature; constant RG temperature or `hk_table` RG mode. These are not the surface-temperature parameters. |
| `SpeciesSensitivity`, `SensitivityH2O`, `SensitivityCO2` | Instrument species-sensitivity correction on/off and its H2O/CO2 factors. |
| `GaugeOffsetSC` | Common gauge position offset, three components in spacecraft-frame meters. |

Changing geometry, shadows, surface boundary condition, acceptance, fixed
gauge temperature, sensitivities, or other cached settings requires a new,
compatible fitting base. Changing a key in the old input does not update the
saved binary kernels. Use distinct output directories for shape, temperature,
or rotation comparisons. Names such as `srccg_table` refer to local model
conventions, not execution of the surrounding srcCG application.

### Change the surface boundary condition

Use the original forward input to rebuild with the other angular law. The
script recomputes a complete static fitting base, not just the pressure table:

```bash
python3 tools/rebuild_cops_boundary.py --input fit_base.in \
  --surface-bc rough --out output_fit_base_rough --threads 4 --dry-run

python3 tools/rebuild_cops_boundary.py --input fit_base.in \
  --surface-bc rough --out output_fit_base_rough --threads 4
```

The first command prints the proposed input and command without writes or
external-file checks. The second executes the existing solver. Use `flat`
instead of `rough` and another fresh output path to rebuild the flat case.

| Option | Meaning |
|---|---|
| `--input FILE` | Source input; default `fit_base.in`. All its ordinary relative data paths are interpreted under the application root, not its own directory. |
| `--surface-bc flat\|rough` | Required target angular emission law. |
| `--out DIR` | Required new output directory; any existing path is refused. |
| `--threads N` | Optional positive OpenMP thread count; otherwise preserve the source input/environment. |
| `--executable FILE` | Existing solver; default application-root `liouville_cops`. No compilation is performed. |
| `--dry-run` | Show configuration and command only. |

Relative `--out` resolves under the output root. Other relative CLI paths
resolve against the application root, with generated-input fallback. The tool
preserves temperatures, source parameters, and selection intervals, sets
`WriteKernels = on`, and saves normalized absolute paths in
`<out>/boundary.in`. It requires observation mode, static transport,
analytical surface temperature, and temperature correction factor 1;
incompatible inputs are rejected, not silently converted.

Results are `<out>/COPS_model.dat`, `elements.bin`, and `kernels/`, plus
`source.in`, `boundary.in`, `rebuild.log`, and `rebuild.json`. It checks the
table's full expected point sequence and version-1 cache sizes/headers after
successful solver exit. Keep `boundary.in` in place: the fitter reads that
recorded input. On failure the tool retains partial output for diagnosis;
do not fit it, and retry with a new destination.

For a fitted model with the new boundary condition, rerun the fitting workflow
in section 9 with `output_fit_base_rough` as the base, for example:

```bash
python3 tools/fit_cops.py output_fit_base_rough \
  --apply-exclude-ng --exclude-ng-saturated --data-only \
  --out output_rough_selection

python3 tools/fit_cops.py output_fit_base_rough \
  --apply-exclude-ng --exclude-ng-saturated --ordered-temperature \
  --lambda 3e-3 --out output_rough_fit
```

Restore the fitter first if absent, as described in section 9. Compare flat
and rough using the same data selection and optimizer settings. `3e-3` is
an illustrative regularization weight, not an established optimum.

The angular kernels `Kn`, `KPhi`, and `KV` must change with the boundary
condition. Surface temperature is not part of those static kernels: it
enters density and the RG velocity correction through `1/sqrt(T_surface)`.
Gauge temperatures enter pressure conversion separately. The tool keeps
temperature settings unchanged; new best-fit temperatures, activities, and
CO2 scale require a new optimization. See README section 6.6 for the full
dependency explanation and why exact-rotation kernels are different.

Tests of this driver (no scientific-data dependency by default):

```bash
python3 -m unittest discover -s tests -p 'test_rebuild_cops_boundary.py' -v
```

### Forward output files

| Under `OutputDirectory` | Contents / use |
|---|---|
| `COPS_model.dat` | Time, geometry, reference and predicted NG/RG pressures, selection flags, gauge temperatures, and other model diagnostics. Plotters read this file. |
| `elements.bin` | Surface-element geometry/source data for fitting, written with kernel output enabled. |
| `kernels/point=*.bin` | Per-observation transport kernels, written with `WriteKernels = on`. |

Fitting requires all three, not just the model table. It loads kernels even
for rows subsequently excluded from the misfit. Static version-1 kernels
are supported by the fitter; exact-rotation version-2 kernels are not.
Keep the input, table, elements, and kernels together as a consistent set.

## 8. Regenerate digitized comparison data

This is optional if the supplied CSV files already exist:

```bash
python3 tools/digitize_cops_figures.py --check
```

| Option | Meaning |
|---|---|
| `--dir DIR` | Input image directory; default is `<COPS_DATA_ROOT>/sources` (or `COPS_SOURCES_ROOT`). |
| `--out-dir DIR` | Generated products directory; default is `data/sources` under the output root. |
| `--check` | Also write the two visual diagnostic PNGs; CSVs are written with or without it. |

Required images are `EGU2019_Tenishev_slide02.jpg`,
`EGU2019_Tenishev_slide16_data.jpg`, and
`ESA_blog_2016-12-15_COPS_EoM.jpg`. Outputs under `--out-dir` are:

- `digitized_EGU2019_Tenishev_slide02_NG.csv` and `_RG.csv`;
- `digitized_EGU2019_Tenishev_slide16_data_NG.csv` and `_RG.csv`;
- `digitized_ESA_blog_2016-12-15_COPS_EoM_NG.csv`;
- with `--check`, `digitized_check.png` and
  `EGU2019_slide02_on_slide16_axes.png`.

The pixel calibrations and color masks are specific to these images. This
is not a general image-to-data converter. Digitized traces are approximate
comparison sources, not recovered raw instrument measurements. Changing
`--dir` changes only image input. Keep the default `--out-dir` for the
collector to discover newly digitized CSVs automatically.

## 9. Data selection, fitting, derivative checks, and L-curves

### Check the fitting entry point

`tools/fit_cops.py` is included with shared-output path handling:

```bash
python3 tools/fit_cops.py --help
```

All listed Python entry points are present as ordinary working files.

### Prerequisites and selection checks

The examples below use `output_fit_base_new` generated in section 7. You may
instead use a complete compatible saved base, but do not mix tables and
kernels from different runs. The fitter checks that its initial evaluation
reproduces the base C++ pressures before optimizing.

Export the selected data without optimization:

```bash
python3 tools/fit_cops.py output_fit_base_new \
  --data-only --exclude-ng-saturated --out output_selection_new
```

Inspect `output_selection_new/fit_data_used.txt`, `fit_data_used.png`, and
`fit_profiles_used.png`. `--data-only` still loads/checks the fit model and
kernels; it is not a lightweight table-only command.

Selection rules matter:

- Without `--data`, the fitter normally re-enables otherwise valid NG
  points in the base input's `ExcludeNG` intervals. Add `--apply-exclude-ng`
  to retain those exclusions.
- Saturated NG plateau points are retained unless `--exclude-ng-saturated`
  is supplied. The HK table is read even when that switch is absent.
- `--data FILE` uses that table's pressure values and selection flags;
  matching observation IDs, UTC strings, and saved RG temperatures are
  required. The explicit exclusion switches can still remove NG points.
- Invalid/nonpositive selected pressures and incompatible tables must be
  corrected rather than treated as genuine observations. Consult README
  before changing missing-data or reconstructed-data handling.

To make the multi-source plot match a deliberate candidate selection and
use exactly its exported table in the fit:

```bash
python3 tools/collect_cops_sources.py \
  --fit-run output_fit_base_new \
  --apply-exclude-ng --exclude-ng-saturated \
  --fit-data data/COPS_fit_data_new.txt \
  --fit-used output_selection_new/fit_data_used.txt \
  --out data/COPS_EoM_sources_new.json --plots data/COPS_EoM_sources_new

python3 tools/fit_cops.py output_fit_base_new \
  --data data/COPS_fit_data_new.txt --data-only \
  --out output_selection_table_new
```

Here the collector still uses its default SHAP5 and M004 geometry runs.
Set `--shap5-run`/`--m004-run` if replacing those runs with matching shapes.
The first command's `--fit-used` is an independently exported selection;
to display the table-based selection instead, rerun the collector pointing
that option to `output_selection_table_new/fit_data_used.txt`.

### One fit

Example using the exported candidate table:

```bash
python3 tools/fit_cops.py output_fit_base_new \
  --data data/COPS_fit_data_new.txt --lambda 3e-3 \
  --ordered-temperature --out output_fit_new
```

Or, independently, fit from the base reference values with explicit NG
selection switches:

```bash
python3 tools/fit_cops.py output_fit_base_new \
  --apply-exclude-ng --exclude-ng-saturated --lambda 3e-3 \
  --ordered-temperature --out output_fit_from_base_new
```

These lambda values are illustrative, not an established optimal choice.
Read `fit_result.txt` for termination status, bound effects, residuals,
conditioning, and parameters. Reaching `--max-nfev` is not convergence.
The default free quantities are the H2O activity coefficients through
degree 4, three analytical surface-temperature parameters, and a CO2
activity multiplier; this is not a fit of arbitrary per-point pressures.

### Derivative check

```bash
python3 tools/fit_cops.py output_fit_base_new \
  --check-jacobian --lambda 1e-3 --out output_jacobian_check_new
```

This compares analytic derivatives with finite differences at the initial
parameters, prints the discrepancies, and stops before optimization. Inspect
the numbers: this mode does not implement an independent pass/fail threshold
for the derivative comparison. `--check-tolerance` controls the prior
pressure-consistency check, not the Jacobian difference.

### Regularization scan

```bash
env OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  python3 tools/fit_cops.py output_fit_base_new \
  --data data/COPS_fit_data_new.txt --ordered-temperature \
  --lambda 1e-4 1e-3 3e-3 1e-2 1e-1 \
  --jobs 2 --out output_lcurve_new

python3 tools/plot_lcurve.py output_lcurve_new \
  --save output_lcurve_new/lcurve.png
```

Each lambda is an independent optimization starting from the base parameters.
`--jobs` controls concurrent Python processes. The environment settings limit
threads for the corresponding numerical backends if used, reducing nested
parallelism. Kernel arrays can consume substantial RAM; use `--jobs 1` when
memory is constrained. The automatic maximum-curvature point is a diagnostic,
not a guarantee of convergence or scientific validity.

### Complete fitting options

| Parameter | Default | Meaning |
|---|---|---|
| Positional `run` | Required | Compatible base-run directory containing table, elements, and kernels. |
| `--out DIR` | `<run-name>_fit` | Results directory under the output root; created if needed, existing same-named products overwritten. |
| `--degree N` | 4 | Fitted harmonic groups through degree N, intended range 0–4; higher coefficients stay at base values. |
| `--lambda L [L ...]` | 0 | Nonnegative regularization weight(s); multiple values create a scan. Zero means no regularization penalty. |
| `--jobs N` | Minimum of number of lambdas and 12 | Concurrent scan fits. |
| `--fix-temperature` | off | Keep all three base surface-temperature parameters. |
| `--ordered-temperature` | off | Enforce Tmax >= T75 >= Tmin through nonnegative increments. |
| `--fix-kappa` | off | Keep the CO2 multiplier at 1. |
| `--tbounds LO HI` | 30 400 | Independent temperature bounds in kelvin, or only Tmin bounds in ordered mode; ordered Tmax/T75 can exceed HI. |
| `--max-nfev N` | 1000 | Maximum objective-function evaluations. |
| `--ftol VALUE` | 1e-10 | Objective-change stopping tolerance. |
| `--xtol VALUE` | 1e-10 | Step-change stopping tolerance. |
| `--gtol VALUE` | 1e-10 | Gradient/optimality stopping tolerance. |
| `--check-tolerance VALUE` | 1e-6 | Maximum relative discrepancy allowed between the prior Python and base C++ pressures. |
| `--check-jacobian` | off | Print derivative comparison and stop. |
| `--data FILE` | Base table observations | Replacement pressure/selection table with matching identifiers and gauge temperatures. |
| `--apply-exclude-ng` | off | Apply base-input NG exclusion intervals. |
| `--exclude-ng-saturated` | off | Remove NG points flagged by the HK saturation check. |
| `--hk-table FILE` | `data/COPS_HK_20160930.txt`, output-root first | HK source for saturation checking, not a way to change cached T_RG. |
| `--data-only` | off | Export selection plots/table after consistency checks, then stop. |

### Fit output locations

| Location | Contents |
|---|---|
| `<out>/fit_data_used.txt`, `fit_data_used.png`, `fit_profiles_used.png` | Actual selection exports, including in data-only and derivative-check modes. |
| Single fit: `<out>/fit_result.txt` | Parameters, termination and misfit/regularization diagnostics. |
| Single fit: `<out>/ActivityFit.txt` | Three species rows of 25 activity coefficients. |
| Single fit: `<out>/cops_fit.in` | Generated input for a full C++ verification run. |
| Single fit: `<out>/COPS_model.dat` | Fitted total NG/RG predictions in the base table layout. |
| Scan: `<out>/lambda=<value>/` | The four optimization products above for each lambda; directory values use `%.3e` formatting. |
| Scan: `<out>/lcurve.txt` | Aggregate scan diagnostics; selection exports remain in `<out>`. |

The Python-fitted table is not a complete forward-model recomputation:
species pressure/density/flux columns are set to -1, while other fields such
as production rates and much of the header are inherited. Use the C++
verification output for recomputed species results and production rates.

## 10. Verify fitted parameters with C++ and test flight-time rotation

For the newly generated single-fit example, inspect `output_fit_new/cops_fit.in`
and its paths, then run:

```bash
env OMP_NUM_THREADS=4 ./liouville_cops output_fit_new/cops_fit.in
```

The generated input normally sets `ActivityModel = file`, points to that
fit's `ActivityFit.txt`, inserts the fitted temperature parameters, and sets
the destination to `output_fit_new_verify`. Check it before running; existing
verification files would be overwritten. For a scan result, use the
`cops_fit.in` inside the chosen `lambda=...` directory; its verification
destination is that fit-directory path with `_verify` appended.

The C++ rerun uses the original observation file and copied input exclusions;
it does not automatically reproduce a replacement `--data` table or the
Python selection switches. `fit_data_used.txt` is the authoritative record
of the data used in the fit.

For an exact-rotation forward check:

```bash
cp -i "${COPS_OUTPUT_ROOT:-../../cops_output}/output_fit_new/cops_fit.in" fit_exact_new.in
```

Edit the copy to use a distinct output destination and exact transport:

```ini
FlightTimeRotation = exact
OutputDirectory = output_fit_new_exact
SpeedNodes = 16
SpeedRangeMax = 5.5
PathTolerance = 0.5
```

Then run it:

```bash
env OMP_NUM_THREADS=4 ./liouville_cops fit_exact_new.in
```

Review numerical convergence before interpreting differences. The current
Python fitter cannot fit the exact-rotation kernels; this is a forward
comparison using the already fitted parameters.

Historical saved `cops_fit.in` files reference an older absolute application
location. Correct `ActivityFile`, `OutputDirectory`, and external input paths
before using them here. A saved input is not automatically portable.

## 11. Plot model runs and fit diagnostics

Compare the base, fitted totals, and full C++ verification:

```bash
python3 tools/plot_cops_model.py \
  output_fit_base_new output_fit_new output_fit_new_verify \
  --ref output_fit_base_new --save COPS_fit_comparison_new.png
```

Plot against altitude instead:

```bash
python3 tools/plot_cops_model_altitude.py \
  output_fit_base_new output_fit_new_verify \
  --ref output_fit_base_new --save COPS_fit_altitude_new.png
```

| Model-plotter argument | Meaning |
|---|---|
| Positional `runs ...` | One or more run directories or explicit `COPS_model.dat` files. Omit to scan output-root `output*/COPS_model.dat`, nonrecursively; if none exist, scan legacy application-root runs. |
| `--x time` / `--x altitude` | Horizontal axis; defaults to time in the main script and altitude in its wrapper. |
| `--save FILE` | Save the figure instead of opening a window. |
| `--theme light` / `--theme dark` | Plot style; default light. |
| `--ref NAME` | Loaded run's directory basename, or its suffix after `output_`, supplying observations/selection and reference geometry. Default is newest loaded table with useNG/useRG flags. |

Set `--ref` explicitly for reproducibility; it is not a path to an unloaded
run. The four panels show combined NG/RG, NG, RG, and their ratio; unused
reference points are hollow. Nested scan directories must be passed
explicitly, for example:

```bash
python3 tools/plot_cops_model.py output_fit_base_new \
  output_lcurve_new/lambda=3.000e-03 \
  --ref output_fit_base_new --save COPS_scan_fit_new.png
```

For the L-curve plotter, the positional argument is a scan directory or its
`lcurve.txt`; `--save FILE` is its only additional option besides help:

```bash
python3 tools/plot_lcurve.py output_lcurve_new/lcurve.txt \
  --save output_lcurve_new/lcurve_diagnostics.png
```

Omit `--save` for an interactive window. Open circles identify scan fits
that have not converged; inspect the status and evaluation counts rather
than assuming the marked curvature point is an accepted fit. Parent directories
for custom plot filenames are created automatically.

All Python entry points provide `--help`, including `cops_hk_table.py`.
The C++ executable requires an input filename;
`--help` is not a special C++ option.

## 12. Troubleshooting and reproducibility checklist

| Symptom | Check |
|---|---|
| Missing `SpiceUsr.h` or `cspice.a` | Toolkit root passed to make; verify its `include/` and `lib/` layout, then rebuild with `-B`. |
| OpenMP compiler/link/runtime error | Compiler's OpenMP support and runtime; keep `-fopenmp` in both compile and link settings. |
| SPICE cannot load a file or resolve geometry | Runtime `KernelDirectory`/`Kernel` paths, exact filenames, time coverage, frames, and instrument IDs; changing build-time `CSPICE` does not fix these. |
| Missing shape, observation, or HK file | Literal paths and launch directory; no shell expansion occurs in `.in` values. |
| Unknown/unused input key | Spelling and selected mode, including conditional HK/activity/shield options. |
| `fit_cops.py` not found | Recover it from source control or another complete LiouvilleCOPS copy; no local RCS archive remains. |
| Fitter missing columns or kernels | Generate a complete current static base with `WriteKernels = on`; an old table alone is insufficient. |
| Fitter consistency check fails | Check paired table/elements/kernels, analytical temperature with correction factor 1, and preserved source/input files. Do not simply loosen tolerance to hide a mismatch. |
| Fit or selection differs from saved results | Check NG exclusion/saturation switches and `--data`; current defaults need not match historical selections. |
| Collector plots appear inconsistent | Check `--fit-run` versus `--fit-used`, SHAP5/M004 geometry inputs, and stale profile PNGs. |
| Plot fails on a headless host | Use `--save` for model/L-curve plots and omit collector `--show`. |
| Run or scan exhausts resources | Reduce runtime threads or fitter `--jobs`; use a separate small forward trial before a full high-resolution run. |

Preserve each run's exact input, source-data provenance, kernel filenames,
build settings, logs, selection exports, and fit options. Keep forward base,
fit, verification, and exact-rotation outputs in separate directories.
Consult [README.md](README.md) before interpreting historical convergence,
NG reconstruction, production rates, or fitted parameter uncertainties.
