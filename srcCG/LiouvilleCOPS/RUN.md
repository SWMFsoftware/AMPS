# LiouvilleCOPS: build and run guide

This guide covers only this standalone application and its subdirectories.
Sections 1–2 describe the multi-source pressure plot; sections 3 onward cover
compilation, paths, the forward model, fitting, verification, and every Python
entry point. See [README.md](README.md) for the model equations, data audit,
file formats, and scientific limitations.

## Output location

For the restored 788k baseline and joint NG-only/cliff fitting workflow, use
[RESTORED_BASELINE.md](RESTORED_BASELINE.md) and
[JOINT_TEMPLATE_FIT.md](JOINT_TEMPLATE_FIT.md). Read
[GRAVITY_REPAIR.md](GRAVITY_REPAIR.md) before rebuilding cliff products: the
original HEAL-derived gravity and cliff masks are invalid and quarantined.
The older degree-8 fixed-baseline examples below are separate experiments.

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

### Run and fit integrity records

New C++ observation runs require a new `OutputDirectory`; an existing path is
refused before mission data are loaded. On success the directory contains:

- `input.snapshot.in`, the exact input bytes used;
- `run_manifest.json`, schema `liouville_cops.run_manifest` version 1;
- `COPS_model.dat` and, when requested, `elements.bin` plus `kernels/`;
- `ActivityEnhancementFacets.dat` only when the optional Stage-A3 residual is enabled.

Before writing the manifest, the solver reads its products back and checks the
point sequence, binary versions and exact lengths, face ranges and duplicates,
and stale extra kernels. The manifest records streamed SHA-256 identities for
the executable, input, shape, observations, SPICE kernels, optional HK/activity
files, and outputs. `process`, `cache`, and `scientific_review` are distinct;
a completed/verified run remains `scientific_review=not_reviewed`.

New fits likewise require a new `--out` directory and write
`fit_manifest.json`. Before creating it, the fitter independently audits the
entire static cache and verifies that its prior reproduces the base C++
pressures. A pre-W01 base without `run_manifest.json` is refused unless
`--allow-legacy-base` is supplied explicitly.

Validate a base without changing it:

```bash
python3 tools/validate_cops_run.py output_fit_base
```

After generating a C++ verification run, register pressure agreement:

```bash
python3 tools/validate_cops_run.py \
  --fit output_fit_new --verification output_fit_new_verify
```

This changes only the fit manifest's C++-agreement state. It never marks the
optimizer converged or the result scientifically reviewed.

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
| `tools/rebuild_cops_boundary.py` | Regenerate a static fitting base for a selected shape and flat/rough emission without overwriting existing runs. | Section 7, "Change the surface boundary condition or shape model"; new table, elements, kernels, input snapshot, log, and metadata. |
| `tools/collect_cops_sources.py` | Compare pressure sources and export candidate fit data. | Section 2; JSON, PNGs, and selection table. |
| `tools/digitize_cops_figures.py` | Regenerate approximate pressure traces from the three supplied images. | Section 8; five CSV files and optional check plots. |
| `tools/fit_cops.py` | Select data, check the fit model/Jacobian, fit once, or scan regularization. | Section 9; included working file with shared output paths. |
| `tools/run_fit_profile.py` | Run one or all versioned data-selection profiles through the same fitter. | Section 9; defaults to `--jobs 14`, distinct outputs, and one numerical thread per fit process. |
| `tools/run_a6_validation.py` | Run the Stage-A6 unit/synthetic gate, C++ self-test, and optional local integration gate. | Section 6; terminal results and optional JSON report. |
| `tools/benchmark_fit_performance.py` | Characterize Stage-A7 production-fitter time, memory, output size, and independent-fit scaling. | Section 6; terminal tables and optional new JSON/Markdown reports. |
| `tools/convert_dsk_to_bdf.py` | Convert archived type-2 DSK meshes to the original SHAP5 NASTRAN surface convention with topology validation. | Section 5; new BDFs and conversion manifests beside external shape data. |
| `tools/benchmark_activity_basis.py` | Build and time the isolated degree-5--8 Stage-A2 basis on a version-1 element cache. | Command below; prints JSON and does not modify the cache. |
| `tools/validate_cops_run.py` | Strictly validate a static cache or register C++ reproduction of fitted pressures. | Run on a base, or with `--fit FIT --verification RUN`; see the integrity section above. |
| `tools/plot_cops_model.py` | Compare modeled and observed NG/RG pressures versus time or altitude. | Section 11; interactive figure or `--save` image. |
| `tools/plot_cops_model_altitude.py` | Same plotter, with altitude as the default horizontal axis. | Section 11; same arguments. |
| `tools/plot_lcurve.py` | Plot the regularization scan and diagnostics. | Sections 9 and 11; interactive figure or `--save` image. |

There is one C++ executable, not separate executables for different shapes,
gauges, or transport modes. The Makefile has no `test`, `install`, MPI, or
Python-fitting target. `__pycache__` files are interpreter caches, not
additional supported launch commands. Saved
`cops_fit.in` files under result directories are inputs to the same executable.

No scientific result directories are currently retained inside LiouvilleCOPS.
C++ runs and fits refuse existing result directories. Preprocessing, collector,
and plot filenames can still be overwritten when explicitly reused, so use
new names when those products must be preserved.
The commands below are examples, not a claim that the existing scientific
results have been regenerated.

### Validate and benchmark the Stage-A2 harmonic basis

The fast test suite compiles a standard-library-only C++ probe and compares it
with the Python evaluator for every mode through degree 8:

```bash
python3 -m unittest tests.test_a2_activity_basis -v
```

Benchmark degrees 5--8 on a real fitting base without changing it:

```bash
python3 tools/benchmark_activity_basis.py \
  output_fit_full_default_base/elements.bin \
  --min-degree 5 --max-degree 8 --repeat 5
```

The optional `--json FILE` writes the report under `<AMPS>/cops_output` (or
`COPS_OUTPUT_ROOT`). The tool validates the version-1 header and exact file
size, memory-maps only the centroid columns, constructs a transient dense
table, reports exact table bytes plus measured maximum-RSS change, and leaves
`elements.bin` untouched. `HARMONIC_BASIS.md` is the normative convention.
Stage A3 connects this basis to an optional C++ H2O forward residual; the
Python fitter still does not estimate it.

Validate the strict Stage-A3 coefficient format, baseline control, a nonzero
single mode, final clipping, and incompatibility failures:

```bash
python3 -m unittest tests.test_a3_activity_enhancement -v
```

To exercise the enabled forward path, copy `cops.in`, choose a fresh
`OutputDirectory`, and use the checked-in zero-residual example:

```ini
ActivityEnhancementModel = real_spherical_harmonics_v1
ActivityEnhancementDegree = 5
ActivityEnhancementFile = examples/activity_enhancement_zero_l5.txt
ActivityEnhancementScale = 1.0
```

```bash
OMP_NUM_THREADS=14 ./liouville_cops my_enhanced_run.in
```

The example is tied to `ActivityModel = shap5` and body frame `67P/C-G_CK`.
It must reproduce baseline pressures while writing
`ActivityEnhancementFacets.dat`. A nonzero file uses the same command. See
`ACTIVITY_ENHANCEMENT_FILE.md` for required metadata and coefficient records.
Stage A3 applies such a file. Stage A4 also lets `fit_cops.py` estimate the
degree-5--8 H2O residual under the fixed-legacy/fixed-temperature/fixed-CO2
contract documented below.

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
| SHAP5 surface mesh | `ShapeModel` in the `.in` file | `/data/vtenishe/rosetta_data_end_of_mission_cops/SHAP5_Volume.bdf` |
| Raw COPS L2 BG files | First argument or `COPS_L2_ROOT` | `<COPS_DATA_ROOT>/data_COPS_L2/mtp35_cops_bg` |
| L4 table and source images | collector `--sources-dir` or `COPS_SOURCES_ROOT` | `<COPS_DATA_ROOT>/sources` |
| Processed HK table | Second argument of the HK script; then `TemperatureRGTable` in the `.in` file | `data/COPS_HK_20160930.txt` |
| Observation times and reference pressures | `ObservationFile` in the `.in` file | `/home/vtenishe/ROS-END-OF-MISSION/AMPS/srcCG/RosinaMeasurements.h` |
| Surface mesh | `ShapeModel` and `ShapeModelFormat` in the `.in` file | `/data/vtenishe/rosetta_data_end_of_mission_cops/SHAP5_Volume.bdf`, `nastran` |
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
`COPS_DATA_ROOT` configures the Python raw/L4/image tools. The C++ solver reads
the literal `ShapeModel` value, so update that entry in `cops.in` and
`fit_base.in` if the external data root is moved.
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
  `TemperatureRGTable`, `ActivityFile`, and `ActivityEnhancementFile` use generated-product lookup above.
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
- Use a fresh `OutputDirectory`. Nested parent directories are created, but
  any existing destination is refused before the run begins.

New runs preserve `input.snapshot.in`; the fitter uses that snapshot for
selection intervals and its generated verification input. Legacy runs lack
the snapshot and fall back to the original input name recorded in the table,
so moving or editing that file makes legacy provenance incomplete.

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

### Convert archived SHAP5 DSKs to the original BDF convention

The original `SHAP5_Volume.bdf` uses metre-unit NASTRAN `GRID*` and `CTRIA3`
cards. Its `CTETRA` volume cards are ignored by LiouvilleCOPS, so the converted
products contain the complete surface cards but do not invent volume cells.
Convert the downloaded PDS models with the authoritative NAIF readers:

```bash
python3 tools/convert_dsk_to_bdf.py \
  /data/vtenishe/rosetta_data_end_of_mission_cops/shape_models/SHAP5_PDS_V2/cg_spc_shap5_199k_cart.bds \
  /data/vtenishe/rosetta_data_end_of_mission_cops/shape_models/SHAP5_PDS_V2/cg_spc_shap5_384k_cart.bds \
  /data/vtenishe/rosetta_data_end_of_mission_cops/shape_models/SHAP5_PDS_V2/cg_spc_shap5_788k_cart.bds \
  /data/vtenishe/rosetta_data_end_of_mission_cops/shape_models/SHAP5_PDS_V2/cg_spc_shap5_001m_cart.bds \
  /data/vtenishe/rosetta_data_end_of_mission_cops/shape_models/SHAP5_PDS_V2/cg_spc_shap5_003m_cart.bds \
  --out-dir /data/vtenishe/rosetta_data_end_of_mission_cops/shape_models/SHAP5_PDS_V2 \
  --frame-kernel /home/vtenishe/SPICE/Kernels/ROSETTA/esa_2022/fk/ROS_V38.TF \
  --cspice-root /home/vtenishe/SPICE/cspice
```

The tool refuses existing BDF/manifest destinations. It requires exactly one
type-2 segment in `67P/C-G_CK`, uses 17-digit `dskexp` vertices, converts km to
m, and preserves the relative order of retained plates. The lower four PDS
reductions contain reversed exact duplicate triangle pairs forming
zero-thickness fins. Both copies of only those exact pairs are removed; any
same-winding duplicate, open edge, non-manifold remainder, or inconsistent
winding is fatal. Each `.bdf.conversion.json` records source/output hashes,
counts, scale, commands, PDS label, and the exact removed input face IDs.

Use the converted surface as follows:

```ini
ShapeModel = /data/vtenishe/rosetta_data_end_of_mission_cops/shape_models/SHAP5_PDS_V2/cg_spc_shap5_199k_cart.bdf
ShapeModelFormat = nastran
BodyFrame = 67P/C-G_CK
```

Choose a fresh `OutputDirectory` and rebuild all geometry-dependent tables and
kernels. Exact retained counts and validation results are in
[A7_PERFORMANCE.md](A7_PERFORMANCE.md).

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

### Stage-A6 implementation validation

Run the complete fast gate after changing the fitter, harmonic basis,
coefficient exchange, selection logic, manifests, or C++ source model:

```bash
python3 tools/run_a6_validation.py --threads 2
```

This discovers every Python test, including deterministic exact/noisy/
unobservable synthetic recovery, Jacobian, regularization, coefficient-file,
basis, selection, and reporting tests. It then runs the compiled C++ analytic
and transport self-test. `--threads` sets `OMP_NUM_THREADS` for C++ work;
BLAS and MKL are held to one thread to avoid multiplying worker counts.

Before accepting an A-stage change, run the external-data integration gate:

```bash
python3 tools/run_a6_validation.py --integration --threads 2 \
  --report a6_validation/validation_YYYYMMDD.json
```

Integration mode sets `COPS_RUN_INTEGRATION=1` and exercises small real C++
runs, a bitwise zero-enhancement regression, a Python degree-5 fit and its
generated C++ replay, the `1e-6` NG/RG comparison gate, plotting, profile
selection, boundary rebuilding, and configured mission/source inputs. Test
products are created in temporary directories and removed by their tests. An
optional report path must be new; a relative path is placed under
`<AMPS>/cops_output`, while an absolute path must remain outside this source
tree. The JSON records exact commands, environment, elapsed times, return
codes, output tails, and hashes of the executable and A6 test source.

Use `--skip-selftest` only when intentionally testing the Python portion; it
does not constitute the full gate. Passing these checks is implementation
consistency, not independent physics, calibration, localization, or
full-descent-fit validation. See [A6_VALIDATION.md](A6_VALIDATION.md) for the
requirement matrix, thresholds, and interpretation.

### Stage-A7 performance characterization

Before a multi-degree/multi-lambda campaign, measure the actual manifested
base on the target host.  Use absolute input paths if the retained base and HK
table are outside the configured output root; doing so does not change where
future relative outputs go:

```bash
env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  python3 tools/benchmark_fit_performance.py \
  /absolute/path/to/output_fit_full_default_base \
  --hk-table /absolute/path/to/COPS_HK_20160930.txt \
  --profile all_available \
  --degrees 5 6 7 8 \
  --repeat 3 \
  --jobs 1 2 4 8 14 \
  --threads-per-process 1 \
  --campaign-profiles 3 \
  --campaign-lambdas 5 \
  --campaign-max-nfev 1000 \
  --campaign-id a7measured \
  --json a7_performance/performance_YYYYMMDD.json \
  --markdown a7_performance/performance_YYYYMMDD.md
```

Relative report paths above are created under `<AMPS>/cops_output` (or
`COPS_OUTPUT_ROOT`), must be new, and do not change the default output root.
The four default output-size probes run a real production fit truncated at
one evaluation in temporary directories and delete it afterward; they are not
science fits.  Use `--skip-output-probe` to omit them.

The measured 100,000-facet base recommended
`--jobs 8 --threads-per-process 1` for at least eight independent lambdas.
Fourteen jobs were faster in aggregate but fell below the declared 60%
parallel-efficiency gate.  Repeat the measurement on another host or shape;
do not treat that recommendation as portable.  Full tables, metric meanings,
the campaign-estimate formula, higher-resolution shape options, and caveats
are in [A7_PERFORMANCE.md](A7_PERFORMANCE.md).

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
| `ActivityEnhancementModel` | `none` by default, or `real_spherical_harmonics_v1` to apply a strict H2O-only residual. |
| `ActivityEnhancementDegree` | Degree 4 for the disabled/zero-length control; 5--8 when enabled, exactly matching the file. |
| `ActivityEnhancementFile`, `ActivityEnhancementScale` | Required coefficient file and optional duplicated scale. The file supplies both degree and scale; any duplicate must agree exactly. |
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

### Change the surface boundary condition or shape model

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

Select a different surface mesh without editing `fit_base.in` by passing the
path and its actual reader format. This example creates the complete static
fitting base for the converted 199k SHAP5 mesh:

```bash
python3 tools/rebuild_cops_boundary.py --input fit_base.in \
  --surface-bc flat \
  --shape-model /data/vtenishe/rosetta_data_end_of_mission_cops/shape_models/SHAP5_PDS_V2/cg_spc_shap5_199k_cart.bdf \
  --shape-model-format nastran \
  --out output_fit_base_shap5_199k --threads 16 --dry-run

python3 tools/rebuild_cops_boundary.py --input fit_base.in \
  --surface-bc flat \
  --shape-model /data/vtenishe/rosetta_data_end_of_mission_cops/shape_models/SHAP5_PDS_V2/cg_spc_shap5_199k_cart.bdf \
  --shape-model-format nastran \
  --out output_fit_base_shap5_199k --threads 16
```

Changing the shape invalidates the old facet table, visibility, shadows,
altitudes, pressures, and every point kernel. The subsequent fitter receives
the shape through this new base; it must not reuse a base made for another
mesh.

| Option | Meaning |
|---|---|
| `--input FILE` | Source input; default `fit_base.in`. All its ordinary relative data paths are interpreted under the application root, not its own directory. |
| `--surface-bc flat\|rough` | Required target angular emission law. |
| `--out DIR` | Required new output directory; any existing path is refused. |
| `--threads N` | Optional positive OpenMP thread count; otherwise preserve the source input/environment. |
| `--shape-model FILE` | Optional surface-mesh override. Relative paths resolve from the application root; the generated input records an absolute path. |
| `--shape-model-format nastran\|dsk` | Optional reader override. No format is inferred from the filename; when omitted, preserve the source value or the solver's `nastran` default. |
| `--executable FILE` | Existing solver; default application-root `liouville_cops`. No compilation is performed. |
| `--dry-run` | Show configuration and command only. |

Relative `--out` resolves under the output root. Other relative CLI paths
resolve against the application root, with generated-input fallback. The tool
preserves temperatures, source parameters, and selection intervals, sets
`WriteKernels = on`, applies any explicit shape overrides, and saves normalized absolute paths in
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
| `ActivityEnhancementFacets.dat` | Enabled-only per-face legacy-linear, residual, total-linear, clipped activity, geometry, and clipping flag. |
| `input.snapshot.in`, `run_manifest.json` | Immutable effective input source and hashed provenance, including the enhancement configuration/file when enabled. |

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

### Versioned selection profiles and the 14-process runner

The three checked-in profiles use one shared implementation in
`tools/cops_selection.py`. The table shows the 100k reference counts and the
current 788k counts; differences come only from the last four points moving
from inside to outside the refined shape:

| Profile | 100k NG/RG | 788k NG/RG | NG treatment |
|---|---:|---:|---|
| `legacy_base_mask` | 839 / 1,146 | 839 / 1,150 | Apply the base gauge-exclusion intervals. |
| `all_available` | 930 / 1,146 | 934 / 1,150 | Retain all valid data, including saturated and post-transition NG. |
| `transition_cut_100m` | 869 / 1,146 | 873 / 1,150 | Exclude NG from `10:24:42` through `10:37:10`, inclusive. |

The audit fixes the 1,150 rows and the pre-inside policy counts, then derives
final counts from the active shape's `inside` flags. It does not relax audits
of pressures, times, intervals, or saturation policy.

Audit all three selections without optimizing:

```bash
python3 tools/run_fit_profile.py \
  --run output_fit_full_default_base \
  --all-profiles \
  --jobs 14 \
  --data-only \
  --out-prefix a1_selection
```

Run all three literal-default fits (degree 4, `lambda=0`, unordered/free
temperature, free `kappa`) with the requested process limit:

```bash
python3 tools/run_fit_profile.py \
  --run output_fit_full_default_base \
  --all-profiles \
  --degree 4 \
  --lambda 0 \
  --jobs 14 \
  --threads-per-process 1 \
  --tbounds 30 400 \
  --max-nfev 1000 \
  --ftol 1e-10 --xtol 1e-10 --gtol 1e-10 \
  --check-tolerance 1e-6 \
  --out-prefix a1_fit
```

The runner executes profiles sequentially to avoid holding three kernel sets
at once. `--jobs` is the number of independent lambda optimizations allowed in
the fitter pool; with the single literal-default `--lambda 0`, only one task
exists and therefore one process performs that optimization. To use as many as
14 fit processes, supply at least 14 distinct, scientifically chosen lambda
values. The runner sets `OMP_NUM_THREADS`, `OPENBLAS_NUM_THREADS`, and
`MKL_NUM_THREADS` to `--threads-per-process` (default 1) to avoid nested
oversubscription. Use `--dry-run` to print every complete command without
creating outputs.

For one profile, use `--profile NAME --out NEW_DIRECTORY`. A profile owns the
selection and cannot be combined with the fitter's legacy
`--apply-exclude-ng` or `--exclude-ng-saturated` switches.

### Controlled degree-5--8 activity-residual fit

The initial A4 experiment changes only the fine-scale H2O activity residual.
The fitter therefore requires `--fix-legacy-activity --fix-temperature
--fix-kappa` whenever `--enhancement-degree` is 5--8. It derives
`A_scale=sqrt(sum(area*A0_linear^2)/sum(area))`, initializes every `b_lm` to
zero, and gives the coefficients open bounds. `--enhancement-scale` may only
assert the derived value; it is not an override or free parameter.

Check all analytic degree-5 derivatives without optimizing:

```bash
python3 tools/run_fit_profile.py \
  --run output_fit_full_default_base \
  --profile all_available \
  --enhancement-degree 5 \
  --fix-legacy-activity --fix-temperature --fix-kappa \
  --lambda 1e-4 --check-jacobian \
  --jobs 14 --threads-per-process 1 \
  --out a4_jacobian_l5_all_available
```

Run a degree-8 regularization scan. Four lambdas mean at most four concurrent
optimizer processes even though the safety ceiling is 14:

```bash
python3 tools/run_fit_profile.py \
  --run output_fit_full_default_base \
  --profile all_available \
  --enhancement-degree 8 \
  --fix-legacy-activity --fix-temperature --fix-kappa \
  --lambda 1e-6 1e-5 1e-4 1e-3 \
  --jobs 14 --threads-per-process 1 \
  --tbounds 30 400 --max-nfev 1000 \
  --ftol 1e-10 --xtol 1e-10 --gtol 1e-10 \
  --check-tolerance 1e-6 \
  --out fit_enh_l8_all_available_scan
```

Each lambda directory contains `ActivityEnhancementFit.txt` and
`ActivityEnhancementFacets.dat` in addition to the normal fit products. It
also contains the universal A5 activity/residual tables and plots listed under
“Fit output locations” below; no separate plotting command is needed to make
those per-fit diagnostics.
`ActivityFit.txt` remains exactly three rows by 25 coefficients. Verify one
selected lambda using that subdirectory's generated input:

```bash
./liouville_cops \
  ../../cops_output/fit_enh_l8_all_available_scan/lambda=1.000e-04/cops_fit.in
python3 tools/validate_cops_run.py \
  --fit fit_enh_l8_all_available_scan/lambda=1.000e-04 \
  --verification fit_enh_l8_all_available_scan/lambda=1.000e-04_verify
```

Do not fit from a base that already contains an enabled enhancement. The
fitter rejects it because its cached H2O activity is no longer the legacy
baseline required by this parameterization. `lambda=0` is allowed for a rank
diagnostic and is explicitly labeled unregularized in the report/manifest.

### Late-NG visibility and effective-slope template fit

This optional controlled experiment adds two nonnegative H2O amplitudes to
the signed degree-5--8 residual: one for facets with an NG response but no RG
response during a declared late interval, and one for facets whose effective
gravitational slope is at least 60 degrees. The latter includes homogeneous
self-gravity and centrifugal acceleration; it is a categorical cliff
hypothesis, not an assumed continuous production law.

Install the pinned gravity package in an isolated environment and build the
SPICE helper. The helper uses the same C++ geometry/SPICE implementation as
the forward solver:

```bash
COPS_GRAVITY_VENV=/absolute/path/outside/LiouvilleCOPS/liouville-cops-gravity-venv
python3 -m venv "$COPS_GRAVITY_VENV"
"$COPS_GRAVITY_VENV/bin/pip" install -r requirements-gravity.txt
make -j16 liouville_cops spice_angular_velocity
```

Generate the 788k base once. `COPS_OUTPUT_ROOT` below is an explicit example;
choose a writable directory outside LiouvilleCOPS:

```bash
export COPS_OUTPUT_ROOT=/data/vtenishe/rosetta_data_end_of_mission_cops/cops_output

python3 tools/rebuild_cops_boundary.py \
  --input fit_base.in --surface-bc flat \
  --shape-model /data/vtenishe/rosetta_data_end_of_mission_cops/shape_models/SHAP5_PDS_V2/cg_spc_shap5_788k_cart.bdf \
  --shape-model-format nastran \
  --out output_fit_shap5_788k_base \
  --threads 10 --executable ./liouville_cops
```

Build the template bundle. The UTC bounds, density, gravity mesh, offset, and
validation tolerances are scientific inputs and are written to the metadata;
do not change them merely to improve the fit:

```bash
env OMP_NUM_THREADS=10 \
  "$COPS_GRAVITY_VENV/bin/python3" tools/build_physical_activity_templates.py \
  output_fit_shap5_788k_base \
  --out physical_templates_shap5_788k.npz \
  --late-start 2016-09-30T10:33:10 \
  --late-stop 2016-09-30T10:39:26 \
  --gravity-shape /data/vtenishe/rosetta_data_end_of_mission_cops/shape_models/SHAP5_PDS_V2/cg_spc_shap5_788k_cart.bdf \
  --density 533 \
  --cliff-threshold 45 --cliff-threshold 60 \
  --gravity-offset 0.1 \
  --gravity-offset-tolerance 0.25 \
  --gravity-resolution-tolerance 5 \
  --gravity-chunk 2048 \
  --omega-executable ./spice_angular_velocity
```

The command writes the immutable `.npz` mask bundle, a JSON audit record, a
body-frame PNG, and `physical_templates_shap5_788k.gravity.npz` under the
output root. The gravity cache preserves per-facet homogeneous self gravity,
the reference centrifugal/effective fields, and effective slope. It is bound
to exact centroid/normal/area ordering and to the gravity-shape hash, density,
offset, and algorithm. To build a different mask bundle on the same modeled
surface without repeating the expensive polyhedral calculation, add:

```bash
  --gravity-cache physical_templates_shap5_788k.gravity.npz
```

The cache loader recomputes the inexpensive centrifugal/effective fields when
the SPICE angular velocity differs and refuses any physical-input or mesh
mismatch. Existing NG-only/cliff masks are already reusable by passing their
bundle directly to `--activity-template`; rebuilding is unnecessary. The main
comparison fits the 60-degree
field; retain the 45-degree field for a separate sensitivity run rather than
fitting two nested masks simultaneously. Run the same positive-lambda grid on
the full profile and the profile that removes the NG transition:

```bash
python3 tools/run_fit_profile.py \
  --run output_fit_shap5_788k_base \
  --profile all_available \
  --enhancement-degree 8 \
  --fix-legacy-activity --fix-temperature --fix-kappa \
  --activity-template physical_templates_shap5_788k.npz late_ng_not_rg \
  --activity-template physical_templates_shap5_788k.npz cliff_ge_60deg \
  --lambda 1e-6 3.162e-6 1e-5 3.162e-5 1e-4 3.162e-4 1e-3 \
  --jobs 2 --threads-per-process 1 \
  --max-nfev 1000 --ftol 1e-10 --xtol 1e-10 --gtol 1e-10 \
  --check-tolerance 1e-6 \
  --out fit_shap5_788k_l8_templates_all_available_scan

python3 tools/run_fit_profile.py \
  --run output_fit_shap5_788k_base \
  --profile transition_cut_100m \
  --enhancement-degree 8 \
  --fix-legacy-activity --fix-temperature --fix-kappa \
  --activity-template physical_templates_shap5_788k.npz late_ng_not_rg \
  --activity-template physical_templates_shap5_788k.npz cliff_ge_60deg \
  --lambda 1e-6 3.162e-6 1e-5 3.162e-5 1e-4 3.162e-4 1e-3 \
  --jobs 2 --threads-per-process 1 \
  --max-nfev 1000 --ftol 1e-10 --xtol 1e-10 --gtol 1e-10 \
  --check-tolerance 1e-6 \
  --out fit_shap5_788k_l8_templates_transition_cut_scan
```

Seven lambda values provide seven independent optimizer tasks. On the measured
788k campaign, each prepared process occupied roughly 66--74 GiB, so only two
are run simultaneously on the approximately 242-GiB host. `--jobs` is a
ceiling, not a request to split one optimizer; remeasure memory before raising
it on another system. Compare these results with a
same-degree harmonic-only scan. A template coefficient at its lower bound has
one-sided inference; its local covariance entry must not be interpreted as a
two-sided detection uncertainty. The exact selected mixed field is saved in
`ActivityEnhancementFit.txt` using the mesh-locked per-facet replay format.

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
| `--out DIR` | `<run-name>_fit` | New results directory under the output root; any existing path is refused. |
| `--degree N` | 4 | Fitted harmonic groups through degree N, intended range 0–4; higher coefficients stay at base values. |
| `--enhancement-degree N` | 4 | Maximum degree of the separate H2O residual: 4 is the no-residual control; 5–8 fit degrees 5 through N. |
| `--fix-legacy-activity` | off | Freeze all 25 legacy H2O coefficients. It is required for residual fits and currently also requires fixed temperature and `kappa`. |
| `--enhancement-scale VALUE` | derived | Assert, but do not replace, the area-weighted RMS scale derived from the base legacy linear activity. |
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
| `--profile NAME` | off | Apply a checked-in versioned selection profile; mutually exclusive with the two legacy NG selection switches. |
| `--data-only` | off | Export selection plots/table after consistency checks, then stop. |
| `--allow-legacy-base` | off | Explicitly accept an unmanifested pre-W01 base after strict structural validation; prefer a regenerated manifested base. |

### Fit output locations

| Location | Contents |
|---|---|
| `<out>/fit_data_used.txt`, `fit_data_used.png`, `fit_profiles_used.png` | Actual selection exports, including per-gauge reason codes in the table and data-only/derivative-check modes. |
| `<out>/selection_profile.json` | Effective named profile, source identity, counts, and reason hash; present only with `--profile`. |
| Single fit: `<out>/fit_result.txt` | Parameters, termination and misfit/regularization diagnostics. |
| Single fit: `<out>/ActivityFit.txt` | Three species rows of 25 activity coefficients. |
| Single fit: `<out>/cops_fit.in` | Generated input for a full C++ verification run. |
| Single fit: `<out>/COPS_model.dat` | Fitted total NG/RG predictions in the base table layout. |
| Single fit: `<out>/fit_manifest.json` | Versioned hashes, cache/selection/prior/settings identities, optimizer result, and independent verification/review states. |
| Every fit result: `fit_report.json` | Canonical identity binding profile/selection, legacy and fitted activity, basis/degree, regularizer, and lambda. |
| Every fit result: `activity_facets.dat` | Centroid, area, selected-kernel contribution flag, legacy/residual/linear-total/clipped activity for every facet. |
| Every fit result: `activity_degree_power.dat`, `activity_region_diagnostics.dat` | Residual power by degree and area-weighted all/seen/unseen diagnostics, including clipped area. |
| Every fit result: `activity_maps.png` | Baseline, residual, and clipped-total activity maps. |
| Every fit result: `pressure_residuals.dat` | Observed/model pressure, selection/saturation/reason fields, and selected-only log residual for both gauges. |
| Every fit result: `pressure_comparison_time.png`, `pressure_comparison_altitude.png`, `pressure_residuals.png` | Selected/excluded/saturated/transition-aware pressure and residual figures. |
| Scan: `<out>/lambda=<value>/` | The four optimization products above for each lambda; directory values use `%.3e` formatting. |
| Scan: `<out>/lcurve.txt` | Aggregate scan diagnostics; selection exports remain in `<out>`. |

The Python-fitted table is not a complete forward-model recomputation:
aggregate NG/RG pressures are recomputed, while production, species
pressure/density/flux, and below-deck diagnostic columns are explicitly
`nan`. They are not inherited fitted quantities. Use the C++ verification
output for complete species results and production rates.

All A5 tables and PNG metadata carry the same `ReportIdentity`. The matching
optimizer record and every output hash are in `fit_manifest.json`. “Seen” in
the activity tables means nonzero contribution to at least one kernel loaded
for a selected NG or RG measurement; it does not by itself establish that the
facet is uniquely constrained.

## 10. Verify fitted parameters with C++ and test flight-time rotation

For the newly generated single-fit example, inspect `output_fit_new/cops_fit.in`
and its paths, then run:

```bash
env OMP_NUM_THREADS=4 ./liouville_cops output_fit_new/cops_fit.in
```

The generated input normally sets `ActivityModel = file`, points to that
fit's `ActivityFit.txt`, inserts the fitted temperature parameters, and sets
the destination to `output_fit_new_verify`. Check it before running; the
solver refuses an existing verification directory. For a scan result, use the
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
  --ref output_fit_new --save COPS_fit_comparison_new.png
```

Plot against altitude instead:

```bash
python3 tools/plot_cops_model_altitude.py \
  output_fit_base_new output_fit_new output_fit_new_verify \
  --ref output_fit_new --save COPS_fit_altitude_new.png
```

| Model-plotter argument | Meaning |
|---|---|
| Positional `runs ...` | One or more run directories or explicit `COPS_model.dat` files. Omit to scan output-root `output*/COPS_model.dat`, nonrecursively; if none exist, scan legacy application-root runs. |
| `--x time` / `--x altitude` | Horizontal axis; defaults to time in the main script and altitude in its wrapper. |
| `--save FILE` | Save the figure instead of opening a window. |
| `--theme light` / `--theme dark` | Plot style; default light. |
| `--ref NAME` | Loaded run's directory basename, or its suffix after `output_`, supplying observations, `useNG`/`useRG` marker masks, and reference geometry. Default is newest loaded table with useNG/useRG flags. Select the Python fit to show the samples actually used by that fit. |

Set `--ref` explicitly for reproducibility; it is not a path to an unloaded
run. The four panels show combined NG/RG, NG, RG, and their ratio; unused
reference points are hollow. For Python fit tables the reader verifies the
table/report hashes and profile, selection, report, and harmonic-basis
identities against `fit_manifest.json`. Labels show profile, enhancement
degree, lambda, and short report ID; an incompatible fit is rejected and an
unmanifested historical table is explicitly labelled legacy. Manifested C++
runs are also checked against their table hash and labelled with basis, degree,
and a short coefficient/activity identity. Nested scan directories must be passed
explicitly, for example:

```bash
python3 tools/plot_cops_model.py output_fit_base_new \
  output_lcurve_new/lambda=3.000e-03 \
  --ref lambda=3.000e-03 --save COPS_scan_fit_new.png
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
