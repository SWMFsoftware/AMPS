# Scripts

This file describes the scripts in the standalone LiouvilleCOPS application.
See [README.md](README.md) for details and [RUN.md](RUN.md) for build/path setup.
Run examples from the LiouvilleCOPS root. Use Python 3 with NumPy/Matplotlib
for plotting and SciPy for fitting/digitization; the housekeeping and boundary
rebuild tools use only the standard library. The rebuild tool requires Python
3.8+ and an already compiled C++ solver. All Python entry points support `--help`.

All new scientific products default to **`<AMPS>/cops_output/`**. Relative output
names below are under that root, not LiouvilleCOPS. Absolute output paths are
honored only outside LiouvilleCOPS; paths resolving inside it are rejected.
`COPS_OUTPUT_ROOT=/absolute/path/outside/LiouvilleCOPS` overrides the root for Python and
the rebuilt C++ executable. Generated-product readers prefer the new root,
then existing local files; explicit absolute or `./` inputs select that location.
Former local products have been removed. Source images/L4 data stay in the application.
See [RUN.md](RUN.md#output-location). `tools/cops_paths.py` is the shared path
helper, not a standalone command.

| Script | Purpose |
|---|---|
| [`tools/cops_hk_table.py`](#cops_hk_tablepy) | make the COPS housekeeping table from the archived L2 files |
| [`tools/rebuild_cops_boundary.py`](#rebuild_cops_boundarypy) | regenerate a static fitting base with flat or rough surface emission, preserving the old run |
| [`tools/plot_cops_model.py`](#plot_cops_modelpy-plot_cops_model_altitudepy) | compare the LiouvilleCOPS runs with the data, vs. time |
| [`tools/plot_cops_model_altitude.py`](#plot_cops_modelpy-plot_cops_model_altitudepy) | the same, vs. altitude |
| [`tools/fit_cops.py`](#fit_copspy) | fit the source rate and temperature to the COPS data |
| [`tools/plot_lcurve.py`](#plot_lcurvepy) | plot the L-curve of the regularized fit |
| [`tools/digitize_cops_figures.py`](#digitize_cops_figurespy) | digitize the published end-of-mission COPS figures |
| [`tools/collect_cops_sources.py`](#collect_cops_sourcespy) | gather the COPS data of all sources into one JSON file and two plots; write and plot the data used in the fit |

Typical order:

```
python3 tools/cops_hk_table.py /data/vtenishe/rosetta_data_end_of_mission_cops/data_COPS_L2/mtp35_cops_bg data/COPS_HK_20160930.txt
./liouville_cops cops.in
python3 tools/plot_cops_model.py
./liouville_cops fit_base.in
./liouville_cops m004.in
python3 tools/fit_cops.py output_fit_base --out output_fit_base_lcurve --lambda 1e-4 1e-3 1e-2 1e-1 1 --jobs 5
python3 tools/plot_lcurve.py output_fit_base_lcurve
python3 tools/digitize_cops_figures.py --check
python3 tools/collect_cops_sources.py
python3 tools/fit_cops.py output_fit_base --data data/COPS_fit_data.txt --lambda 3e-3
```

Units: the archive, the srcCG table and the published figures give pressures in mbar (N2-equivalent, 20 C).
LiouvilleCOPS and all plots of `tools/` use Pa (1 mbar = 100 Pa).

---

## `rebuild_cops_boundary.py`

Regenerates the angular transport kernels and forward pressures when changing
`SurfaceBC`. It runs the existing C++ solver, preserving source temperature,
activity, geometry, and selection settings. It does not run a fit. Requires
`Mode = observations`, `FlightTimeRotation = none`,
`TemperatureModel = analytical`, and `SurfaceTemperatureCorrectionFactor = 1`.

Preview, then execute:

```bash
python3 tools/rebuild_cops_boundary.py --input fit_base.in \
  --surface-bc rough --out output_fit_base_rough --threads 4 --dry-run

python3 tools/rebuild_cops_boundary.py --input fit_base.in \
  --surface-bc rough --out output_fit_base_rough --threads 4
```

Options: `--input FILE` (default `fit_base.in`), required
`--surface-bc flat|rough`, required `--out DIR`, optional `--threads N`,
`--executable FILE` (default `liouville_cops`), and `--dry-run` (no writes or
execution; external files are not checked). Relative `--out` uses the shared
output root. Other relative source paths use the application root, with
generated input/HK/activity lookup; kernel entries use `KernelDirectory`.

The output directory must not exist. It receives `boundary.in`, `source.in`,
`COPS_model.dat`, `elements.bin`, `kernels/`, `rebuild.log`, and `rebuild.json`.
Successful completion includes table/cache completeness checks. Failed or
interrupted output is retained for diagnosis and must not be fitted.

After success, use the new base as the fitter's positional argument, for example:

```bash
python3 tools/fit_cops.py output_fit_base_rough \
  --apply-exclude-ng --exclude-ng-saturated --ordered-temperature \
  --lambda 3e-3 --out output_rough_fit
```

The working fitter is included. The lambda is an
example, not a prescribed optimum. Static kernels depend on the boundary
condition but not surface temperature; the fitter reevaluates thermal speed
factors during optimization. See README section 6.6 for the detailed dependency
table, temperature equations, safeguards, and tests.

## `cops_hk_table.py`

* **Purpose:** make one table of the COPS pressures and housekeeping of 2016-09-30.
* **Input:** the ROSINA COPS level 2 BG files (`bg_*.tab`, ESA PSA data set RO-C-ROSINA-2-EXT3-V2.0, MTP35). The
  local copy is in `/data/vtenishe/rosetta_data_end_of_mission_cops/data_COPS_L2/mtp35_cops_bg`.
* **Output:** `data/COPS_HK_20160930.txt`, one row per file, time = middle of the file interval. Columns: NG and RG
  pressure (mbar), sensor temperature (T_RG) and electronics temperatures, emission currents and ranges, autorange
  flags, data quality. LiouvilleCOPS reads T_RG from it (`TemperatureRGMode = hk_table`).
* **Arguments:** optional raw BG directory, default
  `<COPS_DATA_ROOT>/data_COPS_L2/mtp35_cops_bg`; optional output filename,
  default `data/COPS_HK_20160930.txt` under the output root. `--help` prints usage.
* **Example:**

  ```
  python3 tools/cops_hk_table.py /data/vtenishe/rosetta_data_end_of_mission_cops/data_COPS_L2/mtp35_cops_bg data/COPS_HK_20160930.txt
  ```

## `plot_cops_model.py`, `plot_cops_model_altitude.py`

* **Purpose:** compare the model pressures of one or several LiouvilleCOPS runs (or fits) with the COPS data.
  `plot_cops_model_altitude.py` is the same script with altitude as the default horizontal axis.
* **Input:** `COPS_model.dat` of the runs (`output*/`, or the `<run>_fit/` and `lambda=*/` fit directories). The
  reference data in it are the srcCG table, identical to the archived L2 values.
* **Output:** a window with four panels (NG and RG together, NG, RG, NG/RG). Check boxes switch each run on and off;
  the points not used in the fit are hollow. With `--save`, a PNG file instead.
* **Options:** run directories or files (default: all `output*/`), `--x time|altitude`, `--ref RUN` (the run that
  provides the reference data and their altitudes), `--save FILE`, `--theme light|dark`.
* **Examples:**

  ```
  python3 tools/plot_cops_model.py
  python3 tools/plot_cops_model.py output_fit_base output_fit_base_fit output_rough
  python3 tools/plot_cops_model_altitude.py output_M004 --save m004_altitude.png
  ```

## `fit_cops.py`

* **Purpose:** fit the H2O source map (spherical harmonics up to degree 4), the analytical surface temperature
  (Tmax, T75, Tmin) and the CO2 scale kappa to the COPS NG and RG pressures. The misfit is logarithmic, with equal
  weight per gauge. A Tikhonov term `--lambda` pulls the H2O map toward the srcCG map; several values give an L-curve.
* **Input:** a LiouvilleCOPS run made with `WriteKernels = on` and `FlightTimeRotation = none` (e.g. `fit_base.in`
  -> `output_fit_base`): its geometric kernels `kernels/point=NNNN.bin`, `elements.bin`, `COPS_model.dat` (data, fit
  selection `useNG`/`useRG`, T_RG) and the srcCG coefficients of its `cops.in` (prior). With `--data FILE` the
  pressures and the `useNG`/`useRG` flags come from a table instead (`data/COPS_fit_data.txt` of
  `collect_cops_sources.py`, which can be edited); its points must be points of the run, with the same T_RG.
* **NG data after the emission range change:** the `ExcludeNG` interval of the base run (`fit_base.in`:
  10:24:40-10:40:00, NG filament emission range HIGH->LOW until the impact) is **not applied by default**: its 91 NG
  points with data are used (930 NG points in total). `--apply-exclude-ng` blocks them (839 NG points, as in the fits
  in `output_fit_base_fit` and `output_fit_base_lcurve`).
* **Saturated NG points:** the 13 points 10:23:50-10:24:38, where the NG reading is held at 2.2633e-8 mbar (L2 BG
  records 10:23:46-10:24:38, filament emission current 0), are **used by default**. `--exclude-ng-saturated` removes
  them (identified in the housekeeping table, `--hk-table`, default `data/COPS_HK_20160930.txt`); this works with
  and without `--data`. `fit_result.txt` records the choice.
* **Output:**
  * always, before the fit: `<out>/fit_data_used.txt` and `<out>/fit_data_used.png`, the data that actually enter
    the misfit (after the selection of the run or `--data` and `--exclude-ng-saturated`), taken from the fit problem
    itself: used points filled, other data open grey; and the separate plot `<out>/fit_profiles_used.png` with only
    the used NG and RG profiles. `--data-only` writes them and stops;
  * one lambda: `<run-name>_fit/` under the output root (or `--out`) with `fit_result.txt` (parameters, uncertainties, rms, p_eff, activity
    diagnostics), `ActivityFit.txt` (coefficients), `COPS_model.dat` (model at the fitted parameters) and
    `cops_fit.in` (LiouvilleCOPS input that reruns the fitted model into `<fit dir>_verify`);
  * several lambdas: `<out>/lambda=<value>/` with the same files, and `<out>/lcurve.txt` with the maximum-curvature
    lambda.
* **Options:** `--out DIR`, `--degree N`, `--lambda L [L ...]`, `--jobs N`, `--fix-temperature`,
  `--ordered-temperature`, `--fix-kappa`, `--tbounds TMIN TMAX` (default 30 400 K), `--max-nfev`, `--ftol`, `--xtol`,
  `--gtol`, `--check-jacobian`, `--check-tolerance`, `--data FILE`, `--apply-exclude-ng`, `--exclude-ng-saturated`, `--hk-table FILE`, `--data-only`.
* **Examples:**

  ```
  python3 tools/fit_cops.py output_fit_base --data-only                  # plot the data the fit would use
  python3 tools/fit_cops.py output_fit_base --data-only --exclude-ng-saturated --out output_fit_base_nosat
  python3 tools/collect_cops_sources.py --fit-used output_fit_base_nosat/fit_data_used.txt
  python3 tools/fit_cops.py output_fit_base --check-jacobian
  python3 tools/fit_cops.py output_fit_base --lambda 3e-3
  python3 tools/fit_cops.py output_fit_base --out output_fit_base_lcurve --jobs 8 \
          --lambda 1e-5 1e-4 1e-3 3e-3 1e-2 1e-1 1 10
  python3 tools/fit_cops.py output_fit_base --lambda 3e-3 --fix-temperature --out output_fit_base_fixT
  python3 tools/fit_cops.py output_fit_base --lambda 3e-3 --data data/COPS_fit_data.txt --out output_fit_base_table
  python3 tools/fit_cops.py output_fit_base --lambda 3e-3 --exclude-ng-saturated --out output_fit_base_nosat
  python3 tools/fit_cops.py output_fit_base --lambda 3e-3 --apply-exclude-ng --out output_fit_base_blockNG
  ./liouville_cops output_fit_base_fit/cops_fit.in    # verify a fit with the C++ code
  ```

## `plot_lcurve.py`

* **Purpose:** show the L-curve and how the fitted parameters change with lambda.
* **Input:** `lcurve.txt` written by `fit_cops.py` (the directory or the file).
* **Output:** three panels: the L-curve (misfit vs. deviation from srcCG, maximum curvature marked), Tmax/T75/Tmin
  vs. lambda, and kappa, activity diagnostics and rms vs. lambda. Open circles are fits that reached `--max-nfev`.
  A window, or a PNG with `--save`.
* **Options:** `--save FILE`.
* **Examples:**

  ```
  python3 tools/plot_lcurve.py output_fit_base_lcurve
  python3 tools/plot_lcurve.py output_fit_base_lcurve/lcurve.txt --save output_fit_base_lcurve/lcurve.png
  ```

## `digitize_cops_figures.py`

* **Purpose:** extract the NG and RG curves from the published figures that show end-of-mission COPS data not in
  the archive.
* **Input** (in `<COPS_DATA_ROOT>/sources/` by default):
  * `EGU2019_Tenishev_slide02.jpg`, `EGU2019_Tenishev_slide16_data.jpg`: Tenishev et al., EGU 2019, slides 2 and
    16 (from `~/Downloads/2019-EGU.pptx`), NG and RG vs. altitude;
  * `ESA_blog_2016-12-15_COPS_EoM.jpg`: ESA Rosetta blog 2016-12-15 (credit K. Altwegg), NG vs. UTC.
* **Output:** `<AMPS>/cops_output/data/sources/digitized_<figure>_<NG|RG>.csv` (altitude_m or utc, pressure_mbar); with `--check`
  also `digitized_check.png`, the digitized points drawn over the images, and
  `EGU2019_slide02_on_slide16_axes.png`, slide 2 resampled onto the axes of slide 16 (the two slides are different
  figures: equal above ~4 km, different below ~3 km).
* **Options:** `--check`, `--dir DIR` (input images, default
  `<COPS_DATA_ROOT>/sources`), `--out-dir DIR` (products, default output-root `data/sources`).
* **Example:**

  ```
  python3 tools/digitize_cops_figures.py --check
  ```

## `collect_cops_sources.py`

* **Purpose:** gather the end-of-mission COPS NG and RG pressures of all sources in the same units: one
  machine-readable JSON file and two plots. Also write the pressures that are actually used in the fit as a table
  for `fit_cops.py --data`, and plot them.
* **Input:**
  * `data/COPS_HK_20160930.txt` (archive L2, from `cops_hk_table.py`);
  * `COPS_model.dat` of the SHAP5 run (the srcCG table);
  * `<COPS_DATA_ROOT>/sources/COPS_L4_MTP35.ASC` (archive L4 NG densities, ESA PSA RO-C-ROSINA-4-EXT3-V2.0);
  * output-root `data/sources/digitized_*.csv` (from `digitize_cops_figures.py`);
  * two LiouvilleCOPS runs for the altitudes above SHAP5 and M004;
  * the base run of the fit (`--fit-run`, default `output_fit_base`): its `useNG`/`useRG` flags and its input file
    (`ExcludeNG`, `ExcludeRG`, `FitInterval`) for the reason of each exclusion.
* **Output:**
  * `data/COPS_EoM_sources.json`: datasets `PSA_L2_BG`, `srcCG_table`, `PSA_L4`, `ESA_blog_2016`, `EGU2019_s02`,
    `EGU2019_s16`. Each has source, notes and columns: UTC, altitude above SHAP5 and M004, distance from the centre,
    pressures (Pa), T_RG, NG densities and RG fluxes (N2-equivalent and H2O), flags. Also unit conventions, events
    (NG saturation, NG gap, touchdown) and the sources that could not be used. Missing values are `null`.
  * `data/COPS_EoM_sources_time.png`: NG and RG vs. UTC (whole descent and last 40 min).
  * `data/COPS_EoM_sources_altitude.png`: NG and RG vs. altitude.
  * In both plots the points used in the fit are small red markers; the JSON has them as dataset `fit_data`.
  * `data/COPS_fit_data.txt`: the fit data, one row per point of the fit run: `iPoint UTC altitude_SHAP5_m pNG_Pa
    pRG_Pa useNG useRG T_RG_K ng_saturated reason_NG reason_RG` (Pa, 0 = no data). Edit `useNG`/`useRG` (0/1) to
    change the selection and pass the file to `fit_cops.py --data`. `ng_saturated` = 1 flags the 13 saturated NG
    points (10:23:50-10:24:38); they are included by default (`useNG` = 1), `--exclude-ng-saturated` writes `useNG` = 0
    for them.
  * `data/COPS_fit_data.png`: the fit data vs. UTC and vs. SHAP5 altitude; used points filled, excluded points open
    and coloured by reason; the Exclude intervals and the saturated NG points shaded.
  * `data/COPS_fit_profiles.png`: a separate plot of only the NG and RG profiles actually used in a fit, read from
    the `fit_data_used.txt` of that fit (`--fit-used`, default `output_fit_base_selection/fit_data_used.txt`, the
    default selection; run `fit_cops.py --data-only` first for a new selection).
* **Options:** `--shap5-run DIR` (default `output_fit_base`), `--m004-run DIR` (default `output_M004`),
  `--fit-run DIR` (default `output_fit_base`), `--sources-dir DIR` (default `<COPS_DATA_ROOT>/sources`), `--fit-data FILE` (default `data/COPS_fit_data.txt`),
  `--apply-exclude-ng` (default: the NG data in `ExcludeNG` are not blocked), `--exclude-ng-saturated` (default: the
  saturated NG points are included), `--fit-used FILE` (default `output_fit_base_selection/fit_data_used.txt`),
  `--out FILE.json`, `--plots PREFIX`, `--show`.
* **Clean-checkout prerequisites:** run `cops_hk_table.py`,
  `digitize_cops_figures.py`, `liouville_cops fit_base.in`, and
  `liouville_cops m004.in` first. The supplied `m004.in` writes
  `output_M004/COPS_model.dat` with exact M004/SPICE geometry but zero model
  pressures (`GeometryOnly = on`) and without fitting kernels. It is a geometry
  input for the collector, not a forward-pressure run or fitting base. The collector checks
  all required inputs up front and lists every missing file; it does not
  regenerate these potentially expensive products implicitly.
* **Examples:**

  ```
  python3 tools/collect_cops_sources.py
  python3 tools/collect_cops_sources.py --show
  python3 tools/collect_cops_sources.py --out /tmp/cops.json --plots /tmp/cops
  python3 tools/collect_cops_sources.py --fit-run <other base run> --fit-data data/COPS_fit_data_<name>.txt
  ```

  Excluding the saturated NG points from a fit (default: included), either in the fit or in the table:

  ```
  python3 tools/fit_cops.py output_fit_base --lambda 3e-3 --exclude-ng-saturated
  python3 tools/collect_cops_sources.py --exclude-ng-saturated --fit-data data/COPS_fit_data_nosat.txt
  python3 tools/fit_cops.py output_fit_base --data data/COPS_fit_data_nosat.txt --lambda 3e-3
  ```

  Reading the JSON file:

  ```python
  import json
  d = json.load(open("data/COPS_EoM_sources.json"))
  l2 = next(s for s in d["datasets"] if s["id"] == "PSA_L2_BG")
  utc, pNG, pRG = l2["columns"]["utc"], l2["columns"]["pNG_Pa"], l2["columns"]["pRG_Pa"]
  ```
