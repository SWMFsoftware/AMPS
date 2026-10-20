# Scripts

Restored-fit development: [RESTORED_BASELINE.md](RESTORED_BASELINE.md) documents
the campaign preparation tool; [JOINT_TEMPLATE_FIT.md](JOINT_TEMPLATE_FIT.md)
documents simultaneous baseline/temperature/CO2/template fitting.
[GRAVITY_REPAIR.md](GRAVITY_REPAIR.md) documents the repaired builder and cache
identity. The fitter refuses old cliff columns without repaired orientation
provenance even if the old bundle is renamed; its geometry-only NG column is
unaffected by the HEAL defect.

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
Former local products have been removed. Source images, L4 data, and the SHAP5
mesh are external under `/data/vtenishe/rosetta_data_end_of_mission_cops` by default.
See [RUN.md](RUN.md#output-location). `tools/cops_paths.py` is the shared path
helper, not a standalone command.

| Script | Purpose |
|---|---|
| [`tools/cops_hk_table.py`](#cops_hk_tablepy) | make the COPS housekeeping table from the archived L2 files |
| [`tools/rebuild_cops_boundary.py`](#rebuild_cops_boundarypy) | regenerate a static fitting base with a selected shape and flat or rough surface emission, preserving the old run |
| [`tools/plot_cops_model.py`](#plot_cops_modelpy-plot_cops_model_altitudepy) | compare the LiouvilleCOPS runs with the data, vs. time |
| [`tools/plot_cops_model_altitude.py`](#plot_cops_modelpy-plot_cops_model_altitudepy) | the same, vs. altitude |
| [`tools/fit_cops.py`](#fit_copspy) | fit the source rate and temperature to the COPS data |
| [`tools/run_fit_profile.py`](#run_fit_profilepy) | run one or all audited selection profiles safely |
| [`tools/run_a6_validation.py`](#run_a6_validationpy) | run the Stage-A6 unit, synthetic, C++, and optional integration gate |
| [`tools/benchmark_fit_performance.py`](#benchmark_fit_performancepy) | measure Stage-A7 production-fitter time, memory, output size, and independent-job scaling |
| [`tools/build_physical_activity_templates.py`](#build_physical_activity_templatespy) | build mesh-bound late-NG-only and effective-slope cliff hypotheses |
| [`tools/convert_dsk_to_bdf.py`](#convert_dsk_to_bdfpy) | convert archived type-2 DSK shapes to validated metre-unit NASTRAN surface BDFs |
| [`tools/convert_ply_to_bdf.py`](PLY_CONVERSION.md) | convert binary triangular PLY to metre-unit GRID*/CTRIA3 BDF, preserving geometry and recording the gravity topology audit |
| [`tools/repair_collinear_ply_bdf.py`](PLY_CONVERSION.md#separate-surface-preserving-44m-repair) | repair isolated exactly collinear source triangles by checked edge flips; preserve all BDF vertices and originals, emit a separate mesh and per-face provenance; production validation is still required |
| [`tools/benchmark_activity_basis.py`](#benchmark_activity_basispy) | benchmark the isolated degree-5--8 real-harmonic table on an element cache |
| [`tools/validate_cops_run.py`](#validate_cops_runpy) | validate a fitting cache or register C++ reproduction of fitted pressures |
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

## `benchmark_activity_basis.py`

* **Purpose:** validate a version-1 `elements.bin`, memory-map its facet
  centroids, and measure construction of the independently normalized Stage-A2
  real-harmonic table. It never modifies the cache or runs the source model.
* **Arguments:** optional element file (default
  `output_fit_full_default_base/elements.bin`), `--min-degree` (default 5),
  `--max-degree` (default 8), `--repeat` (default 3), and optional
  `--json FILE` under the configured output root.
* **Output:** JSON on stdout with basis ID, cache header, table shape/exact
  bytes, timing samples/median, maximum-RSS delta, checksum, and maximum
  absolute basis value.
* **Example:**

  ```bash
  python3 tools/benchmark_activity_basis.py \
    output_fit_full_default_base/elements.bin \
    --min-degree 5 --max-degree 8 --repeat 5
  ```

The mathematical source of truth is `HARMONIC_BASIS.md`; Python is implemented
in `tools/activity_basis.py` and C++ in
`src/ActivityEnhancementBasis.{h,cpp}`. Stage A3 makes it an optional C++
H2O forward-model component using the strict format in
`ACTIVITY_ENHANCEMENT_FILE.md`; Stage A4 estimates its coefficients and Stage
A5 reports their facet/degree/pressure diagnostics.

---

## `build_physical_activity_templates.py`

* **Purpose:** make auditable, mesh-bound H2O predictor fields without fitting
  their amplitudes. `late_ng_not_rg` is the union, over an explicit UTC
  interval, of facets with `Kn > 0` and simultaneous `KPhi = KV = 0`.
  `cliff_ge_<angle>deg` uses the angle between the inward facet normal and
  homogeneous-polyhedron self-gravity plus centrifugal acceleration.
* **Inputs:** a manifested static fitting base; explicit late start/stop; a
  closed NASTRAN gravity surface and positive homogeneous density; the compiled
  `spice_angular_velocity` helper; repeatable slope thresholds; and optionally
  a second gravity surface for a sampled resolution check.
* **Safeguards:** the tool records exact input hashes, evaluates the SPICE
  angular velocity at three epochs, checks gravity-offset convergence, can
  compare gravity-shape resolutions, and binds the output to the exact ordered
  centroid/normal/area table. Response magnitudes are diagnostics and do not
  set the NG-only mask.
* **Outputs:** a new `.npz` mask bundle, companion JSON metadata, PNG map, and
  `<bundle-stem>.gravity.npz` under the output root. The mask bundle contains
  only finite, nonnegative template columns and cannot be loaded with another
  mesh/order. The gravity cache preserves per-facet self gravity and reference
  centrifugal/effective-gravity and slope arrays.
* **Reuse:** an existing mask bundle is passed directly to the fitter. For a
  different set of masks on the same exact modeled surface, pass the existing
  gravity file as `--gravity-cache FILE`. The loader requires matching
  centroid/normal/area order, gravity-shape hash, density, offset, and
  algorithm. It can reuse self gravity and cheaply re-derive the
  SPICE-dependent fields for a different angular velocity.
* **Dependency:** install the pinned `polyhedral-gravity` package from
  `requirements-gravity.txt`; `make spice_angular_velocity` builds the helper.
* **Use in a fit:** repeat `fit_cops.py --activity-template BUNDLE NAME` (or
  the equivalent profile-runner option). Every template gets one nonnegative
  amplitude; the controlled fit must also fix legacy activity, temperature,
  and `kappa`. Do not fit the nested 45- and 60-degree masks together.

```bash
COPS_GRAVITY_VENV=/absolute/path/outside/LiouvilleCOPS/liouville-cops-gravity-venv
env OMP_NUM_THREADS=10 \
  "$COPS_GRAVITY_VENV/bin/python3" tools/build_physical_activity_templates.py \
  output_fit_shap5_788k_base \
  --out physical_templates_shap5_788k.npz \
  --late-start 2016-09-30T10:33:10 --late-stop 2016-09-30T10:39:26 \
  --gravity-shape /data/vtenishe/rosetta_data_end_of_mission_cops/shape_models/SHAP5_PDS_V2/cg_spc_shap5_788k_cart.bdf \
  --density 533 --cliff-threshold 45 --cliff-threshold 60 \
  --gravity-offset 0.1 --gravity-offset-tolerance 0.25 \
  --gravity-resolution-tolerance 5 --gravity-chunk 2048 \
  --omega-executable ./spice_angular_velocity
```

`tools/activity_templates.py` is the imported strict bundle I/O module, not a
standalone program. `tools/gravity_cache.py` is the corresponding strict cache
I/O module and is not a command. `tools/spice_angular_velocity.cpp` is also not
a Python entry point; run its compiled target only when diagnosing the SPICE
rate.

---

## `convert_dsk_to_bdf.py`

* **Purpose:** convert an archived one-segment SPICE type-2 DSK to the
  surface-only NASTRAN format consumed from the original SHAP5 BDF. Coordinates
  are converted from kilometres to metres and written as `GRID*`; retained
  triangular plates are written as `CTRIA3` in their original relative order.
* **Topology policy:** reversed exact duplicate triangle pairs in the lower
  PDS SHAP5 reductions are zero-thickness fins. The tool removes both copies
  of only those exact pairs and then proves that every remaining edge occurs
  exactly twice with opposite direction. Other duplicates, open/non-manifold
  edges, bad indices, wrong counts/types/frames, and inconsistent winding are
  fatal. It does not create `CTETRA`; LiouvilleCOPS ignores those volume cards
  in the historical BDF.
* **Inputs:** one or more DSKs, `--out-dir`, `--frame-kernel`, and either
  `--cspice-root` or explicit `--dskbrief`/`--dskexp`. The expected frame
  defaults to `67P/C-G_CK`.
* **Outputs:** one new `.bdf` and `.bdf.conversion.json` per input. Existing
  destinations are refused. The manifest binds source/PDS-label/output hashes,
  DSK and retained counts, conversion scale/commands, and removed face pairs.
* **Example:**

  ```bash
  python3 tools/convert_dsk_to_bdf.py \
    /data/vtenishe/rosetta_data_end_of_mission_cops/shape_models/SHAP5_PDS_V2/cg_spc_shap5_199k_cart.bds \
    --out-dir /data/vtenishe/rosetta_data_end_of_mission_cops/shape_models/SHAP5_PDS_V2 \
    --frame-kernel /home/vtenishe/SPICE/Kernels/ROSETTA/esa_2022/fk/ROS_V38.TF \
    --cspice-root /home/vtenishe/SPICE/cspice
  ```

Select the result with `ShapeModelFormat = nastran`, not `dsk`. See
`A7_PERFORMANCE.md` for the five converted files and exact audit results.

---

## `benchmark_fit_performance.py`

* **Purpose:** characterize the complete Stage-A7 Python fitter on an actual
  manifested static base. It times startup, facet-basis construction,
  production residual and analytic-Jacobian calls, measures peak RSS and
  independent-process scaling, and estimates a declared campaign's capacity
  and storage.
* **Inputs:** positional manifested base; `--hk-table`; one named profile;
  `--degrees` (default 5 6 7 8); timing `--repeat`; candidate `--jobs`
  (default 1 2 4 8 14); and `--threads-per-process` (default 1). Campaign
  counts are explicit planning inputs, not inferred optimizer behavior.
* **Outputs:** resource/scaling/campaign tables on stdout. Optional `--json`
  and `--markdown` paths must be new and are resolved through the normal
  output root. By default, one production fit per degree is deliberately
  stopped after one evaluation to measure a complete A5 package; these
  temporary directories are deleted and are not scientific fits.
* **Parallel meaning:** jobs are independent fit-like workers. The script
  fixes OpenMP/OpenBLAS/MKL to the requested thread count inside each worker;
  it does not parallelize a single SciPy lambda fit.
* **Example:**

  ```bash
  env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
    python3 tools/benchmark_fit_performance.py \
    /absolute/path/to/output_fit_full_default_base \
    --hk-table /absolute/path/to/COPS_HK_20160930.txt \
    --profile all_available --degrees 5 6 7 8 --repeat 3 \
    --jobs 1 2 4 8 14 --threads-per-process 1 \
    --campaign-profiles 3 --campaign-lambdas 5 \
    --campaign-max-nfev 1000 --campaign-id a7measured \
    --json a7_performance/performance_YYYYMMDD.json \
    --markdown a7_performance/performance_YYYYMMDD.md
  ```

Use `--skip-output-probe` when product-size measurement is unnecessary. See
[`A7_PERFORMANCE.md`](A7_PERFORMANCE.md) for measured tables, recommendation,
limits, and the higher-resolution shape-model inventory.

---

## `rebuild_cops_boundary.py`

Regenerates the angular transport kernels and forward pressures when changing
`SurfaceBC` and can replace the surface geometry with `--shape-model` and
`--shape-model-format`. It runs the existing C++ solver, preserving source
temperature, activity, and selection settings. It does not run a fit. Requires
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
`--shape-model FILE`, `--shape-model-format nastran|dsk`,
`--executable FILE` (default `liouville_cops`), and `--dry-run` (no writes or
execution; external files are not checked). Relative `--out` uses the shared
output root. Other relative source paths use the application root, with
generated input/HK/activity lookup; kernel entries use `KernelDirectory`.

The format is never inferred from a suffix. A shape override is written as an
absolute path to `boundary.in` and `rebuild.json`; the C++ run manifest also
hashes the selected mesh. Changing the shape rebuilds all geometry-dependent
products and requires a new fit based on that output directory.

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
  points not selected by the `--ref` run are hollow. Select the Python fit as `--ref` to display its actual fitting
  mask rather than the base input's exclusion flags. Schema-v2 fit inputs are
  accepted only when their table/report hashes and profile, selection, report,
  and basis identities agree with `fit_manifest.json`. Labels expose profile,
  degree, lambda, and report ID; unmanifested inputs are marked legacy. With
  manifested C++ runs, the table hash is verified against `run_manifest.json`
  and the label exposes degree and coefficient/activity identity. With
  `--save`, the PNG metadata lists the loaded fit-report identities.
* **Options:** run directories or files (default: all `output*/`), `--x time|altitude`, `--ref RUN` (the run that
  provides the reference data, `useNG`/`useRG` marker masks, and altitudes), `--save FILE`, `--theme light|dark`.
* **Examples:**

  ```
  python3 tools/plot_cops_model.py
  python3 tools/plot_cops_model.py output_fit_base output_fit_base_fit output_rough
  python3 tools/plot_cops_model_altitude.py output_M004 --save m004_altitude.png
  ```

## `fit_cops.py`

* **Purpose:** fit either the legacy H2O source map (spherical harmonics up to degree 4), analytical surface
  temperature and CO2 scale, or a controlled normalized H2O residual at degrees 5--8. The residual mode fixes all
  legacy activity, temperature, and CO2 parameters. The misfit is logarithmic, with equal weight per gauge. A
  Tikhonov term `--lambda` pulls the legacy map toward its prior or the new unclipped residual toward zero;
  several values give an L-curve.
* **Input:** a LiouvilleCOPS run made with `WriteKernels = on` and `FlightTimeRotation = none` (e.g. `fit_base.in`
  -> `output_fit_base`): its geometric kernels `kernels/point=NNNN.bin`, `elements.bin`, `COPS_model.dat` (data, fit
  selection `useNG`/`useRG`, T_RG) and the srcCG coefficients of its `cops.in` (prior). With `--data FILE` the
  pressures and the `useNG`/`useRG` flags come from a table instead (`data/COPS_fit_data.txt` of
  `collect_cops_sources.py`, which can be edited); its points must be points of the run, with the same T_RG.
  New bases must have a valid `run_manifest.json`. A pre-W01 base requires the explicit `--allow-legacy-base`
  compatibility option and still receives the complete structural cache audit.
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
    `cops_fit.in` (LiouvilleCOPS input that reruns the fitted model into `<fit dir>_verify`), plus
    `fit_manifest.json` (versioned identities and independent process/optimizer/verification/review states).
    Residual fits also write the strict `ActivityEnhancementFit.txt` exchange file and a per-facet
    `ActivityEnhancementFacets.dat`; the legacy `ActivityFit.txt` format is unchanged;
  * every completed result also writes `fit_report.json`, `activity_facets.dat`,
    `activity_degree_power.dat`, `activity_region_diagnostics.dat`,
    `activity_maps.png`, `pressure_residuals.dat`,
    `pressure_comparison_time.png`, `pressure_comparison_altitude.png`, and
    `pressure_residuals.png`. Tables and PNG metadata share one report identity;
    pressure plots distinguish selected, excluded, saturated-NG, transition,
    and model values;
  * several lambdas: `<out>/lambda=<value>/` with the same files, and `<out>/lcurve.txt` with the maximum-curvature
    lambda.
  `--out` must name a new path. In the Python prediction table, only aggregate `pNG`/`pRG` are recomputed;
  unrecomputed production/species/below-deck fields are `nan`, not inherited base results.
* **Options:** `--out DIR`, `--degree N`, `--enhancement-degree 4..8`, `--fix-legacy-activity`,
  `--enhancement-scale VALUE` (assertion of the derived scale), repeatable
  `--activity-template BUNDLE NAME`, `--lambda L [L ...]`, `--jobs N`, `--fix-temperature`,
  `--ordered-temperature`, `--fix-kappa`, `--tbounds TMIN TMAX` (default 30 400 K), `--max-nfev`, `--ftol`, `--xtol`,
  `--gtol`, `--check-jacobian`, `--check-tolerance`, `--data FILE`, `--profile NAME`, `--apply-exclude-ng`, `--exclude-ng-saturated`, `--hk-table FILE`, `--data-only`, `--allow-legacy-base`.
  A named profile owns the selection and cannot be combined with the two
  legacy NG-selection switches. Its audited policy and reason codes are saved
  in `selection_profile.json`, `fit_manifest.json`, and `fit_data_used.txt`.
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
  python3 tools/fit_cops.py output_fit_full_default_base --profile all_available \
          --enhancement-degree 8 --fix-legacy-activity --fix-temperature --fix-kappa \
          --lambda 1e-4 --jobs 14 --out fit_enh_l8_all_available
  python3 tools/fit_cops.py output_fit_shap5_788k_base --profile all_available \
          --enhancement-degree 8 --fix-legacy-activity --fix-temperature --fix-kappa \
          --activity-template physical_templates_shap5_788k.npz late_ng_not_rg \
          --activity-template physical_templates_shap5_788k.npz cliff_ge_60deg \
          --lambda 1e-5 1e-4 1e-3 --jobs 3 --out fit_l8_physical_templates
  ./liouville_cops output_fit_base_fit/cops_fit.in    # verify a fit with the C++ code
  python3 tools/validate_cops_run.py --fit output_fit_base_fit --verification output_fit_base_fit_verify
  ```

## `run_fit_profile.py`

* **Purpose:** validate and execute the checked-in `fit_profiles/*.json`
  policies through `fit_cops.py`, with distinct output directories and exact
  commands printed before execution.
* **Profiles:** on the 100k reference, `legacy_base_mask` has 839 NG/1,146 RG,
  `all_available` 930/1,146, and `transition_cut_100m` 869/1,146. On the 788k
  base the counts are 839/1,150, 934/1,150, and 873/1,150 because its last four
  positions are outside. Profiles audit the 1,150 rows and selection before
  the shape test; an unexplained data/policy change remains fatal.
* **Parallelism:** `--jobs` defaults to 14 and controls independent lambda
  fits. Profiles run sequentially. A single lambda supplies only one optimizer
  task. Numerical-library threads default to one per worker.
* **Examples:**

  ```bash
  python3 tools/run_fit_profile.py --run output_fit_full_default_base \
    --all-profiles --data-only --jobs 14 --out-prefix a1_selection

  python3 tools/run_fit_profile.py --run output_fit_full_default_base \
    --profile transition_cut_100m --lambda 0 --jobs 14 \
    --out a1_fit_transition_cut_100m

  python3 tools/run_fit_profile.py --run output_fit_full_default_base \
    --all-profiles --lambda 0 --jobs 14 --out-prefix a1_fit --dry-run

  python3 tools/run_fit_profile.py --run output_fit_full_default_base \
    --profile all_available --enhancement-degree 8 \
    --fix-legacy-activity --fix-temperature --fix-kappa \
    --lambda 1e-6 1e-5 1e-4 1e-3 --jobs 14 \
    --out fit_enh_l8_all_available_scan
  ```

The runner exposes the current controls for legacy degree, residual degree,
fixed legacy activity, derived-scale assertion, lambda grid,
temperature constraints, `kappa`, optimizer tolerances, replacement data,
housekeeping table, Jacobian/data-only modes, and legacy-base acceptance.
With a named profile, replacement data must include every base point and the
profile overrides its use flags. The runner does not change the physical model.

## `run_a6_validation.py`

* **Purpose:** execute the complete Stage-A6 implementation-consistency gate
  from one command. The default discovers all Python unit/synthetic tests and
  runs `liouville_cops selftest.in`; `--integration` additionally enables the
  local C++/mission-data workflows guarded by `COPS_RUN_INTEGRATION`.
* **Controls:** `--threads N` sets C++ OpenMP threads (default 2),
  `--integration` enables the external-input tests, `--skip-selftest` omits
  the C++ self-test, and `--report FILE` writes a new JSON evidence record.
  Numerical-library threads are fixed at one during the suite.
* **Outputs:** test products use temporary directories. The optional report
  records commands, environment, timings, return codes, output tails, and
  executable/test hashes. Relative report names are under the configured
  output root; existing reports are never overwritten.
* **Examples:**

  ```bash
  python3 tools/run_a6_validation.py --threads 2
  python3 tools/run_a6_validation.py --integration --threads 2 \
    --report a6_validation/validation_YYYYMMDD.json
  ```

The integration gate includes an exact numeric-row and cache-byte zero-
enhancement forward regression
and a degree-5 Python-fit/C++-replay pressure check with an NG/RG relative
tolerance of `1e-6`. Passing does not independently validate the physical
model, instrument calibration, source localization, or a full descent fit.
See [A6_VALIDATION.md](A6_VALIDATION.md) for every test and threshold.

## `validate_cops_run.py`

* **Purpose:** read-only validation of a static fitting base, or comparison of a Python fit prediction with its
  full C++ rerun. It checks table widths and unique point IDs, exact binary lengths and versions, kernel/table point
  correspondence, face bounds/duplicates, hashes, declared states, and fitted aggregate NG/RG agreement. When both
  parent and verification runs are manifested, their non-activity input identities must match.
* **State change:** base validation changes nothing. Fit verification atomically changes only
  `fit_manifest.json`'s `cpp_prediction_agreement` state and records the verification identities/error maxima.
  It never changes `scientific_review` or makes a physical-validity claim.
* **Examples:**

  ```bash
  python3 tools/validate_cops_run.py output_fit_base
  python3 tools/validate_cops_run.py \
    --fit output_fit_new --verification output_fit_new_verify
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
