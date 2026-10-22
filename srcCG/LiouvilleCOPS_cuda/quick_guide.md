# LiouvilleCOPS_cuda: quick guide

Exact commands to set up, run and post-process GPU fits of the COPS model.
Background, evidence and stage status are in `README.md`, `RUN.md`, `STATUS.md` and
`VALIDATION.md`. The CPU workflow in `../LiouvilleCOPS` is unchanged: the GPU path runs
the same `fit_cops.py` (unmodified), and only the pressure/Jacobian evaluation moves to
the GPU.

Shell: tcsh. Two interpreters:
* `.venv/bin/python` (3.11, CuPy): everything in this directory;
* `/usr/bin/python3` (3.8): the CPU tools in `../LiouvilleCOPS`.

## 0. Session setup (every new shell)

```tcsh
cd ~/ROS-END-OF-MISSION/AMPS/srcCG/LiouvilleCOPS_cuda
setenv COPS_OUTPUT_ROOT /data/vtenishe/rosetta_data_end_of_mission_cops/cops_output
setenv OMP_NUM_THREADS 1
setenv OPENBLAS_NUM_THREADS 1
setenv MKL_NUM_THREADS 1
set O = $COPS_OUTPUT_ROOT
set BASE = $O/output_fit_shap5_788k_base
set BUNDLE = $O/cuda_c05/restart_rebuilt/restart_templates_shap5_788k.npz
set CPU_BUNDLE = $O/physical_templates_shap5_788k.npz_heal_defect
```

Notes:
* `COPS_OUTPUT_ROOT` is required: the reference path helpers resolve the base run with it.
* The BLAS variables are set to 1, as in the CPU campaign. The serial optimizer gains
  nothing from more threads.
* `$BUNDLE` is the restart bundle (C05) with `late_ng_not_rg`, `late_rg_not_ng`,
  `cliff_ge_45deg` and `cliff_ge_60deg`, made from the admitted 788k gravity cache.
* `$CPU_BUNDLE` holds the `late_ng_not_rg` column that the CPU campaign's N fit used. Its
  `late_ng_not_rg` is bitwise equal to `$BUNDLE`'s. Its cliff columns are HEAL-defective,
  and the fitter refuses them.

## 1. Build ("compile")

### CUDA environment (once, or after a dependency change)

```tcsh
uv python install 3.11
uv venv --python 3.11 .venv
uv pip sync --python .venv/bin/python --require-hashes requirements-cuda.lock
```

The CUDA kernels (`cops_cuda/kernels/*.cu`) have no separate build step: CuPy compiles
them with NVRTC (`--fmad=false`, FP64) on first use in each process and caches the result.
Check the GPUs, the compilation and the code:

```tcsh
.venv/bin/python tools/cuda_doctor.py
.venv/bin/python -m pytest -q tests          # 142 tests; about 30 s; uses both GPUs
```

### CPU executables (only after a source change in ../LiouvilleCOPS)

These are needed for the C++ replay (section 4) and for SPICE ω (section 2):

```tcsh
make -C ../LiouvilleCOPS liouville_cops spice_angular_velocity
```

## 2. One-time products (only when the base run, mesh or gravity body changes)

These already exist for the 788k base. Skip this section unless something changed.

| Product | Path |
|---|---|
| admitted 788k gravity cache | `$O/cuda_c04/equivalence_20261020T021855Z/admitted_gravity_cache_shap5_788k.npz` |
| restart template bundle | `$BUNDLE` |

To rebuild them:

```tcsh
set SHAPE = /data/vtenishe/rosetta_data_end_of_mission_cops/shape_models/SHAP5_PDS_V2/cg_spc_shap5_788k_cart.bdf
# (a) GPU gravity field at all 787,754 facets (about 53 min on one A6000); prints the new directory
.venv/bin/python tools/c04_gravity_equivalence.py gpu --run $BASE --gravity-shape $SHAPE --device 0
set EQ = <the printed cuda_c04/equivalence_<UTC> directory>
# (b) optional, when a CPU cache of the same body exists: facet-by-facet admission
.venv/bin/python tools/c04_gravity_equivalence.py compare --dir $EQ \
    --cpu-cache <CPU cache .gravity.npz> --cpu-bundle <CPU bundle .npz>
# (c) CPU-format cache with the compare report as evidence
.venv/bin/python tools/c04_write_gpu_cache.py --field-dir $EQ --run $BASE --gravity-shape $SHAPE \
    --out $EQ/admitted_gravity_cache.npz --evidence $EQ/compare_<UTC>.json
# (d) restart bundle: NG-only and RG-only masks (GPU) plus cliff masks from the cache
.venv/bin/python tools/c05_write_bundle.py --run $BASE --gravity-cache $EQ/admitted_gravity_cache.npz \
    --gravity-shape $SHAPE --out $O/cuda_c05/<new dir>/restart_templates.npz
```

## 3. Run fits

All fits use `tools/gpu_fit.py`, which runs the unmodified `fit_cops.py`:

```text
.venv/bin/python tools/gpu_fit.py [driver options] -- <fit_cops.py arguments>
```

| Driver option | Default | Meaning |
|---|---|---|
| `--device N` | 0 | GPU of a single-GPU fit |
| `--devices 0,1` | (none) | one fit split over both GPUs: Jacobian about 1.7× faster, bitwise identical to one GPU |
| `--max-nfev N` | **2000** | optimizer evaluation limit (user decision 2026-10-20: 2000 for fits, 1000 for CPU comparison) |
| `--template-start C` | **1.0** | optimizer start of each template amplitude (lower bound stays 0; the prior stays 0) |
| `--cpu-comparison` | off | exactly the CPU fitter's command: `--max-nfev 1000 --template-start 0` |

Fixed by the user's decisions:
* the temperature bounds stay `--tbounds 30 400` (2026-10-20);
* use one `--lambda` value per process.

Each fit needs about 21 GB on its GPU (kernel data 9 GB plus evaluator arrays 13.2 GiB,
as measured by nvidia-smi). With the C08 fused kernels a B29 fit takes about 16 min
for 2000 evaluations; template fits take a little longer. The fit matrix runs fastest as independent single-GPU fits,
two at a time (C07 benchmark).

The common fitter arguments are:

```tcsh
set FITARGS = "--profile all_available --degree 4 --enhancement-degree 4 --lambda 0 --jobs 1 --tbounds 30 400 --ftol 1e-10 --xtol 1e-10 --gtol 1e-10 --check-tolerance 1e-6"
```

### 3a. Production fits (max_nfev 2000, templates start at 1)

Run two at a time, one per GPU. `nohup … >& log &` keeps a fit running after logout.

```tcsh
# B29: free degree-4 baseline, T, kappa
nohup .venv/bin/python tools/gpu_fit.py --device 0 -- $BASE $FITARGS:x \
    --out $O/fit_788k_gpu_B29_n2000 >& $O/fit_788k_gpu_B29_n2000.log &
# N30: B29 + late NG-only template
nohup .venv/bin/python tools/gpu_fit.py --device 1 -- $BASE $FITARGS:x \
    --out $O/fit_788k_gpu_N30_n2000 --activity-template $BUNDLE late_ng_not_rg \
    >& $O/fit_788k_gpu_N30_n2000.log &
```

When those two have finished:

```tcsh
# C30: B29 + cliff >= 60 deg
nohup .venv/bin/python tools/gpu_fit.py --device 0 -- $BASE $FITARGS:x \
    --out $O/fit_788k_gpu_C30_n2000 --activity-template $BUNDLE cliff_ge_60deg \
    >& $O/fit_788k_gpu_C30_n2000.log &
# NC31: B29 + late NG-only + cliff >= 60 deg
nohup .venv/bin/python tools/gpu_fit.py --device 1 -- $BASE $FITARGS:x \
    --out $O/fit_788k_gpu_NC31_n2000 --activity-template $BUNDLE late_ng_not_rg \
    --activity-template $BUNDLE cliff_ge_60deg >& $O/fit_788k_gpu_NC31_n2000.log &
```

Labelled extensions, run the same way:

```tcsh
# R30: B29 + late RG-only template
nohup .venv/bin/python tools/gpu_fit.py --device 0 -- $BASE $FITARGS:x \
    --out $O/fit_788k_gpu_R30_n2000 --activity-template $BUNDLE late_rg_not_ng \
    >& $O/fit_788k_gpu_R30_n2000.log &
# C30 at 45 deg (cliff sensitivity)
nohup .venv/bin/python tools/gpu_fit.py --device 1 -- $BASE $FITARGS:x \
    --out $O/fit_788k_gpu_C30_45deg_n2000 --activity-template $BUNDLE cliff_ge_45deg \
    >& $O/fit_788k_gpu_C30_45deg_n2000.log &
```

`$FITARGS:x` splits the variable into words (tcsh). Every `--out` must be a **new**
directory, because the fitter refuses to reuse one.

### 3b. Comparison fits against the CPU (max_nfev 1000, templates start at 0)

These are identical to the CPU campaign's B and N commands:

```tcsh
nohup .venv/bin/python tools/gpu_fit.py --device 0 --cpu-comparison -- $BASE $FITARGS:x \
    --out $O/fit_788k_gpu_B29_cpucmp >& $O/fit_788k_gpu_B29_cpucmp.log &
nohup .venv/bin/python tools/gpu_fit.py --device 1 --cpu-comparison -- $BASE $FITARGS:x \
    --out $O/fit_788k_gpu_N30_cpucmp --activity-template $CPU_BUNDLE late_ng_not_rg \
    >& $O/fit_788k_gpu_N30_cpucmp.log &
```

Comparison fits that already exist: `$O/fit_shap5_788k_gpu_B` and `$O/fit_shap5_788k_gpu_N`.
`$O/fit_shap5_788k_gpu2_B` is the same B29 fit on two GPUs, bitwise identical.

### 3c. One fit on both GPUs (when one result is needed quickly)

```tcsh
.venv/bin/python tools/gpu_fit.py --devices 0,1 -- $BASE $FITARGS:x --out $O/<new dir>
```

### Progress and failure

* The log shows the fitter's own messages: the consistency check, the selection,
  `initial cost …`, the TRF summary, then a `gpu_fit:` line with wall time and GPU call
  counts.
* `fit_manifest.json` → `states` gives `process` (`running`, `complete` or `failed`) and
  `optimizer` (`converged` or `not_converged`).
* `compute_backend` records the GPU backend, the devices, the start vector and
  `max_nfev`.
* Reaching the evaluation limit is **not** convergence.

## 4. Post-process the results

Each fit directory contains the fitter's standard products:
* `fit_result.txt`: parameters, σ, termination, bounds, rank and singular values;
* `ActivityFit.txt`, plus `ActivityEnhancementFit.txt` and `ActivityEnhancementFacets.dat`
  when templates are fitted;
* `COPS_model.dat`: fitted NG/RG pressures;
* `pressure_comparison_time.png`, `pressure_comparison_altitude.png` and
  `pressure_residuals.png`, with `.dat` files;
* `activity_maps.png` and the `activity_*.dat` diagnostics;
* `fit_data_used.*` and `fit_profiles_used.png`;
* `cops_fit.in`: the C++ replay input;
* `fit_manifest.json` and `fit_report.json`.

For a template fit, note:
* the "initial" column of `fit_result.txt` is the prior (templates 0);
* the actual optimizer start is in `fit_manifest.json` → `compute_backend.fit_start`.

### 4a. Compare two fits (GPU vs CPU, GPU vs GPU, or two models)

```tcsh
.venv/bin/python tools/compare_fits.py $O/fit_shap5_788k_restored_B $O/fit_shap5_788k_gpu_B
.venv/bin/python tools/compare_fits.py $O/fit_788k_gpu_B29_n2000 $O/fit_788k_gpu_N30_n2000 \
    --json $O/fit_788k_gpu_N30_n2000/compare_vs_B29.json
```

The tool prints:
* whether the inputs, selection, prior and optimizer settings are the same;
* the termination of both fits;
* the final cost and its relative difference;
* which parameters sit on a bound in each fit;
* the parameter differences in units of σ;
* parameters present in only one fit (for example a template amplitude);
* the maximum relative difference of the predicted pressures.

### 4b. Full C++ replay of a fitted state (independent check, CPU)

```tcsh
cd ../LiouvilleCOPS
./liouville_cops $O/fit_788k_gpu_B29_n2000/cops_fit.in            # writes <fit>_verify; Threads=16 is set in the input
/usr/bin/python3 tools/validate_cops_run.py --fit fit_788k_gpu_B29_n2000 \
    --verification fit_788k_gpu_B29_n2000_verify --tolerance 1e-6
cd -
```

The validator compares the fitted aggregate pressures with the C++ forward model and
records the agreement in `fit_manifest.json` (`cpp_prediction_agreement`). Run the replay
when the CPU campaign is not using the same cores.

### 4c. Pressure plots: base, fit and C++ replay against the measurements

```tcsh
cd ../LiouvilleCOPS
/usr/bin/python3 tools/plot_cops_model.py output_fit_shap5_788k_base fit_788k_gpu_B29_n2000 \
    fit_788k_gpu_B29_n2000_verify --ref fit_788k_gpu_B29_n2000 --x time --theme light \
    --save fit_788k_gpu_B29_n2000/comparison_time.png
/usr/bin/python3 tools/plot_cops_model.py output_fit_shap5_788k_base fit_788k_gpu_B29_n2000 \
    fit_788k_gpu_B29_n2000_verify --ref fit_788k_gpu_B29_n2000 --x altitude --theme light \
    --save fit_788k_gpu_B29_n2000/comparison_altitude.png
cd -
```

Run names are relative to `$COPS_OUTPUT_ROOT`. Add more fit directories to overlay
several models in one figure.

### 4c1. Publication figures (.eps and .png): the format of all runs since 2026-10-20

The tools write these figures automatically:

| Tool | Figures |
|---|---|
| `gpu_fit.py` | `<fit>/figures/pressures_{time,altitude}`: base, fit and data (`--no-figures` skips them) |
| `c05_write_bundle.py` | `<bundle stem>_groundtrack_{3d,latlon}`, `_track.dat`, `_impact.json`, and the template maps `<bundle stem>_templates` (+ `_caption.txt`) next to the bundle (`--no-figures` skips them) |
| `restart_tool.py plot` | `<out>_{time,altitude}` |
| `restart_tool.py enhance` | `<out>_{time,altitude}` and the region map `<out>_region` |

Combined figures and regeneration:

```tcsh
# several fits on one figure: --fit LABEL=DIR[:solid|dashed|dotted|dashdot]
.venv/bin/python tools/pub_figures.py pressures --base $BASE --bundle $BUNDLE \
    --fit "N30, 788k=$O/fit_788k_gpu_N30_n2000_r2:dashed" --fit "N30, 1M=$O/fit_001m_gpu_N30_n2000" \
    --out $O/<dir>/N30_788k_vs_1m
# end-of-mission ground track and Rosetta impact point on the base's shape model (GPU, about 30 s)
.venv/bin/python tools/pub_figures.py groundtrack --base $BASE --bundle $BUNDLE --out $O/<dir>/groundtrack --device 0
# activity-template maps (late-visibility, slope, cliff) with ground track and impact (GPU, about 40 s)
.venv/bin/python tools/pub_figures.py templates --base $BASE --bundle $BUNDLE --out $O/<dir>/templates --device 0
```

Format (`cops_cuda/figures.py`):
* Each figure is written as `.eps` (vector; rendered images embedded) and `.png` (300 dpi).
  It uses an Okabe–Ito colour-blind safe palette and no transparency.
* **Pressure figures.**
  * Rows: NG, RG and NG/RG.
  * The time figure has a full-descent column and a zoom on the last 25 min. The altitude
    figure uses a log axis with the descent running left to right.
  * The **late-visibility window** is shaded and labelled with its UTC and altitude ranges.
    The window is taken from `--bundle`'s `late_ng_not_rg` interval, else `--late-start/
    --late-stop`, else c05_write_bundle's default 10:33:10–10:39:26.
* **Ground-track figures.** The track is the sub-spacecraft intercept: the first surface hit
  toward the body-frame origin, as in SPICE `subpnt` INTERCEPT.
  * `_3d` has three panels: the nucleus viewed down the final approach, the template
    regions, and a close-up with 1-minute marks.
  * `_latlon` has the global map and the impact area.
  * The **Rosetta impact point** is the straight-line extrapolation of the last two rows to
    the shape model: 10:39:26.4 UTC, lon 12.23°, lat 31.75° on SHAP5 1M and 788k.
  * The NG-only and RG-only facets of the bundle are coloured.
* **Template maps** (`templates`): the publication form of the builder's bundle map.
  * Rows: late-visibility templates, effective gravitational slope (from the bundle's
    gravity cache), and cliff templates. Left: global map; right: impact area.
  * The header defines "late visibility": NG-only facets have Kn > 0 and KΦ = KV = 0 at
    one of the window's points; RG-only facets are the reverse.
  * `_caption.txt` gives the full definitions (kernels, slope, track, impact).
* Existing figure files are refused unless `--overwrite` is given.

## 4c2. Restart files (C12): save, restart, plot, enhance, translate

Every `gpu_fit.py` fit writes `<out>/fit_restart.npz`: the fitted state, templates,
activity map, mesh identity and predictions.

```tcsh
# continue or redo a fit from a restart (same mesh); templates come from the restart unless given
.venv/bin/python tools/gpu_fit.py --device 0 --restart-from $O/<fit>/fit_restart.npz -- $BASE $FITARGS:x --out $O/<new fit>
# inspect, and compare several states (no recomputation)
.venv/bin/python tools/restart_tool.py show $O/<fit>/fit_restart.npz
.venv/bin/python tools/restart_tool.py plot $O/<fitA>/fit_restart.npz $O/<fitB>/fit_restart.npz \
    --labels A B --out $O/<new dir>/compare
# raise H2O activity on the facets seen by the NG and not by the RG from points below 100 m;
# plots with and without the change; optionally save the edited restart (+ its template bundle)
.venv/bin/python tools/restart_tool.py enhance $O/<fit>/fit_restart.npz --run $BASE \
    --region ng_not_rg --alt-max 100 --mode add --value 2 --out $O/<new dir>/enh \
    --save $O/<new dir>/enh_restart.npz
# carry a restart to another mesh (prefer a bundle recomputed on the target mesh)
.venv/bin/python tools/restart_tool.py translate $O/<fit>/fit_restart.npz --to-run <other base> \
    --save $O/<new dir>/translated_restart.npz [--target-bundle <bundle on that mesh>] --evaluate
```

Regions are `ng_not_rg`, `rg_not_ng`, `ng_seen`, `rg_seen` and `any_seen`, selected with
`--alt-min/--alt-max` and/or `--utc-start/--utc-stop`.

Enhancement modes:
* `add V ≥ 0`: adds V·A_scale to the linear activity. It also acts on facets whose activity
  is clipped.
* `factor V > 0`: multiplies the positive activity on the region by V; V < 1 reduces it.

The region is confined to the window. With `--confine total` (default), the summed change of
the kept facets alters every outside point's NG and RG pressure by at most `--epsilon`
(default 1e-3). With `--confine facet`, each facet alone stays below ε. The summary reports
the change inside and outside the window.

A saved edit can be refit (`--restart-from`) with the new region's amplitude free.

## 4d. Build the fitting base on the GPU (C11)

```tcsh
make -C native                                   # builds native/liouville_cops_cuda (CUDA 12.4, CSPICE)
cd ../LiouvilleCOPS
env COPS_CUDA_DEVICE=0 python3 tools/rebuild_cops_boundary.py --input fit_base.in --surface-bc flat \
    --shape-model <BDF> --shape-model-format nastran --out <new base dir> --threads 16 \
    --executable ../LiouvilleCOPS_cuda/native/liouville_cops_cuda
cd -
.venv/bin/python tools/c11_compare_bases.py <C++ base dir> <GPU base dir>       # optional comparison
```

The output is the same set of products as the C++ build: `elements.bin` is bitwise
identical, and the kernels agree at rounding level. See `native/README.md`.

## 4e. A new surface model, everything on the GPU (SHAP5 1M, 2026-10-20)

The exact sequence used for `cg_spc_shap5_001m_cart.bdf` (1,194,250 facets). The times are
for one A6000.

```tcsh
set SHAPE = /data/vtenishe/rosetta_data_end_of_mission_cops/shape_models/SHAP5_PDS_V2/cg_spc_shap5_001m_cart.bdf
set BASE  = $O/output_fit_shap5_001m_base
# (1) fitting base (C11), 5.5 min; 356,553,303 kernel records, 14 GB
cd ../LiouvilleCOPS
env COPS_CUDA_DEVICE=0 python3 tools/rebuild_cops_boundary.py --input fit_base.in --surface-bc flat \
    --shape-model $SHAPE --shape-model-format nastran --out $BASE --threads 16 \
    --executable ../LiouvilleCOPS_cuda/native/liouville_cops_cuda >& $BASE.build.log
cd -
# (2) optional, no gravity needed: the NG-only/RG-only templates, so N30 can start at once
.venv/bin/python tools/c05_write_bundle.py --run $BASE --visibility-only \
    --out $O/cuda_c05/restart_001m/visibility_templates_shap5_001m.npz
# (3) gravity field at all facets, 2.0 h (1.2e6 x 1.2e6 facet pairs); add --devices 0,1 to split the
#     targets over both GPUs (bitwise the same field, half the time)
.venv/bin/python tools/c04_gravity_equivalence.py gpu --run $BASE --gravity-shape $SHAPE --device 0 \
    >& $O/gravity_001m_gpu.log
set EQ = <the printed cuda_c04/equivalence_<UTC> directory>
# (4) cache (1 min) and full bundle (1 min); no CPU cache of this body exists, so no --evidence
.venv/bin/python tools/c04_write_gpu_cache.py --field-dir $EQ --run $BASE --gravity-shape $SHAPE \
    --out $EQ/admitted_gravity_cache_shap5_001m.npz
.venv/bin/python tools/c05_write_bundle.py --run $BASE --gravity-cache $EQ/admitted_gravity_cache_shap5_001m.npz \
    --gravity-shape $SHAPE --out $O/cuda_c05/restart_001m/restart_templates_shap5_001m.npz
# (5) the four fits exactly as in 3a, with $BASE and the 1M bundle; about 25 min each
```

The GPU evaluator plan of a fit on this base is 19.2 GiB (bind 138 s), within one 48 GB A6000. B29 and N30 can run during step (3) on the other
GPU, using the visibility bundle from (2).

## 4e1. Mixed precision and very large meshes (C14 experiment, SHAP8 12M, 2026-10-21)

Not the validated path (VALIDATION.md covers FP64 only); every product records its precision policy.

```tcsh
set SHAPE = /data/vtenishe/rosetta_data_end_of_mission_cops/shape_models/SHAP8_MIARD/cg-spg_mspcd-shap8-v0.4cmb-12Mfacets.bdf
set BASE  = $O/output_fit_shap8_012m_base                 # (1) base as in 4e: 4.33e9 records, 173 GB, 65 min
# (3) mixed gravity, both GPUs, 18 min (FP64 would take ~4 days); genus 4 admitted explicitly
.venv/bin/python tools/c04_gravity_equivalence.py gpu --run $BASE --gravity-shape $SHAPE --devices 0,1 \
    --precision mixed --near-factor 6 --allow-genus 4 >& $O/gravity_012m_mixed.log
# (4) cache and bundle: as in 4e, each with --allow-genus 4
# (5) fits: FP32 records, FP64 arithmetic, both GPUs (44.9 GB each; bind 12 min)
.venv/bin/python tools/gpu_fit.py --devices 0,1 --headroom-fraction 0.06 --precision mixed -- $BASE <3a args> ...
```

* `--precision mixed` on `gpu_fit.py` sets `COPS_CUDA_RECORD_PRECISION=fp32` (inherited by the GPU workers).
* FP32 record store (C15/P1): `tools/c15_record_store.py --run $BASE --out $O/cuda_c15/store_<name>_fp32`
  once per base, then `gpu_fit.py --precision mixed --record-store <store> ...` (3M: 14x faster loading). Fits
  validate the base once (`--audit-cache auto`, default; `--revalidate` forces a full check).
* Streaming (C15/P2): with `--precision mixed`, records that do not fit the GPUs stay in pinned host memory and are
  streamed each pass (automatic; results bitwise the resident ones). 20M on the local 2 A6000: about 3.7 h per fit.
* Fine tuning of a 12M restart: `tools/fine_tune_gpu.py ... --precision mixed --devices 0,1 --headroom-fraction 0.06`
  (other options as in §4f). With streaming (C15/P2b) the restart and fine-tuning kernels also run on bases larger
  than the GPUs: `restart_tool.py enhance|translate --precision mixed` on one GPU, `fine_tune_gpu.py --precision
  mixed` on 1–2 GPUs (bitwise the resident results; slower by the transfer).
* `gpu_fit.py` validates the base with `cops_cuda/base_validation.py` (identical result to the reference's,
  about 1 min on 12M instead of 1.8 h); `--reference-validation` runs the reference's Python loop.
* Accuracy: gravity slope ≤ 1.1e-5° from FP64 (1M), cliff masks identical (1M, 3M); 1M B29 fit cost within
  1.5e-7 relative, parameters within 8.8e-5 σ of FP64.

## 4e2. One fit over MPI ranks, one GPU per rank (C15/P3, 2026-10-21)

For several nodes (Voyager: 8 A100 per node). Every rank starts `gpu_fit.py`; rank 0 runs the fitter and the first
block of rows, the other ranks serve their blocks. Locally, ranks map onto `--devices` (default: all visible).

```tcsh
# this machine, 2 ranks (environment forwarded with -x)
mpirun -n 2 --bind-to none -x COPS_OUTPUT_ROOT -x OMP_NUM_THREADS .venv/bin/python tools/gpu_fit.py \
    --transport mpi --precision mixed --record-store $STORE -- $BASE <3a args> --out $O/<name>
# a cluster (Slurm): sbatch --nodes=N --ntasks-per-node=<GPUs per node> tools/cluster/fit_mpi.slurm \
#     <base> <store> <out> [fit_cops args]; set its SITE lines first
```

* `--activity device` (MPI default, user decision Q3): each GPU forms the activity of its own facets; checked
  against the host expression at the first evaluation and every 100th (`compute_backend` records the checks).
  `--activity host` keeps the exact host field (broadcast each evaluation), bitwise the local two-GPU fit.
* mpi4py is built against the MPI that runs the job:
  `env MPICC=mpicc uv pip install --python .venv/bin/python --no-binary mpi4py mpi4py==4.0.3`.
* A failure on any rank ends the whole job (`MPI_Abort`); the manifest of rank 0 records `states.process = failed`.

## 4e3. A whole campaign on a cluster: the runner (C15/P0, 2026-10-22; CPU-tested, not yet run on GPUs)

`tools/cluster/cops_pipeline.py` runs base → record store → bundles (gravity) → B29/N30/C30/NC31 → fine tuning for
one surface model and emission law, resumable from one restart directory. Details: `tools/cluster/README.md`.

```tcsh
set R = $SCRATCH/shap8_020m_flat
.venv/bin/python tools/cluster/cops_pipeline.py plan --restart-dir $R --shape <BDF> --surface-bc flat \
    --max-nfev 2000 --fine-tune --allow-genus 4                      # stages, state, exact commands; no GPU
.venv/bin/python tools/cluster/cops_pipeline.py run --restart-dir $R  # inside an allocation; rerun to resume
.venv/bin/python tools/cluster/cops_pipeline.py status --restart-dir $R
.venv/bin/python tools/cluster/cops_pipeline.py submit --restart-dir $R --nodes 3 --wall-time slurm --resubmit 5
# bases/stores copied from this machine: --base DIR --store DIR [--visibility-bundle F --bundle F] (identity-checked)
```

* Hardware is detected at run time (Slurm, nvidia-smi/CUDA, /proc/meminfo; `--nodes --gpus-per-node --gpu-memory-gb
  --host-memory-gb` override). `--precision auto` (default): FP64 when the records fit resident on the fit GPUs, else
  FP32 records with streaming; refused when the pinned host buffers do not fit a node.
* `--gravity-precision mixed` (default, user decision 2026-10-22) or `fp64`, independent of `--precision`.
* `--wall-time HH:MM:SS|slurm` (default none): fits stop at their best state before the deadline (`gpu_fit.py
  --deadline`, exit 75) and continue as the next segment with `--restart-from` and the remaining evaluations; a
  segmented fit is not identical to one continuous run.

## 4f. Late-interval fine tuning of a completed fit (C13, 2026-10-21)

The second fitting stage of `../LiouvilleCOPS/FINE_TUNING_STEP.md`, with the kernel-response preparation on the
GPU (`tools/fine_tune_gpu.py`; `cops_cuda/fine_tuning.py`). It freezes a fit's restart file and changes only the
activity of facets first seen in the last `--minutes`, split into `--intervals` groups (NG-, RG-dominated, mixed).

```tcsh
# the CPU branch's method (default, user decision 2026-10-21): strict support, additive c >= 0,
# --third-group all, CPU group_basis + fit_corrections; earlier pressures bitwise unchanged
.venv/bin/python tools/fine_tune_gpu.py $O/output_fit_shap5_001m_base \
    --restart-from $O/fit_001m_gpu_NC31_n2000/fit_restart.npz --minutes 5 --intervals 5 \
    --devices 0,1 --out $O/finetune_gpu_001m_NC31_m5_n5_cpu            # add --prepare-only to stop before the fit
# variant: multiplier + floor f*A0 + c*A_scale on disjoint groups, support by significance
.venv/bin/python tools/fine_tune_gpu.py ... --method multiplier --confine total --epsilon 1e-3 \
    --gauges joint --out ...                                            # or --gauges independent (NG, RG fits)
```

* Products as the CPU tool's: `fine_tuning.json`, `fine_tuning_prepared.npz`, `fine_tuning_state.npz`,
  `pressure_comparison.txt`, `ActivityFit.txt`, `ActivityEnhancementFit.txt`, `cops_fine_tuning.in` (C++ replay:
  `./liouville_cops <out>/cops_fine_tuning.in`), plus `figures/` in the publication format.
* Checks: restart/base identities, the restart's pressures reproduced (rtol 1e-6), the kernel set against
  `run_manifest.json`, exact zero early response, an independent GPU evaluation of the final field.
* 3M needs `--devices 0,1 --headroom-fraction 0.15`. 1M on two GPUs: about 2.5 min.
* **Result on SHAP5 1M (N30, NC31), last 5 min:** the strict support is **empty** (0 facets first seen in the window;
  the CPU tool gives the same counts). It depends on geometry only, so it is empty for every fit on this base.
  With `--method multiplier --confine total`, ε = 1e-3: 146 facets, late log-rms unchanged (NG 0.6039 → 0.6038);
  ε = 1e-2: 840 facets, NG 0.6039 → 0.6005. See STATUS.md C13.

## 5. Optional numerical re-checks (after a code or driver change)

```tcsh
.venv/bin/python -m pytest -q tests
.venv/bin/python tools/c06_full_run.py   --capture $O/cuda_reference/c00_20261019T064536Z
.venv/bin/python tools/c07_multi_gpu_check.py --capture $O/cuda_reference/c00_20261019T064536Z
```

* `c06_full_run.py`: GPU vs CPU on the full fit set (about 10 min).
* `c07_multi_gpu_check.py`: one vs two GPUs, bitwise (about 8 min); needs both GPUs free.

The approved budgets are in `VALIDATION.md`:
* pressure ≤ 512 eps·M;
* Jacobian ≤ 1024 eps·M_J per entry;
* gravity ≤ 16 eps·M.

## Open decisions and known behaviour

* All fits so far end with Tmax = 30 K and T75 = 400 K on their bounds. The bounds are
  kept (user decision); report the temperatures as boundary-limited.
* With template start 0, the NG-only and RG-only amplitudes cannot leave 0: the fitted
  baseline is negative on those facets, and the single clip zeroes their derivative. The
  production default start of 1.0 addresses that (user decision 2026-10-20).
* The end-to-end CPU-vs-GPU fit budget is still provisional. It is proposed after the CPU
  B fit finishes (`compare_fits.py`, section 4a).
