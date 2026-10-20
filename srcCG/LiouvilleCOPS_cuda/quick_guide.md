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
| `--devices 0,1` | (none) | one fit split over both GPUs: about 1.9× faster per fit, bitwise identical to one GPU |
| `--max-nfev N` | **2000** | optimizer evaluation limit (user decision 2026-10-20: 2000 for fits, 1000 for CPU comparison) |
| `--template-start C` | **1.0** | optimizer start of each template amplitude (lower bound stays 0; the prior stays 0) |
| `--cpu-comparison` | off | exactly the CPU fitter's command: `--max-nfev 1000 --template-start 0` |

Fixed by the user's decisions:
* the temperature bounds stay `--tbounds 30 400` (2026-10-20);
* use one `--lambda` value per process.

Each fit needs about 21 GB on its GPU (kernel data 9 GB plus evaluator arrays 13.2 GiB,
as measured by nvidia-smi), and roughly 18–25 min for 1000 evaluations (about twice
that for 2000). The fit matrix runs fastest as independent single-GPU fits,
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
