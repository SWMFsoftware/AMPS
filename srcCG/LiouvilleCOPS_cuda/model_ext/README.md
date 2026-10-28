# model_ext: dust grain in the nude gauge and the NG = RG crossing constraint

The code used for the 12M grain fits on Voyager (2026-10-27, jobs 120985–120988). These files are byte-identical
copies of `$COPS_OUTPUT_ROOT/model_ext/` (`/data001/highresolu_vt/cops_output/model_ext/`), which is where the jobs
ran them from. Their SHA256 values match the ones recorded in every production fit's
`fit_manifest.json["model_extension"]` and in `<run>/grain/provenance-<job>/model_ext.sha256`.

## Model (user decisions 2026-10-27, handoff §14)

* **Gauges:** NG = surface sources + dust grain; RG = surface sources only.
* **Grain:** grain(t) = max(0, A − B·(t − t_end)) for t ≥ t_end, and 0 before. t_end is the first NG reading after
  the largest NG data gap of the last hour (10:33:10 UTC). A and B are ≥ 0; they are fitted in units of 1e-6 Pa and
  1e-8 Pa per 100 s.
  * The start value of A is the first reading after the gap minus the surface model there, clipped at 0.
  * The start value of B is 0.
* **Crossing:** the model should have the same sign of ln(NG/RG) as the data: NG ≤ RG before the gap and NG ≥ RG
  from t_end on.
  * One-sided penalty rows with weight MU·sqrt(w_NG·w_RG), MU = `--crossing-weight` (10), and an exact Jacobian.
  * The reported costs exclude the penalty.
* **Monotonic penalty** (mono set): it acts on the full NG, surface + grain.
* **Fine tuning:** A and B are refitted together with the late activity groups.
  * Data: the late window plus the NG points from t_end to the window start.
  * Each gauge gets equal total weight.
  * The rows before t_end are asserted unchanged; between t_end and the window, only the grain may change them.

## Files

| File | Role |
|---|---|
| `grain_crossing.py` | `install(fit_cops, crossing_weight, record)`: `Problem` subclass with A and B, a `least_squares` wrapper for the crossing penalty, and a `write_fit` wrapper that writes `grain_term.{txt,json}` and `grain_components.dat` |
| `grain_fit.py` | Driver: loads the unchanged `fit_cops.py`, installs the extension, runs the unchanged `tools/gpu_fit.py`. Adds `--crossing-weight` and `--check-jacobian` |
| `grain_finetune.py` | Driver around the unchanged `tools/fine_tune_gpu.py --method cpu`, with A and B refitted and the earlier-row checks |
| `grain.sbatch` | The 8 refits: `GRAIN_SET=free` (B29, N30, C30, NC31) or `mono` (B29, N30, S30, NS31, `--monotonic-weight 10 --ordered-temperature`); `GRAIN_TEST=1` runs a 4-evaluation test with the finite-difference check |
| `grain_finetune.sbatch` | Fine tuning of every complete `<run>/grain/fits/*_grain_cross10` |

`fit_cops.py`, `cops_cuda`, `gpu_fit.py` and `fine_tune_gpu.py` are not modified. The extension replaces functions
in the running process only.

## Site-specific values (as run)

* `grain_fit.py` and `grain_finetune.py` set `APP = /data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda`.
* The job scripts use Voyager paths, account and partition, and take the code from `$COPS_OUTPUT_ROOT/model_ext`.

These lines were kept as run, so the copies stay identical to the code recorded in the fits. To use the code
elsewhere, change them; the hashes will then differ from the recorded ones.

`../../LiouvilleCOPS` (the CPU workflow) is unchanged. The grain is not in it yet.

---

# How to include or exclude the dust grain

The grain is not a switch in the fitter or in the runner. It is present when a fit runs through `grain_fit.py`, and
absent when it runs through the unchanged `tools/gpu_fit.py`: the runner's `fit:` stage and `jobs/c16.sbatch` both
call `gpu_fit.py`. Both paths use the same base, record store and template bundles, so the two kinds of fit can be
made on one campaign and compared directly.

| | Without grain (standard) | With grain + crossing |
|---|---|---|
| fit driver | `tools/gpu_fit.py` | `model_ext/grain_fit.py` (installs `grain_crossing.py`, then runs `gpu_fit.main()`) |
| free temperatures (B29, N30, C30, NC31) | runner: `cops_pipeline.py submit/run` → `<run>/fits/<FIT>_n2000/` | `sbatch --export=ALL,COPS_RUN=<run>,GRAIN_SET=free model_ext/grain.sbatch` → `<run>/grain/fits/<FIT>_n2000_grain_cross10/` |
| monotonic + ordered T (B29, N30, S30, NS31) | `sbatch --export=ALL,COPS_RUN=<run> jobs/c16.sbatch` → `<run>/c16/fits/<FIT>_n2000_mono10_ordered/` | `... GRAIN_SET=mono model_ext/grain.sbatch` → `<run>/grain/fits/<FIT>_n2000_mono10_ordered_grain_cross10/` |
| fine tuning | runner `finetune:` stage (`tools/fine_tune_gpu.py --method cpu`) → `<run>/finetune/<FIT>_n2000_m5_i5_cpu/` | `sbatch --export=ALL,COPS_RUN=<run> model_ext/grain_finetune.sbatch` (`grain_finetune.py`) → `<run>/grain/finetune/<fit>_m5_i5_cpu/` |
| NG in `COPS_model.dat`, figures | surface only | surface + grain (the split is in `grain_components.dat`) |
| parameters | 25 activity coefficients, temperatures, κ, template amplitudes | the same + `grain_A_uPa`, `grain_B_uPa_per_100s` (after κ) |

**Including the grain** in a fit run by hand: put `model_ext/grain_fit.py` in front of the `gpu_fit.py` command
line. Everything else, including the order of the arguments, stays as it is:

```bash
cd $APP
# standard fit
.venv/bin/python tools/gpu_fit.py --devices 0,1,2,3,4,5,6,7 --precision mixed --record-store $R/store \
    --max-nfev 2000 --wall-time slurm -- $R/base $FITARGS --out <out> [--activity-template <npz> <name> ...]
# the same fit with the grain and the crossing constraint
.venv/bin/python model_ext/grain_fit.py --crossing-weight 10 [--check-jacobian] \
    --devices 0,1,2,3,4,5,6,7 --precision mixed --record-store $R/store \
    --max-nfev 2000 --wall-time slurm -- $R/base $FITARGS --out <out> [--activity-template <npz> <name> ...]
```

The two options belong to `grain_fit.py` and go before `--`:

* `--crossing-weight MU`, default 10. With `0` the crossing rows are all zero, which gives the grain alone. This
  variant has not been run.
* `--check-jacobian` compares the A and B columns of the Jacobian with central finite differences at the start
  and prints the result. It costs two residual evaluations.

The other `gpu_fit.py` options work as before: `--monotonic-weight`, `--ordered-temperature`, templates,
`--restart-from`, `--tbounds`. A restart file of a grain fit contains A and B. Use it only with `grain_fit.py`.
A standard fit's restart has no A and B, so with `grain_fit.py` they start from the default (A = the jump,
B = 0). That combination has not been tested.

**Excluding the grain** means running the standard path: the runner, `jobs/c16.sbatch`, or `tools/gpu_fit.py`
directly. Nothing has to be undone, because nothing outside `model_ext/` is changed. The standard results of the
12M and 20M campaigns were made this way, with no grain, and they stay where they are.

**Fine tuning has to follow its fit:**

* A grain fit is fine-tuned only with `grain_finetune.py` or `grain_finetune.sbatch`. That script needs
  `grain_term.json` next to the restart file and refits A and B together with the late groups.
* A standard fit is fine-tuned with the unchanged `fine_tune_gpu.py`.
* Crossing the two is refused by design:
  * `fine_tune_gpu.py` on a grain fit: the stored NG includes the grain, so the baseline check fails (rtol 1e-6,
    "the restart's stored pressures do not reproduce").
  * `grain_finetune.py` on a standard fit: there is no `grain_term.json`.

**What is fixed in the grain code** (changing any of it changes the model; record the decision in handoff §14):

* The gap: the largest NG data gap in the last hour; t_end is its first NG point after it.
* The form of the grain: A − B·(t − t_end), clipped at 0, with A, B ≥ 0.
* The start values of A and B.
* The units `A_UNIT` and `B_UNIT`.
* The crossing rule: the sign of the data on each side of the gap.
* In the fine tuning: the data selection and the equal total weight per gauge.

# How to change the surface boundary condition (emission law)

`SurfaceBC` sets the angular distribution of the gas leaving each surface element (`../LiouvilleCOPS/src/Liouville.h`):

* `flat`: Lambertian, dF/dΩ = F·cos θ/π, i.e. an isotropic half-Maxwellian. All runs so far use it
  (`runs/shap8_012m_flat_earlier`, `runs/shap8_020m_flat_earlier`; `rebuild.json` `"surface_bc": "flat"`).
* `rough`: uniform over the hemisphere of the external normal, dF/dΩ = F/(2π).

The emission law is part of the geometric kernels Kn, KPhi and KV of the fitting base. The fitter only reads those
kernels, so the law cannot be changed in a fit. **It needs a new base, and therefore a new campaign in a new restart
directory.** The runner refuses to change the emission law, the mesh or the base input in an existing restart
directory, and the existing flat products stay as they are.

The grain is independent of the emission law. On a rough campaign it is included or excluded exactly as above.

## 1. Runner campaign with `--surface-bc rough`

Use the same command as the flat campaign (handoff §8.1), with a new restart directory and `--surface-bc rough`:

```bash
cd $APP
R=$COPS_OUTPUT_ROOT/runs/shap8_012m_rough_earlier          # new directory; never reuse a flat one
SHAPE=$COPS_DATA_ROOT/shape_models/SHAP8_MIARD/cg-spg_mspcd-shap8-v0.4cmb-12Mfacets.bdf
EARLIER="--profile all_available --degree 4 --enhancement-degree 4 --lambda 0 --jobs 1 --tbounds 30 400 --ftol 1e-10 --xtol 1e-10 --gtol 1e-10 --check-tolerance 1e-6"
CAMPAIGN=(--restart-dir $R --shape $SHAPE --surface-bc rough --max-nfev 2000 --allow-genus 4
          --fits B29,N30,C30,NC31 --precision mixed --fit-args "$EARLIER" --output-root $COPS_OUTPUT_ROOT
          --fine-tune --ft-minutes 5 --ft-intervals 5 --ft-method cpu --ft-max-nfev 2000)
.venv/bin/python tools/cluster/cops_pipeline.py plan "${CAMPAIGN[@]}" --gpus-per-node 8 --gpu-memory-gb 80 --host-memory-gb <GB>
.venv/bin/python tools/cluster/cops_pipeline.py submit "${CAMPAIGN[@]}" --nodes 1 --gpus-per-node 8 \
    --slurm-partition gpu --slurm-account general --time-prepare 06:00:00 --time-fit 06:00:00 \
    --wall-time slurm --resubmit 3
$COPS_OUTPUT_ROOT/jobs/submit_campaign.sh $R     # prepare (+ fits, fine tuning chained by the runner) and C16
```

* **What happens:** the runner's `base` stage runs the unchanged `../LiouvilleCOPS/tools/rebuild_cops_boundary.py`
  with `SurfaceBC = rough` and builds the base on one GPU with `native/liouville_cops_cuda`. All later stages are
  rebuilt from that base: record store, visibility, template bundles, then the fits and the fine tuning.
* **Gravity:** gravity and the gravity cache do not depend on the emission law, but the runner rebuilds them in the
  new directory too.
* **Cost on 12M (flat, measured 2026-10-27):** base 48 min, store 3 min, visibility 4 min, gravity 9 min, cache
  3 min, bundle 7 min. About 1 h 15 min before the first fit. Disk: base 157 GB plus store 76 GB.
* **C16 set:** `submit_campaign.sh` also submits `c16.sbatch`. It makes the sin(slope) templates from the new
  base (`$R/c16/bundle/slope_templates.npz`) and runs the C16 fits.
* **Grain fits on the rough base**, once the prepare group and C16 are done:

  ```bash
  sbatch --export=ALL,COPS_RUN=$R,GRAIN_SET=free model_ext/grain.sbatch
  sbatch --export=ALL,COPS_RUN=$R,GRAIN_SET=mono model_ext/grain.sbatch
  sbatch --dependency=afterok:<grain job> --export=ALL,COPS_RUN=$R model_ext/grain_finetune.sbatch
  ```

  `grain.sbatch` checks that `base`, `store`, `bundle` and `c16/bundle` exist in `$R`.
* **Permission:** a rough campaign is a new run with different settings. Per CLAUDE.md §3.5, ask the user before
  submitting it.

## 2. Other places where SurfaceBC is set

* CPU reference (`../LiouvilleCOPS`): `SurfaceBC = flat|rough` in the input file. `fit_base.in` has
  `SurfaceBC = flat` (line 67). The runner comments out that line in the generated `<run>/base/boundary.in` and
  writes its `--surface-bc` value instead. Do not edit `fit_base.in` for this: only its `KernelDirectory` and `ObservationFile` lines may differ from git.
* Surface utility (`../Utility/surface/surface_contributions.py`): add `--set SurfaceBC=rough`. The NG-only and
  RG-only counts then use the rough emission law.
* Restart files, translations and template bundles from a flat campaign belong to the flat base. The runner checks
  `SurfaceBC` when a base, store or bundle made elsewhere is passed (`--base`, `--store`, `--bundle`), and it
  refuses one with the other law.
