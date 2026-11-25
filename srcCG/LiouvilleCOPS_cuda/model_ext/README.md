# model_ext: dust grain in the nude gauge and the NG = RG crossing constraint

The code used for the grain fits on Voyager and for every later model change built on them (2026-10-27 …
2026-11-11). The jobs ran the files from `$COPS_OUTPUT_ROOT/model_ext/`, `$COPS_OUTPUT_ROOT/jobs/`,
`$COPS_OUTPUT_ROOT/figures_cluster/scripts/` and `$COPS_OUTPUT_ROOT/mpi_multinode/`
(`COPS_OUTPUT_ROOT=/data001/highresolu_vt/cops_output`), outside git. This folder keeps them under version control:

* `grain_crossing.py`, `grain_fit.py`, `grain_finetune.py`, `grain.sbatch`, `grain_finetune.sbatch` (committed in
  `af48d01`): byte-identical to the run copies (no header; rows in `RUN_COPIES.sha256` since 2026-11-11). Their
  SHA256 values match the ones recorded in every production fit's `fit_manifest.json["model_extension"]` and in
  `<run>/grain/provenance-<job>/model_ext.sha256`.
* Every other file (added 2026-11-06 and later; sections "Extensions after 2026-10-27", "Job scripts", "Figure
  scripts", "MPI helpers" below) is the run copy **plus a header block of comments** at the top (after the `#!` line):
  what the file is for, which runs used it, how to run it, and the SHA256 of the run copy. Nothing else is changed.
  The block lies between the lines `# ==== version-controlled copy: header added <date>, …` and
  `# ==== end of header ====` (`//` instead of `#` in `kernels/*.cu`; in `PERFORMANCE.md` an HTML comment
  `<!-- ==== version-controlled copy: … ` … `==== end of header ==== -->`).
  * `RUN_COPIES.sha256` lists, for each such file, the SHA256 of its run copy and the run copy's path.
  * `./check_copies.sh` removes the header blocks (all three forms; from 2026-11-11, before only `#`) and checks every
    file against those values (OK / DIFFERS).
  * The SHA256 values are also the ones in each job's `provenance-<job>/model_ext.sha256` (for the job scripts:
    the copy in the provenance folder), so a fit can be matched to the code here.

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

## Extensions after 2026-10-27: when to use which

All of them wrap the unchanged `fit_cops.py` / `gpu_fit.py` / `fine_tune_gpu.py` in the running process, like the
grain code. Each file's header says how to run it; the `jobs/` scripts show complete commands.

| Need | Use | Module(s) | First used |
|---|---|---|---|
| hold some temperatures (Tmax, T75, Tmin) fixed | `grain_fit_fixed.py --fix T75=133 --fix Tmin=133` | `fixed_parameters.py` | 2026-11-04 (older form: `grain_fit_fixed_T.py` + `fixed_temperature.py`, job 121932) |
| hold any parameter fixed, e.g. a template amplitude at 1 | `grain_fit_fixed.py --fix c_<template>=1` | `fixed_parameters.py` | 2026-11-04 |
| use a fine tuning's activity correction as a fixed added term (same or another mesh) | `ft_template.py`, then the fit with `--activity-template <npz> <name> --fix c_<name>=1` | — | 2026-11-04 (20M option 2), 2026-11-05 (12M cycles) |
| start a fit from a fit + its fine-tuned grain | `restart_grain_start.py`, then `--restart-from <new npz> --template-start 0` | — | 2026-11-05 (cycles) |
| fit/fine-tuning cycles (correction templates accumulated) | `jobs/grain_ft_cycle*.sbatch` | the three above | 2026-11-05 |
| give the last N data points a larger weight | `grain_fit_weighted.py --weight-last N --weight-factor F [--weight-profile constant\|linear] [--weight-start S] [--weight-values W1:...:WN] [--weight-gauges ng,rg]` | `late_weights.py` (+ `fixed_parameters.py` for `--fix`) | 2026-11-06 (job 122236: last 5 x10; job 122248: last 6, linear 1 → 20) |
| carry a grain fit to another mesh and fine-tune there | `translate_grain.py`, `grain_translate_finetune.sbatch` | — | 2026-11-04 |
| refine elements of a shape model locally | `refine_mesh.py` (+ `seen_after_gap.py` for the 20M mask), `refined_bundle.py` | — | 2026-11-04 |
| scale the activity of chosen areas without refitting | `area_scaling.py`, `area_scaling.sbatch` | — | 2026-11-04 |
| rerun a finished grain fit exactly | `grain_repro.sbatch` | — | 2026-11-03 |
| evaluate finished fits on a fine time grid of the late window with the COPS time response (10-s running average, 10 s + DPU 5-readout average) and gauge offsets, without refitting | `late_time_grid.py` (observation file with a 1-s grid), `native/liouville_cops_cuda` (bases), `late_response_eval.py`; `jobs/late_response.sbatch` | — | 2026-11-07 (job 122289) |
| fit a strictly positive H2O activity A = exp(Y a + c T_cliff) (c ≥ 0) with a night/shadow ratio, an H2O flux floor F = max(q A, r q_max A_orig), a log-space Tikhonov prior (emission-weighted projection of ln A_orig), the late observations compared with the gauge time response (W10x5) and model-only 1-s grid rows | `positive_fit.py` (options before `--`, gpu_fit.py arguments after it); `jobs/positive_fit.sbatch` (LAMBDA, resume segments); base with the grid: `jobs/prepare_12m_lategrid.sbatch` (`late_time_grid.py --allow-missing-l2`), check `jobs/base_consistency_check.sbatch` | `positive_problem.py`, `positive_activity.py` + `kernels/pressure.cu`, `kernels/activity.cu` (env COPS_POS_ENABLE, COPS_POS_NIGHT_RATIO, COPS_POS_ACTIVITY=exp, COPS_POS_FLOOR_RATIO); tests `positive_activity_check.py`, `test_positive_fit.py` | 2026-11-07 (12M lambda scan, jobs 122307, 122310–122312, 122319, 122320) |
| per-facet fine tuning of a positive-activity fit (A_e = exp(s_e + delta_e), down to the floor, no cap, the fit's lambda, grain refitted), and fit / fine-tuning cycles with the change held in the next fit | `positive_finetune.py` (`--facets late_only --window-min 5`: only the facets seen exclusively in the last 5 min, earlier pressures unchanged; or `all_seen`); `jobs/positive_finetune.sbatch`; cycles: `jobs/positive_cycles.sh`; next fit: `positive_fit.py --start-restart`, `--held-exponent`, `--ratio 0.002` | `kernels/finetune.cu` (J^T u), `positive_problem.py`; test `test_positive_finetune.py` | 2026-11-08 (12M, ratio 0.002, six lambdas x 3 cycles; frozen code `model_ext_frozen/c1_20261108T000327Z` for the fits, `ft_20261108T032615Z` / `ft_20261108T034023Z` for the fine tunings) |
| refit the earlier clipped model (original-T cycle 6) on another base, with or without the gauge time response (diagnostic) | `clip_refit.py --variant D1|D2 --start-restart <c6 restart>`; `jobs/clip_refit.sbatch` (VARIANT D0|D1|D2) | `clip_problem.py` (uses positive_problem's W10x5 functions), `clip_compare.py`; test `test_clip_refit.py`, `jobs/clip_refit_bitwise.sbatch` | 2026-11-08 (D0 122880, D1 122865, D2 122866) |
| ratio 1e-5 with a CO2 flux floor, and an explicit start override | `positive_fit.py --ratio 1e-5` (sets COPS_POS_FLOOR_CO2=1: F_CO2 = max(kappa Q_CO2 ACO2, r qmax_C ACO2)), `--start-set kappa=1`; `jobs/positive_fit.sbatch` RATIO=1e-5, START_SET; `positive_cycles.sh --ratio 1e-5 --c1-start-set kappa=1` | `kernels/pressure.cu` (eval_fused_m_p5f/p14f), `positive_activity.py`, `positive_problem.py`; check `positive_activity_check.py` | 2026-11-09 (six chains 122951-123076; snapshot `r1e-5k_20261109T023540Z`) |
| one activity multiplier per effective-gravitational-slope bin (0-40, 40-50, ..., 80-90, 90-180 deg; 0-40 held at 1) on the 25 signed harmonics, floor 1e-5 A_orig, D2's setup without templates | `slope_mult_fit.py --start srccg` (user 2026-11-10: original srcCG coefficients, multipliers 1, kappa 1, D2's grain); `jobs/slope_mult_fit.sbatch` START, L_THREADS=32 (chunked, bitwise); bin bundle `slope_bins.py` (`jobs/slope_mult_cpu.sbatch STEP=bins`); figures `figures/figures_12m_slope_mult.sh` | `slope_mult.py` (GPU patch of the host-activity path, kernels unchanged), `slope_mult_problem.py` (D2 core of `clip_problem.py`); tests `test_slope_mult.py` (`jobs/slope_mult_cpu.sbatch STEP=tests`), `jobs/slope_mult_threads_test.sbatch` | 2026-11-10 (bins 123624; fit 123786; threads 123891, 123979) |
| Gaussian bumps x slope multipliers (option 3a), fit / fine-tuning cycles (dropped by the user 2026-11-11) | `slope_bumps.py` (bump bundle, NNLS start; CPU job 123981), `slope_mult_fit.py --bump-bundle` (BUMPS, START=nnls / chain), `slope_finetune.py` (`jobs/slope_finetune.sbatch`), chain `jobs/slope_cycles.sh` with guard jobs `jobs/slope_chain_guard.sbatch`; figures `figures/figures_12m_slope_bumps.sh` | `slope_mult.py`, `slope_mult_problem.py`; tests `jobs/slope_bumps_test.sbatch`, `jobs/slope_ft_diag.sbatch` (checks 6/7), `jobs/slope_bumps_threads_test.sbatch` (never submitted) | 2026-11-10 (chain 124008-124086, snapshots v6 `slopemul_20261110T201212Z`, v8 `slopemul_20261110T221713Z`) |
| non-negative activity by construction with the slope multipliers: softplus A = m_b s ln(1 + exp(L/s)), s = 1e-5 A_orig, or squared harmonics A = m_b (A_s (Y b)^2 + A_s sum ft_c*); D2's ft_c1..ft_c6 held at 1; fine tuning with no lower bound | `nonneg_fit.py --model softplus\|squared` (`jobs/nonneg_fit.sbatch`), `nonneg_finetune.py` (`jobs/nonneg_finetune.sbatch`), squared start `nonneg_start.py` (`jobs/nonneg_start.sbatch`); production as one job per chain (`jobs/chain_onejob.sh --kind nonneg`); figures `figures/figures_12m_nonneg.sh` | `nonneg_mult.py`, `nonneg_problem.py` (+ slope_mult helpers); tests `test_nonneg.py` (`jobs/nonneg_cpu.sbatch`), `jobs/nonneg_test.sbatch`, `jobs/nonneg_start_determinism.sbatch` | 2026-11-11 (softplus 124162, squared 124165; snapshot `nonneg_20261111T044920Z`) |
| activity and instantaneous H2O / CO2 source-rate maps of a final state for visualization (ParaView .vtu, Tecplot ASCII), flux exactly as the fitter computes it, sun angle and cast shadow for every facet at chosen UTC times | `jobs/nonneg_maps.sbatch` (CHAINS = chain directories, MODEL_EXT) -> `nonneg_maps.py` (+ `nonneg_maps_geom.cpp`, `nonneg_maps_shadow.cpp`, `nonneg_maps_pvcheck.py`, `test_nonneg_maps.py`, diagnostic `nonneg_maps_diag_dispatch.py`) | `nonneg_maps.py` | 2026-11-11/12 (124322 test; 124324, 124325 production) |
| softplus with the six free slope multipliers bounded below by 1 (user 2026-11-12: "does not allow the multiplier factor decrease below 1; start all of them with 1; they can grow") | `nonneg_fit.py --slope-mult-min 1` (`jobs/nonneg_fit.sbatch` SLOPE_MULT_MIN=1); chains `jobs/chain_onejob.sh --kind nonneg --model softplus --slope-mult-min 1` | `nonneg_problem.py` (bounds override in the existing hook; marker refuses a mixed resume); tests `test_nonneg.py` (multmin), `test_chain_onejob.py` (1c, 3d) | 2026-11-13 (5 cycles 124664, extended to 9 cycles 124772; snapshots `nonneg_20261113T035423Z`, `chain_onejob_20261113T041215Z`; FINDINGS.md A8) |
| per-facet contributions of a final state to the modelled NG / RG signals: P% source areas, source-area occurrence, slope-bin fractions | `jobs/nonneg_contrib.sbatch` (CHAIN, CYCLE, GEOM, MODEL_EXT[, CONTRIB_NAME]) -> `nonneg_contrib.py` (+ `test_nonneg_contrib.py`) | `nonneg_contrib.py` (frozen `nonneg_contrib_20261112T210012Z`) | 2026-11-12 (124549); bounded c5 124668, c9 124776 |
| alternative non-negative activity models (options A, B, D; user 2026-11-14 "a polynomial that is non-negative and would allow fit with more seen surface contributing"): sum of squares A = m_b (‖Lᵀy‖²/A_s + A_s ΣT), background A = m_b (A_bg A_orig + softplus), log-normal spread A = m_b (A_orig e^u + A_s ΣT) with the penalty f C_ref Var_w(u) | `nonneg_fit.py --model sos\|background\|logspread` (`jobs/nonneg_fit.sbatch` MODEL, SOS_RANK, BG_MIN, SPREAD_F, SPREAD_CREF), stage 1 `jobs/nonneg_alt_stage1.sh`, summary `jobs/nonneg_alt_summary.sbatch` | `nonneg_alt.py`, `nonneg_alt_summary.py` (+ nonneg_mult / nonneg_problem); tests `test_nonneg.py`, `jobs/nonneg_alt_test.sbatch` | 2026-11-14 (stage 1: B 125032-125046, A 125101-125105 with the LM start, D 125072-125086) |
| a whole fit / fine-tuning chain as ONE job (records bound once; same products and layout as the split chain) | `jobs/chain_onejob.sh` (dry run; `--submit`) -> `jobs/chain_onejob.sbatch` -> `chain_onejob.py`; KIND positive (default), slope (one cycle only), nonneg | the stage modules of the chosen kind; tests and the split-vs-one-job compare `test_chain_onejob.py` | 2026-11-10 (positive tests 124091-124095, bitwise); 2026-11-11 (nonneg tests 124143/124148, production 124162/124165; frozen `chain_onejob_20261111T060742Z`) |
| one fit split across one or two GPU nodes (MPI, one rank per GPU) | `jobs/mpi_fit_submit.sh 1\|2\|auto <spec>` -> `jobs/mpi_fit.sbatch` (SPEC files in `mpi_multinode/specs/`); probe `jobs/mpi_probe.sbatch` | `mpi_support.py` (per-rank wrapper); comparisons `mpi_multinode/compare_tests.py`, `compare_slope.py` | 2026-11-10 (tests 123623-123649 only; the production fits ran on one node) |
| which configuration (nodes, CPU threads, one job per chain) is fastest | `PERFORMANCE.md` (measured 2026-11-10/11) | — | — |
| retries without starting a stage behind a failed upstream (guard job per stage) | `jobs/slope_cycles.sh`, `jobs/positive_cycles.sh` (from 2026-11-10); guards `jobs/slope_chain_guard.sbatch`, `jobs/positive_chain_guard.sbatch`; layout test `jobs/test_chain_deps.sh` | — | 2026-11-10 (dummy jobs 124035-124058; 3a chain guards 124068-124086; positive guard not yet submitted) |
| non-negative activity from a LOCAL basis (option C): hat functions of a 150 m coarse triangulation, c >= 0, x the 7 slope-bin multipliers (>= 1), D2's templates on/off, smoothness lambda_C x^T Lg x; CPU fit from precomputed design matrices | `jobs/localnn_build.sbatch` (mesh, design, verify, start, tests), `jobs/localnn_fit.sbatch` (TEMPLATES, SMOOTH_FACTOR or SMOOTH_LAMBDA, NFEV, RUNROOT) | `localnn_mesh.py`, `localnn_design.py`, `localnn_problem.py`, `localnn_fit.py`, `localnn_direct.py`; tests `test_localnn.py` | 2026-11-14 (build 125004; timing test 124995; frozen `localnn_20261114T081912Z`) |

Notes:

* **Held parameters** need free temperatures (not `--ordered-temperature`, not `--fix-temperature`), and with held
  templates `--template-start 0` (the start then comes from the restart). `fixed_parameters.py` was fixed on
  2026-11-06: it evaluates the regularization term, its Jacobian and the diagnostic D on the full parameter vector
  (before, held template amplitudes inside the regularized block shifted it; a fit with 7 held parameters crashed,
  job 122118). With λ = 0, which every fit used, cost and parameters are unaffected; only the reported D of earlier
  fits with held templates is wrong.
* **Weights:** W multiplies a point's cost term (its residual and Jacobian row by √W); the crossing penalty uses the
  unweighted residuals; reported costs, rms and the restart file are unweighted, so they compare with other fits;
  the weighted objective is in `fit_manifest.json["model_extension"]["late_weights"]`. Without `--weight-last` the
  driver is the same as `grain_fit_fixed.py`.
  Profiles: `constant` (W = F on all N points), `linear` (from S = `--weight-start`, default 1, at the N-th point
  from the end to F at the last), or explicit `--weight-values` in time order. Example, "from 1 at point 6 from the end
  to 20 at the last point, linearly": `--weight-last 6 --weight-profile linear --weight-start 1 --weight-factor 20`
  (W = 1, 4.8, 8.6, 12.4, 16.2, 20). `--weight-start` and `--weight-values` were added after job 122236; with their
  defaults the weights are unchanged.
* **Positive-activity fits** (user decisions 2026-11-07): mixed precision only (FP64 records and `restart.cu` are refused while COPS_POS_* are active), so the fine-tuning tools (`fine_tune_gpu.py`, `grain_finetune.py`, `late_response_eval.py`, which use max(activity_linear, 0)) are **not valid** for these fits; the restart file's `activity_linear` holds the exponent s. The driver switches the start-up consistency check off (`--check-tolerance inf`), so check each new base once with `jobs/base_consistency_check.sbatch` (unchanged code). Each fit's `ACTIVITY_PARAMETERIZATION.txt` lists which of its products are valid. Grid rows (use 0) are model-only: `figures/positive_table.py` removes them before plotting.
* **The original temperature law** (used by `jobs/grain_ft_cycle_Torig.sbatch`) is the analytical law of
  `../../LiouvilleCOPS/src/Source.cpp` with its defaults Tmax 182.09922, T75 136.11149, Tmin 133 K; to use it in a
  grain fit, hold the three: `--fix Tmax=182.09922 --fix T75=136.11149 --fix Tmin=133`.

## Job scripts (`jobs/`, run copies in `$COPS_OUTPUT_ROOT/jobs/`)

| Script | What it ran | Jobs |
|---|---|---|
| `probe_gpu.sbatch`, `build_check.sbatch`, `hash_check.sbatch` | site probe, builds + self tests, hash check against passion (setup, handoff §3, §5, §9) | 120809, 120810/120814, 120867 |
| `submit_earlier.sh`, `submit_campaign.sh` | runner campaign submission (+ C16 set) with the site environment and provenance | 12M, 20M campaigns |
| `c16.sbatch` | C16 fit set (monotonic penalty + ordered T) | 120870, 120873 |
| `grain_ft_cycle.sbatch` | 12M fit/fine-tuning cycles 1–3 from C30 | 122080–122082 |
| `grain_ft_cycle_T130.sbatch` | cycle 4 with all T ≥ 130 K | 122096 |
| `grain_ft_cycle_Tmin130.sbatch` | cycles 4–6, Tmin = 130 K and Tmax = 231.98 K held | 122117, 122140, 122141 |
| `grain_ft_cycle_Torig.sbatch` | cycles 4–6, the original temperature law held (figures: `figures/figures_12m_cycles_Torig.sh`) | 122219–122221 |
| `grain_late_weights.sbatch` | fit with the last N points weighted (WLAST, WFACTOR, WPROFILE, WSTART, WVALUES), from cycle 3 (START=c3) or the original-temperature cycle 6 (START=torig_c6) | 122236, 122248, 122250 |
| `late_response.sbatch` | late-window forward check: 1-s grid, bases for gauge offsets 0/1.5/3 m along +Z (3 GPUs), time-response windows, five fits + fine tunings (figures: `figures/figures_12m_late_response.sh`) | 122289 |
| `prepare_12m_lategrid.sbatch` | 12M fitting base with the 1-s late grid and gauge offset +1.5 m along +Z (no record store) | 122291 |
| `base_consistency_check.sbatch` | a base against the unchanged GPU forward model (`gpu_fit.py --data-only`, tolerance 1e-6) | 122304 |
| `positive_activity_check.sbatch`, `positive_unit_12m.sbatch`, `positive_fit_test.sbatch` | positive-activity kernel checks, unit tests, smoke test (λ 0, 4 evaluations, Jacobian check) | 122298, —, 122303 |
| `positive_fit.sbatch` | positive-activity fit, LAMBDA required, max_nfev 4000, resume segments; from 2026-11-08 also RATIO, CYCLE, START_RESTART, GRAIN_*, HELD, MODEL_EXT (figures: `figures/figures_12m_positive_scan.sh`) | 122307, 122310–122312, 122319, 122320; cycle-1 fits 122449–122454 |
| `positive_finetune.sbatch`, `positive_cycles.sh`, `positive_held_test.sbatch` | per-facet fine tuning (FACETS required); fit / fine-tuning chains for one lambda with retry attempts; held-term GPU tests (figures: `figures/figures_12m_positive_cycles.sh`) | chains 122596–122703 (λ 0 again 122706–122723); tests 122445, 122446, 122448–122477, 122594, 122595, 122704, 122705 |
| `positive_cycles.sh` (2026-11-10 version), `positive_chain_guard.sbatch` | guard-job layout (user 2026-11-10 "Apply the guard fix"): attempts behind one gate job, a CPU guard job per stage that exits 0 only for a complete stage; settings and products unchanged | not yet used for a submission (dry runs only) |
| `slope_mult_fit.sbatch`, `slope_mult_cpu.sbatch`, `slope_mult_threads_test.sbatch` | slope-multiplier fit (START, L_THREADS, BUMPS, HELD, CHAIN; one node); CPU steps (STEP=bins \| tests); threaded / chunked L test (figures: `figures/figures_12m_slope_mult.sh`) | fit 123786 (earlier script version); bins 123624; CPU tests 123625, 123628, 123636, 123892, 123978, 123983, 124012, 124062; threads 123891, 123979 |
| `slope_cycles.sh`, `slope_finetune.sbatch`, `slope_chain_guard.sbatch` | bump + slope-multiplier chain (option 3a, 3 cycles, guard job per stage); its late-facet fine tuning; the guard (figures: `figures/figures_12m_slope_bumps.sh`) | c1 fit 124008 (by hand), 124064–124086 |
| `slope_bumps_test.sbatch`, `slope_ft_diag.sbatch`, `slope_bumps_threads_test.sbatch`, `test_chain_deps.sh` | 3a GPU tests (STEP=fits \| ft); read-only diagnostic of the fine-tuning checks 6/7; 64-thread timing (never submitted); dummy-job test of the chain dependency layouts | 123984, 123985, 124013, 124063; 124034; —; 124035–124058 |
| `nonneg_fit.sbatch`, `nonneg_finetune.sbatch`, `nonneg_start.sbatch` | softplus / squared fit (MODEL required); its late-facet fine tuning (no lower bound); the squared start (CPU) (figures: `figures/figures_12m_nonneg.sh`) | split tests 124139–124142, 124144–124147; inside nonneg_test.sbatch 124130–124133; start 124129 (test), 124138 (production) |
| `nonneg_cpu.sbatch`, `nonneg_test.sbatch`, `nonneg_start_determinism.sbatch` | CPU tests; GPU tests (STEP=fits \| ft, one MODEL per job); determinism of the squared start per CPU type | 124128; 124130–124133; 124149–124153 |
| `nonneg_maps.sbatch` | activity and source-rate maps of one final state (CPU, normal128cores): tests, geometry build and validation against the base's kernel records, maps (.vtu + Tecplot), verify, ParaView read-back; frozen code `model_ext_frozen/nonneg_maps_20261112T021623Z`; 2026-11-13 output name `MAPS_NAME` (default `maps`) and the `EXPECT_CYCLES` guard | 124318-124322 (development and 3-cycle test maps), 124324 (softplus), 124325 (squared), 124667 (multipliers >= 1, c5), 124775 (the same chain, c9: `maps_c9`) |
| `nonneg_contrib.sbatch` | per-facet contributions of a final state to the modelled NG / RG signals (P% source areas, occurrence, slope-bin fractions; CPU, frozen nonneg_contrib code, NPY_DISABLE_CPU_FEATURES for the bitwise state check); `CYCLE`, output name `CONTRIB_NAME` (2026-11-13, default `contributions`) | 124549 (softplus 5 cycles), 124668 (multipliers >= 1, c5), 124776 (c9: `contributions_c9`) |
| `nonneg_alt_test.sbatch`, `nonneg_alt_stage1.sh`, `nonneg_alt_summary.sbatch`, `nonneg_multmin_test.sbatch` | options A, B, D: GPU test per model (independent start cost, 5 evaluations with the full Jacobian check); stage-1 launcher (3 attempts, guard, summary; run roots `nonneg/<model>/stage1_<setting>_<UTC>/`); CPU summary of one fit; the softplus multiplier-bound GPU test | 124998, 124999, 125028; B 125032-125046, A 125101-125105, D 125072-125086; 125030 (test); 124659 |
| `chain_onejob.sh`, `chain_onejob.sbatch` | ONE job per fit / fine-tuning chain (KIND positive \| slope \| nonneg) | positive tests 123990, 124091–124095 (earlier script); nonneg tests 124143, 124148; production 124162 (softplus), 124165 (squared)  2026-11-15 (v7): nonneg sos \| background \| logspread tokens and `--c1-from` (cycle 1's fit adopted); GPU test 131984-131988; stage-2 chains 131989-131997  2026-11-16 (v7.1): fine-tuning `--ft-delta-*` (FT_DELTA_*); GPU test 135365/135366; reruns 135367-135378; (v7.2) --c1-from with 1 cycle, gate 135861/135862, reruns 135863-135874; 2026-11-17 (v7.3): free-space guard `--min-free-gb` / MIN_FREE_GB (default 20 GiB, exit 5), D 136016-136018, A 136019-136021 |
| `nonneg_stage2_chains.sh` | stage 2 of options A / B / D (user 2026-11-15): one one-job chain per model (background bg-min 0.01, logspread f 0.1 with the fixed-level u0, sos rank 3), 5 cycles, `--slope-mult-min 1`, cycle 1 adopted from the stage-1c fit; chain roots `nonneg/<model>/stage2_<setting>_<UTC>/` with README.md / HOW-TO-RERUN.md; dry run by default | 131989-131991 (B), 131992-131994 (D), 131995-131997 (A), gated afterok:131988 |
| `nonneg_stage2r_chains.sh` | stage-2 RERUNS with the regularized late-facet fine tuning (user 2026-11-16: `--delta-penalty-f 10 --delta-cref 4.2688066575e-03 --delta-cap 4.605170185988092`): softplus reference (start D2, first), background, logspread, sos (cycle 1 adopted from stage1c); 5 cycles, `--slope-mult-min 1`; roots `nonneg/<model>/stage2r_<setting>_<UTC>/`; dry run by default | first submission 135367-135378 cancelled with its failed gate 135365 (v7.1 defect, fixed in v7.2); resubmitted 135863-135865 (softplus), 135866-135868 (B), 135869-135871 (D), 135872-135874 (A), gated afterok:135862 |
| `pub_maps_figures.sbatch`, `pub_maps_v2.sbatch`, `pub_source_areas.sbatch` | full / version-2 figure sets (GPU; inputs from the environment) | 124492-124516, 124565-124595, 124677-124683, 124722 |
| `cops_movies_prep.sbatch`, `cops_movies_render.sbatch`, `cops_movies_encode.sbatch`, `cops_movies_qc.sbatch`, `cops_movies_diag.sbatch` | movies: preparation and validation, render on all GPUs, compose / encode, MP4 check, dark-surface diagnosis; 2026-11-13 also the fixed camera (`CAMERA=fixed`, `FROM_WORK` + step `copy`, step `validate_fixed`), the touchdown data (step `touchdown`), the flag (`TOUCHDOWN`) and `MOVIES_PY` | 124615-124635, 124707-124735, 124786-124788, 124863-124869 |
| `ng_source_map.sbatch`, `pan_candidates.sbatch`, `cops_movies_pan_stills.sbatch` | 2026-11-14 panorama movie: per-facet NG source sums (late and all NG observations; CPU), candidate final camera poses (1 GPU), candidate stills (plan `--pan-to`, mid-path and last frame; 1 GPU per pose); `cops_movies_prep.sbatch` `PAN_TO`/`PAN_START`/`PAN_END` and `cops_movies_encode.sbatch` `NO_SUB_MARK` (defaults unchanged) | 124913, 124915, 124918, 124919, 124920; production 124930, 124931, 124936, 124937 |
| `cops_movies_prep.sbatch` (`SC_BLOCK_FROM`, step `scrows`), `cops_movies_render.sbatch` (`SC_BLOCK_FROM`), `fov_share.sbatch` (`HEMI_SC`, `HPREP`, `HMOVIES`, `DECK`) | 2026-11-14 half-space movies: the fit's RG / NG front half-spaces minus the directions hidden by the spacecraft from the ROSINA deck point; the shares of the modelled signal the movies do not colour (defaults unchanged) | 124942-124945; tests 124940 |
| `frames_base.sbatch`, `frame_contrib.sbatch`; `cops_movies_render.sbatch` (`COLOR_CONTRIB`, `COLOR_RANGE`), `cops_movies_encode.sbatch` (`HOLD_SECONDS`) | 2026-11-14 contribution movies: fitting-base kernels at the 1959 frame times (the lategrid base's executable and boundary.in), the fitted model's contribution per unit area per frame, verification, colour range; render coloured by it; 15 s hold at the end (defaults unchanged) | 124961, 124969, 124972, 124975-124977; tests 124974 |
| `contrib_phase_data.sbatch`, `pub_contrib_density.sbatch` | 2026-11-14 full set of the c9 state: contribution per unit area at snapshot rows and phase means (data, CPU) and their maps (GPU, short); region shares by lobe / distance (`regions_figures.py`, CPU) | 124979, 124990 |
| `options_activity_maps.sbatch` | 2026-11-14 activity maps of the best stage-1 setting of each non-negative option and the bounded softplus c9 (`figures/pub_options_activity_maps.py`; 1 GPU, ~3 min) | 125205; final 131968 |
| `full_set_submit.sh`, `finish_full_set.sh`, `d2_tools.sbatch`; `nonneg_maps.sbatch` (cycle detection skips c<N>.incomplete_*), `contrib_phase_data.sbatch` (optional FT / FIT / SBINS / ADAPTER) | 2026-11-17 full figure sets of the six cases of 12M_stage2_reg_best (softplus, A, B, D regularized, best C, D2): the c9 set's figure jobs parameterized, the last CPU steps, the D2 export / contributions / verify / maps (defaults unchanged) | maps 136152, 136166, 136219, C 136200, D2 136217; figures 136154-136266 (see each set's HOW-TO-RERUN.md) |
| `cliff_effect_data.sbatch`, `cliff_effect_figures.sbatch` | 2026-11-17 slope at equal illumination (user: does the source rate depend on the slope angle at the same illumination?): data per case, figures and answer tables of all six stage-2 cases (CPU) | data 136324-136329, figures 136332, document 136339 |
| `illumination_numbers.sbatch`, `illumination_maps.sbatch` | 2026-11-18 illumination maps at the production-rate map times (user: "i dont see the effect of the shadow or night side decrease of the productions -- is it there? generate illumination maps ..."): numbers and cross-case check (CPU), then the maps of one time in every production-rate layout plus the side-by-side figure (1 GPU, `TIME=`) | numbers 136359, figures 136369 (10:33:10), 136370 (10:39:26); tests 136360, 136367 |
| `source_rate_q002.sbatch` | 2026-11-18 source-rate maps of the runs fitted with the new H2O law (q_min = 0.02 q_max; user: "prepare the surface source rate maps for completed runs for inspection of the source rate"): pub_maps_v2 activity + H2O rate at both map times and the new-vs-old-law inspection figures of one run (GPU; submitted by `figures/source_rate_q002_submit.py`) | D2_q002_warm 136737, 136741 (136740 superseded) |
| `stage2_q002_figures.sbatch` | 2026-11-19 comparison set of the new-law stage-2 round (user: "prepare a full set of comparison plots for the last round of fitting"): stage2_q002_figures.py on CPU | 137248 (test 137230) |
| `plausibility_q002_metrics.sbatch` | 2026-11-21 activity-spread and energy-limit numbers of the ten new-law fits for the plausibility ranking (user: "analyze of the fitting cases and range them by physicall plasuablity ..."); CPU, 1-2 min; output `12M/q002_plausibility/` |
| `q002_top5_figures.sbatch`, `q002_top5_docx.sbatch` | 2026-11-21 figures, numbers and the Word document on the five most plausible new-law fits (user: "prepare a summary word doc for the 5 models in the table ..."); CPU; output `12M/q002_top5_models/` (the docx job waited on the CPU limit and its commands ran on the login node) |
| `q002_seen_gauges.sbatch`, `q002_seen_gauges_gpu.sbatch`, `q002_slope_activity_pass2.sbatch`, `q002_touchdown_maps.sbatch`, `q002_slope_matched.sbatch` | 2026-11-21 per-gauge visibility pass (16 CPUs, 2 min; the gpu-partition variant was only a queue test, cancelled), slope / illumination statistics, touchdown maps, matched slope metric; output `12M/q002_top5_models/slope_activity/` |
| `q002_threshold_options.sbatch` | 2026-11-21 switched-off share of the ten new-law fits with the existing active-area threshold (1 % of the srcCG median) and around it; output `12M/q002_decisions_20261121/` |
| `q002_active_area.sbatch` | 2026-11-21 switched-off share of the ten new-law fits by the user's decision 1b (activity without the imposed floor > 1 % of the srcCG median; FINDINGS A21) |
| `q002_all_fits_summary_docx.sbatch`, `q002_correction_notices.sbatch` | 2026-11-21 decision 2b: rebuild of the all-fits document; notices on the older documents and figure sets |
| `localnn_slopetest_code_test.sbatch`, `localnn_slope_test.sbatch`, `localnn_slope_test_submit.sh` | 2026-11-21 decision 3c: code test (137892, no fit); the eight slope-effect control fits (plan `12M/q002_decisions_20261121/PLAN_3c.md`; submitted 2026-11-21 22:40 CST with the lower bound 0.01, jobs 137898-137913) |
| `q002_slopetest_report.sbatch` (`figures/q002_slopetest_report.py`) | 2026-11-21 decision 3c: results table of the control fits (costs, rms, switched-off share by 1b, ordering) |
| `q002_slopetest_state_check.sbatch` (`figures/q002_slopetest_state_check.py`) | 2026-11-21 decision 3c: a control fit's saved state in the reference's own problem (old frozen code, direct facet sums, paths between the two solutions); no fit |
| `localnn_multistart.sbatch`, `localnn_multistart_submit.sh`, `localnn_multistart_preflight.sbatch` (`localnn_multistart.py`) | 2026-11-22 user "multi-start when ready": the option C reference problem (six multipliers >= 1) refitted with the unchanged frozen fitter from each better end state of the 3c controls with all multipliers >= 1, up to five fits |
| `q002_stationarity_check.sbatch` (`figures/q002_stationarity_check.py`) | 2026-11-22 is an end state a minimum or an early ftol stop: bounded Gauss-Newton step and the objective along it; use it before calling an option C fit converged; with `SOLVER=nnls` (exact; the default trf solve can hit its iteration limit: then inconclusive) |
| `localnn_ac_{fit,chain,check,prepare}.sbatch`, `localnn_ac_submit.sh` (`localnn_ac.py`) | 2026-11-22 follow-up of the slope test (user "as recommended: a, c"): six multipliers >= 0.01 from each of the five solutions; the all-held model from other starts and as a chain (`localnn_chain.sbatch` with the opt-in HOLD_ALL_MULT=1). Plan `12M/q002_decisions_20261121/PLAN_ac.md`; c run 2026-11-22 (FINDINGS A23), a withdrawn |
| `localnn_cliff_{fit,check,code_test}.sbatch`, `localnn_cliff_submit.sh` (`localnn_cliff.py`, `test_localnn_cliff.py`) | 2026-11-22 cliff test (user "prepare the cliff test"): the all-held model with the cliff classes held at 0.01, the moderate slopes as the control (`localnn_fit.py --slope-test --hold-mult BIN=V`). Plan `12M/q002_decisions_20261121/PLAN_cliff.md`; waiting for the user's OK |
| `q002_ac_preview.sbatch` (`figures/q002_ac_preview.py`) | 2026-11-22 what the follow-up a and c can achieve (evaluations only): one bounded step with the bound 0.01 at each of the five solutions; the late facets' share of the late pressure |
| `plausibility_q002_v2.sbatch` (`figures/plausibility_q002_v2.py`) | 2026-11-23 plausibility numbers of all 17 new-law fits with one method (Q(t), seen share, energy limit, decision 1b); ranking in `12M/q002_plausibility/README.md` |
| `full_set_allheld_spec.py`, specs `data/{source_rate_q002_runs,full_sets_q002}_20261123.json` | 2026-11-23 full figure set of the all-held candidate (`12M/full_set_q002_C_f10_allheld/`, pipeline of 3.18 with `full_set_q002.py --spec`, `source_rate_q002_submit.py --runs`, `full_set_q002_finish.sbatch` SPEC). Known limit: `cliff_importance_figures.py` labels slope bins by position, so a partly held state (the cliff-test fits) would be mislabelled |
| `full_set_q002_finish.sbatch`; `full_set_submit.sh` (MEM_CAP, SKIP), `finish_full_set.sh` (BASE, REF_DIR), `pub_maps_v2.sbatch` (BASE, MODEL_TEXT), `pub_source_areas.sbatch`, `pub_maps_figures.sbatch` (BASE, D2REF), `contrib_phase_data.sbatch`, `pub_contrib_density.sbatch` (BASE, MCODE) | 2026-11-20 full figure sets of the ten new-law states (user: "prepare the full set of plots for the last round of fits with the fixed surface source distribution"); the overrides keep every default | figure jobs 137598-137668, finish 137669-137678, D2 contributions 137493 |
| `fov_share.sbatch`, `source_area_compare.sbatch` | signal share outside the IK FOVs; source areas of two fits (CPU) | 124731, 124727; c9: 124789, 124784, 124785 |
| `mpi_fit.sbatch`, `mpi_fit_submit.sh`, `mpi_probe.sbatch` | one fit over 1 or 2 nodes with MPI (SPEC in `mpi_multinode/specs/`, NODES=1\|2\|auto); submit helper; MPI probe | tests 123629–123649 (earlier script versions; the current one not submitted); probe 123623 |
| `grain_rough.sbatch`, `prepare_12m_rough.sbatch` | 12M rough-surface base and the C30 fit on it | 122033, 122034 |
| `prepare_44m.sbatch`, `grain_44m_c30.sbatch` | 44M preparation; C30 with the 12M conditions on 44M | 121996/122020, 122078 |
| `refine_12m.sbatch`, `refine_20m.sbatch` | local mesh refinement tests | 122009, 122011 |
| `localnn_build.sbatch`, `localnn_fit.sbatch` | option C (local non-negative basis): CPU build (coarse mesh, design matrices, D2 / c9 verification, starts, tests) and one CPU fit (wall-time stop / exit 75, RUNROOT resume, direct facet-sum check) | build 125004; timing test 124995 (dev grid); start-rule stage 1 125014-125016 (cancelled, superseded); pilot 125071; anchored stage 1 125107, 125113, 125119 (+ continuations 125108, 125114, 125120) |
| `localnn_gpu_eval.sbatch`, `localnn_tools.sbatch`, `localnn_report.sbatch` | option C after a fit: GPU standard evaluation (gpu002, D2's path, 1 h), the frozen maps / contribution tools through the localnn adapter (CPU), the item-5 numbers and run-root README.md / HOW-TO-RERUN.md (CPU) | tests 125017, 125052, 125019-125031; start-rule stage 1 125053-125064 (cancelled); anchored stage 1 125109-125112, 125115-125118, 125121-125124 |
| `localnn_chain.sbatch`, `localnn_chain_submit.sh`, `localnn_chain_post.sbatch`, `localnn_late.sbatch`; tests `localnn_stage2_tests.sh`, `localnn_ft_verify_softplus.sh`, `localnn_ft_verify_c9_perturb.sh`, `localnn_gpu_solver_{cpuref,bitwise}.sh`, `localnn_gpu_solver_verify.sbatch` | option C stage 2: 5-cycle fit / CPU fine-tuning chains (held term, checkpoints, HANDOVER), their submission and post-fit tools, the late kernels; tests and verifications; FT_DIR in `localnn_gpu_eval.sbatch` / `localnn_tools.sbatch` | late kernels 132008; tests 131980, 132009; verification 131977-131979; chains 132010-132045; GPU solver 132047-132049 |
| `localnn_q002_preflight.sbatch`, `localnn_q002_preflight_compare.py`, `localnn_q002_d2check.sh`, `localnn_q002_posts_maps.sh` | option C for the H2O night law q_min = 0.02 q_max: preflight of the new D2 (selection, templates, base, start files, test_localnn), the design built against the new D2 + verification, the five posts with maps (geometry validated against base_q002, contribution code nonneg_contrib_20261118T192126Z) | preflight 136557, d2check 136558; posts 136694-136698 |

The job scripts in this folder itself (`grain*.sbatch`, `area_scaling.sbatch`) ran from `$COPS_OUTPUT_ROOT/model_ext/`.

## Figure scripts (`figures/`, run copies in `$COPS_OUTPUT_ROOT/figures_cluster/scripts/`)

`model_obs_figures.py` draws every comparison (model vs COPS observations, .eps + .png 300 dpi);
`finetune_table.py` turns a fine tuning into a table it can plot; `how_to_rerun.py` writes HOW-TO-RERUN.md;
`snapshot.sh` copies a campaign's results and restart files next to its plots. `rerun_commands.py` writes the "exact commands" part of HOW-TO-RERUN.md (from 2026-11-06): every command of a fit and its fine tuning with all arguments, the job's exports, the AMPS commit of the run, and for each module and job script the git commit holding an identical copy here; job scripts from `grain_late_weights.sbatch` on also record their resolved commands in `provenance-<job>/commands.sh`. Each `figures_*.sh` makes one figure
folder under `$COPS_OUTPUT_ROOT/figures_cluster/<RES>/` (its header names it) and writes HOW-TO-RERUN.md there; from
2026-11-06 on (`figures_12m_late_weights.sh`) also `code/` with the code that ran and the AMPS commit of each fit.
They run on the login node (no GPU), from `$APP` with `.venv/bin/python`. `late_response_figures.py` + `figures_12m_late_response.sh <job>` (2026-11-07) draw a `late_response.sbatch` job
(model with and without the gauge time response, per gauge offset) into `12M/model_vs_observations_late_response/`.
`figures_12m_positive_scan.sh` (2026-11-07) draws the six positive-activity λ fits into `12M/model_vs_observations_positive_lambda_scan/` (`positive_table.py` removes the model-only grid rows); `scan_summary.py` (run copy in `$COPS_OUTPUT_ROOT/analysis/positive_lambda_scan_12m/`) writes the scan table and the L-curve.
`figures_12m_positive_cycles.sh` + `positive_cycles_table.py` (2026-11-08) draw the fit / fine-tuning cycles into `12M/model_vs_observations_positive_cycles/` (final state per lambda; first scan vs cycle-1 fit vs final; cycles_table.txt).
`RATIO=1e-5 figures_12m_positive_cycles.sh` (2026-11-09) draws the ratio-1e-5 chains into `12M/model_vs_observations_positive_cycles_r1e-5/`.
`figures_12m_clip_vs_positive.sh` (2026-11-09) draws the last round against the clipped-model refit (D0/D1/D2) into `12M/model_vs_observations_clip_vs_positive/`.
`snapshot_positive.sh` (2026-11-09) snapshots the positive-activity work and the clipped refits (restart files hard-linked) into `12M/snapshot_positive_<UTC>/`.
`shadow_classes.py` (2026-11-10; run copy in `$COPS_OUTPUT_ROOT/analysis/shadow_zero_activity_12m/`) classifies the seen facets of a fit by illumination (sunlit / cast shadow / night) over the fit's records, for zero-activity and active facets.
`analysis/overhang_illumination.py` + `analysis/overhang_illumination.sbatch` (2026-11-12; run copies in `$COPS_OUTPUT_ROOT/analysis/overhang_illumination_12m/`, job 124578) split the seen facets of the softplus 5-cycle state by illumination (sunlit / night / cast shadow) per effective-slope bin at the map times (FINDINGS.md A1).
`figures_12m_slope_mult.sh` (2026-11-10) draws the slope-multiplier fit 123786 against D2 and the best positive fit into `12M/model_vs_observations_slope_mult/`.
`figures_12m_slope_bumps.sh` (2026-11-10) draws the bump + slope-multiplier chain (3a, 124008-124086) into `12M/model_vs_observations_slope_bumps/`.
`figures_12m_model_families.sh` (2026-11-11) draws the best fit of each of the four families of 2026-11-08 ... 11-10 (D2, positive, slope multipliers, bumps) into `12M/model_vs_observations_model_families/`.
`figures_12m_nonneg.sh` (2026-11-11; extended 19:26 UTC with every chain stage, cycle 1 vs final and all families) draws the softplus and squared chains (124162, 124165) into `12M/model_vs_observations_nonneg/`.
`figures_12m_nonneg5.sh` (2026-11-12) draws the 5-cycle softplus and squared chains (124310, 124313; cycle-1 fit and the state after each cycle, and both vs D2 and the best positive fit) into `12M/model_vs_observations_nonneg_5cycles/`, with MAPS.txt pointing to the maps of the final states.
`pub_ng_rg_altitude.py` (2026-11-12) draws the publication figure of COPS NG / RG and the softplus model's NG / RG vs altitude (5-cycle final state, legend on the right, gap in the NG data hatched) into `12M/publication_softplus/`.
`pub_maps_figures.py` (2026-11-12; data-driven captions 2026-11-13) draws the full set's activity / production-rate maps with fixed ranges and is the helper of the version-2 scripts; `pub_maps_v2_figures.py` (version 2: seen-only colour, per-panel log classes, `--model-name`), `pub_source_area_figures.py` (P% source areas) and `cliff_importance_figures.py` (slope-bin use; `--groups`, `--fit-dir`, `--maps-dir`) write `12M/full_set_softplus_5cycles_v2/` and `12M/full_set_softplus_smm1_5cycles/`; `collect_captions.py` collects the captions into CAPTIONS.md; `model_data_captions.py` (2026-11-13) writes the model-data captions of a full set.
`figures_12m_nonneg5_smm1.sh` (2026-11-13) draws the softplus refit with slope multipliers >= 1 (chain 124664) vs the unbounded softplus and D2 into `12M/model_vs_observations_softplus_smm1/`; `source_area_compare.py` compares the source areas and slope-bin shares of two fits (`full_set_softplus_smm1_5cycles/comparison_vs_unbounded/`).
`cops_movies.py` + `cops_movies_geom.cpp` (2026-11-12; observable-part masks, spacecraft shadow and the dark-surface diagnosis 2026-11-13) make the one-minute movies of the descent seen from RG / NG (`12M/movies_softplus_5cycles/`, `12M/movies_softplus_smm1_5cycles/`); `fov_share.py` gives the share of the modelled signal from outside the IK fields of view. Since 2026-11-13 (second set) `cops_movies.py` also has the fixed camera (`plan --camera fixed`: frozen in the body frame at the COPS position and attitude of 100 m altitude), the gauges' line of sight from the true gauge position (depth map + exact edge refinement; `validate_fixed`), `touchdown` (the ground-track impact point, its visibility, SPICE check) and `compose --touchdown` (a waving chequered finish flag): `12M/movies_softplus_smm1_9cycles/` (gauge view) and `12M/movies_softplus_smm1_9cycles_fixedcam/`; and (2026-11-13/14) circular fields of view (`cone:<deg>` masks), `contact`, `prep --extend-to-utc`, `sc --hga-ck`, `plan --camera detach` and the detach render with a spacecraft model (bus, wings, HGA), its shadow and the FOV cones (`cops_movies_geom hgageom`): `12M/movies_softplus_smm1_9cycles_detach200/`; `fov_share.py --cone-rg/--cone-ng` gives the share outside circular fields of view. 2026-11-14: `plan --camera detach --pan-to <pose.json>` (the camera moves after the detach to a final pose: smootherstep in the frame index, arc around the pivot, quaternion slerp, focal length geometric) and `compose --no-sub-mark` (no sub-spacecraft +); `ng_source_map.py` (per-facet NG source sums), `ng_source_regions.py` (lobes, shares by lobe / distance / azimuth, relief, map), `ng_lobe_share_vs_altitude.py` (big-lobe share per late observation) and `pan_candidates.py` (candidate final poses): `12M/movies_softplus_smm1_9cycles_detach300_pan/` (the panorama movies with the user's pose B, held from frame 1800; analysis and candidate stills). 2026-11-14: `render/validate_fixed --sc-block-from rosina-deck` (fields of view minus the spacecraft model, exact segment test), the `hemisphere` mask in the detach render (line of sight by cube faces), `scrows`, `fov_share.py --hemi-sc`: `12M/movies_softplus_smm1_9cycles_detach300_pan_hemi/`. 2026-11-14: `render --color-contrib/--color-range` (colour = the fitted model's contribution per unit surface area, `frame_contrib.py` on the frames base `frames_base.sbatch`, range `contrib_range.py`), `compose --hold-seconds`: `12M/movies_softplus_smm1_9cycles_detach300_pan_contrib/`.
`figures_12m_nonneg9_smm1.sh` (2026-11-13) draws the same bounded chain after 9 cycles (jobs 124664, 124772) vs its c5 state, the unbounded softplus and D2, and the cycles c1-c5 / c5-c9, into `12M/model_vs_observations_softplus_smm1_9cycles/`.
`figures_12m_nonneg_options_stage1.sh` (2026-11-14) compares the stage-1 fits of the non-negative activity options (A sum of squares, B background, C local basis, D log spread) with the COPS data, the bounded softplus c9 and D2 in `12M/model_vs_observations_nonneg_options_stage1/` (per option all settings, the best of each, CPU); `options_table.py` (summary.txt / .csv / best.json from the runs' summaries), `options_captions.py` (captions), `options_howto.py` (HOW-TO-RERUN.md, README.txt), `pub_options_activity_maps.py` (activity maps of the best settings and c9; job `jobs/options_activity_maps.sbatch`, 1 GPU). v2 (2026-11-14 evening, `..._stage1_v2/`): wrapper variables `DIAG` (convergence verdicts of the CPU diagnosis, `figures/data/stage1_convergence_20261114.json` -> `options_table.py --diag`), `RESID=1` (`options_residual_figure.py`: ln(model/obs) vs altitude, NG and RG) and `ACT_FROM` (copy unchanged activity maps); defaults unchanged (tested). Final (2026-11-15, `..._stage1_final/`): `RUNS` (candidate list file, `figures/data/stage1_final_runs_20261115.json`: the latest completed fit per setting), `TEXTS` (`figures/data/stage1_final_texts_20261115.json`, `options_captions.py --texts`), `CONV_TERM=1` (convergence from the terminations); continuation rows carry the stage-1 costs they started from (prev_*); activity maps job 131968; defaults unchanged (tested).
`stage2_figures.py` + `stage2_howto.py` (2026-11-17) draw the final stage-2 comparison (regularized and unregularized chains: model/data figures, residuals, cost table with terminations and bounds; spec `figures/data/stage2_final_spec_20261117.json`; CPU, no maps, stops below 25 GB free) into `12M/model_vs_observations_nonneg_options_stage2_final/`. `kappa_table.py` (2026-11-17) lists kappa (value, at the lower bound) of every stage-2 chain, the stage-1 fits, the softplus chains and D2 (`<folder>/kappa_check/`).
`full_set_facts.py` + `full_set_docs.py` (2026-11-17) compute the Facts of one full figure set from its data (cycles, multipliers, kappa, A = 0 shares, Q_H2O and its seen / unseen split, slope-bin production and pressure shares, P50 / P90 source areas, lobe / shell shares) and write its README.txt, HOW-TO-RERUN.md and code/ (spec `figures/data/full_sets_stage2_20261117.json`); `pub_maps_figures.py` (`bin_edges_from_meta`, `multiplier_bounds`), `pub_maps_v2_figures.py`, `cliff_importance_figures.py`, `contrib_phase_data.py` (`--state-adapter`) and `model_data_captions.py` accept the option-C and D2 exports (2026-11-17; the nonneg sets unchanged).
`fullsets_comparison_docx.py` (2026-11-17) writes the Word document that compares the six stage-2 full figure sets plot by plot (`12M/full_sets_stage2_comparison/12M_stage2_full_sets_comparison.docx`: per plot type the six PNGs in one grid, numbers from captions / checks / manifests / facts / summary.csv with their sources in `numbers_used.json`, Appendix A source paths) and validates it (`validate`: counts, page fit, layout previews); runs in the separate venv `$COPS_OUTPUT_ROOT/tools/docx_venv` (python-docx, Pillow), not the project `.venv`; CPU job, ~20 s.
`cliff_effect_data.py` + `cliff_effect_figures.py` + `cliff_effect_docx.py` (2026-11-17) answer whether the source rate depends on the slope angle at the same illumination, for the six stage-2 final states (`12M/cliff_effect_stage2/`): the data script recomputes each kernel record with the fit's arithmetic (equal to `record_contributions` bitwise) and accumulates per illumination class (lit SZA 10-deg classes, shadowed, night) x slope bin the area, A, F, the observation-weight split (1/r^2, cos theta, gauge weight, visible fraction) and the matched-cell and decomposition sums on the gauge-seen surface (snapshot rows, phase means) and, from the maps' .vtu, on the whole surface; the figure script draws P1 (F and A at equal illumination), P2 (observation weight), P3 (matched cells), P4 (exact decomposition) and the SZA-controlled ratio tables (spec `figures/data/cliff_effect_cases_20261117.json`: labels, caveats); the docx script packs each plot of the six cases into one cluster with summaries (imports `fullsets_comparison_docx.py`; docx_venv).
`illumination_maps_figures.py` (2026-11-18) answers whether the night-side / cast-shadow decrease is in the production-rate maps (`12M/illumination_maps_103310_103926/`): `numbers` (CPU job `jobs/illumination_numbers.sbatch`) verifies the reference state's fields against its maps (.vtu, bitwise), checks F_H2O = q_H2O A bitwise and q_min/q_max = 0.2, computes per time the area per illumination class, the mean q/q_max, the night / shadow shares of the non-zero production area and of the production (whole and seen surface) and the class changes between the times, and checks that cos_SZA, illumination and the seen mask are bitwise identical in the maps of the six full-set cases; `figures` (GPU job `jobs/illumination_maps.sbatch`, one time) draws the illumination class, cos SZA (0.1 classes) and the source-rate factor q/q_max (linear and log10) for the whole surface and the seen facets in every layout of the full sets' production-rate figures (imports `pub_maps_v2_figures.py` / `pub_maps_figures.py` unchanged; layouts copied) and the side-by-side figure (production rate | q/q_max on the same facets | q/q_max on the whole surface; classes checked equal to the production-rate manifest); `docs` writes README.txt (answer, numbers), CAPTIONS.md (`collect_captions.py`), HOW-TO-RERUN.md and code/.
`source_rate_q002_inspect.py` + `source_rate_q002_submit.py` + `source_rate_q002_docs.py` + `source_rate_q002_cluster.py` (2026-11-18) make the source-rate figure sets of the runs fitted with the new H2O law (`12M/source_rate_q002/<run>/`): the submit script takes every input of a run from its maps manifest and the runs spec `figures/data/source_rate_q002_runs_20261118.json` (maps, old-law counterpart, model name, the model's own caption text per model); the inspection figure shows the new-law H2O production rate next to the old-law counterpart's (common log classes) and the new q/q_max, with Q_H2O total / seen, night + shadow shares and active area (the old-law state is read from its maps' .vtu: the old maps code's law guard refuses to run under the new reference law); docs writes README / CAPTIONS / HOW-TO-RERUN / code; cluster composes the production-rate panels of all completed runs. `pub_maps_figures.py` / `pub_maps_v2_figures.py` (2026-11-18): `--model-text` (MODEL_TEXT, optional {multipliers}) replaces the softplus description in the map captions (default unchanged).
`stage2_q002_figures.py` + `stage2_q002_docs.py` (2026-11-19) draw the comparison set of the new-law (H2O q_min = 0.02 q_max) stage-2 round (`12M/model_vs_observations_nonneg_options_stage2_q002/`; spec `figures/data/stage2_q002_spec_20261119.json` with each run's old-law counterpart from `stage2_final_spec_20261117.json`): best of each option and the five C chains vs COPS, new vs old law per run (model/data and residuals), the data cost per cycle per chain, summary.txt / .csv with kappa, Q_H2O total / seen at 10:33:10 (source_rate_q002 inspection files, the C report / maps for C f10 1e-4) and the cost change; the readers of `stage2_figures.py` and `kappa_table.py` imported unchanged.
`full_set_q002.py` (2026-11-20) makes the full figure sets of the ten new-law states (`12M/full_set_q002_<key>/`, spec `figures/data/full_sets_q002_20261120.json`): `spec` (inputs from each state's maps manifest), `prepare` (publication NG/RG figure folders `publication_q002_<key>/`; activity/, production_rate/, inspection/ copied from `source_rate_q002/<run>/` with sha256 checks), `submit` (`jobs/full_set_submit.sh` with BASE = base_q002, the new-law snapshots, the model's own `--model-text`, MEM_CAP 24G, SKIP activity,h2o), `finish` (`jobs/full_set_q002_finish.sbatch`), `docs` (`full_set_docs.py` with the new optional spec fields intro / reference / pub / written / extra_dirs).
`q002_comparison_docx.py` (2026-11-19) writes the Word document that clusters the same plot of the ten new-law runs (D2_q002_warm and the nine regularized stage-2 chains, H2O q_min = 0.02 q_max) with per-plot summaries (`12M/stage2_q002_comparison/12M_stage2_q002_comparison.docx`: overview of the multi-run figures, then per plot type the runs in a fixed order with placeholders for runs without the plot; numbers with sources in `numbers_used.json`; Appendix A SHA-256 of every source PNG) and validates it (captions state the law, counts, page fit, previews); derived from `fullsets_comparison_docx.py`, separate venv `$COPS_OUTPUT_ROOT/tools/docx_venv`; CPU job, ~20 s.
`q002_all_fits_summary_docx.py` (2026-11-20) writes ONE Word document summarizing all fits of the new-law round (`12M/q002_all_fits_summary/12M_q002_all_fits_summary.docx`; supersedes the q002 comparison document for this round): runs table new vs old, A18 summary with the full-set facts, the multi-run and per-run comparison-set figures, one cluster per full-set plot type for the ten runs, the source-rate inspection, an old-law appendix; `check` (inputs, gaps, extractor dry run), `build`, `validate`; observations from the images in `figures/data/q002_all_fits_observations_20261120.json`, viewing sheets by `q002_all_fits_sheets.py`; separate venv `$COPS_OUTPUT_ROOT/tools/docx_venv`; CPU job, ~1 min.
`plausibility_q002_metrics.py` (2026-11-21) reads `A_final`, the seen mask and the facet areas of each new-law run and writes the area fractions with zero / graded activity (relative to the trimmed mean), the effective area, f50/f90, the Gini coefficient, and the share of Q_H2O at 10:33:10 above the energy limit S0 / (r_h^2 dH) (`12M/q002_plausibility/metrics.{json,csv}`; ranking in that folder's README.md and FINDINGS A19).
`q002_top5_figures.py`, `q002_top5_sources_figures.py`, `q002_top5_docx.py`, `q002_top5_content.py` (2026-11-21) make `12M/q002_top5_models/12M_q002_top5_models.docx`: the five most plausible new-law fits (ranks 1-5 of `q002_plausibility`): summary, two tables, model/data, residual, criteria and convergence figures, each fit's activity and production-rate maps, and the methods (forward model and assumptions, what the fit minimizes, the activity models A, B, C and D2 and how they differ numerically, results, plausibility; section 4 on the comet and its sources: ground track, source rates, contribution of areas, slopes and cliffs, robust vs model-dependent, from `q002_top5_sources_figures.py`); the method text lives in `q002_top5_content.py` with its sources (OUTPUT_VISUALIZATION 3.20).
`q002_seen_gauges.py`, `q002_slope_activity.py`, `q002_slope_matched.py` and the checks `q002_shutdown_check.py`, `q002_visibility_check.py`, `q002_test_seen_gauges.py`, `q002_spike_check.py` (2026-11-21): visibility per gauge, source rate around the touchdown, activity against the angle to gravity with illumination and location matched, switched-off surface (OUTPUT_VISUALIZATION 3.21, FINDINGS A20).

`q002_active_area.py` (2026-11-21) computes the switched-off share by decision 1b (B's imposed floor not counted) for the ten new-law fits and checks it against the c5 reports; `q002_correction_notices.py` put the notices of decision 2b on the older documents and sets. Use the first one, not `A == 0` and not `q002_shutdown_check.py` (A < k A_u), to say whether a fit switches the surface off (FINDINGS A21).

## MPI helpers (`mpi_multinode/`, run copies in `$COPS_OUTPUT_ROOT/mpi_multinode/`)

Test and comparison scripts of the multi-node work (user decisions 2026-11-10 "one fit split across both nodes", "make
it also runnable on one node"; no physics). The launcher is `jobs/mpi_fit.sbatch` (+ `jobs/mpi_fit_submit.sh`), the
per-rank wrapper `mpi_support.py`; results in `PERFORMANCE.md`.

| File | Role | Jobs |
|---|---|---|
| `mpi_multinode/mpi_probe.py` | MPI probe: rank -> GPU map, ping-pong, Bcast of 96 MB, the transport's request/reply round (run by `jobs/mpi_probe.sbatch`) | 123623 (`mpi_multinode/probe-123623/`) |
| `mpi_multinode/specs/clip_D2.sh` | SPEC for `jobs/mpi_fit.sbatch`: clip_refit.py --variant D2 as job 122866 (test) | 123629–123632, 123638, 123645, 123647, 123649 |
| `mpi_multinode/specs/slope_mult.sh` | SPEC for `jobs/mpi_fit.sbatch`: the slope-multiplier fit over MPI (host activity) | 123646, 123648 (bitwise = 123635) |
| `mpi_multinode/compare_tests.py` | comparison of the launcher's test fits (local / 8 / 16 ranks: costs, SHA-256 of x / residuals / Jacobian, timing, memory) | `clip_refit_tests/mpi_20261110/compare_*` |
| `mpi_multinode/compare_slope.py` | slope-multiplier MPI tests vs the one-node test 123635 | `slope_mult_tests/mpi_*_srccg_20261110T094013Z/compare_vs_123635.json` |

## H2O night-side fix (2026-11-18): what changed here and what it means for old code

User decision 2026-11-18: "(b) -- use 0.02; x10 was not intended; use the same law but with decreased activity". The
H2O night / cast-shadow production is q_min = 0.02 q_max (as srcCG `Comet::GetTotalProductionRateBjornNASTRAN`); it was
0.2 q_max (`10.0*0.02`) in Source.cpp, fit_cops.py, cops_cuda/evaluator.py and build_physical_activity_templates.py, and
**every fit, map and figure up to 2026-11-18 used 0.2**. CO2 unchanged (0.1, r_h^-2; differs from srcCG too, not asked).

* `nonneg_maps.py`: `H2O_QMIN_RATIO = 0.02` / `CO2_QMIN_RATIO = 0.1` are the one model_ext place for the ratio
  (`q_factors`, the law_guard line); tools that need it read it from there (`figures/illumination_maps_figures.py`:
  `QRATIO` = the maps' own `nonneg_maps.H2O_QMIN_RATIO`). `late_response_eval.py` has its own named constant
  (`H2O_QMIN_RATIO`, it imports only fit_cops). `positive_activity.py` `REFERENCE_QMIN_RATIO` H2O 0.02;
  `positive_activity_check.py` reads it. Captions of `figures/pub_maps_figures.py`, `pub_maps_v2_figures.py` take
  q_min / q_max from the maps' law. `test_nonneg_maps.py`: Source.cpp statement 0.02 (36 of 36).
* `slope_mult.py` / `positive_activity.py`: `ORIGINAL_SHA256["evaluator.py"]` = the fixed evaluator.py (its only change
  is `HostConstants.qmin`); without it `nonneg_mult.install` (every nonneg fit and fine tuning) and `positive_activity`
  refuse on the fixed tree.
* New snapshots: `model_ext_frozen/nonneg_20261118T072729Z` (fits / fine tunings), `nonneg_maps_20261118T072729Z`,
  `nonneg_contrib_20261118T072729Z`.
* **Old frozen snapshots** take fit_cops and cops_cuda from the LIVE tree (`$APP`, `$APP/../LiouvilleCOPS/tools`), not
  from a copy: (a) nonneg fit / fine-tuning snapshots before 20261118 refuse (`slope_mult._check_sources`: evaluator.py
  changed); (b) nonneg_maps / nonneg_contrib snapshots before 20261118 refuse wherever `law_guard` runs (maps,
  contributions, `localnn_design.py`, `pub_maps_figures.load_fields`); (c) tools that call `nonneg_maps.source_rates` /
  `q_factors` without the guard (`nonneg_alt_summary.py`, `nonneg_ft_summary.py`, `localnn_report.py`,
  `localnn_chain_summary.py`, `figures/cliff_effect_data.py`, `figures/frame_contrib.py`, `full_set_facts.py`) keep the
  0.2 of the snapshot they load; (d) code that uses the GPU evaluator or fit_cops directly (grain fits, clip_refit / D2,
  the runner, gpu_fit) silently uses 0.02; (e) option C (`localnn_*`) evaluates through the stored design matrices
  (`grid_150m/design/design.npz`, `late_kernels_150m.npz`, built with 0.2 by the frozen contribution code) and has no law
  check: a C fit on them keeps 0.2 silently until the design is rebuilt with `nonneg_contrib_20261118T072729Z`.
* **Fits on an existing base:** `fit_cops.main` (also through gpu_fit) first checks its model at the base's own
  parameters against the base's `COPS_model.dat`, which the C++ wrote with 0.2: on the 12M lategrid base that check now
  gives max rel 28 % (NG) / 25 % (RG) > 1e-6 and stops every fit (old law: 4.8e-11 / 9.9e-10). The base must be rebuilt
  with the rebuilt `native/liouville_cops_cuda` (user decision), see `tests/nightfix_diagnostic_20261118/README.txt`.
* Diagnostic / checks (new): `nightfix_diag.py` (fixed-A re-evaluation of a final state with both laws; `--base-run`:
  the base consistency check), `nightfix_cpu_gpu_check.py` (fit_cops vs the GPU evaluator on real 12M rows at the B29
  state, VALIDATION.md budgets), jobs `jobs/nightfix_*.sbatch` (the build job's copy is `jobs/nightfix_binaries.sbatch`:
  the repository's `*build*` ignore pattern would hide `nightfix_build.sbatch`) / `jobs/nightfix_diag_submit_all*.sh`,
  `jobs/nightfix_getq_probe.cpp`, `nightfix_make_table.py`. Results: `tests/nightfix_diagnostic_20261118/`,
  `tests/nightfix_20261118/`; FINDINGS A18.

* **base_q002 (2026-11-18):** the lategrid base with the new-law COPS_model.dat (hard-linked kernels; no
  run_manifest.json: fits use fit_cops' `--allow-legacy-base`). Maps / contributions: `nonneg_maps.py validate-geometry
  --shape-manifest base_q002/cpp_run_manifest.json --shape-manifest-must-equal base/run_manifest.json` (snapshots
  `nonneg_maps_20261118T192126Z`, `nonneg_contrib_20261118T192126Z`; `jobs/nonneg_maps.sbatch` SHAPE_MANIFEST /
  SHAPE_REFERENCE; `jobs/d2_tools.sbatch` D2DIR / D2_BASE); validated geometry `maps_geometry/job_136691`. Tools:
  `nightfix_base_q002.py`, `jobs/q002_*.sbatch`, `jobs/clip_refit_q002.sbatch` (the D2 refit `clip_refit/D2_q002_warm`).

## Physics findings: `FINDINGS.md`

The physics results of the 12M fits with numbers, evidence (jobs, files) and caveats, for use in publications, plus method facts and open questions: `FINDINGS.md` (kept current; user request 2026-11-12).

## Figures: `OUTPUT_VISUALIZATION.md`

How every figure of the 12M fits is made (model/data comparisons, publication figure, ParaView/Tecplot maps, templates, effective slope, ground tracks, activity and production-rate maps, contribution maps) and the recipe to redo them all after a model rerun: `OUTPUT_VISUALIZATION.md` (kept current with every new or changed figure; user request 2026-11-12).

## How to run the model (fit, fine tuning, cycles, figures)

Added 2026-11-11. All commands are as run on Voyager (paths in "Site-specific values" above). GPU work runs only
inside Slurm jobs; the job scripts write their provenance (commands, git, environment, code copy, module/pip lists)
next to the products.

**Environment** (every job script sets it; for a shell on a compute node):

```bash
export AMPS=/data001/highresolu_vt/AMPS APP=$AMPS/srcCG/LiouvilleCOPS_cuda
export COPS_OUTPUT_ROOT=/data001/highresolu_vt/cops_output
export T=/data001/highresolu_vt/cops_cluster_transfer_20261026 COPS_DATA_ROOT=$T/rosetta_data_end_of_mission_cops
export CUDA_PATH=/usr/local/cuda-12.9 CUDA_HOME=/usr/local/cuda-12.9 CUPY_CACHE_DIR=$COPS_OUTPUT_ROOT/.cupy
export PATH=$CUDA_HOME/bin:$PATH LD_LIBRARY_PATH=$CUDA_HOME/lib64${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1
cd $APP     # interpreter: $APP/.venv/bin/python
```

**What a run needs:**
- a fitting base (`runs/<run>/base`, built by the runner or by a job such as `jobs/prepare_12m_lategrid.sbatch`)
  and its template bundle (`runs/<run>/bundle/restart_templates.npz`);
- a start (a fit's `fit_restart.npz`);
- the code: a **frozen snapshot** of `$COPS_OUTPUT_ROOT/model_ext` for every production run. Copy the files the job
  needs into `model_ext_frozen/<name>_<UTC>/`, write `SHA256SUMS`, make it read-only, and pass it as `MODEL_EXT=`.
  The job copies the snapshot into its provenance folder and runs that copy.

**Three ways to run, from the smallest to the largest unit:**

1. **One stage.** A fit job script (`jobs/positive_fit.sbatch`, `jobs/clip_refit.sbatch`,
   `jobs/slope_mult_fit.sbatch`, `jobs/nonneg_fit.sbatch`) or a fine-tuning job script (`jobs/positive_finetune.sbatch`,
   `jobs/slope_finetune.sbatch`, `jobs/nonneg_finetune.sbatch`), configured with `sbatch --export=ALL,VAR=value,...`.
   The variables are listed in each script's header.
   - A fit stops before its time limit at its best state and exits 75; resubmitting the same command continues it as
     the next segment (`segNN/`, `final` link when complete).
   - A fine tuning runs validation checks 1-8 in its first segment and refuses to optimize if one fails.
2. **A split chain of cycles,** each stage its own job: `jobs/positive_cycles.sh` or `jobs/slope_cycles.sh`. Run it
   without `--submit` first (dry run), then with it. A guard job after every stage starts the next stage only if the
   stage really completed. Each stage has retries that continue after a time-limit stop.
3. **A whole chain in ONE job** (records bound once, bitwise equal to the split chain; about 12 min saved per 3-cycle
   chain): `jobs/chain_onejob.sh --kind positive|nonneg ...` (dry run, then `--submit`), which runs
   `jobs/chain_onejob.sbatch` -> `chain_onejob.py`. A chain that hits the 8-h limit continues in the next attempt and
   skips the completed stages. The bump/slope kind (`--kind slope`) runs one-cycle chains only.

**Speed:** one node, local transport, chunked L step with 32 threads (`L_THREADS=32` / `--l-threads 32`; bitwise
equal to 1 thread). Two nodes (MPI) pay off only for the clipped models with device activity. Details:
`PERFORMANCE.md`.

**Figures:** each `figures/figures_12m_*.sh` (login node, no GPU) writes `figures_cluster/12M/model_vs_observations_*/`
with the figures, the tables, `HOW-TO-RERUN.md` (every stage's exact command, AMPS commit, code commits) and `code/`.

**Bitwise reproducibility.** It holds on the same CPU type at any thread count. Across CPU types the CPU steps that go
through numpy's OpenBLAS differ at rounding level, because OpenBLAS (DYNAMIC_ARCH) picks its kernels by CPU: node023
(EPYC 7763, `Zen`) vs node032 (EPYC 9554, `Cooperlake`) gave 8e-15 in the squared start (jobs 124149-124153).

## Softplus and squared-harmonics models (2026-11-10/11)

**Why** (user, 2026-11-10): "make step back to clipped D2 with 6 templates; propose a global decomposition of the
activity map that is non-negative by construction". Both models replace D2's clip max(0, ·) and its additive
`cliff_ge_60deg` template with an activity that cannot be negative.

**Models** (user decisions 2026-11-10/11, `/data001/highresolu_vt/AMPS/handoff.md` §10). H2O activity per facet e:

- **softplus:** A_e = m_b(e) · s_e ln(1 + exp(L_e / s_e)) · exp(h_e), with
  - L = Σ a_lm Y_lm (25 signed harmonics, degree ≤ 4) + A_s Σ_k T_k (D2's six ft templates ft_c1..ft_c6, held at 1);
  - s_e = 1e-5 · A_orig,e (A_orig = the original srcCG map). Where A_orig,e = 0 (0.71 % of the facets) A_e = max(0, L_e).
  - D2 is the limit s → 0.
- **squared harmonics:** A_e = m_b(e) · (A_s (Σ b_lm Y_lm)² + A_s Σ_k T_k) · exp(h_e). The templates stay outside the
  square.
- **Both:**
  - m_b: one multiplier per effective-gravitational-slope bin [0,40) (held at 1), [40,50), [50,60), [60,70), [70,80),
    [80,90), [90,180] deg (≥ 0; bundle `slope_bins.npz`);
  - no floor, no `cliff_ge_60deg`;
  - h_e: the accumulated per-facet fine-tuning term (0 in cycle 1);
  - everything else as D2: lategrid base (gauge +1.5 m along +Z), W10x5 gauge response from 10:33:10 UTC, original
    temperature law held (Tmax 182.09922, T75 136.11149, Tmin 133 K), κ free, grain + NG = RG crossing (weight 10),
    λ = 0, max_nfev 4000, `--tbounds 30 400`;
  - 34 free parameters: 25 harmonics, 6 multipliers, κ, grain A and B.
- **Starts (cycle 1, user "Both from D2"):**
  - softplus: D2's harmonic coefficients;
  - squared: b by area-weighted least squares so that the activity at m = 1 matches D2's max(0, field)
    (`nonneg_start.py`);
  - multipliers 1, κ and grain from D2.
- **Fine tuning (each cycle):** a factor exp(δ_e) on the facets seen only in the last 5 min (late_only), no lower bound
  and no cap, grain refitted, λ = 0, MAXNFEV 500, validation checks 1-8. Facets with A_e = 0 cannot be tuned. The
  next cycle's fit holds h = the accumulated δ. 3 cycles.

**Code:**
- modules: `nonneg_mult.py` (model, CPU stand-in, GPU patch), `nonneg_problem.py` (fit hooks), `nonneg_fit.py`
  (fit driver), `nonneg_start.py` (squared start), `nonneg_finetune.py`, `test_nonneg.py`;
- jobs: `jobs/nonneg_{fit,finetune,start,cpu,test}.sbatch`, `jobs/chain_onejob.{sh,sbatch}` (kind nonneg);
- snapshots used: `model_ext_frozen/nonneg_20261111T044920Z` (model) and `model_ext_frozen/chain_onejob_20261111T060742Z`
  (driver).

**Results** (2026-11-11, one job per model):

| model | final data cost | rms ln NG / RG | fits stop | κ |
|---|---|---|---|---|
| softplus (job 124162, gpu001, 50 min) | 4.3697e-03 | 0.0710 / 0.0608 | ftol (1517, 134, 291 evaluations) | 3e-36 (bound 0) |
| squared (job 124165, gpu002, 3 h 05 min) | 6.9398e-03 | 0.0902 / 0.0758 | max_nfev 4000 in every cycle | 7e-36 (bound 0) |
| D2, for reference | 4.8864e-03 | 0.0778 / 0.0610 | max_nfev 4000 | 2e-36 |

- **softplus:**
  - With s = 1e-5 A_orig and D2's very negative field, it underflows to exactly 0 on most of the seen surface. At the
    start that is 94.7 % of the seen area, as in D2; in the final fit 96.4 %.
  - Four multipliers end at 0; 60-70 deg is ×26.3 and 90-180 deg ×4.27.
  - 61 % of the seen production comes from slopes > 90 deg.
  - Only 73 of the 271 late facets are tunable.
- **squared:**
  - No facet is at 0.
  - The multipliers are 50-60 deg ×16.7 and 90-180 deg ×3.48; the others ≤ 8e-4.
  - 80 % of the seen production comes from slopes > 90 deg.

Figures: `figures_cluster/12M/model_vs_observations_nonneg/`.

**5-cycle runs and maps** (user 2026-11-11: "repeat softplus and squared runs adding 2 extra cycles of fitting and fine tuning; save the activity and source rate maps for visualization"): fresh runs with `--cycles 5` (jobs 124310 softplus, 124313 squared); cycles 1-3 bitwise equal to the 3-cycle runs. Final data cost: softplus 4.3668e-03 (every fit ftol; cycles 4-5 -0.06 %), squared 6.6491e-03 (every fit max_nfev 4000; cycles 4-5 -4.2 %, not converged). Maps of the final states (`<chain>/maps/`: `<model>_flux.pvd` + 13 `<model>_flux_<UTC>.vtu` with cos SZA, illumination, F_H2O and F_CO2 in m^-2 s^-1 at 00:00 ... 10:00, 10:33:10, 10:39:26 UTC; `<model>_static.vtu` with A_H2O, kappa ACO2, seen, slope bin, area; `<model>_maps.dat` for Tecplot 360; README.txt, maps_manifest.json): `sbatch --export=ALL,CHAINS=<chain>,MODEL_EXT=$COPS_OUTPUT_ROOT/model_ext_frozen/nonneg_maps_20261112T021623Z jobs/nonneg_maps.sbatch` (the chain's last cycle: c<N>/ft/final + c<N>/fit/final); the exact commands are in `figures_cluster/12M/model_vs_observations_nonneg_5cycles/HOW-TO-RERUN.md`. Figures: `figures_cluster/12M/model_vs_observations_nonneg_5cycles/`.

### Running them from scratch

Each step's products are the next step's inputs. `$J=$COPS_OUTPUT_ROOT/jobs`,
`$R=$COPS_OUTPUT_ROOT/runs/shap8_012m_lategrid_dz1.5`.

**0. Earlier inputs, built before 2026-11-07:**
- the 12M base `runs/shap8_012m_flat_earlier/base` and its bundle and gravity cache: the runner campaign
  (`../tools/cluster/README.md`, `../quick_guide.md` §4e3);
- the original-temperature fit / fine-tuning cycles 1-6, which give D2's start (c6 `fit_restart.npz`) and the
  templates ft_c1..ft_c6: `figures_cluster/12M/model_vs_observations_grain_ft_cycles_Torig/HOW-TO-RERUN.md` and
  `jobs/grain_ft_cycle_Torig.sbatch`.

**1. Lategrid base** (the 1150 observations plus a 1-s model grid in the late window, gauge +1.5 m; GPU, about 172 GB;
job of 2026-11-07):

```bash
sbatch --export=ALL,STOP_AFTER=base $J/prepare_12m_lategrid.sbatch     # -> $R/{late_grid.h,late_grid_points.txt,base,bundle}
```

**2. D2** (clipped refit of the c6 model on the lategrid base, W10x5; job 122866). The exact commands are in
`figures_cluster/12M/model_vs_observations_clip_vs_positive/HOW-TO-RERUN.md`:

```bash
sbatch --nodelist=gpu002 --export=ALL,VARIANT=D2,MODEL_EXT=$COPS_OUTPUT_ROOT/model_ext_frozen/clip_20261108T234013Z $J/clip_refit.sbatch
                                                                     # -> $R/clip_refit/D2/final
```

**3. Slope bins** (CPU; job 123624):

```bash
sbatch --export=ALL,STEP=bins,MODEL_EXT=$COPS_OUTPUT_ROOT/model_ext_frozen/slopemul_20261110T091946Z $J/slope_mult_cpu.sbatch
                                                                     # -> $R/slope_mult/bundle/slope_bins.npz
```

**4. Freeze the code.** Copy the files listed in `model_ext_frozen/nonneg_20261111T044920Z/SHA256SUMS` from
`$COPS_OUTPUT_ROOT/model_ext` into `model_ext_frozen/nonneg_<UTC>/`, and `chain_onejob.py` and `test_chain_onejob.py`
into `model_ext_frozen/chain_onejob_<UTC>/`. Write `SHA256SUMS` in both and make them read-only (`chmod -R a-w`).
Below, `NX` and `CX` are these two folders.

**5. Squared start** (CPU, about 1 min; job 124138):

```bash
sbatch --export=ALL,MODEL_EXT=$NX $J/nonneg_start.sbatch     # -> $R/nonneg/bundle/squared_start.npz (+ .json, provenance)
```

**6. Tests** (recommended before production; jobs 124128-124133):
- CPU: `sbatch --export=ALL,MODEL_EXT=$NX $J/nonneg_cpu.sbatch`;
- GPU, one node per model: `sbatch --nodelist=<node> --export=ALL,MODEL_EXT=$NX,TAG=<tag>,MODEL=<model>,STEP=fits
  $J/nonneg_test.sbatch`, then `STEP=ft`. They check the start cost against an independent evaluation, the
  softplus → clip limit, the full Jacobian, fine tuning with 32 threads vs 1 bitwise, and cycle 2 starting at the
  fine-tuning cost.

**7. The two chains, one job per model, concurrently** (jobs 124162 and 124165; 3 cycles; attempts 2 and 3 only
continue after a time-limit stop):

```bash
cd $J
./chain_onejob.sh --kind nonneg --model softplus --chain-root $R/nonneg/softplus/<name> --model-ext $NX --chain-code $CX \
  --cycles 3 --nfev 4000 --check 1 --ft-maxnfev 500 --ft-facets late_only --ft-window-min 5 --l-threads 32 \
  --node gpu001 --time 08:00:00 --attempts 3                # dry run; add --submit
./chain_onejob.sh --kind nonneg --model squared --chain-root $R/nonneg/squared/<name> --model-ext $NX --chain-code $CX \
  --start-file $R/nonneg/bundle/squared_start.npz \
  --cycles 3 --nfev 4000 --check 1 --ft-maxnfev 500 --ft-facets late_only --ft-window-min 5 --l-threads 32 \
  --node gpu002 --time 08:00:00 --attempts 3                # dry run; add --submit
```

Products: `<chain-root>/c<k>/{fit,ft}/` (`segNN/`, `final`, `provenance-<job>-<stage>/commands.sh` reruns that stage
alone) and `<chain-root>/provenance-<job>/`.

**Split alternative (one job per stage, same results bitwise):**

```bash
# cycle 1 fit
sbatch --nodelist=gpu001 --export=ALL,MODEL=softplus,START=d2,NFEV=4000,CHECK=1,L_THREADS=32,FITDIR=<chain>/c1/fit,MODEL_EXT=$NX $J/nonneg_fit.sbatch
#   (squared: MODEL=squared,START=lsq,START_FILE=$R/nonneg/bundle/squared_start.npz)
# its fine tuning
sbatch --nodelist=gpu001 --export=ALL,FIT=<chain>/c1/fit/final,OUT=<chain>/c1/ft,FACETS=late_only,WINDOW_MIN=5,MAXNFEV=500,L_THREADS=32,MODEL_EXT=$NX $J/nonneg_finetune.sbatch
# cycle k > 1 fit: from the previous cycle
sbatch ... --export=ALL,MODEL=<model>,START=chain,START_RESTART=<chain>/c<k-1>/fit/final/fit_restart.npz,HELD=<chain>/c<k-1>/ft/finetune_state.npz,GRAIN_FROM=<chain>/c<k-1>/ft/finetune_state.npz,NFEV=4000,CHECK=1,L_THREADS=32,FITDIR=<chain>/c<k>/fit,MODEL_EXT=$NX $J/nonneg_fit.sbatch
```

Submit each stage only after the previous one is complete. The one-job chain runs exactly these commands, token for
token; its stage commands are in each stage's `provenance-<job>-<stage>/commands.sh`.

**8. Figures** (login node):

```bash
$COPS_OUTPUT_ROOT/figures_cluster/scripts/figures_12m_nonneg.sh   # -> figures_cluster/12M/model_vs_observations_nonneg/
```

## Option C: local non-negative basis (2026-11-14)

**Why** (user 2026-11-14, relayed by the lead): the clipped / softplus models (25 harmonics of degree <= 4 x 7 slope-bin
multipliers) switch off ~97 % of the seen surface, because a global low-degree basis can only make local sources by
clipping; the user wants non-negative alternatives that let more of the seen surface contribute. Options A, B, D are
`nonneg_problem.py --model sos|background|logspread` (another agent); option C is this section.

**Model** (per facet e, slope bin b(e)):

- A_e = m_b(e) (sum_j c_j phi_j(e) + TS_e), c_j >= 0; phi_j = piecewise-linear hat functions of a coarse triangulation
  of the 12M surface (three per facet, phi >= 0, sum_j phi_j = 1 at every facet);
- TS = A_s sum_k ft_ck (D2's six fine-tuning templates held at 1) with `--templates on`, 0 with `--templates off`
  (user 2026-11-14: on, added as is because the sum is >= 0 on every facet; the six bundles are >= 0 everywhere, 265 facets each, so the activity is >= 0
  without a wrapper; the products apply max(0, .) and report 0 clipped facets);
- m_b: the 7 slope bins of `slope_mult/bundle/slope_bins.npz`, m_[0,40) = 1 held, six free >= 1 (start 1), no upper
  bound (as the bounded softplus, user 2026-11-12);
- everything else = D2: lategrid base (1551 rows), D2's selection (934 NG, 1150 RG), W10x5 for the 95 observations from
  10:33:10 UTC (`positive_problem.LateLayout`), grain at the time stamp, crossing weight 10, original T law held, kappa
  >= 0, grain A, B >= 0, ln-residual cost of fit_cops (0.5 sum r^2), scipy TRF (bounds, x_scale jac, exact trust-region
  solver, ftol = xtol = gtol = 1e-10);
- smoothness: lambda_C x_c^T Lg x_c (x_c = c / A_s, A_s = 2.2228299497402522e18; Lg = the P1 stiffness matrix =
  cotangent-weighted graph Laplacian of the coarse mesh, x^T Lg x = integral of |grad sum c phi|^2 over the coarse
  surface), as residual rows sqrt(2 lambda_C) R x_c, R^T R = Lg;
- start (`--start d2`): c = the area-weighted NNLS projection of max(0, D2 field) (minus TS with templates on), m = 1,
  kappa and grain from D2;
- no fine-tuning term h_e yet (stage 1).

**Grid** (`localnn_mesh.py`, `runs/shap8_012m_lategrid_dz1.5/nonneg/localnn/grid_150m`, job 125004): n = round(2 S /
(sqrt 3 h^2)) = 2745 nodes for S = 53.495 km^2, h = 150 m; seeds by area-weighted farthest-point sampling, 6 Lloyd
iterations of the geodesic Voronoi partition of the 12M facets (multi-source Dijkstra on the facet-adjacency graph);
coarse triangles = the dual (5486 triangles, 8229 edges, every edge in exactly 2 triangles; Euler characteristic 2: the
four small tunnels of the genus-4 12M surface are bridged at 150 m); coarse edge length median 141 m (p10 120, p90
171 m); region area 13.6-37.3 x 10^3 m^2. Hat functions: barycentric coordinates of each facet centroid projected on the
plane of the best coarse triangle of its node (2-ring as fallback); 97.1 % of the area inside a triangle, 2.84 % outside
every candidate (negative coordinates set to 0 and renormalized; median -0.012); max |sum phi - 1| = 2.2e-16. Continuity:
L1 jump of phi across fine edges median 0.022 (p99 0.069); jumps > 0.5 on 3627 of 18.5 M fine edges (0.085 % of the
area, mostly steep folded terrain). 2207 nodes have support on seen facets; the other 538 are set by the smoothness term.

**Design matrices** (`localnn_design.py`): for fixed T the pressures are linear in A and kappa, so G^(b)[row, node] =
sum over the row's records with bin b of g_r phi_j(e_r) (g_r = the H2O contribution at A = 1 from the frozen
`nonneg_contrib.record_contributions`: FP32 records, FP64 arithmetic), template and CO2 columns; one pass over the 4.74e9
records, 94 s on 64 cores; design route vs direct facet sums <= 9.3e-15 (every row).

**Verification** (`localnn_fit.py verify`, `test_localnn.py`; job 125004): D2's final activity: data cost
4.8864368867405e-03 (reported 4.886436886740517e-03, rel 8.9e-16), crossing 0, COPS_model.dat pNG / pRG of all 1551 rows
and late_window_model.dat (surface, grain): 0 printed differences; softplus c9 final state: 4.268806657481031e-03
(reported 4.268806657481027e-03, rel 1.1e-15), pressure_comparison.txt (10 digits, all rows): 0 printed differences;
an independent activity at 4 rows by direct facet sums: <= 1.8e-15; Jacobian vs finite differences <= 6.8e-6; bounds;
partition of unity; stiffness (Lg 1 = 0, linear-function energy 4e-16). After a complete fit `localnn_direct.py`
re-evaluates the final per-facet activity by direct facet sums at every row (timing test 124995: cost rel 4.8e-15).

**Lambda rule** (proposed 2026-11-14, accepted by the user the same day; `localnn_fit.py start` writes the values): lambda_C = f data_cost(start) /
x_c0^T Lg x_c0, f = 0.1 / 1 / 10 (weak / medium / strong): templates off: start data cost 1.6089, x_c0^T Lg x_c0 =
3.0513e5, lambda_C = 5.27e-7 / 5.27e-6 / 5.27e-5; templates on: 0.11387, 3.0318e5, 3.76e-8 / 3.76e-7 / 3.76e-6.

**Stage 1 (user / lead decisions 2026-11-14, ~02:35 CST, "implement the recommendations")**:

- templates: A_s sum_k ft_ck is >= 0 on every facet of the bundles in D2's fit manifest (265 facets > 0, 12,364,025 = 0,
  0 < 0; max 7.894e23, smallest positive 1.70e-6; `runs/.../nonneg/localnn/template_signs_20261114.json`), so it is
  added as is: `TEMPLATES=on`, A = m_b (sum_j c_j phi_j + A_s sum_k ft_ck);
- lambda_C (SUPERSEDED ~02:55 CST by the solution-anchored rule below): three fits with the start penalty 0.1, 1, 10 x the
  start data cost (0.11387): lambda_C = 3.7559e-8, 3.7559e-7, 3.7559e-6; start = the NNLS projection of max(0, D2) - A_s sum ft (templates on), m = 1, kappa and grain
  from D2; max_nfev 4000; frozen `model_ext_frozen/localnn_20261114T081912Z`;
- jobs (start rule; cancelled after ~150 evaluations, superseded): f 0.1: 125014 (+125025), f 1: 125015 (+125026),
  f 10: 125016 (+125027); run roots kept as `runs/.../nonneg/localnn/stage1_f<f>_<UTC>.superseded_20261114T090228Z/` (the first three submissions 125007-125009 were
  cancelled after one minute because Slurm put all three on the shared node023; kept as `*.incomplete_20261114T083031Z`);
- after each fit (chained afterok; for the anchored fits from `model_ext_frozen/localnn_anch_20261114T090437Z`): `jobs/localnn_tools.sbatch`
  (export + the frozen maps and contribution tools -> `tools_state/`, `contributions/`, `maps/`),
  `jobs/localnn_report.sbatch` (`analysis/report.{json,txt}`, `fit/final/nonneg_bins.txt`, README.md, HOW-TO-RERUN.md),
  `jobs/localnn_gpu_eval.sbatch` (gpu002: the GPU standard evaluation -> `gpu_eval/fit`, `gpu_eval/gpu_check.json`), then
  README.md once more with the GPU check (`DOCS_ONLY=1`).

**Stage 1, solution-anchored lambda (user decision via the lead, 2026-11-14 ~02:55 CST; replaces the start rule)**:

- lambda_C = f C_ref / R_pilot, f = 0.1 / 1 / 10; C_ref = 4.268806657481027e-03 (the bounded softplus c9 data cost; the
  lead's 4.2688066575e-03); R_pilot = x_c^T Lg x_c (the penalty's normalization) of the pilot state = the templates-on
  start-rule f = 0.1 fit (lambda_C = 3.7558561837e-8) after 200 evaluations: data cost 3.799870121e-03, penalty 1.943515624e-05, 200 evaluations (max_nfev, not converged), job 125071 (`runs/.../nonneg/localnn/pilot_startrule_f0.1_n200_20261114T090209Z`);
  `localnn_anchor.py` -> `runs/shap8_012m_lategrid_dz1.5/nonneg/localnn/lambda_rule_anchored_20261114.json` (given to the
  fits as `SMOOTH_RULE_FILE`, checked and recorded in each fit's lambda_source);
- R_pilot = 517.4627377999832 (= the pilot's penalty / lambda exactly); lambda_C = 8.249495752351284e-07 (f 0.1), 8.249495752351283e-06 (f 1), 8.249495752351283e-05 (f 10); the penalty of the pilot state with these values is 0.1 / 1 / 10 x C_ref;
- the start-rule fits 125014-125016 (and 125025-125027, 125053-125064) were cancelled after ~150 evaluations; their run
  roots are kept as `stage1_f<f>_<UTC>.superseded_20261114T090228Z` (SUPERSEDED.txt); the pilot (job 125071, the same
  command as 125014 with --max-nfev 200, same snapshot and CPU type) reproduces 125014's logged objectives at
  evaluations 50 / 100 / 150 to all printed digits: 5.1129033437e-03, 4.8626753969e-03, 4.2835120963e-03;
- the anchored fits (templates on, start = the NNLS projection of D2 minus the template part, max_nfev 4000; frozen
  `model_ext_frozen/localnn_anch_20261114T090437Z`; EPYC 7763 nodes, 16 threads, 12 h; continuation afterany; the same
  post-fit chain): f 0.1: 125107 (node001; continuation 125108; tools 125109, report 125110, GPU 125111, docs 125112), f 1: 125113 (node002; 125114; 125115, 125116, 125117, 125118), f 10: 125119 (node016; 125120; 125121, 125122, 125123, 125124); submitted 2026-11-14 09:42 UTC (03:42 CST) with `tests/localnn_20261114/submit_anchored.sh <rule.json>`; run roots `runs/shap8_012m_lategrid_dz1.5/nonneg/localnn/stage1_anch_f<f>_<UTC>/`.

**Stage-1 anchored results and incidents (2026-11-14, up to 16:30 CST)**:

| f | lambda_C | data cost | penalty | stop | rms ln NG / RG (all; before the gap; late) | multipliers 40-50 .. 90-180 | c_j at 0 | active seen area | GPU check |
|---|---|---|---|---|---|---|---|---|---|
| 0.1 | 8.2495e-7 | 2.8331230e-3 | 1.494e-5 | max_nfev 4000 (2 segments) | 0.0553 / 0.0510; 0.0412 / 0.0527; 0.1229 / 0.0638 | 1, 1.001, 29.3, 675, 2269, 3.35 | 36 % | 0.136 | passed (125185: cost rel 2.2e-16) |
| 1 | 8.2495e-6 | - | - | resubmitted (125177) | - | - | - | - | chained (125181) |
| 10 | 8.2495e-5 | 2.8306496e-3 | 6.305e-5 | max_nfev 4000 (2 segments) | 0.0558 / 0.0505; 0.0426 / 0.0536; 0.1205 / 0.0501 | 1, 1.504, 209, 1608, 3617, 3.33 | 9.3 % | 0.201 | passed (125183) |

(D2 4.886e-3, bounded softplus c9 4.2688e-3; kappa and grain B at 0 in both; Q_H2O(10:33:10) 1.10e26 / 2.89e25 s^-1;
50 % of the seen production from 0.05 % / 0.09 % of the seen area; `<run root>/analysis/report.txt`.)

- f = 1: the first run (125113) died with node002 (NODE_FAIL at 08:18 CST after ~1400 evaluations; no segment, the fit
  writes its state only at a segment's end); moved aside as `stage1_anch_f1_20261114T094258Z.incomplete_*`; resubmitted
  on node016 (125177, continuation 125178 not pinned, post-fit chain 125179-125182; run root
  `stage1_anch_f1_20261114T222136Z`).
- GPU check: with grain A, B held the D2 path fails (125052, 125123: `fixed_parameters` sets the grain start at x0[-2]
  and `grain_crossing` reads the grain from the last two entries of the reduced vector); fixed in
  `model_ext_frozen/localnn_post2_20261114T222122Z`: grain A, B free (the fitter's start, A = the observed jump, B = 0),
  `--max-nfev 1`, the CPU state evaluated with the GPU's grain for the comparison. f = 10 (125183): data cost GPU
  2.841486039784887e-3 vs CPU 2.8414860397848836e-3 (rel 1.2e-15), 0 printed differences in COPS_model.dat (all rows)
  and late_window_model.dat.
- m_b vs c_j (`localnn_degeneracy.py`, `<run root>/analysis_degeneracy/`): the c-field is nearly bin-independent at
  150 m (<sum c phi> over the bins' seen facets 3.7-6.1e15 for f 10), so the multipliers alone set the steep/flat
  contrast: <A> on the seen facets relative to 0-40 deg = 0.77, 1.05, 134, 979, 2156, 703 (f 10). No node is dominated by
  a steep bin (<= 2 nodes with >= 50 % of their seen support in a bin >= 40 deg), so the trade-off is not node-local:
  it is the global level of the c-field, which the Laplacian does not penalize: lowering every c_j and raising the six
  multipliers by the same factor changes the residual vector by only 2.3 % (f 10) / 4.3 % (f 0.1) of the multiplier step
  alone (only the 0-40 deg bin, m held at 1, and the 90-180 deg template part change); a single multiplier is
  compensated by the c_j to rho = 0.2-0.3 % (full Jacobian with the smoothness rows). The absolute multipliers and the
  flat-terrain activity level are therefore weakly determined; the activity on the 60-90 deg facets (92-94 % of the late
  NG, 86-90 % of the late RG model signal, f 0.1 / 10) is what the data fix.

**Multiplier penalty P_m (user decision 2026-11-14 ~17:30 CST; for stage 2, chosen by the user later)**:

- P_m = f_m C_ref sum_b (ln m_b)^2 over the six free multipliers, f_m = 0.01, C_ref = 4.268806657481027e-03 (the bounded
  softplus c9 data cost; the lead's 4.2688066575e-03): the multipliers stay as small as the data allow, one point of
  the flat direction (A12) is picked. Residual rows sqrt(2 f_m C_ref) ln m_b (0.5 r.r scale), exact Jacobian
  sqrt(2 f_m C_ref) / m_b; reported separately (`final_mult_penalty` in fit_result.txt, `costs.mult_penalty` and
  `mult_penalty` in fit_manifest.json); objective = data + crossing + smoothness + P_m.
- `localnn_fit.py fit --mult-penalty-f F` (`jobs/localnn_fit.sbatch` MULT_PENALTY_F; default off); `localnn_fit.py start
  --start state --state-file <seg>/localnn_state.npz` starts at another fit's end point (START= of the job);
  `jacobian_check`'s node-column step is now scaled by the state's median x_c (at a converged state x_c ~ 1e-3; the old
  step 1e-4 x 1 gave FD truncation errors of 5e-4, with or without P_m). Frozen `model_ext_frozen/localnn_mp_20261115T020233Z`.
- tests (`test_localnn_multpen.py`, `tests/localnn_20261114/multpen/tests_multpen.json`): option absent and f_m = 0 =
  the previous code bitwise (fun, jac and a 3-evaluation TRF run vs `localnn_anch_20261114T090437Z`); P_m rows = f_m C_ref
  sum (ln m)^2 (rel 2e-16); Jacobian vs finite differences <= 8.5e-6 (multiplier columns 3e-9).
- 2026-11-21 (user decision 3c): `localnn_fit.py fit --slope-test` with `--hold-mult slope_40_50,...` (held at 1) or
  `--slope-mult-min L` (0 < L < 1), for control fits only; CPU; without the option the fitter is unchanged (guards as
  before). Frozen `model_ext_frozen/localnn_slopetest_20261122T033334Z`; test `test_localnn_slopetest.py` (job 137892):
  default path = the old code (x0, bounds, residuals, Jacobian, costs identical), Jacobians vs FD <= 3.1e-9.
- check run (user: "one continuation of the f = 10 fit from its end point with P_m on"): job 125360 (node001,
  max_nfev 1000), report 125361, degeneracy 125362; run root
  `runs/shap8_012m_lategrid_dz1.5/nonneg/localnn/stage1_anch_f10_pm0.01_cont_20261115T020243Z`. At the f = 10 end point
  P_m = 6.48e-3 = 2.3 x its data cost (2.8306e-3): f_m = 0.01 is not weak there. Result: ftol after 219 of 1000 evaluations (converged, first-order optimality 3.5e-9). Data cost 3.2655e-3 (+15.4 % vs the f = 10 end point 2.8306e-3; still -23.5 % vs c9), smoothness 1.69e-4, P_m 1.96e-3 (6.48e-3 at the start), objective 5.395e-3 (9.374e-3); rms ln NG / RG 0.0604 / 0.0537 (late 0.130 / 0.063). Multipliers 40-50 1 (bound), 50-60 1 (bound), 60-70 1.010, 70-80 827, 80-90 1.003, 90-180 2.42: the penalty removes the 60-70 and 80-90 deg enhancement and keeps one steep bin; mean activity on the seen facets relative to 0-40 deg 0.87, 0.81, 0.80, 654, 0.79, 190; the 70-80 deg bin carries 86 % of the seen production, 88 % of the late NG and 80 % of the late RG signal; the c-field level rose (<sum c phi> 1.3e16 on every bin, was 3.7-6.1e15) while 63 % of the nodes went to 0; active seen area 24.5 % (f = 10: 20.1 %); 50 % / 90 % of the seen production from 0.10 % / 0.47 % of the seen area; the global level direction now costs 6.2 % of the multiplier step (was 2.3 %), and the multiplier weight in the 12 weakest directions is on 70-80 deg only (0.94); direct facet-sum check passed (7.8e-16); `<run root>/analysis/report.txt`, `analysis_degeneracy/degeneracy.txt`.

**Stage 2 (user decisions 2026-11-15 ~12:30 CST; lead 12:40 CST "option 3")**: four 5-cycle fit / fine-tuning
chains, smoothness f = 1 and f = 10 (lambda_C 8.249495752351283e-06 / 8.249495752351283e-05 from
`lambda_rule_anchored_20261114.json`) x multiplier penalty f_m = 1e-4 and 0.01; templates on; fit max_nfev 4000;
fine tuning as in the softplus chains (late_only facets, WINDOW_MIN 5, no lower bound on delta, grain refitted,
lambda 0, MAXNFEV 500, LSMR 50, ftol/xtol/gtol 1e-10); the accumulated h held in the next fit.

- model: A_e = m_b (sum_j c_j phi_j(e) + A_s TS_e) exp(h_e) (`localnn_problem.Problem(..., held)`): the held facets
  (h != 0, all among the 271 late facets) enter as the exact linear correction (exp(h_e) - 1) m_b (u_e + TS_e) G[:, e]
  with G = dp/dA of the late kernels (`localnn_late.py` -> `runs/.../nonneg/localnn/late_kernels_150m.npz`, job 132008:
  the 9 row masks the frozen `positive_finetune.late_only_facets` asks for, the 271 late facets, G_NG / G_RG of all 1551
  rows); exact Jacobian; held absent / h = 0 = the previous code bitwise.
- CPU fine tuning (`localnn_ft.py run`): the frozen `nonneg_finetune.run` (`model_ext_frozen/nonneg_20261113T035423Z`,
  unchanged) with a CPU engine on the late kernels; held_next = held_prev + delta is held in the next cycle's fit
  (`localnn_fit.py start --start state --grain-from <ft>/finetune_state.npz --held <ft>/finetune_state.npz`). Needs
  `NPY_DISABLE_CPU_FEATURES` (AVX-512 off): np.exp differs in the last bit between the AVX2 and AVX-512 loops and the
  frozen context refuses an activity it cannot reproduce bitwise (131977).
- verification against the GPU fine tuning of the softplus chain `chain_5cycles_smm1_20261113T041215Z`
  (`localnn_ft.py verify-softplus`, 131978, 131979): the same 271 late facets and parameters in every case, start
  costs equal to 8e-16; c9: 4 evaluations each, final cost equal bitwise (4.268806657481027e-03), delta max diff
  5.4e-11 (max |delta| 0.011); c5: 46 / 46 evaluations (ftol), final 4.286063145045e-03 (CPU) vs 4.286063164388e-03
  (GPU, rel 4.5e-9), delta max diff 0.0063 (max |delta| 9.65), iterations 0-3 identical; the same CPU code from the
  c5 start perturbed by 1e-15 relative gives 4.286063190777e-03 (41 evaluations) and 4.286063220094e-03 (42; rel 6.2e-9
  and 1.3e-8 to the GPU, delta diffs 0.040 / 0.019): the CPU-GPU difference is within the rounding amplification of the
  truncated LSMR on the flat ftol end; c1: CPU 4.441499757894e-03 (126 evaluations, ftol) vs GPU 4.4415011703e-03
  (500, max_nfev). Chain consistency: the c2 fit of the f10 / 0.01 chain starts at data cost 2.8885631108e-03 = the c1
  fine tuning's final 2.888563110791e-03.
- checkpoints: `localnn_fit.py fit --checkpoint-every 250` writes `<fit>/checkpoint.npz` (best x, nfev_total, settings)
  atomically; the same command resumes from it (settings checked; kept as `checkpoint_resumed_segNN_<UTC>.npz`).
- tests (`test_localnn_stage2.py`, 132009 from the frozen snapshot, 5 of 5 passed): held h = 0 bitwise; held term vs
  direct facet sums at 4 rows (8.9e-16); Jacobian vs FD with held + P_m (2.8e-7; held-facet nodes 6.1e-7); checkpoint
  resume after SIGKILL at 8 of 14 evaluations (resumes at 8, ends at 14, one segment); `localnn_ft.py run
  --validate-only`: validation checks 1-8 pass. The first run (131980) failed test 3: the FD step came from the median
  of the positive x_c (8e-20: with P_m half of the nodes sit at the bound 0 as 1e-24 .. 1e-19), now the 75th
  percentile (`node_fd_scale`, 1.1e-3 .. 1.5e-3 at the stage-1 end points) and a second-order one-sided difference at
  the bounds (the multipliers at m = 1).
- driver `jobs/localnn_chain.sbatch` (layout `<chain>/c<k>/{fit,ft}`, `chain_events.jsonl` with device / host / job per
  stage, `provenance-<job>/` per attempt with `commands.sh`; a `<chain>/HANDOVER` file naming another device stops a
  job before its next stage with exit 75), submitted by `jobs/localnn_chain_submit.sh` (attempts afternotok, the post
  job `jobs/localnn_chain_post.sbatch`: degeneracy, report, maps / contributions of the final fine-tuned state
  (FT_DIR), chain docs `localnn_chain_docs.py`). Frozen `model_ext_frozen/localnn_s2_20261115T185842Z`.
- chains (2026-11-15 19:00 UTC = 13:00 CST, normal128cores, 16 threads, 12 h attempts, 8 attempts each, first
  attempts on node004-node007, later ones `--exclude=node002,node[019-024]`):
  `runs/shap8_012m_lategrid_dz1.5/nonneg/localnn/stage2_f1_pm1e-4_20261115T190017Z` (132010-132017, post 132018),
  `stage2_f1_pm0.01_20261115T190017Z` (132019-132026, post 132027), `stage2_f10_pm1e-4_20261115T190018Z`
  (132028-132035, post 132036), `stage2_f10_pm0.01_20261115T190019Z` (132037-132044, post 132045). Cycle 1: a
  continuation (P_m on) from the stage-1 end point of its f; for f10 / 0.01 the check run
  `stage1_anch_f10_pm0.01_cont_20261115T020243Z/fit` IS cycle 1's fit (linked as c1/fit). Each chain root has
  README.md (status table) and HOW-TO-RERUN.md.
- GPU solver path (`localnn_gpu_solver.py`, `localnn_fit.py fit --device gpu [--gpu-svd-mode qr|gesvd]`; default cpu
  = the previous path): profile on 16 EPYC threads (f10 / 0.01 state): fun 0.008 s, jac 0.09 s, the TRF SVD of
  J_augmented (8522 x 2754) 9.6 s of ~11 s per evaluation. The design products and the SVD go to one GPU (CuPy FP64;
  qr = QR + SVD(R), U = Q U_R); the rest of TRF is scipy's code on the CPU. Verification (`localnn_gpu_solver.py
  verify`, `jobs/localnn_gpu_solver_verify.sbatch`, tolerances: residuals 1e-13 of max |r|, Jacobian columns 1e-12,
  TRF step 1e-8 (L2, relative), 50-evaluation objective trajectory 1e-10) on the c2 start of the f10 / 0.01 chain
  (held h of its c1 fine tuning; case `tests/localnn_20261114/gpu_solver/case_c2_f10_pm0.01.json`; CPU reference 132047:
  50 evaluations, 11.11 s per evaluation, objective 5.017757127426e-03 -> 4.945852513880e-03); the default cpu path
  with the new option is bitwise the frozen stage-2 code (6-evaluation fit, x equal; 132049). Frozen
  `model_ext_frozen/localnn_s2gpu_20261115T192024Z`; verification 132810 (1 GPU, gpu001 or gpu002, never gpu003 - lead
  2026-11-15 13:30 CST, `scontrol update ExcNodeList=gpu003`; agent 3's 8-GPU A/B/D attempts run first on gpu002); `jobs/localnn_gpu_move.sbatch` 132811 (afterok) submits, if it passes and the GPU needs <= 0.5 x the CPU time
  per evaluation, GPU attempts for the four chains (`jobs/localnn_chain_to_gpu.sh`: 1 GPU, 6 x 8 h, gpu001 or
  gpu002, never gpu003) that take their chain over when they start (they cancel the chain's CPU attempts then; a fit resumes from its
  last checkpoint, <= 250 evaluations of CPU work lost). More CPU threads do not help: the SVD takes 9.6 s on 16, 12.1 s
  on 32, 10.8 s on 64 threads.
- post-fit tools of a fine-tuned state (fix 2026-11-15, frozen `model_ext_frozen/localnn_s2post_20261115T195616Z`):
  `localnn_tools.py export --ft-dir` now writes `fit_state/` (the fit's own late_window_model.dat, nonneg_facets.npz,
  ...) and the tools use it as --fit-dir; before, one folder served both and nonneg_contrib compared the fit state with
  the fine tuning's late_window_model.dat (max rel 2.6e-2; post job 132045 failed in contrib). The post script's
  step list is POST_STEPS (STEPS leaked into localnn_tools.sbatch: 132812 failed). f10 / 0.01 rerun 132818: contrib
  verification NG / RG final and fit state 4.4e-10 (tables print 10 digits), maps verified; the other chains' post
  jobs resubmitted with the fixed snapshot (132814-132816; 132018, 132027, 132036 cancelled), the GPU move 132817
  (132811 cancelled) and the GPU cross-check 132813 (8 GPUs, `#SBATCH --exclude=gpu003`) use it.
- result of the complete chain f10 / 0.01 (all fits ftol; full objective 9.374e-3 -> 4.910e-3, monotone over fits and
  fine tunings; data cost 2.889e-3 after c1, 3.0375e-3 at the end): FINDINGS.md A12.
- after-success attempts of the afternotok chain: when attempt k succeeds, k+1 is cancelled (invalid dependency) and
  k+2 (afternotok of a cancelled job) runs as a no-op (stages complete; an extra provenance folder and
  chain_complete line); harmless.
- fine-tuning regularization tests (lead / user 2026-11-15 ~21:10 CST; agent 1's options in
  `model_ext_frozen/nonneg_20261116T031115Z`: P_delta = f_delta C_ref sum w (h+delta)^2 / sum w over the tuned facets,
  w = area, C_ref 4.2688066575e-03, cap |h+delta| <= ln 100): the c5 fine tuning re-run on CPU (no chain changed) for
  C f1 / 0.01, C f10 / 0.01 (`localnn_ft.py run`), B bg1e-2 and D logspread f0.1 (`localnn_ft.py run-nonneg`: the real
  NonnegFitContext, the fit-state pressures of `localnn_late.py --state`, job 134495), f_delta 0.1, 1, 10, MAXNFEV 500;
  frozen `model_ext_frozen/localnn_s2ftreg_20261116T031038Z`; jobs 134504 (C), 134505 (B, D), summary 134507
  (`localnn_ftreg_summary.py`); folder `tests/localnn_20261114/ftreg/` (`summary_134507.txt` / `.json`). CPU controls
  without the options vs the GPU c5 fine tunings of B / D: final cost rel 6.3e-14 / 4.4e-16, same evaluations (10 / 4),
  start rel 3.3e-16 / 2.2e-16 (134501, 134505). Results: FINDINGS.md A12.
- stage-2 reruns (user decisions 2026-11-15/16 via the lead): the regularized fine tuning (f_delta = 10, C_ref
  4.2688066575e-03, cap ln 100; agent 1's `model_ext_frozen/nonneg_20261116T031115Z`) and D2's templates tested both
  ways on option C, smoothness f = 10 (lambda_C 8.249495752351283e-05), f_m = 0.01, 5 cycles; frozen
  `model_ext_frozen/localnn_s2r_20261116T044801Z` (driver TEMPLATES / C1_START_EXTRA / FT_EXTRA / NONNEG_CODE;
  `localnn_fit.py start --state-drop-templates`; `localnn_chain_summary.py`). Submitted 2026-11-16 04:48 UTC:
  `runs/.../nonneg/localnn/stage2r_f10_pm0.01_tplon_20261116T044825Z` (templates ON; cycle-1 fit = the check run
  `stage1_anch_f10_pm0.01_cont_20261115T020243Z/fit` as in `stage2_f10_pm0.01_*`; jobs 135070-135072, post 135073,
  node007) and `stage2r_f10_pm0.01_tploff_20261116T044825Z` (templates OFF; cycle-1 fit from the stage-1 f = 10 end
  point with the template part removed: the same c, m, kappa, grain, start data cost 5.4582e-03 (2.8306e-03 with the
  templates), late rms ln NG / RG 0.159 / 0.230; jobs 135074-135076, post 135077, node008). Per-cycle numbers:
  `tests/localnn_20261114/stage2r/chain_summary_<job>.txt` (the unregularized f10 / 0.01 chain: 135078).
  Results: ON (23 min, all ftol) data 3.1996e-3 after c5, objective 5.2623e-3; OFF (58 min; c1 fit 238 evaluations,
  then the fine tunings and fits do not move: the late facets carry 0.01 % of the late signal without the templates)
  data 3.6519e-3, objective 6.0628e-3; summaries 135107 (ON), 135114 (OFF), shares `share_c5_*.txt`; FINDINGS.md A12.
- the remaining combinations with templates ON and the regularized fine tuning (user decisions 2026-11-16 ~00:10
  CST; the same settings as `stage2r_f10_pm0.01_tplon`), submitted 2026-11-16 05:59 UTC, 8 x 12 h attempts each:
  `stage2r_f1_pm1e-4_tplon_20261116T055923Z` (135332-135339, post 135340; cycle-1 fit = the c1 fit of the original
  `stage2_f1_pm1e-4_20261115T190017Z`, linked once complete: gate job 135331 waits for its final + direct check),
  `stage2r_f1_pm0.01_tplon_20261116T055923Z` (135341-135348, post 135349, node017; c1 fit of
  `stage2_f1_pm0.01_20261115T190017Z`), `stage2r_f10_pm1e-4_tplon_20261116T055923Z` (135350-135357, post 135358, node018;
  c1 fit of `stage2_f10_pm1e-4_20261115T190018Z`): the c1 fits are the original stage-2 continuations (stage-1 end point,
  P_m on), so only the fine tunings and cycles 2-5 differ. The unregularized chains 132010 / 132028 stay on CPU as
  references (GPU move 132817 cancelled).
- 2026-11-16 (lead 10:05 CST): GPU solver verification 132810 (gpu002, A100 80GB PCIe) PASSED: residuals 1.9e-15 of
  max |r| (tolerance 1e-13), Jacobian columns 1.0e-15 (1e-12), singular values 1.5e-15, one TRF step 1.5e-13 (qr) /
  2.1e-13 (gesvd) (1e-8), 50-evaluation objective trajectory 1.1e-11 (qr) / 2.5e-11 (gesvd) (1e-10), the same 50
  evaluations; 1.33 s per evaluation on the GPU vs 11.1 s on 16 CPU threads (SVD 1.28 s vs 7.9 s). GPU cross-check
  132813 of the f10 / 0.01 final state (c5 fine tuning; D2's path, 8 GPUs): data cost GPU 3.052435634290019e-3 vs CPU
  with the GPU's grain rel 5.6e-16, tables 0 printed differences. The chains stay on CPU (agent 3's GPU reruns have
  gpu002; GPU move 132817 cancelled).
- 2026-11-16: the unregularized f1 / 1e-4 chain failed its c4 Jacobian check (2.89e-5, jobs 132011-132017): FD rounding
  noise at bound nodes (h = 3.8e-9; the error fell as 1/h), not the Jacobian. Fix of the check only (noise floor of
  the step), frozen `model_ext_frozen/localnn_s2jc_20261116T160722Z`; the chain resumed at c4 (135854-135860; check
  7.4e-7); the pending attempts of the three running chains replaced by attempts with that snapshot
  (`jobs/localnn_replace_pending.sh`: stage2r_f10_pm1e-4_tplon 135875-135882, stage2r_f1_pm1e-4_tplon 135883-135890,
  stage2_f10_pm1e-4 135891-135897; the running attempts untouched).
- 2026-11-17: the /data001/highresolu_vt file system filled (0 B free at ~03:10Z); `stage2r_f10_pm1e-4_tplon` (135875) and
  `stage2r_f1_pm1e-4_tplon` (135883) died in their c3 fits, their later attempts and posts failed (exit 0:53). The c3
  checkpoints were written before (03:06:12Z, 3750 evaluations; 02:51:41Z, 3500), zip test and arrays OK; no partial
  products (c3/fit held only checkpoint.npz and start). Resubmitted 06:27Z with the same settings (MODEL_EXT
  localnn_s2jc_20261116T160722Z), the c3 fits resumed from the checkpoints: 136006-136009 + post 136010 (node018),
  136011-136014 + post 136015 (node008). Job scripts: disk guard (localnn_chain.sbatch: >= 20 GB free before every
  stage, exit 4; localnn_chain_post.sbatch: per step, TOOLS_STEPS); the two posts run without the maps
  (TOOLS_STEPS "export contrib": the maps are 15-20 GB per chain; disk budget ~30 GB).
  Both resubmissions then failed the c4 Jacobian check (07:20Z, 08:03Z; 2.15e-5 / 2.49e-5 in c_slope_50_60 /
  c_slope_60_70 at m = 1.0000570 / 1.0000503): for 1 + h/2 <= m < 1 + h the check used the asymmetric two-point
  difference (first-order). Fixed (check only): the second-order one-sided difference whenever the bound clips the lower
  point, the noise-floor cap 1000 x; frozen `model_ext_frozen/localnn_s2jc2_20261117T080730Z` (5 states <= 8.2e-7);
  resubmitted 08:08Z: f10 / 1e-4 136032-136035 + post 136036 (node018), f1 / 1e-4 136037-136040 + post 136041 (node008);
  c4 checks passed (5.9e-7, 8.2e-7); stage-2 tests of the snapshot 136043.
  Device switch (lead 2026-11-17 11:45 CST, gpu002 free): the remaining stages of both chains moved to gpu002 (1 GPU
  each; `jobs/localnn_chain_to_gpu.sh --after-checkpoint`: the GPU attempt cancels the CPU attempt seconds after the c4
  fit's next checkpoint write, the fit resumes from it on the GPU; same settings, check and disk guard): f10 / 1e-4 GPU
  attempts 136092-136093 + post 136094, f1 / 1e-4 136095-136096 + post 136097 (posts without maps). Stages: c1 (linked
  original fits) - c3 and the c4 fit up to the takeover checkpoint on CPU, the rest of c4 and c5 on the GPU
  (chain_events.jsonl device per stage, the fits' meta "device").
  f10 / 1e-4 finished on the GPU (136092; c5 fit data 2.458706e-3, ft 2.45855e-3). Its post 136094 failed in contrib:
  "A_reproduced_bitwise: false" - not the GPU: the post ran on node028 (EPYC 9554, AVX-512) and the tools job script
  did not set NPY_DISABLE_CPU_FEATURES, so np.exp(h) differed in the last bit at 38 of the 271 held facets (rel <=
  3.3e-16); with the setting, and on node018 (EPYC 7763), bitwise (`jobs/localnn_postA_repro.py`, sruns 136119-136122;
  the GPU fit writes its activity on gpu002's CPU, EPYC 7513, no AVX-512). The export's 4.78e-7 is the 7-digit print
  precision of COPS_model.dat (0 printed differences). Fix in the job scripts only (NPY_DISABLE_CPU_FEATURES in
  localnn_tools / chain_post / report / degeneracy .sbatch, as the driver); the check stays bitwise. Rerun 136123 on
  node028: contrib verification 4.7e-10 (final and fit state). kappa ~0 in every C stage-2 state: FINDINGS.md A12.
  f1 / 1e-4 finished on the GPU (136095, 20:37Z; post 136097 on node028 with the fix, contrib 4.2e-10). Stages:
  c1 (linked) - c3 and c4 up to its checkpoint at 3250 evaluations on CPU, the rest of c4 and c5 on gpu002
  (chain_events.jsonl). Per-cycle numbers of both 1e-4 chains and their unregularized references:
  `tests/localnn_20261114/stage2r/chain_summary_136129.txt`; FINDINGS.md A12.
- 2026-11-18 (H2O night law q_min = 0.02 q_max, user decision; FINDINGS.md A18): every option C product above holds
  the 0.2 law. (1) Regression of test_localnn check 5 (Jacobian vs FD at the D2-projection start, 6.0e-3 / 9.4e-3): the
  check, not the Jacobian (bisect `tests/localnn_20261118_jacreg/bisect.out`: the 2026-11-14 code passes, the
  state-scaled node steps from `localnn_mp_20261115T020233Z` on fail; the 2026-11-14 check on the current Problem passes,
  6.8e-6): jacobian_check now tries several steps per column (1e-4 max(|x|, 1), 1e-4 max(|x|, scale), each also
  noise-floored), the smallest discrepancy counts; 2.3e-8 / 2.5e-8 at the check-5 states, <= 5.9e-7 at five converged
  stage-2 states (`model_ext_frozen/localnn_s2jc3_20261118T081235Z`). The stage-2 C fits are not affected (the check
  only reads the model). (2) Design matrices and late kernels of the new law, new folders (the 0.2 ones kept):
  `grid_150m/design_h2oqmin0.02` (137 MB) and `late_kernels_150m_h2oqmin0.02.npz` (158 MB), job 136507 (contribution
  code `nonneg_contrib_20261118T072729Z`, late facets of `nonneg_20261118T072729Z`, the same grid, bins, seen mask and
  D2 templates TS as the 0.2 build; 2 + 1 min). Checks: design vs direct facet sums with the new law 3a / 3b <= 9.4e-15 /
  2.4e-15; held term vs direct sums 2.4e-15; Jacobian 2.4e-8; the same 271 late facets. Expected failures (the law):
  test_localnn 8 (the 0.2-law D2 / c9 tables: data cost 0.0798 / 0.0648 vs 4.89e-3 / 4.27e-3) and test_localnn_stage2
  5 (start cost of a 0.2-law fit). New / 0.2 design columns: G 0.1 - 1 (median 0.79 NG, 0.80 RG), templates median
  0.71, CO2 1.0; the late kernels median 0.1 (most late-facet records are night / shadow). Law check
  (`model_ext_frozen/localnn_s2law_20261118T081952Z`): localnn_problem.load_design and localnn_late.load compare the
  product's H2O q_min ratio (meta law.h2o_qmin_ratio; older products: the fit_cops law line / contribution code of
  their build) with tools/fit_cops.py and refuse a mismatch (LOCALNN_ALLOW_LAW_MISMATCH=1 reads it on purpose, with a
  warning and a record); job-script defaults (DESIGN_DIR, CONTRIB, LATE, NONNEG_CODE) now point to the new-law products.
  Open: the data selection (clip_refit/D2/final) and the D2 templates TS are the 0.2-law D2's until the new D2 exists.
  Preflight against the new D2 (`clip_refit/D2_q002_warm/final`, 2026-11-18 10:36Z, jobs 136557 / 136558; plan
  `cops_output/model_ext/logs/plan_q002_stage2r_C.md`): data selection, references, ft templates, enhancement scale
  and TS equal to the 0.2-law D2's (TS diff 8192 on 7.9e23); base_q002 kernels (1551), elements and observation columns
  equal to base: no rebuild of the design matrices or the late kernels. The design built with the new D2 (test folder
  `tests/localnn_q002/d2check_136558/design_d2q002`) has bitwise the matrices of `design_h2oqmin0.02` and reproduces
  the new D2: data cost 5.6707116011472e-3 (rel 6.7e-16), COPS_model.dat 0 printed differences, late window <= 3.5e-10
  (`localnn_fit verify` compares with the design's STORED direct sums of its build-time D2, so with design_h2oqmin0.02
  (old D2's sums) that check is not meaningful for the new D2). Start files of the new D2: `grid_150m/start_q002_templates_on`
  (data cost 0.168; 0.114 with the 0.2-law D2), `start_q002_templates_off` (1.551; 1.609).
  Lead decision 2026-11-18: that design is the design of the q002 chains, copied to
  `grid_150m/design_h2oqmin0.02_d2q002` (provenance in `design_h2oqmin0.02_d2q002.provenance/`, sha256 in SHA256SUMS)
  and passed with `jobs/localnn_chain_submit.sh --design-dir` (DESIGN_DIR of the attempts and the post job;
  `jobs/localnn_chain_to_gpu.sh` carries it).
  SUBMITTED 2026-11-18 16:24Z (user OK, posts with maps): the five chains
  `nonneg/localnn/stage2r_q002_{f1_pm1e-4_tplon,f10_pm1e-4_tplon,f1_pm0.01_tplon,f10_pm0.01_tplon,f10_pm0.01_tploff}_20261118T162407Z`
  (code `model_ext_frozen/localnn_q002_20261118T162327Z`; CPU attempts 136589-136610, gpu001 takeover of the P_m 1e-4 chains
  136612-136615, posts 136601/136606/136611/136616/136619; job table in `cops_output/model_ext/logs/plan_q002_stage2r_C.md`).
- NFS note (Voyager): `/data001` is mounted with `acregmin=acregmax=600`: for up to 10 min after an edit a compute
  node can see the old size / mtime of a .py file and import its stale `__pycache__` .pyc (seen 2026-11-15: a test ran
  the previous `jacobian_check`). Frozen snapshots are new files and not affected; tests of freshly edited run copies
  set `PYTHONPYCACHEPREFIX` to a per-job folder.

**Post-fit tools**:

- `localnn_gpu_eval.py` + `jobs/localnn_gpu_eval.sbatch`: the final per-facet activity through D2's own GPU path
  (`clip_refit.py --variant D2`, frozen `clip_20261108T234013Z`, unchanged grain_fit_fixed / gpu_fit / fit_cops):
  template localnn_a = A / A_s held at 1, a00..a24 = 0, kappa, grain A, B, T held (exact repr), one free amplitude
  that cannot change any pressure (template probe_unseen = 1 at one facet no fit point sees; fixed_parameters refuses
  to hold every parameter and fit_cops refuses an all-zero template), `--max-nfev 1` (TRF returns its start); the
  fitter's standard products and a comparison with the CPU products;
- `d2_tools.py` + `jobs/d2_tools.sbatch` (2026-11-17): the clipped refit D2 for the maps / contribution tools and the figure scripts: `export` (A = max(0, activity_linear), re-derived bitwise from D2's restart; seen from D2's activity_facets.dat; slope bins of the nonneg fits; the observation-row references from pressure_residuals.dat model_Pa, full precision), `maps` / `contrib` (the frozen tools with state_activity replaced), `verify` (forward model of the export vs D2's per-row model and data cost); `make_state_activity(NM)` for `contrib_phase_data.py --state-adapter`. `nonneg_maps.py` / `nonneg_contrib.py` (2026-11-17): state_activity also reproduces the alternative models sos / background / logspread bitwise (frozen nonneg_maps_20261117T224340Z / nonneg_contrib_20261117T224340Z); softplus / squared unchanged (tested on the c9 state).
- `localnn_tools.py` + `jobs/localnn_tools.sbatch`: the frozen `nonneg_maps.py` / `nonneg_contrib.py` (folder
  `nonneg_contrib_20261112T210012Z`) run UNCHANGED on a localnn state: an export folder in the chain layout (marked: no
  fine tuning; `fit_restart.npz` there is not a gpu_fit restart) and `state_activity` replaced by a localnn version (A
  re-derived bitwise from c, m, the grid and TS). Tested on the 40-evaluation timing state: contributions (the tool's
  own verification vs the tables: every row within 4.8e-10) and maps (VTU bitwise, Tecplot re-parse within 5e-10;
  ~15 GB per state);
- `localnn_report.py` + `jobs/localnn_report.sbatch`: the numbers of the lead's item 5 (with D2 and c9 for comparison)
  and the run-root documents.

**Cost**: CPU only; the TRF step's dense SVD (8516 x 2754) is 5.7 s at 16 OpenBLAS threads (node032; 32/64/128
threads slower); the timing test ran 10.0 s per evaluation on a shared node023 with 32 threads.

**Code and jobs**: `localnn_mesh.py`, `localnn_design.py`, `localnn_problem.py`, `localnn_fit.py` (verify | start | fit),
`localnn_direct.py`, `test_localnn.py`; `jobs/localnn_build.sbatch` (CPU: mesh, design, verify, start, tests),
`jobs/localnn_fit.sbatch` (CPU: one fit, wall-time stop / exit 75, `RUNROOT=` to continue, then the direct check);
snapshot `model_ext_frozen/localnn_20261114T081912Z` (with an unchanged copy of `positive_problem.py`); the contribution
arithmetic comes from `model_ext_frozen/nonneg_contrib_20261112T210012Z`.

**Running it** (`$J=$COPS_OUTPUT_ROOT/jobs`, `$LX` = the frozen snapshot):

```bash
sbatch --export=ALL,MODEL_EXT=$LX $J/localnn_build.sbatch          # -> runs/.../nonneg/localnn/grid_150m (job 125004)
sbatch --export=ALL,MODEL_EXT=$LX,TEMPLATES=off,SMOOTH_FACTOR=1,NFEV=<n> $J/localnn_fit.sbatch
                                                                    # -> runs/.../nonneg/localnn/localnn_toff_f1_<UTC>/fit
sbatch --export=ALL,MODEL_EXT=$LX,TEMPLATES=off,SMOOTH_FACTOR=1,NFEV=<n>,RUNROOT=<run root> $J/localnn_fit.sbatch
                                                                    # continue after exit 75 (wall-time stop)
# 2026-11-14, solution-anchored lambda: the rule file from localnn_anchor.py, one lambda of it per fit (exact repr),
# fit + continuation + post-fit chain for f = 0.1, 1, 10 (as run):
$COPS_OUTPUT_ROOT/tests/localnn_20261114/submit_anchored.sh $R/nonneg/localnn/lambda_rule_anchored_20261114.json
```

Products of a fit (`<run root>/fit/segNN`, `final` -> the segment that ended by an optimizer condition):
`fit_result.txt`, `localnn_state.npz`, `nonneg_facets.npz` (u, TS, A, P, bin, seen, held = 0, area), `COPS_model.dat`,
`late_window_model.dat` (D2's formats; CPU design route, checked by `direct_check_segNN/`), `fit_manifest.json`,
`README.txt`, `HOW-TO-RERUN.md`; `<run root>/provenance-<job>/`. The maps / contribution tools (`nonneg_maps.py`,
`nonneg_contrib.py`) check the softplus / squared reproduction and need an adapter for this model.

## Alternative non-negative activity models (options A-D, 2026-11-14)

**Why** (user, 2026-11-14): "prepare suggestions for modification of the minimization polynomial: the current approach
shuts down almost all the surface, and I need a polynomial that is non-negative and would allow fit with more seen
surface contributing". "All those options need to be switchable in CLI". In the bounded softplus fit (s = 1e-5 A_orig)
A is exactly 0 on about 96 % of the seen area. Then "implement the recommendations" (2026-11-14 ~02:35 CST).

The options are in `nonneg_alt.py`, with `nonneg_mult.py` and `nonneg_problem.py`:
- A: `--model sos`;
- B: `--model background`;
- D: `--model logspread`.

Option C (`--model localnn`, local non-negative hat functions) is a separate module of another agent; it follows the
same template and reporting rules.

**Common scheme (user decisions 2026-11-14):**
- lategrid base, D2's setup: W10x5 from 10:33:10 UTC, grain at the time stamp, NG = RG crossing weight 10, original
  temperature law held (Tmax 182.09922, T75 136.11149, Tmin 133 K), κ free, λ = 0, `--tbounds 30 400`;
- D2's six ft templates ft_c1..ft_c6 held at amplitude 1, no `cliff_ge_60deg`;
- slope-bin multipliers m_b: [0,40) deg held at 1, the six others free with lower bound 1 and start 1
  (`--slope-mult-min 1`);
- max_nfev 4000; κ and grain A, B start from D2, by name.

**Templates (user decision 1):** B puts them inside the softplus argument, exactly as the softplus model does. For A
and D the rule is: if A_s Σ T_k ≥ 0 on every facet, add it as is; otherwise add max(0, A_s Σ T_k). Check of
2026-11-14 (`checks/ft_template_sign_20261114/`):
- each template is > 0 on the same 265 facets and 0 elsewhere (negative on 0 facets);
- A_s Σ T_k: minimum 0, maximum 7.894e23 (3.55e5 A_s); 350.9 m² in total;
- all 265 facets are seen; latitudes 0.5° to 69.0°, longitudes −11.6° to 166.7°.

So the sum is added as is: outside the non-negative map, times m_b exp(h) like the rest.

**Models.** H2O activity per facet e, with y_e = the 25 real harmonics (degree ≤ 4), A_s = D2's enhancement scale,
T = Σ_k T_k (the ft templates), h_e = the accumulated fine-tuning term (0 in cycle 1), A_orig = the srcCG map:

| option | `--model` | activity A_e | extra parameters |
|---|---|---|---|
| A, sum of squares | `sos --sos-rank r` | m_b (‖Lᵀ y_e‖² / A_s + A_s T_e) exp(h_e); L is 25 × r lower-trapezoidal | column 0 of L = a00..a24; the other columns `sosL_<i>_<j>` (i ≥ j ≥ 1). At r = 3: 25 + 24 + 23 = 72 parameters, 81 in the fit |
| B, background | `background --bg-min V` | m_b (A_bg A_orig,e + s_e ln(1 + exp(L_e / s_e))) exp(h_e); L = Y a + A_s T, s = 1e-5 A_orig | `A_bg` (dimensionless, relative to A_orig), lower bound V, start V |
| D, log-normal spread | `logspread --spread-f f --spread-cref C` | m_b (A_orig,e exp(u_e) + A_s T_e) exp(h_e); u = Y a / A_s | none; the spread penalty below |

- **A:** A = 0 only where all r harmonic combinations vanish together. r = 1 is the squared model.
- **B:** every facet with A_orig > 0 keeps at least A_bg·A_orig (× m_b exp(h)). The softplus part is the current
  softplus model.
- **D:** A > 0 wherever A_orig > 0.
- **The other models:** `softplus` and `squared` are unchanged (bitwise against the 2026-11-13 snapshot) and do not
  use `nonneg_alt.py`.

**D penalty (user decision 2):**
- definition: P = f · C_ref · Var_w(u), where Var_w(u) = Σ_seen w_e (u_e − ū)² with w_e = area_e / seen area, over the
  facets seen by the fit;
- values: C_ref = 4.2688066575e-03 (the data cost of the bounded-softplus c9 fit), f = 0.1, 1, 10;
- meaning: an rms spread of 1 neper costs f × the current misfit;
- implementation: P is a term of the fitter's cost 0.5 r·r, on the same scale as the data cost. It enters as 25 rows
  √(2 f C_ref) S b, with SᵀS = the seen-area-weighted covariance of y and b = a / A_s, so that
  0.5 ‖rows‖² = P exactly;
- reporting: P is not part of the data cost. `fit_manifest.json` (`model_extension.nonneg.model_check.alternative.spread`:
  f, C_ref, the value at the start and at the final state) and `NONNEG_MODEL.txt` report the data cost and P apart.

**Starts (cycle 1; user decisions 3 and 4):**
- **A, `--start gram --start-file <npz>`:** the rank-r least-squares fit of D2's activity (user decision 4: "A start:
  least-squares fit of D2's activity"). It is built in two steps:
  - `nonneg_start.py --model sos --rank r` (job 124994): the area-weighted least-squares Gram matrix Q of D2's
    activity (minimum-norm, since the products y_i y_j span only degree ≤ 8), projected onto the PSD matrices; a
    lower-trapezoidal L from its top-r eigenpairs. Relative rms 0.95, total activity × 11.0.
  - `--refine lm` (job 125003, `sos_start_r3_lm.npz`): Levenberg-Marquardt on the rank-r problem itself, from that
    start. Relative rms 0.129, total activity × 1.66; stopped at 300 iterations.
  - The lead chose the LM start on 2026-11-14 (the same decision, implemented better). Data costs by direct facet sums
    on CPU (`checks/sos_start_costs_20261114/`, job 125100): Gram 33.70346689911884 (equal to the GPU start check of
    job 124998 in every printed digit), LM 9.519127043550094; D2 is 4.886e-03.
  - The first A root (Gram start) was moved aside as `sos/stage1_r3_20261114T084711Z.superseded_20261114T091321Z`
    before it ran.
- **B, `--start d2`:** D2's harmonics a00..a24 (as the softplus), with A_bg at its lower bound V.
- **D, `--start uniform`:** u = ln s on every facet (a00 only; every other harmonic 0).
  - s is the unweighted least-squares fit in ln p of the model pressures to D2's modelled pressures, over the fit's
    used NG and RG rows. D2's predictions are matched by UTC as in the start check.
  - The pressures are linear in s, so two evaluations give p(s) = P0 + s P1; a third checks it.
  - The manifest records s, the fit-weighted alternative value and the rms d ln p
    (`model_check.alternative.uniform_start`).
- **Every model, `--start chain`:** the previous cycle, as for softplus.

**Fine tuning:** the same as softplus/squared. It applies A_e = m q_e exp(h_e + δ_e) on the late_only facets with A > 0
(`nonneg_finetune.py`, validation checks 1-8). The extra parameters are read from the fit restart by name.

**CLI** (`nonneg_fit.py`, before `--`):
```
--model sos --sos-rank 3 --start gram --start-file $R/nonneg/bundle/sos_start_r3.npz --slope-mult-min 1
--model background --bg-min 1e-3 --start d2 --slope-mult-min 1
--model logspread --spread-f 1 --spread-cref 4.2688066575e-03 --start uniform --slope-mult-min 1
```
- Every other option is as for softplus (`--ft-templates-from`, `--start-restart`, `--held`, `--grain-from`, ...).
- The model options are refused with the wrong model.
- A resume or chain refuses a changed value (the restart marker records `slope_mult_min`, `bg_min`, `spread_f`,
  `spread_cref`).
- When an option is absent, the records of the other models are unchanged.

**Job variables:**
- `jobs/nonneg_fit.sbatch`:
  - `MODEL=sos|background|logspread`;
  - `SOS_RANK`, `BG_MIN`, `SPREAD_F` + `SPREAD_CREF` (each required for its model);
  - `SLOPE_MULT_MIN=1`;
  - `START` (defaults: gram, d2 and uniform, or chain);
  - `START_FILE` (sos; default `$R/nonneg/bundle/sos_start_r<r>.npz`);
  - `HELD`, `GRAIN_FROM` (chain);
  - `NFEV`, `CHECK`, `L_THREADS`, `FITDIR`, `MODEL_EXT`.
- `jobs/nonneg_start.sbatch`: `MODEL=sos RANK=r [REFINE=lm]`.
- `jobs/nonneg_alt_test.sbatch`: a GPU test per model. It checks the start cost against an independent evaluation,
  then runs 5 evaluations with the full Jacobian check.
- `jobs/nonneg_alt_stage1.sh --model-ext <snapshot> --utc <UTC> --runs bg:1e-3,bg:1e-2,bg:1e-1,sos:3,ls:0.1,ls:1,ls:10
  --spread-cref 4.2688066575e-03 [--submit]`. It submits the stage-1 fits, one 8-GPU job each, in that order on any GPU
  node but gpu003:
  - 3 attempts per fit (attempt a > 1 afternotok of all earlier ones);
  - a guard (`positive_chain_guard.sbatch KIND=fit`) and the summary job;
  - run roots `runs/shap8_012m_lategrid_dz1.5/nonneg/<model>/stage1_<setting>_<UTC>/` with `fit/` (`fit/final`),
    `README.md`, `HOW-TO-RERUN.md` and `summary.txt`.
- `jobs/nonneg_alt_summary.sbatch FIT=<fit/final> OUT=<summary.txt>` (`nonneg_alt_summary.py`) computes:
  - the data cost and the penalties, the termination, nfev;
  - rms ln NG/RG overall / before the gap / late window;
  - the parameters at bounds;
  - Q_H2O at 10:33:10;
  - the active seen area (A > 1 % of the srcCG median A_orig of the seen facets);
  - the seen area holding 50 % / 90 % of the seen production;
  - the NG/RG shares from the active area (the fitted forward model term by term, FP32 records as in the fit).

**GPU tests (gpu002):**

| test | start cost vs independent | Jacobian check (worst) | time per evaluation | A > 0 |
|---|---|---|---|---|
| sos r 3 (124998) | 3.370346689911884e+01, rel 0.0 | 18 columns, 4.7e-7 | fun 0.38 s, jac 2.03 s | all seen facets |
| background 1e-2 (124999) | 1.409815714818658e-01, rel 0.0 | 14 columns, 1.9e-5 | fun 0.29 s, jac 1.17 s | 0.9935 of the seen facets |
| logspread f 1, start uniform (125028) | 2.557169268575762e-01, rel 0.0 | 13 columns, 7.7e-8 | fun 0.32 s, jac 0.63 s | 0.9916 of the seen facets |

- The background worst case (1.9e-5, in the harmonics) is the finite-difference crossing of the softplus kink, as in
  the softplus tests (3.5e-5).
- The D uniform start: s = 1.5503152609e-01 (ln s = −1.8641; the fit-weighted alternative would give 1.586e-01).
  - rms d ln p against D2: NG 0.343, RG 0.616; before 10:33:10 0.464, late window 0.852.
  - Linearity check: 4.4e-15.
- After 5 evaluations: data cost 0.070869, penalty 0.001526 = f C_ref Var_w(u) with Var_w(u) = 0.3575 (rms spread 0.598
  neper). fit_cops' `cost` = data + spread penalty (0.072395); the crossing penalty is reported apart.
- Per-fit time estimates (gpu002, 4000 evaluations, about 0.91 Jacobians per evaluation as in the softplus fit, plus
  about 10 min of binding and checks):
  - B: about 1.8 h;
  - A: about 2.8 h;
  - D: about 1.2 h.

**CPU tests:** `jobs/nonneg_cpu.sbatch`, 41 of 41 passed. They cover:
- softplus/squared bitwise against the previous snapshot;
- each new model's full Jacobian against FD;
- the independent start activity and A ≥ 0;
- the D penalty against an independent f C_ref Var_w(u), and the uniform start;
- fine-tuning checks 1-8, and cycle 2 starting at the fine-tuning cost;
- the sos LM refinement.

**Code:**
- modules: `nonneg_alt.py` (new), `nonneg_alt_summary.py` (new), and `nonneg_mult.py`, `nonneg_problem.py`,
  `nonneg_finetune.py`, `nonneg_start.py` and `test_nonneg.py` (extended);
- snapshots:
  - `model_ext_frozen/nonneg_20261114T084032Z`: the stage-1 fits B and A;
  - `model_ext_frozen/nonneg_20261114T085147Z`: the stage-1 fits D (= 084032Z plus a text fix in NONNEG_MODEL.txt;
    CPU 41 of 41, job 125065);
  - `model_ext_frozen/nonneg_summary_20261114T090853Z`: the summaries (`nonneg_summary_20261114T084417Z` was its test).

**Stage-1 runs (submitted 2026-11-14, `jobs/nonneg_alt_stage1.sh`; `runs/shap8_012m_lategrid_dz1.5/nonneg/`):**

| run root | fit attempts, guard, summary |
|---|---|
| `background/stage1_bg1e-3_20261114T084711Z` | 125032-125034, 125035, 125093 |
| `background/stage1_bg1e-2_20261114T084711Z` | 125037-125039, 125040, 125094 |
| `background/stage1_bg1e-1_20261114T084711Z` | 125042-125044, 125045, 125095 |
| `sos/stage1_r3_20261114T091354Z` (LM start; replaces the Gram-start root `stage1_r3_20261114T084711Z`, moved aside before it ran) | 125101-125103, 125104, 125105 |
| `logspread/stage1_f0.1_20261114T090251Z` | 125072-125074, 125075, 125097 |
| `logspread/stage1_f1_20261114T090251Z` | 125077-125079, 125080, 125098 |
| `logspread/stage1_f10_20261114T090251Z` | 125082-125084, 125085, 125099 |

Each root has README.md, HOW-TO-RERUN.md (the exact commands) and, when the fit is complete, summary.txt.

**Stage-1 convergence diagnosis (2026-11-14 ~16:20-17:20 CST, `checks/stage1_diag_20261114/stage1_diag.py`).** The
check runs on CPU from each fit's products. It reproduces the fit exactly: activity bitwise, data cost to 1e-15,
Jacobian singular values to 1e-9. It then computes the gradient scaled by the Jacobian column norms, the
bound-constrained Gauss-Newton (GN) step, and the actual objective along that step and along steepest descent:
- D f = 1 and f = 10 are genuine minima of data + penalty: scaled gradient about 2e-7, GN-predicted remaining decrease
  1e-10 relative. Their early ftol stops are real.
- D f = 0.1 is a scaling artefact on top of a model degeneracy:
  - the multiplier Jacobian columns shrank 1.2e3 to 1.8e3 times since the start, and scipy's x_scale='jac' keeps the
    largest column norm seen, so the trust region along the multipliers froze;
  - Var_w(u) does not see the uniform level a00, which trades against the free multipliers (6180 and 20637).
- B (all three) and A did not converge at max_nfev:
  - descent is still available only in small steps (B 1e-3 to 1e-2 of the GN step; A has scaled gradient 2.3e-2 and
    lost 3.2 % of its cost over the last 100 evaluations);
  - in B the softplus is active on only 2.6 to 5 % of the seen area.

**Continuation (`--start restart`, user 2026-11-14 ~17:30 CST).** A new fit starts from a COMPLETE fit of the same
model (sos, background, logspread):
- every parameter is taken by name, in exact optimizer coordinates, from its `fit_restart.npz` (grain included), with
  no held term;
- a new scipy least_squares call rebuilds the x_scale='jac' scale from the Jacobian at that point;
- the model settings must equal the source's (marker, s_ratio, bins, templates, sos_rank, bg_min, spread_f,
  spread_cref, multiplier lower bound, parameter names); otherwise the fit is refused. The spread centre may differ;
- the start check compares the start with the source's final cost (equal data cost).

I did not implement `--scale-reset N` (restarting least_squares every N evaluations inside one job). One restart per
continuation rebuilds the scale. A further continuation does it again if the multipliers keep growing.

CLI: `--start restart --start-restart <fit>/final/fit_restart.npz`. Job variable: `START=restart` with `START_RESTART`
(required).

**D penalty about a fixed level (`--spread-u0 U`, user 2026-11-14):**
- definition: P = f · C_ref · Σ_seen w_e (u_e − u0)² = f · C_ref · [Var_w(u) + (ū − u0)²];
- implementation: one more residual row √(2 f C_ref)·(ū − u0), with ū = Σ_seen w y_e · b (26 rows in total);
- default, option absent: the mean-centred penalty, unchanged (bitwise against `nonneg_20261114T085147Z`);
- value used: u0 = −1.864126788612443, the stage-1 uniform start level ln s (`model_check.alternative.uniform_start`);
- records: `model_check.alternative.spread` holds u0, the centre, Var_w(u), ū, (ū − u0)² and both penalty parts, at the
  start and at the final state. NONNEG_MODEL.txt and summary.txt report them apart from the data cost;
- job variable: `SPREAD_U0`.

**Launcher:** `jobs/nonneg_alt_stage1c.sh --model-ext <snapshot> --utc <UTC> --from ROOT[,...] [--spread-u0 U]
[--submit]`.
- It reads the model settings from each source fit and adds 4000 evaluations.
- Queue mechanics are as in stage 1: 3 attempts, guard, summary.
- Run roots are `nonneg/<model>/stage1c_<setting>_<UTC>/`; their README.md and HOW-TO-RERUN.md say "continuation of
  <source root>".

**Tests:**
- CPU (`jobs/nonneg_cpu.sbatch`, 48 of 48 passed):
  - the restart start of each model starts at the source's final data cost and x (bitwise);
  - the u0 penalty at the start and at the end matches an independent f C_ref Σ w (u − u0)²;
  - refusals: another model's restart, `--spread-u0` with another model, restart without a source;
  - the new options absent for sos, background and logspread are bitwise equal to the stage-1 snapshot.
- GPU: `jobs/nonneg_restart_test.sbatch` (job 125365; results below).

**Snapshots:** `model_ext_frozen/nonneg_20261115T020314Z` (code) and `model_ext_frozen/nonneg_summary_20261115T020314Z`
(summary).

**GPU test (job 125365, gpu002):** all three PASS. In each case the start data cost equals the source's final data
cost (rel ≤ 2.2e-16) and x_start equals the source x:

| test | start data cost | Jacobian check, worst column |
|---|---|---|
| background bg1e-2 | 5.378360708129928e-03 | 4.1e-4 |
| logspread f0.1 + u0 | 7.420894283555133e-03 | 1.6e-4 |
| sos r3 | 1.618620951583266e-02 | 6.8e-8 |

- background: the 4.1e-4 is the softplus kink, as in the stage-1 tests.
- logspread: the 1.6e-4 is in multiplier columns at their bound whose bins carry e^-9.7 of the activity
  (finite-difference noise).
- logspread start: Var_w(u) = 3.597 (= the source's final), ū = −9.70, so (ū − u0)² = 61.3. The level penalty is
  0.0262 of the start penalty 0.0277: the stage-1 f0.1 end point had lowered the uniform level by 7.8 neper and
  compensated with the multipliers.

**Stage-1 continuations (submitted 2026-11-15T02:34Z = 2026-11-14 20:34 CST, `jobs/nonneg_alt_stage1c.sh`,
snapshot `nonneg_20261115T020314Z`, +4000 evaluations each; D with u0 = −1.864126788612443; in this order on any
GPU node but gpu003):**

| run root (`runs/shap8_012m_lategrid_dz1.5/nonneg/`) | continuation of | fit attempts, guard, summary |
|---|---|---|
| `background/stage1c_bg1e-2_20261115T023450Z` | `background/stage1_bg1e-2_20261114T084711Z` | 125950-125952, 125953, 125954 |
| `background/stage1c_bg1e-3_20261115T023450Z` | `background/stage1_bg1e-3_20261114T084711Z` | 125955-125957, 125958, 125959 |
| `logspread/stage1c_f0.1_20261115T023450Z` | `logspread/stage1_f0.1_20261114T090251Z` | 125960-125962, 125963, 125964 |
| `logspread/stage1c_f1_20261115T023450Z` | `logspread/stage1_f1_20261114T090251Z` | 125965-125967, 125968, 125969 |
| `logspread/stage1c_f10_20261115T023450Z` | `logspread/stage1_f10_20261114T090251Z` | 125970-125972, 125973, 125974 |
| `sos/stage1c_r3_20261115T023450Z` | `sos/stage1_r3_20261114T091354Z` | 125975-125977, 125978, 125979 |
| `background/stage1c_bg1e-1_20261115T023450Z` | `background/stage1_bg1e-1_20261114T084711Z` | 125980-125982, 125983, 125984 |

**Fine-tuning regularization (`nonneg_finetune.py`, user 2026-11-15 ~21:10 CST; every model, also option C through
`localnn_ft.py`).** Options, off by default (bitwise unchanged when absent):
- `--delta-penalty-f F --delta-cref C` (job variables `DELTA_F`, `DELTA_CREF`): a penalty on the total per-facet log
  factor t_e = h_e + δ_e (h = the fit's held term) over the tuned facets.
  - Definition: P_δ = f · C_ref · Σ_e w_e t_e² / Σ_e w_e, with w_e = facet area and C_ref = 4.2688066575e-03. An rms
    total log factor of 1 neper costs f × the current misfit.
  - Implementation: exact residual rows √(2 f C_ref w_e / Σw) t_e on the 0.5 r·r scale, with their own matrix-free
    Jacobian rows.
- `--delta-cap X` (`DELTA_CAP`): |h_e + δ_e| ≤ X as bounds δ_e ∈ [−X − h_e, X − h_e]. Where |h_e| > X at the start,
  δ_e starts at the nearest bound; the count is in fine_tuning.json `delta_regularization.start_clipped`.

Records (costs_initial / costs_final of the state and fine_tuning.json, and finetune_result.txt):
- cost_delta_penalty, delta_total_rms_area_weighted, delta_total_max_abs;
- facets_at_cap (upper / lower; within 1e-8 max(1, X));
- the settings;
- cost_total = data + crossing + P_δ;
- finetune_state.npz gains `upper_bound` under the cap.

Validation:
- check 1 compares data + crossing with the fit;
- check 3b compares the rows with FD and P_δ with an independent sum;
- under the cap the FD steps of checks 3 and 6 are not clipped.

Tests: CPU 56 of 56 (jobs 134500, 134502). Snapshots: `model_ext_frozen/nonneg_20261116T031115Z`, and
`nonneg_summary_20261116T031115Z` with the new `nonneg_ft_summary.py` (`jobs/nonneg_ft_summary.sbatch`). The summary
tool reports the late rms, P_δ, |h + δ|, the seen area holding 50 / 90 % of the seen production, and the tuned
facets' share of the late signal.

**B / D c5 tests (2026-11-16):**
- Setup: the c5 fine tuning only, rerun from the c5 fit state of the stage-2 chains `background/stage2_bg1e-2_20261115T184949Z`
  and `logspread/stage2_f0.1_20261115T184949Z`; f_δ = 0.1, 1, 10 with the cap ln 100; MAXNFEV 500; late_only, 5-min
  window.
- Runs: on CPU by the option-C agent's engine (jobs 134505 and 134501; `tests/localnn_20261114/ftreg/`).
- Engine check: its unregularized controls reproduce the chains' GPU c5 fine tunings (final data rel 6e-14 B, 4e-16 D).
- Summaries: `checks/ft_delta_reg_20261116/` (jobs 134508-134515).
- At c5 the held term reaches |h| = 12.2 / 12.3 neper; 36 / 37 of the 271 tuned facets start clipped to the cap.

Columns: the data cost; P_δ; late-window rms ln NG / RG; the area-weighted rms and the maximum of |h+δ|; the facets at
the cap; the seen area holding 50 % / 90 % of the seen production at 10:33:10; the tuned facets' share of the seen
production; their mean share of the NG / RG late signal.

| run | data | P_δ | late NG | late RG | rms / max \|h+δ\| | at cap | 50 % / 90 % | tuned share of production | late share NG / RG |
|---|---|---|---|---|---|---|---|---|---|
| B unregularized (chain) | 4.783e-3 | – | 0.139 | 0.040 | 2.49 / 12.2 | – | 1 facet / 0.0110 | 0.57 | 0.071 / 0.046 |
| B f_δ 0.1 | 4.928e-3 | 1.8e-5 | 0.150 | 0.033 | 0.21 / 4.61 | 2 | 1 facet / 0.0085 | 0.65 | 0.055 / 0.042 |
| B f_δ 1 | 4.954e-3 | 9.4e-5 | 0.151 | 0.035 | 0.15 / 2.63 | 0 | 1 facet / 0.0095 | 0.63 | 0.053 / 0.041 |
| B f_δ 10 | 5.239e-3 | 2.3e-5 | 0.167 | 0.045 | 0.023 / 0.95 | 0 | 0.0017 / 0.0200 | 0.14 | 0.034 / 0.044 |
| D unregularized (chain) | 7.214e-3 | – | 0.136 | 0.061 | 2.97 / 12.3 | – | 3 facets / 0.0234 | 0.74 | 0.134 / 0.068 |
| D f_δ 0.1 | 7.737e-3 | 4.4e-5 | 0.158 | 0.093 | 0.32 / 4.61 | 2 | 1 facet / 0.0201 | 0.75 | 0.069 / 0.033 |
| D f_δ 1 | 7.833e-3 | 1.3e-4 | 0.159 | 0.102 | 0.17 / 2.72 | 0 | 1 facet / 0.0238 | 0.74 | 0.066 / 0.032 |
| D f_δ 10 | 8.305e-3 | 5.7e-5 | 0.180 | 0.115 | 0.037 / 1.07 | 0 | 0.0033 / 0.1064 | 0.22 | 0.045 / 0.038 |

The single facet that holds 50 % of the seen production (except at f_δ = 10) is the tuned facet 12232930 (1.03 m²,
on the night side at 10:33:10). Its A is 2e8 A_orig, nearly all from D2's ft templates held at amplitude 1. The δ
regularization only scales it by exp(h + δ).

## Running on another machine

The copies are as run on Voyager. Change, in the files you use: `APP`, `AMPS`, `T`, `COPS_OUTPUT_ROOT`,
`COPS_DATA_ROOT`, `CUDA_PATH`/`CUDA_HOME`, the MPI path in `PATH`, and the `#SBATCH` partition, account,
`--nodelist` and resource lines. Then `./check_copies.sh` reports those files as DIFFERS, as expected; keep the
change in git so that the code of every new run is recorded.

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
