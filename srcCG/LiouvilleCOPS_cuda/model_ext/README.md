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
| a whole fit / fine-tuning chain as ONE job (records bound once; same products and layout as the split chain) | `jobs/chain_onejob.sh` (dry run; `--submit`) -> `jobs/chain_onejob.sbatch` -> `chain_onejob.py`; KIND positive (default), slope (one cycle only), nonneg | the stage modules of the chosen kind; tests and the split-vs-one-job compare `test_chain_onejob.py` | 2026-11-10 (positive tests 124091-124095, bitwise); 2026-11-11 (nonneg tests 124143/124148, production 124162/124165; frozen `chain_onejob_20261111T060742Z`) |
| one fit split across one or two GPU nodes (MPI, one rank per GPU) | `jobs/mpi_fit_submit.sh 1\|2\|auto <spec>` -> `jobs/mpi_fit.sbatch` (SPEC files in `mpi_multinode/specs/`); probe `jobs/mpi_probe.sbatch` | `mpi_support.py` (per-rank wrapper); comparisons `mpi_multinode/compare_tests.py`, `compare_slope.py` | 2026-11-10 (tests 123623-123649 only; the production fits ran on one node) |
| which configuration (nodes, CPU threads, one job per chain) is fastest | `PERFORMANCE.md` (measured 2026-11-10/11) | — | — |
| retries without starting a stage behind a failed upstream (guard job per stage) | `jobs/slope_cycles.sh`, `jobs/positive_cycles.sh` (from 2026-11-10); guards `jobs/slope_chain_guard.sbatch`, `jobs/positive_chain_guard.sbatch`; layout test `jobs/test_chain_deps.sh` | — | 2026-11-10 (dummy jobs 124035-124058; 3a chain guards 124068-124086; positive guard not yet submitted) |

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
| `chain_onejob.sh`, `chain_onejob.sbatch` | ONE job per fit / fine-tuning chain (KIND positive \| slope \| nonneg) | positive tests 123990, 124091–124095 (earlier script); nonneg tests 124143, 124148; production 124162 (softplus), 124165 (squared) |
| `pub_maps_figures.sbatch`, `pub_maps_v2.sbatch`, `pub_source_areas.sbatch` | full / version-2 figure sets (GPU; inputs from the environment) | 124492-124516, 124565-124595, 124677-124683, 124722 |
| `cops_movies_prep.sbatch`, `cops_movies_render.sbatch`, `cops_movies_encode.sbatch`, `cops_movies_qc.sbatch`, `cops_movies_diag.sbatch` | movies: preparation and validation, render on all GPUs, compose / encode, MP4 check, dark-surface diagnosis; 2026-11-13 also the fixed camera (`CAMERA=fixed`, `FROM_WORK` + step `copy`, step `validate_fixed`), the touchdown data (step `touchdown`), the flag (`TOUCHDOWN`) and `MOVIES_PY` | 124615-124635, 124707-124735, 124786-124788, 124863-124869 |
| `ng_source_map.sbatch`, `pan_candidates.sbatch`, `cops_movies_pan_stills.sbatch` | 2026-11-14 panorama movie: per-facet NG source sums (late and all NG observations; CPU), candidate final camera poses (1 GPU), candidate stills (plan `--pan-to`, mid-path and last frame; 1 GPU per pose); `cops_movies_prep.sbatch` `PAN_TO`/`PAN_START`/`PAN_END` and `cops_movies_encode.sbatch` `NO_SUB_MARK` (defaults unchanged) | 124913, 124915, 124918, 124919, 124920; production 124930, 124931, 124936, 124937 |
| `fov_share.sbatch`, `source_area_compare.sbatch` | signal share outside the IK FOVs; source areas of two fits (CPU) | 124731, 124727; c9: 124789, 124784, 124785 |
| `mpi_fit.sbatch`, `mpi_fit_submit.sh`, `mpi_probe.sbatch` | one fit over 1 or 2 nodes with MPI (SPEC in `mpi_multinode/specs/`, NODES=1\|2\|auto); submit helper; MPI probe | tests 123629–123649 (earlier script versions; the current one not submitted); probe 123623 |
| `grain_rough.sbatch`, `prepare_12m_rough.sbatch` | 12M rough-surface base and the C30 fit on it | 122033, 122034 |
| `prepare_44m.sbatch`, `grain_44m_c30.sbatch` | 44M preparation; C30 with the 12M conditions on 44M | 121996/122020, 122078 |
| `refine_12m.sbatch`, `refine_20m.sbatch` | local mesh refinement tests | 122009, 122011 |

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
`cops_movies.py` + `cops_movies_geom.cpp` (2026-11-12; observable-part masks, spacecraft shadow and the dark-surface diagnosis 2026-11-13) make the one-minute movies of the descent seen from RG / NG (`12M/movies_softplus_5cycles/`, `12M/movies_softplus_smm1_5cycles/`); `fov_share.py` gives the share of the modelled signal from outside the IK fields of view. Since 2026-11-13 (second set) `cops_movies.py` also has the fixed camera (`plan --camera fixed`: frozen in the body frame at the COPS position and attitude of 100 m altitude), the gauges' line of sight from the true gauge position (depth map + exact edge refinement; `validate_fixed`), `touchdown` (the ground-track impact point, its visibility, SPICE check) and `compose --touchdown` (a waving chequered finish flag): `12M/movies_softplus_smm1_9cycles/` (gauge view) and `12M/movies_softplus_smm1_9cycles_fixedcam/`; and (2026-11-13/14) circular fields of view (`cone:<deg>` masks), `contact`, `prep --extend-to-utc`, `sc --hga-ck`, `plan --camera detach` and the detach render with a spacecraft model (bus, wings, HGA), its shadow and the FOV cones (`cops_movies_geom hgageom`): `12M/movies_softplus_smm1_9cycles_detach200/`; `fov_share.py --cone-rg/--cone-ng` gives the share outside circular fields of view. 2026-11-14: `plan --camera detach --pan-to <pose.json>` (the camera moves after the detach to a final pose: smootherstep in the frame index, arc around the pivot, quaternion slerp, focal length geometric) and `compose --no-sub-mark` (no sub-spacecraft +); `ng_source_map.py` (per-facet NG source sums), `ng_source_regions.py` (lobes, shares by lobe / distance / azimuth, relief, map), `ng_lobe_share_vs_altitude.py` (big-lobe share per late observation) and `pan_candidates.py` (candidate final poses): `12M/movies_softplus_smm1_9cycles_detach300_pan/` (the panorama movies with the user's pose B, held from frame 1800; analysis and candidate stills).
`figures_12m_nonneg9_smm1.sh` (2026-11-13) draws the same bounded chain after 9 cycles (jobs 124664, 124772) vs its c5 state, the unbounded softplus and D2, and the cycles c1-c5 / c5-c9, into `12M/model_vs_observations_softplus_smm1_9cycles/`.

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
