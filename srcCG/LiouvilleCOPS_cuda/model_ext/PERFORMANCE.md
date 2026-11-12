<!-- ==== version-controlled copy: header added 2026-11-11, not part of the code that ran ====
Run copy: $COPS_OUTPUT_ROOT/model_ext/PERFORMANCE.md (on Voyager: /data001/highresolu_vt/cops_output/model_ext/PERFORMANCE.md)
SHA256 of the run copy = SHA256 of this file without these header lines: b18e042ce600659c34fa03280e9a84177c474dae6880b9d3887cc41065d6e54a
Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
WHAT: document. Measured run configurations and speed of the 12M fits on Voyager (user 2026-11-10: "document which
  configuration results in the most efficient run (the number of nodes, CPU threads)"): the recommendation per model
  family (slope / bump models: 1 node, chunked L step with 32 threads; clipped models: MPI 8 ranks with device activity;
  a whole chain as one job), the timing tables of the slope-multiplier model and of D2 (local / MPI), and the one-job
  section (merged 2026-11-11 01:00 UTC).
USED BY: read by hand (the run choice of the 2026-11-10/11 fits); measurements from jobs 123623-123649, 123891, 123979,
  123986-124105 (named in the document).
==== end of header ==== -->
# Run configuration and speed of the 12M fits (Voyager, 2026-11-10)

User request 2026-11-10: "document which configuration results in the most efficient run (the number of nodes, CPU
threads)". All numbers are measured on gpu001/gpu002 (8 x A100 80 GB, 64 CPU cores, 240 GB, 8 NUMA nodes per node,
200 Gb/s InfiniBand between them), 12M SHAP8 lategrid base (1551 points, 4.74e9 kernel records, mixed precision, no
record store), the clipped refit D2 and the slope-multiplier model; per-evaluation times are means from the fits'
evaluation logs (31 evaluations, 5 Jacobians each). Jobs: 123635, 123646, 123648 (slope model 1 node / MPI 8 / MPI 16),
123891 (threaded OpenBLAS), 123979 (chunked threads), 123629-123632, 123638 (D2 local / MPI), 123623 (MPI probe).

## Recommendation

| Model family | Most efficient configuration | 4000-evaluation fit | Results |
|---|---|---|---|
| slope multipliers (and the bump model built on it): host activity, L = basis x weights on the CPU | **1 node, local transport, chunked L step with 32 threads** (`L_THREADS=32`, chunked mode) | **about 30-35 min** (1 thread: about 55-65 min) | bitwise equal to 1 thread by construction |
| clipped models (D2 type, fit_cops' own activity) | 1 node, `--transport mpi`, 8 ranks, `--activity device`; 2 nodes (16 ranks) another 1.57x | pressures 29x and Jacobians 4.4x faster than the local transport (1 node) | equal to the local fit to rounding only (device summation order); `--activity host` is bitwise but gains nothing |
| positive (exp) model with a held term | 1 node, local transport (MPI refused with a held term) | about 1 h per fit (its L step is also single-threaded; threads not implemented there) | - |
| late-only per-facet fine tuning | 1 node; inside a one-job chain (positive model) | 3.5 min as its own job (bind ~2.3 min); ~1 min inside a one-job chain | one-job = split, bitwise |
| a whole positive fit / fine-tuning chain | 1 node, ONE job (`jobs/chain_onejob.sh`) | ~12 min less per 3-cycle chain than 6 jobs | bitwise equal to the split chain |

Two nodes do NOT help the slope / bump models as implemented: the fitter forms the field L = Y a (12.4 M x 25-32) on
the CPU and broadcasts 98.9 MB to every rank each evaluation; across nodes that costs more than the halved GPU time
(GPUs take only 0.04-0.05 s per evaluation). Computing L on the GPUs (the basis is already resident there) would make
2 nodes about 1.6x faster than 1; not implemented (user 2026-11-10: "One node now").

## Measured: slope-multiplier model (srcCG start, NFEV 5, CHECK 1)

| Configuration | s / evaluation | s / Jacobian | Bitwise vs 1 thread, 1 node |
|---|---|---|---|
| 1 node, local transport, 1 thread (123635) | 0.351 | 0.489 | reference |
| 1 node, MPI 8 ranks, host activity (123646) | 0.369 | 0.507 | yes |
| 2 nodes, MPI 16 ranks, host activity (123648) | 0.465 | 0.559 | yes |
| 1 node, threaded OpenBLAS 8 (123891) | 0.232 | 0.352 | in the test; not guaranteed (2-4 of 12.4 M L values differ by <= 7e-16 rel.) |
| 1 node, threaded OpenBLAS 16 | 0.154 | 0.283 | as above |
| 1 node, threaded OpenBLAS 32 | 0.140 | 0.279 | as above |
| 1 node, chunked 16 (123979) | 0.209 | 0.340 | yes, by construction |
| **1 node, chunked 32 (123979)** | **0.161** | **0.297** | **yes, by construction** |

L step alone (micro-benchmark, 123979): 1 thread 0.238 s (10 GB/s); threaded OpenBLAS 8/16/32/64: 0.111/0.109/0.095/0.092 s
(saturates at ~26 GB/s); chunked 8/16/32/64: 0.081/0.081/0.064/0.051 s (up to 48 GB/s), bitwise equal at every count.
Inside a fit the chunked pool gains less than in the benchmark; 48/64 threads are to be timed in the next GPU test.
Threads apply to the L step of the fitter process only; the 8 GPU worker processes keep OMP/OPENBLAS/MKL_NUM_THREADS=1
(user 2026-11-10 allowed this departure from CLAUDE.md §4 for the fitter process).

## Measured: clipped refit D2 (5 evaluations from c6)

| Configuration | s / pressure evaluation | s / Jacobian | bind (s) | vs local |
|---|---|---|---|---|
| 1 node, local transport | 0.551 | 0.672 | 158 | reference |
| 1 node, MPI 8 ranks, device activity | 0.019 | 0.151 | 130 | rounding (<= 1.7e-15) |
| 2 nodes, MPI 16 ranks, device activity | 0.012 | 0.096 | 77 | = 8 ranks bitwise |
| 2 nodes, MPI 16 ranks, host activity | 0.545 | 0.619 | 91 | bitwise |

The local transport spends ~0.43 s per evaluation forming the activity field on one CPU thread; device activity forms it
on the GPUs. Host memory: peak 36-67 GB per node of 251 GB; inter-node MPI (srun --mpi=pmix, UCX): 24.7 GB/s, 2.2 us.

## One job per chain: bind the records once (option (b), 2026-11-10)

Status: tested 2026-11-10 (jobs 123986-123989 split, 124091 one-job, 124092 recheck, 124093-124095 resume; compare 124099, 124105).

**Configuration.** A whole fit / fine-tuning chain of one lambda (c1 fit -> c1 ft -> ... -> c3 ft) runs as ONE job on
one node (8 GPUs, local transport, host activity, FP32 records), `jobs/chain_onejob.sh` -> `jobs/chain_onejob.sbatch`
-> `model_ext/chain_onejob.py`. The driver binds the 4.74e9 kernel records once and runs every stage in the same
process on that backend; before each later stage only the template basis is reset on the kept records (~1.5 s).
Products, layout and stage settings are those of the split chain (`jobs/positive_cycles.sh`, one job per stage).

**Where the split chain loses time** (ratio 1e-5, lambda 1e-4, jobs 123035-123052, sacct + logs):

| stage | job | elapsed | bind |
|---|---|---|---|
| c1 fit | 123035 | 5037 s | 161.3 s |
| c1 ft | 123038 | 195 s | 129.1 s |
| c2 fit | 123042 | 3784 s | 161.0 s |
| c2 ft | 123045 | 212 s | 142.8 s |
| c3 fit | 123049 | 1735 s | 133.5 s |
| c3 ft | 123052 | 222 s | 148.9 s |
| total | 6 jobs | 11185 s (186.4 min) | 876.6 s (14.6 min, 7.8 %) |

The bind is 66 % of every fine-tuning job (its optimization takes ~15 s) and 3-8 % of a fit job.

**One job:** 1 bind (~130-160 s) + 5 basis resets (~1.5 s each) per 3-cycle chain: about 12 min less per chain
(~7 % of this chain; more for chains whose fits end early, e.g. ftol); the 5 Slurm hand-offs between dependent jobs
also disappear (about 1 s each here, longer when other jobs take the node in between). Measured on the 2-cycle test
chain (ratio 1e-5, lambda 1e-4, fit NFEV 5, fine tuning MAXNFEV 3): split 4 jobs 1322 s with 4 binds 616 s (140-169 s each);
one job 851 s (-36 %) with 1 bind 166 s and basis resets of 0.12-0.13 s; a fine-tuning stage takes 51-57 s in the one job vs
203-209 s as its own job. A resubmission of the completed chain skips every stage (6 s, exit 0).

**Results.** The one-job chain runs the split chain's commands; the test compares split vs one-job bitwise
(fit_result.txt numbers, every fit_restart.npz / activity_exp_facets.npz / finetune_state.npz array, COPS_model.dat,
late_window_model.dat, fine-tuning cost history and grain). Result (124099): c1 fit, c1 ft, c2 fit, c2 ft BITWISE IDENTICAL (checked) to the split chain.
Resume (124093-124095): a job stopped by its time limit ends the running fit at its best evaluated state (exit 75) and
the next attempt continues it as a new segment and completes the chain; like every segmented fit, the result is not
bitwise that of an uninterrupted run (here: final data cost +0.9 %, fits at max_nfev).

**Use.** `jobs/chain_onejob.sh --lambda L --ratio R --ft-maxnfev 500 --ft-facets late_only --ft-window-min 5
--model-ext <frozen> [--source ...] --node gpu002 --time 08:00:00 --attempts 3` (dry run; `--submit`). A chain longer
than the job's limit continues in the next attempt (afternotok): every stage stops at its best state before the
limit (exit 75) and the next attempt skips the completed stages. Not for the slope / bump (3a) models yet: their GPU
patch carries the cycle's held term in its frozen worker configuration (see model_ext/logs/onejob_stage_interface_3a.md).

## How to run

* One node, chunked threads: `sbatch --export=ALL,MODEL_EXT=<snapshot>,L_THREADS=32,L_MODE=chunked jobs/slope_mult_fit.sbatch`
  (exact variable names: the job script's header and model_ext/logs/agent_slopemul_READY.txt).
* Two nodes / auto: `jobs/mpi_fit_submit.sh 2|auto <spec> ...` (model_ext/logs/agent_mpi_READY.txt).
