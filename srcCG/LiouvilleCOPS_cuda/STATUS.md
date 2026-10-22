# Stage status

Stage records follow ROADMAP §5. "Pending" means the evidence does not exist yet. A
synthetic success is never counted as production evidence.

## Plan revision 2 (2026-10-19)

`ROADMAP.md` now incorporates `AMMENDED_ROADMAP.md`; the review comments are in
`AMENDMENT_REVIEW.md` and the original plan is `ROADMAP_v1.md`.

* **Primary target:** the free B29 baseline and the N30/C30/NC31 template families. The
  fixed-baseline degree-5–8 harmonics are secondary.
* **Oscillation diagnostics:** added as reported, non-gating quantities (ROADMAP §3.4).

Reconciled facts, checked 2026-10-19:
* The CPU tree was repaired at 13:02–13:06:
  * `gravity_topology.py` audit;
  * DISABLE in both gravity paths;
  * cache algorithm `…_v2_disable_after_topology_audit`;
  * the joint fitter;
  * `prepare_restored_campaign.py`.
* The five HEAL products are quarantined as `*_heal_defect`.
* The L5 harmonic scan has finished (optimizer state converged, 17:24 UTC); it is a
  secondary-mode result.
* The C02 contract passes a synthetic NC31 configuration (31 columns, joint
  eigen-regularization) bitwise.
* New pending items: C00b re-freeze against the repaired tree, C02 joint extensions, and
  the C04 gravity acceptance route (a user decision).

## C15: scaling to SHAP8 20M/44M and the Voyager cluster (2026-10-21, user request): P1, P2 (with P2b), P3 COMPLETE

Plan and decisions: `SCALING_PLAN.md` (§0a). Phase P1 (data layer), all bitwise-neutral:
* **Template fits keep their records:** `CudaBackend.rebind_basis` and `MultiGpuBackend.rebind_basis`; the registry
  rebinds only the evaluator's basis-sized arrays (event `rebind_basis_only`). A fresh bind before cost 12 min per
  12M template fit. Tests: `test_gpu_model.py::test_basis_rebind_is_bitwise_a_new_backend` (FP64, FP32),
  `test_multi_gpu.py::test_basis_rebind_on_two_devices_is_bitwise`.
* **Face compaction (user decision):** `DeviceDataset.load(compact_faces=True)`, on by default in `CudaBackend`
  (`COPS_CUDA_COMPACT_FACES=0` turns it off). The device holds per-face data only for the faces in its records;
  `kernel.face` is remapped on the device; `gather_faces`/`scatter_faces` translate at every kernel call site
  (evaluator, masks, restart, fine tuning). Measured on 12M: 74.5% of the faces occur in any record; with 8 GPUs the
  largest rank needs 52%. Real-data check, 1M full table: 840,458 of 1,194,250 faces on the device; pressures and
  point changes bitwise equal to uncompacted. `tests/test_face_compaction.py` (5).
* **FP32 record store (user decision Q4):** `cops_cuda/record_store.py`, `tools/c15_record_store.py`; one file per
  point (20 B/record), manifest tied to each kernel file (size, mtime, SHA-256) and to `run_manifest.json`; the FP32
  conversion is one shared function (`device_store.fp32_record_fields`). Loader: `record_store=` (reader threads,
  pinned double buffer); `gpu_fit.py --record-store`. 3M: built in 91 s (17.1 GB); half the rows loaded in 19.8 s
  against 278 s from the kernel files, device arrays bitwise identical. `tests/test_record_store.py` (6).
* **Validate once (user decision Q4):** `base_validation.AuditCache` (`gpu_fit.py --audit-cache auto|off|PATH`,
  `--revalidate`): unchanged files (size, mtime, inode) reuse their first full validation; the reference
  `validate_fit_base` still runs on these results. `test_base_validation.py::test_validate_once_cache`.
* 216 tests pass.

Phase P2 (streaming engine), FP32 records: COMPLETE for fits.
* `CudaBackend` keeps as many records resident as the budget allows (per-face arrays, slots and partial sums
  planned first; aligned to a block end; `COPS_CUDA_MAX_RESIDENT_RECORDS` caps it). The rest stay in pinned host
  memory (`DeviceDataset.host_records`). Each pass runs the resident blocks, then copies chunks
  (`COPS_CUDA_STREAM_CHUNK_RECORDS`, default 2^25 records) into 3 device slots on a copy stream while earlier chunks
  compute, and runs every tile of the request on a chunk while it is on the device (one transfer per pass). The
  kernels see block for block the same records, so results are **bitwise** the resident ones.
* `seen_faces`, audit and device digests cover the streamed part. Restart/fine-tuning kernels refuse a streamed
  backend (`UnsupportedConfiguration`) until P2b.
* `tests/test_streaming.py` (5): caps 0..nrec-1 with 40-record chunks, compaction on/off, record store, two GPUs,
  all bitwise; 221 tests pass.
* 3M, one A6000 (PCIe gen3, 12.3 GB/s), from the store: resident P 0.214 s, J 1.049 s (16.7 GiB); 50% streamed
  P 0.747 s, J 1.159 s (10.6 GiB); 88% streamed P 1.28 s, J 1.311 s (4.5 GiB); all 31 result arrays bitwise equal.
  Pressure passes are transfer-bound; the Jacobian hides most of the transfer.
* 20M chain started (`output_fit_shap8_020m_base` in /tmp/vtenishe; store, gravity, bundles, fits on /data;
  fits with `--headroom-fraction 0.12 --precision mixed --record-store`).
* 20M base: 5.06e9 records (253 per facet, 28% below the projection), 194 GB, 6271 s. Store: 101.1 GB, 354 s.
  The first chain stopped at the visibility masks: their per-facet scratch (uint32, 80 MB at 20M) was outside every
  plan and above the 64 MB slack. `memory.scratch_reservation` now reserves such arrays for the call (up to half
  the slack they stay in it); `tests/test_masks_budget.py` (2). Chain restarted from that step.
* The second chain lost B29 at 22:12 (2026-10-21) to the kernel OOM killer: other jobs had taken about 190 GB of the
  251 GB host RAM. The GPU workers now memory-map `elements.bin` (`formats.map_elements`; RunInputs loads lazily,
  `formats.LazyElementFields`), so a worker holds about 1.3 GB instead of about 15 GB; a 20M fit needs about 55 GB of
  host RAM (11 GB parent, 2 x (1.3 GB + 20.2 GB pinned)). The third chain (user decision: wait for free RAM) starts
  each fit after 10 min with at least 110 GB available, with `oom_score_adj` 1000 so the fit, not another job, is
  killed first.
* **20M fits complete (2026-10-22 13:42)**, `fit_020m_gpu_{B29,N30,C30,NC31}_n2000_mixed`, 3.3–4.0 h each on
  2 A6000 (2.3–2.5 s per pressure call, 4.1–4.4 s per Jacobian, 30% of the records streamed). Comparison:
  `cops_output/comparison_020m_vs_012m_003m_001m/` (README, cost table, figures). B29 0.066997, N30 0.039970, C30
  0.068048: within 0.06–0.3% of 12M, same branches. **NC31 0.033252 is degenerate:** Tmin (87 K) rose above Tmax
  (77 K) and T75 (44 K), so T = max(Tmin, lit law) = Tmin on every facet; Tmax and T75 have zero Jacobian columns
  (data Jacobian rank 29 of 31, condition 1.6e22; their printed sigmas of about 0 are the pseudo-inverse dropping
  them). The optimizer has no gradient out of that region; 12M NC31 reached 0.026600. Rerun from another start and
  an optional Tmin <= T75 <= Tmax constraint: pending user decisions.

Phase P2b (streaming for the restart and fine-tuning kernels): COMPLETE.
* `CudaBackend.record_chunks()` yields the bound records in block order (resident part, then each streamed chunk
  copied into a slot) with the matching block table; every `kernels/restart.cu` call site (facet shares, point
  changes, first structural/significant row, first-interval scores) loops over it. Their outputs are per-facet
  atomic max/min, per-block partials indexed by global block, or chronological per-point sums (one record per facet
  and point), so the results are **bitwise** the resident ones (`tests/test_streaming.py`: caps 0, 50 and 130
  records with 40-record chunks, points split over chunks and over the resident boundary).
* Their per-facet scratch is reserved (`memory.scratch_reservation`), and `StateEvaluator` backends keep 32 B per
  device facet free when sizing the resident records (`scratch_per_face`). Fit backends keep 1 B (seen_faces), plus
  8 B with device activity (the check).
* So `restart_tool.py enhance|translate --precision mixed` and `fine_tune_gpu.py --precision mixed` run on bases
  larger than the GPUs (12M, 20M), streaming the excess. 239 tests pass.

Phase P3 (distribution over MPI): COMPLETE (local `mpirun`; the cluster gates wait for an account).
* `cops_cuda/transport.py`: `MultiGpuBackend` talks to its workers through channels with one interface: spawned
  processes (`transport="local"`, unchanged behaviour), MPI ranks (`transport="mpi"`, mpi4py 4.0.3 built against the
  system Open MPI 4.0.3), and rank 0's own GPU in process (its block runs while the other ranks compute). One worker
  handler serves bind, evaluate, rebasis, seen, identity, close (the rank keeps serving later binds) and exit.
  Template bases travel as chunked raw buffers (below MPI-3's 2^31 count). After a failure, closing drains the
  unanswered replies, so the next bind on the same ranks is clean. A rank's device is its launcher local rank
  (Open MPI, Slurm, MPICH, MVAPICH variables; no collective) modulo the visible devices.
* `tools/gpu_fit.py --transport mpi [--activity host|device]`: every rank starts the tool; ranks >= 1 serve, rank 0
  fits. A failure on rank 0 ends the job with `MPI_Abort`; at the end rank 0 releases the workers.
* **Activity split over the GPUs (user decision Q3):** `CudaBackend(activity="device")` uploads Y and the raw template
  columns of its own facets once per bind and forms act_lin, af and the active set on the device
  (`kernels/activity.cu`, fixed term order, NumPy's clip). Deterministic and the same in every partition; it is not
  the host BLAS rounding, so it is checked against the host expression at the first evaluation and every 100th
  (error ≤ 64 eps of Σ|terms|, no active-set flip outside that bound; `BackendError` otherwise; recorded in the
  identity). Default of the MPI transport; `host` (exact, `MPI_Bcast` of the field) stays available and is the
  local default.
* Partition: contiguous row blocks balanced by records. At 20M the largest block exceeds the ideal by 2.1% at
  24 ranks (the largest point is 3.4% of a rank's share), so the block-granular split is deferred to P4.
* Tests: `tests/test_mpi_transport.py` (8, real `mpirun -n 1..3` jobs, ranks may share a GPU): host activity 2 ranks
  and device activity 1–3 ranks, FP64 and streamed FP32, with a basis rebind, `seen_faces` and a no-basis (B29) bind,
  all BITWISE one GPU; a worker failure is a `BackendError` naming the rank, and the ranks serve the next bind.
  `tests/test_activity_device.py` (5): within the bound, BITWISE across compaction, streaming, two GPUs and a rebind,
  a corrupted device basis is caught. The existing multi-GPU tests pass unchanged on the refactored pipe transport.
  Full suite: 236 tests pass (239 with P2b).
* **Real-data gate, 788k B29** (FP32 records, max_nfev 4): `mpirun -n 2 gpu_fit.py --transport mpi --activity host`
  is BITWISE the local `--devices 0,1` fit (x, a, T, κ, act_lin, predicted NG/RG). With `--activity device`: x within
  8.8e-15 relative, the same final cost (0.8724887), predicted pressures bitwise; activity check 3.16 eps of the
  term sum (bound 64), 0 active-set flips. (The raw ulp distance is large only where the terms cancel to near
  zero, which is why the bound is relative to Σ|terms|.) Each fit 1.5 min wall, 75 s of it the bind.
* Not yet: distributed base validation (validate-once makes it a one-time cost: rank 0 alone); per-rank element
  subsets (every rank reads `elements.bin`: 4 GB of Y per rank at 20M); the Slurm template
  `tools/cluster/fit_mpi.slurm` has SITE placeholders until the Voyager account exists.


Phase P0 preparation (2026-10-22, user request; CPU-tested; developed in a copy while the 20M chain ran, merged into
this tree on 2026-10-22 after the chain, full suite 267 passed including the mpirun tests): the campaign runner for Voyager.
* `tools/cluster/cops_pipeline.py plan|run|submit|status`: the user's five settings (`--shape`, `--surface-bc
  flat|rough`, `--restart-dir`, `--max-nfev`, `--fine-tune`) plus fits B29/N30/C30/NC31 (templates as in the 20M chain),
  `--restart-from FIT=FILE` (template and mesh checks; `--translate-restart` for another mesh), existing products
  (`--base`, `--store`, `--visibility-bundle`, `--bundle`, identity-checked: rebuild.json, run_manifest.json,
  store_manifest.json, bundle sidecars). Stages base → store → visibility → gravity → cache → bundle → fits →
  fine tuning, each with a stamp (inputs, command, code digests, hardware); failed or killed stages are moved aside
  (`.incomplete_<UTC>`) and run again; another mesh, emission law, input, gravity body or precision request is refused.
* Hardware is detected at run time (Slurm variables, CUDA_VISIBLE_DEVICES, nvidia-smi or the CUDA device properties,
  /proc/meminfo; overrides `--nodes --gpus-per-node --gpu-memory-gb --host-memory-gb`); nothing is hardwired. Precision
  `auto` (user decision): FP64 when the FP64 records fit resident on the fit GPUs, else mixed with streaming (streaming
  stays mixed-only, user decision); the pinned host buffers must fit the node (else refused). The memory model
  reproduces the local 20M fit: 3.56e9 records resident / 1.50e9 streamed / 54.9 GB host predicted, against 3.54e9 /
  1.52e9 / about 54 GB measured.
* `--gravity-precision mixed|fp64` (default mixed, user decision 2026-10-22), separate from the record precision; both
  policies are recorded in `pipeline.json`. `plan` estimates the gravity time from recorded c04 runs (20M on 8 GPUs
  of the A6000's rate: mixed about 12 min, fp64 about 70 h), or says unknown.
* Wall time (user decision 2026-10-22; default none): `gpu_fit.py --wall-time HH:MM:SS|slurm | --deadline` stops the
  unchanged fitter's `least_squares` at the best evaluated state before the deadline (`cops_cuda/wall_time.py`),
  writes the restart file with `status = "interrupted"` and a `segments` record, exits 75; the runner continues the
  fit as segment k+1 with `--restart-from` and the remaining evaluations. A segmented fit is not identical to one
  continuous run. `submit --resubmit N` scripts resubmit themselves (afterany) while `status --check` says resumable.
* `submit` writes the task/GPU lines from `--nodes` and `--gpus-per-node` (default `auto`: `sinfo -h -o %G`);
  `tools/cluster/fit_mpi.slurm` no longer hardwires `--ntasks-per-node=8` (given at sbatch time).
* Tests (CPU only, fake tools / fake clock): `tests/test_cops_pipeline.py` (20), `tests/test_wall_time.py` (7);
  `tests/test_gpu_fit.py` unchanged and passing. `plan` read-only against the real 20M products: base, store and both
  bundles accepted as external, the four fits and fine tunings pending, nothing written (8 s, 67 MB).
* PENDING (after the 20M chain, GPUs needed): a real `run` of the runner on this machine (788k, then 20M with
  `--wall-time` segments), a wall-time stop of a real GPU fit, MPI transport through the runner, and the cluster gates.

## C14: mixed precision and the SHAP8 12M mesh (2026-10-21, user request): EXPERIMENT; four fits complete

User request: "switch the gravity calculation and the fitter to single precision as an experiment", "12M, start with
gravity", then "run the full fit of the 12M surface model with mixed precision"; decisions "Admit genus 4
(experiment)", "Yes, implement and run". The FP64 paths are unchanged and stay the validated default.

**Mixed gravity** (`GravityEvaluator(precision="mixed", near_factor=6)`, policy
`mixed_near_fp64_exact_far_fp32_quad3_v1`): faces with |c − P| < 6 circumradii use the exact FP64
`face_contribution`; far faces use an FP32 3-point quadrature of ∫dS/r; FP32 sums of 64-face tiles are added to FP64
running sums (launches aligned to whole tiles, so results do not depend on the launch split).
`tests/test_gravity_mixed.py` (2). 1M: 28 s against 2 h FP64, slope difference ≤ 1.1e-5°, cliff masks identical;
3M: 117 s against 7 h, masks identical.

**Mixed fitter** (`gpu_fit.py --precision mixed`, `COPS_CUDA_RECORD_PRECISION=fp32`, policy
`fp32_records_fp64_arithmetic_v1`): Kn, KPhi, KV and cos stored in FP32 (shadow folded into cos = −2), 20 B/record,
no stored static factors (`static_factors_m` recomputes them in FP64), fused kernels only, FP64 sums. The loader
refuses a value FP32 cannot hold (zero-ness and finiteness must survive) and steps cos back one ulp where FP32
rounding would cross COS75. `tests/test_evaluator_mixed.py`: pressures and every Jacobian column ≤ 1e-6 of the
summed magnitudes. 1M B29 full fit: cost relative difference 1.5e-7, parameters ≤ 8.8e-5 σ (nfev 1892, 1582 s).

**Genus admission** (`GravityBody(allow_genus=g)`, `--allow-genus` in c04/c05 tools): SHAP8 12M has genus 4
(χ = −6); a reference-audit refusal that is due only to genus is tolerated when g matches.

**SHAP8 12M** (`cg-spg_mspcd-shap8-v0.4cmb-12Mfacets.bdf`, 12,364,290 facets, Cheops frame; against SHAP5 3M:
median 3.8 m, p99 44 m, max 103 m):
* Base: 4,334,777,197 records, 173.4 GB, 3906 s (GPU rebuild).
* Gravity (mixed, 2 GPUs): 1065 s, 1.44e11 pairs/s; cache offset-halving max 0.0056°; bundle written (late
  visibility, cliffs ≥ 45°/60°, ground track).
* Fits B29, N30, C30, NC31 (mixed, 2 GPUs, headroom 0.06, 44.9 GB per GPU): bind 715 s.
* **B29** (`fit_012m_gpu_B29_n2000_mixed`): 7520 s wall (bind 715, pressures 1752 for 2004 calls, Jacobian 4758
  for 1557), nfev 2000 (limit). Cost 0.066789, rms NG 0.3218, RG 0.1732, D 15.9; Tmax 30 and T75 400 at the bounds,
  Tmin 40.8, kappa 0.820. The same branch as 788k (0.067692), 1M (0.068711) and 3M (0.067207); active area 15.2%
  (3M 21.3%), unseen activity fraction 0.94 (3M 0.91).
* **N30** (`fit_012m_gpu_N30_n2000_mixed`): 9630 s wall, of which 1390 s is two binds (the fitter's template-free
  consistency evaluation, then the template basis: `rebind_for_new_basis`), nfev 2000 (limit). Cost 0.039995, rms NG
  0.2104, RG 0.1891, D 7.3; c_late_ng_not_rg 7.25 ± 7.1 (1M 6.99, 3M 7.51). **Different temperature branch:** Tmax
  85.3, T75 30 (bound), Tmin 51.2, kappa 0.133, against 1M/3M Tmax 30 (bound), T75 223/206, Tmin 70.0/70.5, kappa
  0.096/0.086 (temperature sigmas 24–116 K; the branches are weakly separated by the data).
* **C30** (`fit_012m_gpu_C30_n2000_mixed`): 9947 s wall (binds 1392 s), nfev 2000 (limit). Cost 0.068196, rms NG
  0.3210, RG 0.1826, D 27.2; c_cliff_ge_60deg 6.71 ± 6.8 (1M 0.22, 3M 0.49); Tmax 30 (bound), T75 382, Tmin 49.0,
  kappa 0.778 (B29 branch). **Not converged:** the cost is above the 12M B29 cost (0.066789) although C30 contains B29
  (c = 0); 3M shows the same (C30 0.067404 > B29 0.067207). 1M C30 0.068664.
* **NC31** (`fit_012m_gpu_NC31_n2000_mixed`): binds 1104 + 1168 s (base partly re-read from disk), optimizer
  13:28–15:50, nfev 2000 (limit), first-order optimality 1.8e-4. **Cost 0.026600, rms NG 0.1817, RG 0.1421: the lowest
  of all fits on any mesh** (12M N30 0.039995). c_late_ng_not_rg 59.5 ± 29.2, c_cliff_ge_60deg 5.84 ± 1.88; Tmax 144,
  T75 30 (bound), Tmin 30 (bound), kappa 0.146; D 33.1, active area 21.8%. 1M/3M NC31 stayed on the B29 branch (cost
  0.068664/0.067368, c_late 1.14/1.12 with σ = 0: the template sat in the clipped region, no gradient).
* Comparison figures: `cops_output/comparison_012m_vs_003m_001m/` (README.md). In the late window 12M N30 and NC31
  raise NG toward the data but drive RG ~4x below the flat RG data below 30 m.

**Two fixes found on the 12M run:**
* `multi_gpu`: the 600 s bind timeout (sized for 788k) failed on 2.2e9 records per worker; the bind allowance is now
  600 s + 1 µs/record (3M binds at 0.23 µs/record).
* `cops_cuda/base_validation.py`: the reference `validate_fit_base` unpacks every record in Python (1.8 h per 12M
  fit). `gpu_fit.py` now runs the reference function with its two per-byte workers swapped for NumPy equivalents
  with threaded prefetch (`--reference-validation` restores the original). `tests/test_base_validation.py` (23): the
  audit is identical and every fault (range, duplicate, non-finite, shadow, header, elements; several faults per
  file) gives the identical first error. 788k real base: identical audit, 6.5 s against 349 s. `cops_io.py` and
  `fit_cops.py` are unchanged.

**Fine tuning and restart evaluation with FP32 records** (so 12M restarts can be fine-tuned): `kernels/restart.cu`
`_m` variants of facet_point_shares, point_change_partials, first_interval_scores and first_significant_row recompute
the H2O factors in FP64 from the FP32 records as `static_factors_m` does; `first_structural_row_f32`.
`restart.record_factor_args` selects them by `backend.mixed`. `fine_tune_gpu.py --precision mixed` and
`restart_tool.py enhance|translate --precision mixed` (one GPU, so not for 12M). Test
`test_fine_tuning_gpu.py::test_fp32_records_match_fp64`: first rows identical, scores, shares and point changes
within 1e-6 of FP64. 201 tests pass. The fine tuning has not yet been run on 12M.

## C13: late-interval fine tuning on the GPU (2026-10-21, user request): IMPLEMENTED

The CPU branch added `tools/fine_tune_activity.py` (`../LiouvilleCOPS/FINE_TUNING_STEP.md`, commit 0ff50ac,
not this work): a second fit on a completed restart, changing only facets first seen in the last m minutes.

**Code** (CUDA tree only): `tools/fine_tune_gpu.py`, `cops_cuda/fine_tuning.py` (`SplitEvaluator`: contiguous row
blocks on 1..N GPUs), `kernels/restart.cu` (`first_structural_row`, `first_interval_scores`,
`first_significant_row`), `StateEvaluator(row_indices=...)`, `DeviceDataset.kernel_entries` (split kernel-set
digest checked against `run_manifest.json`), `figures.window_label` name. Tests: `tests/test_fine_tuning_gpu.py`
(7): first intervals, atoms and atom responses equal the CPU's `prepare_support`/`pressure_matrices`; two-GPU split
bitwise; `--method cpu` equals the CPU pipeline (`group_basis`, `fit_corrections`) for `all` and `mixed`; confine
total; gauge independence.

**User decisions, in order (2026-10-20/21):** multiplier f > 0 (f < 1 lowers) plus an additive floor on disjoint
groups; support by significance, confine total ε = 1e-3; then NG and RG fitted independently with ε = 1e-2; then
both pressures in one fit (no ratio term; the objective never had one), ε = 1e-3; finally **"the same as the CPU,
fit on top of 1M NC31"**: `--method cpu` is the default; the multiplier variant stays as `--method multiplier`.

**Results, SHAP5 1M, last 5 min, 5 intervals:**
* Strict support (CPU definition): **0 facets** first seen in the window, for N30 and NC31 (it depends on geometry
  only): 840,458 facets first seen before the window, 353,792 never seen. The unchanged CPU tool gives the same
  counts (`finetune_cpu_001m_N30_m5_n5_prepare`, 75 min; GPU 2.5 min). So the CPU method returns its unchanged
  control result (`finetune_gpu_001m_NC31_m5_n5_cpu`). Last 20 min: 35 facets, no effect.
* Every facet seen in the last 5 min (85,631) was seen earlier from higher up; one A_scale added to one of them
  changes some earlier pressure by a median 1e-4 (max 1.3e-2).
* `--method multiplier --confine total`: ε = 1e-3, 146 facets (5,200 m², interval 0, mixed), late log-rms NG 0.6039
  → 0.6038, RG 0.5111 → 0.5112 (unconstrained optimum would change earlier pressures by 224%); ε = 1e-2, 840
  facets (29,400 m², patches at lon 55–75°, lat −25–0°, far from the impact), NG-only fit 0.6039 → 0.6005.
* Conclusion offered: the late NG/RG mismatch cannot be fitted by activity confined to surfaces seen only late;
  the surfaces that matter late are also seen earlier.

## C09 on the SHAP5 3M mesh (2026-10-20, user request): COMPLETE; all products and fits on the GPU

The user asked to repeat the full procedure on `cg_spc_shap5_003m_cart.bdf`.

**Mesh:** 3,145,728 facets, 1,572,866 vertices.

**Base:** `output_fit_shap5_003m_base`, built with `liouville_cops_cuda` in 12.4 min.
856,881,542 kernel records, 33 GB; all 1,150 points checked.

**Visibility templates** (`cuda_c05/restart_003m/visibility_templates_shap5_003m.npz`):

| Template | Facets | Area fraction | 1M area fraction |
|---|---|---|---|
| NG-only | 1,130 | 4.90e-4 | 4.91e-4 |
| RG-only | 3,791 | 1.84e-3 | 1.83e-3 |

**Impact point:** 10:39:26.4 UTC, lon 12.227°, lat 31.749°, as on 1M.

The mesh size forced two engineering changes. Both are documented in the code and tested:
* **Gravity on two GPUs.** `c04_gravity_equivalence.py gpu --devices 0,1` splits the targets
  into contiguous blocks, one per GPU, in spawned processes.
  * The split field is bitwise the one-device field (`tests/test_gravity_split.py`): each
    target is a fixed-order sum over all faces.
  * The manifest records the partition.
  * Expected about 7 h, against about 14 h on one GPU.
* **Memory headroom as a run setting.**
  * Every 3M fit needs both GPUs (`--devices 0,1`); one GPU would need about 50 GiB.
  * Even then, device 0 needed 23.40 GiB of evaluator plan on top of 14.80 GiB of kernel
    data. Under the 20% headroom policy only 22.63 GiB were usable, so the first B29 was
    refused (`fit_003m_gpu_B29_n2000.capacity_refused.log`). That cost 22 min of the
    fitter's kernel-file check and no fit time.
  * New `gpu_fit.py --headroom-fraction` (range 0.05–0.5; default 0.20 unchanged). It is
    carried to the workers by `$COPS_CUDA_HEADROOM_FRACTION` and recorded in
    `compute_backend.memory_headroom_fraction`.
  * The 3M fits use 0.15: about 8.7 GiB of each device stays actually free.
  * Results do not depend on the headroom: it only changes what the budget admits, not the
    arrays or the reduction order.
* **Start-up time.** The unchanged fitter's start-up check reads all 33 GB of kernel files
  at about 30 MB/s, about 22 min per fit. This is fixed overhead.

**Chain** (`scratchpad/chain_003m.sh`, 11:53–22:02 CDT): B29 and N30 on both GPUs, then
gravity on both GPUs (6.96 h), then cache and full bundle (with ground-track and template
figures), then C30 and NC31. No step failed after the headroom change.

**Products:** gravity `cuda_c04/equivalence_20261020T182655Z` (field, cache; offset halving
changes the slope by at most 0.0059°); bundle `cuda_c05/restart_003m/restart_templates_shap5_003m.npz`
(its `late_ng_not_rg` column is bitwise the visibility bundle's); fits
`fit_003m_gpu_{B29,C30,N30,NC31}_n2000` with restart files and figures; package
`comparison_003m_vs_001m_788k/` (README, figures, tables, translated states, inputs,
expected hashes, code, PROVENANCE, MANIFEST of 1,497 files). 161 tests pass.

| Fit | Cost 788k / 1M / 3M | Termination 3M | Parameter change vs 1M / 788k | Notes (3M) |
|---|---|---|---|---|
| B29 | 0.067692 / 0.068711 / **0.067207** | max_nfev | ≤ 0.56 σ / ≤ 0.23 σ, same bounds | T 30/400/43.8, κ 0.838 |
| C30 | 0.067711 / 0.068664 / **0.067404** | max_nfev | ≤ 0.46 σ / ≤ 0.13 σ | cliff60 0.49 ± 6.75 (not significant) |
| N30 | 0.042221 / 0.041887 / **0.041294** | max_nfev | ≤ 0.48 σ / ≤ 0.96 σ (κ) | late_ng 7.51 ± 5.14; T75 206.3 K; κ 0.086 |
| NC31 | 0.067662 / 0.068664 / **0.067368** | max_nfev | ≤ 0.46 σ / ≤ 0.11 σ (late_ng inert) | B29 branch again |

**Mesh convergence:**
* Unfitted bases, 3M vs 1M: median 0.04–0.07%, at most 0.19% above 10 m (788k vs 1M was
  0.07%); below 10 m at most 0.85%. 3M vs 788k: at most 0.26% above 10 m.
* 3M fitted states re-evaluated on 1M and 788k (fit selection): B29, C30 and NC31 costs change
  by −0.14% and −0.21%; predictions by at most 1.0% (at 0.4 m) and 0.7%. The 3M B29 state
  costs 0.0671 on 1M, below the 1M fit's 0.0687: the 1–2% cost differences are optimizer
  stopping points.
* N30 transfer: +1.9% (1M) and +3.4% (788k); prediction changes up to 19%/37% (NG/RG) at
  3.9–11 m, at most 1.1–1.9% above 310 m: the per-mesh NG-only region (1,130 / 500 / 379
  facets, area fraction 4.9e-4 on all).
* Cliff area fractions grow with resolution: ≥ 60° 3.77% / 3.92% / 4.10%, ≥ 45° 13.6% /
  13.8% / 14.3% (788k / 1M / 3M). The cliff amplitude is consistent with zero on all meshes.

**Conclusion offered to the user:** the fitted branches, active bounds and template
amplitudes are the same on 788k, 1M and 3M; the forward model changes by ≤ 0.3% above 10 m.

**Provenance note:** `LiouvilleCOPS/tools/fit_cops.py` was changed in the CPU tree at
15:51 CDT (not by this work): optional `--restart-from` and a CPU restart writer. The 3M C30
and NC31 fits record the new sha256 (`eb82bbbd…`), all other fits `5bd6cf12…`. Without
`--restart-from` the optimizer path is unchanged, and the CPU writer defers to `gpu_fit.py`'s
restart file. The package keeps both CPU code snapshots.

## Publication figure format (2026-10-20, user request; default for all future runs)

The user asked for:
* publication-quality plots (.eps and .png);
* the time and altitude range of the late visibility marked;
* the end-of-mission ground track of Rosetta with the touchdown point.

The user's decisions:
* the Rosetta impact only, not Philae's sites;
* both 3D views and a lon/lat map;
* this format for future runs.

**Code:**
* `cops_cuda/figures.py`: style, pressure figures, track and impact geometry, renderers;
* `cops_cuda/kernels/render.cu`: `raster_faces` z-buffer, `ray_first_hit`;
* `tools/pub_figures.py` (`pressures`, `groundtrack`);
* `tests/test_figures.py` (6 tests). The full suite passes: 157 tests.

**Automatic use:**
* `gpu_fit.py` writes `<fit>/figures/pressures_*`;
* `c05_write_bundle.py` writes `<bundle>_groundtrack_*`;
* `restart_tool.py plot` and `enhance` use the format, and `enhance` adds a region lon/lat
  map.

A figure failure never fails a fit; it is recorded in `compute_backend.figures`.

**Definitions:**
* **Late-visibility window:** the bundle's `late_ng_not_rg` interval, 10:33:10–10:39:26 UTC.
  * 1M: 95 points, 307.0–0.38 m altitude.
  * 788k: 307.0–0.34 m.
* **Ground track:** the sub-spacecraft intercept, i.e. the first surface hit toward the
  body-frame origin.
  * The nearest-point track (the altitude definition) jumps between the bumps of the outer
    envelope at high altitude, so it is not plotted.
  * It stays in `_track.dat`. On 1M it differs from the intercept by a median of 55 m (at
    most 100 m) in the late window, and by 0.14 m at the last row.
* **Impact point:** the straight-line extrapolation of the last two rows to the mesh.
  * 1M: 10:39:26.4 UTC, lon 12.227°, lat 31.749°.
  * 788k: 10:39:26.3 UTC, lon 12.228°, lat 31.750°.
  * The last row is 0.38 m above the 1M surface at 0.89 m/s. The nearest surface point
    of the last row is 0.05 m from the impact point.

**Template maps** (follow-up request: publication form of `physical_templates_shap5_788k_corrected.png`,
with "late visibility" defined and the ground track with the touchdown point):
* `figures.templates_figure` and `pub_figures.py templates`;
* rows: late templates, effective gravitational slope, cliffs; global lon/lat map plus a
  zoom on the impact area; the track and impact point on every panel;
* the definition is in the figure header, and the full definitions are in `_caption.txt`;
* `c05_write_bundle.py` writes them next to every new bundle;
* written for 1M and 788k (`figures/templates_{001m,788k}`); 158 tests pass.

**Regenerated:**
* `comparison_001m_vs_788k/figures/`: per-model 788k vs 1M, all 1M models, and ground track
  on both meshes. The old `plot_cops_model` PNGs were removed.
* `fit_001m_gpu_*_n2000/figures/`.
* The package README, PROVENANCE and MANIFEST were updated, and the code snapshot was
  refreshed.

## C09 on the SHAP5 1M mesh (2026-10-20, user request): all products and fits on the GPU

The user asked whether 788k facets are sufficient and requested a full refit on
`cg_spc_shap5_001m_cart.bdf`.

**Mesh:** 1,194,250 facets, 597,127 vertices; audit passed, Euler characteristic 2,
volume 1.876e10 m³.

**Commands:** `quick_guide.md` §4e.

**Products:**
* fitting base `output_fit_shap5_001m_base`: C11, 5.5 min, 356,553,303 records;
* gravity field `cuda_c04/equivalence_20261020T082211Z`: 2.0 h;
* cache `admitted_gravity_cache_shap5_001m.npz` in the same directory;
* bundle `cuda_c05/restart_001m/restart_templates_shap5_001m.npz`;
* fits `fit_001m_gpu_{B29,C30,N30,NC31}_n2000`;
* comparison `comparison_001m_vs_788k/` (README, plots, table, translated states).

**Not CPU-validated on this mesh:** no CPU gravity or C++ base exists for 1M. The GPU codes
were validated on 100k and 788k.

| Fit | Cost 788k → 1M | Termination 1M | Change of parameters | Notes |
|---|---|---|---|---|
| B29 | 0.067692 → 0.068711 | max_nfev | ≤ 0.34 σ, same bounds | T on bounds 30/400, κ 0.855 |
| C30 | 0.067711 → 0.068664 | max_nfev | ≤ 0.34 σ | cliff60 0.22 ± 5.86 (not significant) |
| N30 | 0.042221 → **0.041887** | **converged (ftol), nfev 1439** | ≤ 0.49 σ | late_ng 6.99 ± 4.83; T75 222.8 K; κ 0.096 |
| NC31 | 0.067662 → 0.068664 | max_nfev | ≤ 0.36 σ (late_ng inert) | ends on the B29 branch again |

**Mesh convergence:**
* The unfitted bases differ by at most 0.07% above 10 m altitude (median 0.01–0.02%); at the
  3 points below 10 m they differ by at most 0.52%.
* A 1M fitted state re-evaluated on 788k (`restart_tool.py translate --evaluate`, fit
  selection) changes the B29, C30 and NC31 data costs by 0.05%. The predictions change by at
  most 0.45%.
* The 1–2% cost differences between the meshes are therefore optimizer stopping points
  (non-converged at 2000 evaluations), not resolution.
* N30's larger transfer change (NG 8.5% and RG 20% at the last points) comes from the
  per-mesh NG-only region definition: 500 against 379 facets, same area fraction.

**Conclusion offered to the user:** 788k resolves the smooth forward model. The finer mesh
does not change any fitted branch.

## C12: fit restart files and tools. Status: IMPLEMENTED (2026-10-20, user request)

* **Code:**
  * `cops_cuda/restart.py` (format, name mapping, GPU state evaluation);
  * `tools/gpu_fit.py` (`fit_restart.npz` written at the end of every fit; `--restart-from`;
    `--no-save-restart`);
  * `tools/restart_tool.py` (`show`, `plot`, `enhance`, `translate`);
  * `masks.visibility_masks` (region masks for any set of points);
  * `tests/test_restart.py`. `pytest`: 150 passed.
* **Evidence on the 788k base** (`cops_output/c12_trials/`):
  * **Save:** a 5-evaluation N30 fit wrote `fit_restart.npz`. The data cost recomputed from its
    stored predictions (0.5734177360) matches the optimizer's (0.5734177330) to print precision.
  * **Restart:** `--restart-from` with no template arguments.
    * The template was supplied from the restart's bundle (SHA-256 checked).
    * All 30 parameters were taken from the restart, none clipped.
    * TRF's initial cost (5.7342e-01) equals the restart's final cost: the fit continued
      exactly from the saved state, reaching 0.1422 after 5 more evaluations.
    * The restart history shows two fit steps.
  * **GPU state evaluation:** for the restart state, `StateEvaluator` gives cost
    0.5734177329781822 against the optimizer's 0.5734177329781821.
  * **Enhance** (`ng_not_rg` below 100 m, `add 2`): 29 points; 379 facets, 92.3% of their
    area clipped in this state. The plots, the region map and the saved edited restart with its
    CPU-format bundle were all written.
  * **Translate** of that edited restart to the 100k base:
    * nearest-centroid mapping, distance median 10.8 m, p95 121 m, max 387 m;
    * amplitudes rescaled by A_scale ratio 0.989;
    * templates renamed `xl_*`;
    * predictions evaluated on the 100k base; history fit → enhance → translate.
  * **Plot:** two restarts compared (5 against 5+5 evaluations).
* **Adjustments of 2026-10-20 (user request):**
  * `factor` accepts any V > 0; V < 1 reduces the activity. The amplitude V − 1 is stored, and
    a refit clips a negative amplitude to 0 (reported).
  * **Effect confined to the window.**
    * Strict exclusion of facets seen from outside the window empties the region for every
      low-altitude window on 788k: ≤ 1000 m keeps 19 facets; ≤ 300, 100, 30 and 10 m keep 0.
    * The user chose a significance threshold, ε = 1e-3 by default (`--epsilon`).
    * A per-facet bound (`--confine facet`) still let the region change outside pressures by
      15.7% (NG) on a 24k-facet region.
    * The default `--confine total` therefore bounds the summed change at every outside point.
      It is exact because the pressures are linear in each facet's clipped activity
      (`kernels/restart.cu`, `tests/test_restart.py`). The prediction agrees with the full GPU
      evaluation to 5e-17.
  * **Finding:** with total confinement the change inside the window is no larger than outside
    (788k, ≤ 100 m windows, `c12_trials/v4_*`). The facets that matter for the low points are
    also seen from the higher points, so a static activity map cannot confine a pressure change
    to the late descent.
* **Notes:**
  * Edited and translated restarts keep the source optimizer record as `source_optimizer`,
    with a `state_note`.
  * A translated restart uses the target base table's rows and selection flags.
  * Mapped cliff masks carry no gravity provenance; recompute them on the target mesh
    (`--target-bundle`).

## C11: GPU rebuild of the fitting base. Status: COMPLETE (2026-10-20); validation policy approved by the user

Implemented ahead of C10 at the user's request.

* **Implementation:** `native/` (`gpu_liouville.cu`, `access_shim.cpp`, `Makefile`,
  `README.md`), building `native/liouville_cops_cuda`.
  * It is the unchanged `../LiouvilleCOPS/src` with GPU definitions of
    `cLiouville::EvaluateKernels` (transport kernels) and `cSourceModel::Evaluate` (cast
    shadows). The CPU definitions are renamed by `-D` in their own translation units.
  * SPICE, the shield Monte Carlo, the source model, `ApplySource`, every writer, the output
    audit and the manifest are the original code. The manifest names the CUDA binary as the
    executable.
  * The BVH and the shield table are read from main()'s objects and cross-checked at
    start-up: GPU `FirstHit` against the CPU on 4,096 rays gives 0 face mismatches, and
    `NGWeight` agrees within ≤ 1e-12.
  * Usage: the documented `rebuild_cops_boundary.py … --executable
    native/liouville_cops_cuda` (device: `COPS_CUDA_DEVICE`).
* **Not bitwise by construction:** the C++ binary uses GCC FMA contraction (40 FMAs in
  `EvaluateKernels`, 96 in the BVH routines); the GPU uses FP64 with `--fmad=false`.
* **Evidence:**
  * Self-test (`Mode=selftest`): ALL PASSED, every result line identical to the C++ run.
    Covered by `tests/test_native_liouville.py`.
  * **788k** (`cuda_c11/compare_20261020T063859Z`):
    * `elements.bin` bitwise;
    * 245,515,533 records in both, 1 + 1 flipped facets (2 points, share ≤ 4.5e-13);
    * CosSZA bitwise, 0 shadow mismatches;
    * point sums ≤ 1.1e-14, worst record 9.7e-9;
    * header geometry bitwise except the altitude (150 points, ≤ 1.0e-15);
    * `COPS_model.dat` pressures identical at print precision.
  * **100k** (`SHAP5_Volume.bdf`; `cuda_c11/compare100k_20261020T065426Z`):
    * `elements.bin` bitwise;
    * 31,672,361 against 31,672,363 records (0 + 2 flipped facets, share ≤ 1.7e-12);
    * CosSZA bitwise, 0 shadow mismatches;
    * point sums ≤ 1.2e-14, worst record 1.1e-9;
    * `COPS_model.dat` ≤ 6.6e-11, i.e. the last printed digit.
  * **Fitter:** the unchanged fitter accepts the GPU-built 788k base, with consistency
    figures identical to the C++ base.
  * **A full B29 fit on the GPU-built base is bitwise identical** to the fit on the C++-built
    base (`fit_788k_gpu_B29_n2000_on_gpu_base` vs `fit_788k_gpu_B29_n2000_fused`).
* **Timing,** one A6000 with 16 host threads, including the post-run C++ audit:

  | Base | C++ | GPU |
  |---|---|---|
  | 788k | 292.5 s | **227.5 s** |
  | 100k | 50 s | 38 s |

  * GPU work about 40 s on 788k: rays about 17 s, shadows about 25 s.
  * The rest is unchanged C++ host work that both builds share: the read-back audit
    (`std::set` of all records, SHA-256 of 9.8 GB), `ApplySource`, per-record `fwrite`
    and setup.
  * The remaining gain is therefore limited by code that stays unchanged (user
    constraint).
* **Validation policy approved** by the user on 2026-10-20 (ROADMAP C11 table, VALIDATION.md).
  Both meshes pass with margins of 100× or more.

## C08: FP64 evaluator restructure. Status: COMPLETE (2026-10-20); mixed precision skipped (user)

* **Implemented:**
  * `kernels/pressure.cu`: `eval_fused_p5`, `eval_fused_p14` and `eval_fused_c1`…`c16`
    (template `fused_body<NP, G>`). Fixed-size, fully unrolled accumulators, held in
    registers: ptxas reports no stack frame and no spills; 37–128 registers.
  * The per-record factors and loads are shared by every sum of a launch.
  * Every per-term expression, the per-thread strided fold, the shared-memory tree and the
    partial layout are those of `eval_partial`.
  * `evaluator.CudaBackend(fused=True)` is the default in tree mode. `_fused_plan` maps the
    standard requests onto pressures(+temperatures) plus column groups of at most 16; any
    other code list uses the generic kernel. `identity()["fused_kernels"]` records which
    path ran.
* **Gate: PASS.**
  * Tests: `tests/test_evaluator.py` checks fused against generic **bitwise** for 6
    configurations × 4 request types × 2 block sizes, plus the plan logic. `pytest`: 145
    passed.
  * `tools/c08_fused_check.py` (`cuda_c08/fused_check_20261020T053223Z`): **bitwise** on
    the full 788k set for B29/N30/NC31, all fixture states, 4 request types.
  * `c06_full_run.py` (`cuda_c06/full_20261020T054054Z`): ALL PASS, with the same error
    figures as before.
  * `c07_multi_gpu_check.py` (`cuda_c07/check_20261020T055749Z`): ALL PASS; one vs two GPUs
    bitwise.
  * **Full B29 fit** (`fit_788k_gpu_B29_n2000_fused`): **bitwise identical** to
    `fit_788k_gpu_B29_n2000_r2` (parameters, final cost 6.769214e-02).
* **Speed** (one A6000):

  | | Before | After | Speedup |
  |---|---|---|---|
  | Full Jacobian | 0.684–0.744 s | **0.261–0.279 s** | **2.6–2.7×** |
  | Pressures | 69–77 ms | 54–59 ms | about 1.3× (host matvec 22 ms now dominates) |
  | B29 fit, 2000 evaluations | 1,627 s | **985 s** | **1.65×** |

  * In the fused fit, the Jacobians took 413 s (was 1,035 s) and the pressures 113 s; about
    390 s is the fixed per-fit overhead.
  * Two GPUs, full Jacobian: 0.151–0.161 s, i.e. 1.7× over one GPU.
  * Benchmark: independent per-GPU fits 7.50 Jacobians/s against one split fit 6.60
    Jacobians/s.
* **Remaining costs:** the fixed per-fit overhead (the CPU fitter's base audit, TRF and the
  writers; about 6.5 min of a 16-min fit) and the host activity matvec (22 ms per
  evaluation). Neither depends on GPU precision.

## C08 decision support: evaluator profile (2026-10-20)

`tools/profile_evaluator.py`; report `cops_output/cuda_c08/profile_20261020T052129Z/report.json`.
Full 788k B29 problem, one A6000, tree mode. Hardware counters are unavailable without
admin rights (`RmProfilingAdminOnly=1`), so costs are separated by timing and controlled
experiments.

* **Full Jacobian 0.688 s:**
  * `eval_partial` 0.657 s (96%), in 4 launches;
  * host activity matvec 0.029 s;
  * transfers and assembly about 0.03 s in total (they overlap the matvec).
* **Pressures 0.069 s:** kernel 0.043 s, host matvec 0.022 s.
* **Per sum:** each sum adds 7.4–10.4 ms (0.030–0.043 ns per record and sum) for every code
  family; the per-launch base cost is 7–16 ms.
* **Memory floor:** reading every per-record array once, plus the facet gather, takes
  21 ms (712 GB/s, close to the 768 GB/s peak).
* **FP64 throughput,** measured as independent mul+add chains without FMA: 267 GFLOP/s.
  A pressure term (about 3 FP64 operations per record) would therefore cost about 2.8 ms
  per sum if bound by FP64 arithmetic. The measured cost is 3× that.
* **Root cause (ptxas `-v`, sm_86):** `eval_partial` has a **256-byte stack frame**: the
  runtime-indexed `double acc[MAX_SUMS]`. Every term therefore costs a local-memory load and
  store, plus a load and switch on `codes[k]`. The accumulators are not held in registers.
* **Whole-fit split** (B29, 2000 evaluations, 1,627 s): Jacobians 1,035 s, pressures 138 s,
  bind 69 s, rest 385 s. The rest is the CPU fitter's base audit (hashing 9.8 GB of
  kernels), TRF and the output writers.
* **Conclusion for C08:**
  * An FP64-preserving restructure targets the dominant cost: register accumulators
    through compile-time-specialized kernels for the fixed code sets, and one fused pass
    sharing `u` and the loads.
  * It keeps every per-term expression and the per-thread fold order, so its results stay
    bitwise equal to today's tree mode.
  * Estimate: Jacobian 0.66 s → 0.2–0.25 s; fit about 27 → 13–15 min.
  * Mixed precision would add only about 1.3× on top (fit about 10–11 min), because the
    fixed per-fit overhead of about 6.5 min and the memory floor then dominate. It would
    cost the approved budgets and optimizer reproducibility.

## C09 first production fits (2026-10-20): max_nfev 2000, template start 1.0, rebuilt base and bundle

Outputs: `cops_output/fit_788k_gpu_{B29,C30,N30,NC31}_n2000_r2`. Commands: `quick_guide.md` §3a.
All four stopped at the 2000-evaluation limit (not converged). The template start was
applied (`compute_backend.fit_start.applied`). One GPU per fit, two fits at a time; 27–34 min each.

| Fit | Final cost | rms NG | rms RG | Template amplitude(s) ± σ | Tmax / T75 / Tmin [K] | κ | Data Jacobian rank, condition number |
|---|---|---|---|---|---|---|---|
| B29 | 0.067692 | 0.3217 | 0.1785 | — | 30 (bound) / 400 (bound) / 44.2 | 0.845 | 29/29, 1.2e5 |
| C30 | 0.067711 | 0.3219 | 0.1783 | cliff60: 0.38 ± 8.97 | 30 (bound) / 400 (bound) / 44.2 | 0.845 | 30/30, 1.3e5 |
| N30 | **0.042221** | **0.2238** | 0.1854 | late_ng: **7.21 ± 5.46** | 30 (bound) / **227.4** / 69.7 | **0.105** | 30/30, 9.4e5 |
| NC31 | 0.067662 | 0.3219 | 0.1781 | late_ng: 1.14 (σ 3.9e-13, inert); cliff60: 0.39 ± 8.95 | 30 (bound) / 400 (bound) / 44.3 | 0.845 | 30/31, 5.7e20 |

Observations, for the user's scientific review (no conclusions drawn here):
* B29 reproduces the lost run's final cost exactly (deterministic).
* **N30 reaches a different branch.** Its cost is 38% below B29: NG rms 0.224 against 0.322,
  RG slightly worse. T75 leaves the upper bound, Tmin rises to 70 K and κ falls to 0.105.
* **NC31**, which contains N30 as the special case c_cliff = 0, ends on the B29 branch. Its
  NG-only amplitude is inert again: zero gradient, σ ≈ 4e-13, rank deficiency. Which branch
  a fit reaches therefore depends on the optimizer path. All fits are non-converged at 2000
  evaluations.
* The cliff ≥ 60° amplitude is not significant in C30 or NC31.

## Data loss on 2026-10-20 and rebuild (complete)

* **Incident:** between about 20:41 and 21:00 local, the following were removed from
  `cops_output` by an unidentified actor: every `output_*` run directory, including the 788k
  fitting base, and every `cuda_c02`…`cuda_c07` directory (gravity caches, restart bundle,
  evidence reports). The `fit_*` and `physical_templates_*` products, `cuda_reference` and
  `cuda_doctor` survived. No backup exists (user).
* **Effect:**
  * the GPU B29 production fit (2000 evaluations, final cost 6.7692e-02) could not write
    its outputs, so only its log remains;
  * C30 failed at start;
  * the N30 and NC31 queue was stopped (user).
* **Base rebuilt** with the documented C++ command (`rebuild_cops_boundary.py --input
  fit_base.in --surface-bc flat --shape-model …788k_cart.bdf --threads 16`): 293 s.
  * **Byte-identical** to the lost base, checked against the hashes recorded in the C00
    capture and the fit manifests: `elements.bin`, `COPS_model.dat`, `input.snapshot.in`
    and **all 1,150 kernel files**.
  * Only `run_manifest.json` and `rebuild.json` differ (timestamps).
* **GPU products and evidence regenerated** (2026-10-20, two GPU queues). Every check passes
  again with the same figures, and the evidence paths in this file now point to the
  regenerated reports:
  * C02 load: ALL PASS.
  * C03 gravity: ALL PASS; the full 788k field now lives in the C04 equivalence directory.
  * C04 788k: compare PASS (relative error ≤ 2.45e-13, slopes ≤ 2.49e-11°, masks
    identical); admitted cache written.
  * C04 100k: compare PASS (≤ 7.8e-14); admitted cache written.
  * C05: masks ALL PASS; the restart bundle `cuda_c05/restart_rebuilt/` is bitwise equal to
    the CPU builder in `late_ng_not_rg`, `cliff_ge_45deg` and `cliff_ge_60deg`.
  * C06: check, full run and Jacobian budget all pass (244.7 eps·M_J).
  * C07: ALL PASS (1.90× Jacobian speedup).
  * Not repeated: the CPU builder's reuse test of a GPU cache (C04).
* **Production fits restarted** with suffix `_n2000_r2`. The earlier `fit_788k_gpu_B29_n2000`
  (outputs not written) and `fit_788k_gpu_N30_n2000` (stopped) are left as they are,
  marked failed or incomplete.

## User decisions of 2026-10-20 and the quick guide

* The temperature bounds stay 30–400 K; temperatures on the bounds are reported as
  boundary-limited.
* Template amplitudes start at c = 1 (bound 0): `gpu_fit.py --template-start`, default 1.0.
  `tests/test_gpu_fit.py` checks that only the optimizer start changes and that the prior
  `P.x0`, the regularization and D stay the fitter's.
* `max_nfev` is a driver parameter: 2000 by default, and 1000 with `--cpu-comparison`
  (which also sets template start 0, i.e. the CPU command exactly).
* `quick_guide.md`: exact commands for setup, the one-time products, production and
  comparison fits, fit comparison (`tools/compare_fits.py`), C++ replay, validation and plots.

## C07: two-GPU fitting. Status: COMPLETE (2026-10-20)

* **Implemented:**
  * `cops_cuda/multi_gpu.py`: `MultiGpuBackend` behind the C02 contract.
    * One optimizer on the host and persistent spawned workers, one per device.
    * Each worker owns contiguous rows balanced by record count and rebuilds its model with
      the reference parser (element cache checked by digest).
    * One host activity field per evaluation, in shared memory.
    * Assembly in row order.
    * Failures and dead workers become `BackendError`; `close()` is graceful, then
      terminates.
  * `gpu_model`: `REGISTRY.configure(devices=[0, 1])` and `seen_faces` as the union over
    workers.
  * `tools/gpu_fit.py --devices 0,1`.
  * `CudaBackend.activity_source` hook.
  * `tests/test_multi_gpu.py`:
    * partition;
    * two devices bitwise equal to one in both modes;
    * a SciPy TRF trajectory bitwise on two devices;
    * `seen_faces` union;
    * dead worker gives a clean error.
    * `pytest`: 140 passed.
* **Gate** (`tools/c07_multi_gpu_check.py`; `cops_output/cuda_c07/check_20261020T031417Z`,
  ALL PASS), 788k base with all 1150 points:
  * **one vs two devices bitwise identical** (pressures and every Jacobian column) for
    B29, N30 and NC31, 13 fixture states;
  * repeat evaluations bitwise (deterministic pinned reductions);
  * partition 122,962,701 / 122,552,832 records;
  * memory bounded per device: each worker owns its `DeviceBudget`, about 6.6 GiB per
    device for one split fit.
* **Latency** (warm medians):

  | | One GPU | Two GPUs | Speedup |
  |---|---|---|---|
  | Full Jacobian | 0.683–0.738 s | 0.359–0.388 s | 1.89–1.90× |
  | Pressures | 68–72 ms | 53–58 ms | small (host matvec and IPC dominate) |
  | Bind | 68 s | 37 s | |

* **Benchmark: one split fit vs independent fits (B29):**
  * two concurrent independent fits, one per GPU: 0.682 s per Jacobian each, **2.93
    Jacobians/s** in total, with no contention between devices;
  * one fit split over both GPUs: **2.79 Jacobians/s**, at half the latency.
  * **Recommendation:** run a fit matrix (C09) as independent single-GPU fits, about 5%
    more throughput and no IPC. Use the split mode when one fit's turnaround matters.
* **Full two-GPU fit:** `gpu_fit.py --devices 0,1`, B29 (`fit_shap5_788k_gpu2_B`).
  **Bitwise identical** to the one-GPU fit `fit_shap5_788k_gpu_B`: every parameter and every
  predicted NG/RG pressure (`tools/compare_fits.py`). Wall time 730 s against 1,049 s.
  Jacobians 276 s against 523 s, bind 37 s against 69 s.
* CPU threads: the optimizer is serial and is not parallelized
  through thread settings (ROADMAP C07 item 3); BLAS threads stay at 1 as in the campaign.

## C07 step 1: GPU fit with the unchanged fitter. Status: first B29 fit COMPLETE (2026-10-20); CPU-vs-GPU end-to-end comparison waits for the CPU B fit

* **Implemented:**
  * `cops_cuda/gpu_model.py`: a GPU-backed subclass of the reference `fit_cops.Model`.
    * Host side is header-only.
    * `evaluate` has the reference signature and return format, including the T and κ
      columns when `jac=True`.
    * `seen_faces` runs on the device (`masks.cu: seen_positive`, the reference `> 0` predicate).
    * Backends are bound lazily and shared across models with equal rows; one is rebound
      only for a new enhancement basis.
  * `tools/gpu_fit.py`: runs the unchanged `fit_cops.main()` with `fit_cops.Model` replaced
    in-process; no file in `../LiouvilleCOPS` changes.
    * Writes a `compute_backend` section into the fit manifest on every save.
    * Reproduces the fitter's failure marking.
    * Refuses L-curves (they need fork).
  * `tests/test_gpu_model.py`:
    * sequential `evaluate` and `seen_faces` bitwise equal to the reference;
    * **a full SciPy TRF trajectory (`fit_cops.fit_one`) bitwise identical** to the CPU;
    * tree mode within budget;
    * backend sharing and rebinding. `pytest`: 135 passed.
* **First GPU fit, B29:** the campaign's B command (profile `all_available`, degree 4, λ = 0,
  tbounds 30–400, `max_nfev` 1000, tolerances 1e-10). Output in
  `cops_output/fit_shap5_788k_gpu_B`.
  * The fitter's own checks pass:
    * consistency with LiouvilleCOPS ≤ 4.6e-11 (NG) and 9.9e-10 (RG) against the 1e-6 tolerance;
    * initial cost 3.830468 and seen area 70.5%, identical to the CPU B fit's log.
  * Result:
    * stopped at `max_nfev = 1000` (not converged); final cost 6.7814e-02;
    * rms ln residual NG 0.3218, RG 0.1790; optimality 4.6e-4;
    * data Jacobian rank 29 of 29; condition number 1.3e5.
  * **Temperatures at bounds:** Tmax = 30 K (lower) and T75 = 400 K (upper); Tmin = 44.3 K;
    κ = 0.847. The archived B29 (930/1146 points) shows the same pattern. This is the
    ROADMAP §3.4 trigger, a physical decision for the user.
  * Timing: **1,049 s wall** on one A6000.
    * 768 Jacobians took 523 s, 1,004 pressure calls 69 s, and the bind 69 s.
    * For comparison, the CPU B fit with the same command had been running for more than 5.5 h.
* **N30 GPU fit:** the campaign's N command, with `late_ng_not_rg` from the same
  `physical_templates_shap5_788k.npz_heal_defect` bundle as the CPU N fit. That column is
  geometry-only and unaffected by HEAL. Output in `cops_output/fit_shap5_788k_gpu_N`.
  * Same startup checks as B29. Wall time 1,146 s, including one rebind for the template
    basis (2 × 68 s).
  * Result:
    * `max_nfev` 1000 (not converged); final cost 6.7824e-02;
    * rms NG 0.3219, RG 0.1790;
    * Tmax = 30 K and T75 = 400 K on their bounds, Tmin 44.3 K, κ 0.847, i.e. B29 within 1e-3.
  * **`c_late_ng_not_rg` = 1.05e-10, at its lower bound 0.** The data Jacobian has rank 29 of
    30 (condition number 2.8e21).
  * **Cause (model structure, not numerics):**
    * At the prior, all 379 template facets have positive legacy linear activity.
    * At the final B and N states, none does: the fitted degree-4 baseline is negative
      on every template facet.
    * The single clip is `max(0, legacy + c·A_scale·T)`. With c ≈ 0 the template facets
      are clipped, so ∂p/∂c = 0 exactly (the reference's `[act_lin > 0]` mask), and the
      amplitude cannot leave its bound.
    * It would become active only once c·A_scale exceeds the negative baseline there.
  * This is a scientific finding for the user. With this parameterization a nonnegative
    template on facets that the free baseline drives below zero is unidentifiable from c = 0.
    It concerns N30, NC31 and C30, whichever the cliff facets turn out to be.
  * Area fraction of each restart-bundle template where the B29 final baseline is positive:
    * `late_ng_not_rg` 0.000;
    * `late_rg_not_ng` 0.000;
    * `cliff_ge_45deg` 0.141;
    * `cliff_ge_60deg` 0.076.
    The cliff amplitudes keep a nonzero gradient, but only through those minorities.
* **Pending:**
  * End-to-end comparison with the CPU B fit (`fit_shap5_788k_restored_B`) once it finishes,
    then a proposal for the end-to-end budget. Because both fits stop at `max_nfev`, the
    comparison is between two 1000-evaluation trajectories. TRF amplifies 1e-13 differences
    along them, so a parameter-space criterion may not be meaningful there.
  * The N30 GPU fit, and the C30/NC31 fits with the restart bundle (C09 scope).

## C06: single-GPU FP64 evaluator. Status: COMPLETE (2026-10-19): numerical gate passed on the C00b subset and the full 1150-point set; Jacobian budget approved

* **Implemented:**
  * `cops_cuda/evaluator.py` (`CudaBackend`, behind the C02 contract);
  * `cops_cuda/kernels/pressure.cu`;
  * `tests/test_evaluator.py`, `tools/c06_evaluator_check.py`.
* **Exactness design:**
  * The static factors of `Model.load_kernels`/`prepare` are computed on the device once,
    in the CPU's operation order (`--fmad=false`).
  * Per evaluation, the host computes `act_lin = Y a + A_scale B c`, the clip and the
    active set with the CPU's own NumPy expression, so the clipping decisions are identical.
  * The device computes the per-contribution terms and per-point sums: 5 pressure sums,
    9 temperature sums, and 2 per activity/template column, in tiles of 32 per launch.
  * The host assembles pressures and columns with `Model.evaluate`'s expressions.
* **Modes:**
  * `sequential` (verification): one thread per point, records in file order, the
    `np.bincount` and CSR order.
  * `tree` (production): fixed record blocks, a fixed shared-memory tree, block partials
    added in order.
* **Supported:** `SZAMode=continuous`, all RG models, both heliocentric models, every
  derivative kind. `srccg_table` is refused explicitly.
* **Evidence:** `cops_output/cuda_c06/check_20261020T023349Z/report.json`, ALL PASS, on
  the 47-point C00b subset (8,783,677 contributions); B29/N30/C30/NC31/R30 with 9–10
  states each.
  * **sequential is bitwise identical to the CPU reference** for every state: pressures
    and all 29–31 Jacobian columns. This includes the archived-B29 state (85% of
    contributing facets clipped, Tmax and T75 on their bounds, 1.34M floor-branch
    contributions).
  * The CPU reference in the CUDA environment reproduces the Python 3.8 C00b fixtures
    bitwise in 47 of 47 cases.
  * **tree:** pressure ≤ 234 eps·M against the 512 budget (the CPU's own summation error is
    ≤ 233 eps·M). The derived residual and cost bounds are met in every state. The
    Jacobian is ≤ 6.2e-14 × column max (provisional budget 1e-11). Repeat runs are bitwise
    identical.
  * Warm tree timing on this subset: pressures 20–29 ms; full Jacobian (29–31 columns)
    45–55 ms. Host matvec and transfers are included.
  * `pytest`: 116 passed. `test_evaluator.py` adds bitwise sequential checks for 6
    configurations, plus tree budgets and drop-in use under `fit_cops.Problem` for B29 and
    NC31.
* **Full-set run** (`tools/c06_full_run.py`;
  `cops_output/cuda_c06/full_20261020T024518Z/report.json`, ALL PASS):
  * Host-light model: `cops_cuda/host_model.py`, a header-only subclass of the reference
    `fit_cops.Model`. The unmodified `fit_cops.Problem` is built on all 1150 fit points
    (245,515,533 contributions) in 0.2 s, with no kernel records on the host.
    `tests/test_host_model.py` proves its per-point data, constants and Problem are bitwise
    the reference's, and that a backend bound from it is bitwise one bound from the full model.
  * Verification uses per-point independence: a 23-point CPU subset (`choose_subset`, 16
    spaced plus late rows; 3.65M contributions) is the exact reference for those points of
    the full problem.
  * B29/N30/NC31/R30, 3–5 fixture states each:
    * GPU sequential on the subset is bitwise the CPU;
    * full tree is bitwise the subset tree at the subset points (point-set invariance);
    * full-tree pressures are ≤ 200 eps·M (budget 512);
    * the full-tree Jacobian is ≤ 2.0e-14 × column max (provisional 1e-11).
  * Timing on the full set (warm, tree, one A6000):
    * pressures 52–63 ms;
    * full Jacobian (29–31 columns) 0.66–0.72 s; `Problem.data_jacobian` through the facade is the same;
    * CPU estimate for the same Jacobian 72–78 s (subset time per record × records), about 105×;
    * bind 68 s (disk streaming plus audit); device plan 13.2 GiB.
* **Jacobian budget: approved by the user (2026-10-19):** ≤ 1024 eps·M_J per entry
  (VALIDATION.md). Measured maximum 244.7 eps·M_J
  (`cuda_c06/jacobian_budget_20261020T025421Z`). The first report (329) used KT instead of
  |KT| in the temperature magnitudes and was corrected on 2026-10-20.
  `tests/test_evaluator.py` now asserts the per-entry budget.
  * A 1-vs-2 GPU comparison belongs to C07.

## C02 joint extensions. Status: COMPLETE (2026-10-19)

* **Identities and serialization:**
  * `backend.problem_contract(P)` identifies the parameter layout, template names and basis
    digest, the enhancement scale and the **original prior** (run coefficients, x0,
    regularization operator), separately from any initialization.
  * `ParameterState.to_dict`/`from_dict` and `ColumnRequest.as_dict`/`from_dict`
    round-trip exactly (`float.hex`).
  * `install()` refuses stale or reordered templates and a changed scale.
* **Tests:** the NC31 synthetic configuration passes bitwise through the facade.
  `test_backend_contract.py` has 22 tests.

## C00b: joint-mode re-freeze. Status: COMPLETE (2026-10-19)

* **Stage:** `joint` of `tools/c00_capture_reference.py`, run under the CPU reference
  interpreter. Output: `cuda_reference/c00_20261019T064536Z/joint/`, with the repaired
  reference tree digest `b602a045…` in its stage manifest. The historical `fit` stage is
  unchanged.
* **Subset:** the same 47-point subset (identical to the C00 subset); kernels checked
  against the frozen digests.
* **Templates:**
  * `late_ng_not_rg`: 379 facets, 4.94e-4 of the area;
  * **`late_rg_not_ng`** (new CPU replica, the C05 test oracle): 1,336 facets, 1.83e-3 of
    the area, zero overlap with NG-only;
  * `syn_cliff_standin`: synthetic, standing in for the cliff mask until the C05 bundle.
* **Families:** B29, N30, C30, NC31 and the labelled R30.
* **States:** prior, mild and strong baseline, temperature floor mix and wide, κ scaled,
  mild and strong amplitudes, a template on its bound, and the archived 29-parameter
  optimum (`output_fit_full_current16`).
* **Stored per state:** pressures, exact sums and magnitudes, residuals and Jacobians
  (λ = 0 and λ = 1e-4 with joint regularization), and branch-aware central differences on
  selected states.

## C05: masks and effective fields. Status: COMPLETE (2026-10-19): visibility masks and the restart bundle

* **Implemented:** `cops_cuda/masks.py`, `cops_cuda/kernels/masks.cu`,
  `tools/c05_mask_check.py`.
* **Masks:** the four same-point predicates of the reference builder (`point_ng`,
  `point_rg`), united over the selected points with integer `atomicOr` on one 32-bit
  word per facet, which is deterministic and race-free:
  * `ng_seen`, `rg_seen`;
  * `late_ng_not_rg`;
  * **`late_rg_not_ng`** (user decision of 2026-10-19: facets only RG sees, same late
    interval).
* **Point selection:** the builder's own (`parse_utc` interval, `inside == 0`).
* **Evidence:** `cops_output/cuda_c05/masks_20261020T023337Z/report.json`, ALL PASS.
  * Late interval (95 points):
    * `late_ng_not_rg` (379 facets) is **bit-identical to the unchanged reference
      `late_gauge_template`**;
    * `late_rg_not_ng` (1,336 facets) is bit-identical to the C00b CPU replica and to a
      NumPy replica;
    * the NG-seen and RG-seen counts (72,411 each) and the point list equal the reference.
  * Control interval (10 points): both exclusive masks are empty, as in the reference.
  * GPU 5.9 s against CPU 1.8 s. Loading and hashing 95 kernel files dominates; masks are
    not a performance target.
* **Restart bundle (C05 item 4):** `cops_cuda/restart_bundle.py`, with
  `tools/c05_write_bundle.py` as the glue and `tests/test_restart_bundle.py`.
  * Product: `cops_output/cuda_c05/restart_rebuilt/restart_templates_shap5_788k.{npz,json,png}`,
    made from the admitted 788k cache (C04).
  * Columns:
    * `late_ng_not_rg` (379 facets);
    * `late_rg_not_ng` (1,336 facets);
    * `cliff_ge_45deg` (110,073);
    * `cliff_ge_60deg` (31,212).
  * 0 facets fell in the 1e-6° threshold band.
  * The unchanged fitter's `load_bundle` and `select_templates` accept all four columns,
    including the gravity provenance gate.
  * **Bitwise identical to the unchanged CPU builder's corrected bundle**
    (`physical_templates_shap5_788k_corrected.npz`) for every column the CPU produces
    (`late_ng_not_rg`, `cliff_ge_45deg`, `cliff_ge_60deg`).
  * Restart use: `fit_cops.py ... --activity-template <bundle.npz> <name>`.

## C04: production gravity cache. Status: single-GPU COMPLETE (2026-10-19): 100k and 788k equivalence PASSED, admitted caches written; two-GPU partitioning pending (C07)

* **User decisions (2026-10-19):**
  * test the gravity model at full-field scale with the repaired CPU builder on
    `SHAP5_Volume.bdf` (100k faces);
  * restart reuses the persisted gravity cache and the mask bundle
    (`late_ng_not_rg`, `late_rg_not_ng`, cliff masks);
  * no workflow change in `../LiouvilleCOPS`.
* **Derived body:** `SHAP5_Volume.bdf` is a Gmsh volume mesh (972,652 CTETRA). Its 100,000
  CTRIA3 surface is a closed, consistently oriented manifold (Euler 2) wound **inward**
  (V = −1.909e10 m³), so the repaired builder refuses it.
  `tools/derive_outward_surface.py` wrote
  `shape_models/derived/SHAP5_Volume_surface_outward.bdf` (with `.provenance.json`):
  * vertices byte-identical through the reference reader;
  * faces = source[:, [0,2,1]];
  * V = +1.909e10 m³;
  * passes the CPU reference audit and the CUDA audit.
* **Audit parity:** `GravityBody` now also runs and stores the CPU's
  `require_outward_closed_body` result. The unchanged CPU builder reuses a cache only when
  its `topology` and `orientation_policy` equal these, so a GPU-written cache can carry
  them truthfully.
* **GPU field:** 787,754 targets × 100,000 faces in 406 s (1.94e8 pairs/s), in
  `cuda_c04/equivalence_20261020T022637Z/`.
* **CPU reference:** the unchanged repaired builder, on CPUs 28–35 beside the user's CPU
  campaign: 5,617 s of gravity.
* **Full-field equivalence: PASS** (`cuda_c04/equivalence_20261020T022637Z/compare_*.json`, regenerated 2026-10-20), on all
  787,754 facets:
  * |Δg| ≤ 7.8e-14 relative (median 1.6e-14) and ≤ 0.003 eps·M, with M/|g| up to 1.5e5;
  * slope ≤ 1.7e-11°;
  * 45°/60° cliff masks identical, with **0** facets in the 1e-6° band and 0 disagreements;
  * the CPU bundle's cliff columns equal the cache slopes;
  * offset-halving 0.0053143694° on both sides (difference 5e-13°).
* **By-product: corrected 100k-vs-788k comparison**, same targets; this is physics, not
  numerics.
  * Slope difference median 0.78°, **p95 8.1°**, max 18.4°.
  * Cliff area at ≥45°: 16.2% (100k) against 13.6% (788k); at ≥60°: 5.0% against 3.8%.
  * Classification differs on 2.7% (45°) and 1.3% (60°) of the area.
  * With correct gravity, the 100k surface still fails the historical 5° p95 resolution
    gate. The exact 788k surface remains the right gravity body for production.
* **Cache writer:** `cops_cuda/gravity_cache_io.py` and `tools/c04_write_gpu_cache.py`.
  * They write a GPU field through the CPU's own `make_definition` and `write_cache`.
  * The calculation metadata names the CUDA backend and its evidence. It records the
    reference audit result verbatim as `topology` and the reference `orientation_policy`,
    `integrity_check = "DISABLE after topology audit"`, and the GPU's offset-halving audit
    on the builder's 512-facet sample.
  * No CPU package identity is copied.
  * `tests/test_gravity_cache_io.py` checks loader compatibility, the builder's reuse
    conditions, refusal for another body, and the offset-audit refusal.
  * **The unchanged CPU builder reuses a GPU-written cache:**
    `cuda_c04/equivalence_20261019T185043Z/gpu_gravity_cache_reuse_test.npz` (100k body;
    offset-halving 0.00531°; removed in the 2026-10-20 data loss, test not repeated) was reused with `--gravity-cache` ("reused gravity cache …
    same reference omega"). Its SPICE ω is bit-identical to the tool's.
  * It wrote `cpu_builder_reuse_of_gpu_cache_test.npz`: 379 NG-only facets; cliff ≥45° on
    16.2% of the area and ≥60° on 5.0%.
  * The **unchanged fitter's `select_templates` accepts its `cliff_ge_60deg` column**, and
    the metadata states the CUDA backend.
  * The restart path therefore needs no change in `../LiouvilleCOPS`.
* **788k production equivalence: PASS**
  (`cuda_c04/equivalence_20261020T021855Z/compare_20261020T031257Z.json`). The CPU reference
  is the corrected 788k build of the restored campaign
  (`physical_templates_shap5_788k_corrected.gravity.npz`; 3.6 h on CPUs). On all 787,754
  facets:
  * |Δg| ≤ 2.45e-13 relative and ≤ 0.001 eps·M;
  * slope ≤ 2.5e-11°;
  * 45°/60° cliff masks identical;
  * the CPU bundle's cliff columns equal the cache slopes.
* **Cross-device determinism:** the GPU 1 field (3,191 s, 1.94e8 pairs/s) is **bitwise
  identical** to the C03 field from GPU 0, for self-gravity, term magnitude and potential.
* **Admitted caches** (CPU format, with the compare report as `evidence`):
  * 100k test body:
    `cuda_c04/equivalence_20261020T022637Z/admitted_gravity_cache_shap5vol100k.npz`
    (offset-halving 0.00531°);
  * **788k production:**
    `cuda_c04/equivalence_20261020T021855Z/admitted_gravity_cache_shap5_788k.npz`
    (offset-halving 0.00502°).
* **Pending:** two-GPU partitioning (C07).

## C03: single-GPU exact polyhedral gravity in FP64. Status: COMPLETE (2026-10-19); the full-mask comparison against the corrected CPU reference passed in C04

* **Implemented:**
  * `cops_cuda/gravity_body.py` (project code): BDF through the reference reader, and the
    exact topology audit that replaces HEAL;
  * `cops_cuda/gravity.py` and `cops_cuda/kernels/gravity.cu` (GPLv3, translated from
    polyhedral-gravity 3.3);
  * `tools/c03_gravity_check.py`, `tests/test_gravity.py`;
  * docs: `GRAVITY_REFERENCE.md` §8, the VALIDATION.md C03 amendment, `RUN.md`.
* **Physical assumptions (unchanged from the CPU reference):** homogeneous density 533 kg/m³;
  targets at centroid + 0.1 m·normal; G = 6.67430e-11; centrifugal term at the centroid;
  C00 SPICE omega. The body's winding is evaluated as given, with no HEAL (see the defect in
  `GRAVITY_REFERENCE.md` §5.1). It is admitted only after the topology audit: closed,
  consistently oriented, Euler 2, V > 0, no zero-area or unreferenced item. The 788k body
  passes.
* **Algorithm items:**
  1. Precompute on the CPU: bitwise equal to the package state.
  2. Operation-order translation, compiled with `--fmad=false`.
  3. Fold (the package's order) and fixed-hierarchy split variants; no target×source matrix.
  4. Acceleration and potential only.
  5. Synthetic tests.
  6. Offset-halving audit and difficult samples.
* **Gate evidence:** `cops_output/cuda_c03/gravity_20261020T022050Z/report.json`.
  * **Synthetic** (C00 fixture, campaign cp38 binary): cube, tetrahedron and icospheres L0–4,
    each original, translated, rotated and double-density. Worst case **0.33 eps·M**,
    against a budget of 16 eps·M. The NaN pattern is identical to the package's on every
    exact-coincidence probe.
  * **Precompute:** G, N and n of all 787,754 faces are bitwise equal to
    `GravityEvaluable.__getstate__`.
  * **Wheel equivalence:** the cp311 package is bitwise identical to the cp38 campaign binary
    on 16 788k targets, so it serves as the CPU reference for new samples.
  * **788k C00 targets** (832 at 0.1 m, 832 at 0.05 m, 32 on-surface; fold and split): every
    set passes ≤ 1e-12 relative (max 2.3e-13) and ≤ 16 eps·M, although M/|g| reaches
    1.1e6. Slopes are within 1e-9° and cliff classifications at 45° and 60° are identical.
    The offset-halving audit gives 0.0050239°, equal to the C00 0.00502388° within 1e-9°.
  * **Difficult samples** (new, 160 points on 32 faces: on a vertex, on an edge midpoint,
    +0.1 m above both, +1 mm above the centroid): max relative error 1.4e-13. The single
    package NaN (an exact vertex) is NaN on the GPU too.
  * **Determinism:** split repeats bitwise; fold with different chunking is bitwise; fold on
    GPU-c314803e… equals GPU 0 bitwise.
  * `pytest tests/test_gravity.py`: 21 passed.
* **Performance (FP64, one A6000):**
  * fold: about 2.0e8 face–target pairs/s, at 162 W of 300 W;
  * 832 targets: split 5 s, fold 38 s (the fold has only 7 blocks for so few targets);
  * kernel resources: 162 registers and 72 B stack per thread;
  * the CPU reference on 8 contended cores takes 0.041 s per target;
  * full field: see the next item.
* **Full 788k field (`--full`, fold, GPU 0):**
  * 787,754 targets × 787,754 faces in **3150 s**, 1.97e8 pairs/s, beside the running L5
    CPU fit. The HEAL-affected CPU build took 13,825 s on 36 threads.
  * Bitwise equal at the C00 targets to the separately validated fold runs; all values
    finite.
  * |g_self| = 1.26e-4 to 2.33e-4 m/s².
  * Slope min/median/max 0.03° / 23.1° / 161.1°. Slopes above 90° are overhangs.
  * **≥ 45°: 110,073 facets, 13.56% of the area; ≥ 60°: 31,212 facets, 3.77% of the area.**
    No facet lies within 1e-6° of either threshold, so no CPU re-resolution is needed. The
    HEAL-corrupted CPU bundle (`physical_templates_shap5_788k`) has 71.0% and 59.5% and a
    72° median slope.
  * Written to `full_field.npz` in the report directory. This is C03 evidence validated on
    samples; it is not a campaign template bundle (C04/C05).
* **Limitations and next steps:**
  * The cliff-mask gate ("reproduces valid slope classifications") is shown on the C00 and
    difficult samples only. The CPU builder has since been repaired (DISABLE after the
    topology audit). A full-mask comparison needs either the repaired full CPU 788k field
    or an accepted import route: the C04 decision. The old bundle is quarantined as
    `*_heal_defect`.
  * `gravity_body` should call the CPU's stricter `require_outward_closed_body` audit
    (connectivity, vertex links, unique triangles); planned in ROADMAP C03/C04.
  * C04 splits targets over both GPUs, adds checkpoints and writes `gravity_cache` products.
  * C10 tunes registers and tiles.

## C02: backend contract, memory layout and provenance. Status: COMPLETE (2026-10-19)

* **Implemented:**
  * `cops_cuda/errors.py`, `formats.py`, `precision.py`, `memory.py`, `backend.py`,
    `cuda_sources.py`, `device_store.py`, `manifest.py`, `kernels/store_audit.cu`;
  * `tools/c02_device_load.py`;
  * tests `tests/test_formats_memory.py`, `tests/test_device_store.py`,
    `tests/test_backend_contract.py`, `tests/synthetic_run.py`, `tests/conftest.py`;
  * docs `cops_cuda/README.md` (design) and `cops_cuda/kernels/README.md`.
* **Contract (item 1):** `evaluate(ParameterState, ColumnRequest) -> EvaluationResult`. Columns
  come back in `Model.evaluate`'s `J` order and only requested columns are computed. With
  `ModelFacade`/`install`, the unmodified `fit_cops.Problem` runs on any backend. Through the
  facade, residuals and Jacobians (data and regularized) are **bitwise identical** to direct
  evaluation in 5 configurations: legacy d4 free T/κ, d3 ordered T, d2 fixed T, the L4
  control, and L5 + 2 templates. Unsupported settings are refused before evaluation.
* **Layout (items 2 and 4):**
  * SoA layout: int32 faces (range-checked), FP64 Kn/KPhi/KV/cos, uint8 shadow (checked 0/1),
    **int64** CSR offsets.
  * Mesh vectors are stored component-major: `Y[25, nfaces]`.
  * File record order is preserved, and so is the caller's row order.
  * Device fields are checked against the declared FP64 policy. Reference NVRTC options
    are `--fmad=false`; fast-math is refused.
* **Memory (item 3):**
  * Exact `MemoryPlan` per array, with CuPy's 512 B rounding.
  * Admission against `free − 20 % total − slack`, checked against both the first and a
    live measurement.
  * The pool limit enforces the plan.
  * Tiling splits observations at point boundaries and parameter columns into groups.
* **Provenance (item 5):** the backend manifest records:
  * application and reference tree digests and package versions;
  * GPU UUID/PCI and compute capability;
  * CUDA runtime, driver and NVRTC versions;
  * the precision policy;
  * digest, options and architecture of every compiled CUDA source;
  * the dataset identity: table, elements, mesh identity, selection digests, loaded set
    digest, per-array digests, device audit, plan and timings.
* **Gate evidence:**
  * `cops_output/cuda_c02/load_20261020T021843Z/report.json`: ALL PASS. The run took
    125 s on GPU 0 with host load average 2, beside the running L5 CPU fit.
  * **Subset** (47 C00 points, 8,783,677 records): every file matches its C00 digest, and
    the device arrays are bitwise equal to `fit_cops.read_kernel` output, an independent
    reader.
  * **Load/release cycles** ×3: identical per-array digests every cycle; the pool is empty
    after each release; free-memory drift is **0 bytes**.
  * **Full set** (1150 points, 245,515,533 records):
    * the loaded set digest equals `run_manifest.json` (`908dbbf4…2800`);
    * the device audit is all zero, and the round trip is bitwise;
    * the plan is 9,248,031,744 B, pool in use 9,248,031,744 B, and measured device use
      9,261,023,232 B. nvidia-smi reports 9.53 GB for the process, which includes the
      ~0.28 GB context. After release, free memory equals the value before the load.
    * header scan 1.0 s, stream 65.3 s (150 MB/s, bound by disk and SHA-256), device
      audit 0.013 s, round trip 26.5 s.
  * **Capacity:** a plan one byte larger than usable is refused in 1.0 s, after the header
    scan and before any allocation; the pool stays empty.
  * **Second GPU** (GPU-c314803e…): the subset is bitwise identical to GPU 0.
  * **Synthetic refusals** (pytest), each with no device memory left behind:
    * truncated, appended, missing and wrong-version files;
    * wrong iPoint; out-of-range or duplicate faces, including a duplicate split across
      chunks; non-finite values; invalid shadow flags;
    * header/table mismatch; kernel-set digest mismatch; expected-digest mismatch;
    * a reordered mesh, refused both against the manifest and against an expected mesh
      identity; a basis built for another mesh; an FP32 basis;
    * too little capacity; an allocation outside the plan;
    * device-side corruption: face, NaN, shadow and offsets.
  * `pytest`: 78 passed (CUDA environment). The host-side modules also run under the
    Python 3.8 reference interpreter, which has no pytest; this was checked by import and
    execution.
* **Limitations and next prerequisites:**
  * No pressure kernels yet; C06 implements `PressureBackend` on `DeviceDataset`.
  * `plan_tiles` budgets only the partial sums of the fused C06 kernels. C06 must declare
    any additional scratch in its plan.
  * Streaming is single-threaded and limited by disk and hashing. Revisit in C10 only if
    load time matters.
  * Multi-process context reserves (C07) are declared through the `context` category but
    not yet exercised.

## C01: isolated environment and device diagnostics. Status: COMPLETE (2026-10-19)

* **Implemented:** `requirements-cuda.in`, `requirements-cuda.lock` (hash-locked),
  `cops_cuda/environment.py`, `tools/cuda_doctor.py`, `tests/test_cuda_environment.py`,
  `tests/test_provenance.py`, `tests/test_c00_shapes.py`, `RUN.md`.
* **Environment:** Python 3.11.14 (uv), CuPy 13.6.0 (`cupy-cuda12x`, the only CuPy
  distribution), NumPy 1.24.4 and SciPy 1.10.1 (identical to the CPU reference), system
  CUDA 12.4. Resolved libraries: libcuda 550.163.01, libnvrtc 12.4.127, libcublas 12.4.5.8,
  libcusparse 12.3.1.170. CuPy-linked runtime 12.9 on the 12.4 driver.
* **Devices:** GPU-1dbc0c22-5f6f-535c-a083-41a5492478b9 (PCI 1A:00.0) and
  GPU-c314803e-beea-8763-dad3-d4af458e4686 (PCI 68:00.0). Both RTX A6000, cc 8.6, 84 SMs,
  51.0 GB, ECC off, NVLink NV4.
* **Evidence:**
  * `cops_output/cuda_doctor/doctor_20261019T065452Z.json`: PASS.
  * `doctor_20261019T065601Z.json`: PASS from a **fresh** environment rebuilt from the lock
    file with an empty CuPy kernel cache.
  * `pytest`: 25 passed.
* **Measurements:**
  * FP64 add/mul/div/fma/sqrt are exact; log, exp, atan, atan2 and acos are within 2 ulp.
  * Pinned H2D: 8.9–12.1 GB/s (host loaded).
  * NVLink peer copy: 52.8 GB/s, bit-exact.
  * The fixed-order segmented-sum kernel is bitwise repeatable, with max error
    0.31·eps·Σ|t| against `math.fsum`.
* **Finding that affects C06:** cuSPARSE SpMV/SpMM (`cupyx.scipy.sparse`) are **not**
  bitwise repeatable between identical calls. The pressure and Jacobian reductions use
  custom fixed-order kernels.
* **Limitations:** NCCL is not installed and not needed. Peer access is optional in the
  doctor. The doctor's transfer bandwidths were measured on a heavily loaded host.

## C00: freeze scientific and numerical reference. Status: COMPLETE (2026-10-19)

Capture `cops_output/cuda_reference/c00_20261019T064536Z/`.

* `identities`: done. All 1150 kernels (245,515,533 records) match the manifest set digest
  `908dbbf4…2800`; 0 dependency failures; reference tree and git state recorded (HEAD
  27cca46, 31 dirty/untracked paths); the campaign `polyhedral_gravity` binary is
  byte-identical to the PyPI cp38 wheel; omega =
  (-1.1757e-08, -1.4459e-07, 1.44781e-04) rad/s.
* Campaign inventory:
  * Fit references: the two L4 controls are complete and converged, but have **zero free
    parameters**.
  * The L5 harmonic scan was running at capture time, so it was not a reference. It has
    since finished (converged, 17:24 UTC). It is a fixed-baseline harmonic result, which is
    secondary under plan revision 2.
  * The 788k template build was running at capture time and is affected by the HEAL
    defect. Its products are now quarantined as `*_heal_defect`.
  * The 100k bundle was rejected.
* `masks`: done. Late interval: 95 points, 379 NG-only facets (area fraction 4.94e-4); the
  per-point predicate union equals the reference mask exactly. Control interval: empty mask
  (a negative case).
* `fit`: done.
  * 47 subset rows (38 NG and 47 RG used), 8,783,677 contributions, 6 configurations.
  * Configurations: legacy d4 with free T and κ, d3 with ordered T, d2 with fixed T,
    the L4 control, L6 harmonic, and L5 + 2 synthetic templates.
  * Coverage: clipping from 0.9% to 43% of contributing facets, the T-floor branch on 765k
    contributions, amplitudes on their bounds, and the T75 = Tmin edge.
  * C++ replay ≤ 9.1e-10.
  * CPU summation error ≤ 216 eps·M against `math.fsum`.
  * Central differences ≤ 6.1e-6 on non-crossing columns; branch-crossing columns are
    flagged.
* `gravity`: done.
  * Synthetic shapes: cube, tetrahedron and icospheres L0–4, with translation, rotation,
    density scaling, far-field and refinement checks.
  * NaN cases: the package returns NaN at some exact vertex/edge points. This is recorded
    as reference behavior.
  * 788k body: 832 targets at 0.1 and 0.05 m offsets plus 32 on-surface points. The
    topology audit passed and the HEAL replica flags 46%, so the reference uses
    DISABLE. The offset-halving maximum is 0.0050°, the face-order effect ≤ 9.5e-14
    relative, and evaluation took 0.44 s per target on 6 contended CPUs.
* `VALIDATION.md` and `validation/tolerances.json`: **approved by the user on 2026-10-19**,
  with one change made at sign-off: the objective budget is now derived from the pressure
  budget, |Δcost| ≤ Σ|r_i|δ_i + ½Σδ_i², replacing a fixed 1e-12 relative bound. The C03
  gravity amendment (16 eps·M, plus 1e-12 relative on campaign targets, plus NaN parity) is
  approved too. Still provisional: the Jacobian budget (decided at the end of C06) and the
  end-to-end fit criteria (decided after the first CPU-vs-GPU fit).
* **C00 gate: closed.** Fixture generation is reproducible (exclusive stage directories,
  every input content-hashed) and the comparison policy is approved.
  * The full cliff-mask equivalence (C05) needs a corrected cliff bundle. The CPU builder is
    now repaired; the production gravity acceptance route is the C04 decision.
  * The CPU tree changed after this capture, so a C00b re-freeze with B29/N30/C30/NC31
    states is required (ROADMAP C00b).
  * C06 can start on the fit fixtures.
* **Pending:** production fit vectors, i.e. inspected results of fits with free parameters.
