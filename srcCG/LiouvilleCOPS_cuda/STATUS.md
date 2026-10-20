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
* **Evidence:** `cops_output/cuda_c06/check_20261019T193255Z/report.json`, ALL PASS, on
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
  `cops_output/cuda_c06/full_20261019T221945Z/report.json`, ALL PASS):
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
  (`cuda_c06/jacobian_budget_20261020T000812Z`). The first report (329) used KT instead of
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
* **Evidence:** `cops_output/cuda_c05/masks_20261019T194835Z/report.json`, ALL PASS.
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
  * Product: `cops_output/cuda_c05/restart_20261019/restart_templates_shap5_788k.{npz,json,png}`,
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
  `cuda_c04/equivalence_20261019T185043Z/`.
* **CPU reference:** the unchanged repaired builder, on CPUs 28–35 beside the user's CPU
  campaign: 5,617 s of gravity.
* **Full-field equivalence: PASS** (`equivalence_20261019T185043Z/compare_*.json`), on all
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
    offset-halving 0.00531°) was reused with `--gravity-cache` ("reused gravity cache …
    same reference omega"). Its SPICE ω is bit-identical to the tool's.
  * It wrote `cpu_builder_reuse_of_gpu_cache_test.npz`: 379 NG-only facets; cliff ≥45° on
    16.2% of the area and ≥60° on 5.0%.
  * The **unchanged fitter's `select_templates` accepts its `cliff_ge_60deg` column**, and
    the metadata states the CUDA backend.
  * The restart path therefore needs no change in `../LiouvilleCOPS`.
* **788k production equivalence: PASS**
  (`cuda_c04/equivalence_20261019T221538Z/compare_20261019T230910Z.json`). The CPU reference
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
    `cuda_c04/equivalence_20261019T185043Z/admitted_gravity_cache_shap5vol100k.npz`
    (offset-halving 0.00531°);
  * **788k production:**
    `cuda_c04/equivalence_20261019T221538Z/admitted_gravity_cache_shap5_788k.npz`
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
* **Gate evidence:** `cops_output/cuda_c03/gravity_20261019T171225Z/report.json`.
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
  * `cops_output/cuda_c02/load_20261019T165013Z/report.json`: ALL PASS. The run took
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
