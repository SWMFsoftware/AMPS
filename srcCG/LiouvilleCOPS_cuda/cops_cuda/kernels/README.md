# CUDA C sources

Every `__global__` function the backend runs lives in a `.cu` file in this
directory. It is compiled at run time by NVRTC through
`cops_cuda.cuda_sources.module(name, device_id)`. Sources are never built from
Python strings at the call site, so the SHA-256 recorded in a manifest
identifies exactly the code that ran.

| File | Stage | Functions | Purpose |
|---|---|---|---|
| `store_audit.cu` | C02 | `audit_records`, `audit_offsets` | Recount transport-kernel invariants in device memory after upload |
| `pressure.cu` | C06, C08 | `precompute_static`, `eval_sequential`, `eval_partial`, `eval_finish`; C08 `eval_fused_p5`, `eval_fused_p14`, `eval_fused_c1`…`c16` (register accumulators, bitwise equal to `eval_partial`) | Per-contribution terms of `fit_cops.Model.evaluate` and their per-point sums (pressures, temperature sums, activity/template columns); verification and production reductions |
| `gravity.cu` | C03 | `gravity_fold`, `gravity_partial`, `gravity_reduce` | **GPLv3** translation of polyhedral-gravity 3.3: exact FP64 potential and acceleration, plus the term magnitude M |
| `masks.cu` | C05, C07 | `visibility_union`, `seen_positive` | Late-interval visibility masks (integer OR, bit-identical to the CPU predicates); facets with Kn > 0 for `seen_faces` |
| `restart.cu` | C12, C13 | `facet_point_shares`, `point_change_partials`; C13 `first_structural_row`, `first_interval_scores`, `first_significant_row` | Per-facet influence of an activity change on point pressures (`restart_tool.py enhance`); late-interval fine tuning: earliest structural/significant row per facet (integer atomicMin) and first-interval NG/RG dominance scores (one launch per row, chronological; `tools/fine_tune_gpu.py`) |
| `render.cu` | figures | `raster_faces`, `ray_first_hit` | Figures only, no model quantity: z-buffer rasterization of the shape model and the first ray–mesh hit (the sub-spacecraft intercept) for `cops_cuda/figures.py` |

## Conventions

* **Compiler options** come from `precision.REFERENCE_NVRTC_OPTIONS`:
  `--std=c++14 --fmad=false`. `precision.check_nvrtc_options` refuses
  fast-math, flush-to-zero, approximate division/sqrt and `--fmad=true`. Where
  a kernel deliberately fuses an operation, it calls `fma()` explicitly and
  states why in a comment.
* **64-bit indexing**: record indices and CSR offsets are `long long`. Every
  loop over records is a 64-bit grid-stride loop, so the grid size never has
  to grow with the data.
* **Determinism**:
  * Floating-point results must not depend on thread scheduling. Use
    fixed-order (segmented, two-pass) reductions and no floating-point
    atomics.
  * Integer atomics, which are associative, are allowed for counters and
    Boolean OR (masks, audits).
* **Comments** explain the indexing, the branch logic, the reduction order and
  the precision choices, with a reference to the formula in
  `../LiouvilleCOPS/tools/fit_cops.py` or `GRAVITY_REFERENCE.md`. They do not
  restate expressions.
* **Licensing**: a file translated from GPL source carries the GPLv3 header and
  attribution required in `GRAVITY_REFERENCE.md` §7. All other files are
  original to this project.
