# `cops_cuda`: the LiouvilleCOPS GPU backend package

This package holds the code that moves LiouvilleCOPS data to the GPUs and,
from C03/C06 onward, evaluates gravity, pressures and Jacobians there. The CPU
application in `../LiouvilleCOPS` remains the definition of every quantity; it
is bound explicitly (`reference.py`) and is never copied.

Two words recur with different meanings:

* **transport kernel**: the per-point Kn/KPhi/KV records in
  `kernels/point=NNNN.bin` of a base run. These are handled by
  `device_store.py`.
* **CUDA kernel**: a `__global__` function in `cops_cuda/kernels/*.cu`. These
  are handled by `cuda_sources.py`.

## Modules

| Module | Stage | Interpreter | Role |
|---|---|---|---|
| `errors.py` | C02 | 3.8 + 3.11 | Refusal categories: `CacheFormatError` (malformed bytes), `IdentityError` (wrong product, e.g. another mesh order), `CapacityError` (memory), `UnsupportedConfiguration` |
| `provenance.py` | C00 | 3.8 + 3.11 | File, array and tree digests; exclusive writes; process context |
| `reference.py` | C00 | 3.8 + 3.11 | Explicit binding to `../LiouvilleCOPS/tools` (or `$LIOUVILLE_COPS_REFERENCE`) |
| `formats.py` | C02 | 3.8 + 3.11 | `elements.bin` / `point=NNNN.bin` layouts as NumPy dtypes, with vectorized validation equivalent to `cops_io` |
| `precision.py` | C02 | 3.8 + 3.11 | The FP64 policy: device dtypes per field, NVRTC options (`--fmad=false`, no fast-math) |
| `memory.py` | C02 | plan: both; budget: 3.11 | Memory plans, the 20% headroom admission, pool-limit enforcement, observation/column tiling |
| `backend.py` | C02 | 3.8 + 3.11 | Backend contract, the CPU reference backend, and the facade that runs `fit_cops.Problem` on any backend |
| `cuda_sources.py` | C02 | 3.11 | Compiles `kernels/*.cu` through NVRTC and records source digest, options, NVRTC version and architecture |
| `device_store.py` | C02 | 3.11 | Validated, planned, streamed SoA upload of one run (mesh + transport kernels) to one GPU |
| `manifest.py` | C02 | 3.11 | Software, device, runtime, precision and CUDA-source identity sections for manifests |
| `gravity_body.py` | C03 | 3.8 + 3.11 | Gravity polyhedron: BDF reading through the reference reader, the exact topology audit that replaces HEAL, identity |
| `gravity.py` | C03 | 3.11 | **GPLv3-derived.** Face precompute (`prepare`), bitwise equal to the package; `GravityEvaluator` (fold/split variants, memory plan, launch bound, final scaling) |
| `gravity_cache_io.py` | C04 | 3.11 | Writes a GPU self-gravity field as a CPU-format gravity cache, using the CPU's own `make_definition`/`write_cache` and truthful CUDA metadata, so the unchanged builder and selector reuse it |
| `masks.py` | C05 | 3.11 | GPU late-interval visibility masks (`late_ng_not_rg`, `late_rg_not_ng`, plus the seen diagnostics) by integer OR; bit-identical to the CPU predicates |
| `restart_bundle.py` | C05 | 3.11 | The restart template bundle (NG-only, RG-only, cliff columns from the admitted cache), written with the CPU's bundle writer and checked by the unchanged loader and selector |
| `evaluator.py` | C06 | 3.11 | `CudaBackend`: FP64 pressures and analytic-Jacobian columns behind the C02 contract, for B29/N30/C30/NC31 and the labelled extensions. `"sequential"` mode is bitwise the CPU reference; `"tree"` mode is the deterministic production mode |
| `host_model.py` | C06 | 3.11 | Header-only subclass of the reference `fit_cops.Model`: per-point data and constants bitwise the reference's, with no kernel records on the host. Lets the unmodified `fit_cops.Problem` run on the full 788k selection with a GPU backend installed; CPU evaluation refuses |
| `gpu_model.py` | C06→C07 | 3.11 | GPU-backed subclass of the reference `fit_cops.Model` (header-only host side; `evaluate`, Jacobian and `seen_faces` on the device). `tools/gpu_fit.py` puts it in place of `fit_cops.Model` inside one process, so the unchanged `fit_cops.main()` fits on the GPU. Backends are shared across models with equal rows |
| `multi_gpu.py` | C07 | 3.11 | One fit on several GPUs: `MultiGpuBackend` behind the C02 contract, with persistent spawned workers (one per device, contiguous row blocks balanced by records), one host activity field in shared memory, and assembly in row order. Bitwise equal to one GPU; clean `BackendError` on worker failure |
| `restart.py` | C12 | 3.11 | Fit restart files: write/read (digests, mesh identity), the start vector for any `fit_cops.Problem` by name, state pressures on the GPU (`StateEvaluator`). Used by `gpu_fit.py` and `tools/restart_tool.py` |
| `environment.py` | C01 | 3.11 | CUDA library and device discovery; thread policy for workers |
| `kernels/` | C02+ | — | CUDA C sources; see `kernels/README.md` |

"3.8 + 3.11" modules are also imported by the C00 capture under the CPU
reference interpreter. They must not import CuPy at module level.

## C02 design

### 1. Backend contract (`backend.py`)

```text
backend.evaluate(ParameterState(a[25], T[3], kappa, b[nb] or None),
                 ColumnRequest(basis_index, enhancement_index, temperature, kappa))
    -> EvaluationResult(pNG[npt], pRG[npt], columns{key: (dpNG, dpRG)})
```

* **Column order** is exactly the insertion order of `fit_cops.Model.evaluate`'s
  `J` dictionary: `("a", j)` in the caller's order, `("b", k)` in the caller's
  order, `("T", 0..2)` (Tmax, T75, Tmin), then `("kappa", 0)`. Only the
  requested columns are computed and returned.
* **What stays on the CPU:** the mapping from the optimizer vector `x`
  (`Problem.unpack`), log residuals, gauge weights, regularization, SciPy
  TRF and every diagnostic. A backend returns only FP64 pressures and
  pressure derivatives.
* **`ModelFacade` / `install(problem, backend)`** replace `problem.M` with a
  facade, so the *unmodified* `fit_cops.Problem` evaluates through a backend.
  `tests/test_backend_contract.py` shows that residuals and Jacobians through
  the facade are bitwise identical to direct evaluation, for legacy (free,
  ordered and fixed temperature) and enhancement (control and templates)
  configurations.
* **Support checks:** `capabilities()` and `check_supported()` make a backend
  refuse unsupported settings (SZA mode, heliocentric model, RG pressure model,
  derivative kinds) before evaluating. The run-level settings are parsed by
  `ModelSettings.from_header` with the reference's own header parsers.
* **The C06 oracle:** `ReferenceBackend` wraps `fit_cops.Model` unchanged. C06
  compares the CUDA backend with it through the same interface.

### 2. Device layout (`device_store.py`)

Structure-of-arrays, FP64 for every physical value:

| Array | dtype | Shape | Notes |
|---|---|---|---|
| `kernel.face` | int32 | [nrecords] | checked to lie in `[0, nfaces)`; `nfaces < 2^31` is checked |
| `kernel.kn`, `kernel.kphi`, `kernel.kv`, `kernel.cos` | float64 | [nrecords] | bit-identical to the file |
| `kernel.shadow` | uint8 | [nrecords] | file int32, checked ∈ {0,1} before narrowing |
| `kernel.offsets` | **int64** | [npt+1] | CSR row pointer; 64-bit as ROADMAP C02 requires |
| `kernel.geometry` | float64 | [10, npt] | kernel header: et, V·b_RG, sc xyz, sun xyz, r_h [m], altitude |
| `mesh.Y` | float64 | [25, nfaces] | component-major: coalesced per-facet reads |
| `mesh.centroid`, `mesh.normal` | float64 | [3, nfaces] | component-major |
| `mesh.area`, `mesh.AH2O`, `mesh.ACO2` | float64 | [nfaces] | |
| `basis.enhancement_basis` | float64 | [ncols, nfaces] | harmonics and templates; tied to the mesh identity |

Records keep the kernel-file order within a point and the caller's row order
across points. That is the order in which the reference sums. On the device a
record takes 37 B, so the 245,515,533-record campaign set is 9.08 GB.

### 3. Loading sequence and refusals

1. **Identities.** The table and run manifest are checked, and `elements.bin`
   is fully validated. A SHA-256 differing from `run_manifest.json` means the
   wrong or a reordered mesh. The mesh identity is computed with
   `activity_templates.mesh_identity`, the definition template bundles carry,
   and compared with an optional expected identity.
2. **Header scan.** All files are scanned: version, iPoint, count and exact
   size, and header geometry versus table row (the `cops_io` tolerances).
   This yields exact record counts and offsets.
3. **Plan admission.** If the plan doesn't fit, `CapacityError` is raised with
   the full plan, before any allocation.
4. **Streaming.** Records are read into a reusable buffer, validated (face
   range, duplicates including duplicates split across chunks, finite values,
   shadow flag) and hashed, then de-interleaved into two alternating pinned
   staging sets and copied asynchronously to their final offsets.
5. **Identity of the loaded set.** Per-file SHA-256 are compared with
   optional expected digests (the C00 `identities/kernel_files.json`). For a
   full table load, the set digest must equal `run_manifest.json`
   `kernels.set_sha256`.
6. **Device audit.** `kernels/store_audit.cu` recounts the invariants in
   device memory using race-free integer atomics.
7. **Round trip.** `verify_roundtrip()` downloads every array and compares
   its SHA-256 with the host digest from the upload, bit for bit.

Any failure frees everything already allocated. That includes arrays still
referenced by the exception's traceback; the loader clears those frames. A
partially loaded dataset never exists.

### 4. Memory policy (`memory.py`)

* Every device array is declared in a `MemoryPlan` with its exact size and a
  category (static, basis, partial sums, scratch, slack, modules, context).
  Pool items are rounded to CuPy's 512-byte granule.
* `DeviceBudget` measures free and total memory with `cudaMemGetInfo` after
  context creation, so this process's own context, about 0.3 GB on an A6000,
  is already excluded. A plan is admitted only if it fits
  `free − 20% × total − slack − earlier reservations`. This is checked both
  against the first measurement and against a live measurement, because
  another process may allocate in between.
* The budget sets the CuPy pool limit to the admitted plans plus a declared
  64 MiB slack for small temporaries. An allocation outside the plan raises
  at the allocation site instead of quietly consuming the headroom.
  `allocation_guard` turns it into a `CapacityError` that shows the plan.
* Release frees the pool blocks (`free_all_blocks`), so measured device free
  memory returns to its starting level. The C02 gate checks this across
  repeated cycles.
* CPU RSS is never used as a GPU estimate.
* **Tiling** (`plan_tiles`) splits observations into contiguous point ranges,
  only at point boundaries, and parameter columns into groups. The
  partial-sum scratch of the C06 reduction then fits a declared budget; no
  contribution × parameter array is ever formed. A point larger than the tile
  bound is reduced by several blocks whose partial sums are combined in a
  fixed order.

### 5. Precision and provenance

* `precision.py` declares the FP64 policy, and each array is checked against
  it at allocation. Reference CUDA modules are compiled with `--fmad=false`;
  fast-math options are refused (see the module docstring for the reasons).
* `manifest.backend_manifest()` records:
  * the content digests of this package and of the reference tree;
  * the Python, NumPy, SciPy and CuPy versions;
  * device UUID, PCI bus, compute capability and memory;
  * CUDA runtime, driver and NVRTC versions;
  * the precision policy;
  * the digest, options and architecture of every compiled CUDA module,
    plus the digests of all shipped `.cu` files;
  * the dataset identity: table, elements, mesh identity, selection digests,
    kernel set digest, layout with per-array digests, device audit, memory
    plan and timings.

## Using it

```python
from cops_cuda import device_store
with device_store.DeviceDataset.load(run_dir, rows=M.rows, device_id=0) as ds:
    ds.verify_roundtrip()          # optional bitwise check
    ...                            # C06: CudaBackend(ds, settings, basis)
```

Gravity (C03):

```python
from cops_cuda import gravity, gravity_body
body = gravity_body.GravityBody.from_bdf(bdf_path, density=533.0)   # audited; no HEAL
with gravity.GravityEvaluator(body, device_id=0) as ev:
    V, g, M = ev.evaluate(points, variant="fold", with_magnitude=True)   # SI units
```

The translation details are in `../GRAVITY_REFERENCE.md` §8 and the budget is in
`../VALIDATION.md` (C03 amendment).

Commands and the evidence tools are in `../RUN.md`; results are in `../STATUS.md`.
