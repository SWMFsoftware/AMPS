# LiouvilleCOPS CUDA + CuPy development roadmap

Prepared 2026-10-19. Status: proposed implementation roadmap; no CUDA backend
has been implemented or validated by this document.

## 1. Objective and scope

Port the expensive gravity, masking, pressure, and analytic-Jacobian
calculations of the self-contained `LiouvilleCOPS` application to CUDA,
orchestrated from Python through CuPy. Keep SciPy's bounded trust-region
least-squares optimizer and the existing scientific definitions, data-selection
profiles, regularization, output reports, and C++ replay checks.

The initial target is two NVIDIA RTX A6000 GPUs and the converted SHAP5
`cg_spc_shap5_788k_cart.bdf` surface: 787,754 retained triangles, 1,150
observation points, and 245,515,533 cached point–facet contributions in the
audited base. These are dataset-specific counts, never program constants.

This directory is the development destination. The CPU application in
`../LiouvilleCOPS` is the reference. Its source, manifests, exact input files,
and products must be identified by hashes before comparison. The application
has no dependency on unrelated code in `srcCG`.

Scope includes exact homogeneous-polyhedron self-gravity, effective gravity,
slope and NG-only templates, their persistent caches, and the cached-kernel
fitting evaluator. Initial transport-kernel construction, ray tracing, and
SPICE calculations remain CPU operations. GPU transport/ray tracing would be
a separate later project with its own geometry validation.

## 2. Audited environment and missing dependencies

The following was observed on the current host; re-run discovery on another
machine rather than relying on paths or version assumptions.

| Component | Observed state | Action |
|---|---|---|
| GPUs | Two RTX A6000, compute capability 8.6, approximately 48 GiB each | Use stable GPU UUIDs in manifests; validate both devices |
| Interconnect | `nvidia-smi topo -m` reports NV4 | Measure peer access; do not assume one contiguous 96 GiB allocation |
| Driver | 550.163.01 | Verify runtime compatibility using actual kernel execution |
| CUDA toolkit | 12.4; nvcc 12.4.131 | Available, including CUDA headers |
| CUDA libraries | cudart, NVRTC, cuBLAS, cuSPARSE, cuSOLVER, cuRAND, cuFFT, nvJitLink found | Record resolved library paths and runtime versions |
| Python | 3.8.10 | Create an isolated newer environment |
| NumPy / SciPy | 1.24.4 / 1.10.1 in current Python | Preserve CPU reference environment; explicitly pin GPU environment |
| CuPy | Not installed | Install exactly one matching CUDA wheel package |
| NCCL | Not found in audited library discovery | Optional; unnecessary for independent observation rows returned to CPU |
| Native build | GCC/G++ 9.4, CMake 3.16.3, Python headers, pybind11 2.13.6 | Sufficient starting tools; test actual compilation if native extension is used |
| NVIDIA HPC compiler | `nvc++` not found on PATH | Not needed for CUDA/CuPy |
| Host CPU | i9-10980XE, 18 cores / 36 logical CPUs | User requests 36 threads for appropriate CPU work; avoid nested oversubscription |

The dependencies are therefore **not all installed**. CUDA development tools
are available, but CuPy and a suitable isolated Python environment are missing.
Finding shared libraries does not establish successful kernel compilation,
device execution, or compatibility of a chosen package set.

An explicitly verified documentation candidate is Python 3.11 with
`cupy-cuda12x==13.3.0` and CUDA 12.4. CuPy 13.3 documents Python 3.9–3.12;
installing that package into the present Python 3.8 environment is not the plan.
Treat this as a reproducibility candidate, not a recommendation to ignore
newer maintained releases: stage C01 must select, test, and lock the actual
supported package set. Include NumPy, SciPy, Matplotlib, the CPU gravity
reference, and test dependencies. Do not install PyTorch, JAX, cuDNN, or NCCL
unless a concrete implementation needs them.

## 3. Scientific behavior that must be preserved

### 3.1 Gravity and slope

The current builder uses `polyhedral-gravity==3.3`, its Tsoulis homogeneous
polyhedron formulation, and outward-oriented geometry in metres. Self-gravity
is evaluated at `centroid + offset * outward_normal`; the current campaign
uses density 533 kg/m^3 and offset 0.1 m. Preserve all signs, units, density,
mesh ordering, evaluation position, and gravitational constant from the
pinned reference. Verify those constants from its source before coding.

At a facet centroid r, the current effective field is

```text
g_centrifugal = -omega cross (omega cross r)
g_effective   = g_self_at_offset + g_centrifugal_at_centroid
slope         = acos(clamp(dot(-normal, g_effective) / |g_effective|, -1, 1))
```

The centroid/offset distinction is part of the current algorithm. Preserve it
in equivalence tests; changing it would require a separately documented model
change. The angular velocity comes from SPICE in the mesh body frame. Preserve
the existing sampling/variation audit and record its epoch and kernel identities.
The current constant-omega realization must not silently become a new thermal
or dynamical model.

Cliff templates currently use binary slope thresholds of 45 and 60 degrees.
Use the 60-degree template for the main hypothesis; the 45-degree template is
a separate sensitivity fit. The data do not supply a calibrated continuous
production-versus-slope law. Density and slope thresholds are prescribed
physical/hypothesis settings, not new fitted parameters.

The older 100k gravity surface was rejected for this 788k campaign because
its sampled slope difference against a second surface exceeded the existing
5-degree p95 gate. The accepted approach uses the exact 788k transport surface
as the gravity boundary. Do not reintroduce a coarser boundary as an unnoticed
GPU optimization. Retain the outward-offset halving audit (current 512-point
sample, maximum slope change 0.25 degrees).

### 3.2 Visibility and masks

`late_ng_not_rg` is the union over selected late observation points of facets
with `Kn > 0` and, at that SAME point, `KPhi == 0` and `KV == 0`. It is not
the set that never contributes to RG over the entire interval. The current
campaign uses the closed interval 2016-09-30T10:33:10 to 10:39:26 UTC.

Compute this predicate from the original FP64 kernels before any lossy
conversion. A tiny positive value rounded to zero must never change visibility.
Use integer/Boolean reductions with race-safe OR, not concurrent unsynchronized
Boolean stores. Response-weighted diagnostic maps are separate from the Boolean
mask; they must not become thresholds for defining the source region.

NG-only masks depend on observation geometry and interval as well as shape.
They cannot be reused merely because the BDF filename matches. Cliff masks
depend on gravity, omega, and slope threshold. Self-gravity is epoch-independent
only for fixed geometry, density, and evaluation offset.

### 3.3 Fitting

Preserve `tools/fit_cops.py` as the definition of the evaluator: species
sensitivities, production law, shadow and SZA handling, heliocentric scaling,
analytical temperature branches, gas-mass constants, NG conversion, each RG
pressure model, spacecraft-velocity correction, and CO2 scaling.

Surface activity is the sum of the linear legacy map and scaled harmonic and
template enhancements, followed by ONE nonnegative clip. Derivatives use the
same active-set branch (`linear_activity > 0`). Harmonic coefficients can be
signed; the proposed physical template amplitudes have nonnegative bounds.
The reference supports both legacy fits and fixed-physics enhancement fits;
the GPU implementation must explicitly support or reject each configuration.

Residuals are log pressure ratios with each gauge divided by sqrt(its selected
sample count). Preserve the pressure floor, selection profile, saturation
policy, exclusion intervals, ordered-temperature option, parameter scaling,
bounds, regularization Gram matrix including cross terms, lambda, and analytic
Jacobian. Preserve SciPy TRF, `x_scale='jac'`, stopping tolerances, and evaluation
limits for equivalence runs. Do not replace the objective with a neural-network
optimizer or silently change the treatment of saturated data.

Regularization, parameter vectors, final residuals/Jacobians, SVD/rank checks,
uncertainty diagnostics, and reporting initially remain CPU FP64. Clipping and
temperature floors create nondifferentiable points: finite differences must
distinguish branch crossings from erroneous derivatives.

## 4. Proposed architecture

All names below are proposed deliverables, not existing commands.

```text
LiouvilleCOPS_cuda/
  ROADMAP.md
  README.md, RUN.md, VALIDATION.md, PERFORMANCE.md
  requirements-cuda.lock
  cops_cuda/
    environment.py, backend.py, device_pool.py
    kernel_store.py, model.py, gravity.py, masks.py
    precision.py, provenance.py, checkpoint.py
    kernels/gravity.cu, pressure.cu, jacobian.cu, masks.cu
  tools/
    cuda_doctor.py, validate_cuda.py, benchmark_cuda.py
    build_physical_activity_templates.py, run_fit_profile.py
  tests/
    test_gravity.py, test_masks.py, test_model.py
    test_derivatives.py, test_cache.py, test_multi_gpu.py
```

Choose in C00 how shared CPU readers/profiles/reporting are packaged or copied
with provenance. Avoid accidental imports based on current working directory.
The CUDA application should be deployable with explicit shared dependencies,
without unrelated `srcCG` code. Preserve the original CPU application for replay.

Use CuPy for device ownership, arrays, streams, transfers, and kernel compilation;
use CUDA RawKernel/RawModule implementations for fused indexed reductions.
Start with one device. Introduce persistent spawned workers, one per GPU, for
multi-GPU evaluation; never fork a process after initializing CUDA. Workers
load their own assigned slices and return only small results. Exchange parameters
per iteration; keep static kernels resident where memory permits.

Retain the repository-relative `<AMPS>/cops_output` default and existing
`COPS_OUTPUT_ROOT` override. On this host the campaign override is
`/data/vtenishe/rosetta_data_end_of_mission_cops/cops_output`. Input locations
are configurable; the current L2 and source roots are respectively
`/data/vtenishe/rosetta_data_end_of_mission_cops/data_COPS_L2` and
`/data/vtenishe/rosetta_data_end_of_mission_cops/sources`.
Write generated products outside both source directories, in unique run folders.

## 5. Development stages and acceptance gates

Each stage must deliver commented code, focused tests, documented commands,
machine-readable measurements, and an updated status entry. A stage is complete
only after its acceptance evidence exists; mark unavailable production evidence
as pending rather than inferring success from a synthetic example.

### C00 — Freeze scientific and numerical reference

1. Record source revision and dirty-file hashes, package versions, base manifest,
   BDF/conversion metadata, kernel hashes, SPICE identities, and all fit settings.
2. Save CPU reference pressures, Jacobians, activity fields, gravity vectors,
   slopes, exact masks, regularization residuals, and representative parameter
   vectors. Include prior, fitted vectors when available, active clipping,
   night/lit points, temperature-floor branches, and template bounds.
3. Capture small deterministic fixtures plus a production-subset index list.
   Inspect pinned gravity source and license; document equations, singular
   cases, orientation handling, and attribution obligations before translating.
4. Identify which current campaign products are complete. Running CPU jobs are
   not accepted fit references until termination and diagnostics are inspected.

Gate: reproducible CPU fixture generation and a signed-off numerical comparison
policy. Any unresolved physical convention blocks its dependent kernel stage.

### C01 — Isolated environment and device diagnostics

1. Create a Python 3.11 environment without changing the running CPU environment.
2. Resolve a compatible CuPy CUDA-12 package and pin the full tested environment.
   Record runtime, driver, NVRTC, toolkit, compiler, and dynamic library paths.
3. Implement `cuda_doctor.py`: enumerate UUIDs/capabilities/free memory; compile
   and execute a tiny FP64 RawKernel on each GPU; test transfers, asynchronous
   completion, error reporting, and optional peer access.
4. Test a small CPU↔GPU pressure reduction and sparse operation if used. Detect
   conflicting CuPy packages and library versions explicitly.
5. Document CPU 36-thread preference while preventing 36 BLAS threads per worker.

Gate: both GPUs pass actual kernel tests and a fresh environment reproduces
them. The CUDA toolkit's presence alone is insufficient. NCCL is optional.

### C02 — Backend contract, memory layout, and provenance

1. Define a backend interface that takes parameter vectors and returns pressures
   and requested analytic-Jacobian columns in the exact existing order.
2. Stream audited point files into structure-of-arrays device buffers, preserving
   face indices and point offsets. Use 64-bit offsets for future larger caches;
   validate narrower index ranges explicitly.
3. Budget arrays, basis, partial sums, scratch, CUDA context, compiler cache, and
   CuPy memory pool. CPU RSS is not a GPU memory estimate. Query free memory and
   reserve declared headroom (initial policy: at least 20%).
4. Keep mesh geometry, parameters, and baseline references FP64. Stage all
   contribution fields in FP64 first. Avoid giant contribution-by-parameter
   arrays; tile parameter columns and observations.
5. Add backend/version/precision/kernel-source/device identities to manifests.
   Validate inputs before costly GPU allocations and fail clearly on OOM.

Gate: load/release cycles preserve data and bounded memory; malformed caches,
wrong mesh order, and insufficient capacity are detected without silent changes.

### C03 — Single-GPU exact polyhedral gravity in FP64

1. Precompute geometry-only terms on CPU from the audited closed polyhedron,
   reproducing the reference orientation and any documented healing effects.
2. Translate the reference acceleration formula to CUDA with explicit units and
   equation references. Use stable log/angle expressions and FP64 throughout.
   Do not substitute point masses, a sphere, or a coarse mesh.
3. Tile source faces/edges and target points; reduce partial vector contributions
   with a fixed hierarchy. Benchmark target-per-block and tiled variants. Avoid
   materializing the full target-by-source matrix.
4. Compute only acceleration when sufficient; potential/Hessian work is optional
   unless needed to preserve the chosen formulation or diagnostics.
5. Test tetrahedron/cube symmetry, translated and rotated copies, density scaling,
   far-field total-mass behavior, and near-face/edge/vertex points against the
   pinned CPU reference. An analytic sphere is a convergence test only for
   successively refined sphere meshes, not exact polyhedron equivalence.
6. Reproduce the existing 512-point offset-halving audit and add deterministic
   difficult-location samples. Record vector and slope errors separately.

Gate: FP64 gravity meets C00 numerical budgets and reproduces valid slope
classifications. Exact gravity remains roughly all-target/all-source work;
CUDA changes execution speed, not asymptotic complexity. Report measurements
before promising a production runtime, especially on A6000 FP64 hardware.

### C04 — Two-GPU gravity, checkpoints, and reusable caches

1. Partition target facets between GPUs; replicate source geometry when it fits.
   Balance by measured chunk runtime, retaining original output ordering.
2. Write completed FP64 target chunks atomically with input hashes, target IDs,
   backend identity, and checksums. Resume only complete verified chunks from
   the same numerical configuration; reject mixed/overlapping corrupt chunks.
3. Finalize the existing gravity-cache schema when compatible, adding backend
   provenance without falsely claiming byte-identical CPU results. Version the
   schema if compatibility requires it. Store self-gravity separately from
   omega-dependent derived fields.
4. Preserve exact cache identity: centroid/normal/area ordering, gravity-shape
   content, density, offset, algorithm, and numerical realization. Changed omega
   recomputes derived fields; changed geometry/density/offset invalidates the
   appropriate cached calculation under the established reuse policy.
5. Exercise interrupted write, restart, corrupt chunk, moved input path, changed
   omega, and shape reordering tests. A filename match never establishes reuse.

Gate: one/two-GPU and resumed calculations agree within the approved error budget;
final caches pass CPU loaders and reuse avoids repeating self-gravity work.

### C05 — GPU masking and effective fields

1. Reproduce exact NG/RG predicates from original kernels; compute same-point
   NG-only union using race-safe integer OR/reductions.
2. Derive centrifugal acceleration and slopes in FP64 using CPU/SPICE omega.
3. Generate >=45-degree and >=60-degree masks with explicit threshold convention.
   Flag facets whose numerical slope uncertainty overlaps a threshold; recompute
   those on the CPU reference rather than silently accepting classification drift.
4. Port response-map sums only if useful; their floating reductions need separate
   tolerances. Area fractions, overlap, counts, and mask identities remain audited.
5. Persist mask bundles tied to geometry, kernel identities, time interval, omega,
   density/offset, threshold, and exact definition. Reuse current template loader
   semantics and output maps.

Gate: visibility masks are bit-identical to CPU; cliff classification differences
are resolved and reported; moving paths preserves valid content-based reuse.

### C06 — Single-GPU FP64 pressure and analytic Jacobian

1. Port the facet activity calculation, single clip, temperatures, source factors,
   and each gauge/species contribution exactly. Preserve RG velocity-term signs.
2. Use segmented reductions over contiguous point ranges. Process long rows in
   multiple blocks and reduce their partial sums deterministically.
3. Implement analytic derivatives for all requested legacy/harmonic/template,
   temperature, and CO2 parameters. Compute only requested columns.
4. Reuse pressure/derivative intermediates for identical parameter vectors where
   safe; cache keys must include every setting that affects the result.
5. For fixed-temperature fits, precompute temperature-dependent factors and
   constant CO2 terms once. This is an exact specialization; do not use it when
   temperatures vary. Preserve clipping, which prevents treating the full
   activity problem as globally linear.
6. Return small FP64 arrays to existing CPU residual/regularization/SciPy code.
   Test both residual-only and joint residual/Jacobian evaluation paths.

Gate: prior pressures pass existing C++ replay tolerance (default 1e-6 relative),
and stronger CPU-vs-GPU budgets from C00 pass across representative vectors.
Use central differences away from branch crossings and branch-aware tests at
boundaries. Passing the replay check alone does not validate all derivatives.

### C07 — Two-GPU fitting and resource scheduling

1. Split complete observations by contribution count and measured work, rather
   than splitting solely by number of points. Keep exact global row ordering.
2. Keep device arrays in persistent spawned workers; send parameters and receive
   pressure/Jacobian rows. Propagate worker failures and cancellation cleanly.
3. Compare one fit across two GPUs against one independent lambda fit per GPU.
   Independent fits may require compact storage to fit each full dataset.
4. Enforce explicit ownership of devices so a lambda worker pool cannot duplicate
   the dataset uncontrollably. Schedule gravity and fitting deliberately; default
   to no simultaneous memory-heavy work on the same device.
5. Use host assembly of the small Jacobian initially; require NCCL only if a
   measured alternative needs collective device reductions.

Gate: same objective/data across device layouts, bounded device memory, stable
cleanup, and measured single-fit latency AND campaign throughput. Two devices
do not imply twice the speed or automatically pooled memory.

### C08 — Optional mixed precision

1. Keep FP64 backend as the reference and default until this stage passes.
2. Evaluate FP32 static contribution storage separately from FP32 arithmetic.
   Test per-point scaling if underflow or extreme dynamic range requires it.
3. Keep geometry/near-surface gravity FP64 initially. Keep parameters, activity
   branch decisions, final reductions, residuals, regularization, SVD, and
   uncertainty calculations FP64. Disable fast-math in the reference kernels.
4. Measure cancellation in signed harmonics and RG terms, pressure/log-pressure
   errors, derivative noise, active-set changes, rank changes, and template-bound
   behavior. Test multiple chunk sizes and devices.
5. Accumulating FP32 products in FP64 does not restore information lost in the
   products. Validate both multiplication precision and summation order.
6. Compare complete fits with unchanged CPU tolerances, including current 1e-10
   stopping settings. If mixed precision cannot support them, retain FP64 or
   introduce an explicitly separate validated tolerance profile. Never weaken
   scientific checks solely to claim speed.

Gate: approved end-to-end error budget passes. If unsupported, mixed precision
remains experimental. GPU acceleration does not require accepting FP32 gravity.

### C09 — Production integration, outputs, and validation campaign

1. Add documented backend/device/precision/memory/chunk controls to runners while
   retaining existing profile, lambda, degree, freezing, bounds, and data options.
2. Preserve output refusal/unique-directory behavior, selection exports, manifests,
   gravity/template products, enhancement exchange files, and comparison plots.
3. Test degree-4 control; degrees 5–8; templates alone where supported; harmonic
   plus NG-only plus 60-degree cliff; and separate 45-degree sensitivity.
4. Exercise all three named selection profiles and legacy options using the exact
   same observations for CPU/GPU comparisons. Include the late NG interval,
   saturation, missing data, and shape-inside decisions.
5. Test synthetic recovery and non-identifiability, regularization lambda scans,
   overlapping templates, clipping, and active bounds. Degenerate parameters may
   differ while predictions agree; report null directions and map differences.
6. Run selected full 10-hour production fits; compare objective components,
   pressures, activity maps, coefficients, bound states, gradient/optimality,
   convergence status, singular values, and uncertainty diagnostics.
7. Replay selected GPU solutions in the original C++ forward model using the
   requested 36 CPU threads and existing validation/plotting workflow.

Gate: numerical and scientific report is complete; GPU backend has an explicit
supported-configuration matrix and fails clearly for unsupported modes.

### C10 — Performance characterization and release

1. Benchmark original CPU, optimized CPU where available, one GPU, and two GPUs.
   Distinguish one-fit latency from independent-fit campaign throughput.
2. Measure cold startup/JIT, cache reads, transfers, warm pressure/Jacobian time,
   gravity time, masking, reports, and full fit wall time. Synchronize CUDA events
   correctly; timing asynchronous launches alone is invalid.
3. Record peak host/device memory, PCIe traffic, reduction costs, numerical errors,
   iteration counts, and hardware utilization. Run performance tests without
   competing production jobs or explicitly disclose contention.
4. Compare identical precision and algorithm first, then report mixed-precision
   results separately. Earlier conversational 5–20x/2–8x estimates are unmeasured
   hypotheses, not acceptance evidence or promised outcomes.
5. Tune tiling, shared-memory reuse, parameter-column blocks, pinned host buffers,
   and streams only after profiling identifies the limiting work.
6. Publish installation, validated commands, input/cache identities, failure
   recovery, CPU fallback, reproducibility limits, and benchmark results.

Gate: release only backends that pass validation. Report actual speedup even if
it is modest. Keep a CPU fallback and explicit opt-in for experimental precision.

## 6. Numerical acceptance policy

C00 must write a machine-readable tolerance specification BEFORE tuning. It
must distinguish physical-model uncertainty from numerical implementation error.
Existing 0.25-degree offset and 5-degree cross-shape checks are physical/numerical
model audits, not acceptable GPU roundoff allowances.

Require exact equality for IDs, selection flags, visibility masks, profile
reasons, shape ordering, and cache identity inputs. Use both absolute and relative
vector/pressure errors to avoid division by near-zero references. Evaluate
log-pressure and objective errors directly because those drive this fit.

For gravity, set a slope-error budget much smaller than the offset audit and
resolve every threshold-sensitive facet. For fitting, determine FP64 evaluator
error and repeatability first, then bound additional mixed-precision error so it
does not change scientific conclusions or invalidate stopping diagnostics.
Document chosen numerical values, derivation, and measured maxima. There is
currently no validated mixed-precision budget to copy into this document.

Do not demand bitwise floating equality across arbitrary GPUs, compiler versions,
or block sizes. Offer deterministic reduction ordering for a pinned configuration
and measure cross-configuration differences. Preserve the distinction between
repeatability, numerical equivalence, and physical correctness.

## 7. Execution order and implementation records

Proceed C00 → C01 → C02. Gravity follows C03 → C04 → C05. Fitting follows
C06 → C07 → C08. Integrate both in C09 and characterize in C10. C06 can proceed
using accepted CPU-generated templates while the gravity port is developed.
Full FP64 fitting integration may proceed before optional C08 succeeds.

For each stage record: status, source revision, implemented files, commands,
input hashes, test report, numerical error summary, performance observations,
known limitations, and next-stage prerequisites. Keep TODOs explicit. Document
why each kernel's indexing, branch, reduction, precision, and synchronization
matches the physical formula; comments should explain decisions rather than
repeat expressions.

No package installation, port implementation, or active CPU-job modification is
performed by writing this roadmap. The first concrete implementation task is
C00 reference capture followed by C01 environment validation.

## 8. Reference sources

- CPU implementation: `../LiouvilleCOPS/tools/fit_cops.py`, particularly
  `Model.prepare`, `Model.evaluate`, `Problem`, and `fit_one`.
- Geometry/masks: `../LiouvilleCOPS/tools/build_physical_activity_templates.py`,
  `tools/gravity_cache.py`, and `tools/activity_templates.py`.
- Existing contracts: `../LiouvilleCOPS/HARMONIC_BASIS.md`,
  `ACTIVITY_ENHANCEMENT_FILE.md`, `A6_VALIDATION.md`, `A7_PERFORMANCE.md`,
  `README.md`, and `amendement_development_plan.md`.
- [CuPy 13.3 installation compatibility](https://docs.cupy.dev/en/v13.3.0/install.html).
- [CuPy custom CUDA kernels](https://docs.cupy.dev/en/stable/user_guide/kernel.html).
- [ESA polyhedral gravity documentation](https://esa.github.io/polyhedral-gravity-model/).

Consult the exact pinned library source and its license when translating the
gravity algorithm. These references describe APIs/formulations; they do not
establish a measured speedup or validate this application's future CUDA backend.
