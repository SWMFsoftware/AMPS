# liouville_cops_cuda: the fitting base built on the GPU (ROADMAP C11)

`liouville_cops_cuda` is the **unchanged** LiouvilleCOPS forward model (`../../LiouvilleCOPS/src/*.cpp`,
same compiler flags) with its two expensive loops executed on the GPU:

| Replaced (GPU) | Original (CPU, set aside by `-D` renaming in the Makefile) |
|---|---|
| `cLiouville::EvaluateKernels`: transport kernels Kn, KPhi, KV of every facet | `Liouville.cpp` → `EvaluateKernelsCPU` |
| `cSourceModel::Evaluate`: cast shadows (one ray per sunlit facet) | `Source.cpp` → `EvaluateCPU` |

Everything else is the original C++ code:
* SPICE geometry, input parsing and validation;
* the NG shield-table Monte Carlo, the source model and `ApplySource`;
* every writer (`elements.bin`, `COPS_model.dat`, `kernels/point=NNNN.bin`);
* the read-back output audit and `run_manifest.json`.

The manifest therefore records the CUDA binary as the executable, so it names the build
truthfully.

## Build and run

```tcsh
make -C native                                    # CSPICE=... CUDA=... ARCH=sm_86 to override
cd ../LiouvilleCOPS
env COPS_CUDA_DEVICE=0 ../LiouvilleCOPS_cuda/native/liouville_cops_cuda selftest.in   # analytic self-test
# a fitting base, exactly as the documented C++ rebuild, pointed at the CUDA executable:
env COPS_OUTPUT_ROOT=... COPS_CUDA_DEVICE=0 python3 tools/rebuild_cops_boundary.py --input fit_base.in \
    --surface-bc flat --shape-model <BDF> --shape-model-format nastran --out <new base> --threads 16 \
    --executable ../LiouvilleCOPS_cuda/native/liouville_cops_cuda
```

* `COPS_CUDA_DEVICE` selects the GPU (default 0).
* `COPS_CUDA_TIMING=1` prints cumulative stage times every 100 points.
* Input files are unchanged: an unknown input key would be an error in the reference parser,
  so the device is chosen through the environment.

## Files

| File | Role |
|---|---|
| `gpu_liouville.cu` | GPU implementations; device translations of `BVH.cpp` (FirstHit, AnyHit, RayTriangle, RayBox), `Mesh.cpp` (PointTriangleDistance), `Gauges.cpp` (NGWeight, RGAcceptance, ShieldFraction) and the sample loop of `EvaluateKernels` |
| `access_shim.cpp/.h` | Reads two private members of the unchanged classes: the BVH nodes and face order main() built, and the NG shield table. Cross-checked at start-up against the public methods |
| `Makefile` | Builds from `../../LiouvilleCOPS/src` with the reference flags; the two renames |

## Exactness

* **FP64 throughout,** compiled with `--fmad=false`; expressions are written in the operand
  order of the C++ code.
* **The BVH is the very tree main() built,** traversed in the same order.
* **Samples** are enumerated in the C++ order. A facet's samples are summed in chunks of
  32, and the chunk sums are added in order. Facets with at most 32 samples, i.e. every
  facet farther than a few tens of metres, are summed exactly in the C++ order.
* **Not bitwise:** the reference binary is built with `-O3 -march=native`, and GCC contracts
  multiply-adds into FMA (objdump: 40 in `EvaluateKernels`, 96 in the BVH routines).
  Results agree to rounding level; rare near-threshold visibility decisions can differ.
* **Start-up self-checks** (the run aborts on failure):
  * GPU `FirstHit` against `cBVH::FirstHit` on 4,096 random rays: no face mismatch;
  * GPU `NGWeight` against `cGaugeModel::NGWeight` on 20,000 directions: ≤ 1e-12.
* **Supported:**
  * SurfaceBC `flat`/`rough`;
  * NGAcceptance `srccg_shield`/`hemisphere`/`full`;
  * RGAcceptance `hemisphere`;
  * `FlightTimeRotation=none` (version-1 kernels).

  `ik_fov` is refused with an error; the C++ binary handles it.

## Evidence (788k SHAP5 base, 2026-10-20)

`tools/c11_compare_bases.py`, report `cops_output/cuda_c11/compare_20261020T063859Z/report.json`:
* `elements.bin`: bitwise identical.
* Records: 245,515,533 in both. One facet in each direction differs (2 points), with a share
  of the point totals ≤ 4.5e-13.
* Common records: CosSZA bitwise; 0 shadow mismatches.
* Kn/KPhi/KV: point sums agree to ≤ 1.1e-14 relative; the worst single record differs by
  9.7e-9 (a partially visible subdivided facet).
* Headers: geometry bitwise. The altitude (a minimum point–triangle distance) differs at 150
  points by ≤ 1.0e-15 relative (FMA).
* `COPS_model.dat`: pressures, densities and fluxes identical at the printed precision.

Timing (one A6000, 16 host threads), with the post-run C++ audit included in both:

| Build | Time |
|---|---|
| C++ `liouville_cops` | 292.5 s |
| `liouville_cops_cuda` | 227.5 s |

* In the GPU build, the GPU work is about 40 s (rays about 17 s, shadows about 25 s).
* The rest is unchanged C++ host work that both builds share: the read-back output audit
  (`std::set` of every record, SHA-256 of 9.8 GB), `ApplySource`, per-record `fwrite`, setup.
