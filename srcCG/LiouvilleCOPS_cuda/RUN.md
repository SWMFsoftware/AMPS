# Running the LiouvilleCOPS CUDA development tools

Two interpreters are used on purpose:

| Interpreter | Use | Packages |
|---|---|---|
| `/usr/bin/python3` (3.8.10) | **CPU reference**: C00 capture only. It is the interpreter that produced the campaign. | NumPy 1.24.4, SciPy 1.10.1 (user site), matplotlib 3.1.2; `polyhedral_gravity` 3.3 from `PYTHONPATH=/tmp/liouville_cops_polyhedral_gravity` |
| `.venv/bin/python` (3.11.14) | **CUDA environment**: doctor, tests, and later the GPU backend | `requirements-cuda.lock`: CuPy 13.6.0, NumPy 1.24.4, SciPy 1.10.1, polyhedral-gravity 3.3, pytest |

The shell is tcsh; the commands below use `setenv`.

## Environment (C01)

Create or recreate the CUDA environment. This never touches the system Python:

```tcsh
cd ~/ROS-END-OF-MISSION/AMPS/srcCG/LiouvilleCOPS_cuda
uv python install 3.11
uv venv --python 3.11 .venv
uv pip sync --python .venv/bin/python --require-hashes requirements-cuda.lock
```

To change a dependency, edit `requirements-cuda.in`, then run
`uv pip compile requirements-cuda.in --python-version 3.11 --generate-hashes -o requirements-cuda.lock`,
sync, and rerun the doctor and the tests.

CuPy uses the system CUDA 12.4 toolkit (`/usr/local/cuda-12.4`: NVRTC, cuBLAS, cuSPARSE).
Do **not** install `cupy-cuda12x[ctk]` or `nvidia-*-cu12` wheels. A newer pip NVRTC can emit
PTX that the 550 driver (CUDA 12.4) refuses. CuPy's own linked runtime is 12.9, which runs
on the 12.4 driver through CUDA minor-version compatibility. The doctor verifies that by
executing kernels.

## Device check

```tcsh
setenv COPS_OUTPUT_ROOT /data/vtenishe/rosetta_data_end_of_mission_cops/cops_output
.venv/bin/python tools/cuda_doctor.py            # all GPUs; writes $COPS_OUTPUT_ROOT/cuda_doctor/doctor_<UTC>.json
.venv/bin/python -m pytest -q tests              # provenance, synthetic shapes, per-device checks
```

The exit status is nonzero if a required check fails. Peer access is recorded but optional.

## C02 device data: evidence on the campaign run

```tcsh
setenv COPS_OUTPUT_ROOT /data/vtenishe/rosetta_data_end_of_mission_cops/cops_output
.venv/bin/python tools/c02_device_load.py \
    --capture $COPS_OUTPUT_ROOT/cuda_reference/c00_20261019T064536Z --full
```

The tool writes `$COPS_OUTPUT_ROOT/cuda_c02/load_<UTC>/report.json` and, with `--full`,
`backend_manifest.json`. It loads the 47 C00 fit points and checks them bitwise against
the CPU `fit_cops.read_kernel`. It then runs three load/release cycles. With `--full` it
loads all 1150 points (9.1 GB on the device), verifies them against the manifest's set
digest and runs a round trip. Finally it refuses a plan one byte too large and repeats the
subset on the second GPU. `--device N` selects the GPU and `--cycles N` the cycle count.
Reading the 9.8 GB of kernel files dominates the run time. The tool only reads the run
directory, single-threaded, so it can run beside CPU jobs.

The unit tests in `tests/test_formats_memory.py`, `tests/test_device_store.py` and
`tests/test_backend_contract.py` use small synthetic runs (`tests/synthetic_run.py`) and
need no campaign data. `test_device_store.py` needs a GPU; the others also pass under the
CPU reference interpreter's NumPy.

In code (CUDA environment):

```python
from cops_cuda import device_store
with device_store.DeviceDataset.load(run_dir, rows=model_rows, device_id=0) as ds:
    print(ds.summary()); ds.verify_roundtrip()
```

## C03 gravity: evidence on the 788k body

```tcsh
setenv COPS_OUTPUT_ROOT /data/vtenishe/rosetta_data_end_of_mission_cops/cops_output
.venv/bin/python -u tools/c03_gravity_check.py \
    --capture $COPS_OUTPUT_ROOT/cuda_reference/c00_20261019T064536Z --cpus 28-35 [--full]
```

The tool writes `$COPS_OUTPUT_ROOT/cuda_c03/gravity_<UTC>/report.json`, plus
`full_field.npz` with `--full`. It compares against:
* the C00 fixtures (campaign binary), synthetic shapes and the 832 + 32 targets on the 788k
  body;
* the package state, bit for bit;
* new difficult-location samples evaluated by the cp311 package on the CPUs given by
  `--cpus`.

It also checks determinism, including on the second GPU, and measures timing. `--full`
evaluates all 787,754 targets on one GPU, about 45 min. When the package evaluates
targets far from the body it logs one `isCriticalDifference` warning per face; this is
harmless, so filter its log lines (`| grep -v POLYHEDRAL_GRAVITY_LOGGER`). Unit tests:
`tests/test_gravity.py`.

## C04 full-field gravity equivalence (100k test body)

The user decided on 2026-10-19 to test the gravity model at full-field scale on the smaller
`SHAP5_Volume.bdf`. Its 100k-triangle surface is wound inward, so derive an outward copy
once (exact reversal of every face, with a provenance record):

```tcsh
.venv/bin/python tools/derive_outward_surface.py \
    /data/vtenishe/rosetta_data_end_of_mission_cops/SHAP5_Volume.bdf \
    /data/vtenishe/rosetta_data_end_of_mission_cops/shape_models/derived/SHAP5_Volume_surface_outward.bdf
```

**CPU reference.** Run the unchanged repaired builder from `../LiouvilleCOPS` with this
body. Pin it to a few CPUs when other CPU work is running; TBB follows the affinity mask:

```tcsh
cd ../LiouvilleCOPS
setenv PYTHONPATH /data/vtenishe/rosetta_data_end_of_mission_cops/liouville_cops_polyhedral_gravity
taskset -c 28-35 /usr/bin/python3 tools/build_physical_activity_templates.py \
    $COPS_OUTPUT_ROOT/output_fit_shap5_788k_base \
    --out $COPS_OUTPUT_ROOT/physical_templates_shap5_788k_targets_shap5vol100k_gravity_test.npz \
    --late-start 2016-09-30T10:33:10 --late-stop 2016-09-30T10:39:26 \
    --gravity-shape .../shape_models/derived/SHAP5_Volume_surface_outward.bdf --density 533 \
    --cliff-threshold 45 --cliff-threshold 60 --gravity-offset 0.1 --gravity-offset-tolerance 0.25 \
    --omega-executable ./spice_angular_velocity
```

**GPU field and facet-by-facet comparison:**

```tcsh
.venv/bin/python tools/c04_gravity_equivalence.py gpu --run $COPS_OUTPUT_ROOT/output_fit_shap5_788k_base \
    --gravity-shape .../shape_models/derived/SHAP5_Volume_surface_outward.bdf      # about 7 min on one A6000
.venv/bin/python tools/c04_gravity_equivalence.py compare --dir $COPS_OUTPUT_ROOT/cuda_c04/equivalence_<UTC> \
    --cpu-cache $COPS_OUTPUT_ROOT/physical_templates_shap5_788k_targets_shap5vol100k_gravity_test.gravity.npz \
    --cpu-bundle $COPS_OUTPUT_ROOT/physical_templates_shap5_788k_targets_shap5vol100k_gravity_test.npz \
    --gpu-788k-field $COPS_OUTPUT_ROOT/cuda_c04/equivalence_20261020T021855Z/gpu_field.npz
```

## C05 visibility masks and the C06 evaluator check

```tcsh
.venv/bin/python tools/c05_mask_check.py --capture $COPS_OUTPUT_ROOT/cuda_reference/c00_20261019T064536Z
.venv/bin/python tools/c06_evaluator_check.py --capture $COPS_OUTPUT_ROOT/cuda_reference/c00_20261019T064536Z
```

The C00b fixtures come from `/usr/bin/python3 tools/c00_capture_reference.py joint --capture <C00 capture>`,
run under the CPU reference interpreter (about 20 min on 4 CPUs).

## C04 GPU cache for reuse by the unchanged CPU builder

```tcsh
.venv/bin/python tools/c04_write_gpu_cache.py --field-dir <equivalence dir or C03 report dir> \
    --run $COPS_OUTPUT_ROOT/output_fit_shap5_788k_base --gravity-shape <BDF> --out <new cache.npz> \
    [--evidence <compare report.json>]
# then, in ../LiouvilleCOPS: build_physical_activity_templates.py ... --gravity-cache <new cache.npz>
```

## C05 restart bundle (gravity cache + NG/RG/cliff masks)

```tcsh
.venv/bin/python tools/c05_write_bundle.py --run $COPS_OUTPUT_ROOT/output_fit_shap5_788k_base \
    --gravity-cache <admitted cache.npz> --gravity-shape <the cache's BDF> --out <new bundle.npz>
# columns: late_ng_not_rg late_rg_not_ng cliff_ge_45deg cliff_ge_60deg;
# the unchanged fitter selects each with: fit_cops.py ... --activity-template <bundle.npz> <name>
```

Facets within 1e-6° of a cliff threshold are re-resolved with the CPU reference
package (requires `polyhedral_gravity` in `.venv`); the bundle metadata lists them.

## C06 full-set evaluator run (all 1150 fit points)

```tcsh
env COPS_OUTPUT_ROOT=/data/vtenishe/rosetta_data_end_of_mission_cops/cops_output \
    .venv/bin/python tools/c06_full_run.py --capture $COPS_OUTPUT_ROOT/cuda_reference/c00_20261019T064536Z
```

The full problem uses the header-only model (`cops_cuda/host_model.py`), so host
memory stays small (kernel records live only on the GPU, about 22 GB on device).
A 23-point CPU subset is the exact reference for those points of the full problem.
`COPS_OUTPUT_ROOT` must be set: the reference path helpers resolve the base run with it.

## GPU fit with the unchanged fitter (step 1 of C07)

```tcsh
env COPS_OUTPUT_ROOT=/data/vtenishe/rosetta_data_end_of_mission_cops/cops_output \
    OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
    .venv/bin/python tools/gpu_fit.py --device 1 -- \
    $COPS_OUTPUT_ROOT/output_fit_shap5_788k_base --profile all_available \
    --out $COPS_OUTPUT_ROOT/<new fit directory> --degree 4 --enhancement-degree 4 --lambda 0 \
    --jobs 1 --tbounds 30 400 --max-nfev 1000 --ftol 1e-10 --xtol 1e-10 --gtol 1e-10 \
    --check-tolerance 1e-6 [--activity-template <bundle.npz> <name> ...]
```

* Everything after `--` is passed to `fit_cops.main()` unchanged, so the outputs are the
  fitter's own: `fit_manifest.json`, `fit_result.txt`, `ActivityFit.txt`, `cops_fit.in`, and so on.
* The manifest has an extra `compute_backend` section that names the CUDA backend.
* Use one `--lambda` value per process; an L-curve is refused because the fitter forks.
* `--devices 0,1` splits one fit over both GPUs (bitwise equal to one GPU; `cops_cuda/multi_gpu.py`).
* `--mode sequential` gives the bitwise CPU summation order. It is slow and meant only for
  verification.

## C07 multi-GPU gate and benchmark

```tcsh
env COPS_OUTPUT_ROOT=/data/vtenishe/rosetta_data_end_of_mission_cops/cops_output \
    .venv/bin/python tools/c07_multi_gpu_check.py --capture $COPS_OUTPUT_ROOT/cuda_reference/c00_20261019T064536Z
```

Both GPUs must be free of other fits for the timing to mean anything.

## Threads

Run CPU work that runs alone (C++ replay, CPU gravity reference) with 36 threads. A GPU
worker sets `OMP_NUM_THREADS=OPENBLAS_NUM_THREADS=MKL_NUM_THREADS=1`
(`cops_cuda.environment.worker_thread_environment`), so N workers do not start N
36-thread pools. `polyhedral_gravity` uses TBB, which **ignores** `OMP_NUM_THREADS`. Limit it
with the CPU affinity mask instead: `taskset -c 30-35 …` or the C00 `--cpus 30-35` option.
Otherwise it uses all 36 CPUs.

## C00 reference capture

```tcsh
setenv COPS_OUTPUT_ROOT /data/vtenishe/rosetta_data_end_of_mission_cops/cops_output
setenv PYTHONPATH /tmp/liouville_cops_polyhedral_gravity
set C = `/usr/bin/python3 tools/c00_capture_reference.py new`
/usr/bin/python3 tools/c00_capture_reference.py identities --capture $C   # ~40 min on a loaded host (reads 20 GB)
/usr/bin/python3 tools/c00_capture_reference.py masks      --capture $C
/usr/bin/python3 tools/c00_capture_reference.py fit        --capture $C   # --profile, --points, --lambda
/usr/bin/python3 tools/c00_capture_reference.py gravity    --capture $C --cpus 30-35
```

Each stage creates `$C/<stage>/` (an existing directory is refused) with a
`stage_manifest.json` that lists the digests of its products, of the reference tree and of
the capture tool. The `masks`, `fit` and `gravity` stages check every kernel they read
against `identities/kernel_files.json`.

If `/tmp/liouville_cops_polyhedral_gravity` disappears (for example after a reboot), unpack
`third_party/download_polyhedral_gravity_v3.3/*cp38*.whl` (`unzip -d DIR`) and point
`PYTHONPATH` at DIR. Its `.so` is byte-identical.

## Output locations

All products go below `$COPS_OUTPUT_ROOT` in unique, timestamped directories. They are
never written inside `../LiouvilleCOPS` or this source tree.
