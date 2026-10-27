#!/usr/bin/env python3
"""C01 CUDA environment doctor: prove that every GPU actually executes FP64 work.

Run with the CUDA environment interpreter (``.venv/bin/python``).  The checks
execute real kernels; the presence of a toolkit, a driver or shared libraries
is recorded but never counted as a pass.

Per device
  fp64_kernel      NVRTC-compiled RawKernel (no fast-math): IEEE-exact
                   +,*,/,fma,sqrt and <= 2 ulp log/exp/atan/atan2 versus NumPy
  transfers        bit-exact host->device->host round trips (pageable, pinned)
  async            stream/event ordering: work is pending before synchronize,
                   complete after it, and event timing is positive
  errors           invalid launch, compile error, and a device fault isolated
                   in a spawned process (the fault is sticky for its context)
  reduction        deterministic segmented (CSR) FP64 sum against correctly
                   rounded math.fsum, bounded by k*eps*sum|terms|, and bitwise
                   repeatable
  sparse           cuSPARSE CSR SpMV/SpMM against SciPy and run-to-run equality
  memory           pool allocate/release returns used bytes to zero
Across devices
  peer             peer-access capability and a bit-exact peer copy (optional)
  workers          one spawned persistent worker per GPU (the C07 pattern)

Exit status is nonzero if any required check fails.  The JSON report records
devices by UUID, runtime/driver/NVRTC versions and resolved library paths.
"""

import argparse
import json
import math
import multiprocessing
import os
from pathlib import Path
import sys
import time
import traceback

import numpy as np

APP_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(APP_ROOT))
from cops_cuda import environment  # noqa: E402
from cops_cuda import provenance as pv  # noqa: E402

EPS = np.finfo(np.float64).eps
# Pairwise/tree summation of n terms has error <= ~log2(n)*eps*sum|t|; 64*eps
# covers any segment length used here with margin while still detecting an
# FP32 or dropped-term reduction (those are >= 1e-8 relative).
REDUCTION_BOUND_EPS = 64.0

FP64_SOURCE = r"""
extern "C" __global__
void fp64_ops(const double* a, const double* b, const double* c, double* out, int n)
{
    // One thread per element; each output slot is a separate operation so the
    // host can apply an exact or ulp tolerance per operation.
    int i = blockDim.x * blockIdx.x + threadIdx.x;
    if (i >= n) return;
    double x = a[i], y = b[i], z = c[i];
    double* o = out + 10 * (long long)i;
    o[0] = x + y;
    o[1] = x * y;
    o[2] = x / y;
    o[3] = fma(x, y, z);
    o[4] = sqrt(fabs(x));
    o[5] = log(fabs(x));
    o[6] = exp(z);
    o[7] = atan(x);
    o[8] = atan2(x, y);
    o[9] = acos(fmax(-1.0, fmin(1.0, z)));
}

extern "C" __global__
void spin(double* out, long long iterations)
{
    // Deliberately long dependent loop used to observe asynchronous completion.
    double v = threadIdx.x + 1.0;
    for (long long k = 0; k < iterations; ++k) v = fma(v, 0.999999999, 1e-9);
    out[blockIdx.x * blockDim.x + threadIdx.x] = v;
}

extern "C" __global__
void segmented_sum(const double* values, const long long* indptr, double* out, int rows)
{
    // One block per row; fixed-stride partial sums followed by a fixed-shape
    // shared-memory tree, so the summation order depends only on the block
    // size and the row length: bitwise reproducible for a pinned launch.
    extern __shared__ double partial[];
    int row = blockIdx.x;
    if (row >= rows) return;
    double s = 0.0;
    for (long long k = indptr[row] + threadIdx.x; k < indptr[row + 1]; k += blockDim.x)
        s += values[k];
    partial[threadIdx.x] = s;
    __syncthreads();
    for (int width = blockDim.x / 2; width > 0; width >>= 1) {
        if (threadIdx.x < width) partial[threadIdx.x] += partial[threadIdx.x + width];
        __syncthreads();
    }
    if (threadIdx.x == 0) out[row] = partial[0];
}
"""

FAULT_SOURCE = r"""
extern "C" __global__ void fault(double* p) { p[(long long)1 << 40] = 1.0; }
"""


def ulp_distance(a, b):
    """Units in the last place between float64 arrays (inf for sign/NaN mismatch)."""
    a = np.asarray(a, dtype=np.float64)
    b = np.asarray(b, dtype=np.float64)
    ia = a.view(np.int64).astype(object)
    ib = b.view(np.int64).astype(object)
    # Map the sign-magnitude float layout onto a monotonic integer line.
    ia = np.where(ia < 0, -(ia & 0x7FFFFFFFFFFFFFFF), ia)
    ib = np.where(ib < 0, -(ib & 0x7FFFFFFFFFFFFFFF), ib)
    distance = np.abs((ia - ib)).astype(np.float64)
    distance[np.isnan(a) != np.isnan(b)] = np.inf
    distance[np.isnan(a) & np.isnan(b)] = 0.0
    return distance


def module(cp):
    return cp.RawModule(code=FP64_SOURCE, options=("--std=c++14",))


def check_fp64_kernel(cp, device):
    rng = np.random.default_rng(1234 + device)
    n = 1 << 16
    a = rng.uniform(-1e3, 1e3, n)
    b = rng.uniform(0.5, 2e3, n) * rng.choice((-1.0, 1.0), n)
    c = rng.uniform(-5.0, 5.0, n)
    kernel = module(cp).get_function("fp64_ops")
    out = cp.empty(10 * n)
    kernel(((n + 255) // 256,), (256,), (cp.asarray(a), cp.asarray(b), cp.asarray(c), out, np.int32(n)))
    got = cp.asnumpy(out).reshape(n, 10)
    with np.errstate(all="ignore"):
        expected = np.column_stack((a + b, a * b, a / b, np.array([math.fma(x, y, z) if hasattr(math, "fma") else np.nan
                                                                for x, y, z in zip(a, b, c)]),
                                    np.sqrt(np.abs(a)), np.log(np.abs(a)), np.exp(c), np.arctan(a),
                                    np.arctan2(a, b), np.arccos(np.clip(c, -1.0, 1.0))))
    names = ("add", "mul", "div", "fma", "sqrt", "log", "exp", "atan", "atan2", "acos")
    exact = ("add", "mul", "div", "fma", "sqrt")
    report, ok = {}, True
    for k, name in enumerate(names):
        if name == "fma" and not hasattr(math, "fma"):
            # Python < 3.13 has no math.fma: use exact rational arithmetic.
            from fractions import Fraction
            sample = range(0, n, 64)
            ref = np.array([float(Fraction(a[i]) * Fraction(b[i]) + Fraction(c[i])) for i in sample])
            distance = ulp_distance(got[list(sample), k], ref)
        else:
            distance = ulp_distance(got[:, k], expected[:, k])
        limit = 0.0 if name in exact else 2.0
        report[name] = {"max_ulp": float(np.max(distance)), "limit_ulp": limit}
        ok &= bool(np.max(distance) <= limit)
    return ok, report


def check_transfers(cp, device):
    rng = np.random.default_rng(device)
    host = rng.standard_normal(32 * 1024 * 1024)       # 256 MiB of float64
    t0 = time.perf_counter()
    dev = cp.asarray(host)
    back = cp.asnumpy(dev)
    cp.cuda.Device().synchronize()
    pageable = time.perf_counter() - t0
    pinned_mem = cp.cuda.alloc_pinned_memory(host.nbytes)
    pinned = np.frombuffer(pinned_mem, dtype=np.float64, count=host.size)
    pinned[:] = host
    stream = cp.cuda.Stream(non_blocking=True)
    start, stop = cp.cuda.Event(), cp.cuda.Event()
    with stream:
        start.record()
        dev2 = cp.empty_like(dev)
        dev2.set(pinned, stream=stream)
        stop.record()
    stop.synchronize()
    h2d_ms = cp.cuda.get_elapsed_time(start, stop)
    ok = bool(np.array_equal(back.view(np.int64), host.view(np.int64)) and
              np.array_equal(cp.asnumpy(dev2).view(np.int64), host.view(np.int64)))
    return ok, {"bytes": int(host.nbytes), "pageable_roundtrip_s": pageable,
                "pinned_h2d_GBps": host.nbytes / (h2d_ms * 1e-3) / 1e9}


def check_async(cp, device):
    spin = module(cp).get_function("spin")
    out = cp.zeros(80 * 256)
    stream = cp.cuda.Stream(non_blocking=True)
    start, stop = cp.cuda.Event(), cp.cuda.Event()
    with stream:
        start.record()
        spin((80,), (256,), (out, np.int64(20_000_000)))
        stop.record()
    pending = not stop.done
    stop.synchronize()
    elapsed = cp.cuda.get_elapsed_time(start, stop)
    finished = stop.done and bool(cp.all(out > 0).get())
    return bool(pending and finished and elapsed > 0.0), {
        "pending_before_synchronize": pending, "complete_after": finished,
        "kernel_ms": elapsed}


def _fault_child(device, queue):
    """Run in a spawned process: trigger an illegal address and report the error."""
    try:
        import cupy as cp
        cp.cuda.Device(device).use()
        kernel = cp.RawKernel(FAULT_SOURCE, "fault")
        buffer = cp.zeros(1)
        kernel((1,), (1,), (buffer,))
        cp.cuda.Device(device).synchronize()
        queue.put(("no_error", ""))
    except Exception as exc:  # the expected path
        queue.put((type(exc).__name__, str(exc)[:300]))
    # CuPy's module destructors re-raise the sticky fault while this process
    # exits; that output is expected and would only obscure the doctor report.
    sys.stderr.flush()
    os.dup2(os.open(os.devnull, os.O_WRONLY), 2)


def check_errors(cp, device):
    report = {}
    kernel = module(cp).get_function("spin")
    try:
        kernel((1,), (4096,), (cp.zeros(4096), np.int64(1)))   # > 1024 threads/block
        cp.cuda.Device().synchronize()
        report["invalid_launch"] = "no_error"
    except Exception as exc:
        report["invalid_launch"] = type(exc).__name__ + ": " + str(exc)[:200]
    try:
        cp.RawModule(code="extern \"C\" __global__ void broken( { }").get_function("broken")
        report["compile_error"] = "no_error"
    except Exception as exc:
        report["compile_error"] = type(exc).__name__
    context = multiprocessing.get_context("spawn")
    queue = context.Queue()
    child = context.Process(target=_fault_child, args=(device, queue))
    child.start()
    child.join(120)
    kind, message = queue.get(timeout=5) if not queue.empty() else ("no_report", "")
    report["device_fault_in_child"] = {"exception": kind, "message": message,
                                       "exitcode": child.exitcode}
    # The parent's context must be unaffected by the child's sticky error.
    healthy, _ = check_fp64_kernel(cp, device)
    report["parent_healthy_after_child_fault"] = healthy
    ok = (report["invalid_launch"] != "no_error" and report["compile_error"] != "no_error" and
          kind not in ("no_error", "no_report") and "illegal" in message.lower() and healthy)
    return ok, report


def check_reduction(cp, device):
    rng = np.random.default_rng(42 + device)
    lengths = np.concatenate(([0, 1, 2, 255, 256, 257], rng.integers(1, 400_000, 58)))
    indptr = np.concatenate(([0], np.cumsum(lengths))).astype(np.int64)
    # Mixed-sign terms spanning 12 decades: the RG velocity term is signed,
    # so cancellation must be bounded by the sum of magnitudes, not the sum.
    values = rng.standard_normal(indptr[-1]) * 10.0 ** rng.uniform(-6, 6, indptr[-1])
    kernel = module(cp).get_function("segmented_sum")
    d_values, d_indptr = cp.asarray(values), cp.asarray(indptr)
    results = []
    for _ in range(2):
        out = cp.zeros(len(lengths))
        kernel((len(lengths),), (256,), (d_values, d_indptr, out, np.int32(len(lengths))),
               shared_mem=256 * 8)
        results.append(cp.asnumpy(out))
    exact = np.array([math.fsum(values[indptr[k]:indptr[k + 1]].tolist()) for k in range(len(lengths))])
    magnitude = np.array([math.fsum(np.abs(values[indptr[k]:indptr[k + 1]]).tolist())
                          for k in range(len(lengths))])
    bound = REDUCTION_BOUND_EPS * EPS * magnitude
    error = np.abs(results[0] - exact)
    numpy_error = np.abs(np.add.reduceat(values, indptr[:-1])[lengths > 0] - exact[lengths > 0])
    ok = bool(np.all(error <= bound) and np.array_equal(results[0].view(np.int64), results[1].view(np.int64)))
    return ok, {"segments": int(len(lengths)), "terms": int(indptr[-1]),
                "max_error_over_eps_sum_abs": float(np.max(error / np.maximum(EPS * magnitude, 1e-300))),
                "numpy_reduceat_max_error_over_eps_sum_abs": float(np.max(
                    numpy_error / np.maximum(EPS * magnitude[lengths > 0], 1e-300))),
                "bound_eps_multiple": REDUCTION_BOUND_EPS,
                "bitwise_repeatable": bool(np.array_equal(results[0], results[1]))}


def check_sparse(cp, device):
    import scipy.sparse as sp
    import cupyx.scipy.sparse as csp
    rng = np.random.default_rng(7 + device)
    rows, cols, per_row = 2000, 200_000, 1000
    # Built directly in CSR form: scipy.sparse.random draws unique flat
    # indices from rows*cols and is needlessly slow at this size.
    indices = np.sort(rng.integers(0, cols, (rows, per_row)), axis=1).ravel()
    indptr = np.arange(0, rows * per_row + 1, per_row)
    matrix = sp.csr_matrix((rng.standard_normal(rows * per_row), indices, indptr),
                           shape=(rows, cols))
    vector = rng.standard_normal(cols)
    dense = rng.standard_normal((cols, 8))
    d_matrix = csp.csr_matrix(matrix)
    mv = [cp.asnumpy(d_matrix @ cp.asarray(vector)) for _ in range(2)]
    mm = [cp.asnumpy(d_matrix @ cp.asarray(dense)) for _ in range(2)]
    absolute = abs(matrix)
    bound_mv = REDUCTION_BOUND_EPS * EPS * (absolute @ np.abs(vector))
    bound_mm = REDUCTION_BOUND_EPS * EPS * (absolute @ np.abs(dense))
    ok = bool(np.all(np.abs(mv[0] - matrix @ vector) <= bound_mv) and
              np.all(np.abs(mm[0] - matrix @ dense) <= bound_mm))
    return ok, {"nnz": int(matrix.nnz),
                "spmv_bitwise_repeatable": bool(np.array_equal(mv[0], mv[1])),
                "spmm_bitwise_repeatable": bool(np.array_equal(mm[0], mm[1])),
                "note": "cuSPARSE repeatability is informational; C06 uses its own "
                        "fixed-order reductions where bitwise repeatability is required"}


def check_memory(cp, device):
    pool = cp.get_default_memory_pool()
    pool.free_all_blocks()
    before = pool.used_bytes()
    block = cp.empty(128 * 1024 * 1024)                 # 1 GiB
    during = pool.used_bytes()
    del block
    pool.free_all_blocks()
    after_used, after_total = pool.used_bytes(), pool.total_bytes()
    free, total = cp.cuda.Device(device).mem_info
    ok = during - before >= 1 << 30 and after_used == before and after_total == 0
    return bool(ok), {"pool_used_during_bytes": int(during - before),
                      "pool_used_after_bytes": int(after_used),
                      "pool_total_after_bytes": int(after_total),
                      "device_free_bytes": int(free), "device_total_bytes": int(total),
                      "free_fraction": free / total}


def check_peer(cp, devices):
    if len(devices) < 2:
        return None, {"skipped": "fewer than two devices"}
    runtime = cp.cuda.runtime
    can = bool(runtime.deviceCanAccessPeer(0, 1)) and bool(runtime.deviceCanAccessPeer(1, 0))
    report = {"can_access_peer_0_1_and_1_0": can}
    if not can:
        return None, report
    host = np.random.default_rng(3).standard_normal(64 * 1024 * 1024)   # 512 MiB
    for here, there in ((0, 1), (1, 0)):
        with cp.cuda.Device(here):
            try:
                runtime.deviceEnablePeerAccess(there)
            except cp.cuda.runtime.CUDARuntimeError as exc:
                if "PeerAccessAlreadyEnabled" not in str(exc):
                    raise
    with cp.cuda.Device(0):
        source = cp.asarray(host)
    with cp.cuda.Device(1):
        target = cp.empty_like(source)
        start, stop = cp.cuda.Event(), cp.cuda.Event()
        runtime.memcpyPeer(target.data.ptr, 1, source.data.ptr, 0, source.nbytes)   # warm-up
        start.record()
        runtime.memcpyPeer(target.data.ptr, 1, source.data.ptr, 0, source.nbytes)
        stop.record()
        stop.synchronize()
        report["peer_copy_GBps"] = source.nbytes / (cp.cuda.get_elapsed_time(start, stop) * 1e-3) / 1e9
        report["peer_copy_bit_exact"] = bool(np.array_equal(cp.asnumpy(target).view(np.int64),
                                                            host.view(np.int64)))
    return report["peer_copy_bit_exact"], report


def _worker(device, queue):
    """Persistent-worker pattern: own device, own data, return only small results."""
    os.environ.update(environment.worker_thread_environment(1))
    try:
        import cupy as cp
        cp.cuda.Device(device).use()
        data = cp.arange(10_000_000, dtype=cp.float64)
        queue.put((device, float(cp.sum(data).get()), environment.device_uuid(device)))
    except Exception as exc:
        queue.put((device, None, repr(exc)))


def check_workers(devices):
    context = multiprocessing.get_context("spawn")
    queue = context.Queue()
    processes = [context.Process(target=_worker, args=(d["index"], queue)) for d in devices]
    for p in processes:
        p.start()
    results = sorted(queue.get(timeout=300) for _ in processes)
    for p in processes:
        p.join(60)
    expected = 10_000_000 * (10_000_000 - 1) / 2.0
    ok = all(r[1] == expected for r in results) and all(p.exitcode == 0 for p in processes)
    return bool(ok), {"results": [{"device": r[0], "sum": r[1], "uuid_or_error": r[2]}
                                  for r in results]}


def run_check(name, function, *args):
    t0 = time.perf_counter()
    try:
        ok, detail = function(*args)
    except Exception:
        ok, detail = False, {"exception": traceback.format_exc(limit=4)}
    return {"name": name, "passed": ok, "seconds": time.perf_counter() - t0, "detail": detail}


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--json", help="report path (default: new file under "
                        "$COPS_OUTPUT_ROOT/cuda_doctor/)")
    parser.add_argument("--devices", help="comma-separated device indices (default: all)")
    args = parser.parse_args()

    distributions = environment.cupy_distributions()
    import cupy as cp
    report = {"schema": "liouville_cops_cuda.doctor", "schema_version": 1,
              "process": pv.process_context(), "runtime": environment.runtime_versions(),
              "devices": environment.device_inventory(), "checks": []}
    report["loaded_libraries"] = environment.loaded_cuda_libraries()
    report["checks"].append({"name": "single_cupy_distribution",
                             "passed": len(distributions) == 1,
                             "detail": distributions, "seconds": 0.0})
    selected = [d for d in report["devices"]
                if args.devices is None or str(d["index"]) in args.devices.split(",")]
    for device in selected:
        with cp.cuda.Device(device["index"]):
            for name, function in (("fp64_kernel", check_fp64_kernel),
                                   ("transfers", check_transfers), ("async", check_async),
                                   ("errors", check_errors), ("reduction", check_reduction),
                                   ("sparse", check_sparse), ("memory", check_memory)):
                result = run_check(name, function, cp, device["index"])
                result["device"] = device["uuid"]
                report["checks"].append(result)
                print("%-8s %-12s %s (%.2f s)" % (device["uuid"][:12], name,
                      "PASS" if result["passed"] else "FAIL", result["seconds"]), flush=True)
    peer = run_check("peer", check_peer, cp, selected)
    peer["required"] = False
    report["checks"].append(peer)
    print("%-21s %s" % ("peer", {True: "PASS", False: "FAIL", None: "N/A"}[peer["passed"]]))
    workers = run_check("workers", check_workers, selected)
    report["checks"].append(workers)
    print("%-21s %s" % ("workers", "PASS" if workers["passed"] else "FAIL"))
    report["loaded_libraries_after"] = environment.loaded_cuda_libraries()
    required = [c for c in report["checks"] if c.get("required", True)]
    report["passed"] = all(c["passed"] for c in required)

    if args.json:
        path = Path(args.json)
    else:
        root = Path(os.environ.get("COPS_OUTPUT_ROOT") or APP_ROOT.parents[1] / "cops_output")
        path = root / "cuda_doctor" / ("doctor_%s.json" % pv.utc_stamp())
    path.parent.mkdir(parents=True, exist_ok=True)
    pv.write_json_exclusive(path, report)
    print("report: %s\noverall: %s" % (path, "PASS" if report["passed"] else "FAIL"))
    return 0 if report["passed"] else 1


if __name__ == "__main__":
    sys.exit(main())
