#!/usr/bin/env python3
"""Where does the time of a GPU fit evaluation go? (C08 decision support; CUDA env).

Hardware counters (Nsight Compute) need administrator rights on this host
(``RmProfilingAdminOnly=1``). This tool therefore separates the costs by
timing and controlled experiments; ``nsys`` can trace the same run for exact
kernel and transfer times.

On the full 788k fit set (B29 problem, tree mode, one device):

``breakdown``
    One pressure evaluation and one full Jacobian, split into: host activity
    matvec, uploads, each kernel launch (CUDA events), sum downloads, host
    assembly.
``sums_scaling``
    ``eval_partial`` time against the number of sums per launch (1–32) for each
    code family (pressure, temperature, activity column). The slope is the cost
    per sum per record; the intercept is the per-record base cost (loads and
    ``u = 1/sqrt(T)``).
``memory_baseline``
    A kernel that only reads the per-record arrays ``eval_partial`` reads, plus
    the facet gather, and reduces them. This is the floor set by memory traffic.
``fp64_peak``
    An FP64 multiply-add chain (no FMA, as ``--fmad=false``) with many independent
    accumulators, measuring the FP64 throughput this device delivers.

Writes ``$COPS_OUTPUT_ROOT/cuda_c08/profile_<UTC>/report.json``.
"""

import argparse
import json
import os
from pathlib import Path
import sys
import time

import numpy as np

APP_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(APP_ROOT))
sys.path.insert(0, str(APP_ROOT / "tools"))

DEFAULT_OUTPUT_ROOT = "/data/vtenishe/rosetta_data_end_of_mission_cops/cops_output"

MEMORY_KERNEL = r'''
extern "C" __global__ void read_only(const long long* __restrict__ blocks, const int* __restrict__ face,
    const double* __restrict__ a0, const double* __restrict__ a1, const double* __restrict__ a2,
    const double* __restrict__ a3, const double* __restrict__ a4, const double* __restrict__ a5,
    const unsigned char* __restrict__ lit, const double* __restrict__ af_face,
    const unsigned char* __restrict__ active, double* __restrict__ out)
{
    __shared__ double tree[256];
    const long long b = blockIdx.x, begin = blocks[3*b+1], end = blocks[3*b+2];
    double s = 0.0;
    for (long long r = begin + threadIdx.x; r < end; r += 256) {
        const int f = face[r];
        s += a0[r] + a1[r] + a2[r] + a3[r] + a4[r] + a5[r] + (double)lit[r] + af_face[f] + (double)active[f];
    }
    tree[threadIdx.x] = s; __syncthreads();
    for (int w = 128; w > 0; w >>= 1) { if (threadIdx.x < w) tree[threadIdx.x] += tree[threadIdx.x + w]; __syncthreads(); }
    if (threadIdx.x == 0) out[b] = tree[0];
}
extern "C" __global__ void fp64_chain(double* __restrict__ out, long long iters, double a, double b)
{
    double x[8];
    for (int k = 0; k < 8; ++k) x[k] = threadIdx.x * 1e-9 + k;
    for (long long i = 0; i < iters; ++i)
        for (int k = 0; k < 8; ++k) x[k] = x[k] * a + b;        // mul + add (no FMA with --fmad=false)
    double s = 0.0;
    for (int k = 0; k < 8; ++k) s += x[k];
    out[(long long)blockIdx.x * blockDim.x + threadIdx.x] = s;
}
'''


def log(message):
    print("[%s] %s" % (time.strftime("%H:%M:%S"), message), flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--capture", required=True)
    parser.add_argument("--device", type=int, default=0)
    parser.add_argument("--repeats", type=int, default=5)
    parser.add_argument("--output-root", default=os.environ.get("COPS_OUTPUT_ROOT", DEFAULT_OUTPUT_ROOT))
    args = parser.parse_args()

    import cupy
    from c07_multi_gpu_check import build
    from cops_cuda import backend as bk, cuda_sources, evaluator, manifest, precision, provenance
    capture = Path(args.capture).resolve()
    run, G, problems = build(capture, ["B29"])
    P = problems["B29"]
    out = provenance.make_new_directory(Path(args.output_root) / "cuda_c08" /
                                        ("profile_" + provenance.utc_stamp()))
    report = {"schema": "liouville_cops_cuda.c08_profile", "schema_version": 1,
              "created_utc": provenance.utc_now(), "points": int(G.npt),
              "records": int(G.indptr[-1]), "device": args.device}
    nrec = int(G.indptr[-1])

    with evaluator.CudaBackend(G, run, None, P.enhancement_scale, device_id=args.device) as be:
        dev = cupy.cuda.Device(args.device)
        request = bk.ColumnRequest.from_problem(P)
        state = bk.ParameterState.from_problem(P, P.x0)
        be.evaluate(state, request)                         # warm up (JIT, caches)

        # ---------------------------------------------------------------- breakdown
        def timed(fn):
            dev.synchronize()
            t = time.perf_counter()
            r = fn()
            dev.synchronize()
            return r, time.perf_counter() - t

        events = []
        original = cuda_sources.function

        def traced(module, name, device_id):
            kernel = original(module, name, device_id)

            def launch(grid, block, args_):
                start, stop = cupy.cuda.Event(), cupy.cuda.Event()
                start.record()
                kernel(grid, block, args_)
                stop.record()
                events.append((name, start, stop))
            return launch

        breakdown = {}
        for label, req in (("pressures", bk.PRESSURES_ONLY), ("full_jacobian", request)):
            samples = []
            for _ in range(args.repeats):
                events.clear()
                cuda_sources.function = traced
                try:
                    with dev:
                        _, t_act = timed(lambda: be.activity_linear(state))
                        _, t_total = timed(lambda: be.evaluate(state, req))
                finally:
                    cuda_sources.function = original
                kern = {}
                for name, a, b in events:
                    b.synchronize()
                    kern[name] = kern.get(name, 0.0) + cupy.cuda.get_elapsed_time(a, b) / 1e3
                samples.append({"total_s": t_total, "host_activity_matvec_s": t_act,
                                "kernels_s": kern, "launches": len(events)})
            med = lambda key: float(np.median([s[key] for s in samples]))
            kern_names = samples[0]["kernels_s"].keys()
            k_med = {n: float(np.median([s["kernels_s"][n] for s in samples])) for n in kern_names}
            breakdown[label] = {"total_s": med("total_s"), "host_activity_matvec_s": med("host_activity_matvec_s"),
                                "kernels_s": k_med, "kernel_sum_s": sum(k_med.values()),
                                "other_s": med("total_s") - sum(k_med.values()),
                                "launches": samples[0]["launches"]}
            log("%s: %s" % (label, breakdown[label]))
        report["breakdown"] = breakdown

        # ------------------------------------------------------------ sums scaling
        be._upload_activity(state)
        be._upload_columns(request)
        params = be._params(state)
        families = {"pressure": [0, 1, 2, 3, 4], "temperature": list(range(5, 14)),
                    "activity_column_NG": [1000 + c for c in range(25)],
                    "activity_column_RG": [2000 + c for c in range(25)]}
        scaling = {}
        for fam, base_codes in families.items():
            rows = []
            for n in (1, 2, 4, 8, 16, 32):
                codes = [base_codes[k % len(base_codes)] for k in range(n)]
                ts = []
                for _ in range(args.repeats):
                    events.clear()
                    cuda_sources.function = traced
                    try:
                        with dev:
                            be._run_sums(codes, params)
                    finally:
                        cuda_sources.function = original
                    for name, a, b in events:
                        if name == "eval_partial":
                            b.synchronize()
                            ts.append(cupy.cuda.get_elapsed_time(a, b) / 1e3)
                rows.append({"nsums": n, "eval_partial_s": float(np.median(ts))})
            x = np.array([r["nsums"] for r in rows], float)
            y = np.array([r["eval_partial_s"] for r in rows])
            slope, intercept = np.polyfit(x, y, 1)
            scaling[fam] = {"points": rows, "slope_s_per_sum": float(slope), "intercept_s": float(intercept),
                            "ns_per_record_per_sum": float(slope / nrec * 1e9),
                            "ns_per_record_base": float(intercept / nrec * 1e9)}
            log("%s: slope %.4f s/sum (%.3f ns/record/sum), intercept %.4f s" %
                (fam, slope, slope / nrec * 1e9, intercept))
        report["sums_scaling"] = scaling

        # --------------------------------------------------------- memory baseline
        with dev:
            mod = cupy.RawModule(code=MEMORY_KERNEL, options=precision.REFERENCE_NVRTC_OPTIONS)
            read_only, chain = mod.get_function("read_only"), mod.get_function("fp64_chain")
            a = be.arrays
            nblocks = int(be._first_block[-1])
            red = cupy.empty(nblocks, dtype=np.float64)
            ts = []
            for _ in range(args.repeats):
                s, e = cupy.cuda.Event(), cupy.cuda.Event()
                s.record()
                read_only((nblocks,), (256,), (a["blocks"], be.dataset.arrays["kernel.face"], a["wnH"], a["wvH"],
                                               a["wpH"], a["wnC"], a["wvC"], a["xlit"], a["lit"],
                                               a["af_face"], a["active"], red))
                e.record(); e.synchronize()
                ts.append(cupy.cuda.get_elapsed_time(s, e) / 1e3)
            t_mem = float(np.median(ts))
            bytes_per_record = 4 + 6 * 8 + 1 + 8 + 1
            report["memory_baseline"] = {"read_only_s": t_mem, "bytes_per_record": bytes_per_record,
                                         "effective_GBps": nrec * bytes_per_record / t_mem / 1e9}
            log("memory baseline: %.4f s (%.0f GB/s effective)" % (t_mem, nrec * bytes_per_record / t_mem / 1e9))

            # ------------------------------------------------------------- FP64 peak
            props = cupy.cuda.runtime.getDeviceProperties(args.device)
            nsm = props["multiProcessorCount"]
            grid, block, iters = nsm * 32, 256, 20000
            o = cupy.empty(grid * block, dtype=np.float64)
            chain((grid,), (block,), (o, np.int64(100), np.float64(0.999), np.float64(1e-3)))
            s, e = cupy.cuda.Event(), cupy.cuda.Event()
            s.record()
            chain((grid,), (block,), (o, np.int64(iters), np.float64(0.999), np.float64(1e-3)))
            e.record(); e.synchronize()
            t = cupy.cuda.get_elapsed_time(s, e) / 1e3
            flops = 2.0 * 8 * iters * grid * block
            report["fp64_peak"] = {"seconds": t, "gflops": flops / t / 1e9,
                                   "note": "mul+add chains, 8 independent per thread, --fmad=false"}
            log("FP64 throughput (mul+add, no FMA): %.0f GFLOP/s" % (flops / t / 1e9))

    report["manifest"] = manifest.backend_manifest("cuda_fp64_v1", args.device)
    provenance.write_json_exclusive(out / "report.json", report)
    log("report %s" % (out / "report.json"))
    return 0


if __name__ == "__main__":
    sys.exit(main())
