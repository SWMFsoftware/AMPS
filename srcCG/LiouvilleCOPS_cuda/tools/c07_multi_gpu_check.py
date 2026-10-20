#!/usr/bin/env python3
"""C07 gate: one vs two GPUs on the full 788k fit set, plus the throughput benchmark (CUDA env).

For B29, N30 and NC31 (templates from the C00b capture) on all 1150 fit points:

``bitwise_one_vs_two``
    ``MultiGpuBackend`` on devices 0 and 1 against ``CudaBackend`` on device 0, for every
    selected C00b fixture state and the full Jacobian request. This is the gate:
    one-device and two-device objective equivalence. Because each point is
    independent and the devices are deterministic, the expected result is bitwise
    identity of pressures and every column; the residuals and objective then agree
    exactly.
``repeat_bitwise``
    The same state evaluated twice on two devices (deterministic pinned reductions).
``timing``
    Warm medians for pressures only and the full Jacobian, on one and on two devices,
    including the host matvec and the assembly.

Benchmark (ROADMAP C07 item 2): **one fit on two GPUs versus independent fits, one
per GPU.** Two spawned processes each hold the full B29 dataset on their own
device and time the full Jacobian concurrently. Throughput is then compared with
one fit split over both devices.

Writes ``$COPS_OUTPUT_ROOT/cuda_c07/check_<UTC>/report.json``.
"""

import argparse
import json
import multiprocessing as mp
import os
from pathlib import Path
import sys
import time

import numpy as np

APP_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(APP_ROOT))
sys.path.insert(0, str(APP_ROOT / "tools"))

DEFAULT_OUTPUT_ROOT = "/data/vtenishe/rosetta_data_end_of_mission_cops/cops_output"


def log(message):
    print("[%s] %s" % (time.strftime("%H:%M:%S"), message), flush=True)


def build(capture, configs):
    """The full fit problem per configuration, on the GPU-backed header model."""
    import c00_capture_reference as c00
    from cops_cuda import gpu_model, reference
    root = reference.bind()
    modules = c00.load_reference_modules()
    fit_cops = modules["fit_cops"]
    run = c00.base_directory(modules, c00.DEFAULT_BASE)
    header, names, rows, columns, _, _, _ = c00.fit_selection(modules, run, c00.DEFAULT_PROFILE)
    coef = fit_cops.read_activity_coefficients(str(root), header)
    elements = fit_cops.read_elements(str(run / "elements.bin"))
    G = gpu_model.gpu_model_class()(run, header, columns, rows, elements, coef)
    G.load_kernels(G.sel)
    templates = np.load(str(capture / "joint" / "templates.npz"))
    library = {n: templates[n].astype(float)
               for n in ("late_ng_not_rg", "late_rg_not_ng", "syn_cliff_standin")}
    problems = {}
    for label, template_names, _ in c00.JOINT_CONFIGS:
        if label not in configs:
            continue
        config = {"name": label, "degree": 4, "enhancement_degree": 4,
                  "fix_legacy_activity": False, "fix_temperature": False, "fix_kappa": False,
                  "ordered_temperature": False, "templates": bool(template_names)}
        tpl = (list(template_names),
               np.column_stack([library[n] for n in template_names]) if template_names
               else np.empty((len(elements["area"]), 0)), {})
        problems[label] = c00.make_problem(fit_cops, G, elements, coef, config, tpl)
    return run, G, problems


def bench_worker(capture, device, repeats, barrier, queue):
    """One independent B29 'fit' process: full dataset on one device, timed Jacobians."""
    try:
        from cops_cuda import backend as bk, evaluator
        run, G, problems = build(Path(capture), ["B29"])
        P = problems["B29"]
        request = bk.ColumnRequest.from_problem(P)
        with evaluator.CudaBackend(G, run, None, P.enhancement_scale, device_id=device) as one:
            state = bk.ParameterState.from_problem(P, P.x0)
            one.evaluate(state, request)
            barrier.wait()
            times = []
            for _ in range(repeats):
                t0 = time.perf_counter()
                one.evaluate(state, request)
                times.append(time.perf_counter() - t0)
        queue.put((device, {"full_jacobian_s_median": float(np.median(times)), "samples": times}))
    except BaseException as exc:
        try:
            barrier.abort()
        except Exception:
            pass
        queue.put((device, {"error": repr(exc)}))


def same(u, v):
    return np.asarray(u).tobytes() == np.asarray(v).tobytes()


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--capture", required=True)
    parser.add_argument("--configs", default="B29,N30,NC31")
    parser.add_argument("--vectors", default="prior,legacy_strong_clip,archived_b29_state,"
                                            "enhancement_strong_clip,template_on_bound")
    parser.add_argument("--repeats", type=int, default=5)
    parser.add_argument("--skip-benchmark", action="store_true")
    parser.add_argument("--output-root", default=os.environ.get("COPS_OUTPUT_ROOT", DEFAULT_OUTPUT_ROOT))
    args = parser.parse_args()

    from cops_cuda import backend as bk, evaluator, manifest, provenance
    from cops_cuda.memory import DeviceBudget
    from cops_cuda.multi_gpu import MultiGpuBackend
    capture = Path(args.capture).resolve()
    fixtures = json.loads((capture / "joint" / "joint_fixtures.json").read_text())
    run, G, problems = build(capture, args.configs.split(","))
    out = provenance.make_new_directory(Path(args.output_root) / "cuda_c07" /
                                        ("check_" + provenance.utc_stamp()))
    log("report %s; %d points, %d records" % (out, G.npt, G.indptr[-1]))
    results, all_ok = {}, True
    budget = DeviceBudget(0)
    try:
        for label, P in problems.items():
            arrays = np.load(str(capture / "joint" / ("joint_%s.npz" % label)))
            record = fixtures["configs"][label]
            vectors = [v for v in args.vectors.split(",") if v in record["vectors"]]
            request = bk.ColumnRequest.from_problem(P)
            basis = P.enhancement_basis if P.nb else None
            t0 = time.perf_counter()
            one = evaluator.CudaBackend(G, run, basis, P.enhancement_scale, device_id=0,
                                        budget=budget)
            bind_one = time.perf_counter() - t0
            t0 = time.perf_counter()
            two = MultiGpuBackend(G._spec, G, basis, P.enhancement_scale, devices=(0, 1))
            bind_two = time.perf_counter() - t0
            try:
                per = {}
                for v in vectors:
                    state = bk.ParameterState.from_problem(P, arrays["%s__%s__x" % (label, v)])
                    a, b = one.evaluate(state, request), two.evaluate(state, request)
                    c = two.evaluate(state, request)
                    bit = (same(a.pNG, b.pNG) and same(a.pRG, b.pRG) and
                           all(same(a.columns[k][i], b.columns[k][i])
                               for k in request.keys() for i in (0, 1)))
                    rep = (same(b.pNG, c.pNG) and same(b.pRG, c.pRG) and
                           all(same(b.columns[k][i], c.columns[k][i])
                               for k in request.keys() for i in (0, 1)))
                    per[v] = {"bitwise_one_vs_two": bool(bit), "repeat_bitwise": bool(rep)}
                    all_ok = all_ok and bit and rep
                state = bk.ParameterState.from_problem(P, P.x0)
                timing = {"bind_one_s": bind_one, "bind_two_s": bind_two}
                for name, backend in (("one", one), ("two", two)):
                    backend.evaluate(state, request)
                    tp, tj = [], []
                    for _ in range(args.repeats):
                        t = time.perf_counter(); backend.evaluate(state); tp.append(time.perf_counter() - t)
                        t = time.perf_counter(); backend.evaluate(state, request); tj.append(time.perf_counter() - t)
                    timing["%s_pressures_s_median" % name] = float(np.median(tp))
                    timing["%s_full_jacobian_s_median" % name] = float(np.median(tj))
                timing["two_worker_s_last"] = two.timing.get("worker_s")
                timing["speedup_full_jacobian"] = (timing["one_full_jacobian_s_median"] /
                                                   timing["two_full_jacobian_s_median"])
                partition_info = two.identity()["partition"]
            finally:
                two.close()
                one.release()
            results[label] = {"columns": len(request.keys()), "vectors": per,
                              "timing": timing, "partition": partition_info}
            log("%s: %s; one jac %.3f s, two jac %.3f s (x%.2f); partition %s" % (
                label, {v: (d["bitwise_one_vs_two"], d["repeat_bitwise"]) for v, d in per.items()},
                timing["one_full_jacobian_s_median"], timing["two_full_jacobian_s_median"],
                timing["speedup_full_jacobian"], [p["records"] for p in partition_info]))
    finally:
        budget.close()

    benchmark = None
    if not args.skip_benchmark and "B29" in results:
        log("benchmark: two independent B29 processes, one per device")
        ctx = mp.get_context("spawn")
        barrier, queue = ctx.Barrier(2), ctx.Queue()
        procs = [ctx.Process(target=bench_worker, args=(str(capture), d, args.repeats, barrier, queue))
                 for d in (0, 1)]
        for p in procs:
            p.start()
        got = dict(queue.get(timeout=1800) for _ in procs)
        for p in procs:
            p.join(60)
        if all("error" not in v for v in got.values()):
            independent = [got[d]["full_jacobian_s_median"] for d in (0, 1)]
            split = results["B29"]["timing"]["two_full_jacobian_s_median"]
            benchmark = {"independent_full_jacobian_s_median": independent,
                         "independent_throughput_jacobians_per_s": float(sum(1.0 / t for t in independent)),
                         "split_full_jacobian_s_median": split,
                         "split_throughput_jacobians_per_s": 1.0 / split,
                         "note": "independent = two separate fits, one per GPU, running concurrently; "
                                 "split = one fit over both GPUs (latency per evaluation)"}
        else:
            benchmark = {"errors": got}
        log("benchmark: %s" % benchmark)

    report = {"schema": "liouville_cops_cuda.c07_multi_gpu_check", "schema_version": 1,
              "created_utc": provenance.utc_now(), "capture": str(capture), "run": str(run),
              "points": int(G.npt), "records": int(G.indptr[-1]), "results": results,
              "benchmark": benchmark, "all_pass": bool(all_ok),
              "manifest": [manifest.backend_manifest("cuda_fp64_v1", d) for d in (0, 1)]}
    provenance.write_json_exclusive(out / "report.json", report)
    log("ALL %s; report %s" % ("PASS" if all_ok else "FAIL", out / "report.json"))
    return 0 if all_ok else 1


if __name__ == "__main__":
    sys.exit(main())
