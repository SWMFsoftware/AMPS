#!/usr/bin/env python3
"""C08 gate: fused FP64 kernels against the generic tree kernel on the full 788k set (CUDA env).

For B29, N30 and NC31 (C00b templates) on all 1150 fit points, and every selected
C00b fixture state:

* ``bitwise``: ``CudaBackend(fused=True)`` and ``CudaBackend(fused=False)``, both in tree
  mode with the same block table, give identical pressures and Jacobian columns.
  The request types are the full Jacobian, pressures only, without temperatures,
  and temperatures only.
* ``timing``: warm medians of pressures and of the full Jacobian for both, and
  the speedup.

Writes ``$COPS_OUTPUT_ROOT/cuda_c08/fused_check_<UTC>/report.json``.
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


def log(message):
    print("[%s] %s" % (time.strftime("%H:%M:%S"), message), flush=True)


def same(a, b, request):
    return (a.pNG.tobytes() == b.pNG.tobytes() and a.pRG.tobytes() == b.pRG.tobytes() and
            all(u.tobytes() == v.tobytes() for k in request.keys()
                for u, v in zip(a.columns[k], b.columns[k])))


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--capture", required=True)
    parser.add_argument("--device", type=int, default=0)
    parser.add_argument("--configs", default="B29,N30,NC31")
    parser.add_argument("--repeats", type=int, default=5)
    parser.add_argument("--output-root", default=os.environ.get("COPS_OUTPUT_ROOT", DEFAULT_OUTPUT_ROOT))
    args = parser.parse_args()

    from c07_multi_gpu_check import build
    from cops_cuda import backend as bk, evaluator, manifest, provenance
    from cops_cuda.memory import DeviceBudget
    capture = Path(args.capture).resolve()
    fixtures = json.loads((capture / "joint" / "joint_fixtures.json").read_text())
    run, G, problems = build(capture, args.configs.split(","))
    out = provenance.make_new_directory(Path(args.output_root) / "cuda_c08" /
                                        ("fused_check_" + provenance.utc_stamp()))
    results, all_ok = {}, True
    budget = DeviceBudget(args.device)
    try:
        for label, P in problems.items():
            arrays = np.load(str(capture / "joint" / ("joint_%s.npz" % label)))
            full = bk.ColumnRequest.from_problem(P)
            requests = {"full": full, "pressures": bk.PRESSURES_ONLY,
                        "no_temperature": bk.ColumnRequest(full.basis_index, full.enhancement_index,
                                                           temperature=False, kappa=True),
                        "temperature_only": bk.ColumnRequest((), (), temperature=True, kappa=False)}
            basis = P.enhancement_basis if P.nb else None
            kw = dict(device_id=args.device, budget=budget)
            vectors = record_vectors(fixtures, label)
            states = {v: bk.ParameterState.from_problem(P, arrays["%s__%s__x" % (label, v)]) for v in vectors}
            state0 = bk.ParameterState.from_problem(P, P.x0)
            stored, timing = {}, {}
            # One backend at a time on the same device (two full datasets exceed the budget).
            for name, fused_flag in (("generic", False), ("fused", True)):
                with evaluator.CudaBackend(G, run, basis, P.enhancement_scale, fused=fused_flag, **kw) as be:
                    for v, state in states.items():
                        for rname, r in requests.items():
                            stored[(name, v, rname)] = be.evaluate(state, r)
                    be.evaluate(state0, full)
                    tp, tj = [], []
                    for _ in range(args.repeats):
                        t = time.perf_counter(); be.evaluate(state0); tp.append(time.perf_counter() - t)
                        t = time.perf_counter(); be.evaluate(state0, full); tj.append(time.perf_counter() - t)
                    timing[name] = {"pressures_s": float(np.median(tp)), "full_jacobian_s": float(np.median(tj))}
            per = {v: {rname: bool(same(stored[("generic", v, rname)], stored[("fused", v, rname)], r))
                       for rname, r in requests.items()} for v in vectors}
            all_ok = all_ok and all(all(d.values()) for d in per.values())
            timing["speedup_full_jacobian"] = timing["generic"]["full_jacobian_s"] / timing["fused"]["full_jacobian_s"]
            timing["speedup_pressures"] = timing["generic"]["pressures_s"] / timing["fused"]["pressures_s"]
            results[label] = {"columns": len(full.keys()), "bitwise": per, "timing": timing}
            log("%s: bitwise %s; jacobian %.3f -> %.3f s (x%.2f); pressures %.3f -> %.3f s" % (
                label, all(all(d.values()) for d in per.values()),
                timing["generic"]["full_jacobian_s"], timing["fused"]["full_jacobian_s"],
                timing["speedup_full_jacobian"], timing["generic"]["pressures_s"], timing["fused"]["pressures_s"]))
    finally:
        budget.close()
    report = {"schema": "liouville_cops_cuda.c08_fused_check", "schema_version": 1,
              "created_utc": provenance.utc_now(), "capture": str(capture),
              "points": int(G.npt), "records": int(G.indptr[-1]), "results": results,
              "all_pass": bool(all_ok), "manifest": manifest.backend_manifest("cuda_fp64_v1", args.device)}
    provenance.write_json_exclusive(out / "report.json", report)
    log("ALL %s; report %s" % ("PASS" if all_ok else "FAIL", out / "report.json"))
    return 0 if all_ok else 1


def record_vectors(fixtures, label):
    wanted = ("prior", "legacy_strong_clip", "archived_b29_state", "enhancement_strong_clip",
              "template_on_bound", "temperature_wide")
    return [v for v in wanted if v in fixtures["configs"][label]["vectors"]]


if __name__ == "__main__":
    sys.exit(main())
