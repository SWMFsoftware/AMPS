#!/usr/bin/env python3
# ==== version-controlled copy: header added 2026-11-11, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/mpi_multinode/compare_tests.py (on Voyager: /data001/highresolu_vt/cops_output/mpi_multinode/compare_tests.py)
# SHA256 of the run copy = SHA256 of this file without these header lines: c47648c2620bd883bbb56db98262128ce8b75f012431a4427b11576438f9836e
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: tool (test comparison, no physics). Compares short test fits of jobs/mpi_fit.sbatch (local vs MPI 8 vs 16
#   ranks): start cost, cost history with the SHA-256 of x, residuals and Jacobian per evaluation, predictions bitwise
#   or max rel, timing, memory, bytes per message (2026-11-10).
# RUN: .venv/bin/python compare_tests.py [--json OUT] [--md OUT] REF_DIR DIR ...   (clip D2 tests 123629-123638:
#   runs/shap8_012m_lategrid_dz1.5/clip_refit_tests/mpi_20261110/compare_*.{md,json})
# ==== end of header ====
"""Compare short test fits run through jobs/mpi_fit.sbatch (local vs MPI 8 vs MPI 16 ranks; 2026-11-10; no physics).

    .venv/bin/python compare_tests.py [--json OUT.json] [--md OUT.md] REF_DIR DIR [DIR ...]

Each DIR is a fit directory of the launcher (seg<NN>/, provenance-<job>/, segments.tsv). For every DIR, segment 1 (or
--seg K) and its provenance (the job recorded in segments.tsv) are read:
  * the start cost (clip_refit / positive start check line), the optimizer's cost history and the SHA-256 of x, of the
    residual vector and of the Jacobian at every optimizer evaluation (provenance ranks/optimizer_history.json written
    by model_ext/mpi_support.py), the final cost / nfev / message (fit_manifest.json);
  * the predictions pred_ng, pred_rg and x of fit_restart.npz, compared with REF_DIR: BITWISE (bytes) or the largest
    relative difference;
  * timing: model_ext/mpi_support.py medians per pressure / Jacobian evaluation (after the first), the slowest worker,
    the overhead (communication + host work), the registry totals (bind, pressures, jacobian), the gpu_fit wall line;
  * memory: per node the peak of MemTotal - MemAvailable (ranks/mem_<host>.tsv), per rank the peak RSS;
  * communication per evaluation: request / reply bytes, Bcast bytes per rank (host activity).
"""

import argparse
import glob
import json
import os
import re
import sys
from pathlib import Path

import numpy as np


def segment_job(d, seg):
    rows = [l.rstrip("\n").split("\t") for l in open(os.path.join(d, "segments.tsv")) if not l.startswith("seg\t")]
    rows = [r for r in rows if int(r[0]) == seg]
    if not rows:
        raise SystemExit("%s: no segment %d in segments.tsv" % (d, seg))
    r = rows[-1]
    return {"seg": seg, "job": r[1], "start_utc": r[2], "end_utc": r[3], "nodes": int(r[4]), "node_list": r[5],
            "ranks": int(r[6]), "transport": r[7], "activity": r[8], "status": r[11], "exit": r[12]}


def grep(path, pattern):
    out = []
    try:
        for line in open(path, errors="replace"):
            if re.search(pattern, line):
                out.append(line.rstrip("\n"))
    except OSError:
        pass
    return out


def load_case(d, seg):
    d = os.path.abspath(d)
    info = segment_job(d, seg)
    P = os.path.join(d, "provenance-%s" % info["job"])
    S = os.path.join(d, "seg%02d" % seg)
    log = os.path.join(P, "fit.log")
    info["dir"], info["provenance"] = d, P
    start = grep(log, r"start check \(variant|start check: cost|positive_fit: start|cost at the start")
    info["start_line"] = start[0][:300] if start else None
    m = re.search(r"cost at the start ([0-9.eE+-]+)", start[0]) if start else None
    info["start_cost"] = float(m.group(1)) if m else None
    wall = grep(log, r"^gpu_fit: [0-9.]+ s wall")
    info["gpu_fit_wall_line"] = wall[-1] if wall else None
    hist = json.load(open(os.path.join(P, "ranks", "optimizer_history.json")))
    calls = []
    for c in hist:
        ev = c["evaluations"]
        calls.append({"fun_costs": [e["cost"] for e in ev if e["kind"] == "fun"],
                      "fun_cost_repr": [e["cost_repr"] for e in ev if e["kind"] == "fun"],
                      "f_sha": [e["f_sha256"] for e in ev if e["kind"] == "fun"],
                      "x_sha": [e["x_sha256"] for e in ev],
                      "J_sha": [e.get("J_sha256") for e in ev if e["kind"] == "jac"],
                      "fun_s": [e["seconds"] for e in ev if e["kind"] == "fun"],
                      "jac_s": [e["seconds"] for e in ev if e["kind"] == "jac"],
                      "nfev": c.get("nfev"), "cost": c.get("cost")})
    info["optimizer_calls"] = calls
    man = json.load(open(os.path.join(S, "fit_manifest.json")))
    o = man["optimizer_results"][0]
    info["final"] = {"cost": o["cost"], "data_cost": o.get("data_cost"), "nfev": o["nfev"], "message": o["message"],
                     "status": o["status"]}
    reg = man["compute_backend"]["registry"]
    info["registry"] = {"calls": reg["calls"], "seconds": reg["seconds"], "transport": reg.get("transport"),
                        "activity": reg.get("activity"),
                        "events": [{k: e.get(k) for k in ("event", "rows", "seconds", "records")} for e in reg["events"]]}
    be = reg["backends"][0] if reg.get("backends") else {}
    info["partition"] = be.get("partition")
    info["workers_ready"] = [{k: w.get(k) for k in ("host", "device_id", "rows", "records", "dataset_load_s", "activity")}
                             for w in (be.get("workers_ready") or [])]
    acts = []
    for w in be.get("workers") or []:
        a = (w.get("activity") or {})
        if a.get("checks"):
            acts.append(a["checks"])
    info["activity_checks"] = acts
    z = np.load(os.path.join(S, "fit_restart.npz"), allow_pickle=False)
    info["_arrays"] = {k: np.array(z[k]) for k in ("x", "activity_linear")}
    # full-precision predictions: pressure_residuals.dat model_Pa (%.15e, 16 digits; the restart's pred_ng/pred_rg are
    # read back from COPS_model.dat with 7 digits) and start_check_model.npz (binary, the surface model at the start)
    mp = []
    for line in open(os.path.join(S, "pressure_residuals.dat")):
        if not line.startswith("#") and line.strip():
            mp.append(float(line.split()[7]))
    info["_arrays"]["model_Pa_final"] = np.array(mp)
    sc = os.path.join(S, "start_check_model.npz")
    if os.path.isfile(sc):
        with np.load(sc, allow_pickle=False) as q:
            for k in ("surface_pNG", "surface_pRG", "grain"):
                if k in q.files:
                    info["_arrays"]["start_" + k] = np.array(q[k])
    comm = os.path.join(P, "ranks", "comm_rank000.json")
    info["comm"] = json.load(open(comm))["summary"] if os.path.isfile(comm) else None
    mem = {}
    for f in glob.glob(os.path.join(P, "ranks", "mem_*.tsv")):
        host = os.path.basename(f)[4:-4]
        used = [int(l.split("\t")[4]) for l in open(f) if not l.startswith("utc")]
        mem[host] = {"peak_used_GB": max(used) * 1024 / 1e9 if used else None,
                     "start_used_GB": used[0] * 1024 / 1e9 if used else None, "samples": len(used)}
    info["node_memory"] = mem
    rss = {}
    for f in sorted(glob.glob(os.path.join(P, "ranks", "rank*.end.json"))):
        e = json.load(open(f))
        rss.setdefault(e["host"].split(".")[0], []).append(
            (e["rank"], e["maxrss_self_kB"] * 1024 / 1e9, e["maxrss_children_kB"] * 1024 / 1e9, e["exit"]))
    info["rank_peak_rss_GB"] = rss
    ranks = []
    for f in sorted(glob.glob(os.path.join(P, "ranks", "rank[0-9][0-9][0-9].json"))):
        r = json.load(open(f))
        g = r.get("gpu") or {}
        ranks.append({"rank": r["rank"], "host": r["host"], "local_rank": r["local_rank"],
                      "device": g.get("device"), "uuid": g.get("uuid"), "pci": g.get("pci_bus_id"),
                      "gpus": r.get("gpus")})
    info["rank_records"] = ranks
    return info


def compare(ref, case):
    out = {}
    for k in sorted(set(ref["_arrays"]) & set(case["_arrays"])):
        a, b = ref["_arrays"][k], case["_arrays"][k]
        same = a.tobytes() == b.tobytes()
        with np.errstate(divide="ignore", invalid="ignore"):
            rel = np.abs(a - b) / np.maximum(np.abs(a), 1e-300)
        out[k] = {"bitwise": bool(same), "max_rel": float(np.nanmax(rel)) if a.size else 0.0,
                  "n_differ": int(np.sum(a != b))}
    rc, cc = ref["optimizer_calls"], case["optimizer_calls"]
    out["optimizer_calls"] = []
    for i, (r, c) in enumerate(zip(rc, cc)):
        n = min(len(r["f_sha"]), len(c["f_sha"]))
        first_diff = next((k for k in range(n) if r["f_sha"][k] != c["f_sha"][k]), None)
        rel = [abs(a - b) / abs(a) for a, b in zip(r["fun_costs"], c["fun_costs"])]
        out["optimizer_calls"].append({
            "evaluations_ref": len(r["f_sha"]), "evaluations": len(c["f_sha"]),
            "residuals_bitwise_all": first_diff is None and len(r["f_sha"]) == len(c["f_sha"]),
            "first_residual_difference": first_diff,
            "x_sequence_bitwise": r["x_sha"] == c["x_sha"],
            "jacobians_bitwise_all": r["J_sha"] == c["J_sha"],
            "cost_max_rel_diff": max(rel) if rel else None})
    out["final_cost_rel_diff"] = abs(ref["final"]["cost"] - case["final"]["cost"]) / abs(ref["final"]["cost"])
    if ref["start_cost"] and case["start_cost"]:
        out["start_cost_rel_diff"] = abs(ref["start_cost"] - case["start_cost"]) / abs(ref["start_cost"])
    return out


def summary_row(c):
    cm = c["comm"] or {}
    p, j = cm.get("pressures") or {}, cm.get("jacobian") or {}
    bind = next((e["seconds"] for e in c["registry"]["events"] if e["event"] == "bind"), None)
    return {"case": os.path.basename(c["dir"]), "job": c["job"], "nodes": c["nodes"], "ranks": c["ranks"],
            "transport": c["transport"], "activity": c["activity"], "start_cost": c["start_cost"],
            "final_cost": c["final"]["cost"], "nfev": c["final"]["nfev"],
            "pressures_median_s": p.get("median_s_after_first"), "pressures_worker_max_s": p.get("median_worker_max_s"),
            "pressures_overhead_s": p.get("median_overhead_s"), "pressures_calls": p.get("calls"),
            "jacobian_median_s": j.get("median_s_after_first"), "jacobian_worker_max_s": j.get("median_worker_max_s"),
            "jacobian_overhead_s": j.get("median_overhead_s"), "jacobian_calls": j.get("calls"),
            "request_B": p.get("median_request_bytes"), "reply_B_pressures": p.get("median_reply_bytes"),
            "reply_B_jacobian": j.get("median_reply_bytes"), "bcast_B_per_rank": p.get("bcast_bytes_per_rank"),
            "bind_s": bind, "registry_seconds": c["registry"]["seconds"], "registry_calls": c["registry"]["calls"],
            "node_peak_used_GB": {h: round(v["peak_used_GB"], 1) for h, v in c["node_memory"].items()},
            "rank0_peak_rss_GB": next((round(r[1], 1) for v in c["rank_peak_rss_GB"].values() for r in v if r[0] == 0), None),
            "worker_peak_rss_GB_max": max((r[1] for v in c["rank_peak_rss_GB"].values() for r in v if r[0] != 0),
                                          default=None),
            "gpu_fit_wall": c["gpu_fit_wall_line"]}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("ref")
    ap.add_argument("cases", nargs="+")
    ap.add_argument("--seg", type=int, default=1)
    ap.add_argument("--json")
    ap.add_argument("--md")
    a = ap.parse_args()
    ref = load_case(a.ref, a.seg)
    cases = [load_case(d, a.seg) for d in a.cases]
    result = {"reference": summary_row(ref), "cases": []}
    lines = ["| case | nodes | ranks | transport | activity | start cost | final cost | nfev | p median s | p slowest worker s | "
             "p overhead s | J median s | J slowest worker s | J overhead s | bind s | node peak used GB | vs ref: residuals "
             "bitwise | predictions bitwise (start surface + final model_Pa) | predictions max rel |", "|" + "---|" * 19]

    def row(s, cmp):
        f = lambda v, fmt="%.4f": "-" if v is None else fmt % v
        keys = [k for k in ("model_Pa_final", "start_surface_pNG", "start_surface_pRG") if cmp is not None and k in cmp]
        pb = "-" if cmp is None else ("yes" if all(cmp[k]["bitwise"] for k in keys) else "no")
        pm = "-" if cmp is None else "%.2e" % max(cmp[k]["max_rel"] for k in keys)
        rb = "-" if cmp is None else ("yes" if all(c["residuals_bitwise_all"] for c in cmp["optimizer_calls"]) else
                                      "no (first diff %s)" % [c["first_residual_difference"] for c in cmp["optimizer_calls"]])
        return "| %s | %d | %d | %s | %s | %s | %s | %d | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s |" % (
            s["case"], s["nodes"], s["ranks"], s["transport"], s["activity"], f(s["start_cost"], "%.12e"),
            f(s["final_cost"], "%.12e"), s["nfev"], f(s["pressures_median_s"]), f(s["pressures_worker_max_s"]),
            f(s["pressures_overhead_s"]), f(s["jacobian_median_s"]), f(s["jacobian_worker_max_s"]),
            f(s["jacobian_overhead_s"]), f(s["bind_s"], "%.1f"), s["node_peak_used_GB"], rb, pb, pm)
    lines.append(row(result["reference"], None))
    for c in cases:
        cmp = compare(ref, c)
        s = summary_row(c)
        result["cases"].append({"summary": s, "vs_reference": cmp, "rank_records": c["rank_records"], "partition": c["partition"],
                                "activity_checks": c["activity_checks"], "rank_peak_rss_GB": c["rank_peak_rss_GB"],
                                "cost_history": c["optimizer_calls"][0]["fun_cost_repr"] if c["optimizer_calls"] else None})
        lines.append(row(s, cmp))
    result["reference_cost_history"] = ref["optimizer_calls"][0]["fun_cost_repr"] if ref["optimizer_calls"] else None
    result["reference_rank_records"] = ref["rank_records"]
    text = "\n".join(lines)
    print(text)
    for c in result["cases"]:
        print("\n%s vs %s:" % (c["summary"]["case"], result["reference"]["case"]))
        print(json.dumps(c["vs_reference"], indent=1))
    if a.json:
        Path(a.json).write_text(json.dumps(result, indent=1, default=str))
    if a.md:
        Path(a.md).write_text(text + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
