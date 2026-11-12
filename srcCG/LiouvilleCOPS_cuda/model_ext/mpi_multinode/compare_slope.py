#!/usr/bin/env python3
# ==== version-controlled copy: header added 2026-11-11, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/mpi_multinode/compare_slope.py (on Voyager: /data001/highresolu_vt/cops_output/mpi_multinode/compare_slope.py)
# SHA256 of the run copy = SHA256 of this file without these header lines: e82668278d92d928dd5cf4a4456c6bbd81b452bcde673d341199e3e2f802750f
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: tool (test comparison, no physics). Slope-multiplier fit: the one-node local-transport test against the MPI
#   tests: the cost of every evaluation (bitwise as text), start check, final cost, pressures, restart x, timing, node
#   memory (2026-11-10).
# RUN: .venv/bin/python compare_slope.py REF_SEG_DIR CASE_SEG_DIR ...   (123646 / 123648 vs 123635: BITWISE;
#   compare_vs_123635.json in runs/shap8_012m_lategrid_dz1.5/slope_mult_tests/mpi_*_srccg_20261110T094013Z/)
# ==== end of header ====
"""Slope-multiplier fit: one-node local-transport test vs MPI tests (2026-11-10; no physics).

    .venv/bin/python compare_slope.py REF_SEG_DIR CASE_SEG_DIR [CASE_SEG_DIR ...]

Each argument is a segment directory (seg01) of slope_mult_fit.py; its provenance directory is the one next to it
(the only provenance-* of the fit directory, or the job of segments.tsv). Compared with the reference:
  * evaluations.txt (the driver's per-evaluation log: kind, seconds, data cost repr): the cost of EVERY evaluation,
    as text (bitwise, 17 significant digits) and as the largest relative difference;
  * the start-check cost line, the final cost / nfev / message (fit_manifest.json);
  * pressure_residuals.dat model_Pa (%.15e) and start_check_model.npz (surface model at the start, binary);
  * fit_restart.npz x and template_values (binary);
  * timing: evaluations.txt median seconds per fun / jac (after the first), the registry (calls, seconds, bind),
    mpi_support medians if present (ranks/comm_rank000.json), node memory peaks (ranks/mem_*.tsv or
    memory_used_avail_MB.txt), device-check records of the slope patch (backend identity).
"""

import glob
import json
import os
import re
import sys

import numpy as np


def provenance(seg):
    d = os.path.dirname(os.path.abspath(seg))
    tsv = os.path.join(d, "segments.tsv")
    if os.path.isfile(tsv):
        k = int(os.path.basename(seg)[3:])
        rows = [l.split("\t") for l in open(tsv) if not l.startswith("seg\t")]
        job = [r for r in rows if int(r[0]) == k][-1][1]
        return os.path.join(d, "provenance-%s" % job)
    p = sorted(glob.glob(os.path.join(d, "provenance-*")))
    return p[-1]


def load(seg):
    seg = os.path.abspath(seg)
    P = provenance(seg)
    out = {"seg": seg, "provenance": P}
    ev = [l.split() for l in open(os.path.join(seg, "evaluations.txt")) if not l.startswith("#")]
    out["eval_kind"] = [e[1] for e in ev]
    out["eval_seconds"] = [float(e[2]) for e in ev]
    out["eval_cost"] = [e[3] for e in ev]
    log = os.path.join(P, "fit.log")
    txt = open(log, errors="replace").read()
    m = re.search(r"start check \(variant D2\): cost at the start ([0-9.eE+-]+)", txt)
    out["start_cost"] = m.group(1) if m else None
    m = re.findall(r"^gpu_fit: [0-9.]+ s wall.*$", txt, re.M)
    out["wall_line"] = m[-1] if m else None
    m = re.findall(r"^gpu_model: bound .*$", txt, re.M)
    out["bind_line"] = m[0][:200] if m else None
    man = json.load(open(os.path.join(seg, "fit_manifest.json")))
    o = man["optimizer_results"][0]
    out["final"] = {k: o[k] for k in ("cost", "nfev", "message", "status")}
    reg = man["compute_backend"]["registry"]
    out["registry"] = {"calls": reg["calls"], "seconds": {k: round(v, 2) for k, v in reg["seconds"].items()},
                       "transport": reg.get("transport"), "activity": reg.get("activity")}
    checks = []
    for b in reg.get("backends") or []:
        ext = (b.get("model_extension") or {}).get("slope_multiplier") or {}
        for w in b.get("workers") or []:
            sm = ((w.get("model_extension") or {}).get("slope_multiplier") or {}).get("backend") or {}
            if sm:
                checks.append({k: sm.get(k) for k in ("host_check_max_abs_diff", "host_checks", "active_flips",
                                                       "colbasis_check") if k in sm})
    out["device_checks"] = checks
    mp = [float(l.split()[7]) for l in open(os.path.join(seg, "pressure_residuals.dat"))
          if not l.startswith("#") and l.strip()]
    arr = {"model_Pa_final": np.array(mp)}
    with np.load(os.path.join(seg, "start_check_model.npz"), allow_pickle=False) as q:
        for k in ("surface_pNG", "surface_pRG", "grain"):
            arr["start_" + k] = np.array(q[k])
    with np.load(os.path.join(seg, "fit_restart.npz"), allow_pickle=False) as z:
        for k in ("x", "template_values", "a", "kappa"):
            arr["restart_" + k] = np.array(z[k])
    out["_arrays"] = arr
    comm = os.path.join(P, "ranks", "comm_rank000.json")
    out["mpi_support"] = json.load(open(comm))["summary"] if os.path.isfile(comm) else None
    mem = {}
    for f in glob.glob(os.path.join(P, "ranks", "mem_*.tsv")):
        used = [int(l.split("\t")[4]) for l in open(f) if not l.startswith("utc")]
        mem[os.path.basename(f)[4:-4]] = round(max(used) * 1024 / 1e9, 1)
    f = os.path.join(P, "memory_used_avail_MB.txt")
    if os.path.isfile(f):
        used = [int(l.split()[1]) for l in open(f) if len(l.split()) >= 3]
        mem["node (memory_used_avail_MB)"] = round(max(used) * 1.048576e6 / 1e9, 1) if used else None
    out["node_peak_used_GB"] = mem
    return out


def med(v):
    return float(np.median(v)) if v else None


def summary(c):
    fun = [s for k, s in zip(c["eval_kind"], c["eval_seconds"]) if k == "fun"]
    jac = [s for k, s in zip(c["eval_kind"], c["eval_seconds"]) if k == "jac"]
    s = {"seg": c["seg"], "start_cost": c["start_cost"], "final": c["final"], "registry": c["registry"],
         "fun_median_s_after_first": med(fun[1:]), "jac_median_s_after_first": med(jac[1:]), "n_fun": len(fun),
         "n_jac": len(jac), "wall_line": c["wall_line"], "bind_line": c["bind_line"],
         "node_peak_used_GB": c["node_peak_used_GB"], "device_checks": c["device_checks"][:3],
         "n_device_check_records": len(c["device_checks"])}
    if c["mpi_support"]:
        m = c["mpi_support"]
        s["mpi_support"] = {k: {x: m[k].get(x) for x in ("calls", "median_s_after_first", "median_worker_max_s",
                                                          "median_overhead_s", "median_request_bytes",
                                                          "median_reply_bytes", "bcast_bytes_per_rank")}
                            for k in ("pressures", "jacobian") if k in m}
        s["mpi_support"]["binds"] = [{k: b.get(k) for k in ("seconds", "workers", "activity")} for b in m["binds"]]
    return s


def compare(r, c):
    out = {"evaluations": len(r["eval_cost"]), "evaluations_case": len(c["eval_cost"])}
    n = min(len(r["eval_cost"]), len(c["eval_cost"]))
    same = [r["eval_cost"][k] == c["eval_cost"][k] for k in range(n)]
    out["eval_costs_text_identical_all"] = all(same) and len(r["eval_cost"]) == len(c["eval_cost"])
    out["first_eval_cost_difference"] = next((k + 1 for k in range(n) if not same[k]), None)
    rel = [abs(float(a) - float(b)) / abs(float(a)) for a, b in zip(r["eval_cost"], c["eval_cost"])
           if a != "-" and b != "-"]
    out["eval_cost_max_rel"] = max(rel) if rel else None
    out["start_cost_identical"] = r["start_cost"] == c["start_cost"]
    out["final_cost_rel"] = abs(r["final"]["cost"] - c["final"]["cost"]) / abs(r["final"]["cost"])
    for k in sorted(set(r["_arrays"]) & set(c["_arrays"])):
        a, b = r["_arrays"][k], c["_arrays"][k]
        with np.errstate(divide="ignore", invalid="ignore"):
            d = np.abs(a - b) / np.maximum(np.abs(a), 1e-300)
        out[k] = {"bitwise": a.tobytes() == b.tobytes(), "n_differ": int(np.sum(a != b)),
                  "max_rel": float(np.nanmax(d)) if a.size else 0.0}
    return out


def main():
    ref = load(sys.argv[1])
    res = {"reference": summary(ref), "cases": []}
    for s in sys.argv[2:]:
        c = load(s)
        res["cases"].append({"summary": summary(c), "vs_reference": compare(ref, c)})
    print(json.dumps(res, indent=1, default=str))
    return 0


if __name__ == "__main__":
    sys.exit(main())
