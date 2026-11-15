#!/usr/bin/env python3
# ==== version-controlled copy: header added 2026-11-14, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/model_ext/localnn_direct.py (on Voyager: /data001/highresolu_vt/cops_output/model_ext/localnn_direct.py)
# SHA256 of the run copy = SHA256 of this file without these header lines: 2f1e179ec875f9c7d5b5d19ccf5967173518e51903cacc2e1ad07247c1183793
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: check (option C): a fit segment's final per-facet activity re-evaluated at every base row by direct facet sums
#   (frozen nonneg_contrib arithmetic), compared with the segment's cost and tables (design route).
# RUN: jobs/localnn_fit.sbatch after a complete fit (frozen model_ext_frozen/localnn_20261114T081912Z).
# ==== end of header ====
"""Option C (local non-negative basis): independent re-evaluation of a fit segment's FINAL per-facet activity at every
base row by direct facet sums with the frozen contribution arithmetic (model_ext_frozen/nonneg_contrib_20261112T210012Z
record_contributions: FP32 records, FP64 arithmetic; it reproduces the GPU tables of D2 and of the softplus c9 state to
their printed digits, localnn_fit.py verify), then W10x5 + grain + the fit's cost (localnn_problem.Data).

Compared with the segment's products (design route): the data cost and crossing cost (fit_result.txt / fit_manifest.json),
pNG / pRG of every row of COPS_model.dat (printed digits) and the surface / grain columns of late_window_model.dat.
Written to <fit dir>/direct_check_<segNN>/ (direct_check.json, COPS_model_direct.dat = the same table with the direct
values); exit 1 if the design route and the direct sums differ by more than 1e-11 (relative, any row) or the costs by
more than 1e-10.

CLI: localnn_direct.py --seg <fit>/segNN --base BASE --d2 D2 --point-table PT --design DESIGN --contrib-code FROZEN
     [--procs N]
"""

import argparse
import datetime as dt
import json
import multiprocessing as mp
import os
import socket
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))
import localnn_design as ld        # noqa: E402
import localnn_problem as lp       # noqa: E402

_S = {}


def work(task):
    r, ip = task
    return r, ld.direct_row(_S["CS"], _S["base"], r, ip, _S["A"], _S["kappa"])


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    for k in ("seg", "base", "d2", "point-table", "design", "contrib-code"):
        ap.add_argument("--" + k, required=True)
    ap.add_argument("--procs", type=int, default=int(os.environ.get("SLURM_CPUS_PER_TASK", "8")))
    args = ap.parse_args(argv)
    t0 = time.time()
    seg = Path(args.seg).resolve()
    out = seg.parent / ("direct_check_%s" % seg.name)
    if out.exists():
        raise SystemExit("%s exists" % out)
    with np.load(seg / "nonneg_facets.npz", allow_pickle=False) as z:
        A = np.asarray(z["A"], np.float64)
    with np.load(seg / "localnn_state.npz", allow_pickle=False) as z:
        kappa = float(z["kappa"])
        gA, gB = float(z["grain_A_uPa"]), float(z["grain_B_uPa_per_100s"])
        meta = json.loads(str(z["meta_json"]))
    design = lp.load_design(args.design)
    D = lp.Data(args.base, args.d2, args.point_table, design)
    CS = ld.contrib_setup(args.base, args.contrib_code, lp.USER_T0)
    _S.update({"CS": CS, "base": args.base, "A": A, "kappa": kappa})
    sN, sR = np.zeros(D.nrow), np.zeros(D.nrow)
    tasks = [(r, int(D.ipoint[r])) for r in range(D.nrow)]
    with mp.get_context("fork").Pool(args.procs) as pool:
        for r, (n_, r_) in pool.imap_unordered(work, tasks, chunksize=1):
            sN[r], sR[r] = n_, r_
    rows = D.model_rows(sN, sR, gA, gB)
    cst = D.costs(rows["pNG"], rows["pRG"])
    cm = lp.read_table_cols(seg / "COPS_model.dat", ("pNG", "pRG"))
    lw = lp.read_table_cols(seg / "late_window_model.dat", ("surface_pNG_Pa", "surface_pRG_Pa", "grain_NG_Pa"))
    res = {"created_utc": lp.utc_now(), "host": socket.gethostname(), "slurm_job_id": os.environ.get("SLURM_JOB_ID"),
           "command": sys.argv, "segment": str(seg),
           "data_cost_direct": cst["data_cost"], "data_cost_segment": meta["costs"]["data_cost"],
           "data_cost_rel_diff": abs(cst["data_cost"] / meta["costs"]["data_cost"] - 1.0),
           "crossing_cost_direct": cst["crossing_cost"], "crossing_cost_segment": meta["costs"]["crossing_cost"],
           "COPS_model_pNG": lp.compare(rows["pNG"], {k: v["pNG"] for k, v in cm.items()}, D.ipoint, "%.6e"),
           "COPS_model_pRG": lp.compare(rows["pRG"], {k: v["pRG"] for k, v in cm.items()}, D.ipoint, "%.6e"),
           "late_window_surface_pNG": lp.compare(rows["surface_NG"], {k: v["surface_pNG_Pa"] for k, v in lw.items()},
                                                 D.ipoint, "%.9e"),
           "late_window_surface_pRG": lp.compare(rows["surface_RG"], {k: v["surface_pRG_Pa"] for k, v in lw.items()},
                                                 D.ipoint, "%.9e")}
    cN, cR = lp.all_rows_surface(design, np.asarray(np.load(seg / "localnn_state.npz")["c"]),
                                 np.asarray(np.load(seg / "localnn_state.npz")["m"]), kappa, meta["templates"] == "on")
    res["design_vs_direct_instantaneous_max_rel"] = [float(np.max(np.abs(cN / sN - 1))), float(np.max(np.abs(cR / sR - 1)))]
    res["passed"] = bool(res["data_cost_rel_diff"] < 1e-10 and max(res["design_vs_direct_instantaneous_max_rel"]) < 1e-11
                         and res["COPS_model_pNG"]["max_rel"] <= 5e-7 * (1 + 1e-6)
                         and res["COPS_model_pRG"]["max_rel"] <= 5e-7 * (1 + 1e-6))
    res["wall_seconds"] = time.time() - t0
    stage = out.parent / ("%s.incomplete_%s" % (out.name, dt.datetime.now(dt.timezone.utc).strftime("%Y%m%dT%H%M%SZ")))
    stage.mkdir()
    with open(stage / "direct_check.json", "x") as f:
        json.dump(res, f, indent=1, sort_keys=True, default=float)
    lines = (seg / "COPS_model.dat").read_text().splitlines()
    k_hdr = next(i for i, line in enumerate(lines) if line.startswith("#iPoint"))
    names = lines[k_hdr][1:].split()
    c = {n: i for i, n in enumerate(names)}
    pos = {int(ip): k for k, ip in enumerate(D.ipoint)}
    outl = ["# LocalnnDirect (model_ext/localnn_direct.py, %s): pNG/pRG of every row from DIRECT facet sums of the "
            "segment's per-facet activity (frozen nonneg_contrib.record_contributions), W10x5 and grain as the fit"
            % res["created_utc"]] + lines[:k_hdr + 1]
    for line in lines[k_hdr + 1:]:
        if line.startswith("#") or not line.strip():
            outl.append(line)
            continue
        v = line.split()
        k = pos[int(v[c["iPoint"]])]
        v[c["pNG"]], v[c["pRG"]] = "%.6e" % rows["pNG"][k], "%.6e" % rows["pRG"][k]
        outl.append(" ".join(v))
    with open(stage / "COPS_model_direct.dat", "x") as f:
        f.write("\n".join(outl) + "\n")
    stage.rename(out)
    lp.say("direct check %s: %s" % ("PASSED" if res["passed"] else "FAILED", json.dumps(
        {k: res[k] for k in ("data_cost_direct", "data_cost_segment", "data_cost_rel_diff",
                             "design_vs_direct_instantaneous_max_rel")})))
    return 0 if res["passed"] else 1


if __name__ == "__main__":
    sys.exit(main())
