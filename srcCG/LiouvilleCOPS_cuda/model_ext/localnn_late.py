#!/usr/bin/env python3
# ==== version-controlled copy: header added 2026-11-18, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/model_ext/localnn_late.py (on Voyager: /data001/highresolu_vt/cops_output/model_ext/localnn_late.py)
# SHA256 of the run copy = SHA256 of this file without these header lines: d6b3e8fae6bfdd6ca963657dec6011b4284b6c377bada22418759ebf7eee695d
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: late kernels of the CPU fine tuning (the 9 row masks of the frozen positive_finetune.late_only_facets, the
#   271 late facets, G_NG / G_RG of every row; optional --state pressures); 2026-11-18: meta law record (h2o / co2
#   q_min ratios) and the law check in load().
# RUN 2026-11-18: late_kernels_150m_h2oqmin0.02.npz (136507).
# ==== end of header ====
"""Option C stage 2 (lead / user 2026-11-15): the LATE-FACET KERNELS of the base for the CPU fine tuning of the localnn
chains (localnn_ft.py) and for the held fine-tuning term exp(h_e) of the localnn fits (localnn_problem.py held=).

The fine tuning of the nonneg chains (frozen nonneg_finetune.py / positive_finetune.py, model_ext_frozen/
nonneg_20261113T035423Z, used unchanged) needs from its engine (GpuModelEngine there): the facet masks seen(row mask,
structural) that positive_finetune.late_only_facets asks for, the fit state's surface pressures of every row (forward),
and G = dp/dA of the parameter facets (sensitivities). All three come from the base's kernel records here, with the
frozen contribution arithmetic (nonneg_contrib.record_contributions: FP32 record fields, FP64 arithmetic):
  1. the row masks late_only_facets requests are recorded by calling it with a recording engine (exactly its code);
  2. pass 1 over the records: per row the facets with a record with Kn != 0 | KPhi != 0 | KV != 0 (structural) or
     Kn > 0 | KPhi > 0 | KV > 0 (the '> 0' predicate), OR-ed into the facet mask of every request that selects the row;
  3. late_only_facets is called again with these answers -> the late facets, its record, the support rows;
  4. pass 2: G_NG, G_RG [row, late facet] = the sum over the row's records of the facet of the H2O contribution at
     A_e = 1 (kappa 0) = dp/dA (instantaneous; every base row); optionally the all-row instantaneous surface pressures
     of given per-facet states (--state NAME=facets.npz:KAPPA), e.g. a softplus fit state for the verification.

CLI: localnn_late.py --base BASE --d2 D2 --point-table PT --contrib-code FROZEN_CONTRIB --nonneg-code FROZEN_NONNEG
     --out late_kernels.npz [--window-min 5] [--state NAME=PATH:KAPPA ...] [--procs N]
Output: late_kernels.npz (late, late_idx, G_NG, G_RG [nrow x nlate], requests: row masks + structural flags + answers,
support rows, the late_only record, state pressures P_<NAME>_NG / _RG, meta_json) + .json report.
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
APP = Path(os.environ.get("APP", "/data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda"))
for _p in (str(HERE), str(APP)):
    if _p not in sys.path:
        sys.path.insert(0, _p)
import localnn_design as ld        # noqa: E402  (contrib_setup: the frozen arithmetic)

SCHEMA = "localnn_late_kernels_v1"
_W = {}


def say(text):
    print("localnn_late: [%s] %s" % (dt.datetime.now(dt.timezone.utc).strftime("%H:%M:%SZ"), text), flush=True)


def sha(path):
    return ld.sha256_file(path)


def build_rows(base, d2, point_table, nonneg_code, crossing_weight=10.0):
    """positive_finetune.Rows exactly as nonneg_finetune.NonnegFitContext builds it (D2's selection)."""
    if str(nonneg_code) not in sys.path:
        sys.path.insert(0, str(nonneg_code))
    import positive_finetune as pf
    import positive_problem as pp
    import grain_crossing as gc
    for mod in (pf, pp, gc):
        if Path(mod.__file__).resolve().parent != Path(nonneg_code).resolve():
            raise SystemExit("%s imported from %s, not the frozen %s" % (mod.__name__, mod.__file__, nonneg_code))
    from cops_cuda import reference
    reference.bind()
    fit_cops = reference.load("fit_cops")
    header, names_t, table_rows = fit_cops.read_model_file(str(Path(base) / "COPS_model.dat"))
    cols = {n: i for i, n in enumerate(names_t)}
    ipoint = np.array([int(r[cols["iPoint"]]) for r in table_rows])
    with np.load(Path(d2) / "fit_restart.npz", allow_pickle=False) as z:
        use_ng, use_rg, pid = np.asarray(z["use_ng"]), np.asarray(z["use_rg"]), np.asarray(z["point_ids"])
    if not np.array_equal(pid, ipoint):
        raise SystemExit("D2's rows differ from the base's")
    refNG = np.array([float(r[cols["refNG"]]) for r in table_rows])
    refRG = np.array([float(r[cols["refRG"]]) for r in table_rows])
    utc = [r[cols["UTC"]] for r in table_rows]
    ets = gc.et_by_point(Path(base))
    et = np.array([ets[int(i)] for i in ipoint])
    layout = pp.LateLayout(pp.read_point_table(str(point_table)), pp.read_run_points(Path(base)), ipoint)
    window = pp.WindowOperator(layout, np.arange(len(ipoint)))
    rows = pf.Rows(ipoint, et, refNG, refRG, use_ng, use_rg, window, crossing_weight, gc.GAP_SEARCH_S,
                   (gc.A_UNIT, gc.B_UNIT))
    return {"pf": pf, "pp": pp, "gc": gc, "rows": rows, "layout": layout, "utc": utc, "ipoint": ipoint, "et": et,
            "header": header, "table_rows": table_rows, "cols": cols}


class RecordingEngine:
    def __init__(self, nf):
        self.nf, self.requests = nf, []

    def seen(self, mask, structural=False):
        self.requests.append((np.asarray(mask, bool).copy(), bool(structural)))
        return np.zeros(self.nf, bool)


class AnswerEngine:
    """seen() from precomputed answers (exact request match), for late_only_facets and localnn_ft's CPU engine."""

    def __init__(self, requests, answers):
        self.requests, self.answers = requests, answers

    def seen(self, mask, structural=False):
        mask = np.asarray(mask, bool)
        for (rm, st), ans in zip(self.requests, self.answers):
            if st == bool(structural) and np.array_equal(rm, mask):
                return ans.copy()
        raise SystemExit("seen(): a row mask that was not precomputed (structural %s, %d rows)" % (structural,
                                                                                                  int(mask.sum())))


class Ctx:
    def __init__(self, utc, area, A_orig):
        self.utc, self.area, self.A_orig = utc, area, A_orig


def pass1(task):
    r, ip = task
    NM = _W["CS"]["NM"]
    kp, h, rec = NM.read_kernel_full(Path(_W["base"]) / "kernels" / ("point=%04d.bin" % ip))
    if kp != ip:
        raise RuntimeError("kernel file %d stores %d" % (ip, kp))
    nf = _W["nf"]
    face = rec["face"].astype(np.int64)
    s = np.zeros(nf, bool)
    s[face[(rec["Kn"] != 0) | (rec["KPhi"] != 0) | (rec["KV"] != 0)]] = True
    p = np.zeros(nf, bool)
    p[face[(rec["Kn"] > 0) | (rec["KPhi"] > 0) | (rec["KV"] > 0)]] = True
    return r, np.packbits(s), np.packbits(p)


def pass2(task):
    r, ip = task
    CS = _W["CS"]
    nc, NM = CS["nc"], CS["NM"]
    kp, h, rec = NM.read_kernel_full(Path(_W["base"]) / "kernels" / ("point=%04d.bin" % ip))
    pos = _W["late_pos"][rec["face"].astype(np.int64)]
    sel = pos >= 0
    out = {"row": r}
    nc._G["kappa"] = 0.0
    if sel.any():
        _, ((gN, gR),), _ = nc.record_contributions(rec[sel], h, (_W["ONE"],), "fp32", r)
        out["G_NG"] = np.bincount(pos[sel], weights=gN, minlength=_W["nlate"])
        out["G_RG"] = np.bincount(pos[sel], weights=gR, minlength=_W["nlate"])
    else:
        out["G_NG"] = out["G_RG"] = np.zeros(_W["nlate"])
    for name, A, kappa in _W["states"]:
        nc._G["kappa"] = kappa
        _, ((ng, rg),), _ = nc.record_contributions(rec, h, (A,), "fp32", r)
        out["P_" + name] = (float(np.sum(ng)), float(np.sum(rg)))
    nc._G["kappa"] = 0.0
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    for k in ("base", "d2", "point-table", "contrib-code", "nonneg-code", "out"):
        ap.add_argument("--" + k, required=True)
    ap.add_argument("--window-min", type=float, default=5.0)
    ap.add_argument("--state", action="append", default=[], help="NAME=facets.npz:KAPPA (key A, or A_final)")
    ap.add_argument("--procs", type=int, default=int(os.environ.get("SLURM_CPUS_PER_TASK", "8")))
    args = ap.parse_args(argv)
    t0 = time.time()
    out = Path(args.out)
    if out.exists():
        raise SystemExit("%s exists" % out)
    with np.load(Path(args.d2) / "fit_restart.npz", allow_pickle=False) as z:
        T = tuple(float(v) for v in np.asarray(z["T"]))
    CS = ld.contrib_setup(args.base, args.contrib_code, T)
    el = CS["el"]
    nf = el["nfaces"]
    R = build_rows(args.base, args.d2, args.point_table, args.nonneg_code)
    pf, rows = R["pf"], R["rows"]
    own, late_used = rows.objective_points()
    ctx = Ctx(R["utc"], el["area"], el["AH2O"])
    rec_eng = RecordingEngine(nf)
    try:
        pf.late_only_facets(rec_eng, rows, ctx, args.window_min, own, late_used)
    except SystemExit:
        pass
    reqs = rec_eng.requests
    say("late_only_facets requests %d row masks (%s)" % (len(reqs), [(int(m.sum()), s) for m, s in reqs]))
    _W.update({"CS": CS, "base": args.base, "nf": nf})
    tasks = [(i, int(R["ipoint"][i])) for i in range(len(R["ipoint"]))]
    answers = [np.zeros(nf, bool) for _ in reqs]
    with mp.get_context("fork").Pool(args.procs) as pool:
        for r, ps, pp_ in pool.imap_unordered(pass1, tasks, chunksize=1):
            s = p = None
            for k, (rm, st) in enumerate(reqs):
                if rm[r]:
                    if st:
                        s = np.unpackbits(ps, count=nf).astype(bool) if s is None else s
                        answers[k] |= s
                    else:
                        p = np.unpackbits(pp_, count=nf).astype(bool) if p is None else p
                        answers[k] |= p
    say("pass 1 done (%.0f s)" % (time.time() - t0))
    late, late_rec, support = pf.late_only_facets(AnswerEngine(reqs, answers), rows, ctx, args.window_min, own,
                                                  late_used)
    late_idx = np.flatnonzero(late)
    say("late facets %d (%s)" % (len(late_idx), json.dumps({k: late_rec[k] for k in ("t_start_utc", "t_end_utc",
                                                                                     "points_inside", "facets")})))
    late_pos = np.full(nf, -1, dtype=np.int64)
    late_pos[late_idx] = np.arange(len(late_idx))
    states = []
    for spec in args.state:
        name, rest = spec.split("=", 1)
        path, kappa = rest.rsplit(":", 1)
        with np.load(path, allow_pickle=False) as z:
            A = np.asarray(z["A_final"] if "A_final" in z.files else z["A"], np.float64)
        states.append((name, A, float(kappa)))
    _W.update({"late_pos": late_pos, "nlate": len(late_idx), "ONE": np.ones(nf), "states": states})
    nrow = len(tasks)
    G = {g: np.zeros((nrow, len(late_idx))) for g in ("NG", "RG")}
    P = {name: np.zeros((nrow, 2)) for name, _, _ in states}
    with mp.get_context("fork").Pool(args.procs) as pool:
        for res in pool.imap_unordered(pass2, tasks, chunksize=1):
            r = res["row"]
            G["NG"][r], G["RG"][r] = res["G_NG"], res["G_RG"]
            for name in P:
                P[name][r] = res["P_" + name]
    say("pass 2 done (%.0f s)" % (time.time() - t0))
    meta = {"schema": SCHEMA, "created_utc": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
            "host": socket.gethostname(), "slurm_job_id": os.environ.get("SLURM_JOB_ID"), "command": sys.argv,
            "base": str(Path(args.base).resolve()), "d2": str(Path(args.d2).resolve()),
            "point_table": {"path": str(Path(args.point_table).resolve()), "sha256": sha(args.point_table)},
            "contrib_code": str(Path(args.contrib_code).resolve()), "nonneg_code": str(Path(args.nonneg_code).resolve()),
            "law": CS["law"],
            "nonneg_code_sha256": {m: sha(Path(args.nonneg_code) / m) for m in ("positive_finetune.py",
                                                                              "positive_problem.py", "grain_crossing.py",
                                                                              "nonneg_finetune.py")},
            "window_min": args.window_min, "T": T, "late_only": late_rec, "late_facets": int(len(late_idx)),
            "requests": [{"rows": int(m.sum()), "structural": s, "facets": int(a.sum())} for (m, s), a in
                         zip(reqs, answers)],
            "states": [{"name": n, "spec": s} for (n, _, _), s in zip(states, args.state)],
            "definition": "G_* = d p / d A_e (H2O, instantaneous, every base row) of the late facets; FP32 records, FP64",
            "code_sha256": {"localnn_late.py": sha(__file__), "localnn_design.py": sha(HERE / "localnn_design.py")},
            "wall_seconds": time.time() - t0}
    arrays = {"late": late, "late_idx": late_idx, "G_NG": G["NG"], "G_RG": G["RG"],
              "req_rows": np.stack([m for m, _ in reqs]), "req_structural": np.array([s for _, s in reqs]),
              "support_before": support["before_rows"], "support_inside": support["inside_rows"],
              "support_t": np.array([support["t_start"], support["t_end"]]),
              "meta_json": np.array(json.dumps(meta, sort_keys=True, default=float))}
    for k, a in enumerate(answers):
        arrays["req_answer_%d" % k] = np.flatnonzero(a)
    for name in P:
        arrays["P_%s_NG" % name], arrays["P_%s_RG" % name] = P[name][:, 0], P[name][:, 1]
    tmp = out.parent / (out.name + ".incomplete")
    with open(tmp, "xb") as f:
        np.savez(f, **arrays)
    tmp.rename(out)
    with open(str(out) + ".json", "x") as f:
        json.dump(meta, f, indent=1, sort_keys=True, default=float)
    say("done: %s (%.0f s)" % (out, time.time() - t0))
    return 0


def load(path, nf):
    """(late mask, late_idx, G_NG, G_RG, AnswerEngine, support dict, meta, state pressures)."""
    with np.load(path, allow_pickle=False) as z:
        meta = json.loads(str(z["meta_json"]))
    # 2026-11-18: the production law of the late kernels vs the fitter's (localnn_problem.law_check)
    if "h2o_qmin_ratio" in (meta.get("law") or {}):
        ratio, src = float(meta["law"]["h2o_qmin_ratio"]), "late-kernel meta law.h2o_qmin_ratio"
    else:
        ratio, src = lp.contrib_code_h2o_qmin_ratio(meta["contrib_code"]), "its contrib code %s" % meta["contrib_code"]
    import localnn_problem as lp
    law_rec = lp.law_check("late kernels", Path(path).resolve(), ratio, src)
    with np.load(path, allow_pickle=False) as z:
        rows = np.asarray(z["req_rows"])
        st = np.asarray(z["req_structural"])
        answers = []
        for k in range(len(st)):
            a = np.zeros(nf, bool)
            a[np.asarray(z["req_answer_%d" % k])] = True
            answers.append(a)
        out = {"late": np.asarray(z["late"]), "late_idx": np.asarray(z["late_idx"]), "G_NG": np.asarray(z["G_NG"]),
               "G_RG": np.asarray(z["G_RG"]), "engine": AnswerEngine([(rows[k], bool(st[k])) for k in range(len(st))],
                                                                     answers),
               "support": {"before_rows": np.asarray(z["support_before"]), "inside_rows": np.asarray(z["support_inside"]),
                           "t_start": float(z["support_t"][0]), "t_end": float(z["support_t"][1])},
               "meta": meta, "states": {}, "law_check": law_rec}
        for k in z.files:
            if k.startswith("P_") and k.endswith("_NG"):
                name = k[2:-3]
                out["states"][name] = (np.asarray(z[k]), np.asarray(z["P_%s_RG" % name]))
    return out


if __name__ == "__main__":
    sys.exit(main())
