#!/usr/bin/env python3
# ==== version-controlled copy: header added 2026-11-14, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/model_ext/nonneg_ft_summary.py (on Voyager: /data001/highresolu_vt/cops_output/model_ext/nonneg_ft_summary.py)
# SHA256 of the run copy = SHA256 of this file without these header lines: 5b22e4e04f898a2207e50d52b31af50d77558649f76baf166e5558420cd1e8b9
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: tool. Summary of one late-facet fine tuning (user 2026-11-15 delta regularization tests): costs, late rms
#   ln NG/RG, P_delta, rms / max |h + delta|, facets at the cap, seen area holding 50 / 90 % of the seen production,
#   the tuned facets' share of the late signal (FP32 records, verified against late_window_model.dat).
# USED BY: jobs/nonneg_ft_summary.sbatch (test 134499; snapshot model_ext_frozen/nonneg_summary_20261116T031115Z).
# ==== end of header ====
"""Summary of ONE finished late-facet fine tuning of the non-negative models (nonneg_finetune.py; for the delta
regularization tests of the user's 2026-11-15 decision): from the fine tuning's own products only.

    nonneg_ft_summary.py --ft <ft seg dir> --base <base> --geometry <validated nonneg_maps geometry> --out <summary.txt>
        [--reference-ft <unregularized ft seg dir>] [--procs N]
Writes <out> and <out>.json (never overwrites).

NUMBERS
  costs, late-window rms ln NG / RG (costs_final of finetune_state.npz; the fit's selection; late window = from
    10:33:10 UTC), P_delta, rms (area-weighted) and max |h + delta| over the tuned facets, facets at the cap;
  the seen area holding 50 % / 90 % of the seen production at 10:33:10 (F_H2O(A_final) area over the facets seen by the
    fit, nonneg_facets.npz seen; decreasing F_H2O, ties lower facet index first; as nonneg_alt_summary.py);
  the share of the late signal from the tuned facets: per late observation (used rows from 10:33:10) the surface
    pressure of the tuned facets / that of all facets (fit_cops.Model term by term, the fit's FP32 records,
    nonneg_contrib.record_contributions, H2O + CO2 per facet, W10x5 of the 1-s grid rows), mean / min / ratio of sums;
    verified against the fine tuning's late_window_model.dat surface columns;
  with --reference-ft: the same numbers of the reference (e.g. the chain's unregularized c5 fine tuning) side by side.
"""

import argparse
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

import nonneg_maps as NM            # noqa: E402
import nonneg_contrib as NC         # noqa: E402
import nonneg_alt_summary as NS     # noqa: E402  (row_work, area_holding)
import positive_problem as pp       # noqa: E402

SCHEMA = "nonneg_ft_summary_v1"
LATE_UTC = "2016-09-30T10:33:10"


def one(ft, base, geometry, procs, el, fc, header, cols, rows, helio, sza):
    ft = Path(ft).resolve()
    rec = json.load(open(ft / "fine_tuning.json"))
    fit = Path(rec["fit"]["directory"]).resolve()
    with np.load(ft / "finetune_state.npz", allow_pickle=False) as z:
        meta = json.loads(str(z["meta_json"]))
        A = np.asarray(z["A_final"], float)
        tuned = np.asarray(z["seen"], bool)
        held_next = np.asarray(z["held_next"], float)
        held_prev = np.asarray(z["held_prev"], float)
        gA, gB = float(z["grain_A_uPa"]), float(z["grain_B_uPa_per_100s"])
    with np.load(fit / "nonneg_facets.npz", allow_pickle=False) as z:
        seen = np.asarray(z["seen"], bool)
    with np.load(fit / "fit_restart.npz", allow_pickle=False) as z:
        T = tuple(float(v) for v in np.asarray(z["T"]))
        kappa = float(np.asarray(z["kappa"]).reshape(-1)[0])
        use = {"NG": (np.asarray(z["use_ng"]) > 0) & (np.asarray(z["ref_ng"]) > 0),
               "RG": (np.asarray(z["use_rg"]) > 0) & (np.asarray(z["ref_rg"]) > 0)}
        utc_r = [str(u) for u in z["utc"]]
    area = el["area"]
    nf = el["nfaces"]
    cf = meta["costs_final"]
    ci = meta["costs_initial"]
    t = held_next[tuned]
    w = area[tuned]
    out = {"ft": str(ft), "fit": str(fit), "model": (rec["fit"].get("model") or {}).get("model"),
           "settings": {k: rec["settings"].get(k) for k in ("delta_penalty_f", "delta_cref", "delta_cap", "max_nfev",
                                                            "window_min", "facets")},
           "stop": meta.get("stop"), "nfev_total": meta.get("nfev_total"), "tuned_facets": int(tuned.sum()),
           "costs_initial": {k: ci.get(k) for k in ("cost_total", "cost_data", "cost_crossing", "cost_delta_penalty",
                                                    "late_window_rms_ln_NG", "late_window_rms_ln_RG")},
           "costs_final": {k: cf.get(k) for k in ("cost_total", "cost_data", "cost_crossing", "cost_delta_penalty",
                                                  "rms_ln_NG", "rms_ln_RG", "late_window_rms_ln_NG",
                                                  "late_window_rms_ln_RG", "delta_total_rms_area_weighted",
                                                  "delta_total_max_abs", "facets_at_cap", "facets_at_cap_upper",
                                                  "facets_at_cap_lower")},
           "total_log_factor": {"rms_area_weighted": float(np.sqrt(np.sum(w * t * t) / np.sum(w))),
                                "max_abs": float(np.max(np.abs(t))),
                                "held_prev_max_abs": float(np.max(np.abs(held_prev[tuned]))),
                                "delta_range": [float((held_next - held_prev)[tuned].min()),
                                                float((held_next - held_prev)[tuned].max())]}}
    # production at 10:33:10
    G_ = NM.read_geometry(geometry)
    Tt = G_["times"]["t11"]
    cos, sh = NM.geometry_time(geometry, "t11", nf)
    F_h, F_c, _, _ = NM.source_rates(cos, sh, Tt["rh_m"], helio, sza, A, el["ACO2"], kappa, fc.AU)
    P_e = F_h * area
    hold, tot = NS.area_holding(F_h[seen], P_e[seen], area[seen])
    a_seen = float(area[seen].sum())
    for v in hold.values():
        v["fraction_of_seen_area"] = v["area_m2"] / a_seen
    out["seen_production_103310"] = {"holding": hold, "Q_seen_s-1": tot, "Q_all_s-1": float(P_e.sum()),
                                     "tuned_share_of_seen_production": float(P_e[seen & tuned].sum() / tot)}
    # share of the late signal from the tuned facets
    ipoint = np.array([int(r[cols["iPoint"]]) for r in rows])
    TRG = np.array([float(r[cols["T_RG"]]) for r in rows])
    bm, cNG, cRG, _ = NC.gauge_constants(fc, header, TRG)
    from cops_cuda import device_store
    NC._G.update({"A": A, "ACO2": el["ACO2"], "T": T, "kappa": kappa, "bm": bm, "cNG": cNG, "cRG": cRG,
                  "helio": helio, "AU": float(fc.AU), "COS75": float(fc.COS75), "KT": float(fc.KT),
                  "fp32_record_fields": device_store.fp32_record_fields})
    NS._W.update({"base": str(base), "A": A, "active": tuned})
    pt = rec["fit"]["point_table"]["path"]
    lay = pp.LateLayout(pp.read_point_table(pt), pp.read_run_points(base), ipoint)
    W = lay.W.tocsr()
    need = set(int(r) for r in lay.grid_rows[np.unique(W.indices)])
    tasks = sorted([(r, int(ipoint[r])) for r in need], key=lambda q: -int(float(rows[q[0]][cols["nVisible"]])))
    n = len(ipoint)
    tot_ = {k: np.zeros(n) for k in ("ng", "rg", "ng_t", "rg_t")}
    with mp.get_context("fork").Pool(procs) as pool:
        for r, ng, rg, ngt, rgt, bad, nonf, et in pool.imap_unordered(NS.row_work, tasks, chunksize=1):
            if bad or nonf:
                raise SystemExit("negative / non-finite contributions in row %d" % r)
            tot_["ng"][r], tot_["rg"][r], tot_["ng_t"][r], tot_["rg_t"][r] = ng, rg, ngt, rgt
    lr = lay.late_rows
    late = {k: W @ v[lay.grid_rows] for k, v in tot_.items()}
    # verification against late_window_model.dat (surface columns of the late observations)
    lw = {}
    for line in open(ft / "late_window_model.dat"):
        if line.startswith("#"):
            continue
        v = line.split()
        if v[3] == "late_obs":
            lw[int(v[0])] = (float(v[5]), float(v[6]))
    ip_l = ipoint[lr]
    refN = np.array([lw[int(i)][0] for i in ip_l])
    refR = np.array([lw[int(i)][1] for i in ip_l])
    out["verification_late_surface_max_rel"] = {"NG": float(np.max(np.abs(late["ng"] / refN - 1))),
                                                "RG": float(np.max(np.abs(late["rg"] / refR - 1)))}
    sh_ = {}
    for g, k in (("NG", "ng"), ("RG", "rg")):
        m = use[g][lr]
        s = late[k + "_t"][m] / late[k][m]
        sh_[g] = {"rows": int(m.sum()), "mean": float(np.mean(s)), "min": float(np.min(s)), "max": float(np.max(s)),
                  "ratio_of_sums": float(late[k + "_t"][m].sum() / late[k][m].sum())}
    out["tuned_share_of_late_signal"] = sh_
    out["tuned_share_note"] = "surface pressure (H2O + CO2) of the tuned facets / all facets, late observations (used)"
    return out


def text(S, R=None):
    def g(d, *ks):
        for k in ks:
            d = (d or {}).get(k)
        return d

    def f(v, fmt="%.6g"):
        return "-" if v is None else (fmt % v if isinstance(v, (int, float)) else str(v))
    L = []
    for lab, X in (("run", S), ("reference", R)):
        if X is None:
            continue
        cf = X["costs_final"]
        L.append("%s: %s (stop %s, %s evaluations; settings %s)" % (lab, X["ft"], X["stop"], X["nfev_total"],
                                                                  json.dumps(X["settings"], sort_keys=True)))
        L.append("  data cost %s (start %s); P_delta %s; total %s" % (
            f(cf["cost_data"], "%.10e"), f(g(X, "costs_initial", "cost_data"), "%.10e"),
            f(cf.get("cost_delta_penalty"), "%.6e"), f(cf["cost_total"], "%.10e")))
        L.append("  late-window rms ln NG %s RG %s (start NG %s RG %s); all rows NG %s RG %s" % (
            f(cf["late_window_rms_ln_NG"], "%.5f"), f(cf["late_window_rms_ln_RG"], "%.5f"),
            f(g(X, "costs_initial", "late_window_rms_ln_NG"), "%.5f"), f(g(X, "costs_initial", "late_window_rms_ln_RG"), "%.5f"),
            f(cf["rms_ln_NG"], "%.5f"), f(cf["rms_ln_RG"], "%.5f")))
        tl = X["total_log_factor"]
        L.append("  |h + delta| over %d tuned facets: rms %.4f, max %.4f (held before %.4f); delta range [%.4f, %.4f]; "
                 "facets at the cap %s (upper %s, lower %s)" % (
                     X["tuned_facets"], tl["rms_area_weighted"], tl["max_abs"], tl["held_prev_max_abs"],
                     tl["delta_range"][0], tl["delta_range"][1], f(cf.get("facets_at_cap")),
                     f(cf.get("facets_at_cap_upper")), f(cf.get("facets_at_cap_lower"))))
        h = X["seen_production_103310"]["holding"]
        L.append("  seen area holding 50 %% / 90 %% of the seen production at 10:33:10: %.5f / %.5f; tuned facets' share "
                 "of the seen production %.4f" % (h["50%"]["fraction_of_seen_area"], h["90%"]["fraction_of_seen_area"],
                                                  X["seen_production_103310"]["tuned_share_of_seen_production"]))
        s = X["tuned_share_of_late_signal"]
        L.append("  tuned facets' share of the late signal (95 late observations): NG mean %.4f (min %.4f, max %.4f), RG "
                 "mean %.4f (min %.4f, max %.4f); verification vs late_window_model.dat %s" % (
                     s["NG"]["mean"], s["NG"]["min"], s["NG"]["max"], s["RG"]["mean"], s["RG"]["min"], s["RG"]["max"],
                     json.dumps(X["verification_late_surface_max_rel"])))
    return "\n".join(L) + "\n"


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--ft", required=True)
    ap.add_argument("--reference-ft", default=None)
    ap.add_argument("--base", required=True)
    ap.add_argument("--geometry", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--procs", type=int, default=int(os.environ.get("SLURM_CPUS_PER_TASK", "8")))
    args = ap.parse_args(argv)
    t0 = time.time()
    out = Path(args.out)
    if out.exists():
        raise SystemExit("%s exists" % out)
    base = Path(args.base).resolve()
    gval = json.load(open(Path(args.geometry) / "geometry_validation.json"))
    if not gval.get("passed") or Path(gval["base"]) != base:
        raise SystemExit("geometry not validated against %s" % base)
    fc, header, cols, rows = NM.base_table(base)
    helio, sza = fc.header_value(header, "HeliocentricModel"), fc.header_value(header, "SZAMode")
    el = NM.read_elements(base / "elements.bin")
    S = one(args.ft, base, args.geometry, args.procs, el, fc, header, cols, rows, helio, sza)
    R = one(args.reference_ft, base, args.geometry, args.procs, el, fc, header, cols, rows, helio, sza) \
        if args.reference_ft else None
    doc = {"schema": SCHEMA, "created_utc": NM.utc_now(), "host": socket.gethostname(),
           "slurm_job_id": os.environ.get("SLURM_JOB_ID"), "command": sys.argv, "run": S, "reference": R,
           "code_sha256": {p.name: NM.sha256_file(p) for p in (Path(__file__).resolve(), Path(NS.__file__).resolve(),
                                                               Path(NC.__file__).resolve(), Path(NM.__file__).resolve())},
           "wall_seconds": time.time() - t0}
    with open(out, "x") as f:
        f.write("# %s -- nonneg_ft_summary.py (%s, job %s)\n" % (out.name, doc["created_utc"], doc["slurm_job_id"]))
        f.write(text(S, R))
    with open(str(out) + ".json", "x") as f:
        json.dump(doc, f, indent=1, sort_keys=True, default=float)
    print(text(S, R), flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
