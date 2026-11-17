#!/usr/bin/env python3
# ==== version-controlled copy: header added 2026-11-15, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/model_ext/localnn_chain_summary.py (on Voyager: /data001/highresolu_vt/cops_output/model_ext/localnn_chain_summary.py)
# SHA256 of the run copy = SHA256 of this file without these header lines: 0e449b7eadeae9c1741c33fd5146148e4a470c2a77ccf382264b2a6756fcea45
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: per-cycle summary of option C stage-2 chains (user / lead 2026-11-16): fit and fine-tuning costs (data,
#   smoothness, P_m, P_delta, objective, late rms ln NG / RG), multipliers, Q_H2O(10:33:10), active seen area, seen
#   area holding 50 % / 90 % of the seen production, the facet with the largest activity and its share of the seen
#   production.
# RUN 2026-11-15: 135078 (the unregularized f10 / 0.01 chain); the stage-2r chains.
# ==== end of header ====
"""Per-cycle summary of option C stage-2 chains (user / lead 2026-11-16, the stage-2 reruns templates on / off):

    localnn_chain_summary.py --chain LABEL=CHAIN [--chain ...] --base BASE --geometry GEOM --contrib-code CC
        --out summary.json [--txt summary.txt]

For every chain and cycle k, the fit (c<k>/fit/final) and the fine tuning (c<k>/ft/final):
  fit: data cost, crossing, smoothness penalty, P_m, objective, late-window rms ln NG / RG, evaluations, termination,
       the 7 multipliers (localnn_state.npz meta / m);
  ft:  data cost, P_delta (cost_delta_penalty; 0 without regularization), late-window rms ln NG / RG, evaluations,
       termination (fine_tuning.json costs_final); objective after the fine tuning = ft data + crossing + the fit's
       smoothness + the fit's P_m (+ P_delta reported apart: the next fit does not carry it);
  both states: Q_H2O(10:33:10) and the seen production, the active seen area fraction (A > 0.01 x median AH2O over the
       seen facets), the seen area holding 50 % / 90 % of the seen production (localnn_report: facets in decreasing
       production), the facet with the largest activity A: index, slope bin, area, A, its share of the seen production
       (and the facet with the largest production and its share); validated maps geometry t11 (10:33:10 UTC), the srcCG
       production law (nonneg_maps.source_rates, the fit's kappa).
Light CPU job (one geometry read, one pass per state).
"""

import argparse
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))
import localnn_problem as lp      # noqa: E402
import localnn_design as ld       # noqa: E402
import localnn_report as lr       # noqa: E402


def stats(A, kappa, seen, bins, S):
    F_h, _, _, _ = S["NM"].source_rates(S["cos"], S["sh"], S["rh"], S["helio"], S["sza"], A, S["el"]["ACO2"], kappa,
                                        S["AU"])
    area = S["el"]["area"]
    P_e = F_h * area
    Ps = float(np.sum(P_e[seen]))
    sa = float(np.sum(area[seen]))
    hold = lr.area_holding(P_e[seen], area[seen], (0.5, 0.9))
    ia, ip = int(np.argmax(A)), int(np.argmax(np.where(seen, P_e, -np.inf)))
    return {"Q_H2O_all_s-1": float(np.sum(P_e)), "Q_H2O_seen_s-1": Ps,
            "active_seen_area_fraction": float(np.sum(area[seen & (A > S["thr"])]) / sa),
            "seen_area_holding_production_m2": hold,
            "seen_area_holding_production_fraction": {k: v / sa for k, v in hold.items()},
            "largest_activity_facet": {"facet": ia, "bin": lp.BIN_NAMES[int(bins[ia])], "area_m2": float(area[ia]),
                                       "A": float(A[ia]), "seen": bool(seen[ia]),
                                       "share_of_seen_production": float(P_e[ia] / Ps) if seen[ia] else None},
            "largest_production_facet": {"facet": ip, "bin": lp.BIN_NAMES[int(bins[ip])], "area_m2": float(area[ip]),
                                         "share_of_seen_production": float(P_e[ip] / Ps)}}


def main(argv=None):
    ap = argparse.ArgumentParser()
    for k in ("base", "geometry", "contrib-code", "out"):
        ap.add_argument("--" + k, required=True)
    ap.add_argument("--chain", action="append", required=True)
    ap.add_argument("--txt")
    args = ap.parse_args(argv)
    CS = ld.contrib_setup(args.base, args.contrib_code, lp.USER_T0)
    el, NM, fc = CS["el"], CS["NM"], CS["fc"]
    G = NM.read_geometry(args.geometry)
    Tt = G["times"][lr.GEOM_LABEL]
    cos, sh = NM.geometry_time(args.geometry, lr.GEOM_LABEL, el["nfaces"])
    S = {"NM": NM, "el": el, "cos": cos, "sh": sh, "rh": Tt["rh_m"], "AU": float(fc.AU),
         "helio": fc.header_value(CS["header"], "HeliocentricModel"), "sza": fc.header_value(CS["header"], "SZAMode")}
    res = {"geometry": {"label": lr.GEOM_LABEL, "utc": Tt["utc"]}, "chains": {}}
    for spec in args.chain:
        label, ch = spec.split("=", 1)
        ch = Path(ch).resolve()
        cyc = []
        for k in range(1, 100):
            fd = ch / ("c%d" % k) / "fit" / "final"
            if not fd.exists():
                break
            fd = fd.resolve()
            with np.load(fd / "localnn_state.npz", allow_pickle=False) as z:
                meta = json.loads(str(z["meta_json"]))
                m, kappa = [float(v) for v in z["m"]], float(z["kappa"])
            with np.load(fd / "nonneg_facets.npz", allow_pickle=False) as z:
                A, seen, bins = np.asarray(z["A"], np.float64), np.asarray(z["seen"], bool), np.asarray(z["bin"])
            if "thr" not in S:
                S["thr"] = lr.ACTIVE_FRACTION * float(np.median(el["AH2O"][seen]))
            c = meta["costs"]
            row = {"cycle": k, "fit": {
                "dir": str(fd), "templates": meta.get("templates"), "nfev": meta.get("nfev_total"),
                "termination": meta.get("message"), "data_cost": c["data_cost"], "crossing": c["crossing_cost"],
                "smoothness": c["penalty"], "P_m": c.get("mult_penalty", 0.0), "objective": c["objective"],
                "late_rms_ln_NG": c.get("late_rms_ln_NG"), "late_rms_ln_RG": c.get("late_rms_ln_RG"),
                "multipliers": m, "kappa": kappa, "state": stats(A, kappa, seen, bins, S)}}
            ftd = ch / ("c%d" % k) / "ft" / "final"
            if ftd.exists():
                ftd = ftd.resolve()
                j = json.load(open(ftd / "fine_tuning.json"))
                cf = j["costs_final"]
                with np.load(ftd / "finetune_state.npz", allow_pickle=False) as z:
                    Af, hn, tuned = np.asarray(z["A_final"], np.float64), np.asarray(z["held_next"]), np.asarray(z["seen"], bool)
                pd = float(cf.get("cost_delta_penalty") or 0.0)
                row["ft"] = {"dir": str(ftd), "nfev": j["optimizer"].get("nfev"),
                             "termination": j["optimizer"].get("termination"), "data_cost": cf["cost_data"],
                             "crossing": cf.get("cost_crossing", 0.0), "P_delta": pd,
                             "objective_after_ft": cf["cost_data"] + cf.get("cost_crossing", 0.0) + c["penalty"] +
                             c.get("mult_penalty", 0.0),
                             "late_rms_ln_NG": cf.get("late_window_rms_ln_NG"),
                             "late_rms_ln_RG": cf.get("late_window_rms_ln_RG"),
                             "max_abs_h_next_tuned": float(np.abs(hn[tuned]).max()) if tuned.any() else 0.0,
                             "facets_at_cap": cf.get("facets_at_cap"),
                             "state": stats(Af, kappa, seen, bins, S)}
            cyc.append(row)
            print("%s c%d done" % (label, k), flush=True)
        res["chains"][label] = {"chain": str(ch), "cycles": cyc}
    with open(args.out, "x") as f:
        json.dump(res, f, indent=1, default=float)
    if args.txt:
        L = ["# per chain and cycle (production at %s = %s UTC); fit row then ft row; objective = data + crossing + "
             "smoothness + P_m (ft: with the fit's smoothness and P_m; P_delta apart)" % (lr.GEOM_LABEL, Tt["utc"]),
             "#chain cycle stage nfev stop data smooth P_m P_delta objective lateNG lateRG Q_H2O active seen50_m2 "
             "seen90_m2 maxA_facet maxA_bin maxA_share maxP_facet maxP_share multipliers(40-50..90-180)"]
        for label, cc in res["chains"].items():
            for r in cc["cycles"]:
                for stg in ("fit", "ft"):
                    if stg not in r:
                        continue
                    x, s = r[stg], r[stg]["state"]
                    f = r["fit"]
                    L.append("%s c%d %s %s %s %.6e %.4e %.4e %s %.6e %.4f %.4f %.3e %.3f %.3e %.3e %d %s %s %d %.4f %s" % (
                        label, r["cycle"], stg, x["nfev"], str(x["termination"]).replace(" ", "_")[:24], x["data_cost"],
                        f["smoothness"], f["P_m"], "%.3e" % x["P_delta"] if stg == "ft" else "-",
                        x["objective"] if stg == "fit" else x["objective_after_ft"],
                        x["late_rms_ln_NG"], x["late_rms_ln_RG"], s["Q_H2O_all_s-1"], s["active_seen_area_fraction"],
                        s["seen_area_holding_production_m2"]["50%"], s["seen_area_holding_production_m2"]["90%"],
                        s["largest_activity_facet"]["facet"], s["largest_activity_facet"]["bin"],
                        "%.4f" % s["largest_activity_facet"]["share_of_seen_production"]
                        if s["largest_activity_facet"]["share_of_seen_production"] is not None else "unseen",
                        s["largest_production_facet"]["facet"], s["largest_production_facet"]["share_of_seen_production"],
                        ",".join("%.4g" % v for v in f["multipliers"][1:])))
        Path(args.txt).write_text("\n".join(L) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
