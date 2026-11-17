#!/usr/bin/env python3
# ==== version-controlled copy: header added 2026-11-15, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/model_ext/localnn_ftreg_summary.py (on Voyager: /data001/highresolu_vt/cops_output/model_ext/localnn_ftreg_summary.py)
# SHA256 of the run copy = SHA256 of this file without these header lines: 38fa5335dd39ec8d49b5de936c1097af76031a288594919c83f837a39a25f56b
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: per-run numbers of the fine-tuning regularization tests (lead / user 2026-11-15 ~21:10 CST): data cost,
#   late-window rms ln NG / RG, P_delta, max |h+delta| and facets at the cap, the seen area holding 50 % / 90 % of
#   the seen production and Q_H2O at 10:33:10 (localnn_report.production_stats), the share of the late-window
#   surface pressure from the tuned facets.
# RUN 2026-11-15: reference 134496; the test runs.
# ==== end of header ====
"""Fine-tuning regularization tests (lead / user 2026-11-15 ~21:10 CST): the per-run numbers of the c5 fine tunings
with P_delta = f_delta C_ref sum w (h+delta)^2 / sum w and the cap |h+delta| <= ln 100, against the unregularized c5
fine tuning of the same fit (agent 1's nonneg_finetune options; localnn_ft.py run / run-nonneg on CPU).

    localnn_ftreg_summary.py --late LK --base BASE --d2 D2 --point-table PT --design DESIGN --geometry GEOM
        --contrib-code CC --case LABEL=FT_DIR [--case ...] --out summary.json [--txt summary.txt]

Per case (FT_DIR = a fine-tuning output folder; its fit from fine_tuning.json):
  data cost, late-window rms ln NG / RG, cost_delta_penalty, total, evaluations, termination (fine_tuning.json
  costs_final / optimizer); max |h + delta| over the tuned facets and the facets at the cap (finetune_state.npz held_next,
  seen = the tuned (parameter) facets); the seen area holding 50 % / 90 % of the seen H2O production at 10:33:10 UTC and
  Q_H2O (localnn_report.production_stats on A_final: the validated maps geometry t11, srcCG production law, the fit's
  kappa, seen = the fit's nonneg_facets.npz seen); the share of the late-window model surface pressure (W10x5, 95 late
  observations) from the tuned facets: E (G[:, tuned] A_final[tuned]) / surface_p of the fine tuning's
  late_window_model.dat (mean, min, max, at the last observation), NG and RG. Light CPU job (one pass over the facets).
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
import localnn_late as llate      # noqa: E402
import localnn_report as lr       # noqa: E402


def main(argv=None):
    ap = argparse.ArgumentParser()
    for k in ("late", "base", "d2", "point-table", "design", "geometry", "contrib-code", "out"):
        ap.add_argument("--" + k, required=True)
    ap.add_argument("--case", action="append", required=True)
    ap.add_argument("--txt")
    args = ap.parse_args(argv)
    design = lp.load_design(args.design)
    D = lp.Data(args.base, args.d2, args.point_table, design)
    CS = ld.contrib_setup(args.base, args.contrib_code, lp.USER_T0)
    el, NM, fc = CS["el"], CS["NM"], CS["fc"]
    helio, sza = fc.header_value(CS["header"], "HeliocentricModel"), fc.header_value(CS["header"], "SZAMode")
    nf = el["nfaces"]
    LK = llate.load(args.late, nf)
    pos = np.full(nf, -1, np.int64)
    pos[LK["late_idx"]] = np.arange(len(LK["late_idx"]))
    late_rows = np.asarray(D.late_rows)
    res = {"cases": {}, "late_kernels": str(Path(args.late).resolve()), "geometry_label": lr.GEOM_LABEL}
    for case in args.case:
        label, ftd = case.split("=", 1)
        ftd = Path(ftd).resolve()
        ftj = json.load(open(ftd / "fine_tuning.json"))
        fit = Path(ftj["fit"]["directory"]).resolve()
        cf, ci, opt = ftj["costs_final"], ftj["costs_initial"], ftj["optimizer"]
        with np.load(ftd / "finetune_state.npz", allow_pickle=False) as z:
            A, hn, delta, tuned = (np.asarray(z[k]) for k in ("A_final", "held_next", "delta", "seen"))
        tuned = tuned.astype(bool)
        with np.load(fit / "nonneg_facets.npz", allow_pickle=False) as z:
            seen = np.asarray(z["seen"], bool)
        if (fit / "localnn_state.npz").exists():
            with np.load(fit / "localnn_state.npz", allow_pickle=False) as z:
                kappa = float(z["kappa"])
        else:
            with np.load(fit / "fit_restart.npz", allow_pickle=False) as z:
                kappa = float(np.asarray(z["kappa"]).reshape(-1)[0])
        cap = cf.get("delta_cap") or ftj.get("settings", {}).get("delta_cap")
        ht = np.abs(hn[tuned])
        thr = lr.ACTIVE_FRACTION * float(np.median(el["AH2O"][seen]))
        prod = lr.production_stats(A, kappa, el, seen, args.geometry, NM, helio, sza, float(fc.AU), thr)
        # late-window surface pressure from the tuned facets
        p = pos[np.flatnonzero(tuned)]
        if np.any(p < 0):
            raise SystemExit("%s: tuned facets outside the late kernels" % label)
        At = A[np.flatnonzero(tuned)]
        tN, tR = D.E @ (LK["G_NG"][:, p] @ At), D.E @ (LK["G_RG"][:, p] @ At)
        lw = lp.read_table_cols(ftd / "late_window_model.dat", ("kind", "surface_pNG_Pa", "surface_pRG_Pa"))
        share = {}
        for g, t, col in (("NG", tN, "surface_pNG_Pa"), ("RG", tR, "surface_pRG_Pa")):
            s = np.array([t[r] / float(lw[int(D.ipoint[r])][col]) for r in late_rows])
            share[g] = {"mean": float(s.mean()), "min": float(s.min()), "max": float(s.max()), "last": float(s[-1])}
        res["cases"][label] = {
            "ft_dir": str(ftd), "fit": str(fit), "settings": {k: ftj.get("settings", {}).get(k) for k in (
                "delta_penalty_f", "delta_cref", "delta_cap", "max_nfev", "facets", "window_min")},
            "nfev": opt.get("nfev"), "termination": opt.get("termination"),
            "data_cost": cf.get("cost_data"), "data_cost_initial": ci.get("cost_data"),
            "cost_total": cf.get("cost_total"), "cost_delta_penalty": cf.get("cost_delta_penalty"),
            "late_rms_ln_NG": cf.get("late_window_rms_ln_NG"), "late_rms_ln_RG": cf.get("late_window_rms_ln_RG"),
            "rms_ln_NG": cf.get("rms_ln_NG"), "rms_ln_RG": cf.get("rms_ln_RG"),
            "tuned_facets": int(tuned.sum()), "max_abs_h_plus_delta": float(ht.max()) if ht.size else 0.0,
            "facets_at_cap": int(np.sum(ht >= float(cap) * (1 - 1e-12))) if cap else None,
            "facets_at_cap_recorded": cf.get("facets_at_cap"), "cap": cap,
            "delta_total_rms_area_weighted": cf.get("delta_total_rms_area_weighted"),
            "max_abs_delta": float(np.abs(delta[tuned]).max()) if tuned.any() else 0.0,
            "production_103310": prod, "late_signal_share_tuned_facets": share}
        print("%s: %s" % (label, json.dumps({k: res["cases"][label][k] for k in (
            "data_cost", "late_rms_ln_NG", "late_rms_ln_RG", "cost_delta_penalty", "max_abs_h_plus_delta",
            "facets_at_cap")}, default=float)), flush=True)
    with open(args.out, "x") as f:
        json.dump(res, f, indent=1, default=float)
    if args.txt:
        L = ["# fine-tuning regularization tests: per run (production at %s; share = late-window surface pressure from "
             "the tuned facets, mean over the 95 late observations)" % lr.GEOM_LABEL,
             "#label f_delta cap nfev stop data_cost late_rms_NG late_rms_RG P_delta max|h+d| at_cap tuned "
             "seen50%_m2 seen90%_m2 Q_H2O_s-1 shareNG shareRG shareNG_last shareRG_last"]
        for label, c in res["cases"].items():
            pr = c["production_103310"]
            L.append("%s %s %s %s %s %.6e %.4f %.4f %s %.3f %s %d %.4e %.4e %.4e %.4f %.4f %.4f %.4f" % (
                label, c["settings"]["delta_penalty_f"], c["cap"], c["nfev"], c["termination"], c["data_cost"],
                c["late_rms_ln_NG"], c["late_rms_ln_RG"],
                "%.4e" % c["cost_delta_penalty"] if c["cost_delta_penalty"] is not None else "-",
                c["max_abs_h_plus_delta"], c["facets_at_cap"], c["tuned_facets"],
                pr["seen_area_holding_production_m2"]["50%"], pr["seen_area_holding_production_m2"]["90%"],
                pr["Q_H2O_all_s-1"], c["late_signal_share_tuned_facets"]["NG"]["mean"],
                c["late_signal_share_tuned_facets"]["RG"]["mean"], c["late_signal_share_tuned_facets"]["NG"]["last"],
                c["late_signal_share_tuned_facets"]["RG"]["last"]))
        Path(args.txt).write_text("\n".join(L) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
