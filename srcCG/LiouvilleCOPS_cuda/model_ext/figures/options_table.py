#!/usr/bin/env python3
# ==== version-controlled copy: header added 2026-11-13, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/figures_cluster/scripts/options_table.py (on Voyager: /data001/highresolu_vt/cops_output/figures_cluster/scripts/options_table.py)
# SHA256 of the run copy = SHA256 of this file without these header lines: 2d9af51ce8b46af4bd1fc785b2dab46d39a6ceb0da19935edbea1892c0363e0c
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: table of the stage-1 option fits and references (summary.txt / .csv / best.json) from each run's
#   summary.txt.json or analysis/report.json (2026-11-14); --diag; 2026-11-15: prev_* columns of continuations, C
#   multiplier penalty, --conv-from-termination, setting column widened only for longer settings
# RUN: called by figures_12m_nonneg_options_stage1.sh: options_table.py --spec <folder>/runs.json --out-dir <folder>
#   [--diag <json> | --conv-from-termination]
# ==== end of header ====
"""Table of the stage-1 fits of the non-negative activity options and the references (2026-11-14, user: "make
comparison plot of the results with the data and the current softplus fit"; lead: one row per fit, read from each run's
summary.txt.json (options A, B, D: model_ext/nonneg_alt_summary.py) or analysis/report.json (option C:
model_ext/localnn_report.py), not recomputed).

    options_table.py --spec <runs.json> --out-dir <dir>
runs.json: [{"option": "B", "name": "background", "setting": "bg 1e-2", "root": <run root>, "kind": "summary"|"localnn"}, ...]
References: the bounded softplus c9 (finetune_result.txt: costs, rms overall and late) and D2 (fit_result.txt: cost, rms
overall); their Q_H2O, active seen area and production areas from option C's report.json (production_1033.c9 / .D2, the same
definitions as the summaries). Writes summary.txt, summary.csv, best.json (the lowest-data-cost setting per option).
Option --diag <json> (2026-11-14, lead: convergence column from agent 1's CPU diagnosis): the verdict per fit ("option|setting"
-> code, note, diag files) and, read from the listed diag_*.json, the scaled gradient, the Gauss-Newton predicted relative
decrease and the best trial step; adds the columns DIAG_COLS to summary.csv, a "conv" column and a convergence section to
summary.txt. Without --diag the output is unchanged.
2026-11-15 (final stage-1 comparison): (1) a summary-kind run started from another fit's restart (summary "start" with
start_choice "restart", e.g. the stage-1 continuations) gets that source fit's data cost, penalty, nfev and termination as
prev_* columns (from the source run's summary.txt.json) and a "before the continuation" section; (2) an option-C report
with a multiplier penalty (costs.mult_penalty) gets penalty = smoothness + multiplier penalty and the column penalty_mult;
(3) --conv-from-termination fills the convergence column from each fit's termination (ftol: converged in the optimizer's
sense; max_nfev: not converged). Columns (1) and (2) appear only when such runs are present; without them and without the
options the output is unchanged.
"""
import argparse
import csv
import json
import re
import sys
from pathlib import Path

MULT = ("c_slope_40_50", "c_slope_50_60", "c_slope_60_70", "c_slope_70_80", "c_slope_80_90", "c_slope_90_180")
COLS = ["option", "model", "setting", "data_cost", "penalty", "objective", "termination", "nfev",
        "rms_NG_all", "rms_RG_all", "rms_NG_before_gap", "rms_RG_before_gap", "rms_NG_late", "rms_RG_late",
        "at_bounds", "mult_40_50", "mult_50_60", "mult_60_70", "mult_70_80", "mult_80_90", "mult_90_180",
        "Q_H2O_103310_s-1", "active_seen_area_frac", "seen_area_50pct_prod_frac", "seen_area_90pct_prod_frac",
        "NG_share_active_late_mean", "RG_share_active_late_mean", "NG_share_active_all_mean", "RG_share_active_all_mean",
        "source"]
PREV_COLS = ["prev_data_cost", "prev_penalty", "prev_nfev", "prev_termination", "prev_source"]
DIAG_COLS = ["convergence", "conv_note", "diag_g_scaled_inf", "diag_GN_pred_rel", "diag_best_trial_rel", "diag_best_trial",
             "diag_files"]


def diag_values(spec, r, roots):
    """Verdict + numbers of one fit from the --diag spec (None: fit not listed)."""
    v = spec["fits"].get("%s|%s" % (r["option"], r["setting"]))
    if v is None:
        return None
    out = {"convergence": v["code"], "conv_note": v["note"], "diag_files": " ".join(v["files"])}
    best = None
    for k, fn in enumerate(v["files"]):
        d = json.load(open(Path(spec["diag_dir"]) / fn))
        if not str(d["fit"]).startswith(roots[(r["option"], r["setting"])] + "/"):
            raise SystemExit("%s: diag fit %s is not in the run root %s" % (fn, d["fit"], roots[(r["option"], r["setting"])]))
        if k == 0:
            out["diag_g_scaled_inf"] = d["gradient"]["g_inf_free_scaled"]
            out["diag_GN_pred_rel"] = d["gauss_newton"]["predicted_relative"]
        for t in d.get("trials", []):
            x = t.get("relative_actual_decrease")
            if x is not None and (best is None or x > best[0]):
                best = (x, "%s (%s)" % (t["trial"], fn))
    if best:
        out["diag_best_trial_rel"], out["diag_best_trial"] = best
    return out


def from_summary(e):
    root = Path(e["root"])
    d = json.load(open(root / "summary.txt.json"))
    c, t, r, p = d["costs"], d["termination"], d["rms_ln"]["values"], d["parameters"]
    sp, sh = d["seen_production"]["at_103310"]["holding"], d["shares_from_active_area"]["values"]
    row = {"data_cost": c["data_cost"], "penalty": c.get("spread_penalty") or 0.0, "objective": c["objective_with_crossing"],
           "termination": t["message"], "nfev": t["nfev"],
           "rms_NG_all": r["NG"]["all"][0], "rms_RG_all": r["RG"]["all"][0],
           "rms_NG_before_gap": r["NG"]["before_gap"][0], "rms_RG_before_gap": r["RG"]["before_gap"][0],
           "rms_NG_late": r["NG"]["late_window"][0], "rms_RG_late": r["RG"]["late_window"][0],
           "at_bounds": " ".join(d["at_bounds_all"]),
           "Q_H2O_103310_s-1": d["Q_H2O_103310"]["Q_all_s-1"],
           "active_seen_area_frac": d["active_area"]["active_seen_area_fraction"],
           "seen_area_50pct_prod_frac": sp["50%"]["fraction_of_seen_area"],
           "seen_area_90pct_prod_frac": sp["90%"]["fraction_of_seen_area"],
           "NG_share_active_late_mean": sh["NG"]["late_window"]["mean"], "RG_share_active_late_mean": sh["RG"]["late_window"]["mean"],
           "NG_share_active_all_mean": sh["NG"]["all"]["mean"], "RG_share_active_all_mean": sh["RG"]["all"]["mean"],
           "source": str(root / "summary.txt.json")}
    for m in MULT:
        row["mult_" + m[8:]] = p[m]["x"] if m in p else None
    st = d.get("start") or {}
    if st.get("start_choice") == "restart" and st.get("source_fit_directory"):
        src = Path(st["source_fit_directory"]).parent.parent
        ds = json.load(open(src / "summary.txt.json"))
        row.update({"prev_data_cost": ds["costs"]["data_cost"], "prev_penalty": ds["costs"].get("spread_penalty") or 0.0,
                    "prev_nfev": ds["termination"]["nfev"], "prev_termination": ds["termination"]["message"],
                    "prev_source": str(src / "summary.txt.json")})
    return row


def from_localnn(e):
    root = Path(e["root"])
    d = json.load(open(root / "analysis" / "report.json"))
    c, t, r, b = d["costs"], d["termination"], d["rms_ln"], d["bounds"]
    pr = d["production_1033"]["localnn"]
    sh = d["active_area_shares"]
    ab = [("c_j at 0: %d nodes (%.3f of all)" % (b["c_at_0_nodes"], b["c_at_0_fraction_all_nodes"]))]
    ab += [k.replace("_at_0", "") for k in ("kappa_at_0", "grain_A_at_0", "grain_B_at_0") if b.get(k)]
    ab += ["mult at 1: %s" % ",".join(map(str, b.get("multipliers_at_1", [])))] if b.get("multipliers_at_1") else []
    row = {"data_cost": c["data_cost"], "penalty": c["penalty"] + c.get("mult_penalty", 0.0), "objective": c["objective"],
           "termination": t["message"], "nfev": t.get("nfev_total", t.get("nfev")),
           "rms_NG_all": r["overall"][0], "rms_RG_all": r["overall"][1], "rms_NG_before_gap": r["before_gap"][0],
           "rms_RG_before_gap": r["before_gap"][1], "rms_NG_late": r["late_window"][0], "rms_RG_late": r["late_window"][1],
           "at_bounds": "; ".join(ab), "Q_H2O_103310_s-1": pr["Q_H2O_all_s-1"],
           "active_seen_area_frac": pr["active_seen_area_fraction"],
           "seen_area_50pct_prod_frac": pr["seen_area_holding_production_fraction"]["50%"],
           "seen_area_90pct_prod_frac": pr["seen_area_holding_production_fraction"]["90%"],
           "NG_share_active_late_mean": sh["late_window_95"]["NG"]["mean"], "RG_share_active_late_mean": sh["late_window_95"]["RG"]["mean"],
           "NG_share_active_all_mean": sh["all_used_observations"]["NG"]["mean"],
           "RG_share_active_all_mean": sh["all_used_observations"]["RG"]["mean"],
           "source": str(root / "analysis" / "report.json")}
    if "mult_penalty" in c:
        row["penalty_mult"] = c["mult_penalty"]
    mu = d["multipliers"]
    for i, m in enumerate(MULT):
        row["mult_" + m[8:]] = mu[i + 1] if len(mu) == 7 else None
    return row, d


def refs(spec_refs, rep):
    out = []
    for e in spec_refs:
        row = {"option": e["option"], "model": e["name"], "setting": e["setting"], "source": e["root"]}
        txt = Path(e["result"]).read_text()
        if e["kind"] == "finetune":
            g = lambda k: float(re.search(r"^%s\s+initial \S+\s+final (\S+)" % k, txt, re.M).group(1))   # noqa: E731
            row.update({"data_cost": g("cost_data"), "penalty": 0.0, "objective": g("cost_total"),
                        "rms_NG_all": g("rms_ln_NG"), "rms_RG_all": g("rms_ln_RG"), "rms_NG_late": g("late_window_rms_ln_NG"),
                        "rms_RG_late": g("late_window_rms_ln_RG"),
                        "termination": re.search(r"^message (.*)$", txt, re.M).group(1) + " (fine tuning)",
                        "nfev": int(re.search(r"^nfev_total (\d+)", txt, re.M).group(1))})
        else:
            row.update({"data_cost": float(re.search(r"^final_data_cost (\S+)", txt, re.M).group(1)), "penalty": 0.0,
                        "objective": float(re.search(r"^final_cost (\S+)", txt, re.M).group(1)),
                        "rms_NG_all": float(re.search(r"^rms_ln_residual_NG (\S+)", txt, re.M).group(1)),
                        "rms_RG_all": float(re.search(r"^rms_ln_residual_RG (\S+)", txt, re.M).group(1)),
                        "termination": re.search(r"^optimizer: (.*)$", txt, re.M).group(1)})
        if rep is not None and e["report_key"] in rep["production_1033"]:
            pr = rep["production_1033"][e["report_key"]]
            row.update({"Q_H2O_103310_s-1": pr["Q_H2O_all_s-1"], "active_seen_area_frac": pr["active_seen_area_fraction"],
                        "seen_area_50pct_prod_frac": pr["seen_area_holding_production_fraction"]["50%"],
                        "seen_area_90pct_prod_frac": pr["seen_area_holding_production_fraction"]["90%"]})
            row["source"] += "; production: %s production_1033.%s" % (rep["_path"], e["report_key"])
        for k, v in e.get("multipliers", {}).items():
            row["mult_" + k] = v
        out.append(row)
    return out


def fmt(v):
    if v is None or v == "":
        return "-"
    if isinstance(v, float):
        return "%.4g" % v
    return str(v)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--spec", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--diag", default=None, help="convergence verdicts + diag files (JSON); default: no convergence column")
    ap.add_argument("--conv-from-termination", action="store_true",
                    help="convergence column from each fit's termination (ftol / max_nfev); default: none")
    a = ap.parse_args()
    spec = json.load(open(a.spec))
    out = Path(a.out_dir)
    rows, rep = [], None
    for e in spec["runs"]:
        if e["kind"] == "summary":
            r = from_summary(e)
        else:
            r, rp = from_localnn(e)
            if rep is None:
                rep = dict(rp, _path=str(Path(e["root"]) / "analysis" / "report.json"))
        r.update({"option": e["option"], "model": e["name"], "setting": e["setting"]})
        rows.append(r)
    rows += refs(spec["refs"], rep)
    cols, dspec = COLS, None
    if a.diag:
        dspec = json.load(open(a.diag))
        roots = {(e["option"], e["setting"]): str(Path(e["root"])) for e in spec["runs"]}
        for r in rows:
            if r["option"] != "ref":
                r.update(diag_values(dspec, r, roots) or {"convergence": "not listed"})
        cols = COLS + DIAG_COLS
    if a.conv_from_termination:
        if a.diag:
            raise SystemExit("--diag and --conv-from-termination exclude each other")
        dspec = {"_about": "from each fit's termination message: ftol = converged in the optimizer's sense (relative cost "
                           "change below ftol 1e-10), max_nfev = not converged; nothing else is implied"}
        for r in rows:
            if r["option"] == "ref":
                continue
            msg = r.get("termination") or ""
            if "ftol" in msg:
                r.update({"convergence": "ftol", "conv_note": "converged in the optimizer's sense (ftol)"})
            elif "maximum number of function evaluations" in msg:
                r.update({"convergence": "not converged", "conv_note": "stopped at max_nfev"})
            else:
                r.update({"convergence": "other", "conv_note": msg})
        cols = COLS + DIAG_COLS[:2]
    if any("penalty_mult" in r for r in rows):
        cols = cols[:cols.index("objective")] + ["penalty_mult"] + cols[cols.index("objective"):]
    if any("prev_data_cost" in r for r in rows):
        cols = cols + PREV_COLS
    best = {}
    for r in rows:
        if r["option"] in ("ref",):
            continue
        if r["option"] not in best or r["data_cost"] < best[r["option"]]["data_cost"]:
            best[r["option"]] = r
    for f in ("summary.txt", "summary.csv", "best.json"):
        if (out / f).exists():
            raise SystemExit("%s exists" % (out / f))
    with open(out / "summary.csv", "x", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols, extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, "") for k in cols})
    L = ["# Stage-1 fits of the non-negative activity options (single fits from D2; slope multipliers >= 1; D2 fine-tuning "
         "templates on) and the references; made by options_table.py from each run's summary.txt.json / analysis/report.json",
         "# data cost = 0.5 r.r of the data rows (the fit's cost without penalties); penalty: option D spread penalty, option C "
         "smoothness penalty; rms ln over the used rows (W10x5 late, + grain for NG); before the gap = UTC < 10:33:10",
         "# active seen area = A > 1 % of the srcCG median of A_orig over the seen facets; 50/90 % = seen area holding 50/90 % "
         "of the seen H2O production at 10:33:10 UTC; shares = surface pressure of the active facets / all facets (means)", ""]
    sw = max([10] + [len(r["setting"]) for r in rows])     # 2026-11-15: wider only for longer settings
    rowfmt = "%-3s %-14s %-" + str(sw) + "s %-10s %-9s %-6s %-23s %-23s %-23s %-9s %-7s %-7s %-7s %-7s %-7s"
    hdr = rowfmt % (
        "opt", "model", "setting", "data", "penalty", "nfev", "rms NG all/bgap/late", "rms RG all/bgap/late",
        "multipliers 40..180", "Q_H2O", "active", "A50%", "A90%", "NGact", "RGact")
    if dspec:
        hdr += " conv"
    L.append(hdr)
    for r in rows:
        mults = "/".join(fmt(r.get("mult_" + m[8:])) for m in MULT)
        L.append(rowfmt % (
            r["option"], r["model"][:14], r["setting"][:sw], fmt(r.get("data_cost")), fmt(r.get("penalty")), fmt(r.get("nfev")),
            "/".join(fmt(r.get(k)) for k in ("rms_NG_all", "rms_NG_before_gap", "rms_NG_late")),
            "/".join(fmt(r.get(k)) for k in ("rms_RG_all", "rms_RG_before_gap", "rms_RG_late")),
            mults[:23], fmt(r.get("Q_H2O_103310_s-1")), fmt(r.get("active_seen_area_frac")),
            fmt(r.get("seen_area_50pct_prod_frac")), fmt(r.get("seen_area_90pct_prod_frac")),
            fmt(r.get("NG_share_active_late_mean")), fmt(r.get("RG_share_active_late_mean"))) +
                 ((" " + (r.get("convergence") or "-")) if dspec else ""))
    L += ["", "multipliers (40-50, 50-60, 60-70, 70-80, 80-90, 90-180 deg) and parameters at bounds:"]
    for r in rows:
        L.append("  %s %s %s: %s; at bounds: %s; termination: %s" % (
            r["option"], r["model"], r["setting"], ", ".join(fmt(r.get("mult_" + m[8:])) for m in MULT),
            r.get("at_bounds", "-") or "none", r.get("termination", "-")))
    if dspec:
        L += ["", "convergence (%s):" % dspec["_about"]]
        for r in rows:
            if r["option"] == "ref":
                continue
            L.append("  %s %s %s: %s; %s%s" % (
                r["option"], r["model"], r["setting"], r.get("convergence"), r.get("conv_note", "-"),
                ("; scaled gradient %s, GN predicted relative decrease %s, best trial %s [%s]" % (
                    fmt(r.get("diag_g_scaled_inf")), fmt(r.get("diag_GN_pred_rel")), fmt(r.get("diag_best_trial_rel")),
                    r.get("diag_best_trial"))) if r.get("diag_files") else ""))
    if any("prev_data_cost" in r for r in rows):
        L += ["", "before the continuation (the stage-1 fit each continuation started from; data cost, penalty, nfev, "
                  "termination) -> after:"]
        for r in rows:
            if "prev_data_cost" in r:
                short = lambda m: (m or "-").replace("The maximum number of function evaluations is exceeded.",  # noqa: E731
                                                     "max_nfev").replace("`ftol` termination condition is satisfied.", "ftol")
                L.append("  %s %s %s: %s (penalty %s, %s evaluations, %s) -> %s (penalty %s, %s evaluations, %s); data cost x%.4f" % (
                    r["option"], r["model"], r["setting"], fmt(r["prev_data_cost"]), fmt(r["prev_penalty"]), fmt(r["prev_nfev"]),
                    short(r["prev_termination"]), fmt(r["data_cost"]), fmt(r.get("penalty")), fmt(r.get("nfev")),
                    short(r.get("termination")), r["data_cost"] / r["prev_data_cost"]))
    if any("penalty_mult" in r for r in rows):
        L += ["", "penalty of the option-C variant with the multiplier penalty = smoothness + P_m:"] + [
            "  %s %s %s: P_m %s, smoothness %s" % (r["option"], r["model"], r["setting"], fmt(r["penalty_mult"]),
                                                  fmt(r["penalty"] - r["penalty_mult"])) for r in rows if "penalty_mult" in r]
    L += ["", "best setting per option (lowest data cost): " + "; ".join(
        "%s %s %s (%.4e)" % (k, v["model"], v["setting"], v["data_cost"]) for k, v in sorted(best.items())), "",
          "sources:"] + ["  %s %s %s: %s" % (r["option"], r["model"], r["setting"], r["source"]) for r in rows]
    (out / "summary.txt").write_text("\n".join(L) + "\n")
    json.dump({k: {"model": v["model"], "setting": v["setting"], "data_cost": v["data_cost"]} for k, v in best.items()},
              open(out / "best.json", "x"), indent=1)
    print("\n".join(L[:len(rows) + 6]))
    return 0


if __name__ == "__main__":
    sys.exit(main())
