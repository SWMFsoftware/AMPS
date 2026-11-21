#!/usr/bin/env python3
# ==== version-controlled copy: header added 2026-11-19, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/figures_cluster/scripts/stage2_q002_figures.py (on Voyager: /data001/highresolu_vt/cops_output/figures_cluster/scripts/stage2_q002_figures.py)
# SHA256 of the run copy = SHA256 of this file without these header lines: 5897bf012875f7198a478e6c95b6add9a568dc0691c20eebd8d24bf6c9b3b66a
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: comparison set of the new-law (H2O q_min = 0.02 q_max) stage-2 round vs the data and vs the old law: best of
#   each option, the C chains, new vs old per run, residuals, convergence per chain, summary with kappa, Q_H2O and the
#   cost change (2026-11-19)
# RUN: CPU job jobs/stage2_q002_figures.sbatch (OUT=<new folder>); imports stage2_figures.py and kappa_table.py
# ==== end of header ====
"""Comparison set of the NEW-LAW stage-2 round (H2O q_min = 0.02 q_max, the night-side fix of 2026-11-18) against the
COPS data and against the old-law round (q_min = 0.2 q_max). Lead 2026-11-19 17:05 CST, user: "prepare a full set of
comparison plots for the last round of fitting"; the same set as model_vs_observations_nonneg_options_stage2_final/
(stage2_figures.py, whose readers are imported unchanged), for the new-law states.

    stage2_q002_figures.py --spec data/stage2_q002_spec_20261119.json --out <new folder> [--min-free-gb 25]   (CPU)

Figures (.eps + .png 300 dpi; model/data figures by model_obs_figures.py: _time, _altitude, _altitude_NG, _altitude_RG,
_altitude_NG_RG; late window shaded and labelled with UTC and altitude range):
  12M_q002_best_*            new law: softplus reference, D2 refit, A, B, best C (lowest final data cost), D
  12M_q002_C_*               new law: the five option-C chains + softplus reference
  12M_q002_vs_old_<run>_*    per run: new law (solid) vs old law (dashed): softplus, A, B, D, D2, each C setting
  12M_q002_best_residuals, 12M_q002_C_residuals, 12M_q002_vs_old_<run>_residuals   ln(model/obs) vs altitude
                             (options_residual_figure.py; rms checked against summary.csv)
  12M_q002_convergence_<run> data cost after every cycle's fit and fine tuning, new law vs old law; and
  12M_q002_convergence_all   all new-law chains in one panel
summary.txt / summary.csv: per run (new and old law) final data cost, P_delta, terminations / nfev of the last fit and
fine tuning, rms ln NG / RG (all, late), multipliers, parameters at a bound, kappa (and at its bound), Q_H2O total and of
the seen facets at 10:33:10 UTC (new law: source_rate_q002/<run>/inspection/checks/inspect_103310.json, else the C
report.json production_1033.localnn; old law: the same inspection files' old-law block), the old-law data cost and the
change. Captions state the law in every figure. Refuses an existing --out; stops below --min-free-gb free.
"""
import argparse
import csv
import json
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import stage2_figures as S2  # noqa: E402
import kappa_table as KT  # noqa: E402
from options_captions import SUFFIX  # noqa: E402

ROOT, RUN, NN, APP, PY, BUNDLE = S2.ROOT, S2.RUN, S2.NN, S2.APP, S2.PY, S2.BUNDLE
SRQ = ROOT / "figures_cluster/12M/source_rate_q002"
LAWN, LAWO = "new law q_min 0.02 q_max", "old law q_min 0.2 q_max"
EXTRA = ["law", "kappa", "kappa_at_bound", "Q_H2O_1033", "Q_H2O_seen_1033", "Q_source", "old_key", "old_data_cost",
         "data_cost_change", "old_Q_H2O_1033", "old_Q_H2O_seen_1033"]
COLS = S2.COLS[:5] + ["law"] + [c for c in S2.COLS[5:] if c not in ("note", "root")] + [c for c in EXTRA if c != "law"] + \
    ["note", "root"]


def kappa_of(row, kind):
    root = Path(row["root"])
    if kind == "localnn":
        k = KT.from_report(root / "c5" / "analysis" / "report.json")
    elif kind == "d2":
        k = KT.from_fit(root)
        if "kappa" not in k:
            with np.load(root / "fit_restart.npz", allow_pickle=False) as z:
                k = {"kappa": float(np.asarray(z["kappa"]).reshape(-1)[0]), "at_lower_bound": None,
                     "source": "fit_restart.npz kappa"}
    else:
        k = KT.from_fit(root / "c5" / "fit" / "final")
    return k.get("kappa"), k.get("at_lower_bound")


def q_numbers(srkey, root, kind):
    """(Q new total, Q new seen, source, Q old total, Q old seen) at 10:33:10."""
    f = SRQ / srkey / "inspection" / "checks" / "inspect_103310.json"
    if f.exists():
        r = json.load(open(f))
        return (r["new"]["Q_H2O_s-1"], r["new"]["Q_H2O_seen_s-1"], str(f), r["old"]["Q_H2O_s-1"],
                r["old"]["Q_H2O_seen_s-1"])
    if kind == "localnn":
        # total: the final state's maps (c5/maps_ft manifest, 10:33:10); seen: the C report (production_1033.localnn,
        # computed for the c5 FIT state, i.e. before the last fine tuning: its totals agree with the final maps to
        # 5e-8 ... 1.2e-4 relative in this round)
        rp = root / "c5" / "analysis" / "report.json"
        p = json.load(open(rp))["production_1033"]["localnn"]
        mf = root / "c5" / "maps_ft" / "maps_manifest.json"
        if mf.exists():
            t = [t for t in json.load(open(mf))["times"] if t["utc"].endswith("10:33:10")][0]
            return (t["Q_H2O_s-1"], p["Q_H2O_seen_s-1"], "total: %s; seen: %s production_1033.localnn (c5 fit state)" % (
                mf, rp), "", "")
        return p["Q_H2O_all_s-1"], p["Q_H2O_seen_s-1"], "%s production_1033.localnn (c5 fit state)" % rp, "", ""
    return "", "", "not available", "", ""


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--spec", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--min-free-gb", type=float, default=25.0)
    a = ap.parse_args()
    O = Path(a.out)
    if O.exists():
        raise SystemExit("%s exists" % O)
    gb0 = S2.check_disk(O.parent, a.min_free_gb)
    spec = json.load(open(a.spec))
    old_spec = json.load(open(HERE / "data" / spec["old_spec"]))
    olds = {c["key"]: c for c in old_spec["chains"]}
    O.mkdir()
    (O / "tables").mkdir()
    (O / "code").mkdir()
    log = open(O / "build.log", "x")
    log.write("free at start: %.1f GB\n" % gb0)
    rows, missing, checks = [], [], {}

    def table(row, state_dir):
        S2.run([PY, str(HERE / "positive_table.py"), str(state_dir), str(O / "tables" / (row["key"] + ".dat"))], log)

    for c in spec["chains"]:
        r, why = S2.chain_row(c)
        if r is None:
            missing.append({"key": c["key"], "why": why})
            continue
        r["law"] = "new"
        r["label"] = c["label"].replace("new law", LAWN)
        r["kappa"], r["kappa_at_bound"] = kappa_of(r, c["kind"])
        qn, qs, src, qo, qos = q_numbers(c["srkey"], Path(r["root"]), c["kind"])
        r.update({"Q_H2O_1033": qn, "Q_H2O_seen_1033": qs, "Q_source": src, "old_Q_H2O_1033": qo,
                  "old_Q_H2O_seen_1033": qos})
        if c["kind"] == "localnn" and "inspect_103310" in src:     # cross-check with the C report
            p = json.load(open(Path(r["root"]) / "c5" / "analysis" / "report.json"))["production_1033"]["localnn"]
            checks["Q_inspect_vs_report_" + c["key"]] = {
                "total_rel": abs(qn / p["Q_H2O_all_s-1"] - 1), "seen_rel": abs(qs / p["Q_H2O_seen_s-1"] - 1)}
        table(r, Path(r["root"]) / "c5" / "ft" / "final")
        rows.append(r)
        oc = dict(olds[c["old"]])
        ro, why = S2.chain_row(oc)
        if ro is None:
            missing.append({"key": "old " + c["old"], "why": why})
            continue
        ro.update({"key": "old_" + c["old"], "law": "old", "setting": oc["setting"] + " old law",
                   "label": oc["label"].replace(" (reference)", "") + ", " + LAWO})
        ro["kappa"], ro["kappa_at_bound"] = kappa_of(ro, oc["kind"])
        table(ro, Path(ro["root"]) / ("c%d" % oc["cycles"]) / "ft" / "final")
        rows.append(ro)
        r["old_key"], r["old_data_cost"] = ro["key"], ro["data_cost"]
        r["data_cost_change"] = r["data_cost"] / ro["data_cost"] - 1.0
    d = spec["d2"]
    for law, path, key, lab in (("new", d["path"], d["key"], d["label"].replace("new law", LAWN)),
                                ("old", d["old_path"], "old_D2", "clipped D2, " + LAWO)):
        r = S2.d2_row(RUN / path)
        r.update({"key": key, "law": law, "label": lab, "setting": "D2 %s law" % law,
                  "note": ("the refitted D2 (D2_q002_warm: warm start from D2, one segment of 4000 evaluations)" if law == "new"
                           else "clipped refit D2 (lategrid base, W10x5), no fine tuning")})
        r["kappa"], r["kappa_at_bound"] = kappa_of(r, "d2")
        table(r, RUN / path)
        rows.append(r)
    by = {r["key"]: r for r in rows}
    dq = by[d["key"]]
    qn, qs, src, qo, qos = q_numbers(d["srkey"], RUN / d["path"], "d2")
    dq.update({"Q_H2O_1033": qn, "Q_H2O_seen_1033": qs, "Q_source": src, "old_Q_H2O_1033": qo, "old_Q_H2O_seen_1033": qos,
               "old_key": "old_D2", "old_data_cost": by["old_D2"]["data_cost"],
               "data_cost_change": dq["data_cost"] / by["old_D2"]["data_cost"] - 1.0})
    for r in rows:                                   # old-law Q from the new run's inspection (its old-law block)
        if r["law"] == "old":
            nk = [x for x in rows if x.get("old_key") == r["key"]]
            if nk:
                r["Q_H2O_1033"], r["Q_H2O_seen_1033"] = nk[0].get("old_Q_H2O_1033", ""), nk[0].get("old_Q_H2O_seen_1033", "")
                r["Q_source"] = nk[0]["Q_source"] + " (old-law block)" if "inspect" in nk[0]["Q_source"] else "not available"
    new = [r for r in rows if r["law"] == "new"]
    best = {}
    for o in ("A", "B", "C", "D"):
        cand = [r for r in new if r["option"] == o]
        if cand:
            best[o] = min(cand, key=lambda r: r["data_cost"])
    with open(O / "summary.csv", "x", newline="") as f:
        w = csv.DictWriter(f, fieldnames=COLS, extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, "") for k in COLS})
    # summary.txt
    fmt, short = S2.fmt, S2.short
    L = ["# Comparison of the NEW-LAW stage-2 round (H2O %s; the night-side fix of 2026-11-18) with the old-law round (%s); "
         "spec %s. Final state of each chain = c5/ft/final (D2: final); data cost = 0.5 r.r of the data rows at that state; "
         "P_delta = the regularized fine tuning's penalty; fit = the last cycle's fit; rms ln of the final state (all used "
         "rows; late window W10x5); bounds of the last fit; kappa (CO2 scale) of the last fit; Q_H2O at 10:33:10 UTC." % (
             LAWN, LAWO, a.spec), ""]
    L.append("%-26s %-4s %-10s %-10s %-9s %-9s %-15s %-10s %-9s %-20s %-15s %-15s %-9s %-9s %-9s" % (
        "run", "law", "data", "old data", "change", "P_delta", "ft stop (nfev)", "fit data", "fit pen.", "fit stop (nfev)",
        "rms NG all/late", "rms RG all/late", "kappa", "Q_H2O", "Q seen"))
    for r in sorted(rows, key=lambda r: (r["law"] != "new", r["data_cost"])):
        L.append("%-26s %-4s %-10s %-10s %-9s %-9s %-15s %-10s %-9s %-20s %-15s %-15s %-9s %-9s %-9s" % (
            r["key"], r["law"], fmt(r["data_cost"], "%.4e"), fmt(r.get("old_data_cost"), "%.4e"),
            fmt(100 * r["data_cost_change"], "%+.1f%%") if r.get("data_cost_change", "") != "" else "-",
            fmt(r.get("ft_delta_penalty"), "%.2e"),
            "%s (%s)" % (short(r.get("ft_termination")), r.get("ft_nfev", "-")) if r.get("ft_termination") else "-",
            fmt(r.get("fit_data_cost"), "%.4e"), fmt(r.get("fit_penalty"), "%.2e"),
            "%s (%s%s)" % (short(r.get("fit_termination")), r.get("fit_nfev"),
                           ", seg %s" % r["fit_segments"] if (r.get("fit_segments") or 1) > 1 else ""),
            "%s/%s" % (fmt(r.get("rms_NG_all"), "%.4f"), fmt(r.get("rms_NG_late"), "%.4f")),
            "%s/%s" % (fmt(r.get("rms_RG_all"), "%.4f"), fmt(r.get("rms_RG_late"), "%.4f")),
            fmt(r.get("kappa"), "%.2e"), fmt(r.get("Q_H2O_1033"), "%.3e"), fmt(r.get("Q_H2O_seen_1033"), "%.3e")))
    L += ["", "multipliers of the last fit (40-50, 50-60, 60-70, 70-80, 80-90, 90-180 deg), parameters at a bound, kappa:"]
    for r in rows:
        if r["key"] in ("D2_q002", "old_D2"):
            L.append("  %s: no slope multipliers; kappa %s (at its lower bound: %s)" % (r["key"], fmt(r.get("kappa"), "%.3e"),
                                                                                       r.get("kappa_at_bound")))
            continue
        L.append("  %s: %s; at bounds: %s; kappa %s (at its lower bound: %s); fine tuning: facets at the cap %s, active "
                 "bounds %s" % (r["key"], ", ".join(fmt(r.get("mult_" + m)) for m in S2.MULT), r.get("at_bounds"),
                                fmt(r.get("kappa"), "%.3e"), r.get("kappa_at_bound"), r.get("ft_facets_at_cap") or "-",
                                r.get("ft_active_bounds") or "-"))
    L += ["", "data cost after each cycle (fit | fine tuning):"]
    for r in rows:
        if r.get("history_fit_data"):
            L.append("  %-28s fit %s | ft %s" % (r["key"], r["history_fit_data"], r["history_ft_data"]))
    L += ["", "Q_H2O at 10:33:10 UTC [s^-1] (total, seen facets) and its source:"]
    for r in rows:
        L.append("  %-28s %s, %s  (%s)" % (r["key"], fmt(r.get("Q_H2O_1033"), "%.4e"), fmt(r.get("Q_H2O_seen_1033"), "%.4e"),
                                          r.get("Q_source", "-")))
    L += ["", "best per option, new law (lowest final data cost): " + "; ".join(
        "%s %s (%.4e)" % (o, b["setting"], b["data_cost"]) for o, b in sorted(best.items())),
        "", "checks: %s" % json.dumps(checks), "", "not available / not included: " + (
            "; ".join("%s (%s)" % (m["key"], m["why"]) for m in missing) if missing else "none"),
        "", "notes:"] + ["  %s: %s" % (r["key"], r.get("note")) for r in rows]
    L += ["", "Reaching max_nfev is not convergence. The data costs of both laws are against the same observations "
          "(934 NG, 1150 RG, the same W10x5 late window); the models differ in the H2O night / shadow source.",
          "", "roots:"] + ["  %s: %s" % (r["key"], r["root"]) for r in rows]
    (O / "summary.txt").write_text("\n".join(L) + "\n")
    json.dump({"rows": [r["key"] for r in rows], "missing": missing, "best": {o: b["key"] for o, b in best.items()},
               "checks": checks}, open(O / "included.json", "x"), indent=1)

    # ---- model / data figure groups
    styles = ["dashed", "dashdot", "dotted", "dashed", "dashdot", "dotted"]
    groups = []
    g = [(by["softplus_q002"], "solid"), (dq, "dotted")]
    g += [(best[o], st) for o, st in (("A", "dashed"), ("B", "dashdot"), ("C", "dashed"), ("D", "dashdot")) if o in best]
    groups.append(("12M_q002_best", "SHAP8 12M stage 2, new H2O law: best of each option vs the data", g,
                   "the NEW-LAW final states of the regularized stage-2 chains of options A, B, D and the best option-C chain "
                   "(lowest final data cost), with the new-law softplus reference and the refitted D2"))
    cr = sorted([r for r in new if r["option"] == "C"], key=lambda r: r["data_cost"])
    groups.append(("12M_q002_C", "SHAP8 12M stage 2, new H2O law: the option-C chains vs the data",
                   [(by["softplus_q002"], "solid")] + [(r, styles[i]) for i, r in enumerate(cr)],
                   "the NEW-LAW final states of all five option-C chains (local non-negative basis), with the new-law softplus "
                   "reference"))
    pairs = [r for r in new if r.get("old_key") in by]
    for r in pairs:
        o = by[r["old_key"]]
        nm = "12M_q002_vs_old_%s" % r["key"].replace("_q002", "")
        groups.append((nm, "SHAP8 12M: %s, new vs old H2O law" % r["label"].split(",")[0], [(r, "solid"), (o, "dashed")],
                       "%s: the final state of the NEW-LAW run (solid) and of its OLD-LAW counterpart (dashed)" % (
                           r["label"].split(", new law")[0])))
    plan = []
    for name, title, fits, _ in groups:
        cmd = [PY, "-W", "ignore", str(HERE / "model_obs_figures.py"), "--bundle", str(BUNDLE), "--out", str(O / name),
               "--ref", str(O / "tables" / "softplus_q002.dat"), "--title", title]
        for r, st in fits:
            cmd += ["--fit", "%s, %.4e=%s:%s" % (r["label"], r["data_cost"], O / "tables" / (r["key"] + ".dat"), st)]
        plan.append(cmd)
    res_groups = [("12M_q002_best_residuals", "SHAP8 12M stage 2, new H2O law: log residuals, best of each option",
                   [r for r, _ in groups[0][2]]),
                  ("12M_q002_C_residuals", "SHAP8 12M stage 2, new H2O law: log residuals, option C", [by["softplus_q002"]] + cr)]
    for r in pairs:
        if r["option"] != "C" or r is best.get("C"):
            res_groups.append(("12M_q002_vs_old_%s_residuals" % r["key"].replace("_q002", ""),
                               "SHAP8 12M: %s, log residuals, new vs old H2O law" % r["label"].split(",")[0],
                               [r, by[r["old_key"]]]))
    for nm, title, fits in res_groups:
        cmd = [PY, "-W", "ignore", str(HERE / "options_residual_figure.py"), "--bundle", str(BUNDLE), "--out", str(O / nm),
               "--summary", str(O / "summary.csv"), "--title", title,
               "--caption-intro", "Stage 2 of the NEW H2O law (%s; the night-side fix of 2026-11-18; the old law had %s); "
                                  "the law of each fit is in its label. " % (LAWN.replace("new law ", ""),
                                                                              LAWO.replace("old law ", ""))]
        for r in fits:
            lab = "%s, %.4e" % (r["label"], r["data_cost"])
            cmd += ["--fit", "%s=%s" % (lab, O / "tables" / (r["key"] + ".dat")), "--key", "%s=%s|%s" % (lab, r["option"], r["setting"])]
        plan.append(cmd)
    with open(O / "figure_plan.sh", "x") as f:
        f.write("#!/bin/bash\n# the figure commands of this folder (run from %s)\n" % APP)
        for cmd in plan:
            f.write(" ".join("'%s'" % x.replace("'", "'\\''") for x in cmd) + " || exit 1\n")
    for cmd in plan:
        S2.check_disk(O, a.min_free_gb)
        S2.run(cmd, log)
    # ---- convergence figures (data cost per cycle)
    conv = convergence_figures(O, rows, by, new)
    # ---- captions
    COMMON = ("NEW H2O law: q_H2O = q_min on night (cos SZA <= 0) and cast-shadow facets, q_min + (q_max - q_min) cos SZA on "
              "day facets, with q_min = 0.02 q_max (srcCG's value; the night-side fix of 2026-11-18, FINDINGS A18); the OLD "
              "law used q_min = 0.2 q_max (not intended); the law of every line is in its label. New-law round (fitting base "
              "base_q002): each chain = 5 cycles of a fit (max_nfev 4000; slope-bin multipliers 40-50 ... 90-180 deg >= 1, "
              "D2's fine-tuning templates ft_c1..ft_c6 held at 1 unless stated) and a regularized late-facet fine tuning "
              "(the 271 facets seen only in the last 5 minutes; late_only, max_nfev 500; P_delta = f_delta C_ref mean_w (h + "
              "delta)^2, f_delta 10, C_ref 4.2688e-3 kept from the old-law round, cap |h + delta| <= ln 100); cycle-1 fits "
              "from the refitted D2 (D2_q002_warm: D2 refitted with the new law, warm start, one segment of 4000 "
              "evaluations); final state = the last fine tuning. Option C: local non-negative basis with smoothness and "
              "multiplier penalty P_m = f_m C_ref sum_b (ln m_b)^2 (f_m in the label). Old-law round: the regularized "
              "stage-2 chains of 2026-11-16/17 with the same scheme on the old base, and the old D2.")
    data = ("Black (NG) and grey (RG) dots: the COPS observations used in the fits (934 NG, 1150 RG; NG excluded from "
            "10:24:40 UTC; the 95 late-window observations compared with the gauges' 10-s / 5-readout response averages "
            "W10x5 of the 1-s grid model, the NG model including the fitted grain term). Lines: the final states "
            "(observation rows of COPS_model.dat).")

    def desc(r):
        if r["key"] in ("D2_q002", "old_D2"):
            return "%s: data cost %.4e, rms ln NG %.4f, RG %.4f, fit %s after %s evaluations" % (
                r["label"], r["data_cost"], r["rms_NG_all"], r["rms_RG_all"], short(r["fit_termination"]), r["fit_nfev"])
        return ("%s: final data cost %.4e%s, rms ln NG %.4f (late %.4f), RG %.4f (late %.4f); last fit data cost %.4e%s, "
                "%s after %s evaluations%s; last fine tuning %s after %s evaluations" % (
                    r["label"], r["data_cost"], (", P_delta %.2e" % r["ft_delta_penalty"]) if r.get("ft_delta_penalty") != "" else "",
                    r["rms_NG_all"], r["rms_NG_late"], r["rms_RG_all"], r["rms_RG_late"], r["fit_data_cost"],
                    (", fit penalty %.2e" % r["fit_penalty"]) if r.get("fit_penalty") not in ("", 0.0) else "",
                    short(r["fit_termination"]), r["fit_nfev"],
                    (" (in the last of %d segments)" % r["fit_segments"]) if (r.get("fit_segments") or 1) > 1 else "",
                    short(r["ft_termination"]), r["ft_nfev"]))
    for name, title, fits, what in groups:
        for sfx, w in SUFFIX.items():
            if not (O / (name + sfx + ".png")).exists():
                continue
            cap = ("Modelled COPS pressures of %s, %s. %s %s Fits (first = solid): %s. Reaching max_nfev is not convergence. "
                   "Table: summary.txt." % (what, w, COMMON, data, "; ".join(desc(r) for r, _ in fits)))
            (O / (name + sfx + "_caption.txt")).open("x").write(cap + "\n")
    S2.run([PY, str(HERE / "collect_captions.py"), str(O), "--title",
            "Figure captions: the new-law (q_min = 0.02 q_max) stage-2 round vs the data and vs the old law", "--order", ""],
           log)
    for f in ("stage2_q002_figures.py", "stage2_figures.py", "kappa_table.py", "positive_table.py", "model_obs_figures.py",
              "options_residual_figure.py", "options_captions.py", "collect_captions.py"):
        shutil.copy2(HERE / f, O / "code" / f)
    shutil.copy2(a.spec, O / "code" / Path(a.spec).name)
    shutil.copy2(HERE / "data" / spec["old_spec"], O / "code" / spec["old_spec"])
    git = subprocess.run(["bash", "-c", "git -C /data001/highresolu_vt/AMPS rev-parse HEAD; git -C /data001/highresolu_vt/AMPS "
                                        "status --short"], capture_output=True, text=True).stdout
    (O / "code" / "git.txt").write_text(git)
    subprocess.run(["bash", "-c", "cd '%s' && sha256sum $(ls | grep -v SHA256SUMS) > SHA256SUMS" % (O / "code")], check=True)
    gb1 = S2.check_disk(O, 0)
    log.write("free at end: %.1f GB\n" % gb1)
    print("\n".join(L[:len(rows) + 4]))
    print("missing: %s; convergence figures: %d" % (missing, len(conv)))
    print("free: %.1f GB -> %.1f GB; done: %s" % (gb0, gb1, O))
    return 0


def convergence_figures(O, rows, by, new):
    """Data cost after each cycle's fit (F) and fine tuning (T), new law solid, old law dashed; .eps/.png + captions."""
    sys.path.insert(0, str(APP))
    from cops_cuda import figures as FIG
    plt = FIG.style()
    out = []

    def hist(r):
        f = [float(x) if x != "-" else np.nan for x in r["history_fit_data"].split()]
        t = [float(x) if x != "-" else np.nan for x in r["history_ft_data"].split()]
        xs, ys = [], []
        for k, (a_, b_) in enumerate(zip(f, t)):
            xs += [k + 1.0, k + 1.5]
            ys += [a_, b_]
        return np.array(xs), np.array(ys)

    chains = [r for r in new if r.get("history_fit_data")]
    for r in chains:
        o = by.get(r.get("old_key"))
        fig, ax = plt.subplots(figsize=(4.6, 3.3))
        fig.subplots_adjust(left=0.16, right=0.97, top=0.91, bottom=0.36)
        for rr, ls, lab in ((r, "-", r["label"]), (o, "--", o["label"] if o else None)):
            if rr is None:
                continue
            x, y = hist(rr)
            ax.plot(x, y, ls=ls, marker="o", ms=3, color="k" if ls == "-" else "0.45", label=lab)
            ax.plot(x[1::2], y[1::2], ls="none", marker="s", ms=3.5, mfc="white", mec="k" if ls == "-" else "0.45")
        ax.set_yscale("log")
        ax.set_xticks([1, 2, 3, 4, 5])
        ax.set_xlabel("cycle (circle: after the fit; square: after the fine tuning)", fontsize=7)
        ax.set_ylabel("data cost")
        ax.legend(fontsize=6.3, loc="upper left", bbox_to_anchor=(-0.16, -0.2), frameon=False)
        ax.set_title(r["label"].split(",")[0], fontsize=7.5)
        name = "12M_q002_convergence_%s" % r["key"].replace("_q002", "")
        FIG.save(fig, O / name)
        plt.close(fig)
        cap = ("Convergence of the stage-2 chain %s: data cost (0.5 r.r of the data rows) after the fit (circles) and after "
               "the regularized late-facet fine tuning (squares, at cycle + 0.5) of each of the 5 cycles; solid black: the "
               "NEW H2O law (q_min = 0.02 q_max, the night-side fix of 2026-11-18; history %s | %s); dashed grey: the "
               "OLD-law chain (q_min = 0.2 q_max; %s; history %s | %s). Both against the same observations. Last fit: new law %s "
               "after %s evaluations. Reaching max_nfev is not convergence. Table: summary.txt." % (
                   r["key"], r["history_fit_data"], r["history_ft_data"], o["key"] if o else "-",
                   o["history_fit_data"] if o else "-", o["history_ft_data"] if o else "-", S2.short(r["fit_termination"]),
                   r["fit_nfev"]))
        (O / (name + "_caption.txt")).open("x").write(cap + "\n")
        out.append(name)
    fig, ax = plt.subplots(figsize=(FIG.DOUBLE_COLUMN_IN, 4.6))
    fig.subplots_adjust(left=0.09, right=0.98, top=0.93, bottom=0.36)
    cols = plt.get_cmap("tab10")
    for i, r in enumerate(sorted(chains, key=lambda r: r["data_cost"])):
        x, y = hist(r)
        ax.plot(x, y, "-", marker="o", ms=2.5, color=cols(i % 10), label="%s, %.4e" % (r["label"].replace(", " + LAWN, ""),
                                                                                           r["data_cost"]))
    d2 = by.get("D2_q002")
    if d2:
        ax.axhline(d2["data_cost"], color="0.3", lw=0.8, ls=":", label="refitted D2 (no cycles), %.4e" % d2["data_cost"])
    ax.set_yscale("log")
    ax.set_xticks([1, 2, 3, 4, 5])
    ax.set_xlabel("cycle (fit at the integer, fine tuning at + 0.5)")
    ax.set_ylabel("data cost")
    ax.set_title("SHAP8 12M stage 2, new H2O law (q$_\\mathrm{min}$ / q$_\\mathrm{max}$ 0.02): data cost per cycle", fontsize=8)
    ax.legend(fontsize=6.3, ncol=2, loc="upper left", bbox_to_anchor=(-0.02, -0.16), frameon=False)
    FIG.save(fig, O / "12M_q002_convergence_all")
    plt.close(fig)
    (O / "12M_q002_convergence_all_caption.txt").open("x").write(
        "Data cost after the fit (integer cycle) and after the regularized late-facet fine tuning (cycle + 0.5) of every "
        "cycle of the nine NEW-LAW stage-2 chains (H2O q_min = 0.02 q_max, the night-side fix of 2026-11-18; the old law "
        "had 0.2), legend = run and final data cost; dotted: the refitted D2 (new law, one fit, no cycles). Reaching "
        "max_nfev is not convergence. Table: summary.txt.\n")
    out.append("12M_q002_convergence_all")
    return out


if __name__ == "__main__":
    sys.exit(main())
