#!/usr/bin/env python3
# ==== version-controlled copy: header added 2026-11-13, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/figures_cluster/scripts/stage2_figures.py (on Voyager: /data001/highresolu_vt/cops_output/figures_cluster/scripts/stage2_figures.py)
# SHA256 of the run copy = SHA256 of this file without these header lines: 3aee62eb8b5990ef584252e7bd7750418ff3ea63b5c80df8409e5bc172a26121
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: final stage-2 comparison of the non-negative activity chains (regularized and unregularized late-facet fine
#   tuning): model/data figures, residuals, cost table with terminations and bounds, captions; small products, disk
#   check (2026-11-17)
# RUN: CPU: srun ... stage2_figures.py --spec figures/data/stage2_final_spec_20261117.json --out <new folder>
# ==== end of header ====
"""Final stage-2 comparison of the non-negative activity models (2026-11-17, lead: "prepare the FINAL stage-2 comparison";
same figure style as model_vs_observations_nonneg_options_stage1_final: NG/RG vs data with the late window shaded and
labelled, residuals, a cost table with final cost / termination / active bounds; .eps + .png 300 dpi; small products only).

    stage2_figures.py --spec <spec.json> --out <new folder> [--min-free-gb 25]       (CPU; run from anywhere)

spec.json (figures_cluster/scripts/data/stage2_final_spec_*.json): the chains (key, option, model, setting, label, reg,
kind nonneg|localnn, root relative to runs/shap8_012m_lategrid_dz1.5/nonneg, cycles, note) and D2. A chain is used when its
last cycle's fine tuning is complete (c<N>/ft/finetune_result.txt "complete yes", c<N>/ft/final/COPS_model.dat) and, for
option C, c<N>/analysis/report.json exists; otherwise it is listed as pending. Per used chain, read only:
  * final state = c<N>/ft/final: COPS_model.dat (observation rows -> tables/<key>.dat via positive_table.py), and
    finetune_result.txt: final data cost, delta penalty (regularized), termination, nfev, rms ln NG/RG (all, late window),
    facets at the cap, active bounds of the fine tuning;
  * the last fit c<N>/fit/final: fit_result.txt (data cost, penalty = final cost - data cost, or smoothness + P_m for C,
    optimizer termination and nfev); multipliers and parameters at a bound from fit_restart.npz (A, B, D, softplus: x at
    lo / hi within 1e-9 max(1, |bound|)) or from c<N>/analysis/report.json (C: c_j at 0, multipliers at 1, kappa / grain at 0);
  * the data cost after each cycle's fit and fine tuning (history).
Figures (model_obs_figures.py; the data selection of the regularized softplus reference):
  12M_stage2_reg_best_*      regularized: softplus reference, D2, A, B, best C (lowest final data cost), D
  12M_stage2_C_reg_*         all regularized C chains + softplus reference
  12M_stage2_C_unreg_*       the unregularized C chains + best regularized C + softplus reference
  12M_stage2_<A|B|D>_*       the option's regularized and unregularized chain + softplus reference
  12M_stage2_softplus_*      softplus regularized vs unregularized (5 and 9 cycles)
  12M_stage2_reg_best_residuals, 12M_stage2_unreg_best_residuals (options_residual_figure.py)
plus summary.txt / summary.csv, captions, CAPTIONS.md, code/. Refuses an existing --out; stops when the file system holding
--out has less than --min-free-gb free (before and after every figure call).
"""
import argparse
import csv
import json
import re
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = Path("/data001/highresolu_vt/cops_output")
RUN = ROOT / "runs/shap8_012m_lategrid_dz1.5"
NN = RUN / "nonneg"
APP = Path("/data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda")
PY = str(APP / ".venv/bin/python")
BUNDLE = RUN / "bundle/restart_templates.npz"
MULT = ("40_50", "50_60", "60_70", "70_80", "80_90", "90_180")
sys.path.insert(0, str(HERE))
from options_captions import SUFFIX  # noqa: E402

COLS = ["key", "option", "model", "setting", "label", "regularized", "cycles", "data_cost", "ft_delta_penalty",
        "ft_termination", "ft_nfev", "ft_facets_at_cap", "ft_active_bounds", "fit_data_cost", "fit_penalty",
        "fit_termination", "fit_nfev", "fit_segments", "rms_NG_all", "rms_RG_all", "rms_NG_late", "rms_RG_late",
        "mult_40_50", "mult_50_60", "mult_60_70", "mult_70_80", "mult_80_90", "mult_90_180", "at_bounds",
        "history_fit_data", "history_ft_data", "note", "root"]


def free_gb(p):
    return shutil.disk_usage(p).free / 1e9


def check_disk(p, need):
    f = free_gb(p)
    if f < need:
        raise SystemExit("only %.1f GB free on the file system of %s (< %.0f GB): stopped" % (f, p, need))
    return f


def kv(path):
    """finetune_result.txt: '<name> initial <a> final <b>' rows and '<name> <value>' rows."""
    out = {}
    for line in Path(path).read_text().splitlines():
        if line.startswith("#"):
            continue
        m = re.match(r"^(\S+)\s+initial\s+(\S+)\s+final\s+(\S+)", line)
        if m:
            out[m.group(1)] = m.group(3)
            continue
        p = line.split(None, 1)
        if len(p) == 2:
            out.setdefault(p[0], p[1].strip())
    return out


def fit_result(path):
    t = Path(path).read_text()
    g = lambda k: (float(re.search(r"^%s (\S+)" % k, t, re.M).group(1)) if re.search(r"^%s (\S+)" % k, t, re.M) else None)  # noqa: E731
    m = re.search(r"^optimizer: (.*?) \((.*)\)$", t, re.M)
    nf = re.search(r"nfev (\d+)", m.group(2)) if m else None
    return {"final_cost": g("final_cost"), "final_data_cost": g("final_data_cost"), "final_penalty": g("final_penalty"),
            "final_mult_penalty": g("final_mult_penalty"), "termination": m.group(1) if m else None,
            "nfev": int(nf.group(1)) if nf else None}


def chain_row(c):
    root = NN / c["root"]
    n = c["cycles"]
    ft = root / ("c%d" % n) / "ft"
    res = ft / "final" / "finetune_result.txt"
    if not res.exists() or "complete yes" not in res.read_text() or not (ft / "final" / "COPS_model.dat").exists():
        return None, "c%d fine tuning not complete" % n
    if c["kind"] == "localnn" and not (root / ("c%d" % n) / "analysis" / "report.json").exists():
        return None, "c%d analysis/report.json missing" % n
    k = kv(res)
    fr = fit_result(root / ("c%d" % n) / "fit" / "final" / "fit_result.txt")
    row = {x: c.get(x) for x in ("key", "option", "model", "setting", "label", "cycles", "note")}
    row.update({"regularized": "yes" if c["reg"] else "no", "root": str(root),
                "data_cost": float(k["cost_data"]), "ft_delta_penalty": float(k["cost_delta_penalty"]) if "cost_delta_penalty" in k else "",
                "ft_termination": k.get("termination"), "ft_nfev": k.get("nfev"), "ft_facets_at_cap": k.get("facets_at_cap", ""),
                "ft_active_bounds": k.get("active_bounds", ""),
                "rms_NG_all": float(k["rms_ln_NG"]), "rms_RG_all": float(k["rms_ln_RG"]),
                "rms_NG_late": float(k["late_window_rms_ln_NG"]), "rms_RG_late": float(k["late_window_rms_ln_RG"]),
                "fit_data_cost": fr["final_data_cost"], "fit_termination": fr["termination"], "fit_nfev": fr["nfev"],
                "fit_segments": len([s for s in (root / ("c%d" % n) / "fit").glob("seg[0-9][0-9]") if s.is_dir()])})
    if c["kind"] == "localnn":
        rp = json.load(open(root / ("c%d" % n) / "analysis" / "report.json"))
        row["fit_penalty"] = rp["costs"]["penalty"] + rp["costs"].get("mult_penalty", 0.0)
        b = rp["bounds"]
        ab = ["c_j at 0: %d nodes (%.3f of all)" % (b["c_at_0_nodes"], b["c_at_0_fraction_all_nodes"])]
        ab += [x.replace("_at_0", " at 0") for x in ("kappa_at_0", "grain_A_at_0", "grain_B_at_0") if b.get(x)]
        if b.get("multipliers_at_1"):
            ab.append("multipliers at 1: %s" % ", ".join(m.replace("slope_", "") for m in b["multipliers_at_1"]))
        row["at_bounds"] = "; ".join(ab)
        for i, m in enumerate(MULT):
            row["mult_" + m] = rp["multipliers"][i + 1]
    else:
        row["fit_penalty"] = (fr["final_cost"] - fr["final_data_cost"]) if fr["final_cost"] is not None else ""
        with np.load(root / ("c%d" % n) / "fit" / "final" / "fit_restart.npz", allow_pickle=False) as z:
            names = [str(x) for x in json.loads(str(z["meta_json"]))["names"]]
            x, lo, hi = (np.asarray(z[q], np.float64) for q in ("x", "lo", "hi"))
        ab = []
        for i, nm in enumerate(names):
            if np.isfinite(lo[i]) and x[i] - lo[i] <= 1e-9 * max(1.0, abs(lo[i])):
                ab.append(nm + " (lower)")
            elif np.isfinite(hi[i]) and hi[i] - x[i] <= 1e-9 * max(1.0, abs(hi[i])):
                ab.append(nm + " (upper)")
        row["at_bounds"] = "; ".join(ab) or "none"
        for m in MULT:
            row["mult_" + m] = float(x[names.index("c_slope_" + m)]) if ("c_slope_" + m) in names else ""
    hf, ht = [], []
    for j in range(1, n + 1):
        f = root / ("c%d" % j) / "fit" / "final" / "fit_result.txt"
        hf.append("%.4e" % fit_result(f)["final_data_cost"] if f.exists() else "-")
        r = root / ("c%d" % j) / "ft" / "final" / "finetune_result.txt"
        ht.append("%.4e" % float(kv(r)["cost_data"]) if r.exists() else "-")
    row["history_fit_data"], row["history_ft_data"] = " ".join(hf), " ".join(ht)
    return row, None


def d2_row(d2):
    fr = fit_result(d2 / "fit_result.txt")
    t = (d2 / "fit_result.txt").read_text()
    rn = re.search(r"^rms_ln_residual_NG (\S+)", t, re.M)
    rr = re.search(r"^rms_ln_residual_RG (\S+)", t, re.M)
    return {"key": "D2", "option": "ref", "model": "clipped", "setting": "D2", "label": "clipped D2 (start of stage 1)",
            "regularized": "-", "cycles": "", "data_cost": fr["final_data_cost"], "fit_data_cost": fr["final_data_cost"],
            "fit_termination": fr["termination"], "fit_nfev": fr["nfev"], "fit_penalty": 0.0,
            "rms_NG_all": float(rn.group(1)), "rms_RG_all": float(rr.group(1)), "root": str(d2),
            "note": "clipped refit D2 (lategrid base, W10x5), no fine tuning"}


def fmt(v, f="%.4g"):
    if v in (None, ""):
        return "-"
    return f % v if isinstance(v, float) else str(v)


def short(m):
    return (m or "-").replace("The maximum number of function evaluations is exceeded.", "max_nfev").replace(
        "`ftol` termination condition is satisfied.", "ftol").replace("`xtol` termination condition is satisfied.", "xtol")


def run(cmd, log):
    log.write("$ " + " ".join(cmd) + "\n")
    log.flush()
    r = subprocess.run(cmd, cwd=str(APP), stdout=log, stderr=subprocess.STDOUT)
    if r.returncode:
        raise SystemExit("failed (%d): %s" % (r.returncode, " ".join(cmd)))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--spec", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--min-free-gb", type=float, default=25.0)
    a = ap.parse_args()
    O = Path(a.out)
    if O.exists():
        raise SystemExit("%s exists" % O)
    gb0 = check_disk(O.parent, a.min_free_gb)
    spec = json.load(open(a.spec))
    O.mkdir()
    (O / "tables").mkdir()
    (O / "code").mkdir()
    log = open(O / "build.log", "x")
    log.write("free at start: %.1f GB\n" % gb0)
    rows, pending = [], []
    for c in spec["chains"]:
        r, why = chain_row(c)
        if r is None:
            pending.append({"key": c["key"], "label": c["label"], "root": str(NN / c["root"]), "why": why, "note": c.get("note")})
            continue
        run([PY, str(HERE / "positive_table.py"), str(Path(r["root"]) / ("c%d" % c["cycles"]) / "ft" / "final"),
             str(O / "tables" / (c["key"] + ".dat"))], log)
        rows.append(r)
    d2 = RUN / spec["d2"]
    rows.append(d2_row(d2))
    run([PY, str(HERE / "positive_table.py"), str(d2), str(O / "tables" / "D2.dat")], log)
    by = {r["key"]: r for r in rows}
    if "softplus_reg" not in by:
        raise SystemExit("the regularized softplus reference is not complete")
    reg = [r for r in rows if r["regularized"] == "yes"]
    best = {}
    for o in ("A", "B", "C", "D"):
        for rg in ("yes", "no"):
            cand = [r for r in rows if r["option"] == o and r["regularized"] == rg]
            if cand:
                best[(o, rg)] = min(cand, key=lambda r: r["data_cost"])
    with open(O / "summary.csv", "x", newline="") as f:
        w = csv.DictWriter(f, fieldnames=COLS, extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, "") for k in COLS})
    # summary.txt
    L = ["# Final stage-2 comparison (%s): final state of each chain = its last cycle's fine tuning (c<N>/ft/final); data "
         "cost = 0.5 r.r of the data rows at that state; delta penalty = the regularized fine tuning's P_delta (f_delta C_ref "
         "mean_w (h + delta)^2 over the tuned facets); fit = the last cycle's fit (penalty: D spread, C smoothness + P_m); rms "
         "ln of the final state (all used rows; late window W10x5); bounds of the last fit." % a.spec, ""]
    hdr = "%-26s %-4s %-10s %-9s %-15s %-10s %-9s %-15s %-15s %-15s" % (
        "chain", "reg", "data", "P_delta", "ft stop (nfev)", "fit data", "fit pen.", "fit stop (nfev)", "rms NG all/late",
        "rms RG all/late")
    L.append(hdr)
    for r in sorted(rows, key=lambda r: (r["regularized"] != "yes", r["data_cost"])):
        L.append("%-26s %-4s %-10s %-9s %-15s %-10s %-9s %-15s %-15s %-15s" % (
            r["key"], r["regularized"], fmt(r["data_cost"], "%.4e"), fmt(r.get("ft_delta_penalty"), "%.2e"),
            "%s (%s)" % (short(r.get("ft_termination")), r.get("ft_nfev", "-")) if r.get("ft_termination") else "-",
            fmt(r.get("fit_data_cost"), "%.4e"), fmt(r.get("fit_penalty"), "%.2e"),
            "%s (%s%s)" % (short(r.get("fit_termination")), r.get("fit_nfev"),
                           ", seg %s" % r["fit_segments"] if (r.get("fit_segments") or 1) > 1 else ""),
            "%s/%s" % (fmt(r.get("rms_NG_all"), "%.4f"), fmt(r.get("rms_NG_late"), "%.4f")),
            "%s/%s" % (fmt(r.get("rms_RG_all"), "%.4f"), fmt(r.get("rms_RG_late"), "%.4f"))))
    L += ["", "multipliers of the last fit (40-50, 50-60, 60-70, 70-80, 80-90, 90-180 deg) and parameters at a bound:"]
    for r in rows:
        if r["key"] == "D2":
            continue
        L.append("  %s: %s; at bounds: %s; fine tuning: facets at the cap %s, active bounds %s" % (
            r["key"], ", ".join(fmt(r.get("mult_" + m)) for m in MULT), r.get("at_bounds"), r.get("ft_facets_at_cap") or "-",
            r.get("ft_active_bounds") or "-"))
    L += ["", "data cost after each cycle (fit | fine tuning):"]
    for r in rows:
        if r["key"] != "D2":
            L.append("  %-26s fit %s | ft %s" % (r["key"], r["history_fit_data"], r["history_ft_data"]))
    L += ["", "best per option (lowest final data cost): " + "; ".join(
        "%s %s %s (%.4e)" % (o, "regularized" if rg == "yes" else "unregularized", b["setting"], b["data_cost"])
        for (o, rg), b in sorted(best.items())), "", "notes:"] + ["  %s: %s" % (r["key"], r.get("note")) for r in rows]
    L += ["", "pending (not included): " + ("; ".join("%s (%s; %s)" % (p["key"], p["why"], p["note"]) for p in pending)
                                             if pending else "none"), "", "Reaching max_nfev is not convergence.",
          "", "roots:"] + ["  %s: %s" % (r["key"], r["root"]) for r in rows]
    (O / "summary.txt").write_text("\n".join(L) + "\n")
    json.dump({"rows": [r["key"] for r in rows], "pending": pending,
               "best": {"%s|%s" % k: v["key"] for k, v in best.items()}}, open(O / "included.json", "x"), indent=1)

    # figure groups
    ref = ("softplus_reg", "solid")
    def ent(key, style):
        r = by[key]
        return r, style
    groups = []
    g = [ent("softplus_reg", "solid"), ent("D2", "dotted")]
    g += [ent(best[(o, "yes")]["key"], st) for o, st in (("A", "dashed"), ("B", "dashdot"), ("C", "dashed"), ("D", "dashdot"))
          if (o, "yes") in best]
    groups.append(("12M_stage2_reg_best", "SHAP8 12M stage 2: regularized chains (best of each option) vs the data", g,
                   "the final states of the regularized stage-2 chains of options A, B, D and the best regularized option-C "
                   "chain (lowest final data cost), with the regularized softplus reference and the clipped refit D2"))
    cr = sorted([r for r in reg if r["option"] == "C"], key=lambda r: r["data_cost"])
    styles = ["dashed", "dashdot", "dotted", "dashed", "dashdot", "dotted"]
    groups.append(("12M_stage2_C_reg", "SHAP8 12M stage 2: regularized option-C chains vs the data",
                   [ent("softplus_reg", "solid")] + [(r, styles[i]) for i, r in enumerate(cr)],
                   "the final states of all completed regularized option-C chains (local non-negative basis), with the "
                   "regularized softplus reference"))
    cu = sorted([r for r in rows if r["option"] == "C" and r["regularized"] == "no"], key=lambda r: r["data_cost"])
    if cu:
        groups.append(("12M_stage2_C_unreg", "SHAP8 12M stage 2: unregularized option-C chains vs the data",
                       [ent("softplus_reg", "solid")] + ([(best[("C", "yes")], "dotted")] if ("C", "yes") in best else []) +
                       [(r, styles[i]) for i, r in enumerate(cu)],
                       "the final states of the unregularized option-C chains (unregularized late-facet fine tuning), with "
                       "the best regularized option-C chain and the regularized softplus reference"))
    for o, nm in (("A", "sum of squares rank 3"), ("B", "background bg 1e-2"), ("D", "log spread f 0.1")):
        gg = [ent("softplus_reg", "solid")] + [(r, st) for r, st in
                                              ((by.get("%s_reg" % o), "dashed"), (by.get("%s_unreg" % o), "dotted")) if r]
        groups.append(("12M_stage2_%s" % o, "SHAP8 12M stage 2: option %s (%s), regularized vs unregularized" % (o, nm), gg,
                       "option %s (%s): the final states of its regularized and unregularized stage-2 chains, with the "
                       "regularized softplus reference" % (o, nm)))
    groups.append(("12M_stage2_softplus", "SHAP8 12M stage 2: softplus, regularized vs unregularized",
                   [(by[k], st) for k, st in (("softplus_reg", "solid"), ("softplus_unreg_c5", "dashed"),
                                               ("softplus_unreg_c9", "dotted")) if k in by],
                   "the softplus model (slope multipliers >= 1): the regularized stage-2 chain (5 cycles) and the "
                   "unregularized chain chain_5cycles_smm1 after 5 and 9 cycles (c9 = the stage-1 reference)"))
    plan = []
    for name, title, fits, _ in groups:
        cmd = [PY, "-W", "ignore", str(HERE / "model_obs_figures.py"), "--bundle", str(BUNDLE), "--out", str(O / name),
               "--ref", str(O / "tables" / "softplus_reg.dat"), "--title", title]
        for r, st in fits:
            cmd += ["--fit", "%s, %.4e=%s:%s" % (r["label"], r["data_cost"], O / "tables" / (r["key"] + ".dat"), st)]
        plan.append(cmd)
    for nm, rg, title in (("12M_stage2_reg_best_residuals", "yes", "SHAP8 12M stage 2: log residuals, regularized chains"),
                          ("12M_stage2_unreg_best_residuals", "no", "SHAP8 12M stage 2: log residuals, unregularized chains")):
        fits = [by["softplus_reg"]] + ([by["D2"]] if rg == "yes" else ([by["softplus_unreg_c5"]] if "softplus_unreg_c5" in by else []))
        fits += [best[(o, rg)] for o in ("A", "B", "C", "D") if (o, rg) in best]
        cmd = [PY, "-W", "ignore", str(HERE / "options_residual_figure.py"), "--bundle", str(BUNDLE), "--out", str(O / nm),
               "--summary", str(O / "summary.csv"), "--title", title,
               "--caption-intro", "Final stage 2 (%s). " % ("regularized late-facet fine tuning" if rg == "yes"
                                                             else "unregularized late-facet fine tuning; the regularized "
                                                                  "softplus reference for comparison")]
        for r in fits:
            lab = "%s, %.4e" % (r["label"], r["data_cost"])
            cmd += ["--fit", "%s=%s" % (lab, O / "tables" / (r["key"] + ".dat")), "--key", "%s=%s|%s" % (lab, r["option"], r["setting"])]
        plan.append(cmd)
    with open(O / "figure_plan.sh", "x") as f:
        f.write("#!/bin/bash\n# the figure commands of this folder (run from %s)\n" % APP)
        for cmd in plan:
            f.write(" ".join("'%s'" % x.replace("'", "'\\''") for x in cmd) + " || exit 1\n")
    for cmd in plan:
        check_disk(O, a.min_free_gb)
        run(cmd, log)
    # captions of the model/data figures
    COMMON = ("Stage 2: each chain = 5 cycles of a fit (max_nfev 4000; slope-bin multipliers 40-50 ... 90-180 deg >= 1, D2's "
              "fine-tuning templates ft_c1..ft_c6 held at 1 unless stated) and a late-facet fine tuning (the 271 facets seen "
              "only in the last 5 minutes; late_only, max_nfev 500); final state = the last fine tuning. Regularized chains "
              "(stage2r_*): fine-tuning penalty P_delta = f_delta C_ref mean_w (h + delta)^2 with f_delta = 10, C_ref = "
              "4.2688e-3 and the cap |h + delta| <= ln 100; unregularized chains (stage2_*, and the softplus chain "
              "chain_5cycles_smm1): no penalty on delta. Option C: local non-negative basis on a 150-m grid with smoothness "
              "and multiplier penalty P_m = f_m C_ref sum_b (ln m_b)^2 (f_m in the label), CPU fits; option D with the "
              "fixed-level spread penalty (u0 = -1.8641).")
    data = ("Black (NG) and grey (RG) dots: the COPS observations used in the fits (934 NG, 1150 RG; NG excluded from "
            "10:24:40 UTC; the 95 late-window observations compared with the gauges' 10-s / 5-readout response averages "
            "W10x5 of the 1-s grid model, the NG model including the fitted grain term). Lines: the final states "
            "(observation rows of COPS_model.dat).")

    def desc(r):
        if r["key"] == "D2":
            return "clipped D2: data cost %.4e, rms ln NG %.4f, RG %.4f, fit %s" % (r["data_cost"], r["rms_NG_all"],
                                                                                   r["rms_RG_all"], short(r["fit_termination"]))
        return ("%s: final data cost %.4e%s, rms ln NG %.4f (late %.4f), RG %.4f (late %.4f); last fit data cost %.4e%s, "
                "%s after %s evaluations%s; last fine tuning %s after %s evaluations%s" % (
                    r["label"], r["data_cost"], (", P_delta %.2e" % r["ft_delta_penalty"]) if r.get("ft_delta_penalty") != "" else "",
                    r["rms_NG_all"], r["rms_NG_late"], r["rms_RG_all"], r["rms_RG_late"], r["fit_data_cost"],
                    (", fit penalty %.2e" % r["fit_penalty"]) if r.get("fit_penalty") not in ("", 0.0) else "",
                    short(r["fit_termination"]), r["fit_nfev"],
                    (" (in the last of %d segments; the fit continued after a time limit)" % r["fit_segments"])
                    if (r.get("fit_segments") or 1) > 1 else "", short(r["ft_termination"]), r["ft_nfev"],
                    ("; " + r["note"]) if r.get("note") else ""))
    for name, title, fits, what in groups:
        for sfx, w in SUFFIX.items():
            if not (O / (name + sfx + ".png")).exists():
                continue
            cap = ("Modelled COPS pressures of %s, %s. %s %s Fits (first = solid): %s. Reaching max_nfev is not convergence. "
                   "Table: summary.txt." % (what, w, COMMON, data, "; ".join(desc(r) for r, _ in fits)))
            (O / (name + sfx + "_caption.txt")).open("x").write(cap + "\n")
    if pending:
        (O / "PENDING.txt").write_text("Chains not included yet (to be added when complete):\n" + "\n".join(
            "  %s: %s (%s; %s)" % (p["key"], p["root"], p["why"], p["note"]) for p in pending) + "\n")
    run([PY, str(HERE / "collect_captions.py"), str(O), "--title",
         "Figure captions: final stage-2 comparison of the non-negative activity models", "--order", ""], log)
    for f in ("stage2_figures.py", "positive_table.py", "model_obs_figures.py", "options_residual_figure.py",
              "options_captions.py", "collect_captions.py"):
        shutil.copy2(HERE / f, O / "code" / f)
    shutil.copy2(a.spec, O / "code" / Path(a.spec).name)
    git = subprocess.run(["bash", "-c", "git -C /data001/highresolu_vt/AMPS rev-parse HEAD; git -C /data001/highresolu_vt/AMPS "
                                        "status --short"], capture_output=True, text=True).stdout
    (O / "code" / "git.txt").write_text(git)
    subprocess.run(["bash", "-c", "cd '%s' && sha256sum $(ls | grep -v SHA256SUMS) > SHA256SUMS" % (O / "code")], check=True)
    gb1 = check_disk(O, 0)
    log.write("free at end: %.1f GB\n" % gb1)
    print("\n".join(L[:len(rows) + 4]))
    print("pending: %s" % [p["key"] for p in pending])
    print("free: %.1f GB -> %.1f GB; done: %s" % (gb0, gb1, O))
    return 0


if __name__ == "__main__":
    sys.exit(main())
