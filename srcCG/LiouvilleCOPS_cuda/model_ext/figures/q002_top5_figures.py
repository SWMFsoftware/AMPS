#!/usr/bin/env python3
# ==== version-controlled copy: header added 2026-11-21, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/figures_cluster/scripts/q002_top5_figures.py (on Voyager: /data001/highresolu_vt/cops_output/figures_cluster/scripts/q002_top5_figures.py)
# SHA256 of the run copy = SHA256 of this file without these header lines: b083cd20a24cb3e4cd68ef6cf47da2cf24d8357aa2174ced6b03dcec4301392c
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: figures and numbers.json of the five most plausible new-law fits: model/data and residuals (unchanged
#   model_obs_figures.py / options_residual_figure.py), criteria panels, convergence (2026-11-21)
# RUN: CPU job jobs/q002_top5_figures.sbatch; output figures_cluster/12M/q002_top5_models/
# CHANGED: 2026-11-21 user decisions 1b / 2b (FINDINGS A21)
# ==== end of header ====
"""Figures and numbers of the five most plausible new-law fits (2026-11-21), for the Word document
12M/q002_top5_models/12M_q002_top5_models.docx.

The five fits (ranks 1-5 of 12M/q002_plausibility/README.md): C f10 / P_m 0.01 (templates on), C f1 / P_m 0.01
(templates on), C f10 / P_m 0.01 (templates off), A (sum of squares, rank 3), B (background, bg 1e-2). All with the new
H2O law q_min = 0.02 q_max.

1. Model vs COPS (time, altitude, NG / RG) and the log residuals of the five fits, with the unchanged figure tools of
   the comparison set (model_obs_figures.py, options_residual_figure.py; run copies in this folder) and the model
   tables of 12M/model_vs_observations_nonneg_options_stage2_q002/tables/ (the same final states).
2. 12M_q002_top5_criteria: the plausibility criteria per fit (data cost, late rms ln NG / RG, active area, Q above the
   energy limit, seen share of Q, Q max / min over the map times).
3. 12M_q002_top5_convergence: data cost of the fit and of the fine tuning per chain cycle.
Every figure: .png and .eps at 300 dpi and a _caption.txt. numbers.json: every number used, with its source file.

    q002_top5_figures.py --out <new folder>        (CPU; never overwrites)
"""
import argparse
import csv
import hashlib
import json
import subprocess
import sys
from pathlib import Path

import numpy as np

APP = Path("/data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda")
sys.path.insert(0, str(APP))
import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from cops_cuda import figures  # noqa: E402

ROOT = Path("/data001/highresolu_vt/cops_output")
F12 = ROOT / "figures_cluster/12M"
MV = F12 / "model_vs_observations_nonneg_options_stage2_q002"
PL = F12 / "q002_plausibility"
S = ROOT / "figures_cluster/scripts"
PY = str(APP / ".venv/bin/python")
BUNDLE = ROOT / "runs/shap8_012m_lategrid_dz1.5/bundle/restart_templates.npz"
TEN_LABEL = {"D2_q002": "D2", "softplus_q002": "softplus", "A_q002": "A", "B_q002": "B", "D_q002": "D",
             "C_f1_pm1em4_q002": "C f1 / 1e-4", "C_f10_pm1em4_q002": "C f10 / 1e-4",
             "C_f1_pm0p01_q002": "C f1 / 0.01", "C_f10_pm0p01_q002": "C f10 / 0.01",
             "C_f10_pm0p01_tploff_q002": "C f10 / 0.01, off"}
# key in summary.csv, key in metrics.json / source_rate_q002, full-set key, short label, line style
RUNS = [("C_f10_pm0p01_q002", "C_f10_pm0p01_q002_reg", "C_f10_pm0p01", "C f10, P_m 0.01, templates on", "solid"),
        ("C_f1_pm0p01_q002", "C_f1_pm0p01_q002_reg", "C_f1_pm0p01", "C f1, P_m 0.01, templates on", "dashed"),
        ("C_f10_pm0p01_tploff_q002", "C_f10_pm0p01_tploff_q002_reg", "C_f10_pm0p01_tploff",
         "C f10, P_m 0.01, templates off", "dotted"),
        ("A_q002", "A_q002_reg", "A", "A sum of squares, rank 3", "dashdot"),
        ("B_q002", "B_q002_reg", "B", "B background, bg 1e-2", "solid")]
LAW = "new H2O law q_min = 0.02 q_max"
INK, MUTED, GRID = "#0b0b0b", "#52514e", "#e4e3df"


def sha256(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def run(cmd, log):
    log.write("$ " + " ".join(cmd) + "\n")
    log.flush()
    r = subprocess.run(cmd, cwd=str(APP), stdout=log, stderr=subprocess.STDOUT)
    if r.returncode:
        raise SystemExit("failed (%d): %s" % (r.returncode, " ".join(cmd)))


def style(ax):
    ax.grid(True, axis="x", color=GRID, lw=0.5)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)


def save(fig, out, stem, caption):
    for ext in ("png", "eps"):
        p = out / ("%s.%s" % (stem, ext))
        if p.exists():
            raise SystemExit("%s exists" % p)
        fig.savefig(p, dpi=300)
    (out / ("%s_caption.txt" % stem)).write_text(caption + "\n")
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    O = Path(a.out)
    FIG = O / "figures"
    if FIG.exists():
        raise SystemExit("%s exists (never overwritten)" % FIG)
    FIG.mkdir(parents=True)
    log = open(O / "figures.log", "x")

    rows = {r["key"]: r for r in csv.DictReader(open(MV / "summary.csv"))}
    met = json.load(open(PL / "metrics.json"))["runs"]
    aa_p = F12 / "q002_decisions_20261121/active_area.json"
    aa = json.load(open(aa_p))                          # switched-off area (user decision 1b, 2026-11-21)
    off = aa["runs"]
    num = {"sources": {"summary.csv": [str(MV / "summary.csv"), sha256(MV / "summary.csv")],
                       "metrics.json": [str(PL / "metrics.json"), sha256(PL / "metrics.json")]},
           "law": LAW, "runs": {}}
    num["sources"]["active_area.json"] = [str(aa_p), sha256(aa_p)]
    num["ten_fits"] = []
    for mk, o in off.items():
        key = "D2_q002" if mk.startswith("D2") else mk[:-len("_reg")]
        num["ten_fits"].append(dict(o, key=key, label=TEN_LABEL[key], data_cost=float(rows[key]["data_cost"])))
    for k, mk, fk, lab, _ in RUNS:
        r, m = rows[k], met[mk]
        facts_p = F12 / ("full_set_q002_%s/checks/facts.json" % fk)
        facts = json.load(open(facts_p))
        q = [v[0] for v in facts["Q"].values()]
        num["sources"]["facts_" + fk] = [str(facts_p), sha256(facts_p)]
        num["runs"][k] = {
            "label": lab, "summary_label": r["label"], "option": r["option"], "model": r["model"],
            "setting": r["setting"], "root": r["root"], "note": r["note"],
            "data_cost": float(r["data_cost"]), "P_delta": float(r["ft_delta_penalty"]),
            "fit_penalty": float(r["fit_penalty"]) if r["fit_penalty"] else None,
            "fit_termination": r["fit_termination"], "fit_nfev": int(r["fit_nfev"]),
            "ft_termination": r["ft_termination"], "ft_nfev": int(r["ft_nfev"]),
            "rms_ln_NG_all": float(r["rms_NG_all"]), "rms_ln_RG_all": float(r["rms_RG_all"]),
            "rms_ln_NG_late": float(r["rms_NG_late"]), "rms_ln_RG_late": float(r["rms_RG_late"]),
            "multipliers": {b: float(r["mult_" + b]) for b in ("40_50", "50_60", "60_70", "70_80", "80_90", "90_180")},
            "at_bounds": r["at_bounds"], "kappa": float(r["kappa"]),
            "Q_H2O_1033": float(r["Q_H2O_1033"]), "Q_H2O_seen_1033": float(r["Q_H2O_seen_1033"]),
            "seen_share": float(r["Q_H2O_seen_1033"]) / float(r["Q_H2O_1033"]),
            "old_data_cost": float(r["old_data_cost"]), "data_cost_change": float(r["data_cost_change"]),
            "history_fit_data": [float(x) for x in r["history_fit_data"].split()],
            "history_ft_data": [float(x) for x in r["history_ft_data"].split()],
            "off_seen": off[mk]["off_T"], "off_seen_floor_counted": off[mk]["off_T_floor_counted"],
            "zero_area": m["whole"]["zero_area"], "active_10pct_trim_whole": m["whole"]["ge_10pct_trim"],
            "active_10pct_trim_seen": m["seen"]["ge_10pct_trim"], "active_1pct_trim_whole": m["whole"]["ge_1pct_trim"],
            "active_1pct_trim_seen": m["seen"]["ge_1pct_trim"], "gini_whole": m["whole"]["gini"],
            "Q_share_over_F_E": m["energy"]["Q_share_over_F_E"], "F_max_over_F_E": m["energy"]["F_max_over_F_E"],
            "F_E": m["energy"]["F_E_m-2_s-1"], "rh_AU": m["energy"]["rh_AU"],
            "Q_day_max": max(q), "Q_day_min": min(q), "Q_day_ratio": max(q) / min(q),
            "seen_share_facts": facts["production_seen_unseen"]["10:33:10"]["seen_share"]}

    # ---- 1. model vs data and residuals (unchanged tools, comparison-set tables)
    fits = []
    for k, _, _, lab, st in RUNS:
        fits += ["--fit", "%s, %.4e=%s:%s" % (lab, num["runs"][k]["data_cost"], MV / "tables" / (k + ".dat"), st)]
    run([PY, "-W", "ignore", str(S / "model_obs_figures.py"), "--bundle", str(BUNDLE), "--out",
         str(FIG / "12M_q002_top5"), "--ref", str(MV / "tables" / "softplus_q002.dat"),
         "--title", "SHAP8 12M, %s: the five selected fits vs the data" % LAW] + fits, log)
    res = []
    for k, _, _, lab, _ in RUNS:
        r = num["runs"][k]
        L = "%s, %.4e" % (lab, r["data_cost"])
        res += ["--fit", "%s=%s" % (L, MV / "tables" / (k + ".dat")), "--key", "%s=%s|%s" % (L, r["option"], r["setting"])]
    run([PY, "-W", "ignore", str(S / "options_residual_figure.py"), "--bundle", str(BUNDLE), "--out",
         str(FIG / "12M_q002_top5_residuals"), "--summary", str(MV / "summary.csv"),
         "--title", "SHAP8 12M, %s: log residuals of the five selected fits" % LAW,
         "--caption-intro", "The five selected fits of the new-law round (%s; ranking 12M/q002_plausibility). " % LAW]
        + res, log)

    # ---- 2. criteria
    labs = [lab for _, _, _, lab, _ in RUNS]
    col = [figures.PALETTE[i] for i in range(len(RUNS))]
    R = [num["runs"][k] for k, *_ in RUNS]
    y = np.arange(len(RUNS))[::-1]
    panels = [("(a) final data cost (10$^{-3}$)", [r["data_cost"] * 1e3 for r in R], None, "lower is better"),
              ("(b) late-window rms ln(model/data)", [r["rms_ln_NG_late"] for r in R],
               [r["rms_ln_RG_late"] for r in R], "filled NG, open RG"),
              ("(c) seen surface switched off (%)", [100 * r["off_seen"] for r in R],
               [100 * r["off_seen_floor_counted"] for r in R], "open: imposed floor counted"),
              ("(d) share of Q$_{H_2O}$ above the energy limit (%)", [100 * r["Q_share_over_F_E"] for r in R], None,
               "10:33:10 UTC"),
              ("(e) share of Q$_{H_2O}$ from seen terrain (%)", [100 * r["seen_share"] for r in R], None,
               "10:33:10 UTC"),
              ("(f) Q$_{H_2O}$ max / min over the 13 map times", [r["Q_day_ratio"] for r in R], None,
               "00:00-10:39:26 UTC")]
    fig, axes = plt.subplots(2, 3, figsize=(13.5, 6.4), dpi=300, sharey=True)
    for ax, (title, v1, v2, note) in zip(axes.ravel(), panels):
        for i in range(len(RUNS)):
            ax.plot([0, v1[i]], [y[i], y[i]], color=col[i], lw=2.0, solid_capstyle="round")
            ax.plot(v1[i], y[i], "o", color=col[i], ms=6, zorder=3)
            ax.annotate("%.3g" % v1[i], (v1[i], y[i]), xytext=(6, 4), textcoords="offset points", fontsize=7.2,
                        color=INK)
            if v2 is not None:
                ax.plot(v2[i], y[i] - 0.22, "o", mfc="white", mec=col[i], mew=1.4, ms=6, zorder=3)
                ax.annotate("%.3g" % v2[i], (v2[i], y[i] - 0.22), xytext=(6, -9), textcoords="offset points",
                            fontsize=7.2, color=MUTED)
        hi = max(v1 + (v2 or []))
        if title.startswith("(c)"):
            ax.axvline(50, color=MUTED, lw=0.8, ls="--", zorder=1)
            ax.text(50, len(RUNS) - 0.55, " fail above 50 %", fontsize=7, color=MUTED, va="top")
        ax.set_xlim(0, hi * 1.28 if hi > 0 else 1)
        ax.set_ylim(-1.05, len(RUNS) - 0.5)     # room below the last row for its open marker and the panel note
        ax.set_title(title, loc="left", fontsize=9, color=INK)
        ax.text(0.99, 0.02, note, transform=ax.transAxes, ha="right", va="bottom", fontsize=7, color=MUTED)
        style(ax)
    for ax in axes[:, 0]:
        ax.set_yticks(y)
        ax.set_yticklabels(labs, fontsize=8)
    fig.tight_layout()
    cap = ("The five selected fits of the new-law round (%s), in the order of the withdrawn first ranking (top to bottom). "
           "(a) final data cost (0.5 r.r of the data rows); (b) rms of ln(model/data) in the late window 10:33:10-10:39:26 "
           "UTC (297-5.5 m), NG filled, RG open; (c) share of the seen area that is switched off (user decision 1b): activity without the imposed floor "
           "below T = 2.11e15, 1 %% of the srcCG median (filled; open: with the floor counted; only B has one), "
           "fail above 50 %%; (d) share of Q_H2O at "
           "10:33:10 UTC from facets whose flux exceeds the energy limit F_E = S0 / (r_h^2 dH) = %.3e m^-2 s^-1 at r_h = "
           "%.4f AU (all sunlight at normal incidence into sublimation; an upper bound); (e) share of Q_H2O from terrain "
           "the gauges saw; (f) highest / lowest Q_H2O over the 13 map times 00:00-10:39:26 UTC (activity fixed, "
           "illumination changing). Sources: model_vs_observations_nonneg_options_stage2_q002/summary.csv, "
           "q002_plausibility/metrics.json, full_set_q002_<run>/checks/facts.json." % (LAW, R[0]["F_E"], R[0]["rh_AU"]))
    save(fig, FIG, "12M_q002_top5_criteria", cap)

    # ---- 3. convergence
    fig, ax = plt.subplots(figsize=(7.6, 4.2), dpi=300)
    for i, r in enumerate(R):
        c = np.arange(1, len(r["history_fit_data"]) + 1)
        ax.plot(c, np.array(r["history_fit_data"]) * 1e3, color=col[i], lw=1.4, ls=figures_style(RUNS[i][4]),
                marker="o", ms=3.5, label=labs[i])
        ax.plot(c, np.array(r["history_ft_data"]) * 1e3, color=col[i], lw=0, marker="s", ms=4, mfc="white", mew=1.1)
    ax.set_xlabel("chain cycle")
    ax.set_ylabel("data cost (10$^{-3}$)")
    ax.set_xticks(range(1, 6))
    ax.grid(True, color=GRID, lw=0.5)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.legend(fontsize=7.2, frameon=False, loc="upper right", title="circles: fit, squares: fine tuning",
              title_fontsize=7)
    fig.tight_layout()
    cap = ("Data cost per chain cycle of the five selected fits (%s): circles after the fit of the cycle, squares "
           "after its regularized late-facet fine tuning (the state the next cycle starts from). Source: "
           "model_vs_observations_nonneg_options_stage2_q002/summary.csv (history_fit_data, history_ft_data)." % LAW)
    save(fig, FIG, "12M_q002_top5_convergence", cap)

    (O / "numbers.json").write_text(json.dumps(num, indent=1))
    log.write("done\n")
    print("wrote", FIG, O / "numbers.json")


def figures_style(name):
    return {"solid": "-", "dashed": (0, (4, 2)), "dotted": (0, (1, 1.5)), "dashdot": (0, (5, 1.5, 1, 1.5))}[name]


if __name__ == "__main__":
    main()
