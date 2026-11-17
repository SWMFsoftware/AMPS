#!/usr/bin/env python3
# ==== version-controlled copy: header added 2026-11-13, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/figures_cluster/scripts/options_residual_figure.py (on Voyager: /data001/highresolu_vt/cops_output/figures_cluster/scripts/options_residual_figure.py)
# SHA256 of the run copy = SHA256 of this file without these header lines: 0877407eee37ec947ad8c57509e24b24d5508277a435cab3e2f3111e6d76e97c
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: ln(model / observation) vs altitude of several fits, NG and RG in two panels, rms recomputed and checked
#   against summary.csv (2026-11-14; legend in one column for long labels, 2026-11-17)
# RUN: called by figures_12m_nonneg_options_stage1.sh (RESID=1) and stage2_figures.py (CPU)
# ==== end of header ====
"""Log residuals ln(model / observation) of several fits against altitude, NG and RG in two panels (2026-11-14, lead: "one
residual figure: ln(model/obs) vs altitude for NG and RG (two panels, best of each option + c9 + D2), so the fit differences
are visible directly").

Run from <AMPS>/srcCG/LiouvilleCOPS_cuda with .venv/bin/python (CPU):
    options_residual_figure.py --bundle <bundle> --out <prefix> [--title t] [--summary <summary.csv>]
                               --fit "LABEL=TABLE[:style]" ... [--key "LABEL=option|setting" ...]
TABLE: positive_table.py output (observation rows of a fit's COPS_model.dat). The residual of an observation is
ln(pModel / pObserved) on the rows the fit used (use > 0, observed > 0): for the late-window rows pModel is the W10x5
average of the surface model + grain written into COPS_model.dat, so these are the fit's own log residuals (checked against
pressure_residuals.dat to 4.7e-7 for two fits, 2026-11-14). Every fit is drawn solid in the colour it has in the
model_obs_figures.py figure with the same --fit order (figures.PALETTE); the style in --fit is ignored (noisy residual lines
cannot carry dash patterns). The rms per fit and gauge is recomputed and, with --summary and --key, compared with the
fit's summary (rms_NG_all / rms_RG_all; must agree to 5e-4 absolute, else the script stops). Writes <prefix>.eps/.png,
<prefix>_caption.txt (with the recomputed rms) and <prefix>_check.json. Refuses existing files.
"""
import argparse
import csv
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path.cwd()))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from cops_cuda import figures  # noqa: E402
from model_obs_figures import parse_fit  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bundle", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--title", default=None)
    ap.add_argument("--fit", action="append", required=True)
    ap.add_argument("--summary", default=None)
    ap.add_argument("--key", action="append", default=[])
    ap.add_argument("--caption-intro", default="")
    a = ap.parse_args()
    out = Path(a.out)
    for f in (out.with_name(out.name + "_caption.txt"), out.with_name(out.name + "_check.json")):
        if f.exists():
            raise SystemExit("%s exists" % f)
    plt = figures.style()
    from matplotlib.lines import Line2D
    fits = [(lab, figures.read_table(p), p, figures.PALETTE[k]) for k, (lab, p, _ls) in
            enumerate(parse_fit(f) for f in a.fit)]
    ref = fits[0][1]
    for lab, tab, p, _ in fits[1:]:
        for k in ("iPoint", "refNG", "refRG", "useNG", "useRG", "altitude"):
            if not np.array_equal(tab[k], ref[k]):
                raise SystemExit("%s: column %s differs from the first fit's table" % (p, k))
    window = figures.late_window(ref, bundle=a.bundle)
    alt = ref["altitude"]
    order = np.argsort(-alt)
    keys = dict(k.split("=", 1) for k in a.key)
    summ = {}
    if a.summary:
        summ = {"%s|%s" % (r["option"], r["setting"]): r for r in csv.DictReader(open(a.summary))}
    check, resid = {}, {}
    for lab, tab, p, _ in fits:
        check[lab] = {"table": str(p)}
        for g in ("NG", "RG"):
            use = (tab["use" + g] > 0) & (tab["ref" + g] > 0) & (tab["p" + g] > 0)
            r = np.full(len(alt), np.nan)
            r[use] = np.log(tab["p" + g][use] / tab["ref" + g][use])
            resid[(lab, g)] = r
            rms = float(np.sqrt(np.mean(r[use] ** 2)))
            c = {"n": int(use.sum()), "rms": rms, "mean": float(np.mean(r[use])), "max_abs": float(np.max(np.abs(r[use])))}
            late = use & np.isin(tab["iPoint"], window["points"])
            if late.any():
                c["rms_late"] = float(np.sqrt(np.mean(r[late] ** 2)))
            k = keys.get(lab)
            if k and k in summ and summ[k].get("rms_%s_all" % g) not in (None, "", "-"):
                s = float(summ[k]["rms_%s_all" % g])
                c["rms_summary"] = s
                if abs(s - rms) > 5e-4:
                    raise SystemExit("%s %s: recomputed rms %.6f differs from the summary %.6f" % (lab, g, rms, s))
            check[lab][g] = c
    lim = max(check[lab][g]["max_abs"] for lab, *_ in fits for g in ("NG", "RG"))   # every residual inside the axes
    lim = float(np.ceil(lim * 10 + 1e-9) / 10)
    fig, ax = plt.subplots(2, 1, figsize=(figures.DOUBLE_COLUMN_IN, 5.6), sharex=True,
                           gridspec_kw={"hspace": 0.08})
    nleg = len(fits)
    longest = max(len("%s; rms NG 0.000, RG 0.000" % lab) for lab, *_ in fits)
    ncol = 2 if longest <= 66 else 1      # 2026-11-17: one column for long labels (the stage-1 sets: 2 columns, unchanged)
    fig.subplots_adjust(top=0.95 - 0.028 * (-(-nleg // ncol)), bottom=0.09, left=0.1, right=0.98)
    for i, g in enumerate(("NG", "RG")):
        x = ax[i]
        x.axvspan(window["altitude_min_m"], window["altitude_max_m"], facecolor=figures.LATE_FILL,
                  edgecolor=figures.LATE_EDGE, lw=0.5, zorder=0)
        x.axhline(0.0, color="k", lw=0.6, zorder=1)
        for k, (lab, tab, _p, color) in enumerate(fits):
            r = resid[(lab, g)][order]
            x.plot(alt[order], r, "-", color=color, lw=0.7, zorder=10 + len(fits) - k)
        x.set_xscale("log")
        x.set_xlim(alt.max() * 1.3, max(alt.min(), 1e-3) * 0.7)
        x.set_ylim(-lim, lim)
        x.set_ylabel("ln(model / %s data)" % g)
        x.text(0.01, 0.95, "(%s) %s" % ("ab"[i], g), transform=x.transAxes, va="top", fontsize=8)
    ax[1].set_xlabel("altitude above the shape model [m] (descent →)")
    ax[1].text(0.99, 0.05, "shaded: " + figures.window_label(window), transform=ax[1].transAxes, va="bottom", ha="right",
               fontsize=6.3, color="0.15", bbox={"facecolor": figures.LATE_FILL, "edgecolor": figures.LATE_EDGE, "lw": 0.5,
                                                 "pad": 1.5})
    labels = ["%s; rms NG %.3f, RG %.3f" % (lab, check[lab]["NG"]["rms"], check[lab]["RG"]["rms"]) for lab, *_ in fits]
    fig.legend([Line2D([], [], color=c, lw=1.2) for *_, c in fits], labels, loc="upper center",
               bbox_to_anchor=(0.5, 0.975), ncol=ncol, fontsize=6.5, handlelength=2.2, columnspacing=1.2)
    if a.title:
        fig.suptitle(a.title, y=0.998, fontsize=8)
    written = figures.save(fig, out)
    plt.close(fig)
    cap = ("%sLog residuals ln(model / observation) of %s against the gauge altitude (logarithmic, descent from left to right): "
           "(a) NG, (b) RG; one point per fitted observation (934 NG, 1150 RG; NG not fitted from 10:24:40 UTC), joined in "
           "descent order, every fit solid in its colour of the corresponding model / data figure; zero line black; the late "
           "window shaded (%s). For the late-window observations the model is the gauges' 10-s / 5-readout response average "
           "W10x5 of the 1-s grid model (NG with the fitted grain term), exactly as fitted, so these are the fits' own log "
           "residuals. rms ln over the fitted observations (recomputed from the tables; late window in brackets): %s. "
           "Positive: model above the data; the vertical range (±%.1f) holds every residual." % (
               a.caption_intro, "; ".join(lab for lab, *_ in fits), figures.window_label(window).replace("\n", ", "),
               "; ".join("%s: NG %.4f (%.4f), RG %.4f (%.4f)" % (lab, check[lab]["NG"]["rms"], check[lab]["NG"].get("rms_late", np.nan),
                                                              check[lab]["RG"]["rms"], check[lab]["RG"].get("rms_late", np.nan))
                         for lab, *_ in fits), lim))
    written.append(figures.write_caption(out, cap))
    json.dump({"fits": check, "ylim": lim, "files": written, "command": sys.argv},
              open(out.with_name(out.name + "_check.json"), "x"), indent=1)
    print("written: %s" % written)
    return 0


if __name__ == "__main__":
    sys.exit(main())
