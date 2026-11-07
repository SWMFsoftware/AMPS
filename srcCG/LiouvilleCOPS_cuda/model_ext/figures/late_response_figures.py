#!/usr/bin/env python3
# ==== version-controlled copy: header added 2026-11-07, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/figures_cluster/scripts/late_response_figures.py (on Voyager: /data001/highresolu_vt/cops_output/figures_cluster/scripts/late_response_figures.py)
# SHA256 of the run copy = SHA256 of this file without these header lines: 5cd5a65050a3afb0799540344bd78f7177afd0bb7b5fbdf826753fd287d1ee2c
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: figure script. Model vs COPS data in the late window with the gauge time response (instantaneous 1-s
#   series, W10, W10x5) per gauge offset for a pair of models, and all models on the zero-offset base; .eps + .png
#   300 dpi. User 2026-11-07.
# USED BY: figures/figures_12m_late_response.sh.
# RUN (from <AMPS>/srcCG/LiouvilleCOPS_cuda): .venv/bin/python late_response_figures.py --eval <job>/eval --out <prefix>
#   --pair c3,c3_last5_w10_constant --label c3="cycle 3" ... [--overwrite]
# ==== end of header ====
"""Figures of the late-window forward check (model_ext/late_response_eval.py; user 2026-11-07: COPS time response
"both" windows, gauge-offset sensitivity, "two 12M best fits with and without extra weight of the last point
observations"). Model against the COPS observations only.

Run from <AMPS>/srcCG/LiouvilleCOPS_cuda with .venv/bin/python (cops_cuda.figures: style, palette, save):

    late_response_figures.py --eval <eval dir> --out <prefix> --pair c3,c3_last5_w10_constant \
        --label c3="cycle 3" --label ... [--overwrite]

Writes (.eps and .png, 300 dpi):
* <prefix>_{RG,NG}_{fit,finetuned}: rows = gauge offset along +Z (one per base), columns = the --pair models; in each
  panel the COPS data (black dots), the model at every 1-s grid point (grey line: the value at the time stamp, as the
  fits use it), and the model averaged as the gauge measures: 10-s running average (W10) and 10 s + DPU 5-readout
  average (W10x5), at the observation times; rms of ln(model/obs) in the legend, over the late-window observations
  whose 18-s window starts after the grain onset (10:33:10; the first 5 observations are left out).
* <prefix>_allmodels_dz0_{fit,finetuned}: rows NG, RG; columns inst / W10 / W10x5; every model on the dz0 base.
x axis: altitude of the trajectory point (the dz0 base) at the same times, log, descent to the right.
"""

import argparse
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path.cwd()))
from cops_cuda import figures  # noqa: E402

WLAB = {"inst": "at the time stamp", "W10": "10-s running average", "W10x5": "10 s + DPU 5-readout average"}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--eval", required=True, type=Path)
    ap.add_argument("--out", required=True)
    ap.add_argument("--pair", required=True)
    ap.add_argument("--label", action="append", default=[])
    ap.add_argument("--overwrite", action="store_true")
    a = ap.parse_args()
    labels = dict(s.split("=", 1) for s in a.label)
    bases = sorted([p.stem[len("series_"):] for p in a.eval.glob("series_*.npz")], key=lambda s: float(s[2:]))
    S = {b: np.load(a.eval / ("series_%s.npz" % b), allow_pickle=True) for b in bases}
    z0 = S["dz0"]
    kind = z0["kind"]
    grid, obs = np.nonzero(kind == "grid")[0], z0["obs_index"]
    xg, xo = z0["alt"][grid], z0["alt"][obs]
    models = sorted({k.split("|")[0].rsplit(":", 1)[0] for k in z0.files if "|" in k})
    pair = a.pair.split(",")
    plt = figures.style()
    late = "late visibility: %s–%s UTC, %.0f–%.1f m" % (str(z0["utc"][obs[0]])[11:19], str(z0["utc"][obs[-1]])[11:19], xo[0], xo[-1])

    def axis(ax, ylab=True, gauge="RG"):
        ax.set_xscale("log"); ax.set_yscale("log")
        ax.set_xlim(xo[0] * 1.15, xo[-1] * 0.85)
        if ylab:
            ax.set_ylabel("%s pressure [Pa]" % gauge)

    # rms over the observations whose longest window (18 s) starts after the grain onset (the first post-gap NG point,
    # 10:33:10): before it the model has no grain term, so earlier windows mix grain-free and grain values
    # (an artifact of the step-like grain, NG; RG is treated the same way so both use the same points)
    after = z0["stop_et"] - 18.0 >= z0["et"][obs[0]]

    def rms(model, gauge, w, b):
        r = np.log(S[b]["%s|%s|%s" % (model, gauge, w)] / S[b]["ref" + gauge])
        return np.sqrt(np.nanmean(r[after] ** 2))

    with open(Path(a.out).parent / "summary_after_grain_onset.txt", "w") as f:
        f.write("# rms of ln(model/obs) over the %d late observations from %s UTC (18-s window after the grain onset);\n"
                "# NG also without the readings < 8e-6 Pa (low readings at emission-current jumps)\n"
                "# base model window rms_NG rms_NG_without_low rms_RG\n" % (after.sum(), str(z0["utc"][obs[after][0]])[11:19]))
        for b in bases:
            low = S[b]["refNG"] < 8e-6
            for mdl in models:
                for state in ("fit", "finetuned"):
                    for w in ("inst", "W4", "W10", "W10x5"):
                        nm = "%s:%s" % (mdl, state)
                        rN = np.log(S[b]["%s|NG|%s" % (nm, w)] / S[b]["refNG"])
                        f.write("%-6s %-40s %-6s %.4f %.4f %.4f\n" % (b, nm, w, np.sqrt(np.mean(rN[after] ** 2)),
                                np.sqrt(np.mean(rN[after & ~low] ** 2)), rms(nm, "RG", w, b)))
    for state in ("fit", "finetuned"):
        for gauge in ("RG", "NG"):
            fig, axs = plt.subplots(len(bases), len(pair), figsize=(figures.DOUBLE_COLUMN_IN, 2.0 * len(bases) + 0.6),
                                    sharex=True, sharey=True, squeeze=False)
            for i, b in enumerate(bases):
                for j, mdl in enumerate(pair):
                    ax, nm = axs[i, j], "%s:%s" % (mdl, state)
                    d = S[b]
                    ax.plot(xg, d["%s|%s" % (nm, gauge)][grid], "-", color="0.6", lw=0.8, zorder=2,
                            label="model %s (rms %.3f)" % (WLAB["inst"], rms(nm, gauge, "inst", b)))
                    for w, c, mk in (("W10", figures.PALETTE[0], "s"), ("W10x5", figures.PALETTE[1], "^")):
                        ax.plot(xo, d["%s|%s|%s" % (nm, gauge, w)], "-", marker=mk, ms=2.2, mew=0, lw=0.9, color=c, zorder=4,
                                label="model, %s (rms %.3f)" % (WLAB[w], rms(nm, gauge, w, b)))
                    ok = d["ref" + gauge] > 0
                    ax.plot(xo[ok], d["ref" + gauge][ok], "o", color="k", ms=1.8, mew=0, zorder=5, label="COPS data")
                    axis(ax, j == 0, gauge)
                    ax.set_title("%s, %s; gauge offset +Z %s m" % (labels.get(mdl, mdl), "fit" if state == "fit" else "fine-tuned",
                                                                  b[2:]), fontsize=7)
                    ax.legend(loc="lower left" if gauge == "RG" else "upper left", fontsize=5.5, handlelength=2.0)
                    if i == len(bases) - 1:
                        ax.set_xlabel("altitude of the trajectory point [m] (descent →)")
            fig.suptitle("SHAP8 12M, %s: model vs COPS with the gauge time response (%s)\nrms of ln(model/obs) from %s UTC (windows after the grain onset)"
                         % (gauge, late, str(z0["utc"][obs[after][0]])[11:19]), fontsize=8)
            fig.tight_layout(rect=(0, 0, 1, 0.97))
            print("wrote", figures.save(fig, "%s_%s_%s" % (a.out, gauge, state), overwrite=a.overwrite))
            plt.close(fig)
        fig, axs = plt.subplots(2, 3, figsize=(figures.DOUBLE_COLUMN_IN, 5.4), sharex=True, sharey="row", squeeze=False)
        for i, gauge in enumerate(("NG", "RG")):
            for j, w in enumerate(("inst", "W10", "W10x5")):
                ax = axs[i, j]
                for k, mdl in enumerate(models):
                    nm = "%s:%s" % (mdl, state)
                    ax.plot(xo, z0["%s|%s|%s" % (nm, gauge, w)], "-", color=figures.PALETTE[k], lw=0.9,
                            label=labels.get(mdl, mdl))
                ok = z0["ref" + gauge] > 0
                ax.plot(xo[ok], z0["ref" + gauge][ok], "o", color="k", ms=1.8, mew=0, zorder=5, label="COPS data")
                axis(ax, j == 0, gauge)
                ax.set_title("%s: %s" % (gauge, WLAB[w]), fontsize=7)
                # rms column in the empty corner: NG upper left, RG lower left (header "rms" on top)
                col = [("rms", "0.3")] + [("%.3f" % rms("%s:%s" % (mdl, state), gauge, w, "dz0"), figures.PALETTE[k])
                                          for k, mdl in enumerate(models)]
                for k, (txt, c) in enumerate(col):
                    y = 0.97 - 0.07 * k if gauge == "NG" else 0.03 + 0.07 * (len(col) - 1 - k)
                    ax.text(0.02, y, txt, color=c, transform=ax.transAxes, fontsize=5.5, va="top" if gauge == "NG" else "bottom")
                if i == 1:
                    ax.set_xlabel("altitude of the trajectory point [m]")
        fig.suptitle("SHAP8 12M, all models (%s), gauge at the trajectory point (%s)\nnumbers: rms of ln(model/obs) from %s UTC (windows after the grain onset)"
                     % ("fits" if state == "fit" else "fine-tuned", late, str(z0["utc"][obs[after][0]])[11:19]), fontsize=8)
        h, l = axs[0, 0].get_legend_handles_labels()
        fig.legend(h, l, loc="lower center", ncol=3, fontsize=6.5, frameon=False)
        fig.tight_layout(rect=(0, 0.07, 1, 0.95))
        print("wrote", figures.save(fig, "%s_allmodels_dz0_%s" % (a.out, state), overwrite=a.overwrite))
        plt.close(fig)


if __name__ == "__main__":
    main()
