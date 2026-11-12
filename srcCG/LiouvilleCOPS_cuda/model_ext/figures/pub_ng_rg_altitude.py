#!/usr/bin/env python3
# ==== version-controlled copy: header added 2026-11-12, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/figures_cluster/scripts/pub_ng_rg_altitude.py (on Voyager: /data001/highresolu_vt/cops_output/figures_cluster/scripts/pub_ng_rg_altitude.py)
# SHA256 of the run copy = SHA256 of this file without these header lines: e47c80c3cc586ee5fc85b097dc56cef3653b6b2cc3603945a3db734e14345e9b
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: figure script. Publication figure (eps + png, 300 dpi): COPS NG / RG and one model's NG / RG vs altitude in one
#   panel, legend on the right, model labelled 'Model NG/RG', the gap in the NG data hatched, the model NG = RG crossing
#   (user requests 2026-11-12). Input: a positive_table.py table of a fit's observation rows.
# RUN: pub_ng_rg_altitude.py --table <table .dat> --out <prefix>   (login node; from srcCG/LiouvilleCOPS_cuda with
#   .venv/bin/python) -> figures_cluster/12M/publication_softplus/12M_softplus_NG_RG_altitude.{eps,png}
# ==== end of header ====
"""Publication figure: COPS NG and RG pressures and one model's NG and RG curves as functions of altitude, in one
panel (as panel (c)/(e) "NG and RG" of model_obs_figures.py), legend to the right of the axes.

User requests 2026-11-12: "create publication quality eps and png plot with comparison of RG and NG COPS and softplus
curves as a function of altitude (like the third panel of the data/model comparison plot). place the legend on the
right side of the plot; label the model data as AMPS NG/RG"; then "modify the plot: 1) replace AMPS RG/NG with Model
RG/NG; 2) remove colored 'late visibility' panel; 3) simply mark the gap interval and say that is the gap in NG data".

Run from <AMPS>/srcCG/LiouvilleCOPS_cuda with .venv/bin/python (uses cops_cuda.figures for style, data and saving):

    pub_ng_rg_altitude.py --table <model table .dat> --out <prefix> [--model-label Model] [--overwrite]

<model table .dat>: the observation rows of a fit (figures_cluster/scripts/positive_table.py output: COPS_model.dat
columns, model-only grid rows removed). Writes <prefix>.eps and <prefix>.png (300 dpi).
Content: COPS NG (black) and RG (grey) observations (filled: fitted; open: not fitted); the model NG (solid) and RG
(dashed); the model's NG = RG crossing (open circle); the gap in the NG data (hatched): every interval between two
consecutive NG observations (descent order) whose altitudes differ by more than a factor GAP_RATIO and inside which
RG observations exist. Altitude axis logarithmic, descent from left to right.
"""
import argparse
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path.cwd()))
from cops_cuda import figures  # noqa: E402

RG_DASH = (0, (4, 2))
RG_GREY = "0.55"
NG_COLOR = "#0072B2"     # model NG (blue)
RG_COLOR = "#D55E00"     # model RG (vermillion)
GAP_RATIO = 2.0          # an NG gap: consecutive NG altitudes differ by more than this factor, with RG data inside


def crossings(x, lr):
    out = []
    for i in np.nonzero(np.sign(lr[:-1]) * np.sign(lr[1:]) < 0)[0]:
        f = lr[i] / (lr[i] - lr[i + 1])
        out.append((x[i] + f * (x[i + 1] - x[i]), x[i], x[i + 1]))
    return out


def ng_gaps(tab, alt):
    """[(alt_hi, alt_lo, utc_hi, utc_lo)] of the gaps in the NG data during which RG was measured."""
    o = np.argsort(-alt)
    ng = o[tab["refNG"][o] > 0]
    rg_alt = alt[tab["refRG"] > 0]
    out = []
    for i in range(len(ng) - 1):
        hi, lo = alt[ng[i]], alt[ng[i + 1]]
        if hi / lo > GAP_RATIO and ((rg_alt < hi) & (rg_alt > lo)).any():
            out.append((hi, lo, tab["time"][ng[i]], tab["time"][ng[i + 1]]))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--table", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--model-label", default="Model")
    ap.add_argument("--overwrite", action="store_true")
    a = ap.parse_args()

    plt = figures.style()
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch
    tab = figures.read_table(a.table)
    alt = tab["altitude"]
    order = np.argsort(-alt)                          # descent order

    fig, ax = plt.subplots(figsize=(6.6, 3.6))
    fig.subplots_adjust(left=0.10, right=0.70, bottom=0.14, top=0.95)
    gaps = ng_gaps(tab, alt)
    for hi, lo, *_ in gaps:
        ax.axvspan(lo, hi, facecolor="none", edgecolor="0.45", hatch="////", lw=0, zorder=1)
    # observations
    for gauge, color in (("NG", "k"), ("RG", RG_GREY)):
        obs, use = tab["ref" + gauge], tab["use" + gauge]
        ok = obs > 0
        ax.plot(alt[ok & (use > 0)], obs[ok & (use > 0)], "o", color=color, ms=2.0, mew=0, zorder=3)
        if (ok & (use == 0)).any():
            ax.plot(alt[ok & (use == 0)], obs[ok & (use == 0)], "o", mfc="none", mec=color, ms=2.6, mew=0.5, zorder=3)
    # model
    ax.plot(alt[order], tab["pNG"][order], "-", color=NG_COLOR, lw=1.3, zorder=10)
    ax.plot(alt[order], tab["pRG"][order], ls=RG_DASH, color=RG_COLOR, lw=1.3, zorder=11)
    mc = crossings(np.log(alt[order]), np.log(tab["pNG"][order] / tab["pRG"][order]))
    for c, *_ in mc:
        p = np.exp(np.interp(c, np.log(alt[order][::-1]), np.log(tab["pNG"][order][::-1])))
        ax.plot([np.exp(c)], [p], "o", mfc="white", mec="k", mew=1.0, ms=5.0, zorder=30)

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(alt.max() * 1.3, max(alt.min(), 1e-3) * 0.7)
    ax.set_xlabel("altitude above the surface [m]")
    ax.set_ylabel("pressure [Pa]")

    m = a.model_label
    handles = [Line2D([], [], ls="none", marker="o", color="k", ms=3.5, mew=0),
               Line2D([], [], ls="none", marker="o", color=RG_GREY, ms=3.5, mew=0),
               Line2D([], [], color=NG_COLOR, lw=1.4), Line2D([], [], color=RG_COLOR, lw=1.4, ls=RG_DASH),
               Line2D([], [], ls="none", marker="o", mfc="white", mec="k", mew=1.0, ms=5.0)]
    labels = ["COPS NG", "COPS RG", "%s NG" % m, "%s RG" % m, "%s NG = RG" % m]
    for hi, lo, t0, t1 in gaps:
        handles.append(Patch(facecolor="none", edgecolor="0.45", hatch="////", lw=0))
        labels.append("gap in NG data\n%.0f–%.0f m\n(%s–%s UTC)" % (hi, lo, t0.strftime("%H:%M:%S"), t1.strftime("%H:%M:%S")))
    ax.legend(handles, labels, loc="center left", bbox_to_anchor=(1.02, 0.5), frameon=False, fontsize=8,
              handlelength=2.4, borderaxespad=0.0, labelspacing=0.9)
    figures.save(fig, a.out, overwrite=a.overwrite)
    plt.close(fig)
    print("wrote %s.eps, %s.png; model NG = RG at %s m; NG gaps: %s" % (
        a.out, a.out, ", ".join("%.1f" % np.exp(c[0]) for c in mc) or "none",
        ", ".join("%.1f-%.1f m (%s-%s)" % (hi, lo, t0.strftime("%H:%M:%S"), t1.strftime("%H:%M:%S"))
                  for hi, lo, t0, t1 in gaps) or "none"))


if __name__ == "__main__":
    main()
