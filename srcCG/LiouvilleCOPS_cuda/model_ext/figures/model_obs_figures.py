#!/usr/bin/env python3
# ==== version-controlled copy: header added 2026-11-06, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/figures_cluster/scripts/model_obs_figures.py (on Voyager: /data001/highresolu_vt/cops_output/figures_cluster/scripts/model_obs_figures.py)
# SHA256 of the run copy = SHA256 of this file without these header lines: a81096bfb4e72fefd58c07a58b76c34a67319f1e4221ae26851dc830618427c0
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: the plotting program behind every figures_*.sh: fitted model NG/RG pressures against the COPS observations
#   (time, altitude, NG, RG, NG+RG with the NG = RG crossing), late window shaded and labelled, .eps + .png 300 dpi.
# RUN (from $APP, .venv/bin/python): <this> --bundle <bundle npz> --ref <fit dir> --out <prefix>
#        --fit "LABEL=<fit dir or table>[:solid|dashed|dotted|dashdot]" ... [--title ...] [--overwrite]
#   (a LABEL must not contain '=')
# ==== end of header ====
"""Fitted model pressures vs the COPS observations only (no other machine, no ratio panel; user request 2026-10-27).

Run from <AMPS>/srcCG/LiouvilleCOPS_cuda with .venv/bin/python (uses cops_cuda.figures for style, data and saving):

    model_obs_figures.py --bundle <run>/bundle/restart_templates.npz --ref <fit dir> --out <prefix> --title <t>
                         --fit "LABEL=DIR[:solid|dashed|dotted|dashdot]" ... [--overwrite]

Writes (.eps and .png, 300 dpi, late window shaded and labelled with UTC and altitude):
* <prefix>_time      panels (a) NG, (b) RG, (c) NG and RG together; full descent left, last 25 min right
* <prefix>_altitude  the same three panels against altitude (log, descent left to right)
* <prefix>_altitude_NG, _altitude_RG, _altitude_NG_RG  the panels as separate figures; _NG_RG lists the crossings
In (a) and (b) each fit keeps its line style; in (c) a fit has one colour, NG is solid and RG dashed, the data are
black (NG) and grey (RG) dots, and every crossing NG = RG of a model is marked with an open circle. Crossings: sign
changes of ln(pNG/pRG) along the descent, interpolated (log altitude, or time); for the data on a running median of
9 points where both gauges have a fitted value. A data crossing inside a gap of the data is hatched and given as the
bracketing range, not a value.
"""

import argparse
import datetime as dt
import re
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path.cwd()))
from cops_cuda import figures  # noqa: E402

STYLES = {"solid": "-", "dashed": (0, (4, 2)), "dotted": (0, (1, 1.5)), "dashdot": (0, (5, 1.5, 1, 1.5))}
RG_DASH = (0, (4, 2))
RG_GREY = "0.55"
ZOOM_MIN = 25.0


def parse_fit(text):
    label, _, rest = text.partition("=")
    path, _, style = rest.rpartition(":") if rest.rsplit(":", 1)[-1] in STYLES else (rest, "", "solid")
    return label, path, STYLES[style or "solid"]


def crossings(x, lr):
    """(x at the sign change, x before, x after) of lr along x (descent order); x is a linear coordinate."""
    out = []
    for i in np.nonzero(np.sign(lr[:-1]) * np.sign(lr[1:]) < 0)[0]:
        f = lr[i] / (lr[i] - lr[i + 1])
        out.append((x[i] + f * (x[i + 1] - x[i]), x[i], x[i + 1]))
    return out


def fmt_alt(a):
    return "%.0f m" % a if a >= 10 else "%.1f m" % a


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bundle", required=True)
    ap.add_argument("--ref", required=True, help="fit directory whose data selection is drawn")
    ap.add_argument("--out", required=True)
    ap.add_argument("--title", default=None)
    ap.add_argument("--fit", action="append", required=True)
    ap.add_argument("--overwrite", action="store_true")
    args = ap.parse_args()

    plt = figures.style()
    import matplotlib.dates as mdates
    import matplotlib.transforms as mtransforms
    from matplotlib.lines import Line2D
    fits = [(lab, figures.read_table(p), ls, figures.PALETTE[k]) for k, (lab, p, ls) in
            enumerate(parse_fit(f) for f in args.fit)]
    ref = figures.read_table(args.ref)
    window = figures.late_window(ref, bundle=args.bundle)
    alt, t = ref["altitude"], ref["time"]
    order = np.argsort(-alt)                        # descent order
    tsec = np.array([(x - t[0]).total_seconds() for x in t])
    to_time = lambda s: t[0] + dt.timedelta(seconds=float(s))  # noqa: E731
    span = "%s: %s–%s UTC, %s–%s m" % (window.get("name", "late visibility"), window["start"].split("T")[1],
                                        window["stop"].split("T")[1], figures._hm(window["altitude_max_m"]),
                                        figures._hm(window["altitude_min_m"]))

    # crossings in descent order, in log altitude and in time
    def cross(lr_desc, idx):
        ca = crossings(np.log(alt[idx]), lr_desc)
        ct = crossings(tsec[idx], lr_desc)
        return [(np.exp(a), np.exp(hi), np.exp(lo), to_time(s), to_time(s0), to_time(s1))
                for (a, hi, lo), (s, s0, s1) in zip(ca, ct)]
    both = (ref["refNG"] > 0) & (ref["refRG"] > 0) & (ref["useNG"] > 0) & (ref["useRG"] > 0)
    o = order[both[order]]
    data_cross = []
    if len(o) > 9:
        lr = np.log(ref["refNG"][o] / ref["refRG"][o])
        med = np.array([np.median(lr[max(0, i - 4):i + 5]) for i in range(len(lr))])
        data_cross = cross(med, o)
    model_cross = [cross(np.log(f[1]["pNG"][order] / f[1]["pRG"][order]), order) for f in fits]

    def shade(ax, xaxis):
        lo, hi = ((window["altitude_min_m"], window["altitude_max_m"]) if xaxis == "alt"
                  else (window["t0"], window["t1"]))
        ax.axvspan(lo, hi, facecolor=figures.LATE_FILL, edgecolor=figures.LATE_EDGE, lw=0.5, zorder=0)

    def data(ax, gauge, color, x):
        obs, use = ref["ref" + gauge], ref["use" + gauge]
        ok = obs > 0
        ax.plot(x[ok & (use > 0)], obs[ok & (use > 0)], "o", color=color, ms=1.6, mew=0, zorder=3)
        if (ok & (use == 0)).any():
            ax.plot(x[ok & (use == 0)], obs[ok & (use == 0)], "o", mfc="none", mec=color, ms=2.2, mew=0.4, zorder=3)

    def single(ax, gauge, xaxis):
        x = alt if xaxis == "alt" else t
        data(ax, gauge, "k", x)
        for k, (lab, tab, ls, color) in enumerate(fits):
            ax.plot(x, tab["p" + gauge], ls=ls, color=color, lw=1.0, zorder=10 + len(fits) - k)
        ax.set_yscale("log")
        ax.set_ylabel("%s pressure [Pa]" % gauge)
        shade(ax, xaxis)

    def combined(ax, xaxis, key=True):
        x = alt if xaxis == "alt" else t
        data(ax, "NG", "k", x)
        data(ax, "RG", RG_GREY, x)
        for a, hi, lo, tc, t0, t1 in data_cross:
            gap = hi / lo > 1.2
            if gap:
                ax.axvspan(*((lo, hi) if xaxis == "alt" else (t0, t1)), facecolor="none", edgecolor="0.3",
                           hatch="////", lw=0, zorder=1)
            else:
                ax.axvline(a if xaxis == "alt" else tc, color="k", lw=0.5, ls=(0, (1, 2)), zorder=1)
        for k, (lab, tab, _, color) in enumerate(fits):
            z = 10 + len(fits) - k
            ax.plot(x, tab["pNG"], "-", color=color, lw=1.0, zorder=z)
            ax.plot(x, tab["pRG"], ls=RG_DASH, color=color, lw=1.0, zorder=z)
            for a, *_r, tc, _t0, _t1 in model_cross[k]:
                p = np.exp(np.interp(np.log(a), np.log(alt[order][::-1]), np.log(tab["pNG"][order][::-1])))
                ax.plot([a if xaxis == "alt" else tc], [p], "o", mfc="white", mec=color, mew=1.0, ms=4.5, zorder=30)
        ax.set_yscale("log")
        ax.set_ylabel("NG and RG pressure [Pa]")
        shade(ax, xaxis)
        if not key:
            return
        ax.add_artist(ax.legend([Line2D([], [], color="0.2", lw=1.0), Line2D([], [], color="0.2", lw=1.0, ls=RG_DASH),
                   Line2D([], [], ls="none", marker="o", color="k", ms=2.5, mew=0),
                   Line2D([], [], ls="none", marker="o", color=RG_GREY, ms=2.5, mew=0),
                   Line2D([], [], ls="none", marker="o", mfc="white", mec="0.2", ms=4.5)],
                  ["model NG", "model RG", "NG data", "RG data", "model NG = RG"],
                  fontsize=5.8, ncol=1, frameon=True, framealpha=1.0, edgecolor="0.7", handlelength=2.2,
                  loc="lower center"))   # add_artist: a later ax.legend() must not replace it

    def alt_axis(ax, label=True):
        ax.set_xscale("log")
        ax.set_xlim(alt.max() * 1.3, max(alt.min(), 1e-3) * 0.7)
        if label:
            ax.set_xlabel("altitude above the shape model [m] (descent →)")

    def window_title(ax):
        tr = mtransforms.blended_transform_factory(ax.transData, ax.transAxes)
        ax.text(np.sqrt(window["altitude_min_m"] * window["altitude_max_m"]), 1.02, span, transform=tr,
                ha="center", va="bottom", fontsize=6.5, color="0.15",
                bbox={"facecolor": figures.LATE_FILL, "edgecolor": figures.LATE_EDGE, "lw": 0.5, "pad": 1.5})

    fit_handles = [Line2D([], [], ls="none", marker="o", color="k", ms=2.5, mew=0)]
    fit_labels = ["COPS data (fitted)"]
    fit_handles += [Line2D([], [], color=c, ls=ls, lw=1.2) for _, _, ls, c in fits]
    fit_labels += [lab for lab, _, _, _ in fits]

    shown = max(len(re.sub(r"\$[^$]*\$", "x", s)) for s in fit_labels)   # label length as drawn
    ncol = 1 if shown > 55 else 2
    nrow = -(-len(fit_labels) // ncol)

    def top_legend(fig, y):
        fig.legend(fit_handles, fit_labels, loc="upper center", bbox_to_anchor=(0.5, y), ncol=ncol, fontsize=7,
                   handlelength=2.6, columnspacing=1.4)
        if args.title:
            fig.suptitle(args.title, y=y + 0.025, fontsize=8)

    # ---- time: full descent and the last minutes
    fig, ax = plt.subplots(3, 2, figsize=(figures.DOUBLE_COLUMN_IN, 7.2), sharex="col", sharey="row",
                           gridspec_kw={"width_ratios": [2.2, 1.0], "wspace": 0.06, "hspace": 0.08})
    fig.subplots_adjust(top=0.955 - 0.026 * nrow)
    for col in (0, 1):
        single(ax[0, col], "NG", "time")
        single(ax[1, col], "RG", "time")
        combined(ax[2, col], "time", key=col == 1)   # line-style key only in the zoom panel
    for r in range(3):
        ax[r, 1].set_ylabel("")
        ax[r, 1].set_xlim(window["t1"] - dt.timedelta(minutes=ZOOM_MIN), window["t1"] + dt.timedelta(minutes=0.6))
    ax[0, 0].set_xlim(t.min() - dt.timedelta(minutes=10), window["t1"] + dt.timedelta(minutes=0.6))
    for a in (ax[2, 0], ax[2, 1]):
        a.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
    ax[2, 1].xaxis.set_major_locator(mdates.MinuteLocator(byminute=range(0, 60, 5)))
    ax[2, 0].set_xlabel("UTC, %s" % t[0].strftime("%Y-%m-%d"))
    ax[2, 1].set_xlabel("UTC (last %.0f min)" % ZOOM_MIN)
    ax[1, 1].text(0.04, 0.05, "shaded: " + figures.window_label(window), transform=ax[1, 1].transAxes,
                  va="bottom", ha="left", fontsize=6.3, color="0.15",
                  bbox={"facecolor": figures.LATE_FILL, "edgecolor": figures.LATE_EDGE, "lw": 0.5, "pad": 1.5})
    for k, a in enumerate(ax.T.ravel()):
        a.text(0.02 if k < 3 else 0.04, 0.95, "(%s)" % "acebdf"[k], transform=a.transAxes, va="top", fontsize=8)
    top_legend(fig, 0.975)
    figures.save(fig, args.out + "_time", overwrite=args.overwrite)
    plt.close(fig)

    # ---- altitude, three panels
    fig, ax = plt.subplots(3, 1, figsize=(figures.DOUBLE_COLUMN_IN * 0.72, 7.6), sharex=True,
                           gridspec_kw={"hspace": 0.08})
    fig.subplots_adjust(top=0.935 - 0.025 * nrow)
    single(ax[0], "NG", "alt")
    single(ax[1], "RG", "alt")
    combined(ax[2], "alt")
    for r, a in enumerate(ax):
        alt_axis(a, r == 2)
        a.text(0.01, 0.95, "(%s)" % "abc"[r], transform=a.transAxes, va="top", fontsize=8)
    window_title(ax[0])
    top_legend(fig, 0.975)
    figures.save(fig, args.out + "_altitude", overwrite=args.overwrite)
    plt.close(fig)

    # ---- altitude, the panels as separate figures
    for name, draw in (("NG", lambda a: single(a, "NG", "alt")), ("RG", lambda a: single(a, "RG", "alt")),
                       ("NG_RG", lambda a: combined(a, "alt"))):
        fig, ax = plt.subplots(figsize=(figures.DOUBLE_COLUMN_IN * 0.72, 4.2))
        draw(ax)
        alt_axis(ax)
        window_title(ax)
        if name == "NG_RG":
            notes = ["data: " + (", ".join("between %s and %s (no NG data between)" % (fmt_alt(hi), fmt_alt(lo))
                                           if hi / lo > 1.2 else fmt_alt(a) for a, hi, lo, *_ in data_cross)
                                 or "none")]
            notes += ["%s: %s" % (lab, ", ".join(fmt_alt(c[0]) for c in mc) or "none")
                      for (lab, *_), mc in zip(fits, model_cross)]
            # raise the top of the (log) y axis so that the note box (top left) covers no data or model curve
            lo, hi = ax.get_ylim()
            frac = min(0.6, 0.032 * (len(notes) + 1) + 0.03)
            ax.set_ylim(lo, np.exp(np.log(lo) + (np.log(hi) - np.log(lo)) / (1.0 - frac)))
            ax.text(0.01, 0.98, "NG = RG at (descent order)\n" + "\n".join(notes), transform=ax.transAxes,
                    va="top", ha="left", fontsize=5.8, color="0.15",
                    bbox={"facecolor": "white", "edgecolor": "0.7", "lw": 0.4, "pad": 2}, zorder=40)
        ax.legend(fit_handles, fit_labels, loc="upper center", bbox_to_anchor=(0.5, -0.13),
                  ncol=ncol, fontsize=7, handlelength=2.6, frameon=False)
        if args.title:
            ax.set_title(args.title, fontsize=8, pad=16)
        figures.save(fig, "%s_altitude_%s" % (args.out, name), overwrite=args.overwrite)
        plt.close(fig)


if __name__ == "__main__":
    main()
