#!/usr/bin/env python3
# ==== version-controlled copy: header added 2026-11-21, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/figures_cluster/scripts/q002_top5_sources_figures.py (on Voyager: /data001/highresolu_vt/cops_output/figures_cluster/scripts/q002_top5_sources_figures.py)
# SHA256 of the run copy = SHA256 of this file without these header lines: 1f37e39a2ba7b822d727866ab0783a1f33f6921a346c156dc08cfa847fa68035
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: source figures of the five most plausible new-law fits for section 4 of 12M_q002_top5_models.docx: Q over
#   the day, slope-bin production and pressure shares, lobes / distance / source areas, illumination shares (2026-11-21)
# RUN: python q002_top5_sources_figures.py --out figures_cluster/12M/q002_top5_models (JSON inputs; seconds)
# ==== end of header ====
"""Source and comet figures of the five most plausible new-law fits (2026-11-21), for section 4 of
12M/q002_top5_models/12M_q002_top5_models.docx (user 2026-11-21: "enhance the summary: that also should include ground
tracks, source rates, quantification of the contribution of different area to the observations, importance of the
cliffs, etc -- all useful information about the comet, sources, etc").

Reads only recorded numbers: full_set_q002_<run>/checks/facts.json (Q over the day, seen / unseen production, slope-bin
shares of the seen area and production, pressure shares by slope bin, P50 / P90 source-area unions, lobe and distance
shares) and source_rate_q002/<run>/inspection/checks/inspect_<hhmmss>.json (day / night / cast-shadow shares of Q).
Writes <out>/figures_sources/ (.png + .eps 300 dpi + _caption.txt) and <out>/numbers_sources.json. Never overwrites.

    q002_top5_sources_figures.py --out figures_cluster/12M/q002_top5_models
"""
import argparse
import hashlib
import json
import sys
from pathlib import Path

import numpy as np

APP = Path("/data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda")
sys.path.insert(0, str(APP))
import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from cops_cuda import figures  # noqa: E402

F12 = Path("/data001/highresolu_vt/cops_output/figures_cluster/12M")
SR = F12 / "source_rate_q002"
RUNS = [("C_f10_pm0p01_q002", "C_f10_pm0p01", "C_f10_pm0p01_q002_reg", "C f10, P_m 0.01, templates on"),
        ("C_f1_pm0p01_q002", "C_f1_pm0p01", "C_f1_pm0p01_q002_reg", "C f1, P_m 0.01, templates on"),
        ("C_f10_pm0p01_tploff_q002", "C_f10_pm0p01_tploff", "C_f10_pm0p01_tploff_q002_reg", "C f10, P_m 0.01, templates off"),
        ("A_q002", "A", "A_q002_reg", "A sum of squares, rank 3"),
        ("B_q002", "B", "B_q002_reg", "B background, bg 1e-2")]
BINS = ["[0,40)", "[40,50)", "[50,60)", "[60,70)", "[70,80)", "[80,90)", "[90,180]"]
BIN_LAB = ["0-40°", "40-50°", "50-60°", "60-70°", "70-80°", "80-90°", "90-180°"]
# the colours of the full sets' cliff figures (one per bin)
BIN_COL = ["#9e9e9e", "#E69F00", "#D55E00", "#0b4f8a", "#56B4E9", "#009E73", "#7b3294"]
LAW = "new H2O law q_min = 0.02 q_max"
INK, MUTED, GRID = "#0b0b0b", "#52514e", "#e4e3df"


def sha256(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def style(ax, axis="both"):
    ax.grid(True, axis=axis, color=GRID, lw=0.5)
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
    FIG = O / "figures_sources"
    if FIG.exists() or (O / "numbers_sources.json").exists():
        raise SystemExit("figures_sources/ or numbers_sources.json exists (never overwritten)")
    FIG.mkdir(parents=True)
    num = {"law": LAW, "sources": {}, "runs": {}}
    for key, fk, sk, lab in RUNS:
        fp = F12 / ("full_set_q002_%s/checks/facts.json" % fk)
        f = json.load(open(fp))
        num["sources"]["facts_" + fk] = [str(fp), sha256(fp)]
        ill = {}
        for t in ("103310", "103926"):
            ip = SR / sk / "inspection" / "checks" / ("inspect_%s.json" % t)
            d = json.load(open(ip))["new"]
            num["sources"]["inspect_%s_%s" % (sk, t)] = [str(ip), sha256(ip)]
            ill[t] = {k: d[k] for k in ("Q_H2O_s-1", "seen_share_of_Q", "dark_share_of_Q", "night_share_of_Q",
                                        "cast_shadow_share_of_Q", "dark_share_of_seen_Q", "dark_area_fraction",
                                        "dark_area_fraction_seen")}
        num["runs"][key] = {
            "label": lab,
            "Q_time": {t: v[0] for t, v in f["Q"].items()},
            "production_seen_unseen": f["production_seen_unseen"],
            "production_share_seen_103310": f["production_share_seen_103310"],
            "seen_area_share": f["seen_area_share"],
            "pressure_share_by_bin": f["pressure_share_by_bin"],
            "source_area_unions": f["source_area_unions"],
            "regions": f["regions"], "illumination": ill, "multipliers": f["multipliers"]}
    R = [num["runs"][k] for k, *_ in RUNS]
    labs = [lab for *_, lab in RUNS]
    col = [figures.PALETTE[i] for i in range(len(RUNS))]
    ls = ["-", (0, (4, 2)), (0, (1, 1.5)), (0, (5, 1.5, 1, 1.5)), "-"]

    # ---- 1. Q over the day
    times = list(R[0]["Q_time"])
    def hours(t):
        h, m, s = (int(x) for x in t.split(":"))
        return h + m / 60 + s / 3600
    x = [hours(t) for t in times]
    fig, ax = plt.subplots(figsize=(8.0, 4.3), dpi=300)
    for i, r in enumerate(R):
        ax.plot(x, [r["Q_time"][t] for t in times], color=col[i], ls=ls[i], lw=1.6, marker="o", ms=3.2, label=labs[i])
    ax.set_yscale("log")
    ax.set_xlabel("time on 2016-09-30 (UTC, h)")
    ax.set_ylabel("Q$_{H_2O}$ (s$^{-1}$)")
    ax.axvspan(hours("10:33:10"), hours("10:39:26"), color=figures.LATE_FILL, lw=0)
    style(ax)
    ax.legend(fontsize=7.2, frameon=False, loc="upper left", ncol=2)
    fig.tight_layout()
    save(fig, FIG, "12M_q002_top5_Q_time",
         "Total H2O production Q_H2O of the five selected fits (%s) at the 13 map times of 2016-09-30 (00:00 ... "
         "10:00 UTC hourly, 10:33:10 and 10:39:26 UTC): the final activity of each fit, held fixed, with the "
         "illumination (night, cast shadow, cos SZA) of each time; shaded: the late window 10:33:10-10:39:26 UTC. "
         "Source: full_set_q002_<run>/checks/facts.json (Q)." % LAW)

    # ---- 2. slope bins: seen area vs seen production
    fig, ax = plt.subplots(figsize=(9.2, 4.4), dpi=300)
    xb = np.arange(len(BINS))
    w = 0.13
    ax.bar(xb - 2.5 * w, [100 * R[0]["seen_area_share"][b] for b in BINS], w, color="#bdbdbd", label="share of the seen area")
    for i, r in enumerate(R):
        ax.bar(xb + (i - 1.5) * w, [100 * r["production_share_seen_103310"][b] for b in BINS], w, color=col[i],
               label="production share: " + labs[i])
    ax.set_xticks(xb)
    ax.set_xticklabels(BIN_LAB)
    ax.set_xlabel("effective gravitational slope")
    ax.set_ylabel("share of the seen surface (%)")
    ax.set_yscale("symlog", linthresh=1)
    ax.set_ylim(0, 130)
    ax.set_yticks([0, 1, 2, 5, 10, 20, 50, 100])
    ax.set_yticklabels(["0", "1", "2", "5", "10", "20", "50", "100"])
    style(ax, "y")
    ax.legend(fontsize=6.8, frameon=False, loc="upper right")
    fig.tight_layout()
    save(fig, FIG, "12M_q002_top5_slope_production",
         "Where the seen H2O production comes from, by effective gravitational slope (%s): grey, the share of the seen "
         "surface area in each bin (fit-independent; 'seen' = the facet contributes to at least one fitted NG or RG "
         "observation, 66.5 %% of the area); coloured, each fit's share of the H2O production of the seen surface at "
         "10:33:10 UTC (facet flux times area). Cliffs (>= 60 deg) are 14.6 %% of the seen area. Vertical axis linear "
         "below 1 %%, logarithmic above. Source: full_set_q002_<run>/checks/facts.json (seen_area_share, "
         "production_share_seen_103310)." % LAW)

    # ---- 3. pressure shares by slope bin
    phases = [("pre_mean", "before the NG gap (mean)"), ("late_mean", "late window (mean)"),
              ("last", "last observation (5.5 m)")]
    fig, axes = plt.subplots(2, 3, figsize=(13.0, 5.6), dpi=300, sharey=True)
    y = np.arange(len(RUNS))[::-1]
    for gi, g in enumerate(("NG", "RG")):
        for pi, (ph, title) in enumerate(phases):
            ax = axes[gi, pi]
            for i, r in enumerate(R):
                left = 0.0
                for bi in range(len(BINS)):
                    v = 100 * r["pressure_share_by_bin"][g][ph][bi]
                    ax.barh(y[i], v, left=left, color=BIN_COL[bi], height=0.62, edgecolor="white", lw=0.4)
                    if v >= 9:
                        ax.text(left + v / 2, y[i], "%.0f" % v, ha="center", va="center", fontsize=6.3,
                                color="white" if bi in (2, 3, 6) else INK)
                    left += v
            ax.set_xlim(0, 100)
            ax.set_title("%s, %s" % (g, title), fontsize=8.6, loc="left")
            style(ax, "x")
            if pi == 0:
                ax.set_yticks(y)
                ax.set_yticklabels(labs, fontsize=7.6)
            if gi == 1:
                ax.set_xlabel("share of the modelled surface pressure (%)")
    handles = [plt.Rectangle((0, 0), 1, 1, color=c) for c in BIN_COL]
    fig.legend(handles, BIN_LAB, loc="lower center", ncol=7, frameon=False, fontsize=8, title="effective gravitational "
               "slope (90-180° = overhangs)", title_fontsize=8)
    fig.tight_layout(rect=(0, 0.09, 1, 1))
    save(fig, FIG, "12M_q002_top5_slope_pressure",
         "Which slopes the gauges' modelled signal comes from (%s): share of the modelled surface pressure of NG (top) "
         "and RG (bottom) contributed by each effective-gravitational-slope bin, as the mean over the observation rows "
         "before the NG gap (10:24:38 UTC), the mean over the late window 10:33:10-10:39:26 UTC (297-5.5 m), and at "
         "the last observation (10:39:26 UTC, 5.5 m). Numbers: shares >= 9 %%. Source: full_set_q002_<run>/checks/"
         "facts.json (pressure_share_by_bin; the per-row curves are the full sets' cliffs/12M_cliff_pressure_fraction "
         "figures)." % LAW)

    # ---- 4. regions and source areas
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(12.4, 4.4), dpi=300)
    cats = [("big lobe, NG before late (median)", lambda r: r["regions"]["big_lobe_NG_median_pre"]),
            ("small lobe, NG late (median)", lambda r: r["regions"]["small_lobe_median_late"][0]),
            ("small lobe, RG late (median)", lambda r: r["regions"]["small_lobe_median_late"][1]),
            ("30-100 m from touchdown, NG late", lambda r: r["regions"]["shell_30_100m_median_late"][0]),
            ("30-100 m from touchdown, RG late", lambda r: r["regions"]["shell_30_100m_median_late"][1])]
    yc = np.arange(len(cats))[::-1]
    for i, r in enumerate(R):
        a1.plot([100 * fn(r) for _, fn in cats], yc + (2 - i) * 0.12, "o", color=col[i], ms=5, label=labs[i])
    a1.set_yticks(yc)
    a1.set_yticklabels([c for c, _ in cats], fontsize=7.6)
    a1.set_xlim(0, 100)
    a1.set_xlabel("share of the modelled surface pressure (%)")
    a1.set_title("(a) lobes and distance from the touchdown point", loc="left", fontsize=9)
    style(a1, "x")
    a1.legend(fontsize=6.6, frameon=False, loc="lower right")
    keys = [("NG_P50", "NG, P50"), ("RG_P50", "RG, P50"), ("NG_P90", "NG, P90"), ("RG_P90", "RG, P90")]
    yk = np.arange(len(keys))[::-1]
    for i, r in enumerate(R):
        a2.plot([r["source_area_unions"][k]["area_m2"] for k, _ in keys], yk + (2 - i) * 0.12, "o", color=col[i], ms=5)
    a2.set_xscale("log")
    a2.set_yticks(yk)
    a2.set_yticklabels([l for _, l in keys], fontsize=8)
    a2.set_xlabel("area of the union of the source areas (m$^2$)")
    a2.set_title("(b) P50 / P90 source areas of the 95 late observations", loc="left", fontsize=9)
    style(a2, "x")
    fig.tight_layout()
    save(fig, FIG, "12M_q002_top5_regions_sources",
         "Where the late signal comes from (%s). (a) Median share of the modelled surface pressure from the big lobe "
         "(body) for NG before the late window, from the small lobe (head, where Rosetta touched down) in the late "
         "window for NG and RG, and from facets 30-100 m from the touchdown point in the late window (lobes separated by "
         "the neck plane). (b) Area of the union over the 95 late observations of each observation's P%% source area "
         "(the smallest set of facets, in decreasing order of contribution, that gives P %% of the observation's modelled "
         "signal), P = 50 and 90, NG and RG; logarithmic. Source: full_set_q002_<run>/checks/facts.json (regions, "
         "source_area_unions)." % LAW)

    # ---- 5. illumination classes of the production
    fig, axes = plt.subplots(1, 2, figsize=(11.6, 3.6), dpi=300, sharey=True)
    for ax, t, tl in ((axes[0], "103310", "10:33:10 UTC"), (axes[1], "103926", "10:39:26 UTC")):
        for i, r in enumerate(R):
            d = r["illumination"][t]
            day = 1 - d["dark_share_of_Q"]
            parts = [(100 * day, "#E69F00", "day (lit)"), (100 * d["night_share_of_Q"], "#0b4f8a", "night"),
                     (100 * d["cast_shadow_share_of_Q"], "#56B4E9", "cast shadow")]
            left = 0
            for v, c, _ in parts:
                ax.barh(y[i], v, left=left, color=c, height=0.62, edgecolor="white", lw=0.4)
                if v >= 7:
                    ax.text(left + v / 2, y[i], "%.0f" % v, ha="center", va="center", fontsize=6.6,
                            color="white" if c != "#E69F00" else INK)
                left += v
        ax.set_xlim(0, 100)
        ax.set_title("share of Q$_{H_2O}$ by illumination, %s" % tl, loc="left", fontsize=9)
        ax.set_xlabel("%")
        style(ax, "x")
    axes[0].set_yticks(y)
    axes[0].set_yticklabels(labs, fontsize=7.6)
    hs = [plt.Rectangle((0, 0), 1, 1, color=c) for c in ("#E69F00", "#0b4f8a", "#56B4E9")]
    fig.legend(hs, ["day (lit)", "night", "cast shadow"], loc="lower center", ncol=3, frameon=False, fontsize=8)
    fig.tight_layout(rect=(0, 0.1, 1, 1))
    save(fig, FIG, "12M_q002_top5_illumination",
         "Share of the total H2O production from lit, night and cast-shadow facets at the two map times (%s: night and "
         "cast-shadow facets produce q_min = 0.02 q_max per unit activity, lit facets q_min + (q_max - q_min) cos SZA). "
         "Source: source_rate_q002/<run>/inspection/checks/inspect_<hhmmss>.json (night_share_of_Q, "
         "cast_shadow_share_of_Q, dark_share_of_Q)." % LAW)

    (O / "numbers_sources.json").write_text(json.dumps(num, indent=1))
    print("wrote", FIG, O / "numbers_sources.json")


if __name__ == "__main__":
    main()
