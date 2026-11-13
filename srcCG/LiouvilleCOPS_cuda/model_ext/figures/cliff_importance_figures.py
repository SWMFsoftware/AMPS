#!/usr/bin/env python3
# ==== version-controlled copy: header added 2026-11-13, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/figures_cluster/scripts/cliff_importance_figures.py (on Voyager: /data001/highresolu_vt/cops_output/figures_cluster/scripts/cliff_importance_figures.py)
# SHA256 of the run copy = SHA256 of this file without these header lines: 10d31930ceb89464c12ef2e283386e8267f4ca48078c7f71687c0461ce54b94b
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: figure script: how the fit uses the effective-gravitational-slope bins (areas, multipliers, production
#   shares; stacked pressure shares along the descent); 2026-11-13: --model-name, --fit-dir (bounds), --maps-dir (CO2
#   share), --groups default|all, data-driven caption statements.
# RUN: jobs/pub_source_areas.sbatch STEPS=cliffs (V2 = folder with checks/slope_bins_total_*.json)
# ==== end of header ====
"""Cliff-importance figures of one fitted state (softplus model after 5 cycles): how the fit uses the effective-
gravitational-slope bins (user request 2026-11-12: "as the fitting shows importance of the cliffs, also add a plot that
illustrates that"; lead: the fraction of the modelled NG and RG surface pressure from each slope bin along the descent,
the gap and the late window marked, and per bin the seen-area fraction, the fitted multiplier and the share of the seen
production).

Input: <chain>/contributions/slope_bins.npz (model_ext nonneg_contrib.py, job 124549: per bin multiplier, area, seen
area, production shares at 10:33:10 and 10:39:26; per base row the slope-bin sums and fractions of the modelled NG and
RG surface pressure, kind obs / grid / late_obs, use flags), the fit base's COPS_model.dat (altitude of the rows, the
data for the NG gap), the bundle (late window); cross-check of the production shares against <v2>/checks/
slope_bins_total_*.json (pub_maps_v2_figures.py, total rate from the verified fields).

Figures (each .eps, .png 300 dpi, _caption.txt; refuses existing files):
  <out>/12M_cliff_slope_bins                 (a) area shares and seen fraction, (b) multiplier, (c) production shares
  <out>/12M_cliff_pressure_fraction_{NG,RG,NG_RG}  fraction of the modelled surface pressure from each slope bin vs
                                             altitude along the descent (observation rows; stacked)
Run from <AMPS>/srcCG/LiouvilleCOPS_cuda with .venv/bin/python (CPU; reads small files):
    cliff_importance_figures.py --contrib-dir <chain>/contributions --base <fit base> --bundle <bundle .npz>
        --v2-dir <v2 folder> --out-dir <v2 folder>/cliffs
"""

import argparse
import json
import sys
from pathlib import Path

import numpy as np

APP = Path("/data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda")
sys.path.insert(0, str(APP))
from cops_cuda import figures as FIG  # noqa: E402

# groups of the stacked descent figures: (label, bins, colour, colour name for the caption). Default = the grouping of the
# unbounded softplus fit of 2026-11-12 (40-60 and 70-90 deg switched off there); --groups chooses another (2026-11-13).
GROUPS = (("0–40°", (0,), "#BBBBBB", "grey"), ("60–70°", (3,), "#005A8C", "dark blue"),
          ("90–180° (overhangs)", (6,), "#CC79A7", "purple"), ("40–60° and 70–90°", (1, 2, 4, 5), "#E69F00", "orange"))
GROUP_SETS = {
    "default": GROUPS,
    # every bin on its own (2026-11-13, refit with multipliers >= 1: no bin is switched off)
    "all": (("0–40°", (0,), "#BBBBBB", "grey"), ("40–50°", (1,), "#F0E442", "yellow"),
            ("50–60°", (2,), "#E69F00", "orange"), ("60–70°", (3,), "#005A8C", "dark blue"),
            ("70–80°", (4,), "#56B4E9", "light blue"), ("80–90°", (5,), "#009E73", "green"),
            ("90–180° (overhangs)", (6,), "#CC79A7", "purple")),
}
GAP_RATIO = 2.0


def ng_gaps(tab, alt):
    """Gaps in the NG data during which RG was measured (as figures_cluster/scripts/pub_ng_rg_altitude.py)."""
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
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    for a in ("--contrib-dir", "--base", "--bundle", "--v2-dir", "--out-dir"):
        ap.add_argument(a, required=True)
    ap.add_argument("--model-name", default="the softplus activity model, final state after 5 fit and fine-tuning "
                    "cycles", help="name of the fitted state in the captions (2026-11-13: set per fit)")
    ap.add_argument("--fit-dir", default=None, help="c<N>/fit/final of the state: lower bounds of the multipliers")
    ap.add_argument("--maps-dir", default=None, help="maps of the state (maps_manifest.json): Q_CO2 / Q_H2O")
    ap.add_argument("--groups", default="default", choices=sorted(GROUP_SETS),
                    help="slope-bin groups of the descent figures (default: the 2026-11-12 grouping)")
    args = ap.parse_args()
    global GROUPS
    GROUPS = GROUP_SETS[args.groups]
    cman = json.load(open(Path(args.contrib_dir) / "contrib_manifest.json"))
    cjob = cman.get("slurm_job_id") or "?"
    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    plt = FIG.style()
    from matplotlib.patches import Patch
    z = np.load(Path(args.contrib_dir) / "slope_bins.npz", allow_pickle=False)
    labs = [str(s).replace(",", "–").strip("[]()") for s in z["bin_labels_deg"]]
    labs = ["%s°" % s for s in labs]
    nb = len(labs)
    x = np.arange(nb)
    mult = np.asarray(z["multiplier"])
    times = ("103310", "103926")
    # data-driven statements (2026-11-13): bounds and switched-off bins of the multipliers, CO2 share of the production
    bound_txt = ""
    if args.fit_dir:
        with np.load(Path(args.fit_dir) / "fit_restart.npz", allow_pickle=False) as zr:
            rn = [str(n) for n in json.loads(str(zr["meta_json"]))["names"]]
            rx, rlo = np.asarray(zr["x"], np.float64), np.asarray(zr["lo"], np.float64)
        mb = [(n, rx[i], rlo[i]) for i, n in enumerate(rn) if n.startswith("c_slope_")]
        lows = sorted({lo for _, _, lo in mb})
        at = [labs[1 + k] for k, (n, x, lo) in enumerate(mb) if x - lo <= 1e-9 * max(1.0, abs(lo))]
        bound_txt = (" The %d free multipliers have the lower bound %s; at that bound: %s." % (
            len(mb), " / ".join("%g" % v for v in lows), ", ".join(at) if at else "none"))
    offb = [k for k in range(nb) if mult[k] < 1e-10]
    off_txt = (" The bins %s are switched off by their multipliers (<= %.2g)." % (
        ", ".join(labs[k] for k in offb), mult[offb].max()) if offb else "")
    co2_txt = ""
    if args.maps_dir:
        mt = {t["utc"]: t for t in json.load(open(Path(args.maps_dir) / "maps_manifest.json"))["times"]}
        r = [mt["2016-09-30T10:33:10"]["Q_CO2_s-1"] / mt["2016-09-30T10:33:10"]["Q_H2O_s-1"],
             mt["2016-09-30T10:39:26"]["Q_CO2_s-1"] / mt["2016-09-30T10:39:26"]["Q_H2O_s-1"]]
        co2_txt = "; the CO2 production, %.1e and %.1e of the H2O production at the two times," % tuple(r)
    # cross-check of the H2O production shares against the total-rate shares of the v2 jobs
    xcheck = {}
    for t in times:
        d = json.load(open(Path(args.v2_dir) / "checks" / ("slope_bins_total_%s.json" % t)))["bins"]
        xcheck[t] = {"seen": float(np.max(np.abs(np.array([b["production_share_seen"] for b in d]) -
                                                 z["production_share_seen_%s" % t]))),
                     "all": float(np.max(np.abs(np.array([b["production_share_all"] for b in d]) -
                                                z["production_share_all_%s" % t])))}
    print("production shares, H2O (contributions) vs total (v2 fields): max abs difference %s" % json.dumps(xcheck))
    written = []

    # ---------------- per-bin figure
    fig, ax = plt.subplots(3, 1, figsize=(FIG.SINGLE_COLUMN_IN * 1.5, 7.0), sharex=True, gridspec_kw={"hspace": 0.10})
    w = 0.38
    ax[0].bar(x - w / 2, 100 * z["bin_area_fraction"], w, color="0.6", label="share of the total area")
    ax[0].bar(x + w / 2, 100 * z["bin_share_of_seen_area"], w, color="#0072B2", label="share of the seen area")
    ax[0].plot(x, 100 * z["seen_area_fraction_of_bin"], "o-", color="#E69F00", ms=3.5, lw=1.0,
               label="seen fraction of the bin's area")
    ax[0].set_ylabel("area [%]")
    ax[0].set_ylim(0, 100)
    ax[0].legend(loc="upper right", fontsize=6.3)
    bottom = 10.0 ** np.floor(np.log10(mult.min()) - 1)
    cols = ["0.6" if k == 0 else (FIG.CLIFF60_COLOR if k >= 3 else FIG.CLIFF45_COLOR) for k in range(nb)]
    ax[1].bar(x, mult - bottom, 0.6, bottom=bottom, color=cols)
    ax[1].set_yscale("log")
    ax[1].set_ylim(bottom, 10 ** np.ceil(np.log10(mult.max()) + 5))
    ax[1].axhline(1.0, color="k", lw=0.6, ls=(0, (3, 2)))
    for xi, mi in zip(x, mult):
        ax[1].text(xi, mi * 4, "%.3g" % mi + (" (held)" if xi == 0 else ""), ha="center", va="bottom", fontsize=5.6)
    ax[1].set_ylabel("multiplier $m_b$ of the activity")
    ww = 0.19
    tl = {"103310": "10:33:10", "103926": "10:39:26"}
    tc = {"103310": "#D55E00", "103926": "#CC79A7"}
    for k, t in enumerate(times):
        off = (k - 0.5) * 2 * ww
        ax[2].bar(x + off - ww / 2, 100 * z["production_share_seen_%s" % t], ww, color=tc[t],
                  label="seen surface, %s UTC" % tl[t])
        ax[2].bar(x + off + ww / 2, 100 * z["production_share_all_%s" % t], ww, color="white", edgecolor=tc[t],
                  hatch="//////", lw=0.6, label="whole surface, %s UTC" % tl[t])
    ax[2].set_ylabel("share of the production rate [%]")
    ax[2].set_ylim(0, 100)
    ax[2].legend(loc="upper right", fontsize=6.0, ncol=2)
    ax[2].set_xticks(x)
    ax[2].set_xticklabels(labs, fontsize=6.5)
    ax[2].set_xlabel("effective gravitational slope bin")
    for k, a in enumerate(ax):
        a.text(0.01, 0.95, "(%s)" % "abc"[k], transform=a.transAxes, va="top", fontsize=8)
        a.tick_params(axis="x", which="minor", bottom=False, top=False)
    cl = [k for k in range(nb) if k >= 3]
    cap = ("How the fitted state (%s) uses the effective-"
           "gravitational-slope bins of its multipliers. (a) Share of the total surface area (grey) and of the surface "
           "seen during the fit (blue; 'seen' = the facet contributes to at least one fitted NG or RG observation, the "
           "fit's nonneg_facets.npz flag) in each bin, and the fraction of each bin's area that is seen (orange). (b) The "
           "fitted multiplier m_b of the H2O activity in each bin (log scale; [0,40) deg held at 1): %s; light blue: "
           "40-60 deg, dark blue: >= 60 deg, grey: the held bin. (c) Share of the production rate (facet flux F times area) of the seen "
           "surface (solid) and of the whole surface (hatched) in each bin, at 10:33:10 and 10:39:26 UTC (H2O%s the "
           "total-rate shares differ from these by at most %.1e). The bins >= 60 deg cover %.1f %% of the "
           "seen area and give %.1f %% of its production at 10:33:10 UTC (%s).%s%s Data: %s/slope_bins.npz "
           "(job %s)." % (
               args.model_name, ", ".join("%s %.3g" % (l, m) for l, m in zip(labs, mult)),
               co2_txt or ";", max(max(v.values()) for v in xcheck.values()),
               100 * z["bin_share_of_seen_area"][cl].sum(), 100 * z["production_share_seen_103310"][cl].sum(),
               ", ".join("%s: %.1f %%" % (labs[k], 100 * z["production_share_seen_103310"][k]) for k in cl),
               off_txt, bound_txt, args.contrib_dir, cjob))
    pre = out / "12M_cliff_slope_bins"
    written += FIG.save(fig, pre)
    plt.close(fig)
    written.append(FIG.write_caption(pre, cap))

    # ---------------- descent figures
    tab = FIG.read_table(args.base)
    row = {int(p): i for i, p in enumerate(tab["iPoint"])}
    kind = np.asarray(z["kind"]).astype(str)
    keep = kind != "grid"
    ip = np.asarray(z["iPoint"])[keep]
    idx = np.array([row[int(p)] for p in ip])
    alt = tab["altitude"][idx]
    order = np.argsort(-alt)
    obs = (tab["refNG"] > 0) | (tab["refRG"] > 0)
    otab = {k: (v[obs] if isinstance(v, np.ndarray) and len(v) == len(obs) else v) for k, v in tab.items()}
    gaps = ng_gaps(otab, otab["altitude"])
    window = FIG.late_window(otab, args.bundle, None, None)
    frac = {g: np.asarray(z["fraction_%s" % g])[keep] for g in ("NG", "RG")}
    use = {"NG": np.asarray(z["use_ng"])[keep], "RG": np.asarray(z["use_rg"])[keep]}

    def panel(a, g):
        f = frac[g][order]
        ys = [f[:, list(ix)].sum(axis=1) for _, ix, _, _ in GROUPS]
        a.stackplot(alt[order], *ys, colors=[c for _, _, c, _ in GROUPS], lw=0, zorder=2)
        a.axvspan(window["altitude_min_m"], window["altitude_max_m"], facecolor="none", edgecolor="k", lw=0.6,
                  zorder=3)
        import matplotlib.transforms as mtransforms
        for hi, lo, *_ in gaps:
            if g == "NG":                 # no NG data: the NG panel keeps the gap hatched over
                a.axvspan(lo, hi, facecolor="white", edgecolor="0.3", hatch="////", lw=0, alpha=1.0, zorder=4)
            else:                         # RG was measured in the NG gap: model shares shown, the gap edges dashed
                for x in (lo, hi):        # (user decision 2026-11-12)
                    a.axvline(x, color="k", lw=0.8, ls=(0, (4, 2)), zorder=5)
                a.text(np.sqrt(hi * lo), 1.015, "gap in NG data", ha="center", va="bottom", fontsize=6.0,
                       transform=mtransforms.blended_transform_factory(a.transData, a.transAxes))
        a.set_xscale("log")
        a.set_xlim(alt.max() * 1.15, alt.min() * 0.85)
        a.set_ylim(0, 1)
        a.set_xlabel("gauge altitude [m] (descent →)")
        a.set_title("%s" % g, fontsize=8, loc="left", pad=3)

    from matplotlib.lines import Line2D
    h_slope = [Patch(facecolor=c, edgecolor="none", label="slope %s" % l) for l, _, c, _ in GROUPS]
    h_dash = Line2D([], [], color="k", lw=0.8, ls=(0, (4, 2)), label="edges of the gap in the NG data (RG measured)")
    h_late = Patch(facecolor="none", edgecolor="k", label="late window %s–%s UTC" % (
        window["start"].split("T")[1], window["stop"].split("T")[1]))
    nums = {g: {l: frac[g][:, list(ix)].sum(axis=1) for l, ix, _, _ in GROUPS} for g in ("NG", "RG")}
    for g in ("NG", "RG"):                    # per bin (independent of the grouping)
        for k in range(nb):
            nums[g]["bin%d" % k] = frac[g][:, k]
        nums[g][">=60"] = frac[g][:, 3:].sum(axis=1)
    late = np.isin(ip, window["points"])

    def stat(g, l, sel):
        v = nums[g][l][sel]
        return "%.0f-%.0f %%" % (100 * v.min(), 100 * v.max())
    for gauges in (("NG",), ("RG",), ("NG", "RG")):
        W = FIG.DOUBLE_COLUMN_IN if len(gauges) == 2 else FIG.DOUBLE_COLUMN_IN * 0.6
        fig, axs = plt.subplots(1, len(gauges), figsize=(W, 3.3 if len(gauges) == 2 else 3.6), sharey=True, squeeze=False,
                                gridspec_kw={"wspace": 0.08})
        for a, g in zip(axs[0], gauges):
            panel(a, g)
        axs[0][0].set_ylabel("fraction of the model surface pressure")
        fig.subplots_adjust(bottom=0.33, top=0.92, left=0.11 if len(gauges) == 2 else 0.16, right=0.98)
        h_hatch = Patch(facecolor="white", edgecolor="0.3", hatch="////",
                        label="gap in the NG data (RG measured)" if gauges == ("NG",) else "gap in the NG data (NG panel)")
        handles = h_slope + ([h_hatch] if "NG" in gauges else []) + ([h_dash] if "RG" in gauges else []) + [h_late]
        fig.legend(handles=handles, loc="lower left", bbox_to_anchor=(0.02, 0.0), ncol=3 if len(gauges) == 2 else 2,
                   fontsize=6.0, frameon=False, handlelength=1.6, columnspacing=1.0)
        cap = ("Fraction of the modelled surface pressure of %s contributed by each effective-gravitational-slope bin "
               "along the descent of 2016-09-30 (%s; stacked, %s; multipliers %s);%s the x axis "
               "is the altitude of the COPS gauges above the shape model in the fit's base (1.5 m toward the comet from "
               "the spacecraft reference point), logarithmic, descent from left to right; every observation row (%d; "
               "for the 95 observations from 10:33:10 UTC the response-averaged value of the fit). Box: the late window "
               "%s-%s UTC (gauge altitude %.0f-%.1f m). The gap in the NG data (%s), during which RG was measured: %s. "
               "The NG grain term (a separate model term) is not included. %s Data: %s/slope_bins.npz "
               "(job %s)." % (
                   " and ".join(gauges), args.model_name,
                   ("the bins grouped: " if any(len(ix) > 1 for _, ix, _, _ in GROUPS) else "one colour per bin: ") +
                   ", ".join("%s (%s)" % (l.replace("–", "-").replace("°", " deg").replace(" (overhangs)", " = overhangs"), cn)
                             for l, _, _, cn in GROUPS),
                   ", ".join("%s %.3g" % (l, m) for l, m in zip(labs, mult)),
                   (" bins %s switched off by their multipliers (<= %.2g), invisible;" % (
                       ", ".join(labs[k] for k in offb), mult[offb].max()) if offb else ""),
                   keep.sum(), window["start"].split("T")[1], window["stop"].split("T")[1],
                   window["altitude_max_m"], window["altitude_min_m"],
                   ", ".join("%.0f-%.0f m, %s-%s UTC" % (hi, lo, t0.strftime("%H:%M:%S"), t1.strftime("%H:%M:%S"))
                             for hi, lo, t0, t1 in gaps),
                   "; ".join({"NG": "hatched over in the NG panel (no NG data)",
                              "RG": "in the RG panel the model's RG shares are shown there (the RG observation rows of "
                                    "the gap) and the gap edges are dashed lines labelled 'gap in NG data'; the lower "
                                    "edge (10:33:10 UTC) is also the start of the late window, so it lies on the left "
                                    "edge of the late-window box"}[g]
                             for g in gauges),
                   " ".join("%s: slopes >= 60 deg give %s of the surface pressure before the late window and %s in it "
                            "(%s in it)." % (
                                g, stat(g, ">=60", ~late), stat(g, ">=60", late),
                                ", ".join("%s %s" % (labs[k], stat(g, "bin%d" % k, late)) for k in range(3, nb)))
                            for g in gauges),
                   args.contrib_dir, cjob))
        pre = out / ("12M_cliff_pressure_fraction_%s" % "_".join(gauges))
        written += FIG.save(fig, pre)
        plt.close(fig)
        written.append(FIG.write_caption(pre, cap))
    json.dump({"production_share_crosscheck_max_abs": xcheck, "files": written},
              open(out / "cliff_figures_manifest.json", "x"), indent=1)
    print("\n".join(written))
    return 0


if __name__ == "__main__":
    sys.exit(main())
