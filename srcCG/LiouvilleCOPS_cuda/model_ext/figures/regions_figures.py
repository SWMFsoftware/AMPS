#!/usr/bin/env python3
# ==== version-controlled copy: header added 2026-11-13, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/figures_cluster/scripts/regions_figures.py (on Voyager: /data001/highresolu_vt/cops_output/figures_cluster/scripts/regions_figures.py)
# SHA256 of the run copy = SHA256 of this file without these header lines: ef912e8d018216512b3b03b51b8196cfb536be389fd9e5a288897496efe03a58
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: Stacked share of the modelled surface pressure by lobe and by distance shell from the touchdown point along
#   the descent (cliff-figure style).
# RUN: srun on normal128cores (CPU)
# ==== end of header ====
"""Share of the modelled surface signal by region along the descent (2026-11-14, lead: "stacked share of the modelled
surface signal vs gauge altitude by (i) lobe (head / body, the neck plane of the A10 analysis) and (ii) distance shells
from the touchdown point; NG and RG; NG gap hatched / marked as in the cliff figures; late window boxed").

Data: contrib_phase_data.py regions.npz (per base row: the fitted model's NG and RG surface pressure summed by lobe and by
distance shell of the facet centroid from the touchdown point; the late observations = the W10x5 response average of the
1-s grid rows as in the fit). The style is that of cliff_importance_figures.py (stackplot over the observation rows,
log altitude axis descending, the NG gap hatched over in the NG panel, its edges dashed in the RG panel, the late window
boxed). Figures (.eps, .png 300 dpi, _caption.txt) in <out-dir>/:
  12M_region_share_lobe_{NG,RG,NG_RG}, 12M_region_share_distance_{NG,RG,NG_RG}
Run (CPU): regions_figures.py --data <contrib_phase dir> --base <fit base> --bundle <bundle .npz> --out-dir <dir>
"""

import argparse
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
APP = Path("/data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda")
sys.path.insert(0, str(APP))
from cops_cuda import figures as FIG  # noqa: E402
import cliff_importance_figures as CF  # noqa: E402

LOBE_COLS = ("#E69F00", "#0072B2")
SHELL_COLS = ("#D55E00", "#E69F00", "#F0E442", "#009E73", "#56B4E9", "#0072B2")


def main():
    ap = argparse.ArgumentParser()
    for a in ("--data", "--base", "--bundle", "--out-dir"):
        ap.add_argument(a, required=True)
    ap.add_argument("--model-name", default="the softplus activity model with slope multipliers >= 1, final state after 9 "
                    "fit and fine-tuning cycles")
    args = ap.parse_args()
    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=False)
    plt = FIG.style()
    from matplotlib.patches import Patch
    from matplotlib.lines import Line2D
    import matplotlib.transforms as mtransforms
    z = np.load(Path(args.data) / "regions.npz", allow_pickle=False)
    man = json.load(open(Path(args.data) / "manifest.json"))
    keep = np.asarray(z["obs"]) | np.asarray(z["late"])
    alt = np.asarray(z["alt"])[keep]
    order = np.argsort(-alt)
    tab = FIG.read_table(args.base)
    obsm = (tab["refNG"] > 0) | (tab["refRG"] > 0)
    otab = {k: (v[obsm] if isinstance(v, np.ndarray) and len(v) == len(obsm) else v) for k, v in tab.items()}
    gaps = CF.ng_gaps(otab, otab["altitude"])
    window = FIG.late_window(otab, args.bundle, None, None)
    refng = np.asarray(z["refng"])[keep]
    refrg = np.asarray(z["refrg"])[keep]
    late = np.asarray(z["late"])[keep]
    shells = np.asarray(z["shells_m"])
    shell_labs = ["%s–%s m" % (("%g" % shells[i]), ("%g" % shells[i + 1])) if np.isfinite(shells[i + 1]) else
                  "> %g m" % shells[i] for i in range(len(shells) - 1)]
    groups = {"lobe": (["small lobe (head, with the touchdown point)", "big lobe (body)"], LOBE_COLS, "lobe"),
              "distance": (["distance %s" % s for s in shell_labs], SHELL_COLS, "shell")}
    written, stats = [], {}
    for gname, (labels, cols, key) in groups.items():
        frac = {}
        for g in ("NG", "RG"):
            parts = np.asarray(z["%s_%s_Pa" % (key, g.lower())])[keep]
            tot = np.asarray(z["total_%s_Pa" % g.lower()])[keep]
            with np.errstate(invalid="ignore", divide="ignore"):
                frac[g] = parts / tot[:, None]
        sel = {"NG": refng > 0, "RG": refrg > 0}

        def panel(a, g):
            s = sel[g][order]
            x = alt[order][s]
            f = frac[g][order][s]
            a.stackplot(x, *[f[:, i] for i in range(f.shape[1])], colors=cols, lw=0, zorder=2)
            a.axvspan(window["altitude_min_m"], window["altitude_max_m"], facecolor="none", edgecolor="k", lw=0.6, zorder=3)
            for hi, lo, *_ in gaps:
                if g == "NG":
                    a.axvspan(lo, hi, facecolor="white", edgecolor="0.3", hatch="////", lw=0, alpha=1.0, zorder=4)
                else:
                    for xx in (lo, hi):
                        a.axvline(xx, color="k", lw=0.8, ls=(0, (4, 2)), zorder=5)
                    a.text(np.sqrt(hi * lo), 1.015, "gap in NG data", ha="center", va="bottom", fontsize=6.0,
                           transform=mtransforms.blended_transform_factory(a.transData, a.transAxes))
            a.set_xscale("log")
            a.set_xlim(alt.max() * 1.15, alt.min() * 0.85)
            a.set_ylim(0, 1)
            a.set_xlabel("gauge altitude [m] (descent →)")
            a.set_title(g, fontsize=8, loc="left", pad=3)
        st = {}
        for g in ("NG", "RG"):
            pre_ = sel[g] & ~late
            st[g] = {lab: {"before_late_window": [float(np.nanmin(frac[g][pre_, i])), float(np.nanmax(frac[g][pre_, i])),
                                                  float(np.nanmedian(frac[g][pre_, i]))],
                           "late_window": [float(np.nanmin(frac[g][sel[g] & late, i])), float(np.nanmax(frac[g][sel[g] & late, i])),
                                           float(np.nanmedian(frac[g][sel[g] & late, i]))]}
                     for i, lab in enumerate(labels)}
        stats[gname] = st
        h_reg = [Patch(facecolor=c, edgecolor="none", label=l) for l, c in zip(labels, cols)]
        h_dash = Line2D([], [], color="k", lw=0.8, ls=(0, (4, 2)), label="edges of the gap in the NG data (RG measured)")
        h_late = Patch(facecolor="none", edgecolor="k", label="late window %s–%s UTC" % (
            window["start"].split("T")[1], window["stop"].split("T")[1]))
        for gauges in (("NG",), ("RG",), ("NG", "RG")):
            W = FIG.DOUBLE_COLUMN_IN if len(gauges) == 2 else FIG.DOUBLE_COLUMN_IN * 0.6
            fig, axs = plt.subplots(1, len(gauges), figsize=(W, 3.4 if len(gauges) == 2 else 3.8), sharey=True,
                                    squeeze=False, gridspec_kw={"wspace": 0.08})
            for a, g in zip(axs[0], gauges):
                panel(a, g)
            axs[0][0].set_ylabel("fraction of the model surface pressure")
            fig.subplots_adjust(bottom=0.36, top=0.92, left=0.11 if len(gauges) == 2 else 0.16, right=0.98)
            h_hatch = Patch(facecolor="white", edgecolor="0.3", hatch="////",
                            label="gap in the NG data" + (" (RG measured)" if gauges == ("NG",) else " (NG panel)"))
            handles = h_reg + ([h_hatch] if "NG" in gauges else []) + ([h_dash] if "RG" in gauges else []) + [h_late]
            fig.legend(handles=handles, loc="lower left", bbox_to_anchor=(0.02, 0.0), ncol=3 if len(gauges) == 2 else 2,
                       fontsize=6.0, frameon=False, handlelength=1.6, columnspacing=1.0)
            what = ("by lobe: the small lobe (head, where the touchdown point lies) and the big lobe (body), separated by "
                    "the plane of the neck (the plane of minimum surface area per unit length between the two lobes' "
                    "centres; ng_source_regions.py, FINDINGS A10)" if gname == "lobe" else
                    "by the distance of the facet centroid from the touchdown point (Rosetta impact, lon 12.164°, lat "
                    "31.685°): %s" % ", ".join(shell_labs))
            cap = ("Fraction of the modelled surface pressure of %s contributed %s, along the descent of 2016-09-30 (%s; "
                   "stacked); the x axis is the altitude of the COPS gauges above the shape model in the fit's base, "
                   "logarithmic, descent from left to right; every observation row of the gauge (for the 95 observations "
                   "from 10:33:10 UTC the response-averaged value of the fit, W10x5 of the 1-s grid rows). Box: the late "
                   "window %s-%s UTC (gauge altitude %.0f-%.1f m). The gap in the NG data (%s): %s. Contributions: the "
                   "fit's kernels and arithmetic per facet (the NG grain term, a separate model term, is not included); "
                   "the sums of the regions equal the model surface pressure (max relative difference %.1e NG, %.1e RG) "
                   "and the model equals the contribution data's per-row values to %.1e. %s Data: %s/regions.npz." % (
                       " and ".join(gauges), what, args.model_name, window["start"].split("T")[1],
                       window["stop"].split("T")[1], window["altitude_max_m"], window["altitude_min_m"],
                       ", ".join("%.0f-%.0f m, %s-%s UTC" % (hi, lo, t0.strftime("%H:%M:%S"), t1.strftime("%H:%M:%S"))
                                 for hi, lo, t0, t1 in gaps),
                       "; ".join({"NG": "hatched over in the NG panel (no NG data)",
                                  "RG": "in the RG panel the model's RG shares are shown there and the gap edges are "
                                        "dashed"}[g] for g in gauges),
                       man["checks"]["ng"]["%s_sum_vs_total_max_rel" % key], man["checks"]["rg"]["%s_sum_vs_total_max_rel" % key],
                       max(man["checks"]["ng"]["max_rel_vs_contribution_data"], man["checks"]["rg"]["max_rel_vs_contribution_data"]),
                       " ".join("%s: %s." % (g, "; ".join("%s median %.0f %% before the late window, %.0f %% in it (%.0f-%.0f %%)"
                                                           % (lab, 100 * v["before_late_window"][2], 100 * v["late_window"][2],
                                                              100 * v["late_window"][0], 100 * v["late_window"][1])
                                                           for lab, v in stats[gname][g].items()))
                                for g in gauges),
                       args.data))
            pre = out / ("12M_region_share_%s_%s" % (gname, "_".join(gauges)))
            written += FIG.save(fig, pre)
            plt.close(fig)
            written.append(FIG.write_caption(pre, cap))
    json.dump({"stats": stats, "files": written}, open(out / "regions_figures_manifest.json", "x"), indent=1)
    print("\n".join(written))
    return 0


if __name__ == "__main__":
    sys.exit(main())
