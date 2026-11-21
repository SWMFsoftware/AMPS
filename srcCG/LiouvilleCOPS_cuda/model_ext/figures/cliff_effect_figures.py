#!/usr/bin/env python3
# ==== version-controlled copy: header added 2026-11-17, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/figures_cluster/scripts/cliff_effect_figures.py (on Voyager: /data001/highresolu_vt/cops_output/figures_cluster/scripts/cliff_effect_figures.py)
# SHA256 of the run copy = SHA256 of this file without these header lines: cc95d3372c26a2cafd6d52183be60e1cd597f80130a46fa77d57ab32d310ab30
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: figures P1-P4 and answer tables (SZA-controlled ratios) of the slope-at-equal-illumination analysis of one
#   state (2026-11-17)
# RUN: CPU Slurm job jobs/cliff_effect_figures.sbatch for all cases (HOW-TO-RERUN.md of
#   figures_cluster/12M/cliff_effect_stage2/)
# ==== end of header ====
"""Figures and answer tables of the slope-at-equal-illumination analysis of ONE state (2026-11-17; data:
cliff_effect_data.py). CPU, light.

    cliff_effect_figures.py --data <data/<key>> --out <cliff_effect_stage2/<key>> --label <case label> [--caveat T ...]

Slope categories (the fit's slope bins): 0-40 (reference), 40-50, 50-60, 60-70, 70-80, 80-90, thresholds >40, >60, >80
(everything above, incl. 90-180) and, as an addition (not in the user's list), 90-180. Illumination classes: lit SZA in 10
deg bins, shadowed (cos SZA > 0, cast shadow), night (cos SZA <= 0).
Figures (.eps + .png 300 dpi + _caption.txt):
  P1_<NG|RG>      source rate F and activity A at equal illumination on the gauge-seen surface, the four snapshot times
  P1_whole        the same for the whole surface at the map times 10:33:10 and 10:39:26
  P2_<NG|RG>      observation weight per m^2 g and its factors 1/r^2, cos theta, gauge weight W, visible fraction
  P3_<NG|RG>      matched cells (SZA class x 10 % distance x 0.1 cos theta): contribution-weighted distribution of the
                  ratio of the contribution per m^2, category vs 0-40
  P4_<NG|RG>      signal share / area share = activity x illumination x geometry x covariance (exact)
Tables: answer.csv / answer.json (SZA-controlled ratio category vs 0-40 for F and A: the area-weighted average over the
illumination classes of the per-class ratio, weights = the category's area per class; spread = min / max over the classes
holding >= 1 % of the category's area, and the weighted standard deviation of log10), p3.csv, p4.csv, p1_whole.csv.
"""
import argparse
import csv
import json
import sys
from pathlib import Path

import numpy as np

APP = Path("/data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda")
sys.path.insert(0, str(APP))
from cops_cuda import figures as FIG  # noqa: E402

BINLAB = ["0-40", "40-50", "50-60", "60-70", "70-80", "80-90", "90-180"]
CATS = [("0-40", [0]), ("40-50", [1]), ("50-60", [2]), ("60-70", [3]), ("70-80", [4]), ("80-90", [5]),
        (">40", [1, 2, 3, 4, 5, 6]), (">60", [3, 4, 5, 6]), (">80", [5, 6]), ("90-180 (addition)", [6])]
CLS = ["0-10", "10-20", "20-30", "30-40", "40-50", "50-60", "60-70", "70-80", "80-90", "shad.", "night"]
SNAPS = [(1079, "10:33:10 UTC (297 m)"), (1406, "10:37:31.282 UTC (100 m)"), (1528, "10:39:09.282 UTC (20 m)"),
         (1549, "10:39:26 UTC (5.5 m)")]
PHASES = [("pre_gap", "before the NG gap"), ("gap", "the NG gap (RG)"), ("late_W10x5", "late window (W10x5)")]
STY = {"0-40": dict(color="k", lw=1.6, ls="-"), "40-50": dict(color="#440154", lw=1.0, ls="-"),
       "50-60": dict(color="#3b528b", lw=1.0, ls="-"), "60-70": dict(color="#21918c", lw=1.0, ls="-"),
       "70-80": dict(color="#5ec962", lw=1.0, ls="-"), "80-90": dict(color="#c8b400", lw=1.0, ls="-"),
       ">40": dict(color="#d62728", lw=1.4, ls="--"), ">60": dict(color="#ff7f0e", lw=1.4, ls="--"),
       ">80": dict(color="#e377c2", lw=1.4, ls="--"), "90-180 (addition)": dict(color="0.5", lw=1.0, ls=":")}


def cat_sum(arr, bins):
    """Sum over the slope bins of a category (last axis = bins)."""
    return arr[..., bins].sum(axis=-1)


def safe_div(a, b):
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where(b > 0, a / np.where(b > 0, b, 1.0), np.nan)


def p1_tables(p1):
    """p1 (4, NJ, NB): area, A area, F area, active area -> per category arrays over the classes."""
    out = {}
    for name, bins in CATS:
        a, aA, aF, act = (cat_sum(p1[k], bins) for k in range(4))
        out[name] = {"area": a, "meanA": safe_div(aA, a), "meanF": safe_div(aF, a), "active": safe_div(act, a),
                     "meanF_active": safe_div(aF, act), "meanA_active": safe_div(aA, act)}
    return out


RLIM = (1e-3, 1e4)  # ratio axes of P1; points outside are drawn at the edge as triangles
VSPAN = 1e7         # other P1 log axes: at most 7 decades below the largest value (points below: triangles at the edge)
REFMIN = 0.01      # a class's 0-40 reference is "under-sampled" below 1 % of the reference's area on that surface


def controlled_ratio(t, q, refmin=None):
    """SZA-controlled ratio of quantity q ('meanF' / 'meanA') of every category vs 0-40 (the lead's definition:
    area-weighted average over the illumination classes of the per-class ratio of the means, weights = the category's
    area per class). refmin: sensitivity variant (an addition) comparing only the classes whose 0-40 reference area is >=
    refmin x the reference's area on that surface. weight_fraction_undersampled_ref: the fraction of the averaging weight
    in classes whose reference area is below REFMIN x its total (where the per-class ratio rests on few 0-40 facets)."""
    ref = t["0-40"][q]
    ra = t["0-40"]["area"]
    under = ra < REFMIN * ra.sum()
    res = {}
    for name, _ in CATS[1:]:
        a = t[name]["area"]
        r = safe_div(t[name][q], ref)
        ok = np.isfinite(r) & (a > 0) & (ref > 0)
        if refmin is not None:
            ok &= ra >= refmin * ra.sum()
        if not ok.any() or a.sum() == 0:
            res[name] = None
            continue
        w = a[ok]
        rr = r[ok]
        big = a[ok] >= 0.01 * a.sum()
        lg = np.log10(np.where(rr > 0, rr, np.nan))
        fin = np.isfinite(lg)
        mean = float(np.sum(w * rr) / np.sum(w))
        sd = float(np.sqrt(np.sum(w[fin] * (lg[fin] - np.sum(w[fin] * lg[fin]) / np.sum(w[fin])) ** 2) / np.sum(w[fin]))) \
            if fin.any() else float("nan")
        zero_ref = (a > 0) & ~(ref > 0)
        res[name] = {"ratio": mean, "min": float(np.min(rr[big])) if big.any() else float("nan"),
                     "max": float(np.max(rr[big])) if big.any() else float("nan"), "sd_log10": sd,
                     "classes": int(ok.sum()), "area_fraction_compared": float(w.sum() / a.sum()),
                     "area_fraction_ref_zero": float(a[zero_ref].sum() / a.sum()),
                     "weight_fraction_undersampled_ref": float(a[ok & under].sum() / w.sum()),
                     "area_weighted_mean_ratio_of_sums": float(np.sum(a[ok] * t[name][q][ok]) / np.sum(a[ok] * ref[ok]))}
    return res


def ratios(t):
    return {"F": controlled_ratio(t, "meanF"), "A": controlled_ratio(t, "meanA"),
            "F_refmin": controlled_ratio(t, "meanF", REFMIN), "A_refmin": controlled_ratio(t, "meanA", REFMIN)}


def style():
    return FIG.style()


def plot_p1(tabs, titles, path, label, gauge_text, caption_extra):
    plt = style()
    ncol = len(tabs)
    rows = [("meanF", "mean F [m$^{-2}$ s$^{-1}$]", True), ("ratioF", "F / F(0-40)", True), ("meanA", "mean A", True),
            ("ratioA", "A / A(0-40)", True), ("active", "active-area fraction (A > 0)", True),
            ("meanF_active", "mean F of the active facets", True)]
    fig, ax = plt.subplots(len(rows) + 1, ncol, figsize=(3.6 * ncol + 1.6, 2.15 * (len(rows) + 1)), squeeze=False,
                           gridspec_kw={"hspace": 0.32, "wspace": 0.28})
    x = np.arange(len(CLS))
    for c, (t, ttl) in enumerate(zip(tabs, titles)):
        for r_, (key, lab, logy) in enumerate(rows):
            a = ax[r_, c]
            ys = {}
            for name, _ in CATS:
                if key == "ratioF":
                    y = safe_div(t[name]["meanF"], t["0-40"]["meanF"])
                elif key == "ratioA":
                    y = safe_div(t[name]["meanA"], t["0-40"]["meanA"])
                else:
                    y = t[name][key]
                ys[name] = np.where(np.isfinite(y) & (y > 0), y, np.nan)
            if key.startswith("ratio"):
                lim = RLIM
            else:                                   # log axis over at most VSPAN decades below the largest value
                allv = np.concatenate([v[np.isfinite(v)] for v in ys.values()] + [np.array([1.0])])
                top = float(np.max(allv)) * 2.0
                lim = (max(float(np.min(allv)) / 2.0, top / VSPAN), top)
            for name, _ in CATS:
                y = ys[name]
                out_hi, out_lo = y > lim[1], y < lim[0]
                yc = np.clip(y, lim[0], lim[1])
                a.plot(x, np.where(out_hi | out_lo, np.nan, yc), marker="o", ms=2.2, label=name, **STY[name])
                a.plot(x[out_hi], yc[out_hi], marker="^", ms=3.5, ls="none", color=STY[name]["color"])
                a.plot(x[out_lo], yc[out_lo], marker="v", ms=3.5, ls="none", color=STY[name]["color"])
            ra = t["0-40"]["area"]
            for j in np.nonzero(ra < REFMIN * ra.sum())[0]:
                a.axvspan(j - 0.5, j + 0.5, color="0.88", lw=0, zorder=0)
            if logy:
                a.set_yscale("log")
            if key.startswith("ratio"):
                a.axhline(1.0, color="0.3", lw=0.6, ls=":")
            a.set_ylim(*lim)
            a.set_xticks(x)
            a.set_xticklabels(CLS if r_ == len(rows) - 1 else [], rotation=60, fontsize=6)
            a.tick_params(labelsize=6)
            if c == 0:
                a.set_ylabel(lab, fontsize=7)
            if r_ == 0:
                a.set_title(ttl, fontsize=8)
        # strip: illumination area fractions per category
        a = ax[len(rows), c]
        names = [n for n, _ in CATS]
        lit = [float(t[n]["area"][:9].sum()) for n in names]
        sh = [float(t[n]["area"][9]) for n in names]
        ni = [float(t[n]["area"][10]) for n in names]
        tot = np.array(lit) + np.array(sh) + np.array(ni)
        tot = np.where(tot > 0, tot, 1.0)
        xs = np.arange(len(names))
        a.bar(xs, np.array(lit) / tot, color="#f0c419", label="lit")
        a.bar(xs, np.array(sh) / tot, bottom=np.array(lit) / tot, color="#7f7f7f", label="shadowed")
        a.bar(xs, np.array(ni) / tot, bottom=(np.array(lit) + np.array(sh)) / tot, color="#1f2a44", label="night")
        a.set_xticks(xs)
        a.set_xticklabels(names, rotation=60, fontsize=5.5)
        a.set_ylim(0, 1)
        a.tick_params(labelsize=6)
        if c == 0:
            a.set_ylabel("area fraction", fontsize=7)
            a.legend(fontsize=5.5, loc="lower left", ncol=3, frameon=True)
    h, l = ax[0, 0].get_legend_handles_labels()
    fig.subplots_adjust(top=0.93)
    fig.legend(h, l, loc="upper center", ncol=10, fontsize=7, bbox_to_anchor=(0.5, 0.975))
    fig.suptitle("%s: source rate and activity at equal illumination by slope category, %s" % (label, gauge_text),
                 fontsize=9, y=0.99)
    files = FIG.save(fig, path)
    plt.close(fig)
    return files


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--label", required=True)
    ap.add_argument("--caveat", action="append", default=[])
    a = ap.parse_args()
    D, O = Path(a.data), Path(a.out)
    if O.exists():
        raise SystemExit("%s exists" % O)
    O.mkdir(parents=True)
    z = np.load(D / "cliff_data.npz", allow_pickle=False)
    meta = json.load(open(D / "cliff_data.json"))
    caps, answer, written = {}, {}, []
    common = ("Slope categories: the fit's effective-gravitational-slope bins 0-40 deg (reference), 40-50, 50-60, 60-70, "
              "70-80, 80-90, the thresholds >40, >60, >80 (everything above, including 90-180) and, as an addition not in "
              "the requested list, the fit's 90-180 bin (dotted grey). Illumination classes: lit facets in 10-deg bins of the "
              "solar zenith angle SZA (cos SZA > 0, Sun ray not blocked), shadowed (cos SZA > 0, ray-traced cast shadow: "
              "night-side source), night (cos SZA <= 0). Model: F = A q, q = Qmin + (Qmax - Qmin) cos SZA on the lit side, "
              "Qmin at night and in cast shadow (LiouvilleCOPS Source.cpp; base SZAMode continuous, CastShadows on), so at "
              "equal SZA and shadow state the slope dependence of F is that of the fitted activity A. Area-weighted means; "
              "about 95 % of the facets have A = 0 in the softplus and D2 states, so the active-area fraction and the mean "
              "over the active facets are shown separately. ")
    for g, gt in (("NG", "NG-seen surface (Kn > 0)"), ("RG", "RG-seen surface (KV > 0)")):
        tabs = [p1_tables(z["snap%d_%s_p1" % (r, g)]) for r, _ in SNAPS]
        files = plot_p1(tabs, [t for _, t in SNAPS], O / ("P1_%s" % g), a.label, gt, "")
        written += files
        answer["%s_seen" % g] = {}
        for (r, t), tb in zip(SNAPS, tabs):
            answer["%s_seen" % g][t] = ratios(tb)
        caps["P1_%s" % g] = (
            "Source rate F (instantaneous H2O emitted flux per m^2, the model's) and fitted activity A of %s at equal "
            "illumination, by slope category, on the surface seen by the %s gauge (records with %s > 0 in the fitting "
            "base's kernels) at the four contribution-density snapshot times (columns). Rows: area-weighted mean F per "
            "illumination class; its ratio to 0-40 deg at the same class; mean A; its ratio; the active-area fraction (A > "
            "0); the mean F of the active facets; the strip: each category's lit / shadowed / night area fractions. F, SZA "
            "and the shadow state per facet from the kernel records with the fit's arithmetic (nonneg_contrib.record_"
            "contributions, FP32 records); at 10:33:10 and 10:39:26 they equal the maps' F_H2O_per_m2_s (the production-rate "
            "maps' code path) to %.1e (max relative; zero pattern and illumination classes identical). Grey columns: "
            "classes where the 0-40 deg reference holds less than 1 %% of its seen area (under-sampled reference; there the "
            "per-class ratio rests on a few 0-40 deg facets and can reach any size); ratios beyond 1e-3 ... 1e4, and other values more than 7 decades below a panel's largest value, are drawn at the "
            "axis edge as triangles. %s%s" % (
                a.label, g, "Kn" if g == "NG" else "KV",
                max(meta["checks"]["rows"][str(r)].get("F_records_vs_maps_max_rel", 0) for r in (1079, 1549)), common,
                " ".join("Caveat: %s." % c for c in a.caveat)))
    # whole surface
    wt = [p1_tables(z["whole_%s_all_p1" % lab.replace(":", "")]) for lab in ("10:33:10", "10:39:26")]
    written += plot_p1(wt, ["10:33:10 UTC, whole surface", "10:39:26 UTC, whole surface"], O / "P1_whole", a.label,
                       "whole surface (map times)", "")
    answer["whole"] = {lab: ratios(t) for lab, t in zip(("10:33:10", "10:39:26"), wt)}
    fs = p1_tables(z["whole_103310_fit_seen_p1"])
    answer["fit_seen"] = {"10:33:10": ratios(fs)}
    caps["P1_whole"] = ("As P1 for the whole surface of %s (all 12,364,290 facets) at the map times 10:33:10 and 10:39:26 "
                        "UTC: F, cos SZA and the illumination class (day / night / cast shadow) per facet from the maps' "
                        "flux .vtu (the code path of the production-rate maps). %s" % (a.label, common))
    with open(O / "P1_whole_tables.csv", "x", newline="") as f:
        w = csv.writer(f)
        w.writerow(["time", "category", "class", "area_m2", "meanF", "meanA", "active_fraction", "meanF_active"])
        for lab, t in zip(("10:33:10", "10:39:26"), wt):
            for name, _ in CATS:
                for j, cl in enumerate(CLS):
                    w.writerow([lab, name, cl, "%.6e" % t[name]["area"][j]] + ["%.6e" % t[name][k][j] for k in
                                                                              ("meanF", "meanA", "active", "meanF_active")])
    # P2
    plt = style()
    for g, kname in (("NG", "Kn"), ("RG", "KV")):
        fig, ax = plt.subplots(5, 4, figsize=(16, 13), squeeze=False, gridspec_kw={"hspace": 0.3, "wspace": 0.28})
        labs = [("g", 5, "g = %s/dS [m$^{-2}$]" % kname), ("1/r$^2$", 1, "1/r$^2$ [m$^{-2}$]"), ("cos", 2, "cos θ"),
                ("W", 3, "gauge weight W"), ("fvis", 4, "visible fraction f")]
        for c, (r, ttl) in enumerate(SNAPS):
            p2 = z["snap%d_%s_p2" % (r, g)]
            for rr, (k, idx, lab) in enumerate(labs):
                axx = ax[rr, c]
                for name, bins in CATS:
                    y = safe_div(cat_sum(p2[idx], bins), cat_sum(p2[0], bins))
                    axx.plot(np.arange(len(CLS)), np.where(np.isfinite(y) & (y > 0), y, np.nan), marker="o", ms=2.2,
                             label=name, **STY[name])
                if k in ("g", "1/r$^2$"):
                    axx.set_yscale("log")
                axx.set_xticks(np.arange(len(CLS)))
                axx.set_xticklabels(CLS if rr == len(labs) - 1 else [], rotation=60, fontsize=6)
                axx.tick_params(labelsize=6)
                if c == 0:
                    axx.set_ylabel(lab, fontsize=7)
                if rr == 0:
                    axx.set_title(ttl, fontsize=8)
        h, l = ax[0, 0].get_legend_handles_labels()
        fig.subplots_adjust(top=0.92)
        fig.legend(h, l, loc="upper center", ncol=10, fontsize=7, bbox_to_anchor=(0.5, 0.97))
        fig.suptitle("%s: %s observation weight per m$^2$ by slope category at equal illumination" % (a.label, g),
                     fontsize=9, y=0.99)
        written += FIG.save(fig, O / ("P2_%s" % g))
        plt.close(fig)
        caps["P2_%s" % g] = (
            "%s observation weight per m^2 of surface, g = %s/dS (%s), and its factors at the facet centroid, by slope "
            "category and illumination class, on the %s-seen surface of %s at the four snapshot times: 1/r^2 (r from the "
            "fitting base's gauge position, GaugeOffsetSC 1.5 m), cos theta (emission angle), the gauge weight W (%s; "
            "boresights of the row from SPICE, the movie preparation's prep.npz) and the visible fraction f = g / (cos "
            "theta W / (pi r^2)) (the remainder: partial visibility and the kernel's sub-element sampling; can exceed 1). "
            "Area-weighted means; records whose centroid has cos theta <= 0 or W = 0 are excluded from the split (counts in "
            "p2_summary.json). Kn = sum (dS/r^2)(cos theta/pi) W_NG over the visible sub-elements (LiouvilleCOPS Liouville.h). %s" % (
                g, kname, "NG density kernel" if g == "NG" else "RG flux kernel without the (-l.b) factor", g, a.label,
                "NG: srccg_shield acceptance with the head-on weight 1 - l.b_NG, as the base" if g == "NG" else
                "RG: hemisphere acceptance, 1 for l.b_RG < 0", common))
    # P2 summary: ratio of the area-weighted means over all illumination classes, category / 0-40
    p2s = {}
    for g in ("NG", "RG"):
        p2s[g] = {}
        for ekey, etl in [("snap%d" % r, t) for r, t in SNAPS] + [("phase_%s" % p, t) for p, t in PHASES]:
            p2 = z["%s_%s_p2" % (ekey, g)]
            nn = z["%s_%s_n" % (ekey, g)]
            ref = {f: p2[i][:, 0].sum() / p2[0][:, 0].sum() if p2[0][:, 0].sum() > 0 else float("nan")
                   for f, i in (("inv_r2", 1), ("cos", 2), ("W", 3), ("fvis", 4), ("g", 5))}
            d = {"records_seen": float(nn[0]), "excluded_cos_le_0": float(nn[1]), "excluded_W_0": float(nn[2]),
                 "records_total": float(nn[3])}
            for name, bins in CATS:
                a0 = p2[0][:, bins].sum()
                d[name] = {f: float(p2[i][:, bins].sum() / a0 / ref[f]) if a0 > 0 and ref[f] > 0 else float("nan")
                           for f, i in (("inv_r2", 1), ("cos", 2), ("W", 3), ("fvis", 4), ("g", 5))}
                d[name]["area_m2"] = float(a0)
            p2s[g][etl] = d
    json.dump(p2s, open(O / "p2_summary.json", "x"), indent=1)
    # P3: matched cells
    p3rows = []
    defs = [(">40", [1, 2, 3, 4, 5, 6]), (">60", [3, 4, 5, 6]), (">80", [5, 6]), ("40-50", [1]), ("50-60", [2]),
            ("60-70", [3]), ("70-80", [4]), ("80-90", [5]), ("90-180 (addition)", [6])]
    entries = [("snap%d" % r, t) for r, t in SNAPS] + [("phase_%s" % p, t) for p, t in PHASES]
    edges = np.linspace(-6, 6, 97)
    for g in ("NG", "RG"):
        fig, ax = plt.subplots(3, 3, figsize=(13, 10), squeeze=False, gridspec_kw={"hspace": 0.35, "wspace": 0.25})
        cols_ = ["#1b9e77", "#d95f02", "#7570b3", "#e7298a", "0.2", "0.5", "#66a61e"]
        for k_, (dname, bins) in enumerate(defs):
            axx = ax[k_ // 3, k_ % 3]
            for e_, (ekey, etl) in enumerate(entries):
                if g == "NG" and ekey == "phase_gap":
                    continue
                p3 = z["%s_%s_p3" % (ekey, g)]
                ak, ck = cat_sum(p3[0].T, bins), cat_sum(p3[1].T, bins)
                ar, cr = p3[0][0], p3[1][0]
                both = (ak > 0) & (ar > 0)
                wk = ck
                with np.errstate(divide="ignore", invalid="ignore"):
                    R = (ck / np.where(ak > 0, ak, 1)) / (cr / np.where(ar > 0, ar, 1))
                tot = wk.sum()
                ok = both & (cr > 0) & (ck > 0)
                frac_matched = float(wk[ok].sum() / tot) if tot > 0 else float("nan")
                frac_noref = float(wk[(ck > 0) & ~((ar > 0) & (cr > 0))].sum() / tot) if tot > 0 else float("nan")
                lr = np.log10(R[ok])
                h, _ = np.histogram(np.clip(lr, edges[0], edges[-1]), bins=edges, weights=wk[ok])
                if tot > 0:
                    axx.step(0.5 * (edges[1:] + edges[:-1]), h / tot, where="mid", color=cols_[e_], lw=1.0,
                             label=etl)
                o = np.argsort(lr)
                cw = np.cumsum(wk[ok][o])
                med = float(lr[o][np.searchsorted(cw, 0.5 * cw[-1])]) if len(cw) else float("nan")
                p3rows.append({"gauge": g, "definition": dname, "entry": etl, "weighted_median_log10_ratio": med,
                               "weighted_median_ratio": 10 ** med if np.isfinite(med) else float("nan"),
                               "contribution_fraction_matched": frac_matched,
                               "contribution_fraction_without_flat_reference": frac_noref,
                               "matched_cells": int(ok.sum())})
            axx.axvline(0.0, color="0.3", lw=0.6, ls=":")
            axx.set_title("%s vs 0-40" % dname, fontsize=8)
            axx.tick_params(labelsize=6)
            if k_ // 3 == 2:
                axx.set_xlabel("log$_{10}$ (contribution per m$^2$, category / 0-40) in a matched cell", fontsize=7)
            if k_ % 3 == 0:
                axx.set_ylabel("fraction of the category's contribution", fontsize=7)
        h, l = ax[0, 0].get_legend_handles_labels()
        fig.subplots_adjust(top=0.89)
        fig.legend(h, l, loc="upper center", ncol=4, fontsize=7, bbox_to_anchor=(0.5, 0.965))
        fig.suptitle("%s: %s matched comparison (same SZA class, 10 %% distance, 0.1 cos θ)" % (a.label, g), fontsize=9,
                     y=0.99)
        written += FIG.save(fig, O / ("P3_%s" % g))
        plt.close(fig)
        caps["P3_%s" % g] = (
            "Matched comparison for %s, %s: the base rows' gauge-seen facets are grouped into cells of equal illumination "
            "class (10-deg SZA bins, shadowed, night), distance from the gauge (bins of a factor 1.1) and emission angle (0.1 "
            "in cos theta); in every cell holding both the category and 0-40 deg facets, the ratio of the area-weighted "
            "mean %s contribution per m^2 (category / 0-40) is formed; histograms: the distribution of log10 of that ratio "
            "weighted by the category's contribution in the cell, for the four snapshot rows and the phases pooled over "
            "their rows (before the NG gap; the gap, RG only; the late window = the W10x5 average of the 1-s grid rows). "
            "Contribution with no 0-40 deg facet (or a zero 0-40 deg contribution) in its cell is not in the histograms "
            "(fraction in p3.csv); ratios beyond 1e-6 ... 1e6 are counted in the edge bins (a cell whose 0-40 deg facets have a "
            "near-zero activity gives an arbitrarily large ratio). At matched illumination, distance and emission angle the ratio reflects the activity "
            "and the remaining geometry (gauge weight, visibility). %s" % (a.label, g, g, common))
    with open(O / "p3.csv", "x", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(p3rows[0].keys()))
        w.writeheader()
        for r in p3rows:
            w.writerow(r)
    # P4: exact decomposition
    p4rows = []
    for g in ("NG", "RG"):
        ents = [e for e in entries if not (g == "NG" and e[0] == "phase_gap")]
        fig, ax = plt.subplots(2, 4, figsize=(17, 8.5), squeeze=False, gridspec_kw={"hspace": 0.45, "wspace": 0.25})
        fac = [("share/area", "#000000"), ("activity", "#d62728"), ("illumination", "#f0c419"), ("geometry", "#1f77b4"),
               ("covariance", "#7f7f7f")]
        for e_, (ekey, etl) in enumerate(ents):
            p4 = z["%s_%s_p4" % (ekey, g)]
            tot = p4.sum(axis=1)
            axx = ax[e_ // 4, e_ % 4]
            names = [n for n, _ in CATS]
            vals = np.full((len(names), 5), np.nan)
            for i, (name, bins) in enumerate(CATS):
                s = cat_sum(p4, bins)
                if s[1] <= 0 or tot[0] <= 0:
                    continue
                share, ashare = s[0] / tot[0], s[1] / tot[1]
                mA, mI, mG = s[2] / s[1], s[3] / s[1], s[4] / s[1]
                MA, MI, MG = tot[2] / tot[1], tot[3] / tot[1], tot[4] / tot[1]
                cov_k = (s[0] / s[1]) / (mA * mI * mG) if mA * mI * mG > 0 else np.nan
                cov_t = (tot[0] / tot[1]) / (MA * MI * MG)
                parts = [share / ashare if ashare > 0 else np.nan, mA / MA, mI / MI, mG / MG,
                         cov_k / cov_t if np.isfinite(cov_k) else np.nan]
                vals[i] = parts
                p4rows.append({"gauge": g, "entry": etl, "category": name, "signal_share": share, "area_share": ashare,
                               "share_over_area": parts[0], "activity": parts[1], "illumination": parts[2],
                               "geometry": parts[3], "covariance": parts[4],
                               "product_check": parts[1] * parts[2] * parts[3] * parts[4] / parts[0]
                               if all(np.isfinite(parts)) and parts[0] > 0 else float("nan")})
            xs = np.arange(len(names))
            for k, (fn, col) in enumerate(fac):
                v = np.log10(np.where(vals[:, k] > 0, vals[:, k], np.nan))
                axx.bar(xs + (k - 2) * 0.16, np.where(np.isfinite(v), v, 0), width=0.16, color=col, label=fn)
            axx.axhline(0, color="0.3", lw=0.6)
            axx.set_xticks(xs)
            axx.set_xticklabels(names, rotation=60, fontsize=5.5)
            axx.set_title(etl, fontsize=8)
            axx.tick_params(labelsize=6)
            if e_ % 4 == 0:
                axx.set_ylabel("log$_{10}$ factor", fontsize=7)
        for e_ in range(len(ents), 8):
            ax[e_ // 4, e_ % 4].axis("off")
        h, l = ax[0, 0].get_legend_handles_labels()
        fig.subplots_adjust(top=0.87)
        fig.legend(h, l, loc="upper center", ncol=5, fontsize=7, bbox_to_anchor=(0.5, 0.955))
        fig.suptitle("%s: %s signal share / area share = activity x illumination x geometry x covariance" % (a.label, g),
                     fontsize=9, y=0.99)
        written += FIG.save(fig, O / ("P4_%s" % g))
        plt.close(fig)
        caps["P4_%s" % g] = (
            "Exact decomposition for %s, %s, per slope category on the gauge-seen surface (snapshot rows and the phases "
            "pooled over their rows): the category's share of the %s H2O signal divided by its share of the seen area "
            "(black) = activity factor <A>_k / <A> (red) x illumination factor <I>_k / <I> (yellow; I = q u for NG with u = "
            "1/sqrt(T), I = q for RG) x geometry factor <G>_k / <G> (blue; G = Kn/dS for NG, (KPhi + bm (V.b) u KV)/dS for "
            "RG) x covariance factor (grey; the rest: <A I G> / (<A><I><G>) of the category relative to the whole seen "
            "surface); area-weighted means < >; log10 bars (0 = no difference). The product equals the share ratio "
            "(p4.csv product_check). %s" % (a.label, g, g, common))
    with open(O / "p4.csv", "x", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(p4rows[0].keys()))
        w.writeheader()
        for r in p4rows:
            w.writerow(r)
    # answer table
    rows_ = []
    for surf, d in answer.items():
        for tl, v in d.items():
            for q in ("F", "A", "F_refmin", "A_refmin"):
                for name, res in v[q].items():
                    if res is None:
                        continue
                    rows_.append(dict(surface=surf, time=tl, quantity=q, category=name, **res))
    with open(O / "answer.csv", "x", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows_[0].keys()))
        w.writeheader()
        for r in rows_:
            w.writerow(r)
    json.dump({"answer": answer, "checks": meta["checks"], "label": a.label, "caveats": a.caveat},
              open(O / "answer.json", "x"), indent=1, default=str)
    for k, c in caps.items():
        (O / ("%s_caption.txt" % k)).open("x").write(c + "\n")
    print("written %d files to %s" % (len(written), O))
    return 0


if __name__ == "__main__":
    sys.exit(main())
