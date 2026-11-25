#!/usr/bin/env python3
# ==== version-controlled copy: header added 2026-11-21, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/figures_cluster/scripts/q002_slope_activity.py (on Voyager: /data001/highresolu_vt/cops_output/figures_cluster/scripts/q002_slope_activity.py)
# SHA256 of the run copy = SHA256 of this file without these header lines: 7dcb4cff0d165f618d21e701125be465b3c1065bfa674efa10f0e0f1d3aaffc4
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: touchdown source-rate maps; activity statistics per slope and visibility class (2026-11-21)
# RUN: see OUTPUT_VISUALIZATION.md 3.21
# ==== end of header ====
"""Surface activity against the angle to the local gravity, with the illumination accounted for, and maps of the source
rate around the touchdown point, for the five selected new-law fits (user 2026-11-21: "1) add a plot of the surface
source rate in the vicinity of the touchdawn: one panel is the source rate, another panel angle to the gravity 0-40,
40-50, 50-60, etc. 2) get statistics of the acrivity of the surface depending on the angle to the local gravity; accout
for the affect of changing of the solar illumination conditions (devide by cos of illuminated angle before doing the
snalysis -- propose a better metric ...). do this analyziz separately for faces seen at least by one gauge and then for
faces seen by both gauges concurrently or in different time").

The model's source rate of facet e at time t is F_e(t) = A_e q_e(t), with the illumination factor
    q = q_min (night or cast shadow), q = q_min + (q_max - q_min) cos SZA (lit),  q_min = 0.02 q_max.
Metrics per slope bin (area weighted, within each visibility class):
  M1 F/cos   the user's: F / cos SZA at 10:33:10 on lit facets with cos SZA >= 0.1 (undefined in night / shadow;
             = A [(q_max - q_min) + q_min / cos]: equals A q_max only for cos -> 1, grows at grazing incidence).
  M2 A       the activity A = F / q: the exact, time-independent illumination correction of the model's own law, defined
             on every facet (night and shadow included).
  M3 illum   the illumination each slope class gets: diurnal mean q / q_max over the 13 map times, lit fraction at
             10:33:10 (what any illumination correction has to remove).
  M4 local   A / <A>_patch, <A>_patch = the area-weighted mean A of the same visibility class in the facet's 150-m patch
             (Voronoi region of the option-C mesh): removes the large-scale (regional) variation, so it measures the
             contrast between slopes inside the same neighbourhood; also for the field without the slope multiplier,
             (A / m_b) / <A / m_b>_patch.
  M5 strata  mean A per slope bin inside illumination strata at 10:33:10 (dark; cos 0-0.3; 0.3-0.7; 0.7-1), relative
             to the stratum mean: does the slope dependence persist between facets with the same illumination?
Visibility classes (kernel > 0): any = seen by at least one gauge (the fit's seen mask); NG, RG; both = seen by NG and
by RG at any time (concurrently or at different times); NG only; RG only. From option A's per-phase contributions
(contrib_phase_c5/phase_means.npz: A > 0 on every facet, so a non-zero contribution means a non-zero kernel). With
--seen-gauges (q002_seen_gauges.py) also: both concurrently (same observation) and both only at different times.

The main classes exclude the held-template and late fine-tuned facets (TS != 0 or held != 0 in option A's
nonneg_facets.npz: 275 facets, 367 m2, which hold 11-39 % of the seen activity and would set the overhang means);
the classes *_all include them. "Off" = A below 1e-6 / 1e-2 of A_u, the uniform activity that gives the fit's own seen
production at 10:33:10 (bounded solvers leave switched-off coefficients slightly positive, so A == 0 does not count them).

    q002_slope_activity.py --out <dir> [--seen-gauges seen_gauges.npz] [--maps-only|--stats-only]
"""
import argparse
import csv
import json
import os
import struct
import sys
import time
from pathlib import Path

import numpy as np

APP = Path("/data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda")
sys.path.insert(0, str(APP))
import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.collections import PolyCollection  # noqa: E402
from matplotlib.colors import ListedColormap, LogNorm  # noqa: E402
from cops_cuda import figures  # noqa: E402

ROOT = Path("/data001/highresolu_vt/cops_output")
R12 = ROOT / "runs/shap8_012m_lategrid_dz1.5"
BASE = R12 / "base_q002"
GEO = R12 / "maps_geometry/job_136691/geometry"
GRID = R12 / "nonneg/localnn/grid_150m/grid.npz"
PL = ROOT / "figures_cluster/12M/q002_plausibility/metrics.json"
RUNS = [("C_f10_pm0p01_q002_reg", "C f10, P_m 0.01, templates on", "C_f10_pm0p01"),
        ("C_f1_pm0p01_q002_reg", "C f1, P_m 0.01, templates on", "C_f1_pm0p01"),
        ("C_f10_pm0p01_tploff_q002_reg", "C f10, P_m 0.01, templates off", "C_f10_pm0p01_tploff"),
        ("A_q002_reg", "A sum of squares, rank 3", "A"),
        ("B_q002_reg", "B background, bg 1e-2", "B")]
A_RUN = "A_q002_reg"
BIN_LAB = ["0-40°", "40-50°", "50-60°", "60-70°", "70-80°", "80-90°", "90-180°"]
BIN_COL = ["#9e9e9e", "#E69F00", "#D55E00", "#0b4f8a", "#56B4E9", "#009E73", "#7b3294"]
P = 4.2143229
TD_LON, TD_LAT = 12.164, 31.685
INK, GRID_C = "#0b0b0b", "#e4e3df"
T0 = time.time()
LAW = "new H2O law q_min = 0.02 q_max"


def say(t):
    print("q002_slope_activity [%5.0f s]: %s" % (time.time() - T0, t), flush=True)


def elements():
    path = BASE / "elements.bin"
    with open(path, "rb") as f:
        v, n, nb = struct.unpack("<iii", f.read(12))
    if v != 1 or nb != 25:
        raise SystemExit("elements.bin format")
    mm = np.memmap(path, dtype="<f8", mode="r", offset=12, shape=(n, 34))
    cen = np.array(mm[:, 0:3])
    area = np.array(mm[:, 6])
    del mm
    return n, cen, area


def map_times():
    T = {}
    with open(GEO / "times.txt") as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            t = line.split()
            if t[0].startswith("t"):
                T[t[0]] = {"utc": t[1], "rh_m": float(t[6])}
    return T


def qfac(label, rh_m, nf):
    """q per facet at one map time, q_max, cos SZA and shadow."""
    cos = np.fromfile(GEO / ("%s.cos.f8" % label), dtype="<f8")
    sh = np.fromfile(GEO / ("%s.shadow.u1" % label), dtype=np.uint8) != 0
    if len(cos) != nf:
        raise SystemExit("geometry size")
    rh = rh_m / 1.495978707e11
    qmax = 600.0 * rh ** (-P)
    qmin = 0.02 * qmax
    lit = (cos > 0) & ~sh
    q = np.full(nf, qmin)
    q[lit] = qmin + (qmax - qmin) * cos[lit]
    return q, qmax, cos, sh


def wquant(v, w, ps):
    o = np.argsort(v)
    c = np.cumsum(w[o])
    c /= c[-1]
    return [float(v[o][min(np.searchsorted(c, p), len(v) - 1)]) for p in ps]


def patch_mean(x, a, lab, npatch):
    num = np.bincount(lab, weights=x * a, minlength=npatch)
    den = np.bincount(lab, weights=a, minlength=npatch)
    with np.errstate(invalid="ignore", divide="ignore"):
        m = num / den
    return m


def nanmeanw(x, a):
    ok = np.isfinite(x)
    return float(np.sum(x[ok] * a[ok]) / np.sum(a[ok])) if ok.any() else float("nan")


def stats_for(A, m_e, cls, bin_, area, lab, npatch, cos1, sh1, F1, qbar, qmax1, Au):
    a = area[cls]
    b = bin_[cls]
    Ac = A[cls]
    lc = lab[cls]
    ref = float(np.sum(Ac * a) / a.sum())
    o = np.argsort(-Ac)
    ca = np.cumsum(a[o]) / a.sum()
    keep = o[np.searchsorted(ca, 0.01) + 1:]
    trim = float(np.sum(Ac[keep] * a[keep]) / a[keep].sum())
    contrast = Ac / patch_mean(Ac, a, lc, npatch)[lc]
    Af = Ac / m_e[cls]
    fcontrast = Af / patch_mean(Af, a, lc, npatch)[lc]
    cs = cos1[cls]
    sc = sh1[cls]
    litc = (cs >= 0.1) & ~sc
    M1 = np.where(litc, F1[cls] / np.where(litc, cs, 1.0), np.nan)
    m1ref = nanmeanw(M1, a)
    lit_any = (cs > 0) & ~sc
    strata = [("dark", ~lit_any), ("cos 0-0.3", lit_any & (cs < 0.3)),
              ("cos 0.3-0.7", lit_any & (cs >= 0.3) & (cs < 0.7)), ("cos 0.7-1", lit_any & (cs >= 0.7))]
    out = {"area_m2": float(a.sum()), "mean_A": ref, "trimmed_mean_A": trim, "M1_mean": m1ref,
           "M1_included_area": float(a[litc].sum() / a.sum()), "bins": [],
           "strata_area": {n: float(a[s].sum() / a.sum()) for n, s in strata}}
    smeans = {n: (float(np.sum(Ac[s] * a[s]) / a[s].sum()) if s.any() else float("nan")) for n, s in strata}
    q_rel = qbar[cls] / qmax1
    for k in range(7):
        s = b == k
        if not s.any():
            out["bins"].append(None)
            continue
        ak = a[s]
        Ak = Ac[s]
        med, p16, p84 = wquant(Ak, ak, [0.5, 0.16, 0.84])
        rec = {"area_share": float(ak.sum() / a.sum()), "A_mean": float(np.sum(Ak * ak) / ak.sum()),
               "A_median": med, "A_p16": p16, "A_p84": p84,
               "active_10pct": float(ak[Ak >= 0.1 * trim].sum() / ak.sum()),
               "off_1e-6": float(ak[Ak < 1e-6 * Au].sum() / ak.sum()),
               "off_1pct": float(ak[Ak < 1e-2 * Au].sum() / ak.sum()),
               "local_contrast_mean": nanmeanw(contrast[s], ak),
               "local_contrast_median": wquant(contrast[s], ak, [0.5])[0],
               "field_contrast_mean": nanmeanw(fcontrast[s], ak),
               "field_mean": float(np.sum(Af[s] * ak) / ak.sum()),
               "illum_diurnal_mean_q_over_qmax": float(np.sum(q_rel[s] * ak) / ak.sum()),
               "lit_fraction_103310": float(ak[lit_any[s]].sum() / ak.sum()),
               "M1_F_over_cos_mean": nanmeanw(M1[s], ak),
               "M1_included_fraction": float(ak[litc[s]].sum() / ak.sum()), "strata": {}}
        rec["A_rel"] = rec["A_mean"] / ref
        rec["A_median_rel"] = med / ref
        rec["M1_rel"] = rec["M1_F_over_cos_mean"] / m1ref
        for n, st in strata:
            t = s & st
            rec["strata"][n] = (float(np.sum(Ac[t] * a[t]) / a[t].sum() / smeans[n])
                                if t.any() and a[t].sum() > 1e3 and smeans[n] > 0 else float("nan"))
        out["bins"].append(rec)
    return out


def vicinity(out, cen, area, bin_, seen, Fs, radius=200.0):
    d = np.array([np.cos(np.radians(TD_LAT)) * np.cos(np.radians(TD_LON)),
                  np.cos(np.radians(TD_LAT)) * np.sin(np.radians(TD_LON)), np.sin(np.radians(TD_LAT))])
    r = np.linalg.norm(cen, axis=1)
    c = (cen @ d) / r
    # the facet on the ray to (lon, lat): nearest centroid to the ray within the outer surface layer (fixed 2026-11-21;
    # before, the highest centroid within 0.3 deg = 10 m of the ray was taken, up to 10 m off)
    perp = r * np.sqrt(np.maximum(0.0, 1.0 - c * c))
    cand = np.nonzero((c > 0) & (perp < 4.0))[0]
    cand = cand[r[cand] > r[cand].max() - 5.0]
    td = cand[np.argmin(perp[cand])]
    last = np.array([1.687560e3, 3.653001e2, 1.067887e3])       # spacecraft, last row of base_q002/COPS_model.dat
    say("touchdown facet %d: %.2f m from the lon/lat ray, %.2f m from the last spacecraft position (altitude 5.22 m)"
        % (td, perp[td], np.linalg.norm(cen[td] - last)))
    p0 = cen[td].copy()
    dd = cen - p0
    near = np.nonzero(np.einsum("ij,ij->i", dd, dd) < (1.6 * radius) ** 2)[0]
    nodes = np.memmap(GEO / "mesh_nodes.f8", dtype="<f8", mode="r").reshape(-1, 3)
    faces = np.memmap(GEO / "mesh_faces.i4", dtype="<i4", mode="r").reshape(-1, 3)
    fv = np.array(faces[near])
    if fv.min() == 1 and fv.max() == len(nodes):
        fv -= 1
    tri = np.array(nodes[fv.ravel()]).reshape(-1, 3, 3)
    if np.max(np.abs(tri.mean(axis=1) - cen[near])) > 1.0:
        raise SystemExit("mesh faces do not match the element centroids")
    n = np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0])
    core = np.linalg.norm(cen[near] - p0, axis=1) < 50
    axis = n[core].sum(0)
    axis /= np.linalg.norm(axis)
    if axis @ (p0 / np.linalg.norm(p0)) < 0:
        axis = -axis
    north = np.array([0.0, 0.0, 1.0])
    ey = north - (north @ axis) * axis
    ey /= np.linalg.norm(ey)
    ex = np.cross(ey, axis)
    rel = tri - p0
    xy = np.stack([rel @ ex, rel @ ey], axis=-1)
    inside = np.all(np.abs(xy.mean(axis=1)) < radius * 1.08, axis=1)
    order = np.argsort((rel.mean(axis=1) @ axis)[inside])          # painter's algorithm along the view axis
    sel = near[inside][order]
    xy = xy[inside][order]
    vis = seen[sel]
    info = {"touchdown_facet": int(td), "touchdown_lonlat": [TD_LON, TD_LAT], "window_radius_m": radius,
            "facets_drawn": int(len(sel)), "area_drawn_m2": float(area[sel].sum()),
            "area_share_by_bin": [float(area[sel][bin_[sel] == k].sum() / area[sel].sum()) for k in range(7)]}
    for tag, F in Fs.items():
        allv = np.concatenate([F[key][sel][vis] for key, _, _ in RUNS])
        vmax = float(10 ** np.ceil(np.log10(np.percentile(allv[allv > 0], 99.5))))
        vmin = vmax / 1e4
        fig, axes = plt.subplots(2, 3, figsize=(14.2, 9.8), dpi=300)
        panels = [(None, "(a) slope: angle to the local effective gravity")] + [
            (key, "(%s) source rate, %s" % ("bcdef"[i], lab)) for i, (key, lab, _) in enumerate(RUNS)]
        for ax, (key, title) in zip(axes.ravel(), panels):
            if key is None:
                pc = PolyCollection(xy, array=bin_[sel].astype(float), cmap=ListedColormap(BIN_COL), linewidths=0,
                                    antialiased=False, rasterized=True)
                pc.set_clim(-0.5, 6.5)
                ax.add_collection(pc)
                cb = fig.colorbar(pc, ax=ax, fraction=0.046, pad=0.02, ticks=range(7))
                cb.ax.set_yticklabels(BIN_LAB, fontsize=7)
            else:
                v = F[key][sel]
                colv = np.where(vis, np.clip(v, vmin, vmax), np.nan)
                cmap = plt.get_cmap("magma").copy()
                cmap.set_bad("#c8c8c8")
                pc = PolyCollection(xy, array=np.ma.masked_invalid(colv), cmap=cmap, norm=LogNorm(vmin, vmax),
                                    linewidths=0, antialiased=False, rasterized=True)
                ax.add_collection(pc)
                cb = fig.colorbar(pc, ax=ax, fraction=0.046, pad=0.02, extend="both")
                cb.set_label("F$_{H_2O}$ (m$^{-2}$ s$^{-1}$)", fontsize=8)
                cb.ax.tick_params(labelsize=7)
            ax.plot(0, 0, marker="+", ms=14, mew=2.0, color="white", zorder=5)
            ax.plot(0, 0, marker="+", ms=12, mew=0.9, color="black", zorder=6)
            ax.set_xlim(-radius, radius)
            ax.set_ylim(-radius, radius)
            ax.set_aspect("equal")
            ax.set_title(title, loc="left", fontsize=9)
            ax.tick_params(labelsize=7)
            ax.set_xlabel("m (east)", fontsize=8)
            ax.set_ylabel("m (north)", fontsize=8)
        fig.suptitle("Vicinity of the touchdown point (+), %s UTC 2016-09-30, %s" % (tag, LAW), fontsize=11)
        fig.tight_layout(rect=(0, 0, 1, 0.97))
        stem = "12M_q002_top5_touchdown_source_rate_%s" % tag.replace(":", "")
        for ext in ("png", "eps"):
            fig.savefig(out / ("%s.%s" % (stem, ext)), dpi=300)
        (out / ("%s_caption.txt" % stem)).write_text(
            "Vicinity (%d m x %d m) of Rosetta's touchdown point (+; lon %.3f deg, lat %.3f deg) on the SHAP8 12M shape, "
            "seen along the mean surface normal within 50 m of it (north up; where layers overlap the upper one is "
            "drawn). (a) Effective gravitational slope class (angle between the inward normal and the local effective "
            "gravity). (b-f) H2O source rate F = A q at %s UTC of the five selected fits (%s: night and cast "
            "shadow q = q_min, lit q_min + (q_max - q_min) cos SZA); one logarithmic colour scale %.0e-%.0e m^-2 s^-1 "
            "for all five (values outside clipped to it); grey: facets no fitted observation sees (model "
            "extrapolation only).\n" % (2 * radius, 2 * radius, TD_LON, TD_LAT, tag, LAW, vmin, vmax))
        plt.close(fig)
        say("vicinity figure %s: %d facets" % (tag, len(sel)))
    return info


def line_panels(res, classes, metric, fname, out, ylabel, caption, second=None, names=None):
    fig, axes = plt.subplots(1, len(classes), figsize=(4.9 * len(classes), 4.4), dpi=300, sharey=True)
    axes = np.atleast_1d(axes)
    x = np.arange(7)
    for ax, (cls, title) in zip(axes, classes):
        for i, (key, lab, _) in enumerate(RUNS):
            B = res[key][cls]["bins"]
            ax.plot(x, [b[metric] if b else np.nan for b in B], "-o", ms=3.5, lw=1.4, color=figures.PALETTE[i],
                    label=lab)
            if second:
                ax.plot(x, [b[second] if b else np.nan for b in B], ls=(0, (3, 2)), lw=1.0, marker="s", mfc="white",
                        ms=3.5, color=figures.PALETTE[i])
        ax.axhline(1, color="#52514e", lw=0.7)
        ax.set_yscale("log")
        ax.set_xticks(x)
        ax.set_xticklabels(BIN_LAB, fontsize=7.5, rotation=30)
        ax.set_title(title, loc="left", fontsize=9)
        ax.grid(True, color=GRID_C, lw=0.5)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
    axes[0].set_ylabel(ylabel)
    h, l = axes[0].get_legend_handles_labels()
    if names:
        h += [plt.Line2D([], [], color=INK, ls="-", marker="o", ms=3.5),
              plt.Line2D([], [], color=INK, ls=(0, (3, 2)), marker="s", mfc="white", ms=3.5)]
        l += list(names)
    fig.legend(h, l, fontsize=7, frameon=False, loc="lower center", ncol=4)
    fig.tight_layout(rect=(0, 0.17, 1, 1))
    for ext in ("png", "eps"):
        fig.savefig(out / ("%s.%s" % (fname, ext)), dpi=300)
    (out / ("%s_caption.txt" % fname)).write_text(caption + "\n")
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--out", required=True)
    ap.add_argument("--seen-gauges", default=None)
    ap.add_argument("--maps-only", action="store_true")
    ap.add_argument("--stats-only", action="store_true")
    a = ap.parse_args()
    out = Path(a.out)
    if (out / "slope_activity.json").exists():
        raise SystemExit("slope_activity.json exists in %s" % out)
    out.mkdir(parents=True, exist_ok=True)
    met = json.load(open(PL))["runs"]
    nf, cen, area = elements()
    with np.load(Path(met[A_RUN]["fit_dir"]) / "nonneg_facets.npz", allow_pickle=False) as z:
        seen = np.asarray(z["seen"]).astype(bool)
        special = (np.asarray(z["TS"]) != 0) | (np.asarray(z["held"]) != 0)
        bin_ = np.asarray(z["bin"]).astype(np.int8)
    if len(seen) != nf or bin_.min() < 0 or bin_.max() != 6:
        raise SystemExit("nonneg_facets: unexpected facet or bin data")
    say("elements, seen mask (%.4f of the area), slope classes" % (area[seen].sum() / area.sum()))

    T = map_times()
    lab1 = [k for k, v in T.items() if v["utc"].endswith("10:33:10")][0]
    lab2 = [k for k, v in T.items() if v["utc"].endswith("10:39:26")][0]
    qbar = np.zeros(nf)
    keep = {}
    for lab in sorted(T):
        q, qmax, cos, sh = qfac(lab, T[lab]["rh_m"], nf)
        qbar += q / len(T)
        if lab == lab1:
            keep["1"] = (q, qmax, cos, sh)
        elif lab == lab2:
            keep["2"] = (q, qmax)
        del cos, sh
    q1, qmax1, cos1, sh1 = keep["1"]
    q2, _ = keep["2"]
    say("illumination at %d map times; q_max(10:33:10) %.6f" % (len(T), qmax1))

    res = {"schema": "q002_slope_activity_v1", "law": LAW, "qmax_103310": float(qmax1), "map_times": len(T),
           "definition": __doc__, "runs": {}, "checks": {}}
    A_all = {}
    m_all = {}
    Fs = {"10:33:10": {}, "10:39:26": {}}
    for key, lab, fk in RUNS:
        with np.load(Path(met[key]["ft_dir"]) / "finetune_state.npz", allow_pickle=False) as z:
            A = np.asarray(z["A_final"], np.float64)
        F1 = A * q1
        fp = ROOT / ("figures_cluster/12M/full_set_q002_%s/checks/facts.json" % fk)
        facts = json.load(open(fp))
        Qf = facts["Q"]["10:33:10"][0]
        Q = float(np.sum(F1 * area))
        res["checks"]["Q_103310_" + key] = [Q, Qf, Q / Qf - 1]
        if abs(Q / Qf - 1) > 1e-9:
            raise SystemExit("%s: Q(10:33:10) %.6e from A q differs from the full set's %.6e" % (key, Q, Qf))
        mb = np.ones(7)
        row = {r["key"]: r for r in csv.DictReader(open(ROOT / "figures_cluster/12M/"
               "model_vs_observations_nonneg_options_stage2_q002/summary.csv"))}[key.replace("_reg", "")]
        for k, c in enumerate(("40_50", "50_60", "60_70", "70_80", "80_90", "90_180"), start=1):
            mb[k] = float(row["mult_" + c])
        m_all[key] = mb
        A_all[key] = A
        Fs["10:33:10"][key] = F1
        Fs["10:39:26"][key] = A * q2
        say("%s: A loaded, Q check %.1e" % (key, Q / Qf - 1))

    if not a.stats_only:
        res["vicinity"] = vicinity(out, cen, area, bin_, seen, Fs)
    del cen
    for d in Fs.values():
        if d is not Fs["10:33:10"]:
            d.clear()
    if a.maps_only:
        (out / "slope_activity_maps.json").write_text(json.dumps(res.get("vicinity", {}), indent=1))
        return

    # visibility classes
    pm = Path(str(met[A_RUN]["ft_dir"]).split("/c5/")[0]) / "contrib_phase_c5/phase_means.npz"
    with np.load(pm, allow_pickle=False) as z:
        ng = (np.asarray(z["NG_pre"]) > 0) | (np.asarray(z["NG_late"]) > 0)
        rg = (np.asarray(z["RG_pre"]) > 0) | (np.asarray(z["RG_gap"]) > 0) | (np.asarray(z["RG_late"]) > 0)
    res["checks"]["fit_seen_area"] = float(area[seen].sum() / area.sum())
    res["checks"]["phase_union_area"] = float(area[ng | rg].sum() / area.sum())
    res["checks"]["seen_minus_phase_union_area"] = float(area[seen & ~(ng | rg)].sum() / area.sum())
    res["checks"]["phase_union_minus_seen_area"] = float(area[(ng | rg) & ~seen].sum() / area.sum())
    classes = {"any": seen, "NG": ng, "RG": rg, "both": ng & rg, "NG_only": ng & ~rg, "RG_only": rg & ~ng}
    if a.seen_gauges:
        with np.load(a.seen_gauges, allow_pickle=False) as z:
            U = {k: np.unpackbits(z[k], count=nf).astype(bool) for k in ("ng_obs", "rg_obs", "both_obs")}
        classes["both_concurrent"] = U["both_obs"]
        classes["both_different_times"] = (U["ng_obs"] & U["rg_obs"]) & ~U["both_obs"]
        res["checks"]["kernel_pass_vs_phase"] = {
            "ng_mismatch_area": float(area[U["ng_obs"] != ng].sum() / area.sum()),
            "rg_mismatch_area": float(area[U["rg_obs"] != rg].sum() / area.sum())}
    full = {k + "_all": classes[k] for k in ("any", "both")}
    classes = {k: v & ~special for k, v in classes.items()}
    classes.update(full)
    res["special"] = {"facets": int(special.sum()), "area_m2": float(area[special].sum())}
    wseen = (q1 * area)[seen]
    res["class_area"] = {k: float(area[v].sum() / area.sum()) for k, v in classes.items()}
    say("classes: %s; checks %s" % (res["class_area"], res["checks"]))

    with np.load(GRID, allow_pickle=False) as z:
        patch = np.asarray(z["facet_label"]).astype(np.int32)
    npatch = int(patch.max()) + 1
    for key, lab, _ in RUNS:
        me = m_all[key][bin_]
        F1 = Fs["10:33:10"][key]
        Au = float(np.sum(A_all[key][seen] * wseen) / wseen.sum())
        res["runs"][key] = {"label": lab, "multipliers": m_all[key].tolist(), "A_u": Au}
        for cname, cmask in classes.items():
            res["runs"][key][cname] = stats_for(A_all[key], me, cmask, bin_, area, patch, npatch, cos1, sh1, F1,
                                                qbar, qmax1, Au)
        say("stats %s" % key)
    (out / "slope_activity.json").write_text(json.dumps(res, indent=1))

    cl2 = [("any", "(a) seen by at least one gauge"), ("both", "(b) seen by both gauges (same or different times)")]
    if a.seen_gauges:
        cl2 += [("both_concurrent", "(c) both in the same observation"),
                ("both_different_times", "(d) both, only at different times")]
    R = res["runs"]
    line_panels(R, cl2, "A_rel", "12M_q002_top5_slope_activity", out, "relative to the class mean",
                "Illumination-corrected activity against slope (%s): solid, the activity A = F/q of the model's own "
                "illumination law (defined on every facet), area-weighted mean per effective-gravitational-slope class "
                "divided by the mean over the visibility class; dashed, the same for F/cos SZA at 10:33:10 UTC (lit "
                "facets with cos SZA >= 0.1 only). 1 = no slope dependence. The 275 held-template and late fine-tuned facets (367 m^2) are left out." % LAW, second="M1_rel", names=("A = F/q (all facets)", "F/cos SZA (lit, cos >= 0.1)"))
    line_panels(R, cl2, "local_contrast_mean", "12M_q002_top5_slope_local_contrast", out,
                "A / mean A of the 150-m patch",
                "Slope effect separated from location (%s): solid, the area-weighted mean over each slope class of A "
                "divided by the mean A of the facet's 150-m patch (same visibility class); dashed, the same for the "
                "fitted field without the slope multiplier, A/m_b. 1 = a slope class is as active as its neighbourhood."
                % LAW, second="field_contrast_mean", names=("A", "A / m_b (no multiplier)"))

    fig, axes = plt.subplots(1, 2, figsize=(10.4, 3.8), dpi=300)
    x = np.arange(7)
    r0 = R[RUNS[0][0]]
    for j, (cls, title) in enumerate(cl2[:2]):
        B = r0[cls]["bins"]
        for ax, met_, yl in ((axes[0], "illum_diurnal_mean_q_over_qmax", "diurnal mean q / q$_{max}$"),
                             (axes[1], "lit_fraction_103310", "lit area fraction, 10:33:10 UTC")):
            ax.plot(x, [b[met_] if b else np.nan for b in B], marker="o", ms=3.5, color=("#0072B2", "#D55E00")[j],
                    label=title[4:])
            ax.set_ylabel(yl)
            ax.set_xticks(x)
            ax.set_xticklabels(BIN_LAB, fontsize=7.5, rotation=30)
            ax.grid(True, color=GRID_C, lw=0.5)
            ax.set_ylim(bottom=0)
            for s in ("top", "right"):
                ax.spines[s].set_visible(False)
    axes[0].legend(fontsize=7, frameon=False)
    fig.tight_layout()
    for ext in ("png", "eps"):
        fig.savefig(out / ("12M_q002_top5_slope_illumination.%s" % ext), dpi=300)
    (out / "12M_q002_top5_slope_illumination_caption.txt").write_text(
        "The illumination each slope class receives (geometry only, the same for all fits): left, the area-weighted "
        "mean illumination factor q/q_max over the 13 map times of 2016-09-30; right, the lit (cos SZA > 0, not in "
        "cast shadow) area fraction at 10:33:10 UTC.\n")
    plt.close(fig)

    fig, axes = plt.subplots(1, 5, figsize=(16.5, 3.9), dpi=300, sharey=True)
    sc = {"dark": "#0b4f8a", "cos 0-0.3": "#56B4E9", "cos 0.3-0.7": "#E69F00", "cos 0.7-1": "#D55E00"}
    for ax, (key, lab, _) in zip(axes, RUNS):
        for st, c in sc.items():
            ax.plot(x, [b["strata"][st] if b else np.nan for b in R[key]["any"]["bins"]], "-o", color=c, ms=3,
                    lw=1.2, label=st)
        ax.axhline(1, color="#52514e", lw=0.7)
        ax.set_yscale("log")
        ax.set_title(lab, loc="left", fontsize=8)
        ax.set_xticks(x)
        ax.set_xticklabels(BIN_LAB, fontsize=6.5, rotation=40)
        ax.grid(True, color=GRID_C, lw=0.5)
    axes[0].set_ylabel("mean A / stratum mean")
    axes[0].legend(fontsize=6.6, frameon=False)
    fig.tight_layout()
    for ext in ("png", "eps"):
        fig.savefig(out / ("12M_q002_top5_slope_strata.%s" % ext), dpi=300)
    (out / "12M_q002_top5_slope_strata_caption.txt").write_text(
        "Is the slope dependence an illumination effect? Area-weighted mean activity per slope class inside groups of "
        "facets with the same illumination at 10:33:10 UTC (dark; lit with cos SZA 0-0.3, 0.3-0.7, 0.7-1), divided by "
        "the group's mean; facets seen by at least one gauge (%s).\n" % LAW)
    plt.close(fig)
    say("done")


if __name__ == "__main__":
    main()
