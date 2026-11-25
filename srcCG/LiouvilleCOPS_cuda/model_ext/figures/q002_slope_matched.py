#!/usr/bin/env python3
# ==== version-controlled copy: header added 2026-11-21, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/figures_cluster/scripts/q002_slope_matched.py (on Voyager: /data001/highresolu_vt/cops_output/figures_cluster/scripts/q002_slope_matched.py)
# SHA256 of the run copy = SHA256 of this file without these header lines: efe28cbc9cc21a4022df0e348f6105a4c70158e9b1bce141934474253f8c0262
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: slope effect on activity matched in location and illumination (2026-11-21)
# RUN: see OUTPUT_VISUALIZATION.md 3.21
# ==== end of header ====
"""Slope effect on the surface activity, matched in location and illumination (2026-11-21; user: "propose a better metric
to account for the effect of the illumination and separate the effect of the angle to the local gravity").

For every cell = (150-m patch of the option-C mesh) x (illumination class) that contains slope bin b and 0-40 deg facets:
    r = <X>_b,cell / <X>_0-40,cell          (area-weighted means in the cell)
and per bin the weighted median (and 16 / 84 %) of r over the cells, weight = the smaller of the two areas.
  X = A      illumination class = quintile of the diurnal mean q / q_max (the facet's exposure over the rotation):
             the proposed metric (same neighbourhood, same exposure, illumination law divided out exactly);
  X = A/m_b  the same without the fitted slope multiplier: 1 means the fit has no slope effect other than m_b;
  X = F/cos  the user's metric at 10:33:10, cells = patch x cos SZA class (0.1-0.3, 0.3-0.7, 0.7-1; lit only).
Classes: seen by at least one gauge; by both gauges (any time); by both in the same observation (seen_gauges.npz).
The 275 held-template / late fine-tuned facets are left out. Also writes slope_statistics.csv (every bin statistic of
pass2/slope_activity.json and the matched ratios) and the two-panel figure.

    q002_slope_matched.py --dir <slope_activity dir>      (writes <dir>/matched/, never overwrites)
"""
import argparse
import csv
import json
from pathlib import Path

import numpy as np

import q002_slope_activity as Q

MIN_AREA = 5.0                      # m2 of each slope class in a cell


def matched(X, a, b, cell):
    n = int(cell.max()) + 1
    s0 = b == 0
    n0 = np.bincount(cell[s0], weights=(X * a)[s0], minlength=n)
    d0 = np.bincount(cell[s0], weights=a[s0], minlength=n)
    out = [{"ratio": 1.0}]
    for k in range(1, 7):
        s = b == k
        nk = np.bincount(cell[s], weights=(X * a)[s], minlength=n)
        dk = np.bincount(cell[s], weights=a[s], minlength=n)
        ok = (dk >= MIN_AREA) & (d0 >= MIN_AREA) & (n0 > 0)
        if not ok.any():
            out.append(None)
            continue
        r = (nk[ok] / dk[ok]) / (n0[ok] / d0[ok])
        med, lo, hi = Q.wquant(r, np.minimum(dk[ok], d0[ok]), [0.5, 0.16, 0.84])
        out.append({"ratio": med, "p16": lo, "p84": hi, "cells": int(ok.sum()),
                    "area_coverage": float(dk[ok].sum() / dk.sum())})
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", required=True)
    D = Path(ap.parse_args().dir)
    O = D / "matched"
    if O.exists():
        raise SystemExit("%s exists" % O)
    met = json.load(open(Q.PL))["runs"]
    P2 = json.load(open(D / "pass2/slope_activity.json"))
    nf, cen, area = Q.elements()
    del cen
    with np.load(Path(met[Q.A_RUN]["fit_dir"]) / "nonneg_facets.npz") as z:
        seen = z["seen"].astype(bool)
        bin_ = z["bin"].astype(int)
        special = (z["TS"] != 0) | (z["held"] != 0)
    with np.load(Q.GRID) as z:
        patch = z["facet_label"].astype(np.int64)
    with np.load(D / "seen_gauges.npz") as z:
        U = {k: np.unpackbits(z[k], count=nf).astype(bool)
             for k in ("ng_obs", "rg_obs", "both_obs", "ng_lit", "rg_lit", "both_lit")}
    T = Q.map_times()
    qb = np.zeros(nf)
    for t in sorted(T):
        q, qm, cos, sh = Q.qfac(t, T[t]["rh_m"], nf)
        qb += q / qm / len(T)
        if T[t]["utc"].endswith("10:33:10"):
            q1, c1, s1 = q, cos, sh
    edges = Q.wquant(qb[seen], area[seen], [0.2, 0.4, 0.6, 0.8])
    expo = np.searchsorted(edges, qb)
    lit = (c1 >= 0.1) & ~s1
    cc = np.digitize(c1, [0.3, 0.7])
    classes = {"any": seen & ~special, "both": U["ng_obs"] & U["rg_obs"] & ~special,
               "both_concurrent": U["both_obs"] & ~special}
    tot = area.sum()
    vis = {k: float(area[v].sum() / tot) for k, v in U.items()}
    vis["any"] = float(area[seen].sum() / tot)
    vis["both_any_time"] = float(area[U["ng_obs"] & U["rg_obs"]].sum() / tot)
    vis["both_only_different_times_facets"] = int((U["ng_obs"] & U["rg_obs"] & ~U["both_obs"]).sum())
    vis["seen_lit_at_least_once_of_seen"] = float(area[(U["ng_lit"] | U["rg_lit"])].sum() / area[seen].sum())
    res = {"schema": "q002_slope_matched_v1", "definition": __doc__, "exposure_edges": edges, "visibility": vis,
           "runs": {}}
    O.mkdir()
    rows = []
    for key, lab, _ in Q.RUNS:
        with np.load(Path(met[key]["ft_dir"]) / "finetune_state.npz") as z:
            A = z["A_final"].astype(float)
        m = np.array(P2["runs"][key]["multipliers"])[bin_]
        res["runs"][key] = {"label": lab}
        for cn, s in classes.items():
            cell = patch[s] * 5 + expo[s]
            l = s & lit
            r = {"A": matched(A[s], area[s], bin_[s], cell),
                 "A_over_m": matched((A / m)[s], area[s], bin_[s], cell),
                 "F_over_cos": matched((A * q1)[l] / c1[l], area[l], bin_[l], patch[l] * 3 + cc[l])}
            res["runs"][key][cn] = r
            B = P2["runs"][key][cn]["bins"]
            for k in range(7):
                row = {"fit": lab, "visibility": cn, "slope_bin": Q.BIN_LAB[k].replace("°", " deg"),
                       "multiplier": P2["runs"][key]["multipliers"][k]}
                row.update({n: v for n, v in B[k].items() if n != "strata"})
                row["A_over_0_40"] = B[k]["A_mean"] / B[0]["A_mean"]
                for n, v in r.items():
                    row["matched_" + n] = v[k]["ratio"] if v[k] else ""
                    row["matched_%s_p16" % n] = v[k].get("p16", "") if v[k] else ""
                    row["matched_%s_p84" % n] = v[k].get("p84", "") if v[k] else ""
                rows.append(row)
            if cn != "both":
                print("%-30s %-16s" % (key, cn))
                for n, v in r.items():
                    print("   %-11s" % n, " ".join("%8.3g" % x["ratio"] if x else "   -    " for x in v),
                          "| coverage", " ".join("%.2f" % x["area_coverage"] if x else "-" for x in v[1:]))
    print("visibility", vis)
    (O / "slope_matched.json").write_text(json.dumps(res, indent=1))
    with open(O / "slope_statistics.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)

    plt = Q.plt
    x = np.arange(7)
    fig, axes = plt.subplots(1, 2, figsize=(12.6, 5.2), dpi=300, sharey=True)
    for ax, (cn, title) in zip(axes, (("any", "(a) seen by at least one gauge"),
                                      ("both_concurrent", "(b) seen by both gauges in the same observation"))):
        for i, (key, lab, _) in enumerate(Q.RUNS):
            c = Q.figures.PALETTE[i]
            g = lambda n: [v["ratio"] if v else np.nan for v in res["runs"][key][cn][n]]
            ax.plot(x, g("A"), "-o", color=c, ms=4, lw=1.4, label=lab)
            ax.plot(x, g("F_over_cos"), "s", color=c, mfc="white", ms=6.5, mew=1.0)
            ax.plot(x, g("A_over_m"), ls=(0, (3, 2)), color=c, lw=1.0)
        ax.set_yscale("log")
        ax.axhline(1, color="#52514e", lw=0.6)
        ax.set_xticks(x)
        ax.set_xticklabels(Q.BIN_LAB, fontsize=8)
        ax.set_title(title, loc="left", fontsize=9.5)
        ax.grid(True, color=Q.GRID_C, lw=0.5)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
    axes[0].set_ylabel("ratio to 0-40° facets, same patch and illumination")
    h, l = axes[0].get_legend_handles_labels()
    h += [plt.Line2D([], [], color=Q.INK, marker="o", ms=4),
          plt.Line2D([], [], color=Q.INK, marker="s", mfc="white", ls=""),
          plt.Line2D([], [], color=Q.INK, ls=(0, (3, 2)))]
    l += ["A = F/q (proposed metric)", "F/cos SZA (lit facets only)", "A without the slope multiplier"]
    fig.legend(h, l, fontsize=7.6, frameon=False, loc="lower center", ncol=4)
    fig.tight_layout(rect=(0, 0.16, 1, 1))
    st = "12M_q002_top5_slope_matched"
    for e in ("png", "eps"):
        fig.savefig(O / ("%s.%s" % (st, e)), dpi=300)
    (O / (st + "_caption.txt")).write_text(
        "Effect of the angle to the local gravity on the activity, with location and illumination matched: each slope "
        "class against 0-40 deg facets in the same 150-m patch and the same class of illumination (median over the "
        "cells). Solid: A = F/q, cells matched in exposure over the rotation (proposed metric). Squares: F/cos SZA at "
        "10:33:10 UTC, lit facets with cos SZA >= 0.1, cells matched in cos SZA. Dashed: A without the fitted slope "
        "multiplier; 1 means that the fit has no slope effect other than its multipliers. (a) Seen by at least one "
        "gauge; (b) seen by both gauges in the same observation. No facet with a measurable area is seen by both only "
        "at different times. The 275 template / late fine-tuned facets (367 m^2) are left out.\n")


if __name__ == "__main__":
    main()
