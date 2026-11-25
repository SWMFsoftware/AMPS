#!/usr/bin/env python3
# ==== version-controlled copy: header added 2026-11-21, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/figures_cluster/scripts/q002_slopetest_report.py (on Voyager: /data001/highresolu_vt/cops_output/figures_cluster/scripts/q002_slopetest_report.py)
# SHA256 of the run copy = SHA256 of this file without these header lines: c5ad6cf70b4c9152324afd414c3e600b4341d1322bd75cb5eb3855cdb16d3626
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: table of the slope-effect control fits (decision 3c): costs, rms, switched-off share, ordering
# RUN: jobs/q002_slopetest_report.sbatch
# ==== end of header ====
"""Table of the slope-effect control fits (decision 3c, 2026-11-21; PLAN_3c.md): the finished variants and the
reference: objective, data cost, penalties, rms, termination, multipliers, switched-off share by decision 1b, and the
ordering that holds at global optima (free_min <= reference <= only_X <= none in the objective).

    q002_slopetest_report.py --runroot <slopetest run root> --out <new folder>
"""
import argparse
import csv
import json
from pathlib import Path

import numpy as np

import q002_slope_activity as geo       # q = F / A per facet at 10:33:10 UTC

R = Path("/data001/highresolu_vt/cops_output")
REF = R / ("runs/shap8_012m_lategrid_dz1.5/nonneg/localnn/stage2r_q002_f10_pm0.01_tploff_20261118T162407Z/c1/fit/final")
AA = R / "figures_cluster/12M/q002_decisions_20261121/active_area.json"
RH, F_E = 3.833406592273387, 1.0936266398665254e21      # at 10:33:10 UTC, q002_plausibility/metrics.json
QMAX = 600.0 * RH ** -4.2143229                         # the H2O law's q_max there
BINS = ("40_50", "50_60", "60_70", "70_80", "80_90", "90_180")
ONLY = ["only_" + b for b in BINS]


def q_103310(n):
    if not hasattr(q_103310, "q"):
        t = geo.map_times()
        k = [k for k, v in t.items() if v["utc"].endswith("10:33:10")][0]
        q_103310.q = geo.qfac(k, t[k]["rh_m"], n)[0]
    return q_103310.q


def row(name, seg, T):
    man = json.loads((seg / "fit_manifest.json").read_text())
    c = man["costs"]
    with np.load(seg / "nonneg_facets.npz", allow_pickle=False) as z:
        seen, area, A, bins = z["seen"], z["area"], z["A"], z["bin"].astype(int)
    off = seen & (A <= T)
    r = dict(variant=name, objective=c["objective"], data_cost=c["data_cost"], smoothness=c["penalty"],
             P_m=c.get("mult_penalty", 0.0), crossing=c["crossing_cost"], rms_NG=c["rms_ln_NG"], rms_RG=c["rms_ln_RG"],
             late_rms_NG=c["late_rms_ln_NG"], late_rms_RG=c["late_rms_ln_RG"], termination=man["message"],
             nfev=man["nfev_total"], optimality=man["optimality"], segments=man["segment"],
             c_at_bound=man["active_bounds"]["c_at_lower_bound"],
             m_at_bound=" ".join(man["active_bounds"]["multipliers_at_lower_bound"]),
             off_seen_1b=float(area[off].sum() / area[seen].sum()))
    for tag, f in (("T0p1", 0.1), ("T10", 10.0)):
        r["off_seen_" + tag] = float(area[seen & (A <= f * T)].sum() / area[seen].sum())
    r["passes_1b"] = r["off_seen_1b"] <= 0.5
    r["F_bound_over_F_E"] = float(A.max() * QMAX / F_E)     # max F <= max A q_max: every facet as if subsolar
    r["A_weighted_seen_share"] = float((A * area)[seen].sum() / (A * area).sum())
    F = A * q_103310(len(A))
    Fa = F * area
    r["Q_103310"] = float(Fa.sum())
    r["Q_seen_share"] = float(Fa[seen].sum() / Fa.sum())
    r["Q_share_over_F_E"] = float(Fa[F > F_E].sum() / Fa.sum())
    r["F_max_over_F_E"] = float(F.max() / F_E)
    share = np.array([Fa[seen & (bins == b)].sum() for b in range(7)]) / Fa[seen].sum()
    r["top_class"] = (("00_40",) + BINS)[int(share.argmax())]
    r["top_class_share_of_seen_production"] = float(share.max())
    for b, m in zip(BINS, man["multipliers"][1:]):
        r["m_" + b] = m
    r["off_by_bin"] = " ".join("%.4f" % (area[off & (bins == b)].sum() / area[seen & (bins == b)].sum()) for b in range(7))
    return r


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--runroot", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--multistart", help="run root of the multi-start refits (user 2026-11-22): adds its fits")
    a = ap.parse_args()
    out = Path(a.out)
    out.mkdir()
    T = json.loads(AA.read_text())["T"]
    rows, waiting = [row("reference", REF, T)], []
    for v in ["none"] + ONLY + ["free_min"]:
        seg = Path(a.runroot) / v / "fit/final"
        if seg.exists() and any((Path(a.runroot) / v / "fit").glob("direct_check_seg*")):
            rows.append(row(v, seg.resolve(), T))
        else:
            waiting.append(v)
    skipped = {}
    if a.multistart:                                        # from_<parent>: six free >= 1, started at the parent's end
        for v in ONLY[:5] + ["free_min"]:
            o = Path(a.multistart) / ("from_" + v)
            if (o / "SKIPPED.txt").exists():
                skipped["ms_" + v] = (o / "SKIPPED.txt").read_text().strip()
            elif (o / "fit/final").exists() and any((o / "fit").glob("direct_check_seg*")):
                rows.append(row("ms_" + v, (o / "fit/final").resolve(), T))
            else:
                waiting.append("ms_" + v)
    d = {r["variant"]: r for r in rows}
    pairs = [(x, "none") for x in ONLY] + [("reference", x) for x in ONLY + ["none"]] + [("free_min", "reference")]
    pairs += [("ms_" + x, x) for x in ONLY[:5] + ["free_min"]]      # a refit from x cannot end above x
    pairs += [("free_min", "ms_" + x) for x in ONLY[:5]]            # >= 1 states are feasible for free_min's problem too
    order = [dict(lower=l, upper=u, ok=d[l]["objective"] <= d[u]["objective"] * (1 + 1e-6))
             for l, u in pairs if l in d and u in d]
    with open(out / "slopetest.csv", "w", newline="") as f:
        w = csv.DictWriter(f, list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    (out / "slopetest.json").write_text(json.dumps(dict(T=T, rows=rows, waiting=waiting, skipped=skipped, ordering=order), indent=1) + "\n")
    base = d.get("none")
    for r in rows:
        g = ("  vs none: objective %+.3f %% data %+.3f %%" % (100 * (r["objective"] / base["objective"] - 1),
             100 * (r["data_cost"] / base["data_cost"] - 1))) if base else ""
        print("%-12s obj %.6e data %.6e rms all %.4f/%.4f late %.4f/%.4f off %.1f %% nfev %d %s%s" % (
            r["variant"], r["objective"], r["data_cost"], r["rms_NG"], r["rms_RG"], r["late_rms_NG"], r["late_rms_RG"],
            100 * r["off_seen_1b"], r["nfev"], r["termination"].split()[0], g))
        print("   off by bin", r["off_by_bin"], "| off at T/10 %.1f %% at 10 T %.1f %% | 1b %s | max F/F_E <= %.3f" % (
            100 * r["off_seen_T0p1"], 100 * r["off_seen_T10"], "pass" if r["passes_1b"] else "fail", r["F_bound_over_F_E"]))
    for r in rows:
        print("%-16s top class %s %.1f %% of seen production | Q %.3e seen %.1f %% above F_E %.2f %% max F/F_E %.3f | m %s" % (
            r["variant"], r["top_class"], 100 * r["top_class_share_of_seen_production"], r["Q_103310"],
            100 * r["Q_seen_share"], 100 * r["Q_share_over_F_E"], r["F_max_over_F_E"],
            " ".join("%.4g" % r["m_" + b] for b in BINS)))
    for k, why in skipped.items():
        print("skipped", k, why)
    print("ordering violations:", [o for o in order if not o["ok"]] or "none", "| waiting:", waiting)


if __name__ == "__main__":
    main()
