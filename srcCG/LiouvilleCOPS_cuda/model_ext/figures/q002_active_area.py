#!/usr/bin/env python3
# ==== version-controlled copy: header added 2026-11-21, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/figures_cluster/scripts/q002_active_area.py (on Voyager: /data001/highresolu_vt/cops_output/figures_cluster/scripts/q002_active_area.py)
# SHA256 of the run copy = SHA256 of this file without these header lines: 6e9a46434cc5741406e32935460d2a99311889ba132b2a4bf0859d7523f1a9a7
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: switched-off share of the seen surface by decision 1b (2026-11-21)
# RUN: jobs/q002_active_area.sbatch
# ==== end of header ====
"""Switched-off share of the seen surface of the ten new-law fits, user decision 1b of 2026-11-21.

  active        : A without the imposed floor > T, T = 1 % of the median (one value per facet) of the srcCG map A_orig
                  (elements.bin AH2O) over the seen facets (the threshold of localnn_report.py / nonneg_alt_summary.py)
  imposed floor : option B only, m_b A_bg A_orig exp(h), at the final state. A (squares), C (bounds c_j >= 0), D (spread
                  penalty, no lower bound), D2 and softplus have none. Held templates are not a floor (265 facets).
  fails         : more than 50 % of the seen area switched off

Also: the same with the floor counted, T x 0.1 and x 10, and the share per slope class (without the template and
fine-tuned facets, as in the slope tables). Checks: B's final activity is rebuilt from its terms; with the floor counted
the C shares equal the c5 reports' active seen area.

    q002_active_area.py --out <dir>
"""
import argparse
import csv
import json
from pathlib import Path

import numpy as np

import q002_slope_activity as Q

REPORTED_ACTIVE = {"C_f10_pm0p01_q002_reg": 0.1997, "C_f1_pm0p01_q002_reg": 0.189,
                   "C_f10_pm0p01_tploff_q002_reg": 0.342, "D2_q002_warm": 0.0639}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    out = Path(ap.parse_args().out)
    for n in ("active_area.json", "active_area.csv"):
        if (out / n).exists():
            raise SystemExit("%s exists" % (out / n))
    met = json.load(open(Q.PL))["runs"]
    n, cen, area = Q.elements()
    del cen
    orig = np.array(np.memmap(Q.BASE / "elements.bin", dtype="<f8", mode="r", offset=12, shape=(n, 34))[:, 32])
    with np.load(Path(met[Q.A_RUN]["fit_dir"]) / "nonneg_facets.npz") as z:
        seen, bins = z["seen"].astype(bool), z["bin"].astype(int)
        special = (z["TS"] != 0) | (z["held"] != 0)
    med = float(np.median(orig[seen]))
    T = 0.01 * med
    aseen = area[seen].sum()
    binsel = [seen & ~special & (bins == b) for b in range(7)]
    binarea = [area[s].sum() for s in binsel]
    runs, checks = {}, {}
    for key, v in met.items():
        with np.load(Path(v["ft_dir"]) / "finetune_state.npz") as z:
            A = z["A_final"].astype(float)
            h = z["held_next"].astype(float) if "held_next" in z.files else None
        r = {"data_definition": "final state (finetune_state.npz A_final)", "state": v["ft_dir"]}
        X = A
        if key == "B_q002_reg":
            fit = Path(v["fit_dir"])
            with np.load(fit / "nonneg_facets.npz") as z:
                meta = json.loads(str(z["meta_json"]))
                Afit, u, TS, hfit = z["A"], z["u"], z["TS"], z["held"]
            line = [x.split() for x in (fit / "fit_result.txt").read_text().split("\n") if x.startswith("A_bg")][0]
            Abg = float(line[2])
            m = np.array(meta["multipliers"], float)[bins]
            floor = m * Abg * orig * np.exp(h)
            X = A - floor
            s = 1e-5 * orig
            best = None
            for name, L in (("u + TS", u + TS), ("u", u)):
                g = np.maximum(L, 0.0)
                p = s > 0
                g[p] = s[p] * np.logaddexp(0.0, L[p] / s[p])
                e = np.abs(m * (Abg * orig + g) * np.exp(hfit) - Afit)
                rel = float((e / np.maximum(Afit, T)).max())
                if best is None or rel < best[1]:
                    best = (name, rel)
            checks["B_rebuilt_from_terms"] = {"L": best[0], "max_abs_diff_over_max_A_T": best[1], "A_bg": Abg,
                                              "most_negative_excess_over_T": float(X.min() / T)}
            if best[1] > 1e-6 or X.min() < -1e-6 * T:
                raise SystemExit("B's decomposition does not reproduce its activity: %s" % checks)
        for tag, f in (("T0p1", 0.1), ("T", 1.0), ("T10", 10.0)):
            r["off_%s" % tag] = float(area[seen & (X <= f * T)].sum() / aseen)
            r["off_%s_floor_counted" % tag] = float(area[seen & (A <= f * T)].sum() / aseen)
        r["off_by_bin"] = [float(area[s & (X <= T)].sum() / a) for s, a in zip(binsel, binarea)]
        r["passes"] = r["off_T"] <= 0.5
        r["Q_share_over_F_E"] = v["energy"]["Q_share_over_F_E"]
        r["facets_over_F_E"] = v["energy"]["facets_over_F_E"]
        if key in REPORTED_ACTIVE:
            got = 1.0 - r["off_T_floor_counted"]
            checks[key] = {"active": got, "reported": REPORTED_ACTIVE[key]}
            if abs(got - REPORTED_ACTIVE[key]) > 0.0015:
                raise SystemExit("%s: active %.4f, reported %.4f" % (key, got, REPORTED_ACTIVE[key]))
        runs[key] = r
        print("%-32s off %.4f (floor counted %.4f); x0.1 %.4f x10 %.4f %s" % (
            key, r["off_T"], r["off_T_floor_counted"], r["off_T0p1"], r["off_T10"], "pass" if r["passes"] else "FAIL"))
    rec = {"decision": "1b (user 2026-11-21): active = A without the imposed floor > T; more than 50 % of the seen "
                       "area switched off fails",
           "T": T, "srcCG_median_seen": med, "T_definition": "0.01 median(AH2O of the seen facets, one value per facet)",
           "imposed_floor": "option B only: m_b A_bg A_orig exp(h)", "fail_above": 0.5,
           "off_by_bin": "share of each slope class of the seen area, without the template and fine-tuned facets",
           "checks": checks, "runs": runs}
    (out / "active_area.json").write_text(json.dumps(rec, indent=1, allow_nan=False) + "\n")
    cols = ["off_T", "off_T_floor_counted", "off_T0p1", "off_T10", "passes", "Q_share_over_F_E", "facets_over_F_E"]
    with open(out / "active_area.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["run"] + cols)
        for k, r in runs.items():
            w.writerow([k] + [r[c] for c in cols])
    print("checks", json.dumps(checks))


if __name__ == "__main__":
    main()
