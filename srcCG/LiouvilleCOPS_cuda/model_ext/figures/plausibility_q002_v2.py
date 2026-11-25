#!/usr/bin/env python3
# ==== version-controlled copy: header added 2026-11-23, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/figures_cluster/scripts/plausibility_q002_v2.py (on Voyager: /data001/highresolu_vt/cops_output/figures_cluster/scripts/plausibility_q002_v2.py)
# SHA256 of the run copy = SHA256 of this file without these header lines: 3b98197759505b096ab09f20ddea7f9a919ee5632723d1fe4b5c1a6e2b135ba0
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: plausibility numbers of every new-law fit with one method (2026-11-23)
# RUN: jobs/plausibility_q002_v2.sbatch
# ==== end of header ====
"""Physical plausibility of every new-law (q_min = 0.02 q_max) fit, one method for all (2026-11-23; user: "update the
analysys physical plasuability of each of the fit and range it").

Fits: the ten of the round (q002_plausibility), the four other solutions of option C templates off (FINDINGS A22), the
all-held candidate's final state (A23) and the two cliff-test variants (A24, test fits, not candidate models).

From each fit's final per-facet activity A and the illumination law at the 13 map times of 2016-09-30 (the geometry of
the maps; q002_slope_activity.qfac), with one code for all:
  Q_H2O(t) = sum F a, F = A q;  share of Q from seen terrain; max / min of Q over the 13 times
  energy limit F_E = S0 / (r_h^2 dH): share of Q(10:33:10) from facets with F > F_E; max F / F_E over the 13 times
  switched off (user decision 1b): seen area with A <= T (T = 1 % of the srcCG median of the seen facets); option B's
  imposed floor is not counted (its numbers are taken from q002_decisions_20261121/active_area.json)
plus, from the run files: data cost, rms ln(model/data) overall and in the late window, last-fit termination.
Check: Q(10:33:10) and the share above F_E of the ten agree with q002_plausibility/metrics.json (from the maps).

    plausibility_q002_v2.py --out <new folder>
"""
import argparse
import csv
import json
import sys
from pathlib import Path

import numpy as np

import q002_slope_activity as geo

R = Path("/data001/highresolu_vt/cops_output")
F12 = R / "figures_cluster/12M"
L = R / "runs/shap8_012m_lategrid_dz1.5/nonneg/localnn"
MS = L / "multistart_q002_f10_pm0.01_tploff_20261122T065041Z"
CC = L / "slopetest_c_q002_f10_pm0.01_tploff_20261122T195607Z"
CT = L / "clifftest_q002_f10_tploff_20261122T233834Z"
S0, DH, NA = 1361.0, 51.0e3, 6.02214076e23
AU = 1.495978707e11
TEN = {"D2_q002_warm": "D2_q002", "softplus_q002_reg": "softplus_q002", "A_q002_reg": "A_q002", "B_q002_reg": "B_q002",
       "D_q002_reg": "D_q002", "C_f1_pm1em4_q002_reg": "C_f1_pm1em4_q002", "C_f10_pm1em4_q002_reg": "C_f10_pm1em4_q002",
       "C_f1_pm0p01_q002_reg": "C_f1_pm0p01_q002", "C_f10_pm0p01_q002_reg": "C_f10_pm0p01_q002",
       "C_f10_pm0p01_tploff_q002_reg": "C_f10_pm0p01_tploff_q002"}
LABEL = {"D2_q002_warm": "D2 (clipped harmonics)", "softplus_q002_reg": "softplus", "A_q002_reg": "A (sum of squares)",
         "B_q002_reg": "B (background)", "D_q002_reg": "D (log spread)", "C_f1_pm1em4_q002_reg": "C f1, P_m 1e-4",
         "C_f10_pm1em4_q002_reg": "C f10, P_m 1e-4", "C_f1_pm0p01_q002_reg": "C f1, P_m 0.01",
         "C_f10_pm0p01_q002_reg": "C f10, P_m 0.01", "C_f10_pm0p01_tploff_q002_reg": "C f10, templates off: 40-50 deg solution"}
NEW = [("C_off_50_60", "C f10, templates off: 50-60 deg solution", MS / "from_only_50_60/fit/final", "fit", "other solution (A22)"),
       ("C_off_60_70", "C f10, templates off: 60-70 deg solution", MS / "from_only_60_70/fit/final", "fit", "other solution (A22)"),
       ("C_off_70_80", "C f10, templates off: 70-80 deg solution", MS / "from_only_70_80/fit/final", "fit", "other solution (A22)"),
       ("C_off_80_90", "C f10, templates off: 80-90 deg solution", MS / "from_only_80_90/fit/final", "fit", "other solution (A22)"),
       ("C_allheld", "C f10, templates off, all multipliers 1 (candidate)", CC / "chain/c5/ft/final", "ft", "candidate (A23)"),
       ("C_allheld_cliffs_off", "candidate with cliffs at 0.01 (test)", CT / "cliffs_off_from_held/stage1/fit/final", "fit", "test fit (A24)"),
       ("C_allheld_moderate_off", "candidate with 40-60 deg at 0.01 (control)", CT / "moderate_off_from_held/stage1/fit/final", "fit", "test fit (A24)")]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    out = Path(ap.parse_args().out)
    out.mkdir()
    met = json.load(open(geo.PL))["runs"]
    old = json.load(open(F12 / "q002_plausibility/metrics.json"))["runs"]
    act = json.load(open(F12 / "q002_decisions_20261121/active_area.json"))
    T = act["T"]
    summ = {r["key"]: r for r in csv.DictReader(open(F12 / "model_vs_observations_nonneg_options_stage2_q002/summary.csv"))}
    n, cen, area = geo.elements()
    del cen
    with np.load(Path(met[geo.A_RUN]["fit_dir"]) / "nonneg_facets.npz") as z:
        seen = z["seen"].astype(bool)
    times = geo.map_times()
    labels = sorted(times, key=lambda k: times[k]["utc"])
    t0 = [k for k in labels if times[k]["utc"].endswith("10:33:10")][0]
    q = {k: geo.qfac(k, times[k]["rh_m"], n)[0] for k in labels}
    FE = {k: S0 / ((times[k]["rh_m"] / AU) ** 2 * DH / NA) for k in labels}

    states = []
    for key, s in TEN.items():
        with np.load(Path(met[key]["ft_dir"]) / "finetune_state.npz") as z:
            A = np.asarray(z["A_final"], float)
        r = summ[s]
        states.append(dict(key=key, label=LABEL[key], group="round of 2026-11-18/19", A=A,
                           data_cost=float(r["data_cost"]), rms_NG=float(r["rms_NG_all"]), rms_RG=float(r["rms_RG_all"]),
                           late_NG=float(r["rms_NG_late"]) if r["rms_NG_late"] else None,
                           late_RG=float(r["rms_RG_late"]) if r["rms_RG_late"] else None,
                           termination=r["fit_termination"], off=act["runs"][key]["off_T"],
                           off_T10=act["runs"][key]["off_T10"], templates="off" if "tploff" in key else "on (old-law)"))
    for key, lab, d, kind, group in NEW:
        d = d.resolve()
        if kind == "ft":
            with np.load(d / "finetune_state.npz") as z:
                A = np.asarray(z["A_final"], float)
            c = json.load(open(d / "fine_tuning.json"))["costs_final"]
            man = json.load(open((d.parent.parent / "fit/final").resolve() / "fit_manifest.json"))
            dc, rn, rr, ln, lr = c["cost_data"], c["rms_ln_NG"], c["rms_ln_RG"], c["late_window_rms_ln_NG"], c["late_window_rms_ln_RG"]
        else:
            with np.load(d / "nonneg_facets.npz") as z:
                A = np.asarray(z["A"], float)
            man = json.load(open(d / "fit_manifest.json"))
            cc = man["costs"]
            dc, rn, rr, ln, lr = cc["data_cost"], cc["rms_ln_NG"], cc["rms_ln_RG"], cc["late_rms_ln_NG"], cc["late_rms_ln_RG"]
        states.append(dict(key=key, label=lab, group=group, A=A, data_cost=dc, rms_NG=rn, rms_RG=rr, late_NG=ln, late_RG=lr,
                           termination=man["message"], off=float(area[seen & (A <= T)].sum() / area[seen].sum()),
                           off_T10=float(area[seen & (A <= 10 * T)].sum() / area[seen].sum()), templates="off"))
    rows, checks = [], {}
    for s in states:
        A = s.pop("A")
        Q = {k: float((A * q[k] * area).sum()) for k in labels}
        F0 = A * q[t0]
        Fa = F0 * area
        s.update(Q_103310=Q[t0], Q_seen_share=float(Fa[seen].sum() / Fa.sum()),
                 Q_share_over_F_E=float(Fa[F0 > FE[t0]].sum() / Fa.sum()),
                 max_F_over_F_E=max(float((A * q[k]).max() / FE[k]) for k in labels),
                 Q_max_over_min=max(Q.values()) / min(Q.values()), passes_1b=s["off"] <= 0.5)
        if s["key"] in old:
            e = old[s["key"]]["energy"]
            checks[s["key"]] = dict(Q_rel=s["Q_103310"] / e["Q_H2O_s-1"] - 1.0,
                                    over_F_E_diff=s["Q_share_over_F_E"] - e["Q_share_over_F_E"])
        rows.append(s)
        print("%-30s data %.4e late %s off %.1f %% Q %.3e seen %.1f %% >F_E %.2f %% maxF/F_E %.3g Qmax/min %.2f" % (
            s["key"], s["data_cost"], "-" if s["late_NG"] is None else "%.3f/%.3f" % (s["late_NG"], s["late_RG"]),
            100 * s["off"], s["Q_103310"], 100 * s["Q_seen_share"], 100 * s["Q_share_over_F_E"], s["max_F_over_F_E"],
            s["Q_max_over_min"]), flush=True)
    worst = max(max(abs(c["Q_rel"]) for c in checks.values()), max(abs(c["over_F_E_diff"]) for c in checks.values()))
    print("check against q002_plausibility/metrics.json (from the maps): worst |dQ/Q| or |d share| %.2e" % worst)
    (out / "plausibility_v2.json").write_text(json.dumps(dict(T=T, map_times=[times[k]["utc"] for k in labels],
                                                              rows=rows, check_vs_maps=checks), indent=1) + "\n")
    with open(out / "plausibility_v2.csv", "w", newline="") as f:
        w = csv.DictWriter(f, list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    return 0 if worst < 1e-6 else 1


if __name__ == "__main__":
    sys.exit(main())
