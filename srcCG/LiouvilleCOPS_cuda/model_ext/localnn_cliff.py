#!/usr/bin/env python3
# ==== version-controlled copy: header added 2026-11-22, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/model_ext/localnn_cliff.py (on Voyager: /data001/highresolu_vt/cops_output/model_ext/localnn_cliff.py)
# SHA256 of the run copy = SHA256 of this file without these header lines: a99bcebe9e107b88a621abd275b670ed4aa3a5931492992adeea60a38a87d19d
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: helper of the cliff test: baseline, table, the plan's reading rules
# RUN: jobs/localnn_cliff_*.sbatch
# ==== end of header ====
"""Helper of the cliff test (user 2026-11-22: "prepare the cliff test"; plan figures_cluster/12M/q002_decisions_20261121/
PLAN_cliff.md): option C, templates off, all multipliers held; the cliff classes (>= 60 deg) at 0.01 against the all-held
baseline, with the moderate slopes (40-60 deg) at 0.01 as the control.

  baseline --controls-c RUN      the all-held baseline's fit folder (follow-up c's chain_pick.json)
  states --runroot R             name=state list of the finished fits, for q002_stationarity_check
  summary --runroot R --baseline FIT --out DIR
      table, the cliff classes' share of the seen production and (median over the rows) of the modelled NG / RG signal, and the reading
      rules of the plan
"""
import argparse
import csv
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
R = Path("/data001/highresolu_vt/cops_output")
FITS = ("cliffs_off_from_held", "cliffs_off_from_d2", "moderate_off_from_held", "moderate_off_from_d2")
CLIFF, MODERATE = (3, 4, 5, 6), (1, 2)          # slope-bin numbers: >= 60 deg, 40-60 deg
SAME, NOT_NEEDED, NEEDED, VALID = 1e-3, 0.01, 0.05, 0.01        # PLAN_cliff.md


def last_fit(root, name):
    d = Path(root) / name
    fit = d / "stage2/fit" if (d / "stage2").exists() else d / "stage1/fit"
    if not (fit / "final").exists():
        return None
    seg = (fit / "final").resolve()
    if not (fit / ("direct_check_%s" % seg.name) / "direct_check.json").exists():
        return None
    msg = json.loads((seg / "fit_manifest.json").read_text())["message"]
    if fit.parent.name == "stage1" and "maximum number of function evaluations" in msg:
        return None                                     # the one restart is due
    return fit


def cmd_baseline(a):
    pick = Path(a.controls_c) / "chain_pick.json"
    if not pick.exists():
        print("follow-up c has no chain_pick.json yet", file=sys.stderr)
        return 3
    print(json.loads(pick.read_text())["fit"])
    return 0


def cmd_states(a):
    out = []
    for name in FITS:
        fit = last_fit(a.runroot, name)
        if fit is not None:
            out.append("%s=%s/localnn_state.npz" % (name, (fit / "final").resolve()))
    print(" ".join(out))
    return 0


def class_shares(seg):
    """Share of each slope class of the seen production at 10:33:10; the cliff classes' share of NG / RG (median row)."""
    import localnn_problem as lp
    import q002_slope_activity as geo
    with np.load(seg / "nonneg_facets.npz", allow_pickle=False) as z:
        A, area, seen, bins = np.asarray(z["A"]), np.asarray(z["area"]), np.asarray(z["seen"], bool), np.asarray(z["bin"]).astype(int)
    t = geo.map_times()
    label = [k for k, v in t.items() if v["utc"].endswith("10:33:10")][0]
    F = A * geo.qfac(label, t[label]["rh_m"], len(A))[0] * area
    prod = [float(F[seen & (bins == b)].sum() / F[seen].sum()) for b in range(7)]
    with np.load(seg / "localnn_state.npz", allow_pickle=False) as z:
        c, m, meta = np.asarray(z["c"]), np.asarray(z["m"]), json.loads(str(z["meta_json"]))
    design = lp.load_design(meta["design"]["path"])
    parts = []
    for b in range(7):
        one = np.zeros(7)
        one[b] = m[b]
        parts.append(lp.all_rows_surface(design, c, one, 0.0, False))
    out = {"production": prod}
    for g, j in (("NG", 0), ("RG", 1)):
        P = np.array([p[j] for p in parts])
        out[g] = float(np.median(P[list(CLIFF)].sum(axis=0) / np.maximum(P.sum(axis=0), 1e-300)))
    return out


def cmd_summary(a):
    sys.path.insert(0, a.code)
    sys.path.append(str(HERE))
    import q002_slopetest_report as rep
    out = Path(a.out)
    out.mkdir()
    T = json.loads(rep.AA.read_text())["T"]
    segs = {"baseline": (Path(a.baseline) / "final").resolve()}
    waiting = []
    for name in FITS:
        fit = last_fit(a.runroot, name)
        if fit is None:
            waiting.append(name)
        else:
            segs[name] = (fit / "final").resolve()
    rows = {k: rep.row(k, s, T) for k, s in segs.items()}
    for k, r in rows.items():
        sh = class_shares(segs[k])
        r["cliff_share_of_seen_production"] = sum(sh["production"][b] for b in CLIFF)
        r["moderate_share_of_seen_production"] = sum(sh["production"][b] for b in MODERATE)
        for g in ("NG", "RG"):
            r["cliff_share_of_%s_signal" % g] = sh[g]
        r["data_cost_vs_baseline"] = r["data_cost"] / rows["baseline"]["data_cost"] - 1.0
        print("%-24s obj %.6e data %.6e (%+.2f %%) late %.4f/%.4f off %.1f %% maxF/F_E %.3f cliffs %.2f %% of seen production, "
              "%.2f / %.2f %% of NG / RG %s" % (k, r["objective"], r["data_cost"], 100 * r["data_cost_vs_baseline"],
                                               r["late_rms_NG"], r["late_rms_RG"], 100 * r["off_seen_1b"], r["F_max_over_F_E"],
                                               100 * r["cliff_share_of_seen_production"], 100 * r["cliff_share_of_NG_signal"],
                                               100 * r["cliff_share_of_RG_signal"], r["termination"][:24]))
    verdict = {}
    if not waiting:
        best = {}
        for s in ("cliffs_off", "moderate_off"):
            pair = [rows[s + "_from_held"], rows[s + "_from_d2"]]
            best[s] = min(pair, key=lambda r: r["objective"])
            o = [r["objective"] for r in pair]
            verdict[s + "_starts_agree"] = bool(max(o) / min(o) - 1.0 <= SAME)
        dc, dm = best["cliffs_off"]["data_cost_vs_baseline"], best["moderate_off"]["data_cost_vs_baseline"]
        valid = best["cliffs_off"]["cliff_share_of_seen_production"] < VALID
        if not valid:
            text = "invalid: the held cliffs still give >= 1 % of the seen production"
        elif dc <= NOT_NEEDED:
            text = "the data do not need active cliffs, in this model"
        elif dc >= NEEDED and dc > dm:
            text = "the data need gas from cliffs specifically, in this model"
        elif dc >= NEEDED:
            text = "the data need that much active area; nothing specific to cliffs"
        else:
            text = "undecided"
        verdict.update(delta_cliffs=dc, delta_control=dm, valid=bool(valid), statement=text)
        print("VERDICT (rules of PLAN_cliff.md):", json.dumps(verdict))
    print("waiting:", waiting or "none")
    (out / "cliff_test.json").write_text(json.dumps(dict(rows=rows, verdict=verdict, waiting=waiting), indent=1) + "\n")
    with open(out / "cliff_test.csv", "w", newline="") as f:
        w = csv.DictWriter(f, list(rows["baseline"]))
        w.writeheader()
        w.writerows(rows.values())
    return 0


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("baseline").add_argument("--controls-c", required=True)
    sub.add_parser("states").add_argument("--runroot", required=True)
    p = sub.add_parser("summary")
    for k in ("runroot", "baseline", "out", "code"):
        p.add_argument("--" + k, required=True)
    a = ap.parse_args()
    return {"baseline": cmd_baseline, "states": cmd_states, "summary": cmd_summary}[a.cmd](a)


if __name__ == "__main__":
    sys.exit(main())
