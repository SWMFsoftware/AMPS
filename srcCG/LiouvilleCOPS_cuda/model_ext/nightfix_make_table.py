#!/usr/bin/env python3
# ==== version-controlled copy: header added 2026-11-18, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/tests/nightfix_diagnostic_20261118/code/make_table.py (on Voyager: /data001/highresolu_vt/cops_output/tests/nightfix_diagnostic_20261118/code/make_table.py)
# SHA256 of the run copy = SHA256 of this file without these header lines: ee5eafb60aba2570bb44564d867c6712d8b234ea93ddbadcea8ea666a50a9293
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: table of the night-side fix diagnostic (2026-11-18): results_v2/*.json + the base check -> table.{txt,csv,json};
#   checks v1 = v2 in every common number.
# RUN: cd tests/nightfix_diagnostic_20261118 && $APP/.venv/bin/python code/make_table.py (login node)
# ==== end of header ====
"""H2O night-side fix diagnostic (2026-11-18): collect results_v2/<key>.json (nightfix_diag.py v2) and
results/base_run_consistency.json into table.csv, table.json and table.txt of this folder; also checks that every
number v1 (results/) and v2 (results_v2/) have in common is identical. Standard library + numpy; login node."""
import csv
import json
import sys
from pathlib import Path

D = Path(__file__).resolve().parents[1]
ORDER = [("stage2_softplus_reg", "softplus, regularized stage 2 (c5)"), ("stage2_A_reg", "A sum of squares r3, reg."),
         ("stage2_B_reg", "B background 1e-2, reg."), ("stage2_D_reg", "D log spread f 0.1, reg."),
         ("stage2_C_f1_pm1em4_reg", "C local basis f 1, P_m 1e-4, reg."), ("D2", "D2 clipped refit"),
         ("stage2_C_f10_pm0.01_reg", "C local basis f 10, P_m 0.01, reg."), ("softplus_c9", "softplus c9 (unreg.)")]
LAWS = ("old", "new")


def flat(r):
    row = {"key": r["key"], "model": r["model"], "recorded_data_cost": r.get("recorded_data_cost"),
           "old_cost_vs_recorded_rel": r.get("old_law_cost_vs_recorded_rel"),
           "old_vs_tables_max_rel_NG": r["verification_old_law_vs_tables"]["NG"],
           "old_vs_tables_max_rel_RG": r["verification_old_law_vs_tables"]["RG"]}
    for lab in LAWS:
        L = r["laws"][lab]
        P = L["production_103310"]
        row.update({
            "%s_data_cost" % lab: L["data_cost"],
            "%s_rms_ln_NG" % lab: L["rms_ln_NG"], "%s_rms_ln_RG" % lab: L["rms_ln_RG"],
            "%s_rms_ln_NG_late" % lab: L["rms_ln_NG_late_window"], "%s_rms_ln_RG_late" % lab: L["rms_ln_RG_late_window"],
            "%s_mean_ln_NG" % lab: L["mean_ln_NG"], "%s_mean_ln_RG" % lab: L["mean_ln_RG"],
            "%s_mean_ln_NG_late" % lab: L["mean_ln_NG_late_window"], "%s_mean_ln_RG_late" % lab: L["mean_ln_RG_late_window"],
            "%s_Q_H2O_103310" % lab: P["Q_H2O_s-1"], "%s_Q_H2O_seen_103310" % lab: P["Q_H2O_seen_s-1"],
            "%s_night_shadow_share_Q" % lab: P["night_shadow_share_of_Q_H2O"],
            "%s_night_shadow_share_seen_Q" % lab: P["night_shadow_share_of_seen_production"],
            "%s_late_night_shadow_share_NG" % lab: L["late_window_night_shadow_share_of_surface_mean"]["NG"],
            "%s_late_night_shadow_share_RG" % lab: L["late_window_night_shadow_share_of_surface_mean"]["RG"],
            "%s_all_night_shadow_share_NG" % lab: L["all_used_night_shadow_share_of_surface_mean"]["NG"],
            "%s_all_night_shadow_share_RG" % lab: L["all_used_night_shadow_share_of_surface_mean"]["RG"],
            "%s_uniform_scale_s" % lab: L["uniform_A_scale"]["s"],
            "%s_uniform_scale_cost" % lab: L["uniform_A_scale"]["data_cost_at_s"]})
    row["new_over_old_cost"] = r["new_over_old"]["data_cost"]
    row["new_over_old_Q"] = r["new_over_old"]["Q_H2O"]
    row["new_over_old_Q_seen"] = r["new_over_old"]["Q_H2O_seen"]
    row["slurm_job_id"] = r["slurm_job_id"]
    return row


def main():
    rows, same = [], {}
    for key, label in ORDER:
        r2 = json.load(open(D / "results_v2" / ("%s.json" % key)))
        row = flat(r2)
        row["label"] = label
        rows.append(row)
        p1 = D / "results" / ("%s.json" % key)
        if p1.exists():
            r1 = json.load(open(p1))
            diffs = []
            for lab in LAWS:
                for k, v in r1["laws"][lab].items():
                    if k != "production_103310" and r2["laws"][lab].get(k) != v:
                        diffs.append("%s.%s" % (lab, k))
                for k, v in r1["laws"][lab]["production_103310"].items():
                    if r2["laws"][lab]["production_103310"][k] != v:
                        diffs.append("%s.production.%s" % (lab, k))
            same[key] = {"v1_job": r1["slurm_job_id"], "v2_job": r2["slurm_job_id"], "differing": diffs}
    base = json.load(open(D / "results" / "base_run_consistency.json"))
    cols = list(rows[0].keys())
    with open(D / "table.csv", "x", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader()
        for r in rows:
            w.writerow(r)
    json.dump({"states": rows, "v1_v2_identical": same, "base_run_consistency": {
        k: base["laws"][k] for k in LAWS}, "base_run_job": base["slurm_job_id"]}, open(D / "table.json", "x"), indent=1)
    L = []
    hdr = ("%-34s %-11s %-11s %6s | %-15s %-15s | %-15s %-15s | %-17s %-17s | %-13s %-13s | %s" % (
        "state", "cost old", "cost new", "x", "rms ln NG o/n", "rms ln RG o/n", "late NG o/n", "late RG o/n",
        "Q_H2O old/new", "Q seen old/new", "dark %Q o/n", "dark %seen o/n", "uniform s (new): cost"))
    L.append(hdr)
    L.append("-" * len(hdr))
    for r in rows:
        L.append("%-34s %-11.4e %-11.4e %6.2f | %.4f / %.4f | %.4f / %.4f | %.4f / %.4f | %.4f / %.4f | "
                 "%.2e / %.2e | %.2e / %.2e | %4.1f / %4.1f | %4.1f / %4.1f | %.3f: %.4e" % (
                     r["label"], r["old_data_cost"], r["new_data_cost"], r["new_over_old_cost"],
                     r["old_rms_ln_NG"], r["new_rms_ln_NG"], r["old_rms_ln_RG"], r["new_rms_ln_RG"],
                     r["old_rms_ln_NG_late"], r["new_rms_ln_NG_late"], r["old_rms_ln_RG_late"], r["new_rms_ln_RG_late"],
                     r["old_Q_H2O_103310"], r["new_Q_H2O_103310"], r["old_Q_H2O_seen_103310"],
                     r["new_Q_H2O_seen_103310"], 100 * r["old_night_shadow_share_Q"], 100 * r["new_night_shadow_share_Q"],
                     100 * r["old_night_shadow_share_seen_Q"], 100 * r["new_night_shadow_share_seen_Q"],
                     r["new_uniform_scale_s"], r["new_uniform_scale_cost"]))
    L.append("")
    L.append("late-window share of the modelled surface pressure from night / cast-shadow records (mean over the used "
             "late rows), old / new: " + "; ".join("%s NG %.2f/%.2f RG %.2f/%.2f" % (
                 r["key"], r["old_late_night_shadow_share_NG"], r["new_late_night_shadow_share_NG"],
                 r["old_late_night_shadow_share_RG"], r["new_late_night_shadow_share_RG"]) for r in rows))
    L.append("mean ln(model/obs) new law, all used rows: " + "; ".join("%s NG %+.3f RG %+.3f" % (
        r["key"], r["new_mean_ln_NG"], r["new_mean_ln_RG"]) for r in rows))
    L.append("uniform A scale with the OLD law (should be ~1 for a fit at its optimum): " + "; ".join(
        "%s %.4f" % (r["key"], r["old_uniform_scale_s"]) for r in rows))
    L.append("base consistency (fit_cops check, tolerance 1e-6): old law max rel NG %.3e RG %.3e; new law NG %.3e RG %.3e"
             % (base["laws"]["old"]["max_rel_NG"], base["laws"]["old"]["max_rel_RG"], base["laws"]["new"]["max_rel_NG"],
                base["laws"]["new"]["max_rel_RG"]))
    L.append("v1 (results/) vs v2 (results_v2/) common numbers identical: " + json.dumps(
        {k: (not v["differing"]) for k, v in same.items()}))
    (D / "table.txt").open("x").write("\n".join(L) + "\n")
    print("\n".join(L))
    return 0


if __name__ == "__main__":
    sys.exit(main())
