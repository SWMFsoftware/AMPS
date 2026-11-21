#!/usr/bin/env python3
# ==== version-controlled copy: header added 2026-11-19, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/figures_cluster/scripts/stage2_q002_docs.py (on Voyager: /data001/highresolu_vt/cops_output/figures_cluster/scripts/stage2_q002_docs.py)
# SHA256 of the run copy = SHA256 of this file without these header lines: c22ab47888d057a779a9179084f2f2d66005ea6147bf19c884dad921f7ae9f97
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: README.txt and HOW-TO-RERUN.md of the new-law stage-2 comparison folder (2026-11-19)
# RUN: login node after the figure job
# ==== end of header ====
"""README.txt and HOW-TO-RERUN.md of the new-law stage-2 comparison folder (stage2_q002_figures.py; 2026-11-19).
Reads <folder>/included.json, summary.csv, the job's provenance (<folder>.provenance-<job>/) and the chains' roots
(AMPS commits of their jobs via stage2_howto.commits). Light (login node). Refuses to overwrite.

    stage2_q002_docs.py --folder <folder> [--note "..."]
"""
import argparse
import csv
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from stage2_howto import commits  # noqa: E402


def g(v, f):
    try:
        return f % float(v)
    except (TypeError, ValueError):
        return "-"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--folder", required=True)
    ap.add_argument("--note", action="append", default=[])
    a = ap.parse_args()
    O = Path(a.folder).resolve()
    for n in ("README.txt", "HOW-TO-RERUN.md"):
        if (O / n).exists():
            raise SystemExit("%s exists" % (O / n))
    inc = json.load(open(O / "included.json"))
    rows = list(csv.DictReader(open(O / "summary.csv")))
    by = {r["key"]: r for r in rows}
    provs = sorted(O.parent.glob(O.name + ".provenance-*"))
    figs = sorted({p.name.rsplit(".", 1)[0] for p in O.glob("*.png")})
    L = ["SHAP8 12M: comparison set of the NEW-LAW stage-2 round (H2O q_min = 0.02 q_max) vs the data and vs the old law "
         "(q_min = 0.2 q_max)", "",
         "User 2026-11-19 (relayed by the lead, 17:05 CST): \"prepare a full set of comparison plots for the last round of "
         "fitting\". The last round = the refitted D2 (D2_q002_warm) and the nine regularized stage-2 chains on base_q002 "
         "after the H2O night-side fix of 2026-11-18 (q_min = 0.02 q_max, srcCG's value; FINDINGS A18). Same set as "
         "../model_vs_observations_nonneg_options_stage2_final/ (stage2_figures.py), for the new-law states.", "",
         "FIGURES (.eps + .png 300 dpi + _caption.txt; model/data sets: _time, _altitude, _altitude_NG, _altitude_RG, "
         "_altitude_NG_RG; late window 10:33:10-10:39:26 UTC, 297-5.5 m shaded and labelled)",
         "  12M_q002_best_*                  (1) new law: softplus reference, refitted D2, A, B, best C (%s), D vs COPS" %
         inc["best"].get("C"),
         "  12M_q002_C_*                     (2) new law: the five option-C chains + softplus reference",
         "  12M_q002_vs_old_<run>_*          (3) per run new law (solid) vs old law (dashed): softplus, A, B, D, D2, and the "
         "five C settings",
         "  12M_q002_best_residuals          (4) ln(model/obs) vs altitude, NG and RG, the new-law best set",
         "  12M_q002_C_residuals                 the same for the five C chains + softplus",
         "  12M_q002_vs_old_<run>_residuals      new vs old law: softplus, A, B, D, D2, best C (the other C settings: model/data "
         "figures only)",
         "  12M_q002_convergence_<chain>     (5) data cost after every cycle's fit and fine tuning, new law vs old law, per "
         "chain; 12M_q002_convergence_all: the nine new-law chains + the refitted D2 as a line",
         "  summary.txt / summary.csv        (6) per run (new and old law): final data cost, P_delta, terminations and nfev "
         "of the last fit and fine tuning, rms ln NG / RG all / late, multipliers, parameters at bounds, kappa, Q_H2O total "
         "and seen at 10:33:10 UTC, the old-law cost and the change; cost history per cycle",
         "  CAPTIONS.md, tables/ (the plotted observation rows), figure_plan.sh, included.json, build.log, code/", "",
         "Figures: %d (each with a caption that states the law of every line)." % len(figs), ""]
    L.append("RESULTS (summary.csv; data cost = 0.5 r.r of the 934 NG + 1150 RG observations; change = new / old - 1)")
    L.append("  %-26s %-11s %-11s %-8s %-15s %-15s %-10s %-10s %-10s" % (
        "run (new law)", "data cost", "old law", "change", "rms NG all/late", "rms RG all/late", "Q_H2O", "Q seen",
        "Q old law"))
    for r in sorted([r for r in rows if r["law"] == "new"], key=lambda r: float(r["data_cost"])):
        L.append("  %-26s %-11s %-11s %-8s %-15s %-15s %-10s %-10s %-10s" % (
            r["key"], g(r["data_cost"], "%.4e"), g(r["old_data_cost"], "%.4e"), g(100 * float(r["data_cost_change"]), "%+.1f%%")
            if r["data_cost_change"] else "-", "%s/%s" % (g(r["rms_NG_all"], "%.4f"), g(r["rms_NG_late"], "%.4f")),
            "%s/%s" % (g(r["rms_RG_all"], "%.4f"), g(r["rms_RG_late"], "%.4f")), g(r["Q_H2O_1033"], "%.3e"),
            g(r["Q_H2O_seen_1033"], "%.3e"), g(r["old_Q_H2O_1033"], "%.3e")))
    L += ["  Terminations, multipliers, bounds, kappa, the per-cycle history and the Q sources: summary.txt.", ""]
    L += ["NOT COMPARED / LIMITATIONS",
          "  - C f 10 / P_m 1e-4: its old-law chain has no maps, so there is no old-law Q_H2O. Its new-law Q_H2O total comes "
          "from its final maps (c5/maps_ft manifest), and its seen Q_H2O from the C report (production_1033, c5 FIT state, "
          "i.e. before the last fine tuning). For the other four C chains the report and the inspection files agree to "
          "<= 1.3e-4 (included.json checks).",
          "  - D2 (new and old): one fit, no cycles, so no convergence figure (a line in 12M_q002_convergence_all) and no "
          "late-window rms in the summary (its fit_result has the overall rms only).",
          "  - Residual figures new vs old only for softplus, A, B, D, D2 and the best C (one figure per C setting would repeat "
          "the model/data figures).",
          "  - The old-law Q_H2O values are the old-law blocks of the new runs' inspection files "
          "(source_rate_q002/<run>/inspection/checks/inspect_103310.json).", ""]
    if a.note:
        L += ["NOTES"] + ["  " + n for n in a.note] + [""]
    L.append("JOBS")
    for p in provs:
        tail = (p / "job.log").read_text().strip().splitlines()[-1] if (p / "job.log").exists() else "?"
        L.append("  %s: %s" % (p.name.split("-")[-1], tail))
    L.append("Reaching max_nfev is not convergence; compare fits by final cost, termination and active bounds.")
    (O / "README.txt").open("x").write("\n".join(L) + "\n")
    H = ["# How to rerun: the new-law stage-2 comparison set", "",
         "Folder `%s`; made by `figures_cluster/scripts/stage2_q002_figures.py` (spec "
         "`figures_cluster/scripts/data/stage2_q002_spec_20261119.json`; the old-law spec `stage2_final_spec_20261117.json`), "
         "code as run in `code/` (SHA256SUMS) and in the job's provenance. CPU, about 3 min. The script refuses an existing "
         "folder: move this one aside as `*.superseded_<UTC>` first, never delete it." % O, "",
         "```bash", "cd /data001/highresolu_vt/cops_output/jobs",
         "OUT=%s sbatch --export=ALL stage2_q002_figures.sbatch" % O,
         "# then (login node):",
         "/data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda/.venv/bin/python /data001/highresolu_vt/cops_output/figures_cluster/"
         "scripts/stage2_q002_docs.py --folder %s" % O, "```", ""]
    for p in provs:
        H += ["## The exact commands as run (`%s/commands.sh`)" % p.name, "", "```bash", (p / "commands.sh").read_text().rstrip(),
              "```", ""]
    H += ["The figure commands: `figure_plan.sh`.", "", "## The runs (each root has its own README.md / HOW-TO-RERUN.md)", "",
          "| run | law | final data cost | root | AMPS commit(s) of its jobs (count) |", "|---|---|---|---|---|"]
    for r in rows:
        root = Path(r["root"])
        cm = commits(root if r["key"] not in ("D2_q002", "old_D2") else root.parent)
        H.append("| %s | %s | %s | `%s` | %s |" % (r["key"], r["law"], g(r["data_cost"], "%.4e"), r["root"],
                                                   "; ".join("`%s` (%d)" % (c[:12], n) for c, n in cm.items()) or "-"))
    (O / "HOW-TO-RERUN.md").open("x").write("\n".join(H) + "\n")
    print("wrote README.txt, HOW-TO-RERUN.md in %s (%d figures)" % (O, len(figs)))
    return 0


if __name__ == "__main__":
    sys.exit(main())
