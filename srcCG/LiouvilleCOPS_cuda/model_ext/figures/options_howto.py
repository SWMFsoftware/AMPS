#!/usr/bin/env python3
# ==== version-controlled copy: header added 2026-11-13, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/figures_cluster/scripts/options_howto.py (on Voyager: /data001/highresolu_vt/cops_output/figures_cluster/scripts/options_howto.py)
# SHA256 of the run copy = SHA256 of this file without these header lines: 86e17329b593a25747f75b5f70915bd20f89c4d352410bb8da00cda75ff60721
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: HOW-TO-RERUN.md and README.txt of a stage-1 options comparison folder (figure commands, each fit's full
#   command and AMPS commit, references) (2026-11-14; v2: --env, --maps-copied-from, --note; 2026-11-15: the ftol line
#   of README.txt from summary.csv)
# RUN: login node: options_howto.py --folder <folder> [--wrapper-log <log>] [--maps-job <job>] [--env ...] [--maps-
#   copied-from <folder>] [--note ...]
# ==== end of header ====
"""HOW-TO-RERUN.md and README.txt of a stage-1 options comparison folder (figures_12m_nonneg_options_stage1.sh;
2026-11-14). Reads <folder>/runs.json and, per plotted fit, the run's fit_manifest.json (full fit command), its
provenance git.txt files (AMPS commit) and its own HOW-TO-RERUN.md (path). Writes <folder>/HOW-TO-RERUN.md and
<folder>/README.txt (refuses to overwrite). Light (JSON only): login node.
Usage: options_howto.py --folder <folder> [--wrapper-log <log>] [--maps-job <jobid>] [--env "<VAR=value ...>"]
       [--maps-copied-from <folder>]
2026-11-14 (v2): --env = the wrapper variables of this folder (DIAG, RESID, ACT_FROM), shown in the A1 command;
--maps-copied-from: activity_maps/ was copied (ACT_FROM) instead of a new GPU job; README.txt lists the residual figure and
the convergence columns when they exist. Without these options the output is unchanged.
"""
import argparse
import csv
import glob
import hashlib
import json
import shlex
import subprocess
import sys
from pathlib import Path

ROOT = Path("/data001/highresolu_vt/cops_output")
AMPS = Path("/data001/highresolu_vt/AMPS")


def first_lines(root):
    out = {}
    for g in sorted(glob.glob(str(root / "**" / "provenance*" / "git.txt"), recursive=True)):
        try:
            c = open(g).readline().strip()
        except OSError:
            continue
        if c:
            out.setdefault(c, []).append(str(Path(g).relative_to(root)))
    return out


def fit_command(final):
    m = json.load(open(final / "fit_manifest.json"))
    me = m.get("model_extension")
    if isinstance(me, dict) and me.get("command"):
        return me["command"], (m.get("compute_backend") or {}).get("driver_command")
    return m.get("command"), None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--folder", required=True)
    ap.add_argument("--wrapper-log", default=None)
    ap.add_argument("--maps-job", default=None)
    ap.add_argument("--env", default="")
    ap.add_argument("--maps-copied-from", default=None)
    ap.add_argument("--note", default=None, help="a line added to README.txt and HOW-TO-RERUN.md (e.g. runs not included)")
    a = ap.parse_args()
    O = Path(a.folder)
    d = json.load(open(O / "runs.json"))
    head = subprocess.run(["git", "-C", str(AMPS), "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip()
    L = ["# How to rerun: SHAP8 12M stage-1 fits of the non-negative activity options vs the data and the bounded softplus c9",
         "",
         "Written by `figures_cluster/scripts/options_howto.py` (2026-11-14). Folder: `%s`." % O,
         "Rules (CLAUDE.md): GPU work only inside Slurm jobs; never write into an existing product directory (every script "
         "here refuses); never delete earlier attempts (move them aside as `*.incomplete_<UTC>` / `*.superseded_<UTC>`).",
         "AMPS tree when the figures were made: `%s` (local changes: `code/git.txt`). The version-controlled copies of the "
         "figure scripts and the job: `srcCG/LiouvilleCOPS_cuda/model_ext/figures/` and `model_ext/jobs/` (SHA256 in "
         "`model_ext/RUN_COPIES.sha256`; working-tree copies of 2026-11-14, not committed yet: the user commits)." % head,
         "",
         "## A. The figures (no fit; CPU, then one short GPU job)",
         "",
         "```bash",
         "export COPS_OUTPUT_ROOT=/data001/highresolu_vt/cops_output",
         "cd $COPS_OUTPUT_ROOT",
         "# A1: tables, model/data figures, captions (CPU; OUT must not exist)",
         "OUT=$COPS_OUTPUT_ROOT/figures_cluster/12M/<new folder> " + (a.env + " " if a.env else "") +
         "srun -p normal128cores -A general -N1 -n1 -c 8 --mem=64G "
         "-t 00:50:00 -J cops_options_cmp bash figures_cluster/scripts/figures_12m_nonneg_options_stage1.sh",
         "# A2: activity maps of the best setting per option + c9 (1 GPU on the gpu partition, about 10 min)",
         "sbatch --export=ALL,OPT_DIR=<the same folder> jobs/options_activity_maps.sbatch",
         "# A3: this file + README.txt, then CAPTIONS.md again (it then includes the activity-map caption)",
         "/data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda/.venv/bin/python figures_cluster/scripts/options_howto.py "
         "--folder <folder> [--wrapper-log <A1 log>] [--maps-job <A2 job>]",
         "/data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda/.venv/bin/python figures_cluster/scripts/collect_captions.py "
         "<folder> --title \"Figure captions: stage-1 fits of the non-negative activity options vs the COPS data and the "
         "bounded softplus c9\" --order \"\"   # after moving the A1 CAPTIONS.md aside",
         "```",
         "",
         "Inside A1 (exact; `figure_plan.sh` holds the resolved model_obs_figures.py calls): for each run "
         "`positive_table.py <run>/fit/final tables/<tag>.dat`; the references `positive_table.py <c9 ft/final> "
         "tables/ref_softplus_c9.dat`, `positive_table.py <D2 final> tables/ref_clip_D2.dat`; `options_table.py --spec "
         "runs.json --out-dir <folder>`; `. figure_plan.sh`; `options_captions.py --out-dir <folder>`; `collect_captions.py`. "
         "The activity-map job's exact command: `provenance/activity_maps-<job>/commands.sh`.",
         "Which runs are included: those of the wrapper's candidate list whose summary (A/B/D: `summary.txt.json`; C: "
         "`analysis/report.json`) and `fit/final/COPS_model.dat` exist. Update when the other fits finish: rerun A1-A3 "
         "into a new folder (or move this one aside as `*.superseded_<UTC>` and rerun into the same name).",
         ""]
    if a.wrapper_log:
        L += ["A1 of this folder: log `%s`." % a.wrapper_log, ""]
    if a.maps_job:
        L += ["A2 of this folder: job %s (`provenance/activity_maps-%s/`)." % (a.maps_job, a.maps_job), ""]
    if a.maps_copied_from:
        L += ["A2 of this folder: none - `activity_maps/` was copied by the wrapper (ACT_FROM) from `%s/activity_maps` "
              "(made there by its A2 job; the best settings, their run roots and inputs are unchanged, checks and SHA256 in "
              "`activity_maps/COPIED_FROM.txt`; the job's provenance copied to `provenance/`)." % a.maps_copied_from, ""]
    if a.note:
        L += [a.note, ""]
    L += ["## B. The plotted fits (each run root has its own HOW-TO-RERUN.md with the submission lines, the summary job and "
          "the from-scratch route)", ""]
    for r in d["runs"]:
        root = Path(r["root"])
        final = root / "fit" / "final"
        cmd, drv = fit_command(final)
        commits = first_lines(root)
        L += ["### %s %s, %s (tag %s)" % (r["option"], r["name"], r["setting"], r["tag"]), "",
              "- run root `%s`; its rerun file `%s`" % (root, root / "HOW-TO-RERUN.md"),
              "- AMPS commit(s) of its jobs: %s" % ("; ".join("`%s` (%s)" % (c, ", ".join(v[:3]) + (" ..." if len(v) > 3 else ""))
                                                               for c, v in commits.items()) or "-"),
              "- final fit command (fit/final/fit_manifest.json):", "", "```bash",
              " ".join(shlex.quote(str(x)) for x in cmd) if cmd else "(not recorded)", "```"]
        if drv:
            L += ["- its GPU driver command (compute_backend.driver_command):", "", "```bash",
                  " ".join(shlex.quote(str(x)) for x in drv), "```"]
        L.append("")
    L += ["## C. The references", "",
          "- bounded softplus c9 (slope multipliers >= 1, 9 fit + fine-tuning cycles): `%s`; full rerun commands in "
          "`../model_vs_observations_softplus_smm1_9cycles/HOW-TO-RERUN.md` (chain jobs 124664 / 124772)." % d["refs"][0]["root"],
          "- clipped D2 (the start of every option fit): `%s`; rerun commands in "
          "`../model_vs_observations_clip_vs_positive/HOW-TO-RERUN.md`." % d["refs"][1]["root"], ""]
    (O / "HOW-TO-RERUN.md").open("x").write("\n".join(L) + "\n")
    R = ["SHAP8 12M: stage-1 fits of the four non-negative activity options vs the COPS data and the bounded softplus c9",
         "(user 2026-11-14: \"make comparison plot of the results with the data and the current softplus fit\").",
         "",
         "summary.txt / summary.csv   the table per fit (data cost, penalty, nfev, rms ln NG/RG overall/before the gap/late, "
         "multipliers, parameters at bounds, Q_H2O at 10:33:10, active seen area, seen area holding 50/90 % of the "
         "production, NG/RG shares from the active area), read from each run's summary.txt.json / analysis/report.json",
         "best.json                   the lowest-data-cost setting of each option",
         "12M_options_<A|B|C|D>_*     per option: all its settings vs the data and c9 (time; altitude; NG/RG panels)",
         "12M_options_best_*          the best setting of each option + c9 + D2",
         "activity_maps/              log A of the best setting of each option and c9 (global + impact close-up)" +
         (" - copied from %s (inputs unchanged; COPIED_FROM.txt)" % a.maps_copied_from if a.maps_copied_from else ""),
         "CAPTIONS.md, *_caption.txt  captions; tables/ the plotted observation rows; runs.json the included runs",
         "figure_plan.sh              the resolved figure commands; code/ the scripts as run; provenance/ job records",
         "HOW-TO-RERUN.md             exact commands (figures, each fit, references)"] + (
        ["12M_options_best_residuals  ln(model / observation) vs altitude, NG and RG, best of each option + c9 + D2"]
        if (O / "12M_options_best_residuals.png").exists() else []) + (
        ["summary convergence column  verdict per fit from the CPU diagnosis checks/stage1_diag_20261114 (code/"
         "stage1_convergence_20261114.json; numbers from its diag_*.json)"]
        if "convergence" in open(O / "summary.csv").readline() else []) + [
         "",
         "Included runs:"] + ["  %s %s %s  %s" % (r["option"], r["name"], r["setting"], r["root"]) for r in d["runs"]] + [
         "", "Reaching max_nfev is not convergence; stopped on ftol (evaluations): %s (see the captions)." % (
             ", ".join("%s %s (%s)" % (r["option"], r["setting"], r["nfev"]) for r in csv.DictReader(open(O / "summary.csv"))
                       if r["option"] != "ref" and "ftol" in (r.get("termination") or "")) or "none")] + (
        ["", a.note] if a.note else [])
    (O / "README.txt").open("x").write("\n".join(R) + "\n")
    print("written %s, %s" % (O / "HOW-TO-RERUN.md", O / "README.txt"))
    return 0


if __name__ == "__main__":
    sys.exit(main())
