#!/usr/bin/env python3
# ==== version-controlled copy: header added 2026-11-13, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/figures_cluster/scripts/stage2_howto.py (on Voyager: /data001/highresolu_vt/cops_output/figures_cluster/scripts/stage2_howto.py)
# SHA256 of the run copy = SHA256 of this file without these header lines: 2ad6709900dadf5fe2aaf80a053d8e872cf7fa09d226d72fc205ade88471c32f
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: HOW-TO-RERUN.md and README.txt of the final stage-2 comparison folder (2026-11-17)
# RUN: login node: stage2_howto.py --folder <folder> --spec <spec> --log <build log> [--note ...]
# ==== end of header ====
"""HOW-TO-RERUN.md and README.txt of the final stage-2 comparison folder (stage2_figures.py; 2026-11-17).
Reads <folder>/included.json, summary.csv and code/<spec>.json; per chain: its root, its own HOW-TO-RERUN.md / README.md,
the AMPS commit(s) of its jobs (first line of every provenance*/git.txt below the chain root, depth <= 3). Light (login).
    stage2_howto.py --folder <folder> --spec <spec.json> --log <build log> [--note "..."]
Refuses to overwrite.
"""
import argparse
import csv
import glob
import json
import subprocess
import sys
from pathlib import Path

NN = Path("/data001/highresolu_vt/cops_output/runs/shap8_012m_lategrid_dz1.5/nonneg")


def commits(root):
    out = {}
    for pat in ("provenance*/git.txt", "*/provenance*/git.txt", "*/*/provenance*/git.txt"):
        for g in glob.glob(str(root / pat)):
            try:
                c = open(g).readline().strip()
            except OSError:
                continue
            if c:
                out.setdefault(c, 0)
                out[c] += 1
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--folder", required=True)
    ap.add_argument("--spec", required=True)
    ap.add_argument("--log", required=True)
    ap.add_argument("--note", default=None)
    a = ap.parse_args()
    O = Path(a.folder)
    inc = json.load(open(O / "included.json"))
    spec = json.load(open(a.spec))
    rows = {r["key"]: r for r in csv.DictReader(open(O / "summary.csv"))}
    head = subprocess.run(["git", "-C", "/data001/highresolu_vt/AMPS", "rev-parse", "HEAD"], capture_output=True,
                          text=True).stdout.strip()
    L = ["# How to rerun: SHAP8 12M final stage-2 comparison of the non-negative activity models", "",
         "Folder `%s`; made by `figures_cluster/scripts/stage2_figures.py` (code as run in `code/`, SHA256SUMS; AMPS tree "
         "`%s`, local changes in `code/git.txt`; version-controlled copies in `srcCG/LiouvilleCOPS_cuda/model_ext/figures/`, "
         "SHA256 in `model_ext/RUN_COPIES.sha256`)." % (O, head),
         "Rules (CLAUDE.md): CPU job, no GPU; never write into an existing folder (the script refuses); move earlier "
         "versions aside as `*.superseded_<UTC>`; the script stops below 25 GB free (shared disk).", "",
         "## A. The figures (CPU, about 1 min)", "", "```bash", "export COPS_OUTPUT_ROOT=/data001/highresolu_vt/cops_output",
         "cd $COPS_OUTPUT_ROOT",
         "srun -p normal128cores -A general -N1 -n1 -c 8 --mem=64G -t 00:50:00 -J cops_stage2_final bash -c \"export "
         "OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1 MPLBACKEND=Agg; "
         "/data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda/.venv/bin/python -W ignore "
         "figures_cluster/scripts/stage2_figures.py --spec figures_cluster/scripts/data/%s --out <new folder>\"" %
         Path(a.spec).name,
         "/data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda/.venv/bin/python figures_cluster/scripts/stage2_howto.py "
         "--folder <new folder> --spec figures_cluster/scripts/data/%s --log <build log>" % Path(a.spec).name, "```", "",
         "Exact figure commands: `figure_plan.sh`; build log of this folder: `%s` (and `build.log`)." % a.log,
         "A chain enters when its last cycle's fine tuning is complete (and, for option C, its `analysis/report.json` "
         "exists); the others are listed in `PENDING.txt` / summary.txt. To add them: rerun A into a new folder, or move "
         "this one aside and rerun into the same name.", ""]
    if a.note:
        L += [a.note, ""]
    L += ["## B. The chains (each root has its own README.md and HOW-TO-RERUN.md with the submission lines and the commands "
          "of every stage: `c<k>/{fit,ft}/provenance-<job>-c<k>{fit,ft}/commands.sh` for A, B, D and softplus; "
          "`provenance-<job>*/` for C)", "",
          "| key | included | final data cost | root | AMPS commit(s) of its jobs (count) |", "|---|---|---|---|---|"]
    for c in spec["chains"]:
        root = NN / c["root"]
        r = rows.get(c["key"])
        cm = commits(root)
        L.append("| %s | %s | %s | `%s` | %s |" % (
            c["key"], "yes" if r else "pending", ("%.4e" % float(r["data_cost"])) if r else "-", root,
            "; ".join("`%s` (%d)" % (k[:12], v) for k, v in cm.items()) or "-"))
    L += ["", "Final state of each chain: `c<N>/ft/final` (N = 5; 9 for softplus_unreg_c9); its rerun: the chain's "
          "HOW-TO-RERUN.md. D2: `/data001/highresolu_vt/cops_output/runs/shap8_012m_lategrid_dz1.5/clip_refit/D2/final` "
          "(rerun commands in `../model_vs_observations_clip_vs_positive/HOW-TO-RERUN.md`).", ""]
    (O / "HOW-TO-RERUN.md").open("x").write("\n".join(L) + "\n")
    pend = inc.get("pending", [])
    R = ["SHAP8 12M: final stage-2 comparison of the non-negative activity models (lead 2026-11-17 12:30 CST)", "",
         "12M_stage2_reg_best_*           regularized chains: softplus reference, D2, A, B, best C, D vs the data",
         "12M_stage2_C_reg_*              all completed regularized option-C chains + softplus reference",
         "12M_stage2_C_unreg_*            unregularized option-C chains + best regularized C + softplus reference",
         "12M_stage2_<A|B|D>_*            regularized vs unregularized chain of the option + softplus reference",
         "12M_stage2_softplus_*           softplus regularized vs unregularized (5 and 9 cycles)",
         "12M_stage2_{reg,unreg}_best_residuals   ln(model / observation) vs altitude, NG and RG",
         "  (each model/data set: _time, _altitude, _altitude_NG, _altitude_RG, _altitude_NG_RG; .eps + .png 300 dpi)",
         "summary.txt / summary.csv       final data cost, P_delta, terminations, rms, multipliers, bounds, cost history",
         "CAPTIONS.md, *_caption.txt      captions; tables/ the plotted observation rows; figure_plan.sh; code/; build.log",
         "HOW-TO-RERUN.md                 commands; the chains and their AMPS commits", "",
         "Included chains: " + ", ".join(inc["rows"]),
         "Pending (not included): " + ("; ".join("%s (%s; %s)" % (p["key"], p["why"], p["note"]) for p in pend) or "none"),
         "", "No activity maps (disk: shared, small products only). Reaching max_nfev is not convergence."]
    if a.note:
        R += ["", a.note]
    (O / "README.txt").open("x").write("\n".join(R) + "\n")
    print("written %s, %s" % (O / "HOW-TO-RERUN.md", O / "README.txt"))
    return 0


if __name__ == "__main__":
    sys.exit(main())
