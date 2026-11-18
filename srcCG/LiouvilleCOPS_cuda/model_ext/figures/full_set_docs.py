#!/usr/bin/env python3
# ==== version-controlled copy: header added 2026-11-13, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/figures_cluster/scripts/full_set_docs.py (on Voyager: /data001/highresolu_vt/cops_output/figures_cluster/scripts/full_set_docs.py)
# SHA256 of the run copy = SHA256 of this file without these header lines: e11a01fbc1e5db8cec193bada2ade86aeecfd5191b7f359efda37ffa60ed58c9
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: README.txt, HOW-TO-RERUN.md and code/ of one full figure set (2026-11-17)
# RUN: the full sets of the stage-2 states (jobs/full_set_submit.sh, jobs/finish_full_set.sh)
# ==== end of header ====
"""README.txt, HOW-TO-RERUN.md and code/ of ONE full figure set (2026-11-17, full sets of the six cases of
12M_stage2_reg_best; the layout of ../12M/full_set_softplus_smm1_9cycles/). Light (login node).

    full_set_docs.py --spec <full_sets json> --key <set key>
spec: {"sets": {key: {"title", "model_name", "state", "chain" (or null), "ft", "fit", "maps", "contrib", "phase", "xc",
       "geom", "jobs_state": [ids of the maps / contributions / export jobs], "frozen": {...}, "notes": [...],
       "comparisons": [...]}}}. The figure jobs are read from model_ext/logs/submissions.txt ("FULL SET <key>"), the exact
commands from each job's provenance commands.sh; the facts from checks/facts.txt (full_set_facts.py). Refuses to
overwrite README.txt / HOW-TO-RERUN.md; code/ must not exist.
"""
import argparse
import hashlib
import json
import re
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path("/data001/highresolu_vt/cops_output")
F12 = ROOT / "figures_cluster/12M"
S = ROOT / "figures_cluster/scripts"
SCRIPTS = ("pub_maps_v2_figures.py", "pub_maps_figures.py", "pub_source_area_figures.py", "cliff_importance_figures.py",
           "model_obs_figures.py", "positive_table.py", "model_data_captions.py", "collect_captions.py",
           "contrib_phase_data.py", "frame_contrib.py", "pub_contrib_density_figures.py", "regions_figures.py",
           "full_set_facts.py", "full_set_docs.py", "stage2_figures.py", "options_captions.py", "pub_ng_rg_altitude.py")
JOBS = ("full_set_submit.sh", "finish_full_set.sh", "pub_maps_v2.sbatch", "pub_source_areas.sbatch",
        "pub_maps_figures.sbatch", "contrib_phase_data.sbatch", "pub_contrib_density.sbatch", "nonneg_maps.sbatch",
        "nonneg_contrib.sbatch", "localnn_tools.sbatch", "d2_tools.sbatch")


def sha16(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()[:16]


def du(p):
    r = subprocess.run(["du", "-sh", str(p)], capture_output=True, text=True).stdout.split()
    return r[0] if r else "-"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--spec", required=True)
    ap.add_argument("--key", required=True)
    a = ap.parse_args()
    sp = json.load(open(a.spec))["sets"][a.key]
    O = F12 / ("full_set_%s" % a.key)
    for f in ("README.txt", "HOW-TO-RERUN.md"):
        if (O / f).exists():
            raise SystemExit("%s exists" % (O / f))
    if (O / "code").exists():
        raise SystemExit("%s/code exists" % O)
    sub = [l for l in (ROOT / "model_ext/logs/submissions.txt").read_text().splitlines()
           if ("FULL SET %s:" % a.key) in l or ("FULL SET %s " % a.key) in l]
    jobs = []
    for l in sub:
        m = re.match(r"^(\S+) (\d+) (.*)$", l)
        if m:
            jobs.append((m.group(1), m.group(2), m.group(3)))
    provs = sorted(O.glob("provenance-*"))
    facts = (O / "checks" / "facts.txt").read_text().strip() if (O / "checks" / "facts.txt").exists() else "(facts pending)"
    head = subprocess.run(["git", "-C", "/data001/highresolu_vt/AMPS", "rev-parse", "HEAD"], capture_output=True,
                          text=True).stdout.strip()
    prods = [(n, sp[n]) for n in ("maps", "contrib", "phase", "xc") if sp.get(n)]
    # README
    R = [sp["title"], "=" * len(sp["title"]),
         "Lead 2026-11-17 15:35 CST (user: \"yes, free the 44M store and generate the full sets\"; \"I need the full set\"): a FULL "
         "figure set, exactly like ../full_set_softplus_smm1_9cycles/ (same tools, styles and user decisions), no movies, for "
         "the six cases of 12M_stage2_reg_best. This set: %s." % sp["model_name"],
         "Final state: %s" % sp["state"], "",
         "Every figure: .eps + .png (300 dpi) + _caption.txt; all captions in CAPTIONS.md.", "",
         "visibility/        set 1 (tools/pub_figures.py templates + ground track, unchanged tracked tool) on the observation-only",
         "                   base: fit-independent; the ground track shows the Rosetta impact only; CAPTIONS_COPIED.txt (the",
         "                   PNGs are pixel-identical to the c9 set's, checked; two captions copied)",
         "model_data/        set 4: tables/ (this state and the clipped reference D2), 12M_%s_vs_COPS_*, 12M_%s_pressures_*,"
         % (a.key, a.key),
         "                   captions by model_data_captions.py; publication_NG_RG_altitude/ = copy of",
         "                   ../publication_%s/ (pub_ng_rg_altitude.py)" % a.key.replace("stage2_", "stage2_") if a.key != "D2"
         else "                   ../publication_D2/ (pub_ng_rg_altitude.py)",
         "activity/          final H2O activity A (pub_maps_v2_figures.py --quantity activity --low-cut-area 0.02, flat colours)",
         "production_rate/   H2O production rate F_H2O = q_H2O A at 10:33:10 and 10:39:26 UTC (--quantity h2o, LOWCUT 0.02, relief 0.7)",
         "source_areas/      P% source areas of the 95 post-gap observations, unions and occurrence, NG / RG / NG_RG, P50-90",
         "cliffs/            12M_cliff_slope_bins, 12M_cliff_pressure_fraction_{NG,RG,NG_RG}: one colour per slope bin (--groups all)",
         "contribution_density/  contribution per unit surface area as the fit computes it (snapshots 10:33:10, 10:37:31.282,",
         "                   10:39:09.282, 10:39:26 and phase means), NG / RG / side by side, global + close-up",
         "regions/           share of the modelled surface pressure by lobe and by distance shell from the touchdown point",
         "checks/            class tables, field verification, manifests, slope-bin json, source-area tables, caption json,",
         "                   facts.json / facts.txt (full_set_facts.py), regions / captions logs",
         "code/, provenance-<job>/, HOW-TO-RERUN.md", "",
         "Products next to the state (as for c9): " + "; ".join("%s %s (%s)" % (n, p, du(p)) for n, p in prods) + ".",
         "Not made (c9-specific): comparison_vs_* folders. Possible comparisons (source_area_compare.py, A = this state's",
         "contributions): " + "; ".join(sp.get("comparisons", [])) + ".", ""]
    if sp.get("notes"):
        R += ["Notes:"] + ["* " + n for n in sp["notes"]] + [""]
    R += ["Facts of THIS state (computed from the data by full_set_facts.py; checks/facts.json):", facts, ""]
    (O / "README.txt").open("x").write("\n".join(R) + "\n")
    # code/
    cd = O / "code"
    cd.mkdir()
    (cd / "jobs").mkdir()
    for f in SCRIPTS:
        if (S / f).exists():
            shutil.copy2(S / f, cd / f)
    for f in JOBS:
        if (ROOT / "jobs" / f).exists():
            shutil.copy2(ROOT / "jobs" / f, cd / "jobs" / f)
    shutil.copy2(a.spec, cd / Path(a.spec).name)
    (cd / "git.txt").write_text(subprocess.run(["bash", "-c", "git -C /data001/highresolu_vt/AMPS rev-parse HEAD; git -C "
                                                "/data001/highresolu_vt/AMPS status --short"], capture_output=True,
                                               text=True).stdout)
    subprocess.run(["bash", "-c", "cd '%s' && sha256sum $(find . -type f ! -name SHA256SUMS | sort) > SHA256SUMS" % cd], check=True)
    # HOW-TO
    H = ["# How to rerun: %s" % sp["title"], "",
         "Written 2026-11-17 (UTC) by figures_cluster/scripts/full_set_docs.py. Description and facts: README.txt; captions: "
         "CAPTIONS.md. The state itself: %s (its chain's README.md / HOW-TO-RERUN.md)." % sp["state"],
         "AMPS tree: `%s` (local changes: code/git.txt; not used by the jobs). Code that ran: code/ (SHA256SUMS; the jobs' own "
         "copies in provenance-<job>/code and in each product's provenance). Frozen code: %s." % (
             head, "; ".join("%s `%s`" % kv for kv in sp.get("frozen", {}).items())), "",
         "## A. Submit again (new folders; the tools refuse existing ones)", "", "```bash",
         "cd /data001/highresolu_vt/cops_output/jobs"] + sp.get("state_commands", []) + [
         "# the figure jobs (pub_maps_v2 x5, pub_source_areas sources + cliffs, pub_maps_figures set1/set4/pressures, "
         "contrib_phase_data -> pub_contrib_density):",
         sp["submit_command"],
         "# afterwards, CPU (srun on normal128cores): visibility captions, regions/, model_data captions, facts, CAPTIONS.md",
         sp["finish_command"],
         "# README.txt, HOW-TO-RERUN.md, code/ (login node)",
         "/data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda/.venv/bin/python /data001/highresolu_vt/cops_output/figures_cluster/"
         "scripts/full_set_docs.py --spec %s --key %s" % (a.spec, a.key), "```", "",
         "## B. Jobs as submitted (model_ext/logs/submissions.txt)", ""]
    H += ["| UTC | job | what |", "|---|---|---|"] + ["| %s | %s | %s |" % (t, j, w.replace("|", "/")[:300]) for t, j, w in jobs]
    H += ["", "State jobs (maps / contributions / export): " + ", ".join(sp.get("jobs_state", [])), "",
          "## C. The exact commands as run (provenance-<job>/commands.sh)", ""]
    for p in provs:
        c = p / "commands.sh"
        if c.exists():
            H += ["### %s" % p.name, "", "```bash", c.read_text().strip(), "```", ""]
    (O / "HOW-TO-RERUN.md").open("x").write("\n".join(H) + "\n")
    print("written %s/README.txt, HOW-TO-RERUN.md, code/" % O)
    return 0


if __name__ == "__main__":
    sys.exit(main())
