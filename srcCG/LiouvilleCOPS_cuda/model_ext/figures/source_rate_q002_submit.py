#!/usr/bin/env python3
# ==== version-controlled copy: header added 2026-11-18, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/figures_cluster/scripts/source_rate_q002_submit.py (on Voyager: /data001/highresolu_vt/cops_output/figures_cluster/scripts/source_rate_q002_submit.py)
# SHA256 of the run copy = SHA256 of this file without these header lines: cb248fd920cdc5183b27a5a4979a599012381a786d3f06deb8d1aa5707ffe091
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: submits jobs/source_rate_q002.sbatch for one completed new-law run, inputs from its maps manifest and the
#   runs spec (2026-11-18)
# RUN: login node: source_rate_q002_submit.py --run <key> [--node gpu001] [--dry-run]
# CHANGED: 2026-11-23 full set of the all-held candidate: --runs
# ==== end of header ====
"""Submit jobs/source_rate_q002.sbatch for one completed new-law run (2026-11-18; user: "prepare the surface source rate
maps for completed runs for inspection of the source rate"). Every input of the new-law state is taken from its maps'
manifest (base, geometry, fine-tuning and fit folders; the frozen nonneg_maps.py found by its sha256 among
model_ext_frozen/nonneg_maps_*); the old-law counterpart likewise from its maps (runs spec old_maps; null: no
inspection figure, STEPS=v2); the caption texts (model name, activity model) from the runs spec
(source_rate_q002_runs_20261118.json).

    source_rate_q002_submit.py --run <key of the runs spec> [--node gpu001] [--time 00:45:00] [--dry-run]

Refuses an existing output folder. Logs the submission in model_ext/logs/submissions.txt.
"""
import argparse
import datetime
import hashlib
import json
import os
import subprocess
import sys
from pathlib import Path

ROOT = Path("/data001/highresolu_vt/cops_output")
S = ROOT / "figures_cluster" / "scripts"
RUNS = S / "data" / "source_rate_q002_runs_20261118.json"
JOB = ROOT / "jobs" / "source_rate_q002.sbatch"


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def frozen_maps_dir(code_sha):
    hits = [d for d in sorted((ROOT / "model_ext_frozen").glob("nonneg_maps_*"))
            if (d / "nonneg_maps.py").exists() and sha(d / "nonneg_maps.py") == code_sha]
    if not hits:
        raise SystemExit("no model_ext_frozen/nonneg_maps_* with nonneg_maps.py sha256 %s" % code_sha)
    return hits[0]


def state_of_maps(maps):
    man = json.load(open(Path(maps) / "maps_manifest.json"))
    return {"MAPS": str(Path(maps).resolve()), "BASE": man["base"], "GEOM": man["geometry"],
            "FT": man["state"]["ft_dir_resolved"], "FIT": man["state"]["fit_dir_resolved"],
            "MX": str(frozen_maps_dir(man["code_sha256"]["nonneg_maps.py"]))}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", required=True)
    ap.add_argument("--node", default="gpu001")
    ap.add_argument("--time", default="00:45:00")
    ap.add_argument("--dry-run", action="store_true")
    # 2026-11-23 (all-held candidate's full set): another runs spec (a new dated file; the 2026-11-18 one is not edited)
    ap.add_argument("--runs", default=str(RUNS))
    a = ap.parse_args()
    runs = json.load(open(a.runs))
    r = runs["runs"][a.run]
    maps = Path(r["maps"])
    if not (maps / "maps_manifest.json").exists():
        raise SystemExit("maps not complete: %s/maps_manifest.json missing" % maps)
    out = ROOT / "figures_cluster" / "12M" / "source_rate_q002" / a.run
    if out.exists():
        raise SystemExit("%s exists" % out)
    env = {"KEY": a.run}
    env.update(state_of_maps(maps))
    if r.get("old_maps"):
        o = state_of_maps(r["old_maps"])
        env.update({"OLD_" + k: v for k, v in o.items()})
        env["OLD_NAME"] = r["old_name"]
        env["STEPS"] = "v2 inspect"
    else:
        env["STEPS"] = "v2"
    env["MODEL_NAME"] = r["model_name"]
    env["MODEL_TEXT"] = runs["model_texts"][r["model"]] or ""
    for k in ("FT", "FIT"):
        if not (Path(env[k]) / ("finetune_state.npz" if k == "FT" else "fit_restart.npz")).exists():
            raise SystemExit("%s %s incomplete" % (k, env[k]))
    cmd = ["sbatch", "--parsable", "--nodelist=%s" % a.node, "-t", a.time, "--export=ALL", str(JOB)]
    print("\n".join("export %s=%s" % (k, json.dumps(v)) for k, v in env.items()))
    print(" ".join(cmd))
    if a.dry_run:
        return 0
    e = dict(os.environ)
    e.update(env)
    j = subprocess.run(cmd, env=e, cwd=str(ROOT / "jobs"), check=True, capture_output=True, text=True).stdout.strip()
    with open(ROOT / "model_ext" / "logs" / "submissions.txt", "a") as fh:
        fh.write("%s %s source_rate_q002.sbatch gpu (--nodelist=%s, %s) SOURCE-RATE MAPS q002: KEY=%s maps=%s old=%s "
                 "(source_rate_q002_submit.py)\n" % (datetime.datetime.now(datetime.timezone.utc).strftime(
                     "%Y-%m-%dT%H:%M:%SZ"), j, a.node, a.time, a.run, env["MAPS"], env.get("OLD_MAPS", "-")))
    print("submitted %s" % j)
    return 0


if __name__ == "__main__":
    sys.exit(main())
