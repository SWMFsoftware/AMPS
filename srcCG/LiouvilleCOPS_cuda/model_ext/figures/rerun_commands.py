#!/usr/bin/env python3
# ==== version-controlled copy: header added 2026-11-06, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/figures_cluster/scripts/rerun_commands.py (on Voyager: /data001/highresolu_vt/cops_output/figures_cluster/scripts/rerun_commands.py)
# SHA256 of the run copy = SHA256 of this file without these header lines: e3a788d5b4dbeaec23cfc521e96b041f9b348f151a5af788933ec88c64a50181
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: writes the "exact commands" part of HOW-TO-RERUN.md for one grain fit + fine tuning: the AMPS commit of the run,
#   for every model_ext module and job script used the git commit holding an identical version-controlled copy, the
#   job's exports, and every command with all arguments resolved (the job's commands.sh if it wrote one; otherwise
#   rebuilt from fit_manifest.json, the restart history and fine_tuning.json). Read only. (user request 2026-11-06)
# RUN (from $APP): .venv/bin/python <this> --fit <fit dir> [--finetune <dir>] [--provenance <job dir>]
#        [--restart-start <npz>] [--template NAME=<fine-tuning dir>,<fit_restart.npz>] [--label <text>]
# ==== end of header ====
"""Exact rerun commands and code versions of one grain fit (+ its fine tuning), as Markdown for HOW-TO-RERUN.md
(user request 2026-11-06: "HOW-TO-RERUN.md should contain the exact commands and all arguments that are needed to do
a rerun as the running scripts might change; also it needs to contain the commit id of the code itself").

    .venv/bin/python rerun_commands.py --fit <fit dir> [--finetune <grain fine-tuning dir>] [--provenance <job dir>]
        [--restart-start <npz made by restart_grain_start.py>] [--template NAME=<fine-tuning dir>,<its fit_restart.npz>]
        [--label "<case>"]

Prints, for the case:
* code versions: the AMPS commit and tracked changes the job ran with (provenance git.txt); for every model_ext module
  the job used (provenance model_ext.sha256 = the run copies' SHA256): the git commit that holds a version-controlled
  copy with that content (srcCG/LiouvilleCOPS_cuda/model_ext, header block removed), or "not in git";
* the environment (exports) of the job;
* the commands, every argument resolved: templates (ft_template.py), the restart start (restart_grain_start.py, from
  the restart file's history), the fit (fit_manifest.json model_extension.command, verbatim), the fine tuning
  (grain_finetune.py, from fine_tuning.json); if the job wrote <provenance>/commands.sh, that file is printed instead.
Standard library + numpy only; read only.
"""

import argparse
import hashlib
import json
import shlex
import subprocess
from pathlib import Path

import numpy as np

AMPS = Path("/data001/highresolu_vt/AMPS")
GITDIR = "srcCG/LiouvilleCOPS_cuda/model_ext"
APP = AMPS / "srcCG/LiouvilleCOPS_cuda"
BEGIN = "# ==== version-controlled copy: header added"
END = "# ==== end of header ===="


def strip_header(text):
    out, skip = [], False
    for line in text.split("\n"):
        if line.startswith(BEGIN):
            skip = True
            continue
        if skip and line.startswith(END):
            skip = False
            continue
        if not skip:
            out.append(line)
    return "\n".join(out)


def git(*args):
    return subprocess.run(["git", "-C", str(AMPS)] + list(args), capture_output=True, text=True).stdout


def commit_of(name, sha):
    """newest commit whose version-controlled copy of <name> equals the run copy with this SHA256, else None"""
    for sub in ("", "jobs/", "figures/"):
        path = "%s/%s%s" % (GITDIR, sub, name)
        for c in git("log", "--format=%h", "--", path).split():
            blob = git("show", "%s:%s" % (c, path))
            if blob and hashlib.sha256(strip_header(blob).encode()).hexdigest() == sha:
                return c
    return None


def q(args):
    return " ".join(shlex.quote(str(a)) for a in args)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fit", required=True, type=Path)
    ap.add_argument("--finetune", type=Path)
    ap.add_argument("--provenance", type=Path)
    ap.add_argument("--restart-start", type=Path)
    ap.add_argument("--template", action="append", default=[])
    ap.add_argument("--label", default="")
    a = ap.parse_args()
    m = json.loads((a.fit / "fit_manifest.json").read_text())
    me = m.get("model_extension", {})
    out = ["### %s" % (a.label or a.fit.name), ""]
    # ---- code versions
    if a.provenance and (a.provenance / "git.txt").exists():
        g = (a.provenance / "git.txt").read_text().split("\n")
        changes = [l.strip() for l in g[1:] if l.strip() and not l.startswith("??")]
        out += ["* AMPS commit of the run: `%s`; tracked local changes: %s" % (g[0].strip(), ", ".join(changes) or "none")]
    if a.provenance and (a.provenance / "model_ext.sha256").exists():
        used = set()
        for p in a.provenance.glob("*.py"):
            used.add(p.name)
        for p in a.provenance.glob("*.sbatch"):
            used.add(p.name)
        rows = []
        for line in (a.provenance / "model_ext.sha256").read_text().split("\n"):
            if not line.strip():
                continue
            sha, path = line.split()[:2]
            name = Path(path).name
            if name in used:
                c = commit_of(name, sha)
                rows.append("  | `%s` | `%s` | %s |" % (name, sha[:16], "`%s`" % c if c else "**not in git**"))
        for p in sorted(a.provenance.glob("*.sbatch")):
            sha = hashlib.sha256(p.read_bytes()).hexdigest()
            c = commit_of(p.name, sha)
            rows.append("  | `%s` (job script, provenance copy) | `%s` | %s |" % (p.name, sha[:16], "`%s`" % c if c else "**not in git**"))
        out += ["* Code that ran (run-copy SHA256; git commit holding an identical version-controlled copy in `%s`):" % GITDIR,
                "", "  | file | SHA256 (first 16) | commit |", "  |---|---|---|"] + rows + [""]
    # ---- commands
    out += ["Commands (every argument as run; from `%s`):" % (a.provenance or a.fit), "", "```bash"]
    if a.provenance and (a.provenance / "commands.sh").exists():
        out += [(a.provenance / "commands.sh").read_text().rstrip()]
    else:
        env = []
        # the job script's export lines (static; the provenance environment.txt of jobs before 2026-11-06 23:50 UTC
        # missed CUDA_*/COPS_*/CUPY_* because of a pattern error), in order
        for p in sorted(a.provenance.glob("*.sbatch")) if a.provenance else []:
            env += [l.strip() for l in p.read_text().split("\n")
                    if l.startswith("export ") and "COPS_CUDA_RECORD_STORE" not in l]   # set below, resolved
            break
        out += env + ["cd %s" % APP, 'PY=".venv/bin/python -W ignore"']
        for t in a.template:
            name, rest = t.split("=", 1)
            ftd, rst = rest.split(",")
            tpl_out = Path(m and next((x for x in me.get("command", []) if str(x).endswith("/%s.templates.npz" % name)), "")).parent
            out += ["$PY %s" % q([str(Path(me["driver"]).parent / "ft_template.py"), "--finetune", ftd, "--source-restart", rst,
                                  "--target-restart", rst, "--out", tpl_out, "--name", name])]
        if a.restart_start:
            z = np.load(a.restart_start, allow_pickle=True)
            hist = [h for h in json.loads(str(z["meta_json"])).get("history", []) if h.get("action") == "grain start from the fine tuning"]
            if hist:
                h = hist[-1]
                out += ["$PY %s" % q([str(Path(me["driver"]).parent / "restart_grain_start.py"), "--restart", h["source_restart"]["path"],
                                      "--finetune", h["finetune"], "--out", a.restart_start])]
        cmd = me.get("command") or m["compute_backend"]["driver_command"]
        out += ["bash -c 'echo 1000 > /proc/self/oom_score_adj; exec \"$@\"' _ $PY %s" % q(cmd)]
        if a.finetune:
            f = json.loads((a.finetune / "fine_tuning.json").read_text())
            s = f["settings"]
            method = "cpu" if "--method cpu" in f.get("implementation", "") else "?"
            prec = "mixed" if f.get("record_precision") == "fp32" else "fp64"
            store = m["compute_backend"].get("record_store", {})
            store = store.get("path") if isinstance(store, dict) else store
            if store:
                out += ["export COPS_CUDA_RECORD_STORE=%s" % store]
            out += ["$PY %s" % q([str(Path(me["driver"]).parent / "grain_finetune.py"), f["base_run"], "--restart-from",
                                  f["source_restart"]["path"], "--minutes", "%g" % s["minutes"], "--intervals", s["intervals"],
                                  "--method", method, "--max-nfev", s["max_nfev"], "--devices",
                                  ",".join(str(d) for d in f["devices"]), "--precision", prec, "--out", a.finetune])]
    out += ["```", ""]
    print("\n".join(out))


if __name__ == "__main__":
    main()
