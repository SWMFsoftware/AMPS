#!/usr/bin/env python3
# ==== version-controlled copy: header added 2026-11-06, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/figures_cluster/scripts/how_to_rerun.py (on Voyager: /data001/highresolu_vt/cops_output/figures_cluster/scripts/how_to_rerun.py)
# SHA256 of the run copy = SHA256 of this file without these header lines: 18f5a829a352761492ba088c6f7f87fa6766e9461267d0825806f4c1cbea5473
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: writes HOW-TO-RERUN.md next to figures from the fits' manifests (restart and from-scratch commands).
# RUN: .venv/bin/python <this> --res <RES> --run <run> --figdir <dir> --figure-command "<cmd>" FIT_DIR ...
# ==== end of header ====
"""Write HOW-TO-RERUN.md next to a set of model-vs-observation figures (user request 2026-11-03).

    how_to_rerun.py --res 20M --run <runner restart dir> --figdir <figure dir> --figure-command "<cmd>" FIT_DIR...

For every fit directory plotted (a fit, or a fine tuning whose fit is its source) the file gives the exact commands
  (A) from the available restart files: the figures again (no GPU), resubmission of the unchanged job script
      (finished fits are skipped), a continuation of the fit from its fit_restart.npz, the fine tuning;
  (B) from scratch: the runner campaign (base, record store, visibility, gravity, template bundle, fits, fine
      tuning), the C16 set, the grain sets and their fine tuning, then the figures.
The commands are taken from what ran: <run>/pipeline.json (runner argv), each fit_manifest.json
(compute_backend.driver_command; model_extension.command for grain fits), the job scripts' provenance copies.
Nothing is run. Standard library only.
"""
import argparse
import hashlib
import json
import shlex
from datetime import datetime, timezone
from pathlib import Path

OUT = Path("/data001/highresolu_vt/cops_output")
AMPS = Path("/data001/highresolu_vt/AMPS")
APP = AMPS / "srcCG/LiouvilleCOPS_cuda"
PY = APP / ".venv/bin/python"

ENV = f"""```bash
export AMPS={AMPS} APP={APP}
export T=/data001/highresolu_vt/cops_cluster_transfer_20261026 COPS_OUTPUT_ROOT={OUT}
export COPS_DATA_ROOT=$T/rosetta_data_end_of_mission_cops
export CUDA_PATH=/usr/local/cuda-12.9 CUDA_HOME=/usr/local/cuda-12.9 CUPY_CACHE_DIR=$COPS_OUTPUT_ROOT/.cupy
export PATH=$CUDA_HOME/bin:/usr/mpi/gcc/openmpi-4.1.7rc1/bin:$PATH
export LD_LIBRARY_PATH=$CUDA_HOME/lib64${{LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}}
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1
cd $APP
```"""


def q(argv):
    return " ".join(shlex.quote(str(a)) for a in argv)


def wrap(cmd, width=118):
    """Break a long shell command at option boundaries with backslash continuations."""
    parts, line, out = shlex.split(cmd), "", []
    toks = []
    for p in parts:   # keep an option together with its values
        if p.startswith("--") or not toks or p == "--":
            toks.append([p])
        else:
            toks[-1].append(p)
    for t in toks:
        s = q(t)
        if line and len(line) + 1 + len(s) > width:
            out.append(line + " \\")
            line = "    " + s
        else:
            line = (line + " " + s) if line else s
    out.append(line)
    return "\n".join(out)


def load(p):
    try:
        return json.loads(Path(p).read_text())
    except (OSError, ValueError):
        return None


def replace_opt(argv, opt, nvals, new):
    """argv without `opt` and its values, with `new` (list) appended before '--' if given."""
    a, i = list(argv), 0
    while i < len(a):
        if a[i] == opt:
            del a[i:i + 1 + nvals]
        else:
            i += 1
    if new:
        k = a.index("--") if "--" in a else len(a)
        a[k:k] = new
    return a


def set_out(argv, out):
    a = list(argv)
    if "--out" in a:
        a[a.index("--out") + 1] = str(out)
    return a


def fit_entry(d, run):
    """Markdown for one plotted directory (a fit or a fine tuning)."""
    d = Path(d)
    lines = [f"### `{d.relative_to(run) if str(d).startswith(str(run)) else d}`", ""]
    ft = load(d / "fine_tuning.json")
    if ft is not None:   # a fine tuning: its source fit and the command that made it
        sr = ft.get("source_restart", "")
        src = Path(sr.get("path", "") if isinstance(sr, dict) else sr).parent
        grain = (d / "grain_finetune.json").exists()
        lines += [f"Fine tuning (`{'grain_finetune.py' if grain else 'fine_tune_gpu.py'}`) of the fit `{src}`; "
                  f"settings {json.dumps(ft.get('settings', {}))}.", ""]
        if grain:
            cmd = [PY, "-W", "ignore", OUT / "model_ext/grain_finetune.py", run / "base", "--restart-from",
                   src / "fit_restart.npz", "--minutes", "5", "--intervals", "5", "--method", "cpu", "--max-nfev",
                   "2000", "--devices", "0,1,2,3,4,5,6,7", "--precision", "mixed", "--out", d]
            lines += ["Redo this fine tuning from the fit's restart file, inside a GPU job "
                      "(`grain_finetune.sbatch` does exactly this, with `COPS_CUDA_RECORD_STORE` set). "
                      "Write to a new `--out`; never over this directory:", "",
                      "```bash", f"export COPS_CUDA_RECORD_STORE={run}/store",
                      wrap(q(set_out(cmd, str(d) + "_rerun"))), "```", ""]
        else:
            log = run / "logs" / f"finetune_{src.parent.name if src.name.startswith('seg') else src.name}.log"
            head = log.read_text(errors="replace").splitlines()[:5] if log.exists() else []
            cmd = next((h[4:] for h in head if h.startswith("=== ") and "fine_tune_gpu.py" in h), "")
            env = next((h[4:].split(" ", 1)[1] for h in head if h.startswith("=== ") and "COPS_CUDA_RECORD_STORE" in h), "")
            if cmd:
                argv = set_out(shlex.split(cmd), str(d) + "_rerun")
                lines += ["Redo this fine tuning from the fit's restart file: the command the runner "
                          f"ran (from `{log}`), with a new `--out`, inside a GPU job:", "",
                          "```bash"] + [f"export {e}" for e in shlex.split(env)] + [wrap(q(argv)), "```", ""]
        return lines + ["From scratch: §B (it rebuilds the fit and this fine tuning).", ""]

    m = load(d / "fit_manifest.json")
    if m is None:
        return lines + ["(no fit_manifest.json: not a finished fit)", ""]
    drv = m.get("compute_backend", {}).get("driver_command") or []
    ext = m.get("model_extension")
    res = load(d / "fit_report.json") or {}
    if ext:   # grain fit: grain_fit.py with its own options, then the gpu_fit options
        argv = [PY, "-W", "ignore"] + shlex.split(ext["command"]) if isinstance(ext.get("command"), str) \
            else [PY, "-W", "ignore"] + list(ext.get("command", []))
        kind = "grain fit (`model_ext/grain_fit.py` around `tools/gpu_fit.py`)"
    else:
        argv = list(drv)
        kind = "standard fit (`tools/gpu_fit.py`)"
    # the runner's segment and deadline options belong to one job; a rerun uses --wall-time slurm instead
    for opt, n in (("--segment", 1), ("--deadline", 1), ("--wall-margin", 1), ("--save-seconds", 1)):
        argv = replace_opt(argv, opt, n, None)
    if "--wall-time" not in argv:
        argv = replace_opt(argv, "--wall-time", 1, ["--wall-time", "slurm"])
    rerun = set_out(argv, str(d).rstrip("/") + "_rerun")
    cont = replace_opt(set_out(argv, str(d).rstrip("/") + "_cont"), "--restart-from", 1,
                       ["--restart-from", str(d / "fit_restart.npz")])
    cont = [a for a in cont if a != "--check-jacobian"]
    lines += [f"{kind}; started {m.get('compute_backend', {}).get('started_utc', '?')}; "
              f"process state `{m.get('states', {}).get('process', '?')}`.", "",
              "Continue from this fit's restart file (same base, templates from the restart). "
              "`--max-nfev` is the number of further evaluations. A continuation is not identical to one longer "
              "run (the trust region and scaling restart):", "",
              "```bash", wrap(q(cont)), "```", "",
              "Redo the same fit from its default start, on the existing base and record store "
              "(the exact command that ran, with a new `--out`):", "",
              "```bash", wrap(q(rerun)), "```", ""]
    return lines


def fit_input_option(pipe):
    """'' if the run's fitting input is unchanged, else '--fit-input <as-run copy>' (the runner keys a restart
    directory to the input's SHA-256; fit_base.in changed on 2026-11-03), or a note if no copy matches."""
    ident = (pipe.get("identity") or {}).get("fit_input") or {}
    sha, path = ident.get("sha256"), ident.get("path")
    def digest(f):
        return hashlib.sha256(Path(f).read_bytes()).hexdigest()
    if not sha or (path and Path(path).is_file() and digest(path) == sha):
        return ""
    for copy in sorted((OUT / "runs").glob("fit_base.in.as_run_*")):
        if digest(copy) == sha:
            return " --fit-input %s" % copy
    return " --fit-input <a copy of the fitting input with SHA-256 %s>" % sha


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--res", required=True)
    ap.add_argument("--run", required=True, type=Path)
    ap.add_argument("--figdir", required=True, type=Path)
    ap.add_argument("--figure-command", required=True)
    ap.add_argument("dirs", nargs="+")
    a = ap.parse_args()
    run = a.run.resolve()
    pipe = load(run / "pipeline.json") or {}
    inv = next((i for i in reversed(pipe.get("invocations", [])) if i.get("command") == "submit"), None)  # the last
    dirs = list(dict.fromkeys(str(Path(d)) for d in a.dirs))
    grain = any("grain" in Path(d).name for d in dirs)
    # reproductions (model_ext/grain_repro.sbatch): <run>/grain/repro/<source fit>_repro_<job>[_m5_i5_cpu]
    repro = sorted({Path(d).name.removesuffix("_m5_i5_cpu").rsplit("_repro_", 1)[0]
                    for d in dirs if "/grain/repro/" in d and "_repro_" in Path(d).name})
    c16 = any("/c16/" in d for d in dirs) or grain
    md = [f"# How to rerun: SHAP8 {a.res} figures in this folder", "",
          f"Written {datetime.now(timezone.utc).strftime('%Y-%m-%d %H:%M UTC')} by "
          f"`figures_cluster/scripts/how_to_rerun.py`, from the records of the runs (no command here was run to "
          f"write it). Run `{run}`.", "",
          "Rules (CLAUDE.md): GPU work only inside Slurm jobs; never write into an existing product directory "
          "(use a new `--out` or a new restart directory); never delete earlier attempts.", "",
          "## Environment (every shell and job)", "", ENV, "",
          "## A. From the available restart files", "",
          "### A0. The figures only (login node, no GPU, reads the finished fits)", "",
          "```bash", a.figure_command, "```", "",
          "### A1. Resume an interrupted job (time limit, node failure)", "",
          "Resubmit the same job script unchanged: finished fits are skipped, an unfinished one is moved aside as "
          "`*.incomplete_<UTC>` and restarted.", "", "```bash"]
    fi = fit_input_option(pipe)
    if not grain:
        if fi:
            md += ["# this run's fitting input differs from the current fit_base.in: give the runner the as-run copy;",
                   "# the generated fit.sbatch / finetune.sbatch call the runner without it, so resume with",
                   f"#   {PY} -W ignore tools/cluster/cops_pipeline.py run --restart-dir {run} --group <fit|finetune>{fi}",
                   "# inside a GPU job"]
        md += [f"{PY} -W ignore tools/cluster/cops_pipeline.py status --restart-dir {run}{fi}",
               f"sbatch {run}/slurm/fit.sbatch          # runner fits (it resumes at the first incomplete stage)",
               f"sbatch {run}/slurm/finetune.sbatch     # runner fine tuning"]
    if c16 and not grain:
        md += [f"sbatch --export=ALL,COPS_RUN={run} {OUT}/jobs/c16.sbatch"]
    if grain:
        md += [f"sbatch --export=ALL,COPS_RUN={run},GRAIN_SET=free {OUT}/model_ext/grain.sbatch",
               f"sbatch --export=ALL,COPS_RUN={run},GRAIN_SET=mono {OUT}/model_ext/grain.sbatch",
               f"sbatch --dependency=afterok:<the grain jobs> --export=ALL,COPS_RUN={run} "
               f"{OUT}/model_ext/grain_finetune.sbatch"]
    if repro:
        md += ["# the reproductions plotted here: the source fit's recorded command, a new --out "
               "(<run>/grain/repro/<fit>_repro_<job>), then its fine tuning; the grain code is checked"] + \
              [f"sbatch --export=ALL,COPS_RUN={run},SRC_FIT={s} {OUT}/model_ext/grain_repro.sbatch" for s in repro]
    fixt = any("/grain/fixed_T/" in d for d in dirs)
    if fixt:
        md += ["# the fixed-temperature fit plotted here (T75 = Tmin = 133 K, Tmax free) and its fine tuning; the job "
               "refuses an existing output: move it aside as *.incomplete_<UTC> first",
               f"sbatch --export=ALL,COPS_RUN={run} {OUT}/model_ext/grain_fixT.sbatch"]
    md += ["```", "", "### A2. Per plotted case", "",
           "Each command runs inside a GPU job on one 8-GPU node (`salloc`/`srun`, or put it in a copy of the job "
           "script).", ""]
    for d in dirs:
        md += fit_entry(d, run)
    md += ["## B. From scratch", "",
           "A new restart directory; the existing one is never reused. `fit_base.in` must have the cluster "
           "`KernelDirectory` and `ObservationFile` lines (the only allowed local change).", "",
           "### B1. Runner campaign: base, record store, visibility, gravity, template bundle, fits, fine tuning", "",
           "The exact `cops_pipeline.py submit` that made this run (from `pipeline.json`), with a new "
           "`--restart-dir`. It writes the job scripts and never runs `sbatch`; `submit_campaign.sh` then saves "
           "the provenance and submits the prepare job and the C16 set (the later runner jobs submit themselves).",
           "", "```bash", f"R_NEW={run}_rerun_$(date -u +%Y%m%d)"]
    if inv:
        argv = list(inv["argv"])
        argv[argv.index("--restart-dir") + 1] = "$R_NEW"
        md += [wrap(q([PY, "-W", "ignore"] + argv)).replace("'$R_NEW'", "$R_NEW")]
    else:
        md += ["# pipeline.json has no submit invocation: see handoff_2026-10-26_cluster_setup_and_log.md §8.1"]
    md += [f"{OUT}/jobs/submit_campaign.sh $R_NEW", "```", ""]
    if grain:
        md += ["### B2. Grain sets and their fine tuning (after the prepare job: they need base, store, bundle and "
               "`c16/bundle`)", "", "```bash",
               f"a=$(sbatch --parsable --export=ALL,COPS_RUN=$R_NEW,GRAIN_SET=free {OUT}/model_ext/grain.sbatch)",
               f"b=$(sbatch --parsable --export=ALL,COPS_RUN=$R_NEW,GRAIN_SET=mono {OUT}/model_ext/grain.sbatch)",
               f"sbatch --dependency=afterok:$a:$b --export=ALL,COPS_RUN=$R_NEW {OUT}/model_ext/grain_finetune.sbatch"]
        if repro:
            md += ["# then the reproductions (after the grain jobs: they read the source fits' manifests)"] + \
                  [f"sbatch --dependency=afterok:$a:$b --export=ALL,COPS_RUN=$R_NEW,SRC_FIT={s} "
                   f"{OUT}/model_ext/grain_repro.sbatch" for s in repro]
        if fixt:
            md += ["# the fixed-temperature fit (T75 = Tmin = 133 K) and its fine tuning (after the prepare job)",
                   f"sbatch --export=ALL,COPS_RUN=$R_NEW {OUT}/model_ext/grain_fixT.sbatch"]
        md += ["```", ""]
    md += ["### B3. Figures", "", "```bash", a.figure_command.replace(str(run), "$R_NEW"), "```", "",
           "A rerun from scratch on other hardware, compiler or BLAS agrees at rounding level, not bitwise; "
           "2000 non-converged evaluations can amplify that (CLAUDE.md §3.12).", ""]
    (a.figdir / "HOW-TO-RERUN.md").write_text("\n".join(md))
    print(f"wrote {a.figdir / 'HOW-TO-RERUN.md'}")


if __name__ == "__main__":
    main()
