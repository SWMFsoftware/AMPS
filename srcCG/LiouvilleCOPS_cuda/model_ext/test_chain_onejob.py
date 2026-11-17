#!/usr/bin/env python3
# ==== version-controlled copy: header added 2026-11-13, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/model_ext/test_chain_onejob.py (on Voyager: /data001/highresolu_vt/cops_output/model_ext/test_chain_onejob.py)
# SHA256 of the run copy = SHA256 of this file without these header lines: 64154a0101597f1bb232095d3797b00a965af65c597ffcd520a4217f5410d8a9
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: test. CPU tests of chain_onejob.py (the planned stage commands equal, token by token, the commands.sh of the
#   split chains: positive jobs 123035-123052, the 3a chain 124008-124074, the nonneg job scripts; stage states; module
#   purge / restore; spawn; realfit; tee) and "compare [--generic] A B", the bitwise comparison of two chain roots
#   (split vs one job). User decision 2026-11-10 "yes, implement (b)".
# USED BY: login node / CPU jobs (realfit 124090, nonneg kind 124134); this version ran the nonneg comparisons 124160
#   and 124161 (ALL STAGES BITWISE IDENTICAL).
# 2026-11-13: test 1c softplus case with SLOPE_MULT_MIN=1 (token-equal to jobs/nonneg_fit.sbatch), test 3d (CPU job
#   124663); quick CPU set 87 checks pass.
# 2026-11-15: + 1d (stage-2 commands and the adoption of cycle 1 for sos / background / logspread) and the mode
#   nonnegstage2 (3d-style CPU job 131981); frozen chain_onejob_20261115T184608Z.
# 2026-11-16: + 1e (the fine-tuning delta regularization options, token-for-token; refusals).
# 2026-11-16: the configuration stage for --cycles 1 / 2 (1d, nonnegstage2 incl. a real --cycles 1 parent setup).
# 2026-11-17: + 1f (the free-space guard of the driver and of jobs/chain_onejob.sbatch).
# ==== end of header ====
"""Tests of model_ext/chain_onejob.py (CPU) and the split-vs-one-job comparison (after the GPU test jobs).

    .venv/bin/python -W ignore test_chain_onejob.py [--scratch DIR]          # CPU tests (login node: seconds)
    .venv/bin/python -W ignore test_chain_onejob.py compare <A root> <B root> [--json out.json] [--cycles N]

CPU tests
 1. commands: for the six stages of the reference chain runs/.../positive_r1e-5/lam1e-4 (jobs 123035 ... 123052) the
    driver's planned command equals the command recorded in that job's provenance commands.sh (token by token; the
    code path compared by file name, --out relative to the chain root). The predecessors of each stage are links into
    the reference chain (read only).
 1b. slope commands: the same for the --kind slope port (slope_mult_fit.sbatch CHAIN=1 / slope_finetune.sbatch) against
    the six stages of the completed reference 3a chain slope_bumps/chain_slopemul_20261110T201212Z (jobs 124008-124074).
 1c. nonneg commands: the planned fit / fine-tuning commands (softplus and squared; c1 fit, c1 fine tuning, c2 fit)
    == those of sandboxed copies of agent 1's jobs/nonneg_fit.sbatch / nonneg_finetune.sbatch (run line replaced by a
    print; fake SLURM_JOB_ID and nvidia-smi; snapshot model_ext_frozen/nonneg_20261111T044920Z).
 2. states: the ported positive_fit.sbatch / positive_finetune.sbatch logic on synthetic stage directories: new, complete,
    interrupted (next segment, evaluations left), a fit segment without restart (exit 3, or moved aside after an
    exit-75 event), START_RESTART mismatch (3), fine-tuning segment without state (moved aside), wall_time (resume),
    finished (final link + product links), product check failure (4), held identity of a completed cycle-2 fit.
 3. baseline: a stage module under a purge root patches kept modules (cops_cuda.figures.read_table,
    MultiGpuBackend.activity_linear, gpu_model._CLASS, reference._bound_root), the environment, sys.path, argv, numpy
    error state; Baseline.restore() undoes all of it and purges the stage module (in a subprocess).
 3b. real modules: the purge/restore on the real positive_fit / grain_fit / gpu_fit / fit_cops / positive_finetune
    modules (no GPU); spawn: chain_onejob.py as __mp_main__ of a spawned process (the GPU workers' import path).
 5. realfit (CPU Slurm job, 12M base; `test_chain_onejob.py realfit <old code dir> <new code dir> <scratch>`): the
    real cycle-1 fit stage through fit_cops.main() to GrainProblem's M.evaluate with a fake GPU backend: the old
    version refuses the second model (as job 123990), the new one accepts it; MF vs M differences; refusals.
 4. tee: output of this process and of a child process lands in the job log and, between switches, in the stage log
    (in a subprocess).
compare: bitwise comparison of two chain roots (split vs one-job, or two one-job runs): fit_result.txt numbers, every
array of fit_restart.npz (and meta optimizer/segments), activity_exp_facets.npz, COPS_model.dat, late_window_model.dat,
finetune_state.npz (arrays, meta costs), fine_tuning.json optimizer history and grain, finetune_result.txt numbers.
"""

import argparse
import json
import os
import re
import shlex
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
PY = "/data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda/.venv/bin/python"
ROOT_OUT = Path("/data001/highresolu_vt/cops_output")
RUN = ROOT_OUT / "runs" / "shap8_012m_lategrid_dz1.5"
REF = RUN / "positive_r1e-5" / "lam1e-4"
REF_SOURCE = RUN / "positive_r0.002" / "lam1e-4" / "c1" / "fit" / "final" / "fit_restart.npz"
FROZEN = ROOT_OUT / "model_ext_frozen" / "r1e-5k_20261109T023540Z"
POS_ENV = {"COPS_POS_ENABLE": "1", "COPS_POS_NIGHT_RATIO": "1e-5", "COPS_POS_ACTIVITY": "exp",
           "COPS_POS_FLOOR_RATIO": "1e-5", "COPS_POS_FLOOR_CO2": "1"}
FAILS = []


def check(cond, what):
    print("%s %s" % ("PASS" if cond else "FAIL", what), flush=True)
    if not cond:
        FAILS.append(what)
    return cond


def chain_args(root, extra=()):
    import chain_onejob as co
    return co.parse_args(["--chain-root", str(root), "--lambda", "1e-4", "--ratio", "1e-5", "--cycles", "3",
                          "--source", str(REF_SOURCE), "--nfev", "4000", "--check", "1", "--ft-maxnfev", "500",
                          "--ft-facets", "late_only", "--ft-window-min", "5", "--model-ext", str(FROZEN),
                          "--devices", "0,1,2,3,4,5,6,7"] + list(extra))


def reference_command(stage_dir):
    provs = sorted(p for p in stage_dir.glob("provenance-*") if p.is_dir())
    line = [l for l in (provs[0] / "commands.sh").read_text().splitlines() if l.strip()][-1]
    return shlex.split(line)


def normalize(tokens, root):
    out = []
    k = 0
    while k < len(tokens):
        t = tokens[k]
        if t.endswith(".py") and "/" in t:
            out.append("<code>/" + Path(t).name)
        elif t == "--out":
            out += [t, str(Path(tokens[k + 1]).relative_to(root))]
            k += 1
        else:
            out.append(t)
        k += 1
    return out


# ------------------------------------------------------------------------------------------------ 1 commands
def test_commands(scratch):
    import chain_onejob as co
    os.environ.update(POS_ENV)
    for name in ("c1fit", "c1ft", "c2fit", "c2ft", "c3fit", "c3ft"):
        fake = scratch / ("commands_" + name)
        k, kind = int(name[1]), name[2:]
        for j in range(1, k + (1 if kind == "ft" else 0)):                # completed predecessors: links into REF
            (fake / ("c%d" % j) / "fit").mkdir(parents=True, exist_ok=True)
            os.symlink(str(REF / ("c%d" % j) / "fit" / "seg01"), str(fake / ("c%d" % j) / "fit" / "final"))
        for j in range(1, k):
            (fake / ("c%d" % j) / "ft").mkdir(parents=True, exist_ok=True)
            os.symlink(str(REF / ("c%d" % j) / "ft" / "seg01" / "finetune_state.npz"),
                       str(fake / ("c%d" % j) / "ft" / "finetune_state.npz"))
        if kind == "ft":
            pass
        chain = co.Chain(chain_args(fake))
        stage = next(s for s in chain.stages if s.name == name)
        plan = stage.prepare()
        ok = check(plan.rc is None and not plan.done, "%s: planned (rc %s, notes %s)" % (name, plan.rc, plan.notes[-1:]))
        if not ok:
            continue
        mine = normalize(co.PY + plan.argv, fake)
        ref = normalize(reference_command(REF / ("c%d" % k) / ("fit" if kind == "fit" else "ft")), REF)
        same = mine == ref
        if not same:
            import difflib
            print("\n".join(difflib.unified_diff(ref, mine, "reference", "chain_onejob", lineterm="")))
        check(same, "%s: command == the reference job's commands.sh (%d tokens)" % (name, len(ref)))
        if kind == "fit":
            env_ref = (REF / ("c%d" % k) / "fit").glob("provenance-*/environment.txt")
            env_ref = dict(l.split("=", 1) for l in next(env_ref).read_text().splitlines() if "=" in l)
            for key in ("LAMBDA", "RATIO", "CYCLE", "NFEV", "CHECK", "CHAIN", "HELD", "GRAIN_FROM"):
                if key in env_ref or key in plan.env:
                    want = env_ref.get(key, "")
                    got = plan.env.get(key, "")
                    if key in ("HELD", "GRAIN_FROM"):
                        want = str(Path(want).relative_to(REF)) if want else ""
                        got = str(Path(got).relative_to(fake)) if got else ""
                    check(want == got, "%s: environment %s %r == reference %r" % (name, key, got, want))


# ------------------------------------------------------------------------------------------------ 1b slope commands
REF3A = RUN / "slope_bumps" / "chain_slopemul_20261110T201212Z"
D2_RESTART = RUN / "clip_refit" / "D2" / "final" / "fit_restart.npz"


def slope_chain_args(root, extra=()):
    import chain_onejob as co
    return co.parse_args(["--kind", "slope", "--slope-model", "bumps", "--chain-root", str(root), "--cycles", "3",
                          "--source", str(D2_RESTART), "--nfev", "4000", "--check", "1", "--ft-maxnfev", "500",
                          "--ft-facets", "late_only", "--ft-window-min", "5", "--l-threads", "32",
                          "--model-ext", str(ROOT_OUT / "model_ext_frozen" / "slopemul_20261110T221713Z"),
                          "--devices", "0,1,2,3,4,5,6,7"] + list(extra))


def reference_slope_command(stage_dir):
    provs = sorted(p for p in stage_dir.glob("provenance-*") if p.is_dir())
    line = [l for l in (provs[-1] / "commands.sh").read_text().splitlines() if l.strip()][-1]
    toks = shlex.split(line)
    return toks[toks.index(PY):]                    # drop the "bash -c 'echo 1000 > oom_score_adj; exec' _" prefix


def test_slope_commands(scratch):
    import chain_onejob as co
    if not REF3A.is_dir():
        check(False, "slope commands: reference chain %s missing" % REF3A)
        return
    for name in ("c1fit", "c1ft", "c2fit", "c2ft", "c3fit", "c3ft"):
        fake = scratch / ("slope_commands_" + name)
        k, kind = int(name[1]), name[2:]
        for j in range(1, k + (1 if kind == "ft" else 0)):
            (fake / ("c%d" % j) / "fit").mkdir(parents=True, exist_ok=True)
            os.symlink(str(REF3A / ("c%d" % j) / "fit" / "seg01"), str(fake / ("c%d" % j) / "fit" / "final"))
        for j in range(1, k):
            (fake / ("c%d" % j) / "ft").mkdir(parents=True, exist_ok=True)
            os.symlink(str(REF3A / ("c%d" % j) / "ft" / "seg01" / "finetune_state.npz"),
                       str(fake / ("c%d" % j) / "ft" / "finetune_state.npz"))
        chain = co.Chain(slope_chain_args(fake))
        stage = next(s for s in chain.stages if s.name == name)
        plan = stage.prepare()
        if not check(plan.rc is None and not plan.done, "slope %s: planned (rc %s %s)" % (name, plan.rc, plan.notes[-1:])):
            continue
        mine = normalize(co.PY + plan.argv, fake)
        ref = normalize(reference_slope_command(REF3A / ("c%d" % k) / kind), REF3A)
        if mine != ref:
            import difflib
            print("\n".join(difflib.unified_diff(ref, mine, "reference", "chain_onejob", lineterm="")))
        check(mine == ref, "slope %s: command == the reference 3a job's commands.sh (%d tokens)" % (name, len(ref)))


# ------------------------------------------------------------------------------------------------ 1c nonneg commands
NONNEG_SNAP = ROOT_OUT / "model_ext_frozen" / "nonneg_20261113T035423Z"     # (20261111T044920Z + --slope-mult-min)
NONNEG_TEST_START = RUN / "nonneg_tests" / "start" / "squared_start.npz"
RUN_LINES = {"nonneg_fit.sbatch": '''bash -c 'echo 1000 > /proc/self/oom_score_adj; exec "$@"' _ "${CMD[@]}" > "$P/fit.log" 2>&1; rc=$?''',
             "nonneg_finetune.sbatch": '''bash -c 'echo 1000 > /proc/self/oom_score_adj; exec "$@"' _ "${CMD[@]}" > "$P/finetune.log" 2>&1; rc=$?''',
             "chain_onejob.sbatch": '''bash -c 'echo 1000 > /proc/self/oom_score_adj; exec "$@"' _ "${CMD[@]}" > "$P/chain.log" 2>&1; rc=$?'''}


def emulate_sbatch(name, env_vars, scratch, tag):
    """Run a copy of jobs/<name> in a sandbox (fake SLURM_JOB_ID, fake nvidia-smi listing 8 GPUs) with its run line
    replaced by printing the command; returns the command tokens (from the python interpreter on)."""
    src = ROOT_OUT / "jobs" / name
    txt = src.read_text()
    if txt.count(RUN_LINES[name]) != 1:
        raise RuntimeError("the run line of %s changed" % name)
    out = scratch / ("emu_%s.cmd" % tag)
    txt = txt.replace(RUN_LINES[name], 'printf "%s\\n" "${CMD[@]}" > "$EMU_OUT"; kill ${MEMLOG:-} ${LOADLOG:-} 2>/dev/null; exit 0')
    s = scratch / ("emu_%s_%s" % (tag, name))
    s.write_text(txt)
    fakebin = scratch / "fakebin"
    fakebin.mkdir(exist_ok=True)
    (fakebin / "nvidia-smi").write_text('#!/bin/bash\nif [ "$1" = -L ]; then for i in 0 1 2 3 4 5 6 7; do echo "GPU $i: fake"; done; '
                                        'else echo "fake nvidia-smi $*"; fi\n')
    os.chmod(str(fakebin / "nvidia-smi"), 0o755)
    env = {k: v for k, v in os.environ.items() if not k.startswith(("COPS_POS_", "COPS_SLOPEMUL_", "COPS_NONNEG_",
                                                                     "SLURM_"))}
    env.update(env_vars, SLURM_JOB_ID="990000%03d" % (abs(hash(tag)) % 1000), EMU_OUT=str(out),
               PATH=str(fakebin) + ":" + os.environ.get("PATH", ""), CUDA_VISIBLE_DEVICES="0,1,2,3,4,5,6,7")
    r = subprocess.run(["bash", str(s)], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, universal_newlines=True,
                       env=env, timeout=600)
    if r.returncode != 0 or not out.is_file():
        raise RuntimeError("emulated %s (%s) exit %d: %s" % (name, tag, r.returncode, r.stdout[-1500:]))
    return [l for l in out.read_text().split("\n") if l != ""]


def _norm_roots(tokens, roots):
    out = []
    for t in tokens:
        if t.endswith(".py") and "/" in t:
            t = "<code>/" + Path(t).name
        for r, name in roots:
            t = t.replace(str(r), name)
        out.append(t)
    return out


def _fake_chain(root, cycles_done):
    """c1 fit (final -> seg01 with a restart file) and c1 ft (final + finetune_state.npz link): placeholders."""
    (root / "c1" / "fit" / "seg01").mkdir(parents=True)
    (root / "c1" / "fit" / "seg01" / "fit_restart.npz").write_text("x")
    for f in ("fit_manifest.json", "nonneg_facets.npz", "grain_term.json", "COPS_model.dat", "late_window_model.dat",
              "fit_result.txt"):
        (root / "c1" / "fit" / "seg01" / f).write_text("x")
    os.symlink("seg01", str(root / "c1" / "fit" / "final"))
    if cycles_done:
        (root / "c1" / "ft" / "seg01").mkdir(parents=True)
        (root / "c1" / "ft" / "seg01" / "finetune_state.npz").write_text("x")
        os.symlink("seg01", str(root / "c1" / "ft" / "final"))
        os.symlink("final/finetune_state.npz", str(root / "c1" / "ft" / "finetune_state.npz"))


def test_nonneg_commands(scratch):
    import chain_onejob as co
    if not NONNEG_SNAP.is_dir():
        check(False, "nonneg commands: snapshot %s missing" % NONNEG_SNAP)
        return
    d = scratch / "nonneg_commands"
    d.mkdir()
    common = {"NFEV": "5", "CHECK": "1", "CHAIN": "1", "L_THREADS": "32", "MODEL_EXT": str(NONNEG_SNAP)}
    smm_known = "SLOPE_MULT_MIN" in (ROOT_OUT / "jobs" / "nonneg_fit.sbatch").read_text()
    cases = [("softplus", None), ("squared", None)] + ([("softplus", "1")] if smm_known else [])
    if not smm_known:
        print("SKIP nonneg --slope-mult-min case: jobs/nonneg_fit.sbatch has no SLOPE_MULT_MIN yet (agent 1)")
    for model, smm in cases:
        extra = ["--start-file", str(NONNEG_TEST_START)] if model == "squared" else []
        sf = {"START_FILE": str(NONNEG_TEST_START)} if model == "squared" else {}
        # cycle 1 fit (empty roots), cycle 1 fine tuning and cycle 2 fit (c1 done), both sides the same layout
        smv = {"SLOPE_MULT_MIN": smm} if smm is not None else {}
        sma = ["--slope-mult-min", smm] if smm is not None else []
        tag0 = model + ("_smm%s" % smm if smm is not None else "")
        stages = (("c1fit", None), ("c1ft", 0), ("c2fit", 1)) + ((("c1fit_seg02", "resume"),) if smm else ())
        for stage, done in stages:
            split, mine = d / ("split_%s_%s" % (tag0, stage)), d / ("mine_%s_%s" % (tag0, stage))
            if done == "resume":                          # an interrupted segment 1 of the c1 fit (synthetic marker)
                for r in (split, mine):
                    seg = r / "c1" / "fit" / "seg01"
                    seg.mkdir(parents=True)
                    _npz(seg / "fit_restart.npz", {"status": "interrupted", "segments": [{"nfev": 2}],
                                                    "model_extension": {"fit_chain": {"start_restart": {"path": str(
                                                        (RUN / "clip_refit" / "D2" / "final" / "fit_restart.npz").resolve())}},
                                                        "start_choice": "d2"}})
            elif done is not None:
                _fake_chain(split, done)
                _fake_chain(mine, done)
            else:
                split.mkdir()
                mine.mkdir()
            if stage in ("c1fit", "c1fit_seg02"):
                ref = emulate_sbatch("nonneg_fit.sbatch", dict(common, MODEL=model, FITDIR=str(split / "c1" / "fit"),
                                                               **sf, **smv), d, "%s_%s" % (tag0, stage))
            elif stage == "c2fit":
                ref = emulate_sbatch("nonneg_fit.sbatch", dict(
                    common, MODEL=model, FITDIR=str(split / "c2" / "fit"), START="chain",
                    START_RESTART=str(split / "c1" / "fit" / "final" / "fit_restart.npz"),
                    HELD=str(split / "c1" / "ft" / "finetune_state.npz"),
                    GRAIN_FROM=str(split / "c1" / "ft" / "finetune_state.npz"), **smv), d, "%s_%s" % (tag0, stage))
            else:
                ref = emulate_sbatch("nonneg_finetune.sbatch", {
                    "FIT": str(split / "c1" / "fit" / "final"), "OUT": str(split / "c1" / "ft"), "FACETS": "late_only",
                    "WINDOW_MIN": "5", "MAXNFEV": "3", "L_THREADS": "32", "MODEL_EXT": str(NONNEG_SNAP)}, d,
                    "%s_%s" % (tag0, stage))
            args = co.parse_args(["--kind", "nonneg", "--slope-model", model, "--chain-root", str(mine), "--cycles", "2",
                                  "--nfev", "5", "--check", "1", "--ft-maxnfev", "3", "--ft-facets", "late_only",
                                  "--ft-window-min", "5", "--l-threads", "32", "--model-ext", str(NONNEG_SNAP),
                                  "--devices", "0,1,2,3,4,5,6,7", "--source",
                                  str(RUN / "clip_refit" / "D2" / "final" / "fit_restart.npz")] + extra + sma)
            chain = co.Chain(args)
            st = next(s for s in chain.stages if s.name == stage.split("_")[0])
            plan = st.prepare()
            if not check(plan.rc is None and not plan.done, "nonneg %s %s: planned (rc %s %s)" % (
                    tag0, stage, plan.rc, plan.notes[-1:])):
                continue
            a_ = _norm_roots(co.PY + plan.argv, [(mine, "<root>")])
            b_ = _norm_roots(ref, [(split, "<root>")])
            if a_ != b_:
                import difflib
                print("\n".join(difflib.unified_diff(b_, a_, "nonneg sbatch", "chain_onejob", lineterm="")))
            check(a_ == b_, "nonneg %s %s: command == the emulated %s (%d tokens)" % (
                tag0, stage, "nonneg_finetune.sbatch" if stage.endswith("ft") else "nonneg_fit.sbatch", len(b_)))


# ------------------------------------------------------------------------------------------------ 1d nonneg stage 2
# agent 1's options A, B, D (sos / background / logspread; snapshot nonneg_20261115T020314Z) and the lead's stage-2
# layout (2026-11-15): cycle 1's fit ADOPTED from the complete stage-1 continuation fit (--c1-from), the chain starts
# with the c1 fine tuning, c2.. fit --start chain; --slope-mult-min 1 everywhere.
NONNEG_SNAP2 = ROOT_OUT / "model_ext_frozen" / "nonneg_20261115T020314Z"
STAGE2 = {"background": ({"BG_MIN": "0.01"}, RUN / "nonneg" / "background" / "stage1c_bg1e-2_20261115T023450Z" / "fit"),
          "logspread": ({"SPREAD_F": "0.1", "SPREAD_CREF": "4.2688066575e-03", "SPREAD_U0": "-1.864126788612443"},
                        RUN / "nonneg" / "logspread" / "stage1c_f0.1_20261115T023450Z" / "fit"),
          "sos": ({"SOS_RANK": "3"}, RUN / "nonneg" / "sos" / "stage1c_r3_20261115T023450Z" / "fit")}


def _tok_args(tokens):
    out = []
    for k, v in tokens.items():
        out += ["--" + k.lower().replace("_", "-"), v]
    return out


def _stage2_args(co, model, root, tokens, c1_from=None, source=None, extra=()):
    src = source or (str((Path(c1_from) / "final" / "fit_restart.npz").resolve()) if c1_from
                     else str(RUN / "clip_refit" / "D2" / "final" / "fit_restart.npz"))
    return co.parse_args(["--kind", "nonneg", "--slope-model", model, "--chain-root", str(root), "--cycles", "2",
                          "--nfev", "5", "--check", "1", "--ft-maxnfev", "3", "--ft-facets", "late_only",
                          "--ft-window-min", "5", "--l-threads", "32", "--model-ext", str(NONNEG_SNAP2),
                          "--devices", "0,1,2,3,4,5,6,7", "--source", src, "--slope-mult-min", "1"]
                         + _tok_args(tokens) + (["--c1-from", str(c1_from)] if c1_from else []) + list(extra))


def _adopt_like_split(root, c1):
    """The split equivalent of the adoption: c1/fit/final -> the adopted fit's final segment (absolute)."""
    (root / "c1" / "fit").mkdir(parents=True)
    os.symlink(str((c1 / "final").resolve()), str(root / "c1" / "fit" / "final"))


def _fake_ft(root, k=1):
    (root / ("c%d" % k) / "ft" / "seg01").mkdir(parents=True)
    (root / ("c%d" % k) / "ft" / "seg01" / "finetune_state.npz").write_text("x")
    os.symlink("seg01", str(root / ("c%d" % k) / "ft" / "final"))
    os.symlink("final/finetune_state.npz", str(root / ("c%d" % k) / "ft" / "finetune_state.npz"))


def _compare_plan(co, plan, ref, mine, split, what, script):
    if not check(plan.rc is None and not plan.done and plan.argv, "%s: planned (rc %s %s)" % (what, plan.rc,
                                                                                         plan.notes[-1:])):
        return
    a_ = _norm_roots(co.PY + plan.argv, [(mine, "<root>")])
    b_ = _norm_roots(ref, [(split, "<root>")])
    if a_ != b_:
        import difflib
        print("\n".join(difflib.unified_diff(b_, a_, script, "chain_onejob", lineterm="")))
    check(a_ == b_, "%s: command == the emulated %s (%d tokens)" % (what, script, len(b_)))


def test_nonneg_stage2(scratch):
    import chain_onejob as co
    if not NONNEG_SNAP2.is_dir():
        check(False, "nonneg stage 2: snapshot %s missing" % NONNEG_SNAP2)
        return
    d = scratch / "nonneg_stage2"
    d.mkdir()
    common = {"NFEV": "5", "CHECK": "1", "CHAIN": "1", "L_THREADS": "32", "MODEL_EXT": str(NONNEG_SNAP2),
              "SLOPE_MULT_MIN": "1"}
    ftenv = {"FACETS": "late_only", "WINDOW_MIN": "5", "MAXNFEV": "3", "L_THREADS": "32", "MODEL_EXT": str(NONNEG_SNAP2)}
    for model, (tokens, c1) in STAGE2.items():
        if not (c1 / "final" / "fit_restart.npz").is_file():
            check(False, "nonneg stage 2 %s: the stage-1c fit %s is missing" % (model, c1))
            continue
        c1seg = (c1 / "final").resolve()
        # --- (i) a chain WITHOUT adoption: c1 fit (the model's default start), c1 ft, c2 fit (chain)
        for stage, done in (("c1fit", None), ("c1ft", 0), ("c2fit", 1)):
            split, mine = d / ("split_%s_%s" % (model, stage)), d / ("mine_%s_%s" % (model, stage))
            if done is None:
                split.mkdir(), mine.mkdir()
            else:
                _fake_chain(split, done), _fake_chain(mine, done)
            if stage == "c1fit":
                ref = emulate_sbatch("nonneg_fit.sbatch", dict(common, MODEL=model, FITDIR=str(split / "c1" / "fit"),
                                                               **tokens), d, "s2_%s_%s" % (model, stage))
                script = "nonneg_fit.sbatch"
            elif stage == "c2fit":
                ref = emulate_sbatch("nonneg_fit.sbatch", dict(
                    common, MODEL=model, FITDIR=str(split / "c2" / "fit"), START="chain",
                    START_RESTART=str(split / "c1" / "fit" / "final" / "fit_restart.npz"),
                    HELD=str(split / "c1" / "ft" / "finetune_state.npz"),
                    GRAIN_FROM=str(split / "c1" / "ft" / "finetune_state.npz"), **tokens), d, "s2_%s_%s" % (model, stage))
                script = "nonneg_fit.sbatch"
            else:
                ref = emulate_sbatch("nonneg_finetune.sbatch", dict(ftenv, FIT=str(split / "c1" / "fit" / "final"),
                                                                    OUT=str(split / "c1" / "ft")), d,
                                     "s2_%s_%s" % (model, stage))
                script = "nonneg_finetune.sbatch"
            chain = co.Chain(_stage2_args(co, model, mine, tokens))
            st = next(s for s in chain.stages if s.name == stage)
            _compare_plan(co, st.prepare(), ref, mine, split, "nonneg stage 2 %s %s (no adoption)" % (model, stage),
                          script)
        # --- (ii) the stage-2 chain: c1 adopted (--c1-from), c1 ft, c2 fit (chain from the adopted fit)
        mine, split = d / ("mine_%s_adopted" % model), d / ("split_%s_adopted" % model)
        mine.mkdir()
        _adopt_like_split(split, c1)
        chain = co.Chain(_stage2_args(co, model, mine, tokens, c1_from=c1))
        names = [(s.name, s.kind) for s in chain.stages]
        check(names[:2] == [("c1fit", "nonneg_fit_adopted"), ("c1ft", "nonneg_finetune")]
              and names[2] == ("c2fit", "nonneg_fit"), "nonneg stage 2 %s: layout %s" % (model, names))
        plan = chain.stages[0].prepare()
        link = mine / "c1" / "fit" / "final"
        txt = (mine / "c1" / "fit" / "ADOPTED_C1_FIT.txt").read_text() if \
            (mine / "c1" / "fit" / "ADOPTED_C1_FIT.txt").is_file() else ""
        check(plan.done and plan.adopted == str(c1seg) and link.is_symlink() and os.readlink(str(link)) == str(c1seg)
              and co.sha256_file(c1seg / "fit_restart.npz") in txt and "status complete" in txt,
              "nonneg stage 2 %s: c1 fit adopted: final -> %s (absolute), ADOPTED_C1_FIT.txt with the restart sha256"
              % (model, c1seg))
        again = co.Chain(_stage2_args(co, model, mine, tokens, c1_from=c1)).stages[0].prepare()
        check(again.done and again.adopted == str(c1seg), "nonneg stage 2 %s: adoption again (resubmission): "
                                                          "recognised, nothing changed" % model)
        ref = emulate_sbatch("nonneg_finetune.sbatch", dict(ftenv, FIT=str(split / "c1" / "fit" / "final"),
                                                            OUT=str(split / "c1" / "ft")), d, "s2_%s_adopted_c1ft" % model)
        _compare_plan(co, chain.stages[1].prepare(), ref, mine, split, "nonneg stage 2 %s adopted c1ft" % model,
                      "nonneg_finetune.sbatch")
        # c1 ft done (placeholders); c2 fit from the adopted fit
        for r in (mine, split):
            for e in sorted((r / "c1" / "ft").glob("*")) if (r / "c1" / "ft").is_dir() else []:
                os.rename(str(e), str(e) + ".planned_%s" % co.utc_stamp())
            (r / "c1" / "ft").mkdir(exist_ok=True)
            _fake_ft(r)
        ref = emulate_sbatch("nonneg_fit.sbatch", dict(
            common, MODEL=model, FITDIR=str(split / "c2" / "fit"), START="chain",
            START_RESTART=str(split / "c1" / "fit" / "final" / "fit_restart.npz"),
            HELD=str(split / "c1" / "ft" / "finetune_state.npz"),
            GRAIN_FROM=str(split / "c1" / "ft" / "finetune_state.npz"), **tokens), d, "s2_%s_adopted_c2fit" % model)
        chain = co.Chain(_stage2_args(co, model, mine, tokens, c1_from=c1))
        c2 = chain.stages[2]
        plan = c2.prepare()
        _compare_plan(co, plan, ref, mine, split, "nonneg stage 2 %s adopted c2fit" % model, "nonneg_fit.sbatch")
        if plan.argv:
            check(plan.argv[plan.argv.index("--start-restart") + 1] == str(c1seg / "fit_restart.npz"),
                  "nonneg stage 2 %s: c2 starts from the adopted restart" % model)
        # the parent configuration's argv when c1 is adopted: c2's model tokens, --start restart, no held term
        cargv = c2.config_argv()
        check("--held" not in cargv and cargv[cargv.index("--start") + 1] == "restart"
              and cargv[cargv.index("--start-restart") + 1] == str(c1seg / "fit_restart.npz")
              and all(cargv[cargv.index(o) + 1] == v for o, v in zip(_tok_args(tokens)[::2], _tok_args(tokens)[1::2]))
              and cargv[cargv.index("--slope-mult-min") + 1] == "1",
              "nonneg stage 2 %s: process configuration argv (c1 adopted) = c2's tokens with --start restart "
              "--start-restart <adopted>, no held term" % model)
        # the configuration stage: c2 here; with --cycles 1 (c1 adopted + c1 ft only) the adopted stage (v7.2; 135365)
        check(co.nonneg_config_stage(chain) is c2, "nonneg stage 2 %s: process configuration from c2 (cycles 2)" % model)
        r1c = d / ("cycles1_%s" % model)
        r1c.mkdir()
        a1 = _stage2_args(co, model, r1c, tokens, c1_from=c1, extra=["--cycles", "1"])
        ch1 = co.Chain(a1)
        st1 = co.nonneg_config_stage(ch1)
        ca1 = st1.config_argv()
        check([s_.name for s_ in ch1.stages] == ["c1fit", "c1ft"] and isinstance(st1, co.AdoptedNonnegFitStage)
              and "--held" not in ca1 and "--start-file" not in ca1 and ca1[ca1.index("--start") + 1] == "restart"
              and ca1[ca1.index("--start-restart") + 1] == str(c1seg / "fit_restart.npz")
              and all(ca1[ca1.index(o) + 1] == v for o, v in zip(_tok_args(tokens)[::2], _tok_args(tokens)[1::2])),
              "nonneg stage 2 %s --cycles 1: process configuration from the adopted stage (--start restart "
              "--start-restart <adopted>, tokens, no held / start file)" % model)
        # resume of c2 (segment 2): a synthetic marker with this chain's tokens is continued; another value refused
        fl = {k.lower(): float(v) for k, v in tokens.items()}
        ext = {"model": {"model": model, "sos_rank": int(fl["sos_rank"]) if "sos_rank" in fl else None},
               "slope_mult_min": 1.0, "bg_min": fl.get("bg_min"), "spread_f": fl.get("spread_f"),
               "spread_cref": fl.get("spread_cref"), "spread_u0": fl.get("spread_u0"), "start_choice": "chain",
               "fit_chain": {"start_restart": {"path": str(c1seg / "fit_restart.npz")}}}
        for r in (mine, split):
            _npz(r / "c2" / "fit" / "seg01" / "fit_restart.npz", {"status": "interrupted", "segments": [{"nfev": 2}],
                                                                  "model_extension": ext})
        ref = emulate_sbatch("nonneg_fit.sbatch", dict(
            common, MODEL=model, FITDIR=str(split / "c2" / "fit"), START="chain",
            START_RESTART=str(split / "c1" / "fit" / "final" / "fit_restart.npz"),
            HELD=str(split / "c1" / "ft" / "finetune_state.npz"),
            GRAIN_FROM=str(split / "c1" / "ft" / "finetune_state.npz"), **tokens), d, "s2_%s_c2fit_seg02" % model)
        _compare_plan(co, co.Chain(_stage2_args(co, model, mine, tokens, c1_from=c1)).stages[2].prepare(), ref, mine,
                      split, "nonneg stage 2 %s adopted c2fit resume (segment 2)" % model, "nonneg_fit.sbatch")
        bad = dict(ext, slope_mult_min=None) if model != "sos" else dict(ext, model={"model": "sos", "sos_rank": 2})
        _npz(mine / "c2" / "fit" / "seg01" / "fit_restart.npz", {"status": "interrupted", "segments": [{"nfev": 2}],
                                                                 "model_extension": bad})
        p2 = co.Chain(_stage2_args(co, model, mine, tokens, c1_from=c1)).stages[2].prepare()
        check(p2.rc == 3 and "not this chain's model" in " ".join(p2.notes),
              "nonneg stage 2 %s: a c2 segment with other model tokens: refused (exit 3): %s" % (model, p2.notes[-1:]))
        # adoption refusals
        other = {"background": "logspread", "logspread": "sos", "sos": "background"}[model]
        r1 = d / ("refuse_%s_model" % model)
        r1.mkdir()
        try:
            pr = co.Chain(_stage2_args(co, model, r1, tokens, c1_from=STAGE2[other][1])).stages[0].prepare()
            check(pr.rc == 3 and not (r1 / "c1" / "fit" / "final").exists(),
                  "nonneg stage 2 %s: --c1-from the %s fit refused (exit 3), nothing linked: %s" % (model, other,
                                                                                                  pr.notes[-1:]))
        except SystemExit as exc:
            check(False, "nonneg stage 2 %s: --c1-from another model raised %s" % (model, exc))
        wrong = {"background": {"BG_MIN": "0.02"}, "logspread": dict(tokens, SPREAD_F="0.2"),
                 "sos": {"SOS_RANK": "2"}}[model]
        r2 = d / ("refuse_%s_token" % model)
        r2.mkdir()
        pr = co.Chain(_stage2_args(co, model, r2, wrong, c1_from=c1)).stages[0].prepare()
        check(pr.rc == 3 and not (r2 / "c1" / "fit" / "final").exists(),
              "nonneg stage 2 %s: --c1-from with another token value refused (exit 3): %s" % (model, pr.notes[-1:]))
        r3 = d / ("refuse_%s_segs" % model)
        (r3 / "c1" / "fit" / "seg01").mkdir(parents=True)
        pr = co.Chain(_stage2_args(co, model, r3, tokens, c1_from=c1)).stages[0].prepare()
        check(pr.rc == 3 and not (r3 / "c1" / "fit" / "final").exists(),
              "nonneg stage 2 %s: c1/fit with segments: not adopted (exit 3)" % model)
        r4 = d / ("refuse_%s_link" % model)
        (r4 / "c1" / "fit").mkdir(parents=True)
        os.symlink(str(STAGE2[other][1] / "final"), str(r4 / "c1" / "fit" / "final"))
        pr = co.Chain(_stage2_args(co, model, r4, tokens, c1_from=c1)).stages[0].prepare()
        check(pr.rc == 3, "nonneg stage 2 %s: c1/fit/final pointing elsewhere: refused (exit 3)" % model)
    # option checks
    for what, mdl, tok, extra, kw in (
            ("--c1-from with softplus", "softplus", {}, [], {"c1_from": STAGE2["background"][1]}),
            ("--bg-min with sos", "sos", {"SOS_RANK": "3", "BG_MIN": "0.01"}, [], {}),
            ("logspread without --spread-cref", "logspread", {"SPREAD_F": "0.1"}, [], {}),
            ("sos without --sos-rank", "sos", {}, [], {}),
            ("--slope-mult-min 0.5", "background", {"BG_MIN": "0.01"}, ["--slope-mult-min", "0.5"], {}),
            ("--c1-from with another --source", "background", {"BG_MIN": "0.01"}, [],
             {"c1_from": STAGE2["background"][1], "source": str(RUN / "clip_refit" / "D2" / "final" / "fit_restart.npz")}),
            ("--c1-from with --slope-start", "background", {"BG_MIN": "0.01"}, ["--slope-start", "d2"],
             {"c1_from": STAGE2["background"][1]})):
        try:
            import contextlib, io
            with contextlib.redirect_stderr(io.StringIO()):
                _stage2_args(co, mdl, d / "x", tok, extra=extra, **kw)
            check(False, "nonneg stage 2 options: %s accepted" % what)
        except SystemExit as exc:
            check(exc.code == 2, "nonneg stage 2 options: %s refused (usage error)" % what)
    # jobs/chain_onejob.sbatch -> the driver's options (emulated sbatch, then parse_args of its command)
    for model, (tokens, c1) in STAGE2.items():
        root = d / ("sbatch_%s" % model)
        env = dict(KIND="nonneg", SLOPE_MODEL=model, CYCLES="5", NFEV="4000", CHECK="1", FT_MAXNFEV="500",
                   FT_FACETS="late_only", FT_WINDOW_MIN="5", L_THREADS="32", SLOPE_MULT_MIN="1", C1_FROM=str(c1),
                   CHAIN_ROOT=str(root), MODEL_EXT=str(NONNEG_SNAP2), CHAIN_CODE=str(HERE), **tokens)
        try:
            cmd = emulate_sbatch("chain_onejob.sbatch", env, d, "s2_sbatch_%s" % model)
        except Exception as exc:
            check(False, "nonneg stage 2 %s: emulated chain_onejob.sbatch: %s" % (model, exc))
            continue
        i = next(k for k, t in enumerate(cmd) if t.endswith("chain_onejob.py"))
        a = co.parse_args(cmd[i + 1:])
        want = {k.lower(): v for k, v in tokens.items()}
        got = {k: getattr(a, k) for k in want}
        check(got == want and a.c1_from == str(c1.resolve()) and a.slope_mult_min == "1" and a.cycles == 5
              and a.nfev == 4000 and a.ft_maxnfev == 500 and a.l_threads == 32
              and readlink(a.source) == str((c1 / "final" / "fit_restart.npz").resolve())
              and all(k not in cmd for k in ("--slope-start", "--start-file")),
              "nonneg stage 2 %s: chain_onejob.sbatch passes the tokens verbatim (%s), --c1-from, --source = the "
              "adopted restart, --slope-mult-min 1, 5 cycles, NFEV 4000, MAXNFEV 500, 32 threads" % (model, got))


def readlink(p):
    return os.path.realpath(str(p))


# ------------------------------------------------------------------------------------------------ 1e delta regularization
# agent 1 (2026-11-16, snapshot nonneg_20261116T031115Z): jobs/nonneg_finetune.sbatch DELTA_F / DELTA_CREF / DELTA_CAP ->
# --delta-penalty-f / --delta-cref / --delta-cap after --l-threads (TOPT), before RES and --max-nfev; the lead's settings
# for the stage-2 reruns: 10, 4.2688066575e-03, 4.605170185988092.
NONNEG_SNAP3 = ROOT_OUT / "model_ext_frozen" / "nonneg_20261116T031115Z"
DELTA = {"DELTA_F": "10", "DELTA_CREF": "4.2688066575e-03", "DELTA_CAP": "4.605170185988092"}
DELTA_ARGS = ["--ft-delta-penalty-f", "10", "--ft-delta-cref", "4.2688066575e-03", "--ft-delta-cap", "4.605170185988092"]


def test_nonneg_delta(scratch):
    import chain_onejob as co
    if not NONNEG_SNAP3.is_dir():
        check(False, "nonneg delta: snapshot %s missing" % NONNEG_SNAP3)
        return
    if "DELTA_F" not in (ROOT_OUT / "jobs" / "nonneg_finetune.sbatch").read_text():
        check(False, "nonneg delta: jobs/nonneg_finetune.sbatch has no DELTA_F")
        return
    d = scratch / "nonneg_delta"
    d.mkdir()
    ftenv = {"FACETS": "late_only", "WINDOW_MIN": "5", "MAXNFEV": "3", "L_THREADS": "32", "MODEL_EXT": str(NONNEG_SNAP3)}
    fitenv = {"NFEV": "5", "CHECK": "1", "CHAIN": "1", "L_THREADS": "32", "MODEL_EXT": str(NONNEG_SNAP3),
              "SLOPE_MULT_MIN": "1"}

    def args_for(model, root, tokens=None, c1_from=None, delta=True, extra=()):
        src = str((Path(c1_from) / "final" / "fit_restart.npz").resolve()) if c1_from else \
            str(RUN / "clip_refit" / "D2" / "final" / "fit_restart.npz")
        return co.parse_args(["--kind", "nonneg", "--slope-model", model, "--chain-root", str(root), "--cycles", "2",
                              "--nfev", "5", "--check", "1", "--ft-maxnfev", "3", "--ft-facets", "late_only",
                              "--ft-window-min", "5", "--l-threads", "32", "--model-ext", str(NONNEG_SNAP3),
                              "--devices", "0,1,2,3,4,5,6,7", "--source", src, "--slope-mult-min", "1"]
                             + _tok_args(tokens or {}) + (["--c1-from", str(c1_from)] if c1_from else [])
                             + (DELTA_ARGS if delta else []) + list(extra))
    # softplus (the reference chain of the reruns): c1 ft, its resumed segment 2, c2 fit (no delta in a fit)
    split, mine = d / "split_softplus", d / "mine_softplus"
    _fake_chain(split, 0), _fake_chain(mine, 0)
    ref = emulate_sbatch("nonneg_finetune.sbatch", dict(ftenv, FIT=str(split / "c1" / "fit" / "final"),
                                                        OUT=str(split / "c1" / "ft"), **DELTA), d, "dl_sp_c1ft")
    chain = co.Chain(args_for("softplus", mine))
    plan = chain.stages[1].prepare()
    _compare_plan(co, plan, ref, mine, split, "nonneg delta softplus c1ft", "nonneg_finetune.sbatch")
    if plan.argv:
        k = plan.argv.index("--delta-penalty-f")
        check(plan.argv[k - 2:k + 6] == ["--l-threads", "32", "--delta-penalty-f", "10", "--delta-cref",
                                         "4.2688066575e-03", "--delta-cap", "4.605170185988092"]
              and plan.env.get("DELTA_CAP") == "4.605170185988092",
              "nonneg delta softplus c1ft: the three options right after --l-threads, verbatim; DELTA_* in the stage env")
    fsha = co.sha256_file(mine / "c1" / "fit" / "seg01" / "fit_restart.npz")
    dmeta = {"f": 10.0, "C_ref": 4.2688066575e-03, "cap": 4.605170185988092}
    for r in (split, mine):
        _npz(r / "c1" / "ft" / "seg01" / "finetune_state.npz", {"stop": "wall_time", "nfev_total": 2, "complete": False,
                                                                "fit": {"fit_restart_sha256": fsha},
                                                                "delta_regularization": dmeta})
    ref = emulate_sbatch("nonneg_finetune.sbatch", dict(ftenv, FIT=str(split / "c1" / "fit" / "final"),
                                                        OUT=str(split / "c1" / "ft"), **DELTA), d, "dl_sp_c1ft_seg02")
    plan = co.Chain(args_for("softplus", mine)).stages[1].prepare()
    _compare_plan(co, plan, ref, mine, split, "nonneg delta softplus c1ft resume (segment 2)", "nonneg_finetune.sbatch")
    p2 = co.Chain(args_for("softplus", mine, delta=False)).stages[1].prepare()
    check(p2.rc == 3 and "delta regularization" in " ".join(p2.notes),
          "nonneg delta softplus: a regularized segment continued WITHOUT the options: refused (exit 3)")
    # completed unregularized fine tuning + chain with the options: refused (other setting)
    r2 = d / "mine_softplus_done_plain"
    _fake_chain(r2, 0)
    _npz(r2 / "c1" / "ft" / "seg01" / "finetune_state.npz", {"stop": "max_nfev", "nfev_total": 3, "complete": True,
                                                             "fit": {"fit_restart_sha256": fsha}})
    p4 = co.Chain(args_for("softplus", r2)).stages[1].prepare()
    check(p4.rc == 3, "nonneg delta softplus: a completed UNregularized fine tuning in a regularized chain: refused "
                      "(exit 3)")
    p5 = co.Chain(args_for("softplus", r2, delta=False)).stages[1].prepare()
    check(p5.rc != 3 or "delta" not in " ".join(p5.notes),
          "nonneg delta softplus: the same state without the options: no delta refusal (old behaviour)")
    split2, mine2 = d / "split_softplus_c2", d / "mine_softplus_c2"
    _fake_chain(split2, 1), _fake_chain(mine2, 1)
    ref = emulate_sbatch("nonneg_fit.sbatch", dict(
        fitenv, MODEL="softplus", FITDIR=str(split2 / "c2" / "fit"), START="chain",
        START_RESTART=str(split2 / "c1" / "fit" / "final" / "fit_restart.npz"),
        HELD=str(split2 / "c1" / "ft" / "finetune_state.npz"),
        GRAIN_FROM=str(split2 / "c1" / "ft" / "finetune_state.npz")), d, "dl_sp_c2fit")
    _compare_plan(co, co.Chain(args_for("softplus", mine2)).stages[2].prepare(), ref, mine2, split2,
                  "nonneg delta softplus c2fit (fits take no delta option)", "nonneg_fit.sbatch")
    # the adopted stage-2 chains: c1 ft with the options
    for model, (tokens, c1) in STAGE2.items():
        mine, split = d / ("mine_%s" % model), d / ("split_%s" % model)
        mine.mkdir()
        _adopt_like_split(split, c1)
        chain = co.Chain(args_for(model, mine, tokens, c1_from=c1))
        ad = chain.stages[0].prepare()
        ref = emulate_sbatch("nonneg_finetune.sbatch", dict(ftenv, FIT=str(split / "c1" / "fit" / "final"),
                                                            OUT=str(split / "c1" / "ft"), **DELTA), d, "dl_%s_c1ft" % model)
        check(bool(ad.done), "nonneg delta %s: c1 adopted" % model)
        _compare_plan(co, chain.stages[1].prepare(), ref, mine, split, "nonneg delta %s adopted c1ft" % model,
                      "nonneg_finetune.sbatch")
    # options
    for what, argv in (("--ft-delta-penalty-f without --ft-delta-cref", ["--ft-delta-penalty-f", "10"]),
                       ("--ft-delta-cref without --ft-delta-penalty-f", ["--ft-delta-cref", "1e-3"]),
                       ("--ft-delta-cap not a number", ["--ft-delta-cap", "big"])):
        try:
            import contextlib, io
            with contextlib.redirect_stderr(io.StringIO()):
                args_for("softplus", d / "x", delta=False, extra=argv)
            check(False, "nonneg delta options: %s accepted" % what)
        except SystemExit as exc:
            check(exc.code == 2, "nonneg delta options: %s refused (usage error)" % what)
    try:
        import contextlib, io
        with contextlib.redirect_stderr(io.StringIO()):
            co.parse_args(["--chain-root", str(d / "pos"), "--lambda", "1e-4", "--ratio", "1e-5", "--source", "x",
                           "--ft-maxnfev", "3", "--ft-facets", "late_only", "--ft-window-min", "5", "--model-ext", "x",
                           "--ft-delta-penalty-f", "10", "--ft-delta-cref", "1e-3"])
        check(False, "nonneg delta options: accepted with kind positive")
    except SystemExit as exc:
        check(exc.code == 2, "nonneg delta options: refused with kind positive (usage error)")
    # jobs/chain_onejob.sbatch -> driver options
    root = d / "sbatch_softplus"
    env = dict(KIND="nonneg", SLOPE_MODEL="softplus", CYCLES="5", NFEV="4000", CHECK="1", FT_MAXNFEV="500",
               FT_FACETS="late_only", FT_WINDOW_MIN="5", L_THREADS="32", SLOPE_MULT_MIN="1", CHAIN_ROOT=str(root),
               MODEL_EXT=str(NONNEG_SNAP3), CHAIN_CODE=str(HERE), FT_DELTA_F="10", FT_DELTA_CREF="4.2688066575e-03",
               FT_DELTA_CAP="4.605170185988092")
    cmd = emulate_sbatch("chain_onejob.sbatch", env, d, "dl_sbatch_softplus")
    i = next(k for k, t in enumerate(cmd) if t.endswith("chain_onejob.py"))
    a = co.parse_args(cmd[i + 1:])
    check(a.ft_delta_penalty_f == "10" and a.ft_delta_cref == "4.2688066575e-03" and a.ft_delta_cap == "4.605170185988092"
          and a.slope_model == "softplus" and a.c1_from is None and a.slope_mult_min == "1"
          and readlink(a.source) == readlink(RUN / "clip_refit" / "D2" / "final" / "fit_restart.npz"),
          "nonneg delta: chain_onejob.sbatch passes FT_DELTA_F / FT_DELTA_CREF / FT_DELTA_CAP verbatim (softplus, "
          "source D2 final)")


# ------------------------------------------------------------------------------------------------ 1f free-space guard
# v7.3 (lead 2026-11-17, after the disk-full event of 2026-11-16 21:10 CST): --min-free-gb / MIN_FREE_GB.
DISK_CHILD = r'''
import sys
sys.path.insert(0, sys.argv[1])
import chain_onejob as co
args = co.parse_args(sys.argv[3:])
sys.exit(co.Chain(args).run())
'''


def test_disk_guard(scratch):
    import chain_onejob as co
    d = scratch / "disk_guard"
    d.mkdir()
    root = d / "chain"
    base = ["--kind", "nonneg", "--slope-model", "softplus", "--chain-root", str(root), "--cycles", "2", "--nfev", "5",
            "--check", "1", "--ft-maxnfev", "3", "--ft-facets", "late_only", "--ft-window-min", "5", "--l-threads", "32",
            "--model-ext", str(NONNEG_SNAP3), "--devices", "0", "--source",
            str(RUN / "clip_refit" / "D2" / "final" / "fit_restart.npz"), "--slope-mult-min", "1"]
    a = co.parse_args(base)
    check(a.min_free_gb is None and co.Chain(a).disk_short("x") is None, "disk guard: off by default")
    check(co.Chain(co.parse_args(base + ["--min-free-gb", "0.001"])).disk_short("x") is None,
          "disk guard: enough free space -> no stop")
    msg = co.Chain(co.parse_args(base + ["--min-free-gb", "1e9"])).disk_short("stage c1fit")
    check(bool(msg) and "exit 5" in msg, "disk guard: below the limit -> message %r" % (msg or "")[:120])
    r = subprocess.run([PY, "-W", "ignore", "-c", DISK_CHILD, str(HERE), "-"] + base + ["--min-free-gb", "1e9"],
                       stdout=subprocess.PIPE, stderr=subprocess.STDOUT, universal_newlines=True, timeout=600,
                       cwd="/data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda")
    made = sorted(str(x.relative_to(root)) for x in root.rglob("*")) if root.exists() else []
    check(r.returncode == co.EXIT_DISK and not [m for m in made if m.startswith("c")]
          and "below --min-free-gb" in r.stdout and "process setup (" not in r.stdout,
          "disk guard: run() stops with exit 5 before the process setup, no stage directory written (%s; %s)"
          % (r.returncode, made))
    # a completed stage (final link) is not refused; the first unfinished one is
    root2 = d / "chain2"
    _fake_chain(root2, 1)
    ch = co.Chain(co.parse_args([x if x != str(root) else str(root2) for x in base] + ["--min-free-gb", "1e9"]))
    hints = [st.state_complete_hint() for st in ch.stages]
    check(hints == [True, True, False, False], "disk guard: completed stages recognised by their final link %s" % hints)
    # the job script: below MIN_FREE_GB nothing is written (exit 5); otherwise --min-free-gb (default 20) reaches the driver
    env = dict(KIND="nonneg", SLOPE_MODEL="softplus", CYCLES="5", NFEV="4000", CHECK="1", FT_MAXNFEV="500",
               FT_FACETS="late_only", FT_WINDOW_MIN="5", L_THREADS="32", SLOPE_MULT_MIN="1",
               MODEL_EXT=str(NONNEG_SNAP3), CHAIN_CODE=str(HERE))
    r3 = d / "sbatch_low"
    try:
        emulate_sbatch("chain_onejob.sbatch", dict(env, CHAIN_ROOT=str(r3), MIN_FREE_GB="1000000"), d, "disk_low")
        check(False, "disk guard: chain_onejob.sbatch ran below MIN_FREE_GB")
    except RuntimeError as exc:
        check("exit 5" in str(exc) and "below MIN_FREE_GB" in str(exc) and not r3.exists(),
              "disk guard: chain_onejob.sbatch below MIN_FREE_GB: exit 5, chain root not created")
    cmd = emulate_sbatch("chain_onejob.sbatch", dict(env, CHAIN_ROOT=str(d / "sbatch_default")), d, "disk_default")
    i = next(k for k, t in enumerate(cmd) if t.endswith("chain_onejob.py"))
    a4 = co.parse_args(cmd[i + 1:])
    cmd0 = emulate_sbatch("chain_onejob.sbatch", dict(env, CHAIN_ROOT=str(d / "sbatch_off"), MIN_FREE_GB="0"), d, "disk_off")
    check(a4.min_free_gb == 20.0 and "--min-free-gb" not in cmd0,
          "disk guard: chain_onejob.sbatch passes --min-free-gb 20 by default, nothing with MIN_FREE_GB=0")


# ------------------------------------------------------------------------------------------------ 2 states
def _npz(path, meta, **arrays):
    path.parent.mkdir(parents=True, exist_ok=True)
    np.savez(str(path), meta_json=np.array(json.dumps(meta)), **arrays)


def _fake_run(d):
    base = d / "run" / "base"
    (base / "kernels").mkdir(parents=True)
    for f in ("COPS_model.dat", "elements.bin"):
        (base / f).write_text("x")
    (d / "run" / "late_grid_points.txt").write_text("x")
    (d / "bundle.npz").write_text("x")
    return d / "run"


def test_states(scratch):
    import chain_onejob as co
    os.environ.update(POS_ENV)
    d = scratch / "states"
    run = _fake_run(d)
    src = d / "src" / "fit_restart.npz"
    _npz(src, {"status": "complete", "settings": {"lambda": 1e-4}})
    extra = ["--run", str(run), "--point-table", str(run / "late_grid_points.txt"), "--bundle", str(d / "bundle.npz"),
             "--source", str(src), "--nfev", "100"]
    root = d / "chain"
    chain = co.Chain(chain_args(root, extra))
    fit1, ft1, fit2 = chain.stages[0], chain.stages[1], chain.stages[2]
    p = fit1.prepare()
    check(p.segment == 1 and "--max-nfev" in p.argv and p.argv[p.argv.index("--max-nfev") + 1] == "100"
          and "--check-jacobian" in p.argv and p.argv[p.argv.index("--start-restart") + 1] == str(src),
          "fit new: segment 1, 100 evaluations, --check-jacobian, --start-restart")
    chain_rec = {"start_restart": {"path": str(src)}}
    _npz(root / "c1/fit/seg01/fit_restart.npz", {"status": "interrupted", "segments": [{"nfev": 30}],
                                                 "model_extension": {"fit_chain": chain_rec}})
    p = fit1.prepare()
    check(p.segment == 2 and p.argv[p.argv.index("--max-nfev") + 1] == "70" and "--check-jacobian" not in p.argv
          and p.argv[p.argv.index("--restart-from") + 1] == str(root / "c1/fit/seg01/fit_restart.npz"),
          "fit interrupted: segment 2 with 70 evaluations, --restart-from seg01, no --check-jacobian")
    (root / "c1/fit/seg02").mkdir()
    p = fit1.prepare()
    check(p.rc == 3, "fit segment without restart and no exit-75 event: exit 3 (CHAIN=1)")
    chain.events.append(event="end", stage="c1fit", segment_dir=str(root / "c1/fit/seg02"), rc=75, restart_written=False)
    p = fit1.prepare()
    moved = list((root / "c1/fit").glob("seg02.incomplete_*"))
    check(p.segment == 2 and len(moved) == 1, "fit segment without restart after exit 75: moved aside (%s), rerun"
          % [m.name for m in moved])
    _npz(root / "c1/fit/seg02/fit_restart.npz", {"status": "complete", "segments": [{"nfev": 30}, {"nfev": 70}],
                                                 "model_extension": {"fit_chain": chain_rec}})
    p = fit1.prepare()
    check(p.done and os.readlink(str(root / "c1/fit/final")) == "seg02", "fit complete: done, final -> seg02")
    other = co.Chain(chain_args(root, extra[:-4] + ["--source", str(d / "other.npz"), "--nfev", "100"]))
    _npz(d / "other.npz", {"status": "complete", "settings": {"lambda": 1e-4}})
    check(other.stages[0].prepare().rc == 3, "fit with another START_RESTART: exit 3")
    # --- fine tuning of c1
    fitdir = root / "c1/fit/seg02"
    (fitdir / "fit_manifest.json").write_text(json.dumps({"model_extension": {"positive": {"environment": POS_ENV}}}))
    (fitdir / "activity_exp_facets.npz").write_text("x")
    sha = co.sha256_file(fitdir / "fit_restart.npz")
    p = ft1.prepare()
    check(p.segment == 1 and "--validate" in p.argv and p.argv[p.argv.index("--fit") + 1] == str(fitdir)
          and p.argv[p.argv.index("--floor-ratio") + 1] == "1e-5", "ft new: segment 1, --validate, --fit = the final seg")
    os.environ["COPS_POS_NIGHT_RATIO"] = "0.002"
    check(ft1.prepare().rc == 2, "ft with another COPS_POS_* than the fit: refused (2)")
    os.environ.update(POS_ENV)
    (root / "c1/ft/seg01").mkdir(parents=True)
    p = ft1.prepare()
    check(p.segment == 1 and len(list((root / "c1/ft").glob("seg01.incomplete_*"))) == 1,
          "ft segment without state: moved aside, segment 1 again")
    nf = 5
    h0 = np.zeros(nf)
    dlt = np.array([0.0, 0.1, 0.0, -0.2, 0.0])
    _npz(root / "c1/ft/seg01/finetune_state.npz", {"stop": "wall_time", "nfev_total": 120, "fit_restart_sha256": sha,
                                                   "complete": False},
         held_next=h0 + dlt, delta=dlt, held_prev=h0, grain_A_uPa=np.float64(7.0), grain_B_uPa_per_100s=np.float64(0.0))
    p = ft1.prepare()
    check(p.segment == 2 and p.argv[p.argv.index("--max-nfev") + 1] == "380" and "--validate" not in p.argv
          and p.argv[p.argv.index("--resume-from") + 1] == str(root / "c1/ft/seg01/finetune_state.npz")
          and p.argv[p.argv.index("--nfev-before") + 1] == "120", "ft wall_time: segment 2, 380 left, --resume-from")
    _npz(root / "c1/ft/seg02/finetune_state.npz", {"stop": "max_nfev", "nfev_total": 500, "fit_restart_sha256": sha,
                                                   "complete": True},
         held_next=h0 + 2 * dlt, delta=dlt, held_prev=h0, grain_A_uPa=np.float64(7.0),
         grain_B_uPa_per_100s=np.float64(0.0))
    check(ft1.prepare().rc == 4, "ft product check: held_next != held_prev + delta -> 4")
    _npz(root / "c1/ft/seg02/finetune_state.npz", {"stop": "max_nfev", "nfev_total": 500, "fit_restart_sha256": sha,
                                                   "complete": True},
         held_next=h0 + dlt, delta=dlt, held_prev=h0, grain_A_uPa=np.float64(7.0), grain_B_uPa_per_100s=np.float64(0.0))
    (root / "c1/ft/seg02/finetune_result.txt").write_text("x")
    p = ft1.prepare()
    check(p.done and os.readlink(str(root / "c1/ft/final")) == "seg02"
          and os.readlink(str(root / "c1/ft/finetune_state.npz")) == "final/finetune_state.npz"
          and os.readlink(str(root / "c1/ft/finetune_result.txt")) == "final/finetune_result.txt"
          and not (root / "c1/ft/COPS_model.dat").exists(), "ft finished: final -> seg02, product links")
    # --- cycle 2 fit: start from c1's fit, held = c1's fine tuning
    p = fit2.prepare()
    held = str((root / "c1/ft/seg02/finetune_state.npz").resolve())
    check(p.argv[p.argv.index("--start-restart") + 1] == str(fitdir / "fit_restart.npz")
          and p.argv[p.argv.index("--held-exponent") + 1] == held
          and p.argv[p.argv.index("--grain-start-from") + 1] == held, "fit c2: --start-restart c1 final, held/grain c1 ft")
    _npz(root / "c2/fit/seg01/fit_restart.npz", {"status": "complete", "segments": [{"nfev": 100}],
                                                 "model_extension": {"fit_chain": {"start_restart": {
                                                     "path": str(fitdir / "fit_restart.npz")}},
                                                     "held_exponent": {"sha256": "0" * 64}}})
    check(fit2.prepare().rc == 3, "completed c2 fit whose held term is not c1's fine tuning: exit 3")
    _npz(root / "c2/fit/seg01/fit_restart.npz", {"status": "complete", "segments": [{"nfev": 100}],
                                                 "model_extension": {"fit_chain": {"start_restart": {
                                                     "path": str(fitdir / "fit_restart.npz")}},
                                                     "held_exponent": {"sha256": co.sha256_file(held)}}})
    p = fit2.prepare()
    check(p.done, "completed c2 fit with c1's held term: done")


# ------------------------------------------------------------------------------------------------ 3 baseline
BASELINE_CHILD = r'''
import os, sys, numpy as np
sys.path.insert(0, sys.argv[1])
import chain_onejob as co
from pathlib import Path
root = Path(sys.argv[2])
b = co.Baseline([root])
info = b.capture()
from cops_cuda import figures, multi_gpu, gpu_model, reference, evaluator
orig = (figures.read_table, multi_gpu.MultiGpuBackend.activity_linear, evaluator.CudaBackend.__dict__.get("activity_linear"))
path0, env0, argv0 = list(sys.path), dict(os.environ), list(sys.argv)
sys.path.insert(0, str(root))
import stage_mod_x
stage_mod_x.patch()
assert figures.read_table is not orig[0] and multi_gpu.MultiGpuBackend.activity_linear is not orig[1]
assert gpu_model._CLASS == "patched" and reference._bound_root == "patched" and multi_gpu.MultiGpuBackend.extra == 1
rep = b.restore()
ok = [figures.read_table is orig[0], multi_gpu.MultiGpuBackend.activity_linear is orig[1],
      evaluator.CudaBackend.__dict__.get("activity_linear") is orig[2], gpu_model._CLASS is None,
      reference._bound_root is None, not hasattr(multi_gpu.MultiGpuBackend, "extra"), not hasattr(figures, "added"),
      "stage_mod_x" not in sys.modules, sys.path == path0, dict(os.environ) == env0, sys.argv == argv0,
      np.geterr()["divide"] == "warn"]
print("RESULT", ok, rep["purged"], len(rep["restored"]), len(rep["deleted"]))
sys.exit(0 if all(ok) else 1)
'''
STAGE_MOD = r'''
import os, sys, numpy as np
def patch():
    from cops_cuda import figures, multi_gpu, gpu_model, reference, evaluator
    figures.read_table = lambda p: None
    figures.added = 1
    multi_gpu.MultiGpuBackend.activity_linear = lambda self, s: None
    multi_gpu.MultiGpuBackend.extra = 1
    evaluator.CudaBackend.activity_linear = lambda self, s: None
    gpu_model._CLASS = "patched"
    reference._bound_root = "patched"
    os.environ["CHAIN_TEST_X"] = "1"
    os.environ.pop("HOME", None)
    sys.path.insert(0, "/nonexistent")
    sys.argv[:] = ["x", "y"]
    np.seterr(divide="ignore")
'''


def test_baseline(scratch):
    d = scratch / "baseline"
    d.mkdir()
    (d / "stage_mod_x.py").write_text(STAGE_MOD)
    r = subprocess.run([PY, "-W", "ignore", "-c", BASELINE_CHILD, str(HERE), str(d)], stdout=subprocess.PIPE,
                       stderr=subprocess.STDOUT, universal_newlines=True, cwd=str(d))
    print(r.stdout.strip()[-1500:])
    check(r.returncode == 0, "baseline: kept-module patches, environment, sys.path, argv, numpy restored; stage module "
                             "purged")


# ------------------------------------------------------------------------------------------------ 3b real modules
REAL_CHILD = r'''
import os, sys
from pathlib import Path
import numpy as np
HERE = Path(sys.argv[1]); FROZEN = sys.argv[2]
sys.path.insert(0, str(HERE))
os.chdir("/data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda")
import chain_onejob as co
args = co.parse_args(["--chain-root", sys.argv[3], "--lambda", "1e-4", "--ratio", "1e-5", "--source", "/x",
                      "--ft-maxnfev", "3", "--ft-facets", "late_only", "--ft-window-min", "5", "--model-ext", FROZEN,
                      "--devices", "0"])
os.environ[co.KIND_ENV] = "positive"
co._parent_setup_positive(args)
b = co.Baseline([co.HERE, co.APP / "tools", co.REFERENCE_TOOLS])
b.capture()
from cops_cuda import evaluator, multi_gpu, figures, restart, gpu_model, host_model, reference, transport
base = {"cuda_al": evaluator.CudaBackend.__dict__["activity_linear"],
        "multi_al": multi_gpu.MultiGpuBackend.__dict__["activity_linear"],
        "multi_eval": multi_gpu.MultiGpuBackend.__dict__["_evaluate"],
        "upload": evaluator.CudaBackend.__dict__["_upload_activity"],
        "read_table": figures.read_table, "from_fit": restart.from_fit, "handle": transport.WorkerState.handle}
import positive_activity
def stage_fit(tag):
    sys.argv = [str(HERE / "positive_fit.py"), "--ratio", "1e-5", "--", "x"]
    import importlib
    pf = importlib.import_module("positive_fit")
    pp = sys.modules["positive_problem"]
    import grain_fit, grain_crossing
    import gpu_fit
    fc = reference.load("fit_cops")
    ids = {"positive_fit": id(pf), "positive_problem": id(pp), "fit_cops": id(fc), "gpu_fit": id(gpu_fit),
           "held_before": pp.HELD["h"] is None, "STARTED": gpu_fit.STARTED}
    pp.HELD["h"] = np.zeros(4)
    pp.install_held_activity({})
    pp.install_figure_filter(set(), {})
    pp.install_restart_marker({})
    cls = gpu_model.gpu_model_class()
    ids["class_base_is_this_fit_cops_Model"] = issubclass(cls, fc.Model)
    ids["patched"] = (evaluator.CudaBackend.activity_linear is not base["cuda_al"]
                      and multi_gpu.MultiGpuBackend._evaluate is not base["multi_eval"]
                      and figures.read_table is not base["read_table"] and restart.from_fit is not base["from_fit"])
    os.environ["LAMBDA"] = tag
    return ids
def stage_ft():
    import importlib
    sys.argv = [str(HERE / "positive_finetune.py")]
    m = importlib.import_module("positive_finetune")
    m._process_setup()
    chain = co.Chain(args)
    chain.install_finetune_engine(m, chain.stages[1])
    return {"handle_patched": transport.WorkerState.handle is not base["handle"],
            "engine": m.GpuEngine.__name__, "engine_base": m.GpuEngine.__mro__[1].__name__}
def restored():
    return {"cuda_al": evaluator.CudaBackend.__dict__["activity_linear"] is base["cuda_al"],
            "multi_al": multi_gpu.MultiGpuBackend.__dict__["activity_linear"] is base["multi_al"],
            "multi_eval": multi_gpu.MultiGpuBackend.__dict__["_evaluate"] is base["multi_eval"],
            "upload": evaluator.CudaBackend.__dict__["_upload_activity"] is base["upload"],
            "read_table": figures.read_table is base["read_table"], "from_fit": restart.from_fit is base["from_fit"],
            "handle": transport.WorkerState.handle is base["handle"], "gpu_model._CLASS": gpu_model._CLASS is None,
            "host_model._CLASS": host_model._CLASS is None, "bound_root": reference._bound_root is None,
            "purged": not any(n in sys.modules for n in ("positive_fit", "positive_problem", "grain_fit",
                                                          "grain_crossing", "gpu_fit", "fit_cops", "cops_io",
                                                          "positive_finetune")),
            "flags": not hasattr(positive_activity, "_installed_by_positive_fit")
                     and not hasattr(positive_activity, "_installed_by_positive_finetune"),
            "pa_installed": positive_activity._STATE["installed"], "env": "LAMBDA" not in os.environ}
s1 = stage_fit("1"); rep1 = b.restore(); r1 = restored()
s2 = stage_fit("2"); rep2 = b.restore(); r2 = restored()
s3 = stage_ft(); rep3 = b.restore(); r3 = restored()
fresh = all(s1[k] != s2[k] for k in ("positive_fit", "positive_problem", "fit_cops", "gpu_fit")) and s2["held_before"]
ok = (all(r1.values()) and all(r2.values()) and all(r3.values()) and s1["patched"] and s2["patched"] and fresh
      and s1["class_base_is_this_fit_cops_Model"] and s2["class_base_is_this_fit_cops_Model"] and s3["handle_patched"]
      and s3["engine"] == "ChainGpuEngine" and s3["engine_base"] == "GpuEngine" and s2["STARTED"] > s1["STARTED"])
print("STAGE1", s1); print("RESTORED1", r1, len(rep1["purged"]), "purged:", sorted(rep1["purged"]))
print("STAGE2", s2); print("RESTORED2", r2); print("STAGE3", s3); print("RESTORED3", r3, rep3["restored"])
print("REAL", "PASSED" if ok else "FAILED")
sys.exit(0 if ok else 1)
'''


def test_real_modules(scratch):
    env = dict(os.environ, **POS_ENV, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1")
    r = subprocess.run([PY, "-W", "ignore", "-c", REAL_CHILD, str(HERE), str(FROZEN), str(scratch / "real")],
                       stdout=subprocess.PIPE, stderr=subprocess.STDOUT, universal_newlines=True, env=env, timeout=600)
    print(r.stdout.strip()[-4000:])
    check(r.returncode == 0, "real modules: two fit-like stages and a fine-tuning-like stage import positive_fit / "
                             "positive_problem / grain_fit / gpu_fit / fit_cops / positive_finetune afresh, patch the "
                             "held term, figure filter, restart marker, GpuModel class, worker handler, GpuEngine; "
                             "restore() undoes every patch and purges the modules")


def test_spawn(scratch):
    env = dict(os.environ, **POS_ENV, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1")
    r = subprocess.run([PY, "-W", "ignore", str(HERE / "chain_onejob.py"), "--spawn-selftest", "--chain-root",
                        str(scratch / "spawn"), "--lambda", "1e-4", "--ratio", "1e-5", "--source", "/x", "--ft-maxnfev",
                        "3", "--ft-facets", "late_only", "--ft-window-min", "5", "--model-ext", str(FROZEN)],
                       stdout=subprocess.PIPE, stderr=subprocess.STDOUT, universal_newlines=True, env=env, timeout=600,
                       cwd="/data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda")
    print(r.stdout.strip()[-600:])
    check(r.returncode == 0 and "SPAWN SELFTEST PASSED" in r.stdout,
          "spawn: a spawned process (chain_onejob.py as __mp_main__) installs positive_activity as worker, imports "
          "positive_fit and positive_finetune, has the chain_ft_detach request")


# ------------------------------------------------------------------------------------------------ 3c slope kind
SLOPEKIND_CHILD = r'''
import os, sys, json
from pathlib import Path
HERE, FAKE = Path(sys.argv[1]), Path(sys.argv[2])
sys.path.insert(0, str(HERE))
os.chdir("/data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda")
for k in [k for k in os.environ if k.startswith(("COPS_POS_", "COPS_SLOPEMUL_"))]:
    del os.environ[k]
import test_chain_onejob as t
import chain_onejob as co
chain = co.Chain(t.slope_chain_args(FAKE))
os.environ[co.KIND_ENV] = "slope"
setup = chain.parent_setup()
import slope_mult as sm
from cops_cuda import evaluator, multi_gpu, transport
patched = {"cuda_init": evaluator.CudaBackend.__dict__["__init__"], "pipe": transport.pipe_worker_main}
chain.baseline = co.Baseline([co.HERE, co.APP / "tools", co.REFERENCE_TOOLS], keep_names=chain.keep_names)
chain.baseline.capture()
res = {"installed_in_parent": sm._STATE["installed"], "role": sm._STATE.get("role"),
       "dropped": setup["stage_modules_dropped_after_setup"], "kept_with_patch": setup["kept_with_the_patch"],
       "stage_modules_absent_at_baseline": all(m not in sys.modules for m in ("slope_mult_fit", "slope_mult_problem",
                                                                           "clip_problem", "positive_problem"))}
# a cycle-1 stage: its own install_gpu with the same configuration is a no-op
c1 = chain.stages[0]
cfg1 = co._slope_process_config(c1.config_argv())
res["c1_install_again"] = sm.install_gpu(cfg1)
rep = chain.baseline.restore()
res["after_restore_installed"] = sm._STATE["installed"]
res["after_restore_patch_kept"] = (evaluator.CudaBackend.__dict__["__init__"] is patched["cuda_init"]
                                   and transport.pipe_worker_main is patched["pipe"])
res["purged"] = sorted(rep["purged"])
# cycle 2: its configuration (the held term in the frozen config: requirement A) is refused by install_gpu
c2 = chain.stages[2]
cfg2 = co._slope_process_config(c2.config_argv())
res["c2_config_has_held"] = bool(cfg2.get("held"))
try:
    sm.install_gpu(cfg2); res["c2_install"] = "accepted"
except RuntimeError as exc:
    res["c2_install"] = "refused: %s" % exc
chain.baseline.restore()
print("SLOPEKIND", json.dumps(res, default=str))
'''


def test_slope_kind(scratch):
    fake = scratch / "slope_kind"
    for j in (1,):
        (fake / "c1" / "fit").mkdir(parents=True)
        os.symlink(str(REF3A / "c1" / "fit" / "seg01"), str(fake / "c1" / "fit" / "final"))
        (fake / "c1" / "ft").mkdir(parents=True)
        os.symlink(str(REF3A / "c1" / "ft" / "seg01" / "finetune_state.npz"), str(fake / "c1" / "ft" / "finetune_state.npz"))
    env = {k: v for k, v in os.environ.items() if not k.startswith(("COPS_POS_", "COPS_SLOPEMUL_"))}
    env.update(OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1")
    r = subprocess.run([PY, "-W", "ignore", "-c", SLOPEKIND_CHILD, str(HERE), str(fake)], stdout=subprocess.PIPE,
                       stderr=subprocess.STDOUT, universal_newlines=True, env=env, timeout=900)
    line = [l for l in r.stdout.splitlines() if l.startswith("SLOPEKIND ")]
    if not line:
        print(r.stdout[-3000:])
        check(False, "slope kind: the child did not report")
        return
    res = json.loads(line[-1][10:])
    print("slope kind: %s" % json.dumps(res)[:1500])
    check(res["installed_in_parent"] and res["stage_modules_absent_at_baseline"] and res["c1_install_again"] is False
          and res["after_restore_installed"] and res["after_restore_patch_kept"],
          "slope kind: patch installed in the parent before the baseline, stage modules dropped, a stage's install with "
          "the cycle-1 configuration is a no-op, restore keeps the patch")
    check(res["c2_config_has_held"] and str(res["c2_install"]).startswith("refused"),
          "slope kind: requirement A NOT met by the current slope_mult (the cycle-2 held term is in the frozen config; "
          "install_gpu refuses it): a multi-cycle bumps chain needs agent 1's fix (nonneg has it by design)")


# ------------------------------------------------------------------------------------------------ 3d nonneg kind
NNKIND_CHILD = r'''
import os, sys, json
from pathlib import Path
CODE, ROOT = Path(sys.argv[1]), Path(sys.argv[2])
sys.path.insert(0, str(CODE))
os.chdir("/data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda")
import chain_onejob as co
RUN = co.RUN
args = co.parse_args(["--kind", "nonneg", "--slope-model", sys.argv[3], "--chain-root", str(ROOT), "--cycles", "2",
                      "--nfev", "5", "--check", "1", "--ft-maxnfev", "3", "--ft-facets", "late_only",
                      "--ft-window-min", "5", "--l-threads", "32", "--model-ext", str(CODE), "--devices", "0",
                      "--source", str(RUN / "clip_refit" / "D2" / "final" / "fit_restart.npz")]
                     + (["--start-file", sys.argv[4]] if sys.argv[3] == "squared" else [])
                     + (["--slope-mult-min", sys.argv[5]] if len(sys.argv) > 5 and sys.argv[5] != "-" else []))
chain = co.Chain(args)
os.environ[co.KIND_ENV] = "nonneg"
setup = chain.parent_setup()
import nonneg_mult as nm
from cops_cuda import evaluator, multi_gpu, transport
patched = {"pipe": transport.pipe_worker_main, "multi": dict(vars(multi_gpu.MultiGpuBackend))}
chain.baseline = co.Baseline([co.HERE, co.APP / "tools", co.REFERENCE_TOOLS], keep_names=chain.keep_names)
chain.baseline.capture()
res = {"installed": nm._STATE["installed"], "keep_names": list(chain.keep_names),
       "dropped": setup["stage_modules_dropped_after_setup"], "kept_with_patch": setup["kept_with_the_patch"],
       "threads": nm._STATE.get("threads"), "frozen_env": bool(os.environ.get(nm.ENV_FROZEN))}
import nonneg_fit
c1, c2 = chain.stages[0], chain.stages[2]
cfg1 = nonneg_fit.process_config_from_argv(c1.config_argv())
res["c1_install_again"] = nm.install_gpu(cfg1)
res["L_threads_again"] = nm.set_L_threads(32, "chunked").get("reused")
rep = chain.baseline.restore()
res["after_restore_installed"] = nm._STATE["installed"]
res["after_restore_patch_kept"] = transport.pipe_worker_main is patched["pipe"]
res["slope_mult_kept"] = "slope_mult" in sys.modules
res["purged"] = sorted(rep["purged"])
# cycle 2: its own configuration (with --held) is the same constant configuration (requirement A)
import nonneg_fit as nf2
cfg2 = nf2.process_config_from_argv(c2.config_argv())
res["c2_config_equal"] = nm.canonical(cfg2) == nm.canonical(cfg1)
import nonneg_problem as npb
res["slope_mult_min_in_c1_argv"] = "--slope-mult-min" in c1.config_argv()
res["slope_mult_min_in_c2_argv"] = "--slope-mult-min" in c2.config_argv()
res["rewrite_c1"] = npb.rewrite_arguments(nf2._split_threads(c1.config_argv())[0])[0].get("slope_mult_min")
res["rewrite_c2"] = npb.rewrite_arguments(nf2._split_threads(c2.config_argv())[0])[0].get("slope_mult_min")
try:
    res["c2_install"] = nm.install_gpu(cfg2)
except Exception as exc:
    res["c2_install"] = "refused: %s" % exc
rec = nm.set_stage_held(str(c2.held), "held_next")
res["held_set"] = rec and rec.get("sha256")
nm.set_stage_held(None)
chain.baseline.restore()
print("NNKIND", json.dumps(res, default=str))
'''


NN2_CHILD = r'''
import os, sys, json
from pathlib import Path
CODE, ROOT, MODEL, C1, TOK = Path(sys.argv[1]), Path(sys.argv[2]), sys.argv[3], Path(sys.argv[4]), json.loads(sys.argv[5])
sys.path.insert(0, str(CODE))
os.chdir("/data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda")
import chain_onejob as co
toks = []
for k, v in TOK.items():
    toks += ["--" + k.lower().replace("_", "-"), v]
args = co.parse_args(["--kind", "nonneg", "--slope-model", MODEL, "--chain-root", str(ROOT), "--cycles", "3",
                      "--nfev", "5", "--check", "1", "--ft-maxnfev", "3", "--ft-facets", "late_only",
                      "--ft-window-min", "5", "--l-threads", "32", "--model-ext", str(CODE), "--devices", "0",
                      "--source", str((C1 / "final" / "fit_restart.npz").resolve()), "--slope-mult-min", "1",
                      "--c1-from", str(C1)] + toks)
chain = co.Chain(args)
os.environ[co.KIND_ENV] = "nonneg"
setup = chain.parent_setup()                     # the first fit stage is c2 (c1 adopted): --start restart, no held
import nonneg_mult as nm
from cops_cuda import transport
pipe = transport.pipe_worker_main
chain.baseline = co.Baseline([co.HERE, co.APP / "tools", co.REFERENCE_TOOLS], keep_names=chain.keep_names)
chain.baseline.capture()
res = {"installed": nm._STATE["installed"], "keep_names": list(chain.keep_names),
       "kept_with_patch": setup["kept_with_the_patch"], "frozen_env": bool(os.environ.get(nm.ENV_FROZEN)),
       "threads": nm._STATE.get("threads"), "layout": [(s.name, s.kind) for s in chain.stages]}
import nonneg_fit
c2 = chain.stages[2]
cfg0 = nonneg_fit.process_config_from_argv(c2.config_argv())
res["parent_config_model"] = cfg0.get("model")
res["parent_config_sos_rank"] = cfg0.get("sos_rank")
res["stage_install_again"] = str(nm.install_gpu(cfg0))
rep = chain.baseline.restore()
res["after_restore_installed"] = nm._STATE["installed"]
res["after_restore_patch_kept"] = transport.pipe_worker_main is pipe
res["purged"] = sorted(rep["purged"])
# cycles 2 and 3 (start chain, with --held): the constant configuration is the parent's
import nonneg_fit as nf2
import nonneg_problem as npb
for st in (chain.stages[2], chain.stages[4]):
    argv = st._argv(1, 5, st._start_options(1, None), st.directory / "seg01")[1:]
    cfg = nf2.process_config_from_argv(argv)
    rw = npb.rewrite_arguments(nf2._split_threads(argv)[0])[0]
    res[st.name] = {"config_equal": nm.canonical(cfg) == nm.canonical(cfg0),
                    "rewrite": {k: rw.get(k) for k in ("model", "start_choice", "slope_mult_min", "sos_rank", "bg_min",
                                                       "spread_f", "spread_cref", "spread_u0")},
                    "has_held": "--held" in argv}
    try:
        res[st.name]["install"] = str(nm.install_gpu(cfg))
    except Exception as exc:
        res[st.name]["install"] = "refused: %s" % exc
rec = nm.set_stage_held(str(chain.stages[2].held), "held_next")
res["held_set"] = rec and rec.get("sha256")
nm.set_stage_held(None)
chain.baseline.restore()
print("NN2", json.dumps(res, default=str))
'''


NN2C1_CHILD = r'''
import os, sys, json
from pathlib import Path
CODE, ROOT, MODEL, C1, TOK = Path(sys.argv[1]), Path(sys.argv[2]), sys.argv[3], Path(sys.argv[4]), json.loads(sys.argv[5])
sys.path.insert(0, str(CODE))
os.chdir("/data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda")
import chain_onejob as co
toks = []
for k, v in TOK.items():
    toks += ["--" + k.lower().replace("_", "-"), v]
args = co.parse_args(["--kind", "nonneg", "--slope-model", MODEL, "--chain-root", str(ROOT), "--cycles", "1",
                      "--nfev", "5", "--check", "1", "--ft-maxnfev", "3", "--ft-facets", "late_only",
                      "--ft-window-min", "5", "--l-threads", "32", "--model-ext", str(CODE), "--devices", "0",
                      "--source", str((C1 / "final" / "fit_restart.npz").resolve()), "--slope-mult-min", "1",
                      "--c1-from", str(C1)] + toks)
chain = co.Chain(args)
os.environ[co.KIND_ENV] = "nonneg"
setup = chain.parent_setup()
import nonneg_mult as nm
print("NN2C1", json.dumps({"installed": nm._STATE["installed"], "threads": nm._STATE.get("threads"),
                           "stages": [s.name for s in chain.stages], "config_stage": co.nonneg_config_stage(chain).kind,
                           "keep": list(chain.keep_names)}, default=str))
'''


def nonneg_stage2_kind(code_dir, scratch, ref3a_ft_state):
    """CPU job (3d-style): the nonneg process kind for the stage-2 chains (c1 adopted) of sos, background, logspread
    with the real snapshot code: parent setup from c2's tokens (--start restart, no held), the stage installs."""
    for model, (tokens, c1) in STAGE2.items():
        root = Path(scratch) / ("nn2_" + model)
        (root / "c1" / "fit").mkdir(parents=True)
        os.symlink(str((c1 / "final").resolve()), str(root / "c1" / "fit" / "final"))
        for k in (1, 2):
            (root / ("c%d" % k) / "ft").mkdir(parents=True)
            os.symlink(str(ref3a_ft_state), str(root / ("c%d" % k) / "ft" / "finetune_state.npz"))
        (root / "c2" / "fit" / "seg01").mkdir(parents=True)
        os.symlink("seg01", str(root / "c2" / "fit" / "final"))
        os.symlink(str((c1 / "final" / "fit_restart.npz").resolve()), str(root / "c2" / "fit" / "seg01" / "fit_restart.npz"))
        env = {k: v for k, v in os.environ.items() if not k.startswith(("COPS_POS_", "COPS_SLOPEMUL_", "COPS_NONNEG_"))}
        env.update(OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1")
        r = subprocess.run([PY, "-W", "ignore", "-c", NN2_CHILD, str(code_dir), str(root), model, str(c1),
                            json.dumps(tokens)], stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                           universal_newlines=True, env=env, timeout=2400)
        line = [l for l in r.stdout.splitlines() if l.startswith("NN2 ")]
        if not line:
            print(r.stdout[-4000:])
            check(False, "nonneg stage-2 kind %s: the child did not report (exit %d)" % (model, r.returncode))
            continue
        res = json.loads(line[-1][4:])
        print("nonneg stage-2 kind %s: %s" % (model, json.dumps(res)[:2400]))
        fl = {k.lower(): float(v) for k, v in tokens.items()}
        check(res["installed"] and res["frozen_env"] and "slope_mult" in res["keep_names"]
              and res["after_restore_installed"] and res["after_restore_patch_kept"]
              and "nonneg_fit" in res["purged"] and "nonneg_problem" in res["purged"]
              and res["parent_config_model"] == model and res["threads"] == 32,
              "nonneg stage-2 kind %s: nonneg_mult installed in the parent from c2's tokens (c1 adopted), 32 L threads, "
              "restore keeps the patch, stage modules purged" % model)
        for name in ("c2fit", "c3fit"):
            r_ = res.get(name) or {}
            rw = r_.get("rewrite") or {}
            ok_tok = all(rw.get(k) is not None and float(rw[k]) == v for k, v in fl.items())
            check(r_.get("config_equal") and r_.get("has_held") and not str(r_.get("install")).startswith("refused")
                  and rw.get("start_choice") == "chain" and rw.get("slope_mult_min") == 1.0 and ok_tok,
                  "nonneg stage-2 kind %s %s: --start chain with --held, the parent's configuration, install accepted, "
                  "tokens %s and --slope-mult-min 1 parsed by nonneg_problem" % (model, name, rw))
        check(bool(res.get("held_set")), "nonneg stage-2 kind %s: set_stage_held accepts the held term" % model)
        # --cycles 1 with c1 adopted (c1 ft only): the parent setup from the adopted stage (v7.1 raised StopIteration)
        root1 = Path(scratch) / ("nn2c1_" + model)
        root1.mkdir()
        r = subprocess.run([PY, "-W", "ignore", "-c", NN2C1_CHILD, str(code_dir), str(root1), model, str(c1),
                            json.dumps(tokens)], stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                           universal_newlines=True, env=env, timeout=2400)
        line = [l for l in r.stdout.splitlines() if l.startswith("NN2C1 ")]
        if not line:
            print(r.stdout[-4000:])
            check(False, "nonneg stage-2 kind %s --cycles 1: the child did not report (exit %d)" % (model, r.returncode))
            continue
        res1 = json.loads(line[-1][6:])
        print("nonneg stage-2 kind %s --cycles 1: %s" % (model, json.dumps(res1)))
        check(res1["installed"] and res1["threads"] == 32 and res1["stages"] == ["c1fit", "c1ft"]
              and res1["config_stage"] == "nonneg_fit_adopted" and "slope_mult" in res1["keep"],
              "nonneg stage-2 kind %s --cycles 1 (c1 adopted, c1 ft only): parent setup from the adopted stage, "
              "nonneg_mult installed" % model)


def nonneg_kind(code_dir, scratch, ref3a_ft_state):
    """CPU job: the nonneg process kind with the real snapshot code (both models)."""
    out = {}
    for model, smm in (("softplus", None), ("squared", None), ("softplus", "1")):
        tag = model + ("_smm" + smm if smm else "")
        root = Path(scratch) / ("nnkind_" + tag)
        (root / "c1" / "fit" / "seg01").mkdir(parents=True)
        (root / "c1" / "fit" / "seg01" / "fit_restart.npz").write_text("x")
        os.symlink("seg01", str(root / "c1" / "fit" / "final"))
        (root / "c1" / "ft").mkdir(parents=True)
        os.symlink(str(ref3a_ft_state), str(root / "c1" / "ft" / "finetune_state.npz"))
        env = {k: v for k, v in os.environ.items() if not k.startswith(("COPS_POS_", "COPS_SLOPEMUL_", "COPS_NONNEG_"))}
        env.update(OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1")
        r = subprocess.run([PY, "-W", "ignore", "-c", NNKIND_CHILD, str(code_dir), str(root), model,
                            str(NONNEG_TEST_START), smm or "-"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                           universal_newlines=True, env=env, timeout=2400)
        line = [l for l in r.stdout.splitlines() if l.startswith("NNKIND ")]
        if not line:
            print(r.stdout[-4000:])
            check(False, "nonneg kind %s: the child did not report (exit %d)" % (model, r.returncode))
            continue
        res = json.loads(line[-1][7:])
        out[tag] = res
        print("nonneg kind %s: %s" % (tag, json.dumps(res)[:1600]))
        if smm:
            check(res["slope_mult_min_in_c1_argv"] and res["slope_mult_min_in_c2_argv"]
                  and res["rewrite_c1"] == float(smm) and res["rewrite_c2"] == float(smm),
                  "nonneg kind %s: --slope-mult-min %s in the cycle-1 and cycle-2 fit arguments, accepted by "
                  "nonneg_problem.rewrite_arguments for both (value %s / %s)" % (tag, smm, res["rewrite_c1"],
                                                                                res["rewrite_c2"]))
        model_ = model
        check(res["installed"] and res["frozen_env"] and "slope_mult" in res["keep_names"]
              and res["c1_install_again"] in (False, True) and res["after_restore_installed"]
              and res["after_restore_patch_kept"] and res["slope_mult_kept"]
              and "nonneg_fit" in res["purged"] and "nonneg_problem" in res["purged"],
              "nonneg kind %s: nonneg_mult installed in the parent (frozen env set), slope_mult kept with it, a stage's "
              "install with the cycle-1 configuration accepted, restore keeps the patch, stage modules purged" % tag)
        check(res["c2_config_equal"] and not str(res["c2_install"]).startswith("refused") and res["held_set"],
              "nonneg kind %s: requirement A met (the cycle-2 configuration with --held equals cycle 1's; the held term "
              "is set per stage)" % tag)
    return out


# ------------------------------------------------------------------------------------------------ 4 tee
TEE_CHILD = r'''
import os, subprocess, sys
sys.path.insert(0, sys.argv[1])
import chain_onejob as co
t = co.Tee()
print("before stage")
t.switch(sys.argv[2])
print("inside stage (python)")
sys.stderr.write("inside stage (stderr)\n")
subprocess.run(["/bin/echo", "inside stage (child process)"])
t.switch(None)
print("after stage")
t.close()
print("after close")
'''


def test_tee(scratch):
    d = scratch / "tee"
    d.mkdir()
    job = d / "job.log"
    stage = d / "stage.log"
    with job.open("w") as f:
        r = subprocess.run([PY, "-W", "ignore", "-c", TEE_CHILD, str(HERE), str(stage)], stdout=f,
                           stderr=subprocess.STDOUT, cwd=str(d), timeout=120)
    j, s = job.read_text(), stage.read_text() if stage.exists() else ""
    want_s = ["inside stage (python)", "inside stage (stderr)", "inside stage (child process)"]
    check(r.returncode == 0 and all(w in j for w in want_s + ["before stage", "after stage", "after close"])
          and all(w in s for w in want_s) and "before stage" not in s and "after stage" not in s
          and "\x1e" not in j and "\x1e" not in s, "tee: job log has everything, stage log only the stage (job %r, "
                                                    "stage %r)" % (j[-300:], s[-300:]))


# ------------------------------------------------------------------------------------------------ 5 real fit models
REALFIT_CHILD = r'''
import os, sys, types, copy, json
from pathlib import Path
import numpy as np
CODE, SCRATCH, SOURCE, FROZEN = Path(sys.argv[1]), Path(sys.argv[2]), sys.argv[3], sys.argv[4]
sys.path.insert(0, str(CODE))
os.chdir("/data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda")
import chain_onejob as co
print("chain_onejob from", co.__file__, flush=True)
args = co.parse_args(["--chain-root", str(SCRATCH / "chain"), "--lambda", "1e-4", "--ratio", "1e-5", "--cycles", "1",
                      "--source", SOURCE, "--nfev", "5", "--check", "1", "--ft-maxnfev", "3", "--ft-facets",
                      "late_only", "--ft-window-min", "5", "--model-ext", FROZEN, "--devices", "0,1"])
chain = co.Chain(args)
os.environ[co.KIND_ENV] = "positive"
co.PROCESS_KINDS["positive"]["parent"](args)
chain.baseline = co.Baseline([co.HERE, co.APP / "tools", co.REFERENCE_TOOLS])
chain.baseline.capture()
chain.install_registry()
from cops_cuda import gpu_model, manifest
manifest.backend_manifest = lambda kind, device=0: {"stub": "test_chain_onejob realfit (no GPU)"}
SEEN = []
class StopTest(Exception):
    pass
class FakeBackend:
    def __init__(self, model, basis):
        self.rows = list(model.rows); self.n = len(self.rows); self.nbasis_enhancement = 0 if basis is None else basis.shape[1]
        self.enhancement_basis, self.enhancement_scale, self.enhancement_basis_sha256 = basis, 1.0, None
        self.closed = self.broken = self.released = False; self.block_records = [0]
    def evaluate(self, state, request):
        cols = {k: (np.zeros(self.n), np.zeros(self.n)) for k in request.keys()}
        return types.SimpleNamespace(pNG=np.full(self.n, 1e-9), pRG=np.full(self.n, 1e-9), columns=cols)
    def summary(self): return "FakeBackend(%d rows)" % self.n
    def identity(self): return {"fake": True}
    def release(self): self.closed = True
    close = release
def fake_backend_for(self, model, basis, scale):
    key = (str(model.run), tuple(int(r) for r in model.rows), tuple(self.devices), self.mode)
    if key not in self.backends:
        self.backends[key] = (FakeBackend(model, basis), basis)
    be = self.backends[key][0]
    if basis is not None:
        be.nbasis_enhancement = basis.shape[1]
    if all(m is not model for m in SEEN):
        SEEN.append(model)
        print("FAKE: model %d (%s) reaches the registry; header %d lines" % (len(SEEN), type(model).__name__,
              len(model._spec["header"])), flush=True)
    if len(SEEN) >= 2:
        raise StopTest("the second model passed the chain's verification")
    return be
type(gpu_model.REGISTRY).backend_for = fake_backend_for
stage = chain.stages[0]
chain.current = stage
plan = stage.prepare()
assert plan.rc is None, plan.notes
os.environ.update(plan.env)
rc = stage.call(plan)
print("stage exit %d; models that reached the fake registry: %d" % (rc, len(SEEN)), flush=True)
out = {"code": str(co.__file__), "models_reaching_registry": len(SEEN)}
if len(SEEN) >= 1:
    MF = SEEN[0]
    out["MF_header_lines"] = len(MF._spec["header"])
if len(SEEN) >= 2:
    MF, M = SEEN[0], SEEN[1]
    a, b = MF._spec, M._spec
    cols = a["cols"]
    names = {i: n for n, i in cols.items()}
    diffcols = sorted({names[j] for ra, rb in zip(a["rows"], b["rows"]) for j, (x, y) in enumerate(zip(ra, rb)) if x != y})
    extra = [h for h in b["header"] if h not in a["header"]]
    out.update(M_header_lines=len(b["header"]), differing_columns=diffcols, header_lines_only_in_M=extra,
               rows_equal=[int(r) for r in MF.rows] == [int(r) for r in M.rows],
               run_equal=a["run"] == b["run"], coef_equal=bool(np.array_equal(a["coef"], b["coef"])))
    # v4 only: identity checks on model-like variants of M
    if hasattr(co, "dataset_identity"):
        attrs = ("run", "rows", "ipoint", "Y", "indptr", "rh", "VdotB", "TRG", "bm", "cNG", "cRG", "SZAMode",
                 "helio", "TNG", "RGModel", "S", "T0")
        class Variant:                       # weak-referenceable, like a real model (job 124089 used SimpleNamespace)
            pass
        def variant(change):
            v = Variant()
            for k in attrs:
                if hasattr(M, k):
                    setattr(v, k, getattr(M, k))
            v._spec = copy.deepcopy(M._spec)
            change(v)
            return v
        kept = []                            # keep every variant alive: no id() can be reused between them
        chain.current = types.SimpleNamespace(name="c2fit")
        res = {}
        def tryv(name, change):
            try:
                v = variant(change); kept.append(v)
                chain.verify_model(v); res[name] = "accepted"
            except RuntimeError as exc:
                res[name] = "refused: %s" % str(exc)[-160:]
        tryv("same as M (another stage)", lambda v: None)
        tryv("useNG of row 0 flipped", lambda v: v._spec["rows"][0].__setitem__(cols["useNG"], "1" if v._spec["rows"][0][cols["useNG"]] == "0" else "0"))
        tryv("another run directory", lambda v: v._spec.__setitem__("run", "/data001/highresolu_vt/cops_output/runs/other/base"))
        tryv("T_RG of row 0 changed", lambda v: v._spec["rows"][0].__setitem__(cols["T_RG"], "999"))
        tryv("x of row 5 changed", lambda v: v._spec["rows"][5].__setitem__(cols["x"], "1.0"))
        tryv("another Y", lambda v: setattr(v, "Y", v.Y[:-1]))
        tryv("other rows", lambda v: setattr(v, "rows", v.rows[:-1]))
        tryv("another TemperatureNG", lambda v: setattr(v, "TNG", v.TNG + 1.0))
        # the same refusals with every variant freed right after its check (an id() can be reused: job 124089, v4)
        def tryfreed(name, change):
            try:
                chain.verify_model(variant(change)); res["freed: " + name] = "accepted"
            except RuntimeError as exc:
                res["freed: " + name] = "refused: %s" % str(exc)[-160:]
        tryfreed("same as M (another stage)", lambda v: None)
        for name, change in (("another run directory", lambda v: v._spec.__setitem__("run", "/data001/highresolu_vt/cops_output/runs/other/base")),
                             ("T_RG of row 0 changed", lambda v: v._spec["rows"][0].__setitem__(cols["T_RG"], "999")),
                             ("x of row 5 changed", lambda v: v._spec["rows"][5].__setitem__(cols["x"], "1.0")),
                             ("another Y", lambda v: setattr(v, "Y", v.Y[:-1])),
                             ("other rows", lambda v: setattr(v, "rows", v.rows[:-1])),
                             ("another TemperatureNG", lambda v: setattr(v, "TNG", v.TNG + 1.0))):
            tryfreed(name, change)
        out["variants"] = res
print("REALFIT", json.dumps(out, default=str), flush=True)
'''


def realfit(code_dir, scratch):
    """CPU job (12M base): the real cycle-1 fit stage through fit_cops.main() up to GrainProblem's M.evaluate, with a
    fake GPU backend; reports whether the second model (fit_cops' M) passes the chain's verification (v3: refused, as
    job 123990; v4: accepted) and how MF and M differ; v4: refusals of other bases / records."""
    scratch = Path(scratch)
    copy_dir = scratch / ("code_" + Path(code_dir).resolve().parent.name)
    shutil.copytree(str(code_dir), str(copy_dir), ignore=shutil.ignore_patterns("__pycache__"))
    env = dict(os.environ, **POS_ENV, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", PYTHONUNBUFFERED="1")
    r = subprocess.run([PY, "-W", "ignore", "-c", REALFIT_CHILD, str(copy_dir), str(scratch / copy_dir.name),
                        str(REF_SOURCE), str(FROZEN)], stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                       universal_newlines=True, env=env, timeout=3000)
    keep = [l for l in r.stdout.splitlines() if l.startswith(("chain_onejob", "FAKE", "stage exit", "REALFIT",
                                                               "RuntimeError", "positive_fit: profile"))]
    print("\n".join(l[:600] for l in keep))
    line = [l for l in r.stdout.splitlines() if l.startswith("REALFIT ")]
    return json.loads(line[-1][8:]) if line else None


# ------------------------------------------------------------------------------------------------ compare
NUM = re.compile(r"[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[eE][-+]?\d+)?|nan|inf", re.I)


def numbers(text):
    return [t for t in NUM.findall(text)]


STAMP = re.compile(r"\d{4}-\d\d-\d\dT\d\d:\d\d:\d\d(?:\.\d+)?(?:\+00:00|Z)?")


def text_numbers_diff(a, b, skip_regex=None):
    """Data lines (not starting with '#'): exact text; '#' lines and result files: the numbers, with creation time
    stamps (ISO, e.g. 2026-11-09T06:10:54+00:00) masked. Returns (n lines, list of differing (line no, a, b))."""
    la, lb = a.read_text(errors="replace").splitlines(), b.read_text(errors="replace").splitlines()
    out = []
    if len(la) != len(lb):
        out.append((-1, "%d lines" % len(la), "%d lines" % len(lb)))
    for k, (x, y) in enumerate(zip(la, lb)):
        if skip_regex and (re.search(skip_regex, x) or re.search(skip_regex, y)):
            continue
        if a.suffix == ".dat" and not x.startswith("#"):
            same = x == y
        else:
            same = numbers(STAMP.sub("<utc>", x)) == numbers(STAMP.sub("<utc>", y))
        if not same:
            out.append((k + 1, x[:200], y[:200]))
    return len(la), out


def npz_diff(a, b, skip=("meta_json",)):
    res = {}
    with np.load(str(a), allow_pickle=False) as za, np.load(str(b), allow_pickle=False) as zb:
        keys = sorted(set(za.files) | set(zb.files))
        for k in keys:
            if k in skip:
                continue
            if k not in za.files or k not in zb.files:
                res[k] = "missing in %s" % ("A" if k not in za.files else "B")
                continue
            x, y = za[k], zb[k]
            if x.dtype != y.dtype or x.shape != y.shape:
                res[k] = "dtype/shape %s %s vs %s %s" % (x.dtype, x.shape, y.dtype, y.shape)
                continue
            if x.dtype.kind in "fc":
                same = x.tobytes() == y.tobytes()
                if not same:
                    d = np.abs(x.astype(np.float64) - y.astype(np.float64))
                    sc = np.maximum(np.abs(x), np.abs(y)).astype(np.float64)
                    with np.errstate(divide="ignore", invalid="ignore"):
                        relv = np.where(sc > 0, d / sc, 0.0)
                    res[k] = "DIFF: %d of %d values differ, max abs %.3e, max rel %.3e" % (
                        int(np.sum(x != y)), x.size, float(np.nanmax(d)), float(np.nanmax(relv)))
                else:
                    res[k] = "bitwise"
            else:
                res[k] = "bitwise" if np.array_equal(x, y) else "DIFF (%s)" % x.dtype
        ma = json.loads(str(za["meta_json"])) if "meta_json" in za.files else {}
        mb = json.loads(str(zb["meta_json"])) if "meta_json" in zb.files else {}
    return res, ma, mb


def compare(A, B, cycles=None, out_json=None):
    A, B = Path(A).resolve(), Path(B).resolve()
    report = {"A": str(A), "B": str(B), "stages": {}}
    ks = sorted(int(p.name[1:]) for p in A.glob("c[0-9]") if p.is_dir())
    if cycles:
        ks = [k for k in ks if k <= cycles]
    all_bitwise = True
    for k in ks:
        for kind in ("fit", "ft"):
            name = "c%d%s" % (k, kind)
            da, db = A / ("c%d" % k) / kind / "final", B / ("c%d" % k) / kind / "final"
            rec = {"A": str(da.resolve()) if da.exists() else None, "B": str(db.resolve()) if db.exists() else None}
            report["stages"][name] = rec
            if not (da.exists() and db.exists()):
                rec["status"] = "missing final"
                all_bitwise = False
                continue
            if kind == "fit":
                n, diffs = text_numbers_diff(da / "fit_result.txt", db / "fit_result.txt",
                                             skip_regex=r"(/data001/|written|UTC|utc|seconds| s wall| s\b|time)")
                rec["fit_result.txt"] = {"lines": n, "differing_lines": diffs[:20], "n_differing": len(diffs)}
                arr, ma, mb = npz_diff(da / "fit_restart.npz", db / "fit_restart.npz")
                rec["fit_restart.npz"] = arr
                rec["fit_restart.meta"] = {
                    "optimizer": (ma.get("optimizer") == mb.get("optimizer"), ma.get("optimizer"), mb.get("optimizer")),
                    "status": (ma.get("status"), mb.get("status")),
                    "segments_nfev_cost": ([(s.get("nfev"), s.get("cost")) for s in ma.get("segments", [])],
                                           [(s.get("nfev"), s.get("cost")) for s in mb.get("segments", [])])}
                arr2, _, _ = npz_diff(da / "activity_exp_facets.npz", db / "activity_exp_facets.npz")
                rec["activity_exp_facets.npz"] = arr2
                files = ("COPS_model.dat", "late_window_model.dat")
            else:
                arr, ma, mb = npz_diff(da / "finetune_state.npz", db / "finetune_state.npz")
                rec["finetune_state.npz"] = arr
                rec["finetune_state.meta"] = {"costs": (ma.get("costs") == mb.get("costs"), ma.get("costs"),
                                                        mb.get("costs")),
                                              "stop/nfev": ((ma.get("stop"), ma.get("nfev_total")),
                                                            (mb.get("stop"), mb.get("nfev_total")))}
                ja = json.loads((da / "fine_tuning.json").read_text())
                jb = json.loads((db / "fine_tuning.json").read_text())
                ha = [h["cost"] for h in ja.get("optimizer_history", [])]
                hb = [h["cost"] for h in jb.get("optimizer_history", [])]
                rec["cost_history"] = {"bitwise": ha == hb, "A": ha, "B": hb}
                rec["grain"] = {"A": (ja.get("final") or {}).get("grain_A_uPa"),
                                "B": (jb.get("final") or {}).get("grain_A_uPa"),
                                "A_B": (ja.get("final") or {}).get("grain_B_uPa_per_100s"),
                                "B_B": (jb.get("final") or {}).get("grain_B_uPa_per_100s")}
                rec["check1"] = {"A": (ja.get("start_check") or {}).get("relative_difference_data_plus_reg_vs_fit_cost"),
                                 "B": (jb.get("start_check") or {}).get("relative_difference_data_plus_reg_vs_fit_cost")}
                n, diffs = text_numbers_diff(da / "finetune_result.txt", db / "finetune_result.txt",
                                             skip_regex=r"(/data001/|UTC|utc|seconds| s\b|written|time)")
                rec["finetune_result.txt"] = {"lines": n, "differing_lines": diffs[:20], "n_differing": len(diffs)}
                files = ("COPS_model.dat", "late_window_model.dat")
            for f in files:
                n, diffs = text_numbers_diff(da / f, db / f)
                rec[f] = {"lines": n, "n_differing": len(diffs), "differing_lines": diffs[:10]}
            bit = all(v == "bitwise" for key in ("fit_restart.npz", "finetune_state.npz", "activity_exp_facets.npz")
                      for v in (rec.get(key) or {}).values())
            bit = bit and all(rec[f]["n_differing"] == 0 for f in files)
            if kind == "fit":
                bit = bit and rec["fit_result.txt"]["n_differing"] == 0 and rec["fit_restart.meta"]["optimizer"][0]
            else:
                bit = bit and rec["cost_history"]["bitwise"] and rec["finetune_state.meta"]["costs"][0] and \
                    rec["finetune_result.txt"]["n_differing"] == 0
            rec["status"] = "BITWISE IDENTICAL (checked)" if bit else "DIFFERENT"
            all_bitwise = all_bitwise and bit
            print("%-6s %s" % (name, rec["status"]))
            for key in ("fit_restart.npz", "finetune_state.npz", "activity_exp_facets.npz"):
                if key in rec:
                    bad = {k: v for k, v in rec[key].items() if v != "bitwise"}
                    print("    %s: %d arrays, %s" % (key, len(rec[key]), "all bitwise" if not bad else bad))
            for f in ("fit_result.txt", "finetune_result.txt") + files:
                if f in rec:
                    print("    %s: %d lines, %d with different numbers%s" % (
                        f, rec[f]["lines"], rec[f]["n_differing"],
                        "" if not rec[f]["n_differing"] else "; first %s" % (rec[f]["differing_lines"][:2],)))
            if kind == "ft":
                print("    cost history: %d / %d evaluations, bitwise %s; grain A %r / %r, B %r / %r" % (
                    len(rec["cost_history"]["A"]), len(rec["cost_history"]["B"]), rec["cost_history"]["bitwise"],
                    rec["grain"]["A"], rec["grain"]["B"], rec["grain"]["A_B"], rec["grain"]["B_B"]))
                print("    check 1 (rel. difference to the fit's cost): %r / %r" % (rec["check1"]["A"],
                                                                                    rec["check1"]["B"]))
            else:
                print("    optimizer meta equal %s; segments A %s B %s" % (
                    rec["fit_restart.meta"]["optimizer"][0], rec["fit_restart.meta"]["segments_nfev_cost"][0],
                    rec["fit_restart.meta"]["segments_nfev_cost"][1]))
    report["all_bitwise"] = all_bitwise
    print("ALL STAGES BITWISE IDENTICAL (checked)" if all_bitwise else "NOT all stages bitwise identical")
    if out_json:
        Path(out_json).write_text(json.dumps(report, indent=1, default=str) + "\n")
    return 0 if all_bitwise else 1



TIME_KEY = re.compile(r"(second|_s$|time|utc|wall|elapsed|pid|host|path|command|created|job|provenance|dir|^out$|"
                      r"file)", re.I)                       # paths of the run directories, timings, processes


def strip_times(obj):
    """A record without timing / provenance fields (they differ between any two runs)."""
    if isinstance(obj, dict):
        return {k: strip_times(v) for k, v in obj.items() if not TIME_KEY.search(str(k))}
    if isinstance(obj, list):
        return [strip_times(v) for v in obj]
    return obj


META_KEYS = ("optimizer", "costs", "costs_initial", "costs_final", "stop", "complete", "nfev_total", "parameters",
             "status")


def compare_stage_generic(da, db):
    """Every *.npz (arrays bitwise; meta_json META_KEYS), every *.dat (data lines exact, '#' lines numbers without
    creation stamps), fit_result.txt / finetune_result.txt (numbers), fine_tuning.json (optimizer record, history,
    costs, timing fields stripped)."""
    rec, bit = {}, True
    names = sorted(p.name for p in da.iterdir() if p.is_file())
    for n in names:
        fb = db / n
        if not fb.is_file():
            if n.endswith((".npz", ".dat")) or n in ("fit_result.txt", "finetune_result.txt", "fine_tuning.json"):
                rec[n] = "missing in B"
                bit = False
            continue
        if n.endswith(".npz"):
            arr, ma, mb = npz_diff(da / n, fb)
            meta = {k: (ma.get(k) == mb.get(k)) for k in META_KEYS if k in ma or k in mb}
            rec[n] = {"arrays": arr, "meta_equal": meta}
            bit = bit and all(v == "bitwise" for v in arr.values()) and all(meta.values())
        elif n.endswith(".dat"):
            nl, diffs = text_numbers_diff(da / n, fb)
            rec[n] = {"lines": nl, "n_differing": len(diffs), "differing_lines": diffs[:5]}
            bit = bit and not diffs
        elif n in ("fit_result.txt", "finetune_result.txt"):
            nl, diffs = text_numbers_diff(da / n, fb, skip_regex=r"(/data001/|UTC|utc|seconds| s\b|written|time)")
            rec[n] = {"lines": nl, "n_differing": len(diffs), "differing_lines": diffs[:5]}
            bit = bit and not diffs
        elif n == "fine_tuning.json":
            ja, jb = json.loads((da / n).read_text()), json.loads(fb.read_text())
            sub = {}
            for k in ("optimizer", "optimizer_history", "costs_initial", "costs_final", "initial", "final", "seen",
                      "segment", "stop", "bounds", "start_check", "validation"):
                if k in ja or k in jb:
                    sub[k] = strip_times(ja.get(k)) == strip_times(jb.get(k))
            rec[n] = sub
            bit = bit and all(sub.values())
    rec["status"] = "BITWISE IDENTICAL (checked)" if bit else "DIFFERENT"
    return rec, bit


def compare_generic(A, B, out_json=None):
    A, B = Path(A).resolve(), Path(B).resolve()
    report, allbit = {"A": str(A), "B": str(B), "stages": {}}, True
    for k in sorted(int(p.name[1:]) for p in A.glob("c[0-9]") if p.is_dir()):
        for kind in ("fit", "ft"):
            name = "c%d%s" % (k, kind)
            da, db = A / ("c%d" % k) / kind / "final", B / ("c%d" % k) / kind / "final"
            if not (da.exists() and db.exists()):
                report["stages"][name] = {"status": "missing final"}
                allbit = False
                print("%-6s missing final (A %s, B %s)" % (name, da.exists(), db.exists()))
                continue
            rec, bit = compare_stage_generic(da.resolve(), db.resolve())
            report["stages"][name] = rec
            allbit = allbit and bit
            print("%-6s %s" % (name, rec["status"]))
            for n, v in rec.items():
                if n == "status":
                    continue
                if isinstance(v, dict) and "arrays" in v:
                    bad = {a: x for a, x in v["arrays"].items() if x != "bitwise"}
                    badm = {a: x for a, x in v["meta_equal"].items() if not x}
                    print("    %s: %d arrays %s; meta %s" % (n, len(v["arrays"]), "all bitwise" if not bad else bad,
                                                            "equal" if not badm else "DIFFERENT %s" % badm))
                elif isinstance(v, dict) and "lines" in v:
                    print("    %s: %d lines, %d differing%s" % (n, v["lines"], v["n_differing"],
                                                              "" if not v["n_differing"] else " %s" % v["differing_lines"][:2]))
                else:
                    print("    %s: %s" % (n, v))
    report["all_bitwise"] = allbit
    print("ALL STAGES BITWISE IDENTICAL (checked)" if allbit else "NOT all stages bitwise identical")
    if out_json:
        Path(out_json).write_text(json.dumps(report, indent=1, default=str) + "\n")
    return 0 if allbit else 1


def main():
    if len(sys.argv) > 1 and sys.argv[1] == "nonnegkind":        # nonnegkind <code dir> <scratch> (CPU job)
        nonneg_kind(sys.argv[2], sys.argv[3], REF3A / "c1" / "ft" / "seg01" / "finetune_state.npz")
        print("%d FAILED: %s" % (len(FAILS), FAILS) if FAILS else "NONNEG KIND TESTS PASSED")
        return 1 if FAILS else 0
    if len(sys.argv) > 1 and sys.argv[1] == "nonnegstage2":      # nonnegstage2 <code dir> <scratch> (CPU job)
        nonneg_stage2_kind(sys.argv[2], sys.argv[3], REF3A / "c1" / "ft" / "seg01" / "finetune_state.npz")
        print("%d FAILED: %s" % (len(FAILS), FAILS) if FAILS else "NONNEG STAGE-2 KIND TESTS PASSED")
        return 1 if FAILS else 0
    if len(sys.argv) > 1 and sys.argv[1] == "realfit":           # realfit <v3 code dir> <v4 code dir> <scratch>
        res = {}
        for tag, d in (("v_old", sys.argv[2]), ("v_new", sys.argv[3])):
            print("=== %s: %s" % (tag, d), flush=True)
            sc = Path(sys.argv[4]) / tag
            sc.mkdir(parents=True)
            res[tag] = realfit(d, sc)
        old, new = res["v_old"] or {}, res["v_new"] or {}
        check(old.get("models_reaching_registry") == 1, "realfit: the old version stops at the second model "
                                                         "(reproduces job 123990)")
        check(new.get("models_reaching_registry") == 2, "realfit: the new version passes the second model (fit_cops' M)")
        v = new.get("variants", {})
        check(v.get("same as M (another stage)") == "accepted" and v.get("useNG of row 0 flipped") == "accepted",
              "realfit: M again in another stage and a changed selection column are accepted")
        bad = ("another run directory", "T_RG of row 0 changed", "x of row 5 changed", "another Y", "other rows",
               "another TemperatureNG")
        check(all(str(v.get(k, "")).startswith("refused") for k in bad),
              "realfit: another base, geometry column, T_RG, Y, rows or header setting is refused (variants kept alive)")
        check(all(str(v.get("freed: " + k, "")).startswith("refused") for k in bad)
              and v.get("freed: same as M (another stage)") == "accepted",
              "realfit: the same refusals with every variant freed after its check (no id() reuse; v4's defect)")
        print(json.dumps(res, indent=1, default=str))
        print("%d FAILED: %s" % (len(FAILS), FAILS) if FAILS else "REALFIT TESTS PASSED")
        return 1 if FAILS else 0
    if len(sys.argv) > 1 and sys.argv[1] == "compare":
        ap = argparse.ArgumentParser()
        ap.add_argument("cmd")
        ap.add_argument("A")
        ap.add_argument("B")
        ap.add_argument("--json", default=None)
        ap.add_argument("--cycles", type=int, default=None)
        ap.add_argument("--generic", action="store_true", help="every npz / table / result file (slope, nonneg kinds)")
        a = ap.parse_args()
        if a.generic:
            return compare_generic(a.A, a.B, a.json)
        return compare(a.A, a.B, a.cycles, a.json)
    ap = argparse.ArgumentParser()
    ap.add_argument("--scratch", default=None, help="scratch directory (default: a new one under cops_output/tests)")
    a = ap.parse_args()
    base = Path(a.scratch) if a.scratch else ROOT_OUT / "tests"
    base.mkdir(parents=True, exist_ok=True)
    scratch = Path(tempfile.mkdtemp(prefix="test_chain_onejob_cpu_", dir=str(base)))
    print("scratch %s" % scratch)
    for t in (test_commands, test_slope_commands, test_nonneg_commands, test_nonneg_stage2, test_nonneg_delta, test_disk_guard, test_states, test_baseline, test_real_modules, test_spawn,
              test_slope_kind, test_tee):
        try:
            t(scratch)
        except Exception as exc:
            import traceback
            traceback.print_exc()
            check(False, "%s raised %s" % (t.__name__, exc))
    print("%d FAILED: %s" % (len(FAILS), FAILS) if FAILS else "ALL CPU TESTS PASSED")
    return 1 if FAILS else 0


if __name__ == "__main__":
    sys.exit(main())
