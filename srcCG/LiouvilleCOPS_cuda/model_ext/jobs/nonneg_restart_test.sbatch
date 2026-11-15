#!/bin/bash
# ==== version-controlled copy: header added 2026-11-14, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/jobs/nonneg_restart_test.sbatch (on Voyager: /data001/highresolu_vt/cops_output/jobs/nonneg_restart_test.sbatch)
# SHA256 of the run copy = SHA256 of this file without these header lines: 5f28bd521d75b5a6afb4cdc51b8e00447d59a911c38f2e58392d56e4dbb9555a
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: job. GPU test of --start restart (background, logspread with --spread-u0, sos) from the stage-1 end
#   points: start data cost = source final, x_start = source x, the u0 penalty at the start; NFEV 5, CHECK 1.
# RUN: sbatch --export=ALL,MODEL_EXT=<snapshot>,TAG=<tag> nonneg_restart_test.sbatch (test 125365).
# ==== end of header ====
# GPU test of the stage-1 continuations (user 2026-11-14 ~17:30 CST): --start restart for background, logspread
# (with --spread-u0) and sos, each from its stage-1 fit's final state, NFEV 5 with the full Jacobian check (CHECK 1),
# L_THREADS 32, SLOPE_MULT_MIN 1. Checks: the start data cost equals the source fit's final data cost; x_start equals the
# source x; logspread: the penalty at the start = f C_ref [Var_w(u) + (u_bar - u0)^2] and Var_w(u) at the start = the
# source's final Var_w(u); the Jacobian check incl. the new penalty row.
# Settings: MODEL_EXT (required), TAG (required), U0 (default -1.864126788612).
#SBATCH --job-name=cops12m_nonneg_restart_test
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gpus-per-node=8
#SBATCH --cpus-per-task=64
#SBATCH --mem=0
#SBATCH --exclusive
#SBATCH --time=01:00:00
#SBATCH --partition=gpu
#SBATCH --account=general
#SBATCH --exclude=gpu003
#SBATCH --output=/data001/highresolu_vt/cops_output/jobs/logs/cops12m_nonneg_restart_test-%j.out
set -uo pipefail
export AMPS=/data001/highresolu_vt/AMPS APP=/data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda
export COPS_OUTPUT_ROOT=/data001/highresolu_vt/cops_output
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1
MODEL_EXT=$(readlink -f "${MODEL_EXT:?set MODEL_EXT}"); TAG=${TAG:?set TAG}; U0=${U0:--1.864126788612}
R=$COPS_OUTPUT_ROOT/runs/shap8_012m_lategrid_dz1.5; N=$R/nonneg; J=$COPS_OUTPUT_ROOT/jobs; T_=$R/nonneg_tests
D=$T_/restart_test_$TAG; [ -e "$D" ] && { echo "exists: $D" >&2; exit 2; }
mkdir -p "$D"; P=$D/provenance-$SLURM_JOB_ID; mkdir -p "$P"; cd "$APP" || exit 1
cp -p "$0" "$P/"; cp -p "$J/nonneg_fit.sbatch" "$P/"; cp -p "$MODEL_EXT"/*.py "$MODEL_EXT/SHA256SUMS" "$P/" 2>/dev/null
{ git -C "$AMPS" rev-parse HEAD; git -C "$AMPS" status --short; } > "$P/git.txt" 2>&1
lscpu > "$P/lscpu.txt" 2>&1; nvidia-smi -L > "$P/nvidia_smi_L.txt" 2>&1
{ echo '#!/bin/bash'; echo "# job $SLURM_JOB_ID on $(hostname), $(date -u +%FT%TZ)"; grep '^export ' "$0"; echo "cd $APP"; } > "$P/commands.sh"
echo "$(date -u +%FT%TZ) job $SLURM_JOB_ID on $(hostname): MODEL_EXT=$MODEL_EXT TAG=$TAG U0=$U0"
rc=0
run() {
    local label=$1; shift
    echo "$(date -u +%FT%TZ) $label"
    echo "env $* # $label" >> "$P/commands.sh"
    env "$@" > "$P/$label.log" 2>&1; local r=$?
    echo "  exit $r"
    grep -E "nonneg: (start|alternative)|start check \(start|^lambda|PASS|FAIL|Traceback|Error" "$P/$label.log" | cut -c1-420
    return $r
}
BG=$N/background/stage1_bg1e-2_20261114T084711Z/fit/final/fit_restart.npz
LS=$N/logspread/stage1_f0.1_20261114T090251Z/fit/final/fit_restart.npz
SO=$N/sos/stage1_r3_20261114T091354Z/fit/final/fit_restart.npz
C=(SLOPE_MULT_MIN=1 CHECK=1 NFEV=5 L_THREADS=32 MODEL_EXT=$MODEL_EXT START=restart)
run bg TESTNAME=g_restart_bg_$TAG MODEL=background BG_MIN=1e-2 START_RESTART=$BG "${C[@]}" bash "$J/nonneg_fit.sbatch" || rc=1
run ls TESTNAME=g_restart_ls_$TAG MODEL=logspread SPREAD_F=0.1 SPREAD_CREF=4.2688066575e-03 SPREAD_U0=$U0 START_RESTART=$LS \
    "${C[@]}" bash "$J/nonneg_fit.sbatch" || rc=1
run sos TESTNAME=g_restart_sos_$TAG MODEL=sos SOS_RANK=3 START_RESTART=$SO "${C[@]}" bash "$J/nonneg_fit.sbatch" || rc=1
"$APP/.venv/bin/python" - "$T_" "$TAG" "$BG" "$LS" "$SO" > "$P/summary.log" 2>&1 <<'PY' || rc=1
import json, sys
import numpy as np
T, tag, srcs = sys.argv[1], sys.argv[2], dict(zip(("bg", "ls", "sos"), sys.argv[3:6]))
bad = 0
for lab in ("bg", "ls", "sos"):
    seg = "%s/g_restart_%s_%s/seg01" % (T, lab, tag)
    m = json.load(open(seg + "/fit_manifest.json")); rec = m["model_extension"]["nonneg"]
    src = srcs[lab]
    msrc = json.load(open(src.replace("fit_restart.npz", "fit_manifest.json")))
    c0 = rec["clip_core"]["start_check"]["start"]["data_cost"]; cs = float(msrc["optimizer_results"][0]["data_cost"])
    xs = np.asarray(rec["start_restart"]["x_start"]); xsrc = np.load(src)["x"]
    jc = rec["clip_core"]["jacobian_check"]["columns"]
    worst = max(v["data"] for v in jc.values() if isinstance(v, dict))
    ok = abs(c0 / cs - 1) <= 1e-12 and np.array_equal(xs, xsrc)
    line = "%s: start data cost %.15e, source final %.15e (rel %.1e); x_start == source x %s; Jacobian check worst %.2e over %d columns; cost after 5 %s" % (
        lab, c0, cs, abs(c0 / cs - 1), bool(np.array_equal(xs, xsrc)), worst, len(jc), m["optimizer_results"][0]["cost"])
    if lab == "ls":
        sp = rec["model_check"]["alternative"]["spread"]
        src_sp = msrc["model_extension"]["nonneg"]["model_check"]["alternative"]["spread"]["final"]
        want = sp["f"] * sp["C_ref"] * (sp["var_w_u_at_start"] + sp["dev2_mean_u0_at_start"])
        ok2 = (abs(sp["penalty_cost_at_start"] / want - 1) <= 1e-12 and
               abs(sp["var_w_u_at_start"] / src_sp["var_w_u"] - 1) <= 1e-12 and sp["u0"] is not None)
        ok = ok and ok2
        line += "; u0 %r: start penalty %.12e = f C_ref (Var %.10e + (u_bar - u0)^2 %.10e) (rel %.1e); Var at the start vs the source's final %.10e (rel %.1e); final %s" % (
            sp["u0"], sp["penalty_cost_at_start"], sp["var_w_u_at_start"], sp["dev2_mean_u0_at_start"],
            abs(sp["penalty_cost_at_start"] / want - 1), src_sp["var_w_u"], abs(sp["var_w_u_at_start"] / src_sp["var_w_u"] - 1),
            {k: sp["final"].get(k) for k in ("data_cost", "penalty", "var_w_u", "mean_u", "dev2_mean_u0")})
    print(("PASS " if ok else "FAIL ") + line)
    bad += not ok
sys.exit(1 if bad else 0)
PY
cat "$P/summary.log"
echo "$(date -u +%FT%TZ) exit $rc"
exit $rc
