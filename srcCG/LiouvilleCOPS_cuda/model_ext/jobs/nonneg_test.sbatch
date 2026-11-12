#!/bin/bash
# ==== version-controlled copy: header added 2026-11-11, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/jobs/nonneg_test.sbatch (on Voyager: /data001/highresolu_vt/cops_output/jobs/nonneg_test.sbatch)
# SHA256 of the run copy = SHA256 of this file without these header lines: a1d7ff5fcbb45e14256fa0a7fe634f44a75d5dc11f1c5ff9e553ffca04f447ba
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: job (test). GPU tests of the softplus / squared models (2026-11-11), one MODEL per job: STEP=fits (independent
#   start check against an unpatched process; softplus s = 0 clip check; NFEV-5 fit with the full Jacobian check) and
#   STEP=ft (fine tuning 32 vs 1 thread bitwise, checks 1-8, cycle-2 start cost).
# RUN: sbatch --nodelist=gpu001 --export=ALL,MODEL_EXT=<snapshot>,TAG=<tag>,MODEL=softplus,STEP=fits nonneg_test.sbatch
#   (124130-124133, all PASS)
# ==== end of header ====
# GPU tests of the non-negative activity models (user decisions 2026-11-11), one node, one MODEL per job:
#   STEP=fits: (1) TESTMODE=independent: the driver's start check vs an UNPATCHED process evaluating the independently
#              computed activity (test_nonneg.py --gpu-independent; squared bitwise, softplus rel <= 1e-12);
#              (1b) softplus only, TESTMODE=clip: s = 0, m = 1 == the plain clipped model (bitwise);
#              (2) the short fit (NFEV 5, CHECK 1: full Jacobian check of every parameter) with L_THREADS (32).
#   STEP=ft:   (3) fine tuning (late_only, WINDOW_MIN 5, VALIDATE 1, MAXNFEV FT_NFEV) of that short fit with L_THREADS and
#              with 1 thread -> bitwise comparison; (4) the cycle-2 fit (START=chain, HELD, GRAIN_FROM = that fine tuning,
#              NFEV 3): its start cost = the fine tuning's final data cost (1e-8).
# Settings: MODEL_EXT (frozen snapshot, required), TAG (required), MODEL (softplus|squared), STEP (fits|ft), L_THREADS
#   (default 32), FT_NFEV (default 30), START_FILE (squared; default runs/.../nonneg_tests/start/squared_start.npz).
# Products under runs/shap8_012m_lategrid_dz1.5/nonneg_tests/<name>_<TAG>.
#SBATCH --job-name=cops12m_nonneg_test
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
#SBATCH --output=/data001/highresolu_vt/cops_output/jobs/logs/cops12m_nonneg_test-%j.out
set -uo pipefail
export AMPS=/data001/highresolu_vt/AMPS APP=/data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda
export COPS_OUTPUT_ROOT=/data001/highresolu_vt/cops_output
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1
MODEL_EXT=$(readlink -f "${MODEL_EXT:?set MODEL_EXT}"); TAG=${TAG:?set TAG}; STEP=${STEP:?set STEP=fits or ft}
MODEL=${MODEL:?set MODEL=softplus or squared}; LT=${L_THREADS:-32}; FT_NFEV=${FT_NFEV:-30}
R=$COPS_OUTPUT_ROOT/runs/shap8_012m_lategrid_dz1.5; J=$COPS_OUTPUT_ROOT/jobs; T_=$R/nonneg_tests
START_FILE=$(readlink -f "${START_FILE:-$T_/start/squared_start.npz}")
PY="$APP/.venv/bin/python -W ignore"
D=$T_/test_${MODEL}_${STEP}_$TAG; [ -e "$D" ] && { echo "exists: $D" >&2; exit 2; }
mkdir -p "$D"; P=$D/provenance-$SLURM_JOB_ID; mkdir -p "$P"; cd "$APP" || exit 1
cp -p "$0" "$P/nonneg_test.sbatch"; cp -p "$J/nonneg_fit.sbatch" "$J/nonneg_finetune.sbatch" "$P/"
cp -p "$MODEL_EXT"/*.py "$MODEL_EXT/SHA256SUMS" "$P/" 2>/dev/null
{ git -C "$AMPS" rev-parse HEAD; git -C "$AMPS" status --short; } > "$P/git.txt" 2>&1
nvidia-smi -L > "$P/nvidia_smi_L.txt" 2>&1; df -h /dev/shm > "$P/dev_shm.txt" 2>&1
{ echo '#!/bin/bash'; echo "# job $SLURM_JOB_ID on $(hostname), $(date -u +%FT%TZ)"; grep '^export ' "$0"; echo "cd $APP"; } > "$P/commands.sh"
echo "$(date -u +%FT%TZ) job $SLURM_JOB_ID on $(hostname): MODEL=$MODEL STEP=$STEP MODEL_EXT=$MODEL_EXT TAG=$TAG"
SX=(); [ "$MODEL" = squared ] && SX=(START_FILE=$START_FILE)
rc=0
run() {   # label, env assignments..., script
    local label=$1; shift
    echo "$(date -u +%FT%TZ) $label"
    echo "env $* # $label" >> "$P/commands.sh"
    env "$@" > "$P/$label.log" 2>&1; local r=$?
    echo "  exit $r"
    grep -E "nonneg: (start|model check|held term)|start check \(start|clip_refit: Jacobian check d/d|^lambda|gpu_fit: [0-9]|PASS|FAIL|nonneg_finetune: (check|start|final|late_only)|peak used|complete|Traceback|Error" "$P/$label.log" | cut -c1-420
    return $r
}
if [ "$STEP" = fits ]; then
    run indep TESTNAME=g_${MODEL}_indep_$TAG TESTMODE=independent MODEL=$MODEL "${SX[@]}" CHECK=0 NFEV=5 \
        L_THREADS=$LT MODEL_EXT=$MODEL_EXT bash "$J/nonneg_fit.sbatch" || rc=1
    if [ "$MODEL" = softplus ]; then
        run clip TESTNAME=g_${MODEL}_clip_$TAG TESTMODE=clip MODEL=$MODEL CHECK=0 NFEV=5 L_THREADS=$LT \
            MODEL_EXT=$MODEL_EXT bash "$J/nonneg_fit.sbatch" || rc=1
    fi
    run short TESTNAME=g_${MODEL}_short_$TAG MODEL=$MODEL "${SX[@]}" CHECK=1 NFEV=5 L_THREADS=$LT \
        CHECKPARAMS=$(echo a{00..24} | tr ' ' :):c_slope_40_50:c_slope_50_60:c_slope_60_70:c_slope_70_80:c_slope_80_90:c_slope_90_180:kappa:grain_A_uPa:grain_B_uPa_per_100s \
        MODEL_EXT=$MODEL_EXT bash "$J/nonneg_fit.sbatch" || rc=1
else
    F=$(ls -d "$T_"/g_${MODEL}_short_*/final 2>/dev/null | tail -1)
    [ -n "$F" ] && [ -e "$F" ] || { echo "no short fit of $MODEL found" >&2; exit 2; }
    F=$(readlink -f "$F"); echo "fit: $F"
    run ft_l$LT FIT=$F OUT=$T_/g_${MODEL}_ft_l${LT}_$TAG FACETS=late_only WINDOW_MIN=5 MAXNFEV=$FT_NFEV L_THREADS=$LT \
        VALIDATE=1 MODEL_EXT=$MODEL_EXT bash "$J/nonneg_finetune.sbatch" || rc=1
    run ft_l1 FIT=$F OUT=$T_/g_${MODEL}_ft_l1_$TAG FACETS=late_only WINDOW_MIN=5 MAXNFEV=$FT_NFEV VALIDATE=1 \
        MODEL_EXT=$MODEL_EXT bash "$J/nonneg_finetune.sbatch" || rc=1
    $PY - "$T_/g_${MODEL}_ft_l${LT}_$TAG" "$T_/g_${MODEL}_ft_l1_$TAG" > "$P/ft_compare.log" 2>&1 <<'PY' || rc=1
import json, sys
import numpy as np
a, b = sys.argv[1], sys.argv[2]
za, zb = np.load(a + "/finetune_state.npz"), np.load(b + "/finetune_state.npz")
keys = ["delta", "held_next", "A_final", "grain_A_uPa", "grain_B_uPa_per_100s", "lower_bound"]
same = {k: bool(np.array_equal(za[k], zb[k], equal_nan=True)) for k in keys}
ta = [l for l in open(a + "/COPS_model.dat") if not l.startswith("#")]
tb = [l for l in open(b + "/COPS_model.dat") if not l.startswith("#")]
fa, fb = json.load(open(a + "/fine_tuning.json")), json.load(open(b + "/fine_tuning.json"))
same["COPS_model.dat"] = ta == tb
same["nfev"] = fa["optimizer"]["nfev"] == fb["optimizer"]["nfev"]
ok = all(same.values())
print("fine tuning threads vs 1 thread: bitwise %s %s" % ("PASS" if ok else "FAIL", same))
for f, nm in ((fa, "threads"), (fb, "1 thread")):
    print("  %s: timings %s; host seconds %s; gpu seconds %s; calls %s; wall %.1f s; validation %s" % (
        nm, {k: round(v, 2) for k, v in f["timings"].items()}, {k: round(v, 3) for k, v in f["engine_seconds"]["host"].items()},
        {k: round(v, 2) for k, v in (f["engine_seconds"]["gpu"] or {}).items()}, f["engine_seconds"]["gpu_calls"],
        f["wall_s"], f.get("validation")))
sys.exit(0 if ok else 1)
PY
    cat "$P/ft_compare.log" | cut -c1-1200
    FTS=$T_/g_${MODEL}_ft_l${LT}_$TAG/finetune_state.npz
    run c2 TESTNAME=g_${MODEL}_c2_$TAG MODEL=$MODEL START=chain START_RESTART=$F/fit_restart.npz HELD=$FTS GRAIN_FROM=$FTS \
        CHECK=0 NFEV=3 L_THREADS=$LT MODEL_EXT=$MODEL_EXT bash "$J/nonneg_fit.sbatch" || rc=1
    $PY - "$FTS" "$T_/g_${MODEL}_c2_$TAG/seg01/fit_manifest.json" > "$P/c2_check.log" 2>&1 <<'PY' || rc=1
import json, sys
import numpy as np
with np.load(sys.argv[1]) as z:
    meta = json.loads(str(z["meta_json"]))
ft = meta["costs_final"]["cost_data"]
m = json.load(open(sys.argv[2]))
c2 = m["model_extension"]["nonneg"]["clip_core"]["start_check"]["start"]["data_cost"]
rel = abs(c2 / ft - 1.0)
print("cycle-2 fit start data cost %.15e vs the fine tuning's final data cost %.15e: rel %.2e %s" % (
    c2, ft, rel, "PASS" if rel <= 1e-8 else "FAIL"))
sys.exit(0 if rel <= 1e-8 else 1)
PY
    cat "$P/c2_check.log"
fi
echo "$(date -u +%FT%TZ) STEP=$STEP exit $rc"
exit $rc
