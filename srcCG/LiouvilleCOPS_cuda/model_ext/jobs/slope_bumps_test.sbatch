#!/bin/bash
# ==== version-controlled copy: header added 2026-11-11, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/jobs/slope_bumps_test.sbatch (on Voyager: /data001/highresolu_vt/cops_output/jobs/slope_bumps_test.sbatch)
# SHA256 of the run copy = SHA256 of this file without these header lines: dc231591dc71ceda0d329ef8e07fb0ad74d1b106030707421cb507c090423e2a
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: job (test). GPU tests of the bump + slope-multiplier model (option 3a, 2026-11-10): STEP=fits (floor-0 check
#   against an unpatched independent evaluation; NFEV-5 short fits at L_THREADS 32 / 48 / 64, bitwise) and STEP=ft (fine
#   tuning 32 vs 1 thread bitwise, checks 1-8, cycle-2 start cost).
# RUN: sbatch --nodelist=gpu002 --export=ALL,STEP=ft,MODEL_EXT=<snapshot>,TAG=<tag> slope_bumps_test.sbatch   (123984
#   fits; ft: 123985 v6 and 124013 v7 failed, 124063 v8 PASS)
# ==== end of header ====
# GPU tests of the bump + slope-multiplier model and its fine-tuning cycle (option 3a, user decisions 2026-11-10), one node:
#   STEP=fits: (1) model check: floor ratio 0, all multipliers 1, the NNLS bump weights, stop after the start check, then
#              an UNPATCHED process evaluates the plain clipped model with the 32 bump columns as its only templates
#              (a = 0) -> compared with the start check (test_slope_mult.py --gpu-independent);
#              (2) the short fit (START=nnls from D2's grain, NFEV 5, CHECK 1: Jacobian check of every weight and
#              multiplier, kappa, grain) with L_THREADS = each of THREADS (default "32 48 64"); (3) bitwise comparison of
#              the later thread counts with the first (test_slope_mult.py --compare-runs).
#   STEP=ft:   (4) fine tuning (late_only, WINDOW_MIN 5, VALIDATE 1, MAXNFEV FT_NFEV) of the first short fit with
#              L_THREADS=FT_THREADS and with 1 thread -> bitwise comparison of the products; (5) the cycle-2 fit (START
#              =chain, HELD, GRAIN_FROM = that fine tuning, NFEV 3, CHECK 0): its start cost = the fine tuning's final
#              data cost (1e-8).
# Settings: MODEL_EXT (required, frozen snapshot), TAG (required), STEP (fits | ft), THREADS, FT_THREADS (default 32),
#   FT_NFEV (default 30). Products under runs/shap8_012m_lategrid_dz1.5/slope_mult_tests/<name>_<TAG>.
#SBATCH --job-name=cops12m_bumps_test
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
#SBATCH --output=/data001/highresolu_vt/cops_output/jobs/logs/cops12m_bumps_test-%j.out
set -uo pipefail
export AMPS=/data001/highresolu_vt/AMPS APP=/data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda
export COPS_OUTPUT_ROOT=/data001/highresolu_vt/cops_output
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1
MODEL_EXT=$(readlink -f "${MODEL_EXT:?set MODEL_EXT}"); TAG=${TAG:?set TAG}; STEP=${STEP:?set STEP=fits or ft}
THREADS=${THREADS:-32 48 64}; FT_THREADS=${FT_THREADS:-32}; FT_NFEV=${FT_NFEV:-30}
R=$COPS_OUTPUT_ROOT/runs/shap8_012m_lategrid_dz1.5; J=$COPS_OUTPUT_ROOT/jobs; T_=$R/slope_mult_tests
BUMPS=$R/slope_bumps/bundle/bumps32.npz
PY="$APP/.venv/bin/python -W ignore"
D=$T_/bumps_test_${STEP}_$TAG; [ -e "$D" ] && { echo "exists: $D" >&2; exit 2; }
mkdir -p "$D"; P=$D/provenance-$SLURM_JOB_ID; mkdir -p "$P"; cd "$APP" || exit 1
cp -p "$0" "$P/slope_bumps_test.sbatch"; cp -p "$J/slope_mult_fit.sbatch" "$J/slope_finetune.sbatch" "$P/"
cp -p "$MODEL_EXT"/*.py "$MODEL_EXT/SHA256SUMS" "$P/" 2>/dev/null
{ git -C "$AMPS" rev-parse HEAD; git -C "$AMPS" status --short; } > "$P/git.txt" 2>&1
nvidia-smi -L > "$P/nvidia_smi_L.txt" 2>&1
{ echo '#!/bin/bash'; echo "# job $SLURM_JOB_ID on $(hostname), $(date -u +%FT%TZ)"; grep '^export ' "$0"; echo "cd $APP"; } > "$P/commands.sh"
echo "$(date -u +%FT%TZ) job $SLURM_JOB_ID on $(hostname): STEP=$STEP MODEL_EXT=$MODEL_EXT TAG=$TAG"
rc=0
run() {   # label, env assignments..., script
    local label=$1; shift
    echo "$(date -u +%FT%TZ) $label"
    echo "env $* # $label" >> "$P/commands.sh"
    env "$@" > "$P/$label.log" 2>&1; local r=$?
    echo "  exit $r"; grep -E "slope_mult: (start|model check|L = Y a)|start check \(start|clip_refit: Jacobian check|^lambda|gpu_fit: [0-9]|PASS|FAIL|slope_finetune: (check|start|final|late_only)|peak used|complete" "$P/$label.log" | cut -c1-420
    return $r
}
if [ "$STEP" = fits ]; then
    run floor0 TESTNAME=g_bumps_floor0_$TAG TESTMODE=floor0 BUMPS=$BUMPS START=nnls CHECK=0 NFEV=5 MODEL_EXT=$MODEL_EXT \
        bash "$J/slope_mult_fit.sbatch" || rc=1
    SEGS=()
    for N in $THREADS; do
        run "short_l$N" TESTNAME=g_bumps_short_l${N}_$TAG BUMPS=$BUMPS START=nnls CHECK=1 NFEV=5 L_THREADS=$N \
            MODEL_EXT=$MODEL_EXT bash "$J/slope_mult_fit.sbatch" || rc=1
        SEGS+=("$T_/g_bumps_short_l${N}_$TAG/seg01")
    done
    CCMD=($PY "$MODEL_EXT/test_slope_mult.py" --compare-runs "${SEGS[@]}" --json-out "$D/compare.json")
    printf '%q ' "${CCMD[@]}" >> "$P/commands.sh"; echo >> "$P/commands.sh"
    "${CCMD[@]}" > "$P/compare.log" 2>&1 || rc=1
    grep -E "PASS|FAIL|SUMMARY" "$P/compare.log" | cut -c1-900
else
    F=$T_/g_bumps_short_l32_$TAG/final
    [ -e "$F" ] || F=$(ls -d "$T_"/g_bumps_short_l32_*/final 2>/dev/null | tail -1)
    [ -e "$F" ] || { echo "no 32-thread short fit found" >&2; exit 2; }
    F=$(readlink -f "$F"); echo "fit: $F"
    run ft_l$FT_THREADS FIT=$F OUT=$T_/g_bumps_ft_l${FT_THREADS}_$TAG FACETS=late_only WINDOW_MIN=5 MAXNFEV=$FT_NFEV \
        L_THREADS=$FT_THREADS VALIDATE=1 MODEL_EXT=$MODEL_EXT bash "$J/slope_finetune.sbatch" || rc=1
    run ft_l1 FIT=$F OUT=$T_/g_bumps_ft_l1_$TAG FACETS=late_only WINDOW_MIN=5 MAXNFEV=$FT_NFEV VALIDATE=1 \
        MODEL_EXT=$MODEL_EXT bash "$J/slope_finetune.sbatch" || rc=1
    $PY - "$T_/g_bumps_ft_l${FT_THREADS}_$TAG" "$T_/g_bumps_ft_l1_$TAG" > "$P/ft_compare.log" 2>&1 <<'PY' || rc=1
import json, sys
import numpy as np
a, b = sys.argv[1], sys.argv[2]
za, zb = np.load(a + "/finetune_state.npz"), np.load(b + "/finetune_state.npz")
keys = ["delta", "held_next", "A_final", "grain_A_uPa", "grain_B_uPa_per_100s", "lower_bound"]
same = {k: bool(np.array_equal(za[k], zb[k], equal_nan=True)) for k in keys}
ta = [l for l in open(a + "/COPS_model.dat") if not l.startswith("#")]
tb = [l for l in open(b + "/COPS_model.dat") if not l.startswith("#")]
ra = [l.split()[:4] for l in open(a + "/finetune_result.txt") if l.split() and l.split()[0].startswith("cost")]
rb = [l.split()[:4] for l in open(b + "/finetune_result.txt") if l.split() and l.split()[0].startswith("cost")]
fa, fb = json.load(open(a + "/fine_tuning.json")), json.load(open(b + "/fine_tuning.json"))
same["COPS_model.dat"] = ta == tb
same["finetune_result costs"] = ra == rb
same["nfev"] = fa["optimizer"]["nfev"] == fb["optimizer"]["nfev"]
ok = all(same.values())
print("fine tuning %s threads vs 1 thread: bitwise %s %s" % (a.split("_l")[-1], "PASS" if ok else "FAIL", same))
for f, nm in ((fa, "threads"), (fb, "1 thread")):
    print("  %s: timings %s; host seconds %s; gpu seconds %s; calls %s; wall %.1f s" % (
        nm, {k: round(v, 2) for k, v in f["timings"].items()}, {k: round(v, 3) for k, v in f["engine_seconds"]["host"].items()},
        {k: round(v, 2) for k, v in (f["engine_seconds"]["gpu"] or {}).items()}, f["engine_seconds"]["gpu_calls"], f["wall_s"]))
sys.exit(0 if ok else 1)
PY
    cat "$P/ft_compare.log" | cut -c1-900
    FTS=$T_/g_bumps_ft_l${FT_THREADS}_$TAG/finetune_state.npz
    run c2 TESTNAME=g_bumps_c2_$TAG BUMPS=$BUMPS START=chain START_RESTART=$F/fit_restart.npz HELD=$FTS GRAIN_FROM=$FTS \
        CHECK=0 NFEV=3 L_THREADS=$FT_THREADS MODEL_EXT=$MODEL_EXT bash "$J/slope_mult_fit.sbatch" || rc=1
    $PY - "$FTS" "$T_/g_bumps_c2_$TAG/seg01/fit_manifest.json" > "$P/c2_check.log" 2>&1 <<'PY' || rc=1
import json, sys
import numpy as np
with np.load(sys.argv[1]) as z:
    meta = json.loads(str(z["meta_json"]))
ft = meta["costs_final"]["cost_data"]
m = json.load(open(sys.argv[2]))
c2 = m["model_extension"]["slope_multiplier"]["clip_core"]["start_check"]["start"]["data_cost"]
rel = abs(c2 / ft - 1.0)
print("cycle-2 fit start data cost %.15e vs the fine tuning's final data cost %.15e: rel %.2e %s" % (
    c2, ft, rel, "PASS" if rel <= 1e-8 else "FAIL"))
sys.exit(0 if rel <= 1e-8 else 1)
PY
    cat "$P/c2_check.log"
fi
echo "$(date -u +%FT%TZ) STEP=$STEP exit $rc"
exit $rc
