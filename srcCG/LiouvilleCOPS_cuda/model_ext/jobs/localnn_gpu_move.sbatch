#!/bin/bash
# ==== version-controlled copy: header added 2026-11-15, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/jobs/localnn_gpu_move.sbatch (on Voyager: /data001/highresolu_vt/cops_output/jobs/localnn_gpu_move.sbatch)
# SHA256 of the run copy = SHA256 of this file without these header lines: 5f8feea323b320e6f44a52a28a95d99cebed2bff3cb709ef53f71d17d192fcb7
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: after the GPU solver verification (afterok): if fixed state, one-step and trajectory tests pass and GPU
#   s/eval <= 0.5 x CPU, run localnn_chain_to_gpu.sh for the four stage-2 chains with the fastest passing SVD mode;
#   otherwise nothing (exit 3).
# RUN 2026-11-15: 132811 (afterok:132810).
# ==== end of header ====
# Option C stage 2 (lead 2026-11-15 "option 3": "once tested and a GPU is free, move the remaining cycles there"):
# after the GPU solver verification job VERIFY_JOB (afterok), read its verify.json and, if
#   (1) the fixed-state residuals and Jacobian passed (1e-13 of max |r|; 1e-12 per Jacobian column) and
#   (2) for an SVD mode the one-step test (1e-8) and the 50-evaluation trajectory (1e-10) passed and
#   (3) that mode's GPU seconds per evaluation <= MAX_RATIO x the CPU reference's (default 0.5),
# run jobs/localnn_chain_to_gpu.sh --submit for every chain in CHAINS with the fastest passing mode: GPU attempts that take
# their chain over when they start (they cancel the chain's CPU attempts then). Otherwise nothing is moved (exit 3).
# Settings: VERIFY_JOB (required), CHAINS (required, space separated), MODEL_EXT (the frozen GPU snapshot), MAX_RATIO.
#SBATCH --job-name=lnn_gpu_move
#SBATCH --partition=normal128cores
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:15:00
#SBATCH --output=/data001/highresolu_vt/cops_output/jobs/logs/lnn_gpu_move-%j.out
set -uo pipefail
ROOT=/data001/highresolu_vt/cops_output
PY=/data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda/.venv/bin/python
V=$ROOT/tests/localnn_20261114/gpu_solver/verify_${VERIFY_JOB:?VERIFY_JOB}/verify.json
MEXT=$(readlink -f "${MODEL_EXT:?MODEL_EXT}")
echo "$(date -u +%FT%TZ) job $SLURM_JOB_ID on $(hostname): VERIFY=$V MODEL_EXT=$MEXT CHAINS=${CHAINS:?CHAINS}"
MODE=$("$PY" - "$V" "${MAX_RATIO:-0.5}" <<'EOF'
import json, sys
v = json.load(open(sys.argv[1])); ratio = float(sys.argv[2])
if not v.get("1_passed"):
    print("NONE fixed-state residuals / Jacobian failed"); sys.exit()
cpu = v["3_trajectory"]["cpu_seconds_per_evaluation"]
ok = []
for m, r in v["3_trajectory"]["modes"].items():
    st = v["2_trf_step"]["modes"].get(m, {})
    if r.get("passed") and st.get("passed") and r["gpu_seconds_per_evaluation"] <= ratio * cpu:
        ok.append((r["gpu_seconds_per_evaluation"], m))
print(min(ok)[1] if ok else "NONE no mode passed the step / trajectory tests with GPU s/eval <= %g x CPU %.2f s" % (ratio, cpu))
EOF
)
echo "$(date -u +%FT%TZ) decision: $MODE"
case $MODE in NONE*) echo "nothing moved"; exit 3;; esac
rc=0
for C in $CHAINS; do
  bash "$ROOT/jobs/localnn_chain_to_gpu.sh" --chain "$C" --model-ext "$MEXT" --gpu-svd-mode "$MODE" --submit || rc=$?
done
mkdir -p "$ROOT/model_ext/logs/gpu_move-$SLURM_JOB_ID" && cp -p "$0" "$ROOT/model_ext/logs/gpu_move-$SLURM_JOB_ID/localnn_gpu_move.sbatch.used"
echo "$(date -u +%FT%TZ) done rc $rc"
exit $rc
