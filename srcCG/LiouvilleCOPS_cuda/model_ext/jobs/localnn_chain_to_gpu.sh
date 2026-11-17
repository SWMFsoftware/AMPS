#!/bin/bash
# ==== version-controlled copy: header added 2026-11-15, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/jobs/localnn_chain_to_gpu.sh (on Voyager: /data001/highresolu_vt/cops_output/jobs/localnn_chain_to_gpu.sh)
# SHA256 of the run copy = SHA256 of this file without these header lines: a124d0ee936b9d2f856bf54b1afb6cc4a648650631459ac4676da5ae8b4b2cfb
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: move the remaining cycles of one chain to a GPU: GPU attempts (1 GPU, gpu partition, never gpu003) with
#   CANCEL_JOBS = the chain's CPU attempts; the settings of the chain's last CPU submission including TEMPLATES /
#   FT_EXTRA / NONNEG_CODE (2026-11-17); --script, --no-post, --after-checkpoint (TAKEOVER_AFTER_CHECKPOINT=1).
# RUN 2026-11-15: 2026-11-17 17:45Z: stage2r_f10_pm1e-4_tplon (136092-136093), stage2r_f1_pm1e-4_tplon
#   (136095-136096), gpu002.
# ==== end of header ====
# Option C stage 2 (lead 2026-11-15 "option 3"): move the remaining cycles of ONE localnn chain to a GPU.
#   localnn_chain_to_gpu.sh --chain CHAIN --model-ext FROZEN_GPU_SNAPSHOT [--attempts 6] [--time 08:00:00]
#       [--exclude gpu003] [--nodelist NODE] [--gpu-svd-mode qr] [--cpus 8] [--mem 40G] [--submit]
#       [--script FILE (the driver; default MODEL_EXT/localnn_chain.sbatch.copy)] [--no-post]
#       [--after-checkpoint] (2026-11-17: the GPU attempt cancels the CPU attempt only right after its next checkpoint)
# Submits GPU attempts of the chain (localnn_chain_submit.sh --chain CHAIN --device gpu: partition gpu, account general,
# 1 GPU, gpu001 or gpu002, never gpu003 (lead 2026-11-15 13:30 CST); the settings of CHAIN/submission.json) with CANCEL_JOBS = the chain's CPU attempts + its post job (all
# submitted by this session). Nothing is cancelled now: the first GPU attempt, when it STARTS (a GPU is free), writes
# CHAIN/HANDOVER = gpu, cancels those jobs, waits until they have left the queue and continues the chain (a fit in
# progress resumes from its last checkpoint: <= CHECKPOINT evaluations of CPU work lost; an attempt stopped between a
# fit and its direct check runs the check). So the chain never waits for a GPU while the CPU attempts can run.
# A new post job follows the GPU attempts. Without --submit: dry run.
set -euo pipefail
ROOT=/data001/highresolu_vt/cops_output
CHAIN= MEXT= ATT=6 TIME=08:00:00 EXCL=gpu003 NODE= SVDM=qr CPUS=8 MEM=40G SUBMIT=0 SCR= NOPOST= AFTERCK=
while [ $# -gt 0 ]; do
  case $1 in
    --chain) CHAIN=$2; shift;; --model-ext) MEXT=$2; shift;; --attempts) ATT=$2; shift;; --time) TIME=$2; shift;;
    --exclude) EXCL=$2; shift;; --nodelist) NODE=$2; shift;; --gpu-svd-mode) SVDM=$2; shift;;
    --cpus) CPUS=$2; shift;; --mem) MEM=$2; shift;; --submit) SUBMIT=1;; --script) SCR=$(readlink -f "$2"); shift;;
    --no-post) NOPOST=--no-post;; --after-checkpoint) AFTERCK=",TAKEOVER_AFTER_CHECKPOINT=1";; *) echo "unknown option $1"; exit 2;;
  esac; shift
done
CHAIN=$(readlink -f "${CHAIN:?--chain}"); MEXT=$(readlink -f "${MEXT:?--model-ext}")
[ -f "$MEXT/localnn_gpu_solver.py" ] || { echo "$MEXT has no localnn_gpu_solver.py"; exit 2; }
SCRIPT=${SCR:-$MEXT/localnn_chain.sbatch.copy}         # the driver frozen with the code (or --script)
[ -f "$SCRIPT" ] || { echo "$SCRIPT missing"; exit 2; }
PY=/data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda/.venv/bin/python
eval "$("$PY" - "$CHAIN/submission.json" <<'PYEOF'
import json, sys, shlex
j = json.load(open(sys.argv[1])); s = j["settings"]
name = j["sbatch_lines"][0].split("--job-name=lnn_")[1].split()[0].rsplit("_", 1)[0]
for k in ("LAMBDA", "RULE", "FM", "CYCLES", "NFEV", "FT_MAXNFEV", "WINDOW_MIN", "CHECKPOINT", "LATE", "CHECK", "DEVICE",
          "TEMPLATES", "FT_EXTRA", "NONNEG_CODE"):
    print("S_%s=%s" % (k, shlex.quote(str(s.get(k, "")))))
print("S_C1_FROM=%s S_C1_FIT=%s NAME=%s" % (shlex.quote(s.get("C1_FROM", "")), shlex.quote(s.get("C1_FIT", "")), shlex.quote(name)))
print("JOBS=%s" % shlex.quote(":".join(str(x) for x in j["jobs"])))
PYEOF
)"
[ "$S_DEVICE" = cpu ] || { echo "the last submission of $CHAIN is not a CPU one ($S_DEVICE)"; exit 2; }
grep -q '"stage": "chain_complete", "rc": 0' "$CHAIN/chain_events.jsonl" 2>/dev/null && { echo "$CHAIN is complete"; exit 0; }
C1=(); [ -n "$S_C1_FROM" ] && C1=(--c1-from "$S_C1_FROM"); [ -n "$S_C1_FIT" ] && C1=(--c1-fit "$S_C1_FIT")
X="--account=general --gpus-per-node=1 --cpus-per-task=$CPUS --mem=$MEM"; [ -n "$NODE" ] && X="$X --nodelist=$NODE"
CMD=(bash "$ROOT/jobs/localnn_chain_submit.sh" --chain "$CHAIN" --name "${NAME}_gpu" --lambda "$S_LAMBDA" --rule "$S_RULE"
     --fm "$S_FM" "${C1[@]}" --model-ext "$MEXT" --late "$S_LATE" --device gpu --gpu-svd-mode "$SVDM" --cycles "$S_CYCLES"
     --nfev "$S_NFEV" --ft-maxnfev "$S_FT_MAXNFEV" --window-min "$S_WINDOW_MIN" --checkpoint "$S_CHECKPOINT"
     --check "$S_CHECK" --blas-threads "$CPUS" --attempts "$ATT" --time "$TIME" --partition gpu --exclude "$EXCL"
     --sbatch-extra "$X" --env "CANCEL_JOBS=$JOBS$AFTERCK" --script "$SCRIPT" --templates "${S_TEMPLATES:-on}")
[ -n "$S_FT_EXTRA" ] && CMD+=(--ft-extra "$S_FT_EXTRA"); [ -n "$S_NONNEG_CODE" ] && CMD+=(--nonneg-code "$S_NONNEG_CODE")
[ -n "$NOPOST" ] && CMD+=($NOPOST)
echo "chain $CHAIN ($NAME): the first GPU attempt will cancel $JOBS when it starts"
echo "${CMD[*]} --submit"
if [ $SUBMIT = 1 ]; then
  "${CMD[@]}" --submit
  printf '{"stage": "gpu_attempts_submitted", "rc": 0, "device": "gpu", "host": "%s", "job": "", "utc": "%s", "note": "localnn_chain_to_gpu.sh: takeover of %s at the first GPU attempt start"}\n' \
    "$(hostname)" "$(date -u +%FT%TZ)" "$JOBS" >> "$CHAIN/chain_events.jsonl"
fi
