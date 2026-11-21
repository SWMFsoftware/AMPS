#!/bin/bash
# ==== version-controlled copy: header added 2026-11-18, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/jobs/localnn_chain_submit.sh (on Voyager: /data001/highresolu_vt/cops_output/jobs/localnn_chain_submit.sh)
# SHA256 of the run copy = SHA256 of this file without these header lines: 2848928a9cb411ec0f950052d376b5f11c6f212e4b9dafa411123458e11cf6a2
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: submit one stage-2 chain: attempts (afterok / --dep / afternotok), --sbatch-extra, --env, --script, --device
#   / --gpu-svd-mode, --templates, --c1-start-extra, --ft-extra, --nonneg-code, (2026-11-18) --c1-start d2, --base-
#   dir, --d2-dir, --tools-steps, --design-dir (DESIGN_DIR of the attempts and the post job); the post job
#   (afterany); submission.json, README / HOW-TO-RERUN, submissions.txt.
# RUN 2026-11-18: planned: the q002 option C chains with --design-dir grid_150m/design_h2oqmin0.02_d2q002 (dry run
#   2026-11-18).
# ==== end of header ====
# Option C stage 2 (lead 2026-11-15): create ONE localnn chain root and submit its attempts (jobs/localnn_chain.sbatch)
# and its post-fit job (jobs/localnn_chain_post.sbatch); writes <CHAIN>/submission.json, README.md, HOW-TO-RERUN.md
# (model_ext/localnn_chain_docs.py) and one line per job in model_ext/logs/submissions.txt.
#
#   localnn_chain_submit.sh --name NAME --lambda L --rule RULE --fm FM (--c1-from STATE.npz | --c1-fit FITDIR)
#       [--c1-note TEXT] [--chain DIR (an existing root: resume)] [--attempts 7] [--time 12:00:00]
#       [--partition normal128cores] [--exclude node002,node023] [--first-node NODE] [--after JOBID (afterok)]
#       [--model-ext DIR] [--late FILE] [--device cpu] [--cycles 5] [--nfev 4000] [--ft-maxnfev 500]
#       [--window-min 5] [--checkpoint 250] [--check 1] [--blas-threads 16] [--no-post] [--submit]
#       [--sbatch-extra "OPTIONS"] (added to every attempt's sbatch line, e.g. the GPU options with --device gpu)
#       [--gpu-svd-mode qr|gesvd] [--dep DEPENDENCY (raw, for attempt 1, e.g. afterany:<job>; with --after: both)]
#       [--script FILE (default jobs/localnn_chain.sbatch; e.g. the frozen snapshot's localnn_chain.sbatch.copy)]
#       [--templates on|off] [--c1-start-extra "ARGS"] [--ft-extra "ARGS"] [--nonneg-code DIR] (2026-11-16)
#       [--env "K=V[,K=V]"]
#       [--c1-start d2] [--base-dir DIR] [--d2-dir DIR] [--tools-steps "export contrib"] (2026-11-18: cycle 1 from the
#       D2 projection of D2_DIR; the base / D2 of the new law; the post's tool steps, e.g. without the maps) (added to every attempt's --export, e.g. CANCEL_JOBS=j1:j2:... for a device takeover)
#       [--design-dir DIR] (2026-11-18: DESIGN_DIR of the attempts AND the post job; default: the driver's / post scripts'
#       default grid_150m/design_h2oqmin0.02; the q002 chains use grid_150m/design_h2oqmin0.02_d2q002, lead 2026-11-18)
#
# Attempt 1: --dependency=afterok:<--after> (if given), --nodelist=<--first-node> (if given); attempt k > 1:
# afternotok:<attempt k-1> (it runs only when attempt k-1 ended non-zero, e.g. exit 75 at the wall-time guard), with
# --kill-on-invalid-dep=yes (the remaining attempts are removed once an attempt completes the chain). The post job:
# afterany:<all attempts>. Without --submit: dry run (prints the sbatch lines, writes nothing).
set -euo pipefail
ROOT=/data001/highresolu_vt/cops_output
L=$ROOT/runs/shap8_012m_lategrid_dz1.5/nonneg/localnn
NAME= LAMBDA= RULE= FM= C1_FROM= C1_FIT= C1_NOTE= CHAIN= ATT=7 TIME=12:00:00 PART=normal128cores EXCL=node002,node023
FIRST= AFTER= MEXT=$ROOT/model_ext LATE=$L/late_kernels_150m_h2oqmin0.02.npz DEVICE=cpu CYC=5 NFEV=4000 FTM=500 WMIN=5 CKP=250
CHECK=1 BT=16 POST=1 SUBMIT=0 EXTRA= SVDM=qr DEP= XENV= SCRIPT=$ROOT/jobs/localnn_chain.sbatch TPL=on C1X= FTX= NNC= C1S=state BDIR= DDIR= TSTEPS= DSGN=
CMDLINE="$0 $*"
while [ $# -gt 0 ]; do
  case $1 in
    --name) NAME=$2; shift;; --lambda) LAMBDA=$2; shift;; --rule) RULE=$2; shift;; --fm) FM=$2; shift;;
    --c1-from) C1_FROM=$2; shift;; --c1-fit) C1_FIT=$2; shift;; --c1-note) C1_NOTE=$2; shift;;
    --chain) CHAIN=$2; shift;; --attempts) ATT=$2; shift;; --time) TIME=$2; shift;; --partition) PART=$2; shift;;
    --exclude) EXCL=$2; shift;; --first-node) FIRST=$2; shift;; --after) AFTER=$2; shift;;
    --model-ext) MEXT=$2; shift;; --late) LATE=$2; shift;; --device) DEVICE=$2; shift;; --cycles) CYC=$2; shift;;
    --nfev) NFEV=$2; shift;; --ft-maxnfev) FTM=$2; shift;; --window-min) WMIN=$2; shift;; --checkpoint) CKP=$2; shift;;
    --check) CHECK=$2; shift;; --blas-threads) BT=$2; shift;; --no-post) POST=0;; --submit) SUBMIT=1;;
    --sbatch-extra) EXTRA=$2; shift;; --gpu-svd-mode) SVDM=$2; shift;; --dep) DEP=$2; shift;; --env) XENV=$2; shift;; --templates) TPL=$2; shift;; --c1-start-extra) C1X=$2; shift;;
    --ft-extra) FTX=$2; shift;; --nonneg-code) NNC=$(readlink -f "$2"); shift;;
    --c1-start) C1S=$2; shift;; --base-dir) BDIR=$(readlink -f "$2"); shift;; --d2-dir) DDIR=$(readlink -f "$2"); shift;;
    --tools-steps) TSTEPS=$2; shift;; --design-dir) DSGN=$(readlink -f "$2"); shift;; --script) SCRIPT=$(readlink -f "$2"); shift;;
    *) echo "unknown option $1"; exit 2;;
  esac; shift
done
[ -n "$NAME" ] && [ -n "$LAMBDA" ] && [ -n "$RULE" ] && [ -n "$FM" ] || { echo "--name --lambda --rule --fm are required"; exit 2; }
[ -n "$C1_FROM$C1_FIT" ] || [ "$C1S" = d2 ] || { echo "--c1-from, --c1-fit or --c1-start d2 is required"; exit 2; }
MEXT=$(readlink -f "$MEXT"); RULE=$(readlink -f "$RULE")
U=$(date -u +%Y%m%dT%H%M%SZ)
[ -n "$CHAIN" ] || CHAIN=$L/stage2_${NAME}_$U
EXP="ALL,CHAIN=$CHAIN,LAMBDA=$LAMBDA,RULE=$RULE,FM=$FM,CYCLES=$CYC,NFEV=$NFEV,FT_MAXNFEV=$FTM,WINDOW_MIN=$WMIN,CHECKPOINT=$CKP,LATE=$LATE,DEVICE=$DEVICE,MODEL_EXT=$MEXT,BLAS_THREADS=$BT,CHECK=$CHECK,GPU_SVD_MODE=$SVDM"
[ -n "$C1_FROM" ] && EXP="$EXP,C1_FROM=$(readlink -f "$C1_FROM")"
[ -n "$C1_FIT" ] && EXP="$EXP,C1_FIT=$(readlink -f "$C1_FIT")"
EXP="$EXP,TEMPLATES=$TPL,C1_START=$C1S"
[ -n "$BDIR" ] && EXP="$EXP,BASE_DIR=$BDIR"
[ -n "$DDIR" ] && EXP="$EXP,D2_DIR=$DDIR"
[ -n "$DSGN" ] && { [ -f "$DSGN/design.npz" ] || { echo "--design-dir $DSGN has no design.npz"; exit 2; }; EXP="$EXP,DESIGN_DIR=$DSGN"; }
[ -n "$C1X" ] && EXP="$EXP,C1_START_EXTRA=$C1X"
[ -n "$FTX" ] && EXP="$EXP,FT_EXTRA=$FTX"
[ -n "$NNC" ] && EXP="$EXP,NONNEG_CODE=$NNC"
case "$C1X$FTX" in *,*) echo "--c1-start-extra / --ft-extra must not contain commas"; exit 2;; esac
[ -n "$XENV" ] && EXP="$EXP,$XENV"
COMMON=(--parsable --kill-on-invalid-dep=yes --partition="$PART" --time="$TIME" --export="$EXP")
[ -n "$EXCL" ] && COMMON+=(--exclude="$EXCL")
# shellcheck disable=SC2206
[ -n "$EXTRA" ] && COMMON+=($EXTRA)
LINES=(); JOBS=(); PREV=
for a in $(seq 1 "$ATT"); do
  O=(--job-name="lnn_${NAME}_$a")
  D=
  [ -n "$AFTER" ] && D="afterok:$AFTER"
  [ -n "$DEP" ] && D="${D:+$D,}$DEP"
  [ -n "$PREV" ] && D="${D:+$D,}afternotok:$PREV"
  [ -n "$D" ] && O+=(--dependency="$D")
  [ "$a" = 1 ] && [ -n "$FIRST" ] && O+=(--nodelist="$FIRST")
  LINE="sbatch ${COMMON[*]} ${O[*]} $SCRIPT"
  if [ $SUBMIT = 1 ]; then
    [ "$a" = 1 ] && { mkdir -p "$CHAIN"; }
    J=$(sbatch "${COMMON[@]}" "${O[@]}" "$SCRIPT")
    echo "$(date -u +%FT%TZ) $J $SCRIPT $PART localnn stage-2 chain $NAME attempt $a: ${COMMON[*]:1} ${O[*]}" >> "$ROOT/model_ext/logs/submissions.txt"
  else J="<attempt$a>"; fi
  LINES+=("$LINE"); JOBS+=("$J"); PREV=$J
  echo "attempt $a: $J"; [ $SUBMIT = 1 ] || echo "  $LINE"
done
if [ $POST = 1 ]; then
  DEP="afterany:$(IFS=:; echo "${JOBS[*]}")"
  PO=(--parsable --partition=normal128cores --exclude="$EXCL" --job-name="lnn_${NAME}_post" --dependency="$DEP"
      --export="ALL,CHAIN=$CHAIN,CYCLES=$CYC,MODEL_EXT=$MEXT${BDIR:+,BASE_DIR=$BDIR}${DDIR:+,D2_DIR=$DDIR}${DSGN:+,DESIGN_DIR=$DSGN}${TSTEPS:+,TOOLS_STEPS=$TSTEPS}")
  LINE="sbatch ${PO[*]} $ROOT/jobs/localnn_chain_post.sbatch"
  if [ $SUBMIT = 1 ]; then
    J=$(sbatch "${PO[@]}" "$ROOT/jobs/localnn_chain_post.sbatch")
    echo "$(date -u +%FT%TZ) $J $ROOT/jobs/localnn_chain_post.sbatch normal128cores localnn stage-2 chain $NAME post (degeneracy, report, tools on the final ft state, docs): ${PO[*]:1}" >> "$ROOT/model_ext/logs/submissions.txt"
  else J="<post>"; fi
  LINES+=("$LINE"); JOBS+=("$J"); echo "post: $J"; [ $SUBMIT = 1 ] || echo "  $LINE"
fi
if [ $SUBMIT = 1 ]; then
  /data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda/.venv/bin/python - "$CHAIN/submission.json" "$U" "$CMDLINE" "$C1_NOTE" "$EXP" \
      "$(git -C /data001/highresolu_vt/AMPS rev-parse HEAD)" "${#LINES[@]}" "${LINES[@]}" "${JOBS[@]}" <<'EOF'
import json, sys, os, time
p, u, cmd, note, exp, amps, n = sys.argv[1:8]
n = int(n); lines, jobs = sys.argv[8:8 + n], sys.argv[8 + n:8 + 2 * n]
st = dict(kv.split("=", 1) for kv in exp.split(",")[1:])
old = json.load(open(p)) if os.path.exists(p) else None
rec = {"submitted_utc": u, "command": cmd, "c1_note": note, "settings": st, "amps_commit": amps, "sbatch_lines": lines,
       "jobs": jobs, "earlier_submissions": ([old] + old.get("earlier_submissions", [])) if old else []}
if old:
    old.pop("earlier_submissions", None)
json.dump(rec, open(p + ".tmp", "w"), indent=1)
os.replace(p + ".tmp", p)
EOF
  /data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda/.venv/bin/python "$MEXT/localnn_chain_docs.py" --chain "$CHAIN"
fi
echo "CHAIN=$CHAIN"
