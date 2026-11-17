#!/bin/bash
# ==== version-controlled copy: header added 2026-11-15, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/tests/localnn_20261114/jacfail/replace_pending.sh (on Voyager: /data001/highresolu_vt/cops_output/tests/localnn_20261114/jacfail/replace_pending.sh)
# SHA256 of the run copy = SHA256 of this file without these header lines: 804a92e2321ef376b3079d66bfaf400f8cc1c90bedc679f5845ae1b34375f2b0
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: replace the pending attempts and post job of a running localnn chain by attempts with another MODEL_EXT,
#   afternotok on the running attempt (settings from the chain's submission.json).
# RUN 2026-11-15: 2026-11-16: stage2r_f10_pm1e-4_tplon (135875-135882), stage2r_f1_pm1e-4_tplon (135883-135890),
#   stage2_f10_pm1e-4 (135891-135897).
# ==== end of header ====
# replace the PENDING attempts (and the post job) of a running localnn chain by attempts with another MODEL_EXT
# (the Jacobian-check fix localnn_s2jc), afternotok on the RUNNING attempt; the running attempt is not touched.
#   replace_pending.sh CHAIN NAME MODEL_EXT   (TEMPLATES, FT_EXTRA, NONNEG_CODE as in CHAIN/submission.json)
set -euo pipefail
CHAIN=$1; NAME=$2; MEXT=$3
ROOT=/data001/highresolu_vt/cops_output; PY=/data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda/.venv/bin/python
read -r JOBS LAMBDA RULE FM C1 <<< "$($PY -c "
import json; j=json.load(open('$CHAIN/submission.json')); s=j['settings']
c1 = ('--c1-fit ' + s['C1_FIT']) if s.get('C1_FIT') else ('--c1-from ' + s['C1_FROM'])
print(':'.join(j['jobs']), s['LAMBDA'], s['RULE'], s['FM'], c1.replace(' ', '='))")"
TPL=$($PY -c "import json; print(json.load(open('$CHAIN/submission.json'))['settings'].get('TEMPLATES', 'on'))")
FTX=$($PY -c "import json; print(json.load(open('$CHAIN/submission.json'))['settings'].get('FT_EXTRA', ''))")
NNC=$($PY -c "import json; print(json.load(open('$CHAIN/submission.json'))['settings'].get('NONNEG_CODE', ''))")
X=(--templates "$TPL"); [ -n "$FTX" ] && X+=(--ft-extra "$FTX"); [ -n "$NNC" ] && X+=(--nonneg-code "$NNC")
RUN=; PEND=
for jb in ${JOBS//:/ }; do st=$(squeue -h -j $jb -o %T 2>/dev/null || true)
  case $st in RUNNING) RUN="$RUN $jb";; PENDING) PEND="$PEND $jb";; esac; done
RUN=$(echo $RUN); echo "$CHAIN: running [$RUN] pending [$PEND]"
[ -n "$RUN" ] || { echo "no running attempt: not replaced"; exit 2; }
N=$(echo $PEND | wc -w); NATT=$(( N - 1 )); [ $NATT -ge 1 ] || NATT=1
scancel $PEND
for jb in $PEND; do echo "$(date -u +%FT%TZ) $jb scancel (pending attempt / post of $NAME replaced by attempts with MODEL_EXT=$MEXT: the Jacobian-check noise floor; afternotok:$RUN)" >> $ROOT/model_ext/logs/submissions.txt; done
# shellcheck disable=SC2086
bash $ROOT/jobs/localnn_chain_submit.sh --submit --chain $CHAIN --name $NAME --lambda $LAMBDA --rule $RULE --fm $FM ${C1/=/ } \
  --model-ext $MEXT --script $MEXT/localnn_chain.sbatch.copy --attempts $NATT --dep afternotok:$RUN \
  --exclude "node002,node[019-024]" --c1-note "remaining attempts replaced $(date -u +%FT%TZ): MODEL_EXT $MEXT (the Jacobian check's FD noise floor; the fit arithmetic unchanged); afternotok on the running attempt $RUN" "${X[@]}" | tail -2
