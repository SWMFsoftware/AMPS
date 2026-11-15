#!/bin/bash
# ==== version-controlled copy: header added 2026-11-14, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/tests/localnn_20261114/submit_anchored.sh (on Voyager: /data001/highresolu_vt/cops_output/tests/localnn_20261114/submit_anchored.sh)
# SHA256 of the run copy = SHA256 of this file without these header lines: 0b0a98f70db72966ce82ab543f6c6f5ad82f33208e2db040445c6983cda13fe1
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: submission helper (login node, option C): the three solution-anchored stage-1 fits (f 0.1, 1, 10; lambda from
#   the rule file), each with its continuation and post-fit chain (tools, report, GPU evaluation, docs); logs every
#   submission in model_ext/logs/submissions.txt.
# RUN: 2026-11-14 09:42 UTC -> jobs 125107-125124.
# ==== end of header ====
# Submits the three solution-anchored stage-1 fits (option C) with their continuations and post-fit chains.
# Usage: submit_anchored.sh <rule.json>   (run on the login node; logs every submission in submissions.txt)
set -euo pipefail
RULE=$(readlink -f "$1")
cd /data001/highresolu_vt/cops_output/jobs
AX=/data001/highresolu_vt/cops_output/model_ext_frozen/localnn_anch_20261114T090437Z
L=/data001/highresolu_vt/cops_output/runs/shap8_012m_lategrid_dz1.5/nonneg/localnn
LOG=/data001/highresolu_vt/cops_output/model_ext/logs/submissions.txt
PY=/data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda/.venv/bin/python
set -- node001:0.1 node002:1 node016:10
for spec in "$@"; do
  n=${spec%%:*}; f=${spec#*:}
  LAM=$($PY -c "import json,sys; print(repr(json.load(open('$RULE'))['lambda_C']['$f']))")
  RUN=stage1_anch_f$f; RR=$L/${RUN}_$(date -u +%Y%m%dT%H%M%SZ); mkdir "$RR"
  E="ALL,MODEL_EXT=$AX,TEMPLATES=on,SMOOTH_LAMBDA=$LAM,SMOOTH_RULE_FILE=$RULE,NFEV=4000,RUN=$RUN,RUNROOT=$RR"
  J1=$(sbatch --parsable --time=12:00:00 --nodelist=$n --export=$E localnn_fit.sbatch)
  J2=$(sbatch --parsable --time=12:00:00 --nodelist=$n --dependency=afterany:$J1 --export=$E localnn_fit.sbatch)
  T=$(sbatch --parsable --dependency=afterok:$J2 --export=ALL,MODEL_EXT=$AX,RUNROOT=$RR localnn_tools.sbatch)
  R=$(sbatch --parsable --dependency=afterok:$J2 --export=ALL,MODEL_EXT=$AX,RUNROOT=$RR localnn_report.sbatch)
  G=$(sbatch --parsable --dependency=afterok:$J2 --export=ALL,MODEL_EXT=$AX,SEG=$RR/fit/final localnn_gpu_eval.sbatch)
  D=$(sbatch --parsable --dependency=afterok:$R,afterany:$G --export=ALL,MODEL_EXT=$AX,RUNROOT=$RR,DOCS_ONLY=1 localnn_report.sbatch)
  now=$(date -u +%FT%TZ)
  {
  echo "$now $J1 localnn_fit.sbatch normal128cores --nodelist=$n 16 cpus --time=12:00:00 (option C STAGE-1, SOLUTION-ANCHORED lambda, user ~02:55 CST: TEMPLATES=on SMOOTH_LAMBDA=$LAM (f=$f, rule $RULE) NFEV=4000; frozen localnn_anch_20261114T090437Z) -> $RR"
  echo "$now $J2 localnn_fit.sbatch --nodelist=$n --dependency=afterany:$J1 (continuation segment, same settings; no-op if $J1 completed) -> $RR"
  echo "$now $T localnn_tools.sbatch afterok:$J2 (f=$f: export, contrib, maps) -> $RR/{tools_state,contributions,maps}"
  echo "$now $R localnn_report.sbatch afterok:$J2 (f=$f: item-5 numbers, bins table, README.md, HOW-TO-RERUN.md) -> $RR/analysis"
  echo "$now $G localnn_gpu_eval.sbatch gpu002 afterok:$J2 (f=$f: GPU standard evaluation, D2 path) -> $RR/gpu_eval"
  echo "$now $D localnn_report.sbatch DOCS_ONLY=1 afterok:$R,afterany:$G (f=$f: README.md / HOW-TO-RERUN.md with the GPU check)"
  } >> $LOG
  echo "f=$f lambda=$LAM node=$n fit=$J1 cont=$J2 tools=$T report=$R gpu=$G docs=$D runroot=$RR"
  sleep 1
done
