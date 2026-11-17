#!/bin/bash
# ==== version-controlled copy: header added 2026-11-15, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/tests/localnn_20261114/gpu_solver/cpu_ref.sh (on Voyager: /data001/highresolu_vt/cops_output/tests/localnn_20261114/gpu_solver/cpu_ref.sh)
# SHA256 of the run copy = SHA256 of this file without these header lines: 43b49cbc3148f6a692e748591169157fe4ce7558da4339d3af8940d531a75afa
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: GPU solver path, CPU part: the default cpu path bitwise vs localnn_s2_20261115T185842Z (6-evaluation fit)
#   and the 50-evaluation CPU reference trajectory (localnn_gpu_solver.py trajectory --device cpu).
# RUN 2026-11-15: 132047 (the bitwise part had a wrong start path; rerun as bitwise.sh 132049).
# ==== end of header ====
# GPU solver path, CPU part (job on a CPU node): (a) default cpu path bitwise = the frozen stage-2 snapshot (6-evaluation
# fit from the same start, both codes); (b) the 50-evaluation CPU reference trajectory of the case
set -uo pipefail
export NPY_DISABLE_CPU_FEATURES="AVX512F AVX512CD AVX512_KNL AVX512_KNM AVX512_SKX AVX512_CLX AVX512_CNL AVX512_ICL"
export OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 MKL_NUM_THREADS=16 PYTHONUNBUFFERED=1
export PYTHONPYCACHEPREFIX=/data001/highresolu_vt/cops_output/tests/localnn_20261114/gpu_solver/pycache_$SLURM_JOB_ID
APP=/data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda; PY=$APP/.venv/bin/python; cd $APP
W=/data001/highresolu_vt/cops_output/tests/localnn_20261114/gpu_solver
R=/data001/highresolu_vt/cops_output/runs/shap8_012m_lategrid_dz1.5; L=$R/nonneg/localnn; C=$L/stage2_f10_pm0.01_20261115T190019Z
NEW=/data001/highresolu_vt/cops_output/model_ext; OLD=/data001/highresolu_vt/cops_output/model_ext_frozen/localnn_s2_20261115T185842Z
FT=$(readlink -f $C/c1/ft/final)
ARGS=(--design $L/grid_150m/design --base $R/base --d2 $R/clip_refit/D2/final --point-table $R/late_grid_points.txt --grid $L/grid_150m
      --templates on --slope-mult-min 1 --start-file $C/c2/fit/start --smooth-lambda 8.249495752351283e-05
      --smooth-rule-file $L/lambda_rule_anchored_20261114.json --mult-penalty-f 0.01 --held $FT/finetune_state.npz
      --late-kernels $L/late_kernels_150m.npz --max-nfev 6)
$PY -W ignore $OLD/localnn_fit.py fit "${ARGS[@]}" --out $W/bitwise_old > $W/bitwise_old.log 2>&1; echo "old rc $?"
$PY -W ignore $NEW/localnn_fit.py fit "${ARGS[@]}" --out $W/bitwise_new > $W/bitwise_new.log 2>&1; echo "new rc $?"
$PY -c "
import numpy as np, json
a=np.load('$W/bitwise_old/final/localnn_state.npz'); b=np.load('$W/bitwise_new/final/localnn_state.npz')
ma, mb = json.loads(str(a['meta_json'])), json.loads(str(b['meta_json']))
print('default cpu path bitwise:', np.array_equal(a['x'], b['x']), 'objective', ma['costs']['objective'], mb['costs']['objective'], 'device key in new meta:', 'device' in mb)
"
$PY -W ignore $NEW/localnn_gpu_solver.py trajectory --device cpu --case $W/case_c2_f10_pm0.01.json --nfev 50 --out $W/cpu_ref_50 2>&1 | grep -v "eval " | tail -2
