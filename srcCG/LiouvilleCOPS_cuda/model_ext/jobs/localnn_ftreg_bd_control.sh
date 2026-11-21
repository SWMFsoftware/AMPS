#!/bin/bash
# ==== version-controlled copy: header added 2026-11-16, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/tests/localnn_20261114/ftreg/bd_control.sh (on Voyager: /data001/highresolu_vt/cops_output/tests/localnn_20261114/ftreg/bd_control.sh)
# SHA256 of the run copy = SHA256 of this file without these header lines: b06dd3d63a92c13fa8d23b07a243d9aeb781c3caaba969b384227235b4b89589
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: CPU engine control: the B and D c5 fine tunings unregularized on CPU (nonneg_20261115T020314Z) vs their GPU
#   runs.
# RUN 2026-11-16: 134501 (B rel 6.3e-14, D 4.4e-16; same evaluations).
# ==== end of header ====
# CPU engine on the B and D c5 fits (agent 3's chains), UNREGULARIZED, with the nonneg code that ran those chains
# (model_ext_frozen/nonneg_20261115T020314Z): compare with their GPU c5 fine tuning before the regularization tests
set -uo pipefail
export NPY_DISABLE_CPU_FEATURES="AVX512F AVX512CD AVX512_KNL AVX512_KNM AVX512_SKX AVX512_CLX AVX512_CNL AVX512_ICL"
export OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 MKL_NUM_THREADS=16 PYTHONUNBUFFERED=1
W=/data001/highresolu_vt/cops_output/tests/localnn_20261114/ftreg
export PYTHONPYCACHEPREFIX=$W/pycache_$SLURM_JOB_ID
export LOCALNN_NONNEG_CODE=/data001/highresolu_vt/cops_output/model_ext_frozen/nonneg_20261115T020314Z
APP=/data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda; PY=$APP/.venv/bin/python; cd $APP
M=/data001/highresolu_vt/cops_output/model_ext; N=/data001/highresolu_vt/cops_output/runs/shap8_012m_lategrid_dz1.5/nonneg
P=$W/provenance-$SLURM_JOB_ID-bdcontrol; mkdir -p $P; cp -p "$0" $P/; cp -p $M/localnn_ft.py $M/localnn_late.py $P/
for spec in "B background/stage2_bg1e-2_20261115T184949Z B_c5" "D logspread/stage2_f0.1_20261115T184949Z D_c5"; do
  set -- $spec
  echo "$(date -u +%FT%TZ) $1"
  $PY -W ignore $M/localnn_ft.py run-nonneg --fit $N/$2/c5/fit/seg01 --late $W/late_kernels_BD.npz --state-name $3 \
    --facets late_only --window-min 5 --max-nfev 500 --lsmr-maxiter 50 --ftol 1e-10 --xtol 1e-10 --gtol 1e-10 \
    --x-scale-grain 1.0 --validate --seed 20261110 --save-seconds 900 --out $W/${1}_control_cpu_oldcode > $P/${1}.log 2>&1
  echo "  rc $?"; grep -h "nonneg_finetune: final\|validation" $P/${1}.log | tail -2 | cut -c1-300
  $PY -c "
import json
a=json.load(open('$W/${1}_control_cpu_oldcode/fine_tuning.json')); b=json.load(open('$N/$2/c5/ft/final/fine_tuning.json'))
for k in ('cost_total','cost_data'): print('  ', k, 'CPU', a['costs_final'][k], 'GPU', b['costs_final'][k], 'rel', abs(a['costs_final'][k]/b['costs_final'][k]-1), 'start rel', abs(a['costs_initial'][k]/b['costs_initial'][k]-1))
print('   nfev', a['optimizer']['nfev'], b['optimizer']['nfev'], a['optimizer']['termination'], b['optimizer']['termination'])
"
done
