#!/bin/bash
# ==== version-controlled copy: header added 2026-11-15, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/tests/localnn_20261114/ftreg/ref_summary.sh (on Voyager: /data001/highresolu_vt/cops_output/tests/localnn_20261114/ftreg/ref_summary.sh)
# SHA256 of the run copy = SHA256 of this file without these header lines: 12812580a455d7a95b7f981079ee2a6d589f9e35c8efacf7a64fbae0c088e393
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: reference numbers of the unregularized c5 fine tunings (C f1 / 0.01, C f10 / 0.01, B, D) with
#   localnn_ftreg_summary.py.
# RUN 2026-11-15: 134496.
# ==== end of header ====
# fine-tuning regularization tests: the reference numbers of the UNREGULARIZED c5 fine tunings (C f1/0.01, C f10/0.01,
# B bg1e-2, D logspread f0.1) with localnn_ftreg_summary.py (also a test of the summary script)
set -uo pipefail
export NPY_DISABLE_CPU_FEATURES="AVX512F AVX512CD AVX512_KNL AVX512_KNM AVX512_SKX AVX512_CLX AVX512_CNL AVX512_ICL"
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 PYTHONUNBUFFERED=1
W=/data001/highresolu_vt/cops_output/tests/localnn_20261114/ftreg
export PYTHONPYCACHEPREFIX=$W/pycache_$SLURM_JOB_ID
APP=/data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda; cd $APP
M=${MODEL_EXT:-/data001/highresolu_vt/cops_output/model_ext}
R=/data001/highresolu_vt/cops_output/runs/shap8_012m_lategrid_dz1.5; N=$R/nonneg; L=$N/localnn
CF1=$(ls -d $L/stage2_f1_pm0.01_2*/ | head -1); CF10=$(ls -d $L/stage2_f10_pm0.01_2*/ | head -1)
.venv/bin/python -W ignore $M/localnn_ftreg_summary.py --late $L/late_kernels_150m.npz --base $R/base --d2 $R/clip_refit/D2/final \
  --point-table $R/late_grid_points.txt --design $L/grid_150m/design --geometry $R/maps_geometry/job_124775/geometry \
  --contrib-code /data001/highresolu_vt/cops_output/model_ext_frozen/nonneg_contrib_20261112T210012Z \
  --case C_f1_pm0.01_ref=${CF1}c5/ft/final --case C_f10_pm0.01_ref=${CF10}c5/ft/final \
  --case B_bg1e-2_ref=$N/background/stage2_bg1e-2_20261115T184949Z/c5/ft/final \
  --case D_logspread_f0.1_ref=$N/logspread/stage2_f0.1_20261115T184949Z/c5/ft/final \
  --out $W/summary_ref_$SLURM_JOB_ID.json --txt $W/summary_ref_$SLURM_JOB_ID.txt
