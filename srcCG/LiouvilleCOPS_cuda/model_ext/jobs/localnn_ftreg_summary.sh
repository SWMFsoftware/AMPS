#!/bin/bash
# ==== version-controlled copy: header added 2026-11-15, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/tests/localnn_20261114/ftreg/ftreg_summary.sh (on Voyager: /data001/highresolu_vt/cops_output/tests/localnn_20261114/ftreg/ftreg_summary.sh)
# SHA256 of the run copy = SHA256 of this file without these header lines: 77d9fc9659d056f1eff75548c5bac5933b6983a27a8ae2fabfe43a71e10cbe89
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: the per-run summary of the fine-tuning regularization tests (C, B, D; regularized runs, CPU controls and the
#   chains' c5 fine tunings as references) with localnn_ftreg_summary.py of localnn_s2ftreg.
# RUN 2026-11-15: 134507.
# ==== end of header ====
# the per-run numbers of the fine-tuning regularization tests (localnn_ftreg_summary.py of localnn_s2ftreg)
set -uo pipefail
export NPY_DISABLE_CPU_FEATURES="AVX512F AVX512CD AVX512_KNL AVX512_KNM AVX512_SKX AVX512_CLX AVX512_CNL AVX512_ICL"
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 PYTHONUNBUFFERED=1
W=/data001/highresolu_vt/cops_output/tests/localnn_20261114/ftreg
APP=/data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda; cd $APP
M=/data001/highresolu_vt/cops_output/model_ext_frozen/localnn_s2ftreg_20261116T031038Z
R=/data001/highresolu_vt/cops_output/runs/shap8_012m_lategrid_dz1.5; N=$R/nonneg; L=$N/localnn
CF1=$(ls -d $L/stage2_f1_pm0.01_2*/ | head -1); CF10=$(ls -d $L/stage2_f10_pm0.01_2*/ | head -1)
CASES=(--case C_f1_pm0.01_ref=${CF1}c5/ft/final)
for f in 0.1 1 10; do CASES+=(--case C_f1_pm0.01_fd${f}_cap=$W/C_f1_pm0.01_fd${f}_cap); done
CASES+=(--case C_f10_pm0.01_ref=${CF10}c5/ft/final)
for f in 0.1 1 10; do CASES+=(--case C_f10_pm0.01_fd${f}_cap=$W/C_f10_pm0.01_fd${f}_cap); done
CASES+=(--case B_ref_gpu=$N/background/stage2_bg1e-2_20261115T184949Z/c5/ft/final --case B_control_cpu=$W/B_control_cpu)
for f in 0.1 1 10; do CASES+=(--case B_fd${f}_cap=$W/B_fd${f}_cap); done
CASES+=(--case D_ref_gpu=$N/logspread/stage2_f0.1_20261115T184949Z/c5/ft/final --case D_control_cpu=$W/D_control_cpu)
for f in 0.1 1 10; do CASES+=(--case D_fd${f}_cap=$W/D_fd${f}_cap); done
.venv/bin/python -W ignore $M/localnn_ftreg_summary.py --late $L/late_kernels_150m.npz --base $R/base --d2 $R/clip_refit/D2/final \
  --point-table $R/late_grid_points.txt --design $L/grid_150m/design --geometry $R/maps_geometry/job_124775/geometry \
  --contrib-code /data001/highresolu_vt/cops_output/model_ext_frozen/nonneg_contrib_20261112T210012Z "${CASES[@]}" \
  --out $W/summary_$SLURM_JOB_ID.json --txt $W/summary_$SLURM_JOB_ID.txt
