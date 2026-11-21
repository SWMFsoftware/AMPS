#!/bin/bash
# ==== version-controlled copy: header added 2026-11-16, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/tests/localnn_20261114/stage2r/share_c5.sh (on Voyager: /data001/highresolu_vt/cops_output/tests/localnn_20261114/stage2r/share_c5.sh)
# SHA256 of the run copy = SHA256 of this file without these header lines: 3309d6b85028a36b178ab2064eba299adbc41daddbf077d0d8abaee60024f7b9
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: late-signal share of the tuned facets and production numbers of the c5 fine tunings (unregularized f10 /
#   0.01, stage-2r templates on / off) with localnn_ftreg_summary.py.
# RUN 2026-11-16: the stage-2r comparison (2026-11-16).
# ==== end of header ====
# late-signal share of the tuned facets etc. for the c5 fine tunings of the stage-2r chains (localnn_ftreg_summary.py)
set -uo pipefail
export NPY_DISABLE_CPU_FEATURES="AVX512F AVX512CD AVX512_KNL AVX512_KNM AVX512_SKX AVX512_CLX AVX512_CNL AVX512_ICL"
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 PYTHONUNBUFFERED=1
W=/data001/highresolu_vt/cops_output/tests/localnn_20261114/stage2r
APP=/data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda; cd $APP
M=/data001/highresolu_vt/cops_output/model_ext_frozen/localnn_s2r_20261116T044801Z
R=/data001/highresolu_vt/cops_output/runs/shap8_012m_lategrid_dz1.5; L=$R/nonneg/localnn
.venv/bin/python -W ignore $M/localnn_ftreg_summary.py --late $L/late_kernels_150m.npz --base $R/base --d2 $R/clip_refit/D2/final \
  --point-table $R/late_grid_points.txt --design $L/grid_150m/design --geometry $R/maps_geometry/job_124775/geometry \
  --contrib-code /data001/highresolu_vt/cops_output/model_ext_frozen/nonneg_contrib_20261112T210012Z \
  --case unreg_c5=$L/stage2_f10_pm0.01_20261115T190019Z/c5/ft/final \
  --case tplon_c5=$L/stage2r_f10_pm0.01_tplon_20261116T044825Z/c5/ft/final \
  --case tploff_c5=$L/stage2r_f10_pm0.01_tploff_20261116T044825Z/c5/ft/final \
  --case tploff_c1=$L/stage2r_f10_pm0.01_tploff_20261116T044825Z/c1/ft/final \
  --out $W/share_c5_$SLURM_JOB_ID.json --txt $W/share_c5_$SLURM_JOB_ID.txt
