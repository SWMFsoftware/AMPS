#!/bin/bash
# ==== version-controlled copy: header added 2026-11-14, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/tests/localnn_20261114/degeneracy_run.sh (on Voyager: /data001/highresolu_vt/cops_output/tests/localnn_20261114/degeneracy_run.sh)
# SHA256 of the run copy = SHA256 of this file without these header lines: ebb948ebf40ea54db303ee8005000d95b44abce67217e399f9e12875e4deff10
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: helper (option C): runs localnn_degeneracy.py on the anchored f10 and f0.1 run roots (srun, CPU).
# RUN: 2026-11-14 22:26 UTC (srun job on normal128cores) -> <run root>/analysis_degeneracy.
# ==== end of header ====
export OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 PYTHONUNBUFFERED=1
cd /data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda
for RR in /data001/highresolu_vt/cops_output/runs/shap8_012m_lategrid_dz1.5/nonneg/localnn/stage1_anch_f10_20261114T094259Z /data001/highresolu_vt/cops_output/runs/shap8_012m_lategrid_dz1.5/nonneg/localnn/stage1_anch_f0.1_20261114T094256Z; do
  .venv/bin/python -W ignore /data001/highresolu_vt/cops_output/model_ext/localnn_degeneracy.py --runroot $RR --design /data001/highresolu_vt/cops_output/runs/shap8_012m_lategrid_dz1.5/nonneg/localnn/grid_150m/design --grid /data001/highresolu_vt/cops_output/runs/shap8_012m_lategrid_dz1.5/nonneg/localnn/grid_150m --base /data001/highresolu_vt/cops_output/runs/shap8_012m_lategrid_dz1.5/base --d2 /data001/highresolu_vt/cops_output/runs/shap8_012m_lategrid_dz1.5/clip_refit/D2/final --point-table /data001/highresolu_vt/cops_output/runs/shap8_012m_lategrid_dz1.5/late_grid_points.txt --c9-ft /data001/highresolu_vt/cops_output/runs/shap8_012m_lategrid_dz1.5/nonneg/softplus/chain_5cycles_smm1_20261113T041215Z/c9/ft/final; echo "rc $?"
done
