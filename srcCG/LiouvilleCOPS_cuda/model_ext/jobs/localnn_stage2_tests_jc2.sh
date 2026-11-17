#!/bin/bash
# ==== version-controlled copy: header added 2026-11-15, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/tests/localnn_20261114/stage2/tests_s2_jc2.sh (on Voyager: /data001/highresolu_vt/cops_output/tests/localnn_20261114/stage2/tests_s2_jc2.sh)
# SHA256 of the run copy = SHA256 of this file without these header lines: 63271a7dfff79d02b29ea41839ced5696032d3d962450b8b8a00c3b8be1aa159
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: test_localnn_stage2.py from the frozen snapshot localnn_s2jc2_20261117T080730Z.
# RUN 2026-11-15: 136043.
# ==== end of header ====
export NPY_DISABLE_CPU_FEATURES="AVX512F AVX512CD AVX512_KNL AVX512_KNM AVX512_SKX AVX512_CLX AVX512_CNL AVX512_ICL"
export OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 PYTHONUNBUFFERED=1 APP=/data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda
cd $APP
.venv/bin/python -W ignore /data001/highresolu_vt/cops_output/model_ext_frozen/localnn_s2jc2_20261117T080730Z/test_localnn_stage2.py --grid /data001/highresolu_vt/cops_output/runs/shap8_012m_lategrid_dz1.5/nonneg/localnn/grid_150m --base /data001/highresolu_vt/cops_output/runs/shap8_012m_lategrid_dz1.5/base --d2 /data001/highresolu_vt/cops_output/runs/shap8_012m_lategrid_dz1.5/clip_refit/D2/final --point-table /data001/highresolu_vt/cops_output/runs/shap8_012m_lategrid_dz1.5/late_grid_points.txt --late /data001/highresolu_vt/cops_output/tests/localnn_20261114/stage2/late_kernels_dev.npz --contrib-code /data001/highresolu_vt/cops_output/model_ext_frozen/nonneg_contrib_20261112T210012Z --state /data001/highresolu_vt/cops_output/runs/shap8_012m_lategrid_dz1.5/nonneg/localnn/stage1_anch_f10_pm0.01_cont_20261115T020243Z/fit/seg01 --rule /data001/highresolu_vt/cops_output/runs/shap8_012m_lategrid_dz1.5/nonneg/localnn/lambda_rule_anchored_20261114.json --work /data001/highresolu_vt/cops_output/tests/localnn_20261114/stage2/tests_s2_jc2_work --out /data001/highresolu_vt/cops_output/tests/localnn_20261114/stage2/tests_s2_jc2.json
