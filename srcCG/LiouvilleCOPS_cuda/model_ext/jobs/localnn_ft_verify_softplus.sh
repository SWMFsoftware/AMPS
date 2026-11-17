#!/bin/bash
# ==== version-controlled copy: header added 2026-11-15, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/tests/localnn_20261114/stage2/verify_sp.sh (on Voyager: /data001/highresolu_vt/cops_output/tests/localnn_20261114/stage2/verify_sp.sh)
# SHA256 of the run copy = SHA256 of this file without these header lines: cde9fe607f167dfd37f2e21012996d821b1f46e9da57ce459b277212c1d6aee9
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: CPU fine tuning (localnn_ft.py verify-softplus) vs the GPU nonneg_finetune of the softplus chain
#   chain_5cycles_smm1_20261113T041215Z c5 and c1 (same facets, delta, cost).
# RUN 2026-11-15: 131977 (refused: AVX-512 exp), 131978.
# ==== end of header ====
export NPY_DISABLE_CPU_FEATURES="AVX512F AVX512CD AVX512_KNL AVX512_KNM AVX512_SKX AVX512_CLX AVX512_CNL AVX512_ICL"
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 PYTHONUNBUFFERED=1 APP=/data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda
cd $APP
for k in 5 1; do
  .venv/bin/python -W ignore /data001/highresolu_vt/cops_output/model_ext/localnn_ft.py verify-softplus --fit /data001/highresolu_vt/cops_output/runs/shap8_012m_lategrid_dz1.5/nonneg/softplus/chain_5cycles_smm1_20261113T041215Z/c$k/fit/seg01 --gpu-ft /data001/highresolu_vt/cops_output/runs/shap8_012m_lategrid_dz1.5/nonneg/softplus/chain_5cycles_smm1_20261113T041215Z/c$k/ft/final --late /data001/highresolu_vt/cops_output/tests/localnn_20261114/stage2/late_kernels_dev.npz --state-name softplus_c$k --out /data001/highresolu_vt/cops_output/tests/localnn_20261114/stage2/verify_softplus_c$k; echo "c$k rc $?"
done
