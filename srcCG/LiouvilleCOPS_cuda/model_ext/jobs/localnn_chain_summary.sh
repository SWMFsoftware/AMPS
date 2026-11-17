#!/bin/bash
# ==== version-controlled copy: header added 2026-11-15, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/tests/localnn_20261114/stage2r/chain_summary.sh (on Voyager: /data001/highresolu_vt/cops_output/tests/localnn_20261114/stage2r/chain_summary.sh)
# SHA256 of the run copy = SHA256 of this file without these header lines: 1c0dc3c6214aa3666885c82465fe9b3ff28a29ed58facee59f5c0706980b1e88
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: per-cycle summary job (localnn_chain_summary.py of localnn_s2r); CHAINS = 'LABEL=ROOT ...'.
# RUN 2026-11-15: 135078; the stage-2r comparison.
# ==== end of header ====
# per-cycle summary of option C stage-2 chains (localnn_chain_summary.py of localnn_s2r); CHAINS = "LABEL=ROOT ..."
set -uo pipefail
export NPY_DISABLE_CPU_FEATURES="AVX512F AVX512CD AVX512_KNL AVX512_KNM AVX512_SKX AVX512_CLX AVX512_CNL AVX512_ICL"
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 PYTHONUNBUFFERED=1
W=/data001/highresolu_vt/cops_output/tests/localnn_20261114/stage2r
APP=/data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda; cd $APP
M=/data001/highresolu_vt/cops_output/model_ext_frozen/localnn_s2r_20261116T044801Z
R=/data001/highresolu_vt/cops_output/runs/shap8_012m_lategrid_dz1.5
A=(); for c in ${CHAINS:?CHAINS}; do A+=(--chain "$c"); done
.venv/bin/python -W ignore $M/localnn_chain_summary.py --base $R/base --geometry $R/maps_geometry/job_124775/geometry \
  --contrib-code /data001/highresolu_vt/cops_output/model_ext_frozen/nonneg_contrib_20261112T210012Z "${A[@]}" \
  --out $W/chain_summary_$SLURM_JOB_ID.json --txt $W/chain_summary_$SLURM_JOB_ID.txt
