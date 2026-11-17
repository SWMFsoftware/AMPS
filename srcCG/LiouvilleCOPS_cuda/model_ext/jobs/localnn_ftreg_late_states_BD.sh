#!/bin/bash
# ==== version-controlled copy: header added 2026-11-15, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/tests/localnn_20261114/ftreg/late_states_BD.sh (on Voyager: /data001/highresolu_vt/cops_output/tests/localnn_20261114/ftreg/late_states_BD.sh)
# SHA256 of the run copy = SHA256 of this file without these header lines: 331cb24eccad8390f51582a3e6c44ba03e2e007eeded49420fa20507ec7d5609
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: late kernels + fit-state pressures of the B and D c5 fits for the CPU fine-tuning engine (localnn_late.py
#   --state; G bitwise = late_kernels_150m.npz).
# RUN 2026-11-15: 134495.
# ==== end of header ====
# fine-tuning regularization test (lead / user 2026-11-15 ~21:10 CST): fit-state pressures of every base row for the B and
# D c5 fits (agent 3's chains), for the CPU fine-tuning engine (localnn_late.py --state; frozen localnn_s2post snapshot)
set -uo pipefail
export NPY_DISABLE_CPU_FEATURES="AVX512F AVX512CD AVX512_KNL AVX512_KNM AVX512_SKX AVX512_CLX AVX512_CNL AVX512_ICL"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 PYTHONUNBUFFERED=1
APP=/data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda; cd $APP
M=/data001/highresolu_vt/cops_output/model_ext_frozen/localnn_s2post_20261115T195616Z
R=/data001/highresolu_vt/cops_output/runs/shap8_012m_lategrid_dz1.5; N=$R/nonneg
W=/data001/highresolu_vt/cops_output/tests/localnn_20261114/ftreg
P=$W/provenance-$SLURM_JOB_ID; mkdir -p $P; cp -p "$0" $P/late_states_BD.sh.used
{ git -C /data001/highresolu_vt/AMPS rev-parse HEAD; git -C /data001/highresolu_vt/AMPS status --short; } > $P/git.txt 2>&1
echo "CPU job on $(hostname): no GPU used" > $P/nvidia_smi_q.txt; { gcc --version | head -1; lscpu | grep 'Model name'; } > $P/toolchain.txt
~/.local/bin/uv pip list --python $APP/.venv/bin/python > $P/pip_list.txt 2>&1
bash -c 'echo 1000 > /proc/self/oom_score_adj; exec "$@"' _ .venv/bin/python -W ignore $M/localnn_late.py --base $R/base --d2 $R/clip_refit/D2/final \
  --point-table $R/late_grid_points.txt --contrib-code /data001/highresolu_vt/cops_output/model_ext_frozen/nonneg_contrib_20261112T210012Z \
  --nonneg-code /data001/highresolu_vt/cops_output/model_ext_frozen/nonneg_20261113T035423Z --out $W/late_kernels_BD.npz \
  --state B_c5=$N/background/stage2_bg1e-2_20261115T184949Z/c5/fit/seg01/nonneg_facets.npz:1.6805669286251197e-45 \
  --state D_c5=$N/logspread/stage2_f0.1_20261115T184949Z/c5/fit/seg01/nonneg_facets.npz:7.12261847742288e-37 --procs 64 2>&1 | tail -3
