#!/bin/bash
# ==== version-controlled copy: header added 2026-11-22, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/jobs/localnn_ac_chain.sbatch (on Voyager: /data001/highresolu_vt/cops_output/jobs/localnn_ac_chain.sbatch)
# SHA256 of the run copy = SHA256 of this file without these header lines: fd03520aea80e384bb069cdb661561c9be48a46c6ae0b5bfcac270e5a9734b37
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: five-cycle chain of the all-held candidate
# RUN: through jobs/localnn_ac_submit.sh
# ==== end of header ====
# The five-cycle chain (fit + late-facet fine tuning) of the all-held candidate model, part 'c' of the slope-test
# follow-up (user 2026-11-22). Cycle 1 = the lowest converged all-held fit (localnn_ac.py pick; the existing fit from
# D2 on a tie). The driver is jobs/localnn_chain.sbatch (frozen copy in TOOLS) with HOLD_ALL_MULT=1; settings as the
# templates-off chain of 2026-11-18. Env: RUNROOT, TOOLS; test only: C1_FIT, CYCLES, NFEV, FT_MAXNFEV.
#SBATCH --job-name=cops_ac_chain
#SBATCH --partition=normal128cores
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --output=/data001/highresolu_vt/cops_output/jobs/logs/cops_localnn_chain-%j.out
set -uo pipefail
R=/data001/highresolu_vt/cops_output
RUNROOT=${RUNROOT:?RUNROOT is required}
TOOLS=${TOOLS:?TOOLS is required}
PY=/data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda/.venv/bin/python
(cd "$TOOLS" && sha256sum -c --quiet SHA256SUMS) || { echo "frozen helper changed" >&2; exit 2; }
N=$R/runs/shap8_012m_lategrid_dz1.5
L=$N/nonneg/localnn
if [ -z "${C1_FIT:-}" ]; then
  C1_FIT=$("$PY" "$TOOLS/localnn_ac.py" pick --runroot "$RUNROOT") || { echo "no chain start (see the message above)" >&2; exit 3; }
fi
echo "chain cycle 1: $C1_FIT"
export C1_FIT HOLD_ALL_MULT=1 C1_START=state DEVICE=cpu TEMPLATES=off CHECK=1 BLAS_THREADS=16
export CHAIN=$RUNROOT/chain
export CYCLES=${CYCLES:-5} NFEV=${NFEV:-4000} FT_MAXNFEV=${FT_MAXNFEV:-500} WINDOW_MIN=5 CHECKPOINT=250
export LAMBDA=8.249495752351283e-05 RULE=$L/lambda_rule_anchored_20261114.json FM=0.01
export BASE_DIR=$N/base_q002 D2_DIR=$N/clip_refit/D2_q002_warm/seg01
export GRID_DIR=$L/grid_150m DESIGN_DIR=$L/grid_150m/design_h2oqmin0.02_d2q002
export LATE=$L/late_kernels_150m_h2oqmin0.02.npz
export FT_EXTRA="--delta-penalty-f 10 --delta-cref 4.2688066575e-03 --delta-cap 4.605170185988092"
export MODEL_EXT=$R/model_ext_frozen/localnn_slopetest_20261122T033334Z
export NONNEG_CODE=$R/model_ext_frozen/nonneg_20261118T072729Z
export CONTRIB=$R/model_ext_frozen/nonneg_contrib_20261118T072729Z
mkdir -p "$CHAIN"
cp -p "$0" "$CHAIN/localnn_ac_chain.sbatch.used-$SLURM_JOB_ID"
exec bash "$TOOLS/localnn_chain.sbatch"
