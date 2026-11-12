#!/bin/bash
# ==== version-controlled copy: header added 2026-11-11, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/jobs/slope_ft_diag.sbatch (on Voyager: /data001/highresolu_vt/cops_output/jobs/slope_ft_diag.sbatch)
# SHA256 of the run copy = SHA256 of this file without these header lines: 7c0e56c1bf7a37487c15d6943e0da57c3feb7c19db9dcccf23f1769aafc40e7d
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: job (diagnostic, read-only). Runs checks/slopemul_ftdiag_20261110/ft_diag.py on the two fits whose fine tunings
#   failed checks 6 / 7 (jobs 124013, 124015), one bind; it confirmed G and showed that the failures were rounding of
#   single-facet steps (2026-11-10).
# RUN: sbatch --nodelist=gpu002 slope_ft_diag.sbatch   (job 124034 -> checks/slopemul_ftdiag_20261110/run_124034)
# ==== end of header ====
# Diagnostic GPU job (2026-11-10, after jobs 124013 and 124015 failed the fine tuning's validation checks 6 and 7):
# checks/slopemul_ftdiag_20261110/ft_diag.py on the two fits whose fine tunings failed (the short l32 test fit of job
# 123984 and the production c1 fit 124008), one bind for both. Read-only on every product; it writes only DIAG_OUT.
# See ft_diag.py's docstring for parts A-E. About 15 min.
# Settings: MODEL_EXT (default the v7 snapshot of the failed runs), DIAG_OUT (default checks/slopemul_ftdiag_20261110/
#   run_<jobid>), FITS (space-separated fit segment dirs; default the two fits).
#SBATCH --job-name=cops12m_ftdiag
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gpus-per-node=8
#SBATCH --cpus-per-task=64
#SBATCH --mem=0
#SBATCH --exclusive
#SBATCH --time=00:40:00
#SBATCH --partition=gpu
#SBATCH --account=general
#SBATCH --exclude=gpu003
#SBATCH --output=/data001/highresolu_vt/cops_output/jobs/logs/cops12m_ftdiag-%j.out
set -uo pipefail
export AMPS=/data001/highresolu_vt/AMPS APP=/data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda
export T=/data001/highresolu_vt/cops_cluster_transfer_20261026 COPS_OUTPUT_ROOT=/data001/highresolu_vt/cops_output
export COPS_DATA_ROOT=$T/rosetta_data_end_of_mission_cops
export CUDA_PATH=/usr/local/cuda-12.9 CUDA_HOME=/usr/local/cuda-12.9 CUPY_CACHE_DIR=$COPS_OUTPUT_ROOT/.cupy
export PATH=$CUDA_HOME/bin:/usr/mpi/gcc/openmpi-4.1.7rc1/bin:$PATH
export LD_LIBRARY_PATH=$CUDA_HOME/lib64${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1
UNSET=$(env | grep -o -E '^(COPS_POS_[A-Z_]*|COPS_SLOPEMUL_[A-Z_]*)' | sort | paste -sd' ' -)
[ -n "$UNSET" ] && unset $UNSET
R=$COPS_OUTPUT_ROOT/runs/shap8_012m_lategrid_dz1.5; C=$COPS_OUTPUT_ROOT/checks/slopemul_ftdiag_20261110
export MODEL_EXT=$(readlink -f "${MODEL_EXT:-$COPS_OUTPUT_ROOT/model_ext_frozen/slopemul_20261110T211229Z}")
FITS=${FITS:-"$R/slope_mult_tests/g_bumps_short_l32_20261110T201212Z/seg01 $R/slope_bumps/chain_slopemul_20261110T201212Z/c1/fit/final"}
OUTD=$(readlink -m "${DIAG_OUT:-$C/run_$SLURM_JOB_ID}"); [ -e "$OUTD" ] && { echo "exists: $OUTD" >&2; exit 2; }
mkdir -p "$OUTD"; P=$OUTD/provenance-$SLURM_JOB_ID; mkdir -p "$P"; cd "$APP" || exit 1
cp -p "$0" "$P/slope_ft_diag.sbatch"; cp -p "$C/ft_diag.py" "$P/"; cp -p "$MODEL_EXT"/*.py "$MODEL_EXT/SHA256SUMS" "$P/"
(cd "$P" && sha256sum *.py) > "$P/run_copy.sha256" 2>&1
{ git -C "$AMPS" rev-parse HEAD; git -C "$AMPS" status --short; } > "$P/git.txt" 2>&1
module list > "$P/module_list.txt" 2>&1 || true
nvidia-smi -q > "$P/nvidia_smi_q.txt" 2>&1; nvidia-smi -L > "$P/nvidia_smi_L.txt" 2>&1
{ nvcc --version; gcc --version | head -1; } > "$P/toolchain.txt" 2>&1
~/.local/bin/uv pip list --python "$APP/.venv/bin/python" > "$P/pip_list.txt" 2>&1
env | grep -E '^(CUDA_[A-Z_]*|CUPY_[A-Z_]*|COPS_[A-Z_]*|MODEL_EXT|OMP_[A-Z_]*|OPENBLAS_[A-Z_]*|MKL_[A-Z_]*|PATH|LD_LIBRARY_PATH|SLURM_JOB_NODELIST)=' | sort > "$P/environment.txt"
NG=$(nvidia-smi -L | wc -l); DEV=$(seq -s, 0 $((NG - 1)))
CMD=("$APP/.venv/bin/python" -W ignore "$P/ft_diag.py" --model-ext "$MODEL_EXT" --devices "$DEV" --l-threads 32
     --out "$OUTD/ft_diag.json")
for f in $FITS; do CMD+=(--fit "$(readlink -f "$f")"); done
{ echo '#!/bin/bash'; echo "# job $SLURM_JOB_ID on $(hostname), $(date -u +%FT%TZ); AMPS $(git -C "$AMPS" rev-parse HEAD)"
  grep '^export ' "$0"; echo "cd $APP"
  printf '%q ' bash -c 'echo 1000 > /proc/self/oom_score_adj; exec "$@"' _ "${CMD[@]}"; echo; } > "$P/commands.sh"
echo "$(date -u +%FT%TZ) job $SLURM_JOB_ID on $(hostname): ft_diag on $FITS; MODEL_EXT $MODEL_EXT; out $OUTD"
bash -c 'echo 1000 > /proc/self/oom_score_adj; exec "$@"' _ "${CMD[@]}" > "$P/ft_diag.log" 2>&1; rc=$?
grep -E '^ft_diag|Traceback|Error' "$P/ft_diag.log" | cut -c1-700
echo "$(date -u +%FT%TZ) exit $rc"
exit $rc
