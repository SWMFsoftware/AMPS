#!/bin/bash
# ==== version-controlled copy: header added 2026-11-13, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/jobs/cops_movies_diag.sbatch (on Voyager: /data001/highresolu_vt/cops_output/jobs/cops_movies_diag.sbatch)
# SHA256 of the run copy = SHA256 of this file without these header lines: e2a1cdd8cac979963ecd67ca69e9a7c9fd0235b2500cab7af5bb96b0e85e2c7b
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: job script (GPU): cops_movies.py diagnose (dark-surface classification with the exact shadow ray, annotated
#   stills; 2026-11-13 artifact question).
# RUN: WORK=... GEOM=... OUTD=... BOXES='frame,x0,x1,y0,y1 ...' sbatch --export=ALL cops_movies_diag.sbatch
# ==== end of header ====
# COPS movies: dark-surface diagnosis of rendered frames (figures_cluster/scripts/cops_movies.py diagnose; user question
# 2026-11-13 "starting from about 20 m altitude, the movies have some artifact on the surface ... what is it actually?").
# Environment: WORK (movie work folder), GEOM, FRAMES (default 1500,1600,1700,1799), FRAMES_DIR (default $WORK/frames),
# OUTD (output folder), BOXES (space-separated frame,x0,x1,y0,y1 regions listed facet by facet). GPU: 1.
#SBATCH --job-name=cops_movies_diag
#SBATCH --partition=gpu
#SBATCH --account=general
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gpus-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=24G
#SBATCH --time=01:00:00
#SBATCH --output=/data001/highresolu_vt/cops_output/jobs/logs/cops_movies_diag-%j.out
set -uo pipefail
export AMPS=/data001/highresolu_vt/AMPS APP=/data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda
export COPS_OUTPUT_ROOT=/data001/highresolu_vt/cops_output
export CUDA_PATH=/usr/local/cuda-12.9 CUDA_HOME=/usr/local/cuda-12.9 CUPY_CACHE_DIR=$COPS_OUTPUT_ROOT/.cupy
export PATH=$CUDA_HOME/bin:$PATH LD_LIBRARY_PATH=$CUDA_HOME/lib64${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1 MPLBACKEND=Agg
WORK=${WORK:?}; GEOM=${GEOM:?}; FRAMES=${FRAMES:-1500,1600,1700,1799}; FRAMES_DIR=${FRAMES_DIR:-$WORK/frames}
OUTD=${OUTD:?}; BOXES=${BOXES:-}
S=$COPS_OUTPUT_ROOT/figures_cluster/scripts
PY="$APP/.venv/bin/python -W ignore"
stamp() { echo "$(date -u +%FT%TZ) $*"; }
mkdir -p "$OUTD" || exit 1
P=$OUTD/provenance-$SLURM_JOB_ID; mkdir "$P" || exit 1
cd "$APP" || exit 1
stamp "job $SLURM_JOB_ID on $(hostname)"
for v in WORK GEOM FRAMES FRAMES_DIR OUTD BOXES; do printf '%s=%s\n' "$v" "${!v}"; done > "$P/inputs.txt"
cp -p "$0" "$P/cops_movies_diag.sbatch.used"
mkdir "$P/code" && cp -p "$S/cops_movies.py" "$APP/cops_cuda/figures.py" "$APP/cops_cuda/kernels/render.cu" "$P/code/" && (cd "$P/code" && sha256sum * > SHA256SUMS)
{ git -C "$AMPS" rev-parse HEAD; git -C "$AMPS" status --short; } > "$P/git.txt" 2>&1
{ type module >/dev/null 2>&1 && module list; } > "$P/module_list.txt" 2>&1 || echo "module command not available (no modules loaded; CUDA from $CUDA_HOME)" >> "$P/module_list.txt"
nvidia-smi -q > "$P/nvidia_smi_q.txt" 2>&1
{ nvcc --version; gcc --version | head -1; lscpu | grep -E 'Model name|^CPU\(s\)'; } > "$P/toolchain.txt" 2>&1
~/.local/bin/uv pip list --python "$APP/.venv/bin/python" > "$P/pip_list.txt" 2>&1
{ echo '#!/bin/bash'; grep '^export ' "$0"; echo "cd $APP"; } > "$P/commands.sh"
BA=(); for b in $BOXES; do BA+=(--box "$b"); done
cmd=($PY "$S/cops_movies.py" diagnose --work "$WORK" --geometry "$GEOM" --frames "$FRAMES" --frames-dir "$FRAMES_DIR" --out-dir "$OUTD" "${BA[@]}" --device 0)
printf '%q ' "${cmd[@]}" >> "$P/commands.sh"; echo >> "$P/commands.sh"
bash -c 'echo 1000 > /proc/self/oom_score_adj; exec "$@"' heavy "${cmd[@]}" > "$P/diagnose.log" 2>&1; rc=$?
tail -n 30 "$P/diagnose.log"
stamp "done: rc $rc"
cp -p "$COPS_OUTPUT_ROOT/jobs/logs/cops_movies_diag-$SLURM_JOB_ID.out" "$P/job.log"
exit $rc
