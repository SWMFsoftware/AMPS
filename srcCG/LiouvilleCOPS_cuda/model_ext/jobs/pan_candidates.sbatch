#!/bin/bash
# ==== version-controlled copy: header added 2026-11-13, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/jobs/pan_candidates.sbatch (on Voyager: /data001/highresolu_vt/cops_output/jobs/pan_candidates.sbatch)
# SHA256 of the run copy = SHA256 of this file without these header lines: 997415eaddb3101204dc1c8e32e3538dc2b86c50fa5678a1f1ad847972ba7899
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: Job of pan_candidates.py with provenance.
# RUN: sbatch --nodelist=gpu002 --export=ALL (WORK, MOVIES_PY, SRC, REGIONS, OUTP, AZ, DH, HC, EXTRA)
# ==== end of header ====
# Candidate final camera poses of the panorama movie (figures_cluster/scripts/pan_candidates.py; 2026-11-14, user ~20:05 CST).
# One GPU. Environment: WORK (a detach work folder with prep.npz), MOVIES_PY (cops_movies.py to import), SRC
# (ng_source_map npz), REGIONS (ng_source_regions prefix), OUTP (output prefix), AZ/DH/HC grids (optional).
#SBATCH --job-name=cops_pan_candidates
#SBATCH --partition=gpu
#SBATCH --account=general
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gpus-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=80G
#SBATCH --time=01:00:00
#SBATCH --output=/data001/highresolu_vt/cops_output/jobs/logs/cops_pan_candidates-%j.out
set -uo pipefail
AZ=${AZ:-0:360:15}; DH=${DH:-30,60,120}; HC=${HC:-30,60,100,150,220}; EXTRA=${EXTRA:-}
export AMPS=/data001/highresolu_vt/AMPS APP=/data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda
export COPS_OUTPUT_ROOT=/data001/highresolu_vt/cops_output
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1 MPLBACKEND=Agg
export CUDA_PATH=/usr/local/cuda-12.9 CUDA_HOME=/usr/local/cuda-12.9 CUPY_CACHE_DIR=$COPS_OUTPUT_ROOT/.cupy
export PATH=$CUDA_HOME/bin:$PATH LD_LIBRARY_PATH=$CUDA_HOME/lib64${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}
G=${GEOM:-$COPS_OUTPUT_ROOT/runs/shap8_012m_lategrid_dz1.5/maps_geometry/job_124775/geometry}
S=$COPS_OUTPUT_ROOT/figures_cluster/scripts
PY="$APP/.venv/bin/python -W ignore"
stamp() { echo "$(date -u +%FT%TZ) $*"; }
P=$(dirname "$OUTP")/provenance-$SLURM_JOB_ID; mkdir -p "$P" || exit 1
cd "$APP" || exit 1
stamp "job $SLURM_JOB_ID on $(hostname), CUDA_VISIBLE_DEVICES=${CUDA_VISIBLE_DEVICES:-}"
for v in WORK MOVIES_PY SRC REGIONS OUTP G AZ DH HC EXTRA; do printf '%s=%s\n' "$v" "${!v}"; done > "$P/inputs.txt"
cp -p "$0" "$P/pan_candidates.sbatch.used"
mkdir "$P/code" && cp -p "$S/pan_candidates.py" "$MOVIES_PY" "$P/code/" && (cd "$P/code" && sha256sum * > SHA256SUMS)
{ git -C "$AMPS" rev-parse HEAD; git -C "$AMPS" status --short; } > "$P/git.txt" 2>&1
{ type module >/dev/null 2>&1 && module list; } > "$P/module_list.txt" 2>&1 || echo "module command not available (CUDA from $CUDA_HOME)" >> "$P/module_list.txt"
nvidia-smi -q > "$P/nvidia_smi_q.txt" 2>&1
{ nvcc --version; gcc --version | head -1; } > "$P/toolchain.txt" 2>&1
~/.local/bin/uv pip list --python "$APP/.venv/bin/python" > "$P/pip_list.txt" 2>&1
cmd=($PY "$S/pan_candidates.py" --movies-py "$MOVIES_PY" --work "$WORK" --geometry "$G" --src "$SRC" --regions "$REGIONS"
     --out "$OUTP" --device 0 --az "${AZ:-0:360:15}" --dh "${DH:-30,60,120}" --hc "${HC:-30,60,100,150,220}" ${EXTRA:-})
{ echo '#!/bin/bash'; grep '^export ' "$0"; echo "cd $APP"; printf '%q ' "${cmd[@]}"; echo; } > "$P/commands.sh"
bash -c 'echo 1000 > /proc/self/oom_score_adj; exec "$@"' heavy "${cmd[@]}" > "$P/pan_candidates.log" 2>&1; rc=$?
tail -n 8 "$P/pan_candidates.log"
stamp "done: rc $rc"
cp -p "$COPS_OUTPUT_ROOT/jobs/logs/cops_pan_candidates-$SLURM_JOB_ID.out" "$P/job.log"
exit $rc
