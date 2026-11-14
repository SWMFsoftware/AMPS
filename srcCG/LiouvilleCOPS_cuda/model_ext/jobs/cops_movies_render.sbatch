#!/bin/bash
# ==== version-controlled copy: header added 2026-11-13, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/jobs/cops_movies_render.sbatch (on Voyager: /data001/highresolu_vt/cops_output/jobs/cops_movies_render.sbatch)
# SHA256 of the run copy = SHA256 of this file without these header lines: bdfeea022d7d3d4b487675e0ebe373f7d190286a2f8fc429cd4a6cc315971c86
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: COPS movies step 2 (render on all GPUs); 2026-11-14: SC_BLOCK_FROM, COLOR_CONTRIB / COLOR_RANGE (defaults
#   off).
# RUN: sbatch --export=ALL (see the header)
# ==== end of header ====
# COPS movies, step 2: render the clean view frames (RG tint and NG tint) on all GPUs of the allocation.
# figures_cluster/scripts/cops_movies.py render: one worker per GPU; every worker claims frames dynamically (exclusive
# creation of $WORK/frames/claims/fNNNNN), so several jobs on several nodes share the 1800 frames and the result does not
# depend on the split (each frame is computed independently from prep.npz and plan.npz).
# Environment: WORK (default $CHAIN/movies_work), NGPU (GPUs in this job, default SLURM_GPUS_ON_NODE or 1), WORKER0 (first
#   global worker id), NWORKERS (all workers of all jobs; only orders the frames, claims do the split), RUN, CHAIN, GEOM,
#   FRAMES_DIR (output folder, default $WORK/frames; e.g. a split-check re-render), ONLY (comma list of frames, default all),
#   MASK_RG, MASK_NG (observable part: none = the first movies (default), fit, hemisphere, ik_fov, cone:<deg>; 2026-11-13),
#   TINT_RG, TINT_NG (weight | uniform), NO_BOUNDARY=1 (no mask boundary line), SC_SHADOW=1 (spacecraft shadow, needs
#   $WORK/sc.npz; movies only) with SC_SAMPLES (points on the solar disc, default 64); MOVIES_PY (2026-11-13: the
#   cops_movies.py to run, default figures_cluster/scripts/cops_movies.py). A fixed-camera plan.npz (cops_movies.py plan
#   --camera fixed) makes the render the fixed-camera render (masks ik_fov required; 2026-11-13).
# Shared-node rule: every worker with oom_score_adj 1000 (CLAUDE.md rule 7). Submit e.g.
#   NGPU=8 WORKER0=0 NWORKERS=12 sbatch --nodelist=gpu002 --gpus-per-node=8 --mem=96G --export=ALL cops_movies_render.sbatch
#   NGPU=4 WORKER0=8 NWORKERS=12 sbatch --nodelist=gpu003 --gpus-per-node=4 --mem=26G --export=ALL cops_movies_render.sbatch
#SBATCH --job-name=cops_movies_render
#SBATCH --partition=gpu
#SBATCH --account=general
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gpus-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=24G
#SBATCH --time=02:00:00
#SBATCH --output=/data001/highresolu_vt/cops_output/jobs/logs/cops_movies_render-%j.out
set -uo pipefail
export AMPS=/data001/highresolu_vt/AMPS APP=/data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda
export COPS_OUTPUT_ROOT=/data001/highresolu_vt/cops_output
export CUDA_PATH=/usr/local/cuda-12.9 CUDA_HOME=/usr/local/cuda-12.9 CUPY_CACHE_DIR=$COPS_OUTPUT_ROOT/.cupy
export PATH=$CUDA_HOME/bin:$PATH LD_LIBRARY_PATH=$CUDA_HOME/lib64${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1 MPLBACKEND=Agg
R=${RUN:-$COPS_OUTPUT_ROOT/runs/shap8_012m_lategrid_dz1.5}; RUN=$R
CH=${CHAIN:-$R/nonneg/softplus/chain_5cycles_20261112T012553Z}
GEOM=${GEOM:-$R/maps_geometry/job_124324/geometry}
WORK=${WORK:-$CH/movies_work}
NGPU=${NGPU:-${SLURM_GPUS_ON_NODE:-1}}; WORKER0=${WORKER0:-0}; NWORKERS=${NWORKERS:-$NGPU}
FRAMES_DIR=${FRAMES_DIR:-$WORK/frames}; ONLY=${ONLY:-}
MASK_RG=${MASK_RG:-none}; MASK_NG=${MASK_NG:-none}; TINT_RG=${TINT_RG:-weight}; TINT_NG=${TINT_NG:-weight}; NO_BOUNDARY=${NO_BOUNDARY:-}
MOPTS=(); [ "$MASK_RG$MASK_NG" != nonenone ] && MOPTS=(--mask-rg "$MASK_RG" --mask-ng "$MASK_NG" --tint-rg "$TINT_RG" --tint-ng "$TINT_NG")
[ -n "$NO_BOUNDARY" ] && MOPTS+=(--no-boundary)
SC_SHADOW=${SC_SHADOW:-}; SC_SAMPLES=${SC_SAMPLES:-64}
[ -n "$SC_SHADOW" ] && MOPTS+=(--sc-shadow --sc-samples "$SC_SAMPLES")
SC_BLOCK_FROM=${SC_BLOCK_FROM:-}      # 2026-11-14: render --sc-block-from (default off)
[ -n "$SC_BLOCK_FROM" ] && MOPTS+=(--sc-block-from "$SC_BLOCK_FROM")
COLOR_CONTRIB=${COLOR_CONTRIB:-}; COLOR_RANGE=${COLOR_RANGE:-}    # 2026-11-14: render --color-contrib / --color-range (default off)
[ -n "$COLOR_CONTRIB" ] && MOPTS+=(--color-contrib "$COLOR_CONTRIB" --color-range "$COLOR_RANGE")
S=$COPS_OUTPUT_ROOT/figures_cluster/scripts
MPY=${MOVIES_PY:-$S/cops_movies.py}
PY="$APP/.venv/bin/python -W ignore"
stamp() { echo "$(date -u +%FT%TZ) $*"; }
for f in "$WORK/prep.npz" "$WORK/plan.npz" "$MPY"; do [ -e "$f" ] || { stamp "missing $f"; exit 2; }; done
P=$WORK/provenance-$SLURM_JOB_ID; mkdir "$P" || exit 1
cd "$APP" || exit 1
stamp "job $SLURM_JOB_ID on $(hostname): NGPU=$NGPU WORKER0=$WORKER0 NWORKERS=$NWORKERS CUDA_VISIBLE_DEVICES=${CUDA_VISIBLE_DEVICES:-}"
for v in RUN CH GEOM WORK NGPU WORKER0 NWORKERS FRAMES_DIR ONLY MASK_RG MASK_NG TINT_RG TINT_NG NO_BOUNDARY SC_SHADOW SC_SAMPLES SC_BLOCK_FROM COLOR_CONTRIB COLOR_RANGE MPY; do printf '%s=%s\n' "$v" "${!v}"; done > "$P/inputs.txt"
cp -p "$0" "$P/cops_movies_render.sbatch.used"
mkdir "$P/code" && cp -p "$MPY" "$APP/cops_cuda/figures.py" "$APP/cops_cuda/kernels/render.cu" "$P/code/" && (cd "$P/code" && sha256sum * > SHA256SUMS)
{ git -C "$AMPS" rev-parse HEAD; git -C "$AMPS" status --short; } > "$P/git.txt" 2>&1
{ type module >/dev/null 2>&1 && module list; } > "$P/module_list.txt" 2>&1 || echo "module command not available (no modules loaded; CUDA from $CUDA_HOME)" >> "$P/module_list.txt"
nvidia-smi -q > "$P/nvidia_smi_q.txt" 2>&1
{ nvcc --version; gcc --version | head -1; lscpu | grep -E 'Model name|^CPU\(s\)'; } > "$P/toolchain.txt" 2>&1
~/.local/bin/uv pip list --python "$APP/.venv/bin/python" > "$P/pip_list.txt" 2>&1
free -g > "$P/free_start.txt" 2>&1
{ echo '#!/bin/bash'; grep '^export ' "$0"; echo "cd $APP"; } > "$P/commands.sh"
pids=()
for i in $(seq 0 $((NGPU - 1))); do
    w=$((WORKER0 + i))
    cmd=($PY "$MPY" render --work "$WORK" --geometry "$GEOM" --device "$i" --worker "$w" --nworkers "$NWORKERS" --frames-dir "$FRAMES_DIR" ${ONLY:+--only "$ONLY"} "${MOPTS[@]}")
    printf '%q ' "${cmd[@]}" >> "$P/commands.sh"; echo >> "$P/commands.sh"
    bash -c 'echo 1000 > /proc/self/oom_score_adj; exec "$@"' heavy "${cmd[@]}" > "$P/worker_$w.log" 2>&1 &
    pids+=($!)
done
nfail=0
for pid in "${pids[@]}"; do wait $pid || nfail=$((nfail+1)); done
for i in $(seq 0 $((NGPU - 1))); do tail -n 1 "$P/worker_$((WORKER0 + i)).log"; done
free -g > "$P/free_end.txt" 2>&1
stamp "done: $nfail failed workers; frames with metadata in $FRAMES_DIR: $(ls $FRAMES_DIR/meta 2>/dev/null | wc -l)"
cp -p "$COPS_OUTPUT_ROOT/jobs/logs/cops_movies_render-$SLURM_JOB_ID.out" "$P/job.log"
[ $nfail -eq 0 ]
