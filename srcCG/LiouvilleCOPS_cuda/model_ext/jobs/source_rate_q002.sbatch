#!/bin/bash
# ==== version-controlled copy: header added 2026-11-18, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/jobs/source_rate_q002.sbatch (on Voyager: /data001/highresolu_vt/cops_output/jobs/source_rate_q002.sbatch)
# SHA256 of the run copy = SHA256 of this file without these header lines: cd4a7836b5b9ec79a0fbf76cbf0eb3851525e98e96a3ce4eb87dde2a44961ce7
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: Slurm GPU job: the source-rate figure set (pub_maps_v2 activity + H2O rate at both map times) and the
#   inspection figures of one new-law run (2026-11-18)
# RUN: via source_rate_q002_submit.py, or sbatch --nodelist=<gpu node> --export=ALL with KEY BASE FT FIT MAPS GEOM MX
#   OLD_* MODEL_NAME MODEL_TEXT OLD_NAME STEPS
# ==== end of header ====
# Source-rate figure set of ONE completed new-law run (H2O q_min = 0.02 q_max, the 2026-11-18 night-side fix) for
# inspection (user 2026-11-18: "prepare the surface source rate maps for completed runs for inspection of the source
# rate"): pub_maps_v2_figures.py activity + H2O production rate at 10:33:10 and 10:39:26 UTC (templates, lon/lat and 3-D
# ground tracks, panels; LOWCUT 0.02, relief 0.7 for the rate, 0 for the activity; seen-only colouring) and, per time,
# source_rate_q002_inspect.py (new-law F | old-law F of the counterpart, common classes | new q/q_max). GPU (the
# rasterizer is a CUDA kernel). Shared node (lead 2026-11-18: gpu003, one job at a time, --mem <= 24G, abort if less
# than 30 GB host memory is available, heavy processes with oom_score_adj 1000).
# Environment (sbatch --export=ALL; set MODEL_NAME / MODEL_TEXT / OLD_NAME in the submitting shell, they contain commas):
#   KEY      output folder name under $COPS_OUTPUT_ROOT/figures_cluster/12M/source_rate_q002 (or OUT)
#   BASE FT FIT MAPS GEOM MX       the new-law state (MX: folder of the maps' frozen nonneg_maps.py)
#   OLD_BASE OLD_FT OLD_FIT OLD_MAPS OLD_GEOM OLD_MX   its old-law counterpart (for the inspection figure)
#   MODEL_NAME MODEL_TEXT OLD_NAME  caption texts; STEPS (default "v2 inspect"); TBASE TBUNDLE (track, default as before)
#SBATCH --job-name=cops_srcrate_q002
#SBATCH --partition=gpu
#SBATCH --account=general
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gpus-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=24G
#SBATCH --time=01:00:00
#SBATCH --output=/data001/highresolu_vt/cops_output/jobs/logs/cops_srcrate_q002-%j.out
set -uo pipefail
export AMPS=/data001/highresolu_vt/AMPS APP=/data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda
export T=/data001/highresolu_vt/cops_cluster_transfer_20261026 COPS_OUTPUT_ROOT=/data001/highresolu_vt/cops_output
export COPS_DATA_ROOT=$T/rosetta_data_end_of_mission_cops
export CUDA_PATH=/usr/local/cuda-12.9 CUDA_HOME=/usr/local/cuda-12.9 CUPY_CACHE_DIR=$COPS_OUTPUT_ROOT/.cupy
export PATH=$CUDA_HOME/bin:$PATH LD_LIBRARY_PATH=$CUDA_HOME/lib64${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1 MPLBACKEND=Agg
S=$COPS_OUTPUT_ROOT/figures_cluster/scripts; PY="$APP/.venv/bin/python -W ignore"
KEY=${KEY:?}; O=${OUT:-$COPS_OUTPUT_ROOT/figures_cluster/12M/source_rate_q002/$KEY}
BASE=${BASE:?}; FT=${FT:?}; FIT=${FIT:?}; MAPS=${MAPS:?}; GEOM=${GEOM:?}; MX=${MX:?}
OLD_BASE=${OLD_BASE:-}; OLD_FT=${OLD_FT:-}; OLD_FIT=${OLD_FIT:-}; OLD_MAPS=${OLD_MAPS:-}; OLD_GEOM=${OLD_GEOM:-}; OLD_MX=${OLD_MX:-}
MODEL_NAME=${MODEL_NAME:?}; MODEL_TEXT=${MODEL_TEXT:-}; OLD_NAME=${OLD_NAME:-}; STEPS=${STEPS:-v2 inspect}
TBASE=${TBASE:-$COPS_OUTPUT_ROOT/runs/shap8_012m_flat_earlier/base}
TBUNDLE=${TBUNDLE:-$COPS_OUTPUT_ROOT/runs/shap8_012m_flat_earlier/bundle/restart_templates.npz}
stamp() { echo "$(date -u +%FT%TZ) $*"; }
heavy() { bash -c 'echo 1000 > /proc/self/oom_score_adj; exec "$@"' heavy "$@"; }
avail=$(awk '/MemAvailable/ {print int($2/1048576)}' /proc/meminfo)
stamp "job $SLURM_JOB_ID on $(hostname): host MemAvailable ${avail} GB"
[ "$avail" -ge 30 ] || { stamp "ABORT: less than 30 GB host memory available (lead's condition)"; exit 3; }
for f in "$BASE/COPS_model.dat" "$BASE/elements.bin" "$FT/finetune_state.npz" "$FIT/fit_restart.npz" "$FIT/nonneg_facets.npz" \
         "$MAPS/maps_manifest.json" "$GEOM/times.txt" "$MX/nonneg_maps.py" "$TBASE/elements.bin" "$TBUNDLE"; do
    [ -e "$f" ] || { stamp "missing: $f"; exit 2; }; done
mkdir -p "$O" || exit 1
P=$O/provenance-$SLURM_JOB_ID; mkdir "$P" || exit 1
cd "$APP" || exit 1
for v in KEY O BASE FT FIT MAPS GEOM MX OLD_BASE OLD_FT OLD_FIT OLD_MAPS OLD_GEOM OLD_MX TBASE TBUNDLE STEPS MODEL_NAME MODEL_TEXT OLD_NAME; do
    printf '%s=%s\n' "$v" "${!v}"; done > "$P/inputs.txt"
cp -p "$0" "$P/source_rate_q002.sbatch.used"
mkdir "$P/code" && cp -p "$S/pub_maps_v2_figures.py" "$S/pub_maps_figures.py" "$S/illumination_maps_figures.py" \
    "$S/source_rate_q002_inspect.py" "$APP/cops_cuda/figures.py" "$MX/nonneg_maps.py" "$P/code/" || exit 1
[ -n "$OLD_MX" ] && cp -p "$OLD_MX/nonneg_maps.py" "$P/code/nonneg_maps_old.py"
(cd "$P/code" && sha256sum * > SHA256SUMS)
{ git -C "$AMPS" rev-parse HEAD; git -C "$AMPS" status --short; } > "$P/git.txt" 2>&1
{ type module >/dev/null 2>&1 && module list; } > "$P/module_list.txt" 2>&1 || echo "module command not available (no modules loaded; CUDA from $CUDA_HOME)" >> "$P/module_list.txt"
nvidia-smi -q > "$P/nvidia_smi_q.txt" 2>&1
{ nvcc --version; gcc --version | head -1; grep -m1 CSPICE /home/vtenishev/SPICE/cspice/doc/version.txt; lscpu | grep -E 'Model name|^CPU\(s\)'; } > "$P/toolchain.txt" 2>&1
~/.local/bin/uv pip list --python "$APP/.venv/bin/python" > "$P/pip_list.txt" 2>&1
env | grep -E '^(SLURM_JOB_ID|SLURM_JOB_NODELIST|SLURM_CPUS_PER_TASK|SLURM_MEM_PER_NODE|CUDA_[A-Z_]*|CUPY_[A-Z_]*|OMP_NUM_THREADS|MPLBACKEND)=' | sort > "$P/env.txt"
{ echo '#!/bin/bash'; grep '^export ' "$0"; echo "cd $APP"; } > "$P/commands.sh"
run() { printf '%q ' "$@" >> "$P/commands.sh"; echo >> "$P/commands.sh"; "$@"; }
free -g > "$P/free_start.txt" 2>&1; df -h "$COPS_OUTPUT_ROOT" > "$P/df_start.txt" 2>&1
MT=(); [ -n "$MODEL_TEXT" ] && MT=(--model-text "$MODEL_TEXT")
COMMON=(--base "$BASE" --track-base "$TBASE" --bundle "$TBUNDLE" --ft-dir "$FT" --fit-dir "$FIT" --maps-dir "$MAPS"
        --geometry "$GEOM" --maps-code "$MX/nonneg_maps.py" --out-dir "$O" --classes 15,20 --low-cut-area 0.02
        --model-name "$MODEL_NAME" "${MT[@]}" --device 0)
rc=0
if [[ " $STEPS " == *" v2 "* ]]; then
  run heavy $PY "$S/pub_maps_v2_figures.py" "${COMMON[@]}" --quantity activity --template-relief 0 > "$P/v2_activity.log" 2>&1 || rc=$?
  stamp "v2 activity rc=$rc"; tail -n 3 "$P/v2_activity.log"
  for t in 10:33:10 10:39:26; do
    run heavy $PY "$S/pub_maps_v2_figures.py" "${COMMON[@]}" --quantity h2o --time $t --template-relief 0.7 > "$P/v2_h2o_${t//:/}.log" 2>&1 || rc=$?
    stamp "v2 h2o $t rc=$rc"; tail -n 3 "$P/v2_h2o_${t//:/}.log"
  done
fi
if [[ " $STEPS " == *" inspect "* ]]; then
  for t in 10:33:10 10:39:26; do
    run heavy $PY "$S/source_rate_q002_inspect.py" --time $t --out-dir "$O" \
      --new-base "$BASE" --new-ft "$FT" --new-fit "$FIT" --new-maps "$MAPS" --new-geometry "$GEOM" --new-maps-code "$MX/nonneg_maps.py" \
      --old-base "$OLD_BASE" --old-ft "$OLD_FT" --old-fit "$OLD_FIT" --old-maps "$OLD_MAPS" --old-geometry "$OLD_GEOM" \
      --old-maps-code "$OLD_MX/nonneg_maps.py" --track-base "$TBASE" --bundle "$TBUNDLE" \
      --new-name "$MODEL_NAME" --old-name "$OLD_NAME" --low-cut-area 0.02 --template-relief 0.7 --device 0 > "$P/inspect_${t//:/}.log" 2>&1 || rc=$?
    stamp "inspect $t rc=$rc"; tail -n 3 "$P/inspect_${t//:/}.log"
  done
fi
free -g > "$P/free_end.txt" 2>&1; df -h "$COPS_OUTPUT_ROOT" > "$P/df_end.txt" 2>&1
stamp "done: rc $rc"
cp -p "$COPS_OUTPUT_ROOT/jobs/logs/cops_srcrate_q002-$SLURM_JOB_ID.out" "$P/job.log"
exit $rc
