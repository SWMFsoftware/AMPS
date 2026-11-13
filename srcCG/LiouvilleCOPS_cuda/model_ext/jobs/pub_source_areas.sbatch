#!/bin/bash
# ==== version-controlled copy: header added 2026-11-13, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/jobs/pub_source_areas.sbatch (on Voyager: /data001/highresolu_vt/cops_output/jobs/pub_source_areas.sbatch)
# SHA256 of the run copy = SHA256 of this file without these header lines: 2c8d677ba5ed770d8f4ecd7a7680d033a60b3b2361f5b5319d8bfe52475f6bee
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: job script (GPU): pub_source_area_figures.py and cliff_importance_figures.py; MODEL_NAME, FIT, MAPS,
#   CLIFF_GROUPS added 2026-11-13.
# RUN: STEPS=sources | STEPS=cliffs V2=... sbatch --export=ALL pub_source_areas.sbatch
# ==== end of header ====
# P% source-area maps (union and occurrence, NG and RG, P = 50..90) and the cliff-importance figures of the 12M softplus
# model after 5 cycles (user request 2026-11-12, second round; data: model_ext nonneg_contrib.py, job 124549).
#   1. figures_cluster/scripts/pub_source_area_figures.py (GPU: rasterizer)  -> OUT/source_areas/, OUT/checks/
#   2. figures_cluster/scripts/cliff_importance_figures.py (CPU)            -> OUT/cliffs/
# Every input can be set in the submitting environment (sbatch --export=ALL; values without commas there); defaults =
# the softplus 5-cycle final state of 2026-11-12:
#   RUN (fit base = $RUN/base), CHAIN, CONTRIB (default $CHAIN/contributions), TBASE / TBUNDLE (observation-only base
#   and its bundle: track, impact, late window), MX (frozen nonneg_maps code: elements reader), OUT, V2 (folder with
#   checks/slope_bins_total_*.json for the cross-check; default OUT), STEPS (default sources,cliffs; set in the
#   environment, not in --export), MODEL_NAME (name of the fitted state in the captions; 2026-11-13), NCYC / FIT (the
#   fit of the state: multiplier bounds in the cliff captions; default $CHAIN/c5/fit/final), MAPS (maps of the state:
#   Q_CO2 / Q_H2O; default $CHAIN/maps), CLIFF_GROUPS (slope-bin groups of the cliff descent figures: default | all).
# Provenance in OUT/provenance-<job>/ (inputs.txt lists every input). Shared node: oom_score_adj 1000 (CLAUDE.md rule 7).
#SBATCH --job-name=cops_pub_source_areas
#SBATCH --partition=gpu
#SBATCH --account=general
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gpus-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=01:00:00
#SBATCH --output=/data001/highresolu_vt/cops_output/jobs/logs/cops_pub_source_areas-%j.out
set -uo pipefail
export AMPS=/data001/highresolu_vt/AMPS APP=/data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda
export T=/data001/highresolu_vt/cops_cluster_transfer_20261026 COPS_OUTPUT_ROOT=/data001/highresolu_vt/cops_output
export COPS_DATA_ROOT=$T/rosetta_data_end_of_mission_cops
export CUDA_PATH=/usr/local/cuda-12.9 CUDA_HOME=/usr/local/cuda-12.9 CUPY_CACHE_DIR=$COPS_OUTPUT_ROOT/.cupy
export PATH=$CUDA_HOME/bin:$PATH LD_LIBRARY_PATH=$CUDA_HOME/lib64${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1 MPLBACKEND=Agg
O=${OUT:-$COPS_OUTPUT_ROOT/figures_cluster/12M/full_set_softplus_5cycles_v2}
R=${RUN:-$COPS_OUTPUT_ROOT/runs/shap8_012m_lategrid_dz1.5}; BASE=$R/base; RUN=$R
CH=${CHAIN:-$R/nonneg/softplus/chain_5cycles_20261112T012553Z}; CONTRIB=${CONTRIB:-$CH/contributions}
TBASE=${TBASE:-$COPS_OUTPUT_ROOT/runs/shap8_012m_flat_earlier/base}
TBUNDLE=${TBUNDLE:-$COPS_OUTPUT_ROOT/runs/shap8_012m_flat_earlier/bundle/restart_templates.npz}
MX=${MX:-$COPS_OUTPUT_ROOT/model_ext_frozen/nonneg_maps_20261112T021623Z}; V2=${V2:-$O}
STEPS=${STEPS:-sources,cliffs}
MODEL_NAME=${MODEL_NAME:-the softplus activity model, final state after 5 fit and fine-tuning cycles}
NCYC=${NCYC:-5}; FIT=${FIT:-$CH/c$NCYC/fit/final}; MAPS=${MAPS:-$CH/maps}; CLIFF_GROUPS=${CLIFF_GROUPS:-default}
S=$COPS_OUTPUT_ROOT/figures_cluster/scripts
PY="$APP/.venv/bin/python -W ignore"
stamp() { echo "$(date -u +%FT%TZ) $*"; }
heavy() { bash -c 'echo 1000 > /proc/self/oom_score_adj; exec "$@"' heavy "$@"; }
for f in "$BASE/COPS_model.dat" "$BASE/elements.bin" "$TBASE/COPS_model.dat" "$TBASE/elements.bin" "$TBUNDLE" \
         "$CONTRIB/source_areas_NG.npz" "$CONTRIB/source_areas_RG.npz" "$CONTRIB/slope_bins.npz" "$MX/nonneg_maps.py" \
         "$CONTRIB/contrib_manifest.json" "$FIT/fit_restart.npz" "$MAPS/maps_manifest.json" \
         "$S/pub_source_area_figures.py" "$S/cliff_importance_figures.py" "$S/pub_maps_v2_figures.py" "$S/pub_maps_figures.py"; do
    [ -e "$f" ] || { stamp "missing: $f"; exit 2; }; done
mkdir -p "$O" || exit 1
P=$O/provenance-$SLURM_JOB_ID; mkdir "$P" || exit 1
cd "$APP" || exit 1
stamp "job $SLURM_JOB_ID on $(hostname): STEPS=$STEPS OUT=$O"
for v in RUN BASE CH CONTRIB TBASE TBUNDLE MX O V2 STEPS MODEL_NAME NCYC FIT MAPS CLIFF_GROUPS; do printf '%s=%s\n' "$v" "${!v}"; done > "$P/inputs.txt"
cp -p "$0" "$P/pub_source_areas.sbatch.used"
mkdir "$P/code" && cp -p "$S/pub_source_area_figures.py" "$S/cliff_importance_figures.py" "$S/pub_maps_v2_figures.py" \
    "$S/pub_maps_figures.py" "$APP/cops_cuda/figures.py" "$MX/nonneg_maps.py" "$P/code/" || exit 1
(cd "$P/code" && sha256sum * > SHA256SUMS)
{ git -C "$AMPS" rev-parse HEAD; git -C "$AMPS" status --short; } > "$P/git.txt" 2>&1
{ type module >/dev/null 2>&1 && module list; } > "$P/module_list.txt" 2>&1 || echo "module command not available (no modules loaded; CUDA from $CUDA_HOME)" >> "$P/module_list.txt"
nvidia-smi -q > "$P/nvidia_smi_q.txt" 2>&1
{ nvcc --version; gcc --version | head -1; grep -m1 CSPICE /home/vtenishev/SPICE/cspice/doc/version.txt; lscpu | grep -E 'Model name|^CPU\(s\)'; } > "$P/toolchain.txt" 2>&1
~/.local/bin/uv pip list --python "$APP/.venv/bin/python" > "$P/pip_list.txt" 2>&1 || $PY -m pip list > "$P/pip_list.txt" 2>&1
env | grep -E '^(SLURM_JOB_ID|SLURM_JOB_NODELIST|SLURM_CPUS_PER_TASK|SLURM_MEM_PER_NODE|CUDA_[A-Z_]*|CUPY_[A-Z_]*|STEPS|OUT|RUN|CHAIN|CONTRIB|TBASE|TBUNDLE|MX|V2|OMP_NUM_THREADS|MPLBACKEND)=' | sort > "$P/env.txt"
{ echo '#!/bin/bash'; grep '^export ' "$0"; echo "cd $APP"; } > "$P/commands.sh"
run() { printf '%q ' "$@" >> "$P/commands.sh"; echo >> "$P/commands.sh"; "$@"; }
free -g > "$P/free_start.txt" 2>&1; df -h "$COPS_OUTPUT_ROOT" > "$P/df_start.txt" 2>&1
nfail=0
if [[ ",$STEPS," == *,sources,* ]]; then
    stamp "source areas"
    run heavy $PY "$S/pub_source_area_figures.py" --base "$BASE" --track-base "$TBASE" --bundle "$TBUNDLE" \
        --contrib-dir "$CONTRIB" --maps-code "$MX/nonneg_maps.py" --out-dir "$O" --model-name "$MODEL_NAME" \
        --device 0 > "$P/source_areas.log" 2>&1 \
        || { stamp "source areas FAILED"; nfail=$((nfail+1)); }
    tail -n 12 "$P/source_areas.log"
fi
if [[ ",$STEPS," == *,cliffs,* ]]; then
    stamp "cliff figures"
    run $PY "$S/cliff_importance_figures.py" --contrib-dir "$CONTRIB" --base "$BASE" --bundle "$TBUNDLE" \
        --v2-dir "$V2" --out-dir "$O/cliffs" --model-name "$MODEL_NAME" --fit-dir "$FIT" --maps-dir "$MAPS" \
        --groups "$CLIFF_GROUPS" > "$P/cliffs.log" 2>&1 || { stamp "cliffs FAILED"; nfail=$((nfail+1)); }
    tail -n 14 "$P/cliffs.log"
fi
free -g > "$P/free_end.txt" 2>&1; df -h "$COPS_OUTPUT_ROOT" > "$P/df_end.txt" 2>&1
stamp "done: $nfail failures"
cp -p "$COPS_OUTPUT_ROOT/jobs/logs/cops_pub_source_areas-$SLURM_JOB_ID.out" "$P/job.log"
[ $nfail -eq 0 ]
