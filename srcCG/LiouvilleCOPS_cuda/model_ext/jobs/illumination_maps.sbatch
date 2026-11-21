#!/bin/bash
# ==== version-controlled copy: header added 2026-11-18, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/jobs/illumination_maps.sbatch (on Voyager: /data001/highresolu_vt/cops_output/jobs/illumination_maps.sbatch)
# SHA256 of the run copy = SHA256 of this file without these header lines: e1d35b65e89b03c7a2a5ea75db8725bac342014083a1e156d05081379948e253
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: Slurm job (1 GPU): illumination maps of one production-rate map time in every production-rate layout, and the
#   side-by-side figure (illumination_maps_figures.py figures; 2026-11-18)
# RUN: sbatch --dependency=afterok:<numbers job> --export=ALL,TIME=10:33:10 illumination_maps.sbatch (and TIME=10:39:26)
# ==== end of header ====
# Illumination maps at ONE production-rate map time (TIME=10:33:10 or 10:39:26 UTC, 2016-09-30) in every layout of the
# full sets' production-rate figures, plus the side-by-side figure production rate | illumination factor
# (figures_cluster/scripts/illumination_maps_figures.py figures; 2026-11-18, user: "check the production rate plots --
# i dont see the effect of the shadow or night side decrease of the productions -- is it there? generate illumination
# maps for the time intervals that corresponds to that on the production rate maps"). 1 GPU. Needs the output of
# illumination_numbers.sbatch (submit with --dependency=afterok:<numbers job>).
# Environment (defaults = the stage-2 softplus reg state): TIME (required), RUN, CHAIN, NCYC, FT, FIT, MAPS, GEOM, MX,
# TBASE, TBUNDLE, OUT, PR (its full-set folder: the production-rate classes are checked against its manifest),
# NUMBERS_DIR (default $OUT/checks), TEMPLATE_RELIEF (0.7, as the production-rate templates), GT_RELIEF (0, as the
# production-rate ground tracks), ONLY (test: quantity[:whole|seen],side), MODEL_NAME (set it in the submitting shell).
# Provenance in $OUT/provenance-<job>/. Shared node: heavy processes with oom_score_adj 1000 (CLAUDE.md rule 7).
#SBATCH --job-name=cops_illum_maps
#SBATCH --partition=gpu
#SBATCH --account=general
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gpus-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=01:00:00
#SBATCH --output=/data001/highresolu_vt/cops_output/jobs/logs/cops_illum_maps-%j.out
set -uo pipefail
export AMPS=/data001/highresolu_vt/AMPS APP=/data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda
export T=/data001/highresolu_vt/cops_cluster_transfer_20261026 COPS_OUTPUT_ROOT=/data001/highresolu_vt/cops_output
export COPS_DATA_ROOT=$T/rosetta_data_end_of_mission_cops
export CUDA_PATH=/usr/local/cuda-12.9 CUDA_HOME=/usr/local/cuda-12.9 CUPY_CACHE_DIR=$COPS_OUTPUT_ROOT/.cupy
export PATH=$CUDA_HOME/bin:$PATH LD_LIBRARY_PATH=$CUDA_HOME/lib64${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1 MPLBACKEND=Agg
TIME=${TIME:?set TIME=10:33:10 or TIME=10:39:26}
R=${RUN:-$COPS_OUTPUT_ROOT/runs/shap8_012m_lategrid_dz1.5}; BASE=$R/base
CH=${CHAIN:-$R/nonneg/softplus/stage2r_smm1_20261116T160802Z}; NCYC=${NCYC:-5}
FT=${FT:-$CH/c$NCYC/ft/final}; FIT=${FIT:-$CH/c$NCYC/fit/final}; MAPS=${MAPS:-$CH/maps_c$NCYC}
GEOM=${GEOM:-$R/maps_geometry/job_136152/geometry}; MX=${MX:-$COPS_OUTPUT_ROOT/model_ext_frozen/nonneg_maps_20261112T021623Z}
TBASE=${TBASE:-$COPS_OUTPUT_ROOT/runs/shap8_012m_flat_earlier/base}
TBUNDLE=${TBUNDLE:-$COPS_OUTPUT_ROOT/runs/shap8_012m_flat_earlier/bundle/restart_templates.npz}
O=${OUT:-$COPS_OUTPUT_ROOT/figures_cluster/12M/illumination_maps_103310_103926}
PR=${PR:-$COPS_OUTPUT_ROOT/figures_cluster/12M/full_set_stage2_softplus_reg}
NUMBERS_DIR=${NUMBERS_DIR:-$O/checks}
TEMPLATE_RELIEF=${TEMPLATE_RELIEF:-0.7}; GT_RELIEF=${GT_RELIEF:-0}; ONLY=${ONLY:-}
MODEL_NAME=${MODEL_NAME:-the softplus activity model with slope multipliers >= 1, regularized stage-2 chain, final state after 5 fit and fine-tuning cycles}
S=$COPS_OUTPUT_ROOT/figures_cluster/scripts
PY="$APP/.venv/bin/python -W ignore"
stamp() { echo "$(date -u +%FT%TZ) $*"; }
heavy() { bash -c 'echo 1000 > /proc/self/oom_score_adj; exec "$@"' heavy "$@"; }
for f in "$BASE/COPS_model.dat" "$BASE/elements.bin" "$TBASE/COPS_model.dat" "$TBASE/elements.bin" "$TBUNDLE" \
         "$FT/finetune_state.npz" "$FIT/fit_restart.npz" "$FIT/nonneg_facets.npz" "$MAPS/maps_manifest.json" \
         "$GEOM/times.txt" "$MX/nonneg_maps.py" "$NUMBERS_DIR/illumination_numbers.json" \
         "$NUMBERS_DIR/illumination_crosscase.json" "$PR/checks/h2o_${TIME//:/}_manifest.json" \
         "$S/illumination_maps_figures.py" "$S/pub_maps_v2_figures.py" "$S/pub_maps_figures.py"; do
    [ -e "$f" ] || { stamp "missing: $f"; exit 2; }; done
mkdir -p "$O" || exit 1
P=$O/provenance-$SLURM_JOB_ID; mkdir "$P" || exit 1
cd "$APP" || exit 1
stamp "job $SLURM_JOB_ID on $(hostname): illumination maps TIME=$TIME OUT=$O ONLY=$ONLY"
RUN=$R
for v in RUN BASE CH NCYC FT FIT MAPS GEOM MX TBASE TBUNDLE O PR NUMBERS_DIR TIME TEMPLATE_RELIEF GT_RELIEF ONLY MODEL_NAME; do printf '%s=%s\n' "$v" "${!v}"; done > "$P/inputs.txt"
cp -p "$0" "$P/illumination_maps.sbatch.used"
mkdir "$P/code" && cp -p "$S/illumination_maps_figures.py" "$S/pub_maps_v2_figures.py" "$S/pub_maps_figures.py" \
    "$APP/cops_cuda/figures.py" "$MX/nonneg_maps.py" "$P/code/" || exit 1
(cd "$P/code" && sha256sum * > SHA256SUMS)
# read-only git queries for provenance (--no-optional-locks: no index refresh is written)
{ git -C "$AMPS" rev-parse HEAD; git --no-optional-locks -C "$AMPS" status --short; } > "$P/git.txt" 2>&1
{ type module >/dev/null 2>&1 && module list; } > "$P/module_list.txt" 2>&1 || echo "module command not available (no modules loaded; CUDA from $CUDA_HOME)" >> "$P/module_list.txt"
nvidia-smi -q > "$P/nvidia_smi_q.txt" 2>&1
{ nvcc --version; gcc --version | head -1; grep -m1 CSPICE /home/vtenishev/SPICE/cspice/doc/version.txt; lscpu | grep -E 'Model name|^CPU\(s\)'; } > "$P/toolchain.txt" 2>&1
~/.local/bin/uv pip list --python "$APP/.venv/bin/python" > "$P/pip_list.txt" 2>&1 || $PY -m pip list > "$P/pip_list.txt" 2>&1
env | grep -E '^(SLURM_JOB_ID|SLURM_JOB_NODELIST|SLURM_CPUS_PER_TASK|SLURM_MEM_PER_NODE|CUDA_[A-Z_]*|CUPY_[A-Z_]*|TIME|OUT|RUN|CHAIN|NCYC|FT|FIT|MAPS|GEOM|MX|TBASE|TBUNDLE|PR|NUMBERS_DIR|TEMPLATE_RELIEF|GT_RELIEF|ONLY|OMP_NUM_THREADS|MPLBACKEND)=' | sort > "$P/env.txt"
{ echo '#!/bin/bash'; grep '^export ' "$0"; echo "cd $APP"; } > "$P/commands.sh"
run() { printf '%q ' "$@" >> "$P/commands.sh"; echo >> "$P/commands.sh"; "$@"; }
free -g > "$P/free_start.txt" 2>&1; df -h "$COPS_OUTPUT_ROOT" > "$P/df_start.txt" 2>&1
XARGS=(); [ -n "$ONLY" ] && XARGS+=(--only "$ONLY")
run heavy $PY "$S/illumination_maps_figures.py" figures --out-dir "$O" --base "$BASE" --track-base "$TBASE" \
    --bundle "$TBUNDLE" --ft-dir "$FT" --fit-dir "$FIT" --maps-dir "$MAPS" --geometry "$GEOM" \
    --maps-code "$MX/nonneg_maps.py" --time "$TIME" --pr-dir "$PR" --numbers-dir "$NUMBERS_DIR" \
    --template-relief "$TEMPLATE_RELIEF" --gt-relief "$GT_RELIEF" --model-name "$MODEL_NAME" --device 0 \
    "${XARGS[@]}" > "$P/illumination_${TIME//:/}.log" 2>&1; rc=$?
tail -n 30 "$P/illumination_${TIME//:/}.log"
free -g > "$P/free_end.txt" 2>&1; df -h "$COPS_OUTPUT_ROOT" > "$P/df_end.txt" 2>&1
stamp "done: rc $rc"
cp -p "$COPS_OUTPUT_ROOT/jobs/logs/cops_illum_maps-$SLURM_JOB_ID.out" "$P/job.log"
exit $rc
