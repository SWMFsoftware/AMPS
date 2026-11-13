#!/bin/bash
# ==== version-controlled copy: header added 2026-11-13, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/jobs/pub_maps_v2.sbatch (on Voyager: /data001/highresolu_vt/cops_output/jobs/pub_maps_v2.sbatch)
# SHA256 of the run copy = SHA256 of this file without these header lines: 1d4cef12e40f9cac607cad18ecda8d52680fb3349208d74d90a4a8b2d4a16e33
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: job script (GPU): pub_maps_v2_figures.py, one quantity per job; MODEL_NAME added 2026-11-13.
# RUN: QUANTITY=activity | QUANTITY=h2o TIME=HH:MM:SS [OUT, CHAIN, GEOM, MODEL_NAME, LOWCUT] sbatch --export=ALL
#   pub_maps_v2.sbatch
# ==== end of header ====
# Version 2 of the activity / total-production-rate template and ground-track figures of the 12M softplus model after
# 5 cycles (user request 2026-11-12, second round): colour only on the facets seen during the fit, per-panel discrete
# log10 classes (15-20) covering the whole range of the displayed seen non-zero values, total production rate only,
# templates without the header text and every panel as its own figure.
# figures_cluster/scripts/pub_maps_v2_figures.py, ONE quantity per job (QUANTITY=activity, or QUANTITY=total with
# TIME=HH:MM:SS, one of the map times of the maps job), so that they run in parallel on different GPUs / nodes (user:
# "if you need to run any of that on GPU, use all available nodes"). Provenance in OUT/provenance-<job>/ (inputs.txt
# lists every input). Shared node: heavy processes with oom_score_adj 1000 (CLAUDE.md rule 7).
# Every input can be set in the submitting environment (sbatch --export=ALL; values must not contain commas there);
# defaults = the softplus 5-cycle final state of 2026-11-12:
#   RUN      run directory (fit base = $RUN/base)        CHAIN  chain directory          NCYC  last cycle N (c<N>)
#   FT, FIT  final fine tuning / fit (default $CHAIN/c$NCYC/ft/final, $CHAIN/c$NCYC/fit/final)
#   MAPS     maps of that state ($CHAIN/maps)            GEOM   geometry folder of the maps job
#   MX       frozen nonneg_maps code of the maps job     TBASE, TBUNDLE  observation-only base and its bundle (track)
#   OUT      output folder                               QUANTITY, TIME, CLASSES ("15,20")
#   LOWCUT   area fraction cut at the bottom of each panel's colour range (default 0.02 for every quantity, user decisions
#            2026-11-12 "use 2% cut for activity" and "use 2% cut for the rate too"; 0.001 = the range holding 99.9 % of
#            the seen non-zero area; 0 = the entire range)
#   QUANTITY=h2o: the H2O production rate F_H2O = q_H2O A ("plot H2O only for the rate"); needs TIME
#   TEMPLATE_RELIEF  strength of the relief modulation of the class colours in the template figure and its value panels
#            (default 0.7 for QUANTITY=total and h2o, chosen from test renders 0.45 and 0.7 (jobs 124584, 124585); user request 2026-11-12 "the topology of the nucleus is still seen"; 0 for
#            QUANTITY=activity: flat colours; ground tracks always flat)
#   MODEL_NAME  name of the fitted state in the captions (2026-11-13; default "the softplus activity model, final state
#            after 5 fit and fine-tuning cycles"; set it in the submitting shell: it contains commas)
#SBATCH --job-name=cops_pub_maps_v2
#SBATCH --partition=gpu
#SBATCH --account=general
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gpus-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=01:00:00
#SBATCH --output=/data001/highresolu_vt/cops_output/jobs/logs/cops_pub_maps_v2-%j.out
set -uo pipefail
export AMPS=/data001/highresolu_vt/AMPS APP=/data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda
export T=/data001/highresolu_vt/cops_cluster_transfer_20261026 COPS_OUTPUT_ROOT=/data001/highresolu_vt/cops_output
export COPS_DATA_ROOT=$T/rosetta_data_end_of_mission_cops
export CUDA_PATH=/usr/local/cuda-12.9 CUDA_HOME=/usr/local/cuda-12.9 CUPY_CACHE_DIR=$COPS_OUTPUT_ROOT/.cupy
export PATH=$CUDA_HOME/bin:$PATH LD_LIBRARY_PATH=$CUDA_HOME/lib64${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1 MPLBACKEND=Agg
QUANTITY=${QUANTITY:?set QUANTITY=activity, or QUANTITY=total and TIME=HH:MM:SS}; TIME=${TIME:-}; CLASSES=${CLASSES:-15,20}
O=${OUT:-$COPS_OUTPUT_ROOT/figures_cluster/12M/full_set_softplus_5cycles_v2}
R=${RUN:-$COPS_OUTPUT_ROOT/runs/shap8_012m_lategrid_dz1.5}; BASE=$R/base
CH=${CHAIN:-$R/nonneg/softplus/chain_5cycles_20261112T012553Z}; NCYC=${NCYC:-5}
FT=${FT:-$CH/c$NCYC/ft/final}; FIT=${FIT:-$CH/c$NCYC/fit/final}; MAPS=${MAPS:-$CH/maps}
GEOM=${GEOM:-$R/maps_geometry/job_124324/geometry}; MX=${MX:-$COPS_OUTPUT_ROOT/model_ext_frozen/nonneg_maps_20261112T021623Z}
TBASE=${TBASE:-$COPS_OUTPUT_ROOT/runs/shap8_012m_flat_earlier/base}
TBUNDLE=${TBUNDLE:-$COPS_OUTPUT_ROOT/runs/shap8_012m_flat_earlier/bundle/restart_templates.npz}
LOWCUT=${LOWCUT:-0.02}
MODEL_NAME=${MODEL_NAME:-the softplus activity model, final state after 5 fit and fine-tuning cycles}
if [ "$QUANTITY" = total ] || [ "$QUANTITY" = h2o ]; then TEMPLATE_RELIEF=${TEMPLATE_RELIEF:-0.7}; else TEMPLATE_RELIEF=${TEMPLATE_RELIEF:-0}; fi
QARGS=(--quantity "$QUANTITY" --low-cut-area "$LOWCUT" --template-relief "$TEMPLATE_RELIEF"); [ -n "$TIME" ] && QARGS+=(--time "$TIME"); TAG=$QUANTITY${TIME:+_${TIME//:/}}
S=$COPS_OUTPUT_ROOT/figures_cluster/scripts
PY="$APP/.venv/bin/python -W ignore"
stamp() { echo "$(date -u +%FT%TZ) $*"; }
heavy() { bash -c 'echo 1000 > /proc/self/oom_score_adj; exec "$@"' heavy "$@"; }
for f in "$BASE/COPS_model.dat" "$BASE/elements.bin" "$TBASE/COPS_model.dat" "$TBASE/elements.bin" "$TBUNDLE" \
         "$FT/finetune_state.npz" "$FIT/fit_restart.npz" "$FIT/nonneg_facets.npz" "$MAPS/maps_manifest.json" \
         "$GEOM/times.txt" "$MX/nonneg_maps.py" "$S/pub_maps_v2_figures.py" "$S/pub_maps_figures.py"; do
    [ -e "$f" ] || { stamp "missing: $f"; exit 2; }; done
mkdir -p "$O" || exit 1
P=$O/provenance-$SLURM_JOB_ID; mkdir "$P" || exit 1
cd "$APP" || exit 1
stamp "job $SLURM_JOB_ID on $(hostname): QUANTITY=$QUANTITY TIME=$TIME CLASSES=$CLASSES LOWCUT=$LOWCUT TEMPLATE_RELIEF=$TEMPLATE_RELIEF OUT=$O"
RUN=$R
for v in RUN BASE CH NCYC FT FIT MAPS GEOM MX TBASE TBUNDLE O QUANTITY TIME CLASSES LOWCUT TEMPLATE_RELIEF MODEL_NAME; do printf '%s=%s\n' "$v" "${!v}"; done > "$P/inputs.txt"
cp -p "$0" "$P/pub_maps_v2.sbatch.used"
mkdir "$P/code" && cp -p "$S/pub_maps_v2_figures.py" "$S/pub_maps_figures.py" "$APP/cops_cuda/figures.py" \
    "$MX/nonneg_maps.py" "$P/code/" || exit 1
(cd "$P/code" && sha256sum * > SHA256SUMS)
{ git -C "$AMPS" rev-parse HEAD; git -C "$AMPS" status --short; } > "$P/git.txt" 2>&1
{ type module >/dev/null 2>&1 && module list; } > "$P/module_list.txt" 2>&1 || echo "module command not available (no modules loaded; CUDA from $CUDA_HOME)" >> "$P/module_list.txt"
nvidia-smi -q > "$P/nvidia_smi_q.txt" 2>&1
{ nvcc --version; gcc --version | head -1; grep -m1 CSPICE /home/vtenishev/SPICE/cspice/doc/version.txt; lscpu | grep -E 'Model name|^CPU\(s\)'; } > "$P/toolchain.txt" 2>&1
~/.local/bin/uv pip list --python "$APP/.venv/bin/python" > "$P/pip_list.txt" 2>&1 || $PY -m pip list > "$P/pip_list.txt" 2>&1
env | grep -E '^(SLURM_JOB_ID|SLURM_JOB_NODELIST|SLURM_CPUS_PER_TASK|SLURM_MEM_PER_NODE|CUDA_[A-Z_]*|CUPY_[A-Z_]*|QUANTITY|TIME|CLASSES|LOWCUT|TEMPLATE_RELIEF|OUT|RUN|CHAIN|NCYC|FT|FIT|MAPS|GEOM|MX|TBASE|TBUNDLE|OMP_NUM_THREADS|MPLBACKEND)=' | sort > "$P/env.txt"
{ echo '#!/bin/bash'; grep '^export ' "$0"; echo "cd $APP"; } > "$P/commands.sh"
run() { printf '%q ' "$@" >> "$P/commands.sh"; echo >> "$P/commands.sh"; "$@"; }
free -g > "$P/free_start.txt" 2>&1; df -h "$COPS_OUTPUT_ROOT" > "$P/df_start.txt" 2>&1
run heavy $PY "$S/pub_maps_v2_figures.py" --base "$BASE" --track-base "$TBASE" --bundle "$TBUNDLE" --ft-dir "$FT" \
    --fit-dir "$FIT" --maps-dir "$MAPS" --geometry "$GEOM" --maps-code "$MX/nonneg_maps.py" --out-dir "$O" \
    "${QARGS[@]}" --classes "$CLASSES" --model-name "$MODEL_NAME" --device 0 > "$P/v2_$TAG.log" 2>&1; rc=$?
tail -n 25 "$P/v2_$TAG.log"
free -g > "$P/free_end.txt" 2>&1; df -h "$COPS_OUTPUT_ROOT" > "$P/df_end.txt" 2>&1
stamp "done: rc $rc"
cp -p "$COPS_OUTPUT_ROOT/jobs/logs/cops_pub_maps_v2-$SLURM_JOB_ID.out" "$P/job.log"
exit $rc
