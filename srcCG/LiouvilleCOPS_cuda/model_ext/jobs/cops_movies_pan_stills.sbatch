#!/bin/bash
# ==== version-controlled copy: header added 2026-11-13, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/jobs/cops_movies_pan_stills.sbatch (on Voyager: /data001/highresolu_vt/cops_output/jobs/cops_movies_pan_stills.sbatch)
# SHA256 of the run copy = SHA256 of this file without these header lines: b973184621471a61d58d8c79aa080148adad7514ca55f4e89f2ae719e82e0b33
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: Candidate stills of the panorama movie: per pose a work folder (detach-300 m products copied), plan --pan-to,
#   touchdown, render of the mid-path and last frames, compose --no-sub-mark.
# RUN: sbatch --nodelist=gpu002 --export=ALL (CHAIN, POSES, OUTS, MOVIES_PY)
# ==== end of header ====
# Panorama movie, candidate final poses (2026-11-14, user ~20:05 CST): per candidate pose (one GPU each, in parallel):
# work folder $CHAIN/movies_work_c9_detach300_pan_<id> with the detach-300 m products copied (cp -p + cmp), plan --camera
# detach --detach-altitude 300 --fov-deg 90 --pan-to <pose>, touchdown (flag pixels per frame), render of the mid-path
# frame and the last frame only (cone masks, weight tints, spacecraft shadow), compose of those frames (touchdown flag on
# top, no sub-spacecraft "+"); the composed stills copied to $OUTS/<id>/. Stills only: no movie.
# Environment: CHAIN, POSES (space-separated pose json files; id = basename without .json), OUTS, MOVIES_PY, TABLE, TS_LABEL.
#SBATCH --job-name=cops_pan_stills
#SBATCH --partition=gpu
#SBATCH --account=general
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gpus-per-node=3
#SBATCH --cpus-per-task=24
#SBATCH --mem=150G
#SBATCH --time=01:30:00
#SBATCH --output=/data001/highresolu_vt/cops_output/jobs/logs/cops_pan_stills-%j.out
set -uo pipefail
export AMPS=/data001/highresolu_vt/AMPS APP=/data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda
export COPS_OUTPUT_ROOT=/data001/highresolu_vt/cops_output
export CUDA_PATH=/usr/local/cuda-12.9 CUDA_HOME=/usr/local/cuda-12.9 CUPY_CACHE_DIR=$COPS_OUTPUT_ROOT/.cupy
export PATH=$CUDA_HOME/bin:$PATH LD_LIBRARY_PATH=$CUDA_HOME/lib64${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1 MPLBACKEND=Agg
R=$COPS_OUTPUT_ROOT/runs/shap8_012m_lategrid_dz1.5; GEOM=${GEOM:-$R/maps_geometry/job_124775/geometry}
CH=${CHAIN:?}; FROM=${FROM_WORK:-$CH/movies_work_c9_detach300}; POSES=${POSES:?}; OUTS=${OUTS:?}
S=$COPS_OUTPUT_ROOT/figures_cluster/scripts; MPY=${MOVIES_PY:-$S/cops_movies.py}
F12=$COPS_OUTPUT_ROOT/figures_cluster/12M
TABLE=${TABLE:-$F12/model_vs_observations_softplus_smm1_9cycles/tables/softplus_smm1_c9ft.dat}
TS_LABEL=${TS_LABEL:-'COPS data and the\nsoftplus activity model,\nmultipliers ≥ 1\n(final state, 9 cycles)'}
TRACK_BASE=$R/../shap8_012m_flat_earlier/base; IMPACT_JSON=$F12/full_set_softplus_smm1_5cycles/visibility/12M_groundtrack_visibility_impact.json
PY="$APP/.venv/bin/python -W ignore"
stamp() { echo "$(date -u +%FT%TZ) $*"; }
mkdir -p "$OUTS" || exit 1
P=$OUTS/provenance-$SLURM_JOB_ID; mkdir "$P" || exit 1
cd "$APP" || exit 1
stamp "job $SLURM_JOB_ID on $(hostname), CUDA_VISIBLE_DEVICES=${CUDA_VISIBLE_DEVICES:-}"
for v in CH FROM POSES OUTS MPY TABLE TS_LABEL GEOM TRACK_BASE IMPACT_JSON; do printf '%s=%s\n' "$v" "${!v}"; done > "$P/inputs.txt"
cp -p "$0" "$P/cops_movies_pan_stills.sbatch.used"
mkdir "$P/code" "$P/poses" && cp -p "$MPY" "$APP/cops_cuda/figures.py" "$APP/cops_cuda/kernels/render.cu" "$P/code/" && cp -p $POSES "$P/poses/" \
    && (cd "$P/code" && sha256sum * > SHA256SUMS) || exit 1
{ git -C "$AMPS" rev-parse HEAD; git -C "$AMPS" status --short; } > "$P/git.txt" 2>&1
{ type module >/dev/null 2>&1 && module list; } > "$P/module_list.txt" 2>&1 || echo "module command not available (CUDA from $CUDA_HOME)" >> "$P/module_list.txt"
nvidia-smi -q > "$P/nvidia_smi_q.txt" 2>&1
{ nvcc --version; gcc --version | head -1; lscpu | grep -E 'Model name|^CPU\(s\)'; } > "$P/toolchain.txt" 2>&1
~/.local/bin/uv pip list --python "$APP/.venv/bin/python" > "$P/pip_list.txt" 2>&1
heavy() { bash -c 'echo 1000 > /proc/self/oom_score_adj; exec "$@"' heavy "$@"; }
one() {         # one candidate on device $2
    local pose=$1 dev=$2 id W L
    id=$(basename "$pose" .json); W=$CH/movies_work_c9_detach300_pan_$id; L=$P/$id; mkdir -p "$L"
    [ -e "$W" ] && { stamp "$W exists"; return 1; }
    mkdir "$W" || return 1
    for f in prep.npz prep_checks.json fov.npz fov_checks.json sc.npz sc_checks.json kernels.txt frame_utc.txt rows_utc.txt \
             maps_utc.txt frame_geom.txt rows_geom.txt maps_geom.txt shield_table.txt contact.json; do
        cp -p "$FROM/$f" "$W/$f" && cmp -s "$FROM/$f" "$W/$f" || { stamp "$id: copy of $f failed"; return 1; }
    done
    (cd "$FROM" && sha256sum prep.npz fov.npz sc.npz contact.json) > "$L/copied_from_$(basename "$FROM").sha256"
    local cmds=$L/commands.sh
    { echo '#!/bin/bash'; grep '^export ' "$0"; echo "cd $APP"; } > "$cmds"
    r() { printf '%q ' "$@" >> "$cmds"; echo >> "$cmds"; heavy "$@"; }
    r $PY "$MPY" plan --work "$W" --geometry "$GEOM" --device "$dev" --camera detach --detach-altitude 300 --fov-deg 90 \
        --pan-to "$pose" > "$L/plan.log" 2>&1 || { stamp "$id: plan FAILED"; return 1; }
    r $PY "$MPY" touchdown --work "$W" --geometry "$GEOM" --device "$dev" --track-base "$TRACK_BASE" --impact-json "$IMPACT_JSON" \
        > "$L/touchdown.log" 2>&1 || { stamp "$id: touchdown FAILED"; return 1; }
    local km kl
    km=$($PY -c "import json;print(json.load(open('$W/camera_checks.json'))['pan']['mid_frame'])")
    kl=$($PY -c "import json;print(json.load(open('$W/camera_checks.json'))['frames']-1)")
    r $PY "$MPY" render --work "$W" --geometry "$GEOM" --device "$dev" --worker 0 --nworkers 1 --only "$km,$kl" \
        --mask-rg cone:45 --mask-ng cone:77.5 --tint-rg weight --tint-ng weight --sc-shadow --sc-samples 64 \
        > "$L/render.log" 2>&1 || { stamp "$id: render FAILED"; return 1; }
    r $PY "$MPY" compose --work "$W" --out "$W/composed" --table "$TABLE" --procs 2 --frames "$km,$kl" --ts-label "$TS_LABEL" \
        --touchdown "$W/touchdown.npz" --flag-on-top --no-sub-mark > "$L/compose.log" 2>&1 || { stamp "$id: compose FAILED"; return 1; }
    mkdir -p "$OUTS/$id"
    for k in $km $kl; do for m in 1 2 3; do
        cp -p "$W/composed/movie$m/f$(printf %05d $k).png" "$OUTS/$id/${id}_movie${m}_f$(printf %05d $k).png"; done; done
    cp -p "$W/camera_checks.json" "$OUTS/$id/camera_checks.json"; cp -p "$pose" "$OUTS/$id/"
    stamp "$id: done (frames $km, $kl)"
}
i=0; pids=()
for pose in $POSES; do one "$pose" $i > "$P/one_$(basename "$pose" .json).log" 2>&1 & pids+=($!); i=$((i+1)); done
rc=0; for p_ in "${pids[@]}"; do wait $p_ || rc=1; done
cat "$P"/one_*.log
stamp "done: rc $rc"
cp -p "$COPS_OUTPUT_ROOT/jobs/logs/cops_pan_stills-$SLURM_JOB_ID.out" "$P/job.log"
exit $rc
