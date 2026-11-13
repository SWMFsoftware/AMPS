#!/bin/bash
# ==== version-controlled copy: header added 2026-11-13, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/jobs/pub_maps_figures.sbatch (on Voyager: /data001/highresolu_vt/cops_output/jobs/pub_maps_figures.sbatch)
# SHA256 of the run copy = SHA256 of this file without these header lines: dbbfd55bc1ea4af2f42200208f781e60a8461ffa096dc94aa9a052b1da19c951
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: job script (GPU): full figure set 3.4 (set1 visibility, set4 model/data, pressures, maps with fixed ranges);
#   inputs from the environment since 2026-11-13.
# RUN: STEPS=set1,set4,pressures CHAIN=... PUB=... PUBSTEM=... TAG=... MODEL_LABEL=... sbatch --export=ALL
#   pub_maps_figures.sbatch
# ==== end of header ====
# Full publication figure set of the 12M softplus activity model after 5 fit / fine-tuning cycles (final state).
# User request 2026-11-12: "generate the full set of plots that you have done before (templates, gravity, ground
# tracks, model/data comparison, etc.); i need three sets of the template and groundtrack plots: 1) as now, 2) surface
# activity map instead of visibility maps, 3) production rate [m^-2 s^-1] (H2O, CO2, total) instead of visibility
# maps". Answers: the fit = softplus after 5 cycles (final state); production-rate times 10:33:10 and 10:39:26 UTC;
# colour scale log10 over the non-zero range, zero in a separate grey colour.
#
# Steps (STEPS, comma list; default set1,set4,pressures,maps):
#   set1  visibility/: tools/pub_figures.py templates + groundtrack (UNCHANGED tracked tool) on the observation-only
#         base TBASE (runs/shap8_012m_flat_earlier/base: the 1150 observations, no gauge offset, the template bundle's
#         own base; user decision 2026-11-12, as in the earlier late-visibility figures)
#   set4  model_data/: positive_table.py (observation rows of the fit and of the clipped reference refit),
#         model_obs_figures.py; copy of figures_cluster/12M/publication_softplus
#   pressures  model_data/: tools/pub_figures.py pressures (UNCHANGED tracked tool) with the set4 tables; short
#         labels and title (its legend has 3 columns; long labels widen the figure)
#   maps  activity/, production_rate/, checks/ and the set-1 ground-track captions: figures_cluster/scripts/
#         pub_maps_figures.py (SETS: activity,production), fields of the fit on the lategrid base, ground track and
#         impact from TBASE (the same as set 1), fixed colour ranges (user decision 2026-11-12: log10 A 13..23,
#         F_H2O 12..22, F_CO2 -28..-18, total 12..22; non-zero below the range green); PREVIEW=1 adds the 10-decade
#         percentile proposal (colour_range_proposal/, not used since the decision)
# Products: OUT (default figures_cluster/12M/full_set_softplus_5cycles); provenance in OUT/provenance-<job>/.
# Inputs from the environment (2026-11-13, for other fits; defaults = the softplus 5-cycle fit of 2026-11-12): CHAIN,
#   NCYC (FT, FIT = c<NCYC>/ft/final, fit/final), MAPS, GEOM, PUB (publication-figure folder), PUBSTEM (its figure stem,
#   default 12M_softplus_NG_RG_altitude), PUBTABLE (its table, default softplus_c5ft.dat), TAG (file-name tag, default
#   softplus_5cycles), MODEL_LABEL (legend label of the fit, without its cost), MODEL_SHORT (short legend label of the
#   pressures figure), TITLE_PRESSURES.
# Nothing is overwritten (the tools refuse existing files); GPU: 1 (rasterizer, ray casting), shared node.
#SBATCH --job-name=cops_pub_maps
#SBATCH --partition=gpu
#SBATCH --account=general
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gpus-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=01:00:00
#SBATCH --exclude=gpu003
#SBATCH --output=/data001/highresolu_vt/cops_output/jobs/logs/cops_pub_maps-%j.out
set -uo pipefail
export AMPS=/data001/highresolu_vt/AMPS APP=/data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda
export T=/data001/highresolu_vt/cops_cluster_transfer_20261026 COPS_OUTPUT_ROOT=/data001/highresolu_vt/cops_output
export COPS_DATA_ROOT=$T/rosetta_data_end_of_mission_cops
export CUDA_PATH=/usr/local/cuda-12.9 CUDA_HOME=/usr/local/cuda-12.9 CUPY_CACHE_DIR=$COPS_OUTPUT_ROOT/.cupy
export PATH=$CUDA_HOME/bin:$PATH LD_LIBRARY_PATH=$CUDA_HOME/lib64${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1 MPLBACKEND=Agg
STEPS=${STEPS:-set1,set4,pressures,maps}; SETS=${SETS:-activity,production}; PREVIEW=${PREVIEW:-0}
O=${OUT:-$COPS_OUTPUT_ROOT/figures_cluster/12M/full_set_softplus_5cycles}
R=$COPS_OUTPUT_ROOT/runs/shap8_012m_lategrid_dz1.5; BASE=$R/base; BUNDLE=$R/bundle/restart_templates.npz
CH=${CHAIN:-$R/nonneg/softplus/chain_5cycles_20261112T012553Z}; NCYC=${NCYC:-5}
FT=${FT:-$CH/c$NCYC/ft/final}; FIT=${FIT:-$CH/c$NCYC/fit/final}; MAPS=${MAPS:-$CH/maps}
GEOM=${GEOM:-$R/maps_geometry/job_124324/geometry}; MX=$COPS_OUTPUT_ROOT/model_ext_frozen/nonneg_maps_20261112T021623Z
D2=$R/clip_refit/D2/final; S=$COPS_OUTPUT_ROOT/figures_cluster/scripts
TBASE=$COPS_OUTPUT_ROOT/runs/shap8_012m_flat_earlier/base; TBUNDLE=$COPS_OUTPUT_ROOT/runs/shap8_012m_flat_earlier/bundle/restart_templates.npz
RANGES=(--range A=13,23 --range F_H2O=12,22 --range F_CO2=-28,-18 --range total=12,22)
PUB=${PUB:-$COPS_OUTPUT_ROOT/figures_cluster/12M/publication_softplus}; PUBSTEM=${PUBSTEM:-12M_softplus_NG_RG_altitude}
PUBTABLE=${PUBTABLE:-softplus_c5ft.dat}; TAG=${TAG:-softplus_5cycles}
MODEL_LABEL=${MODEL_LABEL:-softplus activity model after 5 fit and fine-tuning cycles}
MODEL_SHORT=${MODEL_SHORT:-softplus activity model}
TITLE_PRESSURES=${TITLE_PRESSURES:-SHAP8 12M: softplus activity model after 5 cycles vs COPS}
PY="$APP/.venv/bin/python -W ignore"
stamp() { echo "$(date -u +%FT%TZ) $*"; }
heavy() { bash -c 'echo 1000 > /proc/self/oom_score_adj; exec "$@"' heavy "$@"; }   # CLAUDE.md rule 7 (shared node)
for f in "$BASE/COPS_model.dat" "$BASE/elements.bin" "$BUNDLE" "$TBASE/COPS_model.dat" "$TBASE/elements.bin" "$TBUNDLE" "$FT/finetune_state.npz" "$FIT/fit_restart.npz" \
         "$FIT/nonneg_facets.npz" "$MAPS/maps_manifest.json" "$GEOM/times.txt" "$MX/nonneg_maps.py" "$D2/COPS_model.dat" \
         "$S/pub_maps_figures.py" "$S/positive_table.py" "$S/model_obs_figures.py" "$PUB/$PUBSTEM.png"; do
    [ -e "$f" ] || { stamp "missing: $f"; exit 2; }; done
mkdir -p "$O" || exit 1
P=$O/provenance-$SLURM_JOB_ID; mkdir "$P" || exit 1
cd "$APP" || exit 1
stamp "job $SLURM_JOB_ID on $(hostname): STEPS=$STEPS SETS=$SETS PREVIEW=$PREVIEW OUT=$O"

# ---------------------------------------------------------------- provenance (CLAUDE.md rule 8)
cp -p "$0" "$P/pub_maps_figures.sbatch.used"
mkdir "$P/code" && cp -p "$S/pub_maps_figures.py" "$S/positive_table.py" "$S/model_obs_figures.py" "$APP/tools/pub_figures.py" \
    "$APP/cops_cuda/figures.py" "$MX/nonneg_maps.py" "$P/code/" || exit 1
(cd "$P/code" && sha256sum * > SHA256SUMS)
{ git -C "$AMPS" rev-parse HEAD; git -C "$AMPS" status --short; } > "$P/git.txt" 2>&1
{ type module >/dev/null 2>&1 && module list; } > "$P/module_list.txt" 2>&1 || echo "module command not available (no modules loaded; CUDA from $CUDA_HOME)" >> "$P/module_list.txt"
nvidia-smi -q > "$P/nvidia_smi_q.txt" 2>&1
{ nvcc --version; gcc --version | head -1; grep -m1 CSPICE /home/vtenishev/SPICE/cspice/doc/version.txt; lscpu | grep -E 'Model name|^CPU\(s\)'; } > "$P/toolchain.txt" 2>&1
~/.local/bin/uv pip list --python "$APP/.venv/bin/python" > "$P/pip_list.txt" 2>&1 || $PY -m pip list > "$P/pip_list.txt" 2>&1
env | grep -E '^(SLURM_JOB_ID|SLURM_JOB_NODELIST|SLURM_CPUS_PER_TASK|CUDA_[A-Z_]*|CUPY_[A-Z_]*|STEPS|SETS|PREVIEW|OUT|OMP_NUM_THREADS|MPLBACKEND)=' | sort > "$P/env.txt"
for v in CH NCYC FT FIT MAPS GEOM PUB PUBSTEM PUBTABLE TAG MODEL_LABEL MODEL_SHORT TITLE_PRESSURES O STEPS SETS; do printf '%s=%s\n' "$v" "${!v}"; done > "$P/inputs.txt"
{ echo '#!/bin/bash'; grep '^export ' "$0"; echo "cd $APP"; } > "$P/commands.sh"
run() { printf '%q ' "$@" >> "$P/commands.sh"; echo >> "$P/commands.sh"; "$@"; }
free -g > "$P/free_start.txt" 2>&1; df -h "$COPS_OUTPUT_ROOT" > "$P/df_start.txt" 2>&1
nfail=0

# ---------------------------------------------------------------- set 1: visibility (unchanged tool)
if [[ ",$STEPS," == *,set1,* ]]; then
    stamp "set 1: tools/pub_figures.py templates + groundtrack"
    mkdir -p "$O/visibility"
    run heavy $PY tools/pub_figures.py templates --base "$TBASE" --bundle "$TBUNDLE" \
        --out "$O/visibility/12M_templates_visibility" --device 0 > "$P/set1_templates.log" 2>&1 || { stamp "set1 templates FAILED"; nfail=$((nfail+1)); }
    tail -n 5 "$P/set1_templates.log"
    run heavy $PY tools/pub_figures.py groundtrack --base "$TBASE" --bundle "$TBUNDLE" \
        --out "$O/visibility/12M_groundtrack_visibility" --device 0 > "$P/set1_groundtrack.log" 2>&1 || { stamp "set1 groundtrack FAILED"; nfail=$((nfail+1)); }
    tail -n 8 "$P/set1_groundtrack.log"
fi

# ---------------------------------------------------------------- set 4: model / data
if [[ ",$STEPS," == *,set4,* ]]; then
    stamp "set 4: model/data comparison"
    M=$O/model_data; TB=$M/tables; mkdir -p "$TB"
    run $PY "$S/positive_table.py" "$FT" "$TB/${TAG}_final.dat" || nfail=$((nfail+1))
    run $PY "$S/positive_table.py" "$D2" "$TB/clipped_refit_reference.dat" || nfail=$((nfail+1))
    CS=$(awk '/^cost_data/{printf "%.4e", $5}' "$FT/finetune_result.txt")
    CD=$(awk '/^final_data_cost/{printf "%.4e", $2}' "$D2/fit_result.txt")
    LS="$MODEL_LABEL, data cost $CS"
    LD="clipped activity model refitted on the same base (reference), data cost $CD"
    TI="SHAP8 12M, COPS 1.5 m toward the comet, 10-s and 5-readout response from 10:33:10 UTC: model vs COPS"
    run $PY "$S/model_obs_figures.py" --bundle "$BUNDLE" --ref "$TB/${TAG}_final.dat" \
        --out "$M/12M_${TAG}_vs_COPS" --title "$TI" \
        --fit "$LS=$TB/${TAG}_final.dat:solid" --fit "$LD=$TB/clipped_refit_reference.dat:dashed" \
        > "$P/set4_model_obs.log" 2>&1 || { stamp "model_obs_figures FAILED"; nfail=$((nfail+1)); }
    tail -n 12 "$P/set4_model_obs.log"
    PD=$M/publication_NG_RG_altitude; mkdir -p "$PD"
    run cp -p "$PUB/$PUBSTEM.eps" "$PUB/$PUBSTEM.png" "$PUB/$PUBTABLE" "$PD/" || nfail=$((nfail+1))
    run cp -p "$PUB/HOW-TO-RERUN.md" "$PD/HOW-TO-RERUN_$(basename "$PUB").md" || nfail=$((nfail+1))
    run cp -rp "$PUB/code" "$PD/code" || nfail=$((nfail+1))
    (cd "$PD" && sha256sum $PUBSTEM.eps $PUBSTEM.png > COPIED_FROM.sha256
     echo "# copied by job $SLURM_JOB_ID on $(date -u +%FT%TZ) from $PUB (cp -p; not regenerated)" >> COPIED_FROM.sha256
     (cd "$PUB" && sha256sum $PUBSTEM.eps $PUBSTEM.png) | sed 's/^/# source: /' >> COPIED_FROM.sha256)
fi

# ---------------------------------------------------------------- set 4: tools/pub_figures.py pressures
if [[ ",$STEPS," == *,pressures,* ]]; then
    stamp "set 4: tools/pub_figures.py pressures"
    M=$O/model_data; TB=$M/tables
    run $PY tools/pub_figures.py pressures --base "$BASE" --bundle "$BUNDLE" --no-base \
        --ref "$TB/${TAG}_final.dat" --out "$M/12M_${TAG}_pressures" \
        --title "$TITLE_PRESSURES" \
        --fit "$MODEL_SHORT=$TB/${TAG}_final.dat:solid" \
        --fit "clipped model (reference)=$TB/clipped_refit_reference.dat:dashed" \
        > "$P/set4_pressures.log" 2>&1 || { stamp "pub_figures pressures FAILED"; nfail=$((nfail+1)); }
    tail -n 4 "$P/set4_pressures.log"
fi

# ---------------------------------------------------------------- sets 2 and 3: activity and production rate
if [[ ",$STEPS," == *,maps,* ]]; then
    stamp "sets 2 and 3: pub_maps_figures.py (SETS=$SETS, PREVIEW=$PREVIEW)"
    EXTRA=()
    [ -e "$O/visibility/12M_groundtrack_visibility_3d.png" ] && [ ! -e "$O/visibility/12M_groundtrack_visibility_3d_caption.txt" ] && \
        EXTRA+=(--visibility-prefix "$O/visibility/12M_groundtrack_visibility")
    [ "$PREVIEW" = 1 ] && EXTRA+=(--preview-dir "$O/colour_range_proposal" --span-decades 10 --top-percentile 99.99)
    run heavy $PY "$S/pub_maps_figures.py" --base "$BASE" --track-base "$TBASE" --bundle "$TBUNDLE" "${RANGES[@]}" \
        --ft-dir "$FT" --fit-dir "$FIT" \
        --maps-dir "$MAPS" --geometry "$GEOM" --maps-code "$MX/nonneg_maps.py" --out-dir "$O" --sets "$SETS" \
        --device 0 "${EXTRA[@]}" > "$P/maps.log" 2>&1 || { stamp "pub_maps_figures FAILED"; nfail=$((nfail+1)); }
    tail -n 30 "$P/maps.log"
fi
free -g > "$P/free_end.txt" 2>&1; df -h "$COPS_OUTPUT_ROOT" > "$P/df_end.txt" 2>&1
stamp "done: $nfail failures"
cp -p "$COPS_OUTPUT_ROOT/jobs/logs/cops_pub_maps-$SLURM_JOB_ID.out" "$P/job.log"
[ $nfail -eq 0 ]
