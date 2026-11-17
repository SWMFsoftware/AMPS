#!/bin/bash
# ==== version-controlled copy: header added 2026-11-15, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/jobs/localnn_tools.sbatch (on Voyager: /data001/highresolu_vt/cops_output/jobs/localnn_tools.sbatch)
# SHA256 of the run copy = SHA256 of this file without these header lines: 367514075fd9001c4e883130187821ae9b279c491f22e37f3a9b8993d49219e9
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: job script (CPU, option C): export + the frozen contribution and maps tools on one localnn final state.
# RUN: tests 125020, 125021, 125031; stage-1 125053, 125057, 125061 (chained).
# 2026-11-15: 2026-11-17: NPY_DISABLE_CPU_FEATURES (AVX-512 off) as the chain driver: on the AVX-512 node028 (EPYC
#   9554) np.exp differed in the last bit at 38 of the 271 held facets and the bitwise reproduction of A failed
#   (post 136094).
# RUN 2026-11-15: 136123 (the f10 / 1e-4 post again), 136097.
# ==== end of header ====
# Option C (local non-negative basis): the maps and contribution tools of the nonneg chains on ONE localnn final state
# (model_ext/localnn_tools.py: export folder + the frozen nonneg_maps.py / nonneg_contrib.py, unchanged, with their
# state_activity replaced by the localnn one). CPU (normal128cores).
#   1. export  -> $RUNROOT/tools_state (use as --fit-dir and --ft-dir; NOT a fine tuning, see its README.txt)
#   2. contrib -> $RUNROOT/contributions (P% source areas, occurrence, slope-bin fractions; verification vs the tables)
#   3. maps    -> $RUNROOT/maps (activity and source-rate maps, .vtu + Tecplot, 13 times)
# Settings: RUNROOT (required), STEPS (default "export contrib maps"), MODEL_EXT (the frozen localnn snapshot), GRID_DIR,
#   GEOM (default maps_geometry/job_124775/geometry), CONTRIB (frozen nonneg_contrib folder), PROCS,
#   FT_DIR (stage 2, 2026-11-15: a localnn fine tuning of RUNROOT/fit/final, e.g. <chain>/c5/ft/final: the tools run on
#   its final state A_final; the folders get the suffix _ft: tools_state_ft, contributions_ft, maps_ft).
#SBATCH --job-name=cops_localnn_tools
#SBATCH --partition=normal128cores
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=320G
#SBATCH --time=02:00:00
#SBATCH --output=/data001/highresolu_vt/cops_output/jobs/logs/cops_localnn_tools-%j.out
set -uo pipefail
# numpy without AVX-512 kernels, as the chain driver (2026-11-17: on an AVX-512 node (node028, EPYC 9554) np.exp differs
# in the last bit at 38 of the 271 held facets and the tools' bitwise reproduction of A failed, post 136094)
export NPY_DISABLE_CPU_FEATURES="AVX512F AVX512CD AVX512_KNL AVX512_KNM AVX512_SKX AVX512_CLX AVX512_CNL AVX512_ICL"
export AMPS=/data001/highresolu_vt/AMPS APP=/data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda
export COPS_OUTPUT_ROOT=/data001/highresolu_vt/cops_output
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1
UNSET=$(env | grep -o -E '^(COPS_POS_[A-Z_]*|COPS_SLOPEMUL_[A-Z_]*|COPS_NONNEG_[A-Z_]*)' | sort | paste -sd' ' -)
[ -n "$UNSET" ] && unset $UNSET
R=$COPS_OUTPUT_ROOT/runs/shap8_012m_lategrid_dz1.5
BASE=$R/base; D2=$R/clip_refit/D2/final; PT=$R/late_grid_points.txt
export MODEL_EXT=$(readlink -f "${MODEL_EXT:-$COPS_OUTPUT_ROOT/model_ext}")
GRID_DIR=${GRID_DIR:-$R/nonneg/localnn/grid_150m}
GEOM=${GEOM:-$R/maps_geometry/job_124775/geometry}
CONTRIB=${CONTRIB:-$COPS_OUTPUT_ROOT/model_ext_frozen/nonneg_contrib_20261112T210012Z}
RUNROOT=$(readlink -f "${RUNROOT:?RUNROOT is required}")
STEPS=${STEPS:-export contrib maps}
FT_DIR=${FT_DIR:-}; SUF=; FTO=()
[ -n "$FT_DIR" ] && { FT_DIR=$(readlink -f "$FT_DIR"); SUF=_ft; FTO=(--ft-dir "$FT_DIR"); }
PROCS=${PROCS:-${SLURM_CPUS_PER_TASK:-8}}
PY=$APP/.venv/bin/python
stamp() { echo "$(date -u +%FT%TZ) $*"; }
heavy() { bash -c 'echo 1000 > /proc/self/oom_score_adj; exec "$@"' heavy "$@"; }
P=$RUNROOT/provenance-$SLURM_JOB_ID-tools; mkdir "$P" || exit 1
cd "$APP" || exit 1
stamp "job $SLURM_JOB_ID on $(hostname): RUNROOT=$RUNROOT STEPS='$STEPS' FT_DIR=$FT_DIR MODEL_EXT=$MODEL_EXT CONTRIB=$CONTRIB"
cp -p "$0" "$P/localnn_tools.sbatch.used"
cp -p "$MODEL_EXT"/localnn_*.py "$MODEL_EXT/positive_problem.py" "$P/"; mkdir "$P/contrib_code" && cp -p "$CONTRIB"/* "$P/contrib_code/"
( cd "$P" && sha256sum *.py contrib_code/*.py localnn_tools.sbatch.used > code_sha256.txt )
{ git -C "$AMPS" rev-parse HEAD; git -C "$AMPS" status --short; } > "$P/git.txt" 2>&1
{ type module >/dev/null 2>&1 && module list; } > "$P/module_list.txt" 2>&1 || echo "module command not available (no modules loaded)" >> "$P/module_list.txt"
echo "CPU job on $(hostname): no GPU used" > "$P/nvidia_smi_q.txt"
{ gcc --version | head -1; lscpu | grep -E 'Model name|^CPU\(s\)'; } > "$P/toolchain.txt" 2>&1
~/.local/bin/uv pip list --python "$PY" > "$P/pip_list.txt" 2>&1 || "$PY" -m pip list > "$P/pip_list.txt" 2>&1
{ echo '#!/bin/bash'; grep '^export ' "$0"; echo "export MODEL_EXT=$MODEL_EXT"; echo "cd $APP"; } > "$P/commands.sh"
run() { printf '%q ' "$@" >> "$P/commands.sh"; echo >> "$P/commands.sh"; "$@"; }
S=$RUNROOT/tools_state$SUF
rc=0
for step in $STEPS; do
  case $step in
    export)  run heavy "$PY" -W ignore "$MODEL_EXT/localnn_tools.py" export --seg "$RUNROOT/fit/final" --base "$BASE" \
               --d2 "$D2" --point-table "$PT" --design "$GRID_DIR/design" "${FTO[@]}" --out "$S" > "$P/export.log" 2>&1; rc=$?
             tail -n 2 "$P/export.log" | cut -c1-600;;
    contrib) run heavy "$PY" -W ignore "$MODEL_EXT/localnn_tools.py" contrib --state "$S" --contrib-code "$CONTRIB" \
               --base "$BASE" --geometry "$GEOM" --out "$RUNROOT/contributions$SUF" --procs "$PROCS" > "$P/contrib.log" 2>&1; rc=$?
             tail -n 4 "$P/contrib.log" | cut -c1-900;;
    maps)    run heavy "$PY" -W ignore "$MODEL_EXT/localnn_tools.py" maps --state "$S" --contrib-code "$CONTRIB" \
               --base "$BASE" --geometry "$GEOM" --out "$RUNROOT/maps$SUF" --procs "$PROCS" > "$P/maps.log" 2>&1; rc=$?
             tail -n 4 "$P/maps.log" | cut -c1-900;;
    *) stamp "unknown step $step"; rc=2;;
  esac
  [ $rc -eq 0 ] || { stamp "step $step FAILED (rc $rc)"; break; }
done
stamp "done rc $rc"
cp -p "$COPS_OUTPUT_ROOT/jobs/logs/cops_localnn_tools-$SLURM_JOB_ID.out" "$P/job.log" 2>/dev/null
exit $rc
