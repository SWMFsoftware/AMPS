#!/bin/bash
# ==== version-controlled copy: header added 2026-11-13, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/jobs/d2_tools.sbatch (on Voyager: /data001/highresolu_vt/cops_output/jobs/d2_tools.sbatch)
# SHA256 of the run copy = SHA256 of this file without these header lines: 1b47f4674316a7e98b29933b0ec085074b213a29f8f1ce2e8ee4efa25f3d98c3
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: D2 export / contributions / verify / maps (2026-11-17)
# RUN: see model_ext/README.md (full sets of the stage-2 states, 2026-11-17)
# ==== end of header ====
# The clipped refit D2: the maps and contribution tools of the nonneg chains on D2's final state (model_ext/d2_tools.py:
# export folder + the frozen nonneg_maps.py / nonneg_contrib.py, unchanged, with state_activity replaced by the D2 one;
# 2026-11-17, full figure set of D2; pattern of localnn_tools.sbatch). CPU (normal128cores).
#   1. export  -> $D2DIR/tools_state (use as --fit-dir and --ft-dir; NOT a fine tuning, see its README.txt)
#   2. contrib -> $D2DIR/contributions (P% source areas, occurrence, slope-bin fractions; verification vs the tables)
#   3. verify  -> $D2DIR/tools_state_verify.json (forward model of the export vs D2's per-row model and data cost)
#   4. maps    -> $D2DIR/maps (activity and source-rate maps, .vtu + Tecplot, 13 times)
# Settings: STEPS (default "export contrib verify maps"), MODEL_EXT (frozen d2_tools snapshot), GEOM (default
#   maps_geometry/job_124775/geometry), CONTRIB (frozen nonneg_contrib folder), BINS_FROM (nonneg_facets.npz of a nonneg
#   fit on the same base: slope bins; default the c9 fit), PROCS.
#SBATCH --job-name=cops_d2_tools
#SBATCH --partition=normal128cores
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=320G
#SBATCH --time=02:00:00
#SBATCH --output=/data001/highresolu_vt/cops_output/jobs/logs/cops_d2_tools-%j.out
set -uo pipefail
# numpy without AVX-512 kernels (as the chains and the other tools' jobs)
export NPY_DISABLE_CPU_FEATURES="AVX512F AVX512CD AVX512_KNL AVX512_KNM AVX512_SKX AVX512_CLX AVX512_CNL AVX512_ICL"
export AMPS=/data001/highresolu_vt/AMPS APP=/data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda
export COPS_OUTPUT_ROOT=/data001/highresolu_vt/cops_output
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1
R=$COPS_OUTPUT_ROOT/runs/shap8_012m_lategrid_dz1.5
BASE=$R/base; D2DIR=$R/clip_refit/D2; D2=$D2DIR/final
export MODEL_EXT=$(readlink -f "${MODEL_EXT:-$COPS_OUTPUT_ROOT/model_ext}")
GEOM=${GEOM:-$R/maps_geometry/job_124775/geometry}
CONTRIB=${CONTRIB:-$COPS_OUTPUT_ROOT/model_ext_frozen/nonneg_contrib_20261112T210012Z}
BINS_FROM=${BINS_FROM:-$R/nonneg/softplus/chain_5cycles_smm1_20261113T041215Z/c9/fit/final/nonneg_facets.npz}
STEPS=${STEPS:-export contrib verify maps}
PROCS=${PROCS:-${SLURM_CPUS_PER_TASK:-8}}
PY=$APP/.venv/bin/python
stamp() { echo "$(date -u +%FT%TZ) $*"; }
heavy() { bash -c 'echo 1000 > /proc/self/oom_score_adj; exec "$@"' heavy "$@"; }
P=$D2DIR/provenance-$SLURM_JOB_ID-d2tools; mkdir "$P" || exit 1
cd "$APP" || exit 1
stamp "job $SLURM_JOB_ID on $(hostname): STEPS='$STEPS' MODEL_EXT=$MODEL_EXT CONTRIB=$CONTRIB GEOM=$GEOM BINS_FROM=$BINS_FROM"
cp -p "$0" "$P/d2_tools.sbatch.used"
cp -p "$MODEL_EXT/d2_tools.py" "$P/"; mkdir "$P/contrib_code" && cp -p "$CONTRIB"/* "$P/contrib_code/"
( cd "$P" && sha256sum d2_tools.py contrib_code/*.py d2_tools.sbatch.used > code_sha256.txt )
{ git -C "$AMPS" rev-parse HEAD; git -C "$AMPS" status --short; } > "$P/git.txt" 2>&1
{ type module >/dev/null 2>&1 && module list; } > "$P/module_list.txt" 2>&1 || echo "module command not available (no modules loaded)" >> "$P/module_list.txt"
echo "CPU job on $(hostname): no GPU used" > "$P/nvidia_smi_q.txt"
{ gcc --version | head -1; lscpu | grep -E 'Model name|^CPU\(s\)'; } > "$P/toolchain.txt" 2>&1
~/.local/bin/uv pip list --python "$PY" > "$P/pip_list.txt" 2>&1 || "$PY" -m pip list > "$P/pip_list.txt" 2>&1
df -h "$COPS_OUTPUT_ROOT" > "$P/df_start.txt" 2>&1
{ echo '#!/bin/bash'; grep '^export ' "$0"; echo "export MODEL_EXT=$MODEL_EXT"; echo "cd $APP"; } > "$P/commands.sh"
run() { printf '%q ' "$@" >> "$P/commands.sh"; echo >> "$P/commands.sh"; "$@"; }
S=$D2DIR/tools_state
rc=0
for step in $STEPS; do
  case $step in
    export)  run heavy "$PY" -W ignore "$MODEL_EXT/d2_tools.py" export --d2 "$D2" --base "$BASE" --bins-from "$BINS_FROM" \
               --out "$S" > "$P/export.log" 2>&1; rc=$?
             tail -n 3 "$P/export.log" | cut -c1-900;;
    contrib) run heavy "$PY" -W ignore "$MODEL_EXT/d2_tools.py" contrib --state "$S" --contrib-code "$CONTRIB" \
               --base "$BASE" --geometry "$GEOM" --out "$D2DIR/contributions" --procs "$PROCS" > "$P/contrib.log" 2>&1; rc=$?
             tail -n 6 "$P/contrib.log" | cut -c1-900;;
    verify)  run heavy "$PY" -W ignore "$MODEL_EXT/d2_tools.py" verify --state "$S" --contrib "$D2DIR/contributions" \
               --out "$D2DIR/tools_state_verify.json" > "$P/verify.log" 2>&1; rc=$?
             tail -n 2 "$P/verify.log" | cut -c1-1500;;
    maps)    run heavy "$PY" -W ignore "$MODEL_EXT/d2_tools.py" maps --state "$S" --contrib-code "$CONTRIB" \
               --base "$BASE" --geometry "$GEOM" --out "$D2DIR/maps" --procs "$PROCS" > "$P/maps.log" 2>&1; rc=$?
             tail -n 4 "$P/maps.log" | cut -c1-900;;
    *) stamp "unknown step $step"; rc=2;;
  esac
  [ $rc -eq 0 ] || { stamp "step $step FAILED (rc $rc)"; break; }
done
df -h "$COPS_OUTPUT_ROOT" > "$P/df_end.txt" 2>&1
stamp "done rc $rc"
cp -p "$COPS_OUTPUT_ROOT/jobs/logs/cops_d2_tools-$SLURM_JOB_ID.out" "$P/job.log" 2>/dev/null
exit $rc
