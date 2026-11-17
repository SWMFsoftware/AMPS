#!/bin/bash
# ==== version-controlled copy: header added 2026-11-15, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/jobs/localnn_report.sbatch (on Voyager: /data001/highresolu_vt/cops_output/jobs/localnn_report.sbatch)
# SHA256 of the run copy = SHA256 of this file without these header lines: d19388b2b24e030557adae0e5d7e362548d4f53ef0f170c7ae38fb2bb704cafc
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: job script (CPU, option C): localnn_report.py for one run root (DOCS_ONLY=1: only README.md / HOW-TO-RERUN.md).
# RUN: test 125019; stage-1 125054/125056, 125058/125060, 125062/125064 (chained).
# 2026-11-15: 2026-11-17: NPY_DISABLE_CPU_FEATURES (as the chain driver; reproducible on AVX-512 nodes).
# RUN 2026-11-15: the chain posts.
# ==== end of header ====
# Option C (local non-negative basis, lead 2026-11-14 item 5): summary numbers of ONE finished fit, the bins table in
# fit/final and the run root's README.md / HOW-TO-RERUN.md (model_ext/localnn_report.py). CPU (normal128cores): one pass
# over the base's kernel records for the active-area shares, the validated maps geometry (label t11 = 10:33:10 UTC)
# for Q_H2O and the production areas.
# Settings: RUNROOT (required: a run root with fit/final), MODEL_EXT (the frozen localnn snapshot), GRID_DIR (default
#   grid_150m), GEOM (default maps_geometry/job_124775/geometry), CONTRIB (frozen nonneg_contrib), PROCS,
#   DOCS_ONLY=1 (only README.md / HOW-TO-RERUN.md from analysis/report.json, e.g. after the GPU evaluation).
#SBATCH --job-name=cops_localnn_report
#SBATCH --partition=normal128cores
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=200G
#SBATCH --time=01:00:00
#SBATCH --output=/data001/highresolu_vt/cops_output/jobs/logs/cops_localnn_report-%j.out
set -uo pipefail
# numpy without AVX-512 kernels, as the chain driver (2026-11-17: on an AVX-512 node (node028, EPYC 9554) np.exp differs
# in the last bit at 38 of the 271 held facets and the tools' bitwise reproduction of A failed, post 136094)
export NPY_DISABLE_CPU_FEATURES="AVX512F AVX512CD AVX512_KNL AVX512_KNM AVX512_SKX AVX512_CLX AVX512_CNL AVX512_ICL"
export AMPS=/data001/highresolu_vt/AMPS APP=/data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda
export COPS_OUTPUT_ROOT=/data001/highresolu_vt/cops_output
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8 PYTHONUNBUFFERED=1
UNSET=$(env | grep -o -E '^(COPS_POS_[A-Z_]*|COPS_SLOPEMUL_[A-Z_]*|COPS_NONNEG_[A-Z_]*)' | sort | paste -sd' ' -)
[ -n "$UNSET" ] && unset $UNSET
R=$COPS_OUTPUT_ROOT/runs/shap8_012m_lategrid_dz1.5
BASE=$R/base; D2=$R/clip_refit/D2/final; PT=$R/late_grid_points.txt
C9=$R/nonneg/softplus/chain_5cycles_smm1_20261113T041215Z/c9
export MODEL_EXT=$(readlink -f "${MODEL_EXT:-$COPS_OUTPUT_ROOT/model_ext}")
GRID_DIR=${GRID_DIR:-$R/nonneg/localnn/grid_150m}
GEOM=${GEOM:-$R/maps_geometry/job_124775/geometry}
CONTRIB=${CONTRIB:-$COPS_OUTPUT_ROOT/model_ext_frozen/nonneg_contrib_20261112T210012Z}
RUNROOT=$(readlink -f "${RUNROOT:?RUNROOT is required}")
PROCS=${PROCS:-${SLURM_CPUS_PER_TASK:-8}}
PY=$APP/.venv/bin/python
stamp() { echo "$(date -u +%FT%TZ) $*"; }
P=$RUNROOT/provenance-$SLURM_JOB_ID-report; mkdir "$P" || exit 1
cd "$APP" || exit 1
stamp "job $SLURM_JOB_ID on $(hostname): RUNROOT=$RUNROOT MODEL_EXT=$MODEL_EXT"
cp -p "$0" "$P/localnn_report.sbatch.used"
cp -p "$MODEL_EXT"/localnn_*.py "$MODEL_EXT/positive_problem.py" "$P/"
( cd "$P" && sha256sum *.py localnn_report.sbatch.used > code_sha256.txt )
{ git -C "$AMPS" rev-parse HEAD; git -C "$AMPS" status --short; } > "$P/git.txt" 2>&1
{ type module >/dev/null 2>&1 && module list; } > "$P/module_list.txt" 2>&1 || echo "module command not available (no modules loaded)" >> "$P/module_list.txt"
echo "CPU job on $(hostname): no GPU used" > "$P/nvidia_smi_q.txt"
{ gcc --version | head -1; lscpu | grep -E 'Model name|^CPU\(s\)'; } > "$P/toolchain.txt" 2>&1
~/.local/bin/uv pip list --python "$PY" > "$P/pip_list.txt" 2>&1 || "$PY" -m pip list > "$P/pip_list.txt" 2>&1
{ echo '#!/bin/bash'; grep '^export ' "$0"; echo "export MODEL_EXT=$MODEL_EXT"; echo "cd $APP"; } > "$P/commands.sh"
CMD=("$PY" -W ignore "$MODEL_EXT/localnn_report.py" --runroot "$RUNROOT" --base "$BASE" --d2 "$D2" --point-table "$PT" \
     --design "$GRID_DIR/design" --grid "$GRID_DIR" --geometry "$GEOM" --contrib-code "$CONTRIB" --c9-ft "$C9/ft/final" \
     --c9-fit "$C9/fit/final" --procs "$PROCS")
[ "${DOCS_ONLY:-0}" = 1 ] && CMD+=(--docs-only)
printf '%q ' "${CMD[@]}" >> "$P/commands.sh"; echo >> "$P/commands.sh"
bash -c 'echo 1000 > /proc/self/oom_score_adj; exec "$@"' _ "${CMD[@]}" > "$P/report.log" 2>&1; rc=$?
tail -n 4 "$P/report.log" | cut -c1-1500
stamp "done rc $rc"
cp -p "$COPS_OUTPUT_ROOT/jobs/logs/cops_localnn_report-$SLURM_JOB_ID.out" "$P/job.log" 2>/dev/null
exit $rc
