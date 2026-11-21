#!/bin/bash
# ==== version-controlled copy: header added 2026-11-18, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/tests/nightfix_diagnostic_20261118/jobs/nightfix_diag.sbatch (on Voyager: /data001/highresolu_vt/cops_output/tests/nightfix_diagnostic_20261118/jobs/nightfix_diag.sbatch)
# SHA256 of the run copy = SHA256 of this file without these header lines: 5c71964c530c822fa70ab728260ce71a91c918133cd7ca39e20d0bdb272b434f
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: job (2026-11-18). One state of the night-side fix diagnostic (nightfix_diag.py) on a CPU node, or BASE_RUN=1.
# RUN: jobs/nightfix_diag_submit_all.sh (v1 + base check), jobs/nightfix_diag_submit_all_v2.sh (v2).
# ==== end of header ====
#SBATCH --job-name=nightfix_diag
#SBATCH --partition=normal128cores
#SBATCH --account=general
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=96G
#SBATCH --time=02:00:00
#SBATCH --output=/data001/highresolu_vt/cops_output/tests/nightfix_diagnostic_20261118/logs/%x-%j.out
# H2O night-side fix diagnostic (lead task 4, 2026-11-18): one final state re-evaluated with A held fixed under the old
# (0.2) and the new (0.02) H2O q_min/q_max (nightfix_diag.py; no fitting). Settings (sbatch --export=ALL,VAR=...):
#   KEY      the state's key (also the output name <OUT_DIR>/<KEY>.json)
#   FT, FIT  the --ft-dir / --fit-dir the contribution tools used for this state
#   GEOM     the state's validated nonneg_maps geometry
#   ADAPTER  optional state adapter (frozen localnn_tools.py for option C, d2_tools.py for D2)
#   RECORDED optional recorded data cost to reproduce with the old law
#   NPYDIS   1: NPY_DISABLE_CPU_FEATURES as the option-C tool jobs (A re-derived bitwise)
#   BASE_RUN 1: instead of a state, fit_cops' base consistency check for both laws
#   CODE     the diagnostic code dir (default: the frozen snapshot of nightfix_diag.py); CONTRIB the contribution code
set -u
export AMPS=/data001/highresolu_vt/AMPS APP=/data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda
export COPS_OUTPUT_ROOT=/data001/highresolu_vt/cops_output
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1
UNSET=$(env | grep -o -E '^(COPS_POS_[A-Z_]*|COPS_SLOPEMUL_[A-Z_]*|COPS_NONNEG_[A-Z_]*)' | sort | paste -sd' ' -)
[ -n "$UNSET" ] && unset $UNSET
D=$COPS_OUTPUT_ROOT/tests/nightfix_diagnostic_20261118
R=$COPS_OUTPUT_ROOT/runs/shap8_012m_lategrid_dz1.5
BASE=${BASE:-$R/base}
CODE=${CODE:?CODE (frozen nightfix_diag snapshot) required}
CONTRIB=${CONTRIB:-$COPS_OUTPUT_ROOT/model_ext_frozen/nonneg_contrib_20261118T072729Z}
OUT_DIR=${OUT_DIR:-$D/results}
PY=$APP/.venv/bin/python
[ "${NPYDIS:-0}" = 1 ] && export NPY_DISABLE_CPU_FEATURES="AVX512F AVX512CD AVX512_KNL AVX512_KNM AVX512_SKX AVX512_CLX AVX512_CNL AVX512_ICL"
P=$D/provenance/job-$SLURM_JOB_ID; mkdir -p "$P" "$OUT_DIR" || exit 1
cp -p "$0" "$P/nightfix_diag.sbatch.used"
cp -p "$CODE/nightfix_diag.py" "$P/"
{ git -C "$AMPS" rev-parse HEAD; git -C "$AMPS" status --short; } > "$P/git.txt" 2>&1
module list > "$P/module_list.txt" 2>&1
"$PY" -m pip list > "$P/pip_list.txt" 2>&1
sha256sum "$CODE/nightfix_diag.py" "$CONTRIB"/*.py ${ADAPTER:+"$ADAPTER"} "$AMPS/srcCG/LiouvilleCOPS/tools/fit_cops.py" \
    "$APP/cops_cuda/evaluator.py" > "$P/code_sha256.txt"
free -g > "$P/free.txt"; hostname > "$P/host.txt"
cd "$APP" || exit 1
if [ "${BASE_RUN:-0}" = 1 ]; then
  CMD=("$PY" -W ignore "$CODE/nightfix_diag.py" --base-run --base "$BASE" --contrib-code "$CONTRIB"
       --out "$OUT_DIR/base_run_consistency.json" --procs "$SLURM_CPUS_PER_TASK")
else
  CMD=("$PY" -W ignore "$CODE/nightfix_diag.py" --key "$KEY" --ft-dir "$FT" --fit-dir "$FIT" --base "$BASE"
       --geometry "$GEOM" --contrib-code "$CONTRIB" ${ADAPTER:+--adapter "$ADAPTER"}
       ${RECORDED:+--recorded-cost "$RECORDED"} --out "$OUT_DIR/$KEY.json" --procs "$SLURM_CPUS_PER_TASK")
fi
{ echo "cd $APP"; [ -n "${NPY_DISABLE_CPU_FEATURES:-}" ] && echo "export NPY_DISABLE_CPU_FEATURES=\"$NPY_DISABLE_CPU_FEATURES\"";
  echo "export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1"; printf '%q ' "${CMD[@]}"; echo; } > "$P/commands.sh"
echo "$(date -u +%FT%TZ) job $SLURM_JOB_ID on $(hostname): ${CMD[*]}"
bash -c 'echo 1000 > /proc/self/oom_score_adj; exec "$@"' heavy "${CMD[@]}"
rc=$?
echo "$(date -u +%FT%TZ) rc=$rc"
exit $rc
