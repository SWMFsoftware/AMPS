#!/bin/bash
# ==== version-controlled copy: header added 2026-11-21, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/jobs/localnn_slope_test.sbatch (on Voyager: /data001/highresolu_vt/cops_output/jobs/localnn_slope_test.sbatch)
# SHA256 of the run copy = SHA256 of this file without these header lines: 6007ae5c45d7a69ee77760ecf4812c6e7774e96b2c5c85ba5a6d4f45c88d600b
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: one slope-effect control fit (decision 3c; not run yet)
# RUN: through jobs/localnn_slope_test_submit.sh
# ==== end of header ====
# Slope-effect control fit (user 2026-11-21, decision 3c): ONE fit of option C f10, P_m 0.01, templates off, SHAP8 12M,
# new H2O law, with the settings, start and code path of that run's cycle-1 fit (stage2r_q002_f10_pm0.01_tploff
# _20261118T162407Z/c1/fit), except the slope multipliers:
#   VARIANT none                 all six held at 1
#           only_40_50 ... only_90_180  that class free (>= 1), the other five held at 1
#           free_min             all six free, lower bound SLOPE_MIN (< 1; required, no default: the user's decision)
# CPU, FP64 arithmetic on the design matrices (built from FP32 records, as the reference), max_nfev 4000. A segment
# stops at WALL and exits 75; the second job of jobs/localnn_slope_test_submit.sh continues it. After the fit: the
# direct facet-sum check. Submit with that script (--export: RUNROOT, VARIANT, MODEL_EXT, SLOPE_MIN).
#SBATCH --job-name=cops_slope3c
#SBATCH --partition=normal128cores
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --output=/data001/highresolu_vt/cops_output/jobs/logs/cops_slope3c-%j.out
set -uo pipefail
R=/data001/highresolu_vt/cops_output
A=/data001/highresolu_vt/AMPS
export AMPS=$A APP=$A/srcCG/LiouvilleCOPS_cuda COPS_OUTPUT_ROOT=$R
export OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 MKL_NUM_THREADS=16 PYTHONUNBUFFERED=1 PYTHONDONTWRITEBYTECODE=1
export NPY_DISABLE_CPU_FEATURES="AVX512F AVX512CD AVX512_KNL AVX512_KNM AVX512_SKX AVX512_CLX AVX512_CNL AVX512_ICL"
export LOCALNN_NONNEG_CODE=$R/model_ext_frozen/nonneg_20261118T072729Z
VARIANT=${VARIANT:?VARIANT is required}
RUNROOT=${RUNROOT:?RUNROOT is required}
MODEL_EXT=${MODEL_EXT:?MODEL_EXT (the frozen code) is required}
export MODEL_EXT
WALL=${WALL:-41955}
ALL="slope_40_50 slope_50_60 slope_60_70 slope_70_80 slope_80_90 slope_90_180"
case $VARIANT in
  none) VA=(--slope-test --slope-mult-min 1 --hold-mult "${ALL// /,}");;
  only_40_50|only_50_60|only_60_70|only_70_80|only_80_90|only_90_180)
    H=; for b in $ALL; do [ "$b" = "slope_${VARIANT#only_}" ] || H=${H:+$H,}$b; done
    VA=(--slope-test --slope-mult-min 1 --hold-mult "$H");;
  free_min)
    : "${SLOPE_MIN:?SLOPE_MIN is required for free_min (the user decides it)}"
    VA=(--slope-test --slope-mult-min "$SLOPE_MIN");;
  *) echo "unknown VARIANT $VARIANT" >&2; exit 2;;
esac
N=$R/runs/shap8_012m_lategrid_dz1.5
L=$N/nonneg/localnn
REF=$L/stage2r_q002_f10_pm0.01_tploff_20261118T162407Z/c1/fit
DES=$L/grid_150m/design_h2oqmin0.02_d2q002
D2=$N/clip_refit/D2_q002_warm/seg01
PT=$N/late_grid_points.txt
CONTRIB=$R/model_ext_frozen/nonneg_contrib_20261118T072729Z
PY=$APP/.venv/bin/python
O=$RUNROOT/$VARIANT
P=$O/provenance-$SLURM_JOB_ID
mkdir -p "$O" && mkdir "$P" || exit 1
(cd "$MODEL_EXT" && sha256sum -c --quiet SHA256SUMS) || { echo "the frozen code does not match its SHA256SUMS" >&2; exit 2; }
cp -p "$0" "$P/localnn_slope_test.sbatch.used" && cp -p "$MODEL_EXT/SHA256SUMS" "$P/MODEL_EXT_SHA256SUMS" || exit 1
{ git -C "$A" rev-parse HEAD; git -C "$A" status --short; } > "$P/git.txt" 2>&1
{ type module >/dev/null 2>&1 && module list; } > "$P/module_list.txt" 2>&1 || echo "no module command" >> "$P/module_list.txt"
echo "CPU job on $(hostname): no GPU used" > "$P/nvidia_smi_q.txt"
{ gcc --version | head -1; echo "CSPICE, nvcc: not used"; "$PY" --version; lscpu | grep -E 'Model name|^CPU\(s\)'; } > "$P/toolchain.txt" 2>&1
"$PY" -c "import importlib.metadata as m; print('\n'.join(sorted(d.metadata['Name']+' '+d.version for d in m.distributions())))" > "$P/pip_list.txt"
trap 'cp -p "$R/jobs/logs/cops_slope3c-$SLURM_JOB_ID.out" "$P/job.log" 2>/dev/null' EXIT
cd "$APP" || exit 1
{ echo '#!/bin/bash'; export -p | grep -E ' (AMPS|APP|COPS_OUTPUT_ROOT|[A-Z]+_NUM_THREADS|NPY_DISABLE_CPU_FEATURES|LOCALNN_NONNEG_CODE|MODEL_EXT)='; echo "cd $APP"; } > "$P/commands.sh"
run() { printf '%q ' "$@" >> "$P/commands.sh"; echo >> "$P/commands.sh"; bash -c 'echo 1000 > /proc/self/oom_score_adj; exec "$@"' heavy "$@"; }
echo "$(date -u +%FT%TZ) job $SLURM_JOB_ID on $(hostname): $VARIANT ${VA[*]} -> $O"
rc=0
if [ ! -e "$O/fit/final" ]; then
  JC=(); ls -d "$O"/fit/seg[0-9][0-9] >/dev/null 2>&1 || JC=(--check-jacobian)
  run "$PY" -W ignore "$MODEL_EXT/localnn_fit.py" fit --design "$DES" --base "$N/base_q002" --d2 "$D2" \
    --point-table "$PT" --grid "$L/grid_150m" --templates off --start-file "$REF/start" \
    --smooth-lambda 8.249495752351283e-05 --smooth-rule-file "$L/lambda_rule_anchored_20261114.json" \
    --mult-penalty-f 0.01 --max-nfev 4000 --wall-time "$WALL" --checkpoint-every 250 "${JC[@]}" "${VA[@]}" \
    --out "$O/fit" > "$P/fit.log" 2>&1
  rc=$?
  tail -n 4 "$P/fit.log" | cut -c1-600
fi
if [ $rc -eq 0 ] && [ -e "$O/fit/final" ]; then
  S=$(readlink "$O/fit/final")
  if [ ! -e "$O/fit/direct_check_$S" ]; then
    export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
    run "$PY" -W ignore "$MODEL_EXT/localnn_direct.py" --seg "$O/fit/$S" --base "$N/base_q002" --d2 "$D2" \
      --point-table "$PT" --design "$DES" --contrib-code "$CONTRIB" --procs "$SLURM_CPUS_PER_TASK" > "$P/direct.log" 2>&1
    rc=$?
    tail -n 2 "$P/direct.log" | cut -c1-600
  fi
fi
echo "$(date -u +%FT%TZ) done rc $rc"
exit $rc
