#!/bin/bash
# ==== version-controlled copy: header added 2026-11-22, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/jobs/localnn_cliff_fit.sbatch (on Voyager: /data001/highresolu_vt/cops_output/jobs/localnn_cliff_fit.sbatch)
# SHA256 of the run copy = SHA256 of this file without these header lines: a00e8e5ed668b5630af170bedc4c49cf7a177531387a18a2d7bb995a04862bca
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: one fit of the cliff test (not run yet)
# RUN: through jobs/localnn_cliff_submit.sh
# ==== end of header ====
# One fit of the cliff test (user 2026-11-22: "prepare the cliff test"; plan figures_cluster/12M/
# q002_decisions_20261121/PLAN_cliff.md). Option C f10, templates off, ALL multipliers held; P_m is then a constant: f_m 0.
#   FIT=cliffs_off_from_{held,d2}    60-70 ... 90-180 deg held at HELD_VALUE (0.01), 40-60 deg at 1
#   FIT=moderate_off_from_{held,d2}  the control: 40-50, 50-60 deg held at HELD_VALUE, the cliff classes at 1
#   from held: the all-held baseline (follow-up c's chain_pick.json, or BASELINE); from d2: the all-held fit's own start.
# max_nfev 4000, one restart after it (stage2), wall-time stop exit 75. Plumbing = jobs/localnn_ac_fit.sbatch.
# Env: RUNROOT, FIT, TOOLS, CONTROLS_C (the run root of follow-up c); test only: NFEV, BASELINE.
#SBATCH --job-name=cops_cliff
#SBATCH --partition=normal128cores
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --output=/data001/highresolu_vt/cops_output/jobs/logs/cops_cliff-%j.out
set -uo pipefail
R=/data001/highresolu_vt/cops_output
A=/data001/highresolu_vt/AMPS
export AMPS=$A APP=$A/srcCG/LiouvilleCOPS_cuda COPS_OUTPUT_ROOT=$R
export OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 MKL_NUM_THREADS=16 PYTHONUNBUFFERED=1 PYTHONDONTWRITEBYTECODE=1
export NPY_DISABLE_CPU_FEATURES="AVX512F AVX512CD AVX512_KNL AVX512_KNM AVX512_SKX AVX512_CLX AVX512_CNL AVX512_ICL"
export LOCALNN_NONNEG_CODE=$R/model_ext_frozen/nonneg_20261118T072729Z
export MODEL_EXT=${TOOLS:?TOOLS is required}/code
RUNROOT=${RUNROOT:?RUNROOT is required}
FIT=${FIT:?FIT is required}
TOOLS=${TOOLS:?TOOLS is required}
NFEV=${NFEV:-4000}
N=$R/runs/shap8_012m_lategrid_dz1.5
L=$N/nonneg/localnn
COM=(--design "$L/grid_150m/design_h2oqmin0.02_d2q002" --base "$N/base_q002" --d2 "$N/clip_refit/D2_q002_warm/seg01" --point-table "$N/late_grid_points.txt")
REF=$L/stage2r_q002_f10_pm0.01_tploff_20261118T162407Z/c1/fit/final
MULTI=$L/multistart_q002_f10_pm0.01_tploff_20261122T065041Z
CONTRIB=$R/model_ext_frozen/nonneg_contrib_20261118T072729Z
PY=$APP/.venv/bin/python
V=${HELD_VALUE:-0.01}
case $FIT in
  cliffs_off_from_*)   H=slope_40_50,slope_50_60,slope_60_70=$V,slope_70_80=$V,slope_80_90=$V,slope_90_180=$V;;
  moderate_off_from_*) H=slope_40_50=$V,slope_50_60=$V,slope_60_70,slope_70_80,slope_80_90,slope_90_180;;
  *) echo "unknown FIT $FIT" >&2; exit 2;;
esac
MODE=(--slope-test --slope-mult-min 1 --hold-mult "$H")
D2START=
case $FIT in
  *_from_d2) D2START=$(dirname "$REF")/start; FROM=$REF;;
  *_from_held)
    if [ -z "${BASELINE:-}" ]; then
      BASELINE=$("$PY" "$TOOLS/localnn_cliff.py" baseline --controls-c "${CONTROLS_C:?CONTROLS_C is required}") || exit 3
    fi
    FROM=$BASELINE/final;;
  *) echo "unknown FIT $FIT" >&2; exit 2;;
esac
PY=$APP/.venv/bin/python
O=$RUNROOT/$FIT
P=$O/provenance-$SLURM_JOB_ID
mkdir -p "$O" && mkdir "$P" || exit 1
(cd "$MODEL_EXT" && sha256sum -c --quiet SHA256SUMS) && (cd "$TOOLS" && sha256sum -c --quiet SHA256SUMS) || { echo "frozen code or starts changed" >&2; exit 2; }
cp -p "$0" "$P/localnn_cliff_fit.sbatch.used" && cp -p "$MODEL_EXT/SHA256SUMS" "$P/MODEL_EXT_SHA256SUMS" && cp -p "$TOOLS/SHA256SUMS" "$P/TOOLS_SHA256SUMS" || exit 1
{ git -C "$A" rev-parse HEAD; git -C "$A" status --short; } > "$P/git.txt" 2>&1
{ type module >/dev/null 2>&1 && module list; } > "$P/module_list.txt" 2>&1 || echo "no module command" >> "$P/module_list.txt"
echo "CPU job on $(hostname): no GPU used" > "$P/nvidia_smi_q.txt"
{ gcc --version | head -1; echo "CSPICE, nvcc: not used"; "$PY" --version; lscpu | grep -E 'Model name|^CPU\(s\)'; } > "$P/toolchain.txt" 2>&1
"$PY" -c "import importlib.metadata as m; print('\n'.join(sorted(d.metadata['Name']+' '+d.version for d in m.distributions())))" > "$P/pip_list.txt"
trap 'cp -p "$R/jobs/logs/cops_cliff-$SLURM_JOB_ID.out" "$P/job.log" 2>/dev/null' EXIT
cd "$APP" || exit 1
{ echo '#!/bin/bash'; export -p | grep -E ' (AMPS|APP|COPS_OUTPUT_ROOT|[A-Z]+_NUM_THREADS|NPY_DISABLE_CPU_FEATURES|LOCALNN_NONNEG_CODE|MODEL_EXT)='; echo "cd $APP"; } > "$P/commands.sh"
run() { printf '%q ' "$@" >> "$P/commands.sh"; echo >> "$P/commands.sh"; bash -c 'echo 1000 > /proc/self/oom_score_adj; exec "$@"' heavy "$@"; }
secs() { local t=$1 d=0 a b c; [[ $t == *-* ]] && { d=${t%%-*}; t=${t#*-}; }; IFS=: read -r a b c <<< "$t"; [ -z "$c" ] && { c=$b; b=$a; a=0; }; echo $(( d*86400 + 10#$a*3600 + 10#$b*60 + 10#$c )); }
DEADLINE=$(( $(date +%s) + $(secs "$(squeue -h -j "$SLURM_JOB_ID" -o %l)") - 1260 ))
echo "$(date -u +%FT%TZ) job $SLURM_JOB_ID on $(hostname): $FIT ${MODE[*]} NFEV $NFEV -> $O"

for k in 1 2; do
  S=$O/stage$k
  if [ ! -e "$S/fit/final" ]; then
    if [ ! -e "$S/start" ]; then
      mkdir -p "$S"
    fi
    if [ ! -e "$S/start" ] && [ $k -eq 1 ] && [ -n "$D2START" ]; then ln -s "$D2START" "$S/start"; fi
    if [ ! -e "$S/start" ]; then
      run "$PY" -W ignore "$MODEL_EXT/localnn_fit.py" start "${COM[@]}" --grid "$L/grid_150m" --templates off --start state \
        --state-file "$(readlink -f "$FROM")/localnn_state.npz" --out "$S/start" > "$P/stage${k}_start.log" 2>&1 \
        || { tail -n 5 "$P/stage${k}_start.log"; exit 4; }
    fi
    WALL=$(( DEADLINE - $(date +%s) ))
    [ $WALL -gt 300 ] || { echo "no time left for stage $k: exit 75"; exit 75; }
    JC=(); ls -d "$S"/fit/seg[0-9][0-9] >/dev/null 2>&1 || JC=(--check-jacobian)
    run "$PY" -W ignore "$MODEL_EXT/localnn_fit.py" fit "${COM[@]}" --grid "$L/grid_150m" --templates off "${MODE[@]}" \
      --start-file "$S/start" --smooth-lambda 8.249495752351283e-05 --smooth-rule-file "$L/lambda_rule_anchored_20261114.json" \
      --mult-penalty-f 0 --max-nfev "$NFEV" --wall-time "$WALL" --checkpoint-every 250 "${JC[@]}" \
      --out "$S/fit" > "$P/stage${k}_fit.log" 2>&1
    rc=$?
    tail -n 3 "$P/stage${k}_fit.log" | cut -c1-500
    [ $rc -eq 0 ] || exit $rc
  fi
  SEG=$(readlink "$S/fit/final")
  if [ ! -e "$S/fit/direct_check_$SEG" ]; then
    OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 run "$PY" -W ignore "$MODEL_EXT/localnn_direct.py" \
      --seg "$S/fit/$SEG" "${COM[@]}" --contrib-code "$CONTRIB" --procs "$SLURM_CPUS_PER_TASK" > "$P/stage${k}_direct.log" 2>&1
    rc=$?
    tail -n 1 "$P/stage${k}_direct.log" | cut -c1-400
    [ $rc -eq 0 ] || exit $rc
  fi
  "$PY" "$TOOLS/localnn_ac.py" after-stage "$S/fit"; rc=$?
  if [ $rc -eq 0 ]; then echo "$(date -u +%FT%TZ) $FIT finished in stage $k"; exit 0; fi
  [ $rc -eq 10 ] || exit $rc
  echo "$(date -u +%FT%TZ) $FIT stage $k stopped at max_nfev"
  FROM=$S/fit/final
done
echo "$(date -u +%FT%TZ) $FIT: not converged after the restart (reported, not continued)"
exit 0
