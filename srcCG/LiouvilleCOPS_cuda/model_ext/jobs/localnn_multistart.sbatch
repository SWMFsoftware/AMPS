#!/bin/bash
# ==== version-controlled copy: header added 2026-11-22, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/jobs/localnn_multistart.sbatch (on Voyager: /data001/highresolu_vt/cops_output/jobs/localnn_multistart.sbatch)
# SHA256 of the run copy = SHA256 of this file without these header lines: 244c6ae9cff32495d467d0bf06ce5d065fca05910142eb4c02b3628087b65be8
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: one multi-start refit of the option C reference problem
# RUN: through jobs/localnn_multistart_submit.sh
# ==== end of header ====
# Multi-start refit of the option C reference problem (user 2026-11-22: "multi-start when ready"): ONE fit with all six
# slope multipliers free and >= 1, started from the end state of one control fit of decision 3c, IF that state is better
# than the reference and has every multiplier >= 1 (localnn_multistart.py gate). Otherwise it writes SKIPPED.txt and
# ends. Every setting and the fitter are the reference's: model_ext_frozen/localnn_q002_20261118T162327Z, f 10, P_m 0.01,
# templates off, max_nfev 4000, CPU. Submit with jobs/localnn_multistart_submit.sh (two jobs per fit: the second
# continues after a wall-time stop, exit 75, or only confirms). Env: RUNROOT, CONTROLS, PARENT, GATE.
#SBATCH --job-name=cops_multistart
#SBATCH --partition=normal128cores
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --output=/data001/highresolu_vt/cops_output/jobs/logs/cops_multistart-%j.out
set -uo pipefail
R=/data001/highresolu_vt/cops_output
A=/data001/highresolu_vt/AMPS
export AMPS=$A APP=$A/srcCG/LiouvilleCOPS_cuda COPS_OUTPUT_ROOT=$R
export OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 MKL_NUM_THREADS=16 PYTHONUNBUFFERED=1 PYTHONDONTWRITEBYTECODE=1
export NPY_DISABLE_CPU_FEATURES="AVX512F AVX512CD AVX512_KNL AVX512_KNM AVX512_SKX AVX512_CLX AVX512_CNL AVX512_ICL"
export LOCALNN_NONNEG_CODE=$R/model_ext_frozen/nonneg_20261118T072729Z
export MODEL_EXT=$R/model_ext_frozen/localnn_q002_20261118T162327Z
RUNROOT=${RUNROOT:?RUNROOT is required}
CONTROLS=${CONTROLS:?CONTROLS is required}
PARENT=${PARENT:?PARENT is required}
GATE=${GATE:?GATE (the frozen localnn_multistart.py) is required}
WALL=${WALL:-41955}
N=$R/runs/shap8_012m_lategrid_dz1.5
L=$N/nonneg/localnn
DES=$L/grid_150m/design_h2oqmin0.02_d2q002
D2=$N/clip_refit/D2_q002_warm/seg01
PT=$N/late_grid_points.txt
CONTRIB=$R/model_ext_frozen/nonneg_contrib_20261118T072729Z
PY=$APP/.venv/bin/python
O=$RUNROOT/from_$PARENT
P=$O/provenance-$SLURM_JOB_ID
mkdir -p "$O" && mkdir "$P" || exit 1
(cd "$MODEL_EXT" && sha256sum -c --quiet SHA256SUMS) || { echo "the frozen fitter does not match its SHA256SUMS" >&2; exit 2; }
cp -p "$0" "$P/localnn_multistart.sbatch.used" && cp -p "$GATE" "$P/" && cp -p "$MODEL_EXT/SHA256SUMS" "$P/MODEL_EXT_SHA256SUMS" || exit 1
{ git -C "$A" rev-parse HEAD; git -C "$A" status --short; } > "$P/git.txt" 2>&1
{ type module >/dev/null 2>&1 && module list; } > "$P/module_list.txt" 2>&1 || echo "no module command" >> "$P/module_list.txt"
echo "CPU job on $(hostname): no GPU used" > "$P/nvidia_smi_q.txt"
{ gcc --version | head -1; echo "CSPICE, nvcc: not used"; "$PY" --version; lscpu | grep -E 'Model name|^CPU\(s\)'; } > "$P/toolchain.txt" 2>&1
"$PY" -c "import importlib.metadata as m; print('\n'.join(sorted(d.metadata['Name']+' '+d.version for d in m.distributions())))" > "$P/pip_list.txt"
trap 'cp -p "$R/jobs/logs/cops_multistart-$SLURM_JOB_ID.out" "$P/job.log" 2>/dev/null' EXIT
cd "$APP" || exit 1
{ echo '#!/bin/bash'; export -p | grep -E ' (AMPS|APP|COPS_OUTPUT_ROOT|[A-Z]+_NUM_THREADS|NPY_DISABLE_CPU_FEATURES|LOCALNN_NONNEG_CODE|MODEL_EXT)='; echo "cd $APP"; } > "$P/commands.sh"
run() { printf '%q ' "$@" >> "$P/commands.sh"; echo >> "$P/commands.sh"; bash -c 'echo 1000 > /proc/self/oom_score_adj; exec "$@"' heavy "$@"; }
echo "$(date -u +%FT%TZ) job $SLURM_JOB_ID on $(hostname): all six free (>= 1) from $PARENT -> $O"

if [ -e "$O/SKIPPED.txt" ]; then echo "skipped earlier: $(cat "$O/SKIPPED.txt")"; exit 0; fi
if [ ! -e "$O/start" ]; then
  why=$("$PY" -I "$GATE" gate --controls "$CONTROLS" --parent "$PARENT"); g=$?
  echo "$why"
  if [ $g -eq 10 ]; then echo "$why" > "$O/SKIPPED.txt"; exit 0; fi
  [ $g -eq 0 ] || exit $g
  echo "$why" > "$O/GATE.txt"
  SEG=$(readlink -f "$CONTROLS/$PARENT/fit/final")
  run "$PY" -W ignore "$MODEL_EXT/localnn_fit.py" start --design "$DES" --base "$N/base_q002" --d2 "$D2" \
    --point-table "$PT" --grid "$L/grid_150m" --templates off --start state --state-file "$SEG/localnn_state.npz" \
    --out "$O/start" > "$P/start.log" 2>&1 || { tail -n 5 "$P/start.log"; exit 4; }
fi
"$PY" -I "$GATE" check-start --state "$(readlink -f "$CONTROLS/$PARENT/fit/final")/localnn_state.npz" --start "$O/start" || exit 4

rc=0
if [ ! -e "$O/fit/final" ]; then
  JC=(); ls -d "$O"/fit/seg[0-9][0-9] >/dev/null 2>&1 || JC=(--check-jacobian)
  run "$PY" -W ignore "$MODEL_EXT/localnn_fit.py" fit --design "$DES" --base "$N/base_q002" --d2 "$D2" \
    --point-table "$PT" --grid "$L/grid_150m" --templates off --slope-mult-min 1 --start-file "$O/start" \
    --smooth-lambda 8.249495752351283e-05 --smooth-rule-file "$L/lambda_rule_anchored_20261114.json" \
    --mult-penalty-f 0.01 --max-nfev 4000 --wall-time "$WALL" --checkpoint-every 250 "${JC[@]}" \
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
