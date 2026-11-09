#!/bin/bash
# ==== version-controlled copy: header added 2026-11-09, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/jobs/clip_refit.sbatch (on Voyager: /data001/highresolu_vt/cops_output/jobs/clip_refit.sbatch)
# SHA256 of the run copy = SHA256 of this file without these header lines: 44e185b641d514aaf73ebf83b68d7acd1051e12ae8d65311e921ca0be485e62a
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: job. Diagnostic refit of the original-T cycle-6 clipped model: VARIANT=D1 | D2 on the lategrid base (+1.5 m),
#   D0 = control on c6's own base with its record store; MODEL_EXT, NFEV 4000, resume segments. 2026-11-08.
# RUN: sbatch --nodelist=gpu001 --export=ALL,VARIANT=D1,MODEL_EXT=<snapshot> clip_refit.sbatch   (D0 122880, D1 122865, D2 122866)
# ==== end of header ====
# Refit of the earlier CLIPPED-LINEAR model (cycle c6, job 122221: runs/shap8_012m_flat_earlier/grain/cycles_Torig/c6/
# C30_n4000_Torig_grain_cross10_ftc6, final cost 4.781872e-03 on its base) on the lategrid base (user decision
# 2026-11-08: "refit the earlier clipped model (original temperature, with its templates) on the new base, with W10x5
# and the +1.5 m offset"; diagnostic: how much fit quality the gauge treatment costs vs the positive parameterization).
# Model = c6's command with the new base: H2O activity max(0, Y a + A_scale (c_cliff T_cliff + sum_k ft_ck)), 25
# degree-<=4 coefficients, cliff_ge_60deg free >= 0, ft_c1..ft_c6 held at 1, Tmax 182.09922 / T75 136.11149 / Tmin 133 K
# held, kappa free, grain + crossing 10, lambda 0, --tbounds 30 400, max_nfev 4000, the srcCG night/shadow ratios (no
# COPS_POS_*), no flux floor, --precision mixed (no record store on this base: kernels converted at bind).
# Driver model_ext/clip_refit.py (+ clip_problem.py) around the unchanged grain_fit_fixed.py / grain_fit.py / gpu_fit.py.
#   VARIANT=D1: model at the observation time stamps (c6's path; grid rows only in the all-rows model);
#   VARIANT=D2: D1 + W10x5 for the observations from 10:33:10 UTC exactly as positive_problem.py (its own operator).
#   VARIANT=D0 (control, user 2026-11-08 "yes"): c6 continued on its OWN base runs/shap8_012m_flat_earlier/base with its
#     record store, the D1 path (clip_refit.py --variant D1, no point table; = the path of the bitwise test 122854),
#     products runs/shap8_012m_flat_earlier/clip_refit/D0: the like-for-like reference for D1/D2 (same start, 4000
#     fresh evaluations).
# Start: c6's final a, c_cliff, kappa, grain by name (clip_refit.py --start-restart; prior unchanged); the start check
# prints the cost at the start next to c6's final cost and the change of the predictions.
# Base runs/shap8_012m_lategrid_dz1.5/base (1150 observations + 401 model-only 1-s grid rows; gauge +1.5 m), point
# table <run>/late_grid_points.txt.
#
# Settings (sbatch --export=ALL,VAR=value,...; values must not contain commas):
#   VARIANT       D1 or D2 (required)
#   NFEV          total evaluations of the fit over all segments (default 4000)
#   CHECK         1: in segment 1, grain_fit --check-jacobian (as c6) and clip_refit --check-full-jacobian (default 1)
#   START_RESTART the start state (segment 1 only; default the c6 fit_restart.npz)
#   TESTNAME      test runs: products in <run>/clip_refit_tests/<TESTNAME>/seg<NN>
#   MODEL_EXT     the model_ext directory the job copies and runs (default $COPS_OUTPUT_ROOT/model_ext); production runs
#                 use a frozen snapshot, MODEL_EXT=$COPS_OUTPUT_ROOT/model_ext_frozen/clip_<UTC> (never edited)
#   CHAIN         1: a segment directory without restart file stops the job (exit 3) instead of being moved aside
# Products: <run>/clip_refit/<VARIANT>/seg<NN> (+ link final -> seg<NN> when complete); provenance-<job>/ next to them
# (commands.sh with every argument resolved, logs, environment, git, nvidia-smi -q, uv pip list, the model_ext code it
# ran). The fit runs from that provenance copy of model_ext.
#
# Wall time / resume (as jobs/positive_fit.sbatch): the fit stops at its best evaluated state before the job's end,
# writes its products and an interrupted fit_restart.npz, exit 75. Resubmitting the SAME command continues with
# segment k+1 (gpu_fit --restart-from seg<k>/fit_restart.npz; clip_refit.py checks the marker and the start check must
# reproduce seg<k>'s cost to 1e-8), --max-nfev = NFEV - the evaluations of the earlier segments. A segment directory
# without restart file is moved aside as seg<k>.incomplete_<UTC> (never deleted) and run again. A segmented fit is not
# identical to one continuous run. A completed fit is not rerun (exit 0).
#
# Production (one whole node each, in parallel):
#   sbatch --nodelist=gpu001 --export=ALL,VARIANT=D1,MODEL_EXT=<snapshot> clip_refit.sbatch
#   sbatch --nodelist=gpu002 --export=ALL,VARIANT=D2,MODEL_EXT=<snapshot> clip_refit.sbatch
# Disk: every segment writes a full product set, ~3 GB on the share (2 GB activity facet tables, 0.44 GB restart).
#SBATCH --job-name=cops12m_cliprefit
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gpus-per-node=8
#SBATCH --cpus-per-task=64
#SBATCH --mem=0
#SBATCH --exclusive
#SBATCH --time=06:00:00
#SBATCH --partition=gpu
#SBATCH --account=general
#SBATCH --exclude=gpu003
#SBATCH --output=/data001/highresolu_vt/cops_output/jobs/logs/cops12m_cliprefit-%j.out
set -uo pipefail
export AMPS=/data001/highresolu_vt/AMPS APP=/data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda
export T=/data001/highresolu_vt/cops_cluster_transfer_20261026 COPS_OUTPUT_ROOT=/data001/highresolu_vt/cops_output
export COPS_DATA_ROOT=$T/rosetta_data_end_of_mission_cops
export CUDA_PATH=/usr/local/cuda-12.9 CUDA_HOME=/usr/local/cuda-12.9 CUPY_CACHE_DIR=$COPS_OUTPUT_ROOT/.cupy
export PATH=$CUDA_HOME/bin:/usr/mpi/gcc/openmpi-4.1.7rc1/bin:$PATH
export LD_LIBRARY_PATH=$CUDA_HOME/lib64${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1
# the model has the srcCG night/shadow ratios and no floor: the positive-activity GPU patch must be off
POSUNSET=$(env | grep -o '^COPS_POS_[A-Z_]*' | sort | paste -sd' ' -)
[ -n "$POSUNSET" ] && unset $POSUNSET
VARIANT=${VARIANT:-}; case $VARIANT in D0|D1|D2) ;; *) echo "set VARIANT=D0, D1 or D2: sbatch --export=ALL,VARIANT=D1 $0" >&2; exit 2;; esac
NFEV=${NFEV:-4000}; CHECK=${CHECK:-1}; TESTNAME=${TESTNAME:-}; CHAIN=${CHAIN:-0}
export MODEL_EXT=$(readlink -f "${MODEL_EXT:-$COPS_OUTPUT_ROOT/model_ext}")
R=$COPS_OUTPUT_ROOT/runs/shap8_012m_lategrid_dz1.5; E=$COPS_OUTPUT_ROOT/runs/shap8_012m_flat_earlier
X=$MODEL_EXT; PY="$APP/.venv/bin/python -W ignore"
BASE=$R/base; PTS=$R/late_grid_points.txt
# D0: c6's own base and record store, the D1 driver path without a point table (as the bitwise test 122854)
DRV=$VARIANT; PTA=(--point-table "$PTS"); STORE=()
if [ "$VARIANT" = D0 ]; then BASE=$E/base; DRV=D1; PTA=(); STORE=(--record-store "$E/store"); fi
C6=$E/grain/cycles_Torig/c6/C30_n4000_Torig_grain_cross10_ftc6
START_RESTART=$(readlink -f "${START_RESTART:-$C6/fit_restart.npz}")
# c6's templates and held values (fit_manifest.json model_extension.command of job 122221)
TA=(--activity-template "$E/bundle/restart_templates.npz" cliff_ge_60deg)
FIX=(--fix Tmax=182.09922 --fix T75=136.11149 --fix Tmin=133.0)
for j in 1 2 3 4 5 6; do
    if [ $j -le 3 ]; then TD=$E/grain/cycles/c$j; else TD=$E/grain/cycles_Torig/c$j; fi
    TA+=(--activity-template "$TD/template_ft_c$j/ft_c$j.templates.npz" ft_c$j); FIX+=(--fix c_ft_c$j=1)
done
if [ -n "$TESTNAME" ]; then NAME=clip_refit_test_$TESTNAME; D=$R/clip_refit_tests/$TESTNAME
elif [ "$VARIANT" = D0 ]; then NAME=clip_refit_D0; D=$E/clip_refit/D0
else NAME=clip_refit_$VARIANT; D=$R/clip_refit/$VARIANT; fi
P=$D/provenance-$SLURM_JOB_ID
FILES=(clip_refit.py clip_problem.py clip_compare.py grain_fit_fixed.py fixed_parameters.py grain_fit.py grain_crossing.py
       positive_problem.py late_response_eval.py test_clip_refit.py test_positive_fit.py)
for f in "$BASE/COPS_model.dat" "$BASE/elements.bin" "$BASE/kernels" "$PTS" "$START_RESTART" \
         $(printf '%s\n' "${TA[@]}" | grep '\.npz$') $(printf "$X/%s\n" "${FILES[@]}"); do
    [ -e "$f" ] || { echo "missing: $f" >&2; exit 2; }; done
mkdir -p "$D"
# ---- segment state: the last seg<NN>; its restart's status, the evaluations used so far, its first segment's start
state() {  # $1 = directory, $2 = segment (optional: the last one)
$PY - "$@" <<'PY'
import glob, json, os, sys
import numpy as np
d = sys.argv[1]
segs = sorted(int(os.path.basename(p)[3:]) for p in glob.glob(os.path.join(d, "seg[0-9][0-9]")))
k = int(sys.argv[2]) if len(sys.argv) > 2 else (segs[-1] if segs else 0)
if k == 0:
    print(0, "none", 0, "-"); sys.exit()
rs = os.path.join(d, "seg%02d" % k, "fit_restart.npz")
if not os.path.isfile(rs):
    print(k, "norestart", 0, "-"); sys.exit()
with np.load(rs, allow_pickle=False) as z:
    meta = json.loads(str(z["meta_json"]))
used = sum(int(s.get("nfev", 0)) for s in meta.get("segments", []))
chain = (meta.get("model_extension") or {}).get("fit_chain") or {}
print(k, meta.get("status", "unknown"), used, (chain.get("start_restart") or {}).get("path") or "-")
PY
}
read -r LAST LSTATE LUSED LSTART < <(state "$D")
[ -n "${LAST:-}" ] || { echo "cannot read the segment state of $D" >&2; exit 3; }
MOVED=""
if [ "$LSTATE" = norestart ]; then                     # a killed or failed attempt
    if [ "$CHAIN" = 1 ]; then
        echo "$D/seg$(printf %02d "$LAST") has no restart file (a failed attempt): chain stopped (CHAIN=1); inspect it" >&2
        exit 3
    fi
    MOVED="$D/seg$(printf %02d "$LAST").incomplete_$(date -u +%Y%m%dT%H%M%SZ)"   # move aside, run that segment again
    mv "$D/seg$(printf %02d "$LAST")" "$MOVED"
    LAST=$((LAST - 1))
    if [ "$LAST" -ge 1 ]; then read -r LAST LSTATE LUSED LSTART < <(state "$D" "$LAST")
    else LSTATE=none; LUSED=0; LSTART=-; fi
fi
if [ "$LSTATE" != none ] && [ "$LSTART" != "$START_RESTART" ]; then
    echo "START_RESTART=$START_RESTART but $D started from $LSTART: not the same fit" >&2; exit 3
fi
final_link() {   # the completed fit's last segment
    local seg=seg$(printf %02d "$1")
    if [ -L "$D/final" ]; then
        [ "$(readlink "$D/final")" = "$seg" ] || { echo "$D/final -> $(readlink "$D/final"), expected $seg" >&2; return 3; }
    else ln -s "$seg" "$D/final"; fi
    echo "$(date -u +%FT%TZ) $NAME complete: $D/final -> $seg"
}
case $LSTATE in
    none) SEG=1; LEFT=$NFEV;;
    interrupted) SEG=$((LAST + 1)); LEFT=$((NFEV - LUSED));;
    complete) echo "$(date -u +%FT%TZ) $NAME: segment $LAST is complete ($LUSED evaluations); nothing to do"
              final_link "$LAST"; exit $?;;
    *) echo "$D/seg$(printf %02d "$LAST")/fit_restart.npz has status '$LSTATE'; not continued" >&2; exit 3;;
esac
if [ "$LEFT" -le 0 ]; then echo "$(date -u +%FT%TZ) $NAME: $LUSED of $NFEV evaluations used; nothing to do"
    final_link "$LAST"; exit $?; fi
if [ "$SEG" = 1 ]; then START=(--start-restart "$START_RESTART")
else START=(--restart-from "$D/seg$(printf %02d "$LAST")/fit_restart.npz"); fi
O=$D/seg$(printf %02d "$SEG")
[ -e "$O" ] && { echo "exists: $O" >&2; exit 2; }
CHK=(); GCHK=(); [ "$SEG" = 1 ] && [ "$CHECK" = 1 ] && { CHK=(--check-full-jacobian); GCHK=(--check-jacobian); }
# ---- provenance; the fit runs from this copy of model_ext
mkdir -p "$P"; cd "$APP" || exit 1
cp -p "$0" "$P/clip_refit.sbatch"
for f in "${FILES[@]}"; do cp -p "$X/$f" "$P/"; done
sha256sum "$X"/*.py "$PTS" "$START_RESTART" > "$P/model_ext.sha256" 2>&1
(cd "$P" && sha256sum "${FILES[@]}") > "$P/run_copy.sha256" 2>&1
{ git -C "$AMPS" rev-parse HEAD; git -C "$AMPS" status --short; } > "$P/git.txt" 2>&1
module list > "$P/module_list.txt" 2>&1 || true
nvidia-smi -q > "$P/nvidia_smi_q.txt" 2>&1
{ nvcc --version; gcc --version | head -1; grep -m1 CSPICE /home/vtenishev/SPICE/cspice/doc/version.txt; } > "$P/toolchain.txt" 2>&1
~/.local/bin/uv pip list --python "$APP/.venv/bin/python" > "$P/pip_list.txt" 2>&1
env | grep -E '^(CUDA_[A-Z_]*|CUPY_[A-Z_]*|COPS_[A-Z_]*|VARIANT|NFEV|CHECK|TESTNAME|START_RESTART|CHAIN|MODEL_EXT|OMP_[A-Z_]*|OPENBLAS_[A-Z_]*|MKL_[A-Z_]*|PYTHONUNBUFFERED|LD_LIBRARY_PATH|PATH|SLURM_JOB_GPUS|CUDA_VISIBLE_DEVICES)=' | sort > "$P/environment.txt"
# commands.sh: every command of this job with all arguments resolved, runnable as is (HOW-TO-RERUN.md)
{ echo '#!/bin/bash'; echo "# job $SLURM_JOB_ID on $(hostname), $(date -u +%FT%TZ); AMPS $(git -C "$AMPS" rev-parse HEAD)"
  grep '^export ' "$0"; echo "cd $APP"; } > "$P/commands.sh"
rec() { printf '%q ' "$@" >> "$P/commands.sh"; echo >> "$P/commands.sh"; }
free -g > "$P/free_at_start.txt"; df -h "$R" > "$P/df_at_start.txt" 2>&1
{ echo "fit $NAME (variant $VARIANT), directory $D; start restart (segment 1) $START_RESTART"
  echo "model_ext copied and run: $X (MODEL_EXT)"
  [ -n "$POSUNSET" ] && echo "unset (inherited from the submitting shell): $POSUNSET"
  echo "segment $SEG; previous segment $LAST ($LSTATE, $LUSED evaluations used); evaluations for this segment $LEFT of $NFEV"
  [ -n "$MOVED" ] && echo "moved aside: $MOVED"; } | tee "$P/segment.txt"
# devices as CUDA numbers them inside the job
if [ -n "${CUDA_VISIBLE_DEVICES:-}" ]; then NG=$(echo "$CUDA_VISIBLE_DEVICES" | tr , "\n" | wc -l); else NG=$(nvidia-smi -L | wc -l); fi
DEV=$(seq -s, 0 $((NG - 1)))
echo "$(date -u +%FT%TZ) job $SLURM_JOB_ID on $(hostname): $NG GPUs visible ($DEV); fit $NAME segment $SEG ($LEFT evaluations)"
CMD=($PY "$P/clip_refit.py" --variant "$DRV" "${PTA[@]}" "${START[@]}" "${CHK[@]}" "${FIX[@]}" \
    --crossing-weight 10 "${GCHK[@]}" --devices "$DEV" --precision mixed "${STORE[@]}" --max-nfev "$LEFT" --segment "$SEG" \
    --wall-time slurm --template-start 0 -- "$BASE" \
    --profile all_available --degree 4 --enhancement-degree 4 --lambda 0 --jobs 1 --tbounds 30 400 \
    --ftol 1e-10 --xtol 1e-10 --gtol 1e-10 --check-tolerance 1e-6 --out "$O" "${TA[@]}")
rec bash -c 'echo 1000 > /proc/self/oom_score_adj; exec "$@"' _ "${CMD[@]}"
bash -c 'echo 1000 > /proc/self/oom_score_adj; exec "$@"' _ "${CMD[@]}" > "$P/fit.log" 2>&1; rc=$?
echo "$(date -u +%FT%TZ) fit segment $SEG: exit $rc"
grep -E '^(clip_refit|positive_fit|consistency|initial cost|lambda |gpu_model: bound|gpu_fit: |grain_fit: |fit data|selection profile)' "$P/fit.log" | cut -c1-500
tail -6 "$P/fit.log" | cut -c1-250
free -g > "$P/free_at_end.txt"; df -h "$R" > "$P/df_at_end.txt" 2>&1; du -sh "$O" >> "$P/df_at_end.txt" 2>&1
if [ $rc -eq 0 ]; then
    read -r _ st _ _ < <(state "$D" "$SEG")
    [ "$st" = complete ] && { final_link "$SEG" || rc=$?; }
fi
[ $rc -eq 75 ] && echo "stopped before the wall-time limit: resubmit the same command to continue with segment $((SEG + 1))"
exit $rc
