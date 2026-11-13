#!/bin/bash
# ==== version-controlled copy: header added 2026-11-13, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/jobs/nonneg_fit.sbatch (on Voyager: /data001/highresolu_vt/cops_output/jobs/nonneg_fit.sbatch)
# SHA256 of the run copy = SHA256 of this file without these header lines: 9584c5cae6f97b906c35616730678684a45ab0fa88b830f63aacda21a7ae7119
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: job. Softplus / squared fit on the lategrid base (user decisions 2026-11-11): MODEL required, START d2 | lsq |
#   chain, START_FILE, FT_FROM, HELD, GRAIN_FROM, NFEV 4000, CHECK, L_THREADS, MODEL_EXT, TESTMODE; everything else as
#   D2 (job 122866); segments, exit 75 = resubmit.
# RUN: sbatch --nodelist=gpu001 --export=ALL,MODEL=softplus,L_THREADS=32,MODEL_EXT=<snapshot> nonneg_fit.sbatch   (tests
#   124130-124133 inside nonneg_test.sbatch; split tests 124139, 124141, 124144, 124146; the production stages run its
#   logic inside chain_onejob.py)
# 2026-11-13: SLOPE_MULT_MIN (softplus; user 2026-11-12: multipliers >= 1) -> --slope-mult-min; GPU test 124659.
# ==== end of header ====
# Non-negative H2O activity fit on the lategrid base (user decisions 2026-11-11, handoff.md section 10; model_ext/
# nonneg_mult.py, nonneg_problem.py, driver nonneg_fit.py; adapted from jobs/slope_mult_fit.sbatch):
#   MODEL=softplus: A = m_b s ln(1 + exp(L/s)) exp(h), L = Y a + A_s sum ft_c1..ft_c6 (held 1), s = 1e-5 A_orig,
#                   hard clip max(0, L) where A_orig = 0
#   MODEL=squared:  A = m_b (A_s (Y b)^2 + A_s sum ft_c1..ft_c6) exp(h)   (a = A_s b)
#   m_b: 7 slope bins (bundle slope_bins.npz), m_[0,40) = 1 held, six free >= 0; NO floor; NO cliff_ge_60deg.
# Everything else = D2 (jobs/clip_refit.sbatch VARIANT=D2, job 122866): lategrid base, W10x5 from 10:33:10 UTC, grain at the
# time stamp, Tmax 182.09922 / T75 136.11149 / Tmin 133 K held, kappa free, grain + crossing 10, lambda 0, --tbounds 30
# 400, ftol/xtol/gtol 1e-10, --precision mixed, --check-tolerance 1e-6.
# Settings (sbatch --export=ALL,VAR=value,...):
#   MODEL         softplus | squared (required)
#   START         softplus: d2 (default) | chain; squared: lsq (default) | chain
#   START_RESTART the --start-restart file (default D2's final fit_restart.npz; chain: the previous cycle's fit restart)
#   START_FILE    squared, START=lsq: nonneg_start.py's output (default runs/.../nonneg/bundle/squared_start.npz)
#   FT_FROM       D2's fit_manifest.json (its ft_c* templates; default runs/.../clip_refit/D2/final/fit_manifest.json)
#   HELD, GRAIN_FROM  START=chain: the previous fine tuning's finetune_state.npz (held_next; grain A, B)
#   NFEV (default 4000), CHECK (1: segment 1 runs the full Jacobian check, default 1), CHECKPARAMS (':'-separated)
#   L_THREADS     N: every 12.4M-facet host step of the fitter with N chunked threads (bitwise; nonneg_fit.py --l-threads)
#   SLOPE_MULT_MIN  softplus only (user 2026-11-12): lower bound of the six free multipliers (nonneg_fit.py
#                 --slope-mult-min; the decision: 1); unset: the fitter's lower bound 0 (unchanged). A resume and a chain
#                 start must use the same value (the restart marker records it)
#   FITDIR        product directory (default runs/.../nonneg/<MODEL>/fit); TESTNAME: runs/.../nonneg_tests/<TESTNAME>
#   TESTMODE      independent (needs TESTNAME): stop after the start check, then an UNPATCHED process evaluates the
#                 independently computed activity (test_nonneg.py --gpu-independent); clip: the same with softplus
#                 s = 0 (--s-ratio 0 --test-clip)
#   MODEL_EXT     the model_ext copy to run (production: a frozen snapshot); CHAIN=1: a segment without restart stops
# Wall time / resume / move-aside / final link exactly as jobs/slope_mult_fit.sbatch (exit 75 = resubmit).
#SBATCH --job-name=cops12m_nonneg
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
#SBATCH --output=/data001/highresolu_vt/cops_output/jobs/logs/cops12m_nonneg-%j.out
set -uo pipefail
export AMPS=/data001/highresolu_vt/AMPS APP=/data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda
export T=/data001/highresolu_vt/cops_cluster_transfer_20261026 COPS_OUTPUT_ROOT=/data001/highresolu_vt/cops_output
export COPS_DATA_ROOT=$T/rosetta_data_end_of_mission_cops
export CUDA_PATH=/usr/local/cuda-12.9 CUDA_HOME=/usr/local/cuda-12.9 CUPY_CACHE_DIR=$COPS_OUTPUT_ROOT/.cupy
export PATH=$CUDA_HOME/bin:/usr/mpi/gcc/openmpi-4.1.7rc1/bin:$PATH
export LD_LIBRARY_PATH=$CUDA_HOME/lib64${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1
# this model patches the evaluator itself: no positive / slope patch, no inherited frozen config or shm names
UNSET=$(env | grep -o -E '^(COPS_POS_[A-Z_]*|COPS_SLOPEMUL_[A-Z_]*|COPS_NONNEG_[A-Z_]*)' | sort | paste -sd' ' -)
[ -n "$UNSET" ] && unset $UNSET
NFEV=${NFEV:-4000}; CHECK=${CHECK:-1}; CHECKPARAMS=${CHECKPARAMS:-}; TESTNAME=${TESTNAME:-}; TESTMODE=${TESTMODE:-}
CHAIN=${CHAIN:-0}; MODEL=${MODEL:-}; L_THREADS=${L_THREADS:-}; SLOPE_MULT_MIN=${SLOPE_MULT_MIN:-}
HELD=${HELD:-}; GRAIN_FROM=${GRAIN_FROM:-}; FITDIR=${FITDIR:-}
case $MODEL in softplus) START=${START:-d2};; squared) START=${START:-lsq};; *) echo "MODEL must be softplus or squared" >&2; exit 2;; esac
[ -n "$HELD" ] && HELD=$(readlink -f "$HELD"); [ -n "$GRAIN_FROM" ] && GRAIN_FROM=$(readlink -f "$GRAIN_FROM")
case $MODEL:$START in softplus:d2|softplus:chain|squared:lsq|squared:chain) ;; *) echo "START=$START is not a start of $MODEL" >&2; exit 2;; esac
if [ "$START" = chain ] && { [ -z "$HELD" ] || [ -z "$GRAIN_FROM" ]; }; then echo "START=chain needs HELD and GRAIN_FROM" >&2; exit 2; fi
BOPT=(); [ -n "$L_THREADS" ] && BOPT=(--l-threads "$L_THREADS" --l-mode chunked)
export MODEL_EXT=$(readlink -f "${MODEL_EXT:-$COPS_OUTPUT_ROOT/model_ext}")
R=$COPS_OUTPUT_ROOT/runs/shap8_012m_lategrid_dz1.5
X=$MODEL_EXT; PY="$APP/.venv/bin/python -W ignore"
BASE=$R/base; PTS=$R/late_grid_points.txt
BUNDLE=$(readlink -f "${BUNDLE:-$R/slope_mult/bundle/slope_bins.npz}")
FT_FROM=$(readlink -f "${FT_FROM:-$R/clip_refit/D2/final/fit_manifest.json}")
START_FILE=$(readlink -f "${START_FILE:-$R/nonneg/bundle/squared_start.npz}")
START_RESTART=$(readlink -f "${START_RESTART:-$R/clip_refit/D2/final/fit_restart.npz}")
FIX=(--fix Tmax=182.09922 --fix T75=136.11149 --fix Tmin=133.0 --fix c_slope_00_40=1)
case $TESTMODE in
    "") TOPT=();;
    independent) [ -n "$TESTNAME" ] || { echo "TESTMODE needs TESTNAME" >&2; exit 2; }; TOPT=(--stop-after-checks);;
    clip) [ -n "$TESTNAME" ] && [ "$MODEL" = softplus ] || { echo "TESTMODE=clip needs TESTNAME and MODEL=softplus" >&2; exit 2; }
          TOPT=(--s-ratio 0 --test-clip --stop-after-checks);;
    *) echo "TESTMODE must be empty, independent or clip" >&2; exit 2;;
esac
if [ -n "$TESTNAME" ]; then NAME=nonneg_test_$TESTNAME; D=$R/nonneg_tests/$TESTNAME
elif [ -n "$FITDIR" ]; then D=$(readlink -m "$FITDIR"); NAME=nonneg_fit_$(basename "$(dirname "$D")")_$(basename "$D")
else NAME=nonneg_${MODEL}_fit; D=$R/nonneg/$MODEL/fit; fi
MOPT=(--model "$MODEL" --ft-templates-from "$FT_FROM"); [ -n "$HELD" ] && MOPT+=(--held "$HELD")
[ "$START" = lsq ] && MOPT+=(--start-file "$START_FILE")
[ -n "$SLOPE_MULT_MIN" ] && MOPT+=(--slope-mult-min "$SLOPE_MULT_MIN")
P=$D/provenance-$SLURM_JOB_ID
FILES=(nonneg_fit.py nonneg_mult.py nonneg_problem.py nonneg_start.py nonneg_finetune.py test_nonneg.py
       slope_mult.py slope_mult_problem.py test_slope_mult.py positive_finetune.py clip_problem.py clip_refit.py
       grain_fit_fixed.py fixed_parameters.py grain_fit.py grain_crossing.py positive_problem.py late_response_eval.py
       test_clip_refit.py test_positive_fit.py)
LSQ_FILE=(); [ "$START" = lsq ] && LSQ_FILE=("$START_FILE")
for f in "$BASE/COPS_model.dat" "$BASE/elements.bin" "$BASE/kernels" "$PTS" "$START_RESTART" "$BUNDLE" "$FT_FROM" "${LSQ_FILE[@]}" \
         ${HELD:+"$HELD"} ${GRAIN_FROM:+"$GRAIN_FROM"} \
         $(printf "$X/%s\n" "${FILES[@]}"); do
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
    print(0, "none", 0, "-", "-"); sys.exit()
rs = os.path.join(d, "seg%02d" % k, "fit_restart.npz")
if not os.path.isfile(rs):
    print(k, "norestart", 0, "-", "-"); sys.exit()
with np.load(rs, allow_pickle=False) as z:
    meta = json.loads(str(z["meta_json"]))
used = sum(int(s.get("nfev", 0)) for s in meta.get("segments", []))
chain = (meta.get("model_extension") or {}).get("fit_chain") or {}
print(k, meta.get("status", "unknown"), used, (chain.get("start_restart") or {}).get("path") or "-",
      (meta.get("model_extension") or {}).get("start_choice") or chain.get("start_choice") or "-")
PY
}
if [ -n "$TESTMODE" ]; then                        # a check, not a fit: one directory, never resumed
    [ -e "$D/seg01" ] && { echo "exists: $D/seg01 (choose a new TESTNAME)" >&2; exit 2; }
    LAST=0; LSTATE=none; LUSED=0; LSTART=-; LCHOICE=-; MOVED=""
else
read -r LAST LSTATE LUSED LSTART LCHOICE < <(state "$D")
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
    if [ "$LAST" -ge 1 ]; then read -r LAST LSTATE LUSED LSTART LCHOICE < <(state "$D" "$LAST")
    else LSTATE=none; LUSED=0; LSTART=-; LCHOICE=-; fi
fi
if [ "$LSTATE" != none ] && [ "$LSTART" != "$START_RESTART" ]; then
    echo "START_RESTART=$START_RESTART but $D started from $LSTART: not the same fit" >&2; exit 3
fi
if [ "$LSTATE" != none ] && [ "$LCHOICE" != "$START" ]; then
    echo "START=$START but $D started with --start $LCHOICE: not the same fit" >&2; exit 3
fi
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
if [ "$SEG" = 1 ]; then SOPT=(--start "$START" --start-restart "$START_RESTART")
    [ "$START" = chain ] && SOPT+=(--grain-from "$GRAIN_FROM")
else SOPT=(--start "$START" --restart-from "$D/seg$(printf %02d "$LAST")/fit_restart.npz"); fi
O=$D/seg$(printf %02d "$SEG")
[ -e "$O" ] && { echo "exists: $O" >&2; exit 2; }
CHK=(); GCHK=(); CHKP=()
if [ "$SEG" = 1 ] && [ "$CHECK" = 1 ]; then
    CHK=(--check-full-jacobian); GCHK=(--check-jacobian)
    [ -n "$CHECKPARAMS" ] && CHKP=(--check-parameters "${CHECKPARAMS//:/,}")
fi
# ---- provenance; the fit runs from this copy of model_ext
mkdir -p "$P"; cd "$APP" || exit 1
cp -p "$0" "$P/nonneg_fit.sbatch"
for f in "${FILES[@]}"; do cp -p "$X/$f" "$P/"; done
sha256sum "$X"/*.py "$PTS" "$START_RESTART" "$BUNDLE" "$FT_FROM" "${LSQ_FILE[@]}" > "$P/model_ext.sha256" 2>&1
(cd "$P" && sha256sum "${FILES[@]}") > "$P/run_copy.sha256" 2>&1
{ git -C "$AMPS" rev-parse HEAD; git -C "$AMPS" status --short; } > "$P/git.txt" 2>&1
module list > "$P/module_list.txt" 2>&1 || true
nvidia-smi -q > "$P/nvidia_smi_q.txt" 2>&1
nvidia-smi -L > "$P/nvidia_smi_L.txt" 2>&1
{ nvcc --version; gcc --version | head -1; grep -m1 CSPICE /home/vtenishev/SPICE/cspice/doc/version.txt; } > "$P/toolchain.txt" 2>&1
~/.local/bin/uv pip list --python "$APP/.venv/bin/python" > "$P/pip_list.txt" 2>&1
env | grep -E '^(CUDA_[A-Z_]*|CUPY_[A-Z_]*|COPS_[A-Z_]*|NFEV|CHECK|CHECKPARAMS|TESTNAME|TESTMODE|MODEL|START|START_RESTART|START_FILE|FT_FROM|L_THREADS|SLOPE_MULT_MIN|HELD|GRAIN_FROM|FITDIR|BUNDLE|CHAIN|MODEL_EXT|OMP_[A-Z_]*|OPENBLAS_[A-Z_]*|MKL_[A-Z_]*|PYTHONUNBUFFERED|LD_LIBRARY_PATH|PATH|SLURM_JOB_GPUS|SLURM_JOB_NODELIST|SLURM_JOB_NUM_NODES|CUDA_VISIBLE_DEVICES)=' | sort > "$P/environment.txt"
# commands.sh: every command of this job with all arguments resolved, runnable as is (HOW-TO-RERUN.md)
{ echo '#!/bin/bash'; echo "# job $SLURM_JOB_ID on $(hostname), $(date -u +%FT%TZ); AMPS $(git -C "$AMPS" rev-parse HEAD)"
  grep '^export ' "$0"; echo "cd $APP"; } > "$P/commands.sh"
rec() { printf '%q ' "$@" >> "$P/commands.sh"; echo >> "$P/commands.sh"; }
free -g > "$P/free_at_start.txt"; df -h "$R" > "$P/df_at_start.txt" 2>&1
{ echo "fit $NAME, model $MODEL, directory $D; start $START; start restart (segment 1) $START_RESTART; slope bundle $BUNDLE"
  echo "ft templates from $FT_FROM; start file ${LSQ_FILE[*]:-none}; held ${HELD:-none}; grain from (segment 1) ${GRAIN_FROM:-the start restart}"
  echo "host steps of the fitter process: ${L_THREADS:-1} thread(s) (chunked, bitwise); OPENBLAS_NUM_THREADS=$OPENBLAS_NUM_THREADS"
  echo "multiplier lower bound: ${SLOPE_MULT_MIN:-0 (the fitter default, unchanged)}"
  echo "model_ext copied and run: $X (MODEL_EXT); test mode: ${TESTMODE:-none}"
  echo "node $(hostname) (SLURM_JOB_NODELIST ${SLURM_JOB_NODELIST:-?}, nodes ${SLURM_JOB_NUM_NODES:-1}); transport local"
  [ -n "$UNSET" ] && echo "unset (inherited from the submitting shell): $UNSET"
  echo "segment $SEG; previous segment $LAST ($LSTATE, $LUSED evaluations used); evaluations for this segment $LEFT of $NFEV"
  [ -n "$MOVED" ] && echo "moved aside: $MOVED"; } | tee "$P/segment.txt"
# devices as CUDA numbers them inside the job
if [ -n "${CUDA_VISIBLE_DEVICES:-}" ]; then NG=$(echo "$CUDA_VISIBLE_DEVICES" | tr , "\n" | wc -l); else NG=$(nvidia-smi -L | wc -l); fi
DEV=$(seq -s, 0 $((NG - 1)))
echo "$(date -u +%FT%TZ) job $SLURM_JOB_ID on $(hostname): $NG GPUs visible ($DEV); $NAME segment $SEG ($LEFT evaluations)"
( while true; do echo "$(date -u +%T) $(free -m | awk '/Mem/{print $3, $7}')"; sleep 20; done ) > "$P/memory_used_avail_MB.txt" 2>&1 &
MEMLOG=$!
( while true; do echo "== $(date -u +%T) loadavg $(cat /proc/loadavg)"
    ps -u "$(id -un)" -o pid,ppid,nlwp,pcpu,rss,etimes,comm --sort=-pcpu | awk 'NR == 1 || $7 ~ /python/' | head -12
    sleep 10; done ) > "$P/process_load.txt" 2>&1 &
LOADLOG=$!
CMD=($PY "$P/nonneg_fit.py" --point-table "$PTS" --slope-bundle "$BUNDLE" "${MOPT[@]}" "${SOPT[@]}" "${BOPT[@]}" "${CHK[@]}" "${CHKP[@]}" \
    "${TOPT[@]}" "${FIX[@]}" --crossing-weight 10 "${GCHK[@]}" --devices "$DEV" --precision mixed --max-nfev "$LEFT" \
    --segment "$SEG" --wall-time slurm --template-start 0 -- "$BASE" \
    --profile all_available --degree 4 --enhancement-degree 4 --lambda 0 --jobs 1 --tbounds 30 400 \
    --ftol 1e-10 --xtol 1e-10 --gtol 1e-10 --check-tolerance 1e-6 --out "$O")
rec bash -c 'echo 1000 > /proc/self/oom_score_adj; exec "$@"' _ "${CMD[@]}"
bash -c 'echo 1000 > /proc/self/oom_score_adj; exec "$@"' _ "${CMD[@]}" > "$P/fit.log" 2>&1; rc=$?
echo "$(date -u +%FT%TZ) $NAME segment $SEG: exit $rc"
grep -E '^(nonneg|clip_refit|consistency|initial cost|lambda |gpu_model: bound|gpu_fit: |grain_fit: |fit data|selection profile)' "$P/fit.log" | cut -c1-600
tail -6 "$P/fit.log" | cut -c1-250
if [ -n "$TESTMODE" ] && [ $rc -eq 0 ]; then
    ICMD=($PY "$P/test_nonneg.py" --gpu-independent "$O" --devices "$DEV")
    rec bash -c 'echo 1000 > /proc/self/oom_score_adj; exec "$@"' _ "${ICMD[@]}"
    bash -c 'echo 1000 > /proc/self/oom_score_adj; exec "$@"' _ "${ICMD[@]}" > "$P/independent.log" 2>&1; rc=$?
    echo "$(date -u +%FT%TZ) independent evaluation: exit $rc"
    grep -E 'PASS|FAIL|SUMMARY|gpu_model: bound|consistency' "$P/independent.log" | cut -c1-700
fi
kill $MEMLOG $LOADLOG 2>/dev/null
free -g > "$P/free_at_end.txt"; df -h "$R" > "$P/df_at_end.txt" 2>&1; du -sh "$O" >> "$P/df_at_end.txt" 2>&1
echo "peak used host memory: $(sort -k2 -n "$P/memory_used_avail_MB.txt" | tail -1) (MB used, available)"
if [ $rc -eq 0 ] && [ -z "$TESTMODE" ]; then
    read -r _ st _ _ _ < <(state "$D" "$SEG")
    [ "$st" = complete ] && { final_link "$SEG" || rc=$?; }
fi
[ $rc -eq 75 ] && echo "stopped before the wall-time limit: resubmit the same command to continue with segment $((SEG + 1))"
exit $rc
