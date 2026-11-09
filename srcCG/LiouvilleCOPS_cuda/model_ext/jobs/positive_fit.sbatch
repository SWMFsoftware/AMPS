#!/bin/bash
# ==== version-controlled copy: header added 2026-11-09, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/jobs/positive_fit.sbatch (on Voyager: /data001/highresolu_vt/cops_output/jobs/positive_fit.sbatch)
# SHA256 of the run copy = SHA256 of this file without these header lines: 21fa8494c2653966712372cce2abf061e1fc42a66fd69481b23e19222cc8232f
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: job. Positive-activity fit (user decisions 2026-11-07/08); LAMBDA required; RATIO (0.02 default | 0.002), CYCLE, START_RESTART,
#   GRAIN_A_UPA / GRAIN_B / GRAIN_FROM, HELD, MODEL_EXT (frozen code; the job runs from its provenance copy), TESTNAME,
#   CHAIN; NFEV 4000; resume segments seg<NN>; products positive_r<R>/lam<L>/c<k>/fit/seg<NN> (+ final link).
# RUN: sbatch --nodelist=gpu001 --export=ALL,MODEL_EXT=<snapshot>,LAMBDA=1e-3,RATIO=0.002,CYCLE=1,START_RESTART=<npz> positive_fit.sbatch
#   (cycle-1 fits 122449-122454; later stages via positive_cycles.sh)
# 2026-11-09 update: ratio 1e-5 (user 2026-11-08), CO2 flux floor (COPS_POS_FLOOR_CO2), start override --start-set (kappa = 1).
# ==== end of header ====
# Production 12M fits with the positive H2O activity (user decisions 2026-11-07), one Tikhonov lambda per job:
# A = exp(Y a + c_cliff T_cliff) (cliff_ge_60deg inside the exponent), start and prior = the emission-weighted
# (w = area A_orig) LS projection of ln A_orig, original temperature law held, kappa free, log-space Tikhonov
# D^2 = sum area (s - s0)^2 / sum area, the COPS W10x5 time response in the late window (>= 10:33:10 UTC; 1-s grid rows
# of the lategrid base), grain + crossing 10; GPU patch model_ext/positive_activity.py (night/shadow ratio = H2O flux
# floor ratio RATIO; env COPS_POS_*). model_ext/positive_fit.py around the unchanged grain_fit.py / gpu_fit.py /
# fit_cops.py. Base: runs/shap8_012m_lategrid_dz1.5/base (no record store: --precision mixed converts the kernels at
# bind), point table <run>/late_grid_points.txt, the flat_earlier bundle.
#
# Settings (sbatch --export=ALL,VAR=value,...; values must not contain commas):
#   LAMBDA        required, one value per job; used verbatim in the directory name (1e-2, not 0.01, for the cycles)
#   NFEV          total evaluations of the fit over all segments (default 4000)
#   CHECK         1: --check-jacobian in segment 1 (default 1)
#   RATIO         night/shadow ratio = H2O floor ratio: 0.02 (default; the behaviour before 2026-11-07 evening) or
#                 0.002 (user 2026-11-07, second decision) or 1e-5 (user 2026-11-08, third decision). Other values are
#                 refused (positive_problem.RATIO_DECISIONS). RATIO=1e-5 also sets the CO2 flux floor with the same
#                 ratio (user 2026-11-08 "use the same floor for CO2": COPS_POS_FLOOR_CO2=1, F_CO2 = max(kappa Q_CO2
#                 ACO2, 1e-5 q_max,CO2 ACO2); positive_problem.CO2_FLOOR_DECISIONS); for 0.02 / 0.002 it is unset
#   CYCLE         k of the fit/fine-tuning cycles: products in <run>/positive_r<RATIO>/lam<LAMBDA>/c<k>/fit/seg<NN>.
#                 Required whenever RATIO != 0.02 or START_RESTART / GRAIN_* / HELD is given (no collision with the
#                 first scan)
#   START_RESTART <fit_restart.npz> of another positive fit (segment 1 only): the optimizer starts from its a, c, kappa,
#                 grain (prior s0 unchanged); positive_fit.py's start check compares the start cost with the source's
#                 final cost (same ratio/lambda: must reproduce it to 1e-8, else the job fails; otherwise reported)
#   START_SET     (2026-11-08; user "kappa = 1 where stuck") NAME=VALUE[:NAME=VALUE...] (colon-separated; sbatch
#                 --export splits commas), e.g. START_SET=kappa=1: positive_fit.py --start-set in segment 1 (needs
#                 START_RESTART): the optimizer start of the named parameters is VALUE, everything else from the start
#                 restart; recorded in the manifest, restart marker (fit_chain), ACTIVITY_PARAMETERIZATION.txt, the start
#                 check (NOT the source's state), commands.sh, environment.txt, segment.txt. A resume must name the same
#                 START_SET or none (checked against the restart's fit_chain). Absent: nothing changes
#   GRAIN_A_UPA, GRAIN_B   explicit grain start (1e-6 Pa; 1e-6 Pa per 100 s), e.g. from a fine tuning; they win over
#                 the start restart's grain. Default: the start restart's grain, else A = 6.96024299190473, B = 0
#   GRAIN_FROM    <finetune_state.npz>: the grain start from its keys grain_A_uPa, grain_B_uPa_per_100s (read at run
#                 time; excludes GRAIN_A_UPA/GRAIN_B)
#   HELD          <npz> with key held_next (the fine tuning's finetune_state.npz) or held: the per-facet term held fixed
#                 in the exponent, s = Y a + c T + h (cycles 2, 3); D measures the total ln A from s0. Every segment of
#                 the fit must have the same HELD (a resume checks the sha256 recorded in the restart)
#   CHECK_HELD_SPAN 1: diagnostic --check-held-span (tests with a held term in the span of the fitted columns)
#   TESTNAME      test runs: products in <run>/positive_tests/<TESTNAME>/seg<NN> (instead of the CYCLE/first-scan name)
#   MODEL_EXT     the model_ext directory the job copies and runs (default $COPS_OUTPUT_ROOT/model_ext); production
#                 runs use a frozen snapshot, e.g. MODEL_EXT=$COPS_OUTPUT_ROOT/model_ext_frozen/c1_<UTC> (never edited),
#                 so later edits of model_ext cannot reach a job that is still pending. Recorded in provenance
#                 (environment.txt, segment.txt, model_ext.sha256, run_copy.sha256)
#   CHAIN         1 (set by jobs/positive_cycles.sh): a segment directory without restart file (a failed attempt) stops
#                 the job (exit 3) instead of being moved aside and rerun, so a dependency chain stops on real errors
# Without CYCLE/TESTNAME (and RATIO 0.02, no START_RESTART/GRAIN_*): the first scan's name, unchanged:
#   <run>/positive/C30_n<NFEV>_exp_cliff_w10x5_grain_cross10_lam<LAMBDA>/seg<NN>
# Every job: <fit dir>/provenance-<job>/ (commands.sh, logs, environment, the model_ext code it ran). The fit runs from
# that provenance copy of model_ext (positive_fit.py, positive_problem.py, positive_activity.py, kernels/, grain_*.py),
# so later edits of model_ext do not reach a running job or its spawned GPU workers.
#
# Wall time / resume (the gpu_fit --wall-time slurm convention, as tools/cluster/cops_pipeline.py does it): the fit
# stops at its best evaluated state before the job's end, writes its products and an interrupted fit_restart.npz
# (meta status "interrupted", cumulative "segments") and the job exits 75. Resubmitting the SAME command continues: the
# next segment seg<k+1> starts with gpu_fit --restart-from seg<k>/fit_restart.npz (positive_fit.py accepts only its own
# restart marker with this job's ratio), --max-nfev = NFEV - the evaluations of all earlier segments, --segment k+1;
# the prior is unchanged. START_RESTART applies to segment 1 only (a resume must name the same START_RESTART or none).
# A segment directory without a restart file (a killed or failed attempt) is moved aside as seg<k>.incomplete_<UTC>
# (never deleted) and that segment runs again from seg<k-1> (with CHAIN=1: exit 3 instead). A segmented fit is not
# identical to one continuous run. A completed fit (restart status "complete", or no evaluations left) is not rerun
# (exit 0); it gets the link <fit dir>/final -> seg<NN> (the fine tuning and the next cycle read <fit dir>/final).
#
# Examples (each job takes one whole GPU node):
#   first scan (unchanged): sbatch --export=ALL,LAMBDA=0 positive_fit.sbatch
#   cycle 1 at ratio 0.002 from the first scan's lambda 1e-2 fit:
#     sbatch --export=ALL,LAMBDA=1e-2,RATIO=0.002,CYCLE=1,START_RESTART=<run>/positive/C30_n4000_exp_cliff_w10x5_grain_cross10_lam1e-2/seg01/fit_restart.npz positive_fit.sbatch
#   resume after exit 75: the same sbatch command again.
# Disk: every segment writes a full product set, ~6 GB apparent (~3.3 GB on the compressed share); check df first.
#SBATCH --job-name=cops12m_posfit
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gpus-per-node=8
#SBATCH --cpus-per-task=64
#SBATCH --mem=0
#SBATCH --exclusive
#SBATCH --time=06:00:00
#SBATCH --partition=gpu
#SBATCH --account=general
#SBATCH --output=/data001/highresolu_vt/cops_output/jobs/logs/cops12m_posfit-%j.out
set -uo pipefail
export AMPS=/data001/highresolu_vt/AMPS APP=/data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda
export T=/data001/highresolu_vt/cops_cluster_transfer_20261026 COPS_OUTPUT_ROOT=/data001/highresolu_vt/cops_output
export COPS_DATA_ROOT=$T/rosetta_data_end_of_mission_cops
export CUDA_PATH=/usr/local/cuda-12.9 CUDA_HOME=/usr/local/cuda-12.9 CUPY_CACHE_DIR=$COPS_OUTPUT_ROOT/.cupy
export PATH=$CUDA_HOME/bin:/usr/mpi/gcc/openmpi-4.1.7rc1/bin:$PATH
export LD_LIBRARY_PATH=$CUDA_HOME/lib64${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1
[ -n "${LAMBDA:-}" ] || { echo "set LAMBDA (one value per job): sbatch --export=ALL,LAMBDA=0 $0" >&2; exit 2; }
NFEV=${NFEV:-4000}; CHECK=${CHECK:-1}; RATIO=${RATIO:-0.02}; CYCLE=${CYCLE:-}; TESTNAME=${TESTNAME:-}
START_RESTART=${START_RESTART:-}; GRAIN_A_UPA=${GRAIN_A_UPA:-}; GRAIN_B=${GRAIN_B:-}; CHAIN=${CHAIN:-0}
GRAIN_FROM=${GRAIN_FROM:-}; HELD=${HELD:-}; CHECK_HELD_SPAN=${CHECK_HELD_SPAN:-0}; START_SET=${START_SET:-}
case $START_SET in *[!A-Za-z0-9_=.:+-]*) echo "START_SET=$START_SET: NAME=VALUE[:NAME=VALUE...] only" >&2; exit 2;; esac
[ -z "$START_SET" ] || [ -n "$START_RESTART" ] || { echo "START_SET needs START_RESTART (it overrides that start)" >&2; exit 2; }
case $RATIO in 0.02|0.002|1e-5) ;; *) echo "RATIO=$RATIO: the user decided 0.02, 0.002 and 1e-5 (write it exactly so)" >&2; exit 2;; esac
# the GPU patch's configuration (user 2026-11-07): one ratio for night/shadow and the H2O floor; positive_fit.py
# checks these values (--ratio below) and refuses undecided ones
export COPS_POS_ENABLE=1 COPS_POS_NIGHT_RATIO=$RATIO COPS_POS_ACTIVITY=exp COPS_POS_FLOOR_RATIO=$RATIO
# the CO2 flux floor is part of the 1e-5 decision only (user 2026-11-08; positive_problem.CO2_FLOOR_DECISIONS refuses
# any other combination)
if [ "$RATIO" = 1e-5 ]; then export COPS_POS_FLOOR_CO2=1; CO2FLOOR="on (COPS_POS_FLOOR_CO2=1, ratio $RATIO)"
else unset COPS_POS_FLOOR_CO2; CO2FLOOR="off (COPS_POS_FLOOR_CO2 unset)"; fi
export MODEL_EXT=$(readlink -f "${MODEL_EXT:-$COPS_OUTPUT_ROOT/model_ext}")
R=$COPS_OUTPUT_ROOT/runs/shap8_012m_lategrid_dz1.5; X=$MODEL_EXT; PY="$APP/.venv/bin/python -W ignore"
BASE=$R/base; PTS=$R/late_grid_points.txt
BUNDLE=$COPS_OUTPUT_ROOT/runs/shap8_012m_flat_earlier/bundle/restart_templates.npz
if [ -n "$TESTNAME" ]; then
    NAME=test_$TESTNAME; D=$R/positive_tests/$TESTNAME
elif [ -n "$CYCLE" ]; then
    case $CYCLE in [1-9]|[1-9][0-9]) ;; *) echo "CYCLE=$CYCLE must be a positive integer" >&2; exit 2;; esac
    NAME=positive_r${RATIO}_lam${LAMBDA}_c${CYCLE}; D=$R/positive_r$RATIO/lam$LAMBDA/c$CYCLE/fit
else
    if [ "$RATIO" != 0.02 ] || [ -n "$START_RESTART$GRAIN_A_UPA$GRAIN_B$GRAIN_FROM$HELD" ]; then
        echo "RATIO $RATIO / START_RESTART / GRAIN_* / HELD: set CYCLE=<k> (products in positive_r<RATIO>/lam<L>/c<k>/fit)" >&2; exit 2
    fi
    NAME=C30_n${NFEV}_exp_cliff_w10x5_grain_cross10_lam${LAMBDA}; D=$R/positive/$NAME
fi
P=$D/provenance-$SLURM_JOB_ID
for f in "$BASE/COPS_model.dat" "$BASE/elements.bin" "$BASE/kernels" "$PTS" "$BUNDLE" "$X/positive_fit.py" \
         "$X/positive_problem.py" "$X/positive_activity.py" "$X/kernels" ${START_RESTART:+"$START_RESTART"} \
         ${HELD:+"$HELD"} ${GRAIN_FROM:+"$GRAIN_FROM"}; do
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
      ",".join("%s=%r" % (e["name"], float(e["value"])) for e in chain.get("start_set") or []) or "-")
PY
}
read -r LAST LSTATE LUSED LSTART LSSET < <(state "$D")
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
    if [ "$LAST" -ge 1 ]; then read -r LAST LSTATE LUSED LSTART LSSET < <(state "$D" "$LAST")
    else LSTATE=none; LUSED=0; LSTART=-; LSSET=-; fi
fi
if [ "$LSTATE" != none ] && [ -n "$START_RESTART" ] && [ "$LSTART" != "$(readlink -f "$START_RESTART")" ]; then
    echo "START_RESTART=$START_RESTART but $D started from $LSTART: not the same fit" >&2; exit 3
fi
if [ -n "$START_SET" ]; then          # canonical NAME=VALUE list, as the restart's fit_chain records it
    SSET=$($PY -c 'import sys; print(",".join("%s=%r" % (n.strip(), float(v)) for n, v in (i.split("=") for i in sys.argv[1].split(":") if i)))' "$START_SET") || { echo "START_SET=$START_SET: cannot parse" >&2; exit 2; }
    if [ "$LSTATE" != none ] && [ "$LSSET" != "$SSET" ]; then
        echo "START_SET=$START_SET ($SSET) but $D started with --start-set ${LSSET}: not the same fit" >&2; exit 3
    fi
fi
final_link() {   # the completed fit's last segment, for the fine tuning and the next cycle
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
if [ "$SEG" = 1 ]; then
    START=(); [ -n "$START_RESTART" ] && START=(--start-restart "$(readlink -f "$START_RESTART")")
    [ -n "$START_SET" ] && START+=(--start-set "${START_SET//:/,}")
else
    START=(--restart-from "$D/seg$(printf %02d "$LAST")/fit_restart.npz")
fi
GRAIN=(); [ -n "$GRAIN_A_UPA" ] && GRAIN+=(--grain-start-uPa "$GRAIN_A_UPA"); [ -n "$GRAIN_B" ] && GRAIN+=(--grain-start-B "$GRAIN_B")
[ -n "$GRAIN_FROM" ] && GRAIN+=(--grain-start-from "$(readlink -f "$GRAIN_FROM")")
[ "$SEG" = 1 ] && [ -z "$START_RESTART" ] && [ -z "$GRAIN_A_UPA$GRAIN_FROM" ] && GRAIN+=(--grain-start-uPa 6.96024299190473)
HOLD=(); [ -n "$HELD" ] && HOLD=(--held-exponent "$(readlink -f "$HELD")")
[ "$CHECK_HELD_SPAN" = 1 ] && [ "$SEG" = 1 ] && HOLD+=(--check-held-span)
O=$D/seg$(printf %02d "$SEG")
[ -e "$O" ] && { echo "exists: $O" >&2; exit 2; }
CHK=(); [ "$SEG" = 1 ] && [ "$CHECK" = 1 ] && CHK=(--check-jacobian)
# ---- provenance (as grain_late_weights.sbatch); the fit runs from this copy of model_ext
mkdir -p "$P"; cd "$APP" || exit 1
cp -p "$0" "$P/positive_fit.sbatch"
cp -p "$X"/positive_fit.py "$X"/positive_problem.py "$X"/positive_activity*.* "$X"/grain_fit.py "$X"/grain_crossing.py \
    "$X"/late_response_eval.py "$X"/late_time_grid.py "$X"/test_positive_fit.py "$P/"
cp -rp "$X/kernels" "$P/kernels"          # positive_activity.py's patched CUDA sources (model_ext/kernels/)
sha256sum "$X"/*.py "$X"/kernels/* "$PTS" > "$P/model_ext.sha256" 2>&1
(cd "$P" && sha256sum positive_fit.py positive_problem.py positive_activity.py grain_fit.py grain_crossing.py kernels/*) \
    > "$P/run_copy.sha256" 2>&1
{ git -C "$AMPS" rev-parse HEAD; git -C "$AMPS" status --short; } > "$P/git.txt" 2>&1
module list > "$P/module_list.txt" 2>&1 || true
nvidia-smi -q > "$P/nvidia_smi_q.txt" 2>&1
{ nvcc --version; gcc --version | head -1; grep -m1 CSPICE /home/vtenishev/SPICE/cspice/doc/version.txt; } > "$P/toolchain.txt" 2>&1
~/.local/bin/uv pip list --python "$APP/.venv/bin/python" > "$P/pip_list.txt" 2>&1
env | grep -E '^(CUDA_[A-Z_]*|CUPY_[A-Z_]*|COPS_[A-Z0-9_]*|LAMBDA|NFEV|CHECK|CHECK_HELD_SPAN|RATIO|CYCLE|TESTNAME|START_RESTART|START_SET|GRAIN_[A-Z_]*|HELD|CHAIN|MODEL_EXT|OMP_[A-Z_]*|OPENBLAS_[A-Z_]*|MKL_[A-Z_]*|PYTHONUNBUFFERED|LD_LIBRARY_PATH|PATH|SLURM_JOB_GPUS|CUDA_VISIBLE_DEVICES)=' | sort > "$P/environment.txt"
# commands.sh: every command of this job with all arguments resolved, runnable as is (HOW-TO-RERUN.md)
{ echo '#!/bin/bash'; echo "# job $SLURM_JOB_ID on $(hostname), $(date -u +%FT%TZ); AMPS $(git -C "$AMPS" rev-parse HEAD)"
  grep '^export ' "$0" | grep -v COPS_POS_
  echo "export COPS_POS_ENABLE=1 COPS_POS_NIGHT_RATIO=$RATIO COPS_POS_ACTIVITY=exp COPS_POS_FLOOR_RATIO=$RATIO"
  if [ "$RATIO" = 1e-5 ]; then echo "export COPS_POS_FLOOR_CO2=1"; else echo "unset COPS_POS_FLOOR_CO2"; fi
  echo "cd $APP"; } > "$P/commands.sh"
rec() { printf '%q ' "$@" >> "$P/commands.sh"; echo >> "$P/commands.sh"; }
free -g > "$P/free_at_start.txt"; df -h "$R" > "$P/df_at_start.txt" 2>&1
{ echo "fit $NAME, directory $D; ratio $RATIO; cycle ${CYCLE:--}; start restart ${START_RESTART:--}; grain A ${GRAIN_A_UPA:--} B ${GRAIN_B:--} from ${GRAIN_FROM:--}; held ${HELD:--}"
  echo "CO2 flux floor: $CO2FLOOR"
  [ -n "$START_SET" ] && echo "start set: $START_SET (segment 1: --start-set; the fit's chain: ${LSSET:--})"
  echo "model_ext copied and run: $X (MODEL_EXT)"
  [ -n "$HELD" ] && echo "held file sha256 $(sha256sum "$HELD" | cut -c1-64)"
  echo "segment $SEG; previous segment $LAST ($LSTATE, $LUSED evaluations used); evaluations for this segment $LEFT of $NFEV"
  [ -n "$MOVED" ] && echo "moved aside: $MOVED"; } | tee "$P/segment.txt"
# devices as CUDA numbers them inside the job (on a shared node nvidia-smi lists all GPUs; job 122226 failed with them)
if [ -n "${CUDA_VISIBLE_DEVICES:-}" ]; then NG=$(echo "$CUDA_VISIBLE_DEVICES" | tr , "\n" | wc -l); else NG=$(nvidia-smi -L | wc -l); fi
DEV=$(seq -s, 0 $((NG - 1)))
echo "$(date -u +%FT%TZ) job $SLURM_JOB_ID on $(hostname): $NG GPUs visible ($DEV); fit $NAME segment $SEG (lambda $LAMBDA, ratio $RATIO, CO2 floor $CO2FLOOR, $LEFT evaluations)"
CMD=($PY "$P/positive_fit.py" --point-table "$PTS" --ratio "$RATIO" "${START[@]}" "${GRAIN[@]}" "${HOLD[@]}" --crossing-weight 10 \
    "${CHK[@]}" --devices "$DEV" --precision mixed --max-nfev "$LEFT" --segment "$SEG" --wall-time slurm \
    --template-start 0 -- "$BASE" \
    --profile all_available --degree 4 --enhancement-degree 4 --lambda "$LAMBDA" --jobs 1 --fix-temperature \
    --tbounds 30 400 --ftol 1e-10 --xtol 1e-10 --gtol 1e-10 --check-tolerance inf --out "$O" \
    --activity-template "$BUNDLE" cliff_ge_60deg)
rec "${CMD[@]}"
bash -c 'echo 1000 > /proc/self/oom_score_adj; exec "$@"' _ "${CMD[@]}" > "$P/fit.log" 2>&1; rc=$?
echo "$(date -u +%FT%TZ) fit segment $SEG: exit $rc"
grep -E '^(positive_fit|consistency|initial cost|prior cost|lambda |gpu_model: bound|gpu_fit: )' "$P/fit.log" | cut -c1-400
tail -6 "$P/fit.log" | cut -c1-250
df -h "$R" > "$P/df_at_end.txt" 2>&1; du -sh "$O" >> "$P/df_at_end.txt" 2>&1
if [ $rc -eq 0 ]; then
    read -r _ st _ _ < <(state "$D" "$SEG")
    [ "$st" = complete ] && { final_link "$SEG" || rc=$?; }
fi
[ $rc -eq 75 ] && echo "stopped before the wall-time limit: resubmit the same command to continue with segment $((SEG + 1))"
exit $rc
