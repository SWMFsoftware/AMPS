#!/bin/bash
# ==== version-controlled copy: header added 2026-11-11, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/jobs/slope_mult_fit.sbatch (on Voyager: /data001/highresolu_vt/cops_output/jobs/slope_mult_fit.sbatch)
# SHA256 of the run copy = SHA256 of this file without these header lines: 10b6db5711ff7cd5819186c9372869b901b401438b65f99a619e0376dd45d1a5
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: job. Slope-multiplier H2O activity fit on the lategrid base (user decisions 2026-11-10: A = max(1e-5 A_orig, m_b
#   L), 7 slope bins, 0-40 deg held at 1; everything else as D2, job 122866): START srccg | d2 | nnls | chain, L_THREADS
#   (chunked), BUMPS (option 3a), HELD, GRAIN_FROM, TESTNAME, TESTMODE, CHAIN; segments, exit 75 = resubmit; one node,
#   local transport.
# RUN: sbatch --nodelist=gpu002 --export=ALL,START=srccg,MODEL_EXT=<snapshot>,L_THREADS=32 slope_mult_fit.sbatch   (this
#   version: 3a fits 124008, 124069, 124078 and tests 123984-124063; the harmonic fit 123786 ran the earlier version
#   .before_blasthreads_*)
# ==== end of header ====
# Slope-multiplier H2O activity fit on the lategrid base (user decisions 2026-11-10; model_ext/slope_mult.py,
# slope_mult_problem.py, driver slope_mult_fit.py, bins bundle slope_bins.py):
#   A = max(1e-5 A_orig, m_b L),  L = the 25 signed degree-<=4 harmonics, one multiplier m_b >= 0 per effective-
#   gravitational-slope bin [40,50) [50,60) [60,70) [70,80) [80,90) [90,180] deg, m_[0,40) = 1 held; F = 1e-5 x A_orig
#   (elements.bin AH2O); flux q(SZA) A unchanged (srcCG night/shadow ratios, no flux floor, CO2 = kappa x run map).
# Everything else = the diagnostic refit D2 (jobs/clip_refit.sbatch VARIANT=D2, job 122866): base
# runs/shap8_012m_lategrid_dz1.5/base (1150 observations + 401 model-only 1-s grid rows, gauge +1.5 m), W10x5 for the
# observations from 10:33:10 UTC (positive_problem's operator through clip_problem.install_core), grain unaveraged at the
# time stamp, Tmax 182.09922 / T75 136.11149 / Tmin 133 K held, kappa free, grain + crossing 10, lambda 0, --tbounds 30 400,
# ftol/xtol/gtol 1e-10, --precision mixed (no record store on this base). NO cliff_ge_60deg, NO ft_c* templates.
# Start (segment 1), START (user decision 2026-11-10, no silent default in the driver: the job passes --start):
#   START=srccg (default, production): the ORIGINAL srcCG harmonic coefficients (the run's map = the fitter's prior x0),
#     every free multiplier 1 (0-40 deg held at 1), kappa 1, grain A and B from D2's final state by name;
#   START=d2: D2's final a, kappa, grain by name, every free multiplier 1 (the first version, tests a/b of 2026-11-10).
#   A resume passes the same --start; the driver refuses it if the fit started with another one.
# ONE node, local transport (8 GPUs; spawned workers get the GPU patch from COPS_SLOPEMUL_FROZEN_CONFIG). The multi-node
# launcher (MPI) is jobs/mpi_fit.sbatch of the second agent; the driver installs its patch in every process it runs in.
#
# Settings (sbatch --export=ALL,VAR=value,...; values must not contain commas):
#   NFEV          evaluations of the fit over all segments (default 4000)
#   CHECK         1: segment 1 runs grain_fit --check-jacobian and slope_mult --check-full-jacobian (default 1)
#   CHECKPARAMS   parameters of the full Jacobian check, ':'-separated (default: the driver's: the 6 multipliers, a00,
#                 a02, a07, a12, kappa, grain_A_uPa)
#   START         srccg (default) or d2: the optimizer start (see above)
#   START_RESTART the --start-restart file (default D2's final fit_restart.npz): grain (srccg) or a, kappa, grain (d2)
#   BUNDLE        the slope-bin bundle (default runs/shap8_012m_lategrid_dz1.5/slope_mult/bundle/slope_bins.npz)
#   TESTNAME      test runs: products in <run>/slope_mult_tests/<TESTNAME>/seg<NN>
#   TESTMODE      floor0 (needs TESTNAME): GPU check (a): floor ratio 0, all multipliers 1, stop after the start check;
#                 then the independent evaluation of the plain clipped harmonic model (test_slope_mult.py --gpu-independent,
#                 an UNPATCHED process) is compared with the start check
#   L_THREADS     N (optional; user 2026-11-10 "Bitwise chunked version"): the FITTER process computes L = Y a with N
#                 Python threads, each Y[chunk] @ a with OpenBLAS at 1 thread (4096-row-aligned chunks), bitwise the
#                 1-thread L (slope_mult_fit.py --l-threads N). Unset: unchanged (one call, 1 thread). OpenBLAS stays at
#                 1 thread everywhere (OMP/OPENBLAS/MKL_NUM_THREADS=1).
#   L_MODE        chunked (default with L_THREADS) or openblas (explicit, non-default, NOT bitwise: threaded OpenBLAS for
#                 that call). BLAS_THREADS (the earlier test setting) is refused.
#   BUMPS         the bump bundle (slope_bumps.py; option 3a, user 2026-11-10): L = A_scale sum_k w_k phi_k, a00..a24
#                 held at 0 (slope_mult_fit.py --bump-bundle); START must then be nnls or chain. Unset: harmonics.
#   HELD          the previous fine tuning's finetune_state.npz: its held_next is held in A = max(F, m L exp(h))
#                 (--held; every segment); START=chain needs it
#   GRAIN_FROM    the previous fine tuning's finetune_state.npz: grain A, B of the start (--grain-from; segment 1 of
#                 START=chain)
#   FITDIR        product directory of the fit (segments FITDIR/seg<NN>, FITDIR/final); default: runs/.../slope_mult/fit
#                 (harmonics) or runs/.../slope_bumps/fit (BUMPS); TESTNAME overrides both
#   MODEL_EXT     the model_ext directory the job copies and runs (default $COPS_OUTPUT_ROOT/model_ext); production runs
#                 use a frozen snapshot MODEL_EXT=$COPS_OUTPUT_ROOT/model_ext_frozen/slopemul_<UTC> (never edited)
#   CHAIN         1: a segment directory without restart file stops the job (exit 3) instead of being moved aside
# Products: <run>/slope_mult/fit/seg<NN> (+ link final -> seg<NN> when complete); provenance-<job>/ next to them
# (commands.sh with every argument resolved, logs, environment, git, nvidia-smi -q, uv pip list, the model_ext code it
# ran). The fit runs from that provenance copy of model_ext.
# Wall time / resume (as jobs/clip_refit.sbatch): the fit stops at its best evaluated state before the job's end,
# writes its products and an interrupted fit_restart.npz, exit 75; resubmitting the SAME command continues with segment
# k+1 (gpu_fit --restart-from seg<k>/fit_restart.npz; the driver checks the marker and the start check must reproduce
# seg<k>'s cost to 1e-8), --max-nfev = NFEV - the evaluations of the earlier segments. A segment directory without
# restart file is moved aside as seg<k>.incomplete_<UTC> and run again. A completed fit is not rerun (exit 0).
# Disk: ~3.5 GB per segment (2 GB fitter activity tables -- not valid for this model, see SLOPE_MULTIPLIER_MODEL.txt --,
# 0.44 GB restart, ~0.3 GB slope_mult_facets.npz).
#SBATCH --job-name=cops12m_slopemul
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
#SBATCH --output=/data001/highresolu_vt/cops_output/jobs/logs/cops12m_slopemul-%j.out
set -uo pipefail
export AMPS=/data001/highresolu_vt/AMPS APP=/data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda
export T=/data001/highresolu_vt/cops_cluster_transfer_20261026 COPS_OUTPUT_ROOT=/data001/highresolu_vt/cops_output
export COPS_DATA_ROOT=$T/rosetta_data_end_of_mission_cops
export CUDA_PATH=/usr/local/cuda-12.9 CUDA_HOME=/usr/local/cuda-12.9 CUPY_CACHE_DIR=$COPS_OUTPUT_ROOT/.cupy
export PATH=$CUDA_HOME/bin:/usr/mpi/gcc/openmpi-4.1.7rc1/bin:$PATH
export LD_LIBRARY_PATH=$CUDA_HOME/lib64${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1
# this model patches the evaluator itself: the positive-activity patch must be off; no inherited frozen slope config
UNSET=$(env | grep -o -E '^(COPS_POS_[A-Z_]*|COPS_SLOPEMUL_[A-Z_]*)' | sort | paste -sd' ' -)
[ -n "$UNSET" ] && unset $UNSET
NFEV=${NFEV:-4000}; CHECK=${CHECK:-1}; CHECKPARAMS=${CHECKPARAMS:-}; TESTNAME=${TESTNAME:-}; TESTMODE=${TESTMODE:-}
CHAIN=${CHAIN:-0}; START=${START:-srccg}; L_THREADS=${L_THREADS:-}; L_MODE=${L_MODE:-}
BUMPS=${BUMPS:-}; HELD=${HELD:-}; GRAIN_FROM=${GRAIN_FROM:-}; FITDIR=${FITDIR:-}
[ -n "$BUMPS" ] && BUMPS=$(readlink -f "$BUMPS"); [ -n "$HELD" ] && HELD=$(readlink -f "$HELD")
[ -n "$GRAIN_FROM" ] && GRAIN_FROM=$(readlink -f "$GRAIN_FROM")
if [ -n "${BLAS_THREADS:-}" ]; then
    echo "BLAS_THREADS was replaced (user 2026-11-10 'Bitwise chunked version'): use L_THREADS=N (chunked, bitwise) or L_THREADS=N L_MODE=openblas" >&2; exit 2
fi
[ -n "$L_MODE" ] && [ -z "$L_THREADS" ] && { echo "L_MODE needs L_THREADS" >&2; exit 2; }
BOPT=(); [ -n "$L_THREADS" ] && BOPT=(--l-threads "$L_THREADS" --l-mode "${L_MODE:-chunked}")
if [ -n "$BUMPS" ]; then case $START in nnls|chain) ;; *) echo "BUMPS: START must be nnls or chain" >&2; exit 2;; esac
else case $START in srccg|d2|chain) ;; *) echo "START must be srccg, d2 or chain" >&2; exit 2;; esac; fi
if [ "$START" = chain ] && { [ -z "$HELD" ] || [ -z "$GRAIN_FROM" ]; }; then echo "START=chain needs HELD and GRAIN_FROM" >&2; exit 2; fi
export MODEL_EXT=$(readlink -f "${MODEL_EXT:-$COPS_OUTPUT_ROOT/model_ext}")
R=$COPS_OUTPUT_ROOT/runs/shap8_012m_lategrid_dz1.5
X=$MODEL_EXT; PY="$APP/.venv/bin/python -W ignore"
BASE=$R/base; PTS=$R/late_grid_points.txt
BUNDLE=$(readlink -f "${BUNDLE:-$R/slope_mult/bundle/slope_bins.npz}")
START_RESTART=$(readlink -f "${START_RESTART:-$R/clip_refit/D2/final/fit_restart.npz}")
FIX=(--fix Tmax=182.09922 --fix T75=136.11149 --fix Tmin=133.0 --fix c_slope_00_40=1)
case $TESTMODE in
    "") TOPT=();;
    floor0) [ -n "$TESTNAME" ] || { echo "TESTMODE=floor0 needs TESTNAME" >&2; exit 2; }
            TOPT=(--floor-ratio 0 --test-floor-zero --stop-after-checks);;
    *) echo "TESTMODE must be empty or floor0" >&2; exit 2;;
esac
if [ -n "$TESTNAME" ]; then NAME=slope_mult_test_$TESTNAME; D=$R/slope_mult_tests/$TESTNAME
elif [ -n "$FITDIR" ]; then D=$(readlink -m "$FITDIR"); NAME=slope_fit_$(basename "$(dirname "$D")")_$(basename "$D")
elif [ -n "$BUMPS" ]; then NAME=slope_bumps_fit; D=$R/slope_bumps/fit
else NAME=slope_mult_fit; D=$R/slope_mult/fit; fi
MOPT=(); [ -n "$BUMPS" ] && MOPT+=(--bump-bundle "$BUMPS"); [ -n "$HELD" ] && MOPT+=(--held "$HELD")
P=$D/provenance-$SLURM_JOB_ID
FILES=(slope_mult_fit.py slope_mult.py slope_mult_problem.py slope_bins.py slope_bumps.py slope_finetune.py test_slope_mult.py
       positive_finetune.py clip_problem.py clip_refit.py
       grain_fit_fixed.py fixed_parameters.py grain_fit.py grain_crossing.py positive_problem.py late_response_eval.py
       test_clip_refit.py test_positive_fit.py)
for f in "$BASE/COPS_model.dat" "$BASE/elements.bin" "$BASE/kernels" "$PTS" "$START_RESTART" "$BUNDLE" ${BUMPS:+"$BUMPS"} \
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
cp -p "$0" "$P/slope_mult_fit.sbatch"
for f in "${FILES[@]}"; do cp -p "$X/$f" "$P/"; done
sha256sum "$X"/*.py "$PTS" "$START_RESTART" "$BUNDLE" > "$P/model_ext.sha256" 2>&1
(cd "$P" && sha256sum "${FILES[@]}") > "$P/run_copy.sha256" 2>&1
{ git -C "$AMPS" rev-parse HEAD; git -C "$AMPS" status --short; } > "$P/git.txt" 2>&1
module list > "$P/module_list.txt" 2>&1 || true
nvidia-smi -q > "$P/nvidia_smi_q.txt" 2>&1
nvidia-smi -L > "$P/nvidia_smi_L.txt" 2>&1
{ nvcc --version; gcc --version | head -1; grep -m1 CSPICE /home/vtenishev/SPICE/cspice/doc/version.txt; } > "$P/toolchain.txt" 2>&1
~/.local/bin/uv pip list --python "$APP/.venv/bin/python" > "$P/pip_list.txt" 2>&1
env | grep -E '^(CUDA_[A-Z_]*|CUPY_[A-Z_]*|COPS_[A-Z_]*|NFEV|CHECK|CHECKPARAMS|TESTNAME|TESTMODE|START|START_RESTART|L_THREADS|L_MODE|BUMPS|HELD|GRAIN_FROM|FITDIR|BUNDLE|CHAIN|MODEL_EXT|OMP_[A-Z_]*|OPENBLAS_[A-Z_]*|MKL_[A-Z_]*|PYTHONUNBUFFERED|LD_LIBRARY_PATH|PATH|SLURM_JOB_GPUS|SLURM_JOB_NODELIST|SLURM_JOB_NUM_NODES|CUDA_VISIBLE_DEVICES)=' | sort > "$P/environment.txt"
# commands.sh: every command of this job with all arguments resolved, runnable as is (HOW-TO-RERUN.md)
{ echo '#!/bin/bash'; echo "# job $SLURM_JOB_ID on $(hostname), $(date -u +%FT%TZ); AMPS $(git -C "$AMPS" rev-parse HEAD)"
  grep '^export ' "$0"; echo "cd $APP"; } > "$P/commands.sh"
rec() { printf '%q ' "$@" >> "$P/commands.sh"; echo >> "$P/commands.sh"; }
free -g > "$P/free_at_start.txt"; df -h "$R" > "$P/df_at_start.txt" 2>&1
{ echo "fit $NAME, directory $D; start $START; start restart (segment 1) $START_RESTART; slope bundle $BUNDLE"
  echo "bump bundle ${BUMPS:-none (harmonics)}; held ${HELD:-none}; grain from (segment 1) ${GRAIN_FROM:-the start restart}"
  if [ -n "$L_THREADS" ]; then echo "L = Y a of the fitter process: ${L_MODE:-chunked}, $L_THREADS threads; OPENBLAS_NUM_THREADS=$OPENBLAS_NUM_THREADS"
  else echo "L = Y a of the fitter process: one call, 1 thread (unchanged); OPENBLAS_NUM_THREADS=$OPENBLAS_NUM_THREADS"; fi
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
CMD=($PY "$P/slope_mult_fit.py" --point-table "$PTS" --slope-bundle "$BUNDLE" "${MOPT[@]}" "${SOPT[@]}" "${BOPT[@]}" "${CHK[@]}" "${CHKP[@]}" \
    "${TOPT[@]}" "${FIX[@]}" --crossing-weight 10 "${GCHK[@]}" --devices "$DEV" --precision mixed --max-nfev "$LEFT" \
    --segment "$SEG" --wall-time slurm --template-start 0 -- "$BASE" \
    --profile all_available --degree 4 --enhancement-degree 4 --lambda 0 --jobs 1 --tbounds 30 400 \
    --ftol 1e-10 --xtol 1e-10 --gtol 1e-10 --check-tolerance 1e-6 --out "$O")
rec bash -c 'echo 1000 > /proc/self/oom_score_adj; exec "$@"' _ "${CMD[@]}"
bash -c 'echo 1000 > /proc/self/oom_score_adj; exec "$@"' _ "${CMD[@]}" > "$P/fit.log" 2>&1; rc=$?
echo "$(date -u +%FT%TZ) $NAME segment $SEG: exit $rc"
grep -E '^(slope_mult|clip_refit|consistency|initial cost|lambda |gpu_model: bound|gpu_fit: |grain_fit: |fit data|selection profile)' "$P/fit.log" | cut -c1-600
tail -6 "$P/fit.log" | cut -c1-250
if [ "$TESTMODE" = floor0 ] && [ $rc -eq 0 ]; then
    ICMD=($PY "$P/test_slope_mult.py" --gpu-independent "$O" --devices "$DEV")
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
