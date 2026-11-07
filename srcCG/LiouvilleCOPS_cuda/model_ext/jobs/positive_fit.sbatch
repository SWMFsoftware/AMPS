#!/bin/bash
# ==== version-controlled copy: header added 2026-11-07, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/jobs/positive_fit.sbatch (on Voyager: /data001/highresolu_vt/cops_output/jobs/positive_fit.sbatch)
# SHA256 of the run copy = SHA256 of this file without these header lines: d4a8fdcef60986c51d6dfddbf13d52f2e99be81a92c81aa65d49c12e3c716c43
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: job. Production positive-activity fit (user decisions 2026-11-07); LAMBDA required, NFEV 4000 by default, CHECK=1 adds the
#   Jacobian check; resume segments seg<NN> (continues an interrupted segment from its restart file).
# RUN: sbatch --export=ALL,LAMBDA=1e-3 positive_fit.sbatch   (jobs 122307, 122310-122312, 122319, 122320)
# ==== end of header ====
# Production 12M fits with the positive H2O activity (user decisions 2026-11-07), one Tikhonov lambda per job:
# A = exp(Y a + c_cliff T_cliff) (cliff_ge_60deg inside the exponent), start and prior = the emission-weighted
# (w = area A_orig) LS projection of ln A_orig, original temperature law held, kappa free, log-space Tikhonov
# D^2 = sum area (s - s0)^2 / sum area, the COPS W10x5 time response in the late window (>= 10:33:10 UTC; 1-s grid rows
# of the lategrid base), grain + crossing 10 (grain start 6.96024299190473e-06 Pa, B = 0); GPU patch
# model_ext/positive_activity.py (night ratio 0.02, H2O flux floor 0.02 q_max A_orig; env COPS_POS_*).
# model_ext/positive_fit.py around the unchanged grain_fit.py / gpu_fit.py / fit_cops.py. The smoke test is
# jobs/positive_fit_test.sbatch (run it first). Base: runs/shap8_012m_lategrid_dz1.5/base (no record store: --precision
# mixed converts the kernels at bind, ~12 min), point table <run>/late_grid_points.txt, the flat_earlier bundle.
#
# Settings: LAMBDA (required, one value per job), NFEV (total evaluations of the fit over all segments, default 4000 as
# the recent 12M cycles; to be confirmed by the user), CHECK (1: --check-jacobian in segment 1, default 1).
# Products: <run>/positive/C30_n<NFEV>_exp_cliff_w10x5_grain_cross10_lam<LAMBDA>/seg<NN>/ (one fit directory per
# wall-time segment) and provenance-<job>/ (commands.sh, logs, environment) per job.
#
# Wall time / resume (the gpu_fit --wall-time slurm convention, as tools/cluster/cops_pipeline.py does it): the fit
# stops at its best evaluated state before the job's end, writes its products and an interrupted fit_restart.npz
# (meta status "interrupted", cumulative "segments") and the job exits 75. Resubmitting the SAME command continues: the
# next segment seg<k+1> starts with gpu_fit --restart-from seg<k>/fit_restart.npz (positive_fit.py accepts only its own
# restart marker), --max-nfev = NFEV - the evaluations of all earlier segments, --segment k+1; the prior is unchanged
# (the projection is recomputed identically). A segment directory without a restart file (a killed or failed attempt) is
# moved aside as seg<k>.incomplete_<UTC> (never deleted) and that segment runs again from seg<k-1>. A segmented fit is
# not identical to one continuous run (each least_squares run restarts the trust region and the scaling).
# A completed fit (restart status "complete", or no evaluations left) is not rerun (exit 0).
#
# Examples (each job takes one whole GPU node; the user picks the scan values after the lambda = 0 test):
#   sbatch --export=ALL,LAMBDA=0 positive_fit.sbatch
#   for L in 0 L1 L2 L3 L4 L5; do sbatch --export=ALL,LAMBDA=$L positive_fit.sbatch; done
#       (L1..L5: 5 log-spaced values to be chosen, e.g. 1e-4 1e-3 1e-2 1e-1 1 -- placeholders, not a decision)
#   resume after exit 75: the same sbatch command again (e.g. sbatch --export=ALL,LAMBDA=0 positive_fit.sbatch)
# Disk: every segment writes a full product set, ~6 GB apparent (~3.3 GB on the compressed share): activity_facets.dat,
# ActivityEnhancement*, activity_exp_facets.npz; check df before a scan (about 34 GB free on 2026-11-07).
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
# the GPU patch's configuration (user 2026-11-07); positive_fit.py sets the same values and refuses different ones
export COPS_POS_ENABLE=1 COPS_POS_NIGHT_RATIO=0.02 COPS_POS_ACTIVITY=exp COPS_POS_FLOOR_RATIO=0.02
[ -n "${LAMBDA:-}" ] || { echo "set LAMBDA (one value per job): sbatch --export=ALL,LAMBDA=0 $0" >&2; exit 2; }
NFEV=${NFEV:-4000}; CHECK=${CHECK:-1}
R=$COPS_OUTPUT_ROOT/runs/shap8_012m_lategrid_dz1.5; X=$COPS_OUTPUT_ROOT/model_ext; PY="$APP/.venv/bin/python -W ignore"
BASE=$R/base; PTS=$R/late_grid_points.txt
BUNDLE=$COPS_OUTPUT_ROOT/runs/shap8_012m_flat_earlier/bundle/restart_templates.npz
NAME=C30_n${NFEV}_exp_cliff_w10x5_grain_cross10_lam${LAMBDA}
D=$R/positive/$NAME; P=$D/provenance-$SLURM_JOB_ID
for f in "$BASE/COPS_model.dat" "$BASE/elements.bin" "$BASE/kernels" "$PTS" "$BUNDLE" "$X/positive_fit.py" \
         "$X/positive_problem.py" "$X/positive_activity.py" "$X/kernels"; do [ -e "$f" ] || { echo "missing: $f" >&2; exit 2; }; done
mkdir -p "$D"
# ---- segment state: the last seg<NN>; its restart's status and the evaluations used so far
read -r LAST LSTATE LUSED < <($PY - "$D" <<'PY'
import glob, json, os, sys
import numpy as np
segs = sorted(int(os.path.basename(p)[3:]) for p in glob.glob(os.path.join(sys.argv[1], "seg[0-9][0-9]")))
if not segs:
    print(0, "none", 0); sys.exit()
k = segs[-1]
rs = os.path.join(sys.argv[1], "seg%02d" % k, "fit_restart.npz")
if not os.path.isfile(rs):
    print(k, "norestart", 0); sys.exit()
with np.load(rs, allow_pickle=False) as z:
    meta = json.loads(str(z["meta_json"]))
used = sum(int(s.get("nfev", 0)) for s in meta.get("segments", []))
print(k, meta.get("status", "unknown"), used)
PY
)
[ -n "${LAST:-}" ] || { echo "cannot read the segment state of $D" >&2; exit 3; }
MOVED=""
if [ "$LSTATE" = norestart ]; then                     # a killed or failed attempt: move aside, run that segment again
    MOVED="$D/seg$(printf %02d "$LAST").incomplete_$(date -u +%Y%m%dT%H%M%SZ)"
    mv "$D/seg$(printf %02d "$LAST")" "$MOVED"
    LAST=$((LAST - 1))
    if [ "$LAST" -ge 1 ]; then
        read -r LAST LSTATE LUSED < <($PY - "$D" "$LAST" <<'PY'
import json, os, sys
import numpy as np
k = int(sys.argv[2])
with np.load(os.path.join(sys.argv[1], "seg%02d" % k, "fit_restart.npz"), allow_pickle=False) as z:
    meta = json.loads(str(z["meta_json"]))
print(k, meta.get("status", "unknown"), sum(int(s.get("nfev", 0)) for s in meta.get("segments", [])))
PY
)
    else LSTATE=none; LUSED=0; fi
fi
case $LSTATE in
    none) SEG=1; RESTART=(); LEFT=$NFEV;;
    interrupted) SEG=$((LAST + 1)); RESTART=(--restart-from "$D/seg$(printf %02d "$LAST")/fit_restart.npz"); LEFT=$((NFEV - LUSED));;
    complete) echo "$(date -u +%FT%TZ) $NAME: segment $LAST is complete ($LUSED evaluations); nothing to do"; exit 0;;
    *) echo "$D/seg$(printf %02d "$LAST")/fit_restart.npz has status '$LSTATE'; not continued" >&2; exit 3;;
esac
if [ "$LEFT" -le 0 ]; then echo "$(date -u +%FT%TZ) $NAME: $LUSED of $NFEV evaluations used; nothing to do"; exit 0; fi
O=$D/seg$(printf %02d "$SEG")
[ -e "$O" ] && { echo "exists: $O" >&2; exit 2; }
CHK=(); [ "$SEG" = 1 ] && [ "$CHECK" = 1 ] && CHK=(--check-jacobian)
# ---- provenance (as grain_late_weights.sbatch)
mkdir -p "$P"; cd "$APP" || exit 1
cp -p "$0" "$P/positive_fit.sbatch"
cp -p "$X"/positive_fit.py "$X"/positive_problem.py "$X"/positive_activity*.* "$X"/grain_fit.py "$X"/grain_crossing.py \
    "$X"/late_response_eval.py "$X"/late_time_grid.py "$X"/test_positive_fit.py "$P/"
cp -rp "$X/kernels" "$P/kernels"          # positive_activity.py's patched CUDA sources (model_ext/kernels/)
sha256sum "$X"/*.py "$X"/kernels/* "$PTS" > "$P/model_ext.sha256" 2>&1
{ git -C "$AMPS" rev-parse HEAD; git -C "$AMPS" status --short; } > "$P/git.txt" 2>&1
module list > "$P/module_list.txt" 2>&1 || true
nvidia-smi -q > "$P/nvidia_smi_q.txt" 2>&1
{ nvcc --version; gcc --version | head -1; grep -m1 CSPICE /home/vtenishev/SPICE/cspice/doc/version.txt; } > "$P/toolchain.txt" 2>&1
~/.local/bin/uv pip list --python "$APP/.venv/bin/python" > "$P/pip_list.txt" 2>&1
env | grep -E '^(CUDA_[A-Z_]*|CUPY_[A-Z_]*|COPS_[A-Z_]*|LAMBDA|NFEV|CHECK|OMP_[A-Z_]*|OPENBLAS_[A-Z_]*|MKL_[A-Z_]*|PYTHONUNBUFFERED|LD_LIBRARY_PATH|PATH|SLURM_JOB_GPUS|CUDA_VISIBLE_DEVICES)=' | sort > "$P/environment.txt"
# commands.sh: every command of this job with all arguments resolved, runnable as is (HOW-TO-RERUN.md)
{ echo '#!/bin/bash'; echo "# job $SLURM_JOB_ID on $(hostname), $(date -u +%FT%TZ); AMPS $(git -C "$AMPS" rev-parse HEAD)"; grep '^export ' "$0"; echo "cd $APP"; } > "$P/commands.sh"
rec() { printf '%q ' "$@" >> "$P/commands.sh"; echo >> "$P/commands.sh"; }
free -g > "$P/free_at_start.txt"; df -h "$R" > "$P/df_at_start.txt" 2>&1
{ echo "segment $SEG; previous segment $LAST ($LSTATE, $LUSED evaluations used); evaluations for this segment $LEFT of $NFEV"
  [ -n "$MOVED" ] && echo "moved aside: $MOVED"; } | tee "$P/segment.txt"
# devices as CUDA numbers them inside the job (on a shared node nvidia-smi lists all GPUs; job 122226 failed with them)
if [ -n "${CUDA_VISIBLE_DEVICES:-}" ]; then NG=$(echo "$CUDA_VISIBLE_DEVICES" | tr , "\n" | wc -l); else NG=$(nvidia-smi -L | wc -l); fi
DEV=$(seq -s, 0 $((NG - 1)))
echo "$(date -u +%FT%TZ) job $SLURM_JOB_ID on $(hostname): $NG GPUs visible ($DEV); fit $NAME segment $SEG (lambda $LAMBDA, $LEFT evaluations)"
CMD=($PY "$X/positive_fit.py" --point-table "$PTS" --grain-start-uPa 6.96024299190473 --crossing-weight 10 "${CHK[@]}" \
    --devices "$DEV" --precision mixed --max-nfev "$LEFT" --segment "$SEG" --wall-time slurm "${RESTART[@]}" \
    --template-start 0 -- "$BASE" \
    --profile all_available --degree 4 --enhancement-degree 4 --lambda "$LAMBDA" --jobs 1 --fix-temperature \
    --tbounds 30 400 --ftol 1e-10 --xtol 1e-10 --gtol 1e-10 --check-tolerance inf --out "$O" \
    --activity-template "$BUNDLE" cliff_ge_60deg)
rec "${CMD[@]}"
bash -c 'echo 1000 > /proc/self/oom_score_adj; exec "$@"' _ "${CMD[@]}" > "$P/fit.log" 2>&1; rc=$?
echo "$(date -u +%FT%TZ) fit segment $SEG: exit $rc"
grep -E '^(positive_fit|consistency|initial cost|prior cost|lambda |gpu_model: bound|gpu_fit: )' "$P/fit.log" | cut -c1-300
tail -6 "$P/fit.log" | cut -c1-250
df -h "$R" > "$P/df_at_end.txt" 2>&1; du -sh "$O" >> "$P/df_at_end.txt" 2>&1
[ $rc -eq 75 ] && echo "stopped before the wall-time limit: resubmit the same command to continue with segment $((SEG + 1))"
exit $rc
