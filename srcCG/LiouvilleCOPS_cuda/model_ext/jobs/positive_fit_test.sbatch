#!/bin/bash
# ==== version-controlled copy: header added 2026-11-07, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/jobs/positive_fit_test.sbatch (on Voyager: /data001/highresolu_vt/cops_output/jobs/positive_fit_test.sbatch)
# SHA256 of the run copy = SHA256 of this file without these header lines: b073dcfa9ada79e1acbe0bd852c9b3102a2df3d3413e27090a8caa2920698141
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: job. Smoke test of the positive-activity fit (lambda 0, max_nfev 4, --check-jacobian); user decisions 2026-11-07.
# RUN: sbatch positive_fit_test.sbatch   (job 122303)
# ==== end of header ====
# 12M fit with the positive H2O activity A = exp(Y a + c_cliff T_cliff), start and prior = the original map (emission-
# weighted, w = area A_orig, LS projection of ln A_orig; user 2026-11-07), original temperature law held, kappa free,
# log-space Tikhonov (lambda), the COPS W10x5 time response in the
# late window (>= 10:33:10 UTC; 1-s grid rows of the lategrid base), grain + crossing 10 (grain start
# 6.96024299190473e-06 Pa, B = 0); GPU patch model_ext/positive_activity.py (night ratio 0.02, flux floor
# 0.02 q_max A_orig; env COPS_POS_*). User decisions 2026-11-07.
# model_ext/positive_fit.py (positive_problem.py around the unchanged grain_fit.py / gpu_fit.py / fit_cops.py).
# Base: runs/shap8_012m_lategrid_dz1.5/base (observations + late 1-s grid; no record store: --precision mixed converts
# the kernels at bind, ~12 min), point table <run>/late_grid_points.txt, cliff_ge_60deg of the flat_earlier bundle
# (same mesh).
# Smoke test only: max_nfev 4, --check-jacobian, lambda 0, into <run>/positive/test-<jobid>/fit.
#   sbatch positive_fit_test.sbatch
# Production fits (one lambda per job, wall-time segments): jobs/positive_fit.sbatch (2026-11-07 review: the TEST=0
# branch of this script was removed; backup positive_fit_test.sbatch.before_review_20261107T080307Z).
# Disk: one fit writes ~6 GB apparent (~3.3 GB on the compressed share): activity_facets.dat, ActivityEnhancement*.
#SBATCH --job-name=cops12m_pos
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gpus-per-node=8
#SBATCH --cpus-per-task=64
#SBATCH --mem=0
#SBATCH --exclusive
#SBATCH --time=02:00:00
#SBATCH --partition=gpu
#SBATCH --account=general
#SBATCH --output=/data001/highresolu_vt/cops_output/jobs/logs/cops12m_pos-%j.out
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
TEST=1; LAMBDA=0; NFEV=4
R=$COPS_OUTPUT_ROOT/runs/shap8_012m_lategrid_dz1.5; X=$COPS_OUTPUT_ROOT/model_ext; PY="$APP/.venv/bin/python -W ignore"
BASE=$R/base; PTS=$R/late_grid_points.txt
BUNDLE=$COPS_OUTPUT_ROOT/runs/shap8_012m_flat_earlier/bundle/restart_templates.npz
NAME=test-$SLURM_JOB_ID; CHECK=(--check-jacobian)
D=$R/positive/$NAME; O=$D/fit; P=$D/provenance-$SLURM_JOB_ID
for f in "$BASE/COPS_model.dat" "$BASE/elements.bin" "$BASE/kernels" "$PTS" "$BUNDLE" "$X/positive_fit.py" \
         "$X/positive_problem.py" "$X/positive_activity.py" "$X/kernels"; do [ -e "$f" ] || { echo "missing: $f" >&2; exit 2; }; done
[ -e "$O" ] && { echo "exists: $O" >&2; exit 2; }
mkdir -p "$P"; cd "$APP" || exit 1
cp -p "$0" "$P/positive_fit_test.sbatch"
cp -p "$X"/positive_fit.py "$X"/positive_problem.py "$X"/positive_activity*.* "$X"/grain_fit.py "$X"/grain_crossing.py \
    "$X"/late_response_eval.py "$X"/late_time_grid.py "$X"/test_positive_fit.py "$P/"
cp -rp "$X/kernels" "$P/kernels"          # positive_activity.py's patched CUDA sources (model_ext/kernels/)
sha256sum "$X"/*.py "$X"/kernels/* "$PTS" > "$P/model_ext.sha256" 2>&1
{ git -C "$AMPS" rev-parse HEAD; git -C "$AMPS" status --short; } > "$P/git.txt" 2>&1
module list > "$P/module_list.txt" 2>&1 || true
nvidia-smi -q > "$P/nvidia_smi_q.txt" 2>&1
{ nvcc --version; gcc --version | head -1; grep -m1 CSPICE /home/vtenishev/SPICE/cspice/doc/version.txt; } > "$P/toolchain.txt" 2>&1
~/.local/bin/uv pip list --python "$APP/.venv/bin/python" > "$P/pip_list.txt" 2>&1
env | grep -E '^(CUDA_[A-Z_]*|CUPY_[A-Z_]*|COPS_[A-Z_]*|TEST|LAMBDA|NFEV|OMP_[A-Z_]*|OPENBLAS_[A-Z_]*|MKL_[A-Z_]*|PYTHONUNBUFFERED|LD_LIBRARY_PATH|PATH|SLURM_JOB_GPUS|CUDA_VISIBLE_DEVICES)=' | sort > "$P/environment.txt"
# commands.sh: every command of this job with all arguments resolved, runnable as is (HOW-TO-RERUN.md)
{ echo '#!/bin/bash'; echo "# job $SLURM_JOB_ID on $(hostname), $(date -u +%FT%TZ); AMPS $(git -C "$AMPS" rev-parse HEAD)"; grep '^export ' "$0"; echo "cd $APP"; } > "$P/commands.sh"
rec() { printf '%q ' "$@" >> "$P/commands.sh"; echo >> "$P/commands.sh"; }
free -g > "$P/free_at_start.txt"; df -h "$R" > "$P/df_at_start.txt" 2>&1
# devices as CUDA numbers them inside the job (on a shared node nvidia-smi lists all GPUs; job 122226 failed with them)
if [ -n "${CUDA_VISIBLE_DEVICES:-}" ]; then NG=$(echo "$CUDA_VISIBLE_DEVICES" | tr , "\n" | wc -l); else NG=$(nvidia-smi -L | wc -l); fi
DEV=$(seq -s, 0 $((NG - 1)))
echo "$(date -u +%FT%TZ) job $SLURM_JOB_ID on $(hostname): $NG GPUs visible ($DEV); fit $NAME (TEST=$TEST, max_nfev $NFEV, lambda $LAMBDA)"
CMD=($PY "$X/positive_fit.py" --point-table "$PTS" --grain-start-uPa 6.96024299190473 --crossing-weight 10 "${CHECK[@]}" \
    --devices "$DEV" --precision mixed --max-nfev "$NFEV" --wall-time slurm --template-start 0 -- "$BASE" \
    --profile all_available --degree 4 --enhancement-degree 4 --lambda "$LAMBDA" --jobs 1 --fix-temperature \
    --tbounds 30 400 --ftol 1e-10 --xtol 1e-10 --gtol 1e-10 --check-tolerance inf --out "$O" \
    --activity-template "$BUNDLE" cliff_ge_60deg)
rec "${CMD[@]}"
bash -c 'echo 1000 > /proc/self/oom_score_adj; exec "$@"' _ "${CMD[@]}" > "$P/fit.log" 2>&1; rc=$?
echo "$(date -u +%FT%TZ) fit: exit $rc"
grep -E '^(positive_fit|consistency|initial cost|lambda |gpu_model: bound|gpu_fit: [0-9])' "$P/fit.log" | cut -c1-300
tail -6 "$P/fit.log" | cut -c1-250
df -h "$R" > "$P/df_at_end.txt" 2>&1; du -sh "$O" >> "$P/df_at_end.txt" 2>&1
exit $rc
