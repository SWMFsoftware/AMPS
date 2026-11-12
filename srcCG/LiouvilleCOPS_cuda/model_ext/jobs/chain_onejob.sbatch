#!/bin/bash
# ==== version-controlled copy: header added 2026-11-11, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/jobs/chain_onejob.sbatch (on Voyager: /data001/highresolu_vt/cops_output/jobs/chain_onejob.sbatch)
# SHA256 of the run copy = SHA256 of this file without these header lines: 74a6a72f8dd4e911031c8a3ded78b25831c768b7eebbc009a60178838f923c85
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: job. ONE job for a whole fit / fine-tuning chain (user decision 2026-11-10 "yes, implement (b)"): environment,
#   job-level provenance and code copy, then model_ext/chain_onejob.py; KIND=positive | slope | nonneg (v6, 2026-11-11),
#   CYCLES, NFEV, CHECK, FT_*, L_THREADS, MODEL_EXT, CHAIN_CODE, CHAIN_ROOT; exit 75 = resubmit (completed stages are
#   skipped).
# RUN: through jobs/chain_onejob.sh (dry run by default). This version: nonneg tests 124143, 124148, production 124162
#   (softplus) and 124165 (squared); the positive-kind tests 123990, 124091-124095 ran the version before KIND
#   (.before_kinds_*).
# ==== end of header ====
# ONE job for a whole positive-activity fit / fine-tuning chain of ONE lambda (user decision 2026-11-10, "yes,
# implement (b)"): one driver process (model_ext/chain_onejob.py) binds the 12M lategrid base's records on the node's
# GPUs once and runs every stage on that backend:
#     c1 fit -> c1 ft -> c2 fit -> c2 ft -> ...  (CYCLES cycles)
# with the products, directory layout and results of the split chain (jobs/positive_cycles.sh = one job per stage:
# jobs/positive_fit.sbatch CHAIN=1, jobs/positive_finetune.sbatch). Each stage's logic (segments, resume, move-aside,
# START_RESTART / HELD / GRAIN_FROM, final links, product checks) is ported in chain_onejob.py; each stage runs the
# command its split job script would run, in this process (see the driver's docstring for how every stage starts from
# a fresh process's module state).
#
# Settings (sbatch --export=ALL,VAR=value,...; values without commas):
#   LAMBDA          required, verbatim (directory name lam<LAMBDA>)
#   RATIO           required: 0.02 | 0.002 | 1e-5 (the user's decisions; 1e-5 includes the CO2 flux floor)
#   SOURCE          required: the cycle-1 fit's START_RESTART (a positive fit's fit_restart.npz, status complete, this
#                   lambda); jobs/chain_onejob.sh picks positive_cycles.sh's default
#   CYCLES          number of cycles (default 3)
#   NFEV, CHECK     each fit's total evaluations (default 4000) and --check-jacobian in its segment 1 (default 1)
#   C1_START_SET    NAME=VALUE[:NAME=VALUE...] for the cycle-1 fit only (positive_fit.py --start-set; e.g. kappa=1)
#   FT_MAXNFEV      required: each fine tuning's MAXNFEV (user 2026-11-08: 500)
#   FT_FACETS       required: late_only | all_seen; FT_WINDOW_MIN with late_only (user 2026-11-08: 5)
#   FT_LSMR_MAXITER (50), FT_CHECKPOINT_MIN (30), FT_VALIDATE (1), FT_FLOOR_RATIO (default: the fit's ratio)
#   MODEL_EXT       the model_ext the positive code is copied from (default $COPS_OUTPUT_ROOT/model_ext; production:
#                   a frozen snapshot model_ext_frozen/<name>)
#   CHAIN_CODE      directory holding chain_onejob.py (default $COPS_OUTPUT_ROOT/model_ext)
#   CHAIN_ROOT      the chain's directory (default <run>/positive_r<RATIO>/lam<LAMBDA>, as positive_cycles.sh);
#                   tests: a directory under $COPS_OUTPUT_ROOT/tests
#   MIN_FIT_MINUTES, MIN_FT_MINUTES (25), BIND_MINUTES (4): a stage starts only with that much time left
# KIND (2026-11-11, v6): positive (default; everything above) | slope | nonneg. slope / nonneg: the slope-multiplier
#   families (their own GPU patch, slope_mult / nonneg_mult; COPS_POS_* / COPS_SLOPEMUL_* / COPS_NONNEG_* are unset):
#   SLOPE_MODEL     required: bumps (KIND slope: slope_mult_fit.py / slope_finetune.py, the split chain jobs/slope_cycles.sh)
#                   | softplus | squared (KIND nonneg: nonneg_fit.py / nonneg_finetune.py, jobs/nonneg_*.sbatch)
#   CHAIN_ROOT      required (no default for these kinds)
#   SOURCE          the cycle-1 --start-restart (default D2's final fit_restart.npz)
#   SLOPE_START     the cycle-1 --start (default: bumps nnls, softplus d2, squared lsq); cycles > 1: chain
#   START_FILE      squared with --start lsq: nonneg_start.py's start file
#   SLOPE_BUNDLE    the slope-bin bundle (default runs/.../slope_mult/bundle/slope_bins.npz); BUMPS (bumps: bumps32.npz);
#   FT_TEMPLATES_FROM (nonneg: D2's final fit_manifest.json, the ft_c1..ft_c6 templates)
#   L_THREADS, L_MODE (chunked): the fits' and fine tunings' chunked L threads (user 2026-11-10/11: 32); CHECKPARAMS
#   FT_FACETS must be late_only; LAMBDA / RATIO / C1_START_SET / FT_FLOOR_RATIO are positive-only
# Products: CHAIN_ROOT/c<k>/fit/seg<NN>, c<k>/fit/final, c<k>/ft/seg<NN>, c<k>/ft/final + product links; per stage
# CHAIN_ROOT/c<k>/{fit,ft}/provenance-<job>-c<k>{fit,ft}/; per job CHAIN_ROOT/provenance-<job>/ (this script, the code
# copy that ran (code/), chain.log = the full output, chain_summary.json, git/module/nvidia-smi/toolchain/pip records);
# CHAIN_ROOT/chain_onejob_events.jsonl (stage start/end per job).
# Wall time: every stage stops at its best state before the job's end (--wall-time slurm) and the job exits 75; a
# stage is not started with less than MIN_*_MINUTES left (exit 75). Resubmitting the SAME command skips the completed
# stages (checked) and continues the current one (next segment). Exit codes: 0 complete, 75 resubmit, else the failing
# stage's code (2 usage, 3 inconsistent state, 4 product check, 1 error). Nothing is deleted.
# Example (dry run first: jobs/chain_onejob.sh):
#   sbatch --nodelist=gpu002 --time=08:00:00 --export=ALL,LAMBDA=1e-4,RATIO=1e-5,SOURCE=<src>,FT_MAXNFEV=500,\
#       FT_FACETS=late_only,FT_WINDOW_MIN=5,MODEL_EXT=<frozen> chain_onejob.sbatch
#SBATCH --job-name=cops12m_chain
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gpus-per-node=8
#SBATCH --cpus-per-task=64
#SBATCH --mem=0
#SBATCH --exclusive
#SBATCH --time=08:00:00
#SBATCH --partition=gpu
#SBATCH --account=general
#SBATCH --exclude=gpu003
#SBATCH --output=/data001/highresolu_vt/cops_output/jobs/logs/cops12m_chain-%j.out
set -uo pipefail
export AMPS=/data001/highresolu_vt/AMPS APP=/data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda
export T=/data001/highresolu_vt/cops_cluster_transfer_20261026 COPS_OUTPUT_ROOT=/data001/highresolu_vt/cops_output
export COPS_DATA_ROOT=$T/rosetta_data_end_of_mission_cops
export CUDA_PATH=/usr/local/cuda-12.9 CUDA_HOME=/usr/local/cuda-12.9 CUPY_CACHE_DIR=$COPS_OUTPUT_ROOT/.cupy
export PATH=$CUDA_HOME/bin:/usr/mpi/gcc/openmpi-4.1.7rc1/bin:$PATH
export LD_LIBRARY_PATH=$CUDA_HOME/lib64${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1
die() { echo "chain_onejob.sbatch: $*" >&2; exit 2; }
KIND=${KIND:-positive}
case $KIND in positive|slope|nonneg) ;; *) die "KIND=$KIND: positive (default), slope or nonneg";; esac
SOURCE=${SOURCE:-}; CYCLES=${CYCLES:-3}; NFEV=${NFEV:-4000}; CHECK=${CHECK:-1}
FT_MAXNFEV=${FT_MAXNFEV:-}; FT_FACETS=${FT_FACETS:-}; FT_WINDOW_MIN=${FT_WINDOW_MIN:-}
FT_LSMR_MAXITER=${FT_LSMR_MAXITER:-50}; FT_CHECKPOINT_MIN=${FT_CHECKPOINT_MIN:-30}; FT_VALIDATE=${FT_VALIDATE:-1}
FT_FLOOR_RATIO=${FT_FLOOR_RATIO:-}; MIN_FIT_MINUTES=${MIN_FIT_MINUTES:-25}; MIN_FT_MINUTES=${MIN_FT_MINUTES:-25}
BIND_MINUTES=${BIND_MINUTES:-4}
export MODEL_EXT=$(readlink -f "${MODEL_EXT:-$COPS_OUTPUT_ROOT/model_ext}")
CHAIN_CODE=$(readlink -f "${CHAIN_CODE:-$COPS_OUTPUT_ROOT/model_ext}")
R=$COPS_OUTPUT_ROOT/runs/shap8_012m_lategrid_dz1.5; X=$MODEL_EXT; PY="$APP/.venv/bin/python -W ignore"
PTS=$R/late_grid_points.txt
[ -n "$FT_MAXNFEV" ] || die "set FT_MAXNFEV (the fine tuning's MAXNFEV; no default here)"
case $FT_FACETS in
    late_only) [ -n "$FT_WINDOW_MIN" ] || die "FT_FACETS=late_only needs FT_WINDOW_MIN";;
    all_seen) [ "$KIND" = positive ] || die "FT_FACETS=all_seen: positive kind only"
              [ -z "$FT_WINDOW_MIN" ] || die "FT_WINDOW_MIN applies to late_only only";;
    *) die "set FT_FACETS=late_only (with FT_WINDOW_MIN) or all_seen";;
esac
if [ "$KIND" = positive ]; then
[ -n "${LAMBDA:-}" ] || die "set LAMBDA"
RATIO=${RATIO:-}; C1_START_SET=${C1_START_SET:-}
case $RATIO in 0.02|0.002|1e-5) ;; *) die "RATIO=$RATIO: the user decided 0.02, 0.002 and 1e-5 (write it exactly so)";; esac
[ -n "$SOURCE" ] || die "set SOURCE (the cycle-1 fit's START_RESTART)"
case $C1_START_SET in *[!A-Za-z0-9_=.:+-]*) die "C1_START_SET=$C1_START_SET: NAME=VALUE[:NAME=VALUE...] only";; esac
# the GPU patch's configuration (as jobs/positive_fit.sbatch; the fine tunings check it against their fit)
export COPS_POS_ENABLE=1 COPS_POS_NIGHT_RATIO=$RATIO COPS_POS_ACTIVITY=exp COPS_POS_FLOOR_RATIO=$RATIO
if [ "$RATIO" = 1e-5 ]; then export COPS_POS_FLOOR_CO2=1; else unset COPS_POS_FLOOR_CO2; fi
BUNDLE=$COPS_OUTPUT_ROOT/runs/shap8_012m_flat_earlier/bundle/restart_templates.npz
CHAIN_ROOT=$(readlink -m "${CHAIN_ROOT:-$R/positive_r$RATIO/lam$LAMBDA}")
KFILES=("$X/positive_fit.py" "$X/positive_finetune.py" "$X/positive_problem.py" "$X/positive_activity.py" "$X/kernels" "$BUNDLE")
else
# slope (bumps, slope_mult_fit / slope_finetune) or nonneg (softplus / squared, nonneg_fit / nonneg_finetune): their own
# GPU patch (slope_mult / nonneg_mult); the positive patch must be off, no inherited frozen configuration (as
# jobs/slope_mult_fit.sbatch)
UNSET=$(env | grep -o -E '^(COPS_POS_[A-Z_]*|COPS_SLOPEMUL_[A-Z_]*|COPS_NONNEG_[A-Z_]*)' | sort | paste -sd' ' -)
[ -n "$UNSET" ] && unset $UNSET
SLOPE_MODEL=${SLOPE_MODEL:-}; SLOPE_START=${SLOPE_START:-}; START_FILE=${START_FILE:-}; L_THREADS=${L_THREADS:-}
L_MODE=${L_MODE:-chunked}; CHECKPARAMS=${CHECKPARAMS:-}
SOURCE=${SOURCE:-$R/clip_refit/D2/final/fit_restart.npz}
SLOPE_BUNDLE=$(readlink -f "${SLOPE_BUNDLE:-$R/slope_mult/bundle/slope_bins.npz}")
FT_TEMPLATES_FROM=$(readlink -f "${FT_TEMPLATES_FROM:-$R/clip_refit/D2/final/fit_manifest.json}")
BUMPS=$(readlink -f "${BUMPS:-$R/slope_bumps/bundle/bumps32.npz}")
[ "$FT_FACETS" = late_only ] || die "KIND=$KIND: FT_FACETS=late_only"
[ -n "${CHAIN_ROOT:-}" ] || die "KIND=$KIND: set CHAIN_ROOT"
CHAIN_ROOT=$(readlink -m "$CHAIN_ROOT")
case $KIND/$SLOPE_MODEL in
    slope/bumps) KFILES=("$X/slope_mult_fit.py" "$X/slope_finetune.py" "$X/slope_mult.py" "$BUMPS");;
    nonneg/softplus) KFILES=("$X/nonneg_fit.py" "$X/nonneg_finetune.py" "$X/nonneg_mult.py" "$FT_TEMPLATES_FROM");;
    nonneg/squared) KFILES=("$X/nonneg_fit.py" "$X/nonneg_finetune.py" "$X/nonneg_mult.py" "$FT_TEMPLATES_FROM")
                    [ -n "$START_FILE" ] || [ -n "$SLOPE_START" ] || die "SLOPE_MODEL=squared: set START_FILE (nonneg_start.py)"
                    [ -n "$START_FILE" ] && START_FILE=$(readlink -f "$START_FILE") && KFILES+=("$START_FILE");;
    *) die "KIND=$KIND SLOPE_MODEL=$SLOPE_MODEL: slope/bumps, nonneg/softplus or nonneg/squared";;
esac
KFILES+=("$SLOPE_BUNDLE")
fi
case $SOURCE in /*) ;; *) die "SOURCE must be an absolute path";; esac
case $CHAIN_ROOT in "$R"/*|"$COPS_OUTPUT_ROOT"/tests/*) ;; *) die "CHAIN_ROOT=$CHAIN_ROOT: under $R or $COPS_OUTPUT_ROOT/tests";; esac
for f in "$R/base/COPS_model.dat" "$R/base/elements.bin" "$PTS" "$SOURCE" "$CHAIN_CODE/chain_onejob.py" "${KFILES[@]}"; do
    [ -e "$f" ] || die "missing: $f"; done
mkdir -p "$CHAIN_ROOT"; P=$CHAIN_ROOT/provenance-$SLURM_JOB_ID; mkdir -p "$P/code"; cd "$APP" || exit 1
# ---- job-level provenance and the code copy the driver runs (later edits of MODEL_EXT / CHAIN_CODE cannot reach it)
cp -p "$0" "$P/chain_onejob.sbatch"
if [ "$KIND" = positive ]; then
cp -p "$X"/positive_fit.py "$X"/positive_problem.py "$X"/positive_activity*.* "$X"/grain_fit.py "$X"/grain_crossing.py \
    "$X"/late_response_eval.py "$X"/late_time_grid.py "$X"/test_positive_fit.py "$X"/positive_finetune.py \
    "$X"/test_positive_finetune.py "$P/code/" || exit 1
else
cp -p "$X"/*.py "$P/code/" || exit 1                 # a frozen snapshot: every module of it
fi
[ -d "$X/kernels" ] && { cp -rp "$X/kernels" "$P/code/kernels" || exit 1; }
cp -p "$CHAIN_CODE/chain_onejob.py" "$P/code/" || exit 1
[ -f "$CHAIN_CODE/test_chain_onejob.py" ] && cp -p "$CHAIN_CODE/test_chain_onejob.py" "$P/code/"
{ sha256sum "$X"/*.py $( [ -d "$X/kernels" ] && echo "$X"/kernels/*) "$PTS"; sha256sum "$CHAIN_CODE/chain_onejob.py"; } > "$P/model_ext.sha256" 2>&1
(cd "$P/code" && sha256sum *.py $( [ -d kernels ] && echo kernels/*)) > "$P/code.sha256" 2>&1
{ git -C "$AMPS" rev-parse HEAD; git -C "$AMPS" status --short; } > "$P/git.txt" 2>&1
module list > "$P/module_list.txt" 2>&1 || true
nvidia-smi -q > "$P/nvidia_smi_q.txt" 2>&1
{ nvcc --version; gcc --version | head -1; grep -m1 CSPICE /home/vtenishev/SPICE/cspice/doc/version.txt; } > "$P/toolchain.txt" 2>&1
~/.local/bin/uv pip list --python "$APP/.venv/bin/python" > "$P/pip_list.txt" 2>&1
env | grep -E '^(CUDA_[A-Z_]*|CUPY_[A-Z_]*|COPS_[A-Z0-9_]*|KIND|LAMBDA|RATIO|SOURCE|CYCLES|NFEV|CHECK|CHECKPARAMS|C1_START_SET|FT_[A-Z_]*|SLOPE_[A-Z_]*|START_FILE|L_THREADS|L_MODE|BUMPS|MODEL_EXT|CHAIN_[A-Z_]*|MIN_[A-Z_]*|BIND_MINUTES|OMP_[A-Z_]*|OPENBLAS_[A-Z_]*|MKL_[A-Z_]*|PYTHONUNBUFFERED|LD_LIBRARY_PATH|PATH|SLURM_JOB_GPUS|CUDA_VISIBLE_DEVICES)=' | sort > "$P/environment.txt"
free -g > "$P/free_at_start.txt"; df -h "$CHAIN_ROOT" > "$P/df_at_start.txt" 2>&1
if [ -n "${CUDA_VISIBLE_DEVICES:-}" ]; then NG=$(echo "$CUDA_VISIBLE_DEVICES" | tr , "\n" | wc -l); else NG=$(nvidia-smi -L | wc -l); fi
DEV=$(seq -s, 0 $((NG - 1)))
CMD=($PY "$P/code/chain_onejob.py" --kind "$KIND" --chain-root "$CHAIN_ROOT" --cycles "$CYCLES"
     --source "$SOURCE" --nfev "$NFEV" --check "$CHECK" --ft-maxnfev "$FT_MAXNFEV" --ft-facets "$FT_FACETS"
     --ft-lsmr-maxiter "$FT_LSMR_MAXITER" --ft-checkpoint-min "$FT_CHECKPOINT_MIN" --ft-validate "$FT_VALIDATE"
     --model-ext "$X" --job-provenance "$P" --job-script "$P/chain_onejob.sbatch" --devices "$DEV"
     --point-table "$PTS" --run "$R"
     --min-fit-minutes "$MIN_FIT_MINUTES" --min-ft-minutes "$MIN_FT_MINUTES" --bind-minutes "$BIND_MINUTES")
[ -n "$FT_WINDOW_MIN" ] && CMD+=(--ft-window-min "$FT_WINDOW_MIN")
if [ "$KIND" = positive ]; then
    CMD+=(--lambda "$LAMBDA" --ratio "$RATIO" --bundle "$BUNDLE")
    [ -n "$FT_FLOOR_RATIO" ] && CMD+=(--ft-floor-ratio "$FT_FLOOR_RATIO")
    [ -n "$C1_START_SET" ] && CMD+=(--c1-start-set "$C1_START_SET")
else
    CMD+=(--slope-model "$SLOPE_MODEL" --slope-bundle "$SLOPE_BUNDLE")
    [ "$KIND" = slope ] && CMD+=(--bumps "$BUMPS")
    [ "$KIND" = nonneg ] && CMD+=(--ft-templates-from "$FT_TEMPLATES_FROM")
    [ -n "$SLOPE_START" ] && CMD+=(--slope-start "$SLOPE_START")
    [ -n "$START_FILE" ] && CMD+=(--start-file "$START_FILE")
    [ -n "$L_THREADS" ] && CMD+=(--l-threads "$L_THREADS" --l-mode "$L_MODE")
    [ -n "$CHECKPARAMS" ] && CMD+=(--checkparams "$CHECKPARAMS")
fi
{ echo '#!/bin/bash'; echo "# job $SLURM_JOB_ID on $(hostname), $(date -u +%FT%TZ); AMPS $(git -C "$AMPS" rev-parse HEAD)"
  grep '^export ' "$0" | grep -v COPS_POS_
  if [ "$KIND" = positive ]; then
      echo "export COPS_POS_ENABLE=1 COPS_POS_NIGHT_RATIO=$RATIO COPS_POS_ACTIVITY=exp COPS_POS_FLOOR_RATIO=$RATIO"
      if [ "$RATIO" = 1e-5 ]; then echo "export COPS_POS_FLOOR_CO2=1"; else echo "unset COPS_POS_FLOOR_CO2"; fi
  else echo "for v in \$(env | grep -o -E '^(COPS_POS_[A-Z_]*|COPS_SLOPEMUL_[A-Z_]*|COPS_NONNEG_[A-Z_]*)'); do unset \$v; done"; fi
  echo "export MODEL_EXT=$X"; echo "cd $APP"; printf '%q ' "${CMD[@]}"; echo; } > "$P/commands.sh"
echo "$(date -u +%FT%TZ) job $SLURM_JOB_ID on $(hostname): $NG GPUs ($DEV); $KIND chain $CHAIN_ROOT (${SLOPE_MODEL:-lambda ${LAMBDA:-} ratio ${RATIO:-}}), $CYCLES cycles; code $P/code (MODEL_EXT $X, chain code $CHAIN_CODE); log $P/chain.log"
bash -c 'echo 1000 > /proc/self/oom_score_adj; exec "$@"' _ "${CMD[@]}" > "$P/chain.log" 2>&1; rc=$?
echo "$(date -u +%FT%TZ) chain driver: exit $rc"
grep -E '^(chain_onejob: |positive_fit: (start check|Jacobian)|lambda |gpu_model: bound|gpu_fit: [0-9.]+ s wall|positive_finetune: (bound|final|check 1|validation|complete|stopped|ERROR)|(slope_finetune|nonneg_finetune|nonneg): (final|validation|start)|finetune_state.npz:|stopped before|.* complete: )' "$P/chain.log" | cut -c1-300
free -g > "$P/free_at_end.txt"; df -h "$CHAIN_ROOT" > "$P/df_at_end.txt" 2>&1
sacct -j "$SLURM_JOB_ID" --format=JobID,State,Elapsed,MaxRSS,ExitCode > "$P/sacct.txt" 2>&1 || true
[ $rc -eq 75 ] && echo "stopped before the wall-time limit: resubmit the same command to continue the chain"
exit $rc
