#!/bin/bash
# ==== version-controlled copy: header added 2026-11-11, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/jobs/slope_mult_threads_test.sbatch (on Voyager: /data001/highresolu_vt/cops_output/jobs/slope_mult_threads_test.sbatch)
# SHA256 of the run copy = SHA256 of this file without these header lines: 2529c2fccacdb2f79ce94592c303519bc0d71438ce6322ab509fad7e9f9e0e63
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: job (test). Threaded test fit of the slope-multiplier model (user 2026-11-10 "run the threaded test fit", then
#   "Bitwise chunked version"): L = Y a micro-benchmark (threaded OpenBLAS and the chunked pool, bitwise against 1
#   thread), the reference short fit (123635) at L_THREADS 16 / 32, bitwise and timing comparison.
# RUN: sbatch --export=ALL,MODEL_EXT=<snapshot>,TAG=<tag> slope_mult_threads_test.sbatch   (this version: 123979; 123891
#   ran the earlier BLAS_THREADS version .before_chunkedL_*)
# ==== end of header ====
# Threaded test fit of the slope-multiplier model (user 2026-11-10: "can CPU calculation be parallelized with OpenMP?"
# -> "run the threaded test fit"; then "Bitwise chunked version": L_THREADS = chunked L = Y a, OpenBLAS at 1 thread;
# L_MODE=openblas reproduces the first test's threaded OpenBLAS). One job on one node:
#  1. micro-benchmark: L = Y a (the fitter's 12.4M x 25 matvec, fit_cops.read_elements' arrays) with numpy's OpenBLAS at
#     THREADS_BENCH (default 1 8 16 32 64) threads (test_slope_mult.py --blas-bench: threaded OpenBLAS and the chunked
#     thread pool), bitwise vs 1 thread;
#  2. the one-node reference short fit (job 123635: START=srccg, NFEV 5, CHECK 1, local transport) with
#     L_THREADS = each of THREADS (default 16 32), L_MODE (default chunked): jobs/slope_mult_fit.sbatch run in this
#     allocation (one bind each), products runs/shap8_012m_lategrid_dz1.5/slope_mult_tests/<PREFIX><N>_<TAG>/seg01
#     (PREFIX f_lchunk for chunked, e_blas for openblas);
#  3. comparison with the reference (test_slope_mult.py --compare-runs): every evaluation's cost, x, predictions,
#     COPS_model.dat, start cost, Jacobian check (bitwise), time per evaluation / Jacobian, the L step, the worker
#     processes' BLAS threads, host CPU load (process_load.txt of each run).
# Settings: MODEL_EXT (required, a frozen snapshot), TAG (required, names the directories), THREADS, THREADS_BENCH, REF.
#SBATCH --job-name=cops12m_slopemul_blas
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gpus-per-node=8
#SBATCH --cpus-per-task=64
#SBATCH --mem=0
#SBATCH --exclusive
#SBATCH --time=01:00:00
#SBATCH --partition=gpu
#SBATCH --account=general
#SBATCH --exclude=gpu003
#SBATCH --output=/data001/highresolu_vt/cops_output/jobs/logs/cops12m_slopemul_blas-%j.out
set -uo pipefail
export AMPS=/data001/highresolu_vt/AMPS APP=/data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda
export COPS_OUTPUT_ROOT=/data001/highresolu_vt/cops_output
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1
MODEL_EXT=$(readlink -f "${MODEL_EXT:?set MODEL_EXT to a frozen snapshot}"); TAG=${TAG:?set TAG}
THREADS=${THREADS:-16 32}; THREADS_BENCH=${THREADS_BENCH:-1,8,16,32,64}; L_MODE=${L_MODE:-chunked}
case $L_MODE in chunked) PREFIX=f_lchunk;; openblas) PREFIX=e_blas;; *) echo "L_MODE chunked|openblas" >&2; exit 2;; esac
R=$COPS_OUTPUT_ROOT/runs/shap8_012m_lategrid_dz1.5; J=$COPS_OUTPUT_ROOT/jobs
REF=${REF:-$R/slope_mult_tests/d_srccg_20261110T094013Z/seg01}
PY="$APP/.venv/bin/python -W ignore"
D=$R/slope_mult_tests/blas_$TAG; [ -e "$D" ] && { echo "exists: $D" >&2; exit 2; }
mkdir -p "$D"; P=$D/provenance-$SLURM_JOB_ID; mkdir -p "$P"; cd "$APP" || exit 1
cp -p "$0" "$P/slope_mult_threads_test.sbatch"; cp -p "$J/slope_mult_fit.sbatch" "$P/"
cp -p "$MODEL_EXT"/*.py "$MODEL_EXT/SHA256SUMS" "$P/" 2>/dev/null
{ git -C "$AMPS" rev-parse HEAD; git -C "$AMPS" status --short; } > "$P/git.txt" 2>&1
nvidia-smi -L > "$P/nvidia_smi_L.txt" 2>&1; lscpu > "$P/lscpu.txt" 2>&1
echo "$(date -u +%FT%TZ) job $SLURM_JOB_ID on $(hostname): MODEL_EXT=$MODEL_EXT TAG=$TAG THREADS='$THREADS' REF=$REF"
{ echo '#!/bin/bash'; echo "# job $SLURM_JOB_ID on $(hostname), $(date -u +%FT%TZ)"; grep '^export ' "$0"; echo "cd $APP"; } > "$P/commands.sh"
BCMD=($PY "$MODEL_EXT/test_slope_mult.py" --blas-bench "$R/base" --threads "$THREADS_BENCH" --json-out "$D/blas_bench.json")
printf '%q ' "${BCMD[@]}" >> "$P/commands.sh"; echo >> "$P/commands.sh"
"${BCMD[@]}" > "$P/blas_bench.log" 2>&1; echo "micro-benchmark exit $?"; grep blas_bench "$P/blas_bench.log"
SEGS=()
for N in $THREADS; do
    echo "$(date -u +%FT%TZ) short fit with L_THREADS=$N L_MODE=$L_MODE"
    echo "TESTNAME=${PREFIX}${N}_$TAG L_THREADS=$N L_MODE=$L_MODE START=srccg CHECK=1 NFEV=5 MODEL_EXT=$MODEL_EXT bash $J/slope_mult_fit.sbatch" >> "$P/commands.sh"
    TESTNAME=${PREFIX}${N}_$TAG L_THREADS=$N L_MODE=$L_MODE START=srccg CHECK=1 NFEV=5 MODEL_EXT=$MODEL_EXT bash "$J/slope_mult_fit.sbatch" \
        > "$P/fit_l$N.log" 2>&1
    echo "  exit $?"; grep -E "L = Y a of the fitter|start check \(start|^lambda|gpu_fit: [0-9]|peak used" "$P/fit_l$N.log" | cut -c1-400
    SEGS+=("$R/slope_mult_tests/${PREFIX}${N}_$TAG/seg01")
done
CCMD=($PY "$MODEL_EXT/test_slope_mult.py" --compare-runs "$REF" "${SEGS[@]}" --json-out "$D/compare.json")
printf '%q ' "${CCMD[@]}" >> "$P/commands.sh"; echo >> "$P/commands.sh"
"${CCMD[@]}" > "$P/compare.log" 2>&1; rc=$?
grep -E "PASS|FAIL|SUMMARY" "$P/compare.log" | cut -c1-900
echo "$(date -u +%FT%TZ) comparison exit $rc"
exit $rc
