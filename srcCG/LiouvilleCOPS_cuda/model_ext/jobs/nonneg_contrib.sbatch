#!/bin/bash
# ==== version-controlled copy: header added 2026-11-13, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/jobs/nonneg_contrib.sbatch (on Voyager: /data001/highresolu_vt/cops_output/jobs/nonneg_contrib.sbatch)
# SHA256 of the run copy = SHA256 of this file without these header lines: 0ebab40a10b614d6586003a31a0c328d943811764d4cf7f87825db2ee2a80bbd
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: job script (CPU, normal128cores). Per-facet contributions of a final state to the modelled NG / RG signals
#   (P% source areas, occurrence, slope-bin fractions) with the frozen nonneg_contrib code; the bitwise state check
#   needs NPY_DISABLE_CPU_FEATURES (AVX-512 features off, set in the script); products <chain>/<CONTRIB_NAME> (default
#   contributions; 2026-11-13 option CONTRIB_NAME, e.g. contributions_c9), CYCLE selects c<CYCLE>/ft/final.
# RUN: 124549 (softplus 5 cycles), 124668 (softplus, multipliers >= 1, c5); 124776 (the same chain, c9:
#   CONTRIB_NAME=contributions_c9 CYCLE=9)
# ==== end of header ====
# Per-facet contributions of the fitted model to the NG and RG surface pressures, P% source areas of the 95 post-gap
# observations (P = 50..90) with the per-facet source-area occurrence, and the slope-bin fractions of every fit point,
# for one final state of the non-negative models: model_ext/nonneg_contrib.py. CPU node (normal128cores), no GPU:
# the work is reading the 4.74e9 kernel records (169 GB) and elementwise FP64 arithmetic on the fit's FP32 records.
# User request 2026-11-12 (relayed by the lead): "determine the area that contributed more than 50% to the RG
# measurement ... 60, 70, 80, 90% ... how many times the surface element was in the most important area ... repeat the
# same for NG"; "add a plot that illustrates [the importance of the cliffs]". Data only (no plots).
# Steps: provenance; tests (test_nonneg_contrib.py); nonneg_contrib.py -> <CHAIN>/contributions (written as
#   contributions.incomplete_<UTC>, renamed on success; never overwritten).
# Environment: CHAIN (chain root), CYCLE (default: the last c<N>), GEOM (a validated nonneg_maps geometry dir, labels
#   t11/t12), MODEL_EXT (code dir; production: a frozen snapshot), PROCS (default SLURM_CPUS_PER_TASK), CONTRIB_NAME
#   (output folder name inside the chain, default contributions; 2026-11-13: contributions_c9 for the extended chain).
#SBATCH --job-name=cops_nonneg_contrib
#SBATCH --partition=normal128cores
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=320G
#SBATCH --time=04:00:00
#SBATCH --output=/data001/highresolu_vt/cops_output/jobs/logs/cops_nonneg_contrib-%j.out
set -uo pipefail
export AMPS=/data001/highresolu_vt/AMPS APP=/data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda
export COPS_OUTPUT_ROOT=/data001/highresolu_vt/cops_output
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1
UNSET=$(env | grep -o -E '^(COPS_POS_[A-Z_]*|COPS_SLOPEMUL_[A-Z_]*|COPS_NONNEG_[A-Z_]*)' | sort | paste -sd' ' -)
[ -n "$UNSET" ] && unset $UNSET
# numpy without its AVX-512 kernels: the state check (A_fit reproduced bitwise from nonneg_facets.npz) needs the AVX2
# path of the GPU nodes where the fit ran (see jobs/nonneg_maps.sbatch, maps_geometry/job_124318/diag_cpu_dispatch.txt)
export NPY_DISABLE_CPU_FEATURES="AVX512F AVX512CD AVX512_KNL AVX512_KNM AVX512_SKX AVX512_CLX AVX512_CNL AVX512_ICL"
R=$COPS_OUTPUT_ROOT/runs/shap8_012m_lategrid_dz1.5
BASE=$R/base
export MODEL_EXT=$(readlink -f "${MODEL_EXT:-$COPS_OUTPUT_ROOT/model_ext}")
PROCS=${PROCS:-${SLURM_CPUS_PER_TASK:-8}}
CHAIN=${CHAIN%/}
GEOM=${GEOM%/}
PY=$APP/.venv/bin/python
stamp() { echo "$(date -u +%FT%TZ) $*"; }
heavy() { bash -c 'echo 1000 > /proc/self/oom_score_adj; exec "$@"' heavy "$@"; }
MODEL=$(basename "$(dirname "$CHAIN")")
N=${CYCLE:-$(ls -d "$CHAIN"/c[0-9]* 2>/dev/null | sed 's#.*/c##' | sort -n | tail -1)}
CONTRIB_NAME=${CONTRIB_NAME:-contributions}
case $CONTRIB_NAME in */*|"") echo "CONTRIB_NAME must be a plain folder name"; exit 2;; esac
FT=$CHAIN/c$N/ft/final; FIT=$CHAIN/c$N/fit/final; OUT=$CHAIN/$CONTRIB_NAME
J=$R/contrib_jobs/job_$SLURM_JOB_ID; P=$J/provenance-$SLURM_JOB_ID
mkdir -p "$R/contrib_jobs" || exit 1
mkdir "$J" "$P" || exit 1
cd "$APP" || exit 1
stamp "job $SLURM_JOB_ID on $(hostname): MODEL=$MODEL CHAIN=$CHAIN c$N GEOM=$GEOM MODEL_EXT=$MODEL_EXT PROCS=$PROCS"
[ -e "$OUT" ] && { stamp "$OUT exists: not overwritten"; exit 1; }
last=$(tail -n 1 "$CHAIN/chain_onejob_events.jsonl")
echo "$last" | grep -q "\"stage\": \"c${N}ft\"" && echo "$last" | grep -q '"rc": 0' || { stamp "chain last event is not c${N}ft rc 0"; exit 1; }
grep -q '^complete yes$' "$FT/finetune_result.txt" || { stamp "$FT not complete"; exit 1; }

# provenance (CLAUDE.md rule 8)
cp -p "$0" "$P/nonneg_contrib.sbatch.used"
for f in nonneg_contrib.py test_nonneg_contrib.py nonneg_maps.py; do cp -p "$MODEL_EXT/$f" "$P/" || exit 1; done
[ -f "$MODEL_EXT/SHA256SUMS" ] && cp -p "$MODEL_EXT/SHA256SUMS" "$P/MODEL_EXT_SHA256SUMS"
( cd "$P" && sha256sum nonneg_contrib.py test_nonneg_contrib.py nonneg_maps.py nonneg_contrib.sbatch.used > code_sha256.txt )
{ git -C "$AMPS" rev-parse HEAD; git -C "$AMPS" status --short; } > "$P/git.txt" 2>&1
{ type module >/dev/null 2>&1 && module list; } > "$P/module_list.txt" 2>&1 || echo "module command not available (no modules loaded)" >> "$P/module_list.txt"
{ command -v nvidia-smi >/dev/null 2>&1 && nvidia-smi -q; } > "$P/nvidia_smi_q.txt" 2>&1 || echo "CPU job on $(hostname): no GPU used (nvidia-smi not available)" >> "$P/nvidia_smi_q.txt"
{ gcc --version | head -1; echo "CSPICE: not used"; echo "nvcc: not used"; lscpu | grep -E 'Model name|^CPU\(s\)'; } > "$P/toolchain.txt" 2>&1
~/.local/bin/uv pip list --python "$PY" > "$P/pip_list.txt" 2>&1 || "$PY" -m pip list > "$P/pip_list.txt" 2>&1
env | grep -E '^(SLURM_JOB_ID|SLURM_JOB_NODELIST|SLURM_CPUS_PER_TASK|MODEL_EXT|CHAIN|CYCLE|GEOM|PROCS|OMP_NUM_THREADS|NPY_DISABLE_CPU_FEATURES|CONTRIB_NAME)=' > "$P/env.txt"
echo "CONTRIB_NAME=$CONTRIB_NAME OUT=$OUT" >> "$P/env.txt"
{ echo '#!/bin/bash'; grep '^export ' "$0"; echo "export MODEL_EXT=$MODEL_EXT"; echo "cd $APP"; } > "$P/commands.sh"
run() { printf '%q ' "$@" >> "$P/commands.sh"; echo >> "$P/commands.sh"; "$@"; }
free -g > "$P/free_start.txt" 2>&1; df -h "$COPS_OUTPUT_ROOT" > "$P/df_start.txt" 2>&1

stamp "1. tests"
run "$PY" -W ignore "$MODEL_EXT/test_nonneg_contrib.py" > "$J/tests.log" 2>&1; rc=$?
tail -n 2 "$J/tests.log"
[ $rc -eq 0 ] || { stamp "tests FAILED"; exit 1; }
stamp "2. nonneg_contrib $MODEL c$N -> $OUT"
run heavy "$PY" -W ignore "$MODEL_EXT/nonneg_contrib.py" --ft-dir "$FT" --fit-dir "$FIT" --base "$BASE" --geometry "$GEOM" \
  --out "$OUT" --expect-model "$MODEL" --procs "$PROCS" > "$J/contrib.log" 2>&1; rc=$?
tail -n 30 "$J/contrib.log" | cut -c1-1500
free -g > "$P/free_end.txt" 2>&1; df -h "$COPS_OUTPUT_ROOT" > "$P/df_end.txt" 2>&1
if [ $rc -eq 0 ]; then
  cp -rp "$P" "$OUT/provenance-$SLURM_JOB_ID"; cp -p "$J/contrib.log" "$J/tests.log" "$OUT/provenance-$SLURM_JOB_ID/"
  du -sh "$OUT"; ls -l "$OUT"
fi
stamp "done rc $rc"
[ -d "$OUT/provenance-$SLURM_JOB_ID" ] && cp -p "$COPS_OUTPUT_ROOT/jobs/logs/cops_nonneg_contrib-$SLURM_JOB_ID.out" "$OUT/provenance-$SLURM_JOB_ID/job.log"
cp -p "$COPS_OUTPUT_ROOT/jobs/logs/cops_nonneg_contrib-$SLURM_JOB_ID.out" "$P/job.log"
exit $rc
