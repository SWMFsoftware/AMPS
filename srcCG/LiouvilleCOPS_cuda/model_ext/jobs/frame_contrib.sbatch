#!/bin/bash
# ==== version-controlled copy: header added 2026-11-13, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/jobs/frame_contrib.sbatch (on Voyager: /data001/highresolu_vt/cops_output/jobs/frame_contrib.sbatch)
# SHA256 of the run copy = SHA256 of this file without these header lines: 7da76aac81bd17ce85d7ba69eaf0a75b01d5e46245ead3032b74e3a50bbaaf21
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: Job of frame_contrib.py with provenance.
# RUN: sbatch --export=ALL (CHAIN, FB, FB_TAIL, OUTC)
# ==== end of header ====
# Per-facet contribution per unit area at every movie frame (figures_cluster/scripts/frame_contrib.py; 2026-11-14), the
# c9 state on the frames base; verification and the fixed colour range. CPU job.
# Environment: CHAIN, NCYC (default 9), FB (frames base), FB_TAIL (optional second build), OUTC (output folder, must not exist), CX (frozen nonneg_contrib).
#SBATCH --job-name=cops_frame_contrib
#SBATCH --partition=normal128cores
#SBATCH --account=general
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=320G
#SBATCH --time=02:00:00
#SBATCH --output=/data001/highresolu_vt/cops_output/jobs/logs/cops_frame_contrib-%j.out
set -uo pipefail
export AMPS=/data001/highresolu_vt/AMPS APP=/data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda
export COPS_OUTPUT_ROOT=/data001/highresolu_vt/cops_output
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1
export NPY_DISABLE_CPU_FEATURES="AVX512F AVX512CD AVX512_KNL AVX512_KNM AVX512_SKX AVX512_CLX AVX512_CNL AVX512_ICL"
R=$COPS_OUTPUT_ROOT/runs/shap8_012m_lategrid_dz1.5; CH=${CHAIN:?}; NCYC=${NCYC:-9}; FB=${FB:?}; OUTC=${OUTC:?}; FB_TAIL=${FB_TAIL:-}
CX=${CX:-$COPS_OUTPUT_ROOT/model_ext_frozen/nonneg_contrib_20261112T210012Z}; PROCS=${PROCS:-48}
S=$COPS_OUTPUT_ROOT/figures_cluster/scripts; PY="$APP/.venv/bin/python -W ignore"
[ -e "$OUTC" ] && { echo "$OUTC exists" >&2; exit 2; }
P=$OUTC.provenance-$SLURM_JOB_ID; mkdir -p "$P" || exit 1
cd "$APP" || exit 1
for v in R CH NCYC FB FB_TAIL OUTC CX PROCS; do printf '%s=%s\n' "$v" "${!v}"; done > "$P/inputs.txt"
cp -p "$0" "$P/frame_contrib.sbatch.used"
mkdir "$P/code" && cp -p "$S/frame_contrib.py" "$CX/nonneg_contrib.py" "$CX/nonneg_maps.py" "$P/code/" && (cd "$P/code" && sha256sum * > SHA256SUMS)
{ git -C "$AMPS" rev-parse HEAD; git -C "$AMPS" status --short; } > "$P/git.txt" 2>&1
echo "CPU job: no GPU" > "$P/nvidia_smi_q.txt"; { gcc --version | head -1; lscpu | grep -E 'Model name|^CPU\(s\)'; } > "$P/toolchain.txt" 2>&1
~/.local/bin/uv pip list --python "$APP/.venv/bin/python" > "$P/pip_list.txt" 2>&1
cmd=($PY "$S/frame_contrib.py" --contrib-code "$CX" --ft-dir "$CH/c$NCYC/ft/final" --fit-dir "$CH/c$NCYC/fit/final" --base "$R/base"
     --frames-base "$FB" ${FB_TAIL:+--frames-base-tail "$FB_TAIL"} --out "$OUTC" --slope-bins "$CH/contributions_c$NCYC/slope_bins.npz" --procs "$PROCS")
{ echo '#!/bin/bash'; grep '^export ' "$0"; echo "cd $APP"; printf '%q ' "${cmd[@]}"; echo; } > "$P/commands.sh"
bash -c 'echo 1000 > /proc/self/oom_score_adj; exec "$@"' heavy "${cmd[@]}" > "$P/frame_contrib.log" 2>&1; rc=$?
tail -n 60 "$P/frame_contrib.log"
echo "$(date -u +%FT%TZ) done: rc $rc"
cp -p "$COPS_OUTPUT_ROOT/jobs/logs/cops_frame_contrib-$SLURM_JOB_ID.out" "$P/job.log"
exit $rc
