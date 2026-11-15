#!/bin/bash
# ==== version-controlled copy: header added 2026-11-13, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/jobs/contrib_phase_data.sbatch (on Voyager: /data001/highresolu_vt/cops_output/jobs/contrib_phase_data.sbatch)
# SHA256 of the run copy = SHA256 of this file without these header lines: 32a04d43b2c107f443562152ca0712f65678bef387abcf406a429feae5fcd48f
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: Job of contrib_phase_data.py with provenance.
# RUN: sbatch --export=ALL (CHAIN, OUTD)
# ==== end of header ====
# Data of the contribution-density maps and region shares (figures_cluster/scripts/contrib_phase_data.py; 2026-11-14). CPU.
# Environment: CHAIN, NCYC (default 9), OUTD (output folder, must not exist), LOBE, FRAME (A10 analysis files), CX.
#SBATCH --job-name=cops_contrib_phase
#SBATCH --partition=normal128cores
#SBATCH --account=general
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=320G
#SBATCH --time=02:00:00
#SBATCH --output=/data001/highresolu_vt/cops_output/jobs/logs/cops_contrib_phase-%j.out
set -uo pipefail
export AMPS=/data001/highresolu_vt/AMPS APP=/data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda
export COPS_OUTPUT_ROOT=/data001/highresolu_vt/cops_output
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1
export NPY_DISABLE_CPU_FEATURES="AVX512F AVX512CD AVX512_KNL AVX512_KNM AVX512_SKX AVX512_CLX AVX512_CNL AVX512_ICL"
R=$COPS_OUTPUT_ROOT/runs/shap8_012m_lategrid_dz1.5; CH=${CHAIN:?}; NCYC=${NCYC:-9}; OUTD=${OUTD:?}
A10=$COPS_OUTPUT_ROOT/figures_cluster/12M/movies_softplus_smm1_9cycles_detach300_pan/analysis
LOBE=${LOBE:-$A10/ng_source_regions_c9_lobe.npy}; FRAME=${FRAME:-$A10/ng_source_regions_c9_frame.npz}
CX=${CX:-$COPS_OUTPUT_ROOT/model_ext_frozen/nonneg_contrib_20261112T210012Z}; PROCS=${PROCS:-32}
S=$COPS_OUTPUT_ROOT/figures_cluster/scripts; PY="$APP/.venv/bin/python -W ignore"
[ -e "$OUTD" ] && { echo "$OUTD exists" >&2; exit 2; }
P=$OUTD.provenance-$SLURM_JOB_ID; mkdir -p "$P" || exit 1
cd "$APP" || exit 1
for v in R CH NCYC OUTD LOBE FRAME CX PROCS; do printf '%s=%s\n' "$v" "${!v}"; done > "$P/inputs.txt"
cp -p "$0" "$P/contrib_phase_data.sbatch.used"
mkdir "$P/code" && cp -p "$S/contrib_phase_data.py" "$S/frame_contrib.py" "$CX/nonneg_contrib.py" "$CX/nonneg_maps.py" "$P/code/" && (cd "$P/code" && sha256sum * > SHA256SUMS)
{ git -C "$AMPS" rev-parse HEAD; git -C "$AMPS" status --short; } > "$P/git.txt" 2>&1
echo "CPU job: no GPU" > "$P/nvidia_smi_q.txt"; { gcc --version | head -1; lscpu | grep -E 'Model name|^CPU\(s\)'; } > "$P/toolchain.txt" 2>&1
~/.local/bin/uv pip list --python "$APP/.venv/bin/python" > "$P/pip_list.txt" 2>&1
cmd=($PY "$S/contrib_phase_data.py" --contrib-code "$CX" --ft-dir "$CH/c$NCYC/ft/final" --fit-dir "$CH/c$NCYC/fit/final" --base "$R/base"
     --lobe "$LOBE" --frame "$FRAME" --slope-bins "$CH/contributions_c$NCYC/slope_bins.npz" --out "$OUTD"
     --frame-contrib-py "$S/frame_contrib.py" --procs "$PROCS")
{ echo '#!/bin/bash'; grep '^export ' "$0"; echo "cd $APP"; printf '%q ' "${cmd[@]}"; echo; } > "$P/commands.sh"
bash -c 'echo 1000 > /proc/self/oom_score_adj; exec "$@"' heavy "${cmd[@]}" > "$P/contrib_phase_data.log" 2>&1; rc=$?
tail -n 40 "$P/contrib_phase_data.log"
echo "$(date -u +%FT%TZ) done: rc $rc"
cp -p "$COPS_OUTPUT_ROOT/jobs/logs/cops_contrib_phase-$SLURM_JOB_ID.out" "$P/job.log"
exit $rc
