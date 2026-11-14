#!/bin/bash
# ==== version-controlled copy: header added 2026-11-13, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/jobs/ng_source_map.sbatch (on Voyager: /data001/highresolu_vt/cops_output/jobs/ng_source_map.sbatch)
# SHA256 of the run copy = SHA256 of this file without these header lines: b21fb7503a1d532779790eb81f4f0ff0db902fb5c7904b42f69f4cbe100c2125
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: Job of ng_source_map.py with provenance.
# RUN: sbatch --export=ALL (CHAIN, OUTP)
# ==== end of header ====
# Per-facet sums of the modelled NG surface signal (late 95 observations, W10x5; all NG observations of the descent):
# figures_cluster/scripts/ng_source_map.py (lead / user 2026-11-13 ~20:05 CST, panorama camera path). CPU job, read only.
# Environment: CHAIN, NCYC (default 9), CONTRIB (contribution data folder, default $CHAIN/contributions_c$NCYC), CX (frozen
# nonneg_contrib code), OUTP (output prefix), PROCS. Provenance in $(dirname OUTP)/provenance-<job>/.
#SBATCH --job-name=cops_ng_source_map
#SBATCH --partition=normal128cores
#SBATCH --account=general
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=200G
#SBATCH --time=01:00:00
#SBATCH --output=/data001/highresolu_vt/cops_output/jobs/logs/cops_ng_source_map-%j.out
set -uo pipefail
export AMPS=/data001/highresolu_vt/AMPS APP=/data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda
export COPS_OUTPUT_ROOT=/data001/highresolu_vt/cops_output
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1
export NPY_DISABLE_CPU_FEATURES="AVX512F AVX512CD AVX512_KNL AVX512_KNM AVX512_SKX AVX512_CLX AVX512_CNL AVX512_ICL"
R=${RUN:-$COPS_OUTPUT_ROOT/runs/shap8_012m_lategrid_dz1.5}
CH=${CHAIN:?}; NCYC=${NCYC:-9}; CONTRIB=${CONTRIB:-$CH/contributions_c$NCYC}
CX=${CX:-$COPS_OUTPUT_ROOT/model_ext_frozen/nonneg_contrib_20261112T210012Z}; OUTP=${OUTP:?}; PROCS=${PROCS:-32}
S=$COPS_OUTPUT_ROOT/figures_cluster/scripts
PY="$APP/.venv/bin/python -W ignore"
stamp() { echo "$(date -u +%FT%TZ) $*"; }
mkdir -p "$(dirname "$OUTP")" || exit 1
P=$(dirname "$OUTP")/provenance-$SLURM_JOB_ID; mkdir "$P" || exit 1
cd "$APP" || exit 1
stamp "job $SLURM_JOB_ID on $(hostname)"
for v in R CH NCYC CONTRIB CX OUTP PROCS; do printf '%s=%s\n' "$v" "${!v}"; done > "$P/inputs.txt"
cp -p "$0" "$P/ng_source_map.sbatch.used"
mkdir "$P/code" && cp -p "$S/ng_source_map.py" "$CX/nonneg_contrib.py" "$CX/nonneg_maps.py" "$P/code/" && (cd "$P/code" && sha256sum * > SHA256SUMS)
{ git -C "$AMPS" rev-parse HEAD; git -C "$AMPS" status --short; } > "$P/git.txt" 2>&1
echo "CPU job: no GPU" > "$P/nvidia_smi_q.txt"
{ gcc --version | head -1; lscpu | grep -E 'Model name|^CPU\(s\)'; } > "$P/toolchain.txt" 2>&1
~/.local/bin/uv pip list --python "$APP/.venv/bin/python" > "$P/pip_list.txt" 2>&1
{ echo '#!/bin/bash'; grep '^export ' "$0"; echo "cd $APP"; } > "$P/commands.sh"
cmd=($PY "$S/ng_source_map.py" --contrib-code "$CX" --ft-dir "$CH/c$NCYC/ft/final" --fit-dir "$CH/c$NCYC/fit/final" \
     --base "$R/base" --late-npz "$CONTRIB/source_areas_NG.npz" --out "$OUTP" --procs "$PROCS")
printf '%q ' "${cmd[@]}" >> "$P/commands.sh"; echo >> "$P/commands.sh"
bash -c 'echo 1000 > /proc/self/oom_score_adj; exec "$@"' heavy "${cmd[@]}" > "$P/ng_source_map.log" 2>&1; rc=$?
tail -n 12 "$P/ng_source_map.log"
stamp "done: rc $rc"
cp -p "$COPS_OUTPUT_ROOT/jobs/logs/cops_ng_source_map-$SLURM_JOB_ID.out" "$P/job.log"
exit $rc
