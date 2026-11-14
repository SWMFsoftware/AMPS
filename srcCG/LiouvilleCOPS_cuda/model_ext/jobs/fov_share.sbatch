#!/bin/bash
# ==== version-controlled copy: header added 2026-11-13, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/jobs/fov_share.sbatch (on Voyager: /data001/highresolu_vt/cops_output/jobs/fov_share.sbatch)
# SHA256 of the run copy = SHA256 of this file without these header lines: aaeb787f7c3d76971673860a605ec013be3724281b970b2c7d485dcd7fc3e5a2
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: Job of fov_share.py; 2026-11-14: HEMI_SC, HPREP, HMOVIES, DECK for --hemi-sc (default off).
# RUN: sbatch --export=ALL (see the header)
# ==== end of header ====
# Share of the modelled NG / RG surface signal from directions outside each gauge's IK FOV (figures_cluster/scripts/
# fov_share.py; lead / user 2026-11-13, with the observable-part movies). CPU job (normal128cores), read only.
# Environment (submitting shell, --export=ALL): CHAIN (chain dir), NCYC (default 5), RUN, GEOM (maps geometry of the
# chain), FOV (cops_movies.py fov output, default $CHAIN/movies_work/fov.npz), CX (frozen nonneg_contrib code), OUTP
# (output prefix), PROCS (default SLURM_CPUS_PER_TASK); 2026-11-13: CONE_RG, CONE_NG (half-angles [deg] of circular
# fields of view instead of the IK rectangles; default empty = the IK rectangles). Provenance in $(dirname OUTP)/provenance-<job>/.
#SBATCH --job-name=cops_fov_share
#SBATCH --partition=normal128cores
#SBATCH --account=general
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=320G
#SBATCH --time=02:00:00
#SBATCH --output=/data001/highresolu_vt/cops_output/jobs/logs/cops_fov_share-%j.out
set -uo pipefail
export AMPS=/data001/highresolu_vt/AMPS APP=/data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda
export COPS_OUTPUT_ROOT=/data001/highresolu_vt/cops_output
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1
# as jobs/nonneg_contrib.sbatch: the state checks reproduce the fit's activity bitwise with the SIMD features of the GPU nodes
export NPY_DISABLE_CPU_FEATURES="AVX512F AVX512CD AVX512_KNL AVX512_KNM AVX512_SKX AVX512_CLX AVX512_CNL AVX512_ICL"
R=${RUN:-$COPS_OUTPUT_ROOT/runs/shap8_012m_lategrid_dz1.5}
CH=${CHAIN:?}; NCYC=${NCYC:-5}; GEOM=${GEOM:?}; FOV=${FOV:-$CH/movies_work/fov.npz}
CX=${CX:-$COPS_OUTPUT_ROOT/model_ext_frozen/nonneg_contrib_20261112T210012Z}; OUTP=${OUTP:?}
PROCS=${PROCS:-${SLURM_CPUS_PER_TASK:-16}}; CONE_RG=${CONE_RG:-}; CONE_NG=${CONE_NG:-}
HEMI_SC=${HEMI_SC:-}; HPREP=${HPREP:-}; HMOVIES=${HMOVIES:-}; DECK=${DECK:-rosina-deck}   # 2026-11-14: --hemi-sc (default off)
S=$COPS_OUTPUT_ROOT/figures_cluster/scripts
PY="$APP/.venv/bin/python -W ignore"
stamp() { echo "$(date -u +%FT%TZ) $*"; }
mkdir -p "$(dirname "$OUTP")" || exit 1
P=$(dirname "$OUTP")/provenance-$SLURM_JOB_ID; mkdir "$P" || exit 1
cd "$APP" || exit 1
stamp "job $SLURM_JOB_ID on $(hostname)"
for v in R CH NCYC GEOM FOV CX OUTP PROCS CONE_RG CONE_NG HEMI_SC HPREP HMOVIES DECK; do printf '%s=%s\n' "$v" "${!v}"; done > "$P/inputs.txt"
cp -p "$0" "$P/fov_share.sbatch.used"
mkdir "$P/code" && cp -p "$S/fov_share.py" "$CX/nonneg_contrib.py" "$CX/nonneg_maps.py" "$P/code/" && (cd "$P/code" && sha256sum * > SHA256SUMS)
[ -n "$HMOVIES" ] && cp -p "$HMOVIES" "$P/code/cops_movies.py" && (cd "$P/code" && sha256sum * > SHA256SUMS)
{ git -C "$AMPS" rev-parse HEAD; git -C "$AMPS" status --short; } > "$P/git.txt" 2>&1
{ type module >/dev/null 2>&1 && module list; } > "$P/module_list.txt" 2>&1 || echo "module command not available (no modules loaded)" >> "$P/module_list.txt"
echo "CPU job: no GPU" > "$P/nvidia_smi_q.txt"
{ gcc --version | head -1; lscpu | grep -E 'Model name|^CPU\(s\)'; } > "$P/toolchain.txt" 2>&1
~/.local/bin/uv pip list --python "$APP/.venv/bin/python" > "$P/pip_list.txt" 2>&1
{ echo '#!/bin/bash'; grep '^export ' "$0"; echo "cd $APP"; } > "$P/commands.sh"
cmd=($PY "$S/fov_share.py" --contrib-code "$CX" --ft-dir "$CH/c$NCYC/ft/final" --fit-dir "$CH/c$NCYC/fit/final" \
     --base "$R/base" --geometry "$GEOM" --fov "$FOV" --out "$OUTP" --procs "$PROCS" ${CONE_RG:+--cone-rg "$CONE_RG"} ${CONE_NG:+--cone-ng "$CONE_NG"} \
     ${HEMI_SC:+--hemi-sc "$HEMI_SC" --prep "$HPREP" --movies-py "$HMOVIES" --deck "$DECK"})
printf '%q ' "${cmd[@]}" >> "$P/commands.sh"; echo >> "$P/commands.sh"
free -g > "$P/free_start.txt" 2>&1
bash -c 'echo 1000 > /proc/self/oom_score_adj; exec "$@"' heavy "${cmd[@]}" > "$P/fov_share.log" 2>&1; rc=$?
tail -n 20 "$P/fov_share.log"
free -g > "$P/free_end.txt" 2>&1
stamp "done: rc $rc"
cp -p "$COPS_OUTPUT_ROOT/jobs/logs/cops_fov_share-$SLURM_JOB_ID.out" "$P/job.log"
exit $rc
