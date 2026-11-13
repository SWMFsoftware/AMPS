#!/bin/bash
# ==== version-controlled copy: header added 2026-11-13, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/jobs/source_area_compare.sbatch (on Voyager: /data001/highresolu_vt/cops_output/jobs/source_area_compare.sbatch)
# SHA256 of the run copy = SHA256 of this file without these header lines: 42ee2d858637eb3c9b2e35aedcb60bd384def63ed1740d8e8ba473db4516b622
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: job script (CPU): source_area_compare.py with provenance.
# RUN: A=... A_LABEL=... B=... B_LABEL=... OUTP=... sbatch --export=ALL source_area_compare.sbatch
# ==== end of header ====
# Source-area locations and slope-bin pressure shares of two fitted states (figures_cluster/scripts/source_area_compare.py;
# 2026-11-13, report items for the softplus refit with slope multipliers >= 1). CPU job; read only; writes OUTP.txt/.json.
# Environment (submitting shell, --export=ALL): A, A_LABEL, B, B_LABEL (contribution folders and labels), OUTP (output
# prefix), RUN (fit base $RUN/base), MX (frozen nonneg_maps code). Provenance in $(dirname OUTP)/provenance-<job>/.
#SBATCH --job-name=cops_sa_compare
#SBATCH --partition=normal128cores
#SBATCH --account=general
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=00:30:00
#SBATCH --output=/data001/highresolu_vt/cops_output/jobs/logs/cops_sa_compare-%j.out
set -uo pipefail
export AMPS=/data001/highresolu_vt/AMPS APP=/data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda
export COPS_OUTPUT_ROOT=/data001/highresolu_vt/cops_output
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1 MPLBACKEND=Agg
R=${RUN:-$COPS_OUTPUT_ROOT/runs/shap8_012m_lategrid_dz1.5}
MX=${MX:-$COPS_OUTPUT_ROOT/model_ext_frozen/nonneg_maps_20261112T021623Z}
A=${A:?}; B=${B:?}; A_LABEL=${A_LABEL:?}; B_LABEL=${B_LABEL:?}; OUTP=${OUTP:?}
S=$COPS_OUTPUT_ROOT/figures_cluster/scripts
PY="$APP/.venv/bin/python -W ignore"
stamp() { echo "$(date -u +%FT%TZ) $*"; }
mkdir -p "$(dirname "$OUTP")" || exit 1
P=$(dirname "$OUTP")/provenance-$SLURM_JOB_ID; mkdir "$P" || exit 1
cd "$APP" || exit 1
stamp "job $SLURM_JOB_ID on $(hostname)"
for v in R MX A A_LABEL B B_LABEL OUTP; do printf '%s=%s\n' "$v" "${!v}"; done > "$P/inputs.txt"
cp -p "$0" "$P/source_area_compare.sbatch.used"
mkdir "$P/code" && cp -p "$S/source_area_compare.py" "$APP/cops_cuda/figures.py" "$MX/nonneg_maps.py" "$P/code/" && (cd "$P/code" && sha256sum * > SHA256SUMS)
{ git -C "$AMPS" rev-parse HEAD; git -C "$AMPS" status --short; } > "$P/git.txt" 2>&1
{ type module >/dev/null 2>&1 && module list; } > "$P/module_list.txt" 2>&1 || echo "module command not available (no modules loaded)" >> "$P/module_list.txt"
echo "CPU job: no GPU" > "$P/nvidia_smi_q.txt"
{ gcc --version | head -1; lscpu | grep -E 'Model name|^CPU\(s\)'; } > "$P/toolchain.txt" 2>&1
~/.local/bin/uv pip list --python "$APP/.venv/bin/python" > "$P/pip_list.txt" 2>&1
{ echo '#!/bin/bash'; grep '^export ' "$0"; echo "cd $APP"; } > "$P/commands.sh"
cmd=($PY "$S/source_area_compare.py" --a "$A" --a-label "$A_LABEL" --b "$B" --b-label "$B_LABEL" --base "$R/base" --maps-code "$MX/nonneg_maps.py" --out "$OUTP")
printf '%q ' "${cmd[@]}" >> "$P/commands.sh"; echo >> "$P/commands.sh"
bash -c 'echo 1000 > /proc/self/oom_score_adj; exec "$@"' heavy "${cmd[@]}" > "$P/compare.log" 2>&1; rc=$?
tail -n 40 "$P/compare.log"
stamp "done: rc $rc"
cp -p "$COPS_OUTPUT_ROOT/jobs/logs/cops_sa_compare-$SLURM_JOB_ID.out" "$P/job.log"
exit $rc
