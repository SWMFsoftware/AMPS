#!/bin/bash
# ==== version-controlled copy: header added 2026-11-19, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/jobs/stage2_q002_figures.sbatch (on Voyager: /data001/highresolu_vt/cops_output/jobs/stage2_q002_figures.sbatch)
# SHA256 of the run copy = SHA256 of this file without these header lines: 2052997f48cb9e88ead71e2775b22f34e209e8e3674a6fe4b7000a51873d1741
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: Slurm CPU job: the new-law stage-2 comparison set (stage2_q002_figures.py) with provenance (2026-11-19)
# RUN: OUT=<new folder> sbatch --export=ALL stage2_q002_figures.sbatch
# ==== end of header ====
# Comparison set of the new-law (H2O q_min = 0.02 q_max) stage-2 round vs the data and vs the old law
# (figures_cluster/scripts/stage2_q002_figures.py; lead 2026-11-19, user: "prepare a full set of comparison plots for the
# last round of fitting"). CPU. Environment: OUT (new folder), SPEC (default the q002 spec).
#SBATCH --job-name=cops_stage2_q002
#SBATCH --partition=normal128cores
#SBATCH --account=general
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=01:30:00
#SBATCH --output=/data001/highresolu_vt/cops_output/jobs/logs/cops_stage2_q002-%j.out
set -uo pipefail
export AMPS=/data001/highresolu_vt/AMPS APP=/data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda
export COPS_OUTPUT_ROOT=/data001/highresolu_vt/cops_output
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1 MPLBACKEND=Agg
S=$COPS_OUTPUT_ROOT/figures_cluster/scripts; PY="$APP/.venv/bin/python -W ignore"
OUT=${OUT:?}; SPEC=${SPEC:-$S/data/stage2_q002_spec_20261119.json}
[ -e "$OUT" ] && { echo "$OUT exists" >&2; exit 2; }
P=$OUT.provenance-$SLURM_JOB_ID; mkdir -p "$P/code" || exit 1
cd "$APP" || exit 1
printf 'OUT=%s\nSPEC=%s\n' "$OUT" "$SPEC" > "$P/inputs.txt"
cp -p "$0" "$P/stage2_q002_figures.sbatch.used"
cp -p "$S/stage2_q002_figures.py" "$S/stage2_figures.py" "$S/kappa_table.py" "$S/positive_table.py" "$S/model_obs_figures.py" \
      "$S/options_residual_figure.py" "$S/options_captions.py" "$S/collect_captions.py" "$SPEC" "$P/code/" && (cd "$P/code" && sha256sum * > SHA256SUMS)
{ git -C "$AMPS" rev-parse HEAD; git -C "$AMPS" status --short; } > "$P/git.txt" 2>&1
{ type module >/dev/null 2>&1 && module list; } > "$P/module_list.txt" 2>&1 || echo "module command not available" >> "$P/module_list.txt"
echo "CPU job: no GPU" > "$P/nvidia_smi_q.txt"; { gcc --version | head -1; lscpu | grep -E 'Model name|^CPU\(s\)'; } > "$P/toolchain.txt" 2>&1
~/.local/bin/uv pip list --python "$APP/.venv/bin/python" > "$P/pip_list.txt" 2>&1
cmd=($PY "$S/stage2_q002_figures.py" --spec "$SPEC" --out "$OUT")
{ echo '#!/bin/bash'; grep '^export ' "$0"; echo "cd $APP"; printf '%q ' "${cmd[@]}"; echo; } > "$P/commands.sh"
bash -c 'echo 1000 > /proc/self/oom_score_adj; exec "$@"' heavy "${cmd[@]}" > "$P/stage2_q002_figures.log" 2>&1; rc=$?
tail -n 30 "$P/stage2_q002_figures.log"
echo "$(date -u +%FT%TZ) done: rc $rc"
cp -p "$COPS_OUTPUT_ROOT/jobs/logs/cops_stage2_q002-$SLURM_JOB_ID.out" "$P/job.log"
exit $rc
