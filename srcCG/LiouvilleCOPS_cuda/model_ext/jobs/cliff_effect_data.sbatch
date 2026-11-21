#!/bin/bash
# ==== version-controlled copy: header added 2026-11-17, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/jobs/cliff_effect_data.sbatch (on Voyager: /data001/highresolu_vt/cops_output/jobs/cliff_effect_data.sbatch)
# SHA256 of the run copy = SHA256 of this file without these header lines: 4e424d896d7901293d761662266706e5d0a3b0f0c393a8ab55beb5f0869b8c5b
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: Slurm job: data of the slope-at-equal-illumination analysis of one state (KEY, OUTD) (2026-11-17)
# RUN: KEY=<set key> OUTD=<folder> sbatch --export=ALL cliff_effect_data.sbatch
# ==== end of header ====
# Data of the slope-at-equal-illumination analysis of ONE state (figures_cluster/scripts/cliff_effect_data.py; 2026-11-17,
# user: "find whether at the same illuminations condition the source rate depends on the slope angle"). CPU.
# Environment: KEY (set key of the full-sets spec), OUTD (output folder, must not exist), SPEC (default the full-sets spec).
#SBATCH --job-name=cops_cliff_data
#SBATCH --partition=normal128cores
#SBATCH --account=general
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=320G
#SBATCH --time=01:00:00
#SBATCH --output=/data001/highresolu_vt/cops_output/jobs/logs/cops_cliff_data-%j.out
set -uo pipefail
export AMPS=/data001/highresolu_vt/AMPS APP=/data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda
export COPS_OUTPUT_ROOT=/data001/highresolu_vt/cops_output
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1
export NPY_DISABLE_CPU_FEATURES="AVX512F AVX512CD AVX512_KNL AVX512_KNM AVX512_SKX AVX512_CLX AVX512_CNL AVX512_ICL"
S=$COPS_OUTPUT_ROOT/figures_cluster/scripts; PY="$APP/.venv/bin/python -W ignore"
KEY=${KEY:?}; OUTD=${OUTD:?}; SPEC=${SPEC:-$S/data/full_sets_stage2_20261117.json}
[ -e "$OUTD" ] && { echo "$OUTD exists" >&2; exit 2; }
P=$OUTD.provenance-$SLURM_JOB_ID; mkdir -p "$P" || exit 1
cd "$APP" || exit 1
for v in KEY OUTD SPEC; do printf '%s=%s\n' "$v" "${!v}"; done > "$P/inputs.txt"
cp -p "$0" "$P/cliff_effect_data.sbatch.used"
mkdir "$P/code" && cp -p "$S/cliff_effect_data.py" "$SPEC" "$P/code/" && (cd "$P/code" && sha256sum * > SHA256SUMS)
{ git -C "$AMPS" rev-parse HEAD; git -C "$AMPS" status --short; } > "$P/git.txt" 2>&1
{ type module >/dev/null 2>&1 && module list; } > "$P/module_list.txt" 2>&1 || echo "module command not available" >> "$P/module_list.txt"
echo "CPU job: no GPU" > "$P/nvidia_smi_q.txt"; { gcc --version | head -1; lscpu | grep -E 'Model name|^CPU\(s\)'; } > "$P/toolchain.txt" 2>&1
~/.local/bin/uv pip list --python "$APP/.venv/bin/python" > "$P/pip_list.txt" 2>&1
df -h "$COPS_OUTPUT_ROOT" > "$P/df_start.txt"
cmd=($PY "$S/cliff_effect_data.py" --spec "$SPEC" --key "$KEY" --out "$OUTD" --procs 32)
{ echo '#!/bin/bash'; grep '^export ' "$0"; echo "cd $APP"; printf '%q ' "${cmd[@]}"; echo; } > "$P/commands.sh"
bash -c 'echo 1000 > /proc/self/oom_score_adj; exec "$@"' heavy "${cmd[@]}" > "$P/cliff_effect_data.log" 2>&1; rc=$?
tail -n 15 "$P/cliff_effect_data.log"
echo "$(date -u +%FT%TZ) done: rc $rc"
cp -p "$COPS_OUTPUT_ROOT/jobs/logs/cops_cliff_data-$SLURM_JOB_ID.out" "$P/job.log"
exit $rc
