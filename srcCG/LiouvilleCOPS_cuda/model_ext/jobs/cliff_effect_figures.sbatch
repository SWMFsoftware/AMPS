#!/bin/bash
# ==== version-controlled copy: header added 2026-11-17, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/jobs/cliff_effect_figures.sbatch (on Voyager: /data001/highresolu_vt/cops_output/jobs/cliff_effect_figures.sbatch)
# SHA256 of the run copy = SHA256 of this file without these header lines: 6b64501bb233bc90c18f375bf262ef9156a866359d505aec452af2069f3c0506
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: Slurm job: figures and answer tables of the slope-at-equal-illumination analysis for all cases of the cases
#   json (2026-11-17)
# RUN: ROOT=<cliff_effect_stage2 folder> sbatch --export=ALL cliff_effect_figures.sbatch
# ==== end of header ====
# Figures and answer tables of the slope-at-equal-illumination analysis, all cases of CASES (cliff_effect_figures.py;
# 2026-11-17). CPU. Environment: ROOT (cliff_effect_stage2 folder with data/<key>/), CASES (default the cases json).
#SBATCH --job-name=cops_cliff_fig
#SBATCH --partition=normal128cores
#SBATCH --account=general
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=00:40:00
#SBATCH --output=/data001/highresolu_vt/cops_output/jobs/logs/cops_cliff_fig-%j.out
set -uo pipefail
export AMPS=/data001/highresolu_vt/AMPS APP=/data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda
export COPS_OUTPUT_ROOT=/data001/highresolu_vt/cops_output
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1
S=$COPS_OUTPUT_ROOT/figures_cluster/scripts; PY="$APP/.venv/bin/python -W ignore"
ROOT=${ROOT:?}; CASES=${CASES:-$S/data/cliff_effect_cases_20261117.json}
P=$ROOT/provenance/figures-$SLURM_JOB_ID; mkdir -p "$P/code" || exit 1
cd "$APP" || exit 1
cp -p "$0" "$P/cliff_effect_figures.sbatch.used"
cp -p "$S/cliff_effect_figures.py" "$CASES" "$P/code/" && (cd "$P/code" && sha256sum * > SHA256SUMS)
{ git -C "$AMPS" rev-parse HEAD; git -C "$AMPS" status --short; } > "$P/git.txt" 2>&1
{ type module >/dev/null 2>&1 && module list; } > "$P/module_list.txt" 2>&1 || echo "module command not available" >> "$P/module_list.txt"
echo "CPU job: no GPU" > "$P/nvidia_smi_q.txt"; { gcc --version | head -1; lscpu | grep -E 'Model name|^CPU\(s\)'; } > "$P/toolchain.txt" 2>&1
~/.local/bin/uv pip list --python "$APP/.venv/bin/python" > "$P/pip_list.txt" 2>&1
rc=0; echo '#!/bin/bash' > "$P/commands.sh"; grep '^export ' "$0" >> "$P/commands.sh"; echo "cd $APP" >> "$P/commands.sh"
n=$($PY -c "import json,sys; print(len(json.load(open(sys.argv[1]))['cases']))" "$CASES")
for ((i=0; i<n; i++)); do
  mapfile -t A < <($PY -c "import json,sys; c=json.load(open(sys.argv[1]))['cases'][int(sys.argv[2])]; print(c['key']); print(c['label']); [print(x) for x in c['caveats']]" "$CASES" $i)
  key=${A[0]}; lab=${A[1]}; cav=()
  for ((j=2; j<${#A[@]}; j++)); do cav+=(--caveat "${A[$j]}"); done
  cmd=($PY "$S/cliff_effect_figures.py" --data "$ROOT/data/$key" --out "$ROOT/$key" --label "$lab" "${cav[@]}")
  printf '%q ' "${cmd[@]}" >> "$P/commands.sh"; echo >> "$P/commands.sh"
  bash -c 'echo 1000 > /proc/self/oom_score_adj; exec "$@"' heavy "${cmd[@]}" > "$P/figures_$key.log" 2>&1 || rc=$?
  echo "$key rc=$rc: $(tail -1 "$P/figures_$key.log")"
done
echo "$(date -u +%FT%TZ) done: rc $rc"
cp -p "$COPS_OUTPUT_ROOT/jobs/logs/cops_cliff_fig-$SLURM_JOB_ID.out" "$P/job.log"
exit $rc
