#!/bin/bash
# ==== version-controlled copy: header added 2026-11-20, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/jobs/full_set_q002_finish.sbatch (on Voyager: /data001/highresolu_vt/cops_output/jobs/full_set_q002_finish.sbatch)
# SHA256 of the run copy = SHA256 of this file without these header lines: a23ad35628a270486a21c9a70b82888e2176a3f2338d83344ca257a359499da2
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: Slurm CPU job: finish_full_set.sh of one new-law full set with provenance (2026-11-20)
# RUN: KEY=<set key> JOB_SET4=<set1/set4 job> sbatch --dependency=afterok:<the set's jobs> --export=ALL
#   full_set_q002_finish.sbatch
# CHANGED: 2026-11-23 full set of the all-held candidate: SPEC (default the 2026-11-20 spec)
# ==== end of header ====
# Last CPU steps of ONE new-law full set (finish_full_set.sh via full_set_q002.py finish: visibility captions, regions/,
# model_data captions, facts, CAPTIONS.md) with provenance in <set>/provenance-<job>-finish/ (2026-11-20). Environment:
# KEY (set key of data/full_sets_q002_20261120.json), JOB_SET4 (its set1/set4/pressures job).
#SBATCH --job-name=cops_fullset_finish
#SBATCH --partition=normal128cores
#SBATCH --account=general
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=00:45:00
#SBATCH --output=/data001/highresolu_vt/cops_output/jobs/logs/cops_fullset_finish-%j.out
set -uo pipefail
export AMPS=/data001/highresolu_vt/AMPS APP=/data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda
export COPS_OUTPUT_ROOT=/data001/highresolu_vt/cops_output
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1 MPLBACKEND=Agg
S=$COPS_OUTPUT_ROOT/figures_cluster/scripts; KEY=${KEY:?}; JOB_SET4=${JOB_SET4:?}
# 2026-11-23 (full set of the all-held candidate): SPEC = the full-set spec file (default the 2026-11-20 one, unchanged)
SPEC=${SPEC:-$S/data/full_sets_q002_20261120.json}
O=$COPS_OUTPUT_ROOT/figures_cluster/12M/full_set_$KEY; P=$O/provenance-$SLURM_JOB_ID-finish; mkdir -p "$P/code" || exit 1
cd "$APP" || exit 1
printf 'KEY=%s\nJOB_SET4=%s\nSPEC=%s\n' "$KEY" "$JOB_SET4" "$SPEC" > "$P/inputs.txt"
cp -p "$0" "$P/full_set_q002_finish.sbatch.used"
cp -p "$S/full_set_q002.py" "$S/regions_figures.py" "$S/model_data_captions.py" "$S/full_set_facts.py" "$S/collect_captions.py" \
      "$SPEC" "$COPS_OUTPUT_ROOT/jobs/finish_full_set.sh" "$P/code/" && (cd "$P/code" && sha256sum * > SHA256SUMS)
{ git -C "$AMPS" rev-parse HEAD; git -C "$AMPS" status --short; } > "$P/git.txt" 2>&1
{ type module >/dev/null 2>&1 && module list; } > "$P/module_list.txt" 2>&1 || echo "module command not available" >> "$P/module_list.txt"
echo "CPU job: no GPU" > "$P/nvidia_smi_q.txt"; { gcc --version | head -1; lscpu | grep -E 'Model name|^CPU\(s\)'; } > "$P/toolchain.txt" 2>&1
~/.local/bin/uv pip list --python "$APP/.venv/bin/python" > "$P/pip_list.txt" 2>&1
cmd=("$APP/.venv/bin/python" -W ignore "$S/full_set_q002.py" finish --spec "$SPEC" --key "$KEY" --job-set4 "$JOB_SET4")
{ echo '#!/bin/bash'; grep '^export ' "$0"; echo "cd $APP"; printf '%q ' "${cmd[@]}"; echo; } > "$P/commands.sh"
bash -c 'echo 1000 > /proc/self/oom_score_adj; exec "$@"' heavy "${cmd[@]}" > "$P/finish.log" 2>&1; rc=$?
tail -n 40 "$P/finish.log"
echo "$(date -u +%FT%TZ) done: rc $rc"
cp -p "$COPS_OUTPUT_ROOT/jobs/logs/cops_fullset_finish-$SLURM_JOB_ID.out" "$P/job.log"
exit $rc
