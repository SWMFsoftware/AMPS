#!/bin/bash
# ==== version-controlled copy: header added 2026-11-11, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/jobs/slope_mult_cpu.sbatch (on Voyager: /data001/highresolu_vt/cops_output/jobs/slope_mult_cpu.sbatch)
# SHA256 of the run copy = SHA256 of this file without these header lines: e1917cd51e832a5105d1ef9944d05c182abbe3c2b20eea5d972d731e170301d3
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: job (CPU). CPU steps of the slope-multiplier model (user decisions 2026-11-10): STEP=bins (slope_bins.py, the
#   12M slope-bin bundle) | STEP=tests (test_slope_mult.py), from MODEL_EXT.
# RUN: sbatch --export=ALL,STEP=tests,MODEL_EXT=<snapshot> slope_mult_cpu.sbatch   (bins 123624; tests 123625, 123628,
#   123636, 123892, 123978, 123983, 124012, 124062)
# ==== end of header ====
# CPU steps of the slope-multiplier model (user decisions 2026-11-10; model_ext/slope_bins.py, test_slope_mult.py):
#   STEP=bins   the 12M slope-bin bundle runs/shap8_012m_lategrid_dz1.5/slope_mult/bundle/slope_bins.npz (7 indicators
#               of the effective gravitational slope of the admitted gravity cache named by the earlier base's bundle;
#               checks: cache sha256, mesh identity, (slope >= 60) == cliff_ge_60deg, the lead's bin counts)
#   STEP=tests  the CPU tests (test_slope_mult.py: bins, floor, CPU model Jacobian vs FD, absent option, arguments, end
#               to end with the real driver chain on a synthetic base) in checks/slopemul_cpu_<job>
# MODEL_EXT: the model_ext directory the job copies and runs (a frozen snapshot model_ext_frozen/slopemul_<UTC>).
# Provenance in <product dir>/provenance-<job>/ (commands.sh, git, pip list, the code copy).
#SBATCH --job-name=cops_slopemul_cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=01:00:00
#SBATCH --partition=normal128cores
#SBATCH --output=/data001/highresolu_vt/cops_output/jobs/logs/cops_slopemul_cpu-%j.out
set -uo pipefail
export AMPS=/data001/highresolu_vt/AMPS APP=/data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda
export T=/data001/highresolu_vt/cops_cluster_transfer_20261026 COPS_OUTPUT_ROOT=/data001/highresolu_vt/cops_output
export COPS_DATA_ROOT=$T/rosetta_data_end_of_mission_cops
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1
UNSET=$(env | grep -o -E '^(COPS_POS_[A-Z_]*|COPS_SLOPEMUL_[A-Z_]*)' | sort | paste -sd' ' -)
[ -n "$UNSET" ] && unset $UNSET
STEP=${STEP:-}; export MODEL_EXT=$(readlink -f "${MODEL_EXT:-$COPS_OUTPUT_ROOT/model_ext}")
PY="$APP/.venv/bin/python -W ignore"
R=$COPS_OUTPUT_ROOT/runs/shap8_012m_lategrid_dz1.5; E=$COPS_OUTPUT_ROOT/runs/shap8_012m_flat_earlier
case $STEP in
    bins) D=$R/slope_mult/bundle; OUT=$D/slope_bins.npz
          [ -e "$OUT" ] && { echo "exists: $OUT" >&2; exit 2; }
          CMD=($PY "$MODEL_EXT/slope_bins.py" --run "$R/base" --bundle "$E/bundle/restart_templates.npz"
               --gravity-cache "$E/cache/admitted_gravity_cache.npz" --out "$OUT" --expect-counts lead);;
    tests) D=$COPS_OUTPUT_ROOT/checks/slopemul_cpu_$SLURM_JOB_ID
           CMD=($PY "$MODEL_EXT/test_slope_mult.py" --scratch "$D/scratch");;
    *) echo "set STEP=bins or STEP=tests" >&2; exit 2;;
esac
mkdir -p "$D"; P=$D/provenance-$SLURM_JOB_ID; mkdir -p "$P"; cd "$APP" || exit 1
cp -p "$0" "$P/slope_mult_cpu.sbatch"
cp -p "$MODEL_EXT"/*.py "$P/" 2>/dev/null; [ -f "$MODEL_EXT/SHA256SUMS" ] && cp -p "$MODEL_EXT/SHA256SUMS" "$P/"
{ git -C "$AMPS" rev-parse HEAD; git -C "$AMPS" status --short; } > "$P/git.txt" 2>&1
~/.local/bin/uv pip list --python "$APP/.venv/bin/python" > "$P/pip_list.txt" 2>&1
{ echo '#!/bin/bash'; echo "# job $SLURM_JOB_ID on $(hostname), $(date -u +%FT%TZ); AMPS $(git -C "$AMPS" rev-parse HEAD)"
  grep '^export ' "$0"; echo "cd $APP"; printf '%q ' "${CMD[@]}"; echo; } > "$P/commands.sh"
echo "$(date -u +%FT%TZ) job $SLURM_JOB_ID on $(hostname): STEP=$STEP MODEL_EXT=$MODEL_EXT"
"${CMD[@]}" > "$P/run.log" 2>&1; rc=$?
grep -v '^kernels' "$P/run.log" | cut -c1-400 | tail -60
echo "$(date -u +%FT%TZ) exit $rc"
exit $rc
