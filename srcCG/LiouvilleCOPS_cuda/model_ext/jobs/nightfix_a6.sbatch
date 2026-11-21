#!/bin/bash
# ==== version-controlled copy: header added 2026-11-18, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/tests/nightfix_20261118/jobs/nightfix_a6.sbatch (on Voyager: /data001/highresolu_vt/cops_output/tests/nightfix_20261118/jobs/nightfix_a6.sbatch)
# SHA256 of the run copy = SHA256 of this file without these header lines: 31cd747ab994a1bbd9e664c97e360aa5cc869c6bb92cd75df6b2b8bb661c1336
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: job (2026-11-18). LiouvilleCOPS A6 gate and --integration with TMPDIR on /data001 (compute-node /tmp is noexec).
# RUN: sbatch tests/nightfix_20261118/jobs/nightfix_a6.sbatch (job 136429).
# ==== end of header ====
#SBATCH --job-name=nightfix_a6
#SBATCH --partition=normal128cores
#SBATCH --account=general
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=96G
#SBATCH --time=02:00:00
#SBATCH --output=/data001/highresolu_vt/cops_output/tests/nightfix_20261118/jobs/%x-%j.out
# H2O night-side fix (2026-11-18): the LiouvilleCOPS A6 gate (unittest suite + liouville_cops selftest.in) and its
# --integration variant with the rebuilt liouville_cops (job 136406), rerun with TMPDIR on /data001: on the Voyager
# compute nodes /tmp is a tmpfs mounted noexec, so every test that compiles or copies an executable into a temporary
# directory failed in job 136406 with PermissionError / "not executable" (site artifact, not the law change).
set -u
export AMPS=/data001/highresolu_vt/AMPS APP=/data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda
export T=/data001/highresolu_vt/cops_cluster_transfer_20261026
export COPS_DATA_ROOT=$T/rosetta_data_end_of_mission_cops COPS_OUTPUT_ROOT=/data001/highresolu_vt/cops_output
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1
W=$COPS_OUTPUT_ROOT/tests/nightfix_20261118
P=$W/a6/job-$SLURM_JOB_ID; mkdir -p "$P" "$P/tmp" || exit 1
export TMPDIR=$P/tmp TMP=$P/tmp TEMP=$P/tmp
cp -p "$0" "$P/nightfix_a6.sbatch.used"
{ git -C "$AMPS" rev-parse HEAD; git -C "$AMPS" status --short; } > "$P/git.txt" 2>&1
{ hostname; lscpu | grep 'Model name'; mount | grep ' /tmp '; } > "$P/node.txt" 2>&1
sha256sum "$AMPS/srcCG/LiouvilleCOPS/liouville_cops" > "$P/executable_sha256.txt"
"$APP/.venv/bin/python" -m pip list > "$P/pip_list.txt" 2>&1
cd "$AMPS/srcCG/LiouvilleCOPS" || exit 1
echo "$(date -u +%FT%TZ) A6 gate (TMPDIR=$TMPDIR)"
"$APP/.venv/bin/python" tools/run_a6_validation.py --threads 16 --report "$P/a6_validation.json" > "$P/a6_validation.log" 2>&1; rc1=$?
echo "$(date -u +%FT%TZ) A6 gate --integration"
timeout 5400 "$APP/.venv/bin/python" tools/run_a6_validation.py --integration --threads 16 --report "$P/a6_validation_integration.json" > "$P/a6_validation_integration.log" 2>&1; rc2=$?
echo "$(date -u +%FT%TZ) summary: a6=$rc1 a6_integration=$rc2" | tee "$P/summary.txt"
exit 0
