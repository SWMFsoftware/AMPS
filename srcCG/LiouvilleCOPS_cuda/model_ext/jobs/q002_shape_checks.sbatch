#!/bin/bash
# ==== version-controlled copy: header added 2026-11-18, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/tests/q002_base_20261118/jobs/q002_shape_checks.sbatch (on Voyager: /data001/highresolu_vt/cops_output/tests/q002_base_20261118/jobs/q002_shape_checks.sbatch)
# SHA256 of the run copy = SHA256 of this file without these header lines: 9aecdafc9de461571380abe2ae8a2eb35a4b781707cb742213f9fbe3ad783ffa
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: job (2026-11-18). Refusal checks of nonneg_maps.py validate-geometry --shape-manifest on base_q002 (no option;
#   the old base's manifest; --skip-bdf-hash) on symlinks of a validated geometry without geometry_validation.*, and
#   test_nonneg_contrib.py of nonneg_contrib_20261118T192126Z (v1 job 136693 stopped at the existing validation file;
#   v2 job 136699).
# ==== end of header ====
#SBATCH --job-name=q002_shape_checks
#SBATCH --partition=normal128cores
#SBATCH --account=general
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=96G
#SBATCH --time=01:00:00
#SBATCH --output=/data001/highresolu_vt/cops_output/tests/q002_base_20261118/logs/%x-%j.out
# 2026-11-18 (lead 13:30 CST): checks of nonneg_maps.py validate-geometry --shape-manifest (snapshot
# nonneg_maps_20261118T192126Z) on base_q002 with the geometry of the maps job GEOM_JOB (validated there WITH
# --shape-manifest base_q002/cpp_run_manifest.json --shape-manifest-must-equal base/run_manifest.json):
#   N1 without --shape-manifest: must fail (base_q002 has no run_manifest.json; no silent fallback)
#   N2 --shape-manifest base/run_manifest.json (the OLD base's: records another COPS_model.dat): must be refused
#   N3 --shape-manifest ... --skip-bdf-hash: must be refused
#   T  test_nonneg_contrib.py of the new contribution snapshot nonneg_contrib_20261118T192126Z
# v2: on symlinks of the geometry without geometry_validation.* (see below). Each negative check's process is expected to EXIT NON-ZERO with the named message; the job exits 0 only if all hold.
set -u
export AMPS=/data001/highresolu_vt/AMPS APP=/data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda
export COPS_OUTPUT_ROOT=/data001/highresolu_vt/cops_output
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1
R=$COPS_OUTPUT_ROOT/runs/shap8_012m_lategrid_dz1.5; W=$COPS_OUTPUT_ROOT/tests/q002_base_20261118
MX=$COPS_OUTPUT_ROOT/model_ext_frozen/nonneg_maps_20261118T192126Z
CX=$COPS_OUTPUT_ROOT/model_ext_frozen/nonneg_contrib_20261118T192126Z
GEOM=$R/maps_geometry/job_${GEOM_JOB:?GEOM_JOB required}/geometry
PY=$APP/.venv/bin/python
P=$W/provenance/job-$SLURM_JOB_ID; mkdir -p "$P" || exit 1
cp -p "$0" "$P/q002_shape_checks.sbatch.used"
{ git -C "$AMPS" rev-parse HEAD; git -C "$AMPS" status --short; } > "$P/git.txt" 2>&1
cd "$APP" || exit 1
# v2 (2026-11-18): validate-geometry never overwrites geometry_validation.json, so the negative checks run on a folder of
# symlinks to the validated geometry's files WITHOUT geometry_validation.* (job 136693 stopped at "... exists" in v1)
NG=$W/geom_negative_$SLURM_JOB_ID; mkdir "$NG" || exit 1
for f in "$GEOM"/*; do case $(basename "$f") in geometry_validation.*) ;; *) ln -s "$f" "$NG/";; esac; done
GEOM=$NG
bad=0
neg() {  # $1 name, $2 expected message fragment, rest: the arguments
  local name=$1 want=$2; shift 2
  "$PY" -W ignore "$MX/nonneg_maps.py" validate-geometry --geometry "$GEOM" --base "$R/base_q002" "$@" > "$P/$name.log" 2>&1
  local rc=$?
  if [ $rc -ne 0 ] && grep -q -- "$want" "$P/$name.log"; then echo "PASS $name: rc $rc, '$want'"; else
    echo "FAIL $name: rc $rc"; tail -3 "$P/$name.log"; bad=1; fi
}
neg N1_no_option "run_manifest.json"
neg N2_old_base_manifest "not this base's manifest" --shape-manifest "$R/base/run_manifest.json"
neg N3_skip_hash "requires the BDF hash check" --shape-manifest "$R/base_q002/cpp_run_manifest.json" --skip-bdf-hash
"$PY" -W ignore "$CX/test_nonneg_contrib.py" > "$P/test_nonneg_contrib.log" 2>&1 && echo "PASS T test_nonneg_contrib: $(tail -1 "$P/test_nonneg_contrib.log")" || {
  echo "FAIL T test_nonneg_contrib"; tail -3 "$P/test_nonneg_contrib.log"; bad=1; }
echo "$(date -u +%FT%TZ) checks done, bad=$bad"
exit $bad
