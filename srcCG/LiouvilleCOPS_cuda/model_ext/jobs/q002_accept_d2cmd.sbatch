#!/bin/bash
# ==== version-controlled copy: header added 2026-11-18, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/tests/q002_base_20261118/jobs/q002_accept_d2cmd.sbatch (on Voyager: /data001/highresolu_vt/cops_output/tests/q002_base_20261118/jobs/q002_accept_d2cmd.sbatch)
# SHA256 of the run copy = SHA256 of this file without these header lines: 46dd54a5da355ed970e8b0cde596fc7787a6776293e231a52beb572757387fa4
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: job (2026-11-18, Phase 0 step 3, third try). Runs EXACTLY jobs/clip_refit_q002.sbatch (D2, START=warm, base_q002,
#   --allow-legacy-base) with only NFEV=1 and TESTNAME; exits non-zero unless fit rc 0, the legacy warning, the 1551-point
#   consistency <= 1e-6 and the computed kernel set sha256 = the old base manifest's (the D2 refit depends on it).
# ==== end of header ====
#SBATCH --job-name=q002_accept_d2cmd
#SBATCH --partition=gpu
#SBATCH --account=general
#SBATCH --nodelist=gpu002
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gpus-per-node=8
#SBATCH --cpus-per-task=64
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=02:00:00
#SBATCH --output=/data001/highresolu_vt/cops_output/tests/q002_base_20261118/logs/%x-%j.out
# Phase 0 step 3, third try (lead 02:56 CST): job 136539 accepted the legacy base and passed the consistency check
# (1551 points, NG 2.16e-09, RG 4.12e-09) but a PLAIN gpu_fit cannot select the 1551-row lategrid base ("selection
# profile all_available does not match the audited dataset: rows=1551 (expected 1150)"): D2's driver (clip_refit.py +
# positive_problem's selection wrapper for the 401 model-only grid rows) is needed. This job runs EXACTLY the production
# script jobs/clip_refit_q002.sbatch (VARIANT=D2 START=warm, MODEL_EXT = the clip_q002 snapshot, base_q002 with
# --allow-legacy-base) with only NFEV=1 (--max-nfev 1) and TESTNAME (products in clip_refit_tests/<TESTNAME>/, not in the
# production folder), then EXITS NON-ZERO unless: the fit rc is 0; the legacy warning is in fit.log; the consistency
# line covers 1551 points with max rel <= 1e-6; seg01/fit_manifest.json base_cache.kernels.set_sha256 equals the old
# base manifest's outputs.kernels.set_sha256. The production D2 job depends on it with afterok. No check is loosened.
set -u
export COPS_OUTPUT_ROOT=/data001/highresolu_vt/cops_output
R=$COPS_OUTPUT_ROOT/runs/shap8_012m_lategrid_dz1.5
W=$COPS_OUTPUT_ROOT/tests/q002_base_20261118
TN=q002_accept_warm_nfev1_$SLURM_JOB_ID
mkdir -p "$W/provenance/job-$SLURM_JOB_ID" && cp -p "$0" "$W/provenance/job-$SLURM_JOB_ID/q002_accept_d2cmd.sbatch.used"
cp -p "$COPS_OUTPUT_ROOT/jobs/clip_refit_q002.sbatch" "$W/provenance/job-$SLURM_JOB_ID/clip_refit_q002.sbatch.used"
echo "$(date -u +%FT%TZ) running jobs/clip_refit_q002.sbatch VARIANT=D2 START=warm NFEV=1 TESTNAME=$TN"
env VARIANT=D2 START=warm NFEV=1 TESTNAME=$TN \
    MODEL_EXT=$COPS_OUTPUT_ROOT/model_ext_frozen/clip_q002_20261118T081552Z \
    bash "$COPS_OUTPUT_ROOT/jobs/clip_refit_q002.sbatch"; rc=$?
echo "$(date -u +%FT%TZ) clip_refit_q002.sbatch rc=$rc"
D=$R/clip_refit_tests/$TN
/data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda/.venv/bin/python - "$D" "$R/base/run_manifest.json" "$rc" "$SLURM_JOB_ID" <<'PY' | tee "$W/accept/accept_d2cmd_verdict_$SLURM_JOB_ID.txt"
import json, re, sys
d, old_manifest, rc, job = sys.argv[1], sys.argv[2], int(sys.argv[3]), sys.argv[4]
want = json.load(open(old_manifest))["outputs"]["kernels"]["set_sha256"]
try:
    log = open("%s/provenance-%s/fit.log" % (d, job)).read()
except OSError as exc:
    log = ""
    print("fit.log: %s" % exc)
got, legacy_recorded = None, None
try:
    man = json.load(open(d + "/seg01/fit_manifest.json"))
    got = man["base_cache"]["kernels"]["set_sha256"]
    legacy_recorded = man["base_cache"].get("legacy_manifest")
except Exception as exc:
    print("fit_manifest.json: %s" % exc)
m = re.search(r"consistency with LiouvilleCOPS \((\d+) points\): max rel\. diff NG (\S+), RG (\S+)", log)
start = [l for l in log.splitlines() if l.startswith("clip_refit") and "start" in l][:4]
checks = {"fit_rc_0": rc == 0, "legacy_warning_in_fit_log": "Warning: accepting an unmanifested legacy base" in log,
          "consistency_points_1551": bool(m) and m.group(1) == "1551",
          "consistency_le_1e-6": bool(m) and max(float(m.group(2)), float(m.group(3))) <= 1e-6,
          "kernel_set_sha256_equal_old_manifest": got == want, "fit_manifest_records_legacy": legacy_recorded is True}
print("kernel set sha256 computed %s, old manifest %s" % (got, want))
if m:
    print("consistency (%s points): max rel NG %s RG %s" % m.groups())
for l in start:
    print(l[:400])
print(json.dumps(checks))
print("ACCEPT PASSED" if all(checks.values()) else "ACCEPT FAILED")
sys.exit(0 if all(checks.values()) else 1)
PY
v=${PIPESTATUS[0]}
echo "$(date -u +%FT%TZ) verdict rc=$v"
exit $v
