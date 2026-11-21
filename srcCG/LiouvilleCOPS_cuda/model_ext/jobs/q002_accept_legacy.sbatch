#!/bin/bash
# ==== version-controlled copy: header added 2026-11-18, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/tests/q002_base_20261118/jobs/q002_accept_legacy.sbatch (on Voyager: /data001/highresolu_vt/cops_output/tests/q002_base_20261118/jobs/q002_accept_legacy.sbatch)
# SHA256 of the run copy = SHA256 of this file without these header lines: 204235a851c3edace9b956afa957010e602c14b162e19bd11253f89696dfcfc0
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: job (2026-11-18, Phase 0 step 3, second try). gpu_fit (fit_cops.main) on base_q002 with --allow-legacy-base,
#   --max-nfev 1, --check-tolerance 1e-6; exits non-zero unless rc 0, the legacy warning, the 1551-point consistency
#   line <= 1e-6 and the computed kernel set sha256 = the old base manifest's (the D2 refit depends on it with afterok).
# ==== end of header ====
#SBATCH --job-name=q002_accept_legacy
#SBATCH --partition=gpu
#SBATCH --account=general
#SBATCH --nodelist=gpu002
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gres=gpu:8
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=01:00:00
#SBATCH --output=/data001/highresolu_vt/cops_output/tests/q002_base_20261118/logs/%x-%j.out
# Phase 0 step 3, second try (lead 02:50 CST): job 136534's gpu_fit on base_q002 stopped at fit_cops.main's W01 check
# ("the base has no run_manifest.json ... use --allow-legacy-base"). Lead's decision: no hand-written manifest; the
# official fit_cops flag --allow-legacy-base (strict structural validation of every file, the limitation recorded in the
# fit manifest). This job runs gpu_fit (= the unchanged fit_cops.main) on base_q002 with that flag, --max-nfev 1,
# --check-tolerance 1e-6, and EXITS NON-ZERO unless all of these hold:
#   gpu_fit rc 0; the legacy warning printed; "consistency with LiouvilleCOPS (1551 points)" printed and the run went on
#   (fit_cops exits above 1e-6); the kernel set sha256 that validate_fit_base computed (fit_manifest.json base_cache)
#   equals the old base manifest's outputs.kernels.set_sha256
# so that the D2 refit can depend on it with afterok.
set -u
export AMPS=/data001/highresolu_vt/AMPS APP=/data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda
export T=/data001/highresolu_vt/cops_cluster_transfer_20261026 COPS_OUTPUT_ROOT=/data001/highresolu_vt/cops_output
export COPS_DATA_ROOT=$T/rosetta_data_end_of_mission_cops
export CUDA_PATH=/usr/local/cuda-12.9 CUDA_HOME=/usr/local/cuda-12.9 CUPY_CACHE_DIR=$COPS_OUTPUT_ROOT/.cupy
export PATH=$CUDA_HOME/bin:$PATH LD_LIBRARY_PATH=$CUDA_HOME/lib64${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1
POSUNSET=$(env | grep -o '^COPS_POS_[A-Z_]*' | sort | paste -sd' ' -); [ -n "$POSUNSET" ] && unset $POSUNSET
R=$COPS_OUTPUT_ROOT/runs/shap8_012m_lategrid_dz1.5
W=$COPS_OUTPUT_ROOT/tests/q002_base_20261118
PY=$APP/.venv/bin/python
P=$W/provenance/job-$SLURM_JOB_ID; mkdir -p "$P" "$W/accept" || exit 1
cp -p "$0" "$P/q002_accept_legacy.sbatch.used"
{ git -C "$AMPS" rev-parse HEAD; git -C "$AMPS" status --short; } > "$P/git.txt" 2>&1
nvidia-smi -q > "$P/nvidia_smi_q.txt" 2>&1; "$PY" -m pip list > "$P/pip_list.txt" 2>&1
cd "$APP" || exit 1
O=$W/accept/gpu_fit_q002_legacy
CMD=("$PY" -W ignore tools/gpu_fit.py --devices 0,1,2,3,4,5,6,7 --precision mixed --max-nfev 1 --no-figures \
     -- "$R/base_q002" --profile all_available --degree 4 --enhancement-degree 4 --lambda 0 --jobs 1 --tbounds 30 400 \
     --ftol 1e-10 --xtol 1e-10 --gtol 1e-10 --check-tolerance 1e-6 --allow-legacy-base --out "$O")
printf '%q ' "${CMD[@]}" > "$P/commands.sh"; echo >> "$P/commands.sh"
echo "$(date -u +%FT%TZ) gpu_fit on base_q002 with --allow-legacy-base"
bash -c 'echo 1000 > /proc/self/oom_score_adj; exec "$@"' heavy "${CMD[@]}" > "$W/accept/gpu_fit_q002_legacy.log" 2>&1; rc=$?
echo "$(date -u +%FT%TZ) gpu_fit rc=$rc"
grep -E 'Warning: accepting an unmanifested|activity factors reproduced|consistency with LiouvilleCOPS|Error' "$W/accept/gpu_fit_q002_legacy.log" | cut -c1-300
"$PY" - "$O" "$R/base/run_manifest.json" "$W/accept/gpu_fit_q002_legacy.log" "$rc" <<'PY' | tee "$W/accept/accept_legacy_verdict.txt"
import json, re, sys
out, old_manifest, log, rc = sys.argv[1], sys.argv[2], open(sys.argv[3]).read(), int(sys.argv[4])
want = json.load(open(old_manifest))["outputs"]["kernels"]["set_sha256"]
got = None
try:
    man = json.load(open(out + "/fit_manifest.json"))
    got = man["base_cache"]["kernels"]["set_sha256"]
except Exception as exc:
    print("fit_manifest.json: %s" % exc)
m = re.search(r"consistency with LiouvilleCOPS \((\d+) points\): max rel\. diff NG (\S+), RG (\S+)", log)
checks = {"gpu_fit_rc_0": rc == 0, "legacy_warning_printed": "Warning: accepting an unmanifested legacy base" in log,
          "consistency_line": bool(m), "consistency_points_1551": bool(m) and m.group(1) == "1551",
          "consistency_le_1e-6": bool(m) and max(float(m.group(2)), float(m.group(3))) <= 1e-6,
          "kernel_set_sha256_equal_old_manifest": got == want}
print("kernel set sha256 computed %s, old manifest %s" % (got, want))
if m:
    print("consistency (%s points): max rel NG %s RG %s" % m.groups())
print(json.dumps(checks))
print("ACCEPT PASSED" if all(checks.values()) else "ACCEPT FAILED")
sys.exit(0 if all(checks.values()) else 1)
PY
v=${PIPESTATUS[0]}
echo "$(date -u +%FT%TZ) verdict rc=$v"
exit $v
