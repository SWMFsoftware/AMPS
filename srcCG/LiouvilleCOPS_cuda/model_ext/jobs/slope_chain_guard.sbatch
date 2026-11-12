#!/bin/bash
# ==== version-controlled copy: header added 2026-11-11, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/jobs/slope_chain_guard.sbatch (on Voyager: /data001/highresolu_vt/cops_output/jobs/slope_chain_guard.sbatch)
# SHA256 of the run copy = SHA256 of this file without these header lines: 52dd2960d0377f9ca64a39b1d3c4fe2f502684ff770cc0d508eb3ceafae73e9a
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: job (CPU, seconds). Guard job between two stages of a slope / bump chain (jobs/slope_cycles.sh; 2026-11-10,
#   after the dependency flaw 124014 -> 124015): depends on afterok of any attempt of the stage and exits 0 only if the
#   stage's completion marker is present (KIND=fit | ft | gate); layout tested with dummy jobs 124035-124058 (24 of 24
#   as expected).
# RUN: submitted by jobs/slope_cycles.sh (3a chain guards 124068, 124072, 124077, 124081, 124086).
# ==== end of header ====
# Guard job between two stages of a chain (jobs/slope_cycles.sh; 2026-11-10, after the dependency flaw 124014 -> 124015:
# a next-stage attempt a > 1 used only afternotok of the earlier attempts, and an attempt cancelled by
# --kill-on-invalid-dep counts as "not ok"). Submitted with --dependency=afterok:A1?...?afterok:An (any attempt of
# the stage succeeded) and --kill-on-invalid-dep=yes; the next stage's attempts depend on afterok of THIS job only
# (attempt 1: afterok:G; attempt a > 1: afterok:G,afternotok:B1:...:B(a-1)), so when every attempt of the stage failed
# this job is cancelled and so is every later job (layout tested with dummy jobs 124035-124058,
# checks/chain_deps_20261110T214653Z: 24 of 24 as expected).
# It exits 0 only if the stage's completion marker is present:
#   KIND=fit   TARGET=<fit dir>: TARGET/final -> seg<NN> and TARGET/final/fit_restart.npz with status "complete"
#   KIND=ft    TARGET=<ft dir>:  TARGET/final -> seg<NN> and TARGET/finetune_state.npz with the product check of
#              slope_finetune.sbatch's final_link (complete, stop tolerance|max_nfev, held_next = held_prev + delta,
#              finite grain >= 0)
#   KIND=gate  no marker (the guard of an --after list: any of the listed jobs succeeded)
# CPU only, seconds; writes only its log.
#SBATCH --job-name=slope_chain_guard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:05:00
#SBATCH --partition=normal128cores
#SBATCH --output=/data001/highresolu_vt/cops_output/jobs/logs/slope_chain_guard-%j.out
set -uo pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
PY=/data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda/.venv/bin/python
KIND=${KIND:?set KIND=fit|ft|gate}; TARGET=${TARGET:-}
echo "$(date -u +%FT%TZ) guard job ${SLURM_JOB_ID:-local} on $(hostname): KIND=$KIND TARGET=$TARGET"
case $KIND in
    gate) echo "gate: an upstream job succeeded"; exit 0;;
    fit|ft) [ -n "$TARGET" ] || { echo "set TARGET" >&2; exit 2; }
            [ -L "$TARGET/final" ] || { echo "no completion link $TARGET/final: stage not complete" >&2; exit 1; };;
    *) echo "KIND must be fit, ft or gate" >&2; exit 2;;
esac
$PY - "$KIND" "$TARGET" <<'PY'
import json, os, sys
import numpy as np
kind, d = sys.argv[1], sys.argv[2]
if kind == "fit":
    p = os.path.join(d, "final", "fit_restart.npz")
    if not os.path.isfile(p):
        print("missing %s" % p); sys.exit(1)
    with np.load(p, allow_pickle=False) as z:
        meta = json.loads(str(z["meta_json"]))
    ok = meta.get("status") == "complete"
    print("fit %s -> %s: status %s: %s" % (os.path.join(d, "final"), os.readlink(os.path.join(d, "final")),
                                          meta.get("status"), "OK" if ok else "NOT COMPLETE"))
else:
    p = os.path.join(d, "finetune_state.npz")
    if not os.path.isfile(p):
        print("missing %s" % p); sys.exit(1)
    with np.load(p, allow_pickle=False) as z:
        meta = json.loads(str(z["meta_json"]))
        h, dl, hp = z["held_next"], z["delta"], z["held_prev"]
        A, B = float(z["grain_A_uPa"]), float(z["grain_B_uPa_per_100s"])
    ok = (meta.get("complete") is True and meta.get("stop") in ("tolerance", "max_nfev") and h.shape == dl.shape == hp.shape
          and bool(np.all(np.isfinite(h))) and bool(np.array_equal(h, hp + dl)) and np.isfinite(A) and np.isfinite(B)
          and A >= 0 and B >= 0)
    print("fine tuning %s -> %s: stop %s, complete %s, %d evaluations, grain A %.6g B %.6g: %s" % (
        os.path.join(d, "final"), os.readlink(os.path.join(d, "final")), meta.get("stop"), meta.get("complete"),
        meta.get("nfev_total", -1), A, B, "OK" if ok else "NOT COMPLETE"))
sys.exit(0 if ok else 1)
PY
rc=$?
echo "$(date -u +%FT%TZ) guard exit $rc"
exit $rc
