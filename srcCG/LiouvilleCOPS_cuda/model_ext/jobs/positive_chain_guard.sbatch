#!/bin/bash
# ==== version-controlled copy: header added 2026-11-11, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/jobs/positive_chain_guard.sbatch (on Voyager: /data001/highresolu_vt/cops_output/jobs/positive_chain_guard.sbatch)
# SHA256 of the run copy = SHA256 of this file without these header lines: 6da392229638065a2bb4af63d0fcee295cb2c527e884156201efdcd9fb63d473
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: job (CPU, seconds). Guard job between two stages of a positive-activity chain (user 2026-11-10 "Apply the guard
#   fix"; the layout of jobs/slope_cycles.sh): depends on afterok of any attempt of the stage and exits 0 only if the
#   stage is complete (KIND=fit: final link and status complete or all NFEV evaluations used; KIND=ft: final link and
#   the product check; KIND=gate), so the next stage's attempts start only behind it.
# RUN: submitted by jobs/positive_cycles.sh. Not yet submitted to Slurm (local runs and dry runs only, 2026-11-10; no
#   positive chain since).
# ==== end of header ====
# Guard job between two stages of a positive-activity chain (jobs/positive_cycles.sh; user 2026-11-10 "Apply the guard
# fix"; the same layout as jobs/slope_cycles.sh / slope_chain_guard.sbatch). Submitted with
# --dependency=afterok:A1?...?afterok:An (any attempt of the stage succeeded) and --kill-on-invalid-dep=yes; the next
# stage's attempts depend on afterok of THIS job only (attempt 1: afterok:G; attempt a > 1:
# afterok:G,afternotok:B1:...:B(a-1)), so when every attempt of the stage failed this job is cancelled and so is every
# later job (layout tested with dummy jobs 124035-124058, checks/chain_deps_20261110T214653Z: 24 of 24 as expected).
# It exits 0 only if the stage is complete as the stage's own job script defines it:
#   KIND=fit   TARGET=<fit dir> (positive_fit.sbatch): TARGET/final -> seg<NN> (the script links it only for a
#              complete fit) and TARGET/final/fit_restart.npz with status "complete" or, with NFEV set, all NFEV
#              evaluations used (positive_fit.sbatch's second completion rule: "no evaluations left")
#   KIND=ft    TARGET=<ft dir> (positive_finetune.sbatch): TARGET/final -> seg<NN> and TARGET/finetune_state.npz with
#              final_link's product check (complete, stop tolerance|max_nfev, held_next = held_prev + delta, finite
#              grain >= 0)
#   KIND=gate  no marker (the guard of an --after list: any of the listed jobs succeeded)
# CPU only, seconds; writes only its log.
#SBATCH --job-name=positive_chain_guard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:05:00
#SBATCH --partition=normal128cores
#SBATCH --output=/data001/highresolu_vt/cops_output/jobs/logs/positive_chain_guard-%j.out
set -uo pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
PY=/data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda/.venv/bin/python
KIND=${KIND:?set KIND=fit|ft|gate}; TARGET=${TARGET:-}; NFEV=${NFEV:-}
echo "$(date -u +%FT%TZ) guard job ${SLURM_JOB_ID:-local} on $(hostname): KIND=$KIND TARGET=$TARGET"
case $KIND in
    gate) echo "gate: an upstream job succeeded"; exit 0;;
    fit|ft) [ -n "$TARGET" ] || { echo "set TARGET" >&2; exit 2; }
            [ -L "$TARGET/final" ] || { echo "no completion link $TARGET/final: stage not complete" >&2; exit 1; };;
    *) echo "KIND must be fit, ft or gate" >&2; exit 2;;
esac
$PY - "$KIND" "$TARGET" "$NFEV" <<'PY'
import json, os, sys
import numpy as np
kind, d, nfev = sys.argv[1], sys.argv[2], sys.argv[3]
if kind == "fit":
    p = os.path.join(d, "final", "fit_restart.npz")
    if not os.path.isfile(p):
        print("missing %s" % p); sys.exit(1)
    with np.load(p, allow_pickle=False) as z:
        meta = json.loads(str(z["meta_json"]))
    used = sum(int(s.get("nfev", 0)) for s in meta.get("segments", []))
    ok = meta.get("status") == "complete" or (nfev != "" and used >= int(nfev))
    print("fit %s -> %s: status %s, %d evaluations used (NFEV %s): %s" % (
        os.path.join(d, "final"), os.readlink(os.path.join(d, "final")), meta.get("status"), used, nfev or "-",
        "OK" if ok else "NOT COMPLETE"))
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
