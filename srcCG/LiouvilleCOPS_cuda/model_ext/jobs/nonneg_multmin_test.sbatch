#!/bin/bash
# ==== version-controlled copy: header added 2026-11-14, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/jobs/nonneg_multmin_test.sbatch (on Voyager: /data001/highresolu_vt/cops_output/jobs/nonneg_multmin_test.sbatch)
# SHA256 of the run copy = SHA256 of this file without these header lines: 4fd2ae6d1fed24dfb631d6f6c111a8111ea4a516492b28cd1505fd6bf2cab6e0
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: job. GPU test of the softplus multiplier lower bound 1 (user 2026-11-12; --slope-mult-min 1).
# RUN: sbatch --nodelist=gpu002 --export=ALL,MODEL_EXT=<snapshot> nonneg_multmin_test.sbatch (test 124659).
# ==== end of header ====
# GPU test of --slope-mult-min 1 (user 2026-11-12, softplus): the NFEV-5 softplus fit from D2 with SLOPE_MULT_MIN=1
# (CHECK 1, L_THREADS 32) through jobs/nonneg_fit.sbatch, then: the optimizer's lower bounds of the six free multipliers
# are 1 (fit_manifest.json parameters.lower_bounds), every final multiplier >= 1, the restart marker carries the value,
# and the start check cost equals the 5-cycle chain's c1 start (REF_FIT; multipliers 1 at the start in both).
# Settings: MODEL_EXT (required), TAG (required), REF_FIT (default the 5-cycle softplus chain's c1/fit/seg01).
#SBATCH --job-name=cops12m_nonneg_multmin
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gpus-per-node=8
#SBATCH --cpus-per-task=64
#SBATCH --mem=0
#SBATCH --exclusive
#SBATCH --time=00:45:00
#SBATCH --partition=gpu
#SBATCH --account=general
#SBATCH --exclude=gpu003
#SBATCH --output=/data001/highresolu_vt/cops_output/jobs/logs/cops12m_nonneg_multmin-%j.out
set -uo pipefail
export AMPS=/data001/highresolu_vt/AMPS APP=/data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda
export COPS_OUTPUT_ROOT=/data001/highresolu_vt/cops_output
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1
MODEL_EXT=$(readlink -f "${MODEL_EXT:?set MODEL_EXT}"); TAG=${TAG:?set TAG}
R=$COPS_OUTPUT_ROOT/runs/shap8_012m_lategrid_dz1.5; J=$COPS_OUTPUT_ROOT/jobs; T_=$R/nonneg_tests
REF_FIT=$(readlink -f "${REF_FIT:-$R/nonneg/softplus/chain_5cycles_20261112T012553Z/c1/fit/seg01}")
D=$T_/test_multmin_$TAG; [ -e "$D" ] && { echo "exists: $D" >&2; exit 2; }
mkdir -p "$D"; P=$D/provenance-$SLURM_JOB_ID; mkdir -p "$P"; cd "$APP" || exit 1
cp -p "$0" "$P/"; cp -p "$J/nonneg_fit.sbatch" "$P/"; cp -p "$MODEL_EXT"/*.py "$MODEL_EXT/SHA256SUMS" "$P/" 2>/dev/null
{ git -C "$AMPS" rev-parse HEAD; git -C "$AMPS" status --short; } > "$P/git.txt" 2>&1
lscpu > "$P/lscpu.txt" 2>&1; nvidia-smi -L > "$P/nvidia_smi_L.txt" 2>&1
echo "$(date -u +%FT%TZ) job $SLURM_JOB_ID on $(hostname) ($(grep -m1 'Model name' "$P/lscpu.txt" | sed 's/  */ /g')): MODEL_EXT=$MODEL_EXT TAG=$TAG REF_FIT=$REF_FIT"
NAME=g_softplus_multmin_$TAG
echo "env TESTNAME=$NAME MODEL=softplus START=d2 SLOPE_MULT_MIN=1 CHECK=1 NFEV=5 L_THREADS=32 MODEL_EXT=$MODEL_EXT bash $J/nonneg_fit.sbatch" > "$P/commands.sh"
env TESTNAME=$NAME MODEL=softplus START=d2 SLOPE_MULT_MIN=1 CHECK=1 NFEV=5 L_THREADS=32 MODEL_EXT=$MODEL_EXT \
    bash "$J/nonneg_fit.sbatch" > "$P/fit.log" 2>&1; rc=$?
echo "fit exit $rc"
grep -E "nonneg: (multiplier lower bound|start|model check)|start check \(start|^lambda|gpu_fit: [0-9]|Traceback|Error" "$P/fit.log" | cut -c1-400
"$APP/.venv/bin/python" - "$T_/$NAME/seg01" "$REF_FIT" > "$P/check.log" 2>&1 <<'PY' || rc=1
import json, sys
import numpy as np
seg, ref = sys.argv[1], sys.argv[2]
m = json.load(open(seg + "/fit_manifest.json")); mr = json.load(open(ref + "/fit_manifest.json"))
rec, recr = m["model_extension"]["nonneg"], mr["model_extension"]["nonneg"]
names, lo = m["parameters"]["names"], m["parameters"]["lower_bounds"]
free = [n for n in names if n.startswith("c_slope_")]
lbs = {n: lo[names.index(n)] for n in free}
z = np.load(seg + "/fit_restart.npz")
meta = json.loads(str(z["meta_json"]))
x = dict(zip(meta["names"], z["x"]))
xs = dict(zip(rec["model_check"]["layout"], rec["start_restart"]["x_start"]))
c, cr = rec["clip_core"]["start_check"]["start"]["cost"], recr["clip_core"]["start_check"]["start"]["cost"]
ok_b = all(v == 1.0 for v in lbs.values()) and all(x[n] >= 1.0 for n in free) and all(xs[n] == 1.0 for n in free)
ok_m = meta["model_extension"].get("slope_mult_min") == 1.0
ok_c = c == cr
print("lower bounds %s; start %s; final %s" % (lbs, {n: xs[n] for n in free}, {n: float(x[n]) for n in free}))
print("bound holds (lower bounds 1, start 1, every final multiplier >= 1): %s" % ("PASS" if ok_b else "FAIL"))
print("restart marker slope_mult_min = %r: %s" % (meta["model_extension"].get("slope_mult_min"), "PASS" if ok_m else "FAIL"))
print("start check cost %.16e vs the 5-cycle chain's c1 %.16e: %s (rel %.1e)" % (c, cr, "PASS bitwise" if ok_c else "FAIL",
                                                                               abs(c / cr - 1)))
jc = rec["clip_core"]["jacobian_check"]["columns"]
print("Jacobian check worst data %.2e over %d columns" % (max(v["data"] for v in jc.values() if isinstance(v, dict)), len(jc)))
sys.exit(0 if ok_b and ok_m and ok_c else 1)
PY
cat "$P/check.log"
echo "$(date -u +%FT%TZ) exit $rc"
exit $rc
