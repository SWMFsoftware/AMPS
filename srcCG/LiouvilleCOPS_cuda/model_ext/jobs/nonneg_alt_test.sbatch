#!/bin/bash
# ==== version-controlled copy: header added 2026-11-14, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/jobs/nonneg_alt_test.sbatch (on Voyager: /data001/highresolu_vt/cops_output/jobs/nonneg_alt_test.sbatch)
# SHA256 of the run copy = SHA256 of this file without these header lines: 5bfb2f3c7ee7aff381c54710b3badde29d2a29e97cccf5fb3bd02389445597a6
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: job. GPU test of one alternative model (A, B, D): start cost vs an unpatched independent evaluation, then
#   5 evaluations with the full Jacobian check, L_THREADS 32, SLOPE_MULT_MIN 1.
# RUN: sbatch --nodelist=gpu002 --export=ALL,MODEL_EXT=<snapshot>,TAG=<tag>,MODEL=sos|background|logspread
#   nonneg_alt_test.sbatch (tests 124998 sos, 124999 background, 125028 logspread).
# ==== end of header ====
# GPU tests of the alternative non-negative models (options A, B, D; 2026-11-14), one MODEL per job, one node:
#   (1) TESTMODE=independent: the driver's start check vs an UNPATCHED process evaluating the independently computed
#       activity (test_nonneg.py --gpu-independent; rel <= 1e-12);
#   (2) the short fit (NFEV 5, CHECK 1: the full Jacobian check of the default parameters incl. the extra block) with
#       L_THREADS 32 and SLOPE_MULT_MIN 1 (the common scheme), the starts gram / d2 / uniform.
# Settings: MODEL_EXT (frozen snapshot, required), TAG (required), MODEL (sos | background | logspread), SOS_RANK (3),
#   BG_MIN (1e-2), SPREAD_F (1), SPREAD_CREF (4.2688066575e-03), START_FILE (sos; default .../sos_start_r<r>.npz).
#SBATCH --job-name=cops12m_nonneg_alt_test
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gpus-per-node=8
#SBATCH --cpus-per-task=64
#SBATCH --mem=0
#SBATCH --exclusive
#SBATCH --time=01:00:00
#SBATCH --partition=gpu
#SBATCH --account=general
#SBATCH --exclude=gpu003
#SBATCH --output=/data001/highresolu_vt/cops_output/jobs/logs/cops12m_nonneg_alt_test-%j.out
set -uo pipefail
export AMPS=/data001/highresolu_vt/AMPS APP=/data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda
export COPS_OUTPUT_ROOT=/data001/highresolu_vt/cops_output
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1
MODEL_EXT=$(readlink -f "${MODEL_EXT:?set MODEL_EXT}"); TAG=${TAG:?set TAG}; MODEL=${MODEL:?set MODEL}
R=$COPS_OUTPUT_ROOT/runs/shap8_012m_lategrid_dz1.5; J=$COPS_OUTPUT_ROOT/jobs; T_=$R/nonneg_tests
case $MODEL in
    sos) OPTS=(SOS_RANK=${SOS_RANK:-3} START_FILE=$(readlink -f "${START_FILE:-$R/nonneg/bundle/sos_start_r${SOS_RANK:-3}.npz}"));;
    background) OPTS=(BG_MIN=${BG_MIN:-1e-2});;
    logspread) OPTS=(SPREAD_F=${SPREAD_F:-1} SPREAD_CREF=${SPREAD_CREF:-4.2688066575e-03});;
    *) echo "MODEL must be sos, background or logspread" >&2; exit 2;;
esac
D=$T_/alt_test_${MODEL}_$TAG; [ -e "$D" ] && { echo "exists: $D" >&2; exit 2; }
mkdir -p "$D"; P=$D/provenance-$SLURM_JOB_ID; mkdir -p "$P"; cd "$APP" || exit 1
cp -p "$0" "$P/"; cp -p "$J/nonneg_fit.sbatch" "$P/"; cp -p "$MODEL_EXT"/*.py "$MODEL_EXT/SHA256SUMS" "$P/" 2>/dev/null
{ git -C "$AMPS" rev-parse HEAD; git -C "$AMPS" status --short; } > "$P/git.txt" 2>&1
lscpu > "$P/lscpu.txt" 2>&1; nvidia-smi -L > "$P/nvidia_smi_L.txt" 2>&1
{ echo '#!/bin/bash'; echo "# job $SLURM_JOB_ID on $(hostname), $(date -u +%FT%TZ)"; grep '^export ' "$0"; echo "cd $APP"; } > "$P/commands.sh"
echo "$(date -u +%FT%TZ) job $SLURM_JOB_ID on $(hostname): MODEL=$MODEL ${OPTS[*]} MODEL_EXT=$MODEL_EXT TAG=$TAG"
rc=0
run() {
    local label=$1; shift
    echo "$(date -u +%FT%TZ) $label"
    echo "env $* # $label" >> "$P/commands.sh"
    env "$@" > "$P/$label.log" 2>&1; local r=$?
    echo "  exit $r"
    grep -E "nonneg: (start|model check|alternative|multiplier lower)|start check \(start|^lambda|gpu_fit: [0-9]|PASS|FAIL|Traceback|Error" "$P/$label.log" | cut -c1-420
    return $r
}
run indep TESTNAME=g_alt_${MODEL}_indep_$TAG TESTMODE=independent MODEL=$MODEL "${OPTS[@]}" SLOPE_MULT_MIN=1 CHECK=0 NFEV=5 \
    L_THREADS=32 MODEL_EXT=$MODEL_EXT bash "$J/nonneg_fit.sbatch" || rc=1
run short TESTNAME=g_alt_${MODEL}_short_$TAG MODEL=$MODEL "${OPTS[@]}" SLOPE_MULT_MIN=1 CHECK=1 NFEV=5 L_THREADS=32 \
    MODEL_EXT=$MODEL_EXT bash "$J/nonneg_fit.sbatch" || rc=1
"$APP/.venv/bin/python" - "$T_/g_alt_${MODEL}_short_$TAG/seg01" > "$P/summary.log" 2>&1 <<'PY' || rc=1
import json, sys
import numpy as np
seg = sys.argv[1]
m = json.load(open(seg + "/fit_manifest.json")); rec = m["model_extension"]["nonneg"]
jc = rec["clip_core"]["jacobian_check"]["columns"]
ev = [l.split() for l in open(seg + "/evaluations.txt") if not l.startswith("#")]
fun = [float(r[2]) for r in ev if r[1] == "fun"]; jac = [float(r[2]) for r in ev if r[1] == "jac"]
z = np.load(seg + "/nonneg_facets.npz"); A, seen = z["A"], z["seen"].astype(bool)
alt = rec["model_check"].get("alternative", {})
print("start cost %.15e (data %.15e); cost after 5 evaluations %s" % (rec["clip_core"]["start_check"]["start"]["cost"],
      rec["clip_core"]["start_check"]["start"]["data_cost"], m["optimizer_results"][0]["cost"]))
print("Jacobian check: worst data %.2e over %d columns: %s" % (max(v["data"] for v in jc.values() if isinstance(v, dict)),
      len(jc), {k: "%.1e" % v["data"] for k, v in jc.items() if isinstance(v, dict) and v["data"] > 1e-6}))
print("per evaluation: fun median %.3f s, jac median %.3f s; A >= 0 %s; A > 0 on %.4f of the seen facets" % (
      np.median(fun), np.median(jac), bool(np.all(A >= 0)), float(np.mean(A[seen] > 0))))
if alt.get("spread"):
    sp = alt["spread"]
    print("spread: f %r C_ref %r; at the start penalty %r Var_w(u) %r; final %s" % (sp.get("f"), sp.get("C_ref"),
          sp.get("penalty_cost_at_start"), sp.get("var_w_u_at_start"), sp.get("final")))
if alt.get("uniform_start"):
    u = alt["uniform_start"]
    print("uniform start: s %.10e (ln s %.10f), fit-weighted alternative s %.6e; rms d ln p %s; linearity %.1e" % (
          u["s"], u["ln_s"], u["fit_weighted_alternative"]["s"],
          {k: (v["rms_dln"] if v["rms_dln"] is None else round(v["rms_dln"], 5)) for k, v in u["dln_stats"].items()},
          u["linearity_check_max_rel"]))
PY
cat "$P/summary.log"
echo "$(date -u +%FT%TZ) exit $rc"
exit $rc
