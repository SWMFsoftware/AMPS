#!/bin/bash
# ==== version-controlled copy: header added 2026-11-08, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/jobs/positive_held_test.sbatch (on Voyager: /data001/highresolu_vt/cops_output/jobs/positive_held_test.sbatch)
# SHA256 of the run copy = SHA256 of this file without these header lines: 16bcabad1675dfb4e905db5dac2afc4c49619d7443f324df61a6144e7ad444c1
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: job (test). GPU tests of the held term of positive_fit.py (span check, held_next + grain from file + start
#   restart, resume); jobs 122445, 122446. user decisions 2026-11-07/08.
# ==== end of header ====
# GPU tests of the held exponent term (model_ext/positive_fit.py --held-exponent; 2026-11-08) on the 12M lategrid base,
# through jobs/positive_fit.sbatch (TESTNAME), RATIO 0.002, LAMBDA 1e-2, NFEV 4, CHECK 1, START_RESTART = the first
# scan's lambda 1e-2 restart (runs/shap8_012m_lategrid_dz1.5/positive/C30_n4000_exp_cliff_w10x5_grain_cross10_lam1e-2/
# seg01/fit_restart.npz). Products: runs/shap8_012m_lategrid_dz1.5/positive_tests/p2_<MODE>_<UTC>/ (inputs/ = the
# synthetic held files, seg01[, seg02], provenance-<job>/).
#   MODE=span    held h = Y z_Y + z_T T_cliff (z_Y ~ N(0, 0.02), z_T = 0.1; in the span of the fitted columns) with
#                --check-held-span: the data residuals and Jacobian with h at the start x must equal those without h at
#                x + h_B (both evaluated on the GPUs), D(x) == D(x + h_B); Jacobian check; D direct vs fitter.
#   MODE=cycle2  a synthetic fine-tuning state file (held_next = 0.3 N(0,1) on 5 % of the facets, zeros elsewhere,
#                mesh_identity; grain_A_uPa 7.5, grain_B_uPa_per_100s 0.1) as HELD and GRAIN_FROM with START_RESTART
#                (cycle 2 of jobs/positive_cycles.sh); then a resume of that fit (positive_fit.py --restart-from seg01,
#                the same held file, --max-nfev 2, --segment 2): its start check must reproduce seg01's final cost.
#   sbatch --export=ALL,MODE=span positive_held_test.sbatch ; sbatch --export=ALL,MODE=cycle2 positive_held_test.sbatch
#SBATCH --job-name=cops12m_heldtest
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gpus-per-node=8
#SBATCH --cpus-per-task=64
#SBATCH --mem=0
#SBATCH --exclusive
#SBATCH --time=01:00:00
#SBATCH --partition=gpu
#SBATCH --account=general
#SBATCH --output=/data001/highresolu_vt/cops_output/jobs/logs/cops12m_heldtest-%j.out
set -uo pipefail
export AMPS=/data001/highresolu_vt/AMPS APP=/data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda
export T=/data001/highresolu_vt/cops_cluster_transfer_20261026 COPS_OUTPUT_ROOT=/data001/highresolu_vt/cops_output
export COPS_DATA_ROOT=$T/rosetta_data_end_of_mission_cops
export CUDA_PATH=/usr/local/cuda-12.9 CUDA_HOME=/usr/local/cuda-12.9 CUPY_CACHE_DIR=$COPS_OUTPUT_ROOT/.cupy
export PATH=$CUDA_HOME/bin:/usr/mpi/gcc/openmpi-4.1.7rc1/bin:$PATH
export LD_LIBRARY_PATH=$CUDA_HOME/lib64${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1
MODE=${MODE:-}
case $MODE in span|cycle2) ;; *) echo "MODE=span or MODE=cycle2" >&2; exit 2;; esac
R=$COPS_OUTPUT_ROOT/runs/shap8_012m_lategrid_dz1.5; PY="$APP/.venv/bin/python -W ignore"
BASE=$R/base; PTS=$R/late_grid_points.txt
BUNDLE=$COPS_OUTPUT_ROOT/runs/shap8_012m_flat_earlier/bundle/restart_templates.npz
SRC=$R/positive/C30_n4000_exp_cliff_w10x5_grain_cross10_lam1e-2/seg01/fit_restart.npz
TN=p2_${MODE}_$(date -u +%Y%m%dT%H%M%SZ); TD=$R/positive_tests/$TN; IN=$TD/inputs
mkdir -p "$IN"; cp -p "$0" "$TD/positive_held_test.sbatch"
echo "$(date -u +%FT%TZ) job $SLURM_JOB_ID on $(hostname): held test $MODE in $TD"
cd "$APP" || exit 1
# ---- synthetic held files (12M mesh order; the cliff column of the bundle as the fitter selects it)
$PY - "$MODE" "$IN" "$BASE" "$BUNDLE" <<'PY' 2>&1 | tee "$IN/make_inputs.log"
import json, sys
import numpy as np
from pathlib import Path
sys.path[:0] = ["/data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda"]
from cops_cuda import reference
reference.bind()
fit_cops, at = reference.load("fit_cops"), reference.load("activity_templates")
mode, out, base, bundle = sys.argv[1], Path(sys.argv[2]), sys.argv[3], sys.argv[4]
el = fit_cops.read_elements(base + "/elements.bin")
mesh = at.mesh_identity(el["centroid"], el["area"])
nf = len(el["area"])
rng = np.random.default_rng(20261108)
if mode == "span":
    T = np.ascontiguousarray(at.select_templates(at.load_bundle(bundle, el["centroid"], el["area"]), ["cliff_ge_60deg"]))
    zY, zT = rng.normal(0.0, 0.02, 25), 0.1
    h = el["Y"] @ zY + zT * T[:, 0]
    np.savez(out / "held_span.npz", held=h, z_Y=zY, z_T=np.float64(zT), mesh_identity=np.array(mesh))
    print("held_span.npz: h = Y z_Y + %.3f T_cliff, rms %.4f, range [%.4f, %.4f]; z_Y %s"
          % (zT, float(np.sqrt(np.mean(h * h))), h.min(), h.max(), np.round(zY, 4).tolist()))
else:
    h = np.zeros(nf)
    idx = rng.choice(nf, nf // 20, replace=False)
    h[idx] = 0.3 * rng.normal(size=len(idx))
    meta = {"kind": "synthetic fine-tuning state for the held-term GPU test", "mesh_identity": mesh}
    np.savez(out / "finetune_state.npz", held_next=h, grain_A_uPa=np.float64(7.5),
             grain_B_uPa_per_100s=np.float64(0.1), meta_json=np.array(json.dumps(meta)))
    print("finetune_state.npz: held_next on %d of %d facets (0.3 N(0,1)), grain 7.5 uPa, 0.1 uPa/100 s"
          % (len(idx), nf))
PY
[ ${PIPESTATUS[0]} -eq 0 ] || { echo "making the inputs failed" >&2; exit 1; }
sha256sum "$IN"/*.npz | tee "$IN/sha256.txt"
FIT=$COPS_OUTPUT_ROOT/jobs/positive_fit.sbatch
if [ "$MODE" = span ]; then
    TESTNAME=$TN RATIO=0.002 LAMBDA=1e-2 NFEV=4 CHECK=1 START_RESTART=$SRC HELD=$IN/held_span.npz CHECK_HELD_SPAN=1 \
        bash "$FIT"; rc=$?
    echo "$(date -u +%FT%TZ) span test fit: exit $rc"; exit $rc
fi
TESTNAME=$TN RATIO=0.002 LAMBDA=1e-2 NFEV=4 CHECK=1 START_RESTART=$SRC HELD=$IN/finetune_state.npz \
    GRAIN_FROM=$IN/finetune_state.npz bash "$FIT"; rc=$?
echo "$(date -u +%FT%TZ) cycle2 test fit (seg01): exit $rc"
[ $rc -eq 0 ] || exit $rc
# ---- resume of that fit (the gpu_fit --restart-from path), same held term: must reproduce seg01's final cost
export COPS_POS_ENABLE=1 COPS_POS_NIGHT_RATIO=0.002 COPS_POS_ACTIVITY=exp COPS_POS_FLOOR_RATIO=0.002
P=$TD/provenance-$SLURM_JOB_ID
if [ -n "${CUDA_VISIBLE_DEVICES:-}" ]; then NG=$(echo "$CUDA_VISIBLE_DEVICES" | tr , "\n" | wc -l); else NG=$(nvidia-smi -L | wc -l); fi
DEV=$(seq -s, 0 $((NG - 1)))
CMD=($PY "$P/positive_fit.py" --point-table "$PTS" --ratio 0.002 --restart-from "$TD/seg01/fit_restart.npz" \
    --held-exponent "$IN/finetune_state.npz" --crossing-weight 10 --devices "$DEV" --precision mixed --max-nfev 2 \
    --segment 2 --wall-time slurm --template-start 0 -- "$BASE" --profile all_available --degree 4 \
    --enhancement-degree 4 --lambda 1e-2 --jobs 1 --fix-temperature --tbounds 30 400 --ftol 1e-10 --xtol 1e-10 \
    --gtol 1e-10 --check-tolerance inf --out "$TD/seg02" --activity-template "$BUNDLE" cliff_ge_60deg)
{ echo "export COPS_POS_ENABLE=1 COPS_POS_NIGHT_RATIO=0.002 COPS_POS_ACTIVITY=exp COPS_POS_FLOOR_RATIO=0.002"
  echo "cd $APP"; printf '%q ' "${CMD[@]}"; echo; } > "$P/commands_resume_seg02.sh"
bash -c 'echo 1000 > /proc/self/oom_score_adj; exec "$@"' _ "${CMD[@]}" > "$P/resume_seg02.log" 2>&1; rc=$?
echo "$(date -u +%FT%TZ) cycle2 test resume (seg02): exit $rc"
grep -E '^(positive_fit|initial cost|lambda |gpu_fit: )' "$P/resume_seg02.log" | cut -c1-400
exit $rc
