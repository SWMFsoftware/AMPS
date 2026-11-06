#!/bin/bash
# ==== version-controlled copy: header added 2026-11-06, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/jobs/refine_20m.sbatch (on Voyager: /data001/highresolu_vt/cops_output/jobs/refine_20m.sbatch)
# SHA256 of the run copy = SHA256 of this file without these header lines: 36c311a6440ad5190f8fe3e4efb76861b9f697366b3eecdb28cb59bf747ab187
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: 20M mesh refined over the surface seen after the gap (seen_after_gap.py + refine_mesh.py), base and store,
#   bundle, the 12M T75 = Tmin = 133 K fit translated and fine-tuned there. Job 122011.
# RUN: sbatch refine_20m.sbatch
# ==== end of header ====
# 20M mesh refined over the entire surface seen by NG or RG after the gap, its base, and the 12M fixed-T grain fit
# fine-tuned on it (user request 2026-11-04: "do the same for 20M ... over the entire part of the surface seen after
# the gap by either NG or RG ... transfer the best fit from 12M, and run the final fitting phase with the locally
# refined 20M"; answers: 1 level (edges halved, area /4), the 12M fit with T75 = Tmin = 133 K fixed).
# 0. model_ext/seen_after_gap.py (done before this job, srun on gpu001): $R/seen_after_gap.npz, 1,466,454 facets
#    with a record at any of the 95 points from 10:33:10 UTC (points 1055-1149) in the 20M record store
# 1. model_ext/refine_mesh.py --seed-mask, margin 0 (conforming red-green refinement, geometry unchanged) -> $R/mesh/
# 2. runner prepare group with --fits B29 (base + record store only; templates from step 3), the 20M campaign
#    settings, --allow-genus 4 (topology unchanged)
# 3. model_ext/refined_bundle.py: templates of the 20M bundles carried to the children (gravity provenance kept)
# 4. model_ext/translate_grain.py (the 12M fixed-T fit's restart -> the refined base, as job 121974 did for 20M;
#    cliff template from the step-3 bundle), then grain_finetune.py (5 min, 5 intervals, cpu, mixed, max_nfev 2000)
#SBATCH --job-name=cops20m_refine
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gpus-per-node=8
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=06:00:00
#SBATCH --partition=gpu
#SBATCH --account=general
#SBATCH --output=/data001/highresolu_vt/cops_output/jobs/logs/cops20m_refine-%j.out
set -uo pipefail
export AMPS=/data001/highresolu_vt/AMPS APP=/data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda
export T=/data001/highresolu_vt/cops_cluster_transfer_20261026 COPS_OUTPUT_ROOT=/data001/highresolu_vt/cops_output
export COPS_DATA_ROOT=$T/rosetta_data_end_of_mission_cops
export CUDA_PATH=/usr/local/cuda-12.9 CUDA_HOME=/usr/local/cuda-12.9 CUPY_CACHE_DIR=$COPS_OUTPUT_ROOT/.cupy
export PATH=$CUDA_HOME/bin:/usr/mpi/gcc/openmpi-4.1.7rc1/bin:$PATH
export LD_LIBRARY_PATH=$CUDA_HOME/lib64${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1
R12=$COPS_OUTPUT_ROOT/runs/shap8_012m_flat_earlier; R20=$COPS_OUTPUT_ROOT/runs/shap8_020m_flat_earlier
R=$COPS_OUTPUT_ROOT/runs/shap8_020m_refined_late
X=$COPS_OUTPUT_ROOT/model_ext; PY="$APP/.venv/bin/python -W ignore"
BDF20=$COPS_DATA_ROOT/shape_models/SHAP8_MIARD/cg-spg_mspcd-shap8-v0.4cmb-20Mfacets.bdf
SRC=C30_n2000_grain_cross10_T75_133K_Tmin_133K; SF=$R12/grain/fixed_T/$SRC
C20=$R20/grain/translated/${SRC}_from_12M/fit_restart.npz       # a 20M restart (centroids, areas)
P=$R/provenance/refine-$SLURM_JOB_ID; mkdir -p "$P"; cd "$APP" || exit 1
cp -p "$0" "$X"/seen_after_gap.py "$X"/refine_mesh.py "$X"/refined_bundle.py "$X"/translate_grain.py "$X"/grain_finetune.py \
    "$X"/grain_crossing.py "$P/"
sha256sum "$X"/*.py > "$P/model_ext.sha256"
{ git -C "$AMPS" rev-parse HEAD; git -C "$AMPS" status --short; } > "$P/git.txt" 2>&1
module list > "$P/module_list.txt" 2>&1 || true
nvidia-smi -q > "$P/nvidia_smi_q.txt" 2>&1
{ nvcc --version; gcc --version | head -1; grep -m1 CSPICE /home/vtenishev/SPICE/cspice/doc/version.txt; } > "$P/toolchain.txt" 2>&1
~/.local/bin/uv pip list --python "$APP/.venv/bin/python" > "$P/pip_list.txt" 2>&1
sha256sum "$APP/native/liouville_cops_cuda" > "$P/executables.sha256"
env | grep -E '^(CUDA_|CUPY_|COPS_|OMP_|OPENBLAS_|MKL_|PYTHONUNBUFFERED|LD_LIBRARY_PATH|PATH)=' | sort > "$P/environment.txt"
DEV=$(nvidia-smi --query-gpu=index --format=csv,noheader | paste -sd, -)
[ -e "$R/seen_after_gap.npz" ] || { echo "missing $R/seen_after_gap.npz"; exit 1; }
M=$R/mesh; SHAPE=$M/$(basename "$BDF20" .bdf)_refined.bdf
if [ ! -e "$M/refine.json" ]; then
    echo "$(date -u +%FT%TZ) refine: start"
    $PY "$X/refine_mesh.py" --bdf "$BDF20" --centroids "$C20" --seed-mask "$R/seen_after_gap.npz" --margin 0 \
        --out "$M" > "$P/refine.log" 2>&1; rc=$?
    echo "$(date -u +%FT%TZ) refine: exit $rc"; tail -17 "$P/refine.log"; [ $rc -eq 0 ] || exit $rc
fi
OPTS=(--shape "$SHAPE" --surface-bc flat --max-nfev 2000 --allow-genus 4 --fits B29 --precision mixed
      "--fit-args=--profile all_available --degree 4 --enhancement-degree 4 --lambda 0 --jobs 1 --tbounds 30 400 --ftol 1e-10 --xtol 1e-10 --gtol 1e-10 --check-tolerance 1e-6"
      --output-root "$COPS_OUTPUT_ROOT" --fine-tune --ft-minutes 5 --ft-intervals 5 --ft-method cpu --ft-max-nfev 2000)
echo "$(date -u +%FT%TZ) runner prepare (base + store): start"
bash -c 'echo 1000 > /proc/self/oom_score_adj; exec "$@"' _ $PY tools/cluster/cops_pipeline.py run --restart-dir "$R" \
    --group prepare --wall-time slurm "${OPTS[@]}" > "$P/runner_prepare.log" 2>&1; rc=$?
echo "$(date -u +%FT%TZ) runner prepare: exit $rc"; grep -v "^\[base\] .*point " "$P/runner_prepare.log" | tail -6; [ $rc -eq 0 ] || exit $rc
B=$R/bundle_parent/templates_parent.npz
if [ ! -e "$B" ]; then
    mkdir -p "$(dirname "$B")"
    $PY "$X/refined_bundle.py" --base "$R/base" --parent "$M/parent.npy" --bundle "$R20/bundle/restart_templates.npz" \
        --bundle "$R20/c16/bundle/slope_templates.npz" --source-restart "$C20" --out "$B" > "$P/bundle.log" 2>&1 \
        || { cat "$P/bundle.log"; exit 1; }
    cat "$P/bundle.log"
fi
export COPS_CUDA_RECORD_STORE=$R/store
N=${SRC}_from_12M; XD=$R/grain/translated/$N; FD=$R/grain/translated/finetune/${N}_m5_i5_cpu
mkdir -p "$R/grain/translated/finetune"
if [ ! -e "$XD/fit_restart.npz" ]; then
    [ -e "$XD" ] && mv "$XD" "$XD.incomplete_$(date -u +%Y%m%dT%H%M%SZ)"
    echo "$(date -u +%FT%TZ) translate $SRC: start"
    $PY "$X/translate_grain.py" --restart "$SF/fit_restart.npz" --to-run "$R/base" --target-bundle "$B" \
        --out "$XD" --devices "$DEV" --precision mixed > "$P/translate.log" 2>&1; rc=$?
    echo "$(date -u +%FT%TZ) translate $SRC: exit $rc"; tail -1 "$P/translate.log" | cut -c1-250
    [ $rc -eq 0 ] || exit $rc
fi
# the translated predictions against the unrefined 20M translation (job 121974) and the 12M fit (same points)
for REF in "$R20/grain/translated/$N/fit_restart.npz" "$SF/fit_restart.npz"; do
$PY - "$REF" "$XD/fit_restart.npz" <<'PY' | tee -a "$P/compare.txt"
import sys, numpy as np
a, b = np.load(sys.argv[1]), np.load(sys.argv[2])
print("reference", sys.argv[1])
for g in ("ng", "rg"):
    pa, pb = a["pred_" + g], b["pred_" + g]
    rel = np.abs(pb / pa - 1)
    print("%s: refined 20M vs reference prediction, max rel %.3e at point %d, median %.3e, last 95 points max %.3e" % (
        g.upper(), rel.max(), int(a["point_ids"][rel.argmax()]), np.median(rel), rel[-95:].max()))
PY
done
[ -e "$FD" ] && { echo "fine tuning exists, skipped"; exit 0; }
echo "$(date -u +%FT%TZ) fine tuning $SRC: start"
$PY "$X/grain_finetune.py" "$R/base" --restart-from "$XD/fit_restart.npz" --minutes 5 --intervals 5 \
    --method cpu --max-nfev 2000 --devices "$DEV" --precision mixed --out "$FD" > "$P/finetune.log" 2>&1; rc=$?
echo "$(date -u +%FT%TZ) fine tuning $SRC: exit $rc"; tail -4 "$P/finetune.log" | cut -c1-250
exit $rc
