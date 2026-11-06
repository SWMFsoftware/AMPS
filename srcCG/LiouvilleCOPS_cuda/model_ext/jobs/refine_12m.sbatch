#!/bin/bash
# ==== version-controlled copy: header added 2026-11-06, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/jobs/refine_12m.sbatch (on Voyager: /data001/highresolu_vt/cops_output/jobs/refine_12m.sbatch)
# SHA256 of the run copy = SHA256 of this file without these header lines: ad0ea6f08f0cc4a042e51f51152bfc7b521fb7892790fbf7b2966329c797cf83
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: 12M mesh refined near the last observations (refine_mesh.py), its base and store (runner prepare), the bundle
#   (refined_bundle.py), the two best 12M grain fits translated (translate_grain.py) and fine-tuned there. Job 122009.
# RUN: sbatch refine_12m.sbatch (skips the steps whose products exist)
# ==== end of header ====
# 12M mesh refined near the last observations, its base, and the 12M grain fine tunings redone on it (user request
# 2026-11-04: test whether a source-rate gradient across a surface element improves the fit; answers: edges halved
# (1 -> 4), the elements first seen in the 12M fine-tuning window + a 10 m margin, the two best grain fits, run now).
# 1. model_ext/refine_mesh.py (conforming red-green refinement, geometry unchanged) -> $R/mesh/
# 2. runner prepare group with --fits B29 (base + record store only; no visibility, no gravity: the templates are the
#    parent's, step 3), the 12M campaign settings, --allow-genus 4 (topology unchanged)
# 3. model_ext/refined_bundle.py: templates of the 12M bundles carried to the children (slope from the 12M mesh)
# 4. per fit: model_ext/translate_grain.py (the 12M fit's restart -> the refined base: coefficients, T, kappa, grain;
#    templates from the step-3 bundle; predictions on the refined base; the fit's selection), then
#    grain_finetune.py (5 min, 5 intervals, cpu, mixed), as on 12M
# 5. comparison of the translated predictions with the 12M fit's (geometry unchanged: they should agree closely)
#SBATCH --job-name=cops12m_refine
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gpus-per-node=8
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=06:00:00
#SBATCH --partition=gpu
#SBATCH --account=general
#SBATCH --output=/data001/highresolu_vt/cops_output/jobs/logs/cops12m_refine-%j.out
set -uo pipefail
export AMPS=/data001/highresolu_vt/AMPS APP=/data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda
export T=/data001/highresolu_vt/cops_cluster_transfer_20261026 COPS_OUTPUT_ROOT=/data001/highresolu_vt/cops_output
export COPS_DATA_ROOT=$T/rosetta_data_end_of_mission_cops
export CUDA_PATH=/usr/local/cuda-12.9 CUDA_HOME=/usr/local/cuda-12.9 CUPY_CACHE_DIR=$COPS_OUTPUT_ROOT/.cupy
export PATH=$CUDA_HOME/bin:/usr/mpi/gcc/openmpi-4.1.7rc1/bin:$PATH
export LD_LIBRARY_PATH=$CUDA_HOME/lib64${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1
R12=$COPS_OUTPUT_ROOT/runs/shap8_012m_flat_earlier; R=$COPS_OUTPUT_ROOT/runs/shap8_012m_refined_late
X=$COPS_OUTPUT_ROOT/model_ext; PY="$APP/.venv/bin/python -W ignore"
BDF12=$COPS_DATA_ROOT/shape_models/SHAP8_MIARD/cg-spg_mspcd-shap8-v0.4cmb-12Mfacets.bdf
FITS="C30_n2000_grain_cross10 NS31_n2000_mono10_ordered_grain_cross10"
FT12=$R12/grain/finetune/C30_n2000_grain_cross10_m5_i5_cpu
P=$R/provenance/refine-$SLURM_JOB_ID; mkdir -p "$P"; cd "$APP" || exit 1
cp -p "$0" "$X"/refine_mesh.py "$X"/refined_bundle.py "$X"/translate_grain.py "$X"/grain_finetune.py "$X"/grain_crossing.py "$P/"
sha256sum "$X"/*.py > "$P/model_ext.sha256"
{ git -C "$AMPS" rev-parse HEAD; git -C "$AMPS" status --short; } > "$P/git.txt" 2>&1
module list > "$P/module_list.txt" 2>&1 || true
nvidia-smi -q > "$P/nvidia_smi_q.txt" 2>&1
{ nvcc --version; gcc --version | head -1; grep -m1 CSPICE /home/vtenishev/SPICE/cspice/doc/version.txt; } > "$P/toolchain.txt" 2>&1
~/.local/bin/uv pip list --python "$APP/.venv/bin/python" > "$P/pip_list.txt" 2>&1
sha256sum "$APP/native/liouville_cops_cuda" > "$P/executables.sha256"
env | grep -E '^(CUDA_|CUPY_|COPS_|OMP_|OPENBLAS_|MKL_|PYTHONUNBUFFERED|LD_LIBRARY_PATH|PATH)=' | sort > "$P/environment.txt"
DEV=$(nvidia-smi --query-gpu=index --format=csv,noheader | paste -sd, -)
M=$R/mesh; SHAPE=$M/$(basename "$BDF12" .bdf)_refined.bdf
if [ ! -e "$M/refine.json" ]; then
    echo "$(date -u +%FT%TZ) refine: start"
    $PY "$X/refine_mesh.py" --bdf "$BDF12" --centroids "$R12/grain/fits/C30_n2000_grain_cross10/fit_restart.npz" \
        --first-interval "$FT12/fine_tuning_prepared.npz" --margin 10 --out "$M" > "$P/refine.log" 2>&1; rc=$?
    echo "$(date -u +%FT%TZ) refine: exit $rc"; tail -16 "$P/refine.log"; [ $rc -eq 0 ] || exit $rc
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
    $PY "$X/refined_bundle.py" --base "$R/base" --parent "$M/parent.npy" --bundle "$R12/bundle/restart_templates.npz" \
        --bundle "$R12/c16/bundle/slope_templates.npz" --source-restart "$R12/grain/fits/C30_n2000_grain_cross10/fit_restart.npz" \
        --out "$B" > "$P/bundle.log" 2>&1 || { cat "$P/bundle.log"; exit 1; }
    cat "$P/bundle.log"
fi
export COPS_CUDA_RECORD_STORE=$R/store
failed=0
for FIT in $FITS; do
    XD=$R/grain/translated/${FIT}_refined; FD=$R/grain/translated/finetune/${FIT}_refined_m5_i5_cpu
    mkdir -p "$R/grain/translated/finetune"
    if [ ! -e "$XD/fit_restart.npz" ]; then
        [ -e "$XD" ] && mv "$XD" "$XD.incomplete_$(date -u +%Y%m%dT%H%M%SZ)"
        echo "$(date -u +%FT%TZ) translate $FIT: start"
        $PY "$X/translate_grain.py" --restart "$R12/grain/fits/$FIT/fit_restart.npz" --to-run "$R/base" --target-bundle "$B" \
            --out "$XD" --devices "$DEV" --precision mixed > "$P/translate_$FIT.log" 2>&1; rc=$?
        echo "$(date -u +%FT%TZ) translate $FIT: exit $rc"; tail -1 "$P/translate_$FIT.log" | cut -c1-250
        [ $rc -eq 0 ] || { failed=1; continue; }
    fi
    # the translated predictions against the 12M fit's (same observation points)
    $PY - "$R12/grain/fits/$FIT/fit_restart.npz" "$XD/fit_restart.npz" <<'PY' | tee "$P/compare_$FIT.txt"
import sys, numpy as np
a, b = np.load(sys.argv[1]), np.load(sys.argv[2])
for g in ("ng", "rg"):
    pa, pb = a["pred_" + g], b["pred_" + g]
    rel = np.abs(pb / pa - 1)
    print("%s: refined vs 12M prediction, max rel %.3e at point %d, median %.3e, last 95 points max %.3e" % (
        g.upper(), rel.max(), int(a["point_ids"][rel.argmax()]), np.median(rel), rel[-95:].max()))
PY
    if [ -e "$FD" ]; then echo "$FIT: fine tuning exists, skipped"; continue; fi
    echo "$(date -u +%FT%TZ) fine tuning $FIT: start"
    $PY "$X/grain_finetune.py" "$R/base" --restart-from "$XD/fit_restart.npz" --minutes 5 --intervals 5 \
        --method cpu --max-nfev 2000 --devices "$DEV" --precision mixed --out "$FD" > "$P/finetune_$FIT.log" 2>&1; rc=$?
    echo "$(date -u +%FT%TZ) fine tuning $FIT: exit $rc"; tail -4 "$P/finetune_$FIT.log" | cut -c1-250
    [ $rc -eq 0 ] || failed=1
done
exit $failed
