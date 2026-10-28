#!/bin/bash
# Fine tuning of the 12M grain fits with the dust-grain term refitted (user 2026-10-27; model_ext/grain_finetune.py
# around the unchanged tools/fine_tune_gpu.py --method cpu, the runner's settings: --minutes 5 --intervals 5
# --max-nfev 2000, mixed precision with the record store, all GPUs of the node).
#   sbatch --dependency=afterok:<fit jobs> --export=ALL,COPS_RUN=<runner restart dir> grain_finetune.sbatch
# Tunes every complete fit in <run>/grain/fits/ (both sets); products <run>/grain/finetune/<fit>_m5_i5_cpu/.
# A finished one (fine_tuning.json state complete) is skipped; an unfinished one is moved aside as *.incomplete_<UTC>.
#SBATCH --job-name=cops_grain_ft
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gpus-per-node=8
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=04:00:00
#SBATCH --partition=gpu
#SBATCH --account=general
#SBATCH --output=/data001/highresolu_vt/cops_output/model_ext/logs/cops_grain_ft-%j.out
set -uo pipefail
export AMPS=/data001/highresolu_vt/AMPS APP=/data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda
export T=/data001/highresolu_vt/cops_cluster_transfer_20261026 COPS_OUTPUT_ROOT=/data001/highresolu_vt/cops_output
export COPS_DATA_ROOT=$T/rosetta_data_end_of_mission_cops
export CUDA_PATH=/usr/local/cuda-12.9 CUDA_HOME=/usr/local/cuda-12.9 CUPY_CACHE_DIR=$COPS_OUTPUT_ROOT/.cupy
export PATH=$CUDA_HOME/bin:/usr/mpi/gcc/openmpi-4.1.7rc1/bin:$PATH
export LD_LIBRARY_PATH=$CUDA_HOME/lib64${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1
R=${COPS_RUN:?set COPS_RUN}; X=$COPS_OUTPUT_ROOT/model_ext
export COPS_CUDA_RECORD_STORE=$R/store
PY="$APP/.venv/bin/python -W ignore"
G=$R/grain; F=$G/fits; O=$G/finetune; L=$G/logs; P=$G/provenance-ft-$SLURM_JOB_ID
mkdir -p "$O" "$L" "$P"
cd "$APP" || exit 1
cp -p "$0" "$P/grain_finetune.sbatch"; cp -p "$X"/grain_crossing.py "$X"/grain_finetune.py "$P/"; sha256sum "$X"/*.py > "$P/model_ext.sha256"
{ git -C "$AMPS" rev-parse HEAD; git -C "$AMPS" status --short; } > "$P/git.txt" 2>&1
module list > "$P/module_list.txt" 2>&1 || true
nvidia-smi -q > "$P/nvidia_smi_q.txt" 2>&1
{ nvcc --version; gcc --version | head -1; } > "$P/toolchain.txt" 2>&1
~/.local/bin/uv pip list --python "$APP/.venv/bin/python" > "$P/pip_list.txt" 2>&1
env | grep -E '^(CUDA_|CUPY_|COPS_|OMP_|OPENBLAS_|MKL_|PYTHONUNBUFFERED|LD_LIBRARY_PATH|PATH)=' | sort > "$P/environment.txt"
DEV=$(nvidia-smi --query-gpu=index --format=csv,noheader | paste -sd, -)
stamp() { date -u +%Y%m%dT%H%M%SZ; }
done_ft() { $PY -c "import json,sys; sys.exit(json.load(open('$1/fine_tuning.json'))['state'] != 'complete')" 2>/dev/null; }
done_fit() { $PY -c "import json,sys; sys.exit(json.load(open('$1/fit_manifest.json'))['states']['process'] != 'complete')" 2>/dev/null; }
failed=0
for d in "$F"/*_grain_cross10; do
    [ -d "$d" ] || continue
    name=$(basename "$d"); out=$O/${name}_m5_i5_cpu
    done_fit "$d" || { echo "$(date -u +%FT%TZ) $name: fit not complete, skipped"; failed=1; continue; }
    if done_ft "$out"; then echo "$(date -u +%FT%TZ) $name: fine tuning complete, skipped"; continue; fi
    [ -e "$out" ] && mv "$out" "$out.incomplete_$(stamp)"
    [ -e "$L/ft_$name.log" ] && mv "$L/ft_$name.log" "$L/ft_$name.log.incomplete_$(stamp)"
    echo "$(date -u +%FT%TZ) $name: fine tuning start"
    $PY "$X/grain_finetune.py" "$R/base" --restart-from "$d/fit_restart.npz" --minutes 5 --intervals 5 \
        --method cpu --max-nfev 2000 --devices "$DEV" --precision mixed --out "$out" > "$L/ft_$name.log" 2>&1
    rc=$?
    echo "$(date -u +%FT%TZ) $name: exit $rc"
    [ $rc -eq 0 ] || failed=1
done
[ $failed -eq 0 ] && echo "grain_ft: all complete" || echo "grain_ft: something failed or was skipped" >&2
exit $failed
