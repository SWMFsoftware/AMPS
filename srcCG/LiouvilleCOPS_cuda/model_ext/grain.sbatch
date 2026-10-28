#!/bin/bash
# 12M refits with the dust-grain term in NG and the NG = RG crossing constraint (user request 2026-10-27;
# model_ext/grain_crossing.py, driver model_ext/grain_fit.py; the git tree is unchanged).
#   sbatch --export=ALL,COPS_RUN=<runner restart dir>,GRAIN_SET=free|mono[,GRAIN_TEST=1] grain.sbatch
# GRAIN_SET=free: the free-temperature fits (no template, late NG-only, cliff >=60, late NG-only + cliff >=60);
# GRAIN_SET=mono: monotonic penalty 10 (on the full NG) + ordered temperatures (no template, late NG-only,
#   sin(slope), late NG-only + sin(slope)). Both: grain A, B (start A = observed jump, B = 0), crossing weight 10,
#   mixed precision, --tbounds 30 400, max_nfev 2000, default start, reusing the run's base, store and templates.
# GRAIN_TEST=1: one fit (late NG-only + sin(slope) with the mono settings, else late NG-only + cliff) with
#   max_nfev 4 and a finite-difference check of the grain Jacobian columns, into <run>/grain/test-<jobid>/.
# Products: <run>/grain/fits/<fit>/ (+ grain_term.{txt,json}, grain_components.dat), logs/, provenance-<jobid>/.
# A finished fit is skipped; an unfinished one is moved aside as *.incomplete_<UTC> and restarted.
#SBATCH --job-name=cops_grain
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gpus-per-node=8
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=08:00:00
#SBATCH --partition=gpu
#SBATCH --account=general
#SBATCH --output=/data001/highresolu_vt/cops_output/model_ext/logs/cops_grain-%j.out
set -uo pipefail
export AMPS=/data001/highresolu_vt/AMPS APP=/data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda
export T=/data001/highresolu_vt/cops_cluster_transfer_20261026 COPS_OUTPUT_ROOT=/data001/highresolu_vt/cops_output
export COPS_DATA_ROOT=$T/rosetta_data_end_of_mission_cops
export CUDA_PATH=/usr/local/cuda-12.9 CUDA_HOME=/usr/local/cuda-12.9 CUPY_CACHE_DIR=$COPS_OUTPUT_ROOT/.cupy
export PATH=$CUDA_HOME/bin:/usr/mpi/gcc/openmpi-4.1.7rc1/bin:$PATH
export LD_LIBRARY_PATH=$CUDA_HOME/lib64${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1
R=${COPS_RUN:?set COPS_RUN}; SET=${GRAIN_SET:?set GRAIN_SET=free or mono}; TEST=${GRAIN_TEST:-0}
X=$COPS_OUTPUT_ROOT/model_ext
PY="$APP/.venv/bin/python -W ignore"
G=$R/grain; O=$G/fits; L=$G/logs; P=$G/provenance-$SLURM_JOB_ID
[ "$TEST" = 1 ] && { O=$G/test-$SLURM_JOB_ID; L=$O/logs; }
mkdir -p "$O" "$L" "$P"
cd "$APP" || exit 1
cp -p "$0" "$P/grain.sbatch"; cp -p "$X/grain_crossing.py" "$X/grain_fit.py" "$P/"; sha256sum "$X"/*.py > "$P/model_ext.sha256"
{ git -C "$AMPS" rev-parse HEAD; git -C "$AMPS" status --short; } > "$P/git.txt" 2>&1
module list > "$P/module_list.txt" 2>&1 || true
nvidia-smi -q > "$P/nvidia_smi_q.txt" 2>&1
{ nvcc --version; gcc --version | head -1; grep -m1 CSPICE /home/vtenishev/SPICE/cspice/doc/version.txt; } > "$P/toolchain.txt" 2>&1
~/.local/bin/uv pip list --python "$APP/.venv/bin/python" > "$P/pip_list.txt" 2>&1
env | grep -E '^(CUDA_|CUPY_|COPS_|GRAIN_|OMP_|OPENBLAS_|MKL_|PYTHONUNBUFFERED|LD_LIBRARY_PATH|PATH)=' | sort > "$P/environment.txt"
for s in base store bundle c16/bundle; do
    [ -e "$R/$s" ] || { echo "grain: $R/$s missing" >&2; exit 2; }
done
DEV=$(nvidia-smi --query-gpu=index --format=csv,noheader | paste -sd, -)
B=$R/bundle/restart_templates.npz; S=$R/c16/bundle/slope_templates.npz
NFEV=2000; CHECK=--check-jacobian; [ "$TEST" = 1 ] && NFEV=4   # the check costs two residual evaluations
FITARGS="--profile all_available --degree 4 --enhancement-degree 4 --lambda 0 --jobs 1 --tbounds 30 400 --ftol 1e-10 --xtol 1e-10 --gtol 1e-10 --check-tolerance 1e-6"
GA="--crossing-weight 10 $CHECK --devices $DEV --precision mixed --record-store $R/store --max-nfev $NFEV --wall-time slurm"
if [ "$SET" = mono ]; then GA="$GA --monotonic-weight 10"; FITARGS="$FITARGS --ordered-temperature"; fi
stamp() { date -u +%Y%m%dT%H%M%SZ; }
done_fit() { $PY -c "import json,sys; sys.exit(json.load(open('$1/fit_manifest.json'))['states']['process'] != 'complete')" 2>/dev/null; }
failed=0
fit() {   # fit <name> [template arguments]
    local name=$1; shift
    if done_fit "$O/$name"; then echo "$(date -u +%FT%TZ) $name: complete, skipped"; return 0; fi
    [ -e "$O/$name" ] && mv "$O/$name" "$O/$name.incomplete_$(stamp)"
    [ -e "$L/$name.log" ] && mv "$L/$name.log" "$L/$name.log.incomplete_$(stamp)"
    echo "$(date -u +%FT%TZ) $name: start"
    $PY "$X/grain_fit.py" $GA -- "$R/base" $FITARGS --out "$O/$name" "$@" > "$L/$name.log" 2>&1
    local rc=$?
    echo "$(date -u +%FT%TZ) $name: exit $rc"
    [ $rc -eq 0 ] || failed=1
}
if [ "$SET" = mono ]; then
    if [ "$TEST" = 1 ]; then
        fit NS31_n4_mono10_ordered_grain_cross10 --activity-template "$S" late_ng_not_rg --activity-template "$S" slope_sin
    else
        fit B29_n2000_mono10_ordered_grain_cross10
        fit N30_n2000_mono10_ordered_grain_cross10 --activity-template "$S" late_ng_not_rg
        fit S30_n2000_mono10_ordered_grain_cross10 --activity-template "$S" slope_sin
        fit NS31_n2000_mono10_ordered_grain_cross10 --activity-template "$S" late_ng_not_rg --activity-template "$S" slope_sin
    fi
else
    if [ "$TEST" = 1 ]; then
        fit NC31_n4_grain_cross10 --activity-template "$B" late_ng_not_rg --activity-template "$B" cliff_ge_60deg
    else
        fit B29_n2000_grain_cross10
        fit N30_n2000_grain_cross10 --activity-template "$B" late_ng_not_rg
        fit C30_n2000_grain_cross10 --activity-template "$B" cliff_ge_60deg
        fit NC31_n2000_grain_cross10 --activity-template "$B" late_ng_not_rg --activity-template "$B" cliff_ge_60deg
    fi
fi
[ $failed -eq 0 ] && echo "grain: all fits complete ($SET)" || echo "grain: a fit failed ($SET)" >&2
exit $failed
