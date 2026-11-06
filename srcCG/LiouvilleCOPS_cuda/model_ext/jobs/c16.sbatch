#!/bin/bash
# ==== version-controlled copy: header added 2026-11-06, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/jobs/c16.sbatch (on Voyager: /data001/highresolu_vt/cops_output/jobs/c16.sbatch)
# SHA256 of the run copy = SHA256 of this file without these header lines: 4ab2d2826966de7a100df18df1f2e0454e2d2a7bbfb878a328b3ad728a26920b
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: the C16 fit set (B29, N30, S30, NS31 with --monotonic-weight 10 --ordered-temperature --tbounds 30 400, mixed,
#   max_nfev 2000) on a finished runner campaign; builds the sin(slope) templates first. Jobs 120870, 120873.
# RUN: sbatch --export=ALL,COPS_RUN=<run> c16.sbatch (submit_campaign.sh does this)
# ==== end of header ====
# C16 fit set (handoff §8.2; user decision 2026-10-27: part of every campaign from now on) on a runner campaign
# whose prepare group is complete: B29, N30, S30, NS31 with the passion C16 settings (--monotonic-weight 10
# --ordered-temperature --tbounds 30 400, mixed, max_nfev 2000). The runner has no C16 options, so this job runs
# tools/gpu_fit.py directly on the runner's base, record store and bundle.
#   sbatch --export=ALL,COPS_RUN=<runner restart dir> c16.sbatch      (submit_campaign.sh does this)
# Products: $COPS_RUN/c16/{bundle/slope_templates.npz, fits/<fit>/, logs/, provenance-<jobid>/}.
# A finished fit (fit_manifest.json states.process == complete) is skipped, so a resubmission resumes;
# an unfinished one is moved aside as *.incomplete_<UTC> and restarted from the default start.
#SBATCH --job-name=cops_c16
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gpus-per-node=8
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=06:00:00
#SBATCH --partition=gpu
#SBATCH --account=general
#SBATCH --output=/data001/highresolu_vt/cops_output/jobs/cops_c16-%j.out
set -uo pipefail
APP=/data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda
R=${COPS_RUN:?set COPS_RUN to the runner restart dir}
PY="$APP/.venv/bin/python -W ignore"
C=$R/c16; B=$C/bundle/slope_templates.npz; O=$C/fits; L=$C/logs; P=$C/provenance-$SLURM_JOB_ID
mkdir -p "$C/bundle" "$O" "$L" "$P"
cd "$APP" || exit 1
cp -p "$0" "$P/c16.sbatch"
{ git -C "$APP" rev-parse HEAD; git -C "$APP" status --short; } > "$P/git.txt" 2>&1
module list > "$P/module_list.txt" 2>&1 || true
nvidia-smi -q > "$P/nvidia_smi_q.txt" 2>&1
{ nvcc --version; gcc --version | head -1; grep -m1 CSPICE /home/vtenishev/SPICE/cspice/doc/version.txt; } > "$P/toolchain.txt" 2>&1
~/.local/bin/uv pip list --python "$APP/.venv/bin/python" > "$P/pip_list.txt" 2>&1
env | grep -E '^(CUDA_|CUPY_|COPS_|OMP_|OPENBLAS_|MKL_|PYTHONUNBUFFERED|LD_LIBRARY_PATH|PATH)=' | sort > "$P/environment.txt"
for s in base store bundle; do
    [ -d "$R/$s" ] || { echo "c16: $R/$s missing (prepare group not complete)" >&2; exit 2; }
done
DEV=$(nvidia-smi --query-gpu=index --format=csv,noheader | paste -sd, -)
stamp() { date -u +%Y%m%dT%H%M%SZ; }

# C16 bundle: slope_sin from the runner bundle's gravity cache, late_ng_not_rg copied from the runner bundle
if [ ! -f "$B" ]; then
    echo "$(date -u +%FT%TZ) slope template: start"
    $PY tools/c16_slope_template.py --run "$R/base" --bundle "$R/bundle/restart_templates.npz" --out "$B" \
        > "$L/slope_template.log" 2>&1 || { echo "c16: slope template failed, see $L/slope_template.log" >&2; exit 1; }
    echo "$(date -u +%FT%TZ) slope template: done"
fi

FITARGS="--profile all_available --degree 4 --enhancement-degree 4 --lambda 0 --jobs 1 --tbounds 30 400 --ftol 1e-10 --xtol 1e-10 --gtol 1e-10 --check-tolerance 1e-6 --ordered-temperature"
G="--devices $DEV --precision mixed --record-store $R/store --monotonic-weight 10 --max-nfev 2000"
done_fit() { $PY -c "import json,sys; sys.exit(json.load(open('$1/fit_manifest.json'))['states']['process'] != 'complete')" 2>/dev/null; }
failed=0
fit() {   # fit <name> [template arguments]
    local name=$1; shift
    if done_fit "$O/$name"; then echo "$(date -u +%FT%TZ) $name: complete, skipped"; return 0; fi
    [ -e "$O/$name" ] && mv "$O/$name" "$O/$name.incomplete_$(stamp)"
    [ -e "$L/$name.log" ] && mv "$L/$name.log" "$L/$name.log.incomplete_$(stamp)"
    echo "$(date -u +%FT%TZ) $name: start"
    $PY tools/gpu_fit.py $G -- "$R/base" $FITARGS --out "$O/$name" "$@" > "$L/$name.log" 2>&1
    local rc=$?
    echo "$(date -u +%FT%TZ) $name: exit $rc"
    [ $rc -eq 0 ] || failed=1
}
fit S30_n2000_mono10_ordered  --activity-template "$B" slope_sin
fit B29_n2000_mono10_ordered
fit N30_n2000_mono10_ordered  --activity-template "$B" late_ng_not_rg
fit NS31_n2000_mono10_ordered --activity-template "$B" late_ng_not_rg --activity-template "$B" slope_sin
[ $failed -eq 0 ] && echo "c16: all fits complete" || echo "c16: a fit failed" >&2
exit $failed
